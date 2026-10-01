/* This file is part of msolve.
 *
 * msolve is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * msolve is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with msolve.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Jérémy Berthomieu
 * Christian Eder
 * Mohab Safey El Din */

/* Linear algebra for prime fields with 2^32 <= p < 2^52.
 *
 * Products of coefficients do not fit into 64 bits anymore, so the dense
 * rows cannot accumulate unreduced products as for the smaller primes.
 * Instead, all entries of the dense rows are kept in [0,p) and each update
 * d = d - a*m mod p is reduced at once. Here m is the multiplier of the
 * pivot row, it is fixed while the pivot row is applied, so it is
 * prepared once per pivot row:
 *
 * - Scalar code uses Montgomery reduction with R = 2^64, the multiplier is
 *   stored as m*R mod p.
 * - AVX512 IFMA resp. AVX IFMA use Montgomery reduction with R = 2^52 on
 *   the 52 bit multipliers of vpmadd52luq and vpmadd52huq.
 * - AVX2, AVX512F and NEON compute in double precision with FMA: for
 *   a, m < p < 2^52 the product a*m = h + l is computed exactly as the
 *   rounded product h and the error l = fma(a, m, -h). With the quotient
 *   q = round(h/p) the remainder h - q*p + l is exact again, its absolute
 *   value is < 2p. Rounding once more brings it to at most p/2+2 < 2^51.
 */

#include "data.h"
#include "../msolve/streams.h"

#if defined HAVE_AVX2 || defined HAVE_AVX512_F
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

/* rounding mode for _mm*_round*_pd(), clang needs it as a constant
 * expression, a const int variable is not enough */
#define FF64_RN (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC)

/* data for the arithmetic modulo p */
typedef struct {
    uint64_t p;
    uint64_t pi;    /* p^(-1) mod 2^64 */
    uint64_t r128;  /* 2^128 mod p */
    uint64_t r104;  /* 2^104 mod p */
    double pinv;    /* 1/p */
} ff64_t;

static inline ff64_t ff64_init(
        const uint64_t p
        )
{
    ff64_t f;
    /* p is odd, so p*p = 1 mod 8, each Newton step doubles the number
     * of correct bits of the inverse */
    uint64_t x = p;
    for (int i = 0; i < 5; ++i) {
        x *= 2 - p * x;
    }
    f.p     = p;
    f.pi    = x;
    const uint64_t r64  = (UINT64_MAX % p + 1) % p;
    f.r128  = mul_mod_64(r64, r64, p);
    const uint64_t r52  = ((uint64_t)1 << 52) % p;
    f.r104  = mul_mod_64(r52, r52, p);
    f.pinv  = 1.0 / (double)p;
    return f;
}

/* returns t = a*b*2^(-64) mod p with -p < t < p, needs a*b < p*2^64 */
static inline int64_t ff64_redc(
        const uint64_t a,
        const uint64_t b,
        const ff64_t * const f
        )
{
    const unsigned __int128 t = (unsigned __int128)a * b;
    const uint64_t m = (uint64_t)t * f->pi;
    return (int64_t)((uint64_t)(t >> 64)
            - (uint64_t)(((unsigned __int128)m * f->p) >> 64));
}

/* returns a*2^64 mod p for a < p, the multiplier for ff64_mul() and
 * ff64_submul() */
static inline uint64_t ff64_mont(
        const uint64_t a,
        const ff64_t * const f
        )
{
    const int64_t t = ff64_redc(a, f->r128, f);
    return (uint64_t)(t + ((t >> 63) & (int64_t)f->p));
}

/* returns a*m mod p for a < p, mm = ff64_mont(m) */
static inline uint64_t ff64_mul(
        const uint64_t a,
        const uint64_t mm,
        const ff64_t * const f
        )
{
    const int64_t t = ff64_redc(a, mm, f);
    return (uint64_t)(t + ((t >> 63) & (int64_t)f->p));
}

/* returns d - a*m mod p for d, a < p, mm = ff64_mont(m) */
static inline int64_t ff64_submul(
        const int64_t d,
        const uint64_t a,
        const uint64_t mm,
        const ff64_t * const f
        )
{
    const int64_t p = (int64_t)f->p;
    int64_t x = d - ff64_redc(a, mm, f);
    x +=  (x >> 63) & p;
    x -=  p;
    x +=  (x >> 63) & p;
    return x;
}

/* returns a*2^52 mod p for a < p, the multiplier for the IFMA kernels */
static inline uint64_t ff64_mont52(
        const uint64_t a,
        const ff64_t * const f
        )
{
    const uint64_t mask = ((uint64_t)1 << 52) - 1;
    const unsigned __int128 t = (unsigned __int128)a * f->r104;
    const uint64_t m  = ((uint64_t)t * f->pi) & mask;
    const int64_t r   = (int64_t)((uint64_t)(t >> 52)
            - (uint64_t)(((unsigned __int128)m * f->p) >> 52));
    return (uint64_t)(r + ((r >> 63) & (int64_t)f->p));
}

/* Vector arithmetic: ff64_v holds FF64_VL entries of a dense row resp.
 * coefficients in [0,p) as 64 bit integers, ff64_m holds the multipliers
 * in the form needed by ff64_v_submul(), ff64_d holds FF64_VL doubles. */
#if defined HAVE_AVX512_F
#define FF64_VL 8
typedef __m512i ff64_v;
typedef __m512d ff64_d;
#if defined HAVE_AVX512_IFMA
#define FF64_IFMA 1
typedef __m512i ff64_m;
#else
typedef __m512d ff64_m;
#endif
#elif defined HAVE_AVX2 && (defined __AVXIFMA__ || defined HAVE_FMA3)
#define FF64_VL 4
typedef __m256i ff64_v;
typedef __m256d ff64_d;
#if defined __AVXIFMA__
#define FF64_IFMA 1
typedef __m256i ff64_m;
#else
typedef __m256d ff64_m;
#endif
#elif defined __aarch64__
#define FF64_VL 2
typedef int64x2_t ff64_v;
typedef float64x2_t ff64_m;
typedef float64x2_t ff64_d;
#endif

#ifdef FF64_VL
/* constants for the vector arithmetic */
typedef struct {
    ff64_v p;
#if defined FF64_IFMA
    ff64_v pi;      /* p^(-1) mod 2^52 */
    ff64_v r104;
#else
    ff64_m pd;      /* p as double */
    ff64_m pinv;    /* 1/p */
#endif
#if defined HAVE_AVX2 && !defined HAVE_AVX512_F
    ff64_v pm1;     /* p-1 for comparisons */
#endif
} ff64_vc_t;
#endif

#if defined HAVE_AVX512_F
static inline ff64_vc_t ff64_vc_init(
        const ff64_t * const f
        )
{
    ff64_vc_t c;
    c.p     = _mm512_set1_epi64((int64_t)f->p);
#if defined FF64_IFMA
    c.pi    = _mm512_set1_epi64((int64_t)(f->pi & (((uint64_t)1 << 52) - 1)));
    c.r104  = _mm512_set1_epi64((int64_t)f->r104);
#else
    c.pd    = _mm512_set1_pd((double)f->p);
    c.pinv  = _mm512_set1_pd(f->pinv);
#endif
    return c;
}

static inline ff64_v ff64_v_load(
        const void * const p
        )
{
    return _mm512_loadu_si512(p);
}

static inline void ff64_v_store(
        void *p,
        const ff64_v v
        )
{
    _mm512_storeu_si512(p, v);
}

static inline ff64_v ff64_v_set1(
        const uint64_t a
        )
{
    return _mm512_set1_epi64((int64_t)a);
}

static inline ff64_v ff64_v_gather(
        const int64_t * const dr,
        const hm_t * const ds
        )
{
    return _mm512_i32gather_epi64(
            _mm256_loadu_si256((__m256i *)ds), (const void *)dr, 8);
}

static inline void ff64_v_scatter(
        int64_t *dr,
        const hm_t * const ds,
        const ff64_v v
        )
{
    _mm512_i32scatter_epi64((void *)dr,
            _mm256_loadu_si256((__m256i *)ds), v, 8);
}

static inline __m512d ff64_cvt_pd(
        const __m512i a
        )
{
#if defined HAVE_AVX512_DQ
    return _mm512_cvtepi64_pd(a);
#else
    /* exact for 0 <= a < 2^52 */
    const __m512i mg  = _mm512_set1_epi64(0x4330000000000000);
    return _mm512_sub_pd(_mm512_castsi512_pd(_mm512_or_si512(a, mg)),
            _mm512_castsi512_pd(mg));
#endif
}

static inline __m512i ff64_cvt_epi64(
        const __m512d r
        )
{
#if defined HAVE_AVX512_DQ
    return _mm512_cvtpd_epi64(r);
#else
    /* exact for integral r with |r| <= 2^51 */
    const __m512d mg  = _mm512_set1_pd(6755399441055744.0);
    return _mm512_sub_epi64(_mm512_castpd_si512(_mm512_add_pd(r, mg)),
            _mm512_castpd_si512(mg));
#endif
}

#if defined FF64_IFMA
/* multipliers m < p in Montgomery form w.r.t. 2^52 */
static inline ff64_m ff64_m_from(
        const ff64_v m,
        const ff64_vc_t * const c
        )
{
    const __m512i z   = _mm512_setzero_si512();
    const __m512i tl  = _mm512_madd52lo_epu64(z, m, c->r104);
    const __m512i th  = _mm512_madd52hi_epu64(z, m, c->r104);
    const __m512i q   = _mm512_madd52lo_epu64(z, tl, c->pi);
    const __m512i x   = _mm512_sub_epi64(th, _mm512_madd52hi_epu64(z, q, c->p));
    return _mm512_min_epu64(x, _mm512_add_epi64(x, c->p));
}

static inline ff64_m ff64_m_set1(
        const uint64_t m,
        const ff64_t * const f
        )
{
    return _mm512_set1_epi64((int64_t)ff64_mont52(m, f));
}

/* d - a*m mod p: with t = a*mm = th*2^52 + tl and q = tl/p mod 2^52 the
 * low 52 bits of t and q*p agree, so a*m = th - hi(q*p) mod p */
static inline ff64_v ff64_v_submul(
        const ff64_v d,
        const ff64_v a,
        const ff64_m mm,
        const ff64_vc_t * const c
        )
{
    const __m512i z   = _mm512_setzero_si512();
    const __m512i tl  = _mm512_madd52lo_epu64(z, a, mm);
    const __m512i th  = _mm512_madd52hi_epu64(z, a, mm);
    const __m512i q   = _mm512_madd52lo_epu64(z, tl, c->pi);
    /* -p < x < 2p */
    __m512i x = _mm512_sub_epi64(_mm512_madd52hi_epu64(d, q, c->p), th);
    x = _mm512_min_epu64(x, _mm512_add_epi64(x, c->p));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, c->p));
}
#else
static inline ff64_m ff64_m_from(
        const ff64_v m,
        const ff64_vc_t * const c
        )
{
    (void)c;
    return ff64_cvt_pd(m);
}

static inline ff64_m ff64_m_set1(
        const uint64_t m,
        const ff64_t * const f
        )
{
    (void)f;
    return _mm512_set1_pd((double)m);
}

static inline ff64_v ff64_v_submul(
        const ff64_v d,
        const ff64_v a,
        const ff64_m m,
        const ff64_vc_t * const c
        )
{
    const __m512d ad  = ff64_cvt_pd(a);
    const __m512d h   = _mm512_mul_pd(ad, m);
    const __m512d l   = _mm512_fmsub_pd(ad, m, h);
    const __m512d q   = _mm512_roundscale_pd(_mm512_mul_pd(h, c->pinv), FF64_RN);
    /* r = a*m - q*p, |r| < 2p */
    __m512d r = _mm512_add_pd(_mm512_fnmadd_pd(q, c->pd, h), l);
    const __m512d k   = _mm512_roundscale_pd(_mm512_mul_pd(r, c->pinv), FF64_RN);
    r = _mm512_fnmadd_pd(k, c->pd, r);
    /* -p < x < 2p */
    __m512i x = _mm512_sub_epi64(d, ff64_cvt_epi64(r));
    x = _mm512_min_epu64(x, _mm512_add_epi64(x, c->p));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, c->p));
}
#endif

#elif defined FF64_VL && defined HAVE_AVX2
static inline ff64_vc_t ff64_vc_init(
        const ff64_t * const f
        )
{
    ff64_vc_t c;
    c.p     = _mm256_set1_epi64x((int64_t)f->p);
    c.pm1   = _mm256_set1_epi64x((int64_t)f->p - 1);
#if defined FF64_IFMA
    c.pi    = _mm256_set1_epi64x((int64_t)(f->pi & (((uint64_t)1 << 52) - 1)));
    c.r104  = _mm256_set1_epi64x((int64_t)f->r104);
#else
    c.pd    = _mm256_set1_pd((double)f->p);
    c.pinv  = _mm256_set1_pd(f->pinv);
#endif
    return c;
}

static inline ff64_v ff64_v_load(
        const void * const p
        )
{
    return _mm256_loadu_si256((const __m256i *)p);
}

static inline void ff64_v_store(
        void *p,
        const ff64_v v
        )
{
    _mm256_storeu_si256((__m256i *)p, v);
}

static inline ff64_v ff64_v_set1(
        const uint64_t a
        )
{
    return _mm256_set1_epi64x((int64_t)a);
}

static inline ff64_v ff64_v_gather(
        const int64_t * const dr,
        const hm_t * const ds
        )
{
    return _mm256_i32gather_epi64((const long long *)dr,
            _mm_loadu_si128((const __m128i *)ds), 8);
}

static inline void ff64_v_scatter(
        int64_t *dr,
        const hm_t * const ds,
        const ff64_v v
        )
{
    int64_t t[4] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i *)t, v);
    dr[ds[0]] = t[0];
    dr[ds[1]] = t[1];
    dr[ds[2]] = t[2];
    dr[ds[3]] = t[3];
}

/* returns x mod p for -p < x < 2p */
static inline __m256i ff64_v_normalize(
        __m256i x,
        const ff64_vc_t * const c
        )
{
    x = _mm256_add_epi64(x, _mm256_and_si256(
                _mm256_cmpgt_epi64(_mm256_setzero_si256(), x), c->p));
    return _mm256_sub_epi64(x, _mm256_and_si256(
                _mm256_cmpgt_epi64(x, c->pm1), c->p));
}

static inline __m256d ff64_cvt_pd(
        const __m256i a
        )
{
    /* exact for 0 <= a < 2^52 */
    const __m256i mg  = _mm256_set1_epi64x(0x4330000000000000);
    return _mm256_sub_pd(_mm256_castsi256_pd(_mm256_or_si256(a, mg)),
            _mm256_castsi256_pd(mg));
}

static inline __m256i ff64_cvt_epi64(
        const __m256d r
        )
{
    /* exact for integral r with |r| <= 2^51 */
    const __m256d mg  = _mm256_set1_pd(6755399441055744.0);
    return _mm256_sub_epi64(_mm256_castpd_si256(_mm256_add_pd(r, mg)),
            _mm256_castpd_si256(mg));
}

#if defined FF64_IFMA
static inline ff64_m ff64_m_from(
        const ff64_v m,
        const ff64_vc_t * const c
        )
{
    const __m256i z   = _mm256_setzero_si256();
    const __m256i tl  = _mm256_madd52lo_avx_epu64(z, m, c->r104);
    const __m256i th  = _mm256_madd52hi_avx_epu64(z, m, c->r104);
    const __m256i q   = _mm256_madd52lo_avx_epu64(z, tl, c->pi);
    const __m256i x   = _mm256_sub_epi64(th,
            _mm256_madd52hi_avx_epu64(z, q, c->p));
    return _mm256_add_epi64(x, _mm256_and_si256(_mm256_cmpgt_epi64(z, x), c->p));
}

static inline ff64_m ff64_m_set1(
        const uint64_t m,
        const ff64_t * const f
        )
{
    return _mm256_set1_epi64x((int64_t)ff64_mont52(m, f));
}

/* see the AVX512 IFMA version */
static inline ff64_v ff64_v_submul(
        const ff64_v d,
        const ff64_v a,
        const ff64_m mm,
        const ff64_vc_t * const c
        )
{
    const __m256i z   = _mm256_setzero_si256();
    const __m256i tl  = _mm256_madd52lo_avx_epu64(z, a, mm);
    const __m256i th  = _mm256_madd52hi_avx_epu64(z, a, mm);
    const __m256i q   = _mm256_madd52lo_avx_epu64(z, tl, c->pi);
    return ff64_v_normalize(_mm256_sub_epi64(
                _mm256_madd52hi_avx_epu64(d, q, c->p), th), c);
}
#else
static inline ff64_m ff64_m_from(
        const ff64_v m,
        const ff64_vc_t * const c
        )
{
    (void)c;
    return ff64_cvt_pd(m);
}

static inline ff64_m ff64_m_set1(
        const uint64_t m,
        const ff64_t * const f
        )
{
    (void)f;
    return _mm256_set1_pd((double)m);
}

static inline ff64_v ff64_v_submul(
        const ff64_v d,
        const ff64_v a,
        const ff64_m m,
        const ff64_vc_t * const c
        )
{
    const __m256d ad  = ff64_cvt_pd(a);
    const __m256d h   = _mm256_mul_pd(ad, m);
    const __m256d l   = _mm256_fmsub_pd(ad, m, h);
    const __m256d q   = _mm256_round_pd(_mm256_mul_pd(h, c->pinv), FF64_RN);
    /* r = a*m - q*p, |r| < 2p */
    __m256d r = _mm256_add_pd(_mm256_fnmadd_pd(q, c->pd, h), l);
    const __m256d k   = _mm256_round_pd(_mm256_mul_pd(r, c->pinv), FF64_RN);
    r = _mm256_fnmadd_pd(k, c->pd, r);
    return ff64_v_normalize(_mm256_sub_epi64(d, ff64_cvt_epi64(r)), c);
}
#endif

#elif defined FF64_VL && defined __aarch64__
static inline ff64_vc_t ff64_vc_init(
        const ff64_t * const f
        )
{
    ff64_vc_t c;
    c.p     = vdupq_n_s64((int64_t)f->p);
    c.pd    = vdupq_n_f64((double)f->p);
    c.pinv  = vdupq_n_f64(f->pinv);
    return c;
}

static inline ff64_v ff64_v_load(
        const void * const p
        )
{
    return vld1q_s64((const int64_t *)p);
}

static inline void ff64_v_store(
        void *p,
        const ff64_v v
        )
{
    vst1q_s64((int64_t *)p, v);
}

static inline ff64_v ff64_v_set1(
        const uint64_t a
        )
{
    return vdupq_n_s64((int64_t)a);
}

static inline ff64_v ff64_v_gather(
        const int64_t * const dr,
        const hm_t * const ds
        )
{
    return vcombine_s64(vld1_s64(dr+ds[0]), vld1_s64(dr+ds[1]));
}

static inline void ff64_v_scatter(
        int64_t *dr,
        const hm_t * const ds,
        const ff64_v v
        )
{
    vst1q_lane_s64(dr+ds[0], v, 0);
    vst1q_lane_s64(dr+ds[1], v, 1);
}

static inline ff64_m ff64_m_from(
        const ff64_v m,
        const ff64_vc_t * const c
        )
{
    (void)c;
    return vcvtq_f64_s64(m);
}

static inline ff64_m ff64_m_set1(
        const uint64_t m,
        const ff64_t * const f
        )
{
    (void)f;
    return vdupq_n_f64((double)m);
}

/* see the AVX2 version */
static inline ff64_v ff64_v_submul(
        const ff64_v d,
        const ff64_v a,
        const ff64_m m,
        const ff64_vc_t * const c
        )
{
    const float64x2_t ad  = vcvtq_f64_s64(a);
    const float64x2_t h   = vmulq_f64(ad, m);
    /* nl = h - a*m */
    const float64x2_t nl  = vfmsq_f64(h, ad, m);
    const float64x2_t q   = vrndnq_f64(vmulq_f64(h, c->pinv));
    /* r = a*m - q*p, |r| < 2p */
    float64x2_t r = vsubq_f64(vfmsq_f64(h, q, c->pd), nl);
    const float64x2_t k   = vrndnq_f64(vmulq_f64(r, c->pinv));
    r = vfmsq_f64(r, k, c->pd);
    /* -p < x < 2p */
    int64x2_t x = vsubq_s64(d, vcvtq_s64_f64(r));
    x = vaddq_s64(x, vandq_s64(vreinterpretq_s64_u64(vcltzq_s64(x)), c->p));
    return vsubq_s64(x, vandq_s64(
                vreinterpretq_s64_u64(vcgeq_s64(x, c->p)), c->p));
}
#endif

/* used unless AVX512 IFMA is available, with AVX IFMA the 256 bit
 * integer version needs as many instructions as the one below */
#if defined FF64_VL && !(defined HAVE_AVX512_F && defined FF64_IFMA) \
    && (defined HAVE_AVX512_F || defined HAVE_FMA3 || defined __aarch64__)
/* Arithmetic in double precision on the interleaved rows of
 * reduce_dense_rows_by_known_pivots_ff_64(): the entries x are integers
 * stored as doubles with |x| <= p/2+2. For an update x - a*m with a < p
 * and |m| <= p/2+2, the exact product a*m = h + l is split as for
 * ff64_v_submul(). Since |a*m/p| < 2^51, the quotient q = round(h/p)
 * differs by at most 1 from a*m/p (|l|/p and the error of 1/p are both
 * < 1/4), so r = a*m - q*p satisfies |r| <= p. Thus x - r is exact
 * (<= 3p/2+2 < 2^53) and subtracting the closest multiple of p yields
 * |x| <= p/2+2 again. On x86 round(y/p) is computed as
 * fma(y, 1/p, 1.5*2^52) - 1.5*2^52, the fma rounds to an integer since
 * |y/p| < 2^51. This avoids the conversions and corrections of
 * ff64_v_submul(). */
#define FF64_DOUBLE_ROWS 1

typedef struct {
    ff64_d pd;
    ff64_d pinv;
#if defined HAVE_AVX2 || defined HAVE_AVX512_F
    ff64_d mg;      /* 1.5*2^52, fma(y, 1/p, mg) - mg is round(y/p) */
#endif
} ff64_dc_t;

#if defined HAVE_AVX512_F
static inline ff64_dc_t ff64_dc_init(
        const ff64_t * const f
        )
{
    ff64_dc_t c;
    c.pd    = _mm512_set1_pd((double)f->p);
    c.pinv  = _mm512_set1_pd(f->pinv);
    c.mg    = _mm512_set1_pd(6755399441055744.0);
    return c;
}

static inline ff64_d ff64_d_load(
        const int64_t * const p
        )
{
    return _mm512_loadu_pd((const void *)p);
}

static inline void ff64_d_store(
        int64_t *p,
        const ff64_d v
        )
{
    _mm512_storeu_pd((void *)p, v);
}

static inline ff64_d ff64_d_set1(
        const double a
        )
{
    return _mm512_set1_pd(a);
}

static inline int ff64_d_nonzero(
        const ff64_d v
        )
{
    return _mm512_cmp_pd_mask(v, _mm512_setzero_pd(), _CMP_NEQ_OQ) != 0;
}

/* integers 0 <= x < p to doubles |x| <= p/2+2 */
static inline ff64_d ff64_d_from_int(
        const ff64_v v,
        const ff64_dc_t * const c
        )
{
    const __m512d x   = ff64_cvt_pd(v);
    const __m512d k   = _mm512_sub_pd(_mm512_fmadd_pd(x, c->pinv, c->mg), c->mg);
    return _mm512_fnmadd_pd(k, c->pd, x);
}

/* doubles |x| <= p/2+2 to integers in [0,p) */
static inline ff64_v ff64_d_to_int(
        const ff64_d x,
        const ff64_vc_t * const c
        )
{
    const __m512i v   = ff64_cvt_epi64(x);
    return _mm512_min_epu64(v, _mm512_add_epi64(v, c->p));
}

static inline ff64_d ff64_d_submul(
        const ff64_d x,
        const ff64_d a,
        const ff64_d m,
        const ff64_dc_t * const c
        )
{
    const __m512d h   = _mm512_mul_pd(a, m);
    const __m512d l   = _mm512_fmsub_pd(a, m, h);
    const __m512d q   = _mm512_sub_pd(_mm512_fmadd_pd(h, c->pinv, c->mg), c->mg);
    const __m512d r   = _mm512_add_pd(_mm512_fnmadd_pd(q, c->pd, h), l);
    const __m512d y   = _mm512_sub_pd(x, r);
    const __m512d k   = _mm512_sub_pd(_mm512_fmadd_pd(y, c->pinv, c->mg), c->mg);
    return _mm512_fnmadd_pd(k, c->pd, y);
}
#elif defined HAVE_AVX2
static inline ff64_dc_t ff64_dc_init(
        const ff64_t * const f
        )
{
    ff64_dc_t c;
    c.pd    = _mm256_set1_pd((double)f->p);
    c.pinv  = _mm256_set1_pd(f->pinv);
    c.mg    = _mm256_set1_pd(6755399441055744.0);
    return c;
}

static inline ff64_d ff64_d_load(
        const int64_t * const p
        )
{
    return _mm256_loadu_pd((const double *)p);
}

static inline void ff64_d_store(
        int64_t *p,
        const ff64_d v
        )
{
    _mm256_storeu_pd((double *)p, v);
}

static inline ff64_d ff64_d_set1(
        const double a
        )
{
    return _mm256_set1_pd(a);
}

static inline int ff64_d_nonzero(
        const ff64_d v
        )
{
    return _mm256_movemask_pd(
            _mm256_cmp_pd(v, _mm256_setzero_pd(), _CMP_NEQ_OQ)) != 0;
}

static inline ff64_d ff64_d_from_int(
        const ff64_v v,
        const ff64_dc_t * const c
        )
{
    const __m256d x   = ff64_cvt_pd(v);
    const __m256d k   = _mm256_sub_pd(_mm256_fmadd_pd(x, c->pinv, c->mg), c->mg);
    return _mm256_fnmadd_pd(k, c->pd, x);
}

static inline ff64_v ff64_d_to_int(
        const ff64_d x,
        const ff64_vc_t * const c
        )
{
    const __m256i v   = ff64_cvt_epi64(x);
    return _mm256_add_epi64(v, _mm256_and_si256(
                _mm256_cmpgt_epi64(_mm256_setzero_si256(), v), c->p));
}

static inline ff64_d ff64_d_submul(
        const ff64_d x,
        const ff64_d a,
        const ff64_d m,
        const ff64_dc_t * const c
        )
{
    const __m256d h   = _mm256_mul_pd(a, m);
    const __m256d l   = _mm256_fmsub_pd(a, m, h);
    const __m256d q   = _mm256_sub_pd(_mm256_fmadd_pd(h, c->pinv, c->mg), c->mg);
    const __m256d r   = _mm256_add_pd(_mm256_fnmadd_pd(q, c->pd, h), l);
    const __m256d y   = _mm256_sub_pd(x, r);
    const __m256d k   = _mm256_sub_pd(_mm256_fmadd_pd(y, c->pinv, c->mg), c->mg);
    return _mm256_fnmadd_pd(k, c->pd, y);
}
#elif defined __aarch64__
static inline ff64_dc_t ff64_dc_init(
        const ff64_t * const f
        )
{
    ff64_dc_t c;
    c.pd    = vdupq_n_f64((double)f->p);
    c.pinv  = vdupq_n_f64(f->pinv);
    return c;
}

static inline ff64_d ff64_d_load(
        const int64_t * const p
        )
{
    return vreinterpretq_f64_s64(vld1q_s64(p));
}

static inline void ff64_d_store(
        int64_t *p,
        const ff64_d v
        )
{
    vst1q_s64(p, vreinterpretq_s64_f64(v));
}

static inline ff64_d ff64_d_set1(
        const double a
        )
{
    return vdupq_n_f64(a);
}

static inline int ff64_d_nonzero(
        const ff64_d v
        )
{
    const uint64x2_t z = vceqzq_f64(v);
    return (vgetq_lane_u64(z, 0) & vgetq_lane_u64(z, 1)) == 0;
}

static inline ff64_d ff64_d_from_int(
        const ff64_v v,
        const ff64_dc_t * const c
        )
{
    const float64x2_t x = vcvtq_f64_s64(v);
    return vfmsq_f64(x, vrndnq_f64(vmulq_f64(x, c->pinv)), c->pd);
}

static inline ff64_v ff64_d_to_int(
        const ff64_d x,
        const ff64_vc_t * const c
        )
{
    const int64x2_t v = vcvtq_s64_f64(x);
    return vaddq_s64(v, vandq_s64(vreinterpretq_s64_u64(vcltzq_s64(v)), c->p));
}

static inline ff64_d ff64_d_submul(
        const ff64_d x,
        const ff64_d a,
        const ff64_d m,
        const ff64_dc_t * const c
        )
{
    const float64x2_t h = vmulq_f64(a, m);
    /* nl = h - a*m */
    const float64x2_t nl  = vfmsq_f64(h, a, m);
    const float64x2_t q   = vrndnq_f64(vmulq_f64(h, c->pinv));
    const float64x2_t r   = vsubq_f64(vfmsq_f64(h, q, c->pd), nl);
    const float64x2_t y   = vsubq_f64(x, r);
    const float64x2_t k   = vrndnq_f64(vmulq_f64(y, c->pinv));
    return vfmsq_f64(y, k, c->pd);
}
#endif
#endif

/* d[j] = d[j] - m*a[j] mod p for 0 <= j < n */
static inline void ff64_submul_dense(
        int64_t *d,
        const cf64_t * const a,
        const len_t n,
        const uint64_t m,
        const ff64_t * const f
        )
{
    len_t j = 0;
#ifdef FF64_VL
    const ff64_vc_t c   = ff64_vc_init(f);
    const ff64_m mv     = ff64_m_set1(m, f);
    for (; j + FF64_VL <= n; j += FF64_VL) {
        ff64_v_store(d+j, ff64_v_submul(ff64_v_load(d+j),
                    ff64_v_load(a+j), mv, &c));
    }
#endif
    if (j < n) {
        const uint64_t mm = ff64_mont(m, f);
        for (; j < n; ++j) {
            d[j] = ff64_submul(d[j], a[j], mm, f);
        }
    }
}

/* dr[ds[j]] = dr[ds[j]] - m*a[j] mod p for 0 <= j < n */
static inline void ff64_submul_sparse(
        int64_t *dr,
        const hm_t * const ds,
        const cf64_t * const a,
        const len_t n,
        const uint64_t m,
        const ff64_t * const f
        )
{
    len_t j = 0;
#ifdef FF64_VL
    const ff64_vc_t c   = ff64_vc_init(f);
    const ff64_m mv     = ff64_m_set1(m, f);
    for (; j + FF64_VL <= n; j += FF64_VL) {
        ff64_v_scatter(dr, ds+j, ff64_v_submul(ff64_v_gather(dr, ds+j),
                    ff64_v_load(a+j), mv, &c));
    }
#endif
    if (j < n) {
        const uint64_t mm = ff64_mont(m, f);
        for (; j < n; ++j) {
            dr[ds[j]] = ff64_submul(dr[ds[j]], a[j], mm, f);
        }
    }
}

static inline cf64_t *normalize_dense_matrix_row_ff_64(
        cf64_t *row,
        const hm_t len,
        const uint64_t fc
        )
{
    len_t i;

    const ff64_t f      = ff64_init(fc);
    const uint64_t inv  = ff64_mont(mod_p_inverse_64((int64_t)row[0],
                (int64_t)fc), &f);

    for (i = 1; i < len; ++i) {
        row[i]  = ff64_mul(row[i], inv, &f);
    }
    row[0]  = 1;

    return row;
}

static inline cf64_t *normalize_sparse_matrix_row_ff_64(
        cf64_t *row,
        const len_t os,
        const len_t len,
        const uint64_t fc
        )
{
    (void)os;
    return normalize_dense_matrix_row_ff_64(row, len, fc);
}

/* Reduces the dense row dr by the pivots pivs[dpiv],...,pivs[ncols-1].
 * If rba is not NULL, the pivots used from the left part of the matrix
 * are marked in rba. Returns the number of nonzero entries left. */
static len_t reduce_dense_row_by_known_pivots_ff_64(
        int64_t *dr,
        const mat_t * const mat,
        hm_t * const *pivs,
        const hi_t dpiv,
        rba_t *rba,
        uint64_t *nops,
        const ff64_t * const f
        )
{
    hi_t i;
    len_t k = 0;
    const len_t ncols = mat->nc;
    const len_t ncl   = mat->ncl;
    cf64_t * const * const mcf  = mat->cf_64;

    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            k++;
            continue;
        }
        if (rba != NULL && i < ncl) {
            rba[i/32] |= 1U << (i % 32);
        }
        const hm_t * const dts  = pivs[i];
        const len_t len         = dts[LENGTH];
        ff64_submul_sparse(dr, dts+OFFSET, mcf[dts[COEFFS]], len,
                (uint64_t)dr[i], f);
        dr[i] = 0;
        *nops += len;
    }
    return k;
}

/* stores the nonzero entries of dr in the columns ncl,...,ncols-1 as new
 * row of mat, k is an upper bound for their number */
static hm_t *dense_row_to_sparse_row_ff_64(
        int64_t *dr,
        mat_t *mat,
        const len_t k,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi
        )
{
    hi_t i;
    len_t j = 0;
    const len_t ncols = mat->nc;

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf64_t *cf  = (cf64_t *)malloc((uint64_t)(k) * sizeof(cf64_t));
    hm_t *rs    = row + OFFSET;
    for (i = mat->ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf64_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_64[tmp_pos] = cf;

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_ff_64(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        )
{
    (void)bs;
    uint64_t nops = 0;
    const ff64_t f  = ff64_init(st->fc);
    rba_t *rba      = tr > 0 ? mat->rba[tmp_pos] : NULL;

    const len_t k = reduce_dense_row_by_known_pivots_ff_64(
            dr, mat, pivs, dpiv, rba, &nops, &f);
    if (k == 0) {
        return NULL;
    }
    return dense_row_to_sparse_row_ff_64(dr, mat, k, tmp_pos, mh, bi);
}

static hm_t *trace_reduce_dense_row_by_known_pivots_sparse_ff_64(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        md_t *st
        )
{
    (void)bs;
    uint64_t nops = 0;
    const ff64_t f  = ff64_init(st->fc);

    const len_t k = reduce_dense_row_by_known_pivots_ff_64(
            dr, mat, pivs, dpiv, rba, &nops, &f);
    st->trace_nr_mult +=  nops / 1000.0;
    st->trace_nr_add  +=  nops / 1000.0;
    if (k == 0) {
        return NULL;
    }
    return dense_row_to_sparse_row_ff_64(dr, mat, k, tmp_pos, mh, bi);
}

/* Reduces the dense row dr by the new dense pivots pivs[*pc],...,
 * pivs[ncr-1], each one starting at its pivot column. Returns the
 * normalized remainder starting at its first nonzero entry whose column
 * is stored in *pc, NULL and *pc = -1 if the remainder is zero. */
static cf64_t *reduce_dense_row_by_dense_new_pivots_ff_64(
        int64_t *dr,
        len_t *pc,
        cf64_t * const * const pivs,
        const len_t ncr,
        const uint64_t fc
        )
{
    hi_t i;
    len_t np  = -1;
    const ff64_t f  = ff64_init(fc);

    for (i = *pc; i < ncr; ++i) {
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == (len_t)-1) {
                np  = i;
            }
            continue;
        }
        ff64_submul_dense(dr+i, pivs[i], ncr-i, (uint64_t)dr[i], &f);
    }
    if (np == (len_t)-1) {
        *pc = -1;
        return NULL;
    }

    cf64_t *row = (cf64_t *)malloc((uint64_t)(ncr-np) * sizeof(cf64_t));
    for (i = np; i < ncr; ++i) {
        row[i-np]  = (cf64_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_64(row, ncr-np, fc);
    }
    *pc = np;

    return row;
}

/* Reduces the dense row dr by the known sparse pivots pivs[*pc],...,
 * pivs[ncl-1] and afterwards by the new dense pivots dpivs, see
 * reduce_dense_row_by_dense_new_pivots_ff_64(). The column of the
 * remainder is returned in *pc relative to ncl. */
static cf64_t *reduce_dense_row_by_all_pivots_ff_64(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t * const * const pivs,
        cf64_t * const * const dpivs,
        const uint64_t fc
        )
{
    hi_t i;
    const len_t ncl   = mat->ncl;
    const len_t ncr   = mat->ncr;
    const ff64_t f    = ff64_init(fc);

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
        if (dr[i] == 0 || pivs[i] == NULL) {
            continue;
        }
        ff64_submul_sparse(dr, pivs[i]+OFFSET, bs->cf_64[pivs[i][COEFFS]],
                pivs[i][LENGTH], (uint64_t)dr[i], &f);
        dr[i] = 0;
    }
    /* step 2: reduce by new dense pivots */
    len_t np = 0;
    cf64_t *row = reduce_dense_row_by_dense_new_pivots_ff_64(
            dr+ncl, &np, dpivs, ncr, fc);
    *pc = np;
    return row;
}

/* see la_ff_common.c */
#define LA_WIDE 1
#define LA_CF_T cf64_t
#define LA_CF cf_64
#define LA_FF(name) name##_ff_64
/* the two draws are sequenced so that the multiplier does not depend on
 * the evaluation order chosen by the compiler */
static inline int64_t random_multiplier_ff_64(msolve_rng_t *rng, uint64_t fc)
{
    const uint64_t hi = (uint64_t)msolve_rng_rand(rng);
    const uint64_t lo = (uint64_t)msolve_rng_rand(rng);
    return (int64_t)(((hi << 31) ^ lo) % fc);
}
#define LA_RANDOM_MULTIPLIER(rng, fc) random_multiplier_ff_64((rng), (fc))

/* maximal number of rows reduced together by the known pivots, see
 * LA_FF(reduce_dense_rows_by_known_pivots) in la_ff_common.c */
#define LA_NR 8

/* 64 bit version of LA_FF(reduce_dense_rows_by_known_pivots) in
 * la_ff_common.c: the entries of the interleaved rows are in [0,p) */
static void reduce_dense_rows_by_known_pivots_ff_64(
        int64_t *dr,
        cf64_t **drs,
        rba_t * const * const rbas,
        const mat_t * const mat,
        hm_t * const * const pivs,
        const len_t sc,
        const uint64_t fc
        )
{
    len_t i, j, r;

    const len_t ncl     = mat->ncl;
    const len_t ncols   = mat->nc;
    const ff64_t f      = ff64_init(fc);
    cf64_t * const * const mcf = mat->cf_64;
#if defined FF64_DOUBLE_ROWS
    const ff64_vc_t c   = ff64_vc_init(&f);
    const ff64_dc_t dc  = ff64_dc_init(&f);
    const ff64_d zero   = ff64_d_set1(0.0);

    /* the rows are only changed in the columns >= sc, they are kept as
     * doubles during the reduction */
    for (i = sc; i < ncols; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        for (r = 0; r < LA_NR; r += FF64_VL) {
            ff64_d_store(d+r, ff64_d_from_int(ff64_v_load(d+r), &dc));
        }
    }
    for (i = sc; i < ncl; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        ff64_d mv[LA_NR/FF64_VL];
        int nz = 0;
        for (r = 0; r < LA_NR/FF64_VL; ++r) {
            mv[r] = ff64_d_load(d+FF64_VL*r);
            nz |= ff64_d_nonzero(mv[r]);
        }
        if (nz == 0 || pivs[i] == NULL) {
            continue;
        }
        if (rbas != NULL) {
            for (r = 0; r < LA_NR; ++r) {
                double m;
                memcpy(&m, d+r, sizeof(double));
                if (m != 0.0) {
                    rbas[r][i/32] |= 1U << (i % 32);
                }
            }
        }
        const hm_t * const piv  = pivs[i];
        const cf64_t * const cfs = mcf[piv[COEFFS]];
        const len_t len         = piv[LENGTH];
        const hm_t * const ds   = piv + OFFSET;
        for (j = 0; j < len; ++j) {
            int64_t *t        = dr + (uint64_t)LA_NR * ds[j];
            const ff64_d cfv  = ff64_d_set1((double)cfs[j]);
            for (r = 0; r < LA_NR/FF64_VL; ++r) {
                ff64_d_store(t+FF64_VL*r, ff64_d_submul(
                            ff64_d_load(t+FF64_VL*r), cfv, mv[r], &dc));
            }
        }
    }
    /* the columns of the known pivots are zero now */
    for (i = sc; i < ncl; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        for (r = 0; r < LA_NR; r += FF64_VL) {
            ff64_d_store(d+r, zero);
        }
    }
    for (i = ncl; i < ncols; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        for (r = 0; r < LA_NR; r += FF64_VL) {
            ff64_v_store(d+r, ff64_d_to_int(ff64_d_load(d+r), &c));
        }
    }
#else
#ifdef FF64_VL
    const ff64_vc_t c   = ff64_vc_init(&f);
#endif
    for (i = sc; i < ncl; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        int64_t nz = 0;
        for (r = 0; r < LA_NR; ++r) {
            nz |= d[r];
        }
        if (nz == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            for (r = 0; r < LA_NR; ++r) {
                d[r] = 0;
            }
            continue;
        }
        if (rbas != NULL) {
            for (r = 0; r < LA_NR; ++r) {
                if (d[r] != 0) {
                    rbas[r][i/32] |= 1U << (i % 32);
                }
            }
        }
        const hm_t * const piv  = pivs[i];
        const cf64_t * const cfs = mcf[piv[COEFFS]];
        const len_t len         = piv[LENGTH];
        const hm_t * const ds   = piv + OFFSET;
#ifdef FF64_VL
        ff64_m mv[LA_NR/FF64_VL];
        for (r = 0; r < LA_NR/FF64_VL; ++r) {
            mv[r] = ff64_m_from(ff64_v_load(d+FF64_VL*r), &c);
        }
        for (j = 0; j < len; ++j) {
            int64_t *t        = dr + (uint64_t)LA_NR * ds[j];
            const ff64_v cfv  = ff64_v_set1(cfs[j]);
            for (r = 0; r < LA_NR/FF64_VL; ++r) {
                ff64_v_store(t+FF64_VL*r, ff64_v_submul(
                            ff64_v_load(t+FF64_VL*r), cfv, mv[r], &c));
            }
        }
#else
        uint64_t mm[LA_NR];
        for (r = 0; r < LA_NR; ++r) {
            mm[r] = ff64_mont((uint64_t)d[r], &f);
        }
        for (j = 0; j < len; ++j) {
            int64_t *t  = dr + (uint64_t)LA_NR * ds[j];
            for (r = 0; r < LA_NR; ++r) {
                t[r] = ff64_submul(t[r], cfs[j], mm[r], &f);
            }
        }
#endif
        for (r = 0; r < LA_NR; ++r) {
            d[r] = 0;
        }
    }
#endif
    if (drs == NULL) {
        return;
    }
    const len_t ncr = ncols - ncl;
    for (r = 0; r < LA_NR; ++r) {
        cf64_t *row = NULL;
        for (i = ncl; i < ncols; ++i) {
            int64_t *d = dr + (uint64_t)LA_NR * i + r;
            if (*d == 0) {
                continue;
            }
            if (row == NULL) {
                row = (cf64_t *)calloc((uint64_t)ncr, sizeof(cf64_t));
            }
            row[i-ncl] = (cf64_t)*d;
            *d = 0;
        }
        drs[r] = row;
    }
}

/* number of rows reduced together in LA_FF(interreduce_dense_matrix) */
#define LA_NRB 2

/* 64 bit version of LA_FF(reduce_dense_rows_by_dense_block) in
 * la_ff_common.c: the coefficients of the rows are in [0,p) and stored
 * as 64 bit integers, so the rows are reduced in place and dr is not
 * used. */
static void reduce_dense_rows_by_dense_block_ff_64(
        cf64_t **rows,
        const len_t * const rc,
        len_t nr,
        const len_t * const bpc,
        const len_t nb,
        cf64_t * const * const dm,
        const len_t ncr,
        const uint64_t fc,
        int64_t *dr
        )
{
    (void)dr;
    len_t i, j, r;
    int64_t *rs[LA_NRB];

    /* skip the rows which are zero in all pivot columns of the block */
    len_t n = 0;
    for (r = 0; r < nr; ++r) {
        cf64_t * const row = rows[r] - rc[r];
        for (i = 0; i < nb; ++i) {
            if (row[bpc[i]] != 0) {
                break;
            }
        }
        if (i < nb) {
            rs[n++] = (int64_t *)row;
        }
    }
    if (n == 0) {
        return;
    }
    nr = n;
    const ff64_t f      = ff64_init(fc);
#ifdef FF64_VL
    const ff64_vc_t c   = ff64_vc_init(&f);
#endif

    for (i = 0; i < nb; ++i) {
        const len_t cl = bpc[i];
        uint64_t mul[LA_NRB] = {0};
        int64_t nz = 0;
        for (r = 0; r < nr; ++r) {
            mul[r]  = (uint64_t)rs[r][cl];
            nz      |= rs[r][cl];
        }
        if (nz == 0) {
            continue;
        }
        const cf64_t * const piv = dm[cl] - cl;
        j = cl;
#ifdef FF64_VL
        if (nr == LA_NRB) {
            ff64_m mv[LA_NRB];
            for (r = 0; r < LA_NRB; ++r) {
                mv[r] = ff64_m_set1(mul[r], &f);
            }
            for (; j + FF64_VL <= ncr; j += FF64_VL) {
                const ff64_v cfv = ff64_v_load(piv+j);
                for (r = 0; r < LA_NRB; ++r) {
                    ff64_v_store(rs[r]+j, ff64_v_submul(
                                ff64_v_load(rs[r]+j), cfv, mv[r], &c));
                }
            }
        }
#endif
        for (r = 0; r < nr; ++r) {
            if (mul[r] != 0) {
                ff64_submul_dense(rs[r]+j, piv+j, ncr-j, mul[r], &f);
            }
        }
    }
}

#include "la_ff_common.c"
