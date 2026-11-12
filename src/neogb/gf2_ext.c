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

/* Arithmetic in the binary extension fields
 *
 *   GF(16)  = GF(2)[t]/(t^4+t+1),
 *   GF(256) = GF(2)[t]/(t^8+t^4+t^3+t+1).
 *
 * An element is stored in one byte, bit i is the coefficient of t^i, so
 * the elements of GF(16) are the bytes < 16. Addition and subtraction
 * are both XOR. Multiplication is done via lookup tables:
 *
 *   mul[a][b] = a*b,
 *
 * so multiplying many elements by the same a only needs the 256 byte row
 * mul[a]. Since b -> a*b is GF(2)-linear,
 *
 *   a*b = lo[a][b & 15] ^ hi[a][b >> 4]
 *
 * with the 16 byte tables lo[a][n] = a*n and hi[a][n] = a*(n << 4).
 * These are the tables of the byte shuffles (pshufb on x86, tbl on
 * aarch64) the vectorized kernels use, each shuffle multiplies 16, 32 or
 * 64 bytes at once. For GF(16) all hi[a] are zero, the kernels skip
 * them. */

#if defined HAVE_AVX2 || defined HAVE_AVX512_BW || defined HAVE_SSSE3
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

/* number of rows reduced together by the known pivots, the entries of
 * the rows in one column fill one or two vector registers */
#if defined HAVE_AVX512_BW
#define GF_NR 64
#else
#define GF_NR 32
#endif

typedef struct {
    uint8_t mul[256][256];
    uint8_t lo[256][16];
    uint8_t hi[256][16];
    uint8_t inv[256];
} gf2x_t;

/* tables for GF(16) and GF(256) */
static gf2x_t gf2x_tables[2] __attribute__((aligned(64)));
static int gf2x_tables_ready = 0;

/* multiplication of a and b modulo poly of degree deg, only used to
 * fill the tables */
static uint8_t gf2x_mul_slow(
        uint8_t a,
        uint8_t b,
        const uint32_t poly,
        const uint32_t deg
        )
{
    uint32_t p = 0, x = a;
    while (b != 0) {
        if (b & 1) {
            p ^= x;
        }
        b >>= 1;
        x <<= 1;
        if (x & (1u << deg)) {
            x ^= poly;
        }
    }
    return (uint8_t)p;
}

static void gf2x_fill_table(
        gf2x_t *t,
        const uint32_t poly,
        const uint32_t deg
        )
{
    const uint32_t q = 1u << deg;
    uint32_t a, b;

    memset(t, 0, sizeof(gf2x_t));
    for (a = 0; a < q; ++a) {
        for (b = 0; b < q; ++b) {
            t->mul[a][b] = gf2x_mul_slow(a, b, poly, deg);
            if (t->mul[a][b] == 1) {
                t->inv[a] = (uint8_t)b;
            }
        }
        for (b = 0; b < 16; ++b) {
            t->lo[a][b] = t->mul[a][b];
            if ((b << 4) < q) {
                t->hi[a][b] = t->mul[a][b << 4];
            }
        }
    }
}

/* fills the tables once, the caller may run in parallel to other
 * computations, e.g. from Julia */
static void gf2x_init_tables(
        void
        )
{
    if (__atomic_load_n(&gf2x_tables_ready, __ATOMIC_ACQUIRE)) {
        return;
    }
#pragma omp critical(gf2x_init_tables)
    {
        if (!__atomic_load_n(&gf2x_tables_ready, __ATOMIC_RELAXED)) {
            gf2x_fill_table(gf2x_tables, 0x13, 4);
            gf2x_fill_table(gf2x_tables + 1, 0x11B, 8);
            __atomic_store_n(&gf2x_tables_ready, 1, __ATOMIC_RELEASE);
        }
    }
}

static inline const gf2x_t *gf2x_table(
        const uint64_t fc
        )
{
    return gf2x_tables + (fc == 256);
}

/* Vector helpers: GF_VEC(lo, hi, x) multiplies the bytes of x by the
 * element whose tables lo and hi are broadcast into the registers lo and
 * hi, GF_VEC16 does the same if x is in GF(16). */
#if defined HAVE_AVX512_BW
static inline __m512i gf2x_mul_512(
        const __m512i tl,
        const __m512i th,
        const __m512i x
        )
{
    const __m512i m = _mm512_set1_epi8(0x0F);
    return _mm512_xor_si512(
            _mm512_shuffle_epi8(tl, _mm512_and_si512(x, m)),
            _mm512_shuffle_epi8(th, _mm512_and_si512(_mm512_srli_epi16(x, 4), m)));
}
#endif
#if defined HAVE_AVX2
static inline __m256i gf2x_mul_256(
        const __m256i tl,
        const __m256i th,
        const __m256i x
        )
{
    const __m256i m = _mm256_set1_epi8(0x0F);
    return _mm256_xor_si256(
            _mm256_shuffle_epi8(tl, _mm256_and_si256(x, m)),
            _mm256_shuffle_epi8(th, _mm256_and_si256(_mm256_srli_epi16(x, 4), m)));
}
#endif
#if defined HAVE_SSSE3
static inline __m128i gf2x_mul_128(
        const __m128i tl,
        const __m128i th,
        const __m128i x
        )
{
    const __m128i m = _mm_set1_epi8(0x0F);
    return _mm_xor_si128(
            _mm_shuffle_epi8(tl, _mm_and_si128(x, m)),
            _mm_shuffle_epi8(th, _mm_and_si128(_mm_srli_epi16(x, 4), m)));
}
#elif defined __aarch64__
static inline uint8x16_t gf2x_mul_neon(
        const uint8x16_t tl,
        const uint8x16_t th,
        const uint8x16_t x
        )
{
    return veorq_u8(vqtbl1q_u8(tl, vandq_u8(x, vdupq_n_u8(0x0F))),
            vqtbl1q_u8(th, vshrq_n_u8(x, 4)));
}
#endif

/* dst[i] = c*src[i] (axpy == 0) resp. dst[i] ^= c*src[i] (axpy != 0)
 * for 0 <= i < n, dst and src may be equal */
static inline void gf2x_mul_region(
        uint8_t *dst,
        const uint8_t *src,
        const len_t n,
        const uint8_t c,
        const gf2x_t * const t,
        const int axpy
        )
{
    len_t i = 0;
#if defined HAVE_AVX512_BW
    {
        const __m512i tl = _mm512_broadcast_i32x4(
                _mm_load_si128((const __m128i *)t->lo[c]));
        const __m512i th = _mm512_broadcast_i32x4(
                _mm_load_si128((const __m128i *)t->hi[c]));
        for (; i + 64 <= n; i += 64) {
            __m512i p = gf2x_mul_512(tl, th, _mm512_loadu_si512(src+i));
            if (axpy) {
                p = _mm512_xor_si512(p, _mm512_loadu_si512(dst+i));
            }
            _mm512_storeu_si512(dst+i, p);
        }
    }
#endif
#if defined HAVE_AVX2
    {
        const __m256i tl = _mm256_broadcastsi128_si256(
                _mm_load_si128((const __m128i *)t->lo[c]));
        const __m256i th = _mm256_broadcastsi128_si256(
                _mm_load_si128((const __m128i *)t->hi[c]));
        for (; i + 32 <= n; i += 32) {
            __m256i p = gf2x_mul_256(tl, th,
                    _mm256_loadu_si256((const __m256i *)(src+i)));
            if (axpy) {
                p = _mm256_xor_si256(p,
                        _mm256_loadu_si256((const __m256i *)(dst+i)));
            }
            _mm256_storeu_si256((__m256i *)(dst+i), p);
        }
    }
#endif
#if defined HAVE_SSSE3
    {
        const __m128i tl = _mm_load_si128((const __m128i *)t->lo[c]);
        const __m128i th = _mm_load_si128((const __m128i *)t->hi[c]);
        for (; i + 16 <= n; i += 16) {
            __m128i p = gf2x_mul_128(tl, th,
                    _mm_loadu_si128((const __m128i *)(src+i)));
            if (axpy) {
                p = _mm_xor_si128(p, _mm_loadu_si128((const __m128i *)(dst+i)));
            }
            _mm_storeu_si128((__m128i *)(dst+i), p);
        }
    }
#elif defined __aarch64__
    {
        const uint8x16_t tl = vld1q_u8(t->lo[c]);
        const uint8x16_t th = vld1q_u8(t->hi[c]);
        for (; i + 16 <= n; i += 16) {
            uint8x16_t p = gf2x_mul_neon(tl, th, vld1q_u8(src+i));
            if (axpy) {
                p = veorq_u8(p, vld1q_u8(dst+i));
            }
            vst1q_u8(dst+i, p);
        }
    }
#endif
    const uint8_t * const mr = t->mul[c];
    if (axpy) {
        for (; i < n; ++i) {
            dst[i] ^= mr[src[i]];
        }
    } else {
        for (; i < n; ++i) {
            dst[i] = mr[src[i]];
        }
    }
}

/* returns the smallest index i with start <= i < n and dr[i] != 0, n if
 * there is none */
static inline len_t gf2x_next_nonzero(
        const uint8_t * const dr,
        len_t i,
        const len_t n
        )
{
#if defined HAVE_AVX512_BW
    for (; i + 64 <= n; i += 64) {
        const __m512i v = _mm512_loadu_si512(dr+i);
        const uint64_t m = _mm512_test_epi8_mask(v, v);
        if (m != 0) {
            return i + (len_t)__builtin_ctzll(m);
        }
    }
#elif defined HAVE_AVX2
    const __m256i z = _mm256_setzero_si256();
    for (; i + 32 <= n; i += 32) {
        const __m256i v = _mm256_loadu_si256((const __m256i *)(dr+i));
        const uint32_t m = ~(uint32_t)_mm256_movemask_epi8(
                _mm256_cmpeq_epi8(v, z));
        if (m != 0) {
            return i + (len_t)__builtin_ctz(m);
        }
    }
#elif defined HAVE_SSSE3
    const __m128i z = _mm_setzero_si128();
    for (; i + 16 <= n; i += 16) {
        const __m128i v = _mm_loadu_si128((const __m128i *)(dr+i));
        const uint32_t m = 0xFFFFu & ~(uint32_t)_mm_movemask_epi8(
                _mm_cmpeq_epi8(v, z));
        if (m != 0) {
            return i + (len_t)__builtin_ctz(m);
        }
    }
#elif defined __aarch64__
    for (; i + 16 <= n; i += 16) {
        if (vmaxvq_u8(vld1q_u8(dr+i)) != 0) {
            break;
        }
    }
#else
    /* eight bytes at once */
    for (; i + 8 <= n; i += 8) {
        uint64_t v;
        memcpy(&v, dr+i, sizeof(uint64_t));
        if (v != 0) {
            break;
        }
    }
#endif
    for (; i < n; ++i) {
        if (dr[i] != 0) {
            return i;
        }
    }
    return n;
}

/* returns 1 if the GF_NR bytes at d are all zero */
static inline int gf2x_lanes_zero(
        const uint8_t * const d
        )
{
#if defined HAVE_AVX512_BW
    const __m512i v = _mm512_load_si512(d);
    return _mm512_test_epi8_mask(v, v) == 0;
#elif defined HAVE_AVX2
    const __m256i v = _mm256_load_si256((const __m256i *)d);
    return _mm256_testz_si256(v, v);
#else
    uint64_t v[GF_NR/8];
    uint64_t o = 0;
    memcpy(v, d, GF_NR);
    for (len_t r = 0; r < GF_NR/8; ++r) {
        o |= v[r];
    }
    return o == 0;
#endif
}

/* d[r] ^= m[r]*c for 0 <= r < GF_NR, i.e. adds the multiple of a pivot
 * with coefficient c in some column to the interleaved rows. m holds the
 * multipliers of the rows, it is split in nibbles: ml = m & 15 and
 * mh = m >> 4, both as GF_NR bytes, for GF(16) mh is not used. */
static inline void gf2x_lanes_axpy(
        uint8_t * const d,
        const uint8_t * const ml,
        const uint8_t * const mh,
        const uint8_t c,
        const gf2x_t * const t,
        const int q16
        )
{
#if defined HAVE_AVX512_BW
    const __m512i tl = _mm512_broadcast_i32x4(
            _mm_load_si128((const __m128i *)t->lo[c]));
    __m512i p = _mm512_shuffle_epi8(tl, _mm512_load_si512(ml));
    if (!q16) {
        const __m512i th = _mm512_broadcast_i32x4(
                _mm_load_si128((const __m128i *)t->hi[c]));
        p = _mm512_xor_si512(p, _mm512_shuffle_epi8(th, _mm512_load_si512(mh)));
    }
    _mm512_store_si512(d, _mm512_xor_si512(p, _mm512_load_si512(d)));
#elif defined HAVE_AVX2
    const __m256i tl = _mm256_broadcastsi128_si256(
            _mm_load_si128((const __m128i *)t->lo[c]));
    __m256i p = _mm256_shuffle_epi8(tl, _mm256_load_si256((const __m256i *)ml));
    if (!q16) {
        const __m256i th = _mm256_broadcastsi128_si256(
                _mm_load_si128((const __m128i *)t->hi[c]));
        p = _mm256_xor_si256(p, _mm256_shuffle_epi8(th,
                    _mm256_load_si256((const __m256i *)mh)));
    }
    _mm256_store_si256((__m256i *)d,
            _mm256_xor_si256(p, _mm256_load_si256((const __m256i *)d)));
#elif defined HAVE_SSSE3
    const __m128i tl = _mm_load_si128((const __m128i *)t->lo[c]);
    __m128i p0 = _mm_shuffle_epi8(tl, _mm_load_si128((const __m128i *)ml));
    __m128i p1 = _mm_shuffle_epi8(tl, _mm_load_si128((const __m128i *)(ml+16)));
    if (!q16) {
        const __m128i th = _mm_load_si128((const __m128i *)t->hi[c]);
        p0 = _mm_xor_si128(p0, _mm_shuffle_epi8(th,
                    _mm_load_si128((const __m128i *)mh)));
        p1 = _mm_xor_si128(p1, _mm_shuffle_epi8(th,
                    _mm_load_si128((const __m128i *)(mh+16))));
    }
    _mm_store_si128((__m128i *)d,
            _mm_xor_si128(p0, _mm_load_si128((const __m128i *)d)));
    _mm_store_si128((__m128i *)(d+16),
            _mm_xor_si128(p1, _mm_load_si128((const __m128i *)(d+16))));
#elif defined __aarch64__
    const uint8x16_t tl = vld1q_u8(t->lo[c]);
    uint8x16_t p0 = vqtbl1q_u8(tl, vld1q_u8(ml));
    uint8x16_t p1 = vqtbl1q_u8(tl, vld1q_u8(ml+16));
    if (!q16) {
        const uint8x16_t th = vld1q_u8(t->hi[c]);
        p0 = veorq_u8(p0, vqtbl1q_u8(th, vld1q_u8(mh)));
        p1 = veorq_u8(p1, vqtbl1q_u8(th, vld1q_u8(mh+16)));
    }
    vst1q_u8(d, veorq_u8(p0, vld1q_u8(d)));
    vst1q_u8(d+16, veorq_u8(p1, vld1q_u8(d+16)));
#else
    const uint8_t * const mr = t->mul[c];
    (void)q16;
    for (len_t r = 0; r < GF_NR; ++r) {
        d[r] ^= mr[ml[r] | (mh[r] << 4)];
    }
#endif
}

/* splits the GF_NR bytes m into ml = m & 15 and mh = m >> 4 */
static inline void gf2x_lanes_split(
        uint8_t * const ml,
        uint8_t * const mh,
        const uint8_t * const m
        )
{
    for (len_t r = 0; r < GF_NR; ++r) {
        ml[r] = m[r] & 0x0F;
        mh[r] = m[r] >> 4;
    }
}
