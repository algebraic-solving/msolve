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

/* Linear algebra over the integers for strong Groebner bases.
 *
 * The rows of a matrix generate a lattice. The reduction keeps this
 * lattice and computes for each column at most one pivot whose lead
 * coefficient is the gcd of what the rows reach in this column. If the
 * lead coefficient a of a pivot does not divide the lead coefficient b of
 * a row, both are replaced by the unimodular combination
 *
 *   pivot' = u * pivot + v * row,        u * a + v * b = g = gcd(a, b),
 *   row'   = b/g * pivot - a/g * row.
 *
 * pivot' has lead coefficient g, row' is zero in this column. This gcd
 * step adds no columns and no rows to the matrix. All other entries in a
 * pivot column are reduced into the symmetric residue system (-a/2, a/2]
 * modulo the lead coefficient a of the pivot. Lead coefficients are
 * positive. Content is never removed: over the integers x and 2x generate
 * different ideals.
 *
 * The reduced rows are unique: a lattice has only one echelon basis whose
 * rows are reduced modulo the pivots, apart from the reducer rows from the
 * basis which are kept unless a gcd step replaces them. So the result does
 * not depend on the order of the reductions nor on the number of threads.
 *
 * Coefficients: rows whose coefficients fit into 32 bits store them as
 * int32_t, a dense row is reduced in int64_t as long as a bound on its
 * entries stays below 2^62; these reductions use AVX512, AVX2 or NEON.
 * Otherwise, and for gcd steps, GMP integers are used.
 *
 * Parallelization: the rows are reduced in parallel, new pivots and the
 * pivots of gcd steps are set by compare and swap. A replaced pivot may
 * still be read by other threads, it is freed after the reduction.
 *
 * Coefficient growth: a new pivot is reduced by the pivots present when it
 * is created. A pivot with unreduced entries in the columns of later
 * pivots lets the entries of the rows it reduces grow fast. So each column
 * stores the generation in which its pivot was set, each new pivot the
 * generation in which its reduction started, and a pivot with an entry in
 * a column with a newer pivot is reduced again before it is used. The new
 * pivot of a gcd step is reduced right away: the entries of the row right
 * of its lead column are not reduced yet when the gcd step is done.
 *
 * The reducer rows from the basis must have positive lead coefficients.
 * A reducer row replaced by a gcd step stays in the basis, the new pivot
 * is returned as a new row. */

#include "data.h"
#include "../msolve/streams.h"

#if defined HAVE_AVX2
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

/* bound for the entries of dense rows in 64 bits */
#define ZZ_LIM ((uint64_t)1 << 62)

/* a pivot row: coefficients in c32 if all of them fit into 32 bits,
 * otherwise in cf; reducer rows from the basis borrow row and cf */
typedef struct zz_row_t zz_row_t;
struct zz_row_t
{
    hm_t *row;      /* column indices, with the usual row header */
    int32_t *c32;   /* coefficients if they fit into 32 bits, else NULL */
    mpz_t *cf;      /* coefficients otherwise, else NULL */
    uint32_t mx;    /* maximal absolute value of the coefficients in c32 */
    int8_t own;     /* 1 for a new row, 0 for a reducer row from the basis */
    uint64_t gen;   /* generation in which its reduction started */
};

/* a dense row of one thread: entries in d64 while a bound bd on their
 * absolute values stays below 2^62, otherwise (big = 1) in dz; the other
 * array is zero */
typedef struct zz_dense_t zz_dense_t;
struct zz_dense_t
{
    int64_t *d64;
    mpz_t *dz;
    mpz_t *dp;          /* new pivot of a gcd step, zero between two steps */
    uint64_t bd;
    int big;
    len_t lo;           /* entries before column lo are zero */
    mpz_t a, q, g, u, v, s, t;
    zz_row_t **ret;     /* replaced pivots, freed after the reduction */
    len_t nret;
    len_t sret;
    int64_t ngcd;       /* number of gcd steps */
    uint64_t *gen;      /* current generation, shared by all threads */
    uint64_t *cgen;     /* generations in which the pivots were set */
    uint64_t g0;        /* generation when the current reduction started */
    struct zz_dense_t *aux; /* for reducing pivots again, NULL if not */
    /* very sparse rows with 64 bit entries (sp = 1): the columns of
     * possibly nonzero entries are kept in the heap hp of size nh, mk marks
     * them, the columns which stay nonzero are collected in rem */
    int sp;
    hm_t *hp;
    len_t nh;
    uint8_t *mk;
    hm_t *rem;
    len_t nrem;
};

static inline zz_row_t *zz_load_piv(
        zz_row_t * const *pivs,
        const len_t c
        )
{
    return __atomic_load_n(&pivs[c], __ATOMIC_ACQUIRE);
}

static inline int zz_cas_piv(
        zz_row_t **pivs,
        const len_t c,
        zz_row_t *old,
        zz_row_t *new
        )
{
    return __atomic_compare_exchange_n(&pivs[c], &old, new, 0,
            __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE);
}

static inline int zz_fits32(
        const mpz_t x
        )
{
    return mpz_cmpabs_ui(x, INT32_MAX) <= 0;
}

static void zz_dense_init(
        zz_dense_t *d,
        const len_t nc
        )
{
    len_t i;

    d->d64 = (int64_t *)calloc((unsigned long)nc, sizeof(int64_t));
    d->dz  = (mpz_t *)malloc((unsigned long)nc * sizeof(mpz_t));
    d->dp  = (mpz_t *)malloc((unsigned long)nc * sizeof(mpz_t));
    for (i = 0; i < nc; ++i) {
        mpz_init(d->dz[i]);
        mpz_init(d->dp[i]);
    }
    mpz_inits(d->a, d->q, d->g, d->u, d->v, d->s, d->t, NULL);
    d->bd   = 0;
    d->big  = 0;
    d->lo   = 0;
    d->ret  = NULL;
    d->nret = 0;
    d->sret = 0;
    d->ngcd = 0;
    d->gen  = NULL;
    d->cgen = NULL;
    d->g0   = 0;
    d->aux  = NULL;
    d->sp   = 0;
    d->nh   = 0;
    d->nrem = 0;
    d->hp   = (hm_t *)malloc((unsigned long)nc * sizeof(hm_t));
    d->mk   = (uint8_t *)calloc((unsigned long)nc, sizeof(uint8_t));
    d->rem  = (hm_t *)malloc((unsigned long)nc * sizeof(hm_t));
}

static void zz_row_free(
        zz_row_t *r
        )
{
    len_t j;

    free(r->c32);
    if (r->own != 0) {
        if (r->cf != NULL) {
            for (j = 0; j < r->row[LENGTH]; ++j) {
                mpz_clear(r->cf[j]);
            }
            free(r->cf);
        }
        free(r->row);
        free(r);
    }
}

/* frees the dense row and the pivots it replaced */
static void zz_dense_free(
        zz_dense_t *d,
        const len_t nc
        )
{
    len_t i;

    for (i = 0; i < d->nret; ++i) {
        zz_row_free(d->ret[i]);
    }
    free(d->ret);
    for (i = 0; i < nc; ++i) {
        mpz_clear(d->dz[i]);
        mpz_clear(d->dp[i]);
    }
    free(d->dz);
    free(d->dp);
    free(d->d64);
    free(d->hp);
    free(d->mk);
    free(d->rem);
    mpz_clears(d->a, d->q, d->g, d->u, d->v, d->s, d->t, NULL);
}

/* adds the column c to the heap of a sparse dense row */
static inline void zz_heap_push(
        zz_dense_t *d,
        const hm_t c
        )
{
    if (d->mk[c] != 0) {
        return;
    }
    d->mk[c] = 1;
    hm_t *hp = d->hp;
    len_t k  = d->nh++;
    while (k > 0 && hp[(k-1)/2] > c) {
        hp[k] = hp[(k-1)/2];
        k     = (k-1)/2;
    }
    hp[k] = c;
}

/* removes and returns the smallest column of the heap */
static inline hm_t zz_heap_pop(
        zz_dense_t *d
        )
{
    hm_t *hp = d->hp;
    const hm_t c = hp[0];
    const hm_t x = hp[--d->nh];
    const len_t n = d->nh;
    len_t k = 0;
    for (;;) {
        len_t m = 2*k+1;
        if (m >= n) {
            break;
        }
        if (m+1 < n && hp[m+1] < hp[m]) {
            m++;
        }
        if (hp[m] >= x) {
            break;
        }
        hp[k] = hp[m];
        k     = m;
    }
    hp[k]    = x;
    d->mk[c] = 0;
    return c;
}

/* back to scanning all columns of the dense row */
static inline void zz_leave_sparse(
        zz_dense_t *d
        )
{
    while (d->nh > 0) {
        d->mk[d->hp[--d->nh]] = 0;
    }
    d->sp = 0;
}

static void zz_retire(
        zz_dense_t *d,
        zz_row_t *p
        )
{
    if (p->own == 0) {
        return;
    }
    if (d->nret == d->sret) {
        d->sret = d->sret == 0 ? 16 : 2 * d->sret;
        d->ret  = realloc(d->ret, (unsigned long)d->sret * sizeof(zz_row_t *));
    }
    d->ret[d->nret++] = p;
}

/* sets the 32 bit coefficients of r if all coefficients in cf fit */
static void zz_set_c32(
        zz_row_t *r,
        mpz_t *cf
        )
{
    len_t j;

    const len_t len = r->row[LENGTH];
    for (j = 0; j < len; ++j) {
        if (!zz_fits32(cf[j])) {
            r->c32 = NULL;
            return;
        }
    }
    r->c32 = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));
    r->mx  = 0;
    for (j = 0; j < len; ++j) {
        r->c32[j] = (int32_t)mpz_get_si(cf[j]);
        const uint32_t a = (uint32_t)abs(r->c32[j]);
        r->mx = a > r->mx ? a : r->mx;
    }
}

/* dr -= q * row for |q| < 2^31, the caller checks that the result fits */
static inline void zz_submul_64(
        int64_t *dr,
        const int64_t q,
        const hm_t * const row,
        const int32_t * const cf
        )
{
    len_t j = 0;

    const len_t len       = row[LENGTH];
    const hm_t * const ds = row + OFFSET;
#if defined HAVE_AVX512_F
    const __m512i qv = _mm512_set1_epi64(q);
    for (; j + 8 <= len; j += 8) {
        const __m256i iv = _mm256_loadu_si256((const __m256i *)(ds+j));
        const __m512i cv = _mm512_cvtepi32_epi64(
                _mm256_loadu_si256((const __m256i *)(cf+j)));
        __m512i drv = _mm512_i32gather_epi64(iv, dr, 8);
        drv = _mm512_sub_epi64(drv, _mm512_mul_epi32(qv, cv));
        _mm512_i32scatter_epi64(dr, iv, drv, 8);
    }
#elif defined HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    const __m256i qv = _mm256_set1_epi64x(q);
    for (; j + 4 <= len; j += 4) {
        const __m128i iv = _mm_loadu_si128((const __m128i *)(ds+j));
        const __m256i cv = _mm256_cvtepi32_epi64(
                _mm_loadu_si128((const __m128i *)(cf+j)));
        __m256i drv = _mm256_i32gather_epi64((const long long *)dr, iv, 8);
        drv = _mm256_sub_epi64(drv, _mm256_mul_epi32(qv, cv));
        _mm256_store_si256((__m256i *)res, drv);
        dr[ds[j]]   = res[0];
        dr[ds[j+1]] = res[1];
        dr[ds[j+2]] = res[2];
        dr[ds[j+3]] = res[3];
    }
#elif defined __aarch64__
    const int32x2_t qv = vdup_n_s32((int32_t)q);
    for (; j + 4 <= len; j += 4) {
        const int32x4_t cv = vld1q_s32(cf+j);
        int64x2_t lv = vcombine_s64(vcreate_s64((uint64_t)dr[ds[j]]),
                vcreate_s64((uint64_t)dr[ds[j+1]]));
        int64x2_t hv = vcombine_s64(vcreate_s64((uint64_t)dr[ds[j+2]]),
                vcreate_s64((uint64_t)dr[ds[j+3]]));
        lv = vmlsl_s32(lv, vget_low_s32(cv), qv);
        hv = vmlsl_s32(hv, vget_high_s32(cv), qv);
        dr[ds[j]]   = vgetq_lane_s64(lv, 0);
        dr[ds[j+1]] = vgetq_lane_s64(lv, 1);
        dr[ds[j+2]] = vgetq_lane_s64(hv, 0);
        dr[ds[j+3]] = vgetq_lane_s64(hv, 1);
    }
#endif
    for (; j < len; ++j) {
        dr[ds[j]] -= q * cf[j];
    }
}

/* dz += q * p (sub = 0) resp. dz -= q * p (sub = 1) */
static inline void zz_addmul_mpz(
        mpz_t *dz,
        const mpz_t q,
        const zz_row_t * const p,
        const int sub
        )
{
    len_t j;

    const len_t len       = p->row[LENGTH];
    const hm_t * const ds = p->row + OFFSET;
    if (p->c32 != NULL) {
        for (j = 0; j < len; ++j) {
            const int32_t c = p->c32[j];
            if ((c > 0) != (sub != 0)) {
                mpz_addmul_ui(dz[ds[j]], q, (unsigned long)abs(c));
            } else {
                mpz_submul_ui(dz[ds[j]], q, (unsigned long)abs(c));
            }
        }
    } else {
        if (sub != 0) {
            for (j = 0; j < len; ++j) {
                mpz_submul(dz[ds[j]], q, p->cf[j]);
            }
        } else {
            for (j = 0; j < len; ++j) {
                mpz_addmul(dz[ds[j]], q, p->cf[j]);
            }
        }
    }
}

static inline void zz_lead_coefficient(
        mpz_t a,
        const zz_row_t * const p
        )
{
    if (p->c32 != NULL) {
        mpz_set_si(a, p->c32[0]);
    } else {
        mpz_set(a, p->cf[0]);
    }
}

/* q such that b - q * a lies in (-a/2, a/2], a > 0, r is temporary */
static inline void zz_quotient(
        mpz_t q,
        mpz_t r,
        const mpz_t b,
        const mpz_t a
        )
{
    mpz_fdiv_qr(q, r, b, a);
    mpz_mul_2exp(r, r, 1);
    if (mpz_cmp(r, a) > 0) {
        mpz_add_ui(q, q, 1);
    }
}

static inline int64_t zz_quotient_64(
        const int64_t b,
        const int64_t a
        )
{
    int64_t q = b / a;
    int64_t r = b % a;
    if (r < 0) {
        r += a;
        q--;
    }
    return 2 * r > a ? q + 1 : q;
}

/* checks that the entries of d stay below 2^62 when subtracting q times a
 * row whose coefficients are at most mx in absolute value, recomputes
 * the bound from the entries if needed; updates the bound */
static inline int zz_room(
        zz_dense_t *d,
        const int64_t q,
        const uint32_t mx,
        const len_t nc
        )
{
    len_t j;

    const uint64_t aq = q < 0 ? -(uint64_t)q : (uint64_t)q;
    if (aq > INT32_MAX) {
        return 0;
    }
    const uint64_t add = aq * mx;
    if (d->bd + add < ZZ_LIM) {
        d->bd += add;
        return 1;
    }
    uint64_t m = 0;
    for (j = d->lo; j < nc; ++j) {
        const uint64_t a = d->d64[j] < 0 ? -(uint64_t)d->d64[j]
            : (uint64_t)d->d64[j];
        m = a > m ? a : m;
    }
    d->bd = m + add;
    return d->bd < ZZ_LIM;
}

static void zz_promote(
        zz_dense_t *d,
        const len_t nc
        )
{
    len_t j;

    zz_leave_sparse(d);
    for (j = d->lo; j < nc; ++j) {
        if (d->d64[j] != 0) {
            mpz_set_si(d->dz[j], d->d64[j]);
            d->d64[j] = 0;
        }
    }
    d->big = 1;
}

/* loads a row with coefficients cf (if c32 is NULL) or c32 into the zero
 * dense row d */
static void zz_dense_load(
        zz_dense_t *d,
        const hm_t * const row,
        const int32_t * const c32,
        mpz_t * const cf
        )
{
    len_t j;

    const len_t len       = row[LENGTH];
    const hm_t * const ds = row + OFFSET;
    d->lo = ds[0];
    d->sp = 0;
    if (c32 != NULL) {
        uint64_t m = 0;
        for (j = 0; j < len; ++j) {
            d->d64[ds[j]] = c32[j];
            const uint64_t a = (uint64_t)abs(c32[j]);
            m = a > m ? a : m;
        }
        d->bd  = m;
        d->big = 0;
    } else {
        for (j = 0; j < len; ++j) {
            mpz_set(d->dz[ds[j]], cf[j]);
        }
        d->big = 1;
    }
}

/* loads a row with coefficients from the basis, as 64 bit entries if all
 * of them fit into 32 bits */
static void zz_dense_load_basis_row(
        zz_dense_t *d,
        const hm_t * const row,
        mpz_t * const cf
        )
{
    len_t j;

    const len_t len = row[LENGTH];
    for (j = 0; j < len; ++j) {
        if (!zz_fits32(cf[j])) {
            zz_dense_load(d, row, NULL, cf);
            return;
        }
    }
    const hm_t * const ds = row + OFFSET;
    uint64_t m = 0;
    d->lo = ds[0];
    for (j = 0; j < len; ++j) {
        d->d64[ds[j]] = mpz_get_si(cf[j]);
        const uint64_t a = (uint64_t)labs(d->d64[ds[j]]);
        m = a > m ? a : m;
    }
    d->bd  = m;
    d->big = 0;
}

/* moves the entries from column lc on of d64 (big = 0) or dz (big = 1)
 * into a new row with positive lead coefficient, the dense entries are
 * zero afterwards */
static zz_row_t *zz_new_row(
        int64_t *d64,
        mpz_t *dz,
        const int big,
        const len_t lc,
        const len_t nc
        )
{
    len_t i, len = 0;
    int fits = 1;

    for (i = lc; i < nc; ++i) {
        if (big != 0) {
            if (mpz_sgn(dz[i]) != 0) {
                len++;
                fits &= zz_fits32(dz[i]);
            }
        } else {
            if (d64[i] != 0) {
                len++;
                fits &= d64[i] >= -INT32_MAX && d64[i] <= INT32_MAX;
            }
        }
    }
    zz_row_t *r = (zz_row_t *)malloc(sizeof(zz_row_t));
    r->row = (hm_t *)malloc((unsigned long)(len+OFFSET) * sizeof(hm_t));
    r->own = 1;
    r->gen = 0;
    r->mx  = 0;
    r->c32 = NULL;
    r->cf  = NULL;
    const int neg = big != 0 ? mpz_sgn(dz[lc]) < 0 : d64[lc] < 0;
    if (fits != 0) {
        r->c32 = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));
    } else {
        r->cf  = (mpz_t *)malloc((unsigned long)len * sizeof(mpz_t));
    }
    len = 0;
    for (i = lc; i < nc; ++i) {
        if (big != 0) {
            if (mpz_sgn(dz[i]) == 0) {
                continue;
            }
            if (neg != 0) {
                mpz_neg(dz[i], dz[i]);
            }
            if (fits != 0) {
                r->c32[len] = (int32_t)mpz_get_si(dz[i]);
                mpz_set_ui(dz[i], 0);
            } else {
                mpz_init(r->cf[len]);
                mpz_swap(r->cf[len], dz[i]);
            }
        } else {
            if (d64[i] == 0) {
                continue;
            }
            const int64_t e = neg != 0 ? -d64[i] : d64[i];
            if (fits != 0) {
                r->c32[len] = (int32_t)e;
            } else {
                mpz_init_set_si(r->cf[len], e);
            }
            d64[i] = 0;
        }
        if (fits != 0) {
            const uint32_t a = (uint32_t)abs(r->c32[len]);
            r->mx = a > r->mx ? a : r->mx;
        }
        r->row[OFFSET+len] = i;
        len++;
    }
    r->row[PRELOOP] = len % UNROLL;
    r->row[LENGTH]  = len;
    return r;
}

/* the new row of a sparse dense row, its entries are in the columns rem,
 * the dense entries are zero afterwards */
static zz_row_t *zz_new_row_sparse(
        zz_dense_t *d
        )
{
    len_t i;
    int fits = 1;

    int64_t *d64 = d->d64;
    const len_t len = d->nrem;
    const hm_t * const rem = d->rem;
    for (i = 0; i < len; ++i) {
        fits &= d64[rem[i]] >= -INT32_MAX && d64[rem[i]] <= INT32_MAX;
    }
    zz_row_t *r = (zz_row_t *)malloc(sizeof(zz_row_t));
    r->row = (hm_t *)malloc((unsigned long)(len+OFFSET) * sizeof(hm_t));
    r->own = 1;
    r->gen = 0;
    r->mx  = 0;
    r->c32 = NULL;
    r->cf  = NULL;
    const int neg = d64[rem[0]] < 0;
    if (fits != 0) {
        r->c32 = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));
    } else {
        r->cf  = (mpz_t *)malloc((unsigned long)len * sizeof(mpz_t));
    }
    for (i = 0; i < len; ++i) {
        const int64_t e = neg != 0 ? -d64[rem[i]] : d64[rem[i]];
        if (fits != 0) {
            r->c32[i] = (int32_t)e;
            const uint32_t a = (uint32_t)abs(r->c32[i]);
            r->mx = a > r->mx ? a : r->mx;
        } else {
            mpz_init_set_si(r->cf[i], e);
        }
        d64[rem[i]] = 0;
        r->row[OFFSET+i] = rem[i];
    }
    r->row[PRELOOP] = len % UNROLL;
    r->row[LENGTH]  = len;
    d->nrem = 0;
    d->sp   = 0;
    return r;
}

/* the pivot of column c was set, rows with entries in c are reduced again
 * before they are used as pivots */
static inline void zz_new_generation(
        zz_dense_t *d,
        const len_t c
        )
{
    const uint64_t g = __atomic_add_fetch(d->gen, 1, __ATOMIC_ACQ_REL);
    __atomic_store_n(&d->cgen[c], g, __ATOMIC_RELEASE);
}

/* 1 if the pivot p has an entry in a column whose pivot was set after the
 * reduction of p started */
static inline int zz_stale(
        const zz_dense_t * const d,
        const zz_row_t * const p
        )
{
    len_t j;

    const len_t len       = p->row[LENGTH];
    const hm_t * const ds = p->row + OFFSET;
    for (j = 1; j < len; ++j) {
        if (__atomic_load_n(&d->cgen[ds[j]], __ATOMIC_ACQUIRE) > p->gen) {
            return 1;
        }
    }
    return 0;
}

static int64_t zz_reduce_dense_row(zz_dense_t *d, zz_row_t **pivs,
        const len_t sc, int64_t lc, const len_t nc);

/* reduces the new pivot p of column c again and replaces it; returns the
 * new pivot, or NULL if another thread changed the pivot of c */
static zz_row_t *zz_refresh_pivot(
        zz_dense_t *d,
        zz_row_t **pivs,
        zz_row_t *p,
        const len_t c,
        const len_t nc
        )
{
    zz_dense_t *x = d->aux;
    const uint64_t g = __atomic_load_n(d->gen, __ATOMIC_ACQUIRE);
    zz_dense_load(x, p->row, p->c32, p->cf);
    zz_reduce_dense_row(x, pivs, c+1, c, nc);
    zz_row_t *r = zz_new_row(x->d64, x->dz, x->big, c, nc);
    x->big = 0;
    r->gen = g;
    if (zz_cas_piv(pivs, c, p, r)) {
        zz_retire(d, p);
        return r;
    }
    zz_row_free(r);
    return NULL;
}

/* gcd step of the dense row d (big) whose lead coefficient in column c is
 * not divisible by the lead coefficient d->a of the pivot p: tries to
 * replace p by u * p + v * row, on success d becomes b/g * p - a/g * row,
 * which is zero in column c; returns 0 if another thread changed the
 * pivot of c in between */
static int zz_gcd_step(
        zz_dense_t *d,
        zz_row_t **pivs,
        zz_row_t *p,
        const len_t c,
        const len_t nc
        )
{
    len_t j;

    mpz_t *dz = d->dz;
    mpz_gcdext(d->g, d->u, d->v, d->a, dz[c]);
    for (j = c; j < nc; ++j) {
        if (mpz_sgn(dz[j]) != 0) {
            mpz_mul(d->dp[j], d->v, dz[j]);
        }
    }
    zz_addmul_mpz(d->dp, d->u, p, 0);
    zz_row_t *np = zz_new_row(NULL, d->dp, 1, c, nc);
    np->gen      = d->g0;
    /* the entries of the row right of c are not reduced yet, so the new
     * pivot is reduced there before it is used */
    zz_dense_t *x = d->aux;
    if (x != NULL) {
        zz_dense_load(x, np->row, np->c32, np->cf);
        zz_row_free(np);
        zz_reduce_dense_row(x, pivs, c+1, c, nc);
        np      = zz_new_row(x->d64, x->dz, x->big, c, nc);
        x->big  = 0;
        np->gen = x->g0;
    }
    if (!zz_cas_piv(pivs, c, p, np)) {
        zz_row_free(np);
        return 0;
    }
    zz_new_generation(d, c);
    mpz_divexact(d->s, dz[c], d->g);
    mpz_divexact(d->t, d->a, d->g);
    mpz_neg(d->t, d->t);
    for (j = c; j < nc; ++j) {
        if (mpz_sgn(dz[j]) != 0) {
            mpz_mul(dz[j], dz[j], d->t);
        }
    }
    zz_addmul_mpz(dz, d->s, p, 0);
    zz_retire(d, p);
    d->ngcd++;
    return 1;
}

/* reduces the dense row d from column sc on by the pivots. Without a lead
 * column (lc = -1) its lead coefficient is eliminated by an exact multiple
 * of a pivot or by a gcd step until the row is zero or its lead column has
 * no pivot. The other entries are reduced into the residue systems of the
 * pivots. Returns the lead column, -1 if the row is zero. */
static int64_t zz_reduce_dense_row(
        zz_dense_t *d,
        zz_row_t **pivs,
        const len_t sc,
        int64_t lc,
        const len_t nc
        )
{
    len_t i;

    d->g0   = __atomic_load_n(d->gen, __ATOMIC_ACQUIRE);
    d->nrem = 0;
    /* a sparse row takes its columns from the heap, it may become dense */
    i = sc;
    for (;;) {
        if (d->sp != 0) {
            if (d->nh == 0) {
                break;
            }
            i = zz_heap_pop(d);
        } else if (i >= nc) {
            break;
        }
        if (d->big == 0 ? d->d64[i] == 0 : mpz_sgn(d->dz[i]) == 0) {
            goto next;
        }
        zz_row_t *p = zz_load_piv(pivs, i);
        if (p == NULL) {
            if (lc == -1) {
                lc = i;
            }
            goto next;
        }
        if (d->aux != NULL && p->own != 0 && zz_stale(d, p)) {
            p = zz_refresh_pivot(d, pivs, p, i, nc);
            if (p == NULL) {
                /* the pivot changed, redo this column */
                goto redo;
            }
        }
        if (d->big == 0) {
            if (p->c32 != NULL) {
                const int64_t a = p->c32[0];
                const int64_t b = d->d64[i];
                int64_t q;
                if (lc == -1) {
                    q = b % a == 0 ? b / a : 0;
                } else {
                    q = zz_quotient_64(b, a);
                    if (q == 0) {
                        goto next;
                    }
                }
                if (q != 0 && zz_room(d, q, p->mx, nc)) {
                    zz_submul_64(d->d64, q, p->row, p->c32);
                    if (d->sp != 0) {
                        const len_t len = p->row[LENGTH];
                        const hm_t * const ds = p->row + OFFSET;
                        for (len_t j = 1; j < len; ++j) {
                            zz_heap_push(d, ds[j]);
                        }
                    }
                    goto next;
                }
            }
            zz_promote(d, nc);
        }
        zz_lead_coefficient(d->a, p);
        if (lc == -1) {
            if (mpz_divisible_p(d->dz[i], d->a) != 0) {
                mpz_divexact(d->q, d->dz[i], d->a);
                zz_addmul_mpz(d->dz, d->q, p, 1);
            } else {
                if (!zz_gcd_step(d, pivs, p, i, nc)) {
                    /* the pivot changed, redo this column */
                    goto redo;
                }
            }
        } else {
            zz_quotient(d->q, d->t, d->dz[i], d->a);
            if (mpz_sgn(d->q) != 0) {
                zz_addmul_mpz(d->dz, d->q, p, 1);
            }
        }
next:
        if (d->sp != 0) {
            if (d->d64[i] != 0) {
                d->rem[d->nrem++] = i;
            }
        } else {
            i++;
        }
        continue;
redo:
        if (d->sp != 0) {
            zz_heap_push(d, i);
        }
    }
    return lc;
}

/* reduces the dense row d by the pivots and makes it the pivot of its lead
 * column; returns 1 for a new pivot, 0 if the row reduces to zero */
static int zz_insert_dense_row(
        zz_dense_t *d,
        zz_row_t **pivs,
        const len_t nc
        )
{
    len_t sc = d->lo;

    for (;;) {
        const int64_t lc = zz_reduce_dense_row(d, pivs, sc, -1, nc);
        if (lc == -1) {
            d->big = 0;
            d->sp  = 0;
            return 0;
        }
        zz_row_t *r = d->sp != 0 ? zz_new_row_sparse(d)
            : zz_new_row(d->d64, d->dz, d->big, lc, nc);
        d->big = 0;
        r->gen = d->g0;
        if (zz_cas_piv(pivs, lc, NULL, r)) {
            zz_new_generation(d, lc);
            return 1;
        }
        /* another thread set this pivot in between */
        zz_dense_load(d, r->row, r->c32, r->cf);
        zz_row_free(r);
        sc = lc;
    }
}

/* reduces the entries of each new pivot in the other pivot columns into
 * their residue systems, from the last column to the first one */
static void zz_interreduce_new_pivots(
        zz_dense_t *d,
        zz_row_t **pivs,
        const len_t nc
        )
{
    len_t k;

    /* the pivots right of k are final, no need to reduce them again */
    zz_dense_t *aux = d->aux;
    d->aux = NULL;
    for (k = nc; k-- > 0; ) {
        zz_row_t *p = pivs[k];
        if (p == NULL || p->own == 0) {
            continue;
        }
        zz_dense_load(d, p->row, p->c32, p->cf);
        zz_row_free(p);
        zz_reduce_dense_row(d, pivs, k+1, k, nc);
        pivs[k] = zz_new_row(d->d64, d->dz, d->big, k, nc);
        d->big  = 0;
    }
    d->aux = aux;
}

/* hands the new pivots over to mat->tr and mat->cf_qq, the coefficients
 * of a row are stored at the position of its lead column, the rows are
 * ordered by increasing (inc != 0) or decreasing lead column */
static len_t zz_collect_new_pivots(
        mat_t *mat,
        zz_row_t **pivs,
        const len_t nc,
        const int inc
        )
{
    len_t i, j, k, np = 0;

    for (i = 0; i < nc; ++i) {
        np += pivs[i] != NULL && pivs[i]->own != 0;
    }
    free(mat->cf_qq);
    mat->cf_qq = (mpz_t **)calloc((unsigned long)nc, sizeof(mpz_t *));
    mat->tr    = realloc(mat->tr, (unsigned long)np * sizeof(hm_t *));
    k = 0;
    for (i = 0; i < nc; ++i) {
        const len_t c = inc != 0 ? i : nc-1-i;
        zz_row_t *p   = pivs[c];
        if (p == NULL || p->own == 0) {
            continue;
        }
        if (p->c32 != NULL) {
            const len_t len = p->row[LENGTH];
            p->cf = (mpz_t *)malloc((unsigned long)len * sizeof(mpz_t));
            for (j = 0; j < len; ++j) {
                mpz_init_set_si(p->cf[j], p->c32[j]);
            }
            free(p->c32);
        }
        p->row[COEFFS] = c;
        mat->cf_qq[c]  = p->cf;
        mat->tr[k++]   = p->row;
        free(p);
        pivs[c] = NULL;
    }
    mat->np = mat->nr = mat->sz = np;
    return np;
}

/* reduces the rows rows[0], ..., rows[nr-1] with coefficients from the
 * basis by the pivots and among each other, in parallel; frees the rows;
 * returns the number of rows reduced to zero */
static len_t zz_reduce_rows(
        zz_row_t **pivs,
        hm_t **rows,
        const len_t nr,
        const bs_t * const bs,
        zz_dense_t *ds,
        const len_t nc,
        const uint64_t apl,
        const int nthrds
        )
{
    len_t i;
    len_t nzero = 0;

#pragma omp parallel for num_threads(nthrds) schedule(dynamic) \
    reduction(+:nzero)
    for (i = 0; i < nr; ++i) {
        zz_dense_t *d = ds + 2 * omp_get_thread_num();
        zz_dense_load_basis_row(d, rows[i], bs->cf_qq[rows[i][COEFFS]]);
        /* very sparse rows are reduced via a heap of their columns, as
         * over prime fields */
        const len_t len   = rows[i][LENGTH];
        const uint64_t rc = nc - rows[i][OFFSET];
        if (d->big == 0 && (uint64_t)len * 32 < rc
                && (uint64_t)len * apl * 64 < rc) {
            d->sp = 1;
            for (len_t j = 0; j < len; ++j) {
                zz_heap_push(d, rows[i][OFFSET+j]);
            }
        }
        free(rows[i]);
        rows[i] = NULL;
        if (zz_insert_dense_row(d, pivs, nc) == 0) {
            nzero++;
        }
    }
    return nzero;
}

/* two dense rows per thread, the second one for reducing pivots again,
 * sharing the generations */
static zz_dense_t *zz_init_dense_rows(
        uint64_t *gen,
        uint64_t *cgen,
        const len_t nc,
        const int nthrds
        )
{
    int i;

    zz_dense_t *ds = (zz_dense_t *)malloc(
            (unsigned long)(2 * nthrds) * sizeof(zz_dense_t));
    for (i = 0; i < 2 * nthrds; ++i) {
        zz_dense_init(ds + i, nc);
        ds[i].gen  = gen;
        ds[i].cgen = cgen;
    }
    for (i = 0; i < nthrds; ++i) {
        ds[2*i].aux = ds + 2*i + 1;
    }
    return ds;
}

static void zz_free_dense_rows(
        zz_dense_t *ds,
        md_t *st,
        const len_t nc
        )
{
    int i;

    for (i = 0; i < 2 * st->nthrds; ++i) {
        st->num_gcdred += ds[i].ngcd;
        zz_dense_free(ds + i, nc);
    }
    free(ds);
}

/* reduces the rows mat->tr by the reducer rows mat->rr and among each
 * other, the new pivots are returned in mat->tr, ordered by decreasing
 * lead column, returns the number of rows reduced to zero */
static len_t exact_sparse_reduced_echelon_form_zz(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i;

    const len_t ncols = mat->nc;
    const len_t nru   = mat->nru;
    const int nthrds  = st->nthrds > 0 ? st->nthrds : 1;
    st->nthrds        = nthrds;

    zz_row_t **pivs = (zz_row_t **)calloc((unsigned long)ncols,
            sizeof(zz_row_t *));
    zz_row_t *kn    = (zz_row_t *)calloc((unsigned long)nru, sizeof(zz_row_t));
    for (i = 0; i < nru; ++i) {
        kn[i].row = mat->rr[i];
        kn[i].cf  = bs->cf_qq[mat->rr[i][COEFFS]];
        zz_set_c32(kn + i, kn[i].cf);
        pivs[mat->rr[i][OFFSET]] = kn + i;
    }
    uint64_t gen    = 0;
    uint64_t *cgen  = (uint64_t *)calloc((unsigned long)ncols, sizeof(uint64_t));
    zz_dense_t *ds  = zz_init_dense_rows(&gen, cgen, ncols, nthrds);

    /* average length of the known pivots */
    uint64_t apl = 1;
    if (nru > 0) {
        for (i = 0; i < nru; ++i) {
            apl += mat->rr[i][LENGTH];
        }
        apl = apl / nru + 1;
    }
    const len_t nzero = zz_reduce_rows(pivs, mat->tr, mat->nrl, bs, ds,
            ncols, apl, nthrds);
    zz_interreduce_new_pivots(ds, pivs, ncols);
    zz_collect_new_pivots(mat, pivs, ncols, 0);

    /* the reducer rows, replaced by a gcd step or not, are not needed
     * anymore, their coefficients belong to the basis */
    for (i = 0; i < nru; ++i) {
        free(kn[i].c32);
        free(mat->rr[i]);
        mat->rr[i] = NULL;
    }
    free(kn);
    zz_free_dense_rows(ds, st, ncols);
    free(cgen);
    free(pivs);

    return nzero;
}

static void exact_sparse_linear_algebra_zz(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t nzero = exact_sparse_reduced_echelon_form_zz(mat, bs, st);
    st->np = mat->np;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += nzero;
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%7d new %7d zero", mat->np, nzero);
        fflush(VERBSTREAM);
    }
}

/* computes the reduced echelon form of the rows mat->rr, which may have
 * equal lead columns, the result is returned in mat->tr, ordered by
 * increasing lead column */
static void interreduce_matrix_rows_zz(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        int free_basis
        )
{
    const len_t ncols = mat->nc;
    const int nthrds  = st->nthrds > 0 ? st->nthrds : 1;
    st->nthrds        = nthrds;

    /* adjust displaying timings for statistic printout */
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "                        ");
    }

    zz_row_t **pivs = (zz_row_t **)calloc((unsigned long)ncols,
            sizeof(zz_row_t *));
    uint64_t gen    = 0;
    uint64_t *cgen  = (uint64_t *)calloc((unsigned long)ncols, sizeof(uint64_t));
    zz_dense_t *ds  = zz_init_dense_rows(&gen, cgen, ncols, nthrds);
    /* average length of the rows */
    uint64_t apl = 1;
    if (mat->nr > 0) {
        for (len_t i = 0; i < mat->nr; ++i) {
            apl += mat->rr[i][LENGTH];
        }
        apl = apl / mat->nr + 1;
    }
    zz_reduce_rows(pivs, mat->rr, mat->nr, bs, ds, ncols, apl, nthrds);
    zz_interreduce_new_pivots(ds, pivs, ncols);

    if (free_basis != 0) {
        /* free now all polynomials in the basis and reset bs->ld to 0. */
        free_basis_elements(bs);
    }
    free(mat->rr);
    mat->rr = NULL;
    zz_collect_new_pivots(mat, pivs, ncols, 1);
    zz_free_dense_rows(ds, st, ncols);
    free(cgen);
    free(pivs);
}
