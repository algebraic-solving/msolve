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


/* Reduction kernels for small prime fields.
 *
 * As long as the field characteristic fc is < 2^18 we can accumulate
 * products of multipliers and coefficients in int64_t without any
 * intermediate reduction. This file is a template: it is included by
 * la_ff_8.c, la_ff_16.c and la_ff_32.c (for primes < 2^18), see
 * la_ff_common.c for the macros LA_CF_T, LA_CF and LA_FF. Additionally,
 * the macro LA_SMALL(name) gives the name of the generated kernels,
 * e.g. foo_ff_8 or foo_17_bit. LA_SMALL is undefined at the end of this
 * file.
 *
 * If LA_DENSE_BITS is defined,
 * LA_SMALL(reduce_dense_row_by_dense_new_pivots) works on entries of
 * LA_DENSE_BITS bits instead of int64_t: 16 bits in la_ff_8.c and 32 bits
 * in la_ff_16.c, i.e. p^2 < 2^LA_DENSE_BITS. The entries are kept in
 * [0,p^2] as in LA_FF(reduce_dense_rows_by_known_pivots), see
 * la_ff_common.c, and the updates are vectorized. LA_DENSE_BITS is
 * undefined at the end of this file. */

#ifdef LA_DENSE_BITS
#if defined HAVE_AVX2 || defined HAVE_AVX512_F
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif
#if LA_DENSE_BITS == 16
#define LA_DENSE_T uint16_t
#else
#define LA_DENSE_T uint32_t
#endif
#endif

static hm_t *LA_SMALL(trace_reduce_dense_row_by_known_pivots_sparse)(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        md_t *st
        )
{
    hi_t i, j, k;
    hm_t *dts;
    LA_CF_T *cfs;
    int64_t np = -1;
    /* statistics are accumulated locally, updating st for each reducer
     * makes all threads compete for the same cache line */
    uint64_t nops = 0, nred = 0;
    const int64_t mod           = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    LA_CF_T * const * const mcf  = mat->LA_CF;

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->LA_CF[dts[COEFFS]];
            /* set corresponding bit of reducer in reducer bit array */
            rba[i/32] |= 1U << (i % 32);
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds  = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] +=  mul * cfs[j];
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   +=  mul * cfs[j];
            dr[ds[j+1]] +=  mul * cfs[j+1];
            dr[ds[j+2]] +=  mul * cfs[j+2];
            dr[ds[j+3]] +=  mul * cfs[j+3];
        }
        dr[i] = 0;
        nops  +=  len;
        nred++;
    }
    st->trace_nr_mult +=  nops / 1000.0;
    st->trace_nr_add  +=  nops / 1000.0;
    st->trace_nr_red  +=  nred;

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    LA_CF_T *cf  = (LA_CF_T *)malloc((uint64_t)(k) * sizeof(LA_CF_T));
    j = 0;
    hm_t *rs = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (LA_CF_T)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->LA_CF[tmp_pos]  = cf;

    return row;
}

static LA_CF_T *LA_SMALL(reduce_dense_row_by_all_pivots)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const * const pivs,
        LA_CF_T *const * const dpivs,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    const int64_t mod = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);
    len_t np  = -1;
    LA_CF_T *red;

    const len_t ncl   = mat->ncl;
    const len_t ncols = mat->nc;

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        const LA_CF_T *cfs = bs->LA_CF[pivs[i][COEFFS]];
        const len_t os    = pivs[i][PRELOOP];
        const len_t len   = pivs[i][LENGTH];
        const hm_t * const ds  = pivs[i] + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] +=  mul * cfs[j];
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   +=  mul * cfs[j];
            dr[ds[j+1]] +=  mul * cfs[j+1];
            dr[ds[j+2]] +=  mul * cfs[j+2];
            dr[ds[j+3]] +=  mul * cfs[j+3];
        }
        dr[i] = 0;
    }
    k = 0;
    /* step 2: reduce by new dense pivots */
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (dpivs[i-ncl] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        red = dpivs[i-ncl];
        const int64_t mul = mod - dr[i];
        const len_t os    = (ncols - i) % UNROLL;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] +=  mul * red[l];
        }
        for (; j < ncols; l += UNROLL, j += UNROLL) {
            dr[j]   +=  mul * red[l];
            dr[j+1] +=  mul * red[l+1];
            dr[j+2] +=  mul * red[l+2];
            dr[j+3] +=  mul * red[l+3];
        }
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    LA_CF_T *row = (LA_CF_T *)calloc((uint64_t)(ncols-np), sizeof(LA_CF_T));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        row[i-np]  = (LA_CF_T)dr[i];
    }
    if (row[0] != 1) {
        row = LA_FF(normalize_dense_matrix_row)(row, ncols-np, fc);
    }

    *pc = np - ncl;
    return row;
}

#ifdef LA_DENSE_BITS
/* t[j] = t[j] - m*c[j] for 0 <= j < len with a borrow corrected by adding
 * p^2, i.e. the entries stay in [0,p^2], see
 * LA_FF(reduce_dense_rows_by_known_pivots) */
static inline void LA_SMALL(dense_submul)(
        LA_DENSE_T *t,
        const LA_CF_T * const c,
        const len_t len,
        const LA_DENSE_T m,
        const LA_DENSE_T mod2
        )
{
    len_t j = 0;
#if LA_DENSE_BITS == 16 && defined HAVE_AVX512_BW
    const __m512i mv    = _mm512_set1_epi16((int16_t)m);
    const __m512i mod2v = _mm512_set1_epi16((int16_t)mod2);
    for (; j + 32 <= len; j += 32) {
        const __m512i prodv = _mm512_mullo_epi16(mv,
                _mm512_cvtepu8_epi16(_mm256_loadu_si256((__m256i *)(c+j))));
        const __m512i tv    = _mm512_loadu_si512((void *)(t+j));
        const __mmask32 bv  = _mm512_cmple_epu16_mask(tv, prodv);
        const __m512i resv  = _mm512_sub_epi16(tv, prodv);
        _mm512_storeu_si512((void *)(t+j),
                _mm512_mask_add_epi16(resv, bv, resv, mod2v));
    }
#elif LA_DENSE_BITS == 16 && (defined HAVE_AVX2 || defined HAVE_AVX512_F)
    const __m256i mv    = _mm256_set1_epi16((int16_t)m);
    const __m256i mod2v = _mm256_set1_epi16((int16_t)mod2);
    for (; j + 16 <= len; j += 16) {
        const __m256i prodv = _mm256_mullo_epi16(mv,
                _mm256_cvtepu8_epi16(_mm_loadu_si128((__m128i *)(c+j))));
        const __m256i tv    = _mm256_loadu_si256((__m256i *)(t+j));
        const __m256i bv    = _mm256_cmpeq_epi16(
                _mm256_min_epu16(tv, prodv), tv);
        _mm256_storeu_si256((__m256i *)(t+j), _mm256_add_epi16(
                    _mm256_sub_epi16(tv, prodv),
                    _mm256_and_si256(bv, mod2v)));
    }
#elif LA_DENSE_BITS == 16 && defined __aarch64__
    const uint16x8_t mv     = vdupq_n_u16(m);
    const uint16x8_t mod2v  = vdupq_n_u16(mod2);
    for (; j + 8 <= len; j += 8) {
        const uint16x8_t prodv  = vmulq_u16(mv, vmovl_u8(vld1_u8(c+j)));
        const uint16x8_t tv     = vld1q_u16(t+j);
        const uint16x8_t bv     = vcleq_u16(tv, prodv);
        vst1q_u16(t+j, vaddq_u16(vsubq_u16(tv, prodv), vandq_u16(bv, mod2v)));
    }
#elif LA_DENSE_BITS == 32 && defined HAVE_AVX512_F
    const __m512i mv    = _mm512_set1_epi32((int32_t)m);
    const __m512i mod2v = _mm512_set1_epi32((int32_t)mod2);
    for (; j + 16 <= len; j += 16) {
        const __m512i prodv = _mm512_mullo_epi32(mv,
                _mm512_cvtepu16_epi32(_mm256_loadu_si256((__m256i *)(c+j))));
        const __m512i tv    = _mm512_loadu_si512((void *)(t+j));
        const __mmask16 bv  = _mm512_cmple_epu32_mask(tv, prodv);
        const __m512i resv  = _mm512_sub_epi32(tv, prodv);
        _mm512_storeu_si512((void *)(t+j),
                _mm512_mask_add_epi32(resv, bv, resv, mod2v));
    }
#elif LA_DENSE_BITS == 32 && defined HAVE_AVX2
    const __m256i mv    = _mm256_set1_epi32((int32_t)m);
    const __m256i mod2v = _mm256_set1_epi32((int32_t)mod2);
    for (; j + 8 <= len; j += 8) {
        const __m256i prodv = _mm256_mullo_epi32(mv,
                _mm256_cvtepu16_epi32(_mm_loadu_si128((__m128i *)(c+j))));
        const __m256i tv    = _mm256_loadu_si256((__m256i *)(t+j));
        const __m256i bv    = _mm256_cmpeq_epi32(
                _mm256_min_epu32(tv, prodv), tv);
        _mm256_storeu_si256((__m256i *)(t+j), _mm256_add_epi32(
                    _mm256_sub_epi32(tv, prodv),
                    _mm256_and_si256(bv, mod2v)));
    }
#elif LA_DENSE_BITS == 32 && defined __aarch64__
    const uint32x4_t mv     = vdupq_n_u32(m);
    const uint32x4_t mod2v  = vdupq_n_u32(mod2);
    for (; j + 4 <= len; j += 4) {
        const uint32x4_t prodv  = vmulq_u32(mv, vmovl_u16(vld1_u16(c+j)));
        const uint32x4_t tv     = vld1q_u32(t+j);
        const uint32x4_t bv     = vcleq_u32(tv, prodv);
        vst1q_u32(t+j, vaddq_u32(vsubq_u32(tv, prodv), vandq_u32(bv, mod2v)));
    }
#endif
    for (; j < len; ++j) {
        const LA_DENSE_T prod = (LA_DENSE_T)(m * c[j]);
        t[j] = (LA_DENSE_T)(t[j] - prod + (t[j] <= prod ? mod2 : 0));
    }
}

/* Reduces the dense row dr from column *pc on by the dense new pivots
 * pivs, see the int64_t version below. The row is copied to entries of
 * LA_DENSE_BITS bits, reduced there and copied back to dr, reduced
 * modulo fc, so dr can be reduced again from the new *pc on if another
 * thread added a pivot in this column in the meantime. */
static LA_CF_T *LA_SMALL(reduce_dense_row_by_dense_new_pivots)(
        int64_t *dr,
        len_t *pc,
        LA_CF_T * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        )
{
    len_t i;
    len_t np  = -1;
    hi_t k    = 0;
    const int64_t mod       = (int64_t)fc;
    const uint64_t bm       = barrett_multiplier(mod);
    const LA_DENSE_T mod2   = (LA_DENSE_T)(fc * fc);
    const len_t sc          = *pc;

    /* t[i] is the entry of dr in column sc+i */
    LA_DENSE_T *t = (LA_DENSE_T *)malloc(
            (uint64_t)(ncr-sc) * sizeof(LA_DENSE_T));
    for (i = sc; i < ncr; ++i) {
        t[i-sc] = dr[i] != 0 ?
            (LA_DENSE_T)barrett_reduce((uint64_t)dr[i], mod, bm) : 0;
    }
    for (i = sc; i < ncr; ++i) {
        LA_DENSE_T * const ti = t + (i-sc);
        if (*ti == 0) {
            continue;
        }
        *ti = (LA_DENSE_T)barrett_reduce((uint64_t)*ti, mod, bm);
        if (*ti == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }
        /* the pivot starts with 1 in column i */
        LA_SMALL(dense_submul)(ti+1, pivs[i]+1, ncr-i-1, *ti, mod2);
        *ti = 0;
    }
    for (i = sc; i < ncr; ++i) {
        dr[i] = t[i-sc] != 0 ?
            (int64_t)barrett_reduce((uint64_t)t[i-sc], mod, bm) : 0;
    }
    free(t);
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    LA_CF_T *row = (LA_CF_T *)malloc((uint64_t)(ncr-np) * sizeof(LA_CF_T));
    for (i = np; i < ncr; ++i) {
        row[i-np]  = (LA_CF_T)dr[i];
    }
    if (row[0] != 1) {
        row = LA_FF(normalize_dense_matrix_row)(row, ncr-np, fc);
    }
    *pc = np;
    return row;
}
#else
static LA_CF_T *LA_SMALL(reduce_dense_row_by_dense_new_pivots)(
        int64_t *dr,
        len_t *pc,
        LA_CF_T * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    len_t np  = -1;
    const int64_t mod = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);

    for (k = 0, i = *pc; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        const int64_t mul = mod - dr[i];
        const len_t os    = (ncr - i) % UNROLL;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] +=  mul * pivs[i][l];
        }
        for (; j < ncr; l += UNROLL, j += UNROLL) {
            dr[j]   +=  mul * pivs[i][l];
            dr[j+1] +=  mul * pivs[i][l+1];
            dr[j+2] +=  mul * pivs[i][l+2];
            dr[j+3] +=  mul * pivs[i][l+3];
        }
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    LA_CF_T *row = (LA_CF_T *)calloc((uint64_t)(ncr-np), sizeof(LA_CF_T));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        row[i-np]  = (LA_CF_T)dr[i];
    }
    if (row[0] != 1) {
        row = LA_FF(normalize_dense_matrix_row)(row, ncr-np, fc);
    }
    *pc = np;
    return row;
}
#endif

#undef LA_SMALL
#ifdef LA_DENSE_BITS
#undef LA_DENSE_T
#undef LA_DENSE_BITS
#endif
