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

#include "data.h"
#include "../msolve/streams.h"

#if defined(HAVE_AVX2) || defined(HAVE_AVX512_F)
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

static inline cf32_t *normalize_dense_matrix_row_ff_32(
        cf32_t *row,
        const hm_t len,
        const uint32_t fc
        )
{
    len_t i;

    const hm_t os       = len % UNROLL;
    const uint32_t inv  = mod_p_inverse_32(row[0], fc);

    for (i = 1; i < os; ++i) {
        row[i]  = (cf32_t)(((uint64_t)row[i] * inv) % fc);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    = (cf32_t)(((uint64_t)row[i] * inv) % fc);
        row[i+1]  = (cf32_t)(((uint64_t)row[i+1] * inv) % fc);
        row[i+2]  = (cf32_t)(((uint64_t)row[i+2] * inv) % fc);
        row[i+3]  = (cf32_t)(((uint64_t)row[i+3] * inv) % fc);
    }
    row[0]  = 1;

    return row;
}

static inline void normalize_dense_matrix_row_kernel_ff_32(
        cf32_t *row,
        const len_t start,
        const len_t nc,
        const uint32_t fc
        )
{
    len_t i;

    const len_t len     = nc;
    const len_t os      = len % UNROLL;
    const int64_t inv  =  (int64_t)mod_p_inverse_32(row[start], fc);

    for (i = start; i < os; ++i) {
        row[i]  = (cf32_t)(((uint64_t)row[i] * inv) % fc);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    = (cf32_t)(((uint64_t)row[i] * inv) % fc);
        row[i+1]  = (cf32_t)(((uint64_t)row[i+1] * inv) % fc);
        row[i+2]  = (cf32_t)(((uint64_t)row[i+2] * inv) % fc);
        row[i+3]  = (cf32_t)(((uint64_t)row[i+3] * inv) % fc);
    }
}

/* unused at the moment */
#if 0
static int reduce_dense_row_by_dense_pivots_ff_32(
        int64_t *dr,
        cf32_t *dm,
        cf32_t **pivs,
        const len_t start,
        const len_t nc,
        const len_t fc
        )
{
    len_t i, j;

    const int64_t mod   = (int64_t)fc;

    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2  = (int64_t)fc * fc;

    for (i = start; i < nc; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            /* add new pivot to dm */
            /* for (j = 0; j < i; ++j) {
             *     if (dr[j] != 0) {
             *         printf("dr =/= 0 at pos %u / %u || %u\n", j, i, start);
             *         break;
             *     }
             * } */
            cf32_t *np  = dm + i*nc;
            for (j = i; j < nc; ++j) {
                dr[j] = (int64_t)barrett_reduce((uint64_t)dr[j], mod, bm);
                np[j] = (cf32_t)dr[j];
            }
            /* normalize new pivot row */
            normalize_dense_matrix_row_kernel_ff_32(
                    np, i, nc, fc);

            pivs[i] = np;

            return i;
        }

        const int64_t mul = (int64_t)dr[i];
        const cf32_t *red = pivs[i];

        const len_t os  = ((nc-i) % 4) + i;

        /* for (j = 0; j < nc; ++j) {
         *     dr[j] -=  mul * red[j];
         *     dr[j] +=  (dr[j] >> 63) & mod2;
         * } */
        for (j = i; j < os; ++j) {
            dr[j] -=  mul * red[j];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        for (; j < nc; j += 4) {
            dr[j]   -=  mul * red[j];
            dr[j]   +=  (dr[j] >> 63) & mod2;
            dr[j+1] -=  mul * red[j+1];
            dr[j+1] +=  (dr[j+1] >> 63) & mod2;
            dr[j+2] -=  mul * red[j+2];
            dr[j+2] +=  (dr[j+2] >> 63) & mod2;
            dr[j+3] -=  mul * red[j+3];
            dr[j+3] +=  (dr[j+3] >> 63) & mod2;
        }
    }
    return -1;
}
#endif

/* is_kernel_trivial() is unused at the moment */
#if 0
/* If the matrix has size m (rows) x n (columns) with m << n then
 * we do the following: Generate a dense m x m matrix where the entries
 * in the mth column is a random linear combination of column entries
 * m up to n from the correspdoning row. Next we check if the kernel
 * for this m x m matrix is trivial or not. */
static int is_kernel_trivial(
        bs_t *sat,
        hm_t **upivs, /* already sorted by increasing pivot */
        const len_t ncl, /* number of columns left, all of them are zero */
        const len_t ncr, /* number of columns right */
        const md_t * const st
        )
{
    len_t i, j;

    len_t max_col = 0;
    for (i = 0; i < sat->ld; ++i) {
        max_col = max_col < sat->hm[i][OFFSET] ? sat->hm[i][OFFSET] : max_col;
    }
    int ret = 0;
    const int64_t fc   = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;
    /* make a dense matrix first */
    const len_t dim = max_col - ncl;
    /* allocate memory */
    cf32_t *dm    = (cf32_t *)calloc((uint64_t)(sat->ld*dim), sizeof(cf32_t));
    int64_t *dr   = (int64_t *)calloc((uint64_t)dim, sizeof(int64_t));
    int64_t *mul  = (int64_t *)calloc((uint64_t)ncr, sizeof(int64_t));
    cf32_t **pivs = (cf32_t **)calloc((uint64_t)dim, sizeof(cf32_t *));

    /* fill random value array */
    for (i = 0; i < ncr; ++i) {
        do {
            mull[j] = (int64_t)rand() % fc;
        } while (mull[j] == 0);
    }
    for (i = 0; i < sat->ld; ++i) {
        memset(dr, 0, (uint64_t)dim * sizeof(int64_t));
        hm_t *npiv  = upivs[i];
        const len_t len = upivs[i][LENGTH];
        cf32_t *cf  = sat->cf_32[npiv[COEFFS]];
        hm_t *cols  = npiv + OFFSET;
        j = 0;
        /* while (j < len) {
         *     dr[cols[j]-ncl] = cf[j];
         *     j++;
         * } */
        while (j < len && cols[j] < dim+ncl) {
            dr[cols[j]-ncl] = cf[j];
            j++;
        }
        while (j < len) {
            dr[dim-1] -=  mul[cols[j]-ncl] * cf[j];
            dr[dim-1] +=  (dr[dim-1] >> 63) & mod2;
            j++;
        }
        ret = reduce_dense_row_by_dense_pivots_ff_32(
                dr, dm, pivs, cols[0]-ncl, dim, fc);
                /* dr, dm, pivs, cols[0]-ncl, ncr, fc); */
        if (ret == -1) {
            free(dm);
            free(dr);
            free(mul);
            free(pivs);
            return 0;
        }
    }
    free(dm);
    free(dr);
    free(mul);
    free(pivs);

    return 1;
}
#endif

static inline cf32_t *normalize_sparse_matrix_row_ff_32(
        cf32_t *row,
        const len_t os,
        const len_t len,
        const uint32_t fc
        )
{
    len_t i;

    const uint32_t inv  = mod_p_inverse_32(row[0], fc);

    for (i = 0; i < os; ++i) {
        row[i]  =  (cf32_t)(((uint64_t)row[i] * inv) % fc);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    = (cf32_t)(((uint64_t)row[i] * inv) % fc);
        row[i+1]  = (cf32_t)(((uint64_t)row[i+1] * inv) % fc);
        row[i+2]  = (cf32_t)(((uint64_t)row[i+2] * inv) % fc);
        row[i+3]  = (cf32_t)(((uint64_t)row[i+3] * inv) % fc);
    }
    row[0]  = 1;

    return row;
}

static inline cf32_t *adjust_multiplier_sparse_matrix_row_ff_32(
        cf32_t *row,
        const cf32_t * const cfs,
        const len_t os,
        const len_t len,
        const uint32_t fc
        )
{
    len_t i;

    const uint32_t inv  = mod_p_inverse_32(cfs[0], fc);

    for (i = 0; i < os; ++i) {
        row[i]  =  (cf32_t)(((uint64_t)row[i] * inv) % fc);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    = (cf32_t)(((uint64_t)row[i] * inv) % fc);
        row[i+1]  = (cf32_t)(((uint64_t)row[i+1] * inv) % fc);
        row[i+2]  = (cf32_t)(((uint64_t)row[i+2] * inv) % fc);
        row[i+3]  = (cf32_t)(((uint64_t)row[i+3] * inv) % fc);
    }

    return row;
}

static inline cf32_t *multiply_sparse_matrix_row_ff_32(
        cf32_t *row,
        const cf32_t mul,
        const len_t os,
        const len_t len,
        const uint32_t fc
        )
{
    len_t i;
    for (i = 0; i < os; ++i) {
        row[i]  =  (cf32_t)(((uint64_t)row[i] * mul) % fc);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    = (cf32_t)(((uint64_t)row[i] * mul) % fc);
        row[i+1]  = (cf32_t)(((uint64_t)row[i+1] * mul) % fc);
        row[i+2]  = (cf32_t)(((uint64_t)row[i+2] * mul) % fc);
        row[i+3]  = (cf32_t)(((uint64_t)row[i+3] * mul) % fc);
    }

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        )
{
    hi_t i, j, k;
    hm_t *dts;
    cf32_t *cfs;
    int64_t np  = -1;
    const int64_t mod           = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;

    rba_t *rba;
    if (tr > 0) {
        rba = mat->rba[tmp_pos];
    } else {
        rba = NULL;
    }
#if defined HAVE_AVX512_F
    int64_t res[8] __attribute__((aligned(64)));
    __m512i redv, mulv, prodv, drv, resv;
#elif defined HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i redv, mulv, prodv, drv, resv;
#elif defined __aarch64__
    uint64_t tmp[2] __attribute__((aligned(32)));
    uint32x4_t redv;
    uint64x2_t drv, resv;
#endif

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
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        }
        cfs   = mcf[dts[COEFFS]];
#if defined HAVE_AVX512_F
        const len_t len = dts[LENGTH];
        const len_t os  = len % 16;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (int32_t)(mod - dr[i]);
        mulv  = _mm512_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 16) {
            redv  = _mm512_loadu_si512((__m512i*)(cfs+j));
            drv   = _mm512_setr_epi64(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]],
                dr[ds[j+9]],
                dr[ds[j+11]],
                dr[ds[j+13]],
                dr[ds[j+15]]);
            /* first four mult-adds -- lower */
            prodv = _mm512_mul_epu32(mulv, _mm512_srli_epi64(redv, 32));
            resv  = _mm512_add_epi64(drv, prodv);
            _mm512_store_si512((__m512*)(res), resv);
            dr[ds[j+1]]  = res[0];
            dr[ds[j+3]]  = res[1];
            dr[ds[j+5]]  = res[2];
            dr[ds[j+7]]  = res[3];
            dr[ds[j+9]]  = res[4];
            dr[ds[j+11]] = res[5];
            dr[ds[j+13]] = res[6];
            dr[ds[j+15]] = res[7];
            /* second four mult-adds -- higher */
            prodv = _mm512_mul_epu32(mulv, redv);
            drv   = _mm512_setr_epi64(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]],
                dr[ds[j+8]],
                dr[ds[j+10]],
                dr[ds[j+12]],
                dr[ds[j+14]]);
            resv  = _mm512_add_epi64(drv, prodv);
            _mm512_store_si512((__m512i*)(res), resv);
            dr[ds[j]]    = res[0];
            dr[ds[j+2]]  = res[1];
            dr[ds[j+4]]  = res[2];
            dr[ds[j+6]]  = res[3];
            dr[ds[j+8]]  = res[4];
            dr[ds[j+10]] = res[5];
            dr[ds[j+12]] = res[6];
            dr[ds[j+14]] = res[7];
        }
#elif defined HAVE_AVX2
        const len_t len = dts[LENGTH];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (int32_t)(mod - dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 8) {
            redv  = _mm256_lddqu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_add_epi64(drv, prodv);
            _mm256_store_si256((__m256i*)(res), resv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_add_epi64(drv, prodv);
            _mm256_store_si256((__m256i*)(res), resv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
#elif defined __aarch64__
        const len_t len       = dts[LENGTH];
        const len_t os        = len % 16;
        const hm_t * const ds = dts + OFFSET;
        const cf32_t mul32   = (cf32_t)(mod - dr[i]);
        const uint32x2_t mulv = vmov_n_u32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 16) {
            tmp[0] = (uint64_t)dr[ds[j]];
            tmp[1] = (uint64_t)dr[ds[j+1]];
            drv  = vld1q_u64(tmp);
            redv = vld1q_u32((cf32_t *)(cfs)+j);
            resv = vmlal_u32(drv, vget_low_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j]]   = (int64_t)tmp[0];
            dr[ds[j+1]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+2]];
            tmp[1] = (uint64_t)dr[ds[j+3]];
            drv  = vld1q_u64(tmp);
            resv = vmlal_u32(drv, vget_high_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+2]] = (int64_t)tmp[0];
            dr[ds[j+3]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+4]];
            tmp[1] = (uint64_t)dr[ds[j+5]];
            drv  = vld1q_u64(tmp);
            redv = vld1q_u32((cf32_t *)(cfs)+j+4);
            resv = vmlal_u32(drv, vget_low_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+4]] = (int64_t)tmp[0];
            dr[ds[j+5]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+6]];
            tmp[1] = (uint64_t)dr[ds[j+7]];
            drv  = vld1q_u64(tmp);
            resv = vmlal_u32(drv, vget_high_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+6]] = (int64_t)tmp[0];
            dr[ds[j+7]] = (int64_t)tmp[1];

            tmp[0] = (uint64_t)dr[ds[j+8]];
            tmp[1] = (uint64_t)dr[ds[j+9]];
            drv  = vld1q_u64(tmp);
            redv = vld1q_u32((cf32_t *)(cfs)+j+8);
            resv = vmlal_u32(drv, vget_low_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+8]] = (int64_t)tmp[0];
            dr[ds[j+9]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+10]];
            tmp[1] = (uint64_t)dr[ds[j+11]];
            drv  = vld1q_u64(tmp);
            resv = vmlal_u32(drv, vget_high_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+10]] = (int64_t)tmp[0];
            dr[ds[j+11]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+12]];
            tmp[1] = (uint64_t)dr[ds[j+13]];
            drv  = vld1q_u64(tmp);
            redv = vld1q_u32((cf32_t *)(cfs)+j+12);
            resv = vmlal_u32(drv, vget_low_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+12]] = (int64_t)tmp[0];
            dr[ds[j+13]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+14]];
            tmp[1] = (uint64_t)dr[ds[j+15]];
            drv  = vld1q_u64(tmp);
            resv = vmlal_u32(drv, vget_high_u32(redv), mulv);
            vst1q_u64(tmp, resv);
            dr[ds[j+14]] = (int64_t)tmp[0];
            dr[ds[j+15]] = (int64_t)tmp[1];
        }
#else
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
#endif
        dr[i] = 0;
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_32[tmp_pos]  = cf;

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_up_to_ff_31_bit(
        int64_t *dr,
        bs_t *sat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const len_t cf_idx,
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t end,   /* column index up to which we reduce */
        const len_t ncols,
        md_t *st
        )
{
    hi_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int64_t np = -1;
    const int64_t mod   = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2  = (int64_t)st->fc * st->fc;
#ifdef HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov= _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#endif

    for (i = dpiv; i < end; ++i) {
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
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        dts = pivs[i];
        cfs = bs->cf_32[dts[COEFFS]];
#ifdef HAVE_AVX2
        const len_t len = dts[LENGTH];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (uint32_t)(dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
#else
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
#endif
        dr[i] = 0;
    }

    k = ncols-end;
    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = end; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    if (j == 0) {
        free(row);
        free(cf);

        return NULL;
    }
    row[COEFFS]   = cf_idx;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    row = realloc(row, (uint64_t)(j+OFFSET) * sizeof(hm_t));
    cf  = realloc(cf, (uint64_t)j * sizeof(cf32_t));
    sat->cf_32[cf_idx]  = cf;

    return row;
}

static hm_t *sba_reduce_dense_row_by_known_pivots_sparse_31_bit(
        int64_t *dr,
        smat_t *smat,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t sm,      /* signature monomial of row reduced */
        const len_t si,     /* signature index of row reduced */
        const len_t ri,     /* index of row in matrix */
        md_t *st
        )
{
    len_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int32_t fnzc       = -1; /* first non zero column */
    const int64_t mod  = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2 = (int64_t)st->fc * st->fc;
    const len_t nc     = smat->nc;
#ifdef HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov= _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#endif

    k = 0;
    for (i = dpiv; i < nc; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (fnzc == -1) {
                if (i == dpiv) {
                    for (j = dpiv; j < nc; ++j) {
                        if (dr[j] != 0) {
                            k++;
                        }
                    }
                    fnzc = i;
                    break;
                }
                fnzc = i;
            }
            k++;
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        dts = pivs[i];
        cfs = smat->cc32[dts[SM_CFS]];

#ifdef HAVE_AVX2
        const len_t len = dts[SM_LEN];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + SM_OFFSET;
        const uint32_t mul32 = (uint32_t)(dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
#else
        const len_t os  = dts[SM_PRE];
        const len_t len = dts[SM_LEN];
        const hm_t * const ds = dts + SM_OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
#endif
        dr[i] = 0;
    }

    if (k == 0) {
        free(smat->cr[ri]);
        smat->cr[ri] = NULL;

        return smat->cr[ri];
    }

    /* printf("k %d | ri %u | cr[%u] = %p\n", k, ri, ri, smat->cr[ri]);
     * printf("%lu\n", (uint64_t)(k+SM_OFFSET)); */
    smat->cr[ri] = realloc(smat->cr[ri],
            (uint64_t)(k+SM_OFFSET) * sizeof(hm_t));
    cf32_t *cf   = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));

    j = 0;
    hm_t *rs  = smat->cr[ri] + SM_OFFSET;
    for (i = fnzc; i < nc; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    smat->cr[ri][SM_SMON] = sm;
    smat->cr[ri][SM_SIDX] = si;
    smat->cr[ri][SM_CFS]  = ri;
    smat->cr[ri][SM_PRE]  = j % UNROLL;
    smat->cr[ri][SM_LEN]  = j;

    smat->cc32[ri] = cf;

    return smat->cr[ri];
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        )
{
    hi_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int64_t np = -1;
    const int64_t mod           = (int64_t)st->fc;
    const int64_t mod2          = (int64_t)st->fc * st->fc;
    const uint64_t bm           = barrett_multiplier(st->fc);
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;

    rba_t *rba;
    if (tr > 0) {
        rba = mat->rba[tmp_pos];
    } else {
        rba = NULL;
    }
#if defined HAVE_AVX512_F
    int64_t res[8] __attribute__((aligned(64)));
    __m512i redv, drv, mulv, prodv, resv, rresv;
    __mmask8 cmpv;
    __m512i zerov = _mm512_set1_epi64(0);
    __m512i mod2v = _mm512_set1_epi64(mod2);
#elif defined HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov= _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#elif defined __aarch64__
    const int64x2_t mod2v = vmovq_n_s64(mod2);
    int64_t tmp[2] __attribute__((aligned(32)));
    int32x4_t redv;
    int64x2_t drv, mask, resv;
#endif

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
        const int64_t mul = (int64_t)dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            /* cfs   = bs->cf_32[dts[COEFFS]]; */
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        }
        cfs   = mcf[dts[COEFFS]];
#if defined HAVE_AVX512_F
        const len_t len = dts[LENGTH];
        const len_t os  = len % 16;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (int32_t)(dr[i]);
        mulv  = _mm512_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -= mul * cfs[j];
            dr[ds[j]] += (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 16) {
            redv  = _mm512_loadu_si512((__m512i*)(cfs+j));
            drv   = _mm512_setr_epi64(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]],
                dr[ds[j+9]],
                dr[ds[j+11]],
                dr[ds[j+13]],
                dr[ds[j+15]]);
            /* first four mult-adds -- lower */
            prodv = _mm512_mul_epu32(mulv, _mm512_srli_epi64(redv, 32));
            resv  = _mm512_sub_epi64(drv, prodv);
            cmpv  = _mm512_cmpgt_epi64_mask(zerov, resv);
            rresv = _mm512_mask_add_epi64(resv, cmpv, resv, mod2v);
            _mm512_store_si512((__m512*)(res), rresv);
            dr[ds[j+1]]  = res[0];
            dr[ds[j+3]]  = res[1];
            dr[ds[j+5]]  = res[2];
            dr[ds[j+7]]  = res[3];
            dr[ds[j+9]]  = res[4];
            dr[ds[j+11]] = res[5];
            dr[ds[j+13]] = res[6];
            dr[ds[j+15]] = res[7];
            /* second four mult-adds -- higher */
            prodv = _mm512_mul_epu32(mulv, redv);
            drv   = _mm512_setr_epi64(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]],
                dr[ds[j+8]],
                dr[ds[j+10]],
                dr[ds[j+12]],
                dr[ds[j+14]]);
            resv  = _mm512_sub_epi64(drv, prodv);
            cmpv  = _mm512_cmpgt_epi64_mask(zerov, resv);
            rresv = _mm512_mask_add_epi64(resv, cmpv, resv, mod2v);
            _mm512_store_si512((__m512i*)(res), rresv);
            dr[ds[j]]    = res[0];
            dr[ds[j+2]]  = res[1];
            dr[ds[j+4]]  = res[2];
            dr[ds[j+6]]  = res[3];
            dr[ds[j+8]]  = res[4];
            dr[ds[j+10]] = res[5];
            dr[ds[j+12]] = res[6];
            dr[ds[j+14]] = res[7];
        }
#elif defined HAVE_AVX2
        const len_t len = dts[LENGTH];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (uint32_t)(dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
#elif defined __aarch64__
        const len_t len       = dts[LENGTH];
        const len_t os        = len % 8;
        const hm_t * const ds = dts + OFFSET;
        const int32_t mul32   = (int32_t)(dr[i]);
        const int32x2_t mulv  = vmov_n_s32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            tmp[0] = dr[ds[j]];
            tmp[1] = dr[ds[j+1]];
            drv  = vld1q_s64(tmp);
            redv = vld1q_s32((int32_t *)(cfs)+j);
            /* multiply and subtract */
            resv = vmlsl_s32(drv, vget_low_s32(redv), mulv);
            mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(tmp, resv);
            dr[ds[j]]   = tmp[0];
            dr[ds[j+1]] = tmp[1];
            tmp[0] = dr[ds[j+2]];
            tmp[1] = dr[ds[j+3]];
            drv  = vld1q_s64(tmp);
            resv = vmlsl_s32(drv, vget_high_s32(redv), mulv);
            mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(tmp, resv);
            dr[ds[j+2]] = tmp[0];
            dr[ds[j+3]] = tmp[1];
            tmp[0] = dr[ds[j+4]];
            tmp[1] = dr[ds[j+5]];
            drv  = vld1q_s64(tmp);
            redv = vld1q_s32((int32_t *)(cfs)+j+4);
            /* multiply and subtract */
            resv = vmlsl_s32(drv, vget_low_s32(redv), mulv);
            mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(tmp, resv);
            dr[ds[j+4]] = tmp[0];
            dr[ds[j+5]] = tmp[1];
            tmp[0] = dr[ds[j+6]];
            tmp[1] = dr[ds[j+7]];
            drv  = vld1q_s64(tmp);
            resv = vmlsl_s32(drv, vget_high_s32(redv), mulv);
            mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(tmp, resv);
            dr[ds[j+6]] = tmp[0];
            dr[ds[j+7]] = tmp[1];
        }

#else
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
#endif
        dr[i] = 0;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_32[tmp_pos]  = cf;

    return row;
}


static hm_t *reduce_dense_row_by_known_pivots_sparse_sat_ff_31_bit(
        int64_t *dr,
        int64_t *drm,
        cf32_t **pivcf,
        hm_t **mulh,
        cf32_t **mulcf,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs arrays */
        const len_t sat_ld, /* number of current multipliers */
        const len_t ncols,  /* number of columns */
        md_t *st,
        const len_t ncl
        )
{
    hi_t i, j, k;
    cf32_t *cfs, *cfsm;
    hm_t *dts, *dtsm;
    int64_t np = -1;
    const int64_t mod           = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2          = (int64_t)st->fc * st->fc;
#ifdef HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov= _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#endif

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
            break;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        dts   = pivs[i];
        dtsm  = mulh[dts[MULT]];
        cfsm  = mulcf[dtsm[COEFFS]];
        cfs   = pivcf[dts[COEFFS]];
#ifdef HAVE_AVX2
        const len_t len = dts[LENGTH];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (uint32_t)(dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
        const len_t lenm = dtsm[LENGTH];
        const len_t osm  = lenm % 8;
        const hm_t * const dsm  = dtsm + OFFSET;
        /* const uint32_t mulm32 = (uint32_t)(drm[i]); */
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < osm; ++j) {
            drm[dsm[j]] -=  mul * cfsm[j];
            drm[dsm[j]] +=  (drm[dsm[j]] >> 63) & mod2;
        }
        for (; j < lenm; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfsm+j));
            drv   = _mm256_setr_epi64x(
                drm[dsm[j+1]],
                drm[dsm[j+3]],
                drm[dsm[j+5]],
                drm[dsm[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            drm[dsm[j+1]] = res[0];
            drm[dsm[j+3]] = res[1];
            drm[dsm[j+5]] = res[2];
            drm[dsm[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                drm[dsm[j]],
                drm[dsm[j+2]],
                drm[dsm[j+4]],
                drm[dsm[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            drm[dsm[j]]   = res[0];
            drm[dsm[j+2]] = res[1];
            drm[dsm[j+4]] = res[2];
            drm[dsm[j+6]] = res[3];
        }
#else
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
        const len_t osm   = dtsm[PRELOOP];
        const len_t lenm  = dtsm[LENGTH];
        const hm_t * const dsm = dtsm + OFFSET;
        for (j = 0; j < osm; ++j) {
            drm[dsm[j]] -=  mul * cfsm[j];
            drm[dsm[j]] +=  (drm[dsm[j]] >> 63) & mod2;
        }
        for (; j < lenm; j += UNROLL) {
            drm[dsm[j]]   -=  mul * cfsm[j];
            drm[dsm[j+1]] -=  mul * cfsm[j+1];
            drm[dsm[j+2]] -=  mul * cfsm[j+2];
            drm[dsm[j+3]] -=  mul * cfsm[j+3];
            drm[dsm[j]]   +=  (drm[dsm[j]] >> 63) & mod2;
            drm[dsm[j+1]] +=  (drm[dsm[j+1]] >> 63) & mod2;
            drm[dsm[j+2]] +=  (drm[dsm[j+2]] >> 63) & mod2;
            drm[dsm[j+3]] +=  (drm[dsm[j+3]] >> 63) & mod2;
        }
#endif
        dr[i] = 0;
    }

    /* we always have to write the multiplier entry for the kernel
     * construction, especially if we computed a zero row. */
    mulh[tmp_pos]   =
        (hm_t *)malloc((uint64_t)(sat_ld+OFFSET) * sizeof(hm_t));
    mulcf[tmp_pos]  =
        (cf32_t *)malloc((uint64_t)(sat_ld) * sizeof(cf32_t));
    j = 0;
    hm_t *rsm   = mulh[tmp_pos] + OFFSET;
    cf32_t*cfm  = mulcf[tmp_pos];
    for (i = 0; i < sat_ld; ++i) {
        if (drm[i] != 0) {
            drm[i]  = drm[i] % mod;
            if (drm[i] != 0) {
                rsm[j]  = (hm_t)i;
                cfm[j] = (cf32_t)drm[i];
                j++;
            }
        }
    }
    mulh[tmp_pos][PRELOOP]  = j % UNROLL;
    mulh[tmp_pos][LENGTH]   = j;
    mulh[tmp_pos][COEFFS]   = tmp_pos;

    mulh[tmp_pos] = realloc(mulh[tmp_pos],
         (uint64_t)(j+OFFSET) * sizeof(hm_t));
    mulcf[tmp_pos]  = realloc(mulcf[tmp_pos],
            (uint64_t)j * sizeof(cf32_t));


    if (k == 0) {
        return NULL;
    }

    k = ncols - np;
    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    row[COEFFS]     = tmp_pos;
    row[PRELOOP]    = j % UNROLL;
    row[LENGTH]     = j;
    row[MULT]       = tmp_pos;
    pivcf[tmp_pos]  = cf;

    return row;
}

static hm_t *trace_reduce_dense_row_by_known_pivots_sparse_31_bit(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        md_t *st
        )
{
    hi_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int64_t np = -1;
    /* statistics are accumulated locally, updating st for each reducer
     * makes all threads compete for the same cache line */
    uint64_t nops = 0, nred = 0;
    const int64_t mod           = (int64_t)st->fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2          = (int64_t)st->fc * st->fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;
#ifdef HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov = _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#endif

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
        const int64_t mul = (int64_t)dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_32[dts[COEFFS]];
            /* set corresponding bit of reducer in reducer bit array */
            rba[i/32] |= 1U << (i % 32);
        } else {
            cfs   = mcf[dts[COEFFS]];
        }

#ifdef HAVE_AVX2
        const len_t len = dts[LENGTH];
        const len_t os  = len % 8;
        const hm_t * const ds  = dts + OFFSET;
        const uint32_t mul32 = (uint32_t)(dr[i]);
        mulv  = _mm256_set1_epi32(mul32);
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 8) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+3]],
                dr[ds[j+5]],
                dr[ds[j+7]]);
            /* first four mult-adds -- lower */
            prodv = _mm256_mul_epu32(mulv, _mm256_srli_epi64(redv, 32));
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j+1]] = res[0];
            dr[ds[j+3]] = res[1];
            dr[ds[j+5]] = res[2];
            dr[ds[j+7]] = res[3];
            /* second four mult-adds -- higher */
            prodv = _mm256_mul_epu32(mulv, redv);
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+2]],
                dr[ds[j+4]],
                dr[ds[j+6]]);
            resv  = _mm256_sub_epi64(drv, prodv);
            cmpv  = _mm256_cmpgt_epi64(zerov, resv);
            rresv = _mm256_add_epi64(resv, _mm256_and_si256(cmpv, mod2v));
            _mm256_store_si256((__m256i*)(res), rresv);
            dr[ds[j]]   = res[0];
            dr[ds[j+2]] = res[1];
            dr[ds[j+4]] = res[2];
            dr[ds[j+6]] = res[3];
        }
#else
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
#endif
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
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_32[tmp_pos]  = cf;

    return row;
}

/* unused at the moment */
#if 0
static hm_t *reduce_dense_row_by_known_pivots_sparse_32_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        )
{
    hi_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int64_t np = -1;
    const uint64_t mod          = (uint64_t)st->fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;
    rba_t *rba;
    if (tr > 0) {
        rba = mat->rba[tmp_pos];
    } else {
        rba = NULL;
    }
    const uint64_t mask = (uint32_t)0xFFFFFFFF;
    uint64_t RED_32     = ((uint64_t)2<<31) % st->fc;
    uint64_t RED_64     = ((uint64_t)1<<63) % st->fc;
    RED_64              = (RED_64*2) % st->fc;
    uint64_t drlow[ncols];
    uint64_t drhigh[ncols];
    uint64_t udr[ncols];
    uint64_t prod;
    for (i = 0; i < ncols; ++i) {
      drlow[i]  = dr[i] & mask;
      drhigh[i] = dr[i] >> 32;
    }

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        udr[i] =   ((drhigh[i] >> 32) * RED_64) % mod;
        udr[i] +=  ((drhigh[i] & (uint64_t)0xFFFFFFFF) * RED_32) % mod;
        udr[i] += drlow[i];
        udr[i] %=  mod;
        if (udr[i] == 0) {
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
        const uint64_t mul = mod-udr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_32[dts[COEFFS]];
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        /* loop unrolling is slower */
        for (j = 0; j < len; ++j) {
            prod  =  mul * cfs[j];
            drhigh[ds[j]] +=  prod >> 32;
            drlow[ds[j]]  +=  (prod & mask);
        }
        udr[i]  = 0;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = np; i < ncols; ++i) {
        if (udr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)udr[i];
            j++;
        }
    }
    row[BINDEX]         = bi;
    row[MULT]           = mh;
    row[COEFFS]         = tmp_pos;
    row[PRELOOP]        = j % UNROLL;
    row[LENGTH]         = j;
    mat->cf_32[tmp_pos] = cf;

    return row;
}
#endif

/* unused at the moment */
#if 0
static hm_t *trace_reduce_dense_row_by_known_pivots_sparse_32_bit(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        md_t *st
        )
{
    hi_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    int64_t np = -1;
    const uint64_t mod          = (uint64_t)st->fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;
    const uint64_t mask = (uint32_t)0xFFFFFFFF;
    uint64_t RED_32     = ((uint64_t)2<<31) % st->fc;
    uint64_t RED_64     = ((uint64_t)1<<63) % st->fc;
    RED_64              = (RED_64*2) % st->fc;
    uint64_t drlow[ncols];
    uint64_t drhigh[ncols];
    uint64_t udr[ncols];
    uint64_t prod;
    for (i = 0; i < ncols; ++i) {
      drlow[i]  = dr[i] & mask;
      drhigh[i] = dr[i] >> 32;
    }

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        udr[i] =   ((drhigh[i] >> 32) * RED_64) % mod;
        udr[i] +=  ((drhigh[i] & (uint64_t)0xFFFFFFFF) * RED_32) % mod;
        udr[i] += drlow[i];
        udr[i] %=  mod;
        if (udr[i] == 0) {
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
        const uint64_t mul = mod-udr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_32[dts[COEFFS]];
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
        const len_t len = dts[LENGTH];
        const hm_t * const ds = dts + OFFSET;
        /* loop unrolling is slower */
        for (j = 0; j < len; ++j) {
            prod  =  mul * cfs[j];
            drhigh[ds[j]] +=  prod >> 32;
            drlow[ds[j]]  +=  (prod & mask);
        }
        udr[i]  = 0;
        st->trace_nr_mult +=  len / 1000.0;
        st->trace_nr_add  +=  len / 1000.0;
        st->trace_nr_red++;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((uint64_t)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = np; i < ncols; ++i) {
        if (udr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf32_t)udr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_32[tmp_pos]  = cf;

    return row;
}
#endif

static cf32_t *reduce_dense_row_by_all_pivots_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t * const * const pivs,
        cf32_t * const * const dpivs,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    const int64_t mod   = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2  = (int64_t)fc * fc;
    len_t np  = -1;
    cf32_t *red;

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
        const int64_t mul = (int64_t)dr[i];
        const cf32_t *cfs = bs->cf_32[pivs[i][COEFFS]];
        const len_t os    = pivs[i][PRELOOP];
        const len_t len   = pivs[i][LENGTH];
        const hm_t * const ds = pivs[i] + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
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
        const int64_t mul = (int64_t) dr[i];
        const len_t os    = (ncols - i) % UNROLL;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * red[l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        for (; j < ncols; l+=4, j += UNROLL) {
            dr[j]   -=  mul * red[l];
            dr[j+1] -=  mul * red[l+1];
            dr[j+2] -=  mul * red[l+2];
            dr[j+3] -=  mul * red[l+3];
            dr[j]   +=  (dr[j] >> 63) & mod2;
            dr[j+1] +=  (dr[j+1] >> 63) & mod2;
            dr[j+2] +=  (dr[j+2] >> 63) & mod2;
            dr[j+3] +=  (dr[j+3] >> 63) & mod2;
        }
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((uint64_t)(ncols-np), sizeof(cf32_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncols-np, fc);
    }

    *pc = np - ncl;

    return row;
}

static cf32_t *reduce_dense_row_by_dense_new_pivots_31_bit(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    len_t np  = -1;
    const int64_t mod = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);
    const int64_t mod2  = (int64_t)fc * fc;

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

        const int64_t mul = (int64_t)dr[i];
#if defined HAVE_AVX512_F
        const len_t os    = (ncr - i) % 8;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * pivs[i][l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        const __m512i mulv  = _mm512_set1_epi64(mul);
        const __m512i zerov = _mm512_setzero_si512();
        const __m512i mod2v = _mm512_set1_epi64(mod2);
        for (; j < ncr; l += 8, j += 8) {
            const __m512i prodv = _mm512_mul_epu32(mulv, _mm512_cvtepu32_epi64(
                        _mm256_loadu_si256((__m256i *)(pivs[i]+l))));
            __m512i resv  = _mm512_sub_epi64(
                    _mm512_loadu_si512((void *)(dr+j)), prodv);
            const __mmask8 cmpv = _mm512_cmpgt_epi64_mask(zerov, resv);
            resv  = _mm512_mask_add_epi64(resv, cmpv, resv, mod2v);
            _mm512_storeu_si512((void *)(dr+j), resv);
        }
#elif defined HAVE_AVX2
        const len_t os    = (ncr - i) % 8;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * pivs[i][l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        const __m256i mulv  = _mm256_set1_epi32((uint32_t)mul);
        const __m256i zerov = _mm256_setzero_si256();
        const __m256i mod2v = _mm256_set1_epi64x(mod2);
        for (; j < ncr; l += 8, j += 8) {
            __m256i prodv = _mm256_mul_epu32(mulv, _mm256_cvtepu32_epi64(
                        _mm_loadu_si128((__m128i *)(pivs[i]+l))));
            __m256i resv  = _mm256_sub_epi64(
                    _mm256_loadu_si256((__m256i *)(dr+j)), prodv);
            resv  = _mm256_add_epi64(resv, _mm256_and_si256(
                        _mm256_cmpgt_epi64(zerov, resv), mod2v));
            _mm256_storeu_si256((__m256i *)(dr+j), resv);
            prodv = _mm256_mul_epu32(mulv, _mm256_cvtepu32_epi64(
                        _mm_loadu_si128((__m128i *)(pivs[i]+l+4))));
            resv  = _mm256_sub_epi64(
                    _mm256_loadu_si256((__m256i *)(dr+j+4)), prodv);
            resv  = _mm256_add_epi64(resv, _mm256_and_si256(
                        _mm256_cmpgt_epi64(zerov, resv), mod2v));
            _mm256_storeu_si256((__m256i *)(dr+j+4), resv);
        }
#elif defined __aarch64__
        const len_t os    = (ncr - i) % 4;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * pivs[i][l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        const uint32x2_t mulv = vdup_n_u32((uint32_t)mul);
        const int64x2_t mod2v = vdupq_n_s64(mod2);
        for (; j < ncr; l += 4, j += 4) {
            const uint32x4_t cfv  = vld1q_u32(pivs[i]+l);
            int64x2_t resv  = vsubq_s64(vld1q_s64(dr+j),
                    vreinterpretq_s64_u64(vmull_u32(mulv, vget_low_u32(cfv))));
            int64x2_t mask  = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv  = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(dr+j, resv);
            resv  = vsubq_s64(vld1q_s64(dr+j+2),
                    vreinterpretq_s64_u64(vmull_u32(mulv, vget_high_u32(cfv))));
            mask  = vreinterpretq_s64_u64(vcltzq_s64(resv));
            resv  = vaddq_s64(resv, vandq_s64(mask, mod2v));
            vst1q_s64(dr+j+2, resv);
        }
#else
        const len_t os    = (ncr - i) % UNROLL;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * pivs[i][l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        for (; j < ncr; l+=4, j += UNROLL) {
            dr[j]   -=  mul * pivs[i][l];
            dr[j+1] -=  mul * pivs[i][l+1];
            dr[j+2] -=  mul * pivs[i][l+2];
            dr[j+3] -=  mul * pivs[i][l+3];
            dr[j]   +=  (dr[j] >> 63) & mod2;
            dr[j+1] +=  (dr[j+1] >> 63) & mod2;
            dr[j+2] +=  (dr[j+2] >> 63) & mod2;
            dr[j+3] +=  (dr[j+3] >> 63) & mod2;
        }
#endif
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((uint64_t)(ncr-np), sizeof(cf32_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = (int64_t)barrett_reduce((uint64_t)dr[i], mod, bm);
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncr-np, fc);
    }
    *pc = np;

    return row;
}

/* see la_ff_common.c */
#define LA_CF_T cf32_t
#define LA_CF cf_32
#define LA_FF(name) name##_ff_32
/* restrict random multipliers to < 2^15, see
 * probabilistic_sparse_reduced_echelon_form_ff_32() */
#define LA_RANDOM_MULTIPLIER(fc) ((int64_t)rand() & ((1 << 15) - 1))

#define LA_SMALL(name) name##_17_bit
#include "la_ff_small_prime.c"

static void sba_echelon_form_ff_32(
        smat_t *smat,
        crit_t *syz,
        md_t *st,
        const ht_t * const ht
        )
{
    len_t i, j;

    /* row index, might differ from i if we encounter zero reductions */
    len_t ri;

    const len_t nc  = smat->nc;
    const len_t nr  = smat->cld;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs = (hm_t **)calloc((uint64_t)nc, sizeof(hm_t *));

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)nc * sizeof(int64_t));

    /* TODO: At the moment I do not see how to make this parallel due to
     * reduction dependencies on the signatures. */
    for (ri = 0, i = 0; i < nr; ++i) {
        hm_t *npiv      = smat->cr[i];
        cf32_t *cfs     = smat->pc32[npiv[SM_CFS]];
        const hm_t sm   = npiv[SM_SMON];
        const len_t si  = npiv[SM_SIDX];
        const len_t os  = npiv[SM_PRE];
        const len_t len = npiv[SM_LEN];
        const hm_t * const ds = npiv + SM_OFFSET;
        memset(dr, 0, (uint64_t)nc * sizeof(int64_t));
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]    = (int64_t)cfs[j];
            dr[ds[j+1]]  = (int64_t)cfs[j+1];
            dr[ds[j+2]]  = (int64_t)cfs[j+2];
            dr[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        const len_t offset = npiv[SM_OFFSET];
        free(npiv);
        npiv        = NULL;
        smat->cr[i] = NULL;
        /* printf("reducing row %u || %u | %u -- ", i, sm, si);
         * for (int ii = 0; ii < ht->evl; ++ii) {
         *     printf("%u ", ht->ev[sm][ii]);
         * }
         * printf("\n"); */
        npiv = sba_reduce_dense_row_by_known_pivots_sparse_ff_32(
                dr, smat, pivs, offset, sm, si, ri, st);
        if (!npiv) {
            /* row s-reduced to zero, add syzygy and go on with next row */
            add_syzygy_schreyer(syz, sm, si, ht);
            continue;
        }

        ri++;
        /* normalize coefficient array
         * NOTE: this has to be done here, otherwise the reduction may
         * lead to wrong results in a parallel computation since other
         * threads might directly use the new pivot once it is synced. */
        if (smat->cc32[npiv[SM_CFS]][0] != 1) {
            normalize_sparse_matrix_row_ff_32(
                    smat->cc32[npiv[SM_CFS]], npiv[SM_PRE],
                    npiv[SM_LEN], st->fc);
        }
        pivs[npiv[SM_OFFSET]] = npiv;
    }

    /* free initial coefficients coming from previous degree matrix */
    for (i = 0; i < smat->pld; ++i) {
        free(smat->pc32[i]);
        smat->pc32[i] = NULL;
    }
    /* get number of zero reductions and adjust number of
     * rows stored in matrix */
    smat->nz  = smat->cld - ri;
    smat->cld = ri;

    free(pivs);
    pivs = NULL;
    free(dr);
    dr   = NULL;
}

static void exact_sparse_reduced_echelon_form_sat_ff_32(
        bs_t *sat,
        mat_t *mat,
        bs_t *kernel,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    double rt, ct;
    rt = realtime();
    ct = cputime();

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncl   = mat->nru;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    /* saturation elements we have to reduce with the known pivots first */
    hm_t **upivs  = sat->hm;

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)ncols * st->nthrds * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < sat->ld; ++i) {
        int64_t *drl    = dr + (omp_get_thread_num() * (uint64_t)ncols);
        hm_t *npiv      = upivs[i];
        len_t mult      = npiv[MULT];
        /* we only saturate w.r.t. one element at the moment */
        len_t cf_idx    = npiv[COEFFS];
        cf32_t *cfs     = sat->cf_32[cf_idx];
        const len_t os  = npiv[PRELOOP];
        const len_t len = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += UNROLL) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        sc  = 0;
        while (drl[sc] == 0) {
            sc++;
        }
        free(npiv);
        upivs[i]  = NULL;
        free(cfs);
        sat->cf_32[cf_idx]  = NULL;
        /* npiv  = reduce_dense_row_by_known_pivots_sparse_ff_32(
         *         drl, mat, bs, pivs, sc, i, st); */
        npiv  = reduce_dense_row_by_known_pivots_sparse_up_to_ff_31_bit(
                drl, sat, bs, pivs, sc, cf_idx, i, ncl, ncols, st);
        if (!npiv) {
#pragma omp critical
            {
            sat->hm[i]  = NULL;
            kernel->hm[kernel->ld]    = (hm_t *)malloc(
                    (uint64_t)(OFFSET+1) * sizeof(hm_t));
            kernel->cf_32[kernel->ld] = (cf32_t *)malloc(
                    (uint64_t)1 * sizeof(cf32_t));
            kernel->hm[kernel->ld][OFFSET]  = mult;
            kernel->hm[kernel->ld][LENGTH]  = 1;
            kernel->hm[kernel->ld][PRELOOP] = 1;
            kernel->hm[kernel->ld][COEFFS]  = kernel->ld;
            kernel->cf_32[kernel->ld][0]    = 1;
            kernel->ld++;
            }
            continue;
        }
        /* write non zero normal form to sat for reusage in a later
         * saturation step of f4sat */
        npiv[MULT]  = mult;
        sat->hm[i]  = npiv;
    }

    mat->np = mat->nr = mat->sz = nrl;
    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    dr        = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    upivs     = NULL;
    upivs     = (hm_t **)calloc((uint64_t)sat->ld, sizeof(hm_t *));
    len_t ctr = 0;

    if (st->info_level > 1) {
      fprintf(VERBSTREAM, "        normal form time");
    }
    print_sat_nf_round_timings(VERBSTREAM, st, rt, ct);
    /* compute kernel */
    for (i = 0; i < sat->ld; ++i) {
        if (sat->hm[i] != NULL) {
                upivs[ctr++]  = sat->hm[i];
        }
    }
    sort_matrix_rows_mult_increasing(upivs, ctr);
    upivs = realloc(upivs, (uint64_t)ctr * sizeof(hm_t *));

    /* first test if kernel is trivial */
/*     if (is_kernel_trivial(sat, upivs, ncl, mat->ncr, st)) {
 *         [> we do not need the old pivots anymore <]
 *         for (i = 0; i < ncols; ++i) {
 *             free(pivs[i]);
 *             pivs[i] = NULL;
 *         }
 *         free(upivs);
 *         upivs = NULL;
 *         free(pivs);
 *         pivs  = NULL;
 *         free(dr);
 *         dr    = NULL;
 *         return;
 *     }
 *
 *     printf("not trivial! "); */
    /* if kernel is not trivial really compute it */
    /* dense row for multipliers */
    hm_t **mulh     = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    cf32_t **mulcf  = (cf32_t **)calloc((uint64_t)ncols, sizeof(cf32_t *));
    cf32_t **pivcf  = (cf32_t **)calloc((uint64_t)ncols, sizeof(cf32_t *));
    int64_t *drm    = calloc((uint64_t)sat->ld, sizeof(int64_t));
    /* now we do a ususal F4 reduction with updated pivs, but we have
     * to track the reduction steps when reducing each row from upivs */
    for (i = 0; i < ctr; ++i) {
        int64_t *drl          = dr;
        hm_t *npiv            = upivs[i];
        len_t tmp_pos         = npiv[COEFFS];
        cf32_t *cfs           = sat->cf_32[tmp_pos];
        const len_t os        = npiv[PRELOOP];
        const len_t len       = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        k = 0;
        memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));
        memset(drm, 0, (uint64_t)sat->ld * sizeof(int64_t));
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += UNROLL) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        drm[upivs[i][MULT]] = 1;
        do {
            sc    = npiv[OFFSET];
            npiv  = reduce_dense_row_by_known_pivots_sparse_sat_ff_31_bit(
                    drl, drm, pivcf, mulh, mulcf, pivs, sc, tmp_pos,
                    sat->ld, ncols, st, ncl);
            if (!npiv) {
                /* normalize new kernel element */
                normalize_sparse_matrix_row_ff_32(
                        mulcf[mulh[tmp_pos][COEFFS]],
                        mulh[tmp_pos][PRELOOP],
                        mulh[tmp_pos][LENGTH], st->fc);
                kernel->hm[kernel->ld]          = mulh[tmp_pos];
                kernel->cf_32[kernel->ld]       = mulcf[mulh[tmp_pos][COEFFS]];
                mulcf[mulh[tmp_pos][COEFFS]]    = NULL;
                mulh[tmp_pos]                   = NULL;
                kernel->hm[kernel->ld][COEFFS]  = kernel->ld;
                kernel->ld++;
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (pivcf[npiv[COEFFS]][0] != 1) {
                /* adjust mul entry correspondingly */
                adjust_multiplier_sparse_matrix_row_ff_32(
                        mulcf[mulh[npiv[MULT]][COEFFS]],
                        pivcf[npiv[COEFFS]],
                        mulh[npiv[MULT]][PRELOOP],
                        mulh[npiv[MULT]][LENGTH],
                        st->fc);
                /* for (int kk = 0; kk < mulh[npiv[MULT]][LENGTH]; ++kk) {
                 *     printf("%u ", mulcf[mulh[npiv[MULT]][COEFFS]][kk]);
                 * }
                 * printf("\n"); */
                normalize_sparse_matrix_row_ff_32(
                        pivcf[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);

            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
        } while (!k);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncols; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }
    for (i = 0; i < ncols; ++i) {
        free(mulcf[i]);
        mulcf[i]  = NULL;
        free(pivcf[i]);
        pivcf[i]  = NULL;
        free(mulh[i]);
        mulh[i]   = NULL;
    }

    free(pivcf);
    pivcf = NULL;
    free(mulcf);
    mulcf = NULL;
    free(mulh);
    mulh = NULL;
    free(upivs);
    upivs = NULL;
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr    = NULL;
    free(drm);
    drm   = NULL;
}

static void sba_linear_algebra_ff_32(
        smat_t *smat,
        crit_t *syz,
        md_t *st,
        const ht_t * const ht
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    smat->cc32 = realloc(smat->cc32,
            (uint64_t)smat->cld * sizeof(cf32_t *));

    sba_echelon_form_ff_32(smat, syz, st, ht);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += smat->nz;
}

static void copy_kernel_to_matrix(
        mat_t *mat,
        bs_t *kernel,
        const len_t nc
        )
{
    len_t i;

    sort_matrix_rows_increasing(kernel->hm, kernel->ld);

    mat->tr = (hm_t **)malloc((uint64_t)kernel->ld * sizeof(hm_t *));

    mat->nr   = kernel->ld;
    mat->nc   = nc;
    mat->nru  = 0;
    mat->nrl  = kernel->ld;
    mat->ncl  = 0;
    mat->ncr  = mat->nc;

    for (i = 0; i < kernel->ld; ++i) {
        mat->tr[i]  = kernel->hm[i];
    }
}

static void compute_kernel_sat_ff_32(
        bs_t *sat,
        mat_t *mat,
        bs_t *kernel,
        const bs_t * const bs,
        md_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    check_enlarge_basis(kernel, sat->ld, st);

    exact_sparse_reduced_echelon_form_sat_ff_32(
            sat, mat, kernel, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;
}

#include "la_ff_common.c"
