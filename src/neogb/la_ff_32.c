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

/* That's also enough if AVX512 is avaialable on the system */
#if defined HAVE_AVX2
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
    const int64_t mod2  = (int64_t)fc * fc;

    for (i = start; i < nc; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
                dr[j] = dr[j] % mod;
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
    printf("sat->ld %u | max_col %u | ncl %u || %u\n", sat->ld, max_col, ncl, max_col - ncl);
    int ret = 0;
    const int64_t fc   = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;
    /* make a dense matrix first */
    const len_t dim = max_col - ncl;
    /* allocate memory */
    cf32_t *dm    = (cf32_t *)calloc((unsigned long)(sat->ld*dim), sizeof(cf32_t));
    int64_t *dr   = (int64_t *)calloc((unsigned long)dim, sizeof(int64_t));
    int64_t *mul  = (int64_t *)calloc((unsigned long)ncr, sizeof(int64_t));
    cf32_t **pivs = (cf32_t **)calloc((unsigned long)dim, sizeof(cf32_t *));

    /* fill random value array */
    for (i = 0; i < ncr; ++i) {
        mul[i]  = (int64_t)rand() & fc;
    }
    for (i = 0; i < sat->ld; ++i) {
        memset(dr, 0, (unsigned long)dim * sizeof(int64_t));
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
            dr[i] = dr[i] % mod;
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
            cfs   = bs->cf_32[dts[COEFFS]];
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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

static hm_t *trace_reduce_dense_row_by_known_pivots_sparse_17_bit(
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
    cf32_t *cfs;
    int64_t np = -1;
    const int64_t mod           = (int64_t)st->fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
            cfs   = bs->cf_32[dts[COEFFS]];
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
        st->trace_nr_mult +=  len / 1000.0;
        st->trace_nr_add  +=  len / 1000.0;
        st->trace_nr_red++;
    }
    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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
    const int64_t mod2  = (int64_t)st->fc * st->fc;
#ifdef HAVE_AVX2
    int64_t res[4] __attribute__((aligned(32)));
    __m256i cmpv, redv, drv, mulv, prodv, resv, rresv;
    __m256i zerov= _mm256_set1_epi64x(0);
    __m256i mod2v = _mm256_set1_epi64x(mod2);
#endif

    for (i = dpiv; i < end; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }

    k = ncols-end;
    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = end; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
    row = realloc(row, (unsigned long)(j+OFFSET) * sizeof(hm_t));
    cf  = realloc(cf, (unsigned long)j * sizeof(cf32_t));
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
            dr[i] = dr[i] % mod;
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }

    if (k == 0) {
        free(smat->cr[ri]);
        smat->cr[ri] = NULL;

        return smat->cr[ri];
    }

    /* printf("k %d | ri %u | cr[%u] = %p\n", k, ri, ri, smat->cr[ri]);
     * printf("%lu\n", (unsigned long)(k+SM_OFFSET)); */
    smat->cr[ri] = realloc(smat->cr[ri],
            (unsigned long)(k+SM_OFFSET) * sizeof(hm_t));
    cf32_t *cf   = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));

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
            dr[i] = dr[i] % mod;
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
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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
            dr[i] = dr[i] % mod;
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }

    /* we always have to write the multiplier entry for the kernel
     * construction, especially if we computed a zero row. */
    mulh[tmp_pos]   =
        (hm_t *)malloc((unsigned long)(sat_ld+OFFSET) * sizeof(hm_t));
    mulcf[tmp_pos]  =
        (cf32_t *)malloc((unsigned long)(sat_ld) * sizeof(cf32_t));
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
         (unsigned long)(j+OFFSET) * sizeof(hm_t));
    mulcf[tmp_pos]  = realloc(mulcf[tmp_pos],
            (unsigned long)j * sizeof(cf32_t));


    if (k == 0) {
        return NULL;
    }

    k = ncols - np;
    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + OFFSET;
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
    const int64_t mod           = (int64_t)st->fc;
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
            dr[i] = dr[i] % mod;
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
        st->trace_nr_mult +=  len / 1000.0;
        st->trace_nr_add  +=  len / 1000.0;
        st->trace_nr_red++;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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
        st->application_nr_mult +=  len / 1000.0;
        st->application_nr_add  +=  len / 1000.0;
        st->application_nr_red++;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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

    hm_t *row   = (hm_t *)malloc((unsigned long)(k+OFFSET) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(k) * sizeof(cf32_t));
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

static cf32_t *reduce_dense_row_by_all_pivots_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const * const pivs,
        cf32_t *const * const dpivs,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    const int64_t mod = (int64_t)fc;
    len_t np  = -1;
    cf32_t *red;

    const len_t ncl   = mat->ncl;
    const len_t ncols = mat->nc;

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        const cf32_t *cfs = bs->cf_32[pivs[i][COEFFS]];
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
            dr[i] = dr[i] % mod;
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

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncols-np), sizeof(cf32_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncols-np, fc);
    }

    *pc = np - ncl;
    return row;
}

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
    const int64_t mod2  = (int64_t)fc * fc;
    len_t np  = -1;
    cf32_t *red;

    const len_t ncl   = mat->ncl;
    const len_t ncols = mat->nc;

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
            dr[i] = dr[i] % mod;
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

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncols-np), sizeof(cf32_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncols-np, fc);
    }

    *pc = np - ncl;

    return row;
}

static cf32_t *reduce_dense_row_by_old_pivots_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const uint32_t fc
        )
{
    hi_t i, j;
    const int64_t mod = (int64_t)fc;
    const len_t ncols = mat->nc;
    const len_t ncl   = mat->ncl;
    const len_t ncr   = mat->ncr;

    for (i = dpiv; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        const cf32_t *cfs = bs->cf_32[pivs[i][COEFFS]];
        const len_t os    = pivs[i][PRELOOP];
        const len_t len   = pivs[i][LENGTH];
        const hm_t * const ds = pivs[i] + OFFSET;
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

    /* store a dense row for further dense gaussian elimination */
    cf32_t *row  = (cf32_t *)calloc(
            (unsigned long)(ncr), sizeof(cf32_t));
    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf32_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }
    return row;
}

static cf32_t *reduce_dense_row_by_old_pivots_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,  /* pivot of dense row at the beginning */
        const uint32_t fc
        )
{
    hi_t i, j;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;
    const len_t ncols   = mat->nc;
    const len_t ncl     = mat->ncl;
    const len_t ncr     = mat->ncr;

    for (i = dpiv; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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

    /* store a dense row for further dense gaussian elimination */
    cf32_t *row  = (cf32_t *)calloc(
            (unsigned long)(ncr), sizeof(cf32_t));

    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf32_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }

    return row;
}

static cf32_t *reduce_dense_row_by_dense_new_pivots_17_bit(
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

    for (k = 0, i = *pc; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncr-np), sizeof(cf32_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncr-np, fc);
    }
    *pc = np;
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
    const int64_t mod2  = (int64_t)fc * fc;

    for (k = 0, i = *pc; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
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
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncr-np), sizeof(cf32_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncr-np, fc);
    }
    *pc = np;

    return row;
}

static void probabilistic_sparse_reduced_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k, l, m;

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    const uint32_t fc   = st->fc;
    /* Why do we generate the random linear combinations so strangely
     * compared to, e.g. la_ff_8.c or la_ff_16.c?
     * We can have fc > 2^31, so shifting the result of random multiplier
     * and coefficient may lead to elements of size > 2^63 which cannot
     * be represented by int64_t anymore. Thus we restrict the multiplier
     * via masking with mask to < 2^31. Then, after subtracting the
     * multiplied values we have to check if at some point the result gets
     * negative as int64_t, i.e. if drl[i] >> 63 is 1. If so, we add a
     * corresponding multiple of fc. Again, if fc > 2^31 we cannot just
     * add fc * fc, but we have to restrict to fc/2 * fc.
     * Note: We cannot apply this trick to the reduction itself, there we
     * cannot restrict the multiplier to be < 2^31, so we have to handle the
     * fc > 2^31 case differently. */
    uint64_t tmp  = (uint64_t)fc * fc;
    while (tmp > pow(2,63)) {
      tmp -=  (uint64_t)(fc/2) * fc;
    }
    const int64_t mod2  = tmp;
    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;
    /* const int64_t mask  = pow(2,(uint32_t)(ceil(log((double)st->max_uht_size)/log(2))))-1; */
    const int64_t mask  = pow(2,15)-1;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) \
    schedule(dynamic)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            hm_t *npiv  = NULL;
            cf32_t *cfs;
            /* starting column, offset, coefficient array position in tmpcf */
            hm_t sc, cfp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                cfp = bctr + i*rpb;
                sc  = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() & mask;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    cfs   = bs->cf_32[npiv[COEFFS]];
                    const len_t os  = npiv[PRELOOP];
                    const len_t len = npiv[LENGTH];
                    const hm_t * const ds = npiv + OFFSET;
                    sc    = sc < ds[0] ? sc : ds[0];
                    for (l = 0; l < os; ++l) {
                        drl[ds[l]]  -=  mull[k] * cfs[l];
                        drl[ds[l]]  +=  (drl[ds[l]] >> 63) & mod2;
                    }
                    for (; l < len; l += UNROLL) {
                        drl[ds[l]]    -=  mull[k] * cfs[l];
                        drl[ds[l]]    +=  (drl[ds[l]] >> 63) & mod2;
                        drl[ds[l+1]]  -=  mull[k] * cfs[l+1];
                        drl[ds[l+1]]  +=  (drl[ds[l+1]] >> 63) & mod2;
                        drl[ds[l+2]]  -=  mull[k] * cfs[l+2];
                        drl[ds[l+2]]  +=  (drl[ds[l+2]] >> 63) & mod2;
                        drl[ds[l+3]]  -=  mull[k] * cfs[l+3];
                        drl[ds[l+3]]  +=  (drl[ds[l+3]] >> 63) & mod2;
                    }
                }
                k     = 0;
                cfs   = NULL;
                npiv  = NULL;
                /* do the reduction */
                do {
                    free(cfs);
                    cfs = NULL;
                    free(npiv);
                    npiv  = NULL;
                    npiv  = reduce_dense_row_by_known_pivots_sparse_ff_32(
                            drl, mat, bs, pivs, sc, cfp, 0, 0, 0, st);
                    if (!npiv) {
                        bctr  = nrbl;
                        break;
                    }
                    /* normalize coefficient array
                    * NOTE: this has to be done here, otherwise the reduction may
                    * lead to wrong results in a parallel computation since other
                    * threads might directly use the new pivot once it is synced. */
                    if (mat->cf_32[npiv[COEFFS]][0] != 1) {
                        normalize_sparse_matrix_row_ff_32(
                                mat->cf_32[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    cfs = mat->cf_32[npiv[COEFFS]];
                    sc  = npiv[OFFSET];
                    k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(upivs[j]);
                upivs[j]  = NULL;
            }
        }
    }
    free(mul);
    mul   = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (unsigned long)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf32_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc, cfp;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = mat->cf_32[pivs[k][COEFFS]];
            cfp = pivs[k][COEFFS];
            const len_t bi  = pivs[k][BINDEX];
            const len_t mh  = pivs[k][MULT];
            const len_t os  = pivs[k][PRELOOP];
            const len_t len = pivs[k][LENGTH];
            const hm_t * const ds = pivs[k] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                dr[ds[j]]   = (int64_t)cfs[j];
                dr[ds[j+1]] = (int64_t)cfs[j+1];
                dr[ds[j+2]] = (int64_t)cfs[j+2];
                dr[ds[j+3]] = (int64_t)cfs[j+3];
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cfp, mh, bi, 0, st);
        }
    }
    free(mat->rr);
    mat->rr = NULL;

    free(pivs);
    pivs  = NULL;

    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}

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
    hm_t **pivs = (hm_t **)calloc((unsigned long)nc, sizeof(hm_t *));

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)nc * sizeof(int64_t));

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
        memset(dr, 0, (unsigned long)nc * sizeof(int64_t));
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

static void exact_sparse_reduced_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    len_t bad_prime = 0;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        if (bad_prime == 0) {
            int64_t *drl  = dr + (omp_get_thread_num() * ncols);
            hm_t *npiv      = upivs[i];
            cf32_t *cfs     = tbr->cf_32[npiv[COEFFS]];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            const hm_t * const ds = npiv + OFFSET;
            k = 0;
            memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
            for (j = 0; j < os; ++j) {
                drl[ds[j]]  = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                drl[ds[j]]    = (int64_t)cfs[j];
                drl[ds[j+1]]  = (int64_t)cfs[j+1];
                drl[ds[j+2]]  = (int64_t)cfs[j+2];
                drl[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            cfs = NULL;
            do {
                sc  = npiv[OFFSET];
                free(npiv);
                free(cfs);
                npiv  = mat->tr[i] = reduce_dense_row_by_known_pivots_sparse_ff_32(
                        drl, mat, bs, pivs, sc, i, mh, bi, st->trace_level == LEARN_TRACER, st);
                if (st->nf > 0) {
                    if (!npiv) {
                        mat->tr[i]  = NULL;
                        break;
                    }
                    mat->tr[i]  = npiv;
                    cfs = mat->cf_32[npiv[COEFFS]];
                    break;
                } else {
                    if (!npiv) {
                        if (st->trace_level == APPLY_TRACER) {
                            bad_prime = 1;
                        }
                        break;
                    }
                    /* normalize coefficient array
                     * NOTE: this has to be done here, otherwise the reduction may
                     * lead to wrong results in a parallel computation since other
                     * threads might directly use the new pivot once it is synced. */
                    if (mat->cf_32[npiv[COEFFS]][0] != 1) {
                        normalize_sparse_matrix_row_ff_32(
                                mat->cf_32[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                    cfs = mat->cf_32[npiv[COEFFS]];
                }
            } while (!k);
        }
    }

    if (bad_prime == 1) {
        for (i = 0; i < ncl+ncr; ++i) {
            free(pivs[i]);
            pivs[i] = NULL;
        }
        mat->np = 0;
        if (st->info_level > 0) {
            fprintf(stderr, "Zero reduction while applying tracer, bad prime.\n");
        }
        return;
    }

    /* construct the trace */
    if (st->trace_level == LEARN_TRACER) {
        construct_trace(st->tr, mat);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    if (st->nf == 0) {
        dr      = realloc(dr, (unsigned long)ncols * sizeof(int64_t));
        mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

        /* interreduce new pivots */
        cf32_t *cfs;
        hm_t cf_array_pos;
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
                cfs = mat->cf_32[pivs[k][COEFFS]];
                cf_array_pos    = pivs[k][COEFFS];
                const len_t os  = pivs[k][PRELOOP];
                const len_t len = pivs[k][LENGTH];
                const len_t bi  = pivs[k][BINDEX];
                const len_t mh  = pivs[k][MULT];
                const hm_t * const ds = pivs[k] + OFFSET;
                sc  = ds[0];
                for (j = 0; j < os; ++j) {
                    dr[ds[j]] = (int64_t)cfs[j];
                }
                for (; j < len; j += UNROLL) {
                    dr[ds[j]]    = (int64_t)cfs[j];
                    dr[ds[j+1]]  = (int64_t)cfs[j+1];
                    dr[ds[j+2]]  = (int64_t)cfs[j+2];
                    dr[ds[j+3]]  = (int64_t)cfs[j+3];
                }
                free(pivs[k]);
                free(cfs);
                pivs[k] = NULL;
                pivs[k] = mat->tr[npivs++] =
                    reduce_dense_row_by_known_pivots_sparse_ff_32(
                            dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
            }
        }
        mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
        st->np = mat->np = mat->nr = mat->sz = npivs;
    } else {
        st->np = mat->np = mat->nr = mat->sz = nrl;
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;
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
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* saturation elements we have to reduce with the known pivots first */
    hm_t **upivs  = sat->hm;

    /* temporary storage for elements to be saturated */
    /* bs_t *tmp = initialize_basis(sat->ld); */
    /* temporary storage for corresponding multipliers
     * which then generate the kernel elements */
    /* bs_t *mul = initialize_basis(sat->ld); */

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < sat->ld; ++i) {
        int64_t *drl    = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv      = upivs[i];
        len_t mult      = npiv[MULT];
        /* we only saturate w.r.t. one element at the moment */
        len_t cf_idx    = npiv[COEFFS];
        cf32_t *cfs     = sat->cf_32[cf_idx];
        const len_t os  = npiv[PRELOOP];
        const len_t len = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
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
                    (unsigned long)(OFFSET+1) * sizeof(hm_t));
            kernel->cf_32[kernel->ld] = (cf32_t *)malloc(
                    (unsigned long)1 * sizeof(cf32_t));
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

    dr        = realloc(dr, (unsigned long)ncols * sizeof(int64_t));
    upivs     = NULL;
    upivs     = (hm_t **)calloc((unsigned long)sat->ld, sizeof(hm_t *));
    len_t ctr = 0;

    if (st->info_level > 1) {
      printf("        normal form time");
    }
    print_sat_nf_round_timings(stdout, st, rt, ct);
    /* compute kernel */
    for (i = 0; i < sat->ld; ++i) {
        if (sat->hm[i] != NULL) {
                upivs[ctr++]  = sat->hm[i];
        }
    }
    sort_matrix_rows_mult_increasing(upivs, ctr);
    upivs = realloc(upivs, (unsigned long)ctr * sizeof(hm_t *));

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
    hm_t **mulh     = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    cf32_t **mulcf  = (cf32_t **)calloc((unsigned long)ncols, sizeof(cf32_t *));
    cf32_t **pivcf  = (cf32_t **)calloc((unsigned long)ncols, sizeof(cf32_t *));
    int64_t *drm    = calloc((unsigned long)sat->ld, sizeof(int64_t));
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
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        memset(drm, 0, (unsigned long)sat->ld * sizeof(int64_t));
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

static void exact_trace_sparse_reduced_echelon_form_ff_32(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv      = upivs[i];
        rba_t *rba      = mat->rba[i];
        cf32_t *cfs     = bs->cf_32[npiv[COEFFS]];
        const len_t bi  = npiv[BINDEX];
        const len_t mh  = npiv[MULT];
        const len_t os  = npiv[PRELOOP];
        const len_t len = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        k = 0;
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += UNROLL) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        cfs = NULL;
        do {
            sc  = npiv[OFFSET];
            free(npiv);
            free(cfs);
            npiv  = mat->tr[i]  = trace_reduce_dense_row_by_known_pivots_sparse_ff_32(
                    rba, drl, mat, bs, pivs, sc, i, mh, bi, st);
            if (!npiv) {
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mat->cf_32[npiv[COEFFS]][0] != 1) {
                normalize_sparse_matrix_row_ff_32(
                        mat->cf_32[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                st->trace_nr_mult  +=  npiv[LENGTH] / 1000.0;
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
            cfs = mat->cf_32[npiv[COEFFS]];
        } while (!k);
    }

    /* construct the trace */
    construct_trace(trace, mat);

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (unsigned long)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf32_t *cfs;
    hm_t cf_array_pos;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = mat->cf_32[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
            const len_t os  = pivs[k][PRELOOP];
            const len_t len = pivs[k][LENGTH];
            const len_t bi  = pivs[k][BINDEX];
            const len_t mh  = pivs[k][MULT];
            const hm_t * const ds = pivs[k] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                dr[ds[j]]    = (int64_t)cfs[j];
                dr[ds[j+1]]  = (int64_t)cfs[j+1];
                dr[ds[j+2]]  = (int64_t)cfs[j+2];
                dr[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}

static int exact_application_sparse_reduced_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
    int flag  = 1;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        if (flag == 1) {
            int64_t *drl    = dr + (omp_get_thread_num() * ncols);
            hm_t *npiv      = upivs[i];
            cf32_t *cfs     = bs->cf_32[npiv[COEFFS]];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            const hm_t * const ds = npiv + OFFSET;
            k = 0;
            memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
            for (j = 0; j < os; ++j) {
                drl[ds[j]]  = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                drl[ds[j]]    = (int64_t)cfs[j];
                drl[ds[j+1]]  = (int64_t)cfs[j+1];
                drl[ds[j+2]]  = (int64_t)cfs[j+2];
                drl[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            cfs = NULL;
            do {
                sc  = npiv[OFFSET];
                free(npiv);
                free(cfs);
                npiv  = mat->tr[i]  = reduce_dense_row_by_known_pivots_sparse_ff_32(
                        drl, mat, bs, pivs, sc, i, mh, bi, 0, st);
                if (!npiv) {
                    fprintf(stderr, "Unlucky prime detected, row reduced to zero.");
                    flag  = 0;
                    break;
                }

                /* normalize coefficient array
                 * NOTE: this has to be done here, otherwise the reduction may
                 * lead to wrong results in a parallel computation since other
                 * threads might directly use the new pivot once it is synced. */
                if (mat->cf_32[npiv[COEFFS]][0] != 1) {
                    normalize_sparse_matrix_row_ff_32(
                            mat->cf_32[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    st->application_nr_mult  +=  npiv[LENGTH] / 1000.0;
                }
                k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                cfs = mat->cf_32[npiv[COEFFS]];
            } while (!k);
        }
    }
    /* unlucky prime found */
    if (flag == 0) {
        return 1;
    }
    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (unsigned long)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf32_t *cfs;
    hm_t cf_array_pos;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = mat->cf_32[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
            const len_t os  = pivs[k][PRELOOP];
            const len_t len = pivs[k][LENGTH];
            const len_t bi  = pivs[k][BINDEX];
            const len_t mh  = pivs[k][MULT];
            const hm_t * const ds = pivs[k] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                dr[ds[j]]    = (int64_t)cfs[j];
                dr[ds[j+1]]  = (int64_t)cfs[j+1];
                dr[ds[j+2]]  = (int64_t)cfs[j+2];
                dr[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;

    return 0;
}

static cf32_t **sparse_AB_CD_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * bs,
        md_t *st
        )
{
    len_t i = 0, j;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncl   = mat->ncl;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* dense rows representing updated D part;
     * after reducing CD part with AB */
    cf32_t **drs    = (cf32_t **)calloc((unsigned long)nrl, sizeof(cf32_t *));

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        cf32_t *cfs = NULL;
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv    = upivs[i];
        /* do the reduction */
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        cfs = bs->cf_32[npiv[COEFFS]];
        const len_t os  = npiv[PRELOOP];
        const len_t len = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += UNROLL) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        sc  = ds[0];
        free(npiv);
        drs[i]  = reduce_dense_row_by_old_pivots_ff_32(
                drl, mat, bs, pivs, sc, st->fc);
    }
    free(dr);
    dr  = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }
    free(pivs);
    pivs  = NULL;

    /* remove NULL dense rows */
    len_t npivs = 0; /* number of new pivots */
    for (i = 0; i < nrl; ++i) {
        if (drs[i] != NULL) {
            drs[npivs++]  = drs[i];
        }
    }
    if (npivs == 0) {
        free(drs);
        drs = NULL;
    }
    st->np = mat->np = npivs;

    return drs;
}

static cf32_t **interreduce_dense_matrix_ff_32(
    cf32_t **dm,
    const len_t ncr,
    const uint32_t fc
    )
{
    len_t i, j, k, l;
    int64_t *dr = malloc((unsigned long)ncr * sizeof(int64_t));

    for (i = 0; i < ncr; ++i) {
        k = ncr-1-i;
        if (dm[k]) {
            memset(dr, 0, (unsigned long)ncr * sizeof(int64_t));
            const len_t npc  = ncr - k;
            const len_t os   = npc % UNROLL;
            for (j = k, l = 0; l < os; ++j, ++l) {
                dr[j] = (int64_t)dm[k][l];
            }
            for (; l < npc; j += UNROLL, l += UNROLL) {
                dr[j]   = (int64_t)dm[k][l];
                dr[j+1] = (int64_t)dm[k][l+1];
                dr[j+2] = (int64_t)dm[k][l+2];
                dr[j+3] = (int64_t)dm[k][l+3];
            }
            free(dm[k]);
            dm[k] = NULL;
            /* start with previous pivot the reduction process, so keep the
             * pivot element as it is */
            dm[k] = reduce_dense_row_by_dense_new_pivots_ff_32(
                    dr, &k, dm, ncr, fc);
        }
    }
    free(dr);
    return dm;
}

static cf32_t **exact_dense_linear_algebra_ff_32(
        cf32_t **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, l, npivs;

    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));
    /* rows to be further reduced */
    cf32_t **tbr  = (cf32_t **)calloc((unsigned long)nrows, sizeof(cf32_t *));
    int64_t *dr   = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncr) * sizeof(int64_t));

    /* separate rows already representing new pivots and rows to
     * be further reduced by these new pivots */
    j     = 0;
    npivs = 0;
    for (i = 0; i < nrows; ++i) {
        if (dm[i] != NULL) {
            k = 0;
            while (dm[i][k] == 0) {
                ++k;
            }
            if (nps[k] == NULL) {
                /* we have a pivot, cut the dense row down to start
                 * at the first nonzero entry */
                memmove(dm[i], dm[i]+k, (unsigned long)(ncr-k) * sizeof(cf32_t));
                dm[i] = realloc(dm[i], (unsigned long)(ncr-k) * sizeof(cf32_t));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_32(nps[k], ncr-k, st->fc);
                }
            } else {
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (unsigned long)ntr * sizeof(cf32_t *));
    /* offset modulo 4 for loop unrolling, +1 due to storing the first
     * nonzero entry at the first position */

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(dynamic)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncr);
        memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));
        hm_t npc  = 0;
        hm_t os   = 0;
        cf32_t *npiv  = tbr[i];
        os   = (ncr-npc) % UNROLL;
        for (l = 0, j = npc; l < os; ++l, ++j) {
            drl[j]  = (int64_t)npiv[l];
        }
        for (; j < ncr; l += UNROLL, j += UNROLL) {
            drl[j]    = (int64_t)npiv[l];
            drl[j+1]  = (int64_t)npiv[l+1];
            drl[j+2]  = (int64_t)npiv[l+2];
            drl[j+3]  = (int64_t)npiv[l+3];
        }
        do {
            free(npiv);
            npiv = NULL;
            npiv = reduce_dense_row_by_dense_new_pivots_ff_32(
                    drl, &npc, nps, mat->ncr, st->fc);
            if (npc == -1) {
                break;
            }
            k = __sync_bool_compare_and_swap(&nps[npc], NULL, npiv);
            /* some other thread has already added a pivot so we have to
             * recall the dense reduction process */
        } while (!k);
    }
    /* count number of pivots */
    const len_t os  = ncr % UNROLL;
    for (i = 0; i < os; ++i) {
        if (nps[i] != NULL) {
            npivs++;
        }
    }
    for (; i < ncr; i += UNROLL) {
        if (nps[i] != NULL) {
            npivs++;
        }
        if (nps[i+1] != NULL) {
            npivs++;
        }
        if (nps[i+2] != NULL) {
            npivs++;
        }
        if (nps[i+3] != NULL) {
            npivs++;
        }
    }
    st->np = mat->np = npivs;

    free(tbr);
    free(dr);

    return nps;
}

static cf32_t **probabilistic_dense_linear_algebra_ff_32(
        cf32_t **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, l, m, npivs;

    const uint32_t fc = st->fc;
    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncols = mat->nc;
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));
    /* rows to be further reduced */
    cf32_t **tbr  = (cf32_t **)calloc((unsigned long)nrows, sizeof(cf32_t *));

    /* separate rows already representing new pivots and rows to
     * be further reduced by these new pivots */
    j     = 0;
    npivs = 0;
    for (i = 0; i < nrows; ++i) {
        if (dm[i] != NULL) {
            k = 0;
            while (dm[i][k] == 0) {
                ++k;
            }
            if (nps[k] == NULL) {
                /* we have a pivot, cut the dense row down to start
                 * at the first nonzero entry */
                memmove(dm[i], dm[i]+k, (unsigned long)(ncr-k) * sizeof(cf32_t));
                dm[i] = realloc(dm[i], (unsigned long)(ncr-k) * sizeof(cf32_t));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_32(nps[k], ncr-k, st->fc);
                }
            } else {
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (unsigned long)ntr * sizeof(cf32_t *));

    /* Why do we generate the random linear combinations so strangely
     * compared to, e.g. la_ff_8.c or la_ff_16.c?
     * We can have fc > 2^31, so shifting the result of random multiplier
     * and coefficient may lead to elements of size > 2^63 which cannot
     * be represented by int64_t anymore. Thus we restrict the multiplier
     * via masking with mask to < 2^31. Then, after subtracting the
     * multiplied values we have to check if at some point the result gets
     * negative as int64_t, i.e. if drl[i] >> 63 is 1. If so, we add a
     * corresponding multiple of fc. Again, if fc > 2^31 we cannot just
     * add fc * fc, but we have to restrict to fc/2 * fc.
     * Note: We cannot apply this trick to the reduction itself, there we
     * cannot restrict the multiplier to be < 2^31, so we have to handle the
     * fc > 2^31 case differently. */
    uint64_t tmp  = (uint64_t)fc * fc;
    while (tmp > pow(2,63)) {
      tmp -=  (uint64_t)(fc/2) * fc;
    }
    const int64_t mod2  = tmp;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(ntr/3)))+1;
    const len_t rem = (ntr % nb == 0) ? 0 : 1;
    const len_t rpb = (ntr / nb) + rem;
    /* const int64_t mask  = pow(2,(uint32_t)(ceil(log((double)st->max_uht_size)/log(2))))-1; */
    const int64_t mask  = pow(2,15)-1;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(dynamic)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncr);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (ntr > (i+1)*rpb ? (i+1)*rpb : ntr);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);

        if (nrbl > 0) {
            hm_t npc;
            hm_t os;
            cf32_t *tmp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;
                os  = ncr % UNROLL;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() & mask;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    tmp  = tbr[m];
                    for (l = 0, j = npc; l < os; ++l, ++j) {
                        drl[j]  -=  mull[k] * tmp[l];
                        drl[j]  +=  (drl[j] >> 63) & mod2;
                    }
                    for (; j < ncr; l += UNROLL, j += UNROLL) {
                        drl[j]    -=  mull[k] * tmp[l];
                        drl[j]    +=  (drl[j] >> 63) & mod2;
                        drl[j+1]  -=  mull[k] * tmp[l+1];
                        drl[j+1]  +=  (drl[j+1] >> 63) & mod2;
                        drl[j+2]  -=  mull[k] * tmp[l+2];
                        drl[j+2]  +=  (drl[j+2] >> 63) & mod2;
                        drl[j+3]  -=  mull[k] * tmp[l+3];
                        drl[j+3]  +=  (drl[j+3] >> 63) & mod2;
                    }
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = reduce_dense_row_by_dense_new_pivots_ff_32(
                            drl, &npc, nps, mat->ncr, st->fc);
                    if (npc == -1) {
                        bctr  = nrbl;
                        break;
                    }
                    k = __sync_bool_compare_and_swap(&nps[npc], NULL, tmp);
                    /* some other thread has already added a pivot so we have to
                    * recall the dense reduction process */
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(tbr[j]);
                tbr[j]  = NULL;
            }
        }
    }
    /* count number of pivots */
    const len_t os  = ncr % UNROLL;
    for (i = 0; i < os; ++i) {
        if (nps[i] != NULL) {
            npivs++;
        }
    }
    for (; i < ncr; i += UNROLL) {
        if (nps[i] != NULL) {
            npivs++;
        }
        if (nps[i+1] != NULL) {
            npivs++;
        }
        if (nps[i+2] != NULL) {
            npivs++;
        }
        if (nps[i+3] != NULL) {
            npivs++;
        }
    }
    st->np = mat->np = npivs;

    free(mul);
    free(tbr);
    free(dr);

    return nps;
}

static cf32_t **probabilistic_sparse_dense_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k, l, m, npivs;

    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncols = mat->nc;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));

    const uint32_t fc   = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) shared(nps) \
    schedule(dynamic)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl > 0) {
            hm_t *npiv;
            cf32_t *tmp;
            hm_t npc;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    tmp   = bs->cf_32[npiv[COEFFS]];
                    const len_t os  = npiv[PRELOOP];
                    const len_t len = npiv[LENGTH];
                    const hm_t * const ds = npiv + OFFSET;
                    for (l = 0; l < os; ++l) {
                        drl[ds[l]]  -=  mull[k] * tmp[l];
                        drl[ds[l]]  +=  (drl[ds[l]] >> 63) & mod2;
                    }
                    for (; l < len; l += UNROLL) {
                        drl[ds[l]]    -=  mull[k] * tmp[l];
                        drl[ds[l]]    +=  (drl[ds[l]] >> 63) & mod2;
                        drl[ds[l+1]]  -=  mull[k] * tmp[l+1];
                        drl[ds[l+1]]  +=  (drl[ds[l+1]] >> 63) & mod2;
                        drl[ds[l+2]]  -=  mull[k] * tmp[l+2];
                        drl[ds[l+2]]  +=  (drl[ds[l+2]] >> 63) & mod2;
                        drl[ds[l+3]]  -=  mull[k] * tmp[l+3];
                        drl[ds[l+3]]  +=  (drl[ds[l+3]] >> 63) & mod2;
                    }
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = reduce_dense_row_by_all_pivots_ff_32(
                            drl, mat, bs, &npc, pivs, nps, st->fc);
                    if (npc == -1) {
                        bctr  = nrbl;
                        break;
                    }
                    k = __sync_bool_compare_and_swap(&nps[npc], NULL, tmp);
                    /* some other thread has already added a pivot so we have to
                    * recall the dense reduction process */
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(upivs[j]);
                upivs[j]  = NULL;
            }
        }
    }
    npivs = 0;
    /* count number of pivots */
    const len_t os  = ncr % UNROLL;
    for (i = 0; i < os; ++i) {
        if (nps[i] != NULL) {
            npivs++;
        }
    }
    for (; i < ncr; i += UNROLL) {
        if (nps[i] != NULL) {
            npivs++;
        }
        if (nps[i+1] != NULL) {
            npivs++;
        }
        if (nps[i+2] != NULL) {
            npivs++;
        }
        if (nps[i+3] != NULL) {
            npivs++;
        }
    }
    st->np = mat->np = npivs;


    for (i = 0; i < nru; ++i) {
        free(pivs[i]);
    }
    free(pivs);
    pivs  = NULL;
    free(mul);
    mul = NULL;
    free(dr);
    dr  = NULL;

    return nps;
}

static void convert_to_sparse_matrix_rows_ff_32(
        mat_t *mat,
        cf32_t * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, j, k, l, m;
    cf32_t *cfs;
    hm_t *dts, *dss;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->tr     = realloc(mat->tr, (unsigned long)mat->np * sizeof(hm_t *));
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)mat->np * sizeof(cf32_t *));

    l = 0;
    for (i = 0; i < ncr; ++i) {
        m = ncr-1-i;
        if (dm[m] != NULL) {
            cfs = malloc((unsigned long)(ncr-m) * sizeof(cf32_t));
            dts = malloc((unsigned long)(ncr-m+OFFSET) * sizeof(hm_t));
            const hm_t len    = ncr-m;
            const hm_t os     = len % UNROLL;
            const hm_t shift  = ncl+m;
            dss = dts + OFFSET;

            for (k = 0, j = 0; j < os; ++j) {
                if (dm[m][j] != 0) {
                    cfs[k]    = dm[m][j];
                    dss[k++]  = j+shift;
                }
            }
            for (; j < len; j += UNROLL) {
                if (dm[m][j] != 0) {
                    cfs[k]    = dm[m][j];
                    dss[k++]  = j+shift;
                }
                if (dm[m][j+1] != 0) {
                    cfs[k]    = dm[m][j+1];
                    dss[k++]  = j+1+shift;
                }
                if (dm[m][j+2] != 0) {
                    cfs[k]    = dm[m][j+2];
                    dss[k++]  = j+2+shift;
                }
                if (dm[m][j+3] != 0) {
                    cfs[k]    = dm[m][j+3];
                    dss[k++]  = j+3+shift;
                }
            }

            /* store meta data in first entries */
            dts[COEFFS]   = l; /* position of coefficient array in tmpcf */
            dts[PRELOOP]  = k % UNROLL;
            dts[LENGTH]   = k;

            /* adjust memory usage */
            dts = realloc(dts, (unsigned long)(k+OFFSET) * sizeof(hm_t));
            cfs = realloc(cfs, (unsigned long)k * sizeof(cf32_t));

            /* link to basis */
            mat->tr[l]    = dts;
            mat->cf_32[l] = cfs;
            l++;
        }
    }
}

/* NOTE: this note is about the different linear algebra implementations:
 * exact and probabilistic linear algebra differ only in the last,
 * dense reduction step: the reduction of CD via AB is sparse and
 * the same for both. this generates a dense D' part which is then
 * either reduced via exact linear algebra or via probabilistic
 * linear algebra */
static void probabilistic_sparse_linear_algebra_ff_32(
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

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32 = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    probabilistic_sparse_reduced_echelon_form_ff_32(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
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
            (unsigned long)smat->cld * sizeof(cf32_t *));

    sba_echelon_form_ff_32(smat, syz, st, ht);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += smat->nz;
}

/* In f4: tbr == bs
in nf: tbr are the polynomials to be reduced w.r.t. bs */
static void exact_sparse_linear_algebra_ff_32(
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

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    exact_sparse_reduced_echelon_form_ff_32(mat, tbr, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void copy_kernel_to_matrix(
        mat_t *mat,
        bs_t *kernel,
        const len_t nc
        )
{
    len_t i;

    sort_matrix_rows_increasing(kernel->hm, kernel->ld);

    mat->tr = (hm_t **)malloc((unsigned long)kernel->ld * sizeof(hm_t *));

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

static int exact_application_sparse_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    int ret;

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32 = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    ret = exact_application_sparse_reduced_echelon_form_ff_32(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }

    return ret;
}

static void exact_trace_sparse_linear_algebra_ff_32(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    exact_trace_sparse_reduced_echelon_form_ff_32(trace, mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void exact_sparse_dense_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_32(mat, bs, st);
    if (mat->np > 0) {
        dm  = exact_dense_linear_algebra_ff_32(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_32(dm, ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

    /* free dm */
    if (dm) {
        for (i = 0; i < ncr; ++i) {
            free(dm[i]);
        }
        free(dm);
        dm  = NULL;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void probabilistic_sparse_dense_linear_algebra_ff_32_2(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_32(mat, bs, st);
    if (mat->np > 0) {
        dm  = probabilistic_dense_linear_algebra_ff_32(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_32(dm, mat->ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

    /* free dm */
    if (dm) {
        for (i = 0; i < ncr; ++i) {
            free(dm[i]);
        }
        free(dm);
        dm  = NULL;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void probabilistic_sparse_dense_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm = NULL;
    mat->np = 0;
    dm      = probabilistic_sparse_dense_echelon_form_ff_32(mat, bs, st);
    dm      = interreduce_dense_matrix_ff_32(dm, mat->ncr, st->fc);

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

    /* free dm */
    if (dm) {
        for (i = 0; i < ncr; ++i) {
            free(dm[i]);
        }
        free(dm);
        dm  = NULL;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void interreduce_matrix_rows_ff_32(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        const int free_basis
        )
{
    len_t i, j, k, l;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;

    /* adjust displaying timings for statistic printout */
    if (st->info_level > 1) {
        printf("                          ");
    }

    mat->tr = realloc(mat->tr, (unsigned long)ncols * sizeof(hm_t *));

    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)ncols * sizeof(cf32_t *));
    memset(mat->cf_32, 0, (unsigned long)ncols * sizeof(cf32_t *));
    hm_t **pivs = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    /* copy coefficient arrays from basis in matrix, maybe
     * several rows need the same coefficient arrays, but we
     * cannot share them here. */
    for (i = 0; i < nrows; ++i) {
        pivs[mat->rr[i][OFFSET]]  = mat->rr[i];
    }

    int64_t *dr = (int64_t *)malloc((unsigned long)ncols * sizeof(int64_t));
    /* interreduce new pivots */
    cf32_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc;
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = bs->cf_32[pivs[l][COEFFS]];
            const len_t os  = pivs[l][PRELOOP];
            const len_t len = pivs[l][LENGTH];
            const len_t bi  = pivs[l][BINDEX];
            const len_t mh  = pivs[l][MULT];
            const hm_t * const ds = pivs[l] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += UNROLL) {
                dr[ds[j]]   = (int64_t)cfs[j];
                dr[ds[j+1]] = (int64_t)cfs[j+1];
                dr[ds[j+2]] = (int64_t)cfs[j+2];
                dr[ds[j+3]] = (int64_t)cfs[j+3];
            }
            free(pivs[l]);
            pivs[l] = NULL;
            pivs[l] = mat->tr[k--] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, l, mh, bi, 0,  st);
        }
    }
    if (free_basis != 0) {
        /* free now all polynomials in the basis and reset bs->ld to 0. */
        free_basis_elements(bs);
    }
    free(mat->rr);
    mat->rr = NULL;
    mat->np = nrows;
    free(pivs);
    free(dr);
}
