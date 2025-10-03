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

static inline cf16_t *normalize_dense_matrix_row_ff_16(
        cf16_t *row,
        const hm_t len,
        const uint32_t fc
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;
    const hm_t os       = len % UNROLL;
    const uint16_t fc16  = (uint16_t)fc;

    const uint16_t inv = mod_p_inverse_16(row[0], fc16);

    for (i = 1; i < os; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc16;
        tmp1    +=  (tmp1 >> 63) & fc16;
        row[i]  =   (cf16_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        tmp1      =   ((int64_t)row[i] * inv) % fc16;
        tmp2      =   ((int64_t)row[i+1] * inv) % fc16;
        tmp3      =   ((int64_t)row[i+2] * inv) % fc16;
        tmp4      =   ((int64_t)row[i+3] * inv) % fc16;
        tmp1      +=  (tmp1 >> 63) & fc16;
        tmp2      +=  (tmp2 >> 63) & fc16;
        tmp3      +=  (tmp3 >> 63) & fc16;
        tmp4      +=  (tmp4 >> 63) & fc16;
        row[i]    =   (cf16_t)tmp1;
        row[i+1]  =   (cf16_t)tmp2;
        row[i+2]  =   (cf16_t)tmp3;
        row[i+3]  =   (cf16_t)tmp4;
    }
    row[0]  = 1;

    return row;
}

static inline cf16_t *normalize_sparse_matrix_row_ff_16(
        cf16_t *row,
        const len_t os,
        const len_t len,
        const uint32_t fc
        )
{
    len_t i;

    const uint16_t fc16  = (uint16_t)fc;
    const uint16_t inv   = mod_p_inverse_16(row[0], fc16);

    for (i = 0; i < os; ++i) {
        row[i]  = (cf16_t)(((uint32_t)row[i] * inv) % fc16);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]   = (cf16_t)(((uint32_t)row[i] * inv) % fc16);
        row[i+1] = (cf16_t)(((uint32_t)row[i+1] * inv) % fc16);
        row[i+2] = (cf16_t)(((uint32_t)row[i+2] * inv) % fc16);
        row[i+3] = (cf16_t)(((uint32_t)row[i+3] * inv) % fc16);
    }
    row[0]  = 1;

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_ff_16(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        const uint32_t fc
        )
{
    hi_t i, j, k;
    hm_t *dts;
    cf16_t *cfs;
    int64_t np = -1;
    const int64_t mod           = (int64_t)fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf16_t * const * const mcf  = mat->cf_16;

    rba_t *rba;
    if (tr > 0) {
        rba = mat->rba[tmp_pos];
    } else {
        rba = NULL;
    }
#if defined HAVE_AVX512_F
    uint32_t mone32   = (uint32_t)0xFFFFFFFF;
    uint16_t mone16   = (uint16_t)0xFFFF;
    uint32_t mone16h  = (uint32_t)0xFFFF0000;
    __m512i mask32    = _mm512_set1_epi64(mone32);
    __m512i mask16    = _mm512_set1_epi32(mone16);
    __m512i mask16h   = _mm512_set1_epi32(mone16h);

    int64_t res[8] __attribute__((aligned(64)));
    __m512i redv, mulv, prodh, prodl, prod, drv, resv;
#elif defined HAVE_AVX2
    uint32_t mone32   = (uint32_t)0xFFFFFFFF;
    uint16_t mone16   = (uint16_t)0xFFFF;
    uint32_t mone16h  = (uint32_t)0xFFFF0000;
    __m256i mask32    = _mm256_set1_epi64x(mone32);
    __m256i mask16    = _mm256_set1_epi32(mone16);
    __m256i mask16h   = _mm256_set1_epi32(mone16h);

    int64_t res[4] __attribute__((aligned(32)));
    __m256i redv, mulv, prodh, prodl, prod, drv, resv;
#elif defined __aarch64__
    uint64_t tmp[2] __attribute__((aligned(32)));
    uint32x4_t prodv;
    uint16x8_t redv;
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
        const uint32_t mul = (uint32_t)(fc - dr[i]);
        dts   = pivs[i];
        if (i < ncl) {
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        }
        cfs   = mcf[dts[COEFFS]];
#if defined HAVE_AVX512_F
        const uint16_t mul16 = (uint16_t)(fc - dr[i]);
        mulv  = _mm512_set1_epi16(mul16);
        const len_t len = dts[LENGTH];
        const len_t os  = len % 32;
        const hm_t * const ds  = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 32) {
            redv  = _mm512_loadu_si512((__m512i*)(cfs+j));
            prodh = _mm512_mulhi_epu16(mulv, redv);
            prodl = _mm512_mullo_epi16(mulv, redv);
            prod  = _mm512_xor_si512(
                _mm512_and_si512(prodh, mask16h), _mm512_srli_epi32(prodl, 16));
            drv   = _mm512_setr_epi64(
                dr[ds[j+1]],
                dr[ds[j+5]],
                dr[ds[j+9]],
                dr[ds[j+13]],
                dr[ds[j+17]],
                dr[ds[j+21]],
                dr[ds[j+25]],
                dr[ds[j+29]]);
            resv  = _mm512_add_epi64(drv, _mm512_and_si512(prod, mask32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+1]]   = res[0];
            dr[ds[j+5]]   = res[1];
            dr[ds[j+9]]   = res[2];
            dr[ds[j+13]]  = res[3];
            dr[ds[j+17]]  = res[4];
            dr[ds[j+21]]  = res[5];
            dr[ds[j+25]]  = res[6];
            dr[ds[j+29]]  = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+3]],
                dr[ds[j+7]],
                dr[ds[j+11]],
                dr[ds[j+15]],
                dr[ds[j+19]],
                dr[ds[j+23]],
                dr[ds[j+27]],
                dr[ds[j+31]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(prod, 32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+3]]   = res[0];
            dr[ds[j+7]]   = res[1];
            dr[ds[j+11]]  = res[2];
            dr[ds[j+15]]  = res[3];
            dr[ds[j+19]]  = res[4];
            dr[ds[j+23]]  = res[5];
            dr[ds[j+27]]  = res[6];
            dr[ds[j+31]]  = res[7];
            prod  = _mm512_xor_si512(
                _mm512_slli_epi32(prodh, 16), _mm512_and_si512(prodl, mask16));
            drv   = _mm512_setr_epi64(
                dr[ds[j+0]],
                dr[ds[j+4]],
                dr[ds[j+8]],
                dr[ds[j+12]],
                dr[ds[j+16]],
                dr[ds[j+20]],
                dr[ds[j+24]],
                dr[ds[j+28]]);
            resv  = _mm512_add_epi64(drv, _mm512_and_si512(prod, mask32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+0]]   = res[0];
            dr[ds[j+4]]   = res[1];
            dr[ds[j+8]]   = res[2];
            dr[ds[j+12]]  = res[3];
            dr[ds[j+16]]  = res[4];
            dr[ds[j+20]]  = res[5];
            dr[ds[j+24]]  = res[6];
            dr[ds[j+28]]  = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+2]],
                dr[ds[j+6]],
                dr[ds[j+10]],
                dr[ds[j+14]],
                dr[ds[j+18]],
                dr[ds[j+22]],
                dr[ds[j+26]],
                dr[ds[j+30]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(prod, 32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+2]]   = res[0];
            dr[ds[j+6]]   = res[1];
            dr[ds[j+10]]  = res[2];
            dr[ds[j+14]]  = res[3];
            dr[ds[j+18]]  = res[4];
            dr[ds[j+22]]  = res[5];
            dr[ds[j+26]]  = res[6];
            dr[ds[j+30]]  = res[7];
        }
#elif defined HAVE_AVX2
        const uint16_t mul16 = (uint16_t)(fc - dr[i]);
        mulv  = _mm256_set1_epi16(mul16);
        const len_t len = dts[LENGTH];
        const len_t os  = len % 16;
        const hm_t * const ds  = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 16) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            prodh = _mm256_mulhi_epu16(mulv, redv);
            prodl = _mm256_mullo_epi16(mulv, redv);
            prod  = _mm256_xor_si256(
                _mm256_and_si256(prodh, mask16h), _mm256_srli_epi32(prodl, 16));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+5]],
                dr[ds[j+9]],
                dr[ds[j+13]]);
            resv  = _mm256_add_epi64(drv, _mm256_and_si256(prod, mask32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+1]]   = res[0];
            dr[ds[j+5]]   = res[1];
            dr[ds[j+9]]   = res[2];
            dr[ds[j+13]]  = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+3]],
                dr[ds[j+7]],
                dr[ds[j+11]],
                dr[ds[j+15]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(prod, 32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+3]]   = res[0];
            dr[ds[j+7]]   = res[1];
            dr[ds[j+11]]  = res[2];
            dr[ds[j+15]]  = res[3];
            prod  = _mm256_xor_si256(
                _mm256_slli_epi32(prodh, 16), _mm256_and_si256(prodl, mask16));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+0]],
                dr[ds[j+4]],
                dr[ds[j+8]],
                dr[ds[j+12]]);
            resv  = _mm256_add_epi64(drv, _mm256_and_si256(prod, mask32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+0]]   = res[0];
            dr[ds[j+4]]   = res[1];
            dr[ds[j+8]]   = res[2];
            dr[ds[j+12]]  = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+2]],
                dr[ds[j+6]],
                dr[ds[j+10]],
                dr[ds[j+14]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(prod, 32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+2]]   = res[0];
            dr[ds[j+6]]   = res[1];
            dr[ds[j+10]]  = res[2];
            dr[ds[j+14]]  = res[3];
        }
#elif defined __aarch64__
        const len_t len       = dts[LENGTH];
        const len_t os        = len % 8;
        const hm_t * const ds = dts + OFFSET;
        const cf16_t mul16   = (cf16_t)(mod - dr[i]);
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 8) {
            tmp[0] = (uint64_t)dr[ds[j]];
            tmp[1] = (uint64_t)dr[ds[j+1]];
            drv  = vld1q_u64(tmp);
            redv = vld1q_u16((cf16_t *)(cfs)+j);

            prodv = vmull_n_u16(vget_low_u16(redv), mul16);
            resv  = vaddw_u32(drv, vget_low_u32(prodv));
            vst1q_u64(tmp, resv);
            dr[ds[j]]   = (int64_t)tmp[0];
            dr[ds[j+1]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+2]];
            tmp[1] = (uint64_t)dr[ds[j+3]];
            drv  = vld1q_u64(tmp);
            resv  = vaddw_u32(drv, vget_high_u32(prodv));
            vst1q_u64(tmp, resv);
            dr[ds[j+2]] = (int64_t)tmp[0];
            dr[ds[j+3]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+4]];
            tmp[1] = (uint64_t)dr[ds[j+5]];
            drv  = vld1q_u64(tmp);

            prodv = vmull_n_u16( vget_high_u16(redv), mul16);
            resv  = vaddw_u32(drv, vget_low_u32(prodv));
            vst1q_u64(tmp, resv);
            dr[ds[j+4]] = (int64_t)tmp[0];
            dr[ds[j+5]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+6]];
            tmp[1] = (uint64_t)dr[ds[j+7]];
            drv  = vld1q_u64(tmp);
            resv  = vaddw_u32(drv, vget_high_u32(prodv));
            vst1q_u64(tmp, resv);
            dr[ds[j+6]] = (int64_t)tmp[0];
            dr[ds[j+7]] = (int64_t)tmp[1];
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
    cf16_t *cf  = (cf16_t *)malloc((uint64_t)(k) * sizeof(cf16_t));
    j = 0;
    hm_t *rs = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf16_t)dr[i];
            j++;
        }
    }
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_16[tmp_pos]  = cf;

    return row;
}

static hm_t *trace_reduce_dense_row_by_known_pivots_sparse_ff_16(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const uint32_t fc
        )
{
    hi_t i, j, k;
    hm_t *dts;
    cf16_t *cfs;
    int64_t np = -1;
    const int64_t mod           = (int64_t)fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf16_t * const * const mcf  = mat->cf_16;

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
        const uint32_t mul = (uint32_t)(fc - dr[i]);
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_16[dts[COEFFS]];
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
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf16_t *cf  = (cf16_t *)malloc((uint64_t)(k) * sizeof(cf16_t));
    j = 0;
    hm_t *rs = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf16_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_16[tmp_pos]  = cf;

    return row;
}

static cf16_t *reduce_dense_row_by_all_pivots_ff_16(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const * const pivs,
        cf16_t *const * const dpivs,
        const uint32_t fc
        )
{
    hi_t i, j, k, l;
    const int64_t mod = (int64_t)fc;
    len_t np  = -1;
    cf16_t *red;

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
        const uint32_t mul    = (uint32_t)(fc - dr[i]);
        const cf16_t *cfs     = bs->cf_16[pivs[i][COEFFS]];
        const len_t os        = pivs[i][PRELOOP];
        const len_t len       = pivs[i][LENGTH];
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
        const uint32_t mul  = (uint32_t)(fc - dr[i]);
        const len_t os      = (ncols - i) % UNROLL;
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

    cf16_t *row = (cf16_t *)calloc((uint64_t)(ncols-np), sizeof(cf16_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf16_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_16(row, ncols-np, fc);
    }

    *pc = np - ncl;
    return row;
}

static cf16_t *reduce_dense_row_by_old_pivots_ff_16(
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
        const uint32_t mul    = (uint32_t)(fc - dr[i]);
        const cf16_t *cfs     = bs->cf_16[pivs[i][COEFFS]];
        const len_t os        = pivs[i][PRELOOP];
        const len_t len       = pivs[i][LENGTH];
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
    cf16_t *row  = (cf16_t *)calloc(
            (uint64_t)(ncr), sizeof(cf16_t));
    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf16_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }
    return row;
}

static cf16_t *reduce_dense_row_by_dense_new_pivots_ff_16(
        int64_t *dr,
        len_t *pc,
        cf16_t * const * const pivs,
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

        const uint32_t mul  = (uint32_t)(fc - dr[i]);
        const len_t os      = (ncr - i) % UNROLL;
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

    cf16_t *row = (cf16_t *)calloc((uint64_t)(ncr-np), sizeof(cf16_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf16_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_16(row, ncr-np, fc);
    }
    *pc = np;
    return row;
}

static void probabilistic_sparse_reduced_echelon_form_ff_16(
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
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_16[j]      = bs->cf_16[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    const uint32_t fc   = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (uint64_t)ncols * st->nthrds * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (uint64_t)rpb * st->nthrds * sizeof(int64_t));

    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) \
    schedule(dynamic)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            hm_t *npiv  = NULL;
            cf16_t *cfs;
            /* starting column, offset, coefficient array position in tmpcf */
            hm_t sc, cfp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                cfp = bctr + i*rpb;
                sc  = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    cfs   = bs->cf_16[npiv[COEFFS]];
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
                    npiv  = reduce_dense_row_by_known_pivots_sparse_ff_16(
                            drl, mat, bs, pivs, sc, cfp, 0, 0, 0, st->fc);
                    if (!npiv) {
                        bctr  = nrbl;
                        break;
                    }
                    /* normalize coefficient array
                    * NOTE: this has to be done here, otherwise the reduction may
                    * lead to wrong results in a parallel computation since other
                    * threads might directly use the new pivot once it is synced. */
                    if (mat->cf_16[npiv[COEFFS]][0] != 1) {
                        normalize_sparse_matrix_row_ff_16(
                                mat->cf_16[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    cfs = mat->cf_16[npiv[COEFFS]];
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

    /* construct the trace */
    if (st->trace_level == LEARN_TRACER && st->in_final_reduction_step == 0) {
        construct_trace(st->tr, mat);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf16_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc, cfp;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = mat->cf_16[pivs[k][COEFFS]];
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
                reduce_dense_row_by_known_pivots_sparse_ff_16(
                        dr, mat, bs, pivs, sc, cfp, mh, bi, 0, st->fc);
        }
    }
    free(pivs);
    pivs  = NULL;

    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}

static void exact_trace_sparse_reduced_echelon_form_ff_16(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc    = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    const int32_t nthrds = st->in_final_reduction_step == 1 ? 1 : st->nthrds;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_16[j]      = bs->cf_16[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }


    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)ncols * nthrds * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        int64_t *drl    = dr + (omp_get_thread_num() * (uint64_t)ncols);
        hm_t *npiv      = upivs[i];
        rba_t *rba      = mat->rba[i];
        cf16_t *cfs     = bs->cf_16[npiv[COEFFS]];
        const len_t bi  = npiv[BINDEX];
        const len_t mh  = npiv[MULT];
        const len_t os  = npiv[PRELOOP];
        const len_t len = npiv[LENGTH];
        const hm_t * const ds = npiv + OFFSET;
        k = 0;
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
        cfs = NULL;
        do {
            sc  = npiv[OFFSET];
            free(npiv);
            npiv  = NULL;
            free(cfs);
            npiv  = mat->tr[i]  = trace_reduce_dense_row_by_known_pivots_sparse_ff_16(
                    rba, drl, mat, bs, pivs, sc, i, mh, bi, st->fc);
            if (!npiv) {
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mat->cf_16[npiv[COEFFS]][0] != 1) {
                normalize_sparse_matrix_row_ff_16(
                        mat->cf_16[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
            cfs = mat->cf_16[npiv[COEFFS]];
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

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf16_t *cfs;
    hm_t cf_array_pos;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = mat->cf_16[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
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
                dr[ds[j]]    = (int64_t)cfs[j];
                dr[ds[j+1]]  = (int64_t)cfs[j+1];
                dr[ds[j+2]]  = (int64_t)cfs[j+2];
                dr[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_16(
                        dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st->fc);
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}


static void exact_sparse_reduced_echelon_form_ff_16(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc    = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    const int32_t nthrds = st->in_final_reduction_step == 1 ? 1 : st->nthrds;

    len_t bad_prime = 0;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    if (st->in_final_reduction_step == 0) {
        memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    } else {
        for (i = 0;  i < mat->nru; ++i) {
            pivs[mat->rr[i][OFFSET]] = mat->rr[i];
        }
    }
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_16[j]      = bs->cf_16[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }


    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)ncols * nthrds * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        if (bad_prime == 0) {
            int64_t *drl    = dr + (omp_get_thread_num() * (uint64_t)ncols);
            hm_t *npiv      = upivs[i];
            cf16_t *cfs     = tbr->cf_16[npiv[COEFFS]];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const hm_t * const ds = npiv + OFFSET;
            k = 0;
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
            cfs = NULL;
            do {
                /* If we do normal form computations the first monomial in the polynomial might not
                be a known pivot, thus setting it to npiv[OFFSET] can lead to wrong results. */
                sc  = st->nf == 0 ? npiv[OFFSET] : 0;
                free(npiv);
                npiv  = NULL;
                free(cfs);
                cfs = NULL;
                npiv  = mat->tr[i] = reduce_dense_row_by_known_pivots_sparse_ff_16(
                        drl, mat, bs, pivs, sc, i, mh, bi, st->trace_level == LEARN_TRACER, st->fc);
                if (st->nf > 0) {
                    if (!npiv) {
                        mat->tr[i]  = NULL;
                        break;
                    }
                    mat->tr[i]  = npiv;
                    cfs = mat->cf_16[npiv[COEFFS]];
                    break;
                } else {
                    if (!npiv) {
                        break;
                    }
                    /* normalize coefficient array
                     * NOTE: this has to be done here, otherwise the reduction may
                     * lead to wrong results in a parallel computation since other
                     * threads might directly use the new pivot once it is synced. */
                    if (mat->cf_16[npiv[COEFFS]][0] != 1) {
                        normalize_sparse_matrix_row_ff_16(
                                mat->cf_16[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                    cfs = mat->cf_16[npiv[COEFFS]];
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
    if (st->trace_level == LEARN_TRACER && st->in_final_reduction_step == 0) {
        construct_trace(st->tr, mat);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    if (st->nf == 0 && st->in_final_reduction_step == 0) {
        dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
        mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

        /* interreduce new pivots */
        cf16_t *cfs;
        hm_t cf_array_pos;
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
                cfs = mat->cf_16[pivs[k][COEFFS]];
                cf_array_pos    = pivs[k][COEFFS];
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
                    dr[ds[j]]    = (int64_t)cfs[j];
                    dr[ds[j+1]]  = (int64_t)cfs[j+1];
                    dr[ds[j+2]]  = (int64_t)cfs[j+2];
                    dr[ds[j+3]]  = (int64_t)cfs[j+3];
                }
                free(pivs[k]);
                free(cfs);
                pivs[k] = NULL;
                pivs[k] = mat->tr[npivs++] =
                    reduce_dense_row_by_known_pivots_sparse_ff_16(
                            dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st->fc);
            }
        }
        mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
        st->np = mat->np = mat->nr = mat->sz = npivs;
    } else {
        st->np = mat->np = mat->nr = mat->sz = nrl;
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;
}

static int exact_application_sparse_reduced_echelon_form_ff_16(
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

    const int32_t nthrds = st->in_final_reduction_step == 1 ? 1 : st->nthrds;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_16[j]      = bs->cf_16[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }


    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)ncols * nthrds * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
    int flag  = 1;
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        if (flag == 1) {
            int64_t *drl    = dr + (omp_get_thread_num() * (uint64_t)ncols);
            hm_t *npiv      = upivs[i];
            cf16_t *cfs     = bs->cf_16[npiv[COEFFS]];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const hm_t * const ds = npiv + OFFSET;
            k = 0;
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
            cfs = NULL;
            do {
                sc  = npiv[OFFSET];
                free(npiv);
                free(cfs);
                npiv  = mat->tr[i]  = reduce_dense_row_by_known_pivots_sparse_ff_16(
                        drl, mat, bs, pivs, sc, i, mh, bi, 0, st->fc);
                if (!npiv) {
                    fprintf(stderr, "Unlucky prime detected, row reduced to zero.");
                    flag  = 0;
                }

                /* normalize coefficient array
                 * NOTE: this has to be done here, otherwise the reduction may
                 * lead to wrong results in a parallel computation since other
                 * threads might directly use the new pivot once it is synced. */
                if (mat->cf_16[npiv[COEFFS]][0] != 1) {
                    normalize_sparse_matrix_row_ff_16(
                            mat->cf_16[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                }
                k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                cfs = mat->cf_16[npiv[COEFFS]];
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

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf16_t *cfs;
    hm_t cf_array_pos;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = mat->cf_16[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
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
                dr[ds[j]]    = (int64_t)cfs[j];
                dr[ds[j+1]]  = (int64_t)cfs[j+1];
                dr[ds[j+2]]  = (int64_t)cfs[j+2];
                dr[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_16(
                        dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st->fc);
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;

    return 0;
}

static cf16_t **sparse_AB_CD_linear_algebra_ff_16(
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
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* dense rows representing updated D part;
     * after reducing CD part with AB */
    cf16_t **drs    = (cf16_t **)calloc((uint64_t)nrl, sizeof(cf16_t *));

    int64_t *dr  = (int64_t *)malloc(
            (uint64_t)ncols * st->nthrds * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        cf16_t *cfs = NULL;
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols);
        hm_t *npiv  = upivs[i];
        /* do the reduction */
        memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));
        cfs = bs->cf_16[npiv[COEFFS]];
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
        sc  = 0; /* row columns are not sorted, so we have to start at 0 */
        free(npiv);
        drs[i]  = reduce_dense_row_by_old_pivots_ff_16(
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

static cf16_t **interreduce_dense_matrix_ff_16(
    cf16_t **dm,
    const len_t ncr,
    const uint32_t fc
    )
{
    len_t i, j, k, l;
    int64_t *dr = malloc((uint64_t)ncr * sizeof(int64_t));

    for (i = 0; i < ncr; ++i) {
        k = ncr-1-i;
        if (dm[k]) {
            /* printf("interreduce dm[%d]\n", i); */
            memset(dr, 0, (uint64_t)ncr * sizeof(int64_t));
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
            dm[k] = reduce_dense_row_by_dense_new_pivots_ff_16(
                    dr, &k, dm, ncr, fc);
        }
    }
    free(dr);
    return dm;
}

static cf16_t **exact_dense_linear_algebra_ff_16(
        cf16_t **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, l, npivs;

    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    cf16_t **nps  = (cf16_t **)calloc((uint64_t)ncr, sizeof(cf16_t *));
    /* rows to be further reduced */
    cf16_t **tbr  = (cf16_t **)calloc((uint64_t)nrows, sizeof(cf16_t *));
    int64_t *dr   = (int64_t *)malloc(
            (uint64_t)ncr * st->nthrds * sizeof(int64_t));

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
                /* printf("ncr - k = %d - %d\n", ncr, k); */
                memmove(dm[i], dm[i]+k, (uint64_t)(ncr-k) * sizeof(cf16_t));
                dm[i] = realloc(dm[i], (uint64_t)(ncr-k) * sizeof(cf16_t));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_16(nps[k], ncr-k, st->fc);
                }
            } else {
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (uint64_t)ntr * sizeof(cf16_t *));
    /* offset modulo 4 for loop unrolling, +1 due to storing the first
     * nonzero entry at the first position */

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(dynamic)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncr);
        memset(drl, 0, (uint64_t)ncr * sizeof(int64_t));
        hm_t npc  = 0;
        hm_t os   = 0;
        cf16_t *npiv  = tbr[i];
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
            npiv = reduce_dense_row_by_dense_new_pivots_ff_16(
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

static cf16_t **probabilistic_dense_linear_algebra_ff_16(
        cf16_t **dm,
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
    cf16_t **nps  = (cf16_t **)calloc((uint64_t)ncr, sizeof(cf16_t *));
    /* rows to be further reduced */
    cf16_t **tbr  = (cf16_t **)calloc((uint64_t)nrows, sizeof(cf16_t *));

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
                memmove(dm[i], dm[i]+k, (uint64_t)(ncr-k) * sizeof(cf16_t));
                dm[i] = realloc(dm[i], (uint64_t)(ncr-k) * sizeof(cf16_t));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_16(nps[k], ncr-k, st->fc);
                }
            } else {
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (uint64_t)ntr * sizeof(cf16_t *));

    const int64_t mod2  = (int64_t)st->fc * st->fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(ntr/3)))+1;
    const len_t rem = (ntr % nb == 0) ? 0 : 1;
    const len_t rpb = (ntr / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (uint64_t)ncols * st->nthrds * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (uint64_t)rpb * st->nthrds * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(dynamic)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncr);
        int64_t *mull = mul + (omp_get_thread_num() * (uint64_t)rpb);
        const int32_t nbl   = (int32_t) (ntr > (i+1)*rpb ? (i+1)*rpb : ntr);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);

        if (nrbl > 0) {
            hm_t npc;
            hm_t os;
            cf16_t *tmp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;
                os  = ncr % UNROLL;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (uint64_t)ncr * sizeof(int64_t));

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
                    tmp = reduce_dense_row_by_dense_new_pivots_ff_16(
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

static cf16_t **probabilistic_sparse_dense_echelon_form_ff_16(
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
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* rows already representing new pivots */
    cf16_t **nps  = (cf16_t **)calloc((uint64_t)ncr, sizeof(cf16_t *));

    const uint32_t fc   = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (uint64_t)ncols * st->nthrds * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (uint64_t)rpb * st->nthrds * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) shared(nps) \
    schedule(dynamic)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols);
        int64_t *mull = mul + (omp_get_thread_num() * (uint64_t)rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl > 0) {
            hm_t *npiv;
            cf16_t *tmp;
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
                memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    tmp   = bs->cf_16[npiv[COEFFS]];
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
                    tmp = reduce_dense_row_by_all_pivots_ff_16(
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

static void convert_to_sparse_matrix_rows_ff_16(
        mat_t *mat,
        cf16_t * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, j, k, l, m;
    cf16_t *cfs;
    hm_t *dts, *dss;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->tr     = realloc(mat->tr, (uint64_t)mat->np * sizeof(hm_t *));
    mat->cf_16  = realloc(mat->cf_16,
            (uint64_t)mat->np * sizeof(cf16_t *));

    l = 0;
    for (i = 0; i < ncr; ++i) {
        m = ncr-1-i;
        if (dm[m] != NULL) {
            cfs = malloc((uint64_t)(ncr-m) * sizeof(cf16_t));
            dts = malloc((uint64_t)(ncr-m+OFFSET) * sizeof(hm_t));
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
            dts = realloc(dts, (uint64_t)(k+OFFSET) * sizeof(hm_t));
            cfs = realloc(cfs, (uint64_t)k * sizeof(cf16_t));

            /* link to basis */
            mat->tr[l]    = dts;
            mat->cf_16[l] = cfs;
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
static void probabilistic_sparse_linear_algebra_ff_16(
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
     * coefficients of all pivot rows */
    mat->cf_16 = realloc(mat->cf_16,
            (uint64_t)mat->nr * sizeof(cf16_t *));
    probabilistic_sparse_reduced_echelon_form_ff_16(mat, bs, st);

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

/* In f4: tbr == bs
in nf: tbr are the polynomials to be reduced w.r.t. bs */
static void exact_sparse_linear_algebra_ff_16(
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
     * coefficients of all pivot rows */
    mat->cf_16  = realloc(mat->cf_16,
            (uint64_t)mat->nr * sizeof(cf16_t *));
    exact_sparse_reduced_echelon_form_ff_16(mat, tbr, bs, st);

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

static int exact_application_sparse_linear_algebra_ff_16(
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
     * coefficients of all pivot rows */
    mat->cf_16  = realloc(mat->cf_16,
            (uint64_t)mat->nr * sizeof(cf16_t *));
    ret = exact_application_sparse_reduced_echelon_form_ff_16(mat, bs, st);

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

static void exact_trace_sparse_linear_algebra_ff_16(
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
     * coefficients of all pivot rows */
    mat->cf_16  = realloc(mat->cf_16,
            (uint64_t)mat->nr * sizeof(cf16_t *));
    exact_trace_sparse_reduced_echelon_form_ff_16(trace, mat, bs, st);

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

static void exact_sparse_dense_linear_algebra_ff_16(
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
    cf16_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_16(mat, bs, st);
    if (mat->np > 0) {      
        dm  = exact_dense_linear_algebra_ff_16(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_16(dm, ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_16(mat, dm);

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

static void probabilistic_sparse_dense_linear_algebra_ff_16_2(
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
    cf16_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_16(mat, bs, st);
    if (mat->np > 0) {      
        dm  = probabilistic_dense_linear_algebra_ff_16(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_16(dm, mat->ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_16(mat, dm);

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

static void probabilistic_sparse_dense_linear_algebra_ff_16(
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
    cf16_t **dm = NULL;
    mat->np = 0;
    dm      = probabilistic_sparse_dense_echelon_form_ff_16(mat, bs, st);
    dm      = interreduce_dense_matrix_ff_16(dm, mat->ncr, st->fc);

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_16(mat, dm);

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

static void interreduce_matrix_rows_ff_16(
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

    /* for interreduction steps like the final basis reduction we
    need to allocate memory for rba here, even so we do not use
    it at all */
    mat->rba  = (rba_t **)malloc((uint64_t)ncols * sizeof(rba_t *));
    const uint64_t len = ncols / 32 + ((ncols % 32) != 0);
    for (i = 0; i < ncols; ++i) {
        mat->rba[i] = (rba_t *)calloc(len, sizeof(rba_t));
    }

    mat->tr = realloc(mat->tr, (uint64_t)ncols * sizeof(hm_t *));

    mat->cf_16  = realloc(mat->cf_16,
            (uint64_t)ncols * sizeof(cf16_t *));
    memset(mat->cf_16, 0, (uint64_t)ncols * sizeof(cf16_t *));
    hm_t **pivs = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    /* copy coefficient arrays from basis in matrix, maybe
     * several rows need the same coefficient arrays, but we
     * cannot share them here. */
    for (i = 0; i < nrows; ++i) {
        pivs[mat->rr[i][OFFSET]]  = mat->rr[i];
    }

    int64_t *dr = (int64_t *)malloc((uint64_t)ncols * sizeof(int64_t));
    /* interreduce new pivots */
    cf16_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc;
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = bs->cf_16[pivs[l][COEFFS]];
            const len_t bi  = pivs[l][BINDEX];
            const len_t mh  = pivs[l][MULT];
            const len_t os  = pivs[l][PRELOOP];
            const len_t len = pivs[l][LENGTH];
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
                reduce_dense_row_by_known_pivots_sparse_ff_16(
                        dr, mat, bs, pivs, sc, l, mh, bi, 0, st->fc);
        }
    }
    for (i = 0; i < ncols; ++i) {
        free(mat->rba[i]);
        mat->rba[i] = NULL;
    }
    if (free_basis != 0) {
        /* free now all polynomials in the basis and reset bs->ld to 0. */
        free_basis_elements(bs);
    }
    free(mat->rr);
    mat->rr = NULL;
    st->np = mat->np = nrows;
    free(pivs);
    free(dr);
}
