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

static inline cf8_t *normalize_dense_matrix_row_ff_8(
        cf8_t *row,
        const hm_t len,
        const uint32_t fc
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;
    const hm_t os     = len % UNROLL;
    const uint8_t fc8 = (uint8_t)fc;

    const uint8_t inv = mod_p_inverse_8(row[0], fc8);

    for (i = 1; i < os; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc8;
        tmp1    +=  (tmp1 >> 63) & fc8;
        row[i]  =   (cf8_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        tmp1      =   ((int64_t)row[i] * inv) % fc8;
        tmp2      =   ((int64_t)row[i+1] * inv) % fc8;
        tmp3      =   ((int64_t)row[i+2] * inv) % fc8;
        tmp4      =   ((int64_t)row[i+3] * inv) % fc8;
        tmp1      +=  (tmp1 >> 63) & fc8;
        tmp2      +=  (tmp2 >> 63) & fc8;
        tmp3      +=  (tmp3 >> 63) & fc8;
        tmp4      +=  (tmp4 >> 63) & fc8;
        row[i]    =   (cf8_t)tmp1;
        row[i+1]  =   (cf8_t)tmp2;
        row[i+2]  =   (cf8_t)tmp3;
        row[i+3]  =   (cf8_t)tmp4;
    }
    row[0]  = 1;

    return row;
}

static inline cf8_t *normalize_sparse_matrix_row_ff_8(
        cf8_t *row,
        const len_t os,
        const len_t len,
        const uint32_t fc
        )
{
    len_t i;

    const uint8_t fc8 = (uint8_t)fc;
    const uint8_t inv  = mod_p_inverse_8(row[0], fc8);

    for (i = 0; i < os; ++i) {
        row[i]    =   (cf8_t)(((uint16_t)row[i] * inv) % fc8);
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += UNROLL) {
        row[i]    =   (cf8_t)(((uint16_t)row[i] * inv) % fc8);
        row[i+1]  =   (cf8_t)(((uint16_t)row[i+1] * inv) % fc8);
        row[i+2]  =   (cf8_t)(((uint16_t)row[i+2] * inv) % fc8);
        row[i+3]  =   (cf8_t)(((uint16_t)row[i+3] * inv) % fc8);
    }
    row[0]  = 1;

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_ff_8(
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
    cf8_t *cfs;
    int64_t np = -1;
    const uint32_t fc           = st->fc;
    const int64_t mod           = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf8_t * const * const mcf  = mat->cf_8;

    rba_t *rba;
    if (tr > 0) {
        rba = mat->rba[tmp_pos];
    } else {
        rba = NULL;
    }
#if defined HAVE_AVX512_F
    __m512i mask1 = _mm512_set1_epi64(0x000000000000FFFF);
    __m512i mask2 = _mm512_set1_epi64(0x00000000FFFF0000);
    __m512i mask3 = _mm512_set1_epi64(0x0000FFFF00000000);
    __m512i mask4 = _mm512_set1_epi64(0xFFFF000000000000);
    __m512i mask8 = _mm512_set1_epi16(0x00FF);
    int64_t res[8] __attribute__((aligned(64)));
    __m512i redv, mulv, prod, drv, resv;
#elif defined HAVE_AVX2
    __m256i mask1 = _mm256_set1_epi64x(0x000000000000FFFF);
    __m256i mask2 = _mm256_set1_epi64x(0x00000000FFFF0000);
    __m256i mask3 = _mm256_set1_epi64x(0x0000FFFF00000000);
    __m256i mask4 = _mm256_set1_epi64x(0xFFFF000000000000);
    __m256i mask8 = _mm256_set1_epi16(0x00FF);
    int64_t res[4] __attribute__((aligned(32)));
    __m256i redv, mulv, prod, drv, resv;
#elif defined __aarch64__
    uint64_t tmp[2] __attribute__((aligned(32)));
    uint16x8_t prodv;
    uint32x4_t prodvl, prodvh;
    uint8x16_t redv;
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
        const uint32_t mul= (uint32_t)(fc - dr[i]);
        dts   = pivs[i];
        if (i < ncl) {
            /* set corresponding bit of reducer in reducer bit array */
            if (tr > 0) {
                rba[i/32] |= 1U << (i % 32);
            }
        }
        cfs   = mcf[dts[COEFFS]];
#if defined HAVE_AVX512_F
        const uint16_t mul8 = (uint16_t)(fc - dr[i]);
        mulv  = _mm512_set1_epi16(mul8);
        const len_t len = dts[LENGTH];
        const len_t os  = len % 64;
        const hm_t * const ds  = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 64) {
            redv  = _mm512_loadu_si512((__m512i*)(cfs+j));
            /* first 16 mult-adds -- lower */
            prod = _mm512_mullo_epi16(mulv,_mm512_and_si512(redv, mask8));
            drv   = _mm512_setr_epi64(
                dr[ds[j]],
                dr[ds[j+8]],
                dr[ds[j+16]],
                dr[ds[j+24]],
                dr[ds[j+32]],
                dr[ds[j+40]],
                dr[ds[j+48]],
                dr[ds[j+56]]);
            resv  = _mm512_add_epi64(drv, _mm512_and_si512(prod, mask1));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j]]    = res[0];
            dr[ds[j+8]]  = res[1];
            dr[ds[j+16]] = res[2];
            dr[ds[j+24]] = res[3];
            dr[ds[j+32]] = res[4];
            dr[ds[j+40]] = res[5];
            dr[ds[j+48]] = res[6];
            dr[ds[j+56]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+2]],
                dr[ds[j+10]],
                dr[ds[j+18]],
                dr[ds[j+26]],
                dr[ds[j+34]],
                dr[ds[j+42]],
                dr[ds[j+50]],
                dr[ds[j+58]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask2), 16));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+2]]  = res[0];
            dr[ds[j+10]] = res[1];
            dr[ds[j+18]] = res[2];
            dr[ds[j+26]] = res[3];
            dr[ds[j+34]] = res[4];
            dr[ds[j+42]] = res[5];
            dr[ds[j+50]] = res[6];
            dr[ds[j+58]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+4]],
                dr[ds[j+12]],
                dr[ds[j+20]],
                dr[ds[j+28]],
                dr[ds[j+36]],
                dr[ds[j+44]],
                dr[ds[j+52]],
                dr[ds[j+60]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask3), 32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+4]]  = res[0];
            dr[ds[j+12]] = res[1];
            dr[ds[j+20]] = res[2];
            dr[ds[j+28]] = res[3];
            dr[ds[j+36]] = res[4];
            dr[ds[j+44]] = res[5];
            dr[ds[j+52]] = res[6];
            dr[ds[j+60]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+6]],
                dr[ds[j+14]],
                dr[ds[j+22]],
                dr[ds[j+30]],
                dr[ds[j+38]],
                dr[ds[j+46]],
                dr[ds[j+54]],
                dr[ds[j+62]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask4), 48));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+6]]  = res[0];
            dr[ds[j+14]] = res[1];
            dr[ds[j+22]] = res[2];
            dr[ds[j+30]] = res[3];
            dr[ds[j+38]] = res[4];
            dr[ds[j+46]] = res[5];
            dr[ds[j+54]] = res[6];
            dr[ds[j+62]] = res[7];
            /* second 16 mult-adds -- higher */
            prod = _mm512_mullo_epi16(mulv, _mm512_srli_epi16(redv, 8));
            drv   = _mm512_setr_epi64(
                dr[ds[j+1]],
                dr[ds[j+9]],
                dr[ds[j+17]],
                dr[ds[j+25]],
                dr[ds[j+33]],
                dr[ds[j+41]],
                dr[ds[j+49]],
                dr[ds[j+57]]);
            resv  = _mm512_add_epi64(drv, _mm512_and_si512(prod, mask1));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+1]]  = res[0];
            dr[ds[j+9]]  = res[1];
            dr[ds[j+17]] = res[2];
            dr[ds[j+25]] = res[3];
            dr[ds[j+33]] = res[4];
            dr[ds[j+41]] = res[5];
            dr[ds[j+49]] = res[6];
            dr[ds[j+57]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+3]],
                dr[ds[j+11]],
                dr[ds[j+19]],
                dr[ds[j+27]],
                dr[ds[j+35]],
                dr[ds[j+43]],
                dr[ds[j+51]],
                dr[ds[j+59]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask2), 16));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+3]]  = res[0];
            dr[ds[j+11]] = res[1];
            dr[ds[j+19]] = res[2];
            dr[ds[j+27]] = res[3];
            dr[ds[j+35]] = res[4];
            dr[ds[j+43]] = res[5];
            dr[ds[j+51]] = res[6];
            dr[ds[j+59]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+5]],
                dr[ds[j+13]],
                dr[ds[j+21]],
                dr[ds[j+29]],
                dr[ds[j+37]],
                dr[ds[j+45]],
                dr[ds[j+53]],
                dr[ds[j+61]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask3), 32));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+5]]  = res[0];
            dr[ds[j+13]] = res[1];
            dr[ds[j+21]] = res[2];
            dr[ds[j+29]] = res[3];
            dr[ds[j+37]] = res[4];
            dr[ds[j+45]] = res[5];
            dr[ds[j+53]] = res[6];
            dr[ds[j+61]] = res[7];
            drv   = _mm512_setr_epi64(
                dr[ds[j+7]],
                dr[ds[j+15]],
                dr[ds[j+23]],
                dr[ds[j+31]],
                dr[ds[j+39]],
                dr[ds[j+47]],
                dr[ds[j+55]],
                dr[ds[j+63]]);
            resv  = _mm512_add_epi64(drv, _mm512_srli_epi64(_mm512_and_si512(prod, mask4), 48));
            _mm512_store_si512((__m512i*)(res),resv);
            dr[ds[j+7]]  = res[0];
            dr[ds[j+15]] = res[1];
            dr[ds[j+23]] = res[2];
            dr[ds[j+31]] = res[3];
            dr[ds[j+39]] = res[4];
            dr[ds[j+47]] = res[5];
            dr[ds[j+55]] = res[6];
            dr[ds[j+63]] = res[7];
        }
#elif defined HAVE_AVX2
        const uint16_t mul8 = (uint16_t)(fc - dr[i]);
        mulv  = _mm256_set1_epi16(mul8);
        const len_t len = dts[LENGTH];
        const len_t os  = len % 32;
        const hm_t * const ds  = dts + OFFSET;
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 32) {
            redv  = _mm256_loadu_si256((__m256i*)(cfs+j));
            /* first 16 mult-adds -- lower */
            prod = _mm256_mullo_epi16(mulv,_mm256_and_si256(redv, mask8));
            drv   = _mm256_setr_epi64x(
                dr[ds[j]],
                dr[ds[j+8]],
                dr[ds[j+16]],
                dr[ds[j+24]]);
            resv  = _mm256_add_epi64(drv, _mm256_and_si256(prod, mask1));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j]]    = res[0];
            dr[ds[j+8]]  = res[1];
            dr[ds[j+16]] = res[2];
            dr[ds[j+24]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+2]],
                dr[ds[j+10]],
                dr[ds[j+18]],
                dr[ds[j+26]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask2), 16));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+2]]  = res[0];
            dr[ds[j+10]] = res[1];
            dr[ds[j+18]] = res[2];
            dr[ds[j+26]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+4]],
                dr[ds[j+12]],
                dr[ds[j+20]],
                dr[ds[j+28]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask3), 32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+4]]  = res[0];
            dr[ds[j+12]] = res[1];
            dr[ds[j+20]] = res[2];
            dr[ds[j+28]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+6]],
                dr[ds[j+14]],
                dr[ds[j+22]],
                dr[ds[j+30]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask4), 48));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+6]]  = res[0];
            dr[ds[j+14]] = res[1];
            dr[ds[j+22]] = res[2];
            dr[ds[j+30]] = res[3];
            /* second 16 mult-adds -- higher */
            prod = _mm256_mullo_epi16(mulv, _mm256_srli_epi16(redv, 8));
            drv   = _mm256_setr_epi64x(
                dr[ds[j+1]],
                dr[ds[j+9]],
                dr[ds[j+17]],
                dr[ds[j+25]]);
            resv  = _mm256_add_epi64(drv, _mm256_and_si256(prod, mask1));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+1]]  = res[0];
            dr[ds[j+9]]  = res[1];
            dr[ds[j+17]] = res[2];
            dr[ds[j+25]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+3]],
                dr[ds[j+11]],
                dr[ds[j+19]],
                dr[ds[j+27]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask2), 16));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+3]]  = res[0];
            dr[ds[j+11]] = res[1];
            dr[ds[j+19]] = res[2];
            dr[ds[j+27]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+5]],
                dr[ds[j+13]],
                dr[ds[j+21]],
                dr[ds[j+29]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask3), 32));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+5]]  = res[0];
            dr[ds[j+13]] = res[1];
            dr[ds[j+21]] = res[2];
            dr[ds[j+29]] = res[3];
            drv   = _mm256_setr_epi64x(
                dr[ds[j+7]],
                dr[ds[j+15]],
                dr[ds[j+23]],
                dr[ds[j+31]]);
            resv  = _mm256_add_epi64(drv, _mm256_srli_epi64(_mm256_and_si256(prod, mask4), 48));
            _mm256_store_si256((__m256i*)(res),resv);
            dr[ds[j+7]]  = res[0];
            dr[ds[j+15]] = res[1];
            dr[ds[j+23]] = res[2];
            dr[ds[j+31]] = res[3];
        }
#elif defined __aarch64__
        const len_t len       = dts[LENGTH];
        const len_t os        = len % 16;
        const hm_t * const ds = dts + OFFSET;
        const cf8_t mul8      = (cf8_t)(mod - dr[i]);
        const uint8x8_t mulv  = vmov_n_u8(mul8);
        for (j = 0; j < os; ++j) {
            dr[ds[j]]  +=  mul * cfs[j];
        }
        for (; j < len; j += 16) {
            redv = vld1q_u8(cfs+j);

            prodv  = vmull_u8(vget_low_u8(redv), mulv);
            prodvl = vmovl_u16(vget_low_u16(prodv));
            tmp[0] = (uint64_t)dr[ds[j]];
            tmp[1] = (uint64_t)dr[ds[j+1]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_low_u32(prodvl));
            vst1q_u64(tmp, resv);
            dr[ds[j]]   = (int64_t)tmp[0];
            dr[ds[j+1]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+2]];
            tmp[1] = (uint64_t)dr[ds[j+3]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_high_u32(prodvl));
            vst1q_u64(tmp, resv);
            dr[ds[j+2]] = (int64_t)tmp[0];
            dr[ds[j+3]] = (int64_t)tmp[1];
            prodvh = vmovl_u16(vget_high_u16(prodv));
            tmp[0] = (uint64_t)dr[ds[j+4]];
            tmp[1] = (uint64_t)dr[ds[j+5]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_low_u32(prodvh));
            vst1q_u64(tmp, resv);
            dr[ds[j+4]] = (int64_t)tmp[0];
            dr[ds[j+5]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+6]];
            tmp[1] = (uint64_t)dr[ds[j+7]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_high_u32(prodvh));
            vst1q_u64(tmp, resv);
            dr[ds[j+6]] = (int64_t)tmp[0];
            dr[ds[j+7]] = (int64_t)tmp[1];

            prodv  = vmull_u8(vget_high_u8(redv), mulv);
            prodvl = vmovl_u16(vget_low_u16(prodv));
            tmp[0] = (uint64_t)dr[ds[j+8]];
            tmp[1] = (uint64_t)dr[ds[j+9]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_low_u32(prodvl));
            vst1q_u64(tmp, resv);
            dr[ds[j+8]]   = (int64_t)tmp[0];
            dr[ds[j+9]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+10]];
            tmp[1] = (uint64_t)dr[ds[j+11]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_high_u32(prodvl));
            vst1q_u64(tmp, resv);
            dr[ds[j+10]] = (int64_t)tmp[0];
            dr[ds[j+11]] = (int64_t)tmp[1];
            prodvh = vmovl_u16(vget_high_u16(prodv));
            tmp[0] = (uint64_t)dr[ds[j+12]];
            tmp[1] = (uint64_t)dr[ds[j+13]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_low_u32(prodvh));
            vst1q_u64(tmp, resv);
            dr[ds[j+12]] = (int64_t)tmp[0];
            dr[ds[j+13]] = (int64_t)tmp[1];
            tmp[0] = (uint64_t)dr[ds[j+14]];
            tmp[1] = (uint64_t)dr[ds[j+15]];
            drv    = vld1q_u64(tmp);
            resv   = vaddw_u32(drv, vget_high_u32(prodvh));
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
    cf8_t *cf  = (cf8_t *)malloc((uint64_t)(k) * sizeof(cf8_t));
    j = 0;
    hm_t *rs = row + OFFSET;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            rs[j] = (hm_t)i;
            cf[j] = (cf8_t)dr[i];
            j++;
        }
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_8[tmp_pos]  = cf;

    return row;
}

/* see la_ff_common.c */
#define LA_CF_T cf8_t
#define LA_CF cf_8
#define LA_FF(name) name##_ff_8
#define LA_RANDOM_MULTIPLIER(fc) ((int64_t)rand() % (fc))

#define LA_SMALL(name) name##_ff_8
#include "la_ff_small_prime.c"

#include "la_ff_common.c"
