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
        md_t *st
        )
{
    hi_t i, j, k;
    hm_t *dts;
    cf16_t *cfs;
    int64_t np = -1;
    const uint32_t fc           = st->fc;
    const int64_t mod           = (int64_t)fc;
    const uint64_t bm = barrett_multiplier(mod);
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
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_16[tmp_pos]  = cf;

    return row;
}

/* see la_ff_common.c */
#define LA_CF_T cf16_t
#define LA_CF cf_16
#define LA_FF(name) name##_ff_16
#define LA_RANDOM_MULTIPLIER(rng, fc) ((int64_t)msolve_rng_rand(rng) % (fc))

#define LA_SMALL(name) name##_ff_16
#include "la_ff_small_prime.c"

#include "la_ff_common.c"
