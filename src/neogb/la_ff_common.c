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


/* Field independent part of the finite field linear algebra.
 *
 * This file is a template: it is not compiled on its own, but included
 * once at the end of each of la_ff_8.c, la_ff_16.c, la_ff_32.c and
 * la_ff_64.c. Those files contain the reduction kernels specific to the
 * corresponding coefficient size and define the following macros before
 * including this file:
 *
 * LA_CF_T                  coefficient type, e.g. cf8_t
 * LA_CF                    name of the coefficient arrays in bs_t and
 *                          mat_t, e.g. cf_8
 * LA_FF(name)              appends the field suffix to name, e.g.
 *                          LA_FF(foo) expands to foo_ff_8
 * LA_RANDOM_MULTIPLIER(fc) random multiplier for the probabilistic
 *                          linear algebra
 *
 * Moreover, the following functions have to be defined:
 *
 * LA_FF(normalize_dense_matrix_row)
 * LA_FF(normalize_sparse_matrix_row)
 * LA_FF(reduce_dense_row_by_known_pivots_sparse)
 * LA_FF(trace_reduce_dense_row_by_known_pivots_sparse)
 * LA_FF(reduce_dense_row_by_all_pivots)
 * LA_FF(reduce_dense_row_by_dense_new_pivots)
 *
 * If LA_WIDE is defined (la_ff_64.c), the products of coefficients do
 * not fit into 64 bits: the dense rows are always reduced to [0,p), the
 * helpers ff64_init(), ff64_mont(), ff64_submul(), ff64_submul_dense()
 * and ff64_submul_sparse() are used for the updates, and the functions
 * LA_FF(reduce_dense_rows_by_known_pivots) and
 * LA_FF(reduce_dense_rows_by_dense_block) have to be defined as well.
 *
 * All macros are undefined at the end of this file. */

#if defined HAVE_AVX2 || defined HAVE_AVX512_F
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

/* maximal number of rows reduced together by the known pivots, see
 * LA_FF(reduce_dense_rows_by_known_pivots) */
#define LA_NR 8

/* Number of rows reduced together by the known pivots when reducing nrl
 * rows with nthrds threads. The rows of a group are interleaved with
 * LA_NR rows, so a group with fewer rows costs as much as a full one.
 * Thus we only use full groups and only if each thread gets at least two
 * of them, otherwise matrices with few rows to be reduced would only use
 * some of the threads. Returns 0 if the rows should be reduced one by
 * one. */
static inline len_t LA_FF(rows_per_group)(
        const len_t nrl,
        const int32_t nthrds
        )
{
    return (uint64_t)nrl >= 2 * (uint64_t)nthrds * LA_NR ? LA_NR : 0;
}

#ifndef LA_WIDE
/* Reduces LA_NR rows at once by the known pivots pivs[sc],...,
 * pivs[ncl-1]. The dense rows are interleaved in dr, i.e. dr[LA_NR*c+r]
 * is the entry of row r in column c, so applying a pivot updates all
 * rows with one vector operation per term of the pivot. The remainders
 * of the rows live in the columns ncl,...,ncols-1, they are written to
 * drs[0],...,drs[LA_NR-1] as dense rows of length ncr, NULL for zero
 * remainders. If drs is NULL, the remainders are left in dr. If rbas is
 * not NULL, the pivots used for row r are marked in the reducer bit
 * array rbas[r] for the tracer, rbas[r] may be NULL if row r is zero. */
static void LA_FF(reduce_dense_rows_by_known_pivots)(
        int64_t *dr,
        LA_CF_T **drs,
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
    const int64_t mod2  = (int64_t)fc * fc;
    const uint64_t bm   = barrett_multiplier(fc);
    LA_CF_T * const * const mcf = mat->LA_CF;
#if defined HAVE_AVX512_F
    const __m512i zerov = _mm512_setzero_si512();
    const __m512i mod2v = _mm512_set1_epi64(mod2);
#elif defined HAVE_AVX2
    const __m256i zerov = _mm256_setzero_si256();
    const __m256i mod2v = _mm256_set1_epi64x(mod2);
#elif defined __aarch64__
    const int64x2_t mod2v = vdupq_n_s64(mod2);
#endif

    for (i = sc; i < ncl; ++i) {
        int64_t *d = dr + (uint64_t)LA_NR * i;
        int64_t nz = 0;
        for (r = 0; r < LA_NR; ++r) {
            if (d[r] != 0) {
                d[r] = (int64_t)barrett_reduce((uint64_t)d[r], fc, bm);
            }
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
        const LA_CF_T * const cfs = mcf[piv[COEFFS]];
        const len_t len         = piv[LENGTH];
        const hm_t * const ds   = piv + OFFSET;
#if defined HAVE_AVX512_F
        __m512i mulv[LA_NR/8];
        for (r = 0; r < LA_NR/8; ++r) {
            mulv[r] = _mm512_loadu_si512((void *)(d+8*r));
        }
        for (j = 0; j < len; ++j) {
            int64_t *t        = dr + (uint64_t)LA_NR * ds[j];
            const __m512i cfv = _mm512_set1_epi64((int64_t)cfs[j]);
            for (r = 0; r < LA_NR/8; ++r) {
                __m512i prodv = _mm512_mul_epu32(mulv[r], cfv);
                __m512i resv  = _mm512_sub_epi64(
                        _mm512_loadu_si512((void *)(t+8*r)), prodv);
                const __mmask8 cmpv = _mm512_cmpgt_epi64_mask(zerov, resv);
                resv  = _mm512_mask_add_epi64(resv, cmpv, resv, mod2v);
                _mm512_storeu_si512((void *)(t+8*r), resv);
            }
        }
#elif defined HAVE_AVX2
        __m256i mulv[LA_NR/4];
        for (r = 0; r < LA_NR/4; ++r) {
            mulv[r] = _mm256_loadu_si256((__m256i *)(d+4*r));
        }
        for (j = 0; j < len; ++j) {
            int64_t *t        = dr + (uint64_t)LA_NR * ds[j];
            const __m256i cfv = _mm256_set1_epi64x((int64_t)cfs[j]);
            for (r = 0; r < LA_NR/4; ++r) {
                __m256i prodv = _mm256_mul_epu32(mulv[r], cfv);
                __m256i resv  = _mm256_sub_epi64(
                        _mm256_loadu_si256((__m256i *)(t+4*r)), prodv);
                resv  = _mm256_add_epi64(resv, _mm256_and_si256(
                            _mm256_cmpgt_epi64(zerov, resv), mod2v));
                _mm256_storeu_si256((__m256i *)(t+4*r), resv);
            }
        }
#elif defined __aarch64__
        uint32x2_t mulv[LA_NR/2];
        for (r = 0; r < LA_NR/2; ++r) {
            mulv[r] = vmovn_u64(vld1q_u64((const uint64_t *)(d+2*r)));
        }
        for (j = 0; j < len; ++j) {
            int64_t *t            = dr + (uint64_t)LA_NR * ds[j];
            const uint32x2_t cfv  = vdup_n_u32((uint32_t)cfs[j]);
            for (r = 0; r < LA_NR/2; ++r) {
                int64x2_t resv  = vsubq_s64(vld1q_s64(t+2*r),
                        vreinterpretq_s64_u64(vmull_u32(mulv[r], cfv)));
                const int64x2_t mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
                resv  = vaddq_s64(resv, vandq_s64(mask, mod2v));
                vst1q_s64(t+2*r, resv);
            }
        }
#else
        int64_t mul[LA_NR];
        for (r = 0; r < LA_NR; ++r) {
            mul[r] = d[r];
        }
        for (j = 0; j < len; ++j) {
            int64_t *t  = dr + (uint64_t)LA_NR * ds[j];
            const int64_t c = (int64_t)cfs[j];
            for (r = 0; r < LA_NR; ++r) {
                t[r] -= mul[r] * c;
                t[r] += (t[r] >> 63) & mod2;
            }
        }
#endif
        for (r = 0; r < LA_NR; ++r) {
            d[r] = 0;
        }
    }
    if (drs == NULL) {
        return;
    }
    const len_t ncr = ncols - ncl;
    for (r = 0; r < LA_NR; ++r) {
        LA_CF_T *row = NULL;
        for (i = ncl; i < ncols; ++i) {
            int64_t *d = dr + (uint64_t)LA_NR * i + r;
            if (*d == 0) {
                continue;
            }
            const LA_CF_T c = (LA_CF_T)barrett_reduce((uint64_t)*d, fc, bm);
            *d = 0;
            if (c != 0) {
                if (row == NULL) {
                    row = (LA_CF_T *)calloc((uint64_t)ncr, sizeof(LA_CF_T));
                }
                row[i-ncl] = c;
            }
        }
        drs[r] = row;
    }
}

/* load 4 resp. 8 coefficients as 64 bit integers */
#if defined HAVE_AVX512_F
static inline __m512i LA_FF(load_cfs_8)(
        const LA_CF_T * const p
        )
{
    if (sizeof(LA_CF_T) == 4) {
        return _mm512_cvtepu32_epi64(_mm256_loadu_si256((__m256i *)p));
    }
    if (sizeof(LA_CF_T) == 2) {
        return _mm512_cvtepu16_epi64(_mm_loadu_si128((__m128i *)p));
    }
    return _mm512_cvtepu8_epi64(_mm_loadl_epi64((__m128i *)p));
}
#elif defined HAVE_AVX2
static inline __m256i LA_FF(load_cfs_4)(
        const LA_CF_T * const p
        )
{
    if (sizeof(LA_CF_T) == 4) {
        return _mm256_cvtepu32_epi64(_mm_loadu_si128((__m128i *)p));
    }
    if (sizeof(LA_CF_T) == 2) {
        return _mm256_cvtepu16_epi64(_mm_loadl_epi64((__m128i *)p));
    }
    int32_t v;
    memcpy(&v, p, sizeof(int32_t));
    return _mm256_cvtepu8_epi64(_mm_cvtsi32_si128(v));
}
#endif

/* number of rows reduced together in LA_FF(interreduce_dense_matrix) */
#define LA_NRB 2

/* Reduces the nr <= LA_NRB dense rows rows[0],...,rows[nr-1], whose pivots
 * are in the columns rc[0],...,rc[nr-1], by the nb pivots
 * dm[bpc[0]],...,dm[bpc[nb-1]] with rc[r] < bpc[0] < ... < bpc[nb-1].
 * The rows of dm are dense and start at their pivot column, i.e. rows[r]
 * has ncr-rc[r] entries. The pivots have to be zero in the pivot columns
 * of each other, so the coefficients of the rows in all pivot columns
 * can be used as multipliers without updating them in between. Each
 * loaded part of a pivot is applied to all rows. dr is a buffer for
 * LA_NRB*ncr entries. */
static void LA_FF(reduce_dense_rows_by_dense_block)(
        LA_CF_T **rows,
        const len_t * const rc,
        len_t nr,
        const len_t * const bpc,
        const len_t nb,
        LA_CF_T * const * const dm,
        const len_t ncr,
        const uint64_t fc,
        int64_t *dr
        )
{
    len_t i, j, r;
    LA_CF_T *rs[LA_NRB];
    int64_t *drs[LA_NRB];

    /* skip the rows which are zero in all pivot columns of the block */
    len_t n = 0;
    for (r = 0; r < nr; ++r) {
        LA_CF_T * const row = rows[r] - rc[r];
        for (i = 0; i < nb; ++i) {
            if (row[bpc[i]] != 0) {
                break;
            }
        }
        if (i < nb) {
            rs[n]   = row;
            drs[n]  = dr + (uint64_t)n * ncr;
            n++;
        }
    }
    if (n == 0) {
        return;
    }
    nr = n;
    const len_t fb      = bpc[0];
    const int64_t mod2  = (int64_t)fc * fc;
    const uint64_t bm   = barrett_multiplier(fc);

    for (r = 0; r < nr; ++r) {
        for (j = fb; j < ncr; ++j) {
            drs[r][j] = (int64_t)rs[r][j];
        }
    }
    /* the loops below always handle LA_NRB rows, missing rows get the
     * multiplier zero and a scratch row in dr */
    for (r = nr; r < LA_NRB; ++r) {
        drs[r]  = dr + (uint64_t)r * ncr;
    }
#if defined HAVE_AVX512_F
    const __m512i zerov = _mm512_setzero_si512();
    const __m512i mod2v = _mm512_set1_epi64(mod2);
#elif defined HAVE_AVX2
    const __m256i zerov = _mm256_setzero_si256();
    const __m256i mod2v = _mm256_set1_epi64x(mod2);
#elif defined __aarch64__
    const int64x2_t mod2v = vdupq_n_s64(mod2);
#endif
    for (i = 0; i < nb; ++i) {
        const len_t c = bpc[i];
        int64_t mul[LA_NRB] = {0};
        int64_t nz = 0;
        for (r = 0; r < nr; ++r) {
            mul[r]  = (int64_t)rs[r][c];
            nz      |= mul[r];
        }
        if (nz == 0) {
            continue;
        }
        const LA_CF_T * const piv = dm[c] - c;
        j = c;
#if defined HAVE_AVX512_F
        __m512i mulv[LA_NRB];
        for (r = 0; r < LA_NRB; ++r) {
            mulv[r] = _mm512_set1_epi64(mul[r]);
        }
        for (; j + 8 <= ncr; j += 8) {
            const __m512i cfv = LA_FF(load_cfs_8)(piv+j);
            for (r = 0; r < LA_NRB; ++r) {
                __m512i resv  = _mm512_sub_epi64(
                        _mm512_loadu_si512((void *)(drs[r]+j)),
                        _mm512_mul_epu32(mulv[r], cfv));
                const __mmask8 cmpv = _mm512_cmpgt_epi64_mask(zerov, resv);
                resv  = _mm512_mask_add_epi64(resv, cmpv, resv, mod2v);
                _mm512_storeu_si512((void *)(drs[r]+j), resv);
            }
        }
#elif defined HAVE_AVX2
        __m256i mulv[LA_NRB];
        for (r = 0; r < LA_NRB; ++r) {
            mulv[r] = _mm256_set1_epi64x(mul[r]);
        }
        for (; j + 4 <= ncr; j += 4) {
            const __m256i cfv = LA_FF(load_cfs_4)(piv+j);
            for (r = 0; r < LA_NRB; ++r) {
                __m256i resv  = _mm256_sub_epi64(
                        _mm256_loadu_si256((__m256i *)(drs[r]+j)),
                        _mm256_mul_epu32(mulv[r], cfv));
                resv  = _mm256_add_epi64(resv, _mm256_and_si256(
                            _mm256_cmpgt_epi64(zerov, resv), mod2v));
                _mm256_storeu_si256((__m256i *)(drs[r]+j), resv);
            }
        }
#elif defined __aarch64__
        uint32x2_t mulv[LA_NRB];
        for (r = 0; r < LA_NRB; ++r) {
            mulv[r] = vdup_n_u32((uint32_t)mul[r]);
        }
        for (; j + 2 <= ncr; j += 2) {
            const uint32x2_t cfv  = vcreate_u32(
                    (uint64_t)piv[j] | ((uint64_t)piv[j+1] << 32));
            for (r = 0; r < LA_NRB; ++r) {
                int64x2_t resv  = vsubq_s64(vld1q_s64(drs[r]+j),
                        vreinterpretq_s64_u64(vmull_u32(mulv[r], cfv)));
                const int64x2_t mask = vreinterpretq_s64_u64(vcltzq_s64(resv));
                resv  = vaddq_s64(resv, vandq_s64(mask, mod2v));
                vst1q_s64(drs[r]+j, resv);
            }
        }
#endif
        for (; j < ncr; ++j) {
            for (r = 0; r < LA_NRB; ++r) {
                drs[r][j] -=  mul[r] * piv[j];
                drs[r][j] +=  (drs[r][j] >> 63) & mod2;
            }
        }
    }
    for (r = 0; r < nr; ++r) {
        for (j = fb; j < ncr; ++j) {
            rs[r][j] = (LA_CF_T)barrett_reduce((uint64_t)drs[r][j], fc, bm);
        }
    }
}

#endif /* LA_WIDE */

/* Interreduces the dense rows dm[0],...,dm[ncr-1] (NULL if there is no
 * pivot in this column), each row starts at its pivot column and has a
 * leading coefficient of 1. The pivots are handled in blocks from the
 * right to the left: the rows of a block are interreduced, afterwards
 * all rows left of the block are reduced by the pivots of the block in
 * parallel. The result is the reduced row echelon form, so it does not
 * depend on the number of threads. */
static LA_CF_T **LA_FF(interreduce_dense_matrix)(
        LA_CF_T **dm,
        const len_t ncr,
        const uint64_t fc,
        const int32_t nthrds
        )
{
    len_t i, n = 0;
    int64_t t;
    const len_t bs = 64;

    len_t *pc = (len_t *)malloc((uint64_t)ncr * sizeof(len_t));
    for (i = 0; i < ncr; ++i) {
        if (dm[i] != NULL) {
            pc[n++] = i;
        }
    }
    int64_t *dr = (int64_t *)malloc(
            (uint64_t)ncr * LA_NRB * nthrds * sizeof(int64_t));

    len_t e = n;
    while (e > 1) {
        const len_t s = e > bs ? e - bs : 0;
        /* interreduce the rows of the block */
        for (t = (int64_t)e-2; t >= (int64_t)s; --t) {
            LA_FF(reduce_dense_rows_by_dense_block)(dm+pc[t], pc+t, 1,
                    pc+t+1, e-1-(len_t)t, dm, ncr, fc, dr);
        }
        /* reduce all rows left of the block by the block, LA_NRB rows
         * at once */
        const int64_t ng = ((int64_t)s + LA_NRB - 1) / LA_NRB;
#pragma omp parallel for num_threads(nthrds) private(t) schedule(dynamic, 2)
        for (t = 0; t < ng; ++t) {
            int64_t *drl = dr + (uint64_t)omp_get_thread_num() * ncr * LA_NRB;
            const len_t f = (len_t)t * LA_NRB;
            const len_t l = f + LA_NRB < s ? LA_NRB : s - f;
            LA_CF_T *rows[LA_NRB];
            for (len_t r = 0; r < l; ++r) {
                rows[r] = dm[pc[f+r]];
            }
            LA_FF(reduce_dense_rows_by_dense_block)(rows, pc+f, l,
                    pc+s, e-s, dm, ncr, fc, drl);
        }
        e = s;
    }
    free(dr);
    free(pc);
    return dm;
}

/* Interreduces the new pivots pivs[ncl],...,pivs[ncl+ncr-1] as dense rows
 * if they are dense enough and stores them in mat->tr from the right to
 * the left. Returns 0 if the pivots are too sparse, they are not changed
 * then. */
static int LA_FF(interreduce_new_pivots_dense)(
        mat_t *mat,
        hm_t **pivs,
        const len_t ncl,
        const len_t ncr,
        const uint64_t fc,
        const int32_t nthrds,
        len_t *npivsp
        )
{
    len_t i, j;
    uint64_t nnz = 0, dsz = 0;
    len_t n = 0;

    for (i = 0; i < ncr; ++i) {
        if (pivs[ncl+i] != NULL) {
            n++;
            nnz += pivs[ncl+i][LENGTH];
            dsz += ncr - i;
        }
    }
    if (n < 2 || 4 * nnz < dsz) {
        return 0;
    }
    LA_CF_T **dm  = (LA_CF_T **)calloc((uint64_t)ncr, sizeof(LA_CF_T *));
    hm_t *hdr     = (hm_t *)malloc((uint64_t)ncr * OFFSET * sizeof(hm_t));
    for (i = 0; i < ncr; ++i) {
        hm_t *row = pivs[ncl+i];
        if (row == NULL) {
            continue;
        }
        memcpy(hdr + (uint64_t)i * OFFSET, row, OFFSET * sizeof(hm_t));
        LA_CF_T *cfs  = mat->LA_CF[row[COEFFS]];
        dm[i]         = (LA_CF_T *)calloc((uint64_t)(ncr-i), sizeof(LA_CF_T));
        const hm_t * const ds = row + OFFSET;
        for (j = 0; j < row[LENGTH]; ++j) {
            dm[i][ds[j]-ncl-i] = cfs[j];
        }
        free(cfs);
        mat->LA_CF[row[COEFFS]] = NULL;
        free(row);
        pivs[ncl+i] = NULL;
    }
    dm = LA_FF(interreduce_dense_matrix)(dm, ncr, fc, nthrds);

    len_t npivs = 0;
    for (i = ncr; i > 0; --i) {
        LA_CF_T *d = dm[i-1];
        if (d == NULL) {
            continue;
        }
        const len_t len = ncr-i+1;
        len_t k = 0;
        for (j = 0; j < len; ++j) {
            k += d[j] != 0;
        }
        hm_t *row     = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
        LA_CF_T *cfs  = (LA_CF_T *)malloc((uint64_t)k * sizeof(LA_CF_T));
        memcpy(row, hdr + (uint64_t)(i-1) * OFFSET, OFFSET * sizeof(hm_t));
        hm_t *rs = row + OFFSET;
        k = 0;
        for (j = 0; j < len; ++j) {
            if (d[j] != 0) {
                rs[k]   = ncl+i-1+j;
                cfs[k]  = d[j];
                k++;
            }
        }
        row[PRELOOP]  = k % UNROLL;
        row[LENGTH]   = k;
        mat->LA_CF[row[COEFFS]] = cfs;
        pivs[ncl+i-1] = mat->tr[npivs++] = row;
        free(d);
    }
    free(dm);
    free(hdr);
    *npivsp = npivs;
    return 1;
}

static void LA_FF(probabilistic_sparse_reduced_echelon_form)(
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
    /* known pivot rows in one block in the order of their lead terms */
    hm_t *kpivs   = NULL;
    if (mat->rr != NULL) {
        kpivs = reducer_rows_to_block(mat);
        memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    }
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->LA_CF[j]      = bs->LA_CF[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    const uint64_t fc   = st->fc;
    /* Why do we generate the random linear combinations so strangely?
     * We can have fc > 2^31, so shifting the result of random multiplier
     * and coefficient may lead to elements of size > 2^63 which cannot
     * be represented by int64_t anymore. Thus, for 32 bit primes,
     * LA_RANDOM_MULTIPLIER restricts the multiplier to < 2^15. Then,
     * after subtracting the multiplied values we have to check if at some
     * point the result gets negative as int64_t, i.e. if drl[i] >> 63 is 1.
     * If so, we add a corresponding multiple of fc. Again, if fc > 2^31 we cannot just
     * add fc * fc, but we have to restrict to fc/2 * fc.
     * Note: We cannot apply this trick to the reduction itself, there we
     * cannot restrict the multiplier to be < 2^31, so we have to handle the
     * fc > 2^31 case differently. */
#ifdef LA_WIDE
    const ff64_t f  = ff64_init(fc);
#else
    uint64_t tmp  = (uint64_t)fc * fc;
    while (tmp > pow(2,63)) {
      tmp -=  (uint64_t)(fc/2) * fc;
    }
    const int64_t mod2  = tmp;
#endif
    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (uint64_t)(ncols) * st->nthrds * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (uint64_t)(rpb) * st->nthrds * sizeof(int64_t));

    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) \
    schedule(dynamic)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)(ncols));
        int64_t *mull = mul + (omp_get_thread_num() * (uint64_t)(rpb));
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            hm_t *npiv  = NULL;
            LA_CF_T *cfs;
            /* starting column, offset, coefficient array position in tmpcf */
            hm_t sc, cfp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                cfp = bctr + i*rpb;
                sc  = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    do {
                        mull[j] = LA_RANDOM_MULTIPLIER(fc);
                    } while (mull[j] == 0);
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    cfs   = bs->LA_CF[npiv[COEFFS]];
                    const len_t os  = npiv[PRELOOP];
                    const len_t len = npiv[LENGTH];
                    const hm_t * const ds = npiv + OFFSET;
                    sc    = sc < ds[0] ? sc : ds[0];
#ifdef LA_WIDE
                    (void)os;
                    (void)l;
                    ff64_submul_sparse(drl, ds, cfs, len, (uint64_t)mull[k], &f);
#else
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
#endif
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
                    npiv  = LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                            drl, mat, bs, pivs, sc, cfp, 0, 0, 0, st);
                    if (!npiv) {
                        bctr  = nrbl;
                        break;
                    }
                    /* normalize coefficient array
                    * NOTE: this has to be done here, otherwise the reduction may
                    * lead to wrong results in a parallel computation since other
                    * threads might directly use the new pivot once it is synced. */
                    if (mat->LA_CF[npiv[COEFFS]][0] != 1) {
                        LA_FF(normalize_sparse_matrix_row)(
                                mat->LA_CF[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    cfs = mat->LA_CF[npiv[COEFFS]];
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
        pivs[i] = NULL;
    }
    free(kpivs);
    kpivs = NULL;

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    LA_CF_T *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc, cfp;
    if (!LA_FF(interreduce_new_pivots_dense)(mat, pivs, ncl, ncr,
                st->fc, st->nthrds, &npivs)) {
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
                cfs = mat->LA_CF[pivs[k][COEFFS]];
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
                    LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                            dr, mat, bs, pivs, sc, cfp, mh, bi, 0, st);
            }
        }
    }
    free(mat->rr);
    mat->rr = NULL;

    free(pivs);
    pivs  = NULL;

    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}

/* Reduces the row given by the columns cols and the coefficients cfs by
 * the pivots like LA_FF(reduce_dense_row_by_known_pivots_sparse), but
 * only visits the columns that are nonzero at some point: they are kept
 * in a min heap and processed in increasing order. For very sparse rows
 * this avoids clearing and scanning the dense row over all columns.
 * dr and mark have to be zero on all ncols columns and are zero again on
 * return, heap and rem need space for ncols entries. The products are
 * accumulated with the correction by fc^2 known from the 31 bit kernels,
 * which works for all supported field characteristics. */
static hm_t *LA_FF(reduce_sparse_row_by_pivots)(
        int64_t *dr,
        uint8_t *mark,
        hm_t *heap,
        hm_t *rem,
        const hm_t * const cols,
        const LA_CF_T * const cfs,
        const len_t len,
        mat_t *mat,
        hm_t * const *pivs,
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        )
{
    len_t j, nh = 0, k = 0;
#ifdef LA_WIDE
    const ff64_t f      = ff64_init(st->fc);
#else
    const int64_t mod   = (int64_t)st->fc;
    const int64_t mod2  = (int64_t)st->fc * st->fc;
    const uint64_t bm   = barrett_multiplier(st->fc);
#endif
    const len_t ncl     = mat->ncl;
    LA_CF_T * const * const mcf = mat->LA_CF;
    rba_t *rba = tr > 0 ? mat->rba[tmp_pos] : NULL;

    for (j = 0; j < len; ++j) {
        dr[cols[j]]   = (int64_t)cfs[j];
        mark[cols[j]] = 1;
        push_column_to_heap(heap, &nh, cols[j]);
    }
    while (nh > 0) {
        const hm_t c  = pop_column_from_heap(heap, &nh);
        mark[c] = 0;
        if (dr[c] == 0) {
            continue;
        }
#ifdef LA_WIDE
        const int64_t mul = dr[c];
        const uint64_t mm = ff64_mont((uint64_t)mul, &f);
#else
        const int64_t mul = (int64_t)barrett_reduce((uint64_t)dr[c], mod, bm);
#endif
        if (mul == 0) {
            dr[c] = 0;
            continue;
        }
        if (pivs[c] == NULL) {
            dr[c]    = mul;
            rem[k++] = c;
            continue;
        }
        if (rba != NULL && c < ncl) {
            rba[c/32] |= 1U << (c % 32);
        }
        /* pivots are monic, their first column is c */
        const hm_t * const ds = pivs[c] + OFFSET;
        const LA_CF_T * const pcfs = mcf[pivs[c][COEFFS]];
        const len_t plen = pivs[c][LENGTH];
        dr[c] = 0;
        for (j = 1; j < plen; ++j) {
            const hm_t d = ds[j];
            if (mark[d] == 0) {
                mark[d] = 1;
                push_column_to_heap(heap, &nh, d);
            }
#ifdef LA_WIDE
            dr[d] = ff64_submul(dr[d], pcfs[j], mm, &f);
#else
            dr[d] -= mul * pcfs[j];
            dr[d] += (dr[d] >> 63) & mod2;
#endif
        }
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    LA_CF_T *cf = (LA_CF_T *)malloc((uint64_t)k * sizeof(LA_CF_T));
    hm_t *rs    = row + OFFSET;
    for (j = 0; j < k; ++j) {
        rs[j]       = rem[j];
        cf[j]       = (LA_CF_T)dr[rem[j]];
        dr[rem[j]]  = 0;
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = k % UNROLL;
    row[LENGTH]   = k;
    mat->LA_CF[tmp_pos] = cf;

    return row;
}

static void LA_FF(exact_sparse_reduced_echelon_form)(
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

    const int32_t nthrds = st->nthrds;

    len_t bad_prime = 0;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->LA_CF[j]      = bs->LA_CF[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }
    /* copy the known pivot rows into one block in the order of their
     * lead terms, the reduction applies them in this order */
    hm_t *kpivs = reducer_rows_to_block(mat);
    if (mat->nru > 0) {
        memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    }

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* very sparse rows are reduced via a heap of their nonzero columns,
     * see LA_FF(reduce_sparse_row_by_pivots) */
    const int final       = st->in_final_reduction_step;
    const int sparse_rows = st->nf == 0;
    uint64_t apl = 1; /* average length of the known pivots */
    if (sparse_rows && mat->nru > 0) {
        for (i = 0; i < mat->nru; ++i) {
            apl += mat->rr[i][LENGTH];
        }
        apl = apl / mat->nru + 1;
    }

    int64_t *dr  = (int64_t *)calloc(
            (uint64_t)ncols * nthrds, sizeof(int64_t));
    uint8_t *mark = NULL;
    hm_t *heap    = NULL;
    hm_t *rem     = NULL;
    if (sparse_rows) {
        mark = (uint8_t *)calloc((uint64_t)ncols * nthrds, sizeof(uint8_t));
        heap = (hm_t *)malloc((uint64_t)ncols * nthrds * sizeof(hm_t));
        rem  = (hm_t *)malloc((uint64_t)ncols * nthrds * sizeof(hm_t));
    }
    /* the rows which fill at least 1/32 of their remaining columns are
     * first reduced by the known pivots in groups of LA_NR rows, see
     * LA_FF(reduce_dense_rows_by_known_pivots), afterwards each remainder
     * is reduced by the new pivots. For sparser rows the reductions of
     * the rows of a group rarely hit the same pivots, then the grouping
     * does not pay off. */
    const len_t gsz = LA_FF(rows_per_group)(nrl, nthrds);
    const int multi_rows = st->nf == 0 && gsz > 0;
    /* when learning the trace the kernel marks the used pivots of each
     * row in its reducer bit array */
    const int learn = st->trace_level == LEARN_TRACER
        && st->in_final_reduction_step == 0;
    /* the buffers for the groups are allocated by the threads when
     * they need them */
    int64_t **drm = NULL;
    if (multi_rows) {
        drm = (int64_t **)calloc((uint64_t)nthrds, sizeof(int64_t *));
    }
    /* without groups each row is handled on its own */
    const len_t gl  = multi_rows ? gsz : 1;
    const len_t ngr = (nrl + gl - 1) / gl;
    len_t g;
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(nthrds) \
    private(g, i, j, k, sc) \
    schedule(dynamic)
    for (g = 0; g < ngr; ++g) {
        if (bad_prime != 0) {
            continue;
        }
        const len_t gs  = g * gl;
        const len_t ge  = gs + gl < nrl ? gs + gl : nrl;
        len_t grp[LA_NR];
        len_t ng        = 0;
        int64_t *dm     = NULL;
        if (multi_rows) {
            sc  = ncols;
            for (i = gs; i < ge; ++i) {
                const hm_t * const r = upivs[i];
                const hm_t * const ds = r + OFFSET;
                const hm_t c0     = final ? first_column(r) : ds[0];
                const uint64_t rc = ncols - c0;
                const len_t len   = r[LENGTH];
                if ((uint64_t)len * 32 < rc || (final && c0 >= ncl)) {
                    continue;
                }
                if (dm == NULL) {
                    const int tid = omp_get_thread_num();
                    if (drm[tid] == NULL) {
                        drm[tid] = (int64_t *)calloc(
                                (uint64_t)ncols * LA_NR, sizeof(int64_t));
                    }
                    dm = drm[tid];
                }
                const LA_CF_T * const cfs = tbr->LA_CF[r[COEFFS]];
                for (j = 0; j < len; ++j) {
                    dm[(uint64_t)LA_NR * ds[j] + ng] = (int64_t)cfs[j];
                }
                if (c0 < sc) {
                    sc = c0;
                }
                grp[ng++] = i;
            }
            if (ng > 0) {
                rba_t *rbas[LA_NR];
                if (learn) {
                    for (j = 0; j < LA_NR; ++j) {
                        rbas[j] = j < ng ? mat->rba[grp[j]] : NULL;
                    }
                }
                LA_FF(reduce_dense_rows_by_known_pivots)(
                        dm, NULL, learn ? rbas : NULL, mat, pivs, sc, st->fc);
            }
        }
        len_t lane = 0;
        for (i = gs; i < ge && bad_prime == 0; ++i) {
            /* is the row already reduced by the known pivots? */
            const int pre = lane < ng && grp[lane] == i;
            /* in the final reduction step a basis element without terms
             * in the columns of the known pivots is not changed by them,
             * it directly becomes a new pivot */
            if (final && pre == 0 && first_column(upivs[i]) >= ncl) {
                hm_t *row       = upivs[i];
                const len_t len = row[LENGTH];
                const hm_t opos = row[COEFFS];
                LA_CF_T *cf     = (LA_CF_T *)malloc(
                        (uint64_t)len * sizeof(LA_CF_T));
                memcpy(cf, tbr->LA_CF[opos], (uint64_t)len * sizeof(LA_CF_T));
                if (cf[0] != 1) {
                    LA_FF(normalize_sparse_matrix_row)(
                            cf, row[PRELOOP], len, st->fc);
                }
                mat->LA_CF[i] = cf;
                row[COEFFS]   = i;
                if (__sync_bool_compare_and_swap(&pivs[row[OFFSET]], NULL, row)) {
                    continue;
                }
                /* the lead term is already a pivot, reduce the row */
                row[COEFFS]   = opos;
                mat->LA_CF[i] = NULL;
                free(cf);
            }
            const uint64_t tn = (uint64_t)omp_get_thread_num() * ncols;
            int64_t *drl  = dr + tn;
            hm_t *npiv      = upivs[i];
            LA_CF_T *cfs     = tbr->LA_CF[npiv[COEFFS]];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            const hm_t * const ds = npiv + OFFSET;
            const uint64_t rc = ncols - ds[0];
            const int sparse  = pre == 0 && sparse_rows
                && (uint64_t)len * 32 < rc && (uint64_t)len * apl * 64 < rc;
            k = 0;
            if (pre) {
                /* take the remainder of the row from its lane, it only
                 * has entries in the columns of new pivots */
                memset(drl, 0, (uint64_t)ncl * sizeof(int64_t));
                for (j = ncl; j < ncols; ++j) {
                    drl[j]  = dm[(uint64_t)LA_NR * j + lane];
                    dm[(uint64_t)LA_NR * j + lane] = 0;
                }
                lane++;
                cfs = NULL;
            } else if (sparse == 0) {
                /* the dense row is zero after a sparse reduction */
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
            }
            int first = 1;
            do {
                /* If we do normal form computations the first monomial in the polynomial might not
                be a known pivot, thus setting it to npiv[OFFSET] can lead to wrong results. */
                sc  = st->nf == 0 ? npiv[OFFSET] : 0;
                if (pre && first) {
                    sc = ncl;
                } else if (final && first) {
                    sc = first_column(npiv);
                }
                hm_t *prev      = npiv;
                LA_CF_T *pcfs   = cfs;
                if (sparse) {
                    npiv = mat->tr[i] = LA_FF(reduce_sparse_row_by_pivots)(
                            drl, mark + tn, heap + tn, rem + tn, prev + OFFSET,
                            pcfs, prev[LENGTH], mat, pivs, i, mh, bi,
                            st->trace_level == LEARN_TRACER, st);
                    /* the first time cfs belongs to tbr */
                    if (first) {
                        pcfs = NULL;
                    }
                } else {
                    npiv  = mat->tr[i] = LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                            drl, mat, bs, pivs, sc, i, mh, bi, st->trace_level == LEARN_TRACER, st);
                }
                first = 0;
                free(prev);
                free(pcfs);
                if (st->nf > 0) {
                    if (!npiv) {
                        mat->tr[i]  = NULL;
                        break;
                    }
                    mat->tr[i]  = npiv;
                    cfs = mat->LA_CF[npiv[COEFFS]];
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
                    if (mat->LA_CF[npiv[COEFFS]][0] != 1) {
                        LA_FF(normalize_sparse_matrix_row)(
                                mat->LA_CF[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                    }
                    k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                    cfs = mat->LA_CF[npiv[COEFFS]];
                }
            } while (!k);
            /* a dense reduction leaves the new pivot's coefficients in the
             * dense row, the sparse reductions need it to be zero */
            if (sparse_rows && sparse == 0 && npiv != NULL) {
                for (j = 0; j < npiv[LENGTH]; ++j) {
                    drl[npiv[OFFSET+j]] = 0;
                }
            }
        }
    }
    if (drm != NULL) {
        for (i = 0; i < (len_t)nthrds; ++i) {
            free(drm[i]);
        }
        free(drm);
    }

    if (bad_prime == 1) {
        for (i = kpivs != NULL ? ncl : 0; i < ncl+ncr; ++i) {
            free(pivs[i]);
            pivs[i] = NULL;
        }
        free(kpivs);
        free(dr);
        free(mark);
        free(heap);
        free(rem);
        mat->np = 0;
        if (st->info_level > 0) {
            fprintf(ERRSTREAM, "Zero reduction while applying tracer, bad prime.\n");
        }
        return;
    }

    /* construct the trace */
    if (st->trace_level == LEARN_TRACER && st->in_final_reduction_step == 0) {
        construct_trace(st->tr, mat);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        if (kpivs == NULL) {
            free(pivs[i]);
        }
        pivs[i] = NULL;
    }
    free(kpivs);
    kpivs = NULL;

    len_t npivs = 0; /* number of new pivots */

    if (st->nf == 0) {
        dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
        mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

        /* interreduce new pivots */
        LA_CF_T *cfs;
        hm_t cf_array_pos;
        int clean = 0; /* is dr zero? */
        if (!LA_FF(interreduce_new_pivots_dense)(mat, pivs, ncl, ncr,
                    st->fc, nthrds, &npivs)) {
            for (i = 0; i < ncr; ++i) {
                k = ncols-1-i;
                if (pivs[k]) {
                    cfs = mat->LA_CF[pivs[k][COEFFS]];
                    cf_array_pos    = pivs[k][COEFFS];
                    const len_t os  = pivs[k][PRELOOP];
                    const len_t len = pivs[k][LENGTH];
                    const len_t bi  = pivs[k][BINDEX];
                    const len_t mh  = pivs[k][MULT];
                    const hm_t * const ds = pivs[k] + OFFSET;
                    sc  = ds[0];
                    const uint64_t rc = ncols - sc;
                    if ((uint64_t)len * 32 < rc && (uint64_t)len * apl * 64 < rc) {
                        if (clean == 0) {
                            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
                            clean = 1;
                        }
                        hm_t *prev  = pivs[k];
                        pivs[k]     = NULL;
                        pivs[k] = mat->tr[npivs++] =
                            LA_FF(reduce_sparse_row_by_pivots)(
                                    dr, mark, heap, rem, prev + OFFSET, cfs, len,
                                    mat, pivs, cf_array_pos, mh, bi, 0, st);
                        free(prev);
                        free(cfs);
                        continue;
                    }
                    clean = 0;
                    memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
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
                        LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                                dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
                }
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
    free(mark);
    free(heap);
    free(rem);
}

static void LA_FF(exact_trace_sparse_reduced_echelon_form)(
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

    const int32_t nthrds = st->in_final_reduction_step == 1 ? 1 : st->nthrds;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    /* known pivot rows in one block in the order of their lead terms */
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

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
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols);
        hm_t *npiv      = upivs[i];
        rba_t *rba      = mat->rba[i];
        LA_CF_T *cfs     = bs->LA_CF[npiv[COEFFS]];
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
            npiv  = mat->tr[i]  = LA_FF(trace_reduce_dense_row_by_known_pivots_sparse)(
                    rba, drl, mat, bs, pivs, sc, i, mh, bi, st);
            if (!npiv) {
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mat->LA_CF[npiv[COEFFS]][0] != 1) {
                LA_FF(normalize_sparse_matrix_row)(
                        mat->LA_CF[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                st->trace_nr_mult  +=  npiv[LENGTH] / 1000.0;
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
            cfs = mat->LA_CF[npiv[COEFFS]];
        } while (!k);
    }

    /* construct the trace */
    construct_trace(trace, mat);

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        pivs[i] = NULL;
    }
    free(kpivs);
    kpivs = NULL;

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    LA_CF_T *cfs;
    hm_t cf_array_pos;
    if (!LA_FF(interreduce_new_pivots_dense)(mat, pivs, ncl, ncr,
                st->fc, nthrds, &npivs)) {
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
                cfs = mat->LA_CF[pivs[k][COEFFS]];
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
                    LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                            dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
            }
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;
}

static int LA_FF(exact_application_sparse_reduced_echelon_form)(
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
    /* known pivot rows in one block in the order of their lead terms */
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

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
            LA_CF_T *cfs     = bs->LA_CF[npiv[COEFFS]];
            const len_t os  = npiv[PRELOOP];
            const len_t len = npiv[LENGTH];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
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
                npiv  = mat->tr[i]  = LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                        drl, mat, bs, pivs, sc, i, mh, bi, 0, st);
                if (!npiv) {
                    fprintf(ERRSTREAM, "Unlucky prime detected, row reduced to zero.");
                    flag  = 0;
                    break;
                }

                /* normalize coefficient array
                 * NOTE: this has to be done here, otherwise the reduction may
                 * lead to wrong results in a parallel computation since other
                 * threads might directly use the new pivot once it is synced. */
                if (mat->LA_CF[npiv[COEFFS]][0] != 1) {
                    LA_FF(normalize_sparse_matrix_row)(
                            mat->LA_CF[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH], st->fc);
                }
                k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                cfs = mat->LA_CF[npiv[COEFFS]];
            } while (!k);
        }
    }
    /* unlucky prime found */
    if (flag == 0) {
        return 1;
    }
    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        pivs[i] = NULL;
    }
    free(kpivs);
    kpivs = NULL;

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (uint64_t)ncols * sizeof(int64_t));
    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    LA_CF_T *cfs;
    hm_t cf_array_pos;
    if (!LA_FF(interreduce_new_pivots_dense)(mat, pivs, ncl, ncr,
                st->fc, nthrds, &npivs)) {
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
                cfs = mat->LA_CF[pivs[k][COEFFS]];
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
                    LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                            dr, mat, bs, pivs, sc, cf_array_pos, mh, bi, 0, st);
            }
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

/* Reduces the rows to be reduced by the known pivots, returns their
 * nonzero remainders as dense rows on the columns of the D part. If rid
 * is not NULL, the trace is learned: the known pivots used for each row
 * are marked in its reducer bit array, the rows of mat->tr are kept for
 * constructing the trace and *rid maps the returned dense rows to the
 * rows of mat->tr. */
static LA_CF_T **LA_FF(sparse_AB_CD_linear_algebra)(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st,
        len_t **rid
        )
{
    len_t i = 0, j;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncl   = mat->ncl;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    mat->LA_CF    = realloc(mat->LA_CF,
            (uint64_t)(nrl + mat->nru) * sizeof(LA_CF_T *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->LA_CF[j]      = bs->LA_CF[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }
    /* known pivot rows in one block in the order of their lead terms */
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* dense rows representing updated D part;
     * after reducing CD part with AB */
    LA_CF_T **drs    = (LA_CF_T **)calloc((uint64_t)nrl, sizeof(LA_CF_T *));

    /* the rows are reduced in groups of gs rows, see
     * LA_FF(reduce_dense_rows_by_known_pivots) */
    const len_t gs = LA_FF(rows_per_group)(nrl, st->nthrds);
    int64_t *dr  = NULL;
    if (gs > 0) {
        dr  = (int64_t *)calloc(
                (uint64_t)ncols * LA_NR * st->nthrds, sizeof(int64_t));
        const len_t ngr = (nrl + gs - 1) / gs;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
        for (i = 0; i < ngr; ++i) {
            int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols * LA_NR);
            const len_t f = i * gs;
            const len_t l = f + gs < nrl ? gs : nrl - f;
            sc = ncols;
            for (len_t r = 0; r < l; ++r) {
                hm_t *npiv              = upivs[f+r];
                const LA_CF_T * const cfs = tbr->LA_CF[npiv[COEFFS]];
                const len_t len         = npiv[LENGTH];
                const hm_t * const ds   = npiv + OFFSET;
                for (j = 0; j < len; ++j) {
                    drl[(uint64_t)LA_NR * ds[j] + r] = (int64_t)cfs[j];
                }
                if (ds[0] < sc) {
                    sc = ds[0];
                }
                if (rid == NULL) {
                    free(npiv);
                    upivs[f+r] = NULL;
                }
            }
            rba_t *rbas[LA_NR];
            if (rid != NULL) {
                for (len_t r = 0; r < LA_NR; ++r) {
                    rbas[r] = r < l ? mat->rba[f+r] : NULL;
                }
            }
            LA_CF_T *rows[LA_NR];
            /* the known pivots have their lead terms left of ncl, so the
             * remainders live in the columns of the D part; drl is zero
             * again afterwards */
            LA_FF(reduce_dense_rows_by_known_pivots)(
                    drl, rows, rid != NULL ? rbas : NULL, mat, pivs, sc,
                    st->fc);
            for (len_t r = 0; r < l; ++r) {
                drs[f+r] = rows[r];
            }
        }
    } else {
        const len_t ncr = mat->ncr;
        dr  = (int64_t *)malloc(
                (uint64_t)ncols * st->nthrds * sizeof(int64_t));
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
        for (i = 0; i < nrl; ++i) {
            LA_CF_T *cfs  = NULL;
            int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncols);
            hm_t *npiv    = upivs[i];
            memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = tbr->LA_CF[npiv[COEFFS]];
            const len_t len = npiv[LENGTH];
            const hm_t * const ds = npiv + OFFSET;
            for (j = 0; j < len; ++j) {
                drl[ds[j]]  = (int64_t)cfs[j];
            }
            sc  = ds[0];
            if (rid == NULL) {
                free(npiv);
                upivs[i] = NULL;
            }
            npiv = LA_FF(reduce_dense_row_by_known_pivots_sparse)(
                    drl, mat, bs, pivs, sc, i, 0, 0, rid != NULL, st);
            if (npiv != NULL) {
                cfs = mat->LA_CF[i];
                LA_CF_T *row = (LA_CF_T *)calloc((uint64_t)ncr, sizeof(LA_CF_T));
                const hm_t * const rs = npiv + OFFSET;
                for (j = 0; j < npiv[LENGTH]; ++j) {
                    row[rs[j]-ncl] = cfs[j];
                }
                free(cfs);
                free(npiv);
                mat->LA_CF[i] = NULL;
                drs[i]  = row;
            }
        }
    }
    free(dr);
    dr  = NULL;

    /* we do not need the old pivots anymore, when learning the trace
     * the caller frees them after constructing the trace, which needs
     * the basis indices and multipliers of the known pivots */
    for (i = 0; i < ncl; ++i) {
        pivs[i] = NULL;
    }
    if (rid != NULL) {
        mat->rrb = kpivs;
    } else {
        free(kpivs);
    }
    kpivs = NULL;
    free(pivs);
    pivs  = NULL;

    /* remove NULL dense rows */
    len_t npivs = 0; /* number of new pivots */
    len_t *ri   = NULL;
    if (rid != NULL) {
        ri = (len_t *)malloc((uint64_t)nrl * sizeof(len_t));
    }
    for (i = 0; i < nrl; ++i) {
        if (drs[i] != NULL) {
            if (ri != NULL) {
                ri[npivs] = i;
            }
            drs[npivs++]  = drs[i];
        }
    }
    if (npivs == 0) {
        free(drs);
        drs = NULL;
    }
    if (rid != NULL) {
        *rid = ri;
    }
    st->np = mat->np = npivs;

    return drs;
}

/* Echelonizes the mat->np dense rows dm. If rid is not NULL, rid[i] is
 * the row of mat->tr the dense row dm[i] comes from and kept[rid[i]] is
 * set to 1 if dm[i] becomes a pivot, i.e. does not reduce to zero. */
static LA_CF_T **LA_FF(exact_dense_linear_algebra)(
        LA_CF_T **dm,
        mat_t *mat,
        md_t *st,
        const len_t * const rid,
        uint8_t *kept
        )
{
    len_t i, j, k, l, npivs;

    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    LA_CF_T **nps  = (LA_CF_T **)calloc((uint64_t)ncr, sizeof(LA_CF_T *));
    /* rows to be further reduced */
    LA_CF_T **tbr  = (LA_CF_T **)calloc((uint64_t)nrows, sizeof(LA_CF_T *));
    /* rows of mat->tr the rows to be further reduced come from */
    len_t *tid     = NULL;
    if (rid != NULL) {
        tid = (len_t *)malloc((uint64_t)nrows * sizeof(len_t));
    }
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
                memmove(dm[i], dm[i]+k, (uint64_t)(ncr-k) * sizeof(LA_CF_T));
                dm[i] = realloc(dm[i], (uint64_t)(ncr-k) * sizeof(LA_CF_T));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = LA_FF(normalize_dense_matrix_row)(nps[k], ncr-k, st->fc);
                }
                if (rid != NULL) {
                    kept[rid[i]] = 1;
                }
            } else {
                if (rid != NULL) {
                    tid[j] = rid[i];
                }
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (uint64_t)ntr * sizeof(LA_CF_T *));
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
        LA_CF_T *npiv  = tbr[i];
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
            npiv = LA_FF(reduce_dense_row_by_dense_new_pivots)(
                    drl, &npc, nps, mat->ncr, st->fc);
            if (npc == -1) {
                break;
            }
            k = __sync_bool_compare_and_swap(&nps[npc], NULL, npiv);
            /* some other thread has already added a pivot so we have to
             * recall the dense reduction process */
        } while (!k);
        if (npc != -1 && rid != NULL) {
            kept[tid[i]] = 1;
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

    free(tbr);
    free(tid);
    free(dr);

    return nps;
}

static LA_CF_T **LA_FF(probabilistic_dense_linear_algebra)(
        LA_CF_T **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, l, m, npivs;

    const uint64_t fc = st->fc;
    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncols = mat->nc;
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    LA_CF_T **nps  = (LA_CF_T **)calloc((uint64_t)ncr, sizeof(LA_CF_T *));
    /* rows to be further reduced */
    LA_CF_T **tbr  = (LA_CF_T **)calloc((uint64_t)nrows, sizeof(LA_CF_T *));

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
                memmove(dm[i], dm[i]+k, (uint64_t)(ncr-k) * sizeof(LA_CF_T));
                dm[i] = realloc(dm[i], (uint64_t)(ncr-k) * sizeof(LA_CF_T));
                nps[k] = dm[i];
                if (nps[k][0] != 1) {
                    nps[k]  = LA_FF(normalize_dense_matrix_row)(nps[k], ncr-k, st->fc);
                }
            } else {
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (uint64_t)ntr * sizeof(LA_CF_T *));

    /* Why do we generate the random linear combinations so strangely?
     * We can have fc > 2^31, so shifting the result of random multiplier
     * and coefficient may lead to elements of size > 2^63 which cannot
     * be represented by int64_t anymore. Thus, for 32 bit primes,
     * LA_RANDOM_MULTIPLIER restricts the multiplier to < 2^15. Then,
     * after subtracting the multiplied values we have to check if at some
     * point the result gets negative as int64_t, i.e. if drl[i] >> 63 is 1.
     * If so, we add a corresponding multiple of fc. Again, if fc > 2^31 we cannot just
     * add fc * fc, but we have to restrict to fc/2 * fc.
     * Note: We cannot apply this trick to the reduction itself, there we
     * cannot restrict the multiplier to be < 2^31, so we have to handle the
     * fc > 2^31 case differently. */
#ifdef LA_WIDE
    const ff64_t f  = ff64_init(fc);
#else
    uint64_t tmp  = (uint64_t)fc * fc;
    while (tmp > pow(2,63)) {
      tmp -=  (uint64_t)(fc/2) * fc;
    }
    const int64_t mod2  = tmp;
#endif

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
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * (uint64_t)ncr);
        int64_t *mull = mul + (omp_get_thread_num() * (uint64_t)rpb);
        const int32_t nbl   = (int32_t) (ntr > (i+1)*rpb ? (i+1)*rpb : ntr);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);

        if (nrbl > 0) {
            hm_t npc;
            hm_t os;
            LA_CF_T *tmp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;
#ifndef LA_WIDE
                os  = ncr % UNROLL;
#endif

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    do {
                        mull[j] = LA_RANDOM_MULTIPLIER(fc);
                    } while (mull[j] == 0);
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (uint64_t)ncr * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    tmp  = tbr[m];
#ifdef LA_WIDE
                    (void)os;
                    (void)l;
                    ff64_submul_dense(drl, tmp, ncr, (uint64_t)mull[k], &f);
#else
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
#endif
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = LA_FF(reduce_dense_row_by_dense_new_pivots)(
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

static LA_CF_T **LA_FF(probabilistic_sparse_dense_echelon_form)(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k, l, m, npivs;

    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncols = mat->nc;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    /* known pivot rows in one block in the order of their lead terms */
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* rows already representing new pivots */
    LA_CF_T **nps  = (LA_CF_T **)calloc((uint64_t)ncr, sizeof(LA_CF_T *));

    const uint64_t fc   = st->fc;
#ifdef LA_WIDE
    const ff64_t f      = ff64_init(fc);
#else
    const int64_t mod2  = (int64_t)fc * fc;
#endif

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
            LA_CF_T *tmp;
            hm_t npc;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    do {
                        mull[j] = (int64_t)rand() % fc;
                    } while (mull[j] == 0);
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (uint64_t)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    tmp   = bs->LA_CF[npiv[COEFFS]];
                    const len_t os  = npiv[PRELOOP];
                    const len_t len = npiv[LENGTH];
                    const hm_t * const ds = npiv + OFFSET;
#ifdef LA_WIDE
                    (void)os;
                    (void)l;
                    ff64_submul_sparse(drl, ds, tmp, len, (uint64_t)mull[k], &f);
#else
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
#endif
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = LA_FF(reduce_dense_row_by_all_pivots)(
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


    free(kpivs);
    kpivs = NULL;
    free(pivs);
    pivs  = NULL;
    free(mul);
    mul = NULL;
    free(dr);
    dr  = NULL;

    return nps;
}

static void LA_FF(convert_to_sparse_matrix_rows)(
        mat_t *mat,
        LA_CF_T * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, j, k, l, m;
    LA_CF_T *cfs;
    hm_t *dts, *dss;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->tr     = realloc(mat->tr, (uint64_t)mat->np * sizeof(hm_t *));
    mat->LA_CF  = realloc(mat->LA_CF,
            (uint64_t)mat->np * sizeof(LA_CF_T *));

    l = 0;
    for (i = 0; i < ncr; ++i) {
        m = ncr-1-i;
        if (dm[m] != NULL) {
            cfs = malloc((uint64_t)(ncr-m) * sizeof(LA_CF_T));
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
            cfs = realloc(cfs, (uint64_t)k * sizeof(LA_CF_T));

            /* link to basis */
            mat->tr[l]    = dts;
            mat->LA_CF[l] = cfs;
            l++;
        }
    }
}

static void LA_FF(probabilistic_sparse_linear_algebra)(
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
    mat->LA_CF = realloc(mat->LA_CF,
            (uint64_t)mat->nr * sizeof(LA_CF_T *));
    LA_FF(probabilistic_sparse_reduced_echelon_form)(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

/* In f4: tbr == bs
in nf: tbr are the polynomials to be reduced w.r.t. bs */
static void LA_FF(exact_sparse_linear_algebra)(
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

    /* if the rows to be reduced are dense on the columns of the D part,
     * their remainders are echelonized as dense rows */
    uint64_t nnz = 0;
    for (len_t i = 0; i < mat->nrl; ++i) {
        nnz += mat->tr[i][LENGTH];
    }
    if (st->nf == 0 && st->in_final_reduction_step == 0 && mat->nrl > 0
            && 2 * nnz >= (uint64_t)mat->nrl * mat->ncr) {
        /* when learning the trace we need to know which rows of mat->tr
         * reduce to new pivots and which known pivots they used */
        const int learn = st->trace_level == LEARN_TRACER;
        const len_t nrl = mat->nrl;
        len_t *rid      = NULL;
        uint8_t *kept   = NULL;
        LA_CF_T **dm    = LA_FF(sparse_AB_CD_linear_algebra)(
                mat, tbr, bs, st, learn ? &rid : NULL);
        if (learn) {
            kept = (uint8_t *)calloc((uint64_t)nrl, sizeof(uint8_t));
        }
        if (mat->np > 0) {
            dm  = LA_FF(exact_dense_linear_algebra)(dm, mat, st, rid, kept);
            dm  = LA_FF(interreduce_dense_matrix)(dm, mat->ncr, st->fc, st->nthrds);
        }
        if (learn) {
            for (len_t i = 0; i < nrl; ++i) {
                if (kept[i] == 0) {
                    free(mat->tr[i]);
                    mat->tr[i] = NULL;
                }
            }
            construct_trace(st->tr, mat);
            free(mat->rrb);
            mat->rrb = NULL;
            for (len_t i = 0; i < nrl; ++i) {
                free(mat->tr[i]);
                mat->tr[i] = NULL;
            }
            free(kept);
            free(rid);
        }
        LA_FF(convert_to_sparse_matrix_rows)(mat, dm);
        if (dm) {
            for (len_t i = 0; i < mat->ncr; ++i) {
                free(dm[i]);
            }
            free(dm);
        }
        st->np = mat->nr = mat->sz = mat->np;
    } else {
        /* allocate temporary storage space for sparse
         * coefficients of all pivot rows */
        mat->LA_CF  = realloc(mat->LA_CF,
                (uint64_t)mat->nr * sizeof(LA_CF_T *));
        LA_FF(exact_sparse_reduced_echelon_form)(mat, tbr, bs, st);
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static int LA_FF(exact_application_sparse_linear_algebra)(
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
    mat->LA_CF = realloc(mat->LA_CF,
            (uint64_t)mat->nr * sizeof(LA_CF_T *));
    ret = LA_FF(exact_application_sparse_reduced_echelon_form)(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }

    return ret;
}

static void LA_FF(exact_trace_sparse_linear_algebra)(
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
    mat->LA_CF  = realloc(mat->LA_CF,
            (uint64_t)mat->nr * sizeof(LA_CF_T *));
    LA_FF(exact_trace_sparse_reduced_echelon_form)(trace, mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static void LA_FF(exact_sparse_dense_linear_algebra)(
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
    LA_CF_T **dm;
    dm  = LA_FF(sparse_AB_CD_linear_algebra)(mat, tbr, bs, st, NULL);
    if (mat->np > 0) {
        dm  = LA_FF(exact_dense_linear_algebra)(dm, mat, st, NULL, NULL);
        dm  = LA_FF(interreduce_dense_matrix)(dm, ncr, st->fc, st->nthrds);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    LA_FF(convert_to_sparse_matrix_rows)(mat, dm);

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
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static void LA_FF(probabilistic_sparse_dense_linear_algebra_2)(
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
    LA_CF_T **dm;
    dm  = LA_FF(sparse_AB_CD_linear_algebra)(mat, tbr, bs, st, NULL);
    if (mat->np > 0) {
        dm  = LA_FF(probabilistic_dense_linear_algebra)(dm, mat, st);
        dm  = LA_FF(interreduce_dense_matrix)(dm, mat->ncr, st->fc, st->nthrds);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    LA_FF(convert_to_sparse_matrix_rows)(mat, dm);

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
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static void LA_FF(probabilistic_sparse_dense_linear_algebra)(
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
    LA_CF_T **dm = NULL;
    mat->np = 0;
    dm      = LA_FF(probabilistic_sparse_dense_echelon_form)(mat, bs, st);
    dm      = LA_FF(interreduce_dense_matrix)(dm, mat->ncr, st->fc, st->nthrds);

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    LA_FF(convert_to_sparse_matrix_rows)(mat, dm);

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
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static void LA_FF(interreduce_matrix_rows)(
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
        fprintf(VERBSTREAM, "                          ");
    }

    mat->tr = realloc(mat->tr, (uint64_t)ncols * sizeof(hm_t *));

    mat->LA_CF  = realloc(mat->LA_CF,
            (uint64_t)ncols * sizeof(LA_CF_T *));
    memset(mat->LA_CF, 0, (uint64_t)ncols * sizeof(LA_CF_T *));
    hm_t **pivs = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    /* copy coefficient arrays from basis in matrix, maybe
     * several rows need the same coefficient arrays, but we
     * cannot share them here. */
    for (i = 0; i < nrows; ++i) {
        pivs[mat->rr[i][OFFSET]]  = mat->rr[i];
    }

    int64_t *dr = (int64_t *)malloc((uint64_t)ncols * sizeof(int64_t));
    /* interreduce new pivots */
    LA_CF_T *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc;
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            memset(dr, 0, (uint64_t)ncols * sizeof(int64_t));
            cfs = bs->LA_CF[pivs[l][COEFFS]];
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
                LA_FF(reduce_dense_row_by_known_pivots_sparse)(
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

#undef LA_CF_T
#undef LA_CF
#undef LA_FF
#undef LA_RANDOM_MULTIPLIER
#undef LA_WIDE
