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

static inline mpz_t *remove_content_of_sparse_matrix_row_qq(
        mpz_t *row,
        const len_t os,
        const len_t len
        )
{
    len_t i;
    long remove = 1;

    mpz_t content;
    /* compute content, i.e. gcd of all coefficients */
    mpz_init_set(content, row[0]);
    for (i = 1; i < len; ++i) {
        mpz_gcd(content, content, row[i]);
        if (mpz_cmp_si(content, 1) == 0) {
            remove = 0;
            break;
        }
    }
    if (remove == 1) {
        /* remove content */
        for (i = 0; i < os; ++i) {
            mpz_divexact(row[i], row[i], content);
        }
        for (; i < len; i += UNROLL) {
            mpz_divexact(row[i], row[i], content);
            mpz_divexact(row[i+1], row[i+1], content);
            mpz_divexact(row[i+2], row[i+2], content);
            mpz_divexact(row[i+3], row[i+3], content);
        }
    }
    mpz_clear(content);

    /* make lead coefficient positive */
    if (mpz_sgn(row[0]) < 0) {
        for (i = 0; i < os; ++i) {
            mpz_neg(row[i], row[i]);
        }
        for (; i < len; i += UNROLL) {
            mpz_neg(row[i], row[i]);
            mpz_neg(row[i+1], row[i+1]);
            mpz_neg(row[i+2], row[i+2]);
            mpz_neg(row[i+3], row[i+3]);
        }
    }
    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_only_ab_qq(
        mpz_t *dr,
        mat_t *mat,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hi_t i, j;
    hm_t *dts;
    mpz_t *cfs;
    int64_t np  = -1;
    const len_t ncols         = mat->nc;
    const len_t ncl           = mat->ncl;

    hm_t *row = NULL;
    mpz_t *cf = NULL;
    len_t rlen  = 0;

    mpz_t mul1, mul2;
    mpz_inits(mul1, mul2, NULL);
    for (i = dpiv; i < ncl; ++i) {
        /* uses mpz_sgn for checking if dr[i] = 0 */
        if (mpz_sgn(dr[i]) == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                row = (hm_t *)malloc(
                        (unsigned long)(ncols-i+OFFSET) * sizeof(hm_t));
                cf  = (mpz_t *)malloc(
                        (unsigned long)(ncols-i) * sizeof(mpz_t));
                np  = i;
            }
            mpz_init(cf[rlen]);
            mpz_swap(cf[rlen], dr[i]);
            row[rlen+OFFSET] = i;
            rlen++;
            continue;
        }
        /* found reducer row, get multiplier */
        dts = pivs[i];
        cfs = mat->cf_ab_qq[dts[COEFFS]];
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds  = dts + OFFSET;

        /* check if lead coefficient of dr is multiple of lead coefficient
         * of cfs, generate corresponding multipliers respectively */
        if (mpz_divisible_p(dr[i], cfs[0]) != 0) {
            mpz_divexact(mul2, dr[i], cfs[0]);
        } else {
            mpz_lcm(mul1, dr[i], cfs[0]);
            mpz_divexact(mul2, mul1, cfs[0]);
            mpz_divexact(mul1, mul1, dr[i]);
            for (j = 0; j < rlen; ++j) {
                mpz_mul(cf[j], cf[j], mul1);
            }
            for (j = i+1; j < ncols; ++j) {
                if (mpz_sgn(dr[j]) != 0) {
                    mpz_mul(dr[j], dr[j], mul1);
                }
            }
        }
        for (j = 0; j < os; ++j) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
            mpz_submul(dr[ds[j+1]], mul2, cfs[j+1]);
            mpz_submul(dr[ds[j+2]], mul2, cfs[j+2]);
            mpz_submul(dr[ds[j+3]], mul2, cfs[j+3]);
        }
    }
    if (rlen != 0) {
        for (i = ncl; i < ncols; ++i) {
            if (mpz_sgn(dr[i]) == 0) {
                continue;
            }
            mpz_init(cf[rlen]);
            mpz_swap(cf[rlen], dr[i]);
            row[rlen+OFFSET] = i;
            rlen++;
        }
        row     = realloc(row, (unsigned long)(rlen+OFFSET) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)rlen * sizeof(mpz_t));
        row[COEFFS]   = tmp_pos;
        row[PRELOOP]  = rlen % UNROLL;
        row[LENGTH]   = rlen;
        mat->cf_ab_qq[tmp_pos]  = cf;
    }
    mpz_clears(mul1, mul2, NULL);
    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_ab_first_qq(
        mpz_t *dr,
        mat_t *mat,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hi_t i, j;
    hm_t *dts;
    mpz_t *cfs;
    int64_t np  = -1;
    const len_t ncols         = mat->nc;
    const len_t ncl           = mat->ncl;
    mpz_t * const * const mcf = mat->cf_qq;

    hm_t *row = NULL;
    mpz_t *cf = NULL;
    len_t rlen  = 0;

    mpz_t mul1, mul2;
    mpz_inits(mul1, mul2, NULL);
    for (i = dpiv; i < ncols; ++i) {
        /* uses mpz_sgn for checking if dr[i] = 0 */
        if (mpz_sgn(dr[i]) == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                row = (hm_t *)malloc(
                        (unsigned long)(ncols-i+OFFSET) * sizeof(hm_t));
                cf  = (mpz_t *)malloc(
                        (unsigned long)(ncols-i) * sizeof(mpz_t));
                np  = i;
            }
            mpz_init(cf[rlen]);
            mpz_swap(cf[rlen], dr[i]);
            row[rlen+OFFSET] = i;
            rlen++;
            continue;
        }
        /* found reducer row, get multiplier */
        dts = pivs[i];
        if (i < ncl) {
            cfs   = mat->cf_ab_qq[dts[COEFFS]];
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds  = dts + OFFSET;

        /* check if lead coefficient of dr is multiple of lead coefficient
         * of cfs, generate corresponding multipliers respectively */
        if (mpz_divisible_p(dr[i], cfs[0]) != 0) {
            mpz_divexact(mul2, dr[i], cfs[0]);
        } else {
            mpz_lcm(mul1, dr[i], cfs[0]);
            mpz_divexact(mul2, mul1, cfs[0]);
            mpz_divexact(mul1, mul1, dr[i]);
            for (j = 0; j < rlen; ++j) {
                mpz_mul(cf[j], cf[j], mul1);
            }
            for (j = i+1; j < ncols; ++j) {
                if (mpz_sgn(dr[j]) != 0) {
                    mpz_mul(dr[j], dr[j], mul1);
                }
            }
        }
        for (j = 0; j < os; ++j) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
            mpz_submul(dr[ds[j+1]], mul2, cfs[j+1]);
            mpz_submul(dr[ds[j+2]], mul2, cfs[j+2]);
            mpz_submul(dr[ds[j+3]], mul2, cfs[j+3]);
        }
    }
    if (rlen != 0) {
        row     = realloc(row, (unsigned long)(rlen+OFFSET) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)rlen * sizeof(mpz_t));
        row[COEFFS]   = tmp_pos;
        row[PRELOOP]  = rlen % UNROLL;
        row[LENGTH]   = rlen;
        mat->cf_qq[tmp_pos]  = cf;
    }
    mpz_clears(mul1, mul2, NULL);
    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_qq(
        mpz_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hi_t i, j;
    hm_t *dts;
    mpz_t *cfs;
    int64_t np  = -1;
    const len_t ncols         = mat->nc;
    const len_t ncl           = mat->ncl;
    mpz_t * const * const mcf = mat->cf_qq;

    hm_t *row = NULL;
    mpz_t *cf = NULL;
    len_t rlen  = 0;

    mpz_t mul1, mul2;
    mpz_inits(mul1, mul2, NULL);
    for (i = dpiv; i < ncols; ++i) {
        /* uses mpz_sgn for checking if dr[i] = 0 */
        if (mpz_sgn(dr[i]) == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                row = (hm_t *)malloc(
                        (unsigned long)(ncols-i+OFFSET) * sizeof(hm_t));
                cf  = (mpz_t *)malloc(
                        (unsigned long)(ncols-i) * sizeof(mpz_t));
                np  = i;
            }
            mpz_init(cf[rlen]);
            mpz_swap(cf[rlen], dr[i]);
            row[rlen+OFFSET] = i;
            rlen++;
            continue;
        }
        /* found reducer row, get multiplier */
        dts = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_qq[dts[COEFFS]];
        } else {
            cfs   = mcf[dts[COEFFS]];
        }
        const len_t os  = dts[PRELOOP];
        const len_t len = dts[LENGTH];
        const hm_t * const ds  = dts + OFFSET;

        /* check if lead coefficient of dr is multiple of lead coefficient
         * of cfs, generate corresponding multipliers respectively */
        if (mpz_divisible_p(dr[i], cfs[0]) != 0) {
            mpz_divexact(mul2, dr[i], cfs[0]);
        } else {
            mpz_lcm(mul1, dr[i], cfs[0]);
            mpz_divexact(mul2, mul1, cfs[0]);
            mpz_divexact(mul1, mul1, dr[i]);
            for (j = 0; j < rlen; ++j) {
                mpz_mul(cf[j], cf[j], mul1);
            }
            for (j = i+1; j < ncols; ++j) {
                if (mpz_sgn(dr[j]) != 0) {
                    mpz_mul(dr[j], dr[j], mul1);
                }
            }
        }
        for (j = 0; j < os; ++j) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
            mpz_submul(dr[ds[j+1]], mul2, cfs[j+1]);
            mpz_submul(dr[ds[j+2]], mul2, cfs[j+2]);
            mpz_submul(dr[ds[j+3]], mul2, cfs[j+3]);
        }
    }
    if (rlen != 0) {
        row     = realloc(row, (unsigned long)(rlen+OFFSET) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)rlen * sizeof(mpz_t));
        row[COEFFS]   = tmp_pos;
        row[PRELOOP]  = rlen % UNROLL;
        row[LENGTH]   = rlen;
        mat->cf_qq[tmp_pos]  = cf;
    }
    mpz_clears(mul1, mul2, NULL);
    return row;
}

static void exact_sparse_reduced_echelon_form_ab_first_qq(
        mat_t *mat,
        const bs_t * const bs,
        const md_t * const st
        )
{
    len_t i = 0, j, k;
    hi_t sc    = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    mpz_t *cfs;
    hm_t cf_array_pos;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    memcpy(pivs, mat->rr, (unsigned long)mat->nru * sizeof(hm_t *));

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;
    mpz_t *dr  = (mpz_t *)malloc(
            (unsigned long)ncols * sizeof(mpz_t));
    for (i = 0; i < ncols; ++i) {
        mpz_init(dr[i]);
    }

    mat->cf_ab_qq[nru-1]  = (mpz_t *)malloc(
            (unsigned long)pivs[nru-1][LENGTH] * sizeof(mpz_t));
    for (i = 0 ; i < pivs[nru-1][LENGTH]; ++i) {
        mpz_init_set(mat->cf_ab_qq[nru-1][i], bs->cf_qq[pivs[nru-1][COEFFS]][i]);
    }
    pivs[nru-1][COEFFS]        = nru-1;
    for (i = 0; i < nru-1; ++i) {
        k = nru-2-i;
        for (j = 0; j < ncols; ++j) {
            mpz_set_si(dr[j], 0);
        }
        cfs = bs->cf_qq[pivs[k][COEFFS]];
        const len_t os  = pivs[k][PRELOOP];
        const len_t len = pivs[k][LENGTH];
        const hm_t * const ds = pivs[k] + OFFSET;
        sc  = ds[0];
        for (j = 0; j < os; ++j) {
            mpz_set(dr[ds[j]], cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_set(dr[ds[j]], cfs[j]);
            mpz_set(dr[ds[j+1]], cfs[j+1]);
            mpz_set(dr[ds[j+2]], cfs[j+2]);
            mpz_set(dr[ds[j+3]], cfs[j+3]);
        }
        free(pivs[k]);
        cfs = NULL;;
        pivs[k] = NULL;
        pivs[k] =
            reduce_dense_row_by_known_pivots_sparse_only_ab_qq(
                    dr, mat, pivs, sc, k);
        remove_content_of_sparse_matrix_row_qq(
                mat->cf_ab_qq[pivs[k][COEFFS]], pivs[k][PRELOOP], pivs[k][LENGTH]);
    }

    const len_t drlen = st->nthrds * ncols;
    dr  = realloc(dr,
            (unsigned long)drlen * sizeof(mpz_t));
    for (i = ncols; i < drlen; ++i) {
        mpz_init(dr[i]);
    }
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) private(i, j, k, sc) schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        mpz_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv  = upivs[i];
        mpz_t *cfs  = bs->cf_qq[npiv[COEFFS]];
        len_t os    = npiv[PRELOOP];
        len_t len   = npiv[LENGTH];
        hm_t * ds   = npiv + OFFSET;
        k = 0;
        /* reset entries to zero */
        for (j = 0; j < ncols; ++j) {
            mpz_set_si(drl[j], 0);
        }
        for (j = 0; j < os; ++j) {
            mpz_set(drl[ds[j]], cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_set(drl[ds[j]], cfs[j]);
            mpz_set(drl[ds[j+1]], cfs[j+1]);
            mpz_set(drl[ds[j+2]], cfs[j+2]);
            mpz_set(drl[ds[j+3]], cfs[j+3]);
        }
        cfs = NULL;
        k   = 1;
        do {
            sc  = npiv[OFFSET];
            os  = npiv[PRELOOP];
            len = npiv[LENGTH];
            ds  = npiv + OFFSET;
            if (k == 0) {
                /* if we redo the row than we only handle elements
                 * with column index >= sc, so we do not need to
                 * reset the other entries in the row to zero */
                for (j = sc; j < ncols; ++j) {
                    mpz_set_si(drl[j], 0);
                }
                for (j = 0; j < os; ++j) {
                    mpz_swap(drl[ds[j]], cfs[j]);
                    mpz_clear(cfs[j]);
                }
                for (; j < len; j += UNROLL) {
                    mpz_swap(drl[ds[j]], cfs[j]);
                    mpz_clear(cfs[j]);
                    mpz_swap(drl[ds[j+1]], cfs[j+1]);
                    mpz_clear(cfs[j+1]);
                    mpz_swap(drl[ds[j+2]], cfs[j+2]);
                    mpz_clear(cfs[j+2]);
                    mpz_swap(drl[ds[j+3]], cfs[j+3]);
                    mpz_clear(cfs[j+3]);
                }
            }
            free(cfs);
            free(npiv);
            npiv  = reduce_dense_row_by_known_pivots_sparse_ab_first_qq(
                    drl, mat, pivs, sc, i);
            if (!npiv) {
                break;
            }
            /* remove content of coefficient array for better usage as reducer
             * later on.
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mpz_cmp_si(mat->cf_qq[npiv[COEFFS]][0], 1) != 0) {
                remove_content_of_sparse_matrix_row_qq(
                        mat->cf_qq[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH]);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
            cfs = mat->cf_qq[npiv[COEFFS]];
        } while (k == 0);
        cfs = NULL;
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        for (j = 0; j < pivs[i][LENGTH]; ++j) {
            mpz_clear(mat->cf_ab_qq[pivs[i][COEFFS]][j]);
        }
        free(mat->cf_ab_qq[pivs[i][COEFFS]]);
        mat->cf_ab_qq[pivs[i][COEFFS]] = NULL;
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    for (i = ncols; i < drlen; ++i) {
        mpz_clear(dr[i]);
    }
    dr      = realloc(dr, (unsigned long)ncols * sizeof(mpz_t));
    mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            for (j = 0; j < ncols; ++j) {
                mpz_set_si(dr[j], 0);
            }
            cfs = mat->cf_qq[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
            const len_t os  = pivs[k][PRELOOP];
            const len_t len = pivs[k][LENGTH];
            const hm_t * const ds = pivs[k] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                mpz_swap(dr[ds[j]], cfs[j]);
                mpz_clear(cfs[j]);
            }
            for (; j < len; j += UNROLL) {
                mpz_swap(dr[ds[j]], cfs[j]);
                mpz_clear(cfs[j]);
                mpz_swap(dr[ds[j+1]], cfs[j+1]);
                mpz_clear(cfs[j+1]);
                mpz_swap(dr[ds[j+2]], cfs[j+2]);
                mpz_clear(cfs[j+2]);
                mpz_swap(dr[ds[j+3]], cfs[j+3]);
                mpz_clear(cfs[j+3]);
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs] =
                reduce_dense_row_by_known_pivots_sparse_qq(
                        dr, mat, bs, pivs, sc, cf_array_pos);
            remove_content_of_sparse_matrix_row_qq(
                    mat->cf_qq[mat->tr[npivs][COEFFS]],
                    mat->tr[npivs][PRELOOP],
                    mat->tr[npivs][LENGTH]);
            npivs++;
        }
    }
    free(pivs);
    pivs  = NULL;
    for (i = 0; i < ncols; ++i) {
        mpz_clear(dr[i]);
    }
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
    mat->np = mat->nr = mat->sz = npivs;
}

static void exact_sparse_reduced_echelon_form_qq(
        mat_t *mat,
        const bs_t * const bs,
        const md_t * const st
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

    const len_t drlen = st->nthrds * ncols;
    mpz_t *dr  = (mpz_t *)malloc(
            (unsigned long)drlen * sizeof(mpz_t));
    for (i = 0; i < drlen; ++i) {
        mpz_init(dr[i]);
    }
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) private(i, j, k, sc) \
    schedule(dynamic)
    for (i = 0; i < nrl; ++i) {
        mpz_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv  = upivs[i];
        mpz_t *cfs  = bs->cf_qq[npiv[COEFFS]];
        len_t os    = npiv[PRELOOP];
        len_t len   = npiv[LENGTH];
        hm_t * ds   = npiv + OFFSET;
        k = 0;
        /* reset entries to zero */
        for (j = 0; j < ncols; ++j) {
            mpz_set_si(drl[j], 0);
        }
        for (j = 0; j < os; ++j) {
            mpz_set(drl[ds[j]], cfs[j]);
        }
        for (; j < len; j += UNROLL) {
            mpz_set(drl[ds[j]], cfs[j]);
            mpz_set(drl[ds[j+1]], cfs[j+1]);
            mpz_set(drl[ds[j+2]], cfs[j+2]);
            mpz_set(drl[ds[j+3]], cfs[j+3]);
        }
        cfs = NULL;
        k   = 1;
        do {
            sc  = npiv[OFFSET];
            os  = npiv[PRELOOP];
            len = npiv[LENGTH];
            ds  = npiv + OFFSET;
            if (k == 0) {
                /* if we redo the row than we only handle elements
                 * with column index >= sc, so we do not need to
                 * reset the other entries in the row to zero */
                for (j = sc; j < ncols; ++j) {
                    mpz_set_si(drl[j], 0);
                }
                for (j = 0; j < os; ++j) {
                    mpz_swap(drl[ds[j]], cfs[j]);
                    mpz_clear(cfs[j]);
                }
                for (; j < len; j += UNROLL) {
                    mpz_swap(drl[ds[j]], cfs[j]);
                    mpz_clear(cfs[j]);
                    mpz_swap(drl[ds[j+1]], cfs[j+1]);
                    mpz_clear(cfs[j+1]);
                    mpz_swap(drl[ds[j+2]], cfs[j+2]);
                    mpz_clear(cfs[j+2]);
                    mpz_swap(drl[ds[j+3]], cfs[j+3]);
                    mpz_clear(cfs[j+3]);
                }
            }
            free(cfs);
            free(npiv);
            npiv  = reduce_dense_row_by_known_pivots_sparse_qq(
                    drl, mat, bs, pivs, sc, i);
            if (!npiv) {
                break;
            }
            /* remove content of coefficient array for better usage as reducer
             * later on.
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mpz_cmp_si(mat->cf_qq[npiv[COEFFS]][0], 1) != 0) {
                remove_content_of_sparse_matrix_row_qq(
                        mat->cf_qq[npiv[COEFFS]], npiv[PRELOOP], npiv[LENGTH]);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
            cfs = mat->cf_qq[npiv[COEFFS]];
        } while (k == 0);
        cfs = NULL;
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    for (i = ncols; i < drlen; ++i) {
        mpz_clear(dr[i]);
    }
    dr      = realloc(dr, (unsigned long)ncols * sizeof(mpz_t));
    mat->tr = realloc(mat->tr, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    mpz_t *cfs;
    hm_t cf_array_pos;
    for (i = 0; i < ncr; ++i) {
        k = ncols-1-i;
        if (pivs[k]) {
            for (j = 0; j < ncols; ++j) {
                mpz_set_si(dr[j], 0);
            }
            cfs = mat->cf_qq[pivs[k][COEFFS]];
            cf_array_pos    = pivs[k][COEFFS];
            const len_t os  = pivs[k][PRELOOP];
            const len_t len = pivs[k][LENGTH];
            const hm_t * const ds = pivs[k] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                mpz_swap(dr[ds[j]], cfs[j]);
                mpz_clear(cfs[j]);
            }
            for (; j < len; j += UNROLL) {
                mpz_swap(dr[ds[j]], cfs[j]);
                mpz_clear(cfs[j]);
                mpz_swap(dr[ds[j+1]], cfs[j+1]);
                mpz_clear(cfs[j+1]);
                mpz_swap(dr[ds[j+2]], cfs[j+2]);
                mpz_clear(cfs[j+2]);
                mpz_swap(dr[ds[j+3]], cfs[j+3]);
                mpz_clear(cfs[j+3]);
            }
            free(pivs[k]);
            free(cfs);
            pivs[k] = NULL;
            pivs[k] = mat->tr[npivs] =
                reduce_dense_row_by_known_pivots_sparse_qq(
                        dr, mat, bs, pivs, sc, cf_array_pos);
            remove_content_of_sparse_matrix_row_qq(
                    mat->cf_qq[mat->tr[npivs][COEFFS]],
                    mat->tr[npivs][PRELOOP],
                    mat->tr[npivs][LENGTH]);
            npivs++;
        }
    }
    free(pivs);
    pivs  = NULL;
    for (i = 0; i < ncols; ++i) {
        mpz_clear(dr[i]);
    }
    free(dr);
    dr  = NULL;

    mat->tr = realloc(mat->tr, (unsigned long)npivs * sizeof(hi_t *));
    mat->np = mat->nr = mat->sz = npivs;
}

static void exact_sparse_linear_algebra_ab_first_qq(
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
    mat->cf_qq  = realloc(mat->cf_qq,
            (unsigned long)mat->nrl * sizeof(mpz_t *));
    mat->cf_ab_qq  = realloc(mat->cf_ab_qq,
            (unsigned long)mat->nru * sizeof(mpz_t *));
    exact_sparse_reduced_echelon_form_ab_first_qq(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}
static void exact_sparse_linear_algebra_qq(
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
    mat->cf_qq  = realloc(mat->cf_qq,
            (unsigned long)mat->nrl * sizeof(mpz_t *));
    exact_sparse_reduced_echelon_form_qq(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void interreduce_matrix_rows_qq(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        int free_basis
        )
{
    len_t i, j, k, l;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;

    /* adjust displaying timings for statistic printout */
    if (st->info_level > 1) {
        printf("                        ");
    }
    mat->tr = realloc(mat->tr, (unsigned long)ncols * sizeof(hm_t *));

    mat->cf_qq  = realloc(mat->cf_qq,
            (unsigned long)ncols * sizeof(mpz_t *));
    memset(mat->cf_qq, 0, (unsigned long)ncols * sizeof(mpz_t *));
    hm_t **pivs = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    /* copy coefficient arrays from basis in matrix, maybe
     * several rows need the same coefficient arrays, but we
     * cannot share them here. */
    for (i = 0; i < nrows; ++i) {
        pivs[mat->rr[i][OFFSET]]  = mat->rr[i];
    }

    mpz_t *dr = (mpz_t *)malloc((unsigned long)ncols * sizeof(mpz_t));
    for (i = 0; i < ncols; ++i) {
        mpz_init(dr[i]);
    }
    /* interreduce new pivots */
    mpz_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc;
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            for (j = 0; j < ncols; ++j) {
                mpz_set_si(dr[j], 0);
            }
            cfs = bs->cf_qq[pivs[l][COEFFS]];
            const len_t os  = pivs[l][PRELOOP];
            const len_t len = pivs[l][LENGTH];
            const hm_t * const ds = pivs[l] + OFFSET;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                mpz_swap(dr[ds[j]], cfs[j]);
            }
            for (; j < len; j += UNROLL) {
                mpz_swap(dr[ds[j]], cfs[j]);
                mpz_swap(dr[ds[j+1]], cfs[j+1]);
                mpz_swap(dr[ds[j+2]], cfs[j+2]);
                mpz_swap(dr[ds[j+3]], cfs[j+3]);
            }
            free(pivs[l]);
            pivs[l] = NULL;
            pivs[l] = mat->tr[k--] =
                reduce_dense_row_by_known_pivots_sparse_qq(
                        dr, mat, bs, pivs, sc, l);
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
    for (i = 0; i < ncols; ++i) {
        mpz_clear(dr[i]);
    }
    free(dr);
}
