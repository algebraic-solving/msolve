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

/* after calling this procedure we have column indices instead of exponent
 * hashes in the polynomials resp. rows. moreover, we have sorted each row
 * by pivots / non-pivots. thus we get already an A|B splicing of the
 * initial matrix. this is a first step for receiving a full GBLA matrix. */
static void convert_multipliers_to_columns(
        hi_t **hcmp,
        bs_t *sat,
        stat_t *st,
        ht_t *ht
        )
{
    hl_t i;

    hi_t *hcm = *hcmp;
    /* clear ht-ev[0] */
    memset(ht->ev[0], 0, (unsigned long)ht->nv * sizeof(exp_t));

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* all elements in the sht hash table represent
     * exactly one column of the matrix */
    hcm = realloc(hcm, (unsigned long)sat->ld * sizeof(hi_t));
    for (i = 0; i < sat->ld; ++i) {
        hcm[i]  = sat->hm[i][MULT];
    }
    sort_r(hcm, (unsigned long)sat->ld, sizeof(hi_t), hcm_cmp, ht);

    /* printf("hcmm\n");
     * for (int ii=0; ii<sat->ld; ++ii) {
     *     printf("hcmm[%d] = %d | idx %u | ", ii, ht->hd[hcm[ii]].idx, hcm[ii]);
     *     for (int jj = 0; jj < ht->nv; ++jj) {
     *         printf("%d ", ht->ev[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    /* store the other direction (hash -> column) */
    for (i = 0; i < sat->ld; ++i) {
        ht->hd[hcm[i]].idx  = (hi_t)i;
    }

    /* map column positions to mul entries*/
    for (i = 0; i < sat->ld; ++i) {
        sat->hm[i][MULT]  =  ht->hd[sat->hm[i][MULT]].idx;
    }
    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    *hcmp = hcm;
}

static void convert_hashes_to_columns_sat(
        hi_t **hcmp,
        mat_t *mat,
        bs_t *sat,
        stat_t *st,
        ht_t *sht
        )
{
    hl_t i;
    hi_t j, k;
    hm_t *row;
    int64_t nterms = 0;

    hi_t *hcm = *hcmp;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t hi;

    const len_t mnr = mat->nr;
    const hl_t esld = sht->eld;
    hd_t *hds       = sht->hd;
    hm_t **rrows    = mat->rr;

    /* all elements in the sht hash table represent
     * exactly one column of the matrix */
    hcm = realloc(hcm, (esld-1) * sizeof(hi_t));
    for (k = 0, j = 0, i = 1; i < esld; ++i) {
        hi  = hds[i].idx;

        hcm[j++]  = i;
        if (hi == 2) {
            k++;
        }
    }
    sort_r(hcm, (unsigned long)j, sizeof(hi_t), hcm_cmp, sht);

    /* printf("hcm\n");
     * for (int ii=0; ii<j; ++ii) {
     *     printf("hcm[%d] = %d | idx %u | deg %u |", ii, hcm[ii], hds[hcm[ii]].idx, sht->ev[hcm[ii]][DEG];
     *     for (int jj = 0; jj < sht->nv; ++jj) {
     *         printf("%d ", sht->ev[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    mat->ncl  = k;
    mat->ncr  = (len_t)esld - 1 - mat->ncl;

    st->num_rowsred +=  sat->ld;

    /* store the other direction (hash -> column) */
    const hi_t ld = (hi_t)(esld - 1);
    for (k = 0; k < ld; ++k) {
        hds[hcm[k]].idx  = (hi_t)k;
    }

    /* map column positions to reducer matrix */
#pragma omp parallel for num_threads(st->nthrds) private(k, j)
    for (k = 0; k < mat->nru; ++k) {
        const len_t os  = rrows[k][PRELOOP];
        const len_t len = rrows[k][LENGTH];
        row = rrows[k] + OFFSET;
        for (j = 0; j < os; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < len; j += UNROLL) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
    }
    for (k = 0; k < mat->nru; ++k) {
        nterms  +=  rrows[k][LENGTH];
    }
    /* map column positions to saturation elements */
#pragma omp parallel for num_threads(st->nthrds) private(k, j)
    for (k = 0; k < sat->ld; ++k) {
        const len_t os  = sat->hm[k][PRELOOP];
        const len_t len = sat->hm[k][LENGTH];
        row = sat->hm[k] + OFFSET;
        for (j = 0; j < os; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < len; j += UNROLL) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
    }
    for (k = 0; k < mat->nrl; ++k) {
        nterms  +=  sat->hm[k][LENGTH];
    }

    /* next we sort each row by the new colum order due
     * to known / unkown pivots */

    /* NOTE: As strange as it may sound, we do not need to sort the rows.
     * When reducing, we copy them to dense rows, there we copy the coefficients
     * at the right place and reduce then. For the reducers itself it is not
     * important in which order the terms are represented as long as the first
     * term is the lead term, which is always true. Once a row is finally reduced
     * it is copied back to a sparse representation, now in the correct term
     * order since it is coming from the correctly sorted dense row. So all newly
     * added elements have all their terms sorted correctly w.r.t. the given
     * monomial order. */

    /* compute density of matrix */
    nterms  *=  100; /* for percentage */
    double density = (double)nterms / (double)mnr / (double)mat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.2f%%", mat->nr + sat->ld, mat->nc, density);
        fflush(stdout);
    }
    *hcmp = hcm;
}


static void convert_hashes_to_columns(
        hi_t **hcmp,
        mat_t *mat,
        stat_t *st,
        ht_t *sht
        )
{
    hl_t i;
    hi_t j, k;
    hm_t *row;
    int64_t nterms = 0;

    hi_t *hcm = *hcmp;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t hi;

    const len_t mnr = mat->nr;
    const hl_t esld = sht->eld;
    hd_t *hds       = sht->hd;
    hm_t **rrows    = mat->rr;
    hm_t **trows    = mat->tr;

    /* all elements in the sht hash table represent
     * exactly one column of the matrix */
    hcm = realloc(hcm, (esld-1) * sizeof(hi_t));
    for (k = 0, j = 0, i = 1; i < esld; ++i) {
        hi  = hds[i].idx;

        hcm[j++]  = i;
        if (hi == 2) {
            k++;
        }
    }
    sort_r(hcm, (unsigned long)j, sizeof(hi_t), hcm_cmp, sht);

    /* printf("hcm\n");
     * for (int ii=0; ii<j; ++ii) {
     *     printf("hcm[%d] = %d | idx %u | deg %u |", ii, hcm[ii], hds[hcm[ii]].idx, sht->ev[hcm[ii]][DEG]);
     *     for (int jj = 0; jj < sht->nv; ++jj) {
     *         printf("%d ", sht->ev[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    mat->ncl  = k;
    mat->ncr  = (len_t)esld - 1 - mat->ncl;

    st->num_rowsred +=  mat->nrl;

    /* store the other direction (hash -> column) */
    const hi_t ld = (hi_t)(esld - 1);
    for (k = 0; k < ld; ++k) {
        hds[hcm[k]].idx  = (hi_t)k;
    }

    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(st->nthrds) private(k, j)
    for (k = 0; k < mat->nru; ++k) {
        const len_t os  = rrows[k][PRELOOP];
        const len_t len = rrows[k][LENGTH];
        row = rrows[k] + OFFSET;
        for (j = 0; j < os; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < len; j += UNROLL) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
    }
    for (k = 0; k < mat->nru; ++k) {
        nterms  +=  rrows[k][LENGTH];
    }
#pragma omp parallel for num_threads(st->nthrds) private(k, j)
    for (k = 0; k < mat->nrl; ++k) {
        const len_t os  = trows[k][PRELOOP];
        const len_t len = trows[k][LENGTH];
        row = trows[k] + OFFSET;
        for (j = 0; j < os; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < len; j += UNROLL) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
    }
    for (k = 0; k < mat->nrl; ++k) {
        nterms  +=  trows[k][LENGTH];
    }

    /* next we sort each row by the new colum order due
     * to known / unkown pivots */

    /* NOTE: As strange as it may sound, we do not need to sort the rows.
     * When reducing, we copy them to dense rows, there we copy the coefficients
     * at the right place and reduce then. For the reducers itself it is not
     * important in which order the terms are represented as long as the first
     * term is the lead term, which is always true. Once a row is finally reduced
     * it is copied back to a sparse representation, now in the correct term
     * order since it is coming from the correctly sorted dense row. So all newly
     * added elements have all their terms sorted correctly w.r.t. the given
     * monomial order. */

    /* compute density of matrix */
    nterms  *=  100; /* for percentage */
    double density = (double)nterms / (double)mnr / (double)mat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.2f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }
    *hcmp = hcm;
}

static void convert_columns_to_hashes(
        bs_t *bs,
        const hi_t * const hcm,
        const hi_t * const hcmm
        )
{
    len_t i, j;

    for (i = 0; i < bs->ld; ++i) {
        if (bs->hm[i] != NULL) {
            for (j = OFFSET; j < bs->hm[i][LENGTH]+OFFSET; ++j) {
                bs->hm[i][j]  = hcm[bs->hm[i][j]];
            }
            bs->hm[i][MULT] = hcmm[bs->hm[i][MULT]];
        }
    }
}

static void add_kernel_elements_to_basis(
        bs_t *sat,
        bs_t *bs,
        bs_t *kernel,
        const ht_t * const ht,
        const hi_t * const hcm,
        stat_t *st
        )
{
    len_t *terms  = (len_t *)calloc((unsigned long)sat->ld, sizeof(len_t));
    len_t nterms  = 0;
    len_t i, j, k;

    len_t ctr       = 0;
    const len_t bld = bs->ld;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, kernel->ld);

    /* we need to sort the kernel elements first (in order to track
     * redundancy correctly) */
    hm_t **rows = (hm_t **)calloc((unsigned long)kernel->ld, sizeof(hm_t *));
    k = 0;
    for (i = 0; i < kernel->ld; ++i) {
        /* printf("kernel[%u] = %p\n", i, kernel->hm[i]); */
        rows[k++]     = kernel->hm[i];
        kernel->hm[i] = NULL;
    }
    sort_matrix_rows_increasing(rows, k);
    /* only for 32 bit at the moment */
    for (i = 0; i < kernel->ld; ++i) {
        bs->cf_32[bld+ctr]              = kernel->cf_32[rows[i][COEFFS]];
        kernel->cf_32[rows[i][COEFFS]]  = NULL;
        bs->hm[bld+ctr]                 = rows[i];
        bs->hm[bld+ctr][COEFFS]         = bld+ctr;
        for (j = OFFSET; j < bs->hm[bld+ctr][LENGTH]+OFFSET; ++j) {
            bs->hm[bld+ctr][j] = hcm[bs->hm[bld+ctr][j]];
            if (nterms != 0) {
                for (int kk = 0; kk < nterms; ++kk) {
                    if (terms[kk] == bs->hm[bld+ctr][j]) {
                        goto next_j;
                    }
                }
            }
            terms[nterms] = bs->hm[bld+ctr][j];
            nterms++;
next_j:
        }
        if (ht->hd[bs->ev[bld+ctr][OFFSET]][DEG] == 0) {
            bs->constant  = 1;
        }
        /* printf("new element from kernel (%u): length %u | ", bld+ctr, bs->hm[bld+ctr][LENGTH]);
         * for (int kk=0; kk<bs->hm[bld+ctr][LENGTH]; ++kk) {
         *     printf("%u | ", bs->cf_32[bld+ctr][kk]);
         *     printf("%u | ", ht->ev[bs->hm[bld+ctr][OFFSET+kk]][DEG]);
         *     for (int jj=0; jj < ht->nv; ++jj) {
         *         printf("%u ", ht->ev[bs->hm[bld+ctr][OFFSET+kk]][jj]);
         *     }
         *     printf(" || ");
         * }
         * printf("\n"); */
        ctr++;
    }
    /* printf("%u of %u terms are used for kernel elements, %.2f\n", nterms, sat->ld, (float)nterms / (float)sat->ld); */
    free(terms);
    free(rows);
    rows  = NULL;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}

static void return_normal_forms_to_basis(
        mat_t *mat,
        bs_t *bs,
        ht_t *bht,
        const ht_t * const sht,
        const hi_t * const hcm,
        stat_t *st
        )
{
    len_t i;

    const len_t np  = mat->np;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np);

    hm_t **rows = mat->tr;

    /* only for 32 bit at the moment */
    for (i = 0; i < np; ++i) {
        if (rows[i] != NULL) {
            insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
            bs->cf_32[bs->ld] = mat->cf_32[rows[i][COEFFS]];
            rows[i][COEFFS]   = bs->ld;
            bs->hm[bs->ld]    = rows[i];
        } else {
            bs->cf_32[bs->ld] = NULL;
            bs->hm[bs->ld]    = NULL;
        }
        bs->lmps[bs->ld]  = bs->ld;
        bs->lml++;
        bs->ld++;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}

static void convert_sparse_matrix_rows_to_basis_elements(
        mat_t *mat,
        bs_t *bs,
        ht_t *bht,
        const ht_t * const sht,
        const hi_t * const hcm,
        stat_t *st
        )
{
    len_t i;

    const len_t bl  = bs->ld;
    const len_t np  = mat->np;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np);

    hm_t **rows = mat->tr;

    switch (st->ff_bits) {
        case 0:
            for (i = 0; i < np; ++i) {
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
                if (bht->ev[rows[i][OFFSET]][DEG] == 0) {
                    bs->constant  = 1;
                }
                bs->cf_qq[bl+i] = mat->cf_qq[rows[i][COEFFS]];

                rows[i][COEFFS] = bl+i;
                bs->hm[bl+i]    = rows[i];
            }
            break;
        case 8:
            for (i = 0; i < np; ++i) {
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
                if (bht->ev[rows[i][OFFSET]][DEG] == 0) {
                    bs->constant  = 1;
                }
                bs->cf_8[bl+i]  = mat->cf_8[rows[i][COEFFS]];
                rows[i][COEFFS] = bl+i;
                bs->hm[bl+i]    = rows[i];
            }
            break;
        case 16:
            for (i = 0; i < np; ++i) {
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
                if (bht->ev[rows[i][OFFSET]][DEG] == 0) {
                    bs->constant  = 1;
                }
                bs->cf_16[bl+i] = mat->cf_16[rows[i][COEFFS]];
                rows[i][COEFFS] = bl+i;
                bs->hm[bl+i]    = rows[i];
            }
            break;
        case 32:
            int ii;
            long power_t;
            exp_t *e  = bht->ev[0];
            for (i = 0; i < np; ++i) {
                power_t = 0;
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
/*                 for (ii = OFFSET; ii < rows[i][LENGTH]+OFFSET; ++ii) {
 *                     if (bht->ev[rows[i][ii]][bht->nv-1] > 0) {
 *                         power_t = power_t > bht->ev[rows[i][ii]][bht->nv-1] ? power_t : bht->ev[rows[i][ii]][bht->nv-1];
 *                     }
 *                     else {
 *                         power_t  =  -1;
 *                         break;
 *                     }
 *                 }
 *                 if (power_t != -1) {
 *                     for (ii = OFFSET; ii < rows[i][LENGTH]+OFFSET; ++ii) {
 *                         memcpy(e, bht->ev[rows[i][ii]], (unsigned long)bht->nv * sizeof(exp_t));
 *                         e[bht->nv-1] -= power_t;
 *                         rows[i][ii] = insert_in_hash_table(e, bht);
 *
 *                     }
 *                     printf("power_t %ld\n", power_t);
 *                 } */
                if (bht->ev[rows[i][OFFSET]][DEG] == 0) {
                    bs->constant  = 1;
                }
                bs->cf_32[bl+i] = mat->cf_32[rows[i][COEFFS]];
                rows[i][COEFFS] = bl+i;
                bs->hm[bl+i]    = rows[i];
                /* printf("new element (%u): length %u | ", bl+i, bs->hm[bl+i][LENGTH]);
                 * for (int kk=0; kk<bs->hm[bl+i][LENGTH]; ++kk) {
                 *     printf("%u | ", bs->cf_32[bl+i][kk]);
                 * for (int jj=0; jj < bht->nv; ++jj) {
                 *     printf("%u ", bht->ev[bs->hm[bl+i][OFFSET+kk]][jj]);
                 * }
                 * printf(" || ");
                 * }
                 * printf("\n"); */
            }
            break;
        default:
            for (i = 0; i < np; ++i) {
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
                if (bht->ev[rows[i][OFFSET]][DEG] == 0) {
                    bs->constant  = 1;
                }
                bs->cf_32[bl+i] = mat->cf_32[rows[i][COEFFS]];
                rows[i][COEFFS] = bl+i;
                bs->hm[bl+i]    = rows[i];
            }
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}

static void convert_sparse_matrix_rows_to_basis_elements_use_sht(
        mat_t *mat,
        bs_t *bs,
        const hi_t * const hcm,
        stat_t *st
        )
{
    len_t i, j;
    hm_t *row;

    const len_t bl  = bs->ld;
    const len_t np  = mat->np;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np);

    hm_t **rows = mat->tr;

    switch (st->ff_bits) {
        case 0:
            for (i = 0; i < np; ++i) {
                row = rows[i];
                const len_t len = rows[i][LENGTH]+OFFSET;
                for (j = OFFSET; j < len; ++j) {
                    row[j]  = hcm[row[j]];
                }
                bs->cf_qq[bl+i] = mat->cf_qq[row[COEFFS]];
                row[COEFFS]     = bl+i;
                bs->hm[bl+i]    = row;
            }
            break;
        case 8:
            for (i = 0; i < np; ++i) {
                row = rows[i];
                const len_t len = rows[i][LENGTH]+OFFSET;
                for (j = OFFSET; j < len; ++j) {
                    row[j]  = hcm[row[j]];
                }
                bs->cf_8[bl+i]  = mat->cf_8[row[COEFFS]];
                row[COEFFS]     = bl+i;
                bs->hm[bl+i]    = row;
            }
            break;
        case 16:
            for (i = 0; i < np; ++i) {
                row = rows[i];
                const len_t len = rows[i][LENGTH]+OFFSET;
                for (j = OFFSET; j < len; ++j) {
                    row[j]  = hcm[row[j]];
                }
                bs->cf_16[bl+i] = mat->cf_16[row[COEFFS]];
                row[COEFFS]     = bl+i;
                bs->hm[bl+i]    = row;
            }
            break;
        case 32:
            for (i = 0; i < np; ++i) {
                row = rows[i];
                const len_t len = rows[i][LENGTH]+OFFSET;
                for (j = OFFSET; j < len; ++j) {
                    row[j]  = hcm[row[j]];
                }
                bs->cf_32[bl+i] = mat->cf_32[row[COEFFS]];
                row[COEFFS]     = bl+i;
                bs->hm[bl+i]    = row;
            }
            break;
        default:
            for (i = 0; i < np; ++i) {
                row = rows[i];
                const len_t len = rows[i][LENGTH]+OFFSET;
                for (j = OFFSET; j < len; ++j) {
                    row[j]  = hcm[row[j]];
                }
                bs->cf_32[bl+i] = mat->cf_32[row[COEFFS]];
                row[COEFFS]     = bl+i;
                bs->hm[bl+i]    = row;
            }
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
