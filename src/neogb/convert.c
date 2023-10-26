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
        md_t *st,
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
        mat_t *mat,
        bs_t *sat,
        md_t *st,
        ht_t *sht
        )
{
    hl_t i;
    hi_t j, k;
    hm_t *row;
    int64_t nterms = 0;

    hi_t *hcm = st->hcm;

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
     *     printf("hcm[%d] = %d | idx %u | deg %u |", ii, hcm[ii], hds[hcm[ii]].idx, sht->ev[hcm[ii]][DEG]+sht->ev[hcm[ii]][sht->ebl]);
     *     for (int jj = 0; jj < sht->evl; ++jj) {
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
#pragma omp parallel for num_threads(st->nthrds) private(k, j, row)
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
#pragma omp parallel for num_threads(st->nthrds) private(k, j, row)
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
    st->hcm = hcm;
}


static void sba_convert_hashes_to_columns(
        hi_t **hcmp,
        smat_t *smat,
        md_t *st,
        ht_t *ht
        )
{
    len_t i, j, k;

    hm_t *row;
    int64_t nterms = 0;

    hi_t *hcm = *hcmp;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t nr = smat->cld;
    const hl_t eld = ht->eld;
    hd_t *hd       = ht->hd;
    hm_t **cr      = smat->cr;

    hcm = realloc(hcm, (unsigned long)eld * sizeof(hi_t));
    k = 0;
    for (i = 0; i < nr; ++i) {
        const len_t len = SM_OFFSET + cr[i][SM_LEN];
        for (j = SM_OFFSET; j < len; ++j) {
            if (hd[cr[i][j]].idx == 0) {
                hd[cr[i][j]].idx = 1;
                hcm[k++] = cr[i][j];
            }
        }
    }

    hcm = realloc(hcm, (unsigned long)k * sizeof(hi_t));
    sort_r(hcm, (unsigned long)k, sizeof(hi_t), hcm_cmp, ht);

    smat->nc = k;

    /* printf("hcm\n");
     * for (int ii=0; ii<j; ++ii) {
     *     printf("hcm[%d] = %d | idx %u | deg %u |", ii, hcm[ii], hds[hcm[ii]].idx, sht->ev[hcm[ii]][DEG]+sht->ev[hcm[ii]][sht->ebl]);
     *     for (int jj = 0; jj < sht->evl; ++jj) {
     *         printf("%d ", sht->ev[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    /* store the other direction (hash -> column) */
    const hi_t ld = k;
    for (i = 0; i < ld; ++i) {
        hd[hcm[i]].idx = (hi_t)i;
    }

    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(st->nthrds) private(k, j, row)
    for (i = 0; i < nr; ++i) {
        const len_t os  = cr[i][SM_PRE];
        const len_t len = cr[i][SM_LEN];
        row = cr[i] + SM_OFFSET;
        for (j = 0; j < os; ++j) {
            row[j]  = hd[row[j]].idx;
        }
        for (; j < len; j += UNROLL) {
            row[j]    = hd[row[j]].idx;
            row[j+1]  = hd[row[j+1]].idx;
            row[j+2]  = hd[row[j+2]].idx;
            row[j+3]  = hd[row[j+3]].idx;
        }
        nterms += len;
    }

    /* compute density of matrix */
    nterms  *=  100; /* for percentage */
    double density = (double)nterms / (double)nr / (double)smat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf("%4d    %7d x %-7d %8.2f%%", smat->cd, smat->cld, smat->nc, density);
        fflush(stdout);
    }
    *hcmp = hcm;
}


static void convert_hashes_to_columns(
        mat_t *mat,
        md_t *st,
        ht_t *sht
        )
{
    hl_t i;
    hi_t j, k;
    hm_t *row;
    int64_t nterms = 0;

    hi_t *hcm = st->hcm;

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
    for (int ii=0; ii<j; ++ii) {
        printf("hcm[%d] = %d | idx %u | deg %u |", ii, hcm[ii], hds[hcm[ii]].idx, sht->ev[hcm[ii]][DEG]+sht->ev[hcm[ii]][sht->ebl]);
        for (int jj = 0; jj < sht->evl; ++jj) {
            printf("%d ", sht->ev[hcm[ii]][jj]);
        }
        printf("\n");
    } */

    mat->ncl  = k;
    mat->ncr  = (len_t)esld - 1 - mat->ncl;

    st->num_rowsred +=  mat->nrl;

    /* store the other direction (hash -> column) */
    const hi_t ld = (hi_t)(esld - 1);
    for (k = 0; k < ld; ++k) {
        hds[hcm[k]].idx  = (hi_t)k;
    }


    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(st->nthrds) private(k, j, row)
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
#pragma omp parallel for num_threads(st->nthrds) private(k, j, row)
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
    double density = (double)nterms / (double)mnr / (double)ld;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.2f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }
    if ((int64_t)mat->nr * mat->nc > st->mat_max_nrows * st->mat_max_ncols) {
        st->mat_max_nrows = mat->nr;
        st->mat_max_ncols = mat->nc;
        st->mat_max_density = density;
    }
    st->hcm = hcm;
}

static void sba_convert_columns_to_hashes(
        smat_t *smat,
        const hi_t * const hcm
        )
{
    len_t i, j;

    for (i = 0; i < smat->cld; ++i) {
        const len_t len = smat->cr[i][SM_LEN] + SM_OFFSET;
        for (j = SM_OFFSET; j < len; ++j) {
            smat->cr[i][j] = hcm[smat->cr[i][j]];
        }
    }
}

static void convert_columns_to_hashes(
        bs_t *bs,
        const md_t *md,
        const hi_t * const hcmm
        )
{
    len_t i, j;

    const hi_t *hcm = md->hcm;

    for (i = 0; i < bs->ld; ++i) {
        if (bs->hm[i] != NULL) {
            for (j = OFFSET; j < bs->hm[i][LENGTH]+OFFSET; ++j) {
                bs->hm[i][j]  = hcm[bs->hm[i][j]];
            }
            bs->hm[i][MULT] = hcmm[bs->hm[i][MULT]];
        }
    }
}

/* add_kernel_elements_to_basis() is unused at the moment */
#if 0
static void add_kernel_elements_to_basis(
        bs_t *sat,
        bs_t *bs,
        bs_t *kernel,
        const ht_t * const ht,
        const hi_t * const hcm,
        md_t *st
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
        j = OFFSET;
next_j:
        for (; j < bs->hm[bld+ctr][LENGTH]+OFFSET; ++j) {
            bs->hm[bld+ctr][j] = hcm[bs->hm[bld+ctr][j]];
            if (nterms != 0) {
                for (int kk = 0; kk < nterms; ++kk) {
                    if (terms[kk] == bs->hm[bld+ctr][j]) {
                        j++;
                        goto next_j;
                    }
                }
            }
            terms[nterms] = bs->hm[bld+ctr][j];
            nterms++;
        }
        if (ht->ev[bs->hm[bld+ctr][OFFSET]][DEG] == 0) {
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
#endif

static void return_normal_forms_to_basis(
        mat_t *mat,
        bs_t *bs,
        ht_t *bht,
        const ht_t * const sht,
        const hi_t * const hcm,
        md_t *st
        )
{
    len_t i;

    const len_t np  = mat->np;

    /* timings */

    double ct = cputime();
    double rt = realtime();

    free_basis_elements(bs);

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np, st);

    hm_t **rows = mat->tr;

    /* only for 32 bit at the moment */
    for (i = 0; i < np; ++i) {
        if (rows[i] != NULL) {
            insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm, st);
            switch (st->ff_bits) {
                case 0:
                    bs->cf_qq[bs->ld] = mat->cf_qq[rows[i][COEFFS]];
                    break;
                case 8:
                    bs->cf_8[bs->ld]  = mat->cf_8[rows[i][COEFFS]];
                    break;
                case 16:
                    bs->cf_16[bs->ld] = mat->cf_16[rows[i][COEFFS]];
                    break;
                case 32:
                    bs->cf_32[bs->ld] = mat->cf_32[rows[i][COEFFS]];
                    break;
                default:
                    bs->cf_32[bs->ld] = mat->cf_32[rows[i][COEFFS]];
                    break;
            }
            rows[i][COEFFS]   = bs->ld;
            bs->hm[bs->ld]    = rows[i];
        } else {
            switch (st->ff_bits) {
                case 0:
                    bs->cf_qq[bs->ld] = NULL;
                    break;
                case 8:
                    bs->cf_8[bs->ld]  = NULL;
                    break;
                case 16:
                    bs->cf_16[bs->ld] = NULL;
                    break;
                case 32:
                    bs->cf_32[bs->ld] = NULL;
                    break;
                default:
                    bs->cf_32[bs->ld] = NULL;
                    break;
            }
            bs->hm[bs->ld]    = NULL;
        }
        bs->lmps[bs->ld]  = bs->ld;
        bs->lml++;
        bs->ld++;
    }

    /* timings */
    st->convert_ctime +=  cputime() - ct;
    st->convert_rtime +=  realtime() - rt;
}

static void convert_sparse_matrix_rows_to_basis_elements(
        const int sort,
        mat_t *mat,
        bs_t *bs,
        ht_t *bht,
        const ht_t * const sht,
        md_t *st
        )
{
    len_t i, j, k;
    deg_t deg;

    const len_t bl = bs->ld;
    const len_t np = mat->np;
    hi_t *hcm      = st->hcm;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np, st);

    hm_t **rows = mat->tr;
    switch_hcm_data_to_basis_hash_table(hcm, bht, mat, sht);
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k)
    for (k = 0; k < np; ++k) {
        /* We first insert the highest leading monomial element to the basis
         * for a better Gebauer-Moeller application when updating the pair
         * set later on. */
        if (sort == -1) {
            i   =   np - 1 - k;
        } else {
            i = k;
        }
        const len_t len = rows[i][LENGTH]+OFFSET;
        for (j = OFFSET; j < len; ++j) {
            rows[i][j] = hcm[rows[i][j]];
        }
        deg = bht->hd[rows[i][OFFSET]].deg;
        if (st->nev > 0) {
            const len_t len = rows[i][LENGTH]+OFFSET;
            for (j = OFFSET+1; j < len; ++j) {
                if (deg < bht->hd[rows[i][j]].deg) {
                    deg = bht->hd[rows[i][j]].deg;
                }
            }
        }
        switch (st->ff_bits) {
            case 0:
                bs->cf_qq[bl+k] = mat->cf_qq[rows[i][COEFFS]];
                break;
            case 8:
                bs->cf_8[bl+k]  = mat->cf_8[rows[i][COEFFS]];
                break;
            case 16:
                bs->cf_16[bl+k] = mat->cf_16[rows[i][COEFFS]];
                break;
            case 32:
                bs->cf_32[bl+k] = mat->cf_32[rows[i][COEFFS]];
                break;
            default:
                bs->cf_32[bl+k] = mat->cf_32[rows[i][COEFFS]];
                break;
        }
        rows[i][COEFFS]   = bl+k;
        bs->hm[bl+k]      = rows[i];
        bs->hm[bl+k][DEG] = deg;
        if (deg == 0) {
            bs->constant  = 1;
        }
#if 0
        if (st->ff_bits == 32) {
            printf("new element (%u): length %u | degree %d | ", bl+k, bs->hm[bl+k][LENGTH], bs->hm[bl+k][DEG]);
            int kk = 0;
            for (int kk=0; kk<bs->hm[bl+k][LENGTH]; ++kk) {
            printf("%u | ", bs->cf_32[bl+k][kk]);
            for (int jj=0; jj < bht->evl; ++jj) {
                printf("%u ", bht->ev[bs->hm[bl+k][OFFSET+kk]][jj]);
            }
            printf(" || ");
            }
            printf("\n");
        }
        if (st->ff_bits == 16) {
            printf("new element (%u): length %u | degree %d (difference %d) | ", bl+k, bs->hm[bl+k][LENGTH], bs->hm[bl+k][DEG],
                    bs->hm[bl+k][DEG] - bht->hd[bs->hm[bl+k][OFFSET]].deg);
            int kk = 0;
            for (int kk=0; kk<bs->hm[bl+k][LENGTH]; ++kk) {
            printf("%u | ", bs->cf_16[bl+k][kk]);
            for (int jj=0; jj < bht->evl; ++jj) {
                printf("%u ", bht->ev[bs->hm[bl+k][OFFSET+kk]][jj]);
            }
            printf(" || ");
            }
            printf("\n");
        }
        if (st->ff_bits == 8) {
            printf("new element (%u): length %u | degree %d (difference %d) | ", bl+k, bs->hm[bl+k][LENGTH], bs->hm[bl+k][DEG],
                    bs->hm[bl+k][DEG] - bht->hd[bs->hm[bl+k][OFFSET]].deg);
            int kk = 0;
            for (int kk=0; kk<bs->hm[bl+k][LENGTH]; ++kk) {
            printf("%u | ", bs->cf_8[bl+k][kk]);
            for (int jj=0; jj < bht->evl; ++jj) {
                printf("%u ", bht->ev[bs->hm[bl+k][OFFSET+kk]][jj]);
            }
            printf(" || ");
            }
            printf("\n");
        }
#endif
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}

static void convert_sparse_matrix_rows_to_basis_elements_use_sht(
        const int sort,
        mat_t *mat,
        bs_t *bs,
        const ht_t * const sht,
        md_t *st
        )
{
    len_t i, j, k;
    deg_t deg;
    hm_t *row;

    const len_t bl  = bs->ld;
    const len_t np  = mat->np;
    const hi_t * const hcm = st->hcm;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np, st);

    hm_t **rows = mat->tr;

    for (k = 0; k < np; ++k) {
        /* We first insert the highest leading monomial element to the basis
         * for a better Gebauer-Moeller application when updating the pair
         * set later on. */
        if (sort == -1) {
            i   =   np - 1 - k;
        } else {
            i = k;
        }
        row = rows[i];
        deg = sht->hd[hcm[rows[i][OFFSET]]].deg;
        const len_t len = rows[i][LENGTH]+OFFSET;
        if (st->nev ==  0) {
            for (j = OFFSET; j < len; ++j) {
                row[j]  = hcm[row[j]];
            }
        } else {
            for (j = OFFSET; j < len; ++j) {
                row[j]  = hcm[row[j]];
                if (deg < sht->hd[row[j]].deg) {
                    deg = sht->hd[row[j]].deg;
                }
            }
        }
        switch (st->ff_bits) {
            case 0:
                bs->cf_qq[bl+k] = mat->cf_qq[rows[i][COEFFS]];
                break;
            case 8:
                bs->cf_8[bl+k]  = mat->cf_8[rows[i][COEFFS]];
                break;
            case 16:
                bs->cf_16[bl+k] = mat->cf_16[rows[i][COEFFS]];
                break;
            case 32:
                bs->cf_32[bl+k] = mat->cf_32[rows[i][COEFFS]];
                break;
            default:
                bs->cf_32[bl+k] = mat->cf_32[rows[i][COEFFS]];
                break;
        }
        rows[i][COEFFS]   = bl+k;
        bs->hm[bl+k]      = rows[i];
        bs->hm[bl+k][DEG] = deg;
        if (deg == 0) {
            bs->constant  = 1;
        }
#if 0
        if (st->ff_bits == 32) {
            printf("new element (%u): length %u | degree %d | ", bl+k, bs->hm[bl+k][LENGTH], bs->hm[bl+k][DEG]);
            int kk = 0;
            /* for (int kk=0; kk<bs->hm[bl+k][LENGTH]; ++kk) { */
            printf("%u | ", bs->cf_32[bl+k][kk]);
            for (int jj=0; jj < sht->evl; ++jj) {
                printf("%u ", sht->ev[bs->hm[bl+k][OFFSET+kk]][jj]);
            }
            /* printf(" || ");
             * } */
            printf("\n");
        }
#endif
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
