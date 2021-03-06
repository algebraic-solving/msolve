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


#include "nf.h"
/* 31-bit implementation only at the moment */
void get_normal_form_matrix(
        const bs_t * const tbr,
        ht_t *bht,
        const len_t start,
        stat_t *st,
        ht_t **shtp,
        hi_t **hcmp,
        mat_t **matp
        )
{
    hi_t *hcm   = *hcmp;
    ht_t *sht   = *shtp;
    mat_t*mat   = *matp;

    /* mul is (0,...,0) */
    exp_t *mul  = (exp_t *)calloc(bht->nv, sizeof(exp_t));

    select_tbr(tbr, mul, start, mat, st, sht, bht, NULL);
    /* mat->nc usually computed inside symbolic preprocessing,
     * here we have to set it by hand; same holds for mat->nrl */
    mat->nrl  = mat->nr;
    mat->nc   = sht->eld-1;
    convert_hashes_to_columns(&hcm, mat, st, sht);
    sort_matrix_rows_decreasing(mat->rr, mat->nru);

    *hcmp = hcm;
    *shtp = sht;
    *matp = mat;
}


int core_nf(
        bs_t **tbrp,
        ht_t **bhtp,
        stat_t **stp,
        const exp_t * const mul,
        const bs_t * const bs
        )
{
    double rt0, rt1;
    rt0 = realtime();

    bs_t *tbr   = *tbrp;
    ht_t *bht   = *bhtp;
    stat_t *st  = *stp;

    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ht_t *sht = initialize_secondary_hash_table(bht, st);

    select_tbr(tbr, mul, 0, mat, st, sht, bht, NULL);

    symbolic_preprocessing(mat, bs, st, sht, NULL, bht);
    if (st->info_level > 1) {
        printf("nf computation data");
    }
    convert_hashes_to_columns(&hcm, mat, st, sht);
    sort_matrix_rows_decreasing(mat->rr, mat->nru);

    /* linear algebra, depending on choice, see set_function_pointers() */
    exact_sparse_linear_algebra_nf_ff_32(mat, tbr, bs, st);
    /* columns indices are mapped back to exponent hashes */
    return_normal_forms_to_basis(
            mat, tbr, bht, sht, hcm, st);

    /* all rows in mat are now polynomials in the basis,
     * so we do not need the rows anymore */
    clear_matrix(mat);

    rt1 = realtime();
    if (st->info_level > 1) {
        printf("%13.2f sec\n", rt1-rt0);
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    /* free and clean up */
    free(hcm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    *tbrp = tbr;
    *bhtp = bht;
    *stp  = st;

    return 1;
}
