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
        md_t *st,
        ht_t **shtp,
        hi_t **hcmp,
        mat_t **matp
        )
{
    /* hi_t *hcm   = *hcmp; */
    ht_t *sht   = *shtp;
    mat_t*mat   = *matp;

    /* mul is (0,...,0) */
    exp_t *mul  = (exp_t *)calloc(bht->nv, sizeof(exp_t));

    select_tbr(tbr, mul, start, mat, st, sht, bht, NULL);
    /* mat->nc usually computed inside symbolic preprocessing,
     * here we have to set it by hand; same holds for mat->nrl */
    mat->nrl  = mat->nr;
    mat->nc   = sht->eld-1;
    convert_hashes_to_columns(mat, st, sht);
    sort_matrix_rows_decreasing(mat->rr, mat->nru);

    /* *hcmp = hcm; */
    *shtp = sht;
    *matp = mat;
}


bs_t *core_nf(
        bs_t *tbr,
        md_t *md,
        const exp_t * const mul,
        bs_t *bs,
        int32_t *errp
        )
{
    double ct = cputime();
    double rt = realtime();

    ht_t *bht = bs->ht;

    /* reset to exact linear algebra for normal form computation */
    md->laopt = 2;
    set_function_pointers(md);

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    md->hcm = (hi_t *)malloc(sizeof(hi_t));
    md->ht  = initialize_secondary_hash_table(bht, md);

    md->nf = 1;
    select_tbr(tbr, mul, 0, mat, md, md->ht, bht, NULL);

    symbolic_preprocessing(mat, bs, md);
    convert_hashes_to_columns(mat, md, md->ht);
    sort_matrix_rows_decreasing(mat->rr, mat->nru);

    /* linear algebra, depending on choice, see set_function_pointers() */
    linear_algebra(mat, tbr, bs, md);
    /* columns indices are mapped back to exponent hashes */
    return_normal_forms_to_basis(
            mat, tbr, bht, md->ht, md->hcm, md);

    /* all rows in mat are now polynomials in the basis,
     * so we do not need the rows anymore */
    clear_matrix(mat);

    print_round_timings(stdout, md, rt, ct);
    print_round_information_footer(stdout, md);

    /* free and clean up */
    free(md->hcm);
    if (md->ht != NULL) {
        free_hash_table(&(md->ht));
    }
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    *errp = 0;

    return tbr;
}

int64_t export_nf(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *nf_ld,   /* basis load */
        int32_t **nf_len, /* length of each poly in basis */
        int32_t **nf_exp, /* basis exponent vectors */
        void **nf_cf,     /* coefficients of basis elements */
        /* input values */
        const int32_t nr_tbr_gens,
        const int32_t *tbr_lens,
        const int32_t *tbr_exps,
        const void *tbr_cfs,
        const int32_t nr_bs_gens,
        const int32_t *bs_lens,
        const int32_t *bs_exps,
        const void *bs_cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t bs_is_gb,
        const int32_t nr_threads,
        const int32_t info_level
        )
{
    /* timings */
    double ct = cputime();
    double rt = realtime();

    /* data structures for basis, hash table and statistics */
    bs_t *bs  = NULL;
    bs_t *tbr = NULL;
    md_t *md  = NULL;
    ht_t *bht = NULL;

    int32_t err = 0;
    int success = 0;

    /* get polys w.r.t. which we reduce */
    success = initialize_gba_input_data(&bs, &bht, &md,
            bs_lens, bs_exps, bs_cfs, field_char, mon_order, elim_block_len,
            nr_vars, nr_bs_gens, 0, 17,
            nr_threads, 0, 0, 44, 0,
            1, 0, info_level);

    /* all input generators are invalid */
    if (success == -1) {
        return_zero(nf_ld, nf_len, nf_exp, nf_cf, nr_vars, field_char, mallocp);
        return 1;
    }
    if (success == 0) {
        printf("Bad input data, stopped computation.\n");
        exit(1);
    }
    if (bs_is_gb == 1) {
        for (len_t k = 0; k < bs->ld; ++k) {
            bs->lmps[k] = k;
            bs->lm[k]   = bht->hd[bs->hm[k][OFFSET]].sdm;
            bs->lml     = bs->ld;
        }
    } else {
        /* compute a gb for initial generators */
        bs = core_gba(bs, md, &err, md->fc);

        if (err) {
            printf("Problem with F4, stopped computation.\n");
            exit(1);
        }
    }
    /* initialize data for elements to be reduced,
     * NOTE: Don't initialize BEFORE running core_f4, bht may
     * change, so hash values of tbr may become wrong. */
    tbr = initialize_basis(md);
    tbr->ht = bht;
    import_input_data(tbr, md, 0, nr_tbr_gens,
            tbr_lens, tbr_exps, (void *)tbr_cfs, NULL);
    tbr->ld = tbr->lml  =  nr_tbr_gens;

    /* generate array for storing multiplier for polynomial
     * to be reduced by basis */
    exp_t *mul  = (exp_t *)calloc(bht->evl, sizeof(exp_t));
    /* compute normal form of last element in tbr */
    tbr = core_nf(tbr, md, mul, bs, &err);

    if (err) {
        printf("Problem with normalform, stopped computation.\n");
        exit(1);
    }

    int64_t nterms  = export_results_from_f4(nf_ld, nf_len, nf_exp,
            nf_cf, mallocp, &tbr, &bht, &md);

    /* timings */
    md->nf_ctime = cputime() - ct;
    md->nf_rtime = realtime() - rt;

    get_and_print_final_statistics(stderr, md, tbr);

    /* free and clean up */
    free_shared_hash_data(bht);
    if (tbr != NULL) {
        free_basis_without_hash_table(&tbr);
    }
    if (bs != NULL) {
        free_basis(&bs);
    }
    free(md);
    md    = NULL;

    return nterms;
}
