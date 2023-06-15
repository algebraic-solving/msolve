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


#include "f4.h"


void static final_remove_redundant_elements(
        bs_t *bs,
        const ht_t * const ht
        )
{
    len_t i, j;
    for (i = 0; i < bs->lml; ++i) {
        hm_t nch = bs->hm[bs->lmps[i]][OFFSET];
        deg_t dd = bs->hm[bs->lmps[i]][DEG] - ht->hd[nch].deg;
        for (j = 0; j < i; ++j) {
            if (bs->red[bs->lmps[j]] == 0
                    && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                    ) {
                bs->red[bs->lmps[i]]  =   1;
                break;
            }
        }
        for (j = i+1; j < bs->lml; ++j) {
            if (bs->red[bs->lmps[j]] == 0
                    && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                    ) {
                bs->red[bs->lmps[i]]  =   1;
                break;
            }
        }
    }
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;
}

/* The parameters themselves are handled by julia, thus we only
 * free what they are pointing to, julia's garbage collector then
 * takes care of everything leftover. */
void free_f4_julia_result_data(
        void (*freep) (void *),
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,      /* coefficients of basis elements */
        const int64_t ngens,
        const int64_t field_char
        )
{
    /* lengths resp. nterms */
    int32_t *lens  = *blen;

    /* int64_t i;
     * int64_t len = 0;
     * for (i = 0; i < ngens; ++i) {
     *     len += (int64_t)lens[i];
     * } */

    (*freep)(lens);
    lens  = NULL;
    *blen = lens;

    /* exponent vectors */
    int32_t *exps = *bexp;
    (*freep)(exps);
    exps  = NULL;
    *bexp = exps;

    /* coefficients */
    if (field_char == 0) {
        /* mpz_t **cfs = (mpz_t **)bcf;
         * for (i = 0; i < len; ++i) {
         *     mpz_clear((*cfs)[i]);
         * }
         * (*freep)(*cfs);
         * *cfs  = NULL; */
    } else {
        if (field_char > 0) {
            int32_t *cfs  = *((int32_t **)bcf);
            (*freep)(cfs);
            cfs = NULL;
        }
    }
    *bcf  = NULL;
}

static void clear_matrix(
        mat_t *mat
        )
{
    len_t i;
    for (i = 0; i < mat->rbal; ++i) {
        free(mat->rba[i]);
    }
    free(mat->rba);
    mat->rba  = NULL;
    free(mat->rr);
    mat->rr = NULL;
    free(mat->tr);
    mat->tr  = NULL;
    free(mat->cf_8);
    mat->cf_8 = NULL;
    free(mat->cf_16);
    mat->cf_16  = NULL;
    free(mat->cf_32);
    mat->cf_32  = NULL;
    free(mat->cf_qq);
    mat->cf_qq  = NULL;
    free(mat->cf_ab_qq);
    mat->cf_ab_qq  = NULL;
}

#if 0
static void intermediate_reduce_basis(
        bs_t *bs,
        mat_t *mat,
        hi_t **hcmp,
        ht_t **bhtp,
        ht_t **shtp,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j;

    ht_t *bht   = *bhtp;
    ht_t *sht   = *shtp;
    hi_t *hcm   = *hcmp;
    exp_t *etmp = bht->ev[0];
    memset(etmp, 0, (unsigned long)(bht->evl) * sizeof(exp_t));

    mat->rr = (hm_t **)malloc((unsigned long)bs->lml * 2 * sizeof(hm_t *));
    mat->nr = 0;
    mat->sz = 2 * bs->lml;

    /* add all non-redundant basis elements as matrix rows */
    for (i = 0; i < bs->lml; ++i) {
        mat->rr[mat->nr] = multiplied_poly_to_matrix_row(
                sht, bht, 0, etmp, bs->hm[bs->lmps[i]]);
        sht->hd[mat->rr[mat->nr][OFFSET]].idx  = 1;
        mat->nr++;
    }
    mat->nc = mat->nr; /* needed for correct counting in symbol */
    symbolic_preprocessing(mat, bs, st, sht, NULL, bht);
    /* no known pivots, we need mat->ncl = 0, so set all indices to 1 */
    for (i = 0; i < sht->eld; ++i) {
        sht->hd[i].idx = 1;
    }

    /* generate hash <-> column mapping */
    if (st->info_level > 1) {
        printf("reduce intermediate basis ");
        fflush(stdout);
    }
    convert_hashes_to_columns(&hcm, mat, st, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    /* do the linear algebra reduction, do NOT free basis data */
    interreduce_matrix_rows(mat, bs, st, 0);
    /* remap rows to basis elements (keeping their position in bs) */
    convert_sparse_matrix_rows_to_basis_elements(mat, bs, bht, sht, hcm, st);

    clear_matrix(mat);
    clean_hash_table(sht);

    /* printf("bs->lml %u | bs->ld %u | mat->np %u\n", bs->lml, bs->ld, mat->np); */
    /* we may have added some multiples of reduced basis polynomials
     * from the matrix, so we get rid of them. */
    i = 0;
    for (i = 0; i < bs->ld; ++i) {
        for (j = 0; j < mat->np; ++j) {
            if (bs->hm[i][OFFSET] == bs->hm[bs->ld+j][OFFSET]) {
                free(bs->hm[i]);
                bs->hm[i] = bs->hm[bs->ld+j];
                bs->hm[i][COEFFS] = i;
                switch(st->ff_bits) {
                    case 8:
                        free(bs->cf_8[i]);
                        bs->cf_8[i] = bs->cf_8[bs->ld+j];
                        break;
                    case 16:
                        free(bs->cf_16[i]);
                        bs->cf_16[i] = bs->cf_16[bs->ld+j];
                        break;
                    case 32:
                        free(bs->cf_32[i]);
                        bs->cf_32[i] = bs->cf_32[bs->ld+j];
                        break;
                }
            }
        }
    }
    for (i = bs->ld; i < bs->ld+mat->np; ++i) {
        bs->hm[i] = NULL;
        switch(st->ff_bits) {
            case 8:
                bs->cf_8[i] = NULL;
                break;
            case 16:
                bs->cf_16[i] = NULL;
                break;
            case 32:
                bs->cf_32[i] = NULL;
                break;
        }
    }
    *hcmp = hcm;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->reduce_gb_ctime = ct1 - ct0;
    st->reduce_gb_rtime = rt1 - rt0;
    if (st->info_level > 1) {
        printf("%13.2f sec\n", rt1-rt0);
    }

    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
}
#endif

static void reduce_basis(
        bs_t *bs,
        mat_t *mat,
        hi_t **hcmp,
        ht_t **bhtp,
        ht_t **shtp,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k;

    ht_t *bht   = *bhtp;
    ht_t *sht   = *shtp;
    hi_t *hcm   = *hcmp;
    exp_t *etmp = bht->ev[0];
    memset(etmp, 0, (unsigned long)(bht->evl) * sizeof(exp_t));

    mat->rr = (hm_t **)malloc((unsigned long)bs->lml * 2 * sizeof(hm_t *));
    mat->nr = 0;
    mat->sz = 2 * bs->lml;

    /* add all non-redundant basis elements as matrix rows */
    for (i = 0; i < bs->lml; ++i) {
        mat->rr[mat->nr] = multiplied_poly_to_matrix_row(
                sht, bht, 0, etmp, bs->hm[bs->lmps[i]]);
        sht->hd[mat->rr[mat->nr][OFFSET]].idx  = 1;
        mat->nr++;
    }
    mat->nc = mat->nr; /* needed for correct counting in symbol */
    symbolic_preprocessing(mat, bs, st, sht, NULL, bht);
    /* no known pivots, we need mat->ncl = 0, so set all indices to 1 */
    for (i = 0; i < sht->eld; ++i) {
        sht->hd[i].idx = 1;
    }

    /* free data from bht, we use sht later on */
    free_hash_table(&bht);

    /* generate hash <-> column mapping */
    if (st->info_level > 1) {
        printf("reduce final basis ");
        fflush(stdout);
    }
    convert_hashes_to_columns(&hcm, mat, st, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    /* do the linear algebra reduction and free basis data afterwards */
    interreduce_matrix_rows(mat, bs, st, 1);
    /* remap rows to basis elements (keeping their position in bs) */
    convert_sparse_matrix_rows_to_basis_elements_use_sht(1, mat, bs, sht, hcm, st);

    /* bht becomes sht, so we do not have to convert the hash entries */
    bht   = sht;
    *bhtp = bht;

    /* set sht = NULL, otherwise we might run in a double-free
     * of sht and bht at the end */
    sht   = NULL;
    *shtp = sht;

    bs->ld  = mat->np;

    /* clean_hash_table(sht); */
    clear_matrix(mat);

    /* we may have added some multiples of reduced basis polynomials
     * from the matrix, so we get rid of them. */
    k = 0;
    i = 0;
start:
    for (; i < bs->ld; ++i) {
        for (j = 0; j < k; ++j) {
            if (check_monomial_division(
                        bs->hm[bs->ld-1-i][OFFSET],
                        bs->hm[bs->lmps[j]][OFFSET], bht)) {
                ++i;
                goto start;
            }
        }
        bs->lmps[k] = bs->ld-1-i;
        bs->lm[k++] = bht->hd[bs->hm[bs->ld-1-i][OFFSET]].sdm;
    }
    bs->lml = k;

    *hcmp = hcm;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->reduce_gb_ctime = ct1 - ct0;
    st->reduce_gb_rtime = rt1 - rt0;
    if (st->info_level > 1) {
        printf("%13.2f sec\n", rt1-rt0);
    }

    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
}

int core_f4(
        bs_t **bsp,
        ht_t **bhtp,
        stat_t **stp
        )
{
    bs_t *bs    = *bsp;
    ht_t *bht   = *bhtp;
    stat_t *st  = *stp;

    /* timings for one round */
    double rrt0, rrt1;

    /* initialize update hash table, symbolic hash table */
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t *ps = initialize_pairset();

    int32_t round, i, j;

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis_f4(ps, bs, bht, st, st->ngens, 1);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    for (round = 1; ps->ld > 0; ++round) {
      if (round % st->reset_ht == 0) {
        reset_hash_table(bht, bs, ps, st);
        st->num_rht++;
      }
      rrt0  = realtime();
      st->max_bht_size  = st->max_bht_size > bht->esz ?
        st->max_bht_size : bht->esz;
      st->current_rd  = round;

      /* preprocess data for next reduction round */
      select_spairs_by_minimal_degree(mat, bs, ps, st, sht, bht, NULL);
      symbolic_preprocessing(mat, bs, st, sht, NULL, bht);
      convert_hashes_to_columns(&hcm, mat, st, sht);
      sort_matrix_rows_decreasing(mat->rr, mat->nru);
      sort_matrix_rows_increasing(mat->tr, mat->nrl);
      /* print pbm files of the matrices */
      if (st->gen_pbm_file != 0) {
        write_pbm_file(mat, st);
      }
      /* linear algebra, depending on choice, see set_function_pointers() */
      linear_algebra(mat, bs, st);
      /* columns indices are mapped back to exponent hashes */
      if (mat->np > 0) {
        convert_sparse_matrix_rows_to_basis_elements(
            -1, mat, bs, bht, sht, hcm, st);
      }
      clean_hash_table(sht);
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      /* check redundancy only if input is not homogeneous */
      update_basis_f4(ps, bs, bht, st, mat->np, 1-st->homogeneous);

      /* if we found a constant we are done, so remove all remaining pairs */
      if (bs->constant  == 1) {
          ps->ld  = 0;
      }
      rrt1 = realtime();
      if (st->info_level > 1) {
        printf("%13.2f sec\n", rrt1-rrt0);
      }
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    /* remove possible redudant elements */
    final_remove_redundant_elements(bs, bht);

    /* At the moment we do not directly remove the eliminated polynomials from
     * the resulting basis. */
#if 0
    if (st->nev > 0) {
        j = 0;
        for (i = 0; i < bs->lml; ++i) {
            if (bht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
                bs->lm[j]   = bs->lm[i];
                bs->lmps[j] = bs->lmps[i];
                ++j;
            }
        }
        bs->lml = j;
    }
#endif

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        /* note: bht will become sht, and sht will become NULL,
         * thus we need pointers */
        reduce_basis(bs, mat, &hcm, &bht, &sht, st);
    }

    *bsp  = bs;
    *bhtp = bht;
    *stp  = st;

    /* free and clean up */
    free(hcm);
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }

    return 1;
}

int64_t export_results_from_f4(
    /* return values */
    int32_t *bld,   /* basis load */
    int32_t **blen, /* length of each poly in basis */
    int32_t **bexp, /* basis exponent vectors */
    void **bcf,     /* coefficients of basis elements */
    void *(*mallocp) (size_t),
    bs_t **bsp,
    ht_t **bhtp,
    stat_t **stp
    )
{

    bs_t *bs    = *bsp;
    ht_t *bht   = *bhtp;
    stat_t *st  = *stp;

    st->nterms_basis  = export_julia_data(
        bld, blen, bexp, bcf, mallocp, bs, bht, st->fc);
    st->size_basis    = *bld;

    return st->nterms_basis;
}

/* we get from julia the generators as three arrays:
 * 0.  a pointer to an int32_t array for returning the basis to julia
 * 1.  an array of the lengths of each generator
 * 2.  an array of all coefficients of all generators in the order:
 *     first all coefficients of generator 1, then all of generator 2, ...
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *
 *  RETURNs the length of the jl_basis array */
int64_t f4_julia(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* data structures for basis, hash table and statistics */
    bs_t *bs    = NULL;
    ht_t *bht   = NULL;
    stat_t *st  = NULL;

    int success = 0;

    const int32_t use_signatures    =   0;
    success = initialize_gba_input_data(&bs, &bht, &st,
            lens, exps, cfs, field_char, mon_order, elim_block_len,
            nr_vars, nr_gens, 0 /* # normal forms */, ht_size,
            nr_threads, max_nr_pairs, reset_ht, la_option, use_signatures,
            reduce_gb, pbm_file, info_level);

    /* all input generators are invalid */
    if (success == -1) {
        return_zero(bld, blen, bexp, bcf, nr_vars, field_char, mallocp);
        return 1;
    }
    if (success == 0) {
        printf("Bad input data, stopped computation.\n");
        exit(1);
    }

    success = core_f4(&bs, &bht, &st);

    if (!success) {
        printf("Problem with F4, stopped computation.\n");
        exit(1);
    }

    int64_t nterms  = export_results_from_f4(bld, blen, bexp,
            bcf, mallocp, &bs, &bht, &st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 1) {
      print_final_statistics(stderr, st);
    }

    /* free and clean up */
    free_shared_hash_data(bht);
    if (bht != NULL) {
        free_hash_table(&bht);
    }

    if (bs != NULL) {
        free_basis(&bs);
    }

    free(st);
    st    = NULL;

    return nterms;
}
