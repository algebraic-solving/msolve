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


#include "modular.h"

trace_t *initialize_trace(
        void
        )
{
    trace_t *tr = (trace_t *)calloc(1, sizeof(trace_t));

    tr->sz  = 8;
    tr->ld  = 0;
    tr->td  = calloc((unsigned long)tr->sz, sizeof(td_t));
    /* rounds stuff for f4sat */
    tr->rsz = 8;
    tr->rld = 0;
    tr->rd  = calloc((unsigned long)tr->rsz, sizeof(len_t));
    tr->deg = calloc((unsigned long)tr->rsz, sizeof(deg_t));

    return tr;
}

void free_trace(
        trace_t **trp
        )
{
    trace_t *tr = *trp;
    len_t i, j;
    for (i = 0; i < tr->ld; ++i) {
        free(tr->td[i].tri);
        free(tr->td[i].rri);
        for (j = 0; j < tr->td[i].tld/2; ++j) {
            free(tr->td[i].rba[j]);
        }
        free(tr->td[i].rba);
        free(tr->td[i].lms);
    }
    free(tr->lm);
    free(tr->lmps);
    free(tr->td);
    free(tr->rd);
    free(tr->deg);
    free(tr);
    tr    = NULL;
    *trp  = tr;
}

void free_lucky_primes(
        primes_t **lpp
        )
{
    primes_t *lp  = *lpp;
    free(lp->p);
    free(lp);
    lp    = NULL;
    *lpp  = lp;
}

/* keeps bht since we know most of the basis already,
 * unlike the reduce_basis() function in f4.c */
void reduce_basis_no_hash_table_switching(
        bs_t *bs,
        mat_t *mat,
        hi_t **hcmp,
        ht_t *bht,
        ht_t *sht,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k;

    hi_t *hcm   = *hcmp;
    exp_t *etmp = bht->ev[0];
    memset(etmp, 0, (unsigned long)(bht->evl) * sizeof(exp_t));

    mat->rr = (hm_t **)malloc((unsigned long)bs->lml * 2 * sizeof(hm_t *));
    mat->nr = mat->nc = mat->ncl  = mat->ncr  = 0;
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
        printf("reduce basis       ");
        fflush(stdout);
    }
    convert_hashes_to_columns(&hcm, mat, st, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    /* do the linear algebra reduction */
    interreduce_matrix_rows(mat, bs, st);
    /* remap rows to basis elements (keeping their position in bs) */
    convert_sparse_matrix_rows_to_basis_elements(
        mat, bs, bht, sht, hcm, st);

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
        bs->lmps[k++] = bs->ld-1-i;
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

/* applies tracer, if prime is unlucky NULL is returned */
bs_t *f4_trace_application_phase(
        const trace_t * const trace,  /* trace of the F4 Algorithm */
        const ht_t * const tht,       /* trace hash table for multipliers */
        const bs_t * const ggb,       /* global basis */
        ht_t *lbht,                   /* local basis hash table, not shared */
        stat_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t ret, round, i;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));

    /* set routines corresponding to prime size */
    reset_trace_function_pointers(fc);

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    ht_t *bht   = lbht;

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash table */
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* reset bs->ld for first update process */
    bs->ld  = st->ngens;

    if(st->info_level>1){
      printf("Application phase with prime p = %d, overall there are %u rounds\n",
             fc, trace->ld);
    }
    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\nround   deg          mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    for (round = 0; round < trace->ld; ++round) {
      rrt0  = realtime();
      st->max_bht_size  = st->max_bht_size > bht->esz ?
        st->max_bht_size : bht->esz;
      st->current_rd  = round;

      /* generate matrix out of tracer data, rows are then already
       * sorted correspondingly */
      generate_matrix_from_trace(mat, trace, round, bs, st, sht, bht, tht);
        if (st->info_level > 1) {
            printf("%5d", round+1);
            printf("%6u ", sht->ev[mat->tr[0][OFFSET]][DEG]);
            fflush(stdout);
        }
      convert_hashes_to_columns(&hcm, mat, st, sht);
      /* linear algebra, depending on choice, see set_function_pointers() */
      ret = application_linear_algebra(mat, bs, st);
      if (ret != 0) {
          goto stop;
      }

      /* columns indices are mapped back to exponent hashes */
      if (mat->np > 0) {
          if (mat->np != trace->td[round].nlm) {
              fprintf(stderr, "Wrong number of new elements when applying tracer.");
              ret = 1;
              goto stop;
          }
          convert_sparse_matrix_rows_to_basis_elements(
                  mat, bs, bht, sht, hcm, st);
          for (i = 0; i < mat->np; ++i) {
              if (bs->hm[bs->ld+i][OFFSET] != trace->td[round].lms[i]) {
                  fprintf(stderr, "Wrong leading term for new element %u/%u.",
                          i, mat->np);
                  ret = 1;
                  goto stop;
              }
          }
      }
      bs->ld  +=  mat->np;
      clean_hash_table(sht);
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      rrt1 = realtime();
      if (st->info_level > 1) {
        printf("%13.2f sec\n", rrt1-rrt0);
      }
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }

    /* apply non-redundant basis data from trace to basis
     * before interreduction */
    bs->lml  = trace->lml;
    free(bs->lmps);
    bs->lmps = (bl_t *)calloc((unsigned long)bs->lml,
            sizeof(bl_t));
    memcpy(bs->lmps, trace->lmps,
            (unsigned long)bs->lml * sizeof(bl_t));
    free(bs->lm);
    bs->lm   = (sdm_t *)calloc((unsigned long)bs->lml,
            sizeof(sdm_t));
    memcpy(bs->lm, trace->lm,
            (unsigned long)bs->lml * sizeof(sdm_t));

    /* reduce final basis */
    /* note: bht will become sht, and sht will become NULL,
     * thus we need pointers */
    reduce_basis_no_hash_table_switching(
            bs, mat, &hcm, bht, sht, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    /* get basis meta data */
    st->size_basis  = bs->lml;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    if (st->info_level > 0) {
      print_final_statistics(stderr, st);
    }

stop:
    /* free and clean up */
    free(hcm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    gst->application_nr_add   = st->application_nr_add;
    gst->application_nr_mult  = st->application_nr_mult;
    gst->application_nr_red   = st->application_nr_red;
    free(st);

    if (ret != 0) {
        free_basis(&bs);
    }

    return bs;
}

bs_t *f4sat_trace_application_test_phase(
        const trace_t * const trace,  /* trace of the F4 Algorithm */
        const ht_t * const tht,       /* trace hash table for multipliers */
        const bs_t * const ggb,       /* global basis */
        const bs_t * const gsat,      /* global saturation element */
        ht_t *lbht,                   /* local basis hash table, not shared */
        stat_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    ps_t * ps   = initialize_pairset();
    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;
    int32_t round, ctr, i, j;
    ctr = 0;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* hashes-to-columns maps for multipliers in saturation step */
    hi_t *hcmm  = (hi_t *)malloc(sizeof(hi_t));

    /* set routines corresponding to prime size */
    reset_trace_function_pointers(fc);

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    bs_t *sat   = copy_basis_mod_p(gsat, st);
    ht_t *bht   = lbht;

    /* initialize multiplier of first element in sat to be the hash of
     * the all-zeroes exponent vector. */
    memset(bht->ev[0], 0, (unsigned long)(bht->nv+1) * sizeof(exp_t));
    sat->hm[0][MULT]  = insert_in_hash_table(bht->ev[0], bht);
    sat->ld = 1;
    len_t sat_deg = 0;

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash table */
    ht_t *sht = initialize_secondary_hash_table(bht, st);
    ht_t *uht = initialize_secondary_hash_table(bht, st);

    /* elements of kernel in saturation step, to be added to basis bs */
    bs_t *kernel  = initialize_basis(10);

    /* reset bs->ld for first update process */
    bs->ld  = 0;
    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    if(st->info_level>1){
        printf("Application phase with prime p = %d, overall there are %u rounds\n",
                fc, trace->ld);
    }
    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    round = 1;
    for (; ps->ld > 0; ++round) {
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
        exact_sparse_linear_algebra_ff_32(mat, bs, st);
        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            convert_sparse_matrix_rows_to_basis_elements(
                    mat, bs, bht, sht, hcm, st);
        }
        /* all rows in mat are now polynomials in the basis,
         * so we do not need the rows anymore */
        clear_matrix(mat); // does not reset mat->np

        update_basis(ps, bs, bht, uht, st, mat->np, 1);

        /* if we found a constant we are done, so remove all remaining pairs */
        rrt1 = realtime();
        if (st->info_level > 1) {
            printf("%13.2f sec\n", rrt1-rrt0);
        }
        if (bs->constant  == 1) {
            printf("basis is constant\n");
            ps->ld  = 0;
            break;
        }
        clean_hash_table(sht);

        /* saturation step starts here */
        if (trace->rd[ctr]  ==  round-1) {
            sat_deg = trace->deg[ctr++];
            rrt0  = realtime();
            /* printf("sat->deg %u\n", sat_deg); */
            update_multipliers(&qb, &bht, &sht, sat, st, bs, sat_deg);
            /* check for monomial multiples of elements from saturation list */
            select_saturation(sat, mat, st, sht, bht);

            symbolic_preprocessing(mat, bs, st, sht, NULL, bht);

            /* It may happen that there is no reducer at all for the
             * saturation elements, then nothing has to be done. */
            if (mat->nru > 0) {
                if (st->info_level > 1) {
                    /* printf("kernel computation "); */
                    printf("%3u  compute kernel", sat_deg);
                }
                /* int ctr = 0;
                 * for (int ii = 1; ii < sat->ld; ++ii) {
                 *     for (int jj = 0; jj < ii; jj++) {
                 *         if (sat->hm[ii][MULT] == sat->hm[jj][MULT]) {
                 *             printf("MULT %d == %d\n", ii, jj);
                 *             ctr++;
                 *         }
                 *     }
                 * } */
                /* int ctr  = 0;
                 * for (int ii = 0; ii<sat->ld; ++ii) {
                 *     if (sht->hd[sat->hm[ii][OFFSET]].idx == 2) {
                 *         sat->hm[ctr]  = sat->hm[ii];
                 *     } else {
                 *         free(sat->hm[ii]);
                 *         sat->hm[ii] = NULL;
                 *     }
                 * }
                 * sat->ld = ctr; */
                convert_hashes_to_columns_sat(&hcm, mat, sat, st, sht);
                convert_multipliers_to_columns(&hcmm, sat, st, bht);
                sort_matrix_rows_decreasing(mat->rr, mat->nru);

                compute_kernel_sat_ff_32(sat, mat, kernel, bs, st);

                if (kernel->ld > 0) {
                    if (st->info_level > 1) {
                        printf("\n                                               ");
                    }
                    clear_matrix(mat);
                    /* interreduce kernel */
                    copy_kernel_to_matrix(mat, kernel, sat->ld);
                    /* linear algebra, depending on choice, see set_function_pointers() */
                    exact_sparse_linear_algebra_ff_32(mat, kernel, st);
                    /* columns indices are mapped back to exponent hashes */
                    if (mat->np > 0) {
                        convert_sparse_matrix_rows_to_basis_elements_use_sht(
                                mat, bs, hcmm, st);
                    }
                    st->nr_kernel_elts  +=  kernel->ld;
                    free_kernel_coefficients(kernel);
                    update_basis(ps, bs, bht, uht, st, mat->np, 1);
                    kernel->ld  = 0;
                    if (st->info_level > 1) {
                        printf("   ");
                    }
                }
                /* columns indices are mapped back to exponent hashes */
                /* return_normal_forms_to_basis(
                 *         mat, tbr, bht, sht, hcm, st); */

                /* all rows in mat are now polynomials in the basis,
                 * so we do not need the rows anymore */
                convert_columns_to_hashes(sat, hcm, hcmm);
                for (i = 0; i < sat->ld; ++i) {
                    bht->hd[hcmm[i]].idx = 0;
                }
            }
            clear_matrix(mat);

            /* move hashes for sat entries from sht back to bht */
            for (i = 0; i < sat->ld; ++i) {
                if (sat->hm[i] != NULL) {
                    while (bht->esz - bht->eld < sat->hm[i][LENGTH]) {
                        enlarge_hash_table(bht);
                    }
                    for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
                        sat->hm[i][j] = insert_in_hash_table(
                                sht->ev[sat->hm[i][j]], bht);
                    }
                }
            }
            clean_hash_table(sht);

            rrt1 = realtime();
            if (st->info_level > 1) {
                printf("%10.2f sec\n", rrt1-rrt0);
            }
        }
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
                ----------------------------------------\n");
    }
    /* remove possible redudant elements */
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;

    /* apply non-redundant basis data from trace to basis
     * before interreduction */
    /* bs->lml  = trace->lml;
     * free(bs->lmps);
     * bs->lmps = (bl_t *)calloc((unsigned long)bs->lml,
     *         sizeof(bl_t));
     * memcpy(bs->lmps, trace->lmps,
     *         (unsigned long)bs->lml * sizeof(bl_t));
     * free(bs->lm);
     * bs->lm   = (sdm_t *)calloc((unsigned long)bs->lml,
     *         sizeof(sdm_t));
     * memcpy(bs->lm, trace->lm,
     *         (unsigned long)bs->lml * sizeof(sdm_t)); */

    /* reduce final basis */
    /* note: bht will become sht, and sht will become NULL,
     * thus we need pointers */
    reduce_basis_no_hash_table_switching(
            bs, mat, &hcm, bht, sht, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    /* get basis meta data */
    st->size_basis  = bs->lml;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    if (st->info_level > 0) {
        print_final_statistics(stderr, st);
    }

    /* free and clean up */
    free(hcm);
    free(hcmm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    free_basis_elements(sat);
    free_basis(&sat);
    free_basis(&kernel);
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    gst->application_nr_add   = st->application_nr_add;
    gst->application_nr_mult  = st->application_nr_mult;
    gst->application_nr_red   = st->application_nr_red;
    free(st);

    return bs;
}

bs_t *f4sat_trace_application_phase(
        const trace_t * const trace,  /* trace of the F4 Algorithm */
        const ht_t * const tht,       /* trace hash table for multipliers */
        const bs_t * const ggb,       /* global basis */
        const bs_t * const gsat,      /* global saturation element */
        ht_t *lbht,                   /* local basis hash table, not shared */
        stat_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;
    int32_t ret, round, ctr, i, j;
    ctr = 0;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* hashes-to-columns maps for multipliers in saturation step */
    hi_t *hcmm  = (hi_t *)malloc(sizeof(hi_t));

    /* set routines corresponding to prime size */
    reset_trace_function_pointers(fc);

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    bs_t *sat   = copy_basis_mod_p(gsat, st);
    ht_t *bht   = lbht;

    /* initialize multiplier of first element in sat to be the hash of
     * the all-zeroes exponent vector. */
    memset(bht->ev[0], 0, (unsigned long)(bht->nv+1) * sizeof(exp_t));
    sat->hm[0][MULT]  = insert_in_hash_table(bht->ev[0], bht);
    sat->ld = 1;
    len_t sat_deg = 0;

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash table */
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* elements of kernel in saturation step, to be added to basis bs */
    bs_t *kernel  = initialize_basis(10);

    bs->ld  = st->ngens;

    if(st->info_level>1){
        printf("Application phase with prime p = %d, overall there are %u rounds\n",
                fc, trace->ld);
    }
    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    round = 0;
    for (; round < trace->ld; ++round) {
        rrt0  = realtime();
        st->max_bht_size  = st->max_bht_size > bht->esz ?
            st->max_bht_size : bht->esz;
        st->current_rd  = round;

        /* generate matrix out of tracer data, rows are then already
         * sorted correspondingly */
        generate_matrix_from_trace(mat, trace, round, bs, st, sht, bht, tht);
        if (st->info_level > 1) {
            printf("%5d", round+1);
            printf("%6u ", sht->ev[mat->tr[0][OFFSET]][DEG]);
            fflush(stdout);
        }
        convert_hashes_to_columns(&hcm, mat, st, sht);
        /* linear algebra, depending on choice, see set_function_pointers() */
        ret = application_linear_algebra(mat, bs, st);
        if (ret != 0) {
            goto stop;
        }

        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            if (mat->np != trace->td[round].nlm) {
                fprintf(stderr, "Wrong number of new elements when applying tracer.");
                ret = 1;
                goto stop;
            }
            convert_sparse_matrix_rows_to_basis_elements(
                    mat, bs, bht, sht, hcm, st);
            for (i = 0; i < mat->np; ++i) {
                if (bs->hm[bs->ld+i][OFFSET] != trace->td[round].lms[i]) {
                    fprintf(stderr, "Wrong leading term for new element %u/%u.",
                            i, mat->np);
                    ret = 1;
                    goto stop;
                }
            }
            bs->ld  +=  mat->np;
            update_lm(bs, bht, st);
        }
        clean_hash_table(sht);
        /* all rows in mat are now polynomials in the basis,
         * so we do not need the rows anymore */
        clear_matrix(mat);

        rrt1 = realtime();
        if (st->info_level > 1) {
            printf("%13.2f sec\n", rrt1-rrt0);
        }
        /* saturation step starts here */
        if (trace->rd[ctr]  ==  round) {
            sat_deg = trace->deg[ctr++];
            /* check for new elements to be tested for adding saturation
             * information to the intermediate basis */
            rrt0  = realtime();
            update_multipliers(&qb, &bht, &sht, sat, st, bs, sat_deg);
            /* check for monomial multiples of elements from saturation list */
            select_saturation(sat, mat, st, sht, bht);

            symbolic_preprocessing(mat, bs, st, sht, NULL, bht);

            /* It may happen that there is no reducer at all for the
             * saturation elements, then nothing has to be done. */
            if (mat->nru > 0) {
                if (st->info_level > 1) {
                    /* printf("kernel computation "); */
                    printf("%3u  compute kernel", sat_deg);
                }
                /* int ctr = 0;
                 * for (int ii = 1; ii < sat->ld; ++ii) {
                 *     for (int jj = 0; jj < ii; jj++) {
                 *         if (sat->hm[ii][MULT] == sat->hm[jj][MULT]) {
                 *             printf("MULT %d == %d\n", ii, jj);
                 *             ctr++;
                 *         }
                 *     }
                 * } */
                /* int ctr  = 0;
                 * for (int ii = 0; ii<sat->ld; ++ii) {
                 *     if (sht->hd[sat->hm[ii][OFFSET]].idx == 2) {
                 *         sat->hm[ctr]  = sat->hm[ii];
                 *     } else {
                 *         free(sat->hm[ii]);
                 *         sat->hm[ii] = NULL;
                 *     }
                 * }
                 * sat->ld = ctr; */
                convert_hashes_to_columns_sat(&hcm, mat, sat, st, sht);
                convert_multipliers_to_columns(&hcmm, sat, st, bht);
                sort_matrix_rows_decreasing(mat->rr, mat->nru);

                compute_kernel_sat_ff_32(sat, mat, kernel, bs, st);

                if (kernel->ld == 0) {
                    fprintf(stderr, "Trivial kernel when applying tracer.");
                    ret = 1;
                    goto stop;
                }
                clear_matrix(mat);
                /* interreduce kernel */
                copy_kernel_to_matrix(mat, kernel, sat->ld);
                /* linear algebra, depending on choice, see set_function_pointers() */
                exact_sparse_linear_algebra_ff_32(mat, kernel, st);
                /* columns indices are mapped back to exponent hashes */
                if (mat->np > 0) {
                    convert_sparse_matrix_rows_to_basis_elements_use_sht(
                            mat, bs, hcmm, st);
                    bs->ld  +=  mat->np;
                    update_lm(bs, bht, st);
                }
                st->nr_kernel_elts  +=  kernel->ld;
                free_kernel_coefficients(kernel);
                kernel->ld  = 0;
                /* all rows in mat are now polynomials in the basis,
                 * so we do not need the rows anymore */
                convert_columns_to_hashes(sat, hcm, hcmm);
                for (i = 0; i < sat->ld; ++i) {
                    bht->hd[hcmm[i]].idx = 0;
                }
            }
            clear_matrix(mat);

            /* move hashes for sat entries from sht back to bht */
            for (i = 0; i < sat->ld; ++i) {
                if (sat->hm[i] != NULL) {
                    while (bht->esz - bht->eld < sat->hm[i][LENGTH]) {
                        enlarge_hash_table(bht);
                    }
                    for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
                        sat->hm[i][j] = insert_in_hash_table(
                                sht->ev[sat->hm[i][j]], bht);
                    }
                }
            }
            clean_hash_table(sht);

            rrt1 = realtime();
            if (st->info_level > 1) {
                printf("%10.2f sec\n", rrt1-rrt0);
            }
        }
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
                ----------------------------------------\n");
    }

    /* apply non-redundant basis data from trace to basis
     * before interreduction */
    bs->lml  = trace->lml;
    free(bs->lmps);
    bs->lmps = (bl_t *)calloc((unsigned long)bs->lml,
            sizeof(bl_t));
    memcpy(bs->lmps, trace->lmps,
            (unsigned long)bs->lml * sizeof(bl_t));
    free(bs->lm);
    bs->lm   = (sdm_t *)calloc((unsigned long)bs->lml,
            sizeof(sdm_t));
    memcpy(bs->lm, trace->lm,
            (unsigned long)bs->lml * sizeof(sdm_t));

    /* reduce final basis */
    /* note: bht will become sht, and sht will become NULL,
     * thus we need pointers */
    reduce_basis_no_hash_table_switching(
            bs, mat, &hcm, bht, sht, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    /* get basis meta data */
    st->size_basis  = bs->lml;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    if (st->info_level > 0) {
        print_final_statistics(stderr, st);
    }

stop:
    /* free and clean up */
    free(hcm);
    free(hcmm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    free_basis_elements(sat);
    free_basis(&sat);
    free_basis(&kernel);
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    gst->application_nr_add   = st->application_nr_add;
    gst->application_nr_mult  = st->application_nr_mult;
    gst->application_nr_red   = st->application_nr_red;
    free(st);

    if (ret != 0) {
        free_basis(&bs);
    }

    return bs;
}

/* In the learning phase we generate
 * 1. a groebner basis mod fc,
 * 2. a trace of the F4 run including a trace hash table tht
 * 3. a basis hash table (this should be nearly the same for
 *    all further runs of F4, so the subsequent runs just copy
 *    it and we will terms to the global basis hash table if
 *    needed on the fly) */
bs_t *f4_trace_learning_phase(
        trace_t *trace,           /* trace of the F4 Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        ht_t *gbht,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        stat_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t round, i, j;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    /* set routines corresponding to prime size */
    reset_trace_function_pointers(fc);

    ps_t * ps   = initialize_pairset();
    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    ht_t *bht   = gbht;

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash tables */
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
      printf("Learning phase with prime p = %d\n", fc);
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    for (round = 1; ps->ld > 0; ++round) {
      rrt0  = realtime();
      st->max_bht_size  = st->max_bht_size > bht->esz ?
        st->max_bht_size : bht->esz;
      st->current_rd  = round;

      /* preprocess data for next reduction round */
      select_spairs_by_minimal_degree(mat, bs, ps, st, sht, bht, tht);
      symbolic_preprocessing(mat, bs, st, sht, tht, bht);
      convert_hashes_to_columns(&hcm, mat, st, sht);
      sort_matrix_rows_decreasing(mat->rr, mat->nru);
      sort_matrix_rows_increasing(mat->tr, mat->nrl);
      /* linear algebra, depending on choice, see set_function_pointers() */
      trace_linear_algebra(trace, mat, bs, st);
      /* columns indices are mapped back to exponent hashes */
      if (mat->np > 0) {
        convert_sparse_matrix_rows_to_basis_elements(
            mat, bs, bht, sht, hcm, st);
      }
      clean_hash_table(sht);
      /* add lead monomials to trace, stores hashes in basis hash
       * table which is used in all upcoming F4 runs */
      if (mat->np > 0) {
          add_lms_to_trace(trace, bs, mat->np);
          trace->ld++;
      }
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      /* check redundancy only if input is not homogeneous */
      update_basis(ps, bs, bht, uht, st, mat->np, 1-st->homogeneous);

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
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;

    /* store information in trace */
    trace->lml  = bs->lml;
    trace->lmps = (bl_t *)calloc((unsigned long)trace->lml,
            sizeof(bl_t));
    memcpy(trace->lmps, bs->lmps,
            (unsigned long)trace->lml * sizeof(bl_t));
    trace->lm   = (sdm_t *)calloc((unsigned long)trace->lml,
            sizeof(sdm_t));
    memcpy(trace->lm, bs->lm,
            (unsigned long)trace->lml * sizeof(sdm_t));

    /* reduce final basis */
    /* note: bht will become sht, and sht will become NULL,
     * thus we need pointers */
    reduce_basis_no_hash_table_switching(bs, mat, &hcm, bht, sht, st);

    /* get basis meta data */
    st->size_basis  = bs->lml;
    len_t bsctr = 0;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
        if (bht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
            bsctr++;
        }
        /* printf("bs[%u] = ", i);
         * for (int ii = 0; ii < bht->evl; ++ii) {
         * printf("%d ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][ii]);
         * }
         * printf("\n"); */
    }
    if (st->nev > 0 && st->info_level > 0) {
        printf("eliminated basis -> %u\n", bsctr);
    }
    if (st->info_level > 0) {
      print_final_statistics(stderr, st);
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 0) {
      fflush(stdout);
      print_final_statistics(stderr, st);
      fflush(stderr);
    }

    /* free and clean up
     * note: we keep the basis hash table bht for all upcoming runs.
     *       this also means that we do not remove the shared data */
    free(hcm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    if (uht != NULL) {
        free_hash_table(&uht);
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    /* fix size of trace data */
    trace->td = realloc(trace->td, (unsigned long)trace->ld * sizeof(td_t));

    /* Note that we have to add also the numbers of the application counters since
     * for the interreduction of the matrices we use also during the learning phase
     * the reduction function which adds to the application counter (since we do not
     * use any tracing for the interreduction at all).
     * The application counters will be reset at the start of each application phase,
     * so we are not counting too many operations in the application runs. */
    gst->trace_nr_add   = st->trace_nr_add + st->application_nr_add;
    gst->trace_nr_mult  = st->trace_nr_mult + st->application_nr_mult;
    gst->trace_nr_red   = st->trace_nr_red + st->application_nr_red;

    free(st);

    return bs;
}

bs_t *f4sat_trace_learning_phase(
        trace_t *trace,           /* trace of the F4sat Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        const bs_t * const gsat,  /* global saturation element */
        ht_t **gbhtp,             /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        stat_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t round, i, j;
    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;

    len_t next_deg  = 0;
    /* global saturation data */
    len_t sat_test  = 0;
    deg_t sat_deg   = 0;
    /* int sat_done    = 0; */

    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* hashes-to-columns maps for multipliers in saturation step */
    hi_t *hcmm  = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    /* set routines corresponding to prime size */
    reset_trace_function_pointers(fc);

    ps_t * ps   = initialize_pairset();
    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    bs_t *sat   = copy_basis_mod_p(gsat, st);
    ht_t *bht   = *gbhtp;

    /* initialize multiplier of first element in sat to be the hash of
     * the all-zeroes exponent vector. */
    memset(bht->ev[0], 0, (unsigned long)(bht->nv+1) * sizeof(exp_t));
    sat->hm[0][MULT]  = insert_in_hash_table(bht->ev[0], bht);
    sat->ld = 1;

    next_deg  = 2*bht->ev[sat->hm[0][OFFSET]][DEG];

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash tables */
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* elements of kernel in saturation step, to be added to basis bs */
    bs_t *kernel  = initialize_basis(10);

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("Learning phase with prime p = %d\n", fc);
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    round = 1;
end_sat_step:
    for (; ps->ld > 0; ++round) {
        rrt0  = realtime();
        st->max_bht_size  = st->max_bht_size > bht->esz ?
            st->max_bht_size : bht->esz;
        st->current_rd  = round;

        /* preprocess data for next reduction round */
        select_spairs_by_minimal_degree(mat, bs, ps, st, sht, bht, tht);
        symbolic_preprocessing(mat, bs, st, sht, tht, bht);
        convert_hashes_to_columns(&hcm, mat, st, sht);
        sort_matrix_rows_decreasing(mat->rr, mat->nru);
        sort_matrix_rows_increasing(mat->tr, mat->nrl);
        /* linear algebra, depending on choice, see set_function_pointers() */
        trace_linear_algebra(trace, mat, bs, st);
        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            convert_sparse_matrix_rows_to_basis_elements(
                    mat, bs, bht, sht, hcm, st);
            sat_test++;
        }
        clean_hash_table(sht);
        /* add lead monomials to trace, stores hashes in basis hash
         * table which is used in all upcoming F4 runs */
        if (mat->np > 0) {
            add_lms_to_trace(trace, bs, mat->np);
            trace->ld++;
        }
        /* all rows in mat are now polynomials in the basis,
         * so we do not need the rows anymore */
        clear_matrix(mat);

        /* check redundancy only if input is not homogeneous */
        update_basis(ps, bs, bht, uht, st, mat->np, 1-st->homogeneous);

        /* if we found a constant we are done, so remove all remaining pairs */
        rrt1 = realtime();
        if (st->info_level > 1) {
            printf("%13.2f sec\n", rrt1-rrt0);
        }

        if (bs->constant  == 1) {
            printf("basis is constant\n");
            ps->ld  = 0;
            break;
        }
        clean_hash_table(sht);

        /* saturation step starts here */
        if (sat_test != 0 || ps->ld == 0) {
        /* if (bs->mltdeg - next_deg == 0 || ps->ld == 0) { */
        /* if (sat_done == 0 && (sat_test % 3 == 0 || ps->ld == 0)) { */
            if (st->nr_kernel_elts > 0 &&
                    ps->ld == 0 &&
                    is_zero_dimensional(bs, bht) &&
                    is_already_saturated(
                        bs, sat, mat, &hcm, &bht, &sht, &uht, st)) {
                /* sat_done  = 1; */
                goto end_sat_step;
            }
            /* check for new elements to be tested for adding saturation
             * information to the intermediate basis */
            if (ps->ld != 0) {
                sat_deg = 2*bs->mltdeg/3;
            } else {
                sat_deg = bs->mltdeg;
            }
            len_t bld = bs->ld;
            for (deg_t ii = next_deg; ii < sat_deg; ++ii) {
                rrt0  = realtime();
                /* printf("sat->deg %u\n", sat_deg); */
                update_multipliers(&qb, &bht, &sht, sat, st, bs, ii);
                /* check for monomial multiples of elements from saturation list */
                select_saturation(sat, mat, st, sht, bht);

                symbolic_preprocessing(mat, bs, st, sht, NULL, bht);

                /* It may happen that there is no reducer at all for the
                 * saturation elements, then nothing has to be done. */
                if (mat->nru > 0) {
                    if (st->info_level > 1) {
                        /* printf("kernel computation "); */
                        printf("%3u  compute kernel", sat_deg);
                    }
                    /* int ctr = 0;
                     * for (int ii = 1; ii < sat->ld; ++ii) {
                     *     for (int jj = 0; jj < ii; jj++) {
                     *         if (sat->hm[ii][MULT] == sat->hm[jj][MULT]) {
                     *             printf("MULT %d == %d\n", ii, jj);
                     *             ctr++;
                     *         }
                     *     }
                     * } */
                    /* int ctr  = 0;
                     * for (int ii = 0; ii<sat->ld; ++ii) {
                     *     if (sht->hd[sat->hm[ii][OFFSET]].idx == 2) {
                     *         sat->hm[ctr]  = sat->hm[ii];
                     *     } else {
                     *         free(sat->hm[ii]);
                     *         sat->hm[ii] = NULL;
                     *     }
                     * }
                     * sat->ld = ctr; */
                    convert_hashes_to_columns_sat(&hcm, mat, sat, st, sht);
                    convert_multipliers_to_columns(&hcmm, sat, st, bht);
                    sort_matrix_rows_decreasing(mat->rr, mat->nru);

                    compute_kernel_sat_ff_32(sat, mat, kernel, bs, st);

                    if (kernel->ld > 0) {
                        if (st->info_level > 1) {
                            printf("\n                                               ");
                        }
                        clear_matrix(mat);
                        /* interreduce kernel */
                        copy_kernel_to_matrix(mat, kernel, sat->ld);
                        /* linear algebra, depending on choice, see set_function_pointers() */
                        probabilistic_sparse_linear_algebra_ff_32(mat, kernel, st);
                        /* linear_algebra(mat, kernel, st); */
                        /* columns indices are mapped back to exponent hashes */
                        if (mat->np > 0) {
                            convert_sparse_matrix_rows_to_basis_elements_use_sht(
                                    mat, bs, hcmm, st);
                        }
                        /* track round in which kernel computation is not trivial */
                        if (trace->rld == trace->rsz) {
                            trace->rsz  *=  2;
                            trace->rd = realloc(
                                    trace->rd,
                                    (unsigned long)trace->rsz * sizeof(len_t));
                        }
                        /* in the application phase round counting
                         * starts at 0, so that's OK */
                        trace->rd[trace->rld]     = trace->ld-1;
                        trace->deg[trace->rld++]  = ii;

                        st->nr_kernel_elts  +=  kernel->ld;
                        sat_test  = 0;
                        free_kernel_coefficients(kernel);
                        update_basis(ps, bs, bht, uht, st, mat->np, 1);
                        kernel->ld  = 0;
                        if (st->info_level > 1) {
                            printf("   ");
                        }
                    }
                    /* columns indices are mapped back to exponent hashes */
                    /* return_normal_forms_to_basis(
                     *         mat, tbr, bht, sht, hcm, st); */

                    /* all rows in mat are now polynomials in the basis,
                     * so we do not need the rows anymore */
                    convert_columns_to_hashes(sat, hcm, hcmm);
                    for (i = 0; i < sat->ld; ++i) {
                        bht->hd[hcmm[i]].idx = 0;
                    }
                }
                clear_matrix(mat);

                /* move hashes for sat entries from sht back to bht */
                for (i = 0; i < sat->ld; ++i) {
                    if (sat->hm[i] != NULL) {
                        while (bht->esz - bht->eld < sat->hm[i][LENGTH]) {
                            enlarge_hash_table(bht);
                        }
                        for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
                            sat->hm[i][j] = insert_in_hash_table(
                                    sht->ev[sat->hm[i][j]], bht);
                        }
                    }
                }
                clean_hash_table(sht);

                rrt1 = realtime();
                if (st->info_level > 1) {
                    printf("%10.2f sec\n", rrt1-rrt0);
                }
                if (bld != bs->ld) {
                    next_deg  = ii;
                    goto end_sat_step;
                }
            }
            next_deg  = sat_deg;
        }
    }

    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    /* remove possible redudant elements */
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;

    /* store information in trace */
    trace->lml  = bs->lml;
    trace->lmps = (bl_t *)calloc((unsigned long)trace->lml,
            sizeof(bl_t));
    memcpy(trace->lmps, bs->lmps,
            (unsigned long)trace->lml * sizeof(bl_t));
    trace->lm   = (sdm_t *)calloc((unsigned long)trace->lml,
            sizeof(sdm_t));
    memcpy(trace->lm, bs->lm,
            (unsigned long)trace->lml * sizeof(sdm_t));

    /* reduce final basis */
    /* note: bht will become sht, and sht will become NULL,
     * thus we need pointers */
    reduce_basis_no_hash_table_switching(bs, mat, &hcm, bht, sht, st);

/*     printf("basis has  %u elements.\n", bs->lml);
 *
 *     for (i = 0; i < bs->lml; ++i) {
 *         for (j = 0; j < bht->nv; ++j) {
 *             printf("%u ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][j]);
 *         }
 *         printf("\n");
 *     } */
    /* get basis meta data */
    st->size_basis  = bs->lml;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 0) {
      fflush(stdout);
      print_final_statistics(stderr, st);
      fflush(stderr);
    }

    /* free and clean up
     * note: we keep the basis hash table bht for all upcoming runs.
     *       this also means that we do not remove the shared data */
    free(hcm);
    free(hcmm);
    free(qb);
    *gbhtp = bht;

    if (sht != NULL) {
        free_hash_table(&sht);
    }
    if (uht != NULL) {
        free_hash_table(&uht);
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }
    free_basis_elements(sat);
    free_basis(&sat);
    free_basis(&kernel);
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    /* fix size of trace data */
    trace->td = realloc(trace->td, (unsigned long)trace->ld * sizeof(td_t));

    /* Note that we have to add also the numbers of the application counters since
     * for the interreduction of the matrices we use also during the learning phase
     * the reduction function which adds to the application counter (since we do not
     * use any tracing for the interreduction at all).
     * The application counters will be reset at the start of each application phase,
     * so we are not counting too many operations in the application runs. */
    gst->trace_nr_add   = st->trace_nr_add + st->application_nr_add;
    gst->trace_nr_mult  = st->trace_nr_mult + st->application_nr_mult;
    gst->trace_nr_red   = st->trace_nr_red + st->application_nr_red;

    free(st);

    return bs;
}


int64_t f4_trace_julia(
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
        const uint32_t prime_start,
        const int32_t nr_primes,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    /* only for computations over the rationals */
    if (field_char != 0) {
        fprintf(stderr, "Tracer only for computations over Q. Call\n");
        fprintf(stderr, "standard F4 Algorithm for computations over\n");
        fprintf(stderr, "finite fields.\n");
        return 1;
    }

    len_t i;

    /* ps is only needed for storing meta data */
    ps_t *ps  = initialize_pairset();

    /* lucky primes */
    primes_t *lp  = (primes_t *)calloc(1, sizeof(primes_t));

    /* initialize stuff */
    stat_t *st  = initialize_statistics();

    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */
    if (check_and_set_meta_data_trace(st, lens, exps, cfs, field_char,
                mon_order, elim_block_len, nr_vars, nr_gens, ht_size, nr_threads,
                max_nr_pairs, reset_ht, la_option, reduce_gb, prime_start,
                nr_primes, pbm_file, info_level)) {
        return 0;
    }

    /*******************
    * initialize basis
    *******************/
    bs_t *bs_qq = initialize_basis(st->ngens);
    /* initialize basis hash table, update hash table, symbolic hash table */
    ht_t *bht = initialize_basis_hash_table(st);
    /* hash table to store the hashes of the multiples of
     * the basis elements stored in the trace */
    ht_t *tht = initialize_secondary_hash_table(bht, st);
    /* read in ideal, move coefficients to integers */
    import_julia_data(bs_qq, bht, st, lens, exps, cfs);

    if (st->info_level > 0) {
      print_initial_statistics(stderr, st);
    }

    /* for faster divisibility checks, needs to be done after we have
     * read some input data for applying heuristics */
    calculate_divmask(bht);

    /* sort initial elements, smallest lead term first */
    sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *),
            initial_input_cmp, bht);
    remove_content_of_initial_basis(bs_qq);

    /* generate lucky prime numbers */
    generate_lucky_primes(lp, bs_qq, st->prime_start, st->nprimes);

    /* generate array to store modular bases */
    bs_t **bs = (bs_t **)calloc((unsigned long)st->nprimes, sizeof(bs_t *));

    /* initialize tracer */
    trace_t *trace  = initialize_trace();

    /* learning phase */
    bs[0] = f4_trace_learning_phase(trace, tht, bs_qq, bht, st, lp->p[0]);

    /* tracing phase */
#pragma omp parallel for num_threads(st->nthrds) \
    private(i) schedule(dynamic)
    for (i = 1; i < st->nprimes; ++i) {
        bs[i] = f4_trace_application_phase(
                trace, tht, bs_qq, bht, st, lp->p[i]);
    }

    /* reconstruction phase */

    /* testing phase */

    /* free and clean up */
    free_trace(&trace);
    free_shared_hash_data(bht);
    free_hash_table(&bht);
    free_pairset(&ps);
    for (i = 0; i < st->nprimes; ++i) {
        free_basis(&(bs[i]));
    }
    free(bs);
    free_lucky_primes(&lp);
    free(st);

    return 0;
}

/* modular f4 call for usual multi-modular f4, no tracing */
bs_t *modular_f4(
        const bs_t * const ggb,       /* global basis */
        ht_t * gbht,                  /* global basis hash table, shared */
        stat_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t round, i, j;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));

    /* set routines corresponding to prime size */
    reset_function_pointers(fc, gst->laopt);

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t * ps   = initialize_pairset();

    /* copy global data as input */
    stat_t *st  = copy_statistics(gst, fc);
    bs_t *bs    = copy_basis_mod_p(ggb, st);
    ht_t *bht   = gbht;

    /* normalize the copied basis */
    normalize_initial_basis(bs, fc);

    /* initialize specialized hash table */
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    for (round = 1; ps->ld > 0; ++round) {
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
            mat, bs, bht, sht, hcm, st);
      }
      clean_hash_table(sht);
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      /* check redundancy only if input is not homogeneous */
      update_basis(ps, bs, bht, uht, st, mat->np, 1-st->homogeneous);

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
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        reduce_basis_no_hash_table_switching(
                bs, mat, &hcm, bht, sht, st);
        /* reduce_basis_(bs, mat, &hcm, &bht, &sht, st); */
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    /* get basis meta data */
    st->size_basis  = bs->lml;
    len_t bsctr = 0;
    for (i = 0; i < bs->lml; ++i) {
        st->nterms_basis +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
        if (bht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
            bsctr++;
        }
        /* printf("bs[%u] = ", i);
         * for (int ii = 0; ii < bht->evl; ++ii) {
         * printf("%d ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][ii]);
         * }
         * printf("\n"); */
    }
    if (st->nev > 0 && st->info_level > 0) {
        printf("eliminated basis -> %u\n", bsctr);
    }
    if (st->info_level > 0) {
      print_final_statistics(stderr, st);
    }

    /* free and clean up */
    free(hcm);
    if (sht != NULL) {
        free_hash_table(&sht);
    }
    if (uht != NULL) {
        free_hash_table(&uht);
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    free(st);

    return bs;
}
