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


#include "sba.h"

static inline void initialize_signatures_schreyer(
        bs_t *bs
        )
{
    for (len_t i = 0; i < bs->ld; ++i) {
        bs->si[i]   =   i;
        bs->sm[i]   =   bs->hm[i][OFFSET];;
    }
}

static inline void initialize_signatures_not_schreyer(
        bs_t *bs
        )
{
    for (len_t i = 0; i < bs->ld; ++i) {
        bs->si[i]   =   i;
        bs->sm[i]   =   0;
    }
}

int core_sba_schreyer(
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

    int32_t round, i, j;

    /* initialize update hash table, symbolic hash table */
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t *ps = initialize_pairset();

    initialize_signatures_schreyer(bs);

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis_sba_schreyer(ps, bs, bht, uht, st, st->ngens, 1);

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
            mat, bs, bht, sht, hcm, st);
      }
      clean_hash_table(sht);
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      /* check redundancy only if input is not homogeneous */
      update_basis_sba_schreyer(ps, bs, bht, uht, st, mat->np, 1-st->homogeneous);

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

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        /* note: bht will become sht, and sht will become NULL,
         * thus we need pointers */
        reduce_basis(bs, mat, &hcm, &bht, &sht, st);
    }
    /* for (i = 0; i < bs->lml; ++i) { */
    /*     for (j = 0; j < bht->evl; ++j) { */
    /*         printf("%d ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][j]); */
    /*     } */
    /*     printf("\n"); */
    /* } */

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
    if (uht != NULL) {
        free_hash_table(&uht);
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }

    return 1;
}

