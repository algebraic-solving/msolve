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


#include "f4sat.h"

static inline unsigned long sum(
        const len_t * const ind,
        const len_t len)
{

    len_t i;
    unsigned long sum = 0;

    for (i = 0; i < len; i++){
        sum += ind[i];
    }
    return sum;
}

static inline len_t generate_new_basis_elements(
        ht_t **htp,
        hm_t *cqb,
        const hm_t * const oqb,
        const len_t oqb_dim,
        const len_t * const ind,
        const bs_t * const bs
        )
{
    len_t i, j;

    ht_t *ht        = *htp;
    const len_t nv  = ht->nv;
    len_t ctr       = 0;
    exp_t *etmp     = calloc((unsigned long)nv, sizeof(exp_t));

    for (i = nv; i > 0; --i) {
        while (ht->esz - ht->eld < oqb_dim-ind[nv-i]) {
            enlarge_hash_table(ht);
        }
        for (j = ind[nv-i]; j < oqb_dim; ++j) {
            memcpy(etmp, ht->ev[oqb[j]], (unsigned long)nv*sizeof(exp_t));;
            etmp[i-1]++;
            cqb[ctr]  = check_lm_divisibility_and_insert_in_hash_table(
                            etmp, ht, bs
                            );
            if (cqb[ctr] != 0) {
                ctr++;
            }
        }
    }
    *htp  = ht;
    free(etmp);

    return ctr;
}

static inline void update_indices(
        len_t *ind,
        const hm_t * const qb,
        const len_t oqb_dim,
        const len_t cqb_dim,
        const ht_t * const ht
        )
{
    len_t i, j, k;

    const len_t nv  = ht->nv;
    const len_t dim = oqb_dim + cqb_dim;
    ind[0]  = oqb_dim;

    for (i = nv-1; i > 0; --i) {
        for (j = ind[nv-1-i]; j < dim; ++j) {
            if (ht->ev[qb[j]][i] == 0) {
                ind[nv-i] = j;
                break;
            }
            for (k = i; k > 0; --k) {
                ind[nv-k] = dim;
            }
        }
    }
}

static len_t quotient_basis(
        hm_t **qbp,
        ht_t **htp,
        const bs_t * const bs,
        const deg_t max_deg
        )
{
    if (bs->constant == 1) {
        return 0;
    }
    hm_t *qb  = *qbp;
    /* "qb" denotes the quotient basis storing the hash table entries.
     * "qb" is the basis itself from the last step, "nqb" are the new"
     * elements we generate in the current step. */
    hm_t *nqb = NULL;
    /* dimensions of qb and nqb */
    len_t qbd   = 1;
    len_t nqbd  = 1; // for entering first round in for loop
    deg_t deg   = 0;

    const len_t nv  = (*htp)->nv;

    /* clear ht-ev[0] */
    memset((*htp)->ev[0], 0, (unsigned long)nv * sizeof(exp_t));


    len_t *ind      = calloc((unsigned long)nv, sizeof(len_t));

    qb  = (hm_t *)calloc((unsigned long)1, sizeof(hm_t));
    qb[0] = check_lm_divisibility_and_insert_in_hash_table(
                            (*htp)->ev[0], (*htp), bs);

    while (nqbd > 0 && deg < max_deg) {
        nqb = (hm_t *)calloc(sum(ind, nv) + nv, sizeof(hm_t));
        nqbd = generate_new_basis_elements(htp, nqb, qb, qbd, ind, bs);
        qb = realloc(qb, (unsigned long)(qbd + nqbd) * sizeof(hm_t));
        /* printf("qb before adding\n");
         * for (len_t ii = 0; ii < qbd; ++ii) {
         *     printf("pos %u --> ", ii);
         *     for (len_t jj = 0; jj < nv; ++jj) {
         *         printf("%d ", (*htp)->ev[qb[ii]][jj]);
         *     }
         *     printf("\n");
         * } */
        /* printf("nqb before adding\n");
         * for (len_t ii = 0; ii < nqbd; ++ii) {
         *     printf("pos %u --> ", ii);
         *     for (len_t jj = 0; jj < nv; ++jj) {
         *         printf("%d ", (*htp)->ev[nqb[ii]][jj]);
         *     }
         *     printf("\n");
         * } */
        memcpy(qb+qbd, nqb, (unsigned long)nqbd * sizeof(hm_t));
        update_indices(ind, qb, qbd, nqbd, *htp);
        qbd   +=  nqbd;
        deg++;
        /* printf("qb after adding\n");
         * for (len_t ii = 0; ii < qbd; ++ii) {
         *     printf("pos %u --> ", ii);
         *     for (len_t jj = 0; jj < nv; ++jj) {
         *         printf("%d ", (*htp)->ev[qb[ii]][jj]);
         *     }
         *     printf("\n");
         * } */
    }
    free(nqb);
    free(ind);

    *qbp  = qb;
    return qbd;
}

static void update_multipliers(
        hm_t **qdp,
        ht_t **bhtp,
        ht_t **shtp,
        bs_t *sat,
        stat_t *st,
        const bs_t *bs,
        const deg_t max_deg
        )
{
    len_t i;

    len_t qdim  = quotient_basis(qdp, bhtp, bs, max_deg);
    ht_t *bht   = *bhtp;
    ht_t *sht   = *shtp;

    bl_t ctr  = 0;

    check_enlarge_basis(sat, qdim);

    hm_t *qb  = *qdp;

    for (i = 0; i < sat->ld; ++i) {
        while (i < sat->ld && qb[ctr] != sat->hm[i][MULT]) {
            free(sat->hm[i]);
            sat->hm[i]    = NULL;
            free(sat->cf_32[i]);
            sat->cf_32[i] = NULL;
            i++;
        }
        sat->hm[ctr]          = sat->hm[i];
        sat->cf_32[ctr]       = sat->cf_32[i];
        sat->hm[ctr][COEFFS]  = ctr;
        ctr++;
    }

    /* TODO: Here we could apply the simplify idea from F4:
    * Any of these new elements is a multiple of an element
    * of lower degree already handled. Thus we can use this
    * "divisor" instead of the initial saturation element sat[0]. */
    const hm_t * const b      = sat->hm[0];
    const cf32_t * const cfs  = sat->cf_32[0];
    /* new saturation elements from new quotient monomials */
    for (i = ctr; i < qdim; ++i) {
        const hm_t m      = qb[i];
        const hi_t h      = bht->hd[m].val;
        const deg_t d     = bht->hd[m].deg;
        sat->hm[i]        = multiplied_poly_to_matrix_row(
                sht, bht, h, d, bht->ev[m], b);
       sat->hm[ctr][MULT] = qb[i];
       sat->cf_32[ctr]    = (cf32_t *)malloc(
               (unsigned long)b[LENGTH] * sizeof(cf32_t));
       memcpy(sat->cf_32[ctr], cfs,
               (unsigned long)b[LENGTH] * sizeof(cf32_t));
       sat->hm[ctr][COEFFS] = ctr;
       ctr++;
    }
    sat->ld = ctr;
    st->new_multipliers = sat->ld - sat->lo;
    /* for (i = 0; i < sat->ld; ++i) {
     *     printf("%3u -> ", i);
     *     for (len_t j = 0; j < ht->nv; ++j) {
     *         printf("%d ", ht->ev[sat->hm[i][OFFSET]][j]);
     *     }
     *     printf("\n");
     * } */

    if (sat->mltdeg < max_deg) {
        sat->mltdeg = max_deg;
    }
    *bhtp = bht;
    *shtp = sht;
}

int core_f4sat(
        bs_t **bsp,
        bs_t **satp,
        ht_t **bhtp,
        stat_t **stp
        )
{
    bs_t *bs    = *bsp;
    ht_t *bht   = *bhtp;
    stat_t *st  = *stp;

    /* elements to saturate input ideal with */
    bs_t *sat   = *satp;
    /* initialize multiplier of first element in sat to be the hash of
     * the all-zeroes exponent vector. */
    memset(bht->ev[0], 0, (unsigned long)bht->nv * sizeof(exp_t));
    sat->hm[0][MULT]  = insert_in_hash_table(bht->ev[0], bht);
    sat->ld = 1;

    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;
    /* elements tracking reduction steps of tbr, being
     * added to bs if we find a linear dependency on
     * elements in tbr. They represent the multipliers
     * of the elements we saturate with. */
    /* timings for one round */
    double rrt0, rrt1;

    /* initialize update hash table, symbolic hash table */
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm   = (hi_t *)malloc(sizeof(hi_t));
    /* hashes-to-columns maps for multipliers in saturation step */
    hi_t *hcmm  = (hi_t *)malloc(sizeof(hi_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t *ps = initialize_pairset();

    int32_t round, i, j;

    /* elements of kernel in saturation step, to be added to basis bs */
    bs_t *kernel  = initialize_basis(10);
    printf("kernel %p | ld %u\n", kernel, kernel->ld);

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
        if (mat->np > 0) {
            /* check for new elements to be tested for adding saturation
             * information to the intermediate basis */
            rrt0  = realtime();
            update_multipliers(&qb, &bht, &sht, sat, st, bs, bs->mltdeg);
            /* check for monomial multiples of elements from saturation list */
            select_saturation(sat, mat, st, sht, bht);

            symbolic_preprocessing(mat, bs, st, sht, NULL, bht);

            /* It may happen that there is no reducer at all for the
             * saturation elements, then nothing has to be done. */
            if (mat->nru > 0) {
                if (st->info_level > 1) {
                    printf("kernel computation ");
                }
                convert_hashes_to_columns_sat(&hcm, mat, sat, st, sht);
                convert_multipliers_to_columns(&hcmm, sat, st, bht);
                sort_matrix_rows_decreasing(mat->rr, mat->nru);

                compute_kernel_sat_ff_32(sat, mat, kernel, bs, st);

                if (kernel->ld > 0) {
                    add_kernel_elements_to_basis(
                            bs, kernel, bht, hcmm, st);
                    update_basis(ps, bs, bht, uht, st, kernel->ld, 1);

                }
                for (i = 0; i < sat->ld; ++i) {
                    bht->hd[hcmm[i]].idx = 0;
                }
                /* columns indices are mapped back to exponent hashes */
                /* return_normal_forms_to_basis(
                 *         mat, tbr, bht, sht, hcm, st); */

                /* all rows in mat are now polynomials in the basis,
                 * so we do not need the rows anymore */
                convert_columns_to_hashes(sat, hcm);
            }
            clear_matrix(mat);
            clean_hash_table(sht);

            /* move hashes for sat entries from sht back to bht */
            for (i = 0; i < sat->ld; ++i) {
                for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
                    sat->hm[i][j] = insert_in_hash_table(
                            sht->ev[sat->hm[i][j]], bht);
                }
            }

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

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        /* note: bht will become sht, and sht will become NULL,
         * thus we need pointers */
        reduce_basis(bs, mat, &hcm, &bht, &sht, st);
    }
    printf("basis has  %u elements.\n", bs->lml);

    for (i = 0; i < bs->lml; ++i) {
        for (j = 0; j < bht->nv; ++j) {
            printf("%u ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][j]);
        } 
        printf("\n");
    }

    *bsp  = bs;
    *satp = sat;
    *bhtp = bht;
    *stp  = st;

    /* free and clean up */
    free(hcm);
    free(hcmm);
    free(qb);

    free_basis_elements(sat);
    free_basis(&sat);
    free_basis(&kernel);
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
