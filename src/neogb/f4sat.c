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

static inline void free_kernel_coefficients(
        bs_t *kernel
        )
{
    len_t i;

    for (i = 0; i < kernel->ld; ++i) {
        free(kernel->cf_32[i]);
    }
}

static inline int is_pure_power(
        const exp_t * const ev,
        const ht_t * const ht
        )
{
    len_t i;
    len_t ctr = 0;

    const len_t ebl = ht->ebl;
    const len_t evl = ht->evl;
    const len_t nv  = ht->nv;
    
    for (i = 1; i < ebl; ++i) {
        if (ev[i] == 0) {
            ctr++;
        }
    }
    for (i = ebl+1; i < evl; ++i) {
        if (ev[i] == 0) {
            ctr++;
        }
    }

    if (ctr == nv-1) {
        return 1;
    }
    return 0;
}

static inline int is_zero_dimensional(
        const bs_t * const bs,
        const ht_t * const ht
        )
{
    len_t i;

    len_t ctr = 0;
    
    const len_t nv  = ht->nv;
    const len_t lml = bs->lml;

    for (i = 0; i < lml; ++i) {
        if (is_pure_power(ht->ev[bs->hm[bs->lmps[i]][OFFSET]], ht)) {
            ctr ++;
        }
    }
    
    if (ctr == nv) {
        return 1;
    }
    return 0;
}

static int is_already_saturated(
        bs_t *bs,
        const bs_t * const sat,
        mat_t *mat,
        ht_t **bhtp,
        ht_t **shtp,
        md_t *st
        )
{
    printf("testing if system is already saturated: ");
    double rrt0, rrt1;
    rrt0  = realtime();

    len_t i;

    /* hi_t *hcm = *hcmp; */
    ht_t *bht = *bhtp;
    ht_t *sht = *shtp;

    /* add phi to basis and generate pairs with phi */
    check_enlarge_basis(bs, 1, st);

    int is_constant = 0;

    /* copy old data from bs to restore after test */
    const bl_t bld  = bs->ld;
    const bl_t blo  = bs->lo;
    const bl_t bsc  = bs->constant;
    const len_t lml = bs->lml;
    sdm_t *lm       = (sdm_t *)malloc((unsigned long)lml * sizeof(sdm_t));
    memcpy(lm, bs->lm, (unsigned long)lml * sizeof(sdm_t));
    bl_t *lmps      = (bl_t *)malloc((unsigned long)lml * sizeof(bl_t));
    memcpy(lmps, bs->lmps, (unsigned long)lml * sizeof(bl_t));
    int8_t *red     = (int8_t *)malloc((unsigned long)bs->sz * sizeof(int8_t));
    memcpy(red, bs->red, (unsigned long)bs->sz * sizeof(int8_t));
    
    ps_t *ps = initialize_pairset();

    cf32_t *cf  = (cf32_t *)malloc(
            (unsigned long)sat->hm[0][LENGTH] * sizeof(cf32_t));
    memcpy(cf, sat->cf_32[sat->hm[0][COEFFS]],
            (unsigned long)sat->hm[0][LENGTH] * sizeof(cf32_t));
    hm_t *hm  = (hm_t *)malloc(
            (unsigned long)(sat->hm[0][LENGTH]+OFFSET) * sizeof(hm_t));
    memcpy(hm, sat->hm[0],
            (unsigned long)(sat->hm[0][LENGTH]+OFFSET) * sizeof(hm_t));

    bs->cf_32[bs->ld] = cf;
    hm[COEFFS]        = bs->ld;
    bs->hm[bs->ld]    = hm;

    update_basis_f4(ps, bs, bht, st, 1);
    
    /* suppress infolevel printing in the test step */
    int32_t infolevel = st->info_level;
    st->info_level  = 0;
    while (ps->ld > 0) {

        select_spairs_by_minimal_degree(mat, bs, st);
        /* select_all_spairs(mat, bs, ps, st, sht, bht, NULL); */
        symbolic_preprocessing(mat, bs, st);
        convert_hashes_to_columns(mat, st, sht);
        sort_matrix_rows_decreasing(mat->rr, mat->nru);
        sort_matrix_rows_increasing(mat->tr, mat->nrl);
        /* linear algebra, depending on choice, see set_function_pointers() */
        probabilistic_sparse_linear_algebra_ff_32(mat, bs, bs, st);

        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            convert_sparse_matrix_rows_to_basis_elements(
                    -1, mat, bs, bht, sht, st);
        }
        /* all rows in mat are now polynomials in the basis,
         * so we do not need the rows anymore */
        clear_matrix(mat); // does not reset mat->np
        clean_hash_table(sht);

        update_basis_f4(ps, bs, bht, st, mat->np);

        /* if we found a constant we are done, so remove all remaining pairs */
        if (bs->constant  == 1) {
            ps->ld  = 0;
            break;
        }

    }

    is_constant = bs->constant;

    /* reset basis data to state before starting test */
    for (i = bld; i < bs->ld; ++i) {
        free(bs->cf_32[bs->hm[i][COEFFS]]);
        bs->cf_32[bs->hm[i][COEFFS]]  = NULL;
        free(bs->hm[i]);
        bs->hm[i] = NULL;
    }
    if (ps != NULL) {
        free_pairset(&ps);
    }
    bs->ld  = bld;
    bs->lo  = blo;

    bs->constant  = bsc;

    /* reset infolevel */
    st->info_level  = infolevel;

    free(bs->lm);
    bs->lm    = lm;
    free(bs->lmps);
    bs->lmps  = lmps;
    bs->lml   = lml;
    free(bs->red);
    bs->red   = red;

    /* *hcmp = hcm; */
    *bhtp = bht;
    *shtp = sht;

    if (is_constant == 1) {
        printf("yes.");
    } else {
        printf("no.");
    }

    rrt1 = realtime();
    if (st->info_level > 1) {
        printf("%40.2f sec\n", rrt1-rrt0);
    }

    return is_constant;
}

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
    exp_t *etmp     = calloc((unsigned long)(nv+1), sizeof(exp_t));

    for (i = nv; i > 0; --i) {
        while (ht->esz - ht->eld < oqb_dim-ind[nv-i]) {
            enlarge_hash_table(ht);
        }
        for (j = ind[nv-i]; j < oqb_dim; ++j) {
            memcpy(etmp, ht->ev[oqb[j]], (unsigned long)(nv+1)*sizeof(exp_t));;
            etmp[i]++;
            etmp[DEG]++;
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
            if (ht->ev[qb[j]][i+1] == 0) {
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

    exp_t ev[nv+1];
    memset(ev, 0, (unsigned long)(nv+1) * sizeof(exp_t));

    len_t *ind      = calloc((unsigned long)nv, sizeof(len_t));

    qb  = (hm_t *)calloc((unsigned long)1, sizeof(hm_t));
    qb[0] = check_lm_divisibility_and_insert_in_hash_table(
                            ev, (*htp), bs);

    while (nqbd > 0 && deg < max_deg) {
        nqb = realloc(nqb, (sum(ind, nv) + nv) * sizeof(hm_t));
        memset(nqb, 0, (sum(ind, nv) + nv) * sizeof(hm_t));
        nqbd = generate_new_basis_elements(htp, nqb, qb, qbd, ind, bs);
        qb = realloc(qb, (unsigned long)(qbd + nqbd) * sizeof(hm_t));
        /* printf("qb before adding\n");
         * for (len_t ii = 0; ii < qbd; ++ii) {
         *     printf("pos %u --> ", ii);
         *     for (len_t jj = 0; jj <= nv; ++jj) {
         *         printf("%d ", (*htp)->ev[qb[ii]][jj]);
         *     }
         *     printf("\n");
         * }
         * printf("nqb before adding\n");
         * for (len_t ii = 0; ii < nqbd; ++ii) {
         *     printf("pos %u --> ", ii);
         *     for (len_t jj = 0; jj <= nv; ++jj) {
         *         printf("%d ", (*htp)->ev[nqb[ii]][jj]);
         *     }
         *     printf("\n");
         * } */
        memcpy(qb+qbd, nqb, (unsigned long)nqbd * sizeof(hm_t));
        update_indices(ind, qb, qbd, nqbd, *htp);
        qbd   +=  nqbd;
        deg++;
        /* printf("qb after adding\n");
        for (len_t ii = 0; ii < qbd; ++ii) {
            printf("pos %u --> ", ii);
            for (len_t jj = 0; jj < (*htp)->evl; ++jj) {
                printf("%d ", (*htp)->ev[qb[ii]][jj]);
            }
            printf("\n");
        } */
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
        md_t *st,
        const bs_t *bs,
        const deg_t max_deg
        )
{
    len_t i, j, k;

    len_t qdim  = quotient_basis(qdp, bhtp, bs, max_deg);
    ht_t *bht   = *bhtp;
    ht_t *sht   = *shtp;

    bl_t ctr  = 0;

    check_enlarge_basis(sat, qdim, st);

    hm_t *qb  = *qdp;
    exp_t etmp[bht->evl];
    memset(etmp, 0, (unsigned long)bht->evl * sizeof(exp_t));

    /* remove elements that moved directly to kernel, i.e. monomials
     * added to the basis during the first stage of saturation
     * linear algebra */
    for (i = 0; i < sat->ld; ++i) {
        if (sat->hm[i] != NULL) {
            sat->hm[ctr] = sat->hm[i];
            sat->cf_32[ctr] = sat->cf_32[i];
            sat->hm[ctr][COEFFS]  = ctr;
            ctr++;
        }
    }
    sat->ld = ctr;

    ctr = 0;

    for (i = 0; i < sat->ld; ++i) {
        /* printf("qb[%u] = %u !=! %u (%u)\n", ctr, qb[ctr], i, sat->hm[i][MULT]); */
        while (i < sat->ld && ctr < qdim && qb[ctr] != sat->hm[i][MULT]) {
            free(sat->hm[i]);
            sat->hm[i]    = NULL;
            free(sat->cf_32[i]);
            sat->cf_32[i] = NULL;
            i++;
        }
        if (i < sat->ld) {
            sat->hm[ctr]          = sat->hm[i];
            sat->cf_32[ctr]       = sat->cf_32[i];
            sat->hm[ctr][COEFFS]  = ctr;
            ctr++;
        }
    }
    sat->lo = ctr;

    /* Here we apply the simplify idea from F4:
    * Any of these new elements is a multiple of an element
    * of lower degree already handled. Thus we can use this
    * "divisor" instead of the initial saturation element sat[0]. */
    for (i = ctr; i < qdim; ++i) {
        const hm_t m    = qb[i];
        const sdm_t ns  = ~bht->hd[qb[i]].sdm;
        j = sat->lo-1;
        /* printf("sat->lo %u | j %u\n", sat->lo, j); */
sat_restart:
        while (j > 0 && bht->hd[sat->hm[j][MULT]].sdm & ns) {
            j--;
        }
        for (k = 0; k < bht->evl; ++k) {
            if (bht->ev[qb[i]][k] < bht->ev[sat->hm[j][MULT]][k]) {
                j--;
                goto sat_restart;
            }
            etmp[k] = bht->ev[qb[i]][k] - bht->ev[sat->hm[j][MULT]][k];
        }
        const hi_t h      = bht->hd[m].val - bht->hd[sat->hm[j][MULT]].val;
        sat->hm[i]        = multiplied_poly_to_matrix_row(
                sht, bht, h, etmp, sat->hm[j]);
        sat->hm[i][MULT] = qb[i];
        deg_t deg = bht->hd[sat->hm[i][OFFSET]].deg;
        if (st->nev > 0) {
            const len_t len = sat->hm[i][LENGTH]+OFFSET;
            for (j = OFFSET+1; j < len; ++j) {
                if (deg < bht->hd[sat->hm[i][j]].deg) {
                    deg = bht->hd[sat->hm[i][j]].deg;
                }
            }
        }
        sat->hm[i][DEG] = deg;
        sat->cf_32[i]    = (cf32_t *)malloc(
                (unsigned long)sat->hm[j][LENGTH] * sizeof(cf32_t));
        memcpy(sat->cf_32[i], sat->cf_32[sat->hm[j][COEFFS]],
                (unsigned long)sat->hm[j][LENGTH] * sizeof(cf32_t));
        sat->hm[i][COEFFS] = i;
    }
    /* AFTER mapping the new sat elements to sht we can map the old ones.
     * New ones may be generated by the old ones, so we have to keep the
     * monomials of the old ones in bht until all new ones are generated! */
    for (i = 0; i < sat->lo; ++i) {
        while (sht->esz - sht->eld < sat->hm[i][LENGTH]) {
            enlarge_hash_table(sht);
        }
        for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
            sat->hm[i][j] = insert_in_hash_table(bht->ev[sat->hm[i][j]], sht);
        }
    }
    /* for (i = 0; i < sht->eld; ++i) {
        printf("sht[%u] = ", i);
        for (j = 0 ; j < sht->evl; ++j) {
            printf("%u ", sht->ev[i][j]);
        }
        printf("\n");
    } */
    sat->ld = qdim;
    st->new_multipliers = sat->ld - sat->lo;

    if (sat->mltdeg < max_deg) {
        sat->mltdeg = max_deg;
    }
    *bhtp = bht;
    *shtp = sht;
}

int core_f4sat(
        bs_t *gbs,
        bs_t *gsat,
        md_t *gmd,
        int32_t*errp
        )
{
    double ct = cputime();
    double rt = realtime();

    bs_t *bs  = gbs;
    ht_t *bht = bs->ht;
    md_t *st  = gmd;

    /* global saturation data */
    len_t sat_test  = 0;
    deg_t next_deg  = 0;
    deg_t sat_deg   = 0;

    /* elements to saturate input ideal with */
    bs_t *sat   = gsat;
    /* initialize multiplier of first element in sat to be the hash of
     * the all-zeroes exponent vector. */
    memset(bht->ev[0], 0, (unsigned long)bht->evl * sizeof(exp_t));
    sat->hm[0][MULT]  = insert_in_hash_table(bht->ev[0], bht);
    sat->ld = 1;

    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;
    /* elements tracking reduction steps of tbr, being
     * added to bs if we find a linear dependency on
     * elements in tbr. They represent the multipliers
     * of the elements we saturate with. */
    /* timings for one round */
    double rrt, crt;

    /* initialize update hash table, symbolic hash table */
    ht_t *sht = initialize_secondary_hash_table(bht, st);
    st->ht    = sht;

    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t *ps = initialize_pairset();
    st->ps   = ps;
    st->hcm  = (hi_t *)malloc(sizeof(hi_t));
    /* hashes-to-columns maps for multipliers in saturation step */
    hi_t *hcmm  = (hi_t *)malloc(sizeof(hi_t));

    int32_t round, i, j;

    st->max_gb_degree = INT32_MAX;

    /* elements of kernel in saturation step, to be added to basis bs */
    bs_t *kernel  = initialize_basis(st);

    /* reset bs->ld for first update process */
    bs->ld  = 0;
    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis_f4(ps, bs, bht, st, st->ngens);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    print_round_information_header(stdout, st);
    round = 1;
end_sat_step:
    for (; ps->ld > 0; ++round) {
        rrt = realtime();
        crt = cputime();
        st->max_bht_size  = st->max_bht_size > bht->esz ?
            st->max_bht_size : bht->esz;
        st->current_rd  = round;

        /* preprocess data for next reduction round */
        select_spairs_by_minimal_degree(mat, bs, st);
        symbolic_preprocessing(mat, bs, st);
        convert_hashes_to_columns(mat, st, sht);
        sort_matrix_rows_decreasing(mat->rr, mat->nru);
        sort_matrix_rows_increasing(mat->tr, mat->nrl);
        /* linear algebra, depending on choice, see set_function_pointers() */
        linear_algebra(mat, bs, bs, st);
        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            convert_sparse_matrix_rows_to_basis_elements(
                    -1, mat, bs, bht, sht, st);
            sat_test++;
        }
        clean_hash_table(sht);
        /* add lead monomials to trace, stores hashes in basis hash
         * table which is used in all upcoming F4 runs */
        /* if (mat->np > 0) {
         *     add_lms_to_trace(trace, bs, mat->np);
         *     trace->ltd++;
         * } */
        /* all rows in mat are now polynomials in the basis,
         * so we do not need the rows anymore */
        clear_matrix(mat);

        /* check redundancy only if input is not homogeneous */
        update_basis_f4(ps, bs, bht, st, mat->np);

        if (bs->constant  == 1) {
            printf("basis is constant\n");
            ps->ld  = 0;
            break;
        }
        clean_hash_table(sht);
        print_round_timings(stdout, st, rrt, crt);

        /* saturation step starts here */
        if ((bs->mltdeg >= sat->hm[0][DEG] && sat_test != 0) || ps->ld == 0) {
        /* if (bs->mltdeg - next_deg == 0 || ps->ld == 0) { */
        /* if (sat_done == 0 && (sat_test % 3 == 0 || ps->ld == 0)) { */
            if (st->nr_kernel_elts > 0 &&
                    ps->ld == 0 &&
                    is_zero_dimensional(bs, bht) &&
                    is_already_saturated(
                        bs, sat, mat, &bht, &sht, st)) {
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
                rrt  = realtime();
                crt  = cputime();
                /* printf("sat->deg %u\n", sat_deg); */
                update_multipliers(&qb, &bht, &sht, sat, st, bs, ii);
                /* check for monomial multiples of elements from saturation list */
                select_saturation(sat, mat, st, sht, bht);

                symbolic_preprocessing(mat, bs, st);

                /* It may happen that there is no reducer at all for the
                 * saturation elements, then nothing has to be done. */
                if (mat->nru > 0) {
                    if (st->info_level > 1) {
                        /* printf("kernel computation "); */
                        printf("%3u  compute kernel", sat_deg);
                    }
                    convert_hashes_to_columns_sat(mat, sat, st, sht);
                    convert_multipliers_to_columns(&hcmm, sat, st, bht);
                    sort_matrix_rows_decreasing(mat->rr, mat->nru);

                    compute_kernel_sat_ff_32(sat, mat, kernel, bs, st);

                    if (st->info_level > 1) {
                        printf("%56d new kernel elements", kernel->ld);
                        fflush(stdout);
                    }

                    if (kernel->ld > 0) {
                        if (st->info_level > 1) {
                            printf("\n                                               ");
                        }
                        clear_matrix(mat);
                        /* interreduce kernel */
                        copy_kernel_to_matrix(mat, kernel, sat->ld);
                        /* linear algebra, depending on choice, see set_function_pointers() */
                        probabilistic_sparse_linear_algebra_ff_32(mat, kernel, kernel, st);
                        /* linear_algebra(mat, kernel, st); */
                        /* columns indices are mapped back to exponent hashes */
                        if (mat->np > 0) {
                            /* we need to use the right hcmm */
                            hi_t *tmp = st->hcm;
                            st->hcm = hcmm;
                            convert_sparse_matrix_rows_to_basis_elements_use_sht(
                                    -1, mat, bs, bs->ht, st);
                            st->hcm = tmp;
                        }

                        st->nr_kernel_elts  +=  kernel->ld;
                        sat_test  = 0;
                        free_kernel_coefficients(kernel);
                        update_basis_f4(ps, bs, bht, st, mat->np);
                        kernel->ld  = 0;
                        if (st->info_level > 1) {
                            printf("   ");
                        }
                    }
                    /* all rows in mat are now polynomials in the basis,
                     * so we do not need the rows anymore */
                    convert_columns_to_hashes(sat, st, hcmm);
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
                        deg_t deg = bht->hd[sat->hm[i][OFFSET]].deg;
                        if (st->nev > 0) {
                            const len_t len = sat->hm[i][LENGTH]+OFFSET;
                            for (j = OFFSET+1; j < len; ++j) {
                                if (deg < bht->hd[sat->hm[i][j]].deg) {
                                    deg = bht->hd[sat->hm[i][j]].deg;
                                }
                            }
                        }
                        sat->hm[i][DEG] = deg;
                    }
                }
                clean_hash_table(sht);

                print_sat_round_timings(stdout, st, rrt, crt);
                if (bld != bs->ld) {
                    next_deg  = ii;
                    goto end_sat_step;
                }
            }
            next_deg  = sat_deg;
        }
    }
    print_round_information_footer(stdout, st);
    /* remove possible redudant elements */
    final_remove_redundant_elements(bs, st, bht);

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        /* note: bht will become sht, and sht will become NULL,
         * thus we need pointers */
        reduce_basis(bs, mat, st);
    }
    st->f4_rtime = realtime() - rt;
    st->f4_ctime = cputime() - ct;

    get_and_print_final_statistics(stdout, st, bs);
/*     printf("basis has  %u elements.\n", bs->lml);
 *
 *     for (i = 0; i < bs->lml; ++i) {
 *         for (j = 0; j < bht->nv; ++j) {
 *             printf("%u ", bht->ev[bs->hm[bs->lmps[i]][OFFSET]][j]);
 *         }
 *         printf("\n");
 *     } */

    /* free and clean up */
    free(hcmm);
    free(qb);

    free_basis_elements(sat);
    free_basis_without_hash_table(&sat);
    free_basis(&kernel);
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
