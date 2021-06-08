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
    len_t i, j, k;

    ht_t *ht  = *htp;
    const len_t nv  = ht->nv;
    len_t ctr = 0;
    hm_t htmp = 0;
    exp_t *etmp  = calloc((unsigned long)nv, sizeof(exp_t));

    printf("oqb_dim %u\n", oqb_dim);
    printf("oqb[0] %u\n", oqb[0]);
    for (i = nv; i > 0; --i) {
        while (ht->esz - ht->eld < oqb_dim-ind[nv-i]) {
            enlarge_hash_table(ht);
        }
        for (j = ind[nv-i]; j < oqb_dim; ++j) {
            memcpy(etmp, ht->ev[oqb[j]], (unsigned long)nv*sizeof(exp_t));;
            printf("before etmp %u || ", j);
            for (len_t ii = 0; ii < nv; ii++) {
                printf("%d ", etmp[ii]);
            }
            printf("\n");
            etmp[i-1]++;
            printf("etmp changing position %u -> ", i);
            for (len_t ii = 0; ii < nv; ii++) {
                printf("%d ", etmp[ii]);
            }
            printf("\n");
            printf("ctr is at the moment %u\n", ctr);
            cqb[ctr]  = check_lm_divisibility_and_insert_in_hash_table(
                            etmp, ht, bs
                            );
            if (cqb[ctr] != 0) {
                printf(" -- added at pos %u !\n", ctr);
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

    printf("LML\n");
    for (len_t i = 0; i < bs->lml; ++i) {
        for (len_t j = 0; j < nv; ++j) {
            printf("%d ", (*htp)->ev[bs->hm[bs->lmps[i]][OFFSET]][j]);
        }
        printf("\n");
    }
    /* clear ht-ev[0] */
    memset((*htp)->ev[0], 0, (unsigned long)nv * sizeof(exp_t));

    len_t *ind      = calloc((unsigned long)nv, sizeof(len_t));

    printf("max degree: %d\n", max_deg);
    qb  = (hm_t *)calloc((unsigned long)1, sizeof(hm_t));
    printf("nqb allocated %lu\n",sum(ind, nv) + 1);

    while (nqbd > 0 && deg < max_deg) {
        nqb = (hm_t *)calloc(sum(ind, nv) + nv, sizeof(hm_t));
        nqbd = generate_new_basis_elements(htp, nqb, qb, qbd, ind, bs);
        qb = realloc(qb, (unsigned long)qbd * nqbd * sizeof(hm_t));
        printf("qb before adding\n");
        for (len_t ii = 0; ii < qbd; ++ii) {
            printf("pos %u --> ", ii);
            for (len_t jj = 0; jj < nv; ++jj) {
                printf("%d ", (*htp)->ev[qb[ii]][jj]);
            }
            printf("\n");
        }
        printf("nqb before adding\n");
        for (len_t ii = 0; ii < nqbd; ++ii) {
            printf("pos %u --> ", ii);
            for (len_t jj = 0; jj < nv; ++jj) {
                printf("%d ", (*htp)->ev[nqb[ii]][jj]);
            }
            printf("\n");
        }
        memcpy(qb+qbd, nqb, (unsigned long)nqbd * sizeof(hm_t));
        update_indices(ind, qb, qbd, nqbd, *htp);
        qbd   +=  nqbd;
        deg++;
        printf("qb after adding\n");
        for (len_t ii = 0; ii < qbd; ++ii) {
            printf("pos %u --> ", ii);
            for (len_t jj = 0; jj < nv; ++jj) {
                printf("%d ", (*htp)->ev[qb[ii]][jj]);
            }
            printf("\n");
        }
    }
    free(nqb);
    free(ind);

    *qbp  = qb;
    return qbd;
}

static void update_multipliers(
        hm_t **qdp,
        bs_t **mulp,
        ht_t **htp,
        stat_t *st,
        const bs_t *bs
        )
{
    len_t i;

    len_t qdim  = quotient_basis(qdp, htp, bs, bs->mltdeg);
    ht_t *ht  = *htp;
    bs_t *mul = *mulp;

    mul->lo   = mul->ld;
    bl_t ctr  = mul->lo;

    check_enlarge_basis(mul, qdim);

    hm_t *qb  = *qdp;

    /* quotient monomials */
    for (i = 0; i < qdim; ++i) {
       mul->hm[ctr]           = realloc(mul->hm[ctr], (1+OFFSET)*sizeof(hm_t));
       mul->hm[ctr][LENGTH]   = 1;
       mul->hm[ctr][PRELOOP]  = 0;
       mul->hm[ctr][MULT]     = ctr;
       mul->hm[ctr][OFFSET]   = qb[i];
       mul->cf_32[ctr]        = realloc(mul->cf_32[ctr], sizeof(cf32_t));
       mul->cf_32[ctr][0]     = 1;
       ctr++;
    }
    mul->ld = ctr;
    st->new_multipliers = mul->ld - mul->lo;
    for (i = 0; i < mul->ld; ++i) {
        printf("%3u -> ", i);
        for (len_t j = 0; j < ht->nv; ++j) {
            printf("%d ", ht->ev[mul->hm[i][OFFSET]][j]);
        }
        printf("\n");
    }

    *mulp = mul;
    *htp  = ht;
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
    /* current quotient basis up to max lm degree in intermediate basis */
    hm_t *qb    = NULL;
    /* elements to compute the normal forms of */
    bs_t *tbr   = initialize_basis_ff_32(10);
    /* elements tracking reduction steps of tbr, being
     * added to bs if we find a linear dependency on
     * elements in tbr. They represent the multipliers
     * of the elements we saturate with. */
    bs_t *mul   = initialize_basis_ff_32(10);
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
        /* check for new elements to be tested for adding saturation
         * information to the intermediate basis */
        update_multipliers(&qb, &mul, &bht, st, bs);
        /* check for monomial multiples of elements from saturation list */
        select_saturation(sat, mul, mat, st, sht, bht);

        symbolic_preprocessing(mat, bs, st, sht, NULL, bht);
        if (st->info_level > 1) {
            printf("nf computation data");
        }
        convert_hashes_to_columns(&hcm, mat, st, sht);
        convert_multipliers_to_columns(&hcmm, mul, st, bht);
        sort_matrix_rows_decreasing(mat->rr, mat->nru);

        /* linear algebra, depending on choice, see set_function_pointers() */
        exact_sparse_linear_algebra_sat_ff_32(mat, mul, sat, bs, st);
        /* columns indices are mapped back to exponent hashes */
        return_normal_forms_to_basis(
                mat, tbr, bht, sht, hcm, st);

        /* all rows in mat are now polynomials in the basis,
        * so we do not need the rows anymore */
        clear_matrix(mat);
        clean_hash_table(sht);

        /* free multiplier list
         * todo: keep already reduced, nonzero elements for next round */
        free_basis_elements(mul);

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

    *bsp  = bs;
    *satp = sat;
    *bhtp = bht;
    *stp  = st;

    /* free and clean up */
    free(tbr);
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
