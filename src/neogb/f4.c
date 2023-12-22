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
        md_t *md,
        const ht_t * const ht
        )
{
    len_t i, j;
    for (i = 0; i < bs->lml; ++i) {
        hm_t nch = bs->hm[bs->lmps[i]][OFFSET];
        for (j = 0; j < i; ++j) {
            if (bs->red[bs->lmps[j]] == 0
                    && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                    ) {
                bs->red[bs->lmps[i]]  =   1;
                md->num_redundant++;
                break;
            }
        }
        for (j = i+1; j < bs->lml; ++j) {
            if (bs->red[bs->lmps[j]] == 0
                    && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                    ) {
                bs->red[bs->lmps[i]]  =   1;
                md->num_redundant++;
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
        md_t *st
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
    exp_t etmp[bht->evl];
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
    symbolic_preprocessing(mat, bs, st, sht, bht);
    /* no known pivots, we need mat->ncl = 0, so set all indices to 1 */
    for (i = 0; i < sht->eld; ++i) {
        sht->hd[i].idx = 1;
    }

    /* generate hash <-> column mapping */
    if (st->info_level > 1) {
        printf("reduce intermediate basis ");
        fflush(stdout);
    }
    convert_hashes_to_columns(mat, st, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    /* do the linear algebra reduction, do NOT free basis data */
    interreduce_matrix_rows(mat, bs, st, 0);
    /* remap rows to basis elements (keeping their position in bs) */
    convert_sparse_matrix_rows_to_basis_elements(mat, bs, bht, sht, st);

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
        md_t *md
        )
{
    /* timings */
    double ct, rt;
    ct = cputime();
    rt = realtime();

    len_t i, j, k;

    ht_t *bht   = bs->ht;
    ht_t *sht   = md->ht;
    exp_t etmp[bht->evl];
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
    symbolic_preprocessing(mat, bs, md);
    /* no known pivots, we need mat->ncl = 0, so set all indices to 1 */
    for (i = 0; i < sht->eld; ++i) {
        sht->hd[i].idx = 1;
    }

    /* generate hash <-> column mapping */
    if (md->info_level > 1) {
        printf("reduce final basis ");
        fflush(stdout);
    }
    convert_hashes_to_columns(mat, md, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    /* do the linear algebra reduction and free basis data afterwards */
    interreduce_matrix_rows(mat, bs, md, 1);

    convert_sparse_matrix_rows_to_basis_elements(
            1, mat, bs, bht, sht, md);

    /* set sht = NULL, otherwise we might run in a double-free
     * of sht and bht at the end */
    /* sht     = NULL;
    md->ht  = sht; */

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

    print_round_timings(stdout, md, rt, ct);
    print_round_information_footer(stdout, md);
}

static int32_t initialize_f4(
        bs_t **lbsp,
        md_t **lmdp,
        mat_t **matp,
        md_t *gmd,
        bs_t *gbs,
        len_t fc
        )
{
    bs_t *bs     = *lbsp;
    mat_t *mat   = *matp;
    int32_t done = 0;
    md_t *md     = *lmdp;

    md      = copy_meta_data(gmd, fc);
    md->fc  = fc;
    md->hcm = (hi_t *)malloc(sizeof(hi_t));

    if (gmd->fc != fc) {
        reset_function_pointers(fc, md->laopt);
        bs = copy_basis_mod_p(gbs, md);
        if (md->laopt < 40) {
            if (md->trace_level != APPLY_TRACER) {
                md->trace_level = LEARN_TRACER;
            }
        }
    } else {
        bs = gbs;
        md->trace_level = NO_TRACER;
    }
    normalize_initial_basis(bs, fc);
    md->ht = initialize_secondary_hash_table(bs->ht, md);

    /* matrix holding sparse information generated
       during symbolic preprocessing */
    mat = (mat_t *)calloc(1, sizeof(mat_t));

    if (md->trace_level != APPLY_TRACER) {
        /* pair set */
        md->ps = initialize_pairset();

        /* reset bs->ld for first update process */
        bs->ld  = 0;
    } else {
        bs->ld = md->ngens;
    }

    /* TODO: make this a command line argument */
    md->max_gb_degree = INT32_MAX;

    /* link tracer into basis */
    if (md->trace_level == LEARN_TRACER) {
        md->tr = initialize_trace(bs, md);
    }


    /* move input generators to basis and generate first spairs.
       always check redundancy since input generators may be redundant
       even so they are homogeneous. */
    if (md->trace_level != APPLY_TRACER) {
        md->np = md->ngens;
        done   = update(bs, md);
    }

    *lbsp  = bs;
    *matp  = mat;
    *lmdp  = md;
    
    return done;
}

static int32_t compute_new_elements(
    mat_t *mat,
    bs_t *bs,
    md_t *md,
    int32_t *errp
    )
{
    len_t i;

    ht_t *ht  = bs->ht;
    ht_t *sht = md->ht;

    convert_hashes_to_columns(mat, md, sht);
    sort_matrix_rows_decreasing(mat->rr, mat->nru);
    linear_algebra(mat, bs, bs, md);

    /* check for bad prime */
    if (md->trace_level == APPLY_TRACER) {
        if (mat->np != md->tr->td[md->trace_rd].nlm) {
            if (md->info_level > 0) {
                fprintf(stderr, "Wrong number of new elements, bad prime.");
            }
            *errp = 1;
            return 1;
        }
    }
    /* columns indices are mapped back to exponent hashes */
    if (mat->np > 0) {
        convert_sparse_matrix_rows_to_basis_elements(
                -1, mat, bs, ht, sht, md);
    }
    clean_hash_table(sht);
    /* all rows in mat are now polynomials in the basis,
     * so we do not need the rows anymore */
    clear_matrix(mat);

    /* check for bad prime */
    if (md->trace_level == APPLY_TRACER) {
        for (i = 0; i < md->np; ++i) {
            if (bs->hm[bs->ld+i][OFFSET] != md->tr->td[md->trace_rd].nlms[i]) {
                if (md->info_level > 0) {
                fprintf(stderr, "Wrong leading term for new element %u/%u, bad prime.",
                        i, mat->np);
                }
                *errp = 2;
                return 1;
            }
        }
    }
    if (md->trace_level == LEARN_TRACER && md->np > 0) {
        add_lms_to_trace(md->tr, bs, md->np);
        md->tr->ltd++;
    }
    if (md->trace_level == APPLY_TRACER) {
        bs->ld += md->np;
        md->trace_rd++;
        if (md->trace_rd >= md->tr->ltd) {
            return 1;
        }
    }

    return 0;
}

static void process_redundant_elements(
        bs_t *bs,
        md_t *md
        )
{
    len_t i, j;

    ht_t *ht = bs->ht;

    if (md->trace_level != APPLY_TRACER) {
        for (i = 0; i < bs->lml; ++i) {
            hm_t nch = bs->hm[bs->lmps[i]][OFFSET];
            for (j = 0; j < i; ++j) {
                if (bs->red[bs->lmps[j]] == 0
                        && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                        ) {
                    bs->red[bs->lmps[i]]  =   1;
                    md->num_redundant++;
                    break;
                }
            }
            for (j = i+1; j < bs->lml; ++j) {
                if (bs->red[bs->lmps[j]] == 0
                        && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                        ) {
                    bs->red[bs->lmps[i]]  =   1;
                    md->num_redundant++;
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

        /* At the moment we do not directly remove the eliminated polynomials from
         * the resulting basis. */
#if 0
        if (md->nev > 0) {
            j = 0;
            for (i = 0; i < bs->lml; ++i) {
                if (ht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
                    bs->lm[j]   = bs->lm[i];
                    bs->lmps[j] = bs->lmps[i];
                    ++j;
                }
            }
            bs->lml = j;
        }
#endif
    }
    if (md->trace_level == APPLY_TRACER) {
        /* apply non-redundant basis data from trace to basis
         * before interreduction */
        bs->lml  = md->tr->lml;
        free(bs->lmps);
        bs->lmps = (bl_t *)calloc((unsigned long)bs->lml,
                sizeof(bl_t));
        memcpy(bs->lmps, md->tr->lmps,
                (unsigned long)bs->lml * sizeof(bl_t));
        free(bs->lm);
        bs->lm   = (sdm_t *)calloc((unsigned long)bs->lml,
                sizeof(sdm_t));
        memcpy(bs->lm, md->tr->lm,
                (unsigned long)bs->lml * sizeof(sdm_t));
    }

    if (md->trace_level == LEARN_TRACER) {
        /* store information in trace */
        md->tr->lml  = bs->lml;
        md->tr->lmps = (bl_t *)calloc((unsigned long)md->tr->lml,
                sizeof(bl_t));
        memcpy(md->tr->lmps, bs->lmps,
                (unsigned long)md->tr->lml * sizeof(bl_t));
        md->tr->lm   = (sdm_t *)calloc((unsigned long)md->tr->lml,
                sizeof(sdm_t));
        memcpy(md->tr->lm, bs->lm,
                (unsigned long)md->tr->lml * sizeof(sdm_t));
        /* do not track the final reduction step */
    }
}

static void reduce_final_basis(
        bs_t *bs,
        mat_t *mat,
        md_t *md
        )
{
    if (md->reduce_gb) {
        reduce_basis(bs, mat, md);
    }
}

static void free_local_data(
        mat_t **matp,
        md_t **mdp
        )
{
    free_meta_data(mdp);

    free(*matp);
    *matp = NULL;
}

static void finalize_f4(
        md_t *gmd,
        bs_t *gbs,
        bs_t **bsp,
        md_t **lmdp,
        mat_t **matp,
        int32_t err
        )
{
    if (err > 0) {
        free_basis_and_only_local_hash_table_data(bsp);
    }

    if ((*lmdp)->trace_level == LEARN_TRACER) {
        gmd->tr = (*lmdp)->tr;
        gmd->trace_level = APPLY_TRACER;
    }
    free_local_data(matp, lmdp);
}

bs_t *core_f4(
        bs_t *gbs,
        md_t *gmd,
        int32_t *errp,
        const len_t fc
        )
{
    double ct = cputime();
    double rt = realtime();

    bs_t *bs   = NULL;
    md_t *md   = NULL;
    mat_t *mat = NULL;

    /* marker for end of computation */
    int32_t done = 0;

    /* timings for one round */
    double rrt, crt;

    done = initialize_f4(&bs, &md, &mat, gmd, gbs, fc);


    /* let's start the f4 rounds, we are done when no more spairs
       are left in the pairset or if we found a constant in the basis. */
    print_round_information_header(stdout, md);
    
    /* reset error */
    *errp = 0;
    while (!done) {
        rrt = realtime();
        crt = cputime();
        md->max_bht_size = md->max_bht_size > bs->ht->esz ?
            md->max_bht_size : bs->ht->esz;

        done = preprocessing(mat, bs, md);

        if (!done) {
            done = compute_new_elements(mat, bs, md, errp);
        }
        if (!done && md->trace_level != APPLY_TRACER) {
            done = update(bs, md);
        }

        print_round_timings(stdout, md, rrt, crt);
    }
    if (*errp > 0) {
        free_basis_and_only_local_hash_table_data(&bs);
    } else {
        print_round_information_footer(stdout, md);

        /* remove possible redudant elements */
        process_redundant_elements(bs, md);

        /* reduce final basis? */
        reduce_final_basis(bs, mat, md);

        md->f4_rtime = realtime() - rt;
        md->f4_ctime = cputime() - ct;

        get_and_print_final_statistics(stdout, md, bs);

        finalize_f4(gmd, gbs, &bs, &md, &mat, *errp);
    }
    return bs;
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
    md_t **stp
    )
{

    bs_t *bs  = *bsp;
    ht_t *bht = *bhtp;
    md_t *st  = *stp;

    st->nterms_basis  = export_data(
        bld, blen, bexp, bcf, mallocp, bs, bht, st);
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
int64_t export_f4(
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
    bs_t *bs  = NULL;
    ht_t *bht = NULL;
    md_t *md  = NULL;

    int success = 0;

    const int32_t use_signatures    =   0;
    success = initialize_gba_input_data(&bs, &bht, &md,
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

    int err = 0;
    bs = core_f4(bs, md, &err, field_char);

    if (err) {
        printf("Problem with F4, stopped computation.\n");
        exit(1);
    }

    int64_t nterms  = export_results_from_f4(bld, blen, bexp,
            bcf, mallocp, &bs, &bht, &md);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->f4_ctime = ct1 - ct0;
    md->f4_rtime = rt1 - rt0;

    get_and_print_final_statistics(stderr, md, bs);

    /* free and clean up */
    free_shared_hash_data(bht);
    if (bs != NULL) {
        free_basis(&bs);
    }

    free(md);
    md    = NULL;

    return nterms;
}
