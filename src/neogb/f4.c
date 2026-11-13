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
#include "../msolve/streams.h"

static void final_remove_redundant_elements(
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

    int64_t i;
    int64_t len = 0;
    for (i = 0; i < ngens; ++i) {
        len += (int64_t)lens[i];
    }

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
        mpz_t **cfs = (mpz_t **)bcf;
        for (i = 0; i < len; ++i) {
            mpz_clear((*cfs)[i]);
        }
        (*freep)(*cfs);
        *cfs = NULL;
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
    free(mat->rrb);
    mat->rrb = NULL;
    free(mat->tr);
    mat->tr  = NULL;
    free(mat->cf_8);
    mat->cf_8 = NULL;
    free(mat->cf_16);
    mat->cf_16  = NULL;
    free(mat->cf_32);
    mat->cf_32  = NULL;
    free(mat->cf_64);
    mat->cf_64  = NULL;
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
        fprintf(VERBSTREAM, "reduce intermediate basis ");
        fflush(VERBSTREAM);
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
                    case 64:
                        free(bs->cf_64[i]);
                        bs->cf_64[i] = bs->cf_64[bs->ld+j];
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
            case 64:
                bs->cf_64[i] = NULL;
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
        fprintf(VERBSTREAM, "%13.2f sec\n", rt1-rt0);
    }

    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "-------------------------------------------------\
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

    md->in_final_reduction_step = 1;

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
        fprintf(VERBSTREAM, "reduce final basis ");
        fflush(VERBSTREAM);
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

    md->in_final_reduction_step = 0;

    print_round_timings(VERBSTREAM, md, rt, ct);
    print_round_information_footer(VERBSTREAM, md);
}

static int32_t initialize_f4(
        bs_t **lbsp,
        md_t **lmdp,
        mat_t **matp,
        md_t *gmd,
        bs_t *gbs,
        uint64_t fc
        )
{
    bs_t *bs     = *lbsp;
    mat_t *mat   = *matp;
    int32_t done = 0;
    md_t *md     = *lmdp;

    md      = copy_meta_data(gmd, fc);
    md->fc  = fc;
    md->hcm = (hi_t *)malloc(sizeof(hi_t));
    md->mat_nterms = 0;
    md->max_mat_nterms = 0;

    if (gmd->fc != fc) {
        if (gmd->f4_qq_round != 2) {
            reset_function_pointers(fc, md->laopt);
        }
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
        start_trace_cache(md, bs);
    }

    /* TODO: make this a command line argument */
    md->max_gb_degree = INT32_MAX;

    /* link tracer into basis */
    if (md->trace_level == LEARN_TRACER) {
        md->tr = initialize_trace(bs, md);
        md->min_deg_in_first_deg_fall = INT32_MAX;
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

/* frees the mat->np new rows of mat->tr and their coefficient arrays
 * (used when they are not moved to the basis) */
static void free_new_matrix_rows(
    mat_t *mat,
    const md_t * const md
    )
{
    len_t i;
    if (mat->tr == NULL) {
        return;
    }
    for (i = 0; i < mat->np; ++i) {
        hm_t *row = mat->tr[i];
        if (row == NULL) {
            continue;
        }
        switch (md->ff_bits) {
            case 8:
                free(mat->cf_8[row[COEFFS]]);
                mat->cf_8[row[COEFFS]] = NULL;
                break;
            case 16:
                free(mat->cf_16[row[COEFFS]]);
                mat->cf_16[row[COEFFS]] = NULL;
                break;
            case 32:
                free(mat->cf_32[row[COEFFS]]);
                mat->cf_32[row[COEFFS]] = NULL;
                break;
            case 64:
                free(mat->cf_64[row[COEFFS]]);
                mat->cf_64[row[COEFFS]] = NULL;
                break;
            default:
                break;
        }
        free(row);
        mat->tr[i] = NULL;
    }
    mat->np = 0;
}

/* Over the integers the rows of a matrix generate a lattice which the
 * linear algebra keeps, see la_zz.c. Returns the lead columns and basis
 * indices of the rows of mat, sorted, as (column << 32) | index. */
static uint64_t *zz_matrix_rows(
    const mat_t * const mat,
    len_t *nrp
    )
{
    len_t i;

    const len_t nr = mat->nru + mat->nrl;
    uint64_t *r = (uint64_t *)malloc((unsigned long)(nr+1) * sizeof(uint64_t));
    for (i = 0; i < mat->nru; ++i) {
        r[i] = ((uint64_t)mat->rr[i][OFFSET] << 32) | mat->rr[i][BINDEX];
    }
    for (i = 0; i < mat->nrl; ++i) {
        r[mat->nru+i] = ((uint64_t)mat->tr[i][OFFSET] << 32) | mat->tr[i][BINDEX];
    }
    qsort(r, (unsigned long)nr, sizeof(uint64_t), cmp_uint64);
    *nrp = nr;
    return r;
}

/* stores in md->zz_roff and md->zz_rbi for each new pivot of mat the basis
 * indices of the rows r of the matrix with its lead column; the new
 * elements become basis elements in the order of decreasing positions in
 * mat->tr, see convert_sparse_matrix_rows_to_basis_elements */
static void zz_rows_of_new_elements(
    md_t *md,
    const mat_t * const mat,
    const uint64_t * const r,
    const len_t nr
    )
{
    len_t i, k;

    const len_t np = mat->np;
    md->zz_roff = (len_t *)malloc((unsigned long)(np+1) * sizeof(len_t));
    len_t *lo = (len_t *)malloc((unsigned long)(np+1) * sizeof(len_t));
    md->zz_roff[0] = 0;
    for (k = 0; k < np; ++k) {
        const uint64_t c = (uint64_t)mat->tr[np-1-k][OFFSET];
        /* first row with column c */
        len_t a = 0, b = nr;
        while (a < b) {
            const len_t m = a + (b - a) / 2;
            if ((r[m] >> 32) < c) {
                a = m + 1;
            } else {
                b = m;
            }
        }
        lo[k] = a;
        for (i = a; i < nr && (r[i] >> 32) == c; ++i);
        md->zz_roff[k+1] = md->zz_roff[k] + (i - a);
    }
    md->zz_rbi = (bl_t *)malloc((unsigned long)(md->zz_roff[np]+1) * sizeof(bl_t));
    for (k = 0; k < np; ++k) {
        for (i = 0; i < md->zz_roff[k+1] - md->zz_roff[k]; ++i) {
            md->zz_rbi[md->zz_roff[k]+i] = (bl_t)(r[lo[k]+i] & 0xFFFFFFFF);
        }
    }
    free(lo);
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

    if (md->tc_state != 1) {
        convert_hashes_to_columns(mat, md, sht);
        sort_matrix_rows_decreasing(mat->rr, mat->nru);
        if (md->tc_state == 2) {
            store_matrix_in_trace_cache(md->tr, md->trace_rd,
                    mat, sht, md->hcm);
        }
    }
    /* the rows of the matrix with the lead monomials of the new elements,
     * see update_basis_zz */
    uint64_t *zr = NULL;
    len_t nzr    = 0;
    free(md->zz_roff);
    free(md->zz_rbi);
    md->zz_roff = NULL;
    md->zz_rbi  = NULL;
    if (md->zz != 0 && md->tc_state != 1) {
        zr = zz_matrix_rows(mat, &nzr);
    }
    linear_algebra(mat, bs, bs, md);
    if (zr != NULL) {
        zz_rows_of_new_elements(md, mat, zr, nzr);
        free(zr);
    }

    /* a lead coefficient other than +-1 ends a computation which has to be
     * the image of an integer computation, see md->zz_unit */
    if (md->zz_unit != 0 && md->zz_nonunit != 0) {
        /* the linear algebra has stored this round in the trace, it is
         * freed with the trace */
        if (md->trace_level == LEARN_TRACER) {
            md->tr->ltd++;
        }
        free_new_matrix_rows(mat, md);
        clear_matrix(mat);
        md->tc_state = 0;
        *errp = 3;
        return 1;
    }

    /* check for bad prime */
    if (md->trace_level == APPLY_TRACER) {
        if (mat->np != md->tr->td[md->trace_rd].nlm) {
            if (md->info_level > 0) {
                fprintf(ERRSTREAM, "Wrong number of new elements, bad prime.");
            }
            /* new rows (and their coefficient arrays) are owned by the
             * matrix at this point, release them before bailing out */
            free_new_matrix_rows(mat, md);
            clear_matrix(mat);
            md->tc_state = 0;
            *errp = 1;
            return 1;
        }
    }
    /* columns indices are mapped back to exponent hashes */
    if (mat->np > 0) {
        convert_sparse_matrix_rows_to_basis_elements(
                -1, mat, bs, ht, sht, md);
    }
    md->tc_state = 0;
    clean_hash_table(sht);
    /* all rows in mat are now polynomials in the basis,
     * so we do not need the rows anymore */
    clear_matrix(mat);

    /* check for bad prime */
    if (md->trace_level == APPLY_TRACER) {
        for (i = 0; i < md->np; ++i) {
            if (bs->hm[bs->ld+i][OFFSET] != md->tr->td[md->trace_rd].nlms[i]) {
                if (md->info_level > 0) {
                fprintf(ERRSTREAM, "Wrong leading term for new element %u/%u, bad prime.",
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
                        && (md->zz == 0
                            || zz_lt_divides(bs, bs->lmps[j], bs->lmps[i], ht))
                        ) {
                    bs->red[bs->lmps[i]]  =   1;
                    md->num_redundant++;
                    break;
                }
            }
            for (j = i+1; j < bs->lml; ++j) {
                if (bs->red[bs->lmps[j]] == 0
                        && check_monomial_division(nch, bs->hm[bs->lmps[j]][OFFSET], ht)
                        && (md->zz == 0
                            || zz_lt_divides(bs, bs->lmps[j], bs->lmps[i], ht))
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
        md->in_final_reduction_step = 1;

        /* timings */
        double ct, rt;
        ct = cputime();
        rt = realtime();

        len_t i;

        ht_t *bht   = bs->ht;
        ht_t *sht   = md->ht;

        /* add all non-redundant basis elements as matrix rows, their
         * terms are inserted into sht in parallel (multiplied by 1) */
        const len_t lml = bs->lml;
        mat->tr = (hm_t **)malloc((unsigned long)lml * sizeof(hm_t *));
        hm_t **polys  = (hm_t **)malloc((unsigned long)lml * sizeof(hm_t *));
        val_t *hs     = (val_t *)calloc((unsigned long)lml, sizeof(val_t));
        exp_t *ems    = (exp_t *)calloc(
                (unsigned long)lml * bht->evl, sizeof(exp_t));
        for (i = 0; i < lml; ++i) {
            polys[i] = bs->hm[bs->lmps[i]];
        }
        multiplied_polys_to_matrix_rows(mat->tr, lml, hs, ems, polys,
                sht, bht, md->nthrds);
        for (i = 0; i < lml; ++i) {
            sht->hd[mat->tr[i][OFFSET]].idx  = 1;
        }
        free(polys);
        free(hs);
        free(ems);
        mat->nr = mat->nrl = mat->sz = bs->lml;
        mat->nc = 0;

        symbolic_preprocessing(mat, bs, md);
        /* over finite fields the reducers are the known pivots on the
         * left side of the matrix as in any other round, the basis
         * elements are the rows to be reduced. They are reduced by the
         * known pivots in parallel, afterwards the interreduction of the
         * new pivots reduces them by each other. The rational linear
         * algebra needs mat->ncl = 0, so we set all indices to 1. The
         * integer linear algebra uses the known pivots as well. */
        if (md->ff_bits == 0 && md->zz == 0) {
            for (i = 0; i < sht->eld; ++i) {
                sht->hd[i].idx = 1;
            }
        }

        /* generate hash <-> column mapping */
        if (md->info_level > 1) {
            fprintf(VERBSTREAM, "reduce final basis ");
            fflush(VERBSTREAM);
        }
        convert_hashes_to_columns(mat, md, sht);
        mat->nc = mat->ncl + mat->ncr;

        sort_matrix_rows_decreasing(mat->rr, mat->nru);
        sort_matrix_rows_increasing(mat->tr, mat->nrl);

        exact_linear_algebra(mat, bs, bs, md);

        /* the basis keeps the order of the rows before the reduction */
        sort_matrix_rows_increasing(mat->tr, mat->np);

        free_basis_elements(bs);

        convert_sparse_matrix_rows_to_basis_elements(
            0, mat, bs, bht, sht, md);

        bs->ld = bs->lml = mat->np;
        clear_matrix(mat);

        for (i = 0; i < bs->ld; ++i) {
            bs->lmps[i] = i;
            bs->lm[i]   = bht->hd[bs->hm[i][OFFSET]].sdm;
        }

        md->in_final_reduction_step = 0;

        /* timings */
        print_round_timings(VERBSTREAM, md, rt, ct);
        print_round_information_footer(VERBSTREAM, md);
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
        initialize_trace_cache(gmd->tr, *lmdp);
    }
    /* not changed when applying the trace, which may be done for several
     * primes in parallel */
    if ((*lmdp)->trace_level != APPLY_TRACER) {
        gmd->min_deg_in_first_deg_fall = (*lmdp)->min_deg_in_first_deg_fall;
        gmd->mat_nterms = (*lmdp)->mat_nterms;
        gmd->max_mat_nterms = gmd->max_mat_nterms > (*lmdp)->max_mat_nterms ?
            gmd->max_mat_nterms : (*lmdp)->max_mat_nterms;
    }
    free_local_data(matp, lmdp);
}

static bs_t *run_f4(
        bs_t *gbs,
        md_t *gmd,
        int32_t *errp,
        const uint64_t fc
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
    print_round_information_header(VERBSTREAM, md);

    /* reset error */
    *errp = 0;
    /* number of rounds done, see md->zz_max_rd */
    int32_t rd = 0;
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
            rd++;
            if (md->zz_max_rd > 0 && rd >= md->zz_max_rd) {
                done = 1;
            }
        }

        print_round_timings(VERBSTREAM, md, rrt, crt);
    }
    if (md->zz_nonunit != 0) {
        gmd->zz_nonunit = 1;
        gmd->zz_rd      = rd;
    }
    if (*errp > 0) {
        /* when learning the trace the basis shares the hash table of the
         * global basis, see copy_basis_mod_p */
        if (bs->ht == gbs->ht) {
            free_basis_without_hash_table(&bs);
        } else {
            free_basis_and_only_local_hash_table_data(&bs);
        }
        /* a trace learnt so far is not used */
        if (md->trace_level == LEARN_TRACER) {
            free_trace(&md->tr);
        }
        /* also release the local meta data (symbolic hash table, pair
         * set, ...) and the matrix: finalize_f4 is not called here */
        free_local_data(&mat, &md);
    } else {
        print_round_information_footer(VERBSTREAM, md);

        /* remove possible redudant elements */
        process_redundant_elements(bs, md);

        /* reduce final basis? */
        reduce_final_basis(bs, mat, md);

        md->f4_rtime = realtime() - rt;
        md->f4_ctime = cputime() - ct;

        get_and_print_final_statistics(VERBSTREAM, md, bs);

        finalize_f4(gmd, gbs, &bs, &md, &mat, *errp);
    }
    return bs;
}

/* The exponent vectors w.r.t. an elimination order with nev eliminated
 * variables are [deg_1, x_1, ..., x_nev, deg_2, x_nev+1, ..., x_nv],
 * w.r.t. DRL they are [deg, x_1, ..., x_nv], see data.h. */
static inline void block_to_drl_exponents(
        exp_t *f,
        const exp_t * const e,
        const len_t nev,
        const len_t nv
        )
{
    len_t i;

    f[0] = e[0] + e[nev+1];
    for (i = 1; i <= nev; ++i) {
        f[i] = e[i];
    }
    for (i = nev+1; i <= nv; ++i) {
        f[i] = e[i+1];
    }
}

static inline void drl_to_block_exponents(
        exp_t *e,
        const exp_t * const f,
        const len_t nev,
        const len_t nv
        )
{
    len_t i;

    e[0]     = 0;
    e[nev+1] = 0;
    for (i = 1; i <= nev; ++i) {
        e[i]  = f[i];
        e[0] += f[i];
    }
    for (i = nev+1; i <= nv; ++i) {
        e[i+1]    = f[i];
        e[nev+1] += f[i];
    }
}

/* index of the monomial a in ht, 0 if it is not contained in ht */
static inline hi_t find_in_hash_table(
        const exp_t * const a,
        const ht_t * const ht
        )
{
    len_t j;
    hi_t k;
    val_t h = 0;

    for (j = 0; j < ht->evl; ++j) {
        h +=  ht->rn[j] * a[j];
    }
    return is_contained_in_hash_table(a, ht, h, &k) ? k : 0;
}

/* Copies the elements src->hm[idx[0]], ..., src->hm[idx[n-1]] to the
 * empty basis dst whose hash table may use another monomial order,
 * nev is the number of eliminated variables of the elimination order
 * of src resp. dst. The terms are sorted w.r.t. the order of dst, the
 * coefficients are not normalized. If lookup is set the hash table of
 * dst is not changed: the terms are only looked up, if one of them is
 * missing 0 is returned and dst holds the elements copied so far
 * (without coefficients). Otherwise 1 is returned. */
static int copy_basis_to_other_order(
        bs_t *dst,
        const bs_t * const src,
        const bl_t * const idx,
        const len_t n,
        const len_t nev,
        const md_t * const md,
        const int lookup
        )
{
    len_t i, j;

    const ht_t * const sht = src->ht;
    ht_t *dht              = dst->ht;
    const len_t nv         = dht->nv;
    const int to_drl       = dht->ebl == 0;

    exp_t *e = (exp_t *)malloc((unsigned long)dht->evl * sizeof(exp_t));

    check_enlarge_basis(dst, n, md);

    for (i = 0; i < n; ++i) {
        const hm_t * const s = src->hm[idx[i]];
        const len_t len      = s[LENGTH];
        while (!lookup && len >= dht->esz - dht->eld) {
            enlarge_hash_table(dht);
        }
        hm_t *d     = (hm_t *)malloc(((unsigned long)len+OFFSET) * sizeof(hm_t));
        d[COEFFS]   = i;
        d[PRELOOP]  = len % UNROLL;
        d[LENGTH]   = len;
        for (j = 0; j < len; ++j) {
            if (to_drl) {
                block_to_drl_exponents(e, sht->ev[s[OFFSET+j]], nev, nv);
            } else {
                drl_to_block_exponents(e, sht->ev[s[OFFSET+j]], nev, nv);
            }
            if (!lookup) {
                d[OFFSET+j] = insert_in_hash_table(e, dht);
            } else if ((d[OFFSET+j] = find_in_hash_table(e, dht)) == 0) {
                free(d);
                dst->ld = i;
                free(e);
                return 0;
            }
        }
        dst->hm[i]  = d;
        dst->red[i] = 0;

        switch (md->ff_bits) {
            case 0:
                /* the source may also be a basis over a prime field */
                dst->cf_qq[i] = (mpz_t *)malloc((unsigned long)len * sizeof(mpz_t));
                for (j = 0; j < len; ++j) {
                    if (src->cf_qq != NULL) {
                        mpz_init_set(dst->cf_qq[i][j], src->cf_qq[s[COEFFS]][j]);
                    } else if (src->cf_32 != NULL) {
                        mpz_init_set_ui(dst->cf_qq[i][j], src->cf_32[s[COEFFS]][j]);
                    } else if (src->cf_16 != NULL) {
                        mpz_init_set_ui(dst->cf_qq[i][j], src->cf_16[s[COEFFS]][j]);
                    } else if (src->cf_8 != NULL) {
                        mpz_init_set_ui(dst->cf_qq[i][j], src->cf_8[s[COEFFS]][j]);
                    } else {
                        mpz_init_set_ui(dst->cf_qq[i][j], src->cf_64[s[COEFFS]][j]);
                    }
                }
                sort_terms_qq(&(dst->cf_qq[i]), &(dst->hm[i]), dht);
                break;
            case 8:
                dst->cf_8[i] = (cf8_t *)malloc((unsigned long)len * sizeof(cf8_t));
                memcpy(dst->cf_8[i], src->cf_8[s[COEFFS]],
                        (unsigned long)len * sizeof(cf8_t));
                sort_terms_ff_8(&(dst->cf_8[i]), &(dst->hm[i]), dht);
                break;
            case 16:
                dst->cf_16[i] = (cf16_t *)malloc((unsigned long)len * sizeof(cf16_t));
                memcpy(dst->cf_16[i], src->cf_16[s[COEFFS]],
                        (unsigned long)len * sizeof(cf16_t));
                sort_terms_ff_16(&(dst->cf_16[i]), &(dst->hm[i]), dht);
                break;
            case 32:
                dst->cf_32[i] = (cf32_t *)malloc((unsigned long)len * sizeof(cf32_t));
                memcpy(dst->cf_32[i], src->cf_32[s[COEFFS]],
                        (unsigned long)len * sizeof(cf32_t));
                sort_terms_ff_32(&(dst->cf_32[i]), &(dst->hm[i]), dht);
                break;
            case 64:
                dst->cf_64[i] = (cf64_t *)malloc((unsigned long)len * sizeof(cf64_t));
                memcpy(dst->cf_64[i], src->cf_64[s[COEFFS]],
                        (unsigned long)len * sizeof(cf64_t));
                sort_terms_ff_64(&(dst->cf_64[i]), &(dst->hm[i]), dht);
                break;
            default:
                exit(1);
        }
        /* degree as set when importing input data, see import_input_data */
        hm_t *h = dst->hm[i];
        deg_t deg = dht->hd[h[OFFSET]].deg;
        if (!to_drl) {
            for (j = OFFSET+1; j < len+OFFSET; ++j) {
                if (deg < dht->hd[h[j]].deg) {
                    deg = dht->hd[h[j]].deg;
                }
            }
        }
        h[DEG] = deg;
    }
    dst->ld = n;

    free(e);
    return 1;
}

/* Do we compute a DRL Groebner basis first and the Groebner basis w.r.t.
 * the elimination order from it? For an elimination order F4 tends to
 * reach very high degrees: the tails of the polynomials have a much
 * higher degree than their leading terms, and F4 runs through long
 * chains of S-pairs, e.g. m*x^k, m*x^(k-1), ..., m for a variable x of the
 * second block. Starting from the DRL Groebner basis avoids most of
 * them. This pays off for 0-dimensional ideals, for positive-dimensional
 * ones the DRL Groebner basis is often as expensive as the elimination.
 * By Krull's principal ideal theorem the ideal is positive-dimensional
 * if there are less generators than variables. If only one variable is
 * eliminated the elimination order is close to DRL and the direct
 * computation is usually faster. */
static int use_drl_for_elimination(
        const md_t * const md
        )
{
    if (md->nev == 0 || md->elim_drl == 0 || md->boolean
            || md->use_signatures > 0) {
        return 0;
    }
    if (md->elim_drl == 1) {
        return 1;
    }
    return md->nev > 1 && md->ngens >= md->nvars;
}

/* meta data for the DRL Groebner basis computation from the global meta
 * data gmd of the elimination order */
static md_t *initialize_drl_meta_data(
        const md_t * const gmd
        )
{
    md_t *md = (md_t *)malloc(sizeof(md_t));
    memcpy(md, gmd, sizeof(md_t));

    md->nev         = 0;
    md->elim_drl    = 0;
    md->round_deg   = 0;
    md->reduce_gb   = 1;
    md->drl_md      = NULL;
    md->drl_bs      = NULL;
    md->tr          = NULL;
    md->trace_level = NO_TRACER;
    md->trace_rd    = 0;
    md->tc_state    = 0;
    md->tcoff       = 0;
    md->tcok        = NULL;
    md->ht          = NULL;
    md->hcm         = NULL;
    md->ps          = NULL;
    md->elim_md     = NULL;
    md->elim_bs     = NULL;

    return md;
}

/* A basis hash table for the terms of the elements bs->hm[idx[0]], ...,
 * bs->hm[idx[n-1]], at most of the initial size given by md. Over the
 * rationals it is copied for each prime, so it should not be much larger
 * than needed, it is enlarged during the computation if necessary. */
static ht_t *initialize_small_basis_hash_table(
        md_t *md,
        const bs_t * const bs,
        const bl_t * const idx,
        const len_t n
        )
{
    len_t i;
    int64_t nterms = 0;

    for (i = 0; i < n; ++i) {
        nterms += bs->hm[idx[i]][LENGTH];
    }
    const int32_t hts = md->init_hts;
    int32_t sz = 10;
    while (sz < hts && ((int64_t)1 << sz) < 8 * nterms) {
        sz++;
    }
    md->init_hts = sz < hts ? sz : hts;
    ht_t *ht = initialize_basis_hash_table(md);
    md->init_hts = hts;

    return ht;
}

/* the input generators gbs w.r.t. DRL, using the new hash table of
 * the basis */
static bs_t *initialize_drl_basis(
        const bs_t * const gbs,
        md_t *md,
        const len_t nev
        )
{
    len_t i;

    bl_t *idx = (bl_t *)malloc((unsigned long)gbs->ld * sizeof(bl_t));
    for (i = 0; i < gbs->ld; ++i) {
        idx[i] = i;
    }
    bs_t *bs = initialize_basis(md,
            initialize_small_basis_hash_table(md, gbs, idx, gbs->ld));
    copy_basis_to_other_order(bs, gbs, idx, gbs->ld, nev, md, 0);
    free(idx);

    calculate_divmask(bs->ht);
    sort_r(bs->hm, (unsigned long)bs->ld, sizeof(hm_t *),
            initial_input_cmp, bs->ht);

    return bs;
}

void free_elim_drl_data(
        md_t *md
        )
{
    if (md->drl_bs != NULL) {
        free_basis(&(md->drl_bs));
    }
    if (md->drl_md != NULL) {
        if (md->drl_md->tr != NULL) {
            free_trace(&(md->drl_md->tr));
        }
        free(md->drl_md);
        md->drl_md = NULL;
    }
    if (md->elim_bs != NULL) {
        /* its hash table is the one of the input */
        free_basis_without_hash_table(&(md->elim_bs));
    }
    if (md->elim_md != NULL) {
        if (md->elim_md->tr != NULL) {
            free_trace(&(md->elim_md->tr));
        }
        free(md->elim_md);
        md->elim_md = NULL;
    }
}

/* do the elements of a and b have the same terms? */
static int same_supports(
        const bs_t * const a,
        const bs_t * const b
        )
{
    len_t i;

    if (a->ld != b->ld) {
        return 0;
    }
    for (i = 0; i < a->ld; ++i) {
        const len_t len = a->hm[i][LENGTH];
        if (len != b->hm[i][LENGTH]
                || memcmp(a->hm[i]+OFFSET, b->hm[i]+OFFSET,
                    (unsigned long)len * sizeof(hm_t)) != 0) {
            return 0;
        }
    }
    return 1;
}

/* Computes the Groebner basis of gbs w.r.t. the elimination order of
 * gmd via the DRL Groebner basis. Over the rationals each prime fc
 * computes the DRL Groebner basis using the tracer learnt with the first
 * prime, the DRL meta data and the input generators w.r.t. DRL are kept
 * in gmd for the next primes, see free_elim_drl_data. */
static bs_t *core_f4_elim_via_drl(
        bs_t *gbs,
        md_t *gmd,
        int32_t *errp,
        const uint64_t fc
        )
{
    len_t i;

    const int modular = gmd->fc != fc;

    if (gmd->drl_md == NULL) {
        gmd->drl_md = initialize_drl_meta_data(gmd);
        gmd->drl_bs = initialize_drl_basis(gbs, gmd->drl_md, gmd->nev);
    }
    md_t *dmd = gmd->drl_md;
    /* these may be changed by the multi-modular computation */
    dmd->f4_qq_round = gmd->f4_qq_round;
    dmd->info_level  = gmd->info_level;
    dmd->nthrds      = gmd->nthrds;

    if (gmd->info_level > 0) {
        fprintf(VERBSTREAM, "\nGroebner basis w.r.t. DRL first, ");
        fprintf(VERBSTREAM, "elimination order from it.\n");
    }
    /* DRL Groebner basis */
    bs_t *dbs = run_f4(gmd->drl_bs, dmd, errp, fc);
    if (!modular) {
        /* run_f4 has computed in place */
        gmd->drl_bs = NULL;
    }
    if (dbs == NULL || *errp > 0) {
        return NULL;
    }

    /* its elements are the input w.r.t. the elimination order */
    md_t *emd = NULL;
    if (modular) {
        /* over the rationals the elimination step is traced with the
         * first prime and the tracer is applied for the next ones, so
         * the input is given with integer coefficients */
        if (gmd->elim_md == NULL) {
            gmd->elim_md = initialize_drl_meta_data(gmd);
            gmd->elim_md->nev       = gmd->nev;
            gmd->elim_md->round_deg = 1;
            gmd->elim_md->reduce_gb = gmd->reduce_gb;
        }
        emd = gmd->elim_md;
        emd->f4_qq_round = gmd->f4_qq_round;
        emd->info_level  = gmd->info_level;
        emd->nthrds      = gmd->nthrds;
    } else {
        emd = initialize_drl_meta_data(gmd);
        emd->nev       = gmd->nev;
        emd->round_deg = 1;
        emd->reduce_gb = gmd->reduce_gb;
        emd->fc        = fc;
        set_ff_bits(emd, fc);
    }

    /* over the rationals the terms are already in the hash table of
     * gbs after the first prime if the prime is good, they are only
     * looked up, so this hash table is not changed when applying the
     * tracer. It cannot be replaced by a smaller one: the hash table of
     * the basis of the first prime is used by the caller. */
    const int lookup = modular && gmd->elim_bs != NULL;
    bs_t *ebs = initialize_basis(emd, gbs->ht);
    int ok = copy_basis_to_other_order(ebs, dbs, dbs->lmps, dbs->lml,
            gmd->nev, emd, lookup);
    if (ok) {
        sort_r(ebs->hm, (unsigned long)ebs->ld, sizeof(hm_t *),
                initial_input_cmp, ebs->ht);
        for (i = 0; i < ebs->ld; ++i) {
            ebs->lmps[i] = i;
        }
    }

    if (!modular) {
        free_basis(&dbs);
        free_elim_drl_data(gmd);
    } else {
        if (gmd->f4_qq_round != 1) {
            free_basis_and_only_local_hash_table_data(&dbs);
        } else {
            free_basis_without_hash_table(&dbs);
        }
    }

    if (!modular) {
        emd->ngens = ebs->ld;
        /* Groebner basis w.r.t. the elimination order */
        bs_t *bs = run_f4(ebs, emd, errp, fc);
        free(emd);

        if (bs != NULL) {
            /* as run_f4 the basis is returned instead of the input */
            free_basis_without_hash_table(&gbs);
        }
        return bs;
    }

    if (gmd->elim_bs == NULL) {
        emd->ngens = ebs->ld;
    } else if (!ok || !same_supports(ebs, gmd->elim_bs)) {
        /* the DRL Groebner basis differs from the one of the first
         * prime, so the tracer cannot be applied: bad prime */
        free_basis_without_hash_table(&ebs);
        *errp = 1;
        return NULL;
    }

    /* Groebner basis w.r.t. the elimination order, ebs uses the hash
     * table of gbs */
    bs_t *bs = run_f4(ebs, emd, errp, fc);

    if (gmd->elim_bs == NULL) {
        gmd->elim_bs = ebs;
    } else {
        free_basis_without_hash_table(&ebs);
    }
    return bs;
}

bs_t *core_f4(
        bs_t *gbs,
        md_t *gmd,
        int32_t *errp,
        const uint64_t fc
        )
{
    if (use_drl_for_elimination(gmd)) {
        return core_f4_elim_via_drl(gbs, gmd, errp, fc);
    }
    return run_f4(gbs, gmd, errp, fc);
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
 *     For primes < 2^31 the coefficients are int32_t, for primes
 *     2^31 <= p < 2^52 they are int64_t, also for the returned basis.
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *
 *  RETURNs the length of the jl_basis array */
int64_t export_f4_64(
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
        const uint64_t field_char,
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
            reduce_gb, pbm_file, 0 /*truncate_lifting*/, info_level);

    /* all input generators are invalid */
    if (success == -1) {
        return_zero(bld, blen, bexp, bcf, nr_vars, field_char, mallocp);
        return 1;
    }
    if (success == 0) {
        fprintf(ERRSTREAM,"Bad input data, stopped computation.\n");
        exit(1);
    }

    int err = 0;
    bs = core_f4(bs, md, &err, field_char);

    if (err) {
        fprintf(ERRSTREAM,"Problem with F4, stopped computation.\n");
        exit(1);
    }

    int64_t nterms  = export_results_from_f4(bld, blen, bexp,
            bcf, mallocp, &bs, &bht, &md);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->f4_ctime = ct1 - ct0;
    md->f4_rtime = rt1 - rt0;

    get_and_print_final_statistics(ERRSTREAM, md, bs);

    /* free and clean up */
    free_shared_hash_data(bht);
    if (bs != NULL) {
        free_basis(&bs);
    }

    free(md);
    md    = NULL;

    return nterms;
}

/* same as export_f4_64(), but with a 32 bit field characteristic and
 * int32_t coefficients, thus only for primes < 2^31, this is the
 * interface used by AlgebraicSolving.jl */
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
    if (field_char > INT32_MAX) {
        fprintf(ERRSTREAM, "Field characteristic must be 0 or < 2^31.\n");
        exit(1);
    }
    return export_f4_64(mallocp, bld, blen, bexp, bcf, lens, exps, cfs,
            (uint64_t)field_char, mon_order, elim_block_len, nr_vars,
            nr_gens, ht_size, nr_threads, max_nr_pairs, reset_ht,
            la_option, reduce_gb, pbm_file, info_level);
}
