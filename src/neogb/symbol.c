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


#include "data.h"

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

/* select_all_pairs() is unused at the moment */
#if 0 
static void select_all_spairs(
        mat_t *mat,
        const bs_t * const bs,
        ps_t *psl,
        md_t *st,
        ht_t *sht,
        ht_t *bht,
        ht_t *tht
        )
{
    len_t i, j, k, l, nps, npd, nrr = 0, ntr = 0;
    hm_t *b;
    len_t load  = 0;
    hi_t lcm;
    len_t *gens;
    exp_t *elcm, *eb;
    exp_t etmp[bht->evl];

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps     = psl->p;
    const len_t nv  = bht->nv;

    /* sort pair set */
    sort_r(ps, (unsigned long)psl->ld, sizeof(spair_t), spair_degree_cmp, bht);

    /* select pairs of this degree respecting maximal selection size mnsel */
    npd  = psl->ld;
    sort_r(ps, (unsigned long)npd, sizeof(spair_t), spair_cmp, bht);
    /* now do maximal selection if it applies */

    nps = psl->ld;
    
    if (st->info_level > 1) {
        printf("%3d  %6d %7d", 0, nps, psl->ld);
        fflush(stdout);
    }
    /* statistics */
    st->num_pairsred  +=  nps;
    /* list for generators */
    gens  = (len_t *)malloc(2 * (unsigned long)nps * sizeof(len_t));
    /* preset matrix meta data */
    mat->rr       = (hm_t **)malloc(2 * (unsigned long)nps * sizeof(hm_t *));
    hm_t **rrows  = mat->rr;
    mat->tr       = (hm_t **)malloc(2 * (unsigned long)nps * sizeof(hm_t *));
    hm_t **trows  = mat->tr;
    mat->sz = 2 * nps;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = 0;

    int ctr = 0;


    i = 0;

    while (i < nps) {
        /* ncols initially counts number of different lcms */
        mat->nc++;
        load  = 0;
        lcm   = ps[i].lcm;
        j = i;

        while (j < nps && ps[j].lcm == lcm) {
            gens[load++] = ps[j].gen1;
            gens[load++] = ps[j].gen2;
            ++j;
        }
        /* sort gens set */
        qsort(gens, (unsigned long)load, sizeof(len_t), gens_cmp);

        len_t prev  = -1;

        /* first element with given lcm goes into reducer part of matrix,
         * all remaining ones go to to be reduced part */
        prev  = gens[0];
        /* printf("prev %u / %u\n", prev, bs->ld); */
        /* ev might change when enlarging the hash table during insertion of a new
            * row in the matrix, thus we have to reset elcm inside the for loop */
        elcm  = bht->ev[lcm];
        b     = bs->hm[prev];
        eb    = bht->ev[b[OFFSET]];
        for (l = 0; l <= nv; ++l) {
            etmp[l]   =   (exp_t)(elcm[l] - eb[l]);
        }
        const hi_t h    = bht->hd[lcm].val - bht->hd[b[OFFSET]].val;
        /* note that we use index mat->nc and not mat->nr since for each new
         * lcm we add exactly one row to mat->rr */
        rrows[nrr]  = multiplied_poly_to_matrix_row(sht, bht, h, etmp, b);
        /* track trace information ? */
        if (tht != NULL) { 
           rrows[nrr][BINDEX]  = prev;
            if (tht->eld == tht->esz-1) {
                enlarge_hash_table(tht);
            }
            rrows[nrr][MULT]    = insert_in_hash_table(etmp, tht);
        }

        /* mark lcm column as lead term column */
        sht->hd[rrows[nrr++][OFFSET]].idx = 2; 
        /* still we have to increase the number of rows */
        mat->nr++;
        for (k = 1; k < load; ++k) {
            /* check sorted list for doubles */
            if (gens[k] ==  prev) {
                continue;
            }
            prev  = gens[k];
            /* ev might change when enlarging the hash table during insertion of a new
             * row in the matrix, thus we have to reset elcm inside the for loop */
            elcm  = bht->ev[lcm];
            if (elcm[0] > 0) {
                /* printf("pair with lcm ");
                 * for (int ii = 0; ii < nv; ++ii) {
                 *     printf("%u ", elcm[ii]);
                 * }
                 * printf("\n"); */
            }
            b     = bs->hm[prev];
            eb    = bht->ev[b[OFFSET]];
            for (l = 0; l <= nv; ++l) {
                etmp[l]   =   (exp_t)(elcm[l] - eb[l]);
            }
            const hi_t h  = bht->hd[lcm].val - bht->hd[b[OFFSET]].val;
            trows[ntr] = multiplied_poly_to_matrix_row(sht, bht, h, etmp, b);
            /* track trace information ? */
            if (tht != NULL) {
                trows[ntr][BINDEX]  = prev;
                if (tht->eld == tht->esz-1) {
                    enlarge_hash_table(tht);
                }
                trows[ntr][MULT]    = insert_in_hash_table(etmp, tht);
            }
            /* mark lcm column as lead term column */
            sht->hd[trows[ntr++][OFFSET]].idx = 2;
            mat->nr++;
        }
        ctr++;
        i = j;
    }
    /* printf("nc %u | nr %u || %u\n", mat->nc, mat->nr, sht->eld); */
    /* printf("%u pairs in degree %u\n", ctr, md); */
    /* fix rows to be reduced */
    mat->tr = realloc(mat->tr, (unsigned long)(mat->nr - mat->nc) * sizeof(hm_t *));

    st->num_rowsred +=  mat->nr - mat->nc;
    st->current_deg =   etmp[DEG];

    free(gens);

    /* remove selected spairs from pairset */
    memmove(ps, ps+nps, (unsigned long)(psl->ld-nps) * sizeof(spair_t));
    psl->ld -=  nps;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->select_ctime  +=  ct1 - ct0;
    st->select_rtime  +=  rt1 - rt0;
}
#endif

/* selection of spairs, at the moment only selection
by minial degree of the spairs is supported

NOTE: The pair list is already sorted! */
static int32_t select_spairs_by_minimal_degree(
        mat_t *mat,
        bs_t *bs,
        md_t *md
        )
{
    len_t i, j, k, l, mdeg, nps, npd, nrr = 0, ntr = 0;
    hm_t *b;
    len_t load = 0;
    hi_t lcm;
    len_t *gens;
    exp_t *elcm, *eb;
    ht_t *bht   = bs->ht;
    exp_t etmp[bht->evl];
    ps_t *psl   = md->ps;
    ht_t *sht   = md->ht;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps     = psl->p;
    const len_t evl = bht->evl;

    /* sort pair set */
    sort_r(ps, (unsigned long)psl->ld, sizeof(spair_t), spair_cmp, bht);
    /* get minimal degree */
    mdeg  = ps[0].deg;

    /* compute a truncated GB? Check maximal degree. */
    if (md->max_gb_degree < mdeg) {
        return 1; 
    }

    /* select pairs of this degree respecting maximal selection size mnsel */
#if 0
    printf("pair set sorted for symbolic preprocessing\n");
    int pctr = 0;
    deg_t degtest = 0;
    for (i = 0; i < psl->ld; ++i) {
        if (degtest != ps[i].deg) {
            printf("%d elements\n", pctr);
            printf("-- degree %d --\n", ps[i].deg);
            degtest = ps[i].deg;
            pctr = 0;
        }
        printf("%d --> deg %d --> [%u,%u]", i, ps[i].deg, ps[i].gen1, ps[i].gen2);
        for (int jj = 0; jj < evl; ++jj) {
            printf("%d ", bht->ev[ps[i].lcm][jj]);
        }
        pctr++;
        printf("\n");
    }
    printf("\n");
#endif
    for (i = 0; i < psl->ld; ++i) {
        if (ps[i].deg > mdeg) {
            break;
        }
    }
    npd  = i;
    /* printf("npd %d\n", npd); */
    /* sort_r(ps, (unsigned long)npd, sizeof(spair_t), spair_cmp, bht); */
    /* now do maximal selection if it applies */
    
    /* if we stopped due to maximal selection size we still get the following
     * pairs of the same lcm in this matrix */
    if (npd > md->mnsel) {
        nps = md->mnsel;
        lcm = ps[nps].lcm;
        while (nps < npd && ps[nps+1].lcm == lcm) {
            nps++;
        }
    } else {
        nps = npd;
    }
    if (md->info_level > 1) {
        printf("%3d  %6d %7d", mdeg, nps, psl->ld);
        fflush(stdout);
    }
    /* statistics */
    md->num_pairsred  +=  nps;
    /* list for generators */
    gens  = (len_t *)malloc(2 * (unsigned long)nps * sizeof(len_t));
    /* preset matrix meta data */
    mat->rr       = (hm_t **)malloc(2 * (unsigned long)nps * sizeof(hm_t *));
    hm_t **rrows  = mat->rr;
    mat->tr       = (hm_t **)malloc(2 * (unsigned long)nps * sizeof(hm_t *));
    hm_t **trows  = mat->tr;
    mat->sz = 2 * nps;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = 0;

    i = 0;

    while (i < nps) {
        /* ncols initially counts number of different lcms */
        mat->nc++;
        load  = 0;
        lcm   = ps[i].lcm;
        j = i;

        while (j < nps && ps[j].lcm == lcm) {
            gens[load++] = ps[j].gen1;
            gens[load++] = ps[j].gen2;
            ++j;
        }
        /* sort gens set */
        qsort(gens, (unsigned long)load, sizeof(len_t), gens_cmp);

        len_t prev  = -1;

        /* first element with given lcm goes into reducer part of matrix,
         * all remaining ones go to to be reduced part */
        prev  = gens[0];
        /* ev might change when enlarging the hash table during insertion of a new
            * row in the matrix, thus we have to reset elcm inside the for loop */
        elcm  = bht->ev[lcm];
        b     = bs->hm[prev];
        eb    = bht->ev[b[OFFSET]];
        for (l = 0; l < evl; ++l) {
            etmp[l]   =   (exp_t)(elcm[l] - eb[l]);
        }
        const hi_t h    = bht->hd[lcm].val - bht->hd[b[OFFSET]].val;
        /* note that we use index mat->nc and not mat->nr since for each new
         * lcm we add exactly one row to mat->rr */
        rrows[nrr]  = multiplied_poly_to_matrix_row(sht, bht, h, etmp, b);
        /* track trace information ? */
        if (md->trace_level == LEARN_TRACER) { 
           rrows[nrr][BINDEX]  = prev;
            if (bs->ht->eld == bs->ht->esz-1) {
                enlarge_hash_table(bs->ht);
            }
#if PARALLEL_HASHING
            rrows[nrr][MULT]    = check_insert_in_hash_table(etmp, h, bs->ht);
#else
            rrows[nrr][MULT]    = insert_in_hash_table(etmp, bs->ht);
#endif
        }

        /* mark lcm column as lead term column */
        sht->hd[rrows[nrr++][OFFSET]].idx = 2; 
        /* still we have to increase the number of rows */
        mat->nr++;
        for (k = 1; k < load; ++k) {
            /* check sorted list for doubles */
            if (gens[k] ==  prev) {
                continue;
            }
            prev  = gens[k];
            /* ev might change when enlarging the hash table during insertion of a new
             * row in the matrix, thus we have to reset elcm inside the for loop */
            elcm  = bht->ev[lcm];
            if (elcm[0] > 0) {
            }
            b     = bs->hm[prev];
            eb    = bht->ev[b[OFFSET]];
            for (l = 0; l < evl; ++l) {
                etmp[l]   =   (exp_t)(elcm[l] - eb[l]);
            }
            const hi_t h  = bht->hd[lcm].val - bht->hd[b[OFFSET]].val;
            trows[ntr] = multiplied_poly_to_matrix_row(sht, bht, h, etmp, b);
            /* track trace information ? */
            if (md->trace_level == LEARN_TRACER) {
                trows[ntr][BINDEX]  = prev;
                if (bs->ht->eld == bs->ht->esz-1) {
                    enlarge_hash_table(bs->ht);
                }
#if PARALLEL_HASHING
                trows[ntr][MULT]    = check_insert_in_hash_table(etmp, h, bs->ht);
#else
                trows[ntr][MULT]    = insert_in_hash_table(etmp, bs->ht);
#endif
            }
            /* mark lcm column as lead term column */
            sht->hd[trows[ntr++][OFFSET]].idx = 2;
            mat->nr++;
        }
        i = j;
    }
    /* printf("%u pairs in degree %u\n", ctr, mdeg); */
    /* fix rows to be reduced */
    mat->tr = realloc(mat->tr, (unsigned long)(mat->nr - mat->nc) * sizeof(hm_t *));

    md->num_rowsred +=  mat->nr - mat->nc;
    md->current_deg =   mdeg;
    mat->cd         =   mdeg;

    free(gens);

    /* remove selected spairs from pairset */
    memmove(ps, ps+nps, (unsigned long)(psl->ld-nps) * sizeof(spair_t));
    psl->ld -=  nps;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->select_ctime  +=  ct1 - ct0;
    md->select_rtime  +=  rt1 - rt0;

    return 0;
}

/* write elements straight to sat, not to a matrix */
static void select_saturation(
        bs_t *sat,
        mat_t *mat,
        md_t *st,
        ht_t *sht,
        ht_t *bht
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();


    /* preset matrix meta data */
    mat->rr = (hm_t **)malloc(100 * sizeof(hm_t *));
    mat->tr = NULL;

    mat->sz = 100;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = 0;

    /* for (i=0; i < sat->hm[0][LENGTH]; ++i) {
     *     printf("%u | ", sat->cf_32[sat->hm[0][COEFFS]][i]);
     *     for (len_t j = 0; j < bht->nv; ++j) {
     *         printf("%u ", bht->ev[sat->hm[0][OFFSET+i]][j]);
     *     }
     *     printf(" ||| ");
     * }
     * printf("\n"); */
    /* move hashes of sat data from bht to sht for linear algebra */
    /* for (i = 0; i < sat->ld; ++i) {
     *     for (j = OFFSET; j < sat->hm[i][LENGTH]+OFFSET; ++j) {
     *         sat->hm[i][j] = insert_in_hash_table(
     *                 bht->ev[sat->hm[i][j]], sht);
     *     }
     * } */

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->select_ctime  +=  ct1 - ct0;
    st->select_rtime  +=  rt1 - rt0;
}


static void select_tbr(
        const bs_t * const tbr,
        const exp_t * const mul,
        const len_t start,
        mat_t *mat,
        md_t *st,
        ht_t *sht,
        ht_t *bht,
        ht_t *tht
        )
{
    len_t i;

    len_t ntr = 0;

    /* preset matrix meta data */
    mat->rr       = (hm_t **)malloc(100 * sizeof(hm_t *));
    mat->tr       = (hm_t **)malloc((unsigned long)tbr->ld * sizeof(hm_t *));
    hm_t **trows  = mat->tr;

    mat->sz = 100;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = 0;

    /* always take all elements in tbr and
     * multiply them by the given multiple */
    for (i = start; i < tbr->ld; ++i) {
        const hm_t *b   = tbr->hm[i];
        /* remove the multiplier business for the moment, no need
         * and it corrupts a bit the sht size for efficient matrix
         * generation */
        /* const hi_t mulh = insert_in_hash_table(mul, sht);
         * const hi_t h    = sht->hd[mulh].val;
         * const deg_t d   = sht->hd[mulh].deg; */
        const hi_t h    = 0;
        trows[ntr++]    = multiplied_poly_to_matrix_row(
                sht, bht, h, mul, b);
        mat->nr++;
    }
}


static inline void find_multiplied_reducer(
        bs_t *bs,
        const hm_t m,
        len_t *nr,
        hm_t **rows,
        ht_t *sht,
        const md_t * const md
        )
{
    len_t i, k;

    ht_t *bht = bs->ht;

    const len_t rr  = *nr;

    const len_t evl = bht->evl;

    const exp_t * const e  = sht->ev[m];

    const hd_t hdm    = sht->hd[m];
    const len_t lml   = bs->lml;
    const sdm_t ns    = ~hdm.sdm;

    const sdm_t * const lms = bs->lm;
    const bl_t * const lmps = bs->lmps;

    exp_t etmp[bht->evl];
    const hd_t * const hdb  = bht->hd;
    exp_t * const * const evb = bht->ev;

    i = 0;
start:
    while (i < lml && lms[i] & ns) {
        i++;
    }
    if (i < lml) {
        const hm_t *b = bs->hm[lmps[i]];
        const exp_t * const f = evb[b[OFFSET]];
        for (k=0; k < evl; ++k) {
            if (e[k] < f[k]) {
                i++;
                goto start;
            }
            etmp[k] = (exp_t)(e[k]-f[k]);
        }

        const hi_t h  = hdm.val - hdb[b[OFFSET]].val;
        rows[rr]  = multiplied_poly_to_matrix_row(sht, bht, h, etmp, b);
        /* track trace information ? */
        if (md->trace_level == LEARN_TRACER) {
            rows[rr][BINDEX]  = lmps[i];
            if (bht->eld == bht->esz-1) {
                enlarge_hash_table(bht);
            }
#if PARALLEL_HASHING
            rows[rr][MULT]    = check_insert_in_hash_table(etmp, h, bht);
#else
            rows[rr][MULT]    = insert_in_hash_table(etmp, bht);
#endif
        }
        sht->hd[m].idx  = 2;
        *nr             = rr + 1;
    }
}

static void symbolic_preprocessing(
        mat_t *mat,
        bs_t *bs,
        md_t *md
        )
{
    hl_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* at the moment we have as many reducers as we have different lcms */
    len_t nrr = mat->nc;

    ht_t *sht = md->ht;

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following. */

    const hl_t oesld = sht->eld;
    const len_t onrr  = mat->nc;
    i = 1;
    /* we only have to check if idx is set for the elements already set
     * when selecting spairs, afterwards (second for loop) we do not
     * have to do this check */
    while (mat->sz <= nrr + oesld) {
        mat->sz *=  2;
        mat->rr =   realloc(mat->rr, (unsigned long)mat->sz * sizeof(hm_t *));
    }
    for (; i < oesld; ++i) {
        if (!sht->hd[i].idx) {
            sht->hd[i].idx = 1;
            mat->nc++;
            find_multiplied_reducer(bs, i, &nrr, mat->rr, sht, md);
        }
    }
    for (; i < sht->eld; ++i) {
        if (mat->sz == nrr) {
            mat->sz *=  2;
            mat->rr  =  realloc(mat->rr, (unsigned long)mat->sz * sizeof(hm_t *));
        }
        sht->hd[i].idx = 1;
        mat->nc++;
        find_multiplied_reducer(bs, i, &nrr, mat->rr, sht, md);
    }
    /* realloc to real size */
    mat->rr   =   realloc(mat->rr, (unsigned long)nrr * sizeof(hm_t *));
    mat->nr   +=  nrr - onrr;
    mat->nrl  =   mat->nr - nrr;
    mat->nru  =   nrr;
    mat->sz   =   mat->nr;
    mat->rbal =   mat->nrl;

    /* initialize memory for reducer bit arrays for tracing information */
    mat->rba  = (rba_t **)malloc((unsigned long)mat->rbal * sizeof(rba_t *));
    const unsigned long len = nrr / 32 + ((nrr % 32) != 0);
    for (i = 0; i < mat->nrl; ++i) {
        mat->rba[i] = (rba_t *)calloc(len, sizeof(rba_t));
    }

    /* statistics */
    md->max_sht_size  = md->max_sht_size > sht->esz ?
        md->max_sht_size : sht->esz;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->symbol_ctime  +=  ct1 - ct0;
    md->symbol_rtime  +=  rt1 - rt0;
}

static void generate_matrix_from_trace(
        mat_t *mat,
        const bs_t * const bs,
        md_t *md
        )
{
    /* timings */
    double ct, rt;
    ct = cputime();
    rt = realtime();

    len_t i, nr;
    hm_t *b;
    exp_t *emul;
    hi_t h;

    const len_t idx = md->trace_rd;

    td_t td   = md->tr->td[idx];
    ht_t *bht = bs->ht;
    ht_t *sht = md->ht;

    mat->rr       = (hm_t **)malloc((unsigned long)td.rld * sizeof(hm_t *));
    hm_t **rrows  = mat->rr;
    mat->tr       = (hm_t **)malloc((unsigned long)td.tld * sizeof(hm_t *));
    hm_t **trows  = mat->tr;
    mat->rba      = (rba_t **)malloc((unsigned long)td.tld * sizeof(rba_t *));
    rba_t **rba   = mat->rba;

    /* reducer rows, i.e. AB part */
    i   = 0;
    nr  = 0;
    while (i < td.rld) {
        b     = bs->hm[td.rri[i++]];
        emul  = bht->ev[td.rri[i]];
        h     = bht->hd[td.rri[i++]].val;


        rrows[nr] = multiplied_poly_to_matrix_row(sht, bht, h, emul, b);
        sht->hd[rrows[nr][OFFSET]].idx = 2;
        ++nr;

    }
    /* to be reduced rows, i.e. CD part */
    i   = 0;
    nr  = 0;
    while (i < td.tld) {
        b     = bs->hm[td.tri[i++]];
        emul  = bht->ev[td.tri[i]];
        h     = bht->hd[td.tri[i]].val;
        trows[nr] = multiplied_poly_to_matrix_row(sht, bht, h, emul, b);
        /* At the moment rba is unused */
        rba[nr]   = td.rba[i/2];
        i++;
        nr++;
    }
    /* meta data for matrix */
    mat->nru  = td.rld/2;
    mat->nrl  = td.tld/2;
    mat->nr   = mat->sz = mat->nru + mat->nrl;
    mat->nc   = sht->eld-1;

    /* statistics */
    md->max_sht_size  = md->max_sht_size > sht->esz ?
        md->max_sht_size : sht->esz;

    /* timings */
    md->tracer_ctime += cputime() - ct;
    md->tracer_rtime += realtime() - rt;

    print_current_trace_meta_data(md);
}

static void generate_saturation_reducer_rows_from_trace(
        mat_t *mat,
        const trace_t * const trace,
        const len_t idx,
        const bs_t * const bs,
        md_t *st,
        ht_t *sht,
        const ht_t * const bht
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, nr;
    hm_t *b;
    exp_t *emul;
    hi_t h;

    ts_t ts       = trace->ts[idx];
    mat->rr       = (hm_t **)malloc((unsigned long)ts.rld * sizeof(hm_t *));
    hm_t **rrows  = mat->rr;

    /* reducer rows, i.e. AB part */
    i   = 0;
    nr  = 0;
    while (i < ts.rld) {
        b     = bs->hm[ts.rri[i++]];
        emul  = bht->ev[ts.rri[i]];
        h     = bht->hd[ts.rri[i++]].val;

        rrows[nr] = multiplied_poly_to_matrix_row(sht, bht, h, emul, b);
        sht->hd[rrows[nr][OFFSET]].idx = 2;
        ++nr;

    }
    /* meta data for matrix */
    mat->nru  = ts.rld/2;
    mat->nr   = mat->sz = mat->nru + mat->nrl;
    mat->nc   = sht->eld-1;

    /* statistics */
    st->max_sht_size  = st->max_sht_size > sht->esz ?
        st->max_sht_size : sht->esz;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->symbol_ctime  +=  ct1 - ct0;
    st->symbol_rtime  +=  rt1 - rt0;
}

static int preprocessing(
        mat_t *mat,
        bs_t *bs,
        md_t *md
        )
{
    if (md->trace_level != APPLY_TRACER) {
        if (select_spairs_by_minimal_degree(mat, bs, md)) {
            return 1;
        }
        symbolic_preprocessing(mat, bs, md);
    } else {
        generate_matrix_from_trace(mat, bs, md);
    }
    return 0;
}
