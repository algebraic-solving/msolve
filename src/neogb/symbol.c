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

#include "../msolve/streams.h"

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
        fprintf(VERBSTREAM, "%3d  %6d %7d", 0, nps, psl->ld);
        fflush(VERBSTREAM);
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
    ps_t *psl   = md->ps;
    ht_t *sht   = md->ht;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps     = psl->p;
    const len_t evl = bht->evl;

    /* get minimal degree */
    mdeg  = ps[0].deg;
    for (i = 1; i < psl->ld; ++i) {
        if (ps[i].deg < mdeg) {
            mdeg  = ps[i].deg;
        }
    }
    /* only pairs of minimal degree are selected, so we move them to the
     * front and sort only those, the degree is the first criterion of
     * all spair orders */
    for (i = 0, j = 0; i < psl->ld; ++i) {
        if (ps[i].deg == mdeg) {
            spair_t tmp = ps[i];
            ps[i]       = ps[j];
            ps[j++]     = tmp;
        }
    }
    sort_r(ps, (unsigned long)j, sizeof(spair_t), spair_cmp, bht);

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
        fprintf(VERBSTREAM, "%3d  %6d %7d", mdeg, nps, psl->ld);
        fflush(VERBSTREAM);
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

    /* we first collect all rows, i.e. the generators and their multipliers,
     * and then generate all rows at once, see
     * multiplied_polys_to_matrix_rows() */
    const len_t mnr = 2 * nps;
    hm_t **rows   = (hm_t **)malloc((unsigned long)mnr * sizeof(hm_t *));
    hm_t **polys  = (hm_t **)malloc((unsigned long)mnr * sizeof(hm_t *));
    val_t *hs     = (val_t *)malloc((unsigned long)mnr * sizeof(val_t));
    exp_t *ems    = (exp_t *)malloc((unsigned long)mnr * evl * sizeof(exp_t));
    len_t *bis    = (len_t *)malloc((unsigned long)mnr * sizeof(len_t));
    hi_t *mhs     = (hi_t *)malloc((unsigned long)mnr * sizeof(hi_t));
    /* is the row the first one of its lcm, i.e. a reducer row? */
    int8_t *isrr  = (int8_t *)malloc((unsigned long)mnr * sizeof(int8_t));
    len_t nr      = 0;

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

        /* first element with given lcm goes into reducer part of matrix,
         * all remaining ones go to to be reduced part */
        for (k = 0; k < load; ++k) {
            /* check sorted list for doubles */
            if (k > 0 && gens[k] == gens[k-1]) {
                continue;
            }
            /* ev might change when enlarging the hash table during
             * insertion of the multiplier when tracing */
            elcm  = bht->ev[lcm];
            b     = bs->hm[gens[k]];
            eb    = bht->ev[b[OFFSET]];
            exp_t *em = ems + (unsigned long)nr * evl;
            for (l = 0; l < evl; ++l) {
                em[l]   =   (exp_t)(elcm[l] - eb[l]);
            }
            hs[nr]    = bht->hd[lcm].val - bht->hd[b[OFFSET]].val;
            polys[nr] = b;
            bis[nr]   = gens[k];
            isrr[nr]  = (k == 0);
            /* track trace information ? */
            if (md->trace_level == LEARN_TRACER) {
                if (bs->ht->eld == bs->ht->esz-1) {
                    enlarge_hash_table(bs->ht);
                }
#if PARALLEL_HASHING
                mhs[nr] = check_insert_in_hash_table(em, hs[nr], bs->ht);
#else
                mhs[nr] = insert_in_hash_table(em, bs->ht);
#endif
            }
            nr++;
        }
        i = j;
    }

    multiplied_polys_to_matrix_rows(rows, nr, hs, ems, polys, sht, bht,
            md->nthrds);

    for (k = 0; k < nr; ++k) {
        if (md->trace_level == LEARN_TRACER) {
            rows[k][BINDEX] = bis[k];
            rows[k][MULT]   = mhs[k];
        }
        /* mark lcm column as lead term column */
        sht->hd[rows[k][OFFSET]].idx = 2;
        if (isrr[k]) {
            rrows[nrr++]  = rows[k];
        } else {
            trows[ntr++]  = rows[k];
        }
        mat->nr++;
    }
    free(rows);
    free(polys);
    free(hs);
    free(ems);
    free(bis);
    free(mhs);
    free(isrr);
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
        val_t h = 0;
        for (len_t j = 0; j < bht->evl; ++j) {
            h +=  bht->rn[j] * mul[j];
        }
        trows[ntr++]    = multiplied_poly_to_matrix_row(
                sht, bht, h, mul, b);
        mat->nr++;
    }
}


/* returns the position of the reducer of the monomial m in bs->lmps
 * plus one, 0 if m is not divisible by any lead term of bs. Among all
 * divisors we choose the one of shortest length, the newest one among the
 * shortest ones: this leads to sparser matrices.
 * Does not change any data, so it can be called in parallel. */
static inline len_t find_reducer(
        const bs_t * const bs,
        const hm_t m,
        const ht_t * const sht
        )
{
    len_t i, k;

    const ht_t * const bht  = bs->ht;
    const len_t evl         = bht->evl;
    const exp_t * const e   = sht->ev[m];
    const sdm_t ns          = ~sht->hd[m].sdm;
    const len_t lml         = bs->lml;

    const sdm_t * const lms   = bs->lm;
    const bl_t * const lmps   = bs->lmps;
    exp_t * const * const evb = bht->ev;

    len_t bi    = 0;
    len_t blen  = 0;
    for (i = 0; i < lml; ++i) {
        if (lms[i] & ns) {
            continue;
        }
        const hm_t *b = bs->hm[lmps[i]];
        if (blen != 0 && b[LENGTH] > blen) {
            continue;
        }
        const exp_t * const f = evb[b[OFFSET]];
        for (k=0; k < evl; ++k) {
            if (e[k] < f[k]) {
                break;
            }
        }
        if (k == evl && (blen == 0 || b[LENGTH] < blen || lmps[i] > lmps[bi])) {
            bi    = i;
            blen  = b[LENGTH];
        }
    }
    return blen > 0 ? bi + 1 : 0;
}

static void symbolic_preprocessing(
        mat_t *mat,
        bs_t *bs,
        md_t *md
        )
{
    hl_t i;
    len_t j, k;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* at the moment we have as many reducers as we have different lcms */
    len_t nrr = mat->nc;

    ht_t *sht = md->ht;
    ht_t *bht = bs->ht;

    const len_t evl   = bht->evl;
    const len_t onrr  = mat->nc;
    /* for the columns already present we have to check whether they are
     * already handled (e.g. lcms of the selected pairs), the columns added
     * here are always new */
    const hl_t oesld  = sht->eld;
    while (mat->sz <= nrr + oesld) {
        mat->sz *=  2;
        mat->rr =   realloc(mat->rr, (unsigned long)mat->sz * sizeof(hm_t *));
        if (mat->rr == NULL) {
            fprintf(ERRSTREAM, "Allocating memory for matrix failed,\n");
            fprintf(ERRSTREAM, "segmentation fault will follow.\n");
        }
    }

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following.
     *
     * we handle the columns in waves: all columns not yet handled are
     * collected, their reducers are searched in parallel, then the
     * reducer rows are generated, adding new columns for the next wave.
     * the rows and new columns are generated in the same order as when
     * handling one column after the other. */
    len_t sz    = 0;
    hm_t *cols  = NULL;
    len_t *red  = NULL;
    hm_t **polys  = NULL;
    val_t *hs     = NULL;
    exp_t *ems    = NULL;

    hl_t lo = 1;
    while (lo < sht->eld) {
        const hl_t hi = sht->eld;
        if (sz < hi - lo) {
            sz    = hi - lo;
            cols  = realloc(cols, (unsigned long)sz * sizeof(hm_t));
            red   = realloc(red, (unsigned long)sz * sizeof(len_t));
            polys = realloc(polys, (unsigned long)sz * sizeof(hm_t *));
            hs    = realloc(hs, (unsigned long)sz * sizeof(val_t));
            ems   = realloc(ems, (unsigned long)sz * evl * sizeof(exp_t));
        }
        len_t nc = 0;
        for (i = lo; i < hi; ++i) {
            if (i >= oesld || !sht->hd[i].idx) {
                sht->hd[i].idx = 1;
                mat->nc++;
                cols[nc++] = i;
            }
        }
#pragma omp parallel for num_threads(md->nthrds) \
        private(j) schedule(dynamic, 32) if (md->nthrds > 1 && nc > 64)
        for (j = 0; j < nc; ++j) {
            red[j]  = find_reducer(bs, cols[j], sht);
        }
        len_t nr = 0;
        for (j = 0; j < nc; ++j) {
            if (red[j] == 0) {
                continue;
            }
            const hm_t m        = cols[j];
            const bl_t bi       = bs->lmps[red[j]-1];
            const hm_t *b       = bs->hm[bi];
            const exp_t *e      = sht->ev[m];
            const exp_t *f      = bht->ev[b[OFFSET]];
            exp_t *em           = ems + (unsigned long)nr * evl;
            for (k = 0; k < evl; ++k) {
                em[k] = (exp_t)(e[k]-f[k]);
            }
            hs[nr]    = sht->hd[m].val - bht->hd[b[OFFSET]].val;
            polys[nr] = (hm_t *)b;
            /* we need the basis index for tracing */
            red[nr]   = bi;
            sht->hd[m].idx  = 2;
            nr++;
        }
        while (mat->sz <= nrr + nr) {
            mat->sz *=  2;
            mat->rr =   realloc(mat->rr, (unsigned long)mat->sz * sizeof(hm_t *));
            if (mat->rr == NULL) {
                fprintf(ERRSTREAM, "Allocating memory for matrix failed,\n");
                fprintf(ERRSTREAM, "segmentation fault will follow.\n");
            }
        }
        /* track trace information ? */
        hi_t *mhs = NULL;
        if (md->trace_level == LEARN_TRACER) {
            mhs = (hi_t *)malloc((unsigned long)nr * sizeof(hi_t));
            for (j = 0; j < nr; ++j) {
                if (bht->eld == bht->esz-1) {
                    enlarge_hash_table(bht);
                }
#if PARALLEL_HASHING
                mhs[j]  = check_insert_in_hash_table(
                        ems + (unsigned long)j * evl, hs[j], bht);
#else
                mhs[j]  = insert_in_hash_table(
                        ems + (unsigned long)j * evl, bht);
#endif
            }
        }
        multiplied_polys_to_matrix_rows(mat->rr + nrr, nr, hs, ems, polys,
                sht, bht, md->nthrds);
        if (md->trace_level == LEARN_TRACER) {
            for (j = 0; j < nr; ++j) {
                mat->rr[nrr+j][BINDEX]  = red[j];
                mat->rr[nrr+j][MULT]    = mhs[j];
            }
            free(mhs);
        }
        nrr +=  nr;
        lo  =   hi;
    }
    free(cols);
    free(red);
    free(polys);
    free(hs);
    free(ems);
    /* realloc to real size */
    mat->rr   =   realloc(mat->rr, (unsigned long)nrr * sizeof(hm_t *));
    mat->nr   +=  nrr - onrr;
    mat->nrl  =   mat->nr - nrr;
    mat->nru  =   nrr;
    mat->sz   =   mat->nr;
    mat->rbal =   mat->nrl;

    /* initialize memory for reducer bit arrays for tracing information */
    mat->rba  = (rba_t **)malloc((uint64_t)mat->rbal * sizeof(rba_t *));
    const uint64_t len = nrr / 32 + ((nrr % 32) != 0);
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
    if (mat->nc == 0) {
        return 1;
    } else {
        return 0;
    }
}
