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

static int is_signature_needed(
        const smat_t * const smat,
        const crit_t * const syz,
        const crit_t * const rew,
        const smat_t *psmat,
        const len_t idx,
        const len_t var_idx,
        ht_t *ht
        )
{
    /* get exponent vector and increment entry for var_idx */
    exp_t *ev   =   ht->ev[0];
    ev          =   ht->ev[psmat->cols[idx][SM_SMON]];
    len_t shift =   var_idx > ht->ebl ? 2 : 1;
    ev[var_idx+shift]++;

    const len_t sig_idx =   psmat->cols[idx][SM_SIDX];
    const hm_t hm       =   insert_in_hash_table(ev, ht);
    const sdm_t nsdm    =   ~ht->hd[hm].sdm;

    const len_t evl     =   ht->evl;

    /* syzygy criterion */
    const crit_t syz_idx =   syz[sig_idx];
    for (len_t i = 0; i < syz_idx.ld; ++i) {
        if (nsdm & syz_idx.sdm[i]) {
            continue;
        }
        const exp_t *sev    =   ht->ev[syz_idx.hm[i]];
        for (len_t j = 0; j < evl; ++j) {
            if (sev[j] > ev[j]) {
                continue;
            }
        }
        return 1;
    }
    /* rewrite criterion */
    const crit_t rew_idx =   rew[sig_idx];
    for (len_t i = 0; i < rew_idx.ld; ++i) {
        if (nsdm & rew_idx.sdm[i]) {
            continue;
        }
        const exp_t *sev    =   ht->ev[rew_idx.hm[i]];
        for (len_t j = 0; j < evl; ++j) {
            if (sev[j] > ev[j]) {
                continue;
            }
        }
        return 1;
    }
    return 0;
}


static void add_multiples_of_previous_degree_row(
        smat_t *smat,
        smat_t *psmat,
        const len_t idx,
        const crit_t * const syz,
        const crit_t * const rew,
        ht_t *ht,
        stat_t *st)
{
    const len_t nv  =   ht->nv;

    for (len_t i = 0; i < nv; ++i) {
        /* check syzygy and rewrite criterion */
        if (is_signature_needed(smat, syz, rew, psmat, idx, i, ht) == 1) {
        }
    }
}

static inline void add_row_with_signature(
        smat_t *mat,
        const bs_t * const bs,
        const len_t pos
        )
{
    const len_t nr          =   mat->nr;
    const unsigned long len =   bs->hm[pos][LENGTH];
    mat->cols[nr]           =   (hm_t *)malloc(
            (len + SM_OFFSET) * sizeof(hm_t));
    /* copy polynomial data */
    memcpy(mat->cols[nr]+SM_CFS,bs->hm[pos]+COEFFS,
            len + OFFSET - COEFFS);
    mat->prev_cf32[nr]      =   bs->cf_32[bs->hm[pos][COEFFS]];
    mat->cols[nr][SM_CFS]   =   nr;
    /* store also signature data */
    mat->cols[nr][SM_SMON]  =   bs->sm[pos];
    mat->cols[nr][SM_SIDX]  =   bs->si[pos];
    mat->nr++;
}

static inline crit_t *initialize_signature_criteria(
        const stat_t * const st
        )
{
    crit_t *crit    =   calloc((unsigned long)st->ngens, sizeof(crit_t));

    return crit;
}

static inline void free_signature_criteria_data(
        crit_t *crit,
        const stat_t * const st
        )
{
    for (len_t i = 0; i < st->ngens; ++i) {
        free(crit[i].sdm);
        crit[i].sdm =   NULL;
        free(crit[i].hm);
        crit[i].hm  =   NULL;
        crit[i].ld  =   0;
    }
}

static inline void free_signature_criteria(
        crit_t **critp,
        const stat_t * const st
        )
{
    crit_t *crit    =   *critp;
    free_signature_criteria_data(crit, st);
    free(crit);
    crit    =   NULL;

    *critp  =   crit;
}

static inline crit_t *initialize_syzygies_schreyer(
        const bs_t * const bs,
        ht_t *ht
        )
{
    /* when initializing syzygies we assume that bs->ld == st->ngens */
    crit_t *syz =   calloc((unsigned long)bs->ld, sizeof(crit_t));
    syz[0].ld   =   0;
    syz[0].sz   =   0;
    for (len_t i = 1; i < bs->ld; ++i) {
        syz[i].hm   =   calloc((unsigned long)i, sizeof(hm_t));
        syz[i].sdm  =   calloc((unsigned long)i, sizeof(sdm_t));
        syz[i].ld   =   0;
        syz[i].sz   =   i;
        for (len_t j = 0; j < i; ++j) {
            syz[i].hm[j]    = insert_multiplied_signature_in_hash_table(
                    bs->hm[j][OFFSET], bs->sm[i], ht);
            syz[i].sdm[j]   = ht->hd[syz[i].hm[j]].sdm;
        }
    }
    return syz;
}

static inline void initialize_signatures_schreyer(
        bs_t *bs
        )
{
    for (len_t i = 0; i < bs->ld; ++i) {
        bs->si[i]   =   i;
        bs->sm[i]   =   bs->hm[i][OFFSET];
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
        ht_t **htp,
        stat_t **stp
        )
{
    bs_t *in    = *bsp;
    ht_t *ht    = *htp;
    stat_t *st  = *stp;

    /* timings for one round */
    double rrt0, rrt1;

    int try_termination =   0;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hi_t *hcm = (hi_t *)malloc(sizeof(hi_t));

    /* signature matrix and previous degree signature matrix */
    smat_t *smat = NULL, *psmat = NULL;

    /* initialize signature related information */
    initialize_signatures_schreyer(in);
    crit_t *syz =   initialize_syzygies_schreyer(in, ht);
    crit_t *rew =   initialize_signature_criteria(st);

    /* initialize an empty basis for keeping the real basis elements */
    bs_t *bs    =   initialize_basis(st);

    /* sort initial elements, highest lead term first */
    sort_r(in->hm, (unsigned long)in->ld, sizeof(hm_t *),
            initial_input_cmp_sig, ht);

    int32_t next_degree =   in->hm[in->ld-1][DEG];
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
                new data             time(rd)\n");
        printf("-------------------------------------------------\
                ----------------------------------------\n");
    }
    st->current_rd  =   0;

    while (!try_termination) {
        rrt0  = realtime();
        st->max_bht_size    = st->max_bht_size > ht->esz ?
            st->max_bht_size : ht->esz;
        st->current_rd++;

        /* prepare signature matrix for next degree */
        psmat   =   smat;
        smat    =   (smat_t *)calloc(1, sizeof(smat_t));

        smat->sz        =   psmat->nr * st->nvars + in->ld;
        smat->cols      =   (hm_t **)calloc(
                (unsigned long)smat->sz,sizeof(hm_t *));
        smat->prev_cf32 =   (cf32_t **)calloc(
                (unsigned long)psmat->nr,sizeof(cf32_t *));

        /* check if we have initial generators not handled in lower degree
         * until now */
        while (in->ld > 0 && in->hm[in->ld-1][DEG] == next_degree) {
            add_row_with_signature(smat, in, in->ld-1);
            in->ld--;
        }
        /* generate rows from previous degree matrix */
        for (len_t i = 0; i < psmat->nr; ++i) {
            add_multiples_of_previous_degree_row(
                    smat, psmat, i, syz, rew, ht, st);
        }

        /* map hashes to columns */

        /* signature-reduce matrix */

        /* flag rows with new leading terms, i.e. non-multiples of
         * already known leading terms */

        /* fully reduce elements with new leading terms */

        /* maps columns to hashes */

        /* add new elements to basis */




        /* if we found a constant we are done, so remove all remaining pairs */
        if (bs->constant  == 1) {
            try_termination =   1;
        }
        rrt1 = realtime();
        if (st->info_level > 1) {
            printf("%13.2f sec\n", rrt1-rrt0);
        }

        /* TODO: termination check, like no new elements the last rounds, etc. */
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
                ----------------------------------------\n");
    }
    if (st->nev > 0) {
        len_t j = 0;
        for (len_t i = 0; i < bs->lml; ++i) {
            if (ht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
                bs->lm[j]   = bs->lm[i];
                bs->lmps[j] = bs->lmps[i];
                ++j;
            }
        }
        bs->lml = j;
    }

    *bsp    = bs;
    *htp    = ht;
    *stp    = st;

    /* free and clean up */
    free(hcm);

    return 1;
}

