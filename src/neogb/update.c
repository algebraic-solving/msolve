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


#include "update.h" 

ps_t *initialize_pairset(
        void
        )
{
    ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
    ps->ld  = 0;
    ps->sz  = 192;
    ps->p = (spair_t *)calloc((unsigned long)ps->sz, sizeof(spair_t));
    return ps;
}

static inline void check_enlarge_pairset(
        ps_t *ps,
        len_t added
        )
{
    if (ps->ld+added >= ps->sz) {
        ps->sz  = ps->sz*2 > ps->ld+added ? ps->sz*2 : ps->ld+added;
        ps->p   = realloc(ps->p, (unsigned long)ps->sz * sizeof(spair_t));
        memset(ps->p+ps->ld, 0,
                (unsigned long)(ps->sz-ps->ld) * sizeof(spair_t));
    }
}

void free_pairset(
        ps_t **psp
        )
{
    ps_t *ps  = *psp;
    if (ps->p) {
        free(ps->p);
        ps->p   = NULL;
        ps->ld  = 0;
        ps->sz  = 0;
    }
    free(ps);
    ps  = NULL;
    *psp  = ps;
}

static void insert_and_update_spairs(
        ps_t *psl,
        bs_t *bs,
        ht_t *bht,
        ht_t *uht,
        stat_t *st,
        const int32_t check_redundancy
        )
{
    len_t i, j, l;

    spair_t *ps = psl->p;

#ifdef HAVE_OPENMP
    const int max_nthrds = 4 <= st->nthrds ? 4 : st->nthrds;
#endif

    const len_t pl  = psl->ld;
    const len_t bl  = bs->ld;

    const hm_t nch = bs->hm[bl][OFFSET];

    deg_t ndeg  = bht->hd[nch].deg;
    bs->mltdeg  = bs->mltdeg > ndeg ?
        bs->mltdeg : ndeg;

    reinitialize_hash_table(uht, bl);
    /* statistics */
    st->max_uht_size  = st->max_uht_size > uht->esz ?
        st->max_uht_size : uht->esz;

    /* only other lead terms from the matrix may render
     * the current element useless */
    if (check_redundancy == 1) {
        for (i = bs->lo; i < bl; ++i) {
            if (bs->red[i]) {
                continue;
            }
            if (check_monomial_division(nch, bs->hm[i][OFFSET], bht)) {
                ps[pl].gen1 = i;
                ps[pl].gen2 = bl;
                ps[pl].lcm  = get_lcm(bs->hm[i][OFFSET], nch, bht, bht);
                bs->red[bl] = 1;
                st->num_redundant++;
                bs->ld++;
                psl->ld++;
                return;
            }
        }
    }

    hi_t *plcm  = (hi_t *)malloc((unsigned long)(bl+1) * sizeof(hi_t));
    spair_t *pp = ps+pl;

    /* create all possible new pairs */
    if (check_redundancy == 1) {
        for (i = 0; i < bl; ++i) {
            plcm[i]   = get_lcm(bs->hm[i][OFFSET], nch, bht, uht);
            pp[i].deg = uht->hd[plcm[i]].deg;
            if (bs->red[i] == 0) {
                pp[i].gen1  = i;
                pp[i].gen2  = bl;
                pp[i].lcm   = plcm[i];
            }
        }
    } else {
        for (i = 0; i < bl; ++i) {
            plcm[i]     =  get_lcm(bs->hm[i][OFFSET], nch, bht, uht);
            pp[i].deg   = uht->hd[plcm[i]].deg;
            pp[i].gen1  = i;
            pp[i].gen2  = bl;
            pp[i].lcm   = plcm[i];
        }
    }

    len_t nl  = pl+bl;
    /* Gebauer-Moeller: check old pairs first */
    /* note: old pairs are sorted by the given spair order */
#pragma omp parallel for num_threads(max_nthrds) \
    private(i, j,  l)
    for (i = 0; i < pl; ++i) {
        j = ps[i].gen1;
        l = ps[i].gen2;
        const int32_t m = pp[l].deg > pp[j].deg ? pp[l].deg : pp[j].deg;
        if (check_monomial_division(ps[i].lcm, nch, bht)
                && ps[i].deg > m
           ) {
            ps[i].lcm = 0;
        }
    }
    /* check new pairs for redundancy */
    j = 0;
    for (i = 0; i < bl; ++i) {
        if (bs->red[i] == 0) {
            pp[j++] = pp[i];
        }
    }
    /* sort new pairs by increasing lcm, earlier polys coming first */
    sort_r(pp, (unsigned long)j, sizeof(spair_t), spair_cmp, uht);
    for (i = 0; i < j; ++i) {
        plcm[i] = pp[i].lcm;
    }
    plcm[j]  = 0;
    const len_t pc  = j;

#pragma omp parallel for num_threads(max_nthrds) \
    private(j)
    for (j = 0; j < pc; ++j) {
        if (plcm[j] > 0) {
            const hi_t plcmj = plcm[j];
            check_monomial_division_in_update(plcm, j, pc, plcmj, uht);
        }
    }

    /* remove useless pairs from pairset */
    j = 0;
    /* old pairs */
    for (i = 0; i < psl->ld; ++i) {
        if (ps[i].lcm == 0) {
            continue;
        }
        ps[j++] = ps[i];
    }
    if (bht->esz - bht->eld <= pc) {
        enlarge_hash_table(bht);
    }
    /* new pairs, wee need to add the lcm to the basis hash table */
    insert_plcms_in_basis_hash_table(psl, pp, bht, uht, bs, plcm, j, pc);
    free(plcm);

    const bl_t lml          = bs->lml;
    const bl_t * const lmps = bs->lmps;

    if (bs->mltdeg > ndeg) {
        /* mark redundant elements in basis */
        for (i = 0; i < lml; ++i) {
            if (bs->red[lmps[i]] == 0
                    && check_monomial_division(bs->hm[lmps[i]][OFFSET], nch, bht)) {
                bs->red[lmps[i]]  = 1;
                st->num_redundant++;
            }
        }
    }

    st->num_gb_crit +=  nl - psl->ld;

    bs->ld++;
}

static void update_lm(
        bs_t *bs,
        const ht_t * const bht,
        stat_t *st
        )
{
    len_t i, j, k, l;

    const bl_t * const lmps = bs->lmps;

    j = bs->lo;
nextj:
    for (; j < bs->ld; ++j) {
        k = 0;
        for (l = bs->lo; l < j; ++l) {
            if (bs->red[l]) {
                continue;
            }
            if (check_monomial_division(bs->hm[j][OFFSET], bs->hm[l][OFFSET], bht)) {
                bs->red[j]  = 1;
                st->num_redundant++;
                j++;
                goto nextj;
            }
        }
        for (i = 0; i < bs->lml; ++i) {
            if (bs->red[lmps[i]] == 0
                    && check_monomial_division(bs->hm[lmps[i]][OFFSET], bs->hm[j][OFFSET], bht)) {
                bs->red[lmps[i]]  = 1;
                st->num_redundant++;
            }
        }
        const sdm_t *lms  = bs->lm;
        for (i = 0; i < bs->lml; ++i) {
            if (bs->red[lmps[i]] == 0) {
                bs->lm[k]   = lms[i];
                bs->lmps[k] = lmps[i];
                k++;
            }
        }
        bs->lml = k;
        k = bs->lml;
        if (bs->red[j] == 0) {
            bs->lm[k]   = bht->hd[bs->hm[j][OFFSET]].sdm;
            bs->lmps[k] = j;
            k++;
        }
        bs->lml = k;
    }
    bs->lo  = bs->ld;

    st->num_redundant_old = st->num_redundant;
}

static void update_basis(
        ps_t *ps,
        bs_t *bs,
        ht_t *bht,
        ht_t *uht,
        stat_t *st,
        const len_t npivs,
        const int32_t check_redundancy
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* compute number of new pairs we need to handle at most */
    len_t np  = bs->ld * npivs;
    for (i = 1; i < npivs; ++i) {
        np  = np + i;
    }
    check_enlarge_pairset(ps, np);

    for (i = 0; i < npivs; ++i) {
        insert_and_update_spairs(ps, bs, bht, uht, st, check_redundancy);
    }

    const bl_t lml          = bs->lml;
    const bl_t * const lmps = bs->lmps;

    len_t k = 0;
    if (st->mo == 0 && st->num_redundant_old < st->num_redundant) {
        const sdm_t *lms  = bs->lm;
        for (i = 0; i < lml; ++i) {
            if (bs->red[lmps[i]] == 0) {
                bs->lm[k]   = lms[i];
                bs->lmps[k] = lmps[i];
                k++;
            }
        }
        bs->lml = k;
    }
    k = bs->lml;
    for (i = bs->lo; i < bs->ld; ++i) {
        if (bs->red[i] == 0) {
            bs->lm[k]   = bht->hd[bs->hm[i][OFFSET]].sdm;
            bs->lmps[k] = i;
            k++;
        }
    }
    bs->lml = k;
    bs->lo  = bs->ld;

    st->num_redundant_old = st->num_redundant;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->update_ctime  +=  ct1 - ct0;
    st->update_rtime  +=  rt1 - rt0;
}
