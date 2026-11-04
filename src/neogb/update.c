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

/* exponent vector of the lcm of the monomials with exponents ea and eb */
static inline void lcm_exponents(
        exp_t *etmp,
        const exp_t * const ea,
        const exp_t * const eb,
        const ht_t * const ht
        )
{
    len_t i;
    const len_t evl = ht->evl;
    const len_t ebl = ht->ebl;

    for (i = 1; i < evl; ++i) {
        etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
    }
    etmp[0]   = 0;
    etmp[ebl] = 0;
    for (i = 1; i < ebl; ++i) {
        etmp[0]  += etmp[i];
    }
    for (i = ebl+1; i < evl; ++i) {
        etmp[ebl] += etmp[i];
    }
}

/* is the lcm of the lead monomials of g and the new element n different
 * from the monomial lcm, kl being the hash index of the first lcm or 0 if
 * it was not in the hash table when the pairs were generated */
static inline int pair_lcm_differs(
        const hi_t kl,
        const bi_t g,
        const hm_t n,
        const hi_t lcm,
        const bs_t * const bs,
        const ht_t * const bht
        )
{
    if (kl != 0) {
        return kl != lcm;
    }
    exp_t etmp[bht->evl];
    lcm_exponents(etmp, bht->ev[bs->hm[g][OFFSET]], bht->ev[n], bht);
    return memcmp(etmp, bht->ev[lcm], (uint64_t)bht->evl * sizeof(exp_t)) != 0;
}

/* data for comparing the new pairs of one new element whose lcms are
 * stored in a local hash table view */
typedef struct {
    const spair_t *p;
    const ht_t *ht;
} spair_cmp_local_t;

/* same order as spair_cmp_update(), for positions of pairs in a local
 * hash table view where equal lcms may have different positions */
static int spair_cmp_update_local(
        const void *a,
        const void *b,
        void *ctx
        )
{
    const spair_cmp_local_t *c = (spair_cmp_local_t *)ctx;
    const len_t ta = *((len_t *)a);
    const len_t tb = *((len_t *)b);
    const spair_t *sa = c->p + ta;
    const spair_t *sb = c->p + tb;
    const int mc = monomial_cmp(ta, tb, c->ht);
    if (mc != 0) {
        return mc;
    }
    if (sa->deg != sb->deg) {
        return (sa->deg < sb->deg) ? -1 : 1;
    }
    if (sa->gen1 != sb->gen1) {
        return (sa->gen1 < sb->gen1) ? -1 : 1;
    }
    return 0;
}

/* Adds the pairs of the npivs new basis elements bs->hm[bs->ld], ...,
 * bs->hm[bs->ld+npivs-1] to the pairset and applies the Gebauer-Moeller
 * criteria. The result is the same as adding the new elements one after
 * the other (the new pairs of one element, the removal of the old pairs
 * by it, the removal of the multiples among its new pairs and the
 * marking of redundant elements), including the order of the pairs, but
 * all elements are handled in bulk steps that are run in parallel. Most
 * new pairs are removed, so only the lcms of the remaining ones are
 * inserted into the hash table. */
static void insert_and_update_spairs(
        ps_t *psl,
        bs_t *bs,
        ht_t *bht,
        md_t *st,
        const len_t npivs
        )
{
    len_t i, k;
    uint64_t j;
    const int nthrds  = st->nthrds;
    const len_t bl    = bs->ld;
    const len_t pl    = psl->ld;
    const len_t evl   = bht->evl;
    const len_t ebl   = bht->ebl;
    const bl_t lml    = bs->lml;
    const bl_t * const lmps = bs->lmps;

    if (npivs == 0) {
        return;
    }

    /* data of the new elements: the maximal degree up to the k-th one
     * and if it marks old elements as redundant */
    hm_t *nch     = (hm_t *)malloc((uint64_t)npivs * sizeof(hm_t));
    deg_t *ndeg   = (deg_t *)malloc((uint64_t)npivs * sizeof(deg_t));
    int *marks    = (int *)malloc((uint64_t)npivs * sizeof(int));
    /* position of the pairs of the k-th new element */
    uint64_t *off = (uint64_t *)malloc((uint64_t)(npivs+1) * sizeof(uint64_t));
    off[0] = 0;
    for (k = 0; k < npivs; ++k) {
        nch[k]    = bs->hm[bl+k][OFFSET];
        ndeg[k]   = bs->hm[bl+k][DEG];
        bs->mltdeg  = bs->mltdeg > ndeg[k] ? bs->mltdeg : ndeg[k];
        marks[k]  = bs->mltdeg > ndeg[k];
        off[k+1]  = off[k] + bl + k;
    }
    const uint64_t np = off[npivs];

    /* 1 + the index of the new element marking an old one as redundant
     * first, 0 if it is already redundant, npivs + 1 if it stays */
    len_t *mby  = (len_t *)malloc((uint64_t)bl * sizeof(len_t));
#pragma omp parallel for num_threads(nthrds) private(i)
    for (i = 0; i < bl; ++i) {
        mby[i] = bs->red[i] != 0 ? 0 : npivs + 1;
    }
    len_t nred = 0;
#pragma omp parallel for num_threads(nthrds) private(i, k) \
    reduction(+:nred) schedule(dynamic, 64)
    for (i = 0; i < lml; ++i) {
        const bl_t b  = lmps[i];
        const hm_t lm = bs->hm[b][OFFSET];
        if (bs->red[b] != 0) {
            continue;
        }
        const deg_t d = bs->hm[b][DEG] - bht->hd[lm].deg;
        for (k = 0; k < npivs; ++k) {
            if (marks[k] && check_monomial_division(lm, nch[k], bht)
                    && d >= ndeg[k] - bht->hd[nch[k]].deg) {
                /* pairs with b of later elements are useless */
                mby[b] = k + 1;
                nred++;
                break;
            }
        }
    }
    st->num_redundant += nred;

    /* generate the new pairs, lcm is the hash index of the lcm if it is
     * already in the hash table, 0 otherwise */
    spair_t *pp = (spair_t *)malloc(np * sizeof(spair_t));
#pragma omp parallel for num_threads(nthrds) private(i, k) \
    schedule(dynamic)
    for (k = 0; k < npivs; ++k) {
        exp_t etmp[evl];
        const exp_t * const eb = bht->ev[nch[k]];
        for (i = 0; i < bl+k; ++i) {
            spair_t *p  = pp + off[k] + i;
            len_t l;
            val_t h     = 0;
            hi_t pos    = 0;
            lcm_exponents(etmp, bht->ev[bs->hm[i][OFFSET]], eb, bht);
            for (l = 0; l < evl; ++l) {
                h +=  bht->rn[l] * etmp[l];
            }
            p->lcm  = is_contained_in_hash_table(etmp, bht, h, &pos) ? pos : 0;
            p->gen1 = i;
            p->gen2 = bl+k;
            p->type = S_PAIR;
            if ((i < bl && mby[i] <= k) || (i >= bl && bs->red[i] != 0)) {
                p->deg = -1;
            } else {
                if (prime_monomials(bs->hm[i][OFFSET], nch[k], bht)) {
                    p->deg = -2;
                } else {
                    /* compute total degree of pair, not trivial if block order is chosen */
                    const deg_t ldeg  = etmp[0] + (ebl > 0 ? etmp[ebl] : 0);
                    if (st->nev == 0) {
                        p->deg = ldeg;
                    } else {
                        const deg_t deg1  = ldeg
                            - bht->hd[bs->hm[i][OFFSET]].deg + bs->hm[i][DEG];
                        const deg_t deg2  = ldeg
                            - bht->hd[nch[k]].deg + bs->hm[bl+k][DEG];
                        p->deg = deg1 > deg2 ? deg1 : deg2;
                    }
                }
            }
        }
    }
    /* the pairs with the k-th new element before they are sorted, needed
     * for the criteria applied to the older pairs */
    hi_t *plcm  = (hi_t *)malloc(np * sizeof(hi_t));
    deg_t *pdeg = (deg_t *)malloc(np * sizeof(deg_t));
#pragma omp parallel for num_threads(nthrds) private(j)
    for (j = 0; j < np; ++j) {
        plcm[j] = pp[j].lcm;
        pdeg[j] = pp[j].deg;
    }

    /* Gebauer-Moeller: check old pairs first, the missing lcms of the new
     * pairs differ from the ones of the old pairs */
    spair_t *ps = psl->p;
    for (k = 0; k < npivs; ++k) {
        const hi_t *kl  = plcm + off[k];
        const deg_t *kd = pdeg + off[k];
#pragma omp parallel for num_threads(nthrds) private(i) \
        schedule(static, 4096)
        for (i = 0; i < pl; ++i) {
            const bi_t g1 = ps[i].gen1;
            const bi_t g2 = ps[i].gen2;
            if (ps[i].deg >= 0 && kl[g1] != ps[i].lcm && kl[g2] != ps[i].lcm
                    && kd[g1] <= ps[i].deg && kd[g2] <= ps[i].deg
                    && check_monomial_division(ps[i].lcm, nch[k], bht)) {
                ps[i].deg   =   -1;
            }
        }
    }

    /* Gebauer-Moeller: sort the new pairs of each new element by
     * increasing lcm, earlier polys coming first, and remove multiples;
     * the lcms are compared in a local view of the hash table */
#pragma omp parallel for num_threads(nthrds) private(i, k) \
    schedule(dynamic)
    for (k = 0; k < npivs; ++k) {
        spair_t *p    = pp + off[k];
        const len_t n = bl + k;
        len_t l;
        exp_t *e      = (exp_t *)malloc((uint64_t)n * evl * sizeof(exp_t));
        exp_t **ev    = (exp_t **)malloc((uint64_t)n * sizeof(exp_t *));
        hd_t *hd      = (hd_t *)calloc((uint64_t)n, sizeof(hd_t));
        len_t *ord    = (len_t *)malloc((uint64_t)n * sizeof(len_t));
        spair_t *sp   = (spair_t *)malloc((uint64_t)n * sizeof(spair_t));
        ht_t lht      = *bht;
        lht.ev  = ev;
        lht.hd  = hd;
        for (i = 0; i < n; ++i) {
            ev[i]   = e + (uint64_t)i * evl;
            lcm_exponents(ev[i], bht->ev[bs->hm[i][OFFSET]], bht->ev[nch[k]], bht);
            hd[i].sdm = generate_short_divmask(ev[i], bht);
            ord[i]  = i;
        }
        spair_cmp_local_t ctx = {p, &lht};
        sort_r(ord, (unsigned long)n, sizeof(len_t), spair_cmp_update_local, &ctx);
        for (i = 0; i < n; ++i) {
            sp[i] = p[ord[i]];
        }
        for (i = 1; i < n; ++i) {
            if (sp[i].deg < 0) {
                continue;
            }
            for (l = 0; l < i; ++l) {
                if (sp[l].deg == -1) {
                    continue;
                }
                if (check_monomial_division(ord[i], ord[l], &lht)) {
                    sp[i].deg   =   -1;
                    break;
                }
            }
        }
        memcpy(p, sp, (uint64_t)n * sizeof(spair_t));
        free(e);
        free(ev);
        free(hd);
        free(ord);
        free(sp);
    }

    /* remove the useless new pairs, soff[k] is the number of remaining
     * pairs of the first k new elements */
    uint64_t *soff  = (uint64_t *)malloc((uint64_t)(npivs+1) * sizeof(uint64_t));
    j = 0;
    soff[0] = 0;
    for (k = 0; k < npivs; ++k) {
        for (i = 0; i < bl+k; ++i) {
            if (pp[off[k]+i].deg >= 0) {
                pp[j++] = pp[off[k]+i];
            }
        }
        soff[k+1] = j;
    }
    /* insert the missing lcms of the remaining new pairs */
    uint64_t nmiss = 0;
    for (j = 0; j < soff[npivs]; ++j) {
        if (pp[j].lcm == 0) {
            nmiss++;
        }
    }
    while (bht->esz - bht->eld < nmiss) {
        enlarge_hash_table(bht);
    }
    for (j = 0; j < soff[npivs]; ++j) {
        if (pp[j].lcm == 0) {
            pp[j].lcm = get_lcm(bs->hm[pp[j].gen1][OFFSET],
                    bs->hm[pp[j].gen2][OFFSET], bht, bht);
        }
    }

    /* Gebauer-Moeller: the old pair check of the later new elements for
     * the pairs of the earlier ones */
    for (k = 1; k < npivs; ++k) {
        const hi_t *kl  = plcm + off[k];
        const deg_t *kd = pdeg + off[k];
#pragma omp parallel for num_threads(nthrds) private(j) \
        schedule(static, 4096)
        for (j = 0; j < soff[k]; ++j) {
            const bi_t g1 = pp[j].gen1;
            const bi_t g2 = pp[j].gen2;
            if (pp[j].deg >= 0 && kd[g1] <= pp[j].deg && kd[g2] <= pp[j].deg
                    && check_monomial_division(pp[j].lcm, nch[k], bht)
                    && pair_lcm_differs(kl[g1], g1, nch[k], pp[j].lcm, bs, bht)
                    && pair_lcm_differs(kl[g2], g2, nch[k], pp[j].lcm, bs, bht)) {
                pp[j].deg   =   -1;
            }
        }
    }

    /* remove useless pairs from pairset, old ones first, then the new
     * ones in the order of the new elements */
    j = 0;
    for (i = 0; i < pl; ++i) {
        if (ps[i].deg >= 0) {
            ps[j++] = ps[i];
        }
    }
    for (uint64_t m = 0; m < soff[npivs]; ++m) {
        if (pp[m].deg >= 0) {
            ps[j++] = pp[m];
        }
    }
    psl->ld = j;

    /* mark redundant elements in basis */
    for (i = 0; i < lml; ++i) {
        if (mby[lmps[i]] > 0 && mby[lmps[i]] <= npivs && bs->red[lmps[i]] == 0) {
            bs->red[lmps[i]] = 1;
        }
    }

    st->num_gb_crit +=  pl + np - psl->ld;

    bs->ld += npivs;

    free(nch);
    free(ndeg);
    free(marks);
    free(off);
    free(soff);
    free(mby);
    free(pp);
    free(plcm);
    free(pdeg);
}

static void update_lm(
        bs_t *bs,
        const ht_t * const bht,
        md_t *st
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

static void update_basis_f4(
        ps_t *ps,
        bs_t *bs,
        ht_t *bht,
        md_t *st,
        const len_t npivs
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

    insert_and_update_spairs(ps, bs, bht, st, npivs);

    const bl_t lml          = bs->lml;
    const bl_t * const lmps = bs->lmps;

    len_t k = 0;

    /* Check new elements on redundancy:
     * Only elements coming from the same matrix are possible leading
     * monomial divisors, thus we only check down to bs->lo */
/* #pragma omp parallel for num_threads(st->nthrds)
    for (int l = bs->lo; l < bs->ld; ++l) {
        hm_t lm  = bs->hm[l][OFFSET];
        deg_t dd = bs->hm[l][DEG] - bht->hd[lm].deg;
        for (int m = bs->lo; m < l; ++m) {
            if (check_monomial_division(lm, bs->hm[m][OFFSET], bht) == 1
                && dd >= (bs->hm[m][DEG] - bht->hd[bs->hm[m][OFFSET]].deg)) {
                bs->red[l]  =   1;
                st->num_redundant++;
                break;
            }
        }
    } */
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

static int32_t update(
        bs_t *bs,
        md_t *md
        )
{
        update_basis_f4(md->ps, bs, bs->ht, md, md->np);
        if (bs->constant) {
            return 1;
        } else {
            return md->ps->ld == 0;
        }
}

/* not needed right now, maybe in a later iteration of sba implementations */
#if 0
static void update_basis_sba_schreyer(
        ps_t *ps,
        bs_t *bs,
        ht_t *bht,
        ht_t *uht,
        md_t *st,
        const len_t npivs
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
        insert_and_update_spairs(ps, bs, bht, uht, st);
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
            j>lmps[k] = i;
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
#endif
