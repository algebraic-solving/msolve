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
        stat_t *st,
        const int32_t check_redundancy
        )
{
    int i, j, l;
    deg_t deg1, deg2;

    spair_t *ps = psl->p;

#ifdef _OPENMP
    const int nthrds = st->nthrds;
#endif

    const int pl  = psl->ld;
    const int bl  = bs->ld;

    const hm_t nch = bs->hm[bl][OFFSET];

    deg_t ndeg  = bs->hm[bl][DEG];

    bs->mltdeg  = bs->mltdeg > ndeg ?
        bs->mltdeg : ndeg;

    spair_t *pp = ps+pl;

    while (bht->esz - bht->eld < bl) {
        enlarge_hash_table(bht);
    }
#if PARALLEL_HASHING
#pragma omp parallel for num_threads(nthrds) \
    private(i)
#endif
    for (i = 0; i < bl; ++i) {
        pp[i].lcm   =  get_lcm(bs->hm[i][OFFSET], nch, bht, bht);
        pp[i].gen1  = i;
        pp[i].gen2  = bl;
        if (bs->red[i] != 0) {
            pp[i].deg   =   -1;
        } else {
            if (prime_monomials(bs->hm[pp[i].gen1][OFFSET], bs->hm[pp[i].gen2][OFFSET], bht)) {
                pp[i].deg   =   -2;
            } else {
                /* compute total degree of pair, not trivial if block order is chosen */
                if (st->nev == 0) {
                    pp[i].deg = bht->hd[pp[i].lcm].deg;
                } else {
                    deg1  = bht->hd[pp[i].lcm].deg - bht->hd[bs->hm[i][OFFSET]].deg + bs->hm[i][DEG];
                    deg2  = bht->hd[pp[i].lcm].deg - bht->hd[nch].deg + bs->hm[bl][DEG];
                    pp[i].deg = deg1 > deg2 ? deg1 : deg2;
                }
            }
        }
    }

    len_t nl  = pl+bl;
    /* Gebauer-Moeller: check old pairs first */
    /* note: old pairs are sorted by the given spair order */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j,  l)
    for (i = 0; i < pl; ++i) {
        j = ps[i].gen1;
        l = ps[i].gen2;
        if (pp[j].lcm != ps[i].lcm && pp[l].lcm != ps[i].lcm
                && pp[j].deg <= ps[i].deg && pp[l].deg <= ps[i].deg
                && check_monomial_division(ps[i].lcm, nch, bht)) {
            ps[i].deg   =   -1;
        }
    }
    /* sort new pairs by increasing lcm, earlier polys coming first */
    sort_r(pp, (unsigned long)bl, sizeof(spair_t), spair_cmp_update, bht);

    /* Gebauer-Moeller: remove real multiples of new spairs */
    for (i = pl; i < nl; ++i) {
        if (ps[i].deg < 0) {
            continue;
        }
        for (j = pl; j < i; ++j) {
            if (i == j || ps[j].deg == -1) {
                continue;
            }
            if (ps[i].lcm != ps[j].lcm
                    && ps[i].deg >= ps[j].deg
                    && check_monomial_division(ps[i].lcm, ps[j].lcm, bht)) {
                ps[i].deg   =   -1;
                break;
            }
        }
    }


    /* Gebauer-Moeller: remove same lcm spairs from the new ones */
    for (i = pl; i < nl; ++i) {
        if (ps[i].deg == -1) {
            continue;
        }
        /* try to remove all others if product criterion applies */
        if (ps[i].deg == -2) {
            for (j = pl; j < nl; ++j) {
                if (ps[j].lcm == ps[i].lcm) {
                    ps[j].deg   =   -1;
                }
            }
            /* try to eliminate this spair with earlier ones */
        } else { 
            for (j = i-1; j >= pl; --j) {
                if (ps[j].deg != -1
                        && ps[j].deg <= ps[i].deg
                        && ps[i].lcm == ps[j].lcm) {
                    ps[i].deg   =   -1;
                    break;
                }
            }
        }
    }


    /* remove useless pairs from pairset */
    j = 0;
    /* old pairs */
    for (i = 0; i < nl; ++i) {
        if (ps[i].deg < 0) {
            continue;
        }
        ps[j++] = ps[i];
    }

    psl->ld =   j;

    const bl_t lml          = bs->lml;
    const bl_t * const lmps = bs->lmps;

    /* mark redundant elements in basis */
    deg_t dd = ndeg - bht->hd[nch].deg;
    if (bs->mltdeg > ndeg) {
#if PARALLEL_HASHING
#pragma omp parallel for num_threads(nthrds) \
    private(i)
#endif
        for (i = 0; i < lml; ++i) {
            hm_t lm = bs->hm[lmps[i]][OFFSET];
            if (bs->red[lmps[i]] == 0
                    && check_monomial_division(lm, nch, bht)
                    && bs->hm[lmps[i]][DEG]-bht->hd[lm].deg >= dd) {
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

static void update_basis_f4(
        ps_t *ps,
        bs_t *bs,
        ht_t *bht,
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
        insert_and_update_spairs(ps, bs, bht, st, check_redundancy);
    }

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

/* not needed right now, maybe in a later iteration of sba implementations */
#if 0
static void update_basis_sba_schreyer(
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
#endif
