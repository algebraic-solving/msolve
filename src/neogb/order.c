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

static int initial_input_cmp_lex(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = ((hm_t **)a)[0][OFFSET];
    const hm_t hb  = ((hm_t **)b)[0][OFFSET];

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* lexicographical */
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
}

static int initial_input_cmp_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = ((hm_t **)a)[0][OFFSET];
    const hm_t hb  = ((hm_t **)b)[0][OFFSET];

    const deg_t da = ht->hd[ha].deg;
    const deg_t db = ht->hd[hb].deg;

    /* DRL */
    if (da < db) {
        return -1;
    } else {
        if (da != db) {
            return 1;
        }
    }

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* note: reverse lexicographical */
    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
}

static int initial_gens_cmp_lex(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = **(hm_t **)a;
    const hm_t hb  = **(hm_t **)b;

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* lexicographical */
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
}

static int initial_gens_cmp_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = **(hm_t **)a;
    const hm_t hb  = **(hm_t **)b;

    const deg_t da = ht->hd[ha].deg;
    const deg_t db = ht->hd[hb].deg;

    /* DRL */
    if (da < db) {
        return 1;
    } else {
        if (da != db) {
            return -1;
        }
    }

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* note: reverse lexicographical */
    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static int matrix_row_cmp_decreasing(
        const void *a,
        const void *b
        )
{
    hm_t va, vb;
    /* compare pivot resp. column index */
    va  = ((hm_t **)a)[0][OFFSET];
    vb  = ((hm_t **)b)[0][OFFSET];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    /* same column index => compare density of row */
    va  = ((hm_t **)a)[0][LENGTH];
    vb  = ((hm_t **)b)[0][LENGTH];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    return 0;
}

static int matrix_row_cmp_increasing(
        const void *a,
        const void *b
        )
{
    hm_t va, vb;
    /* compare pivot resp. column index */
    va  = ((hm_t **)a)[0][OFFSET];
    vb  = ((hm_t **)b)[0][OFFSET];
    if (va > vb) {
        return -1;
    }
    if (va < vb) {
        return 1;
    }
    /* same column index => compare density of row */
    va  = ((hm_t **)a)[0][LENGTH];
    vb  = ((hm_t **)b)[0][LENGTH];
    if (va > vb) {
        return -1;
    }
    if (va < vb) {
        return 1;
    }
    return 0;
}

static inline void sort_matrix_rows_decreasing(
        hm_t **rows,
        const len_t nrows
        )
{
    qsort(rows, (unsigned long)nrows, sizeof(hm_t *),
            &matrix_row_cmp_decreasing);
}

static inline void sort_matrix_rows_increasing(
        hm_t **rows,
        const len_t nrows
        )
{
    qsort(rows, (unsigned long)nrows, sizeof(hm_t *),
            &matrix_row_cmp_increasing);
}

static int dense_matrix_row_cmp(
        const void *a,
        const void *b
        )
{
    const cf32_t pa = ((cf32_t **)a)[0][0];
    const cf32_t pb = ((cf32_t **)b)[0][0];

    return pa-pb;
}

static inline cf32_t **sort_dense_matrix_rows(
        cf32_t **dm,
        const len_t nr
        )
{
    qsort(dm, (unsigned long)nr, sizeof(cf32_t *), &dense_matrix_row_cmp);
    return dm;
}

static int monomial_cmp_pivots_drl(
        const hi_t a,
        const hi_t b,
        const ht_t * const ht
        )
{
    len_t i;

    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
#if ORDER_COLUMNS
    /* first known pivots vs. tail terms */
    if (ha.idx != hb.idx) {
        if (ha.idx < hb.idx) {
            return 1;
        } else {
            return -1;
        }
    }
#endif

    /* then DRL */
    if (ha.deg > hb.deg) {
        return -1;
    } else {
        if (ha.deg != hb.deg) {
            return 1;
        }
    }

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* note: reverse lexicographical */
    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static int monomial_cmp_pivots_lex(
        const hi_t a,
        const hi_t b,
        const ht_t * const ht
        )
{
    len_t i;

    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
#if ORDER_COLUMNS
    /* first known pivots vs. tail terms */
    if (ha.idx != hb.idx) {
        if (ha.idx < hb.idx) {
            return 1;
        } else {
            return -1;
        }
    }
#endif

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* lexicographical */
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return eb[i] - ea[i];
}

static inline int monomial_cmp_drl(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    len_t i;

    const deg_t da = ht->hd[a].deg;
    const deg_t db = ht->hd[b].deg;

    /* DRL */
    if (da > db) {
        return 1;
    } else {
        if (da != db) {
            return -1;
        }
    }

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
}

static inline int monomial_cmp_lex(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    len_t i;

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
}

/* comparison for spair generators */
static int gens_cmp(
        const void *a,
        const void *b
        )
{
    const len_t ga = *((len_t *)a);
    const len_t gb = *((len_t *)b);

    return (ga - gb);
}

/* comparison for sparse rows in preparation for generation of pbm files */
static int pbm_cmp(
        const void *a,
        const void *b
        )
{
    const hm_t ca = *((hm_t *)a);
    const hm_t cb = *((hm_t *)b);

    return (ca - cb);
}

/* comparison for hash-column-maps */
static int hcm_cmp_pivots_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht  = (ht_t *)htp;
    const hi_t ma  = ((hi_t *)a)[0];
    const hi_t mb  = ((hi_t *)b)[0];

    return monomial_cmp_pivots_drl(ma, mb, ht);
}

static int hcm_cmp_pivots_lex(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht  = (ht_t *)htp;
    const hi_t ma   = ((hi_t *)a)[0];
    const hi_t mb   = ((hi_t *)b)[0];

    return monomial_cmp_pivots_lex(ma, mb, ht);
}

/* comparison for s-pairs once their lcms are in the global hash table */
static int spair_cmp_deglex(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hi_t la   = ((spair_t *)a)->lcm;
    const hi_t lb   = ((spair_t *)b)->lcm;
    const ht_t *ht  = (ht_t *)htp;

    if (ht->hd[la].deg != ht->hd[lb].deg) {
        return (ht->hd[la].deg < ht->hd[lb].deg) ? -1 : 1;
    } else {
        return (int)monomial_cmp(la, lb, ht);
    }
}

static int spair_cmp_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hi_t la   = ((spair_t *)a)->lcm;
    const hi_t lb   = ((spair_t *)b)->lcm;
    const ht_t *ht  = (ht_t *)htp;

    return (int)monomial_cmp(la, lb, ht);
}

static int spair_degree_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hd_t *hd  = ((ht_t *)htp)->hd;
    const deg_t da  = hd[((spair_t *)a)->lcm].deg;
    const deg_t db  = hd[((spair_t *)b)->lcm].deg;

    return (da-db);
}
