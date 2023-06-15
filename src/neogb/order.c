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

/*
 * IMPLEMENTATIONS OF ORDER FUNCTIONIONALITY INDEPENDENT
 * OF CHOSEN MONOMIAL ORDER
 *  */

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

static int matrix_row_mult_cmp_increasing(
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
    /* compare multiplier resp. column index */
    va  = ((hm_t **)a)[0][MULT];
    vb  = ((hm_t **)b)[0][MULT];
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

static inline void sort_matrix_rows_mult_increasing(
        hm_t **rows,
        const len_t nrows
        )
{
    qsort(rows, (unsigned long)nrows, sizeof(hm_t *),
            &matrix_row_mult_cmp_increasing);
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

#if 0
static int spair_degree_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht  = htp;
    const hi_t la   = ((spair_t *)a)->lcm;
    const hi_t lb   = ((spair_t *)b)->lcm;
    const exp_t ea  = ht->hd[la].deg;
    const exp_t eb  = ht->hd[lb].deg;

    if (ea-eb != 0) {
        return (ea-eb);
    } else {
        const deg_t da  = ((spair_t *)a)->deg;
        const deg_t db  = ((spair_t *)b)->deg;
        return (da - db);
    }
}
#endif



/*
 * IMPLEMENTATIONS FOR LEXICOGRAPHICAL ORDER
 *  */

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
    const len_t evl = ht->evl;

    i = 1;
    while(i < evl-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
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
    const len_t evl = ht->evl;

    i = 1;
    while(i < evl-1 && ea[i] == eb[i]) {
        ++i;
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
    const len_t evl  = ht->evl;

    i = 1;
    while(i < evl-1 && ea[i] == eb[i]) {
        ++i;
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
    const len_t evl  = ht->evl;

    i = 1;
    while(i < evl-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
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

    if (ht->ev[la][DEG] != ht->ev[lb][DEG]) {
        return (ht->ev[la][DEG] < ht->ev[lb][DEG]) ? -1 : 1;
    } else {
        return (int)monomial_cmp(la, lb, ht);
    }
}




/*
 * IMPLEMENTATIONS FOR DEGREE REVERSE LEXICOGRAPHICAL ORDER
 *  */

#if 0
static int initial_input_cmp_drl_reverse(
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

    /* DRL */
    if (ea[DEG] > eb[DEG]) {
        return 1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return -1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
}
#endif

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

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* DRL */
    if (ea[DEG] < eb[DEG]) {
        return -1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return 1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
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

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];
    /* DRL */
    if (ea[DEG] < eb[DEG]) {
        return 1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return -1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
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

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* then DRL */
    if (ea[DEG] > eb[DEG]) {
        return -1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return 1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static inline int monomial_cmp_drl(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    len_t i;

    if (a == b) {
        return 0;
    }

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* DRL */
    if (ea[DEG] > eb[DEG]) {
        return 1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return -1;
        }
    }

    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
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

static int spair_cmp_update(
        const void *a,
        const void *b,
        void *htp)
{
    const spair_t *sa   =   ((spair_t *)a);
    const spair_t *sb   =   ((spair_t *)b);
    const ht_t *ht      =   (ht_t *)htp;
    if (sa->lcm != sb->lcm) {
        return (int)monomial_cmp(sa->lcm, sb->lcm, ht);
    } else {
        if (sa->deg != sb->deg) {
            return (sa->deg < sb->deg) ? -1 : 1;
        } else {
            if (sa->gen1 != sb->gen1) {
                return (sa->gen1 < sb->gen1) ? -1 : 1;
            } else {
                return 0;
            }
        }
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

    int mc = (int)monomial_cmp(la, lb, ht);
    if (mc != 0) {
        if (mc < 0)
            return -1;
        else
            return 1;

        /* return mc; */
    }
    /* if (((spair_t *)a)->gen1 < ((spair_t *)b)->gen1) {
     *     return -1;
     * }
     * if (((spair_t *)a)->gen1 > ((spair_t *)b)->gen1) {
     *     return 1;
     * }
     * if (((spair_t *)a)->gen2 < ((spair_t *)b)->gen2) {
     *     return -1;
     * }
     * if (((spair_t *)a)->gen2 > ((spair_t *)b)->gen2) {
     *     return 1;
     * } */
    return 0;
}




/*
 * IMPLEMENTATIONS FOR BLOCK ELIMINATION ORDER:
 * 2 blocks, each block handled by the degree
 * reverse lexicographical order
 *  */

static int initial_input_cmp_be(
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

    /* first block */
    if (ea[0] < eb[0]) {
        return -1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return 1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->ebl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    if (eb[i] - ea[i] != 0) {
        return eb[i] - ea[i];
    } else {
        /* second block */
        if (ea[ht->ebl] < eb[ht->ebl]) {
            return -1;
        } else {
            if (ea[ht->ebl] != eb[ht->ebl]) {
                return 1;
            }
        }

        /* note: reverse lexicographical */
        i = ht->evl-1;
        while (i > ht->ebl && ea[i] == eb[i]) {
            --i;
        }
        return eb[i] - ea[i];
    }
}

static int initial_gens_cmp_be(
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
    /* first block */
    if (ea[0] < eb[0]) {
        return 1;
    } else {
        if (ea[0] != eb[0]) {
            return -1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->ebl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    if (ea[i] - eb[i] != 0) {
        return ea[i] - eb[i];
    } else {
        /* second block */
        if (ea[ht->ebl] < eb[ht->ebl]) {
            return 1;
        } else {
            if (ea[ht->ebl] != eb[ht->ebl]) {
                return -1;
            }
        }

        /* note: reverse lexicographical */
        i = ht->evl-1;
        while (i > ht->ebl && ea[i] == eb[i]) {
            --i;
        }
        return ea[i] - eb[i];
    }
}

static int monomial_cmp_pivots_be(
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

    /* first block */
    if (ea[0] > eb[0]) {
        return -1;
    } else {
        if (ea[0] != eb[0]) {
            return 1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->ebl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    if (ea[i] - eb[i] != 0) {
        return ea[i] - eb[i];
    } else {
        /* second block */
        if (ea[ht->ebl] > eb[ht->ebl]) {
            return -1;
        } else {
            if (ea[ht->ebl] != eb[ht->ebl]) {
                return 1;
            }
        }

        /* note: reverse lexicographical */
        i = ht->evl-1;
        while (i > ht->ebl && ea[i] == eb[i]) {
            --i;
        }
        return ea[i] - eb[i];
    }
}

static inline int monomial_cmp_be(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    len_t i;

    if (a == b) {
        return 0;
    }

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* first block */
    if (ea[0] > eb[0]) {
        return 1;
    } else {
        if (ea[0] != eb[0]) {
            return -1;
        }
    }

    i = ht->ebl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    if (eb[i] - ea[i] != 0) {
        return eb[i] - ea[i];
    } else {
        /* second block */
        if (ea[ht->ebl] > eb[ht->ebl]) {
            return 1;
        } else {
            if (ea[ht->ebl] != eb[ht->ebl]) {
                return -1;
            }
        }

        i = ht->evl-1;
        while (i > ht->ebl && ea[i] == eb[i]) {
            --i;
        }
        return eb[i] - ea[i];
    }
}

/* comparison for hash-column-maps */
static int hcm_cmp_pivots_be(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht  = (ht_t *)htp;
    const hi_t ma  = ((hi_t *)a)[0];
    const hi_t mb  = ((hi_t *)b)[0];

    return monomial_cmp_pivots_be(ma, mb, ht);
}

static int spair_cmp_be(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hi_t da   = ((spair_t *)a)->deg;
    const hi_t db   = ((spair_t *)b)->deg;

    if (da - db != 0) {
        return da - db;
    } else {

        const hi_t la   = ((spair_t *)a)->lcm;
        const hi_t lb   = ((spair_t *)b)->lcm;
        const ht_t *ht  = (ht_t *)htp;

        int mc = (int)monomial_cmp(la, lb, ht);
        if (mc != 0) {
            if (mc < 0) {
                return -1;
            } else {
                return 1;
            }
        }
    }
    return 0;
}

/*
 * IMPLEMENTATIONS FOR SIGNATURE BASED ALGORITHMS
 *  */
static int initial_input_cmp_sig(
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

    /* DRL */
    if (ea[DEG] < eb[DEG]) {
        return 1;
    } else {
        if (ea[DEG] != eb[DEG]) {
            return -1;
        }
    }

    /* note: reverse lexicographical */
    i = ht->evl-1;
    while (i > 1 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static int matrix_row_cmp_by_increasing_signature(
        const void *a,
        const void *b,
        void *htp
        )
{
    const ht_t *ht  = (ht_t *)htp;

    hm_t sig_a, sig_b;
    /* compare pivot resp. column index */
    sig_a   = ((hm_t **)a)[0][SM_SMON];
    sig_b   = ((hm_t **)b)[0][SM_SMON];

    int diff = (int)monomial_cmp(sig_a, sig_b, ht);
    if (diff != 0) {
        return diff;
    } else {
        int sidx_a, sidx_b;
        sidx_a  = (int)((hm_t **)a)[0][SM_SIDX];
        sidx_b  = (int)((hm_t **)b)[0][SM_SIDX];

        return sidx_b - sidx_a;
    }
}

static inline void sort_matrix_rows_by_increasing_signature(
        smat_t *smat,
        void *htp
        )
{
    hm_t **cr = smat->cr;
    sort_r(cr, (unsigned long)smat->cld, sizeof(hm_t *),
            &matrix_row_cmp_by_increasing_signature, htp);
}

