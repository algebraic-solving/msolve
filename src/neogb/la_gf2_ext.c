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

/* Linear algebra over GF(16) and GF(256), see gf2_ext.c for the
 * arithmetic. The coefficients are stored as cf8_t in the arrays cf_8
 * like for primes < 2^8, so the matrix and basis handling is the same,
 * only the reduction differs. Since addition is XOR, no reduction modulo
 * p is needed and the dense rows are byte arrays. Only exact linear
 * algebra for Groebner basis computations is supported, no tracer and
 * no normal forms. */

/* returns memory of size bytes aligned to 64 bytes, *raw is the pointer
 * to be freed */
static uint8_t *gf2x_alloc_aligned(
        void **raw,
        const uint64_t size
        )
{
    *raw = calloc(size + 64, sizeof(uint8_t));
    return (uint8_t *)(((uintptr_t)*raw + 63) & ~(uintptr_t)63);
}

/* makes the sparse row with coefficients cf of length len monic */
static inline void normalize_sparse_row_gf2x(
        cf8_t *cf,
        const len_t len,
        const gf2x_t * const t
        )
{
    if (cf[0] != 1) {
        gf2x_mul_region(cf, cf, len, t->inv[cf[0]], t, 0);
    }
}

/* dr[ds[j]] ^= c*cf[j] for 0 <= j < len */
static inline void gf2x_scatter_axpy(
        uint8_t *dr,
        const hm_t * const ds,
        const cf8_t * const cf,
        const len_t len,
        const uint8_t c,
        const gf2x_t * const t
        )
{
    len_t j = 0;
    const len_t os = len % UNROLL;
    if (c == 1) {
        for (; j < os; ++j) {
            dr[ds[j]] ^= cf[j];
        }
        for (; j < len; j += UNROLL) {
            dr[ds[j]]   ^= cf[j];
            dr[ds[j+1]] ^= cf[j+1];
            dr[ds[j+2]] ^= cf[j+2];
            dr[ds[j+3]] ^= cf[j+3];
        }
        return;
    }
    const uint8_t * const mr = t->mul[c];
    for (; j < os; ++j) {
        dr[ds[j]] ^= mr[cf[j]];
    }
    for (; j < len; j += UNROLL) {
        dr[ds[j]]   ^= mr[cf[j]];
        dr[ds[j+1]] ^= mr[cf[j+1]];
        dr[ds[j+2]] ^= mr[cf[j+2]];
        dr[ds[j+3]] ^= mr[cf[j+3]];
    }
}

/* loads the sparse row with columns ds and coefficients cf of length len
 * into the dense row dr which is zero */
static inline void gf2x_load_row(
        uint8_t *dr,
        const hm_t * const ds,
        const cf8_t * const cf,
        const len_t len
        )
{
    for (len_t j = 0; j < len; ++j) {
        dr[ds[j]] = cf[j];
    }
}

/* Reduces the dense row dr by the pivots pivs starting at column sc. The
 * remainder is returned as sparse row, its coefficients are stored in
 * mat->cf_8[tmp_pos], NULL is returned if the remainder is zero. The
 * remainder is not normalized. Afterwards dr is zero. */
static hm_t *reduce_dense_row_by_pivots_gf2x(
        uint8_t *dr,
        mat_t *mat,
        hm_t * const * const pivs,
        const hi_t sc,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi,
        const gf2x_t * const t
        )
{
    len_t i, j, k = 0;
    hi_t np = mat->nc;

    const len_t ncols = mat->nc;
    cf8_t * const * const mcf = mat->cf_8;

    for (i = gf2x_next_nonzero(dr, sc, ncols); i < ncols;
            i = gf2x_next_nonzero(dr, i+1, ncols)) {
        const hm_t * const piv = pivs[i];
        if (piv == NULL) {
            if (np == ncols) {
                np = i;
            }
            k++;
            continue;
        }
        /* pivots are monic, this also zeroes dr[i] */
        gf2x_scatter_axpy(dr, piv + OFFSET, mcf[piv[COEFFS]],
                piv[LENGTH], dr[i], t);
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    cf8_t *cf = (cf8_t *)malloc((uint64_t)k * sizeof(cf8_t));
    hm_t *rs  = row + OFFSET;
    j = 0;
    for (i = gf2x_next_nonzero(dr, np, ncols); i < ncols;
            i = gf2x_next_nonzero(dr, i+1, ncols)) {
        rs[j] = (hm_t)i;
        cf[j] = dr[i];
        dr[i] = 0;
        j++;
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_8[tmp_pos] = cf;

    return row;
}

/* Reduces the GF_NR interleaved dense rows dm, dm[GF_NR*c+r] being the
 * entry of row r in column c, by the known pivots pivs[sc],...,
 * pivs[ncl-1]. For each pivot the multipliers of all rows are split into
 * nibbles once, then each term of the pivot updates all rows with a few
 * vector operations. Afterwards the rows only have entries in the
 * columns ncl,...,ncols-1. */
static void reduce_dense_rows_by_known_pivots_gf2x(
        uint8_t *dm,
        const mat_t * const mat,
        hm_t * const * const pivs,
        const len_t sc,
        const gf2x_t * const t
        )
{
    len_t i, j;
    uint8_t ml[GF_NR] __attribute__((aligned(64)));
    uint8_t mh[GF_NR] __attribute__((aligned(64)));

    const len_t ncl = mat->ncl;
    const int q16   = t == gf2x_tables;
    cf8_t * const * const mcf = mat->cf_8;

    for (i = sc; i < ncl; ++i) {
        uint8_t *d = dm + (uint64_t)GF_NR * i;
        if (gf2x_lanes_zero(d)) {
            continue;
        }
        gf2x_lanes_split(ml, mh, d);
        const hm_t * const piv = pivs[i];
        const hm_t * const ds  = piv + OFFSET;
        const cf8_t * const cf = mcf[piv[COEFFS]];
        const len_t len        = piv[LENGTH];
        /* the pivot is monic, its lead term zeroes the column */
        memset(d, 0, GF_NR);
        if (q16) {
            for (j = 1; j < len; ++j) {
                gf2x_lanes_axpy(dm + (uint64_t)GF_NR * ds[j], ml, mh,
                        cf[j], t, 1);
            }
        } else {
            for (j = 1; j < len; ++j) {
                gf2x_lanes_axpy(dm + (uint64_t)GF_NR * ds[j], ml, mh,
                        cf[j], t, 0);
            }
        }
    }
}

/* Reduces the dense row dr of length ncr starting at column sc by the
 * dense pivots nps, nps[c] is monic and starts at column c. Stops at the
 * first column without pivot and returns it, ncr if dr becomes zero. */
static len_t reduce_dense_row_by_dense_pivots_gf2x(
        uint8_t *dr,
        const len_t sc,
        uint8_t * const * const nps,
        const len_t ncr,
        const gf2x_t * const t
        )
{
    len_t i;
    for (i = gf2x_next_nonzero(dr, sc, ncr); i < ncr;
            i = gf2x_next_nonzero(dr, i+1, ncr)) {
        const uint8_t *piv = nps[i];
        if (piv == NULL) {
            return i;
        }
        gf2x_mul_region(dr+i, piv, ncr-i, dr[i], t, 1);
    }
    return ncr;
}

/* Interreduces the dense pivots dm, dm[c] is monic, starts at column c
 * and has length ncr-c, or is NULL. The pivots are handled in blocks
 * from the right to the left: the rows of a block are interreduced,
 * afterwards all rows left of the block are reduced by the block in
 * parallel. */
static void interreduce_dense_matrix_gf2x(
        uint8_t **dm,
        const len_t ncr,
        const int32_t nthrds,
        const gf2x_t * const t
        )
{
    len_t i, n = 0;
    int64_t r;
    const len_t bs = 64;

    len_t *pc = (len_t *)malloc((uint64_t)ncr * sizeof(len_t));
    for (i = 0; i < ncr; ++i) {
        if (dm[i] != NULL) {
            pc[n++] = i;
        }
    }
    len_t e = n;
    while (e > 1) {
        const len_t s = e > bs ? e - bs : 0;
        /* interreduce the rows of the block */
        for (r = (int64_t)e-2; r >= (int64_t)s; --r) {
            uint8_t *row  = dm[pc[r]];
            for (i = r+1; i < e; ++i) {
                const len_t c = pc[i] - pc[r];
                if (row[c] != 0) {
                    gf2x_mul_region(row+c, dm[pc[i]], ncr-pc[i], row[c], t, 1);
                }
            }
        }
        /* reduce all rows left of the block by the block */
#pragma omp parallel for num_threads(nthrds) private(i) schedule(dynamic, 8)
        for (r = 0; r < (int64_t)s; ++r) {
            uint8_t *row  = dm[pc[r]];
            for (i = s; i < e; ++i) {
                const len_t c = pc[i] - pc[r];
                if (row[c] != 0) {
                    gf2x_mul_region(row+c, dm[pc[i]], ncr-pc[i], row[c], t, 1);
                }
            }
        }
        e = s;
    }
    free(pc);
}

/* Interreduces the new pivots pivs[ncl],...,pivs[ncl+ncr-1] as dense rows
 * if they are dense enough and stores them in mat->tr from the right to
 * the left. Returns 0 if the pivots are too sparse, they are not changed
 * then. */
static int interreduce_new_pivots_dense_gf2x(
        mat_t *mat,
        hm_t **pivs,
        const int32_t nthrds,
        len_t *npivsp,
        const gf2x_t * const t
        )
{
    len_t i, j;
    uint64_t nnz = 0, dsz = 0;
    len_t n = 0;

    const len_t ncl = mat->ncl;
    const len_t ncr = mat->ncr;

    for (i = 0; i < ncr; ++i) {
        if (pivs[ncl+i] != NULL) {
            n++;
            nnz += pivs[ncl+i][LENGTH];
            dsz += ncr - i;
        }
    }
    if (n < 2 || 4 * nnz < dsz) {
        return 0;
    }
    uint8_t **dm  = (uint8_t **)calloc((uint64_t)ncr, sizeof(uint8_t *));
    hm_t *hdr     = (hm_t *)malloc((uint64_t)ncr * OFFSET * sizeof(hm_t));
    for (i = 0; i < ncr; ++i) {
        hm_t *row = pivs[ncl+i];
        if (row == NULL) {
            continue;
        }
        memcpy(hdr + (uint64_t)i * OFFSET, row, OFFSET * sizeof(hm_t));
        cf8_t *cfs  = mat->cf_8[row[COEFFS]];
        dm[i]       = (uint8_t *)calloc((uint64_t)(ncr-i), sizeof(uint8_t));
        const hm_t * const ds = row + OFFSET;
        for (j = 0; j < row[LENGTH]; ++j) {
            dm[i][ds[j]-ncl-i] = cfs[j];
        }
        free(cfs);
        mat->cf_8[row[COEFFS]] = NULL;
        free(row);
        pivs[ncl+i] = NULL;
    }
    interreduce_dense_matrix_gf2x(dm, ncr, nthrds, t);

    len_t npivs = 0;
    for (i = ncr; i > 0; --i) {
        uint8_t *d = dm[i-1];
        if (d == NULL) {
            continue;
        }
        const len_t len = ncr-i+1;
        len_t k = 0;
        for (j = 0; j < len; ++j) {
            k += d[j] != 0;
        }
        hm_t *row   = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
        cf8_t *cfs  = (cf8_t *)malloc((uint64_t)k * sizeof(cf8_t));
        memcpy(row, hdr + (uint64_t)(i-1) * OFFSET, OFFSET * sizeof(hm_t));
        hm_t *rs = row + OFFSET;
        k = 0;
        for (j = gf2x_next_nonzero(d, 0, len); j < len;
                j = gf2x_next_nonzero(d, j+1, len)) {
            rs[k]   = ncl+i-1+j;
            cfs[k]  = d[j];
            k++;
        }
        row[PRELOOP]  = k % UNROLL;
        row[LENGTH]   = k;
        mat->cf_8[row[COEFFS]] = cfs;
        pivs[ncl+i-1] = mat->tr[npivs++] = row;
        free(d);
    }
    free(dm);
    free(hdr);
    *npivsp = npivs;
    return 1;
}

/* Number of rows reduced together by the known pivots, 0 if the rows
 * shall be reduced one by one. A group costs as much as GF_NR rows, but
 * applying a pivot to all rows of a group is much cheaper than applying
 * it to each row on its own, even if only some rows need it. Thus we use
 * smaller groups if there are not enough rows for each thread. */
static inline len_t rows_per_group_gf2x(
        const len_t nrl,
        const int32_t nthrds
        )
{
    const uint64_t gl = ((uint64_t)nrl + nthrds - 1) / nthrds;
    if (gl >= GF_NR) {
        return GF_NR;
    }
    return gl >= 4 ? (len_t)gl : 0;
}

static void exact_sparse_reduced_echelon_form_gf2x(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    const int32_t nthrds  = st->nthrds;
    const int final       = st->in_final_reduction_step;
    const gf2x_t * const t = gf2x_table(st->fc);

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_8[j]       = bs->cf_8[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }
    hm_t *kpivs = reducer_rows_to_block(mat);
    if (mat->nru > 0) {
        memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    }

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* the dense rows are zero whenever a row is started */
    uint8_t *dr = (uint8_t *)calloc((uint64_t)ncols * nthrds, sizeof(uint8_t));

    /* the rows are first reduced by the known pivots in groups, see
     * reduce_dense_rows_by_known_pivots_gf2x() */
    const len_t gsz = rows_per_group_gf2x(nrl, nthrds);
    void **drm_raw  = NULL;
    uint8_t **drm   = NULL;
    if (gsz > 0) {
        drm_raw = (void **)calloc((uint64_t)nthrds, sizeof(void *));
        drm     = (uint8_t **)calloc((uint64_t)nthrds, sizeof(uint8_t *));
    }
    const len_t gl  = gsz > 0 ? gsz : 1;
    const len_t ngr = (nrl + gl - 1) / gl;
    len_t g;
#pragma omp parallel for num_threads(nthrds) \
    private(g, i, j, k, sc) \
    schedule(dynamic)
    for (g = 0; g < ngr; ++g) {
        const len_t gs  = g * gl;
        const len_t ge  = gs + gl < nrl ? gs + gl : nrl;
        const int tid   = omp_get_thread_num();
        len_t grp[GF_NR];
        len_t ng        = 0;
        uint8_t *dm     = NULL;
        if (gsz > 0) {
            sc  = ncols;
            for (i = gs; i < ge; ++i) {
                const hm_t * const r = upivs[i];
                const hm_t c0     = final ? first_column(r) : r[OFFSET];
                const len_t len   = r[LENGTH];
                if (c0 >= ncl) {
                    continue;
                }
                if (dm == NULL) {
                    if (drm[tid] == NULL) {
                        drm[tid] = gf2x_alloc_aligned(drm_raw + tid,
                                (uint64_t)ncols * GF_NR);
                    }
                    dm = drm[tid];
                }
                const cf8_t * const cfs = tbr->cf_8[r[COEFFS]];
                const hm_t * const ds   = r + OFFSET;
                for (j = 0; j < len; ++j) {
                    dm[(uint64_t)GF_NR * ds[j] + ng] = cfs[j];
                }
                if (c0 < sc) {
                    sc = c0;
                }
                grp[ng++] = i;
            }
            if (ng > 0) {
                reduce_dense_rows_by_known_pivots_gf2x(dm, mat, pivs, sc, t);
            }
        }
        uint8_t *drl  = dr + (uint64_t)tid * ncols;
        len_t lane    = 0;
        for (i = gs; i < ge; ++i) {
            hm_t *npiv      = upivs[i];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            if (lane < ng && grp[lane] == i) {
                /* the remainder of the row is in its lane, it only has
                 * entries in the columns of the new pivots */
                for (j = ncl; j < ncols; ++j) {
                    drl[j]  = dm[(uint64_t)GF_NR * j + lane];
                    dm[(uint64_t)GF_NR * j + lane] = 0;
                }
                lane++;
                sc  = ncl;
            } else {
                gf2x_load_row(drl, npiv + OFFSET, tbr->cf_8[npiv[COEFFS]],
                        npiv[LENGTH]);
                sc  = final ? first_column(npiv) : npiv[OFFSET];
            }
            free(npiv);
            upivs[i] = NULL;
            do {
                npiv = reduce_dense_row_by_pivots_gf2x(
                        drl, mat, pivs, sc, i, mh, bi, t);
                if (npiv == NULL) {
                    break;
                }
                /* normalize before the pivot becomes visible to the other
                 * threads */
                cf8_t *cf = mat->cf_8[i];
                normalize_sparse_row_gf2x(cf, npiv[LENGTH], t);
                k = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                if (k == 0) {
                    /* another thread added a pivot with the same lead
                     * term, reduce the row further */
                    gf2x_load_row(drl, npiv + OFFSET, cf, npiv[LENGTH]);
                    sc  = npiv[OFFSET];
                    free(npiv);
                    free(cf);
                    mat->cf_8[i] = NULL;
                }
            } while (!k);
        }
    }
    if (drm != NULL) {
        for (i = 0; i < (len_t)nthrds; ++i) {
            free(drm_raw[i]);
        }
        free(drm_raw);
        free(drm);
    }

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        if (kpivs == NULL) {
            free(pivs[i]);
        }
        pivs[i] = NULL;
    }
    free(kpivs);
    kpivs = NULL;

    len_t npivs = 0; /* number of new pivots */

    mat->tr = realloc(mat->tr, (uint64_t)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    if (!interreduce_new_pivots_dense_gf2x(mat, pivs, nthrds, &npivs, t)) {
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                hm_t *prev      = pivs[k];
                cf8_t *cfs      = mat->cf_8[prev[COEFFS]];
                const hm_t pos  = prev[COEFFS];
                gf2x_load_row(dr, prev + OFFSET, cfs, prev[LENGTH]);
                pivs[k] = NULL;
                pivs[k] = mat->tr[npivs++] = reduce_dense_row_by_pivots_gf2x(
                        dr, mat, pivs, k, pos, prev[MULT], prev[BINDEX], t);
                free(prev);
                free(cfs);
            }
        }
    }
    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;

    free(pivs);
    free(dr);
}

/* Reduces the rows mat->tr by the known pivots, returns the remainders as
 * dense rows of length ncr on the columns of the D part, NULL rows are
 * removed, mat->np is the number of remaining rows. */
static uint8_t **sparse_AB_CD_linear_algebra_gf2x(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    len_t i = 0, j;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncl   = mat->ncl;
    const len_t ncr   = mat->ncr;
    const int32_t nthrds  = st->nthrds;
    const gf2x_t * const t = gf2x_table(st->fc);

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    mat->cf_8     = realloc(mat->cf_8,
            (uint64_t)(nrl + mat->nru) * sizeof(cf8_t *));
    j = nrl;
    for (i = 0; i < mat->nru; ++i) {
        mat->cf_8[j]       = bs->cf_8[mat->rr[i][COEFFS]];
        mat->rr[i][COEFFS] = j;
        ++j;
    }
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    hm_t **upivs  = mat->tr;

    uint8_t **drs = (uint8_t **)calloc((uint64_t)nrl, sizeof(uint8_t *));

    const len_t gsz = rows_per_group_gf2x(nrl, nthrds);
    const len_t gl  = gsz > 0 ? gsz : 1;
    const len_t ngr = (nrl + gl - 1) / gl;
    /* the rows of a group are interleaved with GF_NR lanes */
    const len_t ls  = gsz > 0 ? GF_NR : 1;
    void *raw       = NULL;
    uint8_t *dr     = gf2x_alloc_aligned(&raw,
            (uint64_t)ncols * ls * nthrds);
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < ngr; ++i) {
        uint8_t *drl  = dr + (uint64_t)omp_get_thread_num() * ncols * ls;
        const len_t f = i * gl;
        const len_t l = f + gl < nrl ? gl : nrl - f;
        sc = ncols;
        for (len_t r = 0; r < l; ++r) {
            hm_t *npiv              = upivs[f+r];
            const cf8_t * const cfs = tbr->cf_8[npiv[COEFFS]];
            const len_t len         = npiv[LENGTH];
            const hm_t * const ds   = npiv + OFFSET;
            for (j = 0; j < len; ++j) {
                drl[(uint64_t)ls * ds[j] + r] = cfs[j];
            }
            if (ds[0] < sc) {
                sc = ds[0];
            }
            free(npiv);
            upivs[f+r] = NULL;
        }
        if (gsz > 0) {
            reduce_dense_rows_by_known_pivots_gf2x(drl, mat, pivs, sc, t);
        } else {
            /* only the known pivots are in pivs, so the reduction stops
             * at column ncl */
            for (len_t c = gf2x_next_nonzero(drl, sc, ncl); c < ncl;
                    c = gf2x_next_nonzero(drl, c+1, ncl)) {
                const hm_t * const piv = pivs[c];
                gf2x_scatter_axpy(drl, piv + OFFSET, mat->cf_8[piv[COEFFS]],
                        piv[LENGTH], drl[c], t);
            }
        }
        for (len_t r = 0; r < l; ++r) {
            uint8_t *row  = NULL;
            for (j = ncl; j < ncols; ++j) {
                const uint8_t v = drl[(uint64_t)ls * j + r];
                if (v != 0) {
                    if (row == NULL) {
                        row = (uint8_t *)calloc((uint64_t)ncr, sizeof(uint8_t));
                    }
                    row[j-ncl] = v;
                    drl[(uint64_t)ls * j + r] = 0;
                }
            }
            drs[f+r]  = row;
        }
    }
    free(raw);

    /* we do not need the old pivots anymore */
    free(kpivs);
    free(pivs);

    /* remove NULL dense rows */
    len_t npivs = 0;
    for (i = 0; i < nrl; ++i) {
        if (drs[i] != NULL) {
            drs[npivs++]  = drs[i];
        }
    }
    if (npivs == 0) {
        free(drs);
        drs = NULL;
    }
    st->np = mat->np = npivs;

    return drs;
}

/* Echelonizes the mat->np dense rows dm of length ncr. Returns the new
 * pivots: nps[c] is monic, starts at column c and has length ncr-c, or
 * is NULL. */
static uint8_t **exact_dense_linear_algebra_gf2x(
        uint8_t **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, npivs;

    const len_t nrows = mat->np;
    const len_t ncr   = mat->ncr;
    const gf2x_t * const t = gf2x_table(st->fc);

    uint8_t **nps = (uint8_t **)calloc((uint64_t)ncr, sizeof(uint8_t *));
    uint8_t **tbr = (uint8_t **)calloc((uint64_t)nrows, sizeof(uint8_t *));

    /* rows with a new lead term become pivots directly, the others have
     * to be reduced */
    j = 0;
    for (i = 0; i < nrows; ++i) {
        k = gf2x_next_nonzero(dm[i], 0, ncr);
        if (nps[k] == NULL) {
            memmove(dm[i], dm[i]+k, (uint64_t)(ncr-k) * sizeof(uint8_t));
            gf2x_mul_region(dm[i], dm[i], ncr-k, t->inv[dm[i][0]], t, 0);
            nps[k] = dm[i];
        } else {
            tbr[j++]  = dm[i];
        }
    }
    free(dm);
    const len_t ntr = j;

#pragma omp parallel for num_threads(st->nthrds) \
    private(i, k) schedule(dynamic)
    for (i = 0; i < ntr; ++i) {
        uint8_t *dr = tbr[i];
        len_t c = 0;
        do {
            c = reduce_dense_row_by_dense_pivots_gf2x(dr, c, nps, ncr, t);
            if (c == ncr) {
                break;
            }
            uint8_t *piv = (uint8_t *)malloc((uint64_t)(ncr-c) * sizeof(uint8_t));
            gf2x_mul_region(piv, dr+c, ncr-c, t->inv[dr[c]], t, 0);
            k = __sync_bool_compare_and_swap(&nps[c], NULL, piv);
            if (k == 0) {
                /* another thread was faster, continue reducing */
                free(piv);
            }
        } while (!k);
        free(dr);
    }
    free(tbr);

    npivs = 0;
    for (i = 0; i < ncr; ++i) {
        npivs += nps[i] != NULL;
    }
    st->np = mat->np = npivs;

    return nps;
}

/* converts the dense pivots dm, see exact_dense_linear_algebra_gf2x(),
 * to sparse rows in mat->tr from the right to the left */
static void convert_to_sparse_matrix_rows_gf2x(
        mat_t *mat,
        uint8_t * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, j, k, l, m;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->tr   = realloc(mat->tr, (uint64_t)mat->np * sizeof(hm_t *));
    mat->cf_8 = realloc(mat->cf_8, (uint64_t)mat->np * sizeof(cf8_t *));

    l = 0;
    for (i = 0; i < ncr; ++i) {
        m = ncr-1-i;
        if (dm[m] != NULL) {
            const len_t len = ncr-m;
            const uint8_t * const d = dm[m];
            k = 0;
            for (j = 0; j < len; ++j) {
                k += d[j] != 0;
            }
            cf8_t *cfs = (cf8_t *)malloc((uint64_t)k * sizeof(cf8_t));
            hm_t *dts  = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
            hm_t *dss  = dts + OFFSET;
            k = 0;
            for (j = gf2x_next_nonzero(d, 0, len); j < len;
                    j = gf2x_next_nonzero(d, j+1, len)) {
                cfs[k]  = d[j];
                dss[k]  = ncl+m+j;
                k++;
            }
            dts[BINDEX]   = 0;
            dts[MULT]     = 0;
            dts[COEFFS]   = l;
            dts[PRELOOP]  = k % UNROLL;
            dts[LENGTH]   = k;
            mat->tr[l]    = dts;
            mat->cf_8[l]  = cfs;
            l++;
        }
    }
}

static void exact_sparse_linear_algebra_gf2x(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* if the rows to be reduced are dense on the columns of the D part,
     * their remainders are echelonized as dense rows */
    uint64_t nnz = 0;
    for (len_t i = 0; i < mat->nrl; ++i) {
        nnz += mat->tr[i][LENGTH];
    }
    if (st->in_final_reduction_step == 0 && mat->nrl > 0
            && 2 * nnz >= (uint64_t)mat->nrl * mat->ncr) {
        uint8_t **dm = sparse_AB_CD_linear_algebra_gf2x(mat, tbr, bs, st);
        if (mat->np > 0) {
            dm  = exact_dense_linear_algebra_gf2x(dm, mat, st);
            interreduce_dense_matrix_gf2x(dm, mat->ncr, st->nthrds,
                    gf2x_table(st->fc));
        }
        convert_to_sparse_matrix_rows_gf2x(mat, dm);
        if (dm) {
            for (len_t i = 0; i < mat->ncr; ++i) {
                free(dm[i]);
            }
            free(dm);
        }
        st->np = mat->nr = mat->sz = mat->np;
    } else {
        mat->cf_8 = realloc(mat->cf_8, (uint64_t)mat->nr * sizeof(cf8_t *));
        exact_sparse_reduced_echelon_form_gf2x(mat, tbr, bs, st);
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "%9d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(VERBSTREAM);
    }
}

static void interreduce_matrix_rows_gf2x(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        const int free_basis
        )
{
    len_t i, k, l;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;
    const gf2x_t * const t = gf2x_table(st->fc);

    /* adjust displaying timings for statistic printout */
    if (st->info_level > 1) {
        fprintf(VERBSTREAM, "                          ");
    }

    mat->tr = realloc(mat->tr, (uint64_t)ncols * sizeof(hm_t *));

    mat->cf_8 = realloc(mat->cf_8, (uint64_t)ncols * sizeof(cf8_t *));
    memset(mat->cf_8, 0, (uint64_t)ncols * sizeof(cf8_t *));
    hm_t **pivs = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    for (i = 0; i < nrows; ++i) {
        pivs[mat->rr[i][OFFSET]]  = mat->rr[i];
    }

    uint8_t *dr = (uint8_t *)calloc((uint64_t)ncols, sizeof(uint8_t));
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            hm_t *prev  = pivs[l];
            gf2x_load_row(dr, prev + OFFSET, bs->cf_8[prev[COEFFS]],
                    prev[LENGTH]);
            pivs[l] = NULL;
            pivs[l] = mat->tr[k--] = reduce_dense_row_by_pivots_gf2x(
                    dr, mat, pivs, l, l, prev[MULT], prev[BINDEX], t);
            free(prev);
        }
    }
    if (free_basis != 0) {
        /* free now all polynomials in the basis and reset bs->ld to 0. */
        free_basis_elements(bs);
    }
    free(mat->rr);
    mat->rr = NULL;
    mat->np = nrows;
    free(pivs);
    free(dr);
}

static void normalize_initial_basis_gf2x(
        bs_t *bs,
        const uint64_t fc
        )
{
    const gf2x_t * const t = gf2x_table(fc);
    for (len_t i = 0; i < bs->ld; ++i) {
        normalize_sparse_row_gf2x(bs->cf_8[bs->hm[i][COEFFS]],
                bs->hm[i][LENGTH], t);
    }
}
