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

/* Linear algebra over GF(2). All nonzero coefficients are 1, they are
 * still stored as cf8_t in the arrays cf_8 like for the other primes
 * < 2^8, so the matrix and basis handling is unchanged, but the linear
 * algebra only works on the columns of the rows. Dense rows are bit
 * vectors, adding a row is XOR.
 *
 * The rows reduced by the known pivots are handled in groups of up to
 * GF2_NR rows which are bit sliced: the GF2_LW words of column c hold
 * the entries of all rows of the group in column c, one bit per row.
 * Reducing all rows of the group by a pivot then costs one vector XOR
 * per term of the pivot. Only exact linear algebra for Groebner basis
 * computations and normal forms is supported, no tracer. */

#if defined HAVE_AVX512_F || defined HAVE_AVX2 || defined HAVE_SSE2
#include <immintrin.h>
#elif defined __aarch64__
#include <arm_neon.h>
#endif

/* number of 64 bit words per column of a group of rows */
#if defined HAVE_AVX512_F
#define GF2_LW 8
#elif defined HAVE_AVX2
#define GF2_LW 4
#elif defined HAVE_SSE2 || defined __aarch64__
#define GF2_LW 2
#else
#define GF2_LW 1
#endif
/* maximal number of rows of a group */
#define GF2_NR (64*GF2_LW)
/* maximal memory for the dense rows of the D part in bytes */
#define GF2_MAX_DENSE_D ((uint64_t)1 << 30)

/* returns zeroed memory of nw words aligned to 64 bytes, *raw is the
 * pointer to be freed */
static uint64_t *gf2_alloc_aligned(
        void **raw,
        const uint64_t nw
        )
{
    *raw = calloc(nw * sizeof(uint64_t) + 64, 1);
    return (uint64_t *)(((uintptr_t)*raw + 63) & ~(uintptr_t)63);
}

/* d ^= s for the GF2_LW words of one column of a group, both aligned */
static inline void gf2_lanes_xor(
        uint64_t *d,
        const uint64_t * const s
        )
{
#if GF2_LW == 8
    _mm512_store_si512((void *)d, _mm512_xor_si512(
                _mm512_load_si512((const void *)d),
                _mm512_load_si512((const void *)s)));
#elif GF2_LW == 4
    _mm256_store_si256((__m256i *)d, _mm256_xor_si256(
                _mm256_load_si256((const __m256i *)d),
                _mm256_load_si256((const __m256i *)s)));
#elif GF2_LW == 2 && defined __aarch64__
    vst1q_u64(d, veorq_u64(vld1q_u64(d), vld1q_u64(s)));
#elif GF2_LW == 2
    _mm_store_si128((__m128i *)d, _mm_xor_si128(
                _mm_load_si128((const __m128i *)d),
                _mm_load_si128((const __m128i *)s)));
#else
    d[0] ^= s[0];
#endif
}

/* returns 1 if the GF2_LW words of one column of a group are zero */
static inline int gf2_lanes_zero(
        const uint64_t * const d
        )
{
#if GF2_LW == 8
    const __m512i v = _mm512_load_si512((const void *)d);
    return _mm512_test_epi64_mask(v, v) == 0;
#elif GF2_LW == 4
    const __m256i v = _mm256_load_si256((const __m256i *)d);
    return _mm256_testz_si256(v, v);
#elif GF2_LW == 2
    return (d[0] | d[1]) == 0;
#else
    return d[0] == 0;
#endif
}

/* d ^= s for nw words, no alignment needed */
static inline void gf2_xor_region(
        uint64_t *d,
        const uint64_t * const s,
        const len_t nw
        )
{
    len_t i = 0;
#if defined HAVE_AVX512_F
    for (; i + 8 <= nw; i += 8) {
        _mm512_storeu_si512((void *)(d+i), _mm512_xor_si512(
                    _mm512_loadu_si512((const void *)(d+i)),
                    _mm512_loadu_si512((const void *)(s+i))));
    }
#endif
#if defined HAVE_AVX2
    for (; i + 4 <= nw; i += 4) {
        _mm256_storeu_si256((__m256i *)(d+i), _mm256_xor_si256(
                    _mm256_loadu_si256((const __m256i *)(d+i)),
                    _mm256_loadu_si256((const __m256i *)(s+i))));
    }
#endif
#if defined __aarch64__
    for (; i + 2 <= nw; i += 2) {
        vst1q_u64(d+i, veorq_u64(vld1q_u64(d+i), vld1q_u64(s+i)));
    }
#elif defined HAVE_SSE2
    for (; i + 2 <= nw; i += 2) {
        _mm_storeu_si128((__m128i *)(d+i), _mm_xor_si128(
                    _mm_loadu_si128((const __m128i *)(d+i)),
                    _mm_loadu_si128((const __m128i *)(s+i))));
    }
#endif
    for (; i < nw; ++i) {
        d[i] ^= s[i];
    }
}

/* returns the first set bit of the bit vector d at position >= i, n if
 * there is none; bits at positions >= n are zero */
static inline len_t gf2_next_bit(
        const uint64_t * const d,
        const len_t i,
        const len_t n
        )
{
    if (i >= n) {
        return n;
    }
    len_t w        = i >> 6;
    const len_t nw = (n + 63) >> 6;
    uint64_t x     = d[w] & (~(uint64_t)0 << (i & 63));
    while (x == 0) {
        if (++w == nw) {
            return n;
        }
        x = d[w];
    }
    return (w << 6) + (len_t)__builtin_ctzll(x);
}

static inline int gf2_get_bit(
        const uint64_t * const d,
        const len_t i
        )
{
    return (d[i >> 6] >> (i & 63)) & 1;
}

/* flips the bits of the columns ds[0],...,ds[len-1] in d */
static inline void gf2_flip_bits(
        uint64_t *d,
        const hm_t * const ds,
        const len_t len
        )
{
    len_t j = 0;
    const len_t os = len % UNROLL;
    for (; j < os; ++j) {
        d[ds[j] >> 6] ^= (uint64_t)1 << (ds[j] & 63);
    }
    for (; j < len; j += UNROLL) {
        d[ds[j] >> 6]   ^= (uint64_t)1 << (ds[j] & 63);
        d[ds[j+1] >> 6] ^= (uint64_t)1 << (ds[j+1] & 63);
        d[ds[j+2] >> 6] ^= (uint64_t)1 << (ds[j+2] & 63);
        d[ds[j+3] >> 6] ^= (uint64_t)1 << (ds[j+3] & 63);
    }
}

/* returns a coefficient array of len ones */
static inline cf8_t *gf2_ones(
        const len_t len
        )
{
    cf8_t *cf = (cf8_t *)malloc((uint64_t)(len > 0 ? len : 1) * sizeof(cf8_t));
    memset(cf, 1, (uint64_t)len * sizeof(cf8_t));
    return cf;
}

/* Reduces the bit vector dr by the pivots pivs starting at column sc.
 * The remainder is returned as sparse row, its coefficients are stored in
 * mat->cf_8[tmp_pos], NULL is returned if the remainder is zero.
 * Afterwards dr is zero. */
static hm_t *reduce_dense_row_by_pivots_gf2(
        uint64_t *dr,
        mat_t *mat,
        hm_t * const * const pivs,
        const hi_t sc,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi
        )
{
    len_t i, j, k = 0;
    hi_t np = mat->nc;

    const len_t ncols = mat->nc;

    for (i = gf2_next_bit(dr, sc, ncols); i < ncols;
            i = gf2_next_bit(dr, i+1, ncols)) {
        const hm_t * const piv = pivs[i];
        if (piv == NULL) {
            if (np == ncols) {
                np = i;
            }
            k++;
            continue;
        }
        /* this also clears the bit of column i */
        gf2_flip_bits(dr, piv + OFFSET, piv[LENGTH]);
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    hm_t *rs  = row + OFFSET;
    j = 0;
    for (i = gf2_next_bit(dr, np, ncols); i < ncols;
            i = gf2_next_bit(dr, i+1, ncols)) {
        rs[j++] = (hm_t)i;
        dr[i >> 6] ^= (uint64_t)1 << (i & 63);
    }
    row[BINDEX]   = bi;
    row[MULT]     = mh;
    row[COEFFS]   = tmp_pos;
    row[PRELOOP]  = j % UNROLL;
    row[LENGTH]   = j;
    mat->cf_8[tmp_pos] = gf2_ones(j);

    return row;
}

/* Reduces the rows of the group dm, the GF2_LW words starting at
 * dm[GF2_LW*c] hold the entries of the rows in column c, by the known
 * pivots pivs[sc],...,pivs[ncl-1]. Afterwards the rows only have entries
 * in the columns ncl,...,ncols-1. */
static void reduce_dense_rows_by_known_pivots_gf2(
        uint64_t *dm,
        const mat_t * const mat,
        hm_t * const * const pivs,
        const len_t sc
        )
{
    len_t i, j;
    uint64_t m[GF2_LW] __attribute__((aligned(64)));

    const len_t ncl = mat->ncl;

    for (i = sc; i < ncl; ++i) {
        uint64_t *d = dm + (uint64_t)GF2_LW * i;
        if (gf2_lanes_zero(d)) {
            continue;
        }
        /* the rows with an entry in column i are reduced by the pivot,
         * its lead term clears the column */
        memcpy(m, d, sizeof(m));
        memset(d, 0, sizeof(m));
        const hm_t * const piv = pivs[i];
        const hm_t * const ds  = piv + OFFSET;
        const len_t len        = piv[LENGTH];
        for (j = 1; j < len; ++j) {
            gf2_lanes_xor(dm + (uint64_t)GF2_LW * ds[j], m);
        }
    }
}

/* Extracts the ng rows of the group dm from the columns ncl,...,ncols-1
 * as lists of columns: rcs[r] are the rcl[r] columns of row r. The group
 * is zero afterwards. */
static void extract_rows_from_group_gf2(
        uint64_t *dm,
        const len_t ncl,
        const len_t ncols,
        const len_t ng,
        hm_t **rcs,
        len_t *rcl
        )
{
    len_t j, w;

    memset(rcl, 0, (uint64_t)ng * sizeof(len_t));
    for (j = ncl; j < ncols; ++j) {
        const uint64_t * const d = dm + (uint64_t)GF2_LW * j;
        for (w = 0; w < GF2_LW; ++w) {
            uint64_t x = d[w];
            while (x != 0) {
                rcl[64*w + __builtin_ctzll(x)]++;
                x &= x - 1;
            }
        }
    }
    for (j = 0; j < ng; ++j) {
        rcs[j] = rcl[j] > 0 ?
            (hm_t *)malloc((uint64_t)rcl[j] * sizeof(hm_t)) : NULL;
        rcl[j] = 0;
    }
    for (j = ncl; j < ncols; ++j) {
        uint64_t *d = dm + (uint64_t)GF2_LW * j;
        for (w = 0; w < GF2_LW; ++w) {
            uint64_t x = d[w];
            while (x != 0) {
                const len_t r = 64*w + __builtin_ctzll(x);
                rcs[r][rcl[r]++] = j;
                x &= x - 1;
            }
            d[w] = 0;
        }
    }
}

/* Reduces the dense row dr of ncr bits starting at column sc by the
 * dense pivots nps, nps[c] has its lead term in column c. Stops at the
 * first column without pivot and returns it, ncr if dr becomes zero. */
static len_t reduce_dense_row_by_dense_pivots_gf2(
        uint64_t *dr,
        const len_t sc,
        uint64_t * const * const nps,
        const len_t ncr
        )
{
    len_t i;
    const len_t nw = (ncr + 63) >> 6;
    for (i = gf2_next_bit(dr, sc, ncr); i < ncr;
            i = gf2_next_bit(dr, i+1, ncr)) {
        const uint64_t *piv = nps[i];
        if (piv == NULL) {
            return i;
        }
        gf2_xor_region(dr + (i >> 6), piv + (i >> 6), nw - (i >> 6));
    }
    return ncr;
}

/* Interreduces the dense pivots dm of ncr bits, dm[c] has its lead term
 * in column c or is NULL. The pivots are handled in blocks from the right
 * to the left: the rows of a block are interreduced, afterwards all rows
 * left of the block are reduced by the block in parallel. */
static void interreduce_dense_matrix_gf2(
        uint64_t **dm,
        const len_t ncr,
        const int32_t nthrds
        )
{
    len_t i, n = 0;
    int64_t r;
    const len_t bs = 64;
    const len_t nw = (ncr + 63) >> 6;

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
            uint64_t *row = dm[pc[r]];
            for (i = r+1; i < e; ++i) {
                const len_t c = pc[i];
                if (gf2_get_bit(row, c)) {
                    gf2_xor_region(row + (c >> 6), dm[c] + (c >> 6),
                            nw - (c >> 6));
                }
            }
        }
        /* reduce all rows left of the block by the block */
#pragma omp parallel for num_threads(nthrds) private(i) schedule(dynamic, 8)
        for (r = 0; r < (int64_t)s; ++r) {
            uint64_t *row = dm[pc[r]];
            for (i = s; i < e; ++i) {
                const len_t c = pc[i];
                if (gf2_get_bit(row, c)) {
                    gf2_xor_region(row + (c >> 6), dm[c] + (c >> 6),
                            nw - (c >> 6));
                }
            }
        }
        e = s;
    }
    free(pc);
}

/* Converts the dense row d of ncr bits, starting at column sc, to a
 * sparse row with header hdr whose columns are shifted by ncl. Its
 * coefficients are stored in mat->cf_8[pos]. */
static hm_t *dense_to_sparse_row_gf2(
        mat_t *mat,
        const uint64_t * const d,
        const len_t sc,
        const len_t ncr,
        const len_t ncl,
        const hm_t * const hdr,
        const hm_t pos
        )
{
    len_t j, k = 0;
    const len_t nw = (ncr + 63) >> 6;
    for (j = sc >> 6; j < nw; ++j) {
        k += (len_t)__builtin_popcountll(d[j]);
    }
    hm_t *row = (hm_t *)malloc((uint64_t)(k+OFFSET) * sizeof(hm_t));
    if (hdr != NULL) {
        memcpy(row, hdr, OFFSET * sizeof(hm_t));
    } else {
        row[BINDEX] = 0;
        row[MULT]   = 0;
    }
    hm_t *rs = row + OFFSET;
    k = 0;
    for (j = gf2_next_bit(d, sc, ncr); j < ncr; j = gf2_next_bit(d, j+1, ncr)) {
        rs[k++] = ncl + j;
    }
    row[COEFFS]   = pos;
    row[PRELOOP]  = k % UNROLL;
    row[LENGTH]   = k;
    mat->cf_8[pos] = gf2_ones(k);
    return row;
}

/* Interreduces the new pivots pivs[ncl],...,pivs[ncl+ncr-1] as dense rows
 * if they are dense enough and stores them in mat->tr from the right to
 * the left. Returns 0 if the pivots are too sparse, they are not changed
 * then. */
static int interreduce_new_pivots_dense_gf2(
        mat_t *mat,
        hm_t **pivs,
        const int32_t nthrds,
        len_t *npivsp
        )
{
    len_t i, j;
    uint64_t nnz = 0, dsz = 0;
    len_t n = 0;

    const len_t ncl = mat->ncl;
    const len_t ncr = mat->ncr;
    const len_t nw  = (ncr + 63) >> 6;

    for (i = 0; i < ncr; ++i) {
        if (pivs[ncl+i] != NULL) {
            n++;
            nnz += pivs[ncl+i][LENGTH];
            dsz += ncr - i;
        }
    }
    /* a dense row costs ncr/64 words instead of nnz columns */
    if (n < 2 || 64 * nnz < dsz) {
        return 0;
    }
    uint64_t **dm = (uint64_t **)calloc((uint64_t)ncr, sizeof(uint64_t *));
    hm_t *hdr     = (hm_t *)malloc((uint64_t)ncr * OFFSET * sizeof(hm_t));
    for (i = 0; i < ncr; ++i) {
        hm_t *row = pivs[ncl+i];
        if (row == NULL) {
            continue;
        }
        memcpy(hdr + (uint64_t)i * OFFSET, row, OFFSET * sizeof(hm_t));
        dm[i] = (uint64_t *)calloc((uint64_t)nw, sizeof(uint64_t));
        const hm_t * const ds = row + OFFSET;
        for (j = 0; j < row[LENGTH]; ++j) {
            const len_t c = ds[j] - ncl;
            dm[i][c >> 6] |= (uint64_t)1 << (c & 63);
        }
        free(mat->cf_8[row[COEFFS]]);
        mat->cf_8[row[COEFFS]] = NULL;
        free(row);
        pivs[ncl+i] = NULL;
    }
    interreduce_dense_matrix_gf2(dm, ncr, nthrds);

    len_t npivs = 0;
    for (i = ncr; i > 0; --i) {
        uint64_t *d = dm[i-1];
        if (d == NULL) {
            continue;
        }
        const hm_t *h = hdr + (uint64_t)(i-1) * OFFSET;
        pivs[ncl+i-1] = mat->tr[npivs++] = dense_to_sparse_row_gf2(
                mat, d, i-1, ncr, ncl, h, h[COEFFS]);
        free(d);
    }
    free(dm);
    free(hdr);
    *npivsp = npivs;
    return 1;
}

/* Number of rows reduced together by the known pivots, 0 if the rows
 * shall be reduced one by one. A group costs as much as GF2_NR rows, but
 * applying a pivot to all rows of a group is much cheaper than applying
 * it to each row on its own, even if only some rows need it. Thus we use
 * smaller groups if there are not enough rows for each thread. */
static inline len_t rows_per_group_gf2(
        const len_t nrl,
        const int32_t nthrds
        )
{
    const uint64_t gl = ((uint64_t)nrl + nthrds - 1) / nthrds;
    if (gl >= GF2_NR) {
        return GF2_NR;
    }
    return gl >= 4 ? (len_t)gl : 0;
}

/* loads the columns ds[0],...,ds[len-1] into the bit vector dr which is
 * zero */
static inline void gf2_load_row(
        uint64_t *dr,
        const hm_t * const ds,
        const len_t len
        )
{
    for (len_t j = 0; j < len; ++j) {
        dr[ds[j] >> 6] |= (uint64_t)1 << (ds[j] & 63);
    }
}

static void exact_sparse_reduced_echelon_form_gf2(
        mat_t *mat,
        md_t *st
        )
{
    len_t i = 0, j, k;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;
    const len_t nw    = (ncols + 63) >> 6;

    const int32_t nthrds  = st->nthrds;
    const int final       = st->in_final_reduction_step;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    hm_t *kpivs   = reducer_rows_to_block(mat);
    if (mat->nru > 0) {
        memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));
    }

    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = mat->tr;

    /* the dense rows are zero whenever a row is started */
    void *dr_raw  = NULL;
    uint64_t *dr  = gf2_alloc_aligned(&dr_raw, (uint64_t)nw * nthrds);

    /* the rows are first reduced by the known pivots in groups, see
     * reduce_dense_rows_by_known_pivots_gf2() */
    const len_t gsz = rows_per_group_gf2(nrl, nthrds);
    void **drm_raw  = NULL;
    uint64_t **drm  = NULL;
    if (gsz > 0) {
        drm_raw = (void **)calloc((uint64_t)nthrds, sizeof(void *));
        drm     = (uint64_t **)calloc((uint64_t)nthrds, sizeof(uint64_t *));
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
        len_t grp[GF2_NR];
        hm_t *rcs[GF2_NR];
        len_t rcl[GF2_NR];
        len_t ng        = 0;
        uint64_t *dm    = NULL;
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
                        drm[tid] = gf2_alloc_aligned(drm_raw + tid,
                                (uint64_t)ncols * GF2_LW);
                    }
                    dm = drm[tid];
                }
                const hm_t * const ds   = r + OFFSET;
                const uint64_t b        = (uint64_t)1 << (ng & 63);
                for (j = 0; j < len; ++j) {
                    dm[(uint64_t)GF2_LW * ds[j] + (ng >> 6)] |= b;
                }
                if (c0 < sc) {
                    sc = c0;
                }
                grp[ng++] = i;
            }
            if (ng > 0) {
                reduce_dense_rows_by_known_pivots_gf2(dm, mat, pivs, sc);
                extract_rows_from_group_gf2(dm, ncl, ncols, ng, rcs, rcl);
            }
        }
        uint64_t *drl = dr + (uint64_t)tid * nw;
        len_t lane    = 0;
        for (i = gs; i < ge; ++i) {
            hm_t *npiv      = upivs[i];
            const len_t bi  = npiv[BINDEX];
            const len_t mh  = npiv[MULT];
            if (lane < ng && grp[lane] == i) {
                /* the remainder of the row only has entries in the
                 * columns of the new pivots */
                gf2_load_row(drl, rcs[lane], rcl[lane]);
                free(rcs[lane]);
                lane++;
                sc  = ncl;
            } else {
                gf2_load_row(drl, npiv + OFFSET, npiv[LENGTH]);
                sc  = final ? first_column(npiv) : npiv[OFFSET];
            }
            free(npiv);
            upivs[i] = NULL;
            do {
                npiv = reduce_dense_row_by_pivots_gf2(
                        drl, mat, pivs, sc, i, mh, bi);
                if (npiv == NULL) {
                    break;
                }
                k = __sync_bool_compare_and_swap(&pivs[npiv[OFFSET]], NULL, npiv);
                if (k == 0) {
                    /* another thread added a pivot with the same lead
                     * term, reduce the row further */
                    gf2_load_row(drl, npiv + OFFSET, npiv[LENGTH]);
                    sc  = npiv[OFFSET];
                    free(npiv);
                    free(mat->cf_8[i]);
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
    if (!interreduce_new_pivots_dense_gf2(mat, pivs, nthrds, &npivs)) {
        for (i = 0; i < ncr; ++i) {
            k = ncols-1-i;
            if (pivs[k]) {
                hm_t *prev      = pivs[k];
                const hm_t pos  = prev[COEFFS];
                gf2_load_row(dr, prev + OFFSET, prev[LENGTH]);
                free(mat->cf_8[pos]);
                mat->cf_8[pos] = NULL;
                pivs[k] = NULL;
                pivs[k] = mat->tr[npivs++] = reduce_dense_row_by_pivots_gf2(
                        dr, mat, pivs, k, pos, prev[MULT], prev[BINDEX]);
                free(prev);
            }
        }
    }
    mat->tr = realloc(mat->tr, (uint64_t)npivs * sizeof(hi_t *));
    st->np = mat->np = mat->nr = mat->sz = npivs;

    free(pivs);
    free(dr_raw);
}

/* Reduces the rows mat->tr by the known pivots, returns the remainders as
 * dense rows of ncr bits on the columns of the D part, NULL rows are
 * removed, mat->np is the number of remaining rows. */
static uint64_t **sparse_AB_CD_linear_algebra_gf2(
        mat_t *mat,
        md_t *st
        )
{
    len_t i = 0, j;
    hi_t sc = 0;    /* starting column */

    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t ncl   = mat->ncl;
    const len_t ncr   = mat->ncr;
    const len_t nwr   = (ncr + 63) >> 6;
    const len_t nw    = (ncols + 63) >> 6;
    const int32_t nthrds  = st->nthrds;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((uint64_t)ncols, sizeof(hm_t *));
    hm_t *kpivs   = reducer_rows_to_block(mat);
    memcpy(pivs, mat->rr, (uint64_t)mat->nru * sizeof(hm_t *));

    hm_t **upivs  = mat->tr;

    uint64_t **drs = (uint64_t **)calloc((uint64_t)nrl, sizeof(uint64_t *));

    const len_t gsz = rows_per_group_gf2(nrl, nthrds);
    const len_t gl  = gsz > 0 ? gsz : 1;
    const len_t ngr = (nrl + gl - 1) / gl;
    /* a group needs GF2_LW words per column, a single row a bit vector */
    const uint64_t tsz = gsz > 0 ? (uint64_t)ncols * GF2_LW : nw;
    void *raw       = NULL;
    uint64_t *dr    = gf2_alloc_aligned(&raw, tsz * nthrds);
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, sc) \
    schedule(dynamic)
    for (i = 0; i < ngr; ++i) {
        uint64_t *drl = dr + (uint64_t)omp_get_thread_num() * tsz;
        const len_t f = i * gl;
        const len_t l = f + gl < nrl ? gl : nrl - f;
        if (gsz > 0) {
            hm_t *rcs[GF2_NR];
            len_t rcl[GF2_NR];
            sc = ncols;
            for (len_t r = 0; r < l; ++r) {
                hm_t *npiv              = upivs[f+r];
                const len_t len         = npiv[LENGTH];
                const hm_t * const ds   = npiv + OFFSET;
                const uint64_t b        = (uint64_t)1 << (r & 63);
                for (j = 0; j < len; ++j) {
                    drl[(uint64_t)GF2_LW * ds[j] + (r >> 6)] |= b;
                }
                if (ds[0] < sc) {
                    sc = ds[0];
                }
                free(npiv);
                upivs[f+r] = NULL;
            }
            reduce_dense_rows_by_known_pivots_gf2(drl, mat, pivs, sc);
            extract_rows_from_group_gf2(drl, ncl, ncols, l, rcs, rcl);
            for (len_t r = 0; r < l; ++r) {
                uint64_t *row = NULL;
                if (rcl[r] > 0) {
                    row = (uint64_t *)calloc((uint64_t)nwr, sizeof(uint64_t));
                    for (j = 0; j < rcl[r]; ++j) {
                        const len_t c = rcs[r][j] - ncl;
                        row[c >> 6] |= (uint64_t)1 << (c & 63);
                    }
                }
                free(rcs[r]);
                drs[f+r]  = row;
            }
        } else {
            hm_t *npiv  = upivs[f];
            gf2_load_row(drl, npiv + OFFSET, npiv[LENGTH]);
            sc  = npiv[OFFSET];
            free(npiv);
            upivs[f] = NULL;
            /* only the known pivots are in pivs, so the reduction stops
             * at column ncl */
            for (len_t c = gf2_next_bit(drl, sc, ncl); c < ncl;
                    c = gf2_next_bit(drl, c+1, ncl)) {
                const hm_t * const piv = pivs[c];
                gf2_flip_bits(drl, piv + OFFSET, piv[LENGTH]);
            }
            uint64_t *row = NULL;
            for (len_t c = gf2_next_bit(drl, ncl, ncols); c < ncols;
                    c = gf2_next_bit(drl, c+1, ncols)) {
                if (row == NULL) {
                    row = (uint64_t *)calloc((uint64_t)nwr, sizeof(uint64_t));
                }
                row[(c-ncl) >> 6] |= (uint64_t)1 << ((c-ncl) & 63);
                drl[c >> 6] ^= (uint64_t)1 << (c & 63);
            }
            drs[f]  = row;
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

/* Echelonizes the mat->np dense rows dm of ncr bits. Returns the new
 * pivots: nps[c] has its lead term in column c or is NULL. */
static uint64_t **exact_dense_linear_algebra_gf2(
        uint64_t **dm,
        mat_t *mat,
        md_t *st
        )
{
    len_t i, j, k, npivs;

    const len_t nrows = mat->np;
    const len_t ncr   = mat->ncr;

    uint64_t **nps = (uint64_t **)calloc((uint64_t)ncr, sizeof(uint64_t *));
    uint64_t **tbr = (uint64_t **)calloc((uint64_t)nrows, sizeof(uint64_t *));

    /* rows with a new lead term become pivots directly, the others have
     * to be reduced */
    j = 0;
    for (i = 0; i < nrows; ++i) {
        k = gf2_next_bit(dm[i], 0, ncr);
        if (nps[k] == NULL) {
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
        uint64_t *dr = tbr[i];
        len_t c = 0;
        do {
            c = reduce_dense_row_by_dense_pivots_gf2(dr, c, nps, ncr);
            if (c == ncr) {
                free(dr);
                break;
            }
            k = __sync_bool_compare_and_swap(&nps[c], NULL, dr);
        } while (!k);
    }
    free(tbr);

    npivs = 0;
    for (i = 0; i < ncr; ++i) {
        npivs += nps[i] != NULL;
    }
    st->np = mat->np = npivs;

    return nps;
}

/* converts the dense pivots dm, see exact_dense_linear_algebra_gf2(), to
 * sparse rows in mat->tr from the right to the left */
static void convert_to_sparse_matrix_rows_gf2(
        mat_t *mat,
        uint64_t * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, l = 0;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->tr   = realloc(mat->tr, (uint64_t)mat->np * sizeof(hm_t *));
    mat->cf_8 = realloc(mat->cf_8, (uint64_t)mat->np * sizeof(cf8_t *));

    for (i = ncr; i > 0; --i) {
        if (dm[i-1] != NULL) {
            mat->tr[l] = dense_to_sparse_row_gf2(mat, dm[i-1], i-1, ncr, ncl,
                    NULL, l);
            l++;
        }
    }
}

static void exact_sparse_linear_algebra_gf2(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    (void)tbr;
    (void)bs;
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* the remainders of the rows reduced by the known pivots only have
     * entries in the columns of the D part, as bit vectors they are so
     * small that echelonizing them as dense rows is much faster than
     * reducing them one by one by the sparse new pivots. we only take care
     * that the dense rows do not need too much memory. */
    const uint64_t dsz = (uint64_t)mat->nrl * ((mat->ncr + 63) >> 6)
        * sizeof(uint64_t);
    if (st->in_final_reduction_step == 0 && mat->nrl > 0
            && dsz <= GF2_MAX_DENSE_D) {
        uint64_t **dm = sparse_AB_CD_linear_algebra_gf2(mat, st);
        if (mat->np > 0) {
            dm  = exact_dense_linear_algebra_gf2(dm, mat, st);
            interreduce_dense_matrix_gf2(dm, mat->ncr, st->nthrds);
        }
        convert_to_sparse_matrix_rows_gf2(mat, dm);
        if (dm) {
            for (len_t i = 0; i < mat->ncr; ++i) {
                free(dm[i]);
            }
            free(dm);
        }
        st->np = mat->nr = mat->sz = mat->np;
    } else {
        mat->cf_8 = realloc(mat->cf_8, (uint64_t)mat->nr * sizeof(cf8_t *));
        memset(mat->cf_8, 0, (uint64_t)mat->nr * sizeof(cf8_t *));
        exact_sparse_reduced_echelon_form_gf2(mat, st);
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

static void interreduce_matrix_rows_gf2(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        const int free_basis
        )
{
    len_t i, k, l;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;

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

    void *raw     = NULL;
    uint64_t *dr  = gf2_alloc_aligned(&raw, (uint64_t)(ncols + 63) >> 6);
    k = nrows - 1;
    for (i = 0; i < ncols; ++i) {
        l = ncols-1-i;
        if (pivs[l] != NULL) {
            hm_t *prev  = pivs[l];
            gf2_load_row(dr, prev + OFFSET, prev[LENGTH]);
            pivs[l] = NULL;
            pivs[l] = mat->tr[k--] = reduce_dense_row_by_pivots_gf2(
                    dr, mat, pivs, l, l, prev[MULT], prev[BINDEX]);
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
    free(raw);
}
