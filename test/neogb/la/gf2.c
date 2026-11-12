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

/* Checks the vectorized bit operations of the linear algebra over GF(2)
 * in la_gf2.c, compiled with the same SIMD flags as the library, against
 * scalar computations. */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../src/neogb/data.h"
#include "../../../src/neogb/tools.h"
#include "../../../src/msolve/streams.h"

/* stubs for the functions of the library la_gf2.c needs which are not
 * exported */
static void free_basis_elements(bs_t *bs)
{
    (void)bs;
}

static hm_t *reducer_rows_to_block(mat_t *mat)
{
    (void)mat;
    return NULL;
}

static hm_t first_column(const hm_t * const row)
{
    return row[OFFSET];
}

#include "../../../src/neogb/la_gf2.c"

static uint32_t rs = 12345;

static uint32_t rnd(void)
{
    rs = rs * 1103515245u + 12345u;
    return rs >> 8;
}

static int fails = 0;

static void check(const int ok, const char *what)
{
    if (!ok) {
        fails++;
        fprintf(stderr, "failed: %s\n", what);
    }
}

static int bit(const uint64_t *d, const len_t i)
{
    return (d[i >> 6] >> (i & 63)) & 1;
}

static void check_kernels(void)
{
    uint64_t a[64] __attribute__((aligned(64)));
    uint64_t b[64] __attribute__((aligned(64)));
    uint64_t c[64];

    for (int it = 0; it < 1000; ++it) {
        for (len_t i = 0; i < 64; ++i) {
            a[i] = ((uint64_t)rnd() << 40) ^ ((uint64_t)rnd() << 20) ^ rnd();
            b[i] = ((uint64_t)rnd() << 40) ^ ((uint64_t)rnd() << 20) ^ rnd();
            c[i] = a[i] ^ b[i];
        }
        /* lanes */
        gf2_lanes_xor(a, b);
        check(memcmp(a, c, GF2_LW * sizeof(uint64_t)) == 0, "lanes xor");
        memset(a, 0, sizeof(a));
        check(gf2_lanes_zero(a), "zero lanes");
        a[rnd() % GF2_LW] = (uint64_t)1 << (rnd() % 64);
        check(!gf2_lanes_zero(a), "nonzero lanes");

        /* regions of any length and offset */
        for (len_t i = 0; i < 64; ++i) {
            a[i] = ((uint64_t)rnd() << 40) ^ ((uint64_t)rnd() << 20) ^ rnd();
            c[i] = a[i];
        }
        const len_t o = rnd() % 8;
        const len_t n = rnd() % (64 - o);
        gf2_xor_region(a + o, b + o, n);
        for (len_t i = o; i < o + n; ++i) {
            c[i] ^= b[i];
        }
        check(memcmp(a, c, sizeof(a)) == 0, "region xor");

        /* next set bit */
        memset(a, 0, sizeof(a));
        const len_t nb = 1 + rnd() % (64 * 64);
        const len_t k = rnd() % 4;
        for (len_t i = 0; i < k; ++i) {
            const len_t p = rnd() % nb;
            a[p >> 6] |= (uint64_t)1 << (p & 63);
        }
        const len_t s = rnd() % (nb + 1);
        len_t e = s;
        while (e < nb && !bit(a, e)) {
            e++;
        }
        check(gf2_next_bit(a, s, nb) == e, "next bit");
    }
}

/* reduces random rows as one group and one by one by random pivots */
static void check_group_reduction(void)
{
    const len_t ncl = 300, ncols = 400;
    mat_t mat;
    memset(&mat, 0, sizeof(mat));
    mat.ncl = ncl;
    mat.nc  = ncols;

    hm_t **pivs = (hm_t **)calloc(ncols, sizeof(hm_t *));
    for (len_t i = 0; i < ncl; ++i) {
        len_t len = 1;
        hm_t cols[64];
        cols[0] = i;
        for (len_t c = i+1; c < ncols && len < 64; ++c) {
            if (rnd() % 16 == 0) {
                cols[len++] = c;
            }
        }
        pivs[i] = (hm_t *)malloc((len + OFFSET) * sizeof(hm_t));
        pivs[i][LENGTH] = len;
        memcpy(pivs[i] + OFFSET, cols, len * sizeof(hm_t));
    }

    void *raw = NULL;
    uint64_t *dm = gf2_alloc_aligned(&raw, (uint64_t)ncols * GF2_LW);
    const len_t nw = (ncols + 63) / 64;
    uint64_t *rows = (uint64_t *)calloc((uint64_t)GF2_NR * nw, sizeof(uint64_t));
    for (len_t r = 0; r < GF2_NR; ++r) {
        for (len_t c = 0; c < ncols; ++c) {
            if (rnd() % 8 == 0) {
                rows[r * nw + (c >> 6)] |= (uint64_t)1 << (c & 63);
                dm[GF2_LW * c + (r >> 6)] |= (uint64_t)1 << (r & 63);
            }
        }
        /* reference: one row at a time */
        uint64_t *d = rows + r * nw;
        for (len_t c = gf2_next_bit(d, 0, ncl); c < ncl;
                c = gf2_next_bit(d, c+1, ncl)) {
            gf2_flip_bits(d, pivs[c] + OFFSET, pivs[c][LENGTH]);
        }
    }
    reduce_dense_rows_by_known_pivots_gf2(dm, &mat, pivs, 0);
    hm_t *rcs[GF2_NR];
    len_t rcl[GF2_NR];
    extract_rows_from_group_gf2(dm, ncl, ncols, GF2_NR, rcs, rcl);
    int ok = 1;
    for (len_t r = 0; r < GF2_NR; ++r) {
        uint64_t d[8] = {0};
        gf2_load_row(d, rcs[r], rcl[r]);
        ok &= memcmp(d, rows + r * nw, nw * sizeof(uint64_t)) == 0;
        free(rcs[r]);
    }
    check(ok, "group reduction");
    for (len_t c = 0; c < ncols * GF2_LW; ++c) {
        ok &= dm[c] == 0;
    }
    check(ok, "group is zero after extraction");

    for (len_t i = 0; i < ncl; ++i) {
        free(pivs[i]);
    }
    free(pivs);
    free(rows);
    free(raw);
}

int main(void)
{
    check_kernels();
    check_group_reduction();

    return fails != 0;
}
