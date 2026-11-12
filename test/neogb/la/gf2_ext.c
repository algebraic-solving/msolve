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

/* Checks the vectorized arithmetic in GF(16) and GF(256) of gf2_ext.c,
 * compiled with the same SIMD flags as the library, against the scalar
 * multiplication tables, and the validation of input data over these
 * fields. */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../src/neogb/data.h"
#include "../../../src/neogb/io.h"
#include "../../../src/neogb/gf2_ext.c"

static uint32_t rs = 12345;

static uint8_t rnd8(void)
{
    rs = rs * 1103515245u + 12345u;
    return (uint8_t)(rs >> 16);
}

static int fails = 0;

static void check(const int ok, const char *what)
{
    if (!ok) {
        fails++;
        fprintf(stderr, "failed: %s\n", what);
    }
}

/* multiplication in GF(2)[x]/(poly) without tables */
static uint8_t mul_slow(uint8_t a, uint8_t b, const uint32_t poly,
        const uint32_t deg)
{
    uint32_t p = 0, x = a;
    while (b != 0) {
        if (b & 1) {
            p ^= x;
        }
        b >>= 1;
        x <<= 1;
        if (x & (1u << deg)) {
            x ^= poly;
        }
    }
    return (uint8_t)p;
}

static void check_field(const uint64_t fc, const uint32_t poly,
        const uint32_t deg)
{
    const gf2x_t *t = gf2x_table(fc);
    const uint32_t q = (uint32_t)fc;
    uint8_t a[1000], b[1000], c[1000];
    uint8_t m[GF_NR] __attribute__((aligned(64)));
    uint8_t ml[GF_NR] __attribute__((aligned(64)));
    uint8_t mh[GF_NR] __attribute__((aligned(64)));
    uint8_t d[4*GF_NR] __attribute__((aligned(64)));
    uint8_t r[4*GF_NR];

    for (uint32_t x = 0; x < q; ++x) {
        if (x > 0) {
            check(t->mul[x][t->inv[x]] == 1, "inverse");
        }
        for (uint32_t y = 0; y < q; ++y) {
            check(t->mul[x][y] == mul_slow(x, y, poly, deg), "table");
            check((t->lo[x][y & 15] ^ t->hi[x][y >> 4]) == t->mul[x][y],
                    "nibble tables");
        }
    }
    for (int it = 0; it < 500; ++it) {
        const len_t n   = rnd8() * 3 + rnd8() % 7;
        const uint8_t e = (uint8_t)(rnd8() % q);
        for (len_t i = 0; i < n; ++i) {
            a[i] = (uint8_t)(rnd8() % q);
            b[i] = c[i] = (uint8_t)(rnd8() % q);
        }
        gf2x_mul_region(b, a, n, e, t, 1);
        for (len_t i = 0; i < n; ++i) {
            check(b[i] == (c[i] ^ t->mul[e][a[i]]), "axpy");
        }
        memcpy(b, a, n);
        gf2x_mul_region(b, b, n, e, t, 0);
        for (len_t i = 0; i < n; ++i) {
            check(b[i] == t->mul[e][a[i]], "scaling");
        }
        memset(a, 0, sizeof(a));
        const len_t p = rnd8() % (n + 1);
        if (p < n) {
            a[p] = (uint8_t)(1 + rnd8() % (q-1));
        }
        const len_t s = rnd8() % (p + 1);
        check(gf2x_next_nonzero(a, s, n) == p, "next nonzero entry");

        for (len_t i = 0; i < GF_NR; ++i) {
            m[i] = rnd8() % 3 ? (uint8_t)(rnd8() % q) : 0;
        }
        for (len_t i = 0; i < 4*GF_NR; ++i) {
            d[i] = r[i] = (uint8_t)(rnd8() % q);
        }
        gf2x_lanes_split(ml, mh, m);
        gf2x_lanes_axpy(d + GF_NR, ml, mh, e, t, q == 16);
        for (len_t i = 0; i < GF_NR; ++i) {
            r[GF_NR + i] ^= t->mul[e][m[i]];
        }
        check(memcmp(d, r, sizeof(r)) == 0, "lanes axpy");
        memset(d, 0, sizeof(d));
        check(gf2x_lanes_zero(d + GF_NR), "zero lanes");
        d[GF_NR + rnd8() % GF_NR] = 1;
        check(!gf2x_lanes_zero(d + GF_NR), "nonzero lanes");
    }
}

/* zero coefficients over GF(256) are removed, only generators without
 * nonzero coefficient are invalid; the number of generators is not
 * changed */
static void check_validation(void)
{
    const int32_t lens[]    =   {2, 2, 1};
    const int32_t cfs[]     =   {256, -512, 3, 0x100, -7};
    uint32_t field_char     =   256;
    int32_t mon_order       =   0;
    int32_t elim_block_len  =   0;
    int32_t nr_vars         =   2;
    int32_t nr_gens         =   3;
    int32_t nr_nf           =   0;
    int32_t ht_size         =   17;
    int32_t nr_threads      =   1;
    int32_t max_nr_pairs    =   0;
    int32_t reset_ht        =   0;
    int32_t la_option       =   2;
    int32_t use_signatures  =   1; /* not supported, disabled */
    int32_t reduce_gb       =   1;
    int32_t truncate_lifting =  0;
    int32_t info_level      =   0;
    int *invalid_gens       =   NULL;

    int res = validate_input_data(&invalid_gens, cfs, lens,
            &field_char, &mon_order, &elim_block_len, &nr_vars,
            &nr_gens, &nr_nf, &ht_size, &nr_threads, &max_nr_pairs,
            &reset_ht, &la_option, &use_signatures, &reduce_gb,
            &truncate_lifting, &info_level);

    check(res == 1, "validation result");
    check(invalid_gens[0] == 1 && invalid_gens[1] == 0
            && invalid_gens[2] == 0, "invalid generators");
    check(nr_gens == 3, "number of generators");
    check(use_signatures == 0, "signatures disabled");
    free(invalid_gens);

    /* normal forms are not supported */
    nr_nf = 1;
    res = validate_input_data(&invalid_gens, cfs, lens,
            &field_char, &mon_order, &elim_block_len, &nr_vars,
            &nr_gens, &nr_nf, &ht_size, &nr_threads, &max_nr_pairs,
            &reset_ht, &la_option, &use_signatures, &reduce_gb,
            &truncate_lifting, &info_level);
    check(res == 0, "normal forms rejected");
}

int main(void)
{
    gf2x_init_tables();
    check_field(16, 0x13, 4);
    check_field(256, 0x11B, 8);
    check_validation();

    return fails != 0;
}
