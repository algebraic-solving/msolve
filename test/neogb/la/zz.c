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

/* Checks the linear algebra over the integers in la_zz.c on random
 * matrices: the rows of the result generate the same lattice as the input
 * rows (compared via Hermite normal forms computed independently here),
 * they are in echelon form with the gcds of the columns as lead
 * coefficients, and their entries are reduced into the residue systems of
 * the pivots. */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../../src/neogb/data.h"
#include "../../../src/neogb/tools.h"
#include "../../../src/msolve/streams.h"

/* stub for the function of the library la_zz.c needs which is not
 * exported */
static void free_basis_elements(bs_t *bs)
{
    (void)bs;
}

#include "../../../src/neogb/la_zz.c"

static void free_test_row(hm_t *row, mpz_t *cf)
{
    for (len_t j = 0; j < row[LENGTH]; ++j) {
        mpz_clear(cf[j]);
    }
    free(cf);
    free(row);
}

static uint32_t rs = 12345;

static uint32_t rnd(void)
{
    rs = rs * 1103515245u + 12345u;
    return rs >> 8;
}

static int fails = 0;

static void check(const int ok, const char *what, const int it)
{
    if (!ok) {
        fails++;
        fprintf(stderr, "failed: %s (instance %d)\n", what, it);
    }
}

/* canonical Hermite normal form of the n x nc matrix M (row major): lead
 * coefficients positive, entries above a pivot in [0, pivot). The rows are
 * inserted one by one into an echelon basis H, by Euclid's algorithm on
 * the pivot row and the inserted row; after each insertion all entries of
 * H above its pivots are reduced, which keeps the entries small. Returns
 * the rank, the first rank rows of M hold the form. */
static len_t hnf(mpz_t *M, const len_t n, const len_t nc)
{
    len_t i, j, c, r;
    mpz_t q;
    mpz_init(q);
    mpz_t *H = (mpz_t *)malloc((unsigned long)nc * nc * sizeof(mpz_t));
    for (i = 0; i < nc * nc; ++i) {
        mpz_init(H[i]);
    }
    int *has = (int *)calloc(nc, sizeof(int));
    for (i = 0; i < n; ++i) {
        mpz_t *x = M + i*nc;
        for (c = 0; c < nc; ++c) {
            if (mpz_sgn(x[c]) == 0) {
                continue;
            }
            mpz_t *p = H + c*nc;
            if (has[c] == 0) {
                for (j = c; j < nc; ++j) {
                    mpz_swap(p[j], x[j]);
                }
                has[c] = 1;
                break;
            }
            while (mpz_sgn(x[c]) != 0) {
                mpz_fdiv_q(q, p[c], x[c]);
                for (j = c; j < nc; ++j) {
                    mpz_submul(p[j], q, x[j]);
                    mpz_swap(p[j], x[j]);
                }
            }
        }
        /* positive pivots, reduced entries above them */
        for (c = 0; c < nc; ++c) {
            if (has[c] == 0) {
                continue;
            }
            mpz_t *p = H + c*nc;
            if (mpz_sgn(p[c]) < 0) {
                for (j = c; j < nc; ++j) {
                    mpz_neg(p[j], p[j]);
                }
            }
            for (r = 0; r < c; ++r) {
                if (has[r] == 0 || mpz_sgn(H[r*nc+c]) == 0) {
                    continue;
                }
                mpz_fdiv_q(q, H[r*nc+c], p[c]);
                for (j = c; j < nc; ++j) {
                    mpz_submul(H[r*nc+j], q, p[j]);
                }
            }
        }
    }
    r = 0;
    for (c = 0; c < nc; ++c) {
        if (has[c] == 0) {
            continue;
        }
        for (j = 0; j < nc; ++j) {
            mpz_swap(M[r*nc+j], H[c*nc+j]);
        }
        r++;
    }
    for (i = r * nc; i < n * nc; ++i) {
        mpz_set_ui(M[i], 0);
    }
    for (i = 0; i < nc * nc; ++i) {
        mpz_clear(H[i]);
    }
    free(H);
    free(has);
    mpz_clear(q);
    return r;
}

static mpz_t *dense_alloc(const len_t n, const len_t nc)
{
    mpz_t *M = (mpz_t *)malloc((unsigned long)n * nc * sizeof(mpz_t));
    for (len_t i = 0; i < n * nc; ++i) {
        mpz_init(M[i]);
    }
    return M;
}

static void dense_free(mpz_t *M, const len_t n, const len_t nc)
{
    for (len_t i = 0; i < n * nc; ++i) {
        mpz_clear(M[i]);
    }
    free(M);
}

static void dense_set_row(mpz_t *d, const hm_t *row, mpz_t *cf)
{
    for (len_t j = 0; j < row[LENGTH]; ++j) {
        mpz_set(d[row[OFFSET+j]], cf[j]);
    }
}

/* number of bits the random coefficients are shifted by, coefficients
 * beyond 32 bits are reduced with GMP */
static unsigned zz_shift = 0;

/* v * 2^zz_shift plus some noise of the same sign */
static void set_coefficient(mpz_t c, const long v)
{
    mpz_init_set_si(c, v);
    if (zz_shift > 0) {
        mpz_mul_2exp(c, c, zz_shift);
        if (v > 0) {
            mpz_add_ui(c, c, rnd() % 1000);
        } else {
            mpz_sub_ui(c, c, rnd() % 1000);
        }
    }
}

/* random sparse row with lead column lc, lead coefficient in [1, ml]
 * (positive) or in [-ml, ml] \ {0}, other entries in [-mt, mt], all of
 * them shifted by zz_shift bits */
static hm_t *random_row(mpz_t **cfp, const len_t lc, const len_t nc,
        const int pos, const long ml, const long mt, const uint32_t dens)
{
    hm_t *row = (hm_t *)malloc((unsigned long)(nc-lc+OFFSET) * sizeof(hm_t));
    mpz_t *cf = (mpz_t *)malloc((unsigned long)(nc-lc) * sizeof(mpz_t));
    len_t len = 0;
    long v    = 1 + (long)(rnd() % ml);
    if (pos == 0 && rnd() % 2 == 0) {
        v = -v;
    }
    set_coefficient(cf[len], v);
    row[OFFSET+len++] = lc;
    for (len_t c = lc+1; c < nc; ++c) {
        if (rnd() % dens != 0) {
            continue;
        }
        v = (long)(rnd() % (2*mt+1)) - mt;
        if (v != 0) {
            set_coefficient(cf[len], v);
            row[OFFSET+len++] = c;
        }
    }
    row[PRELOOP] = len % UNROLL;
    row[LENGTH]  = len;
    *cfp = cf;
    return row;
}

/* hash of the rows, for comparing the results of different numbers of
 * threads */
static uint64_t hash_rows(hm_t **rows, mpz_t **cfs, const len_t n)
{
    uint64_t h = 1469598103934665603ULL;
    for (len_t i = 0; i < n; ++i) {
        for (len_t j = 0; j < rows[i][LENGTH]; ++j) {
            h = (h ^ rows[i][OFFSET+j]) * 1099511628211ULL;
            char *str = mpz_get_str(NULL, 16, cfs[i][j]);
            for (const char *c = str; *c != '\0'; ++c) {
                h = (h ^ (uint64_t)*c) * 1099511628211ULL;
            }
            free(str);
        }
        h = (h ^ 0xff) * 1099511628211ULL;
    }
    return h;
}

/* checks that the rows (rows[i], cfs[i]), i < n, have distinct lead
 * columns, positive lead coefficients equal to the pivots of the Hermite
 * normal form H of rank rk, and that the rows in red are reduced: their
 * entries in the lead columns of the other rows lie in (-a/2, a/2] */
static void check_echelon(hm_t **rows, mpz_t **cfs, const int8_t *red,
        const len_t n, mpz_t *H, const len_t rk, const len_t nc,
        const int it)
{
    len_t i, j, k;
    mpz_t **piv = (mpz_t **)calloc(nc, sizeof(mpz_t *));
    int ok = 1;
    for (i = 0; i < n; ++i) {
        const hm_t c = rows[i][OFFSET];
        ok &= piv[c] == NULL && mpz_sgn(cfs[i][0]) > 0;
        piv[c] = cfs[i];
    }
    check(ok, "distinct lead columns, positive lead coefficients", it);
    check(n == rk, "rank", it);
    for (i = 0, k = 0; i < nc && k < rk; ++i) {
        if (piv[i] != NULL) {
            ok &= mpz_cmp(piv[i][0], H[k*nc+i]) == 0;
            k++;
        }
    }
    check(ok, "lead coefficients are the pivots of the HNF", it);
    mpz_t t;
    mpz_init(t);
    for (i = 0; i < n; ++i) {
        if (red[i] == 0) {
            continue;
        }
        for (j = 1; j < rows[i][LENGTH]; ++j) {
            mpz_t *p = piv[rows[i][OFFSET+j]];
            if (p == NULL) {
                continue;
            }
            /* -a < 2e <= a */
            mpz_mul_2exp(t, cfs[i][j], 1);
            ok &= mpz_cmp(t, p[0]) <= 0;
            mpz_neg(t, t);
            ok &= mpz_cmp(t, p[0]) < 0;
        }
    }
    check(ok, "entries reduced modulo the pivots", it);
    mpz_clear(t);
    free(piv);
}

/* reduces nrl random rows by ncl reducer rows with random positive lead
 * coefficients in the first ncl columns */
static uint64_t check_echelon_form(const int it, const len_t nc,
        const len_t ncl, const len_t nrl, const long ml, const long mt,
        const uint32_t dens, const int nthrds, int64_t *ngcd)
{
    len_t i;
    const len_t nr = ncl + nrl;

    bs_t bs;
    mat_t mat;
    md_t md;
    memset(&bs, 0, sizeof(bs));
    memset(&mat, 0, sizeof(mat));
    memset(&md, 0, sizeof(md));
    md.nthrds = nthrds;
    bs.cf_qq = (mpz_t **)malloc(nr * sizeof(mpz_t *));
    mat.rr   = (hm_t **)malloc(ncl * sizeof(hm_t *));
    mat.tr   = (hm_t **)malloc(nrl * sizeof(hm_t *));
    mat.nc   = nc;
    mat.ncl  = ncl;
    mat.ncr  = nc - ncl;
    mat.nru  = ncl;
    mat.nrl  = nrl;

    mpz_t *in   = dense_alloc(nr, nc);
    mpz_t *known = dense_alloc(ncl, nc);
    for (i = 0; i < ncl; ++i) {
        mat.rr[i] = random_row(&bs.cf_qq[i], i, nc, 1, ml, mt, dens);
        mat.rr[i][COEFFS] = i;
        dense_set_row(in + i*nc, mat.rr[i], bs.cf_qq[i]);
        dense_set_row(known + i*nc, mat.rr[i], bs.cf_qq[i]);
    }
    for (i = 0; i < nrl; ++i) {
        mat.tr[i] = random_row(&bs.cf_qq[ncl+i], rnd() % nc, nc, 0, ml, mt,
                dens);
        mat.tr[i][COEFFS] = ncl+i;
        dense_set_row(in + (ncl+i)*nc, mat.tr[i], bs.cf_qq[ncl+i]);
    }
    /* some rows are multiples or combinations of others */
    for (i = 0; i + 1 < nrl; i += 5) {
        hm_t *r   = mat.tr[i];
        mpz_t *cf = bs.cf_qq[r[COEFFS]];
        for (len_t j = 0; j < r[LENGTH]; ++j) {
            mpz_mul_si(cf[j], cf[j], 2 + (long)(rnd() % 3));
        }
        mpz_t *d = in + (ncl+i)*nc;
        for (len_t j = 0; j < nc; ++j) {
            mpz_set_ui(d[j], 0);
        }
        dense_set_row(d, r, cf);
    }

    /* the input rows are freed, keep the lengths of their coefficients */
    len_t *blen = (len_t *)malloc(nr * sizeof(len_t));
    for (i = 0; i < ncl; ++i) {
        blen[i] = mat.rr[i][LENGTH];
    }
    for (i = 0; i < nrl; ++i) {
        blen[ncl+i] = mat.tr[i][LENGTH];
    }

    const len_t nzero = exact_sparse_reduced_echelon_form_zz(&mat, &bs, &md);
    *ngcd += md.num_gcdred;

    /* the result together with the reducer rows which were not replaced */
    const len_t np = mat.np;
    int8_t *isnew  = (int8_t *)calloc(nc, sizeof(int8_t));
    int ok = 1;
    for (i = 0; i < np; ++i) {
        isnew[mat.tr[i][OFFSET]] = 1;
        ok &= i == 0 || mat.tr[i][OFFSET] < mat.tr[i-1][OFFSET];
        ok &= mat.tr[i][COEFFS] == mat.tr[i][OFFSET];
    }
    check(ok, "new rows ordered by decreasing lead column", it);

    hm_t **rows  = (hm_t **)malloc((np + ncl) * sizeof(hm_t *));
    mpz_t **cfs  = (mpz_t **)malloc((np + ncl) * sizeof(mpz_t *));
    int8_t *red  = (int8_t *)calloc(np + ncl, sizeof(int8_t));
    len_t n = 0;
    for (i = 0; i < np; ++i) {
        rows[n]  = mat.tr[i];
        cfs[n]   = mat.cf_qq[mat.tr[i][COEFFS]];
        red[n++] = 1;
    }
    /* reducer rows were freed, rebuild the remaining ones from the copies */
    for (i = 0; i < ncl; ++i) {
        if (isnew[i] != 0) {
            continue;
        }
        rows[n] = (hm_t *)malloc((nc+OFFSET) * sizeof(hm_t));
        cfs[n]  = (mpz_t *)malloc(nc * sizeof(mpz_t));
        len_t len = 0;
        for (len_t j = i; j < nc; ++j) {
            if (mpz_sgn(known[i*nc+j]) != 0) {
                mpz_init_set(cfs[n][len], known[i*nc+j]);
                rows[n][OFFSET+len++] = j;
            }
        }
        rows[n][LENGTH] = len;
        n++;
    }
    mpz_t *out = dense_alloc(n, nc);
    for (i = 0; i < n; ++i) {
        dense_set_row(out + i*nc, rows[i], cfs[i]);
    }
    const len_t rin  = hnf(in, nr, nc);
    const len_t rout = hnf(out, n, nc);
    ok = rin == rout;
    for (i = 0; ok && i < rin * nc; ++i) {
        ok &= mpz_cmp(in[i], out[i]) == 0;
    }
    check(ok, "same lattice", it);
    check(nzero + np >= nrl, "number of zero rows", it);
    check_echelon(rows, cfs, red, n, in, rin, nc, it);
    const uint64_t h = hash_rows(rows, cfs, np);

    for (i = 0; i < n; ++i) {
        free_test_row(rows[i], cfs[i]);
    }
    free(rows);
    free(cfs);
    free(red);
    free(isnew);
    dense_free(in, nr, nc);
    dense_free(out, n, nc);
    dense_free(known, ncl, nc);
    for (i = 0; i < nr; ++i) {
        for (len_t j = 0; j < blen[i]; ++j) {
            mpz_clear(bs.cf_qq[i][j]);
        }
        free(bs.cf_qq[i]);
    }
    free(blen);
    free(mat.cf_qq);
    free(mat.tr);
    free(mat.rr);
    free(bs.cf_qq);
    return h;
}

/* interreduces n random rows, several of them with equal lead columns */
static uint64_t check_interreduction(const int it, const len_t nc,
        const len_t n, const long ml, const long mt, const uint32_t dens,
        const int nthrds, int64_t *ngcd)
{
    len_t i;

    bs_t bs;
    mat_t mat;
    md_t md;
    memset(&bs, 0, sizeof(bs));
    memset(&mat, 0, sizeof(mat));
    memset(&md, 0, sizeof(md));
    md.nthrds = nthrds;
    bs.cf_qq = (mpz_t **)malloc(n * sizeof(mpz_t *));
    mat.rr   = (hm_t **)malloc(n * sizeof(hm_t *));
    mat.nr   = n;
    mat.nc   = nc;

    mpz_t *in   = dense_alloc(n, nc);
    len_t *blen = (len_t *)malloc(n * sizeof(len_t));
    for (i = 0; i < n; ++i) {
        mat.rr[i] = random_row(&bs.cf_qq[i], rnd() % (nc/2 + 1), nc, 0, ml,
                mt, dens);
        mat.rr[i][COEFFS] = i;
        blen[i] = mat.rr[i][LENGTH];
        dense_set_row(in + i*nc, mat.rr[i], bs.cf_qq[i]);
    }

    interreduce_matrix_rows_zz(&mat, &bs, &md, 0);
    *ngcd += md.num_gcdred;

    const len_t np = mat.np;
    int ok = mat.rr == NULL;
    for (i = 0; i < np; ++i) {
        ok &= i == 0 || mat.tr[i][OFFSET] > mat.tr[i-1][OFFSET];
    }
    check(ok, "interreduced rows ordered by increasing lead column", it);
    mpz_t **cfs = (mpz_t **)malloc((np + 1) * sizeof(mpz_t *));
    int8_t *red = (int8_t *)malloc(np + 1);
    mpz_t *out  = dense_alloc(np, nc);
    for (i = 0; i < np; ++i) {
        cfs[i] = mat.cf_qq[mat.tr[i][COEFFS]];
        red[i] = 1;
        dense_set_row(out + i*nc, mat.tr[i], cfs[i]);
    }
    const len_t rin  = hnf(in, n, nc);
    const len_t rout = hnf(out, np, nc);
    ok = rin == rout;
    for (i = 0; ok && i < rin * nc; ++i) {
        ok &= mpz_cmp(in[i], out[i]) == 0;
    }
    check(ok, "interreduction keeps the lattice", it);
    check_echelon(mat.tr, cfs, red, np, in, rin, nc, it);
    const uint64_t h = hash_rows(mat.tr, cfs, np);

    for (i = 0; i < np; ++i) {
        free_test_row(mat.tr[i], cfs[i]);
    }
    for (i = 0; i < n; ++i) {
        for (len_t j = 0; j < blen[i]; ++j) {
            mpz_clear(bs.cf_qq[i][j]);
        }
        free(bs.cf_qq[i]);
    }
    free(blen);
    free(cfs);
    free(red);
    dense_free(in, n, nc);
    dense_free(out, np, nc);
    free(mat.cf_qq);
    free(mat.tr);
    free(bs.cf_qq);
    return h;
}

int main(void)
{
    int64_t ngcd = 0;
    int it = 0;

    /* small and large instances, with coefficients of up to 32 bits and
     * beyond, each with 1 and 8 threads; the results must agree. Random
     * lattices of lower rank have huge entries in their non-pivot columns,
     * so the large instances have more rows than columns. */
    for (int k = 0; k < 300; ++k, ++it) {
        const int large = k % 4 == 3;
        const len_t nc  = 1 + rnd() % (large ? 120 : 40);
        const len_t ncl = rnd() % (nc + 1);
        const len_t nrl = large ? nc + rnd() % 150 : rnd() % 30;
        const long ml   = 1 + rnd() % (k % 3 == 0 ? 60 : 1000000);
        const long mt   = 1 + rnd() % (k % 3 == 0 ? 30 : 1000000);
        const uint32_t dens = 1 + rnd() % 6;
        zz_shift = k % 5 == 4 ? 20 + rnd() % 60 : 0;
        const uint32_t seed = rs;
        const uint64_t h1 = check_echelon_form(it, nc, ncl, nrl, ml, mt, dens,
                1, &ngcd);
        rs = seed;
        const uint64_t h8 = check_echelon_form(it, nc, ncl, nrl, ml, mt, dens,
                8, &ngcd);
        check(h1 == h8, "same result for 1 and 8 threads", it);
    }
    for (int k = 0; k < 300; ++k, ++it) {
        const int large = k % 4 == 3;
        const len_t nc  = 1 + rnd() % (large ? 120 : 40);
        const len_t n   = large ? 2 * nc + rnd() % 150 : 1 + rnd() % 30;
        const long ml   = 1 + rnd() % (k % 3 == 0 ? 60 : 1000000);
        const long mt   = 1 + rnd() % (k % 3 == 0 ? 30 : 1000000);
        const uint32_t dens = 1 + rnd() % 6;
        zz_shift = k % 5 == 4 ? 20 + rnd() % 60 : 0;
        const uint32_t seed = rs;
        const uint64_t h1 = check_interreduction(it, nc, n, ml, mt, dens, 1,
                &ngcd);
        rs = seed;
        const uint64_t h8 = check_interreduction(it, nc, n, ml, mt, dens, 8,
                &ngcd);
        check(h1 == h8, "same interreduction for 1 and 8 threads", it);
    }
    check(ngcd > 0, "gcd steps occur", it);

    return fails != 0;
}
