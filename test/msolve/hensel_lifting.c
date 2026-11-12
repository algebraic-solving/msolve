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

/* tests of the p-adic lifting of hensel.c: the products of polynomial
 * matrices by transforms, the inversion of J modulo (w, p), and the lifting
 * of parametrizations of systems with known solutions, recombined by CRT,
 * reconstructed over Q and compared with the parametrization modulo another
 * prime */

#include "../../src/msolve/libmsolve.c"

#include <flint/fmpz_mpoly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

static int nfail = 0;

/* c modulo p, in [0, p) */
static mp_limb_t mod_si(const long c, const mp_limb_t p) {
  return c >= 0 ? (mp_limb_t)c % p : (p - (mp_limb_t)(-c) % p) % p;
}

#define CHECK(cond, ...)                                                    \
  do {                                                                      \
    if (!(cond)) {                                                          \
      fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__);                  \
      fprintf(stderr, __VA_ARGS__);                                         \
      fprintf(stderr, "\n");                                                \
      nfail++;                                                              \
    }                                                                       \
  } while (0)

/* a fixed generator, independent of the version of FLINT */
static uint64_t seed = 0x853c49e6748fea9bULL;

static uint64_t rnd(void) {
  seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
  return seed >> 11;
}

/* ------------------------------------------------------------------ */
/* products by transforms                                              */

#if HENSEL_NTT
/* a random polynomial of length at most D with nonnegative coefficients of
 * at most bits bits, zero with probability 1/8 */
static void rnd_poly(fmpz_poly_t a, const long D, const long bits) {
  fmpz_poly_zero(a);
  if (rnd() % 8 == 0) {
    return;
  }
  fmpz_t c;
  fmpz_init(c);
  for (long k = 0; k < D; k++) {
    fmpz_zero(c);
    for (long b = 0; b < bits; b += 32) {
      fmpz_mul_2exp(c, c, 32);
      fmpz_add_ui(c, c, rnd() & 0xffffffffULL);
    }
    fmpz_fdiv_r_2exp(c, c, bits);
    fmpz_poly_set_coeff_fmpz(a, k, c);
  }
  fmpz_clear(c);
}

/* C = A B by hensel_ntt_mul against sums of fmpz_poly_mul */
static void test_ntt_mul(const long r, const long c, const long n,
                         const long D, const long bA, const long bB,
                         hensel_ntt_t *T) {
  fmpz_poly_struct *A = malloc(r * c * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *B = malloc(c * n * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *C = malloc(r * n * sizeof(fmpz_poly_struct));
  for (long i = 0; i < r * c; i++) {
    fmpz_poly_init(A + i);
    rnd_poly(A + i, D, bA);
  }
  for (long i = 0; i < c * n; i++) {
    fmpz_poly_init(B + i);
    rnd_poly(B + i, D, bB);
  }
  for (long i = 0; i < r * n; i++) {
    fmpz_poly_init(C + i);
  }
  hensel_ntt_mat_t M;
  const int ok = hensel_ntt_mat_set(&M, A, r, c, bB, D, T);
  CHECK(ok, "ntt: no transforms for %ld x %ld bits", bA, bB);
  if (ok) {
    hensel_ntt_mul(C, &M, B, n, D, T);
    hensel_ntt_mat_clear(&M);
    fmpz_poly_t s, pr;
    fmpz_poly_init(s);
    fmpz_poly_init(pr);
    for (long i = 0; i < r; i++) {
      for (long l = 0; l < n; l++) {
        fmpz_poly_zero(s);
        for (long j = 0; j < c; j++) {
          fmpz_poly_mul(pr, A + i * c + j, B + j * n + l);
          fmpz_poly_add(s, s, pr);
        }
        CHECK(fmpz_poly_equal(s, C + i * n + l),
              "ntt: product %ld x %ld x %ld, D = %ld, %ld x %ld bits, "
              "entry (%ld, %ld)", r, c, n, D, bA, bB, i, l);
      }
    }
    fmpz_poly_clear(s);
    fmpz_poly_clear(pr);
  }
  for (long i = 0; i < r * c; i++) {
    fmpz_poly_clear(A + i);
  }
  for (long i = 0; i < c * n; i++) {
    fmpz_poly_clear(B + i);
  }
  for (long i = 0; i < r * n; i++) {
    fmpz_poly_clear(C + i);
  }
  free(A);
  free(B);
  free(C);
}

static void test_ntt(void) {
  hensel_ntt_t T;
  hensel_ntt_init(&T);
  /* the transforms are kept by T from one size to the next */
  test_ntt_mul(3, 3, 1, 100, 60, 60, &T);
  test_ntt_mul(4, 4, 4, 200, 500, 30, &T);
  test_ntt_mul(2, 9, 3, 64, 1, 1000, &T);
  test_ntt_mul(5, 5, 1, 300, 1200, 1200, &T);
  test_ntt_mul(3, 3, 2, 100, 60, 60, &T);
  /* coefficients beyond the HENSEL_NTT_MAXP primes */
  fmpz_poly_t a;
  fmpz_poly_init(a);
  fmpz_t c;
  fmpz_init(c);
  fmpz_one(c);
  fmpz_mul_2exp(c, c, 3000);
  fmpz_poly_set_coeff_fmpz(a, 1, c);
  hensel_ntt_mat_t M;
  CHECK(!hensel_ntt_mat_set(&M, a, 1, 1, 3000, 100, &T),
        "ntt: transforms for 3000 x 3000 bits");
  hensel_ntt_mat_clear(&M);
  fmpz_clear(c);
  fmpz_poly_clear(a);
  hensel_ntt_clear(&T);
}
#endif

/* ------------------------------------------------------------------ */
/* inversion modulo (f, p)                                             */

/* hensel_inv_mod_p on A (kept), checked against the determinant: it
 * succeeds if and only if det A is invertible modulo f, and then A B = I
 * modulo f; returns its result */
static int check_inv(const nmod_poly_struct *A, const nmod_poly_t f,
                     const long nv, const char *what) {
  const mp_limb_t p = f->mod.n;
  nmod_poly_struct *Ac = malloc(nv * nv * sizeof(nmod_poly_struct));
  nmod_poly_struct *B = malloc(nv * nv * sizeof(nmod_poly_struct));
  for (long i = 0; i < nv * nv; i++) {
    nmod_poly_init(Ac + i, p);
    nmod_poly_init(B + i, p);
    nmod_poly_set(Ac + i, A + i);
  }
  const int ok = hensel_inv_mod_p(B, Ac, f, nv);
  nmod_poly_mat_t M;
  nmod_poly_mat_init(M, nv, nv, p);
  for (long i = 0; i < nv; i++) {
    for (long j = 0; j < nv; j++) {
      nmod_poly_set(nmod_poly_mat_entry(M, i, j), A + i * nv + j);
    }
  }
  nmod_poly_t d, g, s, t;
  nmod_poly_init(d, p);
  nmod_poly_init(g, p);
  nmod_poly_init(s, p);
  nmod_poly_init(t, p);
  nmod_poly_mat_det(d, M);
  nmod_poly_rem(d, d, f);
  nmod_poly_gcd(g, d, f);
  const int inv = !nmod_poly_is_zero(d) && nmod_poly_degree(g) == 0;
  CHECK(ok == inv, "inv %s: returns %d, det %s", what, ok,
        inv ? "invertible" : "not invertible");
  if (ok) {
    for (long i = 0; i < nv; i++) {
      for (long j = 0; j < nv; j++) {
        nmod_poly_zero(s);
        for (long k = 0; k < nv; k++) {
          nmod_poly_mulmod(t, A + i * nv + k, B + k * nv + j, f);
          nmod_poly_add(s, s, t);
        }
        CHECK(nmod_poly_is_zero(s) == (i != j)
              && (i != j || nmod_poly_is_one(s)),
              "inv %s: A B != I at (%ld, %ld)", what, i, j);
      }
    }
  }
  nmod_poly_clear(d);
  nmod_poly_clear(g);
  nmod_poly_clear(s);
  nmod_poly_clear(t);
  nmod_poly_mat_clear(M);
  for (long i = 0; i < nv * nv; i++) {
    nmod_poly_clear(Ac + i);
    nmod_poly_clear(B + i);
  }
  free(Ac);
  free(B);
  return ok;
}

/* the entry from a list of coefficients, from the constant one */
static void set_nmod(nmod_poly_t a, const long n, const long *c) {
  nmod_poly_zero(a);
  for (long k = 0; k < n; k++) {
    nmod_poly_set_coeff_ui(a, k, mod_si(c[k], a->mod.n));
  }
}

static void test_inv(void) {
  const mp_limb_t p = 1073741827;
  nmod_poly_t f, l;
  nmod_poly_init(f, p);
  nmod_poly_init(l, p);
  nmod_poly_struct A[16];
  for (long i = 0; i < 16; i++) {
    nmod_poly_init(A + i, p);
  }
  /* f = t (t - 1): [[t, 1], [t - 1, 1]] has determinant 1, and no entry of
   * its first column is invertible modulo f (jorge square1.sZ) */
  set_nmod(f, 3, (long[]){0, -1, 1});
  set_nmod(A + 0, 2, (long[]){0, 1});
  set_nmod(A + 1, 1, (long[]){1});
  set_nmod(A + 2, 2, (long[]){-1, 1});
  set_nmod(A + 3, 1, (long[]){1});
  CHECK(check_inv(A, f, 2, "split 2 x 2"), "inv: split 2 x 2 not inverted");
  /* the same in the middle of a 3 x 3 matrix, after a first pivot */
  for (long i = 0; i < 9; i++) {
    nmod_poly_zero(A + i);
  }
  set_nmod(A + 0, 1, (long[]){1});
  set_nmod(A + 1, 2, (long[]){3, 5});
  set_nmod(A + 2, 1, (long[]){7});
  set_nmod(A + 4, 2, (long[]){0, 1});
  set_nmod(A + 5, 1, (long[]){1});
  set_nmod(A + 7, 2, (long[]){-1, 1});
  set_nmod(A + 8, 1, (long[]){1});
  CHECK(check_inv(A, f, 3, "split 3 x 3"), "inv: split 3 x 3 not inverted");
  /* [[t, t], [1, 1]] is singular */
  set_nmod(A + 0, 2, (long[]){0, 1});
  set_nmod(A + 1, 2, (long[]){0, 1});
  set_nmod(A + 2, 1, (long[]){1});
  set_nmod(A + 3, 1, (long[]){1});
  CHECK(!check_inv(A, f, 2, "singular"), "inv: singular matrix inverted");
  /* [[t, 0], [0, 1]] is singular modulo the factor t only */
  set_nmod(A + 0, 2, (long[]){0, 1});
  nmod_poly_zero(A + 1);
  nmod_poly_zero(A + 2);
  set_nmod(A + 3, 1, (long[]){1});
  CHECK(!check_inv(A, f, 2, "singular mod t"),
        "inv: matrix singular modulo t inverted");
  /* random 4 x 4 matrices modulo f = (t - 1) .. (t - 8), each entry
   * vanishing modulo some factors, so that the pivots are often made of
   * combinations of rows */
  nmod_poly_one(f);
  for (long r = 1; r <= 8; r++) {
    set_nmod(l, 2, (long[]){-r, 1});
    nmod_poly_mul(f, f, l);
  }
  long ninv = 0;
  for (long it = 0; it < 200; it++) {
    for (long i = 0; i < 16; i++) {
      nmod_poly_zero(A + i);
      for (long k = 0; k < 8; k++) {
        nmod_poly_set_coeff_ui(A + i, k, rnd() % p);
      }
      /* times a product of some factors of f */
      for (long r = 1; r <= 8; r++) {
        if (rnd() % 3 == 0) {
          set_nmod(l, 2, (long[]){-r, 1});
          nmod_poly_mulmod(A + i, A + i, l, f);
        }
      }
    }
    ninv += check_inv(A, f, 4, "random");
  }
  CHECK(ninv > 0, "inv: no random matrix inverted");
  for (long i = 0; i < 16; i++) {
    nmod_poly_clear(A + i);
  }
  nmod_poly_clear(f);
  nmod_poly_clear(l);
}

/* ------------------------------------------------------------------ */
/* lifting                                                              */

/* A triangular system in s and x_0, .., x_{m-1}: G(s) = 0 and c_k x_k =
 * A_k(s), with the separating form t = s + sum_k l_k x_k. The system lifted
 * is in x_0, .., x_{m-1}, t (nv = m + 1): G(t - sum l x) = 0 and c_k x_k -
 * A_k(t - sum l x) = 0, whose parametrization by t has the degree d of G */
typedef struct {
  long m;
  fmpz_poly_t G;
  fmpz_poly_struct *A;
  long *c, *l;
} tri_t;

static void tri_init(tri_t *S, const long m, const long dG, const long *G) {
  S->m = m;
  fmpz_poly_init(S->G);
  for (long k = 0; k <= dG; k++) {
    fmpz_poly_set_coeff_si(S->G, k, G[k]);
  }
  S->A = malloc(m * sizeof(fmpz_poly_struct));
  S->c = calloc(m, sizeof(long));
  S->l = calloc(m, sizeof(long));
  for (long k = 0; k < m; k++) {
    fmpz_poly_init(S->A + k);
  }
}

static void tri_set(tri_t *S, const long k, const long c, const long l,
                    const long dA, const long *A) {
  S->c[k] = c;
  S->l[k] = l;
  for (long i = 0; i <= dA; i++) {
    fmpz_poly_set_coeff_si(S->A + k, i, A[i]);
  }
}

static void tri_clear(tri_t *S) {
  fmpz_poly_clear(S->G);
  for (long k = 0; k < S->m; k++) {
    fmpz_poly_clear(S->A + k);
  }
  free(S->A);
  free(S->c);
  free(S->l);
}

/* the terms of an fmpz_mpoly in the nv variables */
static void poly_from_mpoly(hensel_poly_t *P, const fmpz_mpoly_t a,
                            const fmpz_mpoly_ctx_t ctx, const long nv) {
  P->mon = NULL;
  P->nt = fmpz_mpoly_length(a, ctx);
  P->c = _fmpz_vec_init(P->nt);
  P->e = malloc(FLINT_MAX(P->nt, 1) * nv * sizeof(long));
  ulong *e = malloc(nv * sizeof(ulong));
  for (long t = 0; t < P->nt; t++) {
    fmpz_mpoly_get_term_coeff_fmpz(P->c + t, a, t, ctx);
    fmpz_mpoly_get_term_exp_ui(e, a, t, ctx);
    for (long k = 0; k < nv; k++) {
      P->e[t * nv + k] = e[k];
    }
  }
  free(e);
}

/* a = f(s) by Horner's rule */
static void mpoly_compose(fmpz_mpoly_t a, const fmpz_poly_t f,
                          const fmpz_mpoly_t s, const fmpz_mpoly_ctx_t ctx) {
  fmpz_mpoly_zero(a, ctx);
  for (long k = fmpz_poly_degree(f); k >= 0; k--) {
    fmpz_mpoly_mul(a, a, s, ctx);
    fmpz_mpoly_add_fmpz(a, a, f->coeffs + k, ctx);
  }
}

/* the equations of S in hs->all, then extra products of the first one by
 * powers of x_0 + t, so that the system is overdetermined; the equations
 * of Newton's iteration are left to hensel_sys_set */
static void sys_from_tri(hensel_sys_t *hs, const tri_t *S, const long extra) {
  const long m = S->m, nv = m + 1;
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_ctx_init(ctx, nv, ORD_LEX);
  fmpz_mpoly_t s, x, a, b;
  fmpz_mpoly_init(s, ctx);
  fmpz_mpoly_init(x, ctx);
  fmpz_mpoly_init(a, ctx);
  fmpz_mpoly_init(b, ctx);
  /* s = t - sum l_k x_k */
  fmpz_mpoly_gen(s, m, ctx);
  for (long k = 0; k < m; k++) {
    fmpz_mpoly_gen(x, k, ctx);
    fmpz_mpoly_scalar_mul_si(x, x, S->l[k], ctx);
    fmpz_mpoly_sub(s, s, x, ctx);
  }
  hs->nv = nv;
  hs->nall = nv + extra;
  hs->all = malloc(hs->nall * sizeof(hensel_poly_t));
  hs->sys = malloc(nv * sizeof(hensel_poly_t));
  hs->jac = malloc(nv * nv * sizeof(hensel_poly_t));
  hs->maxe = calloc(nv, sizeof(long));
  for (long k = 0; k < m; k++) {
    mpoly_compose(a, S->A + k, s, ctx);
    fmpz_mpoly_gen(x, k, ctx);
    fmpz_mpoly_scalar_mul_si(x, x, S->c[k], ctx);
    fmpz_mpoly_sub(a, x, a, ctx);
    poly_from_mpoly(hs->all + k, a, ctx, nv);
  }
  mpoly_compose(a, S->G, s, ctx);
  poly_from_mpoly(hs->all + m, a, ctx, nv);
  /* (x_0 + t)^e times the first equation */
  fmpz_mpoly_gen(x, 0, ctx);
  fmpz_mpoly_gen(b, m, ctx);
  fmpz_mpoly_add(b, b, x, ctx);
  mpoly_compose(a, S->A, s, ctx);
  fmpz_mpoly_scalar_mul_si(x, x, S->c[0], ctx);
  fmpz_mpoly_sub(a, x, a, ctx);
  for (long i = 0; i < extra; i++) {
    fmpz_mpoly_mul(a, a, b, ctx);
    poly_from_mpoly(hs->all + nv + i, a, ctx, nv);
  }
  for (long i = 0; i < hs->nall; i++) {
    const hensel_poly_t *P = hs->all + i;
    for (long t = 0; t < P->nt * nv; t++) {
      hs->maxe[t % nv] = FLINT_MAX(hs->maxe[t % nv], P->e[t]);
    }
  }
  fmpz_mpoly_clear(s, ctx);
  fmpz_mpoly_clear(x, ctx);
  fmpz_mpoly_clear(a, ctx);
  fmpz_mpoly_clear(b, ctx);
  fmpz_mpoly_ctx_clear(ctx);
}

/* the parametrization of S modulo p, in P as hensel_add_param sets it:
 * in (Z/p)[s]/(G), t^0, .., t^(d-1) form a basis since t is separating,
 * w is the minimal polynomial of t and X_k the coordinates of x_k in this
 * basis; returns 0 if p divides a denominator or t does not separate */
static int param_mod_p(hensel_prime_t *P, const tri_t *S, const mp_limb_t p) {
  const long m = S->m, nv = m + 1, d = fmpz_poly_degree(S->G);
  memset(P, 0, sizeof(hensel_prime_t));
  P->p = p;
  P->ng = -1;
  nmod_poly_init(P->w, p);
  P->X = malloc(nv * sizeof(nmod_poly_struct));
  for (long k = 0; k < nv; k++) {
    nmod_poly_init(P->X + k, p);
  }
  nmod_t mod;
  nmod_init(&mod, p);
  nmod_poly_t G, t, x, pw;
  nmod_poly_init(G, p);
  nmod_poly_init(t, p);
  nmod_poly_init(x, p);
  nmod_poly_init(pw, p);
  fmpz_poly_get_nmod_poly(G, S->G);
  int ok = nmod_poly_degree(G) == d;
  /* t = s + sum l_k A_k / c_k */
  nmod_poly_set_coeff_ui(t, 1, 1);
  for (long k = 0; k < m && ok; k++) {
    const mp_limb_t c = mod_si(S->c[k], p);
    ok = c != 0;
    if (ok) {
      fmpz_poly_get_nmod_poly(x, S->A + k);
      nmod_poly_rem(x, x, G);
      nmod_poly_scalar_mul_nmod(x, x, nmod_mul(nmod_inv(c, mod),
                                               mod_si(S->l[k], p),
                                               mod));
      nmod_poly_add(t, t, x);
    }
  }
  /* the columns of M are the coordinates of t^0, .., t^(d-1) */
  nmod_mat_t M, Mi;
  nmod_mat_init(M, d, d, p);
  nmod_mat_init(Mi, d, d, p);
  nmod_poly_one(pw);
  for (long j = 0; j < d && ok; j++) {
    for (long i = 0; i < d; i++) {
      nmod_mat_entry(M, i, j) = nmod_poly_get_coeff_ui(pw, i);
    }
    nmod_poly_mulmod(pw, pw, t, G);
  }
  ok = ok && nmod_mat_inv(Mi, M);
  if (ok) {
    /* w = t^d - sum y_i t^i, y = M^-1 t^d */
    for (long i = 0; i < d; i++) {
      mp_limb_t y = 0;
      for (long j = 0; j < d; j++) {
        y = nmod_add(y, nmod_mul(nmod_mat_entry(Mi, i, j),
                                 nmod_poly_get_coeff_ui(pw, j), mod), mod);
      }
      nmod_poly_set_coeff_ui(P->w, i, nmod_neg(y, mod));
    }
    nmod_poly_set_coeff_ui(P->w, d, 1);
    for (long k = 0; k < m; k++) {
      fmpz_poly_get_nmod_poly(x, S->A + k);
      nmod_poly_rem(x, x, G);
      nmod_poly_scalar_mul_nmod(x, x, nmod_inv(mod_si(S->c[k], p), mod));
      for (long i = 0; i < d; i++) {
        mp_limb_t y = 0;
        for (long j = 0; j < d; j++) {
          y = nmod_add(y, nmod_mul(nmod_mat_entry(Mi, i, j),
                                   nmod_poly_get_coeff_ui(x, j), mod), mod);
        }
        nmod_poly_set_coeff_ui(P->X + k, i, y);
      }
    }
    nmod_poly_set_coeff_ui(P->X + m, 1, 1);
  }
  nmod_mat_clear(M);
  nmod_mat_clear(Mi);
  nmod_poly_clear(G);
  nmod_poly_clear(t);
  nmod_poly_clear(x);
  nmod_poly_clear(pw);
  return ok;
}

/* the equations of hs->all at the lifted X modulo (f, p^prec) vanish */
static int lift_vanishes(const hensel_prime_t *P, const hensel_sys_t *hs) {
  const hensel_lift_t *L = P->nfL;
  const fmpz_mod_ctx_struct *ctx = P->nfctx;
  const long nv = hs->nv;
  fmpz_mod_poly_t v, acc;
  fmpz_mod_poly_init(v, ctx);
  fmpz_mod_poly_init(acc, ctx);
  int ok = 1;
  for (long i = 0; i < hs->nall && ok; i++) {
    const hensel_poly_t *Q = hs->all + i;
    fmpz_mod_poly_zero(v, ctx);
    for (long m = 0; m < Q->nt; m++) {
      fmpz_mod_poly_set_fmpz(acc, Q->c + m, ctx);
      for (long k = 0; k < nv; k++) {
        for (long e = 0; e < Q->e[m * nv + k]; e++) {
          fmpz_mod_poly_mulmod(acc, acc, L->X + k, L->f, ctx);
        }
      }
      fmpz_mod_poly_add(v, v, acc, ctx);
    }
    ok = fmpz_mod_poly_is_zero(v, ctx);
  }
  fmpz_mod_poly_clear(v, ctx);
  fmpz_mod_poly_clear(acc, ctx);
  return ok;
}

/* clears T primes of a parametrization, lifted or not */
static void prime_free(hensel_prime_t *P, const long T, const long nv) {
  hensel_data_t hd;
  hd.nprimes = T;
  hd.alloc = T;
  hd.pr = malloc(T * sizeof(hensel_prime_t));
  memcpy(hd.pr, P, T * sizeof(hensel_prime_t));
  hensel_data_clear(&hd, nv);
}

/* lifts the parametrization of S modulo the T primes after p0 to p^N, J^-1
 * being updated up to p^(N / 2^skip), with the equations of Newton's
 * iteration chosen by hensel_select (sel = 1), as the combinations of all
 * the equations (sel = 0) or with the products by transforms (sel = 2);
 * if rec is set, recombines them by CRT, reconstructs them over Q and
 * compares them with the parametrization modulo another prime. P receives
 * the images of the lifted primes */
static void test_lift(const char *name, const tri_t *S, const long extra,
                      const int sel, const long T, const mp_limb_t p0,
                      const long N, const long skip, const int rec,
                      hensel_prime_t *P) {
  hensel_sys_t hs;
  sys_from_tri(&hs, S, extra);
  const long nv = hs.nv, d = fmpz_poly_degree(S->G);
  mp_limb_t p = p0;
  for (long j = 0; j < T; j++) {
    p = n_nextprime(p, 1);
    CHECK(param_mod_p(P + j, S, p), "%s: no parametrization modulo %lu",
          name, (unsigned long)p);
  }
  long *sub = NULL;
  if (sel == 1) {
    sub = hensel_select(&hs, P);
    CHECK(sub != NULL, "%s: no square subsystem selected", name);
  }
  hs.square = sub != NULL;
  hs.ntt = sel == 2;
  hensel_sys_set(&hs, sub);
  free(sub);
  for (long j = 0; j < T; j++) {
    hensel_nf_start(P + j, &hs);
    CHECK(P[j].nfok, "%s: J not invertible modulo (w, p)", name);
    while (P[j].prec < N) {
      hensel_nf_step(P + j, &hs, N, N >> skip);
    }
    CHECK(lift_vanishes(P + j, &hs), "%s: F(X) != 0 modulo (w, p^%ld)",
          name, N);
    hensel_nf_images(P + j, &hs);
    CHECK(P[j].gdeg[0] == d, "%s: degree %ld instead of %ld", name,
          P[j].gdeg[0], d);
  }
  if (rec) {
    fmpz_t one;
    fmpz_init_set_ui(one, 1);
    hensel_crt_t C;
    hensel_crt_init(&C, one, P, T, d, nv);
    split_group_t G;
    G.deg = d;
    G.n = d * nv;
    G.num = _fmpz_vec_init(G.n);
    G.den = _fmpz_vec_init(nv);
    long nrec = 0;
    CHECK(hensel_reconstruct(&G, &C, NULL, &nrec), "%s: not reconstructed "
          "(%ld of %ld coefficients of w)", name, nrec, d);
    /* the parametrization modulo another prime: w and B_v = X_v w' mod w */
    const mp_limb_t q = n_nextprime(UWORD(1) << 29, 1);
    hensel_prime_t Q;
    CHECK(param_mod_p(&Q, S, q), "%s: no parametrization modulo %lu", name,
          (unsigned long)q);
    CHECK(hensel_check(&G, &hs, q) == 1, "%s: check modulo %lu fails",
          name, (unsigned long)q);
    nmod_poly_t wd, b;
    nmod_poly_init(wd, q);
    nmod_poly_init(b, q);
    nmod_poly_derivative(wd, Q.w);
    for (long i = 0; i < d * nv; i++) {
      mp_limb_t e;
      if (i < d) {
        e = nmod_poly_get_coeff_ui(Q.w, i);
      } else {
        nmod_poly_mulmod(b, Q.X + (i - d) / d, wd, Q.w);
        e = nmod_poly_get_coeff_ui(b, (i - d) % d);
      }
      const mp_limb_t den = fmpz_fdiv_ui(G.den + i / d, q);
      CHECK(den != 0 && fmpz_fdiv_ui(G.num + i, q)
            == n_mulmod2(e, den, q), "%s: coefficient %ld differs modulo %lu",
            name, i, (unsigned long)q);
    }
    nmod_poly_clear(wd);
    nmod_poly_clear(b);
    prime_free(&Q, 1, nv);
    _fmpz_vec_clear(G.num, G.n);
    _fmpz_vec_clear(G.den, nv);
    hensel_crt_clear(&C);
    fmpz_clear(one);
  }
  for (long j = 0; j < T; j++) {
    hensel_nf_clear(P + j, &hs);
  }
  hensel_sys_clear(&hs);
}

static void test_lifts(void) {
  hensel_prime_t P[2];
  /* 2 x - s^2 - 1 = 3 s^3 - s - 7 = 0, t = s + x */
  tri_t S;
  tri_init(&S, 1, 3, (long[]){-7, -1, 0, 3});
  tri_set(&S, 0, 2, 1, 2, (long[]){1, 0, 1});
  test_lift("square", &S, 0, 0, 1, UWORD(1) << 30, 16, 0, 1, P);
  prime_free(P, 1, 2);
  test_lift("square, two primes", &S, 0, 0, 2, UWORD(1) << 30, 8, 0, 1, P);
  prime_free(P, 2, 2);
  /* with two products of the first equation: a square subsystem, then
   * random combinations of the equations */
  test_lift("subsystem", &S, 2, 1, 1, UWORD(1) << 30, 16, 0, 1, P);
  prime_free(P, 1, 2);
  test_lift("combinations", &S, 2, 0, 1, UWORD(1) << 30, 16, 0, 1, P);
  prime_free(P, 1, 2);
  /* a small prime, and the last two steps without update of J^-1, with
   * the corrections of hensel_solve */
  test_lift("small prime, corrections", &S, 0, 0, 1, 1000, 64, 2, 1, P);
  prime_free(P, 1, 2);
  tri_clear(&S);
  /* three variables, 5 x_0 = s^4 - 3 s + 1, 3 x_1 = s^2 - 2 and G of
   * degree 70, t = s + x_0: the products by transforms (degree at
   * least HENSEL_NTT_MIND), in the updates of J^-1 and in the corrections,
   * give the same lifting as the others */
  long G[71] = {0};
  G[0] = -11;
  G[1] = 2;
  G[3] = -5;
  G[17] = 4;
  G[70] = 3;
  tri_init(&S, 2, 70, G);
  tri_set(&S, 0, 5, 1, 4, (long[]){1, -3, 0, 0, 1});
  tri_set(&S, 1, 3, 0, 2, (long[]){-2, 0, 1});
  test_lift("degree 70", &S, 0, 0, 1, UWORD(1) << 30, 16, 2, 0, P);
#if HENSEL_NTT
  hensel_prime_t Pn[1];
  test_lift("degree 70, transforms", &S, 0, 2, 1, UWORD(1) << 30, 16, 2, 0,
            Pn);
  CHECK(fmpz_equal(P->q, Pn->q) && fmpz_poly_equal(P->g, Pn->g)
        && fmpz_poly_equal(P->B, Pn->B)
        && fmpz_poly_equal(P->B + 1, Pn->B + 1),
        "degree 70: transforms give another lifting");
  prime_free(Pn, 1, 3);
#endif
  prime_free(P, 1, 3);
  tri_clear(&S);
}

int main(void) {
#if HENSEL_NTT
  test_ntt();
#endif
  test_inv();
  test_lifts();
  if (nfail > 0) {
    fprintf(stderr, "%d checks failed\n", nfail);
    return 1;
  }
  return 0;
}
