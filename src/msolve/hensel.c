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

/* Parametrization of the solution set by p-adic lifting.

   For a square system with simple solutions, the parametrization modulo a
   prime p, x_v = X_v(t) mod w(t), is lifted to p^N by Newton's iteration
   (the variable t is the last one):
     dX    = J^-1 F(X)
     Delta = t - X_t + dX_t
     w    <- w + (w' Delta mod w)
     X_v  <- X_v - dX_v + (X_v' Delta mod w)
     J^-1 <- J^-1 (2 I - J(X) J^-1)
   all in (Z/p^2N)[t]/(w), at one prime per thread. The parametrization,
   described by w and B_v = X_v w' mod w as in splitting.c, is recombined
   by CRT over the primes and reconstructed over Q, step by step as the
   precision grows, then checked by evaluating the system at it modulo a
   prime that was not used.

   A cost model, fed by the first modular runs, chooses N and decides
   whether this lifting is used instead of the multi-modular one (option
   --lifting of msolve). The options --hensel-precision, --hensel-cost-scale
   and --hensel-memory (fields hensel_* of data_gens_ff_t) fix N, scale the
   cost of the modular runs and fix the memory limit. */

#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>

/* the products of polynomial matrices by transforms (see hensel_ntt_mul)
 * use the interface of fft_small of FLINT 3.2.0 and later; HENSEL_NTT 0,
 * or older versions, keep the products of fmpz_poly_mul */
#ifndef HENSEL_NTT
#if __FLINT_RELEASE >= 30200 && defined(FLINT_HAVE_FFT_SMALL) \
  && FLINT_HAVE_FFT_SMALL
#define HENSEL_NTT 1
#else
#define HENSEL_NTT 0
#endif
#endif
#if HENSEL_NTT
#include <flint/fmpz.h>
#include <flint/nmod.h>
#include <flint/fft_small.h>
#endif

#define HENSEL_RRMARGIN 20        /* bits kept by rational reconstruction */

typedef struct {
  long nt;
  fmpz *c;
  long *e;                        /* nt * nv exponents */
  long *mon;                      /* index of each term in the monomials */
} hensel_poly_t;

/* the monomials of the equations and of their derivatives, closed under
 * m -> m / x_var[m] = pred[m]: each one is evaluated with one product;
 * ord[w] lists those needed by the equations (w = 0) or by the Jacobian
 * (w = 1), by increasing degree; nlive[w] values of them are held at once
 * by hensel_eval_set. The nleaf[w] products of degree at least 2 without
 * successor (eco 12: 56 of the 79 monomials of F, katsura 12: 77 of 90)
 * are only added to the nred[w] polynomials of their tleaf[w] terms: if
 * delay[w] is set, they are not reduced modulo w, and each of these
 * polynomials is reduced once, which saves the reduction of nleaf[w] -
 * nred[w] products (eco 12, J: 66 products for 95 polynomials) */
typedef struct {
  long nmon, alloc;
  long *e;                        /* nmon * nv exponents */
  long *deg, *pred, *var;
  long nord[2], nlive[2];
  long nleaf[2], nred[2], tleaf[2];
  int delay[2];
  long *ord[2];
} hensel_mon_t;

typedef struct {
  long nv, nall;
  hensel_poly_t *all;             /* nall >= nv equations with integer
                                     coefficients, for the checks */
  hensel_poly_t *sys;             /* nv of them or nv combinations of
                                     them, for Newton */
  hensel_poly_t *jac;             /* jac[i * nv + j] = d sys[i] / d x_j */
  long *maxe;                     /* largest exponent of each variable */
  hensel_mon_t mt;
  long N;                         /* precision p^N */
  long skip;                      /* last steps without update of J^-1 */
  int nauto;                      /* N chosen by the cost model */
  double cnewton, cmod;           /* model costs of one prime (bit ops) */
  double gain;                    /* precision per cost relative to the
                                     modular runs alone */
  int square;                     /* sys is a subsystem of all */
  int ntt;                        /* products of hensel_solve and
                                     hensel_matmul by transforms (see
                                     hensel_ntt_mul), chosen by
                                     hensel_choose and hensel_nf_nmem */
  int nvar;                       /* 2 if they may be used, else 1 */
  long vN[2], vskip[2];           /* N and skip without (0) and with them */
  double tsel;                    /* time of its choice (sec) */
} hensel_sys_t;

typedef struct {
  mp_limb_t p;
  nmod_poly_t w;                  /* monic eliminating polynomial mod p */
  nmod_poly_struct *X;            /* x_v = X_v(t) mod w, X_{nv-1} = t */
  /* images at the current precision */
  fmpz_t q;                       /* p^prec */
  long ng;                        /* 1 once they are formed, -1 before */
  long *gdeg;
  fmpz_poly_struct *g;            /* monic g = w, coefficients in [0, q) */
  fmpz_poly_struct *B;            /* nv - 1 polynomials B_v = X_v g' mod g */
  double tnewton;
  double tinv, tf, tmv, tj, tmm, tup; /* parts of Newton's iteration */
  /* lifting, kept between the steps */
  void *nfL;                      /* hensel_lift_t, NULL if not started */
  fmpz_mod_ctx_struct nfctx[1];
  long prec;                      /* current precision p^prec */
  int nfok;                       /* J is invertible modulo (w, p) */
} hensel_prime_t;

typedef struct {
  long nprimes, alloc;
  hensel_prime_t *pr;
} hensel_data_t;

#if HENSEL_NTT
/* the products A B over Z[t] of the matrices of hensel_solve and
 * hensel_matmul are formed from the transforms of their entries modulo np
 * primes of 50 bits (those of fft_small, whose product exceeds the
 * coefficients of A B), so that those of A are computed once for all the
 * columns of B and for all the products by A (J^-1 is the same in all the
 * steps without update, J in all the rounds of hensel_solve): a product
 * of an nv x nv matrix by a vector costs nv transforms, nv^2 products of
 * transforms and nv inverse transforms, instead of nv^2 products of
 * polynomials (3 nv^2 transforms); katsura 12, D = 2046, 124 x 500 bits:
 * 27 ms instead of 145 ms. The transforms have the length 2^L of the
 * products, their coefficients are doubles in [-p/2, p/2], and those of A
 * are multiplied by 2^-L, the factor of the inverse transform */
#define HENSEL_NTT_MAXP 64        /* at most 64 primes, about 3190 bits */
#define HENSEL_NTT_MIND 64        /* degree from which they are used */

typedef struct {
  fmpz_comb_t comb;
  fmpz_comb_temp_t ct;
  ulong *cof;                     /* np x np limbs: P / p_k */
  ulong *P;                       /* the product P of the primes */
  ulong inv[HENSEL_NTT_MAXP];     /* (P / p_k)^-1 mod p_k */
} hensel_ntt_crt_t;

typedef struct {
  long np, L;
  ulong N, trunc;                 /* 2^L, and the length used, a multiple
                                     of BLK_SZ */
  ulong pr[HENSEL_NTT_MAXP];
  sd_fft_ctx_struct Q[HENSEL_NTT_MAXP];
  nmod_t mod[HENSEL_NTT_MAXP];
  double scale[HENSEL_NTT_MAXP];  /* 2^-L mod p */
  double lbits[HENSEL_NTT_MAXP];  /* log2 of the products of the primes */
  hensel_ntt_crt_t *crt[HENSEL_NTT_MAXP + 1]; /* for the first np primes */
  double *buf;                    /* N doubles */
} hensel_ntt_t;

/* transforms of the r x c entries of a matrix, modulo its np primes */
typedef struct {
  long r, c, np;
  double *d;                      /* entry (i, j), prime k: d + ((i c + j)
                                     np + k) trunc */
} hensel_ntt_mat_t;

static ulong hensel_ntt_next(ulong p) {
  /* the primes of fft_small (next_fft_number in fft_small/mpn_mul.c) */
  const ulong bits = n_nbits(p), l = n_trailing_zeros(p - 1);
  const ulong q = p - (UWORD(2) << l);
  if (n_nbits(q) == bits) {
    return q;
  }
  if (l < 5) {
    return (UWORD(1) << (bits - 2)) + 1;
  }
  return (UWORD(1) << bits) - (UWORD(1) << (l - 1)) + 1;
}

static void hensel_ntt_init(hensel_ntt_t *T) {
  memset(T, 0, sizeof(hensel_ntt_t));
}

static void hensel_ntt_clear(hensel_ntt_t *T) {
  for (long k = 0; k < T->np; k++) {
    sd_fft_ctx_clear(T->Q + k);
  }
  for (long n = 0; n <= HENSEL_NTT_MAXP; n++) {
    hensel_ntt_crt_t *C = T->crt[n];
    if (C != NULL) {
      fmpz_comb_temp_clear(C->ct);
      fmpz_comb_clear(C->comb);
      free(C->cof);
      free(C->P);
      free(C);
    }
  }
  flint_aligned_free(T->buf);
  hensel_ntt_init(T);
}

/* the number of primes for coefficients of at most bits bits, 0 if more
 * than HENSEL_NTT_MAXP, after the transforms of length the products of
 * polynomials of length D are set up */
static long hensel_ntt_fit(hensel_ntt_t *T, const long D, const double bits) {
  const long L = FLINT_MAX(n_clog2(2 * D - 1), LG_BLK_SZ);
  if (T->buf == NULL || T->L != L) {
    for (long k = 0; k < T->np; k++) {
      sd_fft_ctx_fit_depth(T->Q + k, L);
      T->scale[k] = (double)nmod_inv(nmod_pow_ui(2, L, T->mod[k]),
                                     T->mod[k]);
    }
    T->L = L;
    T->N = UWORD(1) << L;
    T->trunc = FLINT_MIN(T->N, (ulong)((2 * D - 1 + BLK_SZ - 1) / BLK_SZ)
                         * BLK_SZ);
    flint_aligned_free(T->buf);
    T->buf = flint_aligned_alloc(64, T->N * sizeof(double));
  }
  long np = 0;
  while (np < HENSEL_NTT_MAXP && (np == 0 || T->lbits[np - 1] <= bits + 1)) {
    if (np == T->np) {
      ulong p = np == 0 ? UWORD(0x0003f00000000001) : T->pr[np - 1];
      do {
        p = hensel_ntt_next(p);
      } while (!n_is_prime(p) || n_trailing_zeros(p - 1) < 24);
      T->pr[np] = p;
      sd_fft_ctx_init_prime(T->Q + np, p);
      sd_fft_ctx_fit_depth(T->Q + np, L);
      nmod_init(T->mod + np, p);
      T->scale[np] = (double)nmod_inv(nmod_pow_ui(2, L, T->mod[np]),
                                      T->mod[np]);
      T->lbits[np] = (np > 0 ? T->lbits[np - 1] : 0) + log2((double)p);
      T->np++;
    }
    np++;
  }
  if (T->lbits[np - 1] <= bits + 1) {
    return 0;
  }
  if (T->crt[np] == NULL) {
    hensel_ntt_crt_t *C = malloc(sizeof(hensel_ntt_crt_t));
    fmpz_comb_init(C->comb, T->pr, np);
    fmpz_comb_temp_init(C->ct, C->comb);
    C->P = calloc(np + 1, sizeof(ulong));
    C->cof = calloc(np * np, sizeof(ulong));
    C->P[0] = 1;
    long pl = 1;
    for (long k = 0; k < np; k++) {
      C->P[pl] = mpn_mul_1(C->P, C->P, pl, T->pr[k]);
      pl += C->P[pl] != 0;
    }
    for (long k = 0; k < np; k++) {
      mpn_divrem_1(C->cof + k * np, 0, C->P, np, T->pr[k]);
      C->inv[k] = nmod_inv(mpn_mod_1(C->cof + k * np, np, T->pr[k]),
                           T->mod[k]);
    }
    T->crt[np] = C;
  }
  return np;
}

/* a b mod p in [-p/2, p/2] for doubles |a|, |b| < p, four at a time */
FLINT_FORCE_INLINE vec4d hensel_ntt_mulmod(vec4d a, vec4d b, vec4d p,
                                           vec4d pinv) {
  const vec4d h = vec4d_mul(a, b), l = vec4d_fmsub(a, b, h);
  const vec4d r = vec4d_add(vec4d_fnmadd(vec4d_round(vec4d_mul(h, pinv)),
                                         p, h), l);
  return vec4d_fnmadd(vec4d_round(vec4d_mul(r, pinv)), p, r);
}

FLINT_FORCE_INLINE vec4d hensel_ntt_red(vec4d r, vec4d p, vec4d pinv) {
  return vec4d_fnmadd(vec4d_round(vec4d_mul(r, pinv)), p, r);
}

/* the transforms of a, modulo the first np primes, to out (np blocks of
 * trunc doubles), multiplied by 2^-L if scale is set */
static void hensel_ntt_fwd(double *out, const fmpz_poly_t a, const long np,
                           const int scale, hensel_ntt_t *T) {
  const long len = a->length;
  const ulong N = T->N, trunc = T->trunc;
  hensel_ntt_crt_t *C = T->crt[np];
  ulong *res = malloc(np * sizeof(ulong));
  for (long k = 0; k < np; k++) {
    memset(out + k * trunc, 0, trunc * sizeof(double));
  }
  for (long c = 0; c < len; c++) {
    if (fmpz_is_zero(a->coeffs + c)) {
      continue;
    }
    fmpz_multi_mod_ui(res, a->coeffs + c, C->comb, C->ct);
    for (long k = 0; k < np; k++) {
      const double x = (double)res[k], p = T->Q[k].p;
      out[k * trunc + c] = x > p / 2 ? x - p : x;
    }
  }
  free(res);
  for (long k = 0; k < np && len > 0; k++) {
    double *d = T->buf;
    memcpy(d, out + k * trunc, trunc * sizeof(double));
    if (N > trunc) {
      memset(d + trunc, 0, (N - trunc) * sizeof(double));
    }
    sd_fft_trunc(T->Q + k, d, T->L, (len + BLK_SZ - 1) / BLK_SZ * BLK_SZ,
                 trunc);
    const vec4d p = vec4d_set_d(T->Q[k].p), pinv = vec4d_set_d(T->Q[k].pinv);
    const vec4d s = vec4d_set_d(T->scale[k]);
    double *o = out + k * trunc;
    for (ulong x = 0; x < trunc; x += 4) {
      const vec4d v = vec4d_load_aligned(d + x);
      vec4d_store_aligned(o + x, scale ? hensel_ntt_mulmod(v, s, p, pinv)
                          : hensel_ntt_red(v, p, pinv));
    }
  }
}

/* the transforms of the r x c matrix A, for products by matrices with
 * coefficients of at most bB bits and c rows; 0 if they would need more
 * than HENSEL_NTT_MAXP primes */
static int hensel_ntt_mat_set(hensel_ntt_mat_t *M, const fmpz_poly_struct *A,
                              const long r, const long c, const double bB,
                              const long D, hensel_ntt_t *T) {
  long bA = 0;
  for (long i = 0; i < r * c; i++) {
    bA = FLINT_MAX(bA, FLINT_ABS(_fmpz_vec_max_bits(A[i].coeffs,
                                                    A[i].length)));
  }
  const long np = hensel_ntt_fit(T, D, bA + bB + log2((double)c * D) + 1);
  M->r = r;
  M->c = c;
  M->np = np;
  M->d = NULL;
  if (np == 0) {
    return 0;
  }
  M->d = flint_aligned_alloc(64, r * c * np * T->trunc * sizeof(double));
  for (long i = 0; i < r * c; i++) {
    hensel_ntt_fwd(M->d + i * np * T->trunc, A + i, np, 1, T);
  }
  return 1;
}

static void hensel_ntt_mat_clear(hensel_ntt_mat_t *M) {
  flint_aligned_free(M->d);
  M->d = NULL;
  M->np = 0;
}

/* x = the integer in [0, P) with the residues res modulo the np primes */
static void hensel_ntt_crt(fmpz_t x, const ulong *res, const long np,
                           const hensel_ntt_t *T, ulong *t) {
  const hensel_ntt_crt_t *C = T->crt[np];
  memset(t, 0, (np + 1) * sizeof(ulong));
  for (long k = 0; k < np; k++) {
    t[np] += mpn_addmul_1(t, C->cof + k * np, np,
                          nmod_mul(res[k], C->inv[k], T->mod[k]));
  }
  /* t < np P */
  while (t[np] != 0 || mpn_cmp(t, C->P, np) >= 0) {
    t[np] -= mpn_sub_n(t, t, C->P, np);
  }
  long n = np;
  while (n > 0 && t[n - 1] == 0) {
    n--;
  }
  fmpz_set_ui_array(x, t, FLINT_MAX(n, 1));
}

/* C = A B over Z[t] for the transforms M of A (r x c) and B (c x n): the
 * transforms of each column of B are multiplied by those of A and added */
static void hensel_ntt_mul(fmpz_poly_struct *C, const hensel_ntt_mat_t *M,
                           const fmpz_poly_struct *B, const long n,
                           const long D, hensel_ntt_t *T) {
  const long r = M->r, c = M->c, np = M->np;
  const ulong trunc = T->trunc, zl = 2 * D - 1;
  double *Bh = flint_aligned_alloc(64, c * np * trunc * sizeof(double));
  double *out = flint_aligned_alloc(64, np * trunc * sizeof(double));
  ulong *res = malloc(np * sizeof(ulong)), *t = malloc((np + 1)
                                                         * sizeof(ulong));
  for (long l = 0; l < n; l++) {
    for (long j = 0; j < c; j++) {
      hensel_ntt_fwd(Bh + j * np * trunc, B + j * n + l, np, 0, T);
    }
    for (long i = 0; i < r; i++) {
      for (long k = 0; k < np; k++) {
        const vec4d p = vec4d_set_d(T->Q[k].p);
        const vec4d pinv = vec4d_set_d(T->Q[k].pinv);
        double *acc = T->buf;
        memset(acc, 0, T->N * sizeof(double));
        for (long j = 0; j < c; j++) {
          const double *a = M->d + ((i * c + j) * np + k) * trunc;
          const double *b = Bh + (j * np + k) * trunc;
          /* |acc| < 8 p / 2 < 2^53 before each reduction */
          const int red = (j & 7) == 7 || j == c - 1;
          for (ulong x = 0; x < trunc; x += 4) {
            vec4d s = vec4d_add(vec4d_load_aligned(acc + x),
                                hensel_ntt_mulmod(vec4d_load_aligned(a + x),
                                                  vec4d_load_aligned(b + x),
                                                  p, pinv));
            vec4d_store_aligned(acc + x, red ? hensel_ntt_red(s, p, pinv) : s);
          }
        }
        sd_ifft_trunc(T->Q + k, acc, T->L, trunc);
        const double pk = T->Q[k].p, pki = T->Q[k].pinv;
        for (ulong x = 0; x < zl; x++) {
          double y = acc[x] - pk * rint(acc[x] * pki);
          out[k * trunc + x] = y < 0 ? y + pk : y;
        }
      }
      fmpz_poly_struct *R = C + i * n + l;
      fmpz_poly_fit_length(R, zl);
      for (ulong x = 0; x < zl; x++) {
        for (long k = 0; k < np; k++) {
          res[k] = (ulong)out[k * trunc + x];
        }
        hensel_ntt_crt(R->coeffs + x, res, np, T, t);
      }
      _fmpz_poly_set_length(R, zl);
      _fmpz_poly_normalise(R);
    }
  }
  free(t);
  free(res);
  flint_aligned_free(out);
  flint_aligned_free(Bh);
}
#endif

/* lifting modulo q: w, its inverse for mulmod_preinv, X and J^-1 */
typedef struct {
  fmpz_mod_poly_t f, finv;
  fmpz_mod_poly_struct *X;
  fmpz_mod_poly_struct *Jinv;
  long jprec;                     /* J^-1 is known modulo p^jprec */
#if HENSEL_NTT
  hensel_ntt_t *ntt;
  hensel_ntt_mat_t Jh;            /* transforms of J^-1, while it is not
                                     updated (Jh.d NULL if not formed) */
#endif
  /* timings: F(X), J^-1 F(X) with the corrections, J(X), the two
   * products of the update of J^-1, the update of X and w */
  double tf, tmv, tj, tmm, tup;
} hensel_lift_t;

static void hensel_poly_clear(hensel_poly_t *P) {
  _fmpz_vec_clear(P->c, P->nt);
  free(P->e);
  free(P->mon);
}

/* index of the monomial e in the open-addressing table h of size hsz
 * (a power of 2), inserted if absent */
static long hensel_mon_find(hensel_mon_t *mt, long *h, const long hsz,
                            const long *e, const long nv) {
  uint64_t k = 0;
  for (long i = 0; i < nv; i++) {
    k = (k ^ (uint64_t)e[i]) * 0x100000001b3ULL;
  }
  for (uint64_t j = k & (hsz - 1);; j = (j + 1) & (hsz - 1)) {
    if (h[j] < 0) {
      if (mt->nmon == mt->alloc) {
        mt->alloc *= 2;
        mt->e = realloc(mt->e, mt->alloc * nv * sizeof(long));
      }
      memcpy(mt->e + mt->nmon * nv, e, nv * sizeof(long));
      h[j] = mt->nmon++;
      return h[j];
    }
    if (memcmp(mt->e + h[j] * nv, e, nv * sizeof(long)) == 0) {
      return h[j];
    }
  }
}

static int hensel_mon_lookup(const hensel_mon_t *mt, const long *h,
                             const long hsz, const long *e, const long nv) {
  uint64_t k = 0;
  for (long i = 0; i < nv; i++) {
    k = (k ^ (uint64_t)e[i]) * 0x100000001b3ULL;
  }
  for (uint64_t j = k & (hsz - 1); h[j] >= 0; j = (j + 1) & (hsz - 1)) {
    if (memcmp(mt->e + h[j] * nv, e, nv * sizeof(long)) == 0) {
      return 1;
    }
  }
  return 0;
}

static void hensel_mon_init(hensel_sys_t *hs) {
  const long nv = hs->nv;
  hensel_mon_t *mt = &hs->mt;
  long nterms = 0, maxdeg = 0;
  for (long i = 0; i < nv; i++) {
    nterms += hs->sys[i].nt;
    for (long m = 0; m < hs->sys[i].nt; m++) {
      long d = 0;
      for (long k = 0; k < nv; k++) {
        d += hs->sys[i].e[m * nv + k];
      }
      maxdeg = MAX(maxdeg, d);
    }
  }
  for (long i = 0; i < nv * nv; i++) {
    nterms += hs->jac[i].nt;
  }
  /* at most maxdeg monomials are added for each term */
  long hsz = 16;
  while (hsz < 2 * nterms * (maxdeg + 1)) {
    hsz *= 2;
  }
  long *h = malloc(hsz * sizeof(long));
  for (long j = 0; j < hsz; j++) {
    h[j] = -1;
  }
  mt->nmon = 0;
  mt->alloc = 64;
  mt->e = malloc(mt->alloc * nv * sizeof(long));
  for (long w = 0; w < 2; w++) {
    const long np = w == 0 ? nv : nv * nv;
    hensel_poly_t *Q = w == 0 ? hs->sys : hs->jac;
    for (long i = 0; i < np; i++) {
      hensel_poly_t *P = Q + i;
      P->mon = malloc(MAX(P->nt, 1) * sizeof(long));
      for (long m = 0; m < P->nt; m++) {
        P->mon[m] = hensel_mon_find(mt, h, hsz, P->e + m * nv, nv);
      }
    }
  }
  /* predecessors, preferably among the monomials already there */
  long *e = malloc(nv * sizeof(long));
  mt->pred = NULL;
  mt->var = NULL;
  long palloc = 0;
  for (long done = 0; done < mt->nmon; done++) {
    if (palloc < mt->alloc) {
      palloc = mt->alloc;
      mt->pred = realloc(mt->pred, palloc * sizeof(long));
      mt->var = realloc(mt->var, palloc * sizeof(long));
    }
    memcpy(e, mt->e + done * nv, nv * sizeof(long));
    long kbest = -1;
    for (long k = 0; k < nv && kbest < 0; k++) {
      if (e[k] > 0) {
        e[k]--;
        if (hensel_mon_lookup(mt, h, hsz, e, nv)) {
          kbest = k;
        }
        e[k]++;
      }
    }
    for (long k = 0; k < nv && kbest < 0; k++) {
      if (e[k] > 0) {
        kbest = k;
      }
    }
    mt->var[done] = kbest;
    mt->pred[done] = -1;
    if (kbest >= 0) {
      e[kbest]--;
      mt->pred[done] = hensel_mon_find(mt, h, hsz, e, nv);
    }
  }
  free(e);
  free(h);
  mt->deg = malloc(mt->nmon * sizeof(long));
  for (long m = 0; m < mt->nmon; m++) {
    mt->deg[m] = 0;
    for (long k = 0; k < nv; k++) {
      mt->deg[m] += mt->e[m * nv + k];
    }
  }
  /* the monomials needed by each set, with their predecessors, by
   * increasing degree */
  long maxd = 0;
  for (long m = 0; m < mt->nmon; m++) {
    maxd = MAX(maxd, mt->deg[m]);
  }
  char *need = malloc(mt->nmon);
  for (long w = 0; w < 2; w++) {
    memset(need, 0, mt->nmon);
    const long np = w == 0 ? nv : nv * nv;
    const hensel_poly_t *Q = w == 0 ? hs->sys : hs->jac;
    for (long i = 0; i < np; i++) {
      for (long m = 0; m < Q[i].nt; m++) {
        need[Q[i].mon[m]] = 1;
      }
    }
    for (long d = maxd; d > 0; d--) {
      for (long m = 0; m < mt->nmon; m++) {
        if (need[m] && mt->deg[m] == d) {
          need[mt->pred[m]] = 1;
        }
      }
    }
    mt->ord[w] = malloc(mt->nmon * sizeof(long));
    mt->nord[w] = 0;
    for (long d = 0; d <= maxd; d++) {
      for (long m = 0; m < mt->nmon; m++) {
        if (need[m] && mt->deg[m] == d) {
          mt->ord[w][mt->nord[w]++] = m;
        }
      }
    }
    /* values held by hensel_eval_set: each one until its last successor
     * is computed */
    long *nsucc = calloc(mt->nmon, sizeof(long));
    for (long j = 0; j < mt->nord[w]; j++) {
      const long m = mt->ord[w][j];
      if (mt->pred[m] >= 0) {
        nsucc[mt->pred[m]]++;
      }
    }
    mt->nleaf[w] = mt->nred[w] = mt->tleaf[w] = 0;
    for (long j = 0; j < mt->nord[w]; j++) {
      const long m = mt->ord[w][j];
      mt->nleaf[w] += nsucc[m] == 0 && mt->deg[m] >= 2;
    }
    for (long i = 0; i < np; i++) {
      long c = 0;
      for (long m = 0; m < Q[i].nt; m++) {
        const long k = Q[i].mon[m];
        c += nsucc[k] == 0 && mt->deg[k] >= 2;
      }
      mt->tleaf[w] += c;
      mt->nred[w] += c > 0;
    }
    mt->delay[w] = mt->nleaf[w] > mt->nred[w];
    long live = 0;
    mt->nlive[w] = 0;
    for (long j = 0; j < mt->nord[w]; j++) {
      const long m = mt->ord[w][j];
      live++;
      mt->nlive[w] = MAX(mt->nlive[w], live);
      if (mt->pred[m] >= 0 && --nsucc[mt->pred[m]] == 0) {
        live--;
      }
      if (nsucc[m] == 0) {
        live--;
      }
    }
    free(nsucc);
  }
  free(need);
}

static void hensel_mon_clear(hensel_mon_t *mt) {
  free(mt->e);
  free(mt->deg);
  free(mt->pred);
  free(mt->var);
  free(mt->ord[0]);
  free(mt->ord[1]);
}

/* sorts the terms of P by exponents and adds those with the same
 * monomial, removing the zero ones */
static long hensel_merge_nv;
static int hensel_merge_cmp(const void *a, const void *b) {
  const long *ea = *(long *const *)a, *eb = *(long *const *)b;
  for (long k = 0; k < hensel_merge_nv; k++) {
    if (ea[k] != eb[k]) {
      return ea[k] < eb[k] ? -1 : 1;
    }
  }
  return 0;
}

static void hensel_poly_merge(hensel_poly_t *P, const long nv) {
  long **idx = malloc(MAX(P->nt, 1) * sizeof(long *));
  for (long m = 0; m < P->nt; m++) {
    idx[m] = P->e + m * nv;
  }
  hensel_merge_nv = nv;
  qsort(idx, P->nt, sizeof(long *), hensel_merge_cmp);
  fmpz *c = _fmpz_vec_init(MAX(P->nt, 1));
  long *e = malloc(MAX(P->nt, 1) * nv * sizeof(long));
  long n = 0;
  for (long m = 0; m < P->nt; m++) {
    const long t = (idx[m] - P->e) / nv;
    if (n > 0 && memcmp(e + (n - 1) * nv, idx[m], nv * sizeof(long)) == 0) {
      fmpz_add(c + n - 1, c + n - 1, P->c + t);
    } else {
      if (n > 0 && fmpz_is_zero(c + n - 1)) {
        n--;
      }
      fmpz_set(c + n, P->c + t);
      memcpy(e + n * nv, idx[m], nv * sizeof(long));
      n++;
    }
  }
  if (n > 0 && fmpz_is_zero(c + n - 1)) {
    n--;
  }
  _fmpz_vec_clear(P->c, P->nt);
  free(P->e);
  free(idx);
  P->c = c;
  P->e = e;
  P->nt = n;
}

/* the system of gens with integer coefficients, the equations of
 * Newton's iteration being set by hensel_sys_set; returns 0 if there are
 * fewer equations than variables */
static int hensel_sys_init(hensel_sys_t *hs, const data_gens_ff_t *gens) {
  const long nv = gens->nvars, nall = gens->ngens;
  if (nall < nv || gens->field_char != 0) {
    return 0;
  }
  hs->nv = nv;
  hs->nall = nall;
  hs->all = malloc(nall * sizeof(hensel_poly_t));
  hs->sys = malloc(nv * sizeof(hensel_poly_t));
  hs->jac = malloc(nv * nv * sizeof(hensel_poly_t));
  hs->maxe = calloc(nv, sizeof(long));
  fmpz_t l, num, den;
  fmpz_init(l);
  fmpz_init(num);
  fmpz_init(den);
  long pos = 0;
  for (long i = 0; i < nall; i++) {
    hensel_poly_t *P = hs->all + i;
    P->mon = NULL;
    P->nt = gens->lens[i];
    P->c = _fmpz_vec_init(P->nt);
    P->e = malloc(P->nt * nv * sizeof(long));
    /* l = lcm of the denominators */
    fmpz_one(l);
    for (long m = 0; m < P->nt; m++) {
      fmpz_set_mpz(den, *(gens->mpz_cfs[2 * (pos + m) + 1]));
      fmpz_lcm(l, l, den);
    }
    for (long m = 0; m < P->nt; m++) {
      fmpz_set_mpz(num, *(gens->mpz_cfs[2 * (pos + m)]));
      fmpz_set_mpz(den, *(gens->mpz_cfs[2 * (pos + m) + 1]));
      fmpz_divexact(den, l, den);
      fmpz_mul(P->c + m, num, den);
      for (long k = 0; k < nv; k++) {
        P->e[m * nv + k] = gens->exps[(pos + m) * nv + k];
        hs->maxe[k] = MAX(hs->maxe[k], P->e[m * nv + k]);
      }
    }
    pos += P->nt;
  }
  fmpz_clear(l);
  fmpz_clear(num);
  fmpz_clear(den);
  return 1;
}

/* the equations of Newton's iteration, their Jacobian and their monomials:
 * the equations sub[0], .., sub[nv - 1] of the system or, if sub is NULL,
 * the nv combinations F_i + sum_j r_ij F_j (j >= nv) with random r_ij,
 * whose Jacobian is invertible at the simple solutions of F for generic
 * r_ij */
static void hensel_sys_set(hensel_sys_t *hs, const long *sub) {
  const long nv = hs->nv, nall = hs->nall;
  /* the combinations, with a fixed generator so that the random numbers
   * of msolve are not consumed */
  uint64_t seed = 0x9e3779b97f4a7c15ULL;
  for (long i = 0; i < nv; i++) {
    hensel_poly_t *P = hs->sys + i;
    P->mon = NULL;
    if (sub != NULL) {
      const hensel_poly_t *A = hs->all + sub[i];
      P->nt = A->nt;
      P->c = _fmpz_vec_init(P->nt);
      P->e = malloc(P->nt * nv * sizeof(long));
      _fmpz_vec_set(P->c, A->c, P->nt);
      memcpy(P->e, A->e, P->nt * nv * sizeof(long));
      continue;
    }
    P->nt = hs->all[i].nt;
    for (long j = nv; j < nall; j++) {
      P->nt += hs->all[j].nt;
    }
    P->c = _fmpz_vec_init(P->nt);
    P->e = malloc(P->nt * nv * sizeof(long));
    long m = 0;
    for (long j = i; j < nall; j = (j == i ? nv : j + 1)) {
      const hensel_poly_t *A = hs->all + j;
      ulong r = 1;
      if (j != i) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        r = 1 + (seed >> 54);
      }
      for (long l = 0; l < A->nt; l++, m++) {
        fmpz_mul_ui(P->c + m, A->c + l, r);
        memcpy(P->e + m * nv, A->e + l * nv, nv * sizeof(long));
      }
    }
    /* the equations share monomials (cp_d_3_n_7_p_2: 24987 terms for 330
     * monomials, phuoc1: 88202 terms of about 2000 bits for 2147) */
    if (nall > nv) {
      hensel_poly_merge(P, nv);
    }
  }
  for (long i = 0; i < nv; i++) {
    const hensel_poly_t *S = hs->sys + i;
    for (long j = 0; j < nv; j++) {
      hensel_poly_t *P = hs->jac + i * nv + j;
      P->mon = NULL;
      long cnt = 0;
      for (long m = 0; m < S->nt; m++) {
        cnt += S->e[m * nv + j] != 0;
      }
      P->nt = cnt;
      P->c = _fmpz_vec_init(cnt);
      P->e = malloc(MAX(cnt, 1) * nv * sizeof(long));
      cnt = 0;
      for (long m = 0; m < S->nt; m++) {
        if (S->e[m * nv + j] == 0) {
          continue;
        }
        fmpz_mul_ui(P->c + cnt, S->c + m, S->e[m * nv + j]);
        for (long k = 0; k < nv; k++) {
          P->e[cnt * nv + k] = S->e[m * nv + k] - (k == j);
        }
        cnt++;
      }
    }
  }
  hensel_mon_init(hs);
}

/* clears what hensel_sys_set builds */
static void hensel_sys_unset(hensel_sys_t *hs) {
  for (long i = 0; i < hs->nv; i++) {
    hensel_poly_clear(hs->sys + i);
  }
  hensel_mon_clear(&hs->mt);
  for (long i = 0; i < hs->nv * hs->nv; i++) {
    hensel_poly_clear(hs->jac + i);
  }
}

/* the row of the Jacobian of A at x = X modulo (w, p), with pw[v][k] =
 * X_v^k mod w */
static void hensel_jac_row(nmod_poly_struct *row, const hensel_poly_t *A,
                           nmod_poly_struct **pw, const nmod_poly_t w,
                           const long nv) {
  const mp_limb_t p = w->mod.n;
  nmod_poly_t m;
  nmod_poly_init(m, p);
  long *sup = malloc(nv * sizeof(long));
  for (long j = 0; j < nv; j++) {
    nmod_poly_zero(row + j);
  }
  for (long t = 0; t < A->nt; t++) {
    const long *e = A->e + t * nv;
    long s = 0;
    for (long v = 0; v < nv; v++) {
      if (e[v] > 0) {
        sup[s++] = v;
      }
    }
    const mp_limb_t c = fmpz_fdiv_ui(A->c + t, p);
    for (long a = 0; a < s; a++) {
      const long j = sup[a];
      nmod_poly_set(m, pw[j] + e[j] - 1);
      for (long b = 0; b < s; b++) {
        if (b != a) {
          nmod_poly_mulmod(m, m, pw[sup[b]] + e[sup[b]], w);
        }
      }
      nmod_poly_scalar_addmul_nmod(row + j, m,
                                   n_mulmod2_preinv(c, e[j] % p, p,
                                                    w->mod.ninv));
    }
  }
  free(sup);
  nmod_poly_clear(m);
}

/* the indices, in increasing order, of nv equations whose Jacobian is
 * invertible modulo (w, p) at the parametrization P, chosen from the one
 * with the fewest terms: an equation is kept if its row, reduced by those
 * of the equations kept, has an entry invertible modulo w in a column
 * without pivot, which becomes its pivot; NULL if there are not nv such
 * equations */
static long *hensel_select(const hensel_sys_t *hs, const hensel_prime_t *P) {
  const long nv = hs->nv, nall = hs->nall;
  const mp_limb_t p = P->p;
  long *ord = malloc(nall * sizeof(long));
  for (long i = 0; i < nall; i++) {
    long k = i;
    for (; k > 0 && hs->all[ord[k - 1]].nt > hs->all[i].nt; k--) {
      ord[k] = ord[k - 1];
    }
    ord[k] = i;
  }
  nmod_poly_struct **pw = malloc(nv * sizeof(nmod_poly_struct *));
  for (long v = 0; v < nv; v++) {
    pw[v] = malloc((hs->maxe[v] + 1) * sizeof(nmod_poly_struct));
    for (long k = 0; k <= hs->maxe[v]; k++) {
      nmod_poly_init(pw[v] + k, p);
      if (k == 0) {
        nmod_poly_one(pw[v]);
      } else {
        nmod_poly_mulmod(pw[v] + k, pw[v] + k - 1, P->X + v, P->w);
      }
    }
  }
  /* the rows kept, reduced, with a pivot 1 in the column piv[k] */
  nmod_poly_struct *R = malloc(nv * nv * sizeof(nmod_poly_struct));
  for (long i = 0; i < nv * nv; i++) {
    nmod_poly_init(R + i, p);
  }
  long *piv = malloc(nv * sizeof(long));
  char *used = calloc(nv, 1);
  long *sub = malloc(nv * sizeof(long));
  nmod_poly_t g, tmp;
  nmod_poly_init(g, p);
  nmod_poly_init(tmp, p);
  long nk = 0;
  for (long i = 0; i < nall && nk < nv; i++) {
    nmod_poly_struct *r = R + nk * nv;
    hensel_jac_row(r, hs->all + ord[i], pw, P->w, nv);
    for (long k = 0; k < nk; k++) {
      nmod_poly_set(g, r + piv[k]);
      if (nmod_poly_is_zero(g)) {
        continue;
      }
      for (long j = 0; j < nv; j++) {
        nmod_poly_mulmod(tmp, g, R + k * nv + j, P->w);
        nmod_poly_sub(r + j, r + j, tmp);
      }
    }
    long c = -1;
    for (long j = 0; j < nv && c < 0; j++) {
      if (used[j] || nmod_poly_is_zero(r + j)) {
        continue;
      }
      nmod_poly_gcd(g, r + j, P->w);
      if (nmod_poly_degree(g) == 0) {
        c = j;
      }
    }
    if (c < 0) {
      continue;
    }
    nmod_poly_invmod(g, r + c, P->w);
    for (long j = 0; j < nv; j++) {
      nmod_poly_mulmod(r + j, r + j, g, P->w);
    }
    used[c] = 1;
    piv[nk] = c;
    sub[nk++] = ord[i];
  }
  nmod_poly_clear(g);
  nmod_poly_clear(tmp);
  for (long i = 0; i < nv * nv; i++) {
    nmod_poly_clear(R + i);
  }
  free(R);
  for (long v = 0; v < nv; v++) {
    for (long k = 0; k <= hs->maxe[v]; k++) {
      nmod_poly_clear(pw[v] + k);
    }
    free(pw[v]);
  }
  free(pw);
  free(piv);
  free(used);
  free(ord);
  if (nk < nv) {
    free(sub);
    return NULL;
  }
  for (long i = 1; i < nv; i++) {
    for (long k = i; k > 0 && sub[k - 1] > sub[k]; k--) {
      const long t = sub[k];
      sub[k] = sub[k - 1];
      sub[k - 1] = t;
    }
  }
  return sub;
}

static void hensel_sys_clear(hensel_sys_t *hs) {
  hensel_sys_unset(hs);
  for (long i = 0; i < hs->nall; i++) {
    hensel_poly_clear(hs->all + i);
  }
  free(hs->all);
  free(hs->sys);
  free(hs->jac);
  free(hs->maxe);
}

static void hensel_data_init(hensel_data_t *hd) {
  hd->nprimes = 0;
  hd->alloc = 0;
  hd->pr = NULL;
}

/* adds the normalized parametrization modulo a prime (see
 * normalize_nmod_param): x_v = -coords[v] / denom mod elim */
static void hensel_add_param(hensel_data_t *hd, const param_t *param) {
  if (hd->nprimes == hd->alloc) {
    hd->alloc = MAX(8, 2 * hd->alloc);
    hd->pr = realloc(hd->pr, hd->alloc * sizeof(hensel_prime_t));
  }
  hensel_prime_t *P = hd->pr + hd->nprimes;
  memset(P, 0, sizeof(hensel_prime_t));
  const mp_limb_t p = param->charac;
  const long nv = param->nvars;
  P->p = p;
  P->ng = -1;
  nmod_poly_init(P->w, p);
  nmod_poly_make_monic(P->w, param->elim);
  P->X = malloc(nv * sizeof(nmod_poly_struct));
  nmod_poly_t inv;
  nmod_poly_init(inv, p);
  nmod_poly_invmod(inv, param->denom, P->w);
  nmod_poly_neg(inv, inv);
  for (long v = 0; v < nv - 1; v++) {
    nmod_poly_init(P->X + v, p);
    nmod_poly_mulmod(P->X + v, param->coords[v], inv, P->w);
  }
  nmod_poly_init(P->X + nv - 1, p);
  nmod_poly_set_coeff_ui(P->X + nv - 1, 1, 1);
  nmod_poly_clear(inv);
  hd->nprimes++;
}

static void hensel_data_clear(hensel_data_t *hd, const long nv) {
  for (long j = 0; j < hd->nprimes; j++) {
    hensel_prime_t *P = hd->pr + j;
    nmod_poly_clear(P->w);
    for (long v = 0; v < nv; v++) {
      nmod_poly_clear(P->X + v);
    }
    free(P->X);
    if (P->ng > 0) {
      for (long r = 0; r < P->ng; r++) {
        fmpz_poly_clear(P->g + r);
        for (long v = 0; v < nv - 1; v++) {
          fmpz_poly_clear(P->B + r * (nv - 1) + v);
        }
      }
      free(P->g);
      free(P->B);
      free(P->gdeg);
      fmpz_clear(P->q);
    }
  }
  free(hd->pr);
  hensel_data_init(hd);
}

/* ------------------------------------------------------------------ */
/* Newton's iteration                                                  */

static void hensel_set_finv(hensel_lift_t *L, const fmpz_mod_ctx_t ctx) {
  fmpz_mod_poly_t r;
  fmpz_mod_poly_init(r, ctx);
  fmpz_mod_poly_reverse(r, L->f, L->f->length, ctx);
  fmpz_mod_poly_inv_series(L->finv, r, L->f->length, ctx);
  fmpz_mod_poly_clear(r, ctx);
}

#define HENSEL_MULMOD(r, a, b, L, ctx) \
  fmpz_mod_poly_mulmod_preinv(r, a, b, (L)->f, (L)->finv, ctx)

/* res[i] = Q_i(X) mod f for the equations (w = 0) or the entries of the
 * Jacobian (w = 1), from the values of their monomials, computed by
 * increasing degree; the terms of each monomial are added to the res[i]
 * once its value is known, and the value is freed once its successors are
 * computed, so that about two degrees of them are held at once (sot1: 12328
 * monomials in J, at most mt->nlive[1] values). If mt->delay[w] is set,
 * the products without successor are not reduced modulo f: the res[i]
 * then have degree up to 2 deg f - 2 and are reduced at the end */
static void hensel_eval_set(fmpz_mod_poly_struct *res, const hensel_sys_t *hs,
                            const int w, const hensel_lift_t *L,
                            const fmpz_mod_ctx_t ctx) {
  const hensel_mon_t *mt = &hs->mt;
  const long nv = hs->nv, np = w == 0 ? nv : nv * nv;
  const hensel_poly_t *Q = w == 0 ? hs->sys : hs->jac;
  /* the terms of each monomial: equation and coefficient */
  long *start = calloc(mt->nmon + 1, sizeof(long));
  for (long i = 0; i < np; i++) {
    for (long t = 0; t < Q[i].nt; t++) {
      start[Q[i].mon[t] + 1]++;
    }
  }
  for (long m = 0; m < mt->nmon; m++) {
    start[m + 1] += start[m];
  }
  long *teq = malloc(start[mt->nmon] * sizeof(long));
  const fmpz **tcf = malloc(start[mt->nmon] * sizeof(fmpz *));
  long *pos = malloc(mt->nmon * sizeof(long));
  memcpy(pos, start, mt->nmon * sizeof(long));
  for (long i = 0; i < np; i++) {
    for (long t = 0; t < Q[i].nt; t++) {
      const long k = pos[Q[i].mon[t]]++;
      teq[k] = i;
      tcf[k] = Q[i].c + t;
    }
  }
  long *nsucc = pos;
  memset(nsucc, 0, mt->nmon * sizeof(long));
  for (long j = 0; j < mt->nord[w]; j++) {
    const long m = mt->ord[w][j];
    if (mt->pred[m] >= 0) {
      nsucc[mt->pred[m]]++;
    }
  }
  /* the terms are added with their integer coefficients, which are small,
   * and the sums are reduced once: reduced modulo q, a negative
   * coefficient would have the size of q */
  const slong dl = fmpz_mod_poly_degree(L->f, ctx);
  const int delay = mt->delay[w];
  const slong rl = delay ? 2 * dl - 1 : dl;
  fmpz *pr = delay ? _fmpz_vec_init(rl) : NULL;
  for (long i = 0; i < np; i++) {
    fmpz_mod_poly_fit_length(res + i, rl, ctx);
    _fmpz_vec_zero(res[i].coeffs, rl);
  }
  fmpz_mod_poly_struct *val = malloc(mt->nmon * sizeof(fmpz_mod_poly_struct));
  for (long j = 0; j < mt->nord[w]; j++) {
    const long m = mt->ord[w][j], pm = mt->pred[m];
    fmpz_mod_poly_init(val + m, ctx);
    if (pm < 0) {
      fmpz_mod_poly_one(val + m, ctx);
    } else if (mt->deg[m] == 1) {
      fmpz_mod_poly_set(val + m, L->X + mt->var[m], ctx);
    } else if (delay && nsucc[m] == 0) {
      const fmpz_mod_poly_struct *a = val + pm, *b = L->X + mt->var[m];
      if (a->length > 0 && b->length > 0) {
        const slong ln = a->length + b->length - 1;
        if (a->length >= b->length) {
          _fmpz_poly_mul(pr, a->coeffs, a->length, b->coeffs, b->length);
        } else {
          _fmpz_poly_mul(pr, b->coeffs, b->length, a->coeffs, a->length);
        }
        for (long k = start[m]; k < start[m + 1]; k++) {
          _fmpz_vec_scalar_addmul_fmpz(res[teq[k]].coeffs, pr, ln, tcf[k]);
        }
      }
      if (--nsucc[pm] == 0) {
        fmpz_mod_poly_clear(val + pm, ctx);
      }
      fmpz_mod_poly_clear(val + m, ctx);
      continue;
    } else {
      HENSEL_MULMOD(val + m, val + pm, L->X + mt->var[m], L, ctx);
    }
    if (pm >= 0 && --nsucc[pm] == 0) {
      fmpz_mod_poly_clear(val + pm, ctx);
    }
    for (long k = start[m]; k < start[m + 1]; k++) {
      _fmpz_vec_scalar_addmul_fmpz(res[teq[k]].coeffs, val[m].coeffs,
                                   val[m].length, tcf[k]);
    }
    if (nsucc[m] == 0) {
      fmpz_mod_poly_clear(val + m, ctx);
    }
  }
  fmpz_mod_poly_t qq, rr;
  fmpz_mod_poly_init(qq, ctx);
  fmpz_mod_poly_init(rr, ctx);
  for (long i = 0; i < np; i++) {
    _fmpz_vec_scalar_mod_fmpz(res[i].coeffs, res[i].coeffs, rl,
                              fmpz_mod_ctx_modulus(ctx));
    _fmpz_mod_poly_set_length(res + i, rl);
    _fmpz_mod_poly_normalise(res + i);
    if (res[i].length > dl) {
      fmpz_mod_poly_divrem_newton_n_preinv(qq, rr, res + i, L->f, L->finv,
                                           ctx);
      fmpz_mod_poly_swap(res + i, rr, ctx);
    }
  }
  fmpz_mod_poly_clear(qq, ctx);
  fmpz_mod_poly_clear(rr, ctx);
  if (delay) {
    _fmpz_vec_clear(pr, rl);
  }
  free(val);
  free(nsucc);
  free(tcf);
  free(teq);
  free(start);
}

/* C = A B for nv x nv matrices over (Z/q)[t]/(f): product over Z[t], then
 * reduction */
static void hensel_matmul(fmpz_mod_poly_struct *C,
                          const fmpz_mod_poly_struct *A,
                          const fmpz_mod_poly_struct *B, const long nv,
                          const hensel_lift_t *L, const fmpz_mod_ctx_t ctx) {
#if HENSEL_NTT
  const long D = fmpz_mod_poly_degree(L->f, ctx);
  if (L->ntt != NULL && D >= HENSEL_NTT_MIND) {
    fmpz_poly_struct *Az = calloc(3 * nv * nv, sizeof(fmpz_poly_struct));
    fmpz_poly_struct *Bz = Az + nv * nv, *Cz = Bz + nv * nv;
    long bB = 0;
    for (long i = 0; i < 3 * nv * nv; i++) {
      fmpz_poly_init(Az + i);
    }
    for (long i = 0; i < nv * nv; i++) {
      fmpz_mod_poly_get_fmpz_poly(Az + i, A + i, ctx);
      fmpz_mod_poly_get_fmpz_poly(Bz + i, B + i, ctx);
      bB = FLINT_MAX(bB, FLINT_ABS(_fmpz_vec_max_bits(Bz[i].coeffs,
                                                      Bz[i].length)));
    }
    hensel_ntt_mat_t M;
    const int ok = hensel_ntt_mat_set(&M, Az, nv, nv, bB, D, L->ntt);
    if (ok) {
      hensel_ntt_mul(Cz, &M, Bz, nv, D, L->ntt);
      hensel_ntt_mat_clear(&M);
      fmpz_mod_poly_t r;
      fmpz_mod_poly_init(r, ctx);
      for (long i = 0; i < nv * nv; i++) {
        fmpz_mod_poly_set_fmpz_poly(r, Cz + i, ctx);
        fmpz_mod_poly_rem(C + i, r, L->f, ctx);
      }
      fmpz_mod_poly_clear(r, ctx);
    }
    for (long i = 0; i < 3 * nv * nv; i++) {
      fmpz_poly_clear(Az + i);
    }
    free(Az);
    if (ok) {
      return;
    }
  }
#endif
  fmpz_poly_mat_t Az, Bz, Cz;
  fmpz_poly_mat_init(Az, nv, nv);
  fmpz_poly_mat_init(Bz, nv, nv);
  fmpz_poly_mat_init(Cz, nv, nv);
  for (long i = 0; i < nv; i++) {
    for (long j = 0; j < nv; j++) {
      fmpz_mod_poly_get_fmpz_poly(fmpz_poly_mat_entry(Az, i, j),
                                  A + i * nv + j, ctx);
      fmpz_mod_poly_get_fmpz_poly(fmpz_poly_mat_entry(Bz, i, j),
                                  B + i * nv + j, ctx);
    }
  }
  fmpz_poly_mat_mul(Cz, Az, Bz);
  fmpz_mod_poly_t r;
  fmpz_mod_poly_init(r, ctx);
  for (long i = 0; i < nv; i++) {
    for (long j = 0; j < nv; j++) {
      fmpz_mod_poly_set_fmpz_poly(r, fmpz_poly_mat_entry(Cz, i, j), ctx);
      fmpz_mod_poly_rem(C + i * nv + j, r, L->f, ctx);
    }
  }
  fmpz_mod_poly_clear(r, ctx);
  fmpz_poly_mat_clear(Az);
  fmpz_poly_mat_clear(Bz);
  fmpz_poly_mat_clear(Cz);
}

/* B = A^-1 modulo (f, p) by Gauss-Jordan elimination with invertible
 * pivots, made of random combinations of rows when no entry is; returns 0
 * if J is not invertible (or if four combinations failed, each with a
 * probability of about deg f / p) */
static int hensel_inv_mod_p(nmod_poly_struct *B, nmod_poly_struct *A,
                            const nmod_poly_t f, const long nv) {
  const mp_limb_t p = f->mod.n;
  int ok = 1;
  nmod_poly_t g, iv, tmp, fac;
  nmod_poly_init(g, p);
  nmod_poly_init(iv, p);
  nmod_poly_init(tmp, p);
  nmod_poly_init(fac, p);
  /* a fixed generator, so that the random numbers of msolve are not
   * consumed */
  uint64_t seed = 0x2545f4914f6cdd1dULL;
  for (long i = 0; i < nv; i++) {
    for (long j = 0; j < nv; j++) {
      if (i == j) {
        nmod_poly_one(B + i * nv + j);
      } else {
        nmod_poly_zero(B + i * nv + j);
      }
    }
  }
  for (long k = 0; k < nv && ok; k++) {
    long piv = -1;
    for (long i = k; i < nv && piv < 0; i++) {
      if (nmod_poly_is_zero(A + i * nv + k)) {
        continue;
      }
      nmod_poly_gcd(g, A + i * nv + k, f);
      if (nmod_poly_degree(g) == 0) {
        piv = i;
      }
    }
    /* no entry of the column is invertible modulo f, which happens when
     * f splits, each entry vanishing modulo a factor of f while the matrix
     * is invertible (jorge square1.sZ): the row k plus random multiples
     * of the rows below has an invertible entry, for random multiples,
     * if the rows k to nv - 1 have rank nv - k */
    for (int a = 0; a < 4 && piv < 0; a++) {
      for (long i = k + 1; i < nv; i++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        const mp_limb_t r = n_mod2_preinv(seed >> 1, p, f->mod.ninv);
        for (long j = 0; j < nv; j++) {
          nmod_poly_scalar_addmul_nmod(A + k * nv + j, A + i * nv + j, r);
          nmod_poly_scalar_addmul_nmod(B + k * nv + j, B + i * nv + j, r);
        }
      }
      nmod_poly_gcd(g, A + k * nv + k, f);
      if (!nmod_poly_is_zero(A + k * nv + k) && nmod_poly_degree(g) == 0) {
        piv = k;
      }
    }
    if (piv < 0) {
      ok = 0;
      break;
    }
    if (piv != k) {
      for (long j = 0; j < nv; j++) {
        nmod_poly_swap(A + k * nv + j, A + piv * nv + j);
        nmod_poly_swap(B + k * nv + j, B + piv * nv + j);
      }
    }
    nmod_poly_invmod(iv, A + k * nv + k, f);
    for (long j = 0; j < nv; j++) {
      nmod_poly_mulmod(A + k * nv + j, A + k * nv + j, iv, f);
      nmod_poly_mulmod(B + k * nv + j, B + k * nv + j, iv, f);
    }
    for (long i = 0; i < nv; i++) {
      if (i == k || nmod_poly_is_zero(A + i * nv + k)) {
        continue;
      }
      nmod_poly_set(fac, A + i * nv + k);
      for (long j = 0; j < nv; j++) {
        nmod_poly_mulmod(tmp, fac, A + k * nv + j, f);
        nmod_poly_sub(A + i * nv + j, A + i * nv + j, tmp);
        nmod_poly_mulmod(tmp, fac, B + k * nv + j, f);
        nmod_poly_sub(B + i * nv + j, B + i * nv + j, tmp);
      }
    }
  }
  nmod_poly_clear(g);
  nmod_poly_clear(iv);
  nmod_poly_clear(tmp);
  nmod_poly_clear(fac);
  return ok;
}

/* r = A v modulo (f, q) for an nv x nv matrix A and a vector v of
 * polynomials of degree less than deg f: the sums are formed over Z[t] and
 * reduced once; r has coefficients in [0, q) */
/* r = acc mod (q, f), with coefficients in [0, q) */
static void hensel_matvec_red(fmpz_poly_t r, const fmpz_poly_t acc,
                              fmpz_mod_poly_t t, fmpz_mod_poly_t qq,
                              fmpz_mod_poly_t rr, const hensel_lift_t *L,
                              const fmpz_mod_ctx_t ctx) {
  fmpz_mod_poly_set_fmpz_poly(t, acc, ctx);
  if (t->length >= L->f->length) {
    fmpz_mod_poly_divrem_newton_n_preinv(qq, rr, t, L->f, L->finv, ctx);
    fmpz_mod_poly_get_fmpz_poly(r, rr, ctx);
  } else {
    fmpz_mod_poly_get_fmpz_poly(r, t, ctx);
  }
}

static void hensel_matvec(fmpz_poly_struct *r, const fmpz_poly_struct *A,
                          const fmpz_poly_struct *v, const long nv,
                          const hensel_lift_t *L, const fmpz_mod_ctx_t ctx) {
  fmpz_poly_t acc, pr;
  fmpz_poly_init(acc);
  fmpz_poly_init(pr);
  fmpz_mod_poly_t t, qq, rr;
  fmpz_mod_poly_init(t, ctx);
  fmpz_mod_poly_init(qq, ctx);
  fmpz_mod_poly_init(rr, ctx);
  for (long i = 0; i < nv; i++) {
    fmpz_poly_zero(acc);
    for (long j = 0; j < nv; j++) {
      if (!fmpz_poly_is_zero(A + i * nv + j) && !fmpz_poly_is_zero(v + j)) {
        fmpz_poly_mul(pr, A + i * nv + j, v + j);
        fmpz_poly_add(acc, acc, pr);
      }
    }
    hensel_matvec_red(r + i, acc, t, qq, rr, L, ctx);
  }
  fmpz_poly_clear(acc);
  fmpz_poly_clear(pr);
  fmpz_mod_poly_clear(t, ctx);
  fmpz_mod_poly_clear(qq, ctx);
  fmpz_mod_poly_clear(rr, ctx);
}

#if HENSEL_NTT
/* the same, A given by its transforms M */
static void hensel_matvec_ntt(fmpz_poly_struct *r, const hensel_ntt_mat_t *M,
                              const fmpz_poly_struct *v, const long nv,
                              hensel_ntt_t *T, const hensel_lift_t *L,
                              const fmpz_mod_ctx_t ctx) {
  fmpz_poly_struct *acc = malloc(nv * sizeof(fmpz_poly_struct));
  for (long i = 0; i < nv; i++) {
    fmpz_poly_init(acc + i);
  }
  hensel_ntt_mul(acc, M, v, 1, fmpz_mod_poly_degree(L->f, ctx), T);
  fmpz_mod_poly_t t, qq, rr;
  fmpz_mod_poly_init(t, ctx);
  fmpz_mod_poly_init(qq, ctx);
  fmpz_mod_poly_init(rr, ctx);
  for (long i = 0; i < nv; i++) {
    hensel_matvec_red(r + i, acc + i, t, qq, rr, L, ctx);
    fmpz_poly_clear(acc + i);
  }
  free(acc);
  fmpz_mod_poly_clear(t, ctx);
  fmpz_mod_poly_clear(qq, ctx);
  fmpz_mod_poly_clear(rr, ctx);
}
#endif

/* dX = J(X)^-1 F(X) modulo p^n, where F(X) = 0 modulo p^o, o < n <= 2o,
 * so that dX = p^o y with J y = F(X) / p^o modulo p^(n - o); J^-1 is
 * known modulo p^jprec: if this is not enough, y is refined by
 * corrections y += p^e J^-1 s, s = (s - J y_e) / p^jprec, each one at a
 * precision lower by jprec, with J(X) evaluated modulo p^(n - o) */
static void hensel_solve(fmpz_mod_poly_struct *dX,
                         const fmpz_mod_poly_struct *Fv,
                         hensel_lift_t *L, const hensel_sys_t *hs,
                         const long o, const long n,
                         const fmpz_mod_ctx_t ctx) {
  const long nv = hs->nv, h = n - o, jp = L->jprec;
  double rt = realtime();
  fmpz_t pr, pe, qm;
  fmpz_init(pr);
  fmpz_init(pe);
  fmpz_init(qm);
  /* p, from the modulus p^n */
  fmpz_root(pr, fmpz_mod_ctx_modulus(ctx), n);
  fmpz_poly_t fz, tz;
  fmpz_poly_init(fz);
  fmpz_poly_init(tz);
  fmpz_mod_poly_get_fmpz_poly(fz, L->f, ctx);
  fmpz_poly_struct *Jinvz = malloc(nv * nv * sizeof(fmpz_poly_struct));
  for (long i = 0; i < nv * nv; i++) {
    fmpz_poly_init(Jinvz + i);
    fmpz_mod_poly_get_fmpz_poly(Jinvz + i, L->Jinv + i, ctx);
  }
  /* s = F(X) / p^o */
  fmpz_pow_ui(pe, pr, o);
  fmpz_poly_struct *s = malloc(nv * sizeof(fmpz_poly_struct));
  for (long i = 0; i < nv; i++) {
    fmpz_poly_init(s + i);
    fmpz_mod_poly_get_fmpz_poly(s + i, Fv + i, ctx);
    fmpz_poly_scalar_fdiv_fmpz(s + i, s + i, pe);
    fmpz_mod_poly_zero(dX + i, ctx);
  }
  /* J(X) modulo p^h, only if corrections are needed */
  fmpz_poly_struct *Jz = NULL;
  if (jp < h) {
    fmpz_mod_ctx_t ch;
    fmpz_pow_ui(qm, pr, h);
    fmpz_mod_ctx_init(ch, qm);
    hensel_lift_t Lh[1];
    fmpz_mod_poly_init(Lh->f, ch);
    fmpz_mod_poly_init(Lh->finv, ch);
    fmpz_mod_poly_set_fmpz_poly(Lh->f, fz, ch);
    hensel_set_finv(Lh, ch);
    Lh->X = malloc(nv * sizeof(fmpz_mod_poly_struct));
    for (long k = 0; k < nv; k++) {
      fmpz_mod_poly_init(Lh->X + k, ch);
      fmpz_mod_poly_get_fmpz_poly(tz, L->X + k, ctx);
      fmpz_mod_poly_set_fmpz_poly(Lh->X + k, tz, ch);
    }
    fmpz_mod_poly_struct *Je = malloc(nv * nv * sizeof(fmpz_mod_poly_struct));
    for (long i = 0; i < nv * nv; i++) {
      fmpz_mod_poly_init(Je + i, ch);
    }
    hensel_eval_set(Je, hs, 1, Lh, ch);
    Jz = malloc(nv * nv * sizeof(fmpz_poly_struct));
    for (long i = 0; i < nv * nv; i++) {
      fmpz_poly_init(Jz + i);
      fmpz_mod_poly_get_fmpz_poly(Jz + i, Je + i, ch);
      fmpz_mod_poly_clear(Je + i, ch);
    }
    free(Je);
    for (long k = 0; k < nv; k++) {
      fmpz_mod_poly_clear(Lh->X + k, ch);
    }
    free(Lh->X);
    fmpz_mod_poly_clear(Lh->f, ch);
    fmpz_mod_poly_clear(Lh->finv, ch);
    fmpz_mod_ctx_clear(ch);
  }
  L->tj += realtime() - rt;
  rt = realtime();
#if HENSEL_NTT
  /* the transforms of J^-1, formed at the first step with corrections and
   * kept until it is updated, and those of J for the corrections: y_e and
   * s have at most jprec digits */
  const hensel_ntt_mat_t *Jh = NULL;
  hensel_ntt_mat_t Jzh;
  int useJ = 0;
  const long D = fz->length - 1;
  if (L->ntt != NULL && jp < h && D >= HENSEL_NTT_MIND) {
    fmpz_pow_ui(pe, pr, jp);
    const long bv = fmpz_bits(pe);
    if (L->Jh.d == NULL) {
      hensel_ntt_mat_set(&L->Jh, Jinvz, nv, nv, bv, D, L->ntt);
    }
    if (L->Jh.d != NULL) {
      Jh = &L->Jh;
    }
    useJ = hensel_ntt_mat_set(&Jzh, Jz, nv, nv, bv, D, L->ntt);
  }
#endif
  fmpz_mod_poly_t b;
  fmpz_mod_poly_init(b, ctx);
  fmpz_poly_struct *yz = malloc(nv * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *tv = malloc(nv * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *Am = malloc(nv * nv * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *Bm = malloc(nv * nv * sizeof(fmpz_poly_struct));
  for (long i = 0; i < nv * nv; i++) {
    fmpz_poly_init(Am + i);
    fmpz_poly_init(Bm + i);
  }
  for (long i = 0; i < nv; i++) {
    fmpz_poly_init(yz + i);
    fmpz_poly_init(tv + i);
  }
  /* y_e = J^-1 s modulo p^a, a = min(jprec, h - e), is added to dX with
   * the factor p^(o + e), as in Dixon's lifting: J^-1 is only known
   * modulo p^jprec, so that s and y_e are truncated to a digits, and only
   * J y_e is formed at the precision h - e of s */
  for (long e = 0; e < h; e += jp) {
    const long m = h - e, a = MIN(jp, m);
    fmpz_mod_ctx_t ca;
    fmpz_pow_ui(qm, pr, a);
    fmpz_mod_ctx_init(ca, qm);
    hensel_lift_t La[1];
    fmpz_mod_poly_init(La->f, ca);
    fmpz_mod_poly_init(La->finv, ca);
    fmpz_mod_poly_set_fmpz_poly(La->f, fz, ca);
    hensel_set_finv(La, ca);
    for (long i = 0; i < nv; i++) {
      fmpz_poly_scalar_mod_fmpz(tv + i, s + i, qm);
    }
#if HENSEL_NTT
    if (Jh != NULL) {
      /* J^-1 is not reduced modulo p^a: the result is */
      hensel_matvec_ntt(yz, Jh, tv, nv, L->ntt, La, ca);
    } else
#endif
    {
      for (long i = 0; i < nv * nv; i++) {
        if (jp > a) {
          fmpz_poly_scalar_mod_fmpz(Am + i, Jinvz + i, qm);
        } else {
          fmpz_poly_set(Am + i, Jinvz + i);
        }
      }
      hensel_matvec(yz, Am, tv, nv, La, ca);
    }
    fmpz_mod_poly_clear(La->f, ca);
    fmpz_mod_poly_clear(La->finv, ca);
    fmpz_mod_ctx_clear(ca);
    fmpz_pow_ui(pe, pr, o + e);
    for (long i = 0; i < nv; i++) {
      fmpz_poly_scalar_mul_fmpz(tz, yz + i, pe);
      fmpz_mod_poly_set_fmpz_poly(b, tz, ctx);
      fmpz_mod_poly_add(dX + i, dX + i, b, ctx);
    }
    if (e + jp < h) {
      /* s = (s - J y_e) / p^jprec modulo p^(m - jprec) */
      fmpz_mod_ctx_t cm;
      fmpz_pow_ui(qm, pr, m);
      fmpz_mod_ctx_init(cm, qm);
      hensel_lift_t Lm[1];
      fmpz_mod_poly_init(Lm->f, cm);
      fmpz_mod_poly_init(Lm->finv, cm);
      fmpz_mod_poly_set_fmpz_poly(Lm->f, fz, cm);
      hensel_set_finv(Lm, cm);
#if HENSEL_NTT
      if (useJ) {
        hensel_matvec_ntt(tv, &Jzh, yz, nv, L->ntt, Lm, cm);
      } else
#endif
      {
        for (long i = 0; i < nv * nv; i++) {
          fmpz_poly_scalar_mod_fmpz(Bm + i, Jz + i, qm);
        }
        hensel_matvec(tv, Bm, yz, nv, Lm, cm);
      }
      fmpz_pow_ui(pe, pr, jp);
      for (long i = 0; i < nv; i++) {
        fmpz_poly_sub(s + i, s + i, tv + i);
        fmpz_poly_scalar_mod_fmpz(s + i, s + i, qm);
        fmpz_poly_scalar_fdiv_fmpz(s + i, s + i, pe);
      }
      fmpz_mod_poly_clear(Lm->f, cm);
      fmpz_mod_poly_clear(Lm->finv, cm);
      fmpz_mod_ctx_clear(cm);
    }
  }
  for (long i = 0; i < nv; i++) {
    fmpz_poly_clear(yz + i);
    fmpz_poly_clear(tv + i);
  }
  free(yz);
  free(tv);
  for (long i = 0; i < nv * nv; i++) {
    fmpz_poly_clear(Am + i);
    fmpz_poly_clear(Bm + i);
  }
  free(Am);
  free(Bm);
  fmpz_mod_poly_clear(b, ctx);
  for (long i = 0; i < nv; i++) {
    fmpz_poly_clear(s + i);
  }
  free(s);
  for (long i = 0; i < nv * nv; i++) {
    fmpz_poly_clear(Jinvz + i);
    if (Jz != NULL) {
      fmpz_poly_clear(Jz + i);
    }
  }
  free(Jinvz);
  free(Jz);
#if HENSEL_NTT
  if (useJ) {
    hensel_ntt_mat_clear(&Jzh);
  }
#endif
  fmpz_poly_clear(fz);
  fmpz_poly_clear(tz);
  fmpz_clear(pe);
  fmpz_clear(qm);
  fmpz_clear(pr);
  L->tmv += realtime() - rt;
}

/* one Newton step from the precision p^o to p^n, o < n <= 2o, the context
 * has the modulus p^n; J^-1 is updated when upd is set */
static void hensel_step(hensel_lift_t *L, const hensel_sys_t *hs,
                        const int upd, const long o, const long n,
                        const fmpz_mod_ctx_t ctx) {
  const long nv = hs->nv;
  fmpz_mod_poly_t tmp, Delta, fd;
  fmpz_mod_poly_init(tmp, ctx);
  fmpz_mod_poly_init(Delta, ctx);
  fmpz_mod_poly_init(fd, ctx);
  hensel_set_finv(L, ctx);
  fmpz_mod_poly_struct *Fv = malloc(nv * sizeof(fmpz_mod_poly_struct));
  fmpz_mod_poly_struct *dX = malloc(nv * sizeof(fmpz_mod_poly_struct));
  for (long i = 0; i < nv; i++) {
    fmpz_mod_poly_init(Fv + i, ctx);
    fmpz_mod_poly_init(dX + i, ctx);
  }
  double rt = realtime();
  hensel_eval_set(Fv, hs, 0, L, ctx);
  L->tf += realtime() - rt;
  hensel_solve(dX, Fv, L, hs, o, n, ctx);
  rt = realtime();
  /* Delta = t - X_t + dX_t, where t - X_t vanishes only modulo the
   * previous precision */
  fmpz_mod_poly_zero(Delta, ctx);
  fmpz_mod_poly_set_coeff_ui(Delta, 1, 1, ctx);
  fmpz_mod_poly_sub(Delta, Delta, L->X + nv - 1, ctx);
  fmpz_mod_poly_add(Delta, Delta, dX + nv - 1, ctx);
  for (long i = 0; i < nv; i++) {
    fmpz_mod_poly_derivative(fd, L->X + i, ctx);
    HENSEL_MULMOD(tmp, fd, Delta, L, ctx);
    fmpz_mod_poly_sub(L->X + i, L->X + i, dX + i, ctx);
    fmpz_mod_poly_add(L->X + i, L->X + i, tmp, ctx);
  }
  fmpz_mod_poly_derivative(fd, L->f, ctx);
  HENSEL_MULMOD(tmp, fd, Delta, L, ctx);
  fmpz_mod_poly_add(L->f, L->f, tmp, ctx);
  hensel_set_finv(L, ctx);
  L->tup += realtime() - rt;
  if (upd) {
    rt = realtime();
    fmpz_mod_poly_struct *J = malloc(nv * nv * sizeof(fmpz_mod_poly_struct));
    fmpz_mod_poly_struct *T = malloc(nv * nv * sizeof(fmpz_mod_poly_struct));
    for (long i = 0; i < nv * nv; i++) {
      fmpz_mod_poly_init(J + i, ctx);
      fmpz_mod_poly_init(T + i, ctx);
    }
    hensel_eval_set(J, hs, 1, L, ctx);
    L->tj += realtime() - rt;
    rt = realtime();
    hensel_matmul(T, J, L->Jinv, nv, L, ctx);
    fmpz_t c;
    fmpz_init(c);
    for (long i = 0; i < nv * nv; i++) {
      fmpz_mod_poly_neg(T + i, T + i, ctx);
    }
    for (long i = 0; i < nv; i++) {
      fmpz_mod_poly_get_coeff_fmpz(c, T + i * nv + i, 0, ctx);
      fmpz_mod_add_ui(c, c, 2, ctx);
      fmpz_mod_poly_set_coeff_fmpz(T + i * nv + i, 0, c, ctx);
    }
    fmpz_clear(c);
    hensel_matmul(J, L->Jinv, T, nv, L, ctx);
    L->tmm += realtime() - rt;
    L->jprec = MIN(2 * L->jprec, n);
#if HENSEL_NTT
    hensel_ntt_mat_clear(&L->Jh);
#endif
    for (long i = 0; i < nv * nv; i++) {
      fmpz_mod_poly_swap(L->Jinv + i, J + i, ctx);
      fmpz_mod_poly_clear(J + i, ctx);
      fmpz_mod_poly_clear(T + i, ctx);
    }
    free(J);
    free(T);
  }
  for (long i = 0; i < nv; i++) {
    fmpz_mod_poly_clear(Fv + i, ctx);
    fmpz_mod_poly_clear(dX + i, ctx);
  }
  free(Fv);
  free(dX);
  fmpz_mod_poly_clear(tmp, ctx);
  fmpz_mod_poly_clear(Delta, ctx);
  fmpz_mod_poly_clear(fd, ctx);
}

/* ------------------------------------------------------------------ */
/* reconstruction                                                      */

/* residues of the coefficients of g and of the B_v (n = dg nv of them)
 * modulo M Q, Q the product of the moduli q_j of the T primes pr, formed
 * coefficient by coefficient from their residues res modulo M and from
 * the images at the primes, so that no vector modulo M Q is stored */
typedef struct {
  const ulong *res;       /* modulo M, packed, L limbs each; unused if M = 1 */
  slong L;
  const hensel_prime_t *pr;
  long T, dg, nv;
  fmpz *Mj;               /* M q_0 .. q_(j-1) */
  fmpz *inv;              /* its inverse modulo q_j */
  fmpz_t MQ;
} hensel_crt_t;

/* M is not kept: the products Mj start from it */
static void hensel_crt_init(hensel_crt_t *C, const fmpz_t M,
                            const hensel_prime_t *pr, const long T,
                            const long dg, const long nv) {
  C->res = NULL;
  C->L = 0;
  C->pr = pr;
  C->T = T;
  C->dg = dg;
  C->nv = nv;
  C->Mj = _fmpz_vec_init(T);
  C->inv = _fmpz_vec_init(T);
  fmpz_init_set(C->MQ, M);
  for (long j = 0; j < T; j++) {
    fmpz_set(C->Mj + j, C->MQ);
    if (!fmpz_is_one(C->MQ)) {
      fmpz_invmod(C->inv + j, C->MQ, pr[j].q);
    }
    fmpz_mul(C->MQ, C->MQ, pr[j].q);
  }
}

static void hensel_crt_clear(hensel_crt_t *C) {
  _fmpz_vec_clear(C->Mj, C->T);
  _fmpz_vec_clear(C->inv, C->T);
  fmpz_clear(C->MQ);
}

/* the image of the coefficient i at the prime P */
static const fmpz *hensel_crt_image(const hensel_prime_t *P, const long i,
                                    const long dg) {
  const fmpz_poly_struct *f = i < dg ? P->g : P->B + (i - dg) / dg;
  const long k = i < dg ? i : (i - dg) % dg;
  return k < f->length ? f->coeffs + k : NULL;
}

/* x = the residue of the coefficient i modulo M Q, from x = its residue
 * modulo M (ignored if M = 1); t is a temporary */
static void hensel_crt_update(fmpz_t x, const hensel_crt_t *C, const long i,
                              fmpz_t t) {
  for (long j = 0; j < C->T; j++) {
    const fmpz *b = hensel_crt_image(C->pr + j, i, C->dg);
    if (fmpz_is_one(C->Mj + j)) {
      if (b != NULL) {
        fmpz_set(x, b);
      } else {
        fmpz_zero(x);
      }
      continue;
    }
    /* x += M_j ((b - x) / M_j mod q_j) */
    if (b != NULL) {
      fmpz_sub(t, b, x);
    } else {
      fmpz_neg(t, x);
    }
    fmpz_mul(t, t, C->inv + j);
    fmpz_mod(t, t, C->pr[j].q);
    fmpz_addmul(x, t, C->Mj + j);
  }
}

static void hensel_crt_get(fmpz_t x, const hensel_crt_t *C, const long i,
                           fmpz_t t) {
  if (C->T == 0 || !fmpz_is_one(C->Mj)) {
    fmpz_set_ui_array(x, C->res + i * C->L, C->L);
  }
  hensel_crt_update(x, C, i, t);
}

/* evaluates the system at x_v = B_v / g' mod g (x_t = t) modulo p;
 * returns 1 if it vanishes, 0 if not, -1 if p divides a denominator */
static int hensel_check(const fmpq *rec, const long dg, const hensel_sys_t *hs,
                        const mp_limb_t p) {
  const long nv = hs->nv;
  nmod_t mod;
  nmod_init(&mod, p);
  nmod_poly_t g, gd, ig, tmp, acc, val;
  nmod_poly_init(g, p);
  nmod_poly_init(gd, p);
  nmod_poly_init(ig, p);
  nmod_poly_init(tmp, p);
  nmod_poly_init(acc, p);
  nmod_poly_init(val, p);
  nmod_poly_struct *X = malloc(nv * sizeof(nmod_poly_struct));
  for (long v = 0; v < nv; v++) {
    nmod_poly_init(X + v, p);
  }
  int ok = 1;
  /* g, then the B_v */
  for (long i = 0; i < dg * nv && ok == 1; i++) {
    const mp_limb_t den = fmpz_fdiv_ui(fmpq_denref(rec + i), p);
    if (den == 0) {
      ok = -1;
      break;
    }
    const mp_limb_t u = nmod_mul(fmpz_fdiv_ui(fmpq_numref(rec + i), p),
                                 nmod_inv(den, mod), mod);
    if (i < dg) {
      nmod_poly_set_coeff_ui(g, i, u);
    } else {
      nmod_poly_set_coeff_ui(X + (i - dg) / dg, (i - dg) % dg, u);
    }
  }
  if (ok == 1) {
    nmod_poly_set_coeff_ui(g, dg, 1);
    nmod_poly_derivative(gd, g);
    nmod_poly_gcdinv(tmp, ig, gd, g);
    if (nmod_poly_degree(tmp) != 0) {
      ok = -1;
    }
  }
  if (ok == 1) {
    for (long v = 0; v < nv - 1; v++) {
      nmod_poly_mulmod(X + v, X + v, ig, g);
    }
    nmod_poly_zero(X + nv - 1);
    nmod_poly_set_coeff_ui(X + nv - 1, 1, 1);
    nmod_poly_rem(X + nv - 1, X + nv - 1, g);
    for (long i = 0; i < hs->nall && ok == 1; i++) {
      const hensel_poly_t *P = hs->all + i;
      nmod_poly_zero(val);
      for (long m = 0; m < P->nt; m++) {
        nmod_poly_one(acc);
        for (long k = 0; k < nv; k++) {
          for (long e = 0; e < P->e[m * nv + k]; e++) {
            nmod_poly_mulmod(acc, acc, X + k, g);
          }
        }
        nmod_poly_scalar_mul_nmod(acc, acc, fmpz_fdiv_ui(P->c + m, p));
        nmod_poly_add(val, val, acc);
      }
      if (!nmod_poly_is_zero(val)) {
        ok = 0;
      }
    }
  }
  for (long v = 0; v < nv; v++) {
    nmod_poly_clear(X + v);
  }
  free(X);
  nmod_poly_clear(g);
  nmod_poly_clear(gd);
  nmod_poly_clear(ig);
  nmod_poly_clear(tmp);
  nmod_poly_clear(acc);
  nmod_poly_clear(val);
  return ok;
}

/* reconstruction of a = c d modulo M with the bounds of
 * hensel_reconstruct: M^(1/64), then M^(1/16), then balanced */
static int hensel_rr(fmpq_t r, const fmpz_t a, const fmpz_t M,
                     const slong *eb, fmpz_t Nb, fmpz_t Db) {
  for (int t = 0; t < 3; t++) {
    if (eb[t] < 0) {
      return fmpq_reconstruct_fmpz(r, a, M);
    }
    /* 2 N D < M / 2^HENSEL_RRMARGIN: a random residue is accepted with
     * probability about 2^-HENSEL_RRMARGIN */
    fmpz_one(Db);
    fmpz_mul_2exp(Db, Db, eb[t]);
    fmpz_fdiv_q(Nb, M, Db);
    fmpz_fdiv_q_2exp(Nb, Nb, 2 + HENSEL_RRMARGIN);
    if (fmpq_reconstruct_fmpz_2(r, a, M, Nb, Db)) {
      return 1;
    }
  }
  return 0;
}

/* reconstructs the block blk (g, then the B_v) from the highest degree,
 * with the running denominator d; returns the number of coefficients
 * reconstructed before the first failure */
static long hensel_rec_block(fmpq *rec, const hensel_crt_t *C, const long blk,
                             fmpz_t d, const slong *eb, fmpz *tmp) {
  const long dg = C->dg;
  const fmpz *M = C->MQ;
  fmpz *a = tmp, *c = tmp + 1, *t = tmp + 2, *Nb = tmp + 3, *Db = tmp + 4;
  fmpq_t r;
  fmpq_init(r);
  long k = 0;
  for (; k < dg; k++) {
    const long i = blk * dg + dg - 1 - k;
    hensel_crt_get(c, C, i, t);
    fmpz_mul(a, d, c);
    fmpz_mod(a, a, M);
    if (!hensel_rr(r, a, M, eb, Nb, Db)) {
      break;
    }
    fmpz_mul(fmpq_denref(rec + i), fmpq_denref(r), d);
    fmpz_set(fmpq_numref(rec + i), fmpq_numref(r));
    fmpq_canonicalise(rec + i);
    fmpz_mul(d, d, fmpq_denref(r));
  }
  fmpq_clear(r);
  return k;
}

/* rational reconstruction of the n = dg nv coefficients of g and of the
 * B_v modulo M = C->MQ, from the highest degrees down, with a running
 * common denominator d: d c is reconstructed with a denominator bound
 * M^(1/64), then M^(1/16), then a balanced one. d starts from d0 for g (1
 * if d0 is NULL). Each B_v starts from the denominator of g, then from 1
 * if it fails: the B_v often share the denominator of g (random2_9_9:
 * 39000 bits, which a start from 1 has to reconstruct again with each
 * B_v, 81000 bits instead of 52000), but not always (katsura 13: 27800
 * bits from the denominator of g instead of 20500). nrec receives the
 * number of coefficients of g reconstructed before the first failure.
 * Returns 1 if all succeed */
static int hensel_reconstruct(fmpq *rec, const hensel_crt_t *C,
                              const fmpz *d0, long *nrec) {
  const slong b = fmpz_bits(C->MQ);
  const slong eb[3] = {b / 64, b / 16, -1};
  fmpz *tmp = _fmpz_vec_init(5);
  fmpz_t d, dG;
  fmpz_init(d);
  fmpz_init(dG);
  if (d0 != NULL) {
    fmpz_set(d, d0);
  } else {
    fmpz_one(d);
  }
  const long k = hensel_rec_block(rec, C, 0, d, eb, tmp);
  if (nrec != NULL) {
    *nrec = k;
  }
  int ok = k == C->dg;
  fmpz_set(dG, d);
  for (long blk = 1; blk < C->nv && ok; blk++) {
    fmpz_set(d, dG);
    ok = hensel_rec_block(rec, C, blk, d, eb, tmp) == C->dg;
    if (!ok && !fmpz_is_one(dG)) {
      fmpz_one(d);
      ok = hensel_rec_block(rec, C, blk, d, eb, tmp) == C->dg;
    }
  }
  fmpz_clear(d);
  fmpz_clear(dG);
  _fmpz_vec_clear(tmp, 5);
  return ok;
}


/* ------------------------------------------------------------------ */
/* cost model                                                          */

/* The costs are counted in bit operations of the products of two
   polynomials of degree less than D with b1- and b2-bit coefficients,
   done by Kronecker substitution and an FFT product of L bits:
   L log2(L) with L = D (b1 + b2 + log2 D). A product modulo w with its
   preinverse costs three of them. */
static double hensel_mul_cost(const long D, const double b1, const double b2) {
  const double L = D * (b1 + b2 + log2((double)D));
  return L * log2(L);
}

static double hensel_mm_cost(const long D, const double b) {
  return 3 * hensel_mul_cost(D, b, b);
}

/* number of products of the evaluation of F (w = 0) or of J (w = 1) */
static long hensel_nprod(const hensel_sys_t *hs, const int w) {
  long c = 0;
  for (long j = 0; j < hs->mt.nord[w]; j++) {
    c += hs->mt.deg[hs->mt.ord[w][j]] >= 2;
  }
  return c;
}

/* number of terms of F (w = 0) or of J (w = 1) */
static long hensel_nterms(const hensel_sys_t *hs, const int w) {
  const long np = w == 0 ? hs->nv : hs->nv * hs->nv;
  const hensel_poly_t *Q = w == 0 ? hs->sys : hs->jac;
  long c = 0;
  for (long i = 0; i < np; i++) {
    c += Q[i].nt;
  }
  return c;
}

/* cost of the terms of an evaluation of F or J (see hensel_eval_set): D
 * coefficients of b bits multiplied by a small integer and added, for
 * each term, about 5.5 ns + 0.008 ns per bit for a coefficient (measured
 * for D = 32 to 2048 and 60 to 3900 bits, in units of 5e-11 s, see
 * hensel_choose); this matters on dense systems (cp_d_3_n_7_p_2: 24987
 * terms for 330 monomials in F, 63498 terms in J) */
#define HENSEL_TERM_C0 110.0
#define HENSEL_TERM_C1 0.16
static double hensel_term_cost(const hensel_sys_t *hs, const int w,
                               const long D, const double b) {
  return (double)hensel_nterms(hs, w) * D
    * (HENSEL_TERM_C0 + HENSEL_TERM_C1 * b);
}

/* cost of an evaluation of F (w = 0) or of J (w = 1) at b bits (see
 * hensel_eval_set): a product modulo w is a product and a reduction,
 * which costs two products; if the products without successor are not
 * reduced, their terms have 2 D coefficients of about 2 b bits, and nred
 * polynomials are reduced */
static double hensel_eval_cost(const hensel_sys_t *hs, const int w,
                               const long D, const double b) {
  const hensel_mon_t *mt = &hs->mt;
  double c = hensel_nprod(hs, w) * hensel_mm_cost(D, b)
    + hensel_term_cost(hs, w, D, b);
  if (mt->delay[w]) {
    const double red = hensel_mm_cost(D, b) - hensel_mul_cost(D, b, b);
    c += (mt->nred[w] - mt->nleaf[w]) * red
      + (double)mt->tleaf[w] * D * (HENSEL_TERM_C0 + 3 * HENSEL_TERM_C1 * b);
  }
  return c;
}

/* model costs of the products by transforms (see hensel_ntt_mul),
 * measured for nv = 12, D = 2046 and 31 to 1000 bits (2 to 42 primes):
 * per coefficient and prime, its residue (about 9 ns); per element and
 * level, a transform (0.15 ns); per element, a product added (0.86 ns);
 * per coefficient, the CRT (11.7 ns np + 0.465 ns np^2, within 10% of the
 * measures: 25 ns for 2 primes, 1.31 us for 42) */
#define HENSEL_NTT_CCONV 180.0
#define HENSEL_NTT_CFFT 3.0
#define HENSEL_NTT_CPT 17.0
#define HENSEL_NTT_CCRT1 234.0
#define HENSEL_NTT_CCRT2 9.3

/* the number of primes of the transforms for the products of an nv x nv
 * matrix with coefficients of bA bits by one with coefficients of bB bits,
 * 0 if they are not used */
static long hensel_ntt_np(const hensel_sys_t *hs, const long D,
                          const double bA, const double bB) {
#if HENSEL_NTT
  const long np = (long)ceil((bA + bB + log2((double)hs->nv * D) + 2) / 49.9);
  return hs->ntt && D >= HENSEL_NTT_MIND && np <= HENSEL_NTT_MAXP ? np : 0;
#else
  (void)hs;
  (void)D;
  (void)bA;
  (void)bB;
  return 0;
#endif
}

/* transform, inverse transform (with the CRT) and products added of the
 * transforms of one entry, with np primes */
static double hensel_ntt_cost(const long D, const long np, const int what) {
  const double L = MAX(ceil(log2(2.0 * D - 1)), 8), N = pow(2, L);
  if (what == 0) {
    return np * (D * HENSEL_NTT_CCONV + N * L * HENSEL_NTT_CFFT);
  }
  if (what == 1) {
    return np * N * L * HENSEL_NTT_CFFT
      + (2.0 * D - 1) * np * (HENSEL_NTT_CCRT1 + HENSEL_NTT_CCRT2 * np);
  }
  return np * 2.0 * D * HENSEL_NTT_CPT;
}

/* the products of an nv x nv matrix by k vectors, with coefficients of bA
 * and bB bits: nv^2 k products, or, by transforms if they are used with
 * coefficients of bT bits for the matrix, the transforms of the vectors
 * and the inverse ones of the products, those of the matrix if set */
static double hensel_mv_cost(const hensel_sys_t *hs, const long D,
                             const double bA, const double bB, const double bT,
                             const double k, const int set) {
  const double nv = hs->nv;
  const long np = hensel_ntt_np(hs, D, bT, bB);
  if (np == 0) {
    return nv * nv * k * hensel_mul_cost(D, bA, bB);
  }
  return (set ? nv * nv * hensel_ntt_cost(D, np, 0) : 0)
    + k * nv * (hensel_ntt_cost(D, np, 0) + hensel_ntt_cost(D, np, 1)
                + nv * hensel_ntt_cost(D, np, 2));
}

/* model cost of one Newton step from p^o to p^n, with J^-1 known modulo
 * p^jp and updated if upd, following hensel_step and hensel_solve; lp =
 * log2(p); jt is set if the transforms of J^-1 are formed */
static double hensel_model_step(const hensel_sys_t *hs, const long D,
                                const double lp, const long o, const long n,
                                const long jp, const int upd, const int jt) {
  const double nv = hs->nv;
  const long h = n - o;
  const double b = n * lp;
  /* F(X) and the two preinverses of w */
  double c = hensel_eval_cost(hs, 0, D, b) + 2 * hensel_mm_cost(D, b);
  if (jp < h) {
    c += hensel_eval_cost(hs, 1, D, h * lp);
  }
  /* the rounds of hensel_solve: y_e = J^-1 s modulo p^a, then J y_e
   * modulo p^(h - e), with the preinverses and reductions of each one; with
   * corrections, the products may be formed by transforms, those of J^-1
   * once for the steps without update, those of J once for the rounds */
  const int corr = jp < h;
  const double bj = jp * lp;
  for (long e = 0; e < h; e += jp) {
    const double bm = (h - e) * lp, ba = MIN(jp, h - e) * lp;
    c += hensel_mm_cost(D, ba) + 2 * nv * hensel_mul_cost(D, ba, ba)
      + (corr ? hensel_mv_cost(hs, D, ba, ba, bj, 1, e == 0 && !jt)
         : nv * nv * hensel_mul_cost(D, ba, ba));
    if (e + jp < h) {
      c += hensel_mm_cost(D, bm) + 2 * nv * hensel_mul_cost(D, bm, bm)
        + hensel_mv_cost(hs, D, bm, ba, h * lp, 1, e == 0);
    }
  }
  /* update of X and w */
  c += (nv + 1) * hensel_mm_cost(D, b);
  if (upd) {
    c += hensel_eval_cost(hs, 1, D, b)
      + 2 * (hensel_mv_cost(hs, D, b, b, b, nv, 1)
             + 2 * nv * nv * hensel_mul_cost(D, b, b));
  }
  return c;
}

/* whether J^-1 is updated in the step from p^o to p^n: while it is
 * current and the precision does not exceed jcap */
static int hensel_upd(const long jp, const long o, const long n,
                      const long jcap) {
  return jp >= o && n <= jcap;
}

/* model cost of the lifting of one prime from p^o (J^-1 known modulo
 * p^jp) to p^N by doublings, the last one possibly shorter */
static double hensel_model_from(const hensel_sys_t *hs, const long D,
                                const double lp, long o, long jp,
                                const long N, const long jcap) {
  double c = 0;
  int jt = 0;
  while (o < N) {
    const long n = MIN(2 * o, N);
    const int upd = hensel_upd(jp, o, n, jcap);
    c += hensel_model_step(hs, D, lp, o, n, jp, upd, jt);
    jt = !upd && (jt || jp < n - o);
    if (upd) {
      jp = MIN(2 * jp, n);
    }
    o = n;
  }
  return c;
}

/* model cost of the lifting of one prime to p^N with J^-1 updated up to
 * the precision N / 2^skip, including J mod p, its inverse by Gauss-Jordan
 * elimination and B_v = X_v w' mod w */
static double hensel_model(const hensel_sys_t *hs, const long D,
                           const double lp, const long N, const long skip) {
  const double nv = hs->nv;
  return nv * nv * nv * hensel_mm_cost(D, lp) + hensel_eval_cost(hs, 1, D, lp)
    + hensel_model_from(hs, D, lp, 1, 1, N, N >> skip)
    + nv * hensel_mm_cost(D, N * lp);
}

/* costs of a modular run in units of the model, a bit operation of the
 * products of Newton's iteration, about 5e-11 s (measured 4.3e-11 to
 * 6.7e-11 s for the products modulo w for D = 32 to 4096 and 60 to 3900
 * bits): a multiply-add of F4 by a known pivot, 1.5e-10 s, one of the
 * estimate of the multiply-adds of F4 among the new rows, 2.6e-10 s, and a
 * term of the F4 matrices, 2.2e-9 s, all recorded in the trace (see
 * construct_trace), a multiply-add of the products by the dense rows of
 * the multiplication matrix in FGLM, 7.5e-11 s, and one of the products
 * by the dense matrix of the multiplication by the separating form on
 * the dense route of separating.c, 3.8e-10 s (fitted on the modular runs
 * of 48 systems, the dense route on 4 of them); --hensel-cost-scale
 * multiplies them */
#define HENSEL_KAPPA_F4 3.0
#define HENSEL_KAPPA_F4NEW 5.2
#define HENSEL_KAPPA_F4TERM 44.0
#define HENSEL_KAPPA_FGLM 1.5
#define HENSEL_KAPPA_DENSE 7.6


/* the lifting state across the batches of primes */
#define HENSEL_MAX_MISMATCH 3
#define HENSEL_NHIST 256
#define HENSEL_NMAX 1024          /* largest precision p^N */
#define HENSEL_ATTEMPT_GROWTH 10  /* see hensel_state_step */
#define HENSEL_MIN_GAIN 1.2       /* below, the usual lifting is used */

typedef struct {
  hensel_sys_t hs;
  long ng;                        /* 1 once started, -1 before */
  split_data_t sp;                /* one group, see splitting.c */
  fmpz *M;                        /* modulus of the CRT */
  long *need;                     /* primes used when it was checked */
  long nused, nmismatch;
  mp_limb_t pchk;                 /* prime of the checks */
  /* bits of the modulus and fraction of the coefficients of g
   * reconstructed at each attempt, and the estimated number of bits */
  long nh;
  double hb[HENSEL_NHIST], hf[HENSEL_NHIST];
  double best;
  double lp;                      /* log2 of the primes */
  fmpz_t vden;                    /* denominator of g from the vector
                                     witnesses, 0 if none */
  long nmem;                      /* largest precision of one prime within
                                     the memory limit */
  double memref;                  /* estimate of the memory of the usual
                                     computation with one thread, in
                                     bytes */
  double memf4;                   /* its part for the F4 data */
  int nthr;                       /* threads, one prime each */
  /* residues of g and the B_v modulo M, packed in one block of
   * resL limbs each, which grows as a whole (each residue stored apart
   * would grow at each batch, which fragments the heap) */
  ulong *res;
  slong resL;
  double lastb;                   /* bits at the last attempt */
  double vbits;                   /* bits of the last attempt of LLL */
  double kappa;                   /* factor on the cost of a modular run */
  double memuser;                 /* memory limit in bytes, 0 if computed */
} hensel_state_t;

/* returns 0 if there are fewer equations than variables; the options
 * hensel_* of gens are read (see main.c) and param is the parametrization
 * at the first prime. With more equations than variables, Newton's
 * iteration uses a square subsystem of the sparsest equations whose
 * Jacobian is invertible modulo (w, p) at param (see hensel_select), if
 * there is one: its solutions contain those of the system, at which it is
 * then invertible over Q, and it is lifted as the combinations of all the
 * equations are, with fewer terms when the equations share few monomials
 * (eco 12 times the ideal of a point: 89 terms instead of 5600) */
static int hensel_state_init(hensel_state_t *st, const data_gens_ff_t *gens,
                             const param_t *param) {
  if (!hensel_sys_init(&st->hs, gens)) {
    return 0;
  }
  hensel_sys_t *hs = &st->hs;
  long *sub = NULL;
  double rt = realtime();
  if (hs->nall > hs->nv) {
    hensel_data_t hd;
    hensel_data_init(&hd);
    hensel_add_param(&hd, param);
    sub = hensel_select(hs, hd.pr);
    hensel_data_clear(&hd, hs->nv);
  }
  hs->square = sub != NULL;
  hs->ntt = 0;
  hensel_sys_set(hs, sub);
  free(sub);
  hs->tsel = realtime() - rt;
  const long N = gens->hensel_prec;
  st->kappa = gens->hensel_cost_scale > 0 ? gens->hensel_cost_scale : 1.0;
  st->memuser = gens->hensel_memory > 0 ? gens->hensel_memory * 1048576.0
    : 0;
  /* N <= 0: chosen by the cost model after the first modular runs */
  hs->nauto = N <= 0;
  hs->N = N > 0 ? N : 8;
  hs->nvar = 1;
  hs->cnewton = hs->cmod = 0;
  st->ng = -1;
  st->M = NULL;
  st->need = NULL;
  st->nused = 0;
  st->nmismatch = 0;
  st->nh = 0;
  st->best = -1;
  st->lp = 30;
  fmpz_init(st->vden);
  st->res = NULL;
  st->resL = 0;
  st->vbits = 0;
  st->nmem = HENSEL_NMAX;
  st->memref = 0;
  st->memf4 = 0;
  st->nthr = 1;
  st->lastb = 0;
  st->pchk = n_nextprime(UWORD(1) << 29, 1);
  hs->skip = 0;
  return 1;
}

/* clears the state; the groups in st->sp are cleared if clear_groups */
static void hensel_state_clear(hensel_state_t *st, const int clear_groups) {
  fmpz_clear(st->vden);
  free(st->res);
  if (st->ng > 0) {
    _fmpz_vec_clear(st->M, st->ng);
    free(st->need);
    if (clear_groups) {
      split_data_clear(&st->sp);
    }
  }
  hensel_sys_clear(&st->hs);
}

/* the group of the parametrization, from the images at the first prime */
static void hensel_state_start(hensel_state_t *st, const hensel_prime_t *P) {
  const long ng = P->ng, nv = st->hs.nv;
  st->ng = ng;
  st->M = _fmpz_vec_init(ng);
  st->need = calloc(ng, sizeof(long));
  split_data_t *sp = &st->sp;
  sp->nc = nv - 1;
  sp->nwords = (sp->nc * (sp->nc - 1) / 2 + 63) / 64;
  sp->ngroups = ng;
  sp->groups = malloc(ng * sizeof(split_group_t));
  fmpz_init(sp->modulus);
  sp->nprimes = 0;
  sp->ndone = 0;
  for (long r = 0; r < ng; r++) {
    split_group_t *G = sp->groups + r;
    G->deg = P->gdeg[r];
    G->sig = calloc(MAX(sp->nwords, 1), sizeof(uint64_t));
    G->n = G->deg * nv;
    G->crt = _fmpz_vec_init(G->n);
    G->rec = _fmpq_vec_init(G->n);
    G->state = 0;
    fmpz_one(st->M + r);
  }
}

/* ------------------------------------------------------------------ */
/* lifting of the parametrization, step by step                         */

/* starts the lifting of P: data modulo p and J^-1 modulo (w, p) */
static void hensel_nf_start(hensel_prime_t *P, const hensel_sys_t *hs) {
  const long nv = hs->nv;
  const mp_limb_t p = P->p;
  fmpz_t q;
  fmpz_init_set_ui(q, p);
  fmpz_mod_ctx_init(P->nfctx, q);
  fmpz_clear(q);
  fmpz_mod_ctx_struct *ctx = P->nfctx;
  hensel_lift_t *L = malloc(sizeof(hensel_lift_t));
  P->nfL = L;
  L->tf = L->tmv = L->tj = L->tmm = L->tup = 0;
  L->jprec = 1;
#if HENSEL_NTT
  L->ntt = NULL;
  if (hs->ntt) {
    L->ntt = malloc(sizeof(hensel_ntt_t));
    hensel_ntt_init(L->ntt);
  }
  L->Jh.d = NULL;
  L->Jh.np = 0;
#endif
  P->prec = 1;
  P->tnewton = 0;
  double rt = realtime();
  fmpz_mod_poly_init(L->f, ctx);
  fmpz_mod_poly_init(L->finv, ctx);
  fmpz_mod_poly_set_nmod_poly(L->f, P->w);
  L->X = malloc(nv * sizeof(fmpz_mod_poly_struct));
  L->Jinv = malloc(nv * nv * sizeof(fmpz_mod_poly_struct));
  for (long i = 0; i < nv; i++) {
    fmpz_mod_poly_init(L->X + i, ctx);
    fmpz_mod_poly_set_nmod_poly(L->X + i, P->X + i);
  }
  hensel_set_finv(L, ctx);
  nmod_poly_struct *A = malloc(nv * nv * sizeof(nmod_poly_struct));
  nmod_poly_struct *B = malloc(nv * nv * sizeof(nmod_poly_struct));
  fmpz_mod_poly_struct *e = malloc(nv * nv * sizeof(fmpz_mod_poly_struct));
  for (long i = 0; i < nv * nv; i++) {
    fmpz_mod_poly_init(e + i, ctx);
  }
  hensel_eval_set(e, hs, 1, L, ctx);
  fmpz_poly_t ez;
  fmpz_poly_init(ez);
  for (long i = 0; i < nv * nv; i++) {
    nmod_poly_init(A + i, p);
    nmod_poly_init(B + i, p);
    fmpz_mod_poly_get_fmpz_poly(ez, e + i, ctx);
    fmpz_poly_get_nmod_poly(A + i, ez);
    fmpz_mod_poly_clear(e + i, ctx);
  }
  free(e);
  fmpz_poly_clear(ez);
  double ri = realtime();
  P->nfok = hensel_inv_mod_p(B, A, P->w, nv);
  P->tinv = realtime() - ri;
  for (long i = 0; i < nv * nv; i++) {
    fmpz_mod_poly_init(L->Jinv + i, ctx);
    fmpz_mod_poly_set_nmod_poly(L->Jinv + i, B + i);
    nmod_poly_clear(A + i);
    nmod_poly_clear(B + i);
  }
  free(A);
  free(B);
  P->tnewton += realtime() - rt;
}

/* one step towards p^target; J^-1 is updated up to the precision jcap */
static void hensel_nf_step(hensel_prime_t *P, const hensel_sys_t *hs,
                           const long target, const long jcap) {
  hensel_lift_t *L = P->nfL;
  const long o = P->prec, n = MIN(2 * o, target);
  double rt = realtime();
  fmpz_t q;
  fmpz_init_set_ui(q, P->p);
  fmpz_pow_ui(q, q, n);
  fmpz_mod_ctx_set_modulus(P->nfctx, q);
  fmpz_clear(q);
  hensel_step(L, hs, hensel_upd(L->jprec, o, n, jcap), o, n, P->nfctx);
  P->prec = n;
  P->tnewton += realtime() - rt;
  P->tf = L->tf;
  P->tmv = L->tmv;
  P->tj = L->tj;
  P->tmm = L->tmm;
  P->tup = L->tup;
}

/* the images of the parametrization at the current precision, as one
 * group: g = w and B_v = X_v w' mod w modulo q = p^prec */
static void hensel_nf_images(hensel_prime_t *P, const hensel_sys_t *hs) {
  const long nv = hs->nv;
  hensel_lift_t *L = P->nfL;
  fmpz_mod_ctx_struct *ctx = P->nfctx;
  if (P->ng < 0) {
    P->ng = 1;
    P->gdeg = calloc(1, sizeof(long));
    P->g = malloc(sizeof(fmpz_poly_struct));
    fmpz_poly_init(P->g);
    P->B = malloc((nv - 1) * sizeof(fmpz_poly_struct));
    for (long v = 0; v < nv - 1; v++) {
      fmpz_poly_init(P->B + v);
    }
    fmpz_init(P->q);
  }
  fmpz_set(P->q, fmpz_mod_ctx_modulus(ctx));
  P->gdeg[0] = fmpz_mod_poly_degree(L->f, ctx);
  fmpz_mod_poly_get_fmpz_poly(P->g, L->f, ctx);
  fmpz_mod_poly_t wd, xr;
  fmpz_mod_poly_init(wd, ctx);
  fmpz_mod_poly_init(xr, ctx);
  fmpz_mod_poly_derivative(wd, L->f, ctx);
  for (long v = 0; v < nv - 1; v++) {
    HENSEL_MULMOD(xr, L->X + v, wd, L, ctx);
    fmpz_mod_poly_get_fmpz_poly(P->B + v, xr, ctx);
  }
  fmpz_mod_poly_clear(wd, ctx);
  fmpz_mod_poly_clear(xr, ctx);
}

static void hensel_nf_clear(hensel_prime_t *P, const hensel_sys_t *hs) {
  hensel_lift_t *L = P->nfL;
  if (L == NULL) {
    return;
  }
  fmpz_mod_ctx_struct *ctx = P->nfctx;
  for (long i = 0; i < hs->nv; i++) {
    fmpz_mod_poly_clear(L->X + i, ctx);
  }
  for (long i = 0; i < hs->nv * hs->nv; i++) {
    fmpz_mod_poly_clear(L->Jinv + i, ctx);
  }
  free(L->X);
  free(L->Jinv);
#if HENSEL_NTT
  hensel_ntt_mat_clear(&L->Jh);
  if (L->ntt != NULL) {
    hensel_ntt_clear(L->ntt);
    free(L->ntt);
  }
#endif
  fmpz_mod_poly_clear(L->f, ctx);
  fmpz_mod_poly_clear(L->finv, ctx);
  fmpz_mod_ctx_clear(ctx);
  free(L);
  P->nfL = NULL;
}

/* vector witnesses (see vec_wit_try in msolve.c): the coefficients of
 * the monic g share a common denominator d; k of them, u_1, .., u_k
 * modulo M, give the short vector (d, n_1, .., n_k), n_j = d u_j mod M,
 * found by LLL with about bits(n) + bits(d) / k bits instead of bits(n) +
 * bits(d). The candidate is accepted under the same bound as in msolve,
 * then used as the initial denominator of the reconstruction of g, which
 * fails if it is wrong. Returns 1 if a candidate is found */
#define HENSEL_VW_K 4
static int hensel_vecwit(fmpz_t d, const hensel_crt_t *C) {
  const long dg = C->dg;
  const fmpz *M = C->MQ;
  long idx[HENSEL_VW_K], k = 0;
  fmpz *c;
  fmpz_t t;
  const long cand[HENSEL_VW_K] = {dg - 1, dg / 2, dg / 4, dg / 16};
  if (dg < 4 * HENSEL_VW_K) {
    return 0;
  }
  c = _fmpz_vec_init(HENSEL_VW_K);
  fmpz_init(t);
  for (long j = 0; j < HENSEL_VW_K; j++) {
    long i = cand[j];
    int again = 1;
    while (again && i < dg) {
      hensel_crt_get(c + k, C, i, t);
      again = fmpz_is_zero(c + k);
      for (long l = 0; l < k; l++) {
        again = again || idx[l] == i;
      }
      i += again;
    }
    if (i < dg) {
      idx[k++] = i;
    }
  }
  fmpz_clear(t);
  if (k < 2) {
    _fmpz_vec_clear(c, HENSEL_VW_K);
    return 0;
  }
  fmpz_mat_t B;
  fmpz_mat_init(B, k + 1, k + 1);
  fmpz_one(fmpz_mat_entry(B, 0, 0));
  for (long j = 0; j < k; j++) {
    fmpz_set(fmpz_mat_entry(B, 0, j + 1), c + j);
    fmpz_set(fmpz_mat_entry(B, j + 1, j + 1), M);
  }
  fmpz_lll_t fl;
  fmpz_lll_context_init_default(fl);
  fmpz_lll(B, NULL, fl);
  const slong bits = (fmpz_bits(M) * k) / (k + 1) - HENSEL_RRMARGIN;
  int ok = !fmpz_is_zero(fmpz_mat_entry(B, 0, 0)) && bits > 0;
  for (long j = 0; j <= k && ok; j++) {
    ok = fmpz_bits(fmpz_mat_entry(B, 0, j)) <= bits;
  }
  if (ok) {
    fmpz_abs(d, fmpz_mat_entry(B, 0, 0));
  }
  fmpz_mat_clear(B);
  _fmpz_vec_clear(c, HENSEL_VW_K);
  return ok;
}

/* estimate of the number of bits needed, from the attempts: the bits
 * needed by the coefficients of g, from the highest degree down, increase
 * up to a peak near the middle; the fraction f(b) reconstructed with b
 * bits is extrapolated linearly from the last two attempts to f = 1/2,
 * or a bit further when it is already reached. Once g is reconstructed,
 * the B_v need a few percent more bits (4% on coplanarity): 4% and two
 * digits are added. -1 if unknown */
static double hensel_nf_estimate(const hensel_state_t *st) {
  const long nh = st->nh;
  if (nh >= 1 && st->hf[nh - 1] >= 1) {
    return 1.04 * st->hb[nh - 1] + 2 * st->lp;
  }
  if (nh < 2) {
    return -1;
  }
  const double b1 = st->hb[nh - 1], f1 = st->hf[nh - 1];
  long j = nh - 2;
  while (j >= 0 && st->hf[j] >= f1) {
    j--;
  }
  if (j < 0 || f1 <= 0) {
    return -1;
  }
  const double b0 = st->hb[j], f0 = st->hf[j];
  const double fpk = f1 < 0.5 ? 0.5 : MIN(1.0, f1 + 0.1);
  return b1 + (fpk - f1) * (b1 - b0) / (f1 - f0);
}

/* largest precision of one prime within the memory limit, with J^-1
 * updated up to the precision N / 2^skip: the lifting holds at full
 * precision the values of the monomials of F (mt.nlive[0], see
 * hensel_eval_set) and about 2 nv^2 + 4 nv polynomials of degree D (X,
 * F(X), dX, J^-1 F(X) and the corrections), and at most at half precision
 * (or at full precision without skip) those of J, J itself, J^-1 and the
 * operands and products of hensel_matmul (about 8 nv^2 polynomials); the
 * polynomials of F(X) or of J(X) whose products are not reduced (see
 * hensel_eval_set) have twice as many coefficients of twice as many bits,
 * about 3 more polynomials each, and the product 4 more; the transforms
 * of a matrix (see hensel_ntt_mul) take 16 D bytes per prime of about 50
 * bits and entry, for the bits of its products plus about 25, without the
 * headers of the coefficients: those of J and J^-1 in the steps without
 * update, with (N/2 + N/2^skip) and 2 N/2^skip digits, and those of the
 * products of the update, with twice the digits of J; a
 * coefficient of b bits takes 1.2 b / 8 + HENSEL_MEM_COEF bytes (fmpz,
 * mpz and malloc headers; sot1, D = 92: 132 MB at p^12 for 49 MB counted
 * before). The limit is given by --hensel-memory in megabytes, by default
 * st->memref (see hensel_set_memref), whose F4 part is counted once per
 * thread when all is set, since the usual computation holds the F4 data
 * of one prime per thread, plus HENSEL_MEM_CRT times the size
 * of the residues of the parametrization at the estimated number of
 * bits, which the usual lifting would also hold (chandra 11: 267 MB, 103
 * MB of residues), at least HENSEL_MEM_MIN megabytes. The estimate is
 * bounded by twice the bits bnow of the modulus, since the first ones are
 * often too high (katsura 13: 165000 bits for 20600), so that the limit
 * grows with the batches. It is not taken from the machine, so that the
 * computation does not depend on it: higher precisions save few primes
 * (katsura 13: 686 s at p^43, 684 s at p^87; eco 13: 86 s at p^44, 88 s
 * at p^87, 79 s at p^160). The limit for all the primes of a batch (all
 * set) is shared among them; that for one prime alone (all not set) does
 * not depend on the number of threads, so that neither does the choice of
 * the lifting (see hensel_choose) */
#define HENSEL_MEM_MIN 32
#define HENSEL_MEM_CRT 2
#define HENSEL_MEM_COEF 48
static double hensel_nf_mem(const hensel_state_t *st, const long D,
                            const double bnow, const int all) {
  if (st->memuser > 0) {
    return all ? st->memuser : st->memuser / st->nthr;
  }
  const double est = hensel_nf_estimate(st);
  const double B = est >= 0 ? MIN(est, 2 * bnow) : bnow;
  const double mem = st->memref + (all ? (st->nthr - 1) * st->memf4 : 0)
    + HENSEL_MEM_CRT * (double)D * st->hs.nv * B / 8;
  return MAX(mem, HENSEL_MEM_MIN * 1048576.0);
}

static long hensel_nf_nmem_skip(const hensel_state_t *st, const long D,
                                const double mem, const long skip) {
  const hensel_sys_t *hs = &st->hs;
  const double nv = hs->nv;
  const double npf = hs->mt.nlive[0] + 2 * nv * nv + 4 * nv
    + (hs->mt.delay[0] ? 3 * nv + 4 : 0);
  const double npj = (hs->mt.nlive[1] + 8 * nv * nv
                      + (hs->mt.delay[1] ? 3 * nv * nv + 4 : 0))
    / (skip > 0 ? 2 : 1);
  /* bytes per coefficient of the transforms, per bit of p^N and constant */
  const double tb = hensel_ntt_np(hs, D, 1, 1) > 0 ? 16 * nv * nv / 49.9 : 0;
  const double tfn = skip > 0 ? tb * (0.5 + 3 * ldexp(1, -skip)) : 0;
  const double tfc = skip > 0 ? tb * 2 * 25 : 0;
  const double tjn = tb * 2 / (skip > 0 ? 2 : 1), tjc = tb * 25;
  const double c = 1.2 / 8, lp = st->lp;
  const double Nf = (mem / D - HENSEL_MEM_COEF * npf - tfc)
    / (lp * (c * npf + tfn));
  const double Nj = (mem / D - HENSEL_MEM_COEF * npj - tjc)
    / (lp * (c * npj + tjn));
  return MAX(2, (long)MIN(Nf, Nj));
}

/* the largest precision for the T primes of a batch within the memory
 * limit, with or without the products by transforms, whichever gives the
 * better gain at this limit: hensel_choose decides with the limit for one
 * prime before the lifting, which grows with the bits of the residues
 * (random2_10_10: p^7 with them, p^12 without at the start, then p^128
 * and p^32 at most, so that they save a third of the primes) */
#define HENSEL_NTT_SWITCH 5
static long hensel_nf_nmem(hensel_state_t *st, const long D,
                           const double bnow, const long T,
                           const int info_level) {
  hensel_sys_t *hs = &st->hs;
  const double mem = hensel_nf_mem(st, D, bnow, 1);
  if (hs->nauto && hs->nvar == 2) {
    const int old = hs->ntt;
    double eff[2];
    for (int v = 0; v < 2; v++) {
      hs->ntt = v;
      /* at the best precision up to the limit: a step to a precision just
       * above a power of 2 costs about a full step (KLYZ13: 24.3 s for
       * p^130, 18 s for p^128) */
      const long Nv = MAX(2, MIN(hs->vN[v],
                                 hensel_nf_nmem_skip(st, D, mem,
                                                     hs->vskip[v]) / T));
      eff[v] = 0;
      for (long N0 = MAX(2, Nv / 2); N0 <= Nv; N0++) {
        eff[v] = MAX(eff[v], N0 / (hs->cmod + hensel_model(hs, D, st->lp, N0,
                                                           hs->vskip[v])));
      }
    }
    /* the other one only if it gains HENSEL_NTT_SWITCH percent, since the
     * limit also decreases with the estimate of the bits
     * (cp_d_3_n_7_p_7: 217, 241, 212, 281 MB) */
    const int bv = eff[!old] > (1 + HENSEL_NTT_SWITCH / 100.0) * eff[old]
      ? !old : old;
    hs->ntt = bv;
    hs->N = hs->vN[bv];
    hs->skip = hs->vskip[bv];
    if (info_level && bv != old) {
      fprintf(VERBSTREAM, "products by J and J^-1 %s from now on (memory "
              "limit %.0f MB)\n", bv ? "by transforms" : "of polynomials",
              mem / 1048576);
    }
  }
  return hensel_nf_nmem_skip(st, D, mem, hs->skip);
}

/* chooses N (a power of 2) and skip maximizing the precision per cost of
 * a prime, N / (cost of a modular run + cost of the lifting), where a
 * modular run costs the f4ops multiply-adds by the known pivots, the
 * f4new estimated multiply-adds among the new rows and the f4terms terms
 * of the F4 matrices, and the fglmops operations of FGLM or the denseops
 * ones of the dense route, weighted as above; N and skip given by the
 * user are kept. The gain, which decides whether the lifting is used, is
 * the best precision per cost with N at most the precision allowed by the
 * memory limit at the start (see hensel_nf_nmem), which grows with the
 * estimated number of bits: the run aims at N (random2_10_10: 319 primes
 * at about p^12, the limit at the start, 117 at p^32) */
static void hensel_choose(hensel_state_t *st, const long D, const double lp,
                          const double f4ops, const double f4new,
                          const double f4terms, const double fglmops,
                          const double denseops, const int info_level) {
  hensel_sys_t *hs = &st->hs;
  st->lp = lp;
  const double mem = hensel_nf_mem(st, D, 0, 0);
  const double kappa = st->kappa;
  hs->cmod = kappa * (HENSEL_KAPPA_F4 * f4ops + HENSEL_KAPPA_F4NEW * f4new
                       + HENSEL_KAPPA_F4TERM * f4terms
                       + HENSEL_KAPPA_FGLM * fglmops
                       + HENSEL_KAPPA_DENSE * denseops);
  /* with and without the products by transforms, which take memory and
   * may lower the precision within the memory limit (eco 12: p^7 instead
   * of p^12): the one with the best gain within the memory limit */
  double best = -1, bestc = -1, cnewton = 0;
  long bN = hs->N, bs = hs->skip, bNc = hs->N;
  int bntt = 0;
  const long Nin = hs->N;
  for (int ntt = 0; ntt < 2; ntt++) {
    hs->ntt = ntt;
    if (ntt && hensel_ntt_np(hs, D, 1, 1) == 0) {
      break;
    }
    double bt = -1, btc = -1, ct = 0;
    long tN = Nin, ts = hs->skip, tNc = Nin;
    for (long N0 = 2; N0 <= 1024; N0 *= 2) {
      if (!hs->nauto && N0 != Nin) {
        continue;
      }
      for (long sk = 0; (1L << sk) <= N0; sk++) {
        const double c = hensel_model(hs, D, lp, N0, sk);
        const double eff = N0 / (hs->cmod + c);
        if (eff > bt) {
          bt = eff;
          tN = N0;
          ts = sk;
          ct = c;
        }
        const long Nc = hs->nauto
          ? MIN(N0, hensel_nf_nmem_skip(st, D, mem, sk)) : N0;
        const double effc = Nc / (hs->cmod + hensel_model(hs, D, lp, Nc, sk));
        if (effc > btc) {
          btc = effc;
          tNc = Nc;
        }
      }
    }
    hs->vN[ntt] = tN;
    hs->vskip[ntt] = ts;
    hs->nvar = ntt + 1;
    if (btc > bestc) {
      best = bt;
      bestc = btc;
      bN = tN;
      bs = ts;
      bNc = tNc;
      cnewton = ct;
      bntt = ntt;
    }
  }
  hs->ntt = bntt;
  hs->cnewton = cnewton;
  hs->N = bN;
  hs->skip = bs;
  hs->gain = bestc * hs->cmod;
  if (info_level) {
    /* the costs are in the unit of the model of Newton's iteration; the
     * gain is the precision per cost of the lifting, relative to that of
     * the modular runs (one digit per run) */
    long nt[2] = {0, 0};
    for (long i = 0; i < hs->nv; i++) {
      nt[0] += hs->sys[i].nt;
    }
    for (long i = 0; i < hs->nv * hs->nv; i++) {
      nt[1] += hs->jac[i].nt;
    }
    fprintf(VERBSTREAM, "\n-------------- COST MODEL ---------------\n");
    fprintf(VERBSTREAM, "%-27s %10ld in %ld variables\n", "equations",
            hs->nall, hs->nv);
    if (hs->nall > hs->nv) {
      fprintf(VERBSTREAM, "%-27s %s (%.2f sec)\n", "Newton's iteration on",
              hs->square ? "a square subsystem of them" : "combinations of "
              "all of them, no square subsystem with an invertible Jacobian "
              "found", hs->tsel);
    }
    for (int w = 0; w < 2; w++) {
      fprintf(VERBSTREAM, "%-27s %10ld terms, %ld monomials, at most %ld "
              "values held\n", w == 0 ? "F (Newton's iteration)"
              : "Jacobian of F", nt[w], hs->mt.nord[w], hs->mt.nlive[w]);
      fprintf(VERBSTREAM, "%-27s %10ld products without successor, %s "
              "%ld sums\n", "", hs->mt.nleaf[w], hs->mt.delay[w]
              ? "reduced in" : "each reduced, not in", hs->mt.nred[w]);
    }
#if HENSEL_NTT
    fprintf(VERBSTREAM, "%-27s %s\n", "products by J and J^-1",
            hs->ntt ? "by transforms" : D < HENSEL_NTT_MIND
            ? "of polynomials (low degree)"
            : "of polynomials (better within the memory limit)");
#endif
    fprintf(VERBSTREAM, "%-27s %10.3g operations, %.3g among new rows, "
            "%.3g matrix terms\n", "F4 of a modular run", f4ops, f4new,
            f4terms);
    if (denseops > 0) {
      fprintf(VERBSTREAM, "%-27s %10.3g operations\n",
              "dense route of a modular run", denseops);
    } else {
      fprintf(VERBSTREAM, "%-27s %10.3g operations\n",
              "FGLM of a modular run", fglmops);
    }
    fprintf(VERBSTREAM, "%-27s %10.3g\n", "cost of a modular run", hs->cmod);
    fprintf(VERBSTREAM, "%-27s %10.3g to p^%ld", "cost of lifting one prime",
            hs->cnewton, hs->N);
    if (hs->skip > 0) {
      fprintf(VERBSTREAM, ", last %ld step%s without update of J^-1",
              hs->skip, hs->skip > 1 ? "s" : "");
    }
    fprintf(VERBSTREAM, "\n");
    fprintf(VERBSTREAM, "%-27s %10.2f at p^%ld\n", "gain", best * hs->cmod,
            hs->N);
    fprintf(VERBSTREAM, "%-27s %10.2f at p^%ld (threshold %.2f)\n",
            "gain within memory limit", hs->gain, bNc, HENSEL_MIN_GAIN);
  }
}

/* estimate of the memory of the usual computation, the reference for the
 * memory used by the lifting: the F4 data, about HENSEL_MEM_TERM bytes
 * per term of the largest matrix, and the coefficients preallocated by
 * the usual lifting, 28 D^2 bytes (see initialize_mpz_param in msolve.c),
 * not by this one. The peak of the usual computation is 0.9 to 2.3 times
 * this estimate on katsura 10 and 13, eco 12 and 13 and coplanarity, with
 * one thread; nthr threads hold the F4 data of nthr primes */
#define HENSEL_MEM_TERM 8
static void hensel_set_memref(hensel_state_t *st, const int64_t nterms,
                              const long D, const int nthr) {
  st->memf4 = HENSEL_MEM_TERM * (double)nterms;
  st->memref = st->memf4 + 28.0 * D * D;
  st->nthr = MAX(1, nthr);
}

/* target precision of the T primes of the current batch, at precision
 * prec, given the bits b_done of the primes already used: the remaining
 * bits are split over as few batches as possible with at most hs->N
 * (the most efficient precision of the cost model) per prime; returns
 * prec if the batch should stop there */
static long hensel_nf_target(const hensel_state_t *st, const long T,
                             const long prec, const double b_done) {
  const hensel_sys_t *hs = &st->hs;
  const long Neff = MIN(hs->N, st->nmem / T);
  const double B = hensel_nf_estimate(st);
  if (B < 0) {
    return MAX(prec, Neff);
  }
  long need = (long)ceil((B - b_done) / (T * st->lp)) + 1;
  /* the estimate often grows at each attempt, before g is reconstructed
   * (an extrapolation) and after (the B_v may need many more bits than g):
   * steps to slightly higher precisions would each cost a full step of
   * Newton's iteration, whose cost hardly depends on the gain of precision
   * (sot1: p^16, 19, 23, 27, 32, 2.1 times the cost of the doublings;
   * KLYZ12: p^33, 77, 120, 175, 246, 271, 297), so the precision is at
   * least doubled; only the last step of the lifting may overshoot */
  need = MAX(need, MIN(2 * prec, Neff));
  if (need <= prec) {
    /* estimate too low: a quarter more */
    need = prec + MAX(1, prec / 4);
  }
  const long k = (need + Neff - 1) / Neff;
  long target = (need + k - 1) / k;
  target = MIN(target, MIN(HENSEL_NMAX, MAX(st->nmem / T, prec + 1)));
  if (k > 1 && target <= prec) {
    return prec;
  }
  return MAX(target, prec + 1);
}

/* the legend of the lines printed by hensel_state_step, whose columns
 * span HENSEL_LINE characters */
#define HENSEL_LINE 79
static void hensel_print_legend(void) {
  fprintf(VERBSTREAM, "\nLegend for the p-adic lifting\n");
  fprintf(VERBSTREAM, "--------------------------------------------------"
          "------\n");
  fprintf(VERBSTREAM, "prec      precision p^k of the primes of the batch\n");
  fprintf(VERBSTREAM, "bits      bits of the modulus (primes used and "
          "batch)\n");
  fprintf(VERBSTREAM, "est       estimated number of bits needed\n");
  fprintf(VERBSTREAM, "rec       coefficients of w reconstructed, out of its "
          "degree\n");
  fprintf(VERBSTREAM, "vw        bits of the denominator given by the vector "
          "witnesses\n");
  fprintf(VERBSTREAM, "maxprec   largest precision within the memory limit\n");
  fprintf(VERBSTREAM, "mem       memory limit in MB\n");
  fprintf(VERBSTREAM, "crt       time of the CRT and rational reconstruction "
          "(sec)\n");
  fprintf(VERBSTREAM, "chk       time of the check modulo another prime "
          "(sec)\n");
  fprintf(VERBSTREAM, "Newton    time of Newton's iteration at each prime "
          "(sec), of which\n");
  fprintf(VERBSTREAM, "          inv  J^-1 mod p          F    F(X)\n");
  fprintf(VERBSTREAM, "          JF   J^-1 F(X)           J    J(X)\n");
  fprintf(VERBSTREAM, "          pr   products for J^-1   up   update of X "
          "and w\n");
  fprintf(VERBSTREAM, "--------------------------------------------------"
          "------\n");
  fprintf(VERBSTREAM, "\n  prec      bits       est         rec      vw "
          "maxprec     mem     crt     chk\n");
  for (int i = 0; i < HENSEL_LINE; i++) {
    fputc('-', VERBSTREAM);
  }
  fputc('\n', VERBSTREAM);
  fflush(VERBSTREAM);
}

/* lifts the batch hd step by step, trying the reconstruction after each
 * step with the primes already used; returns 1 when the parametrization
 * is reconstructed and checked, 0 if more primes are needed, -1 if the
 * lifting fails; the elapsed times of Newton's iteration, of the CRT and
 * rational reconstruction and of the checks are added to tnewton, trec
 * and tchk */
static int hensel_state_step(hensel_state_t *st, hensel_data_t *hd,
                             const int nthreads, double *tnewton,
                             double *trec, double *tchk,
                             const int info_level) {
  hensel_sys_t *hs = &st->hs;
  const long nv = hs->nv;
  long T = hd->nprimes;
  double rt = realtime();
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for (long j = 0; j < T; j++) {
    hensel_nf_start(hd->pr + j, hs);
  }
  /* primes where J is singular are dropped */
  long T2 = 0;
  for (long j = 0; j < T; j++) {
    if (hd->pr[j].nfok) {
      if (T2 != j) {
        hensel_prime_t tmp = hd->pr[T2];
        hd->pr[T2] = hd->pr[j];
        hd->pr[j] = tmp;
      }
      T2++;
    } else {
      st->nmismatch++;
    }
  }
  T = T2;
  if (T == 0) {
    *tnewton += realtime() - rt;
    return st->nmismatch >= HENSEL_MAX_MISMATCH ? -1 : 0;
  }
  if (info_level) {
    /* a line of the width of the header of the columns */
    char line[HENSEL_LINE + 1];
    int l;
    if (T == 1) {
      l = snprintf(line, sizeof(line), "---- prime %ld: p = %lu ",
                   st->nused + 1, (unsigned long)hd->pr[0].p);
    } else {
      l = snprintf(line, sizeof(line), "---- primes %ld to %ld: p = %lu to "
                   "%lu ", st->nused + 1, st->nused + T,
                   (unsigned long)hd->pr[0].p, (unsigned long)hd->pr[T - 1].p);
    }
    for (l = MIN(l, HENSEL_LINE); l < HENSEL_LINE; l++) {
      line[l] = '-';
    }
    line[HENSEL_LINE] = '\0';
    fprintf(VERBSTREAM, "%s\n", line);
    fflush(VERBSTREAM);
  }
  st->lp = log2((double)hd->pr[0].p);
  if (st->ng < 0) {
    hensel_nf_images(hd->pr, hs);
    hensel_state_start(st, hd->pr);
    st->nmem = hensel_nf_nmem(st, st->sp.groups->deg, 0, T, info_level);
  }
  *tnewton += realtime() - rt;
  split_group_t *G = st->sp.groups;
  const long dg = G->deg, n = G->n;
  hensel_crt_t C;
  int done = 0, attempted = 0;
  for (;;) {
    const double b_done = fmpz_is_one(st->M) ? 0 : fmpz_bits(st->M);
    const long prec = hd->pr[0].prec;
    long target = hensel_nf_target(st, T, prec, b_done);
    if (target <= prec) {
      if (attempted) {
        break;
      }
      target = prec + 1;
    }
    attempted = 1;
    const long jcap = target >> hs->skip;
    rt = realtime();
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (long j = 0; j < T; j++) {
      hensel_nf_step(hd->pr + j, hs, target, jcap);
      hensel_nf_images(hd->pr + j, hs);
    }
    *tnewton += realtime() - rt;
    /* attempt with the primes used and the batch at this precision, when
     * the modulus has grown by HENSEL_ATTEMPT_GROWTH percent since the
     * last attempt, or at the end of the batch if it reaches the estimated
     * number of bits: with many small batches, an attempt at the end of
     * each one would make the cost of the attempts quadratic (335 attempts
     * for 326 primes at p^8 on random2_9_9) */
    const double bnow = b_done + T * hd->pr[0].prec * st->lp;
    const long tnext = hensel_nf_target(st, T, hd->pr[0].prec, b_done);
    const double est = hensel_nf_estimate(st);
    if (bnow < (1 + HENSEL_ATTEMPT_GROWTH / 100.0) * st->lastb
        && (tnext > hd->pr[0].prec || est < 0 || bnow < est)) {
      continue;
    }
    st->lastb = bnow;
    /* the residues modulo the primes used and the batch are formed by the
     * reconstruction, which stops at the first failure */
    double ar = 0, ac = 0;
    rt = realtime();
    hensel_crt_init(&C, st->M, hd->pr, T, dg, nv);
    C.res = st->res;
    C.L = st->resL;
    const fmpz *M2 = C.MQ;
    long nrec = 0, nrec2 = 0;
    done = hensel_reconstruct(G->rec, &C, NULL, &nrec);
    ar += realtime() - rt;
    if (done) {
      rt = realtime();
      done = hensel_check(G->rec, dg, hs, st->pchk) == 1;
      ac += realtime() - rt;
    }
    rt = realtime();
    /* vector witnesses, tried when the modulus has grown by 5% */
    if (!done && fmpz_is_zero(st->vden) && nrec < dg
        && fmpz_bits(M2) >= 1.05 * st->vbits) {
      st->vbits = fmpz_bits(M2);
      hensel_vecwit(st->vden, &C);
    }
    if (!done && !fmpz_is_zero(st->vden)) {
      done = hensel_reconstruct(G->rec, &C, st->vden, &nrec2);
      ar += realtime() - rt;
      if (done) {
        rt = realtime();
        done = hensel_check(G->rec, dg, hs, st->pchk) == 1;
        ac += realtime() - rt;
      }
      rt = realtime();
      if (!done && nrec2 < nrec) {
        /* probably a wrong candidate */
        fmpz_zero(st->vden);
      }
      nrec = MAX(nrec, nrec2);
    }
    if (!done) {
      /* the partial reconstruction would stay scattered in the heap */
      for (long i = 0; i < n; i++) {
        fmpq_zero(G->rec + i);
      }
    }
    if (st->nh < HENSEL_NHIST) {
      st->hb[st->nh] = fmpz_bits(M2);
      st->hf[st->nh] = (double)nrec / dg;
      st->nh++;
    }
    /* the limit follows the estimated number of bits */
    st->nmem = hensel_nf_nmem(st, dg, bnow, T, info_level);
    const long bM2 = fmpz_bits(M2);
    hensel_crt_clear(&C);
    ar += realtime() - rt;
    *trec += ar;
    *tchk += ac;
    if (info_level) {
      char sp[24], se[24], sr[48], sv[24], sm[24];
      const double est = hensel_nf_estimate(st);
      snprintf(sp, sizeof(sp), "p^%ld", hd->pr[0].prec);
      if (est < 0) {
        snprintf(se, sizeof(se), "-");
      } else {
        snprintf(se, sizeof(se), "%.0f", est);
      }
      snprintf(sr, sizeof(sr), "%ld/%ld", done ? dg : nrec, dg);
      if (fmpz_is_zero(st->vden)) {
        snprintf(sv, sizeof(sv), "-");
      } else {
        snprintf(sv, sizeof(sv), "%ld", (long)fmpz_bits(st->vden));
      }
      snprintf(sm, sizeof(sm), "p^%ld", st->nmem);
      fprintf(VERBSTREAM, "%6s %9ld %9s %11s %7s %7s %7.0f %7.2f %7.2f%s\n",
              sp, bM2, se, sr, sv, sm,
              hensel_nf_mem(st, dg, bnow, 1) / 1048576.0, ar, ac,
              done ? "  reconstructed" : "");
      fflush(VERBSTREAM);
    }
    if (done) {
      break;
    }
  }
  if (info_level) {
    for (long j = 0; j < T; j++) {
      const hensel_prime_t *P = hd->pr + j;
      if (T > 1) {
        fprintf(VERBSTREAM, "p = %lu: ", (unsigned long)P->p);
      }
      fprintf(VERBSTREAM, "Newton p^%ld %.2f: inv %.2f, F %.2f, JF %.2f, "
              "J %.2f, pr %.2f, up %.2f\n", P->prec, P->tnewton, P->tinv,
              P->tf, P->tmv, P->tj, P->tmm, P->tup);
    }
  }
  /* the batch is used: the block of residues is enlarged, then they are
   * updated from the last one, so that the residue i, read at i L, is
   * written at i L' >= i L only after the residues above it */
  rt = realtime();
  hensel_crt_init(&C, st->M, hd->pr, T, dg, nv);
  if (!done) {
    const slong L = st->resL;
    const slong L2 = (fmpz_bits(C.MQ) + FLINT_BITS - 1) / FLINT_BITS;
    st->res = realloc(st->res, n * L2 * sizeof(ulong));
    C.res = st->res;
    C.L = L;
    fmpz_t x, t;
    fmpz_init(x);
    fmpz_init(t);
    for (long i = n - 1; i >= 0; i--) {
      hensel_crt_get(x, &C, i, t);
      fmpz_get_ui_array(st->res + i * L2, L2, x);
    }
    fmpz_clear(x);
    fmpz_clear(t);
    st->resL = L2;
  }
  fmpz_set(st->M, C.MQ);
  hensel_crt_clear(&C);
  if (done) {
    /* the residues are no longer needed */
    free(st->res);
    st->res = NULL;
    st->resL = 0;
  }
  *trec += realtime() - rt;
  st->nused += T;
  st->sp.nprimes = st->nused;
  if (done) {
    G->state = 2;
    st->need[0] = st->nused;
    st->sp.ndone = 1;
  }
  for (long j = 0; j < hd->nprimes; j++) {
    hensel_nf_clear(hd->pr + j, hs);
  }
  /* FLINT keeps the integers freed by the lifting in a cache of each
   * thread, with up to 64 limbs (p^87 on katsura 13 needs 42): 600 MB
   * which would stay scattered in the heap */
#pragma omp parallel num_threads(nthreads)
  _fmpz_cleanup_mpz_content();
  return done;
}
