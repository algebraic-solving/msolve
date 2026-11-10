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

/* Refinement of isolating intervals with Newton's iteration in ball
 * arithmetic.
 *
 * The interval of a root of f is (numer / 2^k, (numer + 1) / 2^k) when
 * k >= 0 and (numer, numer + 2^(-k)) when k < 0. Newton's iteration on
 * f / f' is run from the middle of the interval, on approximations whose
 * precision doubles at each step. Its result is only a candidate: the
 * refined interval (N / 2^K, (N + 1) / 2^K) is certified by the signs of f
 * at its ends, evaluated in ball arithmetic until the balls exclude 0, and
 * by its inclusion in the interval being refined. When Newton's iteration
 * fails, four bisection steps are performed before it is tried again.
 *
 * The evaluations lose bits by cancellation: the terms of f at x are much
 * larger than f(x). The precision used is the number of correct bits
 * wanted plus extra, where extra is the number of bits lost at the last
 * evaluation whose sign was found. It is shared by the roots of a
 * polynomial, since it depends mainly on the size of its coefficients. */

#include "flint/fmpz_vec.h"
#include "flint/fmpq.h"
#include "flint/arb.h"
#include "flint/arb_fmpz_poly.h"

/* sign of f at n / 2^K (K >= 0); the precision is doubled until the ball
 * excludes 0, up to the size of the exact value, beyond which the sign is
 * computed in exact arithmetic */
static int usolve_sign_2exp_arb(const fmpz_poly_t f, const fmpz_t n,
                                const slong K, slong *extra){
  if (fmpz_poly_length(f) == 0) {
    return 0;
  }
  arb_t x, y;
  arb_init(x);
  arb_init(y);
  arb_set_fmpz(x, n);
  arb_mul_2exp_si(x, x, -K);
  /* size of the exact value */
  const slong deg = fmpz_poly_degree(f);
  const slong pex = FLINT_ABS(_fmpz_vec_max_bits(f->coeffs, f->length))
    + deg * FLINT_MAX(K, (slong)fmpz_bits(n)) + FLINT_BIT_COUNT(deg) + 64;
  slong p = K + *extra + 32;
  int s;
  for (;;) {
    arb_fmpz_poly_evaluate_arb(y, f, x, p);
    if (arb_is_positive(y)) {
      s = 1;
      break;
    }
    if (arb_is_negative(y)) {
      s = -1;
      break;
    }
    if (p >= pex) {
      fmpq_t xq, v;
      fmpq_init(xq);
      fmpq_init(v);
      fmpz_set(fmpq_numref(xq), n);
      fmpz_one(fmpq_denref(xq));
      fmpz_mul_2exp(fmpq_denref(xq), fmpq_denref(xq), K);
      fmpq_canonicalise(xq);
      fmpz_poly_evaluate_fmpq(v, f, xq);
      s = fmpq_sgn(v);
      fmpq_clear(xq);
      fmpq_clear(v);
      break;
    }
    p = FLINT_MIN(2 * p, pex);
  }
  if (s) {
    *extra = FLINT_MAX(32, p - K - 32);
  }
  arb_clear(x);
  arb_clear(y);
  return s;
}

/* sign of f just right of n / 2^K */
static int usolve_sign_right_2exp_arb(const fmpz_poly_t f, const fmpz_poly_t df,
                                      const fmpz_t n, const slong K,
                                      slong *extra){
  int s = usolve_sign_2exp_arb(f, n, K, extra);
  if (s) {
    return s;
  }
  /* n / 2^K is a root of f, which is square-free */
  return usolve_sign_2exp_arb(df, n, K, extra);
}

/* Newton's iteration from the middle of (lo / 2^kb, (lo + 1) / 2^kb) to
 * an approximation x of the root with about K + 4 correct bits; on
 * success, returns 1 and N = floor(x 2^K) */
static int usolve_newton_arb(fmpz_t N, const fmpz_poly_t f,
                             const fmpz_poly_t df, const fmpz_t lo,
                             const slong kb, const slong K, slong *extra){
  arb_t x, fx, dfx, step;
  arb_init(x);
  arb_init(fx);
  arb_init(dfx);
  arb_init(step);
  fmpz_t m;
  fmpz_init(m);
  fmpz_mul_2exp(m, lo, 1);
  fmpz_add_ui(m, m, 1);
  arb_set_fmpz(x, m);
  arb_mul_2exp_si(x, x, -(kb + 1));

  int ok = 1;
  slong acc = kb + 1;
  const slong maxit = 2 * FLINT_BIT_COUNT(K) + 16;
  for (slong it = 0; acc < K + 4; it++) {
    if (it == maxit) {
      ok = 0;
      break;
    }
    const slong target = FLINT_MIN(2 * acc, K + 8);
    const slong p = target + *extra + 32;
    arb_fmpz_poly_evaluate_arb(fx, f, x, p);
    arb_fmpz_poly_evaluate_arb(dfx, df, x, p);
    if (arb_contains_zero(dfx)) {
      ok = 0;
      break;
    }
    arb_div(step, fx, dfx, p);
    if (!arb_is_finite(step)) {
      ok = 0;
      break;
    }
    if (arb_contains_zero(step) &&
        mag_cmp_2exp_si(arb_radref(step), -target) > 0) {
      /* f(x) is not known well enough: more precision */
      *extra *= 2;
      continue;
    }
    arb_sub(x, x, step, p);
    arb_get_mid_arb(x, x);
    mag_t a;
    mag_init(a);
    arb_get_mag(a, step);
    if (mag_cmp_2exp_si(a, -kb) > 0) {
      /* the step leaves the interval */
      mag_clear(a);
      ok = 0;
      break;
    }
    /* the error after a step of size 2^-e is about 2^-2e */
    const slong e = mag_is_zero(a) ? -(K + 8) : MAG_EXP(a);
    mag_clear(a);
    acc = FLINT_MIN(FLINT_MAX(acc + 1, -2 * e - 4), target);
  }
  if (ok) {
    arf_t a;
    arf_init(a);
    arf_mul_2exp_si(a, arb_midref(x), K);
    arf_get_fmpz(N, a, ARF_RND_FLOOR);
    arf_clear(a);
  }
  fmpz_clear(m);
  arb_clear(x);
  arb_clear(fx);
  arb_clear(dfx);
  arb_clear(step);
  return ok;
}

/* rt becomes the exact root n / 2^K */
static void usolve_set_exact_root(interval *rt, const fmpz_t n, const slong K){
  fmpz_get_mpz(rt->numer, n);
  rt->k = K;
  rt->isexact = 1;
}

/* Refines the isolating interval rt of a root of the square-free
 * polynomial f, of derivative df, to rt->k >= K, or until the root is
 * found to be exact. extra is the number of bits lost by cancellation in
 * the evaluations (see above), updated by the refinement. */
void usolve_refine_root_arb(interval *rt, const fmpz_poly_t f,
                            const fmpz_poly_t df, const slong K, slong *extra){
  if (rt->isexact == 1 || rt->k >= K) {
    return;
  }
  slong ex = *extra;
  fmpz_t lo, t;
  fmpz_init(lo);
  fmpz_init(t);
  fmpz_set_mpz(lo, rt->numer);
  slong kb = rt->k;

  /* sign of f in the interval, left of the root */
  const int slo = usolve_sign_right_2exp_arb(f, df, lo, FLINT_MAX(kb, 0), &ex);
  /* (lo, lo + 2^-kb) with kb < 0: bisection down to kb = 0 */
  while (kb < 0) {
    fmpz_one(t);
    fmpz_mul_2exp(t, t, -kb - 1);
    fmpz_add(t, lo, t);
    const int s = usolve_sign_2exp_arb(f, t, 0, &ex);
    if (s == 0) {
      usolve_set_exact_root(rt, t, 0);
      goto done;
    }
    if (s == slo) {
      fmpz_swap(lo, t);
    }
    kb++;
  }
  fmpz_t N, u;
  fmpz_init(N);
  fmpz_init(u);
  while (kb < K) {
    if (usolve_newton_arb(N, f, df, lo, kb, K, &ex)) {
      /* the candidate must lie in the current interval */
      fmpz_mul_2exp(t, lo, K - kb);
      fmpz_add_ui(u, lo, 1);
      fmpz_mul_2exp(u, u, K - kb);
      for (int tries = 0; tries < 4; tries++) {
        if (fmpz_cmp(N, t) < 0 || fmpz_cmp(N, u) >= 0) {
          break;
        }
        /* the ends t / 2^K and u / 2^K of the interval may be roots of f
         * outside it: f has the sign slo right of the first one and -slo
         * left of the second one */
        const int s1 = fmpz_equal(N, t) ? slo
          : usolve_sign_2exp_arb(f, N, K, &ex);
        if (s1 == 0) {
          usolve_set_exact_root(rt, N, K);
          goto clear;
        }
        fmpz_add_ui(N, N, 1);
        const int s2 = fmpz_equal(N, u) ? -slo
          : usolve_sign_2exp_arb(f, N, K, &ex);
        if (s2 == 0) {
          usolve_set_exact_root(rt, N, K);
          goto clear;
        }
        fmpz_sub_ui(N, N, 1);
        if (s1 != s2) {
          fmpz_get_mpz(rt->numer, N);
          rt->k = K;
          rt->sign_left = slo;
          goto clear;
        }
        if (s1 == slo) {
          fmpz_add_ui(N, N, 1);
        } else {
          fmpz_sub_ui(N, N, 1);
        }
      }
    }
    /* Newton's iteration failed: bisection steps */
    for (int i = 0; i < 4 && kb < K; i++) {
      fmpz_mul_2exp(lo, lo, 1);
      fmpz_add_ui(t, lo, 1);
      kb++;
      const int s = usolve_sign_2exp_arb(f, t, kb, &ex);
      if (s == 0) {
        usolve_set_exact_root(rt, t, kb);
        goto clear;
      }
      if (s == slo) {
        fmpz_swap(lo, t);
      }
    }
  }
  fmpz_get_mpz(rt->numer, lo);
  rt->k = kb;
  rt->sign_left = slo;
clear:
  fmpz_clear(N);
  fmpz_clear(u);
done:
  *extra = ex;
  fmpz_clear(lo);
  fmpz_clear(t);
}

/* Refines all the non-exact isolating intervals in roots, of the real
 * roots of upol (of degree deg, read in place), to rt->k >= prec, with
 * nthreads threads sharing f and f' (each root gets the same interval
 * whatever the precisions learned on the others, so the output does not
 * depend on the number of threads); the integers cached by FLINT in each
 * thread are then freed */
static void usolve_refine_roots_arb(mpz_t *upol, const unsigned long int deg,
                                    interval *roots,
                                    const unsigned long int nbroots,
                                    const slong prec, const int nthreads){
  /* f on aliases of the coefficients of upol */
  fmpz_poly_struct f[1];
  fmpz_poly_t df;
  f->coeffs = malloc((deg + 1) * sizeof(fmpz));
  for (unsigned long int i = 0; i <= deg; i++) {
    fmpz_init_set_readonly(f->coeffs + i, upol[i]);
  }
  f->alloc = deg + 1;
  f->length = deg + 1;
  _fmpz_poly_normalise(f);
  fmpz_poly_init(df);
  fmpz_poly_derivative(df, f);
#pragma omp parallel num_threads(nthreads)
  {
    slong extra = 64;
#pragma omp for schedule(dynamic, 1)
    for (unsigned long int i = 0; i < nbroots; i++) {
      usolve_refine_root_arb(roots + i, f, df, prec, &extra);
    }
    _fmpz_cleanup_mpz_content();
  }
  for (unsigned long int i = 0; i <= deg; i++) {
    fmpz_clear_readonly(f->coeffs + i);
  }
  free(f->coeffs);
  fmpz_poly_clear(df);
}
