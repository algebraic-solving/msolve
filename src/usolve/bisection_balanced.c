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

/* Descartes bisection with balanced fixed-point polynomials.
 *
 * The real roots of F (degree n, F(0) != 0) in (0, 1) are isolated by
 * bisection, as by bisection_rec. The node (c, k) stands for the interval
 * (c / 2^k, (c + 1) / 2^k) and its polynomial u(x) = F((c + x) / 2^k), up
 * to a positive factor. Its coefficients have sizes spreading over
 * thousands of bits, but this spread is almost linear in the degree: with a
 * slope s, the u_i 2^(s i) have close sizes. A node keeps integers V,
 * truncated to W bits relative to their largest one, error bounds E >= 0
 * and a slope s such that |V_i - lambda u_i 2^(s i)| <= E_i for some
 * lambda > 0. Then
 *   - the left half u(x / 2) has the same V and E, with slope s + 1;
 *   - the right half u((x + 1) / 2) is V(y + 2^-(s + 1)): a Taylor shift
 *     by a dyadic number, computed in fixed point by divide and conquer
 *     with the powers (y + 2^-t)^m truncated to 2^-q, which are cached, the
 *     error bounds following the same computation with upward rounding;
 *     it is then balanced again (new slope) and truncated to W bits.
 * Shifts by 2^-t multiply errors by up to (1 + 2^-t)^n: when the slope is
 * small, near the root of the tree, this is too much. So the top of the
 * tree is investigated by bisection_rec, in exact arithmetic, and the
 * subtree of a node whose best slope reaches USOLVE_BAL_SLOPE is
 * investigated here (usolve_bal_subtree). A node polynomial computed again
 * exactly (see below) stays exact (V = u, E = 0) in the same way.
 *
 * The Descartes test counts the sign variations of the reversed
 * polynomial shifted by 1, computed from V truncated to about n bits
 * relative to its largest coefficient; the error bound of the i-th
 * coefficient is bounded by the largest error bound of the coefficients
 * of index >= i times binom(n + 1, i + 1). A coefficient has a known sign
 * when it exceeds its error bound; the test is conclusive when all the
 * possible signs of the others lead to the same decision (no root, one
 * root, or a bisection, which is always safe). Otherwise, and when the
 * sign of u at 1/2 or 1 is needed but unknown, the node polynomial is
 * computed again exactly from F. W starts at about 1.5 n, which makes
 * these fallbacks rare on most inputs; it is doubled when they become
 * frequent (close roots, large coefficients).
 *
 * As in bisection_rec, the right half is not investigated when the node
 * has at most 2 sign variations and the left half as many (the numbers of
 * sign variations of the halves add up to at most the one of the node),
 * and a node with 2 sign variations whose polynomial has opposite signs
 * at 0 and 1/2 and at 1/2 and 1 gives two isolating intervals. A root at
 * the middle of a node is exact; it is not divided out: Descartes' rule
 * on the open halves is not affected by a root at their ends, whose zero
 * coefficients are handled by the exact fallback. */

#include <math.h>
#include "flint/fmpz_vec.h"

#define USOLVE_BAL_SLOPE 2   /* nodes with a smaller best slope stay exact */
#define USOLVE_BAL_BASE 64   /* Horner's scheme below this length */
#define USOLVE_BAL_EBITS 30  /* bits of the error bounds in the products */
#define USOLVE_BAL_NCACHE 4096

typedef struct {
  fmpz *v;    /* coefficients */
  fmpz *e;    /* error bounds */
  slong s;    /* slope */
  int exact;  /* v is the node polynomial (e = 0, s = 0) */
} usolve_bal_node_t;

typedef struct {
  slong m, t, q;
  fmpz_poly_t c;    /* floor(binom(m, j) 2^(q - t (m - j))) */
  fmpz_poly_t ce;   /* ceil((c + 1) / 2^(q - 64)), or c if q = 0 */
  slong z;          /* c_j = 0 for j < z */
} usolve_bal_power_t;

typedef struct {
  slong n;
  fmpz *f;              /* F, exact */
  slong w;              /* working precision */
  double *lbin;         /* log2 binom(n + 1, i + 1) */
  usolve_bal_power_t *powers;
  long npowers;
  interval *roots;
  unsigned long int *nbr;
  usolve_flags *flags;
  int stop;             /* a root was found and flags->hasrealroots is set */
  long fallbacks;
  long subtrees;        /* subtrees investigated in fixed point */
  slong wmax;           /* largest working precision */
  long nodes;           /* nodes investigated here */
  long recent_nodes;    /* nodes and fallbacks since w last changed */
  long recent_fallbacks;
} usolve_bal_t;

static void usolve_bal_node_init(usolve_bal_node_t *a, const slong n){
  a->v = _fmpz_vec_init(n + 1);
  a->e = _fmpz_vec_init(n + 1);
  a->s = 0;
  a->exact = 0;
}

static void usolve_bal_node_clear(usolve_bal_node_t *a, const slong n){
  _fmpz_vec_clear(a->v, n + 1);
  _fmpz_vec_clear(a->e, n + 1);
}

static void usolve_bal_node_set(usolve_bal_node_t *a,
                                const usolve_bal_node_t *b, const slong n){
  _fmpz_vec_set(a->v, b->v, n + 1);
  _fmpz_vec_set(a->e, b->e, n + 1);
  a->s = b->s;
  a->exact = b->exact;
}

/* index of (y + 2^-t)^m in the cache of B */
static long usolve_bal_power(usolve_bal_t *B, const slong m, const slong t,
                             const slong q){
  for(long i = 0; i < B->npowers; i++){
    if(B->powers[i].m == m && B->powers[i].t == t && B->powers[i].q == q){
      return i;
    }
  }
  if(B->npowers == USOLVE_BAL_NCACHE){
    fprintf(ERRSTREAM, "usolve: too many powers of (x + 2^-t)\n");
    exit(1);
  }
  const long id = B->npowers++;
  usolve_bal_power_t *P = B->powers + id;
  P->m = m;
  P->t = t;
  P->q = q;
  fmpz_poly_init(P->c);
  fmpz_poly_init(P->ce);
  fmpz_poly_fit_length(P->c, m + 1);
  fmpz_poly_fit_length(P->ce, m + 1);
  fmpz_t b;
  fmpz_init_set_ui(b, 1);
  for(slong j = m; j >= 0; j--){
    /* b = binom(m, j) */
    const slong x = q - t * (m - j);
    if(x >= 0){
      fmpz_mul_2exp(P->c->coeffs + j, b, x);
    }
    else{
      fmpz_fdiv_q_2exp(P->c->coeffs + j, b, -x);
    }
    if(q > 64){
      fmpz_add_ui(P->ce->coeffs + j, P->c->coeffs + j, 1);
      fmpz_cdiv_q_2exp(P->ce->coeffs + j, P->ce->coeffs + j, q - 64);
    }
    else{
      fmpz_set(P->ce->coeffs + j, P->c->coeffs + j);
    }
    if(j > 0){
      fmpz_mul_ui(b, b, j);
      fmpz_divexact_ui(b, b, m - j + 1);
    }
  }
  fmpz_clear(b);
  _fmpz_poly_set_length(P->c, m + 1);
  _fmpz_poly_normalise(P->c);
  _fmpz_poly_set_length(P->ce, m + 1);
  _fmpz_poly_normalise(P->ce);
  P->z = 0;
  while(P->z < P->c->length && fmpz_is_zero(P->c->coeffs + P->z)){
    P->z++;
  }
  return id;
}

static void usolve_bal_mul(fmpz *r, const fmpz *a, const slong la,
                           const fmpz *b, const slong lb){
  if(la >= lb){
    _fmpz_poly_mul(r, a, la, b, lb);
  }
  else{
    _fmpz_poly_mul(r, b, lb, a, la);
  }
}

/* v(y + 2^-t), and its error bounds e, in place (len coefficients); when
 * t > 0, products are truncated to 2^-q with q = w + 64 */
static void usolve_bal_shift(usolve_bal_t *B, fmpz *v, fmpz *e,
                             const slong len, const slong t, fmpz_t tmp){
  if(len <= USOLVE_BAL_BASE){
    for(slong i = 0; i + 1 < len; i++){
      for(slong j = len - 2; j >= i; j--){
        if(t >= 0){
          fmpz_fdiv_q_2exp(tmp, v + j + 1, t);
          fmpz_add(v + j, v + j, tmp);
          fmpz_cdiv_q_2exp(tmp, e + j + 1, t);
          fmpz_add(e + j, e + j, tmp);
          if(t > 0){
            fmpz_add_ui(e + j, e + j, 1);
          }
        }
        else{
          fmpz_mul_2exp(tmp, v + j + 1, -t);
          fmpz_add(v + j, v + j, tmp);
          fmpz_mul_2exp(tmp, e + j + 1, -t);
          fmpz_add(e + j, e + j, tmp);
        }
      }
    }
    return;
  }
  slong m = 1;
  while(2 * m < len){
    m *= 2;
  }
  usolve_bal_shift(B, v, e, m, t, tmp);
  usolve_bal_shift(B, v + m, e + m, len - m, t, tmp);
  /* v = v0 + (y + 2^-t)^m v1 */
  const slong q = (t > 0) ? B->w + 64 : 0;
  const usolve_bal_power_t *P = B->powers + usolve_bal_power(B, m, t, q);
  const slong l1 = len - m, lp = l1 + P->c->length - 1;
  fmpz *prod = _fmpz_vec_init(lp);
  fmpz *perr = _fmpz_vec_init(lp);
  if(P->z < P->c->length){
    usolve_bal_mul(prod + P->z, v + m, l1, P->c->coeffs + P->z,
                   P->c->length - P->z);
    /* error bounds of v1 with USOLVE_BAL_EBITS-bit mantissas */
    fmpz *em = _fmpz_vec_init(l1);
    slong mb = _fmpz_vec_max_bits(e + m, l1);
    slong g = (mb > USOLVE_BAL_EBITS) ? mb - USOLVE_BAL_EBITS : 0;
    for(slong i = 0; i < l1; i++){
      fmpz_cdiv_q_2exp(em + i, e + m + i, g);
    }
    usolve_bal_mul(perr + P->z, em, l1, P->ce->coeffs + P->z,
                   P->ce->length - P->z);
    if(g){
      _fmpz_vec_scalar_mul_2exp(perr, perr, lp, g);
    }
    _fmpz_vec_clear(em, l1);
  }
  if(q > 0){
    /* (|v1| + e1) / 2^q bounds the truncation of the powers (including
     * their coefficients below z), 1 the truncation of the products */
    fmpz_t a;
    fmpz_init(a);
    for(slong i = 0; i < l1; i++){
      fmpz_abs(tmp, v + m + i);
      fmpz_add(a, a, tmp);
      fmpz_add(a, a, e + m + i);
    }
    fmpz_cdiv_q_2exp(a, a, q);
    fmpz_add_ui(a, a, 1);
    for(slong i = 0; i < lp; i++){
      fmpz_cdiv_q_2exp(perr + i, perr + i, 64);
      fmpz_add(perr + i, perr + i, a);
      fmpz_fdiv_q_2exp(prod + i, prod + i, q);
    }
    fmpz_clear(a);
  }
  for(slong i = 0; i < len; i++){
    if(i < m){
      fmpz_add(v + i, v + i, prod + i);
      fmpz_add(e + i, e + i, perr + i);
    }
    else{
      fmpz_set(v + i, prod + i);
      fmpz_set(e + i, perr + i);
    }
  }
  _fmpz_vec_clear(prod, lp);
  _fmpz_vec_clear(perr, lp);
}

/* slope d (|d| <= 64) for which the v_i 2^(d i) have the closest sizes */
static slong usolve_bal_best_slope(const fmpz *v, const slong n){
  double best = 1e300;
  slong bd = 0;
  for(slong d = -64; d <= 64; d++){
    double mx = -1e300, mn = 1e300;
    for(slong i = 0; i <= n; i++){
      if(fmpz_is_zero(v + i)){
        continue;
      }
      const double x = fmpz_bits(v + i) + (double)d * i;
      mx = FLINT_MAX(mx, x);
      mn = FLINT_MIN(mn, x);
    }
    if(mx - mn < best){
      best = mx - mn;
      bd = d;
    }
  }
  return bd;
}

/* new slope, then truncation to w bits */
static void usolve_bal_rebalance(usolve_bal_t *B, usolve_bal_node_t *a){
  const slong n = B->n;
  const slong d = usolve_bal_best_slope(a->v, n);
  if(d){
    for(slong i = 0; i <= n; i++){
      if(d > 0){
        fmpz_mul_2exp(a->v + i, a->v + i, d * i);
        fmpz_mul_2exp(a->e + i, a->e + i, d * i);
      }
      else if(i > 0){
        fmpz_fdiv_q_2exp(a->v + i, a->v + i, -d * i);
        fmpz_cdiv_q_2exp(a->e + i, a->e + i, -d * i);
        fmpz_add_ui(a->e + i, a->e + i, 1);
      }
    }
    a->s += d;
  }
  const slong mb = FLINT_ABS(_fmpz_vec_max_bits(a->v, n + 1));
  if(mb > B->w){
    for(slong i = 0; i <= n; i++){
      fmpz_fdiv_q_2exp(a->v + i, a->v + i, mb - B->w);
      fmpz_cdiv_q_2exp(a->e + i, a->e + i, mb - B->w);
      fmpz_add_ui(a->e + i, a->e + i, 1);
    }
  }
}

/* divides v by the largest power of 2 dividing all its coefficients */
static void usolve_bal_remove_2content(fmpz *v, const slong n){
  slong z = -1;
  for(slong i = 0; i <= n; i++){
    if(!fmpz_is_zero(v + i)){
      const slong x = fmpz_val2(v + i);
      if(z < 0 || x < z){
        z = x;
      }
    }
  }
  if(z > 0){
    for(slong i = 0; i <= n; i++){
      fmpz_fdiv_q_2exp(v + i, v + i, z);
    }
  }
}

/* exact polynomial of the node (c, k): 2^(nk) F((c + x) / 2^k) */
static void usolve_bal_exact(usolve_bal_t *B, usolve_bal_node_t *a,
                             const fmpz_t c, const slong k){
  const slong n = B->n;
  for(slong i = 0; i <= n; i++){
    fmpz_mul_2exp(a->v + i, B->f + i, (n - i) * k);
    fmpz_zero(a->e + i);
  }
  _fmpz_poly_taylor_shift(a->v, c, n + 1);
  usolve_bal_remove_2content(a->v, n);
  a->s = 0;
  a->exact = 1;
  B->fallbacks++;
  /* frequent fallbacks: w is too small for the separation of the roots
   * and the size of the coefficients of F, it is doubled (the nodes
   * computed later from exact polynomials use it) */
  B->recent_fallbacks++;
  if(B->recent_fallbacks >= 4 && 8 * B->recent_fallbacks >= B->recent_nodes
     && B->w < B->wmax){
    B->w = FLINT_MIN(2 * B->w, B->wmax);
    B->recent_nodes = 0;
    B->recent_fallbacks = 0;
  }
}

/* sign of the i-th coefficient of a polynomial with error bounds:
 * 2 if unknown */
static inline int usolve_bal_sign(const fmpz_t x, const fmpz_t ex){
  if(fmpz_is_zero(ex)){
    return fmpz_sgn(x);
  }
  if(fmpz_cmpabs(x, ex) > 0){
    return fmpz_sgn(x);
  }
  return 2;
}

/* fewest (vmin) and most (vmax) sign variations of a sequence of signs
 * in {-1, 0, 1, 2}, 2 standing for an unknown sign */
static void usolve_bal_variations(const int *sg, const slong len,
                                  long *vmin, long *vmax){
  long vn = 0, vx = 0, run = 0;
  int last = 0;
  for(slong i = 0; i < len; i++){
    const int s = sg[i];
    if(s == 0){
      continue;
    }
    if(s == 2){
      run++;
      continue;
    }
    if(last == 0){
      vx += run;
    }
    else{
      /* run unknown signs between last and s: at most run + 1 sign
       * variations, with the parity of the variation from last to s */
      const long par = (s != last);
      long m = run + 1;
      vn += par;
      if((m & 1) != par){
        m--;
      }
      vx += m;
    }
    run = 0;
    last = s;
  }
  *vmin = vn;
  *vmax = vx + run;
}

/* Descartes test of the node: 0, 1, 2 (exactly 2 sign variations),
 * 3 (3 or more) or -1 (unknown); trunc = 0 asks for the exact test of an
 * exact node */
static long usolve_bal_test(usolve_bal_t *B, const usolve_bal_node_t *a,
                            const int trunc){
  const slong n = B->n;
  /* u_i = v_i 2^x_i: size mb of the largest coefficient of u */
  slong mb = 0;
  for(slong i = 0; i <= n; i++){
    if(!fmpz_is_zero(a->v + i)){
      const slong x = (a->s >= 0) ? a->s * (n - i) : -a->s * i;
      mb = FLINT_MAX(mb, (slong)fmpz_bits(a->v + i) + x);
    }
  }
  const slong keep = trunc ? n + 1 : mb;
  /* t = reverse of u / 2^(mb - keep) */
  fmpz *t = _fmpz_vec_init(n + 1);
  fmpz *et = _fmpz_vec_init(n + 1);
  for(slong i = 0; i <= n; i++){
    const slong j = n - i;
    const slong x = ((a->s >= 0) ? a->s * (n - j) : -a->s * j) - (mb - keep);
    if(x >= 0){
      fmpz_mul_2exp(t + i, a->v + j, x);
      fmpz_mul_2exp(et + i, a->e + j, x);
    }
    else{
      fmpz_fdiv_q_2exp(t + i, a->v + j, -x);
      fmpz_cdiv_q_2exp(et + i, a->e + j, -x);
      fmpz_add_ui(et + i, et + i, 1);
    }
  }
  fmpz_t one;
  fmpz_init_set_ui(one, 1);
  _fmpz_poly_taylor_shift(t, one, n + 1);
  fmpz_clear(one);
  /* error bound of the i-th coefficient of the shifted polynomial:
   * (max_{j >= i} et_j) binom(n + 1, i + 1), in log2 */
  int *sg = malloc((n + 1) * sizeof(int));
  double smax = -1;
  for(slong i = n; i >= 0; i--){
    if(!fmpz_is_zero(et + i)){
      smax = FLINT_MAX(smax, (double)fmpz_bits(et + i));
    }
    if(smax < 0){
      sg[i] = fmpz_sgn(t + i);
    }
    else if(!fmpz_is_zero(t + i)
            && (double)(fmpz_bits(t + i) - 1) > smax + B->lbin[i] + 1e-6){
      sg[i] = fmpz_sgn(t + i);
    }
    else{
      sg[i] = 2;
    }
  }
  long vmin, vmax;
  usolve_bal_variations(sg, n + 1, &vmin, &vmax);
  free(sg);
  _fmpz_vec_clear(t, n + 1);
  _fmpz_vec_clear(et, n + 1);
  if(vmax == 0){
    return 0;
  }
  if(vmin == 1 && vmax == 1){
    return 1;
  }
  if(vmin >= 2){
    return (vmax == 2) ? 2 : 3;
  }
  return -1;
}

/* sign of u at 2^(s - x): sum v_i 2^(-x i); 2 if unknown */
static int usolve_bal_sign_at(usolve_bal_t *B, const usolve_bal_node_t *a,
                              const slong x){
  const slong n = B->n;
  fmpz_t acc, err, tmp;
  fmpz_init(acc);
  fmpz_init(err);
  fmpz_init(tmp);
  for(slong i = 0; i <= n; i++){
    /* sum v_i 2^(x (n - i)) or sum v_i 2^(-x i) */
    if(x >= 0){
      fmpz_mul_2exp(acc, acc, x);
      fmpz_add(acc, acc, a->v + i);
      fmpz_mul_2exp(err, err, x);
      fmpz_add(err, err, a->e + i);
    }
    else{
      fmpz_mul_2exp(tmp, a->v + i, -x * i);
      fmpz_add(acc, acc, tmp);
      fmpz_mul_2exp(tmp, a->e + i, -x * i);
      fmpz_add(err, err, tmp);
    }
  }
  const int s = usolve_bal_sign(acc, err);
  fmpz_clear(acc);
  fmpz_clear(err);
  fmpz_clear(tmp);
  return s;
}

static void usolve_bal_root(usolve_bal_t *B, const fmpz_t c, const slong k,
                            const int isexact, const int sgnlft){
  usolve_flags *flags = B->flags;
  mpz_t z;
  mpz_init(z);
  fmpz_get_mpz(z, c);
  merge_root(B->roots, z, k, isexact, sgnlft, *(B->nbr),
             flags->bound_pos, flags->bound_neg, flags->sign);
  mpz_clear(z);
  (*(B->nbr))++;
  if(flags->verbose >= 1){
    fprintf(VERBSTREAM, "+");
    if((*(B->nbr)) % 100 == 0){
      fprintf(VERBSTREAM, "[%lu]", *(B->nbr));
    }
  }
  if(flags->hasrealroots == 1){
    B->stop = 1;
  }
}

/* investigates the node (c, k) of polynomial a0, which is not modified
 * (so that the left half can share it); returns its number of sign
 * variations (3 for 3 or more) */
static long usolve_bal_node(usolve_bal_t *B, const usolve_bal_node_t *a0,
                            const fmpz_t c, const slong k){
  const slong n = B->n;
  usolve_flags *flags = B->flags;
  (flags->node_looked)++;
  B->nodes++;
  B->recent_nodes++;
  /* a: a0, or own when the node polynomial is computed again */
  const usolve_bal_node_t *a = a0;
  usolve_bal_node_t own;
  int owned = 0;

  double e_time = realtime();
  long nb = usolve_bal_test(B, a, 1);
  if(nb < 0 && a->exact){
    nb = usolve_bal_test(B, a, 0);
  }
  if(nb < 0){
    usolve_bal_node_init(&own, n);
    owned = 1;
    usolve_bal_exact(B, &own, c, k);
    a = &own;
    nb = usolve_bal_test(B, a, 1);
    if(nb < 0){
      nb = usolve_bal_test(B, a, 0);
    }
  }
  flags->time_desc += (realtime() - e_time);
  if(nb <= 1){
    if(nb == 1){
      const int bot = usolve_bal_sign(a->v + 0, a->e + 0);
      usolve_bal_root(B, c, k, 0, (bot == 2) ? 0 : bot);
    }
    if(owned){
      usolve_bal_node_clear(&own, n);
    }
    return nb;
  }
  const int bot = usolve_bal_sign(a->v + 0, a->e + 0);

  /* signs of u at 1/2 and at 1 */
  int half = usolve_bal_sign_at(B, a, a->s + 1);
  int top = (nb == 2) ? usolve_bal_sign_at(B, a, a->s) : 0;
  if(half == 2 || top == 2 || (nb == 2 && bot == 2)){
    if(!owned){
      usolve_bal_node_init(&own, n);
      owned = 1;
    }
    usolve_bal_exact(B, &own, c, k);
    a = &own;
    half = usolve_bal_sign_at(B, a, 1);
    top = (nb == 2) ? usolve_bal_sign_at(B, a, 0) : 0;
  }
  const int sbot = usolve_bal_sign(a->v + 0, a->e + 0);

  fmpz_t c2;
  fmpz_init(c2);
  fmpz_mul_2exp(c2, c, 1);
  if(nb == 2 && sbot * half < 0 && top * half < 0){
    /* one root in each half */
    (flags->half_done)++;
    usolve_bal_root(B, c2, k + 1, 0, sbot);
    if(!B->stop){
      fmpz_add_ui(c2, c2, 1);
      usolve_bal_root(B, c2, k + 1, 0, half);
    }
    goto clear;
  }

  /* an exact node goes to fixed point once its best slope is large
   * enough for the shifts by 2^-(s + 1) to be well conditioned */
  if(a->exact && usolve_bal_best_slope(a->v, n) >= USOLVE_BAL_SLOPE){
    if(!owned){
      usolve_bal_node_init(&own, n);
      usolve_bal_node_set(&own, a, n);
      owned = 1;
    }
    usolve_bal_rebalance(B, &own);
    own.exact = 0;
    a = &own;
  }

  /* left half: the same data with slope s + 1, or, for an exact node,
   * u(x / 2) */
  long nbleft;
  if(a->exact){
    usolve_bal_node_t L;
    usolve_bal_node_init(&L, n);
    usolve_bal_node_set(&L, a, n);
    for(slong i = 0; i <= n; i++){
      fmpz_mul_2exp(L.v + i, a->v + i, n - i);
    }
    usolve_bal_remove_2content(L.v, n);
    nbleft = usolve_bal_node(B, &L, c2, k + 1);
    usolve_bal_node_clear(&L, n);
  }
  else{
    usolve_bal_node_t L = *a;
    L.s = a->s + 1;
    nbleft = usolve_bal_node(B, &L, c2, k + 1);
  }
  if(B->stop){
    goto clear;
  }

  /* the middle is a root */
  fmpz_add_ui(c2, c2, 1);
  if(half == 0){
    usolve_bal_root(B, c2, k + 1, 1, 0);
    if(B->stop){
      goto clear;
    }
  }

  if(nb <= 2 && nbleft == nb){
    (flags->skipped)++;
    if(flags->verbose >= 1){
      fprintf(VERBSTREAM, "s");
    }
    goto clear;
  }

  /* right half */
  usolve_bal_node_t R;
  usolve_bal_node_init(&R, n);
  usolve_bal_node_set(&R, a, n);
  e_time = realtime();
  if(flags->verbose >= 1){
    fprintf(VERBSTREAM, "*");
  }
  if(a->exact){
    for(slong i = 0; i <= n; i++){
      fmpz_mul_2exp(R.v + i, a->v + i, n - i);
    }
    fmpz_t one;
    fmpz_init_set_ui(one, 1);
    _fmpz_poly_taylor_shift(R.v, one, n + 1);
    fmpz_clear(one);
    usolve_bal_remove_2content(R.v, n);
  }
  else{
    fmpz_t tmp;
    fmpz_init(tmp);
    R.s = a->s + 1;
    usolve_bal_shift(B, R.v, R.e, n + 1, R.s, tmp);
    fmpz_clear(tmp);
    usolve_bal_rebalance(B, &R);
  }
  flags->time_shift += (realtime() - e_time);
  (flags->transl)++;
  usolve_bal_node(B, &R, c2, k + 1);
  usolve_bal_node_clear(&R, n);
clear:
  fmpz_clear(c2);
  if(owned){
    usolve_bal_node_clear(&own, n);
  }
  return nb;
}

static void usolve_bal_init(usolve_bal_t *B, mpz_t *upol,
                            const unsigned long int deg, interval *roots,
                            unsigned long int *nbr, usolve_flags *flags){
  const slong n = deg;
  B->n = n;
  B->f = _fmpz_vec_init(n + 1);
  for(slong i = 0; i <= n; i++){
    fmpz_set_mpz(B->f + i, upol[i]);
  }
  B->w = (3 * (n + 1)) / 2 + 64;
  B->wmax = 64 * (n + 1);
  B->nodes = 0;
  B->recent_nodes = 0;
  B->recent_fallbacks = 0;
  B->lbin = malloc((n + 1) * sizeof(double));
  for(slong i = 0; i <= n; i++){
    B->lbin[i] = (lgamma(n + 2.0) - lgamma(i + 2.0) - lgamma(n - i + 1.0))
      / log(2.0);
  }
  B->powers = malloc(USOLVE_BAL_NCACHE * sizeof(usolve_bal_power_t));
  B->npowers = 0;
  B->roots = roots;
  B->nbr = nbr;
  B->flags = flags;
  B->stop = 0;
  B->fallbacks = 0;
  B->subtrees = 0;
}

static void usolve_bal_clear(usolve_bal_t *B){
  if(B->flags->verbose >= 1){
    fprintf(VERBSTREAM, "\n%ld subtrees in fixed point (%ld nodes), %ld "
            "exact node polynomials computed again, final precision %ld "
            "bits\n", B->subtrees, B->nodes, B->fallbacks, (long)B->w);
  }
  for(long i = 0; i < B->npowers; i++){
    fmpz_poly_clear(B->powers[i].c);
    fmpz_poly_clear(B->powers[i].ce);
  }
  free(B->powers);
  free(B->lbin);
  _fmpz_vec_clear(B->f, B->n + 1);
}

/* called by bisection_rec on the node (c, k), whose exact polynomial is
 * upol, of degree deg: when the node polynomial still has the degree of
 * F (no exact root was divided out), upol(0) != 0, and its best slope is
 * at least USOLVE_BAL_SLOPE, its subtree is investigated in fixed point;
 * then 1 is returned, *nbv is the number of sign variations of the node,
 * and upol is unchanged. Otherwise, 0 is returned. */
static int usolve_bal_subtree(usolve_bal_t *B, mpz_t *upol,
                              const unsigned long int deg, const mpz_t c,
                              const long k, long *nbv){
  const slong n = B->n;
  if((slong)deg != n || mpz_sgn(upol[0]) == 0){
    return 0;
  }
  usolve_bal_node_t a;
  usolve_bal_node_init(&a, n);
  for(slong i = 0; i <= n; i++){
    fmpz_set_mpz(a.v + i, upol[i]);
  }
  if(usolve_bal_best_slope(a.v, n) < USOLVE_BAL_SLOPE){
    usolve_bal_node_clear(&a, n);
    return 0;
  }
  usolve_bal_rebalance(B, &a);
  B->subtrees++;
  fmpz_t cf;
  fmpz_init(cf);
  fmpz_set_mpz(cf, c);
  *nbv = usolve_bal_node(B, &a, cf, k);
  fmpz_clear(cf);
  usolve_bal_node_clear(&a, n);
  return 1;
}
