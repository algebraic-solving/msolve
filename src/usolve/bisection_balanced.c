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
#if defined(__GLIBC__)
#include <malloc.h>
#endif
#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <mach/mach.h>
#endif
#include "flint/fmpz_vec.h"

#define USOLVE_BAL_SLOPE 2   /* nodes with a smaller best slope stay exact */
#define USOLVE_BAL_BASE 64   /* Horner's scheme below this length */
#define USOLVE_BAL_EBITS 30  /* bits of the error bounds in the products */
#define USOLVE_BAL_NCACHE 4096
#define USOLVE_BAL_GRIDMAX 8  /* grids of sign evaluations up to 8 roots */
/* memory, measured on katsura-12 to katsura-14: each thread of
 * usolve_bal_run adds about USOLVE_BAL_THRMEM (n + 1) w bits to the
 * resident memory (mostly the work space of FLINT's products, which each
 * thread keeps), and the exact bisection reaches about USOLVE_BAL_TOPMEM
 * times the size of its largest node polynomial */
#define USOLVE_BAL_THRMEM 55
#define USOLVE_BAL_TOPMEM 18
/* the threads may always use USOLVE_BAL_INMEM times the size of F, about
 * twice the memory reached by the exact bisection on katsura-12 and
 * katsura-13: with a univariate input, nothing else is computed before */
#define USOLVE_BAL_INMEM 64

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

/* the cache of the powers, shared by the threads */
typedef struct {
  usolve_bal_power_t *p;
  long n;
} usolve_bal_cache_t;

typedef struct usolve_bal_s {
  slong n;
  /* F(y) = 2^d g(s 2^b y), with g = f, or f / x when the root 0 of the
   * input f was divided out (o = 1), and s = -1 for the negative roots;
   * f is read in place (fl coefficients, aliases of the input) */
  fmpz *f;
  slong fl, o, b, d;
  int s;
  slong w;              /* working precision */
  double *lbin;         /* log2 binom(n + 1, i + 1) */
  usolve_bal_cache_t *powers;
  /* subtrees left for usolve_bal_run: nodes (rebalanced, packed by
   * usolve_bal_pack), c and k */
  ulong **queue;
  fmpz *qc;
  slong *qk;
  long nqueue, aqueue;
  double smax;          /* bytes of the largest exact node polynomial */
  interval *roots;
  unsigned long int *nbr;
  usolve_flags *flags;
  int stop;             /* a root was found and flags->hasrealroots is set */
  long fallbacks;
  long subtrees;        /* subtrees investigated in fixed point */
  slong sign_extra;     /* bits lost by cancellation in the signs of F */
  long grids, grids_ok, grid_signs; /* grids tried, successful, signs */
  slong wmax;           /* largest working precision */
  long nodes;           /* nodes investigated here */
  long recent_nodes;    /* nodes and fallbacks since w last changed */
  long recent_fallbacks;
  /* in usolve_bal_run: the state where the tasks add their counters */
  struct usolve_bal_s *top;
  long tasks;           /* right halves investigated in their own task */
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
static long usolve_bal_power_new(usolve_bal_cache_t *C, const slong m,
                                 const slong t, const slong q);

/* the entries of the cache are never moved nor modified once computed;
 * the threads look them up and compute them one at a time */
static long usolve_bal_power(usolve_bal_t *B, const slong m, const slong t,
                             const slong q){
  usolve_bal_cache_t *C = B->powers;
  long id = -1;
#pragma omp critical(usolve_bal_power)
  {
    for(long i = 0; i < C->n && id < 0; i++){
      if(C->p[i].m == m && C->p[i].t == t && C->p[i].q == q){
        id = i;
      }
    }
    if(id < 0){
      id = usolve_bal_power_new(C, m, t, q);
    }
  }
  return id;
}

static long usolve_bal_power_new(usolve_bal_cache_t *C, const slong m,
                                 const slong t, const slong q){
  if(C->n == USOLVE_BAL_NCACHE){
    fprintf(ERRSTREAM, "usolve: too many powers of (x + 2^-t)\n");
    exit(1);
  }
  const long id = C->n++;
  usolve_bal_power_t *P = C->p + id;
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
  const usolve_bal_power_t *P = B->powers->p + usolve_bal_power(B, m, t, q);
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
    /* the coefficient of F times 2^((n - i) k), an integer */
    const slong x = B->b * i + B->d + (n - i) * k;
    if(x >= 0){
      fmpz_mul_2exp(a->v + i, B->f + B->o + i, x);
    }
    else{
      fmpz_fdiv_q_2exp(a->v + i, B->f + B->o + i, -x);
    }
    if(B->s < 0 && (i & 1)){
      fmpz_neg(a->v + i, a->v + i);
    }
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
 * 3 (3 or more) or -1 (unknown); *vmaxp is an upper bound on the number
 * of sign variations, hence on the number of roots of the node; trunc = 0
 * asks for the exact test of an exact node */
static long usolve_bal_test(usolve_bal_t *B, const usolve_bal_node_t *a,
                            const int trunc, long *vmaxp){
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
  *vmaxp = vmax;
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
  /* the tasks of usolve_bal_run share the array of the roots */
#pragma omp critical(usolve_bal_root)
  {
    merge_root(B->roots, z, k, isexact, sgnlft, *(B->nbr),
               flags->bound_pos, flags->bound_neg, flags->sign);
    (*(B->nbr))++;
    if(flags->verbose >= 1){
      fprintf(VERBSTREAM, "+");
      if((*(B->nbr)) % 100 == 0){
        fprintf(VERBSTREAM, "[%lu]", *(B->nbr));
      }
    }
  }
  mpz_clear(z);
  if(flags->hasrealroots == 1){
    B->stop = 1;
  }
}

/* sign of F at y = x / 2^K >= 0, from the signs of f (see usolve_bal_t) */
static int usolve_bal_sign_F(usolve_bal_t *B, const fmpz_poly_t f,
                             const fmpz_t x, const slong K){
  if(fmpz_is_zero(x)){
    return fmpz_sgn(B->f + B->o);
  }
  /* t = s 2^b y */
  fmpz_t t;
  fmpz_init(t);
  slong e = K - B->b;
  if(e >= 0){
    fmpz_set(t, x);
  }
  else{
    fmpz_mul_2exp(t, x, -e);
    e = 0;
  }
  if(B->s < 0){
    fmpz_neg(t, t);
  }
  int sg = usolve_sign_2exp_arb(f, t, e, &B->sign_extra);
  /* the sign of t when g = f / x */
  if(B->o && B->s < 0){
    sg = -sg;
  }
  fmpz_clear(t);
  return sg;
}

/* Signs of F at the points (c 2^q + j) / 2^(k + q), 0 <= j <= 2^q, of grids
 * with q from about log2(2 r) to log2(16 r), r = vmax, each point being
 * evaluated once (certified signs, see refine_arb.c). When r cells of a
 * grid have a sign change, each holds one root (the node has at most r):
 * they are given as isolating intervals and 1 is returned. Otherwise, or
 * if F vanishes at a point of a grid, 0 is returned. */
static int usolve_bal_grid(usolve_bal_t *B, const fmpz_t c, const slong k,
                           const long r){
  slong q0 = 1, q1;
  while((1L << q0) < 2 * r){
    q0++;
  }
  q1 = q0 + 3;
  const slong m = 1L << q1;
  /* f as a polynomial, on the coefficients of B */
  fmpz_poly_struct f;
  f.coeffs = B->f;
  f.alloc = B->fl;
  f.length = B->fl;
  int *sg = malloc((m + 1) * sizeof(int));
  for(slong j = 0; j <= m; j++){
    sg[j] = 2;
  }
  fmpz_t x;
  fmpz_init(x);
  int ok = 0;
  B->grids++;
  for(slong q = q0; q <= q1 && !ok; q++){
    const slong step = 1L << (q1 - q);
    long changes = 0;
    int last = 0, zero = 0;
    for(slong j = 0; j <= m; j += step){
      if(sg[j] == 2){
        fmpz_mul_2exp(x, c, q1);
        fmpz_add_ui(x, x, j);
        sg[j] = usolve_bal_sign_F(B, &f, x, k + q1);
        B->grid_signs++;
      }
      if(sg[j] == 0){
        zero = 1;
        break;
      }
      if(last && sg[j] != last){
        changes++;
      }
      last = sg[j];
    }
    if(zero){
      break;
    }
    if(changes == r){
      ok = 1;
      /* one root in each cell with a sign change */
      int prev = sg[0];
      for(slong j = step; j <= m && !B->stop; j += step){
        if(sg[j] != prev){
          fmpz_mul_2exp(x, c, q);
          fmpz_add_ui(x, x, j / step - 1);
          usolve_bal_root(B, x, k + q, 0, prev);
        }
        prev = sg[j];
      }
    }
  }
  if(ok){
    B->grids_ok++;
  }
  fmpz_clear(x);
  free(sg);
  return ok;
}

/* state of a task of usolve_bal_run: that of P (working precision, recent
 * nodes and fallbacks, ...), with its own flags fl (copied from pf) and its
 * own counters */
static void usolve_bal_task_init(usolve_bal_t *T, const usolve_bal_t *P,
                                 const usolve_flags *pf, usolve_flags *fl){
  *fl = *pf;
  fl->node_looked = 0;
  fl->half_done = 0;
  fl->skipped = 0;
  fl->transl = 0;
  fl->time_desc = 0;
  fl->time_shift = 0;
  *T = *P;
  T->flags = fl;
  T->nodes = 0;
  T->fallbacks = 0;
  T->grids = 0;
  T->grids_ok = 0;
  T->grid_signs = 0;
  T->tasks = 0;
}

/* adds the counters of the task T to T->top */
static void usolve_bal_task_done(const usolve_bal_t *T){
  usolve_bal_t *B = T->top;
  usolve_flags *flags = B->flags;
  const usolve_flags *fl = T->flags;
#pragma omp critical(usolve_bal_run)
  {
    B->w = FLINT_MAX(B->w, T->w);
    B->nodes += T->nodes;
    B->fallbacks += T->fallbacks;
    B->grids += T->grids;
    B->grids_ok += T->grids_ok;
    B->grid_signs += T->grid_signs;
    B->tasks += T->tasks;
    flags->node_looked += fl->node_looked;
    flags->half_done += fl->half_done;
    flags->skipped += fl->skipped;
    flags->transl += fl->transl;
    flags->time_desc += fl->time_desc;
    flags->time_shift += fl->time_shift;
  }
}

static long usolve_bal_node(usolve_bal_t *B, const usolve_bal_node_t *a0,
                            const fmpz_t c, const slong k);

/* investigates the right half (c2, k + 1) of a node (c, k) of polynomial
 * a, with c2 = 2 c + 1 */
static void usolve_bal_right(usolve_bal_t *B, const usolve_bal_node_t *a,
                             const fmpz_t c2, const slong k){
  const slong n = B->n;
  usolve_flags *flags = B->flags;
  usolve_bal_node_t R;
  usolve_bal_node_init(&R, n);
  usolve_bal_node_set(&R, a, n);
  double e_time = realtime();
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
  long vmax;
  long nb = usolve_bal_test(B, a, 1, &vmax);
  if(nb < 0 && a->exact){
    nb = usolve_bal_test(B, a, 0, &vmax);
  }
  if(nb < 0){
    usolve_bal_node_init(&own, n);
    owned = 1;
    usolve_bal_exact(B, &own, c, k);
    a = &own;
    nb = usolve_bal_test(B, a, 1, &vmax);
    if(nb < 0){
      nb = usolve_bal_test(B, a, 0, &vmax);
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

  fmpz_t c2, cr;
  fmpz_init(c2);
  fmpz_init(cr);
  int spawn = 0;
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

  /* at most vmax roots: vmax cells of a grid with a sign change each
   * contain one root each, and there is no other one */
  if(vmax >= 2 && vmax <= USOLVE_BAL_GRIDMAX
     && usolve_bal_grid(B, c, k, vmax)){
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

  /* in usolve_bal_run, the right half of a node with at least 3 sign
   * variations, always investigated, is given to a task, which starts from
   * the state of B at this node: the result does not depend on the number
   * of threads */
  if(B->top != NULL && nb >= 3){
    spawn = 1;
    B->tasks++;
    fmpz_mul_2exp(cr, c, 1);
    fmpz_add_ui(cr, cr, 1);
    usolve_bal_t P = *B;
    usolve_flags pf = *flags;
#pragma omp task firstprivate(P, pf, a) shared(cr)
    {
      usolve_flags fl;
      usolve_bal_t T;
      usolve_bal_task_init(&T, &P, &pf, &fl);
      usolve_bal_right(&T, a, cr, k);
      usolve_bal_task_done(&T);
    }
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

  if(spawn){
    goto clear;
  }

  if(nb <= 2 && nbleft == nb){
    (flags->skipped)++;
    if(flags->verbose >= 1){
      fprintf(VERBSTREAM, "s");
    }
    goto clear;
  }

  usolve_bal_right(B, a, c2, k);
clear:
  if(spawn){
    /* the task reads a and cr */
#pragma omp taskwait
  }
  fmpz_clear(c2);
  fmpz_clear(cr);
  if(owned){
    usolve_bal_node_clear(&own, n);
  }
  return nb;
}

static void usolve_bal_init(usolve_bal_t *B, mpz_t *upol0,
                            const unsigned long int deg0, const int o,
                            const int s, const long b, const long z,
                            interval *roots, unsigned long int *nbr,
                            usolve_flags *flags){
  const slong n = deg0 - o;
  B->n = n;
  B->fl = deg0 + 1;
  B->f = malloc(B->fl * sizeof(fmpz));
  for(slong i = 0; i < B->fl; i++){
    fmpz_init_set_readonly(B->f + i, upol0[i]);
  }
  B->o = o;
  B->s = s;
  /* see USOLVEmpz_poly_rescale_normalize_2exp_th */
  B->b = b;
  B->d = ((b > 0) ? 0 : -b * n) - z;
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
  B->powers = malloc(sizeof(usolve_bal_cache_t));
  B->powers->p = malloc(USOLVE_BAL_NCACHE * sizeof(usolve_bal_power_t));
  B->powers->n = 0;
  B->queue = NULL;
  B->qc = NULL;
  B->qk = NULL;
  B->nqueue = 0;
  B->aqueue = 0;
  B->smax = 0;
  B->roots = roots;
  B->nbr = nbr;
  B->flags = flags;
  B->stop = 0;
  B->fallbacks = 0;
  B->subtrees = 0;
  B->sign_extra = 64;
  B->grids = 0;
  B->grids_ok = 0;
  B->grid_signs = 0;
  B->top = NULL;
  B->tasks = 0;
}

static void usolve_bal_clear(usolve_bal_t *B){
  if(B->flags->verbose >= 1){
    fprintf(VERBSTREAM, "\n%ld subtrees in fixed point (%ld nodes), %ld "
            "exact node polynomials computed again, final precision %ld "
            "bits; %ld of %ld grids conclusive (%ld signs); %ld tasks\n",
            B->subtrees, B->nodes, B->fallbacks, (long)B->w, B->grids_ok,
            B->grids, B->grid_signs, B->tasks);
  }
  for(long i = 0; i < B->powers->n; i++){
    fmpz_poly_clear(B->powers->p[i].c);
    fmpz_poly_clear(B->powers->p[i].ce);
  }
  free(B->powers->p);
  free(B->powers);
  for(long i = 0; i < B->nqueue; i++){
    free(B->queue[i]);
  }
  free(B->queue);
  for(long i = 0; i < B->aqueue; i++){
    fmpz_clear(B->qc + i);
  }
  free(B->qc);
  free(B->qk);
  free(B->lbin);
  for(slong i = 0; i < B->fl; i++){
    fmpz_clear_readonly(B->f + i);
  }
  free(B->f);
}

/* the node a in one block of words: s, exact, then for each coefficient
 * of v and e, its number of limbs l and sign (2 l + 1 if negative), then
 * its limbs. A queued node is kept this way while the exact bisection
 * goes on: its integers, scattered in the heap among the large ones of the
 * exact bisection, would keep the memory these free from being reused
 * (and those of a rebalanced node keep the room of the exact
 * coefficients) */
static ulong *usolve_bal_pack(const usolve_bal_node_t *a, const slong n){
  slong len = 2;
  for(slong i = 0; i <= n; i++){
    len += 2 + fmpz_size(a->v + i) + fmpz_size(a->e + i);
  }
  ulong *buf = malloc(len * sizeof(ulong));
  buf[0] = (ulong)a->s;
  buf[1] = (ulong)a->exact;
  slong pos = 2;
  fmpz_t x;
  fmpz_init(x);
  for(slong i = 0; i < 2 * (n + 1); i++){
    const fmpz *y = (i <= n) ? a->v + i : a->e + (i - n - 1);
    const slong l = fmpz_size(y);
    buf[pos++] = 2 * (ulong)l + (fmpz_sgn(y) < 0);
    if(l > 0){
      fmpz_abs(x, y);
      fmpz_get_ui_array(buf + pos, l, x);
      pos += l;
    }
  }
  fmpz_clear(x);
  return buf;
}

static void usolve_bal_unpack(usolve_bal_node_t *a, const ulong *buf,
                              const slong n){
  a->s = (slong)buf[0];
  a->exact = (int)buf[1];
  slong pos = 2;
  for(slong i = 0; i < 2 * (n + 1); i++){
    fmpz *y = (i <= n) ? a->v + i : a->e + (i - n - 1);
    const slong l = buf[pos] / 2;
    const int neg = buf[pos] & 1;
    pos++;
    if(l == 0){
      fmpz_zero(y);
    }
    else{
      fmpz_set_ui_array(y, buf + pos, l);
      if(neg){
        fmpz_neg(y, y);
      }
      pos += l;
    }
  }
}

/* called by bisection_rec on the node (c, k), whose exact polynomial is
 * upol, of degree deg: when the node polynomial still has the degree of
 * F (no exact root was divided out), upol(0) != 0, and its best slope is
 * at least USOLVE_BAL_SLOPE, its subtree is investigated in fixed point,
 * later by usolve_bal_run when the test of the node shows at least two
 * sign variations, at once otherwise; then 1 is returned, *nbv is the
 * number of sign variations of the node, and upol is unchanged.
 * Otherwise, 0 is returned. */
static int usolve_bal_subtree(usolve_bal_t *B, mpz_t *upol,
                              const unsigned long int deg, const mpz_t c,
                              const long k, long *nbv){
  const slong n = B->n;
  double sz = 0;
  for(unsigned long int i = 0; i <= deg; i++){
    sz += (double)mpz_size(upol[i]) * sizeof(mp_limb_t);
  }
  B->smax = FLINT_MAX(B->smax, sz);
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
  /* a subtree with at least 2 sign variations is left for usolve_bal_run,
   * unless the bisection stops at the first root */
  long vmax;
  const long nb = B->flags->hasrealroots ? -1 : usolve_bal_test(B, &a, 1, &vmax);
  if(nb >= 2){
    if(B->nqueue == B->aqueue){
      const long na = FLINT_MAX(16, 2 * B->aqueue);
      B->queue = realloc(B->queue, na * sizeof(ulong *));
      B->qc = realloc(B->qc, na * sizeof(fmpz));
      for(long i = B->aqueue; i < na; i++){
        fmpz_init(B->qc + i);
      }
      B->qk = realloc(B->qk, na * sizeof(slong));
      B->aqueue = na;
    }
    B->queue[B->nqueue] = usolve_bal_pack(&a, n);
    usolve_bal_node_clear(&a, n);
    fmpz_set_mpz(B->qc + B->nqueue, c);
    B->qk[B->nqueue] = k;
    B->nqueue++;
    *nbv = nb;
    return 1;
  }
  fmpz_t cf;
  fmpz_init(cf);
  fmpz_set_mpz(cf, c);
  *nbv = usolve_bal_node(B, &a, cf, k);
  fmpz_clear(cf);
  usolve_bal_node_clear(&a, n);
  return 1;
}

/* peak and current resident memory of the process, in bytes; returns 0
 * if the current one cannot be read (always on systems other than Linux
 * and macOS) */
static int usolve_bal_rss(double *peak, double *cur){
#if defined(__APPLE__) || defined(__linux__)
  struct rusage ru;
  if(getrusage(RUSAGE_SELF, &ru)){
    return 0;
  }
#endif
#if defined(__APPLE__)
  *peak = (double)ru.ru_maxrss;
  struct mach_task_basic_info info;
  mach_msg_type_number_t cnt = MACH_TASK_BASIC_INFO_COUNT;
  if(task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
               (task_info_t)&info, &cnt) != KERN_SUCCESS){
    return 0;
  }
  *cur = (double)info.resident_size;
  return 1;
#elif defined(__linux__)
  *peak = 1024.0 * ru.ru_maxrss;
  long pages, res;
  FILE *fd = fopen("/proc/self/statm", "r");
  if(fd == NULL){
    return 0;
  }
  const int ok = (fscanf(fd, "%ld %ld", &pages, &res) == 2);
  fclose(fd);
  *cur = (double)res * sysconf(_SC_PAGESIZE);
  return ok;
#else
  (void)peak;
  (void)cur;
  return 0;
#endif
}

/* order of the intervals of the roots: by left ends (numer / 2^k, or
 * numer when k <= 0, see merge_root), an exact root before the interval
 * starting at it */
static int usolve_bal_cmp_roots(const void *x, const void *y){
  const interval *a = x, *b = y;
  const long ka = FLINT_MAX(a->k, 0), kb = FLINT_MAX(b->k, 0);
  mpz_t s, t;
  mpz_init(s);
  mpz_init(t);
  if(ka >= kb){
    mpz_set(s, a->numer);
    mpz_mul_2exp(t, b->numer, ka - kb);
  }
  else{
    mpz_mul_2exp(s, a->numer, kb - ka);
    mpz_set(t, b->numer);
  }
  int r = mpz_cmp(s, t);
  mpz_clear(s);
  mpz_clear(t);
  if(r == 0){
    r = (int)b->isexact - (int)a->isexact;
  }
  return (r > 0) - (r < 0);
}

/* investigates the subtrees left by usolve_bal_subtree with the threads of
 * B->flags, as many as fit below the peak memory of the process so far
 * (that of the exact bisection, or of earlier computations, such as the
 * parametrization in msolve), or, when the memory of the process cannot be
 * read, below the memory used by the exact bisection, and at least in
 * USOLVE_BAL_INMEM (n + 1) tau bits, tau being the largest bit size of the
 * coefficients of F: since the memory of a thread grows as n^2, this
 * allows about USOLVE_BAL_INMEM tau / (1.5 USOLVE_BAL_THRMEM n) threads
 * whatever the memory used before. Each one starts from
 * the state of B (working precision, ...) and learns on its own, as does
 * each task of usolve_bal_node, so that the output does not depend on the
 * number of threads; the roots are then sorted back in the order of the
 * sequential bisection (by increasing x for the positive roots, by
 * decreasing x for the negative ones) */
static void usolve_bal_run(usolve_bal_t *B){
  const long nq = B->nqueue;
  if(nq == 0){
    return;
  }
  const slong n = B->n;
  usolve_flags *flags = B->flags;
  const slong w0 = B->w;
  /* the work space of FLINT's products in the exact bisection */
  flint_cleanup();
#if defined(__GLIBC__)
  /* the memory freed so far is given back to the system, so that the
   * resident memory measures what is used */
  malloc_trim(0);
#endif
  const double thr = USOLVE_BAL_THRMEM * (double)(n + 1) * w0 / 8;
  double peak, cur, room;
  if(usolve_bal_rss(&peak, &cur)){
    room = FLINT_MAX(peak - cur, 0);
  }
  else{
    room = USOLVE_BAL_TOPMEM * B->smax + thr;
  }
  const double tau = FLINT_ABS(_fmpz_vec_max_bits(B->f, B->fl));
  room = FLINT_MAX(room, USOLVE_BAL_INMEM * (double)(n + 1) * tau / 8);
  /* the main thread counts too: its work space was freed above */
  int nthreads = (int)FLINT_MIN((double)flags->nthreads, room / thr);
  nthreads = FLINT_MAX(nthreads, 1);
  if(flags->verbose >= 1){
    fprintf(VERBSTREAM, "\n%ld subtrees in fixed point on %d threads "
            "(%.0f MB available, about %.0f MB per thread)\n", nq, nthreads,
            room / 1e6, thr / 1e6);
  }
  /* the state each subtree starts from; B only receives the counters */
  usolve_bal_t S = *B;
  S.top = B;
  S.recent_nodes = 0;
  S.recent_fallbacks = 0;
  const usolve_flags sf = *flags;
#pragma omp parallel num_threads(nthreads)
  {
#pragma omp for schedule(dynamic, 1)
    for(long j = 0; j < nq; j++){
      usolve_flags fl;
      usolve_bal_t T;
      usolve_bal_task_init(&T, &S, &sf, &fl);
      usolve_bal_node_t a;
      usolve_bal_node_init(&a, n);
      usolve_bal_unpack(&a, B->queue[j], n);
      free(B->queue[j]);
      B->queue[j] = NULL;
      usolve_bal_node(&T, &a, B->qc + j, B->qk[j]);
      usolve_bal_node_clear(&a, n);
      usolve_bal_task_done(&T);
    }
    flint_cleanup();
  }
  B->nqueue = 0;
  const unsigned long int nr = *(B->nbr);
  qsort(B->roots, nr, sizeof(interval), usolve_bal_cmp_roots);
  if(flags->sign){
    for(unsigned long int i = 0; i < nr / 2; i++){
      interval x = B->roots[i];
      B->roots[i] = B->roots[nr - 1 - i];
      B->roots[nr - 1 - i] = x;
    }
  }
}
