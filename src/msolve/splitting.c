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
 * Vincent Neiger
 * Mohab Safey El Din */

/* Splitting the solution set along the equalities of coordinates.

   On systems with symmetries, many solutions have equal coordinates
   x_i = x_j, and the eliminating polynomial is the product of many
   factors of small height: its height is mostly the one of this product.
   For each pair (i, j), gcd(elim, N_i - N_j) modulo a prime is the factor
   of elim whose roots are the solutions with x_i = x_j. Refining elim
   along all pairs gives groups of solutions with the same equalities.
   These groups are defined over QQ, hence the same at all good primes
   (they are matched by their sets of equalities), and each is lifted on
   its own. The parametrization of the whole solution set is then
   recovered exactly over QQ, without further primes.

   The parametrizations modulo a prime are normalized (see
   normalize_nmod_param): elim is monic and squarefree, and
   coords[v] = -x_v denom mod elim with denom = elim' / deg(elim).
   A group with monic eliminating polynomial g is described by g and by
   B_v = x_v g' mod g for each coordinate. Over QQ, with E the product of
   the g's and H_g = E / g, the numerators of the parametrization of E are
   R_v = -x_v E' mod E = -sum_g B_v H_g, since E' = g' H_g mod g. */

#ifdef HAVE_FLINT_NMOD_H
#include <flint/nmod.h>
#endif
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_vec.h>

/* The splitting is used, unless disabled by --split 1, when the largest
 * group has at most 1 / SPLIT_MIN_RATIO of the solutions, or whenever there
 * are several groups with --split 2: otherwise, the lifting
 * of the largest group costs about as much as the usual lifting, which has
 * better reconstruction bounds. It is abandoned when the groups differ
 * modulo SPLIT_MAX_MISMATCH primes. */
#define SPLIT_MIN_RATIO 4
#define SPLIT_MAX_MISMATCH 3

typedef struct {
  long deg;
  uint64_t *sig;      /* pairs (i, j) with x_i = x_j on the group, bitset */
  nmod_poly_t g;      /* monic eliminating polynomial of the group */
} split_part_t;

typedef struct {
  long deg;
  uint64_t *sig;
  long n;             /* number of lifted values: deg * (nc + 1) */
  fmpz *crt;          /* g (deg coefficients, monic), then the B_v */
  fmpz *num;          /* reconstructed values, num[i] / den[i / deg]: a
                         common denominator for g and for each B_v, which
                         one rational number per value would repeat */
  fmpz *den;
  int state;          /* 0: not reconstructed, 1: reconstructed, 2: checked */
} split_group_t;

typedef struct {
  long nc;            /* number of coordinates */
  long nwords;        /* words of a bitset of pairs */
  long ngroups;
  split_group_t *groups;
  fmpz_t modulus;
  long nprimes;
  long ndone;
} split_data_t;

static void split_parts_clear(split_part_t *parts, long n) {
  for (long i = 0; i < n; i++) {
    free(parts[i].sig);
    nmod_poly_clear(parts[i].g);
  }
  free(parts);
}

static int split_cmp_sig(const uint64_t *a, const uint64_t *b, long nwords) {
  for (long w = 0; w < nwords; w++) {
    if (a[w] != b[w]) {
      return a[w] < b[w] ? -1 : 1;
    }
  }
  return 0;
}

static long split_cmp_nwords;

static int split_cmp_parts(const void *a, const void *b) {
  const split_part_t *pa = a, *pb = b;
  return split_cmp_sig(pa->sig, pb->sig, split_cmp_nwords);
}

/* refines the eliminating polynomial of param along the equalities of
 * coordinates; returns the number of groups, sorted by their sets of
 * equalities, in *parts_ptr */
static long split_mod(split_part_t **parts_ptr, const param_t *param,
                      const long nc, const long nwords) {
  const mp_limb_t p = param->charac;
  long alloc = 16, n = 1;
  split_part_t *parts = malloc(alloc * sizeof(split_part_t));
  parts[0].sig = calloc(nwords, sizeof(uint64_t));
  nmod_poly_init(parts[0].g, p);
  nmod_poly_make_monic(parts[0].g, param->elim);
  parts[0].deg = nmod_poly_degree(parts[0].g);

  nmod_poly_t d, r, h;
  nmod_poly_init(d, p);
  nmod_poly_init(r, p);
  nmod_poly_init(h, p);
  long k = 0;
  for (long i = 0; i < nc; i++) {
    for (long j = i + 1; j < nc; j++, k++) {
      nmod_poly_sub(d, param->coords[i], param->coords[j]);
      const long n0 = n;
      for (long m = 0; m < n0; m++) {
        nmod_poly_rem(r, d, parts[m].g);
        if (nmod_poly_is_zero(r)) {
          parts[m].sig[k / 64] |= UINT64_C(1) << (k % 64);
          continue;
        }
        nmod_poly_gcd(h, parts[m].g, r);
        if (nmod_poly_degree(h) <= 0) {
          continue;
        }
        if (n == alloc) {
          alloc *= 2;
          parts = realloc(parts, alloc * sizeof(split_part_t));
        }
        parts[n].sig = malloc(nwords * sizeof(uint64_t));
        memcpy(parts[n].sig, parts[m].sig, nwords * sizeof(uint64_t));
        parts[n].sig[k / 64] |= UINT64_C(1) << (k % 64);
        nmod_poly_init(parts[n].g, p);
        nmod_poly_set(parts[n].g, h);
        parts[n].deg = nmod_poly_degree(h);
        nmod_poly_div(parts[m].g, parts[m].g, h);
        parts[m].deg = nmod_poly_degree(parts[m].g);
        n++;
      }
    }
  }
  nmod_poly_clear(d);
  nmod_poly_clear(r);
  nmod_poly_clear(h);
  split_cmp_nwords = nwords;
  qsort(parts, n, sizeof(split_part_t), split_cmp_parts);
  *parts_ptr = parts;
  return n;
}

/* values of a group modulo p: the coefficients of g but the leading one,
 * then x_v g' mod g for each coordinate */
static void split_part_values(mp_limb_t *val, const split_part_t *part,
                              const param_t *param, const long nc) {
  const mp_limb_t p = param->charac;
  const long deg = part->deg;
  nmod_poly_t inv, dg, x;
  nmod_poly_init(inv, p);
  nmod_poly_init(dg, p);
  nmod_poly_init(x, p);
  for (long i = 0; i < deg; i++) {
    val[i] = nmod_poly_get_coeff_ui(part->g, i);
  }
  nmod_poly_rem(inv, param->denom, part->g);
  nmod_poly_invmod(inv, inv, part->g);
  nmod_poly_neg(inv, inv);
  nmod_poly_derivative(dg, part->g);
  nmod_poly_mulmod(inv, inv, dg, part->g);
  for (long v = 0; v < nc; v++) {
    nmod_poly_rem(x, param->coords[v], part->g);
    nmod_poly_mulmod(x, x, inv, part->g);
    for (long i = 0; i < deg; i++) {
      val[(v + 1) * deg + i] = nmod_poly_get_coeff_ui(x, i);
    }
  }
  nmod_poly_clear(inv);
  nmod_poly_clear(dg);
  nmod_poly_clear(x);
}

static void split_data_clear(split_data_t *sp) {
  for (long i = 0; i < sp->ngroups; i++) {
    split_group_t *G = sp->groups + i;
    free(G->sig);
    _fmpz_vec_clear(G->crt, G->n);
    _fmpz_vec_clear(G->num, G->n);
    _fmpz_vec_clear(G->den, sp->nc + 1);
  }
  free(sp->groups);
  fmpz_clear(sp->modulus);
}

/* initializes the groups from the parametrization at the first prime;
 * returns 0 if there is a single group */
static int split_data_init(split_data_t *sp, const param_t *param) {
  sp->nc = param->nvars - 1;
  const long npairs = sp->nc * (sp->nc - 1) / 2;
  sp->nwords = (npairs + 63) / 64;
  split_part_t *parts;
  const long n = split_mod(&parts, param, sp->nc, sp->nwords);
  sp->ngroups = n;
  fmpz_init(sp->modulus);
  sp->groups = NULL;
  if (n < 2) {
    split_parts_clear(parts, n);
    sp->ngroups = 0;
    return 0;
  }
  sp->groups = malloc(n * sizeof(split_group_t));
  for (long i = 0; i < n; i++) {
    split_group_t *G = sp->groups + i;
    G->deg = parts[i].deg;
    G->sig = parts[i].sig;
    parts[i].sig = NULL;
    G->n = G->deg * (sp->nc + 1);
    G->crt = _fmpz_vec_init(G->n);
    G->num = _fmpz_vec_init(G->n);
    G->den = _fmpz_vec_init(sp->nc + 1);
    G->state = 0;
  }
  split_parts_clear(parts, n);
  fmpz_one(sp->modulus);
  sp->nprimes = 0;
  sp->ndone = 0;
  return 1;
}

static int split_check_rec(const split_group_t *G, const mp_limb_t *val,
                           mp_limb_t p) {
  nmod_t mod;
  nmod_init(&mod, p);
  mp_limb_t id = 0;
  for (long i = 0; i < G->n; i++) {
    if (i % G->deg == 0) {
      const mp_limb_t den = fmpz_fdiv_ui(G->den + i / G->deg, p);
      if (den == 0) {
        return 0;
      }
      id = nmod_inv(den, mod);
    }
    if (nmod_mul(fmpz_fdiv_ui(G->num + i, p), id, mod) != val[i]) {
      return 0;
    }
  }
  return 1;
}

/* num[i] / den[b] = r for the value i of the block b, the numerators of
 * the block before i being brought to the new common denominator */
static void split_set_rec(split_group_t *G, const long i, const fmpq_t r,
                          fmpz_t t) {
  const long b = i / G->deg, i0 = b * G->deg;
  fmpz *d = G->den + b;
  if (i == i0) {
    fmpz_one(d);
  }
  if (!fmpz_divisible(d, fmpq_denref(r))) {
    fmpz_lcm(t, d, fmpq_denref(r));
    fmpz_divexact(d, t, d);
    _fmpz_vec_scalar_mul_fmpz(G->num + i0, G->num + i0, i - i0, d);
    fmpz_swap(d, t);
  }
  fmpz_divexact(t, d, fmpq_denref(r));
  fmpz_mul(G->num + i, fmpq_numref(r), t);
}

/* adds the image of the parametrization modulo a new prime; returns 0 if
 * its groups do not match those at the first prime */
static int split_add_image(split_data_t *sp, const param_t *param) {
  const mp_limb_t p = param->charac;
  split_part_t *parts;
  const long n = split_mod(&parts, param, sp->nc, sp->nwords);
  int match = n == sp->ngroups;
  for (long i = 0; i < n && match; i++) {
    match = parts[i].deg == sp->groups[i].deg
      && split_cmp_sig(parts[i].sig, sp->groups[i].sig, sp->nwords) == 0;
  }
  if (!match) {
    split_parts_clear(parts, n);
    return 0;
  }
  for (long i = 0; i < n; i++) {
    split_group_t *G = sp->groups + i;
    if (G->state == 2) {
      continue;
    }
    mp_limb_t *val = malloc(G->n * sizeof(mp_limb_t));
    split_part_values(val, parts + i, param, sp->nc);
    if (G->state == 1) {
      if (split_check_rec(G, val, p)) {
        G->state = 2;
        sp->ndone++;
        free(val);
        continue;
      }
      G->state = 0;
    }
    for (long k = 0; k < G->n; k++) {
      if (fmpz_is_one(sp->modulus)) {
        fmpz_set_ui(G->crt + k, val[k]);
      } else {
        fmpz_CRT_ui(G->crt + k, G->crt + k, sp->modulus, val[k], p, 0);
      }
    }
    free(val);
  }
  split_parts_clear(parts, n);
  fmpz_mul_ui(sp->modulus, sp->modulus, p);
  sp->nprimes++;
  return 1;
}

/* rational reconstruction of the groups which are not reconstructed yet */
static void split_reconstruct(split_data_t *sp) {
  for (long i = 0; i < sp->ngroups; i++) {
    split_group_t *G = sp->groups + i;
    if (G->state != 0) {
      continue;
    }
    fmpq_t r;
    fmpz_t t;
    fmpq_init(r);
    fmpz_init(t);
    long k = 0;
    while (k < G->n && fmpq_reconstruct_fmpz(r, G->crt + k, sp->modulus)) {
      split_set_rec(G, k, r, t);
      k++;
    }
    fmpq_clear(r);
    fmpz_clear(t);
    if (k == G->n) {
      G->state = 1;
    }
  }
}

/* splits the solutions from the parametrization at the first prime and
 * adds its image; returns 0 (sp is then cleared) if the splitting is not
 * used; mode is that of --split (0 chosen, 1 never, 2 always) */
static int split_start(split_data_t *sp, const param_t *param,
                       const int mode, const int info_level) {
  if (mode == 1) {
    return 0;
  }
  double rt = realtime();
  if (!split_data_init(sp, param)) {
    split_data_clear(sp);
    return 0;
  }
  long maxdeg = 0;
  for (long i = 0; i < sp->ngroups; i++) {
    maxdeg = MAX(maxdeg, sp->groups[i].deg);
  }
  const long nsols = nmod_poly_degree(param->elim);
  const int used = mode == 2 || SPLIT_MIN_RATIO * maxdeg <= nsols;
  if (info_level) {
    fprintf(VERBSTREAM, "\nSplitting along equalities of coordinates: "
            "%ld components, largest of degree %ld out of %ld%s (%.2f sec)\n",
            sp->ngroups, maxdeg, nsols, used ? "" : ", not used",
            realtime() - rt);
  }
  if (!used) {
    split_data_clear(sp);
    return 0;
  }
  split_add_image(sp, param);
  split_reconstruct(sp);
  return 1;
}

/* P = product of the g's of groups lo..hi-1, S[v] = sum of B_v P / g;
 * the reconstructed values are freed as they are read if free_rec is set,
 * so that they are not held together with the polynomials (random2_10_10:
 * 547 MB at the end instead of 473 MB) */
static void split_combine(fmpq_poly_t P, fmpq_poly_struct *S,
                          const split_data_t *sp, long lo, long hi,
                          const int free_rec) {
  const long nc = sp->nc;
  if (hi - lo == 1) {
    const split_group_t *G = sp->groups + lo;
    const long deg = G->deg;
    /* g = (sum num_i t^i + den t^deg) / den */
    fmpz_poly_t q;
    fmpz_poly_init2(q, deg + 1);
    for (long b = 0; b <= nc; b++) {
      fmpz_poly_zero(q);
      for (long i = 0; i < deg; i++) {
        fmpz_poly_set_coeff_fmpz(q, i, G->num + b * deg + i);
      }
      if (b == 0) {
        fmpz_poly_set_coeff_fmpz(q, deg, G->den);
      }
      fmpq_poly_struct *R = b == 0 ? P : S + b - 1;
      fmpq_poly_set_fmpz_poly(R, q);
      fmpq_poly_scalar_div_fmpz(R, R, G->den + b);
      if (free_rec) {
        for (long i = 0; i < deg; i++) {
          fmpz_zero(G->num + b * deg + i);
        }
      }
    }
    fmpz_poly_clear(q);
    return;
  }
  const long mid = (lo + hi) / 2;
  fmpq_poly_t PR, tmp;
  fmpq_poly_init(PR);
  fmpq_poly_init(tmp);
  fmpq_poly_struct *SR = malloc(nc * sizeof(fmpq_poly_struct));
  for (long v = 0; v < nc; v++) {
    fmpq_poly_init(SR + v);
  }
  split_combine(P, S, sp, lo, mid, free_rec);
  split_combine(PR, SR, sp, mid, hi, free_rec);
  for (long v = 0; v < nc; v++) {
    fmpq_poly_mul(tmp, SR + v, P);
    fmpq_poly_mul(S + v, S + v, PR);
    fmpq_poly_add(S + v, S + v, tmp);
    fmpq_poly_clear(SR + v);
  }
  fmpq_poly_mul(P, P, PR);
  free(SR);
  fmpq_poly_clear(PR);
  fmpq_poly_clear(tmp);
}

/* out = in, padded with zeros to at least minlen coefficients */
static void split_set_mpz_upoly(mpz_upoly_t out, const fmpz_poly_t in,
                                const long minlen) {
  const long len = MAX(fmpz_poly_length(in), minlen);
  if (out->alloc < len) {
    out->coeffs = realloc(out->coeffs, len * sizeof(mpz_t));
    for (long i = out->alloc; i < len; i++) {
      mpz_init(out->coeffs[i]);
    }
    out->alloc = len;
  }
  for (long i = 0; i < len; i++) {
    if (i < fmpz_poly_length(in)) {
      fmpz_get_mpz(out->coeffs[i], in->coeffs + i);
    } else {
      mpz_set_ui(out->coeffs[i], 0);
    }
  }
  out->length = len;
}

/* the parametrization of the solutions of groups lo..hi-1, with its
 * denominator elim'; their reconstructed values are freed if free_rec is
 * set */
static void split_param_range(mpz_param_t mpz_param, const split_data_t *sp,
                              long lo, long hi, const int free_rec) {
  const long nc = sp->nc;
  fmpq_poly_t E;
  fmpq_poly_init(E);
  fmpq_poly_struct *S = malloc(nc * sizeof(fmpq_poly_struct));
  for (long v = 0; v < nc; v++) {
    fmpq_poly_init(S + v);
  }
  split_combine(E, S, sp, lo, hi, free_rec);

  /* integer primitive eliminating polynomial lambda E, with lambda > 0 */
  fmpz_poly_t num;
  fmpz_poly_init(num);
  fmpq_poly_get_numerator(num, E);
  fmpz_t ct;
  fmpz_init(ct);
  fmpz_poly_content(ct, num);
  fmpz_poly_scalar_divexact_fmpz(num, num, ct);
  split_set_mpz_upoly(mpz_param->elim, num, 0);
  const long deg = fmpz_poly_degree(num);
  fmpq_t lambda;
  fmpq_init(lambda);
  fmpq_set_fmpz_frac(lambda, fmpq_poly_denref(E), ct);
  fmpq_neg(lambda, lambda);
  fmpz_poly_derivative(num, num);
  split_set_mpz_upoly(mpz_param->denom, num, 0);

  /* R_v = -lambda sum_g B_v H_g = coords[v] / cfs[v] */
  for (long v = 0; v < nc; v++) {
    fmpq_poly_scalar_mul_fmpq(S + v, S + v, lambda);
    fmpq_poly_get_numerator(num, S + v);
    /* as in the usual output, coords[v] has deg(elim) coefficients */
    split_set_mpz_upoly(mpz_param->coords[v], num, deg);
    fmpz_get_mpz(mpz_param->cfs[v], fmpq_poly_denref(S + v));
    fmpq_poly_clear(S + v);
  }
  free(S);
  fmpq_clear(lambda);
  fmpz_clear(ct);
  fmpz_poly_clear(num);
  fmpq_poly_clear(E);
}

/* the parametrization of the whole solution set, from the groups, whose
 * reconstructed values are freed if free_rec is set */
static void split_recombine(mpz_param_t mpz_param, const split_data_t *sp,
                            const int free_rec) {
  split_param_range(mpz_param, sp, 0, sp->ngroups, free_rec);
}

static void split_params_clear(mpz_param_array_struct *params) {
  for (long i = 0; i < params->nb; i++) {
    mpz_param_clear(params->params[i]);
  }
  free(params->params);
  params->params = NULL;
  params->nb = 0;
}

/* the parametrizations of the groups, one by one, whose reconstructed
 * values are freed if free_rec is set */
static void split_set_params(mpz_param_array_struct *params,
                             const split_data_t *sp, const int free_rec) {
  const long nc = sp->nc;
  split_params_clear(params);
  params->nb = sp->ngroups;
  params->params = malloc(sp->ngroups * sizeof(mpz_param_t));
  for (long i = 0; i < sp->ngroups; i++) {
    mpz_param_struct *P = params->params[i];
    mpz_param_init(P);
    P->nvars = nc + 1;
    P->nsols = sp->groups[i].deg;
    P->dquot = sp->groups[i].deg;
    P->dim = 0;
    P->coords = malloc(nc * sizeof(mpz_upoly_t));
    P->cfs = malloc(nc * sizeof(mpz_t));
    for (long v = 0; v < nc; v++) {
      P->coords[v]->coeffs = NULL;
      P->coords[v]->alloc = 0;
      P->coords[v]->length = 0;
      mpz_init(P->cfs[v]);
    }
    split_param_range(P, sp, i, i + 1, free_rec);
  }
}
