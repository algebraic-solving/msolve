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

/* Search for a separating element of a zero-dimensional ideal from one
 * reduced DRL Groebner basis modulo a prime.
 *
 * The multiplication matrices by all variables are built from the normal
 * forms of the border monomials, computed by increasing DRL order. The
 * degree of the minimal polynomial of each multiplication matrix is
 * computed with a Wiedemann sequence and Berlekamp-Massey. The target is
 * the degree dgen of the minimal polynomial of a random linear form
 * (dquot for a radical ideal): a linear form reaching it separates the
 * solutions, and has at each of them the largest nilpotency index, which
 * sep_dense_param needs. A variable reaching dgen is separating. If there
 * is none, a linear form sum c_i x_i is built greedily, taking the
 * variables by decreasing degree of their minimal polynomials, and the
 * smallest variable first for equal degrees. For a radical ideal, this
 * degree is the number of values taken on the solutions. The first
 * variable gets coefficient 1. Each next variable x_v gets the smallest
 * c_v >= 1 such that the form reaches the degree of the pair (form, x_v),
 * which is the degree for a coefficient of x_v chosen at random, or 0 if
 * x_v does not increase the degree. For a radical ideal, a coefficient
 * fails only if it merges two solutions separated by (form, x_v), so at
 * most dquot (dquot - 1) / 2 of them fail, and the degree is dquot after
 * the last variable. The construction stops once the degree is dgen, the
 * remaining coefficients being 0. As a safeguard, c_v is at most dquot;
 * if no c_v reaches the target, the best one is kept. The search is
 * deterministic: it only depends on the Groebner basis and the prime. */

#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>


/* hash table of exponent vectors */
typedef struct {
  int32_t nv;
  long nelts;
  long size;      /* power of 2 */
  long *tab;      /* -1 if empty, index in exps otherwise */
  int32_t *exps;  /* nelts * nv, not owned */
} sep_hash_t;

static inline uint64_t sep_hash_exp(const int32_t *e, int32_t nv) {
  uint64_t h = 1469598103934665603ULL;
  for (int32_t k = 0; k < nv; k++) {
    h ^= (uint64_t)(uint32_t)e[k];
    h *= 1099511628211ULL;
  }
  return h ^ (h >> 29);
}

static void sep_hash_init(sep_hash_t *H, int32_t *exps, long nelts,
                          int32_t nv) {
  H->nv = nv;
  H->nelts = nelts;
  H->exps = exps;
  H->size = 16;
  while (H->size < 2 * nelts + 16) {
    H->size *= 2;
  }
  H->tab = malloc(H->size * sizeof(long));
  for (long i = 0; i < H->size; i++) {
    H->tab[i] = -1;
  }
  for (long i = 0; i < nelts; i++) {
    uint64_t h = sep_hash_exp(exps + i * nv, nv) & (H->size - 1);
    while (H->tab[h] != -1) {
      h = (h + 1) & (H->size - 1);
    }
    H->tab[h] = i;
  }
}

static long sep_hash_find(const sep_hash_t *H, const int32_t *e) {
  uint64_t h = sep_hash_exp(e, H->nv) & (H->size - 1);
  while (H->tab[h] != -1) {
    if (memcmp(H->exps + H->tab[h] * H->nv, e, H->nv * sizeof(int32_t)) == 0) {
      return H->tab[h];
    }
    h = (h + 1) & (H->size - 1);
  }
  return -1;
}

static int32_t sep_cmp_nv;

/* increasing DRL order, the first variable being the largest one */
static int sep_cmp_drl(const void *a, const void *b) {
  const int32_t *ea = a, *eb = b;
  int32_t da = 0, db = 0;
  for (int32_t k = 0; k < sep_cmp_nv; k++) {
    da += ea[k];
    db += eb[k];
  }
  if (da != db) {
    return da < db ? -1 : 1;
  }
  for (int32_t k = sep_cmp_nv - 1; k >= 0; k--) {
    if (ea[k] != eb[k]) {
      return ea[k] > eb[k] ? -1 : 1;
    }
  }
  return 0;
}

/* degree of the minimal polynomial of the dense dquot x dquot matrix mat
 * (row major), computed from a Wiedemann sequence of length 2 dquot */
static long sep_minpoly_degree(const mp_limb_t *mat, long D, const mp_limb_t *u,
                               const mp_limb_t *r, nmod_t mod,
                               mp_limb_t *x, mp_limb_t *y, mp_limb_t *seq,
                               mp_limb_t *C, mp_limb_t *B, mp_limb_t *T) {
  const long len = 2 * D;
  memcpy(x, r, D * sizeof(mp_limb_t));
  for (long k = 0; k < len; k++) {
    unsigned __int128 acc = 0;
    for (long i = 0; i < D; i++) {
      acc += (unsigned __int128)u[i] * x[i];
    }
    seq[k] = (mp_limb_t)(acc % mod.n);
    if (k == len - 1) {
      break;
    }
    for (long i = 0; i < D; i++) {
      const mp_limb_t *row = mat + i * D;
      unsigned __int128 a = 0;
      for (long j = 0; j < D; j++) {
        a += (unsigned __int128)row[j] * x[j];
      }
      y[i] = (mp_limb_t)(a % mod.n);
    }
    mp_limb_t *tmp = x;
    x = y;
    y = tmp;
  }
  /* Berlekamp-Massey; the recurrence has degree at most D */
  long L = 0, m = 1;
  mp_limb_t b = 1;
  memset(C, 0, (D + 1) * sizeof(mp_limb_t));
  memset(B, 0, (D + 1) * sizeof(mp_limb_t));
  C[0] = B[0] = 1;
  for (long n = 0; n < len; n++) {
    mp_limb_t d = seq[n];
    for (long k = 1; k <= L; k++) {
      d = nmod_add(d, nmod_mul(C[k], seq[n - k], mod), mod);
    }
    if (d == 0) {
      m++;
      continue;
    }
    const mp_limb_t coef = nmod_mul(d, n_invmod(b, mod.n), mod);
    const int update = (2 * L <= n);
    if (update) {
      memcpy(T, C, (D + 1) * sizeof(mp_limb_t));
    }
    for (long k = 0; k + m <= D; k++) {
      if (B[k]) {
        C[k + m] = nmod_sub(C[k + m], nmod_mul(coef, B[k], mod), mod);
      }
    }
    if (update) {
      L = n + 1 - L;
      memcpy(B, T, (D + 1) * sizeof(mp_limb_t));
      b = d;
      m = 1;
    } else {
      m++;
    }
  }
  return L;
}


/* normal forms of the products of the staircase monomials by the
 * variables, from a reduced DRL Groebner basis modulo a prime */
typedef struct {
  int32_t nv;
  long D;
  nmod_t mod;
  long nbd;          /* number of border monomials */
  long *target;      /* target[i * nv + v]: index of lmb[i] * x_v in the
                        staircase, or -(index in the border) - 1 */
  mp_limb_t **nf;    /* normal forms of the border monomials */
  mp_limb_t *nfdata;
  long idx0;         /* index of 1 in the staircase */
} sep_mulmat_t;

static void sep_mulmat_clear(sep_mulmat_t *mm) {
  free(mm->target);
  free(mm->nf);
  free(mm->nfdata);
  mm->target = NULL;
  mm->nf = NULL;
  mm->nfdata = NULL;
}

/* returns 1 on success, 0 otherwise (mm is then cleared) */
static int sep_mulmat_init(sep_mulmat_t *mm, bs_t *bs, int32_t *lmb, long D,
                           uint32_t fc) {
  const ht_t *ht = bs->ht;
  const int32_t nv = ht->nv;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  mm->nv = nv;
  mm->D = D;
  mm->nbd = 0;
  mm->target = NULL;
  mm->nf = NULL;
  mm->nfdata = NULL;
  mm->idx0 = -1;
  if (bs->cf_32 == NULL || ebl != 0 || D < 2) {
    return 0;
  }
  nmod_init(&mm->mod, fc);
  const nmod_t mod = mm->mod;

  /* Groebner basis elements, monic */
  const long ngb = bs->lml;
  long nterms = 0;
  for (long i = 0; i < ngb; i++) {
    nterms += bs->hm[bs->lmps[i]][LENGTH];
  }
  int32_t *gexps = malloc(nterms * nv * sizeof(int32_t));
  mp_limb_t *gcfs = malloc(nterms * sizeof(mp_limb_t));
  long *gstart = malloc((ngb + 1) * sizeof(long));
  int32_t *lms = malloc(ngb * nv * sizeof(int32_t));
  long pos = 0;
  for (long i = 0; i < ngb; i++) {
    const bl_t bi = bs->lmps[i];
    const long len = bs->hm[bi][LENGTH];
    const hm_t *dt = bs->hm[bi] + OFFSET;
    const cf32_t *cf = bs->cf_32[bs->hm[bi][COEFFS]];
    const mp_limb_t inv = n_invmod(cf[0], fc);
    gstart[i] = pos;
    for (long j = 0; j < len; j++) {
      int32_t *e = gexps + (pos + j) * nv;
      long ce = 0;
      for (len_t k = 1; k < ebl; ++k) {
        e[ce++] = (int32_t)ht->ev[dt[j]][k];
      }
      for (len_t k = ebl + 1; k < evl; ++k) {
        e[ce++] = (int32_t)ht->ev[dt[j]][k];
      }
      gcfs[pos + j] = nmod_mul(cf[j], inv, mod);
    }
    memcpy(lms + i * nv, gexps + pos * nv, nv * sizeof(int32_t));
    pos += len;
  }
  gstart[ngb] = pos;

  sep_hash_t Hb, Hlm, Hbd;
  sep_hash_init(&Hb, lmb, D, nv);
  sep_hash_init(&Hlm, lms, ngb, nv);

  int32_t *e = calloc(nv, sizeof(int32_t));
  mm->idx0 = sep_hash_find(&Hb, e);

  /* border monomials */
  long nbd = 0, abd = D * nv;
  int32_t *bd = malloc(abd * nv * sizeof(int32_t));
  for (long i = 0; i < D; i++) {
    for (int32_t v = 0; v < nv; v++) {
      memcpy(e, lmb + i * nv, nv * sizeof(int32_t));
      e[v]++;
      if (sep_hash_find(&Hb, e) < 0) {
        memcpy(bd + nbd * nv, e, nv * sizeof(int32_t));
        nbd++;
      }
    }
  }
  sep_cmp_nv = nv;
  qsort(bd, nbd, nv * sizeof(int32_t), sep_cmp_drl);
  long k = 0;
  for (long i = 0; i < nbd; i++) {
    if (k == 0 || memcmp(bd + (k - 1) * nv, bd + i * nv, nv * sizeof(int32_t))) {
      memmove(bd + k * nv, bd + i * nv, nv * sizeof(int32_t));
      k++;
    }
  }
  nbd = k;
  mm->nbd = nbd;
  sep_hash_init(&Hbd, bd, nbd, nv);

  long *target = malloc(D * nv * sizeof(long));
  for (long i = 0; i < D; i++) {
    for (int32_t v = 0; v < nv; v++) {
      memcpy(e, lmb + i * nv, nv * sizeof(int32_t));
      e[v]++;
      long t = sep_hash_find(&Hb, e);
      target[i * nv + v] = t >= 0 ? t : -sep_hash_find(&Hbd, e) - 1;
    }
  }
  mm->target = target;

  /* normal forms of the border monomials, by increasing DRL order */
  mp_limb_t **nf = malloc(nbd * sizeof(mp_limb_t *));
  mp_limb_t *nfdata = calloc(nbd * D, sizeof(mp_limb_t));
  mm->nf = nf;
  mm->nfdata = nfdata;
  int ok = mm->idx0 >= 0;
  for (long s = 0; s < nbd && ok; s++) {
    nf[s] = nfdata + s * D;
    const int32_t *t = bd + s * nv;
    long g = sep_hash_find(&Hlm, t);
    if (g >= 0) {
      for (long j = gstart[g] + 1; j < gstart[g + 1]; j++) {
        long idx = sep_hash_find(&Hb, gexps + j * nv);
        if (idx < 0) {
          ok = 0;
          break;
        }
        nf[s][idx] = nmod_neg(gcfs[j], mod);
      }
      continue;
    }
    /* t = x_j t' with t' in the border, t' < t */
    long sp = -1;
    int32_t j;
    for (j = 0; j < nv; j++) {
      if (t[j] == 0) {
        continue;
      }
      memcpy(e, t, nv * sizeof(int32_t));
      e[j]--;
      if (sep_hash_find(&Hb, e) < 0 && (sp = sep_hash_find(&Hbd, e)) >= 0) {
        break;
      }
    }
    if (sp < 0 || sp >= s) {
      ok = 0;
      break;
    }
    for (long i = 0; i < D; i++) {
      const mp_limb_t w = nf[sp][i];
      if (w == 0) {
        continue;
      }
      const long tg = target[i * nv + j];
      if (tg >= 0) {
        nf[s][tg] = nmod_add(nf[s][tg], w, mod);
      } else if (-tg - 1 >= s) {
        ok = 0;
        break;
      } else {
        _nmod_vec_scalar_addmul_nmod(nf[s], nf[-tg - 1], D, w, mod);
      }
    }
  }
  free(e);
  free(bd);
  free(Hb.tab);
  free(Hlm.tab);
  free(Hbd.tab);
  free(lms);
  free(gstart);
  free(gcfs);
  free(gexps);
  if (!ok) {
    sep_mulmat_clear(mm);
  }
  return ok;
}

/* builds the dense matrix of the multiplication by sum c[v] x_v,
 * with 0 <= c[v] < mod.n */
static void sep_form_matrix(mp_limb_t *mat, const int32_t *c,
                            const sep_mulmat_t *mm) {
  const long D = mm->D;
  const int32_t nv = mm->nv;
  const nmod_t mod = mm->mod;
  memset(mat, 0, D * D * sizeof(mp_limb_t));
  for (long i = 0; i < D; i++) {
    for (int32_t v = 0; v < nv; v++) {
      if (c[v] == 0) {
        continue;
      }
      const mp_limb_t cv = (mp_limb_t)c[v];
      const long t = mm->target[i * nv + v];
      if (t >= 0) {
        mat[t * D + i] = nmod_add(mat[t * D + i], cv, mod);
      } else {
        const mp_limb_t *col = mm->nf[-t - 1];
        for (long k = 0; k < D; k++) {
          mat[k * D + i] = nmod_add(mat[k * D + i], nmod_mul(cv, col[k], mod), mod);
        }
      }
    }
  }
}

/* Parametrization of the radical of the ideal w.r.t. the separating
 * linear form t = -sum c[v] x_v, without FGLM (t is the extra variable A
 * of the linear form A + sum c[v] x_v added to the system). With M the
 * dense matrix of the multiplication by t, u a random row vector and
 * w_j = u M^j, the minimal polynomial P of M is the one of
 * s_j = w_j[1] (j < 2 dquot), and the eliminating polynomial is its
 * squarefree part Q. With s_j^(v) = w_j . NF(x_v) (j < dquot) and N_v the
 * polynomial part of P(T) sum_j s_j^(v) T^(-j-1), one has
 * x_v = N_v / N_1 mod Q: at a root of Q of multiplicity m in P, the
 * coefficients of 1 / (T - t(z))^m in the generating series give
 * N_v / N_1 = x_v(z), as soon as (t - t(z))^(m-1) kills the maximal ideal
 * of the local algebra at z, i.e. P has the degree of the minimal
 * polynomial of a generic linear form. param has nv + 1 variables, the
 * last one being t; as with FGLM, its coordinates are the opposites of
 * the x_v and its denominator is 1.
 * Returns 0 on success, 1 if the prime is bad, that is, if P and Q do not
 * have degrees degelim and degsqfr. */
static int sep_dense_param(param_t *param, const sep_mulmat_t *mm,
                           const int32_t *c, long degelim, long degsqfr) {
  const long D = mm->D;
  const int32_t nv = mm->nv;
  const nmod_t mod = mm->mod;
  const long len = 2 * D;
  const long *target = mm->target;

  /* mt[i * D + k]: coefficient of the k-th staircase monomial in
   * NF(t lmb[i]), i.e. M transposed, so that w M is computed with
   * contiguous dot products */
  mp_limb_t *mt = calloc((size_t)D * D, sizeof(mp_limb_t));
  for (long i = 0; i < D; i++) {
    for (int32_t v = 0; v < nv; v++) {
      if (c[v] == 0) {
        continue;
      }
      const mp_limb_t cv = nmod_neg((mp_limb_t)c[v] % mod.n, mod);
      const long t = target[i * nv + v];
      if (t >= 0) {
        mt[i * D + t] = nmod_add(mt[i * D + t], cv, mod);
      } else {
        _nmod_vec_scalar_addmul_nmod(mt + i * D, mm->nf[-t - 1], D, cv, mod);
      }
    }
  }
  /* own generator: do not perturb rand() */
  uint64_t state = 0x2545F4914F6CDD1DULL ^ mod.n;
  mp_limb_t *w = malloc(D * sizeof(mp_limb_t));
  mp_limb_t *w2 = malloc(D * sizeof(mp_limb_t));
  for (long i = 0; i < D; i++) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    w[i] = (state >> 33) % mod.n;
  }
  mp_limb_t *seq = malloc(len * sizeof(mp_limb_t));
  mp_limb_t *seqx = malloc((size_t)nv * D * sizeof(mp_limb_t));
  const long idx0 = mm->idx0;
  dot_params_t params = _nmod_vec_dot_params(D, mod);
  for (long j = 0; j < len; j++) {
    seq[j] = w[idx0];
    if (j < D) {
      for (int32_t v = 0; v < nv; v++) {
        const long t = target[idx0 * nv + v];
        seqx[v * D + j] = t >= 0 ? w[t]
          : _nmod_vec_dot(w, mm->nf[-t - 1], D, mod, params);
      }
    }
    if (j == len - 1) {
      break;
    }
    for (long i = 0; i < D; i++) {
      w2[i] = _nmod_vec_dot(mt + i * D, w, D, mod, params);
    }
    mp_limb_t *tmp = w;
    w = w2;
    w2 = tmp;
  }
  free(mt);
  free(w);
  free(w2);

  int bad = 1;
  nmod_berlekamp_massey_t BM;
  nmod_berlekamp_massey_init(BM, mod.n);
  nmod_berlekamp_massey_add_points(BM, seq, len);
  nmod_berlekamp_massey_reduce(BM);
  nmod_poly_t P, Q, N, inv;
  nmod_poly_init(P, mod.n);
  nmod_poly_init(Q, mod.n);
  nmod_poly_init(N, mod.n);
  nmod_poly_init(inv, mod.n);
  nmod_poly_make_monic(P, nmod_berlekamp_massey_V_poly(BM));
  const long d = nmod_poly_degree(P);
  /* squarefree part, d < mod.n */
  nmod_poly_derivative(N, P);
  nmod_poly_gcd(inv, P, N);
  nmod_poly_div(Q, P, inv);
  mp_limb_t *rev = malloc(D * sizeof(mp_limb_t));
  if (d == degelim && nmod_poly_degree(Q) == degsqfr) {
    fglm_param_set_prime(param, mod.n);
    bad = 0;
    /* N(T) = sum_{m < d} (sum_{j < d - m} p_{m + 1 + j} s_j) T^m */
    for (int32_t v = -1; v < nv && !bad; v++) {
      const mp_limb_t *s = v < 0 ? seq : seqx + v * D;
      nmod_poly_fit_length(N, d);
      for (long m = 0; m < d; m++) {
        for (long j = 0; j < d - m; j++) {
          rev[j] = P->coeffs[m + 1 + j];
        }
        N->coeffs[m] = _nmod_vec_dot(rev, s, d - m, mod, params);
      }
      _nmod_poly_set_length(N, d);
      _nmod_poly_normalise(N);
      nmod_poly_rem(N, N, Q);
      if (v < 0) {
        bad = !nmod_poly_invmod(inv, N, Q);
      } else {
        nmod_poly_mulmod(param->coords[v], N, inv, Q);
        nmod_poly_neg(param->coords[v], param->coords[v]);
      }
    }
    if (!bad) {
      nmod_poly_set(param->elim, Q);
      nmod_poly_one(param->denom);
      param->degelimpol = d;
      param->degsqfrelimpol = degsqfr;
    }
  }
  free(rev);
  nmod_poly_clear(P);
  nmod_poly_clear(Q);
  nmod_poly_clear(N);
  nmod_poly_clear(inv);
  nmod_berlekamp_massey_clear(BM);
  free(seq);
  free(seqx);
  return bad;
}

/* Degrees of the minimal polynomials of linear forms, for the search
 * below. In general, the dense matrix of the form is built, and the degree
 * comes from a Wiedemann sequence: O(dquot^3) per form. For a radical
 * ideal, a random form separates the solutions, and its parametrization
 * x_v = X_v(t) mod Q (see sep_dense_param) is computed once: the degree
 * for sum c[v] x_v is then the one of the minimal polynomial of
 * a = sum c[v] X_v in F_p[t]/(Q), from the sequence l(a^k), k < 2 dquot,
 * with O(dquot) products modulo Q. Both are the number of values of the
 * form at the solutions, up to an unlucky choice of the projections. */
typedef struct {
  const sep_mulmat_t *mm;
  mp_limb_t *mat, *u, *r, *x, *y, *seq, *C, *B, *T;
  param_t *par;          /* NULL: dense matrices */
  nmod_poly_t Qinv;      /* inverse of the reverse of Q, for mulmod_preinv */
  mp_limb_t *l;          /* random projection */
} sep_degree_t;

static long sep_degree(sep_degree_t *sd, const int32_t *c) {
  const sep_mulmat_t *mm = sd->mm;
  const long D = mm->D;
  const nmod_t mod = mm->mod;
  if (sd->par == NULL) {
    sep_form_matrix(sd->mat, c, mm);
    return sep_minpoly_degree(sd->mat, D, sd->u, sd->r, mod, sd->x, sd->y,
                              sd->seq, sd->C, sd->B, sd->T);
  }
  const param_t *par = sd->par;
  const long d = nmod_poly_degree(par->elim);
  nmod_poly_t a, pw;
  nmod_poly_init(a, mod.n);
  nmod_poly_init(pw, mod.n);
  for (int32_t v = 0; v < mm->nv; v++) {
    if (c[v] != 0) {
      /* X_v = -coords[v] */
      nmod_poly_scalar_mul_nmod(pw, par->coords[v],
                                nmod_neg((mp_limb_t)c[v] % mod.n, mod));
      nmod_poly_add(a, a, pw);
    }
  }
  nmod_poly_one(pw);
  for (long k = 0; k < 2 * d; k++) {
    sd->seq[k] = _nmod_vec_dot(sd->l, pw->coeffs, pw->length, mod,
                               _nmod_vec_dot_params(pw->length, mod));
    if (k < 2 * d - 1) {
      nmod_poly_mulmod_preinv(pw, pw, a, par->elim, sd->Qinv);
    }
  }
  nmod_berlekamp_massey_t BM;
  nmod_berlekamp_massey_init(BM, mod.n);
  nmod_berlekamp_massey_add_points(BM, sd->seq, 2 * d);
  nmod_berlekamp_massey_reduce(BM);
  const long deg = nmod_poly_degree(nmod_berlekamp_massey_V_poly(BM));
  nmod_berlekamp_massey_clear(BM);
  nmod_poly_clear(a);
  nmod_poly_clear(pw);
  return deg;
}

/* switches sd to the parametrization of the random form c, which
 * separates the solutions of a radical ideal; keeps the dense matrices if
 * it fails */
static void sep_degree_set_param(sep_degree_t *sd, const int32_t *c) {
  const sep_mulmat_t *mm = sd->mm;
  const long D = mm->D;
  param_t *par = allocate_fglm_param(mm->mod.n, mm->nv + 1);
  if (sep_dense_param(par, mm, c, D, D)) {
    free_fglm_param(par);
    return;
  }
  sd->par = par;
  nmod_poly_t rev;
  nmod_poly_init(rev, mm->mod.n);
  nmod_poly_init(sd->Qinv, mm->mod.n);
  nmod_poly_reverse(rev, par->elim, par->elim->length);
  nmod_poly_inv_series(sd->Qinv, rev, par->elim->length);
  nmod_poly_clear(rev);
  sd->l = malloc(D * sizeof(mp_limb_t));
  /* own generator: do not perturb rand() */
  uint64_t state = 0x6A09E667F3BCC909ULL ^ mm->mod.n;
  for (long i = 0; i < D; i++) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    sd->l[i] = (state >> 33) % mm->mod.n;
  }
}

/* returns 0 if nothing was found, 1 if variable *sep_var is separating,
 * 2 if the linear form sep_form (nv coefficients) is separating */
static int find_separating_element(bs_t *bs, int32_t *lmb, long D,
                                   uint32_t fc, int32_t *sep_var,
                                   int32_t *sep_form, int32_t *nsols,
                                   int32_t info_level) {
  double rt = realtime();
  int ret = 0;

  sep_mulmat_t mm;
  int ok = sep_mulmat_init(&mm, bs, lmb, D, fc);
  if (mm.nbd == 0) {
    return 0;
  }
  const int32_t nv = mm.nv;
  if (info_level > 0) {
    fprintf(VERBSTREAM, "\nSearching for a separating element: %ld border monomials, "
            "multiplication matrices in %.2f sec\n", mm.nbd, realtime() - rt);
  }

  if (ok) {
    const nmod_t mod = mm.mod;
    mp_limb_t *mat = malloc(D * D * sizeof(mp_limb_t));
    mp_limb_t *u = malloc(D * sizeof(mp_limb_t));
    mp_limb_t *r = malloc(D * sizeof(mp_limb_t));
    mp_limb_t *x = malloc(D * sizeof(mp_limb_t));
    mp_limb_t *y = malloc(D * sizeof(mp_limb_t));
    mp_limb_t *seq = malloc(2 * D * sizeof(mp_limb_t));
    mp_limb_t *C = malloc((D + 1) * sizeof(mp_limb_t));
    mp_limb_t *B = malloc((D + 1) * sizeof(mp_limb_t));
    mp_limb_t *T = malloc((D + 1) * sizeof(mp_limb_t));
    /* own generator: do not perturb rand() */
    uint64_t state = 0x9E3779B97F4A7C15ULL ^ fc;
    for (long i = 0; i < D; i++) {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      u[i] = (state >> 33) % fc;
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      r[i] = (state >> 33) % fc;
    }
    int32_t *c = calloc(nv, sizeof(int32_t));
    long deg;

    /* the minimal polynomial of a dense random form gives the degree to
     * reach: dquot for a radical ideal; a form reaching it separates the
     * solutions and has the generic nilpotency index at each of them (see
     * sep_dense_param) */
    for (int32_t v = 0; v < nv; v++) {
      state = state * 6364136223846793005ULL + 1442695040888963407ULL;
      c[v] = 1 + (int32_t)((state >> 33) % (fc - 1));
    }
    sep_form_matrix(mat, c, &mm);
    deg = sep_minpoly_degree(mat, D, u, r, mod, x, y, seq, C, B, T);
    const long dgen = deg;
    {
      /* number of distinct solutions */
      nmod_berlekamp_massey_t BM;
      nmod_poly_t P, dP, G;
      nmod_berlekamp_massey_init(BM, fc);
      nmod_berlekamp_massey_add_points(BM, seq, 2 * D);
      nmod_berlekamp_massey_reduce(BM);
      nmod_poly_init(dP, fc);
      nmod_poly_init(G, fc);
      nmod_poly_init(P, fc);
      nmod_poly_set(P, nmod_berlekamp_massey_V_poly(BM));
      nmod_poly_derivative(dP, P);
      nmod_poly_gcd(G, P, dP);
      *nsols = (int32_t)(nmod_poly_degree(P) - nmod_poly_degree(G));
      nmod_poly_clear(P);
      nmod_poly_clear(dP);
      nmod_poly_clear(G);
      nmod_berlekamp_massey_clear(BM);
    }
    if (info_level > 0) {
      fprintf(VERBSTREAM, "random linear form: minimal polynomial of degree %ld, "
              "squarefree part of degree %d (dimension %ld)\n", dgen, *nsols, D);
    }
    sep_degree_t sd = {&mm, mat, u, r, x, y, seq, C, B, T, NULL, {{0}}, NULL};
    if (dgen == D && *nsols == D) {
      sep_degree_set_param(&sd, c);
    }
    long *degv = calloc(nv, sizeof(long));
    if (dgen > 0) {
      for (int32_t v = 0; v < nv && ret == 0; v++) {
        memset(c, 0, nv * sizeof(int32_t));
        c[v] = 1;
        deg = sep_degree(&sd, c);
        degv[v] = deg;
        if (info_level > 1) {
          fprintf(VERBSTREAM, "  variable %d: degree %ld\n", v + 1, deg);
        }
        if (deg == dgen) {
          *sep_var = v;
          ret = 1;
        }
      }
      long ntested = 0;
      if (ret == 0 && nv > 1) {
        /* greedy construction, variables by decreasing degree and, for
         * equal degrees, from the smallest to the largest one */
        int32_t *ord = malloc(nv * sizeof(int32_t));
        for (int32_t v = 0; v < nv; v++) {
          ord[v] = nv - 1 - v;
        }
        for (int32_t i = 1; i < nv; i++) {
          const int32_t w = ord[i];
          int32_t j = i;
          while (j > 0 && degv[ord[j - 1]] < degv[w]) {
            ord[j] = ord[j - 1];
            j--;
          }
          ord[j] = w;
        }
        memset(c, 0, nv * sizeof(int32_t));
        c[ord[0]] = 1;
        long dcur = degv[ord[0]];
        for (int32_t i = 1; i < nv && dcur < dgen; i++) {
          const int32_t v = ord[i];
          /* number of values of (form, x_v): degree of form + rnd x_v */
          state = state * 6364136223846793005ULL + 1442695040888963407ULL;
          c[v] = 1 + (int32_t)((state >> 33) % (fc - 1));
          const long tgt = sep_degree(&sd, c);
          ntested++;
          /* smallest coefficient reaching tgt; there are at most
           * dquot (dquot - 1) / 2 coefficients that do not, but stop at
           * dquot and keep the best one if it increases the degree,
           * else 0 */
          int32_t best = 0;
          long dbest = dcur;
          for (int32_t cv = 1; cv <= D && tgt > dcur; cv++) {
            c[v] = cv;
            deg = sep_degree(&sd, c);
            ntested++;
            if (deg > dbest) {
              best = cv;
              dbest = deg;
            }
            if (deg >= tgt) {
              break;
            }
          }
          c[v] = best;
          dcur = dbest;
          if (info_level > 1) {
            fprintf(VERBSTREAM, "  coefficient %d on variable %d: degree %ld (target %ld)\n",
                    best, v + 1, dcur, tgt);
          }
        }
        free(ord);
        if (dcur == dgen) {
          memcpy(sep_form, c, nv * sizeof(int32_t));
          ret = 2;
        }
      }
      if (info_level > 0) {
        if (ret == 1) {
          fprintf(VERBSTREAM, "variable %d is separating\n", *sep_var + 1);
        } else if (ret == 2) {
          fprintf(VERBSTREAM, "separating linear form after %ld forms:", ntested);
          for (int32_t v = 0; v < nv; v++) {
            fprintf(VERBSTREAM, " %d", sep_form[v]);
          }
          fprintf(VERBSTREAM, "\n");
        } else {
          fprintf(VERBSTREAM, "no separating element found (%ld forms tested)\n",
                  ntested);
        }
      }
    }
    if (sd.par != NULL) {
      free_fglm_param(sd.par);
      nmod_poly_clear(sd.Qinv);
      free(sd.l);
    }
    free(degv);
    free(c);
    free(mat);
    free(u);
    free(r);
    free(x);
    free(y);
    free(seq);
    free(C);
    free(B);
    free(T);
    sep_mulmat_clear(&mm);
  }
  if (info_level > 0) {
    fprintf(VERBSTREAM, "separating element search: %.2f sec\n", realtime() - rt);
  }
  return ret;
}
