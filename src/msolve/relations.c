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

/* Zero and equal coordinates of the real solutions.

   At a root t of the eliminating polynomial w, the coordinates are
   x_v = -C_v(t) / (c_v D(t)), with C_v = param->coords[v],
   c_v = param->cfs[v] and D = param->denom. The relations between the
   values 0, t and x_0, .., x_{nc - 1} at a root, t = 0, x_v = 0, x_v = t
   and x_i = x_j, are the vanishing at this root of
     t,   C_v,   C_v + c_v t D,   c_j C_i - c_i C_j.
   The roots of w are split into groups with the same relations, so that
   the coordinates which vanish at a real root are given exactly as 0, those
   which are equal are given the same interval, and the other ones are
   refined until their intervals exclude 0.

   Modulo a prime p which does not divide lc(w) and such that w is
   squarefree modulo p, the roots with some relation are those of
   G = gcd(w, prod R mod w), the product being over the relations R which
   do not hold at all the roots (those are checked exactly over QQ). When
   G = 1, which is the usual case, nothing more is needed: the cost is that
   of a few products modulo w and of one gcd modulo p. Otherwise, G is split
   by the relations as in splitting.c, and its groups, of small degree in
   practice, are reconstructed over QQ from their images modulo further
   primes, matched by their sets of relations. They are then certified
   exactly: their product divides w, and on each group, the relations of a
   spanning forest of the classes of equal values (the relations modulo p
   are closed by transitivity) are divisible by the polynomial of the
   group. Since deg gcd(w, R) modulo p is at least its degree over QQ, a
   relation then holds at no other root of w. */

#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/nmod_poly.h>

/* primes of the reconstruction of the groups, at most */
#define REL_MAX_PRIMES 512
/* primes at which the groups differ from those of the first prime, at most */
#define REL_MAX_MISMATCH 8

typedef struct {
  long nc;            /* number of coordinates */
  long ne;            /* number of relations */
  long *ea, *eb;      /* relation e between the elements ea[e] < eb[e]: 0
                         (element 0), t (element 1) and x_v (element 2 + v) */
  long nwords;        /* words of a set of relations */
  uint64_t *ident;    /* relations at all the roots of w */
  long ng;            /* groups of roots with some other relation */
  fmpz_poly_struct *g;  /* their primitive polynomials over ZZ */
  fmpz_poly_t cof;    /* w divided by their product, when they are known */
  uint64_t **sig;     /* their relations */
  long **rep;         /* rep[i][v], for the group i < ng, and rep[ng] for the
                         other roots: -1 if x_v = 0, -2 if x_v = t, otherwise
                         the first coordinate equal to x_v */
  int known;          /* 0 if the relations could not be determined */
  long nprimes;       /* primes used */
  long nident;        /* relations at all the roots of w */
  long nroots;        /* roots of the groups */
  double time;        /* elapsed time of real_rel_init */
} real_rel_t;

/* the parametrization modulo p */
typedef struct {
  nmod_poly_t w, td;  /* w and t D */
  nmod_poly_struct *C;
  mp_limb_t *c;
} rel_mod_t;

static inline int rel_test(const uint64_t *s, const long e) {
  return (s[e / 64] >> (e % 64)) & 1;
}

static inline void rel_set(uint64_t *s, const long e) {
  s[e / 64] |= UINT64_C(1) << (e % 64);
}

static void rel_mpz_to_nmod(nmod_poly_t r, const mpz_upoly_t u) {
  const mp_limb_t p = r->mod.n;
  nmod_poly_fit_length(r, u->length);
  for (long i = 0; i < u->length; i++) {
    r->coeffs[i] = mpz_fdiv_ui(u->coeffs[i], p);
  }
  _nmod_poly_set_length(r, u->length);
  _nmod_poly_normalise(r);
}

static void rel_fmpz_poly(fmpz_poly_t r, const mpz_upoly_t u) {
  fmpz_poly_fit_length(r, u->length);
  for (long i = 0; i < u->length; i++) {
    fmpz_set_mpz(r->coeffs + i, u->coeffs[i]);
  }
  _fmpz_poly_set_length(r, u->length);
  _fmpz_poly_normalise(r);
}

static void rel_mod_init(rel_mod_t *m, const mpz_param_t param,
                         const long nc, const mp_limb_t p) {
  nmod_poly_init(m->w, p);
  nmod_poly_init(m->td, p);
  rel_mpz_to_nmod(m->w, param->elim);
  rel_mpz_to_nmod(m->td, param->denom);
  nmod_poly_shift_left(m->td, m->td, 1);
  m->C = malloc(nc * sizeof(nmod_poly_struct));
  m->c = malloc(nc * sizeof(mp_limb_t));
  for (long v = 0; v < nc; v++) {
    nmod_poly_init(m->C + v, p);
    rel_mpz_to_nmod(m->C + v, param->coords[v]);
    m->c[v] = mpz_fdiv_ui(param->cfs[v], p);
  }
}

static void rel_mod_clear(rel_mod_t *m, const long nc) {
  nmod_poly_clear(m->w);
  nmod_poly_clear(m->td);
  for (long v = 0; v < nc; v++) {
    nmod_poly_clear(m->C + v);
  }
  free(m->C);
  free(m->c);
}

/* the relation e modulo p, reduced modulo w */
static void rel_poly_mod(nmod_poly_t r, const real_rel_t *rel, const long e,
                         const rel_mod_t *m) {
  const long a = rel->ea[e], b = rel->eb[e];
  if (a == 0 && b == 1) {
    nmod_poly_zero(r);
    nmod_poly_set_coeff_ui(r, 1, 1);
  } else if (a == 0) {
    nmod_poly_set(r, m->C + (b - 2));
  } else if (a == 1) {
    const long v = b - 2;
    nmod_poly_scalar_mul_nmod(r, m->td, m->c[v]);
    nmod_poly_add(r, r, m->C + v);
  } else {
    const long i = a - 2, j = b - 2;
    nmod_poly_t s;
    nmod_poly_init_mod(s, m->w->mod);
    nmod_poly_scalar_mul_nmod(r, m->C + i, m->c[j]);
    nmod_poly_scalar_mul_nmod(s, m->C + j, m->c[i]);
    nmod_poly_sub(r, r, s);
    nmod_poly_clear(s);
  }
  nmod_poly_rem(r, r, m->w);
}

/* the relation e over ZZ */
static void rel_poly_zz(fmpz_poly_t r, const real_rel_t *rel, const long e,
                        const mpz_param_t param) {
  const long a = rel->ea[e], b = rel->eb[e];
  fmpz_poly_t s;
  fmpz_poly_init(s);
  fmpz_t c;
  fmpz_init(c);
  if (a == 0 && b == 1) {
    fmpz_poly_zero(r);
    fmpz_poly_set_coeff_ui(r, 1, 1);
  } else if (a <= 1) {
    const long v = b - 2;
    rel_fmpz_poly(r, param->coords[v]);
    if (a == 1) {
      rel_fmpz_poly(s, param->denom);
      fmpz_poly_shift_left(s, s, 1);
      fmpz_set_mpz(c, param->cfs[v]);
      fmpz_poly_scalar_addmul_fmpz(r, s, c);
    }
  } else {
    const long i = a - 2, j = b - 2;
    rel_fmpz_poly(r, param->coords[i]);
    fmpz_set_mpz(c, param->cfs[j]);
    fmpz_poly_scalar_mul_fmpz(r, r, c);
    rel_fmpz_poly(s, param->coords[j]);
    fmpz_set_mpz(c, param->cfs[i]);
    fmpz_poly_scalar_submul_fmpz(r, s, c);
  }
  fmpz_clear(c);
  fmpz_poly_clear(s);
}

/* whether the relation e holds at all the roots of w, exactly */
static int rel_ident_zz(const real_rel_t *rel, const long e,
                        const mpz_param_t param) {
  fmpz_poly_t r, w;
  fmpz_poly_init(r);
  fmpz_poly_init(w);
  rel_poly_zz(r, rel, e, param);
  int ok = fmpz_poly_is_zero(r);
  rel_fmpz_poly(w, param->elim);
  if (!ok && fmpz_poly_degree(r) == fmpz_poly_degree(w)) {
    /* r, of degree at most deg(w), is a multiple of w over QQ */
    fmpz_poly_t s;
    fmpz_poly_init(s);
    fmpz_poly_scalar_mul_fmpz(s, r, fmpz_poly_lead(w));
    fmpz_poly_scalar_submul_fmpz(s, w, fmpz_poly_lead(r));
    ok = fmpz_poly_is_zero(s);
    fmpz_poly_clear(s);
  }
  fmpz_poly_clear(r);
  fmpz_poly_clear(w);
  return ok;
}

typedef struct {
  uint64_t *sig;
  nmod_poly_t g;      /* monic */
} rel_group_t;

static void rel_groups_clear(rel_group_t *G, const long n) {
  for (long i = 0; i < n; i++) {
    free(G[i].sig);
    nmod_poly_clear(G[i].g);
  }
  free(G);
}

static long rel_cmp_nwords;

static int rel_cmp_groups(const void *x, const void *y) {
  const rel_group_t *a = x, *b = y;
  for (long w = 0; w < rel_cmp_nwords; w++) {
    if (a->sig[w] != b->sig[w]) {
      return a->sig[w] < b->sig[w] ? -1 : 1;
    }
  }
  return 0;
}

/* groups of the roots of w with some relation besides those of rel->ident,
 * modulo p, sorted by their relations (rel->ident included) in *out;
 * returns their number, or -1 if p is not usable: p divides lc(w), w is not
 * squarefree modulo p (checked when first is set), or a relation which does
 * not hold at all the roots of w holds at all of them modulo p */
static long rel_groups_mod(rel_group_t **out, const real_rel_t *rel,
                           const mpz_param_t param, const mp_limb_t p,
                           const int first) {
  const long nc = rel->nc;
  rel_mod_t m;
  rel_mod_init(&m, param, nc, p);
  long n = -1;
  rel_group_t *G = NULL;
  nmod_poly_struct *R = NULL;
  long nr = 0;
  nmod_poly_t P, h;
  nmod_poly_init(P, p);
  nmod_poly_init(h, p);
  if (nmod_poly_degree(m.w) != param->elim->length - 1) {
    goto clear;
  }
  if (first) {
    nmod_poly_derivative(h, m.w);
    nmod_poly_gcd(h, m.w, h);
    if (nmod_poly_degree(h) > 0) {
      goto clear;
    }
  }
  /* the relations, and the roots with some of them */
  R = malloc(rel->ne * sizeof(nmod_poly_struct));
  nmod_poly_one(P);
  for (long e = 0; e < rel->ne; e++) {
    if (rel_test(rel->ident, e)) {
      continue;
    }
    nmod_poly_init(R + nr, p);
    rel_poly_mod(R + nr, rel, e, &m);
    if (nmod_poly_is_zero(R + nr)) {
      nr++;
      goto clear;
    }
    nmod_poly_mulmod(P, P, R + nr, m.w);
    nr++;
  }
  nmod_poly_gcd(h, m.w, P);
  n = 0;
  if (nmod_poly_degree(h) <= 0) {
    goto clear;
  }
  /* split along the relations, as split_mod in splitting.c */
  long alloc = 16;
  G = malloc(alloc * sizeof(rel_group_t));
  G[0].sig = malloc(rel->nwords * sizeof(uint64_t));
  memcpy(G[0].sig, rel->ident, rel->nwords * sizeof(uint64_t));
  nmod_poly_init(G[0].g, p);
  nmod_poly_make_monic(G[0].g, h);
  n = 1;
  for (long e = 0, k = 0; e < rel->ne; e++) {
    if (rel_test(rel->ident, e)) {
      continue;
    }
    const long n0 = n;
    for (long i = 0; i < n0; i++) {
      nmod_poly_rem(P, R + k, G[i].g);
      if (nmod_poly_is_zero(P)) {
        rel_set(G[i].sig, e);
        continue;
      }
      nmod_poly_gcd(h, G[i].g, P);
      if (nmod_poly_degree(h) <= 0) {
        continue;
      }
      if (n == alloc) {
        alloc *= 2;
        G = realloc(G, alloc * sizeof(rel_group_t));
      }
      G[n].sig = malloc(rel->nwords * sizeof(uint64_t));
      memcpy(G[n].sig, G[i].sig, rel->nwords * sizeof(uint64_t));
      rel_set(G[n].sig, e);
      nmod_poly_init(G[n].g, p);
      nmod_poly_make_monic(G[n].g, h);
      nmod_poly_div(G[i].g, G[i].g, h);
      n++;
    }
    k++;
  }
  rel_cmp_nwords = rel->nwords;
  qsort(G, n, sizeof(rel_group_t), rel_cmp_groups);
clear:
  for (long k = 0; k < nr; k++) {
    nmod_poly_clear(R + k);
  }
  free(R);
  nmod_poly_clear(P);
  nmod_poly_clear(h);
  rel_mod_clear(&m, nc);
  *out = G;
  return n;
}

static long rel_find(long *parent, long x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}

/* classes of equal values for the relations of s: parent[x] is the
 * smallest element of the class of x */
static void rel_classes(long *parent, const real_rel_t *rel,
                        const uint64_t *s) {
  const long ne = rel->nc + 2;
  for (long x = 0; x < ne; x++) {
    parent[x] = x;
  }
  for (long e = 0; e < rel->ne; e++) {
    if (rel_test(s, e)) {
      long a = rel_find(parent, rel->ea[e]);
      long b = rel_find(parent, rel->eb[e]);
      if (a != b) {
        parent[FLINT_MAX(a, b)] = FLINT_MIN(a, b);
      }
    }
  }
  for (long x = 0; x < ne; x++) {
    parent[x] = rel_find(parent, x);
  }
}

/* index of the relation between the elements a < b */
static long rel_index(const real_rel_t *rel, const long a, const long b) {
  for (long e = 0; e < rel->ne; e++) {
    if (rel->ea[e] == a && rel->eb[e] == b) {
      return e;
    }
  }
  return -1;
}

/* checks exactly the groups rel->g: their product divides w, and each one
 * divides the relations between the elements of its classes and the
 * smallest element of these classes */
static int rel_certify(real_rel_t *rel, const mpz_param_t param) {
  const long nc = rel->nc;
  fmpz_poly_t w, P, q, r;
  fmpz_poly_init(w);
  fmpz_poly_init(P);
  fmpz_poly_init(q);
  fmpz_poly_init(r);
  rel_fmpz_poly(w, param->elim);
  fmpz_poly_one(P);
  for (long i = 0; i < rel->ng; i++) {
    fmpz_poly_mul(P, P, rel->g + i);
  }
  int ok = fmpz_poly_divides(q, w, P);
  if (ok) {
    fmpz_poly_set(rel->cof, q);
  }
  /* the product of the groups on which the relation e is checked */
  fmpz_poly_struct *Pe = malloc(rel->ne * sizeof(fmpz_poly_struct));
  for (long e = 0; e < rel->ne; e++) {
    fmpz_poly_init(Pe + e);
    fmpz_poly_one(Pe + e);
  }
  long *parent = malloc((nc + 2) * sizeof(long));
  for (long i = 0; i < rel->ng && ok; i++) {
    rel_classes(parent, rel, rel->sig[i]);
    for (long x = 0; x < nc + 2; x++) {
      if (parent[x] != x) {
        const long e = rel_index(rel, parent[x], x);
        if (rel_test(rel->ident, e)) {
          continue;
        }
        fmpz_poly_mul(Pe + e, Pe + e, rel->g + i);
      }
    }
  }
  for (long e = 0; e < rel->ne && ok; e++) {
    if (fmpz_poly_degree(Pe + e) > 0) {
      rel_poly_zz(r, rel, e, param);
      ok = fmpz_poly_divides(q, r, Pe + e);
    }
  }
  for (long e = 0; e < rel->ne; e++) {
    fmpz_poly_clear(Pe + e);
  }
  free(Pe);
  free(parent);
  fmpz_poly_clear(w);
  fmpz_poly_clear(P);
  fmpz_poly_clear(q);
  fmpz_poly_clear(r);
  return ok;
}

/* rel->rep from the classes of the groups */
static void rel_set_reps(real_rel_t *rel) {
  const long nc = rel->nc;
  long *parent = malloc((nc + 2) * sizeof(long));
  rel->rep = malloc((rel->ng + 1) * sizeof(long *));
  for (long i = 0; i <= rel->ng; i++) {
    rel_classes(parent, rel, (i < rel->ng) ? rel->sig[i] : rel->ident);
    rel->rep[i] = malloc(nc * sizeof(long));
    for (long v = 0; v < nc; v++) {
      const long x = parent[2 + v];
      rel->rep[i][v] = (x == 0) ? -1 : ((x == 1) ? -2 : x - 2);
    }
  }
  free(parent);
}

/* reconstruction over QQ of the monic polynomial whose images modulo m are
 * the coefficients of u, as a primitive polynomial over ZZ */
static int rel_reconstruct(fmpz_poly_t g, const fmpz_poly_t u,
                           const fmpz_t m) {
  fmpq_poly_t gq;
  fmpq_poly_init(gq);
  fmpq_t c;
  fmpq_init(c);
  fmpz_t a;
  fmpz_init(a);
  int ok = 1;
  for (long i = 0; i < fmpz_poly_length(u) && ok; i++) {
    fmpz_poly_get_coeff_fmpz(a, u, i);
    /* fmpq_reconstruct_fmpz needs a residue in [0, m) */
    fmpz_mod(a, a, m);
    ok = fmpq_reconstruct_fmpz(c, a, m);
    fmpq_poly_set_coeff_fmpq(gq, i, c);
  }
  if (ok) {
    fmpq_poly_get_numerator(g, gq);
    fmpz_poly_primitive_part(g, g);
  }
  fmpz_clear(a);
  fmpq_clear(c);
  fmpq_poly_clear(gq);
  return ok;
}

static void real_rel_clear(real_rel_t *rel) {
  fmpz_poly_clear(rel->cof);
  free(rel->ea);
  free(rel->eb);
  free(rel->ident);
  for (long i = 0; i < rel->ng; i++) {
    fmpz_poly_clear(rel->g + i);
    free(rel->sig[i]);
  }
  free(rel->g);
  free(rel->sig);
  if (rel->rep != NULL) {
    for (long i = 0; i <= rel->ng; i++) {
      free(rel->rep[i]);
    }
    free(rel->rep);
  }
}

/* the relations between the coordinates at the roots of the eliminating
 * polynomial of param (see above); rel->known is 0 if they could not be
 * determined */
static void real_rel_init(real_rel_t *rel, const mpz_param_t param) {
  double rt = realtime();
  const long nc = param->nvars - 1;
  rel->nc = nc;
  rel->ne = 0;
  rel->ng = 0;
  rel->g = NULL;
  rel->sig = NULL;
  rel->rep = NULL;
  rel->known = 0;
  rel->nprimes = 0;
  fmpz_poly_init(rel->cof);
  /* relations between 0, t and the x_v */
  const long nel = nc + 2;
  rel->ea = malloc(nel * nel * sizeof(long));
  rel->eb = malloc(nel * nel * sizeof(long));
  for (long a = 0; a < nel; a++) {
    for (long b = a + 1; b < nel; b++) {
      rel->ea[rel->ne] = a;
      rel->eb[rel->ne] = b;
      rel->ne++;
    }
  }
  rel->nwords = (rel->ne + 63) / 64 + 1;
  rel->ident = calloc(rel->nwords, sizeof(uint64_t));

  /* the first prime, and the relations at all the roots of w */
  mp_limb_t p = n_nextprime(UWORD(1) << 62, 1);
  rel_group_t *G1 = NULL;
  long n1 = -1;
  for (long tries = 0; tries < REL_MAX_MISMATCH && n1 < 0; tries++) {
    rel_mod_t m;
    rel_mod_init(&m, param, nc, p);
    nmod_poly_t r;
    nmod_poly_init(r, p);
    int bad = 0;
    for (long e = 0; e < rel->ne && !bad; e++) {
      if (rel_test(rel->ident, e)) {
        continue;
      }
      rel_poly_mod(r, rel, e, &m);
      if (nmod_poly_is_zero(r)) {
        if (rel_ident_zz(rel, e, param)) {
          rel_set(rel->ident, e);
        } else {
          bad = 1;
        }
      }
    }
    nmod_poly_clear(r);
    rel_mod_clear(&m, nc);
    if (!bad) {
      n1 = rel_groups_mod(&G1, rel, param, p, 1);
    }
    rel->nprimes++;
    p = n_nextprime(p, 1);
  }
  if (n1 < 0) {
    goto done;
  }
  if (n1 == 0) {
    rel->known = 1;
    goto done;
  }

  /* the groups over QQ */
  rel->ng = n1;
  rel->g = malloc(n1 * sizeof(fmpz_poly_struct));
  rel->sig = malloc(n1 * sizeof(uint64_t *));
  fmpz_poly_struct *crt = malloc(n1 * sizeof(fmpz_poly_struct));
  fmpz_poly_struct *prev = malloc(n1 * sizeof(fmpz_poly_struct));
  for (long i = 0; i < n1; i++) {
    fmpz_poly_init(rel->g + i);
    fmpz_poly_init(crt + i);
    fmpz_poly_init(prev + i);
    rel->sig[i] = G1[i].sig;
    G1[i].sig = NULL;
    fmpz_poly_set_nmod_poly_unsigned(crt + i, G1[i].g);
  }
  fmpz_t mod;
  fmpz_init_set_ui(mod, G1[0].g->mod.n);
  rel_groups_clear(G1, n1);
  long mismatch = 0, stable = 0;
  while (rel->nprimes < REL_MAX_PRIMES && mismatch < REL_MAX_MISMATCH) {
    int ok = 1;
    for (long i = 0; i < n1 && ok; i++) {
      ok = rel_reconstruct(rel->g + i, crt + i, mod);
    }
    if (ok) {
      for (long i = 0; i < n1 && ok; i++) {
        ok = fmpz_poly_equal(rel->g + i, prev + i);
      }
      stable = ok ? stable + 1 : 0;
      for (long i = 0; i < n1; i++) {
        fmpz_poly_set(prev + i, rel->g + i);
      }
      if (stable && rel_certify(rel, param)) {
        rel->known = 1;
        break;
      }
    }
    rel_group_t *Gk;
    const long nk = rel_groups_mod(&Gk, rel, param, p, 0);
    rel->nprimes++;
    int match = (nk == n1);
    for (long i = 0; i < n1 && match; i++) {
      match = (memcmp(Gk[i].sig, rel->sig[i], rel->nwords * sizeof(uint64_t)) == 0
               && nmod_poly_degree(Gk[i].g) == fmpz_poly_degree(crt + i));
    }
    if (match) {
      for (long i = 0; i < n1; i++) {
        fmpz_poly_CRT_ui(crt + i, crt + i, mod, Gk[i].g, 0);
      }
      fmpz_mul_ui(mod, mod, p);
    } else {
      mismatch++;
    }
    if (nk > 0) {
      rel_groups_clear(Gk, nk);
    }
    p = n_nextprime(p, 1);
  }
  for (long i = 0; i < n1; i++) {
    fmpz_poly_clear(crt + i);
    fmpz_poly_clear(prev + i);
  }
  free(crt);
  free(prev);
  fmpz_clear(mod);
done:
  if (rel->known) {
    rel_set_reps(rel);
  }
  rel->nident = 0;
  for (long e = 0; e < rel->ne; e++) {
    rel->nident += rel_test(rel->ident, e);
  }
  rel->nroots = 0;
  for (long i = 0; i < rel->ng; i++) {
    rel->nroots += fmpz_poly_degree(rel->g + i);
  }
  rel->time = realtime() - rt;
}

/* prints the relations of the group i (rel->ng for the other roots), as
 * classes of equal values: the coordinates, then t, then 0; vnames are the
 * names of the variables, t being the last one */
static void real_rel_print(FILE *f, const real_rel_t *rel, const long i,
                           char **vnames) {
  const long ne = rel->nc + 2;
  long *parent = malloc(ne * sizeof(long));
  rel_classes(parent, rel, (i < rel->ng) ? rel->sig[i] : rel->ident);
  int any = 0;
  for (long r = 0; r < ne; r++) {
    long n = 0;
    for (long x = r; x < ne; x++) {
      n += (parent[x] == r);
    }
    if (parent[r] != r || n < 2) {
      continue;
    }
    fprintf(f, "%s", any ? ", " : "");
    any = 1;
    long m = 0;
    for (long x = 2; x < ne; x++) {
      if (parent[x] == r) {
        fprintf(f, "%s%s", m++ ? " = " : "", vnames[x - 2]);
      }
    }
    if (parent[1] == r) {
      fprintf(f, "%s%s", m++ ? " = " : "", vnames[rel->nc]);
    }
    if (r == 0) {
      fprintf(f, " = 0");
    }
  }
  if (!any) {
    fprintf(f, "no zero or equal coordinates");
  }
  free(parent);
}
