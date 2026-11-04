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

#include <flint/nmod_mat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Threshold on D for parallelism in fglmc_mulvec and fglmc_lu_reduce:
 * below it, the sequential path is used. Not tuned yet. */
#define FGLMC_PAR_MIN 256

/* Threads for FLINT (products, echelon form, block solving) :
 * activated modulo only one prime */
static int fglmc_flint_threads_push(int t){
  const int old = flint_get_num_threads();
  if(t > 1) {
      flint_set_num_threads(t);
  }
  return old;
}
static void fglmc_flint_threads_pop(int old){
  flint_set_num_threads(old);
}

#define FGLM_ORDER_LEX     0
#define FGLM_ORDER_GREVLEX 1

static inline CF_t fglmc_addmod(CF_t a, CF_t b, CF_t p){
  CF_t s = a + b; return s >= p ? s - p : s;          /* p < 2^31 */
}
static inline CF_t fglmc_negmod(CF_t a, CF_t p){ return a ? p - a : 0; }
static inline CF_t fglmc_mulmod(CF_t a, CF_t b, CF_t p){
  return (CF_t)(((uint64_t)a * b) % p);
}

static int fglmc_cmp(const int32_t *a, const int32_t *b, int32_t nv, int ord){
  if(ord == FGLM_ORDER_LEX){
    for(int32_t i = 0; i < nv; i++)
      if(a[i] != b[i]) return a[i] > b[i] ? 1 : -1;
    return 0;
  }
  int64_t da = 0, db = 0;
  for(int32_t i = 0; i < nv; i++){ da += a[i]; db += b[i]; }
  if(da != db) return da > db ? 1 : -1;
  for(int32_t i = nv - 1; i >= 0; i--)
    if(a[i] != b[i]) return a[i] < b[i] ? 1 : -1;
  return 0;
}

/* delayed accumulation on 64 bits */
static void fglmc_mulvec(CF_t *w, const sp_multimatfglm_t *M, int32_t i,
                         const CF_t *v, uint64_t *acc){
  const szmat_t D = M->ncols;
  const CF_t p = M->charac;
  const uint64_t pp = (uint64_t)p * p;
  const uint64_t Mred = ((UINT64_C(1) << 63) / pp) * pp;
  const szmat_t *sel = M->rowsel + (uint64_t)i * D;
#ifdef _OPENMP
  if(M->nthreads > 1 && D >= FGLMC_PAR_MIN){
    /* parallelisation on the columns */
#pragma omp parallel num_threads(M->nthreads)
    {
      const int tid = omp_get_thread_num(), nth = omp_get_num_threads();
      const szmat_t l0 = (szmat_t)((uint64_t)D * tid / nth);
      const szmat_t l1 = (szmat_t)((uint64_t)D * (tid + 1) / nth);
      memset(acc + l0, 0, (l1 - l0) * sizeof(uint64_t));
      for(szmat_t k = 0; k < D; k++){
        const uint64_t a = v[k];
        if(!a) continue;
        const szmat_t s = sel[k];
        if(s < D){
          if(s >= l0 && s < l1){
            uint64_t t = acc[s] + a;
            acc[s] = t - (Mred & -(t >> 63));
          }
        } else {
          const CF_t *row = M->dense_mat + (uint64_t)(s - D) * D;
          for(szmat_t l = l0; l < l1; l++){
            uint64_t t = acc[l] + a * row[l];
            acc[l] = t - (Mred & -(t >> 63));
          }
        }
      }
      for(szmat_t l = l0; l < l1; l++) w[l] = (CF_t)(acc[l] % p);
    }
    return;
  }
#endif
  memset(acc, 0, D * sizeof(uint64_t));
  for(szmat_t k = 0; k < D; k++){
    const uint64_t a = v[k];
    if(!a) continue;
    const szmat_t s = sel[k];
    if(s < D){
      uint64_t t = acc[s] + a;
      acc[s] = t - (Mred & -(t >> 63));
    } else {
      const CF_t *row = M->dense_mat + (uint64_t)(s - D) * D;
      for(szmat_t l = 0; l < D; l++){
        uint64_t t = acc[l] + a * row[l];
        acc[l] = t - (Mred & -(t >> 63));
      }
    }
  }
  for(szmat_t l = 0; l < D; l++) w[l] = (CF_t)(acc[l] % p);
}

/* Incremental LU decomposition */
typedef struct {
  szmat_t D, r;
  mp_limb_t p;
  nmod_mat_t U;      /* D x D, rows 0..r-1                           */
  nmod_mat_t Pt;     /* D x D, triangular block to reduce            */
  CF_t *L;           /* D x D, lower triangular                      */
  mp_limb_t *dinv;   /* inverse of L[t][t]                           */
  slong *piv;        /* pivots                                       */
  int nthr;          /* number of threads                            */
  uint64_t *tmp;
} fglmc_lu_t;

static int fglmc_lu_init(fglmc_lu_t *E, szmat_t D, CF_t p){
  E->D = D;
  E->r = 0;
  E->p = p;
  E->nthr = 1;
  nmod_mat_init(E->U, D, D, p);
  nmod_mat_init(E->Pt, D, D, p);
  E->L = malloc((uint64_t)D * D * sizeof(CF_t) + sizeof(CF_t));
  E->dinv = malloc((D + 1) * sizeof(mp_limb_t));
  E->piv = malloc((D + 1) * sizeof(slong));
  E->tmp = malloc((D + 1) * sizeof(uint64_t));
  return E->L && E->dinv && E->piv && E->tmp;
}
static void fglmc_lu_clear(fglmc_lu_t *E){
  nmod_mat_clear(E->U);
  nmod_mat_clear(E->Pt);
  free(E->L);
  free(E->dinv);
  free(E->piv);
  free(E->tmp);
}

/* Reduce w by the rows U_j, j0 <= j < r.
 * Returns 1 if w becomes 0 */
static int fglmc_lu_reduce(fglmc_lu_t *E, mp_limb_t *w, slong j0, mp_limb_t *a){
  const slong D = E->D, r = E->r;
  const mp_limb_t p = E->p;
  const uint64_t pp = (uint64_t)p * p;
  const uint64_t Mred = ((UINT64_C(1) << 63) / pp) * pp;
#ifdef _OPENMP
  const int par = E->nthr > 1 && D >= FGLMC_PAR_MIN && r - j0 > 1;
#else
  const int par = 0;
#endif
  if(par){
    for(slong j = j0; j < r; j++){
      uint64_t t = w[E->piv[j]];
      const mp_limb_t *Pj = nmod_mat_entry_ptr(E->Pt, j, 0);
      for(slong i = j0; i < j; i++){
        const mp_limb_t ai = a[i - j0];
        if(ai && Pj[i]){
            t += (uint64_t)(p - ai) * Pj[i];
            t -= Mred & -(t >> 63);
        }
      }
      a[j - j0] = (mp_limb_t)(t % p);
    }
#pragma omp parallel for num_threads(E->nthr) schedule(static)
    for(slong l = 0; l < D; l++){
      uint64_t t = w[l];
      for(slong j = j0; j < r; j++){
        const mp_limb_t aj = a[j - j0];
        const mp_limb_t x = nmod_mat_entry(E->U, j, l);
        if(aj && x){
            t += (uint64_t)(p - aj) * x;
            t -= Mred & -(t >> 63);
        }
      }
      w[l] = (mp_limb_t)(t % p);
    }
  } else {
    for(slong j = j0; j < r; j++){
      const slong pj = E->piv[j];
      const mp_limb_t aj = w[pj];
      a[j - j0] = aj;
      if(!aj) continue;
      const mp_limb_t na = p - aj;
      const mp_limb_t *Uj = nmod_mat_entry_ptr(E->U, j, 0);
      for(slong l = pj; l < D; l++)
        if(Uj[l]) w[l] = (mp_limb_t)((w[l] + (uint64_t)na * Uj[l]) % p);
    }
  }
  for(slong l = 0; l < D; l++){
      if(w[l]){
          return 0;
      }
  }
  return 1;
}

/* Insert w if reduced and non zero */
static void fglmc_lu_insert(fglmc_lu_t *E, const mp_limb_t *w, const mp_limb_t *a){
  const slong D = E->D, r = E->r;
  const mp_limb_t p = E->p;
  slong pv = 0;
  while(!w[pv]) {
      pv++;
  }
  const mp_limb_t d = w[pv];
  const mp_limb_t inv = n_invmod(d, p);
  mp_limb_t *Ur = nmod_mat_entry_ptr(E->U, r, 0);
  for(slong l = 0; l < D; l++){
      Ur[l] = w[l] ? (mp_limb_t)(((uint64_t)w[l] * inv) % p) : 0;
  }
  CF_t *Lr = E->L + (uint64_t)r * D;
  for(slong j = 0; j < r; j++){
      Lr[j] = (CF_t)a[j];
  }
  Lr[r] = (CF_t)d;
  E->dinv[r] = inv;
  mp_limb_t *Pr = nmod_mat_entry_ptr(E->Pt, r, 0);
  for(slong j = 0; j < r; j++){
      Pr[j] = nmod_mat_entry(E->U, j, pv);
  }
  Pr[r] = 1;
  E->piv[r] = pv;
  E->r++;
}

static void fglmc_lu_relation(fglmc_lu_t *E, const mp_limb_t *a, CF_t *c){
  const slong r = E->r, D = E->D;
  const mp_limb_t p = E->p;
  uint64_t *res = E->tmp;
  for(slong k = 0; k < r; k++){
      res[k] = a[k];
  }
  for(slong m = r - 1; m >= 0; m--){
    const mp_limb_t xm = (mp_limb_t)((res[m] % p) * (uint64_t)E->dinv[m] % p);
    c[m] = xm ? (CF_t)(p - xm) : 0;
    if(!xm) continue;
    const mp_limb_t nx = p - xm;
    const CF_t *Lm = E->L + (uint64_t)m * D;
    for(slong k = 0; k < m; k++){
      if(Lm[k]) res[k] = (res[k] + (uint64_t)nx * Lm[k]) % p;
    }
  }
}

static void fglmc_lu_reduce_block(fglmc_lu_t *E, nmod_mat_t B, nmod_mat_t C){
  const slong K = B->r, r = E->r, D = E->D;
  const mp_limb_t p = E->p;
  nmod_mat_init(C, K, r, p);
  if(r == 0 || K == 0) return;
  nmod_mat_t Bt, Ct, P, Ptw, Uw;
  nmod_mat_init(Bt, r, K, p);
  nmod_mat_init(Ct, r, K, p);
  for(slong j = 0; j < r; j++){
    for(slong i = 0; i < K; i++){
      nmod_mat_entry(Bt, j, i) = nmod_mat_entry(B, i, E->piv[j]);
    }
  }
  nmod_mat_window_init(Ptw, E->Pt, 0, 0, r, r);
  nmod_mat_solve_tril(Ct, Ptw, Bt, 1);
  nmod_mat_window_clear(Ptw);
  nmod_mat_transpose(C, Ct);
  nmod_mat_init(P, K, D, p);
  nmod_mat_window_init(Uw, E->U, 0, 0, r, D);
  nmod_mat_mul(P, C, Uw);
  nmod_mat_window_clear(Uw);
  nmod_mat_sub(B, B, P);
  nmod_mat_clear(P); nmod_mat_clear(Bt); nmod_mat_clear(Ct);
}

/* Row rank profile of B */
static slong fglmc_first_dependent_row(const nmod_mat_t B){
  const slong K = B->r, D = B->c;
  nmod_mat_t T;
  nmod_mat_init(T, D, K, B->mod.n);
  nmod_mat_transpose(T, B);
  slong rank = nmod_mat_rref(T);
  slong j0 = 0, col = 0;
  for(slong i = 0; i < rank; i++){
    while(col < K && nmod_mat_entry(T, i, col) == 0) {
      col++;
    }
    if(col != j0) break;   /* row j0 is not a pivot */
    j0++;
    col++;
  }
  nmod_mat_clear(T);
  return j0;
}

/* Candidates */
typedef struct {
    const int32_t *e;
    szmat_t parent;
    int32_t var;
} fglmc_cand_t;
typedef struct {
    fglmc_cand_t *a;
    uint64_t n;
    int32_t nv;
    int ord;
} fglmc_heap_t;

static inline int fglmc_hless(const fglmc_heap_t *H, uint64_t i, uint64_t j){
  return fglmc_cmp(H->a[i].e, H->a[j].e, H->nv, H->ord) < 0;
}
static void fglmc_hpush(fglmc_heap_t *H, fglmc_cand_t c){
  uint64_t i = H->n++;
  H->a[i] = c;
  while(i && fglmc_hless(H, i, (i - 1) / 2)){
    fglmc_cand_t t = H->a[i];
    H->a[i] = H->a[(i - 1) / 2];
    H->a[(i - 1) / 2] = t;
    i = (i - 1) / 2;
  }
}
static fglmc_cand_t fglmc_hpop(fglmc_heap_t *H){
  fglmc_cand_t top = H->a[0];
  H->a[0] = H->a[--H->n];
  uint64_t i = 0;
  for(;;){
    uint64_t l = 2 * i + 1, r = l + 1, m = i;
    if(l < H->n && fglmc_hless(H, l, m)){
        m = l;
    }
    if(r < H->n && fglmc_hless(H, r, m)){
        m = r;
    }
    if(m == i){
        break;
    }
    fglmc_cand_t t = H->a[i];
    H->a[i] = H->a[m];
    H->a[m] = t;
    i = m;
  }
  return top;
}

/* FGLM for lex
 * Monomials m without any xn are scanned in increasing lex order
 * For each of them, we deal with m x_n^k by block
 * Fills lex: the staircase S, the leading monomials g_lm and g_ntail
 * Returns 1 on success, 0 otherwise (in particular if ord is not lex). */
#define FGLMC_KMAX 256   // max size of a block
#define FGLMC_VLEN 16    // vectors handled one by one

/* m x_n^k is added to the staircase: S[t] = m x_n^k, then t++ */
static inline void fglmc_new_stair(const int32_t *m, const int32_t k,
                                   const int32_t nv, int32_t *S,
                                   szmat_t *t){
    int32_t *e = S + (uint64_t)(*t) * nv;
    memcpy(e, m, nv * sizeof(int32_t));
    e[nv - 1] = k;
    (*t)++;
}

/* m x_n^k is a new leading monomial of the lex Gb; its tail lives in
 * s_0, ..., s_{t-1} (the t current elements of the staircase) */
static inline void fglmc_new_head(fglm_lex_staircase_t *lex,
                                  const int32_t *m, const int32_t k,
                                  const int32_t nv, const szmat_t t){
    int32_t *g = lex->g_lm + (uint64_t)lex->ngens * nv;
    memcpy(g, m, nv * sizeof(int32_t));
    g[nv - 1] = k;
    lex->g_ntail[lex->ngens] = t;
    lex->ncoeffs += t;
    lex->ngens++;
}

/* Relation for the current leading monomial: its coefficients against
 * s_0, ..., s_{t-1} are written in coeffs[lex->ncoeffs .. lex->ncoeffs + t),
 * coeffs being enlarged if needed. Returns 0 if the allocation fails. */
static int fglmc_push_relation(fglmc_lu_t *E, const mp_limb_t *al,
                               const fglm_lex_staircase_t *lex,
                               const szmat_t t,
                               CF_t **coeffs, uint64_t *cap_c){
    if(lex->ncoeffs + t > *cap_c){
        uint64_t c = *cap_c;
        while(lex->ncoeffs + t > c) c *= 2;
        CF_t *nc = realloc(*coeffs, c * sizeof(CF_t));
        if(nc == NULL) return 0;
        *coeffs = nc;
        *cap_c = c;
    }
    fglmc_lu_relation(E, al, *coeffs + lex->ncoeffs);
    return 1;
}

static int fglm_classic_learn(const sp_multimatfglm_t *Mt, const int ord,
                              fglm_lex_staircase_t *lex, CF_t **coeffs_out){
  if(ord != FGLM_ORDER_LEX) return 0;
  const szmat_t D = Mt->ncols;
  const int32_t nv = Mt->nvars;
  const int32_t xn = nv - 1;
  const CF_t p = Mt->charac;
  const uint64_t maxc = (uint64_t)nv * D + 1;
  int ok = 1;

  fglmc_lu_t E;
  ok = fglmc_lu_init(&E, D, p);
  E.nthr = Mt->nthreads > 1 ? Mt->nthreads : 1;
  /* D + FGLMC_KMAX rows: rows t, t+1, ... are buffers */
  CF_t *Vs    = malloc((uint64_t)(D + FGLMC_KMAX + 1) * D * sizeof(CF_t));
  int32_t *S  = malloc((uint64_t)D * nv * sizeof(int32_t));
  int32_t *pool = malloc(maxc * nv * sizeof(int32_t));
  uint64_t *acc = malloc(D * sizeof(uint64_t));
  mp_limb_t *wl = malloc((D + 1) * sizeof(mp_limb_t));
  mp_limb_t *al = malloc((D + 1) * sizeof(mp_limb_t));
  uint64_t cap_c = 4 * (uint64_t)D + 16;
  CF_t *coeffs = malloc(cap_c * sizeof(CF_t));
  fglmc_heap_t H = { malloc(maxc * sizeof(fglmc_cand_t)), 0, nv, ord };
  lex->g_lm     = malloc(maxc * nv * sizeof(int32_t));
  lex->g_ntail  = malloc(maxc * sizeof(szmat_t));
  lex->ngens = 0; lex->ncoeffs = 0;
  if(!ok || !Vs || !S || !pool || !acc || !wl || !al || !coeffs || !H.a
     || !lex->g_lm || !lex->g_ntail){
      ok = 0;
      goto done;
  }

  uint64_t npool = 1;
  szmat_t t = 0;
  memset(pool, 0, nv * sizeof(int32_t));
  fglmc_hpush(&H,
              (fglmc_cand_t){
                  pool,
                  UINT32_MAX, -1
              }
              );
  const int32_t *last = NULL;

  while(H.n && ok){
    fglmc_cand_t u = fglmc_hpop(&H);
    if(last && !memcmp(last, u.e, nv * sizeof(int32_t))) continue;
    last = u.e;
    /* is it a multiple of a lm? if not bound kmax */
    int skip = 0;
    uint64_t kmax = D + 1;
    for(int32_t l = 0; l < lex->ngens && !skip; l++){
      const int32_t *lm = lex->g_lm + (uint64_t)l * nv;
      int div = 1;
      for(int32_t i = 0; i < xn && div; i++) if(lm[i] > u.e[i]) div = 0;
      if(!div) continue;
      if(lm[xn] == 0){
          skip = 1;
      }
      else if((uint64_t)lm[xn] < kmax) {
          kmax = lm[xn];
      }
    }
    if(skip) continue;

    const szmat_t t_head = t;
    uint64_t k = 0;
    slong K = FGLMC_VLEN;
    int chain_done = 0;
    while(!chain_done && ok){
      if(k < FGLMC_VLEN){
        if(t >= D + 1){
            ok = 0;
            break;
        }
        CF_t *row = Vs + (uint64_t)t * D;
        if(k == 0){
          if(u.parent == UINT32_MAX){
              memset(row, 0, D * sizeof(CF_t));
              row[Mt->one] = 1;
          }
          else {
              fglmc_mulvec(row, Mt, u.var,
                           Vs + (uint64_t)u.parent * D,
                           acc);
          }
        } else {
          fglmc_mulvec(row, Mt, xn, Vs + (uint64_t)(t - 1) * D,
                       acc);
        }
        for(szmat_t l = 0; l < D; l++){
            wl[l] = row[l];
        }
        if(fglmc_lu_reduce(&E, wl, 0, al)){
            if(!fglmc_push_relation(&E, al, lex, t, &coeffs, &cap_c)){
                ok = 0;
                break;
            }
            fglmc_new_head(lex, u.e, (int32_t)k, nv, t);
            chain_done = 1;
        } else {
          if(t == D){
              ok = 0;
              break;
          }
          fglmc_lu_insert(&E, wl, al);
          fglmc_new_stair(u.e, (int32_t)k, nv, S, &t);
          k++;
          if(k >= kmax) {
              chain_done = 1;          /* multiple of a lm */
          }
        }
        continue;
      }
      if((uint64_t)K > kmax - k) {
          K = (slong)(kmax - k);
      }
      for(slong i = 0; i < K; i++){
          fglmc_mulvec(Vs + (uint64_t)(t + i) * D, Mt, xn,
                       Vs + (uint64_t)(t + i - 1) * D, acc);
      }
      nmod_mat_t B, C;
      nmod_mat_init(B, K, D, p);
      for(slong i = 0; i < K; i++){
        const CF_t *row = Vs + (uint64_t)(t + i) * D;
        for(szmat_t l = 0; l < D; l++){
            nmod_mat_entry(B, i, l) = row[l];
        }
      }
      const slong r_old = E.r;
      fglmc_lu_reduce_block(&E, B, C);          /* B -= C U, C : K x r_old */
      slong j0 = fglmc_first_dependent_row(B);
      if(t + j0 > D){
          nmod_mat_clear(B);
          nmod_mat_clear(C);
          ok = 0;
          break;
      }
      /* independent rows */
      const slong nrows_dep = (j0 < K) ? j0 + 1 : j0;
      for(slong i = 0; i < nrows_dep && ok; i++){
          for(szmat_t l = 0; l < D; l++){
              wl[l] = nmod_mat_entry(B, i, l);
          }
          for(slong j = 0; j < r_old; j++){
              al[j] = nmod_mat_entry(C, i, j);
          }
          const int zero = fglmc_lu_reduce(&E, wl, r_old, al + r_old);
          if(i < j0){
              if(zero){
                  ok = 0;
                  break;
              }
              fglmc_lu_insert(&E, wl, al);
              fglmc_new_stair(u.e, (int32_t)(k + i), nv, S, &t);
          } else {
              if(!zero){
                  ok = 0;
                  break;
              }
              if (!fglmc_push_relation(&E, al, lex, t, &coeffs, &cap_c)){
                  ok = 0;
                  break;
                  }
              fglmc_new_head(lex, u.e, (int32_t)(k + j0), nv, t);
          }
      }
      nmod_mat_clear(B);
      nmod_mat_clear(C);
      if(!ok) {
          break;
      }
      if(j0 < K){
          chain_done = 1;
      } else {
          k += K;
          if(k >= kmax) {
              chain_done = 1;
          }
          if(2 * K <= FGLMC_KMAX) {
              K *= 2;
          }
      }
    }
    /* monomial m is in the staircase, next ones are x_i m, i < n */
    if(ok && t > t_head){
        for(int32_t i = 0; i < xn; i++){
            int32_t *e = pool + npool++ * nv;
            memcpy(e, u.e, nv * sizeof(int32_t));
            e[i]++;
            fglmc_hpush(&H, (fglmc_cand_t){ e, t_head, i });
        }
    }
  }
  if(t != D) {
      ok = 0;
  }

done:
  if(ok){
      lex->S = S;
      S = NULL;
      lex->D = D;
      lex->nvars = nv;
  }
  if(ok && coeffs_out != NULL){
      *coeffs_out = coeffs;
      coeffs = NULL;
  }
  fglmc_lu_clear(&E);
  free(Vs);
  free(S);
  free(pool);
  free(acc);
  free(wl);
  free(al);
  free(coeffs);
  free(H.a);
  return ok;
}


/* Reduced lex Gb (sparse form) from the staircase and the dense coefficients
 * given by fglm_classic_learn */
static fglm_gb_t *fglm_gb_from_staircase(const fglm_lex_staircase_t *lex,
                                         const CF_t *coeffs, const CF_t p){
  const int32_t nv = lex->nvars;
  fglm_gb_t *G = calloc(1, sizeof(fglm_gb_t));
  uint64_t nterms = lex->ngens;
  for(uint64_t i = 0; i < lex->ncoeffs; i++){
      if(coeffs[i]){
          nterms++;
      }
  }
  G->charac = p;
  G->nvars = nv;
  G->ngens = lex->ngens;
  G->lens = malloc((lex->ngens + 1) * sizeof(int32_t));
  G->exps = malloc((nterms + 1) * nv * sizeof(int32_t));
  G->cfs  = malloc((nterms + 1) * sizeof(CF_t));
  uint64_t off = 0, pos = 0;
  for(int32_t j = 0; j < lex->ngens; j++){
      const szmat_t nt = lex->g_ntail[j];
      uint64_t start = pos;
      memcpy(G->exps + pos * nv, lex->g_lm + (uint64_t)j * nv,
             nv * sizeof(int32_t));
      G->cfs[pos++] = 1;
      for(szmat_t k = nt; k-- > 0; ){
          if(!coeffs[off + k]) continue;
          memcpy(G->exps + pos * nv, lex->S + (uint64_t)k * nv,
                 nv * sizeof(int32_t));
          G->cfs[pos++] = coeffs[off + k];
      }
      G->lens[j] = (int32_t)(pos - start);
      off += nt;
  }
  return G;
}

/* FGLM over a prime field: the greedy phase (LU elimination) gives at once the
 * staircase, the leading monomials and the relations. */
static fglm_gb_t *fglm_classic(const sp_multimatfglm_t *Mt, const int ord){
  fglm_lex_staircase_t lex;
  memset(&lex, 0, sizeof(lex));
  fglm_gb_t *G = NULL;
  CF_t *coeffs = NULL;
  const int old_thr = fglmc_flint_threads_push(Mt->nthreads);
  if(fglm_classic_learn(Mt, ord, &lex, &coeffs)){
      G = fglm_gb_from_staircase(&lex, coeffs, Mt->charac);
  }
  free(coeffs);
  fglmc_flint_threads_pop(old_thr);
  free_fglm_lex_staircase(&lex);
  return G;
}

