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

void print_fglm_classic_data(
        FILE *file,
        const md_t * const st,
        const sp_multimatfglm_t *matrix,
        fglm_lex_staircase_t *lex
                     )
{
  if (st->info_level > 0) {
    fprintf(file, "\n---------------- TIMINGS ----------------\n");
    fprintf(file, "overall(elapsed) %11.2f sec\n", st->fglm_rtime);
    fprintf(file, "overall(cpu) %15.2f sec\n", st->fglm_ctime);
    fprintf(file, "-----------------------------------------\n");
    fprintf(file, "\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(file, "degree of ideal    %16lu\n", (unsigned long)matrix->ncols);
    fprintf(file, "#dense rows        %16lu\n", (unsigned long)matrix->nnfs);
    fprintf(file, "#normal forms      %16lu\n", (unsigned long)(matrix->nnfs - matrix->nlms));
    /* fprintf(file, "total density                   %5.1f%%\n", 100*matrix->totaldensity); */
    /* fprintf(file, "density of the free part        %5.1f%%\n", 100*matrix->freepartdensity); */
    /* if(matrix->nnfs){ */
    /*   fprintf(file, "density of the nonfree part     %5.1f%%\n", 100*matrix->nonfreepartdensity); */
    /* } */
    /* fprintf(file, "deg. elim. pol.    %16lu\n", (unsigned long)param->degelimpol); */
    /* fprintf(file, "deg. sqfr. elim. pol. %13lu\n", (unsigned long)param->degsqfrelimpol); */
    fprintf(file, "size of basis      %16lu\n", (unsigned long)(lex->ngens));
    fprintf(file, "#terms in basis    %16lu\n", (unsigned long)(lex->nterms));
    fprintf(file, "-----------------------------------------\n\n");
  }
}


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

/* At -O2, GCC leaves the inner loops of fglmc_mulvec scalar (very cheap cost
 * model); this cost model lets it vectorize them. */
#if defined(__GNUC__) && !defined(__clang__)
#define FGLMC_VECTORIZE __attribute__((optimize("vect-cost-model=dynamic")))
#else
#define FGLMC_VECTORIZE
#endif

/* delayed accumulation on 64 bits; if a dense row which is not computed
 * (M->have_row) is needed, *missing is set to 1 (if missing != NULL) */
FGLMC_VECTORIZE
static void fglmc_mulvec(CF_t *w, const sp_multimatfglm_t *M, int32_t i,
                         const CF_t *v, uint64_t *acc, int *missing){
  const szmat_t D = M->ncols;
  const CF_t p = M->charac;
  const uint64_t pp = (uint64_t)p * p;
  const uint64_t Mred = ((UINT64_C(1) << 63) / pp) * pp;
  /* dense rows are added m at a time (products < (p-1)^2) to accumulators
   * kept below 2^63, then folded */
  const uint64_t q = (uint64_t)(p - 1) * (p - 1);
  uint32_t m = (uint32_t)(((UINT64_C(1) << 63) - 1) / (q ? q : 1));
  if(m < 1) m = 1;
  if(m > 8) m = 8;
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
      uint32_t nd = 0;
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
          if(M->have_row != NULL && !M->have_row[s - D]){
            if(missing != NULL){
#pragma omp atomic write
              *missing = 1;
            }
            continue;
          }
          const CF_t *row = M->dense_mat + (uint64_t)(s - D) * D;
          for(szmat_t l = l0; l < l1; l++) acc[l] += a * row[l];
          if(++nd == m){
            for(szmat_t l = l0; l < l1; l++) acc[l] -= Mred & -(acc[l] >> 63);
            nd = 0;
          }
        }
      }
      for(szmat_t l = l0; l < l1; l++) w[l] = (CF_t)(acc[l] % p);
    }
    return;
  }
#endif
  memset(acc, 0, D * sizeof(uint64_t));
  uint32_t nd = 0;
  for(szmat_t k = 0; k < D; k++){
    const uint64_t a = v[k];
    if(!a) continue;
    const szmat_t s = sel[k];
    if(s < D){
      uint64_t t = acc[s] + a;
      acc[s] = t - (Mred & -(t >> 63));
    } else {
      if(M->have_row != NULL && !M->have_row[s - D]){
        if(missing != NULL) *missing = 1;
        continue;
      }
      const CF_t *row = M->dense_mat + (uint64_t)(s - D) * D;
      for(szmat_t l = 0; l < D; l++) acc[l] += a * row[l];
      if(++nd == m){
        for(szmat_t l = 0; l < D; l++) acc[l] -= Mred & -(acc[l] >> 63);
        nd = 0;
      }
    }
  }
  for(szmat_t l = 0; l < D; l++) w[l] = (CF_t)(acc[l] % p);
}

/* Incremental LU decomposition */
typedef struct {
  szmat_t D, r;
  mp_limb_t p;
  CF_t *U;           /* D x D, rows 0..r-1                           */
  CF_t *PtT;         /* D x D: PtT[j][i] = U_j[piv[i]] for i >= j     */
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
  E->U = malloc((uint64_t)D * D * sizeof(CF_t) + sizeof(CF_t));
  E->PtT = malloc((uint64_t)D * D * sizeof(CF_t) + sizeof(CF_t));
  E->L = malloc((uint64_t)D * D * sizeof(CF_t) + sizeof(CF_t));
  E->dinv = malloc((D + 1) * sizeof(mp_limb_t));
  E->piv = malloc((D + 1) * sizeof(slong));
  E->tmp = malloc((D + 1) * sizeof(uint64_t));
  return E->U && E->PtT && E->L && E->dinv && E->piv && E->tmp;
}
static void fglmc_lu_clear(fglmc_lu_t *E){
  free(E->U);
  free(E->PtT);
  free(E->L);
  free(E->dinv);
  free(E->piv);
  free(E->tmp);
}

static void fglmc_gemm_sub(mp_limb_t *B, const slong ldb,
                           const mp_limb_t *C, const slong ldc,
                           const CF_t *U, const slong ldu,
                           const slong K, const slong n, const slong w,
                           const mp_limb_t p, const int nthr);

/* Reduce w by the rows U_j, j0 <= j < r.
 * Returns 1 if w becomes 0 */
FGLMC_VECTORIZE
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
      for(slong i = j0; i < j; i++){
        const mp_limb_t ai = a[i - j0];
        const mp_limb_t pji = E->PtT[(uint64_t)i * D + j];
        if(ai && pji){
            t += (uint64_t)(p - ai) * pji;
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
        const mp_limb_t x = E->U[(uint64_t)j * D + l];
        if(aj && x){
            t += (uint64_t)(p - aj) * x;
            t -= Mred & -(t >> 63);
        }
      }
      w[l] = (mp_limb_t)(t % p);
    }
  } else {
    /* coefficients: forward substitution on the pivots, PtT row by row
     * (t_i = w[piv_i] - sum_{j<i} a_j U_j[piv_i]), delayed reduction */
    const slong n = r - j0;
    const uint64_t q = (uint64_t)(p - 1) * (p - 1);
    uint32_t m = (uint32_t)(((UINT64_C(1) << 63) - 1) / (q ? q : 1));
    if(m < 1) m = 1;
    if(m > 8) m = 8;
    uint64_t *t = E->tmp;
    for(slong i = 0; i < n; i++) t[i] = w[E->piv[j0 + i]];
    uint32_t nd = 0;
    for(slong j = 0; j < n; j++){
      const mp_limb_t aj = (mp_limb_t)(t[j] % p);
      a[j] = aj;
      if(!aj) continue;
      const uint64_t na = p - aj;
      const CF_t *Pj = E->PtT + (uint64_t)(j0 + j) * D + j0;
      for(slong i = j + 1; i < n; i++) t[i] += na * Pj[i];
      if(++nd == m){
        for(slong i = j + 1; i < n; i++) t[i] -= Mred & -(t[i] >> 63);
        nd = 0;
      }
    }
    /* w -= a U[j0..r) */
    fglmc_gemm_sub(w, D, a, n, E->U + (uint64_t)j0 * D, D, 1, n, D, p, 1);
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
  CF_t *Ur = E->U + (uint64_t)r * D;
  for(slong l = 0; l < D; l++){
      Ur[l] = w[l] ? (CF_t)(((uint64_t)w[l] * inv) % p) : 0;
  }
  CF_t *Lr = E->L + (uint64_t)r * D;
  for(slong j = 0; j < r; j++){
      Lr[j] = (CF_t)a[j];
  }
  Lr[r] = (CF_t)d;
  E->dinv[r] = inv;
  for(slong j = 0; j < r; j++){
      E->PtT[(uint64_t)j * D + r] = E->U[(uint64_t)j * D + pv];
  }
  E->PtT[(uint64_t)r * D + r] = 1;
  E->piv[r] = pv;
  E->r++;
}

FGLMC_VECTORIZE
static void fglmc_lu_relation(fglmc_lu_t *E, const mp_limb_t *a, CF_t *c){
  const slong r = E->r, D = E->D;
  const mp_limb_t p = E->p;
  const uint64_t pp = (uint64_t)p * p;
  const uint64_t Mred = ((UINT64_C(1) << 63) / pp) * pp;
  const uint64_t q = (uint64_t)(p - 1) * (p - 1);
  uint32_t nf = (uint32_t)(((UINT64_C(1) << 63) - 1) / (q ? q : 1));
  if(nf < 1) nf = 1;
  if(nf > 8) nf = 8;
  /* back substitution in L, delayed reduction */
  uint64_t *res = E->tmp;
  for(slong k = 0; k < r; k++){
      res[k] = a[k];
  }
  uint32_t nd = 0;
  for(slong m = r - 1; m >= 0; m--){
    const mp_limb_t xm = (mp_limb_t)((res[m] % p) * (uint64_t)E->dinv[m] % p);
    c[m] = xm ? (CF_t)(p - xm) : 0;
    if(!xm) continue;
    const uint64_t nx = p - xm;
    const CF_t *Lm = E->L + (uint64_t)m * D;
    for(slong k = 0; k < m; k++) res[k] += nx * Lm[k];
    if(++nd == nf){
      for(slong k = 0; k < m; k++) res[k] -= Mred & -(res[k] >> 63);
      nd = 0;
    }
  }
}

/* B -= C U mod p, with B (K x w), C (K x n) and U (n x w) given by their first
 * entry and their row strides; entries of B and C are reduced modulo p, U
 * has 32-bit entries. Tiles of FGLMC_RB rows and FGLMC_WB columns; the
 * products are accumulated in 64-bit integers kept below 2^63, m at a time
 * before folding (m depends on p: 8 for p < 2^30). */
#ifndef FGLMC_RB
#define FGLMC_RB 32
#endif
#define FGLMC_WB 256
FGLMC_VECTORIZE
static void fglmc_gemm_sub(mp_limb_t *B, const slong ldb,
                           const mp_limb_t *C, const slong ldc,
                           const CF_t *U, const slong ldu,
                           const slong K, const slong n, const slong w,
                           const mp_limb_t p, const int nthr){
  if(K == 0 || n == 0 || w == 0) return;
  const uint64_t pp = (uint64_t)p * p;
  const uint64_t Mred = ((UINT64_C(1) << 63) / pp) * pp;
  const uint64_t q = (uint64_t)(p - 1) * (p - 1);
  slong m = (slong)(((UINT64_C(1) << 63) - 1) / (q ? q : 1));
  if(m < 1) m = 1;
  if(m > 8) m = 8;
  (void)nthr;
  /* column tiles in parallel; in a column tile, the rows j0..j1-1 of U are
   * used by FGLMC_RB rows of C */
#pragma omp parallel for num_threads(nthr) schedule(static) \
  if(nthr > 1 && (uint64_t)K * n * w >= (UINT64_C(1) << 20))
  for(slong l0 = 0; l0 < w; l0 += FGLMC_WB){
    const slong wl = (w - l0 < FGLMC_WB) ? w - l0 : FGLMC_WB;
    uint64_t acc[FGLMC_RB][FGLMC_WB];
    for(slong i0 = 0; i0 < K; i0 += FGLMC_RB){
      const slong nr = (K - i0 < FGLMC_RB) ? K - i0 : FGLMC_RB;
      for(slong a = 0; a < nr; a++) memset(acc[a], 0, wl * sizeof(uint64_t));
      for(slong j0 = 0; j0 < n; j0 += m){
        const slong j1 = (j0 + m < n) ? j0 + m : n;
        for(slong a = 0; a < nr; a++){
          uint64_t *ac = acc[a];
          const mp_limb_t *cr = C + (uint64_t)(i0 + a) * ldc;
          for(slong j = j0; j < j1; j++){
            const uint64_t c = cr[j];
            if(!c) continue;
            const CF_t *u = U + (uint64_t)j * ldu + l0;
            for(slong l = 0; l < wl; l++) ac[l] += c * u[l];
          }
          for(slong l = 0; l < wl; l++) ac[l] -= Mred & -(ac[l] >> 63);
        }
      }
      for(slong a = 0; a < nr; a++){
        mp_limb_t *bl = B + (uint64_t)(i0 + a) * ldb + l0;
        for(slong l = 0; l < wl; l++){
          const mp_limb_t sv = (mp_limb_t)(acc[a][l] % p);
          bl[l] = (bl[l] >= sv) ? bl[l] - sv : bl[l] + p - sv;
        }
      }
    }
  }
}

/* row stride of a FLINT matrix (windows included) */
static inline slong fglmc_stride(const nmod_mat_t A){
  return (A->r > 1) ? (slong)(nmod_mat_entry_ptr(A, 1, 0) - nmod_mat_entry_ptr(A, 0, 0))
                    : A->c;
}

#define FGLMC_NB 64      // panel width of the triangular solve

/* Reduces the rows of B against the rows a..b-1 of E, B being already reduced
 * against the rows 0..a-1: B -= C U[a..b), C (K x (b - a)) receiving the
 * coefficients. */
static void fglmc_lu_reduce_range(fglmc_lu_t *E, nmod_mat_t B, nmod_mat_t C,
                                  const slong a, const slong b){
  const slong K = B->r, n = b - a, D = E->D;
  const mp_limb_t p = E->p;
  nmod_mat_init(C, K, n, p);
  if(n == 0 || K == 0) return;
  const slong ldb = fglmc_stride(B), ldc = fglmc_stride(C);
  mp_limb_t *Cp = nmod_mat_entry_ptr(C, 0, 0);
  /* Bp = B[:, piv[a..b)] */
  mp_limb_t *Bp = malloc((uint64_t)K * n * sizeof(mp_limb_t));
  for(slong i = 0; i < K; i++){
    const mp_limb_t *bi = nmod_mat_entry_ptr(B, i, 0);
    for(slong j = 0; j < n; j++) Bp[(uint64_t)i * n + j] = bi[E->piv[a + j]];
  }
  /* C T = Bp with T = PtT[a..b, a..b] unit upper triangular: by panels of
   * FGLMC_NB columns, substitution inside a panel, then update of the next
   * columns of Bp */
  const CF_t *T = E->PtT + (uint64_t)a * D + a;
  for(slong jb = 0; jb < n; jb += FGLMC_NB){
    const slong je = (jb + FGLMC_NB < n) ? jb + FGLMC_NB : n;
    for(slong i = 0; i < K; i++){
      mp_limb_t *bp = Bp + (uint64_t)i * n;
      mp_limb_t *c = Cp + (uint64_t)i * ldc;
      for(slong l = jb; l < je; l++){
        const mp_limb_t cl = bp[l];
        c[l] = cl;
        if(!cl) continue;
        const uint64_t ncl = p - cl;
        const CF_t *tr = T + (uint64_t)l * D;
        for(slong j = l + 1; j < je; j++){
          if(tr[j]) bp[j] = (mp_limb_t)((bp[j] + ncl * tr[j]) % p);
        }
      }
    }
    if(je < n){
      fglmc_gemm_sub(Bp + je, n, Cp + jb, ldc, T + (uint64_t)jb * D + je, D,
                     K, je - jb, n - je, p, E->nthr);
    }
  }
  free(Bp);
  /* B -= C U[a..b) */
  fglmc_gemm_sub(nmod_mat_entry_ptr(B, 0, 0), ldb, Cp, ldc,
                 E->U + (uint64_t)a * D, D, K, n, D, p, E->nthr);
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
/* block sizes; can be set at compile time (-DFGLMC_KMAX=...) */
#ifndef FGLMC_KMAX
#define FGLMC_KMAX 256   // max size of a block
#endif
#ifndef FGLMC_KSUB
#define FGLMC_KSUB 64    // sub-blocks reduced together inside a block
#endif
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

/* Parent of a new element m x_n^k of the staircase or of a new lm (trace
 * over Q): x_var times the parent of the prefix for k = 0, x_n times the
 * previous element of the chain (s_{t-1}) otherwise */
static inline void fglmc_set_parent(szmat_t *parent, int32_t *var,
                                    const fglmc_cand_t *u, const uint64_t k,
                                    const szmat_t t, const int32_t xn){
    if(k == 0){
        *parent = u->parent;
        *var = u->var;
    } else {
        *parent = t - 1;
        *var = xn;
    }
}

/* If tr != NULL, the parents of the elements of the staircase and of the
 * leading monomials are recorded in tr (trace over Q). */
/* Dense rows computed on demand (normal forms by core_nf, see hilbert.c):
 * fill computes the dense rows rows[0..n) of M, returns 0 on failure. */
typedef struct {
  int (*fill)(void *ctx, sp_multimatfglm_t *M, const szmat_t *rows,
              szmat_t n, const int info_level);
  void *ctx;
} fglmc_rows_t;

/* computes, if needed, the rows of the matrix of x_i used by the product
 * with v. Returns 0 on failure. */
static int fglmc_need_rows(const fglmc_rows_t *rp, sp_multimatfglm_t *M,
                           const int32_t i, const CF_t *v,
                           const int info_level){
  if(rp == NULL || rp->fill == NULL || M->have_row == NULL) return 1;
  const szmat_t D = M->ncols;
  const szmat_t *sel = M->rowsel + (uint64_t)i * D;
  szmat_t *rows = malloc((D + 1) * sizeof(szmat_t));
  szmat_t n = 0;
  for(szmat_t k = 0; k < D; k++){
    if(v[k] && sel[k] >= D && !M->have_row[sel[k] - D]){
      rows[n++] = sel[k] - D;
    }
  }
  int ok = 1;
  if(n > 0) ok = rp->fill(rp->ctx, M, rows, n, info_level);
  free(rows);
  return ok;
}

/* 1 if the product of v by the matrix of x_i needs a dense row which is
 * not computed yet */
static int fglmc_rows_missing(const sp_multimatfglm_t *M, const int32_t i,
                              const CF_t *v){
  if(M->have_row == NULL) return 0;
  const szmat_t D = M->ncols;
  const szmat_t *sel = M->rowsel + (uint64_t)i * D;
  for(szmat_t k = 0; k < D; k++){
    if(v[k] && sel[k] >= D && !M->have_row[sel[k] - D]) return 1;
  }
  return 0;
}

/* Bound on the chain of the prefix e: the chain stops before x_n^kmax.
 * Returns 0 if e is a multiple of a lm free of x_n (no chain). */
static uint64_t fglmc_kmax(const fglm_lex_staircase_t *lex, const int32_t *e,
                           const int32_t nv, const szmat_t D){
  const int32_t xn = nv - 1;
  uint64_t kmax = (uint64_t)D + 1;
  for(int32_t l = 0; l < lex->ngens; l++){
    const int32_t *lm = lex->g_lm + (uint64_t)l * nv;
    int div = 1;
    for(int32_t i = 0; i < xn && div; i++) {
        if(lm[i] > e[i]) {
            div = 0;
        }
    }
    if(!div) continue;
    if(lm[xn] == 0){
        return 0;
    }
    if((uint64_t)lm[xn] < kmax) {
        kmax = lm[xn];
    }
  }
  return kmax;
}

#ifndef FGLMC_PRE
#define FGLMC_PRE 32     // vectors reduced together ahead of the candidates
#endif

/* Candidates taken ahead from the heap (increasing order). The first
 * vectors of their chains (NF(e x_n^k), k < len) are reduced together
 * against the rows 0..r0-1 of the LU, U being read once. */
typedef struct {
  fglmc_cand_t *c;
  uint32_t *off;     // first row of c[i] in V and B
  uint32_t *len;     // number of vectors of c[i]
  uint32_t n, h;     // number of candidates, next one
  CF_t *V;           // FGLMC_PRE x D: vectors, not reduced
  nmod_mat_t B, C;   // reduced vectors, coefficients (K x r0)
  slong r0;
  int init;
} fglmc_pre_t;

static void fglmc_pre_clear(fglmc_pre_t *Q){
  if(Q->init){
      nmod_mat_clear(Q->B);
      nmod_mat_clear(Q->C);
      Q->init = 0;
  }
  Q->n = Q->h = 0;
}

/* Number of vectors taken ahead for a chain whose length is at most b:
 * h[a] counts the previous chains with this bound which used a vectors.
 * A vector reduced ahead costs about half of a vector reduced alone, hence
 * the largest a reached by half of these chains. */
static uint64_t fglmc_predict(const uint32_t *h, const uint64_t b){
  uint64_t n = 0;
  for(uint64_t a = 1; a <= FGLMC_VLEN; a++) n += h[a];
  if(n < 8) return b;
  uint64_t s = 0;
  for(uint64_t a = FGLMC_VLEN; a > 1; a--){
    s += h[a];
    if(2 * s >= n) return a < b ? a : b;
  }
  return 1;
}

/* Fills Q with the next candidates of H (the unused ones of Q go back to
 * H first). clen[s] is the length of the
 * chain starting at the element s of the staircase: the chain of x_i m is
 * not longer than the one of m, hence at most clen[parent] + 1 vectors
 * (the last one may give a relation). */
static void fglmc_prefetch(fglmc_pre_t *Q, fglmc_heap_t *H, fglmc_lu_t *E,
                           const sp_multimatfglm_t *Mt,
                           const fglm_lex_staircase_t *lex, const CF_t *Vs,
                           const uint32_t *clen,
                           const uint32_t (*chist)[FGLMC_VLEN + 1],
                           const int32_t *last, uint64_t *acc){
  const szmat_t D = Mt->ncols;
  const int32_t nv = Mt->nvars, xn = nv - 1;
  /* candidates not used yet go back to the heap */
  for(uint32_t i = Q->h; i < Q->n; i++){
      fglmc_hpush(H, Q->c[i]);
  }
  fglmc_pre_clear(Q);
  uint32_t K = 0;
  while(H->n && Q->n < FGLMC_PRE){
    const fglmc_cand_t u = H->a[0];
    const int32_t *prev = Q->n ? Q->c[Q->n - 1].e : last;
    if(prev && !memcmp(prev, u.e, nv * sizeof(int32_t))){
      /* same monomial from another parent: tighter bound on the chain */
      fglmc_hpop(H);
      uint32_t *c = Q->len + Q->n - 1;
      if(Q->n && u.parent != UINT32_MAX && *c > clen[u.parent] + 1){
        K -= *c - (clen[u.parent] + 1);
        *c = clen[u.parent] + 1;
      }
      /* (the bound only decreases: the prediction stays an upper bound) */
      continue;
    }
    uint64_t c = 0;
    if(u.parent != UINT32_MAX){
      c = fglmc_kmax(lex, u.e, nv, D);
      if(c > (uint64_t)clen[u.parent] + 1) c = clen[u.parent] + 1;
      if(c > FGLMC_VLEN) c = FGLMC_VLEN;
      if(c) c = fglmc_predict(chist[c], c);
      if(c > FGLMC_PRE) c = FGLMC_PRE;
      /* normal forms are left to the main loop (computed only if needed) */
      if(c && fglmc_rows_missing(Mt, u.var, Vs + (uint64_t)u.parent * D)){
          c = 0;
      }
    }
    if(K + c > FGLMC_PRE && Q->n > 0) break;
    fglmc_hpop(H);
    Q->c[Q->n] = u;
    Q->len[Q->n] = (uint32_t)c;
    Q->n++;
    K += (uint32_t)c;
  }
  /* vectors NF(e x_n^k), k < len */
  K = 0;
  for(uint32_t i = 0; i < Q->n; i++){
    const fglmc_cand_t *u = Q->c + i;
    Q->off[i] = K;
    for(uint32_t k = 0; k < Q->len[i]; k++){
      CF_t *row = Q->V + (uint64_t)(K + k) * D;
      if(k == 0){
        fglmc_mulvec(row, Mt, u->var, Vs + (uint64_t)u->parent * D, acc, NULL);
      } else {
        fglmc_mulvec(row, Mt, xn, row - D, acc, NULL);
      }
    }
    K += Q->len[i];
  }
  Q->r0 = E->r;
  if(K){
    nmod_mat_init(Q->B, K, D, Mt->charac);
    for(uint32_t i = 0; i < K; i++){
      const CF_t *row = Q->V + (uint64_t)i * D;
      mp_limb_t *b = nmod_mat_entry_ptr(Q->B, i, 0);
      for(szmat_t l = 0; l < D; l++) b[l] = row[l];
    }
    fglmc_lu_reduce_range(E, Q->B, Q->C, 0, E->r);
    Q->init = 1;
  }
}

/* If rp != NULL, the dense rows which are not computed yet are computed on
 * demand (only the rows of the matrix of x_n are needed beforehand). */
static int fglm_classic_learn(sp_multimatfglm_t *Mt, const int ord,
                              fglm_lex_staircase_t *lex, fglm_trace_t *tr,
                              const fglmc_rows_t *rp, CF_t **coeffs_out,
                              md_t * st){
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
  uint32_t *clen = malloc((D + 1) * sizeof(uint32_t));
  /* chist[b][a]: chains with bound b (see fglmc_predict) using a vectors */
  uint32_t chist[FGLMC_VLEN + 1][FGLMC_VLEN + 1];
  memset(chist, 0, sizeof(chist));
  fglmc_pre_t Q;
  memset(&Q, 0, sizeof(Q));
  Q.c   = malloc(FGLMC_PRE * sizeof(fglmc_cand_t));
  Q.off = malloc(FGLMC_PRE * sizeof(uint32_t));
  Q.len = malloc(FGLMC_PRE * sizeof(uint32_t));
  Q.V   = malloc((uint64_t)FGLMC_PRE * D * sizeof(CF_t));
  lex->g_lm     = malloc(maxc * nv * sizeof(int32_t));
  lex->g_ntail  = malloc(maxc * sizeof(szmat_t));
  lex->ngens = 0;
  lex->ncoeffs = 0;
  if(tr != NULL){
      tr->s_parent = malloc((D + 1) * sizeof(szmat_t));
      tr->s_var    = malloc((D + 1) * sizeof(int32_t));
      tr->g_parent = malloc(maxc * sizeof(szmat_t));
      tr->g_var    = malloc(maxc * sizeof(int32_t));
      if(!tr->s_parent || !tr->s_var || !tr->g_parent || !tr->g_var){
          ok = 0;
      }
  }
  if(!ok || !Vs || !S || !pool || !acc || !wl || !al || !coeffs || !H.a
     || !lex->g_lm || !lex->g_ntail || !clen || !Q.c || !Q.off || !Q.len
     || !Q.V){
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

  double rt = realtime();
  double ct = cputime();
  lex->npols = calloc(nv, sizeof(uint64_t)); // number of polynomials by variables
  if(lex->npols == NULL){
      ok = 0;
  }
  uint64_t npols = 0;
  if (st->info_level > 1) {
      fprintf(VERBSTREAM, "elimination polynomial                  ");
  }
  int32_t curvar = xn;
  while((H.n || Q.h < Q.n) && ok){
    /* new candidates ahead when Q is empty or when too many rows were
     * inserted since its reduction */
    if(Q.h == Q.n || E.r - Q.r0 >= FGLMC_PRE){
        fglmc_prefetch(&Q, &H, &E, Mt, lex, Vs, clen,
                       (const uint32_t (*)[FGLMC_VLEN + 1])chist, last, acc);
    }
    /* next candidate: the smallest one of Q and H */
    fglmc_cand_t u;
    int64_t qi = -1;
    if(Q.h < Q.n && (H.n == 0 || fglmc_cmp(Q.c[Q.h].e, H.a[0].e, nv, ord) <= 0)){
        qi = Q.h++;
        u = Q.c[qi];
    } else {
        u = fglmc_hpop(&H);
    }
    if(last && !memcmp(last, u.e, nv * sizeof(int32_t))) continue;
    last = u.e;
    int32_t degu = 0;
    for(int32_t i = 0; i < xn; i++) {
        degu += u.e[i];
        if (degu > 1) break;
    }
    if (degu == 1) { // new variable for the chain
        rt = realtime() - rt;
        ct = cputime() - ct;
        if (st->info_level > 1) {
            if (curvar != xn){
                fprintf(VERBSTREAM,
                        "variable %3d                            ",
                        curvar);
            }
            fprintf(VERBSTREAM,
                    "%16lu new           ",npols);
            fprintf (VERBSTREAM, "%15.2f | %-13.2f\n",rt,ct);
            fflush(VERBSTREAM);
        }
        rt = realtime();
        ct = cputime();
        lex->npols[curvar] = npols;
        curvar--;
        npols = 0;
    }
    /* is it a multiple of a lm? if not bound kmax */
    const uint64_t kmax = fglmc_kmax(lex, u.e, nv, D);
    if(kmax == 0) continue;
    const int32_t ng0 = lex->ngens;

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
        int zero;
        if(qi >= 0 && k < Q.len[qi]){
          /* vector reduced ahead against the rows 0..Q.r0-1 */
          const uint64_t ib = Q.off[qi] + k;
          memcpy(row, Q.V + ib * D, D * sizeof(CF_t));
          const mp_limb_t *b = nmod_mat_entry_ptr(Q.B, ib, 0);
          for(szmat_t l = 0; l < D; l++){
              wl[l] = b[l];
          }
          for(slong j = 0; j < Q.r0; j++){
              al[j] = nmod_mat_entry(Q.C, ib, j);
          }
          zero = fglmc_lu_reduce(&E, wl, Q.r0, al + Q.r0);
        } else {
        if(k == 0){
          if(u.parent == UINT32_MAX){
              memset(row, 0, D * sizeof(CF_t));
              row[Mt->one] = 1;
          }
          else {
              if(!fglmc_need_rows(rp, Mt, u.var,
                                  Vs + (uint64_t)u.parent * D,
                                  st->info_level)){
                  ok = 0;
                  break;
              }
              fglmc_mulvec(row, Mt, u.var,
                           Vs + (uint64_t)u.parent * D,
                           acc, NULL);
          }
        } else {
          fglmc_mulvec(row, Mt, xn, Vs + (uint64_t)(t - 1) * D,
                       acc, NULL);
        }
        for(szmat_t l = 0; l < D; l++){
            wl[l] = row[l];
        }
        zero = fglmc_lu_reduce(&E, wl, 0, al);
        }
        if(zero){
            if(!fglmc_push_relation(&E, al, lex, t, &coeffs, &cap_c)){
                ok = 0;
                break;
            }
            if(tr != NULL){
                fglmc_set_parent(tr->g_parent + lex->ngens,
                                 tr->g_var + lex->ngens, &u, k, t, xn);
            }
            fglmc_new_head(lex, u.e, (int32_t)k, nv, t);
            npols++;
            chain_done = 1;
        } else {
          if(t == D){
              ok = 0;
              break;
          }
          fglmc_lu_insert(&E, wl, al);
          if(tr != NULL){
              fglmc_set_parent(tr->s_parent + t, tr->s_var + t, &u, k, t, xn);
          }
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
                       Vs + (uint64_t)(t + i - 1) * D, acc, NULL);
      }
      nmod_mat_t B, C;
      nmod_mat_init(B, K, D, p);
      for(slong i = 0; i < K; i++){
        const CF_t *row = Vs + (uint64_t)(t + i) * D;
        for(szmat_t l = 0; l < D; l++){
            nmod_mat_entry(B, i, l) = row[l];
        }
      }
      /* the block is reduced against the rows of E (one product), then by
       * sub-blocks of FGLMC_KSUB rows: each one against the rows inserted by
       * the previous sub-blocks (one product), then row by row. The first
       * row which reduces to zero ends the chain. */
      const slong r_old = E.r;
      fglmc_lu_reduce_range(&E, B, C, 0, r_old);   /* C : K x r_old */
      slong j0 = K;
      for(slong i0 = 0; i0 < K && ok && j0 == K; i0 += FGLMC_KSUB){
          const slong i1 = (i0 + FGLMC_KSUB < K) ? i0 + FGLMC_KSUB : K;
          const slong rb = E.r;
          nmod_mat_t Bs, Cs;
          nmod_mat_window_init(Bs, B, i0, 0, i1, D);
          fglmc_lu_reduce_range(&E, Bs, Cs, r_old, rb);
          for(slong i = i0; i < i1; i++){
              for(szmat_t l = 0; l < D; l++){
                  wl[l] = nmod_mat_entry(B, i, l);
              }
              for(slong j = 0; j < r_old; j++){
                  al[j] = nmod_mat_entry(C, i, j);
              }
              for(slong j = r_old; j < rb; j++){
                  al[j] = nmod_mat_entry(Cs, i - i0, j - r_old);
              }
              const int zero = fglmc_lu_reduce(&E, wl, rb, al + rb);
              if(!zero){
                  if(t == D){
                      ok = 0;
                      break;
                  }
                  fglmc_lu_insert(&E, wl, al);
                  if(tr != NULL){
                      fglmc_set_parent(tr->s_parent + t, tr->s_var + t, &u,
                                       k + i, t, xn);
                  }
                  fglmc_new_stair(u.e, (int32_t)(k + i), nv, S, &t);
              } else {
                  if(!fglmc_push_relation(&E, al, lex, t, &coeffs, &cap_c)){
                      ok = 0;
                      break;
                  }
                  if(tr != NULL){
                      fglmc_set_parent(tr->g_parent + lex->ngens,
                                       tr->g_var + lex->ngens, &u, k + i, t, xn);
                  }
                  fglmc_new_head(lex, u.e, (int32_t)(k + i), nv, t);
                  npols++;
                  j0 = i;
                  break;
              }
          }
          nmod_mat_clear(Cs);
          nmod_mat_window_clear(Bs);
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
    /* statistics for fglmc_predict */
    if(ok && u.parent != UINT32_MAX){
        uint64_t b = (uint64_t)clen[u.parent] + 1;
        if(b > kmax) b = kmax;
        if(b > FGLMC_VLEN) b = FGLMC_VLEN;
        uint64_t a = (uint64_t)(t - t_head) + (lex->ngens > ng0);
        if(a > FGLMC_VLEN) a = FGLMC_VLEN;
        chist[b][a]++;
    }
    /* monomial m is in the staircase, next ones are x_i m, i < n */
    if(ok && t > t_head){
        clen[t_head] = t - t_head;
        for(int32_t i = 0; i < xn; i++){
            int32_t *e = pool + npool++ * nv;
            memcpy(e, u.e, nv * sizeof(int32_t));
            e[i]++;
            fglmc_hpush(&H, (fglmc_cand_t){ e, t_head, i });
        }
    }
  }
  // largest variable is done
  rt = realtime() - rt;
  ct = cputime() - ct;
  if(st->info_level > 1) {
      fprintf(VERBSTREAM,
              "variable %3d                            ", curvar);
      fprintf(VERBSTREAM,
              "%16lu new           ",npols);
      fprintf (VERBSTREAM, "%15.2f | %-13.2f\n",rt,ct);
      fflush(VERBSTREAM);
  }
  rt = realtime();
  ct = cputime();
  if(lex->npols != NULL){
      lex->npols[curvar] = npols;
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
      lex->nterms = lex->ngens;
      for(uint64_t i = 0; i < lex->ncoeffs; i++){
          if(coeffs[i]){
              lex->nterms++;
          }
      }

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
  fglmc_pre_clear(&Q);
  free(Q.c);
  free(Q.off);
  free(Q.len);
  free(Q.V);
  free(clen);
  if (st->info_level > 1){
      fprintf(VERBSTREAM,
              "-------------------------------------------------\
-----------------------------------------------------\n");
  }
  return ok;
}


/* Replays the trace on a new prime.
 * V (D x D): row k = NF(s_k), computed from the row of its parent;
 * W (ngens x D): row j = NF(g_j);
 * X such that X V = -W gives g_j = lm + sum_k X[j][k] s_k: the transposed
 * system V^T X^T = -W^T is solved. The trace is valid modulo p iff V is
 * invertible and X[j][k] = 0 for k >= g_ntail[j].
 * coeffs (tr->lex.ncoeffs) receives, polynomial after polynomial, the
 * coefficients of s_0, ..., s_{g_ntail[j]-1}.
 * Returns 0 on success, 1 if V is singular, 2 if the support differs, 3 if
 * a dense row not computed on the first prime is needed. */
static int fglm_classic_apply(CF_t *coeffs, const sp_multimatfglm_t *Mt,
                              const fglm_trace_t *tr){
  const szmat_t D = Mt->ncols;
  const CF_t p = Mt->charac;
  const fglm_lex_staircase_t *lex = &tr->lex;
  const int32_t ng = lex->ngens;
  int ret = 0;
  int missing = 0;
  CF_t *Vs = malloc((uint64_t)D * D * sizeof(CF_t));
  CF_t *w  = malloc(D * sizeof(CF_t));
  uint64_t *acc = malloc(D * sizeof(uint64_t));
  nmod_mat_t Vt, Bt, Xt;
  nmod_mat_init(Vt, D, D, p);
  nmod_mat_init(Bt, D, ng, p);
  nmod_mat_init(Xt, D, ng, p);

  for(szmat_t k = 0; k < D; k++){
      CF_t *row = Vs + (uint64_t)k * D;
      if(tr->s_parent[k] == UINT32_MAX){
          memset(row, 0, D * sizeof(CF_t));
          row[Mt->one] = 1;
      } else {
          fglmc_mulvec(row, Mt, tr->s_var[k],
                       Vs + (uint64_t)tr->s_parent[k] * D, acc, &missing);
      }
      for(szmat_t l = 0; l < D; l++){
          nmod_mat_entry(Vt, l, k) = row[l];
      }
  }
  for(int32_t j = 0; j < ng; j++){
      fglmc_mulvec(w, Mt, tr->g_var[j], Vs + (uint64_t)tr->g_parent[j] * D,
                   acc, &missing);
      for(szmat_t l = 0; l < D; l++){
          nmod_mat_entry(Bt, l, j) = w[l] ? p - w[l] : 0;
      }
  }
  if(missing){
      ret = 3;
      goto done;
  }
  if(!nmod_mat_solve(Xt, Vt, Bt)){
      ret = 1;
      goto done;
  }
  uint64_t pos = 0;
  for(int32_t j = 0; j < ng && !ret; j++){
      for(szmat_t k = lex->g_ntail[j]; k < D; k++){
          if(nmod_mat_entry(Xt, k, j)){
              ret = 2;
              break;
          }
      }
      for(szmat_t k = 0; k < lex->g_ntail[j]; k++){
          coeffs[pos++] = (CF_t)nmod_mat_entry(Xt, k, j);
      }
  }
done:
  nmod_mat_clear(Vt);
  nmod_mat_clear(Bt);
  nmod_mat_clear(Xt);
  free(Vs);
  free(w);
  free(acc);
  return ret;
}

/* Reduced lex Gb (sparse form) from the staircase and the dense coefficients
 * given by fglm_classic_learn */
static fglm_gb_t *fglm_gb_from_staircase(const fglm_lex_staircase_t *lex,
                                         const CF_t *coeffs, const CF_t p){
  const int32_t nv = lex->nvars;
  fglm_gb_t *G = calloc(1, sizeof(fglm_gb_t));
  G->charac = p;
  G->nvars = nv;
  G->ngens = lex->ngens;
  G->nterms = lex->nterms;
  G->lens = malloc((lex->ngens + 1) * sizeof(int32_t));
  G->exps = malloc((lex->nterms + 1) * nv * sizeof(int32_t));
  G->cfs  = malloc((lex->nterms + 1) * sizeof(CF_t));
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
static fglm_gb_t *fglm_classic(sp_multimatfglm_t *Mt, const int ord,
                               const fglmc_rows_t *rp, md_t * st){
  fglm_lex_staircase_t lex;
  memset(&lex, 0, sizeof(lex));
  fglm_gb_t *G = NULL;
  CF_t *coeffs = NULL;
  const int old_thr = fglmc_flint_threads_push(Mt->nthreads);
  if(fglm_classic_learn(Mt, ord, &lex, NULL, rp, &coeffs, st)){
      G = fglm_gb_from_staircase(&lex, coeffs, Mt->charac);
  }
  st->fglm_rtime = realtime() - st->fglm_rtime;
  st->fglm_ctime = cputime() - st->fglm_ctime;
  print_fglm_classic_data (VERBSTREAM, st, Mt, /* G,*/ &lex);
  free(coeffs);
  fglmc_flint_threads_pop(old_thr);
  free_fglm_lex_staircase(&lex);
  return G;
}

