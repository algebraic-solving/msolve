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

#include<flint/fmpz.h>
#include<flint/fmpz_vec.h>
#include<flint/fmpz_mat.h>
#include<flint/fmpz_lll.h>
#include "config.h"
#include <inttypes.h>
#include "streams.h"
#include "../neogb/gmp_u64.h"

#define NBCHECK 2
/* number of bits of slack kept in rational reconstruction bounds */
#define RR_MARGIN 16
/* witnesses are reconstructed every rr primes, rr being the largest
   power of 2 <= nprimes / RR_PERIOD_DIV */
#define RR_PERIOD_DIV 128
/* Witness coefficients of a polynomial reconstructed simultaneously by LLL
   (see lll_witnesses): at most LLL_WIT_K of them, attempts are made once
   at least LLL_WIT_START primes are used, then each time the number of
   primes has grown by LLL_WIT_GROWTH percent */
#define LLL_WIT_K 8
#define LLL_WIT_START 8
#define LLL_WIT_GROWTH 5
/* number of witnesses when no denominator is guessed: LLL is much more
   expensive on these lattices, whose first column is not scaled */
#define LLL_WIT_K1 4
/* LLL without guess is only tried when, for some lifted polynomial, at
   least this share of the bits of the denominator of the previous one does
   not divide its own */
#define LLL_WIT_INFL 0.125
/* the estimated cost of the calls to LLL is kept below LLL_WIT_BUDGET
   percent of the estimated cost of the primes used */
#define LLL_WIT_BUDGET 10
/* bits kept below the guessed denominator in the lattices of
   lll_common_den */
#define LLL_TRUNC_SLACK 64
#ifndef MIN
#define MIN(x, y) ((x) > (y) ? (y) : (x))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

typedef struct{
  uint32_t len; /* length of the encoded polynomial */
  /* Coefficients modulo several primes, stored as one contiguous block of
   * len * alloc entries in prime-major order: the coefficient of index j
   * modulo the k-th prime is cf_32[k * len + j] (use MODCF below).
   * Adding primes only grows the block at its end, and each new prime
   * fills one contiguous row.
   * Warning: this is free-ed dynamically (and set to NULL) as soon as the
   * polynomial has been lifted and checked. */
  uint32_t *cf_32;
  /* same, used instead of cf_32 when the primes have more than 32 bits */
  uint64_t *cf_u64;
  /* len integer numerators of the lifted coefficients; they all share the
   * common denominator lm (hence the polynomial is lm * LM + sum cf_qq[j]*m_j
   * up to the order of the terms). */
  mpz_t *cf_qq;
  mpz_t lm; /* stores the leading coefficient (hence all numerators of cf_qq should be divided by lm) */
} modpolys_struct;

typedef modpolys_struct modpolys_t[1];

/* coefficient of index j of polynomial p modulo the k-th prime */
#define MODCF_IDX(p, j, k) ((size_t)(k) * (p)->len + (j))
#define MODCF(p, j, k) ((p)->cf_u64 != NULL                         \
                        ? (p)->cf_u64[MODCF_IDX(p, j, k)]            \
                        : (uint64_t)(p)->cf_32[MODCF_IDX(p, j, k)])
/* nonzero while the modular images of polynomial p are stored */
#define HAS_MODCF(p) ((p)->cf_32 != NULL || (p)->cf_u64 != NULL)

typedef struct {
  uint32_t alloc; /* alloc -> max number of primes */
  uint32_t nprimes; /* number of primes */
  mp_limb_t *primes; /* array of prime numbers encoded with uint64_t to ensure
                       compatibility with flint */
  mp_limb_t *pinvs; /* their inverses for n_mulmod2_preinv */
  int wide; /* 1 if the primes have more than 32 bits (images in cf_u64) */
  mp_limb_t *cf_64; /* array of length equal to number of primes which will be used
                    to copy coefficients (hence ensuring compatibility with
                    flint) */
  uint32_t ld; /* number of polynomials */
  int nv; /* number of variables */
  int32_t *ldm; /* lead monomials */
  ht_t *bht; /* hash table */
  hm_t **hm; /* hashed monomials representing exponents */
  bl_t *lmps; /* position of non-redundant lead monomials in basis */
  modpolys_t *modpolys; /* array of polynomials modulo primes */
} gb_modpoly_array_struct;

typedef gb_modpoly_array_struct gb_modpoly_t[1];

#define NEWGBLIFT 1

/* F4 computes with 64 bit integers modulo primes of at least 31 bits (see
   set_ff_bits), the modular images are then stored in cf_u64 */
#define GB_PRIME_IS_WIDE(p) ((uint64_t)(p) >= ((uint64_t)1 << 31))


typedef struct{
  int32_t npol; /* number of polynomials to be lifted */
  int32_t rr; /* number of primes before activating rational reconstruction */
  int64_t nall; /* number of primes whose trace has been applied, bad ones
                   included (see witness_attempt_due) */
  int32_t nsteps; /* number of steps for lifting GB (per degree) */
  int32_t *steps; /* array of length nsteps ; the sum of the entries should
                     equal npol */
  int32_t cstep; /* current step lifting GB */
  /* liftings are performed on ranges of polynomials (depending on their degrees) */
  int32_t lstart; /* index of first polynomial to be lifted */
  int32_t lend; /* index of last polynomial to be lifted */

  int crt_mult; /* indicates if multi-mod flint structures need to be
                initialized */
  mpz_t *crt; /* stores current CRT */
  int recon; /* equals 1 when some rational number can be lifted, else 0 */
  int32_t *coef; /* array of indices to lift */
  mpz_t *num; /* lifted numerator */
  mpz_t *den; /* lifted denominator */
  mpz_t gden; /* guessed denominator */
  mpz_t tmp;
  int32_t start; /* smallest index of poly whose witness coef has been lifted but not checked*/
  int32_t end; /* largest index of poly whose witness coef has been lifted but not checked*/
  int *check1; /* tells whether lifted data are ok with one more prime */
  int *check2; /* tells whether lifted data are ok with two more primes */
  int32_t S;
  int32_t V; /* the CRTs of the witness coefficients are up to date for the
                polynomials of index S <= k < V (see witness_crt_window) */
  int32_t *coef2; /* second witness coefficient */
  mpz_t *crt2;
  mpz_t *num2;
  mpz_t *den2;
  mpz_t *cden; /* lcm of the witness denominators */
  mpz_t pgden; /* denominator guessed from the previous polynomial */
  int32_t *coef3; /* third witness coefficient */
  mpz_t *crt3;
  mpz_t *num3;
  mpz_t *den3;
  int32_t *lwit; /* indices of the LLL witnesses, LLL_WIT_K per polynomial */
  int *nlwit; /* number of LLL witnesses of each polynomial */
  int nthreads; /* number of threads for the CRT of the coefficients */
  int32_t lll_next; /* number of primes of the next LLL attempt */
  int32_t lll_fnext; /* same for the full lifting (see lll_full_block) */
  int32_t lll_ntries; /* number of calls to LLL */
  int32_t lll_nsucc; /* number of candidates found by LLL */
  double lll_time; /* elapsed time spent in lll_witnesses */
  double lll_cost; /* estimated cost of the calls to LLL (see
                      lll_witnesses_cost) */
  double infl; /* largest share of the bits of the denominator of a lifted
                  polynomial which do not divide the denominator of the
                  next one (see lll_update_inflation) */
} data_lift_struct;

typedef data_lift_struct data_lift_t[1];


static inline void data_lift_init(data_lift_t dl,
                                  int32_t npol,
                                  int32_t *steps, int32_t nsteps){
  dl->npol = npol;
  dl->rr = 1;
  dl->nall = 0;
  dl->lstart = 0;
  dl->nsteps = nsteps;
  dl->S = 0;
  dl->V = 0;
  mpz_init_set_ui(dl->pgden, 1);
  int32_t i;

  dl->steps = (int32_t *)calloc(nsteps, sizeof(int32_t));
  for(i = 0; i < nsteps; i++){
    dl->steps[i] = steps[i];
  }
  dl->cstep = 0;
  dl->lend = npol;
  dl->crt_mult = 0;
  dl->crt = (mpz_t *)malloc(sizeof(mpz_t) * dl->npol);
  for(int32_t i = 0; i < dl->npol; i++){
    mpz_init(dl->crt[i]);
  }
  dl->recon = 0;
  dl->coef =(int32_t *) calloc(npol, sizeof(int32_t) );

  dl->num = malloc(sizeof(mpz_t) * npol);
  for(i = 0; i < npol; i++){
    mpz_init(dl->num[i]);
  }
  dl->den = malloc(sizeof(mpz_t) * npol);
  for(i = 0; i < npol; i++){
    mpz_init(dl->den[i]);
  }
  mpz_init_set_ui(dl->gden, 1);
  mpz_init(dl->tmp);
  dl->start = 0;
  dl->end = npol;
  dl->check1 = calloc(npol, sizeof(int));
  dl->check2 = calloc(npol, sizeof(int));
  dl->coef2 = calloc(npol, sizeof(int32_t));
  dl->crt2 = malloc(sizeof(mpz_t) * npol);
  dl->num2 = malloc(sizeof(mpz_t) * npol);
  dl->den2 = malloc(sizeof(mpz_t) * npol);
  dl->cden = malloc(sizeof(mpz_t) * npol);
  dl->coef3 = calloc(npol, sizeof(int32_t));
  dl->crt3 = malloc(sizeof(mpz_t) * npol);
  dl->num3 = malloc(sizeof(mpz_t) * npol);
  dl->den3 = malloc(sizeof(mpz_t) * npol);
  for(i = 0; i < npol; i++){
    mpz_init(dl->crt2[i]);
    mpz_init(dl->num2[i]);
    mpz_init(dl->den2[i]);
    mpz_init(dl->cden[i]);
    mpz_init(dl->crt3[i]);
    mpz_init(dl->num3[i]);
    mpz_init(dl->den3[i]);
  }
  dl->lwit = calloc((size_t)npol * LLL_WIT_K, sizeof(int32_t));
  dl->nlwit = calloc(npol, sizeof(int));
  dl->nthreads = 1;
  dl->lll_next = LLL_WIT_START;
  dl->lll_fnext = LLL_WIT_START;
  dl->lll_ntries = 0;
  dl->lll_nsucc = 0;
  dl->lll_time = 0;
  dl->lll_cost = 0;
  dl->infl = 0;
}

static inline void data_lift_clear(data_lift_t dl){
  for(int32_t i = 0; i < dl->npol; i++){
    mpz_clear(dl->crt[i]);
  }
  free(dl->crt);

  free(dl->steps);
  free(dl->coef);

  for(int32_t i = 0; i < dl->npol; i++){
    mpz_clear(dl->num[i]);
  }
  free(dl->num);

  for(int32_t i = 0; i < dl->npol; i++){
    mpz_clear(dl->den[i]);
  }
  free(dl->den);

  mpz_clear(dl->gden);
  mpz_clear(dl->tmp);
  mpz_clear(dl->pgden);
  free(dl->check1);
  free(dl->check2);
  for(int32_t i = 0; i < dl->npol; i++){
    mpz_clear(dl->crt2[i]);
    mpz_clear(dl->num2[i]);
    mpz_clear(dl->den2[i]);
    mpz_clear(dl->cden[i]);
    mpz_clear(dl->crt3[i]);
    mpz_clear(dl->num3[i]);
    mpz_clear(dl->den3[i]);
  }
  free(dl->crt2);
  free(dl->num2);
  free(dl->den2);
  free(dl->cden);
  free(dl->coef2);
  free(dl->crt3);
  free(dl->num3);
  free(dl->den3);
  free(dl->coef3);
  free(dl->lwit);
  free(dl->nlwit);

}


static inline void gb_modpoly_init(gb_modpoly_t modgbs,
                                   uint32_t alloc, int32_t *lens,
                                   bs_t *bs,
                                   int nv, uint32_t ld,
                                   int32_t *lm,
                                   md_t *st, const int wide){
  modgbs->alloc = alloc;
  modgbs->nprimes = 0;
  modgbs->wide = wide;
  modgbs->primes = (mp_limb_t *)calloc(alloc, sizeof(mp_limb_t));
  modgbs->pinvs = (mp_limb_t *)calloc(alloc, sizeof(mp_limb_t));
  modgbs->cf_64 = (mp_limb_t *)calloc(alloc, sizeof(mp_limb_t));
  modgbs->ld = ld;
  modgbs->nv = nv;
  modgbs->modpolys = (modpolys_t *)malloc(sizeof(modpolys_t) * ld);

  modgbs->ldm = (int32_t *)calloc(nv*ld, sizeof(int32_t));
  modgbs->bht = bs->ht;
  modgbs->hm = (hm_t **)malloc(modgbs->ld * sizeof(hm_t *));
  modgbs->lmps = (bl_t *)calloc(modgbs->ld, sizeof(bl_t));

  len_t idx;
  for(len_t i = 0; i < modgbs->ld; i++){
      idx = bs->lmps[i];
      modgbs->hm[i] = malloc((1+lens[i] + OFFSET)*sizeof(hm_t));
      for(len_t j = 0; j<lens[i]+1; j++){
          modgbs->hm[i][j+OFFSET] = bs->hm[idx][j+OFFSET];
      }
  }
  for(bl_t i = 0; i < modgbs->ld; i++){
      modgbs->lmps[i] = bs->lmps[i];
  }
  for(int32_t i = 0; i < ld; i++){
    for(int j = 0; j < nv; j++){
      modgbs->ldm[i*nv+j] = lm[i*nv+j];
    }
  }
  for(uint32_t i = 0; i < ld; i++){
    modgbs->modpolys[i]->len = lens[i];
    modgbs->modpolys[i]->cf_32 = NULL;
    modgbs->modpolys[i]->cf_u64 = NULL;
    modgbs->modpolys[i]->cf_qq = NULL;
    if(lens[i] > 0){
      if(wide){
        modgbs->modpolys[i]->cf_u64 = (uint64_t *)calloc((size_t)lens[i] * alloc,
                                                         sizeof(uint64_t));
      }
      else{
        modgbs->modpolys[i]->cf_32 = (uint32_t *)calloc((size_t)lens[i] * alloc,
                                                        sizeof(uint32_t));
      }
      modgbs->modpolys[i]->cf_qq = (mpz_t *)malloc(sizeof(mpz_t)*lens[i]);
    }
    for(uint32_t j = 0; j < lens[i]; j++){
      mpz_init(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_init(modgbs->modpolys[i]->lm);
    mpz_set_ui(modgbs->modpolys[i]->lm, 1);
  }
}


/* Makes room for newalloc more primes. Polynomials whose modular data have
 * already been free-ed (see HAS_MODCF) are skipped. */
static inline void gb_modpoly_realloc(gb_modpoly_t modgbs,
                                      uint32_t newalloc){

  uint32_t oldalloc = modgbs->alloc;
  modgbs->alloc += newalloc;

  mp_limb_t *newprimes = (mp_limb_t *)realloc(modgbs->primes,
                                            modgbs->alloc * sizeof(mp_limb_t));

  if(newprimes == NULL){
    fprintf(ERRSTREAM, "Problem when reallocating modgbs (primes)\n");
    exit(1);
  }
  modgbs->primes = newprimes;
  for(uint32_t i = oldalloc; i < modgbs->alloc; i++){
    modgbs->primes[i] = 0;
  }
  mp_limb_t *newpinvs = (mp_limb_t *)realloc(modgbs->pinvs,
                                           modgbs->alloc * sizeof(mp_limb_t));
  if(newpinvs == NULL){
    fprintf(ERRSTREAM, "Problem when reallocating modgbs (primes)\n");
    exit(1);
  }
  modgbs->pinvs = newpinvs;

  mp_limb_t *ncf_64 = (mp_limb_t *)realloc(modgbs->cf_64,
                                         modgbs->alloc * sizeof(mp_limb_t));
  if(ncf_64 == NULL){
    fprintf(ERRSTREAM, "Problem when reallocating modgbs (cfs)\n");
    exit(1);
  }
  modgbs->cf_64 = ncf_64;
  for(uint32_t i = oldalloc; i < modgbs->alloc; i++){
    modgbs->cf_64[i] = 0;
  }
  for(uint32_t i = 0; i < modgbs->ld; i++){
    const size_t len = modgbs->modpolys[i]->len;
    if(!HAS_MODCF(modgbs->modpolys[i]) || len == 0){
      continue;
    }
    const size_t sz = modgbs->wide ? sizeof(uint64_t) : sizeof(uint32_t);
    void *old = modgbs->wide ? (void *)modgbs->modpolys[i]->cf_u64
      : (void *)modgbs->modpolys[i]->cf_32;
    char *newcfs_pol = (char *)realloc(old, len * modgbs->alloc * sz);
    if(newcfs_pol == NULL){
      fprintf(ERRSTREAM, "Problem when reallocating modgbs (cfs_pol)\n");
      exit(1);
    }
    /* prime-major storage: the new primes are appended at the end */
    memset(newcfs_pol + len * oldalloc * sz, 0,
           len * (modgbs->alloc - oldalloc) * sz);
    if(modgbs->wide){
      modgbs->modpolys[i]->cf_u64 = (uint64_t *)newcfs_pol;
    }
    else{
      modgbs->modpolys[i]->cf_32 = (uint32_t *)newcfs_pol;
    }
  }
}

/* stores the modular images on 64 bits, so that primes of more than 32
   bits can be used for the next images */
static inline void gb_modpoly_widen(gb_modpoly_t modgbs){
  if(modgbs->wide){
    return;
  }
  for(uint32_t i = 0; i < modgbs->ld; i++){
    modpolys_t *pol = modgbs->modpolys + i;
    if(pol[0]->cf_32 == NULL){
      continue;
    }
    const size_t sz = (size_t)pol[0]->len * modgbs->alloc;
    pol[0]->cf_u64 = (uint64_t *)malloc(sz * sizeof(uint64_t));
    for(size_t j = 0; j < sz; j++){
      pol[0]->cf_u64[j] = pol[0]->cf_32[j];
    }
    free(pol[0]->cf_32);
    pol[0]->cf_32 = NULL;
  }
  modgbs->wide = 1;
}

static inline void display_gbmodpoly_cf_32(FILE *file,
                                     gb_modpoly_t modgbs){
  fprintf(file, "alloc = %d\n", modgbs->alloc);
  fprintf(file, "nprimes = %d\n", modgbs->nprimes);
  fprintf(file, "primes = [");
  for(uint32_t i = 0; i < modgbs->alloc-1; i++){
    fprintf(file, "%" PRIu64 ", ", (uint64_t)modgbs->primes[i]);
  }
  fprintf(file, "%" PRIu64 "]\n", (uint64_t)modgbs->primes[modgbs->alloc -1]);
  fprintf(file, "numpolys = %d\n", modgbs->ld);
  fprintf(file, "[\n");
  for(uint32_t i = 0; i < modgbs->ld; i++){
    uint32_t len = modgbs->modpolys[i]->len;
    fprintf(file, "[%d, ", len);
    if(!HAS_MODCF(modgbs->modpolys[i])){
      fprintf(file, "(lifted)],\n");
      continue;
    }
    for(uint32_t j = 0; j < len; j++){
      fprintf(file, "[");
      for(uint32_t k = 0; k < modgbs->alloc-1; k++){
        fprintf(file, "%" PRIu64 ", ", (uint64_t)MODCF(modgbs->modpolys[i], j, k));
      }
      if(j < len - 1){
        fprintf(file, "%" PRIu64 "], ", (uint64_t)MODCF(modgbs->modpolys[i], j, modgbs->alloc-1));
      }
      else{
        fprintf(file, "%" PRIu64 "]\n", (uint64_t)MODCF(modgbs->modpolys[i], j, modgbs->alloc-1));
      }
    }
    fprintf(file, "],\n");
  }
  fprintf(file, "]:\n");
}

/* decimal strings of the coefficients cf_qq[0], .., cf_qq[len-1] of the
 * pos-th polynomial and of its leading coefficient lm (index len),
 * converted by nthreads threads: the conversion is most of the time spent
 * writing large bases */
static char **modpoly_strings(gb_modpoly_t modgbs, const int32_t pos,
                              const int nthreads){
  modpolys_t *pol = modgbs->modpolys + pos;
  const int64_t len = pol[0]->len;
  char **str = (char **)malloc(sizeof(char *) * (len + 1));
  /* FLINT's conversion, as mpz_out_str_10 (see iofiles.c) */
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 16)
  for(int64_t i = 0; i <= len; i++){
    fmpz_t f;
    fmpz_init_set_readonly(f, i < len ? pol[0]->cf_qq[i] : pol[0]->lm);
    str[i] = fmpz_get_str(NULL, 10, f);
    fmpz_clear_readonly(f);
  }
  return str;
}

static void modpoly_strings_clear(char **str, const int64_t len){
  for(int64_t i = 0; i <= len; i++){
    flint_free(str[i]);
  }
  free(str);
}

/* str holds the decimal strings of the coefficients (see modpoly_strings) */
static inline void display_modpoly(FILE *file,
                                   gb_modpoly_t modgbs,
                                   int32_t pos,
                                   data_gens_ff_t *gens,
                                   char **str){

  if(modgbs->modpolys[pos]->len == 0){
    display_monomial_single(file, gens, pos, &modgbs->ldm);
    return;
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->lm, 1) != 0){
    fputs(str[modgbs->modpolys[pos]->len], file);
    fprintf(file, "*");
  }
  display_monomial_single(file, gens, pos, &modgbs->ldm);

  len_t i, k;
  hm_t *hm    = NULL;
  /* modgbs->hm[pos] is the copy of the monomials of the pos-th element
   * (see gb_modpoly_init), hence it is indexed by pos, not by lmps[pos] */
  hm  = modgbs->hm[pos]+OFFSET;

  ht_t *ht = modgbs->bht;
  const len_t nv  = ht->nv;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  int *evi    =   (int *)malloc((unsigned long)(ht->nv + 1) * sizeof(int));
  if (ebl == 0) {
    for (i = 1; i < evl; ++i) {
      evi[i-1]    =   i;
    }
  } else {
    for (i = 1; i < ebl; ++i) {
      evi[i-1]    =   i;
    }
    for (i = ebl+1; i < evl; ++i) {
      evi[i-2]    =   i;
    }
  }
  len_t len = modgbs->modpolys[pos]->len + 1;
  /* once lifted, all coefficients are integers (common denominator is
   * stored in lm), cf_qq[i] is the coefficient of index i */
  mpz_t *cf = modgbs->modpolys[pos]->cf_qq;
  for(i = modgbs->modpolys[pos]->len-1; i > 0 ; i--){
    if(mpz_cmp_ui(cf[i], 0)>0){
      fprintf(file, "+");
    }
    //coef is -1
    if(mpz_cmp_si(cf[i], -1) == 0){
        fprintf(file, "-");
    }
    int star = 0;
    //coef is neither 1 nor -1
    if(mpz_cmp_ui(cf[i], 1) != 0 && mpz_cmp_si(cf[i], -1) != 0){
        fputs(str[i], file);
        fprintf(file, "*");
    }
    if(mpz_cmp_ui(cf[i], 0) != 0){
      for (k = 0; k < nv; ++k) {
         if (ht->ev[hm[len-i-1]][evi[k]] == 1) {
            if(star == 1){
                fprintf(file, "*");
            }
            fprintf(file, "%s",gens->vnames[k]);
            star = 1;
         }
         if (ht->ev[hm[len-i-1]][evi[k]] > 1) {
            if(star == 1){
                fprintf(file, "*");
            }
            fprintf(file, "%s^%u",gens->vnames[k], ht->ev[hm[len-i-1]][evi[k]]);
            star = 1;
         }
     }
    }
    fflush(file);
  }
  if(mpz_cmp_ui(cf[0], 0) > 0){
    fprintf(file, "+");
  }
  int deg = 0;
  for (k = 0; k < nv; ++k) {
     deg += ht->ev[hm[len-1]][evi[k]];
  }
  if(deg == 0){
    fputs(str[0], file);
  }
  else{
    if(mpz_cmp_si(cf[0], -1) == 0){
        fprintf(file, "-");
    }
    int star = 0;
    //coef is neither 1 nor -1
    if(mpz_cmp_ui(cf[0], 1) != 0 && mpz_cmp_si(cf[0], -1) != 0){
        fputs(str[0], file);
        fprintf(file, "*");
    }
    fflush(file);
    if(mpz_cmp_ui(cf[0], 0) != 0){
      for (k = 0; k < nv; ++k) {
         if (ht->ev[hm[len-1]][evi[k]] == 1) {
            if(star == 1){
                fprintf(file, "*");
            }
            fprintf(file, "%s",gens->vnames[k]);
            star = 1;
         }
         if (ht->ev[hm[len-1]][evi[k]] > 1) {
            if(star == 1){
                fprintf(file, "*");
            }
            fprintf(file, "%s^%u",gens->vnames[k], ht->ev[hm[len-1]][evi[k]]);
            star = 1;
         }
     }
    }
  }

  free(evi);
}

static inline void display_gbmodpoly_cf_qq(FILE *file,
                                           gb_modpoly_t modgbs,
                                           data_gens_ff_t *gens,
                                           const int nthreads){
  int32_t p = modgbs->ld ;
  if(p==0){
    fprintf(file, "[0]:\n");
    return;
  }
  fprintf(file, "[");
  for(int i = 0; i < p; i++){
    char **str = modpoly_strings(modgbs, i, nthreads);
    display_modpoly(file, modgbs, i, gens, str);
    modpoly_strings_clear(str, modgbs->modpolys[i]->len);
    fprintf(file, i < p - 1 ? ", \n" : "]:\n");
  }
}

static inline void display_lm_gbmodpoly_cf_qq(FILE *file,
                                              gb_modpoly_t modgbs,
                                              data_gens_ff_t *gens){
  int32_t p = modgbs->ld ;

  if(p==0){
    fprintf(file, "[0]:\n");
    return;
  }

  fprintf(file, "[");
  for(int i = 0; i < p-1; i++){
    display_monomial_single(file, gens, i, &modgbs->ldm);
    fprintf(file, ", \n");
  }
  display_monomial_single(file, gens, p-1, &modgbs->ldm);
  fprintf(file, "]:\n");
}

/* frees everything but the hash table modgbs->bht */
static inline void gb_modpoly_without_hash_table_clear(gb_modpoly_t modgbs){
  free(modgbs->primes);
  modgbs->primes = NULL;
  free(modgbs->pinvs);
  modgbs->pinvs = NULL;
  free(modgbs->cf_64);
  modgbs->cf_64 = NULL;
  free(modgbs->ldm);
  modgbs->ldm = NULL;
  for(uint32_t i = 0; i < modgbs->ld; i++){
    /* cf_32 and cf_u64 are NULL if already free-ed after lifting */
    free(modgbs->modpolys[i]->cf_32);
    modgbs->modpolys[i]->cf_32 = NULL;
    free(modgbs->modpolys[i]->cf_u64);
    modgbs->modpolys[i]->cf_u64 = NULL;
    for(uint32_t j = 0; j < modgbs->modpolys[i]->len; j++){
      mpz_clear(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_clear(modgbs->modpolys[i]->lm);
    free(modgbs->modpolys[i]->cf_qq);
    modgbs->modpolys[i]->cf_qq = NULL;
    free(modgbs->hm[i]);
  }
  free(modgbs->lmps);
  modgbs->lmps = NULL;
  free(modgbs->hm);
  modgbs->hm = NULL;
  free(modgbs->modpolys);
  modgbs->modpolys = NULL;
}

static inline void gb_modpoly_clear(gb_modpoly_t modgbs){
  gb_modpoly_without_hash_table_clear(modgbs);
  full_free_hash_table(&modgbs->bht);
}

static inline int32_t degree(int32_t *mon, int nv){
  int32_t deg = 0;
  for(int i = 0; i < nv; i++){
    deg += mon[i];
  }
  return deg;
}

/*
  bexp_lm is a list of monomials involving nv variables of increasing degrees
  return an array of length *nb containing at the i-th position the number of
  monomials of the i-th degree (hence nb is the number of different degrees).
 */
static inline int32_t *array_nbdegrees(int32_t *bexp_lm, int len,
                                     int nv, int * nb){
  *nb = 1;
  int32_t deg = degree(bexp_lm, nv);

  for(int32_t i = 1; i < len; i++){
    int32_t newdeg = degree(bexp_lm + i * nv, nv);

    if(deg != newdeg){
      (*nb)++;
      deg = newdeg;
    }
  }
  int32_t *ldeg = calloc(*nb, sizeof(int32_t));
  deg = degree(bexp_lm, nv);
  ldeg[0] = 1;
  int32_t i = 0, j = 1;
  while(j < len){
    int32_t newdeg = degree(bexp_lm + j * nv, nv);
    if(deg == newdeg){
      ldeg[i]++;
    }
    else{
      i++;
      ldeg[i] = 1;
      deg = newdeg;
    }
    j++;
  }
  return ldeg;
}

static inline int grevlex_is_less_than(int nv, int32_t* m1, int32_t *m2){
  int32_t deg1 = 0, deg2 = 0;
  for(int i = 0; i < nv; i++){
    deg1 += m1[i];
  }
  for(int i = 0; i < nv; i++){
    deg2 += m2[i];
  }
  if(deg1 < deg2){
    return 1;
  }
  if(deg1 > deg2){
    return 0;
  }
  for(int i = 0; i < nv; i++){
    if(m1[i] < m2[i]){
      return 0;
    }
  }
  return 1;
}

static inline int32_t compute_length(int32_t *mon, int nv,
                                     int32_t *basis, int dquot){

  for(int i = dquot - 1; i >= 0; i--){
    if(!grevlex_is_less_than(nv, mon, basis + i * nv)){
      return i + 1;
    }
  }
  return 0;
}

/*
 * bexp_lm is an array of len  monomials with nv variables
 * basis is an array of dquot monomials with nv variables
 * (basis of some quotient defined with bexp_lm, all with grevlex order)
 * basis is sorted increasingly, as for bexp_lm
 *
 * returns an array of length len containing the maximum possible lengths
 * of polynomials with leading terms in bexp_lm (length excluding the
 * leading term).
 */
static inline int32_t *array_of_lengths(int32_t *bexp_lm, len_t len,
        const bs_t * const bs, int nv){
  int32_t *lens = calloc(len, sizeof(int32_t));
  for(len_t i = 0; i < len; i++){
    len_t idx = bs->lmps[i];
    lens[i] = bs->hm[idx][LENGTH]-1;
  }
  return lens;
}

/* returns 0 in case of failure else returns 1 */
static inline int modpgbs_set(gb_modpoly_t modgbs,
                              const bs_t *bs, const ht_t * const ht,
                              const uint64_t fc,
                              int32_t start, const long elim){
  if(modgbs->nprimes >= modgbs->alloc-1){
    fprintf(ERRSTREAM, "Not enough space in modgbs\n");
    exit(1);
  }
  modgbs->primes[modgbs->nprimes] = fc;
  modgbs->pinvs[modgbs->nprimes] = n_preinvert_limb(fc);

  len_t i, j, idx;

  /*****************************************************
   * It should be checked if modbgs->ht needs to be updated as well as
   * modgbs->lmps
   * *****************************************************/

  len_t len   = 0;
  /* hm_t *hm    = NULL; */

  for(i = start; i < modgbs->ld; i++){
    idx = bs->lmps[i];
    if (bs->hm[idx] == NULL) {
      fprintf(ERRSTREAM, " poly is 0\n");
      exit(1);
    } else {
      /* hm  = bs->hm[idx]+OFFSET; */
      len = bs->hm[idx][LENGTH];
    }
    if(!HAS_MODCF(modgbs->modpolys[i])){
      /* already lifted, modular data are no longer needed */
      continue;
    }
    int32_t bc = modgbs->modpolys[i]->len - 1;
    const size_t off = MODCF_IDX(modgbs->modpolys[i], 0, modgbs->nprimes);
    if(modgbs->wide){
      const uint64_t * const cfs = bs->cf_64[bs->hm[idx][COEFFS]];
      uint64_t *dst = modgbs->modpolys[i]->cf_u64 + off;
      for (j = 1; j < len; ++j) {
        dst[bc] = cfs[j];
        bc--;
      }
    }
    else{
      const uint32_t * const cfs = bs->cf_32[bs->hm[idx][COEFFS]];
      uint32_t *dst = modgbs->modpolys[i]->cf_32 + off;
      for (j = 1; j < len; ++j) {
        dst[bc] = cfs[j];
        bc--;
      }
    }
  }

  modgbs->nprimes++;

  return 1;
}

static inline int32_t maxbitsize_gens(data_gens_ff_t *gens, len_t ngens){
  if(gens->field_char != 0){
    return -1;
  }
  const int32_t *lens = gens->lens;
  mpz_t **cfs = gens->mpz_cfs;
  int32_t off = 0;
  int32_t mbs = 0;
  for(int32_t i = 0; i < ngens; i++){
    for(int32_t j = off; j < off + lens[i]; ++j){
      mbs = MAX(mbs,
                mpz_sizeinbase(*(cfs[2*j]), 2) + mpz_sizeinbase(*(cfs[2*j+1]), 2));
    }
    off += lens[i];
  }
  return mbs;
}

static inline int32_t maxbitsize_generators(bs_t *bs){
    int32_t max = 0;
    int32_t ngens = bs->ld;
    for(int32_t i = 0; i < ngens; i++){
        int32_t len = bs->hm[i][LENGTH];
        mpz_t *cfs = bs->cf_qq[bs->hm[i][COEFFS]];
        for(int32_t j = 0; j < len; j++){
            max = MAX(max,
                    mpz_sizeinbase(cfs[j], 2)
                    );
        }
    }
    return max;
}

/* returns the index of the largest polynomial in the basis of the elimination
 * ideal */
static inline int32_t compute_num_gb(int32_t *bexp_lm, int32_t len, int nv, int nev){
  if(nev){
    for(int32_t i = 0; i < len; i++){
      for(int32_t j = 0; j < nev; j++){
        if(bexp_lm[i*nv+j]){
          return i;
        }
      }
    }
    return len;
  }
  return len;
}

static int32_t gb_modular_trace_learning(gb_modpoly_t modgbs,
                                           int32_t *mgb,
                                           int32_t *num_gb,
                                           int32_t **leadmons,
                                           trace_t *trace,
                                           bs_t *bs_qq,
                                           md_t *st,
                                           const uint64_t fc,
                                           int info_level,
                                           int print_gb,
                                           int truncate_lifting,
                                           int32_t start,
                                           int32_t maxbitsize,
                                           int *success)
{
    double ca0/*, rt*/;
    ca0 = realtime();

    bs_t *bs = NULL;
    int32_t err     = 0;
    st->f4_qq_round = 1;
    bs = core_gba(bs_qq, st, &err, fc);
    if (err) {
      fprintf(ERRSTREAM,"Problem with F4, stopped computation.\n");
      exit(1);
    }

    /* rt = realtime()-ca0; */
    st->learning_rtime = realtime()-ca0;

    ht_t *bht = bs->ht;

    /* Leading monomials from Grobner basis */
    int32_t *bexp_lm = get_lm_from_bs(bs, bht);
    /* may be non NULL if the computation has been restarted */
    free(leadmons[0]);
    leadmons[0] = bexp_lm;

    int32_t len = bs->lml;
    num_gb[0] = compute_num_gb(bexp_lm, len, bht->nv, st->nev);
    int32_t *bexp_lm2 = NULL;
    if(st->nev){
      bexp_lm2 = calloc(num_gb[0]*(bht->nv - st->nev), sizeof(int32_t));
      for(int32_t i = 0; i < num_gb[0]; i++){
        for(int j = 0; j < bht->nv - st->nev; j++){
          bexp_lm2[i*(bht->nv - st->nev) + j] = bexp_lm[i*bht->nv + st->nev + j];
        }
      }
      leadmons[0] = bexp_lm2;
    }
    /************************************************/

    /**
    * int32_t *lmb;
    * if(st->nev){
    *   lmb = NULL; //monomial_basis_enlarged(num_gb[0], bht->nv - st->nev,  bexp_lm2, &dquot);
    * }
    * else{
    *   lmb = NULL; //monomial_basis_enlarged(num_gb[0], bht->nv,  bexp_lm, &dquot);
    * }
    */

    /************************************************/

    int32_t *lens = array_of_lengths(leadmons[0], num_gb[0], bs, bht->nv - st->nev);

    if(truncate_lifting != 0 && truncate_lifting < num_gb[0]){
      gb_modpoly_init(modgbs, 2, lens, bs, bht->nv - st->nev, truncate_lifting, leadmons[0], st,
                      GB_PRIME_IS_WIDE(fc));
    }
    else{
      gb_modpoly_init(modgbs, 2, lens, bs, bht->nv - st->nev, num_gb[0], leadmons[0], st,
                      GB_PRIME_IS_WIDE(fc));
    }
    free(lens);
    modpgbs_set(modgbs, bs, bht, fc, start, st->nev);
    int is_empty = 0;
    if(bs->lml == 1){
        if(info_level){
            fprintf(VERBSTREAM, "Grobner basis has a single element\n");
        }
        is_empty = 1;
        for(int i = 0; i < bht->nv; i++){
            if(bexp_lm[i]!=0){
                is_empty = 0;
            }
        }
        if(is_empty){
            if(bexp_lm2 != NULL){
              free(bexp_lm);
            }
            free_basis_without_hash_table(&(bs));
            return is_empty;
        }
    }

    /**************************************************/

    if(bexp_lm2 != NULL){
      /* leadmons[0] is bexp_lm2, bexp_lm is no longer needed */
      free(bexp_lm);
    }
    free_basis_without_hash_table(&(bs));
    return is_empty;
}


/* Applies the trace modulo prime (F4 with the tracer learned on bs_qq)
 * and checks that the basis has the shape of the learned one (lengths,
 * leading monomials). Returns the basis modulo prime, or NULL with *bad set
 * to 1 or 2 (see core_groebner_qq) when the prime is bad. Shared data are
 * only read (st->nthrds must be 1), so that several primes can be handled
 * at the same time by different threads, each with its own leadmons_cur. */
static bs_t *gb_modular_trace_apply(gb_modpoly_t modgbs,
                                    const int32_t num_gb,
                                    const int32_t *leadmons_ori,
                                    int32_t *leadmons_cur,
                                    bs_t *bs_qq,
                                    md_t *st,
                                    const uint64_t prime,
                                    int *bad)
{
  *bad = 0;
  int32_t error = 0;
  bs_t *bs = core_gba(bs_qq, st, &error, prime);
  /* core_gba returns NULL (and has already free-ed its data) when the
   * prime has been detected as bad while applying the tracer */
  if (bs == NULL || error) {
      *bad = 1;
      if (bs != NULL) {
          free_basis_and_only_local_hash_table_data(&bs);
      }
      return NULL;
  }
  ht_t *bht = bs->ht;
  for(len_t i = 0; i < modgbs->ld; i++){
        if(modgbs->modpolys[i]->len < bs->hm[bs->lmps[i]][LENGTH]-1){
            *bad = 2;
            free_basis_and_only_local_hash_table_data(&bs);
            return NULL;
        }
        if(modgbs->modpolys[i]->len > bs->hm[bs->lmps[i]][LENGTH]-1){
            *bad = 1;
            free_basis_and_only_local_hash_table_data(&bs);
            return NULL;
        }
  }
  int32_t lml = bs->lml;
  if (st->nev > 0) {
      int32_t j = 0;
      for (len_t i = 0; i < bs->lml; ++i) {
          if (bht->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
              bs->lm[j]   = bs->lm[i];
              bs->lmps[j] = bs->lmps[i];
              ++j;
          }
      }
      lml = j;
  }
  if (lml != num_gb) {
      *bad = 2;
      free_basis_and_only_local_hash_table_data(&bs);
      return NULL;
  }

  if(st->nev){
    get_lm_from_bs_trace_elim(bs, bht, leadmons_cur, num_gb);
  }
  else{
    get_lm_from_bs_trace(bs, bht, leadmons_cur);
  }

  if(!equal_staircase(leadmons_cur, (int32_t *)leadmons_ori,
                      num_gb, num_gb, bht->nv - st->nev)){
    *bad = 1;
    free_basis_and_only_local_hash_table_data(&bs);
    return NULL;
  }
  return bs;
}

/* LLL witnesses of the i-th polynomial: the three witnesses used for
   rational reconstruction (so that the candidates found by LLL are checked
   as the reconstructed ones), then nonzero coefficients spread over the
   polynomial, all distinct and nonzero modulo the first prime */
static inline void choose_lll_witnesses(gb_modpoly_t modgbs, data_lift_t dl,
                                        int32_t i){
  int32_t *w = dl->lwit + (size_t)i * LLL_WIT_K;
  const int32_t len = modgbs->modpolys[i]->len;
  int n = 0;
  int32_t nnz = 0;
  for(int32_t l = 0; l < len; l++){
    nnz += (MODCF(modgbs->modpolys[i], l, 0) != 0);
  }
  /* the j-th candidate is the first coefficient of the j-th group of
     nonzero coefficients, the three first ones are the RR witnesses */
  for(int32_t j = 0; j < LLL_WIT_K + 3 && n < LLL_WIT_K && nnz > 0; j++){
    int32_t c;
    if(j < 3){
      c = (j == 0) ? dl->coef[i] : ((j == 1) ? dl->coef2[i] : dl->coef3[i]);
      if(!MODCF(modgbs->modpolys[i], c, 0)){
        continue;
      }
    }
    else{
      const int32_t target = ((int64_t)nnz * (j - 3)) / LLL_WIT_K;
      int32_t cnt = -1;
      c = -1;
      for(int32_t l = 0; l < len; l++){
        cnt += (MODCF(modgbs->modpolys[i], l, 0) != 0);
        if(cnt == target){
          c = l;
          break;
        }
      }
    }
    int dup = (c < 0);
    for(int m = 0; m < n && !dup; m++){
      dup = (w[m] == c);
    }
    if(!dup){
      w[n++] = c;
    }
  }
  dl->nlwit[i] = n;
}

/* three witness coefficients are chosen per polynomial: the last nonzero
   one and the nonzero ones found scanning downwards from 2*len/3 and
   len/2 */
static inline void choose_coef_to_lift(gb_modpoly_t modgbs, data_lift_t dlift){
  uint32_t ld = modgbs->ld;
  for(int32_t i = 0; i < ld; i++){
    len_t len = modgbs->modpolys[i]->len;
    len_t d = len - 1;
    while(d > 0 && len > 0){
      if(MODCF(modgbs->modpolys[i], d, 0)){
        dlift->coef[i] = d;
        break;
      }
      else{
        d--;
      }
    }
    len_t d2 = 2 * len / 3;
    while(d2 > 0 && !MODCF(modgbs->modpolys[i], d2, 0)){
      d2--;
    }
    dlift->coef2[i] = d2;
    len_t d3 = len / 2;
    while(d3 > 0 && !MODCF(modgbs->modpolys[i], d3, 0)){
      d3--;
    }
    dlift->coef3[i] = d3;
    choose_lll_witnesses(modgbs, dlift, i);
  }
}


/* Incremental CRT (called once FLINT multi_CRT has been called) */
/* mod is the current modulus */
static inline void incremental_dlift_crt(gb_modpoly_t modgbs, data_lift_t dlift,
                                         int32_t *coef, mpz_t *mod_p, mpz_t *prod_p,
                                         int thrds){

  /* all primes are assumed to be good primes */
  mpz_mul_u64(prod_p[0], mod_p[0], modgbs->primes[modgbs->nprimes - 1 ]);
  for(int32_t k = dlift->lstart; k <= dlift->lend; k++){
      if(modgbs->modpolys[k]->len){
        uint64_t c = MODCF(modgbs->modpolys[k], coef[k], modgbs->nprimes  - 1);
        mpz_CRT_ui(dlift->crt[k], dlift->crt[k], mod_p[0],
                   c, modgbs->primes[modgbs->nprimes - 1 ],
                   prod_p[0], dlift->tmp, 1);
      }
  }
  mpz_set(mod_p[0], prod_p[0]);
}


/* Incremental CRT on the whole array of witness coefficients */
/* mod is the current modulus */
static inline void incremental_dlift_crt_full(gb_modpoly_t modgbs, data_lift_t dl,
                                              int32_t *coef, mpz_t mod_p, mpz_t prod_p,
                                              int thrds){

  mp_limb_t newprime = modgbs->primes[modgbs->nprimes - 1 ];
  /* all primes are assumed to be good primes */
  mpz_mul_u64(prod_p, mod_p, newprime);
  /* inverse of mod_p modulo the new prime, the same for all the CRTs */
  mp_limb_t c = n_invmod(mpz_fdiv_u64(mod_p, newprime), newprime);
  const mp_limb_t ninv = n_preinvert_limb(newprime);
  const uint32_t np = modgbs->nprimes - 1;
  /* only the polynomials whose witnesses have already been reconstructed
     once are kept up to date, the next ones are computed when they are
     reached (see witness_crt_window); these CRTs are cheap compared to
     moving their results between the caches of several threads (the
     witnesses are then reconstructed by this thread), hence they are not
     split among the threads */
  dl->V = MAX(dl->V, dl->S);
  for(int32_t k = dl->S; k < dl->V; k++){
      if(modgbs->modpolys[k]->len){
        uint64_t r = MODCF(modgbs->modpolys[k], coef[k], np);
        _mpz_CRT_ui_precomp(dl->crt[k], dl->crt[k], mod_p,
                            r, newprime, ninv, prod_p, c, 1);
        r = MODCF(modgbs->modpolys[k], dl->coef2[k], np);
        _mpz_CRT_ui_precomp(dl->crt2[k], dl->crt2[k], mod_p,
                            r, newprime, ninv, prod_p, c, 1);
        r = MODCF(modgbs->modpolys[k], dl->coef3[k], np);
        _mpz_CRT_ui_precomp(dl->crt3[k], dl->crt3[k], mod_p,
                            r, newprime, ninv, prod_p, c, 1);
      }
  }
  mpz_set(mod_p, prod_p);
}

/* The witness coefficients of the i-th polynomial are reconstructed while
 * its CRTs are not up to date (i >= V): those of the polynomials of index
 * V <= k < min(end, i + nthreads) are computed from their images modulo all
 * the primes, in parallel, and then kept up to date prime by prime (see
 * incremental_dlift_crt_full). The CRTs of the polynomials which are not
 * reached yet are thus not updated with each new prime, the reconstruction
 * of the witnesses stopping at the first failure. The results are the
 * symmetric remainders modulo the product of the primes, as the
 * incremental CRTs. */
static void witness_crt_window(gb_modpoly_t modgbs, data_lift_t dl,
                               const int32_t i, const int32_t end){
  const int32_t v0 = MAX(dl->V, dl->S);
  if(i < v0){
    return;
  }
  const int32_t v1 = MIN(end, i + MAX(1, dl->nthreads));
  const uint32_t np = modgbs->nprimes;
  fmpz_comb_t comb;
  fmpz_comb_init(comb, modgbs->primes, np);
  /* the team keeps the size of the other parallel regions of the lifting:
     changing it makes the OpenMP runtime rebuild it */
#pragma omp parallel num_threads(dl->nthreads)
  {
    fmpz_comb_temp_t comb_temp;
    fmpz_comb_temp_init(comb_temp, comb);
    fmpz_t y;
    fmpz_init(y);
    mp_limb_t *res = (mp_limb_t *)malloc(sizeof(mp_limb_t) * np);
#pragma omp for schedule(dynamic, 1)
    for(int32_t k = v0; k < v1; k++){
      modpolys_t *pol = modgbs->modpolys + k;
      if(pol[0]->len == 0){
        continue;
      }
      const int32_t c[3] = {dl->coef[k], dl->coef2[k], dl->coef3[k]};
      mpz_t *z[3] = {dl->crt + k, dl->crt2 + k, dl->crt3 + k};
      for(int j = 0; j < 3; j++){
        for(uint32_t l = 0; l < np; l++){
          res[l] = MODCF(pol[0], c[j], l);
        }
        fmpz_multi_CRT_ui(y, res, comb, comb_temp, 1);
        fmpz_get_mpz(*z[j], y);
      }
    }
    free(res);
    fmpz_clear(y);
    fmpz_comb_temp_clear(comb_temp);
  }
  fmpz_comb_clear(comb);
  dl->V = v1;
}


/* balanced bounds, used when no denominator is guessed; RR_MARGIN bits of
   slack are kept (N*D <= mod_p / 2^(1+RR_MARGIN)) so that spurious
   reconstructions become unlikely */
static inline void set_recdata(rrec_data_t rd, mpz_t mod_p){
  mpz_fdiv_q_2exp(rd->N, mod_p, 1 + RR_MARGIN);
  mpz_sqrt(rd->N, rd->N);
  mpz_set(rd->D, rd->N);
}

/* bounds used to reconstruct u * gden, where gden is a guessed denominator:
   the numerator is then expected to be about gden times the remaining
   denominator, so that with M = mod_p / 2^(1+RR_MARGIN), the bounds are
   D = sqrt(M / gden) and N = M / D */
static inline void set_recdata_guess(rrec_data_t rd, const mpz_t mod_p,
                                     const mpz_t gden){
  mpz_fdiv_q_2exp(rd->N, mod_p, 1 + RR_MARGIN);
  mpz_fdiv_q(rd->D, rd->N, gden);
  mpz_sqrt(rd->D, rd->D);
  if(mpz_sgn(rd->D) == 0){
    mpz_set_ui(rd->D, 1);
  }
  mpz_fdiv_q(rd->N, rd->N, rd->D);
}

/*
  For the polynomials of index in [start, end) whose witness coefficient has
  been lifted and checked (check1 > 1), lifts all coefficients by CRT
  followed by rational reconstruction, and stores the integer numerators in
  cf_qq and the common denominator in lm.

  The CRT images are only temporaries: each coefficient is lifted with
  multi-CRT into a single scratch integer which is immediately rationally
  reconstructed, so that no array of CRT images is kept alive in memory.
  Denominators are also kept in a scratch array, the size of the longest
  polynomial handled here.

  returns the index of the poly if some of its coeff could not be lifted
  else returns end
 */
static int lll_common_den(gb_modpoly_t modgbs, data_lift_t dl, int32_t i,
                          const int32_t *w, const int k, const uint32_t np,
                          const mpz_t mod_p, const mpz_t guess,
                          mpz_t den, mpz_t *num);
static inline double lll_guess_cost(const int k, const double bits);
static inline int32_t lll_next_attempt(const int32_t nprimes,
                                       const double cost, const double pcost,
                                       const int nthreads);

/*
  Coefficients of the k-th polynomial whose rational reconstruction failed
  with the guessed denominator dl->tmp (at most LLL_WIT_K of them, of
  indices fail[0], .., fail[nfail-1], with CRT images zf[0], ..): they are
  reconstructed simultaneously by LLL (see lll_common_den) with this guess,
  together with already reconstructed coefficients (indices in done[0],
  .., done[ndone-1]) when there are less than LLL_WIT_K failing ones. On
  success the guess is multiplied by the common factor found and they are
  reconstructed again with it, else 0 is returned.

  The guess is wrong if LLL returned a spurious vector: the lifted
  polynomial is then rejected by verif_lifted_basis.
 */
static int lll_full_block(gb_modpoly_t modgbs, data_lift_t dl, int32_t k,
                          int32_t *fail, int nfail, mpz_t *zf,
                          const int32_t *done, int32_t ndone,
                          mpz_t mod_p, rrec_data_t rd1, mpz_t rden,
                          mpz_t *dens, const double pcost){
  modpolys_t *polys = modgbs->modpolys;
  const uint32_t np = modgbs->nprimes - 1;
  const double bits = mpz_sizeinbase(mod_p, 2);
  int32_t w[LLL_WIT_K];
  int kk = 0;
  for(int j = 0; j < nfail; j++){
    w[kk++] = fail[j];
  }
  for(int32_t j = ndone - 1; j >= 0 && kk < LLL_WIT_K; j--){
    w[kk++] = done[j];
  }
  const double cost = lll_guess_cost(kk, bits);
  if(kk < 2 || (int32_t)np < dl->lll_fnext
     || (pcost > 0 && dl->lll_cost + cost
         > LLL_WIT_BUDGET / 100.0 * modgbs->nprimes * pcost)){
    return 0;
  }
  dl->lll_cost += cost;
  mpz_t d, num[LLL_WIT_K];
  mpz_init(d);
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_init(num[j]);
  }
  int ok = lll_common_den(modgbs, dl, k, w, kk, np, mod_p, dl->tmp, d, num);
  dl->lll_fnext = ok ? (int32_t)np
    : lll_next_attempt(np, cost, pcost, dl->nthreads);
  if(ok){
    mpz_mul(dl->tmp, dl->tmp, d);
    set_recdata_guess(rd1, mod_p, dl->tmp);
    for(int j = 0; j < nfail && ok; j++){
      const int32_t l = fail[j];
      ok = lift_ratreconwden(polys[k]->cf_qq[l], rden, zf[j], mod_p,
                             dl->tmp, rd1);
      if(ok){
        mpz_mul(dens[l], rden, dl->tmp);
        if(mpz_cmp_ui(rden, 1)){
          mpz_mul(dl->tmp, dl->tmp, rden);
          set_recdata_guess(rd1, mod_p, dl->tmp);
        }
      }
    }
  }
  mpz_clear(d);
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_clear(num[j]);
  }
  return ok;
}

/* number of coefficients per thread in a chunk of crt_coefs_chunk */
#define CRT_CHUNK 32
/* number of consecutive coefficients whose residues are read together */
#define CRT_BLOCK 4

/* For the coefficients l0 <= l < l1 of the polynomial pol, CRT of the
 * residues of u g, where gp holds the guess g modulo the np first primes,
 * as a symmetric remainder modulo their product (of mbits bits): when it
 * has RR_MARGIN bits less than the modulus, it is the integer u g, stored
 * in cf_qq[l] (isint[l - l0] = 1), else the CRT of u itself is stored in
 * zraw[l - l0] (isint[l - l0] = 0) for a rational reconstruction (see
 * crt_ratrecon_lift_modgbs). The coefficients are independent: they are
 * split among nthreads threads, each with its own CRT workspace. */
/* workspace of a thread for crt_coefs_chunk, kept for all the chunks of a
 * call to crt_ratrecon_lift_modgbs (allocating the CRT temporaries for each
 * chunk would serialize the threads in the memory allocator) */
typedef struct{
  fmpz_comb_temp_t comb_temp;
  fmpz_t y;
  mp_limb_t *res; /* residues of CRT_BLOCK coefficients */
  mp_limb_t *raw;
} crt_ws_struct;

static crt_ws_struct *crt_ws_init(const fmpz_comb_t comb, const uint32_t np,
                                  const int nthreads){
  crt_ws_struct *ws = (crt_ws_struct *)malloc(sizeof(crt_ws_struct) * nthreads);
  for(int t = 0; t < nthreads; t++){
    fmpz_comb_temp_init(ws[t].comb_temp, comb);
    fmpz_init(ws[t].y);
    ws[t].res = (mp_limb_t *)malloc(sizeof(mp_limb_t) * np * CRT_BLOCK);
    ws[t].raw = (mp_limb_t *)malloc(sizeof(mp_limb_t) * np);
  }
  return ws;
}

static void crt_ws_clear(crt_ws_struct *ws, const int nthreads){
  for(int t = 0; t < nthreads; t++){
    fmpz_comb_temp_clear(ws[t].comb_temp);
    fmpz_clear(ws[t].y);
    free(ws[t].res);
    free(ws[t].raw);
  }
  free(ws);
}

static void crt_coefs_chunk(gb_modpoly_t modgbs, modpolys_t *pol,
                            const int32_t l0, const int32_t l1,
                            const mp_limb_t *gp, const uint32_t np,
                            const fmpz_comb_t comb, const size_t mbits,
                            int *isint, mpz_t *zraw, crt_ws_struct *ws,
                            const int nthreads){
  /* the images are stored prime by prime: the residues of CRT_BLOCK
     consecutive coefficients are read together */
  const int32_t nblocks = (l1 - l0 + CRT_BLOCK - 1) / CRT_BLOCK;
  const size_t len = pol[0]->len;
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
  for(int32_t blk = 0; blk < nblocks; blk++){
    crt_ws_struct *w = ws + omp_get_thread_num();
    mp_limb_t *res = w->res;
    const int32_t b0 = l0 + blk * CRT_BLOCK;
    const int32_t nb = MIN(CRT_BLOCK, l1 - b0);
    if(modgbs->wide){
      for(uint32_t i = 0; i < np; i++){
        const uint64_t *cf = pol[0]->cf_u64 + (size_t)i * len + b0;
        for(int32_t j = 0; j < nb; j++){
          res[(size_t)j * np + i] = n_mulmod2_preinv(cf[j], gp[i],
                                                     modgbs->primes[i],
                                                     modgbs->pinvs[i]);
        }
      }
    }
    else{
      for(uint32_t i = 0; i < np; i++){
        const uint32_t *cf = pol[0]->cf_32 + (size_t)i * len + b0;
        for(int32_t j = 0; j < nb; j++){
          res[(size_t)j * np + i] = ((uint64_t)cf[j] * gp[i]) % modgbs->primes[i];
        }
      }
    }
    for(int32_t j = 0; j < nb; j++){
      const int32_t l = b0 + j;
      fmpz_multi_CRT_ui(w->y, res + (size_t)j * np, comb, w->comb_temp, 1);
      isint[l - l0] = (fmpz_bits(w->y) + 1 + RR_MARGIN <= mbits);
      if(isint[l - l0]){
        fmpz_get_mpz(pol[0]->cf_qq[l], w->y);
      }
      else{
        for(uint32_t i = 0; i < np; i++){
          w->raw[i] = MODCF(pol[0], l, i);
        }
        fmpz_multi_CRT_ui(w->y, w->raw, comb, w->comb_temp, 1);
        fmpz_get_mpz(zraw[l - l0], w->y);
      }
    }
  }
}

/* For the coefficients l0 <= l < l1 of the chunk starting at cl0 (see
 * crt_coefs_chunk), computed with the guess g while the guess is now
 * g' = g * ratio: symmetric remainder modulo mod_p (of mbits bits, half =
 * mod_p / 2) of u g', obtained from u g when it is an integer (isint) or
 * else from the CRT of u (zraw), stored in fixv; fixok tells whether it is
 * an integer with RR_MARGIN bits less than mod_p, as the CRT of the
 * residues of u g' would. When it is not and u g was an integer, the CRT
 * of u is computed in zraw for the rational reconstruction. The
 * coefficients are independent and split among nthreads threads. */
static void crt_update_chunk(gb_modpoly_t modgbs, modpolys_t *pol,
                             const int32_t l0, const int32_t l1,
                             const int32_t cl0, const int *isint, mpz_t *zraw,
                             mpz_t *fixv, int *fixok, const mpz_t ratio,
                             const mpz_t guess, const mpz_t mod_p,
                             const mpz_t half, const size_t mbits,
                             const uint32_t np, const fmpz_comb_t comb,
                             crt_ws_struct *ws, const int nthreads){
#pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
  for(int32_t l = l0; l < l1; l++){
    const int32_t j = l - cl0;
    if(isint[j]){
      mpz_mul(fixv[j], pol[0]->cf_qq[l], ratio);
    }
    else{
      mpz_mul(fixv[j], zraw[j], guess);
    }
    mpz_mod(fixv[j], fixv[j], mod_p);
    if(mpz_cmp(fixv[j], half) > 0){
      mpz_sub(fixv[j], fixv[j], mod_p);
    }
    fixok[j] = (mpz_sizeinbase(fixv[j], 2) + 1 + RR_MARGIN <= mbits);
    if(!fixok[j] && isint[j]){
      /* CRT of u itself for the rational reconstruction */
      crt_ws_struct *w = ws + omp_get_thread_num();
      for(uint32_t i = 0; i < np; i++){
        w->raw[i] = MODCF(pol[0], l, i);
      }
      fmpz_multi_CRT_ui(w->y, w->raw, comb, w->comb_temp, 1);
      fmpz_get_mpz(zraw[j], w->y);
    }
  }
}

static inline int crt_ratrecon_lift_modgbs(gb_modpoly_t modgbs, data_lift_t dl,
                                           int32_t start, int32_t end,
                                           mpz_t mod_p, rrec_data_t rd1, rrec_data_t rd2,
                                           double *st_crt, double *st_rrec,
                                           const double pcost){

  modpolys_t *polys = modgbs->modpolys;

  /* scratch data */
  uint32_t maxlen = 0;
  for(int32_t k = start; k < end; k++){
    if(dl->check1[k] > 1 && polys[k]->len > maxlen){
      maxlen = polys[k]->len;
    }
  }
  mpz_t rden, zz;
  mpz_init(rden);
  mpz_init(zz);
  mpz_t *dens = NULL;
  int32_t *done = NULL;
  if(maxlen > 0){
    dens = (mpz_t *)malloc(sizeof(mpz_t) * maxlen);
    for(uint32_t l = 0; l < maxlen; l++){
      mpz_init(dens[l]);
    }
    done = (int32_t *)malloc(sizeof(int32_t) * maxlen);
  }
  /* coefficients whose reconstruction failed, see lll_full_block */
  int32_t fail[LLL_WIT_K];
  mpz_t zf[LLL_WIT_K];
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_init(zf[j]);
  }

  /* Data needed by multi CRT functions, mod_p is the product of the
     np first primes (the last one is used for checking) */
  const uint32_t np = modgbs->nprimes - 1;
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;
  fmpz_t y;
  /* the guess dl->tmp modulo the np first primes, see below */
  mp_limb_t *gp = NULL;
  int gvalid = 0;
  const size_t mbits = mpz_sizeinbase(mod_p, 2);
  /* results of crt_coefs_chunk for the coefficients cl0 <= l < cl1 */
  const int32_t chunk = CRT_CHUNK * dl->nthreads;
  int *isint = NULL;
  mpz_t *zraw = NULL;
  crt_ws_struct *ws = NULL;
  /* updates of the chunk for the guess of version fixver (see
     crt_update_chunk); gver changes with the guess */
  mpz_t *fixv = NULL;
  int *fixok = NULL;
  long gver = 0, fixver = -1;
  int32_t cl0 = 0, cl1 = 0;
  /* guess with which the chunk was computed, and current guess / it */
  mpz_t gchunk, ratio, half;
  mpz_init(gchunk);
  mpz_init(ratio);
  mpz_init(half);
  mpz_fdiv_q_2exp(half, mod_p, 1);
  if(maxlen > 0){
    fmpz_comb_init(comb, modgbs->primes, np);
    fmpz_comb_temp_init(comb_temp, comb);
    fmpz_init(y);
    gp = (mp_limb_t *)malloc(sizeof(mp_limb_t) * np);
    ws = crt_ws_init(comb, np, dl->nthreads);
    isint = (int *)malloc(sizeof(int) * chunk);
    fixok = (int *)malloc(sizeof(int) * chunk);
    zraw = (mpz_t *)malloc(sizeof(mpz_t) * chunk);
    fixv = (mpz_t *)malloc(sizeof(mpz_t) * chunk);
    for(int32_t j = 0; j < chunk; j++){
      mpz_init(zraw[j]);
      mpz_init(fixv[j]);
    }
  }

  set_recdata(rd2, mod_p);

  int32_t ret = end;
  double rt;

  for(int32_t k = start; k < end; k++){
    if(dl->check1[k] > 1 && polys[k]->len > 0){
      /* dl->tmp is the lcm of the denominators found so far, starting
         from the witness ones; it is used as a guess for the next
         coefficient, so that only the missing factor is reconstructed */
      mpz_set(dl->tmp, dl->cden[k]);
      set_recdata_guess(rd1, mod_p, dl->tmp);
      gvalid = 0; gver++;
      cl0 = cl1 = 0;
      /* the guess mostly changes at the first coefficients of a
         polynomial: the chunks start with one coefficient per thread and
         double up to their maximal size, so that few coefficients are
         computed with a guess about to change */
      int32_t csize = dl->nthreads;
      int ok = 1;
      int nfail = 0;
      int32_t ndone = 0;
      int used_lll = 0;
      for(int32_t l = 0; l < polys[k]->len; l++){
        rt = realtime();
        /* When the denominator of the coefficient u divides the guess g,
           u g is an integer, given by the CRT of the residues of u g
           modulo the primes (symmetric remainder modulo mod_p): it is
           accepted when it has RR_MARGIN bits less than mod_p, which
           avoids both the product u g and a rational reconstruction for
           most coefficients. Else u is lifted by CRT and rationally
           reconstructed with the guess g. A wrong integer (probability
           about 2^(-RR_MARGIN)) is rejected when checking the lifted
           polynomial with the next primes. */
        if(l >= cl1){
          /* the CRT of the next chunk of coefficients with the current
             guess g, in parallel */
          if(!gvalid){
#pragma omp parallel for num_threads(dl->nthreads) schedule(static)
            for(uint32_t i = 0; i < np; i++){
              gp[i] = mpz_fdiv_u64(dl->tmp, modgbs->primes[i]);
            }
            gvalid = 1;
          }
          mpz_set(gchunk, dl->tmp);
          fixver = -1;
          cl0 = l;
          cl1 = MIN((int32_t)polys[k]->len, l + csize);
          csize = MIN(2 * csize, chunk);
          crt_coefs_chunk(modgbs, polys + k, cl0, cl1, gp, np, comb, mbits,
                          isint, zraw, ws, dl->nthreads);
        }
        int b = isint[l - cl0];
        if(mpz_cmp(dl->tmp, gchunk)){
          /* the guess has become g' = g r since the chunk was computed:
             the symmetric remainder of u g' modulo mod_p, which the CRT of
             the residues of u g' would give, is obtained from u g (integer)
             or from the CRT of u, so that the decisions are those of a
             coefficient by coefficient loop */
          if(fixver != gver){
            /* all the remaining coefficients of the chunk are updated for
               the current guess, in parallel */
            mpz_divexact(ratio, dl->tmp, gchunk);
            crt_update_chunk(modgbs, polys + k, l, cl1, cl0, isint, zraw,
                             fixv, fixok, ratio, dl->tmp, mod_p, half, mbits,
                             np, comb, ws, dl->nthreads);
            fixver = gver;
          }
          b = fixok[l - cl0];
          if(b){
            mpz_swap(polys[k]->cf_qq[l], fixv[l - cl0]);
          }
          else{
            mpz_set(zz, zraw[l - cl0]);
          }
        }
        else if(!b){
          mpz_set(zz, zraw[l - cl0]);
        }
        if(b){
          mpz_set_ui(rden, 1);
        }
        *st_crt += realtime() - rt;
        rt = realtime();
        if(!b){
          b = lift_ratreconwden(polys[k]->cf_qq[l], rden, zz, mod_p,
                                dl->tmp, rd1);
        }
        if(b){
          mpz_mul(dens[l], rden, dl->tmp);
          if(mpz_cmp_ui(rden, 1)){
            mpz_mul(dl->tmp, dl->tmp, rden);
            set_recdata_guess(rd1, mod_p, dl->tmp);
            gvalid = 0; gver++;
          }
          if(mpz_sgn(polys[k]->cf_qq[l])){
            done[ndone++] = l;
          }
        }
        else{
          /* kept for a simultaneous reconstruction by LLL */
          mpz_set(zf[nfail], zz);
          fail[nfail++] = l;
          b = 1;
        }
        if(nfail == LLL_WIT_K
           || (nfail > 0 && l == (int32_t)polys[k]->len - 1)){
          b = lll_full_block(modgbs, dl, k, fail, nfail, zf, done, ndone,
                             mod_p, rd1, rden, dens, pcost);
          used_lll = 1;
          gvalid = 0; gver++;
          for(int j = 0; j < nfail && b; j++){
            if(mpz_sgn(polys[k]->cf_qq[fail[j]])){
              done[ndone++] = fail[j];
            }
          }
          nfail = 0;
        }
        *st_rrec += realtime() - rt;
        if(!b){
          /* the witnesses of this poly are still valid: their common
             denominator remains the best guess for reconstructing them */
          mpz_set(dl->gden, dl->cden[k]);
          for(int32_t i = k; i < modgbs->ld; i++){
            dl->check2[i] = 0;
          }
          ok = 0;
          break;
        }
      }
      if(!ok){
        ret = k;
        break;
      }
      mpz_set(polys[k]->lm, dl->tmp);
      for(int32_t l = 0; l < polys[k]->len; l++){
        /* num * (lcm / den) rather than (num * lcm) / den: same value, but
         * the limb array of the (long-lived) result is not over-allocated
         * by the size of lcm (mpz_divexact does not shrink its output) */
        mpz_divexact(dens[l], dl->tmp, dens[l]);
        mpz_mul(polys[k]->cf_qq[l], polys[k]->cf_qq[l], dens[l]);
      }
      if(k > 0 && polys[k - 1]->len > 0 && mpz_cmp_ui(polys[k - 1]->lm, 1)){
        /* share of the denominator of the previous polynomial, which is
           the guess for this one, not dividing its denominator */
        mpz_gcd(rden, polys[k - 1]->lm, dl->tmp);
        const double b = mpz_sizeinbase(polys[k - 1]->lm, 2);
        const double f = (b - mpz_sizeinbase(rden, 2)) / b;
        dl->infl = MAX(dl->infl, f);
      }
      if(used_lll){
        /* the factor found by LLL may be a multiple of the missing part of
           the denominator: the content is removed */
        mpz_set(rden, polys[k]->lm);
        for(int32_t l = 0; l < polys[k]->len && mpz_cmp_ui(rden, 1); l++){
          mpz_gcd(rden, rden, polys[k]->cf_qq[l]);
        }
        if(mpz_cmp_ui(rden, 1)){
          mpz_divexact(polys[k]->lm, polys[k]->lm, rden);
          for(int32_t l = 0; l < polys[k]->len; l++){
            mpz_divexact(polys[k]->cf_qq[l], polys[k]->cf_qq[l], rden);
          }
        }
      }
      /* lifted, not checked yet: a polynomial lifted again after the
         failed check of a previous one may have been lifted before, but
         its new lift is unchecked (see verif_lifted_basis) */
      dl->check2[k] = 1;
    }
    if(dl->check1[k] < 1 && polys[k]->len > 0){
      ret = k;
      break;
    }
    if(polys[k]->len == 0){
        dl->check2[k] += 1;
    }
  }

  if(maxlen > 0){
    for(int32_t j = 0; j < chunk; j++){
      mpz_clear(zraw[j]);
      mpz_clear(fixv[j]);
    }
    free(zraw);
    free(fixv);
    free(isint);
    free(fixok);
    crt_ws_clear(ws, dl->nthreads);
    free(gp);
    fmpz_clear(y);
    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);
    for(uint32_t l = 0; l < maxlen; l++){
      mpz_clear(dens[l]);
    }
    free(dens);
    free(done);
  }
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_clear(zf[j]);
  }
  mpz_clear(gchunk);
  mpz_clear(ratio);
  mpz_clear(half);
  mpz_clear(zz);
  mpz_clear(rden);

  return ret;
}


/* returns (coef == num / den mod prime) */
static inline int verif_coef(mpz_t num, mpz_t den, const mp_limb_t prime,
                             const mp_limb_t pinv, const uint64_t coef){
  mp_limb_t lc = mpz_fdiv_u64(den, prime);
  lc = n_invmod(lc, prime);

  mp_limb_t c = mpz_fdiv_u64(num, prime);
  c = n_mulmod2_preinv(c, lc, prime, pinv);

  return (c==coef);
}


static inline int verif_lifted_basis(gb_modpoly_t modgbs, data_lift_t dl,
                                     int thrds){
  for(int32_t k = 0; k < modgbs->ld; k++){
    if(dl->check1[k]>=1 && dl->check2[k]>0 && dl->check2[k] < NBCHECK ){
      modpolys_t *pol = modgbs->modpolys + k;
      if(pol[0]->len == 0){
        dl->check2[k]++;
      }
      else{
        for(int i = 0; i < thrds; i++){
          const uint32_t np = modgbs->nprimes - (thrds - i);
          const mp_limb_t prime = modgbs->primes[np];
          const mp_limb_t pinv = modgbs->pinvs[np];
          /* all coefficients share the denominator lm */
          mp_limb_t ilc = mpz_fdiv_u64(pol[0]->lm, prime);
          ilc = n_invmod(ilc, prime);
          for(int32_t c = 0; c < pol[0]->len; c++){
            mp_limb_t v = mpz_fdiv_u64(pol[0]->cf_qq[c], prime);
            v = n_mulmod2_preinv(v, ilc, prime, pinv);
            if(v != MODCF(pol[0], c, np)){
              dl->check2[k] = 0;
              return 0;
            }
          }
          dl->check2[k]++;
          if(dl->check2[k] >= NBCHECK){
            dl->S = k+1;
            dl->lstart = dl->S;
            /* modular images are no longer needed */
            free(pol[0]->cf_32);
            pol[0]->cf_32 = NULL;
            free(pol[0]->cf_u64);
            pol[0]->cf_u64 = NULL;
            break;
          }
        }
      }
    }
  }
  return 1;
}

static inline int verif_lifted_rational_wcoef(gb_modpoly_t modgbs, data_lift_t dl,
                                              int thrds){
  for(int32_t k = dl->lstart; k < dl->lend; k++){
    if(dl->check1[k]==0){
      /* too early to perform the verification */
      return k;
    }
    if(modgbs->modpolys[k]->len == 0){
        dl->start = MIN(dl->start + 1, dl->lend);
        dl->check1[k]++;
    }
    else{
      for(int i = 0; i < thrds; i++){

        mp_limb_t prime = modgbs->primes[modgbs->nprimes - (thrds - i) ];
        mp_limb_t pinv = modgbs->pinvs[modgbs->nprimes - (thrds - i) ];
        uint64_t coef = MODCF(modgbs->modpolys[k], dl->coef[k], modgbs->nprimes  - (thrds - i));
        int boo = verif_coef(dl->num[k], dl->den[k], prime, pinv, coef);
        if(boo){
          coef = MODCF(modgbs->modpolys[k], dl->coef2[k], modgbs->nprimes  - (thrds - i));
          boo = verif_coef(dl->num2[k], dl->den2[k], prime, pinv, coef);
        }
        if(boo){
          coef = MODCF(modgbs->modpolys[k], dl->coef3[k], modgbs->nprimes  - (thrds - i));
          boo = verif_coef(dl->num3[k], dl->den3[k], prime, pinv, coef);
        }

        if(!boo){
          for(int32_t kk = k; kk < dl->end; kk++){
            dl->check1[kk] = 0;
          }
          dl->start = k;
          return k;
        }
        dl->start = MIN(dl->start + 1, dl->lend);
        dl->check1[k]++;
      }
    }
  }
  return -1;
}

/* returns 0 iff rational number could not be lifted */
static inline int reconstructcoeff(data_lift_t dl, int32_t i, mpz_t mod_p,
                                   rrec_data_t recdata1, rrec_data_t recdata2){

  set_recdata_guess(recdata1, mod_p, dl->gden);
  int b = lift_ratreconwden(dl->num[i], dl->den[i],
                            dl->crt[i], mod_p, dl->gden, recdata1);
  if(b){
    mpz_mul(dl->den[i], dl->den[i], dl->gden);
    mpz_gcd(dl->tmp, dl->den[i], dl->num[i]);

    mpz_divexact(dl->num[i], dl->num[i], dl->tmp);
    mpz_divexact(dl->den[i], dl->den[i], dl->tmp);

    mpz_set(dl->gden, dl->den[i]);
  }
  else{

    b = lift_ratrecon(dl->num[i], dl->den[i],
                      dl->crt[i], mod_p, recdata2);
    if(b){
      mpz_set(dl->gden, dl->den[i]);
    }
  }
  return b;
}

/* reconstructs the second witness using as a guess the lcm of the
   denominator of the first one and of the one guessed from the previous
   polynomial; cden[i] is set to the lcm of both witness denominators */
static inline int reconstructcoeff2(data_lift_t dl, int32_t i, mpz_t mod_p,
                                    rrec_data_t recdata1, rrec_data_t recdata2){
  mpz_lcm(dl->cden[i], dl->den[i], dl->pgden);
  set_recdata_guess(recdata1, mod_p, dl->cden[i]);
  if(!lift_ratreconwden(dl->num2[i], dl->den2[i],
                        dl->crt2[i], mod_p, dl->cden[i], recdata1)){
    /* same fallback as for the first witness: balanced bounds, no guess */
    if(!lift_ratrecon(dl->num2[i], dl->den2[i], dl->crt2[i], mod_p, recdata2)){
      return 0;
    }
    mpz_lcm(dl->cden[i], dl->den[i], dl->den2[i]);
    mpz_set(dl->gden, dl->cden[i]);
    return 1;
  }
  mpz_mul(dl->den2[i], dl->den2[i], dl->cden[i]);
  mpz_gcd(dl->tmp, dl->den2[i], dl->num2[i]);
  mpz_divexact(dl->num2[i], dl->num2[i], dl->tmp);
  mpz_divexact(dl->den2[i], dl->den2[i], dl->tmp);
  mpz_lcm(dl->cden[i], dl->den[i], dl->den2[i]);
  mpz_set(dl->gden, dl->cden[i]);
  return 1;
}

/* reconstructs the third witness using the lcm of the denominators of the
   first two witnesses as a guess; cden[i] is updated accordingly */
static inline int reconstructcoeff3(data_lift_t dl, int32_t i, mpz_t mod_p,
                                    rrec_data_t recdata1, rrec_data_t recdata2){
  set_recdata_guess(recdata1, mod_p, dl->cden[i]);
  if(!lift_ratreconwden(dl->num3[i], dl->den3[i],
                        dl->crt3[i], mod_p, dl->cden[i], recdata1)){
    /* same fallback as for the first witness: balanced bounds, no guess */
    if(!lift_ratrecon(dl->num3[i], dl->den3[i], dl->crt3[i], mod_p, recdata2)){
      return 0;
    }
    mpz_lcm(dl->cden[i], dl->cden[i], dl->den3[i]);
    mpz_set(dl->gden, dl->cden[i]);
    return 1;
  }
  mpz_mul(dl->den3[i], dl->den3[i], dl->cden[i]);
  mpz_gcd(dl->tmp, dl->den3[i], dl->num3[i]);
  mpz_divexact(dl->num3[i], dl->num3[i], dl->tmp);
  mpz_divexact(dl->den3[i], dl->den3[i], dl->tmp);
  mpz_lcm(dl->cden[i], dl->cden[i], dl->den3[i]);
  mpz_set(dl->gden, dl->cden[i]);
  return 1;
}

/* sets num / den to the reduced fraction n / d */
static inline void set_reduced_fraction(mpz_t num, mpz_t den,
                                        const mpz_t n, const mpz_t d,
                                        mpz_t tmp){
  mpz_gcd(tmp, n, d);
  mpz_divexact(num, n, tmp);
  mpz_divexact(den, d, tmp);
}

/*
  Reconstructs simultaneously the k coefficients of index w[0], .., w[k-1]
  of the i-th polynomial, which share the denominator g * d, where g is a
  guessed denominator and d is unknown, from their images modulo the np
  first primes, whose product is mod_p.

  With u_1, .., u_k the CRT images of the coefficients modulo m, the vector
  (2^s d, n_1, .., n_k) with n_j = g d u_j mod m is a short vector of the
  lattice spanned by (2^s, g u_1, .., g u_k) and m e_1, .., m e_k. Since
  n_j / (g d) is a coefficient of the polynomial, n_j has about the bitsize
  of g d times the size of the coefficient: the first column is scaled by
  2^s, s = bits(g) - 1, so that all entries have about the same size. LLL
  then finds this vector, heuristically, once bits(m) exceeds about
  bits(n) + bits(d) / k, instead of bits(n) + 2 bits(d) when a single
  coefficient is reconstructed with the guess g.

  Since n_j has about s + bits(d) bits, the t = s - LLL_TRUNC_SLACK lowest
  bits of the entries are irrelevant: LLL is run on the lattice spanned by
  (2^(s-t), x_1, .., x_k) and m' e_1, .., m' e_k, with x_j and m' the
  entries g u_j mod m and m shifted right by t bits. It contains
  (2^(s-t) d, n_1 / 2^t + e_1, ..), with |e_j| <= 2 d + 1, which is still the
  short vector, while its entries have bits(m) - t bits instead of bits(m):
  this makes LLL several times faster for a large guess.

  The n_j are then recomputed from d as g d u_j mod m (symmetric remainder)
  and a candidate is accepted when 2^s d and all the n_j are at most
  2^((s + k bits(m)) / (k + 1) - RR_MARGIN); d (positive) and the n_j are
  then stored in den and num[0], .., num[k-1], so that the j-th coefficient
  is num[j] / (g den).

  returns 1 if a candidate has been found, else 0.
 */
static int lll_common_den(gb_modpoly_t modgbs, data_lift_t dl, int32_t i,
                          const int32_t *w, const int k, const uint32_t np,
                          const mpz_t mod_p, const mpz_t guess,
                          mpz_t den, mpz_t *num){
  if(k < 2){
    return 0;
  }
  double rt = realtime();
  modpolys_t *pol = modgbs->modpolys + i;

  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;
  fmpz_comb_init(comb, modgbs->primes, np);
  fmpz_comb_temp_init(comb_temp, comb);
  fmpz_t y, m, g;
  fmpz_init(y);
  fmpz_init(m);
  fmpz_init(g);
  fmpz_set_mpz(m, mod_p);
  fmpz_set_mpz(g, guess);
  const slong s = fmpz_bits(g) > 1 ? fmpz_bits(g) - 1 : 0;
  const slong t = s > LLL_TRUNC_SLACK ? s - LLL_TRUNC_SLACK : 0;

  /* x_j = g u_j mod m */
  fmpz *x = _fmpz_vec_init(k);
  fmpz_mat_t B;
  fmpz_mat_init(B, k + 1, k + 1);
  fmpz_one(fmpz_mat_entry(B, 0, 0));
  fmpz_mul_2exp(fmpz_mat_entry(B, 0, 0), fmpz_mat_entry(B, 0, 0), s - t);
  for(int j = 0; j < k; j++){
    for(uint32_t l = 0; l < np; l++){
      modgbs->cf_64[l] = MODCF(pol[0], w[j], l);
    }
    fmpz_multi_CRT_ui(y, modgbs->cf_64, comb, comb_temp, 0);
    fmpz_mul(y, y, g);
    fmpz_mod(x + j, y, m);
    fmpz_fdiv_q_2exp(fmpz_mat_entry(B, 0, j + 1), x + j, t);
    fmpz_fdiv_q_2exp(fmpz_mat_entry(B, j + 1, j + 1), m, t);
  }
  fmpz_lll_t fl;
  fmpz_lll_context_init_default(fl);
  /* the reduction is not certified: on these small lattices, fmpz_lll
     spends most of its time beyond the floating-point reduction (up to 10
     times more with FLINT 3.1), while the candidates found are checked
     anyway, by their sizes below and then with new primes */
  if(fmpz_lll_d_heuristic(B, NULL, fl) == -1){
    fmpz_lll(B, NULL, fl);
  }

  const slong bits = (s + (slong)k * (slong)fmpz_bits(m)) / (k + 1) - RR_MARGIN;
  fmpz *d = fmpz_mat_entry(B, 0, 0);
  int ok = bits > 0 && !fmpz_is_zero(d)
    && (slong)fmpz_val2(d) >= s - t;
  if(ok){
    fmpz_fdiv_q_2exp(d, d, s - t);
    fmpz_abs(d, d);
    ok = ((slong)fmpz_bits(d) + s <= bits);
  }
  for(int j = 0; j < k && ok; j++){
    fmpz_mul(y, d, x + j);
    fmpz_smod(x + j, y, m);
    ok = ((slong)fmpz_bits(x + j) <= bits);
  }
  if(ok){
    fmpz_get_mpz(den, d);
    for(int j = 0; j < k; j++){
      fmpz_get_mpz(num[j], x + j);
    }
    dl->lll_nsucc++;
  }

  _fmpz_vec_clear(x, k);
  fmpz_mat_clear(B);
  fmpz_clear(y);
  fmpz_clear(m);
  fmpz_clear(g);
  fmpz_comb_temp_clear(comb_temp);
  fmpz_comb_clear(comb);
  dl->lll_ntries++;
  dl->lll_time += realtime() - rt;
  return ok;
}

/*
  Reconstructs simultaneously (see lll_common_den) the LLL witnesses of the
  i-th polynomial, at most kmax of them, with the guessed denominator
  guess. This matters for the polynomials whose denominator has a large
  factor which is not in the previous ones (e.g. the first polynomial, for
  which no denominator is guessed).

  The candidate is stored as the reconstructed values of the three RR
  witnesses (and their lcm in cden), so that it is checked with the next
  prime as these.

  returns 1 if a candidate has been found, else 0.
 */
static int lll_witnesses(gb_modpoly_t modgbs, data_lift_t dl, int32_t i,
                         const mpz_t mod_p, const mpz_t guess, const int kmax){
  const int k = MIN(dl->nlwit[i], kmax);
  const int32_t *w = dl->lwit + (size_t)i * LLL_WIT_K;
  mpz_t den, gd, num[LLL_WIT_K];
  mpz_init(den);
  mpz_init(gd);
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_init(num[j]);
  }
  const int ok = lll_common_den(modgbs, dl, i, w, k, modgbs->nprimes,
                                mod_p, guess, den, num);
  if(ok){
    mpz_mul(gd, den, guess);
    mpz_set_ui(dl->cden[i], 1);
    for(int j = 0; j < k; j++){
      /* stores the witnesses used by rational reconstruction */
      if(w[j] == dl->coef[i]){
        set_reduced_fraction(dl->num[i], dl->den[i], num[j], gd, dl->tmp);
      }
      if(w[j] == dl->coef2[i]){
        set_reduced_fraction(dl->num2[i], dl->den2[i], num[j], gd, dl->tmp);
      }
      if(w[j] == dl->coef3[i]){
        set_reduced_fraction(dl->num3[i], dl->den3[i], num[j], gd, dl->tmp);
      }
      /* the denominator of num[j] / (g d) */
      mpz_gcd(dl->tmp, num[j], gd);
      mpz_divexact(dl->tmp, gd, dl->tmp);
      mpz_lcm(dl->cden[i], dl->cden[i], dl->tmp);
    }
    mpz_set(dl->gden, dl->cden[i]);
  }
  mpz_clear(den);
  mpz_clear(gd);
  for(int j = 0; j < LLL_WIT_K; j++){
    mpz_clear(num[j]);
  }
  return ok;
}

/* estimated elapsed time of a call to lll_common_den with k coefficients
   and a guessed denominator, modulo an integer of the given bitsize
   (fitted with FLINT's fmpz_lll on random2_7_7, random2_8_8 and
   minrank_d_2_p_3_q_7_r_3, 3000 to 72000 bits; the cost grows about as
   (k + 1)^2.37). It was fitted before the lattices were truncated (see
   lll_common_den), which makes LLL 2 to 4 times faster for a large guess:
   it is kept as an upper bound, so that the attempts, hence the number of
   primes, do not change. */
static inline double lll_guess_cost(const int k, const double bits){
  return 7.6e-10 * pow(k + 1, 2.37) * pow(bits, 1.3);
}

/* same for an attempt of reconstructing the witnesses: LLL_WIT_K of them
   with the guessed denominator, then, if noguess is nonzero, LLL_WIT_K1 of
   them without guess, a lattice on which LLL is much slower */
static inline double lll_witnesses_cost(const double bits, const int noguess){
  return lll_guess_cost(LLL_WIT_K, bits)
    + (noguess ? 5.0e-11 * pow(LLL_WIT_K1 + 1, 2.37) * pow(bits, 1.63) : 0);
}

/* estimated elapsed time of the application of the trace for one prime,
   from the operations recorded in the trace, with the weights of the cost
   model of hensel.c (see HENSEL_KAPPA_F4); 0 if no trace is learned */
static inline double lll_prime_cost(const md_t *st){
  if(st == NULL || st->tr == NULL){
    return 0;
  }
  return 1.5e-10 * (double)st->tr->nmadd + 2.6e-10 * (double)st->tr->nmadd_new
    + 2.2e-9 * (double)st->tr->nterms;
}

/* number of primes of the next attempt of LLL after a failed one with
   nprimes primes: LLL costs about A primes (A is the ratio of the cost of
   the attempt to the cost of a prime), and trying again after a growth of
   g percent of the number of primes loses up to g * nprimes primes when
   LLL succeeds, hence the growth is about sqrt(2 A / nprimes), at least
   LLL_WIT_GROWTH percent. With nthreads threads, the trace is applied to
   nthreads primes at once while LLL runs on one thread: in elapsed time,
   an attempt then costs nthreads times more primes. It only depends on the
   trace, on the number of primes and on the number of threads, so that the
   number of primes used is reproducible for a given number of threads. */
static inline int32_t lll_next_attempt(const int32_t nprimes,
                                       const double cost, const double pcost,
                                       const int nthreads){
  double g = LLL_WIT_GROWTH / 100.0;
  if(pcost > 0){
    const double A = MAX(1, nthreads) * cost / pcost;
    g = MAX(g, sqrt(2 * A / nprimes));
  }
  return MAX(nprimes + 1, (int32_t)(nprimes * (1 + g)));
}

/* above this number of threads, the witnesses are reconstructed once per
   batch of primes (see witness_attempt_due) */
#define WIT_BATCH_MIN 4
/* number of primes of a batch left after an attempt, see
   witness_attempt_due */
#define WIT_BATCH_SLACK 2

/* whether the witnesses are reconstructed with the current prime: every
   rr primes. With nthreads > WIT_BATCH_MIN threads, the trace is applied to
   batches of nthreads primes at once, and an attempt in the middle of a
   batch saves no time, since the whole batch is computed anyway: the
   attempts are then made at most once per batch, with WIT_BATCH_SLACK
   primes of the batch left to check a reconstruction found (see
   verif_lifted_rational_wcoef and verif_lifted_basis), so that its
   checks do not need another batch. Most attempts fail, and they are made
   by a single thread while the others wait. */
static inline int witness_attempt_due(const gb_modpoly_t modgbs,
                                      const data_lift_t dl){
  if(dl->nthreads <= WIT_BATCH_MIN){
    return modgbs->nprimes % dl->rr == 0;
  }
  const int64_t period = MAX(dl->rr, dl->nthreads);
  return (dl->nall + WIT_BATCH_SLACK) % period == 0;
}

static void ratrecon_gb(gb_modpoly_t modgbs, data_lift_t dl,
                        mpz_t mod_p, mpz_t prod_p,
                        rrec_data_t recdata1, rrec_data_t recdata2,
                        int thrds, double *st_crt, double *st_rrec, double *st_wit,
                        const double pcost){

  double st = realtime();
  verif_lifted_rational_wcoef(modgbs, dl, thrds);
  int  b = verif_lifted_basis(modgbs, dl, thrds);
  if(b == 0){
    for(int32_t i = 0; i < modgbs->ld; i++){
      if(dl->check2[i] == 0){
        dl->lstart = i;
        break;
      }
    }
  }
  int stop = 1;
  for(int32_t i = 0; i < modgbs->ld; i++){
    if(dl->check2[i] < NBCHECK){
      stop = 0;
      break;
    }
  }
  if(stop){
    return;
  }
  *st_rrec += realtime()-st;
  if(dl->lstart != dl->start){
    dl->start = crt_ratrecon_lift_modgbs(modgbs, dl, dl->lstart, dl->start,
                                         mod_p, recdata1, recdata2,
                                         st_crt, st_rrec, pcost);
  }
  dl->lstart = dl->start;

  st = realtime();
  incremental_dlift_crt_full(modgbs, dl,
                             dl->coef, mod_p, prod_p,
                             thrds);
  *st_crt += realtime() - st;

  /********************************************************/
  /*            RATRECON WITNESS COEFFS                   */
  /********************************************************/

  if(dl->lstart == 0){
    mpz_set_ui(dl->gden, 1);
  }

  set_recdata(recdata2, mod_p);

  st = realtime();
  if(witness_attempt_due(modgbs, dl)){
    for(int32_t i = dl->lstart; i < dl->lend; i++){
      if(modgbs->modpolys[i]->len==0){
          dl->check1[i]++;
      }
      else{
        /* guess the denominator from the previous polynomial (its lifted
           denominator, or the one of its witnesses once they are
           verified), so that a spurious reconstruction never becomes the
           guess */
        if(i > 0){
          if(dl->check2[i-1] > 0 && modgbs->modpolys[i-1]->len > 0){
            mpz_set(dl->gden, modgbs->modpolys[i-1]->lm);
          }
          else if(dl->check1[i-1] > 1 && modgbs->modpolys[i-1]->len > 0){
            mpz_set(dl->gden, dl->cden[i-1]);
          }
        }
        mpz_set(dl->pgden, dl->gden);
        double rt = realtime();
        witness_crt_window(modgbs, dl, i, dl->lend);
        rt = realtime() - rt;
        *st_crt += rt;
        st += rt; /* not counted as reconstruction */
        int b = reconstructcoeff(dl, i, mod_p,
                               recdata1, recdata2);
        if(b){
          b = reconstructcoeff2(dl, i, mod_p, recdata1, recdata2);
        }
        if(b){
          b = reconstructcoeff3(dl, i, mod_p, recdata1, recdata2);
        }
        /* the numerators found by LLL are multiples of the guess times
           the coefficients, hence LLL is not tried before the modulus has
           the bitsize of the guess */
        const double lcost = lll_witnesses_cost(mpz_sizeinbase(mod_p, 2),
                                                dl->infl >= LLL_WIT_INFL);
        if(!b && (int32_t)modgbs->nprimes >= dl->lll_next
           && mpz_sizeinbase(mod_p, 2) >= mpz_sizeinbase(dl->pgden, 2) + RR_MARGIN
           && (pcost == 0 || dl->lll_cost + lcost
               <= LLL_WIT_BUDGET / 100.0 * modgbs->nprimes * pcost)){
          dl->lll_cost += lcost;
          /* with the guessed denominator, then without it, since the
             guess inflates the numerators when it shares little with the
             denominator of this polynomial */
          b = lll_witnesses(modgbs, dl, i, mod_p, dl->pgden, LLL_WIT_K);
          /* the attempt without guess only pays when the guess carries
             factors foreign to the denominator, which has not been seen
             so far on many systems (all denominators then divide about
             the same integer) */
          if(!b && mpz_cmp_ui(dl->pgden, 1) && dl->infl >= LLL_WIT_INFL){
            mpz_t one;
            mpz_init_set_ui(one, 1);
            b = lll_witnesses(modgbs, dl, i, mod_p, one, LLL_WIT_K1);
            mpz_clear(one);
          }
          dl->lll_next = b ? (int32_t)modgbs->nprimes + 1
            : lll_next_attempt(modgbs->nprimes, lcost, pcost, dl->nthreads);
        }
        if(!b){
          dl->check1[i] = 0;
          break;
        }
        else{
          dl->check1[i]++;
        }
      }
    }
  }
  *st_rrec += realtime()-st;
  *st_wit += realtime()-st;

}

long max_bit_size_gb(gb_modpoly_t modgbs){
  long nb = 0;
  for(uint32_t i = 0; i < modgbs->ld; i++){
    for(uint32_t j = 0; j < modgbs->modpolys[i]->len; j++){
      nb = MAX(nb, mpz_sizeinbase(modgbs->modpolys[i]->cf_qq[j], 2));
    }
    nb = MAX(nb, mpz_sizeinbase(modgbs->modpolys[i]->lm, 2));
  }
  return nb;
}


/*

  - returns 0 if the computation went ok

  - returns 1 in case of failure

  - returns -3 if meta data were not ok

  - returns -4 if there are too many bad primes
*/

/* Frees the data which are (re)built by a learning phase, so that a new
 * learning phase can start from scratch after a failed attempt. */
static void reset_learning_data(mstrace_t msd, md_t *st){
  for(int i = 0; i < st->nthrds; i++){
    free(msd->leadmons_current[i]);
    msd->leadmons_current[i] = NULL;
    if(i > 0){
      free(msd->leadmons_ori[i]);
      msd->leadmons_ori[i] = NULL;
      if(msd->btrace[i] != NULL){
        free_trace(&(msd->btrace[i]));
      }
    }
  }
  /* the tracer learnt with the previous prime must not be applied when
   * learning again */
  if(st->tr != NULL){
    free_trace(&(st->tr));
  }
  st->trace_level = NO_TRACER;
  /* the same holds for the DRL Groebner basis computed first for
   * elimination orders */
  free_elim_drl_data(st);
}

/* frees the bases modulo the primes of the current batch which have not
 * been used (see core_groebner_qq) and the arrays of the batch */
static void gb_free_batch(bs_t **bbs, uint64_t *bprimes, int *bbad,
                          double *bstf4, const int nthreads){
  for(int i = 0; i < nthreads; i++){
    if(bbs[i] != NULL){
      free_basis_and_only_local_hash_table_data(&(bbs[i]));
    }
  }
  free(bbs);
  free(bprimes);
  free(bbad);
  free(bstf4);
}

/* default bitsize of the primes of the multi-modular computation of
   core_groebner_qq, 31 or 52 (F4 then computes with 64 bit integers, its
   vectorized kernels need primes below 2^52): with AVX-512 IFMA, F4 costs
   less per bit of the modulus with primes of 52 bits */
#if defined HAVE_AVX512_IFMA
#define GB_PRIME_BITS_DEFAULT 52
#else
#define GB_PRIME_BITS_DEFAULT 31
#endif

/* a random prime of pbits bits */
static uint64_t gb_random_prime(const int pbits){
  if(pbits == 52){
    const uint64_t r = (((uint64_t)rand() << 31) ^ (uint64_t)rand())
      & (((uint64_t)1 << 50) - 1);
    return n_nextprime(((uint64_t)1 << 51) + r, 1);
  }
  return next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
}

/* the prime following prime among the primes of pbits bits, cyclically */
static uint64_t gb_next_prime(const uint64_t prime, const int pbits){
  if(pbits == 52){
    uint64_t p = n_nextprime(prime, 1);
    if(p >= ((uint64_t)1 << 52) - ((uint64_t)1 << 32)){
      p = n_nextprime((uint64_t)1 << 51, 1);
    }
    return p;
  }
  const uint32_t lprime = 1303905299;
  uint64_t p = next_prime(prime);
  if(p >= lprime){
    p = next_prime(1<<30);
  }
  return p;
}

/* returns 1 if prime divides a coefficient of the input, i.e. it is not
   a lucky prime (same as is_lucky_prime_ui, for primes of 64 bits) */
static int gb_is_unlucky_prime(const uint64_t prime, const bs_t * const bs){
  for(len_t i = 0; i < bs->ld; ++i){
    mpz_t *cf = bs->cf_qq[bs->hm[i][COEFFS]];
    for(len_t j = 0; j < bs->hm[i][LENGTH]; ++j){
      if(mpz_fdiv_u64(cf[j], prime) == 0){
        return 1;
      }
    }
  }
  return 0;
}

gb_modpoly_t *core_groebner_qq(
        gb_modpoly_t *modgbsp,
        bs_t *bs,
        mstrace_t msd,
        md_t *st,
        int32_t *errp,
        const len_t fc,
        const int print_gb,
        const int prime_bits
        )
{

restart:

  *errp = 0;

  int32_t info_level = st->info_level;
  int32_t truncate_lifting = st->truncate_lifting;

  if(fc == 0){
    /*generate lucky prime numbers */
    generate_lucky_primes(msd->lp, bs, st->prime_start, st->nthrds);
  }
  else{
    msd->lp->old = 0;
    msd->lp->ld = 1;
    msd->lp->p = calloc(1, sizeof(uint32_t));
    normalize_initial_basis(bs, st->gfc);
  }

  /* the trace is learned modulo a prime of 31 bits, the size of the
     next ones is chosen afterwards */
  int pbits = 31;
  uint64_t prime = 0;
  uint64_t primeinit = 0;

  prime = gb_random_prime(pbits);
  while(fc == 0 && gb_is_unlucky_prime(prime, bs)){
    prime = gb_random_prime(pbits);
  }

  primeinit = prime;
  if(fc){
    primeinit = fc;
  }
  if(!GB_PRIME_IS_WIDE(primeinit)){
    msd->lp->p[0] = primeinit;
  }
  if(info_level){
      fprintf(VERBSTREAM, "Initial prime = %" PRIu64 "\n", (uint64_t)primeinit);
  }

  int success = 1;

  int32_t maxbitsize = maxbitsize_generators(bs);

  int learn = 1, apply = 1, nprimes = 0;
  double stf4 = 0;

  rrec_data_t recdata1, recdata2;
  initialize_rrec_data(recdata1);
  initialize_rrec_data(recdata2);

  data_lift_t dlift;

  /* indicates that dlift has been already initialized */
  int dlinit = 0;

  double st_crt = 0;
  double st_rrec = 0;
  double st_wit = 0;

  uint32_t nbadprimes = 0;

  while(learn){

    int32_t is_empty = gb_modular_trace_learning((*modgbsp),
                                                 msd->mgb,
                                                 msd->num_gb, msd->leadmons_ori,
                                                 msd->btrace[0],
                                                 bs, st,
                                                 primeinit,
                                                 info_level,
                                                 print_gb,
                                                 truncate_lifting,
                                                 0,
                                                 maxbitsize,
                                                 &success);
    if(is_empty == 1 || print_gb == 1 || (*modgbsp)->ld == 0){
      if(dlinit){
        data_lift_clear(dlift);
      }
      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      return modgbsp;
    }

    apply = 1;

    gb_modpoly_realloc((*modgbsp), 1);

#ifdef DEBUGGBLIFT
    display_gbmodpoly_cf_32(ERRSTREAM, (*modgbsp));
#endif

    if(!dlinit){
      int nb = 0;
      int32_t *ldeg = array_nbdegrees((*msd->leadmons_ori), msd->num_gb[0],
                                      bs->ht->nv - st->nev, &nb);
      data_lift_init(dlift, (*modgbsp)->ld, ldeg, nb);
      dlift->nthreads = st->nthrds;
      choose_coef_to_lift((*modgbsp), dlift);
      free(ldeg);
      dlinit = 1;
    }

    if(info_level){
      fprintf(VERBSTREAM,"\n---------- COMPUTATIONAL DATA -----------\n");
      int s= 0;
      for(int i = 0; i < dlift->nsteps; i++){
        fprintf(VERBSTREAM, "[%d]", dlift->steps[i]);
        fflush(VERBSTREAM);
        s+=dlift->steps[i];
      }
      fprintf(VERBSTREAM, "\n");
      fprintf(VERBSTREAM,
	      "#polynomials to lift %14lu\n",
	      (unsigned long) s);
      fprintf(VERBSTREAM, "-----------------------------------------\n");
    }

    if(is_empty || success == 0 || fc) {

      apply = 0;
      if(dlinit){
        data_lift_clear(dlift);
      }

      gb_modpoly_without_hash_table_clear((*modgbsp));
      reset_learning_data(msd, st);

      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      fprintf(VERBSTREAM, "Something went wrong in the learning phase, msolve restarts.");
      /* return core_groebner_qq(modgbsp, bs, msd, st, errp, fc, print_gb); */
      goto restart;
    }
    /* duplicate data for multi-threaded multi-mod computation */
    duplicate_data_mthread_gbtrace(st->nthrds, bs, st, msd->num_gb,
                                   msd->leadmons_ori, msd->leadmons_current,
                                   msd->btrace);


    learn = 0;
    if(fc == 0){
      pbits = (prime_bits == 31 || prime_bits == 52) ? prime_bits
        : GB_PRIME_BITS_DEFAULT;
    }
    if(pbits == 52){
      gb_modpoly_widen((*modgbsp));
    }
    if(info_level){
      fprintf(VERBSTREAM, "Primes of %d bits\n", pbits);
    }
    prime = gb_random_prime(pbits);
    if(info_level){
        fprintf(VERBSTREAM, "New prime = %" PRIu64 "\n", (uint64_t)prime);
    }
    /* The trace is applied to batches of nthreads primes in parallel (one
     * prime per thread); the images are then copied and the lifting steps
     * applied prime by prime, in the order of the primes, exactly as with
     * one thread, so that the primes used and the output do not depend on
     * the number of threads (only the end of the last batch is wasted). */
    const int nthreads = st->nthrds;
    uint64_t *bprimes = (uint64_t *)malloc(nthreads * sizeof(uint64_t));
    bs_t **bbs = (bs_t **)calloc(nthreads, sizeof(bs_t *));
    int *bbad = (int *)calloc(nthreads, sizeof(int));
    double *bstf4 = (double *)calloc(nthreads, sizeof(double));
    int bpos = nthreads;
    while(apply){

      if(bpos == nthreads){
        /* generate lucky prime numbers */
        for(int i = 0; i < nthreads; i++){
          prime = gb_next_prime(prime, pbits);
          while(gb_is_unlucky_prime(prime, bs) || prime==primeinit){
            prime = gb_next_prime(prime, pbits);
          }
          bprimes[i] = prime;
        }
        /* F4 is run by each thread on a single thread */
        st->info_level  = 0;
        st->f4_qq_round = 2;
        st->nthrds = 1;
#pragma omp parallel for num_threads(nthreads) schedule(static, 1)
        for(int i = 0; i < nthreads; i++){
          const double rt = realtime();
          bbs[i] = gb_modular_trace_apply((*modgbsp), msd->num_gb[i],
                                          msd->leadmons_ori[i],
                                          msd->leadmons_current[i],
                                          bs, st, bprimes[i], bbad + i);
          bstf4[i] = realtime() - rt;
        }
        st->nthrds = nthreads;
        bpos = 0;
      }
      const int bi = bpos++;
      prime = bprimes[bi];
      stf4 = bstf4[bi];
      msd->bad_primes[0] = bbad[bi];
      int nthrds = 1; /* the lifting steps are applied prime by prime */

      if((*modgbsp)->alloc <= nprimes + 2){
        /* Room for the images modulo the next primes, which are copied one
         * prime at a time. The increment is kept small (at most ~1/16 of
         * the current number of primes, but at least 4) since the unused
         * slots are allocated for every term of the basis which is not
         * lifted yet, i.e. for most of the basis at the end of the
         * computation, when memory consumption peaks. */
        gb_modpoly_realloc((*modgbsp), MAX(4, (uint32_t)nprimes / 16));
      }

      if(bbs[bi] != NULL){
        /* copy of data for multi-mod computation */
        modpgbs_set((*modgbsp), bbs[bi], bbs[bi]->ht, prime, dlift->S, st->nev);
        free_basis_and_only_local_hash_table_data(&(bbs[bi]));
      }

      /* nprimes += st->nthrds; */
      nprimes += 1; /* at the moment, multi-mod comp is mono-threaded */
      dlift->nall = nprimes;
      if(nprimes == 1){
        /* if(info_level>1){ */
        /*   fprintf(ERRSTREAM, "Application phase %.2f Gops/sec\n", */
        /*           (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4)); */
        /*   fprintf(ERRSTREAM, "Elapsed time: %.2f\n", stf4); */
        /* } */
	  if(info_level){
	    fprintf(VERBSTREAM,
		    "\n---------------- TIMINGS ----------------\n");
	    fprintf(VERBSTREAM,
		    "multi-mod overall(elapsed) %9.2f sec\n",
		    stf4);
	    /* no trace is learned with probabilistic linear algebra */
	    if (info_level > 1 && st->tr != NULL){
	      /* multiply-adds by the known pivots, recorded in the trace,
	       * for the prime of this step */
	      fprintf(VERBSTREAM,
		      "application phase          %9.2f Gops/sec\n",
		      2.0*st->tr->nmadd/1000.0/1000.0/1000.0/(stf4));
	    }
	    fprintf(VERBSTREAM,
		    "-----------------------------------------\n");
        }
      if (info_level) {
	      fprintf(VERBSTREAM,
		      "\nmulti-modular steps\n");
	      fprintf(VERBSTREAM, "-------------------------------------------------\
-----------------------------------------------------\n");
        }
      }
      int bad = 0;
      for(int i = 0; i < nthrds/* st->nthrds */; i++){
        if(msd->bad_primes[i] == 1){
            bad = 1;
            if(info_level > 1){
                fprintf(VERBSTREAM, "[!]");
            }
            nbadprimes++;
            msd->bad_primes[i] = 0;
        }
        if(msd->bad_primes[i] == 2){
            bad = 1;
            if(info_level > 1){
                fprintf(VERBSTREAM, "[!!]");
            }
            /* nbadprimes = nprimes; */
            nbadprimes++;
            msd->bad_primes[i] = 0;
        }
      }

      if(nbadprimes >= nprimes){
        if(info_level){
          fprintf(VERBSTREAM, "\nToo many bad primes, computation will restart\n");
          fprintf(VERBSTREAM, "-------------------------------------------------\
-----------------------------------------------------\n");
        }
        if(dlinit){
          data_lift_clear(dlift);
        }

        gb_modpoly_without_hash_table_clear((*modgbsp));
        reset_learning_data(msd, st);
        free_rrec_data(recdata1);
        free_rrec_data(recdata2);
        /* reset info_level to previous setting */
        st->info_level = info_level;
        gb_free_batch(bbs, bprimes, bbad, bstf4, nthreads);
        /* return core_groebner_qq(modgbsp, bs, msd, st, errp, fc, print_gb); */
        goto restart;

      }

      int lstart = dlift->lstart;
      st_wit = 0;

      if(!bad){
        ratrecon_gb((*modgbsp), dlift, msd->mod_p, msd->prod_p, recdata1, recdata2,
                    nthrds/* st->nthrds */, &st_crt, &st_rrec, &st_wit,
                    lll_prime_cost(st));
      }
      /* the period depends only on the number of primes (not on timings)
         so that the number of primes used is reproducible */
      int32_t rr = 1;
      while(2 * rr <= nprimes / RR_PERIOD_DIV){
        rr *= 2;
      }
      if(rr != dlift->rr){
        dlift->rr = rr;
        if(info_level){
          fprintf(VERBSTREAM, "(->%d)", dlift->rr);
	        fflush(VERBSTREAM);
        }
      }
      if(info_level){
        if(!(nprimes & (nprimes - 1))){
          fprintf(VERBSTREAM, "{%d}", nprimes);
	      fflush(VERBSTREAM);
        }
      }
      apply = 0;
      for(len_t i = 0; i < (*modgbsp)->ld; i++){
        if(dlift->check2[i] < NBCHECK){
          apply = 1;
          break;
        }
      }
      if(!bad && print_gb == 1){
          apply = 0;
      }
      if(dlift->lstart != lstart){
        if(info_level){
          fprintf(VERBSTREAM, "<%.2f%%>", 100* (float)MIN((dlift->lstart + 1), (*modgbsp)->ld)/(*modgbsp)->ld);
	      fflush(VERBSTREAM);
        }
        lstart = dlift->lstart;
      }
      /* this is where learn could be reset to 1 */
      /* but then duplicated datas and others should be free-ed */
    }
    gb_free_batch(bbs, bprimes, bbad, bstf4, nthreads);
    if (info_level){
      fprintf(VERBSTREAM, " \n-------------------------------------------------\
-----------------------------------------------------\n");
    }
  }
  if(info_level){
    long nbits = max_bit_size_gb((*modgbsp));
    fprintf(VERBSTREAM,"\n\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(VERBSTREAM, "Max coeff. bitsize %16lu\n", (unsigned long) nbits);
    fprintf(VERBSTREAM, "#primes            %16lu\n", (unsigned long) nprimes);
    fprintf(VERBSTREAM, "#bad primes        %16lu\n", (unsigned long) nbadprimes);
    fprintf(VERBSTREAM, "-----------------------------------------\n");
    fprintf(VERBSTREAM, "\n---------------- TIMINGS ----------------\n");
    fprintf(VERBSTREAM, "CRT     (elapsed)         %10.2f sec\n", st_crt);
    fprintf(VERBSTREAM, "ratrecon(elapsed)         %10.2f sec\n", st_rrec);
    if(dlinit && dlift->lll_ntries > 0){
      fprintf(VERBSTREAM, "  incl. LLL, %4d calls (%d found) %6.2f sec\n",
              dlift->lll_ntries, dlift->lll_nsucc, dlift->lll_time);
    }
    fprintf(VERBSTREAM, "-----------------------------------------\n");
  }

  if(dlinit){
    data_lift_clear(dlift);
  }

  free_rrec_data(recdata1);
  free_rrec_data(recdata2);

  return modgbsp;
}

uint64_t export_results_from_groebner_qq(
        /* return values */
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        void *(*mallocp) (size_t),
        const int32_t elim_block_len,
        const int32_t print_gb, /* if this is 1, only leading monomials are exported */
        gb_modpoly_t gb
        )
{
    int64_t nelts = gb->ld;
    /* Over QQ, if eliminating, these variables are no longer included in */
    /* gb->nv, thus we need to help us with adding elim_block_len */
    /* correspondingly when exporting the basis. */
    int32_t nv    = gb->nv;
    int32_t nve   = gb->bht->nv;

    *bld  = nelts;

    int32_t *lens = (int32_t *)(*mallocp)(((uint64_t)nelts) * sizeof(int32_t));

    int64_t nterms = 0;

    for(int32_t i = 0; i < nelts; i++){
        int32_t nlen = 0;
        if (print_gb == 2) {
            int32_t len = gb->modpolys[i]->len;

            for(int32_t j = len-1; j >= 0; j--){
                if (mpz_cmp_ui(gb->modpolys[i]->cf_qq[j], 0) != 0) {
                    nlen++;
                }
            }
        }
        nlen++; /* to take into account the lm */
        lens[i] = nlen;
        nterms += nlen;
    }

    int32_t *exp = (int32_t *)(*mallocp)(
            ((uint64_t)nterms) * ((uint64_t) nve) * sizeof(int32_t));

    memset(exp, 0, (unsigned long)nterms * nve * sizeof(int32_t));

    mpz_t *cf_qq = (mpz_t *)(*mallocp)(
            nterms *sizeof(mpz_t));
    for(int64_t i = 0; i < nterms; i++){
        mpz_init(cf_qq[i]);
    }

    hm_t *hm = NULL;
    ht_t *ht = gb->bht;
    int *evi = NULL;
    if (print_gb == 2) {
        const len_t ebl = ht->ebl;
        const len_t evl = ht->evl;
        evi = (int *)malloc((unsigned long)(ht->nv + 1) * sizeof(int));
        if (ebl == 0) {
            for (len_t i = 1; i <= evl; ++i) {
                evi[i-1]    =   i;
            }
        } else {
            for (len_t i = 1; i < ebl; ++i) {
                evi[i-1]    =   i;
            }
            for (len_t i = ebl + 1; i < evl; ++i) {
                evi[i - 2] = i;
            }
        }
    }

    int64_t term = 0;
    for(int64_t p = 0; p < nelts; p++){
        hm  = gb->hm[p]+OFFSET;
        int32_t l = gb->modpolys[p]->len;
        for(int32_t n = 0; n < nv; n++){
            exp[term * nve + n + elim_block_len] = gb->ldm[p * nv + n];
        }
        mpz_set(cf_qq[term], gb->modpolys[p]->lm);

        term++;
        if (print_gb == 2) {
            for(int32_t i = l-1; i >= 0; i--) {
                if(mpz_cmp_ui(gb->modpolys[p]->cf_qq[i], 0) != 0) {
                  for (int32_t n = 0; n < nve; n++) {
                    exp[term * nve + n] = ht->ev[hm[l - i]][evi[n]];
                  }
                  mpz_set(cf_qq[term], gb->modpolys[p]->cf_qq[i]);

                  term++;
                }
            }
        }
    }

    *blen = lens;
    *bexp = exp;
    *bcf  = cf_qq;

    free(evi);
    return nterms;
}

gb_modpoly_t *groebner_qq(
        data_gens_ff_t *gens,
        msflags_t flags)
{
  int32_t ht_size = flags->ht_size;
  int32_t nr_threads = flags->nr_threads;
  int32_t max_nr_pairs = flags->max_nr_pairs;
  int32_t elim_block_len = flags->elim_block_len;
  int32_t reset_ht = flags->reset_ht;
  int32_t la_option = flags->la_option;
  const int32_t use_signatures = flags->use_signatures;
  int32_t info_level = flags->info_level;
  int32_t pbm_file = flags->pbm_file;
  int32_t print_gb = flags->print_gb;
  int32_t truncate_lifting = flags->truncate_lifting;
  int mon_order = 0;

  /* input data */
  uint32_t field_char = gens->field_char;
  const void *cfs = gens->mpz_cfs;
  if(gens->field_char){
    cfs = gens->cfs;
  }
  else{
    cfs = gens->mpz_cfs;
  }
  int32_t *exps = gens->exps;
  int32_t *lens = gens->lens;

  int32_t nr_vars = gens->nvars;
  int32_t nr_gens = gens->ngens;
  int reduce_gb = 1;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* data structures for basis, hash table and statistics */
  bs_t *bs  = NULL;
  ht_t *bht = NULL;
  md_t *md  = NULL; /* meta data */

  int success = 0;

  /*msolve trace data */
  mstrace_t msd;
  success = initialize_gba_input_data(&bs, &bht, &md,
          lens, exps, cfs, field_char, mon_order, elim_block_len,
          nr_vars, nr_gens, 0 /* # normal forms */, ht_size,
          nr_threads, max_nr_pairs, reset_ht, la_option, use_signatures,
          reduce_gb, pbm_file, truncate_lifting, info_level);

  /* all input generators are invalid */
  if (success == -1) {
      return NULL;
  }
  if (success == 0) {
      fprintf(ERRSTREAM,"Bad input data, stopped computation.\n");
      exit(1);
  }

  md->elim_drl = gens->elim_drl;
  initialize_mstrace(msd, md, bs);
  int err = 0;

  gb_modpoly_t *modgbsp = malloc(sizeof(gb_modpoly_t));
  modgbsp = core_groebner_qq(modgbsp, bs, msd, md, &err, field_char, print_gb,
                             gens->gb_prime_bits);
  if (err) {
      fprintf(ERRSTREAM,"Problem with groebner_qq, stopped computation (%d).\n", err);
      exit(1);
  }
  free_elim_drl_data(md);

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->f4_ctime = ct1 - ct0;
  md->f4_rtime = rt1 - rt0;

  get_and_print_final_statistics(VERBSTREAM, md, bs);

  /* free and clean up */
  free_mstrace(msd, md);
  free_basis_without_hash_table(&bs);
  free_trace(&md->tr);
  free(md);
  md    = NULL;

  return modgbsp;

}

/*
 Function which is called by core_msolve
 */

void print_msolve_gbtrace_qq(data_gens_ff_t *gens,
                            msflags_t flags)
{
  gb_modpoly_t *modgbsp;

  modgbsp = groebner_qq(gens, flags);

  output_timer_t ot;
  output_timer_start(&ot);
  FILE *ofile;
  if (flags->files->out_file != NULL) {
    ofile = fopen(flags->files->out_file, "wb+");
  } else {
    ofile = OUTSTREAM;
  }
  if (flags->print_gb == 1) {
    fprintf(ofile, "#Leading ideal data\n");
  } else {
    if (flags->print_gb > 1) {
      if(flags->truncate_lifting>0){
        fprintf(ofile, "#Truncated reduced Groebner basis data\n");
      }
      else{
        fprintf(ofile, "#Reduced Groebner basis data\n");
      }
    }
  }
  fprintf(ofile, "#---\n");
  fprintf(ofile, "#field characteristic: 0\n");
  fprintf(ofile, "#variable order:       ");
  for (int i = gens->elim; i < gens->nvars-1; ++i) {
    fprintf(ofile, "%s, ", gens->vnames[i]);
  }
  fprintf(ofile, "%s\n", gens->vnames[gens->nvars-1]);
  fprintf(ofile, "#monomial order:       graded reverse lexicographical\n");
  if ((*modgbsp)->ld == 1) {
    fprintf(ofile, "#length of basis:      1 element\n");
  } else {
    fprintf(ofile, "#length of basis:      %u elements sorted by increasing leading monomials\n", (*modgbsp)->ld);
  }
  fprintf(ofile, "#---\n");
  if (flags->files->out_file != NULL) {
    fclose(ofile);
  }

  if(flags->print_gb > 1){

    if(flags->files->out_file != NULL){
      FILE *ofile = fopen(flags->files->out_file, "ab+");
      display_gbmodpoly_cf_qq(ofile, (*modgbsp), gens, flags->nr_threads);
      fclose(ofile);
    }
    else{
      display_gbmodpoly_cf_qq(OUTSTREAM, (*modgbsp), gens, flags->nr_threads);
    }
  }
  if(flags->print_gb == 1){
    if(flags->files->out_file != NULL){
      FILE *ofile = fopen(flags->files->out_file, "ab+");
      display_lm_gbmodpoly_cf_qq(ofile, (*modgbsp), gens);
      fclose(ofile);
    }
    else{
      display_lm_gbmodpoly_cf_qq(OUTSTREAM, (*modgbsp), gens);
    }
  }
  output_timer_stop(flags->files, &ot);
  gb_modpoly_clear((*modgbsp));
  free(modgbsp);
}

/* we get from julia the generators as three arrays:
 * 0.  a pointer to an int32_t array for returning the basis to julia
 * 1.  an array of the lengths of each generator
 * 2.  an array of all coefficients of all generators in the order:
 *     first all coefficients of generator 1, then all of generator 2, ...
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *
 *  RETURNs the length of the jl_basis array */
int64_t export_groebner_qq(
                void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t print_gb,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t truncate_lifting,
        const int32_t info_level
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* data structures for basis, hash table and statistics */
    bs_t *bs  = NULL;
    ht_t *bht = NULL;
    md_t *md  = NULL; /* meta data */

    int success = 0;

    /*msolve trace data*/
    mstrace_t msd;
    success = initialize_gba_input_data(&bs, &bht, &md,
            lens, exps, cfs, field_char, mon_order, elim_block_len,
            nr_vars, nr_gens, 0 /* # normal forms */, ht_size,
            nr_threads, max_nr_pairs, reset_ht, la_option, 0 /*use_signatures*/,
            reduce_gb, pbm_file, truncate_lifting, info_level);

    /* all input generators are invalid */
    if (success == -1) {
        return_zero(bld, blen, bexp, bcf, nr_vars, field_char, mallocp);
        return 1;
    }
    if (success == 0) {
        fprintf(ERRSTREAM, "Bad input data, stopped computation.\n");
        exit(1);
    }

    initialize_mstrace(msd, md, bs);
    int err = 0;

    gb_modpoly_t *modgbsp = malloc(sizeof(gb_modpoly_t));
    modgbsp =
    core_groebner_qq(modgbsp, bs, msd, md, &err, field_char,
                     print_gb, 0);
    if (err) {
        fprintf(ERRSTREAM, "Problem with groebner_qq, stopped computation.\n");
        exit(1);
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->f4_ctime = ct1 - ct0;
    md->f4_rtime = rt1 - rt0;

    get_and_print_final_statistics(VERBSTREAM, md, bs);

    int64_t nterms  = export_results_from_groebner_qq(bld, blen, bexp,
                                                      bcf, mallocp, elim_block_len, print_gb, (*modgbsp));

    /* free and clean up */
    free_mstrace(msd, md);
    free_basis_without_hash_table(&bs);
    free_trace(&md->tr);
    free(md);
    md = NULL;
    gb_modpoly_clear((*modgbsp));
    free(modgbsp);

    return nterms;

}
