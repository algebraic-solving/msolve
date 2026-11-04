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
#include<flint/fmpq.h>
#include "streams.h"

/* rational reconstruction of u mod m with bounds N, D through FLINT
   (half-gcd based) instead of the quadratic Euclidean algorithm;
   u is assumed to be reduced modulo m */
static int flint_rrec(mpz_t n, mpz_t d, const mpz_t u, const mpz_t m,
                      const mpz_t N, const mpz_t D){
  fmpz_t fu, fm, fN, fD;
  fmpq_t q;
  fmpz_init(fu); fmpz_init(fm); fmpz_init(fN); fmpz_init(fD);
  fmpq_init(q);
  fmpz_set_mpz(fu, u);
  fmpz_set_mpz(fm, m);
  fmpz_set_mpz(fN, N);
  fmpz_set_mpz(fD, D);
  int b = fmpq_reconstruct_fmpz_2(q, fu, fm, fN, fD);
  if(b){
    fmpz_get_mpz(n, fmpq_numref(q));
    fmpz_get_mpz(d, fmpq_denref(q));
  }
  fmpq_clear(q);
  fmpz_clear(fu); fmpz_clear(fm); fmpz_clear(fN); fmpz_clear(fD);
  return b;
}

/* FLINT requires 2ND < m */
static int flint_rrec_bounds_ok(const mpz_t m, rrec_data_t rd){
  mpz_mul(rd->tmp, rd->N, rd->D);
  mpz_mul_2exp(rd->tmp, rd->tmp, 1);
  return mpz_cmp(rd->tmp, m) < 0;
}

/* same as ratrecon, through FLINT when the bounds allow it */
static int lift_ratrecon(mpz_t n, mpz_t d, mpz_t u, const mpz_t mod,
                         rrec_data_t rd){
  if(!flint_rrec_bounds_ok(mod, rd)){
    return ratrecon(n, d, u, mod, rd);
  }
  mpz_mod(u, u, mod);
  return flint_rrec(n, d, u, mod, rd->N, rd->D);
}

/* same as ratreconwden, through FLINT when the bounds allow it */
static int lift_ratreconwden(mpz_t n, mpz_t d, mpz_t u, const mpz_t mod,
                             const mpz_t gden, rrec_data_t rd){
  if(!flint_rrec_bounds_ok(mod, rd)){
    return ratreconwden(n, d, u, mod, gden, rd);
  }
  mpz_mod(u, u, mod);
  mpz_mul(rd->r1, u, gden);
  mpz_mod(rd->r1, rd->r1, mod);
  return flint_rrec(n, d, rd->r1, mod, rd->N, rd->D);
}

#define NBCHECK 2
/* number of bits of slack kept in rational reconstruction bounds */
#define RR_MARGIN 16
/* witnesses are reconstructed every rr primes, rr being the largest
   power of 2 <= nprimes / RR_PERIOD_DIV */
#define RR_PERIOD_DIV 128
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
  /* len integer numerators of the lifted coefficients; they all share the
   * common denominator lm (hence the polynomial is lm * LM + sum cf_qq[j]*m_j
   * up to the order of the terms). */
  mpz_t *cf_qq;
  mpz_t lm; /* stores the leading coefficient (hence all numerators of cf_qq should be divided by lm) */
} modpolys_struct;

typedef modpolys_struct modpolys_t[1];

/* coefficient of index j of polynomial p modulo the k-th prime */
#define MODCF(p, j, k) ((p)->cf_32[(size_t)(k) * (p)->len + (j)])

typedef struct {
  uint32_t alloc; /* alloc -> max number of primes */
  uint32_t nprimes; /* number of primes */
  mp_limb_t *primes; /* array of prime numbers encoded with uint64_t to ensure
                       compatibility with flint */
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


typedef struct{
  int32_t npol; /* number of polynomials to be lifted */
  int32_t rr; /* number of primes before activating rational reconstruction */
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
} data_lift_struct;

typedef data_lift_struct data_lift_t[1];


static inline void data_lift_init(data_lift_t dl,
                                  int32_t npol,
                                  int32_t *steps, int32_t nsteps){
  dl->npol = npol;
  dl->rr = 1;
  dl->lstart = 0;
  dl->nsteps = nsteps;
  dl->S = 0;
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

}


static inline void gb_modpoly_init(gb_modpoly_t modgbs,
                                   uint32_t alloc, int32_t *lens,
                                   bs_t *bs,
                                   int nv, uint32_t ld,
                                   int32_t *lm,
                                   md_t *st){
  modgbs->alloc = alloc;
  modgbs->nprimes = 0;
  modgbs->primes = (mp_limb_t *)calloc(alloc, sizeof(mp_limb_t));
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
    if(lens[i] > 0){
      modgbs->modpolys[i]->cf_32 = (uint32_t *)calloc((size_t)lens[i] * alloc,
                                                      sizeof(uint32_t));
      modgbs->modpolys[i]->cf_qq = (mpz_t *)malloc(sizeof(mpz_t)*lens[i]);
    }
    else{
      modgbs->modpolys[i]->cf_32 = NULL;
      modgbs->modpolys[i]->cf_qq = NULL;
    }
    for(uint32_t j = 0; j < lens[i]; j++){
      mpz_init(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_init(modgbs->modpolys[i]->lm);
    mpz_set_ui(modgbs->modpolys[i]->lm, 1);
  }
}


/* Makes room for newalloc more primes. Polynomials whose modular data have
 * already been free-ed (cf_32 == NULL) are skipped. */
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
    if(modgbs->modpolys[i]->cf_32 == NULL || len == 0){
      continue;
    }
    uint32_t *newcfs_pol = (uint32_t *)realloc(modgbs->modpolys[i]->cf_32,
                                               len * modgbs->alloc * sizeof(uint32_t));
    if(newcfs_pol == NULL){
      fprintf(ERRSTREAM, "Problem when reallocating modgbs (cfs_pol)\n");
      exit(1);
    }
    /* prime-major storage: the new primes are appended at the end */
    memset(newcfs_pol + len * oldalloc, 0,
           len * (modgbs->alloc - oldalloc) * sizeof(uint32_t));
    modgbs->modpolys[i]->cf_32 = newcfs_pol;
  }
}


static inline void display_gbmodpoly_cf_32(FILE *file,
                                     gb_modpoly_t modgbs){
  fprintf(file, "alloc = %d\n", modgbs->alloc);
  fprintf(file, "nprimes = %d\n", modgbs->nprimes);
  fprintf(file, "primes = [");
  for(uint32_t i = 0; i < modgbs->alloc-1; i++){
    fprintf(file, "%lu, ", (unsigned long)modgbs->primes[i]);
  }
  fprintf(file, "%lu]\n", (unsigned long)modgbs->primes[modgbs->alloc -1]);
  fprintf(file, "numpolys = %d\n", modgbs->ld);
  fprintf(file, "[\n");
  for(uint32_t i = 0; i < modgbs->ld; i++){
    uint32_t len = modgbs->modpolys[i]->len;
    fprintf(file, "[%d, ", len);
    if(modgbs->modpolys[i]->cf_32 == NULL){
      fprintf(file, "(lifted)],\n");
      continue;
    }
    for(uint32_t j = 0; j < len; j++){
      fprintf(file, "[");
      for(uint32_t k = 0; k < modgbs->alloc-1; k++){
        fprintf(file, "%d, ", MODCF(modgbs->modpolys[i], j, k));
      }
      if(j < len - 1){
        fprintf(file, "%d], ", MODCF(modgbs->modpolys[i], j, modgbs->alloc-1));
      }
      else{
        fprintf(file, "%d]\n", MODCF(modgbs->modpolys[i], j, modgbs->alloc-1));
      }
    }
    fprintf(file, "],\n");
  }
  fprintf(file, "]:\n");
}

static inline void display_modpoly(FILE *file,
                                   gb_modpoly_t modgbs,
                                   int32_t pos,
                                   data_gens_ff_t *gens){

  if(modgbs->modpolys[pos]->len == 0){
    display_monomial_single(file, gens, pos, &modgbs->ldm);
    return;
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->lm, 1) != 0){
    mpz_out_str(file, 10, modgbs->modpolys[pos]->lm);
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
        mpz_out_str(file, 10, cf[i]);
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
    mpz_out_str(file, 10, cf[0]);
  }
  else{
    if(mpz_cmp_si(cf[0], -1) == 0){
        fprintf(file, "-");
    }
    int star = 0;
    //coef is neither 1 nor -1
    if(mpz_cmp_ui(cf[0], 1) != 0 && mpz_cmp_si(cf[0], -1) != 0){
        mpz_out_str(file, 10, cf[0]);
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
                                           data_gens_ff_t *gens){
  int32_t p = modgbs->ld ;
  if(p==0){
    fprintf(file, "[0]:\n");
    return;
  }
  fprintf(file, "[");
  for(int i = 0; i < p-1; i++){
    display_modpoly(file, modgbs, i, gens);
    fprintf(file, ", \n");
  }
  display_modpoly(file, modgbs, p-1, gens);
  fprintf(file, "]:\n");
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
  free(modgbs->cf_64);
  modgbs->cf_64 = NULL;
  free(modgbs->ldm);
  modgbs->ldm = NULL;
  for(uint32_t i = 0; i < modgbs->ld; i++){
    /* cf_32 is NULL if already free-ed after lifting */
    free(modgbs->modpolys[i]->cf_32);
    modgbs->modpolys[i]->cf_32 = NULL;
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
                              const int32_t fc,
                              int32_t start, const long elim){
  if(modgbs->nprimes >= modgbs->alloc-1){
    fprintf(ERRSTREAM, "Not enough space in modgbs\n");
    exit(1);
  }
  modgbs->primes[modgbs->nprimes] = fc;

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
    if(modgbs->modpolys[i]->cf_32 == NULL){
      /* already lifted, modular data are no longer needed */
      continue;
    }
    int32_t bc = modgbs->modpolys[i]->len - 1;
    const uint32_t * const cfs = bs->cf_32[bs->hm[idx][COEFFS]];
    for (j = 1; j < len; ++j) {
      MODCF(modgbs->modpolys[i], bc, modgbs->nprimes) = cfs[j];
      bc--;
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

/* returns 1 if only the basis of the elimination ideal is lifted, i.e. the
 * eliminated variables are dropped, else returns 0 */
static inline int lift_elim_ideal_only(const md_t *st){
  return st->nev > 0 && !st->elim_full_basis;
}

/* returns the number of variables of the lifted leading monomials */
static inline int lifted_nvars(const md_t *st){
  return lift_elim_ideal_only(st) ? st->nvars - st->nev : st->nvars;
}

static int32_t gb_modular_trace_learning(gb_modpoly_t modgbs,
                                           int32_t *mgb,
                                           int32_t *num_gb,
                                           int32_t **leadmons,
                                           trace_t *trace,
                                           bs_t *bs_qq,
                                           md_t *st,
                                           const int32_t fc,
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

    if(info_level > 2){
        fprintf(VERBSTREAM, "------------------------------------------\n");
        fprintf(VERBSTREAM, "#ADDITIONS       %13lu\n", (unsigned long)st->trace_nr_add * 1000);
        fprintf(VERBSTREAM, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->trace_nr_mult * 1000);
        fprintf(VERBSTREAM, "#REDUCTIONS      %13lu\n", (unsigned long)st->trace_nr_red);
        fprintf(VERBSTREAM, "------------------------------------------\n");
    }

    /* Leading monomials from Grobner basis */
    int32_t *bexp_lm = get_lm_from_bs(bs, bht);
    /* may be non NULL if the computation has been restarted */
    free(leadmons[0]);
    leadmons[0] = bexp_lm;

    int32_t len = bs->lml;
    const int nvl = lifted_nvars(st);
    num_gb[0] = compute_num_gb(bexp_lm, len, bht->nv,
                               lift_elim_ideal_only(st) ? st->nev : 0);
    int32_t *bexp_lm2 = NULL;
    if(lift_elim_ideal_only(st)){
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

    int32_t *lens = array_of_lengths(leadmons[0], num_gb[0], bs, nvl);

    if(truncate_lifting != 0 && truncate_lifting < num_gb[0]){
      gb_modpoly_init(modgbs, 2, lens, bs, nvl, truncate_lifting, leadmons[0], st);
    }
    else{
      gb_modpoly_init(modgbs, 2, lens, bs, nvl, num_gb[0], leadmons[0], st);
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


static void gb_modular_trace_application(gb_modpoly_t modgbs,
                                         int32_t *mgb,
                                         int32_t *num_gb,
                                         int32_t **leadmons_ori,
                                         int32_t **leadmons_current,
                                         trace_t **btrace,
                                         bs_t *bs_qq,
                                         md_t *st,
                                         int info_level,
                                         primes_t *lp,
                                         int32_t start,
                                         double *stf4,
                                         int *bad_primes)
{

  double rt = realtime();
  st->info_level  = 0;
  st->f4_qq_round = 2;
  /* tracing phase */

  /* F4 and FGLM are run using a single thread */
  /* st->nthrds is reset to its original value afterwards */
  /*at the moment multi-threading is not supprted here*/
  memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));

  /* current load of the hash table */
  /* hl_t eld = bs_qq->ht->eld; */
  bs_t *bs = NULL;
  int32_t error = 0;
  bs = core_gba(bs_qq, st, &error, lp->p[0]);
  *stf4 = realtime()-rt;
  /* core_gba returns NULL (and has already free-ed its data) when the
   * prime has been detected as bad while applying the tracer */
  if (bs == NULL || error) {
      bad_primes[0] = 1;
      if (bs != NULL) {
          free_basis_and_only_local_hash_table_data(&bs);
      }
      return;
  }
  ht_t **bht = &(bs->ht);
  if(bs->lml < modgbs->ld){
      bad_primes[0] = 2;
      free_basis_and_only_local_hash_table_data(&bs);
      return;
  }
  for(len_t i = 0; i < modgbs->ld; i++){
        if(modgbs->modpolys[i]->len < bs->hm[bs->lmps[i]][LENGTH]-1){
            bad_primes[0] = 2;
            free_basis_and_only_local_hash_table_data(&bs);
            return;
        }
        if(modgbs->modpolys[i]->len > bs->hm[bs->lmps[i]][LENGTH]-1){
            bad_primes[0]=1;
            free_basis_and_only_local_hash_table_data(&bs);
            return;
        }
  }
  int32_t lml = bs->lml;
  if (lift_elim_ideal_only(st)) {
      int32_t j = 0;
      for (len_t i = 0; i < bs->lml; ++i) {
          if ((*bht)->ev[bs->hm[bs->lmps[i]][OFFSET]][0] == 0) {
              bs->lm[j]   = bs->lm[i];
              bs->lmps[j] = bs->lmps[i];
              ++j;
          }
      }
      lml = j;
  }
  if (lml != num_gb[0]) {
      if (bs != NULL) {
          bad_primes[0] = 2;
          free_basis_and_only_local_hash_table_data(&bs);
      }
      return;
  }

  if(lift_elim_ideal_only(st)){
    get_lm_from_bs_trace_elim(bs, bht[0], leadmons_current[0], num_gb[0]);
  }
  else{
    get_lm_from_bs_trace(bs, bht[0], leadmons_current[0]);
  }

  if(!equal_staircase(leadmons_current[0], leadmons_ori[0],
                      num_gb[0], num_gb[0], lifted_nvars(st))){
    bad_primes[0] = 1;
  }

  if(!bad_primes[0] && bs != NULL){
    /* copy of data for multi-mod computation */
    modpgbs_set(modgbs, bs, bht[0], lp->p[0], start, st->nev);
  }

  if (bs != NULL) {
    free_basis_and_only_local_hash_table_data(&bs);
  }

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
  }
}


/* Incremental CRT (called once FLINT multi_CRT has been called) */
/* mod is the current modulus */
static inline void incremental_dlift_crt(gb_modpoly_t modgbs, data_lift_t dlift,
                                         int32_t *coef, mpz_t *mod_p, mpz_t *prod_p,
                                         int thrds){

  /* all primes are assumed to be good primes */
  mpz_mul_ui(prod_p[0], mod_p[0], modgbs->primes[modgbs->nprimes - 1 ]);
  for(int32_t k = dlift->lstart; k <= dlift->lend; k++){
      if(modgbs->modpolys[k]->len){
        uint32_t c = MODCF(modgbs->modpolys[k], coef[k], modgbs->nprimes  - 1);
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
  mpz_mul_ui(prod_p, mod_p, (uint32_t)newprime);
  for(int32_t k = dl->S; k < dl->end; k++){
      if(modgbs->modpolys[k]->len){
        uint64_t c = MODCF(modgbs->modpolys[k], coef[k], modgbs->nprimes  - 1);
        mpz_CRT_ui(dl->crt[k], dl->crt[k], mod_p,
                   c, newprime, prod_p, dl->tmp, 1);
        c = MODCF(modgbs->modpolys[k], dl->coef2[k], modgbs->nprimes  - 1);
        mpz_CRT_ui(dl->crt2[k], dl->crt2[k], mod_p,
                   c, newprime, prod_p, dl->tmp, 1);
        c = MODCF(modgbs->modpolys[k], dl->coef3[k], modgbs->nprimes  - 1);
        mpz_CRT_ui(dl->crt3[k], dl->crt3[k], mod_p,
                   c, newprime, prod_p, dl->tmp, 1);
      }
  }
  mpz_set(mod_p, prod_p);
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
static inline int crt_ratrecon_lift_modgbs(gb_modpoly_t modgbs, data_lift_t dl,
                                           int32_t start, int32_t end,
                                           mpz_t mod_p, rrec_data_t rd1, rrec_data_t rd2,
                                           double *st_crt, double *st_rrec){

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
  if(maxlen > 0){
    dens = (mpz_t *)malloc(sizeof(mpz_t) * maxlen);
    for(uint32_t l = 0; l < maxlen; l++){
      mpz_init(dens[l]);
    }
  }

  /* Data needed by multi CRT functions */
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;
  fmpz_t y;
  if(maxlen > 0){
    fmpz_comb_init(comb, modgbs->primes, modgbs->nprimes);
    fmpz_comb_temp_init(comb_temp, comb);
    fmpz_init(y);
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
      int ok = 1;
      for(int32_t l = 0; l < polys[k]->len; l++){
        rt = realtime();
        for(uint32_t i = 0; i < modgbs->nprimes-1; i++){
          modgbs->cf_64[i] = MODCF(polys[k], l, i);
        }
        fmpz_multi_CRT_ui(y, modgbs->cf_64,
                          comb, comb_temp, 1);
        fmpz_get_mpz(zz, y);
        *st_crt += realtime() - rt;
        rt = realtime();
        if(lift_ratreconwden(polys[k]->cf_qq[l], rden, zz, mod_p, dl->tmp, rd1)){
          mpz_mul(dens[l], rden, dl->tmp);
          if(mpz_cmp_ui(rden, 1)){
            mpz_mul(dl->tmp, dl->tmp, rden);
            set_recdata_guess(rd1, mod_p, dl->tmp);
          }
          *st_rrec += realtime() - rt;
        }
        else{
          *st_rrec += realtime() - rt;
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
      dl->check2[k] += 1;
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
    fmpz_clear(y);
    fmpz_comb_temp_clear(comb_temp);
    fmpz_comb_clear(comb);
    for(uint32_t l = 0; l < maxlen; l++){
      mpz_clear(dens[l]);
    }
    free(dens);
  }
  mpz_clear(zz);
  mpz_clear(rden);

  return ret;
}


/* returns (coef == num / den mod prime) */
static inline int verif_coef(mpz_t num, mpz_t den, uint32_t prime, uint32_t coef){
  uint32_t lc = mpz_fdiv_ui(den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(num, prime);
  c *= lc;
  c = c % prime;

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
          /* all coefficients share the denominator lm */
          uint64_t ilc = mpz_fdiv_ui(pol[0]->lm, prime);
          ilc = mod_p_inverse_32(ilc, prime);
          for(int32_t c = 0; c < pol[0]->len; c++){
            uint64_t v = mpz_fdiv_ui(pol[0]->cf_qq[c], prime);
            v = (v * ilc) % prime;
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
        uint32_t coef = MODCF(modgbs->modpolys[k], dl->coef[k], modgbs->nprimes  - (thrds - i));
        int boo = verif_coef(dl->num[k], dl->den[k], prime, coef);
        if(boo){
          coef = MODCF(modgbs->modpolys[k], dl->coef2[k], modgbs->nprimes  - (thrds - i));
          boo = verif_coef(dl->num2[k], dl->den2[k], prime, coef);
        }
        if(boo){
          coef = MODCF(modgbs->modpolys[k], dl->coef3[k], modgbs->nprimes  - (thrds - i));
          boo = verif_coef(dl->num3[k], dl->den3[k], prime, coef);
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

static void ratrecon_gb(gb_modpoly_t modgbs, data_lift_t dl,
                        mpz_t mod_p, mpz_t prod_p,
                        rrec_data_t recdata1, rrec_data_t recdata2,
                        int thrds, double *st_crt, double *st_rrec, double *st_wit){

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
                                         st_crt, st_rrec);
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
  if(modgbs->nprimes % dl->rr == 0){
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
        int b = reconstructcoeff(dl, i, mod_p,
                               recdata1, recdata2);
        if(b){
          b = reconstructcoeff2(dl, i, mod_p, recdata1, recdata2);
        }
        if(b){
          b = reconstructcoeff3(dl, i, mod_p, recdata1, recdata2);
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
}

gb_modpoly_t *core_groebner_qq(
        gb_modpoly_t *modgbsp,
        bs_t *bs,
        mstrace_t msd,
        md_t *st,
        int32_t *errp,
        const len_t fc,
        const int print_gb
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

  uint32_t prime = 0;
  uint32_t primeinit = 0;
  uint32_t lprime = 1303905299;

  prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  while(fc == 0 && is_lucky_prime_ui(prime, bs)){
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  }

  primeinit = prime;
  msd->lp->p[0] = primeinit;
  if(fc){
    msd->lp->p[0] = fc;
    primeinit = fc;
  }
  if(info_level){
      fprintf(VERBSTREAM, "Initial prime = %d\n", msd->lp->p[0]);
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
                                                 msd->lp->p[0],
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
                                      lifted_nvars(st), &nb);
      data_lift_init(dlift, (*modgbsp)->ld, ldeg, nb);
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
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
    if(info_level){
        fprintf(VERBSTREAM, "New prime = %d\n", prime);
    }
    while(apply){

      prime = next_prime(prime);
      if(prime >= lprime){
        prime = next_prime(1<<30);
      }
      /* generate lucky prime numbers */
      msd->lp->p[0] = prime;
      while(is_lucky_prime_ui(prime, bs) || prime==primeinit){
        prime = next_prime(prime);
        if(prime >= lprime){
          prime = next_prime(1<<30);
        }
        msd->lp->p[0] = prime;
      }
      int nthrds = 1; /* mono-threaded mult-mid comp */
      for(len_t i = 1; i < nthrds/* st->nthrds */; i++){
        prime = next_prime(prime);
        if(prime >= lprime){
          prime = next_prime(1<<30);
        }
        msd->lp->p[i] = prime;
        while(is_lucky_prime_ui(prime, bs) || prime==primeinit){
          prime = next_prime(prime);
          if(prime >= lprime){
            prime = next_prime(1<<30);
          }
          msd->lp->p[i] = prime;
        }
      }
      prime = msd->lp->p[nthrds /* st->nthrds */ - 1];

      if((*modgbsp)->alloc <= nprimes + 2){
        /* Room for the images modulo the next primes. The multi-modular
         * computation is mono-threaded (one prime per step), so this must
         * not scale with st->nthrds. The increment is kept small (at most
         * ~1/16 of the current number of primes, but at least 4) since the
         * unused slots are allocated for every term of the basis which is
         * not lifted yet, i.e. for most of the basis at the end of the
         * computation, when memory consumption peaks. */
        gb_modpoly_realloc((*modgbsp), MAX(4, (uint32_t)nprimes / 16));
      }

      gb_modular_trace_application((*modgbsp), msd->mgb,
                                   msd->num_gb,
                                   msd->leadmons_ori,
                                   msd->leadmons_current,
                                   msd->btrace,
                                   bs,  st,
                                   0, /* info_level, */
                                   msd->lp,
                                   dlift->S, &stf4, msd->bad_primes);

      /* nprimes += st->nthrds; */
      nprimes += 1; /* at the moment, multi-mod comp is mono-threaded */
      if(nprimes == 1){
        if(info_level>2){
          fprintf(VERBSTREAM, "------------------------------------------\n");
          fprintf(VERBSTREAM, "#ADDITIONS       %13lu\n", (unsigned long)st->application_nr_add * 1000);
          fprintf(VERBSTREAM, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->application_nr_mult * 1000);
          fprintf(VERBSTREAM, "#REDUCTIONS      %13lu\n", (unsigned long)st->application_nr_red);
          fprintf(VERBSTREAM, "------------------------------------------\n");
        }
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
	    if (info_level > 1){
	      fprintf(VERBSTREAM,
		      "learning phase             %9.2f Gops/sec\n",
		      (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/(st->learning_rtime));
	      fprintf(VERBSTREAM,
		      "application phase          %9.2f Gops/sec\n",
		      (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4));
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
        /* return core_groebner_qq(modgbsp, bs, msd, st, errp, fc, print_gb); */
        goto restart;

      }

      int lstart = dlift->lstart;
      st_wit = 0;

      if(!bad){
        ratrecon_gb((*modgbsp), dlift, msd->mod_p, msd->prod_p, recdata1, recdata2,
                    nthrds/* st->nthrds */, &st_crt, &st_rrec, &st_wit);
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
  md->elim_full_basis = gens->elim_full_basis;

  initialize_mstrace(msd, md, bs);
  int err = 0;

  gb_modpoly_t *modgbsp = malloc(sizeof(gb_modpoly_t));
  modgbsp = core_groebner_qq(modgbsp, bs, msd, md, &err, field_char, print_gb);
  if (err) {
      fprintf(ERRSTREAM,"Problem with groebner_qq, stopped computation (%d).\n", err);
      exit(1);
  }

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
  if (gens->elim_full_basis && flags->elim_block_len > 0) {
    if (flags->elim_block_len == 1) {
      fprintf(ofile, "#monomial order:       eliminating first variable, blocks: graded reverse lexicographical\n");
    } else {
      fprintf(ofile, "#monomial order:       eliminating first %d variables, blocks: graded reverse lexicographical\n", flags->elim_block_len);
    }
  } else {
    fprintf(ofile, "#monomial order:       graded reverse lexicographical\n");
  }
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
      display_gbmodpoly_cf_qq(ofile, (*modgbsp), gens);
      fclose(ofile);
    }
    else{
      display_gbmodpoly_cf_qq(OUTSTREAM, (*modgbsp), gens);
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
                     print_gb);
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
