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

#define NBCHECK 2
#ifndef MIN
#define MIN(x, y) ((x) > (y) ? (y) : (x))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

typedef struct{
  uint32_t len; /* length of the encoded polynomial */
  uint32_t **cf_32; /* array of arrays of coefficients
                       * modulo several primes
                       */
  mpz_t *cf_zz; /* array which stores CRT lifting of
                    the coefficiels */
  mpz_t *cf_qq; /* array which stores rational coefficients
                  being lifted, numerators and denominators
                  are given given as mpz_t
               */
  mpz_t lm; /* stores the leading coefficient (hence all numerators of cf_qq should be divided by lm) */
} modpolys_struct;

typedef modpolys_struct modpolys_t[1];

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
  free(dl->check1);
  free(dl->check2);

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

  ht_t *ht = modgbs->bht;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  len_t idx;
  for(len_t i = 0; i < modgbs->ld; i++){
      idx = bs->lmps[i];
      if(i==idx){
          modgbs->hm[i] = malloc((1+lens[i] + OFFSET)*sizeof(hm_t));
          for(len_t j = 0; j<lens[i]+1; j++){
              modgbs->hm[i][j+OFFSET] = bs->hm[i][j+OFFSET];
          }
      }
      else{
          modgbs->hm[i] = NULL;
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
    modgbs->modpolys[i]->cf_32 = (uint32_t **)malloc(sizeof(uint32_t *)*lens[i]);
    modgbs->modpolys[i]->cf_zz = (mpz_t *)malloc(sizeof(mpz_t)*lens[i]);
    modgbs->modpolys[i]->cf_qq = (mpz_t *)malloc(sizeof(mpz_t)*2*lens[i]);
    for(uint32_t j = 0; j < lens[i]; j++){
      modgbs->modpolys[i]->cf_32[j] = calloc(alloc, sizeof(uint32_t));
      mpz_init(modgbs->modpolys[i]->cf_zz[j]);
    }
    for(uint32_t j = 0; j < 2 * lens[i]; j++){
      mpz_init(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_init(modgbs->modpolys[i]->lm);
    mpz_set_ui(modgbs->modpolys[i]->lm, 1);
  }
}


static inline void gb_modpoly_realloc(gb_modpoly_t modgbs,
                                      uint32_t newalloc,
                                      int32_t start){


  uint32_t oldalloc = modgbs->alloc;
  modgbs->alloc += newalloc;

  mp_limb_t *newprimes = (mp_limb_t *)realloc(modgbs->primes,
                                            modgbs->alloc * sizeof(mp_limb_t));

  if(newprimes == NULL){
    fprintf(stderr, "Problem when reallocating modgbs (primes)\n");
    exit(1);
  }
  modgbs->primes = newprimes;
  for(uint32_t i = oldalloc; i < modgbs->alloc; i++){
    modgbs->primes[i] = 0;
  }

  mp_limb_t *ncf_64 = (mp_limb_t *)realloc(modgbs->cf_64,
                                         modgbs->alloc * sizeof(mp_limb_t));
  if(ncf_64 == NULL){
    fprintf(stderr, "Problem when reallocating modgbs (cfs)\n");
    exit(1);
  }
  modgbs->cf_64 = ncf_64;
  for(uint32_t i = oldalloc; i < modgbs->alloc; i++){
    modgbs->cf_64[i] = 0;
  }
  for(uint32_t i = start; i < modgbs->ld; i++){
    for(uint32_t j = 0; j < modgbs->modpolys[i]->len; j++){
      uint32_t *newcfs_pol = (uint32_t *)realloc(modgbs->modpolys[i]->cf_32[j],
                                                 modgbs->alloc * sizeof(uint32_t));
      if(newcfs_pol == NULL){
        fprintf(stderr, "Problem when reallocating modgbs (cfs_pol)\n");
      }
      modgbs->modpolys[i]->cf_32[j] = newcfs_pol;
      for(uint32_t k = oldalloc; k < modgbs->alloc; k++){
        modgbs->modpolys[i]->cf_32[j][k] = 0;
      }
    }
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
    for(uint32_t j = 0; j < len; j++){
      fprintf(file, "[");
      for(uint32_t k = 0; k < modgbs->alloc-1; k++){
        fprintf(file, "%d, ", modgbs->modpolys[i]->cf_32[j][k]);
      }
      if(j < len - 1){
        fprintf(file, "%d], ", modgbs->modpolys[i]->cf_32[j][modgbs->alloc-1]);
      }
      else{
        fprintf(file, "%d]\n", modgbs->modpolys[i]->cf_32[j][modgbs->alloc-1]);
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

  len_t idx, i, j, k;
  hm_t *hm    = NULL;
  idx = modgbs->lmps[pos];
  hm  = modgbs->hm[idx]+OFFSET;

  ht_t *ht = modgbs->bht;
  const len_t nv  = ht->nv;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  int *evi    =   (int *)malloc((unsigned long)ht->nv * sizeof(int));
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
  for(i = modgbs->modpolys[pos]->len-1; i > 0 ; i--){
    if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0)>0){
      fprintf(file, "+");
    }
    //coef is -1
    if(mpz_cmp_si(modgbs->modpolys[pos]->cf_qq[2*i], -1) == 0 && mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) == 0){
        fprintf(file, "-");
    }
    int star = 0;
    //coef is neither 1 nor -1
    if((mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 1) != 0 || mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) != 0) && (mpz_cmp_si(modgbs->modpolys[pos]->cf_qq[2*i], -1) != 0 || mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) != 0)){
        mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[2*i]);
        if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) != 0){
            fprintf(file, "/");
            mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[2*i + 1]);
        }
        fprintf(file, "*");
    }
    if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0) != 0){
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
  if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[0], 0) > 0){
    fprintf(file, "+");
  }
  int deg = 0;
  for (k = 0; k < nv; ++k) {
     deg += ht->ev[hm[len-1]][evi[k]];
  }
  if(deg == 0){
    mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[0]);
    if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) != 0){
      fprintf(file, "/");
      mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[1]);
    }
  }
  else{
    if(mpz_cmp_si(modgbs->modpolys[pos]->cf_qq[0], -1) == 0 && mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) == 0){
        fprintf(file, "-");
    }
    int star = 0;
    //coef is neither 1 nor -1
    if((mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[0], 1) != 0 || mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) != 0) && (mpz_cmp_si(modgbs->modpolys[pos]->cf_qq[0], -1) != 0 || mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) != 0)){
        mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[0]);
        if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) != 0){
            fprintf(file, "/");
            mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[1]);
        }
        fprintf(file, "*");
    }
    fflush(file);
    if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[0], 0) != 0){
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
  fprintf(file, "\n");
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
  fprintf(file, "\n");
  fprintf(file, "]:\n");
}

static inline void gb_modpoly_without_hash_table_clear(gb_modpoly_t modgbs){
  free(modgbs->primes);
  free(modgbs->ldm);
  for(uint32_t i = 0; i < modgbs->ld; i++){
    for(uint32_t j = 0; j < modgbs->modpolys[i]->len; j++){
      free(modgbs->modpolys[i]->cf_32[j]);
      mpz_clear(modgbs->modpolys[i]->cf_zz[j]);
    }
    for(uint32_t j = 0; j < 2 * modgbs->modpolys[i]->len; j++){
      mpz_clear(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_clear(modgbs->modpolys[i]->lm);
    free(modgbs->modpolys[i]->cf_32);
    free(modgbs->modpolys[i]->cf_zz);
    free(modgbs->modpolys[i]->cf_qq);
  }
  free(modgbs->lmps);
  free(modgbs->hm);
  free(modgbs->modpolys);
}

static inline void gb_modpoly_clear(gb_modpoly_t modgbs){
  free(modgbs->primes);
  free(modgbs->ldm);
  for(uint32_t i = 0; i < modgbs->ld; i++){
    for(uint32_t j = 0; j < modgbs->modpolys[i]->len; j++){
      free(modgbs->modpolys[i]->cf_32[j]);
      mpz_clear(modgbs->modpolys[i]->cf_zz[j]);
    }
    for(uint32_t j = 0; j < 2 * modgbs->modpolys[i]->len; j++){
      mpz_clear(modgbs->modpolys[i]->cf_qq[j]);
    }
    mpz_clear(modgbs->modpolys[i]->lm);
    free(modgbs->modpolys[i]->cf_32);
    free(modgbs->modpolys[i]->cf_zz);
    free(modgbs->modpolys[i]->cf_qq);
  }
  free(modgbs->lmps);
  free(modgbs->hm);
  //free_hash_table(&(modgbs->bht));
  free_shared_hash_data(modgbs->bht);
  free(modgbs->modpolys);
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
  int32_t *ldeg = calloc(sizeof(int32_t), *nb);
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
    fprintf(stderr, "Not enough space in modgbs\n");
    exit(1);
  }
  modgbs->primes[modgbs->nprimes] = fc;

  len_t i, j, k, idx;

  /*****************************************************
   * It should be checked if modbgs->ht needs to be updated as well as
   * modgbs->lmps
   * *****************************************************/

  len_t len   = 0;
  hm_t *hm    = NULL;

  const len_t nv  = ht->nv;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  for(i = start; i < modgbs->ld; i++){
    idx = bs->lmps[i];
    if (bs->hm[idx] == NULL) {
      fprintf(stderr, " poly is 0\n");
      exit(1);
    } else {
      hm  = bs->hm[idx]+OFFSET;
      len = bs->hm[idx][LENGTH];
    }
    int32_t bc = modgbs->modpolys[i]->len - 1;
    for (j = 1; j < len; ++j) {
      uint32_t c = bs->cf_32[bs->hm[idx][COEFFS]][j];
      modgbs->modpolys[i]->cf_32[bc][modgbs->nprimes] = c;
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
      printf("Problem with F4, stopped computation.\n");
      exit(1);
    }

    /* rt = realtime()-ca0; */
    st->learning_rtime = realtime()-ca0;

    ht_t *bht = bs->ht;

    if(info_level > 2){
        fprintf(stdout, "------------------------------------------\n");
        fprintf(stdout, "#ADDITIONS       %13lu\n", (unsigned long)st->trace_nr_add * 1000);
        fprintf(stdout, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->trace_nr_mult * 1000);
        fprintf(stdout, "#REDUCTIONS      %13lu\n", (unsigned long)st->trace_nr_red);
        fprintf(stdout, "------------------------------------------\n");
    }

    /* Leading monomials from Grobner basis */
    int32_t *bexp_lm = get_lm_from_bs(bs, bht);
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
    /************************************************/

    int32_t *lmb;
    if(st->nev){
      lmb = NULL; //monomial_basis_enlarged(num_gb[0], bht->nv - st->nev,  bexp_lm2, &dquot);
    }
    else{
      lmb = NULL; //monomial_basis_enlarged(num_gb[0], bht->nv,  bexp_lm, &dquot);
    }

    /************************************************/
    /************************************************/

    int32_t *lens = array_of_lengths(leadmons[0], num_gb[0], bs, bht->nv - st->nev);

    if(truncate_lifting != 0 && truncate_lifting < num_gb[0]){
      gb_modpoly_init(modgbs, 2, lens, bs, bht->nv - st->nev, truncate_lifting, leadmons[0], st);
    }
    else{
      gb_modpoly_init(modgbs, 2, lens, bs, bht->nv - st->nev, num_gb[0], leadmons[0], st);
    }
    modpgbs_set(modgbs, bs, bht, fc, start, st->nev);
    int is_empty = 0;
    if(bs->lml == 1){
        if(info_level){
            fprintf(stdout, "Grobner basis has a single element\n");
        }
        is_empty = 1;
        for(int i = 0; i < bht->nv; i++){
            if(bexp_lm[i]!=0){
                is_empty = 0;
            }
        }
        if(is_empty){
            free_basis_without_hash_table(&(bs));
            return is_empty;
        }
    }

    /**************************************************/

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
  hl_t eld = bs_qq->ht->eld;
  bs_t *bs = NULL;
  int32_t error = 0;
  bs = core_gba(bs_qq, st, &error, lp->p[0]);
  *stf4 = realtime()-rt;
  ht_t **bht = &(bs->ht);
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
  if (bs == NULL) {
      bad_primes[0] = 1;
      return;
  }
  int32_t lml = bs->lml;
  if (st->nev > 0) {
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
        free_basis_and_only_local_hash_table_data(&bs);
      }
      return;
  }

  if(st->nev){
    get_lm_from_bs_trace_elim(bs, bht[0], leadmons_current[0], num_gb[0]);
  }
  else{
    get_lm_from_bs_trace(bs, bht[0], leadmons_current[0]);
  }

  if(!equal_staircase(leadmons_current[0], leadmons_ori[0],
                      num_gb[0], num_gb[0], bht[0]->nv - st->nev)){
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

static inline void choose_coef_to_lift(gb_modpoly_t modgbs, data_lift_t dlift){
  uint32_t ld = modgbs->ld;
  for(int32_t i = 0; i < ld; i++){
    len_t d = 0;
    len_t len = modgbs->modpolys[i]->len;
    while(d < len - 1 && len > 0){
      if(modgbs->modpolys[i]->cf_32[d][0]){
        dlift->coef[i] = d;
        break;
      }
      else{
        d++;
      }
    }
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
        uint32_t c = modgbs->modpolys[k]->cf_32[coef[k]][modgbs->nprimes  - 1 ];
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
  for(int32_t k = 0; k < dl->end; k++){
      if(modgbs->modpolys[k]->len){
        uint64_t c = modgbs->modpolys[k]->cf_32[coef[k]][modgbs->nprimes  - 1 ];
        mpz_CRT_ui(dl->crt[k], dl->crt[k], mod_p,
                   c, newprime, prod_p, dl->tmp, 1);
      }
  }
  mpz_set(mod_p, prod_p);
}


static inline void crt_lift_modgbs(gb_modpoly_t modgbs, data_lift_t dlift,
                                   int32_t start, int32_t end){
  /* Data needed by multi CRT functions */
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;

  fmpz_comb_init(comb, modgbs->primes, modgbs->nprimes);
  fmpz_comb_temp_init(comb_temp, comb);
  fmpz_t y;
  fmpz_init(y);

  modpolys_t *polys = modgbs->modpolys;

  for(int32_t k = start; k < end; k++){
    if(dlift->check1[k] >= 1){
      for(int32_t l = 0; l < polys[k]->len; l++){
        for(uint32_t i = 0; i < modgbs->nprimes-1; i++){
          modgbs->cf_64[i] = polys[k]->cf_32[l][i];
        }
        fmpz_multi_CRT_ui(y, modgbs->cf_64,
                          comb, comb_temp, 1);
        fmpz_get_mpz(polys[k]->cf_zz[l], y);
      }
    }
  }

  fmpz_clear(y);
  fmpz_comb_temp_clear(comb_temp);
  fmpz_comb_clear(comb);

}


static inline void set_recdata(data_lift_t dl, rrec_data_t rd1, rrec_data_t rd2, mpz_t mod_p){
  mpz_fdiv_q_2exp(rd1->N, mod_p, 1);

  if(dl->lstart){

    mpz_set(rd2->N, rd1->N);

    mpz_sqrt(rd2->N, rd2->N);
    mpz_set(rd2->D, rd2->N);

    mpz_root(rd1->D, rd1->N, 3);
    mpz_fdiv_q(rd1->N, rd1->N, rd1->D);

  }
  else{
    mpz_sqrt(rd1->N, rd1->N);
    mpz_set(rd1->D, rd1->N);

    mpz_set(rd2->N, rd1->N);
    mpz_set(rd2->D, rd1->D);

  }
}

/*
  returns the index of the poly if some of its coeff could not be lifted
  else returns -1
 */
static inline int ratrecon_lift_modgbs(gb_modpoly_t modgbs, data_lift_t dl,
                                       int32_t start, int32_t end,
                                       mpz_t mod_p, rrec_data_t rd1, rrec_data_t rd2){
 
  mpz_t rnum, rden;
  mpz_init(rnum);
  mpz_init(rden);

  modpolys_t *polys = modgbs->modpolys;
  set_recdata(dl, rd1, rd2, mod_p);

  for(int32_t k = start; k < end; k++){
    if(dl->check1[k] > 1 && polys[k]->len > 0){
      mpz_set_ui(dl->tmp, 1);
      for(int32_t l = 0; l < polys[k]->len; l++){
        if(ratreconwden(rnum, rden, polys[k]->cf_zz[l], mod_p, dl->den[k], rd1)){
          mpz_set(polys[k]->cf_qq[2*l], rnum);
          mpz_set(polys[k]->cf_qq[2*l + 1], rden);
          mpz_lcm(dl->tmp, dl->tmp, rden);
        }
        else{
          mpz_set_ui(dl->gden, 1);
          dl->check2[k] = 0;
          mpz_clear(rnum);
          mpz_clear(rden);
          return k;
        }
      }
      mpz_set(polys[k]->lm, dl->den[k]);
      for(int32_t l = 0; l < polys[k]->len; l++){
        mpz_mul(polys[k]->cf_qq[2*l], polys[k]->cf_qq[2*l], dl->tmp);
        mpz_divexact(polys[k]->cf_qq[2*l], polys[k]->cf_qq[2*l], polys[k]->cf_qq[2*l+1]);
        mpz_set_ui(polys[k]->cf_qq[2*l+1], 1);
      }
      mpz_mul(polys[k]->lm, polys[k]->lm, dl->tmp);
      dl->check2[k] = 1;
    }
    if(dl->check1[k] < 1 && polys[k]->len > 0){
      mpz_clear(rnum);
      mpz_clear(rden);
      return k;
    }
    if(polys[k]->len == 0){
        dl->check2[k] = 1;
    }
  }
  mpz_clear(rnum);
  mpz_clear(rden);
  return end;
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
  /* verification of the basis is performed at the very end of the computation */
  if(dl->check1[dl->end-1] == 0){
      return 1;
  }
  int b = 1;
  mpz_t den;
  mpz_init(den);
  for(int32_t k = 0; k < modgbs->ld; k++){
    if(dl->check1[k]>=1 && dl->check2[k]>0 && dl->check2[k] < NBCHECK ){
      if(modgbs->modpolys[k]->len == 0){
        dl->check2[k]++;
      }
      else{
        for(int i = 0; i < thrds; i++){
          mp_limb_t prime = modgbs->primes[modgbs->nprimes - (thrds - i) ];
          for(int32_t c = 0; c < modgbs->modpolys[k]->len; c++){
            mpz_mul(den, modgbs->modpolys[k]->lm, modgbs->modpolys[k]->cf_qq[2*c+1]);
            uint32_t coef = modgbs->modpolys[k]->cf_32[c][modgbs->nprimes   - (thrds - i) ];
            int b = verif_coef(modgbs->modpolys[k]->cf_qq[2*c], den, prime, coef);
            if(!b){
              dl->check2[k] = 0;
              mpz_set_ui(dl->gden, 1);
              b = 0;
              mpz_clear(den);
              return b;
            }
          }
        dl->check2[k]++;
        }
      }
    }
  }
  mpz_clear(den);
  return b;
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
        uint32_t coef = modgbs->modpolys[k]->cf_32[dl->coef[k]][modgbs->nprimes  - (thrds - i) ];
        int boo = verif_coef(dl->num[k], dl->den[k], prime, coef);

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

  int b = ratreconwden(dl->num[i], dl->den[i],
                       dl->crt[i], mod_p, dl->gden, recdata1);
  if(b){
    mpz_mul(dl->den[i], dl->den[i], dl->gden);
    mpz_gcd(dl->tmp, dl->den[i], dl->num[i]);

    mpz_divexact(dl->num[i], dl->num[i], dl->tmp);
    mpz_divexact(dl->den[i], dl->den[i], dl->tmp);

    mpz_set(dl->gden, dl->den[i]);
  }
  else{

    b = ratrecon(dl->num[i], dl->den[i],
                 dl->crt[i], mod_p, recdata2);
    if(b){
      mpz_set(dl->gden, dl->den[i]);
    }
  }
  return b;
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
  *st_rrec += realtime()-st;
  if(dl->lstart != dl->start){
    crt_lift_modgbs(modgbs, dl, dl->lstart, dl->start);
    *st_crt += realtime() - st;
    st = realtime();
    dl->start = ratrecon_lift_modgbs(modgbs, dl, dl->lstart, dl->start, mod_p, recdata1, recdata2);
    *st_rrec += realtime()-st;
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

  set_recdata(dl, recdata1, recdata2, mod_p);

  st = realtime();
  if(modgbs->nprimes % dl->rr == 0){
    for(int32_t i = dl->lstart; i < dl->lend; i++){
      if(modgbs->modpolys[i]->len==0){
          dl->check1[i]++;
      }
      else{
        int b = reconstructcoeff(dl, i, mod_p,
                               recdata1, recdata2);
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
      nb = MAX(nb, mpz_sizeinbase(modgbs->modpolys[i]->cf_qq[2*j], 2));
      nb = MAX(nb, mpz_sizeinbase(modgbs->modpolys[i]->cf_qq[2*j + 1], 2));
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
  double st0 = realtime();

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
  srand(time(0));

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
      fprintf(stdout, "Initial prime = %d\n", msd->lp->p[0]);
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

    int32_t is_empty = gb_modular_trace_learning(*modgbsp,
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
    if(is_empty == 1 || print_gb == 1){
      if(dlinit){
        data_lift_clear(dlift);
      }
      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      return modgbsp;
    }

    apply = 1;

    gb_modpoly_realloc((*modgbsp), 1, dlift->S);

#ifdef DEBUGGBLIFT
    display_gbmodpoly_cf_32(stderr, (*modgbsp));
#endif

    if(!dlinit){
      int nb = 0;
      int32_t *ldeg = array_nbdegrees((*msd->leadmons_ori), msd->num_gb[0],
                                      bs->ht->nv - st->nev, &nb);
      data_lift_init(dlift, (*modgbsp)->ld, ldeg, nb);
      choose_coef_to_lift((*modgbsp), dlift);
      free(ldeg);
      dlinit = 1;
    }

    if(info_level){
      fprintf(stdout,"\n---------- COMPUTATIONAL DATA -----------\n");
      int s= 0;
      for(int i = 0; i < dlift->nsteps; i++){
        fprintf(stdout, "[%d]", dlift->steps[i]);
	fflush(stdout);
        s+=dlift->steps[i];
      }
      fprintf(stdout, "\n");
      fprintf(stdout,
	      "#polynomials to lift %14lu\n",
	      (unsigned long) s);
      fprintf(stdout, "-----------------------------------------\n");
    }

    if(is_empty || success == 0 || fc) {

      apply = 0;
      if(dlinit){
        data_lift_clear(dlift);
      }

      gb_modpoly_without_hash_table_clear((*modgbsp));

      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      fprintf(stdout, "Something went wrong in the learning phase, msolve restarts.");
      return core_groebner_qq(modgbsp, bs, msd, st, errp, fc, print_gb); 
    }
    /* duplicate data for multi-threaded multi-mod computation */
    duplicate_data_mthread_gbtrace(st->nthrds, bs, st, msd->num_gb,
                                   msd->leadmons_ori, msd->leadmons_current,
                                   msd->btrace);


    learn = 0;
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
    if(info_level){
        fprintf(stdout, "New prime = %d\n", prime);
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
        gb_modpoly_realloc((*modgbsp), 16*st->nthrds, dlift->S);
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
          fprintf(stdout, "------------------------------------------\n");
          fprintf(stdout, "#ADDITIONS       %13lu\n", (unsigned long)st->application_nr_add * 1000);
          fprintf(stdout, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->application_nr_mult * 1000);
          fprintf(stdout, "#REDUCTIONS      %13lu\n", (unsigned long)st->application_nr_red);
          fprintf(stdout, "------------------------------------------\n");
        }
        /* if(info_level>1){ */
        /*   fprintf(stderr, "Application phase %.2f Gops/sec\n", */
        /*           (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4)); */
        /*   fprintf(stderr, "Elapsed time: %.2f\n", stf4); */
        /* } */
	if(info_level){
	  fprintf(stdout,
		  "\n---------------- TIMINGS ----------------\n");
	  fprintf(stdout,
		  "multi-mod overall(elapsed) %9.2f sec\n",
		  stf4);
	  if (info_level > 1){
	    fprintf(stdout,
		    "learning phase             %9.2f Gops/sec\n",
		    (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/(st->learning_rtime));
	    fprintf(stdout,
		    "application phase          %9.2f Gops/sec\n",
		    (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4));
	  }
	  fprintf(stdout,
		  "-----------------------------------------\n");
      }
      if (info_level) {
	  fprintf(stdout,
		  "\nmulti-modular steps\n");
	  fprintf(stdout, "-------------------------------------------------\
-----------------------------------------------------\n");
      }

      }
      int bad = 0;
      for(int i = 0; i < nthrds/* st->nthrds */; i++){
        if(msd->bad_primes[i] == 1){
          bad = 1;
          if(info_level > 1){
              fprintf(stdout, "[!]");
          }
          nbadprimes++;
          msd->bad_primes[i] = 0;
        }
        if(msd->bad_primes[i] == 2){
            bad = 1;
            if(info_level > 1){
                fprintf(stdout, "[!!]");
            }
            nbadprimes = nprimes;
            msd->bad_primes[i] = 0;
        }
      }

      if(nbadprimes >= nprimes){
        if(info_level){
          fprintf(stdout, "Too many bad primes, computation will restart\n");
        }
        if(dlinit){
          data_lift_clear(dlift);
        }

        gb_modpoly_without_hash_table_clear((*modgbsp));
        free_rrec_data(recdata1);
        free_rrec_data(recdata2);
        return core_groebner_qq(modgbsp, bs, msd, st, errp, fc, print_gb); 

      }

      int lstart = dlift->lstart;
      double ost_rrec = st_rrec;
      double ost_crt = st_crt;
      st_wit = 0;

      if(!bad){
        ratrecon_gb((*modgbsp), dlift, msd->mod_p, msd->prod_p, recdata1, recdata2,
                    nthrds/* st->nthrds */, &st_crt, &st_rrec, &st_wit);
      }
      if(st_wit > 2*dlift->rr * stf4){
        dlift->rr = 2*dlift->rr;
        if(info_level){
          fprintf(stdout, "(->%d)", dlift->rr);
	      fflush(stdout);
        }
      }
      if(info_level){
        if(!(nprimes & (nprimes - 1))){
          fprintf(stdout, "{%d}", nprimes);
	      fflush(stdout);
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
          fprintf(stdout, "<%.2f%%>", 100* (float)MIN((dlift->lstart + 1), (*modgbsp)->ld)/(*modgbsp)->ld);
	  fflush(stdout);
        }
        lstart = dlift->lstart;
      }
      /* this is where learn could be reset to 1 */
      /* but then duplicated datas and others should be free-ed */
    }
    if (info_level){
      fprintf(stdout, " \n-------------------------------------------------\
-----------------------------------------------------\n");
    }
  }
  /* if(info_level){ */
  /*   fprintf(stderr, "\nCRT time = %.2f, Rational reconstruction time = %.2f\n", st_crt, st_rrec); */
  /* } */
  if(info_level){
    long nbits = max_bit_size_gb((*modgbsp));
    /* fprintf(stderr, "Maximum bit size of the coefficients: %ld\n", nbits); */
    /* fprintf(stderr, "%d primes used. \nElapsed time: %.2f\n", nprimes, realtime()-st0); */
    fprintf(stdout,"\n\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(stdout, "Max coeff. bitsize %16lu\n", (unsigned long) nbits);
    fprintf(stdout, "#primes            %16lu\n", (unsigned long) nprimes);
    fprintf(stdout, "#bad primes        %16lu\n", (unsigned long) nbadprimes);
    fprintf(stdout, "-----------------------------------------\n");
    fprintf(stdout, "\n---------------- TIMINGS ----------------\n");
    fprintf(stdout, "CRT     (elapsed)         %10.2f sec\n", st_crt);
    fprintf(stdout, "ratrecon(elapsed)         %10.2f sec\n", st_rrec);
    /* fprintf(stdout, "CRT and ratrecon(elapsed) %10.2f sec\n", realtime()-st0); */
    fprintf(stdout, "-----------------------------------------\n");
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
        gb_modpoly_t gb
        )
{
    int64_t nelts = gb->ld;
    /* Over QQ, if eliminating, these variables are no longer included in */
    /* gb->nv, thus we need to help us with adding elim_block_len */
    /* correspondingly when exporting the basis. */
    int32_t nv    = gb->nv;
    int32_t nve   = gb->nv + elim_block_len;

    *bld  = nelts;

    int32_t *lens = (int32_t *)(*mallocp)(((uint64_t)nelts) * sizeof(int32_t));

    int64_t nterms = 0;

    for(int32_t i = 0; i < nelts; i++){
        int32_t nlen = 0;
        int32_t len = gb->modpolys[i]->len;

        for(int32_t j = len-1; j >= 0; j--){
            if(mpz_cmp_ui(gb->modpolys[i]->cf_qq[2*j], 0) != 0){
                nlen++;
            }
        }
        nlen++; /* to take into account the lm */
        lens[i] = nlen;
        nterms += nlen;
    }

    int32_t *exp = (int32_t *)(*mallocp)( 
            ((uint64_t)nterms) * ((uint64_t) nve) * sizeof(int32_t));

    memset(exp, 0, (unsigned long)nterms * nve * sizeof(int32_t));

    mpz_t *cf_qq = (mpz_t *)malloc(
            nterms *sizeof(mpz_t));
    for(int64_t i = 0; i < nterms; i++){
        mpz_init(cf_qq[i]);
    }

    hm_t *hm    = NULL;
    ht_t *ht = gb->bht;
    const len_t ebl = ht->ebl;
    const len_t evl = ht->evl;
    int *evi    =   (int *)malloc((unsigned long)ht->nv * sizeof(int));
    if (ebl == 0) {
      for (len_t i = 1; i < evl; ++i) {
        evi[i-1]    =   i;
      }
    } else {
      for (len_t i = 1; i < ebl; ++i) {
        evi[i-1]    =   i;
      }
      for (len_t i = ebl+1; i < evl; ++i) {
        evi[i-2]    =   i;
      }
    }

    int64_t term = 0;
    for(int64_t p = 0; p < nelts; p++){

        len_t idx = gb->lmps[p];
        hm  = gb->hm[idx]+OFFSET;
        int32_t l = gb->modpolys[p]->len;
        for(int32_t n = 0; n < nv; n++){
            exp[term * nve + n + elim_block_len] = gb->ldm[p * nv + n];
        }
        mpz_set(cf_qq[term], gb->modpolys[p]->lm);

        term++;
        for(int32_t i = l-1; i >= 0; i--){
            
            if(mpz_cmp_ui(gb->modpolys[p]->cf_qq[2*i], 0) != 0){
                for(int32_t n = 0 ; n < nv; n++){
                    exp[term * nve + n + elim_block_len] = ht->ev[hm[l-i]][evi[n]];
                }
                mpz_set(cf_qq[term], gb->modpolys[p]->cf_qq[2*i]);

                term++;
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
  int32_t *dim_ptr = &flags->dim;
  int64_t *dquot_ptr = &flags->dquot;
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
  int32_t nr_nf = 0;
  const uint32_t prime_start = pow(2, 30);

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
      printf("Bad input data, stopped computation.\n");
      exit(1);
  }

  initialize_mstrace(msd, md, bs);
  int err = 0;

  gb_modpoly_t *modgbsp = malloc(sizeof(gb_modpoly_t));
  modgbsp = core_groebner_qq(modgbsp, bs, msd, md, &err, field_char, print_gb);
  if (err) {
      printf("Problem with groebner_qq, stopped computation (%d).\n", err);
      exit(1);
  }

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->f4_ctime = ct1 - ct0;
  md->f4_rtime = rt1 - rt0;

  get_and_print_final_statistics(stderr, md, bs);

  /* free and clean up */
  free_mstrace(msd, md);

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
    ofile = stdout;
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
      display_gbmodpoly_cf_qq(ofile, (*modgbsp), gens);
      fclose(ofile);
    }
    else{
      display_gbmodpoly_cf_qq(stdout, (*modgbsp), gens);
    }
  }
  if(flags->print_gb == 1){
    if(flags->files->out_file != NULL){
      FILE *ofile = fopen(flags->files->out_file, "ab+");
      display_lm_gbmodpoly_cf_qq(ofile, (*modgbsp), gens);
      fclose(ofile);
    }
    else{
      display_lm_gbmodpoly_cf_qq(stdout, (*modgbsp), gens);
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
        printf("Bad input data, stopped computation.\n");
        exit(1);
    }

    initialize_mstrace(msd, md, bs);
    int err = 0;

    gb_modpoly_t *modgbsp = malloc(sizeof(gb_modpoly_t));
    modgbsp = core_groebner_qq(modgbsp, bs, msd, md, &err, field_char, 
            2/* if set to 1, only the LM of the Gbs are correct */);
    if (err) {
        printf("Problem with groebner_qq, stopped computation.\n");
        exit(1);
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    md->f4_ctime = ct1 - ct0;
    md->f4_rtime = rt1 - rt0;

    get_and_print_final_statistics(stderr, md, bs);

    /* free and clean up */
    free_mstrace(msd, md);
    free(md);
    md    = NULL;

    int64_t nterms  = export_results_from_groebner_qq(bld, blen, bexp,
            bcf, mallocp, elim_block_len, (*modgbsp));

    gb_modpoly_clear((*modgbsp));
    
    return nterms;

}
