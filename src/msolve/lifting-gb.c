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

#define MAX(a,b) (((a)>(b))?(a):(b))

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
  uint64_t *primes; /* array of prime numbers encoded with uint64_t to ensure
                       compatibility with flint */
  uint64_t *cf_64; /* array of length equal to number of primes which will be used
                    to copy coefficients (hence ensuring compatibility with
                    flint) */
  uint32_t ld; /* number of polynomials */
  int nv; /* number of variables */
  int32_t *mb; /* monomial basis enlarged */
  int32_t *ldm; /* lead monomials */
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
  int32_t start; /* indicates smallest index of poly whose coef has been lifted but not checked*/
  int32_t end; /* indicates largest index of poly whose coef has been lifted but not checked*/
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

  dl->steps = calloc(nsteps, sizeof(int32_t));
  for(i = 0; i < nsteps; i++){
    dl->steps[i] = steps[i];
  }
  dl->cstep = 0;
  dl->lend = npol;
  dl->crt_mult = 0;
  dl->crt = malloc(sizeof(mpz_t) * dl->npol);
  for(int32_t i = 0; i < dl->npol; i++){
    mpz_init(dl->crt[i]);
  }
  dl->recon = 0;
  dl->coef = calloc(npol, sizeof(mpz_t) );

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
  dl->end = 0;
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
                                   int nv, uint32_t ld,
                                   int32_t *lm, int32_t *basis){
  modgbs->alloc = alloc;
  modgbs->nprimes = 0;
  modgbs->primes = calloc(alloc, sizeof(uint64_t));
  modgbs->cf_64 = calloc(alloc, sizeof(uint64_t));
  modgbs->ld = ld;
  modgbs->nv = nv;
  modgbs->modpolys = malloc(sizeof(modpolys_struct) * ld);

  modgbs->mb = basis;
  modgbs->ldm = calloc(nv*ld, sizeof(int32_t));
  for(int32_t i = 0; i < ld; i++){
    for(int j = 0; j < nv; j++){
      modgbs->ldm[i*nv+j] = lm[i*nv+j];
    }
  }
  for(uint32_t i = 0; i < ld; i++){
    modgbs->modpolys[i]->len = lens[i];
    modgbs->modpolys[i]->cf_32 = malloc(sizeof(uint32_t **)*lens[i]);
    modgbs->modpolys[i]->cf_zz = malloc(sizeof(mpz_t)*lens[i]);
    modgbs->modpolys[i]->cf_qq = malloc(sizeof(mpz_t)*2*lens[i]);
    for(uint32_t j = 0; j < lens[i]; j++){
      modgbs->modpolys[i]->cf_32[j] = calloc(sizeof(uint32_t), alloc);
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

  uint64_t *newprimes = (uint64_t *)realloc(modgbs->primes,
                                            modgbs->alloc * sizeof(uint64_t));

  if(newprimes == NULL){
    fprintf(stderr, "Problem when reallocating modgbs (primes)\n");
    exit(1);
  }
  modgbs->primes = newprimes;
  for(uint32_t i = oldalloc; i < modgbs->alloc; i++){
    modgbs->primes[i] = 0;
  }

  uint64_t *ncf_64 = (uint64_t *)realloc(modgbs->cf_64,
                                         modgbs->alloc * sizeof(uint64_t));
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
  fprintf(stderr, "primes = [");
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
      fprintf(stderr, "[");
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
    display_monomial(file, gens, pos, &modgbs->ldm);
    return;
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->lm, 1) != 0){
    mpz_out_str(file, 10, modgbs->modpolys[pos]->lm);
    fprintf(file, "*");
  }
  display_monomial_single(file, gens, pos, &modgbs->ldm);
  for(int32_t i = modgbs->modpolys[pos]->len-1; i > 0 ; i--){
    if((mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0) != 0 && mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 1) != 0) || mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) != 0){
      if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0)>0){
        fprintf(file, "+");
      }
      mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[2*i]);
      if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i + 1], 1) != 0){
        fprintf(file, "/");
        mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[2*i + 1]);
      }
      fprintf(file, "*");
    }
    else{
      if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0) != 0){
        fprintf(file, "+");
      }
    }
    if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[2*i], 0) != 0){
      display_monomial_single(file, gens, i, &modgbs->mb);
    }
    fflush(file);
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[0], 0) > 0){
    fprintf(file, "+");
    mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[0]);
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[0], 0) < 0){
    mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[0]);
  }
  if(mpz_cmp_ui(modgbs->modpolys[pos]->cf_qq[1], 1) != 0){
    fprintf(file, "/");
    mpz_out_str(file, 10, modgbs->modpolys[pos]->cf_qq[1]);
  }

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
    if(modgbs->modpolys[i]->len == 0){
      display_monomial(file, gens, i, &modgbs->ldm);
    }
    else{
      display_monomial_single(file, gens, i, &modgbs->ldm);
    }
    fprintf(file, ", \n");
  }
  if(modgbs->modpolys[p-1]->len == 0){
    display_monomial(file, gens, p-1, &modgbs->ldm);
  }
  else{
    display_monomial_single(file, gens, p-1, &modgbs->ldm);
  }
  fprintf(file, "\n");
  fprintf(file, "]:\n");
}

static inline void gb_modpoly_clear(gb_modpoly_t modgbs){
  free(modgbs->primes);
  free(modgbs->mb);
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
static inline int32_t *array_of_lengths(int32_t *bexp_lm, int len,
                                        int32_t *basis, int dquot, int nv){
  int32_t *lens = calloc(sizeof(int32_t), len);
  for(int i = 0; i < len; i++){
    lens[i] = compute_length(bexp_lm + (i * nv), nv, basis, dquot);
  }
  return lens;
}

/* returns 0 in case of failure else returns 1 */
static inline int modpgbs_set(gb_modpoly_t modgbs,
                              const bs_t *bs, const ht_t * const ht,
                              const int32_t fc,
                              int32_t *basis, const int dquot,
                              int *mgb, int32_t start, const long elim){
  if(modgbs->nprimes >= modgbs->alloc-1){
    fprintf(stderr, "Not enough space in modgbs\n");
    exit(1);
  }
  modgbs->primes[modgbs->nprimes] = fc;

  len_t i, j, k, idx;

  len_t len   = 0;
  hm_t *hm    = NULL;

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
  for(i = start; i < modgbs->ld; i++){
    idx = bs->lmps[i];
    if (bs->hm[idx] == NULL) {
      fprintf(stderr, " poly is 0\n");
      free(evi);
      exit(1);
    } else {
      hm  = bs->hm[idx]+OFFSET;
      len = bs->hm[idx][LENGTH];
    }
    int bc = modgbs->modpolys[i]->len - 1;
    for (j = 1; j < len; ++j) {
      uint32_t c = bs->cf_32[bs->hm[idx][COEFFS]][j];
      for (k = 0; k < nv; ++k) {
          mgb[k] = ht->ev[hm[j]][evi[k]];
      }
      while(!is_equal_exponent_elim(mgb, basis + (bc * (nv - elim)), nv, elim)){
        bc--;
      }
      modgbs->modpolys[i]->cf_32[bc][modgbs->nprimes] = c;
      bc--;
    }
  }

  modgbs->nprimes++;

  free(evi);
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

static int32_t * gb_modular_trace_learning(gb_modpoly_t modgbs,
                                           int32_t *mgb,
                                           int32_t *num_gb,
                                           int32_t **leadmons,
                                           trace_t *trace,
                                           ht_t *tht,
                                           bs_t *bs_qq,
                                           ht_t *bht,
                                           md_t *st,
                                           const int32_t fc,
                                           int info_level,
                                           int print_gb,
                                           int *dim,
                                           long *dquot_ori,
                                           int32_t start,
                                           data_gens_ff_t *gens,
                                           int32_t maxbitsize,
                                           files_gb *files,
                                           int *success)
{
    double ca0, rt;
    ca0 = realtime();

    bs_t *bs = NULL;
    /* if(gens->field_char){ */
      int32_t err = 0;
      bs = core_gba(bs_qq, st, &err, fc);
      if (err) {
        printf("Problem with F4, stopped computation.\n");
        exit(1);
      }
      /*
      free_shared_hash_data(bht);
    }
    else{
      if(st->laopt > 40){
        bs = modular_f4(bs_qq, bht, st, fc);
      }
      else{
        bs = gba_trace_learning_phase(trace, tht, bs_qq, bht, st, fc);
      }
    } */
    rt = realtime()-ca0;

    if(info_level > 1){
        fprintf(stderr, "Learning phase %.2f Gops/sec\n",
                (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/rt);
    }
    if(info_level > 2){
        fprintf(stderr, "------------------------------------------\n");
        fprintf(stderr, "#ADDITIONS       %13lu\n", (unsigned long)st->trace_nr_add * 1000);
        fprintf(stderr, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->trace_nr_mult * 1000);
        fprintf(stderr, "#REDUCTIONS      %13lu\n", (unsigned long)st->trace_nr_red);
        fprintf(stderr, "------------------------------------------\n");
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

    long dquot = 0;
    int32_t *lmb;
    if(st->nev){
      lmb = monomial_basis_enlarged(num_gb[0], bht->nv - st->nev, 
                                    bexp_lm2, &dquot);
    }
    else{
      lmb = monomial_basis_enlarged(num_gb[0], bht->nv, 
                                    bexp_lm, &dquot);
    }

    /************************************************/
    /************************************************/

    /* int32_t *lens = array_of_lengths(bexp_lm, num_gb[0], lmb, dquot, bht->nv); */
    int32_t *lens = array_of_lengths(leadmons[0], num_gb[0], lmb, dquot, bht->nv - st->nev);

    gb_modpoly_init(modgbs, 2, lens, bht->nv - st->nev, num_gb[0], leadmons[0], lmb);

    modpgbs_set(modgbs, bs, bht, fc, lmb, dquot, mgb, start, st->nev);
    if(bs->lml == 1){
        if(info_level){
            fprintf(stderr, "Grobner basis has a single element\n");
        }
        int is_empty = 1;
        for(int i = 0; i < bht->nv; i++){
            if(bexp_lm[i]!=0){
                is_empty = 0;
            }
        }
        if(is_empty){
            *dquot_ori = 0;
            *dim = 0;
            /* if(info_level){ */
            /*   fprintf(stderr, "No solution\n"); */
            /* } */
            /* print_ff_basis_data( */
            /*                     files->out_file, "a", bs, bht, st, gens, print_gb); */
            free_basis_without_hash_table(&(bs));
            return NULL;
        }
    }

    /**************************************************/

    free_basis_without_hash_table(&(bs));
    return lmb;
}


static void gb_modular_trace_application(gb_modpoly_t modgbs,
                                         int32_t *mgb,
                                         int32_t *num_gb,
                                         int32_t **leadmons_ori,
                                         int32_t **leadmons_current,

                                         trace_t **btrace,
                                         ht_t **btht,
                                         bs_t *bs_qq,
                                         ht_t **bht,
                                         md_t *st,
                                         const int32_t fc,
                                         int info_level,
                                         bs_t **obs,
                                         int32_t *lmb_ori,
                                         int32_t dquot_ori,
                                         primes_t *lp,
                                         int32_t start,
                                         data_gens_ff_t *gens,
                                         double *stf4,
                                         int *bad_primes){

  double rt = realtime();
  st->info_level = 0;
  /* tracing phase */

  /* F4 and FGLM are run using a single thread */
  /* st->nthrds is reset to its original value afterwards */
  /*at the moment multi-threading is not supprted here*/
  memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));

  bs_t *bs = NULL;
  int32_t error = 0;
  bs = core_gba(bs_qq, st, &error, lp->p[0]);
  *stf4 = realtime()-rt;
  /* if(st->laopt > 40){
    bs = modular_f4(bs_qq, bht[0], st, lp->p[0]);
  }
  else{
    bs = gba_trace_application_phase(btrace[0], btht[0], bs_qq, bht[0], st, lp->p[0]);
  }
  */
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
        free_basis_without_hash_table(&bs);
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
    modpgbs_set(modgbs, bs, bht[0], lp->p[0], lmb_ori, dquot_ori, mgb, start, st->nev);
  }

  if (bs != NULL) {
    free_basis_without_hash_table(&bs);
  }

}

static inline void choose_coef_to_lift(gb_modpoly_t modgbs, data_lift_t dlift){
  uint32_t ld = modgbs->ld;
  for(int32_t i = 0; i < ld; i++){
    uint32_t d = 0;
    uint32_t len = modgbs->modpolys[i]->len;
    while(d < len - 1){
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
    uint32_t c = modgbs->modpolys[k]->cf_32[coef[k]][modgbs->nprimes  - 1 ];
    mpz_CRT_ui(dlift->crt[k], dlift->crt[k], mod_p[0],
               c, modgbs->primes[modgbs->nprimes - 1 ],
               prod_p[0], dlift->tmp, 1);
  }
  mpz_set(mod_p[0], prod_p[0]);
}


/* Incremental CRT on the whole array of witness coefficients */
/* mod is the current modulus */
static inline void incremental_dlift_crt_full(gb_modpoly_t modgbs, data_lift_t dl,
                                              int32_t *coef, mpz_t mod_p, mpz_t prod_p,
                                              int thrds){

  uint64_t newprime = modgbs->primes[modgbs->nprimes - 1 ];

  /* all primes are assumed to be good primes */
  mpz_mul_ui(prod_p, mod_p, (uint32_t)newprime);
  for(int32_t k = dl->lstart; k < modgbs->ld; k++){
    uint64_t c = modgbs->modpolys[k]->cf_32[coef[k]][modgbs->nprimes  - 1 ];

    mpz_CRT_ui(dl->crt[k], dl->crt[k], mod_p,
               c, newprime, prod_p, dl->tmp, 1);

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
    if(dlift->check1[k]){
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
    if(dl->check1[k]){
      mpz_set_ui(dl->tmp, 1);
      for(int32_t l = 0; l < polys[k]->len; l++){

        if(ratreconwden(rnum, rden, polys[k]->cf_zz[l], mod_p, dl->den[k], rd1)){
          mpz_set(polys[k]->cf_qq[2*l], rnum);
          mpz_set(polys[k]->cf_qq[2*l + 1], rden);
          mpz_lcm(dl->tmp, dl->tmp, rden);
        }
        else{
          /* fprintf(stderr, "[%d/%d]", k, modgbs->ld - 1); */
          mpz_set_ui(dl->gden, 1);
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
      /* dl->S++; */
    }
    else{
      mpz_clear(rnum);
      mpz_clear(rden);
      return k;
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

static inline int verif_lifted_rational_wcoef(gb_modpoly_t modgbs, data_lift_t dl,
                                              int thrds){
  for(int32_t k = dl->lstart; k < dl->lend; k++){
    if(!dl->check1[k]){
      /* too early to perform the verification */
      return k;
    }
    for(int i = 0; i < thrds; i++){

      uint32_t prime = modgbs->primes[modgbs->nprimes - (thrds - i) ];
      uint32_t coef = modgbs->modpolys[k]->cf_32[dl->coef[k]][modgbs->nprimes  - (thrds - i) ];
      int b = verif_coef(dl->num[k], dl->den[k], prime, coef);

      if(!b){
        dl->check1[k] = 0;
        return k;
      }
    }
    
    dl->start++;
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
                        int thrds, double *st_crt, double *st_rrec){

  verif_lifted_rational_wcoef(modgbs, dl, thrds);
  if(dl->lstart != dl->start){
    double st = realtime();
    crt_lift_modgbs(modgbs, dl, dl->lstart, dl->start);
    *st_crt += realtime() - st;
    st = realtime();
    dl->start = ratrecon_lift_modgbs(modgbs, dl, dl->lstart, dl->start, mod_p, recdata1, recdata2);
    *st_rrec += realtime()-st;
  }
  dl->lstart = dl->start;

  double st = realtime();

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
      int b = reconstructcoeff(dl, i, mod_p,
                               recdata1, recdata2);
      if(!b){
        break;
      }
      else{
        dl->check1[i] = 1;
      }
    }
  }
  *st_rrec += realtime()-st;

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

int msolve_gbtrace_qq(
                      gb_modpoly_t modgbs,
                      data_gens_ff_t *gens,
                      msflags_t flags){

  double st0 = realtime();

  int *dim_ptr = &flags->dim;
  long *dquot_ptr = &flags->dquot;
  int32_t ht_size = flags->ht_size;
  int32_t nr_threads = flags->nr_threads;
  int32_t max_nr_pairs = flags->max_nr_pairs;
  int32_t elim_block_len = flags->elim_block_len;
  int32_t reset_ht = flags->reset_ht;
  int32_t la_option = flags->la_option;
  int32_t use_signatures = flags->use_signatures;
  int32_t info_level = flags->info_level;
  int32_t pbm_file = flags->pbm_file;
  int32_t print_gb = flags->print_gb;
  files_gb *files = flags->files;

  uint32_t field_char = gens->field_char;
  const void *cfs = gens->mpz_cfs;
  if(gens->field_char){
    cfs = gens->cfs;
  }
  else{
    cfs = gens->mpz_cfs;
  }
  int mon_order = 0;
  int32_t nr_vars = gens->nvars;
  int32_t nr_gens = gens->ngens;
  int reduce_gb = 1;
  int32_t nr_nf = 0;
  const uint32_t prime_start = pow(2, 30);

  /* initialize stuff */
  md_t *st  = allocate_meta_data();

  int *invalid_gens   =   NULL;
  int res = validate_input_data(&invalid_gens, cfs, gens->lens, &field_char, &mon_order,
                                &elim_block_len, &nr_vars, &nr_gens, &nr_nf, &ht_size, &nr_threads,
                                &max_nr_pairs, &reset_ht, &la_option, &use_signatures, &reduce_gb,
                                &info_level);

  /* all data is corrupt */
  if (res == -1) {
    fprintf(stderr, "Invalid input generators, msolve now terminates.\n");
    free(invalid_gens);
    return -3;
  }
  /* checks and set all meta data. if a nonzero value is returned then
   * some of the input data is corrupted. */

  if (check_and_set_meta_data_trace(st, gens->lens, gens->exps, cfs, invalid_gens,
                                    field_char, mon_order, elim_block_len, nr_vars, nr_gens,
                                    nr_nf, ht_size, nr_threads, max_nr_pairs, reset_ht, la_option,
                                    use_signatures, reduce_gb, prime_start,
                                    nr_threads /* nr_primes */,
                                    pbm_file, info_level)) {
    fprintf(stderr, "Something went wrong when checking and setting meta data, msolve now terminates.\n");
    free(st);
    return -3;
  }

  mstrace_t msd;
  initialize_mstrace(msd, st);

  /* read in ideal, move coefficients to integers */
  import_input_data(msd->bs_qq, st, 0, st->ngens_input, gens->lens, gens->exps, cfs, invalid_gens);
  free(invalid_gens);
  invalid_gens  =   NULL;

  print_initial_statistics(stderr, st);

  /* for faster divisibility checks, needs to be done after we have
    * read some input data for applying heuristics */
  calculate_divmask(msd->bht);

  /* sort initial elements, smallest lead term first */
  sort_r(msd->bs_qq->hm, (unsigned long)msd->bs_qq->ld, sizeof(hm_t *),
          initial_input_cmp, msd->bht);

  if(gens->field_char == 0){
    remove_content_of_initial_basis(msd->bs_qq);
    /* generate lucky prime numbers */
    generate_lucky_primes(msd->lp, msd->bs_qq, st->prime_start, st->nthrds);
  }
  else{
    msd->lp->old = 0;
    msd->lp->ld = 1;
    msd->lp->p = calloc(1, sizeof(uint32_t));
    normalize_initial_basis(msd->bs_qq, st->gfc);
  }

  uint32_t prime = 0;
  uint32_t primeinit = 0;
  uint32_t lprime = 1303905299;
  srand(time(0));

  prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  while(gens->field_char==0 && is_lucky_prime_ui(prime, msd->bs_qq)){
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  }

  primeinit = prime;
  msd->lp->p[0] = primeinit;
  if(gens->field_char){
    msd->lp->p[0] = gens->field_char;
    primeinit = gens->field_char;
  }

  int success = 1;

  int32_t maxbitsize = maxbitsize_gens(gens, st->ngens);

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

  uint32_t nbadprimes = 0;

  while(learn){

    int32_t *lmb_ori = gb_modular_trace_learning(modgbs,
                                                 msd->mgb,
                                                 msd->num_gb, msd->leadmons_ori,
                                                 msd->btrace[0],
                                                 msd->tht, msd->bs_qq, msd->bht, st,
                                                 msd->lp->p[0],
                                                 info_level,
                                                 print_gb,
                                                 dim_ptr, dquot_ptr,
                                                 0,
                                                 gens, maxbitsize,
                                                 files,
                                                 &success);
    /*lmb_ori can be NULL when gb = [1]*/
    if(lmb_ori == NULL || print_gb == 1){
      if(dlinit){
        data_lift_clear(dlift);
      }
      free_mstrace(msd, st);
      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      free(st);
      return 0;
    }

    apply = 1;

    gb_modpoly_realloc(modgbs, 1, dlift->S);

#ifdef DEBUGGBLIFT
    display_gbmodpoly_cf_32(stderr, modgbs);
#endif

    if(!dlinit){
      int nb = 0;
      int32_t *ldeg = array_nbdegrees((*msd->leadmons_ori), msd->num_gb[0],
                                      msd->bht->nv - st->nev, &nb);
      data_lift_init(dlift, modgbs->ld, ldeg, nb);
      choose_coef_to_lift(modgbs, dlift);
      free(ldeg);
      dlinit = 1;
    }

    if(info_level){
      int s= 0;
      for(int i = 0; i < dlift->nsteps; i++){
        fprintf(stderr, "[%d]", dlift->steps[i]);
        s+=dlift->steps[i];
      }
      fprintf(stderr, "\n");
      if(s > 1){
        fprintf(stderr, "%d polynomials to lift\n", s);
      }
    }

    if(lmb_ori == NULL || success == 0 || gens->field_char) {

      apply = 0;
      if(dlinit){
        data_lift_clear(dlift);
      }

      gb_modpoly_clear(modgbs);

      free_mstrace(msd, st);
      free_rrec_data(recdata1);
      free_rrec_data(recdata2);

      free(st);
      fprintf(stderr, "Something went wrong in the learning phase, msolve restarts.");
      return msolve_gbtrace_qq(modgbs, gens, flags);

    }
    /* duplicate data for multi-threaded multi-mod computation */
    duplicate_data_mthread_gbtrace(st->nthrds, msd->bs_qq, st, msd->num_gb,
                                   msd->leadmons_ori, msd->leadmons_current,
                                   msd->btrace);

    /* copy of hash tables for tracer application */
    msd->blht[0] = msd->bht;
    for(int i = 1; i < st->nthrds; i++){
      ht_t *lht = copy_hash_table(msd->bht, st);
      msd->blht[i] = lht;
    }

    if(info_level){
      fprintf(stderr, "\nStarts multi-modular computations\n");
    }

    learn = 0;
    while(apply){

      prime = next_prime(prime);
      if(prime >= lprime){
        prime = next_prime(1<<30);
      }
      /* generate lucky prime numbers */
      msd->lp->p[0] = prime;
      while(is_lucky_prime_ui(prime, msd->bs_qq) || prime==primeinit){
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
        while(is_lucky_prime_ui(prime, msd->bs_qq) || prime==primeinit){
          prime = next_prime(prime);
          if(prime >= lprime){
            prime = next_prime(1<<30);
          }
          msd->lp->p[i] = prime;
        }
      }
      prime = msd->lp->p[nthrds /* st->nthrds */ - 1];

      if(modgbs->alloc <= nprimes + 2){
        gb_modpoly_realloc(modgbs, 16*st->nthrds, dlift->S);
      }

      gb_modular_trace_application(modgbs, msd->mgb,
                                   msd->num_gb,
                                   msd->leadmons_ori,
                                   msd->leadmons_current,
                                   msd->btrace,
                                   msd->btht, msd->bs_qq, msd->blht, st,
                                   field_char, 0, /* info_level, */
                                   msd->bs, lmb_ori, *dquot_ptr, msd->lp,
                                   dlift->S, gens, &stf4, msd->bad_primes);


      /* nprimes += st->nthrds; */
      nprimes += 1; /* at the moment, multi-mod comp is mono-threaded */

      if(nprimes == 1){
        if(info_level>2){
          fprintf(stderr, "------------------------------------------\n");
          fprintf(stderr, "#ADDITIONS       %13lu\n", (unsigned long)st->application_nr_add * 1000);
          fprintf(stderr, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->application_nr_mult * 1000);
          fprintf(stderr, "#REDUCTIONS      %13lu\n", (unsigned long)st->application_nr_red);
          fprintf(stderr, "------------------------------------------\n");
        }
        if(info_level>1){
          fprintf(stderr, "Application phase %.2f Gops/sec\n",
                  (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4));
          fprintf(stderr, "Elapsed time: %.2f\n", stf4);
        }
      }
      int bad = 0;
      for(int i = 0; i < nthrds/* st->nthrds */; i++){
        if(msd->bad_primes[i] == 1){
          bad = 1;
          nbadprimes++;
        }
      }

      if(nbadprimes == nprimes){
        fprintf(stderr, "Too many bad primes, computation will restart\n");
        free_mstrace(msd, st);
        if(dlinit){
          data_lift_clear(dlift);
        }

        free_rrec_data(recdata1);
        free_rrec_data(recdata2);

        free(st);
        return msolve_gbtrace_qq(modgbs, gens, flags);

      }

      int lstart = dlift->lstart;
      double ost_rrec = st_rrec;

      if(!bad){
        ratrecon_gb(modgbs, dlift, msd->mod_p, msd->prod_p, recdata1, recdata2,
                    nthrds/* st->nthrds */, &st_crt, &st_rrec);
      }
      if(/* (st_crt -ost_crt) + */ (st_rrec - ost_rrec) > dlift->rr * stf4){
        dlift->rr = 2*dlift->rr;
        if(info_level){
          fprintf(stderr, "(->%d)", dlift->rr);
        }
      }
      if(info_level){
        if(!(nprimes & (nprimes - 1))){
          fprintf(stderr, "{%d}", nprimes);
        }
      }
      if(dlift->lstart != lstart && dlift->lstart < modgbs->ld - 1){
        if(info_level){
          fprintf(stderr, "<%.2f%%>", 100* (float)(dlift->lstart + 1)/modgbs->ld);
        }
        lstart = dlift->lstart;
      }
      if(dlift->lstart >= modgbs->ld){
        if(info_level){
          fprintf(stderr, "<100%%>\n");
          fprintf(stderr, "CRT time = %.2f, Rational reconstruction time = %.2f\n", st_crt, st_rrec);
        }
        apply = 0;
      }
      /* this is where learn could be reset to 1 */
      /* but then duplicated datas and others should be free-ed */
    }
  }
  if(info_level){
    long nbits = max_bit_size_gb(modgbs);
    fprintf(stderr, "Maximum bit size of the coefficients: %ld\n", nbits);
    fprintf(stderr, "%d primes used. \nElapsed time: %.2f\n", nprimes, realtime()-st0);
  }
  free_mstrace(msd, st);
  if(dlinit){
    data_lift_clear(dlift);
  }

  free_rrec_data(recdata1);
  free_rrec_data(recdata2);

  free(st);

  return 0;
}

/*
 Function which is called by core_msolve
 */

void print_msolve_gbtrace_qq(data_gens_ff_t *gens,
                            msflags_t flags){
  gb_modpoly_t modgbs;

  msolve_gbtrace_qq(modgbs, gens, flags);

  FILE *ofile;
  if (flags->files->out_file != NULL) {
    ofile = fopen(flags->files->out_file, "w+");
  } else {
    ofile = stdout;
  }
  if (flags->print_gb == 1) {
    fprintf(ofile, "#Leading ideal data\n");
  } else {
    if (flags->print_gb > 1) {
      fprintf(ofile, "#Reduced Groebner basis data\n");
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
  if (modgbs->ld == 1) {
    fprintf(ofile, "#length of basis:      1 element\n");
  } else {
    fprintf(ofile, "#length of basis:      %u elements sorted by increasing leading monomials\n", modgbs->ld);
  }
  fprintf(ofile, "#---\n");
  if (flags->files->out_file != NULL) {
    fclose(ofile);
  }

  if(flags->print_gb > 1){

    if(flags->files->out_file != NULL){
      FILE *ofile = fopen(flags->files->out_file, "a+");
      display_gbmodpoly_cf_qq(ofile, modgbs, gens);
      fclose(ofile);
    }
    else{
      display_gbmodpoly_cf_qq(stdout, modgbs, gens);
    }
  }
  if(flags->print_gb == 1){
    if(flags->files->out_file != NULL){
      FILE *ofile = fopen(flags->files->out_file, "a+");
      display_lm_gbmodpoly_cf_qq(ofile, modgbs, gens);
      fclose(ofile);
    }
    else{
      display_lm_gbmodpoly_cf_qq(stdout, modgbs, gens);
    }
  }
  gb_modpoly_clear(modgbs);

}
