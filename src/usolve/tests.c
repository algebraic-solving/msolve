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

#include<stdio.h>
#include<stdlib.h>
#include <unistd.h>
#include<gmp.h>
#include<math.h>
#include<time.h>
#ifdef _OPENMP
#include<omp.h>
#endif


#include"taylor_shift.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/nmod_poly.h"

#undef THRESHOLDSHIFT
#define THRESHOLDSHIFT 512
#define BITSPRIME 64

static inline void compute_shift_pwx(mpz_t **shift_pwx,
                                            const unsigned long int deg,
                                            const unsigned long int npwr,
                                            unsigned long int pwx,
                                            const unsigned int nthreads){
  taylorshift1_naive(shift_pwx[0], pwx);
  for(unsigned long int i=1; i<npwr; i++){
    /* mpz_poly_mul(shift_pwx[i], */
    /*              shift_pwx[i-1], pwx, */
    /*              shift_pwx[i-1], pwx, */
    /*              nthreads); */
    //    mpz_poly_print_maple(shift_pwx[i-1], pwx);
    mpz_poly_mul(shift_pwx[i],
                 shift_pwx[i-1], pwx,
                 shift_pwx[i-1], pwx,
                 nthreads);
    //    mpz_poly_print_maple(shift_pwx[i], 2*pwx);

    pwx = 2 * pwx;
  }
}

static void mpz_random_dense_poly(mpz_t *pol, unsigned long int deg, mp_bitcnt_t bits){
  unsigned long i;
  mpz_t coeff;
  mpz_init(coeff);

  clock_t t=clock();
  srand((int)t);

  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, (int)t);

  int sgn;
  for(i=0;i<=deg;i++){
    mpz_urandomb(coeff, state , bits);
    if (rand()%2==0)
      sgn =1;
    else
      sgn =-1;
    mpz_mul_si(coeff, coeff, sgn);
    mpz_set(pol[i], coeff);
  }
  mpz_clear(coeff);
}

//a renseigner
static unsigned long int compute_degpower(const unsigned long int deg){
  unsigned long int degpower=deg;

  while(degpower-1>=THRESHOLDSHIFT){
    degpower = degpower / 2;
  }

  return degpower;
}


//On alloue de quoi stocker X^alpha, X^(2alpha), X^(4*alpha) etc.
//avec alpha = THRESHOLDSHIFT
static inline void allocate_shift_pwx(mpz_t **shift_pwx,
                                             unsigned long int npwr,
                                             unsigned long int pwx){
  //pwx c'est le degre du premier element de shift_pwx
  unsigned long int newpwx=pwx;
  mp_bitcnt_t nbits=LOG2(newpwx);
  for(int i=0;i<npwr;i++){
    shift_pwx[i]=(mpz_t *)malloc(sizeof(mpz_t)*(newpwx+1));
    for(int j=0;j<=newpwx;j++){
      mpz_init2(shift_pwx[i][j], nbits);
    }
    newpwx=2*newpwx;
    nbits=LOG2(newpwx);
  }
  mpz_set_ui(shift_pwx[0][pwx],1);
}


static void initialize_heap_flags(usolve_flags *flags, unsigned long int deg){
  if(flags->classical_algo==0){
    flags->cur_deg = deg;
    flags->pwx = compute_degpower(flags->cur_deg);
    if(flags->pwx>=flags->cur_deg){
      flags->nblocks = 0;
      flags->npwr=0;
    }
    else{
      flags->nblocks = 1<<LOG2(deg/flags->pwx);
      flags->npwr=LOG2(flags->nblocks);
    }
    if(flags->npwr>0){
      flags->shift_pwx=(mpz_t **)malloc(sizeof(mpz_t *)*(flags->npwr));
      allocate_shift_pwx(flags->shift_pwx,
                                flags->npwr,
                                flags->pwx);

      compute_shift_pwx(flags->shift_pwx,
                               flags->cur_deg,
                               flags->npwr,
                               flags->pwx,
                               flags->nthreads);
    }
    else{
      flags->shift_pwx = NULL;
    }
    flags->tmpol = (mpz_t *)(malloc(sizeof(mpz_t)*(deg+1)));
    for(int i=0; i<=deg; i++){
      mpz_init(flags->tmpol[i]);
    }
    flags->tmpol_desc = (mpz_t *)(malloc(sizeof(mpz_t)*(deg+1)));
    for(int i=0; i<=deg; i++){
      mpz_init(flags->tmpol_desc[i]);
    }
  }
}

static void free_heap_flags(usolve_flags *flags, unsigned long int deg){
  if(flags->classical_algo == 0){
    for(unsigned long int i=0; i<=deg;i++){
      mpz_clear(flags->tmpol_desc[i]);
      mpz_clear(flags->tmpol[i]);
    }
  }
}

static void initialize_flags(usolve_flags *flags){
  flags->search = 0;
  flags->prec_isole = 16;
  flags->precision_loss = 0;

  flags->bound_pos = 0;
  flags->bound_neg = 0;
  flags->sign = 0;
  flags->revert = 1;

  flags->hasrealroots = 0;

  flags->transl = 0;
  flags->node_looked = 0;
  flags->half_done = 0;

  flags->cur_deg = 0;
  flags->pwx = 0;
  flags->nblocks = 0;
  flags->npwr = 0;
  flags->shift_pwx = NULL;
  flags->tmpol = NULL;
  flags->tmpol_desc = NULL;
  flags->time_desc = 0;
  flags->time_shift = 0;
  flags->nthreads = 1;
  flags->verbose = 0;
  flags->binary_file = 0;
  flags->classical_algo = 0;
  flags->print_stats = 0;
  flags->debug = 0;
}

/* static unsigned long int mpz_poly_max_bsize_coeffs(mpz_t *upoly, unsigned long int deg){ */
/*   unsigned long int max = 0, bs; */
/*   for(int i=0 ; i<=deg; i++){ */
/*     bs = ilog2_mpz(upoly[i]); */
/*     if(bs>max){ */
/*       max = bs; */
/*     } */
/*   } */
/*   return max; */
/* } */

void
mpz_multi_mod_ui(mp_limb_t * out, const mpz_t in, mp_srcptr primes,
                           long num_primes)
{
  long i;
  for (i = 0; i < num_primes; i++)
    {
      out[i] = mpz_fdiv_ui(in, primes[i]);
    }
}

static inline
mp_limb_t mynmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
{
  const mp_limb_t neg = mod.n - a;
  return ((neg > b) ? (a + b) : b - neg);
  /* if (neg > b) */
  /*   return a + b; */
  /* else  */
  /*   return b - neg; */
}

void nmod_poly_taylor_shift_naive(mp_ptr residue, unsigned long int deg, nmod_t mod){
  long int i, j;
  for(i = 0; i <=deg -1; i++){
    for(j = deg - 1; j >= i; j--){
      residue[j] = mynmod_add(residue[j], residue[j+1], mod);
    }
  }
}

void nmod_upoly_taylor_shift(mp_ptr residue, unsigned long int deg,
                            mp_ptr tmp,
                            mp_ptr *residue_shift_pwx,
                            unsigned long int pwx,
                            nmod_t mod){

  int l;
  long int nblocks;
  if(deg<=pwx){
    nblocks = 0;
  }
  else{
    l =LOG2(deg / pwx);
    nblocks = 1<<l;
  }
  if(nblocks<=1){
    nmod_poly_taylor_shift_naive(residue, deg, mod);
  }

  unsigned long int fdegree = deg - (nblocks - 1)*pwx;
  int i;

  for(long int b=0;b<nblocks-1;b++){
      nmod_poly_taylor_shift_naive(&residue[b * pwx], pwx-1, mod);
  }
  nmod_poly_taylor_shift_naive(residue + (nblocks-1) * (pwx), fdegree, mod);

  unsigned long int npowers = LOG2(nblocks);
  mp_ptr cur_shifted_pX = NULL;

  for(i = 0; i < npowers - 1; i++){
    cur_shifted_pX = residue_shift_pwx[i];
    if(nblocks>1){
      fdegree = deg - (nblocks -1) * pwx;
    }
    else{
      fdegree = deg - pwx;
    }
    nblocks = nblocks / 2;
    for(int cnt = 1; cnt <= nblocks; cnt++){
      unsigned long int ndeg;
      if(cnt==nblocks){
        ndeg = fdegree;
      }
      else{
        ndeg = pwx - 1;
      }
      _nmod_poly_mul(&tmp[2*(cnt-1)*pwx],
                     &residue[2*(cnt-1)*pwx + pwx],
                     ndeg + 1,
                     cur_shifted_pX,
                     pwx + 1,
                     mod);
      _nmod_poly_add(&residue[2*(cnt-1)*pwx],
                     &residue[2*(cnt-1)*pwx], pwx,
                     &tmp[2*(cnt-1)*pwx], pwx, mod);
    }
    pwx = 2 * pwx;
  }
  return;
}


//NON FINIE CAR UNE ANALYSE DE COMPLEXITE MONTRE QUE CE N'EST PAS UTILE
void multi_mod_taylor_shift(mpz_t *pol, unsigned long int deg,
                            mpz_t **shift_pwx, unsigned long int pwx,
                            int nthreads){

  unsigned long int nbits = USOLVEmpz_poly_max_bsize_coeffs(pol, deg);
  fprintf(stderr, "nbits = %lu\n", nbits);
  nbits = nbits + deg + 1 + LOG2(deg + 1);
  fprintf(stderr, "nbits = %lu\n", nbits);

  /* gets primes between 2^63 and 2^64 */
  unsigned long int nprimes = (nbits + (BITSPRIME - 1) - 1) / (BITSPRIME - 1);
  mp_ptr primes = malloc(sizeof(mp_limb_t) * nprimes);
  fprintf(stderr, "%lu\n", UWORD(1) << (BITSPRIME - 1));
  primes[0] = n_nextprime(UWORD(1) << (BITSPRIME - 1), 1);
  for (unsigned long int i = 1; i < nprimes; i++)
    primes[i] = n_nextprime(primes[i-1], 1);

  //allocates space for pol modulo the primes
  mp_ptr *residues = malloc(sizeof(mp_ptr) * nprimes);
  for(unsigned long int i = 0; i < nprimes; i++){
    residues[i] = malloc(sizeof(mp_limb_t) * (deg + 1));
  }

  fprintf(stderr, "nprimes = %lu\n", nprimes);
  //allocates space for shiftedpowers of X
  unsigned long int nblocks = 1<<LOG2(deg/pwx);
  unsigned long int npowers=LOG2(nblocks);
  mp_ptr **residues_shift_pwx = malloc(sizeof(mp_ptr *) * nprimes);
  for(long int i=0; i < nprimes; i++){
    residues_shift_pwx[i] = malloc(sizeof(mp_ptr) * npowers);
    unsigned long int pX = pwx;
    for(unsigned long int j = 0; j < npowers; j++){
      residues_shift_pwx[i][j] = malloc(sizeof(mp_limb_t) * (pX + 1));
      for(unsigned long int k = 0; k <= pX; k++){
        residues_shift_pwx[i][j][k] = mpz_fdiv_ui(shift_pwx[j][k], primes[i]);
      }
      pX = 2 * pX;
    }
  }
  fprintf(stderr, "Modular residues (shifted powers of X) obtained\n");
  mp_ptr tmp;

  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;

  tmp = malloc(sizeof(mp_limb_t) * nprimes);
  fmpz_comb_init(comb, primes, nprimes);
  fmpz_comb_temp_init(comb_temp, comb);

  mp_ptr tmp_pol = malloc(sizeof(mp_limb_t) * (deg + 1));
  unsigned long int length = deg + 1;
  fprintf(stderr, "Allocation done\n");
  double e = realtime();
  for(unsigned long int j = 0; j < length; j++){
    //    fmpz_multi_mod_ui(tmp, pol + i, comb, comb_temp);
    //    mpz_multi_mod_ui(tmp, pol[i], comb->primes, comb->num_primes);
    for (unsigned long i = 0; i < nprimes; i++)
      //      residues[i][j] = mpz_fdiv_ui(pol[j], primes[i]);
      residues[i][j] = mpz_fdiv_ui(pol[j], primes[i]);
  }
  fprintf(stderr, "Residues loaded (%f sec)\n", realtime() - e);
  for(unsigned long int i = 0; i < nprimes; i++){
    nmod_t mod;
    mp_limb_t p;

    p = primes[i];
    nmod_init(&mod, p);
    nmod_upoly_taylor_shift(residues[i], deg, tmp_pol, residues_shift_pwx[i], pwx, mod);
  }
  fprintf(stderr, "Time multimod = %f\n", realtime() - e);
  free(tmp);
  for(unsigned long int i = 0; i < nprimes; i++){
    free(residues[i]);
  }
  free(residues);
  free(primes);

  return;
}

/*

  pol = pols[0] + 2^sb pols [1] + 2^(2*sb) pols[2] + ... + 2^(sb * (nthreads - 1)) * pols[nthreads - 1]

*/
void split_coeffs_poly(mpz_t **pols, mpz_t * pol, unsigned long int deg, long int nbits, int nthreads){
  for(int j = nthreads - 1; j >= 0; j--){
    long int trunc =  j * nbits / nthreads ;
#pragma omp parallel for num_threads(nthreads)
    for(unsigned long int i = 0; i <= deg; i++){
      mpz_tdiv_q_2exp(pols[j][i], pol[i], trunc);
      mpz_tdiv_r_2exp(pol[i], pol[i], trunc);
    }
  }
}


//tbr == to be removed
//Ici, on recupere la taille binaire max "nbits" de pol
//ensuite on ecrit pol = pol0+ 2^q * pol1 + ...+ 2^(q * nthreads) * pol_n
//ou q = nbits / nthreads
void tbr_taylor_shift(mpz_t *pol, unsigned long int deg, mpz_t **tmp, mpz_t **pols,
                  mpz_t **shift_pwx, unsigned long int pwx,
                  int nthreads){
  long int nbits = USOLVEmpz_poly_max_bsize_coeffs(pol, deg);
  //  mpz_poly_print_maple(pol, deg);
  //  long int sb = 2 * nbits / (nthreads * deg);
  //  fprintf(stderr, "[sb = %ld]\n", sb);
  if(nthreads >= 1){
    //    fprintf(stderr, "New algo starts\n");
    split_coeffs_poly(pols, pol, deg, nbits, nthreads);
    //double e = realtime();
#pragma omp parallel for num_threads(nthreads)
    for(int j = 0; j < nthreads; j++){
      //      fprintf(stderr, "bs = %ld\n", mpz_poly_max_bsize_coeffs(pols[j], deg));
      //      mpz_poly_print_maple(pols[j], deg);
      taylorshift1_dac(pols[j], deg, tmp[j], shift_pwx, pwx, 0);
    }
    //    fprintf(stderr, "time %f\n", realtime() - e);
    //    e = realtime();
#pragma omp parallel for num_threads(nthreads)
    for(int j = 0; j < nthreads; j++){
      int trunc = j * nbits / nthreads;
      for(unsigned long int i = 0; i <= deg; i++){
        mpz_mul_2exp(pols[j][i], pols[j][i], trunc);
      }
    }
    //    fprintf(stderr, "Rescaling done (time = %f)\n", realtime() - e);
#pragma omp parallel for num_threads(nthreads)
    for(unsigned long int i = 0; i <= deg; i++){
      mpz_set_ui(pol[i], 0);
      for(int j = 0; j < nthreads; j++){
        mpz_add(pol[i], pol[i], pols[j][i]);
      }
    }
/*     for(unsigned long int i = 0; i <= deg; i++){ */
/*       mpz_set_ui(pol[i], 0); */
/*     } */
/*     for(int j = 0; j < nthreads; j++){ */
/* #pragma omp parallel for num_threads(nthreads) */
/*       for(unsigned long int i = 0; i <= deg; i++) */
/*         mpz_add(pol[i], pol[i], pols[j][i]); */
/*       } */
    }
  else{
    taylorshift1_dac(pol, deg, tmp[0], shift_pwx, pwx, nthreads);
  }
}

static inline void decompose_bits_poly(mpz_t **pols, mpz_t *pol, unsigned long int deg,
                         long int cutoff, long int nbits, long int nblocks){
  for(long int j = nblocks - 1; j >= 0; j--){
    long int trunc = nbits - cutoff * (nblocks - j);
    if(trunc < 0)
      trunc = 0;


    //#pragma omp parallel for num_threads(nblocks)
    for(unsigned long int i = 0; i <= deg; i++){
      mpz_tdiv_q_2exp(pols[j][i], pol[i], trunc);
      mpz_tdiv_r_2exp(pol[i], pol[i], trunc);
    }
  }
}

void taylor_shift(mpz_t *pol, unsigned long int deg, mpz_t **tmp, mpz_t **pols,
                  mpz_t **shift_pwx, unsigned long int pwx,
                  int nthreads){
  long int nbits = USOLVEmpz_poly_max_bsize_coeffs(pol, deg);
  long int cutoff = 2 * deg ;

  long int nblocks = nbits / cutoff + 1;
  double e = 0, tmp_e;
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_set_ui(tmp[nthreads][i], 0);
  }

  if(nthreads >= 1 && nblocks > 0){

    if(nblocks <= nthreads){
      fprintf(stderr, "nblocks (%ld) <= nthreads (%d)\n", nblocks, nthreads);
      //      long int nt = nthreads / nblocks;
      //pol = pols[0] + 2^(2*deg) * pols[1] + 2^(4*deg) * pols[2] + 2^(6*deg) *pols[3] +...+2^(nblocks*2*deg)*pols[nblocks]
      decompose_bits_poly(pols, pol, deg, cutoff, nbits, nblocks);

      //      e = realtime();

      //#pragma omp parallel
      {
        //#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
#pragma omp parallel for schedule(dynamic)
      for(int j = 0; j < nblocks; j++){
      //        fprintf(stderr, "bs = %ld\n", mpz_poly_max_bsize_coeffs(pols[j], deg));
        mpz_t *local_pol = pols[j];
        mpz_t *local_tmp = tmp[j];
        taylorshift1_dac(local_pol, deg, local_tmp, shift_pwx, pwx, 1);
      }
      }
      //      fprintf(stderr, "time %f\n", realtime() - e);
      //      mpz_poly_print_maple(pols[0], deg);
      //      e = realtime();
      for(long int j = nblocks-1; j >=0; j--){
        long int trunc = nbits - cutoff * (nblocks - j);
        if(trunc < 0)
          trunc = 0;
        for(unsigned long int i = 0; i <= deg; i++){
          mpz_mul_2exp(pols[j][i], pols[j][i], trunc);
        }
      }

      mpz_t *local_tmp = tmp[nthreads];
      for(int j = 0; j < nblocks; j++){
        //#pragma omp parallel for num_threads(nthreads)
        for(unsigned long int i = 0; i <= deg; i++){
          mpz_add(local_tmp[i], local_tmp[i], pols[j][i]);
        }
      }
      for(unsigned long int i = 0; i <= deg; i++){
        mpz_set(pol[i], local_tmp[i]);
      }
      return;
    }
    else{//nblocks > nthreads
      fprintf(stderr, "nblocks (%ld) > nthreads (%d)\n", nblocks, nthreads);
      //nblocks = nbits / cutoff
      long int q = nblocks / nthreads ;//pb quand nthreads = 1 (doit prendre q+1)
      if(nthreads==1){
        //        q++;
      }
      while(q){

        tmp_e = realtime();
        decompose_bits_poly(pols, pol, deg, cutoff, nbits, nthreads);

        e += (realtime() - tmp_e);
        //        tmp_e = realtime();
#pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < nthreads; j++){
          taylorshift1_dac(pols[j], deg, tmp[j], shift_pwx, pwx, 1);
        }

        for(long int j = nthreads - 1; j >= 0; j--){
          long int trunc = nbits -  cutoff * (nthreads - j);
          if(trunc < 0)
            trunc = 0;
#pragma omp parallel for num_threads(nthreads)
          for(unsigned long int i = 0; i <= deg; i++){
            mpz_mul_2exp(pols[j][i], pols[j][i], trunc);
          }
        }
#pragma omp parallel for num_threads(nthreads)
        for(unsigned long int i = 0; i <= deg; i++){
          for(int j = 0; j < nthreads; j++){
            mpz_add(tmp[nthreads][i], tmp[nthreads][i], pols[j][i]);
          }
        }
        //        e+=realtime() - tmp_e;

        nbits -= cutoff * nthreads;
        q--;
      }
      //      fprintf(stderr, "time = %f\n", e);

      long rth = nblocks % nthreads + 1;

      if(rth==0){
        for(unsigned long int i = 0; i <= deg; i++){
          mpz_swap(pol[i], tmp[nthreads][i]);
        }
        return;
      }

      tmp_e = realtime();
      decompose_bits_poly(pols, pol, deg, cutoff, nbits, rth);
      e += (realtime() - tmp_e);
      //      fprintf(stderr, "time in decompose_bits_poly = %f\n", e);
      //      e = realtime();
      //#pragma omp parallel for num_threads(rth)
      for(int j = 0; j < rth; j++){
        taylorshift1_dac(pols[j], deg, tmp[j], shift_pwx, pwx, nthreads);
      }

      for(int j = rth - 1; j >= 0; j--){
        int trunc = nbits - cutoff * (rth - j);
        if(trunc < 0)
          trunc = 0;
        for(unsigned long int i = 0; i <= deg; i++){
          mpz_mul_2exp(pols[j][i], pols[j][i], trunc);
        }
      }
#pragma omp parallel for num_threads(nthreads)
      for(unsigned long int i = 0; i <= deg; i++){
        for(int j = 0; j < rth; j++){
          mpz_add(tmp[nthreads][i], tmp[nthreads][i], pols[j][i]);
        }
        mpz_swap(pol[i], tmp[nthreads][i]);
      }
      return;
    }
  }
  else{
    taylorshift1_dac(pol, deg, tmp[0], shift_pwx, pwx, nthreads);
  }
}


int main(int argc, char **argv){
  if(argc!=5){
    fprintf(stderr, "Runs taylor shifts and measures timings\n");
    fprintf(stderr, "Usage is %s <degree> <nbits> <nthreads> <nloops>\n", argv[0]);
    exit(1);
  }
  unsigned long int deg = atoi(argv[1]);
  unsigned long int nbits = atoi(argv[2]);
  int nthreads = atoi(argv[3]);

  mpz_t * pol1 = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
  mpz_t * pol2 = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
  mpz_t * pol3 = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
  mpz_t * tmpol = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_init(pol1[i]);
    mpz_init(pol2[i]);
    mpz_init(pol3[i]);
    mpz_init2(tmpol[i], deg + nbits);
  }
  usolve_flags *flags = (usolve_flags*)(malloc(sizeof(usolve_flags)));
  fprintf(stderr, "Starting initialization of data\n");
  double t = realtime();
  initialize_flags(flags);

  initialize_heap_flags(flags, deg);
  fprintf(stderr, "Time for intialization = %f\n", realtime()-t);
  mpz_random_dense_poly(pol1, deg, nbits);
  //  USOLVEmpz_poly_print_maple(pol1, deg);

  for(unsigned long int i = 0; i <= deg; i++){
    mpz_set(pol2[i], pol1[i]);
    mpz_set(pol3[i], pol1[i]);
  }
  int nloops = atoi(argv[4]);
  fprintf(stderr, "nloops = %d\n\n", nloops);

  double e_time = realtime();

  for(unsigned long int i = 0; i <= deg; i++){
    mpz_set(pol1[i], pol2[i]);
  }

  e_time =  realtime();
  for(int i = 0; i < nloops; i++)
    taylorshift1_dac(pol1, deg, tmpol, flags->shift_pwx, flags->pwx, nthreads);
  fprintf(stdout, "Elapsed time %f (Taylor shift Usolve nthreads = %d)\n", realtime() - e_time, nthreads);

  mpz_t **tmp = malloc(sizeof(mpz_t *) * (nthreads + 1));
  mpz_t **pols = malloc(sizeof(mpz_t *) * nthreads);
  fprintf(stderr, "Allocation starts\n");
  for(int i = 0; i <= nthreads; i++){
    tmp[i] = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
    for(unsigned long int j = 0; j <= deg; j++){
      mpz_init2(tmp[i][j], 2 * deg);
    }
  }
  for(int i = 0; i < nthreads; i++){
    pols[i] = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));
    for(unsigned long int j = 0; j <= deg; j++){
      mpz_init2(pols[i][j], 2 * deg);
    }
  }
  fprintf(stderr, "Allocation done\n\n");

  e_time = realtime();
  for(int i=0; i < nloops; i++)
    tbr_taylor_shift(pol2, deg, tmp, pols, flags->shift_pwx, flags->pwx, nthreads);
  fprintf(stdout, "Elapsed time (New TBR_Taylor shift with nthreads = %d) %f\n\n", nthreads, realtime() - e_time);


  e_time = realtime();
  for(int i=0; i < nloops; i++)
    taylor_shift(pol3, deg, tmp, pols, flags->shift_pwx, flags->pwx, nthreads);
  fprintf(stdout, "Elapsed time (New Taylor shift not in Usolve yet with nthreads = %d) %f\n\n", nthreads, realtime() - e_time);

  for(int i = 0; i < nthreads; i++){
    for(unsigned long int j = 0; j <= deg; j++){
      mpz_clear(tmp[i][j]);
      mpz_clear(pols[i][j]);
    }
    free(tmp[i]);
    free(pols[i]);
  }
  free(tmp);
  free(pols);

  for(unsigned long i=0; i <= deg; i++){
    if(mpz_cmp(pol1[i], pol2[i])!=0){
      fprintf(stderr, "BUG in TBR\n\n");
      return 1;
    }
    if(mpz_cmp(pol1[i], pol3[i])!=0){
      fprintf(stderr, "BUG in TL\n\n");
      return 1;
    }
  }

  for(unsigned long int i = 0; i <= deg; i++){
    mpz_clear(pol1[i]);
    mpz_clear(pol2[i]);
    mpz_clear(pol3[i]);
    mpz_clear(tmpol[i]);
  }
  free_heap_flags(flags, deg);

  free(pol1);
  free(pol2);
  free(pol3);
  free(tmpol);
  return 1;
}



