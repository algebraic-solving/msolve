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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include<stdio.h>
#include<stdlib.h>
#include<gmp.h>
#ifdef _OPENMP
#include<omp.h>
#endif

/* for timing functions */
#include "../neogb/tools.h"

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"

#include "../msolve/msolve-data.h"

#define THRESHOLDSHIFT 256
#define POWER_HACK 1

#define ilog2(a) mpz_sizeinbase(a,2)

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define ilog2_mpz(a) mpz_sizeinbase(a,2)

#define min(a, b) (((a) < (b)) ? (a) : (b))
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#include "data_usolve.c"
#include "utils.c"
#include "taylor_shift.c"
#include "descartes.c"
#include "evaluate.c"
#include "print_usolve.c"
#include "refine.c"

static long bisection_rec(mpz_t *, unsigned long *,
                          mpz_t, long,
                          interval *,
                          unsigned long int *,
                          usolve_flags *,
                          mpz_t);

/* Kioustelidis's bound : 2*max(abs(-ai/an)^(1/(n-i)),
   the maximum being taken over the i with sgn(a_i) != sgn(an) */

/* Computes M such that 2^M is greater than Kioustelidis's bound */
/* which is 2*max(abs(-ci/cn)^(1/(n-i)) for those i's such that sign(ci) !=
   sign(c[deg])*/
/* (the ci's are the coefficients of upol) */
/* returns -1 if there is obviously no positive root */
static long bound_roots(mpz_t *upol, const unsigned long deg){

  int b = 1;

  /* 2^l < c[deg] < 2^(l+1)*/
  const long l = ilog2(upol[deg]) - 1;
  long M = -l;

  for(unsigned long int i = 0; i < deg; i++){
    long pow, pow2;
    if (mpz_sgn(upol[deg]) != mpz_sgn(upol[i])){
      b = 0;
      pow = ilog2(upol[i]) ;
      pow -= l; /* 2^currpow >= abs(-ai/an) */

      if(pow > 0){
        pow2 = pow / (deg - i);
      }
      else{
        pow2 = -((-pow) / (deg - i));
      }

      if(pow2 * ((long) (deg - i)) != pow){
        pow2++;
      }

      if(pow2 > M){
        M = pow2;
      }
    }
  }
  /* there can't be positive roots */
  if (b == 1) return -1;

  M++;
  return M;
}

static inline void compute_shift_pwx(mpz_t **shifted,
                                     const unsigned long int deg,
                                     const unsigned long int npwr,
                                     unsigned long int pwx,
                                     const unsigned int nthreads){
  taylorshift1_naive(shifted[0], pwx);
  for(unsigned long int i=1; i < npwr; i++){
    mpz_poly_mul(shifted[i],
                 shifted[i-1], pwx, shifted[i-1], pwx,
                 nthreads);
    pwx = 2 * pwx;
  }
}

static unsigned long int compute_degpower(const unsigned long int deg){
  unsigned long int pwr=deg;

  while(pwr-1 >= THRESHOLDSHIFT){
    pwr = pwr / 2;
  }
  return pwr;
}


static void assign_sign_left(interval *roots,
                             const unsigned long int nbneg,
                             const unsigned long int nbpos,
                             const int b){
  int s = b;
  for(unsigned long int i = 0; i < nbneg; i++){
    (roots+i)->sign_left = s;
    s=-s;
  }
  for(unsigned long int i = nbneg; i < nbpos; i++){
    (roots+i)->sign_left = s;
    s=-s;
  }
}

static void merge_root(interval *roots,
                       mpz_t c, long k,
                       const unsigned int isexact,
                       const int sgnlft,
                       const unsigned long int nb,
                       const int b_pos,
                       const int b_neg,
                       const int sgn){
  int b;
  if(sgn >0){
    b = b_neg;
  }
  else{
    b = b_pos;
  }
  mpz_init(roots[nb].numer);

  if (k <= b)
    {
      if (sgn){
          mpz_neg(roots[nb].numer, c);
          if (!isexact){
            mpz_sub_ui(roots[nb].numer, roots[nb].numer, 1);
          }
          mpz_mul_2exp(roots[nb].numer, roots[nb].numer, b-k);
        }
      else{
        mpz_mul_2exp(roots[nb].numer, c, b-k);
      }

      roots[nb].isexact = isexact;
      if(roots[nb].isexact==1){
        roots[nb].k = 0;
      }
      else{
        roots[nb].k = k - b;
      }
      roots[nb].sign_left = sgnlft;

      return;
    }
  else
    {
      if (sgn){
        mpz_neg(roots[nb].numer, c);
        if (!isexact)
          mpz_sub_ui(roots[nb].numer, roots[nb].numer, 1);
      }
      else{
        mpz_set(roots[nb].numer, c);
      }
      roots[nb].k = k - b;
      roots[nb].isexact = isexact;
      roots[nb].sign_left = sgnlft;
    }
  return;
}


static inline void allocate_shift_pwx(mpz_t **shifted,
                                      const unsigned long int npwr,
                                      const unsigned long int pwx){

  unsigned long int newpwx=pwx;
  mp_bitcnt_t nbits=LOG2(newpwx);

  for(int i = 0; i < npwr; i++){
    shifted[i] = (mpz_t *)malloc(sizeof(mpz_t)*(newpwx+1));
    for(int j = 0; j <= newpwx; j++){
      mpz_init2(shifted[i][j], nbits);
    }
    newpwx = 2 * newpwx;
    nbits = LOG2(newpwx);
  }

  mpz_set_ui(shifted[0][pwx],1);

}

static inline void unallocate_shift_pwx(mpz_t **shifted,
                                        const unsigned long int npwr,
                                        unsigned long int pwx){
  for(int i = 0; i < npwr; i++){
    for(int j = 0;j <= pwx; j++){
      mpz_clear(shifted[i][j]);
    }
    pwx=2*pwx;
  }
  for(int i = 0; i < npwr; i++){
    free(shifted[i]);
  }
}

static inline long mpz_poly_sign_at_one(mpz_t *upol,
                                        const unsigned long int deg,
                                        mpz_t *vals){
  mpz_set_ui(vals[0], 0);
  for(unsigned long int i=0; i <= deg; i++){
    mpz_add(vals[0], vals[0],upol[i]);
  }
  return mpz_sgn(vals[0]);
}


/* assumes upol[0 != 0] */
static inline long mpz_poly_sgn_variations_coeffs_is_zero(mpz_t* upol,
                                                          const unsigned long deg){
  unsigned long int i;
  long nb = 0;
  int s = mpz_sgn(upol[deg]);

  for(i = deg - 1; i > 0; i--){
    if( mpz_sgn(upol[i]) * s < 0 ){
      nb = nb + 1;
      s = mpz_sgn(upol[i]);
      if(nb >= 1){
        return nb;
      }
    }
  }
  if(s * mpz_sgn(upol[0]) < 0){
    nb = nb + 1;
  }
  return nb;
}


/* if 0 is a root, adds c / 2^k to the array of roots and update upol */
static void is_zero_root(mpz_t *upol, mpz_t c, long k,
                         unsigned long *deg,
                         interval *roots,
                         unsigned long int *nbroot,
                         usolve_flags *flags){
  unsigned long int i, j;

  if(mpz_cmp_ui(upol[0], 0) == 0){
    i = 1;
    while(mpz_cmp_ui(upol[i], 0) == 0){
      i++;
    }

    if(i >= 2){
      fprintf(stderr, "error: the polynomial is not square-free\n");
      exit(1);
    }

    for(j = 0; j < i; j++){
      merge_root(roots, c, k, 1, 0, *nbroot,
                 flags->bound_pos, flags->bound_neg,
                 flags->sign);
      (*nbroot) ++;
    }

    *deg -= i;
    for (j = 0; j <= *deg; j++, i++){
      mpz_set(upol[j], upol[i]);
    }

    unsigned long int newpwx = compute_degpower(flags->cur_deg);

    if(newpwx != flags->pwx){
      unallocate_shift_pwx(flags->shift_pwx,
                           flags->npwr,
                           flags->pwx);
      free(flags->shift_pwx);

      flags->pwx =  newpwx;
      if(newpwx >= flags->cur_deg){
        flags->nblocks = 0;
        flags->npwr= 0;
      }
      else{
        flags->nblocks = 1<<LOG2(flags->cur_deg/flags->pwx);
        flags->npwr=LOG2(flags->nblocks);
      }
      if(flags->npwr>0){
        flags->shift_pwx=(mpz_t **)malloc(sizeof(mpz_t *)*(flags->npwr));

        allocate_shift_pwx(flags->shift_pwx,
                           flags->npwr, flags->pwx);

        compute_shift_pwx(flags->shift_pwx,
                          flags->cur_deg,
                          flags->npwr,
                          flags->pwx,
                          flags->nthreads);
      }
      else{
        flags->shift_pwx = NULL;
      }
    }
  }
}

/* it returns the sign of upoly(1/2) */
static int sgn_mpz_upoly_eval_onehalf(mpz_t *upol,
                                      unsigned long int deg,
                                      usolve_flags *flags){
  long i;

  mpz_set(flags->Values[0], upol[deg]);

  for (i = deg - 1; i >= 0; i--){
      mpz_mul_2exp(flags->Values[1], upol[i], deg - i);
      mpz_add(flags->Values[0], flags->Values[0], flags->Values[1]);
  }

  i = mpz_sgn(flags->Values[0]);
  return i;
}


static void numer_quotient_uint(mpz_t *upol, unsigned long int *deg,
                                unsigned int c, long k){

  for(unsigned long int i = 0; i <= *deg; i++){
    mpz_mul_2exp(upol[i], upol[i], (*deg - 1) * k );
  }
  mpz_t tmp;
  mpz_init(tmp);

  for(long int i = *deg - 1; i >=1; i--){
    mpz_div_2exp(tmp, upol[i+1], k);
    mpz_mul_ui(tmp, tmp, c);
    mpz_add(upol[i], upol[i], tmp);
  }

  for(long int i = 0 ; i <= *deg - 1; i++){
    mpz_set(upol[i], upol[i+1]);
  }

  mpz_poly_remove_binary_content(upol, *deg-1);
  *deg = *deg - 1;
  mpz_clear(tmp);
}


static long int manage_root_at_one_half(mpz_t *upol,
                                        mpz_t c, long k,
                                        unsigned long *deg,
                                        mpz_t tmp,
                                        unsigned long int *nbr,
                                        usolve_flags *flags){
  long int sh = sgn_mpz_upoly_eval_onehalf(upol, *deg, flags);

  while (sh == 0) {
    mpz_set(tmp, c);
    mpz_mul_2exp(tmp, tmp, 1);
    mpz_add_ui(tmp, tmp, 1);

    numer_quotient_uint(upol, deg, 1, 1);

    sh = sgn_mpz_upoly_eval_onehalf(upol, *deg, flags);

    flags->cur_deg = (*deg);
    if(sh == 0){
      fprintf(stderr, "Input polynomial is not square-free\n");
      exit(1);
    }

    unsigned long int newpwx = compute_degpower(flags->cur_deg);

    if(flags->classical_algo == 0 && newpwx !=flags->pwx){
      unallocate_shift_pwx(flags->shift_pwx,
                           flags->npwr,
                           flags->pwx);
      free(flags->shift_pwx);

      flags->pwx =  newpwx;
      if(newpwx>=flags->cur_deg){
        flags->nblocks = 0;
        flags->npwr= 0;
      }
      else{
        flags->nblocks = 1<<LOG2(flags->cur_deg/flags->pwx);
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
    }
  }
  return sh;
}


static long nb_0_case_in_bisection_rec(const mpz_t c, const long k,
                                       mpz_t tmp,
                                       interval *roots, unsigned long int *nbr,
                                       usolve_flags *flags, const int is_half_root){
  if(is_half_root){
    mpz_set(tmp, c);
    mpz_mul_2exp(tmp, tmp, 1);
    mpz_add_ui(tmp, tmp, 1);
    merge_root(roots, tmp, k + 1, 1, 0, *nbr,
               flags->bound_pos, flags->bound_neg, flags->sign);
    (*nbr)++;

    if(flags->hasrealroots==1 && (*nbr)>0){
      mpz_clear(tmp);
      return -1;
    }

    if(flags->verbose >= 1){
      fprintf(stderr,"+");
      if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
    }
  }
  if(flags->verbose >= 1){
    fprintf(stderr,"!");
  }
  mpz_clear(tmp);
  return k;
}


static long nb_1_case_in_bisection_rec(mpz_t c, const long k, 
                                       mpz_t tmp,
                                       interval *roots, unsigned long int *nbr,
                                       usolve_flags *flags, const int is_half_root,
                                       const long bot, const long top,
                                       const long half){
    if(is_half_root){
      if(bot * half <0){
        mpz_set(tmp, c);
        mpz_mul_2exp(tmp, tmp, 1);
        merge_root(roots, tmp, k+1, 0, bot, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr) ++;

        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots == 1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
        /* adds 1 / 2 */
        mpz_add_ui(tmp, tmp, 1);
        merge_root(roots, tmp, k+1, 1, 0, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr)++;

        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
      }
      else{
        /* adds 1/ 2 */
        mpz_set(tmp, c);
        mpz_mul_2exp(tmp, tmp, 1);
        mpz_add_ui(tmp, tmp, 1);
        merge_root(roots, tmp, k+1, 1, 0, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr)++;
        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots==1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
        if(top != 0){
          /* 1 is not a root */
          merge_root(roots, tmp, k+1, 0, half, *nbr,
                     flags->bound_pos, flags->bound_neg, flags->sign);
          (*nbr)++;
          if(flags->verbose >= 1){
            fprintf(stderr,"+");
            if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
          }
        }
        else{
          /* 1 is a root */
          mpz_add_ui(tmp, tmp, 1);
          merge_root(roots, tmp, k+1, 1, 0, *nbr,
                     flags->bound_pos, flags->bound_neg, flags->sign);
          (*nbr) ++;
          if(flags->verbose>=1){
            fprintf(stderr,"+");
            if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
          }
          if(flags->hasrealroots == 1 && (*nbr)>0){
            mpz_clear(tmp);
            return -1;
          }
        }
      }
    }
    else{
      if(top != 0){
        merge_root(roots, c, k, 0, bot, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr)++;
        if(flags->verbose>=1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots==1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
      }
      else{
        mpz_set(tmp, c);
        mpz_add_ui(tmp, tmp, 1);
        merge_root(roots, tmp, k, 1, 0, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr)++;
        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots == 1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
      }
    }
    mpz_clear(tmp);
	  return k;
}


static long nb_2_case_in_bisection_rec(mpz_t c, const long k,
                                       mpz_t tmp,
                                       interval *roots, unsigned long int *nbr,
                                       usolve_flags *flags, const int is_half_root,
                                       const long bot, const long top,
                                       const long shalf){
      //attention t peut etre egale a 0 mais 1/2 ne peut etre racine du polynome
      //courant (voir ci-dessus)
      (flags->half_done)++;
      mpz_set(tmp, c);

      mpz_mul_2exp(tmp, tmp, 1);
      merge_root(roots, tmp, k+1, 0, bot, *nbr,
                 flags->bound_pos, flags->bound_neg, flags->sign);
      (*nbr)++;

      if(flags->verbose>=1){
        fprintf(stderr,"+");
        if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
      }
      if(flags->hasrealroots==1 && (*nbr)>0){
        mpz_clear(tmp);
        return -1;
      }

      mpz_add_ui(tmp, tmp, 1);
      if(is_half_root){
        merge_root(roots, tmp, k+1, 1, 0, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr)++;
        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots==1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
      }
      if(top != 0){
        merge_root(roots, tmp, k+1, 0, shalf, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr) ++;
        if(flags->verbose >= 1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots == 1 && (*nbr) > 0){
          mpz_clear(tmp);
          return -1;
        }
      }
      else{//1 est racine
        //BUG Ici il faut diviser P par x-1

        mpz_add_ui(tmp, tmp, 1);
        merge_root(roots, tmp, k+1, 1, 0, *nbr,
                   flags->bound_pos, flags->bound_neg, flags->sign);
        (*nbr) ++;
        if(flags->verbose>=1){
          fprintf(stderr,"+");
          if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
        }
        if(flags->hasrealroots==1 && (*nbr)>0){
          mpz_clear(tmp);
          return -1;
        }
      }
      /* if(flags->debug==1){ */
      /*   long bo=mpz_sgn(upoly[0]); */
      /*   long to=mpz_poly_sign_at_one(upol, (*deg), flags->Values); */
      /*   fprintf(stderr,"[d:%ld, %ld, %ld]",bo,shalf,to); */
      /* } */
      mpz_clear(tmp);
      return k;
}


static long nb_default_case_in_bisection_rec(mpz_t *upol, unsigned long int *deg,
                                             mpz_t c,
                                             const long k, mpz_t tmp,
                                             interval *roots,
                                             unsigned long int *nbr,
                                             usolve_flags *flags, int is_half_root,
                                             mpz_t tmp_half){
  long oldk, nb;
  double e_time;
  mpz_set(tmp, c);

  int branch_left = 1;
  int branch_right = 1;

  if(flags->verbose>=1){
    fprintf(stderr,"-");
  }
  mpz_mul_2exp(tmp, tmp, 1);
  USOLVEmpz_poly_rescale_normalize_2exp_th(upol, -1,
                                           *deg, flags->nthreads);

  if(branch_left==1){
    oldk = bisection_rec(upol, deg, tmp, k+1, roots, nbr,
                         flags, tmp_half);
  }
  else{
    if(flags->verbose >= 1) fprintf(stderr, "!");
    oldk = k + 1;
  }

  if(flags->hasrealroots == 1 && (*nbr)>0){
    mpz_clear(tmp);
    return -1;
  }
  mpz_add_ui(tmp, tmp, 1);

  if(is_half_root){
    merge_root(roots, tmp, k + 1, 1, 0, *nbr,
               flags->bound_pos, flags->bound_neg,
               flags->sign);
    (*nbr)++;
    if(flags->verbose>=1){
      fprintf(stderr,"+");
      if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
    }
    if(flags->hasrealroots == 1 && (*nbr)>0){
      mpz_clear(tmp);
      return -1;
    }
  }

  if (oldk == -1){
    mpz_clear(tmp);
    return -1;
  }

  e_time =  realtime();
  if(flags->verbose>=1){
    fprintf(stderr,"*");
  }
  if(flags->classical_algo==1){
    taylorshift1_naive(upol, *deg);
  }
  else{
    taylorshift1_dac(upol, *deg,
                     flags->tmpol, flags->shift_pwx, flags->pwx,
                     flags->nthreads);
  }
  flags->time_shift += (realtime()-e_time);
  (flags->transl)++;

  if (oldk > k + 1){
    USOLVEmpz_poly_rescale_normalize_2exp_th(upol, oldk - (k + 1), *deg,
                                             flags->nthreads);
    }
  if(branch_right==1){
    oldk = bisection_rec(upol, deg, tmp, k+1, roots, nbr,
                         flags, tmp_half);
  }
  else{
    oldk = k + 1;
    nb = mpz_poly_sgn_variations_coeffs_is_zero(upol, (*deg));
    if(flags->verbose >= 1) fprintf(stderr, "!");
    if(nb == 0) oldk = -1;
  }

  if (oldk == -1){
    mpz_clear(tmp);
    return -1;
  }

  if(flags->hasrealroots==1 && (*nbr)>0){
    mpz_clear(tmp);
    return -1;
  }
  mpz_clear(tmp);
  return oldk;
}



/* investigates  (c / 2^k, (c+1) / 2^k) */
/* returns k */
static long bisection_rec(mpz_t *upol, unsigned long *deg,
                          mpz_t c, long k,
                          interval *roots,
                          unsigned long int *nbr,
                          usolve_flags *flags,
                          mpz_t tmp_half){

  /* display_roots_system(stderr, roots, *nbr); */
  long nb;
  long shalf;
  long bsgn = 0;
  mpz_t tmp;
  double e_time;
  mpz_init(tmp);

  (flags->node_looked)++;

  if(flags->verbose == 4){
    fprintf(stderr,"[");
    mpz_out_str(stderr, 10, c);
    fprintf(stderr,",%lu]", k);
  }
  if(flags->verbose >= 5){
    fprintf(stderr,"[");
    mpz_out_str(stderr, 10, c);
    fprintf(stderr,",%lu][bs=%lu]", k,
            mpz_poly_max_bsize_coeffs(upol, *deg));
  }

  is_zero_root(upol, c, k, deg, roots, nbr,
               flags);
  unsigned int olddeg = *deg;

  if(flags->hasrealroots == 1 && (*nbr)>0){
    return -1;
  }


  /* if 1 / 2 is a root, upol is divided by (x-1/2) */
  /* shalf takes the sign of the new obtained poly at 1/2 */
  shalf = manage_root_at_one_half(upol, c, k, deg,
                                  tmp_half,
                                  nbr,
                                  flags);
  if(flags->hasrealroots == 1 && (*nbr)>0){
    mpz_clear(tmp);
    return -1;
  }

  int is_half_root = 0;

  if(olddeg != *deg){
    is_half_root = 1;
  }


  nb = mpz_poly_sgn_variations_coeffs_is_zero(upol, (*deg));

  if(nb == 0){

    if(is_half_root){
      mpz_set(tmp, c);
      mpz_mul_2exp(tmp, tmp, 1);
      mpz_add_ui(tmp, tmp, 1);

      merge_root(roots, tmp, k+1, 1, 0, *nbr,
                 flags->bound_pos, flags->bound_neg, flags->sign);
      (*nbr)++;

      if(flags->verbose>=1){
        fprintf(stderr,"+");
        if((*nbr) % 100 == 0) fprintf(stderr, "[%lu]",(*nbr));
      }

      if(flags->hasrealroots==1 && (*nbr)>0){
        mpz_clear(tmp);
        return -1;
      }
    }

    if(flags->verbose>=1){
      fprintf(stderr,"!");
    }
    mpz_clear(tmp);
    return -1;
  }


  if(flags->verbose>=2){
    fprintf(stderr,"c");
  }

  e_time = realtime();
  if(flags->classical_algo==1){
    nb = descartes_classical(upol, flags->tmpol_desc, *deg, shalf, &bsgn);
  }
  else{
    nb = descartes(upol, flags->tmpol_desc, *deg, shalf, &bsgn, flags);
  }
  flags->time_desc += (realtime()-e_time);

  if(flags->verbose>=3){
    fprintf(stderr,"[nb=%lu]",nb);
  }

  if (bsgn == -1){
	  mpz_clear(tmp);
	  return -1;
	}

  long bot = mpz_sgn(upol[0]);
  long top = mpz_poly_sign_at_one(upol, (*deg), flags->Values);

  switch (nb){

  case 0:

    /* no root */
    return nb_0_case_in_bisection_rec(c, k, tmp,
                                      roots, nbr,
                                      flags, is_half_root);

	case 1:

    /* one root exactly */

    return nb_1_case_in_bisection_rec(c, k, tmp,
                                      roots, nbr,
                                      flags, is_half_root,
                                      bot, top, shalf);


	case 2:
	  if ((bsgn && (flags->classical_algo)==1) ||
        (flags->classical_algo == 0 && bot * shalf<0 && top * shalf<=0)){
      return nb_2_case_in_bisection_rec(c, k, tmp,
                                        roots, nbr,
                                        flags, is_half_root,
                                        bot, top, shalf);

    }
	default: /* recursive call on each half of the interval */
    return nb_default_case_in_bisection_rec(upol, deg, c, k, tmp,
                                            roots, nbr,
                                            flags, is_half_root,
                                            tmp_half);
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
  flags->Values = NULL;
  flags->tmp_threads = NULL;
  flags->pols_threads = NULL;

  flags->time_desc = 0;
  flags->time_shift = 0;
  flags->nthreads = 1;
  flags->verbose = 0;
  flags->bfile = 0;
  flags->classical_algo = 0;
  flags->print_stats = 0;
  flags->debug = 0;
}


static void initialize_heap_flags(usolve_flags *flags,
                                  const unsigned long int deg){

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

    if(flags->nthreads>=1 && 0==1){
      /* used for parallel Taylor shift */
      flags->tmp_threads = malloc(sizeof(mpz_t *) * (flags->nthreads + 1));
      flags->pols_threads = malloc(sizeof(mpz_t *) * flags->nthreads);
      mpz_t ** tmp = flags->tmp_threads;
      mpz_t ** pols = flags->pols_threads;

      for(int i = 0; i <= flags->nthreads; i++){
        tmp[i] = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));

        for(unsigned long int j = 0; j <= deg; j++){
          mpz_init2(tmp[i][j], 2 * deg);
        }

      }

      for(int i = 0; i < flags->nthreads; i++){
        pols[i] = (mpz_t *)malloc((deg + 1) * sizeof(mpz_t));

        for(unsigned long int j = 0; j <= deg; j++){
          mpz_init2(pols[i][j], 2 * deg);
        }

      }
    }
  }

  flags->Values = malloc(sizeof(mpz_t) * 2);
  mpz_init(flags->Values[0]);
  mpz_init(flags->Values[1]);
}


/* warning: does not free flags->shift_pwx */

static inline void free_heap_flags(usolve_flags *flags,
                                   const unsigned long int deg){

  if(flags->classical_algo == 0){

    for(unsigned long int i = 0; i <= deg; i++){
      mpz_clear(flags->tmpol_desc[i]);
      mpz_clear(flags->tmpol[i]);
    }
    mpz_clear(flags->Values[0]);
    mpz_clear(flags->Values[1]);
    free(flags->Values);

    if(flags->nthreads>=1 &&0==1){
      for(int i = 0; i < flags->nthreads; i++){

        for(unsigned long int j = 0; j <= deg; j++){
          mpz_clear(flags->tmp_threads[i][j]);
          mpz_clear(flags->pols_threads[i][j]);
        }

        free(flags->tmp_threads[i]);
        free(flags->pols_threads[i]);
      }
      free(flags->tmp_threads);
      free(flags->pols_threads);
    }
  }
}

static void display_stats(usolve_flags *flags){

  fprintf(stderr,"\n");
  fprintf(stderr,"Number of nodes : %lu\n", flags->node_looked);
  fprintf(stderr,"Number of shifts : %lu\n", flags->transl);
  fprintf(stderr,"Number of half splits : %lu\n", flags->half_done);
  fprintf(stderr,"Time in Descartes (elapsed): %.2f sec\n", flags->time_desc);
  fprintf(stderr,"Time in Taylor shifts (elapsed): %.2f sec\n", flags->time_shift);
  fprintf(stderr,"\n");

}



/* warning: does not check that upol is square-free (should be done outside) */

interval *bisection_Uspensky(mpz_t *upol0, unsigned long deg,
                             unsigned long int *nb_pos_roots,
                             unsigned long int *nb_neg_roots,
                             usolve_flags *flags){

  interval *pos_roots = (interval *)malloc(deg * sizeof(interval));
  interval *neg_roots = (interval *)malloc(deg * sizeof(interval));
  unsigned int nb_positive_roots = 0, nb_negative_roots = 0;
  mpz_t tmp_half;
  mpz_init(tmp_half);

  unsigned long deg0 = deg;

  mpz_t e;
  mpz_init_set_ui(e, 0);

  *nb_pos_roots = 0;
  *nb_neg_roots = 0;

  int zero_root = 0;

  if(mpz_sgn(upol0[0]) == 0){
    merge_root(pos_roots, e, 0, 1, 0, *nb_pos_roots,
               0, 0, 1);
    (*nb_pos_roots) ++;
    nb_positive_roots++;
    zero_root = 1;
  }

  if(mpz_sgn(upol0[0]) == 0 && mpz_sgn(upol0[1])==0){

    fprintf(stderr,
            "0 is a multiple root ; the input polynomial must be square-free\n");

    free(pos_roots);
    free(neg_roots);
    mpz_clear(e);
    mpz_clear(tmp_half);
    exit(1);

  }

  if(zero_root==0){
    deg = deg0;
  }
  else{
    deg = (deg0 - 1);
  }
  mpz_t *upol = (mpz_t *) malloc ((deg + 1) * sizeof(mpz_t));

  if(zero_root == 0){
    for (long i = 0; i <= deg; i++){
      mpz_init_set(upol[i], upol0[i]);
    }
  }
  else{
    for (long i = 0; i <= deg; i++){
      mpz_init_set(upol[i], upol0[1 + i]);
    }
  }

  if(flags->search >= 0 && deg > 0){

    flags->sign = 0;
    flags->bound_pos = bound_roots(upol, deg);

    if(flags->verbose>=1){
      fprintf(stderr, "Bound for positive roots: %ld\n\n", flags->bound_pos);
    }

    USOLVEmpz_poly_rescale_normalize_2exp_th(upol,
                                             flags->bound_pos, deg,
                                             flags->nthreads);
    initialize_heap_flags(flags, deg);

    unsigned olddeg = deg;
    bisection_rec(upol, &deg, e, 0,
                  pos_roots, nb_pos_roots,
                  flags, tmp_half);
    nb_positive_roots = *nb_pos_roots;

    free_heap_flags(flags, olddeg);
    unallocate_shift_pwx(flags->shift_pwx,
                         flags->npwr, flags->pwx);
  }

  /* replaces upol(x) by upol(-x) => negative roots */

  if(zero_root == 0){
    deg = deg0;
  }
  else{
    deg = (deg0 - 1);
  }
  if(zero_root==0){
    for(long i = 0; i <= deg; i++){
        if (i % 2 == 1){
          mpz_neg(upol[i], upol0[i]);
        }
        else{
          mpz_set(upol[i], upol0[i]);
        }
      }
  }
  else{
    for(long i = 0; i <= deg; i++)
      {
        if (i % 2 == 1){
          mpz_neg(upol[i], upol0[i + 1]);
        }
        else{
          mpz_set(upol[i], upol0[i + 1]);
        }
      }
  }
  if(flags->search<=0 && deg > 0){

    flags->bound_neg = bound_roots(upol, deg);

    USOLVEmpz_poly_rescale_normalize_2exp_th(upol,
                                             flags->bound_neg, deg,
                                             flags->nthreads);

    if(flags->verbose>=1){
      fprintf(stderr, "\nBound for negative roots: %ld\n\n", flags->bound_neg);
    }

    mpz_set_ui(e, 0);
    flags->sign = 1;

    initialize_heap_flags(flags, deg);
    unsigned long int olddeg = deg;
    bisection_rec(upol, &deg, e, 0,
                  neg_roots, nb_neg_roots,
                  flags, tmp_half);
    nb_negative_roots = (*nb_neg_roots);
    free_heap_flags(flags, olddeg);
    unallocate_shift_pwx(flags->shift_pwx,
                         flags->npwr, flags->pwx);
  }

  unsigned long int nbroots = nb_positive_roots + nb_negative_roots;
  interval *roots = (interval *)malloc(nbroots * sizeof(interval));
  for(unsigned int k = 0; k < nb_negative_roots; k++){
    roots[k] = neg_roots[nb_negative_roots - 1 -k];
  }
  for(unsigned int k = nb_negative_roots; k < nbroots; k++){
    roots[k] = pos_roots[k-nb_negative_roots];
  }

  /* needs to assign signs since this info may have been corrupted */
  /* (because of exact roots) */
  int b = mpz_sgn(upol0[deg0]);
  if(deg0 % 2 == 1){
    b = -b;
  }
  if(nbroots>0){
    assign_sign_left(roots, nb_negative_roots, nb_positive_roots, b);
  }

  if(zero_root == 0){
    deg = deg0;
  }
  else{
    deg = (deg0 - 1);
  }
  for(long i = 0; i <= deg; i++){
    mpz_clear(upol[i]);
  }

  free(upol);
  free(pos_roots);
  free(neg_roots);
  mpz_clear(e);
  mpz_clear(tmp_half);
  return roots;

}


interval *real_roots(mpz_t *upoly, unsigned long deg,
                     unsigned long int *nb_pos_roots,
                     unsigned long int *nb_neg_roots,
                     const int32_t precision,
                     int nthrds,
                     int info_level){
  usolve_flags *flags = (usolve_flags*)(malloc(sizeof(usolve_flags)));
  initialize_flags(flags);
  flags->cur_deg = deg;
  flags->prec_isole = precision;
  if(info_level){
    fprintf(stderr, "Real root isolation starts at precision %d\n",
            precision);
  }
  if (info_level > 0) {
    flags->verbose = info_level - 1;
  } else {
    flags->verbose = 0;
  }
  if(info_level>1){
    flags->print_stats = 1;
  }
  flags->nthreads = nthrds;
  if(flags->verbose>=1 || flags->print_stats == 1){
    fprintf(stderr, "Degree = %ld \t Max bit size = %lu Min bit size = %lu \n",
            flags->cur_deg,
            mpz_poly_max_bsize_coeffs(upoly, deg),
            mpz_poly_min_bsize_coeffs(upoly, deg));
    fprintf(stderr, "nthreads = %d\n", flags->nthreads);
  }
  double e_time = realtime ( );

  interval *roots = bisection_Uspensky(upoly, deg, nb_pos_roots, nb_neg_roots,
                                       flags);

  unsigned long int nbroots = *nb_pos_roots + *nb_neg_roots;
  for(unsigned long int i = 0; i < nbroots; i++){
    if(roots[i].isexact){
      if(roots[i].k < 0){
        roots[i].k = 0;
      }
    }
  }

  /* display_root(stderr, roots); */

  e_time = realtime ( ) - e_time;

  if(flags->verbose>=1){
    fprintf(stderr, "\n");
  }

  if((flags->verbose>=1) || (flags->print_stats>=1)){
    fprintf(stderr,"Time for isolation (elapsed): %.2f sec\n", e_time);
  }

  double refine_time = realtime();
  double step = (e_time+1) / (deg) * 1000 * LOG2(flags->prec_isole) * 2;

  if(nbroots > 0 && flags->prec_isole >= 0){
    if(flags->classical_algo > 0){
      refine_all_roots_naive(upoly,deg, roots, nbroots,
                             flags->prec_isole, flags->classical_algo, flags->debug);
    }
    else{
      refine_QIR_roots_adaptative(upoly, &deg, roots, *nb_neg_roots, *nb_pos_roots,
                                  flags->prec_isole, flags->verbose, step, flags->nthreads);
    }
  }
  refine_time = realtime() - refine_time;

  for(unsigned long int i = 0; i < nbroots; i++){
    if(roots[i].isexact){
      if(roots[i].k < 0){
        roots[i].k = 0;
      }
    }
  }

  if(flags->print_stats>=1){
    display_stats(flags);
  }

  if((flags->verbose>=1) || (flags->print_stats>=1)){
    fprintf(stderr,"Time for isolation (elapsed): %.2f sec\n", e_time);
    fprintf(stderr,"Time for refinement (elapsed): %.2f sec\n", refine_time);
  }

  free(flags);
  return roots;
}
