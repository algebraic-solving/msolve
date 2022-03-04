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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>

#define bit_one_index(x) mpz_scan1((x), 0)

#ifdef USEFLINT
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_poly.h"
#endif

#ifndef USOLVE
#define USOLVE
#endif

unsigned long int mpz_poly_max_bsize_coeffs(mpz_t *upol, long int deg){
  if(deg<0) return -1;
  unsigned long int max = 0, bs;
  for(int i=0 ; i<=deg; i++){
    bs = ilog2_mpz(upol[i]);
    if(bs>max){
      max = bs;
    }
  }
  return max;
}

/* returns c s.t. 2^c divides the content of upol */
static unsigned long int mpz_poly_remove_binary_content(mpz_t *upol, unsigned long deg)
{
  unsigned long int c, i, nbit;

  i = 0;
  if(mpz_sgn(upol[deg])==0){
    return 0;
  }
  while (mpz_sgn(upol[i]) == 0)
    i++;

  c = bit_one_index(upol[i]);

  for( ; (i <= deg) && c; i++){
    if (mpz_sgn(upol[i]) != 0)
      {
        nbit = bit_one_index(upol[i]);
        if (nbit < c)
          c = nbit;
      }
  }
  if (c == 0) return 0;

  for (i = 0; i <= deg; i++){
    mpz_fdiv_q_2exp(upol[i], upol[i], c);
  }

  return c;
}

static inline void USOLVEmpz_get_content(mpz_t *upol, unsigned long int deg, mpz_t *g){
  if(deg>=1){
    mpz_set(g[0], upol[0]);
    for(unsigned long int i = 1; i <= deg; i++){
      mpz_gcd(g[0], g[0], upol[i]);
      if(mpz_cmp_ui(g[0], 1)==0){
        return;
      }
    }
  }
}

static inline void USOLVEmpz_remove_content(mpz_t *upol, unsigned long int deg, mpz_t *g){
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_divexact(upol[i], upol[i], g[0]);
  }
}

static inline void USOLVEmpz_restore_content(mpz_t *upol, unsigned long int deg, mpz_t *g){
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_mul(upol[i], upol[i], g[0]);
  }
}


/* computed the numerator of the quotient of the division of upol by (x-c/2^k) */
static void USOLVEnumer_quotient(mpz_t *upol, unsigned long int *deg, mpz_t c, unsigned long int k){

  for(unsigned long int i = 0; i <= *deg; i++){
    mpz_mul_2exp(upol[i], upol[i], (*deg - 1) * k );
  }
  mpz_t tmp;
  mpz_init(tmp);

  for(long int i = *deg-1; i >=1; i--){
    mpz_div_2exp(tmp, upol[i+1], k);
    mpz_mul(tmp, tmp, c);
    mpz_add(upol[i], upol[i], tmp);
  }

  for(long int i = 0 ; i <= *deg - 1; i++){
    mpz_set(upol[i], upol[i+1]);
  }
  mpz_poly_remove_binary_content(upol, *deg - 1);
  *deg = *deg - 1;
  mpz_clear(tmp);
}

/* prints coefficients of pol in increasing degree order */
static inline void USOLVEmpz_poly_print(mpz_t *upol, unsigned long int deg){
  int  i;
  for(i=0;i<=deg;i++){
    gmp_printf("%Zd", upol[i]);
    fprintf(stdout," ");
  }
  fprintf(stdout, "\n");
}

/* prints coefficients of pol in decreasing degree order */
static inline void USOLVEmpz_poly_print_maple(mpz_t *upol, unsigned long int deg){
  int  i;
  for(i=deg;i>=1;i--){
    gmp_printf("(%Zd)", upol[i]);
    fprintf(stdout,"*x^%d+",i);
  }
  gmp_printf("(%Zd)", upol[0]);
  fprintf(stdout, ";\n");
}

/* From the polynomial upol, of degree deg, and b>=0, computes (inplace)
   - when b>0 upoly(X*2^b), remove its content c, and return the content c.
   - when b<=0 numer(upoly(X*2^b)), remove its content c, and
   return the content c.
*/
static inline int USOLVEmpz_poly_rescale_normalize_2exp_th(mpz_t *upol, long int b, unsigned long deg,
                                                           unsigned int nthreads){
  long int i;
  if (b > 0) {
    //    j = b;
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads)
    for(i = 1; i <= deg; i++){
      mpz_mul_2exp(upol[i], upol[i], i*b);
    }
  }
  else{
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads)
      for(i=0; i<deg; i++){
        mpz_mul_2exp(upol[i], upol[i], (i-deg)*b);
    }
  }
  return mpz_poly_remove_binary_content(upol, deg);
}

/* From the polynomial upoly, of degree deg, and b>=0, computes (inplace)
   - when b>0 upoly(X*c*2^b), remove its content, and return the content.
   - when b<=0 numer(upoly(X*c*2^b)), remove its content, and
   return the content.
*/
static inline int USOLVEmpz_poly_rescale_normalize_2exp_th_long(mpz_t *upol,
                                                                unsigned long deg,
                                                                long int b,
                                                                long int c,
                                                                unsigned int nthreads){
  long int i;
  mpz_t coef;mpz_init(coef);mpz_set_si(coef, c);
  if (b > 0) {
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads)
    for(i = 1; i <= deg; i++){
      mpz_mul(upol[i], upol[i], coef);
      mpz_mul_2exp(upol[i], upol[i], i*b);
      mpz_mul_si(coef, coef, c);
    }
  }
  else{
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads)
    for(i=0; i<deg; i++){
      mpz_mul(upol[i], upol[i], coef);
      mpz_mul_2exp(upol[i], upol[i], (i-deg)*b);
      mpz_mul_si(coef, coef, c);
    }
  }
  mpz_clear(coef);
  return mpz_poly_remove_binary_content(upol, deg);
}

unsigned long int mpz_poly_min_bsize_coeffs(mpz_t *upol,
                                            long int deg){
  if(deg < 0) return 1;
  unsigned long int min = ilog2_mpz(upol[deg]), bs;
  for(long i = deg ; i >= 0; i--){
    bs = ilog2_mpz(upol[i]);
    if(bs<min && mpz_cmp_ui(upol[i], 0) != 0){
      min = bs;
    }
  }
  return min;
}


/* One assumes up(0) != 0 */
/* as soon as more that 3 sign variations are found,
   the computation is stopped */
/* takes into accound the loss of precison */
static long USOLVEmpz_poly_sgn_variations_coeffs_bsize(mpz_t* upol,
                                                       unsigned long deg,
                                                       unsigned long int bsize){
  unsigned long int i;
  long nb = 0;
  int s = mpz_sgn(upol[deg]);
  if(s==0){
    return -1;
  }
  int boo = 1;
  unsigned long int N = 1;
  int L = LOG2(bsize);
  for(i = deg-1 ; i > 0; i--){
    int c = mpz_cmp_ui(upol[i],0);
    long int l = ilog2_mpz(upol[i]);
    N = min((bsize - i + 1) * L + 1, min(deg, (i + 1)*L + 1));
    if(l <= N || c==0){
      boo = 0;
    }
    if( mpz_sgn(upol[i]) * s < 0 && l > N && c!=0){
      nb = nb + 1;
      s = mpz_sgn(upol[i]);
      if(nb >= 3){
        return nb;
      }
    }
  }
  N = L + 1;
  int c = mpz_cmp_ui(upol[0],0);
  long int l = ilog2_mpz(upol[0]);
  if(l <= N || c==0){
    boo = 0;
  }
  if(s*mpz_sgn(upol[0]) < 0 && l> N && c!=0){
    nb = nb + 1;
  }
  if(nb>=3 || boo == 1){
    return nb;
  }
  return -1;
}


