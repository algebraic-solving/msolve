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



#ifdef _OPENMP
#include<omp.h>
#endif

#define USEFLINT 1

#ifdef USEFLINT
  #include "flint/flint.h"
  #include "flint/fmpz.h"
//  #include "flint/fft.h"
//  #if __FLINT_VERSION < 3
//      || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR == 0)
//      || (__FLINT_VERSION == 3 && __FLINT_VERSION_MINOR == 1 && __FLINT_VERSION_PATCHLEVEL < 2)
//    #include "flint/fft_tuning.h"
//  #endif
  #include "flint/fmpz_poly.h"
#else
  #include"mpz_upoly_multiply.h"
#endif

#ifdef USEFLINT


void fmpz_poly_set_coeff_mpz2(fmpz_poly_t poly, slong n,
                             const mpz_t x){
  fmpz_t t;
  fmpz_init_set_readonly(t, x);
  fmpz_poly_set_coeff_fmpz(poly, n, t);
  fmpz_clear_readonly(t);
}

/* from mpz_t *poly to fmpz_poly_t  */

static void mpz_2_fmpz_poly(fmpz_poly_t poly_flint, const mpz_t *poly_gmp,
                            const unsigned long int deg,
                            const unsigned int nthreads){


  poly_flint->length = deg + 1;
  poly_flint->alloc = deg + 1;
  unsigned long int i;
#pragma omp parallel for private(i) num_threads(nthreads)
  for(i = 0; i <= deg; i++){
    fmpz_poly_set_coeff_mpz2(poly_flint, ((slong)i),poly_gmp[i]);
  }
}


/* from fmpz_poly_t to mpz_t *poly */
/* poly must have been allocated */

static void fmpz_poly_2_mpz(mpz_t *poly_gmp, const fmpz_poly_t poly_flint,
                            const unsigned long int deg,
                            const unsigned int nthreads){
  unsigned long int i;
#pragma omp parallel for private(i) num_threads(nthreads)
  for(i = 0; i <= deg; i++){
    fmpz_get_mpz(poly_gmp[i],(poly_flint->coeffs + i));
  }

}

#endif


/* res will contain pol1 * pol2 */
/* pol1 has degree deg1 and pol2 has degree deg2 */
/* res must be already allocated */

static void mpz_poly_mul(mpz_t *res,
                  const mpz_t *pol1,
                  const unsigned long int deg1,
                  const mpz_t *pol2,
                  const unsigned long int deg2,
                  const unsigned int nthreads){

#ifdef USEFLINT

  fmpz_poly_t res_fmpz_poly;
  fmpz_poly_t pol1_fmpz_poly;
  fmpz_poly_t pol2_fmpz_poly;

  fmpz_poly_init2(res_fmpz_poly, deg1+deg2+1);
  fmpz_poly_init2(pol1_fmpz_poly, deg1+1);
  fmpz_poly_init2(pol2_fmpz_poly, deg2+1);

  mpz_2_fmpz_poly(pol1_fmpz_poly, pol1, deg1, nthreads);
  mpz_2_fmpz_poly(pol2_fmpz_poly, pol2, deg2, nthreads);
  flint_set_num_threads(nthreads);

  fmpz_poly_mul(res_fmpz_poly, pol2_fmpz_poly, pol1_fmpz_poly);

  fmpz_poly_2_mpz(res, res_fmpz_poly, deg1+deg2, nthreads);

  fmpz_poly_clear(res_fmpz_poly);
  fmpz_poly_clear(pol1_fmpz_poly);
  fmpz_poly_clear(pol2_fmpz_poly);

#else
  fprintf(stderr, "FLINT is missing for univariate polynomial multiplication\n");
  exit(1);
#endif
}


