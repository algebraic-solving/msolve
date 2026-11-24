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

#include "../msolve/streams.h"

/* val / 2^(deg) is up(1/2) */
static inline void mpz_upoly_eval_onehalf(mpz_t *up, unsigned long int deg,
                                             mpz_t *val){
  long j;
  mpz_t y;

  mpz_set(*val, up[deg]);
  mpz_init(y);

  for (j = deg - 1; j >= 0; j--){
      mpz_mul_2exp(y, up[j], deg - j);
      mpz_add(*val, *val, y);
  }
  mpz_clear(y);
}

/* [val_do , val_up ] will be an interval containing
   all the values of the derivative of upoly at [0,1] */
/*
[ f_0 + negative part, f_0 + positive_part]
 */
static inline void mpz_interveval_deriv_01(mpz_t *upol, unsigned long int deg,
                                           mpz_t *vup, mpz_t *vdo){
  mpz_t x;

  mpz_set_ui(*vup, 0);
  mpz_set_ui(*vdo, 0);

  mpz_init_set(x, upol[deg]);
  mpz_mul_ui(x, x, deg);

  if(mpz_sgn(x)==1){
    mpz_add(*vup, *vup, x);
  }
  if(mpz_sgn(x)==-1){
    mpz_add(*vdo, *vdo, x);
  }

  for (unsigned long int j = deg - 1; j >= 2; j--){
      mpz_mul_ui(x, upol[j], j);
      if(mpz_sgn(x)==1){
        mpz_add(*vup, *vup, x);
      }
      if(mpz_sgn(x)==-1){
        mpz_add(*vdo, *vdo, x);
      }
  }

  mpz_add(*vdo, *vdo, upol[1]);
  mpz_add(*vup, *vup, upol[1]);
  mpz_clear(x);
}


static inline void basic_mpz_poly_eval_at_point(mpz_t *up, unsigned long int deg,
                                                mpz_t *c, mpz_t *val){
  mpz_set_ui(*val,0);
  for(unsigned int i = deg; i > 0; i--){
    mpz_add(*val, *val, up[i]);
    mpz_mul(*val, *val, *c);
  }
  mpz_add(*val, *val, up[0]);
}


static inline int sgn_mpz_poly_eval_at_point_naive(mpz_t *upoly, unsigned long int deg, mpz_t *c){
  mpz_t val;
  mpz_init(val);
  basic_mpz_poly_eval_at_point(upoly, deg, c, &val);

  int s = mpz_sgn(val);
  mpz_clear(val);
  return s;
}

/*
  Evaluates upol at c/2^k using a Horner scheme
  In the end, one has
  val/2(k*deg)=up(c / 2^k) */
void mpz_poly_eval_2exp_naive(mpz_t *upol,
                              long int deg,
                              mpz_t *c, const long k,
                              mpz_t *val, mpz_t *tmp){

  if(deg == -1){
    mpz_set_ui(*val, 0);
    return;
  }
  if(deg == 0){
    mpz_set(*val, upol[0]);
    return;
  }
  mpz_set(*val, upol[deg]);
  mpz_mul(*val, *val, *c);
  for(unsigned long int i = deg - 1; i > 0; i--){
    mpz_mul_2exp(*tmp, upol[i], (deg-i)*k);
    mpz_add(*val, *val, *tmp);
    mpz_mul(*val, *val, *c);
  }
  mpz_mul_2exp(*tmp, upol[0], deg*k);
  mpz_add(*val, *val, *tmp);

}

/* Evaluation over an interval using a Horner scheme */
/* assumes a and b have the same sign */


/* [ val_do / 2^prec, val_up / 2^prec ] will contain */
/* up( [ xdo[1], xup[1] ] ) */


/*toto*/

/* [ val_do / 2^prec, val_up / 2^prec ] will contain */
/* up( [ xdo[1], xup[1] ] ) */
/* xup[i] and xdo[i] are given mod 2^corr when i%b == 0 */



static inline int sgn_mpz_poly_eval_at_point_2exp_naive(mpz_t *upol,
                                                        unsigned long int deg,
                                                        mpz_t *c, int k){

  mpz_t val, coeff;
  mpz_init(coeff);
  mpz_init_set(val, upol[deg]);
  mpz_mul(val, val, *c);
  for(unsigned long int i=deg-1; i>0; i--){
    mpz_mul_2exp(coeff, upol[i], (deg-i)*k);
    mpz_add(val, val, coeff);
    mpz_mul(val, val, *c);
  }
  mpz_mul_2exp(coeff, upol[0], deg*k);
  mpz_add(val, val, coeff);

  int s = mpz_sgn(val);
  mpz_clear(val);
  mpz_clear(coeff);
  return s;
}



