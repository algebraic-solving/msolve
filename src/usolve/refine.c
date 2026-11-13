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

#ifndef USOLVE
#define USOLVE
#endif

/* halves the isolating interval rt of a root of upol */
static void refine_root_naive(mpz_t *upol, unsigned long int deg,
                              interval *rt, mpz_t *middle){
  if(rt->isexact == 1){
    return;
  }

  long newk;
  int sgn_middle;
  if((rt->k) >= 0){
    /* interval is (numer, numer + 2^(-k)) */
    mpz_mul_ui(*middle, rt->numer, 2);
    mpz_add_ui(*middle, *middle, 1);
    newk = rt->k + 1;

    sgn_middle = sgn_mpz_poly_eval_at_point_2exp_naive(upol, deg, middle, newk);
  }
  else{
    mpz_set_ui(*middle, 1);
    mpz_mul_2exp(*middle, *middle, - ( (rt->k) + 1));
    mpz_add(*middle, *middle, rt->numer);
    newk = ( (rt->k) + 1 );

    sgn_middle = sgn_mpz_poly_eval_at_point_naive(upol, deg, middle);
  }
  int sign_left = rt->sign_left;

  if(sgn_middle * sign_left<0){
    if(newk>0){
      mpz_mul_ui(rt->numer, rt->numer, 2);
    }
    rt->k = newk;
  }
  else{
    mpz_set(rt->numer, *middle);
    rt->k = newk;
  }
}

/* refines the isolating intervals of the nb roots of upol by bisection, to
 * k >= prec */
void refine_all_roots_naive(mpz_t *upol, unsigned long int deg,
                            interval *roots, unsigned long int nb,
                            unsigned int prec){
  mpz_t *middle=malloc(sizeof(mpz_t));
  mpz_init(middle[0]);

  for(unsigned long int j = 0; j < nb; j++){
    while((roots+j)->k < prec && (roots+j)->isexact == 0){
      refine_root_naive(upol, deg, roots+j, middle);
    }
  }
  mpz_clear(middle[0]);
  free(middle);
}
