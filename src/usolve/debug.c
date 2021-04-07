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

#include <gmp.h>
//#include "utils.h"
//#include "evaluate.h"
//#include "print_usolve.h"

/* returns 1 if up has sign change over the interval encoded by root */
/* return 0 otherwise */
static inline int check_sign_change(mpz_t *upoly, unsigned long int deg, interval *root){
  if(root->isexact==1){
    mpz_t c;
    mpz_init_set(c, root->numer);
    //  mpz_out_str(stderr, 10, c); puts("");
    int s_down = sgn_mpz_poly_eval_at_point_2exp_naive(upoly, deg, &c, root->k);
    if(s_down!=0){
      fprintf(stderr, "Negative (k): Not implemented yet ; skipped but useless with larger precision\n");
      return 0; //    exit(1);
    }
    return 1;
  }
  if(root->k<=0){
    fprintf(stderr, "Negative (k): Not implemented yet ; skipped but useless with larger precision\n");
    return 1; //    exit(1);
  }
  mpz_t c;
  mpz_init_set(c, root->numer);
  //  mpz_out_str(stderr, 10, c); puts("");
  int s_down = sgn_mpz_poly_eval_at_point_2exp_naive(upoly, deg, &c, root->k);
  mpz_add_ui(c, c, 1);
  //  mpz_out_str(stderr, 10, c); puts("");
  int s_up = sgn_mpz_poly_eval_at_point_2exp_naive(upoly, deg, &c, root->k);
  mpz_clear(c);
  if(s_up*s_down<=0){
    return 1;
  }
  fprintf(stderr, "\nBUG: no sign change\n");
  return 0;
}

static int check_all_sign_changes(mpz_t *upoly, unsigned long int deg, interval *root, unsigned long int nbroots){
  int s = 1;
  for(unsigned long int i =0; i < nbroots; i++){
    if((root+i)->isexact==0){
      s = check_sign_change(upoly, deg, root+i);
      if(s==0){
        fprintf(stderr, "Error occurred at root %lu\n", i);
        USOLVEdisplay_roots(stderr, (root+i), 1);
        fprintf(stderr, "\n");
        exit(1);
      }
    }
  }
  return 1;
}

