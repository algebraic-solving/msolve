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


static inline int sgn_mpz_poly_eval_at_point_naive(mpz_t *upoly, unsigned long int deg, mpz_t *c, int classical_algo){
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

/*
  Evaluates upol at c/2^k using a Horner scheme
  In the end, one has 
  val/2(k*deg)=up(c / 2^k) */
void mpz_poly_eval_2exp_naive2(mpz_t *upol,
                              long int deg,
                              mpz_t c, const long k,
                              mpz_t val, mpz_t tmp){

  if(deg == -1){
    mpz_set_ui(val, 0);
    return;
  }
  if(deg == 0){
    mpz_set(val, upol[0]);
    return;
  }
  mpz_set(val, upol[deg]);
  mpz_mul(val, val, c);
  for(unsigned long int i = deg - 1; i > 0; i--){
    mpz_mul_2exp(tmp, upol[i], (deg-i)*k);
    mpz_add(val, val, tmp);
    mpz_mul(val, val, c);
  }
  mpz_mul_2exp(tmp, upol[0], deg*k);
  mpz_add(val, val, tmp);

}

/* Evaluation over an interval using a Horner scheme */
/* assumes a and b have the same sign */
/* returns (mpz_sgn(val_do) != mpz_sgn(val_up)); */
int mpz_poly_eval_interval(mpz_t *up, const long int deg, const long k,
                           mpz_t a, mpz_t b,
                           mpz_t tmp,
                           mpz_t val_do, mpz_t val_up){
  if(deg == -1){
    mpz_set_ui(val_up, 0);
    mpz_set_ui(val_do, 0);
    return 0;
  }
  if(deg == 0){
    mpz_set(val_up, up[0]);
    mpz_set(val_do, up[0]);
    return 0;
  }
  mpz_set(val_up, up[deg]);
  mpz_set(val_do, up[deg]);
  mpz_t s;
  mpz_init(s);

  if(mpz_sgn(a) != mpz_sgn(b) && mpz_sgn(a) != 0 && mpz_sgn(b) != 0){
    fprintf(stderr, "Entries of mpz_poly_eval_interval are incorrect\n");
    fprintf(stderr, "a and b should have same sign");
    mpz_out_str(stderr, 10, a);
    mpz_out_str(stderr, 10, b);
    mpz_clear(s);
    exit(1);
  }
  /* a and b are non-negative */
  if(mpz_sgn(a) >= 0){
    for(long i = deg - 1; i >= 0; i--){

      if(mpz_sgn(val_do) >= 0){
        mpz_mul_2exp(tmp, up[i], k*(deg - i));

        mpz_mul(val_do, val_do, a);
        mpz_add(val_do, val_do, tmp);

        mpz_mul(val_up, val_up, b);
        mpz_add(val_up, val_up, tmp);
      }
      else{
        mpz_mul_2exp(tmp, up[i], k*(deg - i));

        mpz_mul(val_do, val_do, b);
        mpz_add(val_do, val_do, tmp);

        mpz_mul(val_up, val_up, a);
        mpz_add(val_up, val_up, tmp);
      }

      if(mpz_cmp(val_do, val_up) > 0){
        fprintf(stderr, "BUG ici (den_do > den_up)\n");
        exit(1);
      }

    }
  }
  else{
    /* a and b are negative */
    for(long i = deg - 1; i >= 0; i--){
      if(mpz_sgn(val_do) >= 0){
        mpz_set(s, val_do);
        mpz_mul_2exp(tmp, up[i], k*(deg - i));

        mpz_mul(val_do, val_up, a);
        mpz_add(val_do, val_do, tmp);

        mpz_mul(val_up, s, b);
        mpz_add(val_up, val_up, tmp);

      }
      else{
        mpz_set(s, val_do);
        mpz_mul_2exp(tmp, up[i], k*(deg - i));

        mpz_mul(val_do, val_up, b);
        mpz_add(val_do, val_do, tmp);

        mpz_mul(val_up, s, a);
        mpz_add(val_up, val_up, tmp);
      }

      if(mpz_cmp(val_do, val_up) > 0){
        fprintf(stderr, "BUG ici2 (val_do > val_up)\n");
        fprintf(stderr, "=> sign of val_do = %d\n", mpz_sgn(val_do));
        fprintf(stderr, "=> sign of val_up = %d\n", mpz_sgn(val_up));
        exit(1);
      }
    }
  }
  mpz_clear(s);
  return (mpz_sgn(val_do) != mpz_sgn(val_up));
}

/* evaluates up (of degree deg) over the interval
   through scalar product up*xdo and up*xup
   up to precision corr

   Values are stored in val_do, val_up
   => [val_do / 2^corr, val_up/2^corr]

   returns (mpz_sgn(val_do) != mpz_sgn(val_up))
 */
int mpz_scalar_product_interval(mpz_t *up, const long int deg, const long k,
                                mpz_t *xdo, mpz_t *xup,
                                mpz_t tmp,
                                mpz_t val_do, mpz_t val_up,
                                long corr){
  if(deg == -1){
    mpz_set_ui(val_up, 0);
    mpz_set_ui(val_do, 0);
    return 0;
  }
  if(deg == 0){
    mpz_set(val_up, up[0]);
    mpz_set(val_do, up[0]);
    mpz_mul_2exp(val_up, val_up, corr);
    mpz_mul_2exp(val_do, val_do, corr);
    return 0;
  }
  mpz_set_ui(val_up, 0);
  mpz_set_ui(val_do, 0);

  for(long i = 0; i <= deg; i++){
    if(mpz_sgn(up[i]) >= 0){
      mpz_mul(tmp, up[i], xup[i]);
      mpz_add(val_up, val_up, tmp);
      mpz_mul(tmp, up[i], xdo[i]);
      mpz_add(val_do, val_do, tmp);
    }
    else{
      mpz_mul(tmp, up[i], xdo[i]);
      mpz_add(val_up, val_up, tmp);
      mpz_mul(tmp, up[i], xup[i]);
      mpz_add(val_do, val_do, tmp);
    }
  }
  return (mpz_sgn(val_do) != mpz_sgn(val_up));
}



/* [ val_do / 2^prec, val_up / 2^prec ] will contain */
/* up( [ xdo[1], xup[1] ] ) */


/*toto*/

/* [ val_do / 2^prec, val_up / 2^prec ] will contain */
/* up( [ xdo[1], xup[1] ] ) */
/* xup[i] and xdo[i] are given mod 2^corr when i%b == 0 */

int lazy_mpz_poly_eval_interval(mpz_t *up, const unsigned long int deg,
                                const long k,
                                mpz_t *xdo, mpz_t *xup,
                                const long prec, const long corr,
                                const long b,
                                mpz_t tmp,
                                mpz_t val_do, mpz_t val_up){

  if(deg==-1){
    mpz_set_ui(val_up, 0);
    mpz_set_ui(val_do, 0);
    return 0;
  }
  if(deg==0){
    mpz_set(val_up, up[0]);
    mpz_set(val_do, up[0]);
    return 0;
  }

  mpz_set_ui(val_up, 0);
  mpz_set_ui(val_do, 0);
  const long t = (deg) / b;
  long rd = (deg) % b;
  mpz_t fdo, fup;
  mpz_init(fdo);
  mpz_init(fup);

  for(long j = 0; j < t; j++){
    mpz_set_ui(fdo, 0);
    mpz_set_ui(fup, 0);
    for(long i = 0; i < b; i++){
      if(mpz_sgn(up[i+j*b]) >= 0){
        mpz_mul(tmp, up[i+j*b], xup[i]);
        mpz_mul_2exp(tmp, tmp, k*(b - 1 - i));
        mpz_add(fup, fup, tmp);
        mpz_mul(tmp, up[i+j*b], xdo[i]);
        mpz_mul_2exp(tmp, tmp, k*(b - 1 - i));
        mpz_add(fdo, fdo, tmp);
      }
      else{
        mpz_mul(tmp, up[i+j*b], xdo[i]);
        mpz_mul_2exp(tmp, tmp, k*(b - 1 - i));
        mpz_add(fup, fup, tmp);
        mpz_mul(tmp, up[i+j*b], xup[i]);
        mpz_mul_2exp(tmp, tmp, k*(b - 1 - i));
        mpz_add(fdo, fdo, tmp);
      }
    }

    if(mpz_cmp(fdo, fup) > 0){
      fprintf(stderr, "BUG in preprocess eval (fdo > fup)\n");
      mpz_out_str(stderr, 10, fdo); fprintf(stderr, "\n");
      mpz_out_str(stderr, 10, fup); fprintf(stderr, "\n");
      exit(1);
    }

    if(mpz_sgn(fdo) < 0){
      mpz_mul(fdo, fdo, xup[j*b]);
    }
    else{
      mpz_mul(fdo, fdo, xdo[j*b]);
    }
    if(mpz_sgn(fup) < 0){
      mpz_mul(fup, fup, xdo[j*b]);
    }
    else{
      mpz_mul(fup, fup, xup[j*b]);
    }

    mpz_mul_2exp(fdo, fdo, t+prec);
    mpz_mul_2exp(fup, fup, t+prec);
    if(j){
      mpz_fdiv_q_2exp(fdo, fdo, k*(b-1) + corr);
      mpz_cdiv_q_2exp(fup, fup, k*(b-1) + corr);
    }
    else{
      mpz_fdiv_q_2exp(fdo, fdo, k*(b-1) );
      mpz_cdiv_q_2exp(fup, fup, k*(b-1) );

    }

    mpz_add(val_do, val_do, fdo);
    mpz_add(val_up, val_up, fup);


    if(mpz_cmp(fdo, fup) > 0){
      fprintf(stderr, "BUG in preprocess2 eval (fdo > fup)\n");
      mpz_out_str(stderr, 10, xdo[j*b]); fprintf(stderr, "\n");
      mpz_out_str(stderr, 10, xup[j*b]); fprintf(stderr, "\n");
      fprintf(stderr, "cmp = %d\n", mpz_cmp(xdo[j*b], xup[j*b]));
      exit(1);
    }

    if(mpz_cmp(val_do, val_up) > 0){
      fprintf(stderr, "BUG in eval (val_do > val_up)\n");
      mpz_out_str(stderr, 10, val_do); fprintf(stderr, "\n");
      mpz_out_str(stderr, 10, val_up); fprintf(stderr, "\n");
      exit(1);
    }
  }

  if(rd){
    mpz_set_ui(fdo, 0);
    mpz_set_ui(fup, 0);

    for(long i = 0; i <= rd; i++ ){
      if(mpz_sgn(up[i+t*b]) >= 0){
        mpz_mul(tmp, up[i+t*b], xup[i]);
        mpz_mul_2exp(tmp, tmp, k*(rd - i));
        mpz_add(fup, fup, tmp);
        mpz_mul(tmp, up[i+t*b], xdo[i]);
        mpz_mul_2exp(tmp, tmp, k*(rd - i));
        mpz_add(fdo, fdo, tmp);
      }
      else{
        mpz_mul(tmp, up[i+t*b], xdo[i]);
        mpz_mul_2exp(tmp, tmp, k*(rd - i));
        mpz_add(fup, fup, tmp);
        mpz_mul(tmp, up[i+t*b], xup[i]);
        mpz_mul_2exp(tmp, tmp, k*(rd - i));
        mpz_add(fdo, fdo, tmp);
      }
    }

    if(mpz_cmp(fdo, fup) > 0){
      fprintf(stderr, "BUG in preprocess3 init eval (fdo > fup)\n");
      exit(1);
    }

    if(mpz_cmp(val_do, val_up) > 0){
      fprintf(stderr, "BUG in eval (val_do > val_up)\n");
      exit(1);
    }
    if(mpz_cmp(fdo, fup) > 0){
      fprintf(stderr, "BUG in preprocess3 init0 eval (fdo > fup)\n");
      exit(1);
    }

    if(mpz_sgn(fdo) < 0){
      mpz_mul(fdo, fdo, xup[t*b]);
    }
    else{
      mpz_mul(fdo, fdo, xdo[t*b]);
    }
    if(mpz_sgn(fup) < 0){
      mpz_mul(fup, fup, xdo[t*b]);
    }
    else{
      mpz_mul(fup, fup, xup[t*b]);
    }

    mpz_mul_2exp(fdo, fdo, t+prec);
    mpz_mul_2exp(fup, fup, t+prec);

    if(t){
      mpz_cdiv_q_2exp(fup, fup, k*(rd) + corr);
      mpz_fdiv_q_2exp(fdo, fdo, k*(rd) + corr);
    }
    else{
      mpz_cdiv_q_2exp(fup, fup, k*(rd) );
      mpz_fdiv_q_2exp(fdo, fdo, k*(rd) );
    }

    mpz_add(val_do, val_do, fdo);
    mpz_add(val_up, val_up, fup);

  }

  mpz_mul_2exp(val_do, val_do, prec);
  mpz_mul_2exp(val_up, val_up, prec);
  mpz_fdiv_q_2exp(val_do, val_do, t+prec);
  mpz_cdiv_q_2exp(val_up, val_up, t+prec);

  mpz_clear(fdo);
  mpz_clear(fup);

  return (mpz_sgn(val_do) != mpz_sgn(val_up));
}



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



