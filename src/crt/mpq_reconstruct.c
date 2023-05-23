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

/**

   This implementation is a (very) slight modification of the functions in FLINT.

**/

#include <gmp.h>

/* #define ROT(u,v,t)                                            \ */
/*   do { mpz _t = *u; *u = *v; *v = *t; *t = _t; } while (0); */

static inline void ROT(mpz_t u, mpz_t v, mpz_t t){
    mpz_swap(u, v);
    mpz_swap(v, t);
}

int
_mpq_reconstruct_mpz_2(mpz_t n, mpz_t d,
    const mpz_t a, const mpz_t m, const mpz_t N, const mpz_t D)
{

    mpz_t q, r, s, t;
    int success = 0;

    /* Quickly identify small integers */
    if (mpz_cmp(a, N) <= 0)
    {
        mpz_set(n, a);
        mpz_set_ui(d, 1);
        return 1;
    }
    mpz_sub(n, a, m);
    if(mpz_cmp_ui(n, 0)>=0){
      if(mpz_cmp(n, N) <= 0){
        mpz_set_ui(d, 1);
        return 1;
      }
    }
    else{
      mpz_neg(n, n);
      if(mpz_cmp(n, N) <= 0){
        mpz_set_ui(d, 1);
        mpz_neg(n, n);
        return 1;
      }
    }
    /* if (fmpz_cmpabs(n, N) <= 0) */
    /* { */
    /*     fmpz_one(d); */
    /*     fprintf(stderr, "ici?\n"); */
    /*     return 1; */
    /* } */

    mpz_init(q);
    mpz_init(r);
    mpz_init(s);
    mpz_init(t);

    mpz_set(r, m); mpz_set_ui(s, 0);
    mpz_set(n, a); mpz_set_ui(d, 1);

    while (mpz_cmpabs(n, N) > 0)
    {
        mpz_fdiv_q(q, r, n);
        mpz_mul(t, q, n); mpz_sub(t, r, t); ROT(r, n, t);
        mpz_mul(t, q, d); mpz_sub(t, s, t); ROT(s, d, t);
    }

    if (mpz_sgn(d) < 0)
    {
        mpz_neg(n, n);
        mpz_neg(d, d);
    }

    if (mpz_cmp(d, D) <= 0)
    {
        mpz_gcd(t, n, d);
        if(mpz_cmp_ui(t, 1)==0){
          success = 1;
        }
        else{
          success = 0;
        }
    }

    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(s);
    mpz_clear(t);
    return success;
}

int
mpq_reconstruct_mpz_2(mpq_t res, const mpz_t a, const mpz_t m,
                          const mpz_t N, const mpz_t D)
{
  return _mpq_reconstruct_mpz_2(mpq_numref(res),
                                mpq_denref(res), a, m, N, D);
}


int _mpq_reconstruct_mpz(mpz_t n, mpz_t d,
                         const mpz_t a, const mpz_t m)
{
  mpz_t N;
  int result;

  mpz_init(N);
  mpz_fdiv_q_2exp(N, m, 1);
  mpz_sqrt(N, N);
  result = _mpq_reconstruct_mpz_2(n, d, a, m, N, N);
  mpz_clear(N);

  return result;
}

int _mpq_reconstruct_mpz_with_denom(mpz_t n, mpz_t d,
                                    const mpz_t a, const mpz_t m,
                                    mpz_t num, mpz_t den)
{
  //  mpz_t N;
  int result;

  //  mpz_init(N);
  //  mpz_fdiv_q_2exp(N, m, 1);
  //  mpz_sqrt(N, N);

  result = _mpq_reconstruct_mpz_2(n, d, a, m, num, den);
  //  mpz_clear(N);

  return result;
}

int mpq_reconstruct_mpz(mpq_t *res, mpz_t a, const mpz_t m)
{
  if(mpz_cmp_ui(a, 0)>=0){
    return _mpq_reconstruct_mpz(mpq_numref(*res),
                                mpq_denref(*res), a, m);
  }
  else{
    while(mpz_cmp_ui(a, 0) < 0){
      //      mpz_fprint(stderr, a); fprintf(stderr, "\n");
      mpz_add(a, a, m);
    }
    int b = _mpq_reconstruct_mpz(mpq_numref(*res),
                                  mpq_denref(*res), a, m);
    return b;
  }
}

int mpq_reconstruct_mpz_with_denom(mpq_t *res, mpz_t a, const mpz_t m,
                                   mpz_t num, mpz_t den)
{
  if(mpz_cmp_ui(a, 0)>=0){
    return _mpq_reconstruct_mpz_with_denom(mpq_numref(*res),
                                           mpq_denref(*res), a, m,
                                           num, den);
  }
  else{
    while(mpz_cmp_ui(a, 0) < 0){
      mpz_add(a, a, m);
    }
    int b = _mpq_reconstruct_mpz_with_denom(mpq_numref(*res),
                                            mpq_denref(*res), a, m,
                                            num, den);
    return b;
  }
}


/* Rational reconstruction -- classical algorithm
   returns 1 in case of success else returns 0

   assumes that 2ND < mod

   In case of success, n and d are such that
   n/d = u modulo mod
   gcd(n,d)=1
   and |n| < recadata->N, |d| < recdata->D


 */
int ratreconwden(mpz_t n, mpz_t d, /* output numerator and denominator */
                 mpz_t u, const mpz_t mod, const mpz_t gden,
                 rrec_data_t recdata){

  /* while(mpz_cmp_ui(u, 0) < 0){ */
  /*   mpz_add(u, u, mod); */
  /* } */

  mpz_mod(u, u, mod);
  mpz_set(recdata->r0, mod);
  mpz_set_ui(recdata->t0, 0);

  mpz_set(recdata->r1, u);
  mpz_mul(recdata->r1, recdata->r1, gden);
  mpz_mod(recdata->r1, recdata->r1, mod);
  mpz_set_ui(recdata->t1, 1);

  while(mpz_cmp(recdata->r1, recdata->N)>0){
    mpz_fdiv_q(recdata->q, recdata->r0, recdata->r1);

    mpz_mul(recdata->tmp, recdata->q, recdata->r1);
    mpz_sub(recdata->tmp, recdata->r0, recdata->tmp);
    mpz_swap(recdata->r0, recdata->r1);
    mpz_swap(recdata->r1, recdata->tmp);

    mpz_mul(recdata->tmp, recdata->q, recdata->t1);
    mpz_sub(recdata->tmp, recdata->t0, recdata->tmp);
    mpz_swap(recdata->t0, recdata->t1);
    mpz_swap(recdata->t1, recdata->tmp);
  }
  mpz_set(n, recdata->r1);
  mpz_set(d, recdata->t1);

  if(mpz_sgn(d) < 0){
    mpz_neg(n, n);
    mpz_neg(d, d);
  }
  mpz_gcd(recdata->q, n, d);
  if(mpz_cmp(d, recdata->D) <= 0 && mpz_cmp_ui(recdata->q, 1)==0){
    return 1;
  }

  return 0;
}

/* Rational reconstruction -- classical algorithm
   returns 1 in case of success else returns 0

   assumes that 2ND < mod

   In case of success, n and d are such that
   n/d = u modulo mod
   gcd(n,d)=1
   and |n| < recadata->N, |d| < recdata->D


 */
int ratrecon(mpz_t n, mpz_t d, /* output numerator and denominator */
             mpz_t u, const mpz_t mod,
             rrec_data_t recdata){

  while(mpz_cmp_ui(u, 0) < 0){
    mpz_add(u, u, mod);
  }

  mpz_set(recdata->r0, mod);
  mpz_set_ui(recdata->t0, 0);

  mpz_set(recdata->r1, u);
  mpz_set_ui(recdata->t1, 1);

  while(mpz_cmp(recdata->r1, recdata->N)>0){
    mpz_fdiv_q(recdata->q, recdata->r0, recdata->r1);

    mpz_mul(recdata->tmp, recdata->q, recdata->r1);
    mpz_sub(recdata->tmp, recdata->r0, recdata->tmp);
    mpz_swap(recdata->r0, recdata->r1);
    mpz_swap(recdata->r1, recdata->tmp);

    mpz_mul(recdata->tmp, recdata->q, recdata->t1);
    mpz_sub(recdata->tmp, recdata->t0, recdata->tmp);
    mpz_swap(recdata->t0, recdata->t1);
    mpz_swap(recdata->t1, recdata->tmp);
  }
  mpz_set(n, recdata->r1);
  mpz_set(d, recdata->t1);

  if(mpz_sgn(d) < 0){
    mpz_neg(n, n);
    mpz_neg(d, d);
  }
  mpz_gcd(recdata->q, n, d);
  if(mpz_cmp(d, recdata->D) <= 0 && mpz_cmp_ui(recdata->q, 1)==0){
    return 1;
  }

  return 0;
}
