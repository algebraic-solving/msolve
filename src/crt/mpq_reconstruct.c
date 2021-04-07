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
  mpz_t _t;
  mpz_init(_t);
  mpz_swap(_t, u);
  mpz_swap(u,v);
  mpz_swap(v,t);
  mpz_swap(t, _t);
  mpz_clear(_t);
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
      //      mpz_fprint(stderr, a); fprintf(stderr, "\n");
      mpz_add(a, a, m);
    }
    int b = _mpq_reconstruct_mpz_with_denom(mpq_numref(*res),
                                            mpq_denref(*res), a, m,
                                            num, den);
    return b;
  }
}
