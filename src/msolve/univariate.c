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

/* Univariate input over the rationals.

   When the input polynomials f_1, ..., f_m have one variable x, their
   solutions are the roots of g = gcd(f_1, ..., f_m), and the
   parametrization of the radical of the ideal they generate is given by
   the squarefree part w = g / gcd(g, g') of g: x = t with w(t) = 0, the
   denominator being w'. These are computed with FLINT, without Groebner
   bases, and the real roots are then isolated from this parametrization as
   for any other one. The degree of the quotient is the degree of g. When
   all the input polynomials are 0 (their terms cancel out), the solution
   set is the whole line. */

#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>

static void mpz_upoly_init(mpz_upoly_t poly, deg_t alloc);
void mpz_upoly_clear(mpz_upoly_t pol);

/* f: the i-th input polynomial, whose terms start at pos, times the least
 * common multiple of the denominators of its coefficients, divided by its
 * content (several terms may have the same exponent) */
static void univariate_input_poly(fmpz_poly_t f, const data_gens_ff_t *gens,
                                  const long i, const long pos){
  fmpq_poly_t q;
  fmpq_t c, t;
  fmpq_poly_init(q);
  fmpq_init(c);
  fmpq_init(t);
  for(long j = 0; j < gens->lens[i]; j++){
    const slong e = gens->exps[pos + j];
    fmpz_set_mpz(fmpq_numref(c), *(gens->mpz_cfs[2 * (pos + j)]));
    fmpz_set_mpz(fmpq_denref(c), *(gens->mpz_cfs[2 * (pos + j) + 1]));
    fmpq_canonicalise(c);
    fmpq_poly_get_coeff_fmpq(t, q, e);
    fmpq_add(t, t, c);
    fmpq_poly_set_coeff_fmpq(q, e, t);
  }
  fmpq_poly_get_numerator(f, q);
  fmpz_poly_primitive_part(f, f);
  fmpq_poly_clear(q);
  fmpq_clear(c);
  fmpq_clear(t);
}

/* p: the polynomial f, of length len >= 1 */
static void univariate_set_upoly(mpz_upoly_t p, const fmpz_poly_t f,
                                 const slong len){
  mpz_upoly_clear(p);
  mpz_upoly_init(p, len);
  for(slong i = 0; i < len; i++){
    fmpz_get_mpz(p->coeffs[i], f->coeffs + i);
  }
  p->length = len;
}

/* for the input polynomials given by gens (one variable, rational
 * coefficients): sets the parametrization of their solutions in mpz_param
 * (when there are solutions), *dim_ptr to 0 and *dquot_ptr to the degree of
 * their gcd, or, when they are all 0, *dim_ptr to 1 and *dquot_ptr to -1 */
static void univariate_param_qq(mpz_param_t mpz_param, int *dim_ptr,
                               long *dquot_ptr, const data_gens_ff_t *gens,
                               const int32_t info_level){
  fmpz_poly_t f, g, d;
  fmpz_poly_init(f);
  fmpz_poly_init(g);
  fmpz_poly_init(d);
  long pos = 0;
  for(long i = 0; i < gens->ngens; i++){
    univariate_input_poly(f, gens, i, pos);
    pos += gens->lens[i];
    fmpz_poly_gcd(g, g, f);
  }
  if(fmpz_poly_is_zero(g)){
    mpz_param->dim = 1;
    *dim_ptr = 1;
    *dquot_ptr = -1;
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_poly_clear(d);
    return;
  }
  const slong dq = fmpz_poly_degree(g);
  slong ds = 0;
  if(dq > 0){
    /* w = g / gcd(g, g'), primitive with a positive leading coefficient
     * since g is (fmpz_poly_squarefree_part is public since FLINT 3.6) */
#if defined(__FLINT_RELEASE) && __FLINT_RELEASE >= 30600
    fmpz_poly_squarefree_part(f, g);
#else
    fmpz_poly_derivative(d, g);
    fmpz_poly_gcd(d, g, d);
    fmpz_poly_div(f, g, d);
#endif
    ds = fmpz_poly_degree(f);
    univariate_set_upoly(mpz_param->elim, f, ds + 1);
    fmpz_poly_derivative(d, f);
    univariate_set_upoly(mpz_param->denom, d, ds);
  }
  mpz_param->nvars = 1;
  mpz_param->nsols = ds;
  mpz_param->dquot = dq;
  mpz_param->dim = 0;
  *dim_ptr = 0;
  *dquot_ptr = dq;
  if(info_level){
    fprintf(VERBSTREAM, "\nUnivariate input: gcd of degree %ld, squarefree "
            "part of degree %ld\n", (long)dq, (long)ds);
  }
  fmpz_poly_clear(f);
  fmpz_poly_clear(g);
  fmpz_poly_clear(d);
}
