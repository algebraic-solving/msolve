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
#include <gmp.h>
#include <math.h>
#include <time.h>



void mpz_poly_mul(mpz_t *, mpz_t *, unsigned long int, mpz_t *, unsigned long int, unsigned int);
#ifdef USEFLINT
void mpz_poly_mul_allocated_fmpz_poly(mpz_t *, mpz_t *, unsigned long int, mpz_t *, unsigned long int,
                                      fmpz_poly_t, fmpz_poly_t, fmpz_poly_t);

void mpz_upoly_mul(mpz_t *, mpz_t *, mpz_t *, unsigned long int, unsigned long int,
                   fmpz_poly_t, fmpz_poly_t, fmpz_poly_t,
                   int);
#endif
