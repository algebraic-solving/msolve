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

void basic_mpz_poly_eval_at_point(mpz_t *, unsigned long int, mpz_t*, mpz_t *);

void mpz_poly_eval_at_point_2exp(mpz_t *, unsigned long int,
                                 unsigned long int, unsigned long int,
                                 mpz_t, int,
                                 int, mpz_t *, mpz_t *,
                                 double *);

void mpz_poly_eval_at_point_2exp_ui(mpz_t *, unsigned long int,
                                 unsigned long int, unsigned long int,
                                 unsigned long int, int,
                                 int, mpz_t *, mpz_t *,
                                 double *);

int sgn_mpz_poly_eval_at_point_2exp(mpz_t *, unsigned long int,
                                    unsigned long int, unsigned long int,
                                    mpz_t *, int,
                                    mpz_t *, mpz_t *,
                                    int *, double *);

int sgn_mpz_poly_eval_at_point_2exp_ui(mpz_t *, unsigned long int,
                                    unsigned long int, unsigned long int,
                                    unsigned long int, int,
                                    mpz_t *, mpz_t *,
                                    int *, double *);

int sgn_mpz_poly_eval_at_point_naive(mpz_t *, unsigned long int, mpz_t *, int);

void mpz_poly_eval_2exp_naive(mpz_t *, long int,
                                       mpz_t *, const long int, mpz_t *, mpz_t *);

int sgn_mpz_poly_eval_at_point_2exp_naive(mpz_t *, unsigned long int, mpz_t*, int);
int sgn_mpz_poly_eval_at_point_2exp_naive2(mpz_t *, unsigned long int, mpz_t, int);
