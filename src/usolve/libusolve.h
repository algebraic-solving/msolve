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


interval *real_roots(mpz_t *, unsigned long,
                     unsigned long int *,
                     unsigned long int *,
                     const int32_t,
                     int,
                     int);

void display_roots_system(FILE *, interval *, unsigned long int);

void display_root(FILE *, interval *);

unsigned long int mpz_poly_max_bsize_coeffs(mpz_t *, unsigned long int);

unsigned long int mpz_poly_min_bsize_coeffs(mpz_t *, unsigned long int);

void mpz_poly_eval_2exp_naive(mpz_t *,
                              long int,
                              mpz_t *, const long int,
                              mpz_t *, mpz_t *);

void mpz_poly_eval_2exp_naive2(mpz_t *,
                              unsigned long int,
                              mpz_t, const int,
                              mpz_t, mpz_t);

void get_values_at_bounds(mpz_t *, unsigned long int, interval *, mpz_t *);

void refine_QIR_positive_root(mpz_t *, long int *, interval *,
                                     mpz_t *, int, int);

int mpz_scalar_product_interval(mpz_t *, unsigned long int, long,
                                mpz_t *, mpz_t *, mpz_t, mpz_t, mpz_t, long);

int mpz_poly_eval_interval(mpz_t *, unsigned long int, long,
                           mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);

int lazy_mpz_poly_eval_interval(mpz_t *, unsigned long int, long,
                                mpz_t *, mpz_t *, long, long, long,
                                mpz_t, mpz_t, mpz_t);

int lazy_mpz_poly_eval_interval_old(mpz_t *, unsigned long int, long,
                                mpz_t *, mpz_t *, long, long, long,
                                mpz_t, mpz_t, mpz_t);
