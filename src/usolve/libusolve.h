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

#include <flint/fmpz_poly.h>

void usolve_refine_root_arb(interval *, const fmpz_poly_t, const fmpz_poly_t,
                            const slong, slong *);

int usolve_sign_2exp_arb(const fmpz_poly_t, const fmpz_t, const slong,
                         slong *);

int usolve_sign_right_2exp_arb(const fmpz_poly_t, const fmpz_poly_t,
                               const fmpz_t, const slong, slong *);

interval *real_roots(mpz_t *, unsigned long,
                     unsigned long int *,
                     unsigned long int *,
                     const int32_t,
                     int,
                     int,
                     usolve_stats_t *);

unsigned long int mpz_poly_max_bsize_coeffs(mpz_t *, unsigned long int);

void mpz_poly_eval_2exp_naive(mpz_t *,
                              long int,
                              mpz_t *, const long int,
                              mpz_t *, mpz_t *);

