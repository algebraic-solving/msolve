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

#include"descartes.h"

interval *Bisection_Uspensky(mpz_t *, unsigned long,
                             unsigned long int *, unsigned long int *,
                             usolve_flags *);

interval *Bisection_Uspensky_from_root(mpz_t *, unsigned long,
                                       unsigned long int *, unsigned long int *,
                                       usolve_flags *);
