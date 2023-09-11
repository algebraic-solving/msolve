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
 * Vincent Neiger
 * Mohab Safey El Din */

#include "nmod_mat_poly.h"

void nmod_mat_poly_print(const nmod_mat_poly_t matp)
{
	for (slong i = 0; i < matp->length; i++)
		nmod_mat_print(matp->coeffs + i);
}

void nmod_mat_poly_print_pretty(const nmod_mat_poly_t matp)
{
    printf("r: %ld, c: %ld, len: %ld, mod: %ld\n",
           matp->r, matp->c, matp->length, matp->mod.n);
	for (slong i = 0; i < matp->length; i++)
    {
        printf("coeff %ld:\n", i);
        nmod_mat_print_pretty(matp->coeffs + i);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
