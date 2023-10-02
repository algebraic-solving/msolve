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

#include <flint/nmod_mat.h>
#include "nmod_mat_poly.h"

void nmod_mat_poly_mul_coeff(nmod_mat_t coeff,
                             const nmod_mat_poly_t mat1,
                             const nmod_mat_poly_t mat2,
                             slong k)
{
    // only consider indices i such that:
    //     0 <= i <= k
    //     i < mat1->length
    //     k-i < mat2->length
    // so  i < min(k+1,mat1->length)
    // and i >= max(0, k + 1 - mat2->length)
    const slong ubound = FLINT_MIN(k+1, mat1->length);
    const slong lbound = FLINT_MAX(0, k+1 - mat2->length);

    // lbound >= ubound ==> coeff is zero, no term to consider
    if (lbound >= ubound)
    {
        nmod_mat_zero(coeff);
        return;
    }

    // now lbound < ubound
    // first handle i == lbound separately, to avoid wasting time zero-ing `coeff`
    nmod_mat_mul(coeff, mat1->coeffs + lbound, mat2->coeffs + (k - lbound));

    // `if` just here to avoid initializing temp for nothing
    if (lbound + 1 < ubound)
    {
        nmod_mat_t temp;
        nmod_mat_init(temp, mat1->r, mat2->c, mat1->mod.n);
        for (slong i = lbound+1; i < ubound; i++)
        {
            nmod_mat_mul(temp, mat1->coeffs + i, mat2->coeffs + (k - i));
            nmod_mat_add(coeff, coeff, temp);
        }
        nmod_mat_clear(temp);
    }
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
