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
#include <flint/nmod_poly_mat.h>

void nmod_mat_poly_init_set_from_nmod_mat(nmod_mat_poly_t matp,
                                          const nmod_mat_t cmat)
{
    nmod_mat_poly_init2_preinv(matp, cmat->r, cmat->c, cmat->mod.n, cmat->mod.ninv, 1);
    // if cmat is zero: do nothing more
    // otherwise set length to 1 and copy data
    if (! nmod_mat_is_zero(cmat))
    {
        nmod_mat_init_set(matp->coeffs + 0, cmat);
        matp->length = 1;
    }
}

void nmod_mat_poly_set_from_nmod_mat(nmod_mat_poly_t matp,
                                     const nmod_mat_t cmat)
{
    if (nmod_mat_is_zero(cmat))
        nmod_mat_poly_zero(matp);
    else
    {
        nmod_mat_poly_fit_length(matp, 1);
        _nmod_mat_poly_set_length(matp, 1);
        nmod_mat_set(matp->coeffs + 0, cmat);
    }
}

void nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                           const nmod_poly_mat_t pmat,
                                           slong order)
{
    const slong len = nmod_poly_mat_max_length(pmat);
    if (order > len)
        order = len;

    // allocate and init coefficients
    nmod_mat_poly_fit_length(matp, order);
    _nmod_mat_poly_set_length(matp, order);

    // fill data
    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < matp->r; i++)
            for (slong j = 0; j < matp->c; j++)
                nmod_mat_poly_entry(matp, k, i, j) = nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat, i, j), k);

    // normalize (useless if order==len)
    if (order < len)
        _nmod_mat_poly_normalise(matp);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
