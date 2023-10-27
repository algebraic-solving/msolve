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

#include <flint/nmod_poly.h>
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

void nmod_poly_mat_shift_right(nmod_poly_mat_t smat, const nmod_poly_mat_t pmat, slong k)
{
    for (slong i = 0; i < smat->r; i++)
        for (slong j = 0; j < smat->c; j++)
            nmod_poly_shift_right(smat->rows[i] + j, pmat->rows[i] + j, k);
}

// computes x**(-d) (A*B mod x**h)
void nmod_poly_mat_middle_product(nmod_poly_mat_t res,
                                  const nmod_poly_mat_t A,
                                  const nmod_poly_mat_t B,
                                  slong d,
                                  slong h)
{
    nmod_poly_mat_mul(res, A, B);
    nmod_poly_mat_truncate(res, h);
    nmod_poly_mat_shift_right(res, res, d);
}

void nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                           const nmod_mat_poly_t matp,
                                           slong order)
{
    if (order > matp->length)
        order = matp->length;

    // prepare memory
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_fit_length(nmod_poly_mat_entry(pmat, i, j), order);

    // fill data
    for (slong k = 0; k < order; k++)
        for (slong i = 0; i < pmat->r; i++)
            for (slong j = 0; j < pmat->c; j++)
                nmod_poly_mat_entry(pmat, i, j)->coeffs[k] = nmod_mat_poly_entry(matp, k, i, j);

    // normalize
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
        {
            _nmod_poly_set_length(nmod_poly_mat_entry(pmat, i, j), order);
            _nmod_poly_normalise(nmod_poly_mat_entry(pmat, i, j));
        }
}

void nmod_poly_mat_degree_matrix(fmpz_mat_t dmat,
                   const nmod_poly_mat_t mat)
{
    for(slong i = 0; i < mat->r; i++)
        for(slong j = 0; j < mat->c; j++)
            *fmpz_mat_entry(dmat, i, j) = nmod_poly_degree(nmod_poly_mat_entry(mat, i, j));
}

#if __FLINT_VERSION < 3
void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len)
{
    for (slong i = 0; i < pmat->r; i++)
        for (slong j = 0; j < pmat->c; j++)
            nmod_poly_truncate(pmat->rows[i] + j, len);
}
#endif


#if __FLINT_VERSION < 3
void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * var)
{
    slong rdim = mat->r, cdim = mat->c;

    flint_printf("<%wd x %wd matrix over Z/nZ[%s]>\n", mat->r, mat->c, var);
    flint_printf("[");
    for (slong i = 0; i < rdim; i++)
    {
        flint_printf("[");
        for (slong j = 0; j < cdim; j++)
        {
            nmod_poly_print_pretty(nmod_poly_mat_entry(mat, i, j), var);
            if (j+1 < cdim)
                flint_printf(", ");
        }
        if (i != rdim -1)
            flint_printf("],\n");
        else
            flint_printf("]");
    }
    flint_printf("]\n");
}
#endif


void nmod_poly_mat_degree_matrix_print_pretty(const nmod_poly_mat_t mat)
{
    fmpz_mat_t dmat;
    fmpz_mat_init(dmat, mat->r, mat->c);
    nmod_poly_mat_degree_matrix(dmat, mat);
    fmpz_mat_print_pretty(dmat);
    printf("\n");
    fmpz_mat_clear(dmat);
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
