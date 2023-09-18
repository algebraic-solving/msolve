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

#include "nmod_mat_poly.h"  // truncate, shift, middle_product
#include "nmod_poly_mat_utils.h"  // truncate, shift, middle_product
#include "nmod_poly_mat_pmbasis.h"



void nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                          slong * shift,
                          const nmod_poly_mat_t pmat,
                          ulong order)
{
    nmod_mat_poly_t app, matp;
    nmod_mat_poly_init(matp, pmat->r, pmat->c, pmat->modulus);
    // TODO improve: set init
    nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, order);
    nmod_mat_poly_init(app, pmat->r, pmat->r, pmat->modulus);
    nmod_mat_poly_mbasis(app, shift, matp, order);
    // TODO improve: set init
    nmod_poly_mat_set_from_mat_poly(appbas, app);
    nmod_mat_poly_clear(matp);
    nmod_mat_poly_clear(app);
}


void nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                           slong * shift,
                           const nmod_poly_mat_t pmat,
                           slong order)
{
    if (order <= PMBASIS_THRES)
    {
        nmod_poly_mat_mbasis(appbas, shift, pmat, order);
        return;
    }

    const long order1 = order>>1;
    const long order2 = order - order1;
    nmod_poly_mat_t appbas2, residual;

    nmod_poly_mat_init(appbas2, pmat->r, pmat->r, pmat->modulus);
    nmod_poly_mat_init(residual, pmat->r, pmat->c, pmat->modulus);

    nmod_poly_mat_pmbasis(appbas, shift, pmat, order1);

    nmod_poly_mat_middle_product(residual, appbas, pmat, order1, order);

    nmod_poly_mat_pmbasis(appbas2, shift, residual, order2);

    nmod_poly_mat_mul(appbas, appbas2, appbas);

    nmod_poly_mat_clear(appbas2);
    nmod_poly_mat_clear(residual);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
