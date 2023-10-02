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

/**
 * \file nmod_poly_mat_approximant.h
 * \brief Minimal approximant bases
 *
 */

/** \file nmod_poly_mat_approximant.h
 * Definition (approximant basis).
 * -------------------------------
 * Consider:
 *   - an m x n matrix of univariate polynomials F,
 *   - an approximation order d (list of n positive integers).
 *
 * Then an approximant basis for (F,d) is a matrix over the univariate
 * polynomials whose rows form a basis for the following module:
 * { p in K[X]^{1 x m}  |  the column j of p F is 0 modulo X^{d[j]} }.
 * Note that such a matrix is square, m x m, and nonsingular. Its determinant
 * has the form X^k for some 0 <= k <= d[0] + d[1] + .. + d[n-1].
 */

/** \file nmod_poly_mat_approximant.h
 * Definition (shifted minimal approximant basis).
 * -----------------------------------------------
 * Starting from the definition of an approximant basis, consider further:
 *   - a degree shift s (a list of m integers).
 *
 * Then an approximant basis for (F,d) is said to be <em>a shift-minimal</em>
 * (resp. <em>a shift-ordered weak Popov</em>, resp. <em>the shift-Popov</em>)
 * approximant basis if it is in shift-reduced form (resp. in shift-ordered
 * weak Popov form, resp. in shift-Popov form). See nmod_poly_mat_forms.h
 * for definitions of these forms.
 */

/** \file nmod_poly_mat_approximant.h
 * Conventions.
 * ------------
 * Apart from the general interfaces (TODO) which offer the choice between left or
 * right approximants, all other functions compute left approximant bases
 * (approximants operate on the left of the matrix F; the basis elements are
 * the rows of the matrix).
 *
 * Most functions below use the following parameters.
 *
 * \param[out] appbas the output approximant basis (cannot alias `pmat`)
 * \param[in] pmat the input polynomial matrix (no restriction)
 * \param[in] order the input order (list of strictly positive integers, length must be the number of columns of `pmat`)
 * \param[in,out] shift in: the input shift; and out: the output shifted row degree of `appbas` (list of integers, length must be the number of rows of `pmat`)
 *
 * Note that the latter two restrictions on the lengths of the lists are
 * assuming left approximants; for right approximants, they are swapped.
 */

#ifndef NMOD_POLY_MAT_APPROXIMANT_H
#define NMOD_POLY_MAT_APPROXIMANT_H

#define PMBASIS_THRES 32

#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @name M-Basis algorithm (uniform approximant order)
 *
 * The core of these functions is implemented with `nmod_mat_poly_t` type,
 * for efficiency reasons. The bulk of the documentation can be found
 * in .nmod_mat_poly.h .
 *
 */
//@{

// TODO DOC
// appbas must be initialized with right dimensions
// shift=NULL for uniform shift
FLINT_DLL void
nmod_poly_mat_mbasis(nmod_poly_mat_t appbas,
                     slong * shift,
                     const nmod_poly_mat_t pmat,
                     ulong order);

//@} // doxygen group: M-Basis algorithm (uniform approximant order)


/** @name PM-Basis algorithm (uniform approximant order)
 * \anchor pmbasis
 *
 * These functions compute a `shift`-minimal ordered weak Popov approximant
 * basis for `(pmat,orders)` in the case where `orders` is given by a single
 * integer `orders = (order,...order)`. They use a divide and conquer approach,
 * computing a first basis at order `order/2`, finding the so-called _residual
 * matrix_, computing a second basis at order `order/2`, and deducing the
 * sought basis by multiplying the two obtained bases.
 *
 * The first recursive call returns an approximant basis `appbas1` such that
 * `appbas1*pmat = 0 mod x^{order/2}`, and the residual matrix has the same
 * dimensions as `pmat` and is defined by the matrix middle product
 * `(x^{-order/2} appbas1*pmat) mod x^{order/2}`.
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `appbas`, for the input shift. 
 *
 * This is inspired from the algorithm _pmbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014.
 */
//@{

// TODO DOC
// appbas must be initialized with right dimensions
// shift=NULL for uniform shift
FLINT_DLL void
nmod_poly_mat_pmbasis(nmod_poly_mat_t appbas,
                      slong * shift,
                      const nmod_poly_mat_t pmat,
                      slong order);

//@} // doxygen group: PM-Basis algorithm (uniform approximant order)


#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_APPROXIMANT_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
