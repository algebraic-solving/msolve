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

#ifndef NMOD_POLY_MAT_UTILS_H
#define NMOD_POLY_MAT_UTILS_H

/** \brief Basic routines for univariate polynomial matrices over `nmod`
 * \file nmod_poly_mat_utils.h
 * \date 2023-09-08
 *
 * This provides additional functions for FLINT's for matrices over the
 * univariate polynomials, with type `nmod_poly_mat` (coefficients in a finite
 * field Z/pZ for "small" p).
 *
 */

// FLINT includes
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_mat.h>
#include <flint/fmpz_mat.h>

// for nmod_mat_poly_t
#include "nmod_mat_poly.h"


#ifdef __cplusplus
extern "C" {
#endif

// TODO remove once using flint's comp instead
/** apply permutation to vector */
NMOD_POLY_MAT_INLINE void
apply_perm_to_vector(slong *res, const slong *initial_vect,
                          const slong *perm, slong length)
{
    for (slong i = 0; i < length; i++)
        res[perm[i]] = initial_vect[i];
}


/** computes x**(-d) (A*B mod x**h) */
void nmod_poly_mat_middle_product(nmod_poly_mat_t res,
                                  const nmod_poly_mat_t A,
                                  const nmod_poly_mat_t B,
                                  slong d,
                                  slong h);

#if __FLINT_VERSION < 3
/** Truncate `pmat` at order `len` */
void nmod_poly_mat_truncate(nmod_poly_mat_t pmat, long len);
#endif



/** @name Shifts (multiplication by powers of the variable)
 *
 * Left and right shift operations (X is the variable):
 *  - LeftShift by n means multiplication by X^n
 *  - RightShift by n means division by X^n
 *  - shift `n` has to be nonnegative in both cases.
 *
 *  In all the following functions which involve an OUT parameter (`svec` or
 *  `smat`), this parameter may alias the corresponding IN parameter (`pvec` or
 *  `pmat`).
 */
//@{

/** Computes the right `k`-shift `smat` of the polynomial matrix `pmat` */
void nmod_poly_mat_shift_right(nmod_poly_mat_t smat,
                               const nmod_poly_mat_t pmat,
                               slong k);

//@} // doxygen group: Shifts (multiplication by powers of the variable)

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM CONSTANT OR MATP                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Set from constant matrix and matrix polynomial.
 *
 * Provides `set` and `init_set`, from either a constant matrix,
 * or from a `nmod_mat_poly_t`.
 *
 * Two main representations are used for polynomial matrices:
 *    - a polynomial matrix stored as a matrix with polynomial entries, that
 *    is, of type `nmod_poly_mat_t`
 *    - a polynomial matrix stored as a polynomial with matrix coefficients,
 *    that is, of type `nmod_mat_poly_t`
 *
 * In the latter representation, the array of matrices has length
 * `degree(mat)+1`; in particular, the zero matrix may be represented by an
 * array of length 0. Note also the following requirement: in the second
 * representation, all the matrix coefficients have the same row and column
 * dimensions; this is currently assumed to hold, and is not checked by the
 * algorithms.
 *
 * The following functions perform conversions between these types, with two
 * variants: either the degree is deduced from the input, or a user-provided
 * truncation order is used.
 *
 */
//@{

/** Set from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (a nonnegative integer). */
// TODO benchmark and try variants if needed
FLINT_DLL void
nmod_poly_mat_set_trunc_from_mat_poly(nmod_poly_mat_t pmat,
                                      const nmod_mat_poly_t matp,
                                      slong order);

/** Set from polynomial with matrix coefficients `matp`. */
// TODO benchmark and try variants if needed
NMOD_POLY_MAT_INLINE void
nmod_poly_mat_set_from_mat_poly(nmod_poly_mat_t pmat,
                                const nmod_mat_poly_t matp)
{
    nmod_poly_mat_set_trunc_from_mat_poly(pmat, matp, matp->length);
}
//@} // doxygen group: Conversion from constant matrix and matrix polynomial


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* INPUT/OUTPUT                                               */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

#if __FLINT_VERSION < 3
/** Matrix pretty print to standard output */
void nmod_poly_mat_print(const nmod_poly_mat_t mat, const char * var);
#endif


/** Print the degree matrix, see @ref DegreeMatrix */
void nmod_poly_mat_degree_matrix_print_pretty(const nmod_poly_mat_t mat);




#ifdef __cplusplus
}
#endif

#endif // NMOD_POLY_MAT_UTILS_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
