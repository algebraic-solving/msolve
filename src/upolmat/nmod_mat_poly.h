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
 * \file nmod_mat_poly.h
 * \brief Main header for univariate polynomial with matrix coefficients modulo word-size prime
 * \version 0.0
 * \date 2023-02-13
 *
 * This is the main header for functions for univariate polynomials with
 * coefficients that are matrices over a finite field Z/pZ for "small" p
 * (word-size, so that we use a representation with flint's ``nmod``). Only
 * operations for which this `nmod_mat_poly` format is best suited are
 * provided. For other operations it is suggested to rely the conversions to and
 * from `nmod_poly_mat` format and use the relevant functions for that format.
 *
 * \todo benchmark performance
 * \todo test for memory leaks
 * \todo Note: all parameters are supposed init
 *
 * \todo transpose, swap/permute rows/cols
 * \todo set, init_set, swap
 * \todo init2_mod for init from mod with given alloc
 * \todo set_trunc, set_shift, attach_trunc, attach_shift
 * \todo windows?
 */

#ifdef NMOD_MAT_POLY_INLINES_C
#define NMOD_MAT_POLY_INLINE FLINT_DLL
#else
#define NMOD_MAT_POLY_INLINE static __inline__
#endif

#ifndef NMOD_MAT_POLY_H
#define NMOD_MAT_POLY_H

#include <flint/perm.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------*/
/* struct, typedefs                                           */
/*------------------------------------------------------------*/

/** Struct for matrix polynomials.
 *
 * Storage is a dynamic array of matrices `nmod_mat`. The maximum number of
 * coefficients is `alloc` and the actual number of coefficients (which is the
 * degree plus 1) is `length`. The number of rows and columns are `r` and `c`.
 * The modulus is stored as an `nmod_t`. In the provided functions, e.g. for
 * modifying a coefficient, it is not checked that the dimensions of the
 * modified coefficient are indeed `r x c`.
 */
typedef struct
{
    nmod_mat_struct * coeffs; /**< array of coefficients */
    slong alloc;              /**< allocated length */
    slong length;             /**< actual length */
    slong r;                  /**< number of rows */
    slong c;                  /**< number of columns */
    nmod_t mod;               /**< modulus */
} nmod_mat_poly_struct;

/** nmod_mat_poly_t allows reference-like semantics for nmod_mat_poly_struct */
typedef nmod_mat_poly_struct nmod_mat_poly_t[1];

/*------------------------------------------------------------*/
/* memory management                                          */
/*------------------------------------------------------------*/

/** @name Memory management for matrix polynomials
 */
//@{

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. */
FLINT_DLL void nmod_mat_poly_init(nmod_mat_poly_t matp,
                                  slong r, slong c,
                                  mp_limb_t n);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. The caller supplies a precomputed inverse limb generated by
 * `n_preinvert_limb()`. */
FLINT_DLL void nmod_mat_poly_init_preinv(nmod_mat_poly_t matp,
                                         slong r, slong c,
                                         mp_limb_t n, mp_limb_t ninv);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients modulo
 * `n`. Up to alloc coefficients may be stored in `matp`. Implementation note:
 * the `alloc` matrix coefficients are not initialized. */
FLINT_DLL void nmod_mat_poly_init2(nmod_mat_poly_t matp,
                                   slong r, slong c,
                                   mp_limb_t n,
                                   slong alloc);

/** Initialises `matp`. It will have dimensions `r x c` and coefficients
 * modulo~`n`. The caller supplies a precomputed inverse limb generated by
 * n_preinvert_limb(). Up to alloc coefficients may be stored in `matp`.
 * Implementation note: the `alloc` matrix coefficients are not initialized. */
FLINT_DLL void nmod_mat_poly_init2_preinv(nmod_mat_poly_t matp,
                                          slong r, slong c,
                                          mp_limb_t n, mp_limb_t ninv,
                                          slong alloc);

/** Sets `matp1` to equal `matp2`. It is assumed that these matrices
 * have identical dimensions and modulus. */
FLINT_DLL void
nmod_mat_poly_set(nmod_mat_poly_t matp1, const nmod_mat_poly_t matp2);


/** Initialises `matp` using an already initialised modulus `mod`. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_init_mod(nmod_mat_poly_t matp,
                       slong r, slong c,
                       const nmod_t mod)
{
    matp->coeffs = NULL;
    matp->alloc = 0;
    matp->length = 0;
    matp->r = r;
    matp->c = c;
    matp->mod = mod;
}

/** Set the modulus to `mod`. */
NMOD_MAT_POLY_INLINE
void nmod_mat_poly_set_mod(nmod_mat_poly_t matp, const nmod_t mod)
{
    matp->mod = mod;
}

/** Reallocates `matp` to the given length `alloc`. If the current `length` is
 * more than `alloc`, the polynomial is truncated and normalised. If `alloc` is
 * zero, the polynomial is cleared. */
FLINT_DLL void nmod_mat_poly_realloc(nmod_mat_poly_t matp, slong alloc);

/** Clears `matp` and releases any memory it used. The matrix polynomial cannot be
used again until it is initialised. */
FLINT_DLL void nmod_mat_poly_clear(nmod_mat_poly_t matp);

/** Ensures `matp` has space for at least `alloc` coefficients. This function
 * only ever grows the allocated space, so no data loss can occur. */
FLINT_DLL void nmod_mat_poly_fit_length(nmod_mat_poly_t matp, slong alloc);

/** Sets the length of `matp` to `length`. If `length < matp->length` then the
 * matrix coefficients of `matp` beyond `length` are cleared; otherwise, the
 * matrix coefficients of `matp` are initialized up to `length`. Note: `matp`
 * may not be normalized after this. The provided `length` must be less than
 * `matp->alloc`. */
NMOD_MAT_POLY_INLINE void
_nmod_mat_poly_set_length(nmod_mat_poly_t matp, slong length)
{
    if (matp->length > length)
        for (slong i = length; i < matp->length; i++)
            nmod_mat_clear(matp->coeffs + i);
    else
        for (slong i = matp->length; i < length; i++)
            nmod_mat_init(matp->coeffs + i, matp->r, matp->c, matp->mod.n);
    matp->length = length;
}


/** Normalises a matrix polynomial `matp` so that the top coefficient, if there
 * is one at all, is not zero. */
NMOD_MAT_POLY_INLINE void
_nmod_mat_poly_normalise(nmod_mat_poly_t matp)
{
    while (matp->length && nmod_mat_is_zero(matp->coeffs + matp->length - 1))
    {
        nmod_mat_clear(matp->coeffs + matp->length - 1);
        matp->length--;
    }
}

//@} // doxygen group:  Memory management for matrix polynomials

/*------------------------------------------------------------*/
/* Zero and Identity                                          */
/*------------------------------------------------------------*/

/** @name Zero and Identity
 * Functions to set an already initialized matrix polynomial `matp` to zero or
 * to one (this is the identity matrix if `matp` is square; in general this
 * puts `1` on the main diagonal and `0` elsewhere). Functions to test whether
 * a matrix polynomial is zero or one.
 */
//@{

/** Sets `matp` to the zero matrix polynomial */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_zero(nmod_mat_poly_t matp)
{
   _nmod_mat_poly_set_length(matp, 0);
}

/** \def nmod_mat_poly_is_zero(matp)
 * Tests whether `matp` is the zero matrix polynomial
 */
#define nmod_mat_poly_is_zero(matp) \
    ((matp)->length == 0)

/** Sets `matp` to the identity matrix polynomial */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_one(nmod_mat_poly_t matp)
{
    nmod_mat_poly_fit_length(matp, 1);
    _nmod_mat_poly_set_length(matp, 1);
    nmod_mat_one(matp->coeffs + 0);
}

/** Tests whether `matp` is one (i.e., if square, the identity matrix
 * polynomial; otherwise `1` on the main diagonal and `0` elsewhere). */
NMOD_MAT_POLY_INLINE int
nmod_mat_poly_is_one(const nmod_mat_poly_t matp)
{
    return (matp->length) == 1 && (nmod_mat_is_one(matp->coeffs + 0));
}

//@} // doxygen group:  Zero and Identity


/*------------------------------------------------------------*/
/* Accessing struct info and coefficients                     */
/*------------------------------------------------------------*/

/** @name Accessing struct info and matrix coefficients
 * Get the number of rows, the number of columns, the length, or the degree of
 * this matrix polynomial. Get a reference to the leading coefficient matrix,
 * or to the coefficient of degree `k`, or to the entry `i,j` of the
 * coefficient of degree `k`. Get or set the latter entry.
 */
//@{

/** Returns the number of rows of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_nrows(const nmod_mat_poly_t matp)
{
    return matp->r;
}

/** Returns the number of cols of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_ncols(const nmod_mat_poly_t matp)
{
    return matp->c;
}

/** Returns the length of `matp`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_length(const nmod_mat_poly_t matp)
{
    return matp->length;
}

/** Returns the degree of `matp`. By convention the zero matrix polynomial has
 * degree `-1`. */
NMOD_MAT_POLY_INLINE slong
nmod_mat_poly_degree(const nmod_mat_poly_t matp)
{
    return matp->length - 1;
}

/** \def nmod_mat_poly_coeff(matp, k)
 * Returns a reference to the coefficient of degree `k` in the matrix
 * polynomial `matp`. This function is provided so that individual coefficients
 * can be accessed and operated on by functions in the `nmod_mat` module. This
 * function does not make a copy of the data, but returns a reference
 * `nmod_mat_struct *` to the actual coefficient. Returns `NULL` when `k`
 * exceeds the degree of the matrix polynomial.
 */
#define nmod_mat_poly_coeff(matp, k) \
    ((k) < (matp)->length ? (matp)->coeffs + (k) : NULL)

/** Get the coefficient of degree `k` in the matrix polynomial `matp`. Zeroes
 * the output matrix when `k` exceeds the degree of the matrix polynomial. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_get_coeff(nmod_mat_t coeff, const nmod_mat_poly_t matp, slong k)
{
    if (k < matp->length)
        nmod_mat_set(coeff, matp->coeffs + k);
    else
        nmod_mat_zero(coeff);
}

/** \def nmod_mat_poly_lead(const nmod_mat_poly_t poly)
 * Returns a reference to the leading coefficient of the matrix polynomial, as
 * an `nmod_mat_struct *`. This function is provided so that the leading
 * coefficient can be easily accessed and operated on by functions in the
 * `nmod_mat` module. This function does not make a copy of the data, but
 * returns a reference to the actual coefficient.  Returns `NULL` when the
 * polynomial is zero.
 */
#define nmod_mat_poly_lead(matp) \
    ((matp)->length ? (matp)->coeffs + (matp)->length - 1 : NULL)

/** \def nmod_mat_poly_entry(matp,k,i,j)
 * Directly accesses the entry in the coefficient of `matp` of degree `k`, in
 * row `i` and column `j` (indexed from zero). No bounds checking is performed.
 * This macro can be used both for reading and writing coefficients.
 */
#define nmod_mat_poly_entry(matp,k,i,j) \
    (((matp)->coeffs + (k))->rows[(i)][(j)])

/** Get the entry at row `i` and column `j` in the coefficient of
 * degree `k` of the matrix polynomial `matp`. */
NMOD_MAT_POLY_INLINE mp_limb_t
nmod_mat_poly_get_entry(const nmod_mat_poly_t matp,
                        slong k, slong i, slong j)
{
   return (matp->coeffs + k)->rows[i][j];
}

/** Return a pointer to the entry at row `i` and column `j` of the coefficient
 * of degree `k` of the matrix polynomial `matp` */
NMOD_MAT_POLY_INLINE mp_limb_t *
nmod_mat_poly_entry_ptr(const nmod_mat_poly_t matp,
                        slong k, slong i, slong j)
{
   return (matp->coeffs + k)->rows[i] + j;
}

/** Set to `x` the entry at row `i` and column `j` in the coefficient of degree
 * `k` of the matrix polynomial `matp`. */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_set_entry(nmod_mat_poly_t matp,
                        slong k, slong i, slong j,
                        mp_limb_t x)
{
    nmod_mat_poly_entry(matp, k, i, j) = x;
}

//@} // doxygen group:  Accessing struct info and matrix coefficients


/*------------------------------------------------------------*/
/* Truncate, Shift, Reverse, Permute                          */
/*------------------------------------------------------------*/

/** @name Truncate, shift, reverse, permute
 * \todo TODO
 */
//@{

/** Tests whether `matp` is a constant matrix, that is, of degree 0 */
NMOD_MAT_POLY_INLINE int
nmod_mat_poly_is_constant(const nmod_mat_poly_t matp)
{
    return matp->length == 0;
}

/** Truncates `matp` to the given `order` and normalises it. If `order` is
 * greater than or equal to the current length of `matp`, then nothing happens.
 */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_truncate(nmod_mat_poly_t matp, slong order)
{
    if (matp->length > order)
    {
        for (slong i = order; i < matp->length; i++)
            nmod_mat_clear(matp->coeffs + i);
        matp->length = order;
        _nmod_mat_poly_normalise(matp);
    }
}

/** Sets `(smatp, len + n)` to `(matp, len)` shifted left by `n` coefficients.
 * Inserts zero coefficients at the lower end. Assumes that `len` and `n`
 are positive, and that `smatp` fits `len + n` elements. Supports aliasing
 between res and poly. */
FLINT_DLL void
_nmod_mat_poly_shift_left(nmod_mat_struct * smatp,
                          const nmod_mat_struct * matp,
                          slong len,
                          slong n);

// TODO
//FLINT_DLL void
//_nmod_mat_poly_shift_right(nmod_mat_poly_t smatp,
//                          const nmod_mat_poly_t matp,
//                          slong len,
//                          slong n);

/** Sets `smatp` to `matp` shifted left by `n` coeffs. Zero coefficients are
 * inserted. */
FLINT_DLL void
nmod_mat_poly_shift_left(nmod_mat_poly_t smatp,
                         const nmod_mat_poly_t matp,
                         slong n);

// TODO
//FLINT_DLL void
//nmod_mat_poly_shift_right(nmod_mat_poly_t smatp,
//                          const nmod_mat_poly_t matp,
//                          slong n);
//

/** Permute rows of a matrix polynomial `matp` according to `perm_act`, and
 * propagate the action on `perm_store`.
 * That is, performs for each appropriate index `i`, the operations
 * `perm_store[i] <- perm_store[perm_act[i]]`
 * `rows[i] <- rows[perm_act[i]]` */
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_permute_rows(nmod_mat_poly_t matp,
                           const slong * perm_act,
                           slong * perm_store)
{
    slong i;
    mp_limb_t ** mat_tmp = flint_malloc(matp->r * sizeof(mp_limb_t *));

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, matp->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (slong k = 0; k < matp->length; k++)
    {
        for (i = 0; i < matp->r; i++)
            mat_tmp[i] = matp->coeffs[k].rows[perm_act[i]];
        for (i = 0; i < matp->r; i++)
            matp->coeffs[k].rows[i] = mat_tmp[i];
    }

    flint_free(mat_tmp);
}

//@} // doxygen group:  Truncate, shift, reverse, permute


/*------------------------------------------------------------*/
/* Intput / Output                                            */
/*------------------------------------------------------------*/


/** @name Input / Output
 * Printing, writing to a file, reading from a file.
 * \todo
 */
//@{

/** Basic print to standard output */
FLINT_DLL void nmod_mat_poly_print(const nmod_mat_poly_t matp);

/** Pretty print to standard output */
FLINT_DLL void nmod_mat_poly_print_pretty(const nmod_mat_poly_t matp);

//@} // doxygen group:  Input / Output

/*------------------------------------------------------------*/
/* Basic arithmetic                                           */
/*------------------------------------------------------------*/

/** @name Basic arithmetic
 * \todo doc
 */
//@{

/** Compute the coefficient of degree `k` in the product of the two matrix
 * polynomials. Precisely, if `mat1` is some matrix `A = sum^{deg_A}_{i=0} a_i
 * x^i` and `mat2` is some matrix `B = sum^{deg_B}_{i=0} b_i x^i`, then `coeff`
 * will get the coefficient `k` of the product `AB`, which is `c_k =
 * sum^{i=0}_{k} a_i b_{k-i}`. It is not checked that dimensions or moduli are
 * compatible. The output `coeff` must already be initialized with dimensions
 * `mat1->r x mat2->c`.
 */
FLINT_DLL void
nmod_mat_poly_mul_coeff(nmod_mat_t coeff,
                        const nmod_mat_poly_t mat1,
                        const nmod_mat_poly_t mat2,
                        slong k);

//@} // doxygen group:  Basic arithmetic


/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* SET FROM CONSTANT OR PMAT                                  */
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/

/** @name Set from constant matrix and polynomial matrix.
 *
 * Provides `set` and `init_set`, from either a constant matrix,
 * or from a `nmod_poly_mat_t`.
 *
 * Two main representations are used for polynomial matrices:
 *    - a polynomial matrix stored as a matrix with polynomial entries, that
 *    is, of type `nmod_poly_mat_t`
 *    - a polynomial matrix stored as a polynomial with matrix coefficients,
 *    that is, of type `nmod_mat_poly_t`
 *
 * The following functions perform conversions between these types, with two
 * variants: either the degree is deduced from the input, or a user-provided
 * truncation order is used.
 *
 * \todo init-set variants
 */
//@{

/** Initialize a polynomial matrix `pmat` of the same dimensions and modulus as
 * a constant matrix `mat`, and set its constant coefficient to `mat`.
 * \todo to be written
 * careful: probably better to use nmod_poly_init with preinv...
 * but this should have been done in Flint's native poly_mat_init too?
 **/
FLINT_DLL void
nmod_mat_poly_init_set_from_nmod_mat(nmod_mat_poly_t matp,
                                     const nmod_mat_t cmat);

/** Set the polynomial matrix `pmat` to be a constant polynomial matrix whose
 * constant coefficient is a copy of `cmat`. This assumes `pmat` is already
 * initialized with the same modulus and dimensions as `cmat`.  */
FLINT_DLL void
nmod_mat_poly_set_from_nmod_mat(nmod_mat_poly_t matp, const nmod_mat_t cmat);

/** Set from polynomial with matrix coefficients `matp`, truncated at the
 * specified `order` (a nonnegative integer). */
FLINT_DLL void
nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                      const nmod_poly_mat_t pmat,
                                      slong order);

/** Set from polynomial with matrix coefficients `matp`. */
// TODO benchmark and try variants if needed
NMOD_MAT_POLY_INLINE void
nmod_mat_poly_set_from_poly_mat(nmod_mat_poly_t matp, const nmod_poly_mat_t pmat)
{
    nmod_mat_poly_set_trunc_from_poly_mat(matp, pmat, nmod_poly_mat_max_length(pmat));
}

//@} // doxygen group: Conversion from constant matrix and polynomial matrix


/** @name M-Basis algorithm (uniform approximant order)
 * \anchor mbasis
 *
 * See .nmod_poly_mat_approximant.h for general definitions and conventions
 * about approximant bases in PML.
 *
 * The functions here compute a `shift`-minimal ordered weak Popov approximant
 * basis for `(pmat,orders)` in the case where `orders` is given by a single
 * integer `orders = (order,...order)`. They iterate from `1` to `order`,
 * computing at each step a basis at order `1` (see @ref mbasis1) and using it
 * to update the output `appbas`, the so-called _residual matrix_, and the
 * considered shift. At step `d`, we have `appbas*pmat = 0 mod x^{d-1}`, and we
 * want to update `appbas` so that this becomes zero modulo `x^d`.
 *
 * In this context, the residual matrix is a constant matrix with the same
 * dimensions as `pmat` which, at the iteration `d`, is equal to the
 * coefficient of degree `d` of `appbas*pmat` (the coefficients of lower degree
 * being already zero).
 *
 * At the end of the computation, the vector `shift` contains the shifted row
 * degree of `appbas`, for the input shift. 
 *
 * This is inspired from the algorithm _mbasis_ described in
 *  - P. Giorgi, C.-P. Jeannerod, G. Villard. Proceeding ISSAC 2003,
 *  - P. Giorgi, R. Lebreton. Proceedings ISSAC 2014.
 */
//@{

/** Main mbasis function. */
FLINT_DLL void
nmod_mat_poly_mbasis(nmod_mat_poly_t appbas,
                     slong * shift,
                     const nmod_mat_poly_t matp,
                     slong order);

//@} // doxygen group: M-Basis algorithm (uniform approximant order)

#ifdef __cplusplus
}
#endif

#endif /* NMOD_MAT_POLY_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
