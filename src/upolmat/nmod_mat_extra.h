/** \file nmod_mat_extra.h
 *
 */

#ifndef __NMOD_MAT_EXTRA__H
#define __NMOD_MAT_EXTRA__H

#include <flint/flint.h>
#include <flint/perm.h>
#include <flint/nmod_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Left nullspace of A.
 *
 *  Computes a basis X for the left nullspace of A, in reduced row echelon form
 *  with pivots being the rightmost nonzero entries. X should be given
 *  uninitialized, it will be initialized during the call with the right number
 *  of rows.
 *
 *  This first calls nmod_mat_left_nullspace_compact(), and expands the compact
 *  nullspace representation given by this call into the complete dense
 *  nullspace representation.
 *
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \param [in] A input matrix
 * \return nullity of A (i.e. rank of X)
 *
 * @see nmod_mat_left_nullspace_compact
 */
FLINT_DLL slong nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A);

/** Left nullspace of A in compact form.
 *
 *  Computes a basis X for the left nullspace of A, in reduced row echelon form
 *  with pivots being the rightmost nonzero entries. Only the nonpivot columns
 *  of X are stored, in the order they appear in the nullspace basis. The list
 *  permutation contains the concatenation of two lists, each in increasing
 *  order: the positions of the columns without pivots in the rref nullspace
 *  basis, and the positions of the columns with pivots in the rref nullspace
 *  basis (the former being also the row rank profile of A).
 *
 * \param [out] X matrix where the nullspace will be stored (uninitialized)
 * \param [out] permutation list, allocated with A->r elements
 * \param [in] A input matrix
 * \return nullity of A (i.e. rank of X)
 *
 * \todo efficiency is probably not best for small matrices: this uses Flint's
 * right nullspace and matrix transposition
 */
FLINT_DLL slong nmod_mat_left_nullspace_compact(
                                                nmod_mat_t X,
                                                slong * permutation,
                                                const nmod_mat_t A
                                                );

#ifdef __cplusplus
}
#endif

#endif  // __NMOD_MAT_EXTRA__H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
