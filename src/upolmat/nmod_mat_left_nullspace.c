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

#include "nmod_mat_extra.h"
#include <flint/perm.h>

slong nmod_mat_left_nullspace(nmod_mat_t X, const nmod_mat_t A)
{
    // helper lists of nonpivot|pivot columns of X
    slong * permutation = flint_malloc(A->r * sizeof(slong));

    // compute compact form of left nullspace
    nmod_mat_t Y;
    const slong nullity = nmod_mat_left_nullspace_compact(Y,permutation,A);

    // deduce nullspace
    const slong rank = A->r - nullity;
    nmod_mat_init(X, nullity, A->r, A->mod.n);
    for (slong i = 0; i < nullity; ++i)
        nmod_mat_entry(X, i, permutation[rank+i]) = 1;
    for (slong i = 0; i < nullity; ++i)
        for (slong j = 0; j < rank; ++j)
            nmod_mat_entry(X, i, permutation[j]) = nmod_mat_entry(Y, i, j);

    nmod_mat_clear(Y);
    flint_free(permutation);

    return nullity;
}

slong nmod_mat_left_nullspace_compact(
                                      nmod_mat_t X,
                                      slong * permutation,
                                      const nmod_mat_t A
                                     )
{
    // At <- transpose A
    nmod_mat_t At;
    nmod_mat_init(At, A->c, A->r, A->mod.n);
    nmod_mat_transpose(At,A);

    // Xt <- right kernel of At
    nmod_mat_t Xt;
    nmod_mat_init(Xt, At->c, At->c, At->mod.n);
    slong nullity = nmod_mat_nullspace(Xt, At);

    const slong rank = At->c - nullity;
    nmod_mat_init(X, nullity, rank, A->mod.n);

    if (nullity)
    {
        // find row rank profile of A and its complement, store them concatenated
        // in permutation (first rrp, then complement, each in increasing order)
        // -> recall: row rank profile = indices of nonpivots in Xt which is in
        // reduced column echelon form

        // search pivots first
        for (slong j = rank; j < At->c; ++j)
        {
            permutation[j] = Xt->r - 1;
            while (permutation[j] >= 0 && nmod_mat_entry(Xt, permutation[j], j-rank) == 0)
                --permutation[j];
            // Note: we should never arrive at (permutation[j] < 0)
        }

        slong r = 0;

        for (slong i = 0; i < permutation[rank]; i++, r++)
            permutation[r] = i;

        for (slong j = rank; j < At->c -1; ++j)
            for (slong i = permutation[j]+1; i < permutation[j+1]; i++, r++)
                permutation[r] = i;

        for (slong i = permutation[At->c -1]+1; i < Xt->r; i++, r++)
            permutation[r] = i;

        // extract dense part of kernel
        for (slong i = 0; i < nullity; ++i)
            for (slong j = 0; j < rank; ++j)
                nmod_mat_entry(X, i, j) = nmod_mat_entry(Xt, permutation[j], i);
    }
    else
#if __FLINT_VERSION >= 3 && __FLINT_VERSION_MINOR >= 1
        _perm_one(permutation, A->r);
#else
        _perm_set_one(permutation, A->r);
#endif

    // clean
    nmod_mat_clear(At);
    nmod_mat_clear(Xt);

    return nullity;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
