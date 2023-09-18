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




/* The following is a copy-paste from a Flint file added in version 3 */




#if __FLINT_VERSION < 3

/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat_extra.h"

/** Permute rows of a matrix `mat` according to `perm_act`, and propagate the
 * action on `perm_store`.
 * That is, performs for each appropriate index `i`, the operations
 * `perm_store[i] <- perm_store[perm_act[i]]`
 * `rows[i] <- rows[perm_act[i]]` */
void
nmod_mat_permute_rows(nmod_mat_t mat, const slong * perm_act, slong * perm_store)
{
    slong i;
    mp_limb_t ** mat_tmp = (mp_limb_t **) flint_malloc(mat->r * sizeof(mp_limb_t *));

    /* perm_store[i] <- perm_store[perm_act[i]] */
    if (perm_store)
        _perm_compose(perm_store, perm_store, perm_act, mat->r);

    /* rows[i] <- rows[perm_act[i]]  */
    for (i = 0; i < mat->r; i++)
        mat_tmp[i] = mat->rows[perm_act[i]];
    for (i = 0; i < mat->r; i++)
        mat->rows[i] = mat_tmp[i];

    flint_free(mat_tmp);
}

#endif //  __FLINT_VERSION < 3

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
