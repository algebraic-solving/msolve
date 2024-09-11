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

#include <stdlib.h> // for qsort
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/perm.h>
#include <flint/profiler.h>
#include "nmod_mat_poly.h"
#include "nmod_mat_extra.h"

/* type for stable sort while retaining the permutation */
typedef struct
{
    slong value;
    slong index;
} slong_pair;

/* comparator for quicksort, lexicographic (total) order to ensure stable sort */
static inline int _slong_pair_compare(const void * a, const void * b)
{
    slong_pair aa = * (const slong_pair *) a;
    slong_pair bb = * (const slong_pair *) b;
    if (aa.value == bb.value)
    {
        if (aa.index < bb.index)
            return -1;
        else if (aa.index > bb.index)
            return 1;
        else // aa.index == bb.index
            return 0;
    }
    else if (aa.value < bb.value)
        return -1;
    else // aa.value > bb.value
        return 1;
}

/** Creates a permutation from the sorting of a shift
 * After running this, perm is the unique list of integers which sorts
 * the pairs (shift,index) increasingly, i.e.
 * shift[perm[0]] <= shift[perm[1]] < ... < shift[perm[n-1]]
 *
 * \param perm permutation (list of integers), length n
 * \param shift list of integer to be sorted nondecreasingly, length n
 * \param pair_tmp temporary storage, length n
 * \param n length
 *
 */
static inline void _find_shift_permutation(slong * perm,
                                           const slong * shift,
                                           slong n,
                                           slong_pair * pair_tmp)
{
    for (slong i = 0; i < n; i++)
    {
        pair_tmp[i].value = shift[i];
        pair_tmp[i].index = i;
    }

    qsort(pair_tmp, n, sizeof(slong_pair), _slong_pair_compare);

    for (slong i = 0; i < n; i++)
        perm[i] = pair_tmp[i].index;
}

void nmod_mat_poly_mbasis(nmod_mat_poly_t appbas,
                          slong * shift,
                          const nmod_mat_poly_t matp,
                          slong order)
{
    // dimensions of input matrix
    const slong m = matp->r;
    const slong n = matp->c;

    // initialize output approximant basis with identity
    // except when matp == 0: return appbas = x**order * identity
    if (nmod_mat_poly_is_zero(matp))
    {
        nmod_mat_poly_fit_length(appbas, order+1);
        _nmod_mat_poly_set_length(appbas, order+1);
        nmod_mat_one(appbas->coeffs + order);
        return;
    }
    else
        nmod_mat_poly_one(appbas);
    // -> ensures matp->length > 0 in what follows

    // residual matrix: m x n constant matrix, next coefficient of appbas *
    // matp to be annihilated
    nmod_mat_t res;
    nmod_mat_init(res, m, n, matp->mod.n);
    // temporary matrix used during the computation of residuals
    nmod_mat_t res_tmp;
    nmod_mat_init(res_tmp, m, n, res->mod.n);

    // will hold the stable permutation making the shift nondecreasing
    slong * perm = _perm_init(m);
    slong_pair * pair_tmp = (slong_pair *) flint_malloc(m * sizeof(slong_pair));

    // to store the left constant nullspace bases, their ranks,
    // and corresponding pivot information
    slong nullity;
    slong * pivots = (slong *) flint_malloc(m * sizeof(slong));
    nmod_mat_t nsbas;

    // Note: iterations will guarantee that `shift` (initially, this is the
    // input shift) holds the `shift`-shifted row degree of appbas; in
    // particular this is the case when the algorithm returns
    for (slong ord = 0; ord < order; ++ord)
    {
#ifdef MBASIS_DEBUG
printf("\n\nDEBUG -- Start iteration %ld\n", ord);
#endif // MBASIS_DEBUG

        // Letting s0 be the original input shift:
        // start of loop: appbas is an `s0`-ordered weak Popov approximant
        //    basis at order `ord` for `pmat`, and shift = rdeg_s0(appbas)
        // end of loop: appbas is an `s0`-ordered weak Popov approximant
        //    basis at order `ord+1` for `pmat`, and shift = rdeg_s0(appbas)
        if (ord == 0) // initially res = coeffs_matp[0] (note matp->length > 0 here)
            nmod_mat_set(res, matp->coeffs + 0);
        else // res = coefficient of degree ord in appbas*pmat
            nmod_mat_poly_mul_coeff(res, appbas, matp, ord);
        // TODO try another variant with res_tmp to save memory?

        // compute stable permutation which makes the shift nondecreasing
        // --> we will need to permute things, to take into account the
        // "priority" indicated by the shift
        _find_shift_permutation(perm, shift, m, pair_tmp);

        // permute rows of the residual accordingly
        nmod_mat_permute_rows(res, perm, NULL);

        // find the left nullspace basis of the (permuted) residual, in reduced
        // row echelon form and compact storage (see the documentation)
        if (ord >= 1) // nsbas should be uninitialized for nullspace
            nmod_mat_clear(nsbas);
        nullity = nmod_mat_left_nullspace_compact(nsbas,pivots,res);

        if (nullity==0)
        {
            // Exceptional case: the residual matrix has empty left nullspace
            // --> no need to compute more: the final basis is X^(order-ord)*appbas
            break;
        }

        else if (nullity < m)
        {
            // nothing to do if nullity == m, another exceptional case:
            // residual coeff was zero, and nullspace 'nsbas' is identity
            // --> approximant basis is already correct for this order, no need to
            // change it or to change shift, go to next iteration

            // here, we are in the "usual" case, where the left nullspace of the
            // residual has no special shape

            // things will be easier in the permuted world, gathering pivot
            // rows and non-pivot rows
            // --> we are currently working on perm * appbas
            //         i.e. perm*appbas[i] = appbas[perm[i]]
            // --> we will add to this the pivot permutation:
            //         i.e. pivots*perm*appbas[i] = appbas[perm[pivots[i]]]

            // transform pivots[i] into perm[pivots[i]]
            _perm_compose(pivots, perm, pivots, m);

            // Update the shift
            // --> in the permuted world, we should add 1 to each entry
            // corresponding to a nullspace column which contains a pivot
            // (don't move below! pivots will be modified)
            for (slong i = 0; i < m-nullity; i++)
                shift[pivots[i]] += 1;

            // Update the approximant basis: 1/ go to permuted world
            nmod_mat_poly_permute_rows(appbas, pivots, NULL);

            // Update the approximant basis: 2/ perform constant update
            // Submatrix of rows corresponding to pivots are replaced by
            // nsbas*appbas (note: these rows currently have degree
            // at most deg(appbas))
            // FIXME possible small improvement for uniform shift: these rows
            // have degree less than deg(appbas), in this case (and deg(appbas)
            // is reached on the diagonal, among the pivot degrees)
            // TODO see if a convenient way to make this matrix fixed at the beginning;
            // note nullity is most often unchanged along iterations
            // TODO generally try to improve memory handling in here
            nmod_mat_t ns_app, app_win_mul, app_win_add;
            nmod_mat_init(ns_app, nullity, appbas->c, appbas->mod.n);
            for (slong d = 0; d < appbas->length; ++d)
            {
                nmod_mat_window_init(app_win_mul, appbas->coeffs + d, 0, 0, m-nullity, m);
                nmod_mat_window_init(app_win_add, appbas->coeffs + d, m-nullity, 0, m, m);
                nmod_mat_mul(ns_app, nsbas, app_win_mul);
                nmod_mat_add(app_win_add, app_win_add, ns_app);
            }
            nmod_mat_clear(ns_app);
            nmod_mat_clear(app_win_mul);
            nmod_mat_clear(app_win_add);

            // Update the approximant basis: 3/ left-shift relevant rows by 1
            // increase length of appbas if necessary
            for (slong i = 0; i < m-nullity; ++i)
                if (!_nmod_vec_is_zero(appbas->coeffs[appbas->length-1].rows[i], m))
                {
                    nmod_mat_poly_fit_length(appbas, appbas->length+1);
                    _nmod_mat_poly_set_length(appbas, appbas->length+1);
                    i = m - nullity;
                }

            // rows not corresponding to pivots are multiplied by X
            // note: these rows currently have length strictly less than len(appbas)
            mp_limb_t ** save_zero_rows = (mp_limb_t **) flint_malloc((m-nullity) * sizeof(mp_limb_t *));
            for (slong i = 0; i < m-nullity; ++i)
            {
                save_zero_rows[i] = appbas->coeffs[appbas->length-1].rows[i];
                appbas->coeffs[appbas->length-1].rows[i] = appbas->coeffs[appbas->length-2].rows[i];
                // note: arrived at this stage, appbas->length must be >= 2
            }
            for (slong d = appbas->length-2; d > 0; d--)
                for (slong i = 0; i < m-nullity; ++i)
                    appbas->coeffs[d].rows[i] = appbas->coeffs[d-1].rows[i];
            for (slong i = 0; i < m-nullity; ++i)
                //_nmod_vec_zero(appbas->coeffs[0].rows[i], m);
                appbas->coeffs[0].rows[i] = save_zero_rows[i];

            // Update the approximant basis: 4/ come back from permuted world
            _perm_inv(pivots, pivots, m);
            nmod_mat_poly_permute_rows(appbas, pivots, NULL);
        }
    }

    nmod_mat_clear(res);
    nmod_mat_clear(res_tmp);
    _perm_clear(perm);
    flint_free(pair_tmp);
    flint_free(pivots);
    nmod_mat_clear(nsbas);
    return;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
