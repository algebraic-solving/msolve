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
 * Mohab Safey El Din */


#include "data.h"

/* function pointers */
/* bs_t *(*initialize_basis)(
 *         const int32_t ngens
 *         ); */
void (*normalize_initial_basis)(
        bs_t *bs,
        const uint32_t fc
        );

int (*initial_input_cmp)(
        const void *a,
        const void *b,
        void *ht
        );

int (*initial_gens_cmp)(
        const void *a,
        const void *b,
        void *ht
        );

int (*monomial_cmp)(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        );

int (*spair_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

int (*hcm_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

/* linear algebra routines */
void (*sba_linear_algebra)(
        smat_t *smat,
        crit_t *syz,
        md_t *st,
        const ht_t * const ht
        );
void (*linear_algebra)(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        );

int (*application_linear_algebra)(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        );

void (*trace_linear_algebra)(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        );

void (* interreduce_matrix_rows)(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        int free_basis
        );

cf32_t *(*reduce_dense_row_by_old_pivots_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,
        const uint32_t fc
        );

hm_t *(*sba_reduce_dense_row_by_known_pivots_sparse_ff_32)(
        int64_t *dr,
        smat_t *smat,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t sm,      /* signature monomial of row reduced */
        const len_t si,     /* signature index of row reduced */
        const len_t ri,     /* index of row in matrix */
        md_t *st
        );

hm_t *(*reduce_dense_row_by_known_pivots_sparse_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        );

hm_t *(*trace_reduce_dense_row_by_known_pivots_sparse_ff_32)(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi,
        md_t *st
        );

cf32_t *(*reduce_dense_row_by_all_pivots_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const *pivs,
        cf32_t *const *dpivs,
        const uint32_t fc
        );


cf32_t *(*reduce_dense_row_by_dense_new_pivots_ff_32)(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        );
