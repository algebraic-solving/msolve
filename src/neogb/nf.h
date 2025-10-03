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


#ifndef GB_NF_H
#define GB_NF_H

#include "data.h"
#include "f4.h"

void get_normal_form_matrix(
        const bs_t * const tbr,
        ht_t * bht,
        const len_t start,
        md_t *st,
        ht_t **shtp,
        hi_t **hcmp,
        mat_t **matp
        );

bs_t *core_nf(
        bs_t *tbr,
        md_t *md,
        const exp_t * const mul,
        bs_t *bs,
        int32_t *errp
        );

int64_t export_nf(
        void *(*mallocp) (size_t),
        /* return values */
        int32_t *nf_ld,   /* basis load */
        int32_t **nf_len, /* length of each poly in basis */
        int32_t **nf_exp, /* basis exponent vectors */
        void **nf_cf,     /* coefficients of basis elements */
        /* input values */
        const int32_t nr_tbr_gens,
        const int32_t *tbr_lens,
        const int32_t *tbr_exps,
        const void *tbr_cfs,
        const int32_t nr_bs_gens,
        const int32_t *bs_lens,
        const int32_t *bs_exps,
        const void *bs_cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t bs_is_gb,
        const int32_t nr_threads,
        const int32_t info_level
        );
#endif
