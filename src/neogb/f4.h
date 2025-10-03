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


#ifndef GB_F4_H
#define GB_F4_H

#include "data.h"

void free_f4_julia_result_data(
        void (*freep) (void *),
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,      /* coefficients of basis elements */
        const int64_t ngens,
        const int64_t field_char
        );

int64_t export_f4(
        void *(*mallocp) (size_t),
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        );

int64_t export_results_from_f4(
    /* return values */
    int32_t *bld,   /* basis load */
    int32_t **blen, /* length of each poly in basis */
    int32_t **bexp, /* basis exponent vectors */
    void **bcf,     /* coefficients of basis elements */
    void *(*mallocp) (size_t),
    bs_t **bsp,
    ht_t **bhtp,
    md_t **stp
    );

bs_t *core_f4(
        bs_t *gbs,
        md_t *gmd,
        int32_t *errp,
        const len_t fc
        );

bs_t *modular_f4(
        const bs_t * const ggb,       /* global basis */
        ht_t * gbht,                  /* global basis hash table, shared */
        md_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        );
#endif
