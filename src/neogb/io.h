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

#ifndef GB_IO_H
#define GB_IO_H

#include "data.h"

void set_function_pointers(
        const stat_t *st
        );

int32_t check_and_set_meta_data(
        stat_t *st,
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

/* for normal form input data */
void import_julia_data_nf_ff_32(
        bs_t *tbr,
        ht_t *ht,
        stat_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

void import_julia_data_nf_qq(
        bs_t *tbr,
        ht_t *ht,
        stat_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

int32_t check_and_set_meta_data_trace(
        stat_t *st,
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
        const uint32_t prime_start,
        const int32_t nr_primes,
        const int32_t pbm_file,
        const int32_t info_level
        );
#endif
