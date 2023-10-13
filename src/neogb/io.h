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
        const md_t *st
        );

void set_ff_bits(md_t *st, int32_t fc);

void sort_terms_ff_8(
    cf8_t **cfp,
    hm_t **hmp,
    ht_t *ht
    );

void sort_terms_ff_16(
    cf16_t **cfp,
    hm_t **hmp,
    ht_t *ht
    );

void sort_terms_ff_32(
    cf32_t **cfp,
    hm_t **hmp,
    ht_t *ht
    );

void sort_terms_qq(
    mpz_t **cfp,
    hm_t **hmp,
    ht_t *ht
    );

void import_input_data(
        bs_t *bs,
        md_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs,
        const int *invalid_gens
        );

int32_t check_and_set_meta_data(
        md_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int *invalid_gens,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t nr_nf,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        );

int32_t check_and_set_meta_data_trace(
        md_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int *invalid_gens,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t nr_nf,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t reduce_gb,
        const uint32_t prime_start,
        const int32_t nr_primes,
        const int32_t pbm_file,
        const int32_t info_level
        );

/* for normal form input data */
void import_input_data_nf_ff_32(
        bs_t *tbr,
        ht_t *ht,
        md_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

void import_input_data_nf_ff_16(
                                bs_t *tbr,
                                ht_t *ht,
                                md_t *st,
                                const int32_t start,
                                const int32_t stop,
                                const int32_t *lens,
                                const int32_t *exps,
                                const void *vcfs
                                );

void import_input_data_nf_qq(
        bs_t *tbr,
        ht_t *ht,
        md_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

int validate_input_data(
        int **invalid_gensp,
        const void *cfs,
        const int32_t *lens,
        uint32_t *field_charp,
        int32_t *mon_orderp,
        int32_t *elim_block_lenp,
        int32_t *nr_varsp,
        int32_t *nr_gensp,
        int32_t *nr_nfp,
        int32_t *ht_sizep,
        int32_t *nr_threadsp,
        int32_t *max_nr_pairsp,
        int32_t *reset_htp,
        int32_t *la_optionp,
        int32_t *use_signaturesp,
        int32_t *reduce_gbp,
        int32_t *info_levelp
        );
#endif
