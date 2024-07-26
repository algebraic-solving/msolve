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

#ifndef MSOLVE_MSOLVE_H
#define MSOLVE_MSOLVE_H

#include "msolve-data.h"

int msolve_trace_qq(
        mpz_param_t *mpz_paramp,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t unstable_staircase,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t elim_block_len,
        int32_t reset_ht,
        int32_t la_option,
        int32_t use_signatures,
        int32_t lift_matrix,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        files_gb *files,
        int
        );

int msolve_probabilistic_qq(
        mpz_param_t mp_param,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t elim_block_len,
        int32_t reset_ht,
        int32_t la_option,
        int32_t use_signatures,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        files_gb *files,
        int round
        );

#if 0
int msolve_qq(
        mpz_param_t mp_param,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t elim_block_len,
        int32_t reset_ht,
        int32_t la_option,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        files_gb *files,
        int
        );
#endif

int real_msolve_qq(
        mpz_param_t *mpz_paramp,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        long *nb_real_roots_ptr,
        interval **real_roots_ptr,
        real_point_t **real_pts_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t unstable_staircase,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t elim_block_len,
        int32_t reset_ht,
        int32_t la_option,
        int32_t use_signatures,
        int32_t lift_matrix,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        int32_t precision,
        files_gb *files,
        int,
        int32_t
        );

int core_msolve(
        int32_t la_option,
        int32_t use_signatures,
        int32_t nr_threads,
        int32_t info_level,
        int32_t initial_hts,
        int32_t max_pairs,
        int32_t elim_block_len,
        int32_t update_ht,
        int32_t generate_pbm,
        int32_t reduce_gb,
        int32_t print_gb,
        int32_t truncate_lifting,
        int32_t get_param,
        int32_t genericity_handling,
        int32_t unstable_staircase,
        int32_t saturate,
        int32_t colon,
        int32_t normal_form,
        int32_t normal_form_matrix,
        int32_t is_gb,
        int32_t lift_matrix,
        int32_t precision,
        files_gb *files,
        data_gens_ff_t *gens,
        param_t **paramp,
        mpz_param_t *mpz_paramp,
        long *nb_real_roots_ptr,
        interval **real_roots_ptr,
        real_point_t **real_pts_ptr
        );

void msolve_julia(
        void *(*mallocp) (size_t),
        int32_t *rp_ld,
        int32_t *rp_nr_vars,
        int32_t *rp_dim,
        int32_t *rp_dquot,
        int32_t **rp_lens,
        char ***rp_var_namesp,
        void **rp_cfs_linear_form,
        void **rp_cfs,
        int32_t *n_real_sols,
        void **real_sols_num,
        int32_t **real_sols_den,
        int32_t *lens,
        int32_t *exps,
        void *cfs,
        char **var_names,
        char *output_file,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t initial_hts,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t print_gb,
        const int32_t get_param,
        const int32_t genericity_handling,
        const int32_t precision,
        const int32_t info_level
        );

void free_msolve_julia_result_data(
        void (*freep) (void *),
        int32_t **res_len,
        void **res_cf,
        void **sols_num,
        int32_t **sols_den,
        const int64_t res_ld,
        const int64_t nr_sols,
        const int64_t field_char
        );

int64_t export_groebner_qq(
                void *(*mallocp) (size_t),
        /* return values */
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        /* input values */
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
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t truncate_lifting,
        const int32_t info_level
        );
#endif
