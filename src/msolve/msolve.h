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

#include "data.h"

int msolve_ff_alloc(
        param_t **bparam,
        int32_t *bld,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        data_gens_ff_t *gens,
        int32_t initial_hts,
        int32_t nr_threads,
        int32_t max_pairs,
        int32_t update_ht,
        int32_t la_option,
        int32_t info_level,
        int32_t print_gb,
        files_gb *files
        );

int modular_run_msolve(
        param_t **bparam,
        data_gens_ff_t *gens,
        int32_t initial_hts,
        int32_t nr_threads,
        int32_t max_pairs,
        int32_t update_ht,
        int32_t la_option,
        int32_t info_level,
        files_gb *files,
        int32_t prime
        );

int msolve_trace_qq(
        mpz_param_t mpz_param,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t reset_ht,
        int32_t la_option,
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
        int32_t reset_ht,
        int32_t la_option,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        files_gb *files,
        int round
        );

int msolve_qq(
        mpz_param_t mp_param,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t reset_ht,
        int32_t la_option,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        files_gb *files,
        int
        );

int real_msolve_qq(
        mpz_param_t mp_param,
        param_t **nmod_param,
        int *dim_ptr,
        long *dquot_ptr,
        long *nb_real_roots_ptr,
        interval **real_roots_ptr,
        real_point_t **real_pts_ptr,
        data_gens_ff_t *gens,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t reset_ht,
        int32_t la_option,
        int32_t info_level,
        int32_t print_gb,
        int32_t pbm_file,
        int32_t precision,
        files_gb *files,
        int
        );

int core_msolve(
        int32_t la_option,
        int32_t nr_threads,
        int32_t info_level,
        int32_t initial_hts,
        int32_t max_pairs,
        int32_t update_ht,
        int32_t generate_pbm,
        int32_t reduce_gb,
        int32_t print_gb,
        int32_t get_param,
        int32_t genericity_handling,
        int32_t normal_form,
        int32_t normal_form_matrix,
        int32_t is_gb,
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
        int32_t *rp_ld,
        int32_t *rp_dim,
        int32_t *rp_dquot,
        int32_t **rp_lens,
        void **rp_cfs,
        int32_t *lens,
        int32_t *exps,
        void *cfs,
        char **var_names,
        char *output_file,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t initial_hts,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t print_gb,
        int32_t get_param,
        const int32_t genericity_handling,
        const int32_t precision,
        const int32_t info_level
        );
#endif
