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


#ifndef GB_ENGINE_H
#define GB_ENGINE_H

#include "data.h"

int initialize_gba_input_data(
        bs_t **bsp,
        ht_t **bhtp,
        md_t **stp,
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        uint32_t field_char,
        int32_t mon_order,
        int32_t elim_block_len,
        int32_t nr_vars,
        int32_t nr_gens,
        int32_t nr_nf,
        int32_t ht_size,
        int32_t nr_threads,
        int32_t max_nr_pairs,
        int32_t reset_ht,
        int32_t la_option,
        int32_t use_signatures,
        int32_t reduce_gb,
        int32_t pbm_file,
        int32_t truncate_lifting,
        int32_t info_level
        );

bs_t *core_gba(
        bs_t *bs,
        md_t *md,
        int32_t *errp,
        const len_t fc
        );

int64_t export_results_from_gba(
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

bs_t *gba_trace_learning_phase(
        trace_t *trace,           /* trace of the GB Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        ht_t *gbht,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        md_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        );

bs_t *gba_trace_application_phase(
        trace_t *trace,           /* trace of the GB Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        ht_t *lbht,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        md_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        );
#endif
