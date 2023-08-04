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


#ifndef GB_TRACE_H
#define GB_TRACE_H

#include "data.h"

trace_t *initialize_trace(
        const bs_t * const bs,
        const md_t * const md
        );


void free_trace(
        trace_t **trp
        );

void free_lucky_primes(
        primes_t **lpp
        );

static inline int is_lucky_prime_ui(
                          const uint32_t prime,
                          const bs_t * const bs
                          )
{
  len_t i, j;
  mpz_t *cf;

  const len_t bl  = bs->ld;

  for (i = 0; i < bl; ++i) {
    cf  = bs->cf_qq[bs->hm[i][COEFFS]];
    for (j = 0; j < bs->hm[i][LENGTH]; ++j) {
      if (mpz_divisible_ui_p(cf[j], prime) != 0) {
        return 1;
      }
    }
  }
  return 0;
}


static inline int is_lucky_prime(
        const mpz_t prime,
        const bs_t * const bs
        )
{
    len_t i, j;
    mpz_t *cf;

    const len_t bl  = bs->ld;

    for (i = 0; i < bl; ++i) {
        cf  = bs->cf_qq[bs->hm[i][COEFFS]];
        for (j = 0; j < bs->hm[i][LENGTH]; ++j) {
            if (mpz_divisible_p(cf[j], prime) != 0) {
                return 1;
            }
        }
    }
    return 0;
}

void reduce_basis_no_hash_table_switching(
        bs_t *bs,
        mat_t *mat,
        ht_t *bht,
        ht_t *sht,
        md_t *st
        );

static inline void generate_lucky_primes(
        primes_t *lp,
        const bs_t * const bs,
        const uint32_t start,
        const len_t nr_new_primes
        )
{
    len_t i;

    lp->old =   lp->ld;
    lp->ld  +=  nr_new_primes;
    lp->p   =   realloc(lp->p, (unsigned long)(lp->ld) * sizeof(uint32_t));

    mpz_t last_prime;
    mpz_init(last_prime);
    if (lp->old == 0) {
        mpz_set_ui(last_prime, start);
    } else {
        mpz_set_ui(last_prime, lp->p[lp->old-1]);
    }
    mpz_nextprime(last_prime, last_prime);

    i = lp->old;
    while (i < lp->ld) {
        if (is_lucky_prime(last_prime, bs) == 0) {
            lp->p[i++] = (int32_t)mpz_get_ui(last_prime);
        }
        mpz_nextprime(last_prime, last_prime);
    }
    mpz_clear(last_prime);
}

bs_t *f4_trace_learning_phase(
        trace_t *trace,           /* trace of the F4 Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        ht_t *gbht,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        md_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        );

bs_t *f4sat_trace_learning_phase_1(
        trace_t *trace,           /* trace of the F4 Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        const bs_t * const gsat,  /* global saturation elements */
        ht_t **gbhtp,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        md_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        );

bs_t *f4sat_trace_learning_phase_2(
        trace_t *trace,           /* trace of the F4 Algorithm */
        ht_t * tht,               /* trace hash table for multipliers */
        const bs_t * const ggb,   /* global basis */
        const bs_t * const gsat,  /* global saturation elements */
        ht_t **gbhtp,               /* global basis hash table, generated
                                   * in this run, used in upcoming runs */
        md_t *gst,              /* global statistics */
        const int32_t fc          /* characteristic of field */
        );

bs_t *f4_trace_application_phase(
        const trace_t * const trace,  /* trace of the F4 Algorithm */
        const ht_t * const tht,       /* trace hash table for multipliers */
        const bs_t * const ggb,       /* global basis */
        ht_t *lbht,                   /* local basis hash table, not shared */
        md_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        );

bs_t *f4sat_trace_application_phase(
        const trace_t * const trace,  /* trace of the F4 Algorithm */
        const ht_t * const tht,       /* trace hash table for multipliers */
        const bs_t * const ggb,       /* global basis */
        const bs_t * const gsat,      /* global saturation elements */
        ht_t *lbht,                   /* local basis hash table, not shared */
        md_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        );

int64_t f4_trace_julia(
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
        const uint32_t prime_start,
        const int32_t nr_primes,
        const int32_t pbm_file,
        const int32_t info_level
        );

bs_t *modular_f4(
        const bs_t * const ggb,       /* global basis */
        ht_t * gbht,                  /* global basis hash table, shared */
        md_t *gst,                  /* global statistics */
        const uint32_t fc             /* characteristic of field */
        );
#endif
