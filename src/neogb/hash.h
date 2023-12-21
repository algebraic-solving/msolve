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

#ifndef GB_HASH_H
#define GB_HASH_H

#include "data.h"

/* we have three different hash tables:
 * 1. one hash table for elements in the basis (bht)
 * 2. one hash table for the spairs during the update process (uht)
 * 3. one hash table for the multiplied elements during symbolic
 *    preprocessing (sht) */

/* The idea of the structure of the hash table is taken from an
 * implementation by Roman Pearce and Michael Monagan in Maple. */

void reset_hash_table_indices(
        ht_t *ht,
        const hi_t * const hcm,
        const len_t len
        );

ht_t *initialize_basis_hash_table(
    md_t *st
    );

ht_t *copy_hash_table(
    const ht_t *bht
    );

ht_t *initialize_secondary_hash_table(
    const ht_t * const ht,
    const md_t * const md
    );

void free_shared_hash_data(
    ht_t *ht
    );

void free_hash_table(
    ht_t **htp
    );

void full_free_hash_table(
                     ht_t **htp
                     );

void calculate_divmask(
    ht_t *ht
    );
#endif
