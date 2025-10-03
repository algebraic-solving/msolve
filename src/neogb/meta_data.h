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

#ifndef GB_META_DATA_H
#define GB_META_DATA_H


#include "data.h"

md_t *copy_meta_data(
		     const md_t * const gmd,
		     const int32_t prime
		     );

md_t *allocate_meta_data(
                              void
    );

void print_initial_statistics(
                              FILE *,
                              const md_t *md
    );

void print_tracer_statistics(
                              FILE *,
                              const double rt,
                              const md_t *md
    );

void get_and_print_final_statistics(
                            FILE *file, 
                            md_t *md,
                            const bs_t * const bs
        );
#endif
