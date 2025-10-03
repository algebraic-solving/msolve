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


#ifndef GB_SBA_H
#define GB_SBA_H

#include "data.h"

void add_syzygy_schreyer(
        crit_t *syz,
        const hm_t sm,
        const len_t si,
        const ht_t * const ht
        );

int core_sba_schreyer(
        bs_t **bsp,
        ht_t **bhtp,
        md_t **stp
        );
#endif
