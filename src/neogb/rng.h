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


#ifndef MSOLVE_RNG_H
#define MSOLVE_RNG_H

#include <stdint.h>
#include "mt.h"

/* Platform independent pseudo-random number generation based on the
 * Mersenne Twister MT19937 (M. Matsumoto and T. Nishimura, 1998, with the
 * 2002 seeding procedure). The generator is GSL's gsl_rng_mt19937 taken from
 * GSL 1.9 (see mt.c): for a given seed the generated sequence is
 * identical on all platforms and coincides with the one of GSL (including
 * the convention that seed 0 is replaced by 4357).
 *
 * Each thread owns its own generator (see msolve_rng_current()), the one
 * of the main thread is seeded by msolve_srand(). Code drawing random
 * numbers inside an OpenMP parallel region must not depend on the thread
 * scheduling: it should use a local generator seeded via
 * msolve_rng_derive_seed() from a base seed drawn before entering the
 * parallel region and from the loop index. */

/* largest value returned by msolve_rand() and msolve_rng_rand() */
#define MSOLVE_RAND_MAX 0x7fffffff

/* state of the generator, see mt.h */
typedef mt_state_t msolve_rng_t;

/* seeds rng with seed, seed 0 is replaced by 4357 as in GSL */
void msolve_rng_seed(msolve_rng_t *rng, uint32_t seed);

/* returns the next 32 bit output of rng, i.e. the same value as
 * gsl_rng_get() for gsl_rng_mt19937 */
uint32_t msolve_rng_get(msolve_rng_t *rng);

/* returns the next output of rng in [0, MSOLVE_RAND_MAX] (the upper 31 bits
 * of msolve_rng_get()), as a drop-in replacement for rand() */
int32_t msolve_rng_rand(msolve_rng_t *rng);

/* returns the seed base + idx (mod 2^32) for a child generator, depending
 * only on base and idx */
uint32_t msolve_rng_derive_seed(uint32_t base, uint32_t idx);

/* returns the generator of the calling thread; if it was never seeded it
 * is seeded with 0 (i.e. 4357) */
msolve_rng_t *msolve_rng_current(void);

/* seeds the generator of the calling thread */
void msolve_srand(uint32_t seed);

/* draws from the generator of the calling thread: msolve_rand() returns a
 * value in [0, MSOLVE_RAND_MAX], msolve_rand_u32() a full 32 bit value */
int32_t msolve_rand(void);
uint32_t msolve_rand_u32(void);

#endif
