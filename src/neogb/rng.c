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


#include "rng.h"

#include "mt.h"

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define MSOLVE_THREAD_LOCAL _Thread_local
#else
#define MSOLVE_THREAD_LOCAL __thread
#endif

static MSOLVE_THREAD_LOCAL msolve_rng_t msolve_rng_thread;
static MSOLVE_THREAD_LOCAL int msolve_rng_thread_seeded = 0;

void msolve_rng_seed(msolve_rng_t *rng, uint32_t seed)
{
    mt_set(rng, seed);
}

uint32_t msolve_rng_get(msolve_rng_t *rng)
{
    return mt_get(rng);
}

int32_t msolve_rng_rand(msolve_rng_t *rng)
{
    return (int32_t)(msolve_rng_get(rng) >> 1);
}

uint32_t msolve_rng_derive_seed(uint32_t base, uint32_t idx)
{
    return base + idx;
}

msolve_rng_t *msolve_rng_current(void)
{
    if (!msolve_rng_thread_seeded) {
        msolve_rng_seed(&msolve_rng_thread, 0);
        msolve_rng_thread_seeded = 1;
    }
    return &msolve_rng_thread;
}

void msolve_srand(uint32_t seed)
{
    msolve_rng_seed(&msolve_rng_thread, seed);
    msolve_rng_thread_seeded = 1;
}

int32_t msolve_rand(void)
{
    return msolve_rng_rand(msolve_rng_current());
}

uint32_t msolve_rand_u32(void)
{
    return msolve_rng_get(msolve_rng_current());
}
