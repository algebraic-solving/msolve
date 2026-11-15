#include <stdio.h>
#include "../../../src/neogb/rng.h"

/* returns the n-th output of the generator seeded with seed */
static uint32_t nth_output(uint32_t seed, uint32_t n)
{
    msolve_rng_t rng;
    uint32_t i, v = 0;
    msolve_rng_seed(&rng, seed);
    for (i = 0; i < n; ++i) {
        v = msolve_rng_get(&rng);
    }
    return v;
}

int main(void)
{
    int ret = 0;
    /* reference value of GSL's test suite for gsl_rng_mt19937 */
    if (nth_output(4357, 1000) != 1186927261U) {
        fprintf(stderr, "seed 4357: wrong 1000th output\n");
        ret = 1;
    }
    /* seed 0 is replaced by 4357, as in GSL */
    if (nth_output(0, 1000) != 1186927261U) {
        fprintf(stderr, "seed 0: wrong 1000th output\n");
        ret = 1;
    }
    /* reference value required by the C++ standard for std::mt19937 */
    if (nth_output(5489, 10000) != 4123659995U) {
        fprintf(stderr, "seed 5489: wrong 10000th output\n");
        ret = 1;
    }
    /* the thread generator must follow msolve_srand */
    msolve_srand(5489);
    if (msolve_rand_u32() != nth_output(5489, 1)
            || msolve_rand() != (int32_t)(nth_output(5489, 2) >> 1)) {
        fprintf(stderr, "msolve_srand/msolve_rand mismatch\n");
        ret = 1;
    }
    return ret;
}
