#include "../../../src/neogb/f4.h"
#include "../../../src/neogb/io.h"

/* Computes the reduced Groebner basis of katsura-6 over the prime field
 * of characteristic 2^52-47, the largest supported prime. The result was
 * checked against the reduced Groebner basis over the rationals mapped to
 * this prime field, here only a hash of it is compared. */

static uint64_t hash(uint64_t h, const uint64_t x)
{
    for (int i = 0; i < 8; ++i) {
        h ^= (x >> (8*i)) & 0xff;
        h *= 1099511628211ULL;
    }
    return h;
}

int main(void)
{
    const uint64_t field_char = 4503599627370449ULL;
    const int32_t nr_vars     = 6;
    const int32_t nr_gens     = 6;
    const int32_t lens[] = {7, 7, 6, 6, 5, 5};
    const int64_t cfs[] = {
        1, 2, 2, 2, 2, 2, -1, 1, 2, 2,
        2, 2, 2, -1, 2, 2, 2, 2, 2, -1,
        1, 2, 2, 2, 2, -1, 2, 2, 2, 2,
        -1, 1, 2, 2, 2, -1
    };
    const int32_t exps[] = {
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 0, 0,
        2, 0, 0, 0, 0, 0,
        0, 2, 0, 0, 0, 0,
        0, 0, 2, 0, 0, 0,
        0, 0, 0, 2, 0, 0,
        0, 0, 0, 0, 2, 0,
        0, 0, 0, 0, 0, 2,
        1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 0,
        0, 1, 1, 0, 0, 0,
        0, 0, 1, 1, 0, 0,
        0, 0, 0, 1, 1, 0,
        0, 0, 0, 0, 1, 1,
        0, 1, 0, 0, 0, 0,
        0, 2, 0, 0, 0, 0,
        1, 0, 1, 0, 0, 0,
        0, 1, 0, 1, 0, 0,
        0, 0, 1, 0, 1, 0,
        0, 0, 0, 1, 0, 1,
        0, 0, 1, 0, 0, 0,
        0, 1, 1, 0, 0, 0,
        1, 0, 0, 1, 0, 0,
        0, 1, 0, 0, 1, 0,
        0, 0, 1, 0, 0, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 2, 0, 0, 0,
        0, 1, 0, 1, 0, 0,
        1, 0, 0, 0, 1, 0,
        0, 1, 0, 0, 0, 1,
        0, 0, 0, 0, 1, 0
    };
    const int32_t la_options[]  = {1, 2, 42, 44};
    const int32_t nr_threads[]  = {1, 4};

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            int32_t bld = 0, *blen = NULL, *bexp = NULL;
            void *bcf = NULL;
            const int64_t nterms = export_f4(malloc, &bld, &blen, &bexp,
                    &bcf, lens, exps, cfs, field_char, 0, 0, nr_vars,
                    nr_gens, 17, nr_threads[j], 0, 0, la_options[i], 1, 0, 0);
            const int64_t *cf = (int64_t *)bcf;
            uint64_t h = 14695981039346656037ULL;
            h = hash(h, (uint64_t)bld);
            for (int32_t k = 0; k < bld; ++k) {
                h = hash(h, (uint64_t)blen[k]);
            }
            for (int64_t k = 0; k < nterms * nr_vars; ++k) {
                h = hash(h, (uint64_t)bexp[k]);
            }
            for (int64_t k = 0; k < nterms; ++k) {
                h = hash(h, (uint64_t)cf[k]);
            }
            free(blen);
            free(bexp);
            free(bcf);
            if (bld != 22 || nterms != 528 || h != 0x245c51d9cc3e1546ULL) {
                return 1;
            }
        }
    }

    /* primes >= 2^52 are not supported */
    uint64_t fc = (uint64_t)1 << 52;
    int32_t mon_order = 0, elim_block_len = 0, nv = nr_vars, ng = nr_gens;
    int32_t nr_nf = 0, ht_size = 17, nthrds = 1, max_nr_pairs = 0;
    int32_t reset_ht = 0, la_option = 2, use_signatures = 0, reduce_gb = 1;
    int32_t truncate_lifting = 0, info_level = 0;
    int *invalid_gens = NULL;
    if (validate_input_data_64(&invalid_gens, cfs, lens, &fc, &mon_order,
                &elim_block_len, &nv, &ng, &nr_nf, &ht_size, &nthrds,
                &max_nr_pairs, &reset_ht, &la_option, &use_signatures,
                &reduce_gb, &truncate_lifting, &info_level) != 0) {
        return 1;
    }

    return 0;
}
