#include "../../../src/neogb/io.h"

int main(void)
{
    const int32_t lens[]    =   {2, 2};
    const int32_t cfs[]     =   {1, -1, 1, -1};
    uint32_t field_char     =   32003;
    int32_t mon_order       =   0;
    int32_t elim_block_len  =   1;
    int32_t nr_vars         =   2;
    int32_t nr_gens         =   2;
    int32_t nr_nf           =   1; /* must be <= nr_gens */
    int32_t ht_size         =   17;
    int32_t nr_threads      =   1;
    int32_t max_nr_pairs    =   0;
    int32_t reset_ht        =   2;
    int32_t la_option       =   45; /* incorrect: should be fixed to 2 */
    int32_t use_signatures  =   0;
    int32_t reduce_gb       =   -1; /* incorrect: should be fixed to 0 */
    int32_t truncate_lifting =  -1; /* incorrect: should be fixed to 0 */
    int32_t info_level      =   2;
    int *invalid_gens       =   NULL;

    int res = validate_input_data(&invalid_gens, cfs, lens,
            &field_char, &mon_order, &elim_block_len, &nr_vars,
            &nr_gens, &nr_nf, &ht_size, &nr_threads, &max_nr_pairs,
            &reset_ht, &la_option, &use_signatures, &reduce_gb,
            &truncate_lifting, &info_level);

    if (res == -1) return 1;
    if (res == 0) return 1;

    if (invalid_gens[0] != 0 || invalid_gens[1] != 0) return 1; 
    if (field_char != 32003) return 1;
    if (mon_order != 0) return 1;
    if (elim_block_len != 1) return 1;
    if (nr_vars!= 2) return 1;
    if (nr_gens != 2) return 1;
    if (nr_nf != 1) return 1;
    if (ht_size != 17) return 1;
    if (nr_threads != 1) return 1;
    if (max_nr_pairs != 0) return 1;
    if (reset_ht != 2) return 1;
    if (la_option != 2) return 1; /* fixed */
    if (use_signatures != 0) return 1;
    if (reduce_gb != 0) return 1; /* fixed */
    if (truncate_lifting != 0) return 1;
    if (info_level != 2) return 1;

    return 0;
}
