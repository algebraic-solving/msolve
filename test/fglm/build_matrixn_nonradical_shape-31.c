#include "../../src/msolve/libmsolve.c"

int main(void)
{
    int32_t la_option             = 2; // by default
    int32_t use_signatures        = 0;
    int32_t nr_threads            = 1;
    int32_t info_level            = 0;
    int32_t initial_hts           = 17;
    int32_t max_pairs             = 0;
    int32_t elim_block_len        = 0;
    int32_t update_ht             = 0;
    int32_t generate_pbm          = 0;
    int32_t reduce_gb             = 1;
    int32_t print_gb              = 0;
    int32_t genericity_handling   = 2;
    int32_t saturate              = 0;
    int32_t colon                 = 0;
    int32_t normal_form           = 0;
    int32_t normal_form_matrix    = 0;
    int32_t is_gb                 = 0;
    int32_t get_param             = 0;
    int32_t precision             = 128;
    int32_t refine                = 0; /* not used at the moment */
    int32_t isolate               = 0; /* not used at the moment */
    files_gb *files = malloc(sizeof(files_gb));
    files->in_file  = "test/input_files/nonradical_shape-31.ms";
    files->out_file  = NULL;
    FILE *fh  = fopen(files->in_file, "r");

    //  int32_t mon_order   = 0;
    int32_t nr_vars     = 0;
    int32_t field_char  = 9001;
    int32_t nr_gens     = 0;
    data_gens_ff_t *gens = allocate_data_gens();

    get_data_from_file(files->in_file, &nr_vars, &field_char,&nr_gens,gens);

    if (nr_vars != 4) return 1;
    if (field_char != 1073741831) return 1;
    if (nr_gens != 4) return 1;

    gens->rand_linear           = 0;
    gens->random_linear_form = malloc(sizeof(int32_t)*(nr_vars));

    param_t *param  = NULL;
    mpz_param_t mpz_param;
    mpz_param_init(mpz_param);

    long nb_real_roots      = 0;
    interval *real_roots    = NULL;
    real_point_t *real_pts  = NULL;

    param_t ** paramp = &param;
    mpz_param_t *mpz_paramp = &mpz_param;
    long *nb_real_roots_ptr = &nb_real_roots;
    interval **real_roots_ptr = &real_roots;
    real_point_t **real_pts_ptr = &real_pts;

    int32_t *bld    = NULL;
    int32_t **blen  = NULL;
    int32_t **bexp  = NULL;
    void **bcf      = NULL;
    int b           = 0;
    /* counter for randomly chosen linear forms */
    int round = -1;

    bld   = malloc(sizeof(int32_t));
    blen  = malloc(sizeof(int32_t *));
    bexp  = malloc(sizeof(int32_t *));
    bcf   = malloc(sizeof(void *));

    bs_t *bs    = NULL;
    ht_t *bht   = NULL;
    stat_t *st  = NULL;

    int success = 0;

    success = initialize_gba_input_data(&bs, &bht, &st,gens->lens, gens->exps, (void *)gens->cfs,	1073741831, 0 /* DRL order */,elim_block_len, gens->nvars,/* gens->field_char,0 [> DRL order <], gens->nvars, */ gens->ngens, saturate,	initial_hts, nr_threads, max_pairs,	update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,	info_level);
    success = core_gba(&bs, &bht, &st);
    if (!success) {
      printf("Problem with F4, stopped computation.\n");
      return 1;
    }

    export_results_from_gba(bld, blen, bexp,bcf, &malloc, &bs, &bht, &st);

    int32_t *bcf_ff = (int32_t *)(*bcf);
    int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);
    long dquot = 0;
    int32_t *lmb= monomial_basis (bld[0], gens->nvars, bexp_lm,&dquot);
    sp_matfglm_t  *matrix= build_matrixn(lmb, dquot, bld[0], blen, bexp, bcf_ff,bexp_lm, gens->nvars, gens->field_char);

    /* display_fglm_matrix (stdout, matrix); */

    if (matrix->charac != field_char) return 1;
    if (matrix->ncols != 11) return 1;
    if (matrix->nrows != 6) return 1;

    if (matrix->dense_mat[0] != 695220959) return 1;
    if (matrix->dense_mat[1] != 672652169) return 1;
    if (matrix->dense_mat[2] != 1046599061) return 1;
    if (matrix->dense_mat[3] != 260205064) return 1;
    if (matrix->dense_mat[4] != 787559341) return 1;
    if (matrix->dense_mat[5] != 635594395) return 1;
    if (matrix->dense_mat[6] != 82300796) return 1;
    if (matrix->dense_mat[7] != 148803321) return 1;
    if (matrix->dense_mat[8] != 155205489) return 1;
    if (matrix->dense_mat[9] != 638296297) return 1;
    if (matrix->dense_mat[10] != 873613310) return 1;
    if (matrix->dense_mat[11] != 344758724) return 1;
    if (matrix->dense_mat[12] != 836890318) return 1;
    if (matrix->dense_mat[13] != 158248818) return 1;
    if (matrix->dense_mat[14] != 105450484) return 1;
    if (matrix->dense_mat[15] != 676511729) return 1;
    if (matrix->dense_mat[16] != 165475628) return 1;
    if (matrix->dense_mat[17] != 674062047) return 1;
    if (matrix->dense_mat[18] != 386559681) return 1;
    if (matrix->dense_mat[19] != 291442961) return 1;
    if (matrix->dense_mat[20] != 579903534) return 1;
    if (matrix->dense_mat[21] != 141486867) return 1;
    if (matrix->dense_mat[22] != 572023784) return 1;
    if (matrix->dense_mat[23] != 744133717) return 1;
    if (matrix->dense_mat[24] != 348598510) return 1;
    if (matrix->dense_mat[25] != 468551107) return 1;
    if (matrix->dense_mat[26] != 842113423) return 1;
    if (matrix->dense_mat[27] != 659500774) return 1;
    if (matrix->dense_mat[28] != 733019635) return 1;
    if (matrix->dense_mat[29] != 1071125757) return 1;
    if (matrix->dense_mat[30] != 852888838) return 1;
    if (matrix->dense_mat[31] != 981635303) return 1;
    if (matrix->dense_mat[32] != 305324683) return 1;
    if (matrix->dense_mat[33] != 1039333228) return 1;
    if (matrix->dense_mat[34] != 1067654706) return 1;
    if (matrix->dense_mat[35] != 813998442) return 1;
    if (matrix->dense_mat[36] != 829655474) return 1;
    if (matrix->dense_mat[37] != 1054477649) return 1;
    if (matrix->dense_mat[38] != 100193860) return 1;
    if (matrix->dense_mat[39] != 106199953) return 1;
    if (matrix->dense_mat[40] != 1009183349) return 1;
    if (matrix->dense_mat[41] != 74459732) return 1;
    if (matrix->dense_mat[42] != 328387945) return 1;
    if (matrix->dense_mat[43] != 762192905) return 1;
    if (matrix->dense_mat[44] != 137100198) return 1;
    if (matrix->dense_mat[45] != 916451149) return 1;
    if (matrix->dense_mat[46] != 753995230) return 1;
    if (matrix->dense_mat[47] != 639147380) return 1;
    if (matrix->dense_mat[48] != 202835333) return 1;
    if (matrix->dense_mat[49] != 1020314720) return 1;
    if (matrix->dense_mat[50] != 173981585) return 1;
    if (matrix->dense_mat[51] != 37718888) return 1;
    if (matrix->dense_mat[52] != 471998822) return 1;
    if (matrix->dense_mat[53] != 911855634) return 1;
    if (matrix->dense_mat[54] != 877767489) return 1;
    if (matrix->dense_mat[55] != 620591546) return 1;
    if (matrix->dense_mat[56] != 346575422) return 1;
    if (matrix->dense_mat[57] != 5929868) return 1;
    if (matrix->dense_mat[58] != 99564914) return 1;
    if (matrix->dense_mat[59] != 141010962) return 1;
    if (matrix->dense_mat[60] != 1024131791) return 1;
    if (matrix->dense_mat[61] != 824205027) return 1;
    if (matrix->dense_mat[62] != 495749234) return 1;
    if (matrix->dense_mat[63] != 578054943) return 1;
    if (matrix->dense_mat[64] != 143011498) return 1;
    if (matrix->dense_mat[65] != 476740452) return 1;

    if (matrix->triv_idx[0] != 0) return 1;
    if (matrix->triv_idx[1] != 1) return 1;
    if (matrix->triv_idx[2] != 2) return 1;
    if (matrix->triv_idx[3] != 3) return 1;
    if (matrix->triv_idx[4] != 4) return 1;
    if (matrix->triv_pos[0] != 1) return 1;
    if (matrix->triv_pos[1] != 5) return 1;
    if (matrix->triv_pos[2] != 6) return 1;
    if (matrix->triv_pos[3] != 7) return 1;
    if (matrix->triv_pos[4] != 8) return 1;

    if (matrix->dense_idx[0] != 5) return 1;
    if (matrix->dense_idx[1] != 6) return 1;
    if (matrix->dense_idx[2] != 7) return 1;
    if (matrix->dense_idx[3] != 8) return 1;
    if (matrix->dense_idx[4] != 9) return 1;
    if (matrix->dense_idx[5] != 10) return 1;

    if (matrix->dst[0] != 0) return 1;
    if (matrix->dst[1] != 0) return 1;
    if (matrix->dst[2] != 0) return 1;
    if (matrix->dst[3] != 0) return 1;
    if (matrix->dst[4] != 0) return 1;
    if (matrix->dst[5] != 0) return 1;

    return 0;
}