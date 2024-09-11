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
    files->in_file  = "input_files/nonradical_radicalshape-31.ms";
    files->out_file  = NULL;
    FILE *fh  = fopen(files->in_file, "r");

    //  int32_t mon_order   = 0;
    int32_t nr_vars     = 0;
    int32_t field_char  = 9001;
    int32_t nr_gens     = 0;
    data_gens_ff_t *gens = allocate_data_gens();

    get_data_from_file(files->in_file, &nr_vars, &field_char,&nr_gens,gens);

    if (nr_vars != 2) return 101;
    if (field_char != 1073741827) return 102;
    if (nr_gens != 3) return 103;

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

    bs_t *bs  = NULL;
    ht_t *bht = NULL;
    md_t *st  = NULL;

    int error   = 0;
    int success = 0;

    success = initialize_gba_input_data(&bs, &bht, &st,gens->lens, gens->exps, (void *)gens->cfs,	1073741827, 0 /* DRL order */,elim_block_len, gens->nvars,/* gens->field_char,0 [> DRL order <], gens->nvars, */ gens->ngens, saturate,	initial_hts, nr_threads, max_pairs,	update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,	0/*truncate_lifting*/, info_level);
    bs = core_gba(bs, st, &error, 1073741827);
    if (!success || error) {
      printf("Problem with F4, stopped computation.\n");
      return 104;
    }

    export_results_from_gba(bld, blen, bexp,bcf, &malloc, &bs, &bht, &st);

    int32_t *bcf_ff = (int32_t *)(*bcf);
    int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);
    long dquot = 0;
    int32_t *lmb= monomial_basis (bld[0], gens->nvars, bexp_lm,&dquot);
    sp_matfglm_t  *matrix= build_matrixn(lmb, dquot, bld[0], blen, bexp, bcf_ff,bexp_lm, gens->nvars, gens->field_char);

    /* display_fglm_matrix (stdout, matrix); */

    if (matrix->charac != field_char) return 105;
    if (matrix->ncols != 3) return 106;
    if (matrix->nrows != 2) return 107;

    if (matrix->dense_mat[0] != 0) return 201;
    if (matrix->dense_mat[1] != 0) return 202;
    if (matrix->dense_mat[2] != 0) return 203;
    if (matrix->dense_mat[3] != 0) return 204;
    if (matrix->dense_mat[4] != 0) return 205;
    if (matrix->dense_mat[5] != 0) return 206;

    if (matrix->triv_idx[0] != 0) return 207;
    if (matrix->triv_pos[0] != 1) return 208;

    if (matrix->dense_idx[0] != 1) return 209;
    if (matrix->dense_idx[1] != 2) return 210;

    if (matrix->dst[0] != 3) return 211;
    if (matrix->dst[1] != 3) return 212;

    return 0;
}
