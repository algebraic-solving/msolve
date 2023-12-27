
static inline void duplicate_linear_data(int nthreads, int nvars, int nlins,
                                         nvars_t **blinvars, uint32_t **blineqs,
                                         nvars_t **bsquvars){
  for(int i = 1; i < nthreads; i++){

    blineqs[i] = calloc(nlins*(nvars + 1), sizeof(uint64_t));
    for(long j = 0; j < nlins*(nvars + 1); j++){
      blineqs[i][j] = 0;
    }

    blinvars[i] = calloc(nvars, sizeof(uint64_t));
    for(long j = 0; j < nvars; j++){
      blinvars[i][j] = blinvars[0][j];
    }

    bsquvars[i] = calloc(nvars - 1, sizeof(uint64_t));
    for(nvars_t j = 0; j < nvars - 1; j++){
      bsquvars[i][j] = bsquvars[0][j];
    }

  }
}

static inline void duplicate_tracer(
        const int nthreads,
        const bs_t * const bs,
        const md_t * const md,
        trace_t **btrace)
{
    if (btrace[0] != NULL) {
        for(int i = 1; i < nthreads; i++){
            btrace[i]  = initialize_trace(bs, md);

            /* size for trace data */
            btrace[i]->std  = btrace[0]->std;
            /* load of trace data */
            btrace[i]->ltd  = btrace[0]->ltd;
            /* number of lead monomials of non-redundant elements in Gbasis */
            btrace[i]->lml  = btrace[0]->lml;

            /* position of non redundant elements in Gbasis */
            btrace[i]->lmps = (bl_t *)calloc((unsigned long)btrace[0]->lml,
                    sizeof(bl_t));
            memcpy(btrace[i]->lmps, btrace[0]->lmps,
                    (unsigned long)btrace[0]->lml * sizeof(bl_t));

            /* non-redundant lead mon. as short divmask */
            btrace[i]->lm   = (sdm_t *)calloc((unsigned long)btrace[0]->lml,
                    sizeof(sdm_t));
            memcpy(btrace[0]->lm, btrace[i]->lm,
                    (unsigned long)btrace[0]->lml * sizeof(sdm_t));

            /* array of trace data per F4 round */
            btrace[i]->td  = calloc((unsigned long)btrace[0]->ltd, sizeof(td_t));

            for(len_t l = 0; l < btrace[0]->ltd; ++l){

                btrace[i]->td[l].rld = btrace[0]->td[l].rld;
                btrace[i]->td[l].tld = btrace[0]->td[l].tld;
                btrace[i]->td[l].nlm = btrace[0]->td[l].nlm;

                btrace[i]->td[l].rri = calloc(btrace[0]->td[l].rld,
                        sizeof(len_t));
                for(len_t k = 0; k < (btrace[0]->td[l].rld ); ++k){
                    btrace[i]->td[l].rri[k] = btrace[0]->td[l].rri[k] ;
                }

                btrace[i]->td[l].tri = calloc(btrace[0]->td[l].tld,
                        sizeof(len_t));
                for(len_t k = 0; k < (btrace[0]->td[l].tld); ++k){
                    btrace[i]->td[l].tri[k] = btrace[0]->td[l].tri[k] ;
                }

                btrace[i]->td[l].nlms = calloc(btrace[0]->td[l].nlm,
                        sizeof(len_t));
                for(len_t k = 0; k < (btrace[0]->td[l].nlm); ++k){
                    btrace[i]->td[l].nlms[k] = btrace[0]->td[l].nlms[k];
                }

                btrace[i]->td[l].rba = (rba_t **)malloc((unsigned long)btrace[0]->td[l].tld * sizeof(rba_t *));

                const len_t nrr   = btrace[0]->td[l].rld;

                const unsigned long nlrba = nrr / 2 / 32 + (((nrr / 2) % 32) != 0);

                for(len_t j = 0; j < (btrace[0]->td[l].tld / 2); ++j) {
                    btrace[i]->td[l].rba[j]  = calloc(nlrba, sizeof(rba_t));
                }

                for (len_t j = 0; j < (btrace[0]->td[l].tld / 2); ++j) {
                    for(len_t k = 0; k < nlrba; k++){
                        btrace[i]->td[l].rba[j][k]  = btrace[0]->td[l].rba[j][k];
                    }
                }

            }
        }
    }
}

static inline void duplicate_data_mthread_trace(int nthreads,
                                                bs_t *bs,
                                                md_t *st,
                                                int32_t *num_gb,
                                                int32_t **leadmons_ori,
                                                int32_t **leadmons_current,
                                                /* trace_t **btrace, */
                                                fglm_bms_data_t **bdata_bms,
                                                fglm_data_t **bdata_fglm,
                                                int32_t **bstart_cf_gb_xn,
                                                int32_t **blen_gb_xn,
                                                int32_t **bdiv_xn,
                                                sp_matfglm_t **bmatrix,
                                                param_t **nmod_params,
                                                nvars_t nlins,
                                                nvars_t *bnlins,
                                                nvars_t **blinvars,
                                                uint32_t **blineqs,
                                                nvars_t **bsquvars){
  const long len_xn = bmatrix[0]->nrows;
  const long dquot = bmatrix[0]->ncols;
  const long len = num_gb[0] * (st->nvars);

  for(int i = 0; i <nthreads; i++){
    bnlins[i] = nlins;
  }

  for(int i = 0; i < nthreads; i++){
    leadmons_current[i] = (int32_t *)malloc(len*sizeof(int32_t));
  }
  /* leadmons_ori[0] has already been allocated*/
  for(int i = 1; i < nthreads; i++){
    leadmons_ori[i] = (int32_t *)calloc(len, sizeof(int32_t));
    for(long j = 0; j < len; j++){
      leadmons_ori[i][j] = leadmons_ori[0][j];
    }
  }
  for(long i = 1; i < nthreads; i++){
    bstart_cf_gb_xn[i] = malloc(sizeof(int32_t) * len_xn);
    blen_gb_xn[i] = malloc(sizeof(int32_t) * len_xn);
    bdiv_xn[i] = malloc(sizeof(int32_t) * num_gb[0]/* bs->lml */);
    for(long j = 0; j < len_xn; j++){
      bstart_cf_gb_xn[i][j] = bstart_cf_gb_xn[0][j];
      blen_gb_xn[i][j] = blen_gb_xn[0][j];
    }
    for(long j = 0; j < num_gb[0]; j++){
      bdiv_xn[i][j] = bdiv_xn[0][j];
    }
  }
  for(int i=1; i < nthreads; i++){
    num_gb[i] = num_gb[0];
  }
  /* need to duplicate bmatrix */
  for(int i = 1; i < nthreads; i++){

    bmatrix[i] = calloc(1, sizeof(sp_matfglm_t));
    bmatrix[i]->ncols = dquot;
    bmatrix[i]->nrows = len_xn;
    long len1 = dquot * len_xn;
    long len2 = dquot - len_xn;

    sp_matfglm_t *matrix = bmatrix[i];
    if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
      fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len1; j++){
        matrix->dense_mat[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*(dquot - len_xn))){
      fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < (dquot-len_xn); j++){
        matrix->triv_idx[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
      fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len2; j++){
        matrix->triv_pos[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
      fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len_xn; j++){
        matrix->dense_idx[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
      fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len_xn; j++){
        matrix->dst[j] = 0;
      }
    }
  }

  for(int i = 1; i < nthreads; i++){
    bdata_fglm[i] = allocate_fglm_data(len_xn, dquot, (st->nvars));
    bdata_bms[i] = allocate_fglm_bms_data(dquot, 65521);
    nmod_params[i] = allocate_fglm_param(nmod_params[0]->charac, (st->nvars));
    nmod_poly_set(nmod_params[i]->elim, nmod_params[0]->elim);
    nmod_poly_set(nmod_params[i]->denom, nmod_params[0]->denom);
    for(long j = 0; j < (st->nvars) - 2; j++){
      nmod_poly_set(nmod_params[i]->coords[j], nmod_params[0]->coords[j]);
    }
  }

  /* duplicate_tracer(nthreads, bs, st, btrace); */

  duplicate_linear_data(nthreads, st->nvars, nlins,
                        blinvars, blineqs,
                        bsquvars);

}



static inline void duplicate_data_mthread(int nthreads,
                                          nvars_t nv,
                                        int32_t *num_gb,
                                        int32_t **leadmons_ori,
                                        int32_t **leadmons_current,
                                        fglm_bms_data_t **bdata_bms,
                                        fglm_data_t **bdata_fglm,
                                        int32_t **bstart_cf_gb_xn,
                                        int32_t **blen_gb_xn,
                                        int32_t **bdiv_xn,
                                          sp_matfglm_t **bmatrix,
                                          param_t **nmod_params){
  const long len_xn = bmatrix[0]->nrows;
  const long dquot = bmatrix[0]->ncols;
  const long len = num_gb[0] * nv;

  for(int i = 0; i < nthreads; i++){
    leadmons_current[i] = (int32_t *)calloc(len, sizeof(int32_t));
  }
  /* leadmons_ori[0] has already been allocated*/
  for(int i = 1; i < nthreads; i++){
    leadmons_ori[i] = (int32_t *)calloc(len, sizeof(int32_t));
    for(long j = 0; j < len; j++){
      leadmons_ori[i][j] = leadmons_ori[0][j];
    }
  }
  for(long i = 1; i < nthreads; i++){
    bstart_cf_gb_xn[i] = malloc(sizeof(int32_t) * len_xn);
    blen_gb_xn[i] = malloc(sizeof(int32_t) * len_xn);
    bdiv_xn[i] = malloc(sizeof(int32_t) * num_gb[0]/* bs->lml */);
    for(long j = 0; j < len_xn; j++){
      bstart_cf_gb_xn[i][j] = bstart_cf_gb_xn[0][j];
      blen_gb_xn[i][j] = blen_gb_xn[0][j];
    }
    for(long j = 0; j < num_gb[0]; j++){
      bdiv_xn[i][j] = bdiv_xn[0][j];
    }
  }
  for(int i=1; i < nthreads; i++){
    num_gb[i] = num_gb[0];
  }
  /* need to duplicate bmatrix */
  for(int i = 1; i < nthreads; i++){

    bmatrix[i] = calloc(1, sizeof(sp_matfglm_t));
    bmatrix[i]->ncols = dquot;
    bmatrix[i]->nrows = len_xn;
    long len1 = dquot * len_xn;
    long len2 = dquot - len_xn;

    sp_matfglm_t *matrix = bmatrix[i];
    if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
      fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len1; j++){
        matrix->dense_mat[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*(dquot - len_xn))){
      fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < (dquot-len_xn); j++){
        matrix->triv_idx[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
      fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len2; j++){
        matrix->triv_pos[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
      fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len_xn; j++){
        matrix->dense_idx[j] = 0;
      }
    }
    if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
      fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
      exit(1);
    }
    else{
      for(long j = 0; j < len_xn; j++){
        matrix->dst[j] = 0;
      }
    }
  }

  for(int i = 1; i < nthreads; i++){
    bdata_fglm[i] = allocate_fglm_data(len_xn, dquot, nv);
    bdata_bms[i] = allocate_fglm_bms_data(dquot, 65521);
    nmod_params[i] = allocate_fglm_param(nmod_params[0]->charac, nv);
    nmod_poly_set(nmod_params[i]->elim, nmod_params[0]->elim);
    nmod_poly_set(nmod_params[i]->denom, nmod_params[0]->denom);
    for(long j = 0; j < nv - 2; j++){
      nmod_poly_set(nmod_params[i]->coords[j], nmod_params[0]->coords[j]);
    }
  }

  if(nthreads > 1){
    fprintf(stderr, "Duplication of data to be implemented\n");
    exit(1);
  }

}

static inline void duplicate_data_mthread_gbtrace(int nthreads,
                                                  bs_t *bs,
                                                  md_t *st,
                                                  int32_t *num_gb,
                                                  int32_t **leadmons_ori,
                                                  int32_t **leadmons_current,
                                                  trace_t **btrace){


  const len_t len = num_gb[0] * (st->nvars);

  for(int i = 0; i < nthreads; i++){
    leadmons_current[i] = (int32_t *)calloc(len, sizeof(int32_t));
  }
  /* leadmons_ori[0] has already been allocated*/
  for(int i = 1; i < nthreads; i++){
    leadmons_ori[i] = (int32_t *)calloc(len, sizeof(int32_t));
    for(len_t j = 0; j < len; j++){
      leadmons_ori[i][j] = leadmons_ori[0][j];
    }
  }
  for(int i=1; i < nthreads; i++){
    num_gb[i] = num_gb[0];
  }


  duplicate_tracer(nthreads, bs, st, btrace);

}
