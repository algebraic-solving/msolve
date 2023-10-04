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



/* Initialization of sparse fglm matrix for crt */
static inline void crt_mpz_matfglm_initset(crt_mpz_matfglm_t crt_mat,
                                           sp_matfglm_t *mod_mat){
  crt_mat->ncols = mod_mat->ncols;
  crt_mat->nrows = mod_mat->nrows;
  crt_mat->dense_mat = calloc(crt_mat->ncols*crt_mat->nrows,
                              sizeof(mpz_t));
  uint64_t sz = crt_mat->nrows * crt_mat->ncols;
  for(uint64_t i = 0; i < sz; i++){
    mpz_init_set_ui(crt_mat->dense_mat[i], mod_mat->dense_mat[i]);
  }
  long diff = crt_mat->ncols - crt_mat->nrows;
  crt_mat->triv_idx = calloc(diff, sizeof(uint32_t));
  crt_mat->triv_pos = calloc(diff, sizeof(uint32_t));
  crt_mat->dense_idx = calloc(crt_mat->nrows, sizeof(uint32_t));
  crt_mat->dst = calloc(crt_mat->nrows, sizeof(uint32_t));

  for(long i = 0; i < diff; i++){
    crt_mat->triv_idx[i]= mod_mat->triv_idx[i];
    crt_mat->triv_pos[i]= mod_mat->triv_pos[i];
  }
  for(long i = 0; i < crt_mat->nrows; i++){
    crt_mat->dense_idx[i] = mod_mat->dense_idx[i];
    crt_mat->dst[i] = mod_mat->dst[i];
  }
}

/* Initialization of sparse mpz fglm matrix  */
static inline void mpz_matfglm_initset(mpz_matfglm_t mpz_mat,
                                       sp_matfglm_t *mod_mat){
  mpz_mat->ncols = mod_mat->ncols;
  mpz_mat->nrows = mod_mat->nrows;
  mpz_mat->dense_mat = calloc(mpz_mat->ncols*mpz_mat->nrows,
                              sizeof(mpz_t));
  uint64_t sz = mpz_mat->nrows * mpz_mat->ncols;
  for(uint64_t i = 0; i < sz; i++){
    mpz_init_set_ui(mpz_mat->dense_mat[i], 0);
  }
  mpz_mat->denoms = calloc(mpz_mat->nrows,
                           sizeof(mpz_t));
  for(uint64_t i = 0; i < mpz_mat->nrows; i++){
    mpz_init_set_ui(mpz_mat->denoms[i], 1);
  }
  long diff = mpz_mat->ncols - mpz_mat->nrows;
  mpz_mat->triv_idx = calloc(diff, sizeof(uint32_t));
  mpz_mat->triv_pos = calloc(diff, sizeof(uint32_t));
  mpz_mat->dense_idx = calloc(mpz_mat->nrows, sizeof(uint32_t));
  mpz_mat->dst = calloc(mpz_mat->nrows, sizeof(uint32_t));

  for(long i = 0; i < diff; i++){
    mpz_mat->triv_idx[i]= mod_mat->triv_idx[i];
    mpz_mat->triv_pos[i]= mod_mat->triv_pos[i];
  }
  for(long i = 0; i < mpz_mat->nrows; i++){
    mpz_mat->dense_idx[i] = mod_mat->dense_idx[i];
    mpz_mat->dst[i] = mod_mat->dst[i];
  }
}

/* Initialization of sparse fglm matrix for rational reconstruction */
static inline void mpq_matfglm_initset(mpq_matfglm_t mpq_mat,
                                       sp_matfglm_t *mod_mat){
  mpq_mat->ncols = mod_mat->ncols;
  mpq_mat->nrows = mod_mat->nrows;
  mpq_mat->dense_mat = calloc(2*mpq_mat->ncols*mpq_mat->nrows,
                              sizeof(mpz_t));
  uint64_t nc = 2*mpq_mat->ncols;

  for(uint32_t i = 0; i < mpq_mat->nrows; i++){
    uint64_t c = 2*i*mpq_mat->ncols;
    for(uint32_t j = 0; j < nc; j++){
      mpz_init_set_ui(mpq_mat->dense_mat[c+j], 0);
      j++;
      mpz_init_set_ui(mpq_mat->dense_mat[c+j], 1);
    }
  }
  long diff = mpq_mat->ncols - mpq_mat->nrows;
  mpq_mat->triv_idx = calloc(diff, sizeof(uint32_t));
  mpq_mat->triv_pos = calloc(diff, sizeof(uint32_t));
  mpq_mat->dense_idx = calloc(mpq_mat->nrows, sizeof(uint32_t));
  mpq_mat->dst = calloc(mpq_mat->nrows, sizeof(uint32_t));

  for(long i = 0; i < diff; i++){
    mpq_mat->triv_idx[i]= mod_mat->triv_idx[i];
    mpq_mat->triv_pos[i]= mod_mat->triv_pos[i];
  }
  for(long i = 0; i < mpq_mat->nrows; i++){
    mpq_mat->dense_idx[i] = mod_mat->dense_idx[i];
    mpq_mat->dst[i] = mod_mat->dst[i];
  }
}

static inline void trace_det_initset(trace_det_fglm_mat_t trace_det,
                                     uint32_t trace_mod, uint32_t det_mod,
                                     uint32_t tridx, uint32_t detidx){
  mpz_init_set_ui(trace_det->trace_crt, trace_mod);
  mpz_init_set_ui(trace_det->det_crt, det_mod);
  mpz_init_set_ui(trace_det->trace_num, 0);
  mpz_init_set_ui(trace_det->trace_den, 1);
  mpz_init_set_ui(trace_det->det_num, 0);
  mpz_init_set_ui(trace_det->det_den, 1);
  mpz_init(trace_det->tmp);
  trace_det->check_trace = 0;
  trace_det->check_det = 0;
  trace_det->done_trace = 0;
  trace_det->done_det = 0;
  trace_det->trace_idx = tridx;
  trace_det->det_idx = detidx;
}

static inline void trace_det_clear(trace_det_fglm_mat_t trace_det){
  mpz_clear(trace_det->trace_crt);
  mpz_clear(trace_det->det_crt);
  mpz_clear(trace_det->trace_num);
  mpz_clear(trace_det->trace_den);
  mpz_clear(trace_det->det_num);
  mpz_clear(trace_det->det_den);
  mpz_clear(trace_det->tmp);
}

static inline void crt_lift_trace_det(trace_det_fglm_mat_t trace_det,
                                      uint32_t trace_mod, uint32_t det_mod,
                                      mpz_t modulus, mpz_t prod,
                                      uint32_t prime){
  mpz_CRT_ui(trace_det->trace_crt, trace_det->trace_crt,
             modulus, trace_mod, prime, prod, trace_det->tmp, 1);
  mpz_CRT_ui(trace_det->det_crt, trace_det->det_crt,
             modulus, det_mod, prime, prod, trace_det->tmp, 1);
}

static inline void crt_lift_dense_rows(mpz_t *rows, uint32_t *mod_rows,
                                       const uint64_t sz,
                                       mpz_t modulus, 
                                       mpz_t prod,
                                       int32_t prime,
                                       mpz_t tmp,
                                       const int nthrds){
  len_t i;

  for(i = 0; i < sz; i++){

    mpz_CRT_ui(rows[i], rows[i], modulus,
               mod_rows[i], prime, prod, tmp, 1);

  }
}

static inline void crt_lift_mat(crt_mpz_matfglm_t mat, sp_matfglm_t *mod_mat,
                                mpz_t modulus, mpz_t prod_crt,
                                const int32_t prime,
                                mpz_t tmp,
                                const int nthrds){
  /*assumes prod_crt = modulus * prime */
  const uint64_t sz = mat->nrows * mat->ncols;
  crt_lift_dense_rows(mat->dense_mat, mod_mat->dense_mat,
                      sz, modulus, prod_crt, prime, tmp, nthrds);
}

static inline void build_mpz_matrix(mpq_matfglm_t mpq_mat, mpz_matfglm_t mpz_mat){
  fprintf(stderr, "TODO\n");
}

static inline int rat_recon_trace_det(trace_det_fglm_mat_t trace_det,
                                      rrec_data_t recdata, mpz_t modulus,
                                      mpz_t rnum, mpz_t rden){
  int b = ratrecon(rnum, rden, trace_det->trace_crt, modulus, recdata);
  if(b == 1){
    mpz_set(trace_det->trace_num, rnum);
    mpz_set(trace_det->trace_den, rden);
  }
  else
  {
    return 0;
  }
  b = ratrecon(rnum, rden, trace_det->det_crt, modulus, recdata);
  if(b == 1){
    mpz_set(trace_det->det_num, rnum);
    mpz_set(trace_det->det_den, rden);
  }
  else{
    return 0;
  }
  return 1;
}


static inline int check_trace(trace_det_fglm_mat_t trace_det,
                              uint32_t trace_mod,
                              uint32_t prime){

  uint32_t lc = mpz_fdiv_ui(trace_det->trace_den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(trace_det->trace_num, prime);
  c *= lc;
  c = c % prime;

  return (c == trace_mod);
}

static inline int check_det(trace_det_fglm_mat_t trace_det,
                            uint32_t det_mod,
                            uint32_t prime){

  uint32_t lc = mpz_fdiv_ui(trace_det->det_den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(trace_det->det_num, prime);
  c *= lc;
  c = c % prime;

  return (c==det_mod);
}


#define NEW 1

static inline int64_t rat_recon_dense_rows(mpq_matfglm_t mpq_mat,
                                           crt_mpz_matfglm_t crt_mat,
                                           mpz_matfglm_t mpz_mat,
                                           mpz_t modulus, rrec_data_t rdata,
                                           mpz_t rnum, mpz_t rden,
                                           long *matrec){
  const uint32_t nrows = crt_mat->nrows;
  const uint32_t ncols = crt_mat->ncols;
  int64_t cnt = 0;
#ifdef NEW
  mpz_t lcm, coef;
  mpz_init(lcm);
  mpz_init(coef);
#endif
  for(uint32_t i = 0; i < nrows; i++){
    uint64_t c = i*ncols;
#ifdef NEW
    mpz_set_ui(lcm, 1);
#endif
    for(uint32_t j = 0; j < ncols; j++){
      int b = 1;
      if(*matrec <= c+j){
#ifdef NEW
        mpz_mul(coef, crt_mat->dense_mat[c+j], lcm);
        mpz_mod(coef, coef, modulus);
        b = ratrecon(rnum, rden, coef,
                     modulus, rdata);
#else
        b = ratrecon(rnum, rden, crt_mat->dense_mat[c+j],
                     modulus, rdata);
#endif
        if(b == 1){
          mpz_set(mpq_mat->dense_mat[2*c+2*j], rnum);
          mpz_set(mpq_mat->dense_mat[2*c+2*j+1], rden);
          mpz_lcm(lcm, lcm, rden);
          cnt++;
        }
        else{
          if(cnt > *matrec+1){
            fprintf(stderr, "<%.2f%%>", (100*(((double)cnt)/ncols))/nrows );
          }
          *matrec = c;
#ifdef NEW
          mpz_clear(lcm);
          mpz_clear(coef);
#endif
          return c;
        }
      }
      else{
        cnt++;
      }
    }
  }
  fprintf(stderr, "<100.0%%>\n");
  build_mpz_matrix(mpq_mat, mpz_mat);
  *matrec = cnt;
#ifdef NEW
  mpz_clear(coef);
  mpz_clear(lcm);
#endif
  return cnt;
}



void initialize_rrec_data(rrec_data_t recdata){
  mpz_init(recdata->r0);
  mpz_set_ui(recdata->r0, 0);
  mpz_init(recdata->r1);
  mpz_set_ui(recdata->r1, 0);
  mpz_init(recdata->t0);
  mpz_set_ui(recdata->t0, 0);
  mpz_init(recdata->t1);
  mpz_set_ui(recdata->t1, 0);
  mpz_init(recdata->q);
  mpz_set_ui(recdata->q, 0);
  mpz_init(recdata->tmp);
  mpz_set_ui(recdata->tmp, 0);
  mpz_init(recdata->N);
  mpz_set_ui(recdata->N, 0);
  mpz_init(recdata->D);
  mpz_set_ui(recdata->D, 0);
}

void free_rrec_data(rrec_data_t recdata){
  mpz_clear(recdata->r0);
  mpz_clear(recdata->r1);
  mpz_clear(recdata->t0);
  mpz_clear(recdata->t1);
  mpz_clear(recdata->q);
  mpz_clear(recdata->tmp);
  mpz_clear(recdata->N);
  mpz_clear(recdata->D);
}
