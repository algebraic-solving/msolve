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

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
/* Initialization of sparse fglm matrix for crt */
#include "msolve-data.h"
static inline void crt_mpz_matfglm_initset(crt_mpz_matfglm_t crt_mat,
                                           sp_matfglm_t *mod_mat){
  crt_mat->ncols = mod_mat->ncols;
  crt_mat->nrows = mod_mat->nrows;
  uint64_t sz = crt_mat->nrows * crt_mat->ncols;
  crt_mat->dense_mat = (mpz_t *)malloc(sz * sizeof(mpz_t));
  for (uint64_t i = 0; i < sz; i++){
    mpz_init_set_ui(crt_mat->dense_mat[i], mod_mat->dense_mat[i]);
  }
  long diff = crt_mat->ncols - crt_mat->nrows;
  crt_mat->triv_idx = malloc(diff * sizeof(uint32_t));
  crt_mat->triv_pos = malloc(diff * sizeof(uint32_t));
  crt_mat->dense_idx = malloc(crt_mat->nrows * sizeof(uint32_t));
  crt_mat->dst = malloc(crt_mat->nrows * sizeof(uint32_t));

  for (long i = 0; i < diff; i++){
    crt_mat->triv_idx[i] = mod_mat->triv_idx[i];
    crt_mat->triv_pos[i] = mod_mat->triv_pos[i];
  }
  for (long i = 0; i < crt_mat->nrows; i++){
    crt_mat->dense_idx[i] = mod_mat->dense_idx[i];
    crt_mat->dst[i] = mod_mat->dst[i];
  }
}

static inline void crt_mpz_matfglm_clear(crt_mpz_matfglm_t crt_mat){
  uint64_t sz = crt_mat->nrows * crt_mat->ncols;

  for (uint64_t i = 0; i < sz; i++){
    mpz_clear(crt_mat->dense_mat[i]);
  }
  free(crt_mat->dense_mat);

  free(crt_mat->triv_idx);
  free(crt_mat->triv_pos);
  free(crt_mat->dense_idx);
  free(crt_mat->dst);
}

/* Initialization of sparse mpz fglm matrix  */
static inline void mpz_matfglm_initset(mpz_matfglm_t mpz_mat,
                                       sp_matfglm_t *mod_mat) {
  mpz_mat->ncols = mod_mat->ncols;
  mpz_mat->nrows = mod_mat->nrows;
  mpz_mat->dense_mat = malloc(mpz_mat->ncols * mpz_mat->nrows * sizeof(mpz_t));
  uint64_t sz = mpz_mat->nrows * mpz_mat->ncols;
  for (uint64_t i = 0; i < sz; i++) {
    mpz_init(mpz_mat->dense_mat[i]);
  }
  mpz_mat->denoms = malloc(mpz_mat->nrows * sizeof(mpz_t));
  for (uint64_t i = 0; i < mpz_mat->nrows; i++) {
    mpz_init(mpz_mat->denoms[i]);
  }
  long diff = mpz_mat->ncols - mpz_mat->nrows;
  mpz_mat->triv_idx = calloc(diff, sizeof(uint32_t));
  mpz_mat->triv_pos = calloc(diff, sizeof(uint32_t));
  mpz_mat->dense_idx = calloc(mpz_mat->nrows, sizeof(uint32_t));
  mpz_mat->dst = calloc(mpz_mat->nrows, sizeof(uint32_t));

  for (long i = 0; i < diff; i++) {
    mpz_mat->triv_idx[i] = mod_mat->triv_idx[i];
    mpz_mat->triv_pos[i] = mod_mat->triv_pos[i];
  }
  for (long i = 0; i < mpz_mat->nrows; i++) {
    mpz_mat->dense_idx[i] = mod_mat->dense_idx[i];
    mpz_mat->dst[i] = mod_mat->dst[i];
  }
}

static inline void mpz_matfglm_clear(mpz_matfglm_t mpz_mat) {

  uint64_t sz = mpz_mat->nrows * mpz_mat->ncols;
  for (uint64_t i = 0; i < sz; i++) {
    mpz_clear(mpz_mat->dense_mat[i]);
  }
  free(mpz_mat->dense_mat);
  for (uint64_t i = 0; i < mpz_mat->nrows; i++) {
    mpz_clear(mpz_mat->denoms[i]);
  }
  free(mpz_mat->denoms);

  free(mpz_mat->triv_idx);
  free(mpz_mat->triv_pos);
  free(mpz_mat->dense_idx);
  free(mpz_mat->dst);
}

/* Initialization of sparse fglm matrix for rational reconstruction */
static inline void mpq_matfglm_initset(mpq_matfglm_t mpq_mat,
                                       sp_matfglm_t *mod_mat) {
  mpq_mat->ncols = mod_mat->ncols;
  mpq_mat->nrows = mod_mat->nrows;
  mpq_mat->dense_mat =
      malloc(2 * mpq_mat->ncols * mpq_mat->nrows * sizeof(mpz_t));
  mpq_mat->denoms = malloc(mpq_mat->nrows * sizeof(mpz_t));
  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    mpz_init(mpq_mat->denoms[i]);
  }
  uint64_t nc = 2 * mpq_mat->ncols;

  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    uint64_t c = 2 * i * mpq_mat->ncols;
    for (uint32_t j = 0; j < nc; j++) {
      mpz_init(mpq_mat->dense_mat[c + j]);
      j++;
      mpz_init(mpq_mat->dense_mat[c + j]);
    }
  }
  long diff = mpq_mat->ncols - mpq_mat->nrows;
  mpq_mat->triv_idx = calloc(diff, sizeof(uint32_t));
  mpq_mat->triv_pos = calloc(diff, sizeof(uint32_t));
  mpq_mat->dense_idx = calloc(mpq_mat->nrows, sizeof(uint32_t));
  mpq_mat->dst = calloc(mpq_mat->nrows, sizeof(uint32_t));

  for (long i = 0; i < diff; i++) {
    mpq_mat->triv_idx[i] = mod_mat->triv_idx[i];
    mpq_mat->triv_pos[i] = mod_mat->triv_pos[i];
  }
  for (long i = 0; i < mpq_mat->nrows; i++) {
    mpq_mat->dense_idx[i] = mod_mat->dense_idx[i];
    mpq_mat->dst[i] = mod_mat->dst[i];
  }
}

static inline void mpq_matfglm_clear(mpq_matfglm_t mpq_mat) {
  /* entries of mpq_mat->denoms are cleared previously when building mpz_mat */
  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    mpz_clear(mpq_mat->denoms[i]);
  }
  free(mpq_mat->denoms);

  uint64_t nc = 2 * mpq_mat->ncols;
  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    uint64_t c = 2 * i * mpq_mat->ncols;
    for (uint32_t j = 0; j < nc; j++) {
      mpz_clear(mpq_mat->dense_mat[c + j]);
    }
  }
  free(mpq_mat->dense_mat);

  free(mpq_mat->triv_idx);
  free(mpq_mat->triv_pos);
  free(mpq_mat->dense_idx);
  free(mpq_mat->dst);
}

static inline void mpq_matfglm_partial_clear(mpq_matfglm_t mpq_mat) {
  /* entries of mpq_mat->denoms are cleared previously when building mpz_mat */
  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    mpz_clear(mpq_mat->denoms[i]);
  }
  free(mpq_mat->denoms);

  uint64_t nc = mpq_mat->ncols;
  for (uint32_t i = 0; i < mpq_mat->nrows; i++) {
    uint64_t c = 2 * i * mpq_mat->ncols;
    for (uint32_t j = 0; j < nc; j++) {
      mpz_clear(mpq_mat->dense_mat[c + 2 * j]);
    }
  }
  free(mpq_mat->dense_mat);

  free(mpq_mat->triv_idx);
  free(mpq_mat->triv_pos);
  free(mpq_mat->dense_idx);
  free(mpq_mat->dst);
}

static inline void copy_modular_matrix(trace_det_fglm_mat_t trace_det,
        sp_matfglm_t **mod_mat, uint32_t newalloc, mp_limb_t prime){
    deg_t nrows = (*mod_mat)->nrows;
    deg_t ncols = (*mod_mat)->ncols;
    size_t sz = nrows * ncols;
    if(trace_det->num_mat == trace_det->mat_alloc ){
        int32_t old_alloc = trace_det->mat_alloc;
        trace_det->mat_alloc += newalloc;
        size_t sz2 = sz * trace_det->mat_alloc; 

        int32_t *tmp = (int32_t *)realloc(trace_det->modular_matrices, 
                sz2 * sizeof(int32_t));
        if(tmp == NULL){
            fprintf(stderr, "Problem when allocating modular matrices (amount = %ld)\n", trace_det->mat_alloc * sz);
            exit(1);
        }
        for(int64_t i = sz*old_alloc; i < sz2; i++){
            tmp[i] = 0;
        }
        trace_det->modular_matrices = tmp;
        trace_det->primes = flint_realloc(trace_det->primes, 
                trace_det->mat_alloc * sizeof(mp_limb_t));
    }

    for(int64_t i = 0; i < sz; i++){
        trace_det->modular_matrices[(trace_det->num_mat) * sz + i] = (*mod_mat)->dense_mat[i];
    }
    trace_det->num_mat++;
    trace_det->primes[trace_det->num_primes] = prime;
    trace_det->num_primes++;
}


static inline void trace_det_initset(trace_det_fglm_mat_t trace_det,
                                     uint32_t trace_mod, uint32_t det_mod,
                                     uint32_t tridx, uint32_t detidx, 
                                     sp_matfglm_t **mod_mat, 
                                     int lift_matrix, 
                                     int nlins, 
                                     int32_t nv,
                                     uint32_t **lineqs_ptr,
                                     uint32_t newalloc, 
                                     mp_limb_t prime) {
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

  trace_det->lift_matrix = lift_matrix;
  trace_det->nlins = nlins;
  trace_det->nv = nv;
  if(nlins && lift_matrix){
     trace_det->lin_lifted = 0;
  }
  else{
     trace_det->lin_lifted = 2;
  }
  trace_det->crt_linear_forms = allocate_crt_linear_forms(nlins, nv, lineqs_ptr);
  trace_det->mpq_linear_forms = allocate_mpq_linear_forms(nlins, nv);
  trace_det->mpz_linear_forms = mpz_linear_forms_allocate(nlins, nv);

  uint32_t nrows = (*mod_mat)->nrows;
  uint32_t ncols = (*mod_mat)->ncols;

  if(lift_matrix){
      trace_det->triv_idx = (szmat_t*)malloc((ncols - nrows) * sizeof(szmat_t));
      for(szmat_t i = 0; i < ncols - nrows; i++){
          trace_det->triv_idx[i] = (*mod_mat)->triv_idx[i];
      }
      trace_det->triv_pos = (szmat_t*)malloc((ncols - nrows) * sizeof(szmat_t));
      for(szmat_t i = 0; i < ncols - nrows; i++){
          trace_det->triv_pos[i] = (*mod_mat)->triv_pos[i];
      }
      trace_det->dense_idx = (szmat_t*)malloc(nrows * sizeof(szmat_t));
      for(szmat_t i = 0; i < nrows; i++){
          trace_det->dense_idx[i] = (*mod_mat)->dense_idx[i];
      }
      trace_det->dst = (szmat_t*)malloc(nrows * sizeof(szmat_t));
      for(szmat_t i = 0; i < nrows; i++){
          trace_det->dst[i] = (*mod_mat)->dst[i];
      }
  }
  else{
      trace_det->triv_idx = NULL;
      trace_det->dense_idx = NULL;
      trace_det->triv_pos = NULL;
      trace_det->dst = NULL;
  }
  trace_det->nrows = nrows;
  trace_det->ncols = ncols;
  trace_det->nlifted = 0;
  trace_det->w_checked = 0;
  trace_det->matmul_indices = (uint64_t *)malloc(nrows * sizeof(uint64_t));
  for(uint32_t i = 0; i < nrows; i++){
      trace_det->matmul_indices[i] = 0;
      uint64_t tmp = i*((*mod_mat)->ncols);
      for(uint32_t j = 0; j < (*mod_mat)->ncols; j++){
          if((*mod_mat)->dense_mat[tmp+j] != 0){
              trace_det->matmul_indices[i] = tmp + j;
              break;
          }
      }
  }
  trace_det->matmul_wcrt = (mpz_t *)malloc(nrows * sizeof(mpz_t));
  for(uint32_t i = 0; i < nrows; i++){
      mpz_init_set_ui(trace_det->matmul_wcrt[i], (*mod_mat)->dense_mat[trace_det->matmul_indices[i]]);
  }
  trace_det->matmul_wqq = (mpz_t *)malloc(2 * nrows * sizeof(mpz_t));
  for(uint32_t i = 0; i < 2 * nrows; i++){
      mpz_init(trace_det->matmul_wqq[i]);
  }
  trace_det->done_coeffs = (int16_t *)calloc(nrows, sizeof(int16_t));
  trace_det->check_coeffs = (int16_t *)calloc(nrows, sizeof(int16_t));
  if(!lift_matrix){
      trace_det->mat_alloc = 0;
      trace_det->num_mat = 0;
      trace_det->modular_matrices = NULL;
      trace_det->primes = NULL;
      trace_det->num_primes = 0;
      trace_det->residues = NULL;
  }
  else{
      trace_det->mat_alloc = 4;
      trace_det->num_mat = 0;
      trace_det->modular_matrices = (int32_t *)malloc(trace_det->mat_alloc * nrows * ncols * sizeof(int32_t) );
      trace_det->primes = flint_malloc(trace_det->mat_alloc * sizeof(mp_limb_t));
      trace_det->num_primes = 0;
      trace_det->residues = NULL;
      copy_modular_matrix(trace_det, mod_mat, newalloc, prime);
  }
  trace_det->mat_lifted = 0;
  if(lift_matrix){
      trace_det->mat_denoms = (mpz_t *)malloc(nrows * sizeof(mpz_t));
      for(int32_t i = 0; i < nrows; i++){
          mpz_init(trace_det->mat_denoms[i]);
      }
      trace_det->dense_mat = malloc(nrows * ncols * sizeof(mpz_t));
      for(int32_t i = 0; i < nrows; i++){
          int32_t sz = i*ncols;
          for(int32_t j = 0; j < ncols; j++){
              mpz_init(trace_det->dense_mat[sz + j]);
          }
      }
  }
  else{
      trace_det->mat_denoms = NULL;
      trace_det->dense_mat = NULL;
  }
}

static inline void trace_det_clear(trace_det_fglm_mat_t trace_det) {
  int lift_matrix = trace_det->lift_matrix;
  mpz_clear(trace_det->trace_crt);
  mpz_clear(trace_det->det_crt);
  mpz_clear(trace_det->trace_num);
  mpz_clear(trace_det->trace_den);
  mpz_clear(trace_det->det_num);
  mpz_clear(trace_det->det_den);
  mpz_clear(trace_det->tmp);

  int32_t nlins = trace_det->nlins;
  int32_t nv = trace_det->nv;

  if(lift_matrix){
      free(trace_det->triv_idx);
      free(trace_det->dense_idx);
      free(trace_det->triv_pos);
      free(trace_det->dst);
  }
  mpz_linear_forms_clear(trace_det->mpz_linear_forms, nlins, nv);
  crt_linear_forms_clear(trace_det->crt_linear_forms, nlins, nv);
  mpq_linear_forms_clear(trace_det->mpq_linear_forms, nlins, nv);

  uint32_t nrows = trace_det->nrows;
  for(uint32_t i = 0; i < nrows; i++){
      mpz_clear(trace_det->matmul_wcrt[i]);
      mpz_clear(trace_det->matmul_wqq[2 * i]);
      mpz_clear(trace_det->matmul_wqq[2 * i + 1]);
  }
  free(trace_det->matmul_wcrt);
  free(trace_det->matmul_wqq);
  free(trace_det->matmul_indices);
  free(trace_det->done_coeffs);
  free(trace_det->check_coeffs);
  free(trace_det->modular_matrices);
  flint_free(trace_det->primes);
  if(lift_matrix){
      int32_t nrows = trace_det->nrows;
      int32_t ncols = trace_det->ncols;
      for(int32_t i = 0; i < nrows; i++){
          mpz_clear(trace_det->mat_denoms[i]);
      }
      free(trace_det->mat_denoms);
      for(int32_t i = 0; i < nrows; i++){
          int32_t sz = i*ncols;
          for(int32_t j = 0; j < ncols; j++){
              mpz_clear(trace_det->dense_mat[sz + j]);
          }
      }
      free(trace_det->dense_mat);
  }
}

static inline void crt_lift_dense_rows(mpz_t *rows, uint32_t *mod_rows,
                                       const uint64_t start, const uint64_t end,
                                       mpz_t modulus, mpz_t prod, int32_t prime,
                                       mpz_t tmp, const int nthrds) {
  len_t i;
  for (i = start; i < end; i++) {
    mpz_CRT_ui(rows[i], rows[i], modulus, mod_rows[i], prime, prod, tmp, 0);
  }
}

static inline void crt_lift_trace_det(trace_det_fglm_mat_t trace_det,
                                      uint32_t trace_mod, uint32_t det_mod,
                                      sp_matfglm_t *mod_mat,
                                      uint32_t *lineqs, 
                                      param_t *nmod_param,
                                      mpz_t modulus, mpz_t prod,
                                      uint32_t prime, 
                                      int nthrds
                                      ) {
  mpz_CRT_ui(trace_det->trace_crt, trace_det->trace_crt, modulus, trace_mod,
             prime, prod, trace_det->tmp, 0);
  mpz_CRT_ui(trace_det->det_crt, trace_det->det_crt, modulus, det_mod, prime,
             prod, trace_det->tmp, 0);

  int32_t nlins = trace_det->nlins;
  int32_t lin_lifted = trace_det->lin_lifted;

  if (nlins && (lin_lifted) < 2) {
      for (int32_t i = 0; i < nlins; i++) {
          crt_lift_dense_rows(trace_det->crt_linear_forms + i * (nmod_param->nvars + 1),
                          lineqs + i * (nmod_param->nvars + 1), 0,
                          nmod_param->nvars + 1, modulus, prod, prime,
                          trace_det->tmp, nthrds);
      }
  }

  for(uint32_t i = trace_det->w_checked; i < mod_mat->nrows; i++){
      mpz_CRT_ui(trace_det->matmul_wcrt[i], trace_det->matmul_wcrt[i], modulus, 
                 mod_mat->dense_mat[trace_det->matmul_indices[i]], prime, prod, 
                 trace_det->tmp, 0);
  }
}

static inline void crt_lift_mat(crt_mpz_matfglm_t mat, sp_matfglm_t *mod_mat,
                                mpz_t modulus, mpz_t prod_crt,
                                const int32_t prime, mpz_t tmp,
                                const int32_t nrows, const int nthrds) {
  /*assumes prod_crt = modulus * prime */
  const uint64_t sz = mat->nrows * mat->ncols;
  crt_lift_dense_rows(mat->dense_mat, mod_mat->dense_mat, nrows * mat->ncols,
                      sz, modulus, prod_crt, prime, tmp, nthrds);
}

static inline void build_linear_forms(mpz_t *mpz_linear_forms,
                                      mpz_t *mpq_linear_forms, uint32_t sz,
                                      int nlins) {
  mpz_t lcm;
  mpz_init(lcm);
  for (int i = 0; i < nlins; i++) {
    mpz_set_ui(lcm, 1);
    int32_t nc = i * sz;
    int32_t nc2 = i * (sz + 1);
    for (int32_t j = 0; j < sz; j++) {
      mpz_lcm(lcm, lcm, mpq_linear_forms[2 * nc + 2 * j + 1]);
    }
    for (int32_t j = 0; j < sz; j++) {
      mpz_mul(mpz_linear_forms[nc2 + j], mpq_linear_forms[2 * nc + 2 * j], lcm);
      mpz_divexact(mpz_linear_forms[nc2 + j], mpz_linear_forms[nc2 + j],
                   mpq_linear_forms[2 * nc + 2 * j + 1]);
    }
    mpz_set(mpz_linear_forms[nc2 + sz], lcm);
  }
#ifdef DEBUGLIFTMAT
  for (int i = 0; i < nlins; i++) {
    int32_t nc = i * sz;
    int32_t nc2 = i * (sz + 1);
    fprintf(stderr, " ===> (%d, %d) ", nc, sz);
    mpz_out_str(stderr, 10, mpz_linear_forms[nc2 + sz]);
    fprintf(stderr, "\n");
    for (int j = 0; j < sz; j++) {
      mpz_out_str(stderr, 10, mpz_linear_forms[nc2 + j]);
      fprintf(stderr, " / ");
      mpz_out_str(stderr, 10, mpz_linear_forms[nc2 + sz]);
      fprintf(stderr, ", ");
    }
    fprintf(stderr, "\n");
  }
#endif
  mpz_clear(lcm);
}


static inline int old_rat_recon_trace_det(trace_det_fglm_mat_t trace_det,
                                      rrec_data_t recdata, mpz_t modulus,
                                      mpz_t rnum, mpz_t rden, mpz_t gden){
  if(trace_det->done_trace < 2){
    int b = ratrecon(rnum, rden, trace_det->trace_crt, modulus, recdata);
    if(b == 1){
      mpz_set(trace_det->trace_num, rnum);
      mpz_set(trace_det->trace_den, rden);
    }
    else
    {
      return 0;
    }
  }
  if(trace_det->done_det < 2){
    int b = ratrecon(rnum, rden, trace_det->det_crt, modulus, recdata);
    if(b == 1){
      mpz_set(trace_det->det_num, rnum);
      mpz_set(trace_det->det_den, rden);
    }
    else{
      return 0;
    }
  }
  return 1;
}

static inline int rat_recon_array(mpz_t *res, mpz_t *crt, int32_t sz,
                                  mpz_t modulus, rrec_data_t rdata) {
  for (int i = 0; i < sz; i++) {
    int b = ratrecon(res[2 * i], res[2 * i + 1], crt[i], modulus, rdata);
    if (b == 0)
      return 0;
  }
  return 1;
}

static inline int rat_recon_trace_det(trace_det_fglm_mat_t trace_det,
                                      rrec_data_t recdata, mpz_t modulus,
                                      mpz_t rnum, mpz_t rden, mpz_t gden) {
    
    if(trace_det->done_det > 1 && trace_det->done_trace > 1){
        return 1;
    }

    /* rational reconstruction of linear forms */
    mpz_t *mpq_linear_forms = trace_det->mpq_linear_forms;
    mpz_t *crt_linear_forms = trace_det->crt_linear_forms;
    mpz_t *mpz_linear_forms = trace_det->mpz_linear_forms;
    int32_t nv = trace_det->nv;
    if (trace_det->nlins && trace_det->lin_lifted < 2) {
        int boo = 0;
        for (int i = 0; i < trace_det->nlins; i++) {
            boo = rat_recon_array(mpq_linear_forms + (2 * i * (nv + 1)),
                        crt_linear_forms + (i * (nv + 1)),
                        nv + 1, modulus, recdata);
            if (boo == 0)
                break;
        }
        if (boo) {
            build_linear_forms(mpz_linear_forms, mpq_linear_forms,
                       nv + 1, trace_det->nlins);
            trace_det->lin_lifted = 1;
        }
    }
    if(trace_det->nlins && trace_det->lin_lifted < 2){
        return 0;
    }

    /* rational reconstruction of witness coefficients for mul mat */
    /* not applied once the multiplication matrix is lifted */
    mpz_t gcd;
    mpz_init(gcd);
    int nlifted = trace_det->w_checked;
    for(int32_t i = trace_det->w_checked; i < trace_det->nrows; i++){
        int b = ratreconwden(rnum, rden, trace_det->matmul_wcrt[i], modulus, gden, recdata);
        if(b == 1){
            mpz_set(trace_det->matmul_wqq[2*i], rnum);
            mpz_set(trace_det->matmul_wqq[2*i + 1], rden);
            mpz_mul(trace_det->matmul_wqq[2*i + 1], trace_det->matmul_wqq[2*i + 1],
                    gden);
            mpz_gcd(gcd, trace_det->matmul_wqq[2*i], trace_det->matmul_wqq[2*i + 1]);
            mpz_divexact(trace_det->matmul_wqq[2*i + 1], 
                    trace_det->matmul_wqq[2*i + 1], gcd);
            mpz_divexact(trace_det->matmul_wqq[2*i], 
                    trace_det->matmul_wqq[2*i], gcd);
            nlifted++;
        }
        else{
            break;
        }
    }
    trace_det->nlifted = nlifted;

    /* rational reconstruction of witness coefficients in eliminating polynomial
     * */
    int trace_rec = 0;
    if(trace_det->done_trace < 2){
    int b =
        ratreconwden(rnum, rden, trace_det->trace_crt, modulus, gden, recdata);
    if (b == 1) {
      mpz_set(trace_det->trace_num, rnum);
      mpz_set(trace_det->trace_den, rden);

      mpz_mul(trace_det->trace_den, trace_det->trace_den, gden);
      mpz_gcd(gcd, trace_det->trace_num, trace_det->trace_den);
      mpz_divexact(trace_det->trace_num, trace_det->trace_num, gcd);
      mpz_divexact(trace_det->trace_den, trace_det->trace_den, gcd);
      trace_det->done_trace = 1;
      trace_rec = 1;
    } else {
        trace_det->done_trace = 0;
        mpz_clear(gcd);
        return 0;
    }
    }

    if(trace_det->done_det < 2 && trace_det->done_trace){
        mpz_t oldgden;
        if(trace_rec){
            mpz_init(oldgden);
            mpz_set(oldgden, gden);
            mpz_set(gden, trace_det->trace_den);
        }
        int b = ratreconwden(rnum, rden, trace_det->det_crt, modulus, gden, recdata);

        if (b == 1) {
          mpz_set(trace_det->det_num, rnum);
          mpz_set(trace_det->det_den, rden);
          mpz_mul(trace_det->det_den, trace_det->det_den, gden);
          mpz_gcd(gcd, trace_det->det_num, trace_det->det_den);
          mpz_divexact(trace_det->det_num, trace_det->det_num, gcd);
          mpz_divexact(trace_det->det_den, trace_det->det_den, gcd);
          trace_det->done_det = 1;
          if(trace_rec){
              mpz_clear(oldgden);
          }
        } else {
          mpz_set_ui(gden, 1);
          if(trace_rec){
              mpz_clear(oldgden);
          }
          trace_det->done_det = 0;
          mpz_clear(gcd);
          return 0;
         }
    }
    mpz_clear(gcd);
    return 1;
}

static inline int check_trace(trace_det_fglm_mat_t trace_det,
                              uint32_t trace_mod, uint32_t prime) {

  uint32_t lc = mpz_fdiv_ui(trace_det->trace_den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(trace_det->trace_num, prime);
  c *= lc;
  c = c % prime;

  return (c == trace_mod);
}

static inline int check_det(trace_det_fglm_mat_t trace_det, uint32_t det_mod,
                            uint32_t prime) {

  uint32_t lc = mpz_fdiv_ui(trace_det->det_den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(trace_det->det_num, prime);
  c *= lc;
  c = c % prime;

  return (c == det_mod);
}

static inline int check_lifted_coeff(mpz_t num, mpz_t den, 
        uint32_t mod, uint32_t prime){
  uint32_t lc = mpz_fdiv_ui(den, prime);
  lc = mod_p_inverse_32(lc, prime);

  uint64_t c = mpz_fdiv_ui(num, prime);
  c *= lc;
  c = c % prime;

  return (c == mod);
}

static inline int check_trace_det_data(trace_det_fglm_mat_t trace_det, 
        deg_t *maxrec, 
        uint32_t trace_mod, 
        uint32_t det_mod, 
        uint32_t *modular_linear_forms,
        sp_matfglm_t *mod_mat, 
        int32_t prime, 
        int lift_matrix){
  //checks if trace is lifted
  if (trace_det->done_trace && trace_det->done_trace < 2) {
    if (check_trace(trace_det, trace_mod, prime)) {
        trace_det->done_trace++;
        *maxrec = trace_det->det_idx;
    }
    else{
        trace_det->done_trace = 0;
    }
  }
  //checks if det is lifted
  if (trace_det->done_det && trace_det->done_det < 2) {
    if (check_det(trace_det, det_mod, prime)) {
        trace_det->done_det++;
        *maxrec = trace_det->det_idx;
    }
    else{
        trace_det->done_det = 0;
    }
  }
  if(lift_matrix && trace_det->nlins && trace_det->lin_lifted < 2){
      if(trace_det->lin_lifted < 2){
          nvars_t sz = trace_det->nv + 1;
          for(int i = 0; i < trace_det->nlins; i++){
              nvars_t n = i * sz;
              nvars_t n2 = i * (sz + 1);
              for(nvars_t j = 0; j < sz; j++){
                  if(!check_lifted_coeff(trace_det->mpz_linear_forms[n2+j], 
                             trace_det->mpz_linear_forms[n2+sz], 
                             modular_linear_forms[n+j], 
                             prime)){
                      trace_det->lin_lifted = 0;
                      return 0;
                  }
              }
          }
          trace_det->lin_lifted = 2;
      }
  }
  if(lift_matrix && trace_det->w_checked < trace_det->nlifted){
      int64_t sz = trace_det->nrows * trace_det->ncols * (trace_det->num_mat - 1);
      int32_t nr = trace_det->w_checked;
      for(deg_t i = trace_det->w_checked; i < trace_det->nlifted; i++){
          if(check_lifted_coeff(trace_det->matmul_wqq[2*i], 
                  trace_det->matmul_wqq[2*i+1], 
                  trace_det->modular_matrices[sz + trace_det->matmul_indices[i]], prime)){
              trace_det->done_coeffs[i]++;
              if(trace_det->done_coeffs[i] > 0){
                 nr++;
              }
          }
          else{
             trace_det->done_coeffs[i] = 0;
             break;
          }
      }
      return nr;
  }
  return 0;
}

static inline int check_linear_forms(trace_det_fglm_mat_t trace_det, 
        uint32_t *modular_linear_forms,
        int32_t prime){
    int32_t nlins = trace_det->nlins;
    int32_t nv = trace_det ->nv;
    nvars_t sz = nv + 1;
    for(int32_t i = 0; i < nlins; i++){
        nvars_t n = i *sz;
        nvars_t n2 = i * (sz + 1);
        for(int32_t j = 0; j < sz; j++){
            int b = check_lifted_coeff(trace_det->mpz_linear_forms[n2 + j], 
                    trace_det->mpz_linear_forms[n2+sz],
                    modular_linear_forms[n], 
                    prime);
            if(!b) return 0;
        }
    }
}

static inline void check_matrix(trace_det_fglm_mat_t trace_det, 
        sp_matfglm_t *mod_mat, 
        int32_t prime){
    int32_t nrows = trace_det->nrows;
    int32_t ncols = trace_det->ncols;
    for(int32_t row = 0; row < nrows; row++){
        int32_t sz = row * ncols;
        for(int32_t col = 0; col < ncols; col++){
            int b = check_lifted_coeff(trace_det->dense_mat[sz + col], 
                    trace_det->mat_denoms[row], 
                    mod_mat->dense_mat[sz + col], 
                    prime);
            if(!b){
                trace_det->mat_lifted = 0;
                return;
            }
        }
    }
    trace_det->mat_lifted = 2;
}


static inline void mulmat_reconstruct(trace_det_fglm_mat_t trace_det, 
        mpz_t modulus){
    int32_t num_primes = trace_det->num_primes;
    fmpz_comb_init(trace_det->comb, trace_det->primes, num_primes);
    fmpz_comb_temp_init(trace_det->comb_temp, trace_det->comb);

    trace_det->residues = flint_malloc(num_primes * sizeof(mp_limb_t));
    int32_t nrows = trace_det->nrows;
    int32_t ncols = trace_det->ncols;
    int32_t sz = nrows * ncols;

    fmpz_t y; 
    fmpz_init(y);
    fmpz_t fmodulus;
    fmpz_init(fmodulus);
    fmpq_t res;
    fmpq_init(res);

    fmpz_set_mpz(fmodulus, modulus);
    mpz_t lcm;
    mpz_init(lcm);
    mpz_set_ui(lcm, 1);

    mpq_t rr;
    mpq_init(rr);
    int32_t row = 0;
    int32_t col = 0;
    mpz_t *rowdenoms = (mpz_t *)malloc(ncols * sizeof(mpz_t));
    for(int32_t i = 0; i < ncols; i++){
        mpz_init(rowdenoms[i]);
    }
    int boo = 0;
    for(int32_t elt = 0 ; elt < sz; elt++){
        for(int32_t i = 0; i < num_primes; i++){
            trace_det->residues[i] = (mp_limb_t)trace_det->modular_matrices[i*sz + elt];
        }
        fmpz_multi_CRT_ui(y, trace_det->residues, 
                trace_det->comb, trace_det->comb_temp, 0);
        boo = fmpq_reconstruct_fmpz(res, y, fmodulus);
        if(boo == 0) {
            break;
        }
        fmpq_get_mpq(rr, res);
        mpz_lcm(lcm, lcm, mpq_denref(rr));
        mpz_set(trace_det->dense_mat[elt], mpq_numref(rr));
        mpz_set(rowdenoms[col], mpq_denref(rr));
        col++;
        if(col == ncols){
            mpz_set(trace_det->mat_denoms[row], lcm);
            mpz_set_ui(lcm, 1);
            for(int32_t nc = 0; nc < ncols; nc++){
                mpz_divexact(rowdenoms[nc], trace_det->mat_denoms[row], 
                        rowdenoms[nc]);
                mpz_mul(trace_det->dense_mat[row*ncols + nc], 
                        trace_det->dense_mat[row*ncols + nc], 
                        rowdenoms[nc]);
            }
            row++;
            col = 0;
        }
    }

    for(int32_t i = 0; i < nrows; i++){
        mpz_clear(rowdenoms[i]);
    }
    free(rowdenoms);
    fmpz_comb_temp_clear(trace_det->comb_temp);
    fmpz_comb_clear(trace_det->comb);
    flint_free(trace_det->residues);
    fmpz_clear(y);
    fmpq_clear(res);
    fmpz_clear(fmodulus);
    mpz_clear(lcm);
    mpq_clear(rr);
    if(boo){
        trace_det->mat_lifted = 1;
    }
}

static inline void
rat_recon_matfglm(mpq_matfglm_t mpq_mat, crt_mpz_matfglm_t crt_mat,
                  mpz_t modulus, rrec_data_t rdata,
                  mpz_t rnum, mpz_t rden, mpz_t *guessed_den, deg_t *matrec, int *mat_lifted) {
  const uint32_t nrows = crt_mat->nrows;
  const uint32_t ncols = crt_mat->ncols;

  mpz_t lcm, coef;
  mpz_init(lcm);
  mpz_init(coef);
  mpz_t gden;
  mpz_init(gden);
  if((0==1) && (*matrec)>1){
      mpz_set(gden, *guessed_den);
  }
  else{
      mpz_set_ui(gden, 1);
  }
  for (uint32_t i = *matrec; i < nrows; i++) {
    uint64_t c = i * ncols;
    mpz_set_ui(lcm, 1);

    for (uint32_t j = 0; j < ncols; j++) {
      int b = 0;
      mpz_mul(coef, crt_mat->dense_mat[c + j], lcm);
      mpz_mod(coef, coef, modulus);
//      b = ratrecon(rnum, rden, coef, modulus, rdata);
      b = ratreconwden(rnum, rden, coef, modulus, gden, rdata);
      if (b == 1) {
        mpz_set(mpq_mat->dense_mat[2 * c + 2 * j], rnum);
        mpz_set(mpq_mat->dense_mat[2 * c + 2 * j + 1], rden);
        mpz_mul(mpq_mat->dense_mat[2 * c + 2 * j + 1],
                mpq_mat->dense_mat[2 * c + 2 * j + 1], lcm);
        mpz_mul(mpq_mat->dense_mat[2 * c + 2 * j + 1],
                mpq_mat->dense_mat[2 * c + 2 * j + 1], gden);
        mpz_lcm(lcm, lcm, rden);
      } else {
          mpz_clear(gden);
        mpz_clear(lcm);
        mpz_clear(coef);
        return ;
      }
    }
    mpz_set_ui(mpq_mat->denoms[i], 1);
    for (uint32_t j = 0; j < ncols; j++) {
      mpz_lcm(mpq_mat->denoms[i], mpq_mat->denoms[i],
              mpq_mat->dense_mat[2 * c + 2 * j + 1]);
    }

    for(uint32_t j = 0; j < ncols; j++){
        mpz_mul(mpq_mat->dense_mat[2 * c+ 2 * j], mpq_mat->dense_mat[2 * c+ 2 * j], 
                mpq_mat->denoms[i]);
        mpz_divexact(mpq_mat->dense_mat[2 * c + 2 * j], mpq_mat->dense_mat[2 * c + 2 * j], 
                mpq_mat->dense_mat[2 * c + 2 * j + 1]);
    }
    (*matrec)++;
  }
  mpz_clear(coef);
  mpz_clear(lcm);
  mpz_clear(gden);

  *mat_lifted = 1;
}


static inline void check_matrix_and_linear_forms(mpq_matfglm_t mpq_mat, sp_matfglm_t *mod_mat, 
        mpz_t *mpz_lin, uint32_t *mod_lin, const int nlins, const nvars_t nv,
        int *mat_lifted, int* lin_lifted, deg_t *oldmatrec_checked, deg_t *matrec_checked, const uint32_t prime){
    *oldmatrec_checked = *matrec_checked;
    /* checks lifted linear forms */
    nvars_t sz = nv + 1;
    uint32_t cand;
    if(*lin_lifted < 2){
    for(int i = 0; i < nlins; i++){
        nvars_t n = i * sz;
        nvars_t n2 = i * (sz + 1);
        uint64_t inv = mpz_fdiv_ui(mpz_lin[n2 + sz], prime);
        inv = mod_p_inverse_32(inv, prime);
        for(nvars_t j = 0; j < sz; j++){
            cand = (((uint64_t)mpz_fdiv_ui(mpz_lin[n2+j], prime)) * inv) % prime;
            if(cand != mod_lin[n + j]){
                *lin_lifted = 0;
                *matrec_checked = 0;
                return;
            }
        }
    }
    *lin_lifted = 2;
    }

    /*checks multiplication matrix */
    mod_mat->charac = prime;

    for(int32_t i = (*oldmatrec_checked); i < mpq_mat->nrows; i++){
      uint64_t lc = mpz_fdiv_ui(mpq_mat->denoms[i], prime);
      lc = mod_p_inverse_32(lc, prime);
      uint32_t nc = i*mpq_mat->ncols;
      for(uint32_t j = 0 ; j < mpq_mat->ncols; j++){
        uint32_t mod = mpz_fdiv_ui(mpq_mat->dense_mat[2*(nc+j)], prime);
        if(mod_mat->dense_mat[nc+j] != ( ((uint64_t)mod) * lc )% prime){
            *mat_lifted = 0;
            *matrec_checked = MAX(0, i);
            return;
        }
      }
      *matrec_checked = i + 1;
    }
    *mat_lifted = 2;

}

void initialize_rrec_data(rrec_data_t recdata) {
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

void free_rrec_data(rrec_data_t recdata) {
  mpz_clear(recdata->r0);
  mpz_clear(recdata->r1);
  mpz_clear(recdata->t0);
  mpz_clear(recdata->t1);
  mpz_clear(recdata->q);
  mpz_clear(recdata->tmp);
  mpz_clear(recdata->N);
  mpz_clear(recdata->D);
}
