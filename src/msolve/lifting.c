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


static inline void trace_det_initset(trace_det_fglm_mat_t trace_det,
                                     uint32_t trace_mod, uint32_t det_mod,
                                     uint32_t tridx, uint32_t detidx, 
                                     sp_matfglm_t **mod_mat) {
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

  uint32_t nrows = (*mod_mat)->nrows;
  trace_det->nrows = nrows;
  trace_det->nlifted = 0;
  trace_det->matmul_indices = (uint64_t *)malloc(nrows * sizeof(uint64_t));
  for(uint32_t i = 0; i < nrows; i++){
      uint64_t tmp = i*((*mod_mat)->ncols);
      for(uint32_t j = 0; j < (*mod_mat)->ncols; j++){
          if((*mod_mat)->dense_mat[tmp+j] != 0){
              trace_det->matmul_indices[i] = tmp + j;
              break;
          }
      }
  }
  trace_det->matmul_crt = (mpz_t *)malloc(nrows * sizeof(mpz_t));
  for(uint32_t i = 0; i < nrows; i++){
      mpz_init_set_ui(trace_det->matmul_crt[i], (*mod_mat)->dense_mat[trace_det->matmul_indices[i]]);
  }
  trace_det->matmul_cfs_qq = (mpz_t *)malloc(2 * nrows * sizeof(mpz_t));
  for(uint32_t i = 0; i < 2 * nrows; i++){
      mpz_init(trace_det->matmul_cfs_qq[i]);
  }
  trace_det->done_coeffs = (int16_t *)calloc(nrows, sizeof(int16_t));
  trace_det->check_coeffs = (int16_t *)calloc(nrows, sizeof(int16_t));
}

static inline void trace_det_clear(trace_det_fglm_mat_t trace_det) {
  mpz_clear(trace_det->trace_crt);
  mpz_clear(trace_det->det_crt);
  mpz_clear(trace_det->trace_num);
  mpz_clear(trace_det->trace_den);
  mpz_clear(trace_det->det_num);
  mpz_clear(trace_det->det_den);
  mpz_clear(trace_det->tmp);

  uint32_t nrows = trace_det->nrows;
  for(uint32_t i = 0; i < nrows; i++){
      mpz_clear(trace_det->matmul_crt[i]);
      mpz_clear(trace_det->matmul_cfs_qq[2 * i]);
      mpz_clear(trace_det->matmul_cfs_qq[2 * i + 1]);
  }
  free(trace_det->matmul_crt);
  free(trace_det->matmul_cfs_qq);
  free(trace_det->matmul_indices);
  free(trace_det->done_coeffs);
  free(trace_det->check_coeffs);
}

static inline void crt_lift_trace_det(trace_det_fglm_mat_t trace_det,
                                      uint32_t trace_mod, uint32_t det_mod,
                                      mpz_t modulus, mpz_t prod,
                                      uint32_t prime) {
  mpz_CRT_ui(trace_det->trace_crt, trace_det->trace_crt, modulus, trace_mod,
             prime, prod, trace_det->tmp, 1);
  mpz_CRT_ui(trace_det->det_crt, trace_det->det_crt, modulus, det_mod, prime,
             prod, trace_det->tmp, 1);
}

static inline void crt_lift_dense_rows(mpz_t *rows, uint32_t *mod_rows,
                                       const uint64_t start, const uint64_t end,
                                       mpz_t modulus, mpz_t prod, int32_t prime,
                                       mpz_t tmp, const int nthrds) {
  len_t i;
  for (i = start; i < end; i++) {
    mpz_CRT_ui(rows[i], rows[i], modulus, mod_rows[i], prime, prod, tmp, 1);
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

static inline int rat_recon_trace_det(trace_det_fglm_mat_t trace_det,
                                      rrec_data_t recdata, mpz_t modulus,
                                      mpz_t rnum, mpz_t rden, mpz_t gden) {
  int b =
      ratreconwden(rnum, rden, trace_det->trace_crt, modulus, gden, recdata);
  /* int b = ratrecon(rnum, rden, trace_det->trace_crt, modulus, recdata); */
  if (b == 1) {
    mpz_t gcd;
    mpz_init(gcd);
    mpz_set(trace_det->trace_num, rnum);
    mpz_set(trace_det->trace_den, rden);
    mpz_mul(trace_det->trace_den, trace_det->trace_den, gden);
    mpz_gcd(gcd, trace_det->trace_num, trace_det->trace_den);
    mpz_divexact(trace_det->trace_num, trace_det->trace_num, gcd);
    mpz_divexact(trace_det->trace_den, trace_det->trace_den, gcd);
    mpz_clear(gcd);
  } else {
    return 0;
  }
  b = ratreconwden(rnum, rden, trace_det->det_crt, modulus, gden, recdata);

  if (b == 1) {
    mpz_t gcd;
    mpz_init(gcd);
    mpz_set(trace_det->det_num, rnum);
    mpz_set(trace_det->det_den, rden);
    mpz_mul(trace_det->det_den, trace_det->det_den, gden);
    mpz_gcd(gcd, trace_det->det_num, trace_det->det_den);
    mpz_divexact(trace_det->det_num, trace_det->det_num, gcd);
    mpz_divexact(trace_det->det_den, trace_det->det_den, gcd);
    mpz_clear(gcd);
  } else {
    return 0;
  }

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

static inline int rat_recon_array(mpz_t *res, mpz_t *crt, int32_t sz,
                                  mpz_t modulus, rrec_data_t rdata) {
  for (int i = 0; i < sz; i++) {
    int b = ratrecon(res[2 * i], res[2 * i + 1], crt[i], modulus, rdata);
    if (b == 0)
      return 0;
  }
  return 1;
}

static inline void
rat_recon_matfglm(mpq_matfglm_t mpq_mat, crt_mpz_matfglm_t crt_mat,
                  mpz_t modulus, rrec_data_t rdata,
                  mpz_t rnum, mpz_t rden, long *matrec, int *mat_lifted) {
  const uint32_t nrows = crt_mat->nrows;
  const uint32_t ncols = crt_mat->ncols;

  mpz_t lcm, coef;
  mpz_init(lcm);
  mpz_init(coef);

  for (uint32_t i = *matrec; i < nrows; i++) {
    uint64_t c = i * ncols;
    mpz_set_ui(lcm, 1);

    for (uint32_t j = 0; j < ncols; j++) {
      int b = 1;
      mpz_mul(coef, crt_mat->dense_mat[c + j], lcm);
      mpz_mod(coef, coef, modulus);
      b = ratrecon(rnum, rden, coef, modulus, rdata);
      if (b == 1) {
        mpz_set(mpq_mat->dense_mat[2 * c + 2 * j], rnum);
        mpz_set(mpq_mat->dense_mat[2 * c + 2 * j + 1], rden);
        mpz_mul(mpq_mat->dense_mat[2 * c + 2 * j + 1],
                mpq_mat->dense_mat[2 * c + 2 * j + 1], lcm);
        mpz_lcm(lcm, lcm, rden);
      } else {
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
        mpz_clear(mpq_mat->dense_mat[2 * c + 2 * j + 1]);
    }
    (*matrec)++;
  }
  crt_mpz_matfglm_clear(crt_mat);
  mpz_clear(coef);
  mpz_clear(lcm);

  *mat_lifted = 1;
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
