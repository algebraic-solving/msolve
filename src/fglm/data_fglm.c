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

#include<stdint.h>
#include <flint/flint.h>
#include <flint/longlong.h>
#include <flint/mpn_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>

#include "aligned_alloc.h"

static inline void free_sp_mat_fglm(sp_matfglm_t *mat){
  if(mat!=NULL){
    free(mat->dense_mat);
    free(mat->triv_idx);
    free(mat->triv_pos);
    free(mat->dense_idx);
    free(mat->dst);
    free(mat);
  }
}

static inline fglm_data_t *allocate_fglm_data(szmat_t nrows, szmat_t ncols, szmat_t nvars){
  fglm_data_t * data = malloc(sizeof(fglm_data_t));

  szmat_t block_size = nvars; //block size in data->res

  if(posix_memalign((void **)&data->vecinit, 32, ncols*sizeof(CF_t))){
    fprintf(stderr, "posix_memalign failed\n");
    exit(1);
  }

  if(posix_memalign((void **)&data->res, 32, 2 * block_size * ncols * sizeof(CF_t))){
    fprintf(stderr, "posix_memalign failed\n");
    exit(1);
  }

  if(posix_memalign((void **)&data->vecmult, 32, nrows*sizeof(CF_t))){
    fprintf(stderr, "posix_memalign failed\n");
    exit(1);
  }

  if(posix_memalign((void **)&data->vvec, 32, ncols*sizeof(CF_t))){
    fprintf(stderr, "posix_memalign failed\n");
    exit(1);
  }
  data->pts = calloc(ncols * 2, sizeof(mp_limb_t));

  for(szmat_t i = 0; i < 2*block_size*ncols; i++){
    data->res[i] = 0;
  }
  for(szmat_t i = 0; i < nrows; i++){
    data->vecmult[i] = 0;
  }
  for(szmat_t i = 0; i < ncols; i++){
    data->vvec[i] = 0;
    data->vecinit[i] = 0;
  }

  return data;
}


static inline void free_fglm_data(fglm_data_t *data){
  posix_memalign_free(data->vecinit);
  posix_memalign_free(data->res);
  posix_memalign_free(data->vecmult);
  posix_memalign_free(data->vvec);
  free(data->pts);
  free(data);
}

static inline void display_fglm_matrix(FILE *file, sp_matfglm_t *matrix){

  fprintf(file, "%u\n", matrix->charac);
  fprintf(file, "%u\n", matrix->ncols);
  fprintf(file, "%u\n", matrix->nrows);

  szmat_t len1 = (matrix->ncols)*(matrix->nrows);
  for(szmat_t i = 0; i < len1; i++){
    fprintf(file, "%d ", matrix->dense_mat[i]);
  }
  fprintf(file, "\n");
  szmat_t len2 = (matrix->ncols) - (matrix->nrows);
  for(szmat_t i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(szmat_t i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(szmat_t i = 0; i < matrix->nrows; i++){
    fprintf(file, "%d ", matrix->dense_idx[i]);
  }
  fprintf(file, "\n");
}

static inline void display_fglm_colon_matrix(FILE *file, sp_matfglmcol_t *matrix){

  fprintf(file, "%u\n", matrix->charac);
  fprintf(file, "%u\n", matrix->ncols);
  fprintf(file, "%u\n", matrix->nrows);
  fprintf(file, "%u\n", matrix->nzero);

  szmat_t len1 = (matrix->ncols)*(matrix->nrows);
  for(szmat_t i = 0; i < len1; i++){
    fprintf(file, "%d ", matrix->dense_mat[i]);
  }
  fprintf(file, "\n");
  szmat_t len2 = (matrix->ncols) - (matrix->nrows);
  for(szmat_t i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(szmat_t i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(szmat_t i = 0; i < matrix->nrows; i++){
    fprintf(file, "%d ", matrix->dense_idx[i]);
  }
  fprintf(file, "\n");
  for(szmat_t i = 0; i < matrix->nzero; i++){
    fprintf(file, "%d ", matrix->zero_idx[i]);
  }
  fprintf(file, "\n");
}

static inline param_t *allocate_fglm_param(mp_limb_t prime, long nvars){
  param_t *param = malloc(sizeof(param_t));
  if(param==NULL){
    fprintf(stderr, "Pb when calling malloc to allocate param_t\n");
    exit(1);
    return param;
  }
  param->charac = prime;
  param->nvars = nvars;
  nmod_poly_init(param->elim, prime);
  nmod_poly_init(param->denom, prime);
  param->coords = malloc(sizeof(nmod_poly_t) * (nvars-1));
  for(nvars_t i = 0; i < nvars-1; i++){
    nmod_poly_init(param->coords[i], prime);
  }
  return param;
}

static inline void free_fglm_param(param_t *param){
  nmod_poly_clear(param->elim);
  nmod_poly_clear(param->denom);
  for(szmat_t i = 0; i < param->nvars-1; i++){
    nmod_poly_clear(param->coords[i]);
  }
  free(param->coords);
  free(param);
}

static inline fglm_bms_data_t *allocate_fglm_bms_data(szmat_t dim, mp_limb_t prime){

  fglm_bms_data_t * data_bms = (fglm_bms_data_t *)malloc(sizeof(fglm_bms_data_t));
  nmod_poly_init(data_bms->A, prime);

  nmod_poly_init(data_bms->B, prime);
  nmod_poly_init(data_bms->Z1, prime);

  nmod_poly_init2(data_bms->rZ1, prime, dim+1);

  nmod_poly_init(data_bms->Z2, prime);
  nmod_poly_init2(data_bms->rZ2, prime, dim+1);
  nmod_poly_init2(data_bms->V, prime, dim+1);

  nmod_poly_init2(data_bms->param, prime, dim+1);

  for(len_t i = 0; i < dim + 1; i++){
    data_bms->rZ1->coeffs[i] = 0;
    data_bms->rZ2->coeffs[i] = 0;
    data_bms->V->coeffs[i] = 0;
    data_bms->param->coeffs[i] = 0;
  }

  nmod_berlekamp_massey_init(data_bms->BMS, (mp_limb_t)prime);

  nmod_poly_factor_init(data_bms->sqf);
  return data_bms;
}

static inline void nmod_poly_set_prime(nmod_poly_t poly,
                                       mp_limb_t prime){
  mp_limb_t ninv = n_preinvert_limb(prime);
  poly->mod.n = prime;
  poly->mod.ninv = ninv;
#if __FLINT_VERSION < 3
  count_leading_zeros(poly->mod.norm, prime);
#else
  poly->mod.norm = flint_clz(prime);
#endif
}

static inline void fglm_param_set_prime(param_t *param, mp_limb_t prime){
  param->charac = prime;

  nmod_poly_set_prime(param->elim, prime);
  nmod_poly_set_prime(param->denom, prime);
  for(nvars_t i = 0; i < param->nvars-1; i++){
    nmod_poly_set_prime(param->coords[i], prime);
  }

}


static inline void fglm_bms_data_set_prime(fglm_bms_data_t *data_bms,
                                           mp_limb_t prime){
  nmod_poly_set_prime(data_bms->A, prime);
  nmod_poly_set_prime(data_bms->B, prime);
  nmod_poly_set_prime(data_bms->Z1, prime);
  nmod_poly_set_prime(data_bms->rZ1, prime);
  nmod_poly_set_prime(data_bms->Z2, prime);
  nmod_poly_set_prime(data_bms->rZ2, prime);
  nmod_poly_set_prime(data_bms->V, prime);
  nmod_poly_set_prime(data_bms->param, prime);

  nmod_berlekamp_massey_set_prime(data_bms->BMS, prime);

}

static inline void free_fglm_bms_data(fglm_bms_data_t *data_bms){
  nmod_poly_clear(data_bms->A);
  nmod_poly_clear(data_bms->B);
  nmod_poly_clear(data_bms->Z1);
  nmod_poly_clear(data_bms->Z2);
  nmod_poly_clear(data_bms->rZ1);
  nmod_poly_clear(data_bms->rZ2);
  nmod_poly_clear(data_bms->V);
  nmod_poly_clear(data_bms->param);
  nmod_poly_factor_clear(data_bms->sqf);

  nmod_berlekamp_massey_clear(data_bms->BMS);

  free(data_bms);
}

