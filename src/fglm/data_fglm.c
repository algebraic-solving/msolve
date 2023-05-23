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


typedef uint32_t szmat_t;
typedef uint32_t CF_t;
typedef uint64_t CF_l_t;
typedef uint32_t mod_t;
/* typedef __uint128_t CF_L_t; */

/**

M:=Matrix([[0, 1, 0, 0, 0],[ 0, 0, 0, 1, 0],[1, 2, 4, 5, 6],[0, 0, 0, 0, 1],[7, 8, 1, 9, 4]]);

   Une matrice de la forme

   0 1 0 0 0
   0 0 0 1 0
   1 2 4 5 6
   0 0 0 0 1
   7 8 1 9 4

   est code de la forme
   ncols  = 5
   nrows = 2
   dense_mat -> 1 2 4 5 6 7 8 1 9
   triv_idx -> 0 1 3 (lgr donnee par ncols-nrows)
   triv_pos -> 1 3 4
   dense_idx -> 2 4 (lgr donnee par nrows)
 **/

typedef struct{
  CF_t charac;
  szmat_t ncols; //dimension du quotient
  szmat_t nrows; //nbre de lignes non triviales
  CF_t *dense_mat; // matrice nrows lignes et ncols colonnes (elements donnes par lignes)
  szmat_t *triv_idx; //tableau d'indices des lignes ne contenant que des 0 et un 1
  szmat_t *triv_pos; //position des 1
  szmat_t *dense_idx; //position des lignes non triviales (qui constituent donc
                      //dense_mat)
  szmat_t *dst; //pour la gestion des lignes "denses" mais avec un bloc de zero a la fin
} sp_matfglm_t;

typedef struct{
  CF_t charac;
  szmat_t ncols; //dimension du sev du quotient
  szmat_t nrows; //nbre de lignes non triviales
  szmat_t nzero; //nbre de lignes nulles
  CF_t *dense_mat; // matrice nrows lignes et ncols colonnes (elements donnes par lignes)
  szmat_t *triv_idx; //tableau d'indices des lignes ne contenant que des 0 et un 1
  szmat_t *triv_pos; //position des 1
  szmat_t *dense_idx; //position des lignes non triviales (qui constituent donc
                      //dense_mat)
  szmat_t *zero_idx; //tableau d'indices des lignes ne contenant que des 0
  szmat_t *dst; //pour la gestion des lignes "denses" mais avec un bloc de zero a la fin
} sp_matfglmcol_t;


#ifndef ALIGNED32
#define ALIGNED32 __attribute__((aligned(32)))
#endif
typedef struct{
  CF_t *vecinit ALIGNED32; //vecteur initial
  CF_t *res ALIGNED32; //va contenir les termes de suites qui nous interessent
  CF_t *vecmult ALIGNED32; //utilise pour la multiplication
  CF_t *vvec ALIGNED32;
  CF_l_t *vec_cache ALIGNED32; //useless
  mp_limb_t *pts;
} fglm_data_t;

/* typedef struct { */
/*   slong npoints; */
/*   nmod_poly_t R0, R1; */
/*   nmod_poly_t V0, V1; */
/*   nmod_poly_t qt, rt; //temporaries */
/*   nmod_poly_t points; */
/* } nmod_berlekamp_massey_struct; */
/* typedef nmod_berlekamp_massey_struct nmod_berlekamp_massey_t[1]; */


typedef struct{
  nmod_berlekamp_massey_t BMS;
  nmod_poly_t Z1;
  nmod_poly_t Z2;
  nmod_poly_t rZ1;
  nmod_poly_t rZ2;
  nmod_poly_t A;
  nmod_poly_t B;
  nmod_poly_t V;
  nmod_poly_t param;
  nmod_poly_factor_t sqf;
} fglm_bms_data_t;


typedef struct{
  mp_limb_t charac;
  long nvars;
  nmod_poly_t elim;
  nmod_poly_t denom;
  nmod_poly_t *coords;
} param_t;


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

static inline fglm_data_t *allocate_fglm_data(long nrows, long ncols, long nvars){
  fglm_data_t * data = malloc(sizeof(fglm_data_t));

  szmat_t block_size = nvars; //taille de bloc dans data->res

  if(posix_memalign((void **)&data->vecinit, 32, ncols*sizeof(CF_t))){
    fprintf(stderr, "posix_memalign failed\n");
    exit(1);
  }

  if(posix_memalign((void **)&data->res, 32, block_size * ncols * sizeof(CF_t)*2)){
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
  data->pts = malloc(sizeof(mp_limb_t) * ncols * 2);

  memset(data->res, 0, (block_size)*(ncols)*sizeof(CF_t)*2);
  memset(data->vecinit, 0, (ncols)*sizeof(CF_t));
  memset(data->vecmult, 0, (nrows)*sizeof(CF_t));
  memset(data->vvec, 0, (ncols)*sizeof(CF_t));

  return data;
}


static inline void free_fglm_data(fglm_data_t *data){
  free(data->vecinit);
  free(data->res);
  free(data->vecmult);
  free(data->vvec);
  free(data->pts);
  free(data);
}

static inline void display_fglm_matrix(FILE *file, sp_matfglm_t *matrix){

  fprintf(file, "%u\n", matrix->charac);
  fprintf(file, "%u\n", matrix->ncols);
  fprintf(file, "%u\n", matrix->nrows);

  long len1 = (matrix->ncols)*(matrix->nrows);
  for(long i = 0; i < len1; i++){
    fprintf(file, "%d ", matrix->dense_mat[i]);
  }
  fprintf(file, "\n");
  long len2 = (matrix->ncols) - (matrix->nrows);
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < matrix->nrows; i++){
    fprintf(file, "%d ", matrix->dense_idx[i]);
  }
  fprintf(file, "\n");
}

static inline void display_fglm_colon_matrix(FILE *file, sp_matfglmcol_t *matrix){

  fprintf(file, "%u\n", matrix->charac);
  fprintf(file, "%u\n", matrix->ncols);
  fprintf(file, "%u\n", matrix->nrows);
  fprintf(file, "%u\n", matrix->nzero);

  long len1 = (matrix->ncols)*(matrix->nrows);
  for(long i = 0; i < len1; i++){
    fprintf(file, "%d ", matrix->dense_mat[i]);
  }
  fprintf(file, "\n");
  long len2 = (matrix->ncols) - (matrix->nrows);
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", matrix->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < matrix->nrows; i++){
    fprintf(file, "%d ", matrix->dense_idx[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < matrix->nzero; i++){
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
  for(long i = 0; i < nvars-1; i++){
    nmod_poly_init(param->coords[i], prime);
  }
  return param;
}

static inline void free_fglm_param(param_t *param){
  nmod_poly_clear(param->elim);
  nmod_poly_clear(param->denom);
  for(long i = 0; i < param->nvars-1; i++){
    nmod_poly_clear(param->coords[i]);
  }
  free(param->coords);
  free(param);
}

static inline fglm_bms_data_t *allocate_fglm_bms_data(long dim, mp_limb_t prime){

  fglm_bms_data_t * data_bms = (fglm_bms_data_t *)malloc(sizeof(fglm_bms_data_t));
  nmod_poly_init(data_bms->A, prime);

  nmod_poly_init(data_bms->B, prime);
  nmod_poly_init(data_bms->Z1, prime);

  nmod_poly_init2(data_bms->rZ1, prime, dim+1);

  nmod_poly_init(data_bms->Z2, prime);
  nmod_poly_init2(data_bms->rZ2, prime, dim+1);
  nmod_poly_init2(data_bms->V, prime, dim+1);

  nmod_poly_init2(data_bms->param, prime, dim+1);

  for(long i = 0; i < dim + 1; i++){
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
  count_leading_zeros(poly->mod.norm, prime);
  /* poly->mod.norm = flint_clz(prime); */

}

static inline void fglm_param_set_prime(param_t *param, mp_limb_t prime){
  param->charac = prime;

  nmod_poly_set_prime(param->elim, prime);
  nmod_poly_set_prime(param->denom, prime);
  for(long i = 0; i < param->nvars-1; i++){
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

  /* nmod_poly_factor_init(data_bms->sqf); */

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

