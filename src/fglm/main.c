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

/**

   First FGLM implementation

 **/

//#define DEBUGFGLM 1

//Slower with cache functions
//#define CACHE 2


#include "fglm_core.c"


#ifdef CACHE
static void inner_product_cache(CF_t *__restrict__ res, const CF_t *restrict a, const CF_t *restrict b,
                          const szmat_t sz, const mod_t prime,
                          const szmat_t bl, const szmat_t rst,
                          CF_l_t *__restrict__ vec_cache){
  CF_L_t c = 0;
  for(szmat_t i = 0; i < bl; ++i){
    const szmat_t pos = i * CACHE ;
    for(szmat_t j = 0; j < CACHE; ++j){
      szmat_t k = j+pos;
      vec_cache[j] = a[k]*b[k];
      //      c+=vec_cache[j];
    }
    for(int j = 0; j < CACHE; ++j){
      c+=vec_cache[j];
    }
  }
  const szmat_t pos = bl * CACHE;
  for(szmat_t i =0; i < rst; i++){
    const szmat_t k = i + pos;
    c += a[k] * b[k];
  }
  res[0] = (CF_t) (c % prime);
}
#else
static void inner_product(CF_t *__restrict__ res, const CF_t *restrict a, const CF_t *restrict b,
                          const szmat_t sz, const mod_t prime){
  CF_L_t c = 0;
  /* print_vec(stderr, a, sz); */
  /* print_vec(stderr, b, sz); */
  for(szmat_t i = 0; i < sz; ++i){
    c+= (CF_l_t)a[i]*b[i];
  }
  res[0] = (CF_t) (c % prime);
}
#endif

static inline void dense_mat_vec_mul(CF_t *__restrict__ res,
                                     const CF_t *__restrict__ a,
                                     const CF_t *__restrict__ b,
                                     const szmat_t nrows, const szmat_t ncols,
                                     const mod_t prime, CF_l_t *__restrict__ vec_cache){
  #ifndef CACHE
  //  int CACHE = 1;
  #else
  szmat_t bl = ncols / CACHE;
  szmat_t rst = ncols % CACHE;
  #endif
  for(szmat_t i = 0; i < nrows; i++){
    #ifdef CACHE
    inner_product_cache(res, a, b, ncols, prime, bl, rst, vec_cache);
    #else
    inner_product(res, a, b, ncols, prime);
    #endif
    a+=ncols;
    res++;
  }
}



static sp_matfglm_t *read_sp_matfglm(FILE *file, mod_t *prime_ptr){
  szmat_t ncols=0;
  szmat_t nrows=0;
  if(fscanf(file, "%u\n", prime_ptr)){
  }
  else{
    exit(1);
  }
  if(fscanf(file, "%u\n", &ncols)){
  }
  else{
    exit(1);
  }
  if(fscanf(file, "%u", &nrows)){
  }
  else{
    exit(1);
  }
  if(nrows>ncols){
    fprintf(stderr, "Bad data (nrows > ncols)\n");
    exit(1);
  }
  sp_matfglm_t *mat  ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  if(mat==NULL){
    return mat;
  }
  mat->charac = (*prime_ptr);
  szmat_t len = ncols * nrows ;
  //  cf_t * dense_mat = NULL;// = malloc(sizeof(cf_t) * len);
  if(!posix_memalign((void **)&mat->dense_mat, 32, sizeof(CF_t)*len)){
    for(szmat_t i = 0; i < len; ++i){
      if(fscanf(file, "%u", mat->dense_mat+i)){
      }
      else{
        exit(1);
      }
    }
  }
  else{
    fprintf(stderr, "Error when using posix_memalign(dense_mat)\n");
    return NULL;
  }
  if(fscanf(file, "\n")){}
  szmat_t len2 = ncols - nrows;

  if(!posix_memalign((void **)&mat->triv_idx, 32, sizeof(CF_t)*len2)){
    for(szmat_t i = 0; i < len2; ++i){
      if(fscanf(file, "%u", mat->triv_idx+i)){
      }
      else{
        exit(1);
      }
    }
  }
  else{
    fprintf(stderr, "Error when using posix_memalign(triv_idx)\n");
    return NULL;
  }
  if(!posix_memalign((void **)&mat->triv_pos, 32, sizeof(CF_t)*len2)){
    for(szmat_t i = 0; i < len2; ++i){
      if(fscanf(file, "%u", mat->triv_pos+i)==1){
      }
      else{
        exit(1);
      }
    }
  }
  else{
    fprintf(stderr, "Error when using posix_memalign(triv_pos)\n");
    return NULL;
  }
  if(fscanf(file, "\n")){}
  //  szmat_t *dense_idx = NULL;
  if(!posix_memalign((void **)&mat->dense_idx, 32, sizeof(CF_t)*nrows)){
    for(szmat_t i = 0; i < nrows; ++i){
      if(fscanf(file, "%u", mat->dense_idx+i)==1){
      }
      else{
        fprintf(stderr, "Didn't fill dense_idx entries\n");
        break;
        //        exit(1);
      }
    }
    fprintf(stderr, "Warning: dense_idx entries were not given in the file\n");
    szmat_t ct = 0;
    szmat_t cd = 0;
    for(szmat_t i = 0; i < ncols; ++i){
      if(mat->triv_idx[ct] != i){
        mat->dense_idx[cd] = i;
        cd++;
      }
      else{
        ct++;
      }
    }
  }
  else{
    fprintf(stderr, "Error when using posix_memalign(dense_idx)\n");
    return NULL;
  }
  if(posix_memalign((void **)&mat->dst, 32, sizeof(CF_t)*nrows)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < nrows; i++){
      mat->dst[i] = 0;
    }
  }

  mat->ncols = ncols;
  mat->nrows = nrows;

  //Ici on support que les entres de matrix->dst sont initialisees a 0
  for(long i = 0; i < mat->nrows; i++){
    for(long j = mat->ncols - 1; j >= 0; j--){
      if(mat->dense_mat[i*mat->ncols + j] == 0){
        mat->dst[i]++;
      }
      else{
        break;
      }
    }
  }


  #if DEBUGFGLM > 0
  fprintf(stderr, "mat->triv_pos = ");
  print_vec(stderr, mat->triv_pos, ncols-nrows);
  #endif
  return mat;
}






/***

Algorithme.

On se donne EMGCD(R0, R1)

qui retourne (R_i, R_{i+1}) et (V_{i}, V_{i+1}) tels que

R0 U_i + R1 V_i = R_i
R0 U_{i+1} + R1 V_{i+1} = R_{i+1}

avec deg(R_i) >= deg(R0) / 2
deg(R_{i+1}) < deg(R0) / 2

(donc comme Berlekamp-Massey quand R0 = x^(2*dimension))

On pose T la matrice de multiplication par x_n

On pose U=(u_0, ..., u_d) la suite calculee r^t T Vect(1)

[Vect(1) = (1, 0, ..., 0)]

ou d est la dimension du quotient

On veut resoudre

u_0, ..., u_{d-1}
.
.
u_{d-1}... u_{2d-2}

fois P = V

ou V= (v0, ..., v_{d-1})

Donc on a u_i = r^t T Vect(1)

et les v_i = r^t T^i (T_{x1} 1)

On a deja calcule r^t T^i et par ailleurs T_{x_1} 1 est le vecteur correspondant
a x1 dans la base monomiale

On se donne une fonction AD(d, U)

1- R0 = t^{2d-1}
   R1 = u0+u1 t + ... + u_{2d-2} t^{2d-2}

2- (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)

3- si deg(R_{i+1}) < d-1
 -> erreur (matrice singuliere)

4- si V_{i+1}(0) \neq 0

   Z1 = LC(R_{i+1})^{-1} x (V_{i+1})
   R1 = miroir(R1)
   (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)
   Z2 = LC(R_{i+1})^{-1} x (V_{i+1})

   renvoyer Z1, Z2

5. On est dans le cas V_{i+1} = 0

   R0 = t^{2d+1}
   R1 = 1+u0+...+u_{2d-2}t^{2d-1}+t^{2d}
   (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)

   Si deg(R_{i+1}) = d
     Z1 = LC(R_{i+1})^{-1} x (V_{i+1})
     R1=Miroir(R1)
     (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)
     Z2 = LC(R_{i+1})^{-1} x (V_{i+1})
     renvoyer Z1, Z2

  Sinon
  R0 = t^{2d+1}
  R1 = 1+u0+...+u_{2d-2}t^{2d-1}-t^{2d}
  (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)
  Z1 = LC(R_{i+1})^{-1} x (V_{i+1})
  R1=Miroir(R1)
  (R_i, R_{i+1}, V_i, V_{i+1}) = EMGCD(R_0, R_1)
  Z2 = LC(R_{i+1})^{-1} x (V_{i+1})
  renvoyer Z1, Z2


Et a la fin:

Solve (Z1,Z2,v)
//Z1 and Z2 have degree d-1
V = v_0+v_1*t+...+v_{d-1}*t^{d-1} // d is the dimension of the vector space.
A = Z1^r * V^r mod t^d //Z1^r is the mirror of Z1, same for V^r
B = Z2   * V^r mod t^d
Z3= (Z1*B^r - Z2^r*A^r)/Z1(0) mod t^d //Z1(0)!=0 guaranteed
return Z3^r

***/


#ifdef DEBUGFGLM
static inline void display_basic(){
  fprintf(stderr, "size of int %lu\n", sizeof(int));
  fprintf(stderr, "size of long int %lu\n", sizeof(long int));
  fprintf(stderr, "size of long long int %lu\n", sizeof(long long int));
  fprintf(stderr, "size of unsigned int %lu\n", sizeof(unsigned int));
  fprintf(stderr, "size of unsigned long int %lu\n", sizeof(unsigned long int));
  fprintf(stderr, "size of unsigned long long int %lu\n", sizeof(unsigned long long int));
  fprintf(stderr, "size of uint32_t %lu\n", sizeof(uint32_t));
  fprintf(stderr, "size of uint64_t %lu\n", sizeof(uint64_t));
  fprintf(stderr, "size of __uint128_t %lu\n\n", sizeof(__uint128_t));
}
#endif



int main(int argc, char *argv[]){

#ifdef DEBUGFGLM
  display_basic();
#endif

  fprintf(stderr, "%s\n", argv[1]);
  FILE *filename = fopen(argv[1],"r");
  FILE *output_fname = NULL;
  if(argc==3){
    output_fname = fopen(argv[2], "w");
  }
  else {
    output_fname = stdout;
  }
  mod_t prime = 65521;
  sp_matfglm_t *matrix ALIGNED32 = read_sp_matfglm(filename, &prime);

  fprintf(stderr, "Matrix read\n");
  long nvars = 4;
  if(nvars<=1){
    fprintf(stderr, "Problem with number of variables\n");
    exit(1);
  }
  uint64_t * linvars= calloc (nvars,sizeof(uint64_t));
  uint32_t * lineqs = calloc (nvars,sizeof(uint32_t));
  uint64_t * squvars= calloc (nvars,sizeof(uint64_t));
  param_t * param = nmod_fglm_compute(matrix, prime, nvars, nvars, linvars,
				      lineqs, squvars, 1);
  fprintf(stdout, "\n");
  display_fglm_param(output_fname, param);
  fprintf(stderr, "Done.\n");
  free_fglm_param(param);

  if(matrix!=NULL){
    free_sp_mat_fglm(matrix);
  }
  return 0;
}

