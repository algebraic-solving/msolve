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

#ifndef MSOLVE_DATA_H
#define MSOLVE_DATA_H

#define _GNU_SOURCE
#include "../neogb/data.h"
#include<flint/fmpz.h>
#include<flint/fmpq.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include "../neogb/libneogb.h"

#include <flint/nmod_poly.h>
#define MODP(a,b)                               \
	(a) % (b)


typedef struct
{
  mpz_t numer;
  long k;
  unsigned int isexact;
  int sign_left;
} interval;

typedef len_t nelts_t;
typedef int32_t nvars_t;
typedef int64_t bits_t;



typedef struct{
  nvars_t nvars;
  nvars_t elim;
  int32_t ngens;
  int32_t nterms;
  int32_t field_char;
  /* counts change of variable orders:
   * x1 <-> xn
   * x1 <-> xn-1
   * ...
   * for situation when staircase is not generic enough */
  int32_t change_var_order;
  /* base coefficient for linear form
   * sum(i^k*x[k]) k = 1, ..., nvars
   * It is zero if no linear form is active, otherwise != zero. */
  int32_t linear_form_base_coef;
  /* set to 1 if a linear form is chosen randomly */
  int32_t rand_linear;
  int32_t *random_linear_form;
  char **vnames;
  int32_t *lens;
  int32_t *exps;
  int32_t *cfs;    /* int32_t coeffs */
  mpz_t **mpz_cfs; /* mpz_t coeffs */
} data_gens_ff_t;

typedef uint32_t szmat_t;

typedef uint32_t CF_t;
typedef uint64_t CF_l_t;
typedef uint32_t mod_t;
typedef int32_t nvars_t;
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
  szmat_t nnfs; //number of required normal forms
  CF_t *dense_mat; // matrice nrows lignes et ncols colonnes (elements donnes par lignes)
  szmat_t *triv_idx; //tableau d'indices des lignes ne contenant que des 0 et un 1
  szmat_t *triv_pos; //position des 1
  szmat_t *dense_idx; //position des lignes non triviales (qui constituent donc
                      //dense_mat)
  szmat_t *dst; //pour la gestion des lignes "denses" mais avec un bloc de zero a la fin
  double totaldensity;
  double freepartdensity;
  double nonfreepartdensity;
} sp_matfglm_t;

#ifndef ALIGNED32
#define ALIGNED32 __attribute__((aligned(32)))
#endif
typedef struct{
  CF_t *vecinit ALIGNED32; /* init vector used in Wiedemann's implementation */
  CF_t *vecmult ALIGNED32; /* vector used to store the multiplication of the 
                            * dense part of the multiplication matrix with a 
                            * vector in Wiedemann */
  CF_t *vvec ALIGNED32; /* stores the result of matrix vector product in Wiedeman */
  CF_t *res ALIGNED32; /* array storing the term sequences needed after Wiedeman */
  mp_limb_t *pts;
} fglm_data_t;


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

typedef struct{
  mp_limb_t charac;
  nvars_t nvars;
  nmod_poly_t elim;
  nmod_poly_t denom;
  nmod_poly_t *coords;
  szmat_t degelimpol;
  szmat_t degsqfrelimpol;
} param_t;


typedef struct{
  mpz_t r0;
  mpz_t r1;
  mpz_t t0;
  mpz_t t1;
  mpz_t q;
  mpz_t tmp;
  mpz_t N;
  mpz_t D;
} rrec_data_struct_t;

typedef rrec_data_struct_t rrec_data_t[1];


typedef struct{
  deg_t length;
  deg_t alloc;
  mpz_t *coeffs;
} mpz_upoly_struct;

typedef mpz_upoly_struct mpz_upoly_t[1];

typedef struct{
  nvars_t nvars;
  deg_t nsols;
  len_t dquot;
  deg_t dim;
  mpz_upoly_t elim;
  mpz_upoly_t denom;
  mpz_upoly_t *coords;
  mpz_t *cfs;
} mpz_param_struct;
typedef mpz_param_struct mpz_param_t[1];

typedef struct{
  len_t nb;
  mpz_param_t *params;
} mpz_param_array_struct;
typedef mpz_param_array_struct mpz_param_array_t[1];

typedef struct{
  uint32_t ncols; /* dimension of quotient */
  uint32_t nrows; /* number of non trivial lines */
  mpz_t *dense_mat; /*array of nrows*ncols*2 mpz_t coefficients (num, den)*/
  mpz_t *denoms;
  uint32_t *triv_idx; /*array of indices of rows which are unit vectors*/
  uint32_t *triv_pos; /*position of '1' in unit vectors */
  uint32_t *dense_idx; /* array of rows which are NOT unit vectors */
  uint32_t *dst; /* blocks of 0's in non-trivial rows */
} mpq_matfglm_struct;
typedef mpq_matfglm_struct mpq_matfglm_t[1];

typedef struct{
  uint32_t ncols; /* dimension of quotient */
  uint32_t nrows; /* number of non trivial lines */
  mpz_t *dense_mat; /*array of nrows*ncols mpz_t coefficients*/
  mpz_t *denoms; /*denominators for rows which are not unit vectors*/
  uint32_t *triv_idx; /*array of indices of rows which are unit vectors*/
  uint32_t *triv_pos; /*position of '1' in unit vectors */
  uint32_t *dense_idx; /* array of rows which are NOT unit vectors */
  uint32_t *dst; /* blocks of 0's in non-trivial rows */
} mpz_matfglm_struct;
typedef mpz_matfglm_struct mpz_matfglm_t[1];

typedef struct{
  uint32_t ncols; /* dimension of quotient */
  uint32_t nrows; /* number of non trivial lines */
  mpz_t *dense_mat; /*array of nrows*ncols*2 mpz_t coefficients (num, den)*/
  uint32_t *triv_idx; /*array of indices of rows which are unit vectors*/
  uint32_t *triv_pos; /*position of '1' in unit vectors */
  uint32_t *dense_idx; /* array of rows which are NOT unit vectors */
  uint32_t *dst; /* blocks of 0's in non-trivial rows */
} crt_mpz_matfglm_struct;
typedef crt_mpz_matfglm_struct crt_mpz_matfglm_t[1];

typedef struct{
  uint32_t trace_idx;
  uint32_t det_idx;
  mpz_t trace_crt;
  mpz_t det_crt;
  mpz_t trace_num;
  mpz_t trace_den;
  mpz_t det_num;
  mpz_t det_den;
  mpz_t tmp;
  int16_t done_trace;
  int16_t done_det;
  int16_t check_trace;
  int16_t check_det;

  int lift_matrix;
  uint32_t nlins; /*number of linear forms to lift */
  int32_t nv; /* number of variables for linear forms */
  uint32_t lin_lifted; /* tells if linear forms are lifted */
  mpz_t *mpq_linear_forms;
  mpz_t *crt_linear_forms;
  mpz_t *mpz_linear_forms;

  uint32_t nrows; /* number of non trivial rows in matrix multiplication */
  uint32_t ncols; /* number of non trivial columns in matrix multiplication */
  szmat_t *triv_idx; 
  szmat_t *triv_pos; 
  szmat_t *dense_idx; 
  szmat_t *dst; 

  uint32_t nlifted; /* number of rows whose witness coefficients could be lifted */
  uint32_t w_checked; /* number of rows whose witness coefficients is checked */
  uint64_t *matmul_indices; /* indices per row of matrix multiplication used 
                              for lifting (witness coefficient); given index is the one of 
                              the dense matrix format */
  mpz_t *matmul_wcrt; /* crt for witness coefficients */
  mpz_t *matmul_wqq; /* stores reconstruction */
  int32_t mat_alloc; /* number of matrices that can be stored */
  int32_t *modular_matrices; /*array to modular matrices */
  int16_t *done_coeffs; /* indicates coefficients which have been lifted */
  int16_t *check_coeffs; /* indicates coefficients which have been checked */
  int32_t num_mat; /*number of modular matrices stored */

  mp_limb_t *primes;
  int32_t num_primes;
  mp_limb_t *residues;
  fmpz_comb_t comb;
  fmpz_comb_temp_t comb_temp;
  int32_t mat_lifted;

  mpz_t *mat_denoms;
  mpz_t *dense_mat;

} trace_det_fglm_mat_struct;
typedef trace_det_fglm_mat_struct trace_det_fglm_mat_t[1];

typedef struct{
  mpz_t val_up;
  mpz_t val_do;
  deg_t k_up;
  deg_t k_do;
  unsigned int isexact;
} coord_struct;
typedef coord_struct coord_t[1];

typedef struct{
  nvars_t nvars;
  coord_t *coords;
} real_point_struct;
typedef real_point_struct real_point_t[1];

typedef struct{
  char *in_file;
  char *bin_file;
  char *out_file;
  char *bin_out_file;
} files_gb;

/* data structure for tracing algorithms */
typedef struct{
  primes_t *lp; /* array of lucky primes, usually of size st->nthrds */
  bs_t *bs_qq; /* basis_qq */
  ht_t *bht; /* hash table */
  ht_t *tht; /* hash table to store the hashes of the multiples of the basis
                elements stored in the trace */
  bs_t **bs;
  int *bad_primes;
  trace_t **btrace;

  int32_t *num_gb; /* array storing lengths of computed GBs */
  int32_t **leadmons_ori; /* original leading monomials (from learning) */
  int32_t **leadmons_current; /* leading monomials (from tracing) */

  int32_t *mgb; /* array which stores one monomial */

  ht_t **blht;
  ht_t **btht;

  mpz_t mod_p;
  mpz_t prod_p;
} msolvetrace_data_struct;
typedef msolvetrace_data_struct mstrace_t[1];

typedef struct{
  int dim;
  int64_t dquot;
  int32_t ht_size; /* initial_hts */
  int32_t nr_threads;
  int32_t max_nr_pairs;
  int32_t elim_block_len;
  int32_t reset_ht;
  int32_t la_option;
  int32_t use_signatures;
  int32_t info_level;
  int32_t print_gb;
  int32_t truncate_lifting;
  int32_t pbm_file;
  files_gb *files;
} msolveflags_struct;
typedef msolveflags_struct msflags_t[1];
#endif
