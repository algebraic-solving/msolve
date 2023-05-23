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
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <ctype.h>
#include "../neogb/libneogb.h"

#define MODP(a,b)                               \
	(a) % (b)
typedef len_t nelts_t;
typedef int32_t nvars_t;

typedef struct{
  int32_t nvars;
  int32_t elim;
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
  int32_t length;
  int32_t alloc;
  mpz_t *coeffs;
} mpz_upoly_struct;

typedef mpz_upoly_struct mpz_upoly_t[1];

typedef struct{
  long nvars;
  long nsols;
  long dquot;
  int dim;
  mpz_upoly_t elim;
  mpz_upoly_t denom;
  mpz_upoly_t *coords;
  mpz_t *cfs;
} mpz_param_struct;
typedef mpz_param_struct mpz_param_t[1];

typedef struct{
  int32_t nb;
  mpz_param_t *params;
} mpz_param_array_struct;
typedef mpz_param_array_struct mpz_param_array_t[1];

typedef struct{
  uint32_t ncols; /* dimension of quotient */
  uint32_t nrows; /* number of non trivial lines */
  mpz_t *dense_mat; /*array of nrows*ncols*2 mpz_t coefficients (num, den)*/
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
  int done_trace;
  int done_det;
  int check_trace;
  int check_det;
} trace_det_fglm_mat_struct;
typedef trace_det_fglm_mat_struct trace_det_fglm_mat_t[1];

typedef struct{
  mpz_t val_up;
  mpz_t val_do;
  long k_up;
  long k_do;
  unsigned int isexact;
} coord_struct;
typedef coord_struct coord_t[1];

typedef struct{
  long nvars;
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
  long dquot;
  int32_t ht_size; /* initial_hts */
  int32_t nr_threads;
  int32_t max_nr_pairs;
  int32_t elim_block_len;
  int32_t reset_ht;
  int32_t la_option;
  int32_t use_signatures;
  int32_t info_level;
  int32_t print_gb;
  int32_t pbm_file;
  files_gb *files;
} msolveflags_struct;
typedef msolveflags_struct msflags_t[1];
#endif
