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

#include "msolve.h"
#include "duplicate.c"
#include "linear.c"
#include "lifting.c"
#include "lifting-gb.c"

#define LIFTMATRIX 0
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef MIN
#define MIN(x, y) ((x) > (y) ? (y) : (x))
#endif
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define ilog2_mpz(a) mpz_sizeinbase(a,2)

static void mpz_upoly_init(mpz_upoly_t poly, long alloc){
  mpz_t *tmp = NULL;
  if(alloc){
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < alloc; i++){
      mpz_init(tmp[i]);
      mpz_set_ui(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

static void mpz_upoly_init2(mpz_upoly_t poly, long alloc, long nbits){
  mpz_t *tmp = NULL;
  if(alloc){
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    //    tmp = (mpz_t *)(calloc(alloc, sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < alloc; i++){
      mpz_init2(tmp[i], nbits);
      mpz_set_ui(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

void mpz_upoly_clear(mpz_upoly_t pol){

  for(long i = 0; i < pol->alloc; i++){
    mpz_clear(pol->coeffs[i]);
  }
  free(pol->coeffs);
}

static inline void mpz_upoly_out_str(FILE *file, mpz_upoly_t pol){
  fprintf(file, "[");
  if(pol->length>0){
    fprintf(file, "%d, ", pol->length - 1); //degree
    fprintf(file, "[");
    for(long i = 0; i < pol->length - 1; i++){
      mpz_out_str(file, 10, pol->coeffs[i]);fprintf(file, ", ");
    }
    mpz_out_str(file, 10, pol->coeffs[pol->length - 1]);
    fprintf(file, "]");
  }
  else{
    fprintf(file, "-1, [0]");
  }
  fprintf(file, "]");
}

static inline void mpz_param_init(mpz_param_t param){
  param->nvars = 0;
  param->nsols = 0;
  mpz_upoly_init(param->elim, 0);
  mpz_upoly_init(param->denom, 0);
  param->coords=NULL;
  param->cfs=NULL;
}

static inline void mpz_param_clear(mpz_param_t param){
  mpz_upoly_clear(param->elim);
  mpz_upoly_clear(param->denom);
  if(param->coords != NULL){
    for(long i = 0; i < param->nvars - 1; i++){
      mpz_upoly_clear(param->coords[i]);
      mpz_clear(param->cfs[i]);
    }
  }
  free(param->coords);
  free(param->cfs);
  param->nvars  = 0;
  param->nsols  = 0;
}

static inline void mpz_param_out_str(FILE *file, const data_gens_ff_t *gens,
                                     const long dquot, mpz_param_t param,
                                     param_t *mod_param){
  fprintf(file, "[");
  fprintf(file, "%d, \n", gens->field_char); /* field charac */
  fprintf(file, "%ld, \n", param->nvars); //nvars
  fprintf(file, "%ld, \n", dquot); //dim quotient
  /* Print all variables:
   * 1. may include new variable from added linear form,
   * 2. variables may have a different order than in the.
   *    input system
   * Both of the above points are due to genericity handling. */
  fprintf(file, "[");
  for (int i = 0; i < param->nvars-1; ++i) {
      fprintf(file, "'%s', ", gens->vnames[i]);
  }
  fprintf(file, "'%s'],\n", gens->vnames[param->nvars-1]);

  fprintf(file, "[");
  if(gens->rand_linear){
    for(int i = 0; i < param->nvars - 1; i++){
      fprintf(file, "%d,", gens->random_linear_form[i]);
    }
    fprintf(file, "%d", gens->random_linear_form[param->nvars-1]);
  }
  else{
    if (gens->linear_form_base_coef > 0) {
      const int32_t bcf = gens->linear_form_base_coef;
      for (int i = 0; i < param->nvars-1; ++i) {
        fprintf(file, "%d,", (int32_t)(pow(i+1, bcf - 1)));
      }
      fprintf(file, "%d", (int32_t)(1));
    }
    else{
      for (int i = 0; i < param->nvars-1; ++i) {
        fprintf(file, "%d, ", (int32_t)(0));
      }
      fprintf(file, "%d", (int32_t)(1));
    }
  }
  fprintf(file, "],\n");
  fprintf(file, "[1,\n["); /*at the moment, a single param is returned */
  if(gens->field_char){
    display_nmod_poly(file, mod_param->elim);
  }
  else{
    mpz_upoly_out_str(file, param->elim); //elim. poly
  }
  fprintf(file, ",\n");
  if(gens->field_char){
    display_nmod_poly(file, mod_param->denom);
  }
  else{
  mpz_upoly_out_str(file, param->denom); //denom. poly
  }
  fprintf(file, ",\n");
  fprintf(file, "[\n");
  if(gens->field_char){/* positive characteristic */
    if(mod_param->coords != NULL){
      for(int i = 0; i < mod_param->nvars - 1; i++){
        fprintf(file, "[");
        if(gens->field_char){
          display_nmod_poly(file, mod_param->coords[i]);
        }
        if(i == mod_param->nvars - 2){
          fprintf(file, "]\n");
        }
        else{
          fprintf(file, "],\n");
        }
      }
    }
  }
  else{
    if(param->coords != NULL){
      for(int i = 0; i < param->nvars - 1; i++){
        fprintf(file, "[");
        mpz_upoly_out_str(file, param->coords[i]); //param. polys
        fprintf(file, ",\n");
        mpz_out_str(file, 10, param->cfs[i]);
        if(i==param->nvars-2){
          fprintf(file, "]\n");
        }
        else{
          fprintf(file, "],\n");
        }
      }
    }
  }
  /* fprintf(file, "]"); */
  fprintf(file, "]");
  fprintf(file, "]]");

}

static inline void mpz_param_out_str_maple(FILE *file,
        const data_gens_ff_t *gens,const long dquot,
                                           mpz_param_t param, param_t *mod_param){
  mpz_param_out_str(file, gens, dquot, param, mod_param);
  fprintf(file, "]");
}


static inline void display_fglm_crt_matrix(FILE *file,
                                           crt_mpz_matfglm_t mat){

  fprintf(file, "%u\n", mat->ncols);
  fprintf(file, "%u\n", mat->nrows);

  long len1 = (mat->ncols)*(mat->nrows);
  for(long i = 0; i < len1; i++){
    mpz_out_str(file, 10, mat->dense_mat[i]);
    fprintf(file, " ");
  }
  fprintf(file, "\n");
  long len2 = (mat->ncols) - (mat->nrows);
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", mat->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", mat->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < mat->nrows; i++){
    fprintf(file, "%d ", mat->dense_idx[i]);
  }
  fprintf(file, "\n");
}

static inline void display_fglm_mpq_matrix(FILE *file,
                                           mpq_matfglm_t mat){

  fprintf(file, "%u\n", mat->ncols);
  fprintf(file, "%u\n", mat->nrows);

  uint64_t nc = mat->ncols;

  fprintf(file, "[");
  for(long i = 0; i < mat->nrows; i++){
    long c = 2*i*mat->ncols;
    fprintf(file, "[");
    for(long j = 0; j < nc-1; j++){
      mpz_out_str(file, 10, mat->dense_mat[c+2*j]);
      fprintf(file, "/");
      mpz_out_str(file, 10, mat->dense_mat[c+2*j+1]);
      fprintf(file, ", ");
    }
    mpz_out_str(file, 10, mat->dense_mat[c+2*nc-2]);
    fprintf(file, "/");
    mpz_out_str(file, 10, mat->dense_mat[c+2*nc-1]);
    fprintf(file, "]\n");
  }
  fprintf(file, "]");
  fprintf(file, "\n");
  long len2 = (mat->ncols) - (mat->nrows);
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", mat->triv_idx[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < len2; i++){
    fprintf(file, "%d ", mat->triv_pos[i]);
  }
  fprintf(file, "\n");
  for(long i = 0; i < mat->nrows; i++){
    fprintf(file, "%d ", mat->dense_idx[i]);
  }
  fprintf(file, "\n");
}

static inline data_gens_ff_t *allocate_data_gens(){
  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
  gens->lens  = NULL;
  gens->exps  = NULL;
  gens->cfs   = NULL;
  gens->mpz_cfs = NULL;

  gens->elim = 0;
  return gens;
}

static inline void free_data_gens(data_gens_ff_t *gens){
  for(long i = 0; i < gens->nvars; i++){
    free(gens->vnames[i]);
  }
  free(gens->vnames);
  if (gens->field_char == 0) {
      for(long i = 0; i < 2*gens->nterms; i++){
          mpz_clear(*(gens->mpz_cfs[i]));
          free(gens->mpz_cfs[i]);
      }
  }
  free(gens->mpz_cfs);
  free(gens->lens);
  free(gens->cfs);
  free(gens->exps);
  free(gens->random_linear_form);
  free(gens);
}

static inline void display_term(FILE *file, int64_t ind, data_gens_ff_t *gens,
                                int32_t **blen, int32_t *bcf, int32_t **bexp){
  if(bcf[ind] != 0 && bcf[ind] != 1){
    fprintf(file, "%d", bcf[ind]);
    display_monomial(file, gens, ind, bexp);
  }
  else{
    if(bcf[ind] == 1){
      int32_t b = display_monomial_single(file, gens, ind, bexp);
      if(b==0){
        fprintf(file, "1");
      }
    }
  }
}

static inline void display_basis(FILE *file,
                                 int32_t *bld, int32_t **blen, int32_t **bexp,
                                 int32_t *bcf, data_gens_ff_t *gens){
  int32_t npos = 0;
  for(int32_t i = 0; i < bld[0]; i++){
    int32_t j, end = (*blen)[i];
    for(j = 0; j < end; j++){
      display_term(file, npos+j, gens, blen, bcf, bexp);
      if(j < (*blen)[i] - 1){
        fprintf(file, "+");
      }
    }
    //    display_term(file, npos+j, gens, blen, bcf, bexp);
    if(i < bld[0]-1){
      fprintf(file, ",\n");
    }
    else{
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
}

static inline void display_basis_maple(FILE *file,
                                 int32_t *bld, int32_t **blen, int32_t **bexp,
                                 int32_t *bcf, data_gens_ff_t *gens){
  int32_t npos = 0;
  fprintf(file, "[");
  for(int32_t i = 0; i < bld[0]; i++){
    int32_t j, end = (*blen)[i];
    for(j = 0; j < end; j++){
      display_term(file, npos+j, gens, blen, bcf, bexp);
      if(j < (*blen)[i] - 1){
        fprintf(file, "+");
      }
    }
    //    display_term(file, npos+j, gens, blen, bcf, bexp);
    if(i < bld[0]-1){
      fprintf(file, ",\n");
    }
    else{
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
  fprintf(file, "]:\n");
}

static inline void display_lead_monomials_from_gb(FILE *file,
                                                  int32_t *bld, int32_t **blen, int32_t **bexp,
                                                  int32_t *bcf, data_gens_ff_t *gens){
  int32_t npos = 0;
  for(int32_t i = 0; i < bld[0]; i++){
    display_term(file, npos, gens, blen, bcf, bexp);
    if(i < bld[0]-1){
      fprintf(file, ",\n");
    }
    else{
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
}

static inline void display_monomials_from_array(FILE *file, long length,
                                         int32_t *bexp,
                                         const int nv, char **vnames){
  fprintf(file, "[");
  for(long i = 0; i < length-1; i++){
    display_monomial_full(file, nv, vnames, i, bexp);
    fprintf(file, ", ");
  }
  display_monomial_full(file, nv, vnames, length - 1, bexp);
  fprintf(file, "]");
}

static inline void display_monomials_from_array_maple(FILE *file, long length,
                                                      int32_t *bexp, const int nv,
                                                      char **vnames){
  fprintf(file, "[");
  for(long i = 0; i < length-1; i++){
    display_monomial_full(file, nv, vnames, i, bexp);
    fprintf(file, ", ");
  }
  display_monomial_full(file, nv, vnames, length - 1, bexp);
  fprintf(file, "]:");
}

/* Used if staircase is not generic enough, before
 * we add a linear form with a new variable.
 *
 * Returns 1 if change of variable order is done,
 * 0 if all cyclic changes have already been checked
 * and we should go on to add a linear form with a
 * new variable. */
static int undo_variable_order_change(
        data_gens_ff_t *gens
        )
{
    int32_t i, j;
    int32_t len, tmp;
    char *tmp_char      = NULL;
    const int32_t cvo   = gens->change_var_order;
    const int32_t nvars = gens->nvars;
    const int32_t ngens = gens->ngens;

    if (gens->linear_form_base_coef > 0) {
        return 0;
    }
    if (cvo > -1) {
        /* undo last variable change */
        tmp_char = gens->vnames[nvars-1];
        gens->vnames[nvars-1] = gens->vnames[cvo];
        gens->vnames[cvo]  = tmp_char;
        len = 0;
        tmp = 0;
        for (i = 0; i < ngens; ++i) {
            for (j = 0; j < gens->lens[i]; ++j) {
                tmp = gens->exps[len + j * nvars + nvars - 1];
                gens->exps[len + j * nvars + nvars - 1] =
                    gens->exps[len + j * nvars + cvo];
                gens->exps[len + j * nvars + cvo] = tmp;
            }
            len +=  gens->lens[i]*nvars;
        }
    }
    /* printf("undone:  ");
     * for(int i = 0; i < nvars-1; i++){
     *     fprintf(stdout, "%s, ", gens->vnames[i]);
     * }
     * fprintf(stdout, "%s\n", gens->vnames[nvars-1]); */
    /* all cyclic changes already done, stop here, try to add
     * a linear form with additional variable afterwards */
    gens->change_var_order++;
    if (gens->change_var_order == nvars-1) {
        return 0;
    } else {
        return 1;
    }
}
static int change_variable_order_in_input_system(
        data_gens_ff_t *gens,
        int32_t info_level
        )
{
    int32_t i, j;
    int32_t len, tmp;
    char *tmp_char      = NULL;
    const int32_t cvo   = gens->change_var_order;
    const int32_t nvars = gens->nvars;

    if (undo_variable_order_change(gens) == 0) {
        return 0;
    }
    /* do the current variable change */
    tmp_char = gens->vnames[nvars-1];
    gens->vnames[nvars-1]  = gens->vnames[cvo+1];
    gens->vnames[cvo+1]  = tmp_char;
    len = 0;
    tmp = 0;
    for (i = 0; i < gens->ngens; ++i) {
        for (j = 0; j < gens->lens[i]; ++j) {
            tmp = gens->exps[len + j * nvars + nvars-1];
            gens->exps[len + j * nvars + nvars-1] =
                gens->exps[len + j * nvars + cvo + 1];
            gens->exps[len + j * nvars + cvo + 1] = tmp;
        }
        len +=  gens->lens[i]*nvars;
    }
    if (info_level > 0) {
        printf("\nChanging variable order for possibly more generic staircase:\n");
        for(int i = 0; i < nvars-1; i++){
            fprintf(stdout, "%s, ", gens->vnames[i]);
        }
        fprintf(stdout, "%s\n", gens->vnames[nvars-1]);
    }
    return 1;
}

/* Used if staircase is not generic enough, after
 * we tried to change variable order.
 *
 * Returns 1 if addition of linear form is done,
 * 0 if an error occurs. */
static int add_linear_form_to_input_system(
        data_gens_ff_t *gens,
        int32_t info_level
        )
{
    int64_t i, j;
    int32_t k;
    int32_t nvars_old, nvars_new;
    int64_t len_old = 0, len_new;
    if (gens->linear_form_base_coef == 0) {
        nvars_old = gens->nvars;
        nvars_new = nvars_old + 1;
        for (i = 0; i < gens->ngens; ++i) {
            len_old +=  gens->lens[i];
        }
        len_new = len_old + nvars_new;
    } else {
        nvars_old = gens->nvars - 1;
        nvars_new = nvars_old + 1;
        for (i = 0; i < gens->ngens-1; ++i) {
            len_old +=  gens->lens[i];
        }
        len_new = len_old + gens->lens[gens->ngens-1];
    }
    /* for the first run we have to reinitialize gens data with the newly
     * added variable, then we add the linear form at the end. */
    if (gens->linear_form_base_coef == 0) {
        char *extra_var = (char *)malloc(2*sizeof(char));
        strcpy(extra_var, "A");
        gens->nvars++;
        gens->ngens++;
        gens->lens  = realloc(gens->lens,
                (unsigned long)gens->ngens * sizeof(int32_t));
        gens->lens[gens->ngens-1] = (int32_t)nvars_new;
        /* add new dummy variable (only need for correct freeing of vnames at
         * the end */
        gens->vnames  = realloc(
                gens->vnames, (unsigned long)gens->nvars * sizeof(char *));
        gens->vnames[gens->nvars-1] = (char *)malloc(2 * sizeof(char));
        gens->vnames[gens->nvars-1] = extra_var;
        int32_t *old_exps = gens->exps;
        gens->exps  =
            (int32_t *)calloc(
                    (unsigned long)len_new * nvars_new, sizeof(int32_t));
        i = 0;
        j = 0;
        while (i < (len_old * nvars_old)) {
            memcpy(gens->exps+j, old_exps+i,
                    (unsigned long)nvars_old * sizeof(int32_t));
            i +=  nvars_old;
            j +=  nvars_new;
        }
        free(old_exps);
        /* add linear form exponents */
        while (j < (len_new * nvars_new)) {
            gens->exps[j] = 1;
            j += nvars_new+1;
        }
        /* allocate memory for coefficients of linear form */
        if (gens->field_char > 0) {
            gens->cfs = realloc(gens->cfs,
                    (unsigned long)len_new * sizeof(int32_t));
        } else {
            gens->mpz_cfs = realloc(gens->mpz_cfs,
                    (unsigned long)2 * len_new * (sizeof(mpz_t *)));
            for (i = 2*len_old; i < 2*len_new; i += 2) {
                gens->mpz_cfs[i]    = (mpz_t *)malloc(sizeof(mpz_t));
                mpz_init(*(gens->mpz_cfs[i]));
                gens->mpz_cfs[i+1]  = (mpz_t *)malloc(sizeof(mpz_t));
                mpz_init(*(gens->mpz_cfs[i+1]));
                /* set denominators 1 for all coefficients */
                mpz_set_ui(*(gens->mpz_cfs[i+1]), 1);
            }
        }
    }
    gens->linear_form_base_coef++;
    const int32_t bcf = gens->linear_form_base_coef;
    k = 1;
    if (info_level > 0) {
        printf("\nAdding a linear form with an extra variable ");
        printf("(lowest w.r.t. monomial order)\n");
        printf("[coefficients of linear form are k^%d for k looping over variable index 1...n]\n", bcf - 1);
    }
    if (gens->field_char > 0) {
        for (i = len_old; i < len_new - 1; ++i) {
            gens->cfs[i]  = ((int32_t)(pow(k, bcf - 1)) % gens->field_char);
            k++;
        }
        gens->cfs[len_new - 1]  = 1;
        k++;
    } else {
        for (i = 2*len_old; i < 2*len_new; i += 2) {
            mpz_set_ui(*(gens->mpz_cfs[i]), (int32_t)(pow(k, bcf - 1)));
            k++;
        }
        mpz_set_si(*(gens->mpz_cfs[2*(len_new - 1)]), 1);
    }
    return 1;
}

static int add_random_linear_form_to_input_system(
        data_gens_ff_t *gens,
        int32_t info_level
        )
{

    int64_t i, j;
    int32_t nvars_old, nvars_new;
    int64_t len_old = 0, len_new;

    if (gens->linear_form_base_coef == 0) {
        nvars_old = gens->nvars;
        nvars_new = nvars_old + 1;
        for (i = 0; i < gens->ngens; ++i) {
            len_old +=  gens->lens[i];
        }
        len_new = len_old + nvars_new;
    } else {
        nvars_old = gens->nvars - 1;
        nvars_new = nvars_old + 1;
        for (i = 0; i < gens->ngens-1; ++i) {
            len_old +=  gens->lens[i];
        }
        len_new = len_old + gens->lens[gens->ngens-1];
    }

    /* for the first run we have to reinitialize gens data with the newly
     * added variable, then we add the linear form at the end. */
    if (gens->linear_form_base_coef == 0) {
        char *extra_var = (char *)malloc(2*sizeof(char));
        strcpy(extra_var, "A");
        gens->nvars++;
        gens->ngens++;
        gens->lens  = realloc(gens->lens,
                (unsigned long)gens->ngens * sizeof(int32_t));
        gens->lens[gens->ngens-1] = (int32_t)nvars_new;
        /* add new dummy variable (only need for correct freeing of vnames at
         * the end */
        gens->vnames  = realloc(
                gens->vnames, (unsigned long)gens->nvars * sizeof(char *));
        gens->vnames[gens->nvars-1] = (char *)malloc(2 * sizeof(char));
        gens->vnames[gens->nvars-1] = extra_var;
        int32_t *old_exps = gens->exps;
        gens->exps  =
            (int32_t *)calloc(
                    (unsigned long)len_new * nvars_new, sizeof(int32_t));
        i = 0;
        j = 0;
        while (i < (len_old * nvars_old)) {
            memcpy(gens->exps+j, old_exps+i,
                    (unsigned long)nvars_old * sizeof(int32_t));
            i +=  nvars_old;
            j +=  nvars_new;
        }
        free(old_exps);
        /* add linear form exponents */
        while (j < (len_new * nvars_new)) {
            gens->exps[j] = 1;
            j += nvars_new+1;
        }
        /* allocate memory for coefficients of linear form */
        if (gens->field_char > 0) {
            gens->cfs = realloc(gens->cfs,
                    (unsigned long)len_new * sizeof(int32_t));
        } else {
            gens->mpz_cfs = realloc(gens->mpz_cfs,
                    (unsigned long)2 * len_new * (sizeof(mpz_t *)));
            for (i = 2*len_old; i < 2*len_new; i += 2) {
                gens->mpz_cfs[i]    = (mpz_t *)malloc(sizeof(mpz_t));
                mpz_init(*(gens->mpz_cfs[i]));
                gens->mpz_cfs[i+1]  = (mpz_t *)malloc(sizeof(mpz_t));
                mpz_init(*(gens->mpz_cfs[i+1]));
                /* set denominators 1 for all coefficients */
                mpz_set_ui(*(gens->mpz_cfs[i+1]), 1);
            }
        }
    }
    gens->linear_form_base_coef++;
    /* const int32_t bcf = gens->linear_form_base_coef; */
    if (info_level > 0) {
        printf("\nAdding a linear form with an extra variable ");
        printf("(lowest w.r.t. monomial order)\n");
        printf("[coefficients of linear form are randomly chosen]\n");
    }
    srand(time(0));
    /* gens->random_linear_form = malloc(sizeof(int32_t)*(nvars_new)); */
    gens->random_linear_form = realloc(gens->random_linear_form, sizeof(int32_t)*(nvars_new));

    if (gens->field_char > 0) {
      int j = 0;
      for (i = len_old; i < len_new; ++i) {
        gens->random_linear_form[j] = ((int8_t)(rand()) % gens->field_char);

        while(gens->random_linear_form[j] == 0){
            gens->random_linear_form[j] = ((int8_t)(rand()) % gens->field_char);
       }
        gens->cfs[i]  = gens->random_linear_form[j];
        j++;
      }
    }
    else {
      int j = 0;
      for (i = 2*len_old; i < 2*len_new; i += 2) {
        gens->random_linear_form[j] = ((int8_t)(rand()));

        while(gens->random_linear_form[j] == 0){
            gens->random_linear_form[j] = ((int8_t)(rand()));
        }

        mpz_set_si(*(gens->mpz_cfs[i]), gens->random_linear_form[j]);
        mpz_set_ui(*(gens->mpz_cfs[i+1]), 1);

        j++;
      }
    }
    gens->rand_linear = 1;
    return 1;
}

static inline void nmod_param_out_str(FILE *file, const long dquot, param_t *param){
  fprintf(file, "[");
  fprintf(file, "0, \n"); //charac
  fprintf(file, "%ld, \n", param->nvars); //nvars
  fprintf(file, "%ld, \n", dquot); //dim quotient
  display_nmod_poly(file, param->elim); //elim. poly
  fprintf(file, ",\n");
  display_nmod_poly(file, param->denom); //elim. poly
  fprintf(file, ",\n");
  fprintf(file, "[\n");
  for(int i = 0; i < param->nvars - 1; i++){
    fprintf(file, "[");
    display_nmod_poly(file, param->coords[i]); //elim. poly
    fprintf(file, ",\n");
    fprintf(file, "[1]");
    if(i==param->nvars-2){
      fprintf(file, "]\n");
    }
    else{
      fprintf(file, "],\n");
    }
  }
  fprintf(file, "]");
  fprintf(file, "]");
}

static inline void nmod_param_out_str_maple(FILE *file, const long dquot, param_t *param){
  nmod_param_out_str(file, dquot, param);
  fprintf(file, ":\n");
}


static inline void print_msolve_message(FILE * file, int n){
  if(n==0){
    fprintf(file, "The ideal has positive dimension.\n");
  }
  if(n==1){
    fprintf(file, "The staircase of the Grobner basis is not generic enough.\n");
    fprintf(file, "\nYou may add to your system a random linear combination of \nthe variables plus a new variable which you choose to be \nthe smallest for the monomial ordering.\n\n");
  }
}

static inline void initialize_mpz_param(mpz_param_t param, param_t *bparam){

  param->nvars = bparam->nvars;
  param->nsols = bparam->elim->length - 1;

  mpz_upoly_init2(param->elim, bparam->elim->alloc, 2*32*(bparam->elim->length));
  mpz_upoly_init2(param->denom, bparam->elim->alloc - 1, 2*32*(bparam->elim->length));
  param->elim->length = bparam->elim->length;

  param->coords = (mpz_upoly_t *)malloc(sizeof(mpz_upoly_t)*(param->nvars - 1));
  if(param->coords != NULL){
    for(long i = 0; i < param->nvars - 1; i++){

      mpz_upoly_init(param->coords[i], MAX(1,bparam->elim->alloc - 1));
      /* param->coords[i]->length = bparam->coords[i]->length; */
      param->coords[i]->length = bparam->elim->length - 1;

    }
  }
  else{
    fprintf(stderr, "Error when initializing parametrization\n");
    exit(1);
  }
  param->cfs = (mpz_t *)malloc(sizeof(mpz_t) * (param->nvars - 1));
  if(param->cfs != NULL){
    for(int i = 0; i < param->nvars - 1; i++){
      mpz_init(param->cfs[i]);
    }
  }
  else{
    fprintf(stderr, "Error when allocating cfs\n");
    exit(1);
  }
}

static inline void reduce_generators(mpz_t **tmp, long nterms,
                                     int32_t * cfs, int32_t prime){
  for(int32_t i = 0; i < 2*nterms; i += 2){
    cfs[i/2] = (int32_t)mpz_fdiv_ui(tmp[i][0], prime);
  }
}

static inline void check_and_set_linear_poly_non_hashed(long *nlins_ptr,
                                                        uint64_t *linvars,
                                                        uint32_t** lineqs_ptr,
                                                        int32_t *bld,
                                                        int32_t *bexp_lm,
                                                        int32_t **blen,
                                                        int32_t **bexp,
                                                        int32_t *bcf,
                                                        long nvars){
  long nlins = 0;
  uint32_t *coefpos = calloc(nvars, sizeof(uint32_t));
  uint32_t pos = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1.
  */
  for(long i = 0; i < bld[0]; i++){
    long deg = 0;
    for(int j = 0; j < nvars; j++){
      deg+=bexp_lm[i*nvars+j];
    }
    if(deg == 1){
      nlins++;
      for(int k = 0; k < nvars; k++){
        if(bexp_lm[i*nvars+k] == 1){
          linvars[k] = i + 1;
          coefpos[k] = pos;
        }
      }
    }
    pos += (*blen)[i];
  }
  *nlins_ptr = nlins;

  //On recupere les coefficients des formes lineaires
  uint32_t* lineqs = calloc(nlins*(nvars + 1), sizeof(uint32_t));

  int cnt = 0;

  for(int i = 0; i < nvars; i++){
    if(linvars[i] != 0){

      long len = (*blen)[linvars[i] - 1];

      if(len==nvars+1){
        for(long j = 0; j<len; j++){
          uint32_t coef = (uint32_t)((bcf)[j + coefpos[i]]);
          lineqs[cnt*(nvars+1)+j] = coef;
        }
      }
      else{
        //        hm_t *dt = bs->hm[bi] + OFFSET;
        for(long j = 0; j<len; j++){
          uint32_t coef = (uint32_t)((bcf)[j + coefpos[i]]);
          //bs->cf_32[bs->hm[bi][COEFFS]][j];

          int32_t *exp = (*bexp) + (j + coefpos[i]) * nvars; //bht->ev[dt[j]];
          int isvar = 0;
          for(int k = 0; k < nvars; k++){
            if(exp[k]==1){
              lineqs[cnt*(nvars+1)+k] = coef;
              isvar=1;
            }
          }
          if(isvar==0){
            lineqs[cnt*(nvars+1)+nvars] = coef;
          }
        }
        cnt++;
      }
    }
  }
  free(coefpos);
  lineqs_ptr[0] = lineqs;

}

static inline void
check_and_set_vars_squared_in_monomial_basis(uint64_t *squvars,
					     int32_t *lmb,
					     long dquot,
					     long nvars){
  /*
    la i-ieme entree de squvars est a 0 si le monome vars[i]^2 n'est
    pas dans la base monomiale.
    Sinon, on met l'indice du monome. De toute facon, en 0, on aurait
    le monome 1.
  */
  for (long i = 0; i < dquot; i++) {
    long deg = 0;
    for (long j = 0; j < nvars; j++){
      deg+=lmb[i*nvars+j];
    }
    if (deg == 2) {
      for (long j = 0; j < nvars-1; j++) {
	if (lmb[i*nvars+j] == 2) {
	  /* vars[j]^2 is the monomial */
	  /* fprintf (stderr, */
	  /* 	   "vars[%ld]^2 is the %ldth monomial of the monomial basis\n",j,i); */
	  squvars[j]= i;
	  break;
	}
      }
    }
  }
}

static inline void normalize_nmod_param(param_t *nmod_param){
  if(nmod_param != NULL){
    uint32_t prime = nmod_param->charac;
    uint32_t inv;
    ulong inv2;
    n_gcdinv(&inv2, nmod_param->elim->length - 1, prime);
    inv = inv2 % prime;
    nmod_poly_fit_length(nmod_param->denom, nmod_param->elim->length - 1);
    nmod_param->denom->length = nmod_param->elim->length - 1;

    for(long i = 1; i < nmod_param->elim->length; i++){
      nmod_param->denom->coeffs[i-1] = (i * nmod_param->elim->coeffs[i]) % prime;
    }

    for(long i = 0; i < nmod_param->elim->length-1; i++){
      nmod_param->denom->coeffs[i] = (inv * nmod_param->denom->coeffs[i]) % prime;
    }

    for(int j = 0; j < nmod_param->nvars - 1; j++){
      nmod_poly_mul(nmod_param->coords[j], nmod_param->coords[j],
                    nmod_param->denom);

      nmod_poly_rem(nmod_param->coords[j], nmod_param->coords[j],
                    nmod_param->elim);
    }
  }
}

/*
On copie nmod_param dans mpz_param
On suppose que les allocations ont ete faites via initialize_mpz_param avant.
 */
static inline void set_mpz_param_nmod(mpz_param_t mpz_param, param_t *nmod_param){
  mpz_param->elim->length = nmod_param->elim->length;
  for(long i = 0; i <= mpz_param->nsols; i++){
    mpz_set_ui(mpz_param->elim->coeffs[i], (nmod_param)->elim->coeffs[i]);
  }
  mpz_param->elim->length = nmod_param->elim->length;

  for(long i = 0; i < nmod_param->denom->length; i++){
    mpz_set_ui(mpz_param->denom->coeffs[i], nmod_param->denom->coeffs[i]);
  }
  mpz_param->denom->length = nmod_param->denom->length;

  for(int j = 0; j < mpz_param->nvars - 1; j++){

    for(long i = 0 ; i < nmod_param->coords[j]->length; i++){
      mpz_set_ui(mpz_param->coords[j]->coeffs[i],
                 nmod_param->coords[j]->coeffs[i]);
    }
    for(long i = nmod_param->coords[j]->length; i < nmod_param->elim->length - 1; i++){
      mpz_set_ui(mpz_param->coords[j]->coeffs[i], 0);
    }
    mpz_param->coords[j]->length = nmod_param->elim->length - 1;
  }
}

static inline void crt_lift_mpz_upoly(mpz_upoly_t pol, nmod_poly_t nmod_pol,
                                      mpz_t modulus, int32_t prime,
                                      mpz_t prod, mpz_t tmp,
                                      int nthrds){
  long i;

/* #pragma omp parallel for num_threads(nthrds)    \ */
/*   private(i) schedule(static) */
  for(i = 0; i < pol->length; i++){
    mpz_CRT_ui(pol->coeffs[i], pol->coeffs[i], modulus,
               nmod_pol->coeffs[i], prime, prod, tmp, 1);
  }


}


/* assumes that all degrees are the same */
static inline void crt_lift_mpz_param(mpz_param_t mpz_param, param_t *nmod_param,
                                      mpz_t modulus, mpz_t prod_crt,
                                      const int32_t prime, mpz_t tmp, const int nthrds){

  /*assumes prod_crt = modulus * prime */
  crt_lift_mpz_upoly(mpz_param->elim, nmod_param->elim, modulus, prime,
                     prod_crt, tmp, nthrds);
  for(long i = 0; i < mpz_param->nvars - 1; i++){

    crt_lift_mpz_upoly(mpz_param->coords[i], nmod_param->coords[i],
                       modulus, prime, prod_crt, tmp, nthrds);

  }

}





/**

   la sortie est recons / denominator

 **/
#define RR 1

static inline int rational_reconstruction_mpz_ptr(mpz_t *recons,
                                                  mpz_t denominator,
                                                  mpz_t *pol,
                                                  long len,
                                                  mpz_t modulus,
                                                  long *maxrec,
                                                  mpq_t *coef,
                                                  mpz_t rnum,
                                                  mpz_t rden,
                                                  mpz_t *tmp_num,
                                                  mpz_t *tmp_den,
                                                  mpz_t lcm,
                                                  mpz_t guessed_num,
                                                  mpz_t guessed_den,
                                                  rrec_data_t rdata,
                                                  int info_level){


  if(ratrecon(rnum, rden, pol[*maxrec], modulus, rdata) == 0){
    return 0;
  }
  mpz_set(tmp_num[*maxrec], rnum);
  mpz_set(tmp_den[*maxrec], rden);

  for(long i = *maxrec + 1; i < len; i++){
    int b = ratrecon(rnum, rden, pol[i], modulus, rdata);
    if(b == 0){
      *maxrec = i - 1;
      return b;
    }

    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }
  for(long i = 0; i < *maxrec; i++){
    int b = ratrecon(rnum, rden, pol[i], modulus, rdata);

    if(b == 0){
      if(info_level){
        fprintf(stderr, "[*]");
      }
      *maxrec = MAX(i-1, 0);
      return b;
    }

    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }

  mpz_set(lcm, tmp_den[0]);
  for(long i = 1; i < len; i++){
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_divexact(tmp_den[i], lcm, tmp_den[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_mul(tmp_num[i], tmp_num[i], tmp_den[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_set(recons[i], tmp_num[i]);
  }
  mpz_set(denominator, lcm);

  return 1;

}

static inline int rational_reconstruction_mpz_ptr_with_denom(mpz_t *recons,
                                                  mpz_t denominator,
                                                  mpz_t *pol,
                                                  long len,
                                                  mpz_t modulus,
                                                  long *maxrec,
                                                  mpq_t *coef,
                                                  mpz_t rnum,
                                                  mpz_t rden,
                                                  mpz_t *tmp_num,
                                                  mpz_t *tmp_den,
                                                  mpz_t lcm,
                                                  mpz_t guessed_num,
                                                  mpz_t guessed_den,
                                                  rrec_data_t rdata,
                                                  int info_level){

  mpz_set(guessed_num, pol[*maxrec]);

  if(ratreconwden(rnum, rden, guessed_num, modulus, guessed_den, rdata) == 0){
    return 0;
  }

  mpz_set(tmp_num[*maxrec], rnum);
  mpz_set(tmp_den[*maxrec], rden);

  for(long i = *maxrec + 1; i < len; i++){
    mpz_set(guessed_num, pol[i]);
    int b = ratreconwden(rnum, rden, guessed_num, modulus, guessed_den, rdata);

    if(b == 0){
      *maxrec = MAX(0, i - 1);
      return b;
    }
    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }

  mpz_set(lcm, tmp_den[*maxrec]);
  for(long i = *maxrec + 1; i < len; i++){
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }

  mpz_t newlcm;
  mpz_init(newlcm);
  mpz_set(newlcm, lcm);
  mpz_mul(newlcm, newlcm, guessed_den);
  mpz_fdiv_q(rdata->D, rdata->D, lcm);
  mpz_mul(rdata->N, rdata->N, lcm);

  for(long i = *maxrec-1; i >=0; i--){
    mpz_set(guessed_num, pol[i]);
    int b = ratreconwden(tmp_num[i], tmp_den[i],
                         guessed_num, modulus, newlcm, rdata);

      if(b == 0){
        *maxrec = MAX(i + 1, 0);
        mpz_clear(newlcm);
        return b;
      }

      mpz_divexact(rden, newlcm, guessed_den);
      mpz_mul(tmp_den[i], tmp_den[i], rden);

      mpz_lcm(newlcm, newlcm, rden);

  }

  mpz_set(lcm, tmp_den[0]);
  for(long i = 1; i < len; i++){
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }

  for(long i = 0; i < len; i++){
    mpz_divexact(tmp_den[i], lcm, tmp_den[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_mul(tmp_num[i], tmp_num[i], tmp_den[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_set(recons[i], tmp_num[i]);
  }
  mpz_set(denominator, lcm);
  mpz_clear(newlcm);
  return 1;

}


/**

   la sortie est recons / denominator

 **/

static inline int rational_reconstruction_upoly(mpz_upoly_t recons,
                                                mpz_t denominator,
                                                mpz_upoly_t pol,
                                                long len,
                                                mpz_t modulus,
                                                long *maxrec,
                                                mpq_t *coef,
                                                mpz_t rnum,
                                                mpz_t rden,
                                                mpz_upoly_t tmp_num,
                                                mpz_upoly_t tmp_den,
                                                mpz_t lcm,
                                                mpz_t guessed_num,
                                                mpz_t guessed_den,
                                                rrec_data_t rdata,
                                                int info_level){

  return rational_reconstruction_mpz_ptr(recons->coeffs, denominator,
                                         pol->coeffs, len,
                                         modulus, maxrec, coef,
                                         rnum, rden,
                                         tmp_num->coeffs,
                                         tmp_den->coeffs,
                                         lcm, guessed_num, guessed_den,
                                         rdata,
                                         info_level);
}

static inline int rational_reconstruction_upoly_with_denom(mpz_upoly_t recons,
                                                mpz_t denominator,
                                                mpz_upoly_t pol,
                                                long len,
                                                mpz_t modulus,
                                                long *maxrec,
                                                mpq_t *coef,
                                                mpz_t rnum,
                                                mpz_t rden,
                                                mpz_upoly_t tmp_num,
                                                mpz_upoly_t tmp_den,
                                                mpz_t lcm,
                                                mpz_t guessed_num,
                                                mpz_t guessed_den,
                                                rrec_data_t rdata,
                                                int info_level){

  return rational_reconstruction_mpz_ptr_with_denom(recons->coeffs, denominator,
                                         pol->coeffs, len,
                                         modulus, maxrec, coef,
                                         rnum, rden,
                                         tmp_num->coeffs,
                                         tmp_den->coeffs,
                                         lcm, guessed_num, guessed_den,
                                         rdata,
                                         info_level);
}




/* /\** */

/* renvoie 0 si on n'a pas reussi a reconstruire */

/* numer contient le resultat */

/* **\/ */

/* static inline int rational_reconstruction(mpz_param_t mpz_param, */
/*                                           param_t *nmod_param, */
/*                                           mpz_upoly_t numer, */
/*                                           mpz_upoly_t denom, */
/*                                           mpz_t *modulus, int32_t prime, */
/*                                           mpq_t *coef, */
/*                                           mpz_t *guessed_num, mpz_t *guessed_den, */
/*                                           long *maxrec, */
/*                                           const int info_level){ */

/*   crt_lift_mpz_param(mpz_param, nmod_param, *modulus, prime, 1); */
/*   mpz_mul_ui(*modulus, *modulus, prime); */

/*   mpz_fdiv_q_2exp(*guessed_num, *modulus, 1); */
/*   mpz_sqrt(*guessed_num, *guessed_num); */
/*   mpz_set(*guessed_den, *guessed_num); */

/*   long deg = mpz_param->nsols; */
/*   if(mpq_reconstruct_mpz(coef, mpz_param->elim->coeffs[*maxrec], *modulus) == 0){ */
/*     return 0; */
/*   } */

/*   mpz_set(numer->coeffs[*maxrec], mpq_numref(*coef)); */
/*   mpz_set(denom->coeffs[*maxrec], mpq_denref(*coef)); */

/*   for(long i = *maxrec + 1; i <= deg; i++){ */
/*     int b = mpq_reconstruct_mpz_with_denom(coef, mpz_param->elim->coeffs[i], */
/*                                            *modulus, */
/*                                            *guessed_num, *guessed_den); */
/*     if(b == 0){ */
/*       *maxrec = i - 1; */
/*       return b; */
/*     } */
/*     mpz_set(numer->coeffs[i], mpq_numref(*coef)); */
/*     mpz_set(denom->coeffs[i], mpq_denref(*coef)); */
/*   } */
/*   for(long i = 0; i < *maxrec; i++){ */
/*     int b = mpq_reconstruct_mpz_with_denom(coef, mpz_param->elim->coeffs[i], */
/*                                            *modulus, */
/*                                            *guessed_num, *guessed_den); */

/*     if(b == 0){ */
/*       if (info_level) { */
/*         fprintf(stderr, "[%ld]->[%ld]", *maxrec, i); */
/*       } */
/*       *maxrec = MAX(i-1, 0); */
/*       return b; */
/*     } */
/*     mpz_set(numer->coeffs[i], mpq_numref(*coef)); */
/*     mpz_set(denom->coeffs[i], mpq_denref(*coef)); */
/*   } */

/*   mpz_t lcm; */
/*   mpz_init_set(lcm, denom->coeffs[0]); */
/*   for(long i = 1; i <= deg; i++){ */
/*     mpz_lcm(lcm, lcm, denom->coeffs[i]); */
/*   } */
/*   for(long i = 0; i <= deg; i++){ */
/*     mpz_divexact(denom->coeffs[i], lcm, denom->coeffs[i]); */
/*   } */
/*   for(long i = 0; i <= deg; i++){ */
/*     mpz_mul(numer->coeffs[i], numer->coeffs[i], denom->coeffs[i]); */
/*   } */
/*   mpz_clear(lcm); */
/*   return 1; */
/* } */


/**

returns 0 if rational reconstruction failed

 **/

static inline int new_rational_reconstruction(mpz_param_t mpz_param,
                                              mpz_param_t tmp_mpz_param,
                                              param_t *nmod_param,
                                              mpq_matfglm_t mpq_mat,
                                              crt_mpz_matfglm_t crt_mat,
                                              mpz_matfglm_t mpz_mat,
                                              trace_det_fglm_mat_t trace_det,
                                              sp_matfglm_t *mat,
                                              mpz_upoly_t numer,
                                              mpz_upoly_t denom,
                                              mpz_t *modulus, mpz_t prod_crt,
                                              int32_t prime,
                                              mpq_t *coef,
                                              mpz_t rnum, mpz_t rden,
                                              rrec_data_t recdata,
                                              mpz_t *guessed_num,
                                              mpz_t *guessed_den,
                                              long *maxrec,
                                              long *matrec,
                                              int *is_lifted,
                                              int *mat_lifted,
                                              int doit,
                                              int nthrds,
                                              const int info_level){

  mpz_mul_ui(prod_crt, *modulus, prime);
  crt_lift_mpz_param(tmp_mpz_param, nmod_param, *modulus, prod_crt,
                     prime, trace_det->tmp, nthrds);

  uint32_t trace_mod = nmod_param->elim->coeffs[trace_det->trace_idx];
  uint32_t det_mod = nmod_param->elim->coeffs[trace_det->det_idx];
  int b = 1;
  if(trace_det->done_trace == 0){
    b = 0;
    if(check_trace(trace_det, trace_mod, prime)){
      if(trace_det->check_trace == 0){
        trace_det->check_trace = 1;
      }
      else{
        trace_det->done_trace = 1;
        *maxrec = trace_det->det_idx;
      }
    }
  }
  if(trace_det->done_det == 0){
    b = 0;
    if(check_det(trace_det, det_mod, prime)){
      if(trace_det->check_det == 0){
        trace_det->check_det = 1;
      }
      else{
        trace_det->done_det = 1;
        *maxrec = trace_det->det_idx;
      }
    }
  }
  crt_lift_trace_det(trace_det,
                     trace_mod,
                     det_mod,
                     *modulus, prod_crt, prime);

#if LIFTMATRIX == 1
  if(*matrec < crt_mat->nrows*crt_mat->ncols){
    crt_lift_mat(crt_mat, mat, *modulus, prod_crt, prime, trace_det->tmp, nthrds);
  }
#endif
  mpz_mul_ui(*modulus, *modulus, prime);

  if(doit==0){
    return 0;
  }

  mpz_fdiv_q_2exp(*guessed_num, *modulus, 1);
  mpz_sqrt(*guessed_num, *guessed_num);
  mpz_set(*guessed_den, *guessed_num);
  mpz_set(recdata->N, *guessed_num);
  mpz_set(recdata->D, *guessed_den);
#if LIFTMATRIX == 1
  long cnt = 0;
  if(*matrec < crt_mat->nrows*crt_mat->ncols && *mat_lifted == 0){
    cnt = rat_recon_dense_rows(mpq_mat, crt_mat, mpz_mat, *modulus, recdata,
                               rnum, rden, matrec);
  }
  if(cnt == crt_mat->nrows*crt_mat->ncols || *mat_lifted){
    *mat_lifted = 1;
  }
#else
  *mat_lifted = 1;
#endif
  rat_recon_trace_det(trace_det, recdata,*modulus, rnum, rden);
  if(b && trace_det->done_trace == 1 && trace_det->done_det == 1){

    mpz_t denominator;
    mpz_init(denominator);
    mpz_t lcm;
    mpz_init(lcm);
    int b = 0;

    if(is_lifted[0]==0){

      if(trace_det->done_trace == 1){
        mpz_set(*guessed_den, trace_det->trace_den);
        if(trace_det->done_det == 1){
          mpz_lcm(*guessed_den, *guessed_den, trace_det->det_den);
        }
      }
      else{
        mpz_set(*guessed_den, trace_det->det_den);
        if(trace_det->done_trace == 1){
          mpz_lcm(*guessed_den, *guessed_den, trace_det->trace_den);
        }
      }

      mpz_fdiv_q(recdata->D, recdata->D, *guessed_den);
      mpz_mul(recdata->D, recdata->D, recdata->D);
      mpz_fdiv_q(recdata->N, *modulus, recdata->D);
      mpz_fdiv_q_2exp(recdata->N, recdata->N, 1);

      mpz_root(recdata->D, *modulus, 3);
      mpz_fdiv_q(recdata->N, *modulus, recdata->D);
      mpz_fdiv_q_2exp(recdata->N, recdata->N, 1);

      b = rational_reconstruction_upoly_with_denom(mpz_param->elim,
                                        denominator,
                                        tmp_mpz_param->elim,
                                        nmod_param->elim->length,
                                        *modulus,
                                        maxrec,
                                        coef,
                                        rnum,
                                        rden,
                                        numer,
                                        denom,
                                        lcm,
                                        *guessed_num,
                                        *guessed_den,
                                        recdata,
                                        info_level);
      if(b == 0){
        mpz_root(recdata->D, *modulus, 16);
        mpz_fdiv_q(recdata->N, *modulus, recdata->D);
        mpz_fdiv_q_2exp(recdata->N, recdata->N, 1);

        b = rational_reconstruction_upoly_with_denom(mpz_param->elim,
                                                     denominator,
                                                     tmp_mpz_param->elim,
                                                     nmod_param->elim->length,
                                                     *modulus,
                                                     maxrec,
                                                     coef,
                                                     rnum,
                                                     rden,
                                                     numer,
                                                     denom,
                                                     lcm,
                                                     *guessed_num,
                                                     *guessed_den,
                                                     recdata,
                                                     info_level);
        if(b==0){
          is_lifted[0] = 0;
          mpz_clear(denominator);
          mpz_clear(lcm);
          return b;
        }
      }
      is_lifted[0] = 1;
      if(info_level){
        fprintf(stderr, "[0]");
      }
    }

    long nsols = mpz_param->nsols;
    mpz_t lc;
    mpz_init(lc);
    mpz_set(lc, mpz_param->elim->coeffs[nsols]);
    mpz_mul_ui(lc, lc, nsols);
    mpq_t c;
    mpq_init(c);
    mpz_set_ui(mpq_numref(c), 1);
    mpz_set_ui(mpq_denref(c), 1);
    int nc = mpz_param->nvars - 1;

    mpz_set(*guessed_den, lc);

    mpz_fdiv_q_2exp(*guessed_num, *modulus, 1);
    mpz_sqrt(recdata->D, *guessed_num);
    mpz_set(recdata->N, recdata->D);

    for(int i = 0; i < nc; i++){
      *maxrec = MIN(MAX(0, trace_det->det_idx-1), MAX(0,nmod_param->coords[i]->length - 1));

      if(is_lifted[0]>0 && is_lifted[i+1]==0){

        b = rational_reconstruction_upoly_with_denom(mpz_param->coords[i],
                                                     denominator,
                                                     tmp_mpz_param->coords[i],
                                                     nmod_param->coords[i]->length,
                                                     *modulus,
                                                     maxrec,
                                                     coef,
                                                     rnum,
                                                     rden,
                                                     numer,
                                                     denom,
                                                     lcm,
                                                     *guessed_num,
                                                     *guessed_den,
                                                     recdata,
                                                     info_level);

        if(b == 0){
          mpz_set_ui(recdata->D, 1);
          mpz_mul_2exp(recdata->D, recdata->D, nc);
          mpz_fdiv_q_2exp(recdata->N, *modulus, 1);
          mpz_fdiv_q(recdata->N, recdata->N, recdata->D);

          b = rational_reconstruction_upoly_with_denom(mpz_param->coords[i],
                                        denominator,
                                        tmp_mpz_param->coords[i],
                                        nmod_param->coords[i]->length,
                                        *modulus,
                                        maxrec,
                                        coef,
                                        rnum,
                                        rden,
                                        numer,
                                        denom,
                                        lcm,
                                        *guessed_num,
                                        *guessed_den,
                                        recdata,
                                        info_level);

          if(b == 0){
            mpz_fdiv_q_2exp(recdata->N, *modulus, 1);
            mpz_root(recdata->D, recdata->N, 16);
            mpz_fdiv_q(recdata->N, recdata->N, recdata->D);

            b = rational_reconstruction_upoly_with_denom(mpz_param->coords[i],
                                                   denominator,
                                                   tmp_mpz_param->coords[i],
                                                   nmod_param->coords[i]->length,
                                                   *modulus,
                                                   maxrec,
                                                   coef,
                                                   rnum,
                                                   rden,
                                                   numer,
                                                   denom,
                                                   lcm,
                                                   *guessed_num,
                                                   *guessed_den,
                                                   recdata,
                                                   info_level);
            if(b==0){

              mpz_clear(denominator);
              mpz_clear(lcm);
              mpz_clear(lc);
              mpq_clear(c);

              is_lifted[i+1] = 0;
              return b;
            }
          }
        }
      }
      else{
        /* indicates that there is no need to set up below the data as they were already computed */
        /* also denominator = 0 by now */
        b = 0;
      }
      if(info_level && b && is_lifted[i+1] == 0){
        fprintf(stderr, "[%d]", i+1);
      }
      if(b){
        is_lifted[i+1] = 1;

        mpz_set(mpq_denref(c), denominator);
        mpq_canonicalize(c);
        mpz_set(mpz_param->cfs[i], mpq_denref(c));

        for(long j = 0; j < mpz_param->coords[i]->length; j++){
          mpz_mul(mpz_param->coords[i]->coeffs[j], mpz_param->coords[i]->coeffs[j],
                  mpq_numref(c));
        }
      }
    }

    mpz_clear(denominator);
    mpz_clear(lcm);
    mpz_clear(lc);
    mpq_clear(c);

    return 1;

  }
  return 0;
}





/**
   on verifie que mpz_pol / lc(mpz_pol) mod prime = nm_pol
   renvoie 1 si il faut faire le modular check.
**/
static inline int check_unit_mpz_nmod_poly(const long len,
                                           const mpz_upoly_t mpz_pol,
                                           const nmod_poly_t nm_pol,
                                           const int32_t prime){
  uint32_t lc = mpz_fdiv_ui(mpz_pol->coeffs[len - 1], prime);
  lc = mod_p_inverse_32(lc, prime);
  for(long i = 0; i < len; i++){
    uint64_t c = mpz_fdiv_ui(mpz_pol->coeffs[i], prime);
    c *= (uint64_t)lc;
    c = c % prime;
    if(c != nm_pol->coeffs[i]){
      return i + 1;
    }
  }
  return 0;
}

static inline int check_param_nmod_poly(const long len,
                                        const mpz_upoly_t mpz_pol,
                                        const mpz_t den,
                                        const mpz_t lcelim,
                                        const long nbsol,
                                        const nmod_poly_t nm_pol,
                                        const int32_t prime){
  if(len == 0){
    return 0;
  }
  mpz_t binv, bprime;
  mpz_init(binv);
  mpz_init_set_ui(bprime, prime);

  mpz_mul(binv, lcelim, den);
  mpz_mul_ui(binv, binv, nbsol);
  mpz_invert(binv, binv, bprime);

  uint32_t inv = mpz_mod_ui(binv, binv, prime);

  for(long i = 0; i < len; i++){

    mpz_mul_ui(binv, mpz_pol->coeffs[i], inv);
    if(mpz_mod_ui(binv, binv, prime) != (nm_pol->coeffs[i] % prime)){
      return 1;
    }
  }
  mpz_clear(binv);
  mpz_clear(bprime);
  return 0;
}


/**
   renvoie 1 si il faut faire le modular check.
**/

static inline int check_param_modular(const mpz_param_t mp_param,
                                      const param_t *bparam,
                                      const int32_t prime,
                                      int *is_lifted,
                                      trace_det_fglm_mat_t trace_det,
                                      const int info_level){

  long len = mp_param->nsols + 1;

  int c = check_unit_mpz_nmod_poly(len,
                                    mp_param->elim,
                                    bparam->elim,
                                    prime);
  if(c){
    if(info_level){
      fprintf(stderr, "<0,%d>", c);
    }
    is_lifted[0] = 0;
    for(int i = 0; i<mp_param->nvars-1; i++){
      is_lifted[i+1] = 0;
    }
    trace_det->done_trace = 0;
    trace_det->check_trace = 0;
    trace_det->done_det = 0;
    trace_det->check_det = 0;
    return 1;
  }

  for(int i = 0; i <mp_param->nvars-1; i++ ){
    len = mp_param->coords[0]->length;

    if(check_param_nmod_poly(bparam->coords[i]->length,
                             mp_param->coords[i],
                             mp_param->cfs[i],
                             mp_param->elim->coeffs[mp_param->elim->length - 1],
                             mp_param->elim->length - 1,
                             bparam->coords[i],
                             prime)){
      is_lifted[i+1] = 0;
      if(info_level){
        fprintf(stderr, "<%d>", i+1);
      }
      return 1;
    }
  }
  return 0;
}


static inline void get_leading_ideal_information(
        int32_t *num_gb,
        int32_t **lead_mons,
        const int32_t pos,
        const bs_t * const bs
        )
{
    lead_mons[pos] = get_lm_from_bs(bs, bs->ht);
    num_gb[pos]    = bs->lml;
}

static inline void print_groebner_basis(
        files_gb *files,
        const data_gens_ff_t * const gens,
        const bs_t * const bs,
        md_t *md,
        const int32_t fc
        )
{
    if(md->print_gb){
        int32_t gfc = md->gfc;
        md->gfc     = fc;
        print_ff_basis_data(files->out_file, "a", bs, bs->ht, md,
                gens, md->print_gb);
        md->gfc     = gfc;
    }
}

static int32_t check_for_single_element_groebner_basis(
        int *dim,
        long *dquot_ori,
        const bs_t * const bs,
        int32_t **leadmons,
        const int32_t pos,
        const md_t * const md
        )
{
    int32_t empty_solution_set = 1;
    int32_t i;

    if (bs->lml == 1) {
        if (md->info_level > 0) {
            fprintf(stdout, "Grobner basis has a single element\n");
        }
        for (i = 0; i < bs->ht->nv; i++) {
            if (leadmons[pos][i] != 0) {
                empty_solution_set = 0;
                break;
            }
        }
        if (empty_solution_set == 1) {
            *dquot_ori = 0;
            *dim = 0;
            if (md->info_level > 0) {
                fprintf(stdout, "No solution\n");
            }
        }
    } else {
        empty_solution_set = 0;
    }

    return empty_solution_set;
}

static int32_t *initial_modular_step(
        sp_matfglm_t **bmatrix,
        int32_t **bdiv_xn,
        int32_t **blen_gb_xn,
        int32_t **bstart_cf_gb_xn,

        long *nlins_ptr,
        uint64_t *linvars,
        uint32_t** lineqs_ptr,
        uint64_t *squvars,

        fglm_data_t **bdata_fglm,
        fglm_bms_data_t **bdata_bms,

        int32_t *num_gb,
        int32_t **leadmons,
        uint64_t *bsz,
        param_t **bparam,
        bs_t *gbg,
        md_t *md,
        const int32_t fc,
        int print_gb,
        int *dim,
        long *dquot_ori,
        data_gens_ff_t *gens,
        files_gb *files,
        int *success)
{
    double rt = realtime();

    md->print_gb = print_gb;

    int32_t error              = 0;
    int32_t empty_solution_set = 1;
    bs_t *bs = core_gba(gbg, md, &error, fc);

    print_tracer_statistics(stdout, rt, md);

    get_leading_ideal_information(num_gb, leadmons, 0, bs);

    print_groebner_basis(files, gens, bs, md, fc);

    empty_solution_set = check_for_single_element_groebner_basis(dim, dquot_ori,
            bs, leadmons, 0, md);

    if (empty_solution_set == 1) {
        return NULL;
    }

    check_and_set_linear_poly(nlins_ptr, linvars, lineqs_ptr, bs->ht,
            leadmons[0], bs);
    if (has_dimension_zero(bs->lml, bs->ht->nv, leadmons[0])) {
        long dquot = 0;
        int32_t *lmb = monomial_basis(bs->lml, bs->ht->nv, leadmons[0], &dquot);

        if(md->info_level){
            fprintf(stderr, "Dimension of quotient: %ld\n", dquot);
        }
        if(print_gb==0){
            *bmatrix = build_matrixn_from_bs_trace(bdiv_xn,
                    blen_gb_xn,
                    bstart_cf_gb_xn,
                    lmb, dquot, bs, bs->ht,
                    leadmons[0], bs->ht->nv,
                    fc,
                    md->info_level);
            if(*bmatrix == NULL){
                *success = 0;
                *dim = 0;
                *dquot_ori = dquot;
                return NULL;
            }

            *bsz = bs->ht->nv - (*nlins_ptr); //nlins ;

            check_and_set_vars_squared_in_monomial_basis(squvars, lmb,
                    dquot, gens->nvars);
            *bparam = nmod_fglm_compute_trace_data(*bmatrix, fc, bs->ht->nv,
                    *bsz, *nlins_ptr, linvars, lineqs_ptr[0], squvars,
                    md->info_level, bdata_fglm, bdata_bms, success, md);
        }
        free_basis_without_hash_table(&(bs));
        *dim = 0;
        *dquot_ori = dquot;
        return lmb;
    }
    else{
        *dim  = 1;
        *dquot_ori = -1;
        free_basis_without_hash_table(&(bs));
        return NULL;
    }
}



/**

   On appelle F4 en generant la trace du calcul.
   On recupere les monomes dominants.

   Si l'ideal est de dimension > 0, on renvoie un pointeur nul.

   Si la dimension est <= 0, on construit la matrice de multiplication si
   l'escalier est suffisamment generique.

   Dans ce cas, on appelle alors un algorithme de changement d'ordre et on renvoie
   un tableau contenant la base monomiale.

**/

#if 0
static int32_t * modular_trace_learning(sp_matfglm_t **bmatrix,
                                        int32_t **bdiv_xn,
                                        int32_t **blen_gb_xn,
                                        int32_t **bstart_cf_gb_xn,

                                        long *nlins_ptr,
                                        uint64_t *linvars,
                                        uint32_t** lineqs_ptr,
                                        uint64_t *squvars,

                                        fglm_data_t **bdata_fglm,
                                        fglm_bms_data_t **bdata_bms,

                                        int32_t *num_gb,
                                        int32_t **leadmons,
                                        uint64_t *bsz,
                                        param_t **bparam,
                                        trace_t *trace,
                                        bs_t *bs_qq,
                                        md_t *st,
                                        const int32_t fc,
                                        int info_level,
                                        int print_gb,
                                        int *dim,
                                        long *dquot_ori,
                                        data_gens_ff_t *gens,
                                        files_gb *files,
                                        int *success)
{
    double rt = realtime();

    bs_t *bs = NULL;
    int32_t err = 0;
    bs = core_gba(bs_qq, st, &err, fc);
#if 0
    if(gens->field_char){
      if (err) {
        printf("Problem with F4, stopped computation.\n");
        exit(1);
      }
      /* free_shared_hash_data(bht); */
    }
    else{
      if(st->laopt > 40){
        bs = modular_f4(bs_qq, bs->ht, st, fc);
      }
      else{
        bs = gba_trace_learning_phase(trace, st->tr->ht, bs_qq, bs->ht, st, fc);
      }
    }
#endif
    print_tracer_statistics(stdout, rt, st);

    /* Leading monomials from Grobner basis */
    ht_t *bht = bs->ht;
    int32_t *bexp_lm = get_lm_from_bs(bs, bht);
    leadmons[0] = bexp_lm;
    num_gb[0] = bs->lml;

    if(bs->lml == 1){
        if(info_level){
            fprintf(stderr, "Grobner basis has a single element\n");
        }
        int is_empty = 1;
        for(int i = 0; i < bht->nv; i++){
            if(bexp_lm[i]!=0){
                is_empty = 0;
            }
        }
        if(is_empty){
            *dquot_ori = 0;
            *dim = 0;
            if(info_level){
              fprintf(stderr, "No solution\n");
            }
            print_ff_basis_data(
                                files->out_file, "a", bs, bht, st, gens, print_gb);
            return NULL;
        }
    }
    if(print_gb){
      if(st->gfc == 0){
        /* to fix display inconsistency when gens->fc = 0 */
        st->gfc = fc;

        print_ff_basis_data(
                            files->out_file, "a", bs, bht, st, gens, print_gb);
        st->gfc = 0;

      }
      else{
        print_ff_basis_data(
                            files->out_file, "a", bs, bht, st, gens, print_gb);
      }
    }
    check_and_set_linear_poly(nlins_ptr, linvars, lineqs_ptr, bht, bexp_lm, bs);

    if(has_dimension_zero(bs->lml, bht->nv, bexp_lm)){

      long dquot = 0;
      int32_t *lmb = monomial_basis(bs->lml, bht->nv, bexp_lm, &dquot);

        if(info_level){
            fprintf(stderr, "Dimension of quotient: %ld\n", dquot);
        }
        if(print_gb==0){
        *bmatrix = build_matrixn_from_bs_trace(bdiv_xn,
                                               blen_gb_xn,
                                               bstart_cf_gb_xn,
                                               lmb, dquot, bs, bht,
                                               bexp_lm, bht->nv,
                                               fc,
                                               info_level);
        if(*bmatrix == NULL){
            *success = 0;
            *dim = 0;
            *dquot_ori = dquot;
            return NULL;
        }

        *bsz = bht->nv - (*nlins_ptr); //nlins ;

        check_and_set_vars_squared_in_monomial_basis(squvars, lmb,
                dquot, gens->nvars);
        *bparam = nmod_fglm_compute_trace_data(*bmatrix, fc, bht->nv,
                                               *bsz,
                                               *nlins_ptr,
                                               linvars,
                                               lineqs_ptr[0],
                                               squvars,
                                               info_level,
                                               bdata_fglm, bdata_bms,
                                               success,
					       st);
        }
        free_basis(&(bs));
        *dim = 0;
        *dquot_ori = dquot;
        return lmb;
    }
    else{
        *dim  = 1;
        *dquot_ori = -1;
        free_basis(&(bs));
        return NULL;
    }
}
#endif

#if 0
static int32_t * modular_probabilistic_first(sp_matfglm_t **bmatrix,
                                             int32_t **bdiv_xn,
                                             int32_t **blen_gb_xn,
                                             int32_t **bstart_cf_gb_xn,

                                             long *nlins_ptr,
                                             uint64_t *linvars,
                                             uint32_t** lineqs_ptr,
                                             uint64_t * squvars,

                                             fglm_data_t **bdata_fglm,
                                             fglm_bms_data_t **bdata_bms,

                                             int32_t *num_gb,
                                             int32_t **leadmons,
                                             uint64_t *bsz,
                                             param_t **bparam,
                                             const bs_t *bs_qq,
                                             ht_t *bht,
                                             md_t *st,
                                             const int32_t fc,
                                             int info_level,
                                             int32_t print_gb,
                                             int *dim,
                                             long *dquot_ori,
                                             data_gens_ff_t *gens,
                                             files_gb *files,
                                             int *success)
{
    double ca0;
    ca0 = realtime();
    bs_t *bs = modular_f4(bs_qq, bht, st, fc);
    int32_t *bexp_lm = get_lm_from_bs(bs, bht);
    leadmons[0] = bexp_lm;
    num_gb[0] = bs->lml;

    if(bs->lml == 1){
      if(info_level){
        fprintf(stderr, "Grobner basis has a single element\n");
      }
      int is_empty = 1;
      for(int i = 0; i < bht->nv; i++){
        if(bexp_lm[i]!=0){
          is_empty = 0;
        }
      }
      if(is_empty){
        *dquot_ori = 0;
        *dim = 0;
        return NULL;
      }
    }

    print_ff_basis_data(
            files->out_file, "a", bs, bht, st, gens, print_gb);

    check_and_set_linear_poly(nlins_ptr, linvars, lineqs_ptr, bht, bexp_lm, bs);

    if(has_dimension_zero(bs->lml, bht->nv, bexp_lm)){
        long dquot = 0;
        int32_t *lmb = monomial_basis(bs->lml, bht->nv, bexp_lm, &dquot);
        if(info_level){
          fprintf(stderr, "Dimension of quotient: %ld\n", dquot);
        }

        *bmatrix = build_matrixn_from_bs_trace(bdiv_xn,
                                               blen_gb_xn,
                                               bstart_cf_gb_xn,
                                               lmb, dquot, bs, bht,
                                               bexp_lm, bht->nv,
                                               fc,
                                               info_level);

        if(*bmatrix == NULL){
          *success = 0;
          *dim = 0;
          *dquot_ori = dquot;
          return NULL;
        }

        *bsz = bht->nv - (*nlins_ptr); //nlins ;
        check_and_set_vars_squared_in_monomial_basis(squvars, lmb,
                dquot, gens->nvars);
        *bparam = nmod_fglm_compute_trace_data(*bmatrix, fc, bht->nv, *bsz,
                *nlins_ptr, linvars, lineqs_ptr[0],
                squvars, info_level,
		bdata_fglm, bdata_bms, success,st);
        *dim = 0;
        *dquot_ori = dquot;

        if (info_level) {
            fprintf(stderr, "First probabilistic F4 + FGLM time (elapsed): %.2f sec\n",
                    realtime() - ca0);
        }
        return lmb;
    }
    else{
        *dim  = 1;
        return NULL;
    }
}
#endif


#if 0
static void modular_probabilistic_apply(sp_matfglm_t **bmatrix,
                               int32_t **div_xn,
                               int32_t **len_gb_xn,
                               int32_t **start_cf_gb_xn,

                               long nlins,
                               uint64_t *linvars,
                               uint32_t *lineqs,
                               uint64_t *squvars,

                               fglm_data_t **bdata_fglm,
                               fglm_bms_data_t **bdata_bms,

                               int32_t *num_gb,
                               int32_t **leadmons_ori,
                               int32_t **leadmons_current,

                                        uint64_t bsz,
                               param_t **nmod_params,
                               const bs_t *bs_qq,
                               ht_t *bht,
                               md_t *st,
                               const int32_t fc,
                               int info_level,
                               bs_t **bs,
                               int32_t *lmb_ori,
                               int32_t dquot_ori,
                               primes_t *lp,
                               data_gens_ff_t *gens,
                               double *stf4,
                               const long nbsols,
                               int *bad_primes){

st->info_level = 0;
/* tracing phase */
len_t i;
double ca0;


memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));
/* #pragma omp parallel for num_threads(st->nthrds)  \ */
/*   private(i) schedule(dynamic) */
for (i = 0; i < st->nprimes; ++i){
  ca0 = realtime();
  bs[i] = modular_f4(bs_qq, bht, st, lp->p[i]);
  *stf4 = realtime()-ca0;

  if (bs[i] == NULL) {
      return;
  }
  int32_t lml = bs[i]->lml;
  if (st->nev > 0) {
      int32_t j = 0;
      for (len_t k = 0; k < bs[i]->lml; ++k) {
          if ((*bht)->ev[bs[i]->hm[bs[i]->lmps[k]][OFFSET]][0] == 0) {
              bs[i]->lm[j]   = bs[i]->lm[k];
              bs[i]->lmps[j] = bs[i]->lmps[k];
              ++j;
          }
      }
      lml = j;
  }
  if(lml != num_gb[i]){
    /* nmod_params[i] = NULL; */
    bad_primes[i] = 1;
    return;
  }
  get_lm_from_bs_trace(bs[i], bht, leadmons_current[i]);

  if(equal_staircase(leadmons_current[i], leadmons_ori[i],
                      num_gb[i], num_gb[i], bht->nv)){
    set_linear_poly(nlins, lineqs, linvars, bht, leadmons_current[i], bs[i]);

    build_matrixn_from_bs_trace_application(bmatrix[i],
                                            div_xn[i],
                                            len_gb_xn[i],
                                            start_cf_gb_xn[i],
                                            lmb_ori, dquot_ori, bs[i], bht,
                                            leadmons_ori[i], bht->nv,
                                            lp->p[i]);

    if(nmod_fglm_compute_apply_trace_data(bmatrix[i], lp->p[i],
                                          nmod_params[i],
                                          bht->nv,
                                          bsz,
                                          nlins, linvars, lineqs, squvars,
                                          bdata_fglm[i],
                                          bdata_bms[i],
                                          nbsols,
                                          info_level,
					  st)){
      bad_primes[i] = 1;
    }
    free_basis(&(bs[i]));
  }
  else{
      bad_primes[i] = 1;
  }
 }
}
#endif

static void secondary_modular_steps(sp_matfglm_t **bmatrix,
                                   int32_t **div_xn,
                                   int32_t **len_gb_xn,
                                   int32_t **start_cf_gb_xn,

                                   long *bnlins,
                                   uint64_t **blinvars,
                                   uint32_t **blineqs,
                                   uint64_t **bsquvars,

                                   fglm_data_t **bdata_fglm,
                                   fglm_bms_data_t **bdata_bms,

                                   int32_t *num_gb,
                                   int32_t **leadmons_ori,
                                   int32_t **leadmons_current,

                                   uint64_t bsz,
                                   param_t **nmod_params,
                                   /* trace_t **btrace, */
                                   bs_t *bs_qq,
                                   md_t *st,
                                   const int32_t fc,
                                   int info_level,
                                   bs_t **bs,
                                   int32_t *lmb_ori,
                                   int32_t dquot_ori,
                                   primes_t *lp,
                                   data_gens_ff_t *gens,
                                   double *stf4,
                                   const long nbsols,
                                   int *bad_primes)
{
    st->info_level = 0;

    double rt = realtime();
    /* tracing phase */
    len_t i;
    int32_t error = 0;

    /* F4 and FGLM are run using a single thread */
    /* st->nthrds is reset to its original value afterwards */
    const int nthrds = st->nthrds;
    st->nthrds = 1;

    memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));
#pragma omp parallel for num_threads(nthrds)  \
    private(i) schedule(static)
    for (i = 0; i < st->nprimes; ++i){
        bs[i] = core_gba(bs_qq, st, &error, lp->p[i]);
        *stf4 = realtime()-rt;
        /* printf("F4 trace timing %13.2f\n", *stf4); */

        if (error > 0) {
            if (bs[i] != NULL) {
                free(bs[i]);
                bs[i] = NULL;
            }
            nmod_params[i] = NULL;
            bad_primes[i] = 1;
            continue;
        }
        int32_t lml = bs[i]->lml;
        if (st->nev > 0) {
            int32_t j = 0;
            for (len_t k = 0; k < bs[i]->lml; ++k) {
                if (bs[i]->ht->ev[bs[i]->hm[bs[i]->lmps[k]][OFFSET]][0] == 0) {
                    bs[i]->lm[j]   = bs[i]->lm[k];
                    bs[i]->lmps[j] = bs[i]->lmps[k];
                    ++j;
                }
            }
            lml = j;
        }
        if(lml != num_gb[i]){
            if (bs[i] != NULL) {
                free_basis(&(bs[i]));
            }
            /* nmod_params[i] = NULL; */
            bad_primes[i] = 1;
            continue;
            /* return; */
        }
        get_lm_from_bs_trace(bs[i], bs[i]->ht, leadmons_current[i]);

        if(equal_staircase(leadmons_current[i], leadmons_ori[i],
                    num_gb[i], num_gb[i], bs[i]->ht->nv)){

            set_linear_poly(bnlins[i], blineqs[i], blinvars[i], bs[i]->ht,
                    leadmons_current[i], bs[i]);

            build_matrixn_from_bs_trace_application(bmatrix[i],
                    div_xn[i],
                    len_gb_xn[i],
                    start_cf_gb_xn[i],
                    lmb_ori, dquot_ori, bs[i], bs[i]->ht,
                    leadmons_ori[i], bs[i]->ht->nv,
                    lp->p[i]);
            if(nmod_fglm_compute_apply_trace_data(bmatrix[i], lp->p[i],
                        nmod_params[i],
                        bs[i]->ht->nv,
                        bsz,
                        bnlins[i], blinvars[i], blineqs[i],
                        bsquvars[i],
                        bdata_fglm[i],
                        bdata_bms[i],
                        nbsols,
                        info_level,
                        st)){
                bad_primes[i] = 1;
            }
        }
        else{
            bad_primes[i] = 1;
        }
        if (bs[i] != NULL) {
            free_basis_without_hash_table(&(bs[i]));
        }
    }
    st->nthrds = nthrds;
}

#if 0
static void modular_trace_application(sp_matfglm_t **bmatrix,
                                   int32_t **div_xn,
                                   int32_t **len_gb_xn,
                                   int32_t **start_cf_gb_xn,

                                   long *bnlins,
                                   uint64_t **blinvars,
                                   uint32_t **blineqs,
                                   uint64_t **bsquvars,

                                   fglm_data_t **bdata_fglm,
                                   fglm_bms_data_t **bdata_bms,

                                   int32_t *num_gb,
                                   int32_t **leadmons_ori,
                                   int32_t **leadmons_current,

                                   uint64_t bsz,
                                   param_t **nmod_params,
                                   /* trace_t **btrace, */
                                   ht_t **btht,
                                   bs_t *bs_qq,
                                   ht_t **bht,
                                   md_t *st,
                                   const int32_t fc,
                                   int info_level,
                                   bs_t **bs,
                                   int32_t *lmb_ori,
                                   int32_t dquot_ori,
                                   primes_t *lp,
                                   data_gens_ff_t *gens,
                                   double *stf4,
                                   const long nbsols,
                                   int *bad_primes){

  st->info_level = 0;

  /* tracing phase */
  len_t i;
  double ca0;
  int32_t error = 0;

  /* F4 and FGLM are run using a single thread */
  /* st->nthrds is reset to its original value afterwards */
  const int nthrds = st->nthrds;
  st->nthrds = 1 ;

  memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));
  #pragma omp parallel for num_threads(nthrds)  \
    private(i) schedule(static)
  for (i = 0; i < st->nprimes; ++i){
    ca0 = realtime();
    bs[i] = core_gba(bs_qq, st, &error, lp->p[i]);
    /* if(st->laopt > 40){
      bs[i] = modular_f4(bs_qq, bht[i], st, lp->p[i]);
    }
    else{
      bs[i] = gba_trace_application_phase(btrace[i], btht[i], bs_qq, bht[i], st, lp->p[i]);
    } */
    *stf4 = realtime()-ca0;
    /* printf("F4 trace timing %13.2f\n", *stf4); */

    if (bs[i] == NULL) {
        /* nmod_params[i] = NULL; */
        bad_primes[i] = 1;
        continue;
    }
    int32_t lml = bs[i]->lml;
    if (st->nev > 0) {
        int32_t j = 0;
        for (len_t k = 0; k < bs[i]->lml; ++k) {
            if ((*bht)->ev[bs[i]->hm[bs[i]->lmps[k]][OFFSET]][0] == 0) {
                bs[i]->lm[j]   = bs[i]->lm[k];
                bs[i]->lmps[j] = bs[i]->lmps[k];
                ++j;
            }
        }
        lml = j;
    }
    if(lml != num_gb[i]){
      if (bs[i] != NULL) {
        free_basis(&(bs[i]));
      }
      /* nmod_params[i] = NULL; */
      bad_primes[i] = 1;
      continue;
      /* return; */
    }
    get_lm_from_bs_trace(bs[i], bht[i], leadmons_current[i]);

    if(equal_staircase(leadmons_current[i], leadmons_ori[i],
                       num_gb[i], num_gb[i], bht[i]->nv)){

      set_linear_poly(bnlins[i], blineqs[i], blinvars[i], bht[i],
                      leadmons_current[i], bs[i]);

      build_matrixn_from_bs_trace_application(bmatrix[i],
                                              div_xn[i],
                                              len_gb_xn[i],
                                              start_cf_gb_xn[i],
                                              lmb_ori, dquot_ori, bs[i], bht[i],
                                              leadmons_ori[i], bht[i]->nv,
                                              lp->p[i]);
      if(nmod_fglm_compute_apply_trace_data(bmatrix[i], lp->p[i],
                                            nmod_params[i],
                                            bht[i]->nv,
                                            bsz,
                                            bnlins[i], blinvars[i], blineqs[i],
                                            bsquvars[i],
                                            bdata_fglm[i],
                                            bdata_bms[i],
                                            nbsols,
                                            info_level,
					    st)){
        bad_primes[i] = 1;
      }
    }
    else{
      bad_primes[i] = 1;
    }
    if (bs[i] != NULL) {
      free_basis(&(bs[i]));
    }
  }
  st->nthrds = nthrds;
}
#endif

/* sets function pointer */
void set_linear_function_pointer(int32_t fc){
  int nbits = 0;
  if(fc==0){
    nbits = 32;
  }
  else{
    if(fc < pow(2, 8)) {
      nbits = 8;
    }
    else{
      if(fc < pow (2, 16)){
        nbits = 16;
      }
      else nbits = 32;
    }
  }
  switch (nbits) {
  case 8:
    set_linear_poly = set_linear_poly_8;
    check_and_set_linear_poly = check_and_set_linear_poly_8;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_8;
    break;
  case 16:
    set_linear_poly = set_linear_poly_16;
    check_and_set_linear_poly = check_and_set_linear_poly_16;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_16;
    break;
  case 32:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
    break;
  case 0:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
    break;
  default:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
  }
}







/*

  - renvoie 0 si le calcul est ok.
  => GB = [1] dim =0 dquot = 0
  => Positive dimension dim > 0
  => Dimension zero + calcul qui a pu etre fait. dim=0 dquot > 0

  - renvoie 1 si le calcul a echoue
  => Dimension 0 => pas en position generique

  - renvoie 2 si besoin de plus de genericite.
  => (tous les carres ne sont pas sous l'escalier)

  - renvoie -2 si la carac est > 0

  - renvoie -3 si meta data pas bonnes

  - renvoie -4 si bad prime
*/

int msolve_trace_qq(mpz_param_t mpz_param,
                    param_t **nmod_param,
                    int *dim_ptr,
                    long *dquot_ptr,
                    data_gens_ff_t *gens,
                    int32_t ht_size, //initial_hts,
                    int32_t nr_threads,
                    int32_t max_nr_pairs,
                    int32_t elim_block_len,
                    int32_t reset_ht,
                    int32_t la_option,
                    int32_t use_signatures,
                    int32_t info_level,
                    int32_t print_gb,
                    int32_t pbm_file,
                    files_gb *files,
                    int round){

  const int32_t *lens = gens->lens;
  const int32_t *exps = gens->exps;
  uint32_t field_char = gens->field_char;
  const void *cfs = gens->mpz_cfs;
  if(gens->field_char){
    cfs = gens->cfs;
  }
  else{
    cfs = gens->mpz_cfs;
  }
  int mon_order = 0;
  int32_t nr_vars = gens->nvars;
  int32_t nr_gens = gens->ngens;
  int reduce_gb = 1;
  int32_t nr_nf = 0;
  const uint32_t prime_start = pow(2, 30);
  const int32_t nr_primes = nr_threads;

  len_t i;

  /* initialize stuff */
  md_t *st  = allocate_meta_data();

    int *invalid_gens   =   NULL;
    int res = validate_input_data(&invalid_gens, cfs, lens, &field_char, &mon_order,
            &elim_block_len, &nr_vars, &nr_gens, &nr_nf, &ht_size, &nr_threads,
            &max_nr_pairs, &reset_ht, &la_option, &use_signatures, &reduce_gb,
            &info_level);

    /* all data is corrupt */
    if (res == -1) {
        fprintf(stderr, "Invalid input generators, msolve now terminates.\n");
        free(invalid_gens);
        return -3;
    }
    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */

  if (check_and_set_meta_data_trace(st, lens, exps, cfs, invalid_gens,
              field_char, mon_order, elim_block_len, nr_vars, nr_gens,
              nr_nf, ht_size, nr_threads, max_nr_pairs, reset_ht, la_option,
              use_signatures, reduce_gb, prime_start, nr_primes, pbm_file,
              info_level)) {
    free(st);
    return -3;
  }

  /* lucky primes */
  primes_t *lp  = (primes_t *)calloc(st->nthrds, sizeof(primes_t));

  /*******************
  * initialize basis
  *******************/
  bs_t *bs_qq = initialize_basis(st);
  /* read in ideal, move coefficients to integers */
  import_input_data(bs_qq, st, lens, exps, cfs, invalid_gens);
  free(invalid_gens);
  invalid_gens  =   NULL;

  print_initial_statistics(stderr, st);

  /* for faster divisibility checks, needs to be done after we have
    * read some input data for applying heuristics */
  calculate_divmask(bs_qq->ht);

  /* sort initial elements, smallest lead term first */
  sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *),
          initial_input_cmp, bs_qq->ht);
  if(gens->field_char == 0){
    remove_content_of_initial_basis(bs_qq);
    /* generate lucky prime numbers */
    generate_lucky_primes(lp, bs_qq, st->prime_start, st->nthrds);
  }
  else{
    lp->old = 0;
    lp->ld = 1;
    lp->p = calloc(1, sizeof(uint32_t));
    normalize_initial_basis(bs_qq, st->fc);
  }

  /* generate array to store modular bases */
  bs_t **bs = (bs_t **)calloc((unsigned long)st->nthrds, sizeof(bs_t *));

  param_t **nmod_params =  (param_t **)calloc((unsigned long)st->nthrds,
                                              sizeof(param_t *));

  int *bad_primes = calloc((unsigned long)st->nthrds, sizeof(int));

  /* initialize tracers */
  /* trace_t **btrace = (trace_t **)calloc(st->nthrds,
                                       sizeof(trace_t *));
  btrace[0]  = initialize_trace(bs_qq, st); */
  /* initialization of other tracers is done through duplication */

  uint32_t prime = next_prime(1<<30);
  uint32_t primeinit;
  srand(time(0));
  prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  while(gens->field_char==0 && is_lucky_prime_ui(prime, bs_qq)){
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  }

  primeinit = prime;
  lp->p[0] = primeinit;

  if(gens->field_char){
    lp->p[0] = gens->field_char;
  }
  sp_matfglm_t **bmatrix = (sp_matfglm_t **)calloc(st->nthrds,
                                                  sizeof(sp_matfglm_t *));

  int32_t **bdiv_xn = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));
  int32_t **blen_gb_xn = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));
  int32_t **bstart_cf_gb_xn = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));

  fglm_data_t **bdata_fglm = (fglm_data_t **)calloc(st->nthrds,
                                                    sizeof(fglm_data_t *));
  fglm_bms_data_t **bdata_bms = (fglm_bms_data_t **)calloc(st->nthrds,
                                                          sizeof(fglm_bms_data_t *));
  int32_t *num_gb = (int32_t *)calloc(st->nthrds, sizeof(int32_t));
  int32_t **leadmons_ori = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));
  int32_t **leadmons_current = (int32_t**)calloc(st->nthrds, sizeof(int32_t *));

  uint64_t bsz = 0;

  /* data for linear forms */
  long nlins = 0;
  long *bnlins = (long *)calloc(st->nthrds, sizeof(long));
  uint64_t **blinvars = (uint64_t **)calloc(st->nthrds, sizeof(uint64_t *));
  uint64_t *linvars = calloc(bs_qq->ht->nv, sizeof(uint64_t));
  blinvars[0] = linvars;
  uint32_t **lineqs_ptr = calloc(st->nthrds, sizeof(uint32_t *));
  uint64_t **bsquvars = (uint64_t **) calloc(st->nthrds, sizeof(uint64_t *));
  uint64_t *squvars = calloc(nr_vars-1, sizeof(uint64_t));
  bsquvars[0] = squvars;

  set_linear_function_pointer(gens->field_char);

  int success = 1;
  int squares = 1;
#if 0
  int32_t *lmb_ori = modular_trace_learning(bmatrix, bdiv_xn, blen_gb_xn,
                                            bstart_cf_gb_xn,

                                            &nlins, blinvars[0], lineqs_ptr,
                                            squvars,

                                            bdata_fglm, bdata_bms,

                                            num_gb, leadmons_ori,

                                            &bsz, nmod_params,
                                            btrace[0],
                                            bs_qq, st,
                                            lp->p[0], //prime,
                                            st->info_level,
                                            print_gb,
                                            dim_ptr, dquot_ptr,
                                            gens,
                                            files,
                                            &success);
#else
  int32_t *lmb_ori = initial_modular_step(bmatrix, bdiv_xn, blen_gb_xn,
                                            bstart_cf_gb_xn,

                                            &nlins, blinvars[0], lineqs_ptr,
                                            squvars,

                                            bdata_fglm, bdata_bms,

                                            num_gb, leadmons_ori,

                                            &bsz, nmod_params,
                                            bs_qq, st,
                                            lp->p[0], //prime,
                                            print_gb,
                                            dim_ptr, dquot_ptr,
                                            gens,
                                            files,
                                            &success);
#endif
  if(*dim_ptr == 0 && success && *dquot_ptr > 0 && print_gb == 0){
    if(nmod_params[0]->elim->length - 1 != *dquot_ptr){
      for(int i = 0; i < nr_vars - 1; i++){
        if((squvars[i] == 0) && round ){
          squares = 0;
          success = 0;
        }
      }
    }
  }

  mpz_param->dim    = *dim_ptr;
  mpz_param->dquot  = *dquot_ptr;

  if(lmb_ori == NULL || success == 0 || print_gb || gens->field_char) {
      /* print_msolve_message(stderr, 1); */
    /* for(int i = 0; i < st->nthrds; i++){
      free_trace(&btrace[i]);
    }
    free(btrace); */
    if(gens->field_char==0){
      /* free_shared_hash_data(bht);
      if(bht!=NULL){
        free_hash_table(&bht);
      }
      free(bht); */
    }
    /* if(tht!=NULL){
      free_hash_table(&tht);
    }
    free(tht); */
    if (gens->field_char == 0) {
        for (i = 0; i < st->nthrds; ++i) {
            if (bs[i] != NULL) {
                free_basis_without_hash_table(&(bs[i]));
            }
        }
    }
    free(bs);
    if(gens->field_char==0){
      free_basis(&bs_qq);
    }
    //here we should clean nmod_params
    free_lucky_primes(&lp);
    free(bad_primes);
    free(lp);
    free(linvars);
    if(nlins){
      free(lineqs_ptr[0]);
    }
    free(bnlins);
    free(lineqs_ptr);
    free(squvars);
    if(print_gb){
      return 0;
    }
    if(*dim_ptr == 0 && gens->field_char && success){
      /* copy of parametrization */
      if(*dquot_ptr != 0){
        param_t *par = allocate_fglm_param(gens->field_char, st->nvars);
        nmod_poly_set(par->elim, nmod_params[0]->elim);
        nmod_poly_set(par->denom, nmod_params[0]->denom);
        for(long j = 0; j <= st->nvars - 2; j++){
          nmod_poly_set(par->coords[j], nmod_params[0]->coords[j]);
        }
        (*nmod_param) = par;

      }
      free(st);
      return 0;
    }
    free(st);
    free(nmod_params);
    if(*dim_ptr==1){
      if(info_level){
        fprintf(stderr, "Positive dimensional Grobner basis\n");
      }
      return 0;
    }
    if(*dquot_ptr==0){
      return 0;
    }
    if(*dquot_ptr>0){
      if(squares == 0){
        return 2;
      }
      return 1;
    }
  }
  crt_mpz_matfglm_t crt_mat;
  mpq_matfglm_t mpq_mat;
  mpz_matfglm_t mpz_mat;

#if LIFTMATRIX == 1
  crt_mpz_matfglm_initset(crt_mat, *bmatrix);
  mpq_matfglm_initset(mpq_mat, *bmatrix);
  mpz_matfglm_initset(mpz_mat, *bmatrix);
#endif

  /* btrace[0] = st->tr; */

  /* duplicate data for multi-threaded multi-mod computation */
  duplicate_data_mthread_trace(st->nthrds, bs_qq, st, num_gb,
                               leadmons_ori, leadmons_current,
                               /* btrace, */
                               bdata_bms, bdata_fglm,
                               bstart_cf_gb_xn, blen_gb_xn, bdiv_xn, bmatrix,
                               nmod_params, nlins, bnlins,
                               blinvars, lineqs_ptr,
                               bsquvars);

  normalize_nmod_param(nmod_params[0]);

  if(info_level && st->trace_level == APPLY_TRACER){
    fprintf(stderr, "\nStarts trace based multi-modular computations\n");
  }

  mpz_param_t tmp_mpz_param;
  mpz_param_init(tmp_mpz_param);

  initialize_mpz_param(mpz_param, nmod_params[0]);
  initialize_mpz_param(tmp_mpz_param, nmod_params[0]);
  //attention les longueurs des mpz_param sont fixees par nmod_params[0]
  //dans des cas exceptionnels, ca peut augmenter avec un autre premier.

  /* data for rational reconstruction of trace and det of mult. mat. */
  trace_det_fglm_mat_t trace_det;
  uint32_t detidx = 0;
  /* int32_t tridx = nmod_params[0]->elim->length-2; */
  int32_t tridx = 3* (nmod_params[0]->elim->length - 1) / 4 ;
  /* tridx = nmod_params[0]->elim->length - 2; */
  while(nmod_params[0]->elim->coeffs[tridx] == 0 && tridx > 0){
    tridx--;
  }
  detidx = 2 * (nmod_params[0]->elim->length - 1) / 3;
  while(nmod_params[0]->elim->coeffs[detidx] == 0 && detidx < nmod_params[0]->elim->length-2){
    detidx++;
  }

  trace_det_initset(trace_det,
                    nmod_params[0]->elim->coeffs[tridx],
                    nmod_params[0]->elim->coeffs[detidx],
                    tridx, detidx);
  /********************************************************************/

  mpz_t modulus;
  mpz_init_set_ui(modulus, primeinit);
  mpz_t prod_crt;
  mpz_init_set_ui(prod_crt, primeinit);

  mpq_t result, test;
  mpq_init(result);
  mpq_init(test);
  mpz_t rnum, rden; /* num and den of reconstructed rationals */
  mpz_init(rnum);
  mpz_init(rden);
  set_mpz_param_nmod(tmp_mpz_param, nmod_params[0]);


  long nsols = tmp_mpz_param->nsols;

  mpz_upoly_t numer;
  mpz_upoly_init2(numer, (nsols + 1),
                  32*(nsols + 1));
  numer->length = nsols + 1;

  mpz_upoly_t denom;
  mpz_upoly_init2(denom, (nsols + 1), 64*(nsols + 1));
  denom->length = nsols + 1;

  mpz_t guessed_den;
  mpz_init2(guessed_den, 32*nsols);
  mpz_t guessed_num;
  mpz_init2(guessed_num, 32*nsols);

  long maxrec = 0;
  long matrec = 0;

  int rerun = 1, nprimes = 1, mcheck =1;

  long nbadprimes = 0;

  int *is_lifted = malloc(sizeof(int)*nr_vars);
  for(int i = 0; i < nr_vars; ++i){
    is_lifted[i] = 0;
  }
  int mat_lifted = 0;
  int nbdoit = 1;
  int doit = 1;
  int prdone = 0;
  int lpow2 = 0;
  int clog = 0;
  int br = 0;
  prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));

  rrec_data_t recdata;
  initialize_rrec_data(recdata);

  /* measures time spent in rational reconstruction */
  double strat = 0;

  while(rerun == 1 || mcheck == 1){

    /* controls call to rational reconstruction */
    doit = ((prdone%nbdoit) == 0);

    /* generate lucky prime numbers */
    prime = next_prime(prime);
    lp->p[0] = prime;
    while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
      prime = next_prime(prime);
      lp->p[0] = prime;
    }

    for(len_t i = 1; i < st->nthrds; i++){
      prime = next_prime(prime);
      lp->p[i] = prime;
      while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
        prime = next_prime(prime);
        lp->p[i] = prime;
      }
    }
    prime = lp->p[st->nthrds - 1];

    double ca0 = realtime();

    double stf4 = 0;
    secondary_modular_steps(bmatrix,
                              bdiv_xn,
                              blen_gb_xn,
                              bstart_cf_gb_xn,

                              bnlins,
                              blinvars,
                              lineqs_ptr,
                              bsquvars,

                              bdata_fglm,
                              bdata_bms,
                              num_gb,
                              leadmons_ori,
                              leadmons_current,

                              bsz,
                              nmod_params, /* btrace, */
                              bs_qq, st,
                              field_char, 0, /* info_level, */
                              bs, lmb_ori, *dquot_ptr, lp,
                              gens, &stf4, nsols, bad_primes);
    double ca1 = realtime() - ca0;

    if(nprimes==1){
      if(info_level>2){
        fprintf(stderr, "------------------------------------------\n");
        fprintf(stderr, "#ADDITIONS       %13lu\n", (unsigned long)st->application_nr_add * 1000);
        fprintf(stderr, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->application_nr_mult * 1000);
        fprintf(stderr, "#REDUCTIONS      %13lu\n", (unsigned long)st->application_nr_red);
        fprintf(stderr, "------------------------------------------\n");
      }
      if(info_level>1){
        fprintf(stderr, "Application phase %.2f Gops/sec\n",
                (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4));
        fprintf(stderr, "Tracer + fglm time (elapsed): %.2f sec\n",
                (ca1) );
      }
    }

    for(int i = 0; i < st->nthrds; i++){
      if(bad_primes[i] == 0){
        normalize_nmod_param(nmod_params[i]);
      }
    }

    /* scrr measures time spent in ratrecon for modular images */
    double crr = 0, scrr = 0;
    /* CRT + rational reconstruction */
    for(len_t i = 0; i < st->nthrds; i++){
      if(bad_primes[i] == 0){
        if(rerun == 0){
          mcheck = check_param_modular(mpz_param, nmod_params[i], lp->p[i],
                                       is_lifted, trace_det, info_level);
        }
        crr = realtime();
        if(mcheck==1){
          br = new_rational_reconstruction(mpz_param,
                                           tmp_mpz_param,
                                           nmod_params[i],
                                           mpq_mat,
                                           crt_mat,
                                           mpz_mat,
                                           trace_det,
                                           bmatrix[i],
                                           numer, denom,
                                           &modulus, prod_crt,
                                           lp->p[i],
                                           &result,
                                           rnum, rden, recdata,
                                           &guessed_num, &guessed_den,
                                           &maxrec,
                                           &matrec,
                                           is_lifted,
                                           &mat_lifted,
                                           doit,
                                           st->nthrds, info_level);

#if LIFTMATRIX == 1
          if(mat_lifted && display){
            /* display_fglm_mpq_matrix(stderr, mpq_mat); */
            /* exit(1); */
            int maxnum = 0;
            int maxden = 0;
            int maxbits = 0;
            for(long i = 0; i < mpq_mat->nrows; i++){
              long c = 2*i*mpq_mat->ncols;
              for(long j = 0; j < mpq_mat->ncols; j++){
                maxnum = MAX(maxnum, mpz_sizeinbase(mpq_mat->dense_mat[c+2*j], 2));
                maxden = MAX(maxden, mpz_sizeinbase(mpq_mat->dense_mat[c+2*j+1], 2));
                maxbits = MAX(maxbits, mpz_sizeinbase(mpq_mat->dense_mat[c+2*j], 2) + mpz_sizeinbase(mpq_mat->dense_mat[c+2*j+1], 2));
              }
            }
            fprintf(stderr, "BIT SIZE IN MATRIX : %d, %d, %d => %f\n", maxnum, maxden, maxbits, (((double) maxbits))/16);
            fprintf(stderr, "nprimes = %d\n", nprimes);
            display = 0;
          }
#endif
          if(br == 1){
            rerun = 0;
          }
          else{
            rerun = 1;
          }
        }
        scrr += realtime()-crr;
        nprimes++;
      }
      else{
        if(info_level){
          fprintf(stderr, "<bp: %d>\n", lp->p[i]);
        }
        nbadprimes++;
        if(nbadprimes > nprimes){
          free(linvars);
          free(bnlins);
          free(lineqs_ptr[0]);
          free(lineqs_ptr);
          free(squvars);
          free_rrec_data(recdata);
          mpz_clear(prod_crt);
          trace_det_clear(trace_det);
          free_rrec_data(recdata);
          fprintf(stderr, "Many other data should be cleaned\n");
          return -4;
        }
      }
    }
    strat += scrr;

    double t = ((double)nbdoit)*ca1;
    if((t == 0) || (scrr >= 0.2*t && br == 0)){
      nbdoit=2*nbdoit;
      lpow2 = nprimes - lpow2;
      doit = 0;
      if(info_level){
        fprintf(stderr, "\n<Step:%d/%.2f/%.2f>",nbdoit,scrr,t);
      }
      prdone = 0;
    }
    else{
      prdone++;
    }

    if( (LOG2(nprimes) > clog) || (nbdoit != 1 && (nprimes % lpow2 == 0) ) ){
        if(info_level){
          fprintf(stderr, "{%d}", nprimes);
        }
        clog++;
    }

  }

  mpz_param->denom->length = mpz_param->nsols;
  for(long i = 1; i <= mpz_param->nsols; i++){
    mpz_set(mpz_param->denom->coeffs[i-1], mpz_param->elim->coeffs[i]);
    mpz_mul_ui(mpz_param->denom->coeffs[i-1], mpz_param->denom->coeffs[i-1], i);
  }


  if(info_level){
    fprintf(stderr, "\n%d primes used\n", nprimes);
    fprintf(stderr, "Time for CRT + rational reconstruction = %.2f\n", strat);
  }

  mpz_param_clear(tmp_mpz_param);

  mpz_upoly_clear(numer);
  mpz_upoly_clear(denom);
  mpz_clear(guessed_num);
  mpz_clear(guessed_den);
  mpq_clear(test);
  mpq_clear(result);
  mpz_clear(rnum);
  mpz_clear(rden);
  mpz_clear(modulus);
  mpz_clear(prod_crt);
  free_rrec_data(recdata);
  trace_det_clear(trace_det);

  /* free and clean up */
  /* free_shared_hash_data(bht); */
  /* for(int i = 0; i < st->nthrds; i++){
    free_hash_table(blht+i);
    free_hash_table(btht+i);
  } */

  //here we should clean nmod_params

  for(i = 0; i < st->nthrds; ++i){
    if (bs[i] != NULL) {
      free_basis(&(bs[i]));
    }
    free_fglm_bms_data(bdata_bms[i]);
    free_fglm_data(bdata_fglm[i]);
    free(blen_gb_xn[i]);
    free(bstart_cf_gb_xn[i]);
    free(bdiv_xn[i]);
    free(bmatrix[i]->dense_mat);
    free(bmatrix[i]->dense_idx);
    free(bmatrix[i]->triv_idx);
    free(bmatrix[i]->triv_pos);
    free(bmatrix[i]->dst);
    free(bmatrix[i]);
    free(leadmons_ori[i]);
    free(leadmons_current[i]);
    /* free_trace(&btrace[i]); */
    free(nmod_params[i]);

    free(blinvars[i]);
    free(lineqs_ptr[i]);
    free(bsquvars[i]);
  }

  free_basis(&(bs_qq));
  free(bs);
  free(bdata_fglm);
  free(bdata_bms);
  free(bmatrix);
  free(leadmons_ori);
  free(leadmons_current);
  free(nmod_params);

  free_lucky_primes(&lp);
  free_trace(&(st->tr));
  free(st);
  free(bad_primes);
  free(bnlins);
  free(blinvars);
  free(lineqs_ptr);
  free(bsquvars);
  free(is_lifted);
  free(num_gb);
  free(blen_gb_xn);
  free(bstart_cf_gb_xn);
  free(bdiv_xn);
  /* free(btrace); */

  return 0;
}


void real_point_init(real_point_t pt, long nvars){
  pt->nvars = nvars;
  pt->coords = malloc(sizeof(coord_t) * nvars);
  for(long i = 0; i < nvars; i++){
    mpz_init(pt->coords[i]->val_up);
    mpz_init(pt->coords[i]->val_do);
    pt->coords[i]->k_up = 0;
    pt->coords[i]->k_do = 0;
    pt->coords[i]->isexact = 0;
  }
}

void real_point_clear(real_point_t pt){
  for(long i = 0; i < pt->nvars; i++){
    mpz_clear(pt->coords[i]->val_up);
    mpz_clear(pt->coords[i]->val_do);
  }
  free(pt->coords);
}


void display_real_point_middle(FILE *fstream, real_point_t pt){
  mpz_t c;
  mpz_init(c);
  fprintf(fstream, "[");
  for(long i = 0; i < pt->nvars - 1; i++){
    mpz_add(c, pt->coords[i]->val_do, pt->coords[i]->val_up);
    mpz_out_str(fstream, 10, c);
    fprintf(fstream, " / ");
    fprintf(fstream, "2^%ld, ", pt->coords[i]->k_do + 1);
  }
  mpz_add(c, pt->coords[pt->nvars - 1]->val_do, pt->coords[pt->nvars - 1]->val_up);
  mpz_out_str(fstream, 10, c);
  fprintf(fstream, " / ");
  fprintf(fstream, "2^%ld ", pt->coords[pt->nvars - 1]->k_do + 1);
  fprintf(fstream, "]");
  mpz_clear(c);
}

void display_real_points_middle(FILE *fstream, real_point_t *pts, long nb){
  fprintf(fstream, "[");
  for(long i = 0; i < nb - 1; i++){
    display_real_point_middle(fstream, pts[i]);
    fprintf(fstream, ", ");
  }
  /* There might be no real solutions, so nb could be zero and
   * we have to recheck this here again. */
  if (nb > 0) {
    display_real_point_middle(fstream, pts[nb - 1]);
  }
  fprintf(fstream, "]:\n");
}

void display_real_point(FILE *fstream, real_point_t pt){

  fprintf(fstream, "[");
  for(long i = 0; i < pt->nvars - 1; i++){
    fprintf(fstream, "[");
    mpz_out_str(fstream, 10, pt->coords[i]->val_do);
    fprintf(fstream, " / ");
    fprintf(fstream, "2^%ld, ", pt->coords[i]->k_do);
    mpz_out_str(fstream, 10, pt->coords[i]->val_up);
    fprintf(fstream, " / ");
    fprintf(fstream, "2^%ld", pt->coords[i]->k_up);
    fprintf(fstream, "], ");
  }
  fprintf(fstream, "[");
  mpz_out_str(fstream, 10, pt->coords[pt->nvars - 1]->val_do);
  fprintf(fstream, " / ");
  fprintf(fstream, "2^%ld, ", pt->coords[pt->nvars - 1]->k_do);
  mpz_out_str(fstream, 10, pt->coords[pt->nvars - 1]->val_up);
  fprintf(fstream, " / ");
  fprintf(fstream, "2^%ld", pt->coords[pt->nvars - 1]->k_up);
  fprintf(fstream, "]");
  fprintf(fstream, "]");

}

void display_real_points(FILE *fstream, real_point_t *pts, long nb){
  fprintf(fstream, "[1,\n"); /* because at the moment we return a single list */
  fprintf(fstream, "[");
  for(long i = 0; i < nb - 1; i++){
    display_real_point(fstream, pts[i]);
    fprintf(fstream, ", ");
  }
  if(nb){
    display_real_point(fstream, pts[nb - 1]);
  }
  fprintf(fstream, "]\n");
  fprintf(fstream, "]");
}

void single_exact_real_root_param(mpz_param_t param, interval *rt, long nb,
                                  mpz_t *xdo, mpz_t *xup, mpz_t den_up, mpz_t den_do,
                                  mpz_t c, mpz_t tmp, mpz_t val_do, mpz_t val_up,
                                  mpz_t *tab, real_point_t pt, long prec,
                                  int info_level){

  mpz_poly_eval_2exp_naive(param->denom->coeffs, param->denom->length - 1,
                           &rt->numer, rt->k, tab, tab + 1);

  mpz_set(den_up, tab[0]);
  mpz_set(den_do, tab[0]);

  for(long nv = 0; nv < param->nvars - 1; nv++){

    mpz_poly_eval_2exp_naive(param->coords[nv]->coeffs,
                             param->coords[nv]->length - 1,
                             &rt->numer, rt->k, tab, tab + 1);
    mpz_set(val_up, tab[0]);
    mpz_set(val_do, tab[0]);

    mpz_neg(val_do, val_do);
    mpz_neg(val_up, val_up);
    mpz_swap(val_up, val_do);

    long exp = (rt->k) * ((param->denom->length ) -
                          (param->coords[nv]->length));

    mpz_mul_2exp(val_up, val_up, exp + prec);
    mpz_mul_2exp(val_do, val_do, exp + prec);
    mpz_mul(tab[1], den_up, param->cfs[nv]);
    mpz_cdiv_q(val_up, val_up, tab[1]);
    mpz_fdiv_q(val_do, val_do, tab[1]);

    mpz_set(pt->coords[nv]->val_up, val_up);
    mpz_set(pt->coords[nv]->val_do, val_do);
    pt->coords[nv]->k_up = prec;
    pt->coords[nv]->k_do = prec;
    pt->coords[nv]->isexact = 1;
  }

  mpz_set(pt->coords[param->nvars - 1]->val_do, rt->numer);
  mpz_set(pt->coords[param->nvars - 1]->val_up, rt->numer);
  pt->coords[param->nvars - 1]->k_up = rt->k;
  pt->coords[param->nvars - 1]->k_do = rt->k;
  pt->coords[param->nvars - 1]->isexact = 1;

}


/* assumes b is even */

void generate_table_values(interval *rt, mpz_t c,
                           const long ns, const long b,
                           const long corr,
                           mpz_t *xdo, mpz_t *xup){

  mpz_add_ui(c, rt->numer, 1);

  if(mpz_sgn(rt->numer) >= 0){

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for(long i = 1; i < ns; i++){
      if(i <= b){
        mpz_mul(xup[i], xup[i-1], c);
        mpz_mul(xdo[i], xdo[i-1], rt->numer);
      }
      if((i%b)==0 && i > b){
        long q = i / b;
        mpz_mul(xup[i], xup[(q-1)*b], xup[b]);
        mpz_mul(xdo[i], xdo[(q-1)*b], xdo[b]);
      }
    }
  }
  else{

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for(long i = 1; i < ns; i++){
      if(i<=b){
      if((i & 1) == 1){
        mpz_mul(xup[i], xdo[i-1], c);
        mpz_mul(xdo[i], xup[i-1], rt->numer);
      }
      else{
        mpz_mul(xup[i], xdo[i-1], rt->numer);
        mpz_mul(xdo[i], xup[i-1], c);
      }
      }
      if((i%b)==0 && i > b){
        long q = i / b;
        mpz_mul(xup[i], xdo[(q-1)*b], xup[b]);
        mpz_mul(xdo[i], xup[(q-1)*b], xdo[b]);
      }
    }
  }

  long q = (ns - 1) / b;

  for(long i = 1; i <= q; i++){
    mpz_mul_2exp(xup[i*b], xup[i*b], corr);
    mpz_cdiv_q_2exp(xup[i*b], xup[i*b], (rt->k)*i*b);

    mpz_mul_2exp(xdo[i*b], xdo[i*b], corr);
    mpz_fdiv_q_2exp(xdo[i*b], xdo[i*b], (rt->k)*i*b);
  }
}

/*
 xdo[i] and xup[i] are rt^i and c^i at precision 2^corr
[xdo[i]/2^corr, xup[i]/2^corr] contain [rt^i, c^i]
 */
void generate_table_values_full(interval *rt, mpz_t c,
                                const long ns, const long b,
                                const long corr,
                                mpz_t *xdo, mpz_t *xup){

  mpz_add_ui(c, rt->numer, 1);

  if(mpz_sgn(rt->numer) >= 0){

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for(long i = 1; i < ns; i++){

        mpz_mul(xup[i], xup[i-1], c);
        mpz_mul(xdo[i], xdo[i-1], rt->numer);
    }
  }
  else{

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for(long i = 1; i < ns; i++){
      if((i & 1) == 1){
        mpz_mul(xup[i], xdo[i-1], c);
        mpz_mul(xdo[i], xup[i-1], rt->numer);
      }
      else{
        mpz_mul(xup[i], xdo[i-1], rt->numer);
        mpz_mul(xdo[i], xup[i-1], c);
      }
    }
  }
  mpz_mul_2exp(xdo[0], xdo[0], corr);
  mpz_mul_2exp(xup[0], xup[0], corr);
  for(long i = 1; i <ns; i++){
    mpz_mul_2exp(xup[i], xup[i], corr);
    mpz_cdiv_q_2exp(xup[i], xup[i], (rt->k)*i);

    mpz_mul_2exp(xdo[i], xdo[i], corr);
    mpz_fdiv_q_2exp(xdo[i], xdo[i], (rt->k)*i);

  }

}



/*
  quad:=a*x^2+b*x+c;
  rr1:=r/2^k;
  rr2:=(r+1)/2^k;
  numer(expand(subs(x=1/(x+1),subs(x=x*(rr2-rr1),subs(x=x+rr1, quad)))));

  > numer(subs(x=x+r/2^k, quad)) =
  (2^k)^2*a*x^2+((2^k)^2*b+2*2^k*a*r)*x+c*(2^k)^2+2^k*b*r+a*r^2

  return 1 if quad has real roots in [0, 1] else it returns 0
 */
int evalquadric(mpz_t *quad, mpz_t r, long k,
                mpz_t *tmpquad, mpz_t tmp){

  /*We start by computing numer(subs(x=x*(rr2-rr1),subs(x=x+rr1, quad))) */
  /* This is  */
  /*a*x^2+(b*2^k+2*a*r)*x+c*(2^k)^2+2^k*b*r+a*r^2*/
  mpz_set(tmpquad[2], quad[2]);

  mpz_set(tmp, quad[2]);
  mpz_mul(tmp, tmp, r);
  /* tmp = a * r */
  mpz_set(tmpquad[0], tmp);
  mpz_mul(tmpquad[0], tmpquad[0], r);
  /* tmpquad[0] = a*r^2 */
  mpz_mul_2exp(tmp, tmp, 1);
  /* Now tmp = 2*a*r*/

  mpz_set(tmpquad[1], quad[1]);
  mpz_mul_2exp(tmpquad[1], tmpquad[1], k);
  /*tmpquad[1] = b*2^k*/
  mpz_add(tmpquad[1], tmpquad[1], tmp);
  /* Now tmpquad[1] = 2^k * b + 2 * a * r */

  mpz_set(tmp, quad[1]);
  /*tmp = b*/
  mpz_mul(tmp, tmp, r);
  mpz_mul_2exp(tmp, tmp, k);
  mpz_add(tmpquad[0], tmpquad[0], tmp);
  /* tmpquad[0] =  2^k*b*r + a*r^2*/

  mpz_set(tmp, quad[0]);
  mpz_mul_2exp(tmp, tmp, 2*k);
  /*tmp = c*2^(2k)*/
  mpz_add(tmpquad[0], tmpquad[0], tmp);

  /* At this stage, one has computed
     tmpquad = numer(subs(x=x*(rr2-rr1),subs(x=x+rr1, quad)))
     We want to know if it has real roots in [0, 1]
   */
  /* If all coefficients have the same sign, there is no root */
  int s=mpz_sgn(tmpquad[0]);
  if(s==mpz_sgn(tmpquad[1]) && s==mpz_sgn(tmpquad[2])){
    return 0;
  }
  /* One computes numer(subs(x=1/(x+1)), tmpquad) to apply Descartes*/
  /* if tmpquad = a2*x^2+a1*x+a0 computes
     a0*x^2+(2*a0+a1)*x+a0+a1+a2
  */
  mpz_add(tmpquad[1], tmpquad[1], tmpquad[0]);
  mpz_add(tmpquad[2], tmpquad[2], tmpquad[1]);
  mpz_add(tmpquad[1], tmpquad[1], tmpquad[0]);

  /* We apply now Descartes */
  s=mpz_sgn(tmpquad[0]);
  if(s==mpz_sgn(tmpquad[1]) && s==mpz_sgn(tmpquad[2])){
    return 0;
  }
  return 1;
}

/* evaluates denom (which has degree deg)
   at the interval [r/2^k, (r+1)/2^k]
   returns  */
int value_denom(mpz_t *denom, long deg, mpz_t r, long k,
                mpz_t *xdo, mpz_t *xup,
                mpz_t tmp, mpz_t den_do, mpz_t den_up,
                long corr, mpz_t c){

  /* /\* boo is 1 if den_do and den_up have not the same sign */
  /*    else it is 0 */
  mpz_add_ui(c, r, 1);

  /* mpz_poly_eval_2exp_naive2(denom, deg, r, k, den_do, tmp); */
  /* mpz_poly_eval_2exp_naive2(denom, deg, c, k, den_up, tmp); */


  /* if(mpz_sgn(den_do)!=mpz_sgn(den_up)){ */
  /*   return 1; */
  /* } */
  /* MODIFS START HERE */
  /* if(mpz_cmp(den_do, den_up)>0){ */
  /*   mpz_swap(den_do, den_up); */
  /* } */
  /* mpz_mul_2exp(den_do, den_do, corr); */
  /* mpz_mul_2exp(den_up, den_up, corr); */
  /* mpz_fdiv_q_2exp(den_do, den_do, k*deg); */
  /* mpz_cdiv_q_2exp(den_up, den_up, k*deg); */
  /* return 0; */

  int boo = mpz_poly_eval_interval(denom, deg, k,
                               r, c,
                               tmp, den_do, den_up);
  if(mpz_cmp(den_do, den_up) > 0){
    fprintf(stderr, "BUG (den_do > den_up)\n");
    exit(1);
  }
  mpz_mul_2exp(den_do, den_do, corr);
  mpz_mul_2exp(den_up, den_up, corr);
  mpz_fdiv_q_2exp(den_do, den_do, k*deg);
  mpz_cdiv_q_2exp(den_up, den_up, k*deg);

  if(mpz_sgn(den_do)!=mpz_sgn(den_up)){
    return 1;
  }
  return boo;
}




int newvalue_denom(mpz_t *denom, long deg, mpz_t r, long k,
                   mpz_t *xdo, mpz_t *xup,
                   mpz_t tmp, mpz_t den_do, mpz_t den_up,
                   long corr, mpz_t c){

  mpz_add_ui(c, r, 1);
  int boo = mpz_poly_eval_interval(denom, deg, k,
                               r, c,
                               tmp, den_do, den_up);
  if(mpz_cmp(den_do, den_up) > 0){
    fprintf(stderr, "BUG (den_do > den_up)\n");
    exit(1);
  }
  mpz_mul_2exp(den_do, den_do, corr);
  mpz_mul_2exp(den_up, den_up, corr);
  mpz_fdiv_q_2exp(den_do, den_do, k*deg);
  mpz_cdiv_q_2exp(den_up, den_up, k*deg);

  return boo;
}


void lazy_single_real_root_param(mpz_param_t param, mpz_t *polelim,
                                 interval *rt, long nb, interval *pos_root,
                                 mpz_t *xdo, mpz_t *xup, mpz_t den_up, mpz_t den_do,
                                 mpz_t c, mpz_t tmp, mpz_t val_do, mpz_t val_up,
                                 mpz_t *tab, real_point_t pt,
                                 long prec, long nbits, mpz_t s,
                                 int info_level){
  long ns = param->nsols ;
  /* root is exact */
  if(rt->isexact==1){
    single_exact_real_root_param(param, rt, nb,
                                 xdo, xup, den_up, den_do,
                                 c, tmp, val_do, val_up,
                                 tab, pt, prec,
                                 info_level);
    return ;
  }

  long b = 16;
  long corr = 2*(ns + rt->k);


  /* checks whether the abs. value of the root is greater than 1 */

  generate_table_values_full(rt, c, ns, b,
                             corr,
                             xdo, xup);

  while(newvalue_denom(param->denom->coeffs,
                       param->denom->length - 1,
                       rt->numer,
                       rt->k,
                       xdo, xup,
                       tmp, den_do, den_up, corr, s)){

    /* fprintf(stderr, "==> "); mpz_out_str(stderr, 10, rt->numer); */
    /* fprintf(stderr, " / 2^%ld\n", rt->k); */

    /* root is positive */
    if(mpz_sgn(rt->numer)>=0){
      get_values_at_bounds(param->elim->coeffs, ns, rt, tab);
      refine_QIR_positive_root(polelim, &ns, rt, tab, 2* (rt->k),
                               info_level);
    }
    else{
      /* root is positive */
      mpz_add_ui(pos_root->numer, rt->numer, 1);
      mpz_neg(pos_root->numer, pos_root->numer);
      pos_root->k = rt->k;
      pos_root->sign_left = - (rt->sign_left);
      pos_root->isexact = rt->isexact;
      for(long i = 0; i<=ns; i++){
        if((i & 1) == 1){
          mpz_neg(polelim[i], polelim[i]);
        }
      }
      get_values_at_bounds(polelim, ns, pos_root, tab);
      refine_QIR_positive_root(polelim, &ns, pos_root, tab, 2* (pos_root->k) + ns,
                               info_level);
      for(long i = 0; i<=ns; i++){
        if((i & 1) == 1){
          mpz_neg(polelim[i], polelim[i]);
        }
      }
      if(pos_root->isexact!=1){
        rt->k = pos_root->k;
        rt->isexact = pos_root->isexact;
        mpz_add_ui(rt->numer, pos_root->numer, 1);
        mpz_neg(rt->numer, rt->numer);
      }
      else{
        rt->k = pos_root->k;
        if(rt->isexact!=1){
          rt->isexact = pos_root->isexact;
          mpz_set(rt->numer, pos_root->numer);
          mpz_neg(rt->numer, rt->numer);
        }
      }
    }

    if(ns != param->nsols){
      for(long i = 0; i < param->elim->length; i++){
        mpz_set(polelim[i], param->elim->coeffs[i]);
      }
      ns = param->nsols;
    }

    corr *= 2;  /* *((rt->k) + prec); */
    b *= 2;
    /* generate_table_values(rt, c, ns, b, corr, */
    /*                       xdo, xup); */
    generate_table_values_full(rt, c, ns, b, corr,
                               xdo, xup);

    if(info_level){
      fprintf(stderr, "<%ld>", rt->k);
    }

  }


  mpz_t v1, v2;
  mpz_init(v1); mpz_init(v2);

  for(long nv = 0; nv < param->nvars - 1; nv++){

    /* lazy_mpz_poly_eval_interval(param->coords[nv]->coeffs, */
    /*                             param->coords[nv]->length - 1, */
    /*                             rt->k, */
    /*                             xdo, xup, */
    /*                             (rt->k)*(param->coords[nv]->length-1), corr, b, */
    /*                             tmp, val_do, val_up); */

    mpz_scalar_product_interval(param->coords[nv]->coeffs,
                                param->coords[nv]->length - 1,
                                rt->k,
                                xdo, xup,
                                tmp, val_do, val_up, corr);
    mpz_neg(val_do, val_do);
    mpz_neg(val_up, val_up);
    mpz_swap(val_up, val_do);

    /* long delta = param->denom->length - param->coords[nv]->length; */

    /* mpz_mul_2exp(val_do, val_do, corr * delta); */
    /* mpz_mul_2exp(val_up, val_up, corr * delta); */

    long dec = prec ;

    mpz_mul_2exp(val_up, val_up, dec);
    mpz_mul_2exp(val_do, val_do, dec);

    if(mpz_cmp(val_do, val_up) > 0){
      fprintf(stderr, "BUG in real root extractor(2)\n");
      exit(1);
    }

    if(mpz_sgn(den_do) >=0 && mpz_sgn(den_up) >= 0){
      if(mpz_sgn(val_do)>=0 && mpz_sgn(val_up) >= 0){
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_cdiv_q(v2, val_up, tmp);
      }
      if(mpz_sgn(val_do)<=0 && mpz_sgn(val_up) >= 0){
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_cdiv_q(v2, val_up, tmp);
      }
      if(mpz_sgn(val_do)<=0 && mpz_sgn(val_up) <= 0){
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_cdiv_q(v2, val_up, tmp);
      }
    }
    else{
      if(mpz_sgn(val_do)>=0 && mpz_sgn(val_up) >= 0){
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_up, tmp);
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_cdiv_q(v2, val_do, tmp);
      }
      if(mpz_sgn(val_do)<=0 && mpz_sgn(val_up) >= 0){
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_up, tmp);
        mpz_cdiv_q(v2, val_do, tmp);
      }
      if(mpz_sgn(val_do)<=0 && mpz_sgn(val_up) <= 0){
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_fdiv_q(v1, val_up, tmp);
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_cdiv_q(v2, val_do, tmp);
      }
    }
    mpz_set(val_do, v1);
    mpz_set(val_up, v2);

    mpz_set(pt->coords[nv]->val_up, val_up);
    mpz_set(pt->coords[nv]->val_do, val_do);

    pt->coords[nv]->k_up = prec;
    pt->coords[nv]->k_do = prec;
    pt->coords[nv]->isexact = 0;

  }
  mpz_set(pt->coords[param->nvars - 1]->val_do, rt->numer);
  mpz_set(pt->coords[param->nvars - 1]->val_up, rt->numer);
  mpz_add_ui(pt->coords[param->nvars - 1]->val_up,
             pt->coords[param->nvars - 1]->val_up, 1);
  pt->coords[param->nvars - 1]->k_up = rt->k;
  pt->coords[param->nvars - 1]->k_do = rt->k;
  pt->coords[param->nvars - 1]->isexact = 0;

  mpz_clear(v1);
  mpz_clear(v2);

}


void extract_real_roots_param(mpz_param_t param, interval *roots, long nb,
                              real_point_t *pts, long prec, long nbits,
                              double step, int info_level){
  long nsols = param->elim->length - 1;
  mpz_t *xup = malloc(sizeof(mpz_t)*nsols);
  mpz_t *xdo = malloc(sizeof(mpz_t)*nsols);
  mpz_t c, tmp, den_up, den_do, val_up, val_do;
  mpz_init(c);
  mpz_init(tmp);
  mpz_init(den_up);
  mpz_init(den_do);
  mpz_init(val_up);
  mpz_init(val_do);
  for(long i = 0; i < nsols; i++){
    mpz_init_set_ui(xup[i], 1);
    mpz_init_set_ui(xdo[i], 1);
  }
  mpz_t *tab = (mpz_t*)(calloc(8,sizeof(mpz_t)));//table for some intermediate values
  for(int i=0;i<8;i++)mpz_init(tab[i]);

  mpz_t *polelim = calloc(param->elim->length, sizeof(mpz_t));
  for(long i = 0; i < param->elim->length; i++){
    mpz_init_set(polelim[i], param->elim->coeffs[i]);
  }
  interval *pos_root = calloc(1, sizeof(interval));
  mpz_init(pos_root->numer);
  mpz_t s;
  mpz_init(s);

  double et = realtime();

  for(long nc = 0; nc < nb; nc++){
    interval *rt = roots+nc;

    lazy_single_real_root_param(param, polelim, rt, nb, pos_root,
                                xdo, xup, den_up, den_do,
                                c, tmp, val_do, val_up, tab,
                                pts[nc], prec, nbits, s,
                                info_level);

    if(info_level){
      if(realtime() - et >= step){
        fprintf(stderr, "{%.2f%%}", 100*nc/((double) nb));
        et = realtime();
      }
    }

  }

  for(long i = 0; i < nsols; i++){
    mpz_clear(xup[i]);
    mpz_clear(xdo[i]);
  }
  free(xup);
  free(xdo);
  mpz_clear(c);
  mpz_clear(s);
  mpz_clear(tmp);
  mpz_clear(den_up);
  mpz_clear(den_do);
  mpz_clear(val_up);
  mpz_clear(val_do);
  for(int i=0;i<8;i++)mpz_clear(tab[i]);
  free(tab);

  for(long i = 0; i < param->elim->length; i++){
    mpz_clear(polelim[i]);
  }
  free(polelim);
  mpz_clear(pos_root->numer);
  free(pos_root);
}


static real_point_t *isolate_real_roots_param(mpz_param_t param, long *nb_real_roots_ptr,
                                              interval **real_roots_ptr, 
                                              int32_t precision, int32_t nr_threads, int32_t info_level){
  mpz_t *pol = calloc(param->elim->length, sizeof(mpz_t));
  for(long i = 0; i < param->elim->length; i++){
    mpz_init_set(pol[i], param->elim->coeffs[i]);
  }
  long maxnbits = mpz_poly_max_bsize_coeffs(param->elim->coeffs,
                                            param->elim->length - 1);

  for(int i = 0; i < param->nvars - 1; i++){
    long cmax = mpz_poly_max_bsize_coeffs(param->coords[i]->coeffs,
                                          param->coords[i]->length - 1);
    maxnbits = MAX(cmax, maxnbits);
  }
  long prec = MAX(precision, 128 + (maxnbits) / 32 );
  double st = realtime();

  long unsigned int nbpos = 0;
  long unsigned int nbneg = 0;
  interval *roots = real_roots(pol, param->elim->length - 1,
                               &nbpos, &nbneg, prec, nr_threads, info_level );
  long nb = nbpos + nbneg;
  double step = (realtime() - st) / (nb) * 10 * LOG2(precision);

  real_point_t *pts = NULL;
  if(info_level > 0){
    fprintf(stderr, "Number of real roots: %ld\n", nb);
  }
  if(nb){
    /* */
    if(info_level){
      fprintf(stderr, "Starts real root extraction.\n");
    }
    double st = realtime();
    pts = malloc(sizeof(real_point_t) * nb);

    for(long i = 0; i < nb; i++){
      real_point_init(pts[i], param->nvars);
    }

    extract_real_roots_param(param, roots, nb, pts, precision, maxnbits,
                             step, info_level);
    if(info_level){
      fprintf(stderr, "Elapsed time (real root extraction) = %.2f\n",
              realtime() - st);
    }
  }
  *real_roots_ptr = roots;
  *nb_real_roots_ptr  = nb;

  for(long i = 0; i < param->elim->length; i++){
    mpz_clear(pol[i]);
  }
  free(pol);
  return pts;
}

void isolate_real_roots_lparam(mpz_param_array_t lparams, long **lnbr_ptr,
                               interval ***lreal_roots_ptr, real_point_t ***lreal_pts_ptr,
                               int32_t precision, int32_t nr_threads, int32_t info_level){
  long *lnbr = malloc(sizeof(long) * lparams->nb);
  interval **lreal_roots = malloc(sizeof(interval *) * lparams->nb);
  real_point_t **lreal_pts = malloc(sizeof(real_point_t *) * lparams->nb);
  for(int i = 0; i < lparams->nb; i++){
    lreal_roots[i] = NULL;
    lreal_pts[i] = NULL;
  }

  for(int i = 0; i < lparams->nb; i++){
    lreal_pts[i] = isolate_real_roots_param(lparams->params[i], lnbr + i,
                                            lreal_roots + i,
                                            precision, nr_threads, info_level);
  }
  (*lnbr_ptr)        = lnbr;
  (*lreal_roots_ptr) = lreal_roots;
  (*lreal_pts_ptr) = lreal_pts;

}

int real_msolve_qq(mpz_param_t mp_param,
                   param_t **nmod_param,
                   int *dim_ptr,
                   long *dquot_ptr,
                   long *nb_real_roots_ptr,
                   interval **real_roots_ptr,
                   real_point_t **real_pts_ptr,
                   data_gens_ff_t *gens,
                   int32_t ht_size, //initial_hts,
                   int32_t nr_threads,
                   int32_t max_nr_pairs,
                   int32_t elim_block_len,
                   int32_t reset_ht,
                   int32_t la_option,
                   int32_t use_signatures,
                   int32_t info_level,
                   int32_t print_gb,
                   int32_t pbm_file,
                   int32_t precision,
                   files_gb *files,
                   int round,
                   int32_t get_param){

  /*
    0 is comp. is ok
    1 if comp. failed
    2 if more genericity is required
    -2 if charac is > 0
    -3 if meta data are corrupted
    -4 if bad prime
  */

  double ct0 = cputime();
  double rt0 = realtime();
  int b = msolve_trace_qq(mp_param,
                          nmod_param,
                          dim_ptr,
                          dquot_ptr,
                          gens,
                          ht_size, //initial_hts,
                          nr_threads,
                          max_nr_pairs,
                          elim_block_len,
                          reset_ht,
                          la_option,
                          use_signatures,
                          info_level,
                          print_gb,
                          pbm_file,
                          files,
                          round);
  double ct1 = cputime();
  double rt1 = realtime();

  if(get_param>1){
    return b;
  }

  if(print_gb){
    return 0;
  }

  if(info_level){
    fprintf(stderr, "Time for rational param: %13.2f (elapsed) sec / %5.2f sec (cpu)\n\n",
            rt1 - rt0, ct1 - ct0);
  }

  real_point_t *pts = NULL;

  if(b==0 && *dim_ptr == 0 && *dquot_ptr > 0 && gens->field_char == 0){

    pts = isolate_real_roots_param(mp_param, nb_real_roots_ptr, real_roots_ptr,
                                   precision, nr_threads, info_level);
    int32_t nb = *nb_real_roots_ptr;
    if(nb){
      /* If we added a linear form for genericity reasons remove do not
       * return the last (new) variable in the solutions later on */
      if (gens->linear_form_base_coef > 0) {
        for (long i = 0; i < nb; ++i) {
          pts[i]->nvars--;
        }
      }
      /* If we changed the variable order for genericity reasons we have
       * to rechange the entries in the solution points. */
      if (gens->change_var_order != -1 &&
          gens->change_var_order != mp_param->nvars-1) {
        coord_t *tmp = malloc(sizeof(coord_t));
        int32_t lidx  = pts[0]->nvars - 1 - gens->change_var_order;
        for (long i = 0; i < nb; ++i) {
          memcpy(tmp,pts[i]->coords[0], sizeof(coord_t));
          memcpy(pts[i]->coords[0], pts[i]->coords[lidx], sizeof(coord_t));
          memcpy(pts[i]->coords[lidx], tmp, sizeof(coord_t));
        }
        free(tmp);
      }
      *real_pts_ptr = pts;
    }
  }
  return b;
}

void display_arrays_of_real_roots(files_gb *files, int32_t len, real_point_t **lreal_pts, long *lnbr){
  if(files->out_file != NULL){
    FILE *ofile = fopen(files->out_file, "a+");
    fprintf(ofile, "[");
    for(int i = 0; i < len - 1; i++){
      display_real_points(ofile, lreal_pts[i], lnbr[i]);
      fprintf(ofile, ", \n");
    }
    display_real_points(ofile, lreal_pts[len - 1], lnbr[len - 1]);
    fprintf(ofile, "];\n");
    fclose(ofile);
  }
  else{
    fprintf(stdout, "[");
    for(int i = 0; i < len - 1; i++){
      display_real_points(stdout, lreal_pts[i], lnbr[i]);
      fprintf(stdout, ", \n");
    }
    display_real_points(stdout, lreal_pts[len - 1], lnbr[len - 1]);
    fprintf(stdout, "];\n");
  }

}



void display_output(int b, int dim, int dquot,
                    files_gb *files, data_gens_ff_t *gens,
                    param_t *param, mpz_param_t *mpz_paramp, int get_param,
                    long *nb_real_roots_ptr,
                    interval **real_roots_ptr,
                    real_point_t **real_pts_ptr,
                    int info_level){
  if(dquot == 0){
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "a+");
      fprintf(ofile, "[-1]:\n");
      fclose(ofile);
    }
    else{
      fprintf(stdout, "[-1]:\n");
    }
    return;
  }

  if(dim == 0 && dquot >= 0){
    (*mpz_paramp)->nvars  = gens->nvars;
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "a+");
      fprintf(ofile, "[0, ");
      if (get_param >= 1 || gens->field_char) {
        mpz_param_out_str_maple(ofile, gens, dquot, *mpz_paramp, param);
      }
      if(get_param <= 1 && gens->field_char == 0){
        if(get_param){
          fprintf(ofile, ",");
        }
        display_real_points(
                            ofile, *real_pts_ptr, *nb_real_roots_ptr);
      }
      fprintf(ofile, "]:\n");
      fclose(ofile);
    }
    else{
      fprintf(stdout, "[0, ");
      if (get_param >= 1  || gens->field_char) {
        mpz_param_out_str_maple(stdout, gens, dquot, *mpz_paramp, param);
      }
      if(get_param <= 1 && gens->field_char == 0){
        if(get_param){
          fprintf(stdout, ",");
        }
        display_real_points(stdout, *real_pts_ptr,
                            *nb_real_roots_ptr);
      }
      fprintf(stdout, "]:\n");
    }
  }
  if(dim > 0){
    if (info_level > 0) {
      fprintf(stderr, "The ideal has positive dimension\n");
    }
    if(files->out_file != NULL){
      FILE *ofile2 = fopen(files->out_file, "a+");
      //1 because dim is >0
      fprintf(ofile2, "[1, %d, -1, []]:\n", gens->nvars);
      fclose(ofile2);
    }
    else{
      fprintf(stdout, "[1, %d, -1, []]:\n", gens->nvars);
    }
  }
}

void manage_output(int b, int dim, int dquot,
                   files_gb *files, data_gens_ff_t *gens,
                   param_t *param, mpz_param_t *mpz_paramp, int get_param,
                   long *nb_real_roots_ptr,
                   interval **real_roots_ptr,
                   real_point_t **real_pts_ptr,
                   int info_level){
  if(b == 0){
    display_output(b, dim, dquot, files, gens, param, mpz_paramp, get_param,
                   nb_real_roots_ptr,
                   real_roots_ptr,
                   real_pts_ptr,
                   info_level);
  }
  if(b==-2){
    fprintf(stderr, "Characteristic of the field here shouldn't be positive\n");
    (*mpz_paramp)->dim = -2;
  }
  if(b==-3){
    fprintf(stderr, "Problem when checking meta data\n");
    (*mpz_paramp)->dim  = -3;
  }

}


int core_msolve(
  int32_t la_option,
  int32_t use_signatures,
  int32_t nr_threads,
  int32_t info_level,
  int32_t initial_hts,
  int32_t max_pairs,
  int32_t elim_block_len,
  int32_t update_ht,
  int32_t generate_pbm,
  int32_t reduce_gb,
  int32_t print_gb,
  int32_t get_param,
  int32_t genericity_handling,
  int32_t saturate,
  int32_t colon,
  int32_t normal_form,
  int32_t normal_form_matrix,
  int32_t is_gb,
  int32_t precision,
  files_gb *files,
  data_gens_ff_t *gens,
  param_t **paramp,
  mpz_param_t *mpz_paramp,
  long *nb_real_roots_ptr,
  interval **real_roots_ptr,
  real_point_t **real_pts_ptr
  )
{

    param_t *param  = NULL;
    int32_t *bld    = NULL;
    int32_t **blen  = NULL;
    int32_t **bexp  = NULL;
    void **bcf      = NULL;
    int b           = 0;
    /* counter for randomly chosen linear forms */
    int round = -1;

restart:

    param = NULL;
    bld   = malloc(sizeof(int32_t));
    blen  = malloc(sizeof(int32_t *));
    bexp  = malloc(sizeof(int32_t *));
    bcf   = malloc(sizeof(void *));
    b     = 0;

    if(gens->field_char > 0){
        if (use_signatures > 0) {
            /* timings */
            double ct0, ct1, rt0, rt1;
            ct0 = cputime();
            rt0 = realtime();

            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            ht_t *bht   = NULL;
            md_t *st  = NULL;

            /* for (int ii = 0; ii<gens->nvars; ++ii) {
             *     mul[ii] = 1;
             * } */

            int success = 0;

            success = initialize_gba_input_data(&bs, &bht, &st,
                    gens->lens, gens->exps, (void *)gens->cfs,
                    1073741827, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, saturate, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    info_level);

            if (st->homogeneous != 1) {
                fprintf(stderr,
                        "Input system must be homogeneous.\n");
                exit(1);
            }

            st->gfc  = gens->field_char;
            if(info_level){
                fprintf(stderr,
                        "NOTE: Field characteristic is now corrected to %u\n",
                        st->gfc);
            }
            if (!success) {
                printf("Bad input data, stopped computation.\n");
                exit(1);
            }
            /* compute a gb for initial generators */
            success = core_sba_schreyer(&bs, &bht, &st);

            if (!success) {
                printf("Problem with sba, stopped computation.\n");
                exit(1);
            }
            int64_t nb  = export_results_from_gba(bld, blen, bexp,
                    bcf, &malloc, &bs, &bht, &st);

            /* timings */
            ct1 = cputime();
            rt1 = realtime();
            st->f4_ctime = ct1 - ct0;
            st->f4_rtime = rt1 - rt0;

            get_and_print_final_statistics(stderr, st, bs);

            if(nb==0){
                fprintf(stderr, "Something went wrong during the computation\n");
                return -1;
            }
            return 0;
        }

        if (saturate == 1) { /* positive characteristic */
            /* timings */
            double ct0, ct1, rt0, rt1;
            ct0 = cputime();
            rt0 = realtime();

            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            bs_t *sat   = NULL;
            ht_t *bht   = NULL;
            md_t *st  = NULL;

            /* for (int ii = 0; ii<gens->nvars; ++ii) {
             *     mul[ii] = 1;
             * } */

            int32_t error = 0;
            int success   = 0;

            /*             initialize generators of ideal, note the "gens->ngens-normal_form" which
             *             means that we only take the first nr_gens-normal_form generators from
             *             the input file, the last normal_form polynomial in the file will
             *             be reduced w.r.t. the basis
             *
             *             NOTE: There is a little hack here, instead of gens->field_char we
             *             give 1073741827 as parameter, which ensures that all F4 internal
             *             routines are the 32-bit implementations (since nf is at the moment
             *             only implemented for 32-bit elements). Later on we set st-fc by hand
             *             to the correct field characteristic. */
            success = initialize_gba_input_data(&bs, &bht, &st,
                    gens->lens, gens->exps, (void *)gens->cfs,
                    1073741827, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, saturate, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    info_level);

            st->gfc  = gens->field_char;
            set_ff_bits(st, st->gfc);
            if(info_level){
                fprintf(stderr,
                        "NOTE: Field characteristic is now corrected to %u\n",
                        st->gfc);
            }
            if(st->ff_bits < 32){
              fprintf(stderr, "Error: not implemented yet (prime field of too low characteristic)\n");
              return 1;
            }
            if (!success) {
                printf("Bad input data, stopped computation.\n");
                exit(1);
            }

            if (is_gb == 1) {
                for (len_t k = 0; k < bs->ld; ++k) {
                    bs->lmps[k] = k;
                    bs->lm[k]   = bht->hd[bs->hm[k][OFFSET]].sdm;
                    bs->lml     = bs->ld;
                }
            } else {
                sat = initialize_basis(st);
                import_input_data_nf_ff_32(
                                           sat, bht, st, gens->ngens-saturate, gens->ngens,
                                           gens->lens, gens->exps, (void *)gens->cfs);

                sat->ld = sat->lml  =  saturate;
                /* normalize_initial_basis(tbr, st->gfc); */
                for (int k = 0; k < saturate; ++k) {
                    sat->lmps[k]  = k; /* fix input element in tbr */
                }

                /* compute a gb for initial generators */
                success = core_f4sat(bs, sat, st, &error);

                if (!success) {
                    printf("Problem with f4sat, stopped computation.\n");
                    exit(1);
                }
                int64_t nb  = export_results_from_gba(bld, blen, bexp,
                        bcf, &malloc, &bs, &bht, &st);

                /* timings */
                ct1 = cputime();
                rt1 = realtime();
                st->f4_ctime = ct1 - ct0;
                st->f4_rtime = rt1 - rt0;

                get_and_print_final_statistics(stderr, st, bs);

                if(nb==0){
                    fprintf(stderr, "Something went wrong during the computation\n");
                    return -1;
                }
                if (print_gb) {
		  print_ff_basis_data(files->out_file, "a", bs, bht,
                    st, gens, print_gb);
		}
            }
            return 0;
        }
	if (colon == 1) {
	    /* colon is 1 */
            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            bs_t *tbr   = NULL;
            ht_t *bht   = NULL;
            md_t *st  = NULL;

	    double ct0, rt0, ct0p, rt0p, ct1, rt1, ct2, rt2, ct3, rt3, ct4, rt4;
	    ct0 = cputime();
	    rt0 = realtime();

            /* generate array for storing multiplier for polynomial
             * to be reduced by basis */

            exp_t *mul  = (exp_t *)calloc(gens->nvars, sizeof(exp_t));

            /* for (int ii = 0; ii<gens->nvars; ++ii) {
             *     mul[ii] = 1;
             * } */
            int32_t err = 0;

            /* initialize generators of ideal, note the "gens->ngens-1_form" which
             * means that we only take the first nr_gens-1 generators from
             * the input file, the last polynomial in the file will
             * be reduced w.r.t. the basis
             *
             * NOTE: There is a little hack here, instead of gens->field_char we
             * give 1073741827 as parameter, which ensures that all F4 internal
             * routines are the 32-bit implementations (since nf is at the moment
             * only implemented for 32-bit elements). Later on we set st-fc by hand
             * to the correct field characteristic. */
            int success = initialize_gba_input_data(&bs, &bht, &st,
                    gens->lens, gens->exps, (void *)gens->cfs,
                    1073741827, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, 1, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    info_level);

	    st->gfc  = gens->field_char;
            if(info_level){
                fprintf(stderr,
                        "NOTE: Field characteristic is now corrected to %u\n",
                        st->gfc);
            }
            if (!success) {
                printf("Bad input data, stopped computation.\n");
                exit(1);
            }

	    ct0p = cputime();
	    rt0p = realtime();
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "INIT   TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt0p-rt0, ct0p-ct0);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }
            if (is_gb == 1) {
                for (len_t k = 0; k < bs->ld; ++k) {
                    bs->lmps[k] = k;
                    bs->lm[k]   = bht->hd[bs->hm[k][OFFSET]].sdm;
                    bs->lml     = bs->ld;
                }
            } else {

                /* compute a gb for initial generators */
                bs = core_gba(bs, st, &err, gens->field_char);

                if (err) {
                    printf("Problem with F4, stopped computation.\n");
                    exit(1);
                }
            }
	    export_results_from_gba(bld, blen, bexp,
						  bcf, &malloc, &bs, &bht, &st);
            printf("size of basis: %u\n", bs->lml);
	    ct1 = cputime();
	    rt1 = realtime();
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "F4     TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt1-rt0p, ct1-ct0p);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }

            /* initialize data for elements to be reduced,
             * NOTE: Don't initialize BEFORE running core_f4, bht may
             * change, so hash values of tbr may become wrong. */
            tbr = initialize_basis(st);
            import_input_data_nf_ff_32(tbr, bht, st, gens->ngens-1, gens->ngens,
				       gens->lens, gens->exps, (void *)gens->cfs);
            tbr->ld = tbr->lml  =  1;
            /* normalize_initial_basis(tbr, st->gfc); */
            for (int k = 0; k < 1; ++k) {
                tbr->lmps[k]  = k; /* fix input element in tbr */
            }

            /* compute normal form of last element in tbr */

            tbr = core_nf(tbr, st, mul, bs, &err);

            if (err) {
                printf("Problem with normalform, stopped computation.\n");
                exit(1);
            }
            /* print reduced element in tbr, last one is the input element */
	    /* printf ("normal form:\n"); */
            /* print_msolve_polynomials_ff(stdout, 1, tbr->lml, tbr, bht, */
	    /* 				st, gens->vnames, 0); */
	    ct2 = cputime();
	    rt2 = realtime();
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "NF     TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt2-rt1, ct2-ct1);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }
            /* print all reduced elements in tbr, first  one
             * is the input element */
            /* print_msolve_polynomials_ff(stdout, 1, tbr->lml, tbr, bht, */
	    /* 				st, gens->vnames, 0); */
	    /* printf("\n"); */
	    /* list of monomials */
	    /* size of the list */
	    long suppsize= tbr->hm[tbr->lmps[1]][LENGTH]; // bs->hm[bs->lmps[1]][LENGTH]
	    printf("Length of the support of phi: %lu\n",
		   suppsize);
	    
	    /* sht and hcm will store the support of the normal form in tbr. */
	    ht_t *sht   = initialize_secondary_hash_table(bht, st);
	    hi_t *hcm   = (hi_t *)malloc(sizeof(hi_t));
	    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));
	    
	    /* printf("Starts computation of normal form matrix\n"); */
	    get_normal_form_matrix(tbr, bht, 1,
				   st, &sht, &hcm, &mat);
	    printf("Length of union of support of all normal forms: %u\n",
		   mat->nc);
	    
	    /* printf("\nUnion of support, sorted by decreasing monomial order:\n"); */
	    /* for (len_t k = 0; k < mat->nc; ++k) { */
	    /*   for (len_t l = 1; l <= sht->nv; ++l) { */
	    /* 	printf("%2u ", sht->ev[hcm[k]][l]); */
	    /*   } */
	    /*   printf("\n"); */
	    /* } */

	    int32_t *bcf_ff = (int32_t *)(*bcf);
	    int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);
	    
	    long maxdeg = sht->ev[hcm[0]][0]; /* degree of the normal
						 form */
	    for (long i = 0; i < bld[0]; i++) {
	      long degi = 0;
	      for (long k = 0; k < gens->nvars; k++) {
		degi += bexp_lm[i*gens->nvars+k];
	      }
	      maxdeg = MAX(maxdeg,degi);
	    }
	    printf ("Maximal degree of the truncated staircase: %ld\n",maxdeg);
	    long dquot;
#define ZERO 0
#if ZERO
	    int32_t *lmb= monomial_basis_colon (bld[0],
						gens->nvars, bexp_lm,
						&dquot, maxdeg);
#else
	    int32_t *lmb= monomial_basis_colon_no_zero (bld[0],
						 gens->nvars, bexp_lm,
						 &dquot, maxdeg);
#endif
	    /* printf("\nMonomial basis:\n"); */
	    /* for (len_t k = 0; k < dquot; ++k) { */
	    /*   for (len_t l = 0; l < gens->nvars; ++l){ */
	    /* 	printf("%2u ", lmb[k*gens->nvars+l]); */
	    /*   } */
	    /*   printf("\n"); */
	    /* } */
	    printf("Subspace of quotient algebra has dimension: %ld\n",dquot);
	    uint32_t * leftvector = calloc(dquot,sizeof (uint32_t));
	    uint32_t ** leftvectorsparam = malloc(2*(gens->nvars-1)*sizeof (uint32_t *));
	    for (long i = 0; i < 2*(gens->nvars-1); i++) {
	      leftvectorsparam[i] = calloc(dquot,sizeof (uint32_t));
	    }

	    /* we assume that the support of phi is enough to encode
	     * the multiplication matrix */
	    /* we need the nf of sigma x_n for all sigma is this
	       support */
#if ZERO
	    sp_matfglmcol_t  *matrix
	      = build_matrixn_colon(lmb, dquot, bld[0], blen, bexp, bcf_ff,
				    bexp_lm, tbr, bht, st, mul, bs,
				    gens->nvars, gens->field_char,
				    maxdeg, gens,
				    leftvector, leftvectorsparam,
				    suppsize);
#else
	    sp_matfglmcol_t  *matrix
	      = build_matrixn_colon_no_zero(lmb, dquot, bld[0], blen, bexp, bcf_ff,
					    bexp_lm, tbr, bht, st, mul, bs,
					    gens->nvars, gens->field_char,
					    maxdeg, gens,
					    leftvector, leftvectorsparam,
					    suppsize);
#endif
#undef ZERO
	    ct3 = cputime();
	    rt3 = realtime();
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "MATRIX TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt3-rt2, ct3-ct2);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }
	    uint64_t *linvars = calloc(gens->nvars, sizeof(uint64_t));
	    uint32_t *lineqs = calloc(gens->nvars,sizeof(uint32_t));
	    uint64_t *squvars = calloc(gens->nvars-1, sizeof(uint64_t));
	    param_t * param = nmod_fglm_guess_colon(matrix, gens->field_char,
						    leftvector, leftvectorsparam,
						    gens->nvars,
						    0, linvars, lineqs, squvars, 1, st);
	    ct4 = cputime();
	    rt4 = realtime();
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "FGLM  TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt4-rt3, ct4-ct3);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }
	    display_fglm_param(stdout, param);
	    if (info_level) {
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "TOTAL  TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
		      rt4-rt0, ct4-ct0);
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	      fprintf(stderr, "INIT   PERCTG %13.2f%% (REAL) / %5.2f%% (CPU)\n",
		      100*(rt0p-rt0)/(rt4-rt0), 100*(ct0p-ct0)/(ct4-ct0));
	      fprintf(stderr, "F4     PERCTG %13.2f%% (REAL) / %5.2f%% (CPU)\n",
		      100*(rt1-rt0p)/(rt4-rt0), 100*(ct1-ct0p)/(ct4-ct0));
	      fprintf(stderr, "NF     PERCTG %13.2f%% (REAL) / %5.2f%% (CPU)\n",
		      100*(rt2-rt1)/(rt4-rt0), 100*(ct2-ct1)/(ct4-ct0));
	      fprintf(stderr, "MATRIX PERCTG %13.2f%% (REAL) / %5.2f%% (CPU)\n",
		      100*(rt3-rt2)/(rt4-rt0), 100*(ct3-ct2)/(ct4-ct0));
	      fprintf(stderr, "FGLM   PERCTG %13.2f%% (REAL) / %5.2f%% (CPU)\n",
		      100*(rt4-rt3)/(rt4-rt0), 100*(ct4-ct3)/(ct4-ct0));
	      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
	    }
	    free(param);
	    free(squvars);
	    free(lineqs);
	    free(linvars);
	    free(matrix);
	    for (long i = 0; i < 2*(gens->nvars-1); i++) {
	      free(leftvectorsparam[i]);
	    }
	    free(leftvectorsparam);
	    free(leftvector);
	    free(lmb);
	    free(bexp_lm);
	    free(bcf_ff);
	    free(hcm);
	    free (mul);
	    hcm = NULL;
	    if (sht != NULL) {
	      free_hash_table(&sht);
	    }
	    
            /* free and clean up */
            if (bs != NULL) {
	      free_basis(&bs);
            }
            if (tbr != NULL) {
	      free_basis(&tbr);
	    }
            free(st);
            st  = NULL;
            free_shared_hash_data(bht);
            if (bht != NULL) {
	      free_hash_table(&bht);
            }
	    return 0;
	}
	/* no saturate = 0 */
	/* no colon    = 0 */
	if (normal_form == 0) {/* positive characteristic */

	  int dim = - 2;
	  long dquot = -1;

    if(elim_block_len > 0 && print_gb == 0){
      if(info_level){
        fprintf(stderr, "Warning: elim order not available for rational parametrizations\n");
        fprintf(stderr, "Computing Groebner basis\n");
        print_gb=2;
      }
    }
	  b = real_msolve_qq(*mpz_paramp,
                       &param,
                       &dim,
                       &dquot,
                       nb_real_roots_ptr,
                       real_roots_ptr,
                       real_pts_ptr,
                       gens,
                       initial_hts, nr_threads, max_pairs,
                       elim_block_len, update_ht,
                       la_option, use_signatures, info_level, print_gb,
                       generate_pbm, precision, files, round, get_param);
          if(print_gb){
            return 0;
          }

          manage_output(b, dim, dquot, files, gens, param, mpz_paramp, get_param,
                        nb_real_roots_ptr,
                        real_roots_ptr,
                        real_pts_ptr,
                        info_level);


          /* if (b == 0 && gens->field_char > 0) { */
          /*   if(dim == 0){ */
          /*     if(files->out_file != NULL){ */
          /*       FILE *ofile = fopen(files->out_file, "a"); */
          /*       if(dquot == 0){ */
          /*         fprintf(ofile, "[-1]:\n"); */
          /*         return 0; */
          /*       } */
          /*       display_fglm_param_maple(ofile, param); */
          /*       fclose(ofile); */
          /*     } */
          /*     else{ */
          /*       if(dquot == 0){ */
          /*         fprintf(stdout, "[-1]:\n"); */
          /*         return 0; */
          /*       } */
          /*       display_fglm_param_maple(stdout, param); */
          /*     } */
          /*     return 0; */
          /*   } */
          /* } */
          if (b == 1) {
            free(bld);
            bld = NULL;
            free(blen);
            blen  = NULL;
            free(bexp);
            bexp  = NULL;
            free(bcf);
            bcf = NULL;
            free(param);
            param = NULL;
            if (genericity_handling > 0) {
              if (change_variable_order_in_input_system(gens, info_level)) {
                goto restart;
              }
              if (genericity_handling == 2) {
                if (add_linear_form_to_input_system(gens, info_level)) {
                  goto restart;
                }
              }
            }
            fprintf(stderr, "\n=====> Computation failed <=====\n");
            fprintf(stderr, "Try to add a random linear form with ");
            fprintf(stderr, "a new variable\n");
            fprintf(stderr, "(smallest w.r.t. DRL) to the input system. ");
            fprintf(stderr, "This will\n");
            fprintf(stderr, "be done automatically if you run msolve with option\n");
            fprintf(stderr, "\"-c2\" which is the default.\n");
          }
          if(b == 2){
            free(bld);
            bld = NULL;
            free(blen);
            blen  = NULL;
            free(bexp);
            bexp  = NULL;
            free(bcf);
            bcf = NULL;
            free(param);
            param = NULL;
            round++;
            if(gens->change_var_order >= 0){
              undo_variable_order_change(gens);
            }
            if (add_random_linear_form_to_input_system(gens, info_level)) {
              goto restart;
            }
          }
        }
	else {
          /* normal_form is 1 */
        printf("normal form active? %d\n", normal_form);
            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            bs_t *tbr   = NULL;
            ht_t *bht   = NULL;
            md_t *st  = NULL;

            /* generate array for storing multiplier for polynomial
             * to be reduced by basis */

            exp_t *mul  = (exp_t *)calloc(gens->nvars, sizeof(exp_t));

            /* for (int ii = 0; ii<gens->nvars; ++ii) {
             *     mul[ii] = 1;
             * } */

            int32_t err = 0;
            int success = 0;

            /* initialize generators of ideal, note the "gens->ngens-normal_form" which
             * means that we only take the first nr_gens-normal_form generators from
             * the input file, the last normal_form polynomial in the file will
             * be reduced w.r.t. the basis
             *
             * NOTE: There is a little hack here, instead of gens->field_char we
             * give 1073741827 as parameter, which ensures that all F4 internal
             * routines are the 32-bit implementations (since nf is at the moment
             * only implemented for 32-bit elements). Later on we set st-fc by hand
             * to the correct field characteristic. */
            success = initialize_gba_input_data(&bs, &bht, &st,
                    gens->lens, gens->exps, (void *)gens->cfs,
                    1073741827, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, normal_form, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    info_level);

            st->gfc  = gens->field_char;
            if(info_level){
                fprintf(stderr,
                        "NOTE: Field characteristic is now corrected to %u\n",
                        st->gfc);
            }
            if (!success) {
                printf("Bad input data, stopped computation.\n");
                exit(1);
            }

            if (is_gb == 1) {
                for (len_t k = 0; k < bs->ld; ++k) {
                    bs->lmps[k] = k;
                    bs->lm[k]   = bht->hd[bs->hm[k][OFFSET]].sdm;
                    bs->lml     = bs->ld;
                }
            } else {

                /* compute a gb for initial generators */
                bs = core_gba(bs, st, &err, gens->field_char);

                if (err) {
                    printf("Problem with F4, stopped computation.\n");
                    exit(1);
                }
            }
            printf("size of basis %u\n", bs->lml);
            /* initialize data for elements to be reduced,
             * NOTE: Don't initialize BEFORE running core_f4, bht may
             * change, so hash values of tbr may become wrong. */
            tbr = initialize_basis(st);
            import_input_data_nf_ff_32(
                    tbr, bht, st, gens->ngens-normal_form, gens->ngens,
                    gens->lens, gens->exps, (void *)gens->cfs);
            tbr->ld = tbr->lml  =  normal_form;
            /* normalize_initial_basis(tbr, st->gfc); */
            for (int k = 0; k < normal_form; ++k) {
                tbr->lmps[k]  = k; /* fix input element in tbr */
            }
            /* compute normal form of last element in tbr */
            tbr = core_nf(tbr, st, mul, bs, &err);

            if (err) {
                printf("Problem with normalform, stopped computation.\n");
                exit(1);
            }
            /* print all reduced elements in tbr, first normal_form ones
             * are the input elements */
	    printf ("normal form:\n");
            print_msolve_polynomials_ff(stdout, normal_form,
					tbr->lml, tbr, bht, st, gens->vnames, 0);
            if (normal_form_matrix > 0) {
                /* sht and hcm will store the union of the support
                 * of all normal forms in tbr. */
                ht_t *sht   = initialize_secondary_hash_table(bht, st);
                hi_t *hcm   = (hi_t *)malloc(sizeof(hi_t));
                mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

                printf("\nStarts computation of normal form matrix\n");
                get_normal_form_matrix(tbr, tbr->ht, normal_form,
                        st, &sht, &hcm, &mat);

                printf("\n\nLength of union of support of all normal forms: %u\n",
                        mat->nc);

                printf("\nUnion of support, sorted by decreasing monomial order:\n");
                for (len_t k = 0; k < mat->nc; ++k) {
                    for (len_t l = 1; l <= sht->nv; ++l) {
                        printf("%2u ", sht->ev[hcm[k]][l]);
                    }
                    printf("\n");
                }

                /* sparse represented matrix of normal forms, note that the column entries
                 * of the rows are not sorted, but you can do so using any sort algorithm */
                printf("\nMatrix of normal forms (sparse format, note that entries are\n");
                printf("NOT completely sorted by column index):\n");
                int64_t nterms  = 0;
                for (len_t k = 0; k < mat->nr; ++k) {
                    printf("row %u | ", k);
                    for (len_t l = 0; l < mat->tr[k][LENGTH]; ++l) {
                        printf("%u at %u, ",
                                tbr->cf_32[mat->tr[k][COEFFS]][l],
                                mat->tr[k][l+OFFSET]);
                    }
                    printf("\n");
                    nterms  +=  mat->tr[k][LENGTH];
                }
                nterms  *=  100; /* for percentage */
                double density = (double)nterms / (double)mat->nr / (double)mat->nc;
                printf("\nMatrix of normal forms (dense format)\n");
                cf32_t *dr  = (cf32_t *)malloc(
                        (unsigned long)mat->nc * sizeof(cf32_t));
                for (len_t k = 0; k < mat->nr; ++k) {
                    memset(dr, 0, (unsigned long)mat->nc * sizeof(cf32_t));
                    printf("row %u | ", k);
                    for (len_t l = 0; l < mat->tr[k][LENGTH]; ++l) {
                        dr[mat->tr[k][l+OFFSET]]  =
                            tbr->cf_32[mat->tr[k][COEFFS]][l];
                    }
                    for (len_t l = 0; l < mat->nc; ++l) {
                        printf("%u, ", dr[l]);
                    }
                    printf("\n");
                }
                printf("density of matrix: %.2f%%\n", density);
                for (len_t k = 0; k < mat->nr; ++k) {
                    free(mat->tr[k]);
                }
                free(mat);
                mat = NULL;
                free(hcm);
                hcm = NULL;
                if (sht != NULL) {
                    free_hash_table(&sht);
                }
            }
            /* free and clean up */
            if (bs != NULL) {
                free_basis(&bs);
            }
            if (tbr != NULL) {
                free_basis(&tbr);
            }
            free(st);
            st  = NULL;
            free_shared_hash_data(bht);
            if (bht != NULL) {
                free_hash_table(&bht);
            }

        }
        return 0;
    }
    else{/* characteristic is 0 */
#if 0
      if (elim_block_len > 0) { /* characteristic is 0 */
            /* timings */
            double ct0, ct1, rt0, rt1;
            ct0 = cputime();
            rt0 = realtime();
            uint32_t field_char         = gens->field_char;
            const uint32_t prime_start = pow(2, 30);
            const int32_t nr_primes = nr_threads;
            /* initialize stuff */
            md_t *st  = allocate_meta_data();

            int *invalid_gens       =   NULL;
            int32_t monomial_order  =   0;
            int32_t reduce_gb       =   1;
            int res = validate_input_data(&invalid_gens, gens->mpz_cfs,
                    gens->lens, &field_char, &monomial_order, &elim_block_len,
                    &gens->nvars, &gens->ngens, &saturate, &initial_hts,
                    &nr_threads, &max_pairs, &update_ht, &la_option,
                    &use_signatures, &reduce_gb, &info_level);
	    
            /* all data is corrupt */
            if (res == -1) {
                fprintf(stderr, "Invalid input generators, msolve now terminates.\n");
                free(invalid_gens);
                return -3;
            }

            /* checks and set all meta data. if a nonzero value is returned then
             * some of the input data is corrupted. */
            if (check_and_set_meta_data_trace(st, gens->lens, gens->exps,
                        (void *)gens->mpz_cfs, invalid_gens, gens->field_char, 0,
                        elim_block_len, gens->nvars, gens->ngens, saturate,
                        initial_hts, nr_threads, max_pairs, update_ht,
                        la_option, use_signatures, 1, prime_start,
                        nr_primes, 0, info_level)) {
                free(st);
                return -3;
            }

            /* lucky primes */
            primes_t *lp  = (primes_t *)calloc(1, sizeof(primes_t));

            /*******************
             * initialize basis
             *******************/
            bs_t *bs_qq = initialize_basis(st);
            /* initialize basis hash table, update hash table, symbolic hash table */
            ht_t *bht = initialize_basis_hash_table(st);
            /* hash table to store the hashes of the multiples of
             * the basis elements stored in the trace */
            ht_t *tht = initialize_secondary_hash_table(bht, st);
            /* read in ideal, move coefficients to integers */
            import_input_data(bs_qq, bht, st, gens->lens, gens->exps,
                    (void *)gens->mpz_cfs, invalid_gens);

            print_initial_statistics(stderr, st);

            /* for faster divisibility checks, needs to be done after we have
             * read some input data for applying heuristics */
            calculate_divmask(bht);

            /* sort initial elements, smallest lead term first */
            sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *),
                    initial_input_cmp, bht);
            remove_content_of_initial_basis(bs_qq);


            /* generate lucky prime numbers */
            generate_lucky_primes(lp, bs_qq, st->prime_start, st->nprimes);

            /* generate array to store modular bases */
            bs_t **bs = (bs_t **)calloc((unsigned long)st->nprimes, sizeof(bs_t *));

            /* initialize tracer */
            trace_t *trace  = initialize_trace(bs_qq, st);

            srand(time(0));
            uint32_t prime = next_prime(1<<30);
            prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
            while(is_lucky_prime_ui(prime, bs_qq)){
                prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
            }

            uint32_t primeinit = prime;
            lp->p[0] = primeinit;


            if (is_gb == 1) {
                for (len_t k = 0; k < bs_qq->ld; ++k) {
                    bs_qq->lmps[k] = k;
                    bs_qq->lm[k]   = bht->hd[bs_qq->hm[k][OFFSET]].sdm;
                    bs_qq->lml     = bs_qq->ld;
                }
            }

            st->laopt = 44;
            bs_t *bsprob = modular_f4(bs_qq, bht, st, lp->p[0]);


            st->laopt = 2;
            /* compute a gb for initial generators */
            gba_trace_learning_phase(
                    trace,
                    tht,
                    bs_qq,
                    bht,
                    st,
                    lp->p[0]);

            /* int64_t nb  = export_results_from_gba(bld, blen, bexp,
             *         bcf, &bs, &bht, &st); */

            /* timings */
            ct1 = cputime();
            rt1 = realtime();
            st->f4_ctime = ct1 - ct0;
            st->f4_rtime = rt1 - rt0;

            if (st->info_level > 1) {
                print_final_statistics(stderr, st);
            }
            if(info_level){
                fprintf(stderr, "\nStarts trace based multi-modular computations\n");
            }

            int i;

            ht_t *lht = copy_hash_table(bht, st);

            prime = next_prime(1<<30);

            /* while(rerun == 1 || mcheck == 1){ */

                /* generate lucky prime numbers */
                prime = next_prime(prime);
                lp->p[0] = prime;
                while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
                    prime = next_prime(prime);
                    lp->p[0] = prime;
                }

                for(len_t i = 1; i < st->nprimes; i++){
                    prime = next_prime(prime);
                    lp->p[i] = prime;
                    while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
                        prime = next_prime(prime);
                        lp->p[i] = prime;
                    }
                }
                prime = lp->p[st->nprimes - 1];

                double ca0;

                double stf4 = 0;

                for (i = 0; i < st->nprimes; ++i){
                    ca0 = realtime();
                    bs[i] = gba_trace_application_phase(
                            trace,
                            tht,
                            bs_qq,
                            lht,
                            st,
                            lp->p[i]);

                    stf4 = realtime()-ca0;
                    printf("F4 trace timing %13.2f\n", stf4);
                    /* printf("bs[%u]->lml = %u\n", i, bs[i]->lml); */
                }
            /* } */
            return 0;
        }
#endif
      /* characteristic is 0 */
      /* characteristic is 0 and elim_block = 0 */
      if (saturate == 1) {       /* characteristic is 0 and elim_block = 0 */
            /* timings */
            double ct0, ct1, rt0, rt1;
            ct0 = cputime();
            rt0 = realtime();
            uint32_t field_char         = gens->field_char;
            const uint32_t prime_start  = pow(2, 30);
            const int32_t nr_primes     = nr_threads;

            /* data structures for basis, hash table and statistics */
            bs_t *sat_qq   = NULL;

            /* initialize stuff */
            md_t *st  = allocate_meta_data();

            int *invalid_gens       =   NULL;
            int32_t monomial_order  =   0;
            int32_t reduce_gb       =   1;
            int res = validate_input_data(&invalid_gens, gens->mpz_cfs,
                    gens->lens, &field_char, &monomial_order, &elim_block_len,
                    &gens->nvars, &gens->ngens, &saturate, &initial_hts,
                    &nr_threads, &max_pairs, &update_ht, &la_option,
                    &use_signatures, &reduce_gb, &info_level);

            /* all data is corrupt */
            if (res == -1) {
                fprintf(stderr, "Invalid input generators, msolve now terminates.\n");
                free(invalid_gens);
                return -3;
            }

            /* checks and set all meta data. if a nonzero value is returned then
             * some of the input data is corrupted. */
            if (check_and_set_meta_data_trace(st, gens->lens, gens->exps,
                        (void *)gens->mpz_cfs, invalid_gens,
                        field_char, 0, elim_block_len, gens->nvars,
                        gens->ngens, saturate, initial_hts, nr_threads,
                        max_pairs, update_ht, la_option, use_signatures,
                        1, prime_start, nr_primes, 0, info_level)) {
                free(st);
                return -3;
            }

            /* lucky primes */
            primes_t *lp  = (primes_t *)calloc(1, sizeof(primes_t));

            /*******************
             * initialize basis
             *******************/
            bs_t *bs_qq = initialize_basis(st);
            /* initialize basis hash table, update hash table, symbolic hash table */
            ht_t *bht = bs_qq->ht;
            /* hash table to store the hashes of the multiples of
             * the basis elements stored in the trace */
            ht_t *tht = initialize_secondary_hash_table(bht, st);
            /* read in ideal, move coefficients to integers */
            import_input_data(bs_qq, st, gens->lens, gens->exps,
                    (void *)gens->mpz_cfs, invalid_gens);
            free(invalid_gens);
            invalid_gens    =   NULL;

            print_initial_statistics(stderr, st);

            /* for faster divisibility checks, needs to be done after we have
             * read some input data for applying heuristics */
            calculate_divmask(bht);

            /* sort initial elements, smallest lead term first */
            sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *),
                    initial_input_cmp, bht);
            remove_content_of_initial_basis(bs_qq);


            /* generate lucky prime numbers */
            generate_lucky_primes(lp, bs_qq, st->prime_start, st->nprimes);

            /* generate array to store modular bases */
            bs_t **bs = (bs_t **)calloc((unsigned long)st->nprimes, sizeof(bs_t *));

            /* initialize tracer */
            trace_t *trace  = initialize_trace(bs_qq, st);

            st->tr = trace;

            srand(time(0));
            uint32_t prime = next_prime(1<<30);
            prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
            while(is_lucky_prime_ui(prime, bs_qq)){
                prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
            }

            uint32_t primeinit = prime;
            lp->p[0] = primeinit;


            if (is_gb == 1) {
                for (len_t k = 0; k < bs_qq->ld; ++k) {
                    bs_qq->lmps[k] = k;
                    bs_qq->lm[k]   = bht->hd[bs_qq->hm[k][OFFSET]].sdm;
                    bs_qq->lml     = bs_qq->ld;
                }
            }
            sat_qq = initialize_basis(st);
            import_input_data_nf_qq(
                    sat_qq, bht, st, gens->ngens-saturate, gens->ngens,
                    gens->lens, gens->exps, (void *)gens->mpz_cfs);
            sat_qq->ld = sat_qq->lml  =  saturate;
            /* normalize_initial_basis(tbr, st->gfc); */
            for (int k = 0; k < saturate; ++k) {
                sat_qq->lmps[k]  = k; /* fix input element in tbr */
            }

            printf("LEARNING PHASE -- PART 1\n");
            f4sat_trace_learning_phase_1(
                    trace,
                    tht,
                    bs_qq,
                    sat_qq,
                    &bht,
                    st,
                    lp->p[0]);

            /* int64_t nb  = export_results_from_gba(bld, blen, bexp,
             *         bcf, &bs, &bht, &st); */

            /* timings */
            ct1 = cputime();
            rt1 = realtime();
            st->f4_ctime = ct1 - ct0;
            st->f4_rtime = rt1 - rt0;

            get_and_print_final_statistics(stderr, st, bs_qq);

            if(info_level){
                fprintf(stderr, "\nStarts trace based multi-modular computations\n");
            }

            prime = next_prime(1<<30);

            lp->p[0]  = prime;

            printf("LEARNING PHASE -- PART 2\n");
            f4sat_trace_learning_phase_2(
                    trace,
                    tht,
                    bs_qq,
                    sat_qq,
                    &bht,
                    st,
                    lp->p[0]);

            /* int64_t nb  = export_results_from_gba(bld, blen, bexp,
             *         bcf, &bs, &bht, &st); */

            /* timings */
            ct1 = cputime();
            rt1 = realtime();
            st->f4_ctime = ct1 - ct0;
            st->f4_rtime = rt1 - rt0;

            get_and_print_final_statistics(stderr, st, bs_qq);

            if(info_level){
                fprintf(stderr, "\nStarts trace based multi-modular computations\n");
            }

            int i;

            ht_t *lht = copy_hash_table(bht, st);

            prime = next_prime(1<<30);

            /* while(rerun == 1 || mcheck == 1){ */

                /* generate lucky prime numbers */
                prime = next_prime(prime);
                lp->p[0] = prime;
                while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
                    prime = next_prime(prime);
                    lp->p[0] = prime;
                }

                for(len_t i = 1; i < st->nprimes; i++){
                    prime = next_prime(prime);
                    lp->p[i] = prime;
                    while(is_lucky_prime_ui(prime, bs_qq) || prime==primeinit){
                        prime = next_prime(prime);
                        lp->p[i] = prime;
                    }
                }
                prime = lp->p[st->nprimes - 1];

                double ca0;

                double stf4 = 0;

                for (i = 0; i < st->nprimes; ++i){
                    ca0 = realtime();
                    bs[i] = f4sat_trace_application_phase(
                            trace,
                            tht,
                            bs_qq,
                            sat_qq,
                            lht,
                            st,
                            lp->p[i]);

                    stf4 = realtime()-ca0;
                    printf("F4 trace timing %13.2f\n", stf4);
                    /* printf("bs[%u]->lml = %u\n", i, bs[i]->lml); */
                }
            /* } */
            return 0;
      } else {       /* characteristic is 0 and elim_block = 0 and saturate = 0 */

            int dim = - 2;
            long dquot = -1;

            if(elim_block_len && print_gb == 0){
              if(info_level){
                fprintf(stderr, "Warning: elim order not available for rational parametrizations\n");
                fprintf(stderr, "Computing Groebner basis\n");
                print_gb=2;
              }
            }

            if(print_gb){

              msflags_t flags;

              flags->ht_size = initial_hts;
              flags->nr_threads = nr_threads;
              flags->max_nr_pairs = max_pairs;
              flags->elim_block_len = elim_block_len;
              flags->reset_ht = update_ht;
              flags->la_option = la_option;
              flags->use_signatures = use_signatures;
              flags->info_level = info_level;
              flags->pbm_file = generate_pbm;
              flags->print_gb = print_gb;
              flags->files = files;

              print_msolve_gbtrace_qq(gens, flags);
              return 0;
            }

            b = real_msolve_qq(*mpz_paramp,
                    &param,
                    &dim,
                    &dquot,
                    nb_real_roots_ptr,
                    real_roots_ptr,
                    real_pts_ptr,
                    gens,
                    initial_hts, nr_threads, max_pairs,
                    elim_block_len, update_ht,
                    la_option, use_signatures, info_level, print_gb,
                    generate_pbm, precision, files, round, get_param);

            if(print_gb){
              return 0;
            }

            manage_output(b, dim, dquot, files, gens, param, mpz_paramp, get_param,
                          nb_real_roots_ptr,
                          real_roots_ptr,
                          real_pts_ptr,
                          info_level);
            if (b == 1) {
                free(bld);
                bld = NULL;
                free(blen);
                blen  = NULL;
                free(bexp);
                bexp  = NULL;
                free(bcf);
                bcf = NULL;
                free(param);
                param = NULL;
                if (genericity_handling > 0) {
                    if (change_variable_order_in_input_system(gens, info_level)) {
                        goto restart;
                    }
                    if (genericity_handling == 2) {
                        if (add_linear_form_to_input_system(gens, info_level)) {
                            goto restart;
                        }
                    }
                }
                fprintf(stderr, "\n=====> Computation failed <=====\n");
                fprintf(stderr, "Try to add a random linear form with ");
                fprintf(stderr, "a new variable\n");
                fprintf(stderr, "(smallest w.r.t. DRL) to the input system. ");
                fprintf(stderr, "This will\n");
                fprintf(stderr, "be done automatically if you run msolve with option\n");
                fprintf(stderr, "\"-c2\" which is the default.\n");
                (*mpz_paramp)->dim  = -1;
            }
            if(b == -4){
                fprintf(stderr, "Bad prime chosen initially\n");
                goto restart;
                /* (*mpz_paramp)->dim  = -4; */
            }
            if(b == 2){
                free(bld);
                bld = NULL;
                free(blen);
                blen  = NULL;
                free(bexp);
                bexp  = NULL;
                free(bcf);
                bcf = NULL;
                free(param);
                param = NULL;
                round++;
                if(gens->change_var_order >= 0){
                  undo_variable_order_change(gens);
                }
                if (add_random_linear_form_to_input_system(gens, info_level)) {
                    goto restart;
                }
            }
        }

        free(bld);
        free(blen);
        free(bexp);
        free(bcf);

        /* get parametrization */
        *paramp     = param;

        return !(b==0);
    }
}

static void export_julia_rational_parametrization_qq(
        void *(*mallocp) (size_t),
        int32_t *load,
        int32_t *nvars,
        int32_t *dim,
        int32_t *dim_quot,
        int32_t **lens,
        char ***var_namesp,
        void **cfs_linear_form,
        void **cfs,
        void **real_sols_num,
        int32_t **real_sols_den,
        data_gens_ff_t *gens, /* might change vnames, thus not const */
        const mpz_param_t param,
        const long nb_real_roots,
        const real_point_t *real_pts
        )
{
    int32_t i, j;
    int64_t ctr = 0;

    *load     = (int32_t)param->nvars+1;
    *dim      = (int32_t)param->dim;
    *dim_quot = (int32_t)param->dquot;
    *nvars    = (int32_t)gens->nvars;

    /* keep variable names for returning rational parametrization */
    *var_namesp  = gens->vnames;
    gens->vnames = NULL;

    mpz_t *cf_lf = NULL;
    /* check existence of linear form */
    if (gens->linear_form_base_coef > 0) {
        cf_lf = (mpz_t *)(*mallocp)(
                (unsigned long)(gens->nvars) * sizeof(mpz_t));
        int64_t len = 0;
        for (i = 0; i < gens->ngens-1; ++i) {
            len += 2*gens->lens[i]; /* numerator plus denominator, thus "2*" */
        }
        j = 0;
        for (i = 0; i < 2*gens->nvars; i += 2) {
            mpz_init_set(cf_lf[j++], *(gens->mpz_cfs[i+len]));
        }
    }

    if ((param->dim > 0) || (param->dim == 0 && param->dquot == 0)) {
        *lens = NULL;
        *cfs  = NULL;
    } else {
        int32_t *len  = (int32_t *)(*mallocp)(
                (unsigned long)(param->nvars+1) * sizeof(int32_t));

        /* precompute number of all terms of all polynomials */
        int64_t nterms  = 0;
        nterms  += (int64_t)param->elim->length;
        len[0]  =  (int32_t)param->elim->length;
        nterms  += (int64_t)param->denom->length;
        len[1]  =  (int32_t)param->denom->length;
        /* In the following we add one more term since each coords[i] polynomial
         * also gets another coefficient for the correct denominator in the
         * rational parametrization of the solution set. */
        for (i = 0; i < param->nvars-1; ++i) {
            nterms    +=  param->coords[i]->length+1;
            len[i+2]  =   param->coords[i]->length+1;
        }

        mpz_t *cf     = (mpz_t *)(*mallocp)(
                (unsigned long)(nterms) * sizeof(mpz_t));

        /* store elim */
        for (i = 0; i < param->elim->length; ++i) {
            mpz_init_set((cf+ctr)[i], param->elim->coeffs[i]);
        }
        ctr +=  param->elim->length;
        /* store denom */
        for (i = 0; i < param->denom->length; ++i) {
            mpz_init_set((cf+ctr)[i], param->denom->coeffs[i]);
        }
        ctr +=  param->denom->length;

        /* store param */
        for (i = 0; i < param->nvars-1; ++i) {
            for (j = 0; j < param->coords[i]->length; ++j) {
                mpz_init_set((cf+ctr)[j], param->coords[i]->coeffs[j]);
            }
            mpz_init_set((cf+ctr)[j], param->cfs[i]);
            ctr +=  param->coords[i]->length+1;
        }
        *lens            = len;
        *cfs             = (void *)cf;
        *cfs_linear_form = (void *)cf_lf;

        /* if there are no real solutions return the parametrization at least */
        if (nb_real_roots <= 0) {
            return;
        }

        const long nb_real_roots_intervall  = 2 * nb_real_roots;

        mpz_t *sols_num = (mpz_t *)(*mallocp)(
                (unsigned long)nb_real_roots_intervall * real_pts[0]->nvars * sizeof(mpz_t));

        int32_t *sols_den = (int32_t *)(*mallocp)(
                (unsigned long)nb_real_roots_intervall * real_pts[0]->nvars * sizeof(int32_t));


        /* mpz_t tmp;
         * mpz_init(tmp); */

        ctr = 0;
        for (i = 0; i < nb_real_roots; ++i) {
            for (j = 0; j < real_pts[i]->nvars; ++j) {
                /* mpz_add(tmp, real_pts[i]->coords[j]->val_do,
                 *         real_pts[i]->coords[j]->val_up); */
                mpz_init_set(sols_num[ctr], real_pts[i]->coords[j]->val_do);
                sols_den[ctr++] = real_pts[i]->coords[j]->k_do;
                mpz_init_set(sols_num[ctr], real_pts[i]->coords[j]->val_up);
                sols_den[ctr++] = real_pts[i]->coords[j]->k_up;
            }
        }

        *real_sols_num = (void *)sols_num;
        *real_sols_den = sols_den;
    }
}

void msolve_julia(
        void *(*mallocp) (size_t),
        int32_t *rp_ld,
        int32_t *rp_nr_vars,
        int32_t *rp_dim,
        int32_t *rp_dquot,
        int32_t **rp_lens,
        char ***rp_var_namesp,
        void **rp_cfs_linear_form,
        void **rp_cfs,
        int32_t *n_real_sols,
        void **real_sols_num,
        int32_t **real_sols_den,
        int32_t *lens,
        int32_t *exps,
        void *cfs,
        char **var_names,
        char *output_file,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t initial_hts,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t print_gb,
        const int32_t get_param,
        const int32_t genericity_handling,
        const int32_t precision,
        const int32_t info_level
        )
{
    /* timinigs */
    double st0 = cputime();
    double rt0 = realtime();

    len_t i;
    files_gb *files = calloc(1, sizeof(files_gb));

    if (output_file != NULL) {
        files->out_file = output_file;
    }

    data_gens_ff_t *gens = allocate_data_gens();

    unsigned long nterms  = 0;
    for (i = 0; i < nr_gens; ++i) {
        nterms  +=  lens[i];
    }
    gens->nvars                 = nr_vars;
    gens->ngens                 = nr_gens;
    gens->field_char            = field_char;
    gens->change_var_order      = -1;
    gens->linear_form_base_coef = 0;
    /* gens->vnames                = var_names; */
    gens->vnames  = (char **)malloc((unsigned long)nr_vars * sizeof(char *));
    for (i = 0; i < nr_vars; ++i) {
        gens->vnames[i] = calloc((unsigned long)strlen(var_names[i]), sizeof(char));
        memcpy(gens->vnames[i], var_names[i], (unsigned long)strlen(var_names[i]) * sizeof(char));
    }
    /* gens->lens                  = lens; */
    gens->lens  = (int32_t *)malloc((unsigned long)nr_gens * sizeof(int32_t));
    memcpy(gens->lens, lens, (unsigned long)nr_gens * sizeof(int32_t));
    /* gens->exps                  = exps; */
    gens->exps  = (int32_t *)malloc(nterms * nr_vars * sizeof(int32_t));
    memcpy(gens->exps, exps, nterms * nr_vars * sizeof(int32_t));
    gens->rand_linear           = 0;

    if (field_char > 0) {
        gens->cfs = (int32_t *)malloc(nterms * sizeof(int32_t));
        memcpy(gens->cfs, (int32_t *)cfs, nterms * sizeof(int32_t));
        /* gens->cfs = (int32_t *)cfs; */
    } else {
        if (field_char == 0) {
            gens->mpz_cfs = (mpz_t **)malloc(nterms * 2 * sizeof(mpz_t *));
            for (i = 0; i < 2*nterms; ++i) {
                gens->mpz_cfs[i]  = (mpz_t *)malloc(sizeof(mpz_t));
                mpz_init_set(*(gens->mpz_cfs[i]), *(((mpz_t **)cfs)[i]));
            }
        }
    }

    /* data structures for parametrization */
    param_t *param  = NULL;
    mpz_param_t mpz_param;
    mpz_param_init(mpz_param);

    long nb_real_roots      = 0;
    interval *real_roots    = NULL;
    real_point_t *real_pts  = NULL;

    /* main msolve functionality */
    int ret = core_msolve(la_option, use_signatures, nr_threads, info_level,
            initial_hts, max_nr_pairs, elim_block_len, reset_ht,
            0 /* generate pbm */, 1 /* reduce_gb */, print_gb, get_param,
            genericity_handling, 0 /* saturate */, 0 /* colon */,
	    0 /* normal_form */, 0 /* normal_form_matrix */,
	    0 /* is_gb */, precision, files,
            gens, &param, &mpz_param, &nb_real_roots, &real_roots, &real_pts);

    if (ret == -1) {
        exit(1);
    }
    *rp_dim =   mpz_param->dim;

    char **rp_var_names = NULL;
    if (mpz_param->dim != -1) {
        export_julia_rational_parametrization_qq(
                mallocp, rp_ld, rp_nr_vars, rp_dim, rp_dquot, rp_lens,
                &rp_var_names, rp_cfs_linear_form, rp_cfs, real_sols_num,
                real_sols_den, gens, mpz_param, nb_real_roots, real_pts);
    } else {
        *rp_ld  = -1;
    }

    /* clean up data storage, but do not free data handled by julia */

    free(gens);
    gens  = NULL;

    *rp_var_namesp = rp_var_names;

    /* free parametrization */
    free(param);
    mpz_param_clear(mpz_param);

    *n_real_sols = nb_real_roots;

    free(real_roots);

    if (nb_real_roots > 0) {
        for(long i = 0; i < nb_real_roots; i++){
          real_point_clear(real_pts[i]);
        }
        free(real_pts);
    }

    /* timings */
    if (info_level > 0) {
        double st1 = cputime();
        double rt1 = realtime();
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
        fprintf(stderr, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
    }
}

/* The parameters themselves are handled by julia, thus we only
 * free what they are pointing to, julia's garbage collector then
 * takes care of everything leftover. */
void free_msolve_julia_result_data(
        void (*freep) (void *),
        int32_t **res_len,
        void **res_cf,
        void **sols_num,
        int32_t **sols_den,
        const int64_t res_ld,
        const int64_t nr_sols,
        const int64_t field_char
        )
{


    int32_t *lens  = *res_len;

    (*freep)(lens);
    lens      = NULL;
    *res_len  = lens;

    if (field_char > 0) {
        int32_t *numerators = *(int32_t **)sols_num;
        (*freep)(numerators);
        numerators  = NULL;
        int32_t *cfs= *(int32_t **)res_cf;
        (*freep)(cfs);
        cfs = NULL;
    } else {
        /* denominators only exist if char == 0 */
        int32_t *denominators = *sols_den;
        (*freep)(denominators);
        denominators  = NULL;
        *sols_den     = denominators;
        /* mpz_t **numerators  = (mpz_t **)sols_num;
         * for (i = 0; i < nr_sols; ++i) {
         *     (*freep)((*numerators)[i]);
         * }
         * (*freep)(*numerators);
         * *numerators = NULL;
         * mpz_t **cfs = (mpz_t **)res_cf;
         * for (i = 0; i < len; ++i) {
         *     (*freep)((*cfs)[i]);
         * }
         * (*freep)(*cfs);
         * *cfs  = NULL; */
    }
    *sols_num = NULL;
    *res_cf   = NULL;
}
