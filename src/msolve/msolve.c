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

#define MAX(a,b) (((a)>(b))?(a):(b))
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define ilog2_mpz(a) mpz_sizeinbase(a,2)

static void mpz_upoly_init(mpz_upoly_t poly, long length){
  mpz_t *tmp = NULL;
  if(length){
    tmp = (mpz_t *)(malloc(length * sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < length; i++){
      mpz_init(tmp[i]);
      mpz_set_ui(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = length;
  poly->length = -1;
}

static void mpz_upoly_init2(mpz_upoly_t poly, long length, long nbits){
  mpz_t *tmp = NULL;
  if(length){
    tmp = (mpz_t *)(malloc(length * sizeof(mpz_t)));
    //    tmp = (mpz_t *)(calloc(length, sizeof(mpz_t)));
    if(tmp==NULL){
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for(long i = 0; i < length; i++){
      mpz_init2(tmp[i], nbits);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = length;
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
  for(long i = 0; i < param->nvars - 1; i++){
    mpz_upoly_clear(param->coords[i]);
    mpz_clear(param->cfs[i]);
  }
  free(param->coords);
  free(param->cfs);
  param->nvars  = 0;
  param->nsols  = 0;
}

static inline void mpz_param_out_str(FILE *file, const data_gens_ff_t *gens,
        const long dquot, mpz_param_t param){
  fprintf(file, "[");
  fprintf(file, "0, \n"); //dimension of input ideal
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
  }
  fprintf(file, "],\n");
  mpz_upoly_out_str(file, param->elim); //elim. poly
  fprintf(file, ",\n");
  mpz_upoly_out_str(file, param->denom); //denom. poly
  fprintf(file, ",\n");
  fprintf(file, "[\n");
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
  /* fprintf(file, "]"); */
  fprintf(file, "]");
}

static inline void mpz_param_out_str_maple(FILE *file,
        const data_gens_ff_t *gens,const long dquot, mpz_param_t param){
  mpz_param_out_str(file, gens, dquot, param);
  fprintf(file, "],");
}

static inline data_gens_ff_t *allocate_data_gens(){
  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
  gens->lens  = NULL;
  gens->exps  = NULL;
  gens->cfs   = NULL;
  gens->mpz_cfs = NULL;

  return gens;
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
        mpz_set_ui(*(gens->mpz_cfs[2*(len_new - 1)]), 1); 
    }
    return 1;
}

static int add_random_linear_form_to_input_system(
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
    /* const int32_t bcf = gens->linear_form_base_coef; */
    k = 1;
    if (info_level > 0) {
        printf("\nAdding a linear form with an extra variable ");
        printf("(lowest w.r.t. monomial order)\n");
        printf("[coefficients of linear form are randomly chosen]\n");
    }
    srand(time(0));
    gens->random_linear_form = malloc(sizeof(int32_t *)*(nvars_new));
    if (gens->field_char > 0) {
      int j = 0;
      for (i = len_old; i < len_new; ++i) {
        gens->random_linear_form[j] = ((int16_t)(rand()) % gens->field_char);
        gens->cfs[i]  = gens->random_linear_form[j];
        k++;
        j++;
      }
    }
    else {
      int j = 0;

      for (i = 2*len_old; i < 2*len_new; i += 2) {
        gens->random_linear_form[j] = ((int16_t)(rand()));
        mpz_set_ui(*(gens->mpz_cfs[i]), gens->random_linear_form[j]);
        k++;
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
  mpz_upoly_init2(param->elim, bparam->elim->length, 2*32*(bparam->elim->length));
  mpz_upoly_init(param->denom, bparam->elim->length - 1);
  param->elim->length = bparam->elim->length;

  param->coords = (mpz_upoly_t *)malloc(sizeof(mpz_upoly_t)*(param->nvars - 1));
  if(param->coords != NULL){
    for(long i = 0; i < param->nvars - 1; i++){
      mpz_upoly_init(param->coords[i], bparam->elim->length - 1);
      param->coords[i]->length = bparam->coords[i]->length;
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
  uint32_t* lineqs = calloc(nlins*(nvars + 1), sizeof(uint64_t));
  for(long i = 0; i < nlins*(nvars+1); i++){
    lineqs[i] = 0;
  }

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



/**

On fait un calcul modulaire.
donc ici on suppose gens->field_char >0

Si la dimension est 0 et que l'on est en position generique
on renvoie 0 (et une param.).

Si la dimension est > 0, on renvoie un nombre positif (ce sera la dimension a
terme) et on a la BG.

Si il y a eu un pb (pas de position generique ou echec quelconque) on renvoie -1.
 **/

int msolve_ff(param_t **bparam,
              //int32_t *bld, int32_t **blen, int32_t **bexp, void **bcf,
              data_gens_ff_t *gens,
              int32_t initial_hts,
              int32_t nr_threads,
              int32_t max_pairs,
              int32_t update_ht,
              int32_t la_option,
              int32_t info_level,
              files_gb *files){

  int32_t *bld = malloc(sizeof(int32_t)*gens->ngens);
  int32_t **blen = malloc(sizeof(int32_t *));
  int32_t **bexp = malloc(sizeof(int32_t *));
  void **bcf = malloc(sizeof(void *));


  if(info_level > 0){
    fprintf(stderr, "Starts F4 with prime = %d\n", gens->field_char);
  }
  long nlins = 0;
  uint64_t *linvars = calloc(gens->nvars, sizeof(uint64_t));
  uint32_t **lineqs_ptr = malloc(sizeof(uint32_t *));

  int64_t nb = f4_julia(bld, blen, bexp, bcf,
                        gens->lens, gens->exps, (void *)gens->cfs, gens->field_char,
                        0, //mon_order,
                        gens->nvars, gens->ngens, initial_hts,
                        nr_threads, max_pairs, update_ht, la_option,
                        1, //reduce_gb
                        0, //generate_pbm
                        info_level);

  if(nb==0){
    fprintf(stderr, "Something went wrong during the computation\n");
    return -1;
  }

  int32_t *bcf_ff = (int32_t *)(*bcf);

#if DEBUGGB > 0
  if(files->out_file != NULL){
    char fn[200];
    char spr[100];
    sprintf(spr, "%d", gens->field_char);
    strcpy(fn, "/tmp/sys.basis.");
    strcat(fn,spr);
    FILE *ofilebs = fopen(fn, "w");
    display_basis_maple(ofilebs, bld, blen, bexp, bcf_ff, gens);
    fclose(ofilebs);
  }
  else{
    display_basis_maple(stdout, bld, blen, bexp, bcf_ff, gens);
  }
#endif

#if DEBUGGB>0
  if(files->out_file != NULL){
    char fn[200];
    char spr[100];
    sprintf(spr, "%d", gens->field_char);
    strcpy(fn, "/tmp/sys_lead_monomial.");
    strcat(fn,spr);
    FILE *ofilelm = fopen(fn, "w");
    display_lead_monomials_from_gb(ofilelm, bld, blen, bexp, bcf_ff, gens);
    fclose(ofilelm);
  }
  else{
    display_lead_monomials_from_gb(stdout, bld, blen, bexp, bcf_ff, gens);
    fprintf(stdout, "\n");
  }
#endif

  double st_lm = realtime();
  int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);
  check_and_set_linear_poly_non_hashed(&nlins, linvars, lineqs_ptr,
                                       bld, bexp_lm,
                                       blen, bexp, bcf_ff, gens->nvars);

  if(has_dimension_zero(bld[0], gens->nvars, bexp_lm)){
    if(info_level){
      fprintf(stderr, "The ideal has dimension zero\n");
    }
    long dquot = 0;
    int32_t *lmb = monomial_basis(bld[0], gens->nvars, bexp_lm, &dquot);

    if(info_level){
      fprintf(stderr, "Dimension of quotient ring = %ld\n", dquot);
    }
#if DEBUGGB>0
    char fn[200];
    char spr[100];
    sprintf(spr, "%d", gens->field_char);
    strcpy(fn, "/tmp/sys.mon_basis.");
    strcat(fn,spr);
    FILE *ofile = fopen(fn, "w");
    display_monomials_from_array_maple(ofile, dquot, lmb, gens);
    fclose(ofile);
#endif

    if(dquot==1){
#if DEBUGGB>0
      if(files->out_file!=NULL){
        FILE *ofile = fopen(files->out_file, "w");
        display_basis_maple(ofile, bld, blen, bexp, bcf_ff, gens);
        fclose(ofile);
       }
      else{
        display_basis_maple(stdout, bld, blen, bexp, bcf_ff, gens);
      }
#endif
      free(bld);
      free(*blen);
      free(blen);
      free(*bexp);
      free(bexp);
      free(*bcf);
      free(bcf);
      free(linvars);
      if(nlins){
        free(lineqs_ptr[0]);
      }
      free(lineqs_ptr);
      return 0;
    }
    if(info_level > 1){
      fprintf(stderr, "Build monomial basis: %.2f sec.\n", realtime()-st_lm);
    }
    sp_matfglm_t *matrix = build_matrixn(lmb, dquot, bld[0], blen, bexp, bcf_ff,
                                         bexp_lm, gens->nvars, gens->field_char);

    if(matrix==NULL){
      /* fprintf(stderr, "Problem when building the matrix multiplication\n");
       * fprintf(stderr, "Some data should be free-ed\n"); */
      free(bld);
      free(*blen);
      free(blen);
      free(*bexp);
      free(bexp);
      free(*bcf);
      free(bcf);
      free(bexp_lm);
      free(lmb);
      free(linvars);
      if(nlins){
        free(lineqs_ptr[0]);
      }
      free(lineqs_ptr);
      if (dquot > 0) {
          return 1;
      } else {
        return -1;
      }
    }
#if DEBUGGB > 0
    display_fglm_matrix(stdout, matrix);
#endif
    //display_fglm_matrix(stderr, matrix);

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    if(info_level > 1){
      fprintf(stderr, "Starts FGLM\n");
    }
    uint64_t *squvars = calloc(gens->nvars-1, sizeof(uint64_t));
    check_and_set_vars_squared_in_monomial_basis(squvars, lmb,
						 dquot, gens->nvars);
    *bparam = nmod_fglm_compute(matrix, gens->field_char,
                                gens->nvars, nlins, linvars, lineqs_ptr[0],
                                squvars, info_level);

#if DEBUGGB>0
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "w");
      display_fglm_param_maple(ofile, *bparam);
      fclose(ofile);
    }
    else{
      display_fglm_param_maple(stdout, *bparam);
    }
#endif

    free_sp_mat_fglm(matrix);
    free(bexp_lm);
    free(lmb);

    free(bld);
    free(*blen);
    free(blen);
    free(*bexp);
    free(bexp);
    free(*bcf);
    free(bcf);
    free(linvars);
    if(nlins){
      free(lineqs_ptr[0]);
    }
    free(lineqs_ptr);
    free(squvars);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    if (info_level) {
      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
      fprintf(stderr, "FGLM TIMING %13.2f sec (REAL) / %5.2f sec (CPU)\n",
            rt1-rt0, ct1-ct0);
      fprintf(stderr, "-------------------------------------------------\
----------------------------------------\n");
      fprintf(stderr, "Done (matrix is free-ed).\n");
    }

    return 0;
  }
  else{
    fprintf(stderr, "The ideal is not zero-dimensional\n");

    free(bld);
    free(*blen);
    free(blen);
    free(*bexp);
    free(bexp);
    free(*bcf);
    free(bcf);
    free(linvars);
    if(nlins){
      free(lineqs_ptr[0]);
    }
    free(lineqs_ptr);

    return 2;
#if DEBUGGB>0
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "w");
      display_basis_maple(ofile, bld, blen, bexp, bcf_ff, gens);
      fclose(ofile);
    }
    else{
      display_basis_maple(stdout, bld, blen, bexp, bcf_ff, gens);
    }
#endif
  }

  free(bld);
  free(*blen);
  free(blen);
  free(*bexp);
  free(bexp);
  free(*bcf);
  free(bcf);

  return -1;
}



/**

Modular computation:
does as msolve_ff but with bld, blen, etc already allocated.

returns 0 when the ideal has dim. 0 and coordinates are generic enough to obtain the parametrization. 
returns 1 in case of failure.
returns 2 when the dimension is > 0

Small hack here: returns -1 when the dimension of the quotient is 1. 
**/

int msolve_ff_alloc(param_t **bparam,
                    int32_t *bld, int32_t **blen, int32_t **bexp, void **bcf,
                    data_gens_ff_t *gens,
                    int32_t initial_hts,
                    int32_t nr_threads,
                    int32_t max_pairs,
                    int32_t update_ht,
                    int32_t la_option,
                    int32_t info_level,
                    int32_t print_gb,
                    files_gb *files){


    if(info_level){
        fprintf(stderr, "Starts F4 with prime = %d\n", gens->field_char);
    }
    long nlins = 0;
    uint64_t *linvars = calloc(gens->nvars, sizeof(uint64_t));
    uint32_t **lineqs_ptr = malloc(sizeof(uint32_t *));

    /* int64_t nb = f4_julia(bld, blen, bexp, bcf,
     *                       gens->lens, gens->exps, (void *)gens->cfs, gens->field_char,
     *                       0, //mon_order,
     *                       gens->nvars, gens->ngens, initial_hts,
     *                       nr_threads, max_pairs, update_ht, la_option,
     *                       1, //reduce_gb
     *                       0, //generate_pbm
     *                       info_level); */
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* data structures for basis, hash table and statistics */
    bs_t *bs    = NULL;
    ht_t *bht   = NULL;
    stat_t *st  = NULL;

    int success = 0;

    success = initialize_f4_input_data(&bs, &bht, &st,
            gens->lens, gens->exps, (void *)gens->cfs,
            gens->field_char, 0, gens->nvars, gens->ngens,
            initial_hts, nr_threads, max_pairs,
            update_ht, la_option, 1, 0, info_level);

    if (!success) {
        printf("Bad input data, stopped computation.\n");
        exit(1);
    }

    success = core_f4(&bs, &bht, &st);

    if (!success) {
        printf("Problem with F4, stopped computation.\n");
        exit(1);
    }
    int64_t nb  = export_results_from_f4(bld, blen, bexp,
            bcf, &bs, &bht, &st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 1) {
        print_final_statistics(stderr, st);
    }

    if(nb==0){
        fprintf(stderr, "Something went wrong during the computation\n");
        return -1;
    }

    int32_t *bcf_ff = (int32_t *)(*bcf);

#if DEBUGGB > 0
    if(files->out_file != NULL){
        char fn[200];
        char spr[100];
        sprintf(spr, "%d", gens->field_char);
        strcpy(fn, "/tmp/sys.basis.");
        strcat(fn,spr);
        FILE *ofilebs = fopen(fn, "w");
        display_basis_maple(ofilebs, bld, blen, bexp, bcf_ff, gens);
        fclose(ofilebs);
    }
    else{
        display_basis_maple(stdout, bld, blen, bexp, bcf_ff, gens);
    }
#endif

#if DEBUGGB>0
    if(files->out_file != NULL){
        char fn[200];
        char spr[100];
        sprintf(spr, "%d", gens->field_char);
        strcpy(fn, "/tmp/sys_lead_monomial.");
        strcat(fn,spr);
        FILE *ofilelm = fopen(fn, "w");
        display_lead_monomials_from_gb(ofilelm, bld, blen, bexp, bcf_ff, gens);
        fclose(ofilelm);
    }
    else{
        display_lead_monomials_from_gb(stdout, bld, blen, bexp, bcf_ff, gens);
        fprintf(stdout, "\n");
    }
#endif

    double st_lm = realtime();
    int32_t *bexp_lm = get_lead_monomials(bld, blen, bexp, gens);
    check_and_set_linear_poly_non_hashed(&nlins, linvars, lineqs_ptr,
            bld, bexp_lm,
            blen, bexp, bcf_ff, gens->nvars);

    if(has_dimension_zero(bld[0], gens->nvars, bexp_lm)){
        if(info_level > 1){
            fprintf(stderr, "The ideal has dimension zero\n");
        }
        long dquot = 0;
        int32_t *lmb = monomial_basis(bld[0], gens->nvars, bexp_lm, &dquot);

        if(info_level){
            fprintf(stderr, "Dimension of quotient ring = %ld\n", dquot);
        }
        if(dquot==0){
            fprintf(stderr, "\nGrobner basis is [1]\n");

            free(linvars);
            /* free(lineqs_ptr[0]); */
            free(lineqs_ptr);
            return 3;
        }
#if DEBUGGB>0
        char fn[200];
        char spr[100];
        sprintf(spr, "%d", gens->field_char);
        strcpy(fn, "/tmp/sys.mon_basis.");
        strcat(fn,spr);
        FILE *ofile = fopen(fn, "w");
        display_monomials_from_array_maple(ofile, dquot, lmb, gens);
        fclose(ofile);
#endif

        if(dquot==1){
            print_ff_basis_data(
                    files->out_file, "a", bs, bht, st, gens, print_gb);
            return -1;
        }
        if(info_level > 1){
            fprintf(stderr, "Build monomial basis: %.2f sec.\n", realtime()-st_lm);
        }
        sp_matfglm_t *matrix = build_matrixn(lmb, dquot, bld[0], blen, bexp, bcf_ff,
                bexp_lm, gens->nvars, gens->field_char);

        if(matrix==NULL){
            /* fprintf(stderr, "Problem when building the matrix multiplication\n");
             * fprintf(stderr, "Some data should be free-ed\n"); */
            free(bexp_lm);
            free(lmb);
            if (dquot > 0) {
                return 1;
            } else {
                print_ff_basis_data(
                        files->out_file, "a", bs, bht, st, gens, print_gb);
                return -1;
            }
        }
#if DEBUGGB > 0
        display_fglm_matrix(stdout, matrix);
#endif

        /* timings */
        double rt0, rt1;
        rt0 = realtime();

        if(info_level > 1){
            fprintf(stderr, "Starts FGLM\n");
        }
        uint64_t *squvars = calloc(gens->nvars-1, sizeof(uint64_t));
        check_and_set_vars_squared_in_monomial_basis(squvars, lmb,
                dquot, gens->nvars);
        *bparam = nmod_fglm_compute(matrix, gens->field_char, gens->nvars,
                nlins, linvars, lineqs_ptr[0],
                squvars, info_level);
#if DEBUGGB>0
        if(files->out_file != NULL){
            FILE *ofile = fopen(files->out_file, "w");
            display_fglm_param_maple(ofile, *bparam);
            fclose(ofile);
        }
        else{
            display_fglm_param_maple(stdout, *bparam);
        }
#endif

        free_sp_mat_fglm(matrix);
        free(bexp_lm);
        free(lmb);
        free(linvars);
        if(nlins){
            free(lineqs_ptr[0]);
        }
        free(lineqs_ptr);
        free(squvars);

        /* timings */
        rt1 = realtime();
        if (info_level > 0) {
            fprintf(stderr, "FGLM time (elapsed) %.2f sec\n",
                    rt1-rt0);
        }
        print_ff_basis_data(
                files->out_file, "a", bs, bht, st, gens, print_gb);
        if(bparam==NULL){
            return 1;
        }
        return 0;
    }
    /* else{ */
    /*     print_ff_basis_data( */
    /*             files->out_file, "a", bs, bht, st, gens, print_gb); */
    /*     fprintf(stderr, "The ideal is not zero-dimensional\n"); */

    /*     return 1; */
    /* } */
    print_ff_basis_data(
            files->out_file, "a", bs, bht, st, gens, print_gb);
    free(linvars);
    if(nlins){
        free(lineqs_ptr[0]);
    }
    free(lineqs_ptr);

    if(info_level){
        fprintf(stderr, "Positive dimensional Grobner basis\n");
    }
    return 2;
}


/**

   Renvoie 0 si tout s'est bien passe (ideal de dim 0 en position generique)

 **/

int modular_run_msolve(param_t **bparam,
                       //                       int32_t *bld, int32_t **blen, int32_t **bexp, void **bcf,
                       data_gens_ff_t *gens,
                       int32_t initial_hts,
                       int32_t nr_threads,
                       int32_t max_pairs,
                       int32_t update_ht,
                       int32_t la_option,
                       int32_t info_level,
                       files_gb *files,
                       int32_t prime){
  long nterms = 0;
  for(long i = 0; i < gens->ngens; i++){
    nterms += gens->lens[i];
  }

  reduce_generators((gens->mpz_cfs), nterms, gens->cfs, prime);

  gens->field_char = prime;

  int b = msolve_ff(bparam,
                    //bld, blen, bexp, bcf,
                    gens,
                    initial_hts, nr_threads, max_pairs, update_ht, la_option,
                    info_level,
                    files);
  gens->field_char = 0;

  return b;
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
    mpz_param->coords[j]->length = nmod_param->coords[j]->length;
  }
}

static inline void crt_lift_mpz_upoly(mpz_upoly_t pol, nmod_poly_t nmod_pol,
                                      mpz_t modulus, int32_t prime,
                                      mpz_t prod,
                                      int nthrds){
  len_t i;

#pragma omp parallel for num_threads(nthrds)    \
  private(i) schedule(static)
  for(i = 0; i < pol->length; i++){

    mpz_CRT_ui(pol->coeffs[i], pol->coeffs[i], modulus,
               nmod_pol->coeffs[i], prime, prod, 1);

  }
}


/* assumes that all degrees are the same */
static inline void crt_lift_mpz_param(mpz_param_t mpz_param, param_t *nmod_param,
                                      mpz_t modulus, int32_t prime, int nthrds){
  mpz_t prod;
  mpz_init(prod);
  mpz_mul_ui(prod, modulus, prime);
  crt_lift_mpz_upoly(mpz_param->elim, nmod_param->elim, modulus, prime,
                     prod, nthrds);
  for(long i = 0; i < mpz_param->nvars - 1; i++){
    crt_lift_mpz_upoly(mpz_param->coords[i], nmod_param->coords[i],
                       modulus, prime, prod, nthrds);
  }

  mpz_clear(prod);

}

void initialize_rrec_data(rrec_data_t recdata){
  mpz_init(recdata->r0);
  mpz_init(recdata->r1);
  mpz_init(recdata->t0);
  mpz_init(recdata->t1);
  mpz_init(recdata->q);
  mpz_init(recdata->tmp);
}

void free_rrec_data(rrec_data_t recdata){
  mpz_clear(recdata->r0);
  mpz_clear(recdata->r1);
  mpz_clear(recdata->t0);
  mpz_clear(recdata->t1);
  mpz_clear(recdata->q);
  mpz_clear(recdata->tmp);
}



int mqrr(mpz_t n, mpz_t d, mpz_t u, mpz_t mod, mpz_t T,
         rrec_data_t rdata){
  if(mpz_cmp_ui(u, 0) == 0){
    if(mpz_cmp(mod, T) > 0){
      mpz_set_ui(n, 0);
      mpz_set_ui(d, 1);
      return 1;
    }
    else{
      return 0;
    }
  }
  mpz_set_ui(n, 0);
  mpz_set_ui(d, 0);

  mpz_set_ui(rdata->t0, 0);
  mpz_set(rdata->r0, mod);

  mpz_set_ui(rdata->t1, 1);
  mpz_set(rdata->r1, u);

  while(( mpz_cmp_ui(rdata->r1, 0) ) && ( mpz_cmp(rdata->r0, T) > 0 )){
    fprintf(stderr, "*");
    mpz_fdiv_q(rdata->q, rdata->r0, rdata->r1);
    if(mpz_cmp(rdata->q, T) > 0){
      mpz_set(n, rdata->r1);
      mpz_set(d, rdata->t1);
      mpz_set(T, rdata->q);
    }
    mpz_mul(rdata->tmp, rdata->q, rdata->r1);
    mpz_sub(rdata->tmp, rdata->r0, rdata->tmp);

    mpz_set(rdata->r0, rdata->r1);
    mpz_set(rdata->r1, rdata->tmp);

    mpz_mul(rdata->tmp, rdata->q, rdata->t1);
    mpz_sub(rdata->tmp, rdata->t0, rdata->tmp);

    mpz_set(rdata->t0, rdata->t1);
    mpz_set(rdata->t1, rdata->tmp);
  }
  mpz_gcd(rdata->q, n, d);
  if((mpz_cmp_ui(d, 0) == 0) || mpz_cmp_ui(rdata->q, 1) != 0){
    fprintf(stderr, "<<%d>>", mpz_cmp_ui(d, 0));
    return 0;
  }
  if(mpz_sgn(d) == -1){
    mpz_neg(n, n);
    mpz_neg(d, d);
  }
  return 1;
}


int mqrr_rat_recon(mpz_t *num, mpz_t *den, mpz_t *lu, 
                   mpz_t mod, mpz_t T,
                   const long min, const long max,
                   rrec_data_t rdata){
  mpz_set_ui(T, mpz_sizeinbase(mod,2));
  mpz_mul_ui(T, T, mpz_sizeinbase(mod,2));
  mpz_mul_2exp(T, T, 1);
  mpz_add_ui(T, T, 1);
  mpz_out_str(stderr, 10, mod);fprintf(stderr, "->\n");
  mpz_out_str(stderr, 10, T);fprintf(stderr, "->");
  for(long i = min; i < max ; i++){
    int b = mqrr(num[i], den[i], lu[i], mod, T, rdata);
    if(b == 0) {
      fprintf(stderr, "i=%ld[%ld,%ld]", i, min, max);
      return 0;
    }
  }
  return 1;
}

/**

   la sortie est recons / denominator
 **/

static inline int rational_reconstruction_mpz_ptr(mpz_t *recons,
                                                  mpz_t denominator,
                                                  mpz_t *pol,
                                                  long len,
                                                  mpz_t modulus,
                                                  long *maxrec,
                                                  mpq_t *coef,
                                                  mpz_t *tmp_num,
                                                  mpz_t *tmp_den,
                                                  mpz_t lcm,
                                                  mpz_t guessed_num,
                                                  mpz_t guessed_den,
                                                  rrec_data_t rdata,
                                                  int info_level){

  if(mpq_reconstruct_mpz(coef, pol[*maxrec], modulus) == 0){
    return 0;
  }

  /* mpq_canonicalize(*coef); */
  mpz_set(tmp_num[*maxrec], mpq_numref(*coef));
  mpz_set(tmp_den[*maxrec], mpq_denref(*coef));

  /* int boo = mqrr_rat_recon(tmp_num, tmp_den, pol, modulus, lcm, *maxrec, len, rdata); */
  /* fprintf(stderr, "(%d)\n", boo); */
  for(long i = *maxrec + 1; i < len; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, pol[i],
                                           modulus,
                                           guessed_num, guessed_den);
    if(b == 0){
      *maxrec = i - 1;
      return b;
    }
    /* mpq_canonicalize(*coef); */

    mpz_set(tmp_num[i], mpq_numref(*coef));
    mpz_set(tmp_den[i], mpq_denref(*coef));
  }
  for(long i = 0; i < *maxrec; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, pol[i],
                                           modulus,
                                           guessed_num, guessed_den);

    if(b == 0){
      if(info_level){
        fprintf(stderr, "[*]");
      }
      *maxrec = MAX(i-1, 0);
      return b;
    }
    /* mpq_canonicalize(*coef); */

    mpz_set(tmp_num[i], mpq_numref(*coef));
    mpz_set(tmp_den[i], mpq_denref(*coef));
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
                                         tmp_num->coeffs,
                                         tmp_den->coeffs,
                                         lcm, guessed_num, guessed_den,
                                         rdata,
                                         info_level);
}


/**

   la sortie est recons / denominator
 **/

static inline int rational_reconstruction_upoly_old(mpz_upoly_t recons,
                                                mpz_t denominator,
                                                mpz_upoly_t pol,
                                                long len,
                                                mpz_t modulus,
                                                long *maxrec,
                                                mpq_t *coef,
                                                mpz_upoly_t tmp_num,
                                                mpz_upoly_t tmp_den,
                                                mpz_t lcm,
                                                mpz_t guessed_num,
                                                mpz_t guessed_den){
  if(mpq_reconstruct_mpz(coef, pol->coeffs[*maxrec], modulus) == 0){
    return 0;
  }
  mpz_set(tmp_num->coeffs[*maxrec], mpq_numref(*coef));
  mpz_set(tmp_den->coeffs[*maxrec], mpq_denref(*coef));

  for(long i = *maxrec + 1; i < len; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, pol->coeffs[i],
                                           modulus,
                                           guessed_num, guessed_den);
    if(b == 0){
      *maxrec = i - 1;
      return b;
    }
    mpz_set(tmp_num->coeffs[i], mpq_numref(*coef));
    mpz_set(tmp_den->coeffs[i], mpq_denref(*coef));
  }
  for(long i = 0; i < *maxrec; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, pol->coeffs[i],
                                           modulus,
                                           guessed_num, guessed_den);

    if(b == 0){
      fprintf(stderr, "[%ld]->[%ld]", *maxrec, i);
      *maxrec = MAX(i-1, 0);
      return b;
    }
    mpz_set(tmp_num->coeffs[i], mpq_numref(*coef));
    mpz_set(tmp_den->coeffs[i], mpq_denref(*coef));
  }

  mpz_set(lcm, tmp_den->coeffs[0]);
  for(long i = 1; i < len; i++){
    mpz_lcm(lcm, lcm, tmp_den->coeffs[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_divexact(tmp_den->coeffs[i], lcm, tmp_den->coeffs[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_mul(tmp_num->coeffs[i], tmp_num->coeffs[i], tmp_den->coeffs[i]);
  }
  for(long i = 0; i < len; i++){
    mpz_set(recons->coeffs[i], tmp_num->coeffs[i]);
  }
  mpz_set(denominator, lcm);

  return 1;

}


/**

renvoie 0 si on n'a pas reussi a reconstruire

numer contient le resultat

**/

static inline int rational_reconstruction(mpz_param_t mpz_param,
                                          param_t *nmod_param,
                                          mpz_upoly_t numer,
                                          mpz_upoly_t denom,
                                          mpz_t *modulus, int32_t prime,
                                          mpq_t *coef,
                                          mpz_t *guessed_num, mpz_t *guessed_den,
                                          long *maxrec,
                                          const int info_level){

  crt_lift_mpz_param(mpz_param, nmod_param, *modulus, prime, 1);
  mpz_mul_ui(*modulus, *modulus, prime);

  mpz_fdiv_q_2exp(*guessed_num, *modulus, 1);
  mpz_sqrt(*guessed_num, *guessed_num);
  mpz_set(*guessed_den, *guessed_num);

  long deg = mpz_param->nsols;
  if(mpq_reconstruct_mpz(coef, mpz_param->elim->coeffs[*maxrec], *modulus) == 0){
    return 0;
  }

  mpz_set(numer->coeffs[*maxrec], mpq_numref(*coef));
  mpz_set(denom->coeffs[*maxrec], mpq_denref(*coef));

  for(long i = *maxrec + 1; i <= deg; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, mpz_param->elim->coeffs[i],
                                           *modulus,
                                           *guessed_num, *guessed_den);
    if(b == 0){
      *maxrec = i - 1;
      return b;
    }
    mpz_set(numer->coeffs[i], mpq_numref(*coef));
    mpz_set(denom->coeffs[i], mpq_denref(*coef));
  }
  for(long i = 0; i < *maxrec; i++){
    int b = mpq_reconstruct_mpz_with_denom(coef, mpz_param->elim->coeffs[i],
                                           *modulus,
                                           *guessed_num, *guessed_den);

    if(b == 0){
      if (info_level) {
        fprintf(stderr, "[%ld]->[%ld]", *maxrec, i);
      }
      *maxrec = MAX(i-1, 0);
      return b;
    }
    mpz_set(numer->coeffs[i], mpq_numref(*coef));
    mpz_set(denom->coeffs[i], mpq_denref(*coef));
  }

  mpz_t lcm;
  mpz_init_set(lcm, denom->coeffs[0]);
  for(long i = 1; i <= deg; i++){
    mpz_lcm(lcm, lcm, denom->coeffs[i]);
  }
  for(long i = 0; i <= deg; i++){
    mpz_divexact(denom->coeffs[i], lcm, denom->coeffs[i]);
  }
  for(long i = 0; i <= deg; i++){
    mpz_mul(numer->coeffs[i], numer->coeffs[i], denom->coeffs[i]);
  }
  mpz_clear(lcm);
  return 1;
}


/**

returns 0 if rational reconstruction failed

 **/

//#define TRY 1
static inline int new_rational_reconstruction(mpz_param_t mpz_param,
                                              mpz_param_t tmp_mpz_param,
                                              param_t *nmod_param,
                                              mpz_upoly_t numer,
                                              mpz_upoly_t denom,
                                              mpz_t *modulus, int32_t prime,
                                              mpq_t *coef,
                                              rrec_data_t recdata,
                                              mpz_t *guessed_num,
                                              mpz_t *guessed_den,
                                              long *maxrec,
                                              int *is_lifted,
                                              int doit,
                                              int nthrds,
                                              const int info_level){

  crt_lift_mpz_param(tmp_mpz_param, nmod_param, *modulus, prime, nthrds);
  mpz_mul_ui(*modulus, *modulus, prime);

  if(doit==0){
    return 0;
  }

  mpz_fdiv_q_2exp(*guessed_num, *modulus, 1);
  mpz_sqrt(*guessed_num, *guessed_num);
  mpz_set(*guessed_den, *guessed_num);

  mpz_t denominator;
  mpz_init(denominator);
  mpz_t lcm;
  mpz_init(lcm);
  int b = 0;

  double st = realtime();

  if(is_lifted[0]==0){
      b = rational_reconstruction_upoly(mpz_param->elim,
                                        denominator,
                                        tmp_mpz_param->elim,
                                        nmod_param->elim->length,
                                        *modulus,
                                        maxrec,
                                        coef,
                                        numer,
                                        denom,
                                        lcm,
                                        *guessed_num,
                                        *guessed_den,
                                        recdata,
                                        info_level);
      if(b == 0){
        is_lifted[0] = 0;
        mpz_clear(denominator);
        mpz_clear(lcm);
        return b;
      }
      else{
        is_lifted[0] = 1;
      }
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

  for(int i = 0; i < nc; i++){
    *maxrec = 0;

    if(is_lifted[0]>0 && is_lifted[i+1]==0){

      b = rational_reconstruction_upoly(mpz_param->coords[i],
                                        denominator,
                                        tmp_mpz_param->coords[i],
                                        nmod_param->coords[i]->length,
                                        *modulus,
                                        maxrec,
                                        coef,
                                        numer,
                                        denom,
                                        lcm,
                                        *guessed_num,
                                        *guessed_den,
                                        recdata,
                                        info_level);
    }
    if(b == 0){
      mpz_clear(denominator);
      mpz_clear(lcm);
      mpz_clear(lc);
      mpq_clear(c);

      is_lifted[i+1] = 0;
      return b;
    }
    else{
      is_lifted[i+1] = 1;
      if(info_level && b){
        fprintf(stderr, "[%d]", i+1);
      }
    }

    mpz_set(mpq_denref(c), denominator);
    mpz_mul_ui(mpq_numref(c), lc, 1);
    mpq_canonicalize(c);
    mpz_set(mpz_param->cfs[i], mpq_denref(c));

    for(long j = 0; j < mpz_param->coords[i]->length; j++){
      mpz_mul(mpz_param->coords[i]->coeffs[j], mpz_param->coords[i]->coeffs[j],
              mpq_numref(c));
    }
  }
  mpz_clear(denominator);
  mpz_clear(lcm);
  mpz_clear(lc);
  mpq_clear(c);
  /* fprintf(stderr, "RATRECON time = %.2f\n", realtime()-st); */
  return 1;
}



/**
renvoie 1 si il faut faire le modular check.
A terme on va pouvoir enlever cette fonction.
**/

/* static inline int check_param_modular_elim(mpz_param_t mp_param, */
/*                                            mpz_upoly_t numer, */
/*                                            param_t *bparam, */
/*                                            int32_t prime){ */
/*   uint32_t lc = mpz_fdiv_ui(numer->coeffs[mp_param->nsols], prime); */
/*   lc = mod_p_inverse_32(lc, prime); */
/*   for(long i = 0; i <= mp_param->nsols; i++){ */
/*     uint64_t c = mpz_fdiv_ui(numer->coeffs[i], prime); */
/*     c *= (uint64_t)lc; */
/*     c = c % prime; */
/*     if(c != bparam->elim->coeffs[i]){ */
/*       return 1; */
/*     } */
/*   } */
/*   return 0; */
/* } */

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

/* renvoie 0 si c'est bon sinon on renvoie l'indice du coeff problematique + 1  */
static inline int check_proportional_mpz_nmod_poly(const long len,
                                                   const mpz_upoly_t mpz_pol,
                                                   const nmod_poly_t nm_pol,
                                                   const int32_t prime){
  if(len == 0){
    return 0;
  }
  uint32_t lc = mpz_fdiv_ui(mpz_pol->coeffs[len - 1], prime);
  uint32_t nmodlc = nm_pol->coeffs[len - 1] % prime;
  lc = mod_p_inverse_32(lc, prime);
  nmodlc = mod_p_inverse_32(nmodlc, prime);
  for(long i = 0; i < len; i++){
    uint64_t c = mpz_fdiv_ui(mpz_pol->coeffs[i], prime);
    c *= (uint64_t)lc;
    c = c % prime;

    uint64_t nmod_coef = nm_pol->coeffs[i] * nmodlc;
    nmod_coef = nmod_coef % prime;
    if(c != nmod_coef){
      return 1;
    }
  }
  return 0;
}


/**
   renvoie 1 si il faut faire le modular check. 
**/

static inline int check_param_modular(const mpz_param_t mp_param,
                                      const param_t *bparam,
                                      const int32_t prime,
                                      int *is_lifted,
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
    return 1;
  }

  for(int i = 0; i <mp_param->nvars-1; i++ ){
    len = mp_param->coords[0]->length;
    if(check_proportional_mpz_nmod_poly(bparam->coords[i]->length,
                                        mp_param->coords[i],
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


/**

renvoie 0 si qqc plante.

 **/

int msolve_prob_linalg_qq(mpz_param_t mp_param,
                   /* int32_t *bld, */
                   /* int32_t **blen, int32_t **bexp, void **bcf, */
                   param_t **bparam, data_gens_ff_t *gens,
                   int32_t initial_hts,
                   int32_t nr_threads,
                   int32_t max_pairs,
                   int32_t update_ht,
                   int32_t la_option,
                   int32_t info_level,
                   int32_t pbm_file,
                   files_gb *files){

  int32_t prime = next_prime(1<<30);
  fprintf(stderr, "{%d}", prime);
  int b = modular_run_msolve(bparam, //bld, blen, bexp, bcf,
                             gens,
                             initial_hts, nr_threads, max_pairs, update_ht,
                             la_option,
                             info_level, files,
                             prime);
  if(b!=0){
    return 0;
  }
  initialize_mpz_param(mp_param, *bparam);
  mpz_t modulus;
  mpz_init_set_ui(modulus, prime);
  fprintf(stderr, "modulus = "); mpz_out_str(stderr, 10, modulus);fprintf(stderr, "\n");

  mpq_t result, test;
  mpq_init(result);
  mpq_init(test);
  normalize_nmod_param(*bparam);
  set_mpz_param_nmod(mp_param, *bparam);

  mpz_upoly_t numer;
  mpz_upoly_init2(numer, (mp_param->nsols + 1), 32*(mp_param->nsols + 1));
  numer->length = mp_param->nsols + 1;

  mpz_upoly_t denom;
  mpz_upoly_init2(denom, (mp_param->nsols + 1), 64*(mp_param->nsols + 1));
  //  mpz_poly_init(denom, 33);
  denom->length = mp_param->nsols + 1;

  mpz_t guessed_den;
  mpz_init2(guessed_den, 32*mp_param->nsols);
  mpz_t guessed_num;
  mpz_init2(guessed_num, 32*mp_param->nsols);

  long maxrec = 0;

  long degree = mp_param->nsols;
  int rerun = 1, nprimes = 1, mcheck =1;

  double rr = 0, tr = 0;
  while(rerun == 1 || mcheck == 1){
    prime = next_prime(prime + 1);
    b = modular_run_msolve(bparam, //bld, blen, bexp, bcf,
                           gens,
                           initial_hts, nr_threads, max_pairs, update_ht,
                           la_option,
                           0, //info_level,
                           files,
                           prime);
    normalize_nmod_param(*bparam);
    if(degree == (*bparam)->elim->length-1 && b == 0){
      double cr0 = cputime();
      double tr0 = realtime();
      if(rerun == 0){
        mcheck = check_unit_mpz_nmod_poly(mp_param->nsols + 1, numer,
                                          (*bparam)->elim,
                                          prime);
      }
      if(rational_reconstruction(mp_param, *bparam,
                                 numer, denom,
                                 &modulus, prime,
                                 &result,
                                 &guessed_num, &guessed_den,
                                 &maxrec, info_level)){
        rerun = 0;
      }
      double cr1 = cputime();
      double tr1 = realtime();
      rr += (cr1-cr0);
      tr += (tr1-tr0);

      nprimes++;
    }
    else{
        if(info_level){
            fprintf(stderr, "<bp: %d>\n", prime);
        }
    }
    if(info_level){
      if( (nprimes & (nprimes - 1)) == 0){
        fprintf(stderr, "{%d}", nprimes);
      }
    }
  }

  fprintf(stderr, "\n");mpz_out_str(stderr, 10, numer->coeffs[degree]); fprintf(stderr, "\n");
  fprintf(stderr, "%d primes used\n", nprimes);
  fprintf(stderr, "Time for rational reconstruction %13.2f sec (elapsed) / %5.2f sec (cpu)\n", rr, tr);

  mpz_upoly_clear(numer);
  mpz_upoly_clear(denom);
  mpz_clear(guessed_num);
  mpz_clear(guessed_den);
  mpq_clear(test);
  mpq_clear(result);
  mpz_clear(modulus);

  return 1;
}

static inline int32_t *get_lm_from_bs(bs_t *bs, const ht_t *ht){
  hm_t *dt;
  const len_t nelts = bs->lml;
  const int nv = ht->nv;
  int32_t *exp  = (int32_t *)malloc(
                                    (unsigned long)(nelts) * (unsigned long)(nv) * sizeof(int32_t));
  /* counters for lengths, exponents and coefficients */
  int64_t cl = 0, ce = 0;//, cc = 0, ctmp  = 0;;


  for (long i = 0; i < nelts; ++i) {
    const bl_t bi = bs->lmps[i];
    //    len[cl] = bs->hm[bi][LENGTH];

    dt  = bs->hm[bi] + OFFSET;
    for (int k = 0; k < nv; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    //    cc  +=  len[cl];
    cl++;
  }
  return exp;
}


static inline void get_lm_from_bs_trace(bs_t *bs, const ht_t *ht, int32_t *exp){
  hm_t *dt;
  const len_t nelts = bs->lml;
  const int nv = ht->nv;

  /* counters for lengths, exponents and coefficients */
  int64_t cl = 0, ce = 0;//, cc = 0, ctmp  = 0;;

  for (long i = 0; i < nelts; ++i) {
    const bl_t bi = bs->lmps[i];
    //    len[cl] = bs->hm[bi][LENGTH];

    dt  = bs->hm[bi] + OFFSET;
    for (int k = 0; k < nv; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    //    cc  +=  len[cl];
    cl++;
  }
}


static inline void set_linear_poly(long nlins, uint32_t *lineqs, uint64_t *linvars,
                                   ht_t *bht, int32_t *bexp_lm, bs_t *bs){

  for(long i = 0; i < nlins*((bht->nv + 1)); i++){
    lineqs[i] = 0;
  }
  for(long i = 0; i < nlins*(bht->nv+1); i++){
    lineqs[i] = 0;
  }
  int cnt = 0;

  for(int i = 0; i < bht->nv; i++){
    if(linvars[i] != 0){

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if(len==bht->nv+1){
        for(long j = 0; j<len; j++){
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          lineqs[cnt*(bht->nv+1)+j] = coef;
        }
      }
      else{
        hm_t *dt = bs->hm[bi] + OFFSET;
        for(long j = 0; j<len; j++){
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          int isvar = 0;
          for(int k = 0; k < bht->nv; k++){
            if(exp[k]==1){
              lineqs[cnt*(bht->nv+1)+k] = coef;
              isvar=1;
            }
          }
          if(isvar==0){
            lineqs[cnt*(bht->nv+1)+bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
}


static inline void check_and_set_linear_poly(long *nlins_ptr, uint64_t *linvars,
                                             uint32_t** lineqs_ptr,
                                             ht_t *bht, int32_t *bexp_lm, bs_t *bs){
  long nlins = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1. 
  */
  for(long i = 0; i < bs->lml; i++){
    long deg = 0;
    for(int j = 0; j < bht->nv; j++){
      deg+=bexp_lm[i*bht->nv+j];
    }

    if(deg == 1){
      nlins++;
      for(int k = 0; k < bht->nv; k++){
        if(bexp_lm[i*bht->nv+k] == 1){
          linvars[k]=i + 1;
        }
      }
    }
  }

  *nlins_ptr = nlins;

  //On recupere les coefficients des formes lineaires
  uint32_t* lineqs = calloc(nlins*(bht->nv + 1), sizeof(uint64_t));
  for(long i = 0; i < nlins*(bht->nv+1); i++){
    lineqs[i] = 0;
  }
  int cnt = 0;
  for(int i = 0; i < bht->nv; i++){
    if(linvars[i] != 0){

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if(len==bht->nv+1){
        for(long j = 0; j<len; j++){
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          lineqs[cnt*(bht->nv+1)+j] = coef;
        }
      }
      else{
        hm_t *dt = bs->hm[bi] + OFFSET;
        for(long j = 0; j<len; j++){
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          int isvar = 0;
          for(int k = 0; k < bht->nv; k++){
            if(exp[k]==1){
              lineqs[cnt*(bht->nv+1)+k] = coef;
              isvar=1;
            }
          }
          if(isvar==0){
            lineqs[cnt*(bht->nv+1)+bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  lineqs_ptr[0] = lineqs;
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
                                        ht_t *tht,
                                        const bs_t *bs_qq,
                                        ht_t *bht,
                                        stat_t *st,
                                        const int32_t fc,
                                        int info_level,
                                        int print_gb,
                                        int *dim,
                                        long *dquot_ori,
                                        data_gens_ff_t *gens,
                                        files_gb *files,
                                        int *success)
{
    double ca0, rt;
    ca0 = realtime(); 
    bs_t *bs = f4_trace_learning_phase(trace, tht, bs_qq, bht, st, fc);
    rt = realtime()-ca0;

    if(info_level > 1){
        fprintf(stderr, "Learning phase %.2f Gops/sec\n",
                (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/rt);
        fprintf(stderr, "------------------------------------------\n");
        fprintf(stderr, "#ADDITIONS       %13lu\n", (unsigned long)st->trace_nr_add * 1000);
        fprintf(stderr, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->trace_nr_mult * 1000);
        fprintf(stderr, "#REDUCTIONS      %13lu\n", (unsigned long)st->trace_nr_red);
        fprintf(stderr, "------------------------------------------\n");
    }

    /* Leading monomials from Grobner basis */
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
    /* set st->fc to finite characteristic for printing */
    st->fc  = fc;
    print_ff_basis_data(
                        files->out_file, "a", bs, bht, st, gens, print_gb);

    st->fc  = 0;

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
        free_basis(&(bs));
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
                                               success);


        *dim = 0;
        *dquot_ori = dquot;
        return lmb;
    }
    else{
        *dim  = 1;
        *dquot_ori = -1;
        return NULL;
    }
}

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
                                             stat_t *st,
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
                bdata_fglm, bdata_bms, success);
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



static inline int equal_staircase(int32_t *lmb, int32_t *lmb_ori,
                                  long dquot, long dquot_ori,
                                  int32_t nv){
if(dquot != dquot_ori){
  return 0;
}
for(long i = 0; i < dquot; i++){
  for(int32_t j = 0; j < nv; j++){
    if(lmb[i*nv+j] != lmb_ori[i*nv+j]){
      return 0;
    }
  }
}
return 1;
}



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
                               stat_t *st,
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

  if(bs[i]->lml != num_gb[i]){
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
                                          info_level)){
      bad_primes[i] = 1;
    }
    free_basis(&(bs[i]));
  }
  else{
      bad_primes[i] = 1;
  }
 }
}

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
                                   trace_t **btrace,
                                   ht_t **btht,
                                   const bs_t *bs_qq,
                                   ht_t **bht,
                                   stat_t *st,
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

  /* F4 and FGLM are run using a single thread */
  /* st->nthrds is reset to its original value afterwards */
  const int nthrds = st->nthrds;
  st->nthrds = 1 ;

  memset(bad_primes, 0, (unsigned long)st->nprimes * sizeof(int));
  #pragma omp parallel for num_threads(nthrds)  \
    private(i) schedule(static)
  for (i = 0; i < st->nprimes; ++i){
    ca0 = realtime();
    bs[i] = f4_trace_application_phase(btrace[i], btht[i], bs_qq, bht[i], st, lp->p[i]);
    *stf4 = realtime()-ca0;
    /* printf("F4 trace timing %13.2f\n", *stf4); */

    if(bs[i]->lml != num_gb[i]){
      if (bs[i] != NULL) {
        free_basis(&(bs[i]));
      }
      nmod_params[i] = NULL;
      bad_primes[i] = 1;
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
                                            info_level)){
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




/*

  - renvoie 0 si le calcul est ok.
  => GB = [1] dim =0 dquot = 0
  => Positive dimension dim > 0
  => Dimension zero + calcul qui a pu etre fait. dim=0 dquot > 0

  - renvoie 1 si le calcul a echoue.
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
                    int32_t reset_ht,
                    int32_t la_option,
                    int32_t info_level,
                    int32_t print_gb,
                    int32_t pbm_file,
                    files_gb *files,
                    int round){

  const int32_t *lens = gens->lens;
  const int32_t *exps = gens->exps;
  const void *cfs = gens->mpz_cfs;
  const uint32_t field_char = gens->field_char;
  const int mon_order = 0;
  const int32_t nr_vars = gens->nvars;
  const int32_t nr_gens = gens->ngens;
  const int reduce_gb = 1;
  const uint32_t prime_start = pow(2, 30);
  const int32_t nr_primes = nr_threads;

  /* only for computations over the rationals */
  if (field_char != 0) {
      fprintf(stderr, "Tracer only for computations over Q. Call\n");
      fprintf(stderr, "standard F4 Algorithm for computations over\n");
      fprintf(stderr, "finite fields.\n");
      return -2;
  }
  len_t i;

  /* initialize stuff */
  stat_t *st  = initialize_statistics();

  /* checks and set all meta data. if a nonzero value is returned then
    * some of the input data is corrupted. */
  if (check_and_set_meta_data_trace(st, lens, exps, cfs, field_char,
              mon_order, nr_vars, nr_gens, ht_size, nr_threads,
              max_nr_pairs, reset_ht, la_option, reduce_gb, prime_start,
              nr_primes, pbm_file, info_level)) {
    free(st);
    return -3;
  }

  /* lucky primes */
  primes_t *lp  = (primes_t *)calloc(st->nthrds, sizeof(primes_t));

  /*******************
  * initialize basis
  *******************/
  bs_t *bs_qq = initialize_basis(st->ngens);
  /* initialize basis hash table, update hash table, symbolic hash table */
  ht_t *bht = initialize_basis_hash_table(st);
  /* hash table to store the hashes of the multiples of
    * the basis elements stored in the trace */
  ht_t *tht = initialize_secondary_hash_table(bht, st);
  /* read in ideal, move coefficients to integers */
  import_julia_data(bs_qq, bht, st, lens, exps, cfs);

  if (st->info_level > 0) {
    print_initial_statistics(stderr, st);
  }

  /* for faster divisibility checks, needs to be done after we have
    * read some input data for applying heuristics */
  calculate_divmask(bht);

  /* sort initial elements, smallest lead term first */
  sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *),
          initial_input_cmp, bht);
  remove_content_of_initial_basis(bs_qq);

  /* generate lucky prime numbers */
  generate_lucky_primes(lp, bs_qq, st->prime_start, st->nthrds);

  /* generate array to store modular bases */
  bs_t **bs = (bs_t **)calloc((unsigned long)st->nthrds, sizeof(bs_t *));

  param_t **nmod_params =  (param_t **)calloc((unsigned long)st->nthrds,
                                              sizeof(param_t *));

  int *bad_primes = calloc((unsigned long)st->nthrds, sizeof(int));

  /* initialize tracers */
  trace_t **btrace = (trace_t **)calloc(st->nthrds,
                                       sizeof(trace_t *));
  btrace[0]  = initialize_trace();
  /* initialization of other tracers is done through duplication */

  srand(time(0));
  uint32_t prime = next_prime(1<<30);
  prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  while(is_lucky_prime_ui(prime, bs_qq)){
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
  }

  uint32_t primeinit = prime;
  lp->p[0] = primeinit;

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
  long *bnlins = (long *)malloc(sizeof(long) * st->nthrds);
  uint64_t **blinvars = (uint64_t **)malloc(sizeof(uint64_t *) * st->nthrds);
  uint64_t *linvars = calloc(bht->nv, sizeof(uint64_t));
  blinvars[0] = linvars;
  uint32_t **lineqs_ptr = malloc(sizeof(uint32_t *) * st->nthrds);
  uint64_t **bsquvars = (uint64_t **) malloc(sizeof(uint64_t *) * st->nthrds);
  uint64_t *squvars = calloc(nr_vars-1, sizeof(uint64_t));
  bsquvars[0] = squvars;


  int success = 1;
  int squares = 1;
  int32_t *lmb_ori = modular_trace_learning(bmatrix, bdiv_xn, blen_gb_xn,
                                            bstart_cf_gb_xn,

                                            &nlins, blinvars[0], lineqs_ptr,
                                            squvars,

                                            bdata_fglm, bdata_bms,

                                            num_gb, leadmons_ori,

                                            &bsz, nmod_params, btrace[0],
                                            tht, bs_qq, bht, st,
                                            lp->p[0], //prime,
                                            info_level,
                                            print_gb,
                                            dim_ptr, dquot_ptr,
                                            gens,
                                            files,
                                            &success);


  if(*dim_ptr == 0 && success && *dquot_ptr > 0){
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

  if(lmb_ori == NULL || success == 0) {
      /* print_msolve_message(stderr, 1); */
    for(int i = 0; i < st->nthrds; i++){
      free_trace(&btrace[i]);
    }
    free_shared_hash_data(bht);
    free_hash_table(&bht);
    free_hash_table(&tht);

    for (i = 0; i < st->nthrds; ++i) {
      //      free_basis(&(bs[i]));
    }
    free(bs);
    //here we should clean nmod_params
    free_lucky_primes(&lp);
    free(st);
    free(linvars);
    free(nmod_params);
    if(nlins){
      free(lineqs_ptr[0]);
    }
    free(lineqs_ptr);
    free(squvars);
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

  /* duplicate data for multi-threaded multi-mod computation */
  duplicate_data_mthread_trace(st->nthrds, st, num_gb,
                               leadmons_ori, leadmons_current,
                               btrace,
                               bdata_bms, bdata_fglm,
                               bstart_cf_gb_xn, blen_gb_xn, bdiv_xn, bmatrix,
                               nmod_params, nlins, bnlins,
                               blinvars, lineqs_ptr,
                               bsquvars);

  /* copy of hash tables for tracer application */
  ht_t **blht = (ht_t **)malloc((st->nthrds) * sizeof(ht_t *));
  blht[0] = bht;
  for(int i = 1; i < st->nthrds; i++){
    ht_t *lht = copy_hash_table(bht, st);
    blht[i] = lht;
  }
  ht_t **btht = (ht_t **)malloc((st->nthrds) * sizeof(ht_t *));
  btht[0] = tht;
  for(int i = 1; i < st->nthrds; i++){
    btht[i] = copy_hash_table(tht, st);
  }

  normalize_nmod_param(nmod_params[0]);


  if(info_level){
    fprintf(stderr, "\nStarts trace based multi-modular computations\n");
  }

  mpz_param_t tmp_mpz_param;
  mpz_param_init(tmp_mpz_param);

  initialize_mpz_param(mpz_param, nmod_params[0]);
  initialize_mpz_param(tmp_mpz_param, nmod_params[0]);
  //attention les longueurs des mpz_param sont fixees par nmod_params[0]
  //dans des cas exceptionnels, ca peut augmenter avec un autre premier. 

  mpz_t modulus;
  mpz_init_set_ui(modulus, primeinit);

  mpq_t result, test;
  mpq_init(result);
  mpq_init(test);
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

  int rerun = 1, nprimes = 1, mcheck =1;

  long nbadprimes = 0;

  int *is_lifted = malloc(sizeof(int)*nr_vars);
  for(int i = 0; i < nr_vars; ++i){
    is_lifted[i] = 0;
  }
  int nbdoit = 1;
  int doit = 1;
  int prdone = 0;
  int lpow2 = 0;
  int clog = 0;
  int br = 0;
  prime = next_prime(1<<30);

  rrec_data_t recdata;
  initialize_rrec_data(recdata);

  /* measures time spent in rational reconstruction */
  double strat = 0;
  while(rerun == 1 || mcheck == 1){

    /* controls call to rational reconstruction */
    doit = ((prdone %nbdoit) == 0);

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
    modular_trace_application(bmatrix,
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
                              nmod_params, btrace,
                              btht, bs_qq, blht, st,
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
                                       is_lifted, info_level);
        }
        crr = realtime();
        if(mcheck==1){
          br = new_rational_reconstruction(mpz_param,
                                           tmp_mpz_param,
                                           nmod_params[i],
                                           numer, denom,
                                           &modulus, lp->p[i],
                                           &result, recdata,
                                           &guessed_num, &guessed_den,
                                           &maxrec,
                                           is_lifted,
                                           doit,
                                           st->nthrds, info_level);
          if(br == 1){
            rerun = 0;
          }
          else{
            rerun = 1;
          }
        }
        scrr += realtime()-crr;
      }
      else{
        if(info_level){
          fprintf(stderr, "<bp: %d>\n", lp->p[i]);
        }
        nbadprimes++;
        if(nbadprimes > nprimes){
          free(linvars);
          free(lineqs_ptr[0]);
          free(lineqs_ptr);
          free(squvars);
          free_rrec_data(recdata);
          return -4;
        }
      }
      nprimes++;
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
  mpz_clear(modulus);
  free_rrec_data(recdata);

  /* free and clean up */
  free_shared_hash_data(bht);
  for(int i = 0; i < st->nthrds; i++){
    free_hash_table(blht+i);
    free_hash_table(btht+i);
  }
  /* free_hash_table(&bht); */
  /* free_hash_table(&tht); */

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
    free_trace(&btrace[i]);
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
  free(st);
  free(bad_primes);
  free(blinvars);
  free(lineqs_ptr);
  free(bsquvars);
  free(is_lifted);
  free(num_gb);
  free(blen_gb_xn);
  free(bstart_cf_gb_xn);
  free(bdiv_xn);
  free(btrace);

  return 0;
}


/*

  - renvoie 0 si le calcul est ok.
  => GB = [1] dim =0 dquot = 0
  => Positive dimension dim > 0
  => Dimension zero + calcul qui a pu etre fait. dim=0 dquot > 0

  - renvoie 1 si le calcul a echoue.
  => Dimension 0 => pas en position generique

  - renvoie 2 si besoin de plus de genericite.
  => (tous les carres ne sont pas sous l'escalier)

  - renvoie -2 si la carac est > 0

  - renvoie -3 si meta data pas bonnes

  - renvoie -4 si bad prime
*/

int msolve_probabilistic_qq(mpz_param_t mpz_param,
                            param_t **nmod_param,
                            int *dim_ptr,
                            long *dquot_ptr,
                            data_gens_ff_t *gens,
                            int32_t ht_size, //initial_hts,
                            int32_t nr_threads,
                            int32_t max_nr_pairs,
                            int32_t reset_ht,
                            int32_t la_option,
                            int32_t info_level,
                            int32_t print_gb,
                            int32_t pbm_file,
                            files_gb *files,
                            int round){

    const int32_t *lens = gens->lens;
    const int32_t *exps = gens->exps;
    const void *cfs = gens->mpz_cfs;
    const uint32_t field_char = gens->field_char;
    const int mon_order = 0;
    const int32_t nr_vars = gens->nvars;
    const int32_t nr_gens = gens->ngens;
    const int reduce_gb = 1;
    const uint32_t prime_start = pow(2, 30);
    const int32_t nr_primes = nr_threads;

    /* only for computations over the rationals */
    if (field_char != 0) {
        fprintf(stderr, "Modular F4 only for computations over Q. Call\n");
        fprintf(stderr, "standard F4 Algorithm for computations over\n");
        fprintf(stderr, "finite fields.\n");
        return 1;
    }

    len_t i;

    /* initialize stuff */
    stat_t *st  = initialize_statistics();

    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */
    if (check_and_set_meta_data_trace(st, lens, exps, cfs, field_char,
                mon_order, nr_vars, nr_gens, ht_size, nr_threads,
                max_nr_pairs, reset_ht, la_option, reduce_gb, prime_start,
                nr_primes, pbm_file, info_level)) {
        free(st);
        return -3;
    }

    /* lucky primes */
    primes_t *lp  = (primes_t *)calloc(1, sizeof(primes_t));

    /*******************
     * initialize basis
     *******************/
    bs_t *bs_qq = initialize_basis(st->ngens);
    /* initialize basis hash table, update hash table, symbolic hash table */
    ht_t *bht = initialize_basis_hash_table(st);
    /* read in ideal, move coefficients to integers */
    import_julia_data(bs_qq, bht, st, lens, exps, cfs);

    if (st->info_level > 0) {
        print_initial_statistics(stderr, st);
    }

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

    param_t **nmod_params =  (param_t **)calloc((unsigned long)st->nprimes,
            sizeof(param_t *));

    int *bad_primes = calloc((unsigned long)st->nprimes, sizeof(int));

    int dim = -1;
    long dquot_ori  = 0;

    srand(time(0));
    uint32_t prime = next_prime(1<<30);
    prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
    while(is_lucky_prime_ui(prime, bs_qq)){
      prime = next_prime(rand() % (1303905301 - (1<<30) + 1) + (1<<30));
    }
    uint32_t primeinit = prime;
    lp->p[0] = primeinit;

    sp_matfglm_t **bmatrix = (sp_matfglm_t **)calloc(st->nprimes,
            sizeof(sp_matfglm_t *));
    int32_t **bdiv_xn = (int32_t **)calloc(st->nprimes, sizeof(int32_t *));
    int32_t **blen_gb_xn = (int32_t **)calloc(st->nprimes, sizeof(int32_t *));
    int32_t **bstart_cf_gb_xn = (int32_t **)calloc(st->nprimes, sizeof(int32_t *));

    fglm_data_t **bdata_fglm = (fglm_data_t **)calloc(st->nprimes,
                                                      sizeof(fglm_data_t *));
    fglm_bms_data_t **bdata_bms = (fglm_bms_data_t **)calloc(st->nprimes,
            sizeof(fglm_bms_data_t *));
    int32_t *num_gb = (int32_t *)calloc(st->nprimes, sizeof(int32_t));
    int32_t **leadmons_ori = (int32_t **)calloc(st->nprimes, sizeof(int32_t *));
    int32_t **leadmons_current = (int32_t**)calloc(st->nprimes, sizeof(int32_t *));

    double stf4 = 0;
    uint64_t bsz = 0;
    long nlins = 0;
    uint64_t *linvars = calloc(bht->nv, sizeof(uint64_t));
    uint32_t **lineqs_ptr = malloc(sizeof(uint32_t *));
    uint64_t *squvars = calloc(nr_vars-1, sizeof(uint64_t));
    int success = 1;
    int squares = 1;
    int32_t *lmb_ori = modular_probabilistic_first(bmatrix, bdiv_xn, blen_gb_xn,
                                                   bstart_cf_gb_xn,

                                                   &nlins, linvars, lineqs_ptr,
                                                   squvars,

                                                   bdata_fglm, bdata_bms,

                                                   num_gb, leadmons_ori,

                                                   &bsz, nmod_params, bs_qq, bht,
                                                   st,
                                                   lp->p[0], //prime,
                                                   info_level,
                                                   print_gb,
                                                   &dim,
                                                   &dquot_ori,
                                                   gens,
                                                   files,
                                                   &success);

    for(int i = 0; i < nr_vars - 1; i++){
      if((squvars[i] == 0) && round == 1){
        squares = 0;
        success = 0;
      }
    }

    *dim_ptr = dim;
    *dquot_ptr = dquot_ori;

    mpz_param->dim   = dim;
    mpz_param->dquot = dquot_ori;
    if(lmb_ori == NULL || success == 0){
        if(*dim_ptr==1){
          if(info_level){
            fprintf(stderr, "Positive dimensional Grobner basis\n");
          }
          free_shared_hash_data(bht);
          free_hash_table(&bht);

          for (i = 0; i < st->nprimes; ++i) {
            //      free_basis(&(bs[i]));
          }
          free(bs);
          //here we should clean nmod_params
          free_lucky_primes(&lp);
          free(st);
          free(linvars);
          if(nlins){
            free(lineqs_ptr[0]);
          }
          free(lineqs_ptr);
          free(squvars);
          return 0;
        }
        if(*dquot_ptr==0){
            free_shared_hash_data(bht);
            free_hash_table(&bht);
            //sometimes this may crash (for instance of the ideal has positive dimension)
            for (i = 0; i < st->nprimes; ++i) {
                //      free_basis(&(bs[i]));
            }
            //      free(bs);
            //here we should clean nmod_params
            free_lucky_primes(&lp);
            free(st);

            free(linvars);
            /* free(lineqs_ptr[0]); */
            free(squvars);
            free(lineqs_ptr);
            return 0;
        }
        if(*dquot_ptr>0){
            /* print_msolve_message(stderr, 1); */
            free_shared_hash_data(bht);
            free_hash_table(&bht);
            //sometimes this may crash (for instance of the ideal has positive dimension)
            for (i = 0; i < st->nprimes; ++i) {
                //      free_basis(&(bs[i]));
            }
            free(bs);
            //here we should clean nmod_params
            free_lucky_primes(&lp);
            free(st);
            free(linvars);
            if(nlins){
                free(lineqs_ptr[0]); 
            }
            free(lineqs_ptr);
            free(squvars);
            if(squares == 0){
              return 2;
            }
            return 1;
        }
    }

    duplicate_data_mthread(st->nthrds, nr_vars, num_gb,
                           leadmons_ori, leadmons_current,
                           bdata_bms, bdata_fglm,
                           bstart_cf_gb_xn, blen_gb_xn, bdiv_xn,
                           bmatrix, nmod_params);

    if(info_level){
        fprintf(stderr, "\nStarts multi-modular computations based on probabilistic linear algebra\n");
    }

    for(int i = 0; i < st->nprimes; i++){
        normalize_nmod_param(nmod_params[i]);
    }

    mpz_param_t tmp_mpz_param;
    mpz_param_init(tmp_mpz_param);

    initialize_mpz_param(mpz_param, nmod_params[0]);
    initialize_mpz_param(tmp_mpz_param, nmod_params[0]);

    mpz_t modulus;
    mpz_init_set_ui(modulus, primeinit);

    mpq_t result, test;
    mpq_init(result);
    mpq_init(test);
    set_mpz_param_nmod(tmp_mpz_param, nmod_params[0]);
    long nsols = tmp_mpz_param->nsols;

    mpz_upoly_t numer;
    mpz_upoly_init2(numer, (nsols + 1), 32*(nsols + 1));
    numer->length = nsols + 1;

    mpz_upoly_t denom;
    mpz_upoly_init2(denom, (nsols + 1), 64*(nsols + 1));
    //  mpz_poly_init(denom, 33);
    denom->length = nsols + 1;

    mpz_t guessed_den;
    mpz_init2(guessed_den, 32*nsols);
    mpz_t guessed_num;
    mpz_init2(guessed_num, 32*nsols);

    long maxrec = 0;

    int rerun = 1, nprimes = 1, mcheck =1;
    long nbadprimes = 0;

    int *is_lifted = malloc(sizeof(int)*nr_vars);
    for(int i = 0; i < nr_vars; ++i){
      is_lifted[i] = 0;
    }

    prime = next_prime(1<<30);
    int br = 0;
    int doit = 1;
    rrec_data_t recdata;
    initialize_rrec_data(recdata);

    while(rerun == 1 || mcheck == 1){

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

      double ca0 = 0;
      if (nprimes==1) {
        ca0 = realtime();
      }

      modular_probabilistic_apply(bmatrix,
                                bdiv_xn,
                                blen_gb_xn,
                                bstart_cf_gb_xn,
                                nlins,
                                linvars,
                                lineqs_ptr[0],
                                squvars,
                                bdata_fglm,
                                bdata_bms,
                                num_gb,
                                leadmons_ori,
                                leadmons_current,
                                bsz,
                                nmod_params,
                                bs_qq,
                                bht,
                                st,
                                field_char,
                                0, //info_level,
                                bs,
                                lmb_ori,
                                dquot_ori,
                                lp,
                                gens,
                                &stf4,
                                nsols,
                                bad_primes);

      if(info_level){
        if(nprimes==1) {
          fprintf(stderr, "Probabilistic F4 + FGLM time (elapsed): %.2f sec\n",
                  realtime() - ca0);
        }
      }

      for(int i = 0; i < st->nprimes; i++){
        if (bad_primes[i] == 0) {
          normalize_nmod_param(nmod_params[i]);
        }
      }

      for(len_t i = 0; i < st->nprimes; i++){
        if(bad_primes[i] == 0){
          if(rerun == 0){
            mcheck = check_param_modular(mpz_param, nmod_params[i], lp->p[i],
                                         is_lifted, info_level);
          }

          if(mcheck==1){
            br = new_rational_reconstruction(mpz_param,
                                             tmp_mpz_param,
                                             nmod_params[i],
                                             numer,
                                             denom,
                                             &modulus, lp->p[i],
                                             &result,
                                             recdata,
                                             &guessed_num, &guessed_den,
                                             &maxrec,
                                             is_lifted,
                                             doit, st->nthrds,
                                             info_level);
          }


          if(br == 1){
            rerun = 0;
          }
          else{
            rerun = 1;
          }

          /* if(rational_reconstruction(mp_param, nmod_params[i], */
          /*                            numer, denom, */
          /*                            &modulus, lp->p[i], */
          /*                            &result, */
          /*                            &guessed_num, &guessed_den, */
          /*                            &maxrec, info_level)){ */
          /*   rerun = 0; */
          /* } */
          nprimes++;
        }
        else{
          if(info_level){
            fprintf(stderr, "<bp: %d>\n", lp->p[i]);
          }
          nbadprimes++;
          if(nbadprimes > nprimes){
            free(linvars);
            free(lineqs_ptr[0]);
            free(lineqs_ptr);
            free(squvars);
            free_rrec_data(recdata);
            return -4;
          }
        }
        if(info_level){
          if( (nprimes & (nprimes - 1)) == 0){
            fprintf(stderr, "{%d}", nprimes);
          }
        }
      }
    }

    mpz_param->denom->length = mpz_param->nsols;
    for(long i = 1; i <= mpz_param->nsols; i++){
      mpz_set(mpz_param->denom->coeffs[i-1], mpz_param->elim->coeffs[i]);
      mpz_mul_ui(mpz_param->denom->coeffs[i-1], mpz_param->denom->coeffs[i-1], i);
    }

    if(info_level){
      fprintf(stderr, "\n%d used primes\n", nprimes);
    }

    mpz_param_clear(tmp_mpz_param);
    mpz_upoly_clear(numer);
    mpz_upoly_clear(denom);
    mpz_clear(guessed_num);
    mpz_clear(guessed_den);
    mpq_clear(test);
    mpq_clear(result);
    mpz_clear(modulus);

    free_shared_hash_data(bht);
    free_hash_table(&bht);
    free_rrec_data(recdata);
    //sometimes this may crash (for instance of the ideal has positive dimension)
    for (i = 0; i < st->nprimes; ++i) {
        //      free_basis(&(bs[i]));
    }
    free(bs);
    //here we should clean nmod_params
    free_lucky_primes(&lp);
    free(st);
    free(bad_primes);

    for(i = 0; i < st->nprimes; ++i){
      free_fglm_bms_data(bdata_bms[i]);
      free_fglm_data(bdata_fglm[i]);
    }
    free(bdata_fglm);
    free(bdata_bms);

    free(linvars);
    free(lineqs_ptr[0]);
    free(lineqs_ptr);
    free(squvars);
    free(is_lifted);

    return 0;
}

int msolve_qq(mpz_param_t mp_param,
              param_t **nmod_param,
              int *dim_ptr,
              long *dquot_ptr,
              data_gens_ff_t *gens,
              int32_t ht_size, //initial_hts,
              int32_t nr_threads,
              int32_t max_nr_pairs,
              int32_t reset_ht,
              int32_t la_option,
              int32_t info_level,
              int32_t print_gb,
              int32_t pbm_file,
              files_gb *files,
              int round){

  if(la_option == 2 || la_option == 1){
    return msolve_trace_qq(mp_param,
                           nmod_param,
                           dim_ptr,
                           dquot_ptr,
                           gens,
                           ht_size, //initial_hts,
                           nr_threads,
                           max_nr_pairs,
                           reset_ht,
                           la_option,
                           info_level,
                           print_gb,
                           pbm_file,
                           files,
                           round);
  }
  else{
    return msolve_probabilistic_qq(mp_param,
                                   nmod_param,
                                   dim_ptr,
                                   dquot_ptr,
                                   gens,
                                   ht_size, //initial_hts,
                                   nr_threads,
                                   max_nr_pairs,
                                   reset_ht,
                                   la_option,
                                   info_level,
                                   print_gb,
                                   pbm_file,
                                   files,
                                   round);
  }
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
  for(long i = 0; i < nb - 1; i++){
    display_real_point_middle(fstream, pts[i]);
    fprintf(fstream, ", ");
  }
  display_real_point(fstream, pts[nb - 1]);
  fprintf(fstream, "\n");
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



int value_denom(mpz_t *denom, long deg, mpz_t r, long k,
                mpz_t *xdo, mpz_t *xup,
                mpz_t tmp, mpz_t den_do, mpz_t den_up,
                long corr){
  int boo = mpz_scalar_product_interval(denom, deg, k,
                                        xdo, xup,
                                        tmp, den_do, den_up, corr);
  if(mpz_cmp(den_do, den_up) > 0){
    fprintf(stderr, "BUG (den_do > den_up)\n");
    mpz_out_str(stderr, 10, den_do); fprintf(stderr, "\n");
    mpz_out_str(stderr, 10, den_up); fprintf(stderr, "\n");
    exit(1);
  }
  if(boo == 0){
    return boo;
  }
  mpz_t c;
  mpz_init(c);
  mpz_add_ui(c, r, 1);
  boo = mpz_poly_eval_interval(denom, deg, k,
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
  mpz_clear(c);

  return boo;
}

void lazy_single_real_root_param(mpz_param_t param, mpz_t *polelim,
                                 interval *rt, long nb, interval *pos_root,
                                 mpz_t *xdo, mpz_t *xup, mpz_t den_up, mpz_t den_do,
                                 mpz_t c, mpz_t tmp, mpz_t val_do, mpz_t val_up,
                                 mpz_t *tab, real_point_t pt,
                                 long prec, long nbits,
                                 int info_level){

  unsigned long ns = param->nsols ;
  if(rt->isexact==1){
    single_exact_real_root_param(param, rt, nb,
                                 xdo, xup, den_up, den_do,
                                 c, tmp, val_do, val_up,
                                 tab, pt, prec,
                                 info_level);
    return ;
  }

  long b = 16;
  prec = MAX(prec, rt->k);
  long corr = ns + rt->k;

  /* generate_table_values(rt, c, ns, b, */
  /*                       corr, */
  /*                       xdo, xup); */
  generate_table_values_full(rt, c, ns, b,
                             corr,
                             xdo, xup);

  while(/* lazy_mpz_poly_eval_interval(param->denom->coeffs, */
        /*                             param->denom->length - 1, */
        /*                             rt->k, */
        /*                             xdo, xup, */
        /*                             2*(prec+ns+rt->k), corr, b, */
        /*                             tmp, den_do, den_up) */
        value_denom(param->denom->coeffs,
                    param->denom->length - 1,
                    rt->numer,
                    rt->k,
                    xdo, xup,
                    tmp, den_do, den_up, corr)){
    if(mpz_sgn(rt->numer)>=0){
      get_values_at_bounds(param->elim->coeffs, ns, rt, tab);
      refine_QIR_positive_root(polelim, &ns, rt, tab, 2*(rt->k),
                               info_level);
    }
    else{
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
      refine_QIR_positive_root(polelim, &ns, pos_root, tab, 2*(pos_root->k),
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

    prec *= 2;
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


void real_roots_param(mpz_param_t param, interval *roots, long nb,
                      real_point_t *pts, long prec, long nbits,
                      double step,
                      int info_level){
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
  mpz_t *tab = (mpz_t*)(malloc(sizeof(mpz_t)*8));//table for some intermediate values
  for(int i=0;i<8;i++)mpz_init(tab[i]);

  mpz_t *polelim = malloc(sizeof(mpz_t) * param->elim->length);
  for(long i = 0; i < param->elim->length; i++){
    mpz_init_set(polelim[i], param->elim->coeffs[i]);
  }
  interval *pos_root = malloc(sizeof(interval));
  mpz_init(pos_root->numer);
  double et = realtime();
  for(long nc = 0; nc < nb; nc++){
    interval *rt = roots+nc;

    lazy_single_real_root_param(param, polelim, rt, nb, pos_root,
                                xdo, xup, den_up, den_do,
                                c, tmp, val_do, val_up, tab,
                                pts[nc], prec, nbits,
                                info_level);

    if(info_level){
      if(realtime() - et >= step){
        fprintf(stderr, "{%.2f%%}", 100*nc/((double) nb));
        et = realtime();
      }
    }

  }

  if(info_level){
    fprintf(stderr, "\n");
  }
  for(long i = 0; i < param->nsols; i++){
    mpz_clear(xup[i]);
    mpz_clear(xdo[i]);
  }
  free(xup);
  free(xdo);
  mpz_clear(c);
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
                   int32_t reset_ht,
                   int32_t la_option,
                   int32_t info_level,
                   int32_t print_gb,
                   int32_t pbm_file,
                   int32_t precision,
                   files_gb *files,
                   int round){
    if(la_option == 2 || la_option == 1){
        int b = msolve_trace_qq(mp_param,
                                nmod_param,
                                dim_ptr,
                                dquot_ptr,
                                gens,
                                ht_size, //initial_hts,
                                nr_threads,
                                max_nr_pairs,
                                reset_ht,
                                la_option,
                                info_level,
                                print_gb,
                                pbm_file,
                                files,
                                round);
        long unsigned int nbpos = 0;
        long unsigned int nbneg = 0;
        interval *roots   = NULL;
        real_point_t *pts = NULL;
        if(b==0 && *dim_ptr == 0 && *dquot_ptr > 0){

            mpz_t *pol = malloc(sizeof(mpz_t)*mp_param->elim->length);
            for(long i = 0; i < mp_param->elim->length; i++){
                mpz_init_set(pol[i], mp_param->elim->coeffs[i]);
            }
            long maxnbits = mpz_poly_max_bsize_coeffs(mp_param->elim->coeffs,
                                                   mp_param->elim->length - 1);
            long minnbits = mpz_poly_min_bsize_coeffs(mp_param->elim->coeffs,
                                                      mp_param->elim->length - 1);

            long prec = MAX(precision, 64 + 4*(LOG2(mp_param->elim->length) + LOG2(maxnbits)) );
            if(ilog2_mpz(mp_param->elim->coeffs[0]) > ilog2_mpz(mp_param->elim->coeffs[mp_param->nsols])){
              prec += (maxnbits - minnbits) / 2;
            }

            double st = realtime();
            roots = real_roots(pol, mp_param->elim->length - 1,
                               &nbpos, &nbneg, prec, nr_threads, info_level );
            long nb = nbpos + nbneg;
            double step = (realtime() - st) / (nb) * 10 * LOG2(precision);

            /* if(info_level){ */
            /*   display_roots_system(stderr, roots, nb); */
            /* } */
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
                    real_point_init(pts[i], mp_param->nvars);
                }

                real_roots_param(mp_param, roots, nb, pts, precision, maxnbits,
                                 step, info_level);
                if(info_level){
                    fprintf(stderr, "Elapsed time (real root extraction) = %.2f\n",
                            realtime() - st);
                }

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
            }

            for(long i = 0; i < mp_param->elim->length; i++){
                mpz_clear(pol[i]);
            }
            free(pol);
            *real_roots_ptr     = roots;
            *nb_real_roots_ptr  = nb;
            *real_pts_ptr       = pts;
        }

        return b;
    }
    else{
        int b = msolve_probabilistic_qq(mp_param,
                                        nmod_param,
                                        dim_ptr,
                                        dquot_ptr,
                                        gens,
                                        ht_size, //initial_hts,
                                        nr_threads,
                                        max_nr_pairs,
                                        reset_ht,
                                        la_option,
                                        info_level,
                                        print_gb,
                                        pbm_file,
                                        files,
                                        round);
        return b;
    }
    return 0;
}

int core_msolve(
  int32_t la_option,
  int32_t nr_threads,
  int32_t info_level,
  int32_t initial_hts,
  int32_t max_pairs,
  int32_t update_ht,
  int32_t generate_pbm,
  int32_t reduce_gb,
  int32_t print_gb,
  int32_t get_param,
  int32_t genericity_handling,
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
        if (normal_form == 0) {
            b = msolve_ff_alloc(&param, bld, blen, bexp, bcf,
                    gens, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, info_level, print_gb,
                    files);
            if (b == 0) {
              //When dquot = 1 
                if(files->out_file != NULL){
                    FILE *ofile = fopen(files->out_file, "a");
                    display_fglm_param_maple(ofile, param);
                    fclose(ofile);
                }
                else{
                    display_fglm_param_maple(stdout, param);
                }
            }
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
                fprintf(stderr, "\"-g2\" which is the default.\n");
            }
            if (b == 2) {
                fprintf(stderr, "The ideal has positive dimension\n");
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
            if(b == 3){
                if(files->out_file != NULL){
                    FILE *ofile2 = fopen(files->out_file, "a+");
                    fprintf(ofile2, "[0, %d, 0, [0, [1]]]:\n", gens->nvars);
                    fclose(ofile2);
                }
                else{
                    fprintf(stdout, "[0, %d, 0, [0, [1]]]:\n", gens->nvars);
                }
            }
            if (b == 2 || b == -1) {
                if(b==-1) b = 0;
            }
        } else {
            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            bs_t *tbr   = NULL;
            ht_t *bht   = NULL;
            stat_t *st  = NULL;

            /* generate array for storing multiplier for polynomial
             * to be reduced by basis */

            exp_t *mul  = (exp_t *)calloc(gens->nvars, sizeof(exp_t));

            /* for (int ii = 0; ii<gens->nvars; ++ii) {
             *     mul[ii] = 1;
             * } */

            int success = 0;

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
            success = initialize_f4_input_data(&bs, &bht, &st,
                    gens->lens, gens->exps, (void *)gens->cfs,
                    1073741827, 0 /* DRL order */, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens-normal_form, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, 1 /* reduce_gb */, 0,
                    info_level);

            st->fc  = gens->field_char;
            if(info_level){
              fprintf(stderr,
                      "NOTE: Field characteristic is now corrected to %u\n",
                      st->fc);
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
                success = core_f4(&bs, &bht, &st);

                if (!success) {
                    printf("Problem with F4, stopped computation.\n");
                    exit(1);
                }
            }
            /* initialize data for elements to be reduced,
             * NOTE: Don't initialize BEFORE running core_f4, bht may
             * change, so hash values of tbr may become wrong. */
            tbr = initialize_basis(2*normal_form);
            import_julia_data_nf_ff_32(
                    tbr, bht, st, gens->ngens-normal_form, gens->ngens,
                    gens->lens, gens->exps, (void *)gens->cfs);
            tbr->ld = tbr->lml  =  normal_form;
            /* normalize_initial_basis(tbr, st->fc); */
            for (int k = 0; k < normal_form; ++k) {
                tbr->lmps[k]  = k; /* fix input element in tbr */
            }

            /* compute normal form of last element in tbr */
            success = core_nf(&tbr, &bht, &st, mul, bs);

            if (!success) {
                printf("Problem with normalform, stopped computation.\n");
                exit(1);
            }
            /* print all reduced elements in tbr, first normal_form ones
             * are the input elements */
            print_msolve_polynomials_ff_32(
                    stdout, normal_form, tbr->lml, tbr, bht, st, gens->vnames, 0);

            if (normal_form_matrix > 0) {
            /* sht and hcm will store the union of the support
             * of all normal forms in tbr. */
            ht_t *sht   = initialize_secondary_hash_table(bht, st);
            hi_t *hcm   = (hi_t *)malloc(sizeof(hi_t));
            mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

            printf("\nStarts computation of normal form matrix\n");
            get_normal_form_matrix(tbr, bht, normal_form,
                    st, &sht, &hcm, &mat);

            printf("\n\nLength of union of support of all normal forms: %u\n",
                    mat->nc);

            printf("\nUnion of support, sorted by decreasing monomial order:\n");
            for (len_t k = 0; k < mat->nc; ++k) {
                for (len_t l = 0; l < sht->nv; ++l) {
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
    }
    else{
        int dim = - 2;
        long dquot = -1;
        b = real_msolve_qq(*mpz_paramp, 
                           &param,
                           &dim,
                           &dquot,
                           nb_real_roots_ptr,
                           real_roots_ptr,
                           real_pts_ptr,
                           gens,
                           initial_hts, nr_threads, max_pairs, update_ht,
                           la_option, info_level, print_gb,
                           generate_pbm, precision, files, round);

        if(b == 0){
            if(dim == 0 && dquot > 0){
                (*mpz_paramp)->nvars  = gens->nvars;
                if(files->out_file != NULL){
                    FILE *ofile = fopen(files->out_file, "a+");
                    if (get_param == 1) {
                        mpz_param_out_str_maple(ofile, gens, dquot, *mpz_paramp);
                    }
                    display_real_points_middle(
                            ofile, *real_pts_ptr, *nb_real_roots_ptr);
                    fclose(ofile);
                }
                else{
                    if (get_param == 1) {
                        mpz_param_out_str_maple(stdout, gens, dquot, *mpz_paramp);
                    }
                    display_real_points_middle(
                            stdout, *real_pts_ptr, *nb_real_roots_ptr);
                }
            }
            if(dquot == 0){

                if(files->out_file != NULL){
                    FILE *ofile2 = fopen(files->out_file, "a+");
                    if(get_param == 1){
                      fprintf(ofile2, "[0, %d, 0, [0, [1]]],", gens->nvars);
                    }
                    display_real_points_middle(
                                   stdout, *real_pts_ptr, *nb_real_roots_ptr);
                    fclose(ofile2);
                }
                else{
                  if(get_param == 1){
                    fprintf(stdout, "[0, %d, 0, [0, [1]]]:\n", gens->nvars);
                  }
                  display_real_points_middle(
                                   stdout, *real_pts_ptr, *nb_real_roots_ptr);
                }
            }
            if(dim > 0){
                fprintf(stderr, "The ideal has positive dimension\n");
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
            fprintf(stderr, "\"-g2\" which is the default.\n");
            (*mpz_paramp)->dim  = -1;
        }
        if(b==-2){
            fprintf(stderr, "Characteristic of the field here shouldn't be positive\n");
           (*mpz_paramp)->dim = -2;
        }
        if(b==-3){
            fprintf(stderr, "Problem when checking meta data\n");
            (*mpz_paramp)->dim  = -3;
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

static void export_julia_rational_parametrization_qq(
        int32_t *load,
        int32_t *dim,
        int32_t *dim_quot,
        int32_t **lens,
        void **cfs,
        const mpz_param_t param
        )
{
    int32_t i, j;
    int64_t ctr = 0;

    *load     = (int32_t)param->nvars+1;
    *dim      = (int32_t)param->dim;
    *dim_quot = (int32_t)param->dquot;

    if ((param->dim > 0) || (param->dim == 0 && param->dquot == 0)) {
        *lens = NULL;
        *cfs  = NULL;
    } else {
        int32_t *len  = (int32_t *)malloc(
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

        mpz_t *cf     = (mpz_t *)malloc(
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
        *lens = len;
        *cfs  = (void *)cf;
    }
}

void msolve_julia(
        int32_t *rp_ld,
        int32_t *rp_dim,
        int32_t *rp_dquot,
        int32_t **rp_lens,
        void **rp_cfs,
        int32_t *lens,
        int32_t *exps,
        void *cfs,
        char **var_names,
        char *output_file,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t initial_hts,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
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

    files_gb *files = calloc(1, sizeof(files_gb));

    if (output_file != NULL) {
        files->out_file = output_file;
    }

    data_gens_ff_t *gens = allocate_data_gens();

    gens->nvars                 = nr_vars;
    gens->ngens                 = nr_gens;
    gens->field_char            = field_char;
    gens->change_var_order      = -1;
    gens->linear_form_base_coef = 0;
    gens->vnames                = var_names;
    gens->lens                  = lens;
    gens->exps                  = exps;
    gens->rand_linear           = 0;

    if (field_char > 0) {
        gens->cfs = (int32_t *)cfs;
    } else {
        if (field_char == 0) {
            gens->mpz_cfs = (mpz_t **)cfs;
        }
    }

    /* data structures for parametrization */
    param_t *param;
    mpz_param_t mpz_param;
    mpz_param_init(mpz_param);

    long nb_real_roots      = 0;
    interval *real_roots    = NULL;
    real_point_t *real_pts  = NULL;

    /* main msolve functionality */
    core_msolve(la_option, nr_threads, info_level, initial_hts,
            max_nr_pairs, reset_ht, 0 /* generate pbm */, 1 /* reduce_gb */,
            print_gb, get_param, genericity_handling, 0 /* normal_form */,
            0 /* normal_form_matrix */, 0 /* is_gb */, precision,
            files, gens, &param, &mpz_param, &nb_real_roots,
            &real_roots, &real_pts);

    /* clean up data storage, but do not free data handled by julia */
    free(gens);
    gens  = NULL;

    export_julia_rational_parametrization_qq(
            rp_ld, rp_dim, rp_dquot, rp_lens, rp_cfs, mpz_param);
    /* free parametrization */
    free(param);
    mpz_param_clear(mpz_param);

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
