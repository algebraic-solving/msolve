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

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(x, y) ((x) > (y) ? (y) : (x))
#endif
#define LOG2(X)                                                                \
  ((unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((X)) - 1))
#define ilog2_mpz(a) mpz_sizeinbase(a, 2)

static void mpz_upoly_init(mpz_upoly_t poly, deg_t alloc) {
  mpz_t *tmp = NULL;
  if (alloc) {
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    if (tmp == NULL) {
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for (deg_t i = 0; i < alloc; i++) {
      mpz_init(tmp[i]);
      mpz_set_si(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

static void mpz_upoly_init2(mpz_upoly_t poly, deg_t alloc, bits_t nbits) {
  mpz_t *tmp = NULL;
  if (alloc) {
    tmp = (mpz_t *)(malloc(alloc * sizeof(mpz_t)));
    if (tmp == NULL) {
      fprintf(stderr, "Unable to allocate in mpz_upoly_init\n");
      exit(1);
    }
    for (deg_t i = 0; i < alloc; i++) {
      mpz_init2(tmp[i], nbits);
      mpz_set_si(tmp[i], 0);
    }
  }
  poly->coeffs = tmp;
  poly->alloc = alloc;
  poly->length = -1;
}

void mpz_upoly_clear(mpz_upoly_t pol) {

  for (long i = 0; i < pol->alloc; i++) {
    mpz_clear(pol->coeffs[i]);
  }
  free(pol->coeffs);
}

static inline void mpz_upoly_out_str(FILE *file, mpz_upoly_t pol) {
  fprintf(file, "[");
  if (pol->length > 0) {
    fprintf(file, "%d, ", pol->length - 1); // degree
    fprintf(file, "[");
    for (long i = 0; i < pol->length - 1; i++) {
      mpz_out_str(file, 10, pol->coeffs[i]);
      fprintf(file, ", ");
    }
    mpz_out_str(file, 10, pol->coeffs[pol->length - 1]);
    fprintf(file, "]");
  } else {
    fprintf(file, "-1, [0]");
  }
  fprintf(file, "]");
}

void mpz_param_init(mpz_param_t param) {
  param->nvars = 0;
  param->nsols = 0;
  mpz_upoly_init(param->elim, 0);
  mpz_upoly_init(param->denom, 0);
  param->coords = NULL;
  param->cfs = NULL;
}

void mpz_param_clear(mpz_param_t param) {
  mpz_upoly_clear(param->elim);
  mpz_upoly_clear(param->denom);
  if (param->coords != NULL) {
    for (long i = 0; i < param->nvars - 1; i++) {
      mpz_upoly_clear(param->coords[i]);
      mpz_clear(param->cfs[i]);
    }
  }
  free(param->coords);
  free(param->cfs);
  param->nvars = 0;
  param->nsols = 0;
}

static inline void mpz_param_out_str(FILE *file, const data_gens_ff_t *gens,
                                     const long dquot, mpz_param_t param,
                                     param_t *mod_param) {
  fprintf(file, "[");
  fprintf(file, "%d, \n", gens->field_char); /* field charac */
  fprintf(file, "%d, \n", param->nvars);     // nvars
  fprintf(file, "%ld, \n", dquot);           // dim quotient
  /* Print all variables:
   * 1. may include new variable from added linear form,
   * 2. variables may have a different order than in the.
   *    input system
   * Both of the above points are due to genericity handling. */
  fprintf(file, "[");
  for (int i = 0; i < param->nvars - 1; ++i) {
    fprintf(file, "'%s', ", gens->vnames[i]);
  }
  fprintf(file, "'%s'],\n", gens->vnames[param->nvars - 1]);

  fprintf(file, "[");
  if (gens->rand_linear) {
    int32_t sum = 0;
    if (gens->field_char == 0) {
      for (int i = 0; i < param->nvars; i++) {
        sum += abs(gens->random_linear_form[i]) * param->nvars - 1;
      }
    }
    for (int i = 0; i < param->nvars - 1; i++) {
      fprintf(file, "%d", gens->random_linear_form[i]);
      if (gens->field_char == 0) {
        fprintf(file, "/%d", sum);
      }
      fprintf(file, ",");
    }
    fprintf(file, "%d", gens->random_linear_form[param->nvars - 1]);
    if (gens->field_char == 0) {
      fprintf(file, "/%d", sum);
    }
  } else {
    if (gens->linear_form_base_coef > 0) {
      const int32_t bcf = gens->linear_form_base_coef;
      for (int i = 0; i < param->nvars - 1; ++i) {
        fprintf(file, "%d,", (int32_t)(pow(i + 1, bcf - 1)));
      }
      fprintf(file, "%d", (int32_t)(1));
    } else {
      for (int i = 0; i < param->nvars - 1; ++i) {
        fprintf(file, "%d, ", (int32_t)(0));
      }
      fprintf(file, "%d", (int32_t)(1));
    }
  }
  fprintf(file, "],\n");
  fprintf(file, "[1,\n["); /*at the moment, a single param is returned */
  if (gens->field_char) {
    display_nmod_poly(file, mod_param->elim);
  } else {
    mpz_upoly_out_str(file, param->elim); // elim. poly
  }
  fprintf(file, ",\n");
  if (gens->field_char) {
    display_nmod_poly(file, mod_param->denom);
  } else {
    mpz_upoly_out_str(file, param->denom); // denom. poly
  }
  fprintf(file, ",\n");
  fprintf(file, "[\n");
  if (gens->field_char) { /* positive characteristic */
    if (mod_param->coords != NULL) {
      for (int i = 0; i < mod_param->nvars - 1; i++) {
        fprintf(file, "[");
        if (gens->field_char) {
          display_nmod_poly(file, mod_param->coords[i]);
        }
        if (i == mod_param->nvars - 2) {
          fprintf(file, "]\n");
        } else {
          fprintf(file, "],\n");
        }
      }
    }
  } else {
    if (param->coords != NULL) {
      for (int i = 0; i < param->nvars - 1; i++) {
        fprintf(file, "[");
        mpz_upoly_out_str(file, param->coords[i]); // param. polys
        fprintf(file, ",\n");
        mpz_out_str(file, 10, param->cfs[i]);
        if (i == param->nvars - 2) {
          fprintf(file, "]\n");
        } else {
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
                                           const data_gens_ff_t *gens,
                                           const long dquot, mpz_param_t param,
                                           param_t *mod_param) {
  mpz_param_out_str(file, gens, dquot, param, mod_param);
  fprintf(file, "]");
}

static inline void display_fglm_crt_matrix(FILE *file, crt_mpz_matfglm_t mat) {

  fprintf(file, "%u\n", mat->ncols);
  fprintf(file, "%u\n", mat->nrows);

  long len1 = (mat->ncols) * (mat->nrows);
  for (long i = 0; i < len1; i++) {
    mpz_out_str(file, 10, mat->dense_mat[i]);
    fprintf(file, " ");
  }
  fprintf(file, "\n");
  long len2 = (mat->ncols) - (mat->nrows);
  for (long i = 0; i < len2; i++) {
    fprintf(file, "%d ", mat->triv_idx[i]);
  }
  fprintf(file, "\n");
  for (long i = 0; i < len2; i++) {
    fprintf(file, "%d ", mat->triv_pos[i]);
  }
  fprintf(file, "\n");
  for (long i = 0; i < mat->nrows; i++) {
    fprintf(file, "%d ", mat->dense_idx[i]);
  }
  fprintf(file, "\n");
}

static inline void display_fglm_mpq_matrix(FILE *file, mpq_matfglm_t mat) {

  fprintf(file, "%u\n", mat->ncols);
  fprintf(file, "%u\n", mat->nrows);

  uint64_t nc = mat->ncols;

  fprintf(file, "[");
  for (long i = 0; i < mat->nrows; i++) {
    long c = 2 * i * mat->ncols;
    fprintf(file, "[");
    for (long j = 0; j < nc - 1; j++) {
      mpz_out_str(file, 10, mat->dense_mat[c + 2 * j]);
      fprintf(file, "/");
      mpz_out_str(file, 10, mat->dense_mat[c + 2 * j + 1]);
      fprintf(file, ", ");
    }
    mpz_out_str(file, 10, mat->dense_mat[c + 2 * nc - 2]);
    fprintf(file, "/");
    mpz_out_str(file, 10, mat->dense_mat[c + 2 * nc - 1]);
    fprintf(file, "]\n");
  }
  fprintf(file, "]");
  fprintf(file, "\n");
  long len2 = (mat->ncols) - (mat->nrows);
  for (long i = 0; i < len2; i++) {
    fprintf(file, "%d ", mat->triv_idx[i]);
  }
  fprintf(file, "\n");
  for (long i = 0; i < len2; i++) {
    fprintf(file, "%d ", mat->triv_pos[i]);
  }
  fprintf(file, "\n");
  for (long i = 0; i < mat->nrows; i++) {
    fprintf(file, "%d ", mat->dense_idx[i]);
  }
  fprintf(file, "\n");
}

data_gens_ff_t *allocate_data_gens() {
  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
  gens->lens = NULL;
  gens->exps = NULL;
  gens->cfs = NULL;
  gens->mpz_cfs = NULL;
  gens->random_linear_form = NULL;

  gens->elim = 0;
  return gens;
}

void free_data_gens(data_gens_ff_t *gens) {
  for (long i = 0; i < gens->nvars; i++) {
    free(gens->vnames[i]);
  }
  free(gens->vnames);
  if (gens->field_char == 0) {
    for (long i = 0; i < 2 * gens->nterms; i++) {
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
                                int32_t **blen, int32_t *bcf, int32_t **bexp) {
  if (bcf[ind] != 0 && bcf[ind] != 1) {
    fprintf(file, "%d", bcf[ind]);
    display_monomial(file, gens, ind, bexp);
  } else {
    if (bcf[ind] == 1) {
      int32_t b = display_monomial_single(file, gens, ind, bexp);
      if (b == 0) {
        fprintf(file, "1");
      }
    }
  }
}

static inline void display_basis(FILE *file, int32_t *bld, int32_t **blen,
                                 int32_t **bexp, int32_t *bcf,
                                 data_gens_ff_t *gens) {
  int32_t npos = 0;
  for (int32_t i = 0; i < bld[0]; i++) {
    int32_t j, end = (*blen)[i];
    for (j = 0; j < end; j++) {
      display_term(file, npos + j, gens, blen, bcf, bexp);
      if (j < (*blen)[i] - 1) {
        fprintf(file, "+");
      }
    }
    //    display_term(file, npos+j, gens, blen, bcf, bexp);
    if (i < bld[0] - 1) {
      fprintf(file, ",\n");
    } else {
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
}

static inline void display_basis_maple(FILE *file, int32_t *bld, int32_t **blen,
                                       int32_t **bexp, int32_t *bcf,
                                       data_gens_ff_t *gens) {
  int32_t npos = 0;
  fprintf(file, "[");
  for (int32_t i = 0; i < bld[0]; i++) {
    int32_t j, end = (*blen)[i];
    for (j = 0; j < end; j++) {
      display_term(file, npos + j, gens, blen, bcf, bexp);
      if (j < (*blen)[i] - 1) {
        fprintf(file, "+");
      }
    }
    //    display_term(file, npos+j, gens, blen, bcf, bexp);
    if (i < bld[0] - 1) {
      fprintf(file, ",\n");
    } else {
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
  fprintf(file, "]:\n");
}

static inline void display_lead_monomials_from_gb(FILE *file, int32_t *bld,
                                                  int32_t **blen,
                                                  int32_t **bexp, int32_t *bcf,
                                                  data_gens_ff_t *gens) {
  int32_t npos = 0;
  for (int32_t i = 0; i < bld[0]; i++) {
    display_term(file, npos, gens, blen, bcf, bexp);
    if (i < bld[0] - 1) {
      fprintf(file, ",\n");
    } else {
      fprintf(file, "\n");
    }
    npos += (*blen)[i];
  }
}

static inline void display_monomials_from_array(FILE *file, long length,
                                                int32_t *bexp, const int nv,
                                                char **vnames) {
  fprintf(file, "[");
  for (long i = 0; i < length - 1; i++) {
    display_monomial_full(file, nv, vnames, i, bexp);
    fprintf(file, ", ");
  }
  display_monomial_full(file, nv, vnames, length - 1, bexp);
  fprintf(file, "]");
}

static inline void display_monomials_from_array_maple(FILE *file, long length,
                                                      int32_t *bexp,
                                                      const int nv,
                                                      char **vnames) {
  fprintf(file, "[");
  for (long i = 0; i < length - 1; i++) {
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
static int undo_variable_order_change(data_gens_ff_t *gens) {
  int32_t i, j;
  int32_t len, tmp;
  char *tmp_char = NULL;
  const int32_t cvo = gens->change_var_order;
  const int32_t nvars = gens->nvars;
  const int32_t ngens = gens->ngens;

  if (gens->linear_form_base_coef > 0) {
    return 0;
  }
  if (cvo > -1) {
    /* undo last variable change */
    tmp_char = gens->vnames[nvars - 1];
    gens->vnames[nvars - 1] = gens->vnames[cvo];
    gens->vnames[cvo] = tmp_char;
    len = 0;
    tmp = 0;
    for (i = 0; i < ngens; ++i) {
      for (j = 0; j < gens->lens[i]; ++j) {
        tmp = gens->exps[len + j * nvars + nvars - 1];
        gens->exps[len + j * nvars + nvars - 1] =
            gens->exps[len + j * nvars + cvo];
        gens->exps[len + j * nvars + cvo] = tmp;
      }
      len += gens->lens[i] * nvars;
    }
  }
  /* all cyclic changes already done, stop here, try to add
   * a linear form with additional variable afterwards */
  gens->change_var_order++;
  if (gens->change_var_order == nvars - 1) {
    return 0;
  } else {
    return 1;
  }
}
static int change_variable_order_in_input_system(data_gens_ff_t *gens,
                                                 int32_t info_level) {
  int32_t i, j;
  int32_t len, tmp;
  char *tmp_char = NULL;
  const int32_t cvo = gens->change_var_order;
  const int32_t nvars = gens->nvars;

  if (undo_variable_order_change(gens) == 0) {
    return 0;
  }
  /* do the current variable change */
  tmp_char = gens->vnames[nvars - 1];
  gens->vnames[nvars - 1] = gens->vnames[cvo + 1];
  gens->vnames[cvo + 1] = tmp_char;
  len = 0;
  tmp = 0;
  for (i = 0; i < gens->ngens; ++i) {
    for (j = 0; j < gens->lens[i]; ++j) {
      tmp = gens->exps[len + j * nvars + nvars - 1];
      gens->exps[len + j * nvars + nvars - 1] =
          gens->exps[len + j * nvars + cvo + 1];
      gens->exps[len + j * nvars + cvo + 1] = tmp;
    }
    len += gens->lens[i] * nvars;
  }
  if (info_level > 0) {
    printf("\nChanging variable order for possibly more generic staircase:\n");
    for (int i = 0; i < nvars - 1; i++) {
      fprintf(stdout, "%s, ", gens->vnames[i]);
    }
    fprintf(stdout, "%s\n", gens->vnames[nvars - 1]);
  }
  return 1;
}

/* Used if staircase is not generic enough, after
 * we tried to change variable order.
 *
 * Returns 1 if addition of linear form is done,
 * 0 if an error occurs. */
static int add_linear_form_to_input_system(data_gens_ff_t *gens,
                                           int32_t info_level) {
  int64_t i, j;
  int32_t k;
  int32_t nvars_old, nvars_new;
  int64_t len_old = 0, len_new;
  if (gens->linear_form_base_coef == 0) {
    nvars_old = gens->nvars;
    nvars_new = nvars_old + 1;
    for (i = 0; i < gens->ngens; ++i) {
      len_old += gens->lens[i];
    }
    len_new = len_old + nvars_new;
  } else {
    nvars_old = gens->nvars - 1;
    nvars_new = nvars_old + 1;
    for (i = 0; i < gens->ngens - 1; ++i) {
      len_old += gens->lens[i];
    }
    len_new = len_old + gens->lens[gens->ngens - 1];
  }
  /* for the first run we have to reinitialize gens data with the newly
   * added variable, then we add the linear form at the end. */
  if (gens->linear_form_base_coef == 0) {
    char *extra_var = (char *)malloc(2 * sizeof(char));
    strcpy(extra_var, "A");
    gens->nvars++;
    gens->ngens++;
    gens->lens =
        realloc(gens->lens, (unsigned long)gens->ngens * sizeof(int32_t));
    gens->lens[gens->ngens - 1] = (int32_t)nvars_new;
    /* add new dummy variable (only need for correct freeing of vnames at
     * the end */
    gens->vnames =
        realloc(gens->vnames, (unsigned long)gens->nvars * sizeof(char *));
    gens->vnames[gens->nvars - 1] = (char *)malloc(2 * sizeof(char));
    gens->vnames[gens->nvars - 1] = extra_var;
    int32_t *old_exps = gens->exps;
    gens->exps =
        (int32_t *)calloc((unsigned long)len_new * nvars_new, sizeof(int32_t));
    i = 0;
    j = 0;
    while (i < (len_old * nvars_old)) {
      memcpy(gens->exps + j, old_exps + i,
             (unsigned long)nvars_old * sizeof(int32_t));
      i += nvars_old;
      j += nvars_new;
    }
    free(old_exps);
    /* add linear form exponents */
    while (j < (len_new * nvars_new)) {
      gens->exps[j] = 1;
      j += nvars_new + 1;
    }
    /* allocate memory for coefficients of linear form */
    if (gens->field_char > 0) {
      gens->cfs = realloc(gens->cfs, (unsigned long)len_new * sizeof(int32_t));
    } else {
      gens->mpz_cfs = realloc(gens->mpz_cfs,
                              (unsigned long)2 * len_new * (sizeof(mpz_t *)));
      for (i = 2 * len_old; i < 2 * len_new; i += 2) {
        gens->mpz_cfs[i] = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(gens->mpz_cfs[i]));
        /* set denominators 1 for all coefficients */
        mpz_set_ui(*(gens->mpz_cfs[i]), 1);
        gens->mpz_cfs[i + 1] = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(gens->mpz_cfs[i + 1]));
        /* set denominators 1 for all coefficients */
        mpz_set_ui(*(gens->mpz_cfs[i + 1]), 1);
      }
    }
  }
  gens->linear_form_base_coef++;
  const int32_t bcf = gens->linear_form_base_coef;
  k = 1;
  if (info_level > 0) {
    printf("\nAdding a linear form with an extra variable ");
    printf("(lowest w.r.t. monomial order)\n");
    printf("[coefficients of linear form are k^%d for k looping over variable "
           "index 1...n]\n",
           bcf - 1);
  }
  if (gens->field_char > 0) {
    for (i = len_old; i < len_new - 1; ++i) {
      gens->cfs[i] = ((int32_t)(pow(k, bcf - 1)) % gens->field_char);
      k++;
    }
    gens->cfs[len_new - 1] = 1;
    k++;
  } else {
    for (i = 2 * len_old; i < 2 * len_new; i += 2) {
      mpz_set_ui(*(gens->mpz_cfs[i]), (int32_t)(pow(k, bcf - 1)));
      k++;
    }
    mpz_set_si(*(gens->mpz_cfs[2 * (len_new - 1)]), 1);
  }
  return 1;
}

static int add_random_linear_form_to_input_system(data_gens_ff_t *gens,
                                                  int32_t info_level) {

  int64_t i, j;
  int32_t nvars_old, nvars_new;
  int64_t len_old = 0, len_new;

  if (gens->linear_form_base_coef == 0) {
    nvars_old = gens->nvars;
    nvars_new = nvars_old + 1;
    for (i = 0; i < gens->ngens; ++i) {
      len_old += gens->lens[i];
    }
    len_new = len_old + nvars_new;
  } else {
    nvars_old = gens->nvars - 1;
    nvars_new = nvars_old + 1;
    for (i = 0; i < gens->ngens - 1; ++i) {
      len_old += gens->lens[i];
    }
    len_new = len_old + gens->lens[gens->ngens - 1];
  }

  /* for the first run we have to reinitialize gens data with the newly
   * added variable, then we add the linear form at the end. */
  if (gens->linear_form_base_coef == 0) {
    char *extra_var = (char *)malloc(2 * sizeof(char));
    strcpy(extra_var, "A");
    gens->nvars++;
    gens->ngens++;
    gens->lens =
        realloc(gens->lens, (unsigned long)gens->ngens * sizeof(int32_t));
    gens->lens[gens->ngens - 1] = (int32_t)nvars_new;
    /* add new dummy variable (only need for correct freeing of vnames at
     * the end */
    gens->vnames =
        realloc(gens->vnames, (unsigned long)gens->nvars * sizeof(char *));
    gens->vnames[gens->nvars - 1] = (char *)malloc(2 * sizeof(char));
    gens->vnames[gens->nvars - 1] = extra_var;
    int32_t *old_exps = gens->exps;
    gens->exps =
        (int32_t *)calloc((unsigned long)len_new * nvars_new, sizeof(int32_t));
    i = 0;
    j = 0;
    while (i < (len_old * nvars_old)) {
      memcpy(gens->exps + j, old_exps + i,
             (unsigned long)nvars_old * sizeof(int32_t));
      i += nvars_old;
      j += nvars_new;
    }
    free(old_exps);
    /* add linear form exponents */
    while (j < (len_new * nvars_new)) {
      gens->exps[j] = 1;
      j += nvars_new + 1;
    }
    /* allocate memory for coefficients of linear form */
    if (gens->field_char > 0) {
      gens->cfs = realloc(gens->cfs, (unsigned long)len_new * sizeof(int32_t));
    } else {
      gens->mpz_cfs = realloc(gens->mpz_cfs,
                              (unsigned long)2 * len_new * (sizeof(mpz_t *)));
      for (i = 2 * len_old; i < 2 * len_new; i += 2) {
        gens->mpz_cfs[i] = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(gens->mpz_cfs[i]));
        gens->mpz_cfs[i + 1] = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(gens->mpz_cfs[i + 1]));
        /* set denominators 1 for all coefficients */
        mpz_set_ui(*(gens->mpz_cfs[i + 1]), 1);
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
  gens->random_linear_form =
      realloc(gens->random_linear_form, sizeof(int32_t) * (nvars_new));

  if (gens->field_char > 0) {
    int j = 0;
    for (i = len_old; i < len_new; ++i) {
      gens->random_linear_form[j] = ((int8_t)(rand()) % gens->field_char);

      while (gens->random_linear_form[j] == 0) {
        gens->random_linear_form[j] = ((int8_t)(rand()) % gens->field_char);
      }
      gens->cfs[i] = gens->random_linear_form[j];
      j++;
    }
  } else {
    int j = 0;
    int32_t sum = 0;
    for (i = 2 * len_old; i < 2 * len_new; i += 2) {
      gens->random_linear_form[j] = ((int8_t)(rand()));

      while (gens->random_linear_form[j] == 0) {
        gens->random_linear_form[j] = ((int8_t)(rand()));
      }
      if (i < 2 * len_new - 1) {
        sum += nvars_old * abs(gens->random_linear_form[j]);
      } else {
        gens->random_linear_form[j] = sum;
      }
      mpz_set_si(*(gens->mpz_cfs[i]), gens->random_linear_form[j]);
      mpz_set_ui(*(gens->mpz_cfs[i + 1]), 1);

      j++;
    }
  }
  gens->rand_linear = 1;
  return 1;
}

static inline void nmod_param_out_str(FILE *file, const long dquot,
                                      param_t *param) {
  fprintf(file, "[");
  fprintf(file, "0, \n");                // charac
  fprintf(file, "%d, \n", param->nvars); // nvars
  fprintf(file, "%ld, \n", dquot);       // dim quotient
  display_nmod_poly(file, param->elim);  // elim. poly
  fprintf(file, ",\n");
  display_nmod_poly(file, param->denom); // elim. poly
  fprintf(file, ",\n");
  fprintf(file, "[\n");
  for (int i = 0; i < param->nvars - 1; i++) {
    fprintf(file, "[");
    display_nmod_poly(file, param->coords[i]); // elim. poly
    fprintf(file, ",\n");
    fprintf(file, "[1]");
    if (i == param->nvars - 2) {
      fprintf(file, "]\n");
    } else {
      fprintf(file, "],\n");
    }
  }
  fprintf(file, "]");
  fprintf(file, "]");
}

static inline void nmod_param_out_str_maple(FILE *file, const long dquot,
                                            param_t *param) {
  nmod_param_out_str(file, dquot, param);
  fprintf(file, ":\n");
}

static inline void print_msolve_message(FILE *file, int n) {
  if (n == 0) {
    fprintf(file, "The ideal has positive dimension.\n");
  }
  if (n == 1) {
    fprintf(file,
            "The staircase of the Grobner basis is not generic enough.\n");
    fprintf(file, "\nYou may add to your system a random linear combination of "
                  "\nthe variables plus a new variable which you choose to be "
                  "\nthe smallest for the monomial ordering.\n\n");
  }
}

static inline void initialize_mpz_param(mpz_param_t param, param_t *bparam) {

  param->nvars = bparam->nvars;
  param->nsols = bparam->elim->length - 1;
  mpz_upoly_init2(param->elim, bparam->elim->alloc,
                  2 * 32 * (bparam->elim->length));
  mpz_upoly_init2(param->denom, bparam->elim->alloc - 1,
                  2 * 32 * (bparam->elim->length));
  param->elim->length = bparam->elim->length;

  param->coords =
      (mpz_upoly_t *)malloc(sizeof(mpz_upoly_t) * (param->nvars - 1));
  if (param->coords != NULL) {
    for (len_t i = 0; i < param->nvars - 1; i++) {
      mpz_upoly_init(param->coords[i], MAX(1, bparam->elim->alloc - 1));
      param->coords[i]->length = bparam->elim->length - 1;
    }
  } else {
    fprintf(stderr, "Error when initializing parametrization\n");
    exit(1);
  }
  param->cfs = (mpz_t *)malloc(sizeof(mpz_t) * (param->nvars - 1));
  if (param->cfs != NULL) {
    for (len_t i = 0; i < param->nvars - 1; i++) {
      mpz_init(param->cfs[i]);
      mpz_set_ui(param->cfs[i], 1);
    }
  } else {
    fprintf(stderr, "Error when allocating cfs\n");
    exit(1);
  }
}

static inline void reduce_generators(mpz_t **tmp, long nterms, int32_t *cfs,
                                     int32_t prime) {
  for (int32_t i = 0; i < 2 * nterms; i += 2) {
    cfs[i / 2] = (int32_t)mpz_fdiv_ui(tmp[i][0], prime);
  }
}

static inline void
check_and_set_linear_poly_non_hashed(long *nlins_ptr, nvars_t *linvars,
                                     uint32_t **lineqs_ptr, int32_t *bld,
                                     int32_t *bexp_lm, int32_t **blen,
                                     int32_t **bexp, int32_t *bcf, long nvars) {
  long nlins = 0;
  uint32_t *coefpos = calloc(nvars, sizeof(uint32_t));
  uint32_t pos = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1.
  */
  for (long i = 0; i < bld[0]; i++) {
    long deg = 0;
    for (int j = 0; j < nvars; j++) {
      deg += bexp_lm[i * nvars + j];
    }
    if (deg == 1) {
      nlins++;
      for (int k = 0; k < nvars; k++) {
        if (bexp_lm[i * nvars + k] == 1) {
          linvars[k] = i + 1;
          coefpos[k] = pos;
        }
      }
    }
    pos += (*blen)[i];
  }
  *nlins_ptr = nlins;

  // On recupere les coefficients des formes lineaires
  uint32_t *lineqs = calloc(nlins * (nvars + 1), sizeof(uint32_t));

  int cnt = 0;

  for (int i = 0; i < nvars; i++) {
    if (linvars[i] != 0) {

      long len = (*blen)[linvars[i] - 1];

      if (len == nvars + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = (uint32_t)((bcf)[j + coefpos[i]]);
          lineqs[cnt * (nvars + 1) + j] = coef;
        }
      } else {
        //        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = (uint32_t)((bcf)[j + coefpos[i]]);
          // bs->cf_32[bs->hm[bi][COEFFS]][j];

          int32_t *exp = (*bexp) + (j + coefpos[i]) * nvars; // bht->ev[dt[j]];
          int isvar = 0;
          for (int k = 0; k < nvars; k++) {
            if (exp[k] == 1) {
              lineqs[cnt * (nvars + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (nvars + 1) + nvars] = coef;
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
check_and_set_vars_squared_in_monomial_basis(nvars_t *squvars, int32_t *lmb,
                                             long dquot, long nvars) {
  /*
    la i-ieme entree de squvars est a 0 si le monome vars[i]^2 n'est
    pas dans la base monomiale.
    Sinon, on met l'indice du monome. De toute facon, en 0, on aurait
    le monome 1.
  */
  for (deg_t i = 0; i < dquot; i++) {
    deg_t deg = 0;
    for (nvars_t j = 0; j < nvars; j++) {
      deg += lmb[i * nvars + j];
    }
    if (deg == 2) {
      for (nvars_t j = 0; j < nvars - 1; j++) {
        if (lmb[i * nvars + j] == 2) {
          squvars[j] = i;
          break;
        }
      }
    }
  }
}

static inline void normalize_nmod_param(param_t *nmod_param) {
  if (nmod_param != NULL) {
    uint32_t prime = nmod_param->charac;
    uint32_t inv;
    ulong inv2;
    n_gcdinv(&inv2, nmod_param->elim->length - 1, prime);
    inv = inv2 % prime;
    nmod_poly_fit_length(nmod_param->denom, nmod_param->elim->length - 1);
    nmod_param->denom->length = nmod_param->elim->length - 1;

    for (long i = 1; i < nmod_param->elim->length; i++) {
      nmod_param->denom->coeffs[i - 1] =
          (i * nmod_param->elim->coeffs[i]) % prime;
    }

    for (long i = 0; i < nmod_param->elim->length - 1; i++) {
      nmod_param->denom->coeffs[i] =
          (inv * nmod_param->denom->coeffs[i]) % prime;
    }

    for (int j = 0; j < nmod_param->nvars - 1; j++) {
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
static inline void set_mpz_param_nmod(mpz_param_t mpz_param,
                                      param_t *nmod_param) {
  mpz_param->elim->length = nmod_param->elim->length;
  for (long i = 0; i <= mpz_param->nsols; i++) {
    mpz_set_ui(mpz_param->elim->coeffs[i], (nmod_param)->elim->coeffs[i]);
  }
  mpz_param->elim->length = nmod_param->elim->length;

  for (long i = 0; i < nmod_param->denom->length; i++) {
    mpz_set_ui(mpz_param->denom->coeffs[i], nmod_param->denom->coeffs[i]);
  }
  mpz_param->denom->length = nmod_param->denom->length;

  for (int j = 0; j < mpz_param->nvars - 1; j++) {

    for (long i = 0; i < nmod_param->coords[j]->length; i++) {
      mpz_set_ui(mpz_param->coords[j]->coeffs[i],
                 nmod_param->coords[j]->coeffs[i]);
    }
    for (long i = nmod_param->coords[j]->length;
         i < nmod_param->elim->length - 1; i++) {
      mpz_set_ui(mpz_param->coords[j]->coeffs[i], 0);
    }
    mpz_param->coords[j]->length = nmod_param->elim->length - 1;
  }
}

static inline void crt_lift_mpz_upoly(mpz_upoly_t pol, nmod_poly_t nmod_pol,
                                      mpz_t modulus, int32_t prime, mpz_t prod,
                                      mpz_t tmp, int nthrds) {
  len_t i;
  /* #pragma omp parallel for num_threads(nthrds)    \ */
  /*   private(i) schedule(static) */
  for (i = 0; i < pol->length; i++) {
    mpz_CRT_ui(pol->coeffs[i], pol->coeffs[i], modulus, nmod_pol->coeffs[i],
               prime, prod, tmp, 0);
  }
}

/* assumes that all degrees are the same */
static inline void crt_lift_mpz_param(mpz_param_t mpz_param,
                                      param_t *nmod_param, mpz_t modulus,
                                      mpz_t prod_crt, const int32_t prime,
                                      mpz_t tmp, const int nthrds) {

  /*assumes prod_crt = modulus * prime */
  crt_lift_mpz_upoly(mpz_param->elim, nmod_param->elim, modulus, prime,
                     prod_crt, tmp, nthrds);
  for (len_t i = 0; i < mpz_param->nvars - 1; i++) {
    crt_lift_mpz_upoly(mpz_param->coords[i], nmod_param->coords[i], modulus,
                       prime, prod_crt, tmp, nthrds);
  }
}

/**

   la sortie est recons / denominator

 **/
#define RR 1

static inline int rational_reconstruction_mpz_ptr(
    mpz_t *recons, mpz_t denominator, mpz_t *pol, deg_t len, mpz_t modulus,
    deg_t *maxrec, mpq_t *coef, mpz_t rnum, mpz_t rden, mpz_t *tmp_num,
    mpz_t *tmp_den, mpz_t lcm, mpz_t guessed_num, mpz_t guessed_den,
    rrec_data_t rdata, int info_level) {

  if (ratrecon(rnum, rden, pol[*maxrec], modulus, rdata) == 0) {
    return 0;
  }
  mpz_set(tmp_num[*maxrec], rnum);
  mpz_set(tmp_den[*maxrec], rden);

  for (deg_t i = *maxrec + 1; i < len; i++) {
    int b = ratrecon(rnum, rden, pol[i], modulus, rdata);
    if (b == 0) {
      *maxrec = i - 1;
      return b;
    }

    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }
  for (deg_t i = 0; i < *maxrec; i++) {
    int b = ratrecon(rnum, rden, pol[i], modulus, rdata);

    if (b == 0) {
      if (info_level) {
        fprintf(stderr, "[*]");
      }
      *maxrec = MAX(i - 1, 0);
      return b;
    }

    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }

  mpz_set(lcm, tmp_den[0]);
  for (deg_t i = 1; i < len; i++) {
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }
  for (deg_t i = 0; i < len; i++) {
    mpz_divexact(tmp_den[i], lcm, tmp_den[i]);
  }
  for (deg_t i = 0; i < len; i++) {
    mpz_mul(tmp_num[i], tmp_num[i], tmp_den[i]);
  }
  for (deg_t i = 0; i < len; i++) {
    mpz_set(recons[i], tmp_num[i]);
  }
  mpz_set(denominator, lcm);

  return 1;
}

static inline int rational_reconstruction_mpz_ptr_with_denom(
    mpz_t *recons, mpz_t denominator, mpz_t *pol, deg_t len, mpz_t modulus,
    deg_t *maxrec, mpq_t *coef, mpz_t rnum, mpz_t rden, mpz_t *tmp_num,
    mpz_t *tmp_den, mpz_t lcm, mpz_t gnum, mpz_t guessed_den, rrec_data_t rdata,
    int info_level) {

  mpz_set(gnum, pol[*maxrec]);

  if (ratreconwden(rnum, rden, gnum, modulus, guessed_den, rdata) == 0) {
    return 0;
  }

  mpz_set(tmp_num[*maxrec], rnum);
  mpz_set(tmp_den[*maxrec], rden);

  for (deg_t i = *maxrec + 1; i < len; i++) {
    mpz_set(gnum, pol[i]);
    int b = ratreconwden(rnum, rden, gnum, modulus, guessed_den, rdata);

    if (b == 0) {
      *maxrec = MAX(0, i - 1);
      return b;
    }
    mpz_set(tmp_num[i], rnum);
    mpz_set(tmp_den[i], rden);
  }

  mpz_set(lcm, tmp_den[*maxrec]);
  for (deg_t i = *maxrec + 1; i < len; i++) {
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }

  mpz_t newlcm;
  mpz_init(newlcm);
  mpz_set(newlcm, lcm);
  mpz_mul(newlcm, newlcm, guessed_den);
  mpz_fdiv_q(rdata->D, rdata->D, lcm);
  mpz_mul(rdata->N, rdata->N, lcm);

  for (deg_t i = *maxrec - 1; i >= 0; i--) {
    mpz_set(gnum, pol[i]);
    int b = ratreconwden(tmp_num[i], tmp_den[i], gnum, modulus, newlcm, rdata);

    if (b == 0) {
      *maxrec = MAX(i + 1, 0);
      mpz_clear(newlcm);
      return b;
    }

    mpz_divexact(rden, newlcm, guessed_den);
    mpz_mul(tmp_den[i], tmp_den[i], rden);

    mpz_lcm(newlcm, newlcm, rden);
  }

  mpz_set(lcm, tmp_den[0]);
  for (deg_t i = 1; i < len; i++) {
    mpz_lcm(lcm, lcm, tmp_den[i]);
  }

  for (deg_t i = 0; i < len; i++) {
    mpz_divexact(tmp_den[i], lcm, tmp_den[i]);
  }
  for (deg_t i = 0; i < len; i++) {
    mpz_mul(tmp_num[i], tmp_num[i], tmp_den[i]);
  }
  for (deg_t i = 0; i < len; i++) {
    mpz_set(recons[i], tmp_num[i]);
  }
  mpz_set(denominator, lcm);
  mpz_clear(newlcm);
  return 1;
}

/**

   la sortie est recons / denominator

 **/

static inline int rational_reconstruction_upoly(
    mpz_upoly_t recons, mpz_t denominator, mpz_upoly_t pol, long len,
    mpz_t modulus, deg_t *maxrec, mpq_t *coef, mpz_t rnum, mpz_t rden,
    mpz_upoly_t tmp_num, mpz_upoly_t tmp_den, mpz_t lcm, mpz_t guessed_num,
    mpz_t guessed_den, rrec_data_t rdata, int info_level) {

  return rational_reconstruction_mpz_ptr(
      recons->coeffs, denominator, pol->coeffs, len, modulus, maxrec, coef,
      rnum, rden, tmp_num->coeffs, tmp_den->coeffs, lcm, guessed_num,
      guessed_den, rdata, info_level);
}

static inline int rational_reconstruction_upoly_with_denom(
    mpz_upoly_t recons, mpz_t denominator, mpz_upoly_t pol, long len,
    mpz_t modulus, deg_t *maxrec, mpq_t *coef, mpz_t rnum, mpz_t rden,
    mpz_upoly_t tmp_num, mpz_upoly_t tmp_den, mpz_t lcm, mpz_t guessed_num,
    mpz_t guessed_den, rrec_data_t rdata, int info_level) {

  return rational_reconstruction_mpz_ptr_with_denom(
      recons->coeffs, denominator, pol->coeffs, len, modulus, maxrec, coef,
      rnum, rden, tmp_num->coeffs, tmp_den->coeffs, lcm, guessed_num,
      guessed_den, rdata, info_level);
}

/**

returns 0 if rational reconstruction failed

 **/

static inline int new_rational_reconstruction(
    mpz_param_t mpz_param, mpz_param_t tmp_mpz_param, param_t *nmod_param,
    nvars_t nlins, nvars_t *linvars, uint32_t *lineqs, 
    trace_det_fglm_mat_t trace_det, sp_matfglm_t *mat, mpz_upoly_t numer,
    mpz_upoly_t denom, mpz_t modulus, mpz_t prod_crt, int32_t prime,
    mpq_t *coef, mpz_t rnum, mpz_t rden, rrec_data_t recdata,
    mpz_t *guessed_num, mpz_t *guessed_den, deg_t *maxrec, 
    deg_t *matrec, deg_t *oldmatrec_checked, deg_t *matrec_checked,
    int *is_lifted, int *mat_lifted, int *lin_lifted, int doit, 
    int nbdoit,
    int nthrds,
    const int info_level) {

  uint32_t trace_mod = nmod_param->elim->coeffs[trace_det->trace_idx];
  uint32_t det_mod = nmod_param->elim->coeffs[trace_det->det_idx];
  if(trace_det->lift_matrix && trace_det->mat_lifted < 2){
      copy_modular_matrix(trace_det, &mat, nthrds * nbdoit, prime);
  }
  
  int nr = check_trace_det_data(trace_det, maxrec, trace_mod, det_mod, 
          lineqs, mat, prime, trace_det->lift_matrix);
  if(nr > trace_det->w_checked){
      trace_det->w_checked = nr;
      if(info_level){
         fprintf(stderr, "[%.2f%%]", 100*(float)trace_det->w_checked/trace_det->nrows);
      }
  }
  if(trace_det->mat_lifted == 1){
      check_matrix(trace_det, mat, prime);
  }
  mpz_mul_ui(prod_crt, modulus, prime);
  if(trace_det->w_checked == trace_det->nrows && trace_det->mat_lifted == 0){
      mulmat_reconstruct(trace_det, prod_crt);
  }

  /**    CRT PART            **/
  crt_lift_mpz_param(tmp_mpz_param, nmod_param, modulus, prod_crt, prime,
                     trace_det->tmp, nthrds);

  crt_lift_trace_det(trace_det, trace_mod, det_mod, mat,
          lineqs, nmod_param,
          modulus, prod_crt, prime, nthrds);

  *matrec = *matrec_checked;

  mpz_mul_ui(modulus, modulus, prime);

  if (doit == 0) {
    return 0;
  }

  /**     CRT DONE                             **/

  /** RATIONAL RECONSTRUCTIONS                 **/
  mpz_sub_ui(*guessed_num, modulus, 1);
  mpz_fdiv_q_2exp(*guessed_num, *guessed_num, 1);
  mpz_sqrt(*guessed_num, *guessed_num);
  mpz_set(recdata->D, *guessed_num);
  mpz_set(recdata->N, *guessed_num);

  mpz_set_ui(rnum, 0);
  mpz_set_ui(rden, 1);
  if(trace_det->done_trace > 1){
      mpz_set(*guessed_den, trace_det->trace_den);
      if(trace_det->done_det > 1){
          mpz_lcm(*guessed_den, *guessed_den, trace_det->det_den);
      }
  }
  else{
      if(trace_det->done_det > 1){
          mpz_set(*guessed_den, trace_det->det_den);
          if(trace_det->done_trace > 1){
            mpz_lcm(*guessed_den, *guessed_den, trace_det->trace_den);
          }
       }
  }
  mpz_sub_ui(*guessed_num, modulus, 1);
  mpz_fdiv_q_2exp(*guessed_num, *guessed_num, 1);
  mpz_sqrt(*guessed_num, *guessed_num);
  mpz_set(recdata->N, *guessed_num);
  mpz_set(recdata->D, *guessed_num);

  if(trace_det->done_trace < 2 || trace_det->done_det < 2){
      rat_recon_trace_det(trace_det, recdata, modulus, rnum, rden, *guessed_den);
  }
  if (trace_det->done_trace > 1 && trace_det->done_det > 1) {

    mpz_sub_ui(*guessed_num, modulus, 1);
    mpz_fdiv_q_2exp(*guessed_num, *guessed_num, 1);
    mpz_sqrt(*guessed_num, *guessed_num);
    mpz_set(recdata->N, *guessed_num);
    mpz_set(recdata->D, *guessed_num);

    mpz_t denominator;
    mpz_init(denominator);
    mpz_t lcm;
    mpz_init(lcm);
    int b = 0;

    if (is_lifted[0] == 0) {

      mpz_root(recdata->D, modulus, 3);
      mpz_fdiv_q(recdata->N, modulus, recdata->D);
      mpz_fdiv_q_2exp(recdata->N, recdata->N, 1);
      b = rational_reconstruction_upoly_with_denom(
          mpz_param->elim, denominator, tmp_mpz_param->elim,
          nmod_param->elim->length, modulus, maxrec, coef, rnum, rden, numer,
          denom, lcm, *guessed_num, *guessed_den, recdata, info_level);
      if (b == 0) {
        mpz_root(recdata->D, modulus, 16);
        mpz_fdiv_q(recdata->N, modulus, recdata->D);
        mpz_fdiv_q_2exp(recdata->N, recdata->N, 1);

        b = rational_reconstruction_upoly_with_denom(
            mpz_param->elim, denominator, tmp_mpz_param->elim,
            nmod_param->elim->length, modulus, maxrec, coef, rnum, rden, numer,
            denom, lcm, *guessed_num, *guessed_den, recdata, info_level);
        if (b == 0) {
          is_lifted[0] = 0;
          mpz_clear(denominator);
          mpz_clear(lcm);
          return b;
        }
      }
      is_lifted[0] = 1;
      if (info_level) {
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

    mpz_fdiv_q_2exp(*guessed_num, modulus, 1);
    mpz_sqrt(recdata->D, *guessed_num);
    mpz_set(recdata->N, recdata->D);

    for (int i = 0; i < nc; i++) {
      *maxrec = MIN(MAX(0, trace_det->det_idx - 1),
                    MAX(0, nmod_param->coords[i]->length - 1));

      if (is_lifted[0] > 0 && is_lifted[i + 1] == 0) {

        b = rational_reconstruction_upoly_with_denom(
            mpz_param->coords[i], denominator, tmp_mpz_param->coords[i],
            nmod_param->coords[i]->length, modulus, maxrec, coef, rnum, rden,
            numer, denom, lcm, *guessed_num, *guessed_den, recdata, info_level);
        if (b == 0) {
          mpz_set_ui(recdata->D, 1);
          mpz_mul_2exp(recdata->D, recdata->D, nc);
          mpz_fdiv_q_2exp(recdata->N, modulus, 1);
          mpz_fdiv_q(recdata->N, recdata->N, recdata->D);

          b = rational_reconstruction_upoly_with_denom(
              mpz_param->coords[i], denominator, tmp_mpz_param->coords[i],
              nmod_param->coords[i]->length, modulus, maxrec, coef, rnum, rden,
              numer, denom, lcm, *guessed_num, *guessed_den, recdata,
              info_level);

          if (b == 0) {
            mpz_fdiv_q_2exp(recdata->N, modulus, 1);
            mpz_root(recdata->D, recdata->N, 16);
            mpz_fdiv_q(recdata->N, recdata->N, recdata->D);

            b = rational_reconstruction_upoly_with_denom(
                mpz_param->coords[i], denominator, tmp_mpz_param->coords[i],
                nmod_param->coords[i]->length, modulus, maxrec, coef, rnum,
                rden, numer, denom, lcm, *guessed_num, *guessed_den, recdata,
                info_level);
            if (b == 0) {

              mpz_clear(denominator);
              mpz_clear(lcm);
              mpz_clear(lc);
              mpq_clear(c);

              is_lifted[i + 1] = 0;
              return b;
            }
          }
        }
      } else {
        /* indicates that there is no need to set up below the data as they were
         * already computed */
        /* also denominator = 0 by now */
        b = 0;
      }
      if (info_level && b && is_lifted[i + 1] == 0) {
        fprintf(stderr, "[%d]", i + 1);
      }
      if (b) {
        is_lifted[i + 1] = 1;

        mpz_set(mpq_denref(c), denominator);
        mpq_canonicalize(c);
        mpz_set(mpz_param->cfs[i], mpq_denref(c));

        for (long j = 0; j < mpz_param->coords[i]->length; j++) {
          mpz_mul(mpz_param->coords[i]->coeffs[j],
                  mpz_param->coords[i]->coeffs[j], mpq_numref(c));
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
static inline int check_unit_mpz_nmod_poly(const deg_t len,
                                           const mpz_upoly_t mpz_pol,
                                           const nmod_poly_t nm_pol,
                                           const int32_t prime) {
  uint32_t lc = mpz_fdiv_ui(mpz_pol->coeffs[len - 1], prime);
  lc = mod_p_inverse_32(lc, prime);
  for (deg_t i = 0; i < len; i++) {
    uint64_t c = mpz_fdiv_ui(mpz_pol->coeffs[i], prime);
    c *= (uint64_t)lc;
    c = c % prime;
    if (c != nm_pol->coeffs[i]) {
      return i + 1;
    }
  }
  return 0;
}

static inline int
check_param_nmod_poly(const long len, const mpz_upoly_t mpz_pol,
                      const mpz_t den, const mpz_t lcelim, const long nbsol,
                      const nmod_poly_t nm_pol, const int32_t prime) {
  if (len == 0) {
    return 0;
  }
  mpz_t binv, bprime;
  mpz_init(binv);
  mpz_init_set_ui(bprime, prime);

  mpz_mul(binv, lcelim, den);
  mpz_mul_ui(binv, binv, nbsol);
  mpz_invert(binv, binv, bprime);

  uint32_t inv = mpz_mod_ui(binv, binv, prime);

  for (long i = 0; i < len; i++) {

    mpz_mul_ui(binv, mpz_pol->coeffs[i], inv);
    if (mpz_mod_ui(binv, binv, prime) != (nm_pol->coeffs[i] % prime)) {
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
                                      const int32_t prime, int *is_lifted,
                                      trace_det_fglm_mat_t trace_det,
                                      const int info_level) {

  long len = mp_param->nsols + 1;

  int c = check_unit_mpz_nmod_poly(len, mp_param->elim, bparam->elim, prime);
  if (c) {
    if (info_level) {
      fprintf(stderr, "<0,%d>", c);
    }
    is_lifted[0] = 0;
    for (int i = 0; i < mp_param->nvars - 1; i++) {
      is_lifted[i + 1] = 0;
    }
    trace_det->done_trace = 0;
    trace_det->check_trace = 0;
    trace_det->done_det = 0;
    trace_det->check_det = 0;
    return 1;
  }

  for (int i = 0; i < mp_param->nvars - 1; i++) {
    len = mp_param->coords[0]->length;

    if (check_param_nmod_poly(
            bparam->coords[i]->length, mp_param->coords[i], mp_param->cfs[i],
            mp_param->elim->coeffs[mp_param->elim->length - 1],
            mp_param->elim->length - 1, bparam->coords[i], prime)) {
      is_lifted[i + 1] = 0;
      if (info_level) {
        fprintf(stderr, "<%d>", i + 1);
      }
      return 1;
    }
  }
  return 0;
}

static inline void get_leading_ideal_information(int32_t *num_gb,
                                                 int32_t **lead_mons,
                                                 const int32_t pos,
                                                 const bs_t *const bs) {
  lead_mons[pos] = get_lm_from_bs(bs, bs->ht);
  num_gb[pos] = bs->lml;
}

static inline void print_groebner_basis(files_gb *files,
                                        const data_gens_ff_t *gens,
                                        const bs_t *const bs, md_t *md,
                                        const int32_t fc) {
  if (md->print_gb) {
    int32_t gfc = md->gfc;
    md->gfc = fc;
    print_ff_basis_data(files->out_file, "ab", bs, bs->ht, md, gens,
                        md->print_gb);
    md->gfc = gfc;
  }
}


static int32_t check_for_single_element_groebner_basis(
    int *dim, long *dquot_ori, const bs_t *const bs, int32_t **leadmons,
    const int32_t pos, const md_t *const md) {
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
	long **bextra_nf,
	int32_t **blens_extra_nf,
	int32_t **bexps_extra_nf,
	int32_t **bcfs_extra_nf,

                     nvars_t *nlins_ptr, nvars_t *linvars,
                     uint32_t **lineqs_ptr, nvars_t *squvars,

                     fglm_data_t **bdata_fglm, fglm_bms_data_t **bdata_bms,
        int32_t *num_gb,
        int32_t **leadmons,
        uint64_t *bsz,
        param_t **bparam,
        bs_t *gbg,
        md_t *md,
        const int32_t fc,
	const int32_t unstable_staircase,
        int print_gb,
        int *dim,
        long *dquot_ori,
        data_gens_ff_t *gens,
        files_gb *files,
        int *success)
{
    double rt = realtime();

  md->print_gb = print_gb;
  md->f4_qq_round = 1;

  int32_t error = 0;
  int32_t empty_solution_set = 1;
  bs_t *bs = core_gba(gbg, md, &error, fc);

  md->learning_rtime = realtime()-rt;
  print_tracer_statistics(stdout, rt, md);

  get_leading_ideal_information(num_gb, leadmons, 0, bs);

  print_groebner_basis(files, gens, bs, md, fc);

  empty_solution_set = check_for_single_element_groebner_basis(
      dim, dquot_ori, bs, leadmons, 0, md);

  if (empty_solution_set == 1) {
    return NULL;
  }


    check_and_set_linear_poly(nlins_ptr, linvars, lineqs_ptr, bs->ht,
                              leadmons[0], bs);
    if (has_dimension_zero(bs->lml, bs->ht->nv, leadmons[0])) {
        long dquot = 0;
        int32_t *lmb = monomial_basis(bs->lml, bs->ht->nv, leadmons[0], &dquot);

        /* if(md->info_level){ */
        /*     fprintf(stderr, "Dimension of quotient: %ld\n", dquot); */
        /* } */
        if(print_gb==0){
	    /* *bmatrix = build_matrixn_from_bs_trace(bdiv_xn, */
	    /* 					 blen_gb_xn, */
	    /* 					 bstart_cf_gb_xn, */
	    /* 					 lmb, dquot, bs, bs->ht, */
	    /* 					 leadmons[0], bs->ht->nv, */
	    /* 					 fc, */
	    /* 					 md->info_level); */
	    md->fglm_rtime = realtime();
	    md->fglm_ctime = cputime();
	    print_fglm_header (stdout,md);
	    *bmatrix = build_matrixn_unstable_from_bs_trace(bdiv_xn,
							    blen_gb_xn,
							    bstart_cf_gb_xn,
							    bextra_nf,
							    blens_extra_nf,
							    bexps_extra_nf,
							    bcfs_extra_nf,
							    lmb, dquot, bs, bs->ht,
							    leadmons[0], md,
							    bs->ht->nv,
							    fc, unstable_staircase,
							    md->info_level);

	    if(*bmatrix == NULL){
	      *success = 0;
	      *dim = 0;
	      *dquot_ori = dquot;
	      if(md->info_level > 1){
		fprintf (stdout,"------------------------------------------------------------------------------------------------------\n");
	      }
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




static void secondary_modular_steps(sp_matfglm_t **bmatrix,
				    int32_t **bdiv_xn,
				    int32_t **blen_gb_xn,
				    int32_t **bstart_cf_gb_xn,
				    long **bextra_nf,
				    int32_t **blens_extra_nf,
				    int32_t **bexps_extra_nf,
				    int32_t **bcfs_extra_nf,

				    nvars_t *bnlins,
				    nvars_t **blinvars,
				    uint32_t **blineqs,
				    nvars_t **bsquvars,

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
				    const int32_t unstable_staircase,
				    int info_level,
				    bs_t **bs,
				    int32_t *lmb_ori,
				    int32_t dquot_ori,
				    primes_t *lp,
				    data_gens_ff_t *gens,
				    double *stf4,
				    const long nbsols,
				    uint32_t *bad_primes,
                    trace_det_fglm_mat_t trace_det) {
    st->info_level  = 0;
    st->f4_qq_round = 2;

    double rt = realtime();
    /* tracing phase */
    len_t i;
    int32_t error = 0;

    /* F4 and FGLM are run using a single thread */
    /* st->nthrds is reset to its original value afterwards */
    const int nthrds = st->nthrds;
    st->nthrds = 1;
    for(nvars_t i = 0; i < st->nprimes; i++){
      bad_primes[i] = 0;
    }
#pragma omp parallel for num_threads(nthrds)  \
    private(i) schedule(static)
    for (i = 0; i < st->nprimes; ++i){
      if (trace_det->mat_lifted < 2 || trace_det->lin_lifted < 2) {
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
      }
      if (trace_det->mat_lifted == 2 && trace_det->lin_lifted == 2) {
        compute_modular_linear_forms(bnlins[i], bs_qq->ht->nv + 1, blineqs[i],
                                   trace_det->mpz_linear_forms, lp->p[i]);
        compute_modular_matrix(bmatrix[i], trace_det, lp->p[i]);
    } else {
      if (equal_staircase(leadmons_current[i], leadmons_ori[i], num_gb[i],
                          num_gb[i], bs[i]->ht->nv)) {

          set_linear_poly(bnlins[i], blineqs[i], blinvars[i], bs[i]->ht,
                    leadmons_current[i], bs[i]);
	    /* build_matrixn_from_bs_trace_application(bmatrix[i], */
	    /* 					    div_xn[i], */
	    /* 					    len_gb_xn[i], */
	    /* 					    start_cf_gb_xn[i], */
	    /* 					    lmb_ori, dquot_ori, bs[i], bs[i]->ht, */
	    /* 					    leadmons_ori[i], bs[i]->ht->nv, */
	    /* 					    lp->p[i]); */
	       build_matrixn_unstable_from_bs_trace_application(bmatrix[i],
							     bdiv_xn[i],
							     blen_gb_xn[i],
							     bstart_cf_gb_xn[i],
							     bextra_nf[i],
							     blens_extra_nf[i],
							     bexps_extra_nf[i],
							     bcfs_extra_nf[i],
							     lmb_ori, dquot_ori, bs[i], bs[i]->ht,
							     leadmons_ori[i], st, bs[i]->ht->nv,
							     lp->p[i],i);
        }
        else{
            bad_primes[i] = 1;
            }
     }
        if(nmod_fglm_compute_apply_trace_data(bmatrix[i], lp->p[i],
                                              nmod_params[i],
                                              bs_qq->ht->nv,
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
        if (bs[i] != NULL) {
            free_basis_and_only_local_hash_table_data(&(bs[i]));
        }
    }
  st->nthrds = nthrds;
}





/* sets function pointer */
void set_linear_function_pointer(int32_t fc) {
  int nbits = 0;
  if (fc == 0) {
    nbits = 32;
  }
  else{
    if(fc < (int32_t)(1) << 8) {
      nbits = 8;
    }
    else{
      if(fc < (int32_t)(1) << 16){
        nbits = 16;
      } else
        nbits = 32;
    }
  }
  switch (nbits) {
  case 8:
    set_linear_poly = set_linear_poly_8;
    check_and_set_linear_poly = check_and_set_linear_poly_8;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_8;
    copy_nf_in_matrix_from_bs = copy_nf_in_matrix_from_bs_8;
    break;
  case 16:
    set_linear_poly = set_linear_poly_16;
    check_and_set_linear_poly = check_and_set_linear_poly_16;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_16;
    copy_nf_in_matrix_from_bs = copy_nf_in_matrix_from_bs_16;
    break;
  case 32:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
    copy_nf_in_matrix_from_bs = copy_nf_in_matrix_from_bs_32;
    break;
  case 0:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
    copy_nf_in_matrix_from_bs = copy_nf_in_matrix_from_bs_32;
    break;
  default:
    set_linear_poly = set_linear_poly_32;
    check_and_set_linear_poly = check_and_set_linear_poly_32;
    copy_poly_in_matrix_from_bs = copy_poly_in_matrix_from_bs_32;
    copy_nf_in_matrix_from_bs = copy_nf_in_matrix_from_bs_32;
  }
}

static inline int is_lucky_matmul_prime_ui(uint32_t prime,
        trace_det_fglm_mat_t trace_det) {
  if (trace_det->mat_lifted < 2) {
    return 0;
  }
  for (uint32_t i = 0; i < trace_det->nrows; i++) {
    if (mpz_fdiv_ui(trace_det->mat_denoms[i], prime) == 0) {
      return 1;
    }
  }
  return 0;
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

int msolve_trace_qq(mpz_param_t *mpz_paramp,
                    param_t **nmod_param,
                    int *dim_ptr,
                    long *dquot_ptr,
                    data_gens_ff_t *gens,
                    int32_t ht_size, //initial_hts,
                    int32_t unstable_staircase,
                    int32_t nr_threads,
                    int32_t max_nr_pairs,
                    int32_t elim_block_len,
                    int32_t reset_ht,
                    int32_t la_option,
                    int32_t use_signatures,
                    int32_t lift_matrix,
                    int32_t info_level,
                    int32_t print_gb,
                    int32_t pbm_file,
                    files_gb *files,
                    int round){

  const int32_t *lens = gens->lens;
  const int32_t *exps = gens->exps;
  uint32_t field_char = gens->field_char;
  const void *cfs = gens->mpz_cfs;
  if (gens->field_char) {
    cfs = gens->cfs;
  } else {
    cfs = gens->mpz_cfs;
  }
  int mon_order = 0;
  int32_t nr_vars = gens->nvars;
  int32_t nr_gens = gens->ngens;
  int reduce_gb = 1;
  int32_t nr_nf = 0;
  const uint32_t prime_start = (uint32_t)(1) << 30;
  const int32_t nr_primes = nr_threads;

  len_t i;

  /* initialize stuff */
  md_t *st = allocate_meta_data();

  int truncate_lifting = 0;
  int *invalid_gens = NULL;
  int res = validate_input_data(
      &invalid_gens, cfs, lens, &field_char, &mon_order, &elim_block_len,
      &nr_vars, &nr_gens, &nr_nf, &ht_size, &nr_threads, &max_nr_pairs,
      &reset_ht, &la_option, &use_signatures, &reduce_gb, 
      &truncate_lifting, &info_level);

  /* all data is corrupt */
  if (res == -1) {
    fprintf(stderr, "Invalid input generators, msolve now terminates.\n");
    free(invalid_gens);
    return -3;
  }
  /* checks and set all meta data. if a nonzero value is returned then
   * some of the input data is corrupted. */

  if (check_and_set_meta_data_trace(
          st, lens, exps, cfs, invalid_gens, field_char, mon_order,
          elim_block_len, nr_vars, nr_gens, nr_nf, ht_size, nr_threads,
          max_nr_pairs, reset_ht, la_option, use_signatures, reduce_gb,
          prime_start, nr_primes, pbm_file, 0 /*truncate_lifting */, 
          info_level)) {
    free(st);
    return -3;
  }

  /* lucky primes */
  primes_t *lp = (primes_t *)calloc(st->nthrds, sizeof(primes_t));

  /*******************
   * initialize basis
   *******************/
  bs_t *bs_qq = initialize_basis(st);
  /* read in ideal, move coefficients to integers */
  import_input_data(bs_qq, st, 0, st->ngens_input, lens, exps, cfs,
                    invalid_gens);
  free(invalid_gens);
  invalid_gens = NULL;

  print_initial_statistics(stderr, st);

  /* for faster divisibility checks, needs to be done after we have
   * read some input data for applying heuristics */
  calculate_divmask(bs_qq->ht);

  /* sort initial elements, smallest lead term first */
  sort_r(bs_qq->hm, (unsigned long)bs_qq->ld, sizeof(hm_t *), initial_input_cmp,
         bs_qq->ht);
  if (gens->field_char == 0) {
    remove_content_of_initial_basis(bs_qq);
    /* generate lucky prime numbers */
    generate_lucky_primes(lp, bs_qq, st->prime_start, st->nthrds);
  } else {
    lp->old = 0;
    lp->ld = 1;
    lp->p = calloc(1, sizeof(uint32_t));
    normalize_initial_basis(bs_qq, st->fc);
  }

  /* generate array to store modular bases */
  bs_t **bs = (bs_t **)malloc((unsigned long)st->nthrds * sizeof(bs_t *));

  param_t **nmod_params =
      (param_t **)malloc((unsigned long)st->nthrds * sizeof(param_t *));

  uint32_t *bad_primes = calloc((unsigned long)st->nthrds, sizeof(uint32_t));

  uint32_t prime = 0;
  uint32_t primeinit = 0;
  uint32_t lprime = 1303905299;
  srand(time(0));
  prime = next_prime(rand() % (1303905301 - (1 << 30) + 1) + (1 << 30));
  while (gens->field_char == 0 && is_lucky_prime_ui(prime, bs_qq)) {
    prime = next_prime(rand() % (1303905301 - (1 << 30) + 1) + (1 << 30));
  }
  primeinit = prime;
  lp->p[0] = primeinit;

  if (gens->field_char) {
    lp->p[0] = gens->field_char;
  }
  sp_matfglm_t **bmatrix =
      (sp_matfglm_t **)malloc(st->nthrds * sizeof(sp_matfglm_t *));

  int32_t **bdiv_xn = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  int32_t **blen_gb_xn = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  int32_t **bstart_cf_gb_xn = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  long **bextra_nf = (long **)malloc(st->nthrds * sizeof(long *));
  int32_t **blens_extra_nf = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  int32_t **bexps_extra_nf = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  int32_t **bcfs_extra_nf = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));

  fglm_data_t **bdata_fglm =
      (fglm_data_t **)malloc(st->nthrds * sizeof(fglm_data_t *));
  fglm_bms_data_t **bdata_bms =
      (fglm_bms_data_t **)malloc(st->nthrds * sizeof(fglm_bms_data_t *));
  int32_t *num_gb = (int32_t *)calloc(st->nthrds, sizeof(int32_t));
  int32_t **leadmons_ori = (int32_t **)malloc(st->nthrds * sizeof(int32_t *));
  int32_t **leadmons_current =
      (int32_t **)malloc(st->nthrds * sizeof(int32_t *));

  uint64_t bsz = 0;

  int32_t nv = bs_qq->ht->nv;

  /* data for linear forms */
  nvars_t nlins = 0; /*number of linear forms*/
  nvars_t *bnlins = (nvars_t *)calloc(st->nthrds, sizeof(nvars_t));
  nvars_t **blinvars = (nvars_t **)malloc(
      st->nthrds * sizeof(nvars_t *)); /*indicates which variables are linear*/
  nvars_t *linvars = calloc(bs_qq->ht->nv, sizeof(nvars_t));
  blinvars[0] = linvars;
  uint32_t **lineqs_ptr =
      malloc(st->nthrds * sizeof(uint32_t *)); /*coeffs of linear forms*/

  /*data for squared variables*/
  nvars_t **bsquvars = (nvars_t **)malloc(st->nthrds * sizeof(nvars_t *));
  nvars_t *squvars = calloc(nr_vars - 1, sizeof(nvars_t));
  bsquvars[0] = squvars;


  set_linear_function_pointer(gens->field_char);

  int success = 1;
  int squares = 1;

  int32_t *lmb_ori = initial_modular_step(bmatrix, bdiv_xn, blen_gb_xn,
					  bstart_cf_gb_xn,
					  bextra_nf,
					  blens_extra_nf,
					  bexps_extra_nf,
					  bcfs_extra_nf,

					  &nlins, blinvars[0], lineqs_ptr,
					  squvars,

					  bdata_fglm, bdata_bms,

					  num_gb, leadmons_ori,

					  &bsz, nmod_params,
					  bs_qq, st,
					  lp->p[0], //prime,
					  unstable_staircase,
					  print_gb,
					  dim_ptr, dquot_ptr,
					  gens,
					  files,
					  &success);


  if (*dim_ptr == 0 && success && *dquot_ptr > 0 && print_gb == 0) {
    if (nmod_params[0]->elim->length - 1 != *dquot_ptr) {
      for (int i = 0; i < nr_vars - 1; i++) {
        if ((squvars[i] == 0) && round) {
          squares = 0;
          success = 0;
        }
      }
    }
  }

  (*mpz_paramp)->dim = *dim_ptr;
  (*mpz_paramp)->dquot = *dquot_ptr;

  if (lmb_ori == NULL || success == 0 || print_gb || gens->field_char) {
    free(bs);
    if (gens->field_char == 0) {
      free_basis(&bs_qq);
      /*nmod_params[0] should not be cleaned here (change of primitive
       * element_*/
    }
    free_lucky_primes(&lp);
    free(bad_primes);
    free(lp);
    free(linvars);
    if (nlins) {
      free(lineqs_ptr[0]);
    }
    free(bnlins);
    free(lineqs_ptr);
    free(squvars);
    free(lmb_ori);
    if (print_gb) {
      return 0;
    }
    if (*dim_ptr == 0 && gens->field_char && success) {
      /* copy of parametrization */
      if (*dquot_ptr != 0) {
        param_t *par = allocate_fglm_param(gens->field_char, st->nvars);
        nmod_poly_set(par->elim, nmod_params[0]->elim);
        nmod_poly_set(par->denom, nmod_params[0]->denom);
        for (long j = 0; j <= st->nvars - 2; j++) {
          nmod_poly_set(par->coords[j], nmod_params[0]->coords[j]);
        }
        free_fglm_param(nmod_params[0]);
        (*nmod_param) = par;
      }
      free(st);
      return 0;
    }
    free(st);
    if (*dim_ptr == 1) {
      if (info_level) {
        fprintf(stderr, "Positive dimensional Grobner basis\n");
      }
      return 0;
    }
    if (*dquot_ptr == 0) {
      return 0;
    }
    if (*dquot_ptr > 0) {
      if (squares == 0) {
        return 2;
      }
      return 1;
    }
  }

  /* duplicate data for multi-threaded multi-mod computation */
  duplicate_data_mthread_trace(st->nthrds, bs_qq, st, num_gb,
                              leadmons_ori, leadmons_current,
                               bdata_bms, bdata_fglm,
                               bstart_cf_gb_xn, blen_gb_xn, bdiv_xn,
			       bextra_nf,blens_extra_nf,bexps_extra_nf,bcfs_extra_nf,
			       bmatrix,
                               nmod_params, nlins, bnlins,
                               blinvars, lineqs_ptr,
                               bsquvars);
  normalize_nmod_param(nmod_params[0]);

  /* if (info_level) { */
  /*   fprintf(stderr, "\nStarts multi-modular computations\n"); */
  /* } */
  /* print postponed */

  mpz_param_t tmp_mpz_param;
  mpz_param_init(tmp_mpz_param);

  initialize_mpz_param(*mpz_paramp, nmod_params[0]);
  initialize_mpz_param(tmp_mpz_param, nmod_params[0]);

  // attention les longueurs des mpz_param sont fixees par nmod_params[0]
  // dans des cas exceptionnels, ca peut augmenter avec un autre premier.

  /* data for rational reconstruction of trace and det of mult. mat. */
  trace_det_fglm_mat_t trace_det;
  uint32_t detidx = 0;
  /* int32_t tridx = nmod_params[0]->elim->length-2; */
  int32_t tridx = 3 * (nmod_params[0]->elim->length - 1) / 4;
  /* tridx = nmod_params[0]->elim->length - 2; */
  while (nmod_params[0]->elim->coeffs[tridx] == 0 && tridx > 0) {
    tridx--;
  }
  detidx = 2 * (nmod_params[0]->elim->length - 1) / 3;
  while (nmod_params[0]->elim->coeffs[detidx] == 0 &&
         detidx < nmod_params[0]->elim->length - 2) {
    detidx++;
  }
  trace_det_initset(trace_det, nmod_params[0]->elim->coeffs[tridx],
                    nmod_params[0]->elim->coeffs[detidx], tridx, detidx, 
                    bmatrix, lift_matrix, nlins, nv, lineqs_ptr, 4*st->nthrds, 
                    lp->p[0]);
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
  mpz_set_ui(rnum, 0);
  mpz_init(rden);
  mpz_set_ui(rden, 1);
  set_mpz_param_nmod(tmp_mpz_param, nmod_params[0]);

  deg_t nsols = tmp_mpz_param->nsols;

  mpz_upoly_t numer;
  mpz_upoly_init2(numer, (nsols + 1), 32 * (nsols + 1));
  numer->length = nsols + 1;

  mpz_upoly_t denom;
  mpz_upoly_init2(denom, (nsols + 1), 64 * (nsols + 1));
  denom->length = nsols + 1;

  mpz_t guessed_den;
  mpz_init2(guessed_den, 32 * nsols);
  mpz_set_ui(guessed_den, 1);
  mpz_t guessed_num;
  mpz_init2(guessed_num, 32 * nsols);
  mpz_set_ui(guessed_num, 0);

  deg_t maxrec = 0;
  deg_t matrec = 0;
  deg_t oldmatrec_checked = 0;
  deg_t matrec_checked = 0;

  int rerun = 1, nprimes = 1, mcheck = 1;

  long nbadprimes = 0;

  int *is_lifted = calloc(nr_vars, sizeof(int));
  int mat_lifted = 0;
  int lin_lifted = 0;
  if(nlins == 0){
      lin_lifted = 2;
  }
  int nbdoit = 1;
  int doit = 1;
  int prdone = 0;
  int lpow2 = 1;
  int clog = 0;
  int br = 0;

  rrec_data_t recdata;
  initialize_rrec_data(recdata);

  /* measures time spent in rational reconstruction */
  double strat = 0;

  while (rerun == 1 || mcheck == 1) {
    /* controls call to rational reconstruction */
    doit = ((prdone % nbdoit) == 0);

    /* generate lucky prime numbers */
    prime = next_prime(prime);
    if (prime >= lprime) {
      prime = next_prime(1 << 30);
    }
    lp->p[0] = prime;
    while (is_lucky_prime_ui(prime, bs_qq) || prime == primeinit) {
      prime = next_prime(prime);
      if (prime >= lprime) {
        prime = next_prime(1 << 30);
      }
      lp->p[0] = prime;
    }

    for (len_t i = 1; i < st->nthrds; i++) {
      prime = next_prime(prime);
      if (prime >= lprime) {
        prime = next_prime(1 << 30);
      }
      lp->p[i] = prime;
      if(trace_det->lift_matrix){
          while (is_lucky_prime_ui(prime, bs_qq) || prime == primeinit ||
             is_lucky_matmul_prime_ui(prime, trace_det)) {
              prime = next_prime(prime);
              if (prime >= lprime) {
                  prime = next_prime(1 << 30);
              }
              lp->p[i] = prime;
          }
      }
      else{
        while (is_lucky_prime_ui(prime, bs_qq) || prime == primeinit) {
          prime = next_prime(prime);
          if (prime >= lprime) {
            prime = next_prime(1 << 30);
          }
          lp->p[i] = prime;
        }
      }
    }
    prime = lp->p[st->nthrds - 1];

    double ca0 = realtime();
    double stf4 = 0;
    secondary_modular_steps(bmatrix,
			    bdiv_xn,
			    blen_gb_xn,
			    bstart_cf_gb_xn,
			    bextra_nf,
			    blens_extra_nf,
			    bexps_extra_nf,
			    bcfs_extra_nf,

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
			    field_char, unstable_staircase, 0, /* info_level, */
			    bs, lmb_ori, *dquot_ptr, lp,
			    gens, &stf4, nsols, bad_primes,
                trace_det);
    double ca1 = realtime() - ca0;

    if (nprimes == 1) {
      if (info_level > 2) {
        fprintf(stderr, "------------------------------------------\n");
        fprintf(stderr, "#ADDITIONS       %13lu\n",
                (unsigned long)st->application_nr_add * 1000);
        fprintf(stderr, "#MULTIPLICATIONS %13lu\n",
                (unsigned long)st->application_nr_mult * 1000);
        fprintf(stderr, "#REDUCTIONS      %13lu\n",
                (unsigned long)st->application_nr_red);
        fprintf(stderr, "------------------------------------------\n");
      }
      /* if (info_level > 1) { */
      /*   fprintf(stderr, "Application phase %.2f Gops/sec\n", */
      /*           (st->application_nr_add + st->application_nr_mult) / 1000.0 / */
      /*               1000.0 / (stf4)); */
      /*   fprintf(stderr, "Multi-mod time: GB + fglm (elapsed): %.2f sec\n", */
      /*           (ca1) ); */
      /* } */
      if(info_level){
	    fprintf(stdout,
		    "\n---------------- TIMINGS ----------------\n");
	    fprintf(stdout,
		    "multi-mod overall(elapsed) %9.2f sec\n",
		    ca1);
	    fprintf(stdout,
		    "multi-mod F4               %9.2f sec\n",
		    stf4);
	    fprintf(stdout,
		    "multi-mod FGLM             %9.2f sec\n",
		    ca1-stf4);
	    if (info_level > 1){
	      fprintf(stdout,
		      "learning phase             %9.2f Gops/sec\n",
		      (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/(st->learning_rtime));
	      fprintf(stdout,
		      "application phase          %9.2f Gops/sec\n",
		      (st->application_nr_add+st->application_nr_mult)/1000.0/1000.0/(stf4));
	    }
	    fprintf(stdout,
		    "-----------------------------------------\n");
      }
      if (info_level) {
	  fprintf(stdout,
		  "\nmulti-modular steps\n");
	  fprintf(stdout, "-------------------------------------------------\
-----------------------------------------------------\n");
      }

    }
    for (int i = 0; i < st->nthrds; i++) {
      if (bad_primes[i] == 0) {
        normalize_nmod_param(nmod_params[i]);
      }
    }
    /* scrr measures time spent in ratrecon for modular images */
    double crr = 0, scrr = 0;
    /* CRT + rational reconstruction */
    for (len_t i = 0; i < st->nthrds; i++) {
      if (bad_primes[i] == 0) {
        if (rerun == 0) {
          mcheck = check_param_modular(*mpz_paramp, nmod_params[i], lp->p[i],
                                       is_lifted, trace_det, info_level);
        }
        crr = realtime();
        if (mcheck == 1) {
          br = new_rational_reconstruction(
              *mpz_paramp, tmp_mpz_param, nmod_params[i], 
              bnlins[i], blinvars[i], lineqs_ptr[i], 
              trace_det, bmatrix[i], numer,
              denom, modulus, prod_crt, lp->p[i], &result, rnum, rden, recdata,
              &guessed_num, &guessed_den, &maxrec, &matrec, &oldmatrec_checked, 
              &matrec_checked, is_lifted,
              &mat_lifted, &lin_lifted, doit, nbdoit, st->nthrds, info_level);

          if (br == 1) {
            rerun = 0;
          } else {
            rerun = 1;
          }
        }
        scrr += realtime() - crr;
        nprimes++;
      } else {
        if (info_level) {
          fprintf(stdout, "<bp: %d>\n", lp->p[i]);
	  fflush(stdout);
        }
        nbadprimes++;
        if (nbadprimes > nprimes) {
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

    double t = ((double)nbdoit) * ca1;
    if ((t == 0) || (scrr >= 0.2 * t && br == 0)) {
      nbdoit = 2 * nbdoit;
      lpow2 = 2 * nprimes;
      doit = 0;
      if (info_level) {
        fprintf(stdout, "\n<Step:%d/%.2f/%.2f>", nbdoit, scrr, t);
	fflush(stdout);
      }
      prdone = 0;
    } else {
      prdone++;
    }

    if ((LOG2(nprimes) > clog) ||
        (nbdoit != 1 && (nprimes % (lpow2 + 1) == 0))) {
      if (info_level) {
        fprintf(stdout, "{%d}", nprimes);
	fflush(stdout);
      }
      clog++;
      lpow2 = 2 * lpow2;
    }
  }

  (*mpz_paramp)->denom->length = (*mpz_paramp)->nsols;
  for (long i = 1; i <= (*mpz_paramp)->nsols; i++) {
    mpz_set((*mpz_paramp)->denom->coeffs[i - 1],
            (*mpz_paramp)->elim->coeffs[i]);
    mpz_mul_ui((*mpz_paramp)->denom->coeffs[i - 1],
               (*mpz_paramp)->denom->coeffs[i - 1], i);
  }
  if(info_level){
    fprintf(stdout,
	    "\n-------------------------------------------------\
-----------------------------------------------------\n");
  }


  if(info_level){
    /* fprintf(stderr, "\n%d primes used\n", nprimes); */
    /* fprintf(stderr, "Time for CRT + rational reconstruction = %.2f\n", strat); */
    fprintf(stdout,"\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(stdout, "#primes            %16lu\n", (unsigned long) nprimes);
    fprintf(stdout, "#bad primes        %16lu\n", (unsigned long) nbadprimes);
    fprintf(stdout, "-----------------------------------------\n");
    fprintf(stdout, "\n---------------- TIMINGS ----------------\n");
    fprintf(stdout, "CRT and ratrecon(elapsed) %10.2f sec\n", st->fglm_rtime);
    fprintf(stdout, "-----------------------------------------\n");
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

  for (i = 0; i < st->nthrds; ++i) {
    if (bs[i] != NULL) {
      free_basis(&(bs[i]));
    }
    free_fglm_bms_data(bdata_bms[i]);
    free_fglm_data(bdata_fglm[i]);
    
    free_fglm_param(nmod_params[i]);
    free(bcfs_extra_nf[i]);
    free(bexps_extra_nf[i]);
    free(blens_extra_nf[i]);
    free(bextra_nf[i]);
    free(blen_gb_xn[i]);
    free(bstart_cf_gb_xn[i]);
    free(bdiv_xn[i]);
    posix_memalign_free(bmatrix[i]->dense_mat);
    posix_memalign_free(bmatrix[i]->dense_idx);
    posix_memalign_free(bmatrix[i]->triv_idx);
    posix_memalign_free(bmatrix[i]->triv_pos);
    posix_memalign_free(bmatrix[i]->dst);
    free(bmatrix[i]);
    free(leadmons_ori[i]);
    free(leadmons_current[i]);
    /* free_trace(&btrace[i]); */

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
  free(lmb_ori);

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
  free(bcfs_extra_nf);
  free(bexps_extra_nf);
  free(blens_extra_nf);
  free(bextra_nf);
  free(blen_gb_xn);
  free(bstart_cf_gb_xn);
  free(bdiv_xn);
  /* free(btrace); */

  return 0;
}

void real_point_init(real_point_t pt, long nvars) {
  pt->nvars = nvars;
  pt->coords = malloc(sizeof(coord_t) * nvars);
  for (long i = 0; i < nvars; i++) {
    mpz_init(pt->coords[i]->val_up);
    mpz_init(pt->coords[i]->val_do);
    pt->coords[i]->k_up = 0;
    pt->coords[i]->k_do = 0;
    pt->coords[i]->isexact = 0;
  }
}

void real_point_clear(real_point_t pt) {
  for (long i = 0; i < pt->nvars; i++) {
    mpz_clear(pt->coords[i]->val_up);
    mpz_clear(pt->coords[i]->val_do);
  }
  free(pt->coords);
}

void display_real_point(FILE *fstream, real_point_t pt) {
  fprintf(fstream, "[");
  for (long i = 0; i < pt->nvars - 1; i++) {
    fprintf(fstream, "[");
    mpz_out_str(fstream, 10, pt->coords[i]->val_do);
    if (pt->coords[i]->k_do && mpz_sgn(pt->coords[i]->val_do)) {
      fprintf(fstream, " / ");
      fprintf(fstream, "2");
      if (pt->coords[i]->k_do > 1) {
        fprintf(fstream, "^%d", pt->coords[i]->k_do);
      }
    }
    fprintf(fstream, ", ");
    mpz_out_str(fstream, 10, pt->coords[i]->val_up);
    if (pt->coords[i]->k_up && mpz_sgn(pt->coords[i]->val_up)) {
      fprintf(fstream, " / ");
      fprintf(fstream, "2");
      if (pt->coords[i]->k_up > 1) {
        fprintf(fstream, "^%d", pt->coords[i]->k_up);
      }
    }
    fprintf(fstream, "], ");
  }
  fprintf(fstream, "[");
  mpz_out_str(fstream, 10, pt->coords[pt->nvars - 1]->val_do);
  if (pt->coords[pt->nvars - 1]->k_do &&
      mpz_sgn(pt->coords[pt->nvars - 1]->val_do)) {
    fprintf(fstream, " / ");
    fprintf(fstream, "2");
    if (pt->coords[pt->nvars - 1]->k_do > 1) {
      fprintf(fstream, "^%d", pt->coords[pt->nvars - 1]->k_do);
    }
  }
  fprintf(fstream, ", ");

  mpz_out_str(fstream, 10, pt->coords[pt->nvars - 1]->val_up);
  if (pt->coords[pt->nvars - 1]->k_up &&
      mpz_sgn(pt->coords[pt->nvars - 1]->val_up)) {
    fprintf(fstream, " / ");
    fprintf(fstream, "2");
    if (pt->coords[pt->nvars - 1]->k_up > 1) {
      fprintf(fstream, "^%d", pt->coords[pt->nvars - 1]->k_up);
    }
  }

  fprintf(fstream, "]");
  fprintf(fstream, "]");
}

void display_real_points(FILE *fstream, real_point_t *pts, long nb) {
  fprintf(fstream, "[1,\n"); /* because at the moment we return a single list */
  fprintf(fstream, "[");
  for (long i = 0; i < nb - 1; i++) {
    display_real_point(fstream, pts[i]);
    fprintf(fstream, ", ");
  }
  if (nb) {
    display_real_point(fstream, pts[nb - 1]);
  }
  fprintf(fstream, "]\n");
  fprintf(fstream, "]");
}

void single_exact_real_root_param(mpz_param_t param, interval *rt, long nb,
                                  mpz_t *xdo, mpz_t *xup, mpz_t den_up,
                                  mpz_t den_do, mpz_t c, mpz_t tmp,
                                  mpz_t val_do, mpz_t val_up, mpz_t *tab,
                                  real_point_t pt, long prec, int info_level) {

  mpz_poly_eval_2exp_naive(param->denom->coeffs, param->denom->length - 1,
                           &rt->numer, rt->k, tab, tab + 1);

  mpz_set(den_up, tab[0]);
  mpz_set(den_do, tab[0]);

  for (long nv = 0; nv < param->nvars - 1; nv++) {

    mpz_poly_eval_2exp_naive(param->coords[nv]->coeffs,
                             param->coords[nv]->length - 1, &rt->numer, rt->k,
                             tab, tab + 1);
    mpz_set(val_up, tab[0]);
    mpz_set(val_do, tab[0]);

    mpz_neg(val_do, val_do);
    mpz_neg(val_up, val_up);
    mpz_swap(val_up, val_do);

    long exp = (rt->k) * ((param->denom->length) - (param->coords[nv]->length));

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

void generate_table_values(interval *rt, mpz_t c, const long ns, const long b,
                           const long corr, mpz_t *xdo, mpz_t *xup) {

  mpz_add_ui(c, rt->numer, 1);

  if (mpz_sgn(rt->numer) >= 0) {

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for (long i = 1; i < ns; i++) {
      if (i <= b) {
        mpz_mul(xup[i], xup[i - 1], c);
        mpz_mul(xdo[i], xdo[i - 1], rt->numer);
      }
      if ((i % b) == 0 && i > b) {
        long q = i / b;
        mpz_mul(xup[i], xup[(q - 1) * b], xup[b]);
        mpz_mul(xdo[i], xdo[(q - 1) * b], xdo[b]);
      }
    }
  } else {

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for (long i = 1; i < ns; i++) {
      if (i <= b) {
        if ((i & 1) == 1) {
          mpz_mul(xup[i], xdo[i - 1], c);
          mpz_mul(xdo[i], xup[i - 1], rt->numer);
        } else {
          mpz_mul(xup[i], xdo[i - 1], rt->numer);
          mpz_mul(xdo[i], xup[i - 1], c);
        }
      }
      if ((i % b) == 0 && i > b) {
        long q = i / b;
        mpz_mul(xup[i], xdo[(q - 1) * b], xup[b]);
        mpz_mul(xdo[i], xup[(q - 1) * b], xdo[b]);
      }
    }
  }

  long q = (ns - 1) / b;

  for (long i = 1; i <= q; i++) {
    mpz_mul_2exp(xup[i * b], xup[i * b], corr);
    mpz_cdiv_q_2exp(xup[i * b], xup[i * b], (rt->k) * i * b);

    mpz_mul_2exp(xdo[i * b], xdo[i * b], corr);
    mpz_fdiv_q_2exp(xdo[i * b], xdo[i * b], (rt->k) * i * b);
  }
}

/*
 xdo[i] and xup[i] are rt^i and c^i at precision 2^corr
[xdo[i]/2^corr, xup[i]/2^corr] contain [rt^i, c^i]
 */
void generate_table_values_full(interval *rt, mpz_t c, const long ns,
                                const long b, const long corr, mpz_t *xdo,
                                mpz_t *xup) {

  mpz_add_ui(c, rt->numer, 1);

  if (mpz_sgn(rt->numer) >= 0) {

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for (long i = 1; i < ns; i++) {

      mpz_mul(xup[i], xup[i - 1], c);
      mpz_mul(xdo[i], xdo[i - 1], rt->numer);
    }
  } else {

    mpz_set_ui(xup[0], 1);
    mpz_set_ui(xdo[0], 1);
    for (long i = 1; i < ns; i++) {
      if ((i & 1) == 1) {
        mpz_mul(xup[i], xdo[i - 1], c);
        mpz_mul(xdo[i], xup[i - 1], rt->numer);
      } else {
        mpz_mul(xup[i], xdo[i - 1], rt->numer);
        mpz_mul(xdo[i], xup[i - 1], c);
      }
    }
  }
  mpz_mul_2exp(xdo[0], xdo[0], corr);
  mpz_mul_2exp(xup[0], xup[0], corr);
  for (long i = 1; i < ns; i++) {
    mpz_mul_2exp(xup[i], xup[i], corr);
    mpz_cdiv_q_2exp(xup[i], xup[i], (rt->k) * i);

    mpz_mul_2exp(xdo[i], xdo[i], corr);
    mpz_fdiv_q_2exp(xdo[i], xdo[i], (rt->k) * i);
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
int evalquadric(mpz_t *quad, mpz_t r, long k, mpz_t *tmpquad, mpz_t tmp) {

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
  mpz_mul_2exp(tmp, tmp, 2 * k);
  /*tmp = c*2^(2k)*/
  mpz_add(tmpquad[0], tmpquad[0], tmp);

  /* At this stage, one has computed
     tmpquad = numer(subs(x=x*(rr2-rr1),subs(x=x+rr1, quad)))
     We want to know if it has real roots in [0, 1]
   */
  /* If all coefficients have the same sign, there is no root */
  int s = mpz_sgn(tmpquad[0]);
  if (s == mpz_sgn(tmpquad[1]) && s == mpz_sgn(tmpquad[2])) {
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
  s = mpz_sgn(tmpquad[0]);
  if (s == mpz_sgn(tmpquad[1]) && s == mpz_sgn(tmpquad[2])) {
    return 0;
  }
  return 1;
}

/* evaluates denom (which has degree deg)
   at the interval [r/2^k, (r+1)/2^k]
   returns  */
int value_denom(mpz_t *denom, long deg, mpz_t r, long k, mpz_t *xdo, mpz_t *xup,
                mpz_t tmp, mpz_t den_do, mpz_t den_up, long corr, mpz_t c) {

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

  int boo = mpz_poly_eval_interval(denom, deg, k, r, c, tmp, den_do, den_up);
  if (mpz_cmp(den_do, den_up) > 0) {
    fprintf(stderr, "BUG (den_do > den_up)\n");
    exit(1);
  }
  mpz_mul_2exp(den_do, den_do, corr);
  mpz_mul_2exp(den_up, den_up, corr);
  mpz_fdiv_q_2exp(den_do, den_do, k * deg);
  mpz_cdiv_q_2exp(den_up, den_up, k * deg);

  if (mpz_sgn(den_do) != mpz_sgn(den_up)) {
    return 1;
  }
  return boo;
}

int newvalue_denom(mpz_t *denom, long deg, mpz_t r, long k, mpz_t *xdo,
                   mpz_t *xup, mpz_t tmp, mpz_t den_do, mpz_t den_up, long corr,
                   mpz_t c) {

  mpz_add_ui(c, r, 1);
  /*boo = 1 if sgn(den_do) != sgn(den_up) else it is 0*/
  int boo = mpz_poly_eval_interval(denom, deg, k, r, c, tmp, den_do, den_up);
  if (mpz_cmp(den_do, den_up) > 0) {
    fprintf(stderr, "BUG (den_do > den_up)\n");
    exit(1);
  }
  mpz_mul_2exp(den_do, den_do, corr);
  mpz_mul_2exp(den_up, den_up, corr);
  mpz_fdiv_q_2exp(den_do, den_do, k * deg);
  mpz_cdiv_q_2exp(den_up, den_up, k * deg);

  return (boo || (mpz_sgn(den_do)==0) || (mpz_sgn(den_up)==0));
}

void lazy_single_real_root_param(mpz_param_t param, mpz_t *polelim,
                                 interval *rt, long nb, interval *pos_root,
                                 mpz_t *xdo, mpz_t *xup, mpz_t den_up,
                                 mpz_t den_do, mpz_t c, mpz_t tmp, mpz_t val_do,
                                 mpz_t val_up, mpz_t *tab, real_point_t pt,
                                 long prec, long nbits, mpz_t s,
                                 int info_level) {
  long ns = param->nsols;
  /* root is exact */
  if (rt->isexact == 1) {
    single_exact_real_root_param(param, rt, nb, xdo, xup, den_up, den_do, c,
                                 tmp, val_do, val_up, tab, pt, MAX(rt->k, prec),
                                 info_level);
    return;
  }

  int64_t b = 16;
  int64_t corr = 2 * (ns + rt->k);

  /* checks whether the abs. value of the root is greater than 1 */

  generate_table_values_full(rt, c, ns, b, corr, xdo, xup);
  while (rt->isexact == 0 && 
          newvalue_denom(param->denom->coeffs, param->denom->length - 1,
                        rt->numer, rt->k, xdo, xup, tmp, den_do, den_up, corr,
                        s)) {

    /* root is positive */
    if (mpz_sgn(rt->numer) >= 0) {
      get_values_at_bounds(param->elim->coeffs, ns, rt, tab);
      refine_QIR_positive_root(polelim, &ns, rt, tab, 2 * (rt->k), info_level);
    } else {
      /* root is negative */
      mpz_add_ui(pos_root->numer, rt->numer, 1);
      mpz_neg(pos_root->numer, pos_root->numer);
      pos_root->k = rt->k;
      pos_root->sign_left = -(rt->sign_left);
      pos_root->isexact = rt->isexact;
      for (long i = 0; i <= ns; i++) {
        if ((i & 1) == 1) {
          mpz_neg(polelim[i], polelim[i]);
        }
      }
      get_values_at_bounds(polelim, ns, pos_root, tab);
      refine_QIR_positive_root(polelim, &ns, pos_root, tab,
                               2 * (pos_root->k) + ns, info_level);
      for (long i = 0; i <= ns; i++) {
        if ((i & 1) == 1) {
          mpz_neg(polelim[i], polelim[i]);
        }
      }
      if (pos_root->isexact != 1) {
        rt->k = pos_root->k;
        rt->isexact = pos_root->isexact;
        mpz_add_ui(rt->numer, pos_root->numer, 1);
        mpz_neg(rt->numer, rt->numer);
      } else {
        rt->k = pos_root->k;
        if (rt->isexact != 1) {
          rt->isexact = pos_root->isexact;
          mpz_set(rt->numer, pos_root->numer);
          mpz_neg(rt->numer, rt->numer);
        }
      }
    }

    if (ns != param->nsols) {
      for (long i = 0; i < param->elim->length; i++) {
        mpz_set(polelim[i], param->elim->coeffs[i]);
      }
      ns = param->nsols;
    }

    corr *= 2;
    b *= 2;

    if(rt->isexact){
        mpz_poly_eval_2exp_naive(param->denom->coeffs, param->denom->length - 1,
                        &rt->numer, rt->k, xdo, xup); 
        mpz_set(den_up, xdo[0]);
        mpz_set(den_do, xdo[0]);
        corr = (param->denom->length - 1) * rt->k;
    }
    generate_table_values_full(rt, c, ns, b, corr, xdo, xup);

    if (info_level) {
      fprintf(stderr, "<%ld>", rt->k);
    }
  }

  
  if (rt->isexact == 1) {
    single_exact_real_root_param(param, rt, nb, xdo, xup, den_up, den_do, c,
                                 tmp, val_do, val_up, tab, pt, MAX(prec, rt->k),
                                 info_level);
    return;
  }

  mpz_t v1, v2;
  mpz_init(v1);
  mpz_init(v2);

  for (long nv = 0; nv < param->nvars - 1; nv++) {

    mpz_scalar_product_interval(param->coords[nv]->coeffs,
                                param->coords[nv]->length - 1, rt->k, xdo, xup,
                                tmp, val_do, val_up, corr);
    mpz_neg(val_do, val_do);
    mpz_neg(val_up, val_up);
    mpz_swap(val_up, val_do);

    long dec = prec;

    long nbits = mpz_sizeinbase(val_up, 2) - mpz_sizeinbase(den_up, 2) - mpz_sizeinbase(param->cfs[nv], 2);
    if(nbits <= 0){
        dec = prec - nbits;
    }

    mpz_mul_2exp(val_up, val_up, dec);
    mpz_mul_2exp(val_do, val_do, dec);

    if (rt->isexact==0 && mpz_cmp(val_do, val_up) > 0) {
      fprintf(stderr, "BUG in real root extractor(2)\n");
      exit(1);
    }

    if (mpz_sgn(den_do) >= 0 && mpz_sgn(den_up) >= 0) {
      if (mpz_sgn(val_do) >= 0 && mpz_sgn(val_up) >= 0) {
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_cdiv_q(v2, val_up, tmp);
        
      }
      if (mpz_sgn(val_do) <= 0 && mpz_sgn(val_up) >= 0) {
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_cdiv_q(v2, val_up, tmp);
      }
      if (mpz_sgn(val_do) <= 0 && mpz_sgn(val_up) <= 0) {
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_fdiv_q(v1, val_do, tmp);
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_cdiv_q(v2, val_up, tmp);
      }
    } else {
      if (mpz_sgn(val_do) >= 0 && mpz_sgn(val_up) >= 0) {
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_up, tmp);
        mpz_mul(tmp, den_do, param->cfs[nv]);
        mpz_cdiv_q(v2, val_do, tmp);
      }
      if (mpz_sgn(val_do) <= 0 && mpz_sgn(val_up) >= 0) {
        mpz_mul(tmp, den_up, param->cfs[nv]);
        mpz_fdiv_q(v1, val_up, tmp);
        mpz_cdiv_q(v2, val_do, tmp);
      }
      if (mpz_sgn(val_do) <= 0 && mpz_sgn(val_up) <= 0) {
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

    pt->coords[nv]->k_up = dec;
    pt->coords[nv]->k_do = dec;
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

void normalize_points(real_point_t *pts, int64_t nb, int32_t nv) {

  for (int64_t i = 0; i < nb; i++) {
    for (int32_t j = 0; j < nv; j++) {

      int64_t b = 0;
      while (mpz_cmp_ui(pts[i]->coords[j]->val_up, 0) != 0 &&
             mpz_divisible_2exp_p(pts[i]->coords[j]->val_up, b + 1) != 0) {
        b++;
      }
      b = MIN(b, pts[i]->coords[j]->k_up);
      if (b) {
        mpz_tdiv_q_2exp(pts[i]->coords[j]->val_up, pts[i]->coords[j]->val_up,
                        b);
        pts[i]->coords[j]->k_up -= b;
      }

      b = 0;

      while (mpz_cmp_ui(pts[i]->coords[j]->val_do, 0) != 0 &&
             mpz_divisible_2exp_p(pts[i]->coords[j]->val_do, b + 1) != 0) {
        b++;
      }
      b = MIN(b, pts[i]->coords[j]->k_do);
      if (b) {
        mpz_tdiv_q_2exp(pts[i]->coords[j]->val_do, pts[i]->coords[j]->val_do,
                        b);
        pts[i]->coords[j]->k_do -= b;
      }
    }
  }
}

void extract_real_roots_param(mpz_param_t param, interval *roots, long nb,
                              real_point_t *pts, long prec, long nbits,
                              double step, int info_level) {
  long nsols = param->elim->length - 1;
  mpz_t *xup = malloc(sizeof(mpz_t) * nsols);
  mpz_t *xdo = malloc(sizeof(mpz_t) * nsols);
  mpz_t c, tmp, den_up, den_do, val_up, val_do;
  mpz_init(c);
  mpz_init(tmp);
  mpz_init(den_up);
  mpz_init(den_do);
  mpz_init(val_up);
  mpz_init(val_do);
  for (long i = 0; i < nsols; i++) {
    mpz_init_set_ui(xup[i], 1);
    mpz_init_set_ui(xdo[i], 1);
  }
  mpz_t *tab =
      (mpz_t *)(calloc(8, sizeof(mpz_t))); // table for some intermediate values
  for (int i = 0; i < 8; i++) {
    mpz_init(tab[i]);
    mpz_set_ui(tab[i], 0);
  }

  mpz_t *polelim = calloc(param->elim->length, sizeof(mpz_t));
  for (long i = 0; i < param->elim->length; i++) {
    mpz_init_set(polelim[i], param->elim->coeffs[i]);
  }
  interval *pos_root = calloc(1, sizeof(interval));
  mpz_init(pos_root->numer);
  mpz_t s;
  mpz_init(s);

  double et = realtime();

  for (long nc = 0; nc < nb; nc++) {
    interval *rt = roots + nc;

    lazy_single_real_root_param(param, polelim, rt, nb, pos_root, xdo, xup,
                                den_up, den_do, c, tmp, val_do, val_up, tab,
                                pts[nc], prec, nbits, s, info_level);

    if (info_level) {
      if (realtime() - et >= step) {
        fprintf(stderr, "{%.2f%%}", 100 * nc / ((double)nb));
        et = realtime();
      }
    }
  }

  for (long i = 0; i < nsols; i++) {
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
  for (int i = 0; i < 8; i++)
    mpz_clear(tab[i]);
  free(tab);

  for (long i = 0; i < param->elim->length; i++) {
    mpz_clear(polelim[i]);
  }
  free(polelim);
  mpz_clear(pos_root->numer);
  free(pos_root);

  normalize_points(pts, nb, param->nvars);
}

real_point_t *isolate_real_roots_param(mpz_param_t param, long *nb_real_roots_ptr,
                         interval **real_roots_ptr, int32_t precision,
                         int32_t nr_threads, int32_t info_level) {
  mpz_t *pol = malloc(param->elim->length * sizeof(mpz_t));

  for (long i = 0; i < param->elim->length; i++) {
    mpz_init_set(pol[i], param->elim->coeffs[i]);
  }
  long maxnbits =
      mpz_poly_max_bsize_coeffs(param->elim->coeffs, param->elim->length - 1);

  for (int i = 0; i < param->nvars - 1; i++) {
    long cmax = mpz_poly_max_bsize_coeffs(param->coords[i]->coeffs,
                                          param->coords[i]->length - 1);
    maxnbits = MAX(cmax, maxnbits);
  }
  long prec = MAX(precision, 128 + (maxnbits) / 32);
  double st = realtime();

  long unsigned int nbpos = 0;
  long unsigned int nbneg = 0;
  interval *roots = real_roots(pol, param->elim->length - 1, &nbpos, &nbneg,
                               prec, nr_threads, info_level);
  long nb = nbpos + nbneg;
  double step = (realtime() - st) / (nb) * 10 * LOG2(precision);

  real_point_t *pts = NULL;
  if (info_level > 0) {
    fprintf(stderr, "Number of real roots: %ld\n", nb);
  }
  if (nb) {
    /* */
    if (info_level) {
      fprintf(stderr, "Starts real root extraction.\n");
    }
    double st = realtime();
    pts = malloc(sizeof(real_point_t) * nb);

    for (long i = 0; i < nb; i++) {
      real_point_init(pts[i], param->nvars);
    }

    extract_real_roots_param(param, roots, nb, pts, precision, maxnbits, step,
                             info_level);
    if (info_level) {
      fprintf(stderr, "Elapsed time (real root extraction) = %.2f\n",
              realtime() - st);
    }
  }
  *real_roots_ptr = roots;
  *nb_real_roots_ptr = nb;

  for (long i = 0; i < param->elim->length; i++) {
    mpz_clear(pol[i]);
  }
  free(pol);
  return pts;
}

void isolate_real_roots_lparam(mpz_param_array_t lparams, long **lnbr_ptr,
                               interval ***lreal_roots_ptr,
                               real_point_t ***lreal_pts_ptr, int32_t precision,
                               int32_t nr_threads, int32_t info_level) {
  long *lnbr = malloc(sizeof(long) * lparams->nb);
  interval **lreal_roots = malloc(sizeof(interval *) * lparams->nb);
  real_point_t **lreal_pts = malloc(sizeof(real_point_t *) * lparams->nb);
  for (int i = 0; i < lparams->nb; i++) {
    lreal_roots[i] = NULL;
    lreal_pts[i] = NULL;
  }

  for (int i = 0; i < lparams->nb; i++) {
    lreal_pts[i] =
        isolate_real_roots_param(lparams->params[i], lnbr + i, lreal_roots + i,
                                 precision, nr_threads, info_level);
  }
  (*lnbr_ptr) = lnbr;
  (*lreal_roots_ptr) = lreal_roots;
  (*lreal_pts_ptr) = lreal_pts;
}

int real_msolve_qq(mpz_param_t *mpz_paramp, param_t **nmod_param, int *dim_ptr,
                   long *dquot_ptr, long *nb_real_roots_ptr,
                   interval **real_roots_ptr, real_point_t **real_pts_ptr,
                   data_gens_ff_t *gens,
                   int32_t ht_size, //initial_hts,
		   int32_t unstable_staircase,
                   int32_t nr_threads,
                   int32_t max_nr_pairs,
                   int32_t elim_block_len,
                   int32_t reset_ht,
                   int32_t la_option,
                   int32_t use_signatures,
                   int32_t lift_matrix,
                   int32_t info_level,
                   int32_t print_gb,
                   int32_t pbm_file,
                   int32_t precision,
                   files_gb *files,
                   int round,
                   int32_t get_param){

  /*
    0 if comp. is ok
    1 if comp. failed
    2 if more genericity is required
    -2 if charac is > 0
    -3 if meta data are corrupted
    -4 if bad prime
  */

  double ct0 = cputime();
  double rt0 = realtime();
  int b = msolve_trace_qq(mpz_paramp,
                          nmod_param,
                          dim_ptr,
                          dquot_ptr,
                          gens,
                          ht_size, //initial_hts,
			  unstable_staircase,
                          nr_threads,
                          max_nr_pairs,
                          elim_block_len,
                          reset_ht,
                          la_option,
                          use_signatures,
                          lift_matrix,
                          info_level,
                          print_gb,
                          pbm_file,
                          files,
                          round);
  double ct1 = cputime();
  double rt1 = realtime();

  if(info_level && print_gb == 0){
    /* fprintf( */
    /*     stderr, */
    /*     "Time for rational param: %13.2f (elapsed) sec / %5.2f sec (cpu)\n\n", */
    /*     rt1 - rt0, ct1 - ct0); */
    fprintf (stdout,
	     "\n---------------- TIMINGS ----------------\n");
    fprintf(stdout,
	    "rational param(elapsed) %12.2f sec\n",
	    rt1-rt0);
    fprintf(stdout,
	    "rational param(cpu) %16.2f sec\n",
	    ct1-ct0);
    fprintf(stdout,
	    "-----------------------------------------\n");
  }

  if (get_param > 1) {
    return b;
  }

  if (print_gb) {
    return 0;
  }
  real_point_t *pts = NULL;

  if (b == 0 && *dim_ptr == 0 && *dquot_ptr > 0 && gens->field_char == 0) {

    pts =
        isolate_real_roots_param(*mpz_paramp, nb_real_roots_ptr, real_roots_ptr,
                                 precision, nr_threads, info_level);
    int32_t nb = *nb_real_roots_ptr;
    if (nb) {
      /* If we added a linear form for genericity reasons remove do not
       * return the last (new) variable in the solutions later on */
      if (gens->linear_form_base_coef > 0) {
        for (long i = 0; i < nb; ++i) {
          pts[i]->nvars--;
        }
      }
      /* If we changed the variable order for genericity reasons we have
       * to rechange the entries in the solution points. */
      /* This is to be done only when the parametrization is not requested */
      if (get_param == 0 && gens->change_var_order != -1 &&
          gens->change_var_order != (*mpz_paramp)->nvars - 1 &&
          gens->linear_form_base_coef == 0) {
        coord_t *tmp = malloc(sizeof(coord_t));
        const int32_t nvars = gens->nvars;
        int32_t lidx = gens->change_var_order;
        for (long i = 0; i < nb; ++i) {
          memcpy(tmp, pts[i]->coords[nvars - 1], sizeof(coord_t));
          memcpy(pts[i]->coords[nvars - 1], pts[i]->coords[lidx],
                 sizeof(coord_t));
          memcpy(pts[i]->coords[lidx], tmp, sizeof(coord_t));
        }
        free(tmp);
      }
      *real_pts_ptr = pts;
    }
  }
  return b;
}

void display_arrays_of_real_roots(files_gb *files, int32_t len,
                                  real_point_t **lreal_pts, long *lnbr) {
  if (files->out_file != NULL) {
    FILE *ofile = fopen(files->out_file, "ab+");
    fprintf(ofile, "[");
    for (int i = 0; i < len - 1; i++) {
      display_real_points(ofile, lreal_pts[i], lnbr[i]);
      fprintf(ofile, ", \n");
    }
    display_real_points(ofile, lreal_pts[len - 1], lnbr[len - 1]);
    fprintf(ofile, "];\n");
    fclose(ofile);
  } else {
    fprintf(stdout, "[");
    for (int i = 0; i < len - 1; i++) {
      display_real_points(stdout, lreal_pts[i], lnbr[i]);
      fprintf(stdout, ", \n");
    }
    display_real_points(stdout, lreal_pts[len - 1], lnbr[len - 1]);
    fprintf(stdout, "];\n");
  }
}

void display_output(int b, int dim, int dquot, files_gb *files,
                    data_gens_ff_t *gens, param_t *param,
                    mpz_param_t *mpz_paramp, int get_param,
                    long *nb_real_roots_ptr, interval **real_roots_ptr,
                    real_point_t **real_pts_ptr, int info_level) {
  if (dquot == 0) {
    if (files->out_file != NULL) {
      FILE *ofile = fopen(files->out_file, "ab+");
      fprintf(ofile, "[-1]:\n");
      fclose(ofile);
    } else {
      fprintf(stdout, "[-1]:\n");
    }
    return;
  }

  if (dim == 0 && dquot >= 0) {
    (*mpz_paramp)->nvars = gens->nvars;
    if (files->out_file != NULL) {
      FILE *ofile = fopen(files->out_file, "ab+");
      fprintf(ofile, "[0, ");
      if (get_param >= 1 || gens->field_char) {
        mpz_param_out_str_maple(ofile, gens, dquot, *mpz_paramp, param);
      }
      if (get_param <= 1 && gens->field_char == 0) {
        if (get_param) {
          fprintf(ofile, ",");
        }
        display_real_points(ofile, *real_pts_ptr, *nb_real_roots_ptr);
      }
      fprintf(ofile, "]:\n");
      fclose(ofile);
    } else {
      fprintf(stdout, "[0, ");
      if (get_param >= 1 || gens->field_char) {
        mpz_param_out_str_maple(stdout, gens, dquot, *mpz_paramp, param);
      }
      if (get_param <= 1 && gens->field_char == 0) {
        if (get_param) {
          fprintf(stdout, ",");
        }
        display_real_points(stdout, *real_pts_ptr, *nb_real_roots_ptr);
      }
      fprintf(stdout, "]:\n");
    }
  }
  if (dim > 0) {
    if (info_level > 0) {
      fprintf(stderr, "The ideal has positive dimension\n");
    }
    if (files->out_file != NULL) {
      FILE *ofile2 = fopen(files->out_file, "ab+");
      // 1 because dim is >0
      fprintf(ofile2, "[1, %d, -1, []]:\n", gens->nvars);
      fclose(ofile2);
    } else {
      fprintf(stdout, "[1, %d, -1, []]:\n", gens->nvars);
    }
  }
}

void manage_output(int b, int dim, int dquot, files_gb *files,
                   data_gens_ff_t *gens, param_t *param,
                   mpz_param_t *mpz_paramp, int get_param,
                   long *nb_real_roots_ptr, interval **real_roots_ptr,
                   real_point_t **real_pts_ptr, int info_level) {
  if (b == 0) {
    display_output(b, dim, dquot, files, gens, param, mpz_paramp, get_param,
                   nb_real_roots_ptr, real_roots_ptr, real_pts_ptr, info_level);
  }
  if (b == -2) {
    fprintf(stderr, "Characteristic of the field here shouldn't be positive\n");
    (*mpz_paramp)->dim = -2;
  }
  if (b == -3) {
    fprintf(stderr, "Problem when checking meta data\n");
    (*mpz_paramp)->dim = -3;
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
  int32_t truncate_lifting,
  int32_t get_param,
  int32_t genericity_handling,
  int32_t unstable_staircase,
  int32_t saturate,
  int32_t colon,
  int32_t normal_form,
  int32_t normal_form_matrix,
  int32_t is_gb,
  int32_t lift_matrix,
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
                    0 /*truncate_lifting */, info_level);

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

            if(check_ff_bits(gens->field_char) < 32){
              fprintf(stderr, "Error: not implemented yet (prime field of too low characteristic)\n");
              return 1;
            }
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
                    gens->field_char, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, saturate, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    0 /*truncate_lifting */, info_level);

            if (!success) {
                printf("Bad input data, stopped computation.\n");
                exit(1);
            }

            st->gfc = gens->field_char;

            if (is_gb == 1) {
                for (len_t k = 0; k < bs->ld; ++k) {
                    bs->lmps[k] = k;
                    bs->lm[k]   = bht->hd[bs->hm[k][OFFSET]].sdm;
                    bs->lml     = bs->ld;
                }
            } else {
                sat = initialize_basis(st);
                sat->ht = bht;
                import_input_data(sat, st, gens->ngens-saturate, gens->ngens, gens->lens, gens->exps, (void *)gens->cfs, NULL);

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
		  print_ff_basis_data(files->out_file, "ab", bs, bht,
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
                    0 /*truncate_lifting */, info_level);

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
            tbr->ht = bht;
            import_input_data(tbr, st, gens->ngens-1, gens->ngens, gens->lens, gens->exps, (void *)gens->cfs, NULL);
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
	    nvars_t *linvars = calloc(gens->nvars, sizeof(nvars_t));
	    uint32_t *lineqs = calloc(gens->nvars,sizeof(uint32_t));
	    nvars_t *squvars = calloc(gens->nvars-1, sizeof(nvars_t));
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
          fprintf(stderr, "Warning: elim order not available for rational parametrizations\n");
          fprintf(stderr, "Computing Groebner basis\n");
          print_gb=2;
      }
	  b = real_msolve_qq(mpz_paramp,
                       paramp,
                       &dim,
                       &dquot,
                       nb_real_roots_ptr,
                       real_roots_ptr,
                       real_pts_ptr,
                       gens,
		       initial_hts, unstable_staircase, nr_threads, max_pairs,
                       elim_block_len, update_ht,
                       la_option, use_signatures, lift_matrix, info_level, print_gb,
                       generate_pbm, precision, files, round, get_param);
          if(print_gb){
            return 0;
          }

          manage_output(b, dim, dquot, files, gens, (*paramp), mpz_paramp, get_param,
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
	    if (genericity_handling == 2) {
	      if (add_random_linear_form_to_input_system(gens, info_level)) {
		goto restart;
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
        }
	else {
          /* normal_form is 1 */
            /* data structures for basis, hash table and statistics */
            bs_t *bs    = NULL;
            bs_t *tbr   = NULL;
            ht_t *bht   = NULL;
            md_t *st  = NULL;


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
                    gens->field_char, 0 /* DRL order */, elim_block_len, gens->nvars,
                    /* gens->field_char, 0 [> DRL order <], gens->nvars, */
                    gens->ngens, normal_form, initial_hts, nr_threads, max_pairs,
                    update_ht, la_option, use_signatures, 1 /* reduce_gb */, 0,
                    0 /*truncate_lifting */, info_level);

            st->gfc  = gens->field_char;
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
            /* initialize data for elements to be reduced,
             * NOTE: Don't initialize BEFORE running core_f4, bht may
             * change, so hash values of tbr may become wrong. */
            tbr = initialize_basis(st);
            tbr->ht = bht;
            import_input_data(tbr, st, gens->ngens-normal_form, gens->ngens,
                    gens->lens, gens->exps, (void *)gens->cfs, NULL);
            tbr->ld = tbr->lml  =  normal_form;
            /* normalize_initial_basis(tbr, st->gfc); */
            for (int k = 0; k < normal_form; ++k) {
                tbr->lmps[k]  = k; /* fix input element in tbr */
            }

            /* generate array for storing multiplier for polynomial
             * to be reduced by basis */
            exp_t *mul  = (exp_t *)calloc(bht->evl, sizeof(exp_t));
            /* compute normal form of last element in tbr */
            tbr = core_nf(tbr, st, mul, bs, &err);

            if (err) {
                printf("Problem with normalform, stopped computation.\n");
                exit(1);
            }
            /* print all reduced elements in tbr, first normal_form ones
             * are the input elements */
            print_ff_nf_data(files->out_file, "ab", 0, normal_form, tbr, bht, st, gens, 2);
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
                free_basis_without_hash_table(&bs);
            }
            if (tbr != NULL) {
                free_basis(&tbr);
            }
            free(st);
            st  = NULL;

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
            import_input_data(bs_qq, st, 0, st->ngens_input, gens->lens, gens->exps,
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

            ht_t *lht = copy_hash_table(bht);

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
            const uint32_t prime_start  = (uint32_t)(1) << 30;
            const int32_t nr_primes     = nr_threads;

            /* data structures for basis, hash table and statistics */
            bs_t *sat_qq   = NULL;

            /* initialize stuff */
            md_t *st  = allocate_meta_data();

            int *invalid_gens       =   NULL;
            int32_t monomial_order  =   0;
            int32_t reduce_gb       =   1;
            int32_t truncate_lifting =  0;
            int res = validate_input_data(&invalid_gens, gens->mpz_cfs,
                    gens->lens, &field_char, &monomial_order, &elim_block_len,
                    &gens->nvars, &gens->ngens, &saturate, &initial_hts,
                    &nr_threads, &max_pairs, &update_ht, &la_option,
                    &use_signatures, &reduce_gb, &truncate_lifting, &info_level);

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
                        1, prime_start, nr_primes, 0, truncate_lifting, 
                        info_level)) {
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
            import_input_data(bs_qq, st, 0, st->ngens_input, gens->lens, gens->exps,
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
            sat_qq->ht = bht;
            import_input_data(
                    sat_qq, st, gens->ngens-saturate, gens->ngens,
                    gens->lens, gens->exps, (void *)gens->mpz_cfs, NULL);
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

            ht_t *lht = copy_hash_table(bht);

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
                    if (info_level) {
                        printf("F4 trace timing %13.2f\n", stf4);
                    }
                    /* printf("bs[%u]->lml = %u\n", i, bs[i]->lml); */
                }
            /* } */
            return 0;
      } else {       /* characteristic is 0 and elim_block = 0 and saturate = 0 */

            int dim = - 2;
            long dquot = -1;

            if(elim_block_len && print_gb == 0){
                fprintf(stderr, "Warning: elim order not available for rational parametrizations\n");
                fprintf(stderr, "Computing Groebner basis\n");
                print_gb=2;
            }

            if(print_gb){

              msflags_t flags;

              flags->ht_size = initial_hts;
              flags->nr_threads = nr_threads;
              flags->max_nr_pairs = max_pairs;
              flags->elim_block_len = elim_block_len;
              flags->truncate_lifting = truncate_lifting;
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

            b = real_msolve_qq(mpz_paramp,
                    &param,
                    &dim,
                    &dquot,
                    nb_real_roots_ptr,
                    real_roots_ptr,
                    real_pts_ptr,
                    gens,
		    initial_hts, unstable_staircase, nr_threads, max_pairs,
                    elim_block_len, update_ht,
                    la_option, use_signatures, lift_matrix, info_level, print_gb,
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
                if (genericity_handling == 2) {
                  if (add_random_linear_form_to_input_system(gens, info_level)) {
                    goto restart;
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

void export_julia_rational_parametrization_qq(
    void *(*mallocp)(size_t), int32_t *load, int32_t *nvars, int32_t *dim,
    int32_t *dim_quot, int32_t **lens, char ***var_namesp,
    void **cfs_linear_form, void **cfs, void **real_sols_num,
    int32_t **real_sols_den,
    data_gens_ff_t *gens, /* might change vnames, thus not const */
    const mpz_param_t param, const long nb_real_roots,
    const real_point_t *real_pts) {
  int32_t i, j;
  int64_t ctr = 0;

  *load = (int32_t)param->nvars + 1;
  *dim = (int32_t)param->dim;
  *dim_quot = (int32_t)param->dquot;
  *nvars = (int32_t)gens->nvars;

  /* keep variable names for returning rational parametrization */
  *var_namesp = gens->vnames;
  gens->vnames = NULL;

  mpz_t *cf_lf = NULL;
  /* check existence of linear form */
  if (gens->linear_form_base_coef > 0) {
    cf_lf = (mpz_t *)(*mallocp)((unsigned long)(gens->nvars) * sizeof(mpz_t));
    int64_t len = 0;
    for (i = 0; i < gens->ngens - 1; ++i) {
      len += 2 * gens->lens[i]; /* numerator plus denominator, thus "2*" */
    }
    j = 0;
    for (i = 0; i < 2 * gens->nvars; i += 2) {
      mpz_init_set(cf_lf[j++], *(gens->mpz_cfs[i + len]));
    }
  }

  if ((param->dim > 0) || (param->dim == 0 && param->dquot == 0)) {
    *lens = NULL;
    *cfs = NULL;
  } else {
    int32_t *len = (int32_t *)(*mallocp)((unsigned long)(param->nvars + 1) *
                                         sizeof(int32_t));

    /* precompute number of all terms of all polynomials */
    int64_t nterms = 0;
    nterms += (int64_t)param->elim->length;
    len[0] = (int32_t)param->elim->length;
    nterms += (int64_t)param->denom->length;
    len[1] = (int32_t)param->denom->length;
    /* In the following we add one more term since each coords[i] polynomial
     * also gets another coefficient for the correct denominator in the
     * rational parametrization of the solution set. */
    for (i = 0; i < param->nvars - 1; ++i) {
      nterms += param->coords[i]->length + 1;
      len[i + 2] = param->coords[i]->length + 1;
    }

    mpz_t *cf = (mpz_t *)(*mallocp)((unsigned long)(nterms) * sizeof(mpz_t));

    /* store elim */
    for (i = 0; i < param->elim->length; ++i) {
      mpz_init_set((cf + ctr)[i], param->elim->coeffs[i]);
    }
    ctr += param->elim->length;
    /* store denom */
    for (i = 0; i < param->denom->length; ++i) {
      mpz_init_set((cf + ctr)[i], param->denom->coeffs[i]);
    }
    ctr += param->denom->length;

    /* store param */
    for (i = 0; i < param->nvars - 1; ++i) {
      for (j = 0; j < param->coords[i]->length; ++j) {
        mpz_init_set((cf + ctr)[j], param->coords[i]->coeffs[j]);
      }
      mpz_init_set((cf + ctr)[j], param->cfs[i]);
      ctr += param->coords[i]->length + 1;
    }
    *lens = len;
    *cfs = (void *)cf;
    *cfs_linear_form = (void *)cf_lf;

    /* if there are no real solutions return the parametrization at least */
    if (nb_real_roots <= 0) {
      return;
    }

    const long nb_real_roots_intervall = 2 * nb_real_roots;

    mpz_t *sols_num =
        (mpz_t *)(*mallocp)((unsigned long)nb_real_roots_intervall *
                            real_pts[0]->nvars * sizeof(mpz_t));

    int32_t *sols_den =
        (int32_t *)(*mallocp)((unsigned long)nb_real_roots_intervall *
                              real_pts[0]->nvars * sizeof(int32_t));

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
        gens->vnames[i] = calloc(strlen(var_names[i]) + 1, sizeof(char));
        memcpy(gens->vnames[i], var_names[i], (strlen(var_names[i]) + 1) * sizeof(char));
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
                          0 /* generate pbm */, 1 /* reduce_gb */,
                          print_gb, 0 /*truncate_lifting*/, get_param,
			  genericity_handling, 0 /* unstable_staircase -> change to 2?*/,
			  0 /* saturate */, 0 /* colon */,
			  0 /* normal_form */, 0 /* normal_form_matrix */,
			  0 /* is_gb */, 0 /*lift_matrix */, precision, files,
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
    if(param != NULL && gens->field_char){
        free_fglm_param(param);
    }
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
        fprintf(stderr, "\n-------------------------------------------------\
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
void free_msolve_julia_result_data(void (*freep)(void *), int32_t **res_len,
                                   void **res_cf, void **sols_num,
                                   int32_t **sols_den, const int64_t res_ld,
                                   const int64_t nr_sols,
                                   const int64_t field_char) {

  int32_t *lens = *res_len;

  (*freep)(lens);
  lens = NULL;
  *res_len = lens;

  if (field_char > 0) {
    int32_t *numerators = *(int32_t **)sols_num;
    (*freep)(numerators);
    numerators = NULL;
    int32_t *cfs = *(int32_t **)res_cf;
    (*freep)(cfs);
    cfs = NULL;
  } else {
    /* denominators only exist if char == 0 */
    int32_t *denominators = *sols_den;
    (*freep)(denominators);
    denominators = NULL;
    *sols_den = denominators;
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
  *res_cf = NULL;
}
