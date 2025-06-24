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

#include "getdelim.h"

static inline void store_exponent(const char *term, data_gens_ff_t *gens, int32_t pos)
{
    len_t i, j, k;

    len_t op = 0;
    char *var = NULL;
    char *ev  = NULL;
    for (i = 0; i < strlen(term)+1; ++i) {
        if (term[i] == '*' || i == strlen(term)) {
            j = op;
            while (j < i && term[j] != '^') {
                ++j;
            }
            if (term[j-1] == ',') {
                --j;
            }
            while(term[op] == ' ' || term[op] == '+' || term[op] == '-') {
                ++op;
            }
            var = realloc(var, sizeof(char)*(j-op+1));
            memcpy(var, term+op, j-op);
            var[j-op] = '\0';
            if (term[j] == '^') {
                ev = realloc(ev, (sizeof(char)*(i-j)));
                ev = memcpy(ev, term+(j+1), i-j-1);
                ev[i-j-1] = '\0';
            } else {
                ev = realloc(ev, (sizeof(char)*2));
                ev[0] = '1';
                ev[1] ='\0';
            }
            for (k = 0; k < gens->nvars; ++k) {
                if (strcmp(gens->vnames[k], var) == 0) {
                    ((gens->exps) + pos)[k] = strtol(ev, NULL, 10);
                    break;
                }
            }
            op = i+1;
        }
    }
    free(var);
    free(ev);
}



static inline void display_monomial(FILE *file, data_gens_ff_t *gens, int64_t pos,
                                    int32_t **bexp){
  int32_t exp = 0;
  for(int k = 0; k < gens->nvars - gens->elim; k++){
    exp = (*bexp)[(pos)*(gens->nvars - gens->elim) + k];
    if(exp > 0){
      break;
    }
  }
  if(exp == 0){
    fprintf(file, "1");
    return;
  }
  for(int k = 0; k < (gens->nvars - gens->elim); k++){
    exp = (*bexp)[(pos)*(gens->nvars - gens->elim) + k];
    if(exp > 0){
      fprintf(file, "*");
      if(exp==1){
        fprintf(file, "%s", gens->vnames[k + gens->elim]);
      }
      else{
        fprintf(file, "%s^%d",gens->vnames[k + gens->elim], exp);
      }
    }
  }
}

static inline int32_t display_monomial_single(FILE *file, data_gens_ff_t *gens,
                                              int64_t pos, int32_t **bexp){
  int32_t exp = 0;
  int32_t b = 0;
  for(int k = 0; k < gens->nvars - gens->elim; k++){
    exp = (*bexp)[(pos)*(gens->nvars - gens->elim) + k];
    if(exp > 0){
      break;
    }
  }
  if(exp == 0){
    fprintf(file, "1");
    return 0;
  }
  for(int k = 0; k < gens->nvars - gens->elim; k++){
    exp = (*bexp)[(pos)*(gens->nvars - gens->elim) + k];
    if(exp > 0){
      if(exp==1){
        if(b){
          fprintf(file, "*%s", gens->vnames[k + gens->elim]);
        }
        else{
          fprintf(file, "%s", gens->vnames[k + gens->elim]);
        }
      }
      else{
        if(b){
          fprintf(file, "*%s^%d",gens->vnames[k + gens->elim], exp);
        }
        else{
          fprintf(file, "%s^%d",gens->vnames[k + gens->elim], exp);
        }
      }
      b = 1;
    }
  }
  return b;
}


static inline int32_t display_monomial_full(FILE *file, const int nv,
                                            char **vnames,
                                            int64_t pos, int32_t *bexp){
  int32_t exp, b = 0;
  for(int k = 0; k < nv; k++){
    exp = (bexp)[(pos)*nv + k];
    if(exp > 0){
      if(exp==1){
        if(b){
          //          fprintf(file, "*%s", vnames[k]);
          fprintf(file, "*x%d", k+1);
        }
        else{
          //          fprintf(file, "%s", vnames[k]);
          fprintf(file, "x%d", k+1);
        }
      }
      else{
        if(b){
          fprintf(file, "*x%d^%d",k+1, exp);
        }
        else{
          fprintf(file, "x%d^%d",k+1, exp);
        }
      }
      b = 1;
    }
  }
  if(b==0){
    fprintf(file, "1");
  }
  return b;
}

static void print_msolve_polynomials_ff(
        FILE *file,
        const bi_t from,
        const bi_t to,
        const bs_t * const bs,
        const ht_t * const ht,
        const md_t *st,
        char **vnames,
        const int lead_ideal_only,
        const int is_nf
        )
{
    len_t i, j, k, idx;

    len_t len   = 0;
    hm_t *hm    = NULL;

    const len_t nv  = ht->nv;
    const len_t ebl = ht->ebl;
    const len_t evl = ht->evl;

    /* state context if full basis is printed */
    if (is_nf == 0 && from == 0 && to == bs->lml) {
        if (lead_ideal_only != 0) {
            fprintf(file, "#Leading ideal data\n");
        } else {
            fprintf(file, "#Reduced Groebner basis data\n");
        }
        fprintf(file, "#---\n");
        fprintf(file, "#field characteristic: %u\n", st->gfc);
        fprintf(file, "#variable order:       ");
        for (i = 0; i < nv-1; ++i) {
            fprintf(file, "%s, ", vnames[i]);
        }
        fprintf(file, "%s\n", vnames[nv-1]);
        if (st->nev == 0) {
            fprintf(file, "#monomial order:       graded reverse lexicographical\n");
        } else {
            if (st->nev == 1) {
                fprintf(file, "#monomial order:       eliminating first variable, blocks: graded reverse lexicographical\n");
            } else {
                fprintf(file, "#monomial order:       eliminating first %d variables, blocks: graded reverse lexicographical\n", st->nev);
            }
        }
        if (bs->lml == 1) {
            fprintf(file, "#length of basis:      1 element\n");
        } else {
            fprintf(file, "#length of basis:      %u elements sorted by increasing leading monomials\n", bs->lml);
        }
        fprintf(file, "#---\n");
    }

    int *evi    =   (int *)malloc((unsigned long)ht->nv * sizeof(int));
    if (ebl == 0) {
        for (i = 1; i < evl; ++i) {
            evi[i-1]    =   i;
        }
    } else {
        for (i = 1; i < ebl; ++i) {
            evi[i-1]    =   i;
        }
        for (i = ebl+1; i < evl; ++i) {
            evi[i-2]    =   i;
        }
    }

    if (lead_ideal_only != 0) {
        int ctr = 0;
        fprintf(file, "[");
        for (i = from; i < to; ++i) {
            idx = bs->lmps[i];
            if (bs->hm[idx] == NULL) {
                fprintf(file, "0,\n");
            } else {
                hm  = bs->hm[idx]+OFFSET;
                len = bs->hm[idx][LENGTH];
                ctr = 0;
                k = 0;
                while (ctr == 0 && k < nv) {
                    if (ht->ev[hm[0]][evi[k]] > 0) {
                        fprintf(file, "%s^%u",vnames[k], ht->ev[hm[0]][evi[k]]);
                        ctr++;
                    }
                    k++;
                }
                for (;k < nv; ++k) {
                    if (ht->ev[hm[0]][evi[k]] > 0) {
                        fprintf(file, "*%s^%u",vnames[k], ht->ev[hm[0]][evi[k]]);
                    }
                }
                if (i < to-1) {
                    fprintf(file, ",\n");
                } else {
                    fprintf(file, "]:\n");
                }
            }
        }
    } else {
        fprintf(file, "[");
        for (i = from; i < to; ++i) {
            idx = bs->lmps[i];
            if (bs->hm[idx] == NULL) {
                fprintf(file, "0");
            } else {
                hm  = bs->hm[idx]+OFFSET;
                len = bs->hm[idx][LENGTH];
                switch (st->ff_bits) {
                case 0 :
                  fprintf(file, "%u", bs->cf_32[bs->hm[idx][COEFFS]][0]);
                  break;
                case 8:
                  fprintf(file, "%u", bs->cf_8[bs->hm[idx][COEFFS]][0]);
                  break;
                case 16:
                  fprintf(file, "%u", bs->cf_16[bs->hm[idx][COEFFS]][0]);
                  break;
                case 32:
                  fprintf(file, "%u", bs->cf_32[bs->hm[idx][COEFFS]][0]);
                  break;
                }
                for (k = 0; k < nv; ++k) {
                    if (ht->ev[hm[0]][evi[k]] > 0) {
                        fprintf(file, "*%s^%u",vnames[k], ht->ev[hm[0]][evi[k]]);
                    }
                }
                for (j = 1; j < len; ++j) {
		  switch (st->ff_bits) {
		  case 0:
		    fprintf(file, "+%u", bs->cf_32[bs->hm[idx][COEFFS]][j]);
		    break;
		  case 8:
		    fprintf(file, "+%u", bs->cf_8[bs->hm[idx][COEFFS]][j]);
		    break;
		  case 16:
		    fprintf(file, "+%u", bs->cf_16[bs->hm[idx][COEFFS]][j]);
		    break;
		  case 32:
		    fprintf(file, "+%u", bs->cf_32[bs->hm[idx][COEFFS]][j]);
		    break;
		  }
		  for (k = 0; k < nv; ++k) {
		    if (ht->ev[hm[j]][evi[k]] > 0) {
		      fprintf(file, "*%s^%u",vnames[k], ht->ev[hm[j]][evi[k]]);
		    }
		  }
                }
	    }
	    if (i < to-1) {
	      fprintf(file, ",\n");
	    } else {
	      fprintf(file, "]:\n");
	    }
        }
    }
    free(evi);
}

static void print_ff_nf_data(
        const char *fn,
        const char *mode,
        const bi_t from,
        const bi_t to,
        const bs_t *bs,
        const ht_t *ht,
        const md_t *st,
        const data_gens_ff_t *gens,
        const int32_t print_gb)
{
    if (print_gb > 0) {
        if(fn != NULL){
            FILE *ofile = fopen(fn, mode);
            print_msolve_polynomials_ff(ofile, from, to, bs, ht,
                    st, gens->vnames, 2-print_gb, 1);
            fclose(ofile);
        }
        else{
            print_msolve_polynomials_ff(stdout, from, to, bs, ht,
                    st, gens->vnames, 2-print_gb, 1);
        }
    }
}

static void print_ff_basis_data(
        const char *fn,
        const char *mode,
        const bs_t *bs,
        const ht_t *ht,
        const md_t *st,
        const data_gens_ff_t *gens,
        const int32_t print_gb)
{
    if (print_gb > 0) {
        if(fn != NULL){
            FILE *ofile = fopen(fn, mode);
            print_msolve_polynomials_ff(ofile, 0, bs->lml, bs, ht,
                    st, gens->vnames, 2-print_gb, 0);
            fclose(ofile);
        }
        else{
            print_msolve_polynomials_ff(stdout, 0, bs->lml, bs, ht,
                    st, gens->vnames, 2-print_gb, 0);
        }
    }
}

static int32_t get_nvars(const char *fn)
{
    FILE * fh = fopen(fn, "r");
    char * line = NULL;
    char * line2 = NULL;
    size_t len;
    nvars_t nvars = -1; 

    /* number of variables is read from first line, it is 1 + (number of commata) */
    if (getline(&line, &len, fh) != -1)
    {
        nvars = 1;
        line2 = strchr(line, ',');
        while (line2 != NULL)
        {
            nvars++;
            // line points to a comma, which must be followed by one or more characters
            // --> line+1 is valid
            line2 = strchr(line2+1, ',');
        }
    }

    free(line);
    free(line2);
    fclose(fh);

    return nvars;
}

/**
 * Checks if a null-or-delim-terminated string is just empty
 * ("only whitespaces" counts as empty)
 *
 * \return 1 if the line is empty, else 0
 */
static inline int is_line_empty(const char *line, const char delim) {
    while (*line != '\0' && *line != ',') {
        if (!isspace(*line))
            return 0;
        line++;
    }
    return 1;
}

static int32_t get_ngenerators(char *fn)
{
    int32_t ngens = 0;
    char *line  = NULL;
    size_t len;
    FILE *fh  = fopen(fn,"r");

    /* 1st and 2nd lines are ignored; still check all went fine */
    if (   getline(&line, &len, fh) == -1
        || getline(&line, &len, fh) == -1)
        ngens = -1;

    /* go through subsequent lines, not counting empty ones */
    else
        while(getdelim(&line, &len, ',', fh) != -1)
            if (! is_line_empty(line, ','))
                ngens++;

    free(line);
    fclose(fh);

    return ngens;
}


static inline char *get_variable_name(const char *line, char **prev_pos)
{
  const char comma_splicer  = ',';

  char *tmp_var   = (char *)malloc(50 * sizeof(char));
  char *curr_pos  = strchr(*prev_pos, comma_splicer);
  if (curr_pos != NULL) {
    size_t var_diff   = (size_t)(curr_pos - *prev_pos);
    memcpy(tmp_var, *prev_pos, var_diff);
    tmp_var[var_diff] = '\0';
    *prev_pos = curr_pos+1;
  } else { /** we are at the last variable */
    int prev_idx    = (int)(*prev_pos - line);
    int curr_idx    = (int)(strlen(line)+1);
    size_t var_diff = (size_t)(curr_idx - prev_idx);
    memcpy(tmp_var, *prev_pos, var_diff);
    tmp_var[var_diff] = '\0';
  }

  /** trim variable, remove blank spaces */
  char *tmp_var_begin, *tmp_var_end;
  tmp_var_begin = tmp_var_end =  tmp_var;
  while (isspace(*tmp_var_begin))
    tmp_var_begin++;
  if (*tmp_var_begin != 0) {
    while (*tmp_var_end)
      tmp_var_end++;
    do {
      tmp_var_end--;
    } while (isspace(*tmp_var_end));
  }
  int fin_len  = (int)(tmp_var_end - tmp_var_begin + 1);
  size_t final_length = 0;
  if (fin_len > 0)
    final_length  = (size_t)fin_len;
  char *final_var   = (char *)malloc((final_length+1) * sizeof(char));
  memcpy(final_var, tmp_var_begin, final_length);
  final_var[final_length] = '\0';

  free(tmp_var);
  return final_var;
}

static inline nelts_t get_number_of_terms(const char *line)
{
  const char add_splicer        = '+';
  const char minus_splicer      = '-';
  const char whitespace_splicer = ' ';
  char *tmp = NULL;
  nelts_t nterms  = 1;
  /** remove useless whitespaces at the beginning */
  int i = 0;
  while (strncmp(&whitespace_splicer, line+i, 1) == 0)
    i++;
  /** check if first non-whitespace char is "-", set term counter -1 in this case */
  if (strncmp(&minus_splicer, line+i, 1) == 0)
    nterms--;
  /** now count terms */
  tmp = strchr(line, add_splicer);
  while (tmp != NULL) {
    nterms++;
    tmp = strchr(tmp+1, add_splicer);
  }
  tmp = strchr(line, minus_splicer);
  while (tmp != NULL) {
    nterms++;
    tmp = strchr(tmp+1, minus_splicer);
  }
  //MS que faire si dans la ligne on a 0*x+0*y ?
  return nterms;
}


static void get_variables(FILE *fh, char * line, int max_line_size,
                          int32_t *nr_vars, data_gens_ff_t *gens, char **vnames){
  long i;
  char *tmp = NULL;
  if (fgets(line, max_line_size, fh) != NULL) {
    tmp = line;
  } else {
    printf("Bad file format (variable names).\n");
    /* free(vnames);
     * free(line); */
    fclose(fh);
    exit(1);
  }
  for (i=0; i<(*nr_vars); ++i) {
    vnames[i]  = get_variable_name(line, &tmp);
  }
  gens->vnames = vnames;

}

static void get_characteristic(FILE *fh, char * line, int max_line_size,
                               int32_t *field_char, char **vnames){
  *field_char = 0;

  if(fgets(line, max_line_size, fh) != NULL){
    int64_t tmp_mod = atol(line);
    if(tmp_mod >= 0){
      *field_char = (int32_t) tmp_mod;
      if(tmp_mod > 2147483647){
        fprintf(stderr, "Warning: characteristic must be 0 or < 2^31\n");
        free(vnames);
        free(line);
        fclose(fh);
        exit(1);
      }
    }
    else{
      fprintf(stderr, "Bad file format (characteristic)\n");
      free(vnames);
      free(line);
      fclose(fh);
      exit(1);
    }
  }

}

static void trim_polynomial(char *poly) {
    /*
     * WW: The function trims the input, supposedly a polynomial,
     * in a possibly too trivial way, but is left as is for
     * backward compatibility.
     *
     * For example, trim_polynomial converts "x 1, pi/2" to "x1",
     * even when the former is not a valid polynomial.
     * This actually mimics the behavior of msolve v0.7.5,
     * for which the following is a valid input:
     *   x1
     *   0
     *   x 1, pi/2
     *
     * Ideally, in a future version, we should throw an error if the
     * input is not a valid polynomial.
     */
    size_t k = 0;
    for (size_t j = 0; poly[j] != '\0'; j++) {
        if (poly[j] == ',') {
            break;
        }
        if (!isspace(poly[j])) {
            poly[k++] = poly[j];
        }
    }
    poly[k] = '\0';
}

static void remove_newlines(char *line, size_t *len) {
    while (*len && line[*len - 1] == '\0') {
        (*len)--;
    }
    while (*len && isspace(line[*len - 1])) {
        line[--(*len)] = '\0';
    }
    for (size_t i = 0; i < *len; i++) {
        if (line[i] == '\n' || line[i] == '\r') {
            line[i] = ' ';
        }
    }
}

static void remove_trailing_delim(char *line, size_t *len, char delim) {
    if (*len > 0 && line[*len - 1] == delim) {
        line[--(*len)] = '\0';
    }
}

static void get_nterms_and_all_nterms(FILE *fh, 
                                      int max_line_size, data_gens_ff_t *gens,
                                      int32_t *nr_gens, nelts_t *nterms, nelts_t *all_nterms){

    char *line  = NULL; 
    size_t len = 0;
    ssize_t nread;
    for (int32_t i = 0; i < *nr_gens; i++) {
        do {
            nread = getdelim(&line, &len, ',', fh);
        } while (nread != -1 && is_line_empty(line, ','));
        if (i < *nr_gens - 1) {
            remove_trailing_delim(line, &len, ',');
        }
        remove_newlines(line, &len);
        trim_polynomial(line);
        *nterms  = get_number_of_terms(line);
        gens->lens[i] = *nterms;
        *all_nterms += *nterms;
    }
    free(line);
    gens->nterms = *all_nterms;
}


static inline void get_term(const char *line, char **prev_pos,
                            char **term, size_t *term_size)
{
  /** note that maximal term length we handle */
  const char add_splicer    = '+';
  const char minus_splicer  = '-';

  char *start_pos;
  char *curr_pos_add    = strchr(*prev_pos, add_splicer);
  char *curr_pos_minus  = strchr(*prev_pos, minus_splicer);
  if (curr_pos_minus == line)
    curr_pos_minus  = strchr(*prev_pos+1, minus_splicer);

  if (*prev_pos != line)
    start_pos = *prev_pos - 1;
  else
    start_pos = *prev_pos;

  if (curr_pos_add != NULL && curr_pos_minus != NULL) {
    size_t term_diff_add   = (size_t)(curr_pos_add - start_pos);
    size_t term_diff_minus = (size_t)(curr_pos_minus - start_pos);
    /** if "-" is the first char in the line, we have to adjust the
      * if minus is nearer */
    if (term_diff_add > term_diff_minus) {
      if( term_diff_minus > *term_size){
        fprintf(stderr, "Too large input integers, exit...\n");
        exit(1);
      }
      memcpy(*term, start_pos, term_diff_minus);
      (*term)[term_diff_minus]  = '\0';
      *prev_pos                 = curr_pos_minus+1;
      return;
    /** if plus is nearer */
    } else {
      if(term_diff_add > *term_size){
        fprintf(stderr, "Too large input integers, exit...\n");
        exit(1);
      }
      memcpy(*term, start_pos, term_diff_add);
      (*term)[term_diff_add]  = '\0';
      *prev_pos               = curr_pos_add+1;
      return;
    }
  } else {
    if (curr_pos_add != NULL) {
      size_t term_diff_add   = (size_t)(curr_pos_add - start_pos);
      if(term_diff_add > *term_size){
        fprintf(stderr, "Too large input integers, exit...\n");
        exit(1);
      }
      memcpy(*term, start_pos, term_diff_add);
      (*term)[term_diff_add]  = '\0';
      *prev_pos               = curr_pos_add+1;
      return;
    }
    if (curr_pos_minus != NULL) {
      size_t term_diff_minus = (size_t)(curr_pos_minus - start_pos);
      if(term_diff_minus > *term_size){
        fprintf(stderr, "Too large input integers, exit...\n");
        exit(1);
      }
      memcpy(*term, start_pos, term_diff_minus);
      (*term)[term_diff_minus]  = '\0';
      *prev_pos                 = curr_pos_minus+1;
      return;
    }
    if (curr_pos_add == NULL && curr_pos_minus == NULL) {
      size_t prev_idx  = (size_t)(start_pos - line);
      size_t term_diff = strlen(line) + 1 - prev_idx;
      if(term_diff > *term_size){
        fprintf(stderr, "Too large input integers, exit...\n");
        exit(1);
      }
      memcpy(*term, start_pos, term_diff);
      (*term)[term_diff]  = '\0';
      return;
    }
  }
}


/*assumes that coeffs in file fit in word size */
static int get_coefficient_ff_and_term_from_line(char *line, int32_t nterms,
                                          int32_t field_char,
                                          data_gens_ff_t *gens, int32_t pos){
  char *prev_pos = NULL;
  size_t term_size = 50000;
  char *term  = (char *)malloc(term_size * sizeof(char));
  int64_t cf_tmp  = 0; /** temp for coefficient value, possibly coeff is negative. */

  prev_pos = line;
  get_term(line, &prev_pos, &term, &term_size);
  if(term != NULL){
    int32_t iv_tmp  = (int32_t)strtol(term, NULL, 10);
    if (iv_tmp == 0) {
      switch (term[0]) {
          case '0':
            iv_tmp = 0;
            break;
          case '-':
            iv_tmp = -1;
            break;
          default:
            iv_tmp = 1;
            break;
      }
    }
    while (iv_tmp < 0) {
      iv_tmp  +=  field_char; //MS change int -> long int
    }
    gens->cfs[pos]  = (int32_t)iv_tmp;
    store_exponent(term, gens, pos*gens->nvars);
    for(int j = 1; j < nterms; j++){
      get_term(line, &prev_pos, &term, &term_size);
      if (term != NULL) {
        cf_tmp  = (int64_t)strtol(term, NULL, 10);

        if (cf_tmp == 0) {
          if (term[0] == '-') {
            cf_tmp = -1;
          } else {
            cf_tmp = 1;
          }
        }
        while (cf_tmp < 0) {
          cf_tmp  += field_char;
        }
        gens->cfs[pos+j] = (int32_t)(cf_tmp % field_char);
        store_exponent(term, gens, (pos+j)*gens->nvars);
      }
      //      store_exponent(term, basis, ht);
    }
    free(term);
    return 0;
  }
  free(term);
  return 1;
}

static void beginning_strterm_to_mpz(char *str, mpz_t *num, mpz_t *den){
  /* fprintf(stderr, "TERM = %s\n", str); */
  mpq_t tmp;
  mpq_init(tmp);
  mpz_set_ui(*num, 1);
  mpz_set_ui(*den, 1);
  size_t m = 0;
  size_t len = strlen(str);
  for(m = 0; m < len; m++){
    if(str[m]=='*' || str[m]==',') break;
  }
  const size_t sz = m;
  char *buf = (char *)malloc((sz+1) * sizeof(char));
  for(long i = 0; i < sz; i++){
    buf[i]=str[i];
  }
  buf[sz]='\0';

  int b = mpq_set_str(tmp, buf, 10);
  mpq_get_num(*num, tmp);
  mpq_get_den(*den, tmp);
  if(b!=0) {
    if (buf[0] == '-') {
      mpz_set_si(*num, -1);
      mpz_set_si(*den, 1);
    } else {
      mpz_set_si(*num, 1);
      mpz_set_si(*den, 1);
    }
  }
  free(buf);
  mpq_clear(tmp);
}

static void inner_strterm_to_mpz(char *str, mpz_t *num, mpz_t *den){
  mpq_t tmp;
  mpq_init(tmp);
  mpz_set_ui(*num, 1);
  mpz_set_ui(*den, 1);
  size_t m = 0;
  size_t len = strlen(str);
  for(m = 0; m < len; m++){
    if(str[m]=='*' || str[m]==',') break;
  }
  const size_t sz = m;
  char *buf = (char *)malloc((sz+1) * sizeof(char));
  for(long i = 1; i < sz; i++){
    buf[i-1]=str[i];
  }
  buf[sz-1]='\0';

  int b = mpq_set_str(tmp, buf, 10);
  mpq_get_num(*num, tmp);
  mpq_get_den(*den, tmp);
  if(b!=0) {
    mpz_set_si(*num, 1);
    mpz_set_si(*den, 1);
  }
  if(str[0]=='-'){
    mpz_neg(*num, *num);
  }
  free(buf);
  mpq_clear(tmp);
}

static int get_coefficient_mpz_and_term_from_line(char *line, int32_t nterms,
                                          int32_t field_char,
                                          data_gens_ff_t *gens, int32_t pos){
  char *prev_pos = NULL;
  size_t term_size = 50000;
  char *term  = (char *)malloc(term_size * sizeof(char));

  prev_pos = line;
  get_term(line, &prev_pos, &term, &term_size);
  if(term != NULL){

    beginning_strterm_to_mpz(term, gens->mpz_cfs[pos], gens->mpz_cfs[pos+1]);
    store_exponent(term, gens, pos/2*gens->nvars);
    for(int j = 2; j < 2*nterms; j+=2){
      get_term(line, &prev_pos, &term, &term_size);
      inner_strterm_to_mpz(term, gens->mpz_cfs[pos+j], gens->mpz_cfs[pos+j+1]);
      store_exponent(term, gens, ((pos+j)/2)*gens->nvars);
    }
    free(term);
    return 0;
  }
  free(term);

  return 1;
}

static void get_coeffs_and_exponents_ff32(FILE *fh, nelts_t all_nterms,
        int32_t *nr_gens, data_gens_ff_t *gens){
    int32_t pos = 0;
    size_t len = 0;
    ssize_t nread;

    char *line  = NULL; 
    if(getline(&line, &len, fh) !=-1){
    }
    if(getline(&line, &len, fh) !=-1){
    }

    gens->cfs = (int32_t *)(malloc(sizeof(int32_t) * all_nterms));
    gens->exps = (int32_t *)calloc(all_nterms * gens->nvars, sizeof(int32_t));
    for (int32_t i = 0; i < *nr_gens; i++) {
        do {
            nread = getdelim(&line, &len, ',', fh);
        } while (nread != -1 && is_line_empty(line, ','));
        if (i < *nr_gens - 1) {
            remove_trailing_delim(line, &len, ',');
        }
        remove_newlines(line, &len);
        trim_polynomial(line);
        if(get_coefficient_ff_and_term_from_line(line, gens->lens[i], gens->field_char,
                    gens, pos)){
            fprintf(stderr, "Error when reading file (exit but things need to be free-ed)\n");
            free(line);
            fclose(fh);
            exit(1);
        }
        pos += gens->lens[i];
    }
    free(line);
}


static void get_coeffs_and_exponents_mpz(FILE *fh, nelts_t all_nterms,
        int32_t *nr_gens, data_gens_ff_t *gens){
    int32_t pos = 0;
    size_t len = 0;
    ssize_t nread;

    char *line  = NULL; 
    if(getline(&line, &len, fh) !=-1){
    }
    if(getline(&line, &len, fh) !=-1){
    }

    gens->cfs = (int32_t*)(malloc(sizeof(int32_t) * all_nterms));

    gens->mpz_cfs = (mpz_t **)(malloc(sizeof(mpz_t *) * 2 * all_nterms));
    for(long i = 0; i < 2 * all_nterms; i++){
      gens->mpz_cfs[i]  = (mpz_t *)malloc(sizeof(mpz_t));
      mpz_init(*(gens->mpz_cfs[i]));
    }

    gens->exps = (int32_t *)calloc(all_nterms * gens->nvars, sizeof(int32_t));
    for (int32_t i = 0; i < *nr_gens; i++) {
        do {
            nread = getdelim(&line, &len, ',', fh);
        } while (nread != -1 && is_line_empty(line, ','));
        if (i < *nr_gens - 1) {
            remove_trailing_delim(line, &len, ',');
        }
        remove_newlines(line, &len);
        trim_polynomial(line);
        if(get_coefficient_mpz_and_term_from_line(line, gens->lens[i], gens->field_char,
                    gens, pos)){
            fprintf(stderr, "Error when reading file (exit but things need to be free-ed)\n");
            free(line);
            fclose(fh);
            exit(1);
        }
        pos += 2 * gens->lens[i];
    }
    free(line);
}



static inline void initialize_data_gens(int32_t nvars, int32_t ngens, int32_t field_char, data_gens_ff_t *gens){
  gens->nvars = nvars;
  gens->ngens = ngens;
  gens->field_char = field_char;
  gens->change_var_order  = -1;
  gens->linear_form_base_coef = 0;
  gens->rand_linear = 0;
  gens->lens = (int32_t *)malloc(sizeof(int32_t) * ngens);
}

static int duplicate_vnames(char **vnames, int32_t nvars) {
  int32_t i, j;

  for (i = 1; i < nvars; ++i) {
    for (j = 0; j < i; ++j) {
      if (strcmp(vnames[i], vnames[j]) == 0) {
        fprintf(stderr, "Duplicate variable name %s in input file.\n", vnames[i]);
        return 1;
      }
    }
  }

  return 0;
}

//nr_gens is a pointer to the number of generators
static inline void get_data_from_file(char *fn, int32_t *nr_vars,
                                      int32_t *field_char,
                                      int32_t *nr_gens, data_gens_ff_t *gens){
  *nr_vars = get_nvars(fn);
  if (*nr_vars == -1)
    printf("Bad file format (first line).\n");

  *nr_gens = get_ngenerators(fn);
  if (*nr_gens == -1)
    printf("Bad file format (generators).\n");

  const int max_line_size  = 1073741824;
  char *line  = (char *)malloc((nelts_t)max_line_size * sizeof(char));

  FILE *fh  = fopen(fn,"r");

  /** allocate memory for storing variable names */
  char **vnames = (char **)malloc((*nr_vars) * sizeof(char *));
  get_variables(fh, line, max_line_size, nr_vars, gens, vnames);
  if (duplicate_vnames(vnames, *nr_vars) == 1) {
      for(int32_t i = 0; i < *nr_vars; i++){
          free(vnames[i]);
      }
      free(vnames);
      free(line);
      exit(1);
  }
  get_characteristic(fh, line, max_line_size, field_char, vnames);

  initialize_data_gens(*nr_vars, *nr_gens, *field_char, gens);

  nelts_t nterms, all_nterms = 0;
  get_nterms_and_all_nterms(fh, max_line_size, gens, nr_gens,
                            &nterms, &all_nterms);

  fclose(fh);
  fh = fopen(fn, "r");

  if(gens->field_char){
    get_coeffs_and_exponents_ff32(fh, all_nterms, nr_gens, gens);
  }
  else{
    get_coeffs_and_exponents_mpz(fh, all_nterms, nr_gens, gens);
  }

  free(line);
  fclose(fh);

  return;
}

static inline void display_gens_ff(FILE *fh, data_gens_ff_t *gens){
  int64_t pos = 0;
  int32_t c;
  for(long i = 0; i < gens->ngens; i++){
    for(long j = 0; j < gens->lens[i]-1; j++){
      c = gens->cfs[pos+j];
      if(c != 1){
        fprintf(fh, "%d", c);
        display_monomial(fh, gens, pos+j, &gens->exps);
      }
      else{
        display_monomial_single(fh, gens, pos+j, &gens->exps);
      }
      fprintf(fh, "+");
    }
    c = gens->cfs[pos+gens->lens[i]-1];
    if(c!= 1){
      fprintf(fh, "%d*", c);
      display_monomial(fh, gens, pos+gens->lens[i]-1, &gens->exps);
    }
    else{
      int b = display_monomial_single(fh, gens, pos+gens->lens[i]-1, &gens->exps);
      if(b==0){
        fprintf(fh, "1");
      }
    }
    pos+=gens->lens[i];
    if(i < gens->ngens-1){
      fprintf(fh,",\n");
    }
    else{
      fprintf(fh,"\n");
    }
  }
}

static void display_gens_mpz(FILE *fh, data_gens_ff_t *gens){
  long pos = 0, posex = 0;
  for(long i = 0; i < gens->ngens; i++){
    for(long j = 0; j < gens->lens[i]-1; j++){
      if(mpz_cmp_ui(*(gens->mpz_cfs[pos+2*j]), 1) != 0){
        mpz_out_str(fh, 10, *(gens->mpz_cfs[pos+2*j]));
        display_monomial(fh, gens, posex+j, &gens->exps);
      }
      else{
        display_monomial_single(fh, gens, posex+j, &gens->exps);
      }
      if(mpz_cmp_ui(*(gens->mpz_cfs[pos+2*(j+1)]), 0) > 0){
        fprintf(fh, "+");
      }
    }
    if(mpz_cmp_ui(*(gens->mpz_cfs[pos+ 2*(gens->lens[i] -1)]), 1) != 0){
      mpz_out_str(fh, 10, *(gens->mpz_cfs[pos+2*(gens->lens[i]-1)]));
      display_monomial(fh, gens, posex+(gens->lens[i]-1), &gens->exps);
    }
    else{
      int b = display_monomial_single(fh, gens, posex+gens->lens[i]-1, &gens->exps);
      if(b==0){
        fprintf(fh, "1");
      }
    }
    pos+=2*gens->lens[i];
    posex+=gens->lens[i];
    if(i < gens->ngens-1){
      fprintf(fh,",\n");
    }
    else{
      fprintf(fh,"\n");
    }
  }
}

static inline void display_gens(FILE *fh, data_gens_ff_t *gens){
  for(int i = 0; i < gens->nvars-1; i++){
    fprintf(fh, "%s, ", gens->vnames[i]);
  }
  fprintf(fh, "%s\n", gens->vnames[gens->nvars-1]);
  fprintf(fh, "%d\n", gens->field_char);
  if(gens->field_char > 0){
    display_gens_ff(fh, gens);
  }
  else{
    display_gens_mpz(fh, gens);
  }
}

static inline void get_poly_bin(FILE *file, mpz_upoly_t pol){
  if(!fscanf(file, "%d\n", &pol->alloc)){
    fprintf(stderr, "Issue when reading binary file (alloc = %d)\n", pol->alloc);
    exit(1);
  }

  pol->coeffs = malloc(sizeof(mpz_t) * pol->alloc);
  pol->length = pol->alloc;

  for(int32_t i = 0; i < pol->length; i++){
    mpz_init(pol->coeffs[i]);
    if(!mpz_inp_raw(pol->coeffs[i], file)){
      fprintf(stderr, "An error occurred when reading file (i=%d)\n", i);
      exit(1);
    }
  }
}

static inline void get_poly(FILE *file, mpz_upoly_t pol){
  if(!fscanf(file, "%d\n", &pol->alloc)){
    fprintf(stderr, "Issue when reading binary file (alloc = %d)\n", pol->alloc);
    exit(1);
  }

  pol->coeffs = malloc(sizeof(mpz_t) * pol->alloc);
  pol->length = pol->alloc;
  for(int32_t i = 0; i < pol->length; i++){
    mpz_init(pol->coeffs[i]);
    if(!mpz_inp_str(pol->coeffs[i], file, 10)){
      fprintf(stderr, "An error occurred when reading file (i=%d)\n", i);
      exit(1);
    }
  }
}


static inline void get_single_param_from_file_bin(FILE *file, mpz_param_t param){
  get_poly_bin(file, param->elim);

  get_poly_bin(file, param->denom);

  if(!fscanf(file, "%d\n", &param->nvars)){
    fprintf(stderr, "Issue when reading binary file (nvars)\n");
    exit(1);
  }

  param->nsols = param->elim->length - 1;
  param->dquot = param->elim->length - 1;

  param->coords = malloc(sizeof(mpz_upoly_t) * param->nvars);
  param->cfs = malloc(sizeof(mpz_t) * param->nvars);

  for(int32_t i = 0; i < param->nvars - 1; i++){
    get_poly_bin(file, param->coords[i]);

    mpz_init(param->cfs[i]);
    if(!mpz_inp_raw(param->cfs[i], file)){
      fprintf(stderr, "An error occurred when reading file (lcm coord i=%d)\n", i);
      exit(1);
    }

  }

}

static inline void get_single_param_from_file(FILE *file, mpz_param_t param){
  get_poly(file, param->elim);

  get_poly(file, param->denom);
  if(!fscanf(file, "%d\n", &param->nvars)){
    fprintf(stderr, "Issue when reading binary file (nvars)\n");
    exit(1);
  }

  param->nsols = param->elim->length - 1;
  param->dquot = param->elim->length - 1;

  param->coords = malloc(sizeof(mpz_upoly_t) * param->nvars);
  param->cfs = malloc(sizeof(mpz_t) * param->nvars);

  for(int32_t i = 0; i < param->nvars - 1; i++){
    get_poly(file, param->coords[i]);

    mpz_init(param->cfs[i]);

    if(!mpz_inp_str(param->cfs[i], file, 10)){
      fprintf(stderr, "An error occurred when reading file (i=%d)\n", i);
      exit(1);
    }

  }

}

static inline void get_params_from_file_bin(char *fn, mpz_param_array_t lparam){
  FILE *file = fopen(fn,"r");
  int32_t nb = 0;
  if(!fscanf(file, "%d\n", &nb)){
    fprintf(stderr, "Issue when reading binary file (nb = %d)\n", nb);
    exit(1);
  }
  lparam->nb = nb;

  lparam->params = malloc(sizeof(mpz_param_t) * lparam->nb);
  for(int32_t i = 0; i < lparam->nb; i++){
    get_single_param_from_file_bin(file, lparam->params[i]);
  }
  fclose(file);
}

static inline void get_params_from_file(char *fn, mpz_param_array_t lparam){
  FILE *file = fopen(fn,"r");
  int32_t nb = 0;
  if(!fscanf(file, "%d\n", &nb)){
    fprintf(stderr, "Issue when reading binary file (nb = %d)\n", nb);
    exit(1);
  }
  lparam->nb = nb;

  lparam->params = malloc(sizeof(mpz_param_t) * lparam->nb);
  for(int32_t i = 0; i < lparam->nb; i++){
    get_single_param_from_file(file, lparam->params[i]);
  }
  fclose(file);
}
