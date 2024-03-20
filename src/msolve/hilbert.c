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

#include "../fglm/data_fglm.c"
#include "../fglm/libfglm.h"
#include "../neogb/meta_data.h"
#define REDUCTION_ALLINONE 1

static void (*copy_poly_in_matrix_from_bs)(sp_matfglm_t* matrix,
                                           long nrows,
                                           bs_t *bs,
                                           ht_t *ht,
                                           long idx, long len,
                                           long start, long pos,
                                           int32_t *lmb,
                                           const int nv,
                                           const long fc);

static void (*copy_nf_in_matrix_from_bs)(sp_matfglm_t* matrix,
					 long nrows,
					 long pos,
					 int32_t *lmb,
					 const bs_t * const tbr,
					 const ht_t * const bht,
					 int32_t * evi,
					 const md_t *st,
					 const int nv);


static int is_pure_power(const int32_t *bexp, const int nv){
    int cnt = 0;
    for(int i = 0; i < nv; i++){
        if(bexp[i]==0){
            cnt++;
        }
    }
    if(cnt>=nv-1){
        return 1;
    }
    return 0;

}

static int32_t *get_lead_monomials(const int32_t *bld, int32_t **blen, int32_t **bexp,
				   data_gens_ff_t *gens){
  long ngens = bld[0];
  long nvars = gens->nvars;
  int32_t *bexp_lm = malloc(sizeof(int32_t)*ngens*nvars);

  long npos = 0;
  for(long i = 0; i < ngens; i++){
    for(long k = 0; k < nvars; k++){
      (bexp_lm)[i*nvars + k] = (*bexp)[npos*nvars + k];
    }
    npos+=(*blen)[i];
  }

  return bexp_lm;
}


static inline int has_dimension_zero(const long length, const int nv,
                                     int32_t *bexp_lm){
  long ppwr = 0;
  for(long i = 0; i < length; i++){
    int boo = is_pure_power(bexp_lm + i*nv, nv);
    if(boo){
      ppwr++;
    }
    if(ppwr >= nv){
      return 1;
    }
  }
  if(ppwr >= nv){
    return 1;
  }
  return 0;
}

static inline int is_divisible_exp(long nvars, int32_t *exp1, int32_t *exp2){
  for(long i = 0; i < nvars; i++){
    if(exp2[i]>exp1[i]){
      return 0;
    }
  }
  return 1;
}

static inline int is_divisible_lexp(long nvars, long length,
                                    int32_t *exp1, int32_t *bexp){
  for(long i = 0; i < length; i++){
    if(is_divisible_exp(nvars, exp1, (bexp+i*nvars))){
      return 1;
    }
  }
  return 0;
}

static inline int is_divisible_exp_elim(int32_t nvars, int32_t nev,
                                        int32_t *exp1, int32_t *exp2){
  for(int32_t i = nev; i < nvars; i++){
    if(exp2[i]>exp1[i]){
      return 0;
    }
  }
  return 1;
}

static inline int is_divisible_lexp_elim(int32_t nvars, int32_t nev,
                                         long length,
                                         int32_t *exp1, int32_t *bexp){
  for(long i = 0; i < length; i++){
    if(is_divisible_exp_elim(nvars, nev, exp1, (bexp+i*nvars))){
      return 1;
    }
  }
  return 0;
}

static inline int is_divisible_lexp_without_last_variable(long nvars, long length,
							  int32_t
							  *exp1,
							  int32_t *bexp){
  for(long i = 0; i < length; i++){
    if(is_divisible_exp(nvars-1, exp1, (bexp+i*nvars))){
      return 1;
    }
  }
  return 0;
}

/*ind is an array of long of length length
  adds entries of ind  */
static inline int32_t sum(int32_t *ind, int32_t length){
  int32_t s= 0;
  for(int32_t i = 0; i < length; i++){
    s += ind[i];
  }
  return s;
}

/*
  - len1 is the current dimension of the quotient
  - len_lm is the number of polynomials in the GB
  - basis will contain the final monomial basis
  - new_basis
*/
static inline int32_t generate_new_elts_basis(int32_t nvars, int32_t *ind,
                                              len_t len1,
                                              int32_t *basis, int32_t *new_basis,
                                              len_t len_lm,
                                              int32_t *bexp_lm){

  len_t c = 0;
  len_t def = 0;
  for(int64_t n = nvars-1; n >= 0; n--){
    for(len_t i = ind[nvars-1-n]; i < len1; i++){
      for(int64_t k = 0; k < nvars; k++){
        new_basis[c*nvars+k] = basis[i*nvars+k];
      }
      new_basis[c*nvars+n]++;
      if(!is_divisible_lexp(nvars, len_lm, (new_basis+c*nvars), bexp_lm)){
        c++;
      }
      else{
        def++;
      }
    }
  }

  return c;
}

/* takes as input an exponent and an index and returns the degree of the
   variable next to the index */
static inline int degree_prev_var(int32_t *exp, long n){
  return(exp[n+1]);
}

static inline void update_indices(int32_t *ind, int32_t *basis,
                                  long dquot, long new_length, long nvars){

  ind[0] = dquot;
  for(long n = nvars - 2; n >= 0 ; n--){
    for(long i = ind[nvars-1-n-1]; i < dquot + new_length; i++){
      if(degree_prev_var(basis+i*nvars, n)==0){
        ind[nvars-n-1] = i;
        break;
      }
      for(long k = n; k>=0; k--){
        ind[nvars-k-1] = dquot + new_length;
      }
    }
  }
}

/** length is the length of the GB
    nvars is the number of variables
    bexp_lm encodes the leading monomials
    dquot is a pointer to an integer that will be the dimension of the quotient

    assumes bexp_lm is for grevlex ordering
 */
static inline int32_t *monomial_basis(long length, long nvars,
                                      int32_t *bexp_lm, long *dquot){
  /* basis will contain the monomial basis */
  int32_t *basis = calloc(nvars, sizeof(int32_t));
  (*dquot) = 0;

  if(is_divisible_lexp(nvars, length, (basis), bexp_lm)){
    free(basis);
    return NULL;
  }
  else{
    (*dquot)++;
  }
  int32_t *ind = calloc(nvars, sizeof(int32_t));

#ifdef DEBUGHILBERT
  fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

  int32_t *new_basis = malloc(sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
  int32_t new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                               basis, new_basis,
                                               length, bexp_lm);
#ifdef DEBUGHILBERT
  //  display_monomials_from_array(stderr, new_length, new_basis, gens);
  fprintf(stderr, "%ld new elements.\n", new_length);
#endif
  while(new_length>0){
    int32_t *basis2 = realloc(basis,
                              ((*dquot) + new_length) * nvars * sizeof(int32_t *));
    if(basis2==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }

    basis = basis2;
    for(long i = 0; i < new_length; i++){
      for(long k = 0; k < nvars; k++){
        (basis)[((*dquot) + i)*nvars+k] = (new_basis)[i*nvars+k];
      }
    }

    update_indices(ind, basis, *dquot, new_length, nvars);
    (*dquot) += new_length;

#ifdef DEBUGHILBERT
    fprintf(stderr, "Update indices\n");
    for(long i = 0; i < nvars; i++) fprintf(stderr, "%ld ", ind[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

    int32_t *new_basis2 = realloc(new_basis,
                                 sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
    if(new_basis==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }
    new_basis=new_basis2;
    new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                         basis, new_basis,
                                         length, bexp_lm);
#ifdef DEBUGHILBERT
    fprintf(stderr, "%ld new elements.\n", new_length);
#endif
  }

  free(new_basis);
  free(ind);
  return basis;
}



/** nvars is the number of variables
    sht is a secondary htable of the monomials
    dquot is a integer representing the dimension of the
    subspace
*/
static inline int32_t *monomial_basis_colon(long length, long nvars,
					    int32_t *bexp_lm, long *dquot,
					    const long maxdeg){
  int32_t *basis = calloc(nvars, sizeof(int32_t));
  (*dquot) = 0;

  if(is_divisible_lexp(nvars, length, (basis), bexp_lm)){
    fprintf(stderr, "Stop\n");
    free(basis);
    return NULL;
  }
  else{
    (*dquot)++;
  }
  int32_t *ind = calloc(nvars, sizeof(int32_t));

#ifdef DEBUGHILBERT
  fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif


  int32_t *new_basis = malloc(sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
  long new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                            basis, new_basis,
                                            length, bexp_lm);
#ifdef DEBUGHILBERT
  display_monomials_from_array(stderr, new_length, new_basis, gens);
  fprintf(stderr, "%ld new elements.\n", new_length);
#endif
  long deg = 1;
  while(new_length>0 && deg <= maxdeg){
    int32_t *basis2 = realloc(basis, ((*dquot) + new_length) * nvars * sizeof(int32_t *));
    if(basis2==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }

    basis = basis2;
    for(long i = 0; i < new_length; i++){
      for(long k = 0; k < nvars; k++){
        (basis)[((*dquot) + i)*nvars+k] = (new_basis)[i*nvars+k];
      }
    }

    update_indices(ind, basis, *dquot, new_length, nvars);
    (*dquot) += new_length;

#ifdef DEBUGHILBERT
    fprintf(stderr, "Update indices\n");
    for(long i = 0; i < nvars; i++) fprintf(stderr, "%ld ", ind[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

    int32_t *new_basis2 = realloc(new_basis,
                                 sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
    if(new_basis==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }
    new_basis=new_basis2;
    new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                         basis, new_basis,
                                         length, bexp_lm);
#ifdef DEBUGHILBERT
    fprintf(stderr, "%ld new elements.\n", new_length);
#endif
    deg++;
  }

  free(new_basis);
  free(ind);
  return basis;
}


/** nvars is the number of variables
    sht is a secondary htable of the monomials
    dquot is a integer representing the dimension of the
    subspace
*/
static inline int32_t *monomial_basis_colon_no_zero(long length, long nvars,
						    int32_t *bexp_lm, long *dquot,
						    const long maxdeg){
  int32_t *basis = calloc(nvars, sizeof(int32_t));
  (*dquot) = 0;

  if(is_divisible_lexp(nvars, length, (basis), bexp_lm)){
    fprintf(stderr, "Stop\n");
    free(basis);
    return NULL;
  }
  else{
    (*dquot)++;
  }
  int32_t *ind = calloc(nvars, sizeof(int32_t));

#ifdef DEBUGHILBERT
  fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

  int32_t *new_basis = malloc(sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
  long new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                            basis, new_basis,
                                            length, bexp_lm);
#ifdef DEBUGHILBERT
  display_monomials_from_array(stderr, new_length, new_basis, gens);
  fprintf(stderr, "%ld new elements.\n", new_length);
#endif
  long deg = 1;
  while(new_length>0 && deg <= maxdeg){
    int32_t *basis2 = realloc(basis, ((*dquot) + new_length) * nvars * sizeof(int32_t *));
    if(basis2==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }

    basis = basis2;
    for(long i = 0; i < new_length; i++){
      for(long k = 0; k < nvars; k++){
        (basis)[((*dquot) + i)*nvars+k] = (new_basis)[i*nvars+k];
      }
    }

    update_indices(ind, basis, *dquot, new_length, nvars);
    (*dquot) += new_length;

#ifdef DEBUGHILBERT
    fprintf(stderr, "Update indices\n");
    for(long i = 0; i < nvars; i++) fprintf(stderr, "%ld ", ind[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

    int32_t *new_basis2 = realloc(new_basis,
                                 sizeof(int32_t) * nvars * (sum(ind, nvars) + nvars));
    if(new_basis==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }
    new_basis=new_basis2;
    new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                         basis, new_basis,
                                         length, bexp_lm);
#ifdef DEBUGHILBERT
    fprintf(stderr, "%ld new elements.\n", new_length);
#endif
    deg++;
  }
  free(new_basis);
  free(ind);

  /* cleanup by removing monomials that will be sent to 0 after
     iterative multiplication by xn */
  long new_dquot=0;
  for (long i = 0; i < *dquot; i++) {
    if (is_divisible_lexp_without_last_variable(nvars,length,
						(basis+i*nvars),bexp_lm)) {
      /* keep this monomial */
      for (long k = 0; k < nvars; k++) {
	basis[new_dquot*nvars+k]= basis[i*nvars+k];
      }
      new_dquot++;
    }
  }
  if (new_dquot == 0) {
    fprintf(stderr, "Vector space is too small\n");
    exit(1);
  }
  *dquot= new_dquot;
  basis= realloc (basis, new_dquot * nvars * sizeof(int32_t));
  return basis;
}


static inline long get_div_xn(int32_t *bexp_lm, long length, long nvars,
                             int32_t *div_xn){
  long l = 0;
  for(long i = 0; i < length; i++){
    if(bexp_lm[i*nvars + (nvars-1)]!=0){
      div_xn[l] = i;
      l++;
    }
  }
  return l;
}

static inline long get_div_xn_bounded(int32_t *bexp_lm, long length, long nvars,
				      int32_t *div_xn, int32_t *div_not_xn,
				      long* l_not,
				      const long maxdeg){
  long l = 0;
  *l_not = 0;
  for(long i = 0; i < length; i++){
    long deg = 0;
    for (long j = 0; j < nvars; j++) {
      deg += bexp_lm[i*nvars+j];
    }
    if (deg <= maxdeg) {
      if (bexp_lm[i*nvars + (nvars-1)]!=0){
	div_xn[l] = i;
	l++;
      }else{
	div_not_xn[*l_not] = i;
	(*l_not)++;
      }
    }
  }
  return l;
}

/*
One of the exp. encoded is hashed
 */
static inline int is_equal_exponent_dm(bs_t *bs, ht_t *ht, long idx, long pos,
                                       int32_t *exp2, const long nv){
  const bl_t bi = bs->lmps[idx];
  hm_t *dt = bs->hm[bi] + OFFSET;

  for(long k = 0; k < nv-1 ; k++){
    /* in hash table exponent vectors store their degree in the first entry,
     * thus we have to go by "k+1" */
    int e = (int32_t)ht->ev[dt[pos]][k+1];
    if(e != exp2[k]){
      return 0;
    }
  }
  int e = (int32_t)ht->ev[dt[pos]][nv];
  return (e == exp2[nv - 1]);
}

static inline int is_equal_exponent(int32_t *exp1, int32_t *exp2, const long nvars){
  for(long i = 0; i < nvars - 1; i++){
    if(exp1[i]!=exp2[i]){
      return 0;
    }
  }
  return ((exp1[nvars-1]) == exp2[nvars-1]);
}

static inline int is_equal_exponent_elim(int32_t *exp1, int32_t *exp2,
                                         const long nvars, const long elim){
  for(long i = 0; i < nvars - elim - 1; i++){
    if(exp1[i + elim]!=exp2[i]){
      return 0;
    }
  }
  return ((exp1[nvars-1]) == exp2[nvars - elim -1]);
}

static inline int is_equal_exponent_bs(const ht_t * const exp1, int32_t hmj,
				       int32_t *evi,
				       int32_t *exp2,
				       const long nvars){
  /* printf ("\n"); */
  for(long i = 0; i < nvars - 1; i++){
    /* printf ("exp1[%ld]=%d\texp2[%ld]=%d\n",i,exp1->ev[hmj][evi[i]],i,exp2[i]); */
    if(exp1->ev[hmj][evi[i]]!=exp2[i]){
      return 0;
    }
  }
  /* printf ("exp1[%ld]=%d\texp2[%ld]=%d\n",nvars-1,exp1->ev[hmj][evi[nvars-1]], */
  /* 	  nvars-1,exp2[nvars-1]); */
  return ((exp1->ev[hmj][evi[nvars-1]]) == exp2[nvars-1]);
}

/* returns -1 if exp1 is smaller
 * return 0 if they are the same
 * returns 1 if exp1 is larger
 */
static inline int is_larger_exponent(int32_t *exp1, int32_t *exp2, const long nvars){
  int32_t deg1 = 0;
  int32_t deg2 = 0;
  for(long i = 0; i < nvars; i++){
    deg1 += exp1[i];
    deg2 += exp2[i];
  }
  /* printf ("deg1= %d, deg2=%d\n",deg1,deg2); */
  if (deg1 < deg2){
    return -1;
  }
  if (deg1 > deg2){
    return 1;
  }
  for (long i = nvars-1; i>1; i--) {
    /* printf ("x[%ld]: %d, %d\n",i,exp1[i],exp2[i]); */
    if(exp1[i]<exp2[i]){
      return 1;
    }
    if(exp1[i]>exp2[i]){
      return -1;
    }
  }
  return 0;
}

/* returns -1 if exp1 is smaller
 * return 0 if they are the same
 * returns 1 if exp1 is larger
 */
static inline int is_larger_exponent_bs(const ht_t * const exp1, int32_t hmj,
					int32_t *evi,
					int32_t *exp2,
					const long nvars){
  int32_t deg1 = 0;
  int32_t deg2 = 0;
  for(long i = 0; i < nvars; i++){
    deg1 += exp1->ev[hmj][evi[i]];
    deg2 += exp2[i];
  }
  /* printf ("deg1= %d, deg2=%d\n",deg1,deg2); */
  if (deg1 < deg2){
    return -1;
  }
  if (deg1 > deg2){
    return 1;
  }
  for (long i = nvars-1; i>1; i--) {
    /* printf ("x[%ld]: %d, %d\n",i,exp1->ev[hmj][evi[i]],exp2[i]); */
    if(exp1->ev[hmj][evi[i]]<exp2[i]){
      return 1;
    }
    if(exp1->ev[hmj][evi[i]]>exp2[i]){
      return -1;
    }
  }
  return 0;
}

static inline int is_equal_exponent_xxn(int32_t *exp1, int32_t *exp2, const long nvars){
  for(long i = 0; i < nvars - 1; i++){
    if(exp1[i]!=exp2[i]){
      return 0;
    }
  }
  return ((exp1[nvars-1]+1) == exp2[nvars-1]);
}

static inline int is_divisible_exponent_xxn(int32_t *exp1, int32_t *exp2,
					    const long nvars){
  for(long i = 0; i < nvars - 1; i++){
    if(exp1[i]<exp2[i]){
      return 0;
    }
  }
  return ((exp1[nvars-1]+1) >= exp2[nvars-1]);
}

static inline int member_xxn(int32_t *exp, int32_t *bexp_lm, int len,
                             long *pos, const int nv){
  for(long i=1; i < len; i++){
    if(is_equal_exponent_xxn(exp, bexp_lm + (i*nv), nv)){
      *pos = i;
      return 1;
    }
  }
  return 0;
}

static inline void copy_poly_in_matrix_old(data_gens_ff_t *gens,
                                       sp_matfglm_t* matrix,
                                       long nrows,
                                       int32_t *bcf, int32_t **bexp,
                                       int32_t **blen, long start, long pos,
                                       int32_t *lmb){
  int32_t j;
  int32_t end = start + pos;//(*blen)[pos];
  fprintf(stderr, "copy_poly\n");
  /* fprintf(stderr, "\nstart = %ld, end = %ld\n", start, end); */
  /* for(j = start; j < end; j++){ */
  /*   display_term(stderr, j, gens, blen, bcf, bexp); */
  /*   if(j < end - 1){ */
  /*   //    if(j < (*blen)[pos] - 1){ */
  /*     fprintf(stderr, "+"); */
  /*   } */
  /* } */

  long N = nrows * (matrix->ncols) - (start + 1);
  if((end-start) == matrix->ncols + 1){
    for(j = start + 1; j < end; j++){
      matrix->dense_mat[N + j] = gens->field_char - bcf[(end + start) - j];
    }
  }
  else{
    long i;
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    display_monomial_full(stderr, gens->nvars, NULL, 0, (*bexp) + (start+1)*gens->nvars);
    fprintf(stderr, ", ");
    display_monomial_full(stderr, gens->nvars, NULL, 0, lmb + ((end-start-2))*gens->nvars);
    fprintf(stderr, ", ");
    display_monomial_full(stderr, gens->nvars, NULL, 0, lmb + ((end-start-1))*gens->nvars);
    fprintf(stderr, ", ");
    display_monomial_full(stderr, gens->nvars, NULL, 0, lmb + ((end-start-3))*gens->nvars);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    for(i=(matrix->ncols)-1; i >= 0; i--){
      display_monomial_full(stderr, gens->nvars, NULL, 0, (*bexp) + (start+1)*gens->nvars);fprintf(stderr, ", ");
      display_monomial_full(stderr, gens->nvars, NULL, 0, lmb + (i*gens->nvars));
      fprintf(stderr, "\n");
      if(is_equal_exponent((*bexp) + (start+1)*gens->nvars,
                           lmb+(i*gens->nvars),
                           gens->nvars))
        break;
    }

    if(i+1 == end - start - 1){
      for(j = start + 1; j < end; j++){
        matrix->dense_mat[N + j] = gens->field_char - bcf[(end + start) - j];
      }
    }
    else{
      fprintf(stderr, "Some coefficients are null\n");
      fprintf(stderr, "i + 1 = %ld\n", i + 1);
      fprintf(stderr, "end - start - 1 = %ld\n", end - start - 1);
      exit(1);
    }
  }
}

/**

copy of the pol whose coeffs are indexed between start and start + pos

 **/


static inline void copy_poly_in_matrix(sp_matfglm_t* matrix,
                                       long nrows,
                                       int32_t *bcf, int32_t **bexp,
                                       int32_t **blen, long start, long pos,
                                       int32_t *lmb,
                                       const int nv,
                                       const long fc){
  int32_t j;
  long end = start + pos;//(*blen)[pos];

#if DEBUGBUILDMATRIX > 1
  fprintf(stderr, "\nstart = %ld, end = %ld\n", start, end);
  for(j = start; j < end; j++){
    //    display_term(stderr, j, gens, blen, bcf, bexp);
    if(j < end - 1){
    //    if(j < (*blen)[pos] - 1){
      fprintf(stderr, "+");
    }
  }
  fprintf(stderr, "\n");
#endif

  long N = nrows * (matrix->ncols) - (start + 1);

  if((end-start) == matrix->ncols + 1){
    for(j = start + 1; j < end; j++){
      matrix->dense_mat[N + j] = fc - bcf[(end + start) - j];
    }
  }
  else{
    if(is_equal_exponent((*bexp) + (start+1)*nv,
                         lmb+((end-start-2)*nv),
                         nv)){
      for(j = start + 1; j < end; j++){
        matrix->dense_mat[N + j] = fc - bcf[(end + start) - j];
      }
    }
    else{
      long i;

      long N = nrows * matrix->ncols ;
      long k = 0;
      for(i = 0; i < matrix->ncols; i++){

          if(is_equal_exponent((*bexp) + (end - 1 - k) * nv,
                      lmb + i * nv,
                      nv)){
            matrix->dense_mat[N + i] = fc - bcf[end - 1 -  k];
            k++;
        }
      }
    }
  }
}

static inline void
copy_poly_in_matrixcol(sp_matfglmcol_t* matrix, long nrows,
		       int32_t *bcf, int32_t **bexp, int32_t **blen,
		       long start, long pos, int32_t *lmb,
		       const int nv, const long fc){
  int32_t j;
  long end = start + pos;//(*blen)[pos];

#if DEBUGBUILDMATRIX > 0
  fprintf(stderr, "\nstart = %ld, end = %ld\n", start, end);
  for(j = start; j < end; j++){
    //    display_term(stderr, j, gens, blen, bcf, bexp);
    if(j < end - 1){
    //    if(j < (*blen)[pos] - 1){
      fprintf(stderr, "+");
    }
  }
  fprintf(stderr, "\n");
#endif
  long N = nrows * (matrix->ncols) - (start + 1);

  if((end-start) == matrix->ncols + 1){
    for(j = start + 1; j < end; j++){
      matrix->dense_mat[N + j] = fc - bcf[(end + start) - j];
    }
  }
  else{
    if(is_equal_exponent((*bexp) + (start+1)*nv,
                         lmb+((end-start-2)*nv),
                         nv)){
      for(j = start + 1; j < end; j++){
        matrix->dense_mat[N + j] = fc - bcf[(end + start) - j];
      }
    }
    else{
      long i;

      long N = nrows * matrix->ncols ;
      long k = 0;
      for(i = 0; i < matrix->ncols; i++){
	if(is_equal_exponent((*bexp) + (end - 1 - k) * nv,
			     lmb + i * nv, nv)){
	  matrix->dense_mat[N + i] = fc - bcf[end - 1 -  k];
	  k++;
        }
      }
    }
  }
}

static inline void
copy_poly_in_matrixcol_no_zero(sp_matfglmcol_t* matrix, long nrows,
			       int32_t *bcf, int32_t **bexp, int32_t **blen,
			       long start, long pos, int32_t *lmb,
			       const int nv, const long fc){
  long end = start + pos;//(*blen)[pos];

#if DEBUGBUILDMATRIX > 0
  int32_t j;
  fprintf(stderr, "\nstart = %ld, end = %ld\n", start, end);
  for(j = start; j < end; j++){
    //    display_term(stderr, j, gens, blen, bcf, bexp);
    if(j < end - 1){
    //    if(j < (*blen)[pos] - 1){
      fprintf(stderr, "+");
    }
  }
  fprintf(stderr, "\n");
#endif
  long i;
  long N = nrows * matrix->ncols ;
  long k = 0;
  /* printf("["); */
  for(i = 0; i < matrix->ncols; i++){
    if (k < end) {
      int b = is_larger_exponent((*bexp) + (end - 1 - k) * nv, lmb + i * nv, nv);
      while (b < 0 && k < end-2) {
	k++;
	b = is_larger_exponent((*bexp) + (end - 1 - k) * nv, lmb + i *  nv, nv);
      }
      if(!b){ /* equal exponent */
	matrix->dense_mat[N + i] = fc - bcf[end - 1 -  k];
	/* printf("%u, ",matrix->dense_mat[N + i]); */
	k++;
      }
    }
  }
  /* printf("]\n"); */
}

static inline void
copy_extrapoly_in_vector(uint32_t* vector, long ncols, int32_t *lmb,
			 len_t pos, const bs_t * const tbr,
			 const ht_t * const bht, int32_t* evi,
			 const md_t *st, const int nv, const long maxdeg){
  len_t idx = tbr->lmps[pos];
  /* printf ("idx=%d\n",idx); */
  len_t * hm  = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  /* printf ("len=%d\n",len); */
  long i;
  long k = 0;
  /* to remove monomials outside the vector space, we just look at the
   * leading ones */
  uint32_t deglm = 0;
  for (long j = 0; j < nv; j++) {
    deglm += lmb[k*nv + j];
  }
  while (deglm > maxdeg) {
    k++;
    deglm = 0;
    for (long j = 0; j < nv; j++) {
      deglm += lmb[k*nv + j];
    }
  }
  /* printf ("starts at k=%ld with coeff %d\n",k,tbr->cf_32[tbr->hm[idx][COEFFS]][k]); */
  /* printf ("ends at  k=%ld with coeff %d\n",len-1-k, */
  /* 	  tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k]); */
  /* printf ("["); */
  for(i = 0; i < ncols; i++){
    if(is_equal_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv)){
      vector[i] = tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k];
      /* printf ("%u, ",vector[i]); */
      k++;
    }
  }
  /* printf("]\n"); */
}

static inline void
copy_extrapoly_in_vector_no_zero(uint32_t* vector, long ncols,
				 int32_t *lmb, len_t pos, const bs_t * const tbr,
				 const ht_t * const bht, int32_t* evi,
				 const md_t *st, const int nv,
				 const long maxdeg){

  len_t idx = tbr->lmps[pos];
  /* printf ("idx=%d\n",idx); */
  len_t * hm  = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  /* printf ("len=%d\n",len); */
  long i;
  long k = 0;
  /* printf ("["); */
  for(i = 0; i < ncols; i++){
    if (k < len) {
      int b = is_larger_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv);
      while (b < 0 && k < len-2) {
	k++;
	b = is_larger_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv);
      }
      if(!b){ /* equal_exponent */
	vector[i] = tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k];
	/* printf ("%u, ",vector[i]); */
	k++;
      }
    }
  }
  /* printf("]\n"); */
}

static inline void
copy_extrapoly_in_matrixcol(sp_matfglmcol_t* matrix, long nrows,
			    int32_t *lmb, len_t pos, const bs_t * const tbr,
			    const ht_t * const bht, int32_t* evi,
			    const md_t *st, const int nv,
			    const long maxdeg){
  len_t idx = tbr->lmps[pos];
  /* printf ("idx=%d\n",idx); */
  len_t * hm  = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  /* printf ("len=%d\n",len); */
  long i;
  long N = nrows * matrix->ncols ;
  long k = 0;
  /* to remove monomials outside the vector space, we just look at the
   * leading ones */
  uint32_t deglm = 0;
  for (long j = 0; j < nv; j++) {
    deglm += lmb[k*nv + j];
  }
  while (deglm > maxdeg) {
    k++;
    deglm = 0;
    for (long j = 0; j < nv; j++) {
      deglm += lmb[k*nv + j];
    }
  }
  /* printf ("starts at k=%ld with coeff %d\n",k,tbr->cf_32[tbr->hm[idx][COEFFS]][k]); */
  /* printf ("ends at  k=%ld with coeff %d\n",len-1-k, */
  /* 	  tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k]); */
  /* printf ("["); */
  for(i = 0; i < matrix->ncols; i++){
    if(is_equal_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv)){
      matrix->dense_mat[N + i] = tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k];
      /* printf ("%u, ",matrix->dense_mat[N+i]); */
      k++;
    }
  }
  /* printf("]\n"); */
}

static inline void
copy_extrapoly_in_matrixcol_no_zero(sp_matfglmcol_t* matrix, long nrows,
				    int32_t *lmb, len_t pos, const bs_t * const tbr,
				    const ht_t * const bht, int32_t* evi,
				    const md_t *st, const int nv,
				    const long maxdeg){

  len_t idx = tbr->lmps[pos];
  len_t * hm  = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  /* printf ("len=%d\n",len); */
  long i;
  long N = nrows * matrix->ncols ;
  long k = 0;
  /* printf("["); */
  for(i = 0; i < matrix->ncols; i++){
    if (k < len) {
      int b = is_larger_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv);
      while (b < 0 && k < len-2) {
	k++;
	b = is_larger_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv);
      }
      if(!b){ /* equal_exponent */
	matrix->dense_mat[N + i] = tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k];
	/* printf("%u, ",matrix->dense_mat[N + i]); */
	k++;
      }
    }
  }
  /* printf("]\n"); */
}


/**

idx is the position in the GB
len is the length of the pol

 **/

static inline void copy_poly_in_matrix_from_bs_8(sp_matfglm_t* matrix,
                                               long nrows,
                                               bs_t *bs,
                                               ht_t *ht,
                                               long idx, long len,
                                               long start, long pos,
                                               int32_t *lmb,
                                               const int nv,
                                               const long fc){
  int32_t j;
  long end = start + pos;

  long N = nrows * (matrix->ncols) - (start + 1);
  if((len) == matrix->ncols + 1){
    const bl_t bi = bs->lmps[idx];
    long k = 0;
    for(j = start + 1; j < end; j++){
      long ctmp  = bs->cf_8[bs->hm[bi][COEFFS]][len - k - 1];
      k++;
      matrix->dense_mat[N + j] = fc - ctmp;
    }
  }
  else{
    if(1==0 && is_equal_exponent_dm(bs, ht, idx, 1,
                            lmb+((end-start-2)*nv),
                            nv)){
      const bl_t bi = bs->lmps[idx];
      long k = 0;
      for(j = start + 1; j < end; j++){
        long ctmp  = bs->cf_8[bs->hm[bi][COEFFS]][len - k];
        k++;
        matrix->dense_mat[N + j] = fc - ctmp; //bcf[(end + start) - j];
      }
    }
    else{
      long i;
      long N = nrows * matrix->ncols ;
      long k = 0;

      const bl_t bi = bs->lmps[idx];

      for(i = 0; i < matrix->ncols; i++){
        int boo = is_equal_exponent_dm(bs, ht, idx, len - k - 1, //pos-1-k,
                                       lmb + i * nv,
                                       nv);
        if(boo){
            long ctmp  = bs->cf_8[bs->hm[bi][COEFFS]][len - k - 1];
            matrix->dense_mat[N + i] = fc - ctmp; //fc - bcf[end - 1 -  k];
            k++;
        }
      }
    }
  }
}


/**

idx is the position in the GB
len is the length of the pol

 **/

static inline void copy_poly_in_matrix_from_bs_16(sp_matfglm_t* matrix,
                                               long nrows,
                                               bs_t *bs,
                                               ht_t *ht,
                                               long idx, long len,
                                               long start, long pos,
                                               int32_t *lmb,
                                               const int nv,
                                               const long fc){
  int32_t j;
  long end = start + pos;

  long N = nrows * (matrix->ncols) - (start + 1);
  if((len) == matrix->ncols + 1){
    const bl_t bi = bs->lmps[idx];
    long k = 0;
    for(j = start + 1; j < end; j++){
      long ctmp  = bs->cf_16[bs->hm[bi][COEFFS]][len - k - 1];
      k++;
      matrix->dense_mat[N + j] = fc - ctmp;
    }
  }
  else{
    if(1==0 && is_equal_exponent_dm(bs, ht, idx, 1,
                            lmb+((end-start-2)*nv),
                            nv)){
      const bl_t bi = bs->lmps[idx];
      long k = 0;
      for(j = start + 1; j < end; j++){
        long ctmp  = bs->cf_16[bs->hm[bi][COEFFS]][len - k];
        k++;
        matrix->dense_mat[N + j] = fc - ctmp; //bcf[(end + start) - j];
      }
    }
    else{
      long i;
      long N = nrows * matrix->ncols ;
      long k = 0;

      const bl_t bi = bs->lmps[idx];

      for(i = 0; i < matrix->ncols; i++){
        int boo = is_equal_exponent_dm(bs, ht, idx, len - k - 1, //pos-1-k,
                                       lmb + i * nv,
                                       nv);
        if(boo){
            long ctmp  = bs->cf_16[bs->hm[bi][COEFFS]][len - k - 1];
            matrix->dense_mat[N + i] = fc - ctmp; //fc - bcf[end - 1 -  k];
            k++;
        }
      }
    }
  }
}


/**

idx is the position in the GB
len is the length of the pol

 **/

static inline void copy_poly_in_matrix_from_bs_32(sp_matfglm_t* matrix,
                                               long nrows,
                                               bs_t *bs,
                                               ht_t *ht,
                                               long idx, long len,
                                               long start, long pos,
                                               int32_t *lmb,
                                               const int nv,
                                               const long fc){
  int32_t j;
  long end = start + pos;

  long N = nrows * (matrix->ncols) - (start + 1);
  if((len) == matrix->ncols + 1){
    const bl_t bi = bs->lmps[idx];
    long k = 0;
    for(j = start + 1; j < end; j++){
      long ctmp  = bs->cf_32[bs->hm[bi][COEFFS]][len - k - 1];
      k++;
      matrix->dense_mat[N + j] = fc - ctmp;
    }
  }
  else{
    if(1==0 && is_equal_exponent_dm(bs, ht, idx, 1,
                            lmb+((end-start-2)*nv),
                            nv)){
      const bl_t bi = bs->lmps[idx];
      long k = 0;
      for(j = start + 1; j < end; j++){
        long ctmp  = bs->cf_32[bs->hm[bi][COEFFS]][len - k];
        k++;
        matrix->dense_mat[N + j] = fc - ctmp; //bcf[(end + start) - j];
      }
    }
    else{
      long i;
      long N = nrows * matrix->ncols ;
      long k = 0;

      const bl_t bi = bs->lmps[idx];

      for(i = 0; i < matrix->ncols; i++){
        int boo = is_equal_exponent_dm(bs, ht, idx, len - k - 1, //pos-1-k,
                                       lmb + i * nv,
                                       nv);
        if(boo){
            long ctmp  = bs->cf_32[bs->hm[bi][COEFFS]][len - k - 1];
            matrix->dense_mat[N + i] = fc - ctmp; //fc - bcf[end - 1 -  k];
            k++;
        }
      }
    }
  }
}

static inline void copy_nf_in_matrix_from_bs_8(sp_matfglm_t* matrix,
                                               long nrows,
                                               long pos,
                                               int32_t *lmb,
					       const bs_t * const tbr,
					       const ht_t * const bht,
					       int32_t * evi,
					       const md_t *st,
                                               const int nv){
  len_t idx = tbr->lmps[pos];
  len_t * hm = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  long N = nrows * matrix->ncols ;
  long i = 0;
  long k = 0;
  while(k < len) {
    if(is_equal_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv)){
      matrix->dense_mat[N + i] = tbr->cf_8[tbr->hm[idx][COEFFS]][len-1-k];
      k++;
    }
    i++;
  }
}

static inline void copy_nf_in_matrix_from_bs_16(sp_matfglm_t* matrix,
						long nrows,
						long pos,
						int32_t *lmb,
						const bs_t * const tbr,
						const ht_t * const bht,
						int32_t * evi,
						const md_t *st,
						const int nv){
  len_t idx = tbr->lmps[pos];
  len_t * hm = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  long N = nrows * matrix->ncols ;
  long i = 0;
  long k = 0;
  while(k < len) {
    if(is_equal_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv)){
      matrix->dense_mat[N + i] = tbr->cf_16[tbr->hm[idx][COEFFS]][len-1-k];
      k++;
    }
    i++;
  }
}

static inline void copy_nf_in_matrix_from_bs_32(sp_matfglm_t* matrix,
						long nrows,
						long pos,
						int32_t *lmb,
						const bs_t * const tbr,
						const ht_t * const bht,
						int32_t * evi,
						const md_t *st,
						const int nv){
  len_t idx = tbr->lmps[pos];
  len_t * hm = tbr->hm[idx]+OFFSET;
  len_t len = tbr->hm[idx][LENGTH];
  long N = nrows * matrix->ncols ;
  long i = 0;
  long k = 0;
  while(k < len) {
    if(is_equal_exponent_bs(bht,hm[len-1-k],evi,lmb + i * nv,nv)){
      matrix->dense_mat[N + i] = tbr->cf_32[tbr->hm[idx][COEFFS]][len-1-k];
      k++;
    }
    i++;
  }
}


/** length is the length of the GB
    nvars is the number of variables
    bexp_lm encodes the leading monomials
    dquot is a pointer to an integer that will be the dimension of the quotient

    assumes bexp_lm is for grevlex ordering

    compared to monomial_basis this function does not assume
    bexp_lm to be 0-dim.

    It returns all monomials that can appear in the GB
 */
static inline int32_t *monomial_basis_enlarged(long length, long nvars,
                                               int32_t *bexp_lm, long *dquot){

  int32_t maxdeg = sum(bexp_lm+(length-1)*(nvars), nvars);

  /* basis will contain the monomial basis */
  int32_t *basis = calloc(nvars, sizeof(int32_t));
  (*dquot) = 0;

  int32_t deg = 0;

  if(is_divisible_lexp(nvars, length, basis, bexp_lm)){
    free(basis);
    return NULL;
  }
  else{
    (*dquot)++;
  }

  int32_t *ind = calloc(nvars, sizeof(int32_t));

#ifdef DEBUGHILBERT
  fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif
  uint64_t len_newbs = nvars * (sum(ind, nvars) + nvars);

  int32_t *new_basis = malloc(sizeof(int32_t) * len_newbs);
  /* generates monomial basis candidates of degree 1 */
  int32_t new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                               basis, new_basis, length, bexp_lm);

  deg++;
#ifdef DEBUGHILBERT
  //  display_monomials_from_array(stderr, new_length, new_basis, gens);
  fprintf(stderr, "%ld new elements.\n", new_length);
#endif

  while(new_length>0 && deg <= maxdeg){
    uint64_t len_bs = ((*dquot) + new_length) * (nvars);
    int32_t *basis2 = realloc(basis,
                              ((*dquot) + new_length) * (nvars) * sizeof(int32_t));
    if(basis2==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }

    basis = basis2;

    for(long i = 0; i < new_length; i++){
      for(long k = 0; k < nvars; k++){
        (basis)[((*dquot) + i)*(nvars)+k] = (new_basis)[i*(nvars)+k];
      }
    }

    update_indices(ind, basis, *dquot, new_length, nvars);


    (*dquot) += new_length;

#ifdef DEBUGHILBERT
    fprintf(stderr, "Update indices\n");
    for(long i = 0; i < nvars; i++) fprintf(stderr, "%ld ", ind[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "new = %ld \n", sum(ind, nvars) + nvars);
#endif

    len_newbs = nvars * (sum(ind, nvars) + nvars);

    int32_t *new_basis2 = realloc(new_basis,
                                  sizeof(int32_t) * len_newbs);
    if(new_basis==NULL){
      fprintf(stderr, "Issue with realloc\n");
      exit(1);
    }
    new_basis=new_basis2;

    new_length = generate_new_elts_basis(nvars, ind, (*dquot),
                                         basis, new_basis,
                                         length, bexp_lm);
    deg++;

#ifdef DEBUGHILBERT
    fprintf(stderr, "%ld new elements.\n", new_length);
#endif
  }

  free(new_basis);
  free(ind);
  return basis;
}

/**

   lmb is the monomial basis (of the quotient ring) given by ascending order.

   dquo is the dimension of the quotient.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

 **/
static inline sp_matfglm_t * build_matrixn(int32_t *lmb, long dquot, int32_t bld,
                                           int32_t **blen, int32_t **bexp,
                                           int32_t *bcf,
                                           int32_t *bexp_lm,
                                           const int nv, const long fc){


  /* takes monomials in bexp_lm which are reducible by xn */
  /*div_xn contains the indices of those monomials*/
  int32_t *div_xn = calloc(bld, sizeof(int32_t));

  /* len_xn is the number of those monomials. */
  long len_xn = get_div_xn(bexp_lm, bld, nv, div_xn);

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the staircase) which are divisible by x_n = %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
#endif

  /* lengths of the polys which we need to build the matrix */

  int32_t *len_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  long pos = 0, k = 0;
  for(long i = 0; i < bld; i++){
    if(i==div_xn[k]){
      len_gb_xn[k]=(*blen)[i];
      start_cf_gb_xn[k]=pos;
      pos+=(*blen)[i];
      k++;
    }
    else{
      pos+=(*blen)[i];
    }
  }

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "Length of polynomials whose leading terms are divisible by x_n\n");
  for(long i = 0; i < len_xn; i++){
    fprintf(stderr, "%d, ", len_gb_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  sp_matfglm_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nrows = len_xn;
  long len1 = dquot * len_xn;
  long len2 = dquot - len_xn;

  if(posix_memalign((void **)&(matrix->dense_mat), 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  for(long i = 0; i < dquot; i++){

    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif
      /* mult by xn stays in the basis */
      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{
      /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
        copy_poly_in_matrix(matrix, nrows, bcf, bexp, blen,
                            start_cf_gb_xn[count], len_gb_xn[count], lmb,
                            nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
          fprintf(stderr, "One should not arrive here (build_matrix)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix);

          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
          return NULL;
        }
      }
      else{
        fprintf(stderr, "\nStaircase is not generic\n");
        fprintf(stderr, "Multiplication by ");
        display_monomial_full(stderr, nv, NULL, 0, exp);
        fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);
        matrix  = NULL;

        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return matrix;
      }
    }
  }

  /* assumes the entries of matrix->dst are 0 */
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }

  free(len_gb_xn);
  free(start_cf_gb_xn);
  free(div_xn);

  return matrix;
}

/**

   lmb is the monomial basis (of the subscpace of the quotient ring)
   given by ascending order.

   dquo is the dimension of this subspace.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

 **/
static inline sp_matfglmcol_t *
build_matrixn_colon(int32_t *lmb, long dquot, int32_t bld,
		    int32_t **blen, int32_t **bexp, int32_t *bcf, int32_t *bexp_lm,
		    bs_t *tbr, ht_t *bht, md_t *st,
		    const exp_t * const mul, bs_t * bs,
		    const int nv, const long fc, const long maxdeg,
		    const data_gens_ff_t *gens,
		    uint32_t * leftvector,
		    uint32_t ** leftvectorsparam, long suppsize){

  const len_t ebl = bht->ebl;
  const len_t evl = bht->evl;
  int32_t *evi    =  (int *)malloc((unsigned long)nv * sizeof(int));
  if (ebl == 0) {
    for (long i = 1; i < evl; ++i) {
      evi[i-1]    =   i;
    }
  } else {
    for (long i = 1; i < ebl; ++i) {
      evi[i-1]    =   i;
    }
    for (long i = ebl+1; i < evl; ++i) {
      evi[i-2]    =   i;
    }
  }
  copy_extrapoly_in_vector(leftvector, dquot, lmb, 1,
			   tbr, bht, evi, st, nv, maxdeg);

  /* takes monomials in bexp_lm which are reducible by xn */
  /* div_xn contains the indices of those monomials*/
  int32_t *div_xn = calloc(bld, sizeof(int32_t));
  /* div_xn contains the indices of those not reducible by xn*/
  int32_t *div_not_xn = calloc(bld, sizeof(int32_t));

  /* len_xn is the number of those monomials of degree at most maxdeg+1. */
  long len_not_xn = 0;
  long len_xn = get_div_xn_bounded(bexp_lm, bld, nv, div_xn,div_not_xn,
				   &len_not_xn,maxdeg+1);

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the Gb) "
	  "which are divisible by x_n "
	  "and with bounded degree: %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the Gb) "
	  "which are not divisible by x_n "
	  "and with bounded degree: %ld\n", len_not_xn);
  for(long i=0; i < len_not_xn; i++){
    fprintf(stderr, "%d, ", div_not_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  long count_lm = 0;
  /* list of monomials in the staircase that leave the staircase after
     multiplication by xn and land on a multiple of a leading monomial of
     the Gb */
  long *extra_nf= calloc (dquot, sizeof(long));
  long count_not_lm = 0;
  /* list of monomials in the staircase that leave the staircase after
     multiplication by xn and land on zero */
  long *zeronf= calloc (dquot, sizeof(long));
  long count_zero = 0;
  for (long i = 0; i < dquot; i++) {
        long pos = -1;
	int32_t *exp = lmb + (i * nv);
	if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX>0
	  display_monomial_full(stderr, nv, NULL, 0, exp);
	  fprintf(stderr, " => remains in monomial basis\n");
#endif
	}
	else{
	  /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
	  display_monomial_full(stderr, nv, NULL, 0, exp);
	  fprintf(stderr, " => does NOT remain in monomial basis");
#endif
	  if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count_lm])*nv, nv)){
	    count_lm++;
#if DEBUGBUILDMATRIX > 0
	    fprintf(stderr, " => lands on a leading monomial\n");
#endif
	  }
	  else{
	    int is_divisible = 0;
	    for(long j = 0; j < len_xn; j++) {
	      if(is_divisible_exponent_xxn(exp, bexp_lm+(div_xn[j])*nv, nv)){
		extra_nf[count_not_lm]=i;
		count_not_lm++;
#if DEBUGBUILDMATRIX > 0
		fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
		is_divisible = 1;
		break;
	      }
	    }
	    if (!is_divisible) {
#if DEBUGBUILDMATRIX > 0
	      fprintf(stderr, " => lands on 0\n");
#endif
	      zeronf[count_zero]=i;
	      count_zero++;
	    }
	  }
	}
  }

  printf ("Number of extra normal forms for the matrix to compute: %ld\n",count_not_lm);
  printf ("Number of extra normal forms for the vectors to compute: %d\n",2*nv-2);
  /* Computation of the extra normal forms */
  /* count_not_lm monomials to reduce
   * each has length 1 */
  /* 2*(nv-1) shifts of phi to reduce
   * each has the same length as phi */
  long tobereduced = count_not_lm + 2*nv-2;
  int32_t* lens=(int32_t *) (malloc(sizeof(int32_t) * tobereduced));
  int32_t* exps = (int32_t *) (malloc(sizeof(int32_t) * (count_not_lm +
							 suppsize * (2*nv-2)) * nv));
  int32_t* cfs = (int32_t *) (malloc(sizeof(int32_t) * (count_not_lm +
							suppsize * (2*nv-2))));
  /* pure monomials to be reduced */
  for (long i = 0; i < count_not_lm;i++){
    lens[i]=1;
    cfs[i]=1;
    long j= extra_nf[i];
    for (long k = 0; k < nv-1; k++) {
      exps[i*nv+k]=lmb[j*nv+k];
    }
    exps[i*nv+nv-1]=lmb[j*nv+nv-1]+1;
  }

  /* shifts of to be reduced */
  len_t idx = tbr->lmps[1];
  /* printf ("idx=%d\n",idx); */
  len_t * hm  = tbr->hm[idx]+OFFSET;
  /* len_t len = tbr->hm[idx][LENGTH]; */
  for (long i = 0; i < 2*nv-2;i++){
    lens[count_not_lm + i]=suppsize;
    for (long j = 0; j < suppsize; j++) {
      cfs[count_not_lm + i*suppsize + j]=tbr->cf_32[tbr->hm[idx][COEFFS]][j];
      for (long k = 0; k < nv-1; k++) {
	int32_t cpt = 0;
	if (i == k) {
	  cpt = 1;
	}
	else if (i == nv - 1 + k) {
	  cpt = 2;
	}
	exps[(count_not_lm + i*suppsize+j)*nv+k]=bht->ev[hm[j]][evi[k]] + cpt;
      }
      exps[(count_not_lm + i*suppsize+j)*nv+nv-1]=bht->ev[hm[j]][evi[nv-1]];
    }
  }
  tbr = initialize_basis(st);
  tbr->ht = bht;
#if REDUCTION_ALLINONE
  import_input_data(tbr, st, 0, tobereduced, lens, exps, (void *)cfs, NULL);
  tbr->ld = tbr->lml  =  tobereduced;
  /* printf ("%ld imported\n",tobereduced); */
  for (int k = 0; k < tobereduced; ++k) {
    tbr->lmps[k]  = k; /* fix input element in tbr */
  }
  /* printf ("polynomials to be reduced\n"); */
  /* print_msolve_polynomials_ff(stdout, 0, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
  int32_t err = 0;
  tbr = core_nf(tbr, st, mul, bs, &err);
  if (err) {
    printf("Problem with normalform, stopped computation.\n");
    exit(1);
  }
  /* printf ("reductions\n"); */
  /* print_msolve_polynomials_ff(stdout, tobereduced, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
#endif
  printf ("Number of zero normal forms: %ld\n",count_zero);
  /* lengths of the polys which we need to build the matrix */
  int32_t *len_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  long pos = 0, k = 0;
  for(long i = 0; i < bld; i++){
    if(i==div_xn[k]){
      len_gb_xn[k]=(*blen)[i];
      start_cf_gb_xn[k]=pos;
      pos+=(*blen)[i];
      k++;
    }
    else{
      pos+=(*blen)[i];
    }
  }

#if DEBUGBUILDMATRIX > 0
  fprintf(stderr, "Length of polynomials whose leading terms are divisible by x_n\n");
  for(long i = 0; i < len_xn-1; i++){
    fprintf(stderr, "%u, ", len_gb_xn[i]);
  }
  fprintf(stderr, "%u\n", len_gb_xn[len_xn-1]);
#endif
  sp_matfglmcol_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglmcol_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nzero = count_zero;
  matrix->nrows = len_xn + count_not_lm;
  long len1 = dquot * (len_xn + count_not_lm);
  long len2 = dquot - (len_xn + count_not_lm + count_zero);

  if(posix_memalign((void **)&(matrix->dense_mat), 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->zero_idx, 32, sizeof(CF_t)*count_zero)){
    fprintf(stderr, "Problem when allocating matrix->zero_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < count_zero; i++){
      matrix->zero_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32,
		    sizeof(CF_t)*(len_xn + count_not_lm))){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn + count_not_lm; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*(len_xn + count_not_lm))){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn + count_not_lm; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long l_zero = 0;
  long nrows = 0;
  long count = 0;
  long count_nf = 0;
  for(long i = 0; i < dquot; i++){
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
    /* display_monomial_full(stderr, nv, NULL, 0, exp); */
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => remains in monomial basis\n");
#endif
      /* mult by xn stays in the basis */
      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{
      /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => does NOT remain in monomial basis");
#endif
      if (i == zeronf[l_zero]){
#if DEBUGBUILDMATRIX > 0
	fprintf(stderr, " => lands on 0\n");
	printf ("zero row #%ld in row #%ld\n",l_zero,i);
#endif
	matrix->zero_idx[l_zero] = i;
	l_zero++;
      }
      else {
	matrix->dense_idx[l_dens] = i;
	l_dens++;
	if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
#if DEBUGBUILDMATRIX > 0
	  fprintf(stderr, " => lands on a leading monomial\n");
#endif
	  copy_poly_in_matrixcol(matrix, nrows, bcf, bexp, blen,
				 start_cf_gb_xn[count], len_gb_xn[count], lmb,
				 nv, fc);
	  nrows++;
	  count++;
	  if(len_xn < count && i < dquot){
	    fprintf(stderr, "One should not arrive here (build_matrix)\n");
	    free(lens);
	    free(exps);
	    free(cfs);
	    free(matrix->dense_mat);
	    free(matrix->dense_idx);
	    free(matrix->triv_idx);
	    free(matrix->triv_pos);
	    free(matrix->zero_idx);
	    free(matrix->dst);
	    free(matrix);
	    free(evi);
	    free(len_gb_xn);
	    free(start_cf_gb_xn);
	    free(div_xn);
	    free(div_not_xn);
	    return NULL;
	  }
	}
	else if (i == extra_nf[count_nf]){
#if DEBUGBUILDMATRIX > 0
	  fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
#if REDUCTION_ALLINONE
	  copy_extrapoly_in_matrixcol(matrix, nrows, lmb,
				      tobereduced + count_nf,
				      tbr, bht, evi, st, nv, maxdeg);
#else /* reduction one by one */
	  int success;
	  import_input_data(tbr, st, count_nf, count_nf+1, lens, exps, (void *)cfs, NULL);
	  tbr->ld = tbr->lml  =  1;
	  /* printf ("1 imported\n"); */
	  for (int k = 0; k < 1; ++k) {
	    tbr->lmps[k]  = k; /* fix input element in tbr */
	  }
	  /* printf ("mononomial to be reduced\n"); */
	  /* print_msolve_polynomials_ff(stdout, 0, tbr->lml, tbr, bht, */
	  /* 		      st, gens->vnames, 0); */
	  success = core_nf(&tbr, &bht, &st, mul, bs);
	  if (!success) {
	    printf("Problem with normalform, stopped computation.\n");
	    exit(1);
	  }
	  /* printf ("reduction\n"); */
	  /* print_msolve_polynomials_ff(stdout, 1, tbr->lml, tbr, bht, */
	  /* 		      st, gens->vnames, 0); */
	  copy_extrapoly_in_matrixcol(matrix, nrows, lmb,
				      1,
				      tbr, bht, evi, st, nv, maxdeg);
#endif
	  nrows++;
	  count_nf++;
	}
      }
    }
  }
  /* printf ("matrix finished\n[\n"); */
  /* assumes the entries of matrix->dst are 0 */
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
	/* printf ("%d, ",matrix->dst[i]); */
      }
      else{
        break;
      }
    }
  }

#if REDUCTION_ALLINONE
  for (long i = 0; i < 2*nv - 2; i++) {
    copy_extrapoly_in_vector(leftvectorsparam[i], dquot, lmb,
			     tobereduced + count_not_lm + i,
			     tbr, bht, evi, st, nv, maxdeg);
  }
#else /* the shifts of phi are now reduced */
  tbr = initialize_basis(st);
  tbr->ht = bht;
  import_input_data(tbr, st, count_not_lm, tobereduced, lens, exps, (void *)cfs, NULL);
  tbr->ld = tbr->lml  =  2*nv-2;
  /* printf ("%d imported\n",2*nv-2); */
  for (int k = 0; k < 2*nv-2; ++k) {
    tbr->lmps[k]  = k; /* fix input element in tbr */
  }
  /* printf ("shifts of phi to be reduced\n"); */
  /* print_msolve_polynomials_ff(stdout, 0, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
  int success = core_nf(&tbr, &bht, &st, mul, bs);
  if (!success) {
    printf("Problem with normalform, stopped computation.\n");
    exit(1);
  }
  /* printf ("reductions\n"); */
  /* print_msolve_polynomials_ff(stdout, 2*nv-2, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
  for (long i = 0; i < 2*nv - 2; i++) {
    copy_extrapoly_in_vector(leftvectorsparam[i], dquot, lmb,
			     2*nv-2+i,
			     tbr, bht, evi, st, nv, maxdeg);
  }
#endif
  free(lens);
  free(exps);
  free(cfs);
  free(evi);
  free(len_gb_xn);
  free(start_cf_gb_xn);
  free(div_xn);
  free(div_not_xn);
  return matrix;
}

/**

   lmb is the monomial basis (of the subscpace of the quotient ring)
   given by ascending order.

   dquo is the dimension of this subspace.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

 **/
static inline sp_matfglmcol_t *
build_matrixn_colon_no_zero(int32_t *lmb, long dquot, int32_t bld,
			    int32_t **blen, int32_t **bexp,
			    int32_t *bcf, int32_t *bexp_lm, bs_t *tbr,
			    ht_t *bht, md_t *st,
			    const exp_t * const mul, bs_t * bs,
			    const int nv, const long fc, const long maxdeg,
			    const data_gens_ff_t *gens,
			    uint32_t * leftvector,
			    uint32_t ** leftvectorsparam, long suppsize){

  const len_t ebl = bht->ebl;
  const len_t evl = bht->evl;
  int32_t *evi    =  (int *)malloc((unsigned long)nv * sizeof(int));
  if (ebl == 0) {
    for (long i = 1; i < evl; ++i) {
      evi[i-1]    =   i;
    }
  } else {
    for (long i = 1; i < ebl; ++i) {
      evi[i-1]    =   i;
    }
    for (long i = ebl+1; i < evl; ++i) {
      evi[i-2]    =   i;
    }
  }
  for (long i = 1; i < tbr->lml; i++) {
    len_t idx = tbr->lmps[i];
    len_t * hm  = tbr->hm[idx]+OFFSET;
    len_t len = tbr->hm[idx][LENGTH];
    long k = 0;
    for (long j = 0; j < 5; j++) {
      while (!is_equal_exponent_bs (bht,hm[len-1-j],evi,lmb+k*nv,nv)) {
	k++;
      }
    }
  }
  copy_extrapoly_in_vector_no_zero(leftvector, dquot, lmb, 1,
				   tbr, bht, evi, st, nv, maxdeg);
  /* takes monomials in bexp_lm which are reducible by xn */
  /* div_xn contains the indices of those monomials*/
  int32_t *div_xn = calloc(bld, sizeof(int32_t));
  /* div_xn contains the indices of those not reducible by xn*/
  int32_t *div_not_xn = calloc(bld, sizeof(int32_t));

  /* len_xn is the number of those monomials of degree at most maxdeg+1. */
  /* int is_divisible = 0; */
  long len_not_xn = 0;
  long len_xn = get_div_xn_bounded(bexp_lm, bld, nv, div_xn,div_not_xn,
				   &len_not_xn,maxdeg+1);

#if DEBUGBUILDMATRIX > 0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the Gb) "
	  "which are divisible by x_n "
	  "and with bounded degree: %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the Gb) "
	  "which are not divisible by x_n "
	  "and with bounded degree: %ld\n", len_not_xn);
  for(long i=0; i < len_not_xn; i++){
    fprintf(stderr, "%d, ", div_not_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  long count_lm = 0;
  /* list of monomials in the staircase that leave the staircase after
     multiplication by xn and land on a multiple of a leading monomial of
     the Gb */
  long *extra_nf= calloc (dquot, sizeof(long));
  long count_not_lm = 0;
  /* list of monomials in the staircase that leave the staircase after
     multiplication by xn and land on zero */
  /* long *zeronf= calloc (dquot, sizeof(long)); */
  /* long count_zero = 0; */
  for (long i = 0; i < dquot; i++) {
        long pos = -1;
	int32_t *exp = lmb + (i * nv);
	if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
	  display_monomial_full(stderr, nv, NULL, 0, exp);
	  fprintf(stderr, " => remains in monomial basis\n");
#endif
	}
	else{
	  /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
	  display_monomial_full(stderr, nv, NULL, 0, exp);
	  fprintf(stderr, " => does NOT remain in monomial basis");
#endif
	  if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count_lm])*nv, nv)){
	    count_lm++;
#if DEBUGBUILDMATRIX > 0
	    fprintf(stderr, " => lands on a leading monomial\n");
#endif
	  }
	  else{
	    for(long j = 0; j < len_xn; j++) {
	      if(is_divisible_exponent_xxn(exp, bexp_lm+(div_xn[j])*nv, nv)){
		extra_nf[count_not_lm]=i;
		count_not_lm++;
#if DEBUGBUILDMATRIX > 0
		fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
		/* is_divisible = 1; */
		break;
	      }
	    }
/* 	    if (!is_divisible) { */
/* #if DEBUGBUILDMATRIX > 0 */
/* 	      fprintf(stderr, " => land on 0\n"); */
/* #endif */
/* 	      zeronf[count_zero]=i; */
/* 	      count_zero++; */
/* 	    } */
	  }
	}
  }

  printf ("Number of extra normal forms for the matrix to compute: %ld\n",count_not_lm);
  printf ("Number of extra normal forms for the vectors to compute: %d\n",2*nv-2);
  /* Computation of the extra normal forms */
  /* count_not_lm monomials to reduce
   * each has length 1 */
  /* 2*(nv-1) shifts of phi to reduce
   * each has the same length as phi */
  long tobereduced = count_not_lm + 2*nv-2;
  int32_t* lens=(int32_t *) (malloc(sizeof(int32_t) * tobereduced));
  int32_t* exps = (int32_t *) (malloc(sizeof(int32_t) * (count_not_lm +
							 suppsize * (2*nv-2)) * nv));
  int32_t* cfs = (int32_t *) (malloc(sizeof(int32_t) * (count_not_lm +
							suppsize * (2*nv-2))));
  /* pure monomials to be reduced */
  for (long i = 0; i < count_not_lm;i++){
    lens[i]=1;
    cfs[i]=1;
    long j= extra_nf[i];
    for (long k = 0; k < nv-1; k++) {
      exps[i*nv+k]=lmb[j*nv+k];
    }
    exps[i*nv+nv-1]=lmb[j*nv+nv-1]+1;
  }

  /* shifts of to be reduced */
  len_t idx = tbr->lmps[1];
  len_t * hm  = tbr->hm[idx]+OFFSET;
  /* len_t len = tbr->hm[idx][LENGTH]; */
  for (long i = 0; i < 2*nv-2;i++){
    lens[count_not_lm + i]=suppsize;
    for (long j = 0; j < suppsize; j++) {
      cfs[count_not_lm + i*suppsize + j]=tbr->cf_32[tbr->hm[idx][COEFFS]][j];
      for (long k = 0; k < nv-1; k++) {
	int32_t cpt = 0;
	if (i == k) {
	  cpt = 1;
	}
	else if (i == nv - 1 + k) {
	  cpt = 2;
	}
	exps[(count_not_lm + i*suppsize+j)*nv+k]=bht->ev[hm[j]][evi[k]] + cpt;
      }
      exps[(count_not_lm + i*suppsize+j)*nv+nv-1]=bht->ev[hm[j]][evi[nv-1]];
    }
  }
  tbr = initialize_basis(st);
  tbr->ht = bht;
  import_input_data(tbr, st, 0, tobereduced, lens, exps, (void *)cfs, NULL);
  tbr->ld = tbr->lml  =  tobereduced;
  /* printf ("%ld imported\n",tobereduced); */
  for (int k = 0; k < tobereduced; ++k) {
    tbr->lmps[k]  = k; /* fix input element in tbr */
  }
  /* printf ("polynomials to be reduced\n"); */
  /* print_msolve_polynomials_ff(stdout, 0, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
  int32_t err = 0;
  tbr = core_nf(tbr, st, mul, bs, &err);
  if (err) {
    printf("Problem with normalform, stopped computation.\n");
    exit(1);
  }
  /* printf ("reductions\n"); */
  /* print_msolve_polynomials_ff(stdout, tobereduced, tbr->lml, tbr, bht, */
  /* 			      st, gens->vnames, 0); */
  /* printf ("Number of zero normal forms: %ld\n",count_zero); */
  /* lengths of the polys which we need to build the matrix */
  int32_t *len_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  long pos = 0, k = 0;
  for(long i = 0; i < bld; i++){
    if(i==div_xn[k]){
      len_gb_xn[k]=(*blen)[i];
      start_cf_gb_xn[k]=pos;
      pos+=(*blen)[i];
      k++;
    }
    else{
      pos+=(*blen)[i];
    }
  }

#if DEBUGBUILDMATRIX > 0
  fprintf(stderr, "Length of polynomials whose leading terms are divisible by x_n\n");
  for(long i = 0; i < len_xn-1; i++){
    fprintf(stderr, "%u, ", len_gb_xn[i]);
  }
  fprintf(stderr, "%u\n", len_gb_xn[len_xn-1]);
#endif
  sp_matfglmcol_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglmcol_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nzero = 0; /* count_zero; */
  matrix->nrows = len_xn + count_not_lm;
  long len1 = dquot * (len_xn + count_not_lm);
  long len2 = dquot - (len_xn + count_not_lm /* + count_zero*/);

  if(posix_memalign((void **)&(matrix->dense_mat), 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  /* if(posix_memalign((void **)&matrix->zero_idx, 32, sizeof(CF_t)*count_zero)){ */
  /*   fprintf(stderr, "Problem when allocating matrix->zero_idx\n"); */
  /*   exit(1); */
  /* } */
  /* else{ */
  /*   for(long i = 0; i < count_zero; i++){ */
  /*     matrix->zero_idx[i] = 0; */
  /*   } */
  /* } */
  if(posix_memalign((void **)&matrix->dense_idx, 32,
		    sizeof(CF_t)*(len_xn + count_not_lm))){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn + count_not_lm; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*(len_xn + count_not_lm))){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn + count_not_lm; i++){
      matrix->dst[i] = 0;
    }
  }

  for (long i = count_not_lm; i  < tbr->lml; i++) {
    len_t idx = tbr->lmps[i];
    len_t * hm  = tbr->hm[idx]+OFFSET;
    len_t len = tbr->hm[idx][LENGTH];
    long k = 0;
    for (long j = 0; j < 5; j++) {
      while (!is_equal_exponent_bs (bht,hm[len-1-j],evi,lmb+k*nv,nv)) {
	k++;
      }
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  /* long l_zero = 0; */
  long nrows = 0;
  long count = 0;
  long count_nf = 0;
  for(long i = 0; i < dquot; i++){
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
    /* display_monomial_full(stderr, nv, NULL, 0, exp); */
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => remains in monomial basis\n");
#endif
      /* mult by xn stays in the basis */
      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{
      /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => does NOT remain in monomial basis");
#endif
/*       if (i == zeronf[l_zero]){ */
/* #if DEBUGBUILDMATRIX > 0 */
/* 	fprintf(stderr, " => land on 0\n"); */
/* 	printf ("zero row #%ld in row #%ld\n",l_zero,i); */
/* #endif */
/* 	matrix->zero_idx[l_zero] = i; */
/* 	l_zero++; */
/*       } */
/*       else { */
	matrix->dense_idx[l_dens] = i;
	l_dens++;
	if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
#if DEBUGBUILDMATRIX > 0
	  fprintf(stderr, " => lands on a leading monomial\n");
#endif
	  copy_poly_in_matrixcol_no_zero(matrix, nrows, bcf, bexp, blen,
					 start_cf_gb_xn[count],
					 len_gb_xn[count], lmb, nv, fc);
	  nrows++;
	  count++;
	  if(len_xn < count && i < dquot){
	    fprintf(stderr, "One should not arrive here (build_matrix)\n");
	    free(lens);
	    free(exps);
	    free(cfs);
	    free(matrix->dense_mat);
	    free(matrix->dense_idx);
	    free(matrix->triv_idx);
	    free(matrix->triv_pos);
	    free(matrix->dst);
	    /* free(matrix->zero_idx); */
	    free(matrix);
	    free(evi);
	    free(len_gb_xn);
	    free(start_cf_gb_xn);
	    free(div_xn);
	    free(div_not_xn);
	    return NULL;
	  }
	}
	else if (i == extra_nf[count_nf]){
#if DEBUGBUILDMATRIX > 0
	  fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
	  copy_extrapoly_in_matrixcol_no_zero(matrix, nrows, lmb,
					      tobereduced + count_nf,
					      tbr, bht, evi, st, nv, maxdeg);
	  nrows++;
	  count_nf++;
	}
    }
  }
  /* printf ("matrix finished\n[\n"); */
  /* for (long i = 0; i< matrix->nrows; i++){ */
  /*   printf ("["); */
  /*   for (long j = 0; j < matrix->ncols; j++){ */
  /*     printf ("%d, ", matrix->dense_mat[i*matrix->ncols+j]); */
  /*   } */
  /*   printf ("]\n"); */
  /* } */
  /* printf("]\n"); */
  /* assumes the entries of matrix->dst are 0 */
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
	/* printf ("%d, ",matrix->dst[i]); */
      }
      else{
        break;
      }
    }
  }

  for (long i = 0; i < 2*nv - 2; i++) {
    copy_extrapoly_in_vector_no_zero(leftvectorsparam[i], dquot, lmb,
				     tobereduced + count_not_lm + i,
				     tbr, bht, evi, st, nv, maxdeg);
  }
  free(lens);
  free(exps);
  free(cfs);
  free(evi);
  free(len_gb_xn);
  free(start_cf_gb_xn);
  free(div_xn);
  free(div_not_xn);
  return matrix;
}

/**

   lmb is the monomial basis (of the quotient ring) given by ascending order.

   dquo is the dimension of the quotient.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

 **/
static inline sp_matfglm_t * build_matrixn_trace(int32_t **bdiv_xn,
                                                 int32_t **blen_gb_xn,
                                                 int32_t **bstart_cf_gb_xn,
                                                 int32_t *lmb, long dquot,
                                                 int32_t bld,
                                                 int32_t **blen, int32_t **bexp,
                                                 int32_t *bcf,
                                                 int32_t *bexp_lm,
                                                 const int nv, const long fc){


  *bdiv_xn = calloc((unsigned long)bld, sizeof(int32_t));
  int32_t *div_xn = *bdiv_xn;

  long len_xn = get_div_xn(bexp_lm, bld, nv, div_xn);

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the staircase) which are divisible by x_n = %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
#endif


  *blen_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *len_gb_xn = *blen_gb_xn;

  *bstart_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = *bstart_cf_gb_xn;

  long pos = 0, k = 0;
  for(long i = 0; i < bld; i++){
    if(i==div_xn[k]){
      len_gb_xn[k]=(*blen)[i];
      start_cf_gb_xn[k]=pos;
      pos+=(*blen)[i];
      k++;
    }
    else{
      pos+=(*blen)[i];
    }
  }

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "Length ofpolynomials whose leading terms is divisible by x_n\n");
  for(long i = 0; i < len_xn; i++){
    fprintf(stderr, "%d, ", len_gb_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  sp_matfglm_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nrows = len_xn;
  long len1 = dquot * len_xn;
  long len2 = dquot - len_xn;

  if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  for(long i = 0; i < dquot; i++){

    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
    //    fprintf(stderr, "\n");
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{

#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
        copy_poly_in_matrix(matrix, nrows, bcf, bexp, blen,
                            start_cf_gb_xn[count], len_gb_xn[count], lmb,
                            nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
          fprintf(stderr, "One should not arrive here (build_matrix)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
          return NULL;
        }
      }
      else{
        fprintf(stderr, "Staircase is not generic\n");
        fprintf(stderr, "Multiplication by ");
        display_monomial_full(stderr, nv, NULL, 0, exp);
        fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return NULL;
      }
    }
  }

  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }

  return matrix;
}





/**

   lmb is the monomial basis (of the quotient ring) given by ascending order.

   dquo is the dimension of the quotient.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

 **/
static inline sp_matfglm_t * build_matrixn_from_bs(int32_t *lmb, long dquot,
                                                   bs_t *bs,
                                                   ht_t *ht,
                                                   int32_t *bexp_lm,
                                                   const int nv, const long fc){


  int32_t *div_xn = calloc(bs->lml, sizeof(int32_t));

  long len_xn = get_div_xn(bexp_lm, bs->lml, nv, div_xn);

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the staircase) which are divisible by x_n = %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
#endif

  int32_t *len_gb_xn = malloc(sizeof(int32_t) * len_xn);

  int32_t *start_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  long pos = 0, k = 0;
  for(long i = 0; i < bs->lml; i++){
    long len = bs->hm[bs->lmps[i]][LENGTH];
    if(i==div_xn[k]){
      len_gb_xn[k]=len;
      start_cf_gb_xn[k]=pos;
      pos+=len;
      k++;
    }
    else{
      pos+=len;
    }
  }

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "Length ofpolynomials whose leading terms is divisible by x_n\n");
  for(long i = 0; i < len_xn; i++){
    fprintf(stderr, "%d, ", len_gb_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  sp_matfglm_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nrows = len_xn;
  long len1 = dquot * len_xn;
  long len2 = dquot - len_xn;

  if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  for(long i = 0; i < dquot; i++){

    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
    //    fprintf(stderr, "\n");
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{

#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
        copy_poly_in_matrix_from_bs(matrix, nrows, bs, ht, //bcf, bexp, blen,
                                    div_xn[count], len_gb_xn[count],
                                    start_cf_gb_xn[count], len_gb_xn[count], lmb,
                                    nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
          fprintf(stderr, "One should not arrive here (build_matrix with trace)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
          return NULL;
        }
      }
      else{
        fprintf(stderr, "Staircase is not generic\n");
        fprintf(stderr, "Multiplication by ");
#if DEBUGBUILDMATRIX > 0
        display_monomial_full(stderr, nv, NULL, 0, exp);
#endif
        fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return NULL;
      }
    }
  }

  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }

  free(len_gb_xn);
  free(start_cf_gb_xn);
  free(div_xn);

  return matrix;
}

/**

   lmb is the monomial basis (of the quotient ring) given by ascending order.

   dquo is the dimension of the quotient.

   data of gb are given by ascending order.

   bexp_lm is the leading monomials of gb ; there are bld[0] of them.

   The matrix is already allocated.

   div_xn tableau des monomes divisibles par xn
   len_gb_xn + startc_cf_gb_xn permettent d'avoir les longueurs des
   polys dont les monomes de tete sont divisibles par xn
 **/
static inline void build_matrixn_from_bs_trace_application(sp_matfglm_t *matrix,
                                                           int32_t *div_xn,
                                                           int32_t *len_gb_xn,
                                                           int32_t *start_cf_gb_xn,
                                                           int32_t *lmb, long dquot,
                                                           bs_t *bs,
                                                           ht_t *ht,
                                                           int32_t *bexp_lm,
							   const int nv,
                                                           const long fc){

  long len_xn = matrix->nrows; //get_div_xn(bexp_lm, bs->lml, nv, div_xn);

  matrix->charac = fc;
  /* matrix->ncols = dquot; */
  /* matrix->nrows = len_xn; */
  long len1 = dquot * matrix->nrows;
  long len2 = dquot - matrix->nrows;

  for(long i = 0; i < len1; i++){
    matrix->dense_mat[i] = 0;
  }
  for(long i = 0; i < len2; i++){
    matrix->triv_idx[i] = 0;
  }
  for(long i = 0; i < len2; i++){
    matrix->triv_pos[i] = 0;
  }
  for(long i = 0; i < len_xn; i++){
    matrix->dense_idx[i] = 0;
  }
  for(long i = 0; i < len_xn; i++){
    matrix->dst[i] = 0;
  }

  long pos = 0, k = 0;
  for(long i = 0; i < bs->lml; i++){
    long len = bs->hm[bs->lmps[i]][LENGTH];
    if(i==div_xn[k]){
      len_gb_xn[k]=len;
      start_cf_gb_xn[k]=pos;
      pos+=len;
      k++;
    }
    else{
      pos+=len;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  for(long i = 0; i < dquot; i++){

    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{

#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
        copy_poly_in_matrix_from_bs(matrix, nrows, bs, ht, //bcf, bexp, blen,
                                    div_xn[count], len_gb_xn[count],
                                    start_cf_gb_xn[count], len_gb_xn[count], lmb,
                                    nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
          fprintf(stderr, "One should not arrive here (build_matrix with trace)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
          exit(1);
        }
      }
      else{
        fprintf(stderr, "Staircase is not generic\n");
        fprintf(stderr, "Multiplication by ");
        display_monomial_full(stderr, nv, NULL, 0, exp);
        fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return ;
        //        exit(1);
      }
    }
  }
  //Ici on support que les entres de matrix->dst sont initialisees a 0
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }
}

static inline void build_matrixn_unstable_from_bs_trace_application(sp_matfglm_t *matrix,
								    int32_t *div_xn,
								    int32_t *len_gb_xn,
								    int32_t *start_cf_gb_xn,
								    long *extra_nf,
								    int32_t *lens_extra_nf,
								    int32_t *exps_extra_nf,
								    int32_t *cfs_extra_nf,
								    int32_t *lmb, long dquot,
								    bs_t *bs,
								    ht_t *ht,
								    int32_t *bexp_lm,
								    const md_t const *st,
								    const int nv,
								    const long fc,
								    const int thread_number){
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;
  int32_t *evi    =  (int *)malloc((unsigned long)nv * sizeof(int));
  if (ebl == 0) {
    for (long i = 1; i < evl; ++i) {
      evi[i-1]    =   i;
    }
  } else {
    for (long i = 1; i < ebl; ++i) {
      evi[i-1]    =   i;
    }
    for (long i = ebl+1; i < evl; ++i) {
      evi[i-2]    =   i;
    }
  }
  bs_t *tbr;
  long count_not_lm = matrix->nnfs;
  if (count_not_lm) {
    md_t *md = copy_meta_data(st,fc);
    tbr = initialize_basis(md);
    exp_t *mul = (exp_t *)calloc(ht->evl, sizeof(exp_t));
    tbr->ht = ht;
    /* reduction */
    import_input_data(tbr, md, 0, count_not_lm, lens_extra_nf, exps_extra_nf,
		      (void *)cfs_extra_nf, NULL);
    tbr->ld = tbr->lml  =  count_not_lm;
    for (long k = 0; k < count_not_lm; ++k) {
      tbr->lmps[k]  = k; /* fix input element in bs */
    }
    int32_t err = 0;
    tbr = core_nf(tbr, md, mul, bs, &err);
    if (err) {
      printf("Problem with normalform, stopped computation.\n");
      exit(1);
    }
    free(mul);
    free(md);
  }
  long len0 = matrix->nrows;
  long len_xn = len0-count_not_lm; //get_div_xn(bexp_lm, bs->lml, nv, div_xn);

  matrix->charac = fc;
  long len1 = dquot * len0;
  long len2 = dquot - len0;

  for(long i = 0; i < len1; i++){
    matrix->dense_mat[i] = 0;
  }
  for(long i = 0; i < len2; i++){
    matrix->triv_idx[i] = 0;
  }
  for(long i = 0; i < len2; i++){
    matrix->triv_pos[i] = 0;
  }
  for(long i = 0; i < len0; i++){
    matrix->dense_idx[i] = 0;
  }
  for(long i = 0; i < len0; i++){
    matrix->dst[i] = 0;
  }

  long pos = 0, k = 0;
  for(long i = 0; i < bs->lml; i++){
    long len = bs->hm[bs->lmps[i]][LENGTH];
    if(i==div_xn[k]){
      len_gb_xn[k]=len;
      start_cf_gb_xn[k]=pos;
      pos+=len;
      k++;
    }
    else{
      pos+=len;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  long count_nf = 0;
  
  for(long i = 0; i < dquot; i++){
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{

#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
        copy_poly_in_matrix_from_bs(matrix, nrows, bs, ht, //bcf, bexp, blen,
                                    div_xn[count], len_gb_xn[count],
                                    start_cf_gb_xn[count], len_gb_xn[count], lmb,
                                    nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
          fprintf(stderr, "One should not arrive here (build_matrix with trace)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

	  free_basis_without_hash_table(&tbr);
	  free(cfs_extra_nf);
	  free(exps_extra_nf);
	  free(lens_extra_nf);
	  free(extra_nf);
          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
	  free(evi);
          exit(1);
        }
      }
      else if (i == extra_nf[count_nf]){
#if DEBUGBUILDMATRIX > 0
	fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
	copy_nf_in_matrix_from_bs(matrix, nrows, count_nf, lmb,
				  tbr, ht, evi, st, nv);
	nrows++;
	count_nf++;
	if (count_not_lm < count_nf && i < dquot) {
          fprintf(stderr, "One should not arrive here (build_matrix with trace)\n");
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

	  free_basis_without_hash_table(&tbr);
	  free(cfs_extra_nf);
	  free(exps_extra_nf);
	  free(lens_extra_nf);
	  free(extra_nf);
          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
	  free(evi);
          exit(1);
        }
      }
      else{ /* should not arrive here */
        fprintf(stderr, "Staircase is not generic\n");
        fprintf(stderr, "Multiplication by ");
        display_monomial_full(stderr, nv, NULL, 0, exp);
        fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

	free_basis_without_hash_table(&tbr);
	free(cfs_extra_nf);
	free(exps_extra_nf);
	free(lens_extra_nf);
	free(extra_nf);
	free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
	free(evi);
        return ;
        //        exit(1);
      }
    }
  }
  //Ici on suppose que les entres de matrix->dst sont initialisees a 0
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }
  if (count_not_lm) {
    free_basis_without_hash_table(&tbr);
  }
  if(st->info_level){
    fprintf(stderr, "[%lu, %lu], Free / Dense = %.2f%%\n",
            len0, len_xn,
            100*((double)len_xn / (double)len0));
  }
}


static inline sp_matfglm_t * build_matrixn_from_bs_trace(int32_t **bdiv_xn,
                                                         int32_t **blen_gb_xn,
                                                         int32_t **bstart_cf_gb_xn,
                                                         int32_t *lmb, long dquot,
                                                         bs_t *bs,
                                                         ht_t *ht,
                                                         int32_t *bexp_lm,
                                                         const int nv,
                                                         const long fc,
                                                         const int info_level){


  *bdiv_xn = calloc((unsigned long)bs->lml, sizeof(int32_t));
  int32_t *div_xn = *bdiv_xn;

  long len_xn = get_div_xn(bexp_lm, bs->lml, nv, div_xn);


#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the staircase) which are divisible by x_n = %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
#endif

  *blen_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *len_gb_xn = *blen_gb_xn;

  *bstart_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = *bstart_cf_gb_xn;

  long pos = 0, k = 0;
  for(long i = 0; i < bs->lml; i++){
    long len = bs->hm[bs->lmps[i]][LENGTH];
    if(i==div_xn[k]){
      len_gb_xn[k]=len;
      start_cf_gb_xn[k]=pos;
      pos+=len;
      k++;
    }
    else{
      pos+=len;
    }
  }

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "Length of polynomials whose leading terms is divisible by x_n\n");
  for(long i = 0; i < len_xn; i++){
    fprintf(stderr, "%d, ", len_gb_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  sp_matfglm_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nrows = len_xn;
  long len1 = dquot * len_xn;
  long len2 = dquot - len_xn;

  if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i <len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len_xn)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len_xn; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;

  for(long i = 0; i < dquot; i++){
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
    //    fprintf(stderr, "\n");
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      int boo = 0;
      if(count < len_xn){
        boo = is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv);
      }
      if(boo){
        matrix->dense_idx[l_dens] = i;
        l_dens++;
        copy_poly_in_matrix_from_bs(matrix, nrows, bs, ht,
                                    div_xn[count], len_gb_xn[count],
                                    start_cf_gb_xn[count], len_gb_xn[count], lmb,
                                    nv, fc);
        nrows++;
        count++;

        if(len_xn < count && i < dquot){
          if(info_level){
            fprintf(stderr, "Staircase is not generic (1 => explain better)\n");
          }
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
          return NULL;
        }
      }
      else{
        if(info_level){
          fprintf(stderr, "Staircase is not generic\n");
          fprintf(stderr, "Multiplication by ");
          display_monomial_full(stderr, nv, NULL, 0, exp);
          fprintf(stderr, " gets outside the staircase\n");
        }
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return NULL;
      }
    }
  }

  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }

  return matrix;
}

static inline sp_matfglm_t * build_matrixn_unstable_from_bs_trace(int32_t **bdiv_xn,
								  int32_t **blen_gb_xn,
								  int32_t **bstart_cf_gb_xn,
								  long **bextra_nf,
								  int32_t **blens_extra_nf,
								  int32_t **bexps_extra_nf,
								  int32_t **bcfs_extra_nf,
								  int32_t *lmb, long dquot,
								  bs_t *bs, ht_t *ht,
								  int32_t *bexp_lm,
								  const md_t const *st,
								  const int nv,
								  const long fc,
								  const int32_t unstable_staircase,
								  const int info_level){
  *bdiv_xn = calloc((unsigned long)bs->lml, sizeof(int32_t));
  int32_t *div_xn = *bdiv_xn;
  long len_xn = get_div_xn(bexp_lm, bs->lml, nv, div_xn);

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of monomials (in the Gröbner basis) which are divisible by x_n = %ld\n", len_xn);
  for(long i=0; i < len_xn; i++){
    fprintf(stderr, "%d, ", div_xn[i]);
  }
  fprintf(stderr, "\n");
#endif
  *blen_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *len_gb_xn = *blen_gb_xn;

  *bstart_cf_gb_xn = malloc(sizeof(int32_t) * len_xn);
  int32_t *start_cf_gb_xn = *bstart_cf_gb_xn;

  long pos = 0, k = 0;
  for(long i = 0; i < bs->lml; i++){
    long len = bs->hm[bs->lmps[i]][LENGTH];
    if(i==div_xn[k]){
      len_gb_xn[k]=len;
      start_cf_gb_xn[k]=pos;
      pos+=len;
      k++;
    }
    else{
      pos+=len;
    }
  }

#if DEBUGBUILDMATRIX>0
  fprintf(stderr, "Length of polynomials whose leading terms is divisible by x_n\n");
  for(long i = 0; i < len_xn; i++){
    fprintf(stderr, "%d, ", len_gb_xn[i]);
  }
  fprintf(stderr, "\n");
#endif

  long count_lm = 0;
  /* number of monomials whose normal forms are needed */
  long count_not_lm = 0;

  /* list of monomials in the staircase that leave the staircase after
     multiplication by xn and land on a multiple of a leading monomial of
     the Gb */
  /* at most dquot-len_xn new columns to compute */
  *bextra_nf = calloc (dquot-len_xn, sizeof(long));
  long *extra_nf = *bextra_nf;
  for (long i = 0; i < dquot; i++) {
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX>0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => remains in monomial basis\n");
#endif
    }
    else{
      /* we get now outside the basis */
#if DEBUGBUILDMATRIX > 0
      display_monomial_full(stderr, nv, NULL, 0, exp);
      fprintf(stderr, " => does NOT remain in monomial basis");
#endif
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count_lm])*nv, nv)){
!	count_lm++;
#if DEBUGBUILDMATRIX > 0
	fprintf(stderr, " => lands on a leading monomial\n");
#endif
      }
      else{
	int is_divisible = 0;
	for(long j = 0; j < len_xn; j++) {
	  if(is_divisible_exponent_xxn(exp, bexp_lm+(div_xn[j])*nv, nv)){
	    extra_nf[count_not_lm]=i;
	    count_not_lm++;
#if DEBUGBUILDMATRIX > 0
	    fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
	    is_divisible = 1;
	    break;
	  }
	}
	if (!is_divisible) {
#if DEBUGBUILDMATRIX > 0
	  fprintf(stderr, " => BUG\n");
#endif
	}
      }
    }
  }
  long len0 = len_xn + count_not_lm;

#if DEBUGBUILDMATRIX > 0
  printf ("Number of dense columns:      %ld\n",len0);
  printf ("Number of free dense columns: %ld\n",len_xn);
  printf ("Number of normal forms:       %ld\n",count_not_lm);
#endif

  long threshold = dquot-len_xn;
  
  switch(unstable_staircase) {
  case 0:
    threshold = 0;
    break;
  case 1:
    threshold = (threshold+1)/4; /* round to nearest integer */
    break;
  case 2:
    threshold = threshold/2;
    break;
  case 3:
    threshold = (3*threshold+1)/4; /* round to nearest integer */
    break;
  }
  
  if (count_not_lm > threshold) {
    if(info_level){
      fprintf(stderr, "Staircase is not generic\n");
      fprintf(stderr, "and too many normal forms need to be computed\n");
    }
    free(extra_nf);
    free(len_gb_xn);
    free(start_cf_gb_xn);
    free(div_xn);
    return NULL;
  }

  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;
  int32_t *evi    =  (int32_t *)malloc((unsigned long)nv * sizeof(int32_t));
  if (ebl == 0) {
    for (long i = 1; i < evl; ++i) {
      evi[i-1]    =   i;
    }
  } else {
    for (long i = 1; i < ebl; ++i) {
      evi[i-1]    =   i;
    }
    for (long i = ebl+1; i < evl; ++i) {
      evi[i-2]    =   i;
    }
  }

  *blens_extra_nf=(int32_t *) (malloc(sizeof(int32_t) * count_not_lm));
  int32_t* lens_extra_nf= *blens_extra_nf;
  *bexps_extra_nf = (int32_t *) (malloc(sizeof(int32_t) * count_not_lm * nv));
  int32_t* exps_extra_nf= *bexps_extra_nf;
  *bcfs_extra_nf = (int32_t *) (malloc(sizeof(int32_t) * count_not_lm));
  int32_t* cfs_extra_nf= *bcfs_extra_nf;
  /* pure monomials to be reduced */
  for (int32_t i = 0; i < count_not_lm;i++){
    lens_extra_nf[i]=1;
    int32_t j= extra_nf[i];
    for (int k = 0; k < nv-1; k++) {
      exps_extra_nf[i*nv+k]=lmb[j*nv+k];
    }
    exps_extra_nf[i*nv+nv-1]=lmb[j*nv+nv-1]+1;
    cfs_extra_nf[i]=1;
  }
  bs_t* tbr;
  if (count_not_lm) {
    md_t* md = copy_meta_data(st,fc);
    tbr = initialize_basis(md);
    exp_t *mul = (exp_t *)calloc(bs->ht->evl, sizeof(exp_t));
    tbr->ht = ht;
    /* reduction */
    import_input_data(tbr, md, 0, count_not_lm, lens_extra_nf, exps_extra_nf,
		      (void *)cfs_extra_nf, NULL);
    tbr->ld = tbr->lml  =  count_not_lm;
    for (long k = 0; k < count_not_lm; ++k) {
      tbr->lmps[k]  = k; /* fix input element in bs */
    }
    int32_t err = 0;
    tbr = core_nf(tbr, md, mul, bs, &err);
    if (err) {
      printf("Problem with normalform, stopped computation.\n");
      exit(1);
    }
    free(mul);
    free(md);
  }
  sp_matfglm_t *matrix ALIGNED32 = calloc(1, sizeof(sp_matfglm_t));
  matrix->charac = fc;
  matrix->ncols = dquot;
  matrix->nrows = len0;
  matrix->nnfs  = count_not_lm;
  long len1 = dquot * len0;
  long len2 = dquot - len0;

  if(posix_memalign((void **)&matrix->dense_mat, 32, sizeof(CF_t)*len1)){
    fprintf(stderr, "Problem when allocating matrix->dense_mat\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len1; i++){
      matrix->dense_mat[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_idx, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->triv_pos, 32, sizeof(CF_t)*len2)){
    fprintf(stderr, "Problem when allocating matrix->triv_pos\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len2; i++){
      matrix->triv_pos[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dense_idx, 32, sizeof(CF_t)*len0)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len0; i++){
      matrix->dense_idx[i] = 0;
    }
  }
  if(posix_memalign((void **)&matrix->dst, 32, sizeof(CF_t)*len0)){
    fprintf(stderr, "Problem when allocating matrix->dense_idx\n");
    exit(1);
  }
  else{
    for(long i = 0; i < len0; i++){
      matrix->dst[i] = 0;
    }
  }

  long l_triv = 0;
  long l_dens = 0;
  long nrows = 0;
  long count = 0;
  long count_nf = 0;

  for(long i = 0; i < dquot; i++){
    long pos = -1;
    int32_t *exp = lmb + (i * nv);
#if DEBUGBUILDMATRIX > 0
    display_monomial_full(stderr, nv, NULL, 0, exp);
    //    fprintf(stderr, "\n");
#endif
    if(member_xxn(exp, lmb + (i * nv), dquot - i, &pos, nv)){
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => remains in monomial basis\n");
#endif

      matrix->triv_idx[l_triv] = i;
      matrix->triv_pos[l_triv] = pos + i;

      l_triv++;
    }
    else{
#if DEBUGBUILDMATRIX > 0
      fprintf(stderr, " => does NOT remain in monomial basis\n");
#endif
      matrix->dense_idx[l_dens] = i;
      l_dens++;
      if(is_equal_exponent_xxn(exp, bexp_lm+(div_xn[count])*nv, nv)){
#if DEBUGBUILDMATRIX > 0
	fprintf(stderr, " => lands on a leading monomial\n");
#endif
        copy_poly_in_matrix_from_bs(matrix, nrows, bs, ht, //bcf, bexp, blen,
                                    div_xn[count], len_gb_xn[count],
                                    start_cf_gb_xn[count], len_gb_xn[count], lmb,
                                    nv, fc);
        nrows++;
        count++;
        if(len_xn < count && i < dquot){
	  if (info_level){
	    fprintf(stderr, "Staircase is not generic (1 => explain better)\n");
	  }
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

	  free_basis_without_hash_table(&tbr);
	  free(cfs_extra_nf);
	  free(exps_extra_nf);
	  free(lens_extra_nf);
	  free(extra_nf);
          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
	  free(evi);
          return NULL;
        }
      }
      else if (i == extra_nf[count_nf]){
#if DEBUGBUILDMATRIX > 0
	fprintf(stderr, " => lands on a MULTIPLE of a leading monomial\n");
#endif
	copy_nf_in_matrix_from_bs(matrix, nrows, count_nf, lmb,
				  tbr, ht, evi, st, nv);
	nrows++;
	count_nf++;
        if(count_not_lm < count_nf && i < dquot){
	  if (info_level){
	    fprintf(stderr, "Staircase is not generic (1 => explain better)\n");
	  }
          free(matrix->dense_mat);
          free(matrix->dense_idx);
          free(matrix->triv_idx);
          free(matrix->triv_pos);
          free(matrix->dst);
          free(matrix);

	  free_basis_without_hash_table(&tbr);
	  free(cfs_extra_nf);
	  free(exps_extra_nf);
	  free(lens_extra_nf);
	  free(extra_nf);
          free(len_gb_xn);
          free(start_cf_gb_xn);
          free(div_xn);
	  free(evi);
          return NULL;
        }
      }
      else{ /* should not arrive here */
	fprintf(stderr, "Staircase is not generic\n");
	fprintf(stderr, "Multiplication by ");
	display_monomial_full(stderr, nv, NULL, 0, exp);
	fprintf(stderr, " gets outside the staircase\n");
        free(matrix->dense_mat);
        free(matrix->dense_idx);
        free(matrix->triv_idx);
        free(matrix->triv_pos);
        free(matrix->dst);
        free(matrix);

	free_basis_without_hash_table(&tbr);
	free(cfs_extra_nf);
	free(exps_extra_nf);
	free(lens_extra_nf);
	free(extra_nf);
        free(len_gb_xn);
        free(start_cf_gb_xn);
        free(div_xn);
        return NULL;
      }
    }
  }
  for(long i = 0; i < matrix->nrows; i++){
    for(long j = matrix->ncols - 1; j >= 0; j--){
      if(matrix->dense_mat[i*matrix->ncols + j] == 0){
        matrix->dst[i]++;
      }
      else{
        break;
      }
    }
  }

  if (count_not_lm) {
    free_basis_without_hash_table(&tbr);
  }
  if(st->info_level){
    fprintf(stderr, "[%lu, %lu], Free / Dense = %.2f%%\n",
            len0, len_xn,
            100*((double)len_xn / (double)len0));
  }
  return matrix;
}



static inline int32_t *get_lm_from_bs(const bs_t * const bs, const ht_t * const ht){
  hm_t *dt;
  const len_t nelts = bs->lml;
  const int nv    = ht->nv;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;
  int32_t *exp  = (int32_t *)malloc(
                                    (unsigned long)(nelts) * (unsigned long)(nv) * sizeof(int32_t));
  /* counters for lengths, exponents and coefficients */
  int64_t cl = 0, ce = 0;//, cc = 0, ctmp  = 0;;

  for (long i = 0; i < nelts; ++i) {
    const bl_t bi = bs->lmps[i];

    dt  = bs->hm[bi] + OFFSET;
    for (int k = 1; k < ebl; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    for (int k = ebl+1; k < evl; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    cl++;
  }
  return exp;
}


static inline void get_lm_from_bs_trace(bs_t *bs, const ht_t *ht, int32_t *exp){
  hm_t *dt;
  const len_t nelts = bs->lml;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  /* counters for lengths, exponents and coefficients */
  int64_t cl = 0, ce = 0;//, cc = 0, ctmp  = 0;;

  for (long i = 0; i < nelts; ++i) {
    const bl_t bi = bs->lmps[i];
    //    len[cl] = bs->hm[bi][LENGTH];

    dt  = bs->hm[bi] + OFFSET;
    for (int k = 1; k < ebl; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    for (int k = ebl+1; k < evl; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    //    cc  +=  len[cl];
    cl++;
  }
}

static inline void get_lm_from_bs_trace_elim(bs_t *bs, const ht_t *ht, int32_t *exp,
                                             const len_t nelts){
  hm_t *dt;
  const len_t ebl = ht->ebl;
  const len_t evl = ht->evl;

  /* counters for lengths, exponents and coefficients */
  int64_t cl = 0, ce = 0;//, cc = 0, ctmp  = 0;;

  for (long i = 0; i < nelts; ++i) {
    const bl_t bi = bs->lmps[i];
    //    len[cl] = bs->hm[bi][LENGTH];

    dt  = bs->hm[bi] + OFFSET;
    /* for (int k = 1; k < ebl; ++k) { */
    /*   exp[ce++] = (int32_t)ht->ev[dt[0]][k]; */
    /* } */
    for (int k = ebl+1; k < evl; ++k) {
      exp[ce++] = (int32_t)ht->ev[dt[0]][k];
    }
    //    cc  +=  len[cl];
    cl++;
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


#undef REDUCTION_ALLINONE

