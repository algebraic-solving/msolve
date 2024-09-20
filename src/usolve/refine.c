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

#ifndef USOLVE
#define USOLVE
#endif

/* based on Quadratic interval refinement for real roots */
/* J. Abott, ACM communications in computer algebra */

static void make_exact_root(mpz_t *upol, unsigned long int *deg,
                            interval *rt, mpz_t *x, int k){
  mpz_set(rt->numer, *x);
  rt->k = k;
  rt->isexact=1;
  if(k>=0){
    USOLVEnumer_quotient(upol, deg, rt->numer, k);
  }
  else{
    USOLVEnumer_quotient(upol, deg, rt->numer, 0);
  }
}

/* computes 2^logN f(a) / (f(a) - f(b)) where f(a) = vala, f(b) = valb */
static long long int index_linearinterp(mpz_t *vala, mpz_t *valb, mpz_t *q,
                                        long long int logN){
  mpz_sub(*valb, *vala, *valb);
  mpz_mul_2exp(*vala, *vala, logN);
  mpz_tdiv_q(*q, *vala, *valb);

  long int sizeq = ilog2_mpz(*q);
  if(sizeq >= 8 * sizeof(long long int)){
    if(sizeq > logN){
      fprintf(stderr,"Valeur de q = "); mpz_out_str(stderr, 10, *q);puts("");
      fprintf(stderr, "Valeur de Nlog = %lld\n", logN);
      fprintf(stderr, "ilog2(q) = %ld\n", sizeq);
      return -2;
    }
    return -1;
  }
  long long int index = mpz_get_ui(*q);
  return index;
}

static void getx_and_eval_2exp(mpz_t *upol, unsigned long int deg,
                               mpz_t *a, mpz_t *x, mpz_t *value, mpz_t *q,
                               int k, long int index){
  mpz_set(*x, *a);
  mpz_add_ui(*x, *x, index);
  mpz_poly_eval_2exp_naive(upol, deg, x, k, value, q);
}

static void getx_and_eval_2exp_mpzidx(mpz_t *upol, unsigned long int deg,
                                      mpz_t *a, mpz_t *x, mpz_t *value, mpz_t *q,
                                      int k, mpz_t* index_gmp){
  mpz_set(*x, *a);
  mpz_add(*x, *x, *index_gmp);
  mpz_poly_eval_2exp_naive(upol, deg, x, k, value, q);
}

static void getx_and_eval_2expleft(mpz_t *upol, unsigned long int deg,
                                   mpz_t *a, mpz_t *x, mpz_t *value, mpz_t *q,
                                   int k, long int index){
  mpz_set(*x, *a);
  mpz_sub_ui(*x, *x, index);
  mpz_poly_eval_2exp_naive(upol, deg, x, k, value, q);
}

/* assumes k >= -2 */
static void getx_and_eval(mpz_t *upol, unsigned long int deg,
                          mpz_t *a, mpz_t *x, mpz_t *value, mpz_t *q,
                          long int Nlog, long int k, long int index){

  if(k-Nlog>=0){
    mpz_set_ui(*x, 1);
    mpz_mul_2exp(*x, *x, k - Nlog);
    mpz_mul_ui(*x,*x,index);
    mpz_add(*x, *x, *a);

    mpz_poly_eval_2exp_naive(upol, deg, x, 0, value, q);
  }
  else{
    mpz_set(*x, *a);
    mpz_mul_2exp(*x, *x, Nlog - k);
    mpz_add_ui(*x, *x, index);

    mpz_poly_eval_2exp_naive(upol, deg, x, Nlog - k, value, q);
  }
}

static void getx_and_eval_mpzidx(mpz_t *upol, unsigned long int deg,
                                 mpz_t *a, mpz_t *x, mpz_t *value, mpz_t *q,
                                 long int Nlog, long int k, mpz_t *index){
  if(k-Nlog >= 0){
    mpz_set_ui(*x, 1);
    mpz_mul_2exp(*x, *x, k - Nlog);
    mpz_mul(*x, *x, *index);
    mpz_add(*x, *x, *a);

    mpz_poly_eval_2exp_naive(upol, deg, x, 0, value, q);
  }
  else{
    mpz_set(*x, *a);
    mpz_mul_2exp(*x, *x, Nlog - k);
    mpz_add(*x, *x, *index);

    mpz_poly_eval_2exp_naive(upol, deg, x, Nlog - k, value, q);

  }
}

static int right_interval_2exp(mpz_t *upol, unsigned long int *deg_ptr,
                               interval *rt, mpz_t *x, mpz_t *b,
                               mpz_t *vala, mpz_t *valb, mpz_t *q,
                               int k, int newk){

  mpz_swap(*x, *b);
  mpz_swap(*vala, *valb);

  getx_and_eval_2exp(upol, *deg_ptr, x, b, valb, q, newk, 1);

  int sgnb = mpz_sgn(*valb);
  int sgnx = mpz_sgn(*vala);

  if(sgnb==0){
    make_exact_root(upol, deg_ptr, rt, b, newk);
    return 1;
  }
  if(sgnx!=sgnb){
    mpz_set(rt->numer, *x);
    rt->k = newk;
    return 1;
  }
  return 0;
}

static int left_interval_2exp(mpz_t *upol, unsigned long int *deg_ptr,
                                   interval *rt, mpz_t *x, mpz_t *b,
                                   mpz_t *vala, mpz_t *valb, mpz_t *q,
                                   int k, int newk){

  getx_and_eval_2expleft(upol, *deg_ptr, x, b, valb, q,
                         newk, 1);
  mpz_swap(*x, *b);
  mpz_swap(*vala, *valb);
  int sgnb = mpz_sgn(*valb);
  int sgnx = mpz_sgn(*vala);

  if(sgnx==0){
    make_exact_root(upol, deg_ptr, rt, x, newk);
    return 1;
  }
  if(sgnx!=sgnb){
    mpz_set(rt->numer, *x);
    rt->k = newk;
    return 1;
  }
  return 0;
}

static int right_interval(mpz_t *upol, unsigned long int *deg_ptr,
                               interval *rt,
                               mpz_t *x, mpz_t *b,
                               mpz_t *vala, mpz_t *valb, mpz_t *q,
                               int Nlog, int k, int newk){

  mpz_swap(*x, *b);
  mpz_swap(*vala, *valb);

  if(k-Nlog>=0){
    mpz_set_ui(*b, 1);
    mpz_mul_2exp(*b, *b, k-Nlog);
    mpz_add(*b, *x, *b);
    mpz_poly_eval_2exp_naive(upol, *deg_ptr,
                             b, 0, valb, q);
  }
  else{
    mpz_add_ui(*b, *x, 1);
    mpz_poly_eval_2exp_naive(upol, *deg_ptr,
                             b, Nlog - k, valb, q);
  }

  int sgnb = mpz_sgn(*valb);
  int sgnx = mpz_sgn(*vala);

  if(sgnb==0){
    make_exact_root(upol, deg_ptr, rt, b, newk);
    return 1;
  }
  if(sgnx!=sgnb){
    mpz_set(rt->numer, *x);
    if(k-Nlog>=0){
      rt->k = newk;
    }
    else{
      rt->k = Nlog - k;
    }
    return 1;
  }
  return 0;
}

static int left_interval(mpz_t *upoly, unsigned long int *deg_ptr,
                              interval *rt,
                              mpz_t *x, mpz_t *b,
                              mpz_t *vala, mpz_t *valb,
                              mpz_t *q, int Nlog, int k, int newk){

  if(k-Nlog>=0){
    mpz_set_ui(*b, 1);
    mpz_mul_2exp(*b, *b, k-Nlog);
    mpz_sub(*b, *x, *b);
    mpz_poly_eval_2exp_naive(upoly, *deg_ptr,
                             b, 0, valb, q);
  }
  else{
    mpz_sub_ui(*b, *x, 1);
    mpz_poly_eval_2exp_naive(upoly, *deg_ptr,
                             b, Nlog - k, valb, q);
  }

  mpz_swap(*x, *b);
  mpz_swap(*vala, *valb);

  int sgnb = mpz_sgn(*valb);
  int sgnx = mpz_sgn(*vala);
  if(sgnx==0){
    make_exact_root(upoly, deg_ptr, rt, x, newk);
    return 1;
  }
  if(sgnx!=sgnb){
    mpz_set(rt->numer, *x);
    if(k-Nlog>=0){
      rt->k = newk;
    }
    else{
      rt->k = Nlog - k;
    }
    return 1;
  }
  return 0;
}


static void refine_root_by_N_positive_k(mpz_t *upol, unsigned long int *deg_ptr,
                                        interval *rt,
                                        mpz_t *tab, int64_t Nlog,
                                        int *success, int verbose){
  int64_t k = (rt-> k);
  int64_t newk = k+Nlog;
  int sgna, sgnb, sgnx;
  int64_t index;

  mpz_t *vala = (tab);
  mpz_t *valb = (tab+1);
  mpz_t *b = (tab+2);
  mpz_t *q = (tab+3);
  mpz_t *x = (tab+4);
  mpz_t *tmpvala = (tab+5);
  mpz_t *tmpvalb = (tab+6);
  mpz_t *index_gmp = (tab+7);
  mpz_set(*tmpvala, *vala);
  mpz_set(*tmpvalb, *valb);

  *success = 1;
  sgna =  mpz_sgn(*vala);
  sgnb =  mpz_sgn(*valb);
  /* computes 2^logN f(a) / (f(a) - f(b)) where f(a) = vala, f(b) = valb */
  index = index_linearinterp(vala, valb, index_gmp, Nlog);

  mpz_set(*vala, *tmpvala);
  mpz_set(*valb, *tmpvalb);
  int64_t maxindex = (1L<<(Nlog));

  if(index == -2 || index == 0 || (LOG2(index) > Nlog && index > 0) ){
    if(Nlog == 2) index = 2;
    else{
      *success = 0;
      mpz_set(*vala, *tmpvala);
      mpz_set(*valb, *tmpvalb);

      return;
    }
  }

  /* One takes x = (2^Nlog * a + index) / 2^(k+Nlog) */
  /* We will compute f(x) and store the value in vala */
  mpz_mul_2exp(rt->numer, rt->numer, Nlog);

  if(index >= 0){
    getx_and_eval_2exp(upol, *deg_ptr, &(rt->numer),x, vala, q,
                       newk, index);
  }
  else{
    getx_and_eval_2exp_mpzidx(upol, *deg_ptr, &(rt->numer), x, vala, q,
                              newk, index_gmp);
  }

  sgnx = mpz_sgn(*vala);

  if(sgnx==0){
    make_exact_root(upol, deg_ptr, rt, x, newk);
    return;
  }

  if(sgna == sgnx){
    if(index == 3 && Nlog == 2){
      mpz_set(rt->numer, *x);
      rt->k = newk;
      return;
    }

    getx_and_eval_2exp(upol, *deg_ptr, x, b, valb, q, newk, 1);

    sgnb = mpz_sgn(*valb);

    if(sgnb==0){
      make_exact_root(upol, deg_ptr, rt, b, newk);
      return;
    }
    if(sgnx != sgnb){
      mpz_set(rt->numer, *x);
      rt->k = newk;
      return;
    }
    else{
      *success = 0;
      if(Nlog!=2){
        mpz_set(*vala, *tmpvala);
        mpz_set(*valb, *tmpvalb);
        mpz_tdiv_q_2exp(rt->numer, rt->numer, Nlog);
        return;
      }
      while(right_interval_2exp(upol, deg_ptr, rt, x, b,
                                vala, valb, q, k, newk)==0 && index <= maxindex){
        if(verbose>0){
          fprintf(stderr, "->");
        }
        index++;
      }
    }
  }

  if(sgna != sgnx){
    if(index == 1){
      rt->k = newk;
      mpz_swap(*vala, *valb);
      mpz_set(*vala, *tmpvala);
      return;
    }
    /* One computes x - 2^(k-2) */
    getx_and_eval_2expleft(upol, *deg_ptr, x, b, valb, q,
                           newk, 1);
    /* x is now the least element */
    mpz_swap(*x, *b);
    mpz_swap(*vala, *valb);
    sgnb = mpz_sgn(*valb);
    sgnx = mpz_sgn(*vala);
    if(sgnx==0){
      make_exact_root(upol, deg_ptr, rt, x, newk);
      return;
    }
    if(sgnx != sgnb){
      mpz_set(rt->numer, *x);
      rt->k = newk;
      return;
    }
    else{
      *success = 0;
      if(Nlog != 2){
        mpz_set(*vala, *tmpvala);
        mpz_set(*valb, *tmpvalb);
        mpz_tdiv_q_2exp(rt->numer, rt->numer, Nlog);
        return;
      }
      while(left_interval_2exp(upol, deg_ptr, rt,
                               x, b, vala, valb, q, k, newk) == 0){
        if(verbose>0){
          fprintf(stderr, "<-");
        }
      }
    }
  }
}

//Dans ce cas, l'intervalle est code de la maniere suivante :
//On note u  = root->numer et k = root-> k
//On considere alors l'intervalle (a, b) avec
//a = u et b = u+2^(-k) et k < 0 (donc 2^(-k) est un entier)
//tab est un tableau de 7 entiers GMP
//tab[0] et tab[1] contiennent f(a) et f(b) respectivement
//tab[2] contient b
static void refine_root_by_N_negative_k(mpz_t *upol, unsigned long int *deg_ptr,
                                        interval *rt, mpz_t *tab, int64_t Nlog,
                                        int *success, int verbose){
  long newk;
  /* one takes the opposite */
  long k = - (rt-> k);
  int sgna, sgnb, sgnx;
  long long int index;

  mpz_t *vala = (tab);
  mpz_t *valb = (tab+1);
  mpz_t *b = (tab+2);
  mpz_t *q = (tab+3);
  mpz_t *x = (tab+4);
  mpz_t *tmpvala = (tab+5);
  mpz_t *tmpvalb = (tab+6);
  mpz_t *index_gmp = (tab+7);
  mpz_set(*tmpvala, *vala);
  mpz_set(*tmpvalb, *valb);
  *success = 1;
  sgna =  mpz_sgn(*vala);
  sgnb =  mpz_sgn(*valb);
  index = index_linearinterp(vala, valb, index_gmp, Nlog);

  mpz_set(*vala, *tmpvala);
  mpz_set(*valb, *tmpvalb);
  newk = -k + Nlog;

  if(index==-2 || index ==0 || (index>0 && LOG2(index) > Nlog)){
    if(Nlog == 2)index = 2;
    else{
      *success = 0;
      mpz_set(*vala, *tmpvala);
      mpz_set(*valb, *tmpvalb);
      return;
    }
  }

  if(index >= 0){
    getx_and_eval(upol, *deg_ptr, &(rt->numer), x, vala, q,
                  Nlog, k, index);
  }
  else{
    getx_and_eval_mpzidx(upol, *deg_ptr, &(rt->numer),x, vala, q,
                         Nlog, k, index_gmp);
  }

  sgnx = mpz_sgn(*vala);

  if(sgnx==0){
    make_exact_root(upol, deg_ptr, rt, x, newk);
    return;
  }

  if(sgna == sgnx){
    if(index == 3 && Nlog == 2){
      mpz_set(rt->numer, *x);
      rt->k = newk;
      return;
    }
    if(k-Nlog>=0){
      mpz_set_ui(*b, 1);
      mpz_mul_2exp(*b, *b, k-Nlog);
      mpz_add(*b, *x, *b);

      mpz_poly_eval_2exp_naive(upol, *deg_ptr, b, 0, valb, q);
    }
    else{
      mpz_add_ui(*b, *x, 1);
      mpz_poly_eval_2exp_naive(upol, *deg_ptr, b, Nlog - k, valb, q);
    }
    sgnb = mpz_sgn(*valb);
    sgnx = mpz_sgn(*tmpvala);

    if(sgnb==0){
      make_exact_root(upol, deg_ptr, rt, b, newk);
      return;
    }
    if(sgnx != sgnb){
      mpz_set(rt->numer, *x);
      if(k-Nlog>=0){
        rt->k = newk;
      }
      else{
        rt->k = Nlog - k;
      }
      return;
    }
    else{

      *success = 0;
      if(Nlog!=2){
        mpz_set(*vala, *tmpvala);
        mpz_set(*valb, *tmpvalb);
        return;
      }

      while(right_interval(upol, deg_ptr, rt, x, b,
                           vala, valb, q, Nlog, k, newk) == 0){
        if(verbose>0){
          fprintf(stderr, "|->");
        }
      }
    }
  }

  if(sgna != sgnx){
    if(index == 1){
      if(newk > 0){
        mpz_mul_2exp(rt->numer, rt->numer, newk);
      }
      rt->k = newk;
      mpz_swap(*vala, *valb);
      mpz_set(*vala, *tmpvala);
      return;
    }

    if(k-Nlog >= 0){
      mpz_set_ui(*b, 1);
      mpz_mul_2exp(*b, *b, k-Nlog);
      mpz_sub(*b, *x, *b);

      mpz_poly_eval_2exp_naive(upol, *deg_ptr, b, 0, valb, q);
    }
    else{
      mpz_sub_ui(*b, *x, 1);
      mpz_poly_eval_2exp_naive(upol, *deg_ptr, b, Nlog - k, valb, q);
    }

    mpz_swap(*x, *b);
    mpz_swap(*vala, *valb);

    sgnb = mpz_sgn(*valb);
    sgnx = mpz_sgn(*vala);
    if(sgnx==0){
      make_exact_root(upol, deg_ptr, rt, x, newk);
      return;
    }
    if(sgnx != sgnb){
      mpz_set(rt->numer, *x);
      if(k-Nlog>=0){
        rt->k = newk;
      }
      else{
        rt->k = Nlog - k;
      }
      return;
    }
    else{
      *success = 0;
      if(Nlog!=2){
        mpz_set(*vala, *tmpvala);
        mpz_set(*valb, *tmpvalb);
        return;
      }
      while(left_interval(upol, deg_ptr, rt, x, b,
                          vala, valb, q, Nlog, k, newk) == 0){
        if(verbose>0){
          fprintf(stderr, "<-|");
        }
      }
    }
  }
}


/*
  Let f = upol and (a, b) be an isolating interval
  Computes ( 4 f(a) / (f(a) - f(b)) ) and the closest integer k 

  One has w = (b-a) / 4 and a_i = a_i + i x w

  Set x = a_k.

  if f(x) = 0, root is x
  if sgn(f(x))=sgn(f(a)) :
       if f(a_{k+1}) =0 root is a_{k+1}
       if f(x) x f(a_{k+1}) < 0, one gets (x, a_{k+1})
       else success = false

  if sgn(f(x)) = -sgn(f(a))
       if f(a_{k-1}) =0 root is a_{k-1}
       if f(a_{k-1}) x f(x) <0 one gets (a_{k-1}, x)
       else success = false

 */
static void refine_positive_root_by_N(mpz_t *upol, unsigned long int *deg_ptr,
                                      interval *rt, mpz_t *tab,
                                      int64_t Nlog, int *success,
                                      int verbose){
  *success = 1;

  long k = rt-> k;

  if(rt->isexact == 1){
    return;
  }

  if(k < 0){

    /* interval is (rt->numer, rt->numer + 2^(-k)) */
    refine_root_by_N_negative_k(upol, deg_ptr, rt,
                                tab, Nlog, success, verbose);
  }
  else{

    /* interval is (rt->numer/2^k, rt->numer/2^k) */
    refine_root_by_N_positive_k(upol, deg_ptr, rt,
                                tab, Nlog, success, verbose);
  }
}

void refine_QIR_positive_root(mpz_t *upol, unsigned long int *deg_ptr,
                              interval *rt, mpz_t *tab, int prec, int verbose){
  if(rt->isexact==1) return;
  int64_t Nlog = 2;
  int success = 1;
  while(rt->isexact != 1 && rt->k < prec){

    refine_positive_root_by_N(upol, deg_ptr, rt, tab,
                              Nlog, &success, verbose);
    if(rt->isexact==1) return;

    if(mpz_sgn(tab[0]) == mpz_sgn(tab[1])) {
      fprintf(stderr, "BUG in refine_QIR_positive_root");
      exit(1);
    }
    if(success == 1){
      if(prec - rt->k > Nlog)
      Nlog=2*Nlog;
    }
    else{
      if(Nlog>2){
        Nlog = Nlog / 2;
      }
    }
  }
}


void get_values_at_bounds(mpz_t *upol, unsigned long int deg,
                          interval *rt, mpz_t *tab){
  if(rt->k > 0){
    mpz_poly_eval_2exp_naive(upol, deg, &rt->numer, rt->k, tab, tab+5);
    mpz_set(tab[3], rt->numer);
    mpz_add_ui(tab[3], tab[3], 1);
    mpz_poly_eval_2exp_naive(upol, deg, tab+3, rt->k, tab+1, tab+5);
  }
  else{
    basic_mpz_poly_eval_at_point(upol, deg, &(rt->numer), tab);
    mpz_set_ui(tab[3],1);
    mpz_mul_2exp(tab[3], tab[3], - (rt->k));
    mpz_add(tab[3], tab[3], rt->numer);
    basic_mpz_poly_eval_at_point(upol, deg, tab+3, tab+1);
  }
}

static void refine_root_naive(mpz_t *upol, unsigned long int deg,
                              interval *rt, mpz_t *middle, int calgo){
  if(rt->isexact == 1){
    return;
  }

  long newk;
  int sgn_middle; 
  if((rt->k) >= 0){
    /* interval is (numer, numer + 2^(-k)) */
    mpz_mul_ui(*middle, rt->numer, 2);
    mpz_add_ui(*middle, *middle, 1);
    newk = rt->k + 1;

    sgn_middle = sgn_mpz_poly_eval_at_point_2exp_naive(upol, deg, middle, newk);
  }
  else{
    mpz_set_ui(*middle, 1);
    mpz_mul_2exp(*middle, *middle, - ( (rt->k) + 1));
    mpz_add(*middle, *middle, rt->numer);
    newk = ( (rt->k) + 1 );

    sgn_middle = sgn_mpz_poly_eval_at_point_naive(upol, deg, middle, calgo);
  }
  int sign_left = rt->sign_left;

  if(sgn_middle * sign_left<0){
    if(newk>0){
      mpz_mul_ui(rt->numer, rt->numer, 2);
    }
    rt->k = newk;
  }
  else{
    mpz_set(rt->numer, *middle);
    rt->k = newk;
  }
}


static void remove_exact_roots_by_division(mpz_t *upol, unsigned long int *deg,
                                           interval *roots,
                                           unsigned long int nbroots, int nthreads){
  for(unsigned long int i = 0; i < nbroots; i++){
    interval *rt = roots + i;
    if(rt->isexact==1){
      if(rt->k < 0){
        USOLVEnumer_quotient(upol, deg, rt->numer, 0);
      }
      else{
        USOLVEnumer_quotient(upol, deg, rt->numer, rt->k);
      }
    }
  }
}


/* Refinement using Newton-Interval like technique (but replacing Newton with */
/* linear interpolation) */
/* it takes as input a pointer to deg because it may change after performing */
/* divisions when there are exact roots */
void refine_QIR_roots(mpz_t *upol, unsigned long int *deg, interval *roots,
                      int nbneg, int nbpos,
                      int prec, int verbose, double step, int nthreads){
  unsigned long int i;
  /* table for intermediate values */
  mpz_t *tab = (mpz_t *)(malloc(sizeof(mpz_t) * 8));
  for(i=0;i<8;i++){
    mpz_init(tab[i]);
  }

  double e_time = 0, refine_time = realtime();
  int nb = nbneg + nbpos;

  remove_exact_roots_by_division(upol, deg, roots, nb, nthreads);


  interval *pos_rt = (interval *)(malloc(sizeof(interval)));
  mpz_init(pos_rt->numer);
  mpz_t newc;
  mpz_init(newc);

  for(i = 0; i <= *deg; i++){
    if(i%2 == 1){
      mpz_neg(upol[i], upol[i]);
    }
  }

  for(i = 0; i < nbneg; i++){

    interval *rt = roots + i;

    /* display_root(stderr, rt); */

    if(rt->k > 0){
      if(rt->isexact!=1){
        mpz_add_ui(pos_rt->numer, rt->numer, 1);
        mpz_neg(pos_rt->numer, pos_rt->numer);
      }
      pos_rt->k = rt->k;
      pos_rt->sign_left = - (rt->sign_left);
      pos_rt->isexact = rt->isexact;
    }
    else {
      if(rt->isexact!=1){
        mpz_set_ui(newc, 1);
        mpz_mul_2exp(newc, newc, -rt->k);
        mpz_add(pos_rt->numer, rt->numer, newc);
        mpz_neg(pos_rt->numer, pos_rt->numer);
      }
      pos_rt->k = rt->k;
      pos_rt->sign_left = - (rt->sign_left);
      pos_rt->isexact = rt->isexact;
    }

    if(pos_rt->isexact==0){
      get_values_at_bounds(upol, *deg, pos_rt, tab);
      if(mpz_sgn(tab[0])==0 || mpz_sgn(tab[1])==0){
        fprintf(stderr, "Error in refinement (neg. roots): these values should not be zero\n");
        exit(1);
      }
      refine_QIR_positive_root(upol, deg, pos_rt, tab, prec, verbose);

      if(mpz_sgn(tab[0])==mpz_sgn(tab[1])){
        fprintf(stderr, "BUG in refinement (sgn tab[0]==sgn tab[1]) for neg. roots");
        exit(1);
      }
    }

    if(pos_rt->isexact==1){
      if(pos_rt->k < 0){
        pos_rt->k = 0;
      }
    }
    //We assume precision >=0
    if(pos_rt->isexact!=1){
      rt->k = pos_rt->k;
      rt->isexact = pos_rt->isexact;
      mpz_add_ui(rt->numer, pos_rt->numer, 1);
      mpz_neg(rt->numer, rt->numer);
    }
    else{
      rt->k = pos_rt->k;
      if(rt->isexact!=1){
        rt->isexact = pos_rt->isexact;
        mpz_set(rt->numer, pos_rt->numer);
        mpz_neg(rt->numer, rt->numer);
      }
    }

    e_time += realtime() - refine_time;
    if(e_time>=step){
      refine_time = realtime();
      e_time = 0;
      if(verbose>=1){
        fprintf(stderr, "{%.2f%s}", ((double)i / nb) * 100, "%");
      }
    }
  }

  mpz_clear(newc);
  mpz_clear(pos_rt->numer);
  free(pos_rt);

  /* now we refine positive roots */
  for(i = 0; i <= *deg; i++){
    if(i%2 == 1){
      mpz_neg(upol[i], upol[i]);
    }
  }

  for(i=nbneg; i < nb; i++){
    interval *rt = roots + i;

    if(rt->isexact==0){
      get_values_at_bounds(upol, *deg, rt, tab);
      if(mpz_sgn(tab[1])==0 || mpz_sgn(tab[0])==0){
        fprintf(stderr, "Error in refinement (pos. roots): these values should not be zero\n");
        exit(1);
      }
      refine_QIR_positive_root(upol, deg, rt, tab, prec, verbose);
      if(mpz_sgn(tab[0])==mpz_sgn(tab[1])){
        fprintf(stderr,"BUG in refinement (sgn tab[0]=sgn tab[1] for pos. roots)");
        exit(1);
      }
      if(rt->isexact==1){
        if(rt->k < 0){
          rt->k = 0;
        }
      }
    }

    e_time += realtime() - refine_time;
    if(e_time>=step){
      refine_time = realtime();
      e_time = 0;
      if(verbose>=1){
        fprintf(stderr, "{%.2f%s}", ((double)(i) / nb) * 100, "%");
      }
    }

  }
  if(verbose>=1){
    fprintf(stderr, "\n");
  }
  for(i = 0; i < 8; i++){
    mpz_clear(tab[i]);
  }
  free(tab);
}

/* Refinement using Newton-Interval like technique (but replacing Newton with */
/* linear interpolation) */
/* it takes as input a pointer to deg because it may change after performing */
/* divisions when there are exact roots */
/* the precision of the refinement depends on the value of the root to be defined */
void refine_QIR_roots_adaptative(mpz_t *upol, unsigned long int *deg, interval *roots,
                                 int nbneg, int nbpos,
                                 int prec, int verbose, double step, int nthreads){

  unsigned long int i;
  /* table for intermediate values */
  mpz_t *tab = (mpz_t *)(malloc(sizeof(mpz_t) * 8));
  for(i=0;i<8;i++){
    mpz_init(tab[i]);
  }

  double e_time = 0, refine_time = realtime();
  int nb = nbneg + nbpos;

  remove_exact_roots_by_division(upol, deg, roots, nb, nthreads);


  interval *pos_rt = (interval *)(malloc(sizeof(interval)));
  mpz_init(pos_rt->numer);
  mpz_t newc;
  mpz_init(newc);

  for(i = 0; i <= *deg; i++){
    if(i%2 == 1){
      mpz_neg(upol[i], upol[i]);
    }
  }

  for(i = 0; i < nbneg; i++){

    interval *rt = roots + i;

    /* display_root(stderr, rt); */

    if(rt->k > 0){
      if(rt->isexact!=1){
        mpz_add_ui(pos_rt->numer, rt->numer, 1);
        mpz_neg(pos_rt->numer, pos_rt->numer);
      }
      pos_rt->k = rt->k;
      pos_rt->sign_left = - (rt->sign_left);
      pos_rt->isexact = rt->isexact;
    }
    else {
      if(rt->isexact!=1){
        mpz_set_ui(newc, 1);
        mpz_mul_2exp(newc, newc, -rt->k);
        mpz_add(pos_rt->numer, rt->numer, newc);
        mpz_neg(pos_rt->numer, pos_rt->numer);
      }
      pos_rt->k = rt->k;
      pos_rt->sign_left = - (rt->sign_left);
      pos_rt->isexact = rt->isexact;
    }

    if(pos_rt->isexact==0){
      get_values_at_bounds(upol, *deg, pos_rt, tab);
      if(mpz_sgn(tab[0])==0 || mpz_sgn(tab[1])==0){
        fprintf(stderr, "Error in refinement (neg. roots): these values should not be zero\n");
        exit(1);
      }
      long d = 1 + ilog2_mpz(pos_rt->numer) - rt->k;

      /* fprintf(stderr, "[%d, %ld]", prec, */
      /*         prec + ((*deg) * MAX(0, d)) / 32); */
      refine_QIR_positive_root(upol, deg, pos_rt, tab,
                               prec + (((*deg)-1) * MAX(0, d)) / 32, verbose);

      if(mpz_sgn(tab[0])==mpz_sgn(tab[1])){
        fprintf(stderr, "BUG in refinement (sgn tab[0]==sgn tab[1]) for neg. roots");
        exit(1);
      }
    }

    if(pos_rt->isexact==1){
      if(pos_rt->k < 0){
        pos_rt->k = 0;
      }
    }
    //We assume precision >=0
    if(pos_rt->isexact!=1){
      rt->k = pos_rt->k;
      rt->isexact = pos_rt->isexact;
      mpz_add_ui(rt->numer, pos_rt->numer, 1);
      mpz_neg(rt->numer, rt->numer);
    }
    else{
      rt->k = pos_rt->k;
      if(rt->isexact!=1){
        rt->isexact = pos_rt->isexact;
        mpz_set(rt->numer, pos_rt->numer);
        mpz_neg(rt->numer, rt->numer);
      }
    }

    e_time += realtime() - refine_time;
    if(e_time>=step){
      refine_time = realtime();
      e_time = 0;
      if(verbose>=1){
        fprintf(stderr, "{%.2f%s}", ((double)i / nb) * 100, "%");
      }
    }
  }

  mpz_clear(newc);
  mpz_clear(pos_rt->numer);
  free(pos_rt);

  /* now we refine positive roots */
  for(i = 0; i <= *deg; i++){
    if(i%2 == 1){
      mpz_neg(upol[i], upol[i]);
    }
  }

  for(i=nbneg; i < nb; i++){
    interval *rt = roots + i;

    if(rt->isexact==0){
      get_values_at_bounds(upol, *deg, rt, tab);
      if(mpz_sgn(tab[1])==0 || mpz_sgn(tab[0])==0){
        fprintf(stderr, "Error in refinement (pos. roots): these values should not be zero\n");
        exit(1);
      }
      long d = 1 + ilog2_mpz(rt->numer) - rt->k;

      /* fprintf(stderr, "[%d, %ld]", prec, prec + ((*deg) * MAX(0, 1 + d)) / 32); */
      refine_QIR_positive_root(upol, deg, rt, tab,
                               prec + (((*deg) - 1) * MAX(0, 1 + d)) / 32, verbose);
      if(mpz_sgn(tab[0])==mpz_sgn(tab[1])){
        fprintf(stderr,"BUG in refinement (sgn tab[0]=sgn tab[1] for pos. roots)");
        exit(1);
      }
      if(rt->isexact==1){
        if(rt->k < 0){
          rt->k = 0;
        }
      }
    }

    e_time += realtime() - refine_time;
    if(e_time>=step){
      refine_time = realtime();
      e_time = 0;
      if(verbose>=1){
        fprintf(stderr, "{%.2f%s}", ((double)(i) / nb) * 100, "%");
      }
    }

  }
  if(verbose>=1){
    fprintf(stderr, "\n");
  }
  for(i = 0; i < 8; i++){
    mpz_clear(tab[i]);
  }
  free(tab);
}


void refine_all_roots_naive(mpz_t *upol, unsigned long int deg,
                            interval *roots, unsigned long int nb,
                            unsigned int prec, int calgo, int debug){
  mpz_t *middle=malloc(sizeof(mpz_t));
  mpz_init(middle[0]);

  for(unsigned long int j = 0; j < nb; j++){
    while((roots+j)->k < prec && (roots+j)->isexact == 0){
      refine_root_naive(upol, deg, roots+j, middle, calgo);
    }
  }
  mpz_clear(middle[0]);
  free(middle);
}



