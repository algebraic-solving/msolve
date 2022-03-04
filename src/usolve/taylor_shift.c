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

#include "univmultiply.c"

#ifndef USOLVE
#define USOLVE
#endif



static void mpz_poly_add(mpz_t *res,
                         mpz_t *upol1, const unsigned long int deg1,
                         mpz_t *upol2, const unsigned long int deg2){
  if(deg1 > deg2){
    mpz_poly_add(res, upol2, deg2, upol1, deg1);
    return;
  }
  for(unsigned long int i = 0; i <= deg1; i++){
    mpz_add(res[i], upol1[i], upol2[i]);
  }
  for(unsigned long int i = deg1 + 1;i <= deg2; i++){
    mpz_set(res[i],upol2[i]);
  }
}

/*naive taylor shift by 1 (in place)*/
static void taylorshift1_naive(mpz_t *upol, const unsigned long int deg){
  long int i,j;
  for (i = 0; i <= deg-1; i++){
    for (j = deg-1 ; j >= i; j--){
      mpz_add(upol[j], upol[j], upol[j+1]);
    }
  }
}


static inline void mpz_poly_swap_th(mpz_t *res, mpz_t *pol,
                                    const unsigned long int deg,
                                    const unsigned int nthreads){
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif

#pragma omp parallel for num_threads(nthreads)
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_swap(res[i], pol[i]);
  }
}

static inline void mpz_poly_add_th(mpz_t *res, mpz_t *upol1,
                                   const unsigned long int deg1,
                                   mpz_t *upol2, const unsigned long int deg2,
                                   const unsigned int nthreads){
  if(deg1 > deg2){
    mpz_poly_add(res, upol2, deg2, upol1, deg1);
    return;
  }
  unsigned long int i;
#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads) //schedule(static)
  for(i = 0; i <= deg1; i++){
    mpz_add(res[i], upol1[i], upol2[i]);
  }
#pragma omp parallel for num_threads(nthreads) //schedule(static)
  for(i = deg1 + 1; i <= deg2; i++){
    mpz_set(res[i],upol2[i]);
  }
}

static void rescale_upoly_2exp(mpz_t *upol, const unsigned long int deg,
                               const unsigned long int content){
  if(content==0){
    return;
  }
  for(unsigned long int i = 0; i <= deg; i++){
    mpz_mul_2exp(upol[i], upol[i], content);
  }
}

//pol : polynome
//deg : degre(pol)
//tmpol : polynome de degre deg -- sert de tmpol
//shift_pwx : tableau de polynome
//le premier est de degre size, le second 2*size, le troisieme 4*size etc.
//le dernier est de degre 2^l*size avec l le plus grand possible tel que
//2^(l*size) est inferieur a deg. 
//nblocks doit etre la plus petite puissance de 2 qui divise deg/size.

//on a pol = pol1+X^size*pol2+X^(2*size)*pol3+X^(3*size)*pol4+...

//pol0, pol1, etc. sont de degre size-1 sauf le dernier qui est de degre
//final_degree

//Ainsi on a : pol = [pol_1, pol_2, ..., pol_k]
//deg(pol_i)=size-1 pour i<k
//deg(pol_k)=final_degree
//Tout a ete fait pour que k soit une puissance de 2 => c'est nblocks.

//On fait le Taylor shift sur chaque bloc pol0, pol1, etc.

//On note TL la fonction TaylorShift
//On a maintenant
//pol = [TL(pol_1), TL(pol_2), ..., TL(pol_k)]
//ce sont les nouveaux pol_1, pol_2, etc.

//On calcule (X+1)^size*pol_2, (X+1)^size*pol_2, (X+1)^size*pol_4, ...
//tous les blocs pairs

//ce sont des polynomes de degre 2*power-1
//sauf le derner qui est de degre size+final_degree

//tmpol = [tmpol_1, tmpol_2, tmpol_3, tmpol_4, ... tmpol_{k/2}]

//avec tmpol_i=pol_{2*i} * (X+1)^size donc tous de degre 2*size-1 sauf
//le dernier qui est de degre final_degree+size

static void taylorshift1_dac(mpz_t *upol,
                             const unsigned long int deg,
                             mpz_t *tmpol,
                             mpz_t **shift_pwx,
                             unsigned long int sz,
                             const unsigned int nthreads){

  int l;
  unsigned long int nblocks;
  if(deg <= sz){
    nblocks=0;
  }
  else{
    l = LOG2(deg / sz);
    nblocks = 1<<l;
  }
  if(nblocks<=1){
    taylorshift1_naive(upol, deg);
    return;
  }

  unsigned long int fdeg = deg - (nblocks - 1)*sz;

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  for(unsigned long int i = 0; i <= nblocks - 1; i++){
    unsigned long int cont2 = 0;
    if(i < nblocks-1){
      cont2 = mpz_poly_remove_binary_content(&upol[i * sz],
                                             sz-1);
      taylorshift1_naive(&upol[i * sz], sz-1);
      rescale_upoly_2exp(upol + i * sz, sz - 1, cont2);
    }
    else{
      cont2 = mpz_poly_remove_binary_content(upol+(nblocks-1)*sz,
                                             fdeg);
      taylorshift1_naive(upol+(nblocks-1)*sz, fdeg);
      rescale_upoly_2exp(upol + (nblocks-1)*sz, fdeg, cont2);
    }
  }

  unsigned long int cont2 = 0;
  unsigned long int npowers = LOG2(nblocks);
  mpz_t *shifted = NULL;

  for(unsigned long int i = 0; i < npowers - 1; i++){
    shifted = shift_pwx[i];
    if(nblocks>1){
      fdeg = deg - (nblocks-1)*sz;
    }
    else{
      fdeg = deg - sz;
    }

    nblocks = nblocks / 2;
    for(int count = 1; count <= nblocks ; count++){
      unsigned long int newdeg;
      if(count == nblocks){
        newdeg = fdeg;
      }
      else{
        newdeg = sz - 1;
      }

      cont2 = mpz_poly_remove_binary_content(upol + 2*(count-1)*sz+sz,
                                             newdeg);
      mpz_poly_mul(&tmpol[2*(count-1)*sz],
                   shifted,
                   sz,
                   &upol[2*(count-1)*sz+sz],
                   newdeg,
                   nthreads);

      rescale_upoly_2exp(upol + 2*(count-1)*sz + sz, newdeg, cont2);
      rescale_upoly_2exp(tmpol + 2*(count-1)*sz, newdeg + sz, cont2);

      mpz_poly_add_th(&upol[2*(count-1)*sz], &upol[2*(count-1)*sz], sz-1,
                      &tmpol[2*(count-1)*sz], sz-1, nthreads);

      mpz_poly_swap_th(&upol[2*(count-1)*sz+sz],
                       &tmpol[2*(count-1)*sz + sz],
                       newdeg, nthreads);

    }

    sz = 2 * sz;
    fdeg = deg - (nblocks-1)*sz;

  }

  fdeg = deg - sz;
  shifted = shift_pwx[npowers-1];

  cont2 = mpz_poly_remove_binary_content(upol + sz, fdeg);

  mpz_poly_mul(tmpol, shifted, sz,
               upol+sz, fdeg,
               nthreads);

  rescale_upoly_2exp(upol + sz, fdeg, cont2);
  rescale_upoly_2exp(tmpol, fdeg + sz, cont2);

  mpz_poly_add_th(upol, upol, sz - 1, tmpol, sz - 1,nthreads);
  mpz_poly_swap_th(upol+sz, tmpol+sz, fdeg, nthreads);

}



static inline void decompose_bits_poly(mpz_t **upols, mpz_t *upol,
                                       const unsigned long int deg,
                                       const long int cutoff, const long int nbits,
                                       const long int nblocks){
  for(long int j = nblocks - 1; j >= 0; j--){
    long int trunc = nbits - cutoff * (nblocks - j);
    if(trunc < 0)
      trunc = 0;
    /* #pragma omp parallel for num_threads(nblocks) */
    for(unsigned long int i = 0; i <= deg; i++){
      mpz_tdiv_q_2exp(upols[j][i], upol[i], trunc);
      mpz_tdiv_r_2exp(upol[i], upol[i], trunc);
    }
  }
}




static long mpz_poly_sgn_variations_coeffs_bsize_with_index(mpz_t* upol,
                                                            const unsigned long deg,
                                                            const unsigned long int bsize,
                                                            unsigned long int *index){
  int s = mpz_sgn(upol[deg]);
  if(s==0){
    return -1;
  }
  unsigned long int i;
  long nb = 0;
  int boo = 1;
  unsigned long int N = 1;
  int L = LOG2(bsize);
  for(i = deg - 1; i > 0; i--){
    int c = mpz_cmp_ui(upol[i],0);
    long int l = ilog2_mpz(upol[i]);
    N = min((bsize - i +1)*L + 1, min(deg, (i+1)*L + 1));
    if(l <= N || c==0){
      boo = 0;
    }
    if( mpz_sgn(upol[i]) * s < 0 && l > N && c!=0){
      nb = nb + 1;
      s = mpz_sgn(upol[i]);
      if(nb >= 3){
        *index = i;
        return nb;
      }
    }
  }
  /* int oldnb = nb; */
  N = L + 1;
  int c = mpz_cmp_ui(upol[0],0);
  long int l = ilog2_mpz(upol[0]);
  if(l <= N || c==0){
    boo = 0;
  }
  if(s*mpz_sgn(upol[0]) < 0 && l> N && c!=0){
    nb = nb + 1;
  }
  if(nb>=3 || boo == 1){
    *index = i;
    return nb;
  }
  return -1;
}

/* This is called to perform a Taylor shift after truncation
   of coefficients of pol by maxnbits - 2*(deg + 1).

   Hence the first deg bits of the computed Taylor shifts are correct
   (actually, the sgn variation count functions use a more accurate
   estimate)
*/
static long taylorshift1_dac_wsgnvar(mpz_t *pol,
                                     const unsigned long int deg,
                                     mpz_t *tmpol,
                                     mpz_t **shift_pwx,
                                     unsigned long int pwx,
                                     const unsigned int nthreads){


  /* nblocks is the largest power of 2 dividing deg / pwx */
  int l ;
  unsigned long int nblocks ;
  if(deg <= pwx){
    nblocks=0;
  }
  else{
    l = LOG2(deg / pwx);
    nblocks = 1<<l;
  }
  if(nblocks<=1){
    return -1;
  }

  unsigned long int fdeg = deg - (nblocks - 1)*pwx;

  int i;
  long nb;
  taylorshift1_naive(pol+(nblocks-1)*(pwx), fdeg );

#ifdef _OPENMP
  omp_set_num_threads(nthreads);
#endif
#pragma omp parallel for num_threads(nthreads) //schedule(static)
  for(unsigned long int block = 0; block < nblocks-1; block++){
    taylorshift1_naive(&pol[block*(pwx)], pwx - 1);
  }

  nb = USOLVEmpz_poly_sgn_variations_coeffs_bsize(pol, deg, deg);
  if(nb==0) return 0;

  unsigned long int npwr = LOG2(nblocks);
  mpz_t *current_shifted_powerX = NULL;

  for(i=0;i<npwr-1; i++){
    current_shifted_powerX = shift_pwx[i];

    if(nblocks>1){
      fdeg = deg - (nblocks-1)*pwx;
    }
    else{
      fdeg = deg - pwx;
    }

    nblocks = nblocks / 2;

    for(int count=nblocks; count >= 1 ; count--){

      unsigned long int ndeg;

      if(count==nblocks){
        ndeg = fdeg;
      }
      else{
        ndeg = pwx - 1;
      }

      mpz_poly_mul(&tmpol[2*(count-1)*pwx],
                   current_shifted_powerX,
                   pwx,
                   &pol[2*(count-1)*pwx+pwx],
                   ndeg,
                   nthreads);
      mpz_poly_add_th(&pol[2*(count-1)*pwx], &pol[2*(count-1)*pwx],
                      pwx-1, &tmpol[2*(count-1)*pwx],
                      pwx-1, nthreads);
      mpz_poly_swap_th(&pol[2*(count-1)*pwx+pwx], &tmpol[2*(count-1)*pwx + pwx],
                      ndeg, nthreads);


    }
    pwx = 2 * pwx;
    fdeg = deg - (nblocks-1)*pwx;


  }

  fdeg = deg - pwx;
  current_shifted_powerX = shift_pwx[npwr-1];

  mpz_poly_mul(tmpol,
               current_shifted_powerX,
               pwx,
               pol+pwx,
               fdeg,
               nthreads);

  mpz_poly_add_th(pol, pol, pwx-1, tmpol, pwx-1,nthreads);
  mpz_poly_swap_th(pol+pwx, tmpol+pwx, fdeg, nthreads);

  unsigned long int index = 0;
  nb = mpz_poly_sgn_variations_coeffs_bsize_with_index(pol, deg, deg, &index);
  if(nb>=0) {
    return nb;
  }

  return -1;
}

