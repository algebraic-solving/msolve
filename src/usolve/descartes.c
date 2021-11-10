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

#include<stdio.h>
#include<stdlib.h>

#ifndef USOLVE
#define USOLVE
#endif

static unsigned long descartes_classical(const mpz_t *, mpz_t *,
                                         const unsigned long, long, long *);

static long descartes_truncate(mpz_t *, const unsigned long,
                               const unsigned long,
                               long, long *, usolve_flags *);

/* one assumes upol[0] != 0 */
/* stopped as soon as more than 3 sgn variations are found */
static long mpz_poly_sgn_variations_coeffs(mpz_t* upol, unsigned long deg){
  unsigned long int i;
  long nb = 0;
  int s = mpz_sgn(upol[deg]);

  for(i = deg-1; i > 0; i--){
    if( mpz_sgn(upol[i]) * s < 0 ){
      nb = nb + 1;
      s = mpz_sgn(upol[i]);
      if(nb >= 3){
        return nb;
      }
    }
  }
  if(s*mpz_sgn(upol[0]) < 0){
    if(nb==0){
    }
    nb = nb + 1;
  }
  return nb;
}

static unsigned long int descartes(mpz_t *upol1, mpz_t *upol2,
                                   const unsigned long deg,
                                   long sigh, long *flag, usolve_flags *flags){
  unsigned long int i;
#pragma omp parallel for num_threads(flags->nthreads)
  for(i = 0; i <= deg; i++){
    mpz_set(upol2[i],upol1[deg-i]);
  }

  /* Max bit size of coefficients in upol1 */
  const unsigned long nbits = mpz_poly_max_bsize_coeffs(upol1, deg);
  long nb = descartes_truncate(upol2, deg, nbits, sigh, flag, flags);

  if(nb >= 0){
    return nb;
  }

#pragma omp parallel for num_threads(flags->nthreads)
  for(i = 0; i <= deg; i++){
    mpz_set(upol2[i], upol1[deg-i]);
  }

  taylorshift1_dac(upol2, deg, flags->tmpol,
                   flags->shift_pwx, flags->pwx, flags->nthreads);
  nb = mpz_poly_sgn_variations_coeffs(upol2, deg);

  return nb;
}

/* assumes upol1 and upol2 are the same polynomials */
static long descartes_truncate(mpz_t *upol2,
                               const unsigned long deg,
                               const unsigned long nbits, long sigh,
                               long *flag, usolve_flags *flags){
  int i;
  /* /\* Max bit size of coefficients in upol1 *\/ */
  /* const unsigned long int nbits = mpz_poly_max_bsize_coeffs(upol1, deg); */

  const unsigned long int lc = ilog2_mpz(upol2[deg]);
  /* One wants to divide all coefficients by 2^trunc */
  /* Hence only 2*(deg+1) bits are taken into account */
  /* After applying a Taylor shift we'll lose deg+1 bits */
  /* Should be sufficient for sign variations of coefficients
   in most cases*/
  unsigned long int trunc = nbits - 2 * (deg + 1);

  /*Necessary to avoid that lead coef becomes 0 after truncation*/
  if(lc<=trunc){
    trunc = lc - 1;
  }
  /* All coefficients are divided by 2^trunc */
  /* Only 2*(deg + 1) coeffs are meaningful */
#pragma omp parallel for num_threads(flags->nthreads)
  for(i = 0; i <= deg; i++){
    mpz_tdiv_q_2exp(upol2[i],upol2[i],trunc);
  }

  long check = taylorshift1_dac_wsgnvar(upol2, deg,
                                        flags->tmpol,
                                        flags->shift_pwx,
                                        flags->pwx,
                                        flags->nthreads);
  if(check >= 0){
    return check;
  }

  return -1;
}

/* upol -> numer( upol( 1 / (x+1) ) ) */
static unsigned long descartes_classical(const mpz_t *upol, mpz_t *tpol,
                                         const unsigned long deg,
                                         long hsgn, long *bsgn){
  unsigned long nb = 0;

  long j = deg;
  long t = mpz_sgn(upol[j]);
  while (j >= 0 && mpz_sgn(upol[j]) == t){
    j--;
  }
  /* no sgn change in upol */
  if (j < 0){
      *bsgn = -1;
      return nb;
  }

  for (long i = 0; i <= deg; i++) {
    mpz_set(tpol[i], upol[i]);
  }

  for (long i = 0; i <= deg - 1; i++){
    mpz_add(tpol[i+1], tpol[i+1], tpol[i]);
  }

  int s = mpz_sgn(tpol[deg]);

  *bsgn = s && (s == mpz_sgn(upol[0])) && (s == -hsgn);

  for (long i = 1; i <= deg - 1; i++){
    j = deg - i;
    t = s;
    while (t == 0){
      t = mpz_sgn(tpol[j]); j--;
    }
    while (j >= 0 && mpz_sgn(tpol[j]) == t){
      j--;
    }
    if (j < 0){
        return nb;
    }
    for (j = 0; j <= deg - i - 1; j++){
      mpz_add(tpol[j+1], tpol[j+1], tpol[j]);
    }
    if (s == 0){
      s = mpz_sgn(tpol[deg-i]);
    }
    else{
      if (s == -mpz_sgn(tpol[deg-i])){
        if((nb == 1 && !*bsgn) || nb == 2){
          return (nb + 1);
        }
        nb++;
        s = -s;
      }
    }
  }
  if(s == -mpz_sgn(tpol[0])){
    nb++;
  }
  return nb;
}


