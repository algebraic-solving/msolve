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

typedef struct
{
  mpz_t numer;
  long k;
  unsigned int isexact;
  int sign_left;
} interval;

typedef struct{
  int search;/*when >0 (resp. <0, =0) computes positive (resp. negative, all) roots */
  long int bound_pos; /*log2(M) where M dominates the largest positive root*/
  long int bound_neg; /*log2(-M) where M is less than the least negative root*/
  int sign;/* when >=0 (resp. <0), one computes non-negative (resp. negative) roots */
  int revert; /* when set to 1 one reads input coeffs by ascending degree */
  int prec_isole; /* nbr of bits for real root isolation */
  int hasrealroots; /*when set to 1, computation stops when one finds a root*/

  long int precision_loss;

  unsigned long int transl;
  unsigned long int node_looked;
  unsigned long int half_done;

  unsigned long int cur_deg;
  unsigned long int pwx;
  unsigned long int nblocks;
  unsigned long int npwr;
  mpz_t **shift_pwx;
  mpz_t *tmpol;
  mpz_t *tmpol_desc;
  /* used for parallel taylor_shift */
  mpz_t **tmp_threads;
  mpz_t **pols_threads;
  mpz_t *Values;

  float time_desc;
  float time_shift;

  unsigned int nthreads;
  unsigned int verbose;
  unsigned int bfile;
  unsigned int classical_algo;

  unsigned int print_stats;
  int debug;
} usolve_flags;
