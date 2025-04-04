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

static void initialize_mstrace(mstrace_t msd, md_t *st){
  msd->lp  = (primes_t *)calloc(st->nthrds, sizeof(primes_t));

  /*******************
   * initialize basis
   *******************/
  msd->bs_qq = initialize_basis(st);

  msd->tht = NULL;

  /* generate array to store modular bases */
  msd->bs = (bs_t **)calloc((unsigned long)st->nthrds, sizeof(bs_t *));
  msd->bad_primes = calloc((unsigned long)st->nthrds, sizeof(int));

  /* initialize tracers */
  msd->btrace = (trace_t **)calloc(st->nthrds,
                                        sizeof(trace_t *));
  msd->btrace[0]  = initialize_trace(msd->bs_qq, st);
  /* initialization of other tracers is done through duplication */

  msd->num_gb = (int32_t *)calloc(st->nthrds, sizeof(int32_t));
  msd->leadmons_ori = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));
  msd->leadmons_current = (int32_t**)calloc(st->nthrds, sizeof(int32_t *));

  /* array to store one monomial */
  msd->mgb = calloc(st->nvars, sizeof(uint32_t));

  mpz_init(msd->mod_p);
  mpz_set_ui(msd->mod_p, 1);

  mpz_init(msd->prod_p);
  mpz_set_ui(msd->prod_p, 1);

}

static void free_mstrace(mstrace_t msd, md_t *st){
  free_lucky_primes(&msd->lp);
  free(msd->lp);
  /* to be checked if that is to be done when st->ff_bits != 0
     This was previously done only when characteristic is zero
   */

  fprintf(stderr, "removed free_basis in free_mstrace\n");
  //  free_basis(&(msd->bs_qq));
//  free(msd->bs_qq);

  fprintf(stderr, "Check that when st->ff_bits != 0 hash table should be cleared in modgbs\n");
  /***********************************************************
    to be checked if that is to be done when st->ff_bits != 0
     This was previously done only when characteristic is zero
  ************************************************************/

  if(msd->tht!=NULL){
    free_hash_table(&(msd->tht));
  }
  free(msd->tht);

  for(int i = 0; i < st->nthrds; ++i){
    if (msd->bs[i] != NULL) {
      free_basis(&(msd->bs[i]));
    }
  }
  free(msd->bs);

  free(msd->bad_primes);

  for(int i = 0; i < st->nthrds; ++i){
    if(msd->btrace[i] != NULL){
      free_trace(&(msd->btrace[i]));
    }
  }
  free(msd->btrace);

  free(msd->num_gb);

  for(int i = 0; i < st->nthrds; ++i){
    if(msd->leadmons_ori[i] != NULL){
      free(msd->leadmons_ori[i]);
    }
  }
  free(msd->leadmons_ori);

  for(int i = 0; i < st->nthrds; ++i){
    if(msd->leadmons_current[i] != NULL){
      free(msd->leadmons_current[i]);
    }
  }
  free(msd->leadmons_current);

  free(msd->mgb);
  for(int i = 1; i < st->nthrds; i++){
    if(msd->btht[i] != NULL){
      free_hash_table((msd->btht)+i);
    }
  }
  free(msd->btht);

  mpz_clear(msd->mod_p);
  mpz_clear(msd->prod_p);

}
