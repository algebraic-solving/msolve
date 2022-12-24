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

static void initialize_mstrace(mstrace_t msd, stat_t *st){
  msd->lp  = (primes_t *)calloc(st->nthrds, sizeof(primes_t));

  msd->bs_qq = initialize_basis(st);
  msd->bht = initialize_basis_hash_table(st);
  msd->tht = initialize_secondary_hash_table(msd->bht, st);

  msd->bs = (bs_t **)calloc((unsigned long)st->nthrds, sizeof(bs_t *));
  msd->bad_primes = calloc((unsigned long)st->nthrds, sizeof(int));

  msd->btrace = (trace_t **)calloc(st->nthrds,
                                        sizeof(trace_t *));
  msd->btrace[0]  = initialize_trace();

  msd->num_gb = (int32_t *)calloc(st->nthrds, sizeof(int32_t));
  msd->leadmons_ori = (int32_t **)calloc(st->nthrds, sizeof(int32_t *));
  msd->leadmons_current = (int32_t**)calloc(st->nthrds, sizeof(int32_t *));

  msd->mgb = calloc(sizeof(uint32_t), msd->bht->nv);

  msd->blht = (ht_t **)malloc((st->nthrds) * sizeof(ht_t *));
  btht = (ht_t **)malloc((st->nthrds) * sizeof(ht_t *));

}

static void free_mstrace(mstrace_t msd, stat_t *st){
  free(msd->lp);
  /* to be checked if that is to be done when st->ff_bits != 0 */
  free(msd->bs_qq);

  free_shared_hash_data(msd->bht);
  if(msd->bht!=NULL){
    free_hash_table(&(msd->bht));
  }
  free(msd->bht);

  if(msd->tht!=NULL){
    free_hash_table(&(msd->tht));
  }
  free(msd->tht);

  free(msd->bs);
  free(msd->bad_primes);
  free(msd->btrace);

  free(msd->num_gb);
  free(msd->leadmons_ori);
  free(msd->leadmons_current);

  free(msd->mgb);

  for(int i = 0; i < st->nthrds; i++){
    free_hash_table(msd->(blht+i));
    free_hash_table(msd->(btht+i));
  }

}
