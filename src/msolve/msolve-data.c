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
}

static void free_mstrace(mstrace_t msd, stat_t *st){
  free(msd->lp);
  /* to be checked if that is to be done when st->ff_bits != 0 */
  free(msd->bs_qq);
  free(msd->bht);
  free(msd->tht);
  free(msd->bs);
}
