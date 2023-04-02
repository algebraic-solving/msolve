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

param_t *nmod_fglm_compute(sp_matfglm_t *, mod_t, long,
                           long, uint64_t *, uint32_t *, uint64_t*, int);
param_t *nmod_fglm_guess_colon(sp_matfglmcol_t *, mod_t,
			       CF_t *, CF_t **, long,
			       long, uint64_t *, uint32_t *, uint64_t *, int);
param_t *nmod_fglm_compute_trace_data(sp_matfglm_t *, mod_t, long,
                                      long,
                                      long,
                                      uint64_t *,
                                      uint32_t *,
                                      uint64_t*,
                                      int,
                                      fglm_data_t **,
                                      fglm_bms_data_t **,
                                      int *,
				      stat_t *);
int nmod_fglm_compute_apply_trace_data(sp_matfglm_t *,
                                       const mod_t,
                                       param_t *,
                                       const long,
                                       const long,
                                       const long,
                                       uint64_t *,
                                       uint32_t *,
                                       uint64_t*,
                                       fglm_data_t *,
                                       fglm_bms_data_t *,
                                       const long,
                                       const int,
				       const stat_t *);

void display_fglm_param(FILE *, param_t *);
void display_fglm_param_maple(FILE *, param_t *);
void display_nmod_poly(FILE *, nmod_poly_t);
