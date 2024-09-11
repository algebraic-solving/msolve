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
 * Vincent Neiger
 * Mohab Safey El Din */

#include <stdint.h>
#include <stdio.h>
#include "msolve-data.h"

/*---------------------------------*/
/* Integer vector helper functions */
/*---------------------------------*/

void print_vec(FILE * file, uint32_t * vec, uint32_t len)
{
    fprintf(file, "[");
    for(uint32_t i = 0; i < len-1; ++i)
        fprintf(file, "%u, ", vec[i]);
    fprintf(file, "%u]\n",vec[len-1]);
}

/*---------------------------*/
/* FGLM I/O helper functions */
/*---------------------------*/

void display_fglm_param(FILE * file, param_t * param)
{
    fprintf(file, "%ld,\n", param->charac);
    fprintf(file, "%d,\n", param->nvars);

    nmod_poly_fprint(file, param->elim);
    fprintf(file, ",\n");
    nmod_poly_fprint(file, param->denom);
    fprintf(file, ",\n");
    fprintf(file, "[");
    for(int c = param->nvars-2; c >= 0; c--)
    {
        nmod_poly_fprint(file, param->coords[c]);
        fprintf(file, "\n");
    }
    fprintf(file, "]");
}

void display_fglm_param_maple(FILE * file, param_t * param)
{
    fprintf(file, "[%ld, \n", param->charac);
    fprintf(file, "%d, \n", param->nvars);

    nmod_poly_fprint(file, param->elim);
    fprintf(file, ", \n");
    nmod_poly_fprint(file, param->denom);
    fprintf(file, ", \n");

    for(int c = param->nvars-2; c > 0; c--)
    {
        nmod_poly_fprint(file, param->coords[c]);
        fprintf(file, ", \n");
    }
    nmod_poly_fprint(file, param->coords[0]);
    fprintf(file, "]:\n");
}

void print_fglm_data(
        FILE *file,
        const md_t * const st,
	sp_matfglm_t *matrix,
	param_t *param
		     )
{
  if (st->info_level > 0) {
    fprintf(file, "\n---------------- TIMINGS ----------------\n");
    fprintf(file, "overall(elapsed) %11.2f sec\n", st->fglm_rtime);
    fprintf(file, "overall(cpu) %15.2f sec\n", st->fglm_ctime);
    fprintf(file, "-----------------------------------------\n");
    fprintf(file, "\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(file, "degree of ideal    %16lu\n", (unsigned long)matrix->ncols);
    fprintf(file, "#dense rows        %16lu\n", (unsigned long)matrix->nrows);
    fprintf(file, "#normal forms      %16lu\n", (unsigned long)matrix->nnfs);
    fprintf(file, "total density                   %5.1f%%\n", 100*matrix->totaldensity);
    fprintf(file, "density of the free part        %5.1f%%\n", 100*matrix->freepartdensity);
    if(matrix->nnfs){
      fprintf(file, "density of the nonfree part     %5.1f%%\n", 100*matrix->nonfreepartdensity);
    }
    fprintf(file, "deg. elim. pol.    %16lu\n", (unsigned long)param->degelimpol);
    fprintf(file, "deg. sqfr. elim. pol. %13lu\n", (unsigned long)param->degsqfrelimpol);
    fprintf(file, "-----------------------------------------\n\n");
  }
}
