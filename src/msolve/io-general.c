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

