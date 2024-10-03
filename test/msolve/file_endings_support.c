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

#include "../../src/msolve/libmsolve.c"

typedef struct{
  nvars_t nvars;
  nvars_t elim;
  int32_t ngens;
  int32_t nterms;
  int32_t field_char;
  /* counts change of variable orders:
   * x1 <-> xn
   * x1 <-> xn-1
   * ...
   * for situation when staircase is not generic enough */
  int32_t change_var_order;
  /* base coefficient for linear form
   * sum(i^k*x[k]) k = 1, ..., nvars
   * It is zero if no linear form is active, otherwise != zero. */
  int32_t linear_form_base_coef;
  /* set to 1 if a linear form is chosen randomly */
  int32_t rand_linear;
  int32_t *random_linear_form;
  char **vnames;
  int32_t *lens;
  int32_t *exps;
  int32_t *cfs;    /* int32_t coeffs */
  mpz_t **mpz_cfs; /* mpz_t coeffs */
} data_gens_ff_t;

int data_gens_cmp(data_gens_ff_t * gens1, data_gens_ff_t * gens2)
{
    if (gens1->nvars != gens2->nvars) return 1;
    if (gens1->elim != gens2->elim) return 1;
    if (gens1->ngens != gens2->ngens) return 1;
    if (gens1->nterms != gens2->nterms) return 1;
    if (gens1->field_char != gens2->field_char) return 1;
    if (gens1->change_var_order != gens2->change_var_order) return 1;
    if (gens1->linear_form_base_coef != gens2->linear_form_base_coef) return 1;
    if (gens1->rand_linear != gens2->rand_linear) return 1;
    if (gens1->rand_linear)
        for (int k = 0; k < gens1->nvars; k++)
            if (gens1->random_linear_form[k] != gens2->random_linear_form[k])
                return 1;
    for (int k = 0; k < gens1->nvars; k++)
        if (strcmp(gens1->vnames[k], gens2->vnames[k]))
            return 1;

    // TODO
}


#define TEST_IOFILES_NUMBER_FILES 4

/** this tests that input files are interpreted the same way
 * with dos line endings or unix line endings
 **/

int main(void)
{

    char * fn_dos[TEST_IOFILES_NUMBER_FILES] = {
        "in1_dos.ms",
        "in2_dos.ms",
        "in3_dos.ms",
        "in4_dos.ms",
    };

    char * fn_unix[TEST_IOFILES_NUMBER_FILES] = {
        "in1_unix.ms",
        "in2_unix.ms",
        "in3_unix.ms",
        "in4_unix.ms",
    };

    for (int k = 0; k < TEST_IOFILES_NUMBER_FILES; k++)
    {
        char * fn1 = fn_dos[k];
        char * fn2 = fn_unix[k];
        int32_t nr_vars1, nr_vars2, field_char1, field_char2, nr_gens1, nr_gens2;
        data_gens_ff_t * gens1 = allocate_data_gens();
        data_gens_ff_t * gens2 = allocate_data_gens();

        get_data_from_file(fn1, &nr_vars1, &field_char1, &nr_gens1, gens1);
        get_data_from_file(fn2, &nr_vars2, &field_char2, &nr_gens2, gens2);

        if ( (nr_vars1 != nr_vars2) || (field_char1 != field_char2) || (nr_gens1 != nr_gens2) )
            return 1;

        if (data_gens_cmp(gens1, gens2))
            return 1;
    }

    return 0;
}
