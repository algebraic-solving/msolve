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

/** compare two data_gens_ff_t which have been filled with int32_t coefficients
 * returns 0 if equal, returns nonzero otherwise
 **/
int data_gens_cmp_int32(data_gens_ff_t * gens1, data_gens_ff_t * gens2)
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

    long pos = 0;
    int32_t c1, c2;
    const long nv_m_el = gens1->nvars - gens1->elim;
    for(long i = 0; i < gens1->ngens; i++)
    {
        if (gens1->lens[i] != gens2->lens[i])
            return 1;

        for(long j = 0; j < gens1->lens[i]; j++)
        {
            c1 = gens1->cfs[pos+j];
            c2 = gens2->cfs[pos+j];
            if (c1 != c2)
                return 1;

            for(long k = 0; k < nv_m_el; k++)
                if (gens1->exps[pos * nv_m_el + k] != gens2->exps[pos * nv_m_el + k])
                    return 1;
        }

        pos += gens1->lens[i]; // TODO move to loop
    }

    return 0;
}


#define TEST_IOFILES_NUMBER_FILES 4

/** this tests that input files are interpreted the same way
 * with dos line endings or unix line endings
 **/

int main(void)
{

    char * fn_dos[TEST_IOFILES_NUMBER_FILES] = {
        "files/in1_dos.ms",
        "files/in2_dos.ms",
        "files/in3_dos.ms",
        "files/in4_dos.ms",
    };

    char * fn_unix[TEST_IOFILES_NUMBER_FILES] = {
        "files/in1_unix.ms",
        "files/in2_unix.ms",
        "files/in3_unix.ms",
        "files/in4_unix.ms",
    };

    for (long k = 0; k < TEST_IOFILES_NUMBER_FILES; k++)
    {
        printf("testing file %ld\n", k);
        char * fn1 = fn_dos[k];
        char * fn2 = fn_unix[k];
        int32_t nr_vars1, nr_vars2, field_char1, field_char2, nr_gens1, nr_gens2;
        data_gens_ff_t * gens1 = allocate_data_gens();
        data_gens_ff_t * gens2 = allocate_data_gens();

        get_data_from_file(fn1, &nr_vars1, &field_char1, &nr_gens1, gens1);
        get_data_from_file(fn2, &nr_vars2, &field_char2, &nr_gens2, gens2);

        printf("get data ok\n");

        if ((nr_vars1 != nr_vars2) || (field_char1 != field_char2) || (nr_gens1 != nr_gens2))
            return 1;

        if (data_gens_cmp_int32(gens1, gens2))
            return 1;
    }

    return 0;
}
