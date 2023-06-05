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


#include "stat.h"
static stat_t *copy_statistics(
        const stat_t * const gst,
        const int32_t prime
        )
{
    stat_t *st  = (stat_t *)malloc(sizeof(stat_t));
    memcpy(st, gst, sizeof(stat_t));
    st->fc  = prime;
    st->application_nr_mult = 0;
    st->application_nr_add  = 0;
    st->application_nr_red  = 0;

    if (st->fc < pow(2,7)) {
        st->ff_bits = 8;
    } else {
        if (st->fc < pow(2,15)) {
            st->ff_bits = 16;
        } else {
            if (st->fc < pow(2,31)) {
                st->ff_bits = 32;
            }
        }
    }
    set_ff_bits(st, st->fc);
    return st;
}

stat_t *initialize_statistics(
    void
    )
{
    stat_t *st  = (stat_t *)calloc(1, sizeof(stat_t));

    return st;
}

void print_initial_statistics(
                              FILE *file,
                              const stat_t *st
    )
{

    fprintf(file, "\n--------------- INPUT DATA ---------------\n");
    fprintf(file, "#variables             %11d\n", st->nvars);
    fprintf(file, "#equations             %11d\n", st->ngens);
    fprintf(file, "#invalid equations     %11d\n", st->ngens_invalid);
    fprintf(file, "field characteristic   %11u\n", st->fc);
    fprintf(file, "homogeneous input?     %11d\n", st->homogeneous);
    fprintf(file, "signature-based computation %6d\n", st->use_signatures);
    if (st->mo == 0 && st->nev == 0) {
        fprintf(file, "monomial order                 DRL\n");
    }
    if (st->mo == 0 && st->nev > 0) {
        fprintf(file, "monomial order             ELIM(%d)\n", st->nev);
    }
    if (st->mo == 1 && st->nev == 0) {
        fprintf(file, "monomial order                 LEX\n");
    }
    if ((st->mo != 0) && (st->mo != 1)) {
        fprintf(file, "monomial order           DONT KNOW\n");
    }
    if (st->reset_ht == 2147483647) {
        fprintf(file, "basis hash table resetting     OFF\n");
    } else {
        fprintf(file, "basis hash table resetting  %6d\n", st->reset_ht);
    }
    fprintf(file, "linear algebra option  %11d\n", st->laopt);
    fprintf(file, "initial hash table size %10lu (2^%d)\n",
            (unsigned long)pow(2,st->init_hts), st->init_hts);
    if (st->mnsel == 2147483647) {
        fprintf(file, "max pair selection             ALL\n");
    } else {
        fprintf(file, "max pair selection     %11d\n", st->mnsel);
    }
    fprintf(file, "reduce gb              %11d\n", st->reduce_gb);
    fprintf(file, "#threads               %11d\n", st->nthrds);
    fprintf(file, "info level             %11d\n", st->info_level);
    fprintf(file, "generate pbm files     %11d\n", st->gen_pbm_file);
    fprintf(file, "------------------------------------------\n");
}

void print_final_statistics(
                            FILE *file, 
                            const stat_t * const st
        )
{
    fprintf(file, "\n---------------- TIMINGS ---------------\n");
    fprintf(file, "overall(elapsed) %11.2f sec\n", st->overall_rtime);
    fprintf(file, "overall(cpu) %15.2f sec\n", st->overall_ctime);
    fprintf(file, "select       %15.2f sec %5.1f%%\n",
            st->select_rtime,
            (double)100*(double)st->select_rtime
            / (double)(st->overall_rtime));
    fprintf(file, "symbolic prep.       %7.2f sec %5.1f%%\n",
            st->symbol_rtime,
            (double)100*(double)st->symbol_rtime
            / (double)(st->overall_rtime));
    fprintf(file, "update       %15.2f sec %5.1f%%\n",
            st->update_rtime,
            (double)100*(double)st->update_rtime
            / (double)(st->overall_rtime));
    fprintf(file, "convert      %15.2f sec %5.1f%%\n",
            st->convert_rtime,
            (double)100*(double)st->convert_rtime
            / (double)(st->overall_rtime));
    fprintf(file, "linear algebra   %11.2f sec %5.1f%%\n",
            st->la_rtime,
            (double)100*(double)st->la_rtime
            / (double)(st->overall_rtime));
    if (st->reduce_gb == 1) {
        fprintf(file, "reduce gb    %15.2f sec %5.1f%%\n",
                st->reduce_gb_rtime,
                (double)100*(double)st->reduce_gb_rtime
                / (double)(st->overall_rtime));
    }
    if (st->reset_ht != 2147483647) {
        fprintf(file, "rht          %15.2f sec %5.1f%%\n",
                st->rht_rtime,
                (double)100*(double)st->rht_rtime
                / (double)(st->overall_rtime));
    }
    fprintf(file, "-----------------------------------------\n");
    fprintf(file, "\n---------- COMPUTATIONAL DATA -----------\n");
    fprintf(file, "size of basis      %16lu\n", (unsigned long)st->size_basis);
    fprintf(file, "#terms in basis    %16lu\n", (unsigned long)st->nterms_basis);
    fprintf(file, "#pairs reduced     %16lu\n", (unsigned long)st->num_pairsred);
    fprintf(file, "#GM criterion      %16lu\n", (unsigned long)st->num_gb_crit);
    fprintf(file, "#redundant elements      %10lu\n", (unsigned long)st->num_redundant);
    fprintf(file, "#reset basis hash table    %8lu\n", (unsigned long)st->num_rht);
    fprintf(file, "#rows reduced      %16lu\n", (unsigned long)st->num_rowsred);
    fprintf(file, "#zero reductions   %16lu\n", (unsigned long)st->num_zerored);
    fprintf(file, "max. update hash table size    2^%d\n",
            (uint32_t)(ceil(log((double)st->max_uht_size)/log(2))));
    fprintf(file, "max. symbolic hash table size  2^%d\n",
            (int32_t)(ceil(log((double)st->max_sht_size)/log(2))));
    fprintf(file, "max. basis hash table size     2^%d\n",
            (int32_t)(ceil(log((double)st->max_bht_size)/log(2))));
    fprintf(file, "-----------------------------------------\n\n");
}
