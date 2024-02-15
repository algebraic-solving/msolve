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


#include "meta_data.h"
static md_t *copy_meta_data(
        const md_t * const gmd,
        const int32_t prime
        )
{
    md_t *md = (md_t *)malloc(sizeof(md_t));
    memcpy(md, gmd, sizeof(md_t));
    md->fc  = prime;
    md->min_deg_in_first_deg_fall = gmd->min_deg_in_first_deg_fall;
    md->application_nr_mult = 0;
    md->application_nr_add  = 0;
    md->application_nr_red  = 0;

    if (md->fc < pow(2,7)) {
        md->ff_bits = 8;
    } else {
        if (md->fc < pow(2,15)) {
            md->ff_bits = 16;
        } else {
            if (md->fc < pow(2,31)) {
                md->ff_bits = 32;
            }
        }
    }
    set_ff_bits(md, md->fc);
    return md;
}

void free_meta_data(
     md_t **mdp
     )
{
    md_t *md = *mdp;

    if (md->ps != NULL) {
        free_pairset(&(md->ps));
    }
    free(md->hcm);

    ht_t *ht = md->ht;
    if (ht != NULL) {
        free_hash_table(&ht);
    }

    free(md);
    md = NULL;

    *mdp = md;
}

md_t *allocate_meta_data(
    void
    )
{
    md_t *md  = (md_t *)calloc(1, sizeof(md_t));

    return md;
}

void print_tracer_statistics(
    FILE *file,
    const double rt,
    const md_t *st
    )
{
    if (st->trace_level == APPLY_TRACER) {
        if(st->info_level > 1){
            fprintf(stderr, "Learning phase %.2f Gops/sec\n",
                    (st->trace_nr_add+st->trace_nr_mult)/1000.0/1000.0/(realtime()-rt));
        }
        if(st->info_level > 2){
            fprintf(stderr, "------------------------------------------\n");
            fprintf(stderr, "#ADDITIONS       %13lu\n", (unsigned long)st->trace_nr_add * 1000);
            fprintf(stderr, "#MULTIPLICATIONS %13lu\n", (unsigned long)st->trace_nr_mult * 1000);
            fprintf(stderr, "#REDUCTIONS      %13lu\n", (unsigned long)st->trace_nr_red);
            fprintf(stderr, "------------------------------------------\n");
        }
    }
}

void print_initial_statistics(
    FILE *file,
    const md_t *st
    )
{
    if (st->info_level > 0) {
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
}

void print_round_information_header(
        FILE *f,
        const md_t * const st
        )
{
    if (st->info_level > 1) {
        if (st->trace_level != APPLY_TRACER) {
            fprintf(f, "\n");
            fprintf(f, "Legend for f4 information\n");
            fprintf(f, "--------------------------------------------------------\n");
            fprintf(f, "deg       current degree of pairs selected in this round\n");
            fprintf(f, "sel       number of pairs selected in this round\n");
            fprintf(f, "pairs     total number of pairs in pair list\n");
            fprintf(f, "mat       matrix dimensions (# rows x # columns)\n");
            fprintf(f, "density   density of the matrix\n");
            fprintf(f, "new data  # new elements for basis in this round\n");
            fprintf(f, "          # zero reductions during linear algebra\n");
            fprintf(f, "time(rd)  time of the current f4 round in seconds given\n");
            fprintf(f, "          for real and cpu time\n");
            fprintf(f, "--------------------------------------------------------\n");
            fprintf(f, "\ndeg     sel   pairs        mat          density \
           new data         time(rd) in sec (real|cpu)\n");
            fprintf(f, "-------------------------------------------------\
-----------------------------------------------------\n");
        } else {
            fprintf(f, "Legend for f4 information\n");
            fprintf(f, "--------------------------------------------------------\n");
            fprintf(f, "round     # of current tracer round\n");
            fprintf(f, "deg       current degree of pairs selected in this round\n");
            fprintf(f, "mat       matrix dimensions (# rows x # columns)\n");
            fprintf(f, "density   density of the matrix\n");
            fprintf(f, "new data  # new elements for basis in this round\n");
            fprintf(f, "          # zero reductions during linear algebra\n");
            fprintf(f, "time(rd)  time of the current f4 round in seconds given\n");
            fprintf(f, "          for real and cpu time\n");
            fprintf(f, "--------------------------------------------------------\n");
            fprintf(f, "\n    round     deg          mat          density \
           new data         time(rd) in sec (real|cpu)\n");
            fprintf(f, "-------------------------------------------------\
-----------------------------------------------------\n");
        }
    }
}

void print_sat_nf_round_timings(
        FILE *f,
        const md_t * const st,
        const double rrt,
        const double crt
        )
{
    if (st->info_level > 1) {
        printf("%15.2f | %-13.2f\n",
                realtime() - rrt, cputime() - crt);
    }
}

void print_sat_round_timings(
        FILE *f,
        const md_t * const st,
        const double rrt,
        const double crt
        )
{
    if (st->info_level > 1) {
        printf("%10.2f | %-13.2f\n",
                realtime() - rrt, cputime() - crt);
    }
}

void print_round_timings(
        FILE *f,
        const md_t * const st,
        const double rrt,
        const double crt
        )
{
    if (st->info_level > 1) {
        printf("%13.2f | %-13.2f\n",
                realtime() - rrt, cputime() - crt);
    }
}

void print_round_information_footer(
        FILE *f,
        const md_t * const st
        )
{
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
-----------------------------------------------------\n");
    }
}

static void print_current_trace_meta_data(
        md_t *md
        )
{
    len_t rd  = md->trace_rd;
    len_t deg = md->tr->td[rd].deg;

    if (md->info_level > 1) {
        printf("%9d  %6d  ", rd+1, deg);
        fflush(stdout);
    }
}

static void get_final_statistics(
        md_t *md,
        const bs_t * const bs
        )
{
    len_t i = 0;
    int64_t nterms = 0;
    md->size_basis = bs->lml;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->hm[bs->lmps[i]] == NULL) {
            nterms++;
        } else {
            nterms += bs->hm[bs->lmps[i]][LENGTH];
        }
    }
    md->nterms_basis = nterms;

}

void get_and_print_final_statistics(
        FILE *file, 
        md_t *st,
        const bs_t * const bs
        )
{
    get_final_statistics(st, bs);

    if (st->info_level > 0) {
        fprintf(file, "\n---------------- TIMINGS ---------------\n");
        fprintf(file, "overall(elapsed) %11.2f sec\n", st->f4_rtime);
        fprintf(file, "overall(cpu) %15.2f sec\n", st->f4_ctime);
        if(st->trace_level == APPLY_TRACER) {
            fprintf(file, "tracer       %15.2f sec %5.1f%%\n",
                    st->tracer_rtime,
                    (double)100*(double)st->tracer_rtime
                    / (double)(st->f4_rtime));
        } else {
            fprintf(file, "select       %15.2f sec %5.1f%%\n",
                    st->select_rtime,
                    (double)100*(double)st->select_rtime
                    / (double)(st->f4_rtime));
            fprintf(file, "symbolic prep.       %7.2f sec %5.1f%%\n",
                    st->symbol_rtime,
                    (double)100*(double)st->symbol_rtime
                    / (double)(st->f4_rtime));
            fprintf(file, "update       %15.2f sec %5.1f%%\n",
                    st->update_rtime,
                    (double)100*(double)st->update_rtime
                    / (double)(st->f4_rtime));
        }
        fprintf(file, "convert      %15.2f sec %5.1f%%\n",
                st->convert_rtime,
                (double)100*(double)st->convert_rtime
                / (double)(st->f4_rtime));
        fprintf(file, "linear algebra   %11.2f sec %5.1f%%\n",
                st->la_rtime,
                (double)100*(double)st->la_rtime
                / (double)(st->f4_rtime));
        if (st->reduce_gb == 1) {
            fprintf(file, "reduce gb    %15.2f sec %5.1f%%\n",
                    st->reduce_gb_rtime,
                    (double)100*(double)st->reduce_gb_rtime
                    / (double)(st->f4_rtime));
        }
        if (st->reset_ht != 2147483647) {
            fprintf(file, "rht          %15.2f sec %5.1f%%\n",
                    st->rht_rtime,
                    (double)100*(double)st->rht_rtime
                    / (double)(st->f4_rtime));
        }
        fprintf(file, "-----------------------------------------\n");
        fprintf(file, "\n---------- COMPUTATIONAL DATA -----------\n");
        fprintf(file, "size of basis      %16lu\n", (unsigned long)st->size_basis);
        fprintf(file, "#terms in basis    %16lu\n", (unsigned long)st->nterms_basis);
        fprintf(file, "#pairs reduced     %16lu\n", (unsigned long)st->num_pairsred);
        fprintf(file, "#GM criterion      %16lu\n", (unsigned long)st->num_gb_crit);
        fprintf(file, "#redundant elements      %10lu\n", (unsigned long)st->num_redundant);
        fprintf(file, "#rows reduced      %16lu\n", (unsigned long)st->num_rowsred);
        fprintf(file, "#zero reductions   %16lu\n", (unsigned long)st->num_zerored);
        fprintf(file, "max. matrix data   %16ld x %ld (%.3f%%)\n",
                (long)st->mat_max_nrows, (long)st->mat_max_ncols, st->mat_max_density);
        fprintf(file, "max. symbolic hash table size  2^%d\n",
                (int32_t)(ceil(log((double)st->max_sht_size)/log(2))));
        fprintf(file, "max. basis hash table size     2^%d\n",
                (int32_t)(ceil(log((double)st->max_bht_size)/log(2))));
        fprintf(file, "-----------------------------------------\n\n");
    }
}
