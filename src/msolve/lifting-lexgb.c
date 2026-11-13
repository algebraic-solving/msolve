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

/* Reduced lex Gb of a zero-dimensional ideal over Q by multi-modular FGLM:
 *  - learning on a first prime p0: F4, multiplication matrices and FGLM,
 *    recording their trace (fglm_trace_t);
 *  - application on the next primes: traced F4, matrices rebuilt from the
 *    trace (one call to core_nf), FGLM reduced to a linear system;
 *  - lifting: CRT and rational reconstruction with the machinery of
 *    lifting-gb.c (ratrecon_gb), checked on extra primes.
 * Included by msolve.c, after lifting-gb.c. */

#include <flint/fmpz.h>

#define LEXGB_MAX_ATTEMPTS  5   // full restarts (unlucky p0)
#define LEXGB_MAX_BAD      10   // consecutive bad primes before a restart

/* Storage of the lex Gb modulo the primes in the structures of lifting-gb.c,
 * so that its lifting (ratrecon_gb) can be reused: polynomial j has
 * g_ntail[j] coefficients, the one of index k being the coefficient of s_k
 * (index 0 is the smallest monomial, as in modpgbs_set). The fields which
 * describe a DRL basis (hash table, hashed monomials) are not used. */
static void lexgb_modpoly_init(gb_modpoly_t modgbs,
                               const fglm_lex_staircase_t *lex,
                               uint32_t alloc){
  const uint32_t ld = lex->ngens;
  const int nv = lex->nvars;
  modgbs->alloc = alloc;
  modgbs->nprimes = 0;
  modgbs->primes = calloc(alloc, sizeof(mp_limb_t));
  modgbs->cf_64 = calloc(alloc, sizeof(mp_limb_t));
  modgbs->ld = ld;
  modgbs->nv = nv;
  modgbs->modpolys = malloc(sizeof(modpolys_t) * ld);
  modgbs->ldm = malloc((uint64_t)nv * ld * sizeof(int32_t));
  memcpy(modgbs->ldm, lex->g_lm, (uint64_t)nv * ld * sizeof(int32_t));
  modgbs->bht = NULL;
  modgbs->hm = calloc(ld, sizeof(hm_t *));
  modgbs->lmps = calloc(ld, sizeof(bl_t));
  for(uint32_t i = 0; i < ld; i++){
    const uint32_t len = lex->g_ntail[i];
    modgbs->modpolys[i]->len = len;
    modgbs->modpolys[i]->cf_32 = calloc((size_t)len * alloc + 1, sizeof(uint32_t));
    modgbs->modpolys[i]->cf_qq = malloc(sizeof(mpz_t) * (len + 1));
    for(uint32_t j = 0; j < len; j++) mpz_init(modgbs->modpolys[i]->cf_qq[j]);
    mpz_init_set_ui(modgbs->modpolys[i]->lm, 1);
  }
}

/* Coefficients of the lex Gb modulo p (as given by fglm_classic_learn or
 * fglm_classic_apply) */
static void lexgb_modpoly_add(gb_modpoly_t modgbs, const fglm_lex_staircase_t *lex,
                              const CF_t *coeffs, uint32_t p){
  const uint32_t np = modgbs->nprimes;
  if(modgbs->alloc <= np + 2){
    gb_modpoly_realloc(modgbs, MAX(4, np / 16));
  }
  modgbs->primes[np] = p;
  uint64_t off = 0;
  for(int32_t i = 0; i < lex->ngens; i++){
    const uint32_t len = lex->g_ntail[i];
    if(modgbs->modpolys[i]->cf_32 != NULL){
      for(uint32_t k = 0; k < len; k++){
        MODCF(modgbs->modpolys[i], k, np) = coeffs[off + k];
      }
    }
    off += len;
  }
  modgbs->nprimes++;
}

/* Output in the format of print_fglm_gb, with integer coefficients
 * (denominators cleared polynomial by polynomial) */
static void lexgb_print_term(FILE *f, const fmpz_t c, const int32_t *e,
                             const int nv, char **vnames, const int first){
  int one = 1;
  for(int i = 0; i < nv; i++) {
      if(e[i]) {
          one = 0;
      }
  }
  fmpz_t a;
  fmpz_init(a);
  fmpz_abs(a, c);
  if(fmpz_sgn(c) < 0) {
      fprintf(f, "-");
  } else if(!first) {
      fprintf(f, "+");
  }
  if(one){
    fmpz_fprint(f, a);
  } else {
    int star = 0;
    if(!fmpz_is_one(a)){
        fmpz_fprint(f, a);
        star = 1;
    }
    for(int i = 0; i < nv; i++){
      if(!e[i]) continue;
      if(star) fprintf(f, "*");
      if(e[i] == 1) fprintf(f, "%s", vnames[i]);
      else fprintf(f, "%s^%d", vnames[i], e[i]);
      star = 1;
    }
  }
  fmpz_clear(a);
}

static void lexgb_print(FILE *file, const fglm_lex_staircase_t *lex,
                        const gb_modpoly_t modgbs, char **vnames){
  const int nv = lex->nvars;
  fprintf(file, "#Reduced Groebner basis data (change of order by FGLM)\n");
  fprintf(file, "#---\n");
  fprintf(file, "#field characteristic: 0\n");
  fprintf(file, "#variable order:       ");
  for(int i = 0; i < nv - 1; i++) {
      fprintf(file, "%s, ", vnames[i]);
  }
  fprintf(file, "%s\n", vnames[nv - 1]);
  fprintf(file, "#monomial order:       lexicographical\n");
  if(lex->ngens == 1){
    fprintf(file, "#length of basis:      1 element\n");
  } else {
    fprintf(file, "#length of basis:      %d elements ",
            lex->ngens);
    fprintf(file, "sorted by increasing leading monomials\n");
  }
  fprintf(file, "#---\n");
  fmpz_t c;
  fmpz_init(c);
  fprintf(file, "[");
  for(int32_t j = 0; j < lex->ngens; j++){
    const modpolys_struct *P = modgbs->modpolys[j];
    /* lm * LM + sum cf_qq[k] s_k */
    fmpz_set_mpz(c, P->lm);
    lexgb_print_term(file, c, lex->g_lm + (uint64_t)j * nv, nv, vnames, 1);
    for(szmat_t k = lex->g_ntail[j]; k-- > 0; ){
      if(mpz_sgn(P->cf_qq[k]) == 0) continue;
      fmpz_set_mpz(c, P->cf_qq[k]);
      lexgb_print_term(file, c, lex->S + (uint64_t)k * nv, nv, vnames, 0);
    }
    if (j < lex->ngens -1) {
        fprintf(file, ",\n");
    }
  }
  fprintf(file, "]:\n");
  fmpz_clear(c);
}

static FILE *lexgb_open_output(const files_gb *files){
  if(files->out_file != NULL) return fopen(files->out_file, "wb+");
  return OUTSTREAM;
}
static void lexgb_close_output(FILE *f, const files_gb *files){
  if(files->out_file != NULL) fclose(f);
}

static uint32_t lexgb_next_prime(uint32_t prime, const bs_t *bs, uint32_t avoid){
  const uint32_t lprime = 1303905299;
  do{
    prime = next_prime(prime);
    if(prime >= lprime) prime = next_prime(1 << 30);
  } while(is_lucky_prime_ui(prime, bs) || prime == avoid);
  return prime;
}

/* Entry point for -g 3 in characteristic 0.
 * Returns 0 if the basis was printed, -1 otherwise. */
static int print_msolve_lexgb_qq(data_gens_ff_t *gens, msflags_t flags){
  double ca0 = realtime();

  const int32_t info_level = flags->info_level;

  for(int attempt = 0; attempt < LEXGB_MAX_ATTEMPTS; attempt++){
    bs_t *bs  = NULL;
    ht_t *bht = NULL;
    md_t *md  = NULL;
    int success = initialize_gba_input_data(&bs, &bht, &md,
            gens->lens, gens->exps, (void *)gens->mpz_cfs, 0 /* char */,
            0 /* DRL */, flags->elim_block_len, gens->nvars, gens->ngens,
            0, flags->ht_size, flags->nr_threads, flags->max_nr_pairs,
            flags->reset_ht, flags->la_option, flags->use_signatures,
            1 /* reduce_gb */, flags->pbm_file, flags->truncate_lifting,
            info_level);
    if(success == -1){
      fprintf(ERRSTREAM, "All input generators are zero: "
              "the ideal is not zero-dimensional.\n");
      return -1;
    }
    if(success == 0){
      fprintf(ERRSTREAM, "Bad input data, stopped computation.\n");
      exit(1);
    }

    /* learning on p0 */
    uint32_t p0 = lexgb_next_prime(rand() % (1303905301 - (1 << 30) + 1) + (1 << 30),
                                   bs, 0);
    if(info_level){
      fprintf(VERBSTREAM, "Initial prime = %u\n", p0);
    }
    int32_t err = 0;
    md->f4_qq_round = 1;
    bs_t *bsp = core_gba(bs, md, &err, p0);
    if(err || bsp == NULL){
      fprintf(ERRSTREAM, "Problem with F4, stopped computation.\n");
      exit(1);
    }
    const int nv = bsp->ht->nv;
    const long nlm0 = bsp->lml;
    int32_t *lm0 = get_lm_from_bs(bsp, bsp->ht);

    int unit = (nlm0 == 1);
    for(int i = 0; i < nv && unit; i++) {
        if(lm0[i]) {
            unit = 0;
        }
    }
    if(unit || !has_dimension_zero(nlm0, nv, lm0)){
      FILE *file = lexgb_open_output(flags->files);
      if(unit){
        /* the reduced Gb of the unit ideal is {1} for any order; same
         * output as over a prime field */
        fprintf(file, "#Reduced Groebner basis data\n#---\n");
        fprintf(file, "#field characteristic: 0\n");
        fprintf(file, "#variable order:       ");
        for(int i = 0; i < nv - 1; i++) fprintf(file, "%s, ", gens->vnames[i]);
        fprintf(file, "%s\n", gens->vnames[nv - 1]);
        fprintf(file, "#monomial order:       lexicographical\n");
        fprintf(file, "#length of basis:      1 element\n");
        fprintf(file, "#---\n");
        fprintf(file, "[1]:\n");
      }
      lexgb_close_output(file, flags->files);
      if(!unit){
        fprintf(ERRSTREAM, "The ideal is not zero-dimensional: "
                "no change of order performed.\n");
      }
      free(lm0);
      free_basis_without_hash_table(&bsp);
      free_basis_without_hash_table(&bs);
      free_shared_hash_data(bht);
      free_hash_table(&bht);
      free_trace(&md->tr);
      free(md);
      return unit ? 0 : -1;
    }

    md->fglm_rtime = realtime();
    md->fglm_ctime = cputime();
    long dquot = 0;
    int32_t *lmb = monomial_basis(nlm0, nv, lm0, &dquot);
    fglm_trace_t tr;
    memset(&tr, 0, sizeof(tr));
    fglmc_rows_t rp;
    memset(&rp, 0, sizeof(rp));
    print_fglm_header(VERBSTREAM,md);
    sp_multimatfglm_t *M =
        build_multiplication_matrices_from_bs(lmb, dquot, bsp, bsp->ht, lm0,
                                              md, nv, p0, info_level, &tr, &rp);
    /* p0 uses all the threads, each next prime only one */
    M->nthreads = md->nthrds;
    const int old_flint = fglmc_flint_threads_push(md->nthrds);
    CF_t *res = NULL;
    int ok = fglm_classic_learn(M, FGLM_ORDER_LEX, &tr.lex, &tr, &rp, &res, md);
    /* the next primes only compute the rows computed on p0 */
    mm_trace_keep_computed_rows(&tr, M);
    mm_rows_clear(&rp);
    if(res == NULL) {
        res = malloc(sizeof(CF_t));
    }
    fglmc_flint_threads_pop(old_flint);
    md->fglm_rtime = realtime() - md->fglm_rtime;
    md->fglm_ctime = cputime() - md->fglm_ctime;
    md->learning_rtime = realtime()-ca0;
    print_fglm_classic_data(VERBSTREAM, md, M, &tr.lex);
    if(info_level){
      fprintf(VERBSTREAM,"\n---------- COMPUTATIONAL DATA -----------\n");
      if(tr.lex.npols != NULL){
        for(int i = nv-1; i >= 0; i--){
          fprintf(VERBSTREAM, "[%lu]", (unsigned long)tr.lex.npols[i]);
        }
        fprintf(VERBSTREAM, "\n");
      }
      fprintf(VERBSTREAM,
              "#polynomials to lift %14lu\n",
              (unsigned long) tr.lex.ngens);
      fprintf(VERBSTREAM, "-----------------------------------------\n");
    }
    free_sp_multimat_fglm(M);
    free_basis_without_hash_table(&bsp);

    /* lifting with the machinery of lifting-gb.c */
    gb_modpoly_t modgbs;
    data_lift_t dl;
    rrec_data_t rd1, rd2;
    mpz_t mod_p, prod_p;
    mpz_init_set_ui(mod_p, 1);
    mpz_init(prod_p);
    initialize_rrec_data(rd1);
    initialize_rrec_data(rd2);
    lexgb_modpoly_init(modgbs, &tr.lex, 2);
    int32_t step = tr.lex.ngens;
    data_lift_init(dl, tr.lex.ngens, &step, 1);
    double st_crt = 0, st_rrec = 0, st_wit = 0;
    int done = 0;
    if(ok){
      lexgb_modpoly_add(modgbs, &tr.lex, res, p0);
      choose_coef_to_lift(modgbs, dl);
    }

    /* application: as in secondary_modular_steps, nthr primes are handled in
     * parallel (one thread each), then their results are used in order */
    md->info_level = 0;
    const int nthr = md->nthrds > 0 ? md->nthrds : 1;
    md->nthrds = 1;
    md->f4_qq_round = 2;
    uint32_t *bprimes = malloc(nthr * sizeof(uint32_t));
    int *bbad = malloc(nthr * sizeof(int));
    CF_t *bres = malloc((uint64_t)nthr * (tr.lex.ncoeffs + 1) * sizeof(CF_t));
    uint32_t prime = rand() % (1303905301 - (1 << 30) + 1) + (1 << 30);
    int consec_bad = 0;
    uint32_t nbad = 0;
    /* timings of the second prime (first prime of the first batch) */
    double st_f4 = 0, st_fglm = 0, st_prime2 = 0;
    int first_batch = 1;
    /* progress display: {n} = good primes used, p0 included, as for the
     * parametrizations (msolve_trace_qq) */
    int32_t lstart = dl->lstart;
    while(ok && !done && consec_bad < LEXGB_MAX_BAD){
      for(int i = 0; i < nthr; i++){
        prime = lexgb_next_prime(prime, bs, p0);
        bprimes[i] = prime;
      }
      if(info_level && first_batch){
        fprintf(VERBSTREAM, "New prime = %u\n", bprimes[0]);
      }
#pragma omp parallel for num_threads(nthr) schedule(dynamic)
      for(int i = 0; i < nthr; i++){
        CF_t *ri = bres + (uint64_t)i * (tr.lex.ncoeffs + 1);
        const int timed = (first_batch && i == 0);
        double rt_p = timed ? realtime() : 0;
        int32_t lerr = 0;
        bs_t *bsq = core_gba(bs, md, &lerr, bprimes[i]);
        double rt_f = 0;
        if(timed){
          rt_f = realtime();
          st_f4 = rt_f - rt_p;
        }
        int bad = (lerr != 0 || bsq == NULL);
        if(!bad && bsq->lml != nlm0){
            bad = 1;
        }
        if(!bad){
          int32_t *lm = get_lm_from_bs(bsq, bsq->ht);
          if(memcmp(lm, lm0, (uint64_t)nlm0 * nv * sizeof(int32_t))) bad = 1;
          free(lm);
        }
        if(!bad){
          sp_multimatfglm_t *Mq =
              build_multiplication_matrices_from_trace(&tr, lmb, bsq, bsq->ht,
                                                       md, bprimes[i]);
          if(fglm_classic_apply(ri, Mq, &tr)) bad = 1;
          free_sp_multimat_fglm(Mq);
          if(timed) st_fglm = realtime() - rt_f;
        }
        if(bsq != NULL){
            free_basis_and_only_local_hash_table_data(&bsq);
        }
        bbad[i] = bad;
        if(timed){
          st_prime2 = realtime() - rt_p;
          if(info_level){
            fprintf(VERBSTREAM, "\n---------------- TIMINGS ----------------\n");
            fprintf(VERBSTREAM, "multi-mod overall(elapsed) %9.2f sec\n", st_prime2);
            fprintf(VERBSTREAM, "multi-mod F4               %9.2f sec\n", st_f4);
            fprintf(VERBSTREAM, "multi-mod FGLM             %9.2f sec\n", st_fglm);
            fprintf(VERBSTREAM, "-----------------------------------------\n");
            fprintf(VERBSTREAM, "\nmulti-modular steps\n");
            fprintf(VERBSTREAM, "-------------------------------------------------\
-----------------------------------------------------\n");
            fflush(VERBSTREAM);
          }
        }
      }
      first_batch = 0;
      for(int i = 0; i < nthr && !done; i++){
        const CF_t *ri = bres + (uint64_t)i * (tr.lex.ncoeffs + 1);
        if(bbad[i]){
          nbad++;
          consec_bad++;
          if(info_level > 1){
            fprintf(VERBSTREAM, "[!]");
            fflush(VERBSTREAM);
          }
          continue;
        }
        consec_bad = 0;
        lexgb_modpoly_add(modgbs, &tr.lex, ri, bprimes[i]);
        ratrecon_gb(modgbs, dl, mod_p, prod_p, rd1, rd2, nthr,
                    &st_crt, &st_rrec, &st_wit);
        /* same period of reconstruction attempts as in lifting-gb.c */
        int32_t rr = 1;
        while(2 * rr <= (int32_t)modgbs->nprimes / RR_PERIOD_DIV){
            rr *= 2;
        }
        if(rr != dl->rr){
          dl->rr = rr;
          if(info_level){
            fprintf(VERBSTREAM, "(->%d)", dl->rr);
            fflush(VERBSTREAM);
          }
        }
        if(info_level && !(modgbs->nprimes & (modgbs->nprimes - 1))){
          fprintf(VERBSTREAM, "{%u}", modgbs->nprimes);
          fflush(VERBSTREAM);
        }
        if(dl->lstart != lstart){
          if(info_level){
            fprintf(VERBSTREAM, "<%.2f%%>",
                    100 * (float)MIN((dl->lstart + 1), modgbs->ld) / modgbs->ld);
            fflush(VERBSTREAM);
          }
          lstart = dl->lstart;
        }
        done = 1;
        for(int32_t j = 0; j < tr.lex.ngens; j++){
          if(dl->check2[j] < NBCHECK){ done = 0; break; }
        }
      }
    }
    if(info_level && ok){
      fprintf(VERBSTREAM, " \n-------------------------------------------------\
-----------------------------------------------------\n");
    }
    md->nthrds = nthr;
    free(bprimes); free(bbad); free(bres);
    md->info_level = info_level;

    if(done){
      if(info_level){
        fprintf(VERBSTREAM,"\n\n---------- COMPUTATIONAL DATA -----------\n");
        fprintf(VERBSTREAM, "Max coeff. bitsize %16lu\n",
                (unsigned long)max_bit_size_gb(modgbs));
        fprintf(VERBSTREAM, "#primes            %16lu\n", (unsigned long) modgbs->nprimes);
        fprintf(VERBSTREAM, "#bad primes        %16lu\n", (unsigned long) nbad);
        fprintf(VERBSTREAM, "-----------------------------------------\n");
        fprintf(VERBSTREAM, "\n---------------- TIMINGS ----------------\n");
        fprintf(VERBSTREAM, "CRT     (elapsed)         %10.2f sec\n", st_crt);
        fprintf(VERBSTREAM, "ratrecon(elapsed)         %10.2f sec\n", st_rrec);
        fprintf(VERBSTREAM, "-----------------------------------------\n");
      }
      FILE *f = lexgb_open_output(flags->files);
      lexgb_print(f, &tr.lex, modgbs, gens->vnames);
      lexgb_close_output(f, flags->files);
    } else if(info_level){
      fprintf(VERBSTREAM, "Too many bad primes (learning prime %u is probably "
              "unlucky), computation restarts.\n", p0);
    }
    data_lift_clear(dl);
    gb_modpoly_without_hash_table_clear(modgbs);
    free_rrec_data(rd1);
    free_rrec_data(rd2);
    mpz_clear(mod_p);
    mpz_clear(prod_p);
    free(res);
    free_fglm_trace(&tr);
    free(lmb);
    free(lm0);
    free_basis_without_hash_table(&bs);
    free_shared_hash_data(bht);
    free_hash_table(&bht);
    free_trace(&md->tr);
    free(md);
    if(done){
        return 0;
    }
  }
  fprintf(ERRSTREAM, "Lex Groebner basis over Q: too many failures, giving up.\n");
  return -1;
}
