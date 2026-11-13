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

/* -L: lifting of the multiplication matrices. All the T_i share dense_mat
 * and rowsel, which does not depend on the prime: the rows of dense_mat
 * computed on the first prime (tails of the DRL Gb, normal forms) are lifted
 * over Q, each one seen as the monic polynomial m_r - NF(m_r) with D
 * coefficients, with the machinery used for the lex Gb. Once lifted, the
 * matrices modulo a new prime are obtained by reduction: neither F4 nor
 * core_nf is needed any more. */
#define LEXGB_MATB 32           // primes per batch once the matrices are lifted
#define LEXGB_MATBUF (1 << 24)  // max. coefficients of these batches

typedef struct {
  int lifted;
  szmat_t D, nnfs;    // columns, rows of the full dense_mat
  szmat_t nrows;      // lifted rows
  szmat_t *rows;      // their indices in dense_mat
  gb_modpoly_t mod;   // residues and lifting (ratrecon_gb)
  data_lift_t dl;
  rrec_data_t rd1, rd2;
  mpz_t mod_p, prod_p;
  double st_crt, st_rrec, st_wit;
  uint32_t nprimes;   // primes used to lift them
  fmpz *num;          // once lifted: nrows * D numerators
  fmpz *den;          // and one denominator per row
} lexgb_matlift_t;

static void lexgb_matlift_init(lexgb_matlift_t *ML, const fglm_trace_t *tr){
  memset(ML, 0, sizeof(*ML));
  ML->D = tr->D;
  ML->nnfs = tr->nnfs;
  ML->rows = malloc(((uint64_t)tr->nnfs + 1) * sizeof(szmat_t));
  for(szmat_t r = 0; r < tr->nnfs; r++){
    if(tr->have_row == NULL || tr->have_row[r]) ML->rows[ML->nrows++] = r;
  }
  gb_modpoly_array_struct *m = ML->mod;
  const uint32_t ld = ML->nrows, alloc = 2;
  m->alloc = alloc;
  m->nprimes = 0;
  m->primes = calloc(alloc, sizeof(mp_limb_t));
  m->cf_64 = calloc(alloc, sizeof(mp_limb_t));
  m->ld = ld;
  m->nv = tr->nvars;
  m->modpolys = malloc(sizeof(modpolys_t) * (ld + 1));
  m->ldm = calloc((uint64_t)tr->nvars * ld + 1, sizeof(int32_t));
  m->bht = NULL;
  m->hm = calloc(ld + 1, sizeof(hm_t *));
  m->lmps = calloc(ld + 1, sizeof(bl_t));
  for(uint32_t i = 0; i < ld; i++){
    m->modpolys[i]->len = ML->D;
    m->modpolys[i]->cf_32 = calloc((size_t)ML->D * alloc + 1, sizeof(uint32_t));
    m->modpolys[i]->cf_qq = malloc(sizeof(mpz_t) * (ML->D + 1));
    for(szmat_t j = 0; j < ML->D; j++) mpz_init(m->modpolys[i]->cf_qq[j]);
    mpz_init_set_ui(m->modpolys[i]->lm, 1);
  }
  int32_t step = ld;
  data_lift_init(ML->dl, ld, &step, 1);
  initialize_rrec_data(ML->rd1);
  initialize_rrec_data(ML->rd2);
  mpz_init_set_ui(ML->mod_p, 1);
  mpz_init(ML->prod_p);
}

static void lexgb_matlift_clear(lexgb_matlift_t *ML){
  if(ML->num != NULL){
    _fmpz_vec_clear(ML->num, (slong)ML->nrows * ML->D + 1);
    _fmpz_vec_clear(ML->den, (slong)ML->nrows + 1);
  }
  data_lift_clear(ML->dl);
  gb_modpoly_without_hash_table_clear(ML->mod);
  free_rrec_data(ML->rd1);
  free_rrec_data(ML->rd2);
  mpz_clear(ML->mod_p);
  mpz_clear(ML->prod_p);
  free(ML->rows);
}

/* Lifted rows: row r of dense_mat is num[r] / den[r] (den[r] = lm) */
static void lexgb_matlift_finalize(lexgb_matlift_t *ML){
  const szmat_t D = ML->D;
  ML->num = _fmpz_vec_init((slong)ML->nrows * D + 1);
  ML->den = _fmpz_vec_init((slong)ML->nrows + 1);
  for(szmat_t i = 0; i < ML->nrows; i++){
    const modpolys_struct *P = ML->mod->modpolys[i];
    fmpz_set_mpz(ML->den + i, P->lm);
    for(szmat_t k = 0; k < D; k++){
      fmpz_set_mpz(ML->num + (uint64_t)i * D + k, P->cf_qq[k]);
    }
  }
  ML->nprimes = ML->mod->nprimes;
  ML->lifted = 1;
}

/* Rows of dm (full dense_mat modulo p) for the lifting; returns 1 once the
 * rows are lifted and checked (as for the lex Gb: NBCHECK primes) */
static int lexgb_matlift_add(lexgb_matlift_t *ML, const CF_t *dm, uint32_t p){
  if(ML->lifted) return 1;
  gb_modpoly_array_struct *m = ML->mod;
  const uint32_t np = m->nprimes;
  if(m->alloc <= np + 2){
    gb_modpoly_realloc(m, MAX(4, np / 16));
  }
  m->primes[np] = p;
  for(szmat_t i = 0; i < ML->nrows; i++){
    /* residues freed by ratrecon_gb once the row is lifted and checked */
    if(m->modpolys[i]->cf_32 == NULL) continue;
    const CF_t *row = dm + (uint64_t)ML->rows[i] * ML->D;
    for(szmat_t k = 0; k < ML->D; k++){
      MODCF(m->modpolys[i], k, np) = row[k];
    }
  }
  m->nprimes++;
  if(np == 0){
    choose_coef_to_lift(m, ML->dl);
    return 0;
  }
  ratrecon_gb(m, ML->dl, ML->mod_p, ML->prod_p, ML->rd1, ML->rd2, 1,
              &ML->st_crt, &ML->st_rrec, &ML->st_wit);
  int32_t rr = 1;
  while(2 * rr <= (int32_t)m->nprimes / RR_PERIOD_DIV){
    rr *= 2;
  }
  ML->dl->rr = rr;
  for(szmat_t j = 0; j < ML->nrows; j++){
    if(ML->dl->check2[j] < NBCHECK) return 0;
  }
  lexgb_matlift_finalize(ML);
  return 1;
}

/* Lifted rows modulo primes[0..K): row ML->rows[r] of the dense_mat of
 * prime k is written in out + k * stride (stride = nnfs * D + 1, the other
 * rows are not written); bad[k] = 1 if primes[k] divides a denominator. */
static void lexgb_matlift_residues(CF_t *out, const uint64_t stride, int *bad,
                                   const lexgb_matlift_t *ML,
                                   const mp_limb_t *primes, const slong K,
                                   const int nthr){
  const szmat_t D = ML->D, nr = ML->nrows;
  fmpz_comb_t comb;
  fmpz_comb_init(comb, primes, K);
  mp_limb_t *inv = malloc(((uint64_t)nr * K + 1) * sizeof(mp_limb_t));
  for(slong k = 0; k < K; k++) bad[k] = 0;
  fmpz_comb_temp_t temp0;
  fmpz_comb_temp_init(temp0, comb);
  fmpz_t a0;
  fmpz_init(a0);
  for(szmat_t r = 0; r < nr; r++){
    mp_limb_t *iv = inv + (uint64_t)r * K;
    /* residues of |x|, then sign (as for the numerators below) */
    fmpz_abs(a0, ML->den + r);
    fmpz_multi_mod_ui(iv, a0, comb, temp0);
    const int neg = fmpz_sgn(ML->den + r) < 0;
    for(slong k = 0; k < K; k++){
      if(iv[k] == 0){
        bad[k] = 1;
      } else {
        iv[k] = n_invmod(neg ? primes[k] - iv[k] : iv[k], primes[k]);
      }
    }
  }
  fmpz_clear(a0);
  fmpz_comb_temp_clear(temp0);
#pragma omp parallel num_threads(nthr)
  {
    fmpz_comb_temp_t temp;
    fmpz_comb_temp_init(temp, comb);
    mp_limb_t *res = malloc((K + 1) * sizeof(mp_limb_t));
    fmpz_t a;
    fmpz_init(a);
#pragma omp for schedule(dynamic)
    for(slong r = 0; r < (slong)nr; r++){
      const mp_limb_t *iv = inv + (uint64_t)r * K;
      const uint64_t pos = (uint64_t)ML->rows[r] * D;
      for(szmat_t l = 0; l < D; l++){
        const fmpz *x = ML->num + (uint64_t)r * D + l;
        if(fmpz_is_zero(x)){
          for(slong k = 0; k < K; k++) out[k * stride + pos + l] = 0;
          continue;
        }
        fmpz_abs(a, x);
        fmpz_multi_mod_ui(res, a, comb, temp);
        const int neg = fmpz_sgn(x) < 0;
        for(slong k = 0; k < K; k++){
          const mp_limb_t v = (neg && res[k]) ? primes[k] - res[k] : res[k];
          out[k * stride + pos + l] =
              bad[k] ? 0 : (CF_t)(((uint64_t)v * iv[k]) % primes[k]);
        }
      }
    }
    fmpz_clear(a);
    free(res);
    fmpz_comb_temp_clear(temp);
  }
  free(inv);
  fmpz_comb_clear(comb);
}

/* Matrices modulo p with the dense array dm (not copied: the caller sets
 * dense_mat to NULL before free_sp_multimat_fglm) */
static sp_multimatfglm_t *lexgb_matrices_view(const fglm_trace_t *tr, CF_t *dm,
                                              const uint32_t p){
  const szmat_t D = tr->D;
  sp_multimatfglm_t *M = calloc(1, sizeof(sp_multimatfglm_t));
  M->charac = p;
  M->nvars = tr->nvars;
  M->ncols = D;
  M->nnfs = tr->nnfs;
  M->nlms = tr->nnfs - tr->nnf;
  M->one = 0;
  M->dense_mat = dm;
  M->rowsel = malloc((uint64_t)tr->nvars * D * sizeof(szmat_t));
  memcpy(M->rowsel, tr->rowsel, (uint64_t)tr->nvars * D * sizeof(szmat_t));
  M->is_lm = malloc(tr->nnfs + 1);
  memcpy(M->is_lm, tr->is_lm, tr->nnfs + 1);
  if(tr->have_row != NULL){
    M->have_row = malloc(tr->nnfs + 1);
    memcpy(M->have_row, tr->have_row, tr->nnfs + 1);
  }
  return M;
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
    /* -L: the rows of dense_mat computed on p0 are lifted too */
    const int use_ml = flags->lift_matrix && ok;
    lexgb_matlift_t ML;
    if(use_ml){
      lexgb_matlift_init(&ML, &tr);
      lexgb_matlift_add(&ML, M->dense_mat, p0);
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
    /* -L: rows of dense_mat of each prime of a batch (before the lifting of
     * the matrices) or lifted matrices modulo the primes of a batch (after:
     * batches of up to LEXGB_MATB primes) */
    const uint64_t mstride = (uint64_t)tr.nnfs * tr.D + 1;
    int Kmax = nthr;
    if(use_ml){
      uint64_t kb = LEXGB_MATBUF / mstride;
      if(kb > LEXGB_MATB) kb = LEXGB_MATB;
      if(kb > (uint64_t)Kmax) Kmax = (int)kb;
    }
    uint32_t *bprimes = malloc(Kmax * sizeof(uint32_t));
    mp_limb_t *lprimes = malloc(Kmax * sizeof(mp_limb_t));
    int *bbad = malloc(Kmax * sizeof(int));
    int *pbad = malloc(Kmax * sizeof(int));
    CF_t *bres = malloc((uint64_t)Kmax * (tr.lex.ncoeffs + 1) * sizeof(CF_t));
    CF_t *bmat = use_ml ? malloc((uint64_t)Kmax * mstride * sizeof(CF_t)) : NULL;
    uint32_t nskip = 0, nf4 = 0;
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
      const int lifted = use_ml && ML.lifted;
      const int K = lifted ? Kmax : nthr;
      for(int i = 0; i < K; i++){
        prime = lexgb_next_prime(prime, bs, p0);
        bprimes[i] = prime;
        lprimes[i] = prime;
      }
      if(lifted){
        lexgb_matlift_residues(bmat, mstride, pbad, &ML, lprimes, K, nthr);
      } else {
        nf4 += K;
      }
      if(info_level && first_batch){
        fprintf(VERBSTREAM, "New prime = %u\n", bprimes[0]);
      }
#pragma omp parallel for num_threads(nthr) schedule(dynamic)
      for(int i = 0; i < K; i++){
        CF_t *ri = bres + (uint64_t)i * (tr.lex.ncoeffs + 1);
        if(lifted){
          /* lifted matrices: neither F4 nor core_nf */
          int bad = 0;
          if(pbad[i]){
            bad = 3;                   /* p divides a denominator */
          } else {
            sp_multimatfglm_t *Mq =
                lexgb_matrices_view(&tr, bmat + (uint64_t)i * mstride,
                                    bprimes[i]);
            if(fglm_classic_apply(ri, Mq, &tr)) bad = 1;
            Mq->dense_mat = NULL;
            free_sp_multimat_fglm(Mq);
          }
          bbad[i] = bad;
          continue;
        }
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
          if(use_ml && !bad){
            /* rows for the lifting of the matrices */
            CF_t *bm = bmat + (uint64_t)i * mstride;
            for(szmat_t q = 0; q < ML.nrows; q++){
              const uint64_t o = (uint64_t)ML.rows[q] * tr.D;
              memcpy(bm + o, Mq->dense_mat + o, tr.D * sizeof(CF_t));
            }
          }
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
      for(int i = 0; i < K && !done; i++){
        const CF_t *ri = bres + (uint64_t)i * (tr.lex.ncoeffs + 1);
        if(bbad[i] == 3){
          nskip++;
          continue;
        }
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
        if(use_ml && !lifted && !ML.lifted){
          if(lexgb_matlift_add(&ML, bmat + (uint64_t)i * mstride, bprimes[i])
             && info_level > 1){
            fprintf(VERBSTREAM, "[L]");
            fflush(VERBSTREAM);
          }
        }
        lexgb_modpoly_add(modgbs, &tr.lex, ri, bprimes[i]);
        /* one new prime: the checks use the last prime only (thrds = 1) */
        ratrecon_gb(modgbs, dl, mod_p, prod_p, rd1, rd2, 1,
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
    free(lprimes); free(pbad); free(bmat);
    md->info_level = info_level;

    if(done){
      if(info_level){
        fprintf(VERBSTREAM,"\n\n---------- COMPUTATIONAL DATA -----------\n");
        fprintf(VERBSTREAM, "Max coeff. bitsize %16lu\n",
                (unsigned long)max_bit_size_gb(modgbs));
        fprintf(VERBSTREAM, "#primes            %16lu\n", (unsigned long) modgbs->nprimes);
        fprintf(VERBSTREAM, "#bad primes        %16lu\n", (unsigned long) nbad);
        if(use_ml){
          fprintf(VERBSTREAM, "#primes with F4    %16lu\n", (unsigned long) nf4 + 1);
          if(ML.lifted){
            fprintf(VERBSTREAM, "#primes matrices   %16lu\n",
                    (unsigned long) ML.nprimes);
            fprintf(VERBSTREAM, "matrices bitsize   %16lu\n",
                    (unsigned long)FLINT_ABS(_fmpz_vec_max_bits(ML.num,
                                                  (slong)ML.nrows * ML.D)));
          }
          fprintf(VERBSTREAM, "#skipped primes    %16lu\n", (unsigned long) nskip);
        }
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
    if(use_ml){
      lexgb_matlift_clear(&ML);
    }
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
