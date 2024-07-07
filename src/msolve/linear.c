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

void (*set_linear_poly)(nvars_t nlins, uint32_t *lineqs, nvars_t *linvars,
                        ht_t *bht, int32_t *bexp_lm, bs_t *bs);

void (*check_and_set_linear_poly)(nvars_t *nlins_ptr, nvars_t *linvars,
                                  uint32_t **lineqs_ptr, ht_t *bht,
                                  int32_t *bexp_lm, bs_t *bs);

static inline void set_linear_poly_8(nvars_t nlins, uint32_t *lineqs,
                                     nvars_t *linvars, ht_t *bht,
                                     int32_t *bexp_lm, bs_t *bs) {

  const int nv = bht->nv;
  const len_t ebl = bht->ebl;
  const len_t evl = bht->evl;
  len_t ctr = 0;
  exp_t *etmp = (exp_t *)calloc((unsigned long)nv, sizeof(exp_t));
  for (long i = 0; i < nlins * ((nv + 1)); i++) {
    lineqs[i] = 0;
  }
  int cnt = 0;

  for (int i = 0; i < nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_8[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_8[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          /* convert to usual exponent vector without block elimination storage
           * structure */
          ctr = 0;
          for (int k = 1; k < ebl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }
          for (int k = ebl + 1; k < evl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }

          int isvar = 0;
          for (int k = 0; k < nv; k++) {
            if (etmp[k] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  free(etmp);
}

static inline void set_linear_poly_16(nvars_t nlins, uint32_t *lineqs,
                                      nvars_t *linvars, ht_t *bht,
                                      int32_t *bexp_lm, bs_t *bs) {

  const int nv = bht->nv;
  const len_t ebl = bht->ebl;
  const len_t evl = bht->evl;
  len_t ctr = 0;
  exp_t *etmp = (exp_t *)calloc((unsigned long)nv, sizeof(exp_t));
  for (long i = 0; i < nlins * ((nv + 1)); i++) {
    lineqs[i] = 0;
  }
  /* for(long i = 0; i < nlins*(bht->nv+1); i++){
   *   lineqs[i] = 0;
   * } */
  int cnt = 0;

  for (int i = 0; i < nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_16[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_16[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          /* convert to usual exponent vector without block elimination storage
           * structure */
          ctr = 0;
          for (int k = 1; k < ebl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }
          for (int k = ebl + 1; k < evl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }

          int isvar = 0;
          for (int k = 0; k < nv; k++) {
            if (etmp[k] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  free(etmp);
}

static inline void set_linear_poly_32(nvars_t nlins, uint32_t *lineqs,
                                      nvars_t *linvars, ht_t *bht,
                                      int32_t *bexp_lm, bs_t *bs) {

  const int nv = bht->nv;
  const len_t ebl = bht->ebl;
  const len_t evl = bht->evl;
  len_t ctr = 0;
  exp_t *etmp = (exp_t *)calloc((unsigned long)nv, sizeof(exp_t));
  for (long i = 0; i < nlins * ((nv + 1)); i++) {
    lineqs[i] = 0;
  }
  int cnt = 0;

  for (int i = 0; i < nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          /* convert to usual exponent vector without block elimination storage
           * structure */
          ctr = 0;
          for (int k = 1; k < ebl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }
          for (int k = ebl + 1; k < evl; ++k) {
            etmp[ctr++] = (int32_t)exp[k];
          }

          int isvar = 0;
          for (int k = 0; k < nv; k++) {
            if (etmp[k] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  free(etmp);
}

static inline void check_and_set_linear_poly_8(nvars_t *nlins_ptr,
                                               nvars_t *linvars,
                                               uint32_t **lineqs_ptr, ht_t *bht,
                                               int32_t *bexp_lm, bs_t *bs) {
  long nlins = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1.
  */
  for (long i = 0; i < bs->lml; i++) {
    long deg = 0;
    for (int j = 0; j < bht->nv; j++) {
      deg += bexp_lm[i * bht->nv + j];
    }

    if (deg == 1) {
      nlins++;
      for (int k = 0; k < bht->nv; k++) {
        if (bexp_lm[i * bht->nv + k] == 1) {
          linvars[k] = i + 1;
        }
      }
    }
  }

  *nlins_ptr = nlins;

  // On recupere les coefficients des formes lineaires
  uint32_t *lineqs = calloc(nlins * (bht->nv + 1), sizeof(uint32_t));
  int cnt = 0;
  for (int i = 0; i < bht->nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_8[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (bht->nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_8[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          int isvar = 0;
          for (int k = 0; k < bht->nv; k++) {
            /* exponent vectors in hash table store the degree
             * at the first position, thus "+1" */
            if (exp[k + 1] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  lineqs_ptr[0] = lineqs;
}

static inline void check_and_set_linear_poly_16(nvars_t *nlins_ptr,
                                                nvars_t *linvars,
                                                uint32_t **lineqs_ptr,
                                                ht_t *bht, int32_t *bexp_lm,
                                                bs_t *bs) {
  long nlins = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1.
  */
  for (long i = 0; i < bs->lml; i++) {
    long deg = 0;
    for (int j = 0; j < bht->nv; j++) {
      deg += bexp_lm[i * bht->nv + j];
    }

    if (deg == 1) {
      nlins++;
      for (int k = 0; k < bht->nv; k++) {
        if (bexp_lm[i * bht->nv + k] == 1) {
          linvars[k] = i + 1;
        }
      }
    }
  }

  *nlins_ptr = nlins;

  // On recupere les coefficients des formes lineaires
  uint32_t *lineqs = calloc(nlins * (bht->nv + 1), sizeof(uint32_t));
  int cnt = 0;
  for (int i = 0; i < bht->nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_16[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (bht->nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_16[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          int isvar = 0;
          for (int k = 0; k < bht->nv; k++) {
            /* exponent vectors in hash table store the degree
             * at the first position, thus "+1" */
            if (exp[k + 1] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  lineqs_ptr[0] = lineqs;
}

static inline void check_and_set_linear_poly_32(nvars_t *nlins_ptr,
                                                nvars_t *linvars,
                                                uint32_t **lineqs_ptr,
                                                ht_t *bht, int32_t *bexp_lm,
                                                bs_t *bs) {
  long nlins = 0;
  /*
    la i-ieme entree de linvars est a 0 si il n'y a pas de forme lineaire dont
    le terme dominant est vars[i].
    Sinon, on met l'indice du polynome dans la base + 1.
  */
  for (long i = 0; i < bs->lml; i++) {
    long deg = 0;
    for (int j = 0; j < bht->nv; j++) {
      deg += bexp_lm[i * bht->nv + j];
    }

    if (deg == 1) {
      nlins++;
      for (int k = 0; k < bht->nv; k++) {
        if (bexp_lm[i * bht->nv + k] == 1) {
          linvars[k] = i + 1;
        }
      }
    }
  }

  *nlins_ptr = nlins;

  // On recupere les coefficients des formes lineaires
  uint32_t *lineqs = calloc(nlins * (bht->nv + 1), sizeof(uint32_t));
  int cnt = 0;
  for (int i = 0; i < bht->nv; i++) {
    if (linvars[i] != 0) {

      long len = bs->hm[bs->lmps[linvars[i] - 1]][LENGTH];

      const bl_t bi = bs->lmps[linvars[i] - 1];
      if (len == bht->nv + 1) {
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          lineqs[cnt * (bht->nv + 1) + j] = coef;
        }
      } else {
        hm_t *dt = bs->hm[bi] + OFFSET;
        for (long j = 0; j < len; j++) {
          uint32_t coef = bs->cf_32[bs->hm[bi][COEFFS]][j];
          exp_t *exp = bht->ev[dt[j]];
          int isvar = 0;
          for (int k = 0; k < bht->nv; k++) {
            /* exponent vectors in hash table store the degree
             * at the first position, thus "+1" */
            if (exp[k + 1] == 1) {
              lineqs[cnt * (bht->nv + 1) + k] = coef;
              isvar = 1;
            }
          }
          if (isvar == 0) {
            lineqs[cnt * (bht->nv + 1) + bht->nv] = coef;
          }
        }
        cnt++;
      }
    }
  }
  lineqs_ptr[0] = lineqs;
}

static void inline compute_modular_linear_forms(int nlins, nvars_t sz,
                                                uint32_t *mod_linear_forms,
                                                mpz_t *mpz_linear_forms,
                                                uint32_t prime) {
  for (int i = 0; i < nlins; i++) {
    nvars_t n = i * sz;
    nvars_t n2 = i * (sz + 1);
    uint64_t inv = mpz_fdiv_ui(mpz_linear_forms[n2 + sz], prime);
    inv = mod_p_inverse_32(inv, prime);
    for (nvars_t j = 0; j < sz; j++) {
      mod_linear_forms[n + j] = mpz_fdiv_ui(mpz_linear_forms[n2 + j], prime);
      mod_linear_forms[n + j] =
          (((uint64_t)mod_linear_forms[n + j]) * inv) % prime;
    }
  }
#ifdef DEBUGLIFTMAT
  fprintf(stderr, "\nModular linear forms (prime = %u)\n", prime);
  for (int i = 0; i < nlins; i++) {
    nvars_t n = i * sz;
    nvars_t n2 = i * (sz + 1);
    for (int j = 0; j < sz; j++) {
      fprintf(stderr, "%d, ", mod_linear_forms[n + j]);
    }
    fprintf(stderr, "\n");
  }
#endif
}

static inline mpz_t *allocate_crt_linear_forms(int nlins, int nv,
                                               uint32_t **lineqs_ptr) {
  mpz_t *crt_linear_forms = malloc(sizeof(mpz_t) * nlins * (nv + 1));
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < (nv + 1); j++) {
      mpz_init(crt_linear_forms[i * (nv + 1) + j]);
      mpz_set_ui(crt_linear_forms[i * (nv + 1) + j],
                 lineqs_ptr[0][i * (nv + 1) + j]);
    }
  }
  return crt_linear_forms;
}

static inline mpz_t *allocate_mpq_linear_forms(int nlins, int nv) {
  mpz_t *mpq_linear_forms = malloc(2 * sizeof(mpz_t) * nlins * (nv + 1));
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < 2 * (nv + 1); j++) {
      mpz_init(mpq_linear_forms[2 * i * (nv + 1) + j]);
    }
  }
  return mpq_linear_forms;
}

static inline void crt_linear_forms_clear(mpz_t *crt_linear_forms, int nlins,
                                          int nv) {
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < (nv + 1); j++) {
      mpz_clear(crt_linear_forms[i * (nv + 1) + j]);
    }
  }
  free(crt_linear_forms);
}

static inline void mpq_linear_forms_clear(mpz_t *mpq_linear_forms, int nlins,
                                          int nv) {
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < 2 * (nv + 1); j++) {
      mpz_clear(mpq_linear_forms[2 * i * (nv + 1) + j]);
    }
  }
  free(mpq_linear_forms);
}

static inline mpz_t *mpz_linear_forms_allocate(int nlins, int nv) {
  mpz_t *mpz_linear_forms = malloc(sizeof(mpz_t) * nlins * (nv + 2));
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < (nv + 2); j++) {
      mpz_init(mpz_linear_forms[i * (nv + 2) + j]);
    }
  }
  return mpz_linear_forms;
}

static inline void mpz_linear_forms_clear(mpz_t *mpz_linear_forms, int nlins,
                                          int nv) {
  for (int32_t i = 0; i < nlins; i++) {
    for (int32_t j = 0; j < (nv + 2); j++) {
      mpz_clear(mpz_linear_forms[i * (nv + 2) + j]);
    }
  }
  free(mpz_linear_forms);
}
