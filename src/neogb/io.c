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


#include "io.h"

/* See exponent vector description in data.h for more information. */
static inline void set_exponent_vector(
        exp_t *ev,
        const int32_t *iev,  /* input exponent vectors */
        const int32_t idx,
        const ht_t *ht,
        const md_t *st
        )
{
    len_t i;

    const len_t nv  = ht->nv;
    const len_t ebl = ht->ebl;
    const len_t nev = st->nev;
    const len_t off = ebl - nev + 1;

    ev[0]   = 0;
    ev[ebl] = 0;

    for (i = 0; i < nev; ++i) {
        ev[i+1] = (exp_t)(iev+(nv*idx))[i];
        /* degree */
        ev[0]   +=  ev[i+1];
    }
    for (i = nev; i < nv; ++i) {
        ev[i+off] = (exp_t)(iev+(nv*idx))[i];
        /* degree */
        ev[ebl]    +=  ev[i+off];
    }
}

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
void sort_terms_ff_8(
    cf8_t **cfp,
    hm_t **hmp,
    ht_t *ht
    )
{
  cf8_t *cf = *cfp;
  hm_t *hm  = *hmp;
  hm_t *hmo = hm+OFFSET;

  const len_t len = hm[LENGTH];

  len_t i, j, k;

  hm_t tmphm  = 0;
  cf8_t tmpcf = 0;

  /* generate array of pointers to hm entries */
  hm_t *phm[len];
  for (i = 0; i < len; ++i) {
    phm[i]  = &hmo[i];
  }

  /* sort pointers to hm entries -> getting permutations */
  sort_r(phm, (unsigned long)len, sizeof(phm[0]), initial_gens_cmp, ht);

  /* sort cf and hm using permutations stored in phm */
  for (i = 0; i < len; ++i) {
    if (i != phm[i]-hmo) {
      tmpcf = cf[i];
      tmphm = hmo[i];
      k     = i;
      while (i != (j = phm[k]-hmo)) {
        cf[k]   = cf[j];
        hmo[k] = hmo[j];
        phm[k] = &hmo[k];
        k      = j;
      }
      cf[k]   = tmpcf;
      hmo[k]  = tmphm;
      phm[k]  = &hmo[k];
    }
  }

  *cfp  = cf;
  *hmp  = hm;
}

void sort_terms_ff_16(
    cf16_t **cfp,
    hm_t **hmp,
    ht_t *ht
    )
{
  cf16_t *cf  = *cfp;
  hm_t *hm    = *hmp;
  hm_t *hmo   = hm+OFFSET;

  const len_t len = hm[LENGTH];

  len_t i, j, k;

  hm_t tmphm    = 0;
  cf16_t tmpcf  = 0;

  /* generate array of pointers to hm entries */
  hm_t *phm[len];
  for (i = 0; i < len; ++i) {
    phm[i]  = &hmo[i];
  }

  /* sort pointers to hm entries -> getting permutations */
  sort_r(phm, (unsigned long)len, sizeof(phm[0]), initial_gens_cmp, ht);

  /* sort cf and hm using permutations stored in phm */
  for (i = 0; i < len; ++i) {
    if (i != phm[i]-hmo) {
      tmpcf = cf[i];
      tmphm = hmo[i];
      k     = i;
      while (i != (j = phm[k]-hmo)) {
        cf[k]   = cf[j];
        hmo[k]  = hmo[j];
        phm[k]  = &hmo[k];
        k       = j;
      }
      cf[k]   = tmpcf;
      hmo[k]  = tmphm;
      phm[k]  = &hmo[k];
    }
  }

  *cfp  = cf;
  *hmp  = hm;
}

void sort_terms_ff_32(
    cf32_t **cfp,
    hm_t **hmp,
    ht_t *ht
    )
{
  cf32_t *cf  = *cfp;
  hm_t *hm    = *hmp;
  hm_t *hmo   = hm+OFFSET;

  const len_t len = hm[LENGTH];

  len_t i, j, k;

  hm_t tmphm    = 0;
  cf32_t tmpcf  = 0;

  /* generate array of pointers to hm entries */
  hm_t *phm[len];
  for (i = 0; i < len; ++i) {
    phm[i]  = &hmo[i];
  }

  /* sort pointers to hm entries -> getting permutations */
  sort_r(phm, (unsigned long)len, sizeof(phm[0]), initial_gens_cmp, ht);

  /* sort cf and hm using permutations stored in phm */
  for (i = 0; i < len; ++i) {
    if (i != phm[i]-hmo) {
      tmpcf = cf[i];
      tmphm = hmo[i];
      k     = i;
      while (i != (j = phm[k]-hmo)) {
        cf[k]   = cf[j];
        hmo[k]  = hmo[j];
        phm[k]  = &hmo[k];
        k       = j;
      }
      cf[k]   = tmpcf;
      hmo[k]  = tmphm;
      phm[k]  = &hmo[k];
    }
  }

  *cfp  = cf;
  *hmp  = hm;
}

void sort_terms_qq(
    mpz_t **cfp,
    hm_t **hmp,
    ht_t *ht
    )
{
  mpz_t *cf = *cfp;
  hm_t *hm  = *hmp;
  hm_t *hmo = hm+OFFSET;

  const len_t len = hm[LENGTH];

  len_t i, j, k;

  hm_t tmphm  = 0;
  mpz_t tmpcf;
  mpz_init(tmpcf);

  /* generate array of pointers to hm entries */
  hm_t *phm[len];
  for (i = 0; i < len; ++i) {
    phm[i]  = &hmo[i];
  }

  /* sort pointers to hm entries -> getting permutations */
  sort_r(phm, (unsigned long)len, sizeof(phm[0]), initial_gens_cmp, ht);

  /* sort cf and hm using permutations stored in phm */
  for (i = 0; i < len; ++i) {
    if (i != phm[i]-hmo) {
      mpz_swap(tmpcf, cf[i]);
      tmphm = hmo[i];
      k     = i;
      while (i != (j = phm[k]-hmo)) {
        mpz_swap(cf[k], cf[j]);
        hmo[k]  = hmo[j];
        phm[k]  = &hmo[k];
        k       = j;
      }
      mpz_swap(cf[k], tmpcf);
      hmo[k]  = tmphm;
      phm[k]  = &hmo[k];
    }
  }

  *cfp  = cf;
  *hmp  = hm;
}

void import_input_data(
        bs_t *bs,
        md_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs,
        const int *invalid_gens
        )
{
    int32_t i, j;
    len_t k;
    hm_t *hm;
    len_t ctr = 0; /* ctr for valid input elements */

    cf8_t *cf8      =   NULL;
    cf16_t *cf16    =   NULL;
    cf32_t *cf32    =   NULL;
    mpz_t *cfq      =   NULL;
    int32_t *cfs_ff =   NULL;
    mpz_t **cfs_qq  =   NULL;

    ht_t *ht = bs->ht;

    int32_t off       = 0; /* offset in arrays */
    int32_t init_off  = 0;
    const len_t fc    = st->fc;

    len_t ngens = stop - start;

    for (i = 0; i < start; ++i) {
        init_off +=  lens[i];
    }

    /* check basis size first */
    /* check_enlarge_basis(bs, ngens_input, st); */
    check_enlarge_basis(bs, ngens, st);

    /* import monomials */
    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    off = init_off;
    for (i = start; i < stop; ++i) {
        if (invalid_gens == NULL || invalid_gens[i] == 0) {
            while (lens[i] >= ht->esz-ht->eld) {
                enlarge_hash_table(ht);
                e  = ht->ev[0]; /* reset e if enlarging */
            }
            hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
            bs->hm[ctr] = hm;

            hm[COEFFS]  = ctr; /* link to matcf entry */
            hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
            hm[LENGTH]  = lens[i]; /* length */

            bs->red[ctr] = 0;

            for (j = off; j < off+lens[i]; ++j) {
                set_exponent_vector(e, exps, j, ht, st);
                hm[j-off+OFFSET]  =   insert_in_hash_table(e, ht);
            }
            ctr++;
        }
        off +=  lens[i];
    }
    /* import coefficients */
    off = init_off;
    ctr = 0;
    switch (st->ff_bits) {
        case 8:
            cfs_ff  =   (int32_t *)vcfs;
            for (i = start; i < stop; ++i) {
                if (invalid_gens == NULL || invalid_gens[i] == 0) {
                    cf8 = (cf8_t *)malloc((unsigned long)(lens[i]) * sizeof(cf8_t));
                    bs->cf_8[ctr] = cf8;

                    for (j = off; j < off+lens[i]; ++j) {
                        /* make coefficient positive */
                        cfs_ff[j]   +=  (cfs_ff[j] >> 31) & fc;
                        cf8[j-off]  =   (cf8_t)(cfs_ff[j] % fc);
                    }
                    sort_terms_ff_8(&(bs->cf_8[ctr]), &(bs->hm[ctr]), ht);
                    ctr++;
                }
                off +=  lens[i];
            }
            break;
        case 16:
            cfs_ff  =   (int32_t *)vcfs;
            for (i = start; i < stop; ++i) {
                if (invalid_gens == NULL || invalid_gens[i] == 0) {
                    cf16    = (cf16_t *)malloc((unsigned long)(lens[i]) * sizeof(cf16_t));
                    bs->cf_16[ctr] = cf16;

                    for (j = off; j < off+lens[i]; ++j) {
                        /* make coefficient positive */
                        cfs_ff[j]   +=  (cfs_ff[j] >> 31) & fc;
                        cf16[j-off] =   (cf16_t)(cfs_ff[j] % fc);
                    }
                    sort_terms_ff_16(&(bs->cf_16[ctr]), &(bs->hm[ctr]), ht);
                    ctr++;
                }
                off +=  lens[i];
            }
            break;
        case 32:
            cfs_ff  =   (int32_t *)vcfs;
            for (i = start; i < stop; ++i) {
                if (invalid_gens == NULL || invalid_gens[i] == 0) {
                    cf32    = (cf32_t *)malloc((unsigned long)(lens[i]) * sizeof(cf32_t));
                    bs->cf_32[ctr] = cf32;

                    for (j = off; j < off+lens[i]; ++j) {
                        /* make coefficient positive */
                        cfs_ff[j]   +=  (cfs_ff[j] >> 31) & fc;
                        cf32[j-off] =   (cf32_t)(cfs_ff[j] % fc);
                    }
                    sort_terms_ff_32(&(bs->cf_32[ctr]), &(bs->hm[ctr]), ht);
                    ctr++;
                }
                off +=  lens[i];
            }
            break;
        case 0:
            cfs_qq  =   (mpz_t **)vcfs;
            mpz_t prod_den, mul;
            mpz_inits(prod_den, mul, NULL);
            for (i = start; i < stop; ++i) {
                if (invalid_gens == NULL || invalid_gens[i] == 0) {
                    mpz_set_si(prod_den, 1);

                    for (j = off; j < off+lens[i]; ++j) {
                        /* printf("i %u | j %u\n", i, j);
                         * gmp_printf("%Zd\n", *(cfs[2*j+1])); */
                        mpz_mul(prod_den, prod_den, *(cfs_qq[2*j+1]));
                    }
                    cfq = (mpz_t *)malloc((unsigned long)(lens[i]) * sizeof(mpz_t));

                    bs->cf_qq[ctr]  = cfq;

                    for (j = 0; j < lens[i]; ++j) {
                        mpz_init(cfq[j]);
                    }
                    for (j = off; j < off+lens[i]; ++j) {
                        mpz_divexact(mul, prod_den, *(cfs_qq[2*j+1]));
                        mpz_mul(cfq[j-off], mul, *(cfs_qq[2*j]));
                    }
                    sort_terms_qq(&(bs->cf_qq[ctr]), &(bs->hm[ctr]), ht);
                    ctr++;
                }
                off +=  lens[i];
            }
            break;
        default:
            exit(1);
    }
    /* maybe some input elements are invalid, so reset ngens here */
    ngens = ctr;
            
    /* set total degree of input polynomials */
    deg_t deg = 0;
    if (st->nev) {
        for (i = 0; i < ngens; ++i) {
            hm  = bs->hm[i];
            deg = ht->hd[hm[OFFSET]].deg;
            k   = hm[LENGTH] + OFFSET;
            for (j = OFFSET+1; j < k; ++j) {
                if (deg < ht->hd[hm[j]].deg) {
                    deg = ht->hd[hm[j]].deg;
                    st->homogeneous = 1;
                }
            }
            bs->hm[i][DEG]  = deg;
        }
    } else {
        for (i = 0; i < ngens; ++i) {
            hm  = bs->hm[i];
            bs->hm[i][DEG]  = ht->hd[hm[OFFSET]].deg;
        }
    }
    if (st->homogeneous == 0) {
        /* check if input system is homogeneous or not */
        for (i = 0; i < ngens; ++i) {
            hm  = bs->hm[i];
            deg = ht->hd[hm[OFFSET]].deg;
            k   = hm[LENGTH] + OFFSET;
            for (j = OFFSET+1; j < k; ++j) {
                if (deg != ht->hd[hm[j]].deg) {
                    st->homogeneous = 0;
                    goto done;
                }
            }
        }
        st->homogeneous = 1;
    }
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;
}

/* return zero generator if all input generators are invalid */
static void return_zero(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const int32_t nr_vars,
        const uint32_t field_char,
        void *(*mallocp) (size_t)
        )
{

    int32_t *len  = (int32_t *)(*mallocp)(
            (unsigned long)1 * sizeof(int32_t));
    len[0]  =   1;
    int32_t *exp  = (int32_t *)(*mallocp)(
            (unsigned long)1 * (unsigned long)(nr_vars) * sizeof(int32_t));
    memset(exp, 0, (unsigned long)nr_vars * sizeof(int32_t));
    if (field_char > 0) {
    int32_t *cf   = (int32_t *)(*mallocp)(
            (unsigned long)1 * sizeof(int32_t));
    cf[0]   =   0;
    } else {
        fprintf(stderr, "We only support finite fields.\n");
    }
}

static int64_t export_data(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        void *(*mallocp) (size_t),
        const bs_t * const bs,
        const ht_t * const ht,
        const md_t * const md
        )
{
    len_t i, j, k;
    hm_t *dt;

    void *cf;
    mpz_t *tmp_cf_q;

    const len_t nv  = ht->nv;
    const len_t evl = ht->evl;
    const len_t ebl = ht->ebl;
    const len_t lml = bs->lml;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts   = 0; /* # elemnts in basis */

    for (i = 0; i < lml; ++i) {
        if (bs->hm[bs->lmps[i]] != NULL) {
        nterms +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
        } else {
            nterms++; /* allocate one term for 0 polynomial */
        }
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)(*mallocp)(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)(*mallocp)(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    if (md->ff_bits == 0) {
        cf = (mpz_t *)(*mallocp)(
            (unsigned long)(nterms) * sizeof(mpz_t));
    } else {
        cf = (int32_t *)(*mallocp)(
            (unsigned long)(nterms) * sizeof(int32_t));
    }

    /* counters for lengths, exponents and coefficients */
    int64_t cl = 0, ce = 0, cc = 0;

    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        /* polynomial is zero */
        if (bs->hm[bi] == NULL) {
            if (md->ff_bits == 0) {
                mpz_init(((mpz_t *)cf+cc)[0]);
            } else {
                ((int32_t *)cf+cc)[0] = (int32_t)0;
            }
            for (k = 1; k < evl; ++k) {
                exp[ce++] = (int32_t)0;
            }
            cc += 1;
            cl++;
            continue;
        }
        /* polynomial is nonzero */
        len[cl] = bs->hm[bi][LENGTH];
        switch (md->ff_bits) {
            case 8:
                for (j = 0; j < len[cl]; ++j) {
                    ((int32_t *)cf+cc)[j] = (int32_t)bs->cf_8[bs->hm[bi][COEFFS]][j];
                }
                break;
            case 16:
                for (j = 0; j < len[cl]; ++j) {
                    ((int32_t *)cf+cc)[j] = (int32_t)bs->cf_16[bs->hm[bi][COEFFS]][j];
                }
                break;
            case 32:
                for (j = 0; j < len[cl]; ++j) {
                    ((int32_t *)cf+cc)[j] = (int32_t)bs->cf_32[bs->hm[bi][COEFFS]][j];
                }
                break;
            case 0:
                tmp_cf_q =  bs->cf_qq[bs->hm[bi][COEFFS]];
                for (j = 0; j < len[cl]; ++j) {
                    mpz_init_set(((mpz_t *)cf+cc)[j], tmp_cf_q[j]);
                }
                break;
            default:
                exit(1);
        }
        dt  = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 1; k < ebl; ++k) {
                exp[ce++] = (int32_t)ht->ev[dt[j]][k];
            }
            for (k = ebl+1; k < evl; ++k) {
                exp[ce++] = (int32_t)ht->ev[dt[j]][k];
            }
        }
        cc  +=  len[cl];
        cl++;
    }

    *bload  = (int32_t)nelts;
    *blen   = len;
    *bexp   = exp;
    *bcf    = (void *)cf;

    return nterms;
}

void set_ff_bits(md_t *st, int32_t fc){
  if (fc == 0) {
    st->ff_bits = 0;
  } else {
    if (fc < pow(2,8)) {
      st->ff_bits = 8;
    } else {
      if (fc < pow(2,16)) {
        st->ff_bits = 16;
      } else {
        if (fc < pow(2,32)) {
          st->ff_bits = 32;
        }
      }
    }
  }
}

/* return 1 if validation was possible, zero otherwise */
int validate_input_data(
        int **invalid_gensp,
        const void *cfs,
        const int32_t *lens,
        uint32_t *field_charp,
        int32_t *mon_orderp,
        int32_t *elim_block_lenp,
        int32_t *nr_varsp,
        int32_t *nr_gensp,
        int32_t *nr_nfp,
        int32_t *ht_sizep,
        int32_t *nr_threadsp,
        int32_t *max_nr_pairsp,
        int32_t *reset_htp,
        int32_t *la_optionp,
        int32_t *use_signaturesp,
        int32_t *reduce_gbp,
        int32_t *info_levelp
        )
{
    /* biggest prime msovle can handle */
    if (*field_charp > 4294967291) {
        fprintf(stderr, "Field characteristic not valid.\n");
        return 0;
    }
    if (*nr_varsp < 0) {
        fprintf(stderr, "Number of variables not valid.\n");
        return 0;
    }
    if (*nr_gensp < 1) {
        fprintf(stderr, "Number of generators not valid.\n");
        return 0;
    }
    if (*nr_nfp < 0 || *nr_nfp >= *nr_gensp) {
        fprintf(stderr, "Number of normal forms not valid.\n");
        return 0;
    }
    if (*mon_orderp < 0) {
        fprintf(stderr, "Fixes monomial order to DRL.\n");
        *mon_orderp =   0;
    }
    if (*elim_block_lenp < 0) {
        fprintf(stderr, "Fixes elim block order length to 0.\n");
        *elim_block_lenp =   0;
    }
    if (*ht_sizep < 0) {
        fprintf(stderr, "Fixes initial hash table size to 2^17.\n");
        *ht_sizep   =   17;
    }
    if (*nr_threadsp < 0) {
        fprintf(stderr, "Fixes number of threads to 1.\n");
        *nr_threadsp    =   1;
    }
    if (*max_nr_pairsp < 0) {
        fprintf(stderr, "Fixes maximal number of spairs chosen to all possible.\n");
        *max_nr_pairsp  =   0;
    }
    if (*la_optionp != 1 && *la_optionp != 2 &&
            *la_optionp != 42 && *la_optionp != 44) {
        fprintf(stderr, "Fixes linear algebra option to exact sparse.\n");
        *la_optionp =   2;
    }
    if (*use_signaturesp < 0 || *use_signaturesp > 3) {
        fprintf(stderr, "Usage of signature not valid, disabled.\n");
        *use_signaturesp = 0;
    }
    if (*reduce_gbp < 0 || *reduce_gbp > 1) {
        fprintf(stderr, "Fixes reduction of GB to 0 (false).\n");
        *reduce_gbp =   0;
    }
    if (*info_levelp < 0 || *info_levelp > 2) {
        fprintf(stderr, "Fixes info level to no output.\n");
        *info_levelp    =   0;
    }
    const int ngens     =   *nr_gensp;
    int *invalid_gens   =   (int *)calloc((unsigned long)ngens, sizeof(int));
    int ctr             =   0;
    long len            =   0;
    if (*field_charp == 0) {
        mpz_t **cf  =   (mpz_t **)cfs;
        for (long i = 0; i < 2*len; ++i) {
            if (mpz_cmp_si(*(cf[0]), 0) == 0) {
                invalid_gens[i]   =   1;
                ctr++;
            }
        }
    } else {
        int32_t *cf =   (int32_t *)cfs;
        for (int i = 0; i < ngens; ++i) {
            for (int j = 0; j < lens[i]; ++j) {
                if (cf[j+len] == 0) {
                    invalid_gens[i]   =   1;
                    ctr++;
                    break;
                }
            }
            len +=  lens[i];
        }
    }

    *invalid_gensp  =   invalid_gens;

    if (ctr == 0) {
        return 1;
    }

    *nr_gensp   -=  ctr;

    /* recheck number of generators, if the input was only (0) then
     * no more generators are left over. */
    if (*nr_gensp < 1) {
        return -1;
    }

    return 1;
}

int32_t check_and_set_meta_data(
        md_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int *invalid_gens,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t nr_nf,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    if (nr_gens <= 0
            || nr_nf < 0
            || nr_vars <= 0
            || field_char < 0
            || use_signatures < 0
            || lens == NULL
            || cfs == NULL
            || exps == NULL) {
      fprintf(stderr, "Problem with meta data [%d, %d, %d]\n",
              (lens==NULL),(cfs==NULL),(exps==NULL));
      return 1;
    }

    long ctr    =   0;
    for (int i = 0; i < nr_gens; ++i) {
        ctr +=  invalid_gens[i];
    }

    /* number of generators given from input file */
    st->ngens_input     = nr_gens - nr_nf;
    /* number of generators from input which are invalid */
    st->ngens_invalid   = ctr;
    /* number of valid generators */
    st->ngens           = st->ngens_input - ctr;
    st->init_bs_sz      = 2 * nr_gens;

    st->nvars = nr_vars;
    /* note: prime check should be done in julia */
    st->fc    = field_char;

    set_ff_bits(st, st->fc);

    st->use_signatures  =   use_signatures;

    /* monomial order */
    if (mon_order != 0 && mon_order != 1) {
        st->mo  = 0;
    } else {
        st->mo  = mon_order;
    }
    /* elimination block order? If so, store the blocks length */
    st->nev = elim_block_len >= 0 ? elim_block_len : 0;
    if (st->nev >= st->nvars) {
        printf("error: Too large elimination block.\n");
        exit(1);
    }
    /* set hash table size */
    st->init_hts  = ht_size;
    if (st->init_hts <= 0) {
        st->init_hts  = 12;
    }
    /* info level */
    st->info_level  = info_level >= 0 ? info_level : 0;
    if (st->info_level > 2) {
        st->info_level = 2;
    }

    /* generation of pbm files on the fly? */
    st->gen_pbm_file  = pbm_file > 0 ? 1 : 0;

    /* resetting basis hash table */
    st->reset_ht = reset_hash_table > 0 ? reset_hash_table : 2147483647; /* 2^31-1 */;

    /* set number of threads */
    if (nr_threads <= 0) {
        st->nthrds  = 1;
    } else {
        st->nthrds  = nr_threads;
    }

    if (max_nr_pairs <= 0) {
        st->mnsel = 2147483647; /* 2^31-1 */
    } else {
        st->mnsel = max_nr_pairs;
    }

    /* set linear algebra option */
    if (la_option <= 0) {
        st->laopt = 1;
    } else {
        st->laopt = la_option;
    }

    if (reduce_gb < 0 || reduce_gb > 1) {
        st->reduce_gb = 0;
    } else {
        st->reduce_gb = reduce_gb;
    }
    set_function_pointers(st);

    return 0;
}

void set_function_pointers(
        const md_t *st
        )
{
  /* todo: this needs to be generalized for different monomial orders */
    if (st->nev > 0) {
      initial_input_cmp   = initial_input_cmp_be;
      initial_gens_cmp    = initial_gens_cmp_be;
      monomial_cmp        = monomial_cmp_be;
      spair_cmp           = spair_cmp_be;
      hcm_cmp             = hcm_cmp_pivots_be;
    } else {
        switch (st->mo) {
            case 0:
                initial_input_cmp   = initial_input_cmp_drl;
                initial_gens_cmp    = initial_gens_cmp_drl;
                monomial_cmp        = monomial_cmp_drl;
                spair_cmp           = spair_cmp_drl;
                hcm_cmp             = hcm_cmp_pivots_drl;
                break;
            case 1:
                initial_input_cmp   = initial_input_cmp_lex;
                initial_gens_cmp    = initial_gens_cmp_lex;
                monomial_cmp        = monomial_cmp_lex;
                spair_cmp           = spair_cmp_deglex;
                hcm_cmp             = hcm_cmp_pivots_lex;
                break;
            default:
                initial_input_cmp   = initial_input_cmp_drl;
                initial_gens_cmp    = initial_gens_cmp_drl;
                monomial_cmp        = monomial_cmp_drl;
                spair_cmp           = spair_cmp_drl;
                hcm_cmp             = hcm_cmp_pivots_drl;
        }
    }

  /* up to 17 bits we can use one modular operation for reducing a row. this works
   * for matrices with #rows <= 54 million */
  switch (st->ff_bits) {
    case 0:
      switch (st->laopt) {
        case 1:
          linear_algebra  = exact_sparse_linear_algebra_ab_first_qq;
          break;
        case 2:
          linear_algebra  = exact_sparse_linear_algebra_qq;
          break;
        default:
          linear_algebra  = exact_sparse_linear_algebra_qq;
      }
      interreduce_matrix_rows = interreduce_matrix_rows_qq;
      break;

    case 8:
      switch (st->laopt) {
        case 1:
          linear_algebra  = exact_sparse_dense_linear_algebra_ff_8;
          break;
        case 2:
          linear_algebra  = exact_sparse_linear_algebra_ff_8;
          break;
        case 42:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_8;
          break;
        case 43:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_8_2;
          break;
        case 44:
          linear_algebra  = probabilistic_sparse_linear_algebra_ff_8;
          break;
        default:
          linear_algebra  = exact_sparse_linear_algebra_ff_8;
      }
      interreduce_matrix_rows     = interreduce_matrix_rows_ff_8;
      normalize_initial_basis     = normalize_initial_basis_ff_8;
      break;

    case 16:
      switch (st->laopt) {
        case 1:
          linear_algebra  = exact_sparse_dense_linear_algebra_ff_16;
          break;
        case 2:
          linear_algebra  = exact_sparse_linear_algebra_ff_16;
          break;
        case 42:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_16;
          break;
        case 43:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_16_2;
          break;
        case 44:
          linear_algebra  = probabilistic_sparse_linear_algebra_ff_16;
          break;
        default:
          linear_algebra  = exact_sparse_linear_algebra_ff_16;
      }
      interreduce_matrix_rows     = interreduce_matrix_rows_ff_16;
      normalize_initial_basis     = normalize_initial_basis_ff_16;
      break;

    case 32:
      switch (st->laopt) {
        case 1:
          linear_algebra  = exact_sparse_dense_linear_algebra_ff_32;
          break;
        case 2:
          linear_algebra  = exact_sparse_linear_algebra_ff_32;
          break;
        case 42:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32;
          break;
        case 43:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32_2;
          break;
        case 44:
          linear_algebra  = probabilistic_sparse_linear_algebra_ff_32;
          break;
        default:
          linear_algebra  = exact_sparse_linear_algebra_ff_32;
      }
      interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
      normalize_initial_basis     = normalize_initial_basis_ff_32;
      sba_linear_algebra          = sba_linear_algebra_ff_32;

      sba_reduce_dense_row_by_known_pivots_sparse_ff_32 =
        sba_reduce_dense_row_by_known_pivots_sparse_31_bit;
      /* if coeffs are smaller than 17 bit we can optimize reductions */
      if (st->fc < pow(2, 18)) {
        reduce_dense_row_by_all_pivots_ff_32 =
          reduce_dense_row_by_all_pivots_17_bit;
        reduce_dense_row_by_old_pivots_ff_32 =
          reduce_dense_row_by_old_pivots_17_bit;
        reduce_dense_row_by_known_pivots_sparse_ff_32 =
          reduce_dense_row_by_known_pivots_sparse_17_bit;
        reduce_dense_row_by_dense_new_pivots_ff_32  =
          reduce_dense_row_by_dense_new_pivots_17_bit;
      } else {
        if (st->fc < pow(2, 31)) {
          reduce_dense_row_by_all_pivots_ff_32 =
            reduce_dense_row_by_all_pivots_31_bit;
          reduce_dense_row_by_old_pivots_ff_32 =
            reduce_dense_row_by_old_pivots_31_bit;
          reduce_dense_row_by_known_pivots_sparse_ff_32 =
            reduce_dense_row_by_known_pivots_sparse_31_bit;
          reduce_dense_row_by_dense_new_pivots_ff_32  =
            reduce_dense_row_by_dense_new_pivots_31_bit;
        } else {
          reduce_dense_row_by_all_pivots_ff_32 =
            reduce_dense_row_by_all_pivots_31_bit;
          reduce_dense_row_by_old_pivots_ff_32 =
            reduce_dense_row_by_old_pivots_31_bit;
          reduce_dense_row_by_known_pivots_sparse_ff_32 =
            reduce_dense_row_by_known_pivots_sparse_32_bit;
          reduce_dense_row_by_dense_new_pivots_ff_32  =
            reduce_dense_row_by_dense_new_pivots_31_bit;
        }
      }
      break;

    default:
      switch (st->laopt) {
        case 1:
          linear_algebra  = exact_sparse_dense_linear_algebra_ff_32;
          break;
        case 2:
          linear_algebra  = exact_sparse_linear_algebra_ff_32;
          break;
        case 42:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32;
          break;
        case 43:
          linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32_2;
          break;
        case 44:
          linear_algebra  = probabilistic_sparse_linear_algebra_ff_32;
          break;
        default:
          linear_algebra  = exact_sparse_linear_algebra_ff_32;
      }
      interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
      normalize_initial_basis     = normalize_initial_basis_ff_32;

      /* if coeffs are smaller than 17 bit we can optimize reductions */
      if (st->fc < pow(2, 18)) {
        reduce_dense_row_by_all_pivots_ff_32 =
          reduce_dense_row_by_all_pivots_17_bit;
        reduce_dense_row_by_old_pivots_ff_32 =
          reduce_dense_row_by_old_pivots_17_bit;
        reduce_dense_row_by_known_pivots_sparse_ff_32 =
          reduce_dense_row_by_known_pivots_sparse_17_bit;
        reduce_dense_row_by_dense_new_pivots_ff_32  =
          reduce_dense_row_by_dense_new_pivots_17_bit;
      } else {
        if (st->fc < pow(2, 31)) {
          reduce_dense_row_by_all_pivots_ff_32 =
            reduce_dense_row_by_all_pivots_31_bit;
          reduce_dense_row_by_old_pivots_ff_32 =
            reduce_dense_row_by_old_pivots_31_bit;
          reduce_dense_row_by_known_pivots_sparse_ff_32 =
            reduce_dense_row_by_known_pivots_sparse_31_bit;
          reduce_dense_row_by_dense_new_pivots_ff_32  =
            reduce_dense_row_by_dense_new_pivots_31_bit;
        } else {
          reduce_dense_row_by_all_pivots_ff_32 =
            reduce_dense_row_by_all_pivots_31_bit;
          reduce_dense_row_by_old_pivots_ff_32 =
            reduce_dense_row_by_old_pivots_31_bit;
          reduce_dense_row_by_known_pivots_sparse_ff_32 =
            reduce_dense_row_by_known_pivots_sparse_32_bit;
          reduce_dense_row_by_dense_new_pivots_ff_32  =
            reduce_dense_row_by_dense_new_pivots_31_bit;
        }
      }
  }
}

int32_t check_and_set_meta_data_trace(
        md_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int *invalid_gens,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t elim_block_len,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t nr_nf,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t use_signatures,
        const int32_t reduce_gb,
        const uint32_t prime_start,
        const int32_t nr_primes,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    st->prime_start = prime_start;
    if (st->prime_start <= 0) {
        st->prime_start = 32003;
    }
    st->nprimes = nr_primes;
    if (st->nprimes <= 0) {
        st->nprimes = 10;
    }
    return check_and_set_meta_data(st, lens, exps, cfs, invalid_gens,
            field_char, mon_order, elim_block_len, nr_vars, nr_gens,
            nr_nf, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
            la_option, use_signatures, reduce_gb, pbm_file, info_level);
}

static inline void reset_function_pointers(
        const uint32_t prime,
        const uint32_t laopt
        )
{
    if (prime < pow(2,8)) {
        interreduce_matrix_rows     = interreduce_matrix_rows_ff_8;
        normalize_initial_basis     = normalize_initial_basis_ff_8;
        switch (laopt) {
          case 1:
            linear_algebra  = exact_sparse_dense_linear_algebra_ff_8;
            break;
          case 2:
            linear_algebra  = exact_sparse_linear_algebra_ff_8;
            break;
          case 42:
            linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_8;
            break;
          case 43:
            linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_8_2;
            break;
          case 44:
            linear_algebra  = probabilistic_sparse_linear_algebra_ff_8;
            break;
          default:
            linear_algebra  = exact_sparse_linear_algebra_ff_8;
        }
    } else {
        if (prime < pow(2,16)) {
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_16;
            normalize_initial_basis     = normalize_initial_basis_ff_16;
            switch (laopt) {
              case 1:
                linear_algebra  = exact_sparse_dense_linear_algebra_ff_16;
                break;
              case 2:
                linear_algebra  = exact_sparse_linear_algebra_ff_16;
                break;
              case 42:
                linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_16;
                break;
              case 43:
                linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_16_2;
                break;
              case 44:
                linear_algebra  = probabilistic_sparse_linear_algebra_ff_16;
                break;
              default:
                linear_algebra  = exact_sparse_linear_algebra_ff_16;
            }
        } else {
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
            normalize_initial_basis     = normalize_initial_basis_ff_32;
            switch (laopt) {
              case 1:
                linear_algebra  = exact_sparse_dense_linear_algebra_ff_32;
                break;
              case 2:
                linear_algebra  = exact_sparse_linear_algebra_ff_32;
                break;
              case 42:
                linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32;
                break;
              case 43:
                linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_32_2;
                break;
              case 44:
                linear_algebra  = probabilistic_sparse_linear_algebra_ff_32;
                break;
              default:
                linear_algebra  = exact_sparse_linear_algebra_ff_32;
            }
            if (prime < pow(2,18)) {
                reduce_dense_row_by_all_pivots_ff_32 =
                    reduce_dense_row_by_all_pivots_17_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                    reduce_dense_row_by_old_pivots_17_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                    reduce_dense_row_by_known_pivots_sparse_17_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                    reduce_dense_row_by_dense_new_pivots_17_bit;
            } else {
              if (prime < pow(2,31)) {
                reduce_dense_row_by_all_pivots_ff_32 =
                  reduce_dense_row_by_all_pivots_31_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                  reduce_dense_row_by_old_pivots_31_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  reduce_dense_row_by_known_pivots_sparse_31_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                  reduce_dense_row_by_dense_new_pivots_31_bit;
              } else {
                reduce_dense_row_by_all_pivots_ff_32 =
                  reduce_dense_row_by_all_pivots_31_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                  reduce_dense_row_by_old_pivots_31_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  reduce_dense_row_by_known_pivots_sparse_32_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                  reduce_dense_row_by_dense_new_pivots_31_bit;
              }
            }
        }
    }

}
static inline void reset_trace_function_pointers(
        const uint32_t prime
        )
{
    if (prime < pow(2,8)) {
        interreduce_matrix_rows     = interreduce_matrix_rows_ff_8;
        normalize_initial_basis     = normalize_initial_basis_ff_8;
        application_linear_algebra  = exact_application_sparse_linear_algebra_ff_8;
        trace_linear_algebra        = exact_trace_sparse_linear_algebra_ff_8;
    } else {
        if (prime < pow(2,16)) {
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_16;
            normalize_initial_basis     = normalize_initial_basis_ff_16;
            application_linear_algebra  = exact_application_sparse_linear_algebra_ff_16;
            trace_linear_algebra        = exact_trace_sparse_linear_algebra_ff_16;
        } else {
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
            normalize_initial_basis     = normalize_initial_basis_ff_32;
            application_linear_algebra  = exact_application_sparse_linear_algebra_ff_32;
            trace_linear_algebra        = exact_trace_sparse_linear_algebra_ff_32;
            if (prime < pow(2,18)) {
                reduce_dense_row_by_all_pivots_ff_32 =
                    reduce_dense_row_by_all_pivots_17_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                    reduce_dense_row_by_old_pivots_17_bit;
                trace_reduce_dense_row_by_known_pivots_sparse_ff_32 =
                    trace_reduce_dense_row_by_known_pivots_sparse_17_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                    reduce_dense_row_by_known_pivots_sparse_17_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                    reduce_dense_row_by_dense_new_pivots_17_bit;
            } else {
              if (prime < pow(2,31)) {
                reduce_dense_row_by_all_pivots_ff_32 =
                  reduce_dense_row_by_all_pivots_31_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                  reduce_dense_row_by_old_pivots_31_bit;
                trace_reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  trace_reduce_dense_row_by_known_pivots_sparse_31_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  reduce_dense_row_by_known_pivots_sparse_31_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                  reduce_dense_row_by_dense_new_pivots_31_bit;
              } else {
                reduce_dense_row_by_all_pivots_ff_32 =
                  reduce_dense_row_by_all_pivots_31_bit;
                reduce_dense_row_by_old_pivots_ff_32 =
                  reduce_dense_row_by_old_pivots_31_bit;
                trace_reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  trace_reduce_dense_row_by_known_pivots_sparse_32_bit;
                reduce_dense_row_by_known_pivots_sparse_ff_32 =
                  reduce_dense_row_by_known_pivots_sparse_32_bit;
                reduce_dense_row_by_dense_new_pivots_ff_32  =
                  reduce_dense_row_by_dense_new_pivots_31_bit;
              }
            }
        }
    }

}

static void write_pbm_file(
    mat_t *mat,
    const md_t * const st
    )
{
    len_t i, j, k;
    unsigned char b = 0;
    char buffer[512];
    char fn[200];

    hm_t **rows = mat->rr;
    const len_t ncols = mat->nc;
    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;

    snprintf(fn, 200, "%d-%d-%d-%d.pbm", st->current_rd, nru+nrl, ncols, st->current_deg);
    FILE *fh  = fopen(fn, "wb");

    /* magic header */
    snprintf(buffer, 512, "P4\n# matrix size(%u, %u)\n%u %u\n", nru+nrl, ncols, ncols, nru+nrl);

    fwrite(buffer, sizeof(char), strlen(buffer), fh);


    for (i = 0; i < nru; ++i) {
        const len_t len = rows[i][LENGTH];
        hm_t row[len];
        memcpy(row, rows[i]+OFFSET, (unsigned long)len * sizeof(hm_t));
        qsort(row, (unsigned long)len, sizeof(hm_t), pbm_cmp);
        /* the rows may not be sorted by column indices, thus we
         * have to go over them again and again and again */
        k = 0;
        for (j = 0; j < ncols; ++j) {
            if (k < len && row[k] == j) {
                b |=  (unsigned char)(1 << (7 - (j % 8)));
                k++;
            } else {
                b &= (unsigned char) ~(1 << (7 - (j % 8)));
            }
            /* write each byte */
            if (j % 8 == 7) {
                fwrite(&b, sizeof(unsigned char), 1, fh);
                b = 0;
            }
        }
        /* write leftover bits */
        if (j % 8 != 0) {
            fwrite(&b, sizeof(unsigned char), 1, fh);
        }
        fflush(fh);
    }
    rows  = mat->tr;
    for (i = 0; i < nrl; ++i) {
        const len_t len = rows[i][LENGTH];
        hm_t row[len];
        memcpy(row, rows[i]+OFFSET, (unsigned long)len * sizeof(hm_t));
        qsort(row, (unsigned long)len, sizeof(hm_t), pbm_cmp);
        /* the rows may not be sorted by column indices, thus we
         * have to go over them again and again and again */
        k = 0;
        for (j = 0; j < ncols; ++j) {
            if (k < len && row[k] == j) {
                b |=  (unsigned char)(1 << (7 - (j % 8)));
                k++;
            } else {
                b &= (unsigned char) ~(1 << (7 - (j % 8)));
            }
            /* write each byte */
            if (j % 8 == 7) {
                fwrite(&b, sizeof(unsigned char), 1, fh);
                b = 0;
            }
        }
        /* write leftover bits */
        if (j % 8 != 0) {
            fwrite(&b, sizeof(unsigned char), 1, fh);
        }
        fflush(fh);
    }
    fclose(fh);
}
