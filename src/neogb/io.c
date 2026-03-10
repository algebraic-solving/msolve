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
#include "../msolve/streams.h"

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
            mpz_clears(prod_den, mul, NULL);
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
void return_zero(
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
        fprintf(ERRSTREAM, "We only support finite fields.\n");
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
        fprintf(ERRSTREAM, "Basis has more than 2^31 elements, cannot store it.\n");
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
            len[cl] = 1;
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

int32_t check_ff_bits(int32_t fc){
    if (fc == 0) {
        return 0;
    } else {
        if (fc < (int32_t)(1) << 8) {
            return 8;
        } else {
            if (fc < (int32_t)(1) << 16) {
                return 16;
            } else {
                if (fc < (int32_t)(1) << 23) {
                    return 32;
                } else {
                    return 32;
                }
            }
        }
    }
}

void set_ff_bits(md_t *st, int32_t fc){
    if (fc == 0) {
        st->ff_bits = 0;
    } else {
        if (fc < (int32_t)(1) << 8) {
            st->ff_bits = 8;
        } else {
            if (fc < (int32_t)(1) << 16) {
                st->ff_bits = 16;
            } else {
                if (fc < (int32_t)(1) << 23) {
                    st->ff_bits = 32;
                } else {
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
        int32_t *truncate_liftingp,
        int32_t *info_levelp
        )
{
    /* biggest prime msovle can handle */
    if (*field_charp > 4294967291) {
        fprintf(ERRSTREAM, "Field characteristic not valid.\n");
        return 0;
    }
    if (*nr_varsp < 0) {
        fprintf(ERRSTREAM, "Number of variables not valid.\n");
        return 0;
    }
    if (*nr_gensp < 1) {
        fprintf(ERRSTREAM, "Number of generators not valid.\n");
        return 0;
    }
    if (*nr_nfp < 0 || *nr_nfp >= *nr_gensp) {
        fprintf(ERRSTREAM, "Number of normal forms not valid.\n");
        return 0;
    }
    if (*mon_orderp < 0) {
        fprintf(ERRSTREAM, "Fixes monomial order to DRL.\n");
        *mon_orderp =   0;
    }
    if (*elim_block_lenp < 0) {
        fprintf(ERRSTREAM, "Fixes elim block order length to 0.\n");
        *elim_block_lenp =   0;
    }
    if (*ht_sizep < 0) {
        fprintf(ERRSTREAM, "Fixes initial hash table size to 2^17.\n");
        *ht_sizep   =   17;
    }
    if (*nr_threadsp < 0) {
        fprintf(ERRSTREAM, "Fixes number of threads to 1.\n");
        *nr_threadsp    =   1;
    }
    if (*max_nr_pairsp < 0) {
        fprintf(ERRSTREAM, "Fixes maximal number of spairs chosen to all possible.\n");
        *max_nr_pairsp  =   0;
    }
    if (*la_optionp != 1 && *la_optionp != 2 &&
            *la_optionp != 42 && *la_optionp != 44) {
        fprintf(ERRSTREAM, "Fixes linear algebra option to exact sparse.\n");
        *la_optionp =   2;
    }
    if (*use_signaturesp < 0 || *use_signaturesp > 3) {
        fprintf(ERRSTREAM, "Usage of signature not valid, disabled.\n");
        *use_signaturesp = 0;
    }
    if (*reduce_gbp < 0 || *reduce_gbp > 1) {
        fprintf(ERRSTREAM, "Fixes reduction of GB to 0 (false).\n");
        *reduce_gbp =   0;
    }
    if(*truncate_liftingp < 0){
        fprintf(ERRSTREAM, "Removes truncation of lifted Groebner bases\n");
        *truncate_liftingp = 0;
    }
    if (*info_levelp < 0 || *info_levelp > 2) {
        fprintf(ERRSTREAM, "Fixes info level to no output.\n");
        *info_levelp    =   0;
    }
    const int ngens     =   *nr_gensp;
    int *invalid_gens   =   (int *)calloc((unsigned long)ngens, sizeof(int));
    int ctr             =   0;
    long len            =   0;
    if (*field_charp == 0) {
        mpz_t **cf  =   (mpz_t **)cfs;
        for (int i = 0; i < ngens; ++i) {
            for (int j = 0; j < 2*lens[i]; ++j) {
                if (mpz_cmp_si(*(cf[j+len]), 0) == 0) {
                    invalid_gens[i]   =   1;
                    ctr++;
                    break;
                }
            }
            len += 2*lens[i];
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

static inline int use_block_order(
        const ht_t *ht
        )
{
    return ht->ebl > 0;
}

static inline int use_lex_order(
        const ht_t *ht
        )
{
    return ht->ebl == 0 && ht->mo == 1;
}

static inline int use_17_bit_reduction(
        const uint32_t fc
        )
{
    return fc < ((uint32_t)1u << 18);
}

int dispatch_initial_input_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;
    if (use_block_order(ht)) {
        return initial_input_cmp_be(a, b, htp);
    }
    if (use_lex_order(ht)) {
        return initial_input_cmp_lex(a, b, htp);
    }
    return initial_input_cmp_drl(a, b, htp);
}

int dispatch_initial_gens_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;
    if (use_block_order(ht)) {
        return initial_gens_cmp_be(a, b, htp);
    }
    if (use_lex_order(ht)) {
        return initial_gens_cmp_lex(a, b, htp);
    }
    return initial_gens_cmp_drl(a, b, htp);
}

int dispatch_monomial_cmp(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        )
{
    if (use_block_order(ht)) {
        return monomial_cmp_be(a, b, ht);
    }
    if (use_lex_order(ht)) {
        return monomial_cmp_lex(a, b, ht);
    }
    return monomial_cmp_drl(a, b, ht);
}

int dispatch_spair_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;
    if (use_block_order(ht)) {
        return spair_cmp_be(a, b, htp);
    }
    if (use_lex_order(ht)) {
        return spair_cmp_deglex(a, b, htp);
    }
    return spair_cmp_drl(a, b, htp);
}

int dispatch_hcm_cmp(
        const void *a,
        const void *b,
        void *htp
        )
{
    ht_t *ht = (ht_t *)htp;
    if (use_block_order(ht)) {
        return hcm_cmp_pivots_be(a, b, htp);
    }
    if (use_lex_order(ht)) {
        return hcm_cmp_pivots_lex(a, b, htp);
    }
    return hcm_cmp_pivots_drl(a, b, htp);
}

void dispatch_sba_linear_algebra(
        smat_t *smat,
        crit_t *syz,
        md_t *st,
        const ht_t * const ht
        )
{
    sba_linear_algebra_ff_32(smat, syz, st, ht);
}

void dispatch_linear_algebra(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    switch (st->ff_bits) {
        case 0:
            if (st->laopt == 1) {
                exact_sparse_linear_algebra_ab_first_qq(mat, tbr, bs, st);
            } else {
                exact_sparse_linear_algebra_qq(mat, tbr, bs, st);
            }
            return;
        case 8:
            switch (st->laopt) {
                case 1:
                    exact_sparse_dense_linear_algebra_ff_8(mat, tbr, bs, st);
                    return;
                case 2:
                    exact_sparse_linear_algebra_ff_8(mat, tbr, bs, st);
                    return;
                case 42:
                    probabilistic_sparse_dense_linear_algebra_ff_8(mat, tbr, bs, st);
                    return;
                case 43:
                    probabilistic_sparse_dense_linear_algebra_ff_8_2(mat, tbr, bs, st);
                    return;
                case 44:
                    probabilistic_sparse_linear_algebra_ff_8(mat, tbr, bs, st);
                    return;
                default:
                    exact_sparse_linear_algebra_ff_8(mat, tbr, bs, st);
                    return;
            }
        case 16:
            switch (st->laopt) {
                case 1:
                    exact_sparse_dense_linear_algebra_ff_16(mat, tbr, bs, st);
                    return;
                case 2:
                    exact_sparse_linear_algebra_ff_16(mat, tbr, bs, st);
                    return;
                case 42:
                    probabilistic_sparse_dense_linear_algebra_ff_16(mat, tbr, bs, st);
                    return;
                case 43:
                    probabilistic_sparse_dense_linear_algebra_ff_16_2(mat, tbr, bs, st);
                    return;
                case 44:
                    probabilistic_sparse_linear_algebra_ff_16(mat, tbr, bs, st);
                    return;
                default:
                    exact_sparse_linear_algebra_ff_16(mat, tbr, bs, st);
                    return;
            }
        case 32:
        default:
            switch (st->laopt) {
                case 1:
                    exact_sparse_dense_linear_algebra_ff_32(mat, tbr, bs, st);
                    return;
                case 2:
                    exact_sparse_linear_algebra_ff_32(mat, tbr, bs, st);
                    return;
                case 42:
                    probabilistic_sparse_dense_linear_algebra_ff_32(mat, tbr, bs, st);
                    return;
                case 43:
                    probabilistic_sparse_dense_linear_algebra_ff_32_2(mat, tbr, bs, st);
                    return;
                case 44:
                    probabilistic_sparse_linear_algebra_ff_32(mat, tbr, bs, st);
                    return;
                default:
                    exact_sparse_linear_algebra_ff_32(mat, tbr, bs, st);
                    return;
            }
    }
}

void dispatch_exact_linear_algebra(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        )
{
    switch (st->ff_bits) {
        case 0:
            exact_sparse_linear_algebra_qq(mat, tbr, bs, st);
            return;
        case 8:
            exact_sparse_linear_algebra_ff_8(mat, tbr, bs, st);
            return;
        case 16:
            exact_sparse_linear_algebra_ff_16(mat, tbr, bs, st);
            return;
        case 32:
        default:
            exact_sparse_linear_algebra_ff_32(mat, tbr, bs, st);
            return;
    }
}

int dispatch_application_linear_algebra(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    switch (st->ff_bits) {
        case 8:
            return exact_application_sparse_linear_algebra_ff_8(mat, bs, st);
        case 16:
            return exact_application_sparse_linear_algebra_ff_16(mat, bs, st);
        case 32:
        default:
            return exact_application_sparse_linear_algebra_ff_32(mat, bs, st);
    }
}

void dispatch_trace_linear_algebra(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        )
{
    switch (st->ff_bits) {
        case 8:
            exact_trace_sparse_linear_algebra_ff_8(trace, mat, bs, st);
            return;
        case 16:
            exact_trace_sparse_linear_algebra_ff_16(trace, mat, bs, st);
            return;
        case 32:
        default:
            exact_trace_sparse_linear_algebra_ff_32(trace, mat, bs, st);
            return;
    }
}

void dispatch_interreduce_matrix_rows(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        int free_basis
        )
{
    switch (st->ff_bits) {
        case 0:
            interreduce_matrix_rows_qq(mat, bs, st, free_basis);
            return;
        case 8:
            interreduce_matrix_rows_ff_8(mat, bs, st, free_basis);
            return;
        case 16:
            interreduce_matrix_rows_ff_16(mat, bs, st, free_basis);
            return;
        case 32:
        default:
            interreduce_matrix_rows_ff_32(mat, bs, st, free_basis);
            return;
    }
}

void dispatch_normalize_initial_basis(
        bs_t *bs,
        const uint32_t fc
        )
{
    if (fc == 0) {
        return;
    }
    if (fc < ((uint32_t)1u << 8)) {
        normalize_initial_basis_ff_8(bs, fc);
        return;
    }
    if (fc < ((uint32_t)1u << 16)) {
        normalize_initial_basis_ff_16(bs, fc);
        return;
    }
    normalize_initial_basis_ff_32(bs, fc);
}

cf32_t *dispatch_reduce_dense_row_by_old_pivots_ff_32(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,
        const uint32_t fc
        )
{
    if (use_17_bit_reduction(fc)) {
        return reduce_dense_row_by_old_pivots_17_bit(dr, mat, bs, pivs, dpiv, fc);
    }
    return reduce_dense_row_by_old_pivots_31_bit(dr, mat, bs, pivs, dpiv, fc);
}

hm_t *dispatch_sba_reduce_dense_row_by_known_pivots_sparse_ff_32(
        int64_t *dr,
        smat_t *smat,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t sm,
        const len_t si,
        const len_t ri,
        md_t *st
        )
{
    return sba_reduce_dense_row_by_known_pivots_sparse_31_bit(
            dr, smat, pivs, dpiv, sm, si, ri, st);
}

hm_t *dispatch_reduce_dense_row_by_known_pivots_sparse_ff_32(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi,
        const len_t tr,
        md_t *st
        )
{
    if (use_17_bit_reduction(st->fc)) {
        return reduce_dense_row_by_known_pivots_sparse_17_bit(
                dr, mat, bs, pivs, dpiv, tmp_pos, mh, bi, tr, st);
    }
    return reduce_dense_row_by_known_pivots_sparse_31_bit(
            dr, mat, bs, pivs, dpiv, tmp_pos, mh, bi, tr, st);
}

hm_t *dispatch_trace_reduce_dense_row_by_known_pivots_sparse_ff_32(
        rba_t *rba,
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        const len_t mh,
        const len_t bi,
        md_t *st
        )
{
    if (use_17_bit_reduction(st->fc)) {
        return trace_reduce_dense_row_by_known_pivots_sparse_17_bit(
                rba, dr, mat, bs, pivs, dpiv, tmp_pos, mh, bi, st);
    }
    return trace_reduce_dense_row_by_known_pivots_sparse_31_bit(
            rba, dr, mat, bs, pivs, dpiv, tmp_pos, mh, bi, st);
}

cf32_t *dispatch_reduce_dense_row_by_all_pivots_ff_32(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const *pivs,
        cf32_t *const *dpivs,
        const uint32_t fc
        )
{
    if (use_17_bit_reduction(fc)) {
        return reduce_dense_row_by_all_pivots_17_bit(dr, mat, bs, pc, pivs, dpivs, fc);
    }
    return reduce_dense_row_by_all_pivots_31_bit(dr, mat, bs, pc, pivs, dpivs, fc);
}

cf32_t *dispatch_reduce_dense_row_by_dense_new_pivots_ff_32(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        )
{
    if (use_17_bit_reduction(fc)) {
        return reduce_dense_row_by_dense_new_pivots_17_bit(dr, pc, pivs, ncr, fc);
    }
    return reduce_dense_row_by_dense_new_pivots_31_bit(dr, pc, pivs, ncr, fc);
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
        const int32_t truncate_lifting,
        const int32_t info_level
        )
{
    if (nr_gens <= 0
            || nr_nf < 0
            || nr_vars <= 0
            || use_signatures < 0
            || lens == NULL
            || cfs == NULL
            || exps == NULL) {
      fprintf(ERRSTREAM, "Problem with meta data [%d, %d, %d]\n",
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
        fprintf(ERRSTREAM,"error: Too large elimination block.\n");
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

    st->truncate_lifting = truncate_lifting >= 0 ? truncate_lifting : 0;

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
    (void)st;
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
        const int32_t truncate_lifting,
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
            la_option, use_signatures, reduce_gb, pbm_file, truncate_lifting,
            info_level);
}

static inline void reset_function_pointers(
        const uint32_t prime,
        const uint32_t laopt
        )
{
    (void)prime;
    (void)laopt;
}
static inline void reset_trace_function_pointers(
        const uint32_t prime
        )
{
    (void)prime;
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
