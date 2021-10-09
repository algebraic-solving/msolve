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

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
static void sort_terms_ff_8(
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

static void import_julia_data_ff_8(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    cf8_t *cf;
    hm_t *hm;

    int32_t *cfs  = (int32_t *)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;
    const len_t fc    = st->fc;

    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = 0; i < ngens; ++i) {
        while (lens[i] >= ht->esz-ht->eld) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (cf8_t *)malloc((unsigned long)(lens[i]) * sizeof(cf8_t));
        bs->hm[i]   = hm;
        bs->cf_8[i] = cf;

        hm[COEFFS]  = i; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        bs->red[i] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET]  =   insert_in_hash_table(e, ht);
            /* make coefficient positive */
            cfs[j]            +=  (cfs[j] >> 31) & fc;
            cf[j-off]         =   (cf8_t)cfs[j];
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_ff_8(&cf, &hm, ht);
    }
    deg_t deg = 0;
    for (i = 0; i < ngens; ++i) {
        hm  = bs->hm[i];
        deg = ht->ev[hm[OFFSET]][DEG];
        k   = hm[LENGTH] + OFFSET;
        for (j = OFFSET+1; j < k; ++j) {
            if (deg != ht->ev[hm[j]][DEG]) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
    st->homogeneous = 1;
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;
}

static void sort_terms_ff_16(
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

static void import_julia_data_ff_16(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    cf16_t *cf;
    hm_t *hm;

    int32_t *cfs  = (int32_t *)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;
    const len_t fc    = st->fc;

    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = 0; i < ngens; ++i) {
        while (lens[i] >= ht->esz-ht->eld) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (cf16_t *)malloc((unsigned long)(lens[i]) * sizeof(cf16_t));
        bs->hm[i]     = hm;
        bs->cf_16[i]  = cf;

        hm[COEFFS]  = i; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        bs->red[i] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET]  =   insert_in_hash_table(e, ht);
            /* make coefficient positive */
            cfs[j]            +=  (cfs[j] >> 31) & fc;
            cf[j-off]         =   (cf16_t)cfs[j];
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_ff_16(&cf, &hm, ht);
    }
    deg_t deg = 0;
    for (i = 0; i < ngens; ++i) {
        hm  = bs->hm[i];
        deg = ht->ev[hm[OFFSET]][DEG];
        k   = hm[LENGTH] + OFFSET;
        for (j = OFFSET+1; j < k; ++j) {
            if (deg != ht->ev[hm[j]][DEG]) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
    st->homogeneous = 1;
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;
}

static void sort_terms_ff_32(
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

static void import_julia_data_ff_32(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    cf32_t *cf    = NULL;
    int64_t tmpcf = 0;
    hm_t *hm      = NULL;

    int32_t *cfs  = (int32_t *)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;
    const len_t fc    = st->fc;

    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = 0; i < ngens; ++i) {
        while (lens[i] >= ht->esz-ht->eld) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (cf32_t *)malloc((unsigned long)(lens[i]) * sizeof(cf32_t));
        bs->hm[i]     = hm;
        bs->cf_32[i]  = cf;

        hm[COEFFS]  = i; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        bs->red[i] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET]  =   insert_in_hash_table(e, ht);
            /* make coefficient positive */
            tmpcf             =   (int64_t)cfs[j];
            tmpcf             +=  (tmpcf >> 63) & fc;
            cf[j-off]         =   (cf32_t)tmpcf;
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_ff_32(&cf, &hm, ht);
    }
    deg_t deg = 0;
    for (i = 0; i < ngens; ++i) {
        hm  = bs->hm[i];
        deg = ht->ev[hm[OFFSET]][DEG];
        k   = hm[LENGTH] + OFFSET;
        for (j = OFFSET+1; j < k; ++j) {
            if (deg != ht->ev[hm[j]][DEG]) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
    st->homogeneous = 1;
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;
}

void import_julia_data_nf_ff_32(
        bs_t *tbr,
        ht_t *ht,
        stat_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    cf32_t *cf    = NULL;
    int64_t tmpcf = 0;
    hm_t *hm      = NULL;

    int32_t *cfs  = (int32_t *)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t fc    = st->fc;

    for (i = 0; i < start; ++i) {
        off +=  lens[i];
    }
    exp_t *e  = ht->ev[0]; /* use as temporary storage */

    int32_t nterms  = 0;
    for (i = start; i < stop; ++i) {
        nterms  +=  lens[i];
    }
    for (i = start; i < stop; ++i) {
        while (lens[i] >= ht->esz-ht->eld) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (cf32_t *)malloc((unsigned long)(lens[i]) * sizeof(cf32_t));
        tbr->hm[i-start]    = hm;
        tbr->cf_32[i-start] = cf;

        hm[COEFFS]  = i-start; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        tbr->red[i-start] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET]  =   insert_in_hash_table(e, ht);
            /* make coefficient positive */
            tmpcf             =   (int64_t)cfs[j];
            tmpcf             +=  (tmpcf >> 63) & fc;
            cf[j-off]         =   (cf32_t)tmpcf;
        }
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_ff_32(&cf, &hm, ht);
    }
}

static void sort_terms_qq(
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


void import_julia_data_nf_qq(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t start,
        const int32_t stop,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    mpz_t *cf;
    hm_t *hm;
    mpz_t prod_den, mul;
    mpz_inits(prod_den, mul, NULL);

    /* these coefficients are numerator, denominator, numerator, denominator, ...
     * i.e. the array has length 2*nterms */
    mpz_t **cfs  = (mpz_t **)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;

    /* we want to get rid of denominators, i.e. we want to handle
     * the coefficients as integers resp. mpz_t numbers. for this we
     * first get the product of all denominators of a polynomial and
     * then multiply with this product each term. the polynomials are
    * then be made content free by another function. */
    for (i = 0; i < start; ++i) {
        off +=  lens[i];
    }
    int32_t nterms  = 0;
    for (i = start; i < stop; ++i) {
        nterms  +=  lens[i];
    }

    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = start; i < stop; ++i) {
        while (lens[i] >= ht->esz) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        mpz_set_si(prod_den, 1);

        for (j = off; j < off+lens[i]; ++j) {
            mpz_mul(prod_den, prod_den, *(cfs[2*j+1]));
        }

        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (mpz_t *)malloc((unsigned long)(lens[i]) * sizeof(mpz_t));

        bs->hm[i-start]     = hm;
        bs->cf_qq[i-start]  = cf;

        for (j = 0; j < lens[i]; ++j) {
          mpz_init(cf[j]);
        }
        hm[COEFFS]  = i-start; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        bs->red[i-start] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET] = insert_in_hash_table(e, ht);
            mpz_divexact(mul, prod_den, *(cfs[2*j+1]));
            mpz_mul(cf[j-off], mul, *(cfs[2*j]));
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_qq(&cf, &hm, ht);
    }
    mpz_clears(prod_den, mul, NULL);
}


static void import_julia_data_qq(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        )
{
    int32_t i, j;
    len_t k;
    mpz_t *cf;
    hm_t *hm;
    mpz_t prod_den, mul;
    mpz_inits(prod_den, mul, NULL);

    /* these coefficients are numerator, denominator, numerator, denominator, ...
     * i.e. the array has length 2*nterms */
    mpz_t **cfs  = (mpz_t **)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;

    /* we want to get rid of denominators, i.e. we want to handle
     * the coefficients as integers resp. mpz_t numbers. for this we
     * first get the product of all denominators of a polynomial and
     * then multiply with this product each term. the polynomials are
    * then be made content free by another function. */

    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = 0; i < ngens; ++i) {
        while (lens[i] >= ht->esz) {
            enlarge_hash_table(ht);
            e  = ht->ev[0]; /* reset e if enlarging */
        }
        mpz_set_si(prod_den, 1);

        for (j = off; j < off+lens[i]; ++j) {
            /* printf("i %u | j %u\n", i, j);
             * gmp_printf("%Zd\n", *(cfs[2*j+1])); */
            mpz_mul(prod_den, prod_den, *(cfs[2*j+1]));
        }

        hm  = (hm_t *)malloc(((unsigned long)lens[i]+OFFSET) * sizeof(hm_t));
        cf  = (mpz_t *)malloc((unsigned long)(lens[i]) * sizeof(mpz_t));

        bs->hm[i]     = hm;
        bs->cf_qq[i]  = cf;

        for (j = 0; j < lens[i]; ++j) {
          mpz_init(cf[j]);
        }
        hm[COEFFS]  = i; /* link to matcf entry */
        hm[PRELOOP] = (lens[i] % UNROLL); /* offset */
        hm[LENGTH]  = lens[i]; /* length */

        bs->red[i] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            e[DEG]  = 0;
            for (k = 0; k < nv; ++k) {
                e[k+1]  = (exp_t)(exps+(nv*j))[k];
                e[DEG]  +=  e[k+1];
            }
            hm[j-off+OFFSET] = insert_in_hash_table(e, ht);
            mpz_divexact(mul, prod_den, *(cfs[2*j+1]));
            mpz_mul(cf[j-off], mul, *(cfs[2*j]));
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
        /* sort terms in polynomial w.r.t. given monomial order */
        sort_terms_qq(&cf, &hm, ht);
    }
    deg_t deg = 0;
    for (i = 0; i < ngens; ++i) {
        hm  = bs->hm[i];
        deg = ht->ev[hm[OFFSET]][DEG];
        k   = hm[LENGTH] + OFFSET;
        for (j = OFFSET+1; j < k; ++j) {
            if (deg != ht->ev[hm[j]][DEG]) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
    st->homogeneous = 1;
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;

    mpz_clears(prod_den, mul, NULL);
}

static int64_t export_julia_data_ff_8(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht,
        const uint32_t fc
        )
{
    len_t i, j, k;
    hm_t *dt;

    const len_t nv  = ht->nv;
    const len_t lml = bs->lml;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts   = 0; /* # elemnts in basis */

    for (i = 0; i < lml; ++i) {
        nterms +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)malloc(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)malloc(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    int32_t *cf   = (int32_t *)malloc(
            (unsigned long)(nterms) * sizeof(int32_t));

    /* counters for lengths, exponents and coefficients */
    int64_t cl = 0, ce = 0, cc = 0;

    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        len[cl] = bs->hm[bi][LENGTH];
        for (j = 0; j < len[cl]; ++j) {
            (cf+cc)[j] = (int32_t)bs->cf_8[bs->hm[bi][COEFFS]][j];
        }
        dt  = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 0; k < nv; ++k) {
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

static int64_t export_julia_data_ff_16(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht,
        const uint32_t fc
        )
{
    len_t i, j, k;
    hm_t *dt;

    const len_t nv  = ht->nv;
    const len_t lml = bs->lml;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts   = 0; /* # elemnts in basis */

    for (i = 0; i < lml; ++i) {
        nterms +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)malloc(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)malloc(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    int32_t *cf   = (int32_t *)malloc(
            (unsigned long)(nterms) * sizeof(int32_t));

    /* counters for lengths, exponents and coefficients */
    int64_t cl = 0, ce = 0, cc = 0;

    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        len[cl] = bs->hm[bi][LENGTH];
        for (j = 0; j < len[cl]; ++j) {
            (cf+cc)[j] = (int32_t)bs->cf_16[bs->hm[bi][COEFFS]][j];
        }
        dt  = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 0; k < nv; ++k) {
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

static int64_t export_julia_data_ff_32(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht,
        const uint32_t fc
        )
{
    len_t i, j, k;
    hm_t *dt;

    const len_t nv  = ht->nv;
    const len_t lml = bs->lml;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts   = 0; /* # elemnts in basis */

    for (i = 0; i < lml; ++i) {
        nterms +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)malloc(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)malloc(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    int32_t *cf   = (int32_t *)malloc(
            (unsigned long)(nterms) * sizeof(int32_t));

    /* counters for lengths, exponents and coefficients */
    int64_t cl = 0, ce = 0, cc = 0, ctmp  = 0;;

    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        len[cl] = bs->hm[bi][LENGTH];
        /* we may need to make coefficients negative since we
         * output int32_t but cf32_t is uint_32. */
        for (j = 0; j < len[cl]; ++j) {
            ctmp  = (int64_t)bs->cf_32[bs->hm[bi][COEFFS]][j];
            ctmp  -=  (ctmp >> 31) & fc;
            (cf+cc)[j] = (int32_t)ctmp;
        }
        memcpy(cf+cc, bs->cf_32[bs->hm[bi][COEFFS]],
                (unsigned long)(len[cl]) * sizeof(cf32_t));

        dt  = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 0; k < nv; ++k) {
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

static int64_t export_julia_data_qq(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht,
        const uint32_t fc
        )
{
    len_t i, j, k;
    hm_t *dt;

    const len_t nv  = ht->nv;
    const len_t lml = bs->lml;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts   = 0; /* # elemnts in basis */

    for (i = 0; i < lml; ++i) {
        nterms +=  (int64_t)bs->hm[bs->lmps[i]][LENGTH];
    }
    nelts = lml;

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)malloc(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)malloc(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    mpz_t *cf     = (mpz_t *)malloc(
            (unsigned long)(nterms) * sizeof(mpz_t));

    /* counters for lengths, exponents and coefficients */
    int64_t cl = 0, ce = 0, cc = 0;
    for (i = 0; i < lml; ++i) {
        const bl_t bi = bs->lmps[i];
        len[cl] = bs->hm[bi][LENGTH];
        mpz_t *coeffs =  bs->cf_qq[bs->hm[bi][COEFFS]];
        for (j = 0; j < len[cl]; ++j) {
            mpz_init_set((cf+cc)[j], coeffs[j]);
        }

        dt  = bs->hm[bi] + OFFSET;
        for (j = 0; j < len[cl]; ++j) {
            for (k = 0; k < nv; ++k) {
                exp[ce++] = (int32_t)ht->ev[dt[j]][k];
            }
        }
        cc  += len[cl];
        cl++;
    }

    *bload  = (int32_t)nelts;
    *blen   = len;
    *bexp   = exp;
    *bcf    = (void *)cf;

    return nterms;
}

int32_t check_and_set_meta_data(
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    if (nr_gens <= 0
            || nr_vars <= 0
            || field_char < 0
            || lens == NULL
            || cfs == NULL
            || exps == NULL) {
        return 1;
    }

    st->ngens = nr_gens;
    st->nvars = nr_vars;
    /* note: prime check should be done in julia */
    st->fc    = field_char;
    if (st->fc == 0) {
        st->ff_bits = 0;
    } else {
        if (st->fc < pow(2,8)) {
            st->ff_bits = 8;
        } else {
            if (st->fc < pow(2,16)) {
                st->ff_bits = 16;
            } else {
                if (st->fc < pow(2,32)) {
                    st->ff_bits = 32;
                }
            }
        }
    }
    /* monomial order */
    if (mon_order != 0 && mon_order != 1) {
        st->mo  = 0;
    } else {
        st->mo  = mon_order;
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
        const stat_t *st
        )
{
  /* todo: this needs to be generalized for different monomial orders */
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
      initialize_basis        = initialize_basis_qq;
      import_julia_data       = import_julia_data_qq;
      export_julia_data       = export_julia_data_qq;
      check_enlarge_basis     = check_enlarge_basis_qq;
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
      initialize_basis            = initialize_basis_ff_8;
      import_julia_data           = import_julia_data_ff_8;
      export_julia_data           = export_julia_data_ff_8;
      check_enlarge_basis         = check_enlarge_basis_ff_8;
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
      initialize_basis            = initialize_basis_ff_16;
      import_julia_data           = import_julia_data_ff_16;
      export_julia_data           = export_julia_data_ff_16;
      check_enlarge_basis         = check_enlarge_basis_ff_16;
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
      initialize_basis            = initialize_basis_ff_32;
      import_julia_data           = import_julia_data_ff_32;
      export_julia_data           = export_julia_data_ff_32;
      check_enlarge_basis         = check_enlarge_basis_ff_32;
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
      initialize_basis            = initialize_basis_ff_32;
      import_julia_data           = import_julia_data_ff_32;
      export_julia_data           = export_julia_data_ff_32;
      check_enlarge_basis         = check_enlarge_basis_ff_32;
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
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const uint32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
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
    return check_and_set_meta_data(st, lens, exps,
            cfs, field_char, mon_order, nr_vars, nr_gens,
            ht_size, nr_threads, max_nr_pairs, reset_hash_table,
            la_option, reduce_gb, pbm_file, info_level);
}

static inline void reset_function_pointers(
        const uint32_t prime,
        const uint32_t laopt
        )
{
    if (prime < pow(2,8)) {
        copy_basis_mod_p            = copy_basis_mod_p_8;
        interreduce_matrix_rows     = interreduce_matrix_rows_ff_8;
        initialize_basis            = initialize_basis_ff_8;
        import_julia_data           = import_julia_data_ff_8;
        export_julia_data           = export_julia_data_ff_8;
        check_enlarge_basis         = check_enlarge_basis_ff_8;
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
            copy_basis_mod_p            = copy_basis_mod_p_16;
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_16;
            initialize_basis            = initialize_basis_ff_16;
            import_julia_data           = import_julia_data_ff_16;
            export_julia_data           = export_julia_data_ff_16;
            check_enlarge_basis         = check_enlarge_basis_ff_16;
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
            copy_basis_mod_p            = copy_basis_mod_p_32;
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
            initialize_basis            = initialize_basis_ff_32;
            import_julia_data           = import_julia_data_ff_32;
            export_julia_data           = export_julia_data_ff_32;
            check_enlarge_basis         = check_enlarge_basis_ff_32;
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
        copy_basis_mod_p            = copy_basis_mod_p_8;
        interreduce_matrix_rows     = interreduce_matrix_rows_ff_8;
        initialize_basis            = initialize_basis_ff_8;
        import_julia_data           = import_julia_data_ff_8;
        export_julia_data           = export_julia_data_ff_8;
        check_enlarge_basis         = check_enlarge_basis_ff_8;
        normalize_initial_basis     = normalize_initial_basis_ff_8;
        application_linear_algebra  = exact_application_sparse_linear_algebra_ff_8;
        trace_linear_algebra        = exact_trace_sparse_linear_algebra_ff_8;
    } else {
        if (prime < pow(2,16)) {
            copy_basis_mod_p            = copy_basis_mod_p_16;
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_16;
            initialize_basis            = initialize_basis_ff_16;
            import_julia_data           = import_julia_data_ff_16;
            export_julia_data           = export_julia_data_ff_16;
            check_enlarge_basis         = check_enlarge_basis_ff_16;
            normalize_initial_basis     = normalize_initial_basis_ff_16;
            application_linear_algebra  = exact_application_sparse_linear_algebra_ff_16;
            trace_linear_algebra        = exact_trace_sparse_linear_algebra_ff_16;
        } else {
            copy_basis_mod_p            = copy_basis_mod_p_32;
            interreduce_matrix_rows     = interreduce_matrix_rows_ff_32;
            initialize_basis            = initialize_basis_ff_32;
            import_julia_data           = import_julia_data_ff_32;
            export_julia_data           = export_julia_data_ff_32;
            check_enlarge_basis         = check_enlarge_basis_ff_32;
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
    const stat_t * const st
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

    sprintf(fn, "%d-%d-%d-%d.pbm", st->current_rd, nru+nrl, ncols, st->current_deg);
    FILE *fh  = fopen(fn, "wb");

    /* magic header */
    sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", nru+nrl, ncols, ncols, nru+nrl);

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
