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


#include "hash.h"

/* we have three different hash tables:
 * 1. one hash table for elements in the basis (bht)
 * 2. one hash table for the spairs during the update process (uht)
 * 3. one hash table for the multiplied elements during symbolic
 *    preprocessing (sht) */

/* The idea of the structure of the hash table is taken from an
 * implementation by Roman Pearce and Michael Monagan in Maple. */

static val_t pseudo_random_number_generator(
    uint32_t *seed
    )
{
    uint32_t rseed  = *seed;
    rseed ^=  (rseed << 13);
    rseed ^=  (rseed >> 17);
    rseed ^=  (rseed << 5);
    *seed =   rseed;
    return (val_t)rseed;
}

ht_t *initialize_basis_hash_table(
    md_t *st
    )
{
    len_t i;
    hl_t j;

    const len_t nv  = st->nvars;

    ht_t *ht  = (ht_t *)malloc(sizeof(ht_t));
    ht->nv    = nv;
    /* generate map */
    ht->bpv = (len_t)((CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nv);
    if (ht->bpv == 0) {
        ht->bpv++;
    }
    ht->ndv = (unsigned long)nv < (CHAR_BIT * sizeof(sdm_t)) ?
        nv : (len_t)((CHAR_BIT * sizeof(sdm_t)));
    ht->dv  = (len_t *)calloc((unsigned long)ht->ndv, sizeof(len_t));

    ht->hsz   = (hl_t)pow(2, st->init_hts);
    ht->esz   = ht->hsz / 2;
    ht->hmap  = calloc(ht->hsz, sizeof(hi_t));

    if (st->nev == 0) {
        ht->evl = nv + 1; /* store also degree at first position */
        ht->ebl = 0;
        for (i = 1; i <= ht->ndv; ++i) {
            ht->dv[i-1] = i;
        }
    } else {
        ht->evl = nv + 2; /* store also degrees for both blocks, see
                           * data.h for more on exponent vector structure */
        ht->ebl = st->nev + 1; /* store also degree at first position */
        if (st->nev >= ht->ndv) {
            for (i = 1; i <= ht->ndv; ++i) {
                ht->dv[i-1] = i;
            }
        } else {
            len_t ctr = 0;
            for (i = 1; i <= st->nev; ++i) {
                ht->dv[ctr++] = i;
            }
            for (i = ht->ebl+1; i < ht->ndv+2; ++i) {
                ht->dv[ctr++] = i;
            }
        }

    }
    /* generate divmask map */
    ht->dm  = (sdm_t *)calloc(
            (unsigned long)(ht->ndv * ht->bpv), sizeof(sdm_t));

    /* generate random values */
    ht->rsd = 2463534242;
    ht->rn  = calloc((unsigned long)ht->evl, sizeof(val_t));
    for (i = ht->evl; i > 0; --i) {
        /* random values should not be zero */
        ht->rn[i-1] = pseudo_random_number_generator(&(ht->rsd)) | 1;
    }
    /* generate exponent vector */
    /* keep first entry empty for faster divisibility checks */
    ht->eld = 1;
    ht->hd  = (hd_t *)calloc(ht->esz, sizeof(hd_t));
    ht->ev  = (exp_t **)malloc(ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        fprintf(stderr, "Computation needs too much memory on this machine,\n");
        fprintf(stderr, "could not initialize exponent vector for hash table,\n");
        fprintf(stderr, "esz = %lu, segmentation fault will follow.\n", (unsigned long)ht->esz);
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)ht->evl * ht->esz * sizeof(exp_t));
    if (tmp == NULL) {
        fprintf(stderr, "Exponent storage needs too much memory on this machine,\n");
        fprintf(stderr, "initialization failed, esz = %lu,\n", (unsigned long)ht->esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    const hl_t esz  = ht->esz;
    for (j = 0; j < esz; ++j) {
        ht->ev[j]  = tmp + (j*ht->evl);
    }
    st->max_bht_size  = ht->esz;
    return ht;
}

ht_t *copy_hash_table(
    const ht_t *bht
    )
{
    hl_t j;

    ht_t *ht  = (ht_t *)malloc(sizeof(ht_t));

    ht->nv    = bht->nv;
    ht->evl   = bht->evl;
    ht->ebl   = bht->ebl;
    ht->hsz   = bht->hsz;
    ht->esz   = bht->esz;

    ht->hmap  = calloc(ht->hsz, sizeof(hi_t));
    memcpy(ht->hmap, bht->hmap, (unsigned long)ht->hsz * sizeof(hi_t));

    ht->ndv = bht->ndv;
    ht->bpv = bht->bpv;
    ht->dm  = bht->dm;
    ht->rn  = bht->rn;

    ht->dv  = (len_t *)calloc((unsigned long)ht->ndv, sizeof(len_t));
    memcpy(ht->dv, bht->dv, (unsigned long)ht->ndv * sizeof(len_t));

    /* generate exponent vector */
    /* keep first entry empty for faster divisibility checks */
    ht->hd  = (hd_t *)calloc(ht->esz, sizeof(hd_t));

    memcpy(ht->hd, bht->hd, (unsigned long)ht->esz * sizeof(hd_t));
    ht->ev  = (exp_t **)malloc(ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        fprintf(stderr, "Computation needs too much memory on this machine,\n");
        fprintf(stderr, "could not initialize exponent vector for hash table,\n");
        fprintf(stderr, "esz = %lu, segmentation fault will follow.\n", (unsigned long)ht->esz);
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)ht->evl * ht->esz * sizeof(exp_t));
    if (tmp == NULL) {
        fprintf(stderr, "Exponent storage needs too much memory on this machine,\n");
        fprintf(stderr, "initialization failed, esz = %lu,\n", (unsigned long)ht->esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    memcpy(tmp, bht->ev[0], (unsigned long)ht->evl * ht->esz * sizeof(exp_t));
    ht->eld = bht->eld;
    const hl_t esz  = ht->esz;
    for (j = 0; j < esz; ++j) {
        ht->ev[j]  = tmp + (j*ht->evl);
    }
    return ht;
}


ht_t *initialize_secondary_hash_table(
    const ht_t * const bht,
    const md_t * const md
    )
{
    hl_t j;

    ht_t *ht  = (ht_t *)malloc(sizeof(ht_t)); 
    ht->nv    = bht->nv;
    ht->evl   = bht->evl;
    ht->ebl   = bht->ebl;

    /* generate map */
    int32_t min = 3 > md->init_hts-5 ? 3 : md->init_hts-5;
    ht->hsz   = (hl_t)pow(2, min);
    ht->esz   = ht->hsz / 2;
    ht->hmap  = calloc(ht->hsz, sizeof(hi_t));

    /* divisor mask and random number seeds from basis hash table */
    ht->ndv = bht->ndv;
    ht->bpv = bht->bpv;
    ht->dm  = bht->dm;
    ht->rn  = bht->rn;
    ht->dv  = bht->dv;

    /* generate exponent vector */
    /* keep first entry empty for faster divisibility checks */
    ht->eld = 1;
    ht->hd  = (hd_t *)calloc(ht->esz, sizeof(hd_t));
    ht->ev  = (exp_t **)malloc(ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        fprintf(stderr, "Computation needs too much memory on this machine,\n");
        fprintf(stderr, "could not initialize exponent vector for hash table,\n");
        fprintf(stderr, "esz = %lu, segmentation fault will follow.\n", (unsigned long)ht->esz);
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)ht->evl * ht->esz * sizeof(exp_t));
    if (tmp == NULL) {
        fprintf(stderr, "Exponent storage needs too much memory on this machine,\n");
        fprintf(stderr, "initialization failed, esz = %lu,\n", (unsigned long)ht->esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    const hl_t esz  = ht->esz;
    for (j = 0; j < esz; ++j) {
        ht->ev[j]  = tmp + (j*ht->evl);
    }
    return ht;
}

void free_shared_hash_data(
    ht_t *ht
    )
{
    if (ht != NULL) {
        if (ht->rn) {
            free(ht->rn);
            ht->rn = NULL;
        }
        if (ht->dv) {
            free(ht->dv);
            ht->dv = NULL;
        }
        if (ht->dm) {
            free(ht->dm);
            ht->dm = NULL;
        }
    }
}

void free_hash_table(
    ht_t **htp
    )
{
    ht_t *ht  = *htp;
    if (ht->hmap) {
        free(ht->hmap);
        ht->hmap = NULL;
    }
    if (ht->hd) {
        free(ht->hd);
        ht->hd  = NULL;
    }
    if (ht->ev) {
        /* note: memory is allocated as one big block,
            *       so freeing ev[0] is enough */
        free(ht->ev[0]);
        free(ht->ev);
        ht->ev  = NULL;
    }
    free(ht);
    ht    = NULL;
    *htp  = ht;
}

void full_free_hash_table(
                     ht_t **htp
                     )
{
  ht_t *ht  = *htp;
  if (ht->hmap) {
    free(ht->hmap);
    ht->hmap = NULL;
  }
  if (ht->hd) {
    free(ht->hd);
    ht->hd  = NULL;
  }
  if (ht->ev) {
    /* note: memory is allocated as one big block,
     *       so freeing ev[0] is enough */
    free(ht->ev[0]);
    free(ht->ev);
    ht->ev  = NULL;
  }
  if (ht != NULL) {
    if (ht->rn) {
      free(ht->rn);
      ht->rn = NULL;
    }
    if (ht->dv) {
      free(ht->dv);
      ht->dv = NULL;
    }
    if (ht->dm) {
      free(ht->dm);
      ht->dm = NULL;
    }
  }
  free(ht);
  ht    = NULL;
  *htp  = ht;
}

/* we just double the hash table size */
static void enlarge_hash_table(
    ht_t *ht
    )
{
    hl_t i, j;
    val_t h, k;

    ht->esz = 2 * ht->esz;
    const hl_t esz  = ht->esz;
    const hi_t eld  = ht->eld;

    ht->hd    = realloc(ht->hd, esz * sizeof(hd_t));
    memset(ht->hd+eld, 0, (esz-eld) * sizeof(hd_t));
    ht->ev    = realloc(ht->ev, esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        fprintf(stderr, "Enlarging hash table failed for esz = %lu,\n", (unsigned long)esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    /* note: memory is allocated as one big block, so reallocating
     *       memory from ev[0] is enough    */
    ht->ev[0] = realloc(ht->ev[0],
            esz * (unsigned long)ht->evl * sizeof(exp_t));
    if (ht->ev[0] == NULL) {
        fprintf(stderr, "Enlarging exponent vector for hash table failed\n");
        fprintf(stderr, "for esz = %lu, segmentation fault will follow.\n", (unsigned long)esz);
    }
    /* due to realloc we have to reset ALL ev entries,
     * memory might have been moved */
    for (i = 1; i < esz; ++i) {
        ht->ev[i] = ht->ev[0] + (i*ht->evl);
    }

    /* The hash table should be double the size of the exponent space in
     * order to never get a fill in over 50%. If the exponent size is now
     * enlarge to 2^31 elements that's the limit we can go. Thus we cannot
     * enlarge the hash table size any further and have to live with more
     * than 50% fill in. */
    if (ht->hsz < (hl_t)pow(2,32)) {
        ht->hsz = 2 * ht->hsz;
        const hl_t hsz  = ht->hsz;
        ht->hmap  = realloc(ht->hmap, hsz * sizeof(hi_t));
        if (ht->hmap == NULL) {
            fprintf(stderr, "Enlarging hash table failed for hsz = %lu,\n", (unsigned long)hsz);
            fprintf(stderr, "segmentation fault will follow.\n");
        }
        memset(ht->hmap, 0, hsz * sizeof(hi_t));
        const hi_t mod =  (hi_t )(hsz-1);

        /* reinsert known elements */
        for (i = 1; i < eld; ++i) {
            h = ht->hd[i].val;

            /* probing */
            k = h;
            for (j = 0; j < hsz; ++j) {
                k = (k+j) & mod;
                if (ht->hmap[k]) {
                    continue;
                }
                ht->hmap[k] = i;
                break;
            }
        }
    } else {
        if (ht->hsz == (hl_t)pow(2,32)) {
          printf("Exponent space is now 2^32 elements wide, we cannot\n");
          printf("enlarge the hash table any further, thus fill in gets\n");
          printf("over 50%% and performance of hashing may get worse.\n");
        } else {
          printf("Hash table is full, we can no longer enlarge\n");
          printf("Segmentation fault will follow.\n");
          free(ht->hmap);
          ht->hmap  = NULL;
        }
    }
}

static inline sdm_t generate_short_divmask(
    const exp_t * const a,
    const ht_t *ht
    )
{
  len_t i, j;
  int32_t res = 0;
  int32_t ctr = 0;
  const len_t ndv         = ht->ndv;
  const len_t * const dv  = ht->dv;
  const len_t bpv         = ht->bpv;

  for (i = 0; i < ndv; ++i) {
    for (j = 0; j < bpv; ++j) {
      if ((sdm_t)a[dv[i]] >= ht->dm[ctr]) {
        res |= 1 << ctr;
      }
      ctr++;
    }
  }
 
  return res;
}

/* note: we calculate the divmask after reading in the input generators. thoseV
 * are first stored in the local hash table. thus we use the local exponents to
 * generate the divmask */
void calculate_divmask(
    ht_t *ht
    )
{
  hi_t i;
  hl_t k;
  len_t j, steps;
  int32_t ctr = 0;
  const len_t * const dv  = ht->dv;
  exp_t **ev  = ht->ev;

  deg_t *max_exp  = (deg_t *)malloc((unsigned long)ht->ndv * sizeof(deg_t));
  deg_t *min_exp  = (deg_t *)malloc((unsigned long)ht->ndv * sizeof(deg_t));

  exp_t *e  = ev[1];

  /* get initial values from first hash table entry */
  for (i = 0; i < ht->ndv; ++i) {
    max_exp[i]  = min_exp[i]  = e[dv[i]];
  }

  /* get maximal and minimal exponent element entries in hash table */
  for (i = 2; i < ht->eld; ++i) {
    e = ev[i];
    for (j = 0; j < ht->ndv; ++j) {
      if (e[dv[j]] > max_exp[j]) {
        max_exp[j]  = e[dv[j]];
        continue;
      }
      if (e[dv[j]] < min_exp[j]) {
        min_exp[j]  = e[dv[j]];
      }
    }
  }

  /* calculate average values for generating divmasks */
  for (i = 0; i < ht->ndv; ++i) {
    steps = (max_exp[i] - min_exp[i]) / ht->bpv;
    if (steps == 0)
      steps++;
    for (j = 0; j < ht->bpv; ++j) {
      ht->dm[ctr++] = (sdm_t)steps++;
    }
  }

  /* initialize divmasks for elements already added to hash table */
  for (k = 1; k < ht->eld; k++) {
    ht->hd[k].sdm = generate_short_divmask(ev[k], ht);
  }

  free(max_exp);
  free(min_exp);
}

/* returns zero if a is not divisible by b, else 1 is returned */
static inline hi_t check_monomial_division(
    const hi_t a,
    const hi_t b,
    const ht_t *ht
    )
{
  len_t i;

  /* short divisor mask check */
  if (ht->hd[b].sdm & ~ht->hd[a].sdm) {
    return 0;
  }

  const len_t evl = ht->evl;

  const exp_t *const ea = ht->ev[a];
  const exp_t *const eb = ht->ev[b];

  /* printf("! no sdm decision !\n"); */
  /* exponent check */
  for (i = 0; i < evl-1; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    }
  }
  if (ea[evl-1] < eb[evl-1]) {
    return 0;
  }
  return 1;
}

static inline void check_monomial_division_in_update(
    hi_t *a,
    const len_t start,
    const len_t end,
    const hi_t b,
    const ht_t *ht
    )
{
    len_t i, j;
    const len_t evl = ht->evl;

    const sdm_t sb        = ht->hd[b].sdm;
    const exp_t *const eb = ht->ev[b];
    /* pairs are sorted, we only have to search entries
     * above the starting point */
        j = start+1;
restart:
    for (; j < end; ++j) {
        if (a[j] == 0) {
            continue;
        }
        /* short divisor mask check */
        if (~ht->hd[a[j]].sdm & sb) {
            continue;
        }
        const exp_t *const ea = ht->ev[a[j]];
        /* exponent check */
        for (i = 0; i < evl-1; i += 2) {
            if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
                j++;
                goto restart;
            }
        }
        if (ea[evl-1] < eb[evl-1]) {
            continue;
        }
        a[j]  = 0;
    }
}

static inline hi_t check_lm_divisibility_and_insert_in_hash_table(
    const exp_t *a,
    ht_t *ht,
    const bs_t * const bs
    )
{
    hl_t i;
    hi_t k, pos;
    len_t j;
    exp_t *e;
    hd_t *d;
    const len_t lml   = bs->lml;

    const sdm_t * const lms = bs->lm;
    const bl_t * const lmps = bs->lmps;

    const sdm_t nsdm  = ~generate_short_divmask(a, ht);

    val_t h = 0;
    const len_t evl = ht->evl;
    const hl_t hsz  = ht->hsz;
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod = (hi_t)(ht->hsz - 1);

    /* check divisibility w.r.t. current lead monomials */
    i = 0;
start:
    while (i < lml && lms[i] & nsdm) {
        i++;
    }
    if (i < lml) {
        e = ht->ev[bs->hm[lmps[i]][OFFSET]];
        for (j = 0; j < evl; ++j) {
            if (e[j] > a[j]) {
                i++;
                goto start;
            }
        }
        /* divisible by lm */
        return 0;
    }
    /* if we are here then a is not divisible by a current
     * lead monomial and we can add it to the hash table */

    /* generate hash value */
    for (j = 0; j < evl; ++j) {
        h +=  ht->rn[j] * a[j];
    }
    /* probing */
    k = h;
    i = 0;
restart:
    for (; i < hsz; ++i) {
        k = (hi_t)((k+i) & mod);
        const hi_t hm = ht->hmap[k];
        if (!hm) {
            break;
        }
        if (ht->hd[hm].val != h) {
            continue;
        }
        const exp_t * const ehm = ht->ev[hm];
        for (j = 0; j < evl-1; j += 2) {
            if (a[j] != ehm[j] || a[j+1] != ehm[j+1]) {
                i++;
                goto restart;
            }
        }
        if (a[evl-1] != ehm[evl-1]) {
            i++;
            goto restart;
        }
        return hm;
    }

    /* add element to hash table */
    ht->hmap[k]  = pos = (hi_t)ht->eld;
    e   = ht->ev[pos];
    d   = ht->hd + pos;
    memcpy(e, a, (unsigned long)evl * sizeof(exp_t));
    d->sdm  =   generate_short_divmask(e, ht);
    d->deg  =   e[0];
    d->deg  +=  ht->ebl > 0 ? e[ht->ebl] : 0;
    d->val  =   h;

    ht->eld++;

    return pos;
}

static inline hi_t insert_multiplied_signature_in_hash_table(
    const hm_t h1,
    const hm_t h2,
    ht_t *ht
    )
{
    hl_t i;
    hi_t k, pos;
    len_t j;
    exp_t *e;
    exp_t *a = ht->ev[0];
    hd_t *d;
    val_t h = 0;
    const len_t evl = ht->evl;
    const hl_t hsz = ht->hsz;
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod = (hi_t)(ht->hsz - 1);

    h   =   h1 + h2;

    /* generate exponent vector */
    for (j = 0; j < evl; ++j) {
        a[j] = ht->ev[h1][j] + ht->ev[h2][j];
    }
    /* probing */
    k = h;
    i = 0;
restart:
    for (; i < hsz; ++i) {
        k = (hi_t)((k+i) & mod);
        const hi_t hm = ht->hmap[k];
        if (!hm) {
            break;
        }
        if (ht->hd[hm].val != h) {
            continue;
        }
        const exp_t * const ehm = ht->ev[hm];
        for (j = 0; j < evl-1; j += 2) {
            if (a[j] != ehm[j] || a[j+1] != ehm[j+1]) {
                i++;
                goto restart;
            }
        }
        if (a[evl-1] != ehm[evl-1]) {
            i++;
            goto restart;
        }
        return hm;
    }

    /* add element to hash table */
    ht->hmap[k]  = pos = (hi_t)ht->eld;
    e   = ht->ev[pos];
    d   = ht->hd + pos;
    memcpy(e, a, (unsigned long)evl * sizeof(exp_t));
    d->sdm  =   generate_short_divmask(e, ht);
    d->deg  =   e[0];
    d->deg  +=  ht->ebl > 0 ? e[ht->ebl] : 0;
    d->val  =   h;

    ht->eld++;

    return pos;
}

/* If the exponent vector is not contained in the hash table
 * we return 0 and hp is a pointer to the hash value and kp is
 * a pointer to the index of hmap where to store the exponent.
 * If the exponent vector is contained in the hash table we
 * return 1 and kp is the pointer of the index where the
 * exponent vector is stored in the exponents array of the
 * hash table. */
static inline int32_t is_contained_in_hash_table(
        const exp_t *a,
        const ht_t * const ht,
        const val_t h,
        hi_t *kp
        )
{
    hl_t i;
    hi_t k;
    len_t j;
    /* const len_t evl = ht->evl;
     * const hl_t hsz = ht->hsz; */
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod = (hi_t)(ht->hsz - 1);

    /* probing */
    k = h;
    i = 0;
restart:
    for (; i < ht->hsz; ++i) {
        k = (hi_t)((k+i) & mod);
        const hi_t hm = ht->hmap[k];
        if (!hm) {
            *kp = k;
            return 0;
        }
        if (ht->hd[hm].val != h) {
            continue;
        }
        const exp_t * const ehm = ht->ev[hm];
        for (j = 0; j < ht->evl-1; j += 2) {
            if (a[j] != ehm[j] || a[j+1] != ehm[j+1]) {
                i++;
                goto restart;
            }
        }
        if (a[ht->evl-1] != ehm[ht->evl-1]) {
            i++;
            goto restart;
        }
        *kp = hm;
        return 1;
    }
    return -1;
}

/* This function assumes that is_contained_in_hash_table() was
 * called beforehand such that the values for h and k are already
 * precomputed. */
static inline len_t add_to_hash_table(
    const exp_t * const a,
    const val_t h,
    const hi_t k,
    ht_t *ht
    )
{
    /* add element to hash table */
    hi_t pos;
    ht->hmap[k] = pos = (hi_t)ht->eld;
    exp_t *e    = ht->ev[pos];
    hd_t *d     = ht->hd + pos;
    memcpy(e, a, (unsigned long)ht->evl * sizeof(exp_t));
    d->sdm  =   generate_short_divmask(e, ht);
    d->deg  =   e[0];
    d->deg  +=  ht->ebl > 0 ? e[ht->ebl] : 0;
    d->val  =   h;

    ht->eld++;

    return pos;
}

static inline len_t check_insert_in_hash_table(
        const exp_t *a,
        val_t h,
        ht_t *ht
        )
{
    if (h == 0) {
        /* generate hash value */
        for (len_t j = 0; j < ht->evl; ++j) {
            h +=  ht->rn[j] * a[j];
        }
    }

    hi_t k  = 0;

#if 1
    len_t ld = 0;
    while (1) {
        ld = ht->eld;
        if (is_contained_in_hash_table(a, ht, h, &k)) {
            return k;
        } else {
            if (ht->eld == ld) {
#pragma omp critical
                ld = add_to_hash_table(a, h, k, ht);
                return ld;
            }
        }
    }
#else
    return is_contained_in_hash_table(a, ht, h, &k) ?
        k : add_to_hash_table(a, h, k, ht);
#endif
}

static inline hi_t insert_in_hash_table(
    const exp_t *a,
    ht_t *ht
    )
{
    hl_t i;
    hi_t k, pos;
    len_t j;
    exp_t *e;
    hd_t *d;
    val_t h = 0;
    const len_t evl = ht->evl;
    const hl_t hsz = ht->hsz;
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod = (hi_t)(ht->hsz - 1);

    /* generate hash value */
    for (j = 0; j < evl; ++j) {
        h +=  ht->rn[j] * a[j];
    }
    /* probing */
    k = h;
    i = 0;
restart:
    for (; i < hsz; ++i) {
        k = (hi_t)((k+i) & mod);
        const hi_t hm = ht->hmap[k];
        if (!hm) {
            break;
        }
        if (ht->hd[hm].val != h) {
            continue;
        }
        const exp_t * const ehm = ht->ev[hm];
        for (j = 0; j < evl-1; j += 2) {
            if (a[j] != ehm[j] || a[j+1] != ehm[j+1]) {
                i++;
                goto restart;
            }
        }
        if (a[evl-1] != ehm[evl-1]) {
            i++;
            goto restart;
        }
        return hm;
    }

    /* add element to hash table */
    ht->hmap[k]  = pos = (hi_t)ht->eld;
    e   = ht->ev[pos];
    d   = ht->hd + pos;
    memcpy(e, a, (unsigned long)evl * sizeof(exp_t));
    d->sdm  =   generate_short_divmask(e, ht);
    d->deg  =   e[0];
    d->deg  +=  ht->ebl > 0 ? e[ht->ebl] : 0;
    d->val  =   h;

    ht->eld++;

    return pos;
}

static inline void reinitialize_hash_table(
    ht_t *ht,
    const hl_t size
    )
{
    hl_t i;
    /* is there still enough space in the local table? */
    if (size >= (ht->esz)) {
        while (size >= ht->esz) {
            ht->esz = 2 * ht->esz;
            ht->hsz = 2 * ht->hsz;
        }
        const hl_t esz  = ht->esz;
        const hl_t hsz  = ht->hsz;
        const len_t evl = ht->evl;
        ht->hd  = realloc(ht->hd, esz * sizeof(hd_t));
        ht->ev  = realloc(ht->ev, esz * sizeof(exp_t *));
        if (ht->ev == NULL) {
            fprintf(stderr, "Computation needs too much memory on this machine,\n");
            fprintf(stderr, "could not reinitialize exponent vector for hash table,\n");
            fprintf(stderr, "esz = %lu, segmentation fault will follow.\n", (unsigned long)esz);
        }
        /* note: memory is allocated as one big block, so reallocating
         *       memory from evl[0] is enough    */
        ht->ev[0]  = realloc(ht->ev[0],
                esz * (unsigned long)evl * sizeof(exp_t));
        if (ht->ev[0] == NULL) {
            fprintf(stderr, "Exponent storage needs too much memory on this machine,\n");
            fprintf(stderr, "reinitialization failed, esz = %lu\n", (unsigned long)esz);
            fprintf(stderr, "segmentation fault will follow.\n");
        }
        /* due to realloc we have to reset ALL evl entries, memory might be moved */
        for (i = 1; i < esz; ++i) {
            ht->ev[i] = ht->ev[0] + (i*evl);
        }
        ht->hmap  = realloc(ht->hmap, hsz * sizeof(hi_t));
    }
    memset(ht->hd, 0, ht->esz * sizeof(hd_t));
    memset(ht->hmap, 0, ht->hsz * sizeof(hi_t));

    ht->eld  = 1;
}

static inline void clean_hash_table(
        ht_t *ht
    )
{
    memset(ht->hd, 0, ht->esz * sizeof(hd_t));
    memset(ht->hmap, 0, ht->hsz * sizeof(hi_t));

    ht->eld  = 1;
}

static inline int prime_monomials(
    const hi_t a,
    const hi_t b,
    const ht_t *ht
    )
{
    len_t i;

    const exp_t * const ea = ht->ev[a];
    const exp_t * const eb = ht->ev[b];

    const len_t evl = ht->evl;
    const len_t ebl = ht->ebl;

    for (i = 1; i < ebl; ++i) {
        if (ea[i] != 0 && eb[i] != 0) {
            return 0;
        }
    }
    for (i = ebl+1; i < evl; ++i) {
        if (ea[i] != 0 && eb[i] != 0) {
            return 0;
        }
    }
    return 1;
}

static inline void insert_plcms_in_basis_hash_table(
    ps_t *psl,
    spair_t *pp,
    ht_t *bht,
    const ht_t *uht,
    const bs_t * const bs,
    const hi_t * const lcms,
    const len_t start,
    const len_t end
    )
{
    hl_t i;
    hi_t k, pos;
    len_t j, l, m;
    hd_t *d;

    spair_t *ps     = psl->p;
    const len_t evl = bht->evl;
    const hl_t hsz  = bht->hsz;
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod = (hi_t)(hsz - 1);
    hm_t * const * const hm = bs->hm;
    m = start;
    l = 0;
letsgo:
    for (; l < end; ++l) {
        if (lcms[l] == 0) {
            continue;
        }
        if (prime_monomials(
                    hm[pp[l].gen1][OFFSET], hm[pp[0].gen2][OFFSET], bht)) {
            continue;
        }
        ps[m] = pp[l];
        const val_t h = uht->hd[lcms[l]].val;
        memcpy(bht->ev[bht->eld], uht->ev[lcms[l]],
                (unsigned long)evl * sizeof(exp_t));
        const exp_t * const n = bht->ev[bht->eld];
        k = h;
        i = 0;
restart:
        for (; i < hsz; ++i) {
            k = (hi_t)(k+i) & mod;
            const hi_t hm = bht->hmap[k];
            if (!hm) {
                break;
            }
            if (bht->hd[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = bht->ev[hm];
            for (j = 0; j < evl-1; j += 2) {
                if (n[j] != ehm[j] || n[j+1] != ehm[j+1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[evl-1] != ehm[evl-1]) {
                i++;
                goto restart;
            }
            ps[m++].lcm = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        bht->hmap[k] = pos = (hi_t)bht->eld;
        d = bht->hd + bht->eld;
        d->sdm  = uht->hd[lcms[l]].sdm;
        d->deg  = uht->hd[lcms[l]].deg;
        d->val  = h;

        bht->eld++;
        ps[m++].lcm =  pos;
    }
    psl->ld = m;
}

static inline void switch_hcm_data_to_basis_hash_table(
    hi_t *hcm,
    ht_t *bht,
    const mat_t *mat,
    const ht_t * const sht
    )
{
    const len_t start = mat->ncl;
    const len_t end   = mat->nc;

    while (bht->esz - bht->eld < mat->ncr) {
        enlarge_hash_table(bht);
    }

    for (len_t i = start; i < end; ++i) {
#if PARALLEL_HASHING
        hcm[i] = check_insert_in_hash_table(
                sht->ev[hcm[i]], sht->hd[hcm[i]].val, bht);
#else
        hcm[i] = insert_in_hash_table(sht->ev[hcm[i]], bht);
#endif
    }
}

static inline void insert_in_basis_hash_table_pivots(
    hm_t *row,
    ht_t *bht,
    const ht_t * const sht,
    const hi_t * const hcm,
    const md_t * const st
    )
{
    len_t l;

    /* while (bht->esz - bht->eld < row[LENGTH]) {
        enlarge_hash_table(bht);
    } */

    const len_t len = row[LENGTH]+OFFSET;
    const len_t evl = bht->evl;

    const hd_t * const hds    = sht->hd;
    exp_t * const * const evs = sht->ev;
    
    exp_t *evt  = (exp_t *)malloc(
        (unsigned long)(st->nthrds * evl) * sizeof(exp_t));
/* #if PARALLEL_HASHING
#pragma omp parallel for num_threads(st->nthrds) \
    private(l)
#endif */
    for (l = OFFSET; l < len; ++l) {
        exp_t *evtl = evt + (omp_get_thread_num() * evl);
        memcpy(evtl, evs[hcm[row[l]]],
                (unsigned long)evl * sizeof(exp_t));

#if PARALLEL_HASHING
        const val_t h = hds[hcm[row[l]]].val;
        row[l] = check_insert_in_hash_table(evtl, h, bht);
#else
        row[l] = insert_in_hash_table(evtl, bht);
#endif
    }
    free(evt);
}

static inline void insert_multiplied_poly_in_hash_table(
    hm_t *row,
    const val_t h1,
    const exp_t * const ea,
    const hm_t * const b,
    const ht_t * const ht1,
    ht_t *ht2
    )
{
    len_t j, l;
    exp_t *n;

    const len_t len = b[LENGTH]+OFFSET;
    const len_t evl = ht1->evl;

    exp_t * const *ev1      = ht1->ev;
    const hd_t * const hd1  = ht1->hd;
    
    exp_t **ev2     = ht2->ev;

    l = OFFSET;

    for (; l < len; ++l) {
        const exp_t * const eb = ev1[b[l]];

        n = ev2[ht2->eld];
        for (j = 0; j < evl; ++j) {
            n[j]  = (exp_t)(ea[j] + eb[j]);
        }

#if PARALLEL_HASHING
        const val_t h   = h1 + hd1[b[l]].val;
        row[l] = check_insert_in_hash_table(n, h, ht2);
#else
        row[l] = insert_in_hash_table(n, ht2);
#endif
    }
}

static inline void reinsert_in_hash_table(
    hm_t *row,
    exp_t * const *oev,
    ht_t *ht
    )
{
    hl_t i;
    hi_t k, pos;
    len_t j, l;
    exp_t *e;
    hd_t *d;
    val_t h;

    const len_t len = row[LENGTH]+OFFSET;
    const len_t evl = ht->evl;
    const hi_t hsz  = ht->hsz;
    /* ht->hsz <= 2^32 => mod is always uint32_t */
    const hi_t mod  = (hi_t)(hsz - 1);
    l = OFFSET;
letsgo:
    for (; l < len; ++l) {
        const exp_t * const n = oev[row[l]];
        /* generate hash value */
        h = 0;
        for (j = 0; j < evl; ++j) {
            h +=  ht->rn[j] * n[j];
        }
        k = h;
        i = 0;
restart:
        for (; i < hsz; ++i) {
            k = (hi_t)(k+i) & mod;
            const hi_t hm  = ht->hmap[k];
            if (!hm) {
                break;
            }
            if (ht->hd[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = ht->ev[hm];
            for (j = 0; j < evl-1; j += 2) {
                if (n[j] != ehm[j] || n[j+1] != ehm[j+1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[evl-1] != ehm[evl-1]) {
                i++;
                goto restart;
            }
            row[l] = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        ht->hmap[k] = pos = (hi_t)ht->eld;
        e = ht->ev[ht->eld];
        d = ht->hd + ht->eld;
        memcpy(e, n, (unsigned long)evl * sizeof(exp_t));
        d->sdm  =   generate_short_divmask(e, ht);
        d->deg  =   e[0];
        d->deg  +=  ht->ebl > 0 ? e[ht->ebl] : 0;
        d->val  =   h;

        ht->eld++;
        row[l] =  pos;
    }
}

void reset_hash_table_indices(
        ht_t *ht,
        const hi_t * const hcm,
        const len_t len
        )
{
    for (len_t i = 0; i < len; ++i) {
        ht->hd[hcm[i]].idx = 0;
    }
}

static void reset_hash_table(
    ht_t *ht,
    bs_t *bs,
    ps_t *psl,
    md_t *st
    )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i;
    hi_t k;
    exp_t *e;

    spair_t *ps = psl->p;
    exp_t **oev  = ht->ev;

    const len_t evl = ht->evl;
    const hl_t esz  = ht->esz;
    const bl_t bld  = bs->ld;
    const len_t pld = psl->ld;

    ht->ev  = calloc(esz, sizeof(exp_t *));
    if (ht->ev == NULL) {
        fprintf(stderr, "Computation needs too much memory on this machine,\n");
        fprintf(stderr, "cannot reset ht->ev, esz = %lu\n", (unsigned long)esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)evl * esz * sizeof(exp_t));
    if (tmp == NULL) {
        fprintf(stderr, "Computation needs too much memory on this machine,\n");
        fprintf(stderr, "resetting table failed, esz = %lu\n", (unsigned long)esz);
        fprintf(stderr, "segmentation fault will follow.\n");
    }
    for (k = 0; k < esz; ++k) {
        ht->ev[k]  = tmp + k*evl;
    }
    ht->eld = 1;
    memset(ht->hmap, 0, ht->hsz * sizeof(hi_t));
    memset(ht->hd, 0, esz * sizeof(hd_t));

    /* reinsert known elements */
    for (i = 0; i < bld; ++i) {
      if (bs->red[i] < 2) {
        reinsert_in_hash_table(bs->hm[i], oev, ht);
      }
    }
    for (i = 0; i < pld; ++i) {
        e = oev[ps[i].lcm];
#if PARALLEL_HASHING
        ps[i].lcm = check_insert_in_hash_table(e, 0, ht);
#else
        ps[i].lcm = insert_in_hash_table(e, ht);
#endif
    }
    /* note: all memory is allocated as a big block, so it is
     *       enough to free oev[0].       */
    free(oev[0]);
    free(oev);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->rht_ctime  +=  ct1 - ct0;
    st->rht_rtime  +=  rt1 - rt0;
}

/* computes lcm of a and b from ht1 and inserts it in ht2 */
static inline hi_t get_lcm(
    const hi_t a,
    const hi_t b,
    const ht_t *ht1,
    ht_t *ht2
    )
{
    len_t i;

    /* exponents of basis elements, thus from basis hash table */
    const exp_t * const ea = ht1->ev[a];
    const exp_t * const eb = ht1->ev[b];
    exp_t etmp[ht1->evl];
    const len_t evl = ht1->evl;
    const len_t ebl = ht1->ebl;

    /* set degree(s), if ebl == 0, i.e. we do not have an elimination block
     * order then the second for loop is just not executed and the third one
     * computes correctly the full degree of the lcm. */
    for (i = 1; i < evl; ++i) {
        etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
    }
    /* reset degree entries */
    etmp[0]   = 0;
    etmp[ebl] = 0;
    for (i = 1; i < ebl; ++i) {
        etmp[0]  += etmp[i];
    }
    for (i = ebl+1; i < evl; ++i) {
        etmp[ebl] += etmp[i];
    }
    /* printf("lcm -> ");
     * for (int ii = 0; ii < evl; ++ii) {
     *     printf("%d ", etmp[ii]);
     * }
     * printf("\n"); */
#if PARALLEL_HASHING
    return check_insert_in_hash_table(etmp, 0, ht2);
#else
    return insert_in_hash_table(etmp, ht2);
#endif
}

static inline hm_t *multiplied_poly_to_matrix_row(
    ht_t *sht,
    const ht_t *bht,
    const val_t hm,
    const exp_t * const em,
    const hm_t *poly
    )
{
  hm_t *row = (hm_t *)malloc((unsigned long)(poly[LENGTH]+OFFSET) * sizeof(hm_t));
  row[COEFFS]   = poly[COEFFS];
  row[PRELOOP]  = poly[PRELOOP];
  row[LENGTH]   = poly[LENGTH];
  /* hash table product insertions appear only here:
   * we check for hash table enlargements first and then do the insertions
   * without further elargment checks there */
  while (sht->eld+poly[LENGTH] >= sht->esz) {
    enlarge_hash_table(sht);
  }
  insert_multiplied_poly_in_hash_table(row, hm, em, poly, bht, sht);

  return row;
}
