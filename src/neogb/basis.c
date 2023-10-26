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


#include "basis.h"

static void free_basis_elements(
        bs_t *bs
        )
{
    len_t i, j, len;
    if (bs->cf_8) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_8[i]);
            bs->cf_8[i] = NULL;
            free(bs->hm[i]);
            bs->hm[i] = NULL;
        }
    }
    if (bs->cf_16) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_16[i]);
            bs->cf_16[i]  = NULL;
            free(bs->hm[i]);
            bs->hm[i] = NULL;
        }
    }
    if (bs->cf_32) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_32[i]);
            bs->cf_32[i]  = NULL;
            free(bs->hm[i]);
            bs->hm[i] = NULL;
        }
    }
    if (bs->cf_qq) {
        for (i = 0; i < bs->ld; ++i) {
            len = bs->hm[i][LENGTH];
            mpz_t *coeffs =  bs->cf_qq[bs->hm[i][COEFFS]];
            for (j = 0; j < len; ++j) {
                mpz_clear(coeffs[j]);
            }
            free(bs->cf_qq[bs->hm[i][COEFFS]]);
            bs->cf_qq[bs->hm[i][COEFFS]]  = NULL;
            free(bs->hm[i]);
            bs->hm[i] = NULL;
        }
    }
    /* signatures */
    free(bs->sm);
    bs->sm  =   NULL;
    free(bs->si);
    bs->si  =   NULL;

    bs->ld  = bs->lo  = bs->lml = 0;
}

void free_basis(
        bs_t **bsp
        )
{
    full_free_hash_table(&((*bsp)->ht));
    free_basis_without_hash_table(bsp);
}

void free_basis_without_hash_table(
        bs_t **bsp
        )
{
    len_t i, j, len;
    bs_t *bs  = *bsp;
    if (bs->cf_8) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_8[i]);
            free(bs->hm[i]);
        }
        free(bs->cf_8);
        bs->cf_8  = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    if (bs->cf_16) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_16[i]);
            free(bs->hm[i]);
        }
        free(bs->cf_16);
        bs->cf_16 = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    if (bs->cf_32) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_32[i]);
            free(bs->hm[i]);
        }
        free(bs->cf_32);
        bs->cf_32 = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    if (bs->cf_qq) {
        for (i = 0; i < bs->ld; ++i) {
            len = bs->hm[i][LENGTH];
            mpz_t *coeffs =  bs->cf_qq[bs->hm[i][COEFFS]];
            for (j = 0; j < len; ++j) {
                mpz_clear(coeffs[j]);
            }
            free(bs->cf_qq[bs->hm[i][COEFFS]]);
            free(bs->hm[i]);
        }
        free(bs->cf_qq);
        bs->cf_qq  = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    free(bs->lmps);
    bs->lmps  = NULL;
    free(bs->lm);
    bs->lm  = NULL;
    free(bs->red);
    bs->red = NULL;
    /* signatures */
    free(bs->sm);
    bs->sm  =   NULL;
    free(bs->si);
    bs->si  =   NULL;

    free(bs);
    bs    = NULL;
    *bsp  = bs;
}

bs_t *initialize_basis(
        md_t *md
        )
{
    bs_t *bs  = (bs_t *)calloc(1, sizeof(bs_t));
    /* initialize meta data */
    bs->lo        = 0;
    bs->ld        = 0;
    bs->lml       = 0;
    bs->constant  = 0;
    bs->sz        = md->init_bs_sz;
    bs->mltdeg    = 0;
    bs->ht        = initialize_basis_hash_table(md);

    /* initialize basis elements data */
    bs->hm    = (hm_t **)malloc((unsigned long)bs->sz * sizeof(hm_t *));
    bs->lm    = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));
    bs->lmps  = (bl_t *)malloc((unsigned long)bs->sz * sizeof(bl_t));
    bs->red   = (int8_t *)calloc((unsigned long)bs->sz, sizeof(int8_t));
    /* signature-based groebner basis computation? */
    if (md->use_signatures > 0) {
        bs->sm  =   (sm_t *)malloc((unsigned long)bs->sz * sizeof(sm_t));
        bs->si  =   (si_t *)malloc((unsigned long)bs->sz * sizeof(si_t));
    }
    /* initialize coefficients depending on ground field */
    switch (md->ff_bits) {
        case 8:
            bs->cf_8  = (cf8_t **)malloc((unsigned long)bs->sz * sizeof(cf8_t *));
            break;
        case 16:
            bs->cf_16  = (cf16_t **)malloc((unsigned long)bs->sz * sizeof(cf16_t *));
            break;
        case 32:
            bs->cf_32  = (cf32_t **)malloc((unsigned long)bs->sz * sizeof(cf32_t *));
            break;
        case 0:
            bs->cf_qq = (mpz_t **)malloc((unsigned long)bs->sz * sizeof(mpz_t *));
            break;
        default:
            exit(1);
    }

    return bs;
}

void check_enlarge_basis(
        bs_t *bs,
        const len_t added,
        const md_t * const st
        )
{
    if (bs->ld + added >= bs->sz) {
        bs->sz    = bs->sz * 2 > bs->ld + added ? bs->sz * 2 : bs->ld + added;
        bs->hm    = realloc(bs->hm, (unsigned long)bs->sz * sizeof(hm_t *));
        memset(bs->hm+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(hm_t *));
        bs->lm    = realloc(bs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
        memset(bs->lm+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(sdm_t));
        bs->lmps  = realloc(bs->lmps, (unsigned long)bs->sz * sizeof(bl_t));
        memset(bs->lmps+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(bl_t));
        bs->red   = realloc(bs->red, (unsigned long)bs->sz * sizeof(int8_t));
        memset(bs->red+bs->ld, 0,
                (unsigned long)(bs->sz-bs->ld) * sizeof(int8_t));

        switch (st->ff_bits) {
            case 8:
                bs->cf_8  = realloc(bs->cf_8,
                        (unsigned long)bs->sz * sizeof(cf8_t *));
                memset(bs->cf_8+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(cf8_t *));
                break;
            case 16:
                bs->cf_16  = realloc(bs->cf_16,
                        (unsigned long)bs->sz * sizeof(cf16_t *));
                memset(bs->cf_16+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(cf16_t *));
                break;
            case 32:
                bs->cf_32  = realloc(bs->cf_32,
                        (unsigned long)bs->sz * sizeof(cf32_t *));
                memset(bs->cf_32+bs->ld, 0, (unsigned long)(bs->sz-bs->ld) * sizeof(cf32_t *));
                break;
            case 0:
                bs->cf_qq = realloc(bs->cf_qq,
                        (unsigned long)bs->sz * sizeof(mpz_t *));
                break;
            default:
                exit(1);
        }
    }
}

/* finite field stuff  --  8 bit */
static inline void normalize_initial_basis_ff_8(
        bs_t *bs,
        const uint32_t fc
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    cf8_t **cf         = bs->cf_8;
    hm_t * const *hm    = bs->hm;
    const bl_t ld       = bs->ld;
    const int16_t fc8    = (int16_t)fc;

    for (i = 0; i < ld; ++i) {
        cf8_t *row  = cf[hm[i][COEFFS]];

        const uint8_t inv = mod_p_inverse_8((int16_t)row[0], (int16_t)fc8);
        const len_t os    = hm[i][PRELOOP];
        const len_t len   = hm[i][LENGTH];

        for (j = 0; j < os; ++j) {
            tmp1    =   ((int64_t)row[j] * inv) % fc8;
            tmp1    +=  (tmp1 >> 63) & fc8;
            row[j]  =   (cf8_t)tmp1;
        }
        for (j = os; j < len; j += UNROLL) {
            tmp1      =   ((int64_t)row[j] * inv) % fc8;
            tmp2      =   ((int64_t)row[j+1] * inv) % fc8;
            tmp3      =   ((int64_t)row[j+2] * inv) % fc8;
            tmp4      =   ((int64_t)row[j+3] * inv) % fc8;
            tmp1      +=  (tmp1 >> 63) & fc8;
            tmp2      +=  (tmp2 >> 63) & fc8;
            tmp3      +=  (tmp3 >> 63) & fc8;
            tmp4      +=  (tmp4 >> 63) & fc8;
            row[j]    =   (cf8_t)tmp1;
            row[j+1]  =   (cf8_t)tmp2;
            row[j+2]  =   (cf8_t)tmp3;
            row[j+3]  =   (cf8_t)tmp4;
        }
    }
}

/* finite field stuff  --  16 bit */
static inline void normalize_initial_basis_ff_16(
        bs_t *bs,
        const uint32_t fc
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    cf16_t **cf         = bs->cf_16;
    hm_t * const *hm    = bs->hm;
    const bl_t ld       = bs->ld;
    const int32_t fc16  = (int32_t)fc;

    for (i = 0; i < ld; ++i) {
        cf16_t *row = cf[hm[i][COEFFS]];

        const uint16_t inv  = mod_p_inverse_16((int32_t)row[0], (int32_t)fc16);
        const len_t os      = hm[i][PRELOOP];
        const len_t len     = hm[i][LENGTH];

        for (j = 0; j < os; ++j) {
            tmp1    =   ((int64_t)row[j] * inv) % fc16;
            tmp1    +=  (tmp1 >> 63) & fc16;
            row[j]  =   (cf16_t)tmp1;
        }
        for (j = os; j < len; j += UNROLL) {
            tmp1      =   ((int64_t)row[j] * inv) % fc16;
            tmp2      =   ((int64_t)row[j+1] * inv) % fc16;
            tmp3      =   ((int64_t)row[j+2] * inv) % fc16;
            tmp4      =   ((int64_t)row[j+3] * inv) % fc16;
            tmp1      +=  (tmp1 >> 63) & fc16;
            tmp2      +=  (tmp2 >> 63) & fc16;
            tmp3      +=  (tmp3 >> 63) & fc16;
            tmp4      +=  (tmp4 >> 63) & fc16;
            row[j]    =   (cf16_t)tmp1;
            row[j+1]  =   (cf16_t)tmp2;
            row[j+2]  =   (cf16_t)tmp3;
            row[j+3]  =   (cf16_t)tmp4;
        }
    }
}

/* finite field stuff  --  32 bit */
static inline void normalize_initial_basis_ff_32(
        bs_t *bs,
       const uint32_t fc
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    cf32_t **cf       = bs->cf_32;
    hm_t * const *hm  = bs->hm;
    const bl_t ld     = bs->ld;

    for (i = 0; i < ld; ++i) {
        cf32_t *row = cf[hm[i][COEFFS]];

        const uint32_t inv  = mod_p_inverse_32((int64_t)row[0], (int64_t)fc);
        const len_t os      = hm[i][PRELOOP]; 
        const len_t len     = hm[i][LENGTH]; 

        for (j = 0; j < os; ++j) {
            tmp1    =   ((int64_t)row[j] * inv) % fc;
            tmp1    +=  (tmp1 >> 63) & fc;
            row[j]  =   (cf32_t)tmp1;
        }
        for (j = os; j < len; j += UNROLL) {
            tmp1      =   ((int64_t)row[j] * inv) % fc;
            tmp2      =   ((int64_t)row[j+1] * inv) % fc;
            tmp3      =   ((int64_t)row[j+2] * inv) % fc;
            tmp4      =   ((int64_t)row[j+3] * inv) % fc;
            tmp1      +=  (tmp1 >> 63) & fc;
            tmp2      +=  (tmp2 >> 63) & fc;
            tmp3      +=  (tmp3 >> 63) & fc;
            tmp4      +=  (tmp4 >> 63) & fc;
            row[j]    =   (cf32_t)tmp1;
            row[j+1]  =   (cf32_t)tmp2;
            row[j+2]  =   (cf32_t)tmp3;
            row[j+3]  =   (cf32_t)tmp4;
        }
    }
}

/* characteristic zero stuff */
bs_t *copy_basis_mod_p(
        const bs_t * const gbs,
        const md_t *const st
        )
{
    len_t i, j, idx;

    /* set field characteristic */
    unsigned long prime = (unsigned long)st->fc;

    /* initialize basis */
    bs_t *bs        = (bs_t *)calloc(1, sizeof(bs_t));
    bs->lo          = gbs->lo;
    bs->ld          = gbs->ld;
    bs->lml         = gbs->lml;
    bs->sz          = gbs->sz;
    bs->constant    = gbs->constant;
    bs->ht          = gbs->ht;
    bs->hm          = (hm_t **)malloc((unsigned long)bs->sz * sizeof(hm_t *));
    bs->lm          = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));
    bs->lmps        = (bl_t *)malloc((unsigned long)bs->sz * sizeof(bl_t));
    bs->red         = (int8_t *)calloc((unsigned long)bs->sz, sizeof(int8_t));

    /* copy data */
    memcpy(bs->lm, gbs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
    memcpy(bs->lmps, gbs->lmps, (unsigned long)bs->sz * sizeof(bl_t));
    memcpy(bs->red, gbs->red, (unsigned long)bs->sz * sizeof(int8_t));
    if (st->use_signatures > 0) {
        memcpy(bs->sm, gbs->sm, (unsigned long)bs->sz * sizeof(sm_t));
        memcpy(bs->si, gbs->si, (unsigned long)bs->sz * sizeof(si_t));
    }
    /* copy monomials */
    for (i = 0; i < bs->ld; ++i) {
        bs->hm[i] =
            (hm_t *)malloc(((unsigned long)gbs->hm[i][LENGTH]+OFFSET) * sizeof(hm_t));
        memcpy(bs->hm[i], gbs->hm[i],
                ((unsigned long)gbs->hm[i][LENGTH]+OFFSET) * sizeof(hm_t));
    }
    /* copy coefficients */
    switch (st->ff_bits) {
        case 8:
            bs->cf_8    = (cf8_t **)malloc((unsigned long)bs->sz * sizeof(cf8_t *));
            for (i = 0; i < bs->ld; ++i) {
                idx = gbs->hm[i][COEFFS];
                bs->cf_8[idx]  =
                    (cf8_t *)malloc((unsigned long)(gbs->hm[i][LENGTH]) * sizeof(cf8_t));
                for (j = 0; j < gbs->hm[i][LENGTH]; ++j) {
                    bs->cf_8[idx][j] = (cf8_t)mpz_fdiv_ui(gbs->cf_qq[idx][j], prime);
                }
            }
            break;
        case 16:
            bs->cf_16   = (cf16_t **)malloc((unsigned long)bs->sz * sizeof(cf16_t *));
            for (i = 0; i < bs->ld; ++i) {
                idx = gbs->hm[i][COEFFS];
                bs->cf_16[idx]  =
                    (cf16_t *)malloc((unsigned long)(gbs->hm[i][LENGTH]) * sizeof(cf16_t));
                for (j = 0; j < gbs->hm[i][LENGTH]; ++j) {
                    bs->cf_16[idx][j] = (cf16_t)mpz_fdiv_ui(gbs->cf_qq[idx][j], prime);
                }
            }
            break;
        case 32:
            bs->cf_32   = (cf32_t **)malloc((unsigned long)bs->sz * sizeof(cf32_t *));
            for (i = 0; i < bs->ld; ++i) {
                idx = gbs->hm[i][COEFFS];
                bs->cf_32[idx]  =
                    (cf32_t *)malloc((unsigned long)(gbs->hm[i][LENGTH]) * sizeof(cf32_t));
                for (j = 0; j < gbs->hm[i][LENGTH]; ++j) {
                    bs->cf_32[idx][j] = (cf32_t)mpz_fdiv_ui(gbs->cf_qq[idx][j], prime);
                }
            }
            break;
        default:
            exit(1);
    }

    return bs;
}

void remove_content_of_initial_basis(
        bs_t *bs
        )
{
    len_t i, j;

    mpz_t **cf        = bs->cf_qq;
    hm_t * const *hm  = bs->hm;
    const bl_t ld     = bs->ld;

    mpz_t content;
    mpz_init(content);
    /* compute content, i.e. gcd of all coefficients */
    i = 0;
next_poly:
    for (; i < ld; ++i) {
        mpz_t *row = cf[hm[i][COEFFS]];
        mpz_set(content, row[0]);
        const len_t os  = hm[i][PRELOOP];
        const len_t len = hm[i][LENGTH];
        if (mpz_cmp_si(content, 0) != 0) {
            for (j = 1; j < len; ++j) {
                mpz_gcd(content, content, row[j]);
                if (mpz_cmp_si(content, 1) == 0) {
                    i++;
                    goto next_poly;
                }
            }
            /* remove content */
            for (j = 0; j < os; ++j) {
                mpz_divexact(row[j], row[j], content);
            }
            for (; j < len; j += UNROLL) {
                mpz_divexact(row[j], row[j], content);
                mpz_divexact(row[j+1], row[j+1], content);
                mpz_divexact(row[j+2], row[j+2], content);
                mpz_divexact(row[j+3], row[j+3], content);
            }
        }
    }
    mpz_clear(content);

    /* make lead coefficient positive */
    for (i = 0; i < ld; ++i) {
        mpz_t *row = cf[hm[i][COEFFS]];
        const len_t os  = hm[i][PRELOOP];
        const len_t len = hm[i][LENGTH];
        if (mpz_sgn(row[0]) == -1) {
            for (j = 0; j < os; ++j) {
                mpz_neg(row[j], row[j]);
            }
            for (; j < len; j += UNROLL) {
                mpz_neg(row[j], row[j]);
                mpz_neg(row[j+1], row[j+1]);
                mpz_neg(row[j+2], row[j+2]);
                mpz_neg(row[j+3], row[j+3]);
            }
        }
    }
}
