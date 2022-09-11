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

#ifndef GB_LA_H
#define GB_LA_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "config.h"
#include "types.h"

#if 0
static inline void copy_first_row_from_dense(mat_gb_block_t *bl,
    const mat_gb_block_t *obl, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  bl->len[0]  = 0;
  bl->len[1]  = bl->len[0];
  
  for (i=0; i<meta->bs; ++i) {
    if (obl->val[i] != 0) {
      bl->pos[bl->len[1]] = (bs_t)i;
      bl->val[bl->len[1]] = obl->val[i];
      bl->len[1]++;
    }
  }
}

static inline void copy_first_row_from_sparse(mat_gb_block_t *bl,
    const mat_gb_block_t *obl)
{
  bl->len[0]  = obl->len[0];
  bl->len[1]  = obl->len[1];
  memcpy(bl->pos, obl->pos, (bl->len[1]-bl->len[0]) * sizeof(bs_t));
  memcpy(bl->val, obl->val, (bl->len[1]-bl->len[0]) * sizeof(cf_t));
}

static inline void load_dense_row_for_update_from_dense(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  memset(dr, 0, meta->bs * sizeof(bf_t));

  for (i=0; i<meta->bs; ++i)
    dr[i]  = (bf_t)bl->val[idx*meta->bs+i];
}

static inline void load_dense_block_from_sparse(bf_t *db,
    const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  memset(db, 0, (meta->bs*meta->bs) * sizeof(bf_t));

  if (bl->len != NULL) {
    for (i=0; i<bl->nr; ++i) {
      for (j=bl->len[i]; j<bl->len[i+1]; ++j) {
        db[i*meta->bs + bl->pos[j]]  = (bf_t)bl->val[j];
      }
    }
  }
}

static inline void load_dense_block_for_update_full(bf_t *db,
    const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  memset(db, 0, meta->bs * meta->bs *sizeof(bf_t));

  if (bl->len != NULL) {
    for (i = 0; i < bl->nr; ++i) {
      for (j=bl->len[i]; j<bl->len[i+1]; ++j) {
        db[bl->pos[j] + i*meta->bs]  = (bf_t)bl->val[j];
      }
    }
  } else {
    if (bl->val == NULL)
      return;
    for (i = 0; i < bl->nr; ++i) {
      for (j=0; j<meta->bs; ++j) {
        db[j + i*meta->bs]  = (bf_t)bl->val[j + i*meta->bs];
      }
    }
  }
}

static inline void load_dense_block_for_update(bf_t *db,
    const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  memset(db, 0, meta->bs * bl->nr *sizeof(bf_t));

  if (bl->len != NULL) {
    for (i = 0; i < bl->nr; ++i) {
      for (j=bl->len[i]; j<bl->len[i+1]; ++j) {
        db[bl->pos[j] + i*meta->bs]  = (bf_t)bl->val[j];
      }
    }
  } else {
    for (i = 0; i < bl->nr; ++i) {
      for (j=0; j<meta->bs; ++j) {
        db[j + i*meta->bs]  = (bf_t)bl->val[j + i*meta->bs];
      }
    }
  }
}

static inline void load_dense_row_for_update(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  memset(dr, 0, meta->bs * sizeof(bf_t));

  if (bl->len != NULL) {
    for (i=bl->len[idx]; i<bl->len[idx+1]; ++i)
      dr[bl->pos[i]]  = (bf_t)bl->val[i];
  } else {
    if (bl->val != NULL) {
      for (i=0; i<meta->bs; ++i)
        dr[i]  = (bf_t)bl->val[idx*meta->bs+i];
    }
  }
}

static inline void load_dense_row_for_update_from_sparse(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  memset(dr, 0, meta->bs * sizeof(bf_t));

  if (bl->len != NULL) {
    for (i=bl->len[idx]; i<bl->len[idx+1]; ++i)
      dr[bl->pos[i]]  = (bf_t)bl->val[i];
  }
}

static inline void write_updated_block_to_sparse_format_directly(mat_gb_block_t *bl,
    const bf_t *db, const nelts_t lr, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;
  bf_t tmp;

  free(bl->len);
  bl->len = (nelts_t *)calloc((meta->bs+1), sizeof(nelts_t));
  bl->pos = realloc(bl->pos, meta->bs*meta->bs * sizeof(bs_t));
  bl->val= realloc(bl->val, meta->bs*meta->bs * sizeof(cf_t));

  for (i = 0; i < lr; ++i) {
    bl->len[i+1]  = bl->len[i];
    for (j = 0; j < meta->bs; ++j) {
      if (db[j + i*meta->bs] != 0) {
        tmp = db[j + i*meta->bs] % meta->mod;
        if (tmp != 0) {
          bl->pos[bl->len[i+1]] = (bs_t)j;
          bl->val[bl->len[i+1]] = (cf_t)tmp;
          bl->len[i+1]++;
        }
      }
    }
  }
}

static inline void write_updated_block_to_sparse_format(mat_gb_block_t *bl,
    const bf_t *db, const nelts_t lr, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;
  bf_t tmp;

  for (i = 0; i < lr; ++i) {
    bl->len[i+1]  = bl->len[i];
    for (j = 0; j < meta->bs; ++j) {
      if (db[j + i*meta->bs] != 0) {
        tmp = db[j + i*meta->bs] % meta->mod;
        if (tmp != 0) {
          bl->pos[bl->len[i+1]] = (bs_t)j;
          bl->val[bl->len[i+1]] = (cf_t)tmp;
          bl->len[i+1]++;
        }
      }
    }
  }
}

static inline void write_updated_row_to_sparse_format(mat_gb_block_t *bl,
    const bf_t *dr, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;
  bf_t tmp;

  bl->len[idx+1]  = bl->len[idx];
  for (i=0; i<meta->bs; ++i) {
    if (dr[i] != 0) {
      tmp = dr[i] % meta->mod;
      if (tmp != 0) {
        bl->pos[bl->len[idx+1]] = (bs_t)i;
        bl->val[bl->len[idx+1]] = (cf_t)tmp;
        bl->len[idx+1]++;
      }
    }
  }
}

static inline void set_updated_block(mat_gb_block_t **bl, mat_gb_block_t* new_bl)
{

  mat_gb_block_t *blin  = *bl;
  free(blin->len);
  free(blin->pos);
  free(blin->val);

  blin->len  = new_bl->len;
  blin->pos  = new_bl->pos;
  blin->val  = new_bl->val;

  free(new_bl);
}

static inline void update_dense_row_dense_rev_order_dense(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_block_t *mbl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=idx*meta->bs; i<(idx+1)*meta->bs; ++i) {
    if (mbl->val[i] != 0) {
      const bf_t mul    = (bf_t)meta->mod - mbl->val[i];
      const nelts_t ri  = bl->nr -1 - (i - idx*meta->bs);
      const cf_t *red   = bl->val+(ri*meta->bs);
      for (j=0; j<meta->bs; j=j+8) {
        dr[j]   +=  mul * red[j];
        dr[j+1] +=  mul * red[j+1];
        dr[j+2] +=  mul * red[j+2];
        dr[j+OFFSET] +=  mul * red[j+OFFSET];
        dr[j+4] +=  mul * red[j+4];
        dr[j+5] +=  mul * red[j+5];
        dr[j+6] +=  mul * red[j+6];
        dr[j+7] +=  mul * red[j+7];
      }
    }
  }
}
static inline void update_dense_row_dense_rev_order_sparse(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_block_t *mbl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  i=mbl->len[idx];
  
  if (mbl->len[idx+1] > 3) {
    for (; i<mbl->len[idx+1]-3; i=i+4) {
      const bf_t mul1   = (bf_t)meta->mod - mbl->val[i];
      const bf_t mul2   = (bf_t)meta->mod - mbl->val[i+1];
      const bf_t mul3   = (bf_t)meta->mod - mbl->val[i+2];
      const bf_t mul4   = (bf_t)meta->mod - mbl->val[i+OFFSET];
      const nelts_t ri1 = bl->nr -1 - mbl->pos[i];
      const nelts_t ri2 = bl->nr -1 - mbl->pos[i+1];
      const nelts_t ri3 = bl->nr -1 - mbl->pos[i+2];
      const nelts_t ri4 = bl->nr -1 - mbl->pos[i+OFFSET];
      const cf_t *red1  = bl->val+(ri1*meta->bs);
      const cf_t *red2  = bl->val+(ri2*meta->bs);
      const cf_t *red3  = bl->val+(ri3*meta->bs);
      const cf_t *red4  = bl->val+(ri4*meta->bs);
      for (j=0; j<meta->bs; j=j+8) {
        dr[j]   +=
          mul1 * red1[j] + mul2 * red2[j] + mul3 * red3[j] + mul4 * red4[j];
        dr[j+1] +=
          mul1 * red1[j+1] + mul2 * red2[j+1] + mul3 * red3[j+1] + mul4 * red4[j+1];
        dr[j+2] +=  
          mul1 * red1[j+2] + mul2 * red2[j+2] + mul3 * red3[j+2] + mul4 * red4[j+2];
        dr[j+OFFSET] += 
          mul1 * red1[j+OFFSET] + mul2 * red2[j+OFFSET] + mul3 * red3[j+OFFSET] + mul4 * red4[j+OFFSET];
        dr[j+4] += 
          mul1 * red1[j+4] + mul2 * red2[j+4] + mul3 * red3[j+4] + mul4 * red4[j+4];
        dr[j+5] += 
          mul1 * red1[j+5] + mul2 * red2[j+5] + mul3 * red3[j+5] + mul4 * red4[j+5];
        dr[j+6] += 
          mul1 * red1[j+6] + mul2 * red2[j+6] + mul3 * red3[j+6] + mul4 * red4[j+6];
        dr[j+7] += 
          mul1 * red1[j+7] + mul2 * red2[j+7] + mul3 * red3[j+7] + mul4 * red4[j+7];
      }
    }
  }
  for (; i<mbl->len[idx+1]; i=i+1) {
    const bf_t mul1   = (bf_t)meta->mod - mbl->val[i];
    const nelts_t ri1 = bl->nr -1 - mbl->pos[i];
    const cf_t *red1  = bl->val+(ri1*meta->bs);
    for (j=0; j<meta->bs; j=j+8) {
      dr[j]   +=  mul1 * red1[j];
      dr[j+1] +=  mul1 * red1[j+1];
      dr[j+2] +=  mul1 * red1[j+2];
      dr[j+OFFSET] +=  mul1 * red1[j+OFFSET];
      dr[j+4] +=  mul1 * red1[j+4];
      dr[j+5] +=  mul1 * red1[j+5];
      dr[j+6] +=  mul1 * red1[j+6];
      dr[j+7] +=  mul1 * red1[j+7];
    }
  }
}

static inline void update_dense_row_dense(bf_t *dr, const nelts_t idx,
    const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  /* printf("bl->val %p\n", bl->val);
   * for (nelts_t k=0; k< meta->bs; ++k) {
   *   printf("%p -- ",bl->val+(k*meta->bs));
   *   for (nelts_t l=0; l< meta->bs; ++l) {
   *     printf("%u ", bl->val[k*meta->bs+l]);
   *   }
   *   printf("\n");
   * }
   * printf("\n"); */

  /* printf("val %p, pos %p\n", mbl->val, mbl->pos); */
  /* printf("idx %u\n", idx);
   * for (int ii=0; ii<meta->bs; ++ii)
   *   printf("%lu ",dr[ii]);
   * printf("\n"); */
  i = mbl->len[idx];
  /* printf("len[idx+1] %u\n", mbl->len[idx+1]); */
  if (mbl->len[idx+1] > 3) {
#if 1
    /* printf("1\n"); */
    for (; i<mbl->len[idx+1]-3; i=i+4) {
      const bf_t mul1   = (bf_t)meta->mod - mbl->val[i];
      const nelts_t ri1 = mbl->pos[i];
      const bf_t mul2   = (bf_t)meta->mod - mbl->val[i+1];
      const nelts_t ri2 = mbl->pos[i+1];
      const bf_t mul3   = (bf_t)meta->mod - mbl->val[i+2];
      const nelts_t ri3 = mbl->pos[i+2];
      const bf_t mul4   = (bf_t)meta->mod - mbl->val[i+OFFSET];
      const nelts_t ri4 = mbl->pos[i+OFFSET];
      /* const nelts_t ri  = bl->nr -1 - mbl->pos[i]; */
      const cf_t *red1  = bl->val+(ri1*meta->bs);
      const cf_t *red2  = bl->val+(ri2*meta->bs);
      const cf_t *red3  = bl->val+(ri3*meta->bs);
      const cf_t *red4  = bl->val+(ri4*meta->bs);
      /*     printf("red %p\n",red);
       *     printf("ROW INDEX %u = %u - %u\n", ri, bl->nr, mbl->pos[i]);
       *
       *     printf("RED1: ");
       *     for (j=0; j<meta->bs; j++) {
       *       printf("%u ", red[j]);
       *     }
       *     printf("\n");
       *     printf("BS %u\n", meta->bs); */
      for (j=0; j<meta->bs; j=j+8) {
        /* printf("j %u\n", j);
         * printf("dr %lu\n", dr[j]);
         * printf("mul %lu\n", mul);
         * printf("red %u\n", red[j]); */
        dr[j]   +=  mul1 * red1[j] + mul2 * red2[j] + mul3 * red3[j] + mul4 * red4[j];
        dr[j+1] +=  mul1 * red1[j+1] + mul2 * red2[j+1] + mul3 * red3[j+1] + mul4 * red4[j+1];
        dr[j+2] +=  mul1 * red1[j+2] + mul2 * red2[j+2] + mul3 * red3[j+2] + mul4 * red4[j+2];
        dr[j+OFFSET] +=  mul1 * red1[j+OFFSET] + mul2 * red2[j+OFFSET] + mul3 * red3[j+OFFSET] + mul4 * red4[j+OFFSET];
        dr[j+4] +=  mul1 * red1[j+4] + mul2 * red2[j+4] + mul3 * red3[j+4] + mul4 * red4[j+4];
        dr[j+5] +=  mul1 * red1[j+5] + mul2 * red2[j+5] + mul3 * red3[j+5] + mul4 * red4[j+5];
        dr[j+6] +=  mul1 * red1[j+6] + mul2 * red2[j+6] + mul3 * red3[j+6] + mul4 * red4[j+6];
        dr[j+7] +=  mul1 * red1[j+7] + mul2 * red2[j+7] + mul3 * red3[j+7] + mul4 * red4[j+7];
        /* printf("%lu\n", dr[bl->pos[j]]); */
      }
    }
  }
#endif
  for (; i<mbl->len[idx+1]; ++i) {
    const bf_t mul    = (bf_t)meta->mod - mbl->val[i];
    const nelts_t ri  = mbl->pos[i];
    /* const nelts_t ri  = bl->nr -1 - mbl->pos[i]; */
    const cf_t *red   = bl->val+(ri*meta->bs);
/*     printf("red %p\n",red);
 *     printf("ROW INDEX %u = %u - %u\n", ri, bl->nr, mbl->pos[i]);
 *
 *     printf("RED1: ");
 *     for (j=0; j<meta->bs; j++) {
 *       printf("%u ", red[j]);
 *     }
 *     printf("\n");
 *     printf("BS %u\n", meta->bs); */
    for (j=0; j<meta->bs; j=j+8) {
      /* printf("j %u\n", j);
       * printf("dr %lu\n", dr[j]);
       * printf("mul %lu\n", mul);
       * printf("red %u\n", red[j]); */
      dr[j]   +=  mul * red[j];
      dr[j+1] +=  mul * red[j+1];
      dr[j+2] +=  mul * red[j+2];
      dr[j+OFFSET] +=  mul * red[j+OFFSET];
      dr[j+4] +=  mul * red[j+4];
      dr[j+5] +=  mul * red[j+5];
      dr[j+6] +=  mul * red[j+6];
      dr[j+7] +=  mul * red[j+7];
      /* printf("%lu\n", dr[bl->pos[j]]); */
    }
  /*   printf("ri %u\n", ri);
   * for (int ii=0; ii<meta->bs; ++ii)
   *   printf("%lu ",dr[ii]);
   * printf("\n"); */
  }
}

static inline void update_dense_block_sparse(bf_t *db,
    const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
    const nelts_t max, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j, k;

  register bf_t mul;
  register nelts_t ri;
  register nelts_t shift;

  for (i=0; i<max; ++i) {
    for (j=mbl->len[i]; j<mbl->len[i+1]; ++j) {
      mul   = (uint32_t)meta->mod - mbl->val[j];
      ri    = mbl->pos[j];
      shift = (bl->len[ri+1]-bl->len[ri]) % 8;
      for (k=bl->len[ri]; k<bl->len[ri]+shift; ++k) {
      /* for (k=bl->len[ri]; k<bl->len[ri+1]; ++k) { */
        db[i*meta->bs + bl->pos[k]] +=  mul * bl->val[k];
      }
      for (; k<bl->len[ri+1]; k=k+8) {
        db[i*meta->bs + bl->pos[k]]    +=  mul * bl->val[k];
        db[i*meta->bs + bl->pos[k+1]]  +=  mul * bl->val[k+1];
        db[i*meta->bs + bl->pos[k+2]]  +=  mul * bl->val[k+2];
        db[i*meta->bs + bl->pos[k+OFFSET]]  +=  mul * bl->val[k+OFFSET];
        db[i*meta->bs + bl->pos[k+4]]  +=  mul * bl->val[k+4];
        db[i*meta->bs + bl->pos[k+5]]  +=  mul * bl->val[k+5];
        db[i*meta->bs + bl->pos[k+6]]  +=  mul * bl->val[k+6];
        db[i*meta->bs + bl->pos[k+7]]  +=  mul * bl->val[k+7];
      }
    }
  }
}

static inline void update_dense_block_sparse_self(bf_t *db,
    mat_gb_block_t *bl, const mat_gb_block_t *mbl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i, j, k;

  for (i=0; i<bl->nr; ++i) {
    for (j=mbl->len[i]; j<mbl->len[i+1]; ++j) {
      const bf_t mul    = (bf_t)meta->mod - mbl->val[j];
      /* printf("mul %u\n", mul); */
      const nelts_t ri  = bl->nr -1 - mbl->pos[j];
      const nelts_t shift = (bl->len[ri+1]-bl->len[ri]) % 4;
      /* printf("bl->len[ri] %u | bl->len[ri+1] %u | shift %u\n",
      *     bl->len[ri], bl->len[ri+1], shift); */
      for (k=bl->len[ri]; k<bl->len[ri]+shift; ++k) {
        db[i*meta->bs + bl->pos[k]] +=  mul * bl->val[k];
      }
      for (; k<bl->len[ri+1]; k=k+4) {
        db[i*meta->bs + bl->pos[k]]    +=  mul * bl->val[k];
        db[i*meta->bs + bl->pos[k+1]]  +=  mul * bl->val[k+1];
        db[i*meta->bs + bl->pos[k+2]]  +=  mul * bl->val[k+2];
        db[i*meta->bs + bl->pos[k+OFFSET]]  +=  mul * bl->val[k+OFFSET];
      }
    }
    /* for (int ii=0; ii<meta->bs; ++ii)
     *   printf("%lu ",db[i*meta->bs+ii]);
     * printf("\n"); */

    write_updated_row_to_sparse_format(bl, db+i*meta->bs, i, meta);
  }
}

static inline void update_dense_row_sparse_rev_order_dense(bf_t *dr,
  const nelts_t idx, const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
  const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=idx*meta->bs; i<(idx+1)*meta->bs; ++i) {
    if (mbl->val[i] != 0) {
      const bf_t mul      = (bf_t)meta->mod - mbl->val[i];
      const nelts_t ri    = bl->nr -1 - (i - idx*meta->bs);
      const nelts_t shift = (bl->len[ri+1]-bl->len[ri]) % 8;

      /* printf("ri %u | idx %u | shift %u | bl->nr %u | len[ri] %u | len[ri+1] %u\n", ri, idx, */
      for (j=bl->len[ri]; j<bl->len[ri]+shift; ++j) {
        dr[bl->pos[j]] +=  mul * bl->val[j];
      }
      for (; j<bl->len[ri+1]; j=j+8) {
        dr[bl->pos[j]]    +=  mul * bl->val[j];
        dr[bl->pos[j+1]]  +=  mul * bl->val[j+1];
        dr[bl->pos[j+2]]  +=  mul * bl->val[j+2];
        dr[bl->pos[j+OFFSET]]  +=  mul * bl->val[j+OFFSET];
        dr[bl->pos[j+4]]  +=  mul * bl->val[j+4];
        dr[bl->pos[j+5]]  +=  mul * bl->val[j+5];
        dr[bl->pos[j+6]]  +=  mul * bl->val[j+6];
        dr[bl->pos[j+7]]  +=  mul * bl->val[j+7];
      }
    }
  }
}

static inline void update_dense_row_sparse_rev_order_sparse(bf_t *dr,
  const nelts_t idx, const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
  const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=mbl->len[idx]; i<mbl->len[idx+1]; ++i) {
    const bf_t mul      = (bf_t)meta->mod - mbl->val[i];
    const nelts_t ri    = bl->nr -1 - mbl->pos[i];
    const nelts_t shift = (bl->len[ri+1]-bl->len[ri]) % 8;

    for (j=bl->len[ri]; j<bl->len[ri]+shift; ++j) {
      dr[bl->pos[j]] +=  mul * bl->val[j];
    }
    for (; j<bl->len[ri+1]; j=j+8) {
      dr[bl->pos[j]]    +=  mul * bl->val[j];
      dr[bl->pos[j+1]]  +=  mul * bl->val[j+1];
      dr[bl->pos[j+2]]  +=  mul * bl->val[j+2];
      dr[bl->pos[j+OFFSET]]  +=  mul * bl->val[j+OFFSET];
      dr[bl->pos[j+4]]  +=  mul * bl->val[j+4];
      dr[bl->pos[j+5]]  +=  mul * bl->val[j+5];
      dr[bl->pos[j+6]]  +=  mul * bl->val[j+6];
      dr[bl->pos[j+7]]  +=  mul * bl->val[j+7];
    }
  }
}

static inline void update_dense_row_sparse(bf_t *dr, const nelts_t idx,
    const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=mbl->len[idx]; i<mbl->len[idx+1]; ++i) {
    const bf_t mul    = (bf_t)meta->mod - mbl->val[i];
    const nelts_t ri  = mbl->pos[i];
    const nelts_t shift = (bl->len[ri+1]-bl->len[ri]) % 16;
    j=bl->len[ri];
    for (; j<bl->len[ri]+shift; ++j) {
      dr[bl->pos[j]] +=  mul * bl->val[j];
    }
    for (; j<bl->len[ri+1]; j=j+16) {
      dr[bl->pos[j]]    +=  mul * bl->val[j];
      dr[bl->pos[j+1]]  +=  mul * bl->val[j+1];
      dr[bl->pos[j+2]]  +=  mul * bl->val[j+2];
      dr[bl->pos[j+OFFSET]]  +=  mul * bl->val[j+OFFSET];
      dr[bl->pos[j+4]]  +=  mul * bl->val[j+4];
      dr[bl->pos[j+5]]  +=  mul * bl->val[j+5];
      dr[bl->pos[j+6]]  +=  mul * bl->val[j+6];
      dr[bl->pos[j+7]]  +=  mul * bl->val[j+7];
      dr[bl->pos[j+8]]  +=  mul * bl->val[j+8];
      dr[bl->pos[j+9]]  +=  mul * bl->val[j+9];
      dr[bl->pos[j+10]]  +=  mul * bl->val[j+10];
      dr[bl->pos[j+11]]  +=  mul * bl->val[j+11];
      dr[bl->pos[j+12]]  +=  mul * bl->val[j+12];
      dr[bl->pos[j+13]]  +=  mul * bl->val[j+13];
      dr[bl->pos[j+14]]  +=  mul * bl->val[j+14];
      dr[bl->pos[j+15]]  +=  mul * bl->val[j+15];
    }
  }
}

static inline void update_single_block_dense_rev_order(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  if (fbl->len != NULL) {
    for (i=0; i<ubl->nr; ++i) {
      /* load dense row to be updated */
      load_dense_row_for_update_from_dense(dr, i, ubl, meta);

      /* find corresponding row and multiplier */
      update_dense_row_sparse_rev_order_sparse(dr, i, nbl, fbl, meta);

      /* write updated row to new storage holders */
      write_updated_row_to_sparse_format(nbl, dr, i, meta);
    }
  } else {
    for (i=0; i<ubl->nr; ++i) {
      /* load dense row to be updated */
      load_dense_row_for_update_from_dense(dr, i, ubl, meta);

      /* find corresponding row and multiplier */
      update_dense_row_sparse_rev_order_dense(dr, i, nbl, fbl, meta);

      /* write updated row to new storage holders */
      write_updated_row_to_sparse_format(nbl, dr, i, meta);
    }
  }
  free(dr);
  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_dense(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  for (i=1; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_dense(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row_sparse(dr, i, nbl, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  free(dr);
  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_sparse_2(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *db  = (bf_t *)malloc((meta->bs*meta->bs) * sizeof(bf_t));

  load_dense_block_from_sparse(db, ubl, meta);

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);
  update_dense_block_sparse(db, nbl, fbl, ubl->nr, meta);
  free(db);

  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_sparse_rev_order(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  for (i=0; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_sparse(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    /* printf("UPDATE ROW %u\n", i); */
    update_dense_row_sparse_rev_order_sparse(dr, i, nbl, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }
  free(dr);

  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_sparse(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  for (i=0; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_sparse(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    /* printf("UPDATE ROW %u\n", i); */
    update_dense_row_sparse(dr, i, nbl, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }
  free(dr);

  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_rev_order(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  if (mat[idx].len != NULL) {
    update_single_block_sparse_rev_order(mat, shift, idx, meta);
  } else {
    if (mat[idx].val != NULL)
      update_single_block_dense_rev_order(mat, shift, idx, meta);
    else
      return;
  }
}

static inline void update_single_block(mat_gb_block_t *mat,
    const nelts_t shift, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  if (mat[idx].len != NULL) {
    update_single_block_sparse(mat, shift, idx, meta);
  } else {
    if (mat[idx].val != NULL)
      update_single_block_dense(mat, shift, idx, meta);
    else
      return;
  }
}

static inline void update_upper_row_block(mat_gb_block_t *mat,
    const nelts_t shift, const mat_gb_meta_data_t *meta, const int t)
{
  /* printf("mat %p\n", mat);
   * printf("ncb %u\n", meta->ncb); */
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* the first block is used for updated the remaining ones,
       * i.e. we start at block index i=shift+1 */
      /* printf("shift %u +1 < %u ncb?\n", shift, meta->ncb); */
      for (nelts_t i=meta->ncb_AC; i<meta->ncb; ++i) {
#pragma omp task
        update_single_block_rev_order(mat, shift, i, meta);
      }
    }
  }
  /* check density of blocks */
  /* adjust_block_row_types(mat, meta); */
  adjust_block_row_types_including_dense_righthand(mat, meta);
}

static inline void update_block(mat_gb_block_t **bl, const bf_t *db,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;
  const nelts_t bs_square = meta->bs * meta->bs;
  bf_t tmp;

  mat_gb_block_t *blin  = *bl;
  blin->len     = realloc(blin->len, (meta->bs+1) * sizeof(nelts_t));
  blin->pos     = realloc(blin->pos, bs_square * sizeof(bs_t));
  blin->val     = realloc(blin->val, bs_square * sizeof(cf_t));
  blin->len[0]  = 0;

  for (i=0; i<meta->bs; ++i) {
    blin->len[i+1]  = blin->len[i];
    for (j=0; j<meta->bs; ++j) {
      tmp = db[i*meta->bs+j] %  meta->mod;
      if (tmp != 0) {
        blin->pos[blin->len[i+1]] = (bs_t)j;
        blin->val[blin->len[i+1]] = (cf_t)tmp;
        blin->len[i+1]++;
      }
    }
  }
  if (blin->len[meta->bs] == 0) {
    free(blin->len);
    blin->len = NULL;
    free(blin->pos);
    blin->pos = NULL;
    free(blin->val);
    blin->val = NULL;
  } else {
    blin->pos = realloc(blin->pos, blin->len[meta->bs] * sizeof(bs_t));
    blin->val = realloc(blin->val, blin->len[meta->bs] * sizeof(cf_t));
  }
}

static inline void sparse_update_lower_block_by_upper_block_2(mat_gb_block_t *l,
    const mat_gb_block_t *u, const nelts_t shift, const nelts_t rbi,
    const nelts_t cbi, const mat_gb_meta_data_t *meta)
{
  /* first multiplier block used as lookup table for updates */
  const mat_gb_block_t *mbl  = l + rbi*meta->ncb+shift;
  /* block to be updated */
  /* printf("rbi %u | cbi %u | ncb %u | nrb %u | shift %u\n",
   *     rbi, cbi, meta->ncb, meta->nrb_CD, shift); */
  mat_gb_block_t *ubl  = l + rbi*meta->ncb+cbi;
  /* row to be updated in dense format */
  bf_t *db  = (bf_t *)malloc((meta->bs*meta->bs) * sizeof(bf_t));

  load_dense_block_from_sparse(db, ubl, meta);

  update_dense_block_sparse(db, u+cbi, mbl, ubl->nr, meta);

  update_block(&ubl, db, meta);
  free(db);

}

static inline void update_upper_blocks(const mat_gb_block_t *mul,
    const mat_gb_block_t *red, mat_gb_block_t *ubl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  if (red->len != NULL) {
    for (i=0; i<ubl->nr; ++i) {
      /* load dense row to be updated */
      load_dense_row_for_update(dr, i, ubl, meta);

      update_dense_row_sparse_rev_order_sparse(dr, i, red, mul, meta);

      /* write updated row to new storage holders */
      write_updated_row_to_sparse_format(nbl, dr, i, meta);
    }
  } else {
    for (i=0; i<ubl->nr; ++i) {
      /* load dense row to be updated */
      load_dense_row_for_update(dr, i, ubl, meta);

      update_dense_row_dense_rev_order_sparse(dr, i, red, mul, meta);

      /* write updated row to new storage holders */
      write_updated_row_to_sparse_format(nbl, dr, i, meta);
    }
  }

  set_updated_block(&ubl, nbl);

  free(dr);
}

static inline void sparse_update_lower_block_by_upper_block_rev_order(
    mat_gb_block_t *l, const mat_gb_block_t *u, const nelts_t shift,
    const nelts_t rbi, const nelts_t cbi, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first multiplier block used as lookup table for updates */
  const mat_gb_block_t *mbl  = l + rbi*meta->ncb+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = l + rbi*meta->ncb+cbi;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  /* printf("cbi %u | rbi %u\n", cbi, rbi); */
  for (i=0; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update(dr, i, ubl, meta);
    /* load_dense_row_for_update_from_sparse(dr, i, ubl, meta); */

    /* find corresponding row and multiplier */

    /* we know already that u[cbi] is not the empty block,
     * so it is either sparse or dense */
    if (u[cbi].len != NULL) {
      update_dense_row_sparse_rev_order_sparse(dr, i, u+cbi, mbl, meta);
    } else {
      update_dense_row_dense_rev_order_sparse(dr, i, u+cbi, mbl, meta);
    }

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);
  
  free(dr);
}

static inline void sparse_update_lower_block_by_upper_block_at_once(mat_gb_block_t *l,
    const mat_gb_block_t *u, const nelts_t shift, const nelts_t rbi,
    const nelts_t cbi, const mat_gb_meta_data_t *meta)
{
  /* first multiplier block used as lookup table for updates */
  const mat_gb_block_t *mbl  = l + rbi*meta->ncb+shift;
  mat_gb_block_t *ubl  = l + rbi*meta->ncb+cbi;
  /* row to be updated in dense format */
  bf_t *db  = (bf_t *)malloc(meta->bs * ubl->nr * sizeof(bf_t));

  load_dense_block_for_update(db, ubl, meta);

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  /* printf("cbi %u | rbi %u\n", cbi, rbi); */
  if (u[cbi].len != NULL) {
    /* printf("sparse\n"); */
    update_dense_block_sparse(db, u+cbi, mbl, ubl->nr, meta);
#if 0
  } else {
    /* printf("dense\n"); */
    update_dense_block_dense(db, u+cbi, mbl, ubl->nr, meta);
#endif
  }

  /* write updated row to new storage holders */
  write_updated_block_to_sparse_format(nbl, db, ubl->nr, meta);

  set_updated_block(&ubl, nbl);

  free(db);
}

static inline void sparse_update_lower_block_by_upper_block(mat_gb_block_t *l,
    const mat_gb_block_t *u, const nelts_t shift, const nelts_t rbi,
    const nelts_t cbi, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first multiplier block used as lookup table for updates */
  const mat_gb_block_t *mbl  = l + rbi*meta->ncb+shift;
  /* block to be updated */
  mat_gb_block_t *ubl  = l + rbi*meta->ncb+cbi;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  for (i=0; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update(dr, i, ubl, meta);
    /* load_dense_row_for_update_from_sparse(dr, i, ubl, meta); */

    /* find corresponding row and multiplier */

    /* we know already that u[cbi] is not the empty block,
      * so it is either sparse or dense */
    update_dense_row_sparse(dr, i, u+cbi, mbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);

  free(dr);
}

static inline void sparse_update_D_blocks_by_B_blocks(
    mat_gb_block_t *l, mat_gb_block_t **u, const nelts_t rbi,
    const nelts_t cbi, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  for (i = meta->ncb_AC; i < meta->ncb; ++i) {
    if (u[cbi][i].val != NULL) {
      sparse_update_lower_block_by_upper_block_rev_order(
          l, u[cbi], cbi, rbi, i, meta);
    }
  }
  /* free multiplier block from C */
  free_mat_gb_block(l + cbi + rbi*meta->ncb);
  adjust_block_row_types_including_dense_righthand(l + rbi*meta->ncb, meta);
}

static inline void compute_multiples_in_first_block_at_once(bf_t *db,
    const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j, l;

  for (l = 0; l < meta->bs; ++l) {
    for (i=0; i<bl->nr; ++i) {
      if (db[i+l*meta->bs] != 0) {
        db[i+l*meta->bs] = db[i+l*meta->bs] % meta->mod;
        if (db[i+l*meta->bs] != 0) {
          const bf_t mul  = (bf_t)meta->mod - db[i+l*meta->bs];
          const nelts_t ri  = i;
          const nelts_t shift = (bl->len[ri+1]-bl->len[ri]-1) % 4;
          for (j=bl->len[ri]+1; j<bl->len[ri]+shift+1; ++j) {
            db[bl->pos[j]+l*meta->bs] +=  mul * bl->val[j];
          }
          for (; j<bl->len[ri+1]; j=j+4) {
            db[bl->pos[j]+l*meta->bs] +=  mul * bl->val[j];
            db[bl->pos[j+1]+l*meta->bs] +=  mul * bl->val[j+1];
            db[bl->pos[j+2]+l*meta->bs] +=  mul * bl->val[j+2];
            db[bl->pos[j+OFFSET]+l*meta->bs] +=  mul * bl->val[j+OFFSET];
          }
        }
      }
    }
  }
}

static inline void compute_multiples_in_first_block(bf_t *dr,
    const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=0; i<bl->nr; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % meta->mod;
      if (dr[i] != 0) {
        const bf_t mul  = (bf_t)meta->mod - dr[i];
        const nelts_t ri  = i;
        const nelts_t shift = (bl->len[ri+1]-bl->len[ri]-1) % 4;
        for (j=bl->len[ri]+1; j<bl->len[ri]+shift+1; ++j) {
          dr[bl->pos[j]] +=  mul * bl->val[j];
        }
        for (; j<bl->len[ri+1]; j=j+4) {
          dr[bl->pos[j]] +=  mul * bl->val[j];
          dr[bl->pos[j+1]] +=  mul * bl->val[j+1];
          dr[bl->pos[j+2]] +=  mul * bl->val[j+2];
          dr[bl->pos[j+OFFSET]] +=  mul * bl->val[j+OFFSET];
        }
      }
    }
  }
}

static inline void sparse_update_first_block_by_upper_block(mat_gb_block_t *l,
    const mat_gb_block_t *u, const nelts_t shift, const nelts_t rbi,
   const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* block to be updated */
  mat_gb_block_t *ubl  = l + rbi*meta->ncb+shift;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = generate_mat_gb_block(meta, ubl->nr);

  for (i=0; i<ubl->nr; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update(dr, i, ubl, meta);
    /* load_dense_row_for_update_from_sparse(dr, i, ubl, meta); */

    /* find corresponding row and multiplier */

    /* we know already that u[cbi] is not the empty block,
     * so it is either sparse or dense */
    compute_multiples_in_first_block(dr, u+shift, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);
  
  free(dr);
}

static inline void update_lower_by_upper_row_block_not_prereduced(
    mat_gb_block_t *l, const mat_gb_block_t *u, const nelts_t shift,
    const mat_gb_meta_data_t *meta, const int t)
{
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* the first block is used for updated the remaining ones,
       * i.e. we start at block index i=1 */
      for (nelts_t i=0; i<meta->nrb_CD; ++i) {
        /* need to look at the first block in each row in
        * order to decide which algorithm to be chosen */
        if (l[i*meta->ncb+shift].val != NULL) {
          sparse_update_first_block_by_upper_block(l, u, shift, i, meta);
          for (nelts_t j=shift+1; j<meta->ncb; ++j) {
            if (u[j].val != NULL) {
              /* reduce first block */
#pragma omp task
              {
                /* sparse_update_lower_block_by_upper_block_at_once(l, u, shift, i, j, meta); */
                sparse_update_lower_block_by_upper_block(l, u, shift, i, j, meta);
              }
            }
          }
#pragma omp taskwait
          free_mat_gb_block(l+i*meta->ncb+shift);
          adjust_block_row_types_including_dense_righthand(l+i*meta->ncb, meta);
        }
      }
    }
  }
}

static inline void reduce_upper_column_block(mat_gb_block_t **l,
    const nelts_t shift, const mat_gb_meta_data_t *meta, const int t)
{
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* we only have to reduce B */
      for (nelts_t i=0; i<shift; ++i) {
        if (l[i][shift].val != NULL) {
          /* need to look at the first block in each row in
           * order to decide which algorithm to be chosen */
          for (nelts_t j=meta->ncb_AC; j<meta->ncb; ++j) {
            if (l[shift][j].val != NULL) {
#pragma omp task
              {
                update_upper_blocks(l[i]+shift, l[shift]+j, l[i]+j, meta);
              }
            }
          }
#pragma omp taskwait
          free_mat_gb_block(l[i]+shift);
          adjust_block_row_types_including_dense_righthand(l[i], meta);
          /* adjust_block_row_types_including_dense(l[i], meta); */
        }
      }
    }
  }
}

static inline void update_lower_by_unreduced_upper_row_block(mat_gb_block_t *l,
    const mat_gb_block_t *u, const mat_gb_meta_data_t *meta, const int t)
{
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* loop over row blocks in CD */
      for (nelts_t i = 0; i < meta->nrb_CD; ++i) {
#pragma omp task
        {
          /* allocate memory for a dense block:
           * always the block to be updated */
          bf_t *db  = (bf_t *)malloc(meta->bs * meta->bs * sizeof(bf_t));

          /* multiply with all previously updated blocks from A multiplied by C */
          for (nelts_t j = 0; j < meta->ncb_AC; ++j) {
            load_dense_block_for_update_full(db, l+j+i*meta->ncb, meta);
            for (nelts_t k = 0; k < j; ++k) {
              if (u[k*meta->ncb + j].val != NULL && l[k+i*meta->ncb].val != NULL) {
                update_dense_block_sparse(db, u + k*meta->ncb + j,
                    l + i*meta->ncb + k, meta->bs, meta);
              }
            }
            /* finally update block computing needed multiples for updating D later on */
            compute_multiples_in_first_block_at_once(db, u+j*meta->ncb + j, meta);
            write_updated_block_to_sparse_format_directly(l+j+i*meta->ncb,
                db, l[j+i*meta->ncb].nr, meta);
          }

          /* now update D */
          for (nelts_t j = meta->ncb_AC; j < meta->ncb; ++j) {
            load_dense_block_for_update_full(db, l+j+i*meta->ncb, meta);
            for (nelts_t k = 0; k < meta->nrb_AB; ++k) {
              if (u[k*meta->ncb + j].val != NULL && l[k+i*meta->ncb].val != NULL) {
                update_dense_block_sparse(db, u + k*meta->ncb + j,
                    l + i*meta->ncb + k, meta->bs, meta);
              }
            }
            write_updated_block_to_sparse_format_directly(l+j+i*meta->ncb,
                db, l[j+i*meta->ncb].nr, meta);
          }
          /* free_mat_gb_block(l + j + i * meta->ncb); */
          adjust_block_row_types_including_dense_righthand(l + i * meta->ncb, meta);

          free(db);
        }
#pragma omp taskwait
      }
    }
  }
}

static inline void update_lower_by_reduced_upper_row_block(mat_gb_block_t *l,
    mat_gb_block_t **u, const mat_gb_meta_data_t *meta, const int t)
{
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      for (nelts_t i = 0; i < meta->ncb_AC; ++i) {
        for (nelts_t j = 0; j < meta->nrb_CD; ++j) {
          if (l[i + j*meta->ncb].val != 0) {
#pragma omp task
            {
              sparse_update_D_blocks_by_B_blocks(l, u, j, i, meta);
            }
          }
        }
#pragma omp taskwait
      }
    }
  }
}

static inline void update_lower_by_upper_row_block(mat_gb_block_t *l,
    mat_gb_block_t *u, const nelts_t shift,
    const mat_gb_meta_data_t *meta, const int t)
{
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* the first block is used for updated the remaining ones,
       * i.e. we start at block index i=1 */
      for (nelts_t j=shift+1; j<meta->ncb; ++j) {
        if (u[j].val != NULL) {
          for (nelts_t i=0; i<meta->nrb_CD; ++i) {
            /* need to look at the first block in each row in
            * order to decide which algorithm to be chosen */
            if (l[i*meta->ncb+shift].val != NULL) {
#pragma omp task
              {
                sparse_update_lower_block_by_upper_block(l, u, shift, i, j, meta);
              }
            }
          }
        } else {
          continue;
        }
      }
    }
  }
    /* remove the first block in each block row */
    for (nelts_t i=0; i<meta->nrb_CD; ++i) {
      free_mat_gb_block(l+i*meta->ncb+shift);
    }

/* at the moment we only work with sparse blocks */
#if 1
#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      for (nelts_t i=0; i<meta->nrb_CD; ++i) {
#pragma omp task
        {
          /* check density of blocks */
          adjust_block_row_types(l+i*meta->ncb, meta);
        }
      }
    }
  }
#endif
}

static inline void reduce_upper_part(mat_gb_block_t **mat,
    const mat_gb_meta_data_t *meta, const int t)
{
  for (nelts_t j=0; j<meta->nrb_AB; ++j) {
    /* if first block is not already unit matrix block
     * we have to update the block row */
    if (mat[meta->nrb_AB-1-j][meta->nrb_AB-1-j].val != NULL) {
      update_upper_row_block(mat[meta->nrb_AB-1-j], meta->nrb_AB-1-j, meta, t);
    }
    reduce_upper_column_block(mat, meta->nrb_AB-1-j, meta, t);
  }
}
#endif
#endif
