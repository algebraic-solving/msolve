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

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#ifdef HAVE_AVX2
#include <immintrin.h>

#define AVX2LOAD(A) _mm256_load_si256((__m256i*)(A))
#define AVX2LOADU(A) _mm256_loadu_si256((__m256i*)(A))
#define AVX2STORE(res,A) _mm256_store_si256((__m256i*)(res),A);
#define AVX2STOREU(res,A) _mm256_storeu_si256((__m256i*)(res),A);
#define AVX2SETZERO() _mm256_setzero_si256()
#define AVX2SET1_64(A0) _mm256_set1_epi64x(A0)
#define AVX2AND_(A,B) _mm256_and_si256(A,B)
#define AVX2ADD_64(A,B) _mm256_add_epi64(A,B)
#define AVX2MUL(A,B) _mm256_mul_epu32(A,B)
#define AVX2SRLI_64(A,i) _mm256_srli_epi64(A,i)
#endif

/* Matrix-vector product in GF(p), with p on 32 bits. */

/* Euclidean division of LENGTH by 8 */
/* LENGTH/8 */
#define LENGTHQ(l) (l>>3)
/* LENGTH%8 */
#define LENGTHR(l) (l&7)

#define LENGTHQ8(l) (l / 32)
#define LENGTHR8(l) (l%32)

#define LENGTHQ14(l) (l / 56)
#define LENGTHR14(l) (l%56)

//-1 sur 32 bits
#define MONE32 ((uint32_t)0xFFFFFFFF)

#define ALIGNED32_MALLOC(p,type,nmemb,size) p=(type)_mm_malloc(nmemb*size,32);
#define ALIGNED32_CALLOC(p,type,nmemb,size)     \
  ALIGNED32_MALLOC(p,type,nmemb,size);          \
  memset((void*)p,0,nmemb*size);
#define ALIGNED_FREE(p) _mm_free(p);
#define ALIGNED32 __attribute__((aligned(32)))

/* Faster when set to 1 */
//#define IS_ALIGNED_TMP 1
/* Requirement: the vector is aligned on 32 bytes */
//#define IS_ALIGNED_VECTOR 1
/* Requirement: the matrix is aligned on 32 bytes */
//#define IS_ALIGNED_MATRIX 1

#include <inttypes.h>
#define PRINT_U32(a) printf("%"PRIu32,a);


static inline uint32_t mul_mod_barrett(const uint32_t x, const uint32_t y,
                                       const uint32_t p, const uint32_t pi){
  const int n = 32;
  uint32_t c = (x * ((int64_t) pi)) >> n;
  int64_t d = ((int64_t) x) * y - ((int64_t) c) * p;
  return (d >= p) ? d - p : d;
}

static inline void non_avx_matrix_vector_product(uint32_t* vec_res, const uint32_t* mat,
                                         const uint32_t* vec, const uint32_t ncols,
                                         const uint32_t nrows, const uint32_t PRIME,
					 const uint32_t RED_32, const uint32_t RED_64,
					 md_t *st)
{
    uint32_t i, j;
    int64_t prod1, prod2, prod3, prod4;
    const int64_t modsquare = (int64_t)PRIME*PRIME;

    j = 0;
    if (nrows >= 4) {
        for (j = 0; j < nrows-3; j += 4) {
            i = 0;
            prod1 =  0;
            prod2 =  0;
            prod3 =  0;
            prod4 =  0;
            if (ncols >= 8) {
                while (i < ncols-7) {
                    prod1 -=  (int64_t)mat[j*ncols+i] * vec[i];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i] * vec[i];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i] * vec[i];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i] * vec[i];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+1] * vec[i+1];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+1] * vec[i+1];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+1] * vec[i+1];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+1] * vec[i+1];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+2] * vec[i+2];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+2] * vec[i+2];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+2] * vec[i+2];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+2] * vec[i+2];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+3] * vec[i+3];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+3] * vec[i+3];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+3] * vec[i+3];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+3] * vec[i+3];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+4] * vec[i+4];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+4] * vec[i+4];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+4] * vec[i+4];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+4] * vec[i+4];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+5] * vec[i+5];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+5] * vec[i+5];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+5] * vec[i+5];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+5] * vec[i+5];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+6] * vec[i+6];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+6] * vec[i+6];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+6] * vec[i+6];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+6] * vec[i+6];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    prod1 -=  (int64_t)mat[j*ncols+i+7] * vec[i+7];
                    prod2 -=  (int64_t)mat[(j+1)*ncols+i+7] * vec[i+7];
                    prod3 -=  (int64_t)mat[(j+2)*ncols+i+7] * vec[i+7];
                    prod4 -=  (int64_t)mat[(j+3)*ncols+i+7] * vec[i+7];
                    prod1 +=  ((prod1 >> 63)) & modsquare;
                    prod2 +=  ((prod2 >> 63)) & modsquare;
                    prod3 +=  ((prod3 >> 63)) & modsquare;
                    prod4 +=  ((prod4 >> 63)) & modsquare;
                    i     +=  8;
                }
            }
            while (i < ncols) {
                prod1 -=  (int64_t)mat[j*ncols+i] * vec[i];
                prod2 -=  (int64_t)mat[(j+1)*ncols+i] * vec[i];
                prod3 -=  (int64_t)mat[(j+2)*ncols+i] * vec[i];
                prod4 -=  (int64_t)mat[(j+3)*ncols+i] * vec[i];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod2 +=  ((prod2 >> 63)) & modsquare;
                prod3 +=  ((prod3 >> 63)) & modsquare;
                prod4 +=  ((prod4 >> 63)) & modsquare;
                i     +=  1;
            }
            /* ensure prod being positive */
            prod1 =   -prod1;
            prod1 +=  (prod1 >> 63) & modsquare;
            prod2 =   -prod2;
            prod2 +=  (prod2 >> 63) & modsquare;
            prod3 =   -prod3;
            prod3 +=  (prod3 >> 63) & modsquare;
            prod4 =   -prod4;
            prod4 +=  (prod4 >> 63) & modsquare;
            vec_res[j]    = (uint32_t)(prod1 % PRIME);
            vec_res[j+1]  = (uint32_t)(prod2 % PRIME);
            vec_res[j+2]  = (uint32_t)(prod3 % PRIME);
            vec_res[j+3]  = (uint32_t)(prod4 % PRIME);
        }
    }
    for (; j < nrows; ++j) {
        i = 0;
        prod1 =  0;
        if (ncols >= 8) {
            while (i < ncols-7) {
                prod1 -=  (int64_t)mat[j*ncols+i] * vec[i];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+1] * vec[i+1];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+2] * vec[i+2];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+3] * vec[i+3];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+4] * vec[i+4];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+5] * vec[i+5];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+6] * vec[i+6];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 -=  (int64_t)mat[j*ncols+i+7] * vec[i+7];
                prod1 +=  ((prod1 >> 63)) & modsquare;
                prod1 +=  ((prod1 >> 63)) & modsquare;
                i     +=  8;
            }
        }
        while (i < ncols) {
            prod1 -=  (int64_t)mat[j*ncols+i] * vec[i];
            prod1 +=  ((prod1 >> 63)) & modsquare;
            i     +=  1;
        }
        /* ensure prod being positive */
        prod1 =   -prod1;
        prod1 +=  (prod1 >> 63) & modsquare;
        vec_res[j]    = (uint32_t)(prod1 % PRIME);
    }
}

#ifdef HAVE_AVX2
static inline void matrix_vector_product(uint32_t* vec_res, const uint32_t* mat,
                                         const uint32_t* vec, const uint32_t ncols,
                                         const uint32_t nrows, const uint32_t PRIME,
                                         const uint32_t RED_32, const uint32_t RED_64)
{
    __m256i acc_low,acc_high,prod,mask,vec8,mat8;
#if IS_ALIGNED_TMP
    uint64_t acc4x64[8] ALIGNED32;
    //#if LENGTHR(ncols)
    uint32_t vecr[8] ALIGNED32,matr[8] ALIGNED32;
    //#endif
#else
    uint64_t acc4x64[8];
    //#if LENGTHR(ncols)
    uint32_t vecr[8],matr[8];
    //#endif
#endif
    uint64_t acc64;
    const uint32_t *vec_cp;

    unsigned int i,j;

    //mask pour recuperer les parties basses
    mask=AVX2SET1_64(MONE32);

    if(LENGTHR(ncols)){
      for(i=LENGTHR(ncols);i<8;++i)
        {
          matr[i]=0;
          vecr[i]=0;
        }
    }

    /* For each row of the matrix, we compute a dot product */
    for(j=0;j<nrows;++j){
      vec_cp=vec;

      acc_low=AVX2SETZERO();
      acc_high=AVX2SETZERO();

      /* Accumulation */
      for(i=0;i<LENGTHQ(ncols);++i)
        {
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat);
#else
          mat8=AVX2LOADU(mat);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          mat+=8;
          vec_cp+=8;
          /* four 32-bits mul, lower parts */
          prod=AVX2MUL(mat8,vec8);
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          /* four 32-bits mul, higher parts */
          prod=AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
        }

      /* Remainder */
      if(LENGTHR(ncols)){
        for(i=0;i<LENGTHR(ncols);++i)
          {
            matr[i]=mat[i];
            vecr[i]=vec_cp[i];
          }
        mat+=LENGTHR(ncols);

#if IS_ALIGNED_TMP
          mat8=AVX2LOAD(matr);
          vec8=AVX2LOAD(vecr);
#else
          mat8=AVX2LOADU(matr);
          vec8=AVX2LOADU(vecr);
#endif

          /* four 32-bits mul, lower parts */
          prod=AVX2MUL(mat8,vec8);
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          /* Useless when LENGTHR==1, because higher parts are 0 */
          if(LENGTHR(ncols)>1){
            /* four 32-bits mul, higher parts */
            prod=AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
            acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
          }
        }


        #if IS_ALIGNED_TMP
            AVX2STORE(acc4x64,acc_low);
            AVX2STORE(acc4x64+4,acc_high);
        #else
            AVX2STOREU(acc4x64,acc_low);
            AVX2STOREU(acc4x64+4,acc_high);
        #endif

        /* Reduction */
        acc64=0;
        for(i=0;i<4;++i){
          //partie haute du registre haut (2^64 ->2^95)
          acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
          //partie basse du registre haut (2^32->2^63)
          acc4x64[i]+=((acc4x64[i+4]&((uint64_t)0xFFFFFFFF))*RED_32)%PRIME;
          //          acc64+=acc4x64[i]%PRIME;
          acc64+=acc4x64[i]; 
        }
        *vec_res=acc64%PRIME;
        ++vec_res;
    }
}


static inline void _2mul_new_matrix_vector_product(uint32_t* vec_res, const uint32_t* mat,
                                         const uint32_t* vec, const uint32_t ncols,
                                         const uint32_t nrows, const uint32_t PRIME,
                                         const uint32_t RED_32, const uint32_t RED_64)
{
  __m256i acc_low,acc_high,mask,vec8,mat8,prod,prod1,prod2;
#if IS_ALIGNED_TMP
    uint64_t acc4x64[8] ALIGNED32;
    //#if LENGTHR(ncols)
    uint32_t vecr[8] ALIGNED32,matr[8] ALIGNED32;
    //#endif
#else
    uint64_t acc4x64[8];
    //#if LENGTHR(ncols)
    uint32_t vecr[8],matr[8];
    //#endif
#endif
    uint64_t acc64;
    const uint32_t *vec_cp;

    unsigned int i,j;

    //mask pour recuperer les parties basses
    mask=AVX2SET1_64(MONE32);

    if(LENGTHR(ncols)){
      for(i=LENGTHR(ncols);i<8;++i)
        {
          matr[i]=0;
          vecr[i]=0;
        }
    }

    /* For each row of the matrix, we compute a dot product */
    for(j=0;j<nrows;++j){
      vec_cp=vec;

      acc_low=AVX2SETZERO();
      acc_high=AVX2SETZERO();

      /* Accumulation */
      for(i=0;i<LENGTHQ(ncols);++i)
        {
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat);
#else
          mat8=AVX2LOADU(mat);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          mat+=8;
          vec_cp+=8;
          /* four 32-bits mul, lower parts */
          /* four 32-bits mul, higher parts */
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          prod = prod1+prod2;
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          //          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          //          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
        }

      /* Remainder */
      if(LENGTHR(ncols)){
        for(i=0;i<LENGTHR(ncols);++i)
          {
            matr[i]=mat[i];
            vecr[i]=vec_cp[i];
          }
        mat+=LENGTHR(ncols);

#if IS_ALIGNED_TMP
          mat8=AVX2LOAD(matr);
          vec8=AVX2LOAD(vecr);
#else
          mat8=AVX2LOADU(matr);
          vec8=AVX2LOADU(vecr);
#endif

          /* four 32-bits mul, lower parts */
          prod=AVX2MUL(mat8,vec8);
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          /* Useless when LENGTHR==1, because higher parts are 0 */
          if(LENGTHR(ncols)>1){
            /* four 32-bits mul, higher parts */
            prod=AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
            acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
          }
        }


        #if IS_ALIGNED_TMP
            AVX2STORE(acc4x64,acc_low);
            AVX2STORE(acc4x64+4,acc_high);
        #else
            AVX2STOREU(acc4x64,acc_low);
            AVX2STOREU(acc4x64+4,acc_high);
        #endif

        /* Reduction */
        acc64=0;
        for(i=0;i<4;++i){
          //partie haute du registre haut (2^64 ->2^95)
          acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
          //partie basse du registre haut (2^32->2^63)
          acc4x64[i]+=((acc4x64[i+4]&((uint64_t)0xFFFFFFFF))*RED_32)%PRIME;
          //          acc64+=acc4x64[i]%PRIME;
          acc64+=acc4x64[i]; 
        }
        *vec_res=acc64%PRIME;
        ++vec_res;
    }
}

static inline void _4mul_new_matrix_vector_product(uint32_t* vec_res, const uint32_t* mat,
                                         const uint32_t* vec, const uint32_t ncols,
                                         const uint32_t nrows, const uint32_t PRIME,
                                         const uint32_t RED_32, const uint32_t RED_64)
{
  __m256i acc_low,acc_high,mask,vec8,mat8;
  __m256i prod,prod1,prod2,prod3,prod4;
#if IS_ALIGNED_TMP
    uint64_t acc4x64[8] ALIGNED32;
    //#if LENGTHR(ncols)
    uint32_t vecr[8] ALIGNED32,matr[8] ALIGNED32;
    //#endif
#else
    uint64_t acc4x64[8];
    //#if LENGTHR(ncols)
    uint32_t vecr[8],matr[8];
    //#endif
#endif
    uint64_t acc64;
    const uint32_t *vec_cp;

    unsigned int i,j;

    //mask pour recuperer les parties basses
    mask=AVX2SET1_64(MONE32);

    if(LENGTHR(ncols)){
      for(i=LENGTHR(ncols);i<8;++i)
        {
          matr[i]=0;
          vecr[i]=0;
        }
    }

    /* For each row of the matrix, we compute a dot product */
    for(j=0;j<nrows;++j){
      vec_cp=vec;

      acc_low=AVX2SETZERO();
      acc_high=AVX2SETZERO();

      /* Accumulation */
      for(i=0;i<LENGTHQ(ncols);++i)
        {
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat);
#else
          mat8=AVX2LOADU(mat);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          mat+=8;
          vec_cp+=8;
          /* four 32-bits mul, lower parts */
          /* four 32-bits mul, higher parts */
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));

          i++;
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat);
#else
          mat8=AVX2LOADU(mat);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif

          mat+=8;
          vec_cp+=8;
          /* four 32-bits mul, lower parts */
          /* four 32-bits mul, higher parts */
          prod3=AVX2MUL(mat8,vec8);
          prod4= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));

          prod = prod1+prod2+prod3+prod4;
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          //          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          //          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
        }

      /* Remainder */
      if(LENGTHR(ncols)){
        for(i=0;i<LENGTHR(ncols);++i)
          {
            matr[i]=mat[i];
            vecr[i]=vec_cp[i];
          }
        mat+=LENGTHR(ncols);

#if IS_ALIGNED_TMP
          mat8=AVX2LOAD(matr);
          vec8=AVX2LOAD(vecr);
#else
          mat8=AVX2LOADU(matr);
          vec8=AVX2LOADU(vecr);
#endif

          /* four 32-bits mul, lower parts */
          prod=AVX2MUL(mat8,vec8);
          acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));

          /* Useless when LENGTHR==1, because higher parts are 0 */
          if(LENGTHR(ncols)>1){
            /* four 32-bits mul, higher parts */
            prod=AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
            acc_low=AVX2ADD_64(acc_low,AVX2AND_(prod,mask));
            acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(prod,32));
          }
        }


        #if IS_ALIGNED_TMP
            AVX2STORE(acc4x64,acc_low);
            AVX2STORE(acc4x64+4,acc_high);
        #else
            AVX2STOREU(acc4x64,acc_low);
            AVX2STOREU(acc4x64+4,acc_high);
        #endif

        /* Reduction */
        acc64=0;
        for(i=0;i<4;++i){
          //partie haute du registre haut (2^64 ->2^95)
          acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
          //partie basse du registre haut (2^32->2^63)
          acc4x64[i]+=((acc4x64[i+4]&((uint64_t)0xFFFFFFFF))*RED_32)%PRIME;
          //          acc64+=acc4x64[i]%PRIME;
          acc64+=acc4x64[i]; 
        }
        *vec_res=acc64%PRIME;
        ++vec_res;
    }
}


/**
AVX2-based matrix vector product
**/

static inline void _8mul_matrix_vector_product(uint32_t* vec_res,
                                               const uint32_t* mat,
                                               const uint32_t* vec,
                                               const uint32_t *dst,
                                               const uint32_t ncols,
                                               const uint32_t nrows,
                                               const uint32_t PRIME,
                                               const uint32_t RED_32,
                                               const uint32_t RED_64,
                                               const uint32_t preinv,
					       md_t *st){
    //mask pour recuperer les parties basses
    __m256i mask=AVX2SET1_64(MONE32);
    const long quo0 = LENGTHQ8(ncols);
    const long rem0 = LENGTHR8(ncols);

#pragma omp parallel num_threads (st->nthrds)
    {
      unsigned int i,j;
#if IS_ALIGNED_TMP
      uint64_t acc4x64[8] ALIGNED32;
#else
      uint64_t acc4x64[8];
#endif
      uint64_t acc64;
      __m256i acc_low,acc_high,vec8,mat8;
      __m256i prod1,prod2;
      __m256i res1;
      const uint32_t *vec_cp;
      const uint32_t *mat_cp;
    
      /* parallelization of the outer loop (rows) */
#pragma omp for
      /* For each row of the matrix, we compute a dot product */
      for(j = 0; j < nrows; ++j){
	vec_cp=vec;
	mat_cp=mat + j*ncols;
	
	acc_low=AVX2SETZERO();
	acc_high=AVX2SETZERO();
	
	long local_ncols, quo, rem;
	if(dst[j] != 0){
	  local_ncols = ncols - dst[j];
	  quo = LENGTHQ8(local_ncols);
	  rem = LENGTHR8(local_ncols);
	}
	else{
	  local_ncols = ncols;
	  quo = quo0;
	  rem = rem0;
	}
	const long bs =  quo;
	/* Accumulation */
	for(i = 0 ; i < bs; ++i){
	  
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat_cp);
#else
          mat8=AVX2LOADU(mat_cp);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          /* four 32-bits mul, lower parts */
          /* four 32-bits mul, higher parts */
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          res1 = prod1 + prod2;
          mat_cp+=8;
          vec_cp+=8;
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat_cp);
#else
          mat8=AVX2LOADU(mat_cp);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          res1 += prod1 + prod2;
          mat_cp+=8;
          vec_cp+=8;
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat_cp);
#else
          mat8=AVX2LOADU(mat_cp);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          res1 += prod1 + prod2;
          mat_cp+=8;
          vec_cp+=8;
#if IS_ALIGNED_MATRIX
          mat8=AVX2LOAD(mat_cp);
#else
          mat8=AVX2LOADU(mat_cp);
#endif
#if IS_ALIGNED_VECTOR
          vec8=AVX2LOAD(vec_cp);
#else
          vec8=AVX2LOADU(vec_cp);
#endif
          prod1=AVX2MUL(mat8,vec8);
          prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));
          res1 += prod1 + prod2;
          mat_cp+=8;
          vec_cp+=8;

          acc_low=AVX2ADD_64(acc_low,AVX2AND_(res1,mask));
          acc_high=AVX2ADD_64(acc_high,AVX2SRLI_64(res1,32));

	}
	
#if IS_ALIGNED_TMP
	AVX2STORE(acc4x64,acc_low);
	AVX2STORE(acc4x64+4,acc_high);
#else
	AVX2STOREU(acc4x64,acc_low);
	AVX2STOREU(acc4x64+4,acc_high);
#endif
	
	/* Reduction */
	acc64=0;
	for(i=0;i<4;++i){
	  //partie haute du registre haut (2^64 ->2^95)
	  acc4x64[i]+=((acc4x64[i+4]>>32)*RED_64)%PRIME;
	  //partie basse du registre haut (2^32->2^63)
	  acc4x64[i]+=((acc4x64[i+4]&((uint64_t)0xFFFFFFFF))*RED_32)%PRIME;
	  acc64+=acc4x64[i]%PRIME;
	}
	
	/* *vec_res[j]=MODRED32(acc64, PRIME, preinv); */
	vec_res[j]=acc64%PRIME;
	
	long tmp = 0;
	for(long k = 0; k < rem; k++){
	  tmp += ((long)mat_cp[k] * (long)vec_cp[k]) % PRIME;
	}
	vec_res[j] = (vec_res[j] + tmp) % PRIME;
      }
    }
}
#endif
