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

#include <stdint.h>
#include <inttypes.h>

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

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
