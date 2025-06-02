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
#include <flint/nmod.h>

#ifdef HAVE_AVX2
#include <immintrin.h>
#endif

// for dot product, simultaneous dot products, and matrix-vector products,
// vectorized functions below accumulate 8 terms:
//  -> modulus <= DOT2_ACC8_MAX_MODULUS   (about 2**30.5, slightly less)
//  -> len <= DOT2_ACC8_MAX_LEN           (about 2**28.5)
// (see bottom of flint/src/nmod_vec/dot.c for more details)
#define DOT2_ACC8_MAX_MODULUS UWORD(1515531528)
#define DOT2_ACC8_MAX_LEN UWORD(380368697)

// parameters for splitting
#define DOT_SPLIT_BITS 56
#define DOT_SPLIT_MASK 72057594037927935UL // (1UL << DOT_SPLIT_BITS) - 1

/*--------------------------------------*/
/* non-vectorized matrix vector product */
/*--------------------------------------*/
// FIXME what are the constraints (prime bitsize, length, ..) for this non-vectorized version?

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

/*-----------------------------------------*/
/* vectorized (AVX2) matrix vector product */
/*-----------------------------------------*/

#ifdef HAVE_AVX2

// avx2 horizontal sum
static inline
uint64_t _mm256_hsum(__m256i a)
{
    __m256i a_hi = _mm256_shuffle_epi32(a, 14);  // 14 == 0b00001110
    __m256i sum_lo = _mm256_add_epi64(a, a_hi);
    __m128i sum_hi = _mm256_extracti128_si256(sum_lo, 1);
    __m128i sum = _mm_add_epi64(_mm256_castsi256_si128(sum_lo), sum_hi);
    return (uint64_t) _mm_cvtsi128_si64(sum);
}

uint32_t _nmod32_vec_dot_split_avx2(const uint32_t * vec1, const uint32_t * vec2, int64_t len,
                                    nmod_t mod, uint64_t pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();

    int64_t i = 0;

    for ( ; i+31 < len; i+=32)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 0));
        __m256i v1_1 = _mm256_loadu_si256((const __m256i *) (vec1+i+ 8));
        __m256i v1_2 = _mm256_loadu_si256((const __m256i *) (vec1+i+16));
        __m256i v1_3 = _mm256_loadu_si256((const __m256i *) (vec1+i+24));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 0));
        __m256i v2_1 = _mm256_loadu_si256((const __m256i *) (vec2+i+ 8));
        __m256i v2_2 = _mm256_loadu_si256((const __m256i *) (vec2+i+16));
        __m256i v2_3 = _mm256_loadu_si256((const __m256i *) (vec2+i+24));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        __m256i dp_lo1 = _mm256_mul_epu32(v1_1, v2_1);
        __m256i dp_lo2 = _mm256_mul_epu32(v1_2, v2_2);
        __m256i dp_lo3 = _mm256_mul_epu32(v1_3, v2_3);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v1_2 = _mm256_shuffle_epi32(v1_2, 0xB1);
        v1_3 = _mm256_shuffle_epi32(v1_3, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm256_shuffle_epi32(v2_1, 0xB1);
        v2_2 = _mm256_shuffle_epi32(v2_2, 0xB1);
        v2_3 = _mm256_shuffle_epi32(v2_3, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_1, v2_1));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_2, v2_2));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_3, v2_3));

        // gather results in dp_lo0
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo1);
        dp_lo2 = _mm256_add_epi64(dp_lo2, dp_lo3);
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo2);

        // split
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    }

    // the following loop iterates < 4 times,
    // each iteration accumulates 2 terms in dp_lo0
    for ( ; i+7 < len; i+=8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0 = _mm256_loadu_si256((const __m256i *) (vec2+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm256_shuffle_epi32(v2_0, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0));
    }

    dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);

    uint64_t hsum_lo = _mm256_hsum(dp_lo0);
    uint64_t hsum_hi = _mm256_hsum(dp_hi0) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
        hsum_lo += (uint64_t)vec1[i] * vec2[i];

    hsum_hi += (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // the requirement on "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    uint64_t res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return (uint32_t)res;
}

void _nmod32_vec_dot2_split_avx2(uint32_t * res0, uint32_t * res1,
                                 const uint32_t * vec1, const uint32_t * vec2_0, const uint32_t * vec2_1,
                                 int64_t len, nmod_t mod, uint64_t pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);
    __m256i dp_lo0 = _mm256_setzero_si256();
    __m256i dp_lo1 = _mm256_setzero_si256();
    __m256i dp_hi0 = _mm256_setzero_si256();
    __m256i dp_hi1 = _mm256_setzero_si256();

    int64_t i = 0;

    for ( ; i+31 < len; i+=32)
    {
        __m256i v1_0   = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        __m256i v1_1   = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 8));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        __m256i v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        __m256i v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        __m256i dp_lo2 = _mm256_mul_epu32(v1_1, v2_0_1);
        __m256i dp_lo3 = _mm256_mul_epu32(v1_1, v2_1_1);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        v1_0   = _mm256_loadu_si256((const __m256i *) (vec1  +i+16));
        v1_1   = _mm256_loadu_si256((const __m256i *) (vec1  +i+24));
        v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+16));
        v2_0_1 = _mm256_loadu_si256((const __m256i *) (vec2_0+i+24));
        v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+16));
        v2_1_1 = _mm256_loadu_si256((const __m256i *) (vec2_1+i+24));

        // 3rd term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        // 4th term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm256_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm256_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm256_shuffle_epi32(v2_1_1, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm256_add_epi64(dp_lo2, _mm256_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm256_add_epi64(dp_lo3, _mm256_mul_epu32(v1_1, v2_1_1));

        // gather results in dp_lo0 and dp_lo1
        dp_lo0 = _mm256_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm256_add_epi64(dp_lo1, dp_lo3);

        // split
        dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
        dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0 and dp_lo1
    for ( ; i+7 < len; i+=8)
    {
        __m256i v1_0 = _mm256_loadu_si256((const __m256i *) (vec1+i));
        __m256i v2_0_0 = _mm256_loadu_si256((const __m256i *) (vec2_0+i));
        __m256i v2_1_0 = _mm256_loadu_si256((const __m256i *) (vec2_1+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm256_shuffle_epi32(v1_0, 0xB1);
        v2_0_0 = _mm256_shuffle_epi32(v2_0_0, 0xB1);
        v2_1_0 = _mm256_shuffle_epi32(v2_1_0, 0xB1);

        dp_lo0 = _mm256_add_epi64(dp_lo0, _mm256_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm256_add_epi64(dp_lo1, _mm256_mul_epu32(v1_0, v2_1_0));
    }

    dp_hi0 = _mm256_add_epi64(dp_hi0, _mm256_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm256_and_si256(dp_lo0, low_bits);
    dp_hi1 = _mm256_add_epi64(dp_hi1, _mm256_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
    dp_lo1 = _mm256_and_si256(dp_lo1, low_bits);

    uint64_t hsum_lo0 = _mm256_hsum(dp_lo0);
    uint64_t hsum_hi0 = _mm256_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    uint64_t hsum_lo1 = _mm256_hsum(dp_lo1);
    uint64_t hsum_hi1 = _mm256_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
    {
        hsum_lo0 += (uint64_t)vec1[i] * vec2_0[i];
        hsum_lo1 += (uint64_t)vec1[i] * vec2_1[i];
    }
    hsum_hi0 += (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    hsum_hi1 += (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    // the requirement on "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

void _nmod32_vec_dot3_split_avx2(uint32_t * res, const uint32_t * vec1,
                                 const uint32_t * vec2_0, const uint32_t * vec2_1, const uint32_t * vec2_2,
                                 int64_t len, nmod_t mod, uint64_t pow2_precomp)
{
    const __m256i low_bits = _mm256_set1_epi64x(DOT_SPLIT_MASK);

    __m256i dp_lo[3];
    __m256i dp_hi[3];
    __m256i v1;
    __m256i v2[3];

    dp_lo[0] = _mm256_setzero_si256();
    dp_lo[1] = _mm256_setzero_si256();
    dp_lo[2] = _mm256_setzero_si256();

    dp_hi[0] = _mm256_setzero_si256();
    dp_hi[1] = _mm256_setzero_si256();
    dp_hi[2] = _mm256_setzero_si256();

    int64_t i = 0;

    for ( ; i+31 < len; i+=32)
    {
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 0));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 3rd+4th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 8));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 8));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 8));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 8));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 5th+6th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+16));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+16));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+16));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+16));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 7th+8th term
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+24));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+24));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+24));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+24));

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // split
        dp_hi[0] = _mm256_add_epi64(dp_hi[0], _mm256_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
        dp_hi[1] = _mm256_add_epi64(dp_hi[1], _mm256_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
        dp_hi[2] = _mm256_add_epi64(dp_hi[2], _mm256_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
        dp_lo[0] = _mm256_and_si256(dp_lo[0], low_bits);
        dp_lo[1] = _mm256_and_si256(dp_lo[1], low_bits);
        dp_lo[2] = _mm256_and_si256(dp_lo[2], low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms
    for ( ; i+7 < len; i += 8)
    {
        v1    = _mm256_loadu_si256((const __m256i *) (vec1  +i+ 0));
        v2[0] = _mm256_loadu_si256((const __m256i *) (vec2_0+i+ 0));
        v2[1] = _mm256_loadu_si256((const __m256i *) (vec2_1+i+ 0));
        v2[2] = _mm256_loadu_si256((const __m256i *) (vec2_2+i+ 0));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        v1    = _mm256_shuffle_epi32(v1   , 0xB1);
        v2[0] = _mm256_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm256_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm256_shuffle_epi32(v2[2], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm256_add_epi64(dp_lo[0], _mm256_mul_epu32(v1, v2[0]));
        dp_lo[1] = _mm256_add_epi64(dp_lo[1], _mm256_mul_epu32(v1, v2[1]));
        dp_lo[2] = _mm256_add_epi64(dp_lo[2], _mm256_mul_epu32(v1, v2[2]));
    }

    // split
    dp_hi[0] = _mm256_add_epi64(dp_hi[0], _mm256_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
    dp_hi[1] = _mm256_add_epi64(dp_hi[1], _mm256_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
    dp_hi[2] = _mm256_add_epi64(dp_hi[2], _mm256_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
    dp_lo[0] = _mm256_and_si256(dp_lo[0], low_bits);
    dp_lo[1] = _mm256_and_si256(dp_lo[1], low_bits);
    dp_lo[2] = _mm256_and_si256(dp_lo[2], low_bits);

    uint64_t hsum_lo[3];
    uint64_t hsum_hi[3];
    hsum_lo[0] = _mm256_hsum(dp_lo[0]);
    hsum_hi[0] = _mm256_hsum(dp_hi[0]) + (hsum_lo[0] >> DOT_SPLIT_BITS);
    hsum_lo[0] &= DOT_SPLIT_MASK;
    hsum_lo[1] = _mm256_hsum(dp_lo[1]);
    hsum_hi[1] = _mm256_hsum(dp_hi[1]) + (hsum_lo[1] >> DOT_SPLIT_BITS);
    hsum_lo[1] &= DOT_SPLIT_MASK;
    hsum_lo[2] = _mm256_hsum(dp_lo[2]);
    hsum_hi[2] = _mm256_hsum(dp_hi[2]) + (hsum_lo[2] >> DOT_SPLIT_BITS);
    hsum_lo[2] &= DOT_SPLIT_MASK;

    // less than 8 terms remaining, can accumulate
    for (; i < len; i++)
    {
        hsum_lo[0] += (uint64_t)vec1[i] * vec2_0[i];
        hsum_lo[1] += (uint64_t)vec1[i] * vec2_1[i];
        hsum_lo[2] += (uint64_t)vec1[i] * vec2_2[i];
    }

    hsum_hi[0] += (hsum_lo[0] >> DOT_SPLIT_BITS);
    hsum_hi[1] += (hsum_lo[1] >> DOT_SPLIT_BITS);
    hsum_hi[2] += (hsum_lo[2] >> DOT_SPLIT_BITS);
    hsum_lo[0] &= DOT_SPLIT_MASK;
    hsum_lo[1] &= DOT_SPLIT_MASK;
    hsum_lo[2] &= DOT_SPLIT_MASK;

    NMOD_RED(res[0], pow2_precomp * hsum_hi[0] + hsum_lo[0], mod);
    NMOD_RED(res[1], pow2_precomp * hsum_hi[1] + hsum_lo[1], mod);
    NMOD_RED(res[2], pow2_precomp * hsum_hi[2] + hsum_lo[2], mod);
}

static inline void _avx2_matrix_vector_product(uint32_t * vec_res,
                                               const uint32_t * mat,
                                               const uint32_t * vec,
                                               const uint32_t * dst,
                                               const uint32_t ncols,
                                               const uint32_t nrows,
                                               const nmod_t mod,
                                               const uint64_t pow2_precomp,
                                               md_t *st)
{
    slong i = 0;

    for ( ; i+2 < nrows; i+=3)
    {
        int64_t len = ncols - MIN(dst[i], MIN(dst[i+1], dst[i+2]));
        _nmod32_vec_dot3_split_avx2(vec_res+i, vec,
                                    mat + i*ncols,
                                    mat + (i+1)*ncols,
                                    mat + (i+2)*ncols,
                                    len, mod, pow2_precomp);
    }

    if (nrows - i == 2)
    {
        int64_t len = ncols - MIN(dst[i], dst[i+1]);
        _nmod32_vec_dot2_split_avx2(vec_res+i, vec_res+i+1,
                                    vec, mat + i*ncols, mat + (i+1)*ncols,
                                    len, mod, pow2_precomp);
    }
    else if (nrows - i == 1)
        vec_res[i] = _nmod32_vec_dot_split_avx2(vec, mat + i*ncols, ncols - dst[i], mod, pow2_precomp);
}
#endif

/*-------------------------------------------*/
/* vectorized (AVX512) matrix vector product */
/*-------------------------------------------*/

#ifdef HAVE_AVX512_F

// avx512 horizontal sum
FLINT_FORCE_INLINE uint64_t _mm512_hsum(__m512i a)
{
    return _mm512_reduce_add_epi64(a);
}

uint32_t _nmod32_vec_dot_split_avx512(const uint32_t * vec1, const uint32_t * vec2, int64_t len, nmod_t mod, uint64_t pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();

    int64_t i = 0;

    for ( ; i+63 < len; i+=64)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i+ 0));
        __m512i v1_1 = _mm512_loadu_si512((const __m512i *) (vec1+i+16));
        __m512i v1_2 = _mm512_loadu_si512((const __m512i *) (vec1+i+32));
        __m512i v1_3 = _mm512_loadu_si512((const __m512i *) (vec1+i+48));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2+i+ 0));
        __m512i v2_1 = _mm512_loadu_si512((const __m512i *) (vec2+i+16));
        __m512i v2_2 = _mm512_loadu_si512((const __m512i *) (vec2+i+32));
        __m512i v2_3 = _mm512_loadu_si512((const __m512i *) (vec2+i+48));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        __m512i dp_lo1 = _mm512_mul_epu32(v1_1, v2_1);
        __m512i dp_lo2 = _mm512_mul_epu32(v1_2, v2_2);
        __m512i dp_lo3 = _mm512_mul_epu32(v1_3, v2_3);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v1_2 = _mm512_shuffle_epi32(v1_2, 0xB1);
        v1_3 = _mm512_shuffle_epi32(v1_3, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);
        v2_2 = _mm512_shuffle_epi32(v2_2, 0xB1);
        v2_3 = _mm512_shuffle_epi32(v2_3, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm512_srli_epi64(v1_0, 32)

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_1, v2_1));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_2, v2_2));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_3, v2_3));

        // gather results in dp_lo0
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo1);
        dp_lo2 = _mm512_add_epi64(dp_lo2, dp_lo3);
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo2);

        // split
        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0
    for ( ; i+15 < len; i+=16)
    {
        __m512i v1_0 = _mm512_loadu_si512((const __m512i *) (vec1+i));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        // 2nd term
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        __m512i v1_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1+i));
        __m512i v2_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
        // 2nd term
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0));
    }

    // split
    dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);

    // gather 8 terms in single uint64_t
    uint64_t hsum_lo = _mm512_hsum(dp_lo0);
    uint64_t hsum_hi = _mm512_hsum(dp_hi0) + (hsum_lo >> DOT_SPLIT_BITS);
    hsum_lo &= DOT_SPLIT_MASK;

    // the requirement "len <= DOT2_ACC8_MAX_LEN"
    // ensures pow2_precomp * hsum_hi + hsum_lo fits in 64 bits
    uint64_t res;
    NMOD_RED(res, pow2_precomp * hsum_hi + hsum_lo, mod);
    return (uint32_t)res;
}

void _nmod32_vec_dot2_split_avx512(uint32_t * res0, uint32_t * res1,
                                   const uint32_t * vec1, const uint32_t * vec2_0, const uint32_t * vec2_1,
                                   int64_t len, nmod_t mod, uint64_t pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);
    __m512i dp_lo0 = _mm512_setzero_si512();
    __m512i dp_lo1 = _mm512_setzero_si512();
    __m512i dp_hi0 = _mm512_setzero_si512();
    __m512i dp_hi1 = _mm512_setzero_si512();

    int64_t i = 0;

    for ( ; i+63 < len; i+=64)
    {
        __m512i v1_0   = _mm512_loadu_si512((const __m512i *) (vec1  +i+ 0));
        __m512i v1_1   = _mm512_loadu_si512((const __m512i *) (vec1  +i+16));
        __m512i v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+ 0));
        __m512i v2_0_1 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+16));
        __m512i v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+ 0));
        __m512i v2_1_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+16));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        __m512i dp_lo2 = _mm512_mul_epu32(v1_1, v2_0_1);
        __m512i dp_lo3 = _mm512_mul_epu32(v1_1, v2_1_1);

        // 2nd term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm512_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm512_shuffle_epi32(v2_1_1, 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm256_srli_epi64(v1_0, 32)

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        v1_0   = _mm512_loadu_si512((const __m512i *) (vec1  +i+32));
        v1_1   = _mm512_loadu_si512((const __m512i *) (vec1  +i+48));
        v2_0_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+32));
        v2_0_1 = _mm512_loadu_si512((const __m512i *) (vec2_0+i+48));
        v2_1_0 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+32));
        v2_1_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i+48));

        // 3rd term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        // 4th term: high 32 bit word of each 64 bit word
        v1_0 = _mm512_shuffle_epi32(v1_0, 0xB1);
        v1_1 = _mm512_shuffle_epi32(v1_1, 0xB1);
        v2_0_0 = _mm512_shuffle_epi32(v2_0_0, 0xB1);
        v2_0_1 = _mm512_shuffle_epi32(v2_0_1, 0xB1);
        v2_1_0 = _mm512_shuffle_epi32(v2_1_0, 0xB1);
        v2_1_1 = _mm512_shuffle_epi32(v2_1_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1_0, v2_0_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1_0, v2_1_0));
        dp_lo2 = _mm512_add_epi64(dp_lo2, _mm512_mul_epu32(v1_1, v2_0_1));
        dp_lo3 = _mm512_add_epi64(dp_lo3, _mm512_mul_epu32(v1_1, v2_1_1));

        // gather results in dp_lo0 and dp_lo1
        dp_lo0 = _mm512_add_epi64(dp_lo0, dp_lo2);
        dp_lo1 = _mm512_add_epi64(dp_lo1, dp_lo3);

        // split
        dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
        dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
        dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
        dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo0 and dp_lo1
    for ( ; i+15 < len; i+=16)
    {
        __m512i v1   = _mm512_loadu_si512((const __m512i *) (vec1  +i));
        __m512i v2_0 = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        __m512i v2_1 = _mm512_loadu_si512((const __m512i *) (vec2_1+i));

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));

        // 2nd term: high 32 bit word of each 64 bit word
        v1 = _mm512_shuffle_epi32(v1, 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        __m512i v1   = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1  +i));
        __m512i v2_0 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_0+i));
        __m512i v2_1 = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_1+i));

        // 1st term
        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));

        // 2nd term
        v1   = _mm512_shuffle_epi32(v1  , 0xB1);
        v2_0 = _mm512_shuffle_epi32(v2_0, 0xB1);
        v2_1 = _mm512_shuffle_epi32(v2_1, 0xB1);

        dp_lo0 = _mm512_add_epi64(dp_lo0, _mm512_mul_epu32(v1, v2_0));
        dp_lo1 = _mm512_add_epi64(dp_lo1, _mm512_mul_epu32(v1, v2_1));
    }

    // split
    dp_hi0 = _mm512_add_epi64(dp_hi0, _mm512_srli_epi64(dp_lo0, DOT_SPLIT_BITS));
    dp_lo0 = _mm512_and_si512(dp_lo0, low_bits);
    dp_hi1 = _mm512_add_epi64(dp_hi1, _mm512_srli_epi64(dp_lo1, DOT_SPLIT_BITS));
    dp_lo1 = _mm512_and_si512(dp_lo1, low_bits);

    // gather 8 terms in single uint64_t
    uint64_t hsum_lo0 = _mm512_hsum(dp_lo0);
    uint64_t hsum_hi0 = _mm512_hsum(dp_hi0) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    uint64_t hsum_lo1 = _mm512_hsum(dp_lo1);
    uint64_t hsum_hi1 = _mm512_hsum(dp_hi1) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;

    NMOD_RED(*res0, pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(*res1, pow2_precomp * hsum_hi1 + hsum_lo1, mod);
}

void _nmod32_vec_dot3_split_avx512(uint32_t * res,
                                   const uint32_t * vec1, const uint32_t * vec2_0, const uint32_t * vec2_1, const uint32_t * vec2_2,
                                   int64_t len, nmod_t mod, uint64_t pow2_precomp)
{
    const __m512i low_bits = _mm512_set1_epi64(DOT_SPLIT_MASK);

    __m512i dp_lo[3];
    __m512i dp_hi[3];
    __m512i v1[2];
    __m512i v2[6];

    dp_lo[0] = _mm512_setzero_si512();
    dp_lo[1] = _mm512_setzero_si512();
    dp_lo[2] = _mm512_setzero_si512();

    dp_hi[0] = _mm512_setzero_si512();
    dp_hi[1] = _mm512_setzero_si512();
    dp_hi[2] = _mm512_setzero_si512();

    int64_t i = 0;

    for ( ; i+63 < len; i+=64)
    {
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i+ 0));
        v1[1] = _mm512_loadu_si512((const __m512i *) (vec1  +i+16));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+ 0));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+16));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+ 0));
        v2[3] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+16));
        v2[4] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+ 0));
        v2[5] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+16));

        // 1st+3rd term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // 2nd+4th term: high 32 bit word of each 64 bit word
        v1[0] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v1[1] = _mm512_shuffle_epi32(v1[1], 0xB1);
        v2[0] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm512_shuffle_epi32(v2[2], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[3], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[4], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[5], 0xB1);
        // the above uses vpshufd
        // shuffle [0,1,2,3] => [1,0,3,2]    (imm8 = 0b10110001  -->  0xB1)
        // one could also have used vpsrlq, e.g. v1_0 = _mm512_srli_epi64(v1_0, 32)

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // 5th, 6th, 7th, 8th terms
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i+32));
        v1[1] = _mm512_loadu_si512((const __m512i *) (vec1  +i+48));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+32));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_0+i+48));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+32));
        v2[3] = _mm512_loadu_si512((const __m512i *) (vec2_1+i+48));
        v2[4] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+32));
        v2[5] = _mm512_loadu_si512((const __m512i *) (vec2_2+i+48));

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        v1[0] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v1[1] = _mm512_shuffle_epi32(v1[1], 0xB1);
        v2[0] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[1] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[2] = _mm512_shuffle_epi32(v2[2], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[3], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[4], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[5], 0xB1);

        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[2]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[4]));
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[1]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));

        // split
        dp_hi[0] = _mm512_add_epi64(dp_hi[0], _mm512_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
        dp_lo[0] = _mm512_and_si512(dp_lo[0], low_bits);
        dp_hi[1] = _mm512_add_epi64(dp_hi[1], _mm512_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
        dp_lo[1] = _mm512_and_si512(dp_lo[1], low_bits);
        dp_hi[2] = _mm512_add_epi64(dp_hi[2], _mm512_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
        dp_lo[2] = _mm512_and_si512(dp_lo[2], low_bits);
    }

    // the following loop iterates <= 3 times,
    // each iteration accumulates 2 terms in dp_lo
    for ( ; i+15 < len; i+=16)
    {
        v1[0] = _mm512_loadu_si512((const __m512i *) (vec1  +i));
        v2[0] = _mm512_loadu_si512((const __m512i *) (vec2_0+i));
        v2[1] = _mm512_loadu_si512((const __m512i *) (vec2_1+i));
        v2[2] = _mm512_loadu_si512((const __m512i *) (vec2_2+i));
        v1[1] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[2], 0xB1);

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[1]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[4]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));
    }

    // finally, do a last iteration which may be "incomplete": use mask load
    // (at each position we accumulate 0 or 2 terms, so including the above
    // loop this is a total of <= 8 terms)
    if (i < len)
    {
        // mask == 0b0...01...1 with number of 1's == number of remaining terms
        __mmask16 mask = (1 << (len-i)) - 1;
        v1[0] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec1  +i));
        v2[0] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_0+i));
        v2[1] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_1+i));
        v2[2] = _mm512_maskz_loadu_epi32(mask, (const __m512i *) (vec2_2+i));
        v1[1] = _mm512_shuffle_epi32(v1[0], 0xB1);
        v2[3] = _mm512_shuffle_epi32(v2[0], 0xB1);
        v2[4] = _mm512_shuffle_epi32(v2[1], 0xB1);
        v2[5] = _mm512_shuffle_epi32(v2[2], 0xB1);

        // 1st term: low 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[0], v2[0]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[0], v2[1]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[0], v2[2]));

        // 2nd term: high 32 bit word of each 64 bit word
        dp_lo[0] = _mm512_add_epi64(dp_lo[0], _mm512_mul_epu32(v1[1], v2[3]));
        dp_lo[1] = _mm512_add_epi64(dp_lo[1], _mm512_mul_epu32(v1[1], v2[4]));
        dp_lo[2] = _mm512_add_epi64(dp_lo[2], _mm512_mul_epu32(v1[1], v2[5]));
    }

    // split
    dp_hi[0] = _mm512_add_epi64(dp_hi[0], _mm512_srli_epi64(dp_lo[0], DOT_SPLIT_BITS));
    dp_lo[0] = _mm512_and_si512(dp_lo[0], low_bits);
    dp_hi[1] = _mm512_add_epi64(dp_hi[1], _mm512_srli_epi64(dp_lo[1], DOT_SPLIT_BITS));
    dp_lo[1] = _mm512_and_si512(dp_lo[1], low_bits);
    dp_hi[2] = _mm512_add_epi64(dp_hi[2], _mm512_srli_epi64(dp_lo[2], DOT_SPLIT_BITS));
    dp_lo[2] = _mm512_and_si512(dp_lo[2], low_bits);

    // gather 8 terms in single uint64_t
    uint64_t hsum_lo0 = _mm512_hsum(dp_lo[0]);
    uint64_t hsum_hi0 = _mm512_hsum(dp_hi[0]) + (hsum_lo0 >> DOT_SPLIT_BITS);
    hsum_lo0 &= DOT_SPLIT_MASK;
    uint64_t hsum_lo1 = _mm512_hsum(dp_lo[1]);
    uint64_t hsum_hi1 = _mm512_hsum(dp_hi[1]) + (hsum_lo1 >> DOT_SPLIT_BITS);
    hsum_lo1 &= DOT_SPLIT_MASK;
    uint64_t hsum_lo2 = _mm512_hsum(dp_lo[2]);
    uint64_t hsum_hi2 = _mm512_hsum(dp_hi[2]) + (hsum_lo2 >> DOT_SPLIT_BITS);
    hsum_lo2 &= DOT_SPLIT_MASK;

    NMOD_RED(res[0], pow2_precomp * hsum_hi0 + hsum_lo0, mod);
    NMOD_RED(res[1], pow2_precomp * hsum_hi1 + hsum_lo1, mod);
    NMOD_RED(res[2], pow2_precomp * hsum_hi2 + hsum_lo2, mod);
}

static inline void _avx512_matrix_vector_product(uint32_t * vec_res,
                                                 const uint32_t * mat,
                                                 const uint32_t * vec,
                                                 const uint32_t * dst,
                                                 const uint32_t ncols,
                                                 const uint32_t nrows,
                                                 const nmod_t mod,
                                                 const uint64_t pow2_precomp,
                                                 md_t *st)
{
    slong i = 0;

    for ( ; i+2 < nrows; i+=3)
    {
        int64_t len = ncols - MIN(dst[i], MIN(dst[i+1], dst[i+2]));
        _nmod32_vec_dot3_split_avx512(vec_res+i, vec,
                                      mat + i*ncols,
                                      mat + (i+1)*ncols,
                                      mat + (i+2)*ncols,
                                      len, mod, pow2_precomp);
    }

    if (nrows - i == 2)
    {
        int64_t len = ncols - MIN(dst[i], dst[i+1]);
        _nmod32_vec_dot2_split_avx512(vec_res+i, vec_res+i+1,
                                      vec, mat + i*ncols, mat + (i+1)*ncols,
                                      len, mod, pow2_precomp);
    }
    else if (nrows - i == 1)
        vec_res[i] = _nmod32_vec_dot_split_avx512(vec, mat + i*ncols, ncols - dst[i], mod, pow2_precomp);
}

#endif
