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

//-1 sur 32 bits
#define MONE32 ((uint32_t)0xFFFFFFFF)

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

#define AVXINNERPROD(res1, mat8, vec8, va, vtmp, prod1, prod2)        \
  {                                                                   \
    mat8=AVX2LOADU(va);                                           \
    vec8=AVX2LOADU(vtmp);                                         \
    /* four 32-bits mul, lower parts */                           \
    /* four 32-bits mul, higher parts */                          \
    prod1=AVX2MUL(mat8,vec8);                                     \
    prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));    \
    res1 = prod1 + prod2;                                         \
    va+=8;                                                      \
    vtmp+=8;                                                    \
  }                                                               \


#define AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2)  \
  {                                                                 \
    mat8=AVX2LOADU(va);                                             \
    vec8=AVX2LOADU(vtmp);                                           \
    /* four 32-bits mul, lower parts */                             \
    /* four 32-bits mul, higher parts */                            \
    prod1=AVX2MUL(mat8,vec8);                                       \
    prod2= AVX2MUL(AVX2SRLI_64(mat8,32),AVX2SRLI_64(vec8,32));      \
    res1 += prod1 + prod2;                                          \
    va+=8;                                                        \
    vtmp+=8;                                                      \
  }                                                                 \

#endif



/* /\* */
/* Barrett reduction */
/*   2^30 < p < 2^31 */

/*   r=31 and n=32 */

/*   alpha = 2^(n-1) = 2^31, t=n=32, s=r-1=30 */

/*   preinv = q = floor( 2^(s+t) / p ) */

/*   const uint32_t n, const uint32_t p, */
/*   const uint32_t preinv */
static inline uint64_t _MODRED32(const uint64_t n, const uint32_t p, const uint32_t preinv){
  /* const uint32_t s = 30; */
  /* const uint32_t t = 32; */

  /* const uint64_t t1 = n >> s; */
  const uint64_t t1 = n >> 30;
  /* const uint64_t t2 = (t1 * preinv) >> t; */
  const uint64_t t2 = (t1 * preinv) >> 32;
  int64_t r = n - t2 * p;
  (r >= p) ? r -= p : r;
  (r >= p) ? r -= p : r;
  return r;
}

static inline uint64_t MODRED32(const uint64_t n, const uint32_t p, const uint32_t preinv){
  return ((n > ((uint64_t)p)<<31)) ? MODRED32(n-(((uint64_t)p)<<31), p, preinv) : _MODRED32(n, p, preinv);
}

static inline uint64_t ADDMODRED32(uint64_t a, uint64_t b,
                                   const uint32_t p, const uint32_t preinv)
{
  const long neg = p - a;
  return (neg > b) ? MODRED32((a + b), p, preinv) : MODRED32((b - neg), p, preinv);
}

static inline uint64_t SUBMODRED32(uint64_t a, uint64_t b,
                                   const uint32_t p, const uint32_t preinv)
{
  const long neg = a - b;
  return (a < b) ? (p + neg) : (neg);
}

#if HAVE_AVX2
static inline void REDUCE(uint64_t *acc64, uint64_t *acc4x64,
                          __m256i acc_low, __m256i acc_high,
                          const uint32_t fc, const uint32_t preinv,
                          const uint32_t RED_32, const uint32_t RED_64){
    AVX2STOREU(acc4x64,acc_low);
    AVX2STOREU(acc4x64+4,acc_high);
    /* Reduction */
    *acc64=0;
    acc4x64[0]+=MODRED32((acc4x64[4]>>32)*RED_64,fc,preinv);
    /* partie basse du registre haut (2^32->2^63) */
    acc4x64[0]+=MODRED32((acc4x64[4]&((uint64_t)0xFFFFFFFF))*RED_32,fc,preinv);
    *acc64+=MODRED32(acc4x64[0],fc,preinv);

    acc4x64[1]+=MODRED32((acc4x64[5]>>32)*RED_64,fc,preinv);
    /* partie basse du registre haut (2^32->2^63) */
    acc4x64[1]+=MODRED32((acc4x64[5]&((uint64_t)0xFFFFFFFF))*RED_32,fc,preinv);
    *acc64+=MODRED32(acc4x64[1],fc,preinv);

    acc4x64[2]+=MODRED32((acc4x64[6]>>32)*RED_64,fc,preinv);
    /* partie basse du registre haut (2^32->2^63) */
    acc4x64[2]+=MODRED32((acc4x64[6]&((uint64_t)0xFFFFFFFF))*RED_32,fc,preinv);
    *acc64+=MODRED32(acc4x64[2],fc,preinv);

    acc4x64[3]+=MODRED32((acc4x64[7]>>32)*RED_64, fc, preinv);
    /* partie basse du registre haut (2^32->2^63) */
    acc4x64[3]+=MODRED32((acc4x64[7]&((uint64_t)0xFFFFFFFF))*RED_32,fc,preinv);
    *acc64+=MODRED32(acc4x64[3],fc,preinv);
    *acc64=MODRED32(*acc64,fc,preinv);
}
#endif

/* B has l rows, n columns */
/* stores in tmp the entries of B column-wise */
/* tmp = [ B[0][0], B[1][0], ..., B[n][0],
           B[0][1], ..., B[n][1], ... ] */
static __inline__ void specialize_tmp(uint32_t *tmp,
                                      uint32_t ** B,
                                      const long l, const long n){
  for(long i = 0; i < l; i++){
    for(long j = 0; j < n; j++)
      tmp[j*l + i] = B[i][j];
  }
}

static __inline__ void make_zero(uint32_t ** D,
                                long m, long n){
  for(long i = 0; i < m; i++){
    for(long j = 0; j < n; j++){
      D[i][j] = 0;
    }
  }

}

#ifdef HAVE_AVX2
/*
  m = number of rows in A
  l = number of cols in A = number of rows in B
  n = number of cols in B

  D = A*B

*/
void _mod_mat_addmul_transpose_op(uint32_t *D, 
                                  uint32_t *A, uint32_t *B,
                                  const uint32_t m, const uint32_t l,
                                  const uint32_t n,
                                  const uint32_t fc, const uint32_t preinv,
                                  const uint32_t RED_32, const uint32_t RED_64){

  __m256i acc_low, acc_high, mask, vec8, mat8;
  __m256i prod1,prod2;
  __m256i res1;
#if IS_ALIGNED_TMP
  uint64_t acc4x64[8] ALIGNED32;
#else
  uint64_t acc4x64[8];
#endif
  uint64_t acc64 = 0;

  mask=AVX2SET1_64(MONE32);

  const uint32_t *va;
  const uint32_t *vtmp;

  acc_low=AVX2SETZERO();
  acc_high=AVX2SETZERO();

  //quo(l, 2^6) * 2^6
  const uint32_t bl = ((l >> 6) << 6);

  for(uint32_t i = 0; i < m; i++){
    for(uint32_t j = 0; j < n; j++){
      acc64 = 0;

      acc_low=AVX2SETZERO();
      acc_high=AVX2SETZERO();
      for(uint32_t k = 0; k < bl; k+=64){

        /* va = A[i]+k; */
        va = A+(i*l)+k;
        vtmp = B+j*l+k;

        AVXINNERPROD(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);

        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);
        AVXINNERPRODLOOP(res1, mat8, vec8, va, vtmp, prod1, prod2);

        acc_low = AVX2ADD_64(acc_low,AVX2AND_(res1,mask));
        acc_high = AVX2ADD_64(acc_high,AVX2SRLI_64(res1,32));

      }

      REDUCE(&acc64, acc4x64, acc_low, acc_high, fc, preinv,
             RED_32, RED_64);

      uint64_t p = 0, u = 0, d = 0;
      uint64_t r, Shi = 0, Slo = 0;
      const long c2 = j*l;
      const long c1 = i*l;
      for(uint32_t k = bl; k < l; k++){

        /* p = ((uint64_t)A[i][k] * (uint64_t)tmp[c+k]); */
        p = ((uint64_t)A[c1+k] * (uint64_t)B[c2+k]);
        u = p>>32;
        d = p - (u<<32);
        Shi += u;
        Slo += d;

      }
      uint64_t hi, lo;
      hi = Shi>>32;
      lo = Shi - (hi<<32);
      r = MODRED32(hi<<32, fc, preinv);
      r = ADDMODRED32(lo, r, fc, preinv);
      r = MODRED32((r<<32), fc, preinv);
      r = ADDMODRED32(Slo, r, fc, preinv);

      acc64 = ADDMODRED32(acc64, r, fc, preinv);

      D[i*n+j] = ADDMODRED32(D[i*n+j], acc64, fc, preinv);

    }
  }

}
#endif

static inline void sparse_matfglm_mul(CF_t *res, sp_matfglm_t *matxn, CF_t *R,
                                      CF_t *tres,
                                      const int nc,
                                      const mod_t prime, const uint32_t preinv,
                                      const uint32_t RED_32,
                                      const uint64_t RED_64){
  szmat_t ncols = matxn->ncols;
  szmat_t nrows = matxn->nrows;
  szmat_t ntriv = ncols - nrows;

  for(szmat_t j = 0; j < ntriv; j++){
    for(int i = 0; i < nc; i++){
      uint32_t idx = i*(matxn->ncols);
      res[idx + matxn->triv_idx[j]] = R[idx + matxn->triv_pos[j]];
    }
  }


  /* real product */
#ifdef HAVE_AVX2
  _mod_mat_addmul_transpose_op(tres, matxn->dense_mat, R,
                               matxn->nrows, matxn->ncols, nc,
                               prime, preinv,
                               RED_32, RED_64);
#else
  fprintf(stderr, "Not implemented yet\n");
  exit(1);
#endif

  for(szmat_t j = 0; j < nrows; j++){
    for(int i = 0; i < nc; i++){
      uint32_t idx = i*(matxn->ncols);
      res[idx + matxn->dense_idx[j]] = tres[idx + j];
    }
  }
}
