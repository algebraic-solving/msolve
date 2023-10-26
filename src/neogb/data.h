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


#ifndef GB_DATA_H
#define GB_DATA_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <string.h> /* for memset et al. */
#include <limits.h>
#include <math.h>

/* check if OpenMP is available */
#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num(void) { return 0;}
inline omp_int_t omp_get_max_threads(void) { return 1;}
#endif

#define PARALLEL_HASHING 1
#define ORDER_COLUMNS 1
/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial.
 * NOTE: if one changes UNROLL then also code needs to be changed:
 * the unrolled loops are hardcoded to 4 at the moment */
#define UNROLL  4
/* we store some more information in the row arrays,
 * real data starts at index OFFSET */
#define OFFSET  6           /* real data starts at OFFSET */
#define LENGTH  (OFFSET-1)  /* length of the row */
#define PRELOOP (OFFSET-2)  /* length of not unrolled loop part */
#define COEFFS  (OFFSET-3)  /* index of corresponding coefficient vector */
#define MULT    (OFFSET-4)  /* hash of multiplier (for tracing and saturation) */
#define BINDEX  (OFFSET-5)  /* basis index of element (for tracing) */
#define DEG     (OFFSET-6)  /* the first entry in each exponent vector
                             * stores the total degree of the polynomial */

/* there is a different prelude with meta data for signature based matrices */
#define SM_OFFSET  5            /* real data starts at SIGOFFSET for signature
                                 * based comptutations */
#define SM_LEN   (SM_OFFSET-1)  /* signature meta data length of polynomial */
#define SM_PRE   (SM_OFFSET-2)  /* signature meta data preloop of polynomial */
#define SM_CFS   (SM_OFFSET-3)  /* index of corresponding coefficient array */
#define SM_SIDX  (SM_OFFSET-4)  /* index of signautre */
#define SM_SMON  (SM_OFFSET-5)  /* hash value of signature monomial */

/* computational data */
typedef uint8_t cf8_t;   /* coefficient type finite field (8 bit) */
typedef uint16_t cf16_t; /* coefficient type finite field (16 bit) */
typedef uint32_t cf32_t; /* coefficient type finite field (32 bit) */
typedef uint32_t val_t;  /* core values like hashes */
typedef val_t hi_t;      /* index of hash table entries*/
typedef hi_t hm_t;       /* hashed monomials for polynomial entries */
typedef hm_t sm_t;       /* hashed monomial of signature */
typedef uint16_t si_t;   /* index of signature */
typedef uint64_t hl_t;   /* hash table length (maybe >= 2^32) */
/* like exponent hashes, etc. */
typedef uint32_t rba_t;  /* reducer binary array */
typedef uint32_t ind_t;  /* index in hash table structure */
typedef uint32_t sdm_t;  /* short divmask for faster divisibility checks */
typedef uint32_t len_t;  /* length type for different structures */
typedef uint16_t exp_t;  /* exponent type */
typedef int32_t deg_t;   /* (total) degree of polynomial */
typedef len_t bi_t;      /* basis index of element */
typedef len_t bl_t;      /* basis load */
typedef len_t pl_t;      /* pair set load */

/* hash data structure */
typedef struct hd_t hd_t;
struct hd_t
{
    val_t val;
    sdm_t sdm;
    ind_t idx;
    deg_t deg;
};

/* 
 * Exponent vectors look the following for n variables:
 *
 * 1. If we use a non-block monomial order like DRL
 * [deg, exp_v1, ..., exp_vn]
 * -> length is n+1
 *
 * 2. If we use a block elimination order with two blocks
 * of k variables and n-k variables
 * [deg_b1, exp_v1, ..., exp_vk, deg_b2, exp_vk+1, ..., exp_vn]
 * -> length is n+2
 *
 *  In any of the above situations nv will be n.
 *  evl will be nv + 1 + (ebl != 0) where ebl is the number
 *  of variables in the first variable block + 1 (for the degree of
 *  this block) if we use an elimination block order, 0 otherwise.
 *  */

/* hash table data structure */
typedef struct ht_t ht_t;
struct ht_t
{
    exp_t **ev;   /* exponent vector */
    hd_t *hd;     /* hash data */
    hi_t *hmap;   /* hash map */
    len_t elo;    /* load of exponent vector before current step */
    hl_t eld;     /* load of exponent vector */
    hl_t esz;     /* size of exponent vector */
    hl_t hsz;     /* size of hash map, might be 2^32 */
    len_t ebl;    /* elimination block length:
                   * degree + #elimination variables,
                   * 0 if no elimination order */
    len_t nv;     /* number of variables */
    len_t evl;    /* real length of exponent vector,
                   * includes degree (or two degrees
                   * if an elimination order is used) */
    sdm_t *dm;    /* divisor map for divisibility checks */
    len_t *dv;    /* variables for divmask */
    len_t ndv;    /* number of variables for divmask */
    len_t bpv;    /* bits per variable in divmask */
    val_t *rn;    /* random numbers for hash generation */
    uint32_t rsd; /* seed for random number generator */
};

/* S-pair types */
typedef enum {S_PAIR, GCD_PAIR, GEN_PAIR} spt_t;
typedef struct spair_t spair_t;
struct spair_t
{
    hi_t lcm;
    bi_t gen1;
    bi_t gen2;
    deg_t deg;
    spt_t type;
};

typedef struct ps_t ps_t;
struct ps_t
{
    len_t ld;
    len_t sz;
    spair_t *p;
};

/* signature criteria structure,
 * at the moment we only store the lead terms (i e. signatures) */
typedef struct crit_t crit_t;
struct crit_t
{
    sdm_t *sdm; // array of shot divisor mask of signature monomial
    hm_t *hm;   // array of hash value of the of the signature monomial
    len_t ld;   // load of the corresponding arrays
    len_t sz;   // allocated memory / size of the corresponding arrays
};

/* basis stuff */
typedef struct bs_t bs_t;
struct bs_t
{
    bl_t ld;        /* load of basis */
    bl_t sz;        /* size allocated for basis */
    bl_t lo;        /* load before current update */
    bl_t constant;  /* 1 if constant is found in basis */
    deg_t mltdeg;   /* maximal appearing degree in lead term in basis */
    bl_t *lmps;     /* position of non-redundant lead monomials in basis */
    sdm_t *lm;      /* non-redundant lead monomials as short divmask */
    bl_t lml;       /* number of lead monomials of non redundant
                       elements in basis */
    ht_t *ht;       /* hash table for basis elements */
    int8_t *red;    /* tracks redundancy of basis elements */
    hm_t **hm;      /* hashed monomials representing exponents */
    sm_t *sm;       /* signatures for F5-style computations */
    si_t *si;       /* signatures index for F5-style computations */
    cf8_t **cf_8;   /* coefficients for finite fields (8 bit) */
    cf16_t **cf_16; /* coefficients for finite fields (16 bit) */
    cf32_t **cf_32; /* coefficients for finite fields (32 bit) */
    mpz_t **cf_qq;  /* coefficients for rationals (always multiplied such that
                       the denominator is 1) */
};

/* matrix stuff */
typedef struct mat_t mat_t;
struct mat_t
{
    hm_t **tr;          /* rows to be reduced of the matrix, only column */
                        /* entries, coefficients are handled via linking */
                        /* to coefficient arrays */
    rba_t **rba;        /* bit array for each row to be reduced storing */
                        /* if a reducer row is used during reduction. */
                        /* Thus we can reconstruct the trace for the */
                        /* rows not reduced to zero easily by bitwise */
                        /* AND operations on these bit arrays. */
    hm_t **rr;          /* reducer rows of the matrix, only column */
                        /* entries, coefficients are handled via linking */
                        /* to coefficient arrays. */
    cf8_t **cf_8;       /* coefficients for finite fields (8 bit) */
    cf16_t **cf_16;     /* coefficients for finite fields (16 bit) */
    cf32_t **cf_32;     /* coefficients for finite fields (32 bit) */
    mpz_t **cf_qq;      /* coefficients for rationals */
    mpz_t **cf_ab_qq;   /* coefficients for rationals */
    len_t sz;           /* number of rows allocated resp. size */
    len_t np;           /* number of new pivots */
    len_t nr;           /* number of rows set */
    len_t nc;           /* number of columns */
    len_t nru;          /* number of upper rows (in ABCD splicing) */
    len_t nrl;          /* number of lower rows (in ABCD splicing) */
    len_t ncl;          /* number of left columns (in ABCD splicing) */
    len_t ncr;          /* number of right columns (in ABCD splicing) */
    len_t rbal;         /* length of reducer binary array */
    deg_t cd;           /* current degree */
};

/* signature matrix stuff, stores information from previous and current step */
typedef struct smat_t smat_t;
struct smat_t
{
    hm_t **cr;     /* current matrix rows, hashes resp. columns */
    hm_t **pr;     /* previous matrix rows, hashes resp. columns */
    cf32_t **cc32; /* current matrix coefficients */
    cf32_t **pc32; /* previous matrix coefficients */
    deg_t cd;      /* current degree */
    len_t nlm;     /* number of new leading monomials for basis from */
                   /* current matrix */
    len_t csz;     /* number of rows memory is allocated */
                   /* for in current matrix */
    len_t cld;     /* number of rows stored in current matrix */
    len_t pld;     /* number of rows stored in previous matrix */
    len_t nc;      /* number of columns in current matrix */
    len_t nz;      /* number of zero reductions during linear algebra */
                   /* on current matrix */
};

/* tracer stuff */
typedef struct primes_t primes_t;
struct primes_t
{
    uint32_t *p;  /* array of primes */
    len_t old;    /* old load of array */
    len_t ld;     /* current load of array */
};

/* represents the trace data for one saturation step */
typedef struct ts_t ts_t;
struct ts_t
{
    len_t *rri;   /* reducer rows information in the format
                   * basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    len_t *tri;   /* to be reduced rows information in the format */
                  /* basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    len_t lml;    /* number of non-redundant elements in basis */
    deg_t deg;    /* minimal degree to start saturation process */
    len_t f4rd;   /* round of f4 */
    hm_t *nlms;   /* hashes of new leading monomials represented
                   * in basis hash table */
    len_t rld;    /* load of reducer rows information*/
    len_t tld;    /* load of to be reduced rows information*/
    len_t nlm;    /* number of new leading monomials in this step */
};

typedef struct td_t td_t;
struct td_t
{
    len_t *rri;   /* reducer rows information in the format */
                  /* basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    len_t *tri;   /* to be reduced rows information in the format */
                  /* basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    hm_t *nlms;   /* hashes of new leading monomials represented
                   * in basis hash table */
    rba_t **rba;  /* reducer binary array for each to be reduced row */
    deg_t deg;    /* degree of elements in trace */
    len_t rld;    /* load of reducer rows information*/
    len_t tld;    /* load of to be reduced rows information*/
    len_t nlm;    /* number of new leading monomials in this step */
};

/* possible trace levels */
typedef enum {NO_TRACER, LEARN_TRACER, APPLY_TRACER} tl_t;
typedef struct trace_t trace_t;
struct trace_t
{
    td_t *td;     /* array of trace data for each round of F4 */
    ts_t *ts;     /* array of trace data for each saturation step */
    len_t ltd;    /* load of trace data td */
    len_t lts;    /* load of trace data ts */
    len_t std;    /* size allocated for trace data td */
    len_t sts;    /* size allocated for trace data ts */
    sdm_t *lm;    /* final minimal leading ideal represented as
                     short divisor masks */
    bl_t *lmps;   /* minimal basis geneator positions */
    hm_t *lmh;    /* minimal basis leading monomial hashes represented
                     in basis hash table */
    bl_t lml;     /* final number of lead monomials of non redundant
                     elements in basis */
    len_t *rd;    /* rounds in which saturation steps lead to
                   * non-trivial kernels */
    len_t rld;    /* load of rounds stored, i.e. how often do saturate */
    len_t rsz;    /* size of rounds stored */
};


/* meta data stuff */
typedef struct md_t md_t;
struct md_t
{
    /* trace data */
    trace_t *tr;
    tl_t trace_level;
    int32_t trace_rd;

    /* hash table data */
    ht_t *ht;

    len_t np; /* new pivots */

    hi_t *hcm;
    ps_t *ps;

    double round_ctime;
    double select_ctime;
    double symbol_ctime;
    double la_ctime;
    double update_ctime;
    double convert_ctime;
    double f4_ctime;
    double reduce_gb_ctime;
    double tracer_ctime;
    double rht_ctime;
    double nf_ctime;

    double round_rtime;
    double select_rtime;
    double symbol_rtime;
    double la_rtime;
    double update_rtime;
    double convert_rtime;
    double f4_rtime;
    double reduce_gb_rtime;
    double tracer_rtime;
    double rht_rtime;
    double nf_rtime;

    int64_t num_pairsred;
    int64_t num_gb_crit;
    int64_t num_syz_crit;
    int64_t num_rew_crit;
    int64_t num_redundant_old;
    int64_t num_redundant;
    int64_t num_rht;
    int64_t num_rowsred;
    int64_t num_zerored;
    int64_t mat_max_nrows;
    int64_t mat_max_ncols;
    double  mat_max_density;

    int32_t ngens_input;
    int32_t ngens_invalid;
    int32_t ngens;
    int32_t init_bs_sz;
    int32_t nvars;
    int32_t mnsel;
    int32_t homogeneous;
    uint32_t gfc; /* global field characteristic */
    uint32_t fc;
    int32_t nev; /* number of elimination variables */
    int32_t mo; /* monomial ordering: 0=DRL, 1=LEX*/
    int32_t laopt;
    int32_t init_hts;
    int32_t nthrds;
    int32_t reset_ht;
    int32_t current_rd;
    int32_t current_deg;
    deg_t max_gb_degree;
    uint64_t max_bht_size;
    uint64_t max_sht_size;
    uint64_t max_uht_size;
    int64_t nterms_basis;
    int32_t size_basis;
    int32_t ff_bits;
    int32_t nf;
    int32_t use_signatures; /* module monomial ordering:
                               0 = off,
                               1=SCHREYER,
                               2=POT,
                               3=DPOT */
    int32_t reduce_gb;

    uint32_t prime_start;
    int32_t nprimes;

    int32_t info_level;
    int32_t gen_pbm_file;

    double trace_nr_mult;
    double trace_nr_add;
    uint64_t  trace_nr_red;
    double application_nr_mult;
    double application_nr_add;
    uint64_t application_nr_red;

    int32_t print_gb;

    /* for f4sat */
    uint32_t new_multipliers;
    uint32_t nr_kernel_elts;
};

/* function pointers */
/* extern bs_t *(*initialize_basis)(
 *         const int32_t ngens
 *         ); */
extern void (*normalize_initial_basis)(
        bs_t *bs,
        const uint32_t fc
        );

extern int (*initial_input_cmp)(
        const void *a,
        const void *b,
        void *ht
        );

extern int (*initial_gens_cmp)(
        const void *a,
        const void *b,
        void *ht
        );

extern int (*monomial_cmp)(
        const hi_t a,
        const hi_t b,
        const ht_t *ht
        );

extern int (*spair_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

extern int (*hcm_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

/* linear algebra routines */
extern void (*sba_linear_algebra)(
        smat_t *smat,
        crit_t *syz,
        md_t *st,
        const ht_t * const ht
        );

extern void (*linear_algebra)(
        mat_t *mat,
        const bs_t * const tbr,
        const bs_t * const bs,
        md_t *st
        );

extern int (*application_linear_algebra)(
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        );

extern void (*trace_linear_algebra)(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        md_t *st
        );

extern void (* interreduce_matrix_rows)(
        mat_t *mat,
        bs_t *bs,
        md_t *st,
        int free_basis
        );

extern cf32_t *(*reduce_dense_row_by_old_pivots_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,
        const uint32_t fc
        );

extern hm_t *(*sba_reduce_dense_row_by_known_pivots_sparse_ff_32)(
        int64_t *dr,
        smat_t *smat,
        hm_t *const *pivs,
        const hi_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t sm,      /* signature monomial of row reduced */
        const len_t si,     /* signature index of row reduced */
        const len_t ri,     /* index of row in matrix */
        md_t *st
        );

extern hm_t *(*reduce_dense_row_by_known_pivots_sparse_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        const len_t mh,     /* multiplier hash for tracing */
        const len_t bi,     /* basis index of generating element */
        const len_t tr,     /* trace data? */
        md_t *st
        );

extern hm_t *(*trace_reduce_dense_row_by_known_pivots_sparse_ff_32)(
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
        );

extern cf32_t *(*reduce_dense_row_by_all_pivots_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const *pivs,
        cf32_t *const *dpivs,
        const uint32_t fc
        );


extern cf32_t *(*reduce_dense_row_by_dense_new_pivots_ff_32)(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const uint32_t fc
        );

#endif
