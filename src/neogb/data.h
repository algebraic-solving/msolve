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
#ifdef HAVE_OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num(void) { return 0;}
inline omp_int_t omp_get_max_threads(void) { return 1;}
#endif

#define ORDER_COLUMNS 1
/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial.
 * NOTE: if one changes UNROLL then also code needs to be changed:
 * the unrolled loops are hardcoded to 4 at the moment */
#define UNROLL  4
/* we store some more information in the row arrays,
 * real data starts at index OFFSET */
#define OFFSET  5         /* real data starts at OFFSET */
#define LENGTH  OFFSET-1  /* length of the row */
#define PRELOOP OFFSET-2  /* length of not unrolled loop part */
#define COEFFS  OFFSET-3  /* index of corresponding coefficient vector */
#define MULT    OFFSET-4  /* hash of multiplier (for tracing and saturation) */
#define BINDEX  OFFSET-5  /* basis index of element (for tracing) */


/* computational data */
typedef uint8_t cf8_t;   /* coefficient type finite field (8 bit) */
typedef uint16_t cf16_t; /* coefficient type finite field (16 bit) */
typedef uint32_t cf32_t; /* coefficient type finite field (32 bit) */
typedef uint32_t val_t; /* core values like hashes */
typedef val_t hi_t;     /* index of hash table entries*/
typedef hi_t hm_t;      /* hashed monomials for polynomial entries */
typedef uint64_t hl_t;  /* hash table length (maybe >= 2^32) */
/* like exponent hashes, etc. */
typedef uint32_t rba_t; /* reducer binary array */
typedef uint32_t ind_t; /* index in hash table structure */
typedef uint32_t sdm_t;  /* short divmask for faster divisibility checks */
typedef uint32_t len_t; /* length type for different structures */
typedef int16_t exp_t;  /* exponent type */
typedef int32_t deg_t;  /* (total) degree of polynomial */
typedef len_t bi_t;     /* basis index of element */
typedef len_t bl_t;     /* basis load */
typedef len_t pl_t;     /* pair set load */

/* hash data structure */
typedef struct hd_t hd_t;
struct hd_t
{
    val_t val;
    sdm_t sdm;
    deg_t deg;
    ind_t idx;
};

/* hash table data structure */
typedef struct ht_t ht_t;
struct ht_t
{
    exp_t **ev;   /* exponent vector */
    hd_t *hd;     /* hash data */
    hi_t *hmap;   /* hash map */
    hl_t eld;     /* load of exponent vector */
    hl_t esz;     /* size of exponent vector */
    hl_t hsz;     /* size of hash map, might be 2^32 */
    len_t nv;     /* number of variables */
    sdm_t *dm;    /* divisor map for divisibility checks */
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
    spt_t type;
};

typedef struct ps_t ps_t;
struct ps_t
{
    len_t ld;
    len_t sz;
    spair_t *p;
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
    int8_t *red;    /* tracks redundancy of basis elements */
    hm_t **hm;      /* hashed monomials representing exponents */
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
};

/* tracer stuff */
typedef struct primes_t primes_t;
struct primes_t
{
    uint32_t *p;  /* array of primes */
    len_t old;    /* old load of array */
    len_t ld;     /* current load of array */
};

/* represents the trace data of one step of the F4 algorithm */
typedef struct td_t td_t;
struct td_t
{
    len_t *rri;   /* reducer rows information in the format */
                  /* basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    len_t *tri;   /* to be reduced rows information in the format */
                  /* basis index1, multiplier1,
                   * basis index2, multiplier2,... */
    hm_t *lms;    /* hashes of new leading monomials represented */
                  /* in basis hash table */
    rba_t **rba;  /* reducer binary array for each to be reduced row */
    len_t rld;    /* load of reducer rows information*/
    len_t tld;    /* load of to be reduced rows information*/
    len_t nlm;    /* number of new leading monomials in this step */
};

typedef struct trace_t trace_t;
struct trace_t
{
    td_t *td;     /* array of trace data for each round of F4 */
    len_t ld;     /* load of trace data */
    len_t sz;     /* size allocated for trace data */
    bl_t *lmps;   /* position of non-redundant lead monomials in basis */
    sdm_t *lm;    /* non-redundant lead monomials as short divmask */
    bl_t lml;     /* number of lead monomials of non redundant
                     elements in basis */
    len_t *rd;    /* rounds in which saturation steps lead to
                   * non-trivial kernels */
    deg_t *deg;   /* degree for multipliers in saturation step */
    len_t rld;    /* load of rounds stored, i.e. how often do saturate */
    len_t rsz;    /* size of rounds stored */
};


/* statistic stuff */
typedef struct stat_t stat_t;
struct stat_t
{
    double round_ctime;
    double select_ctime;
    double symbol_ctime;
    double la_ctime;
    double update_ctime;
    double convert_ctime;
    double overall_ctime;
    double reduce_gb_ctime;
    double rht_ctime;

    double round_rtime;
    double select_rtime;
    double symbol_rtime;
    double la_rtime;
    double update_rtime;
    double convert_rtime;
    double overall_rtime;
    double reduce_gb_rtime;
    double rht_rtime;

    int64_t num_pairsred;
    int64_t num_gb_crit;
    int64_t num_redundant_old;
    int64_t num_redundant;
    int64_t num_rht;
    int64_t num_rowsred;
    int64_t num_zerored;

    int32_t ngens;
    int32_t nvars;
    int32_t mnsel;
    int32_t homogeneous;
    uint32_t fc;
    int32_t mo;
    int32_t laopt;
    int32_t init_hts;
    int32_t nthrds;
    int32_t reset_ht;
    int32_t current_rd;
    int32_t current_deg;
    uint64_t max_bht_size;
    uint64_t max_sht_size;
    uint64_t max_uht_size;
    int64_t nterms_basis;
    int32_t size_basis;
    int32_t ff_bits;
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

    /* for f4sat */
    uint32_t new_multipliers;
    uint32_t nr_kernel_elts;
};

/* function pointers */
extern bs_t *(*initialize_basis)(
        const int32_t ngens
        );
extern void (*check_enlarge_basis)(
        bs_t *bs,
        const len_t added
        );
extern void (*normalize_initial_basis)(
        bs_t *bs,
        const uint32_t fc
        );
extern bs_t *(*copy_basis_mod_p)(
        const bs_t * const gbs,
        const stat_t * const st
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

extern void (*import_julia_data)(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

extern int64_t (*export_julia_data)(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht,
        const uint32_t fc
        );

/* linear algebra routines */
extern void (*linear_algebra)(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        );

extern int (*application_linear_algebra)(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        );

extern void (*trace_linear_algebra)(
        trace_t *trace,
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        );

extern void (* interreduce_matrix_rows)(
        mat_t *mat,
        bs_t *bs,
        stat_t *st
        );

extern cf32_t *(*reduce_dense_row_by_old_pivots_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hi_t dpiv,
        const uint32_t fc
        );

extern hm_t *(*reduce_dense_row_by_known_pivots_sparse_ff_32)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hi_t dpiv,
        const hm_t tmp_pos,
        stat_t *st
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
        stat_t *st
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
