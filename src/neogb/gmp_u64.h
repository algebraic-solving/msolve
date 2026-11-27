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

#ifndef GMP_U64_H
#define GMP_U64_H

/* GMP functions with an operand of 64 bits. The functions mpz_*_ui of GMP
 * take or return an unsigned long, which has only 32 bits on Windows, so
 * that they cannot be used with the primes of more than 32 bits (e.g. the
 * primes of 52 bits of the Groebner bases over the rationals). The
 * operand is passed as a limb instead, the limbs having 64 bits on all
 * the 64 bit platforms. */

#include <stdint.h>
#include <gmp.h>

#if GMP_NUMB_BITS < 64
#error "msolve needs GMP limbs of at least 64 bits"
#endif

/* r = v */
static inline void mpz_set_u64(mpz_t r, const uint64_t v)
{
    const mp_limb_t l = v;
    mpz_t t;
    mpz_roinit_n(t, &l, v != 0);
    mpz_set(r, t);
}

/* initializes r to v */
static inline void mpz_init_set_u64(mpz_t r, const uint64_t v)
{
    mpz_init(r);
    mpz_set_u64(r, v);
}

/* r = a * v */
static inline void mpz_mul_u64(mpz_t r, const mpz_t a, const uint64_t v)
{
    const mp_limb_t l = v;
    mpz_t t;
    mpz_roinit_n(t, &l, v != 0);
    mpz_mul(r, a, t);
}

/* r = r + a * v */
static inline void mpz_addmul_u64(mpz_t r, const mpz_t a, const uint64_t v)
{
    const mp_limb_t l = v;
    mpz_t t;
    mpz_roinit_n(t, &l, v != 0);
    mpz_addmul(r, a, t);
}

/* a mod d in [0, d), for d > 0 (as mpz_fdiv_ui) */
static inline uint64_t mpz_fdiv_u64(const mpz_t a, const uint64_t d)
{
    const size_t n = mpz_size(a);
    if (n == 0) {
        return 0;
    }
    mp_limb_t r = mpn_mod_1(mpz_limbs_read(a), (mp_size_t)n, (mp_limb_t)d);
    if (mpz_sgn(a) < 0 && r != 0) {
        r = d - r;
    }
    return r;
}

#endif
