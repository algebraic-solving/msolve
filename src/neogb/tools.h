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

#ifndef GB_TOOLS_H
#define GB_TOOLS_H

#include <time.h>
#include "data.h"

/* cpu time */
double cputime(
    void
    );

/* wall time */
double realtime(
    void
    );

static inline uint8_t mod_p_inverse_8(
        const int16_t val,
        const int16_t p
        )
{
    int16_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 15) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 15) & p;

    return d;
}

static inline uint16_t mod_p_inverse_16(
        const int32_t val,
        const int32_t p
        )
{
    int32_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 31) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 31) & p;

    return (uint16_t)d;
}

/* Barrett reduction: returns x mod p for 0 <= x < 2^64, where
 * m = barrett_multiplier(p) = floor((2^64-1) / p). The quotient estimate
 * is at most one too small, so one correction step suffices. This is
 * much cheaper than a 64 bit division in the linear algebra kernels. */
static inline uint64_t barrett_multiplier(
        const uint64_t p
        )
{
    return UINT64_MAX / p;
}

static inline uint64_t barrett_reduce(
        const uint64_t x,
        const uint64_t p,
        const uint64_t m
        )
{
    const uint64_t q  = (uint64_t)(((unsigned __int128)x * m) >> 64);
    uint64_t r        = x - q * p;
    return r >= p ? r - p : r;
}

static inline uint32_t mod_p_inverse_32(
        const int64_t val,
        const int64_t p
        )
{
    int64_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 63) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 63) & p;

    return d;
}

static inline uint64_t mod_p_inverse_64(
        const int64_t val,
        const int64_t p
        )
{
    int64_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 63) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 63) & p;

    return (uint64_t)d;
}

/* a * b mod p for a, b < p < 2^63 */
static inline uint64_t mul_mod_64(
        const uint64_t a,
        const uint64_t b,
        const uint64_t p
        )
{
    return (uint64_t)(((unsigned __int128)a * b) % p);
}
#endif
