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
#include <sys/time.h>
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
#endif
