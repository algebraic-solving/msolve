/*
    Copyright (C) 2006, 2007, 2008, 2009, 2016 William Hart
    Copyright (C) 2008, Peter Shrimpton
    Copyright (C) 2009, Tom Boothby
    Copyright (C) 2010, Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ULONG_EXTRAS_H
#define ULONG_EXTRAS_H

#ifdef ULONG_EXTRAS_INLINES_C
#define ULONG_EXTRAS_INLINE FLINT_DLL
#else
#define ULONG_EXTRAS_INLINE static __inline__
#endif

#include <gmp.h>
#include "longlong.h"

#if FLINT64
#define UWORD_MAX_PRIME UWORD(18446744073709551557)
#else
#define UWORD_MAX_PRIME UWORD(4294967291)
#endif

#define UWORD_HALF (UWORD_MAX / 2 + 1)

static __inline__
double n_precompute_inverse(ulong n)
{
   return (double) 1 / (double) n;
}

static __inline__
ulong n_preinvert_limb(ulong n)
{
   ulong norm, ninv;

   count_leading_zeros(norm, n);
   invert_limb(ninv, n << norm);

   return ninv;
}

/* ulong n_ll_mod_preinv(ulong a_hi, ulong a_lo, */
/*                       ulong n, ulong ninv); */
#include "ulong_extras/ll_mod_preinv.c"

static __inline__ 
ulong n_mulmod2_preinv(ulong a, ulong b, ulong n, ulong ninv)
{
    ulong p1, p2;

    /* FLINT_ASSERT(n != 0); */

    umul_ppmm(p1, p2, a, b);
    return n_ll_mod_preinv(p1, p2, n, ninv);
}

static __inline__ 
ulong n_mulmod2(ulong a, ulong b, ulong n)
{
    ulong p1, p2, ninv;

    /* FLINT_ASSERT(n != 0); */

    ninv = n_preinvert_limb(n);
    umul_ppmm(p1, p2, a, b);
    return n_ll_mod_preinv(p1, p2, n, ninv);
}

static __inline__ 
ulong n_addmod(ulong x, ulong y, ulong n)
{
    /* FLINT_ASSERT(x < n); */
    /* FLINT_ASSERT(y < n); */
    /* FLINT_ASSERT(n != 0); */

    return (n - y > x ? x + y : x + y - n);
}

static __inline__
ulong n_submod(ulong x, ulong y, ulong n)
{
    /* FLINT_ASSERT(x < n); */
    /* FLINT_ASSERT(y < n); */
    /* FLINT_ASSERT(n != 0); */
    return (y > x ? x - y + n : x - y);
}


static __inline__
ulong n_negmod(ulong x, ulong n)
{
    /* FLINT_ASSERT(x < n); */
    /* FLINT_ASSERT(n != 0); */

    return n_submod(0, x, n);
}


/* static __inline__ ulong n_gcdinv(ulong * a, ulong x, ulong y); */
#include "ulong_extras/gcdinv.c"

static __inline__
ulong n_invmod(ulong x, ulong y)
{
   ulong r, g;

   g = n_gcdinv(&r, x, y);
   if (g != 1)
      flint_throw(FLINT_IMPINV, "Cannot invert modulo %wd*%wd\n", g, g/y);

   return r;
}



#endif
