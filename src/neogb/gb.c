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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "libneogb.h"
#include "data.c"
#include "meta_data.c"/* computational meta data */
#include "tools.c"    /* tools like inversion mod p,
                       * tracer construction, timings etc. */
#include "sort_r.h"   /* special quicksort implementation */
#include "hash.c"     /* hash table stuff */
#include "order.c"    /* order and comparison procedures */
#include "basis.c"    /* basis and polynomial handling */
#include "la_ff_8.c"  /* finite field linear algebra (8 bit) */
#include "la_ff_16.c" /* finite field linear algebra (16 bit) */
#include "la_ff_32.c" /* finite field linear algebra (32 bit) */
#include "la_qq.c"    /* rational linear algebra */
#include "update.c"   /* update process and pairset handling */
#include "convert.c"  /* conversion between hashes and column indices*/
#include "symbol.c"   /* symbolic preprocessing */
#include "io.c"       /* input and output data handling */
#include "engine.c"   /* global, shared parts of gb engine */
#include "f4.c"       /* implemenation of f4 algorithm */
#include "sba.c"      /* implemenation of sba algorithm */
#include "nf.c"       /* implemenation of normal form algorithm */
#include "f4sat.c"    /* implemenation of saturation algorithm */
#include "modular.c"  /* implemenation of modular Groebner for F4 */

