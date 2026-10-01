/* This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA

   Original implementation was copyright (C) 1997 Makoto Matsumoto and
   Takuji Nishimura. Coded by Takuji Nishimura, considering the
   suggestions by Topher Cooper and Marc Rieffel in July-Aug. 1997, "A
   C-program for MT19937: Integer version (1998/4/6)"

   This implementation copyright (C) 1998 Brian Gough. I reorganized
   the code to use the module framework of GSL.  The license on this
   implementation was changed from LGPL to GPL, following paragraph 3
   of the LGPL, version 2.

   Update:

   The seeding procedure has been updated to match the 10/99 release
   of MT19937.

   Update:

   The seeding procedure has been updated again to match the 2002
   release of MT19937

   The original code included the comment: "When you use this, send an
   email to: matumoto@math.keio.ac.jp with an appropriate reference to
   your work".

   Makoto Matsumoto has a web page with more information about the
   generator, http://www.math.keio.ac.jp/~matumoto/emt.html. 

   The paper below has details of the algorithm.

   From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
   623-dimensionally equidistributerd uniform pseudorandom number
   generator". ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1 (Jan. 1998), Pages 3-30

   You can obtain the paper directly from Makoto Matsumoto's web page.

   The period of this generator is 2^{19937} - 1.

*/

/* This file is derived from rng/mt.c of GSL 1.9, the last GSL release
   under GPL version 2 or later.

   It was modified for msolve on 2026-10-01 as follows:

   - The GSL includes and the gsl_rng_type framework were removed.
   - The 1998 and 1999 seeding procedures were removed.
   - The state type mt_state_t and the functions mt_get and mt_set are
     declared in mt.h. Hence, mt_get and mt_set are no longer static.
   - The state words, as well as all other unsigned long variables,
     parameters and return values, are uint32_t instead of unsigned long.
     The UL constants are then U constants, and the masks with 0xffffffff
     in mt_set were removed.
   - The local macros are undefined at the end of the file.

*/

#include "mt.h"

#define N 624   /* Period parameters */
#define M 397

/* most significant w-r bits */
static const uint32_t UPPER_MASK = 0x80000000U;   

/* least significant r bits */
static const uint32_t LOWER_MASK = 0x7fffffffU;   

uint32_t
mt_get (void *vstate)
{
  mt_state_t *state = (mt_state_t *) vstate;

  uint32_t k ;
  uint32_t *const mt = state->mt;

#define MAGIC(y) (((y)&0x1) ? 0x9908b0dfU : 0)

  if (state->mti >= N)
    {   /* generate N words at one time */
      int kk;

      for (kk = 0; kk < N - M; kk++)
        {
          uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);
        }
      for (; kk < N - 1; kk++)
        {
          uint32_t y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
          mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
        }

      {
        uint32_t y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
      }

      state->mti = 0;
    }

  /* Tempering */
  
  k = mt[state->mti];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680U;
  k ^= (k << 15) & 0xefc60000U;
  k ^= (k >> 18);

  state->mti++;

  return k;
}

void
mt_set (void *vstate, uint32_t s)
{
  mt_state_t *state = (mt_state_t *) vstate;
  int i;

  if (s == 0)
    s = 4357;   /* the default seed is 4357 */

  state->mt[0]= s;

  for (i = 1; i < N; i++)
    {
      /* See Knuth's "Art of Computer Programming" Vol. 2, 3rd
         Ed. p.106 for multiplier. */

      state->mt[i] =
        (1812433253U * (state->mt[i-1] ^ (state->mt[i-1] >> 30)) + i);
    }

  state->mti = i;
}

#undef MAGIC
#undef M
#undef N
