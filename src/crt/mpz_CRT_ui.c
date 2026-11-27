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

/**

This implementation is a (very) slight modification of the functions in FLINT.

 **/

/* Force the real FLINT ulong_extras.h to be processed first: it and our
 * local ulong_extras.h below share the ULONG_EXTRAS_H include guard, so
 * whichever is seen first "wins". Without this, our copy only works by
 * accident of #include order elsewhere pulling in FLINT's real header
 * first; if it doesn't, our copy's own "#include
 * ulong_extras/ll_mod_preinv.c" fails since that file isn't vendored here. */
#include <flint/ulong_extras.h>
#include "ulong_extras.h"
#include "../msolve/streams.h"
#include "../neogb/gmp_u64.h"

void
_mpz_CRT_ui_precomp(mpz_t out, const mpz_t r1, const mpz_t m1, uint64_t r2,
                    uint64_t m2, mp_limb_t m2inv, const mpz_t m1m2, mp_limb_t c,
                    int sign)
{
  mp_limb_t r1mod, s;
  mpz_t tmp;
  mpz_init(tmp);


  if (mpz_sgn(r1) < 0)
    mpz_add(tmp, r1, m1);
  else
    mpz_set(tmp, r1);

  r1mod = mpz_fdiv_u64(tmp, m2);
  s = n_submod(r2, r1mod, m2);
  s = n_mulmod2_preinv(s, c, m2, m2inv);
  mpz_addmul_u64(tmp, m1, s);

  if (sign)
    {
      mpz_sub(out, tmp, m1m2);
      if (mpz_cmpabs(tmp, out) <= 0)
        mpz_swap(out, tmp);
    }
  else
    {
      mpz_swap(out, tmp);
    }
  mpz_clear(tmp);
}

/**

   returns out s.t. 0 < out < m1 x m2 (if sign=0) or -(m1 x m2) / 2 < out < (m1 x m2) / 2
   out mod m1 = r1 and out mod m2 = r2

 **/

void mpz_CRT_ui(mpz_t out, const mpz_t r1, const mpz_t m1,
                uint64_t r2, uint64_t m2, const mpz_t m1m2,
                mpz_t tmp, int sign)
{
  mp_limb_t c;

  c = mpz_fdiv_u64(m1, m2);
  c = n_invmod(c, m2);

  if (c == 0)
    {
      fprintf(ERRSTREAM, "Exception (fmpz_CRT_ui). m1 not invertible modulo m2.\n");
      exit(1);
    }

  _mpz_CRT_ui_precomp(out, r1, m1, r2, m2, n_preinvert_limb(m2),
                      m1m2, c, sign);

}
