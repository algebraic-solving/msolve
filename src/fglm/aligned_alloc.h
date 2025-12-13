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

#ifndef ALIGNED_ALLOC_HEADER_H
#define ALIGNED_ALLOC_HEADER_H

#ifdef _WIN32

#include <errno.h>
#include <malloc.h>

static inline int posix_memalign(void **__memptr, size_t __alignment, size_t __size)
{
    void *p = _aligned_malloc(__size, __alignment);
    if (!p)
    {
        return ENOMEM;
    }
    *__memptr = p;
    return 0;
}
#endif

static inline void posix_memalign_free(void *__p)
{
#ifdef _WIN32
    _aligned_free(__p);
#else
    free(__p);
#endif
}

#endif /* ALIGNED_ALLOC_HEADER_H */
