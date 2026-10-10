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

#ifndef CPU_FEATURES_H
#define CPU_FEATURES_H

/* SIMD intrinsics and CPU feature detection for the vectorized kernels.
 *
 * x86-64: HAVE_AVX2_KERNELS resp. HAVE_AVX512_KERNELS are defined by
 * configure whenever the compiler is able to build AVX2 resp. AVX-512
 * kernels, independently of the CPU of the build machine. Such kernels are
 * marked with TARGET_AVX2 resp. TARGET_AVX512, the remaining code is
 * compiled for the baseline architecture. Callers choose a kernel at runtime
 * via cpu_has_avx2() resp. cpu_has_avx512(), so that the same binary runs on
 * every x86-64 CPU. If msolve is compiled with these extensions enabled
 * globally (e.g. with -march=native) the checks are resolved at compile time.
 *
 * AArch64: NEON is part of the baseline architecture, so the NEON kernels
 * are selected at compile time via __aarch64__. */

#if defined HAVE_AVX2_KERNELS || defined HAVE_AVX512_KERNELS
#include <immintrin.h>
#endif

#ifdef __aarch64__
#include <arm_neon.h>
#endif

#ifdef HAVE_AVX2_KERNELS
#define TARGET_AVX2 __attribute__((target("avx2")))

static inline int cpu_has_avx2(void)
{
#ifdef __AVX2__
    return 1;
#else
    return __builtin_cpu_supports("avx2");
#endif
}
#endif

#ifdef HAVE_AVX512_KERNELS
/* the 8 and 16 bit kernels need AVX-512 BW for 16 bit multiplications */
#define TARGET_AVX512 __attribute__((target("avx512f,avx512bw")))

static inline int cpu_has_avx512(void)
{
#if defined __AVX512F__ && defined __AVX512BW__
    return 1;
#else
    return __builtin_cpu_supports("avx512f")
        && __builtin_cpu_supports("avx512bw");
#endif
}
#endif

#endif
