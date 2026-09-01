/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink
Arm NEON port contributed by Martin Larralde <martin.larralde@embl.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#ifdef WITH_AVX2
#define HAVE_AVX2(x) x
#else
#define HAVE_AVX2(x)
#endif

#ifdef WITH_AVX512
#define HAVE_AVX512(x) x
#else
#define HAVE_AVX512(x)
#endif

#ifdef WITH_SSE4_1
#define HAVE_SSE4_1(x) x
#else
#define HAVE_SSE4_1(x)
#endif

#ifdef WITH_NEON
#define HAVE_NEON(x) x
#else
#define HAVE_NEON(x)
#endif

#if defined(WITH_NEON) | defined(WITH_AVX512) | defined(WITH_AVX2) | defined(WITH_SSE4_1)
#define HAVE_SIMD(x) x
#else
#define HAVE_SIMD(x)
#endif

#if ARCH_ID == 0

/* The parameter and argument lists are passed as parenthesized blobs so that the
   commas separating them are hidden from the preprocessor. */

#define SIMD_DISPATCH_DECL(ret, name, params)\
HAVE_SSE4_1(namespace ARCH_SSE4_1 { ret name params; })\
HAVE_AVX2(namespace ARCH_AVX2 { ret name params; })\
HAVE_AVX512(namespace ARCH_AVX512 { ret name params; })\
HAVE_NEON(namespace ARCH_NEON { ret name params; })

/* `return f(...);` is well formed in a function returning void as long as f does,
   so this also serves the void variants below. */
#define SIMD_DISPATCH_BODY(name, args)\
HAVE_SIMD(switch(::SIMD::arch()) {)\
HAVE_NEON(case ::SIMD::Arch::NEON: return ARCH_NEON::name args;)\
HAVE_AVX512(case ::SIMD::Arch::AVX512: return ARCH_AVX512::name args;)\
HAVE_AVX2(case ::SIMD::Arch::AVX2: return ARCH_AVX2::name args;)\
HAVE_SSE4_1(case ::SIMD::Arch::SSE4_1: return ARCH_SSE4_1::name args;)\
HAVE_SIMD(default:)\
return ARCH_GENERIC::name args;\
HAVE_SIMD(})

#define SIMD_DISPATCH_FN(ret, name, params, args)\
SIMD_DISPATCH_DECL(ret, name, params)\
ret name params {\
SIMD_DISPATCH_BODY(name, args)\
}

#define DISPATCH_0V(name)\
SIMD_DISPATCH_FN(void, name, (), ())

#define DISPATCH_1V(name, t1, n1)\
SIMD_DISPATCH_FN(void, name, (t1 n1), (n1))

#define DISPATCH_1(ret, name, t1, n1)\
SIMD_DISPATCH_FN(ret, name, (t1 n1), (n1))

#define DISPATCH_2(ret, name, t1, n1, t2, n2)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2), (n1, n2))

#define DISPATCH_3(ret, name, t1, n1, t2, n2, t3, n3)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3), (n1, n2, n3))

#define DISPATCH_3V(name, t1, n1, t2, n2, t3, n3)\
SIMD_DISPATCH_FN(void, name, (t1 n1, t2 n2, t3 n3), (n1, n2, n3))

#define DISPATCH_4(ret, name, t1, n1, t2, n2, t3, n3, t4, n4)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3, t4 n4), (n1, n2, n3, n4))

#define DISPATCH_5(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5), (n1, n2, n3, n4, n5))

#define DISPATCH_6(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6), (n1, n2, n3, n4, n5, n6))

#define DISPATCH_6V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)\
SIMD_DISPATCH_FN(void, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6), (n1, n2, n3, n4, n5, n6))

#define DISPATCH_7(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7), (n1, n2, n3, n4, n5, n6, n7))

#define DISPATCH_7V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)\
SIMD_DISPATCH_FN(void, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7), (n1, n2, n3, n4, n5, n6, n7))

#define DISPATCH_8(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7, t8, n8)\
SIMD_DISPATCH_FN(ret, name, (t1 n1, t2 n2, t3 n3, t4 n4, t5 n5, t6 n6, t7 n7, t8 n8), (n1, n2, n3, n4, n5, n6, n7, n8))

#else

#define DISPATCH_0V(name)
#define DISPATCH_1V(name, t1, n1)
#define DISPATCH_1(ret, name, t1, n1)
#define DISPATCH_2(ret, name, t1, n1, t2, n2)
#define DISPATCH_3(ret, name, t1, n1, t2, n2, t3, n3)
#define DISPATCH_4(ret, name, t1, n1, t2, n2, t3, n3, t4, n4)
#define DISPATCH_5(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5)
#define DISPATCH_6(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)
#define DISPATCH_7(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)
#define DISPATCH_3V(name, t1, n1, t2, n2, t3, n3)
#define DISPATCH_6V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6)
#define DISPATCH_7V(name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7)
#define DISPATCH_8(ret, name, t1, n1, t2, n2, t3, n3, t4, n4, t5, n5, t6, n6, t7, n7, t8, n8)

#endif
