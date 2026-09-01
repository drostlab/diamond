/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
#include <stdint.h>
#ifdef _MSC_VER
#include <intrin.h>
#endif

// High 64 bits of the product of two 64 bit values.
static inline uint64_t mul_hi(const uint64_t a, const uint64_t b)
{
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_ARM64))
	return __umulh(a, b);
#elif defined(__SIZEOF_INT128__)
	return (uint64_t)(((__uint128_t)a * (__uint128_t)b) >> 64);
#else
	const uint64_t a0 = a & 0xffffffffull, a1 = a >> 32, b0 = b & 0xffffffffull, b1 = b >> 32;
	const uint64_t t = a0 * b0, u = a1 * b0 + (t >> 32), v = a0 * b1 + (u & 0xffffffffull);
	return a1 * b1 + (u >> 32) + (v >> 32);
#endif
}

static inline unsigned popcount32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}

static inline unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
	return (unsigned)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
}

static inline int ctz(uint32_t x)
{
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanForward(&i, x);
	return i;
#else
	return __builtin_ctz(x);
#endif
}


static inline int ctz(uint64_t x)
{
#ifdef _MSC_VER
	if (x)
		return (int)__popcnt64((x ^ (x - 1)) >> 1);
	else
		return CHAR_BIT * sizeof(x);
#else
	return __builtin_ctzll(x);
#endif
}

static inline int clz(uint64_t x) {
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanReverse64(&i, x);
	return 63 - (int)i;
#else
	return __builtin_clzll((unsigned long long)x);
#endif
}

static inline int clz(uint32_t x) {
#ifdef _MSC_VER
	unsigned long i;
	unsigned char c = _BitScanReverse(&i, x);
	return 31 - (int)i;
#else
	return __builtin_clz((unsigned int)x);
#endif
}