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

#include <stdexcept>
#include <vector>
#include <stdint.h>
#include "simd.h"
#include "util/string/string.h"

#ifdef HAVE_GETAUXVAL
#include <sys/auxv.h>
#endif

using std::string;

#ifdef _WIN32
#define cpuid(info,x)    __cpuidex(info,x,0)
#ifdef __SSE2__
inline uint64_t xgetbv0() {
	return _xgetbv(0);
}
#endif
#else
#ifdef __SSE2__
inline void cpuid(int CPUInfo[4], int InfoType) {
	__asm__ __volatile__(
		"cpuid":
	"=a" (CPUInfo[0]),
		"=b" (CPUInfo[1]),
		"=c" (CPUInfo[2]),
		"=d" (CPUInfo[3]) :
		"a" (InfoType), "c" (0)
		);
}

inline uint64_t xgetbv0() {
	uint32_t eax, edx;
	__asm__ __volatile__("xgetbv" : "=a" (eax), "=d" (edx) : "c" (0));
	return ((uint64_t)edx << 32) | eax;
}
#endif
#endif

namespace SIMD {

int flags = 0;

Arch init_arch() {
#if defined(__aarch64__) | defined(__arm__)

#if defined(__aarch64__)
	flags |= NEON;
#elif defined(HAVE_MFPU_NEON) & defined(HAVE_GETAUXVAL)
	unsigned long hwcaps = getauxval(AT_HWCAP);
	if ((hwcaps & HWCAP_ARM_NEON) != 0)
		flags |= NEON;
#endif

	if ((flags & NEON) != 0)
		return Arch::NEON;

#else

#ifdef __SSE2__
	int info[4];
	cpuid(info, 0);
	int nids = info[0];
	if (nids >= 1) {
		cpuid(info, 1);
	}
	else
		throw std::runtime_error("Incompatible CPU type. Please try to compile the software from source.");

	if ((info[2] & (1 << 9)) != 0)
		flags |= SSSE3;
	if ((info[2] & (1 << 23)) != 0)
		flags |= POPCNT;
	if ((info[2] & (1 << 19)) != 0)
		flags |= SSE4_1;

	/* Wide registers are only usable if the OS has enabled saving them in XCR0.
	   XCR0 itself may only be read when the CPU reports XSAVE (bit 26) and the OS
	   has enabled it via CR4.OSXSAVE (bit 27). */
	static const uint64_t XCR0_YMM = 0x6, XCR0_ZMM = 0xE6;
	const bool osxsave = (info[2] & (1 << 26)) != 0 && (info[2] & (1 << 27)) != 0;
	const uint64_t xcr0 = osxsave ? xgetbv0() : 0;

	if (nids >= 7) {
		cpuid(info, 7);
		if ((info[1] & (1 << 5)) != 0 && (xcr0 & XCR0_YMM) == XCR0_YMM)
			flags |= AVX2;
#ifdef WITH_AVX512
		/* The AVX512 object library is compiled for x86-64-v4, so the full feature
		   set that implies has to be present, not just F and BW. */
		static const uint32_t AVX512F = 1u << 16, AVX512DQ = 1u << 17, AVX512CD = 1u << 28,
			AVX512BW = 1u << 30, AVX512VL = 1u << 31;
		static const uint32_t V4 = AVX512F | AVX512DQ | AVX512CD | AVX512BW | AVX512VL;
		if (((uint32_t)info[1] & V4) == V4 && (xcr0 & XCR0_ZMM) == XCR0_ZMM)
			flags |= AVX512;
#endif
	}
#endif

#ifdef __SSSE3__
	if ((flags & SSSE3) == 0)
		throw std::runtime_error("CPU does not support SSSE3. Please compile the software from source.");
#endif
#ifdef __POPCNT__
	if ((flags & POPCNT) == 0)
		throw std::runtime_error("CPU does not support POPCNT. Please compile the software from source.");
#endif
#ifdef __SSE4_1__
	if ((flags & SSE4_1) == 0)
		throw std::runtime_error("CPU does not support SSE4.1. Please compile the software from source.");
#endif
#ifdef __AVX2__
	if ((flags & AVX2) == 0)
		throw std::runtime_error("CPU does not support AVX2. Please compile the software from source.");
#endif

	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1) && (flags & AVX2) && (flags & AVX512))
		return Arch::AVX512;
	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1) && (flags & AVX2))
		return Arch::AVX2;
	if ((flags & SSSE3) && (flags & POPCNT) && (flags & SSE4_1))
		return Arch::SSE4_1;
#endif
	return Arch::Generic;
}

Arch arch() {
	static Arch a = Arch::None;
	return a == Arch::None ? (a = init_arch()) : a;
}

string features() {
	init_arch();
	std::vector<string> r;
	if (flags & NEON)
		r.push_back("neon");
	if (flags & SSSE3)
		r.push_back("ssse3");
	if (flags & POPCNT)
		r.push_back("popcnt");
	if (flags & SSE4_1)
		r.push_back("sse4.1");
	if (flags & AVX2)
		r.push_back("avx2");
	if (flags & AVX512)
		r.push_back("avx512");
	return r.empty() ? "None" : join(" ", r.begin(), r.end());
}

}