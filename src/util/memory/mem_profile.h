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

//
// Category-attributed heap profiler.
//
// Answers "which part of the code holds how much memory at the moment the
// process is at its largest", as opposed to "how many bytes did this part
// allocate in total". Every allocation is charged to the innermost enclosing
// MEM_SCOPE, and the charge is reversed when the block is freed, so a category
// always shows *live* (resident-relevant) bytes, never a running total.
//
// Enabled by configuring with -DMEM_PROFILE=ON. Without it every macro below
// compiles to nothing and the global operator new/delete are not replaced.
//
// Usage:
//
//   void foo() {
//       MEM_SCOPE("extend/ungapped-stage");   // RAII, innermost scope wins
//       ...
//   }
//
//   MEM_ADD("align/trace-points", n);         // memory obtained outside the
//   MEM_SUB("align/trace-points", n);         // heap (mmap, VirtualAlloc)
//
//   std::pmr::unsynchronized_pool_resource pool(MEM_POOL("align/pmr-pool"));
//                                             // charge pool chunks to a fixed
//                                             // category rather than to
//                                             // whichever scope grew the pool
//
//   MemProfile::report(std::cerr, "align_queries");
//
// Attribution is exclusive: a nested scope takes the allocation away from its
// parent, so category figures are disjoint and sum to the tracked total.
// Anything allocated with no scope active lands in "(unattributed)", which
// therefore covers the whole rest of the process.
//
// Two numbers are reported per category:
//   peak     - the largest this category was at any point, on its own. Peaks of
//              different categories generally do not coincide, so these do NOT
//              add up to the process peak.
//   at peak  - what the category held at the moment the *tracked total* was at
//              its maximum. These do add up, and this is the column that
//              explains the peak resident set size.
//

#pragma once

#include <cstdint>
#include <iosfwd>
#include "memory_resource.h"

namespace MemProfile {

#ifdef MEM_PROFILE
constexpr bool ENABLED = true;
#else
constexpr bool ENABLED = false;
#endif

constexpr int MAX_CATEGORIES = 128;
constexpr int UNATTRIBUTED = 0;

#ifdef MEM_PROFILE

// Registers a category name (idempotent) and returns its id. The name must have
// static storage duration; it is not copied.
int category(const char* name);

void push(int cat);
void pop();

// Id of the innermost active scope, for objects that have to remember which
// category a later release should be credited to.
int current();

// Accounting for memory that does not come from operator new.
void add(int cat, int64_t bytes);
void sub(int cat, int64_t bytes);

// For raw malloc/realloc wrappers: behaves exactly like operator new/delete,
// charging the innermost active scope and remembering the block so that the
// release is credited back to the same category. Must be paired.
void track(const void* p, int64_t bytes);
void untrack(const void* p);

// Upstream memory resource that charges everything it hands out to `cat`,
// regardless of which scope is active. Use it for pooling resources, whose
// chunks stay resident long after the allocation that triggered them.
std::pmr::memory_resource* upstream(int cat);

void report(std::ostream& out, const char* title);

// Drops peaks and the peak snapshot, keeps live counters. Lets a report cover
// one phase of the run rather than everything since process start.
void reset_peaks();

#else

inline int category(const char*) { return UNATTRIBUTED; }
inline void push(int) {}
inline void pop() {}
inline int current() { return UNATTRIBUTED; }
inline void add(int, int64_t) {}
inline void sub(int, int64_t) {}
inline void track(const void*, int64_t) {}
inline void untrack(const void*) {}
inline std::pmr::memory_resource* upstream(int) { return std::pmr::get_default_resource(); }
inline void report(std::ostream&, const char*) {}
inline void reset_peaks() {}

#endif

struct Scope {
	explicit Scope(int cat) {
		push(cat);
	}
	~Scope() {
		pop();
	}
	Scope(const Scope&) = delete;
	Scope& operator=(const Scope&) = delete;
};

}

#define MEM_PROFILE_CAT2(a, b) a##b
#define MEM_PROFILE_CAT(a, b) MEM_PROFILE_CAT2(a, b)

#ifdef MEM_PROFILE

#define MEM_SCOPE(name) \
	static const int MEM_PROFILE_CAT(mem_profile_id_, __LINE__) = ::MemProfile::category(name); \
	const ::MemProfile::Scope MEM_PROFILE_CAT(mem_profile_scope_, __LINE__)(MEM_PROFILE_CAT(mem_profile_id_, __LINE__))

#define MEM_ADD(name, bytes) do { \
	static const int mem_profile_id_ = ::MemProfile::category(name); \
	::MemProfile::add(mem_profile_id_, (int64_t)(bytes)); } while (0)

#define MEM_SUB(name, bytes) do { \
	static const int mem_profile_id_ = ::MemProfile::category(name); \
	::MemProfile::sub(mem_profile_id_, (int64_t)(bytes)); } while (0)

#define MEM_POOL(name) (::MemProfile::upstream(::MemProfile::category(name)))
#define MEM_TRACK(p, bytes) (::MemProfile::track((p), (int64_t)(bytes)))
#define MEM_UNTRACK(p) (::MemProfile::untrack(p))

#else

#define MEM_SCOPE(name) ((void)0)
#define MEM_ADD(name, bytes) ((void)0)
#define MEM_SUB(name, bytes) ((void)0)
#define MEM_POOL(name) (std::pmr::get_default_resource())
#define MEM_TRACK(p, bytes) ((void)0)
#define MEM_UNTRACK(p) ((void)0)

#endif
