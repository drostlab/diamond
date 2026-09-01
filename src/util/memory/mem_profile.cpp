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

#include "mem_profile.h"

#ifdef MEM_PROFILE

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <new>
#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include "util/system/system.h"

#ifndef _MSC_VER
#include <stdlib.h>
#endif

namespace MemProfile {

namespace {

// ---------------------------------------------------------------------------
// Everything on the allocation path uses malloc/free directly. Anything that
// went through operator new would recurse.
// ---------------------------------------------------------------------------

constexpr int SHARDS = 256;
constexpr int STACK_DEPTH = 32;
constexpr int64_t MIN_SNAPSHOT_STEP = 1 << 20;

struct CatStat {
	std::atomic<int64_t> live;
	std::atomic<int64_t> peak;
	std::atomic<int64_t> live_blocks;
	std::atomic<int64_t> allocs;
	std::atomic<int64_t> total_bytes;
};

// Static storage: zero initialized before any dynamic initialization runs, so
// operator new is safe to call from other translation units static ctors.
CatStat g_cat[MAX_CATEGORIES];
const char* g_name[MAX_CATEGORIES];
std::atomic<int> g_ncat;
std::mutex g_cat_mutex;

std::atomic<int64_t> g_live;
std::atomic<int64_t> g_peak;

std::mutex g_snap_mutex;
std::atomic<int64_t> g_snap_trigger;
int64_t g_snap_bytes[MAX_CATEGORIES];
int64_t g_snap_blocks[MAX_CATEGORIES];
int64_t g_snap_total;
int64_t g_snap_rss;

thread_local int t_stack[STACK_DEPTH];
thread_local int t_depth;
// Suppresses accounting for allocations made by the profiler itself.
thread_local bool t_reentrant;

int current_cat() {
	const int d = std::min(t_depth, STACK_DEPTH);
	return d > 0 ? t_stack[d - 1] : UNATTRIBUTED;
}

// ---------------------------------------------------------------------------
// Live block table: pointer -> (size, category). Needed so a free can be
// credited back to the category that made the allocation, which is what turns
// "bytes allocated" into "bytes held".
// ---------------------------------------------------------------------------

struct Entry {
	const void* p;
	uint64_t size;
	int32_t cat;
};

struct Shard {
	std::mutex mtx;
	Entry* tab;
	uint64_t cap;
	uint64_t n;
};

Shard g_shard[SHARDS];

inline uint64_t hash_ptr(const void* p) {
	uint64_t x = (uint64_t)(uintptr_t)p >> 4;
	x ^= x >> 33;
	x *= 0xff51afd7ed558ccdULL;
	x ^= x >> 33;
	x *= 0xc4ceb9fe1a85ec53ULL;
	x ^= x >> 33;
	return x;
}

// Caller holds s.mtx.
void grow(Shard& s) {
	const uint64_t newcap = s.cap ? s.cap * 2 : 1024;
	Entry* nt = (Entry*)std::calloc((size_t)newcap, sizeof(Entry));
	if (!nt)
		return;
	for (uint64_t i = 0; i < s.cap; ++i) {
		if (!s.tab[i].p)
			continue;
		uint64_t j = (hash_ptr(s.tab[i].p) >> 8) & (newcap - 1);
		while (nt[j].p)
			j = (j + 1) & (newcap - 1);
		nt[j] = s.tab[i];
	}
	std::free(s.tab);
	s.tab = nt;
	s.cap = newcap;
}

void table_insert(const void* p, uint64_t size, int32_t cat) {
	const uint64_t h = hash_ptr(p);
	Shard& s = g_shard[h & (SHARDS - 1)];
	std::lock_guard<std::mutex> lock(s.mtx);
	if ((s.n + 1) * 4 > s.cap * 3)
		grow(s);
	if (!s.cap)
		return;
	const uint64_t mask = s.cap - 1;
	uint64_t j = (h >> 8) & mask;
	while (s.tab[j].p && s.tab[j].p != p)
		j = (j + 1) & mask;
	if (!s.tab[j].p)
		++s.n;
	s.tab[j] = Entry{ p, size, cat };
}

bool table_erase(const void* p, uint64_t& size, int32_t& cat) {
	const uint64_t h = hash_ptr(p);
	Shard& s = g_shard[h & (SHARDS - 1)];
	std::lock_guard<std::mutex> lock(s.mtx);
	if (!s.cap)
		return false;
	const uint64_t mask = s.cap - 1;
	uint64_t j = (h >> 8) & mask;
	while (s.tab[j].p && s.tab[j].p != p)
		j = (j + 1) & mask;
	if (!s.tab[j].p)
		return false;
	size = s.tab[j].size;
	cat = s.tab[j].cat;
	// Backward shift deletion: keeps probe sequences intact without tombstones.
	uint64_t i = j, k = (j + 1) & mask;
	s.tab[i].p = nullptr;
	while (s.tab[k].p) {
		const uint64_t ideal = (hash_ptr(s.tab[k].p) >> 8) & mask;
		if (((k - ideal) & mask) >= ((k - i) & mask)) {
			s.tab[i] = s.tab[k];
			s.tab[k].p = nullptr;
			i = k;
		}
		k = (k + 1) & mask;
	}
	--s.n;
	return true;
}

int64_t table_bytes() {
	int64_t n = 0;
	for (int i = 0; i < SHARDS; ++i) {
		std::lock_guard<std::mutex> lock(g_shard[i].mtx);
		n += (int64_t)g_shard[i].cap * (int64_t)sizeof(Entry);
	}
	return n;
}

// ---------------------------------------------------------------------------
// Accounting
// ---------------------------------------------------------------------------

void take_snapshot(int64_t total) {
	t_reentrant = true;
	{
		std::lock_guard<std::mutex> lock(g_snap_mutex);
		if (total > g_snap_total) {
			const int n = g_ncat.load(std::memory_order_acquire);
			for (int i = 0; i < n; ++i) {
				g_snap_bytes[i] = g_cat[i].live.load(std::memory_order_relaxed);
				g_snap_blocks[i] = g_cat[i].live_blocks.load(std::memory_order_relaxed);
			}
			g_snap_total = total;
			g_snap_rss = (int64_t)getCurrentRSS();
			g_snap_trigger.store(total + std::max(MIN_SNAPSHOT_STEP, total / 64), std::memory_order_relaxed);
		}
	}
	t_reentrant = false;
}

void account(int cat, int64_t bytes, int64_t blocks) {
	CatStat& c = g_cat[cat];
	const int64_t live = c.live.fetch_add(bytes, std::memory_order_relaxed) + bytes;
	c.live_blocks.fetch_add(blocks, std::memory_order_relaxed);
	const int64_t total = g_live.fetch_add(bytes, std::memory_order_relaxed) + bytes;
	if (bytes <= 0)
		return;
	c.allocs.fetch_add(1, std::memory_order_relaxed);
	c.total_bytes.fetch_add(bytes, std::memory_order_relaxed);
	int64_t p = c.peak.load(std::memory_order_relaxed);
	while (live > p && !c.peak.compare_exchange_weak(p, live, std::memory_order_relaxed))
		;
	p = g_peak.load(std::memory_order_relaxed);
	while (total > p && !g_peak.compare_exchange_weak(p, total, std::memory_order_relaxed))
		;
	// The per category breakdown is only re-snapshotted once the total has moved
	// a meaningful step past the last snapshot, otherwise walking all categories
	// would dominate the allocation path. g_peak itself stays exact.
	if (total >= g_snap_trigger.load(std::memory_order_relaxed) && !t_reentrant)
		take_snapshot(total);
}

void note_alloc(void* p, size_t size) {
	if (!p || t_reentrant)
		return;
	const int cat = current_cat();
	t_reentrant = true;
	table_insert(p, size, cat);
	t_reentrant = false;
	account(cat, (int64_t)size, 1);
}

void note_free(void* p) {
	if (!p || t_reentrant)
		return;
	uint64_t size;
	int32_t cat;
	t_reentrant = true;
	const bool found = table_erase(p, size, cat);
	t_reentrant = false;
	if (found)
		account(cat, -(int64_t)size, -1);
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

std::string format_bytes(int64_t n) {
	char buf[64];
	const double x = (double)n;
	if (n < 1024)
		std::snprintf(buf, sizeof(buf), "%lld B", (long long)n);
	else if (n < (1LL << 20))
		std::snprintf(buf, sizeof(buf), "%.1f KB", x / 1024.0);
	else if (n < (1LL << 30))
		std::snprintf(buf, sizeof(buf), "%.1f MB", x / 1048576.0);
	else
		std::snprintf(buf, sizeof(buf), "%.2f GB", x / 1073741824.0);
	return std::string(buf);
}

std::string format_count(int64_t n) {
	char buf[32];
	std::snprintf(buf, sizeof(buf), "%lld", (long long)n);
	std::string s(buf);
	for (int i = (int)s.size() - 3; i > 0; i -= 3)
		s.insert(i, ",");
	return s;
}

std::string pad_left(const std::string& s, size_t w) {
	return s.size() >= w ? s : std::string(w - s.size(), ' ') + s;
}

std::string pad_right(const std::string& s, size_t w) {
	return s.size() >= w ? s : s + std::string(w - s.size(), ' ');
}

struct Row {
	std::string name;
	int64_t peak, at_peak, blocks_at_peak, allocs, live, live_blocks;
};

std::string group_of(const std::string& name) {
	const size_t slash = name.find('/');
	return slash == std::string::npos ? name : name.substr(0, slash);
}

void print_row(std::ostream& out, const std::string& label, const Row& r, int64_t denom) {
	char pct[32];
	std::snprintf(pct, sizeof(pct), "%.1f%%", denom > 0 ? 100.0 * (double)r.at_peak / (double)denom : 0.0);
	out << "  " << pad_right(label, 36)
		<< pad_left(format_bytes(r.peak), 11)
		<< pad_left(format_bytes(r.at_peak), 12)
		<< pad_left(pct, 8)
		<< pad_left(format_count(r.blocks_at_peak), 14)
		<< pad_left(format_count(r.allocs), 16)
		<< '\n';
}

// Charges everything it hands out to a fixed category. Pooling resources hold
// on to their chunks across many scopes, so charging a chunk to whichever scope
// happened to exhaust the pool would be misleading.
class PinnedResource : public std::pmr::memory_resource {
public:
	explicit PinnedResource(int cat) :
		cat_(cat)
	{}
protected:
	void* do_allocate(size_t bytes, size_t alignment) override {
		const Scope s(cat_);
		return std::pmr::new_delete_resource()->allocate(bytes, alignment);
	}
	void do_deallocate(void* p, size_t bytes, size_t alignment) override {
		std::pmr::new_delete_resource()->deallocate(p, bytes, alignment);
	}
	bool do_is_equal(const std::pmr::memory_resource& other) const noexcept override {
		return this == &other;
	}
private:
	const int cat_;
};

PinnedResource* g_pinned[MAX_CATEGORIES];
std::mutex g_pinned_mutex;

}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

int category(const char* name) {
	int n = g_ncat.load(std::memory_order_acquire);
	for (int i = 0; i < n; ++i)
		if (std::strcmp(g_name[i], name) == 0)
			return i;
	std::lock_guard<std::mutex> lock(g_cat_mutex);
	n = g_ncat.load(std::memory_order_relaxed);
	if (n == 0) {
		g_name[0] = "(unattributed)";
		n = 1;
		g_ncat.store(1, std::memory_order_release);
	}
	for (int i = 0; i < n; ++i)
		if (std::strcmp(g_name[i], name) == 0)
			return i;
	if (n >= MAX_CATEGORIES)
		return UNATTRIBUTED;
	g_name[n] = name;
	g_ncat.store(n + 1, std::memory_order_release);
	return n;
}

void push(int cat) {
	if (t_depth < STACK_DEPTH)
		t_stack[t_depth] = cat;
	++t_depth;
}

void pop() {
	--t_depth;
}

int current() {
	return current_cat();
}

// Regions booked by hand are not malloc blocks, so they do not take part in
// the block count that the allocator overhead estimate is based on.
void add(int cat, int64_t bytes) {
	account(cat, bytes, 0);
}

void sub(int cat, int64_t bytes) {
	account(cat, -bytes, 0);
}

void track(const void* p, int64_t bytes) {
	if (!p || t_reentrant)
		return;
	const int cat = current_cat();
	t_reentrant = true;
	table_insert(p, (uint64_t)bytes, cat);
	t_reentrant = false;
	account(cat, bytes, 1);
}

void untrack(const void* p) {
	note_free(const_cast<void*>(p));
}

std::pmr::memory_resource* upstream(int cat) {
	if (cat < 0 || cat >= MAX_CATEGORIES)
		return std::pmr::new_delete_resource();
	std::lock_guard<std::mutex> lock(g_pinned_mutex);
	if (!g_pinned[cat]) {
		void* p = std::malloc(sizeof(PinnedResource));
		if (!p)
			return std::pmr::new_delete_resource();
		g_pinned[cat] = new (p) PinnedResource(cat);
	}
	return g_pinned[cat];
}

void reset_peaks() {
	std::lock_guard<std::mutex> lock(g_snap_mutex);
	const int n = g_ncat.load(std::memory_order_acquire);
	for (int i = 0; i < n; ++i) {
		g_cat[i].peak.store(g_cat[i].live.load(std::memory_order_relaxed), std::memory_order_relaxed);
		g_snap_bytes[i] = 0;
		g_snap_blocks[i] = 0;
	}
	g_peak.store(g_live.load(std::memory_order_relaxed), std::memory_order_relaxed);
	g_snap_total = 0;
	g_snap_rss = 0;
	g_snap_trigger.store(MIN_SNAPSHOT_STEP, std::memory_order_relaxed);
}

void report(std::ostream& out, const char* title) {
	std::vector<Row> rows;
	int64_t snap_total, snap_rss;
	{
		std::lock_guard<std::mutex> lock(g_snap_mutex);
		const int n = g_ncat.load(std::memory_order_acquire);
		rows.reserve(n);
		for (int i = 0; i < n; ++i)
			rows.push_back(Row{ std::string(g_name[i]),
				g_cat[i].peak.load(std::memory_order_relaxed),
				g_snap_bytes[i],
				g_snap_blocks[i],
				g_cat[i].allocs.load(std::memory_order_relaxed),
				g_cat[i].live.load(std::memory_order_relaxed),
				g_cat[i].live_blocks.load(std::memory_order_relaxed) });
		snap_total = g_snap_total;
		snap_rss = g_snap_rss;
	}

	const int64_t peak = g_peak.load(std::memory_order_relaxed);
	const int64_t live = g_live.load(std::memory_order_relaxed);
	int64_t live_blocks = 0, snap_blocks = 0;
	for (const Row& r : rows) {
		live_blocks += r.live_blocks;
		snap_blocks += r.blocks_at_peak;
	}

	out << '\n' << "Memory profile";
	if (title && *title)
		out << " [" << title << ']';
	out << '\n';
	out << "  peak of tracked heap        " << pad_left(format_bytes(peak), 11) << '\n';
	out << "  process RSS at that point   " << pad_left(format_bytes(snap_rss), 11) << '\n';
	out << "  process peak RSS            " << pad_left(format_bytes((int64_t)getPeakRSS()), 11) << '\n';
	out << "  process RSS now             " << pad_left(format_bytes((int64_t)getCurrentRSS()), 11) << '\n';
	out << "  tracked heap live now       " << pad_left(format_bytes(live), 11)
		<< "  in " << format_count(live_blocks) << " blocks" << '\n';
	out << "  allocator overhead at peak  " << pad_left(format_bytes(snap_blocks * 16), 11)
		<< "  (est. 16 B over " << format_count(snap_blocks) << " live blocks)" << '\n';
	out << "  profiler table              " << pad_left(format_bytes(table_bytes()), 11) << '\n';
	if (snap_rss > 0)
		out << "  outside the tracked heap    " << pad_left(format_bytes(snap_rss - snap_total), 11)
			<< "  (statics, thread stacks, mapped files, allocator retention)" << '\n';
	out << '\n';
	out << "  " << pad_right("category", 36) << pad_left("peak", 11) << pad_left("at peak", 12)
		<< pad_left("share", 8) << pad_left("blocks@peak", 14) << pad_left("allocs", 16) << '\n';
	out << "  " << std::string(97, '-') << '\n';

	// Group by the part of the name before the first slash, so a feature made of
	// several scopes can be read off as a single number.
	std::vector<std::string> groups;
	for (const Row& r : rows) {
		const std::string g = group_of(r.name);
		if (std::find(groups.begin(), groups.end(), g) == groups.end())
			groups.push_back(g);
	}
	std::vector<std::pair<int64_t, std::string>> ordered;
	for (const std::string& g : groups) {
		int64_t s = 0;
		for (const Row& r : rows)
			if (group_of(r.name) == g)
				s += r.at_peak;
		ordered.emplace_back(s, g);
	}
	std::sort(ordered.begin(), ordered.end(), [](const std::pair<int64_t, std::string>& a, const std::pair<int64_t, std::string>& b) {
		return a.first > b.first;
		});

	for (const std::pair<int64_t, std::string>& g : ordered) {
		std::vector<Row> members;
		for (const Row& r : rows)
			if (group_of(r.name) == g.second)
				members.push_back(r);
		std::sort(members.begin(), members.end(), [](const Row& a, const Row& b) { return a.at_peak > b.at_peak; });
		if (members.size() > 1) {
			Row sum{ std::string(), 0, 0, 0, 0, 0, 0 };
			for (const Row& r : members) {
				sum.peak += r.peak;
				sum.at_peak += r.at_peak;
				sum.blocks_at_peak += r.blocks_at_peak;
				sum.allocs += r.allocs;
			}
			print_row(out, g.second + "/*", sum, snap_total);
			for (const Row& r : members)
				print_row(out, "  " + r.name, r, snap_total);
		}
		else if (!members.empty())
			print_row(out, members.front().name, members.front(), snap_total);
	}

	out << "  " << std::string(97, '-') << '\n';
	Row total{ std::string(), 0, 0, 0, 0, 0, 0 };
	for (const Row& r : rows) {
		total.at_peak += r.at_peak;
		total.blocks_at_peak += r.blocks_at_peak;
		total.allocs += r.allocs;
	}
	total.peak = peak;
	print_row(out, "total", total, snap_total);
	out << '\n'
		<< "  Peaks cover the window since they were last reset (entry to the phase this" << '\n'
		<< "  report is named after); live figures cover everything resident, including" << '\n'
		<< "  memory held over from earlier phases." << '\n'
		<< '\n'
		<< "  peak    = largest this category ever was on its own. Peaks of different" << '\n'
		<< "            categories need not coincide, so they do not add up." << '\n'
		<< "  at peak = what the category held at the moment the tracked heap was at" << '\n'
		<< "            its maximum. This column is the breakdown of the peak RSS." << '\n'
		<< "  Attribution is exclusive: a nested scope takes bytes from its parent," << '\n'
		<< "  and anything allocated with no scope active lands in (unattributed)." << '\n';
	out.flush();
}

}

// ---------------------------------------------------------------------------
// Global operator new / delete replacement
// ---------------------------------------------------------------------------

namespace {

void* profiled_malloc(size_t size) {
	void* p = std::malloc(size ? size : 1);
	if (p)
		MemProfile::note_alloc(p, size);
	return p;
}

void* profiled_malloc_aligned(size_t size, size_t alignment) {
	if (alignment < sizeof(void*))
		alignment = sizeof(void*);
	void* p;
#ifdef _MSC_VER
	p = _aligned_malloc(size ? size : 1, alignment);
#else
	if (posix_memalign(&p, alignment, size ? size : 1) != 0)
		p = nullptr;
#endif
	if (p)
		MemProfile::note_alloc(p, size);
	return p;
}

void profiled_free(void* p) {
	MemProfile::note_free(p);
	std::free(p);
}

void profiled_free_aligned(void* p) {
	MemProfile::note_free(p);
#ifdef _MSC_VER
	_aligned_free(p);
#else
	std::free(p);
#endif
}

}

void* operator new(size_t size) {
	void* p = profiled_malloc(size);
	if (!p)
		throw std::bad_alloc();
	return p;
}

void* operator new[](size_t size) {
	void* p = profiled_malloc(size);
	if (!p)
		throw std::bad_alloc();
	return p;
}

void* operator new(size_t size, const std::nothrow_t&) noexcept {
	return profiled_malloc(size);
}

void* operator new[](size_t size, const std::nothrow_t&) noexcept {
	return profiled_malloc(size);
}

void* operator new(size_t size, std::align_val_t al) {
	void* p = profiled_malloc_aligned(size, (size_t)al);
	if (!p)
		throw std::bad_alloc();
	return p;
}

void* operator new[](size_t size, std::align_val_t al) {
	void* p = profiled_malloc_aligned(size, (size_t)al);
	if (!p)
		throw std::bad_alloc();
	return p;
}

void* operator new(size_t size, std::align_val_t al, const std::nothrow_t&) noexcept {
	return profiled_malloc_aligned(size, (size_t)al);
}

void* operator new[](size_t size, std::align_val_t al, const std::nothrow_t&) noexcept {
	return profiled_malloc_aligned(size, (size_t)al);
}

void operator delete(void* p) noexcept { profiled_free(p); }
void operator delete[](void* p) noexcept { profiled_free(p); }
void operator delete(void* p, size_t) noexcept { profiled_free(p); }
void operator delete[](void* p, size_t) noexcept { profiled_free(p); }
void operator delete(void* p, const std::nothrow_t&) noexcept { profiled_free(p); }
void operator delete[](void* p, const std::nothrow_t&) noexcept { profiled_free(p); }
void operator delete(void* p, std::align_val_t) noexcept { profiled_free_aligned(p); }
void operator delete[](void* p, std::align_val_t) noexcept { profiled_free_aligned(p); }
void operator delete(void* p, size_t, std::align_val_t) noexcept { profiled_free_aligned(p); }
void operator delete[](void* p, size_t, std::align_val_t) noexcept { profiled_free_aligned(p); }
void operator delete(void* p, std::align_val_t, const std::nothrow_t&) noexcept { profiled_free_aligned(p); }
void operator delete[](void* p, std::align_val_t, const std::nothrow_t&) noexcept { profiled_free_aligned(p); }

#endif
