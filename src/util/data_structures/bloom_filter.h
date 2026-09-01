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
#include <algorithm>
#include <atomic>
#include <cmath>
#include <memory>
#include <stdint.h>
#include <thread>
#include <vector>
#include "../hash_function.h"
#include "../intrin.h"

/* Blocked Bloom filter over 64 bit keys. All bits of a key lie within the same block
   of 512 bits, so a query reads one or two cache lines instead of scattering over the
   whole filter. Insertions are lock free and may run from any number of threads;
   queries are not synchronized against them, the filter is meant to be filled first
   and queried afterwards.

   A query for a key that was inserted always returns true, a query for any other key
   returns true with the probability given by false_positive_rate(). */

struct BloomFilter {

	static constexpr int BLOCK_BITS = 512;
	static constexpr int BLOCK_WORDS = BLOCK_BITS / 64;
	static constexpr int HASHES = 3;
	// 8 bits per key put the false positive rate at ~3% for 3 hash functions.
	static constexpr int BITS_PER_KEY = 8;
	// Number of bits of the hash value that one in block bit index is taken from. The
	// index of the block is taken from the high bits of the same hash value, which
	// HASHES * BIT_INDEX_BITS has to stay clear of.
	static constexpr int BIT_INDEX_BITS = 9;

	/* expected_keys sizes the filter. Exceeding it does not break the filter, it only
	   raises the false positive rate. */
	BloomFilter(int64_t expected_keys, int threads = 1, int bits_per_key = BITS_PER_KEY) :
		blocks_((uint64_t)std::max<int64_t>(1, (expected_keys * bits_per_key + BLOCK_BITS - 1) / BLOCK_BITS)),
		data_(new std::atomic<uint64_t>[blocks_ * BLOCK_WORDS])
	{
		clear(std::max(threads, 1));
	}

	void insert(const uint64_t key) {
		const uint64_t h = MurmurHash()(key);
		std::atomic<uint64_t>* const block = data_.get() + block_index(h) * BLOCK_WORDS;
		for (int i = 0; i < HASHES; ++i) {
			const uint64_t b = bit_index(h, i);
			block[b >> 6].fetch_or(uint64_t(1) << (b & 63), std::memory_order_relaxed);
		}
	}

	bool contains(const uint64_t key) const {
		const uint64_t h = MurmurHash()(key);
		const std::atomic<uint64_t>* const block = data_.get() + block_index(h) * BLOCK_WORDS;
		for (int i = 0; i < HASHES; ++i) {
			const uint64_t b = bit_index(h, i);
			if ((block[b >> 6].load(std::memory_order_relaxed) & (uint64_t(1) << (b & 63))) == 0)
				return false;
		}
		return true;
	}

	// Size of the filter in bytes.
	int64_t size() const {
		return (int64_t)blocks_ * BLOCK_WORDS * (int64_t)sizeof(uint64_t);
	}

	// Fraction of the bits of the filter that are set.
	double load() const {
		uint64_t n = 0;
		const uint64_t words = blocks_ * BLOCK_WORDS;
		for (uint64_t i = 0; i < words; ++i)
			n += popcount64(data_[i].load(std::memory_order_relaxed));
		return (double)n / ((double)blocks_ * BLOCK_BITS);
	}

	// Probability that a key that was never inserted tests positive, estimated from the
	// fraction of bits that are actually set, as returned by load().
	static double false_positive_rate(const double load) {
		return std::pow(load, HASHES);
	}

private:

	void clear(const int threads) {
		const uint64_t n = blocks_ * BLOCK_WORDS;
		std::vector<std::thread> workers;
		for (int i = 0; i < threads; ++i)
			workers.emplace_back([this, n, threads, i] {
				const uint64_t begin = n / threads * i, end = i == threads - 1 ? n : n / threads * (i + 1);
				for (uint64_t j = begin; j < end; ++j)
					data_[j].store(0, std::memory_order_relaxed);
				});
		for (auto& t : workers)
			t.join();
	}

	// Maps a hash value to [0, blocks_) without a division (Lemire's fastrange), which
	// uses its high bits.
	uint64_t block_index(const uint64_t h) const {
		return mul_hi(h, blocks_);
	}

	static uint64_t bit_index(const uint64_t h, const int i) {
		return (h >> (i * BIT_INDEX_BITS)) & (BLOCK_BITS - 1);
	}

	const uint64_t blocks_;
	const std::unique_ptr<std::atomic<uint64_t>[]> data_;

};
