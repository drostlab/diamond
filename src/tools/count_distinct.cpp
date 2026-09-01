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

#include <algorithm>
#include <array>
#include <atomic>
#include <inttypes.h>
#include <iostream>
#include <memory>
#include <vector>
#include "basic/config.h"
#include "basic/sequence.h"
#include "basic/value.h"
#include "data/block/block.h"
#include "data/sequence_file.h"
#include "murmurhash/MurmurHash3.h"
#include "util/algo/hyperloglog.h"
#include "util/io/file.h"
#include "util/log_stream.h"
#include "util/parallel/simple_thread_pool.h"
#include "util/string/string.h"

using std::array;
using std::atomic;
using std::cout;
using std::endl;
using std::runtime_error;
using std::unique_ptr;
using std::vector;

static const int PRECISION = 20;
// Sequences per work packet handed to a hashing thread.
static const BlockId CHUNK_SIZE = 1024;

static uint64_t seq_hash(const Sequence& seq) {
	static const array<uint64_t, 2> SEED = { 0, 0 };
	array<uint64_t, 2> h;
	MurmurHash3_x64_128(seq.data(), (int)seq.length(), (const char*)SEED.data(), h.data());
	return h[0];
}

// Hashes all sequences of the block in parallel, one counter per thread.
static void hash_block(const SequenceSet& seqs, vector<HyperLogLog>& hll, const int threads) {
	const BlockId n = seqs.size();
	atomic<BlockId> next(0);
	SimpleThreadPool pool;
	auto worker = [&](const atomic<bool>& stop, const int thread_id) {
		HyperLogLog& counter = hll[thread_id];
		while (!stop) {
			const BlockId begin = next.fetch_add(CHUNK_SIZE, std::memory_order_relaxed);
			if (begin >= n)
				break;
			const BlockId end = std::min(begin + CHUNK_SIZE, n);
			for (BlockId i = begin; i < end; ++i)
				counter.add_hash(seq_hash(seqs[i]));
		}
		};
	for (int i = 0; i < threads; ++i)
		pool.spawn(worker, i);
	pool.join_all();
}

static HyperLogLog merge_all(const vector<HyperLogLog>& hll) {
	HyperLogLog total(PRECISION);
	for (const HyperLogLog& h : hll)
		total.merge(h);
	return total;
}

void count_distinct() {
	if (config.input_ref_file.size() != 1)
		throw runtime_error("The countdistinct workflow requires exactly one input file (option --in).");
	const int threads = std::max(config.threads_, 1);
	unique_ptr<SequenceFile> file;
	try {
		// Titles are not loaded: only the sequence data is hashed.
		file.reset(SequenceFile::auto_create(config.input_ref_file, SequenceFile::Flags::SEQS, amino_acid_traits));
	}
	catch (const FormatDetectionError& e) {
		throw runtime_error("Error opening input file: " + std::string(e.what()));
	}
	const uint64_t block_size = std::min<uint64_t>(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), SequenceFile::DEFAULT_LOAD_SIZE);
	*message_stream << "Counting distinct sequences. Counter precision=" << PRECISION << " block size=" << block_size << endl;

	vector<HyperLogLog> hll(threads, HyperLogLog(PRECISION));
	int64_t blocks = 0, seqs = 0, letters = 0;
	TaskTimer total_timer;
	for (;;) {
		TaskTimer timer("Loading block");
		unique_ptr<Block> block(file->load_seqs(block_size));
		timer.finish();
		if (block->empty())
			break;
		timer.go("Hashing block");
		hash_block(block->seqs(), hll, threads);
		timer.finish();
		++blocks;
		seqs += block->seqs().size();
		letters += block->seqs().letters();
		*message_stream << "Block=" << blocks << " sequences=" << seqs << " letters=" << letters
			<< " distinct sequences=" << merge_all(hll).estimate() << endl;
	}
	file->close();

	const HyperLogLog total = merge_all(hll);
	if (!config.output_file.empty()) {
		TaskTimer timer("Writing the counter");
		File out(config.output_file, "wb");
		total.serialize(out);
		out.close();
	}

	*message_stream << "Total time = " << total_timer.get() << "s" << endl;
	cout << "Sequences: " << seqs << endl;
	cout << "Distinct sequences (estimate): " << total.estimate() << endl;
}