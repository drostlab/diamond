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

#include <inttypes.h>
#include <algorithm>
#include <memory>
#include "multinode.h"
#include "data/block/block.h"
#include "data/fasta/fasta_file.h"
#include "search/lin_index/lin_index.h"
#include "search/search.h"
#include "search/seed_array/enum_seeds.h"
#include "search/seed_complexity.h"
#include "util/log_stream.h"
#include "util/ptr_vector.h"
#include "multinode.h"
#include "util/algo/hyperloglog.h"

using std::endl;
using std::string;
using std::ofstream;
using std::string;
using std::unique_ptr;
using std::vector;
using std::pair;

static const int MIN_PRECISION = 8;
static const int MAX_PRECISION = 10;


static int hll_precision(const uint64_t max_distinct) {
	int p = MIN_PRECISION;
	while (p < MAX_PRECISION && ((uint64_t)1 << p) < max_distinct)
		++p;
	return p;
}

static uint64_t seed_hash(const uint64_t key) {
	uint64_t x = key + 0x9e3779b97f4a7c15ull;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
	return x ^ (x >> 31);
}

struct CountCallback {

	CountCallback(const SequenceSet& seqs, const double seed_complexity_cut, int shapes, const int precision) :
		seqs(seqs),
		seed_complexity_cut(seed_complexity_cut),
		hll(shapes, HyperLogLog(precision))
	{
	}

	bool operator()(const uint64_t key, const uint64_t pos, uint32_t, const uint64_t shape_id) {
		/*if (!Search::seed_is_complex(seqs.data(pos), shapes[(int)shape_id], seed_complexity_cut))
			return true;*/
		hll[shape_id].add_hash(seed_hash(key));
		return true;
	}

	void finish() {}

	const SequenceSet& seqs;
	const double seed_complexity_cut;
	vector<HyperLogLog> hll;

};

static vector<HyperLogLog> count_block(Block& block, const Search::Config& cfg, const int threads, const int precision) {
	const SequenceSet& seqs = block.seqs();
	const auto partition = seqs.partition(threads);
	// The counters size the hash tables of the seed index, so the seeds have to be
	// enumerated exactly as build_lin_index does it.
	const EnumCfg enum_cfg{ &partition, 0, shapes.count(), Search::LinIndex::encoding(cfg), nullptr, true, false, cfg.seed_complexity_cut,
		cfg.soft_masking, cfg.minimizer_window, false, false, cfg.sketch_size, nullptr };
	PtrVector<CountCallback> v;
	for (int i = 0; i < threads; ++i)
		v.push_back(new CountCallback(seqs, cfg.seed_complexity_cut, shapes.count(), precision));
	enum_seeds(block, v, &no_filter, enum_cfg);
	vector<HyperLogLog> hll(shapes.count(), HyperLogLog(precision));
	for (int i = 0; i < threads; ++i)
		for (int j = 0; j < shapes.count(); ++j)
			hll[j].merge(v[i].hll[j]);
	return hll;
}

pair<vector<vector<HyperLogLog>>, vector<uint64_t>> count_distinct_seeds(Job& job, const VolumedFile& minichunks, const uint64_t letter_count) {
	unique_ptr<vector<BitVector>> no_seed_hits;
	Search::Config cfg(no_seed_hits);
	Search::setup_search(config.sensitivity, cfg);
	const int threads = std::max(config.threads_, 1);
	const int64_t n = (int64_t)minichunks.size();
	const int precision = hll_precision(letter_count);
	job.log("Counting distinct seeds. Counter precision=%i", precision);

	vector<vector<HyperLogLog>> counters;
	vector<uint64_t> letters;
	counters.reserve(n);
	letters.reserve(n);
	TaskTimer total_timer;
	for (int64_t v = 0; v < n; ++v) {
		TaskTimer timer("Loading minichunk");
		unique_ptr<SequenceFile> file;
		try {
			file.reset(new FastaFile({ minichunks[v].path }, SequenceFile::Flags::SEQS, amino_acid_traits));
		}
		catch (const EmptyFileError& e) {
			counters.push_back(vector<HyperLogLog>(shapes.count(), HyperLogLog(precision)));
			letters.push_back(0);
			continue;
		}
		unique_ptr<Block> block(file->load_seqs(INT64_MAX));
		file->close();
		timer.finish();
		timer.go("Counting distinct seeds");
		counters.push_back(count_block(*block, cfg, threads, precision));
		timer.finish();
		std::ostringstream ss;
		for (size_t i = 0; i < counters.back().size(); ++i) {
			ss << counters.back()[i].estimate() << (i + 1 < counters.back().size() ? "," : "");
		}
		letters.push_back(block->seqs().letters());
		job.log("Counting distinct seeds. Minichunk=%lli/%lli letters=%" PRIu64 " distinct seeds=", v + 1, n, (uint64_t)block->seqs().letters());
		job.log_raw(ss.str());
	}
	return std::make_pair(counters, letters);
}