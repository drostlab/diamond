/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink
With contributions from Patrick Ettenhuber Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)

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
#include <string>
#include <vector>
#include "basic/value.h"
#include "data/sequence_file.h"
#include "util/data_structures/flat_array.h"
#include "basic/match.h"
#include "dp/flags.h"
#include "output/output_format.h"

class ClusteringAlgorithm {
public:
	virtual void run() = 0;
	virtual std::string get_description() = 0;
	virtual ~ClusteringAlgorithm(){};
};

namespace Cluster {

const double DEFAULT_MEMBER_COVER = 80.0;
extern const char* const HEADER_LINE;

struct CentroidSorted {};

void realign(const FlatArray<OId>& clusters, const std::vector<OId>& centroids, SequenceFile& db, std::function<void(const HspContext&)>& callback, HspValues hsp_values);
void realign(const std::vector<OId>& clustering, SequenceFile& db, std::function<void(const HspContext&)>& callback, HspValues hsp_values);
template<typename Int>
std::pair<FlatArray<Int>, std::vector<Int>> read(const std::string& file_name, const SequenceFile& db, CentroidSorted);
template<typename Int>
std::vector<Int> read(const std::string& file_name, const SequenceFile& db);
template<typename Int>
std::vector<Int> member2centroid_mapping(const FlatArray<Int>& clusters, const std::vector<Int>& centroids);
template<typename Int>
std::pair<FlatArray<Int>, std::vector<Int>> cluster_sorted(const std::vector<Int>& mapping);
template<typename Int>
std::pair<std::vector<Int>, std::vector<Int>> split(const std::vector<Int>& mapping);
std::vector<SuperBlockId> member_counts(const std::vector<SuperBlockId>& mapping);
void init_thresholds();
double round_value(const std::vector<std::string>& par, const std::string& name, int round, int round_count);

template<typename Int, typename Int2>
std::vector<Int2> convert_mapping(const std::vector<Int>& mapping, Int2) {
	std::vector<Int2> out;
	out.reserve(mapping.size());
	std::transform(mapping.begin(), mapping.end(), std::back_inserter(out), [](Int x) { return (Int2)x; });
	return out;
}

template<typename It, typename It2>
int64_t update_clustering(It clustering, It2 mapping, It2 query_begin, It2 query_end, It2 db_begin) {
	const uint64_t n = query_end - query_begin;
	int64_t k = 0;
	for (OId i = 0; i < n; ++i)
		if (mapping[i] >= 0 && clustering[query_begin[i]] != db_begin[mapping[i]]) {
			clustering[query_begin[i]] = (typename It::value_type)db_begin[mapping[i]];
			++k;
		}
	return k;
}

template<typename It>
std::vector<OId> cluster_members(It begin, It end, const FlatArray<OId>& clusters) {
	std::vector<OId> out;
	for (It i = begin; i != end; ++i) {
		for (auto it = clusters.cbegin(*i); it != clusters.cend(*i); ++it) {
			out.push_back(*it);
		}
	}
	return out;
}

struct Cfg {
	Cfg(HspValues hsp_values, bool lazy_titles, const FlatArray<OId>& clusters, const std::vector<OId>& centroids, SequenceFile* db) :
		hsp_values(hsp_values),
		lazy_titles(lazy_titles),
		clusters(clusters),
		centroids(centroids),
		db(db)
	{
	}
	const HspValues hsp_values;
	const bool lazy_titles;
	const FlatArray<OId>& clusters;
	const std::vector<OId>& centroids;
	SequenceFile* db;
	std::shared_ptr<Block> centroid_block, member_block;
};

File* realign_block_pair(CentroidId begin, CentroidId end, Cfg& cfg);
std::vector<std::string> cluster_steps(double approx_id, bool linear);
std::vector<std::string> default_round_cov(int steps);
bool is_linclust(const std::vector<std::string>& steps);
std::vector<std::string> default_round_approx_id(int steps);
int round_ccd(std::vector<std::string> param, int round, int round_count, bool linear);

}