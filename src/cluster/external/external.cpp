#include <type_traits>
#include "../basic/config.h"
#include "util/tsv/tsv.h"
#include "util/parallel/atomic.h"
#include "util/log_stream.h"
#include "data/sequence_file.h"
#include "basic/reduction.h"
#include "basic/shape_config.h"
#include "basic/seed_iterator.h"
#include "search/search.h"
#include "external.h"
#include "radix_sort.h"
#define _REENTRANT
#include "lib/ips4o/ips4o.hpp"
#include "util/algo/hyperloglog.h"
#include "util/sequence//sequence.h"
#include "util/string/string.h"
#include "../cascaded/cascaded.h"
#include "build_pair_table.h"

using std::unique_ptr;
using std::thread;
using std::string;
using std::vector;
using std::endl;
using std::array;
using std::mutex;
using std::lock_guard;
using std::pair;
using std::atomic;
using Util::String::format;
using std::shared_ptr;

namespace Cluster {

struct ChunkTableEntry {
	ChunkTableEntry():
		oid(),
		chunk()
	{}
	ChunkTableEntry(int64_t oid, int32_t chunk):
		oid(oid),
		chunk(chunk)
	{}
	int64_t key() const {
		return oid;
	}
	bool operator<(const ChunkTableEntry& e) const {
		return oid < e.oid || (oid == e.oid && chunk < e.chunk);
	}
	int64_t oid;
	int32_t chunk;
};

static vector<string> build_seed_table(Job& job, const VolumedFile& volumes, int shape) {
	const int64_t BUF_SIZE = 4096;
	const double SKETCH_SIZE_RATIO = 0.1;
	//Reduction::set_reduction(Search::no_reduction);
	Loc sketch_size = config.sketch_size == 0 ? Search::sensitivity_traits.at(config.sensitivity).sketch_size : config.sketch_size;
	if (sketch_size == 0)
		sketch_size = std::numeric_limits<Loc>::max();
	//const uint64_t shift = shapes[shape].bit_length() - RADIX_BITS;

	const std::string base_dir = job.base_dir() + PATH_SEPARATOR + "seed_table_" + std::to_string(shape) + PATH_SEPARATOR, qpath = base_dir + "queue";
	mkdir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, RADIX_COUNT, job.worker_id()));
	
	Atomic q(qpath);
	atomic<int> volumes_processed(0);
	vector<thread> workers;
	auto worker = [&](int thread_id) {
		BufferArray buffers(*output_files, RADIX_COUNT);
		int64_t v = 0;
		vector<Letter> buf;
		while (v = q.fetch_add(), v < (int64_t)volumes.size()) {
			job.log("Building seed table. Shape=%i/%i Volume=%lli/%lli Records=%s", shape + 1, ::shapes.count(), v + 1, volumes.size(), format(volumes[v].record_count).c_str());
			unique_ptr<OutputFile> oid_out;
			if (job.round() > 0)
				oid_out.reset(new OutputFile(volumes[v].path + ".oid"));
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[v].path }));
			string id;
			vector<Letter> seq;
			int64_t oid = volumes[v].oid_begin;
			while (in->read_seq(seq, id, nullptr)) {
				if (job.round() > 0) {
					const int64_t prev_oid = atoll(id.c_str());
					oid_out->write(&prev_oid, 1);
				}
				Reduction::reduce_seq(Sequence(seq), buf);
				const Shape& sh = shapes[shape];
				if (seq.size() < (size_t)sh.length_) {
					++oid;
					continue;
				}
				//SketchIterator it(buf, sh, std::min(sketch_size, std::max((int)std::round(seq.size() * SKETCH_SIZE_RATIO), 1)));
				SketchIterator it(buf, sh, sketch_size);
				while (it.good()) {
					const uint64_t key = *it, radix = MurmurHash()(key) & (RADIX_COUNT - 1);
					buffers.write(radix, SeedEntry(key, oid, (int32_t)seq.size()));
					++it;
				}
				++oid;
			}
			in->close();
			if (job.round() > 0)
				oid_out->close();
			volumes_processed.fetch_add(1, std::memory_order_relaxed);
		}
		};
	for (int i = 0; i < config.threads_; ++i)
		workers.emplace_back(worker, i);
	for (auto& t : workers)
		t.join();
	const vector<string> buckets = output_files->buckets();
	TaskTimer timer("Closing the output files");
	output_files.reset();
	Atomic finished(base_dir + "finished");
	finished.fetch_add(volumes_processed);
	finished.await(volumes.size());
	return buckets;
}

static vector<string> build_pair_table(Job& job, const vector<string>& seed_table, int64_t db_size, FileArray& output_files) {
	const int64_t BUF_SIZE = 4096, shift = bit_length(db_size - 1) - RADIX_BITS; // promiscuous_cutoff = db_size / config.promiscuous_seed_ratio;
	const string queue_path = base_path(seed_table.front()) + PATH_SEPARATOR + "build_pair_table_queue";
	const bool unid = !config.mutual_cover.present();
	Atomic queue(queue_path);
	int64_t bucket, buckets_processed = 0;
	while (bucket = queue.fetch_add(), bucket < (int64_t)seed_table.size()) {
		VolumedFile file(seed_table[bucket]);
		InputBuffer<SeedEntry> data(file);
		job.log("Building pair table. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, seed_table.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		ips4o::parallel::sort(data.begin(), data.end(), std::less<SeedEntry>(), config.threads_);
		auto worker = [&](int thread_id) {
			BufferArray buffers(output_files, RADIX_COUNT);
			auto it = merge_keys(data.begin(thread_id), data.end(thread_id), SeedEntry::Key());
			while (it.good()) {
				/*if (it.count() >= promiscuous_cutoff) {
					++it;
					continue;
				}*/
				if (unid)
					get_pairs_uni_cov(it, buffers);
				else
					get_pairs_mutual_cov(it, buffers);
				++it;
			}
			};
		vector<thread> workers;
		for (int i = 0; i < data.parts(); ++i)
			workers.emplace_back(worker, i);
		for (auto& t : workers)
			t.join();
		file.remove();
		++buckets_processed;
	}
	const vector<string> buckets = output_files.buckets();
	Atomic finished(base_path(seed_table.front()) + PATH_SEPARATOR + "pair_table_finished");
	finished.fetch_add(buckets_processed);
	finished.await(seed_table.size());
	return buckets;
}

struct SizeCounter {
	void add(int64_t oid, int32_t len) {
		const int64_t x = oid << 17, n = x + int64_t(len + 63) / 64;
		for (int64_t i = x; i < n; ++i)
			hll.add(i);
	}
	HyperLogLog hll;
};

struct Chunk {
	Chunk(Atomic& next_chunk, const string& chunks_path):
		id(next_chunk.fetch_add())
	{
		mkdir(chunks_path + std::to_string(id));
		pairs_out.reset(new OutputFile(chunks_path + std::to_string(id) + PATH_SEPARATOR + "pairs"));
	}
	void write(vector<PairEntryShort>& pairs_buffer, SizeCounter& size) {
		lock_guard<mutex> lock(mtx);
		pairs_out->write(pairs_buffer.size());
		pairs_out->write(pairs_buffer.data(), pairs_buffer.size());
		pairs_buffer.clear();
		this->size.merge(size.hll);
		size.hll = HyperLogLog();
	}
	~Chunk() {
		pairs_out->close();
	}
	const int id;
	unique_ptr<OutputFile> pairs_out;
	HyperLogLog size;
	mutex mtx;
};

static pair<vector<string>, int> build_chunk_table(Job& job, const vector<string>& pair_table, int64_t db_size) {
	const int64_t BUF_SIZE = 4096, shift = bit_length(db_size - 1) - RADIX_BITS, max_chunk_size = Util::String::interpret_number(config.linclust_chunk_size) / 64,
		max_processed = std::max(std::min(INT64_C(262144), max_chunk_size / config.threads_ / 16), INT64_C(1)); // ???
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "chunk_table",
		chunks_path = job.base_dir() + PATH_SEPARATOR + "chunks" + PATH_SEPARATOR;
	mkdir(base_path);
	mkdir(chunks_path);
	unique_ptr<FileArray> output_files(new FileArray(base_path, RADIX_COUNT, job.worker_id()));
	Atomic queue(base_path + PATH_SEPARATOR + "queue");
	Atomic next_chunk(base_path + PATH_SEPARATOR + "next_chunk");
	shared_ptr<Chunk> current_chunk(new Chunk(next_chunk, chunks_path));
	int64_t bucket, total_pairs = 0, total_distinct_pairs = 0, buckets_processed = 0;
	mutex mtx;
	while (bucket = queue.fetch_add(), bucket < (int64_t)pair_table.size()) {
		VolumedFile file(pair_table[bucket]);
		InputBuffer<PairEntry> data(file);
		job.log("Building chunk table. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, pair_table.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		total_pairs += data.size();
		ips4o::parallel::sort(data.begin(), data.end(), std::less<PairEntry>(), config.threads_);
		auto worker = [&](int thread_id) {
			shared_ptr<Chunk> my_chunk(current_chunk);
			BufferArray buffers(*output_files, RADIX_COUNT);
			vector<PairEntryShort> pairs_buffer;
			SizeCounter size;
			int64_t distinct_pairs = 0, processed = 0;
			auto it = merge_keys(data.begin(thread_id), data.end(thread_id), PairEntry::Key());
			while (it.good()) {
				const int64_t rep_oid = it.begin()->rep_oid, radix = rep_oid >> shift; // Hashing?
				buffers.write(radix, ChunkTableEntry(rep_oid, my_chunk->id));
				size.add(rep_oid, it.begin()->rep_len);
				processed += it.begin()->rep_len;
				for (auto j = it.begin(); j < it.end(); ++j) {
					if (j > it.begin() && j->member_oid == (j - 1)->member_oid)
						continue;
					const int64_t radix = j->member_oid >> shift;
					buffers.write(radix, ChunkTableEntry(j->member_oid, my_chunk->id));
					size.add(j->member_oid, j->member_len);
					pairs_buffer.emplace_back(rep_oid, j->member_oid);
					++distinct_pairs;
					processed += j->member_len;
					if (processed >= max_processed) {
						my_chunk->write(pairs_buffer, size);
						processed = 0;
						bool new_chunk = false;
						if (my_chunk != current_chunk) {
							my_chunk = current_chunk;
							new_chunk = true;
						}
						else {
							const int64_t est = my_chunk->size.estimate();
							if (est >= max_chunk_size) {
								lock_guard<mutex> lock(mtx);
								if (my_chunk == current_chunk) {
									log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << est * 64 << endl;
									current_chunk.reset(new Chunk(next_chunk, chunks_path));
									my_chunk = current_chunk;
									new_chunk = true;
								}
								else {
									// should not happen?
								}
							}
						}
						if (new_chunk) {
							buffers.write(rep_oid >> shift, ChunkTableEntry(rep_oid, my_chunk->id));
							size.add(rep_oid, it.begin()->rep_len);
							processed += it.begin()->rep_len;
						}
					}
				}
				++it;
			}
			my_chunk->write(pairs_buffer, size);
			total_distinct_pairs += distinct_pairs;
			};
		vector<thread> workers;
		for (int i = 0; i < data.parts(); ++i)
			workers.emplace_back(worker, i);
		for (auto& t : workers)
			t.join();
		const int64_t est = current_chunk->size.estimate();
		if (est >= max_chunk_size) {
			log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << est * 64 << endl;
			current_chunk.reset(new Chunk(next_chunk, chunks_path));
		}
		file.remove();
		++buckets_processed;
	}
	log_stream << "build_chunk_table chunk=" << current_chunk->id << " est_size=" << current_chunk->size.estimate() << " total_pairs=" << total_pairs << " total_distinct_pairs=" << total_distinct_pairs << endl;
	const vector<string> buckets = output_files->buckets();
	TaskTimer timer("Closing the output files");
	output_files.reset();
	current_chunk.reset();
	timer.go("Waiting for other workers");
	Atomic finished(base_path + PATH_SEPARATOR + "finished");
	finished.fetch_add(buckets_processed);
	finished.await(pair_table.size());
	return { buckets, next_chunk.get() };
}

static void build_chunks(Job& job, const VolumedFile& db, const vector<string>& chunk_table, int chunk_count) {
	const int64_t BUF_SIZE = 64 * 1024;
	const std::string base_path = job.base_dir() + PATH_SEPARATOR + "chunks" + PATH_SEPARATOR,
		queue_path = base_path + "queue";
	unique_ptr<FileArray> output_files(new FileArray(base_path, chunk_count, job.worker_id(), 1024 * 1024 * 1024));
	Atomic queue(queue_path);
	int64_t bucket, buckets_processed = 0;
	atomic<int64_t> oid_counter(0), distinct_oid_counter(0);
	while (bucket = queue.fetch_add(), bucket < (int64_t)chunk_table.size()) {
		VolumedFile file(chunk_table[bucket]);
		InputBuffer<ChunkTableEntry> data(file);
		job.log("Building chunks. Bucket=%lli/%lli Records=%s Size=%s", bucket + 1, chunk_table.size(), format(data.size()).c_str(), format(data.byte_size()).c_str());
		ips4o::parallel::sort(data.begin(), data.end(), std::less<ChunkTableEntry>(), config.threads_);
		const int64_t oid_begin = data.front().oid, oid_end = data.back().oid + 1;
		const pair<vector<Volume>::const_iterator, vector<Volume>::const_iterator> volumes = db.find(oid_begin, oid_end);
		atomic<int64_t> next(0);
		auto worker = [&]() {
			int64_t volume;
			const ChunkTableEntry* table_ptr = data.begin();
			BufferArray output_bufs(*output_files, chunk_count);
			TextBuffer buf;
			while (volume = next.fetch_add(1, std::memory_order_relaxed), volume < volumes.second - volumes.first) {
				const Volume& v = volumes.first[volume];
				while (table_ptr->oid < v.oid_begin)
					++table_ptr;
				unique_ptr<SequenceFile> in(SequenceFile::auto_create({ v.path }));
				string id;
				vector<Letter> seq;
				int64_t file_oid = v.oid_begin;
				while (file_oid < oid_end && in->read_seq(seq, id, nullptr)) {
					if (table_ptr->oid > file_oid) {
						++file_oid;
						continue;
					}
					Util::Seq::format(seq, std::to_string(file_oid).c_str(), nullptr, buf, "fasta", amino_acid_traits);
					const ChunkTableEntry* begin = table_ptr;
					while (table_ptr < data.end() && table_ptr->oid == file_oid) {
						if (table_ptr == begin || table_ptr->chunk != table_ptr[-1].chunk) {
							output_bufs.write(table_ptr->chunk, buf.data(), buf.size());
							oid_counter.fetch_add(1, std::memory_order_relaxed);
						}
						++table_ptr;
					}
					buf.clear();
					distinct_oid_counter.fetch_add(1, std::memory_order_relaxed);
					++file_oid;
				}
				in->close();
			}
		};
		vector<thread> workers;
		for (int i = 0; i < std::min(config.threads_, int(volumes.second - volumes.first)); ++i)
			workers.emplace_back(worker);
		for (auto& t : workers)
			t.join();
		file.remove();
		++buckets_processed;
	}
	TaskTimer timer("Closing the output files");
	output_files.reset();
	timer.go("Waiting for other workers");
	Atomic finished(base_path + "finished");
	finished.fetch_add(buckets_processed);
	finished.await(chunk_table.size());
	timer.finish();
	log_stream << "build_chunks oids=" << oid_counter << '/' << db.records() << " distinct_oids=" << distinct_oid_counter << endl;
	db.remove();
}

string round(Job& job, const VolumedFile& volumes) {
	::shapes = ShapeConfig(Search::shape_codes.at(config.sensitivity), 0);
	job.log("Starting round %i sensitivity %s %i shapes\n", job.round(), to_string(config.sensitivity).c_str(), ::shapes.count());
	job.set_round(volumes.size(), volumes.records());
	const std::string pair_table_base = job.base_dir() + PATH_SEPARATOR + "pair_table";
	mkdir(pair_table_base);
	unique_ptr<FileArray> pair_table_files(new FileArray(pair_table_base, RADIX_COUNT, job.worker_id()));
	vector<string> pair_table;
	for (int shape = 0; shape < ::shapes.count(); ++shape) {
		const vector<string> buckets = build_seed_table(job, volumes, shape);
		const vector<string> sorted_seed_table = radix_sort<SeedEntry>(job, buckets, shapes[0].bit_length() - RADIX_BITS);
		pair_table = build_pair_table(job, sorted_seed_table, volumes.records(), *pair_table_files);
	}
	pair_table_files.reset();
	const vector<string> sorted_pair_table = radix_sort<PairEntry>(job, pair_table, bit_length(volumes.records() - 1) - RADIX_BITS);
	const pair<vector<string>, int> chunk_table = build_chunk_table(job, sorted_pair_table, volumes.records());
	const vector<string> sorted_chunk_table = radix_sort<ChunkTableEntry>(job, chunk_table.first, bit_length(volumes.records() - 1) - RADIX_BITS);
	build_chunks(job, volumes, sorted_chunk_table, chunk_table.second);
	const vector<string> edges = align(job, chunk_table.second, volumes.records());
	if (config.mutual_cover.present()) {
		return cluster_bidirectional(job, edges, volumes);
	}
	else {
		const vector<string> sorted_edges = radix_sort<Edge>(job, edges, bit_length(volumes.records() - 1) - RADIX_BITS);
		//const vector<string> sorted_edges = read_list(job.base_dir() + PATH_SEPARATOR + "alignments" + PATH_SEPARATOR + "radix_sort_out");
		return cluster(job, sorted_edges, volumes);
	}
}

void external() {
	if (config.output_file.empty())
		throw std::runtime_error("Option missing: output file (--out/-o)");
	TaskTimer total;
	Job job;
	VolumedFile volumes(config.database.get_present());
	if (job.worker_id() == 0) {
		if(config.mutual_cover.present())
			job.log("Bi-directional coverage = %f", config.mutual_cover.get_present());
		else
			job.log("Uni-directional coverage = %f", config.member_cover.get(80));
		job.log("Approx. id = %f", config.approx_min_id.get(0));
		job.log("#Volumes = %lli", volumes.size());
		job.log("#Sequences = %lli", volumes.records());
	}
	if (config.mutual_cover.present()) {
		config.min_length_ratio = std::min(config.mutual_cover.get_present() / 100 + 0.05, 1.0);
		config.query_or_target_cover = 0;
		config.query_cover = config.mutual_cover.get_present();
		config.subject_cover = config.mutual_cover.get_present();
	}
	else {
		config.query_or_target_cover = config.member_cover.get(80);
		config.query_cover = 0;
		config.subject_cover = 0;
	}
#ifdef WIN32
	_setmaxstdio(8192);
#endif
	vector<string> steps = Cluster::cluster_steps(config.approx_min_id.get(0), true);
	string reps;
	job.set_round_count((int)steps.size());
	for (size_t i = 0; i < steps.size(); ++i) {
		config.sensitivity = from_string<Sensitivity>(rstrip(steps[i], "_lin"));
		reps = round(job, i == 0 ? volumes : VolumedFile(reps));
		if (i < steps.size() - 1)
			job.next_round();
	}
	Atomic output_lock(job.base_dir() + PATH_SEPARATOR + "output_lock");
	if(output_lock.fetch_add() == 0)
		output(job);
	log_stream << "Total time = " << (double)total.milliseconds() / 1000 << 's' << endl;
}

}