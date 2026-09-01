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

#include <sstream>
#include <utility>
#include "volumed_fasta_file.h"
#include "fasta_file.h"
#include "cluster/volume.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "util/string/string.h"

using std::endl;
using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;

static const char* VOLUME_LIST_EXTENSION = ".tsv";

static const SequenceFile::Flags UNSUPPORTED_FLAGS = SequenceFile::Flags::TAXON_MAPPING
	| SequenceFile::Flags::TAXON_NODES
	| SequenceFile::Flags::TAXON_RANKS
	| SequenceFile::Flags::TAXON_SCIENTIFIC_NAMES
	| SequenceFile::Flags::NEED_EARLY_TAXON_MAPPING
	| SequenceFile::Flags::ACC_TO_OID_MAPPING
	| SequenceFile::Flags::OID_TO_ACC_MAPPING
	| SequenceFile::Flags::NEED_LENGTH_LOOKUP;

namespace {

struct EmptyChunk : public RawChunk {
	virtual bool empty() const override {
		return true;
	}
	virtual OId begin() const noexcept override {
		return 0;
	}
	virtual OId end() const noexcept override {
		return 0;
	}
	virtual DecodedPackage* decode(SequenceFile::Flags flags, const BitVector* filter, std::unordered_map<string, bool>* accs, SequenceType seq_type) const override {
		throw OperationNotSupported();
	}
	virtual size_t letters() const noexcept override {
		return 0;
	}
	virtual size_t bytes() const noexcept override {
		return 0;
	}
};

}

bool VolumedFastaFile::is_volume_list(const string& path) {
	return ends_with(path, VOLUME_LIST_EXTENSION);
}

VolumedFastaFile::VolumedFastaFile(const string& file_name, Flags flags, const ValueTraits& value_traits) :
	SequenceFile(SequenceFile::Type::VOLUMED, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	file_name_(file_name),
	volumes_(new VolumedFile(file_name)),
	current_(0),
	oid_(0),
	raw_chunk_no_(0),
	have_seq_count_(false),
	have_letter_count_(false),
	seqs_(0),
	letters_(0)
{
	if (flag_any(flags, UNSUPPORTED_FLAGS))
		throw runtime_error("The FASTA volume list format does not support taxonomic features or filtering by accession.");
	if (volumes_->empty())
		throw runtime_error("Empty volume list file: " + file_name);

	const bool have_counts = volumes_->have_record_counts();
	disk_size_ = 0;
	files_.reserve(volumes_->size());	
	vector<const Volume*> volumes;
	volumes.reserve(volumes_->size());
	int64_t empty_volumes = 0;
	for (const Volume& v : *volumes_) {
		if (have_counts && v.record_count == std::numeric_limits<OId>::max())
			throw runtime_error("Invalid sequence count for volume " + v.path + " in " + file_name);
		unique_ptr<FastaFile> f;
		try {
			f.reset(new FastaFile({ v.path }, flags, value_traits));
		}
		catch (const EmptyFileError&) {
			if (have_counts && v.record_count > 0)
				throw runtime_error("Volume " + v.path + " is empty but a sequence count of "
					+ std::to_string(v.record_count) + " is given for it in " + file_name);
			++empty_volumes;
			continue;
		}
		if (!f->is_fasta())
			throw runtime_error("Only the FASTA format is supported for volumes of a volume list: " + v.path);
		disk_size_ += f->disk_size();
		files_.push_back(std::move(f));
		volumes.push_back(&v);
	}

	seqs_ = have_counts ? (uint64_t)volumes_->sparse_records() : 0;
	have_seq_count_ = have_counts;
	if (flag_any(flags, Flags::NEED_LETTER_COUNT)) {
		have_seq_count_ = have_letter_count_ = true;
		uint64_t seqs = 0;
		for (int64_t i = 0; i < (int64_t)files_.size(); ++i) {
			const uint64_t n = files_[i]->sequence_count().value_or(0);
			if (have_counts && n != (uint64_t)volumes[i]->record_count)
				throw runtime_error("Sequence count of volume " + volumes[i]->path + " (" + std::to_string(n)
					+ ") does not match the count given in " + file_name + " (" + std::to_string(volumes[i]->record_count) + ")");
			seqs += n;
			letters_ += files_[i]->letters().value_or(0);
		}
		seqs_ = seqs;
	}

	std::ostringstream ss;
	ss << "FASTA volumes: " << files_.size() << ", total size: " << disk_size_ << " bytes" << endl;
	if (empty_volumes > 0)
		ss << "Empty FASTA volumes skipped: " << empty_volumes << endl;
	open_stats_ += ss.str();
	for (auto& f : files_)
		f->rewind();
}

VolumedFastaFile::~VolumedFastaFile() {
	close();
}

int64_t VolumedFastaFile::file_count() const {
	return (int64_t)files_.size();
}

void VolumedFastaFile::close() {
	for (auto& f : files_)
		f->close();
}

void VolumedFastaFile::set_seqinfo_ptr(OId i) {
	if (i == oid_)
		return;
	if (i != 0)
		throw runtime_error("The FASTA volume list format only supports sequential reading, seeking to OId " + std::to_string(i) + " was requested.");
	for (auto& f : files_)
		f->rewind();
	current_ = 0;
	oid_ = 0;
	raw_chunk_no_ = 0;
}

OId VolumedFastaFile::tell_seq() const {
	return oid_;
}

bool VolumedFastaFile::eof() const {
	return current_ >= (int64_t)files_.size()
		|| (current_ == (int64_t)files_.size() - 1 && files_[current_]->eof());
}

void VolumedFastaFile::init_seq_access() {
	set_seqinfo_ptr(0);
}

void VolumedFastaFile::print_info() const {
	*message_stream << "Database: " << file_name_ << ' ';
	*message_stream << "(type: " << to_string(type()) << ", volumes: " << files_.size();
	if (have_seq_count_)
		*message_stream << ", sequences: " << seqs_;
	if (have_letter_count_)
		*message_stream << ", letters: " << letters_;
	*message_stream << ')' << endl;
}

optional<uint64_t> VolumedFastaFile::sequence_count() const {
	return have_seq_count_ ? optional<uint64_t>(seqs_) : nullopt;
}

optional<uint64_t> VolumedFastaFile::letters() const {
	return have_letter_count_ ? optional<uint64_t>(letters_) : nullopt;
}

string VolumedFastaFile::file_name() {
	return file_name_;
}

bool VolumedFastaFile::read_seq(vector<Letter>& seq, string& id, vector<char>* quals) {
	while (current_ < (int64_t)files_.size()) {
		if (files_[current_]->read_seq(seq, id, quals)) {
			++oid_;
			return true;
		}
		++current_;
	}
	return false;
}

RawChunk* VolumedFastaFile::raw_chunk(size_t bytes, Flags flags) {
	while (current_ < (int64_t)files_.size()) {
		RawChunk* c = files_[current_]->raw_chunk(bytes, flags);
		if (!c->empty()) {
			c->no = raw_chunk_no_++;
			return c;
		}
		delete c;
		if (bytes == 0)
			break;
		++current_;
	}
	RawChunk* c = new EmptyChunk();
	c->no = raw_chunk_no_;
	return c;
}

int VolumedFastaFile::raw_chunk_no() const {
	return raw_chunk_no_;
}

void VolumedFastaFile::advance_seq_count(OId n) {
	oid_ += n;
}

double VolumedFastaFile::rel_file_ptr() {
	if (current_ >= (int64_t)files_.size())
		return 1.0;
	return ((double)current_ + files_[current_]->rel_file_ptr()) / (double)files_.size();
}

void VolumedFastaFile::end_random_access(bool dictionary) {
	if (!dictionary)
		return;
	free_dictionary();
}

void VolumedFastaFile::init_seqinfo_access() {
	throw OperationNotSupported();
}

void VolumedFastaFile::seek_chunk(const Chunk& chunk) {
	throw OperationNotSupported();
}

SequenceFile::SeqInfo VolumedFastaFile::read_seqinfo() {
	throw OperationNotSupported();
}

void VolumedFastaFile::putback_seqinfo() {
	throw OperationNotSupported();
}

size_t VolumedFastaFile::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) {
	throw OperationNotSupported();
}

void VolumedFastaFile::seek_offset(size_t p) {
	throw OperationNotSupported();
}

void VolumedFastaFile::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) {
	throw OperationNotSupported();
}

void VolumedFastaFile::read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) {
	throw OperationNotSupported();
}

void VolumedFastaFile::skip_id_data() {
	throw OperationNotSupported();
}

int VolumedFastaFile::db_version() const {
	throw OperationNotSupported();
}

int VolumedFastaFile::program_build_version() const {
	throw OperationNotSupported();
}

int VolumedFastaFile::build_version() {
	throw OperationNotSupported();
}

void VolumedFastaFile::create_partition_balanced(int64_t max_letters) {
	throw OperationNotSupported();
}

void VolumedFastaFile::save_partition(const string& partition_file_name, const string& annotation) {
	throw OperationNotSupported();
}

int VolumedFastaFile::get_n_partition_chunks() {
	throw OperationNotSupported();
}

DbFilter* VolumedFastaFile::filter_by_accession(const string& file_name) {
	throw runtime_error("The FASTA volume list format does not support filtering by accession.");
}

vector<TaxId> VolumedFastaFile::taxids(size_t oid) const {
	throw OperationNotSupported();
}

void VolumedFastaFile::seq_data(size_t oid, vector<Letter>& dst) {
	throw OperationNotSupported();
}

Loc VolumedFastaFile::seq_length(size_t oid) {
	throw OperationNotSupported();
}