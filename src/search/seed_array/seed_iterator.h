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
#include <deque>
#include <type_traits>
#include <vector>
#include "basic/shape.h"
#include "basic/shape_config.h"
#include "basic/sequence.h"
#include "util/hash_function.h"

/* Tag requesting that seeds containing a masked letter be rejected instead of being
   encoded as if the letter were an ordinary residue. */
struct FilterMaskedSeeds { };

template<typename It>
struct SeedIterator
{
	SeedIterator(It begin, It end, const Shape& sh):
		ptr_ (begin),
		end_ (end - sh.length_ + 1)
	{}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed, const Shape &sh)
	{
		return sh.set_seed_reduced(seed, ptr_++);
	}
	SeedIterator& operator++() {
		ptr_++;
		return *this;
	}
private:
	It ptr_;
	const It end_;
};

template<typename It>
struct MinimizerIterator
{
	MinimizerIterator(It begin, It end, const Shape& sh, Loc window) :
		ptr_(begin),
		begin_(begin),
		end_(end - sh.length_ + 1),
		window_(window),
		sh_(sh)
	{
		next();
		if (good())
			min_idx_ = get();
	}
	bool good() const
	{
		return (Loc)seeds_.size() == window_;
	}
	uint64_t operator*() const {
		return seeds_[min_idx_];
	}
	MinimizerIterator& operator++() {
		int m = 0;
		const uint64_t current = **this;
		do {
			seeds_.pop_front();
			hashes_.pop_front();
			pos_.pop_front();
			next();
		} while (good() && seeds_[m = get()] == current);
		min_idx_ = m;
		return *this;
	}
	Loc pos() const {
		return pos_[min_idx_];
	}
private:
	void next() {
		while ((Loc)seeds_.size() < window_ && ptr_ < end_) {
			uint64_t s;
			if (sh_.set_seed_reduced(s, ptr_)) {
				seeds_.push_back(s);
				hashes_.push_back(MurmurHash()(s));
				pos_.push_back(Loc(ptr_ - begin_));
			}
			++ptr_;
		}
	}
	int get() const {
		std::deque<uint64_t>::const_iterator i = hashes_.begin(), j = i;
		uint64_t s = *i++;
		while (i < hashes_.end()) {
			if (*i < s) {
				s = *i;
				j = i;
			}
			++i;
		}
		return int(j - hashes_.begin());
	}
	It ptr_, begin_, end_;
	std::deque<uint64_t> seeds_, hashes_;
	std::deque<Loc> pos_;
	const Loc window_;
	const Shape& sh_;
	int min_idx_;
};

struct SketchIterator
{
	template<typename It>
	SketchIterator(It begin, It end, const Shape& sh, Loc n)
	{
		std::vector<Kmer> v;
		v.reserve((end - begin) - sh.length_ + 1);
		const It end2 = end - sh.length_ + 1;
		uint64_t s;
		for (It p = begin; p < end2; ++p)
			if (sh.set_seed_reduced(s, p))
				v.emplace_back(s, MurmurHash()(s), Loc(p - begin));
		std::sort(v.begin(), v.end());
		data_.insert(data_.end(), v.begin(), v.begin() + std::min(n, (Loc)v.size()));
		it_ = data_.begin();
	}
	bool good() const {
		return it_ < data_.end();
	}
	uint64_t operator*() const {
		return it_->seed;
	}
	Loc pos() const {
		return it_->pos;
	}
	SketchIterator& operator++() {
		++it_;
		return *this;
	}
private:
	struct Kmer {
		Kmer(uint64_t seed, uint64_t hash, Loc pos) : seed(seed), hash(hash), pos(pos) {}
		uint64_t seed, hash;
		Loc pos;
		bool operator<(const Kmer& k) const {
			if (hash != k.hash)
				return hash < k.hash;
			return pos < k.pos;
		}
	};
	std::vector<Kmer> data_;
	std::vector<Kmer>::const_iterator it_;
};

/* Computes the spaced seed of a window by masking instead of by gathering the letters of
   the shape one by one: the reduced letters of the window are kept in the b-bit fields of
   a machine word, which the shape masks down to its match positions in a single AND, and
   the result is hashed. A window is advanced by one shift and one OR, so the cost per seed
   no longer depends on the weight of the shape.

   The seed is the hash value, not the packed letters, so distinct seeds can collide. It is
   therefore only usable where a collision costs an additional candidate that a later stage
   discards, not correctness.

   Requires the b-bit fields of all length_ letters of the shape to fit into a 64 bit word,
   see hashed_seeds_supported(). */
template<uint64_t B, typename Filter = void>
struct HashedSeedIterator
{
	/* Rejecting windows that have a masked letter on a match position of the shape, which
	   is what the letterwise encoding does, costs a second shift register holding one bit
	   per letter of the window. It is compiled out unless it is asked for. */
	static constexpr bool FILTER_MASKED = std::is_same<Filter, FilterMaskedSeeds>::value;

	HashedSeedIterator(Letter* seq, Loc len, const Shape &sh):
		long_mask(sh.long_mask()),
		shape_mask_(sh.rev_mask_),
		ptr_(seq),
		end_(ptr_ + len),
		last_(0),
		masked_(0)
	{
		for (int i = 0; i < sh.length_ && ptr_ < end_; ++i) {
			const Letter l = letter_mask(*(ptr_++));
			push(l, !is_amino_acid(l));
		}
	}
	bool good() const
	{
		return ptr_ <= end_;
	}
	/* True if the seed of the current window is defined. Always true unless
	   FilterMaskedSeeds was requested. */
	bool valid() const {
		return !FILTER_MASKED || (masked_ & shape_mask_) == 0;
	}
	uint64_t operator*() const {
		return MurmurHash()(last_ & long_mask);
	}
	HashedSeedIterator& operator++() {
		while (ptr_ < end_) {
			const Letter l = letter_mask(*(ptr_++));
			const bool masked = !is_amino_acid(l);
			push(l, masked);
			// The last position of a shape is always a match position, so a window ending
			// on a masked letter never has a seed and is skipped outright.
			if (masked)
				continue;
			return *this;
		}
		++ptr_;
		return *this;
	}
	Letter* seq_ptr(const Shape& sh) const {
		return ptr_ - sh.length_;
	}
private:
	/* Shifts one letter into the window. A masked letter contributes an empty field: its
	   reduced value does not fit into b bits and would corrupt the field of its
	   predecessor. */
	void push(const Letter l, const bool masked) {
		last_ = (last_ << B) | (masked ? uint64_t(0) : (uint64_t)Reduction::get_reduction()(l));
		if (FILTER_MASKED)
			masked_ = (masked_ << 1) | (uint32_t)masked;
	}
	const uint64_t long_mask;
	const uint32_t shape_mask_;
	Letter *ptr_, *end_;
	uint64_t last_;
	uint32_t masked_;
};

/* Recomputes the seed that HashedSeedIterator yields for the window starting at seq.
   Returns false if the window has no seed. */
template<uint64_t B>
static inline bool hashed_seed(const Letter* seq, const Shape& sh, uint64_t& seed)
{
	uint64_t last = 0;
	uint32_t masked = 0;
	for (int i = 0; i < sh.length_; ++i) {
		const Letter l = letter_mask(seq[i]);
		const bool m = !is_amino_acid(l);
		last = (last << B) | (m ? uint64_t(0) : (uint64_t)Reduction::get_reduction()(l));
		masked = (masked << 1) | (uint32_t)m;
	}
	if (masked & sh.rev_mask_)
		return false;
	seed = MurmurHash()(last & sh.long_mask());
	return true;
}

/* True if the shapes and the reduction in use allow the hashed seed encoding, i.e. if the
   b-bit fields of a whole shape fit into the 64 bit window of HashedSeedIterator. The
   enumeration is only instantiated for a four bit reduction (see enum_seeds_worker). */
static inline bool hashed_seeds_supported()
{
	const int b = Reduction::get_reduction().bit_size();
	if (b != 4)
		return false;
	for (int i = 0; i < shapes.count(); ++i)
		if ((int64_t)shapes[i].length_ * b > 64)
			return false;
	return true;
}

/* Sketch of the n smallest seeds of a sequence over the hashed seed encoding. Mirrors
   SketchIterator, except that the seed is already a hash value and is therefore its own
   sort key. */
template<uint64_t B, typename Filter = void>
struct HashedSketchIterator
{
	HashedSketchIterator(Letter* seq, Loc len, const Shape& sh, Loc n)
	{
		std::vector<Kmer> v;
		v.reserve(std::max(len - sh.length_ + 1, 0));
		for (HashedSeedIterator<B, Filter> it(seq, len, sh); it.good(); ++it)
			if (it.valid())
				v.emplace_back(*it, Loc(it.seq_ptr(sh) - seq));
		std::sort(v.begin(), v.end());
		data_.insert(data_.end(), v.begin(), v.begin() + std::min(n, (Loc)v.size()));
		it_ = data_.begin();
	}
	bool good() const {
		return it_ < data_.end();
	}
	bool valid() const {
		return true;
	}
	uint64_t operator*() const {
		return it_->seed;
	}
	Loc pos() const {
		return it_->pos;
	}
	HashedSketchIterator& operator++() {
		++it_;
		return *this;
	}
private:
	struct Kmer {
		Kmer(uint64_t seed, Loc pos) : seed(seed), pos(pos) {}
		uint64_t seed;
		Loc pos;
		bool operator<(const Kmer& k) const {
			if (seed != k.seed)
				return seed < k.seed;
			return pos < k.pos;
		}
	};
	std::vector<Kmer> data_;
	std::vector<Kmer>::const_iterator it_;
};

template<int L, uint64_t B, typename Filter>
struct ContiguousSeedIterator
{
	ContiguousSeedIterator(const Sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0)
	{
		for (int i = 0; i < L - 1; ++i)
			last_ = (last_ << B) | Reduction::get_reduction()(letter_mask(*(ptr_++)));
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		for (;;) {
			last_ <<= B;
			last_ &= (uint64_t(1) << (B*L)) - 1;
			const Letter l = letter_mask(*(ptr_++));
			last_ |= Reduction::get_reduction()(l);
			seed = last_;
			return true;
		}
	}
	static int length()
	{
		return L;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
	unsigned mask_;
};

template<int L, uint64_t B>
struct ContiguousSeedIterator<L, B, FilterMaskedSeeds>
{
	ContiguousSeedIterator(const Sequence &seq) :
		ptr_(seq.data()),
		end_(ptr_ + seq.length()),
		last_(0),
		mask_(0)
	{
		for (int i = 0; i < L - 1; ++i) {
			const Letter l = letter_mask(*(ptr_++));
			last_ = (last_ << B) | Reduction::get_reduction()(l);
			if (!is_amino_acid(l))
				mask_ |= 1;
			mask_ <<= 1;
		}
	}
	bool good() const
	{
		return ptr_ < end_;
	}
	bool get(uint64_t &seed)
	{
		for (;;) {
			last_ <<= B;
			last_ &= (uint64_t(1) << (B*L)) - 1;
			mask_ <<= 1;
			mask_ &= (1 << L) - 1;
			const Letter l = letter_mask(*(ptr_++));
			const unsigned r = Reduction::get_reduction()(l);
			last_ |= r;
			seed = last_;
			if (!is_amino_acid(l))
				mask_ |= 1;
			return mask_ == 0;
		}
	}
	static int length()
	{
		return L;
	}
private:
	const Letter *ptr_, *end_;
	uint64_t last_;
	unsigned mask_;
};