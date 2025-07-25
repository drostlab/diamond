/****
DIAMOND protein aligner
Copyright (C) 2024 Max Planck Society for the Advancement of Science e.V.
				   Benjamin Buchfink

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#pragma once
#include <algorithm>
#include <utility>
#include "../io/text_input_file.h"

namespace Util { namespace String {

static inline const char* trim_cr(const char* begin, const char* end) {
	if (end > begin && end[-1] == '\r')
		return end - 1;
	else
		return end;
}

inline std::pair<const char*, int> find_newline(const char* begin, const char* end) {
	const char* p = std::find(begin, end, '\n');
	if (p == end)
		return { p,0 };
	if (p > begin && p[-1] == '\r')
		return { p + 1,2 };
	else
		return { p + 1,1 };
}

inline bool is_whitespace(const std::string& s) {
	for (char c : s)
		if (c != ' ' && c != '\n' && c != '\r')
			return false;
	return true;
}

struct TokenizerBase {
	virtual TokenizerBase* clone() const = 0;
	void reset(const char* begin, const char* end) {
		ptr_ = begin;
		end_ = end;
	}
	bool good() const {
		return ptr_ < end_;
	}
	const char* ptr() const {
		return ptr_;
	}
	virtual ~TokenizerBase() {}
	virtual std::string operator*() const = 0;
	virtual TokenizerBase& operator++() = 0;
	virtual std::pair<bool, std::string> read_record() = 0;
	virtual std::pair<bool, std::string> read_record(TextInputFile& file) = 0;
	const char* ptr_, * end_;
};

struct CharTokenizer : public TokenizerBase {
	CharTokenizer(char delimiter) :
		delimiter_(delimiter)
	{}
	virtual CharTokenizer* clone() const override
	{
		return new CharTokenizer(delimiter_);
	}
	virtual std::string operator*() const override {
		return std::string(ptr_, std::find(ptr_, end_, delimiter_));
	}
	virtual ~CharTokenizer() {}
	virtual CharTokenizer& operator++() override {
		ptr_ = std::find(ptr_, end_, delimiter_);
#ifndef NDEBUG
		if (ptr_ < end_)
#endif
			++ptr_;
		return *this;
	}
	virtual std::pair<bool, std::string> read_record() override {
		if (ptr_ >= end_)
			return { false,std::string() };
		const char* b = ptr_;
		ptr_ = std::find(ptr_, end_, '\n');
		if (ptr_ > b && ptr_[-1] == '\r') {			
			++ptr_;
			return { true, std::string(b, ptr_ - 2) };
		}
		else {
			++ptr_;
			return { true, std::string(b, ptr_ - 1) };
		}
	}
	virtual std::pair<bool, std::string> read_record(TextInputFile&) override {
		throw std::runtime_error("Invalid use of CharTokenizer");
	}
private:
	const char delimiter_;
};

/*struct MultiCharTokenizer : public TokenizerBase {
	MultiCharTokenizer(const char* delimiter) :
		delimiter_(delimiter)
	{}
	virtual MultiCharTokenizer* clone() const override
	{
		return new MultiCharTokenizer(delimiter_.c_str());
	}
	virtual std::string operator*() const override {
		std::cout << std::string(ptr_, std::search(ptr_, end_, delimiter_.begin(), delimiter_.end())) << std::endl;
		return std::string(ptr_, std::search(ptr_, end_, delimiter_.begin(), delimiter_.end()));
	}
	virtual ~MultiCharTokenizer() {}
	virtual MultiCharTokenizer& operator++() override {
		ptr_ = std::search(ptr_, end_, delimiter_.begin(), delimiter_.end());
#ifndef NDEBUG
		if (ptr_ < end_)
#endif
			ptr_ += delimiter_.length();
		return *this;
	}
private:
	const std::string delimiter_;
};

inline TokenizerBase* make_tokenizer(const std::string& delimiter) {
	assert(!delimiter.empty());
	return delimiter.length() == 1 ? static_cast<TokenizerBase*>(new CharTokenizer(delimiter[0])) : static_cast<TokenizerBase*>(new MultiCharTokenizer(delimiter.c_str()));
}*/

template<typename It, char delimiter>
struct TokenIterator {
	TokenIterator(It begin, It end) :
		ptr_(begin),
		end_(end)
	{}
	bool good() const {
		return ptr_ < end_;
	}
	std::string operator*() const {
		return std::string(ptr_, std::find(ptr_, end_, delimiter));
	}
	TokenIterator& operator++() {
		ptr_ = std::find(ptr_, end_, delimiter);
#ifndef NDEBUG
		if (ptr_ < end_)
#endif
			++ptr_;
		return *this;
	}
	It ptr() const {
		return ptr_;
	}
private:
	It ptr_, end_;
};

template<typename It> using TabIterator = TokenIterator<It, '\t'>;
using LineIterator = TokenIterator<const char*, '\n'>;

struct FastaTokenizer : public TokenizerBase
{
	virtual FastaTokenizer* clone() const override {
		return new FastaTokenizer();
	}
	virtual std::string operator*() const override {
		if (*ptr_ == '>') {
			const char* i = std::find(ptr_, end_, '\n');
			return std::string(ptr_ + 1, trim_cr(ptr_ + 1, i));
		}
		else {
			std::string r;
			r.reserve(end_ - ptr_);
			const char* p = ptr_;
			while (p < end_) {
				const char* i = std::find(p, end_, '\n');
				r.append(p, trim_cr(p, i));
				p = i + 1;
			}
			return r;
		}
	}
	virtual FastaTokenizer& operator++() override {
		if (*ptr_ == '>') {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
		}
		else {
			ptr_ = end_;
		}
		return *this;
	}
	virtual std::pair<bool, std::string> read_record() override {
		static const char* delimiter = "\n>";
		if (ptr_ >= end_)
			return { false,std::string() };
		const char* p = ptr_;
		ptr_ = std::search(ptr_, end_, delimiter, delimiter + 2);
		return { true, std::string(p, ptr_) };
	}
	virtual std::pair<bool, std::string> read_record(TextInputFile& file) override {
		if (file.eof())
			return { false,std::string() };
		std::string s;
		file.read_to(std::back_inserter(s), '\n', '>');
		return { !is_whitespace(s), s};
	}
};

struct MalformedFastqRecord : public std::runtime_error {
	MalformedFastqRecord(int64_t line = -1):
		std::runtime_error("Malformed FASTQ record" + (line == -1 ? "" : " on line " + std::to_string(line)))
	{}
};

struct FastqTokenizer : public TokenizerBase
{
	virtual FastqTokenizer* clone() const override {
		return new FastqTokenizer();
	}
	virtual std::string operator*() const override {
		if (*ptr_ == '@') {
			const char* i = std::find(ptr_, end_, '\n');
			return std::string(ptr_ + 1, trim_cr(ptr_ + 1, i));
		}
		else if (*ptr_ == '+') {
			const char* i = std::find(ptr_, end_, '\n');
			if (i == end_)
				return std::string();
			return std::string(i + 1, trim_cr(i + 1, end_));
		}
		else {
			const char* i = std::find(ptr_, end_, '\n');
			return std::string(ptr_, trim_cr(ptr_, i));
		}
	}
	virtual FastqTokenizer& operator++() override {
		if (*ptr_ == '@') {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
		}
		else if (*ptr_ == '+') {
			ptr_ = end_;
		}
		else {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
			if (ptr_ < end_ && *ptr_ != '+')
				throw MalformedFastqRecord();
		}
		return *this;
	}
	virtual std::pair<bool, std::string> read_record() override {
		if (ptr_ == end_)
			return { false,std::string() };
		if (*ptr_ != '@')
			throw MalformedFastqRecord();
		const char* l, *begin = ptr_;
		int n;
		std::tie(l, n) = find_newline(ptr_, end_);
		if (n == 0) throw MalformedFastqRecord();
		ptr_ = l;
		int64_t len = 0;
		for (;;) {
			std::tie(l, n) = find_newline(ptr_, end_);
			if (n == 0) throw MalformedFastqRecord();
			len += l - ptr_ - n;
			ptr_ = l;
			if (ptr_ == end_) throw MalformedFastqRecord();
			if (*ptr_ == '+') {
				std::tie(l, n) = find_newline(ptr_, end_);
				if (n == 0) throw MalformedFastqRecord();
				ptr_ = l;
				break;
			}
		}
		int64_t qlen = 0;
		for (;;) {
			std::tie(l, n) = find_newline(ptr_, end_);
			qlen += l - ptr_ - n;
			ptr_ = l;
			if (ptr_ == end_ || qlen >= len) break;
		}
		return { true, std::string(begin,ptr_) };
	}

	virtual std::pair<bool, std::string> read_record(TextInputFile& file) override {
		if (file.eof())
			return { false,std::string() };
		file.getline();
		if (file.eof() && file.line.empty())
			return { false,std::string() };
		const int64_t line = file.line_count;
		if (file.line.empty() || file.line[0] != '@' || file.eof())
			throw MalformedFastqRecord(line);
		std::string s(file.line + '\n');
		int64_t len = 0;
		for (;;) {
			file.getline();
			if (file.eof() || file.line.empty()) throw MalformedFastqRecord(line);
			if (file.line[0] == '+') {
				s += file.line;
				s += '\n';
				break;
			}
			len += file.line.length();
			s += file.line;
			s += '\n';
		}
		int64_t qlen = 0;
		for (;;) {
			file.getline();
			qlen += file.line.length();
			s += file.line;
			s += '\n';
			if (file.eof() || qlen >= len) break;
		}
		return { true,s };
	}
};

}}