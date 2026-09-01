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
#include "../memory/alignment.h"
#include "../memory/memory_resource.h"

template<typename T>
struct MemBuffer {

	enum { ALIGN = 32 };

	typedef T value_type;

	MemBuffer(std::pmr::memory_resource* resource = nullptr):
		resource_(resource),
		data_(nullptr),
		size_(0),
		alloc_size_(0)
	{}

	MemBuffer(size_t n, std::pmr::memory_resource* resource = nullptr):
		resource_(resource),
		data_((T*)alloc(n)),
		size_(n),
		alloc_size_(n)
	{
	}

	MemBuffer(const MemBuffer&) = delete;
	MemBuffer& operator=(const MemBuffer&) = delete;

	~MemBuffer() {
		free_data();
	}

	void resize(size_t n) {
		if (alloc_size_ < n) {
			free_data();
			data_ = (T*)alloc(n);
			alloc_size_ = n;
		}
		size_ = n;
	}

	size_t size() const {
		return size_;
	}

	T* begin() {
		return data_;
	}

	T* end() {
		return data_ + size_;
	}

	const T* begin() const {
		return data_;
	}

	const T* end() const {
		return data_ + size_;
	}

	T& operator[](size_t i) {
		return data_[i];
	}

	const T& operator[](size_t i) const {
		return data_[i];
	}

private:

	void* alloc(size_t n) {
#ifdef HAVE_MEMORY_RESOURCE
		if (resource_)
			return resource_->allocate(n * sizeof(T), ALIGN);
#endif
		return Util::Memory::aligned_malloc(n * sizeof(T), ALIGN);
	}

	void free_data() {
		if (data_ == nullptr)
			return;
#ifdef HAVE_MEMORY_RESOURCE
		if (resource_)
			resource_->deallocate(data_, alloc_size_ * sizeof(T), ALIGN);
		else
#endif
			Util::Memory::aligned_free(data_);
		data_ = nullptr;
		size_ = alloc_size_ = 0;
	}

	std::pmr::memory_resource* resource_;
	T *data_;
	size_t size_, alloc_size_;

};