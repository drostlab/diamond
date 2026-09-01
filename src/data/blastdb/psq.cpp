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

#include "volume.h"

using std::vector;
using std::out_of_range;
using std::runtime_error;

void decode_protein_sequence(const char* data, size_t len, vector<Letter>& out)
{
    size_t begin = 0, end = len;
    if (begin < end && data[begin] == '\0')
        ++begin;
    if (end > begin && data[end - 1] == '\0')
        --end;

    out.clear();
    for (size_t i = begin; i < end; ++i) {
        const uint8_t aa = data[i];
        if (aa == '\0')
            throw runtime_error("Unexpected null terminator in sequence data");
        if (aa >= sizeof(NCBI_TO_STD) / sizeof(NCBI_TO_STD[0]))
            throw runtime_error("Invalid amino acid code in sequence data");
        out.push_back(NCBI_TO_STD[aa]);
    }
}

vector<Letter> decode_protein_sequence(const char* data, size_t len)
{
    vector<Letter> decoded;
    decoded.reserve(len);
    decode_protein_sequence(data, len, decoded);
    return decoded;
}

vector<Letter> BlastVolume::sequence(uint32_t oid)
{
    if (oid >= index_.num_oids) {
        throw out_of_range("OID exceeds number of sequences in volume");
    }

    const uint32_t start = index_.sequence_index[oid];
    const uint32_t end = index_.sequence_index[oid + 1];

    if (!index_.is_protein) {
        throw runtime_error("Nucleotide sequence decoding is not supported yet");
    }
    if (oid != seq_ptr_)
        psq_mapping_.seek(start, SEEK_SET);
    seq_ptr_ = oid + 1;
    return decode_protein_sequence(psq_mapping_.read_bytes(end - start), end - start);
}

vector<char> BlastVolume::raw_sequence(uint32_t count) {
    const size_t n = index_.sequence_index[seq_ptr_ + count] - index_.sequence_index[seq_ptr_];
    vector<char> v(n);
    psq_mapping_.read(v.data(), n);
    seq_ptr_ += count;
    return v;
}

Loc BlastVolume::length(uint32_t oid) {
    const uint32_t start = index_.sequence_index[oid];
    const uint32_t end = index_.sequence_index[oid + 1];
    return end - start - 1;
}