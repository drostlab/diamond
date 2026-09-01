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

/* Based on tantan by Martin C. Frith. See:
http://cbrc3.cbrc.jp/~martin/tantan/
A new repeat-masking method enables specific detection of homologous sequences, MC Frith, Nucleic Acids Research 2011 39(4):e23. */

#include <vector>
#include "../basic/value.h"
#include "def.h"
#include "../util/simd/dispatch.h"
#include "tantan_step.h"
#include "masking.h"

using std::fill;

namespace Util { namespace tantan { namespace DISPATCH_ARCH {

Mask::Ranges mask(
	Letter *seq,
	int len,
	const float **likelihood_ratio_matrix,
	float p_repeat,
	float p_repeat_end,
	float repeat_growth,
	float p_mask,
	int mask_mode)
{
	constexpr int RESERVE = 50000;

	Mask::Ranges ranges;
	if (len == 0) return ranges;

	alignas(WINDOW_ALIGN) float f[WINDOW];
	alignas(WINDOW_ALIGN) float d[WINDOW];

	const float b2b  = 1.0f - p_repeat;
	const float f2f  = 1.0f - p_repeat_end;
	const float b2f0 = p_repeat * (1.0f - repeat_growth) / (1.0f - std::pow(repeat_growth, (float)WINDOW));

	d[WINDOW - 1] = b2f0;
	for (int i = WINDOW - 2; i >= 0; --i)
		d[i] = d[i + 1] * repeat_growth;

	fill(f, f + WINDOW, 0.0f);
	
	thread_local std::vector<float> pb;
	thread_local std::vector<float> scale;
	thread_local std::vector<float> e[AMINO_ACID_COUNT];

	const int pb_cap = std::max(len, RESERVE);
	const int sc_len = std::max((len - 1) / 16 + 1, (RESERVE - 1) / 16 + 1);

	pb.resize(pb_cap);
	scale.resize(sc_len);

	for (size_t aa = 0; aa < (size_t)AMINO_ACID_COUNT; ++aa) {
		auto &E = e[aa];
		const int need = std::max(RESERVE, len + WINDOW);
		E.resize(need);

		const float* L = likelihood_ratio_matrix[aa];
		float* p = &E[len - 1];
		for (int j = 0; j < len; ++j) {
			const uint8_t idx = static_cast<uint8_t>(letter_mask(seq[j]));
			*(p--) = L[(size_t)idx];
		}
		fill(E.data() + len, E.data() + len + WINDOW, (float)0.0);
	}

	float b = 1.0f;
	float f_sum = 0.0f;

	for (int i = 0; i < len; ++i) {
		const uint8_t ltr = static_cast<uint8_t>(letter_mask(seq[i]));
		const float* e_seg = &e[ltr][len - i];
		float f_sum_new = forward_step(f, d, e_seg, b, f2f, p_repeat_end, b2b, f_sum);
		f_sum = f_sum_new;
		if ((i & 15) == 15) {
			const float s = 1.0f / b;
			scale[(size_t)i / 16] = s;
			b *= s;
			scale_window(f, s);
			f_sum *= s;
		}
		pb[i] = b;
	}

	const float z = b * b2b + sum_window(f) * p_repeat_end;
	const float zinv = 1.0f / z;

	b = b2b;
	fill(f, f + WINDOW, p_repeat_end);
	const Letter mask = value_traits.mask_char;
	
	for (int i = len - 1; i >= 0; --i) {
		const float pf = 1.0f - (pb[i] * b * zinv);

		if ((i & 15) == 15) {
			const float s = scale[(size_t)i / 16];
			b *= s;
			scale_window(f, s);
		}

		const uint8_t ltr = static_cast<uint8_t>(letter_mask(seq[i]));
		const float* e_seg = &e[ltr][len - i];
		backward_step(f, d, e_seg, b, f2f, p_repeat_end, b2b);

		if (pf >= p_mask) {
			if (mask_mode == 1)
				seq[i] = mask;
			else if (mask_mode == 2)
				seq[i] |= Masking::bit_mask;
			ranges.push_front(i);
		}
	}

	return ranges;
}

}

DISPATCH_8(Mask::Ranges, mask, Letter*, seq, int, len, const float**, likelihood_ratio_matrix, float, p_repeat, float, p_repeat_end, float, repeat_decay, float, p_mask, int, mask_mode)

}}