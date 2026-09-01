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
#include "tantan_step.h"

namespace Util { namespace tantan { namespace DISPATCH_ARCH {

static inline float forward_step(float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b, float f_sum_prev)
{
	const float b_old = b;
	float f_sum_new = 0.0f;

	for (int off = 0; off < WINDOW; ++off) {
		const float vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
		f[off] = vf;
		f_sum_new += vf;
	}

	b = b_old * b2b + f_sum_prev * p_repeat_end;
	return f_sum_new;
}

static inline float backward_step(float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b) {
	float tsum = 0.0f;

	for (int off = 0; off < WINDOW; ++off) {
		const float vf = f[off] * e_seg[off];
		tsum += vf * d[off];
		f[off] = vf * f2f + p_repeat_end * b;
	}

	b = b2b * b + tsum;
	return tsum;
}

static inline void scale_window(float* __restrict f, float s) {
	for (int off = 0; off < WINDOW; ++off)
		f[off] *= s;
}

static inline float sum_window(const float* __restrict f) {
	float s = 0.0f;
	for (int off = 0; off < WINDOW; ++off)
		s += f[off];
	return s;
}

}}}