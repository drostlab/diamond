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

static constexpr int LANES = 4;
static constexpr int VEC_END = WINDOW - WINDOW % LANES;

static inline __m128 fmadd4(__m128 a, __m128 b, __m128 c) {
	return _mm_add_ps(_mm_mul_ps(a, b), c);
}

static inline float hsum4(__m128 v) {
	__m128 shuf = _mm_shuffle_ps(v, v, _MM_SHUFFLE(2, 3, 0, 1));
	__m128 sums = _mm_add_ps(v, shuf);
	shuf = _mm_shuffle_ps(sums, sums, _MM_SHUFFLE(1, 0, 3, 2));
	sums = _mm_add_ss(sums, shuf);
	return _mm_cvtss_f32(sums);
}

static inline float forward_step(float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b, float f_sum_prev) {
	const float b_old = b;
	const __m128 vf2f = _mm_set1_ps(f2f);
	const __m128 vb_old = _mm_set1_ps(b_old);

	float f_sum_new = 0.0f;

	for (int off = 0; off < VEC_END; off += LANES) {
		const __m128 vd = _mm_load_ps(d + off);
		const __m128 ve = _mm_loadu_ps(e_seg + off);
		__m128 vf = _mm_load_ps(f + off);
		vf = _mm_mul_ps(fmadd4(vf, vf2f, _mm_mul_ps(vb_old, vd)), ve);
		_mm_store_ps(f + off, vf);
		f_sum_new += hsum4(vf);
	}
	for (int off = VEC_END; off < WINDOW; ++off) {
		const float vf = (f[off] * f2f + b_old * d[off]) * e_seg[off];
		f[off] = vf;
		f_sum_new += vf;
	}

	b = b_old * b2b + f_sum_prev * p_repeat_end;
	return f_sum_new;
}

static inline float backward_step(float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b) {
	const __m128 vf2f = _mm_set1_ps(f2f);
	const __m128 vc = _mm_set1_ps(p_repeat_end * b);

	float tsum = 0.0f;

	for (int off = 0; off < VEC_END; off += LANES) {
		const __m128 ve = _mm_loadu_ps(e_seg + off);
		const __m128 vd = _mm_load_ps(d + off);
		__m128 vf = _mm_mul_ps(_mm_load_ps(f + off), ve);
		const __m128 vt = _mm_mul_ps(vf, vd);
		vf = fmadd4(vf, vf2f, vc);
		_mm_store_ps(f + off, vf);
		tsum += hsum4(vt);
	}
	for (int off = VEC_END; off < WINDOW; ++off) {
		const float vf = f[off] * e_seg[off];
		tsum += vf * d[off];
		f[off] = vf * f2f + p_repeat_end * b;
	}

	b = b2b * b + tsum;
	return tsum;
}

static inline void scale_window(float* __restrict f, float s) {
	const __m128 vs = _mm_set1_ps(s);
	for (int off = 0; off < VEC_END; off += LANES)
		_mm_store_ps(f + off, _mm_mul_ps(_mm_load_ps(f + off), vs));
	for (int off = VEC_END; off < WINDOW; ++off)
		f[off] *= s;
}

static inline float sum_window(const float* __restrict f) {
	__m128 acc = _mm_setzero_ps();
	for (int off = 0; off < VEC_END; off += LANES)
		acc = _mm_add_ps(acc, _mm_load_ps(f + off));
	float s = hsum4(acc);
	for (int off = VEC_END; off < WINDOW; ++off)
		s += f[off];
	return s;
}

}}}