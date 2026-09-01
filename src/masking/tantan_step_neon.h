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

static inline float32x4_t fmadd4(float32x4_t a, float32x4_t b, float32x4_t c) {
#if defined(__aarch64__)
	return vfmaq_f32(c, a, b);
#else
	return vmlaq_f32(c, a, b);
#endif
}

static inline float hsum4(float32x4_t v) {
#if defined(__aarch64__)
	return vaddvq_f32(v);
#else
	float32x2_t sum2 = vadd_f32(vget_low_f32(v), vget_high_f32(v));
	return vget_lane_f32(vpadd_f32(sum2, sum2), 0);
#endif
}

static inline float forward_step(float* __restrict f, const float* __restrict d, const float* __restrict e_seg,
	float& b, float f2f, float p_repeat_end, float b2b, float f_sum_prev)
{
	const float b_old = b;
	const float32x4_t vf2f = vdupq_n_f32(f2f);
	const float32x4_t vb_old = vdupq_n_f32(b_old);

	float f_sum_new = 0.0f;

	for (int off = 0; off < VEC_END; off += LANES) {
		const float32x4_t vd = vld1q_f32(d + off);
		const float32x4_t ve = vld1q_f32(e_seg + off);
		float32x4_t vf = vld1q_f32(f + off);
		vf = vmulq_f32(fmadd4(vf, vf2f, vmulq_f32(vb_old, vd)), ve);
		vst1q_f32(f + off, vf);
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
	float& b, float f2f, float p_repeat_end, float b2b)
{
	const float32x4_t vf2f = vdupq_n_f32(f2f);
	const float32x4_t vc = vdupq_n_f32(p_repeat_end * b);

	float tsum = 0.0f;

	for (int off = 0; off < VEC_END; off += LANES) {
		const float32x4_t ve = vld1q_f32(e_seg + off);
		const float32x4_t vd = vld1q_f32(d + off);
		float32x4_t vf = vmulq_f32(vld1q_f32(f + off), ve);
		const float32x4_t vt = vmulq_f32(vf, vd);
		vf = fmadd4(vf, vf2f, vc);
		vst1q_f32(f + off, vf);
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
	const float32x4_t vs = vdupq_n_f32(s);
	for (int off = 0; off < VEC_END; off += LANES)
		vst1q_f32(f + off, vmulq_f32(vld1q_f32(f + off), vs));
	for (int off = VEC_END; off < WINDOW; ++off)
		f[off] *= s;
}

static inline float sum_window(const float* __restrict f) {
	float32x4_t acc = vdupq_n_f32(0.0f);
	for (int off = 0; off < VEC_END; off += LANES)
		acc = vaddq_f32(acc, vld1q_f32(f + off));
	float s = hsum4(acc);
	for (int off = VEC_END; off < WINDOW; ++off)
		s += f[off];
	return s;
}

}}}