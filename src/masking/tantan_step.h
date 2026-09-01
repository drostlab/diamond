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
#include "util/simd.h"

namespace Util { namespace tantan {

static constexpr int WINDOW = 50;
static constexpr int WINDOW_ALIGN = 64;

}}

#if ARCH_ID == 3 && defined(WITH_AVX512_WIDE)
#include "tantan_step_avx512.h"
#elif ARCH_ID == 2
#include "tantan_step_avx2.h"
#elif defined(__ARM_NEON)
#include "tantan_step_neon.h"
#elif defined(__SSE2__)
#include "tantan_step_sse.h"
#else
#include "tantan_step_generic.h"
#endif