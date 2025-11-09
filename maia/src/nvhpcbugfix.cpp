// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "compiler_config.h"

// TODO: The following is a WAR to avoid free() of invalid pointer on gpu
#if defined(MAIA_NVHPC_COMPILER) && defined(MAIA_PSTL)
void free(void*) {}
#endif
