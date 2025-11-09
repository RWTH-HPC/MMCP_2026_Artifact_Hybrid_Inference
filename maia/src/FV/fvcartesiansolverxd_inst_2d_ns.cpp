// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesiansyseqnns.h"
#include "fvcartesiansolverxd.cpp"

// Explicit instantiation for 2D
template class FvCartesianSolverXD<2, FvSysEqnNS<2>>;
