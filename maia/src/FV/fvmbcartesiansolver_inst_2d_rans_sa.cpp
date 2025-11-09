// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#define ASSERT2 ASSERT
#include "fvmbcartesiansolverxd.cpp"


template class FvMbCartesianSolverXD<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
