// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgcartesiansolver.cpp"
#include "dgcartesiansyseqnacousticperturb.h"

// Prevent instantiation of boundary condition factory
// -> already done in dgboundaryconditionfactoryacousticperturb.cpp
extern template class DgBoundaryConditionFactory<3, DgSysEqnAcousticPerturb<3>>;

/// Explicit instantiation of 3D DG solver
template class DgCartesianSolver<3, DgSysEqnAcousticPerturb<3>>;
