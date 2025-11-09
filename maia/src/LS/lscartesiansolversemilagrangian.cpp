// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lscartesiansolversemilagrangian.h"

template <MInt nDim>
LsCartesianSolverSemiLagrangian<nDim>::LsCartesianSolverSemiLagrangian(const MInt solverId,
                                                                       const MBool* propertiesGroups,
                                                                       GridProxy& gridProxy, Geometry<nDim>& geometry,
                                                                       const MPI_Comm comm)
  : LsCartesianSolver<nDim>(solverId, propertiesGroups, gridProxy, geometry, comm) {}

template class LsCartesianSolverSemiLagrangian<2>;
template class LsCartesianSolverSemiLagrangian<3>;
