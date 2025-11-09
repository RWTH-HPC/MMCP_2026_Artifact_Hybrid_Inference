// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianapesolver3d.h"

using namespace std;

FvApeSolver3D::FvApeSolver3D(MInt solverId, MInt noSpecies, MBool* propertiesGroups, maia::grid::Proxy<3>& gridProxy_,
                             Geometry<3>& geometry_, const MPI_Comm comm)
  : FvCartesianSolverXD<3, FvSysEqnNS<3>>(solverId, noSpecies, propertiesGroups, gridProxy_, geometry_, comm) {
  TRACE();
  TERMM(1, "Not implemented");
}
