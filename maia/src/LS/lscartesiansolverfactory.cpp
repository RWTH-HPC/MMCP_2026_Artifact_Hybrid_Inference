// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lscartesiansolverfactory.h"
#include "lscartesiansolverfreeadvection.h"
#include "lscartesiansolversemilagrangian.h"

template <MInt nDim>
std::unique_ptr<LsCartesianSolver<nDim>>
LsCartesianSolverFactory<nDim>::create(const MInt solverId, const MBool* propertiesGroups, GridProxy& gridProxy,
                                       Geometry<nDim>& geometry, const MPI_Comm comm) {
  MInt solverMethodInt = -1;
  if(Context::propertyExists("solverMethod", solverId)) {
    solverMethodInt = string2enum(Context::getSolverProperty<MString>("solverMethod", solverId, AT_));
  }
  // TODO labels:LS More relevant input to trigger correct built of desired LS solver version is necessary
  const MBool isSemiLagrange =
      (solverMethodInt == MAIA_SEMI_LAGRANGE_LEVELSET || solverMethodInt == MAIA_RUNGE_KUTTA_MB_SEMI_LAGRANGE_LEVELSET
       || solverMethodInt == MAIA_SEMI_LAGRANGE_LEVELSET_LB);
  const MBool isFreeAdvection = (solverMethodInt == MAIA_LEVELSET_SURFACE);
  if(isSemiLagrange) {
    return std::make_unique<LsCartesianSolverSemiLagrangian<nDim>>(solverId, propertiesGroups, gridProxy, geometry,
                                                                   comm);
  } else if(isFreeAdvection) {
    return std::make_unique<LsCartesianSolverFreeAdvection<nDim>>(solverId, propertiesGroups, gridProxy, geometry,
                                                                  comm);
  } else {
    return std::make_unique<LsCartesianSolver<nDim>>(solverId, propertiesGroups, gridProxy, geometry, comm);
  }
}

template class LsCartesianSolverFactory<2>;
template class LsCartesianSolverFactory<3>;
