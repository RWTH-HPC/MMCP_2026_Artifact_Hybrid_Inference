// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSSOLVERSEMILAGRANGIAN_H_
#define LSSOLVERSEMILAGRANGIAN_H_

#include "lscartesiansolver.h"

/** \brief A LS solver based on the semi-lagrangian method
 *  // TODO labels:LS documentation LsCartesianSolverSemiLagrangian
 */
template <MInt nDim>
class LsCartesianSolverSemiLagrangian : public LsCartesianSolver<nDim> {
 public:
  using CartesianSolver = typename maia::CartesianSolver<nDim, LsCartesianSolverSemiLagrangian>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;
  using Geom = Geometry<nDim>;

  // C'tor
  LsCartesianSolverSemiLagrangian<nDim>(const MInt solverId, const MBool* propertiesGroups, GridProxy& gridProxy_,
                                        Geom& geometry_, const MPI_Comm comm);

 protected:
 private:
};

#endif // ifndef LSSOLVERSEMILAGRANGIAN_H_
