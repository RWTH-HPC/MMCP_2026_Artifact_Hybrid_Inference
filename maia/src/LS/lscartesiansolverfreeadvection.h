// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSSOLVERFREEADVECTION_H_
#define LSSOLVERFREEADVECTION_H_

#include "lscartesiansolver.h"

/** \brief A LS solver solving for free advection
 *  // TODO labels:LS documentation LsCartesianSolverFreeAdvection
 */
template <MInt nDim>
class LsCartesianSolverFreeAdvection : public LsCartesianSolver<nDim> {
 public:
  using CartesianSolver = typename maia::CartesianSolver<nDim, LsCartesianSolverFreeAdvection>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;
  using Geom = Geometry<nDim>;

  // C'tor
  LsCartesianSolverFreeAdvection<nDim>(const MInt solverId, const MBool* propertiesGroups, GridProxy& gridProxy_,
                                       Geom& geometry_, const MPI_Comm comm);

 protected:
 private:
};

#endif // ifndef LSSOLVERFREEADVECTION_H_
