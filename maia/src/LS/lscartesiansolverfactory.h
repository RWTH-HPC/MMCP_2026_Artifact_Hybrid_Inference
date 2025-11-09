// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiansolver.h"
#include "lscartesiansolver.h"

template <MInt nDim>
class LsCartesianSolverFactory {
  using CartesianSolver = typename maia::CartesianSolver<nDim, LsCartesianSolver<nDim>>;
  using GridProxy = typename CartesianSolver::GridProxy;

 public:
  /**  \brief  Factory method for LsCartesianSolver
   *   \author Miro Gondrum
   *   \date   26.05.2023
   *   \return LsCartesianSolver
   */
  static std::unique_ptr<LsCartesianSolver<nDim>> create(MInt solverId_, const MBool* propertiesGroups,
                                                         GridProxy& gridProxy_, Geometry<nDim>& geometry_,
                                                         const MPI_Comm comm);
};
