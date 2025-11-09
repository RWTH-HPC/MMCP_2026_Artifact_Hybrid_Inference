// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVAPESOLVER2D_H_
#define FVAPESOLVER2D_H_

#include "fvcartesiansolverxd.h"
#include "fvcartesiansyseqnns.h"


class FvApeSolver2D : public FvCartesianSolverXD<2, FvSysEqnNS<2>> {
 public:
  FvApeSolver2D(MInt, MInt, MBool*, maia::grid::Proxy<2>& gridProxy_, Geometry<2>& geometry_, const MPI_Comm comm);
  void Ausm();
  MBool maxResidual(MInt /*mode*/);
  MFloat m_advectionVelocity[nDim];

 private:
  static constexpr const MInt nDim = 2;
};

// TODO labels:FV
// 1) Add the following methods (copy from original FV code and modify for APE):
// - Ausm
// - computePrimitiveVariables
// - computeConservativeVariables
// - viscousFlux
// - maxResidual (?)
// - lhsBnd (?)
//
// 2) Define new boundary conditions
//
// 3) Define new initial condition for first test
//
// Overall strategy:
// - Use FV solver "as-is", i.e., only modify as little as absolutely necessary
// - Continue to use variable names and nDim+2 variables: u, v, w, p are to be their perturbed
//   counterparts; set rho to be constant everywhere and just ignore it wherever possible
// - Primitive/conservative variables: use same variables everywhere

#endif // ifndef FVAPESOLVER2D_H_
