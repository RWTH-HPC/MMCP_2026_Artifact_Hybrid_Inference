// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONEXACT_H_
#define DGBOUNDARYCONDITIONEXACT_H_

#include "dgcartesianboundarycondition.h"
#include "dgcartesiansyseqnlinearscalaradv.h"

/// \brief Boundary condition which imposes initial condition ("exact" boundary
///        conditions) at the domain boundaries.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014
///
/// Note: The exact BC uses the *internal* nodeVars for both the external as
/// well as the internal state. This is due to the fact that for direct-hybrid
/// simulations, i.e., where the nodeVars are not initialized by analytical
/// initial conditions but from a file, the external nodeVars are not set.
/// If this causes problems, another solution (probably more complex) has to be
/// found.
template <MInt nDim, class SysEqn>
class DgBcExact final : public DgBoundaryCondition<nDim, SysEqn> {
  // Typedefs
 public:
  using Base = DgBoundaryCondition<nDim, SysEqn>;
  using Base::begin;
  using Base::end;
  using Base::flux;
  using Base::surfaces;
  using Base::sysEqn;
  using typename Base::SolverType;

  // Methods
 public:
  DgBcExact(SolverType& solver_, MInt bcId) : Base(solver_, bcId) {}
  MString name() const override { return "exact"; }

  void apply(const MFloat time) override { maia::dg::bc::loop(&DgBcExact::applyAtSurface, this, begin(), end(), time); }

  void applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time));
};


template <MInt nDim, class SysEqn>
void DgBcExact<nDim, SysEqn>::applyAtSurface(const MInt surfaceId, const MFloat time) {
  MFloat* nodeVarsL = &surfaces().nodeVars(surfaceId, 0);
  MFloat* nodeVarsR = &surfaces().nodeVars(surfaceId, 1);
  MFloat* stateL = &surfaces().variables(surfaceId, 0);
  MFloat* stateR = &surfaces().variables(surfaceId, 1);

  // Calculate boundary state from initial condition
  // The boundary side is marked by an elment id of '-1'
  MFloat* boundaryState = (surfaces().internalSideId(surfaceId) == 1) ? stateL : stateR;
  MFloat* boundaryNodeVars = (surfaces().internalSideId(surfaceId) == 1) ? nodeVarsL : nodeVarsR;
  MFloat* innerNodeVars = (surfaces().internalSideId(surfaceId) == 1) ? nodeVarsR : nodeVarsL;

  // Iterate over all nodes and apply the initial condition at the boundary
  const MInt noNodes1D = surfaces().noNodes1D(surfaceId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MInt dirId = surfaces().orientation(surfaceId);
  MFloatTensor u(boundaryState, noNodes1D, noNodes1D3, SysEqn::noVars());
  // Set node vars to minimum 1 to avoid errors in Tensor
  // TODO labels:DG Check if this is really sensible
  MFloatTensor nodeVars(boundaryNodeVars, noNodes1D, noNodes1D3, std::max(SysEqn::noNodeVars(), 1));
  MFloatTensor x(&surfaces().nodeCoords(surfaceId), noNodes1D, noNodes1D3, nDim);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      sysEqn().calcInitialCondition(time, &x(i, j, 0), &nodeVars(i, j, 0), &u(i, j, 0));
    }
  }

  // Calculate Riemann flux
  // Use inner nodeVars for both sides (depending on the initial condition the
  // nodeVars are not set in calcInitialCondition)
  MFloat* f = flux(surfaceId);
  sysEqn().calcRiemann(innerNodeVars, innerNodeVars, stateL, stateR, noNodes1D, dirId, f);
}

#endif // DGBOUNDARYCONDITIONEXACT_H_
