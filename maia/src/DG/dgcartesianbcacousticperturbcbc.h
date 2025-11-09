// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONACOUSTICPERTURBCBC_H_
#define DGBOUNDARYCONDITIONACOUSTICPERTURBCBC_H_

#include "dgcartesianboundarycondition.h"

template <MInt nDim>
class DgBcAcousticPerturbCBC final : public DgBoundaryCondition<nDim, DgSysEqnAcousticPerturb<nDim>> {
  // Typedefs
 public:
  using Base = DgBoundaryCondition<nDim, DgSysEqnAcousticPerturb<nDim>>;
  using Base::begin;
  using Base::end;
  using Base::flux;
  using Base::id;
  using Base::surfaces;
  using Base::sysEqn;
  using typename Base::SolverType;
  using typename Base::SysEqn;

  // Methods
  DgBcAcousticPerturbCBC(SolverType& solver_, MInt bcId) : Base(solver_, bcId) {}
  MString name() const override { return "CBC"; }

  void apply(const MFloat time) override {
    maia::dg::bc::loop(&DgBcAcousticPerturbCBC::applyAtSurface, this, begin(), end(), time);
  }

  void applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time));
};


template <MInt nDim>
void DgBcAcousticPerturbCBC<nDim>::applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time)) {
  // Storing the surface-information and define the inner state
  MFloat* nodeVarsL = &surfaces().nodeVars(surfaceId, 0);
  MFloat* nodeVarsR = &surfaces().nodeVars(surfaceId, 1);
  MFloat* stateL = &surfaces().variables(surfaceId, 0);
  MFloat* stateR = &surfaces().variables(surfaceId, 1);
  // TODO labels:DG,totest This might not work for all cases
  MFloat* innerNodeVars = (surfaces().internalSideId(surfaceId) != 1) ? nodeVarsL : nodeVarsR;
  MFloat* boundaryNodeVars = (surfaces().internalSideId(surfaceId) != 1) ? nodeVarsR : nodeVarsL;
  MFloat* innerState = (surfaces().nghbrElementIds(surfaceId, 0) == -1) ? stateR : stateL;
  MFloat* boundaryState = (surfaces().nghbrElementIds(surfaceId, 0) == -1) ? stateL : stateR;

  // Collect data for flux-computation
  const MInt dirId = surfaces().orientation(surfaceId);
  const MInt noNodes1D = surfaces().noNodes1D(surfaceId);
  const MInt noNodesXD = surfaces().noNodesXD(surfaceId);
  MFloat* f = flux(surfaceId);

  switch(id()) {
    // CBC OutFlow
    case 301: {
      std::fill_n(&boundaryState[0], noNodesXD * Base::SysEqn::noVars(), 0.0);
      std::fill_n(&boundaryNodeVars[0], noNodesXD * Base::SysEqn::noNodeVars(), 0.0);
      // Set default node vars of boundary state that are not zero
      for(MInt i = 0; i < noNodesXD; i++) {
        MFloat* const nodeVars = &boundaryNodeVars[i * Base::SysEqn::noNodeVars()];
        // Set default density
        nodeVars[Base::SysEqn::CV::RHO0] = 1.0;
        // Set default speed of sound
        nodeVars[Base::SysEqn::CV::C0] = 1.0;
      }
      sysEqn().calcRiemannRoe(nodeVarsL, nodeVarsR, stateL, stateR, noNodes1D, dirId, f);
      break;
    }

    // "periodic" extension of domain for 1D test
    case 304: {
      sysEqn().calcRiemannRoe(innerNodeVars, innerNodeVars, innerState, innerState, noNodes1D, dirId, f);
      break;
    }

    default: {
      TERMM(1, "Bad boundary condition id for CBC boundary conditions");
      break;
    }
  }
}

#endif // DGBOUNDARYCONDITIONACOUSTICPERTURBCBC_H_
