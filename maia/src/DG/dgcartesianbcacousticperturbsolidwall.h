// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONSOLIDWALL_H_
#define DGBOUNDARYCONDITIONSOLIDWALL_H_

#include "dgcartesianboundarycondition.h"

/// \brief Solid (slip) wall boundary condition.
///
/// \author Lev Liberson, Hsun-Jen Cheng
/// \date 2014
template <MInt nDim, class SysEqn, MBool slipWall = false>
class DgBcAcousticPerturbSolidWall final : public DgBoundaryCondition<nDim, SysEqn> {
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
  DgBcAcousticPerturbSolidWall(SolverType& solver_, MInt bcId) : Base(solver_, bcId) {}
  MString name() const override {
    return slipWall ? "solid wall (mean flow slip wall)" : "solid wall (mean flow no-slip wall)";
  }

  void apply(const MFloat time) override {
    maia::dg::bc::loop(&DgBcAcousticPerturbSolidWall::applyAtSurface, this, begin(), end(), time);
  }

  void applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time));
};


template <MInt nDim, class SysEqn, MBool slipWall>
void DgBcAcousticPerturbSolidWall<nDim, SysEqn, slipWall>::applyAtSurface(const MInt surfaceId,
                                                                          const MFloat NotUsed(time)) {
  // Storing the surface-information and define the inner state
  MFloat* stateL = &surfaces().variables(surfaceId, 0);
  MFloat* stateR = &surfaces().variables(surfaceId, 1);
  MFloat* innerState = (surfaces().internalSideId(surfaceId) == 0) ? stateL : stateR;
  MFloat* nodeVarsL = &surfaces().nodeVars(surfaceId, 0);
  MFloat* nodeVarsR = &surfaces().nodeVars(surfaceId, 1);
  MFloat* innerNodeVars = (surfaces().nghbrElementIds(surfaceId, 0) == -1) ? nodeVarsR : nodeVarsL;

  // Collect data for flux-computation
  const MInt noNodes1D = surfaces().noNodes1D(surfaceId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MInt dirId = surfaces().orientation(surfaceId);

  MFloatTensor variables(innerState, noNodes1D, noNodes1D3, SysEqn::noVars());
  MFloatTensor nodeVars(innerNodeVars, noNodes1D, noNodes1D3, SysEqn::noNodeVars());
  MFloatTensor f(flux(surfaceId), noNodes1D, noNodes1D3, SysEqn::noVars());

  // Computing the flux at all nodes considering zero-velocities at the wall
  //
  // Flux-components for boundary-surfaces (no-slip wall version):
  //       u'  v'  w'  p'
  // Fx = (p', 0,  0,  0)
  // Fy = (0,  p', 0,  0)
  // Fz = (0,  0,  p', 0)
  //
  // Reference: M. Bauer, Airframe noise prediction using a discontinuous
  //            Galerkin method, PhD thesis, p. 18, 2011
  //
  // Additionally, if slipWall == true, additional velocity terms are added to
  // account for non-zero mean flow in the wall-parallel direction.
  //
  // Flux-components for boundary-surfaces (slip-wall version):
  //       u'                v'                w'                p'
  // Fx = (p' + v'v0 + w'w0, 0,                0,                0)
  // Fy = (0,                p' + u'u0 + w'w0, 0,                0)
  // Fz = (0,                0,                p' + u'u0 + v'v0, 0)

  // Set the flux to zero
  f.set(0.0);

  // Set the x- or y- or z-component of the velocity to p', depending on the
  // surface orientation
  if(slipWall) {
    // If it is a slip wall, add pressure and additional velocity terms
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D3; j++) {
        // This following line is just so that std::inner_product can be used
        variables(i, j, dirId) = 0.0;

        // Set flux normal to wall (see above for equation)
        f(i, j, dirId) =
            std::inner_product(&variables(i, j, 0), &variables(i, j, nDim), &nodeVars(i, j, 0), variables(i, j, nDim));
      }
    }
  } else {
    // If it is a no-slip wall, just add the pressure term and implicity assume
    // that all velocity compontents are zero at the wall
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D3; j++) {
        // Set flux normal to wall (see above for equation)
        f(i, j, dirId) = variables(i, j, nDim);
      }
    }
  }
}

#endif // DGBOUNDARYCONDITIONSOLIDWALL_H_
