// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONSTRAIGHTDUCTEXIT_H_
#define DGBOUNDARYCONDITIONSTRAIGHTDUCTEXIT_H_

#include "dgcartesianbcacousticperturbsolidwall.h"
#include "dgcartesianboundarycondition.h"

/* \brief Outlet of a straight duct in a solid wall (2D).
 *
 * The extend of the duct is determined by the property ductPosition.
 * It gives the interval on the boundary where the boundary condition is
 * applied, the remaining part of this boundary is treated as a solid wall.
 */
template <MInt nDim, class SysEqn>
class DgBcAcousticPerturbStraightDuctExit final : public DgBoundaryCondition<nDim, SysEqn> {
  // Typedefs
 public:
  using Base = DgBoundaryCondition<nDim, SysEqn>;
  using Base::begin;
  using Base::end;
  using Base::flux;
  using Base::solver;
  using Base::surfaces;
  using Base::sysEqn;
  using typename Base::SolverType;

  // Methods
  DgBcAcousticPerturbStraightDuctExit(SolverType& solver_, MInt bcId) : Base(solver_, bcId), m_wallBc(solver_, bcId) {}
  MString name() const override { return "Straight duct"; }

  void init() override;

  void apply(const MFloat time) override {
    maia::dg::bc::loop(&DgBcAcousticPerturbStraightDuctExit::applyAtSurface, this, begin(), end(), time);
  }

  void applyAtSurface(const MInt surfaceId, const MFloat time);

 private:
  // Wall boundary condition
  DgBcAcousticPerturbSolidWall<nDim, SysEqn> m_wallBc;

  // Position of duct exit: interval [p1,p2]
  MFloat m_ductPosition[2]{};
};

/* \brief Initialize boundary condition.
 *
 * Reads the property ductPosition which determines the position of the duct
 * exit.
 */
template <MInt nDim, class SysEqn>
void DgBcAcousticPerturbStraightDuctExit<nDim, SysEqn>::init() {
  // Check that this boundary condition is not accidentally used in 3D
  IF_CONSTEXPR(nDim != 2) { TERMM(1, "This boundary condition is only usable for 2D simulations"); }

  /*! \property
    \page propertyPageDG DG
    \section ductPosition
    <code>MFloat DgBcAcousticPerturbStraightDuctExit::m_ductPosition</code> \n
    Determines the interval on the boundary in which the duct exit boundary
    condition is placed.\n
    Two values are needed which are the lower and upper bound of the duct along
    the coordinate direction of the boundary.\n
    The boundary condition is applied to a surface if its center for this
    coordinate direction lies in the given interval.
    Keywords: <i>DISCONTINUOUS_GALERKIN</i>
  */
  for(MInt i = 0; i < 2; i++) {
    m_ductPosition[i] = Context::getSolverProperty<MFloat>("ductPosition", solver().solverId(), AT_, i);
  }
}


template <MInt nDim, class SysEqn>
void DgBcAcousticPerturbStraightDuctExit<nDim, SysEqn>::applyAtSurface(const MInt surfaceId, const MFloat time) {
  const MFloat* const coordinates = &surfaces().coords(surfaceId, 0);
  const MInt sOrientation = surfaces().orientation(surfaceId);

  const MFloat pos = coordinates[(sOrientation + 1) % 2];
  const MFloat lb = m_ductPosition[0];
  const MFloat ub = m_ductPosition[1];

  // Check position of surface center
  if(lb > pos || pos > ub) {
    // Apply wall boundary condition
    m_wallBc.applyAtSurface(surfaceId, time);
  } else {
    // Calculate flux with the inner state of a surface and boundary state
    const MInt noNodes1D = surfaces().noNodes1D(surfaceId);

    const MFloat* surfaceStateL = &surfaces().variables(surfaceId, 0);
    const MFloat* surfaceStateR = &surfaces().variables(surfaceId, 1);
    const MFloat* outerState = (surfaces().nghbrElementIds(surfaceId, 0) == -1) ? surfaceStateL : surfaceStateR;
    MFloatTensor oS(const_cast<MFloat*>(outerState), noNodes1D, SysEqn::noVars());

    // Check the orientation of the Duct-exit and set the outer state
    if(sOrientation == 1) {
      for(MInt n = 0; n < noNodes1D; n++) {
        // for y-orientation x-velocity is zero
        oS(n, 0) = 0.0;
        oS(n, 1) = sin(2 * PI * time);
        oS(n, 2) = sin(2 * PI * time);
      }
    } else {
      for(MInt n = 0.0; n < noNodes1D; n++) {
        // for x orientation y-velocity is zero
        oS(n, 0) = sin(2 * PI * time);
        oS(n, 1) = 0;
        oS(n, 2) = sin(2 * PI * time);
      }
    }

    // Call regular Riemann solver
    const MFloat* nodeVarsL = &surfaces().nodeVars(surfaceId, 0);
    const MFloat* nodeVarsR = &surfaces().nodeVars(surfaceId, 1);
    const MFloat* innerNodeVars = (surfaces().nghbrElementIds(surfaceId, 0) == -1) ? nodeVarsR : nodeVarsL;
    sysEqn().calcRiemann(innerNodeVars, innerNodeVars, surfaceStateL, surfaceStateR, noNodes1D, sOrientation,
                         flux(surfaceId));
  }
}

#endif // DGBOUNDARYCONDITIONSTRAIGHTDUCTEXIT_H_
