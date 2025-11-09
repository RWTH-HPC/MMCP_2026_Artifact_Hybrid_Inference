// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianapesolver2d.h"

using namespace std;

FvApeSolver2D::FvApeSolver2D(MInt solverId, MInt noSpecies, MBool* propertiesGroups, maia::grid::Proxy<2>& gridProxy_,
                             Geometry<2>& geometry_, const MPI_Comm comm)
  : FvCartesianSolverXD<2, FvSysEqnNS<2>>(solverId, noSpecies, propertiesGroups, gridProxy_, geometry_, comm) {
  m_advectionVelocity[0] =
      Context::getSolverProperty<MFloat>("advectionVelocity", this->m_solverId, AT_, &m_advectionVelocity[0], 0);
  m_advectionVelocity[1] =
      Context::getSolverProperty<MFloat>("advectionVelocity", this->m_solverId, AT_, &m_advectionVelocity[1], 1);
  IF_CONSTEXPR(nDim == 3) {
    m_advectionVelocity[2] =
        Context::getSolverProperty<MFloat>("advectionVelocity", this->m_solverId, AT_, &m_advectionVelocity[2], 2);
  }
  TRACE();
}

void FvApeSolver2D::Ausm() {
  // TRACE();
  const MUint noVars = 2 + nDim;
  const MUint surfaceVarMemory = 2 * noVars;
  const MUint noSurfaces = a_noSurfaces();
  // const MUint noInternalSurfaces = m_bndrySurfacesOffset;
  // static constexpr MFloat epsilon = 0.000000000001; // pow( 10.0, -12.0 );

  MFloat* const RESTRICT fluxes = ALIGNED_MF(&a_surfaceFlux(0, 0));
  const MFloat* const RESTRICT surfaceVars = ALIGNED_F(&a_surfaceVariable(0, 0, 0));
  const MFloat* const RESTRICT area = ALIGNED_F(&a_surfaceArea(0));
  const MInt* const RESTRICT orientations = ALIGNED_I(&a_surfaceOrientation(0));
  // const MInt *const RESTRICT bndryCndIds = ALIGNED_I(srfcs[0].m_bndryCndId);

  // loop over all surfaces
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MUint srfcId = 0; srfcId < noSurfaces; ++srfcId) {
    //   const MUint orientation = orientations[ srfcId ];

    const MUint offset = srfcId * surfaceVarMemory;
    const MFloat* const RESTRICT leftVars = ALIGNED_F(surfaceVars + offset);
    const MFloat* const RESTRICT rightVars = ALIGNED_F(leftVars + noVars);

    // catch the primitive variables u, v, rho and p,
    const MFloat UL = leftVars[PV->U];
    const MFloat VL = leftVars[PV->V];
    const MFloat PL = leftVars[CV->RHO_E];

    const MFloat UR = rightVars[PV->U];
    const MFloat VR = rightVars[PV->V];
    const MFloat PR = rightVars[CV->RHO_E];

    // catch cell surface area A and advection velocity components
    const MFloat A = area[srfcId];

    const MFloat uBar = m_advectionVelocity[0];
    const MFloat vBar = m_advectionVelocity[1];

    const MUint fluxOffset = srfcId * noVars;
    MFloat* const RESTRICT flux = ALIGNED_MF(fluxes + fluxOffset);

    /* // APE fluxes are not yet implemented: here is a first draft of how they can be implemented
    in 2D

    //calc flux function, uses description of APE flux from Bauer 2011: Airframe noise prediction
    using a DG method
      //calc average density
          const MFloat rho0 = 1.0; //not yet implemented
      //calc speed of sound
          const MFloat c0 = 1.0 ; //not yet implemented
      //calc left fluxes
          const MFloat flux_u_l = orientations[srfcId]==0? uBar * UL + vBar * VL + rho0 * PL : 0.0
    ; const MFloat flux_v_l = orientations[srfcId]==1? uBar * UL + vBar * VL + rho0 * PL : 0.0  ;
          const MFloat flux_p_l = (rho0 * c0 * c0 * + PL) * (orientations[srfcId]==0? UL:VL) ;
      //calc right fluxes
          const MFloat flux_u_r = orientations[srfcId]==0? uBar * UR + vBar * VR + rho0 * PR : 0.0
    ; const MFloat flux_v_r = orientations[srfcId]==1? uBar * UR + vBar * VR + rho0 * PR : 0.0  ;
          const MFloat flux_p_r = (rho0 * c0 * c0 * + PR) * (orientations[srfcId]==0? UR:VR) ;

    //calc numerical flux using a Lax-Friedrichs-scheme as in Schlottke-Lakemper et al 2017
          const MFloat lambdaMax = max(abs(orientations[srfcId]==0? UL:VL),
    abs(orientations[srfcId]==0? UR:VR)); flux[FV->RHO_U] = 0.5 * (flux_u_l + flux_u_r + lambdaMax *
    (UL-UR) ); flux[FV->RHO_V] = 0.5 * (flux_v_l + flux_v_r + lambdaMax * (VL-VR) );

        //flux of perturbed pressure (for APE purposes we use the energy density for the perturbed
    pressure) flux[FV->RHO_E] = 0.5 * (flux_p_l + flux_p_r + lambdaMax * (PL-PR) );
    */

    // Temporary LAE: ic 1184 is used as a testcase for flux calculation with fluxes
    // for the 2D linear advection equation in the rhoU variable while the APE fluxes are still in
    // development. The other state variables(RHO_V, RHO_E, etc) are unused, since the LAE are only
    // for a one-dimensional time dependent state. Feel free to remove this as soon as APE fluxes
    // are implemented and well tested.
    if(m_initialCondition == 1184) {
      flux[FV->RHO_U] =
          orientations[srfcId] == 0 ? (uBar > 0.0 ? UL : UR) * A * uBar : (vBar > 0.0 ? UL : UR) * A * vBar;
      flux[FV->RHO_V] = 0.0;
      flux[FV->RHO_E] = 0.0;
      // This silences " ... is unused"-type warnings,
      // because these variables will be used for APE, but are not used for LAE.
      static_cast<void>(PL);
      static_cast<void>(PR);
      static_cast<void>(VL);
      static_cast<void>(VR);
    }
    flux[FV->RHO] = 0.0;
  }
}


// This function got hijacked to calculate an unnormed L_2 error for the domain for
// convergence/correctness analysis, it does not actually calculate a "residual".
// Rename/Restructure this if appropriate.
MBool FvApeSolver2D::maxResidual(MInt /*mode*/) {
  const MInt noVars = 2 + nDim;

  MFloat domainErrorSquared = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : domainErrorSquared)
#endif
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) { // TODO labels:FV iterate only over inner cells
    if(!a_hasProperty(cellId, SolverCell::IsOnCurrentMGLevel)) continue;
    if(a_isPeriodic(cellId)) continue;
    if(a_hasProperty(cellId, SolverCell::IsPeriodicWithRot)) continue;
    if(a_hasProperty(cellId, SolverCell::IsCutOff)) continue;
    if(a_hasProperty(cellId, SolverCell::IsInvalid)) continue;
    MInt gridcellId = cellId;
    if(a_hasProperty(cellId, SolverCell::IsSplitClone)) {
      gridcellId = m_splitChildToSplitCell.find(cellId)->second;
    }
    if(c_noChildren(gridcellId) > 0) continue;
    if(a_bndryId(cellId) > -1 && m_fvBndryCnd->m_bndryCells->a[a_bndryId(cellId)].m_linkedCellId > -1) continue;

    MFloat cellErrorSquared = 0.0;
    for(MInt var = 0; var < noVars; var++) {
      if(var != 0) {
        continue; // for LAE, only RHO_U is used
      }
      // This analytical solution matches the initial condition prescribed at the start of the
      // simulation.
      if(m_initialCondition == 1184) {
        MFloat pulse_radius = 0;
        for(MInt spaceId = 0; spaceId < nDim; spaceId++) {
          pulse_radius += std::pow(a_coordinate(cellId, spaceId) - 0.5 - m_advectionVelocity[spaceId] * m_time, 2.0);
        }
        MFloat analyticalSolutionCell = m_rhoVVInfinity[var] * (1.0 + std::exp(-pulse_radius * 20.0));
        cellErrorSquared += POW2(analyticalSolutionCell - a_variable(cellId, var));
      }
    }
    domainErrorSquared += cellErrorSquared;
  }
  m_log << "time: " << m_time << " domainerror: " << std::sqrt(domainErrorSquared) << std::endl;

  return true;
}
