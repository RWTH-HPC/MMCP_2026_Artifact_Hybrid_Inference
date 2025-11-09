// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "couplerlbfv.h"
#include "coupling.h"

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>::CouplerLbFv(const MInt couplingId, LbSolver* lb, FvCartesianSolver* fv)
  : Coupling(couplingId),
    CouplingFv<nDim, SysEqnFv>(couplingId, fv),
    CouplingLB<nDim, nDist, SysEqnLb>(couplingId, lb) {
  initData();
  readProperties();

  // calculation of cfl to sync the timesteps of the FV to the LB
  const MFloat newCFL = (1.0 + fvSolver().a_Ma()) / sqrt(3.0);
  fvSolver().a_cfl() = newCFL;

  if(fvSolver().domainId() == 0) {
    std::cerr << "CFL number calculated and set to: " << newCFL << std::endl;
    m_log << "CFL number calculated and set to: " << newCFL << std::endl;
  }

  initConversionFactors();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>::initConversionFactors() {
  const MFloat lbL2FvL = fvSolver().m_referenceLength / fvSolver().c_cellLengthAtLevel(fvSolver().maxLevel());
  conversionLbFv.length = lbL2FvL;
  conversionFvLb.length = 1.0 / lbL2FvL;

  const MFloat lbU2FvU = sqrt(3.0 / (1.0 + (fvSolver().m_gamma - 1.0) / 2.0 * fvSolver().a_Ma() * fvSolver().a_Ma()));
  conversionLbFv.velocity = lbU2FvU;
  conversionFvLb.velocity = 1.0 / lbU2FvU;

  MFloat lbP2FvP = pow(fvSolver().m_TInfinity, fvSolver().m_gamma / (fvSolver().m_gamma - 1.0));
  if(std::isnan(lbP2FvP)) {
    const MFloat TInf = fvSolver().m_initialCondition == 465 || fvSolver().m_initialCondition == 9465
                            ? 1.0
                            : 1.0 / (1.0 + 0.5 * (fvSolver().m_gamma - 1.0) * POW2(fvSolver().m_Ma));
    lbP2FvP = pow(TInf, fvSolver().m_gamma / (fvSolver().m_gamma - 1.0));
  }
  conversionLbFv.pressure = lbP2FvP;
  conversionFvLb.pressure = 1.0 / lbP2FvP;

  MFloat lbNu2FvNu = fvSolver().c_cellLengthAtLevel(fvSolver().maxLevel()) * sqrt(fvSolver().m_TInfinity * 3.0)
                     * fvSolver().sysEqn().m_Re0;
  if(std::isnan(lbNu2FvNu)) {
    const MFloat TInf = fvSolver().m_initialCondition == 465 || fvSolver().m_initialCondition == 9465
                            ? 1.0
                            : 1.0 / (1.0 + 0.5 * (fvSolver().m_gamma - 1.0) * POW2(fvSolver().m_Ma));
    lbNu2FvNu = fvSolver().c_cellLengthAtLevel(fvSolver().maxLevel()) * sqrt(TInf * 3.0) * fvSolver().sysEqn().m_Re0;
  }
  conversionLbFv.viscosity = lbNu2FvNu;
  conversionFvLb.viscosity = 1.0 / lbNu2FvNu;
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>::initData() {
  m_fvSolverId = fvSolver().solverId();
  m_lbSolverId = lbSolver().solverId();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>::checkProperties() {}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>::readProperties() {}

template class CouplerLbFv<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>, FvSysEqnNS<3>>;
template class CouplerLbFv<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>, FvSysEqnNS<3>>;
template class CouplerLbFv<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>, FvSysEqnEEGas<3>>;
template class CouplerLbFv<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>, FvSysEqnEEGas<3>>;
