// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IO/context.h"
#include "UTIL/debug.h"
#include "UTIL/parallelfor.h"
#include "enums.h"
#include "lbsolverdxqy.h"
#include "lbsrcterm.h"
#include "lbsrctermcontroller.h"

namespace maia::lb {

//---LbSrcTerm_nonnewtonian--------------------------------------------------
//-declaration
/** \brief  Class to handle non-newtonian fluids
 *  \author Moritz Waldmann
 *  \date   12.10.2023
 */
template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTerm_nonnewtonian : public LbSrcTerm<nDim, nDist, SysEqn> {
 protected:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;

  MFloat m_nu0 = F0;
  MFloat m_nuInf = F0;
  MFloat m_lambda = F0;
  MFloat m_n = F0;
  MFloat m_a = F0;
  MString m_model = "";
  void readProperties() override;

 public:
  //  static const LbRegSrcTerm<nDim, nDist, LbSrcTerm_nonnewtonian<nDim, nDist, SysEqn>> reg;

  LbSrcTerm_nonnewtonian(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver) { readProperties(); };

  void init() override{/*NOT USED*/};
  void apply_preCollision() override;
  void apply_postCollision() override{/*NOT USED*/};
  void apply_postPropagation() override{/*NOT USED*/};

  MFloat (LbSrcTerm_nonnewtonian::*m_getNuFromNonNewtonianModel)(const MFloat gamma_dot);
  MFloat powerLaw(const MFloat gamma_dot);
  MFloat carreau(const MFloat gamma_dot);
};

//-definiton
template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_nonnewtonian<nDim, nDist, SysEqn>::readProperties() {
  TRACE();
  const MInt solverId = m_solver->m_solverId;

  m_model = Context::getSolverProperty<MString>("nonNewtonianModel", solverId, AT_);

  switch(string2enum(m_model)) {
    case CARREAU: {
      m_nu0 = m_solver->m_nu;
      const MFloat nu0 = Context::getSolverProperty<MFloat>("nu0", solverId, AT_);
      const MFloat c_nu = nu0 / m_nu0;

      const MFloat nuInf = Context::getSolverProperty<MFloat>("nuInf", solverId, AT_);
      m_nuInf = nuInf / c_nu;

      const MFloat lambda = Context::getSolverProperty<MFloat>("lambda", solverId, AT_);

      m_lambda = lambda * m_solver->m_referenceLength / (m_solver->m_Ma * LBCS);

      m_n = Context::getSolverProperty<MFloat>("nonNewtonian_n", solverId, AT_);
      m_a = Context::getSolverProperty<MFloat>("nonNewtonian_a", solverId, AT_);
      m_getNuFromNonNewtonianModel = &LbSrcTerm_nonnewtonian::carreau;
      break;
    }
    case POWERLAW: {
      m_n = Context::getSolverProperty<MFloat>("nonNewtonian_n", solverId, AT_);

      m_solver->m_nu = pow(m_solver->m_Ma * LBCS, F2 - m_n) / m_solver->m_Re * pow(m_solver->m_referenceLength, m_n);
      m_nu0 = m_solver->m_nu;
      m_getNuFromNonNewtonianModel = &LbSrcTerm_nonnewtonian::powerLaw;
      break;
    }
    default: {
      mTerm(1, AT_, "This is not an implemented model!!!");
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTerm_nonnewtonian<nDim, nDist, SysEqn>::apply_preCollision() {
  TRACE();
  LbSolverDxQy<nDim, nDist, SysEqn>& s = *(m_solver); // alias for readability

  s.m_nu = s.m_Ma * LBCS / s.m_Re * s.m_referenceLength;
  maia::parallelFor(0, s.m_currentMaxNoCells, [&](MInt index) {
    const MInt pCellId = s.m_activeCellList[index];
    const MInt lvlDiff = s.maxLevel() - s.a_level(pCellId);
    if((globalTimeStep - 1) % IPOW2(lvlDiff) != 0) return;

    constexpr MInt nDimSqr = nDim * nDim;
    std::array<MFloat, nDimSqr> c{};
    const MFloat rho = s.a_variable(pCellId, s.PV->RHO);
#ifndef WAR_NVHPC_PSTL
    const MFloat* const u = &s.a_variable(pCellId, s.PV->U);
    const MFloat* const dist = &s.a_oldDistribution(pCellId, 0);
    s.sysEqn().calcMomentumFlux(rho, u, dist, c.data());
#else
    MFloat u[nDim] = {F0};
    MFloat dist[nDist] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = s.a_variable(pCellId, d);
    }
    for(MInt d = 0; d < nDist; d++) {
      dist[d] = s.a_oldDistribution(pCellId, d);
    }
    s.sysEqn().calcMomentumFlux(rho, u, dist, c.data());
#endif

    // Compute derivative of rate-of-strain tensor
    std::array<MFloat, nDimSqr> localS{};

    MFloat nu = F0;
    MInt noIterations = 0;
    MBool converged = false;
    MFloat nuOld = s.a_oldNu(pCellId);
    const MFloat F1Brho = F1 / rho;
    while(!converged && noIterations < 50) {
      MFloat omega = F2 / (F6 * nuOld * FFPOW2(lvlDiff) + F1);
      for(MInt d = 0; d < nDimSqr; d++) {
        localS[d] = -F3B2 * omega * F1Brho * c[d];
      }

      if(nDim == 3) {
        const MFloat volumeStrains = F1B3 * (localS[0] + localS[4] + localS[8]);
        localS[0] -= volumeStrains;
        localS[4] -= volumeStrains;
        localS[8] -= volumeStrains;
      }

      // Compute second invariant of the rate-of-strain tensor
      MFloat D = F0; // Second invariant of rate-of-strain tensor
      for(MInt d = 0; d < nDimSqr; d++) {
        D += localS[d] * localS[d];
      }
      const MFloat gamma_dot = sqrt(F2 * D);

      nu = (this->*m_getNuFromNonNewtonianModel)(gamma_dot);

      if(fabs(nuOld - nu) < 1e-10) converged = true;

      nuOld = nu;
      noIterations++;
    }

    s.a_nu(pCellId) = nu;
  });
}

template <MInt nDim, MInt nDist, class SysEqn>
MFloat LbSrcTerm_nonnewtonian<nDim, nDist, SysEqn>::powerLaw(const MFloat gamma_dot) {
  TRACE();

  // Compute local viscosity
  const MFloat exp = m_n - F1;
  const MFloat base = gamma_dot;
  const MFloat nu = m_nu0 * (pow(base, exp));

  return nu;
}

template <MInt nDim, MInt nDist, class SysEqn>
MFloat LbSrcTerm_nonnewtonian<nDim, nDist, SysEqn>::carreau(const MFloat gamma_dot) {
  TRACE();

  // Compute local viscosity
  const MFloat exp = (m_n - F1) / m_a;
  const MFloat base = pow(m_lambda * gamma_dot, m_a) + F1;
  const MFloat nu = m_nuInf + (m_nu0 - m_nuInf) * (pow(base, exp));

  return nu;
}

//-registration for LbSrcTermFactory class
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<2, 9, LbSysEqnIncompressible<2, 9>>>
    reg_nonNewtoniand2q9("LB_SRC_NONNEWTONIAN");
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<3, 19, LbSysEqnIncompressible<3, 19>>>
    reg_nonNewtoniand3q19("LB_SRC_NONNEWTONIAN");
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<3, 27, LbSysEqnIncompressible<3, 27>>>
    reg_nonNewtoniand3q27("LB_SRC_NONNEWTONIAN");
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<2, 9, LbSysEqnCompressible<2, 9>>>
    reg_nonNewtoniand2q9C("LB_SRC_NONNEWTONIAN");
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<3, 19, LbSysEqnCompressible<3, 19>>>
    reg_nonNewtoniand3q19C("LB_SRC_NONNEWTONIAN");
static const LbRegSrcTerm<LbSrcTerm_nonnewtonian<3, 27, LbSysEqnCompressible<3, 27>>>
    reg_nonNewtoniand3q27C("LB_SRC_NONNEWTONIAN");
} // namespace maia::lb
