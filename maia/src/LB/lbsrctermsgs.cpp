// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "IO/context.h"
#include "UTIL/debug.h"
#include "UTIL/parallelfor.h"
#include "lbsolverdxqy.h"
#include "lbsrcterm.h"
#include "lbsrctermcontroller.h"

namespace maia::lb {

//---LbSrcTerm_smagorinsky--------------------------------------------------
//-declaration
/** \brief  Class to handle sub-grid scale model (Smagorinsky)
 *  \author Miro Gondrum
 *  \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn, MBool bubble = false>
class LbSrcTerm_smagorinsky : public LbSrcTerm<nDim, nDist, SysEqn> {
 protected:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;
  MFloat m_Cs = 0.1;           ///< Smagorinsky constant
  const MFloat m_deltaX = 1.0; ///< Filter width

  void readProperties() override;

 public:
  //  static const LbRegSrcTerm<nDim, nDist, LbSrcTerm_smagorinsky<nDim, nDist, SysEqn>> reg;

  LbSrcTerm_smagorinsky(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver) { readProperties(); };

  void init() override{/*NOT USED*/};
  void apply_preCollision() override;
  void apply_postCollision() override{/*NOT USED*/};
  void apply_postPropagation() override{/*NOT USED*/};
};

//-definiton
template <MInt nDim, MInt nDist, class SysEqn, MBool bubble>
void LbSrcTerm_smagorinsky<nDim, nDist, SysEqn, bubble>::readProperties() {
  TRACE();
  const MInt solverId = m_solver->m_solverId;
  if(Context::propertyExists("smagorinskyConstant")) {
    m_Cs = Context::getSolverProperty<MFloat>("smagorinskyConstant", solverId, AT_);
  }
}

template <MInt nDim, MInt nDist, class SysEqn, MBool bubble>
void LbSrcTerm_smagorinsky<nDim, nDist, SysEqn, bubble>::apply_preCollision() {
  TRACE();
  LbSolverDxQy<nDim, nDist, SysEqn>* const s = m_solver; // alias for readability

  s->m_nu = s->m_Ma * LBCS / s->m_Re * s->m_referenceLength;
  maia::parallelFor<false>(0, s->m_currentMaxNoCells, [=](MInt index) {
    const MInt pCellId = s->m_activeCellList[index];
    const MInt lvlDiff = s->maxLevel() - s->a_level(pCellId);
    if((globalTimeStep - 1) % IPOW2(lvlDiff) != 0) return;
    // TODO labels:LB daniell: This is critical for setting a_oldNu ?! But without
    // test case 3D_LB-FV_Euler-Euler_bubble_sphere fails
    s->a_nu(pCellId) = s->m_nu;
    // Ref: Hou, Sterling, Chen, and Doolen 1994
    // Also worth reading: Siebler, Krafczyk, Freudiger, and Geier 2011
    constexpr MInt nDimSqr = nDim * nDim;
    MFloat c[nDimSqr]{};
    const MFloat rho = s->a_variable(pCellId, s->PV->RHO);
#ifndef WAR_NVHPC_PSTL
    const MFloat* const u = &s->a_variable(pCellId, s->PV->U);
    const MFloat* const dist = &s->a_oldDistribution(pCellId, 0);
    s->sysEqn().calcMomentumFlux(rho, u, dist, c);
#else
    MFloat u[nDim] = {F0};
    MFloat dist[nDist] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = s->a_variable(pCellId, d);
    }
    for(MInt d = 0; d < nDist; d++) {
      dist[d] = s->a_oldDistribution(pCellId, d);
    }
    s->sysEqn().calcMomentumFlux(rho, u, dist, c);
#endif
    const MFloat Q = sqrt(2.0 * std::inner_product(&c[0], &c[nDimSqr], &c[0], .0));
    MFloat nu_BIT = 0.0;
    if constexpr(bubble) {
      if(s->m_isEELiquid) {
        // BIT (bubble induced turbulence) see Mohammadi 19 / Sato 81 / Milelli 01
        const MFloat u_d = sqrt(POW2(s->a_uOtherPhase(pCellId, 0) - u[0]) + POW2(s->a_uOtherPhase(pCellId, 1) - u[1])
                                + POW2(s->a_uOtherPhase(pCellId, 2) - u[2]));
        nu_BIT = m_Cs * s->a_alphaGasLim(pCellId) * u_d;
      }
    }
    // The stress tensor S_mean depends on tau and tau depends on S_mean.
    // This leads to a quadratic equation, that is solved by this.
    MFloat tauBase = F1B2 + 3.0 * (s->m_nu + nu_BIT) * FFPOW2(s->maxLevel() - s->a_level(pCellId));
    const MFloat tau = F1B2
                       * (tauBase
                          + sqrt(tauBase * tauBase
                                 + 2.0 * m_Cs * m_Cs * m_deltaX * m_deltaX * (F1BCSsq * F1BCSsq) * Q
                                       * FPOW2(s->maxLevel() - s->a_level(pCellId))));

    // eddy viscosity
    if constexpr(bubble) {
      if(s->m_cells.saveOldNu()) {
        s->a_oldNu(pCellId) = s->a_nu(pCellId);
      }
    }
    s->a_nu(pCellId) = (2.0 * tau - 1.0) * CSsq * FPOW2(s->maxLevel() - s->a_level(pCellId)) / 2.0;
    if constexpr(bubble) {
      if(s->m_isEELiquid) {
        if(s->m_cells.saveOldNu()) {
          s->a_oldNuT(pCellId) = s->a_nuT(pCellId);
        }
        s->a_nuT(pCellId) = s->a_nu(pCellId) - s->m_nu - nu_BIT;
      }
    }
  });
}

//-registration for LbSrcTermFactory class
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<2, 9, LbSysEqnIncompressible<2, 9>, false>>
    reg_Smagorinskyd2q9("LB_SRC_SMAGORINSKY");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 19, LbSysEqnIncompressible<3, 19>, false>>
    reg_Smagorinskyd3q19("LB_SRC_SMAGORINSKY");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 27, LbSysEqnIncompressible<3, 27>, false>>
    reg_Smagorinskyd3q27("LB_SRC_SMAGORINSKY");

static const LbRegSrcTerm<LbSrcTerm_smagorinsky<2, 9, LbSysEqnIncompressible<2, 9>, true>>
    reg_SmagorinskydBubble2q9("LB_SRC_SMAGORINSKY_BUBBLE");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 19, LbSysEqnIncompressible<3, 19>, true>>
    reg_SmagorinskydBubble3q19("LB_SRC_SMAGORINSKY_BUBBLE");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 27, LbSysEqnIncompressible<3, 27>, true>>
    reg_SmagorinskydBubble3q27("LB_SRC_SMAGORINSKY_BUBBLE");

static const LbRegSrcTerm<LbSrcTerm_smagorinsky<2, 9, LbSysEqnCompressible<2, 9>, false>>
    reg_Smagorinskyd2q9C("LB_SRC_SMAGORINSKY");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 19, LbSysEqnCompressible<3, 19>, false>>
    reg_Smagorinskyd3q19C("LB_SRC_SMAGORINSKY");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 27, LbSysEqnCompressible<3, 27>, false>>
    reg_Smagorinskyd3q27C("LB_SRC_SMAGORINSKY");

static const LbRegSrcTerm<LbSrcTerm_smagorinsky<2, 9, LbSysEqnCompressible<2, 9>, true>>
    reg_SmagorinskydBubble2q9C("LB_SRC_SMAGORINSKY_BUBBLE");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 19, LbSysEqnCompressible<3, 19>, true>>
    reg_SmagorinskydBubble3q19C("LB_SRC_SMAGORINSKY_BUBBLE");
static const LbRegSrcTerm<LbSrcTerm_smagorinsky<3, 27, LbSysEqnCompressible<3, 27>, true>>
    reg_SmagorinskydBubble3q27C("LB_SRC_SMAGORINSKY_BUBBLE");

} // namespace maia::lb
