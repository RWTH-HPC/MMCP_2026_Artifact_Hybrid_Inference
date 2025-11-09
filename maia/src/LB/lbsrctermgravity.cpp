// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsrctermgravity.h"

#include "UTIL/parallelfor.h"
#include "lbsrctermcontroller.h"
#include "lbsyseqn.h"

namespace maia::lb {

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermGravity<nDim, nDist, SysEqn>::readProperties() {
  auto& s = m_solver; // alias for readability
  const MInt solverId = s->solverId();

  if(Context::propertyExists("volumeAcceleration", solverId)) {
    m_mode = Mode::DIRECT;
    for(MInt i = 0; i < nDim; i++) {
      /*! \page propertyPage1
        \section volumeAcceleration
        default = <code>[0.0, 0.0, 0.0]</code>\n\n
        This property defines the amount of acceleration applied in each Cartesian direction
        The acceleration has to be specified in LB non-dimensional form
        by multiplying by (\delta t ^ 2) / (\delta x) = (\delta x) / (3 * a_inf^2) (see wiki)\n
        Keywords: <i>LATTICE_BOLTZMANN</i>
      */
      m_acceleration[i] = Context::getSolverProperty<MFloat>("volumeAcceleration", solverId, AT_, i);
    }
  }
  /*! \page propertyPage1
    \section Ga
    <code>MFloat LbSolver::m_Ga</code>\n
    default = <code>0.0</code>\n\n
    Galileo number for each Cartesian direction\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  if(Context::propertyExists("Ga", solverId)) {
    m_mode = Mode::GALILEO;
    ASSERT(Context::propertyLength("Ga") == nDim, "Ga should be of size nDim");
    const MFloat nu = s->m_Ma / F1BCS / s->m_Re * s->m_referenceLength;
    for(MInt n = 0; n < nDim; n++) {
      m_Ga[n] = Context::getSolverProperty<MFloat>("Ga", solverId, AT_, m_Ga.data(), n);
      const MFloat gravity = m_Ga[n] * POW2(nu) / POW3(s->m_referenceLength);
      m_acceleration[n] = gravity;
    }
  }

  // Save density gradient for initial condition
  for(MInt n = 0; n < nDim; n++) {
    m_solver->m_volumeAccel[n] = m_acceleration[n];
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermGravity<nDim, nDist, SysEqn>::init() {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  std::array<MFloat, nDim> densityGradient{};
  for(MInt n = 0; n < nDim; n++) {
    densityGradient[n] = m_acceleration[n] * F1BCSsq;
  }

  for(MInt dir = 0; dir < nDim; dir++) {
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_forcing[Ld::nFld(dir, mi)] += Ld::tp(Ld::distType(Ld::nFld(dir, mi))) * -1.0 * densityGradient[dir];
      m_forcing[Ld::pFld(dir, mi)] += Ld::tp(Ld::distType(Ld::pFld(dir, mi))) * 1.0 * densityGradient[dir];
    }
  }

  // Save forcing for BC and interface
  for(MInt j = 0; j < nDist; j++) {
    m_solver->m_Fext[j] = m_forcing[j];
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSrcTermGravity<nDim, nDist, SysEqn>::apply_postPropagation() {
  TRACE();
  const auto& s = m_solver; // alias for readability

  const MInt timeStep = s->getCurrentTimeStep();
  maia::parallelFor<true>(0, s->m_currentMaxNoCells, [=](MInt index) {
    const MInt cellId = s->m_activeCellList[index];
    const MInt lvlDiff = s->maxLevel() - s->a_level(cellId);
    if((timeStep) % IPOW2(lvlDiff) != 0) return;
    for(MInt j = 0; j < nDist - 1; j++) {
      s->a_oldDistribution(cellId, j) += FPOW2(lvlDiff) * m_forcing[j];
    }
  });
}

//---registration for LbSrcTermFactory class------------------------------------
static const LbRegSrcTerm<LbSrcTermGravity<2, 9, LbSysEqnCompressible<2, 9>>> reg_SpongeRhoConstd2q9C("LB_GRAVITY");
static const LbRegSrcTerm<LbSrcTermGravity<2, 9, LbSysEqnIncompressible<2, 9>>> reg_SpongeRhoConstd2q9("LB_GRAVITY");
static const LbRegSrcTerm<LbSrcTermGravity<3, 27, LbSysEqnCompressible<3, 27>>> reg_SpongeRhoConstd3q27C("LB_GRAVITY");
static const LbRegSrcTerm<LbSrcTermGravity<3, 27, LbSysEqnIncompressible<3, 27>>> reg_SpongeRhoConstd3q27("LB_GRAVITY");

} // namespace maia::lb
