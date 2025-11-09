// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbdgape.h"

//--
//--public
//--

template <MInt nDim, MInt nDist, class SysEqn>
void LbDgApe<nDim, nDist, SysEqn>::finalizeCouplerInit() {
  TRACE();
  // TODO labels:COUPLER,LB,DG for now, DG does the same time step as LB on finest level
  this->dgSolver().forceTimeStep(m_conversionLb2Dg.time * 1.0);
}

//--
//--private
//--

template <MInt nDim, MInt nDist, class SysEqn>
void LbDgApe<nDim, nDist, SysEqn>::initConversionFactors() {
  // TODO labels:lb Miro: document why conversion factors are needed for inactive solver
  TRACE();

  // Description: for generic variable phi conversion factor C_phi is defined as
  //              phi_dg = phi_lb * C_phi
  // base units
  constexpr MFloat dgGamma = 1.4;
  MFloat maSq = this->a_Ma() * this->a_Ma(); // Note: Ma is present in inactive LB solver

  // Inactive donor solver does not know the maxLevel
  const MBool hasLbSolver = this->lbSolver().isActive();
  MFloat dxLb =
      (hasLbSolver) ? this->a_cellLengthAtLevel(donorSolver().maxLevel()) : std::numeric_limits<MFloat>::max();
  MPI_Allreduce(MPI_IN_PLACE,
                &dxLb,
                1,
                maia::type_traits<MFloat>::mpiType(),
                MPI_MIN,
                this->dgSolver().grid().raw().mpiComm(),
                AT_,
                "MPI_IN_PLACE",
                "dxLb");

  m_conversionLb2Dg.length = dxLb;
  m_conversionLb2Dg.density = pow(1.0 + 0.5 * (dgGamma - 1.0) * maSq, 1.0 / (1.0 - dgGamma));
  m_conversionLb2Dg.time = dxLb * LBCS * sqrt(1.0 + 0.5 * (dgGamma - 1.0) * maSq); // sqrt-term: a_0 / a_inf
  // derived units
  m_conversionLb2Dg.velocity = m_conversionLb2Dg.length / m_conversionLb2Dg.time;
  m_conversionLb2Dg.vorticity = 1.0 / m_conversionLb2Dg.time;
  m_conversionLb2Dg.lamb = m_conversionLb2Dg.vorticity * m_conversionLb2Dg.velocity;
}

/** \brief  Get velocity and vorticity of donor solver
 *  \author Miro Gondrum
 *  \date   03.08.2021
 *
 *  \param[in]  donorCellIds  collections of cellIds of interest
 *  \param[out] p_velocity    pointer to write velocity field
 *  \param[out] p_vorticity   pointer to write vorticity field
 *
 *  \note   Here, in case of LB solver as donor solver, the units are converted
 *          into the dimensions used by the DG solver.
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbDgApe<nDim, nDist, SysEqn>::getDonorVelocityAndVorticity(const std::vector<MInt>& donorCellIds,
                                                                MFloatScratchSpace& p_velocity,
                                                                MFloatScratchSpace& p_vorticity) {
  TRACE();
  if(this->m_hasDonorCartesianSolver) {
    // Velocity
    for(auto donorId : donorCellIds) {
      const MFloat* cellVars = nullptr;
      this->donorSolver().getSampleVariables(donorId, cellVars);
      for(MInt dirId = 0; dirId < this->noVelocities(); dirId++) {
        p_velocity(donorId, dirId) = cellVars[this->donorSolver().PV->VV[dirId]] * m_conversionLb2Dg.velocity;
      }
    }
    // Vorticity
    for(auto donorId : donorCellIds) {
      MFloat velocityGradient[nDim][nDim];
      this->donorSolver().calculateVelocityDerivative(donorId, velocityGradient);
      if constexpr(nDim == 2) {
        p_vorticity(donorId, 0) = (velocityGradient[1][0] - velocityGradient[0][1]) * m_conversionLb2Dg.vorticity;
      } else if(nDim == 3) {
        p_vorticity(donorId, 0) = (velocityGradient[2][1] - velocityGradient[1][2]) * m_conversionLb2Dg.vorticity;
        p_vorticity(donorId, 1) = (velocityGradient[0][2] - velocityGradient[2][0]) * m_conversionLb2Dg.vorticity;
        p_vorticity(donorId, 2) = (velocityGradient[1][0] - velocityGradient[0][1]) * m_conversionLb2Dg.vorticity;
      }
    }
  }
}

/** \brief  Transform data from donor solver unit system into DG unit system
 *  \author Miro Gondrum
 *  \date   03.08.2021
 *  \note   Here, in case of LB solver as donor solver, the units are converted
 *          into the dimensions used by the DG solver.
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbDgApe<nDim, nDist, SysEqn>::performUnitConversion(const MString& name,
                                                         const MInt count,
                                                         const MInt stride,
                                                         MFloat* data) {
  TRACE();
  // Determine corresponding conversion factor
  std::map<MString, MFloat> conversionMap;
  conversionMap["wxv_x"] = m_conversionLb2Dg.lamb;
  conversionMap["wxv_y"] = m_conversionLb2Dg.lamb;
  conversionMap["wxv_z"] = m_conversionLb2Dg.lamb;
  conversionMap["um"] = m_conversionLb2Dg.velocity;
  conversionMap["vm"] = m_conversionLb2Dg.velocity;
  conversionMap["wm"] = m_conversionLb2Dg.velocity;
  conversionMap["rhom"] = m_conversionLb2Dg.density;
  conversionMap["c0"] = m_conversionLb2Dg.velocity;
  conversionMap["dc0_dx"] = m_conversionLb2Dg.vorticity;
  conversionMap["dc0_dy"] = m_conversionLb2Dg.vorticity;
  conversionMap["dc0_dz"] = m_conversionLb2Dg.vorticity;
  conversionMap["vort_x"] = m_conversionLb2Dg.vorticity;
  conversionMap["vort_y"] = m_conversionLb2Dg.vorticity;
  conversionMap["vort_z"] = m_conversionLb2Dg.vorticity;

  auto search = conversionMap.find(name);
  if(search == conversionMap.end()) {
    std::stringstream ss;
    ss << "ERROR: No unit type found for " << name << " ." << std::endl;
    TERMM(1, ss.str());
  }
  const MFloat conversionFactor = conversionMap[name];
  for(MInt i = 0; i < count; i++) {
    const MInt j = i * stride;
    data[j] *= conversionFactor;
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbDgApe<nDim, nDist, SysEqn>::calcSourceLambNonlinear(const MFloat* const velocity,
                                                           const MFloat* const vorticity,
                                                           MFloat* const sourceTerms) {
  TRACE();
  RECORD_TIMER_START(this->m_timers[BaseDg::Timers::CalcSourceLamb]);
  // Index of first lamb vector component in m_meanVars
  const MInt lambIndex = this->m_meanVarsIndex[BaseDg::MV::LAMB0[0]];
  // Calculate source terms on all 'active' donor leaf cells, i.e. that are
  // mapped to elements with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < this->m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = this->m_calcSourceDonorCells[donorIndex];
    // Create convenience pointers
    const MFloat* const meanVars = &this->m_meanVars[donorIndex * this->noMeanVars()];
    const MFloat* const velo = &velocity[donorId * this->noVelocities()];
    const MFloat* const vort = &vorticity[donorId * this->noVorticities()];
    // Compute instantaneous Lamb vector
    MFloat lambVector[nDim];
    if(nDim == 2) {
      lambVector[0] = -vort[0] * velo[1];
      lambVector[1] = vort[0] * velo[0];
    } else {
      lambVector[0] = vort[1] * velo[2] - vort[2] * velo[1];
      lambVector[1] = vort[2] * velo[0] - vort[0] * velo[2];
      lambVector[2] = vort[0] * velo[1] - vort[1] * velo[0];
    }
    // The complete source term is:
    // S = -L' = -(w x u - bar(w x u))
    const MInt offset = donorId * SysEqnDg::noVars() + BaseDg::CV::UU[0];
    for(MInt dim = 0; dim < nDim; dim++) {
      sourceTerms[offset + dim] += -(lambVector[dim] - meanVars[lambIndex + dim]);
    }
  }
  RECORD_TIMER_STOP(this->m_timers[BaseDg::Timers::CalcSourceLamb]);
}

//--
//--explicit instantiations
//--
template class LbDgApe<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LbDgApe<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
