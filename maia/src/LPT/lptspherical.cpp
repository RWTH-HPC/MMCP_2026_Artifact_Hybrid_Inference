// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptspherical.h"
#include <cmath>
#include <functional>
#include "lpt.h"

using namespace std;
using namespace maia::lpt;

// reference values used for LPT non-dimensionalisation:
template <MInt nDim>
MFloat LPTSpherical<nDim>::s_lengthFactor = 1.0;
// L_ref / L_refGrid (length factor between the particle reference length
//                    and the grid reference length and)
template <MInt nDim>
MFloat LPTSpherical<nDim>::s_Re = 1.0;
// reference Re-number: Re_ref = rho_ref * u_ref * L_ref / mu_ref
template <MInt nDim>
MFloat LPTSpherical<nDim>::s_Pr = 1.0;
// reference Pr-number: Pr_ref = mu_ref * Cp_ref / lamda_ref
template <MInt nDim>
MFloat LPTSpherical<nDim>::s_Sc = 1.0;
// reference Sc-number: Sc_ref = mu_ref / (rho_ref * D_ref)
template <MInt nDim>
MFloat LPTSpherical<nDim>::s_We = 1.0;
// reference We-number: We_ref = rho_ref * u_ref^2 * L_ref / sigma_ref
template <MInt nDim>
std::array<MFloat, nDim> LPTSpherical<nDim>::s_Frm{};
// modified reference Fr-number with : Fr_ref = u_ref / sqrt( g * L_ref)
// the modified version is Frm = 1/ (Fr_ref^2) used to reduce computations

/**
 * \brief constructor of spherical particles, used to invalidate values when debugging!
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
LPTSpherical<nDim>::LPTSpherical() {
  //#ifdef LPT_DEBUG
  const auto nan = std::numeric_limits<MFloat>::quiet_NaN();

  for(MInt i = 0; i < nDim; i++) {
    m_fluidVel[i] = nan;
    m_oldFluidVel[i] = nan;
    m_velocity[i] = nan;
    m_oldVel[i] = nan;
    m_position[i] = nan;
    m_oldPos[i] = nan;
    m_accel[i] = nan;
    m_oldAccel[i] = nan;
  }
  //#endif
}


/**
 * \brief advance to new timeStep
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPTSpherical<nDim>::advanceParticle() {
  firstStep() = false;
  hasCollided() = false;
  hadWallColl() = false;
  // m_creationTime = 0.0;
  for(MInt i = 0; i < nDim; i++) {
    ASSERT(!isnan(m_fluidVel[i]), "");
    m_oldFluidVel[i] = m_fluidVel[i];
  }
  ASSERT(m_fluidDensity > 0, "");
  m_oldFluidDensity = m_fluidDensity;
}


/**
 * \brief single particle coupling terms
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPTSpherical<nDim>::coupling() {
  if(m_cellId < 0) {
    mTerm(1, AT_, "coupling of invalid particle?");
  }

  if(!s_backPtr->m_heatCoupling && !s_backPtr->m_evaporation) {
    // NOTE: set weights for later source term redistribution if not done for the evaporation yet!
    m_redistWeight.clear();
    interpolateVelocityAndDensity(m_cellId, m_position.data(), m_fluidVel.data(), m_fluidDensity, m_redistWeight);
  }

  // check that all variables have been exchanged/set correctly
  for(MInt i = 0; i < nDim; i++) {
    ASSERT(!isnan(m_oldFluidVel[i]), "");
  }
  ASSERT(m_oldFluidDensity > 0, "");

  // check all source terms
  // TODO: add debug-flags!
  if((s_backPtr->m_heatCoupling && std::isnan(m_heatFlux)) || (s_backPtr->m_massCoupling && std::isnan(m_dM))) {
    cerr << "Invalid energy equation!" << m_partId << " " << m_position[0] << " " << m_position[1] << " "
         << m_position[2] << endl;
  }
  for(MInt i = 0; i < nDim; i++) {
    if(std::isnan(m_accel[i])) {
      cerr << "Invalid motion equation!" << m_partId << " " << m_position[0] << " " << m_position[1] << " "
           << m_position[2] << endl;
      isInvalid() = true;
      return;
    }
    if(std::isnan(m_velocity[i])) {
      cerr << "Invalid motion equation2!" << m_partId << " " << m_position[0] << " " << m_position[1] << " "
           << m_position[2] << " " << m_densityRatio << " " << hadWallColl() << endl;
      isInvalid() = true;
      return;
    }
  }

  if(s_backPtr->m_momentumCoupling) momentumCoupling();
  if(s_backPtr->m_heatCoupling) heatCoupling();
  if(s_backPtr->m_massCoupling) massCoupling();
}
/**
 * \brief add the mass flux from the particle to the cell/all surrounding cells
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
void LPTSpherical<nDim>::massCoupling() {
  if(m_creationTime > s_backPtr->m_timeStep) {
    cerr << "Limiting creation time from " << m_creationTime << " to " << s_backPtr->m_timeStep - MFloatEps << endl;
    m_creationTime = s_backPtr->m_timeStep - MFloatEps;
  }

  const MFloat relT = (s_backPtr->m_timeStep - m_creationTime) / s_backPtr->m_timeStep;

  if(m_dM < 0 || relT < 0 || relT > 1) {
    cerr << m_creationTime << " " << m_partId << " " << m_dM << " " << relT << " " << s_backPtr->m_timeStep << endl;
    mTerm(1, AT_ "Invalid mass flux");
  }

  // NOTE: negative massFlux means going from liquid/disperse(LPT) to fluid/continous(FV)
  const MFloat massFlux = -m_noParticles * m_dM * relT;

  s_backPtr->m_sumEvapMass -= (massFlux * s_backPtr->m_timeStep);

  if(!s_backPtr->m_couplingRedist) {
    s_backPtr->a_massFlux(m_cellId) += massFlux;
  } else {
    for(MInt n = 0; n < (signed)m_neighborList.size(); n++) {
      s_backPtr->a_massFlux(m_neighborList[n]) += m_redistWeight.at(n) * massFlux;
    }
  }
}


/**
 * \brief add the momentum flux from the particle to the cell/all surrounding cells
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
void LPTSpherical<nDim>::momentumCoupling() {
  const MFloat relT = (s_backPtr->m_timeStep - m_creationTime) / s_backPtr->m_timeStep;

  // NOTE: as the momentum coupling is applied after the evaporation, the original mass needs
  //      to be used for the force coupling!
  const MFloat mass =
      (s_backPtr->m_evaporation) ? m_dM * relT * s_backPtr->m_timeStep + sphericalMass() : sphericalMass();

  array<MFloat, nDim> force{};
  array<MFloat, 3> evapMom{F0, F0, F0};

  for(MInt i = 0; i < nDim; i++) {
    // momentum flux due to drag force
    const MFloat invDensityRatio = (abs(m_densityRatio) <= MFloatEps) ? 1.0 : 1.0 / m_densityRatio;
    force[i] = m_noParticles * mass * (m_accel[i] - ((1.0 - invDensityRatio) * s_Frm[i])) * relT;

    // momentum flux due to mass change going into the continous/fluid(FV)
    // its energy is considered in the heat-coupling
    if(s_backPtr->m_evaporation) evapMom[i] = -m_noParticles * m_dM * m_velocity[i] * relT;
  }

  // work done by the momentum
  // const MFloat work = std::inner_product(force.begin(), force.end(), m_velocity.begin(), 0.0);
  const MFloat work = std::inner_product(force.begin(), force.end(), m_oldVel.begin(), 0.0);


  if(!s_backPtr->m_couplingRedist) {
    for(MInt i = 0; i < nDim; i++) {
      s_backPtr->a_momentumFlux(m_cellId, i) += (force[i] + evapMom[i]);
    }
    s_backPtr->a_workFlux(m_cellId) += work;

  } else {
    ASSERT(m_redistWeight.size() == m_neighborList.size(),
           "Invalid " + to_string(m_redistWeight.size()) + "!=" + to_string(m_neighborList.size()));
    for(MInt n = 0; n < (signed)m_neighborList.size(); n++) {
      // a weightHeat is used to distribute the source term to multiple cells...
      for(MInt i = 0; i < nDim; i++) {
        s_backPtr->a_momentumFlux(m_neighborList[n], i) += m_redistWeight.at(n) * (force[i] + evapMom[i]);
      }
      s_backPtr->a_workFlux(m_neighborList[n]) += m_redistWeight.at(n) * work;
    }
  }
}


/**
 * \brief add the heat flux from the particle to the cell/all surrounding cells
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
void LPTSpherical<nDim>::heatCoupling() {
  const MFloat relT = (s_backPtr->m_timeStep - m_creationTime) / s_backPtr->m_timeStep;

  // kin. energy change due to evaporation
  MFloat velMagSquared = 0;
  for(MInt i = 0; i < nDim; i++) {
    velMagSquared += POW2(m_oldVel[i]);
  }

  // NOTE: must be negative if going trom the disperse/liquid(LPT) to the continous/fluid (FV)
  const MFloat energySource = m_noParticles * (m_heatFlux - 0.5 * m_dM * velMagSquared);

  if(s_backPtr->m_couplingRedist) {
    for(MInt n = 0; n < (signed)m_neighborList.size(); n++) {
      // a weight is used to distribute the source term to multiple cells
      s_backPtr->a_heatFlux(m_neighborList[n]) += m_redistWeight.at(n) * energySource * relT;
    }
  } else {
    s_backPtr->a_heatFlux(m_cellId) += energySource * relT;
  }
}


template <MInt nDim>
void LPTSpherical<nDim>::motionEquation() {
  const MFloat dt = s_backPtr->m_timeStep - m_creationTime;
  ASSERT(dt > 0, "Timestep < 0");

  // NOTE: at this point old means 1 time-Step ago and the current values will be updated below!
  m_oldPos = m_position;
  m_oldVel = m_velocity;
  m_oldAccel = m_accel;
  m_oldCellId = m_cellId;

#ifdef LPT_DEBUG
  MLong debugPartId = -1;
  if(m_partId == debugPartId) {
    cerr << "BT " << m_velocity[0] << " " << m_velocity[1] << " " << m_velocity[2] << endl;
  }
#endif

  if(firstStep()) {
    for(MInt i = 0; i < nDim; i++) {
      m_accel[i] = s_Frm[i];
      m_oldAccel[i] = s_Frm[i];
    }
  } /*else {
    m_creationTime = 0;
  } */

  for(MInt i = 0; i < nDim; i++) {
    ASSERT(!isnan(m_oldFluidVel[i]), "");
  }
  ASSERT(m_oldFluidDensity > 0, "");

  // reduce variance between runs (round to 9 decimal places)
  // TODO-timw labels:LPT,totest check if rounding is still necessary?
  const MFloat invDensityRatio =
      1.0 / (static_cast<MFloat>(ceil(s_backPtr->m_material->density(m_temperature) / m_oldFluidDensity * 1E9)) / 1E9);

  m_densityRatio = 1 / invDensityRatio;

  /// 1. Predictor step
  array<MFloat, nDim> predictedPos{};
  array<MFloat, nDim> predictedVel{};

  // predict position and veloctity
  for(MInt i = 0; i < nDim; i++) {
    predictedVel[i] = m_oldVel[i] + dt * m_oldAccel[i];
    predictedPos[i] = m_oldPos[i] + dt * m_oldVel[i] * s_lengthFactor + 0.5 * dt * dt * m_oldAccel[i] * s_lengthFactor;
  }

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    cerr << "PT " << m_velocity[0] << " " << m_velocity[1] << " " << m_velocity[2] << endl;
  }
#endif

  /// 2. get velocity, density and velocity magnitude at the predicted position
  //  find the predicted CellId at the predicted position around the previous cellId
  const MInt predictedCellId = s_backPtr->grid().findContainingLeafCell(&predictedPos[0], m_cellId);

  if(predictedCellId < 0 || m_cellId < 0) {
    // particle has left the LPT bounding-box completely
    //->to-be removed!
    isInvalid() = true;
    return;
  }

  ASSERT(m_cellId > -1, "Invalid cell!");
  ASSERT(predictedCellId > -1, "Invalid predicetd cell!");
  ASSERT(s_backPtr->c_isLeafCell(predictedCellId), "Invalid cell!");
  ASSERT(predictedCellId > -1, "Invalid cell!");

  array<MFloat, nDim> fluidVel{};
  MFloat fluidDensity{};
  vector<MFloat> predictedWeights(0);
  interpolateVelocityAndDensity(predictedCellId, predictedPos.data(), fluidVel.data(), fluidDensity, predictedWeights);

  const MFloat T = s_backPtr->a_fluidTemperature(m_cellId);

  MFloat fluidViscosity = s_backPtr->m_material->dynViscosityFun(T);

  /// 2. corrector step
  if(s_backPtr->m_dragModelType == 0) {
    for(MInt i = 0; i < nDim; i++) {
      m_velocity[i] = predictedVel[i];
      m_position[i] = predictedPos[i];
      m_accel[i] = m_oldAccel[i];
    }
  } else {
    switch(s_backPtr->m_motionEquationType) {
      case 3:
      case 0: {
        // see also Crowe et al "Multiphase Flows with Droplets and Particles" 1998
        // motion equation based on analytical solution of the simplified motion equation
        // contained when assuming constant old-drag count!
        // in this case the newton equation takes the form of:
        // du/dt = a + b * u
        // for which an exponential analytical solution can be found and
        // defined by using the old-position at relative time t_old = 0.
        // for further details contact t.wegmann@aia.rwth-aachen
        const MFloat relVel_old = magRelVel(&m_oldFluidVel[0], &m_oldVel[0]);
        const MFloat relVel = magRelVel(&fluidVel[0], &predictedVel[0]);

        //#ifdef LPT_DEBUG
        if(isnan(relVel)) {
          cerr << "relVel " << fluidVel[0] << " " << fluidVel[1] << " " << fluidVel[2] << " " << predictedVel[0] << " "
               << predictedVel[1] << " " << predictedVel[2] << " " << predictedCellId << endl;
          mTerm(1, AT_, "Invalid relative velocity!");
        }
        //#endif

        // get the particle drag based on velocities and density of the old position!
        const MFloat DC_old = dragFactor(particleRe(relVel_old, m_oldFluidDensity, fluidViscosity) * s_Re) / s_Re
                              * fParticleRelTime(fluidViscosity);
        // get the particle drag based on velocities and density of the predicted position
        const MFloat DC = dragFactor(particleRe(relVel, fluidDensity, fluidViscosity) * s_Re) / s_Re
                          * fParticleRelTime(fluidViscosity);

        // average the fluid velocity of old and predicted position
        // increases stability for flow fields with high velocity gradients!
        array<MFloat, nDim> avgFlowVel{};
        for(MInt i = 0; i < nDim; i++) {
          avgFlowVel[i] = 0.5 * (m_oldFluidVel[i] + fluidVel[i]);
        }
        const MFloat avgDC = s_backPtr->m_motionEquationType == 0 ? DC_old : 0.5 * (DC_old + DC);

        // analytical solution assuming constant DC
        for(MInt i = 0; i < nDim; i++) {
          m_velocity[i] =
              avgFlowVel[i] + (1.0 - invDensityRatio) * s_Frm[i] / avgDC
              + exp(-avgDC * dt) * (m_oldVel[i] - avgFlowVel[i] + (invDensityRatio - 1.0) * s_Frm[i] / avgDC);

          m_accel[i] = (m_velocity[i] - m_oldVel[i]) / dt;

          m_position[i] = m_oldPos[i] + 0.5 * (m_oldVel[i] + m_velocity[i]) * dt * s_lengthFactor;
        }
        break;
      }
      case 1: {
        // non-dimensional version based on:
        // Baumgarten, "Mixture Formation in Internal Combustion Engines", 2005
        const MFloat relVel = magRelVel(&fluidVel[0], &predictedVel[0]);

        const MFloat Rep = particleRe(relVel, m_oldFluidDensity, fluidViscosity) * s_Re;

        const MFloat CD = dragFactor(Rep) / s_Re * fParticleRelTime(fluidViscosity);

        for(MInt i = 0; i < nDim; i++) {
          m_accel[i] = CD * (fluidVel[i] - predictedVel[i]) + (1.0 - invDensityRatio) * s_Frm[i];
          m_velocity[i] = m_oldVel[i] + dt * F1B2 * (m_accel[i] + m_oldAccel[i]);
          m_position[i] = m_oldPos[i] + dt * F1B2 * (m_velocity[i] + m_oldVel[i]) * s_lengthFactor;
        }
        break;
      }
      case 2: {
        // particle moves with fluid velocity
        for(MInt i = 0; i < nDim; i++) {
          m_position[i] = m_oldPos[i] + 0.5 * dt * (fluidVel[i] + m_oldVel[i]) * s_lengthFactor;
          m_velocity[i] = fluidVel[i];
          m_accel[i] = (m_velocity[i] - m_oldVel[i]) / dt;
        }
        break;
      }
      default: {
        mTerm(1, AT_, "Unknown particle corrector equation type!");
      }
    }
  }

#ifdef LPT_DEBUG
  if(m_partId == debugPartId) {
    cerr << "AT " << m_velocity[0] << " " << m_velocity[1] << " " << m_velocity[2] << endl;
  }
#endif

  // 3. update cellId based on the corrected position and set the new status!
  checkCellChange(&m_oldPos[0]);

#ifdef LPT_DEBUG
  for(MInt i = 0; i < nDim; i++) {
    if(std::isnan(m_accel[i])) {
      cerr << "Nan acc. for " << s_backPtr->domainId() << " " << m_partId << " " << m_position[0] << " "
           << m_position[1] << " " << m_position[nDim - 1] << " " << s_backPtr->a_isValidCell(m_cellId) << " "
           << m_oldPos[0] << " " << m_oldPos[1] << " " << m_oldPos[nDim - 1] << " " << m_creationTime << " "
           << fluidVel[0] << " " << fluidVel[1] << " " << fluidVel[nDim - 1] << fluidDensity << " " << T << " "
           << fluidViscosity << " " << m_oldVel[0] << " " << m_oldVel[1] << " " << m_oldVel[nDim - 1] << endl;
    }
  }
#endif
}


/// \brief Calculate drag factor of the current particle for the given particle Reynolds number.
///        in the literature this is referred to as: C_d * Re_p / 24
/// \author Sven Berger
/// \date   August 2015
/// \param[in] partRe Particle Reynolds number.
template <MInt nDim>
MFloat LPTSpherical<nDim>::dragFactor(const MFloat partRe) {
  MFloat result = NAN;
  switch(s_backPtr->m_dragModelType) {
    case 0: { // drag force switched off
      result = 0.0;
      break;
    }
    case 1: { // linear Stokes drag
      result = 1.0;
      break;
    }
    case 2: { // nonlinear drag according to Schiller & Naumann
      result = 1.0 + 0.15 * pow(partRe, 0.687);
      break;
    }
    case 3: { // nonlinear drag according to Pinsky & Khain, J. Aerosol Sci. 28(7), 1177-1214 (1997).
      result = 1.0 + 0.17 * pow(partRe, 0.632) + 1.0e-6 * pow(partRe, 2.25);
      break;
    }
    case 4: {
      // drag according to Carsten Baumgarten
      // Mixture Formation in Internal Combustion Engines, 2005
      // originates from Putnam 1961
      // (C_d * Re_p / 24 )
      if(partRe > 1000) {
        // Newton flow regime Re from 1k to 250k
        result = 0.424 * partRe / 24.0;
      } else if(partRe <= 0.1) {
        // Stokes flow
        result = 1.0;
      } else {
        // transition regime
        result = 1.0 + pow(partRe, 2.0 / 3.0) / 6.0;
      }
      break;
    }
    case 5: {
      // drag according to Rowe for interphase momentum exchange, 1961
      if(partRe >= 1000) {
        result = 0.44 * partRe / 24.0;
      } else {
        result = 1.0 + 0.15 * pow(partRe, 0.687);
      }
      // limited drag to 10% fluid fraction
      const MFloat fluidFractionUnlimited = 1 - s_backPtr->a_volumeFraction(m_cellId);
      const MFloat fluidFraction = fluidFractionUnlimited < 0.1 ? 0.1 : fluidFractionUnlimited;
      result *= pow(fluidFraction, -2.65) * (fluidFraction - POW2(fluidFraction));
      break;
    }
    case 6: {
      // drag corrected for fluid fraction
      // limited drag to 10% fluid fraction
      const MFloat fluidFractionUnlimited = 1 - s_backPtr->a_volumeFraction(m_cellId);
      const MFloat fluidFraction = fluidFractionUnlimited < 0.1 ? 0.1 : fluidFractionUnlimited;
      result = pow(fluidFraction, -2.65) + 1.0 / 6.0 * pow(partRe, 2.0 / 3.0) * pow(fluidFraction, -1.78);
      break;
    }
    case 7: {
      /* for testing const c_d of 1.83, according to Maier in Multisclae Simulation with a Two-Way
      Coupled Lattice Boltzmann Methoc and Discrete Element Method (DOI: 10.1002/ceat.201600547) */
      result = std::max(1.83 * partRe / 24, 0.1);
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown particle drag method");
    }
  }
  return result;
}


template <MInt nDim>
void LPTSpherical<nDim>::energyEquation() {
  // on first time step for this particle:
  // determine special timestep size which is used to get continuous particle generation
  // i.e. the particle has not lived for the entire timeStep!
  const MFloat dt = s_backPtr->m_timeStep - m_creationTime;
  ASSERT(dt > 0, "Timestep < 0");

  //#ifdef LPT_DEBUG
  if(isnan(m_dM) || isnan(m_diameter) || isnan(m_temperature)) {
    cerr << globalTimeStep << " " << m_dM << " " << m_diameter << " " << m_temperature << endl;
    mTerm(1, AT_, "Nan input values for evaporation!");
  }
  //#endif

  // sub-scripts:
  // fluid : (continous) fluid phase
  // liquid: (disperse, particle) liquid phase
  // vap   : vapour phase of the originally liquid particle
  // mix   : mixture of vapour and fluid phase at the bndry-state

  // surf  : surface state at the droplet interface
  // bndry : (reference) boundray layer state (including fluid and vapour)


  // NOTE: interpolate all fluid variables for evaporation to the particle position and
  //      already set weights for later source term redistribution!
  // As in "Direct numerical simulation of a confined three-dimensional gas mixing layer
  //           with one evaporating hydrocarbon-droplet-laden stream", R.S. Miller and J. Bellan
  m_redistWeight.clear();
  MFloat p_fluid = -1;
  MFloat T_fluid = -1;
  MFloat Y_fluid = 0;
  if(!s_backPtr->m_massCoupling) {
    interpolateAllFluidVariables(m_cellId, m_position.data(), m_fluidVel.data(), m_fluidDensity, p_fluid, T_fluid,
                                 m_redistWeight);
  } else {
    interpolateAllFluidVariables(m_cellId, m_position.data(), m_fluidVel.data(), m_fluidDensity, p_fluid, T_fluid,
                                 Y_fluid, m_redistWeight);
  }

  // particles with small diameters are evaporating fully
  if(m_diameter <= s_backPtr->m_sizeLimit) {
    m_dM = sphericalMass() / dt;
    m_diameter = 0.0;
    isInvalid() = true;
    fullyEvaporated() = true;
    m_heatFlux = 0.0;
    return;
  }

  if(hadWallColl()) {
    // special post-wall collision treatment:
    m_dM = 0;
    m_heatFlux = 0.0;
    return;
  }

  const MFloat volume = sphericalVolume();
  ASSERT(volume > 0, "Invalid volume!");

  static constexpr MFloat convEvap = 1E-10;
  const MFloat oldDiameter = m_diameter;

  MFloat previousTemp = m_temperature;
  MInt heatStep = 0;

  const MFloat oldMass = sphericalMass();
  const MFloat oldTemperature = m_temperature;
  const MFloat oldSheddedMass = 1.0 / 6.0 * PI * s_backPtr->m_material->density(m_temperature) * POW3(m_shedDiam);

  // heat capacity of the liquid phase
  const MFloat cp_liquid = s_backPtr->m_material->cp(m_temperature);

  static constexpr MBool belanHarstad = true;
  if(belanHarstad) {
    // 0) get material properties:

    // boiling point of the disperse phase
    const MFloat T_Boil = s_backPtr->m_material->boilingPoint();
    // referenece pressure for the boiling point of the disperse phase
    const MFloat p_Boil = s_backPtr->m_material->bpRefPressure();

    // specific gas constant of the vapour of the originally disperse phase
    // NOTE: in the non-dimensional formulation: 1/gamma * M_ref / M_p
    const MFloat gasConstant_vap = s_backPtr->m_material->gasConstant();

    // molar weight ration (disperse/continuous) M_p / M_ref
    const MFloat molWeightRatio = s_backPtr->m_material->molWeightRatio();

    // NOTE: for the dimensional case this is just 1
    const MFloat gammaMinusOne = s_backPtr->m_material->gammaMinusOne();

    do {
      previousTemp = m_temperature;

      // particles with small diameters are evaporating fully
      if(m_diameter <= s_backPtr->m_sizeLimit) {
        m_dM = oldMass / dt;
        m_diameter = 0.0;
        m_temperature = oldTemperature;
        isInvalid() = true;
        fullyEvaporated() = true;
        m_heatFlux = 0.0;
        break;
      }

      // 1) use the 1/3-rule to interpolate boundary layer temperature:
      const MFloat T_bndry = 2.0 / 3.0 * m_temperature + 1.0 / 3.0 * T_fluid;

      // 2) use this reference temperture to calculate all material properties:

      // viscosity of the continous phase
      MFloat viscosity_fluid = s_backPtr->m_material->dynViscosityFun(T_bndry);
      if(viscosity_fluid < 0) {
        viscosity_fluid = numeric_limits<MFloat>::epsilon();
      }

      // diffusion coefficient
      MFloat diffusionC = s_backPtr->m_material->diffusionCoefficient(T_bndry);
      ASSERT(diffusionC > 0, "ERROR: Invalid diffusion coefficient.");

      // latent heat of Evaporation of the liquid/disperse phase
      // MFloat lH_ev = s_backPtr->m_material->latentHeatEvap(T_bndry);
      MFloat lH_ev = s_backPtr->m_material->latentHeatEvap(m_temperature);

      // thermal conductivity of the fluid/continous/gas phase
      MFloat lamda_fluid = s_backPtr->m_material->airThermalConductivity(T_bndry);
      if(lamda_fluid < 0) {
        lamda_fluid = numeric_limits<MFloat>::epsilon();
      }

      // thermal condictivity of the liquid/disperse phase
      const MFloat lamda_vap = s_backPtr->m_material->thermalConductivity(T_bndry);

      // dynamic viscosity of the vapour phase
      const MFloat viscosity_vap = s_backPtr->m_material->dynamicViscosity(T_bndry);

      // density of the vapour from the originally liquid/disperse phase
      const MFloat density_vap = p_fluid / (gasConstant_vap * T_bndry);

      // Schmidt-number Sc: viscous diffusion rate / molecular diffusion rate
      // NOTE: in dimensional form for Ranz-Marshall correlation!
      const MFloat Sc = viscosity_fluid / (m_fluidDensity * diffusionC) * s_Sc;

      // Prandtl-number of the fluid/continous/gas phase
      const MFloat Pr = s_backPtr->m_material->airPrandtl(T_bndry) * s_Pr;

      // 3) surface equilibrium mole fraction

      // surface mole fraction is obtained from equilibrium assumption
      // (partial pressure of disperse-material vapour is equal to the equilibrium vapour temperature
      // at the drop temperature) which is used to solve the Clausius-Clapeyron
      // equation for X_F,s (X_surf)
      MFloat X_surf = p_Boil / p_fluid * exp(lH_ev / gasConstant_vap * (1.0 / T_Boil - 1.0 / m_temperature));

      static constexpr MFloat X_surfLimit = 0.99999;
      if(X_surf > X_surfLimit) {
        X_surf = X_surfLimit;
      } else if(X_surf < 0) {
        X_surf = 1.0 - X_surfLimit;
      }

      // 4) non-equilibrium corrections to the surface mole-fraction
      static constexpr MBool nonEquilibrium = true;
      MFloat beta = NAN;
      if(nonEquilibrium) {
        /*
        // 4.1) the beta-factor
        // NOTE: in dimensional system
        if(heatStep == 0) {
          beta = (cp_fluid * oldDM) / (2.0 * PI * lamda_fluid * m_diameter) * s_Pr * s_Re;
        } else {
          beta = (cp_fluid * m_dM) / (2.0 * PI * lamda_fluid * m_diameter) * s_Pr * s_Re;
        }
        */

        // blowing Re-number
        const MFloat Re_b = m_dM / (PI * viscosity_fluid * m_diameter) * s_Re;
        beta = 0.5 * Pr * Re_b;

        // 4.2) the Langmuir-Knudsen law
        // dimensionless Knudsen-layer parameter ( Lk / diameter)
        const MFloat Lk_d =
            viscosity_fluid * sqrt(2.0 * PI * m_temperature * gasConstant_vap) / (Sc * p_fluid * m_diameter * s_Re);

        // non-equilibrium correction:
        X_surf = X_surf - (2.0 * beta * Lk_d);

        // limiting...
        if(X_surf < 0) {
          X_surf = 1.0 - X_surfLimit;
        }
      }

      ASSERT(X_surf <= 1.0, "ERROR: Invalid corrected fuel surface mole concentration (>100%): " + to_string(X_surf));
      ASSERT(X_surf >= 0.0, "ERROR: Invalid corrected fuel surface mole concentration (<0%): " + to_string(X_surf));

      // 5) surface fuel mass fraction
      MFloat Y_surf = X_surf * molWeightRatio / (X_surf * molWeightRatio + (1.0 - X_surf));

      ASSERT(Y_surf <= 1.0, "ERROR: Invalid fuel surface mass concentration (>100%): " + to_string(Y_surf));
      if(std::isnan(Y_surf)) {
        cerr << "X_surf " << X_surf << " m_temperature " << m_temperature << " " << p_fluid << m_partId << " "
             << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;
      }

      // 6) Spalding mass transfer number BM
      // BM = (Y_F,s - Y_F,infty)/(1-Y_F,s)
      // where Y is the fuel mass fraction at the surface (s) and
      // Y_F,infty is the fuel mass fraction at the the far-field conditions
      // for the droplet particle
      MFloat BM = (Y_surf - Y_fluid) / (1.0 - Y_surf);
      if(BM < MFloatEps) {
        // cerr << "BM " << BM << " Y_surf " << Y_surf << " Y_inf " << Y_fluid << " "
        //     << oldDiameter << " " << beta << " " << X_surf << endl;
        BM = 0;
      }

      // 7) Compute Mixture of disperse material-vapour and continous phase in the boundary layer
      //   based on Yuen and Chen, 1976 i.e. the 1/3-rule / the Wilke Rule (Edwards et al., 1979)

      // boundary layer disperse-material vapour mass fraction values chosen as proposed by
      MFloat Y_bndry = 2.0 / 3.0 * Y_surf + 1.0 / 3.0 * Y_fluid;

      ASSERT(Y_bndry <= 1.0, "ERROR: Invalid fuel boundary layer mass concentration (>100%): " + to_string(Y_bndry));

      if(std::isnan(Y_bndry)) {
        cerr << "Y_surf " << Y_surf << " Y_f " << Y_fluid << endl;
      }

      // boundary layer disperse vapour mole fraction
      MFloat X_bndry = -Y_bndry / (Y_bndry * molWeightRatio - Y_bndry - molWeightRatio);


      // thermal conductivity of vapour+fluid mixture
      MFloat lamda_vapFluid = POW2(1.0 + sqrt(lamda_vap / lamda_fluid) * pow(1 / molWeightRatio, 0.25))
                              / sqrt(8.0 * (1.0 + molWeightRatio));
      if(isnan(lamda_vapFluid)) {
        cerr << "lamda_vap " << lamda_vap << " lamda_fluid " << lamda_fluid << " lamda_vapFluid " << lamda_vapFluid
             << endl;
      }

      MFloat lamda_fluidVap = POW2(1.0 + sqrt(lamda_fluid / lamda_vap) * pow(molWeightRatio, 0.25))
                              / sqrt(8.0 * (1.0 + 1 / molWeightRatio));

      MFloat lamda_mix = X_bndry * lamda_vap / (X_bndry + (1.0 - X_bndry) * lamda_vapFluid)
                         + (1.0 - X_bndry) * lamda_fluid / (X_bndry * lamda_fluidVap + (1.0 - X_bndry));
      if(isnan(lamda_mix)) {
        cerr << "X_bndry " << X_bndry << " lamda_vap " << lamda_vap << " lamda_vapFluid " << lamda_vapFluid << endl;
      }
      ASSERT(!std::isnan(lamda_mix), "ERROR: lamda_mix is NaN!");

      // mixture of dynamic viscosity
      MFloat viscosity_vapFluid = POW2(1.0 + sqrt(viscosity_vap / viscosity_fluid) * pow(1 / molWeightRatio, 0.25))
                                  / sqrt(8.0 * (1.0 + molWeightRatio));

      MFloat viscosity_fluidVap = POW2(1.0 + sqrt(viscosity_fluid / viscosity_vap) * pow(molWeightRatio, 0.25))
                                  / sqrt(8.0 * (1.0 + 1 / molWeightRatio));

      MFloat viscosity_mix = X_bndry * viscosity_vap / (X_bndry + (1.0 - X_bndry) * viscosity_vapFluid)
                             + (1.0 - X_bndry) * viscosity_fluid / (X_bndry * viscosity_fluidVap + (1.0 - X_bndry));

      ASSERT(!std::isnan(viscosity_mix),
             "ERROR: viscosity_mix is NaN!" + to_string(X_bndry) + " " + to_string(viscosity_vapFluid));


      // const MFloat density_mix = 1.0 / (Y_bndry / density_vap + (1.0 - Y_bndry) / m_fluidDensity);
      // ASSERT(density_mix > 0, "ERROR: Invalid mixture density (<0).");
      // debug!
      // const MFloat density_mix = s_backPtr->m_material->density(m_temperature);
      // const MFloat density_mix = m_fluidDensity;
      const MFloat density_mix = Y_bndry * density_vap + (1.0 - Y_bndry) * m_fluidDensity;

      // heat capacity of the mixture
      // const MFloat cp_mix = Y_bndry * cp_liquid + (1-Y_bndry) * cp_fluid;

      // 8) Compute Nu and Sh-number based on emperical correlations

      // particle Reynolds-Number with slip-velocity
      // NOTE: in dimensional formulation, following:
      //"A comparative study of numerical models for Eulerian-Lagrangian simulations
      // of turbulent evaporating sprays" D.I. Kolaitis, M.A. Founti
      // Int. J. of Heat and Fluid Flow 27 (2006) 424-435
      // is using the surface/film viscosity over the gas/fluid viscosity
      MFloat Re = particleRe(m_fluidDensity, magRelVel(m_fluidVel.data(), &m_velocity[0]), viscosity_mix) * s_Re;

      // Sherwood number (Sh) is calculated using the Ranz-Marshall correlation (1952)
      // Sh = Convective Mass transfer / Diffusion rate
      // Limitations of the Ranz-Marshall correlation Re < 200, Pr < 250, Sc < 250
      // NOTE: using dimensional Re/Sc-numbers for the correlation!
      const MFloat Sh = 2.0 + 0.552 * sqrt(Re) * pow(Sc, 1.0 / 3.0);

      //#ifdef LPT_DEBUG
      ASSERT(!isnan(Sh), "Sh is nan!" + to_string(Re) + " " + to_string(Sc));
      //#endif

      // Nusselt number Nu :  Convective heat transfer / Conductive heat transfer
      // NOTE: using dimensional Re/Pr-numbers for the correlation!
      MFloat Nu = 2.0 + 0.552 * sqrt(Re) * pow(Pr, 1.0 / 3.0);

      //#ifdef LPT_DEBUG
      if(isnan(Nu)) cerr << " Re " << Re << " Pr " << Pr << endl;
      //#endif

      // non-equilibrium correction of Nusselt-number
      if(nonEquilibrium) {
        if(beta > MFloatEps) {
          Nu = Nu * (beta / (exp(beta) - 1));
        }
#ifdef LPT_DEBUG
        if(isnan(Nu)) cerr << "beta " << beta << endl;
#endif
      }

      MFloat mass = oldMass;

      // 9) solve evaporation equation (dm/dt)
      if(s_backPtr->m_evaporation) {
        MFloat A = 0;
        ///////////////////////////////////////////////////////////////////////////////////////////////
        /// This is the analytical solution of the following differential equation:
        /// Spalding equation (1953) gives the mass evaporation rate
        ///  dm/dt = - pi * d * density * D_AB * Sh * ln(1 + BM) as also in
        /// "A compararive study of numerical models for Eulerian-Lagrangian simulations of
        ///  turbulent evaporating sprays" D.I. Kolaitis, M.A. Founti J. Heat and Fluid Flow (2006)
        ///  This is the same equation as used by Miller/Bellan in:
        ///  "Direct numerical simulation of a confined three-dimensional gar mixing layer with one
        ///  evaporating hydrocarbon-droplet-laden stream" J. Fluid Mech. (1999) but they have
        ///  rearranged it in the form of:
        ///  dm/dt = - Sh / (3 * Sc) * m / tau * ln(1+BM)
        ///  using particle-relaxation time tau.
        ///  For detais of the analytical solution contact t.wegmann@aia.rwth-aachen.de
        ///////////////////////////////////////////////////////////////////////////////////////////////
        if(BM < MFloatEps) {
          mass = oldMass;
        } else {
          A = -2.0 / 3.0 * density_mix * PI * diffusionC * Sh * log(1 + BM)
              * pow(6.0 / (PI * s_backPtr->m_material->density(m_temperature)), 1.0 / 3.0) * dt / (s_Sc * s_Re);

          const MFloat B = pow(oldMass, 2.0 / 3.0);

          //#ifdef LPT_DEBUG
          if(isnan(B)) {
            cerr << "oldMass " << oldMass << endl;
          }
          if(isnan(A)) {
            cerr << "A " << A << " density_mix " << density_mix << " diffusion Coefficient " << diffusionC << " SH "
                 << Sh << " BM " << BM << endl;
          }
          //#endif

          // limit mass to be > 0
          if(fabs(A) > B) {
            mass = 0.0;
          } else {
            mass = pow(B + A, 3.0 / 2.0);
          }
        }

        ASSERT(!std::isnan(mass), "ERROR: Mass is NaN!");

        // NOTE: the direction from m_dM is switched, instead of beeing negtive for evaporation
        //       the sign is changed here!
        // TODO-timw labels:LPT switch m_dM around to match the description in the literature!
        //           use below and switch signs everywhere!
        // m_dM = (mass - oldMass) / dt;
        m_dM = (oldMass - mass) / dt;

        constexpr static MBool limitCondensation = true;
        if(limitCondensation && m_dM < 0) {
          m_dM = 0;
          mass = oldMass;
        }

        ASSERT(!std::isnan(m_dM), "ERROR: Evaporation rate is NaN!");
        ASSERT(mass - oldMass <= 0 || abs(mass) < numeric_limits<MFloat>::epsilon(),
               "ERROR: droplets are increasing in size due to condensation! m_dM " + to_string(m_dM) + " oldMass "
                   + to_string(oldMass) + " mass " + to_string(mass));


        if(mass > 0) {
          m_diameter = pow(mass / s_backPtr->m_material->density(m_temperature) * 6.0 / PI, 1.0 / 3.0);
          // reduce shedded diameter consistently based on evaporated mass m_dM * dt
          // relevant for KH-secondary break-up
          const MFloat reducedSheddedMass = oldSheddedMass - m_dM * dt;
          m_shedDiam = pow(reducedSheddedMass / s_backPtr->m_material->density(m_temperature) * 6.0 / PI, 1.0 / 3.0);

        } else { // catch case in which particle evaporates completely
          m_dM = oldMass / dt;
          m_diameter = 0.0;
          mass = 0.0;
          isInvalid() = true;
          fullyEvaporated() = true;
          m_temperature = oldTemperature;
          m_heatFlux = 0;
          // cerr << "Evap. Completly " << mass << " " << oldMass << " " << BM
          //     << " " << oldTemperature << " " << m_temperature << " " << density_mix << " " <<
          //     A << " " << Sh << " " << Re << " " << Sc << " " << dt << endl;
          // cerr << m_temperature << " " << oldTemperature << " " << T_bndry
          //     << " " << T_fluid
          //     << endl;
          // cerr << "Fully evaporating particle with mass " << oldMass * m_noParticles << endl;
          break;
        }
      }

      // 10) solve temperature equation (dT/dt)

      /////////////////////////////////////////////////////////////////////
      /// The heat balance equation for a droplet was given by Faeth 1983 to be
      /// (tp1-tp0)/dt = (pi * dp1 * Nu * km * (T_f - tp1) + dm1 * LH_ev)/(mp1 * cp)
      /// again this is the same equation as used by Miller/Bellan in:
      /// "Direct numerical simulation of a confined three-dimensional gar mixing layer with one
      ///  evaporating hydrocarbon-droplet-laden stream" J. Fluid Mech. (1999) but they have
      ///  rearranged it in the form of:
      ///  dT/dt = Nu / (3 * Pr) * cp_F/cp_L * 1/tau_p * (T_f - T_p) + dm/m * L_ev/cp_L
      ///  For detais of the analytical solution contact t.wegmann@aia.rwth-aachen.de
      /////////////////////////////////////////////////////////////////////
      if(Nu > 100 * MFloatEps) {
        const MFloat c1 =
            PI * 0.5 * (m_diameter + oldDiameter) * Nu * lamda_mix * dt / (mass * cp_liquid) / (s_Pr * s_Re);
        MFloat c2 =
            T_fluid
            - m_dM * lH_ev / (Nu * lamda_mix * PI * 0.5 * (m_diameter + oldDiameter)) * (s_Re * s_Pr * gammaMinusOne);
        m_temperature = c2 + (oldTemperature - c2) * exp(-c1);
      } else {
        //////////////////////////////////////////////////////////////////
        /// NOTE: the type of differential equation changes for Nu -> 0
        ///       and is merely: dT/dt = dm/m * L_ev/cp_L
        //////////////////////////////////////////////////////////////////
        // negative sign due to changed sign of evaporation rate!
        m_temperature = oldTemperature - m_dM / (mass * cp_liquid) * lH_ev * dt * gammaMinusOne;
      }

#ifdef LPT_DEBUG
      if(std::isnan(m_temperature)) {
        cerr << "m_diam " << m_diameter << " Nu " << Nu << " thCond " << lamda_mix << endl;
        cerr << "m_t " << m_temperature << " oldT " << oldTemperature << endl;
        cerr << " Re " << Re << " Pr " << Pr << " " << T_bndry << " " << magRelVel(m_fluidVel.data(), &m_velocity[0])
             << " " << viscosity_mix << " " << m_fluidDensity << " " << 2.0 + 0.552 * sqrt(Re) * pow(Pr, 1.0 / 3.0)
             << endl;
        TERM(-1);
      }
#endif

      if(m_temperature < 0) {
        m_temperature = MFloatEps;
      }

      heatStep++;

    } while(abs(m_temperature - previousTemp) > convEvap && heatStep < 100);

  } else {
    // constant evaporation rate used for reference/validation cases
    static MFloat dm_dt = (oldMass / dt) / 100000.0;
    m_dM = dm_dt;
    if(m_dM * dt > oldMass) {
      m_dM = oldMass / dt;
    }
    m_diameter = pow((oldMass - m_dM * dt) / s_backPtr->m_material->density(m_temperature) * 6.0 / PI, 1.0 / 3.0);
  }


  if(m_temperature < 0) {
    TERMM(-1, "Invalid temperature!");
  }

  // store particle energy change in heat-flux
  if(s_backPtr->m_heatCoupling) {
    const MFloat mass = sphericalMass();
    // see "Direct numerical simulation of a confined three-dimensional gas mixing layer with one
    // evaporating hydrocarbon-droplet-laden stream" by R.S. Miller and J. Bellan
    // in J. Fluid Mech. (1999) for reference!
    m_heatFlux = mass * cp_liquid * (m_temperature - oldTemperature) / dt * 1 / s_backPtr->m_material->gammaMinusOne();


    ASSERT(!std::isnan(m_heatFlux) && !std::isinf(m_heatFlux),
           "ERROR: m_heatFlux is NaN or inf! " + to_string(m_temperature) + " " + to_string(oldTemperature) + " "
               + to_string(volume));

    if(s_backPtr->m_evaporation) {
      m_heatFlux -= (m_dM * cp_liquid * m_temperature * 1 / s_backPtr->m_material->gammaMinusOne());
    }


    ASSERT(!std::isnan(m_heatFlux) && !std::isinf(m_heatFlux), "ERROR: m_heatFlux is NaN or inf!");
  }
}


// Explicit instantiations for 2D and 3D
template class LPTSpherical<3>;
