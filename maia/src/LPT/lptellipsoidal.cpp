// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptellipsoidal.h"
#include <cmath>
#include <functional>
#include <numeric>
#include "lpt.h"

using namespace std;
using namespace maia::lpt;

// reference values used for LPT non-dimensionalisation:
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::s_lengthFactor = F1;
// L_ref / L_refGrid (length factor between the particle reference length
//                    and the grid reference length and)
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::s_Re = F1;
// reference Re-number: Re_ref = rho_ref * u_ref * L_ref / mu_ref
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::s_Pr = F1;
// reference Pr-number: Pr_ref = mu_ref * Cp_ref / lamda_ref
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::s_Sc = F1;
// reference Sc-number: Sc_ref = mu_ref / (rho_ref * D_ref)
template <MInt nDim>
std::array<MFloat, nDim> LPTEllipsoidal<nDim>::s_Frm{};
// modified reference Fr-number with : Fr_ref = u_ref / sqrt( g * L_ref)
// the modified version is Frm = 1/ (Fr_ref^2) used to reduce computations

/**
 * \brief constructor of ellipsoidal particles
 * \author Laurent Andre
 */
template <MInt nDim>
LPTEllipsoidal<nDim>::LPTEllipsoidal() {
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
    m_angularVel[i] = nan;
    m_oldAngularVel[i] = nan;
    m_angularAccel[i] = nan;
    m_oldAngularAccel[i] = nan;
  }
  for(MInt i = 0; i < 4; i++) {
    m_quaternion[i] = nan;
    m_oldQuaternion[i] = nan;
  }
}


/**
 * \brief advance to new timeStep
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPTEllipsoidal<nDim>::advanceParticle() {
  firstStep() = false;
  hasCollided() = false;
  hadWallColl() = false;
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
void LPTEllipsoidal<nDim>::coupling() {
  if(m_cellId < 0) {
    mTerm(1, AT_, "coupling of invalid particle?");
  }

  if(!s_backPtr->m_heatCoupling && !s_backPtr->m_evaporation) {
    m_redistWeight.clear();
    interpolateVelocityAndDensity(m_cellId, m_position.data(), m_fluidVel.data(), m_fluidDensity, m_redistWeight);
  }

  // check that all variables have been exchanged/set correctly
  for(MInt i = 0; i < nDim; i++) {
    ASSERT(!isnan(m_oldFluidVel[i]), "");
  }
  ASSERT(m_oldFluidDensity > 0, "");

  // check all source terms
#ifdef LPT_DEBUG
  for(MInt i = 0; i < nDim; i++) {
    if(std::isnan(m_accel[i])) {
      cerr << "Invalid motion equation! " << m_partId << " " << m_position[0] << " " << m_position[1] << " "
           << m_position[2] << endl;
      isInvalid() = true;
      return;
    }
    if(std::isnan(m_velocity[i])) {
      cerr << "Invalid motion equation2! " << m_partId << " " << m_position[0] << " " << m_position[1] << " "
           << m_position[2] << " " << m_densityRatio << " " << hadWallColl() << endl;
      isInvalid() = true;
      return;
    }
  }
#endif

  if(s_backPtr->m_momentumCoupling) momentumCoupling();
  if(s_backPtr->m_heatCoupling) heatCoupling();
  if(s_backPtr->m_massCoupling) massCoupling();
}

/**
 * \brief add the mass flux from the particle to the cell/all surrounding cells
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
void LPTEllipsoidal<nDim>::massCoupling() {
  mTerm(1, AT_, "Mass coupling not implemented for ellipsoidal particles");
}

/**
 * \brief add the momentum flux from the particle to the cell/all surrounding cells
 * \author Sven Berger, Tim Wegmann
 */
template <MInt nDim>
void LPTEllipsoidal<nDim>::momentumCoupling() {
  const MFloat relT = (s_backPtr->m_timeStep - m_creationTime) / s_backPtr->m_timeStep;

  const MFloat mass = particleMass();

  array<MFloat, nDim> force{};

  for(MInt i = 0; i < nDim; i++) {
    // momentum flux due to drag force
    const MFloat invDensityRatio = (abs(m_densityRatio) <= MFloatEps) ? 1.0 : 1.0 / m_densityRatio;
    force[i] = mass * (m_accel[i] - ((1.0 - invDensityRatio) * s_Frm[i])) * relT;
  }

  // work done by the momentum
  const MFloat work = std::inner_product(force.begin(), force.end(), m_oldVel.begin(), 0.0);


  if(!s_backPtr->m_couplingRedist) {
    for(MInt i = 0; i < nDim; i++) {
      s_backPtr->a_momentumFlux(m_cellId, i) += (force[i]);
    }
    s_backPtr->a_workFlux(m_cellId) += work;

  } else {
    ASSERT(m_redistWeight.size() == m_neighborList.size(),
           "Invalid " + to_string(m_redistWeight.size()) + "!=" + to_string(m_neighborList.size()));
    for(MInt n = 0; n < (signed)m_neighborList.size(); n++) {
      // a weightHeat is used to distribute the source term to multiple cells...
      for(MInt i = 0; i < nDim; i++) {
        s_backPtr->a_momentumFlux(m_neighborList[n], i) += m_redistWeight.at(n) * (force[i]);
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
void LPTEllipsoidal<nDim>::heatCoupling() {
  mTerm(1, AT_, "Heat coupling not tested for ellipsoidal particles");
}


template <MInt nDim>
void LPTEllipsoidal<nDim>::motionEquation() {
  const MFloat dt = s_backPtr->m_timeStep - m_creationTime;
  ASSERT(dt > 0, "Timestep < 0");

  // NOTE: at this point old means 1 time-Step ago and the current values will be updated below!
  m_oldPos = m_position;
  m_oldVel = m_velocity;
  m_oldAccel = m_accel;
  m_oldAngularVel = m_angularVel;
  m_oldAngularAccel = m_angularAccel;
  m_oldQuaternion = m_quaternion;
  m_oldCellId = m_cellId;

  if(firstStep()) {
    for(MInt i = 0; i < nDim; i++) {
      m_accel[i] = s_Frm[i];
      m_oldAccel[i] = s_Frm[i];
    }
  }

  for(MInt i = 0; i < nDim; i++) {
    ASSERT(!isnan(m_oldFluidVel[i]), "");
  }
  ASSERT(m_oldFluidDensity > 0, "");

  // reduce variance between runs (round to 9 decimal places)
  // TODO-timw labels:LPT,totest check if rounding is still necessary?
  const MFloat invDensityRatio =
      F1 / (static_cast<MFloat>(ceil(s_backPtr->m_material->density() / m_oldFluidDensity * 1E9)) / 1E9);

  m_densityRatio = 1 / invDensityRatio;

  /// 1. Predictor step
  array<MFloat, nDim> predictedPos{};
  array<MFloat, nDim> predictedVel{};
  array<MFloat, nDim> predictedAngularVel{};
  array<MFloat, nDim + 1> predictedQuaternion{};

  MFloat q[3];
  MFloat q0[3];
  MFloat tq[3];
  MFloat tq0[3];
  MFloat w = m_oldQuaternion[0];
  MFloat x = m_oldQuaternion[1];
  MFloat y = m_oldQuaternion[2];
  MFloat z = m_oldQuaternion[3];
  for(MInt i = 0; i < 3; i++) {
    q0[i] = m_oldAngularVel[i];
  }

  // predict position, velocity, orientation and angular velocity
  for(MInt i = 0; i < nDim; i++) {
    predictedVel[i] = m_oldVel[i] + dt * m_oldAccel[i];
    predictedPos[i] = m_oldPos[i] + dt * (m_oldVel[i]) + F1B2 * POW2(dt) * (m_oldAccel[i]);
  }
  predictedQuaternion[0] = w + F1B2 * dt * (-x * q0[0] - y * q0[1] - z * q0[2]);
  predictedQuaternion[1] = x + F1B2 * dt * (w * q0[0] - z * q0[1] + y * q0[2]);
  predictedQuaternion[2] = y + F1B2 * dt * (z * q0[0] + w * q0[1] - x * q0[2]);
  predictedQuaternion[3] = z + F1B2 * dt * (-y * q0[0] + x * q0[1] + w * q0[2]);
  for(MInt i = 0; i < 3; i++) {
    predictedAngularVel[i] = m_oldAngularVel[i] + dt * m_oldAngularAccel[i];
  }

  /// 2. get velocity, density and velocity magnitude at the predicted position
  //  find the predicted CellId at the predicted position around the previous cellId
  const MInt predictedCellId = s_backPtr->grid().findContainingLeafCell(&predictedPos[0], m_cellId);

  ASSERT(m_cellId > -1, "Invalid cell!");
  ASSERT(predictedCellId > -1, "Invalid cell!");
  ASSERT(s_backPtr->c_isLeafCell(predictedCellId), "Invalid cell!");

  vector<MFloat> predictedWeights(0);
  array<MFloat, nDim> fluidVel{};
  array<MFloat, nDim * nDim> fluidVelGradient{};
  array<MFloat, nDim * nDim> fluidVelGradientHat{};
  std::fill(fluidVel.begin(), fluidVel.end(), std::numeric_limits<MFloat>::quiet_NaN());
  std::fill(fluidVelGradient.begin(), fluidVelGradient.end(), std::numeric_limits<MFloat>::quiet_NaN());
  std::fill(fluidVelGradientHat.begin(), fluidVelGradientHat.end(), std::numeric_limits<MFloat>::quiet_NaN());
  interpolateVelocityAndVelocitySlopes(predictedCellId, predictedPos.data(), fluidVel.data(), fluidVelGradient.data(),
                                       predictedWeights);

  // calculate gradients according to principle axes of ellipsoid
  MFloatScratchSpace R(3, 3, AT_, "R");
  maia::math::computeRotationMatrix(R, &(m_quaternion[0]));
  for(MInt i = 0; i < nDim; i++) {
    maia::math::matrixVectorProduct(&fluidVelGradientHat[nDim * i], R, &fluidVelGradient[nDim * i]);
  }

  m_fluidVelMag = sqrt(std::inner_product(fluidVel.begin(), fluidVel.end(), fluidVel.begin(), F0));

  const MFloat T = s_backPtr->a_fluidTemperature(m_cellId);

  MFloat fluidViscosity = s_backPtr->m_material->dynViscosityFun(T);

  /// 3. corrector step
  if(s_backPtr->m_dragModelType == 0) {
    // no drag or neglible small (are going to be deleted...)
    for(MInt i = 0; i < nDim; i++) {
      m_velocity[i] = predictedVel[i];
      m_position[i] = predictedPos[i];
      m_angularVel[i] = predictedAngularVel[i];
      m_quaternion[i] = predictedQuaternion[i];
      m_accel[i] = m_oldAccel[i];
      m_angularAccel[i] = m_oldAngularAccel[i];
    }
  } else {
    switch(s_backPtr->m_motionEquationType) {
      case 0: {
        // particle moves with fluid velocity, no rotational motion
        for(MInt i = 0; i < nDim; i++) {
          m_position[i] = m_oldPos[i] + F1B2 * dt * (fluidVel[i] + m_oldVel[i]) * s_lengthFactor;
          m_velocity[i] = fluidVel[i];
          m_accel[i] = (m_velocity[i] - m_oldVel[i]) / dt;
          m_angularVel[i] = F0;
          m_angularAccel[i] = F0;
        }
        break;
      }
      case 1: {
        // spherical case copied from lptspherical
        // serves for debugging to have a version independant of the strain, vorticity
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
        // case derived for creeping flow
        const MFloat fdtB2 = F1B2 * dt;
        MFloatScratchSpace K(3, 3, AT_, "K");
        MFloatScratchSpace W(4, 4, AT_, "W");
        MFloatScratchSpace iW(4, 4, AT_, "iW");
        MFloatScratchSpace rhs(4, AT_, "rhs");

        MFloat vrel[3];
        MFloat vrelhat[3];
        MFloat tmp[3];
        MFloat vort[3];
        MFloat strain[3];

        const MFloat beta = m_aspectRatio;
        const MFloat beta2 = POW2(beta);
        const MFloat fTauP = fParticleRelTime(fluidViscosity) / s_Re;

        const MFloat linAccFac = (8.0 / 3.0) * pow(beta, F2B3) * fTauP;
        const MFloat angAccFac = (40.0 / 9.0) * pow(beta, F2B3) * fTauP;

        // integrate angular velocity
        const MInt maxit = 100;
        MFloat delta = F1;
        MInt it = 0;

        // set temporary variables for angular velocity (q) and acceleration (tq)
        for(MInt i = 0; i < 3; i++) {
          tq0[i] = m_oldAngularAccel[i];
          q0[i] = m_oldAngularVel[i];
          tq[i] = tq0[i];
          q[i] = q0[i];
        }

        // compute vorticity and strain
        for(MInt i = 0; i < 3; i++) {
          MInt id0 = (i + 1) % 3;
          MInt id1 = (id0 + 1) % 3;
          strain[i] = F1B2 * (fluidVelGradientHat[3 * id1 + id0] + fluidVelGradientHat[3 * id0 + id1]);
          vort[i] = F1B2 * (fluidVelGradientHat[3 * id1 + id0] - fluidVelGradientHat[3 * id0 + id1]);
        }

        W(0, 0) = F1 + fdtB2 * angAccFac / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
        W(1, 1) = F1 + fdtB2 * angAccFac / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
        W(2, 2) = F1 + fdtB2 * angAccFac / (F2 * m_shapeParams[1]);
        W(2, 0) = F0;
        W(2, 1) = F0;

        // Newton iterations
        while(delta > 1e-10 && it < maxit) {
          W(0, 1) = -fdtB2 * q[2] * (beta2 - F1) / (beta2 + F1);
          W(0, 2) = -fdtB2 * q[1] * (beta2 - F1) / (beta2 + F1);
          W(1, 0) = -fdtB2 * q[2] * (F1 - beta2) / (beta2 + F1);
          W(1, 2) = -fdtB2 * q[0] * (F1 - beta2) / (beta2 + F1);

          maia::math::invert(W, iW, 3, 3);

          // compute angular acceleration tq
          tq[0] = q[1] * q[2] * (beta2 - F1) / (beta2 + F1)
                  + angAccFac * (((F1 - beta2) / (F1 + beta2)) * strain[0] + vort[0] - q[0])
                        / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
          tq[1] = q[0] * q[2] * (F1 - beta2) / (beta2 + F1)
                  + angAccFac * (((beta2 - F1) / (F1 + beta2)) * strain[1] + vort[1] - q[1])
                        / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
          tq[2] = angAccFac * (vort[2] - q[2]) / (F2 * m_shapeParams[1]);

          for(MInt i = 0; i < 3; i++)
            rhs[i] = q[i] - q0[i] - fdtB2 * (tq0[i] + tq[i]);

          delta = F0;
          for(MInt i = 0; i < 3; i++) {
            const MFloat qq = q[i];
            for(MInt j = 0; j < 3; j++) {
              q[i] -= iW(i, j) * rhs(j);
            }
            delta = mMax(delta, fabs(q[i] - qq));
          }
          it++;
        }

        if(it >= maxit || it <= 0) {
          cerr << "Newton iterations did not converge " << m_partId << endl;
        }

        for(MInt i = 0; i < 3; i++) {
          m_angularVel[i] = q[i];
          m_angularAccel[i] = tq[i];
        }

        // integrate quaternions
        w = m_oldQuaternion[0];
        x = m_oldQuaternion[1];
        y = m_oldQuaternion[2];
        z = m_oldQuaternion[3];

        W(0, 0) = F0;
        W(0, 1) = -q[0];
        W(0, 2) = -q[1];
        W(0, 3) = -q[2];
        W(1, 0) = q[0];
        W(1, 1) = F0;
        W(1, 2) = q[2];
        W(1, 3) = -q[1];
        W(2, 0) = q[1];
        W(2, 1) = -q[2];
        W(2, 2) = F0;
        W(2, 3) = q[0];
        W(3, 0) = q[2];
        W(3, 1) = q[1];
        W(3, 2) = -q[0];
        W(3, 3) = F0;

        for(MInt i = 0; i < 4; i++) {
          for(MInt j = 0; j < 4; j++) {
            W(i, j) = -F1B2 * fdtB2 * W(i, j);
          }
          W(i, i) = F1;
        }

        rhs(0) = w + F1B2 * fdtB2 * (-x * q0[0] - y * q0[1] - z * q0[2]);
        rhs(1) = x + F1B2 * fdtB2 * (w * q0[0] - z * q0[1] + y * q0[2]);
        rhs(2) = y + F1B2 * fdtB2 * (z * q0[0] + w * q0[1] - x * q0[2]);
        rhs(3) = z + F1B2 * fdtB2 * (-y * q0[0] + x * q0[1] + w * q0[2]);

        maia::math::invert(W, iW, 4, 4);

        for(MInt i = 0; i < 4; i++) {
          m_quaternion[i] = F0;
          for(MInt j = 0; j < 4; j++) {
            m_quaternion[i] += iW(i, j) * rhs(j);
          }
        }

        MFloat abs = F0;
        for(MInt i = 0; i < 4; i++)
          abs += POW2(m_quaternion[i]);
        for(MInt i = 0; i < 4; i++)
          m_quaternion[i] /= sqrt(abs);

        // integrate linear motion
        K.fill(F0);
        K(0, 0) = F1 / (m_shapeParams[0] / POW2(m_semiMinorAxis) + m_shapeParams[1]);
        K(1, 1) = K(0, 0);
        K(2, 2) = F1 / (m_shapeParams[0] / POW2(m_semiMinorAxis) + beta2 * m_shapeParams[3]);
        maia::math::computeRotationMatrix(R, &(m_quaternion[0]));

        for(MInt i = 0; i < nDim; i++) {
          vrel[i] = fluidVel[i] - predictedVel[i];
        }

        maia::math::matrixVectorProduct(vrelhat, R, vrel); // principal axes frame relative velocity
        maia::math::matrixVectorProduct(tmp, K, vrelhat);
        maia::math::matrixVectorProductTranspose(vrel, R, tmp);

        for(MInt i = 0; i < nDim; i++) {
          m_accel[i] = linAccFac * vrel[i] + (1.0 - invDensityRatio) * s_Frm[i];
          m_velocity[i] = m_oldVel[i] + dt * F1B2 * (m_accel[i] + m_oldAccel[i]);
          m_position[i] = m_oldPos[i] + dt * F1B2 * (m_velocity[i] + m_oldVel[i]) * s_lengthFactor
                          + F1B4 * POW2(dt) * (m_accel[i] + m_oldAccel[i]);
        }
        break;
      }
      case 3: {
        // case using new correlations functions by Konstantin (2020)
        const MFloat fdtB2 = F1B2 * dt;
        MFloatScratchSpace K(3, 3, AT_, "K");
        MFloatScratchSpace W(4, 4, AT_, "W");
        MFloatScratchSpace iW(4, 4, AT_, "iW");
        MFloatScratchSpace rhs(4, AT_, "rhs");

        MFloat vrel[3];
        MFloat vrelhat[3];
        MFloat tmp[3];
        MFloat vort[3];
        MFloat strain[3];

        const MFloat beta = m_aspectRatio;
        const MFloat beta2 = POW2(beta);
        const MFloat eqRadius = F1B2 * m_eqDiameter;
        const MFloat fTauP = fParticleRelTime(fluidViscosity) / s_Re;
        const MFloat magRelVeloc = magRelVel(&fluidVel[0], &predictedVel[0]);
        const MFloat Rep = particleRe(magRelVeloc, m_oldFluidDensity, fluidViscosity) * s_Re;

        const MFloat angAccFac = (40.0 / 9.0) * pow(beta, F2B3) * fTauP;
        const MFloat pitchingTorqueFac = (15.0 / 8.0) * invDensityRatio * pow(beta, F2B3) / POW2(eqRadius);
        const MFloat addFacCd =
            48.0 * fluidViscosity / (s_Re * s_backPtr->m_material->density() * POW2(m_eqDiameter)) * pow(beta, F2B3);
        const MFloat momI[3] = {1.0 + POW2(beta), 1.0 + POW2(beta), 2.0}; // factors for moments of inertia

        maia::math::computeRotationMatrix(R, &(m_quaternion[0]));

        for(MInt i = 0; i < nDim; i++) {
          vrel[i] = fluidVel[i] - predictedVel[i];
        }

        maia::math::matrixVectorProduct(vrelhat, R, vrel); // principal axes frame relative velocity
        maia::math::matrixVectorProduct(tmp, K, vrelhat);
        maia::math::matrixVectorProductTranspose(vrel, R, tmp);

        const MFloat magRelVelHat = sqrt(POW2(vrelhat[0]) + POW2(vrelhat[1]) + POW2(vrelhat[2]));

        const auto nan = std::numeric_limits<MFloat>::quiet_NaN();
        MFloat inclinationAngle = nan;
        std::array<MFloat, 3> accHat = {nan, nan, nan};
        std::array<MFloat, 3> accDragHat = {nan, nan, nan};
        std::array<MFloat, 3> accLiftHat = {nan, nan, nan};
        std::array<MFloat, 3> dirDHat = {nan, nan, nan};
        std::array<MFloat, 3> dirTHat = {nan, nan, nan};
        std::array<MFloat, 3> dirLHat = {nan, nan, nan};
        std::array<MFloat, 3> dirDCrossZ = {nan, nan, nan};
        std::array<MFloat, 3> dirZHat{nan, nan, nan};

        // dirZHat
        for(MInt i = 0; i < nDim; i++) {
          i == m_particleMajorAxis ? dirZHat[i] = F1 : dirZHat[i] = F0;
        }
        // dirDhat
        for(MInt i = 0; i < nDim; i++) {
          dirDHat[i] = vrelhat[i] / magRelVelHat;
        }
        // dirDcrossZ = dirDhat x dirZhat
        maia::math::cross(&dirDHat[0], &dirZHat[0], &dirDCrossZ[0]);
        const MFloat dotProduct = std::inner_product(dirDHat.begin(), dirDHat.end(), dirZHat.begin(), F0);
        inclinationAngle = acos(fabs(dotProduct));
        MInt sgn = (dotProduct > 0 ? 1 : -1);
        for(MInt i = 0; i < nDim; i++) {
          dirTHat[i] = dirDCrossZ[i] / sqrt(POW2(dirDCrossZ[0]) + POW2(dirDCrossZ[1]) + POW2(dirDCrossZ[2])) * sgn;
        }
        // dirLHat = dirDhat x dirThat
        maia::math::cross(&dirDHat[0], &dirTHat[0], &dirLHat[0]);

        // integrate angular velocity
        const MInt maxit = 100;
        MFloat delta = F1;
        MInt it = 0;

        W(0, 0) = F1 + fdtB2 * angAccFac / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
        W(1, 1) = F1 + fdtB2 * angAccFac / (m_shapeParams[1] + beta2 * m_shapeParams[3]);
        W(2, 2) = F1 + fdtB2 * angAccFac / (F2 * m_shapeParams[1]);
        W(2, 0) = F0;
        W(2, 1) = F0;

        // set temporary variables for angular velocity (q) and acceleration (tq)
        for(MInt i = 0; i < 3; i++) {
          tq0[i] = m_oldAngularAccel[i];
          q0[i] = m_oldAngularVel[i];
          tq[i] = tq0[i];
          q[i] = q0[i];
        }

        // compute pitching torque unless torque modelling is turned off
        MFloat pitch[3] = {F0};
        if(s_backPtr->m_torqueModelType > 0 && Rep > 1e-8) {
          MFloat CT = torqueFactor(Rep, beta, inclinationAngle);
          for(MInt i = 0; i < 3; i++) {
            pitch[i] = pitchingTorqueFac * POW2(magRelVelHat) * CT / momI[i] * dirTHat[i];
            if(pitchingTorqueFac <= F0 || CT <= F0 || momI[i] <= F0)
              cout << "Calc pitching torque " << pitchingTorqueFac << " " << CT << " " << momI[i] << endl;
          }
        }

        // compute vorticity and strain
        for(MInt i = 0; i < 3; i++) {
          MInt id0 = (i + 1) % 3;
          MInt id1 = (id0 + 1) % 3;
          strain[i] = F1B2 * (fluidVelGradientHat[3 * id1 + id0] + fluidVelGradientHat[3 * id0 + id1]);
          vort[i] = F1B2 * (fluidVelGradientHat[3 * id1 + id0] - fluidVelGradientHat[3 * id0 + id1]);
        }

        // Newton iterations
        while(delta > 1e-10 && it < maxit) {
          W(0, 1) = -fdtB2 * q[2] * (beta2 - F1) / (beta2 + F1);
          W(0, 2) = -fdtB2 * q[1] * (beta2 - F1) / (beta2 + F1);
          W(1, 0) = -fdtB2 * q[2] * (F1 - beta2) / (beta2 + F1);
          W(1, 2) = -fdtB2 * q[0] * (F1 - beta2) / (beta2 + F1);

          maia::math::invert(W, iW, 3, 3);

          tq[0] = q[1] * q[2] * (beta2 - F1) / (beta2 + F1)
                  + angAccFac * (((F1 - beta2) / (F1 + beta2)) * strain[0] + vort[0] - q[0])
                        / (m_shapeParams[1] + beta2 * m_shapeParams[3])
                  + pitch[0];
          tq[1] = q[0] * q[2] * (F1 - beta2) / (beta2 + F1)
                  + angAccFac * (((beta2 - F1) / (F1 + beta2)) * strain[1] + vort[1] - q[1])
                        / (m_shapeParams[1] + beta2 * m_shapeParams[3])
                  + pitch[1];
          tq[2] = angAccFac * (vort[2] - q[2]) / (F2 * m_shapeParams[1]) + pitch[2];

          for(MInt i = 0; i < 3; i++)
            rhs[i] = q[i] - q0[i] - fdtB2 * (tq0[i] + tq[i]);

          delta = F0;
          for(MInt i = 0; i < 3; i++) {
            const MFloat qq = q[i];
            for(MInt j = 0; j < 3; j++) {
              q[i] -= iW(i, j) * rhs(j);
            }
            delta = mMax(delta, fabs(q[i] - qq));
          }
          it++;
        }

        if(it >= maxit || it <= 0) {
          cerr << "Newton iterations did not converge " << m_partId << endl;
        }

        for(MInt i = 0; i < 3; i++) {
          m_angularVel[i] = q[i];
          m_angularAccel[i] = tq[i];
        }

        // integrate quaternions
        w = m_oldQuaternion[0];
        x = m_oldQuaternion[1];
        y = m_oldQuaternion[2];
        z = m_oldQuaternion[3];

        W(0, 0) = F0;
        W(0, 1) = -q[0];
        W(0, 2) = -q[1];
        W(0, 3) = -q[2];
        W(1, 0) = q[0];
        W(1, 1) = F0;
        W(1, 2) = q[2];
        W(1, 3) = -q[1];
        W(2, 0) = q[1];
        W(2, 1) = -q[2];
        W(2, 2) = F0;
        W(2, 3) = q[0];
        W(3, 0) = q[2];
        W(3, 1) = q[1];
        W(3, 2) = -q[0];
        W(3, 3) = F0;

        for(MInt i = 0; i < 4; i++) {
          for(MInt j = 0; j < 4; j++) {
            W(i, j) = -F1B2 * fdtB2 * W(i, j);
          }
          W(i, i) = F1;
        }

        rhs(0) = w + F1B2 * fdtB2 * (-x * q0[0] - y * q0[1] - z * q0[2]);
        rhs(1) = x + F1B2 * fdtB2 * (w * q0[0] - z * q0[1] + y * q0[2]);
        rhs(2) = y + F1B2 * fdtB2 * (z * q0[0] + w * q0[1] - x * q0[2]);
        rhs(3) = z + F1B2 * fdtB2 * (-y * q0[0] + x * q0[1] + w * q0[2]);

        maia::math::invert(W, iW, 4, 4);

        for(MInt i = 0; i < 4; i++) {
          m_quaternion[i] = F0;
          for(MInt j = 0; j < 4; j++) {
            m_quaternion[i] += iW(i, j) * rhs(j);
          }
        }

        MFloat abs = F0;
        for(MInt i = 0; i < 4; i++)
          abs += POW2(m_quaternion[i]);
        for(MInt i = 0; i < 4; i++)
          m_quaternion[i] /= sqrt(abs);

        // integrate linear motion
        K.fill(F0);
        K(0, 0) = F1 / (m_shapeParams[0] / POW2(m_semiMinorAxis) + m_shapeParams[1]);
        K(1, 1) = K(0, 0);
        K(2, 2) = F1 / (m_shapeParams[0] / POW2(m_semiMinorAxis) + beta2 * m_shapeParams[3]);

        // new drag correlations
        MFloatScratchSpace accDrag(3, AT_, "accDrag");
        accDrag.fill(F0);
        MFloat facCL = liftFactor(Rep, beta, inclinationAngle, K(0, 0), K(2, 2));
        MFloat facCD = dragFactor(Rep, beta, inclinationAngle, K(0, 0), K(2, 2));
        for(MInt i = 0; i < nDim; i++) {
          accLiftHat[i] = addFacCd * facCL * magRelVelHat * dirLHat[i];
          accDragHat[i] = addFacCd * facCD * magRelVelHat * dirDHat[i];
          accHat[i] = accDragHat[i] + accLiftHat[i];
        }
        maia::math::matrixVectorProductTranspose(&accDrag[0], R, &accHat[0]);

        for(MInt i = 0; i < nDim; i++) {
          m_accel[i] = accDrag[i] + (1.0 - invDensityRatio) * s_Frm[i];
          m_velocity[i] = m_oldVel[i] + dt * F1B2 * (m_accel[i] + m_oldAccel[i]);
          m_position[i] = m_oldPos[i] + dt * F1B2 * (m_velocity[i] + m_oldVel[i]) * s_lengthFactor
                          + F1B4 * POW2(dt) * (m_accel[i] + m_oldAccel[i]);
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


/**
 * \brief Calculates the drag factor for current particle with given Re number, aspect ratio and inclination angle
 * \author Laurent Andre
 * \date  August 2022
 */
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::dragFactor(const MFloat partRe, const MFloat beta, const MFloat inclinationAngle,
                                        const MFloat K0, const MFloat K2) {
  MFloat result = NAN;
  switch(s_backPtr->m_dragModelType) {
    case 0: { // drag force switched off
      result = F0;
      break;
    }
    case 1: { // linear Stokes drag
      result = F1;
      break;
    }
    case 2: { // nonlinear drag according to Schiller & Naumann
      result = 1.0 + 0.15 * pow(partRe, 0.687);
      break;
    }
    case 3: { // new correlations by Konstantin (2020)
      std::array<MFloat, 8> dCorrs{-0.007, 1.0, 1.17, -0.07, 0.047, 1.14, 0.7, -0.008};
      MFloat fd0 = F1;
      MFloat fd90 = F1;
      fd0 = 1.0 + 0.15 * pow(partRe, 0.687)
            + dCorrs[0] * pow(std::log(beta), dCorrs[1]) * pow(partRe, dCorrs[2] + dCorrs[3] * std::log(beta));
      fd90 = 1.0 + 0.15 * pow(partRe, 0.687)
             + dCorrs[4] * pow(log(beta), dCorrs[5]) * pow(partRe, dCorrs[6] + dCorrs[7] * log(beta));

      MFloat fd0K2 = fd0 * K2;
      MFloat f90K0 = fd90 * K0;
      result = fd0K2 + (f90K0 - fd0K2) * POW2(sin(inclinationAngle));
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown particle drag method");
    }
  }
  return result;
}

/**
 * \brief Calculates the lift factor for current particle with given Re number, aspect ratio, and inclination angle
 * \author Laurent Andre
 * \date  August 2022
 */
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::liftFactor(const MFloat partRe, const MFloat beta, const MFloat inclinationAngle,
                                        const MFloat K0, const MFloat K2) {
  MFloat result = NAN;
  switch(s_backPtr->m_liftModelType) {
    case 0: { // lift force switched off
      result = F0;
      break;
    }
    case 1: { // linear
      result = F1;
      break;
    }
    case 2: { // new correlations by Konstantin (2020)
      std::array<MFloat, 6> lCorrs{0.01, 0.86, 1.77, 0.34, 0.88, -0.05};
      MFloat flMax = F1;
      MFloat flShift = F1;
      if(partRe > 1) flShift += lCorrs[0] * pow(log(beta), lCorrs[1]) * pow(log(partRe), lCorrs[2]);
      MFloat psi = F1B2 * M_PI * pow(inclinationAngle / M_PI * F2, flShift);
      flMax = F1 + lCorrs[3] * pow(partRe, lCorrs[4] + lCorrs[5]);
      MFloat ClMax = F1B2 * flMax * (K0 - K2);
      MFloat Cl = F2 * sin(psi) * cos(psi) * ClMax;
      result = Cl;
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown particle drag method");
    }
  }
  return result;
}

/**
 * \brief Calculates the torque factor for current particle with given Re number, aspect ratio and inclination angle
 * \author Laurent Andre
 * \date  August 2022
 */
template <MInt nDim>
MFloat LPTEllipsoidal<nDim>::torqueFactor(const MFloat partRe, const MFloat beta, const MFloat inclinationAngle) {
  MFloat result = NAN;
  switch(s_backPtr->m_torqueModelType) {
    case 0: { // torque switched off
      result = F0;
      break;
    }
    case 1: { // linear
      result = F1;
      break;
    }
    case 2: { // new correlations by Konstantin (2020)
      std::array<MFloat, 7> tCorrs{0.931, 0.675, 0.162, 0.657, 2.77, 0.178, 0.177};
      MFloat CTMax = tCorrs[0] * pow(log(beta), tCorrs[1]) / pow(partRe, tCorrs[2])
                     + tCorrs[3] * pow(log(beta), tCorrs[4]) / pow(partRe, tCorrs[5] + tCorrs[6] * log(beta));
      MFloat CT = F2 * sin(inclinationAngle) * cos(inclinationAngle) * CTMax;
      result = CT;
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown particle torque method");
    }
  }
  return result;
}


/**
 * \brief Calcultes shape parameters for current particle (see Siewert, 2014)
 * \author Laurent Andre
 * \date  August 2022
 */
template <MInt nDim>
void LPTEllipsoidal<nDim>::initShapeParams() {
  MFloat beta = m_aspectRatio;
  MFloat beta2 = beta * beta;
  if(fabs(beta - F1) < 1e-14) {
    // spherical
    m_shapeParams[0] = F2 * POW2(m_semiMinorAxis);
    m_shapeParams[1] = F2B3;
    m_shapeParams[2] = F2B3;
    m_shapeParams[3] = F2B3;
  } else if(beta > F1) {
    // prolate
    MFloat kappa = log((beta - sqrt(beta2 - F1)) / (beta + sqrt(beta2 - F1)));
    m_shapeParams[0] = -beta * POW2(m_semiMinorAxis) * kappa / sqrt(beta2 - F1);
    m_shapeParams[1] = beta2 / (beta2 - F1) + beta * kappa / (F2 * pow(beta2 - F1, F3B2));
    m_shapeParams[2] = beta2 / (beta2 - F1) + beta * kappa / (F2 * pow(beta2 - F1, F3B2));
    m_shapeParams[3] = -F2 / (beta2 - F1) - beta * kappa / (pow(beta2 - F1, F3B2));
  } else if(beta < F1) {
    // oblate
    MFloat kappa = F2 * atan2(beta, sqrt(F1 - beta2));
    m_shapeParams[0] = beta * POW2(m_semiMinorAxis) * (PI - kappa) / sqrt(F1 - beta2);
    m_shapeParams[1] = -beta * (kappa - PI + F2 * beta * sqrt(F1 - beta2)) / (F2 * pow(F1 - beta2, F3B2));
    m_shapeParams[1] = -beta * (kappa - PI + F2 * beta * sqrt(F1 - beta2)) / (F2 * pow(F1 - beta2, F3B2));
    m_shapeParams[1] = (beta * kappa - beta * PI + F2 * sqrt(F1 - beta2)) / (pow(F1 - beta2, F3B2));
  } else {
    // general triaxial ellipsoid
    mTerm(1, "Generalize here");
  }
}


template <MInt nDim>
void LPTEllipsoidal<nDim>::energyEquation() {
  mTerm(1, AT_, "Not yet implemented");
}


/**
 * \brief Calculate the orientation of the major axis in the world coordinate system.
 *
 * \param majorAxis - output array to store the major axis
 */
template <MInt nDim>
void LPTEllipsoidal<nDim>::calculateMajorAxisOrientation(MFloat* majorAxis) {
  std::array<MFloat, nDim> majorAxisParticleFixed{}; // major axis in particle fixed coordinate system
  for(MInt i = 0; i < nDim; i++) {
    i == m_particleMajorAxis ? majorAxisParticleFixed[i] = maxParticleRadius() : majorAxisParticleFixed[i] = F0;
  }
  MFloatScratchSpace R(3, 3, AT_, "R");
  maia::math::computeRotationMatrix(R, &m_quaternion[0]);
  maia::math::matrixVectorProductTranspose(majorAxis, R, majorAxisParticleFixed.begin());
}

// Explicit instantiations
template class LPTEllipsoidal<3>;
