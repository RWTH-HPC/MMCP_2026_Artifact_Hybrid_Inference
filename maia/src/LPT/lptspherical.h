// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _PARTICLESPHERICAL_H
#define _PARTICLESPHERICAL_H

#include <cmath>
#include "UTIL/materialstate.h"
#include "lptbase.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim>
class LPTBase;

template <MInt nDim>
class LPTSpherical : public LPTBase<nDim> {

 public:
  using LPTBase<nDim>::m_oldPos;
  using LPTBase<nDim>::m_position;
  using LPTBase<nDim>::m_velocity;
  using LPTBase<nDim>::m_accel;
  using LPTBase<nDim>::m_oldVel;
  using LPTBase<nDim>::m_densityRatio;
  using LPTBase<nDim>::s_backPtr;
  using LPTBase<nDim>::m_cellId;
  using LPTBase<nDim>::m_oldCellId;
  using LPTBase<nDim>::m_partId;
  using LPTBase<nDim>::m_creationTime;

  using LPTBase<nDim>::firstStep;
  using LPTBase<nDim>::isWindow;
  using LPTBase<nDim>::reqSend;
  using LPTBase<nDim>::reqBroadcast;
  using LPTBase<nDim>::isInvalid;
  using LPTBase<nDim>::hasCollided;
  using LPTBase<nDim>::hadWallColl;
  using LPTBase<nDim>::toBeDeleted;
  using LPTBase<nDim>::wasSend;
  using LPTBase<nDim>::toBeRespawn;
  using LPTBase<nDim>::fullyEvaporated;
  using LPTBase<nDim>::particleWallCollision;
  using LPTBase<nDim>::wallParticleCollision;

  using LPTBase<nDim>::interpolateAndCalcWeights;
  using LPTBase<nDim>::getNghbrList;
  using LPTBase<nDim>::checkCellChange;
  using LPTBase<nDim>::updateProperties;
  using LPTBase<nDim>::m_neighborList;

  LPTSpherical();
  ~LPTSpherical() override = default;

  LPTSpherical(const LPTSpherical&) = default;
  LPTSpherical& operator=(const LPTSpherical&) = default;

  /// fluid velocity
  std::array<MFloat, nDim> m_fluidVel{};
  /// fluid density
  MFloat m_fluidDensity = std::numeric_limits<MFloat>::quiet_NaN();

  /// particle acceleration of the last time step
  std::array<MFloat, nDim> m_oldAccel{};
  /// fluid velocity of the last time step
  std::array<MFloat, nDim> m_oldFluidVel{};

  /// old fluid density
  MFloat m_oldFluidDensity = std::numeric_limits<MFloat>::quiet_NaN();

  /// particle diameter
  MFloat m_diameter = std::numeric_limits<MFloat>::quiet_NaN();
  /// temperature
  MFloat m_temperature = std::numeric_limits<MFloat>::quiet_NaN();
  /// time since last breakUp
  MFloat m_breakUpTime = std::numeric_limits<MFloat>::quiet_NaN();
  /// shedding diameter
  MFloat m_shedDiam = std::numeric_limits<MFloat>::quiet_NaN();
  /// parceled particles
  MInt m_noParticles = -1;
  /// mass evaporation rate
  MFloat m_dM = std::numeric_limits<MFloat>::quiet_NaN();
  /// fluid velocity magnitude
  MFloat m_fluidVelMag = std::numeric_limits<MFloat>::quiet_NaN();
  /// heat flux to cell
  MFloat m_heatFlux = std::numeric_limits<MFloat>::quiet_NaN();

  static constexpr MInt s_floatElements = 9 + 4 * nDim;
  static constexpr MInt s_intElements = 1;

  static MFloat s_lengthFactor;
  static MFloat s_Re;
  static MFloat s_Pr;
  static MFloat s_Sc;
  static MFloat s_We;
  static std::array<MFloat, nDim> s_Frm;

  // store neighbor cells including diagonals
  // std::vector<MInt> m_neighborList;
  // weights of the matching neighbors
  std::vector<MFloat> m_redistWeight{};

  void motionEquation() override;
  void energyEquation() override;
  void coupling() override;
  void advanceParticle() override;
  void resetWeights() override {
    m_neighborList.clear();
    this->template interpolateAndCalcWeights<0, 0>(m_cellId, m_position.data(), nullptr, m_redistWeight);
  };

  //  void timeStepRK();
  MFloat dragFactor(const MFloat partRe);

  template <MInt>
  friend std::ostream& operator<<(std::ostream& os, const LPTSpherical& m);

  /// \brief Calculate the current volume
  inline MFloat sphericalVolume() const { return 1.0 / 6.0 * PI * POW3(m_diameter); }

  /// \brief Calculate the current mass
  inline MFloat sphericalMass() const {
    return 1.0 / 6.0 * PI * POW3(m_diameter) * s_backPtr->m_material->density(m_temperature);
  }

  /// \brief Calculate the magnitude of the relative velocity of the two given velocity vectors
  inline MFloat magRelVel(const MFloat* const velocity1, const MFloat* const velocity2) const {
    MFloat magnitude = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      magnitude += POW2(velocity1[i] - velocity2[i]);
    }
    return sqrt(magnitude);
  }

  /// \brief Calculate the magnitude of the relative velocity of the two given velocity vectors
  inline MFloat magVel() const {
    MFloat magnitude = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      magnitude += POW2(m_velocity[i]);
    }
    return sqrt(magnitude);
  }

  /// \brief Calculate the particle reynoldsnumber
  inline MFloat particleRe(const MFloat velocity, const MFloat density, const MFloat dynamicViscosity) {
    return density * velocity * m_diameter / dynamicViscosity;
  }

  /// \brief Calculate the reciprocal of the particle relaxation time
  inline MFloat fParticleRelTime(const MFloat dynamicViscosity) {
    return 18.0 * dynamicViscosity / (m_diameter * m_diameter * s_backPtr->m_material->density(m_temperature));
  }

  /// \brief Calculate the particle Weber-number (radius based)
  inline MFloat WeberNumber(const MFloat density, const MFloat velocity, const MFloat temperature) {
    return density * velocity * 0.5 * m_diameter / s_backPtr->m_material->spraySurfaceTension(temperature);
  }

  /// \brief Returns the relevant particle radius for the wall collision
  inline MFloat effectiveWallCollisionRadius() const override { return F1B2 * m_diameter; }

  void initVelocityAndDensity() {
    std::array<MFloat, nDim + 1> result{};
    std::vector<MFloat> weights;
    this->template interpolateAndCalcWeights<0, nDim + 1>(m_cellId, m_position.data(), result.data(), weights);
    for(MInt i = 0; i < nDim; i++) {
      m_oldFluidVel[i] = result[i];
      m_fluidVel[i] = result[i];
    }
    m_oldFluidDensity = result[nDim];
    m_fluidDensity = result[nDim];
  }

 private:
  void interpolateVelocityAndDensity(const MInt cellId, const MFloat* const position, MFloat* const velocity,
                                     MFloat& density, std::vector<MFloat>& weights) {
    std::array<MFloat, nDim + 1> result{};
    this->template interpolateAndCalcWeights<0, nDim + 1>(cellId, position, result.data(), weights);
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = result[n];
    }
    density = result[nDim];
  }

  void interpolateAllFluidVariables(const MInt cellId, const MFloat* const position, MFloat* const velocity,
                                    MFloat& density, MFloat& pressure, MFloat& temperature, MFloat& species,
                                    std::vector<MFloat>& weights) {
    std::array<MFloat, nDim + 4> result{};
    this->template interpolateAndCalcWeights<0, nDim + 4>(cellId, position, result.data(), weights);
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = result[n];
    }
    density = result[nDim];
    pressure = result[nDim + 1];
    temperature = result[nDim + 2];
    species = result[nDim + 3];
  }
  void interpolateAllFluidVariables(const MInt cellId, const MFloat* const position, MFloat* const velocity,
                                    MFloat& density, MFloat& pressure, MFloat& temperature,
                                    std::vector<MFloat>& weights) {
    std::array<MFloat, nDim + 3> result{};
    this->template interpolateAndCalcWeights<0, nDim + 3>(cellId, position, result.data(), weights);
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = result[n];
    }
    density = result[nDim];
    pressure = result[nDim + 1];
    temperature = result[nDim + 2];
  }

  void momentumCoupling();
  void heatCoupling();
  void massCoupling();
};

template <MInt nDim>
std::ostream& operator<<(std::ostream& os, const LPTSpherical<nDim>& m) {
  return os << m.m_partId;
}

#endif
