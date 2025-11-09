// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _MAIAPARTICLEELLIPSOIDAL_H
#define _MAIAPARTICLEELLIPSOIDAL_H

#include <cmath>
#include "UTIL/materialstate.h"
#include "lptbase.h"

template <MInt nDim, class SysEqn>
class MAIAFvCartesianSolverXD;

template <MInt nDim>
class LPTBase;

template <MInt nDim>
class LPTEllipsoidal : public LPTBase<nDim> {
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
  using LPTBase<nDim>::particleWallCollision;
  using LPTBase<nDim>::wallParticleCollision;

  using LPTBase<nDim>::interpolateAndCalcWeights;
  using LPTBase<nDim>::getNghbrList;
  using LPTBase<nDim>::checkCellChange;
  using LPTBase<nDim>::updateProperties;
  using LPTBase<nDim>::m_neighborList;

  LPTEllipsoidal();
  ~LPTEllipsoidal() override = default;

  LPTEllipsoidal(const LPTEllipsoidal&) = default;
  LPTEllipsoidal& operator=(const LPTEllipsoidal&) = default;

  /// particle angular velocity
  std::array<MFloat, nDim> m_angularVel{};
  /// particle angular acceleration
  std::array<MFloat, nDim> m_angularAccel{};
  /// fluid velocity
  std::array<MFloat, nDim> m_fluidVel{};
  /// fluid density
  MFloat m_fluidDensity = std::numeric_limits<MFloat>::quiet_NaN();

  /// particle acceleration of the last time step
  std::array<MFloat, nDim> m_oldAccel{};
  /// particle angular velocity of the last time step
  std::array<MFloat, nDim> m_oldAngularVel{};
  /// particle angular acceleration of the last time step
  std::array<MFloat, nDim> m_oldAngularAccel{};
  /// fluid velocity of the last time step
  std::array<MFloat, nDim> m_oldFluidVel{};

  /// gradient of fluid velocity in the particle fixed coordinate system
  std::array<MFloat, nDim * nDim> m_velocityGradientFluid{};

  /// old fluid density
  MFloat m_oldFluidDensity = std::numeric_limits<MFloat>::quiet_NaN();

  /// aspect ratio a/c of ellipsoidal: 0 < oblate < 1, =1 spherical, >1 prolate
  MFloat m_aspectRatio = std::numeric_limits<MFloat>::quiet_NaN();
  /// semi minor axis a
  MFloat m_semiMinorAxis = std::numeric_limits<MFloat>::quiet_NaN();
  /// Shape parameters for ellipsoid
  std::array<MFloat, 4> m_shapeParams{std::numeric_limits<MFloat>::quiet_NaN()};

  // Particle orientation using quaternions
  std::array<MFloat, 4> m_quaternion{};
  std::array<MFloat, 4> m_oldQuaternion{};

  /// temperature
  MFloat m_temperature = std::numeric_limits<MFloat>::quiet_NaN();
  /// fluid velocity magnitude
  MFloat m_fluidVelMag = std::numeric_limits<MFloat>::quiet_NaN();
  /// heat flux to cell
  MFloat m_heatFlux = std::numeric_limits<MFloat>::quiet_NaN();

  static constexpr MInt s_floatElements = 7 + 8 * nDim + 4 * 2;
  static constexpr MInt s_intElements = 0;

  MInt m_particleMajorAxis = 2; // index to the major axis, i.e. z-axis

  static MFloat s_lengthFactor;
  static MFloat s_Re;
  static MFloat s_Pr;
  static MFloat s_Sc;
  static MFloat s_We;
  static std::array<MFloat, nDim> s_Frm;

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

  void calculateMajorAxisOrientation(MFloat* majorAxis);

  MFloat dragFactor(const MFloat partRe, const MFloat beta = 1, const MFloat inclinationAngle = 0, const MFloat K0 = 0,
                    const MFloat K2 = 0);
  MFloat liftFactor(const MFloat partRe, const MFloat beta, const MFloat inclinationAngle, const MFloat K0,
                    const MFloat K2);
  MFloat torqueFactor(const MFloat partRe, const MFloat beta, const MFloat inclinationAngle);
  void setParticleFluidVelocities(MInt cellId, MFloat* position, MFloat* quaternions, MFloat* fluidVel,
                                  MFloat* fluidVelGrad);


  template <MInt>
  friend std::ostream& operator<<(std::ostream& os, const LPTEllipsoidal& m);

  /// \brief Returns the largest radius
  inline MFloat maxParticleRadius() const { return m_semiMinorAxis * m_aspectRatio; }
  /// \brief Returns the smallest radius
  inline MFloat minParticleRadius() const { return m_semiMinorAxis; }
  /// \brief Returns the radius in direction d
  inline MFloat particleRadius(MInt direction) const {
    if(direction == m_particleMajorAxis) {
      return maxParticleRadius();
    }
    return minParticleRadius();
  }

  /// \brief Returns the equivalent radius
  inline MFloat equivalentRadius() const { return F1B2 * m_eqDiameter; }

  /// \brief Returns the equivalent diameter
  inline MFloat equivalentDiameter() const { return m_eqDiameter; }

  /// \brief Calculates the current volume
  inline MFloat particleVolume() const { return F4B3 * PI * POW3(m_semiMinorAxis) * m_aspectRatio; }

  /// \brief Calculates the current mass
  inline MFloat particleMass() const {
    return F4B3 * PI * POW3(m_semiMinorAxis) * m_aspectRatio * s_backPtr->m_material->density();
  }

  /// \brief Calculates the magnitude of the relative velocity of the two given velocity vectors
  inline MFloat magRelVel(const MFloat* const velocity1, const MFloat* const velocity2) const {
    MFloat magnitude = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      magnitude += POW2(velocity1[i] - velocity2[i]);
    }
    return sqrt(magnitude);
  }

  /// \brief Calculates the particle Reynolds number
  inline MFloat particleRe(const MFloat velocity, const MFloat density, const MFloat dynamicViscosity) {
    return density * velocity * m_eqDiameter / dynamicViscosity;
  }

  /// \brief Calculates the reciprocal of the particle relaxation time
  inline MFloat fParticleRelTime(const MFloat dynamicViscosity) {
    return 18.0 * dynamicViscosity / (m_eqDiameter * m_eqDiameter * s_backPtr->m_material->density());
  }

  /// \brief Returns the relevant particle radius for the wall collision
  inline MFloat effectiveWallCollisionRadius() const override { return equivalentRadius(); }

  inline void initEllipsoialProperties() {
    ASSERT(!std::isnan(m_aspectRatio), "Aspect Ratio not set!");
    ASSERT(!std::isnan(m_semiMinorAxis), "Semi Minor Axis not set!");
    initEqDiameter();
    initShapeParams();
  }

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
  /// equivalent particle diameter
  MFloat m_eqDiameter = std::numeric_limits<MFloat>::quiet_NaN();

  void initShapeParams();
  void initEqDiameter() { m_eqDiameter = 2 * m_semiMinorAxis * pow(m_aspectRatio, F1B3); }

  void interpolateVelocityAndDensity(const MInt cellId, const MFloat* const position, MFloat* const velocity,
                                     MFloat& density, std::vector<MFloat>& weights) {
    std::array<MFloat, nDim + 1> result{};
    this->template interpolateAndCalcWeights<0, nDim + 1>(cellId, position, result.data(), weights);
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = result[n];
    }
    density = result[nDim];
  }

  void interpolateVelocityAndVelocitySlopes(const MInt cellId, const MFloat* const position, MFloat* const velocity,
                                            MFloat* const velocityGradient, std::vector<MFloat>& weights) {
    std::array<MFloat, nDim> result{};
    std::array<MFloat, nDim * nDim> gradientResult{};
    this->template interpolateAndCalcWeights<0, nDim>(cellId, position, result.data(), weights, gradientResult.data());
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = result[n];
    }
    for(MInt n = 0; n < nDim * nDim; n++) {
      velocityGradient[n] = gradientResult[n];
    }
  }

  void momentumCoupling();
  void heatCoupling();
  void massCoupling();
};

template <MInt nDim>
std::ostream& operator<<(std::ostream& os, const LPTEllipsoidal<nDim>& m) {
  return os << m.m_partId;
}

#endif
