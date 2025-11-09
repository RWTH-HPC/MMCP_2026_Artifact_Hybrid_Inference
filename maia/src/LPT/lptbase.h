// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _ParticleBase_H
#define _ParticleBase_H

#include <array>
#include <vector>
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"
#include "lptbaseproperties.h"

template <MInt nDim>
class LPT;

template <MInt nDim>
class LPTBase {
 public:
  virtual ~LPTBase() = default;

  using BitsetType = maia::lpt::baseProperty::BitsetType;

  static MInt s_interpolationOrder;
  static MInt s_interpolationMethod;
  static MFloat s_distFactorImp;
  static LPT<nDim>* s_backPtr;

  static constexpr MInt s_floatElements = 3 * nDim + 1;
  static constexpr MInt s_intElements = 5;

  std::array<MFloat, nDim> m_position{};
  std::array<MFloat, nDim> m_velocity{};
  std::array<MFloat, nDim> m_accel{};

  /// particle position of the last time step
  std::array<MFloat, nDim> m_oldPos{};
  /// particle velocity of the last time step
  std::array<MFloat, nDim> m_oldVel{};

  // initialised with invalid values
  MFloat m_densityRatio = -2;
  MLong m_partId = -2;
  MInt m_cellId = -2;
  MInt m_oldCellId = -2;

  /// creation time modifier
  MFloat m_creationTime = std::numeric_limits<MFloat>::quiet_NaN();

  BitsetType m_properties;
  // unversal acces
  BitsetType::reference hasProperty(const LptBaseProperty p);
  MBool hasProperty(const LptBaseProperty p) const;

  // access directly to individual base-properties!
  BitsetType::reference isWindow();
  MBool isWindow() const;
  BitsetType::reference reqSend();
  MBool reqSend() const;
  BitsetType::reference reqBroadcast();
  MBool reqBroadcast() const;
  BitsetType::reference isInvalid();
  MBool isInvalid() const;
  BitsetType::reference hasCollided();
  MBool hasCollided() const;
  BitsetType::reference hadWallColl();
  MBool hadWallColl() const;
  BitsetType::reference firstStep();
  MBool firstStep() const;
  BitsetType::reference toBeDeleted();
  MBool toBeDeleted() const;
  BitsetType::reference wasSend();
  MBool wasSend() const;
  BitsetType::reference toBeRespawn();
  MBool toBeRespawn() const;
  BitsetType::reference fullyEvaporated();
  MBool fullyEvaporated() const;

  // store neighbor cells including diagonals
  std::vector<MInt> m_neighborList;

  void getNghbrList(std::vector<MInt>&, const MInt);

  template <MInt a, MInt b>
  void interpolateAndCalcWeights(const MInt cellId, const MFloat* const x, MFloat* const result,
                                 std::vector<MFloat>& weight, MFloat* const gradientResult = nullptr);

  virtual void energyEquation() = 0;
  virtual void coupling() = 0;
  virtual void motionEquation() = 0;
  virtual void advanceParticle() = 0;
  virtual void resetWeights() = 0;


  virtual void particleWallCollision();
  virtual void wallParticleCollision();
  void initProperties();

  void checkCellChange(const MFloat* oldPosition, const MBool allowHaloNonLeaf = false);
  void updateProperties(const MBool init = true);

  virtual MFloat effectiveWallCollisionRadius() const { TERMM(-1, "Effective wall collision radius not implemented!"); }
};

/// Accessor for properties.
template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::hasProperty(const LptBaseProperty p) {
  return m_properties[maia::lpt::baseProperty::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool LPTBase<nDim>::hasProperty(const LptBaseProperty p) const {
  return m_properties[maia::lpt::baseProperty::p(p)];
}

template <MInt nDim>
MBool LPTBase<nDim>::isWindow() const {
  return hasProperty(LptBaseProperty::IsWindow);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::isWindow() {
  return hasProperty(LptBaseProperty::IsWindow);
}
template <MInt nDim>
MBool LPTBase<nDim>::reqSend() const {
  return hasProperty(LptBaseProperty::ReqSend);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::reqSend() {
  return hasProperty(LptBaseProperty::ReqSend);
}
template <MInt nDim>
MBool LPTBase<nDim>::reqBroadcast() const {
  return hasProperty(LptBaseProperty::ReqBroadcast);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::reqBroadcast() {
  return hasProperty(LptBaseProperty::ReqBroadcast);
}
template <MInt nDim>
MBool LPTBase<nDim>::isInvalid() const {
  return hasProperty(LptBaseProperty::IsInvalid);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::isInvalid() {
  return hasProperty(LptBaseProperty::IsInvalid);
}
template <MInt nDim>
MBool LPTBase<nDim>::hasCollided() const {
  return hasProperty(LptBaseProperty::HasCollided);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::hasCollided() {
  return hasProperty(LptBaseProperty::HasCollided);
}
template <MInt nDim>
MBool LPTBase<nDim>::hadWallColl() const {
  return hasProperty(LptBaseProperty::HadWallColl);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::hadWallColl() {
  return hasProperty(LptBaseProperty::HadWallColl);
}

template <MInt nDim>
MBool LPTBase<nDim>::firstStep() const {
  return hasProperty(LptBaseProperty::FirstStep);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::firstStep() {
  return hasProperty(LptBaseProperty::FirstStep);
}

template <MInt nDim>
MBool LPTBase<nDim>::toBeDeleted() const {
  return hasProperty(LptBaseProperty::ToBeDeleted);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::toBeDeleted() {
  return hasProperty(LptBaseProperty::ToBeDeleted);
}

template <MInt nDim>
MBool LPTBase<nDim>::wasSend() const {
  return hasProperty(LptBaseProperty::WasSend);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::wasSend() {
  return hasProperty(LptBaseProperty::WasSend);
}
template <MInt nDim>
MBool LPTBase<nDim>::toBeRespawn() const {
  return hasProperty(LptBaseProperty::ToBeRespawn);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::toBeRespawn() {
  return hasProperty(LptBaseProperty::ToBeRespawn);
}

template <MInt nDim>
MBool LPTBase<nDim>::fullyEvaporated() const {
  return hasProperty(LptBaseProperty::FullyEvaporated);
}

template <MInt nDim>
LPTBase<nDim>::BitsetType::reference LPTBase<nDim>::fullyEvaporated() {
  return hasProperty(LptBaseProperty::FullyEvaporated);
}

template <MInt nDim>
void LPTBase<nDim>::initProperties() {
  isWindow() = false;
  reqSend() = false;
  reqBroadcast() = false;
  isInvalid() = false;
  hasCollided() = false;
  firstStep() = false;
  toBeDeleted() = false;
  wasSend() = false;
  toBeRespawn() = false;
  fullyEvaporated() = false;
  hadWallColl() = false;
}

#endif
