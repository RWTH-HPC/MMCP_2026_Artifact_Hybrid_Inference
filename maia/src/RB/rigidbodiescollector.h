// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RBCOLLECTOR_H_
#define RBCOLLECTOR_H_

#include <algorithm>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computation. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define RBCOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'RBCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks
// for all accessors. It is enabled for build type "extra_debug".
//#ifdef MAIA_EXTRA_DEBUG
// NOTE: enable checks in normal debug mode until everything is working correctly!
#ifndef NDEBUG
#define RBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(RBCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID(id)                                                                                            \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_DIM(id)                                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < nDim,                                                                        \
                          "dim id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(nDim) + ")", AT_);   \
  } while(false)
#define ENSURE_VALID_ROT(id)                                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < nRot,                                                                        \
                          "rot id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(nRot) + ")", AT_);   \
  } while(false)
#define ENSURE_VALID_QUAT(id)                                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < nQuat,                                                                       \
                          "quat id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(nQuat) + ")", AT_); \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(condition, message, AT_);                                                                    \
  } while(false)
#else
#define ENSURE_VALID_ID(id)                                                                                            \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIM(id)                                                                                           \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_ROT(id)                                                                                           \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_QUAT(id)                                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace rb {
namespace collector {

// Type traits for invalid values. These values are used to initialize/erase nodes
template <class T>
struct Invalid {};

// Invalid value for ids is 'INT_MIN'
template <>
struct Invalid<MInt> {
  static constexpr MInt value() { return std::numeric_limits<MInt>::min(); }
};

// Invalid value for longs is 'INT_MIN'
template <>
struct Invalid<MLong> {
  static constexpr MLong value() { return std::numeric_limits<MLong>::min(); }
};

// Invalid value for floats is 'NaN'
template <>
struct Invalid<MFloat> {
  static constexpr MFloat value() {
#ifdef MAIA_PGI_COMPILER
    return std::numeric_limits<MFloat>::quiet_NaN();
#else
    return std::numeric_limits<MFloat>::signaling_NaN();
#endif
  }
};

// Body status. Used for paralellization
enum class Status : MInt { local, remote, invalid };
// Default (or invalidated) value
template <>
struct Invalid<Status> {
  static constexpr Status value() { return Status::invalid; }
};


/// Class that represents DG element collector.
template <MInt nDim>
class RigidBodyCollector : public maia::container::Container<RigidBodyCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<RigidBodyCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<RigidBodyCollector<nDim>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  // Types
  template <class T>
  using Invalid = typename maia::rb::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr RigidBodyCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;


  // Accessors
  MFloat& bodyCenter(const MInt id, const MInt dim);
  MFloat bodyCenter(const MInt id, const MInt dim) const;
  MFloat& bodyVelocity(const MInt id, const MInt dim);
  MFloat bodyVelocity(const MInt id, const MInt dim) const;
  MFloat& bodyAcceleration(const MInt id, const MInt dim);
  MFloat bodyAcceleration(const MInt id, const MInt dim) const;
  MFloat& bodyTemperature(const MInt id);
  MFloat bodyTemperature(const MInt id) const;
  MFloat& bodyCenterOld(const MInt id, const MInt dim);
  MFloat bodyCenterOld(const MInt id, const MInt dim) const;
  MFloat& bodyVelocityOld(const MInt id, const MInt dim);
  MFloat bodyVelocityOld(const MInt id, const MInt dim) const;
  MFloat& bodyAccelerationOld(const MInt id, const MInt dim);
  MFloat bodyAccelerationOld(const MInt id, const MInt dim) const;
  MFloat& bodyTemperatureOld(const MInt id);
  MFloat bodyTemperatureOld(const MInt id) const;
  MFloat& bodyQuaternionT1B2(const MInt id, const MInt dim);
  MFloat bodyQuaternionT1B2(const MInt id, const MInt dim) const;
  MFloat& bodyQuaternionT1(const MInt id, const MInt dim);
  MFloat bodyQuaternionT1(const MInt id, const MInt dim) const;
  MFloat& angularVelocityT1(const MInt id, const MInt dim);
  MFloat angularVelocityT1(const MInt id, const MInt dim) const;
  MFloat& angularVelocityBodyT1(const MInt id, const MInt dim);
  MFloat angularVelocityBodyT1(const MInt id, const MInt dim) const;
  MFloat& angularVelocityT1B2(const MInt id, const MInt dim);
  MFloat angularVelocityT1B2(const MInt id, const MInt dim) const;
  MFloat& angularVelocityBodyT1B2(const MInt id, const MInt dim);
  MFloat angularVelocityBodyT1B2(const MInt id, const MInt dim) const;
  MFloat& angularAccelerationT1(const MInt id, const MInt dim);
  MFloat angularAccelerationT1(const MInt id, const MInt dim) const;
  MFloat& angularAccelerationBody(const MInt id, const MInt dim);
  MFloat angularAccelerationBody(const MInt id, const MInt dim) const;
  MFloat& torqueT1(const MInt id, const MInt dim);
  MFloat torqueT1(const MInt id, const MInt dim) const;
  MFloat& bodyForce(const MInt id, const MInt dim);
  MFloat bodyForce(const MInt id, const MInt dim) const;
  MFloat& bodyHeatFlux(const MInt id);
  MFloat bodyHeatFlux(const MInt id) const;
  MFloat& bodyDensityRatio(const MInt id);
  MFloat bodyDensityRatio(const MInt id) const;
  MFloat& bodyInertia(const MInt id, const MInt dim);
  MFloat bodyInertia(const MInt id, const MInt dim) const;
  MFloat& bodyRadius(const MInt id);
  MFloat bodyRadius(const MInt id) const;
  MFloat& bodyRadii(const MInt id, const MInt dim);
  MFloat bodyRadii(const MInt id, const MInt dim) const;

  Status& status(const MInt id);
  Status status(const MInt id) const;

  // Auxiliary methods
  void advanceBodies();

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

 private:
  // Member variables
  static constexpr MInt nRot = (nDim == 3) ? 3 : 1;
  static constexpr MInt nQuat = (nDim == 3) ? 4 : 0;

  // Data containers

  // Body data current timestep
  Storage<MFloat> m_bodyCenter{};
  Storage<MFloat> m_bodyVelocity{};
  Storage<MFloat> m_bodyAcceleration{};

  Storage<MFloat> m_bodyTemperature{};

  // Body data old timestep
  Storage<MFloat> m_bodyCenterOld{};
  Storage<MFloat> m_bodyVelocityOld{};
  Storage<MFloat> m_bodyAccelerationOld{};

  Storage<MFloat> m_bodyTemperatureOld{};

  // Rotational data
  Storage<MFloat> m_bodyQuaternionT1B2{};
  Storage<MFloat> m_bodyQuaternionT1{};

  Storage<MFloat> m_angularVelocityT1{};
  Storage<MFloat> m_angularVelocityBodyT1{};
  Storage<MFloat> m_angularVelocityT1B2{};
  Storage<MFloat> m_angularVelocityBodyT1B2{};

  Storage<MFloat> m_angularAccelerationT1{};
  Storage<MFloat> m_angularAccelerationBody{};

  Storage<MFloat> m_torqueT1{};

  Storage<MFloat> m_bodyForce{};
  Storage<MFloat> m_bodyHeatFlux{};

  Storage<MFloat> m_bodyDensityRatio{};
  Storage<MFloat> m_bodyInertia{};

  Storage<MFloat> m_bodyRadius{};
  Storage<MFloat> m_bodyRadii{};

  Storage<Status> m_status{};
};

/**
 * \brief Copy the body data for time t to t-1 to prepare the next timestep
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodyCollector<nDim>::advanceBodies() {
  m_bodyCenterOld = m_bodyCenter;
  m_bodyVelocityOld = m_bodyVelocity;
  m_bodyAccelerationOld = m_bodyAcceleration;

  m_bodyTemperatureOld = m_bodyTemperature;
}


/// Reset, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void RigidBodyCollector<nDim>::reset() {
  resetStorage(nDim, m_bodyCenter);
  resetStorage(nDim, m_bodyVelocity);
  resetStorage(nDim, m_bodyAcceleration);
  resetStorage(1, m_bodyTemperature);
  resetStorage(nDim, m_bodyCenterOld);
  resetStorage(nDim, m_bodyVelocityOld);
  resetStorage(nDim, m_bodyAccelerationOld);
  resetStorage(1, m_bodyTemperatureOld);
  resetStorage(nQuat, m_bodyQuaternionT1B2);
  resetStorage(nQuat, m_bodyQuaternionT1);
  resetStorage(nRot, m_angularVelocityT1);
  resetStorage(nRot, m_angularVelocityBodyT1);
  resetStorage(nRot, m_angularVelocityT1B2);
  resetStorage(nRot, m_angularVelocityBodyT1B2);
  resetStorage(nRot, m_angularAccelerationT1);
  resetStorage(nRot, m_angularAccelerationBody);
  resetStorage(nRot, m_torqueT1);
  resetStorage(nDim, m_bodyForce);
  resetStorage(1, m_bodyHeatFlux);
  resetStorage(1, m_bodyDensityRatio);
  resetStorage(nRot, m_bodyInertia);
  resetStorage(1, m_bodyRadius);
  resetStorage(nDim, m_bodyRadii);
  resetStorage(1, m_status);
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void RigidBodyCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  fill_invalid(m_bodyCenter, begin, end, nDim);
  fill_invalid(m_bodyVelocity, begin, end, nDim);
  fill_invalid(m_bodyAcceleration, begin, end, nDim);
  fill_invalid(m_bodyTemperature, begin, end, 1);
  fill_invalid(m_bodyCenterOld, begin, end, nDim);
  fill_invalid(m_bodyVelocityOld, begin, end, nDim);
  fill_invalid(m_bodyAccelerationOld, begin, end, nDim);
  fill_invalid(m_bodyTemperatureOld, begin, end, 1);
  fill_invalid(m_bodyQuaternionT1B2, begin, end, nQuat);
  fill_invalid(m_bodyQuaternionT1, begin, end, nQuat);
  fill_invalid(m_angularVelocityT1, begin, end, nRot);
  fill_invalid(m_angularVelocityBodyT1, begin, end, nRot);
  fill_invalid(m_angularVelocityT1B2, begin, end, nRot);
  fill_invalid(m_angularVelocityBodyT1B2, begin, end, nRot);
  fill_invalid(m_angularAccelerationT1, begin, end, nRot);
  fill_invalid(m_angularAccelerationBody, begin, end, nRot);
  fill_invalid(m_torqueT1, begin, end, nRot);
  fill_invalid(m_bodyForce, begin, end, nDim);
  fill_invalid(m_bodyHeatFlux, begin, end, 1);
  fill_invalid(m_bodyDensityRatio, begin, end, 1);
  fill_invalid(m_bodyInertia, begin, end, nRot);
  fill_invalid(m_bodyRadius, begin, end, 1);
  fill_invalid(m_bodyRadii, begin, end, nDim);
  fill_invalid(m_status, begin, end, 1);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void RigidBodyCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                              const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  copyData(source.m_bodyCenter, m_bodyCenter, c, begin, end, destination, nDim);
  copyData(source.m_bodyVelocity, m_bodyVelocity, c, begin, end, destination, nDim);
  copyData(source.m_bodyAcceleration, m_bodyAcceleration, c, begin, end, destination, nDim);
  copyData(source.m_bodyTemperature, m_bodyTemperature, c, begin, end, destination, 1);
  copyData(source.m_bodyCenterOld, m_bodyCenterOld, c, begin, end, destination, nDim);
  copyData(source.m_bodyVelocityOld, m_bodyVelocityOld, c, begin, end, destination, nDim);
  copyData(source.m_bodyAccelerationOld, m_bodyAccelerationOld, c, begin, end, destination, nDim);
  copyData(source.m_bodyTemperatureOld, m_bodyTemperatureOld, c, begin, end, destination, 1);
  copyData(source.m_bodyQuaternionT1B2, m_bodyQuaternionT1B2, c, begin, end, destination, nQuat);
  copyData(source.m_bodyQuaternionT1, m_bodyQuaternionT1, c, begin, end, destination, nQuat);
  copyData(source.m_angularVelocityT1, m_angularVelocityT1, c, begin, end, destination, nRot);
  copyData(source.m_angularVelocityBodyT1, m_angularVelocityBodyT1, c, begin, end, destination, nRot);
  copyData(source.m_angularVelocityT1B2, m_angularVelocityT1B2, c, begin, end, destination, nRot);
  copyData(source.m_angularVelocityBodyT1B2, m_angularVelocityBodyT1B2, c, begin, end, destination, nRot);
  copyData(source.m_angularAccelerationT1, m_angularAccelerationT1, c, begin, end, destination, nRot);
  copyData(source.m_angularAccelerationBody, m_angularAccelerationBody, c, begin, end, destination, nRot);
  copyData(source.m_torqueT1, m_torqueT1, c, begin, end, destination, nRot);
  copyData(source.m_bodyForce, m_bodyForce, c, begin, end, destination, nDim);
  copyData(source.m_bodyHeatFlux, m_bodyHeatFlux, c, begin, end, destination, 1);
  copyData(source.m_bodyDensityRatio, m_bodyDensityRatio, c, begin, end, destination, 1);
  copyData(source.m_bodyInertia, m_bodyInertia, c, begin, end, destination, nRot);
  copyData(source.m_bodyRadius, m_bodyRadius, c, begin, end, destination, 1);
  copyData(source.m_bodyRadii, m_bodyRadii, c, begin, end, destination, nDim);
  copyData(source.m_status, m_status, c, begin, end, destination, 1);
}

/// Accessor for the current body center
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyCenter(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyCenter[id * nDim + dim];
}
/// Accessor for the current body center (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyCenter(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyCenter[id * nDim + dim];
}

/// Accessor for the old body center
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyCenterOld(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyCenterOld[id * nDim + dim];
}
/// Accessor for the old body center (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyCenterOld(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyCenterOld[id * nDim + dim];
}

/// Accessor for the current body velocity
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyVelocity(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyVelocity[id * nDim + dim];
}
/// Accessor for the current body velocity (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyVelocity(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyVelocity[id * nDim + dim];
}

/// Accessor for the old body velocity
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyVelocityOld(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyVelocityOld[id * nDim + dim];
}
/// Accessor for the old body velocity (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyVelocityOld(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyVelocityOld[id * nDim + dim];
}

/// Accessor for the current body acceleration
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyAcceleration(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyAcceleration[id * nDim + dim];
}
/// Accessor for the current body acceleration (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyAcceleration(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyAcceleration[id * nDim + dim];
}

/// Accessor for the old body acceleration
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyAccelerationOld(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyAccelerationOld[id * nDim + dim];
}
/// Accessor for the old body acceleration (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyAccelerationOld(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyAccelerationOld[id * nDim + dim];
}

/// Accessor for the body force
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyForce(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyForce[id * nDim + dim];
}
/// Accessor for the body force (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyForce(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyForce[id * nDim + dim];
}

/// Accessor for the body temperature
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyTemperature(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyTemperature[id];
}
/// Accessor for the body temperature (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyTemperature(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyTemperature[id];
}

/// Accessor for the old body temperature
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyTemperatureOld(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyTemperatureOld[id];
}
/// Accessor for the old body temperature (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyTemperatureOld(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyTemperatureOld[id];
}

/// Accessor for the body heat flux
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyHeatFlux(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyHeatFlux[id];
}
/// Accessor for the body heat flux (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyHeatFlux(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyHeatFlux[id];
}

/// Accessor for the body density ratio
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyDensityRatio(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyDensityRatio[id];
}
/// Accessor for the body density ratio (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyDensityRatio(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyDensityRatio[id];
}

/// Accessor for the body inertia
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyInertia(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_bodyInertia[id * nRot + dim];
}
/// Accessor for the body inertia (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyInertia(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_bodyInertia[id * nRot + dim];
}

/// Accessor for the body radius
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyRadius(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyRadius[id];
}
/// Accessor for the body radius (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyRadius(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_bodyRadius[id];
}

/// Accessor for the body radii
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyRadii(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyRadii[id * nDim + dim];
}
/// Accessor for the body radii (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyRadii(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_DIM(dim);
  return m_bodyRadii[id * nDim + dim];
}

/// Accessor for the body bodyQuaternion T1
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyQuaternionT1(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_QUAT(dim);
  return m_bodyQuaternionT1[id * nQuat + dim];
}
/// Accessor for the body bodyQuaternion T1 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyQuaternionT1(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_QUAT(dim);
  return m_bodyQuaternionT1[id * nQuat + dim];
}

/// Accessor for the body bodyQuaternion T1B2
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::bodyQuaternionT1B2(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_QUAT(dim);
  return m_bodyQuaternionT1B2[id * nQuat + dim];
}
/// Accessor for the body bodyQuaternion T1B2 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::bodyQuaternionT1B2(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_QUAT(dim);
  return m_bodyQuaternionT1B2[id * nQuat + dim];
}

/// Accessor for the angular velocity T1
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularVelocityT1(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityT1[id * nRot + dim];
}
/// Accessor for the angular velocity T1 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularVelocityT1(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityT1[id * nRot + dim];
}

/// Accessor for the angular velocity T1B2
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularVelocityT1B2(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityT1B2[id * nRot + dim];
}
/// Accessor for the angular velocity T1B2 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularVelocityT1B2(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityT1B2[id * nRot + dim];
}

/// Accessor for the angular velocity of the body T1
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularVelocityBodyT1(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityBodyT1[id * nRot + dim];
}
/// Accessor for the angular velocity of the body T1 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularVelocityBodyT1(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityBodyT1[id * nRot + dim];
}

/// Accessor for the angular velocity of the body T1B2
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularVelocityBodyT1B2(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityBodyT1B2[id * nRot + dim];
}
/// Accessor for the angular velocity of the body T1B2 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularVelocityBodyT1B2(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularVelocityBodyT1B2[id * nRot + dim];
}

/// Accessor for the angular acceleration
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularAccelerationBody(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularAccelerationBody[id * nRot + dim];
}
/// Accessor for the angular acceleration (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularAccelerationBody(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularAccelerationBody[id * nRot + dim];
}

/// Accessor for the angular acceleration T1
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::angularAccelerationT1(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularAccelerationT1[id * nRot + dim];
}
/// Accessor for the angular acceleration T1 (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::angularAccelerationT1(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_angularAccelerationT1[id * nRot + dim];
}

/// Accessor for the torque
template <MInt nDim>
MFloat& RigidBodyCollector<nDim>::torqueT1(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_torqueT1[id * nRot + dim];
}
/// Accessor for the torque (const)
template <MInt nDim>
MFloat RigidBodyCollector<nDim>::torqueT1(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  ENSURE_VALID_ROT(dim);
  return m_torqueT1[id * nRot + dim];
}

/// Accessor for the body status
template <MInt nDim>
Status& RigidBodyCollector<nDim>::status(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_status[id];
}
/// Accessor for the body status (const)
template <MInt nDim>
Status RigidBodyCollector<nDim>::status(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef RBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID(id);
  return m_status[id];
}

} // namespace collector
} // namespace rb
} // namespace maia

// Undefine macros that should not be used outside this file
#undef RBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID
#undef ENSURE_VALID_DIM
#undef ENSURE_VALID_ROT
#undef ENSURE_VALID_QUAT
#undef ENSURE_CONDITION

#endif // ifndef DGELEMENTCOLLECTOR_H_
