// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LPTCOLLECTOR_H_
#define LPTCOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
#include "lptcellproperties.h"


// The macro 'LPTCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity checks
// for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define LPTCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

#ifdef LPTCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

// Sanity-checking macros for accessors
#if defined(MAIA_ASSERT_ACCESSORS) || defined(LPTCOLLECTOR_SANITY_CHECKS_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_SET_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        (id) >= 0 && (id) < maxNoSets(),                                                                               \
        "set id = " + std::to_string(id) + "is out-of-bounds [0, " + std::to_string(maxNoSets()) + ")", AT_);          \
  } while(false)
#define ENSURE_VALID_DIM_ID_ACCESSOR(dim)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(dim >= 0 && dim < nDim,                                                                      \
                          "dim = " + std::to_string(dim) + " is out-of-bounds [0, " + std::to_string(nDim) + ")",      \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != LptCell::NumProperties, "Invalid property", AT_);                                       \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SET_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIM_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dim)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


// Namespace for auxiliary functions/classes
namespace maia {
namespace lpt {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::lpt::cell::BitsetType;

// Type traits for invalid values. These values are used to initialize/erase nodes
template <class T>
struct Invalid {};

// Invalid value for ids is 'INT_MIN'
template <>
struct Invalid<MInt> {
  static constexpr MInt value() { return std::numeric_limits<MInt>::min(); }
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

// Invalid value for BitsetProperties is '0'
template <>
struct Invalid<BitsetType> {
  static constexpr BitsetType value() { return {0}; }
};


/// Class that represents LPT cell collector.
template <MInt nDim>
class LptCells : public maia::container::Container<LptCells<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<LptCells<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<LptCells<nDim>, Invalid>;
  using Base::resetStorage;
  using self = LptCells<nDim>;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  using BitsetType = maia::lpt::cell::BitsetType;

  // Types
  template <class T>
  using Invalid = typename maia::lpt::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr LptCells() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Property-related accessors
  BitsetType::reference hasProperty(const MInt id, const LptCell p);
  MBool hasProperty(const MInt id, const LptCell p) const;
  void resetProperties(const MInt id);
  BitsetType& properties(const MInt id);

  MFloat volumeFraction(const MInt id) const;
  MFloat& volumeFraction(const MInt id);

  MInt noParticles(const MInt id) const;
  MInt& noParticles(const MInt id);

  MInt noEllipsoids(const MInt id) const;
  MInt& noEllipsoids(const MInt id);

  MInt bndryCellId(const MInt id) const;
  MInt& bndryCellId(const MInt id);

  MFloat fluidVariable(const MInt id, const MInt var) const;
  MFloat& fluidVariable(const MInt id, const MInt var);

  MFloat fluidSpecies(const MInt id) const;
  MFloat& fluidSpecies(const MInt id);

  MFloat& massFlux(const MInt id);
  MFloat massFlux(const MInt id) const;

  MFloat& heatFlux(const MInt id);
  MFloat heatFlux(const MInt id) const;

  MFloat& momentumFlux(const MInt id, const MInt dim);
  MFloat momentumFlux(const MInt id, const MInt dim) const;

  MFloat& workFlux(const MInt id);
  MFloat workFlux(const MInt id) const;

  MFloat& velocitySlope(const MInt id, const MInt varId, const MInt dim);
  MFloat velocitySlope(const MInt id, const MInt varId, const MInt dim) const;

  void setLptCollectorNoSpecies(const MInt noSpecies) { m_noSpecies = noSpecies; }

  void setLptCollectorCoupling(const MBool mass, const MBool momentum, const MBool heat);

  void setLptCollectorSlopes() { m_storeSlopes = true; }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  MBool m_massCoupling = false;
  MBool m_momentumCoupling = false;
  MBool m_heatCoupling = false;
  static constexpr MInt m_noVariables = nDim + 3;
  MInt m_noSpecies = 0;
  MBool m_storeSlopes = false;

  constexpr MBool hasMassCoupling() const { return m_massCoupling; }
  constexpr MBool hasMomentumCoupling() const { return m_momentumCoupling; }
  constexpr MBool hasHeatCoupling() const { return m_heatCoupling; }
  constexpr MInt noVars() const { return m_noVariables; };
  constexpr MInt noSpecies() const { return m_noSpecies; };
  constexpr MBool hasSlopes() const { return m_storeSlopes; }

  Storage<BitsetType> m_properties;

  Storage<MInt> m_noParticles{};
  Storage<MInt> m_noEllipsoids{};
  Storage<MFloat> m_volumeFraction{};
  Storage<MInt> m_bndryCellId{};

  Storage<MFloat> m_variables{};

  Storage<MFloat> m_species{};

  Storage<MFloat> m_massFlux{};
  Storage<MFloat> m_heatFlux{};
  Storage<MFloat> m_momentumFlux{};
  Storage<MFloat> m_workFlux{};

  Storage<MFloat> m_velocitySlopes{};
};

/// Set the lpt-collector type
template <MInt nDim>
void LptCells<nDim>::setLptCollectorCoupling(const MBool mass, const MBool momentum, const MBool heat) {
  m_massCoupling = mass;
  m_momentumCoupling = momentum;
  m_heatCoupling = heat;
}
/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void LptCells<nDim>::reset() {
  // resetStorage(#of items per cell, Storage<> name);

  resetStorage(1, m_properties);
  resetStorage(1, m_noParticles);
  resetStorage(1, m_noEllipsoids);
  resetStorage(1, m_volumeFraction);
  resetStorage(1, m_bndryCellId);

  resetStorage(noVars(), m_variables);

  resetStorage(noSpecies(), m_species);

  if(hasMassCoupling()) {
    resetStorage(1, m_massFlux);
  }
  if(hasHeatCoupling()) resetStorage(1, m_heatFlux);
  if(hasMomentumCoupling()) {
    resetStorage(nDim, m_momentumFlux);
    resetStorage(1, m_workFlux);
  }
  if(hasSlopes()) {
    resetStorage(nDim * nDim, m_velocitySlopes);
  }
}

/// Accessor for properties.
template <MInt nDim>
LptCells<nDim>::BitsetType::reference LptCells<nDim>::hasProperty(const MInt id, const LptCell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties.at(id)[maia::lpt::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool LptCells<nDim>::hasProperty(const MInt id, const LptCell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties.at(id)[maia::lpt::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void LptCells<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties.at(id).reset();
}
/// Accessor for properties.
template <MInt nDim>
BitsetType& LptCells<nDim>::properties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}

/// Accessor for volumeFraction
template <MInt nDim>
MFloat& LptCells<nDim>::volumeFraction(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_volumeFraction[id];
}
template <MInt nDim>
MFloat LptCells<nDim>::volumeFraction(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_volumeFraction[id];
}
/// Accessor for noParticles
template <MInt nDim>
MInt& LptCells<nDim>::noParticles(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noParticles[id];
}
template <MInt nDim>
MInt LptCells<nDim>::noParticles(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noParticles[id];
}

/// Accessor for noEllipsoids
template <MInt nDim>
MInt& LptCells<nDim>::noEllipsoids(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noEllipsoids[id];
}
template <MInt nDim>
MInt LptCells<nDim>::noEllipsoids(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noEllipsoids[id];
}

/// Accessor for bndryCellId
template <MInt nDim>
MInt& LptCells<nDim>::bndryCellId(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCellId[id];
}
template <MInt nDim>
MInt LptCells<nDim>::bndryCellId(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCellId[id];
}

/// Accessor for fluidVariables
template <MInt nDim>
MFloat& LptCells<nDim>::fluidVariable(const MInt id, const MInt var) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_variables[id * noVars() + var];
}
template <MInt nDim>
MFloat LptCells<nDim>::fluidVariable(const MInt id, const MInt var) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_variables[id * noVars() + var];
}

/// Accessor for fluidVariables
template <MInt nDim>
MFloat& LptCells<nDim>::fluidSpecies(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_species[id];
}
template <MInt nDim>
MFloat LptCells<nDim>::fluidSpecies(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_species[id];
}

/// Accessor for massFlux
template <MInt nDim>
MFloat& LptCells<nDim>::massFlux(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_massFlux[id];
}
template <MInt nDim>
MFloat LptCells<nDim>::massFlux(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_massFlux[id];
}
/// Accessor for momentumFlux
template <MInt nDim>
MFloat& LptCells<nDim>::momentumFlux(const MInt id, const MInt dim) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_momentumFlux[id * nDim + dim];
}
template <MInt nDim>
MFloat LptCells<nDim>::momentumFlux(const MInt id, const MInt dim) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_momentumFlux[id * nDim + dim];
}
/// Accessor for workFlux
template <MInt nDim>
MFloat& LptCells<nDim>::workFlux(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_workFlux[id];
}
template <MInt nDim>
MFloat LptCells<nDim>::workFlux(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_workFlux[id];
}

/// Accessor for heatFlux
template <MInt nDim>
MFloat& LptCells<nDim>::heatFlux(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_heatFlux[id];
}
template <MInt nDim>
MFloat LptCells<nDim>::heatFlux(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_heatFlux[id];
}

template <MInt nDim>
MFloat& LptCells<nDim>::velocitySlope(const MInt id, const MInt varId, const MInt dir) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dir);
  return m_velocitySlopes[id * nDim * nDim + varId * nDim + dir];
}

template <MInt nDim>
MFloat LptCells<nDim>::velocitySlope(const MInt id, const MInt varId, const MInt dir) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dir);
  return m_velocitySlopes[id * nDim * nDim + varId * nDim + dir];
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void LptCells<nDim>::invalidate(const MInt begin, const MInt end) {
  fill_invalid(m_properties, begin, end, 1, false);
  fill_invalid(m_volumeFraction, begin, end);
  fill_invalid(m_noParticles, begin, end);
  fill_invalid(m_noEllipsoids, begin, end);
  fill_invalid(m_bndryCellId, begin, end);

  fill_invalid(m_variables, begin, end, noVars());
  fill_invalid(m_species, begin, end, noSpecies());

  if(hasMassCoupling()) {
    fill_invalid(m_massFlux, begin, end);
  }
  if(hasHeatCoupling()) fill_invalid(m_heatFlux, begin, end);
  if(hasMomentumCoupling()) {
    fill_invalid(m_workFlux, begin, end);
    fill_invalid(m_momentumFlux, begin, end, nDim);
  }
  if(hasSlopes()) {
    fill_invalid(m_velocitySlopes, begin, end, nDim * nDim);
  }
}

/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void LptCells<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                    const MInt destination) {
  copyData(source.m_properties, m_properties, c, begin, end, destination);
  copyData(source.m_volumeFraction, m_volumeFraction, c, begin, end, destination);
  copyData(source.m_noParticles, m_noParticles, c, begin, end, destination);
  copyData(source.m_noEllipsoids, m_noEllipsoids, c, begin, end, destination);
  copyData(source.m_bndryCellId, m_bndryCellId, c, begin, end, destination);

  copyData(source.m_variables, m_variables, c, begin, end, destination, noVars());
  copyData(source.m_species, m_species, c, begin, end, destination, noSpecies());

  if(hasMassCoupling()) {
    copyData(source.m_massFlux, m_massFlux, c, begin, end, destination);
  }
  if(hasHeatCoupling()) copyData(source.m_heatFlux, m_heatFlux, c, begin, end, destination);
  if(hasMomentumCoupling()) {
    copyData(source.m_workFlux, m_workFlux, c, begin, end, destination);
    copyData(source.m_momentumFlux, m_momentumFlux, c, begin, end, destination, nDim);
  }
  if(hasSlopes()) {
    copyData(source.m_velocitySlopes, m_velocitySlopes, c, begin, end, destination, nDim * nDim);
  }
}

} // namespace collector
} // namespace lpt
} // namespace maia


// Undefine macros that should not be used outside this file
#undef LPTCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_SET_ID_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef LPTCOLLECTOR_H_
