// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSCOLLECTOR_H_
#define LSCOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
#include "lscartesiancellproperties.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define LSCOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'LSCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity checks
// for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define LSCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

#ifdef LSCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

// Sanity-checking macros for accessors
#if defined(MAIA_ASSERT_ACCESSORS) || defined(LSCOLLECTOR_SANITY_CHECKS_ACCESSORS)
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

#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != LsCell::NumProperties, "Invalid property", AT_);                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SET_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_BODY_ID_ACCESSOR(id)                                                                              \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


// Namespace for auxiliary functions/classes
namespace maia {
namespace ls {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::ls::cell::BitsetType;
using BitsetTypeSet = maia::ls::cell::BitsetTypeSet;


// Type traits for invalid values. These values are used to initialize/erase nodes
template <class T>
struct Invalid {};

// Invalid value for ids is 'INT_MIN'
template <>
struct Invalid<MInt> {
  static constexpr MInt value() { return std::numeric_limits<MInt>::min(); }
};

// Invalid value for ids is 'LONG_INT_MIN'
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

// Invalid value for BitsetProperties is '0'
template <>
struct Invalid<BitsetType> {
  static constexpr BitsetType value() { return BitsetType(0); }
};
/* // Invalid value for BitsetProperties is '0' */
/* template <> */
/* struct Invalid<BitsetTypeSet> { */
/* static constexpr BitsetTypeSet value() { return BitsetTypeSet(0); } */
/* }; */
// NOTE: additional declaration only necessary if:
//      NumProperties != NumSetProperties is lscellproperties.h
// Invalid value for BitsetProperties is '0'
// template <> struct Invalid<BitsetTypeSet> {
//  static constexpr BitsetTypeSet value() { return BitsetTypeSet(0); }
//};


/// Class that represents LS cell collector.
template <MInt nDim>
class GCells : public maia::container::Container<GCells<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<GCells<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<GCells<nDim>, Invalid>;
  using Base::resetStorage;
  using self = GCells<nDim>;
  template <class T>
  using Storage = typename Base::template Storage<T>;


 public:
  using BitsetType = maia::ls::cell::BitsetType;
  using BitsetTypeSet = maia::ls::cell::BitsetTypeSet;


  // Types
  template <class T>
  using Invalid = typename maia::ls::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr GCells() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  void resetExtensionVelocity();
  void fillContainingCell();
  void fillContainingDomain();

  // accessors to specify GCell collector type
  constexpr MBool hasBodyId() const { return m_hasBodyId; }
  constexpr MBool hasOldG() const { return m_hasOldG; }
  constexpr MBool hasRotatingLs() const { return m_hasRotatingLs; }
  constexpr MBool hasReconstructOldG() const { return m_hasReconstructOldG; }
  constexpr MBool hasFSRhs() const { return m_hasFSRhs; }
  constexpr MBool hasReinit() const { return m_hasReinit; }
  constexpr MBool hasCorrectedBurningVelocity() const { return m_hasCorrectedBurningVelocity; }

  constexpr MBool hasFExt() const { return m_hasFExt; }
  constexpr MBool hasCurvature() const { return m_hasCurvature; }
  constexpr MBool hasNormalVector() const { return m_hasNormal; }
  constexpr MBool hasGapClosing() const { return m_hasGapClosing; }

  /// Accessors
  // ls-value
  MFloat& gFunction(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_gFunction.at(m_maxNoSets * id + set);
  }
  MFloat gFunction(const MInt id, const MInt set) const { return const_cast<self&>(*this).gFunction(id, set); }
  // old ls-value (for semi-Lagrange levelset!)
  MFloat& oldGFunction(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_oldGFunction.at(m_maxNoSets * id + set);
  }
  MFloat oldGFunction(const MInt id, const MInt set) const { return const_cast<self&>(*this).oldGFunction(id, set); }

  // ls-function slope
  MFloat& levelSetFunctionSlope(const MInt id, const MInt dim, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_levelSetFunctionSlope.at((id * m_maxNoSets + set) * nDim + dim);
  }
  MFloat levelSetFunctionSlope(const MInt id, const MInt dim, const MInt set) const {
    return const_cast<self&>(*this).levelSetFunctionSlope(id, dim, set);
  }
  // ls-rhs
  MFloat& levelSetRHS(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_levelSetRHS.at(m_maxNoSets * id + set);
  }
  MFloat levelSetRHS(const MInt id, const MInt set) const { return const_cast<self&>(*this).levelSetRHS(id, set); }

  // corrected burning velocity
  MFloat& correctedBurningVelocity(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_correctedBurningVelocity.at(m_maxNoSets * id + set);
  }
  MFloat correctedBurningVelocity(const MInt id, const MInt set) const {
    return const_cast<self&>(*this).correctedBurningVelocity(id, set);
  }

  // Containing Cell
  MLong& containingCell(const MInt id, const MInt body) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_BODY_ID_ACCESSOR(body);
    return m_containingCell.at(m_maxBodies * id + body);
  }
  MLong containingCell(const MInt id, const MInt body) const {
    return const_cast<self&>(*this).containingCell(id, body);
  }

  // Containg Domain
  MInt& containingDomain(const MInt id, const MInt body) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_BODY_ID_ACCESSOR(body);
    return m_containingDomain.at(m_maxBodies * id + body);
  }
  MInt containingDomain(const MInt id, const MInt body) const {
    return const_cast<self&>(*this).containingDomain(id, body);
  }

  // curvature
  MFloat& curvature(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_curvature.at(m_maxNoSets * id + set);
  }
  MFloat curvature(const MInt id, const MInt set) const { return const_cast<self&>(*this).curvature(id, set); }
  // extension-velocity
  MFloat& fExt(const MInt id, const MInt dim, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_fExt.at((id * m_maxNoSets + set) * nDim + dim);
  }
  MFloat fExt(const MInt id, const MInt dim, const MInt set) const {
    return const_cast<self&>(*this).fExt(id, dim, set);
  }
  // normal vector
  MFloat& normalVector(const MInt id, const MInt dim, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_normalVector.at((id * m_maxNoSets + set) * nDim + dim);
  }
  MFloat normalVector(const MInt id, const MInt dim, const MInt set) const {
    return const_cast<self&>(*this).normalVector(id, dim, set);
  }
  // bodyId
  MInt& bodyId(const MInt id, const MInt set) {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_bodyId.at(m_maxNoSets * id + set);
  }
  MInt bodyId(const MInt id, const MInt set) const {
    ENSURE_VALID_ID_ACCESSOR(id);
    ENSURE_VALID_SET_ID_ACCESSOR(set);
    return m_bodyId.at(m_maxNoSets * id + set);
  }

  // secondBodyId
  MInt& secondBodyId(const MInt id) {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_secondBodyId.at(id);
  }
  MInt secondBodyId(const MInt id) const {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_secondBodyId.at(id);
  }

  // gap width
  MFloat& gapWidth(const MInt id) {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_gapWidth.at(id);
  }
  MFloat gapWidth(const MInt id) const {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_gapWidth.at(id);
  }

  // potentialGapCell
  MInt& potentialGapCell(const MInt id) {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_potentialGapCell.at(id);
  }
  MInt potentialGapCell(const MInt id) const {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_potentialGapCell.at(id);
  }
  // potentialGapCellClose
  MInt& potentialGapCellClose(const MInt id) {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_potentialGapCellClose.at(id);
  }
  MInt potentialGapCellClose(const MInt id) const {
    ENSURE_VALID_ID_ACCESSOR(id);
    return m_potentialGapCellClose.at(id);
  }

  /// properties:
  MBool regridTrigger(const MInt id) const { return hasProperty(id, LsCell::RegridTrigger); }
  BitsetType::reference regridTrigger(const MInt id) { return hasProperty(id, LsCell::RegridTrigger); }
  MBool nearGap(const MInt id) const { return hasProperty(id, LsCell::NearGap); }
  BitsetType::reference nearGap(const MInt id) { return hasProperty(id, LsCell::NearGap); }
  MBool isBndryG(const MInt id) const { return hasProperty(id, LsCell::IsBndryG); }
  BitsetType::reference isBndryG(const MInt id) { return hasProperty(id, LsCell::IsBndryG); }

  MBool inBand(const MInt id, const MInt set) const { return hasSetProperty(id, set, LsSet::InBand); }
  BitsetTypeSet::reference inBand(const MInt id, const MInt set) { return hasSetProperty(id, set, LsSet::InBand); }
  MBool isGBndryCell(const MInt id, const MInt set) const { return hasSetProperty(id, set, LsSet::IsGBndryCell); }
  BitsetTypeSet::reference isGBndryCell(const MInt id, const MInt set) {
    return hasSetProperty(id, set, LsSet::IsGBndryCell);
  }
  MBool isGZero(const MInt id, const MInt set) const { return hasSetProperty(id, set, LsSet::IsGZero); }
  BitsetTypeSet::reference isGZero(const MInt id, const MInt set) { return hasSetProperty(id, set, LsSet::IsGZero); }
  MBool wasGZero(const MInt id, const MInt set) const { return hasSetProperty(id, set, LsSet::WasGZero); }
  BitsetTypeSet::reference wasGZero(const MInt id, const MInt set) { return hasSetProperty(id, set, LsSet::WasGZero); }
  MBool hasPositiveSign(const MInt id, const MInt set) const { return hasSetProperty(id, set, LsSet::HasPositiveSign); }
  BitsetTypeSet::reference hasPositiveSign(const MInt id, const MInt set) {
    return hasSetProperty(id, set, LsSet::HasPositiveSign);
  }

  // Property-related accessors
 public:
  BitsetType::reference hasProperty(const MInt id, const LsCell p);
  MBool hasProperty(const MInt id, const LsCell p) const;
  void resetProperties(const MInt id);

  BitsetTypeSet::reference hasSetProperty(const MInt id, const MInt set, const LsSet p);
  MBool hasSetProperty(const MInt id, const MInt set, const LsSet p) const;
  void resetSetProperties(const MInt id, const MInt set);


 public:
  /// Return number of level set functions
  MInt maxNoSets() const { return m_maxNoSets; }
  MInt maxBodies() const { return m_maxBodies; }

  void setMaxNoSets(const MInt value) {
    ASSERT(value >= 0, "setMaxNoSets size must be >= 0");
    m_maxNoSets = value;
  }

  void setMaxBodiesToCompute(const MInt noBodies) {
    ASSERT(noBodies >= 0, "setMaxBodiesToCompute size must be >= 0");
    m_maxBodies = noBodies;
  }

  void setRotatingLs(const MBool rotatingLS) { m_hasRotatingLs = rotatingLS; }

  void setGapClosing(const MBool gapClosing) { m_hasGapClosing = gapClosing; }

  void setReconstructOldG(const MBool ReconstructOldG) { m_hasReconstructOldG = ReconstructOldG; }

  void setReinit(const MBool Reinit) { m_hasReinit = Reinit; }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::ls::cell::p(LsCell::NumProperties); }
  static constexpr MInt noSetProperties() { return maia::ls::cell::setP(LsSet::NumSetProperties); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

 public:
  void setLsCollectorType(const MInt mode);

 private:
  /// Maximum number of sets
  MInt m_maxNoSets = 0;
  MInt m_maxBodies = 0;

  /// Bools for ls-collector types
  MBool m_hasBodyId = false;
  MBool m_hasOldG = false;

  MBool m_hasRotatingLs = false;
  MBool m_hasReconstructOldG = false;

  MBool m_hasFSRhs = false;
  MBool m_hasReinit = false;

  MBool m_hasCorrectedBurningVelocity = false;

  MBool m_hasFExt = false;
  MBool m_hasNormal = false;
  MBool m_hasCurvature = false;

  MBool m_hasGapClosing = false;

  // Data containers
  Storage<MFloat> m_gFunction;

  Storage<MLong> m_containingCell; // Has to be MLong because containing cells are converted to globalIds during restart
                                   // and load balance.
  Storage<MInt> m_containingDomain;

  Storage<MFloat> m_curvature;
  Storage<MFloat> m_fExt;
  Storage<MFloat> m_normalVector;

  Storage<MInt> m_bodyId;
  Storage<MInt> m_secondBodyId;
  Storage<MFloat> m_oldGFunction;

  Storage<MFloat> m_levelSetFunctionSlope;
  Storage<MFloat> m_levelSetRHS;

  Storage<MFloat> m_correctedBurningVelocity;

  Storage<BitsetType> m_properties;
  Storage<BitsetTypeSet> m_setProperties;

  Storage<MFloat> m_gapWidth;
  Storage<MInt> m_potentialGapCell;
  Storage<MInt> m_potentialGapCellClose;
};

/// Set the ls-collector type
template <MInt nDim>
void GCells<nDim>::setLsCollectorType(const MInt mode) {
  ASSERT(mode >= 0, "");
  switch(mode) {
    case 0: // combustion:
    {
      m_hasFExt = true;
      m_hasNormal = true;
      m_hasCurvature = true;
      m_hasOldG = true;
      m_hasFSRhs = true;
      m_hasCorrectedBurningVelocity = true;
      ASSERT(!m_hasBodyId, "");
      break;
    }
    case 1: // semi-lagrange moving bndry without reinitialisation
    {
      m_hasBodyId = true;
      m_hasOldG = true;
      ASSERT(!m_hasFExt, "");
      ASSERT(!m_hasNormal, "");
      ASSERT(!m_hasCurvature, "");
      ASSERT(!m_hasFSRhs, "");
      ASSERT(!m_hasCorrectedBurningVelocity, "");
      break;
    }
    case 2: // semi-lagrange moving bndry with reinitialisation
    {
      m_hasBodyId = true;
      m_hasOldG = true;
      m_hasNormal = true;
      m_hasCurvature = true;
      ASSERT(!m_hasFExt, "");
      ASSERT(!m_hasFSRhs, "");
      ASSERT(!m_hasCorrectedBurningVelocity, "");
      break;
    }
    case 3: // moving bndry with flameExt-velocity
    {
      m_hasBodyId = true;
      m_hasNormal = true;
      m_hasCurvature = true;
      m_hasFExt = true;
      m_hasOldG = true;
      m_hasFSRhs = true;
      m_hasCorrectedBurningVelocity = true;
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown lsCollector Type!");
    }
  }
}


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void GCells<nDim>::reset() {
  // resetStorage(#of items per cell, Storage<> name);
  resetStorage(m_maxNoSets, m_gFunction);

  if(hasCurvature()) resetStorage(m_maxNoSets, m_curvature);
  if(hasFExt()) resetStorage(m_maxNoSets * nDim, m_fExt);
  if(hasNormalVector()) resetStorage(m_maxNoSets * nDim, m_normalVector);
  if(hasBodyId()) {
    resetStorage(m_maxNoSets, m_bodyId);
    resetStorage(1, m_secondBodyId);
  }
  if(hasOldG()) resetStorage(m_maxNoSets, m_oldGFunction);

  if(hasFSRhs() || hasReinit()) {
    resetStorage(m_maxNoSets * nDim, m_levelSetFunctionSlope);
    resetStorage(m_maxNoSets, m_levelSetRHS);
  }

  if(hasCorrectedBurningVelocity()) {
    resetStorage(m_maxNoSets, m_correctedBurningVelocity);
  }

  if(hasRotatingLs()) resetStorage(m_maxBodies, m_containingCell);
  if(hasRotatingLs() && !hasReconstructOldG()) resetStorage(m_maxBodies, m_containingDomain);

  if(hasGapClosing()) {
    resetStorage(1, m_gapWidth);
    resetStorage(1, m_potentialGapCell);
    resetStorage(1, m_potentialGapCellClose);
  }

  resetStorage(1, m_properties);
  resetStorage(m_maxNoSets, m_setProperties);
}

template <MInt nDim>
void GCells<nDim>::resetExtensionVelocity() {
  ASSERT(hasFExt(), "Extension velocity not activated in cell collector!");
  fill_invalid(m_fExt, 0, this->size(), m_maxNoSets * nDim, 0.0);
}

template <MInt nDim>
void GCells<nDim>::fillContainingCell() {
  ASSERT(hasRotatingLs(), "non-rotating");
  fill_invalid(m_containingCell, 0, this->size(), m_maxBodies, -1);
}

template <MInt nDim>
void GCells<nDim>::fillContainingDomain() {
  ASSERT(hasRotatingLs(), "non-rotating");
  ASSERT(!hasReconstructOldG(), "ReconstructOldG");
  fill_invalid(m_containingDomain, 0, this->size(), m_maxBodies, -1);
}

/// Accessor for properties.
template <MInt nDim>
GCells<nDim>::BitsetType::reference GCells<nDim>::hasProperty(const MInt id, const LsCell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties.at(id)[maia::ls::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool GCells<nDim>::hasProperty(const MInt id, const LsCell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties.at(id)[maia::ls::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void GCells<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties.at(id).reset();
}

/// Accessor for set properties.
template <MInt nDim>
GCells<nDim>::BitsetTypeSet::reference GCells<nDim>::hasSetProperty(const MInt id, const MInt set, const LsSet p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  // ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_setProperties.at(m_maxNoSets * id + set)[maia::ls::cell::setP(p)];
}
/// Accessor for set properties (const version).
template <MInt nDim>
MBool GCells<nDim>::hasSetProperty(const MInt id, const MInt set, const LsSet p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  // ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_setProperties.at(m_maxNoSets * id + set)[maia::ls::cell::setP(p)];
}
/// Reset all set properties.
template <MInt nDim>
void GCells<nDim>::resetSetProperties(const MInt id, const MInt set) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_setProperties.at(m_maxNoSets * id + set).reset();
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void GCells<nDim>::invalidate(const MInt begin, const MInt end) {
  fill_invalid(m_gFunction, begin, end, m_maxNoSets);

  if(hasCurvature()) fill_invalid(m_curvature, begin, end, m_maxNoSets, 0);
  if(hasFExt()) fill_invalid(m_fExt, begin, end, m_maxNoSets * nDim);
  if(hasNormalVector()) fill_invalid(m_normalVector, begin, end, m_maxNoSets * nDim);
  if(hasBodyId()) {
    fill_invalid(m_bodyId, begin, end, m_maxNoSets, -1);
    fill_invalid(m_secondBodyId, begin, end, 1, -1);
  }
  if(hasOldG()) fill_invalid(m_oldGFunction, begin, end, m_maxNoSets, 0);

  if(hasFSRhs() || hasReinit()) {
    fill_invalid(m_levelSetFunctionSlope, begin, end, m_maxNoSets * nDim);
    fill_invalid(m_levelSetRHS, begin, end, m_maxNoSets);
  }

  if(hasCorrectedBurningVelocity()) {
    fill_invalid(m_correctedBurningVelocity, begin, end, m_maxNoSets);
  }

  if(hasRotatingLs()) fill_invalid(m_containingCell, begin, end, m_maxBodies, -1);
  if(hasRotatingLs() && !hasReconstructOldG()) fill_invalid(m_containingDomain, begin, end, m_maxBodies, -1);

  if(hasGapClosing()) {
    fill_invalid(m_gapWidth, begin, end, 1);
    fill_invalid(m_potentialGapCell, begin, end, 1, 0);
    fill_invalid(m_potentialGapCellClose, begin, end, 1, 0);
  }

  fill_invalid(m_properties, begin, end, 1, false);
  fill_invalid(m_setProperties, begin, end, m_maxNoSets, false);
}

/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void GCells<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                  const MInt destination) {
  copyData(source.m_gFunction, m_gFunction, c, begin, end, destination, m_maxNoSets);

  if(hasCurvature()) {
    copyData(source.m_curvature, m_curvature, c, begin, end, destination, m_maxNoSets);
  }
  if(hasFExt()) {
    copyData(source.m_fExt, m_fExt, c, begin, end, destination, m_maxNoSets * nDim);
  }
  if(hasNormalVector()) {
    copyData(source.m_normalVector, m_normalVector, c, begin, end, destination, m_maxNoSets * nDim);
  }
  if(hasBodyId()) {
    copyData(source.m_bodyId, m_bodyId, c, begin, end, destination, m_maxNoSets);
    copyData(source.m_secondBodyId, m_secondBodyId, c, begin, end, destination, 1);
  }
  if(hasOldG()) {
    copyData(source.m_oldGFunction, m_oldGFunction, c, begin, end, destination, m_maxNoSets);
  }

  if(hasFSRhs() || hasReinit()) {
    copyData(source.m_levelSetFunctionSlope, m_levelSetFunctionSlope, c, begin, end, destination, m_maxNoSets * nDim);
    copyData(source.m_levelSetRHS, m_levelSetRHS, c, begin, end, destination, m_maxNoSets);
  }

  if(hasCorrectedBurningVelocity()) {
    copyData(source.m_correctedBurningVelocity, m_correctedBurningVelocity, c, begin, end, destination, m_maxNoSets);
  }

  if(hasRotatingLs()) copyData(source.m_containingCell, m_containingCell, c, begin, end, destination, m_maxBodies);
  if(hasRotatingLs() && !hasReconstructOldG())
    copyData(source.m_containingDomain, m_containingDomain, c, begin, end, destination, m_maxBodies);

  if(hasGapClosing()) {
    copyData(source.m_gapWidth, m_gapWidth, c, begin, end, destination, 1);
    copyData(source.m_potentialGapCell, m_potentialGapCell, c, begin, end, destination, 1);
    copyData(source.m_potentialGapCellClose, m_potentialGapCellClose, c, begin, end, destination, 1);
  }

  copyData(source.m_properties, m_properties, c, begin, end, destination);
  copyData(source.m_setProperties, m_setProperties, c, begin, end, destination, m_maxNoSets);
}

} // namespace collector
} // namespace ls
} // namespace maia


// Undefine macros that should not be used outside this file
#undef LSCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_SET_ID_ACCESSOR
#undef ENSURE_VALID_BODY_ID_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef LSCOLLECTOR_H_
