// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FCCOLLECTOR_H_
#define FCCOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
#include "fccellproperties.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define FCCOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'FCCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define FCCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(FCCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noVariables(),                                                                                 \
        "variable id = " + std::to_string(id) + " is out-of-bounds [0, " + std::to_string(noVariables()) + ")", AT_);  \
  } while(false)
#define ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(id)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noDistributions(),                                                           \
                          "distribution id = " + std::to_string(id) + " is out-of-bounds [0, "                         \
                              + std::to_string(noDistributions()),                                                     \
                          AT_);                                                                                        \
  } while(false)

#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != FcCell::NumProperties, "Invalid property", AT_);                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(id)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


// Namespace for auxiliary functions/classes
namespace maia {
namespace fc {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::fc::cell::BitsetType;

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
  static constexpr BitsetType value() { return BitsetType(0); }
};

/// Class that represents FC cell collector.
template <MInt nDim>
class FcCellCollector : public maia::container::Container<FcCellCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<FcCellCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<FcCellCollector<nDim>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;
  using BitsetType = maia::fc::cell::BitsetType;


 public:
  // Types
  template <class T>
  using Invalid = typename maia::fc::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr FcCellCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  MFloat& alpha(const MInt id);
  MFloat alpha(const MInt id) const;
  MFloat& poissonRatio(const MInt id);
  MFloat poissonRatio(const MInt id) const;
  MFloat& invJacobian(const MInt id);
  MFloat invJacobian(const MInt id) const;
  MInt& pRfnmnt(const MInt id);
  MInt pRfnmnt(const MInt id) const;
  MInt& maxSubCellLvl(const MInt id);
  MInt maxSubCellLvl(const MInt id) const;
  MInt& noNodesPerCell(const MInt id);
  MInt noNodesPerCell(const MInt id) const;
  MInt& bndId(const MInt id);
  MInt bndId(const MInt id) const;
  MInt& nodeIdsLoc(const MInt id, const MInt eid);
  MInt nodeIdsLoc(const MInt id, const MInt eid) const;
  MInt& nodeIdsGlob(const MInt id, const MInt eid);
  MInt nodeIdsGlob(const MInt id, const MInt eid) const;
  MFloat& elementDisplacements(const MInt id, const MInt eid);
  MFloat elementDisplacements(const MInt id, const MInt eid) const;
  MFloat& elementStrains(const MInt id, const MInt eid);
  MFloat elementStrains(const MInt id, const MInt eid) const;
  MFloat& nodalStrains(const MInt id, const MInt nid, const MInt eid);
  MFloat nodalStrains(const MInt id, const MInt nid, const MInt eid) const;
  MFloat& elementStresses(const MInt id, const MInt eid);
  MFloat elementStresses(const MInt id, const MInt eid) const;
  MFloat& nodalStresses(const MInt id, const MInt nid, const MInt eid);
  MFloat nodalStresses(const MInt id, const MInt nid, const MInt eid) const;
  MFloat& nodePosition(const MInt id, const MInt eid);
  MFloat nodePosition(const MInt id, const MInt eid) const;
  MFloat& deltaGamma(const MInt id, const MInt eid);
  MFloat deltaGamma(const MInt id, const MInt eid) const;
  MFloat& epsilonBarP(const MInt id, const MInt eid);
  MFloat epsilonBarP(const MInt id, const MInt eid) const;


  // Property-related accessors
  BitsetType::reference hasProperty(const MInt id, const FcCell p);
  MBool hasProperty(const MInt id, const FcCell p) const;
  void resetProperties(const MInt id);
  BitsetType& allProperties(const MInt id);

  /// Allow setting whether to support thermal computations
  void setThermal(const MBool isThermal_);

  /// Allow setting whether to support plastic computations
  void setPlasticity(const MBool isPlastic_);

  /// Sets the number of strains
  void setMaxPRfnmnt(const MInt maxPRfnmnt_);

  /// Return if simulation is thermal
  constexpr MBool isThermal() const { return m_isThermal; }

  /// Return if simulation is plastic
  constexpr MBool isPlastic() const { return m_isPlastic; }

  /// Return if simulation is p-refined
  constexpr MBool isPRefined() const { return m_isPRefined; }

  /// Return number of variables
  constexpr MInt noVariables() const { return m_noVariables; }

  /// Return number of strains
  constexpr MInt noStrains() const { return m_noStrainsPerCell; }

  /// Return number of stresses
  constexpr MInt noStresses() const { return m_noStressesPerCell; }

  constexpr MInt noNodes() const { return m_numberOfNodes; }

  /// Return max pRfnmnt
  constexpr MInt maxPRfnmnt() const { return m_maxPRfnmnt; }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::fc::cell::p(FcCell::NumProperties); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  /// Number of variables
  MInt m_noVariables = 1;

  MInt m_noStrainsPerCell = (nDim == 2) ? (nDim + 1) : (nDim * 2);

  MInt m_noStressesPerCell = m_noStrainsPerCell;

  MInt m_maxPRfnmnt = 0;

  MInt m_numberOfNodes = 0;

  /// Use thermal model
  MBool m_isThermal = false;

  /// Use plastic model
  MBool m_isPlastic = false;

  /// Use p-refinement
  MBool m_isPRefined = false;

  // Data containers
  Storage<MFloat> m_alpha{};
  Storage<MFloat> m_poissonRatio{};
  Storage<MFloat> m_invJacobian{};
  Storage<MInt> m_pRfnmnt{};
  Storage<MInt> m_maxSubCellLvl{};
  Storage<MInt> m_noNodesPerCell{};
  Storage<MInt> m_bndId{};
  Storage<MInt> m_nodeIdsLoc{};
  Storage<MInt> m_nodeIdsGlob{};
  Storage<MFloat> m_elementDisplacements{};
  Storage<MFloat> m_elementStrains{};
  Storage<MFloat> m_nodalStrains{};
  Storage<MFloat> m_elementStresses{};
  Storage<MFloat> m_nodalStresses{};
  Storage<MFloat> m_nodePosition{};
  Storage<MFloat> m_deltaGamma{};
  Storage<MFloat> m_epsilonBarP{};
  Storage<BitsetType> m_properties{};
};

/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void FcCellCollector<nDim>::reset() {
  resetStorage(1, m_alpha);
  resetStorage(1, m_poissonRatio);
  resetStorage(1, m_invJacobian);
  resetStorage(1, m_pRfnmnt);
  resetStorage(1, m_maxSubCellLvl);
  resetStorage(1, m_noNodesPerCell);
  resetStorage(1, m_bndId);
  resetStorage(noNodes(), m_nodeIdsLoc);
  resetStorage(noNodes(), m_nodeIdsGlob);
  resetStorage(nDim, m_elementDisplacements);
  resetStorage(noStrains(), m_elementStrains);
  resetStorage(maxPRfnmnt() + 2, m_nodePosition);
  resetStorage(noStresses(), m_elementStresses);
  if(isPlastic()) {
    resetStorage(noNodes(), m_deltaGamma);
    resetStorage(noNodes(), m_epsilonBarP);
    resetStorage(noStrains() * noNodes(), m_nodalStrains);
    resetStorage(noStresses() * noNodes(), m_nodalStresses);
  }

  resetStorage(1, m_properties);
}


/// Accessor for alpha.
template <MInt nDim>
MFloat& FcCellCollector<nDim>::alpha(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_alpha[id];
}
/// Accessor for alpha (const version).
template <MInt nDim>
MFloat FcCellCollector<nDim>::alpha(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_alpha[id];
}

/// Accessor for poissonRatio.
template <MInt nDim>
MFloat& FcCellCollector<nDim>::poissonRatio(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_poissonRatio[id];
}
/// Accessor for poissonRatio (const version).
template <MInt nDim>
MFloat FcCellCollector<nDim>::poissonRatio(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_poissonRatio[id];
}

/// Accessor for inverse jacobian.
template <MInt nDim>
MFloat& FcCellCollector<nDim>::invJacobian(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invJacobian[id];
}
/// Accessor for inverse jacobian (const version).
template <MInt nDim>
MFloat FcCellCollector<nDim>::invJacobian(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invJacobian[id];
}

/// Accessor for pRfnmnt
template <MInt nDim>
MInt& FcCellCollector<nDim>::pRfnmnt(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_pRfnmnt[id];
}
/// Accessor for pRfnmnt (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::pRfnmnt(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_pRfnmnt[id];
}

/// Accessor for maxSubCellLvl
template <MInt nDim>
MInt& FcCellCollector<nDim>::maxSubCellLvl(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_maxSubCellLvl[id];
}
/// Accessor for maxSubCellLvl (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::maxSubCellLvl(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_maxSubCellLvl[id];
}

/// Accessor for number of elements.
template <MInt nDim>
MInt& FcCellCollector<nDim>::noNodesPerCell(const MInt id) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noNodesPerCell[id];
}
/// Accessor for number of elements (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::noNodesPerCell(const MInt id) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noNodesPerCell[id];
}

/// Accessor for bndId.
template <MInt nDim>
MInt& FcCellCollector<nDim>::bndId(const MInt id) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndId[id];
}
/// Accessor for bndId (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::bndId(const MInt id) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndId[id];
}

/// Accessor for local nodeIds.
template <MInt nDim>
MInt& FcCellCollector<nDim>::nodeIdsLoc(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_nodeIdsLoc[id * noNodes() + eid];
}
/// Accessor for local nodeIds (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::nodeIdsLoc(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_nodeIdsLoc[id * noNodes() + eid];
}

/// Accessor for global nodeIds.
template <MInt nDim>
MInt& FcCellCollector<nDim>::nodeIdsGlob(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_nodeIdsGlob[id * noNodes() + eid];
}
/// Accessor for global nodeIds (const version).
template <MInt nDim>
MInt FcCellCollector<nDim>::nodeIdsGlob(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_nodeIdsGlob[id * noNodes() + eid];
}

/// Accessor for element displacements.
template <MInt nDim>
MFloat& FcCellCollector<nDim>::elementDisplacements(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementDisplacements[id * nDim + eid];
}
/// Accessor for element displacements (const version).
template <MInt nDim>
MFloat FcCellCollector<nDim>::elementDisplacements(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementDisplacements[id * nDim + eid];
}

/// Accessor for element strains, i.e., cell based strains.
template <MInt nDim>
MFloat& FcCellCollector<nDim>::elementStrains(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementStrains[id * noStrains() + eid];
}
/// Accessor for element strains, i.e., cell based strains (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::elementStrains(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementStrains[id * noStrains() + eid];
}

/// Accessor for nodal strains, i.e., node based strains
template <MInt nDim>
MFloat& FcCellCollector<nDim>::nodalStrains(const MInt id, const MInt nid, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(nid * noStrains() + eid);
  return m_nodalStrains[id * noStrains() * noNodes() + nid * noStrains() + eid];
}
/// Accessor for nodal strains, i.e., node based strains (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::nodalStrains(const MInt id, const MInt nid, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(nid * noStrains() + eid);
  return m_nodalStrains[id * noStrains() * noNodes() + nid * noStrains() + eid];
}

/// Accessor for element stresses, i.e., cell based stresses
template <MInt nDim>
MFloat& FcCellCollector<nDim>::elementStresses(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementStresses[id * noStresses() + eid];
}
/// Accessor for element stresses, i.e., cell based stresses (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::elementStresses(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_elementStresses[id * noStresses() + eid];
}

/// Accessor for nodal stresses, i.e., node based stresses
template <MInt nDim>
MFloat& FcCellCollector<nDim>::nodalStresses(const MInt id, const MInt nid, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(nid * noStresses() + eid);
  return m_nodalStresses[id * noStresses() * noNodes() + nid * noStresses() + eid];
}
/// Accessor for nodal stresses, i.e., node based stresses (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::nodalStresses(const MInt id, const MInt nid, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(nid * noStresses() + eid);
  return m_nodalStresses[id * noStresses() * noNodes() + nid * noStresses() + eid];
}

/// Accessor for the node position
template <MInt nDim>
MFloat& FcCellCollector<nDim>::nodePosition(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_nodePosition[id * (maxPRfnmnt() + 2) + eid];
}
/// Accessor for the node position (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::nodePosition(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_nodePosition[id * (maxPRfnmnt() + 2) + eid];
}

/// Accessor for delta gamma, a factor in plasticity models
template <MInt nDim>
MFloat& FcCellCollector<nDim>::deltaGamma(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_deltaGamma[id * noNodes() + eid];
}
/// Accessor for delta gamma, a factor in plasticity models (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::deltaGamma(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_deltaGamma[id * noNodes() + eid];
}

/// Accessor for epsilon bar, a factor in plasticity models
template <MInt nDim>
MFloat& FcCellCollector<nDim>::epsilonBarP(const MInt id, const MInt eid) {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_epsilonBarP[id * noNodes() + eid];
}
/// Accessor for epsilon bar, a factor in plasticity models (const version)
template <MInt nDim>
MFloat FcCellCollector<nDim>::epsilonBarP(const MInt id, const MInt eid) const {
  // Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_epsilonBarP[id * noNodes() + eid];
}

/// Accessor for properties.
template <MInt nDim>
FcCellCollector<nDim>::BitsetType::reference FcCellCollector<nDim>::hasProperty(const MInt id, const FcCell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::fc::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool FcCellCollector<nDim>::hasProperty(const MInt id, const FcCell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::fc::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void FcCellCollector<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties[id].reset();
}

/// Accessor for properties.
template <MInt nDim>
BitsetType& FcCellCollector<nDim>::allProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}

/// Set use of thermal model and update number of variables
template <MInt nDim>
void FcCellCollector<nDim>::setThermal(const MBool isThermal_) {
  m_isThermal = isThermal_;
  m_noVariables = isThermal() ? 1 + 1 : 1;
}

/// Set use of plastic model and update number of variables
template <MInt nDim>
void FcCellCollector<nDim>::setPlasticity(const MBool isPlastic_) {
  m_isPlastic = isPlastic_;
}

/// Set use of p-refinement and update number of variables
template <MInt nDim>
void FcCellCollector<nDim>::setMaxPRfnmnt(const MInt maxPRfnmnt_) {
  m_maxPRfnmnt = maxPRfnmnt_;
  m_numberOfNodes = (m_maxPRfnmnt + 2) * (m_maxPRfnmnt + 2);
  if(nDim == 3) m_numberOfNodes *= (m_maxPRfnmnt + 2);

  if(maxPRfnmnt_ > 0) m_isPRefined = true;
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void FcCellCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  fill_invalid(m_alpha, begin, end, 1, 1.);
  fill_invalid(m_poissonRatio, begin, end, 1, 0.);
  fill_invalid(m_invJacobian, begin, end, 1, 0.);
  fill_invalid(m_pRfnmnt, begin, end, 1, 0);
  fill_invalid(m_maxSubCellLvl, begin, end, 1, 0);
  fill_invalid(m_noNodesPerCell, begin, end, 1, -1);
  fill_invalid(m_bndId, begin, end, 1, -1);
  fill_invalid(m_nodeIdsLoc, begin, end, noNodes(), -1);
  fill_invalid(m_nodeIdsGlob, begin, end, noNodes(), -1);
  fill_invalid(m_elementDisplacements, begin, end, nDim, 0.);
  fill_invalid(m_elementStrains, begin, end, noStrains(), 0.);
  fill_invalid(m_nodePosition, begin, end, maxPRfnmnt() + 2, 0.);
  fill_invalid(m_elementStresses, begin, end, noStresses(), 0.);

  if(isPlastic()) {
    fill_invalid(m_deltaGamma, begin, end, noNodes(), 0.);
    fill_invalid(m_epsilonBarP, begin, end, noNodes(), 0.);
    fill_invalid(m_nodalStrains, begin, end, noStrains() * noNodes(), 0.);
    fill_invalid(m_nodalStresses, begin, end, noStresses() * noNodes(), 0.);
  }

  // Properties
  fill_invalid(m_properties, begin, end);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void FcCellCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                           const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef FCCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  copyData(source.m_alpha, m_alpha, c, begin, end, destination);
  copyData(source.m_poissonRatio, m_poissonRatio, c, begin, end, destination);
  copyData(source.m_invJacobian, m_invJacobian, c, begin, end, destination);
  copyData(source.m_pRfnmnt, m_pRfnmnt, c, begin, end, destination);
  copyData(source.m_maxSubCellLvl, m_maxSubCellLvl, c, begin, end, destination);
  copyData(source.m_noNodesPerCell, m_noNodesPerCell, c, begin, end, destination);
  copyData(source.m_bndId, m_bndId, c, begin, end, destination);
  copyData(source.m_nodeIdsLoc, m_nodeIdsLoc, c, begin, end, destination, noNodes());
  copyData(source.m_nodeIdsGlob, m_nodeIdsGlob, c, begin, end, destination, noNodes());
  copyData(source.m_elementDisplacements, m_elementDisplacements, c, begin, end, destination, nDim);
  copyData(source.m_elementStrains, m_elementStrains, c, begin, end, destination, noStrains());
  copyData(source.m_nodePosition, m_nodePosition, c, begin, end, destination, maxPRfnmnt() + 2);
  copyData(source.m_elementStresses, m_elementStresses, c, begin, end, destination, noStresses());

  if(isPlastic()) {
    copyData(source.m_deltaGamma, m_deltaGamma, c, begin, end, destination);
    copyData(source.m_epsilonBarP, m_epsilonBarP, c, begin, end, destination);
    copyData(source.m_nodalStrains, m_nodalStrains, c, begin, end, destination);
    copyData(source.m_nodalStresses, m_nodalStresses, c, begin, end, destination);
  }


  // Properties
  copyData(source.m_properties, m_properties, c, begin, end, destination);
}

} // namespace collector
} // namespace fc
} // namespace maia


// Undefine macros that should not be used outside this file
#undef FCCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef FCCOLLECTOR_H_
