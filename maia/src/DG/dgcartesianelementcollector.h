// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGELEMENTCOLLECTOR_H_
#define DGELEMENTCOLLECTOR_H_

#include <algorithm>
#include <limits>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define DGCOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'DGCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks
// for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define DGCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(DGCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noVars(),                                                                                      \
        "variable id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noVars()) + ")", AT_);            \
  } while(false)
#define ENSURE_VALID_NODE_VARIABLE_ID_ACCESSOR(id)                                                                     \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noNodeVars(),                                                                                  \
        "node variable id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noNodeVars()) + ")", AT_);   \
  } while(false)
#define ENSURE_VALID_SURFACE_ID_ACCESSOR(id)                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noSurfaces(),                                                                \
                          "surface id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noSurfaces()),   \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(condition, message, AT_);                                                                    \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_NODE_VARIABLE_ID_ACCESSOR(id)                                                                     \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SURFACE_ID_ACCESSOR(id)                                                                           \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace dg {
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


/// Class that represents DG element collector.
template <MInt nDim, class SysEqn>
class ElementCollector : public maia::container::Container<ElementCollector<nDim, SysEqn>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<ElementCollector<nDim, SysEqn>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<ElementCollector<nDim, SysEqn>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  // Types
  template <class T>
  using Invalid = typename maia::dg::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr ElementCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  MInt& cellId(const MInt id);
  MInt cellId(const MInt id) const;
  MInt& polyDeg(const MInt id);
  MInt polyDeg(const MInt id) const;
  MInt noNodes1D(const MInt id) const;
  MInt& noNodes1D(const MInt id);
  MInt noNodesXD(const MInt id) const;
  MInt& surfaceIds(const MInt id, const MInt dir);
  MInt surfaceIds(const MInt id, const MInt dir) const;
  MFloat& nodeCoordinates(const MInt id);
  MFloat nodeCoordinates(const MInt id) const;
  MFloat& variables(const MInt id, const MInt pos);
  MFloat& variables(const MInt id) { return variables(id, 0); };
  MFloat variables(const MInt id, const MInt pos) const;
  MFloat variables(const MInt id) const { return variables(id, 0); };
  MFloat& timeIntStorage(const MInt id);
  MFloat timeIntStorage(const MInt id) const;
  MFloat& nodeVars(const MInt id);
  MFloat nodeVars(const MInt id) const;
  MFloat& rightHandSide(const MInt id);
  MFloat rightHandSide(const MInt id) const;
  MFloat& externalSource(const MInt id);
  MFloat externalSource(const MInt id) const;
  MFloat& invJacobian(const MInt id);
  MFloat invJacobian(const MInt id) const;

  /// Return maximum polynomial degree
  MInt maxPolyDeg() const { return m_maxPolyDeg; }
  // Set maximum polynomial degree
  MInt maxPolyDeg(const MInt maxPolyDeg_);

  /// Return maximum number of nodes 1D
  MInt maxNoNodes1D() const { return m_maxNoNodes1D; }
  // Set maximum number of nodes 1D
  void maxNoNodes1D(const MInt maxNoNodesDeg_);
  /// Return maximum number of nodes XD
  MInt maxNoNodesXD() const { return m_maxNoNodesXD; }

  // Set number of node variables
  void noNodeVars(const MInt noNodeVars_) { m_noNodeVars = noNodeVars_; };

  /// Return number of nodes in 1D
  // constexpr MInt noNodes1D() const { return m_noNodes1D; }

  /// Return number of nodes in 2D/3D
  // constexpr MInt noNodesXD() const { return m_noNodesXD; }

  /// Return number of variables
  static constexpr MInt noVars() { return SysEqn::noVars(); }

  /// Return number of node variables
  constexpr MInt noNodeVars() const { return m_noNodeVars; }

  /// Return number of surfaces per element
  static constexpr MInt noSurfaces() { return 2 * nDim; }

  // Auxiliary methods
  MInt getElementByCellId(const MInt cellId) const;
  MInt getElementByCellId(const MInt first, const MInt last, const MInt cellId) const;

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  /// Maximum polynomial degree
  MInt m_maxPolyDeg = -1;

  /// Maximum number of nodes 1D
  MInt m_maxNoNodes1D = -1;

  /// Maximum number if nodes XD
  MInt m_maxNoNodesXD = -1;

  // Number of node variables
  MInt m_noNodeVars = -1;

  // Data containers
  Storage<MInt> m_cellId{};
  Storage<MInt> m_polyDeg{};
  Storage<MInt> m_noNodes1D{};
  Storage<MInt> m_surfaceIds{};
  Storage<MFloat> m_nodeCoordinates{};
  Storage<MFloat> m_variables{};
  Storage<MFloat> m_timeIntStorage{};
  Storage<MFloat> m_nodeVariables{};
  Storage<MFloat> m_rightHandSide{};
  Storage<MFloat> m_externalSource{};
  Storage<MFloat> m_invJacobian{};
};


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim, class SysEqn>
void ElementCollector<nDim, SysEqn>::reset() {
  resetStorage(1, m_cellId);
  resetStorage(1, m_polyDeg);
  resetStorage(1, m_noNodes1D);
  resetStorage(noSurfaces(), m_surfaceIds);
  resetStorage(maxNoNodesXD() * nDim, m_nodeCoordinates);
  resetStorage(maxNoNodesXD() * noVars(), m_variables);
  resetStorage(maxNoNodesXD() * noVars(), m_timeIntStorage);
  resetStorage(maxNoNodesXD() * noNodeVars(), m_nodeVariables);
  resetStorage(maxNoNodesXD() * noVars(), m_rightHandSide);
  resetStorage(maxNoNodesXD() * noVars(), m_externalSource);
  resetStorage(1, m_invJacobian);
}


/// Accessor for cell id.
template <MInt nDim, class SysEqn>
MInt& ElementCollector<nDim, SysEqn>::cellId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellId[id];
}
/// Accessor for cell id (const version).
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::cellId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellId[id];
}


/// Accessor for polynomial degree.
template <MInt nDim, class SysEqn>
MInt& ElementCollector<nDim, SysEqn>::polyDeg(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_polyDeg[id];
}
/// Accessor for polynomial degree (const version).
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::polyDeg(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_polyDeg[id];
}

/// Accessor for number of nodes 1D
template <MInt nDim, class SysEqn>
MInt& ElementCollector<nDim, SysEqn>::noNodes1D(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noNodes1D[id];
}

/// Accessor for number of nodes 1D (const version).
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::noNodes1D(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noNodes1D[id];
}

/// Accessor for number of nodes XD (const version).
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::noNodesXD(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return ipow(m_noNodes1D[id], nDim);
}

/// Accessor for surface ids.
template <MInt nDim, class SysEqn>
MInt& ElementCollector<nDim, SysEqn>::surfaceIds(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SURFACE_ID_ACCESSOR(dir);
  return m_surfaceIds[id * noSurfaces() + dir];
}
/// Accessor for surface ids (const version).
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::surfaceIds(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SURFACE_ID_ACCESSOR(dir);
  return m_surfaceIds[id * noSurfaces() + dir];
}


/// Accessor for node coordinates.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::nodeCoordinates(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nodeCoordinates[id * nDim * maxNoNodesXD()];
}
/// Accessor for node coordinates (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::nodeCoordinates(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nodeCoordinates[id * nDim * maxNoNodesXD()];
}


/// Accessor for variables.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::variables(const MInt id, const MInt pos) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_CONDITION(pos >= 0 && pos < noVars() * maxNoNodesXD(), "Invalid position.");
  return m_variables[id * noVars() * maxNoNodesXD() + pos];
}
/// Accessor for variables (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::variables(const MInt id, const MInt pos) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_CONDITION(pos >= 0 && pos < noVars() * maxNoNodesXD(), "Invalid position.");
  return m_variables[id * noVars() * maxNoNodesXD() + pos];
}


/// Accessor for storage variables.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::timeIntStorage(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_timeIntStorage[id * noVars() * maxNoNodesXD()];
}
/// Accessor for storage variables (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::timeIntStorage(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_timeIntStorage[id * noVars() * maxNoNodesXD()];
}


/// Accessor for node variables.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::nodeVars(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nodeVariables[id * noNodeVars() * maxNoNodesXD()];
}
/// Accessor for node variables (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::nodeVars(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nodeVariables[id * noNodeVars() * maxNoNodesXD()];
}


/// Accessor for right hand side.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::rightHandSide(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_rightHandSide[id * noVars() * maxNoNodesXD()];
}
/// Accessor for right hand side (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::rightHandSide(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_rightHandSide[id * noVars() * maxNoNodesXD()];
}


/// Accessor for external source terms.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::externalSource(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_externalSource[id * noVars() * maxNoNodesXD()];
}
/// Accessor for external source terms (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::externalSource(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_externalSource[id * noVars() * maxNoNodesXD()];
}


/// Accessor for inverse jacobian.
template <MInt nDim, class SysEqn>
MFloat& ElementCollector<nDim, SysEqn>::invJacobian(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invJacobian[id];
}
/// Accessor for inverse jacobian (const version).
template <MInt nDim, class SysEqn>
MFloat ElementCollector<nDim, SysEqn>::invJacobian(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invJacobian[id];
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim, class SysEqn>
void ElementCollector<nDim, SysEqn>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  const MInt coordSize = nDim * maxNoNodesXD();
  const MInt varSize = maxNoNodesXD() * noVars();
  const MInt nodeVarSize = maxNoNodesXD() * noNodeVars();

  // Cell id
  fill_invalid(m_cellId, begin, end);

  // Polynomial degree
  fill_invalid(m_polyDeg, begin, end);

  // Number of nodes 1D
  fill_invalid(m_noNodes1D, begin, end);

  // Surface ids
  fill_invalid(m_surfaceIds, begin, end, noSurfaces());

  // Node coordinates
  fill_invalid(m_nodeCoordinates, begin, end, coordSize);

  // Variables
  fill_invalid(m_variables, begin, end, varSize);

  // Old variables
  fill_invalid(m_timeIntStorage, begin, end, varSize);

  // Node variables
  fill_invalid(m_nodeVariables, begin, end, nodeVarSize);

  // Right hand side
  fill_invalid(m_rightHandSide, begin, end, varSize);

  // External source terms
  fill_invalid(m_externalSource, begin, end, varSize);

  // Inverse jacobian
  fill_invalid(m_invJacobian, begin, end);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim, class SysEqn>
template <class Functor, class T>
void ElementCollector<nDim, SysEqn>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                                    const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  const MInt varSize = maxNoNodesXD() * noVars();

  // Cell id
  copyData(source.m_cellId, m_cellId, c, begin, end, destination);

  // Polynomial degree
  copyData(source.m_polyDeg, m_polyDeg, c, begin, end, destination);

  // Number of nodes 1D
  copyData(source.m_noNodes1D, m_noNodes1D, c, begin, end, destination);

  // Surface ids
  copyData(source.m_surfaceIds, m_surfaceIds, c, begin, end, destination, noSurfaces());

  // Node coordinates
  copyData(source.m_nodeCoordinates, m_nodeCoordinates, c, begin, end, destination, nDim * maxNoNodesXD());

  // Variables
  copyData(source.m_variables, m_variables, c, begin, end, destination, varSize);

  // Old variables
  copyData(source.m_timeIntStorage, m_timeIntStorage, c, begin, end, destination, varSize);

  // Node variables
  copyData(source.m_nodeVariables, m_nodeVariables, c, begin, end, destination, maxNoNodesXD() * noNodeVars());

  // Right hand side
  copyData(source.m_rightHandSide, m_rightHandSide, c, begin, end, destination, varSize);

  // External source terms
  copyData(source.m_externalSource, m_externalSource, c, begin, end, destination, varSize);

  // Inverse jacobian
  copyData(source.m_invJacobian, m_invJacobian, c, begin, end, destination);
}


/// Set maximum polynomial degree
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::maxPolyDeg(const MInt maxPolyDeg_) {
  const MInt oldMaxPolyDeg = m_maxPolyDeg;
  m_maxPolyDeg = maxPolyDeg_;
  return oldMaxPolyDeg;
}

/// Set maximum noNodes
template <MInt nDim, class SysEqn>
void ElementCollector<nDim, SysEqn>::maxNoNodes1D(const MInt maxNoNodes1D_) {
  m_maxNoNodes1D = maxNoNodes1D_;
  m_maxNoNodesXD = ipow(m_maxNoNodes1D, nDim);
}


/// \brief Return element id for a given cell id (or -1 if not found).
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-13
///
/// \param[in] cellId Cell id for which the element is searched.
///
/// \return The element id if found, -1 otherwise.
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::getElementByCellId(const MInt cellId) const {
  return getElementByCellId(0, this->size(), cellId);
}


/// \brief Search for element with cell id `cellId` in [first, last).
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-13
///
/// \param[in] first First element id to consider.
/// \param[in] last Last element id to consider + 1.
/// \param[in] cellId Cell id of the searched element.
///
/// \return The element id if found, -1 otherwise.
template <MInt nDim, class SysEqn>
MInt ElementCollector<nDim, SysEqn>::getElementByCellId(const MInt first, const MInt last, const MInt cellId) const {
  const MInt* const begin = &m_cellId[0];
  const MInt* const end = begin + (last - first);
  const MInt* const low = std::lower_bound(begin, end, cellId);

  // Return not found if searched cell id is smaller than [first, last]
  if(low == begin && *low != cellId) {
    return -1;
  }

  // Return not found if std::lower_bound does not find anything
  if(low == end) {
    return -1;
  }

  // Return not found if the lower bound does not match the searched id (happens
  // if a cell has no corresponding element)
  if(*low != cellId) {
    return -1;
  }

  // Otherwise return found element id
  return first + std::distance(begin, low);
}

} // namespace collector
} // namespace dg
} // namespace maia


// Undefine macros that should not be used outside this file
#undef DGCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_NODE_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_SURFACE_ID_ACCESSOR
#undef ENSURE_CONDITION

#endif // ifndef DGELEMENTCOLLECTOR_H_
