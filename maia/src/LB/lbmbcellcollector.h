// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBMBCOLLECTOR_H_
#define LBMBCOLLECTOR_H_

#include <vector>

#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "compiler_config.h"
#include "lbcellcollector.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define LBMBCOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'LBCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define LBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(LBCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
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
#define ENSURE_VALID_DISTANCE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noDistances(),                                                                                 \
        "distance id = " + std::to_string(id) + " is out-of-bounds [0, " + std::to_string(noDistances()), AT_);        \
  } while(false)
#define ENSURE_VALID_DIM_ID_ACCESSOR(dim)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(dim >= 0 && dim < nDim,                                                                      \
                          "dim = " + std::to_string(dim) + " is out-of-bounds [0, " + std::to_string(nDim) + ")",      \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != LbCell::NumProperties, "Invalid property", AT_);                                        \
  } while(false)

#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DISTANCE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIM_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


// Namespace for auxiliary functions/classes
namespace maia::lb::collector {

// Type traits for invalid values. These values are used to initialize/erase nodes.
// Already defined in LbCellCollector.
// template <class T> using Invalid = maia::lb::collector::Invalid<T>;

/// Class that represents LB cell collector.
template <MInt nDim>
class LbMbCellCollector : public maia::container::Container<LbMbCellCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<LbMbCellCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<LbMbCellCollector<nDim>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  // Types
  // template <class T> using Invalid = typename maia::lb::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr LbMbCellCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  // TODO labels:LB Access distance via bndryCellId, setId and distanceId
  // (Or find another way to distinquish between boundaries of different sets ...)

  MInt& cellId(const MInt id);
  MInt cellId(const MInt id) const;
  MFloat& distance(const MInt id, const MInt did);
  MFloat distance(const MInt id, const MInt did) const;
  MFloat& velocity(const MInt id, const MInt dim);
  MFloat velocity(const MInt id, const MInt dim) const;
  MFloat& force(const MInt id, const MInt did, const MInt dim);
  MFloat force(const MInt id, const MInt did, const MInt dim) const;
  MFloat& surfaceCenter(const MInt id, const MInt did, const MInt dim);
  MFloat surfaceCenter(const MInt id, const MInt did, const MInt dim) const;
  MFloat& cellCenter(const MInt id, const MInt dim);
  MFloat cellCenter(const MInt id, const MInt dim) const;
  MFloat& normal(const MInt id, const MInt dim);
  MFloat normal(const MInt id, const MInt dim) const;
  MFloat& density(const MInt id);
  MFloat density(const MInt id) const;

  // Setter
  void setVelocity(const MInt id, const std::array<MFloat, nDim>& velocity);

  // Adder
  void addVelocity(const MInt id, const std::array<MFloat, nDim>& velocity);

  /// Number of distances
  constexpr MInt noDistances() const { return m_noDistances; }
  void noDistances(const MInt noDistances) { m_noDistances = noDistances; }

  /// Number of forces
  constexpr MInt noForces() const { return m_noDistances; }
  void noForces(const MInt noDistances) { m_noDistances = noDistances; }

  /// Number of dimensions
  constexpr MInt noDimensions() const { return nDim; }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  // Properties
  MInt m_noDistances = 0;

  // Data containers
  Storage<MInt> m_cellId{};
  Storage<MFloat> m_distances{};
  Storage<MFloat> m_velocities{};
  Storage<MFloat> m_forces{};
  Storage<MFloat> m_surfaceCenters{};
  Storage<MFloat> m_cellCenters{};
  Storage<MFloat> m_normal{};
  Storage<MFloat> m_density{};
};

/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void LbMbCellCollector<nDim>::reset() {
  resetStorage(1, m_cellId);
  resetStorage(noDistances(), m_distances);
  resetStorage(nDim, m_velocities);
  resetStorage(nDim * noDistances(), m_forces);
  resetStorage(nDim * noDistances(), m_surfaceCenters);
  resetStorage(nDim, m_cellCenters);
  resetStorage(nDim, m_normal);
  resetStorage(1, m_density);
}

//---------------------------------ACCESSORS------------------------------------

/// Accessor for cellId.
template <MInt nDim>
MInt& LbMbCellCollector<nDim>::cellId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellId[id];
}
/// Accessor for cellId (const version).
template <MInt nDim>
MInt LbMbCellCollector<nDim>::cellId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellId[id];
}

/// Accessor for distances.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::distance(const MInt id, const MInt did) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_distances[id * noDistances() + did];
}
/// Accessor for distances (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::distance(const MInt id, const MInt did) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_distances[id * noDistances() + did];
}

/// Accessor for velocities.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::velocity(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_velocities[id * nDim + dim];
}
/// Accessor for velocities (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::velocity(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_velocities[id * nDim + dim];
}
// Setter for velocities
template <MInt nDim>
void LbMbCellCollector<nDim>::setVelocity(const MInt id, const std::array<MFloat, nDim>& velocity) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  for(MInt n = 0; n < nDim; n++) {
    m_velocities[id * nDim + n] = velocity[n];
  }
}
// Adder for velocities
template <MInt nDim>
void LbMbCellCollector<nDim>::addVelocity(const MInt id, const std::array<MFloat, nDim>& velocity) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  for(MInt n = 0; n < nDim; n++) {
    m_velocities[id * nDim + n] += velocity[n];
  }
}

/// Accessor for forces.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::force(const MInt id, const MInt did, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_forces[id * nDim * noDistances() + nDim * did + dim];
}
/// Accessor for forces (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::force(const MInt id, const MInt did, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_forces[id * nDim * noDistances() + nDim * did + dim];
}

/// Accessor for surfaceCenter.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::surfaceCenter(const MInt id, const MInt did, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_surfaceCenters[id * nDim * noDistances() + nDim * did + dim];
}
/// Accessor for surfaceCenter (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::surfaceCenter(const MInt id, const MInt did, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  ENSURE_VALID_DISTANCE_ID_ACCESSOR(did);
  return m_surfaceCenters[id * nDim * noDistances() + nDim * did + dim];
}

/// Accessor for cellCenter.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::cellCenter(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_cellCenters[id * nDim + dim];
}
/// Accessor for cellCenter (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::cellCenter(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_cellCenters[id * nDim + dim];
}

/// Accessor for surface normal.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::normal(const MInt id, const MInt dim) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_normal[id * nDim + dim];
}
/// Accessor for surface normal (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::normal(const MInt id, const MInt dim) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIM_ID_ACCESSOR(dim);
  return m_normal[id * nDim + dim];
}

/// Accessor for density.
template <MInt nDim>
MFloat& LbMbCellCollector<nDim>::density(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_density[id];
}
/// Accessor for density (const version).
template <MInt nDim>
MFloat LbMbCellCollector<nDim>::density(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_density[id];
}


//---------------------------CRTP COLLECTOR STUFF------------------------------

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void LbMbCellCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  fill_invalid(m_cellId, begin, end);
  fill_invalid(m_distances, begin, end, noDistances());
  fill_invalid(m_velocities, begin, end, nDim);
  fill_invalid(m_forces, begin, end, nDim * noDistances(), 0.0);
  fill_invalid(m_surfaceCenters, begin, end, nDim * noDistances());
  fill_invalid(m_cellCenters, begin, end, nDim);
  fill_invalid(m_normal, begin, end, nDim);
  fill_invalid(m_density, begin, end);
}

/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void LbMbCellCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                             const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef LBMBCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  copyData(source.m_cellId, m_cellId, c, begin, end, destination, 1);
  copyData(source.m_distances, m_distances, c, begin, end, destination, noDistances());
  coptData(source.m_velocities, m_velocities, c, begin, end, destination, nDim);
  copyData(source.m_forces, m_forces, c, begin, end, destination, nDim * noDistances());
  copyData(source.m_surfaceCenters, m_surfaceCenters, c, begin, end, destination, nDim * noDistances());
  copyData(source.m_cellCenters, m_cellCenters, c, begin, end, destination, nDim);
  copyData(source.m_normal, m_normal, c, begin, end, destination, nDim);
  copyData(source.m_density, m_density, c, begin, end, destination, 1);
}

} // namespace maia::lb::collector


// Undefine macros that should not be used outside this file
#undef LBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_DISTANCES_ID_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef LBMBCOLLECTOR_H_
