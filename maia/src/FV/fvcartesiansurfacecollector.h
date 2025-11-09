// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIAFVSURFACECOLLECTOR_H_
#define MAIAFVSURFACECOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <limits>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT

// Sanity-checking macros for accessors
#if defined(FVSURFACECOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noVariables(),                                                                                 \
        "variable id = " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noVariables()) + ")", AT_);     \
  } while(false)
#define ENSURE_VALID_FLUX_VARIABLE_ID_ACCESSOR(id)                                                                     \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noFVariables(),                                                              \
                          "flux variable id = " + std::to_string(id) + " out-of-bounds [0, "                           \
                              + std::to_string(noFVariables()) + ")",                                                  \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_DIR_ACCESSOR(id)                                                                                  \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < 2, "dir = " + std::to_string(id) + " out-of-bounds [0, 2)", AT_);            \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        dir >= 0 && dir < nDim,                                                                                        \
        "direction dir = " + std::to_string(dir) + " out-of-bounds [0, " + std::to_string(nDim) + ")", AT_);           \
  } while(false)
#define ENSURE_VALID_SURFACE_COEFFICIENT_ID_ACCESSOR(id)                                                               \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noSurfaceCoefficients(),                                                     \
                          "surface coefficient id = " + std::to_string(id) + " out-of-bounds [0, "                     \
                              + std::to_string(noSurfaceCoefficients()) + ")",                                         \
                          AT_);                                                                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_FLUX_VARIABLE_ID_ACCESSOR(id)                                                                     \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIR_ACCESSOR(id)                                                                                  \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SURFACE_COEFFICIENT_ID_ACCESSOR(id)                                                               \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia::fv::surface_collector {

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


/// Class that represents FV surface collector.
template <MInt nDim>
class FvSurfaceCollector : public maia::container::Container<FvSurfaceCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<FvSurfaceCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<FvSurfaceCollector<nDim>, Invalid>;
  using Base::resetStorage;
  using Base::resizeStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;


 public:
  // Types
  template <class T>
  using Invalid = typename maia::fv::surface_collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr FvSurfaceCollector() = default;
  void checkVariables();

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;
  using Base::resize;

  // Accessors
  MInt& bndryCndId(const MInt id);
  MInt bndryCndId(const MInt id) const;

  MInt& orientation(const MInt id);
  MInt orientation(const MInt id) const;

  MFloat& factor(const MInt id, const MInt varId);
  MFloat factor(const MInt id, const MInt varId) const;

  MFloat& area(const MInt id);
  MFloat area(const MInt id) const;

  MFloat& coordinate(const MInt id, const MInt dir);
  MFloat coordinate(const MInt id, const MInt dir) const;

  MFloat& deltaX(const MInt id, const MInt varId);
  MFloat deltaX(const MInt id, const MInt varId) const;

  MInt& nghbrCellId(const MInt id, const MInt dir);
  MInt nghbrCellId(const MInt id, const MInt dir) const;

  MFloat& variable(const MInt id, const MInt dir, const MInt varId);
  MFloat variable(const MInt id, const MInt dir, const MInt varId) const;

  MFloat& upwindCoefficient(const MInt id);
  MFloat upwindCoefficient(const MInt id) const;

  MFloat& surfaceCoefficient(const MInt id, const MInt dimCoefficient);
  MFloat surfaceCoefficient(const MInt id, const MInt dimCoefficient) const;

  MFloat& flux(const MInt id, const MInt fVarId);
  MFloat flux(const MInt id, const MInt fVarId) const;

  // Allow setting number of PVs
  void setNoSpecies(const MInt noSpecies_);
  void setNoMaxSrfcs(const MInt noMaxSrfcs_);
  void setNoVariables(const MInt noVariables_);
  void setNoFVariables(const MInt noFVariables_);
  void setNoSurfaceCoefficients(const MInt noSurfaceCoefficients_);
  void setSolverType(const MInt solverType);

  /// Return number of species
  constexpr MInt noSpecies() const { return m_noSpecies; }

  /// Return number of maximum surfaces
  constexpr MInt noMaxSrfcs() const { return m_noMaxSrfcs; }

  /// Return number of variables
  constexpr MInt noVariables() const { return m_noVariables; }

  /// Return number of flux variables
  constexpr MInt noFVariables() const { return m_noFVariables; }

  /// Return number of surface coefficients
  constexpr MInt noSurfaceCoefficients() const { return m_noSurfaceCoefficients; }

  /// Return solver type
  constexpr SolverType solverType() const { return (SolverType)m_solverType; }

 private:
  // Methods required by base class for CRTP
  void reset();
  void resize() override;
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  // TODO:FV,toremove
  /// Number of species (not being used)
  MInt m_noSpecies = -1;

  // Number of maximum surfaces
  MInt m_noMaxSrfcs = -1;

  // Number of variables
  MInt m_noVariables = -1;

  // Number of flux variables
  MInt m_noFVariables = -1;

  // Number of surface coefficients
  MInt m_noSurfaceCoefficients = -1;

  // Solver type
  MInt m_solverType = -1;

  // Data containers
  Storage<MInt> m_bndryCndId{};
  Storage<MInt> m_orientation{};
  Storage<MFloat> m_factor{};
  Storage<MFloat> m_area{};
  Storage<MFloat> m_coordinates{};
  Storage<MFloat> m_deltaX{};
  Storage<MInt> m_nghbrCellIds{};
  Storage<MFloat> m_variables{};
  Storage<MFloat> m_upwindCoefficent{};
  Storage<MFloat> m_flux{};
  Storage<MFloat> m_surfaceCoefficients{};
};


/// Print data size variables
template <MInt nDim>
void FvSurfaceCollector<nDim>::checkVariables() {
  std::cerr << "@FvSurfaceCollector: noVariables() = " << noVariables() << ", noFVariables() = " << noFVariables()
            << ", noSurfaceCoefficients() = " << noSurfaceCoefficients() << ", m_noMaxSrfcs = " << noMaxSrfcs()
            << std::endl;
}

// TODO labels:FV replace reset with resize+invalidate?
/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void FvSurfaceCollector<nDim>::reset() {
  resetStorage(1, m_bndryCndId);
  resetStorage(1, m_orientation);
  resetStorage(2, m_factor);
  resetStorage(1, m_area);
  resetStorage(nDim, m_coordinates);
  resetStorage(2 * nDim, m_deltaX);
  resetStorage(2, m_nghbrCellIds);
  resetStorage(2 * noVariables(), m_variables);
  resetStorage(1, m_upwindCoefficent);
  resetStorage(noSurfaceCoefficients(), m_surfaceCoefficients);
  resetStorage(noFVariables(), m_flux);
}


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void FvSurfaceCollector<nDim>::resize() {
  resizeStorage(1, m_bndryCndId);
  resizeStorage(1, m_orientation);
  resizeStorage(2, m_factor);
  resizeStorage(1, m_area);
  resizeStorage(nDim, m_coordinates);
  resizeStorage(2 * nDim, m_deltaX);
  resizeStorage(2, m_nghbrCellIds);
  resizeStorage(2 * noVariables(), m_variables);
  resizeStorage(1, m_upwindCoefficent);
  resizeStorage(noSurfaceCoefficients(), m_surfaceCoefficients);
  resizeStorage(noFVariables(), m_flux);
}


// TODO:FV Accessor function variables and macros might need to be changed later
// TODO:FV id is used as the surfaceId
//         Which must be defined when the collector is created as
//         noBoundarySurfaces + noInnerSurfaces. However, since this information
//         is not known, the number of cells can be used to estimate the collector size
/// Accessor for bndryCndId.
template <MInt nDim>
MInt& FvSurfaceCollector<nDim>::bndryCndId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCndId[id];
}
/// Accessor for bndryCndId (const version).
template <MInt nDim>
MInt FvSurfaceCollector<nDim>::bndryCndId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCndId[id];
}

/// Accessor for orientation.
template <MInt nDim>
MInt& FvSurfaceCollector<nDim>::orientation(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_orientation[id];
}
/// Accessor for orientation (const version).
template <MInt nDim>
MInt FvSurfaceCollector<nDim>::orientation(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_orientation[id];
}

/// Accessor for factor.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::factor(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_factor[id * 2 + varId];
}
/// Accessor for factor (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::factor(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_factor[id * 2 + varId];
}

/// Accessor for area.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::area(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_area[id];
}
/// Accessor for area (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::area(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_area[id];
}

/// Accessor for coordinate.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::coordinate(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}
/// Accessor for coordinate (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::coordinate(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}

/// Accessor for deltaX.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::deltaX(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_deltaX[id * 2 * nDim + varId];
}
/// Accessor for deltaX (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::deltaX(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_deltaX[id * 2 * nDim + varId];
}

/// Accessor for nghbrCellId.
template <MInt nDim>
MInt& FvSurfaceCollector<nDim>::nghbrCellId(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_nghbrCellIds[id * 2 + dir];
}
/// Accessor for nghbrCellId (const version).
template <MInt nDim>
MInt FvSurfaceCollector<nDim>::nghbrCellId(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_nghbrCellIds[id * 2 + dir];
}

/// Accessor for variable.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::variable(const MInt id, const MInt dir, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIR_ACCESSOR(dir);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_variables[(id * 2 + dir) * noVariables() + varId];
}
/// Accessor for variable (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::variable(const MInt id, const MInt dir, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIR_ACCESSOR(dir);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_variables[(id * 2 + dir) * noVariables() + varId];
}

/// Accessor for upwind coefficient.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::upwindCoefficient(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_upwindCoefficent[id];
}
/// Accessor for upwind coefficient (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::upwindCoefficient(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_upwindCoefficent[id];
}

/// Accessor for surfaceCoefficient.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::surfaceCoefficient(const MInt id, const MInt dimCoefficient) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SURFACE_COEFFICIENT_ID_ACCESSOR(dimCoefficient);
  return m_surfaceCoefficients[id * noSurfaceCoefficients() + dimCoefficient];
}
/// Accessor for surfaceCoefficient (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::surfaceCoefficient(const MInt id, const MInt dimCoefficient) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SURFACE_COEFFICIENT_ID_ACCESSOR(dimCoefficient);
  return m_surfaceCoefficients[id * noSurfaceCoefficients() + dimCoefficient];
}

/// Accessor for flux.
template <MInt nDim>
MFloat& FvSurfaceCollector<nDim>::flux(const MInt id, const MInt fVarId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_FLUX_VARIABLE_ID_ACCESSOR(varId);
  return m_flux[id * noFVariables() + fVarId];
}
/// Accessor for flux (const version).
template <MInt nDim>
MFloat FvSurfaceCollector<nDim>::flux(const MInt id, const MInt fVarId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_FLUX_VARIABLE_ID_ACCESSOR(varId);
  return m_flux[id * noFVariables() + fVarId];
}
/// Set number of species and update number of variables
template <MInt nDim>
void FvSurfaceCollector<nDim>::setNoSpecies(const MInt noSpecies_) {
  m_noSpecies = noSpecies_;
}

/// Set number of maximum surfaces and update number of variables
template <MInt nDim>
void FvSurfaceCollector<nDim>::setNoMaxSrfcs(const MInt noMaxSrfcs_) {
  m_noMaxSrfcs = noMaxSrfcs_;
}

/// Set number of variables and update number of variables
template <MInt nDim>
void FvSurfaceCollector<nDim>::setNoVariables(const MInt noVariables_) {
  m_noVariables = noVariables_;
}

/// Set number of flux variables
template <MInt nDim>
void FvSurfaceCollector<nDim>::setNoFVariables(const MInt noFVariables_) {
  m_noFVariables = noFVariables_;
}

/// Set number of surface coefficients and update number of surface coefficients
template <MInt nDim>
void FvSurfaceCollector<nDim>::setNoSurfaceCoefficients(const MInt noSurfaceCoefficients_) {
  m_noSurfaceCoefficients = noSurfaceCoefficients_;
}

/// Set solver type
template <MInt nDim>
void FvSurfaceCollector<nDim>::setSolverType(const MInt solverType) {
  m_solverType = solverType;
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void FvSurfaceCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // bndryCndId
  fill_invalid(m_bndryCndId, begin, end, 1, -1);
  // orientation
  fill_invalid(m_orientation, begin, end, 1, -1);
  // factors
  fill_invalid(m_factor, begin, end, 2);
  // area
  fill_invalid(m_area, begin, end);
  // coordinates
  fill_invalid(m_coordinates, begin, end, nDim);
  // deltaX
  fill_invalid(m_deltaX, begin, end, 2 * nDim);
  // nghbrCellIds
  fill_invalid(m_nghbrCellIds, begin, end, 2, -1);
  // variables
  fill_invalid(m_variables, begin, end, 2 * noVariables());
  // upwindCoefficient
  fill_invalid(m_upwindCoefficent, begin, end);
  // surfaceCoefficients
  fill_invalid(m_surfaceCoefficients, begin, end, noSurfaceCoefficients());
  // flux
  fill_invalid(m_flux, begin, end, noFVariables());
}

/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void FvSurfaceCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                              const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVSURFACECOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // bndryCndId
  copyData(source.m_bndryCndId, m_bndryCndId, c, begin, end, destination);
  // orientation
  copyData(source.m_orientation, m_orientation, c, begin, end, destination);
  // factors
  copyData(source.m_factor, m_factor, c, begin, end, destination, 2);
  // area
  copyData(source.m_area, m_area, c, begin, end, destination);
  // coordinates
  copyData(source.m_coordinates, m_coordinates, c, begin, end, destination, nDim);
  // deltaX
  copyData(source.m_deltaX, m_deltaX, c, begin, end, destination, 2 * nDim);
  // nghbrCellIds
  copyData(source.m_nghbrCellIds, m_nghbrCellIds, c, begin, end, destination, 2);
  // variables
  copyData(source.m_variables, m_variables, c, begin, end, destination, 2 * noVariables());
  // upwindCoefficient
  copyData(source.m_upwindCoefficent, m_upwindCoefficent, c, begin, end, destination);
  // surfaceCoefficients
  copyData(source.m_surfaceCoefficients, m_surfaceCoefficients, c, begin, end, destination, noSurfaceCoefficients());
  // flux
  copyData(source.m_flux, m_flux, c, begin, end, destination, noFVariables());
}
} // namespace maia::fv::surface_collector


// Undefine macros that should not be used outside this file
#undef FVSURFACECOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_FLUX_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_DIR_ACCESSOR
#undef ENSURE_VALID_COORDINATE_DIR_ACCESSOR
#undef ENSURE_VALID_SURFACE_COEFFICIENT_ID_ACCESSOR

#endif // ifndef MAIAFVSURFACECOLLECTOR_H_
