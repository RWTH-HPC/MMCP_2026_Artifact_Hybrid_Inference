// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ACASURFACEDATACOLLECTOR_H_
#define ACASURFACEDATACOLLECTOR_H_

#include <algorithm>
#include <numeric>
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
// #define ACACOLLECTOR_SOA_MEMORY_LAYOUT

// The macro 'ACACOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks
// for all accessors. It is enabled for build type "extra_debug".
//#ifdef MAIA_EXTRA_DEBUG
// NOTE: enable checks in normal debug mode until everything is working correctly!
#ifndef NDEBUG
#define ACACOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(ACACOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
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
#define ENSURE_VALID_SAMPLE_ID_ACCESSOR(id)                                                                            \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noSamples(),                                                                                   \
        "sample id " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noSamples()) + ")", AT_);           \
  } while(false)
#define ENSURE_VALID_DIR_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < nDim,                                                                        \
                          "direction id = " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(nDim) + ")", \
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
#define ENSURE_VALID_SAMPLE_ID_ACCESSOR(id)                                                                            \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIR_ID_ACCESSOR(id)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace acoustic_analogy {
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
template <MInt nDim>
class SurfaceDataCollector : public maia::container::Container<SurfaceDataCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<SurfaceDataCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<SurfaceDataCollector<nDim>, Invalid>;
  using Base::resetStorage;
  using Base::resizeStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  // Types
  template <class T>
  using Invalid = typename maia::acoustic_analogy::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr SurfaceDataCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;
  using Base::resize;

  // Accessors
  MFloat& surfaceCoordinates(const MInt id);
  MFloat surfaceCoordinates(const MInt id) const;

  MFloat& surfaceArea(const MInt id);
  MFloat surfaceArea(const MInt id) const;

  MFloat& surfaceNormal(const MInt id);
  MFloat surfaceNormal(const MInt id) const;

  MFloat& variables(const MInt id, const MInt var);
  MFloat& variables(const MInt id) { return variables(id, 0); };

  MFloat& variables(const MInt id, const MInt var, const MInt sample);
  MFloat variables(const MInt id, const MInt var, const MInt sample) const;

  MFloat& complexVariables(const MInt id, const MInt var);
  MFloat& complexVariables(const MInt id) { return complexVariables(id, 0); };

  MFloat& complexVariables(const MInt id, const MInt var, const MInt sample, const MInt component);
  MFloat complexVariables(const MInt id, const MInt var, const MInt sample, const MInt component) const;

  // Auxiliary methods

  /// Return number of variables
  MInt noVars() { return m_noVars; }
  /// Set number of variables
  void setNoVars(const MInt noVars_) {
    ENSURE_CONDITION(noVars_ > 0 && noVars_ <= s_maxNoVars, "Invalid number of variables.");
    m_noVars = noVars_;
  };

  /// Return number of complex variables
  MInt noComplexVars() { return m_noComplexVars; }
  /// Set number of complex variables
  void setNoComplexVars(const MInt noComplexVars_) {
    ENSURE_CONDITION(noComplexVars_ > 0 && noComplexVars_ <= s_maxNoComplexVars,
                     "Invalid number of complex variables.");
    m_noComplexVars = noComplexVars_;
  };

  /// Return number of samples (i.e. number of time steps)
  MInt noSamples() { return m_noSamples; }
  /// Set number of samples
  void setNoSamples(const MInt noSamples_) {
    ENSURE_CONDITION(noSamples_ > 0 && noSamples_ <= s_maxNoSamples, "Invalid number of samples.");
    m_noSamples = noSamples_;
  };

 private:
  // Methods required by base class for CRTP
  void reset();
  void resize() override;
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

 private:
  // Member variables

  // Number of variables
  MInt m_noVars = -1;

  // Number of complex variables
  MInt m_noComplexVars = -1;

  // Number of samples
  MInt m_noSamples = -1;

  // Data containers
  Storage<MFloat> m_surfaceCoordinates{};
  Storage<MFloat> m_surfaceArea{};
  Storage<MFloat> m_surfaceNormal{};
  Storage<MFloat> m_variables{};
  Storage<MFloat> m_complexVariables{};

  /// Maximum number of variables (in the time domain)
  static constexpr MInt s_maxNoVars = 14;
  /// Maximum number of complex variables (in the frequency domain)
  static constexpr MInt s_maxNoComplexVars = 5;
  /// Maximum number of samples
  static constexpr MInt s_maxNoSamples = 1048576; // 2^20
};


/// Reset, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void SurfaceDataCollector<nDim>::reset() {
  const MInt varSize = noVars() * noSamples();
  const MInt complexVarSize = 2 * noComplexVars() * noSamples();

  resetStorage(nDim, m_surfaceCoordinates);
  resetStorage(1, m_surfaceArea);
  resetStorage(nDim, m_surfaceNormal);
  resetStorage(varSize, m_variables);
  resetStorage(complexVarSize, m_complexVariables);
}

/// Resize data strucutres reusing values from previous state
template <MInt nDim>
void SurfaceDataCollector<nDim>::resize() {
  const MInt varSize = noVars() * noSamples();
  const MInt complexVarSize = 2 * noComplexVars() * noSamples();

  resizeStorage(nDim, m_surfaceCoordinates);
  resizeStorage(1, m_surfaceArea);
  resizeStorage(nDim, m_surfaceNormal);
  resizeStorage(varSize, m_variables);
  resizeStorage(complexVarSize, m_complexVariables);
}


/// Accessor for surface element coordinates.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::surfaceCoordinates(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceCoordinates[id * nDim];
}
/// Accessor for surface element coordinates (const version).
template <MInt nDim>
MFloat SurfaceDataCollector<nDim>::surfaceCoordinates(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceCoordinates[id * nDim];
}


/// Accessor for surface area (or segment length in 2D).
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::surfaceArea(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceArea[id];
}
/// Accessor for surface area (or segment length in 2D) (const version).
template <MInt nDim>
MFloat SurfaceDataCollector<nDim>::surfaceArea(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceCoordinates[id];
}


/// Accessor for surface normal.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::surfaceNormal(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceNormal[id * nDim];
}
/// Accessor for surface normal (const version).
template <MInt nDim>
MFloat SurfaceDataCollector<nDim>::surfaceNormal(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_surfaceNormal[id * nDim];
}


/// Accessor for variables.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::variables(const MInt id, const MInt var) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  return m_variables[id * noVars() * noSamples() + var * noSamples()];
}


/// NOTE: the following variables accessors might be very inefficient to use in a loop e.g. over all
/// samples and should only be used during the development of the method.
/// Accessor for variables.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::variables(const MInt id, const MInt var, const MInt sample) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  ENSURE_VALID_SAMPLE_ID_ACCESSOR(sample);
  return m_variables[id * noVars() * noSamples() + var * noSamples() + sample];
}
/// Accessor for variables (const version).
template <MInt nDim>
MFloat SurfaceDataCollector<nDim>::variables(const MInt id, const MInt var, const MInt sample) const {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  ENSURE_VALID_SAMPLE_ID_ACCESSOR(sample);
  return m_variables[id * noVars() * noSamples() + var * noSamples() + sample];
}


/// Accessor for complex variables.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::complexVariables(const MInt id, const MInt var) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  return m_complexVariables[id * 2 * noComplexVars() * noSamples() + var * 2 * noSamples()];
}


/// NOTE: the following complexariables accessors might be very inefficient to use in a loop e.g.
/// over all samples and should only be used during the development of the method.
/// Accessor for complex variables.
template <MInt nDim>
MFloat& SurfaceDataCollector<nDim>::complexVariables(const MInt id, const MInt var, const MInt sample,
                                                     const MInt component) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  ENSURE_VALID_SAMPLE_ID_ACCESSOR(sample);
  return m_complexVariables[id * 2 * noComplexVars() * noSamples() + var * 2 * noSamples() + 2 * sample + component];
}
/// Accessor for complex variables (const version).
template <MInt nDim>
MFloat SurfaceDataCollector<nDim>::complexVariables(const MInt id, const MInt var, const MInt sample,
                                                    const MInt component) const {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(var);
  ENSURE_VALID_SAMPLE_ID_ACCESSOR(sample);
  return m_complexVariables[id * 2 * noComplexVars() * noSamples() + var * 2 * noSamples() + 2 * sample + component];
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void SurfaceDataCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  const MInt varSize = noVars() * noSamples();
  const MInt complexVarSize = 2 * noComplexVars() * noSamples();

  // Surface coordinates
  fill_invalid(m_surfaceCoordinates, begin, end, nDim);

  // Surface area
  fill_invalid(m_surfaceArea, begin, end);

  // Surface normal
  fill_invalid(m_surfaceNormal, begin, end, nDim);

  // Variables
  fill_invalid(m_variables, begin, end, varSize);

  // Complex variables
  fill_invalid(m_complexVariables, begin, end, complexVarSize);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void SurfaceDataCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                                const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef ACACOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  const MInt varSize = noVars() * noSamples();
  const MInt complexVarSize = 2 * noComplexVars() * noSamples();

  // Surface coordinates
  copyData(source.m_surfaceCoordinates, m_surfaceCoordinates, c, begin, end, destination, nDim);

  // Surface area
  copyData(source.m_surfaceArea, m_surfaceArea, c, begin, end, destination);

  // Surface normal
  copyData(source.m_surfaceNormal, m_surfaceNormal, c, begin, end, destination, nDim);

  // Variables
  copyData(source.m_variables, m_variables, c, begin, end, destination, varSize);

  // Variables
  copyData(source.m_complexVariables, m_complexVariables, c, begin, end, destination, complexVarSize);
}


} // namespace collector
} // namespace acoustic_analogy
} // namespace maia


// Undefine macros that should not be used outside this file
#undef ACACOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_SAMPLE_ID_ACCESSOR
#undef ENSURE_VALID_DIR_ID_ACCESSOR
#undef ENSURE_CONDITION

#endif // ifndef DGELEMENTCOLLECTOR_H_
