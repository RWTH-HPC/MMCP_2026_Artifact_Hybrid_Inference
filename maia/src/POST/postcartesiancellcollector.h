// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTCELLCOLLECTOR_H_
#define POSTCELLCOLLECTOR_H_

#include "postcellproperties.h"
#include "property.h"

/// Namespace for auxiliary functions/classes
namespace maia {
namespace post {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::post::cell::BitsetType;

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

/// Class that represents Post cell collector.
template <MInt nDim>
class PostCellCollector : public maia::container::Container<PostCellCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<PostCellCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<PostCellCollector<nDim>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;
  using BitsetType = maia::post::cell::BitsetType;

 public:
  // Types
  template <class T>
  using Invalid = typename maia::post::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr PostCellCollector() = default;

  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  MFloat& variable(const MInt id, const MInt dim);
  MFloat variable(const MInt id, const MInt dim) const;

  /// Return number of variables
  constexpr MInt noVariables() const { return m_noVariables; }

  // Property-related accessors
  BitsetType::reference hasProperty(const MInt id, const PostCell p);
  MBool hasProperty(const MInt id, const PostCell p) const;
  void resetProperties(const MInt id);
  BitsetType& properties(const MInt id);

  // Allow setting number of species and rans variables
  void setNoVariables(const MInt noVars_);

 private:
  /// Number of variables
  MInt m_noVariables = 0;

  // Data containers
  Storage<MFloat> m_Variables{};
  Storage<BitsetType> m_properties{};

  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);
};


template <MInt nDim>
void PostCellCollector<nDim>::setNoVariables(const MInt noVars_) {
  m_noVariables = noVars_;
}

/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void PostCellCollector<nDim>::reset() {
  resetStorage(noVariables(), m_Variables);
  resetStorage(1, m_properties);
}

/// Accessor for variables.
template <MInt nDim>
MFloat& PostCellCollector<nDim>::variable(const MInt id, const MInt varId) {
  return m_Variables[id * noVariables() + varId];
}

/// Accessor for variables (const version).
template <MInt nDim>
MFloat PostCellCollector<nDim>::variable(const MInt id, const MInt varId) const {
  return m_Variables[id * noVariables() + varId];
}

/// Accessor for properties.
template <MInt nDim>
BitsetType::reference PostCellCollector<nDim>::hasProperty(const MInt id, const PostCell p) {
  return m_properties[id][maia::post::cell::p(p)];
}

/// Accessor for properties (const version).
template <MInt nDim>
MBool PostCellCollector<nDim>::hasProperty(const MInt id, const PostCell p) const {
  return m_properties[id][maia::post::cell::p(p)];
}

/// Reset all properties.
template <MInt nDim>
void PostCellCollector<nDim>::resetProperties(const MInt id) {
  m_properties.at(id).reset();
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void PostCellCollector<nDim>::invalidate(const MInt begin, const MInt end) {
  // variables
  fill_invalid(m_Variables, begin, end, noVariables());

  // Properties
  fill_invalid(m_properties, begin, end);
}

/// Accessor for properties.
template <MInt nDim>
BitsetType& PostCellCollector<nDim>::properties(const MInt id) {
  // ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}

/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void PostCellCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                             const MInt destination) {
  // variables
  copyData(source.m_Variables, m_Variables, c, begin, end, destination, noVariables());
  // Properties
  copyData(source.m_properties, m_properties, c, begin, end, destination);
}


} // namespace collector
} // namespace post
} // namespace maia

#endif // ifndef POSTCELLCOLLECTOR_H_
