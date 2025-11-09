// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSPONGELEMENTCOLLECTOR_H_
#define DGSPONGELEMENTCOLLECTOR_H_


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
#define ENSURE_VALID_NODE_ID_ACCESSOR(pos)                                                                             \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(pos >= 0 && pos < noNodesXD(),                                                               \
                          "node id " + std::to_string(pos) + " out-of-bounds [0, " + std::to_string(noNodesXD()),      \
                          AT_);                                                                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_NODE_ID_ACCESSOR(pos)                                                                             \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace dg {
namespace collector {

/// Class that represents DG element collector.
template <MInt nDim, class SysEqn>
class SpongeElementCollector : public maia::container::Container<SpongeElementCollector<nDim, SysEqn>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<SpongeElementCollector<nDim, SysEqn>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<SpongeElementCollector<nDim, SysEqn>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

  /// Return number of nodes in 2D/3D
  constexpr MInt noNodesXD() const { return m_noNodesXD; }

 public:
  // Types
  template <class T>
  using Invalid = typename maia::dg::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr SpongeElementCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  MInt& elementId(const MInt id);
  MInt elementId(const MInt id) const;
  MFloat& spongeEta(const MInt id);
  MFloat spongeEta(const MInt id, const MInt pos) const;
  MFloat spongeEta(const MInt id) const { return spongeEta(id, 0); };

  /// Return maximum polynomial degree
  MInt maxPolyDeg() const { return m_maxPolyDeg; }
  // Set maximum polynomial degree
  MInt maxPolyDeg(const MInt maxPolyDeg_);

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  /// Maximum polynomial degree
  MInt m_maxPolyDeg = -1;

  /// Number of nodes in 2D/3D
  MInt m_noNodesXD = -1;

  // Data containers
  Storage<MInt> m_elementId{};
  Storage<MFloat> m_spongeEta{};
};


/// Reset SpongeElementCollector, re-create data structures.
template <MInt nDim, class SysEqn>
void SpongeElementCollector<nDim, SysEqn>::reset() {
  resetStorage(1, m_elementId);
  resetStorage(noNodesXD(), m_spongeEta);
}


/// Accessor for element id.
template <MInt nDim, class SysEqn>
MInt& SpongeElementCollector<nDim, SysEqn>::elementId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_elementId[id];
}
/// Accessor for element id (const version).
template <MInt nDim, class SysEqn>
MInt SpongeElementCollector<nDim, SysEqn>::elementId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_elementId[id];
}


/// Accessor for sponge eta.
template <MInt nDim, class SysEqn>
MFloat& SpongeElementCollector<nDim, SysEqn>::spongeEta(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeEta[id * noNodesXD()];
}
/// Accessor for sponge eta (const version).
template <MInt nDim, class SysEqn>
MFloat SpongeElementCollector<nDim, SysEqn>::spongeEta(const MInt id, const MInt pos) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_NODE_ID_ACCESSOR(pos);
  return m_spongeEta[id * noNodesXD() + pos];
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim, class SysEqn>
void SpongeElementCollector<nDim, SysEqn>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // Element id
  fill_invalid(m_elementId, begin, end);

  // Sponge eta
  fill_invalid(m_spongeEta, begin, end, noNodesXD());
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim, class SysEqn>
template <class Functor, class T>
void SpongeElementCollector<nDim, SysEqn>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin,
                                                          const MInt end, const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // Element id
  copyData(source.m_elementId, m_elementId, c, begin, end, destination);

  // Sponge eta
  copyData(source.m_spongeEta, m_spongeEta, c, begin, end, destination, noNodesXD());
}


/// Set maximum polynomial degree and update number of nodes in 2D/3D
template <MInt nDim, class SysEqn>
MInt SpongeElementCollector<nDim, SysEqn>::maxPolyDeg(const MInt maxPolyDeg_) {
  const MInt oldMaxPolyDeg = m_maxPolyDeg;
  m_maxPolyDeg = maxPolyDeg_;
  m_noNodesXD = ipow(m_maxPolyDeg + 1, nDim);
  return oldMaxPolyDeg;
}

} // namespace collector
} // namespace dg
} // namespace maia


// Undefine macros that should not be used outside this file
#undef DGCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_NODE_ID_ACCESSOR

#endif // ifndef DGSPONGELEMENTCOLLECTOR_H_
