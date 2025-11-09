// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGHELEMENTCOLLECTOR_H_
#define DGHELEMENTCOLLECTOR_H_


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
#define ENSURE_VALID_SURFACE_DIR_ACCESSOR(dir)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        dir >= 0 && dir < noDirs(),                                                                                    \
        "surface direction " + std::to_string(dir) + " out-of-bounds [0, " + std::to_string(noDirs()) + ")", AT_);     \
  } while(false)
#define ENSURE_VALID_SURFACE_POS_ACCESSOR(pos)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(pos >= 0 && pos < noHrefSurfaces(),                                                          \
                          "surface position out-of-bounds [0, " + std::to_string(noHrefSurfaces()) + ")", AT_);        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SURFACE_DIR_ACCESSOR(dir)                                                                         \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SURFACE_POS_ACCESSOR(pos)                                                                         \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace dg {
namespace collector {

/// Class that represents DG element collector.
template <MInt nDim, class SysEqn>
class HElementCollector : public maia::container::Container<HElementCollector<nDim, SysEqn>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<HElementCollector<nDim, SysEqn>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<HElementCollector<nDim, SysEqn>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;


 public:
  // Types
  template <class T>
  using Invalid = typename maia::dg::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr HElementCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  MInt& elementId(const MInt id);
  MInt elementId(const MInt id) const;
  MInt& hrefSurfaceIds(const MInt id, const MInt dir, const MInt pos);
  MInt hrefSurfaceIds(const MInt id, const MInt dir, const MInt pos) const;

  /// Return number of h-refined surfaces
  static constexpr MInt noHrefSurfaces() { return 2 * (nDim - 1); }

  /// Return number of directions
  static constexpr MInt noDirs() { return 2 * nDim; }

  // Return number of h-refined surface ids
  static constexpr MInt noHrefSurfaceIds() { return noHrefSurfaces() * noDirs(); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  // Data containers
  Storage<MInt> m_elementId{};
  Storage<MInt> m_hrefSurfaceIds{};
};


/// Reset HElementCollector, re-create data structures.
template <MInt nDim, class SysEqn>
void HElementCollector<nDim, SysEqn>::reset() {
  resetStorage(1, m_elementId);
  resetStorage(noHrefSurfaceIds(), m_hrefSurfaceIds);
}


/// Accessor for element id.
template <MInt nDim, class SysEqn>
MInt& HElementCollector<nDim, SysEqn>::elementId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_elementId[id];
}
/// Accessor for element id (const version).
template <MInt nDim, class SysEqn>
MInt HElementCollector<nDim, SysEqn>::elementId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_elementId[id];
}


/// Accessor for h-refined surface ids.
template <MInt nDim, class SysEqn>
MInt& HElementCollector<nDim, SysEqn>::hrefSurfaceIds(const MInt id, const MInt dir, const MInt pos) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SURFACE_DIR_ACCESSOR(dir);
  ENSURE_VALID_SURFACE_POS_ACCESSOR(pos);
  return m_hrefSurfaceIds[id * noHrefSurfaceIds() + dir * noHrefSurfaces() + pos];
}
/// Accessor for h-refined surface ids (const version).
template <MInt nDim, class SysEqn>
MInt HElementCollector<nDim, SysEqn>::hrefSurfaceIds(const MInt id, const MInt dir, const MInt pos) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id * noHrefSurfaceIds());
  ENSURE_VALID_SURFACE_DIR_ACCESSOR(dir);
  ENSURE_VALID_SURFACE_POS_ACCESSOR(pos);
  return m_hrefSurfaceIds[id * noHrefSurfaceIds() + dir * noHrefSurfaces() + pos];
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim, class SysEqn>
void HElementCollector<nDim, SysEqn>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // Element id
  fill_invalid(m_elementId, begin, end);

  // H-refined surface ids
  fill_invalid(m_hrefSurfaceIds, begin, end, noHrefSurfaceIds());
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim, class SysEqn>
template <class Functor, class T>
void HElementCollector<nDim, SysEqn>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                                     const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // Element id
  copyData(source.m_elementId, m_elementId, c, begin, end, destination);

  // H-refined surface ids
  copyData(source.m_hrefSurfaceIds, m_hrefSurfaceIds, c, begin, end, destination, noHrefSurfaceIds());
}

} // namespace collector
} // namespace dg
} // namespace maia


// Undefine macros that should not be used outside this file
#undef DGCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_SURFACE_DIR_ACCESSOR
#undef ENSURE_VALID_SURFACE_POS_ACCESSOR

#endif // ifndef DGHELEMENTCOLLECTOR_H_
