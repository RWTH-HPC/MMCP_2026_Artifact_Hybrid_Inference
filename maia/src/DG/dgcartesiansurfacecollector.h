// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_DGSURFACECOLLECTOR_H
#define MAIA_DGSURFACECOLLECTOR_H

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
#define ENSURE_VALID_SIDE_ACCESSOR(side)                                                                               \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(side >= 0 && side < 2, "side " + std::to_string(side) + " out-of-bounds [0, 2)", AT_);       \
  } while(false)
#define ENSURE_VALID_DIR_ACCESSOR(orientation)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        orientation >= 0 && orientation < nDim,                                                                        \
        "orientation " + std::to_string(orientation) + " out-of-bounds [0, " + std::to_string(nDim) + ")", AT_);       \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SIDE_ACCESSOR(side)                                                                               \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIR_ACCESSOR(orientation)                                                                         \
  do {                                                                                                                 \
  } while(false)
#endif

/// Namespace for auxiliary functions/classes
namespace maia {
namespace dg {
namespace collector {

/// Class that represents DG element collector.
template <MInt nDim, class SysEqn>
class SurfaceCollector : public maia::container::Container<SurfaceCollector<nDim, SysEqn>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<SurfaceCollector<nDim, SysEqn>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<SurfaceCollector<nDim, SysEqn>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;

 public:
  // Types
  template <class T>
  using Invalid = typename maia::dg::collector::Invalid<T>;

  // Constructor
  /// Default c'tor does nothing
  constexpr SurfaceCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessor
  MInt& orientation(const MInt srfcId);
  MInt orientation(const MInt srfcId) const;
  MInt& fineCellId(const MInt srfcId);
  MInt fineCellId(const MInt srfcId) const;
  MInt& internalSideId(const MInt srfcId);
  MInt internalSideId(const MInt srfcId) const;
  MInt& polyDeg(const MInt srfcId);
  MInt polyDeg(const MInt srfcId) const;
  MInt& noNodes1D(const MInt srfcId);
  MInt noNodes1D(const MInt srfcId) const;
  MInt noNodesXD(const MInt srfcId) const;
  MLong& globalId(const MInt srfcId);
  MLong globalId(const MInt srfcId) const;
  MFloat& coords(const MInt srfcId, const MInt dir);
  MFloat coords(const MInt srfcId, const MInt dir) const;
  MFloat& nodeCoords(const MInt srfcId);
  MFloat nodeCoords(const MInt srfcId) const;
  MFloat& flux(const MInt srfcId);
  MFloat flux(const MInt srfcId) const;
  MInt& nghbrElementIds(const MInt srfcId, const MInt side);
  MInt nghbrElementIds(const MInt srfcId, const MInt side) const;
  MFloat& variables(const MInt srfcId, const MInt side);
  MFloat& variables(const MInt srfcId, const MInt side) const;
  MFloat& nodeVars(const MInt srfcId, const MInt side);
  MFloat& nodeVars(const MInt srfcId, const MInt side) const;

  /// Return maximum polynomial degree
  MInt maxPolyDeg() const { return m_maxPolyDeg; }

  // Set maximum polynomial degree
  void maxPolyDeg(const MInt maxPolyDeg_) { m_maxPolyDeg = maxPolyDeg_; };

  /// Return maximum number of nodes 1D
  MInt maxNoNodes1D() const { return m_maxNoNodes1D; }

  // Set maximum number of nodes
  void maxNoNodes1D(const MInt maxNoNodes1D_) {
    m_maxNoNodes1D = maxNoNodes1D_;
    m_maxNoNodesXD = ipow(m_maxNoNodes1D, nDim - 1);
  };

  /// Return maximum number if nodes XD
  MInt maxNoNodesXD() const { return m_maxNoNodesXD; }

  // Set number of node variables
  void noNodeVars(const MInt noNodeVars_) { m_noNodeVars = noNodeVars_; };

  /// Return number of node variables
  constexpr MInt noNodeVars() const { return m_noNodeVars; }

  constexpr MInt noValuesVariables() const { return maxNoNodesXD() * SysEqn::noVars(); }

  constexpr MInt noValuesNodeVars() const { return maxNoNodesXD() * noNodeVars(); }

  constexpr MInt noValuesNodeCoordinates() const { return maxNoNodesXD() * nDim; }

  constexpr MInt noValuesFlux() const { return maxNoNodesXD() * SysEqn::noVars(); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  MInt m_maxPolyDeg = -1;
  MInt m_maxNoNodes1D = -1;
  MInt m_maxNoNodesXD = -1;

  // Number of node variables
  MInt m_noNodeVars = -1;

  // Data containers
  Storage<MInt> m_orientation{};     //!< Orientiation of face normal vector
  Storage<MInt> m_nghbrElementIds{}; //!< Nghbr. element ids (or -1 of no neighbor)
  Storage<MInt> m_fineCellId{};      //!< H-refinement fine element id
  Storage<MInt> m_internalSideId{};  //!< Side id of the internal element
  Storage<MInt> m_polyDeg{};         //!< Polynomial degrees of the neighbors
  Storage<MInt> m_noNodes1D{};       //!< Number of nodes 1D of the neighbors

  // Note: The following global id seemed to be a good idea when it was
  //       introduced in July 2013. There are certain advantages of having a
  //       globally unique surface id available permanently, especially for
  //       debugging purposes. The disadvantage, however, is of course thei
  //       increased storage requirements. Should this not be used more
  //       extensively in the future (or only during initialization (as it is
  //       now) where one could just recalculate the id), this should be removed
  //       again.
  // Note: The global surface id is calculated as the lower global cell id of
  //       both adjacent cells (if existing), which is multiplied by 2*nDim and
  //       then added to the direction id (0 = -x, 1 = +x, 2 = -y, etc.)
  //       relative to that cell.
  //       Example: global cell ids: 24, 13; nDim: 2; directionId: 2;
  //                global surface id: 54
  Storage<MLong> m_globalId{};

  Storage<MFloat> m_variables{};       //!< Surface variables
  Storage<MFloat> m_nodeVars{};        //!< Additional variables at each node
  Storage<MFloat> m_coordinates{};     //!< Coordinates of the surface center
  Storage<MFloat> m_nodeCoordinates{}; //!< Integration node coordinates
  Storage<MFloat> m_flux{};            //!< Store the numerical (Riemann) flux
};


/// Accessor for orientation.
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::orientation(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_orientation[srfcId];
}
/// Accessor for orientation (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::orientation(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_orientation[srfcId];
}


/// Accessor for fine cell id.
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::fineCellId(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_fineCellId[srfcId];
}
/// Accessor for fine cell id (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::fineCellId(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_fineCellId[srfcId];
}


/// Accessor for internal side id.
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::internalSideId(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_internalSideId[srfcId];
}
/// Accessor for internal side id (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::internalSideId(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_internalSideId[srfcId];
}


/// Accessor for polynomial degree.
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::polyDeg(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_polyDeg[srfcId];
}
/// Accessor for polynomial degree (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::polyDeg(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_polyDeg[srfcId];
}

/// Accessor for number of nodes 1D
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::noNodes1D(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_noNodes1D[srfcId];
}
/// Accessor for number of nodes 1D (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::noNodes1D(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_noNodes1D[srfcId];
}
/// Accessor for number of nodes XD (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::noNodesXD(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return ipow(m_noNodes1D[srfcId], nDim - 1);
}

/// Accessor for global id.
template <MInt nDim, class SysEqn>
MLong& SurfaceCollector<nDim, SysEqn>::globalId(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_globalId[srfcId];
}
/// Accessor for global id (const version).
template <MInt nDim, class SysEqn>
MLong SurfaceCollector<nDim, SysEqn>::globalId(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_globalId[srfcId];
}


/// Accessor for coordinates.
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::coords(const MInt srfcId, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_DIR_ACCESSOR(dir);
  return m_coordinates[srfcId * nDim + dir];
}
/// Accessor for coordinates (const version).
template <MInt nDim, class SysEqn>
MFloat SurfaceCollector<nDim, SysEqn>::coords(const MInt srfcId, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_DIR_ACCESSOR(dir);
  return m_coordinates[srfcId * nDim + dir];
}


/// Accessor for node coordinates.
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::nodeCoords(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_nodeCoordinates[srfcId * noValuesNodeCoordinates()];
}
/// Accessor for node coordinates (const version).
template <MInt nDim, class SysEqn>
MFloat SurfaceCollector<nDim, SysEqn>::nodeCoords(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_nodeCoordinates[srfcId * noValuesNodeCoordinates()];
}


/// Accessor for flux.
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::flux(const MInt srfcId) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_flux[srfcId * noValuesFlux()];
}
/// Accessor for flux (const version).
template <MInt nDim, class SysEqn>
MFloat SurfaceCollector<nDim, SysEqn>::flux(const MInt srfcId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  return m_flux[srfcId * noValuesFlux()];
}


/// Accessor for neighbor element ids.
template <MInt nDim, class SysEqn>
MInt& SurfaceCollector<nDim, SysEqn>::nghbrElementIds(const MInt srfcId, const MInt side) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_nghbrElementIds[2 * srfcId + side];
}
/// Accessor for neighbor element ids (const version).
template <MInt nDim, class SysEqn>
MInt SurfaceCollector<nDim, SysEqn>::nghbrElementIds(const MInt srfcId, const MInt side) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_nghbrElementIds[2 * srfcId + side];
}


/// Accessor for variables.
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::variables(const MInt srfcId, const MInt side) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_variables[2 * noValuesVariables() * srfcId + side * noValuesVariables()];
}
/// Accessor for variables (const version).
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::variables(const MInt srfcId, const MInt side) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_variables[2 * noValuesVariables() * srfcId + side * noValuesVariables()];
}


/// Accessor for node variables.
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::nodeVars(const MInt srfcId, const MInt side) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_nodeVars[2 * noValuesNodeVars() * srfcId + side * noValuesNodeVars()];
}
/// Accessor for node variables (const version).
template <MInt nDim, class SysEqn>
MFloat& SurfaceCollector<nDim, SysEqn>::nodeVars(const MInt srfcId, const MInt side) const {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(srfcId);
  ENSURE_VALID_SIDE_ACCESSOR(side);
  return m_nodeVars[2 * noValuesNodeVars() * srfcId + side * noValuesNodeVars()];
}


/// Reset SurfaceCollector, re-create data structures.
template <MInt nDim, class SysEqn>
void SurfaceCollector<nDim, SysEqn>::reset() {
  resetStorage(1, m_orientation);
  resetStorage(2, m_nghbrElementIds);
  resetStorage(1, m_fineCellId);
  resetStorage(1, m_internalSideId);
  resetStorage(1, m_polyDeg);
  resetStorage(1, m_noNodes1D);
  resetStorage(1, m_globalId);
  resetStorage(2 * noValuesVariables(), m_variables);
  resetStorage(2 * noValuesNodeVars(), m_nodeVars);
  resetStorage(nDim, m_coordinates);
  resetStorage(noValuesNodeCoordinates(), m_nodeCoordinates);
  resetStorage(noValuesFlux(), m_flux);
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim, class SysEqn>
void SurfaceCollector<nDim, SysEqn>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  fill_invalid(m_orientation, begin, end);

  fill_invalid(m_nghbrElementIds, begin, end, 2);

  fill_invalid(m_fineCellId, begin, end);

  fill_invalid(m_internalSideId, begin, end);

  fill_invalid(m_polyDeg, begin, end);

  fill_invalid(m_noNodes1D, begin, end);

  fill_invalid(m_globalId, begin, end);

  fill_invalid(m_variables, begin, end, 2 * noValuesVariables());

  fill_invalid(m_nodeVars, begin, end, 2 * noValuesNodeVars());

  fill_invalid(m_coordinates, begin, end, nDim);

  fill_invalid(m_nodeCoordinates, begin, end, noValuesNodeCoordinates());

  fill_invalid(m_flux, begin, end, noValuesFlux());
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim, class SysEqn>
template <class Functor, class T>
void SurfaceCollector<nDim, SysEqn>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                                    const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef DGCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  copyData(source.m_orientation, m_orientation, c, begin, end, destination);

  copyData(source.m_nghbrElementIds, m_nghbrElementIds, c, begin, end, destination, 2);

  copyData(source.m_fineCellId, m_fineCellId, c, begin, end, destination);

  copyData(source.m_internalSideId, m_internalSideId, c, begin, end, destination);

  copyData(source.m_polyDeg, m_polyDeg, c, begin, end, destination);

  copyData(source.m_noNodes1D, m_noNodes1D, c, begin, end, destination);

  copyData(source.m_globalId, m_globalId, c, begin, end, destination);

  copyData(source.m_variables, m_variables, c, begin, end, destination, 2 * noValuesVariables());

  copyData(source.m_nodeVars, m_nodeVars, c, begin, end, destination, 2 * noValuesNodeVars());

  copyData(source.m_coordinates, m_coordinates, c, begin, end, destination, nDim);

  copyData(source.m_nodeCoordinates, m_nodeCoordinates, c, begin, end, destination, noValuesNodeCoordinates());

  copyData(source.m_flux, m_flux, c, begin, end, destination, noValuesFlux());
}

} // namespace collector
} // namespace dg
} // namespace maia


// Undefine macros that should not be used outside this file
#undef DGCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_SIDE_ACCESSOR
#undef ENSURE_VALID_DIR_ACCESSOR

#endif // MAIA_DGSURFACECOLLECTOR_H
