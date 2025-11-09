// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDTREE_H_
#define GRIDTREE_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/parallelio.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "cartesiangridcellproperties.h"
#include "compiler_config.h"
#include "config.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define GRIDTREE_SOA_MEMORY_LAYOUT

// The macro 'GRIDTREE_SANITY_CHECKS' enables (potentially expensive) sanity checks for many grid
// operations (but not the accessors). It is enabled for build type "debug"

// The macro 'GRIDTREE_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity checks
// for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define GRIDTREE_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#ifdef GRIDTREE_SANITY_CHECKS_ACCESSORS
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_CHILD_POSITION_ACCESSOR(pos)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(pos >= 0 && pos < noChildrenPerNode(),                                                       \
                          "child position = " + std::to_string(pos) + " out-of-bounds [0, "                            \
                              + std::to_string(noChildrenPerNode()) + ")",                                             \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir)                                                                        \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(dir >= 0 && dir < noNeighborsPerNode(),                                                      \
                          "neighbor direction = " + std::to_string(dir) + " out-of-bounds [0, "                        \
                              + std::to_string(noNeighborsPerNode()) + ")",                                            \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        dir >= 0 && dir < nDim,                                                                                        \
        "coordinate direction " + std::to_string(dir) + " out-of-bounds [0, " + std::to_string(nDim) + ")", AT_);      \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != Cell::NumProperties, "Invalid property", AT_);                                          \
  } while(false)
#define ENSURE_VALID_SOLVER_ID_ACCESSOR(id)                                                                            \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noSolvers(), "Invalid solver id", AT_);                                      \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_CHILD_POSITION_ACCESSOR(pos)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir)                                                                        \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_SOLVER_ID_ACCESSOR(id)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace grid {
namespace tree {

/// Underlying enum type for property access
using Cell = GridCell;

/// Underlying bitset type for property storage
using PropertyBitsetType = maia::grid::cell::BitsetType;

/// Maximum number of supported solvers
/// Note: If you increase this value, you also need to adapt the data type in the grid file
static constexpr MInt MAIA_MULTISOLVER_MAX_NO_SOLVERS = 8;

/// Underlying bitset type for solver use storage (Note: If there are more solvers, change size here)
using SolverBitsetType = std::bitset<MAIA_MULTISOLVER_MAX_NO_SOLVERS>;

// Type traits for invalid values. These values are used to initialize/erase nodes
template <class T>
struct Invalid {};

// Invalid value for ids is 'INT_MIN'
template <>
struct Invalid<MInt> {
  static constexpr MInt value() { return std::numeric_limits<MInt>::min(); }
};

template <>
struct Invalid<MLong> {
  static constexpr MLong value() { return std::numeric_limits<MLong>::min(); }
};

// Invalid value for chars is 'INT_MIN'
template <>
struct Invalid<MChar> {
  static constexpr MInt value() { return std::numeric_limits<MChar>::min(); }
};

// Invalid value for chars is 'INT_MIN'
template <>
struct Invalid<MUchar> {
  static constexpr MInt value() { return std::numeric_limits<MUchar>::min(); }
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

// Invalid value for bitsets is '0'
template <>
struct Invalid<PropertyBitsetType> {
  static constexpr PropertyBitsetType value() { return {0}; }
};
template <>
struct Invalid<SolverBitsetType> {
  static constexpr SolverBitsetType value() { return {0}; }
};


////////////////////////////////////////////////////////////////////////////////////////////////////
// Concept
////////////////////////////////////////////////////////////////////////////////////////////////////
//
// - Only high-level API is exposed to users (i.e., is public)
// - High-level API calls preserve tree consistency
// - Bounds and sanity checking in each high-level API call
// - Low-level API performs only memory operations (i.e., is auxiliary to the high-level calls)
// - Low-level API and does usually not perform bounds checking
//
////////////////////////////////////////////////////////////////////////////////////////////////////


/// Class that represents grid tree and contains all relevant per-node data.
template <MInt nDim>
class Tree : public maia::container::Container<Tree<nDim>, maia::grid::tree::Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<Tree<nDim>, maia::grid::tree::Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<Tree<nDim>, maia::grid::tree::Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;


 public:
  // Types
  template <class T>
  using Invalid = typename maia::grid::tree::Invalid<T>;
  using Cell = maia::grid::tree::Cell;
  using PropertyBitsetType = maia::grid::tree::PropertyBitsetType;
  using SolverBitsetType = maia::grid::tree::SolverBitsetType;

  // Constructors
  /// Default c'tor does nothing
  constexpr Tree() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Parent-child relationship
  MLong& parent(const MInt id);
  MLong parent(const MInt id) const;
  MBool hasParent(const MInt id) const;
  MLong& child(const MInt id, const MInt pos);
  MLong child(const MInt id, const MInt pos) const;
  MBool hasChild(const MInt id, const MInt pos) const;
  MBool hasChildren(const MInt id) const;
  MBool hasChildren(const MInt id, const MInt solverId) const;
  MInt noChildren(const MInt id) const;

  // Neighbors
  MLong& neighbor(const MInt id, const MInt dir);
  MLong neighbor(const MInt id, const MInt dir) const;
  MBool hasNeighbor(const MInt id, const MInt dir) const;
  MBool hasAnyNeighbor(const MInt id, const MInt dir) const;

  // Other data fields
  MLong& globalId(const MInt id);
  MLong globalId(const MInt id) const;
  MInt& level(const MInt id);
  MInt level(const MInt id) const;
  MFloat& coordinate(const MInt id, const MInt dim);
  MFloat coordinate(const MInt id, const MInt dim) const;
  const MFloat* coordinate(const MInt id) const;
  MFloat& weight(const MInt id);
  MFloat weight(const MInt id) const;

  SolverBitsetType::reference solver(const MInt id, const MInt solverId);
  MBool solver(const MInt id, const MInt solverId) const;
  void resetSolver(const MInt id);
  MUlong solverToBits(const MInt id) const;
  void solverFromBits(const MInt id, const MUlong bits);
  MBool cellHasSolver(const MInt cellId, const MInt solverId);
  SolverBitsetType& solverBits(const MInt id);

  MBool isLeafCell(const MInt id) const;
  SolverBitsetType::reference isLeafCell(const MInt id, const MInt solverId);
  MBool isLeafCell(const MInt id, const MInt solverId) const;
  void resetIsLeafCell(const MInt id);
  SolverBitsetType& leafCellBits(const MInt id);

  PropertyBitsetType::reference hasProperty(const MInt id, const Cell p);
  MBool hasProperty(const MInt id, const Cell p) const;
  void resetProperties(const MInt id);
  MUlong propertiesToBits(const MInt id) const;
  void propertiesFromBits(const MInt id, const MUlong bits);
  MString propertiesToString(const MInt id) const;
  PropertyBitsetType& properties(const MInt id);

  // Other data fields (subject to change)
  MInt& noOffsprings(const MInt id);
  MInt noOffsprings(const MInt id) const;
  MFloat& workload(const MInt id);
  MFloat workload(const MInt id) const;

  // Entries per tree node
  /// Return maximum number of children per node
  static constexpr MInt noChildrenPerNode() { return 4 * nDim - 4; }
  /// Return maximum number of same-level neighbors per node
  static constexpr MInt noNeighborsPerNode() { return 2 * nDim; }

  /// Return opposite direction for given neighbor direction
  static constexpr MInt oppositeNeighborDir(const MInt dir) { return dir + 1 - 2 * (dir % 2); }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::grid::cell::p(Cell::NumProperties); }

  /// Return maximum number of supported solvers
  static constexpr MInt maxNoSolvers() { return SolverBitsetType().size(); }

  // Methods related to multi-solver features
  /// Return currently set number of solvers
  constexpr MInt noSolvers() const { return m_noSolvers; }

  // Set number of solvers
  void setNoSolvers(const MInt count);

  MInt noNodesBySolver(const MInt solverId) const;
  MInt nodesBySolver(const MInt solverId, MInt* const ids) const;

  MInt capacity() { return m_parentIds.capacity(); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);
  void deleteConnectivity(const MInt begin, const MInt end);
  void moveConnectivity(const MInt begin, const MInt end, const MInt to);

  // Data containers
  Storage<MLong> m_parentIds{};
  Storage<MLong> m_childIds{};
  Storage<MLong> m_neighborIds{};
  Storage<MLong> m_globalIds{};
  Storage<MInt> m_levels{};
  Storage<MFloat> m_coordinates{};
  Storage<MFloat> m_weight{};
  Storage<SolverBitsetType> m_solver{};
  Storage<SolverBitsetType> m_isLeafCell{};
  Storage<PropertyBitsetType> m_properties{};
  Storage<MInt> m_noOffsprings{};
  Storage<MFloat> m_workload{};

  /// Number of solvers that are actively using this tree
  MInt m_noSolvers = -1;
};


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void Tree<nDim>::reset() {
  resetStorage(1, m_parentIds);
  resetStorage(noChildrenPerNode(), m_childIds);
  resetStorage(noNeighborsPerNode(), m_neighborIds);
  resetStorage(1, m_levels);
  resetStorage(1, m_globalIds);
  resetStorage(nDim, m_coordinates);
  resetStorage(1, m_weight);
  resetStorage(1, m_solver);
  resetStorage(1, m_isLeafCell);
  resetStorage(1, m_properties);
  resetStorage(1, m_noOffsprings);
  resetStorage(1, m_workload);
}


/// Accessor for parent node.
template <MInt nDim>
MLong& Tree<nDim>::parent(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_parentIds[id];
}
/// Accessor for parent node (const version).
template <MInt nDim>
MLong Tree<nDim>::parent(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_parentIds[id];
}


/// Return whether node has parent.
template <MInt nDim>
MBool Tree<nDim>::hasParent(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_parentIds[id] > -1;
}


/// Accessor for child node.
template <MInt nDim>
MLong& Tree<nDim>::child(const MInt id, const MInt pos) {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_CHILD_POSITION_ACCESSOR(pos);
  return m_childIds[id * noChildrenPerNode() + pos];
}
/// Accessor for child node (const version).
template <MInt nDim>
MLong Tree<nDim>::child(const MInt id, const MInt pos) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_CHILD_POSITION_ACCESSOR(pos);
  return m_childIds[id * noChildrenPerNode() + pos];
}


/// Return whether node has child at given position.
template <MInt nDim>
MBool Tree<nDim>::hasChild(const MInt id, const MInt pos) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_CHILD_POSITION_ACCESSOR(pos);
  return m_childIds[id * noChildrenPerNode() + pos] > -1;
}


/// Return whether node has any children.
template <MInt nDim>
MBool Tree<nDim>::hasChildren(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return noChildren(id) > 0;
}

/// Return whether node has any children living on the specified solver (Lennart).
template <MInt nDim>
MBool Tree<nDim>::hasChildren(const MInt id, const MInt solverId) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return std::any_of(&m_childIds[id * noChildrenPerNode()], &m_childIds[id * noChildrenPerNode() + noChildrenPerNode()],
                     [&](MLong c) { return c > -1 && solver(c, solverId); });
}

/// Return number of children of given node.
template <MInt nDim>
MInt Tree<nDim>::noChildren(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return std::count_if(&m_childIds[id * noChildrenPerNode() + 0],
                       &m_childIds[id * noChildrenPerNode() + noChildrenPerNode()],
                       [](const ParallelIo::size_type childId) { return childId > -1; });
}


/// Accessor for neighbor node.
template <MInt nDim>
MLong& Tree<nDim>::neighbor(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir);
  return m_neighborIds[id * noNeighborsPerNode() + dir];
}
/// Accessor for neighbor node (const version).
template <MInt nDim>
MLong Tree<nDim>::neighbor(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir);
  return m_neighborIds[id * noNeighborsPerNode() + dir];
}


/// Return whether node has same-level neighbor in given direction.
template <MInt nDim>
MBool Tree<nDim>::hasNeighbor(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir);
  return m_neighborIds[id * noNeighborsPerNode() + dir] > -1;
}
/// Return whether node or its parent has neighbor in given direction.
template <MInt nDim>
MBool Tree<nDim>::hasAnyNeighbor(const MInt id, const MInt dir) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR(dir);
  return hasNeighbor(id, dir) || (hasParent(id) && hasNeighbor(parent(id), dir));
}


/// Accessor for global id.
template <MInt nDim>
MLong& Tree<nDim>::globalId(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_globalIds[id];
}
/// Accessor for global id (const version).
template <MInt nDim>
MLong Tree<nDim>::globalId(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_globalIds[id];
}


/// Accessor for level.
template <MInt nDim>
MInt& Tree<nDim>::level(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_levels[id];
}
/// Accessor for level (const version).
template <MInt nDim>
MInt Tree<nDim>::level(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_levels[id];
}


/// Accessor for coordinates.
template <MInt nDim>
MFloat& Tree<nDim>::coordinate(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}
/// Accessor for coordinates (const version).
template <MInt nDim>
MFloat Tree<nDim>::coordinate(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}

template <MInt nDim>
const MFloat* Tree<nDim>::coordinate(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return &m_coordinates[id * nDim];
}


/// Accessor for weight.
template <MInt nDim>
MFloat& Tree<nDim>::weight(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_weight[id];
}
/// Accessor for weight (const version).
template <MInt nDim>
MFloat Tree<nDim>::weight(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_weight[id];
}


/// Accessor for noOffsprings.
template <MInt nDim>
MInt& Tree<nDim>::noOffsprings(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noOffsprings[id];
}
/// Accessor for noOffsprings (const version).
template <MInt nDim>
MInt Tree<nDim>::noOffsprings(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noOffsprings[id];
}


/// Accessor for workload.
template <MInt nDim>
MFloat& Tree<nDim>::workload(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_workload[id];
}
/// Accessor for workload (const version).
template <MInt nDim>
MFloat Tree<nDim>::workload(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_workload[id];
}


/// Accessor for solver usage.
template <MInt nDim>
Tree<nDim>::SolverBitsetType::reference Tree<nDim>::solver(const MInt id, const MInt solverId) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SOLVER_ID_ACCESSOR(solverId);
  return m_solver[id][solverId];
}
/// Accessor for solver usage (const version).
template <MInt nDim>
MBool Tree<nDim>::solver(const MInt id, const MInt solverId) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SOLVER_ID_ACCESSOR(solverId);
  return m_solver[id][solverId];
}
/// Reset all solver use.
template <MInt nDim>
void Tree<nDim>::resetSolver(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_solver[id].reset();
}
/// Convert solver usage to bits.
template <MInt nDim>
MUlong Tree<nDim>::solverToBits(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_solver[id].to_ulong();
}
/// Convert solver usage from bits.
template <MInt nDim>
void Tree<nDim>::solverFromBits(const MInt id, const MUlong bits) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_solver[id] = SolverBitsetType(bits);
}
/// Check if solver is contained for cell
template <MInt nDim>
MBool Tree<nDim>::cellHasSolver(const MInt cellId, const MInt solverId) {
  return m_solver[cellId][solverId];
}
/// Accessor for properties.
template <MInt nDim>
SolverBitsetType& Tree<nDim>::solverBits(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_solver[id];
}


/// Accessor for isLeafCell usage (const version).
template <MInt nDim>
MBool Tree<nDim>::isLeafCell(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return (m_isLeafCell[id].any() && m_isLeafCell[id] == m_solver[id]);
}
/// Accessor for isLeafCell usage.
template <MInt nDim>
Tree<nDim>::SolverBitsetType::reference Tree<nDim>::isLeafCell(const MInt id, const MInt solverId) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SOLVER_ID_ACCESSOR(solverId);
  // ASSERT( m_solver[id][solverId], "Invalid cell accessed, not belonging to solver!" );
  return m_isLeafCell[id][solverId];
}
/// Accessor for isLeafCell usage (const version).
template <MInt nDim>
MBool Tree<nDim>::isLeafCell(const MInt id, const MInt solverId) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_SOLVER_ID_ACCESSOR(solverId);
  // ASSERT( m_solver[id][solverId], "Invalid cell accessed, not belonging to solver!" );
  return m_isLeafCell[id][solverId];
}
/// Reset all isLeafCell.
template <MInt nDim>
void Tree<nDim>::resetIsLeafCell(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_isLeafCell[id].reset();
}
/// Accessor for properties.
template <MInt nDim>
SolverBitsetType& Tree<nDim>::leafCellBits(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_isLeafCell[id];
}

/// Accessor for properties.
template <MInt nDim>
Tree<nDim>::PropertyBitsetType::reference Tree<nDim>::hasProperty(const MInt id, const Cell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::grid::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool Tree<nDim>::hasProperty(const MInt id, const Cell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::grid::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void Tree<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties[id].reset();
}
/// Convert properties to bits.
template <MInt nDim>
MUlong Tree<nDim>::propertiesToBits(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id].to_ulong();
}
/// Convert properties from bits.
template <MInt nDim>
void Tree<nDim>::propertiesFromBits(const MInt id, const MUlong bits) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties[id] = PropertyBitsetType(bits);
}
/// Convert properties to string.
template <MInt nDim>
MString Tree<nDim>::propertiesToString(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id].to_string();
}
/// Accessor for properties.
template <MInt nDim>
PropertyBitsetType& Tree<nDim>::properties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}


/// Set number of solvers
template <MInt nDim>
void Tree<nDim>::setNoSolvers(const MInt count) {
  MAIA_CONTAINER_ENSURE(count > 0, "Number of solvers must be greater than zero.", AT_);
  MAIA_CONTAINER_ENSURE(count <= maxNoSolvers(),
                        "Number of solvers out-of-bounds [0, " + std::to_string(maxNoSolvers()) + "].", AT_);
  m_noSolvers = count;
}


/// Return number of nodes for a given solver
template <MInt nDim>
MInt Tree<nDim>::noNodesBySolver(const MInt solverId) const {
  MAIA_CONTAINER_ENSURE(solverId >= 0 && solverId < noSolvers(),
                        "solver id = " + std::to_string(solverId) + " out-of-bounds [0, " + std::to_string(noSolvers())
                            + ")",
                        AT_);

  // Count nodes with solver bit set
  MInt noNodes = 0;
  for(MInt id = 0; id < this->size(); id++) {
    noNodes += solver(id, solverId);
  }

  return noNodes;
}


/// Generate list of node ids that are used by a given solver and return number of used nodes.
template <MInt nDim>
MInt Tree<nDim>::nodesBySolver(const MInt solverId, MInt* const ids) const {
  MAIA_CONTAINER_ENSURE(solverId >= 0 && solverId < noSolvers(),
                        "solver id = " + std::to_string(solverId) + " out-of-bounds [0, " + std::to_string(noSolvers())
                            + ")",
                        AT_);
  MAIA_CONTAINER_ENSURE(ids != nullptr, "Data pointer must not be null", AT_);

  // Loop over all nodes and store ids with solver bit set
  MInt count = 0;
  for(MInt id = 0; id < this->size(); id++) {
    // Skip node if solver bit not set
    if(!solver(id, solverId)) {
      continue;
    }

    // Add id to list and increase count
    ids[count] = id;
    count++;
  }

  return count;
}


/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void Tree<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // Parent
  fill_invalid(m_parentIds, begin, end);

  // Children
  fill_invalid(m_childIds, begin, end, noChildrenPerNode());

  // Neighbors
  fill_invalid(m_neighborIds, begin, end, noNeighborsPerNode());

  // Global id
  fill_invalid(m_globalIds, begin, end);

  // Level
  fill_invalid(m_levels, begin, end);

  // Coordinates
  fill_invalid(m_coordinates, begin, end, nDim);

  // Solver usage
  fill_invalid(m_solver, begin, end);

  // Solver usage
  fill_invalid(m_isLeafCell, begin, end);

  // Weight
  fill_invalid(m_weight, begin, end);

  // Properties
  fill_invalid(m_properties, begin, end);

  // No offsprings
  fill_invalid(m_noOffsprings, begin, end);

  // Workload
  fill_invalid(m_workload, begin, end);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void Tree<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef GRIDTREE_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  // Parent
  copyData(source.m_parentIds, m_parentIds, c, begin, end, destination);

  // Children
  copyData(source.m_childIds, m_childIds, c, begin, end, destination, noChildrenPerNode());

  // Neighbors
  copyData(source.m_neighborIds, m_neighborIds, c, begin, end, destination, noNeighborsPerNode());

  // Global id
  copyData(source.m_globalIds, m_globalIds, c, begin, end, destination);

  // Level
  copyData(source.m_levels, m_levels, c, begin, end, destination);

  // Coordinates
  copyData(source.m_coordinates, m_coordinates, c, begin, end, destination, nDim);

  // Solver usage
  copyData(source.m_solver, m_solver, c, begin, end, destination);

  // Leaf cell flags
  copyData(source.m_isLeafCell, m_isLeafCell, c, begin, end, destination);

  // Weight
  copyData(source.m_weight, m_weight, c, begin, end, destination);

  // Properties
  copyData(source.m_properties, m_properties, c, begin, end, destination);

  // No offsprings
  copyData(source.m_noOffsprings, m_noOffsprings, c, begin, end, destination);

  // Workload
  copyData(source.m_workload, m_workload, c, begin, end, destination);
}


/// Update parent/children/neighbors *before* nodes are erased.
template <MInt nDim>
void Tree<nDim>::deleteConnectivity(const MInt begin, const MInt end) {
  for(MInt i = begin; i < end; i++) {
    // Parent
    if(hasParent(i)) {
      const MInt p = parent(i);
      for(MInt j = 0; j < noChildrenPerNode(); j++) {
        if(child(p, j) == i) {
          child(p, j) = -1;
        }
      }
    }

    // Children
    for(MInt j = 0; j < noChildrenPerNode(); j++) {
      if(hasChild(i, j)) {
        parent(child(i, j)) = -1;
      }
    }

    // Neighbors
    for(MInt j = 0; j < noNeighborsPerNode(); j++) {
      if(hasNeighbor(i, j)) {
        neighbor(neighbor(i, j), oppositeNeighborDir(j)) = -1;
      }
    }
  }
}


/// Update parent/children/neighbors after nodes have moved.
template <MInt nDim>
void Tree<nDim>::moveConnectivity(const MInt begin, const MInt end, const MInt to) {
  // Auxiliary method for checking if a given id is within the original range that was moved
  auto inMovedRange = [begin, end](const MInt id) { return (id >= begin && id < end); };

  // General strategy:
  // 1) Loop over moved nodes and check all tree connections (parents/children/neighbors)
  // 2) If a given connection is to a node that was moved: apply offset to current node
  // 3) If a given connection is to a node that was not moved: change connectivity in other node
  for(MInt from = begin; from < end; from++) {
    const MInt distance = to - begin;
    const MInt destination = from + distance;

    // Parent
    if(hasParent(destination)) {
      const MInt p = parent(destination);
      if(inMovedRange(p)) {
        parent(destination) += distance;
      } else {
        for(MInt j = 0; j < noChildrenPerNode(); j++) {
          if(child(p, j) == from) {
            child(p, j) = destination;
          }
        }
      }
    }

    // Children
    for(MInt j = 0; j < noChildrenPerNode(); j++) {
      if(hasChild(destination, j)) {
        const MInt c = child(destination, j);
        if(inMovedRange(c)) {
          child(destination, j) += distance;
        } else {
          parent(c) = destination;
        }
      }
    }

    // Neighbors
    for(MInt j = 0; j < noNeighborsPerNode(); j++) {
      if(hasNeighbor(destination, j)) {
        const MInt n = neighbor(destination, j);
        if(inMovedRange(n)) {
          neighbor(destination, j) += distance;
        } else {
          neighbor(n, oppositeNeighborDir(j)) = destination;
        }
      }
    }
  }
}

} // namespace tree
} // namespace grid
} // namespace maia


// Undefine macros that should not be used outside this file
#undef GRIDTREE_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_CHILD_POSITION_ACCESSOR
#undef ENSURE_VALID_NEIGHBOR_DIR_ACCESSOR
#undef ENSURE_VALID_COORDINATE_DIR_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef GRIDTREE_H_
