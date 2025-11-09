// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDPROXY_H_
#define GRIDPROXY_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "cartesiangrid.h"
#include "cartesiangridtree.h"
#include "typetraits.h"

#include "GEOM/geometry.h"

namespace maia {
namespace grid {

// Forward declare tree proxy
namespace tree {
template <MInt nDim>
class TreeProxy;
}

/// Class that acts like a fully functional grid for a single solver
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-16
template <MInt nDim>
class Proxy {
 public:
  // Types
  using Grid = CartesianGrid<nDim>;
  using Tree = maia::grid::tree::Tree<nDim>;
  using Cell = maia::grid::tree::Cell;
  using TreeProxy = maia::grid::tree::TreeProxy<nDim>;
  using Geom = Geometry<nDim>;

  // Methods
  // Constructor
  Proxy(const MInt solverId, Grid& grid_, Geometry<nDim>& geometry_);
  Proxy(const MInt solverId, Grid& grid_); // ParaView plugin

  // Destructor
  ~Proxy();

  // Add this accessor in case the full grid access is needed
  Grid& raw() { return m_grid; }
  const Grid& raw() const { return m_grid; }

  // Return proxy grid tree
  const TreeProxy& tree() const { return m_tree; }

  // Other methods
  void update();
  void updateOther();
  void updateLeafCellExchange();
  void initGridMap();
  void updateGridMap();
  void resizeGridMap(const MInt solverSize);
  void correctAzimuthalHaloCells();
  void findDirectNghbrs(const MInt, std::vector<MInt>&);
  void findNeighborHood(const MInt, const MInt, std::vector<MInt>&);
  void getAllLeafChilds(const MInt, std::vector<MInt>&);
  void updatePartitionCellOffsets();
  MInt domainContainingCell(MInt globalId) const;

 private:
  // Internal methods for updating the proxy class and its tree proxy
  void updateParallelizationInfo();
  void checkOffsetConsistency() const;
  void checkNeighborConsistency() const;
  void setupWindowHaloConnectivityOnLeafLvl(std::map<MInt, MInt>&);
  void checkWindowHaloConsistency() const;
  void checkNeighborConsistencyAzimuthal() const;
  void checkWindowHaloConsistencyAzimuthal() const;
  void descendStoreGlobalId(MInt cellId, MInt& localCnt);
  void updateCutOff();
  void updateTreeData();

  void getLocalSameLevelCellIds(std::vector<MInt>& levelCellId, const MInt level);

  void setSolverFlagsForAddedSolver();

 public:
  void updateGridInfo();
  void findEqualLevelNeighborsParDiagonal(MBool idsAreGlobal = true);
  void findCartesianNghbIds();
  MInt getAdjacentGridCells(
      MInt cellId, MInt noLayers, MIntScratchSpace& adjacentCells, MInt level, MInt diagonalNeighbors = 0);
  MInt getNghbrCells(MInt cellId, MInt level, MInt* nghbrs, MInt dir0, MInt dir1 = -1, MInt dir2 = -1);
  MInt findNeighborDomainId(const MLong globalId);
  MInt globalToLocalId(const MLong globalId) const { return tree().localId(globalId); }
  void smoothFilter(const MInt level, MFloat** value);

  MInt findContainingLeafCell(const MFloat* coord) {
    const MInt gridCellId = raw().findContainingLeafCell(coord, nullptr, solverId());
    if(gridCellId < 0) return gridCellId;
    ASSERT(solverFlag(gridCellId, solverId()), "");
    return tree().grid2solver(gridCellId);
  }
  MInt findContainingLeafCell(const MFloat* coord, const MInt startId, const MBool allowNonLeafHalo = false) {
    const MInt startGridId = tree().solver2grid(startId);
    const MInt gridCellId = raw().findContainingLeafCell(coord, startGridId, nullptr, solverId(), allowNonLeafHalo);
    if(gridCellId < 0) return gridCellId;
    ASSERT(solverFlag(gridCellId, solverId()), "");
    return tree().grid2solver(gridCellId);
  }

  /// Return whether the solver is active on the current domain
  MBool isActive() const { return m_isActive; }

  /// Return whether this solver is inactive on one of the global domains
  MBool hasInactiveRanks() const { return m_hasInactiveRanks; }

  /// Return solver id
  MInt solverId() const { return m_solverId; }

  MBool wasAdapted() const { return raw().wasAdapted(); }
  MBool wasBalanced() const { return raw().wasBalanced(); }

  // General grid & parallelization information
  MInt noCells() const { return tree().size(); }
  MInt noInternalCells() const { return m_noInternalCells; }

  MInt noNeighborDomains() const { return m_neighborDomains.size(); }
  MInt neighborDomain(const MInt id) const { return m_neighborDomains[id]; }
  const std::vector<MInt>& neighborDomains() const { return m_neighborDomains; }

  MInt noWindowCells(const MInt domainId) const { return m_windowCells[domainId].size(); }
  const MInt& windowCell(const MInt domainId, const MInt id) const { return m_windowCells[domainId][id]; }
  const std::vector<std::vector<MInt>>& windowCells() const { return m_windowCells; };

  MInt noHaloCells(const MInt domainId) const { return m_haloCells[domainId].size(); }
  const MInt& haloCell(const MInt domainId, const MInt id) const { return m_haloCells[domainId][id]; }
  const std::vector<std::vector<MInt>>& haloCells() const { return m_haloCells; };

  MInt noLeafSendNeighborDomains() const { return m_leafSendNeighborDomains.size(); }
  MInt leafSendNeighborDomain(const MInt id) const { return m_leafSendNeighborDomains[id]; }

  MInt noLeafRecvNeighborDomains() const { return m_leafRecvNeighborDomains.size(); }
  MInt leafRecvNeighborDomain(const MInt id) const { return m_leafRecvNeighborDomains[id]; }

  MInt noLeafWindowCells(const MInt domainId) const { return m_leafWindowCells[domainId].size(); }
  const MInt& leafWindowCell(const MInt domainId, const MInt id) const { return m_leafWindowCells[domainId][id]; }
  const std::vector<std::vector<MInt>>& leafWindowCells() const { return m_leafWindowCells; };

  MInt noLeafHaloCells(const MInt domainId) const { return m_leafHaloCells[domainId].size(); }
  const MInt& leafHaloCell(const MInt domainId, const MInt id) const { return m_leafHaloCells[domainId][id]; }
  const std::vector<std::vector<MInt>>& leafHaloCells() const { return m_leafHaloCells; };

  MInt leafRecSize() const { return m_leafRecvSize; }
  MInt leafSendSize() const { return m_leafSendSize; }

  // Azimuthal periodicity
  MInt noAzimuthalNeighborDomains() const { return m_azimuthalNeighborDomains.size(); }
  MInt azimuthalNeighborDomain(const MInt azimuthalId) const { return m_azimuthalNeighborDomains[azimuthalId]; }
  const std::vector<MInt>& azimuthalNeighborDomains() const { return m_azimuthalNeighborDomains; }

  // Window cells
  MInt noAzimuthalWindowCells(const MInt azimuthalDomainId) const {
    return m_azimuthalWindowCells[azimuthalDomainId].size();
  }
  const MInt& azimuthalWindowCell(const MInt azimuthalDomainId, const MInt id) const {
    return m_azimuthalWindowCells[azimuthalDomainId][id];
  }
  const std::vector<std::vector<MInt>>& azimuthalWindowCells() const { return m_azimuthalWindowCells; };

  // Halo cells
  MInt noAzimuthalHaloCells(const MInt azimuthalDomainId) const {
    return m_azimuthalHaloCells[azimuthalDomainId].size();
  }
  const MInt& azimuthalHaloCell(const MInt azimuthalDomainId, const MInt id) const {
    return m_azimuthalHaloCells[azimuthalDomainId][id];
  }
  const std::vector<std::vector<MInt>>& azimuthalHaloCells() const { return m_azimuthalHaloCells; };

  MInt noAzimuthalUnmappedHaloCells() const { return m_azimuthalUnmappedHaloCells.size(); }
  const MInt& azimuthalUnmappedHaloCell(const MInt id) const { return m_azimuthalUnmappedHaloCells[id]; }
  const MInt& azimuthalUnmappedHaloDomain(const MInt id) const { return m_azimuthalUnmappedHaloDomains[id]; }

  MInt noDomains() const { return m_noDomains; }
  MInt domainId() const { return m_domainId; }
  MPI_Comm mpiComm() const { return m_mpiComm; }

  const MLong& domainOffset(const MInt id) const { return m_domainOffsets[id]; }
  MInt maxLevel() const { return m_maxLevel; }

  // Methods that are just pass-through to the original grid
  MFloat reductionFactor() const { return raw().reductionFactor(); }
  MInt maxRefinementLevel() const { return m_maxRefinementLevel; }
  MInt maxUniformRefinementLevel() const { return raw().maxUniformRefinementLevel(); }
  MFloat lengthLevel0() const { return raw().lengthLevel0(); }
  MInt minLevel() const { return raw().minLevel(); }
  MInt newMinLevel() const { return raw().m_newMinLevel; }
  MInt centerOfGravity(MInt dir) const { return raw().centerOfGravity(dir); }
  MString gridInputFileName() const { return raw().gridInputFileName(); }
  MFloat cellLengthAtLevel(const MInt level) const { return raw().cellLengthAtLevel(level); }
  MFloat cellLengthAtCell(const MInt cellId) const { return raw().cellLengthAtCell(m_tree.solver2grid(cellId)); }
  MFloat halfCellLength(const MInt cellId) const { return raw().halfCellLength(m_tree.solver2grid(cellId)); }
  MFloat cellVolumeAtLevel(const MInt level) const { return raw().cellVolumeAtLevel(level); }
  constexpr MBool hasCutOff() const { return raw().hasCutOff(); }

  static constexpr MInt m_maxNoChilds = IPOW2(nDim);
  MLong bitOffset() const { return raw().bitOffset(); }
  MBool allowInterfaceRefinement() const { return raw().allowInterfaceRefinement(); }
  MLong noCellsGlobal() const { return raw().noCellsGlobal(); }
  MInt minCell(const MInt id) const {
    if(raw().minCell(id) > -1 && solverFlag(raw().minCell(id), solverId())) {
      return m_tree.grid2solver(raw().minCell(id));
    }
    return -1;
  }
  // NOTE: The number of min-cells is the same in the solver and the grid.
  //      But the ids might have changed during adaptation/balance!
  MInt noMinCells() const { return raw().noMinCells(); }
  MBool isPeriodic(const MInt cellId) {
    const MInt gridCellId = tree().solver2grid(cellId);
    return raw().treeb().hasProperty(gridCellId, Cell::IsPeriodic);
  }
  MInt periodicCartesianDir(const MInt dir) const { return raw().periodicCartesianDir(dir); }
  MFloat periodicCartesianLength(const MInt dir) const { return raw().periodicCartesianLength(dir); }
  MInt maxNoCells() const { return m_maxNoCells; }
  MFloat gridCellVolume(const MInt level) const { return raw().gridCellVolume(level); }
  MBool azimuthalPeriodicity() const { return raw().m_azimuthalPer; }
  MInt determineAzimuthalBoundarySide(const MFloat* coords);
  MInt azimuthalDir(MInt dir) { return raw().m_azimuthalPeriodicDir[dir]; }
  MFloat azimuthalCenter() { return raw().m_azimuthalBbox.azimuthalCenter(); }
  MFloat azimuthalAngle() { return raw().m_azimuthalAngle; }
  void rotateCartesianCoordinates(MFloat* coords, MFloat angle) { raw().rotateCartesianCoordinates(coords, angle); }
  MInt neighborList(const MInt cellId, const MInt dir) const { return m_neighborList[cellId][dir]; }
  MInt& neighborList(const MInt cellId, const MInt dir) { return m_neighborList[cellId][dir]; }
  MLong localPartitionCellOffsets(const MInt index) const { return raw().localPartitionCellOffsets(index); }
  MInt noLocalPartitionCells() { return localPartitionCellOffsets(1) - localPartitionCellOffsets(0); }

  MInt localPartitionCellLocalIds(const MInt id) const {
    return m_tree.grid2solver(raw().localPartitionCellLocalIds(id));
  }
  MLong localPartitionCellGlobalIds(const MInt id) const {
    // NOTE: solver based globalId
    return localPartitionCellLocalIds(id) > -1 ? m_tree.globalId(localPartitionCellLocalIds(id)) : -1;
  }

  // restart-file versions
  MLong localPartitionCellOffsetsRestart(const MInt index) const {
    /*if(!raw().updatedPartitionCells() && !hasInactiveRanks()) {
      return raw().m_localPartitionCellOffsetsRestart[index];
    } else if (!hasInactiveRanks()) {
      return localPartitionCellOffsets(index);
    } else {*/
    return m_localPartitionCellOffsets[index];
    //}
  }
  MLong localPartitionCellGlobalIdsRestart(const MInt id) const {
    /*if(!raw().updatedPartitionCells() && !hasInactiveRanks()) {
      return m_tree.grid2solver(raw().m_localPartitionCellLocalIdsRestart[id]) > -1 ?
             m_tree.globalId(m_tree.grid2solver(raw().m_localPartitionCellLocalIdsRestart[id])): -1;
    } else if (!hasInactiveRanks()) {
      return localPartitionCellGlobalIds(id);
    } else {*/
    return m_partitionCellGlobalId[id];
    //}
  }

  MBool isPeriodic(const MInt cellId) const {
    return raw().a_hasProperty(m_tree.solver2grid(cellId), Cell::IsPeriodic);
  }
  MBool solverFlag(const MInt gridId, const MInt solverId) const { return raw().treeb().solver(gridId, solverId); }

  MLong generateHilbertIndex(const MInt cellId, const MInt refLevel = -1) {
    const MInt level = refLevel > -1 ? refLevel : minLevel();
    return raw().generateHilbertIndex(m_tree.solver2grid(cellId), level);
  }

  MInt domainIndex(const MInt id) {
    if(!g_multiSolverGrid && noDomains() > 1) {
      ASSERT(m_neighborDomainIndex[id] == raw().m_nghbrDomainIndex[id], "");
    }
    return m_neighborDomainIndex[id];
  }
  MInt azimuthalDomainIndex(const MInt id) { return m_azimuthalNeighborDomainIndex[id]; }

  MInt a_storeNghbrIds(const MInt id) const { return m_storeNghbrIds[id]; }
  MInt a_identNghbrIds(const MInt id) const { return m_identNghbrIds[id]; }
  MInt a_neighborList(const MInt cellId, const MInt dir) const { return m_neighborList[cellId][dir]; }

  MBool checkNghbrIds() { return m_storeNghbrIds != nullptr && m_identNghbrIds != nullptr; }

  // constexpr MInt noHaloLayers() const { return raw().noHaloLayers(); }
  constexpr MInt noHaloLayers() const { return raw().noHaloLayers(solverId()); }

  // Methods modifying the tree proxy directly (use with care!)
  void setSolver2grid(const MInt solverId, const MInt treeId) { m_tree.setSolver2grid(solverId, treeId); }
  void setGrid2solver(const MInt treeId, const MInt solverId) { m_tree.setGrid2solver(treeId, solverId); }
  void swapSolverIds(const MInt id0, const MInt id1) { m_tree.swapSolverIds(id0, id1); }
  void swapGridIds(const MInt id0, const MInt id1) { m_tree.swapGridIds(id0, id1); }


  /** \brief Exchange Halo/Window data in a general way for ParaView or other external visuliazation
   ** software
   **
   ** @author: Pascal Meysonnat
   ** @date: May 1984
   **/
  template <typename U>
  void exchangeHaloCellsForVisualization(U* data) {
    if(noNeighborDomains() > 0) {
      maia::mpi::exchangeData(m_neighborDomains, m_haloCells, m_windowCells, mpiComm(), data);
    }
  }
  void exchangeHaloCellsForVisualizationDG(MFloat* data, const MInt* const polyDegs, const MInt* const dataOffsets);
  void exchangeHaloCellsForVisualizationSBP(MFloat* data, const MInt* const noNodes1D, const MInt* const dataOffsets);

  // Check wether cell is located outside geometry
  // Necessary for azimuthal periodicity to check if azimuthal halo cell is requiered
  MBool checkOutsideGeometry(MInt gridId);

 private:
  // Data
  /// Solver id
  const MInt m_solverId = -1;

  static constexpr const MInt m_noDirs = 2 * nDim;

  /// Reference to actual grid
  Grid& m_grid;

  /// Tree proxy object
  TreeProxy m_tree;

  /// Reference to solver geometry
  Geom* m_geometry;

  // Parallelization-related data
  MInt m_noInternalCells = -1;
  MPI_Comm m_mpiComm = MPI_COMM_NULL;
  MInt m_domainId = -1;
  MInt m_noDomains = -1;
  MBool m_isActive = false;
  MBool m_hasInactiveRanks = false;
  std::vector<MLong> m_domainOffsets;

  // Neighbor/window/halo data
  std::vector<MInt> m_neighborDomains;
  std::vector<MInt> m_neighborDomainIndex;
  std::vector<std::vector<MInt>> m_windowCells;
  std::vector<std::vector<MInt>> m_haloCells;
  std::vector<MInt> m_leafSendNeighborDomains;
  std::vector<MInt> m_leafRecvNeighborDomains;
  std::vector<std::vector<MInt>> m_leafWindowCells;
  std::vector<std::vector<MInt>> m_leafHaloCells;
  MInt m_leafRecvSize = -1;
  MInt m_leafSendSize = -1;
  std::map<MInt, MInt> m_global2solver;

  // Azimuthal periodicity window/halo data
  std::vector<MInt> m_azimuthalNeighborDomains;
  std::vector<MInt> m_azimuthalNeighborDomainIndex;
  std::vector<std::vector<MInt>> m_azimuthalWindowCells;
  std::vector<std::vector<MInt>> m_azimuthalHaloCells;
  std::vector<MInt> m_azimuthalUnmappedHaloCells;
  std::vector<MInt> m_azimuthalUnmappedHaloDomains;
  std::vector<MInt> m_isOutsideHalo;
  std::vector<MInt> m_isOutsideWindow;

  // Other grid information
  MInt m_maxLevel = -1;
  MInt m_maxRefinementLevel = -1;
  MInt m_maxNoCells = 0;

  // Holds the reverse-direction-Information
  const MInt m_revDir[6] = {1, 0, 3, 2, 5, 4};

  MInt* m_storeNghbrIds = nullptr;
  MInt* m_identNghbrIds = nullptr;
  MInt** m_neighborList{};

  MLong* m_partitionCellGlobalId = nullptr;
  MLong m_localPartitionCellOffsets[3] = {static_cast<MLong>(-1), static_cast<MLong>(-1), static_cast<MLong>(-1)};
};


namespace tree {

/// Class that acts like the original grid tree but only for a single solver
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-16
template <MInt nDim>
class TreeProxy {
  // Friends & family
  // Grid proxy must be declared as template friend as otherwise classes deriving from this tree
  // proxy do not work
  template <MInt nDim_>
  friend class maia::grid::Proxy;

 public:
  // Types
  using Tree = maia::grid::tree::Tree<nDim>;
  using Cell = maia::grid::tree::Cell;

  /// Constructor does nothing but pass arguments to members
  TreeProxy(const MInt solverId_, Tree& tree_) : m_tree(tree_), m_solverId(solverId_) {
    m_solver2grid.reserve(m_tree.capacity());
    m_grid2solver.reserve(m_tree.capacity());
    m_globalIds.reserve(m_tree.capacity());
    m_isCutOffCell.reserve(m_tree.capacity());
  }

  /// Convert solver cell id to grid cell id
  MInt solver2grid(const MInt id) const {
    // TODO labels:GRID c++14
    //#ifndef NDEBUG
    //      ASSERT(id >= 0 && id < m_solver2grid.size(), "id in solver2grid is out of bounds");
    //#endif
    return m_solver2grid[id];
  }
  /// Convert grid cell id to solver cell id (+1 in size is to map "-1" to "-1")
  MInt grid2solver(const MInt id) const {
    // TODO labels:GRID c++14
    //#ifndef NDEBUG
    //    ASSERT(id+1 >= 0 && id+1 < m_grid2solver.size(), "id in grid2solver is out of bounds");
    //#endif
    return m_grid2solver[id + 1];
  }

  void setSolver2grid(const MInt solverId, const MInt treeId) {
    ASSERT(solverId > -1 && solverId < (signed)m_solver2grid.size(),
           "id, size: " << solverId << " " << m_solver2grid.size());
    m_solver2grid[solverId] = treeId;
  }
  void setGrid2solver(const MInt treeId, const MInt solverId) {
    ASSERT((treeId + 1) > -1 && (treeId + 1) < (signed)m_grid2solver.size(),
           "id size: " << treeId + 1 << " " << m_grid2solver.size());
    m_grid2solver[treeId + 1] = solverId;
  }

  void swapSolverIds(const MInt id0, const MInt id1) {
    ASSERT(id0 > -1 && id1 > -1 && id0 < size() && id1 < size(),
           "id0, id1, size: " << id0 << " , " << id1 << " , " << size());
    const MInt g0 = m_solver2grid[id0];
    const MInt g1 = m_solver2grid[id1];
    if(g0 > -1) m_grid2solver[g0 + 1] = id1;
    if(g1 > -1) m_grid2solver[g1 + 1] = id0;
    std::swap(m_solver2grid[id0], m_solver2grid[id1]);
  }

  void swapGridIds(const MInt id0, const MInt id1) {
    ASSERT(id0 > -1 && id1 > -1 && id0 < m_tree.size() && id1 < m_tree.size()
               && (id0 + 1) < (signed)m_grid2solver.size() && (id1 + 1) < (signed)m_grid2solver.size(),
           "");
    const MInt b0 = m_grid2solver[id0 + 1];
    const MInt b1 = m_grid2solver[id1 + 1];
    if(b0 > -1) m_solver2grid[b0] = id1;
    if(b1 > -1) m_solver2grid[b1] = id0;
    std::swap(m_grid2solver[id0 + 1], m_grid2solver[id1 + 1]);
  }

  /// Return tree size (= number of cells)
  MInt size() const { return m_solver2grid.size(); }

  // Parent-child relationship
  MInt parent(const MInt id) const;
  MBool hasParent(const MInt id) const;
  MInt child(const MInt id, const MInt pos) const;
  MBool hasChild(const MInt id, const MInt pos) const;
  MBool hasChildren(const MInt id) const;
  MInt noChildren(const MInt id) const;
  MBool isLeafCell(const MInt id) const;

  // Neighbors
  MInt neighbor(const MInt id, const MInt dir) const;
  MBool hasNeighbor(const MInt id, const MInt dir) const;
  MBool hasAnyNeighbor(const MInt id, const MInt dir) const;

  // // Other data fields
  MLong globalId(const MInt id) const;
  MInt localId(const MLong id) const;
  MInt cutOff(const MInt id) const;
  MInt level(const MInt id) const;
  const MFloat& coordinate(const MInt id, const MInt dir) const;
  // MFloat& weight(const MInt id);
  MFloat weight(const MInt id) const;
  MBool hasProperty(const MInt id, const Cell p) const;
  // MUlong propertiesToBits(const MInt id) const;
  // MString propertiesToString(const MInt id) const;
  // auto properties(const MInt id) -> decltype(tree().properties(solver2grid(id))) & ;

  // // Other data fields (subject to change)
  // MInt& noOffsprings(const MInt id);
  // MInt noOffsprings(const MInt id) const;
  // MFloat& workload(const MInt id);
  // MFloat workload(const MInt id) const;

  // Entries per tree node
  /// Return maximum number of children per node
  static constexpr MInt noChildrenPerNode() { return 4 * nDim - 4; }
  /// Return maximum number of same-level neighbors per node
  static constexpr MInt noNeighborsPerNode() { return 2 * nDim; }

  /// Return opposite direction for given neighbor direction
  static constexpr MInt oppositeNeighborDir(const MInt dir) { return dir + 1 - 2 * (dir % 2); }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::grid::cell::p(Cell::NumProperties); }

 private:
  /// Return solver id
  MInt solverId() const { return m_solverId; }

  /// Default implementation does nothing as everything relevant is set by grid proxy class
  void update(const Proxy<nDim>& NotUsed(gridProxy_)) {}

  /// Create global to local id mapping
  void createGlobalToLocalIdMapping();

  // Data
  /// Reference to actual grid tree
  Tree& m_tree;

  /// Solver id
  const MInt m_solverId = -1;

  // Grid cell <-> solver cell mappings
  std::vector<MInt> m_solver2grid;
  std::vector<MInt> m_grid2solver;
  std::vector<MLong> m_globalIds;
  std::map<MLong, MInt> m_globalToLocalId{};
  std::vector<MInt> m_isCutOffCell;
};


/// Accessor for parent node.
template <MInt nDim>
MInt TreeProxy<nDim>::parent(const MInt id) const {
  return grid2solver(m_tree.parent(solver2grid(id)));
}


/// Return whether node has parent.
template <MInt nDim>
MBool TreeProxy<nDim>::hasParent(const MInt id) const {
  return parent(id) > -1;
}


/// Accessor for child node.
template <MInt nDim>
MInt TreeProxy<nDim>::child(const MInt id, const MInt pos) const {
  return grid2solver(m_tree.child(solver2grid(id), pos));
}


/// Return whether node has child at given position.
template <MInt nDim>
MBool TreeProxy<nDim>::hasChild(const MInt id, const MInt pos) const {
  return child(id, pos) > -1;
}


/// Return whether node has any children.
template <MInt nDim>
MBool TreeProxy<nDim>::hasChildren(const MInt id) const {
  return noChildren(id) > 0;
}


/// Return number of children of given node.
template <MInt nDim>
MInt TreeProxy<nDim>::noChildren(const MInt id) const {
  MInt count = 0;
  for(MInt pos = 0; pos < noChildrenPerNode(); pos++) {
    count += hasChild(id, pos);
  }
  return count;
}

/// Accessor for isLeafCell
template <MInt nDim>
MBool TreeProxy<nDim>::isLeafCell(const MInt id) const {
  return m_tree.isLeafCell(solver2grid(id), solverId());
}


/// Accessor for neighbor node.
template <MInt nDim>
MInt TreeProxy<nDim>::neighbor(const MInt id, const MInt dir) const {
  return grid2solver(m_tree.neighbor(solver2grid(id), dir));
}


/// Return whether node has same-level neighbor in given direction.
template <MInt nDim>
MBool TreeProxy<nDim>::hasNeighbor(const MInt id, const MInt dir) const {
  return neighbor(id, dir) > -1;
}
/// Return whether node or its parent has neighbor in given direction.
template <MInt nDim>
MBool TreeProxy<nDim>::hasAnyNeighbor(const MInt id, const MInt dir) const {
  return hasNeighbor(id, dir) || (hasParent(id) && hasNeighbor(parent(id), dir));
}


/// Accessor for global id.
template <MInt nDim>
MLong TreeProxy<nDim>::globalId(const MInt id) const {
  return m_globalIds[id];
}

/// Accessor for cutOff dir.
template <MInt nDim>
MInt TreeProxy<nDim>::cutOff(const MInt id) const {
  return m_isCutOffCell[id];
}

/// Accessor for level.
template <MInt nDim>
MInt TreeProxy<nDim>::level(const MInt id) const {
  return m_tree.level(solver2grid(id));
}


/// Accessor for coordinates.
template <MInt nDim>
const MFloat& TreeProxy<nDim>::coordinate(const MInt id, const MInt dir) const {
  return m_tree.coordinate(solver2grid(id), dir);
}


// /// Accessor for weight.
// template <MInt nDim> MFloat& TreeProxy<nDim>::weight(const MInt id) {
//   return m_tree.weight(solver2grid(id));
// }
/// Accessor for weight (const version).
template <MInt nDim>
MFloat TreeProxy<nDim>::weight(const MInt id) const {
  return m_tree.weight(solver2grid(id));
}
//
//
// /// Accessor for noOffsprings.
// template <MInt nDim> MInt& TreeProxy<nDim>::noOffsprings(const MInt id) {
//   return m_tree.noOffsprings(solver2grid(id));
// }
// /// Accessor for noOffsprings (const version).
// template <MInt nDim> MInt TreeProxy<nDim>::noOffsprings(const MInt id) const {
//   return m_tree.noOffsprings(solver2grid(id));
// }
//
//
// /// Accessor for workload.
// template <MInt nDim> MFloat& TreeProxy<nDim>::workload(const MInt id) {
//   return m_tree.workload(solver2grid(id));
// }
// /// Accessor for workload (const version).
// template <MInt nDim> MFloat TreeProxy<nDim>::workload(const MInt id) const {
//   return m_tree.workload(solver2grid(id));
// }


/// Accessor for properties.
template <MInt nDim>
MBool TreeProxy<nDim>::hasProperty(const MInt id, const Cell p) const {
  ASSERT(p != Cell::IsWindow, "Query for IsWindow won't give solver-specific information!");
  return m_tree.hasProperty(solver2grid(id), p);
}

// /// Convert properties to string.
// template <MInt nDim> MString TreeProxy<nDim>::propertiesToString(const MInt id) const {
//   return m_tree.propertiesToString(solver2grid(id));
// }
// /// Accessor for properties.
// template <MInt nDim>
// auto TreeProxy<nDim>::properties(const MInt id) ->
// decltype(m_tree.properties(solver2grid(id))) & {
//   return m_tree.properties(solver2grid(id));
// }

/// \brief Create the mapping from global to local cell ids
template <MInt nDim>
void TreeProxy<nDim>::createGlobalToLocalIdMapping() {
  TRACE();

  std::map<MLong, MInt>().swap(m_globalToLocalId);
  std::set<MInt> periodicCells;
  for(MInt cellId = 0; cellId < size(); cellId++) {
    if(hasProperty(cellId, Cell::IsPeriodic)) {
      periodicCells.insert(cellId);
      continue;
    }
    ASSERT(m_globalToLocalId.count(globalId(cellId)) == 0, "Global id already in map");
    m_globalToLocalId[globalId(cellId)] = cellId;
  }

  for(auto it = periodicCells.begin(); it != periodicCells.end(); it++) {
    if(m_globalToLocalId.count(globalId(*it)) == 0) {
      m_globalToLocalId[globalId(*it)] = *it;
    }
  }
}

/// \brief Accesor for local id
template <MInt nDim>
MInt TreeProxy<nDim>::localId(const MLong id) const {
  TRACE();

  if(m_globalToLocalId.count(id) > 0) {
    return m_globalToLocalId.at(id);
  } else {
    return -1;
  }
}

} // namespace tree

} // namespace grid
} // namespace maia

#endif // ifndef GRIDPROXY_H_
