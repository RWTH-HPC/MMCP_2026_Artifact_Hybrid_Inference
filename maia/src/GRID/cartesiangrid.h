// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CARTESIANGRID_H
#define CARTESIANGRID_H
//#define IMPAIR_SUBDIVISION

#include <algorithm>
#include <array>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/timer.h"
#include "cartesiangridtree.h"

template <class T>
class Collector;

// Forward declaration of cartesiangridio.h:
namespace maia {
namespace grid {
template <typename Grid>
class IO;
template <MInt nDim>
class Controller;
} // namespace grid
} // namespace maia

/** Template class for the operations and data structures of the grid
 */
template <MInt nDim>
class CartesianGrid {
  /// [Splitt] The following is part of a first step to splitt CartesianGrid
  /// from the inheritance hierarchy:
  ///
  /// \todo labels:GRID the following friend declarations will be removed in a future commit
  template <MInt nDim_, class SysEqn>
  friend class FvCartesianSolverXD;
  template <MInt nDim_>
  friend class LsCartesianSolver;
  template <MInt nDim_, class SysEqn>
  friend class FvCartesianSolverXD;
  template <MInt nDim_, class SysEqn>
  friend class FvBndryCnd2D;
  template <MInt nDim_, class SysEqn>
  friend class FvBndryCnd3D;
  template <MInt nDim_, class SysEqn>
  friend class FvMbCartesianSolverXD;
  template <MInt nDim_>
  friend class LbSolver;
#if not defined(MAIA_MS_COMPILER)
  template <MInt nDim_, class SysEqn>
  friend class DgCartesianSolver;
#endif
  friend class maia::grid::Controller<nDim>;

 public:
  // cartesiangridio.h API:
  friend class maia::grid::IO<CartesianGrid<nDim>>;
  static constexpr MInt nDim_() { return nDim; }
  using MInt_dyn_array = MInt*;
  using Tree = maia::grid::tree::Tree<nDim>;

  // Type for cell properties
  using Cell = typename Tree::Cell;

  /// \brief Returns weather the grid has a cutOff
  MBool hasCutOff() const { return m_cutOff; }

  // NOTE: the cartesiangrid requires parent/child/neighbor relations as Longs
  //      during initialisation and balance of the grid,
  //      localToGlobalIds and changeGlobalToLocalIds switches between local and global storage!
  /// \brief Returns the parent of the cell \p cellId
  MLong& a_parentId(const MInt cellId) { return m_tree.parent(cellId); }

  /// \brief Returns the parent of the cell \p cellId
  MLong a_parentId(const MInt cellId) const { return m_tree.parent(cellId); }

  /// \brief Returns the parent of the cell \p cellId
  MLong a_parentId(const MInt cellId, const MInt solverId) const {
    return (m_tree.parent(cellId) > -1 && m_tree.solver(m_tree.parent(cellId), solverId)) ? m_tree.parent(cellId) : -1;
  }

  /// \brief Returns if the cell \p cellId has a parent
  MBool a_hasParent(const MInt cellId) const { return m_tree.hasParent(cellId); }

  /// \brief Returns if the cell \p cellId has a parent
  MBool a_hasParent(const MInt cellId, const MInt solverId) const {
    return m_tree.hasParent(cellId) && m_tree.solver(m_tree.parent(cellId), solverId);
  }

  /// \brief Returns the child id of the cell \p cellId at position \p pos
  MLong& a_childId(const MInt cellId, const MInt pos) { return m_tree.child(cellId, pos); }

  /// \brief Returns the child id of the cell \p cellId at position \p pos
  MLong a_childId(const MInt cellId, const MInt pos) const { return m_tree.child(cellId, pos); }

  /// \brief Returns the child id of the cell \p cellId at position \p pos
  MLong a_childId(const MInt cellId, const MInt pos, const MInt solverId) const {
    return (m_tree.child(cellId, pos) > -1 && m_tree.solver(m_tree.child(cellId, pos), solverId))
               ? m_tree.child(cellId, pos)
               : -1;
  }

  MLong noPartitionCellsGlobal() const { return m_noPartitionCellsGlobal; }

  /// \brief Returns the level of the cell \p cellId
  MInt& a_level(const MInt cellId) { return m_tree.level(cellId); }

  /// \brief Returns the level of the cell \p cellId
  MInt a_level(const MInt cellId) const { return m_tree.level(cellId); }

  /// \brief Returns the workload of the cell \p cellId
  MFloat a_workload(const MInt cellId) const { return m_tree.workload(cellId); }

  /// \brief Returns the workload of the cell \p cellId
  MFloat& a_workload(const MInt cellId) { return m_tree.workload(cellId); }

  /// \brief Returns the noOffsprings of the cell \p
  MInt a_noOffsprings(const MInt cellId) const { return m_tree.noOffsprings(cellId); }

  /// \brief Returns the noOffsprings of the cell \p
  MInt& a_noOffsprings(const MInt cellId) { return m_tree.noOffsprings(cellId); }

  /// \brief Returns the globalId of the cell \p cellId in collector \p cells_
  MLong& a_globalId(const MInt cellId) { return m_tree.globalId(cellId); }

  /// \brief Returns the globalId of the cell \p cellId in collector \p cells_
  MLong a_globalId(const MInt cellId) const { return m_tree.globalId(cellId); }

  /// \brief Returns the weight of the cell \p cellId
  MFloat& a_weight(const MInt cellId) { return m_tree.weight(cellId); }

  /// \brief Returns the weight of the cell \p cellId
  MFloat a_weight(const MInt cellId) const { return m_tree.weight(cellId); }

  /// \brief Returns the no. of children of the cell \p cellId
  MInt a_noChildren(const MInt cellId) const { return m_tree.noChildren(cellId); }

  /// \brief Returns if the cell \p cellId has children
  MInt a_hasChildren(const MInt cellId) const { return m_tree.hasChildren(cellId); }

  /// \brief Returns if the cell \p cellId has children
  MInt a_hasChildren(const MInt cellId, const MInt solverId) const { return m_tree.hasChildren(cellId, solverId); }

  /// \brief Returns if the cell \p cellId has a child at position \p pos
  MInt a_hasChild(const MInt cellId, const MInt pos) const { return m_tree.hasChild(cellId, pos); }

  /// \brief Returns if the cell \p cellId has a child at position \p pos
  MInt a_hasChild(const MInt cellId, const MInt pos, const MInt solverId) const {
    return m_tree.hasChild(cellId, pos) && m_tree.solver(m_tree.child(cellId, pos), solverId);
  }

  /// \brief Returns the delete of the cell \p cellId
  maia::grid::cell::BitsetType::reference a_isToDelete(const MInt cellId) {
    return m_tree.hasProperty(cellId, Cell::IsToDelete);
  }

  /// \brief Returns the delete of the cell \p cellId
  MBool a_isToDelete(const MInt cellId) const { return m_tree.hasProperty(cellId, Cell::IsToDelete); }

  /// \brief Returns the neighbor id of the cell \p cellId \p dir
  MLong& a_neighborId(const MInt cellId, const MInt dir) { return m_tree.neighbor(cellId, dir); }

  /// \brief Returns the neighbor id of the cell \p cellId \p dir
  MLong a_neighborId(const MInt cellId, const MInt dir) const { return m_tree.neighbor(cellId, dir); }

  /// \brief Returns the neighbor id of the cell \p cellId \p dir
  MLong a_neighborId(const MInt cellId, const MInt dir, const MInt solverId) const {
    return (m_tree.neighbor(cellId, dir) > -1 && m_tree.solver(m_tree.neighbor(cellId, dir), solverId))
               ? m_tree.neighbor(cellId, dir)
               : -1;
  }

  /// \brief Returns the coordinate of the cell \p cellId for direction \p dir
  MFloat& a_coordinate(const MInt cellId, const MInt dir) { return m_tree.coordinate(cellId, dir); }

  /// \brief Returns the coordinate of the cell \p cellId for direction \p dir
  MFloat a_coordinate(const MInt cellId, const MInt dir) const { return m_tree.coordinate(cellId, dir); }

  /// \brief Returns the coordinate of the cell \p cellId in collector \p cells_ for direction \p dir
  template <class U>
  MFloat a_coordinate(Collector<U>* NotUsed(cells_), const MInt cellId, const MInt dir) const {
    return a_coordinate(cellId, dir);
  }


  /// \brief Returns whether cell is leaf \p cellId
  MBool a_isLeafCell(const MInt cellId) const { return m_tree.isLeafCell(cellId); }

  /// \brief Returns whether cell is leaf \p cellId
  maia::grid::tree::SolverBitsetType::reference a_isLeafCell(const MInt cellId, const MInt solverId) {
    return m_tree.isLeafCell(cellId, solverId);
  }

  /// \brief Returns whether cell is leaf \p cellId
  MBool a_isLeafCell(const MInt cellId, const MInt solverId) const { return m_tree.isLeafCell(cellId, solverId); }

  /// \brief Returns noNeighborIds of the cell \p CellId
  MBool a_hasNeighbor(const MInt cellId, const MInt dir) const { return m_tree.hasNeighbor(cellId, dir); }

  /// \brief Returns noNeighborIds of the cell \p CellId
  MBool a_hasNeighbor(const MInt cellId, const MInt dir, const MInt solverId) const {
    if(treeb().noSolvers() == 1 && m_tree.hasNeighbor(cellId, dir)
       && !m_tree.solver(m_tree.neighbor(cellId, dir), solverId)) {
      std::cerr << "mis " << cellId << " " << dir << " " << m_tree.neighbor(cellId, dir) << " " << solverId << " "
                << m_tree.solver(m_tree.neighbor(cellId, dir), solverId) << std::endl;
      ASSERT(false, "");
    }
    return m_tree.hasNeighbor(cellId, dir) && m_tree.solver(m_tree.neighbor(cellId, dir), solverId);
  }

  /// \brief Returns property \p p of the cell \p cellId
  void a_resetProperties(const MInt cellId) { m_tree.resetProperties(cellId); }

  /// \brief Returns property \p p of the cell \p cellId
  void a_copyProperties(const MInt fromCellId, const MInt toCellId) {
    m_tree.properties(toCellId) = m_tree.properties(fromCellId);
  }


  /// \brief Returns property \p p of the cell \p cellId
  maia::grid::cell::BitsetType::reference a_hasProperty(const MInt cellId, const Cell p) {
    return m_tree.hasProperty(cellId, p);
  }

  /// \brief Returns property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const Cell p) const { return m_tree.hasProperty(cellId, p); }

  MString a_propertiesToString(const MInt cellId) { return m_tree.propertiesToString(cellId); }

  /// \brief Returns property \p p of the cell \p cellId
  MBool a_isHalo(const MInt cellId) const { return m_tree.hasProperty(cellId, Cell::IsHalo); }

  /// \brief Returns if cell \p cellId is used by solver \p solverId
  maia::grid::tree::SolverBitsetType::reference a_solver(const MInt cellId, const MInt solverId) {
    return m_tree.solver(cellId, solverId);
  }

  /// \brief Sets if cell is used by solver
  void setSolver(const MInt cellId, const MInt solverId, MBool flag) { m_tree.solver(cellId, solverId) = flag; }

  MBool addSolverToGrid() { return m_addSolverToGrid; }

  MBool wasAdapted() const { return m_wasAdapted; }
  MBool wasBalanced() const { return m_wasBalanced; }
  MBool wasBalancedAtLeastOnce() const { return m_wasBalancedAtLeastOnce; }

  MBool m_wasAdapted = false;
  MBool m_wasBalanced = false;
  MBool m_wasBalancedAtLeastOnce = false;

  MInt m_noPeriodicCartesianDirs = 0;
  std::array<MInt, 3> m_periodicCartesianDir{};
  MFloat m_periodicCartesianLength[3]{};
  MBool m_azimuthalPer = false;
  std::array<MFloat, 3> m_azimuthalPerCenter{};
  std::array<MInt, 2> m_azimuthalPeriodicDir{};
  MInt m_azimuthalAxialDir = -1;
  MFloat m_azimuthalAngle = F0;
  struct azimuthalBbox {
    MFloat azimuthalAngle = F0;
    MFloat minCoord[nDim];
    MFloat boxLength[nDim];
    MInt nBins[nDim];
    MFloat* minPer3 = nullptr;
    MFloat* maxPer3 = nullptr;
    MFloat center = F0;
    MInt axDir = -1;
    MInt perDir1 = -1;
    MInt perDir2 = -1;
    MInt noHaloLayers = -1;

    void init(MFloat angle, MInt periodicDir1, MInt periodicDir2, MInt noLayers) {
      azimuthalAngle = angle;
      perDir1 = periodicDir1;
      perDir2 = periodicDir2;
      noHaloLayers = noLayers;
      MInt perSize1 = nBins[perDir1];
      MInt perSize2 = nBins[perDir2];
      mDeallocate(minPer3);
      mDeallocate(maxPer3);
      mAlloc(minPer3, perSize1 * perSize2, "minPer3", std::numeric_limits<MFloat>::max(), AT_);
      mAlloc(maxPer3, perSize1 * perSize2, "maxPer3", -std::numeric_limits<MFloat>::max(), AT_);
    }
    MInt azimuthalSide(const MFloat phi) {
      MFloat mid = center;
      if(mid > phi) {
        return -1;
      } else {
        return 1;
      }
    }
    MFloat azimuthalCenter() { return center; }
    MFloat azimuthalBoundary(const MFloat* coords, MInt side) {
      MInt perSize2 = nBins[perDir2];
      MInt ind1 = index(coords, perDir1);
      MInt ind2 = index(coords, perDir2);
      if(side < 0) {
        MFloat minPhi = center - F1B2 * azimuthalAngle; // std::numeric_limits<MFloat>::max();
        for(MInt i = -(noHaloLayers + 1); i < (noHaloLayers + 2); i++) {
          MInt ind1_ = ind1 + i;
          ind1_ = mMax(0, mMin(ind1_, nBins[perDir1] - 1));
          MInt ind2_ = ind2 + i;
          ind2_ = mMax(0, mMin(ind2_, nBins[perDir2] - 1));
          minPhi = mMin(minPhi, minPer3[ind1 * perSize2 + ind2_]);
          minPhi = mMin(minPhi, minPer3[ind1_ * perSize2 + ind2]);
        }
        return minPhi;
      } else if(side > 0) {
        MFloat maxPhi = center + F1B2 * azimuthalAngle; //-std::numeric_limits<MFloat>::max();
        for(MInt i = -(noHaloLayers + 1); i < (noHaloLayers + 2); i++) {
          MInt ind1_ = ind1 + i;
          ind1_ = mMax(0, mMin(ind1_, nBins[perDir1] - 1));
          MInt ind2_ = ind2 + i;
          ind2_ = mMax(0, mMin(ind2_, nBins[perDir2] - 1));
          maxPhi = mMax(maxPhi, maxPer3[ind1 * perSize2 + ind2_]);
          maxPhi = mMax(maxPhi, maxPer3[ind1_ * perSize2 + ind2]);
        }
        return maxPhi;
      } else {
        mTerm(1, AT_, "side = 0?");
      }
    }
    MInt index(const MFloat* coord, MInt dir) {
      MInt ind = -1;
      ind = (MInt)((coord[dir] - minCoord[dir]) * nBins[dir] / boxLength[dir]);
      if(ind < 0)
        ind = 0;
      else if(ind >= nBins[dir])
        ind = nBins[dir] - 1;

      return ind;
    }
  };
  azimuthalBbox m_azimuthalBbox;

  MLong m_noCellsGlobal{};

  // Marks the beginning of halo cells in the collector
  MInt m_noInternalCells{};

  MInt m_noDomains{};

  MInt m_maxUniformRefinementLevel = -1;
  MInt m_maxRfnmntLvl = -1;
  //! The max Level of refinement of the solver length
  MInt m_maxLevel = -1;
  //! The min Level of refinement of the solver length
  MInt m_minLevel = -1;
  //! The new min Level which has been increased at the restart!
  MInt m_newMinLevel = -1;

  MBool m_allowInterfaceRefinement = false;
  MInt m_cutOff = 0;
  MBool m_lowMemAdaptation = true;

  // Grid cell volume
  MFloat* m_gridCellVolume = nullptr;

  // File processing
  MString m_gridInputFileName = "";

  template <class F>
  MInt reduceToLevel(const MInt reductionLevel, F interpolateToParentCells);

  /// \brief Returns cell length at cell level \p level.
  ///
  /// Note cell length h = h_0 / 2^level where h_0 is the length of the grid's
  /// root cell.
  MFloat cellLengthAtLevel(const MInt level) const { return m_lengthLevel0 * FFPOW2(level); }
  /// \brief Returns the cell length of cell \p cellId.
  MFloat cellLengthAtCell(const MInt cellId) const { return cellLengthAtLevel(a_level(cellId)); }
  /// \brief Returns the half cell length of cell \p cellId.
  MFloat halfCellLength(const MInt cellId) const { return F1B2 * cellLengthAtCell(cellId); }
  MFloat cellVolumeAtLevel(const MInt level) const {
    return (nDim == 3) ? POW3(cellLengthAtLevel(level)) : POW2(cellLengthAtLevel(level));
  }

  /// Store the grid bounding box.
  ///
  /// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
  /// \date 2015-03-10
  ///
  /// \param[out] box Location where bounding box information is to be stored.
  ///                 Must be at least of size 2 * dim.
  ///
  /// The bounding box is stored as follows:
  /// 2D: min_x | min_y | max_x | max_y
  /// 3D: min_x | min_y | min_z | max_x | max_y | max_z
  void boundingBox(MFloat* const box) const { std::copy_n(&m_boundingBox[0], 2 * nDim, box); }

  /// return global Bounding box
  const MFloat* globalBoundingBox() const { return &m_boundingBox[0]; }

  /// Store center of gravity coordinates.
  ///
  /// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
  /// \date 2015-03-10
  ///
  /// \param[out] center Location where coordinates are to be stored. Must be at
  ///                    least of size nDim.
  void centerOfGravity(MFloat* const center) const { std::copy_n(m_centerOfGravity, nDim, center); }

  /// \brief Return the length of the level 0 cell.
  ///
  /// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
  /// \date 2015-03-10
  ///
  /// \return Length of the cell at level 0.
  MFloat lengthLevel0() const { return m_lengthLevel0; }

  /// Return number of internal cells (i.e., excluding halo cells)
  MInt noInternalCells() const { return m_noInternalCells; }

  /// Return number of neighbor domains
  MInt noNeighborDomains() const { return (signed)m_nghbrDomains.size(); }

  /// Return neighbor domain
  const MInt& neighborDomain(const MInt id) const { return m_nghbrDomains[id]; }

  /// Return vector of neighbor domains
  const std::vector<MInt>& neighborDomains() const { return m_nghbrDomains; }

  /// Return domain offset
  const MLong& domainOffset(const MInt id) const { return m_domainOffsets[id]; }

  /// Return number of halo cells for given domain
  MInt noHaloCells(const MInt domainId) const { return m_haloCells[domainId].size(); }
  constexpr MBool paraViewPlugin() const { return m_paraViewPlugin; }

  /// Return halo cell id
  MInt haloCell(const MInt domainId, const MInt cellId) const { return m_haloCells[domainId][cellId]; }

  // Returns vector of halo cells
  const std::vector<std::vector<MInt>>& haloCells() const { return m_haloCells; };

  /// Return number of window cells for given domain
  MInt noWindowCells(const MInt domainId) const { return m_windowCells[domainId].size(); }

  /// Return window cell id
  MInt windowCell(const MInt domainId, const MInt cellId) const { return m_windowCells[domainId][cellId]; }

  // Returns vector of window cells
  const std::vector<std::vector<MInt>>& windowCells() const { return m_windowCells; };

  //
  MInt windowLayer(const MInt domainId, const MInt cellId, const MInt solverId) const {
    return m_windowLayer_[domainId].at(cellId).get(solverId);
  }

  MBool isSolverWindowCell(const MInt domainId, const MInt cellId, const MInt solverId) const {
#ifndef NDEBUG
    if(m_windowLayer_[domainId].find(cellId) != m_windowLayer_[domainId].end()
       && m_windowLayer_[domainId].at(cellId).get(solverId) <= m_noSolverHaloLayers[solverId])
      TERMM_IF_COND(!m_tree.solver(cellId, solverId), "");
#endif
    return m_windowLayer_[domainId].find(cellId) != m_windowLayer_[domainId].end()
           && m_windowLayer_[domainId].at(cellId).get(solverId) <= m_noSolverHaloLayers[solverId];
  }

  MInt noSolverWindowCells(const MInt domainId, const MInt solverId) const {
    MInt cnt = 0;
    for(const auto& cellId : m_windowCells[domainId]) {
      ASSERT(m_windowLayer_[domainId].find(cellId) != m_windowLayer_[domainId].end(), "");
      if(m_windowLayer_[domainId].at(cellId).get(solverId) <= m_noSolverHaloLayers[solverId]) {
        ASSERT(m_tree.solver(cellId, solverId), "");
        ++cnt;
      }
    }
    return cnt;
  }

  MInt noSolverHaloCells(const MInt domainId, const MInt solverId, const MBool periodicCells) const {
    MInt cnt = 0;
    for(const auto& cellId : m_haloCells[domainId]) {
      if(m_tree.solver(cellId, solverId)) {
        if(periodicCells || !a_hasProperty(cellId, Cell::IsPeriodic)) ++cnt;
      }
    }
    return cnt;
  }

  /// Return number of azimuthal neighbor domains
  MInt noAzimuthalNeighborDomains() const { return (signed)m_azimuthalNghbrDomains.size(); }

  /// Return azimuthal neighbor domain
  const MInt& azimuthalNeighborDomain(const MInt id) const { return m_azimuthalNghbrDomains[id]; }

  /// Return vector of azimuthal neighbor domains
  const std::vector<MInt>& azimuthalNeighborDomains() const { return m_azimuthalNghbrDomains; }

  /// Return number of azimuthal Halo cells for given domain
  MInt noAzimuthalHaloCells(const MInt domainId) const { return m_azimuthalHaloCells[domainId].size(); }

  /// Return azimuthal Halo cell id
  MInt azimuthalHaloCell(const MInt domainId, const MInt cellId) const {
    return m_azimuthalHaloCells[domainId][cellId];
  }

  // Returns vector of azimuthal Halo cells
  const std::vector<std::vector<MInt>>& azimuthalHaloCells() const { return m_azimuthalHaloCells; };

  // Returns number of unmapped azimuthal halo cells
  MInt noAzimuthalUnmappedHaloCells() const { return m_azimuthalUnmappedHaloCells.size(); };

  // Returns unmapped azimuthal halo cell id
  MInt azimuthalUnmappedHaloCell(const MInt cellId) const { return m_azimuthalUnmappedHaloCells[cellId]; };

  // Returns unmapped azimuthal halo cell id
  MInt azimuthalUnmappedHaloDomain(const MInt cellId) const { return m_azimuthalUnmappedHaloDomains[cellId]; };

  /// Return number of azimuthal window cells for given domain
  MInt noAzimuthalWindowCells(const MInt domainId) const { return m_azimuthalWindowCells[domainId].size(); }

  /// Return azimuthal window cell id
  MInt azimuthalWindowCell(const MInt domainId, const MInt cellId) const {
    return m_azimuthalWindowCells[domainId][cellId];
  }

  // Returns vector of azimuthal window cells
  const std::vector<std::vector<MInt>>& azimuthalWindowCells() const { return m_azimuthalWindowCells; };
  /// Return grid file name
  MString gridInputFileName() const { return m_gridInputFileName; }

  /// Return maximum homogeneously-refined level
  MInt maxUniformRefinementLevel() const { return m_maxUniformRefinementLevel; }

  /// Return the reductionFactor
  MFloat reductionFactor() const { return m_reductionFactor; }

  /// Return maximum possible refinement level
  MInt maxRefinementLevel() const { return m_maxRfnmntLvl; }

  /// Return the center of gravity
  MInt centerOfGravity(MInt dir) const { return m_centerOfGravity[dir]; }

  /// Return the number of Cells on the minLevel
  MInt noMinCells() const { return m_minLevelCells.size(); }

  /// Return min-level cell id
  MInt minCell(const MInt id) const { return m_minLevelCells[id]; }

  /// Return the Global-noCells
  MLong noCellsGlobal() const { return m_noCellsGlobal; }

  /// Return noHalo-a_noHaloLayers
  MInt noHaloLayers() const { return m_noHaloLayers; }
  constexpr MInt noHaloLayers(const MInt solverId) const { return m_noSolverHaloLayers[solverId]; }
  // WH_old
  MInt haloMode() const { return m_haloMode; }

  MLong localPartitionCellOffsets(const MInt index) const { return m_localPartitionCellOffsets[index]; }
  MInt periodicCartesianDir(const MInt dir) const { return m_periodicCartesianDir[dir]; }
  MFloat periodicCartesianLength(const MInt dir) const { return m_periodicCartesianLength[dir]; }
  MFloat gridCellVolume(const MInt level) const { return m_gridCellVolume[level]; }
  MBool checkOutsideAzimuthalDomain(MFloat*);
  void tagAzimuthalHigherLevelExchangeCells(std::vector<MLong>&, /*std::vector<MLong>&,*/ std::vector<MLong>&, MInt);
  void tagAzimuthalUnmappedHaloCells(std::vector<std::vector<MLong>>&, std::vector<MLong>&, std::vector<MInt>&, MInt);
  void correctAzimuthalSolverBits();

  MInt localPartitionCellLocalIds(const MInt id) const { return m_localPartitionCellLocalIds[id]; }
  constexpr MInt noPartitionCells() const { return m_noPartitionCells; }
  constexpr MBool allowInterfaceRefinement() const { return m_allowInterfaceRefinement; }

  void gridSanityChecks();
  void checkWindowHaloConsistency(const MBool fullCheck = false, const MString a = "");
  void checkAzimuthalWindowHaloConsistency();
  void dumpCellData(const MString name);

  MBool updatedPartitionCells() const { return m_updatedPartitionCells; }
  void setUpdatedPartitionCells(const MBool flag) { m_updatedPartitionCells = flag; }

 private:
  //! coordinates of the root cell
  MFloat m_centerOfGravity[3]{};
  MFloat m_boundingBox[6]{};

  // Remove the following line once transition to new grid format is complete
  // (see also setupWindowHaloCellConnectivity())
  MFloat m_boundingBoxBackup[6]{};

  /// Local bounding box for fast checking if a point lies outside the local domain
  MFloat m_localBoundingBox[6]{};
  MBool m_localBoundingBoxInit = false;

  /// Multisolver grid information
  MFloat m_targetGridBoundingBox[6]{};
  MFloat m_targetGridCenterOfGravity[3]{};
  MFloat m_targetGridLengthLevel0 = -1;
  MInt m_targetGridMinLevel = -1;
  MBool m_hasMultiSolverBoundingBox = false;

  //! The initial length of the solver
  MFloat m_lengthLevel0{};

  // Maximum number of grid cells
  MInt m_maxNoCells = -1;

  // Window/Halo related props
  MInt m_haloMode = -1;
  MInt m_noHaloLayers = -1;
  // solver number of haloLayers
  MInt* m_noSolverHaloLayers = nullptr;

  // Holds the ids of the adjacent domains
  std::vector<MInt> m_nghbrDomains;

  // Holds the ids of the adjacent azimuthal domains
  std::vector<MInt> m_azimuthalNghbrDomains;

  // Holds the Backup for to-be-deleted periodic-neighbor-Links
  std::vector<std::tuple<MInt, MInt, MInt>> m_neighborBackup;

  // Holds the reverse-direction-Information
  const MInt m_revDir[6] = {1, 0, 3, 2, 5, 4};

  std::vector<std::vector<MInt>> m_windowCells;

  // Holds the ids of the halo cells according to original adjacent domain
  std::vector<std::vector<MInt>> m_haloCells;

  /* m_windowLayer_[nghbrDomId][cellId][solverId] stores the window layer of this cell w.r.t.
   * solverId. The entries are initialized by the higheset possible value of 15 (corresponding
   * to all bits set to 1), which indicates that cell is no window, so effectively allowed range is
   * noHaloLayers<15. Note that this is not exactly the windowLayer, but an auxiliary value
   * required by the algorithm.
   */
  std::vector<std::unordered_map<MInt, M32X4bit<true>>> m_windowLayer_;
  static constexpr const auto WINDOWLAYER_MAX = decltype(m_windowLayer_)::value_type::mapped_type::MAX;

  // Azimuthal periodicity exchange data
  std::vector<std::vector<MInt>> m_azimuthalHigherLevelConnectivity;
  std::vector<std::vector<MInt>> m_azimuthalWindowCells;
  std::vector<std::vector<MInt>> m_azimuthalHaloCells;
  std::vector<MInt> m_azimuthalUnmappedHaloCells;
  std::vector<MInt> m_azimuthalUnmappedHaloDomains;
  std::vector<MInt> m_gridBndryCells;

  // maximum size of a partitionCell (important for efficient distribution of cells)
  MInt m_partitionCellOffspringThreshold;
  MFloat m_partitionCellWorkloadThreshold;

  // offset array
  MLong* m_domainOffsets = nullptr;

  // save the offsets of the local partitionCells and the globalIds of the local partitionCells for later use (e.g. the
  // particle restarting) chrs 26.11
  MInt m_noPartitionCells{};
  MLong m_noPartitionCellsGlobal{};
  MInt m_noHaloPartitionLevelAncestors{};
  MLong m_localPartitionCellOffsets[3]{static_cast<MLong>(-1), static_cast<MLong>(-1), static_cast<MLong>(-1)};
  MLong* m_localPartitionCellGlobalIds = nullptr;
  MInt* m_localPartitionCellLocalIds = nullptr;
  MBool m_updatePartitionCellsOnRestart = true;

  // Indicates that the partition cells are currently being updated (e.g. to handle the case when
  // there are 0 local partition cells temporarily before balancing)
  MBool m_updatingPartitionCells = false;
  // Indicate if the partition cells have been updated/changed during balancing
  MBool m_updatedPartitionCells = false;

  /// Local number of partition level ancestors on this domain (without halos)
  MInt m_noPartitionLevelAncestors{};

  MLong m_noPartitionLevelAncestorsGlobal{};
  MLong m_32BitOffset = 0;
  MBool m_restart;

  MInt m_noMinLevelCellsGlobal{};

  std::vector<MInt> m_minLevelCells{};

  MInt m_maxPartitionLevelShift{};

  MBool* m_checkRefinementHoles = nullptr;
  MBool* m_diagSmoothing = nullptr;
  MInt m_noIdenticalSolvers = 0;
  MInt* m_identicalSolvers = nullptr;
  MInt* m_identicalSolverMaxLvl = nullptr;
  MInt* m_identicalSolverLvlJumps = nullptr;

 private:
  /// Partition level shift

  /// List of partition level ancestor cell ids (including halos of the missing subtree of the grid
  /// in case of a partition level shift)
  std::vector<MInt> m_partitionLevelAncestorIds{};

  std::vector<MLong> m_partitionLevelAncestorChildIds{};
  std::vector<MInt> m_partitionLevelAncestorNghbrDomains{};
  std::vector<std::vector<MInt>> m_partitionLevelAncestorHaloCells{};
  std::vector<std::vector<MInt>> m_partitionLevelAncestorWindowCells{};


  MInt** m_neighborList = nullptr;

  MFloat m_reductionFactor{};
  MInt m_decisiveDirection{};

  MString m_outputDir = "";
  MString m_restartDir = "";

  MBool m_loadGridPartition = false;

  MBool m_loadPartition = false;

  MBool m_partitionParallelSplit = false;

  MBool m_addSolverToGrid = false;
  MInt m_referenceSolver = -1;

 public:
  // Data structures for neighbor find and walk algorithms
  // Number of possible paths for 1D neighbors
  MInt m_counter1D{};
  // Number of possible paths for 2D neighbors
  MInt m_counter2D{};
  // Number of possible paths for 3D neighbors
  MInt m_counter3D{};

  // Holds all possible paths to 1D neighbors
  MInt* m_paths1D = nullptr;
  // Holds all possible paths to 2D neighbors
  MInt** m_paths2D = nullptr;
  // Holds all possible paths to 3D neighbors
  MInt** m_paths3D = nullptr;

  // 2D: 11 elements, 3D: 43 elements
  std::array<MInt, 11 + (nDim - 2) * 32> m_neighborCode{};

  /// Mapping: global cell id -> local cell id
  std::map<MLong, MInt> m_globalToLocalId{};

  // Status flag for additional LB grid checks
  MBool m_lbGridChecks;
  MFloat m_coarseRatio = 0.2;
  MBool m_allowCoarsening = true;

  std::set<MInt> m_freeIndices;

  // Mapping from domainId to neighborDomainId
  std::vector<MInt> m_nghbrDomainIndex;

  // Mapping from domainId to azimuthalNeighborDomainId
  std::vector<MInt> m_azimuthalNghbrDomainIndex;
  MInt* m_localPartitionCellLocalIdsRestart = nullptr;
  MLong m_localPartitionCellOffsetsRestart[3]{static_cast<MLong>(-1), static_cast<MLong>(-1), static_cast<MLong>(-1)};

 private:
  // Store the MPI communicator to use for all MPI communication
  const MPI_Comm m_mpiComm;
  /// Store the tree
  Tree m_tree;
  const MBool m_paraViewPlugin;
  const MInt m_maxNoNghbrs;
  static constexpr MInt m_noDirs = 2 * nDim;
  static constexpr MInt m_maxNoChilds = IPOW2(nDim);
  // Jannik TODO labels:GRID for now needed for testcase hack -> change tagActiveWindows
  MBool m_zonal;

 public:
  CartesianGrid(MInt maxCells, const MFloat* const bBox, const MPI_Comm comm, const MString& fileName = "");
  ~CartesianGrid();

  /// Full access to tree (for now)
  // Accessor to tree is named "treeb" to indicate that this is the tree of the base grid and not a
  // solver-specific grid tree
  Tree& treeb() { return m_tree; }
  const Tree& treeb() const { return m_tree; }

  static constexpr MInt m_maxNoSensors = 64;

  MInt maxLevel() const { return m_maxLevel; };
  MInt minLevel() const { return m_minLevel; };

  MInt targetGridMinLevel() const { return m_targetGridMinLevel; };
  MInt maxPartitionLevelShift() const { return m_maxPartitionLevelShift; };


  // Return maximum number of cells
  MInt maxNoCells() const { return m_maxNoCells; };

  /// Return the MPI communicator used by this grid
  MPI_Comm mpiComm() const { return m_mpiComm; }

  /// Return the 32-BitOffset
  MLong bitOffset() const { return m_32BitOffset; }

  /// Return the cell Id of the partition cell containing the coords
  template <MBool t_correct = false>
  MInt findContainingPartitionCell(const MFloat* const coord, const MInt solverId = -1,
                                   std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr);

  // Return the cell Id of the partition cell intersecting with a sphere
  MInt intersectingWithPartitioncells(MFloat* const center, MFloat const radius);

  MInt intersectingWithHaloCells(MFloat* const center, MFloat const radius, MInt domainId, MBool onlyPartition);

  MBool intersectingWithLocalBoundingBox(MFloat* const center, MFloat const radius);

  MBool boxSphereIntersection(const MFloat* bMin, const MFloat* bMax, const MFloat* const sphereCenter,
                              MFloat const radius);

  MBool hollowBoxSphereIntersection(const MFloat* bMin, const MFloat* bMax, const MFloat* const sphereCenter,
                                    MFloat const radius);

  /// Return the cell Id of the partition cell that has the given cellId as a descendant
  inline MInt ancestorPartitionCellId(const MInt cellId) const;

  /// Return the cell id of the halo partition cell containing the coords
  inline MInt findContainingHaloPartitionCell(const MFloat* const coors, const MInt solverId = -1);

  /// Return the cell id of the halo partition cell containing the coords
  template <MBool t_correct = false>
  MInt findContainingHaloCell(const MFloat* const coord, const MInt solverId, MInt domainId, MBool onlyPartitionCells,
                              std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr);

  /// Return the cell Id of the leaf cell containing the coords
  template <MBool t_correct = false>
  MInt findContainingLeafCell(const MFloat* coord,
                              std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr,
                              const MInt solverId = -1);

  template <MBool t_correct = false>
  MInt findContainingLeafCell(const std::array<MFloat, nDim>& coord, const MInt startId,
                              std::function<MFloat*(MInt, MFloat* const)>* correctCellCoord = nullptr,
                              const MInt solverId = -1, const MBool allowNonLeafHalo = false) {
    return findContainingLeafCell<t_correct>(&coord[0], startId, correctCellCoord, solverId, allowNonLeafHalo);
  }
  template <MBool t_correct = false>
  MInt findContainingLeafCell(const MFloat* const coord, const MInt startId,
                              std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr,
                              const MInt solverId = -1, const MBool allowNonLeafHalo = false);

  template <MBool t_correct = false, MBool insideLimit>
  inline MBool pointWthCell(const std::array<MFloat, nDim>& coord, const MInt cellId,
                            std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr) const {
    return pointWthCell<t_correct, insideLimit>(&coord[0], cellId, correctCellCoord);
  }

  template <MBool t_correct = false, MBool insideLimit>
  inline MBool pointWthCell(const MFloat* const coord, const MInt cellId,
                            std::function<MFloat*(MInt, MFloat* const)> correctCellCoord = nullptr) const;

  inline MFloat* dummyCorrect(const MInt cellId, MFloat* coords) {
    std::copy(&a_coordinate(cellId, 0), &a_coordinate(cellId, nDim - 1), coords);
    return coords;
  }

  /// Compute the bounding box of all local cells
  void computeLocalBoundingBox(MFloat* const bbox);

  /// Check if the given point lies in the local bounding box
  MBool pointInLocalBoundingBox(const MFloat* coord);

  void propagateDistance(std::vector<MInt>& list, MIntScratchSpace& distMem, MInt dist);
  void propagationStep(MInt cellId, MInt dist, MInt* distMem, MInt endDist);

  template <typename DATATYPE>
  void exchangeNotInPlace(DATATYPE* exchangeVar, DATATYPE* recvMem);
  template <typename DATATYPE>
  void generalExchange(MInt noVars, const MInt* vars, DATATYPE** exchangeVar, MInt noDomSend, MInt* domSend,
                       const MInt* noCellsSendPerDom, const MInt* cellIdsSend, MInt noDomRecv, MInt* domRecv,
                       const MInt* noCellsRecvPerDom, const MInt* cellIdsRecv);

  // Create a 2D slice gird from 3D grid
  void createGridSlice(const MString& axis, const MFloat intercept, const MString& fileName);
  void createGridSlice(const MString& axis, const MFloat intercept, const MString& fileName, const MInt solverId,
                       MInt* const noSliceCellIds, MInt* const sliceCellIds, MInt* const noSliceHilbertIds = nullptr,
                       MInt* const sliceHilbertInfo = nullptr, MInt* const noSliceContHilbertIds = nullptr,
                       MInt* const sliceContiguousHilbertInfo = nullptr);

  void computeGlobalIds();
  void localToGlobalIds();
  void descendNoOffsprings(const MLong cellId, const MLong offset = 0);
  void storeMinLevelCells(const MBool updateMinlevel = false);
  void computeLeafLevel();
  void determineNoPartitionCellsAndOffsets(MLong* const noPartitionCells, MLong* const partitionCellOffsets);

  void determinePartLvlAncestorHaloWindowCells(std::vector<std::vector<MInt>>& partLvlAncestorHaloCells,
                                               std::vector<std::vector<MInt>>& partLvlAncestorWindowCells);

  void balance(const MInt* const noCellsToReceiveByDomain, const MInt* const noCellsToSendByDomain,
               const MInt* const sortedCellId, const MLong* const offset, const MLong* const globalIdOffsets);

  MBool updatePartitionCells(MInt offspringThreshold, MFloat workloadThreshold);
  MLong generateHilbertIndex(const MInt cellId, const MInt targetMinLevel) {
    return hilbertIndexGeneric(&a_coordinate(cellId, 0), m_targetGridCenterOfGravity, m_targetGridLengthLevel0,
                               targetMinLevel);
  }

 private:
  void createGlobalToLocalIdMapping();
  void changeGlobalToLocalIds();
  void descendStoreGlobalId(MInt cellId, MInt& localCnt);
  void refineCell(const MInt cellId, const MLong* const refineChildIds = nullptr, const MBool mayHaveChildren = false,
                  const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& =
                      std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>(),
                  const std::bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()> refineFlag =
                      std::bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()>(1ul));
  void removeChilds(const MInt cellId);
  template <MBool removeChilds>
  void removeCell(const MInt cellId);
  MInt deleteCell(const MInt cellId);
  void createPaths();
  void setLevel();
  void updatePartitionCellInformation();

  void setGridInputFilename(MBool = false);
  void loadGridFile(const MInt, const MInt);

  MLong hilbertIndexGeneric(const MFloat* const coords, const MFloat* const center, const MFloat length0,
                            const MInt baseLevel) const;

  inline MLong hilbertIndexGeneric(const MFloat* const coords) const {
    return hilbertIndexGeneric(coords, m_targetGridCenterOfGravity, m_targetGridLengthLevel0, m_targetGridMinLevel);
  }

  MLong generateHilbertIndex(const MInt cellId) { return hilbertIndexGeneric(&a_coordinate(cellId, 0)); }

  template <class ITERATOR, typename U>
  MInt findIndex(ITERATOR, ITERATOR, const U&);
  template <typename TA, typename TB, typename TC>
  MInt setNeighborDomainIndex(const MInt, std::vector<TA>&, std::vector<TB>&, std::vector<TC>&);
  template <typename TA, typename TB>
  MInt setNeighborDomainIndex(const MInt, std::vector<TA>&, std::vector<TB>&);

  MInt setAzimuthalNeighborDomainIndex(const MInt, std::vector<std::vector<MInt>>&, std::vector<std::vector<MInt>>&);
  MInt setAzimuthalNeighborDomainIndex(const MInt, std::vector<std::vector<MInt>>&, std::vector<std::vector<MLong>>&);

  template <std::size_t N>
  void exchangeSolverBitset(std::bitset<N>* const data, const MBool defaultVal = false);
  void checkWindowLayer(const MString a);

 public:
  void saveGrid(const MChar*, const std::vector<std::vector<MInt>>&, const std::vector<std::vector<MInt>>&,
                const std::vector<std::vector<MInt>>&, const std::vector<MInt>&, const std::vector<std::vector<MInt>>&,
                MInt* recalcIdTree = nullptr);
  void saveGrid(const MChar*, MInt* recalcIdTree);
  void savePartitionFile();
  void loadPartitionFile(const MString& partitionFileName, MLong* partitionCellOffsets);
  void savePartitionCellWorkloadsGridFile();

  void deletePeriodicConnection(const MBool = true);
  void restorePeriodicConnection();
  MInt findNeighborDomainId(const MLong globalId) {
    return findNeighborDomainId(globalId, noDomains(), m_domainOffsets);
  };
  MInt findNeighborDomainId(const MLong globalId, const MInt noDomains, const MLong* domainOffsets);

  void createLeafCellMapping(const MInt donorSolverId, const MInt gridCellId, std::vector<MInt>& mappedLeafCells,
                             MBool allChilds = false);
  MInt globalIdToLocalId(const MLong& globalId, const MBool termIfNotExisting = false);

 private:
  void loadGrid(const MString& fileName);
  void resetCell(const MInt& cellId);
  void initGridMap();
  void createGridMap(const MString& donorGridFileName, const MString& gridMapFileName);
  void saveDonorGridPartition(const MString& gridMapFileName, const MString& gridPartitionFileName);
  void loadDonorGridPartition(const MLong* const partitionCellsId, const MInt noPartitionCells);

  void savePartitionFile(const MString& partitionFileNameBase, const MLong partitionCellOffset);

  template <typename CELLTYPE>
  void calculateNoOffspringsAndWorkload(Collector<CELLTYPE>* input_cells, MInt input_noCells);
  void partitionParallel(const MInt tmpCount, const MLong tmpOffset, const MFloat* const partitionCellsWorkload,
                         const MLong* const partitionCellsGlobalId, const MFloat totalWorkload,
                         MLong* partitionCellOffsets, MLong* globalIdOffsets, MBool computeOnlyPartition = false);

  void updateHaloCellCollectors();
  void tagActiveWindows(std::vector<MLong>&, MInt);
  void tagActiveWindows2_(std::vector<MLong>&, const MInt);
  void tagActiveWindowsOnLeafLvl3(
      const MInt maxLevel, const MBool duringMeshAdaptation,
      const std::vector<std::function<void(const MInt)>>& = std::vector<std::function<void(const MInt)>>());
  void tagActiveWindowsAtMinLevel(const MInt);
  MInt getAdjacentGridCells(MInt, MIntScratchSpace&, MInt, MBool diagonalNeighbors = true);
  template <MBool finer, MBool coarser>
  inline MInt getAdjacentGridCells5(const MInt, MInt* const, const MInt level = -1,
                                    const MBool diagonalNeighbors = true);
  template <MBool finer, MBool coarser>
  inline void getAdjacentGridCells1d5(std::set<MInt>&, const MInt, const MInt, const MInt dimNext = -1,
                                      const MInt dimNextNext = -1, const uint_fast8_t childCodes = (uint_fast8_t)~0);
  void setupWindowHaloCellConnectivity();
  void createMinLevelExchangeCells();
  void createHigherLevelExchangeCells(
      const MInt onlyLevel = -1, const MBool duringMeshAdaptation = false,
      const std::vector<std::function<void(const MInt)>>& = std::vector<std::function<void(const MInt)>>(),
      const std::vector<std::function<void(const MInt)>>& = std::vector<std::function<void(const MInt)>>(),
      const MBool forceLeafLvlCorrection = false);
  void createHigherLevelExchangeCells_old(
      const MInt onlyLevel = -1,
      const std::vector<std::function<void(const MInt)>>& = std::vector<std::function<void(const MInt)>>());
  MInt createAdjacentHaloCell(const MInt, const MInt, const MLong*, std::unordered_multimap<MLong, MInt>&,
                              const MFloat*, MInt* const, std::vector<std::vector<MInt>>&,
                              std::vector<std::vector<MLong>>&);
  MInt createAzimuthalHaloCell(const MInt, const MInt, const MLong*, std::unordered_multimap<MLong, MInt>&,
                               const MFloat*, MInt* const, std::vector<std::vector<MInt>>&,
                               std::vector<std::vector<MLong>>&);
  void meshAdaptation(std::vector<std::vector<MFloat>>&, std::vector<MFloat>&, std::vector<std::bitset<64>>&,
                      std::vector<MInt>&, const std::vector<std::function<void(const MInt)>>&,
                      const std::vector<std::function<void(const MInt)>>&,
                      const std::vector<std::function<void(const MInt)>>&,
                      const std::vector<std::function<void(const MInt, const MInt)>>&,
                      const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>&,
                      const std::vector<std::function<void()>>&);
  void meshAdaptationLowMem(std::vector<std::vector<MFloat>>&, std::vector<MFloat>&, std::vector<std::bitset<64>>&,
                            std::vector<MInt>&, const std::vector<std::function<void(const MInt)>>&,
                            const std::vector<std::function<void(const MInt)>>&,
                            const std::vector<std::function<void(const MInt)>>&,
                            const std::vector<std::function<void(const MInt, const MInt)>>&,
                            const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>&,
                            const std::vector<std::function<void()>>&);

  void meshAdaptationDefault(std::vector<std::vector<MFloat>>&, std::vector<MFloat>&, std::vector<std::bitset<64>>&,
                             std::vector<MInt>&, const std::vector<std::function<void(const MInt)>>&,
                             const std::vector<std::function<void(const MInt)>>&,
                             const std::vector<std::function<void(const MInt)>>&,
                             const std::vector<std::function<void(const MInt, const MInt)>>&,
                             const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>&,
                             const std::vector<std::function<void()>>&);


  void compactCells(const std::vector<std::function<void(const MInt, const MInt)>>& =
                        std::vector<std::function<void(const MInt, const MInt)>>());
  void swapCells(const MInt cellId, const MInt otherId);
  MBool refineCheck(const MInt cellId, const MInt solverId,
                    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver);
  MBool coarsenCheck(const MInt cellId, const MInt solverId);
  void exchangeProperties();
  void setChildSolverFlag(const MInt cellId, const MInt solver,
                          const std::function<MInt(const MFloat*, const MInt, const MInt)>& cellOutsideSolver);

 public:
  /// Return the total number of domains (total number of ranks in current MPI
  /// communicator)
  MInt noDomains() const { return m_noDomains; }

  /// Return the domainId (rank)
  MInt domainId() const { return m_domainId; }

  void cartesianToCylindric(const MFloat*, MFloat*);
  void rotateCartesianCoordinates(MFloat*, MFloat);

 private:
  // holds the domain id
  MInt m_domainId{};

  /// Return true if this is domain zero (a.k.a. as the root domain)
  MBool isMpiRoot() const { return domainId() == 0; }
};


/** \brief Reduces the current grid to a certain level
 *
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * The current grid is reduced by calling interpolateToParentCells.
 *
 * Requires rewriting of certain arrays for the execution of the other subroutines:
 *
 * 1. Save the globalIds of the window/halocells for later identification
 * 2. Reallocate m_windowCells and m_haloCells
 *    (saveGridDomain requires the arrays to be of size m_noDomains, currenty has size noNeighborDomains())
 * 3. Delete all cells above a certain level provided by m_reductionLevel, uses deleteCell from CartesianGrid
 * 4. Update the m_windowCell and m_haloCell arrays, therefore also reallocate these arrays
 *    (required by saveGridDomainPar)
 * 5. Update m_domainOffsets since the sizes of the Collectors for each process has changed,
 *    this requires global communication (this is required by the IO routine of the solver)
 */
template <MInt nDim>
template <class F>
MInt CartesianGrid<nDim>::reduceToLevel(const MInt reductionLevel, F interpolateToParentCells) {
  TRACE();
  // 1. before we delete, generate a copy of the global ids for the domains so
  // that we can identify them later (as the global id will not be updated
  // during delete cell)
  std::vector<std::vector<MInt>> old_globalhaloCells, old_globalwindowCells;

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    std::vector<MInt> tmp_halo;
    std::vector<MInt> tmp_window;
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      tmp_halo.push_back(a_globalId(m_haloCells[i][j]));
    }
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      tmp_window.push_back(a_globalId(m_windowCells[i][j]));
    }
    old_globalhaloCells.push_back(tmp_halo);
    old_globalwindowCells.push_back(tmp_window);
  }

  // 2. we need to realloc m_haloCells and m_windowCells so that we
  // suit the requirements of saveGridDomainPar
  m_haloCells.clear();
  m_windowCells.clear();
  m_haloCells.resize(noDomains());
  m_windowCells.resize(noDomains());

  // 3. Delete all cells that are above a certain level, interpolate all
  // finer cells to coarser cells
  for(MInt level = m_maxLevel; level > reductionLevel; level--) {
    interpolateToParentCells(level - 1);
    for(MInt c = m_tree.size() - 1; c >= 0; c--) {
      if(a_level(c) == level) {
        deleteCell(c);
      }
    }
  }
  m_noInternalCells = m_tree.size();

  // 4. Update window and halo cell arrays

  for(MInt c = 0; c < m_tree.size(); c++) {
    if(a_hasProperty(c, Cell::IsHalo)) {
      MInt nghr = -1;
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(std::size_t j = 0; j < old_globalhaloCells[i].size(); j++)
          if(a_globalId(c) == old_globalhaloCells[i][j]) {
            nghr = i;
            break;
          }
        if(nghr != -1) {
          break;
        }
      }
      m_haloCells[m_nghbrDomains[nghr]].push_back(c);

    } else if(a_hasProperty(c, Cell::IsWindow)) {
      MInt nghr = -1;
      std::vector<MInt> t;
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(std::size_t j = 0; j < old_globalwindowCells[i].size(); j++)
          if(a_globalId(c) == old_globalwindowCells[i][j]) {
            nghr = i;
            t.push_back(i);
            break;
          }
        if(nghr != -1) {
          continue;
        }
      }
      for(MInt& i : t) {
        m_windowCells[m_nghbrDomains[i]].push_back(c);
      }
    }
  }


  // Update the number of internal cells
  for(MInt i = 0; i < noDomains(); i++) {
    m_noInternalCells -= (signed)m_haloCells[i].size();
  }

  // 5. Update m_domainOffsets
  MIntScratchSpace rcvBufRed(noDomains(), AT_, "rcvBufRed");

  // fill buffer and exchange
  MPI_Allgather(&(m_noInternalCells), 1, MPI_INT, rcvBufRed.begin(), 1, MPI_INT, mpiComm(), AT_, "(m_noInternalCells)",
                "rcvBufRed.begin()");

  // put the result into the array
  MInt offset = 0;
  for(MInt i = 0; i < noDomains(); i++) {
    m_domainOffsets[i] = offset;
    offset += rcvBufRed[i];
  }
  // the last element holds the total number
  m_domainOffsets[noDomains()] = offset;

  return rcvBufRed[domainId()];
}

#endif
