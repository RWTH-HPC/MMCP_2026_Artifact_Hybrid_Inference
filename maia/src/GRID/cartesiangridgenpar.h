// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDGENPAR_H
#define GRIDGENPAR_H

#include <bitset>
#include <set>
#include <stack>

#include <tuple>
#include "GEOM/geometry.h"
#include "GEOM/geometry2d.h"
#include "GEOM/geometry3d.h"
#include "GEOM/geometrycontext.h"
#include "GEOM/geometryelement.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/collector.h"
#include "cartesiangridgencell.h"
#include "cartesiannetcdf.h"
#include "globals.h"

class GeometryRoot;

struct SolverRefinement {
  // members
  MInt maxUniformRefinementLevel;
  MInt maxRfnmntLvl;
  MInt maxBoundaryRfnLvl;
  // refinement patch
  MInt noLocalPatchRfnLvls() { return (MInt)localRfnLevelMethods.size(); };
  MInt noPatchesPerLevel(MInt addLevel) { return (MInt)localRfnLevelMethods[addLevel].size(); };
  std::vector<MString> localRfnLevelMethods;
  MInt* localRfnLevelPropertiesOffset = nullptr;
  MInt* noLocalRfnPatchProperties = nullptr;
  MFloat** localRfnPatchProperties = nullptr;
  // boundary refinement
  MInt noLocalBndRfnLvls;
  MInt localBndRfnMethod;
  MInt noLocalRfnBoundaryIds;
  MInt* localRfnBoundaryIds = nullptr;
  MInt* localMinBoundaryThreshold = nullptr;
  MInt* smoothDistance = nullptr;
  MInt* localRfnLvlDiff = nullptr;
  MInt* localBndRfnMinLvlDiff = nullptr;
  MFloat* localBndRfnDistance = nullptr;
  // cut off
  MInt cutOff;
  std::vector<MString> cutOffMethods;
  MFloat** cutOffCoordinates = nullptr;
  MInt** cutOffNmbrLayers = nullptr;
};

template <MInt nDim>
class GridgenPar {
 public:
  /// \brief Returns the parent of the cell \p cellId
  MLong& a_parentId(const MInt cellId) { return *m_pCells[cellId].m_parentId_; }
  /// \brief Returns the parent of the cell \p cellId
  const MLong& a_parentId(const MInt cellId) const { return *m_pCells[cellId].m_parentId_; };
  /// \brief Returns the level of the cell \p cellId
  MInt& a_level(const MInt cellId) { return *m_pCells /**/[cellId].m_level_; }
  /// \brief Returns the level of the cell \p cellId
  const MInt& a_level(const MInt cellId) const { return *m_pCells /**/[cellId].m_level_; };
  /// \brief Returns the solid layer of the cell \p cellId
  MInt& a_noSolidLayer(const MInt cellId, const MInt solver) { return m_pCells /**/[cellId].m_noSolidLayer_[solver]; };
  /// \brief Returns the solid layer of the cell \p cellId
  const MInt& a_noSolidLayer(const MInt cellId, const MInt solver) const {
    return m_pCells /**/[cellId].m_noSolidLayer_[solver];
  };
  /// \brief Returns the globalId of the cell \p cellId
  MLong& a_globalId(const MInt cellId) { return *m_pCells[cellId].m_globalId_; }
  /// \brief Returns the globalId of the cell \p cellId
  const MLong& a_globalId(const MInt cellId) const { return *m_pCells[cellId].m_globalId_; };
  /// \brief Returns the no. of children of the cell \p cellId
  MInt& a_noChildren(const MInt cellId) { return *m_pCells[cellId].m_noChildIds_; }
  /// \brief Returns the no. of children of the cell \p cellId
  const MInt& a_noChildren(const MInt cellId) const { return *m_pCells[cellId].m_noChildIds_; };
  /// \brief Returns the coordinate of the cell \p cellId for dimension \p dim
  MFloat& a_coordinate(const MInt cellId, const MInt dim) { return m_pCells /**/[cellId].m_coordinates_[dim]; }
  /// \brief Returns the coordinate of the cell \p cellId for dimension \p dim
  const MFloat& a_coordinate(const MInt cellId, const MInt dim) const {
    return m_pCells /**/[cellId].m_coordinates_[dim];
  };
  /// \brief Returns the child id of the cell \p cellId \p position
  MLong& a_childId(const MInt cellId, const MInt position) { return m_pCells /**/[cellId].m_childIds_[position]; }
  /// \brief Returns the child id of the cell \p cellId \p position
  const MLong& a_childId(const MInt cellId, const MInt position) const {
    return m_pCells /**/[cellId].m_childIds_[position];
  }
  /// \brief Returns property \p p of the cell \p cellId
  MChar& a_hasProperty(const MInt cellId, const MInt p) { return m_pCells /**/[cellId].b_properties_[p]; }
  /// \brief Returns property \p p of the cell \p cellId
  MChar a_hasProperty(const MInt cellId, const MInt p) const { return m_pCells /**/[cellId].b_properties_[p]; };
  /// \brief Returns rfnDistance of the cell \p cellId \p position
  MInt& a_refinementDistance(const MInt cellId) { return m_pCells /**/[cellId].m_rfnDistance_; }

  /// \brief Returns rfnDistance of the cell \p cellId \p position
  const MInt& a_refinementDistance(const MInt cellId) const { return m_pCells /**/[cellId].m_rfnDistance_; }

  /// \brief Returns the neighbor id of the cell \p cellId \p position
  MLong& a_neighborId(const MInt cellId, const MInt position) { return m_pCells /**/[cellId].m_nghbrIds_[position]; }

  /// \brief Returns the neighbor id of the cell \p cellId \p position
  const MLong& a_neighborId(const MInt cellId, const MInt position) const {
    return m_pCells /**/[cellId].m_nghbrIds_[position];
  }

  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar& a_isInSolver(const MInt cellId, const MInt solver) {
    return m_pCells /**/[cellId].b_solverAffiliation_[solver];
  }
  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar a_isInSolver(const MInt cellId, const MInt solver) const {
    return m_pCells /**/[cellId].b_solverAffiliation_[solver];
  };

  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar& a_isSolverBoundary(const MInt cellId, const MInt solver) {
    return m_pCells /**/[cellId].b_solverBoundary_[solver];
  }
  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar a_isSolverBoundary(const MInt cellId, const MInt solver) const {
    return m_pCells /**/[cellId].b_solverBoundary_[solver];
  };

  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar& a_isToRefineForSolver(const MInt cellId, const MInt solver) {
    return m_pCells /**/[cellId].b_solverToRefine_[solver];
  }
  /// \brief Multisolver grid: does a cell belong to a certain solver
  MChar a_isToRefineForSolver(const MInt cellId, const MInt solver) const {
    return m_pCells /**/[cellId].b_solverToRefine_[solver];
  };

  GridgenPar(const MPI_Comm comm, const MInt noSolvers);
  template <class T>
  static void quickSort(T* globalIdArray, MInt* lookup, MInt startindex, MInt endindex);


  MPI_Comm mpiComm() const { return m_mpiComm; }
  MInt domainId() const { return m_domainId; }
  MInt noDomains() const { return m_noDomains; }

 protected:
  void initTimers();
  void initMembers();
  void readProperties();
  void readSolverProperties(MInt solver);
  void initGeometry();
  void gridAlignCutOff();
  void createInitialGrid();
  void createStartGrid();
  void finalizeGrid();
  void createComputationalMultisolverGrid();
  void checkLBRefinementValidity();
  MInt getAdjacentGridCells(MInt cellId, MInt* adjacentCells);
  void refineComputationalGrid(MInt lvl);
  void saveGrid();
  void saveGridDomain(MInt level_, MInt tag);

  void checkMemoryAvailability(MInt stage, MInt level_);
  void refineGrid(MInt** offsets, MInt level_, MBool halo);
  void refineGridPatch(MInt** offsets, MInt level_, MBool halo);
  void refineCell(MInt id, MInt* currentChildId);
  void deleteOutsideCellsSerial(MInt level_);
  void deleteCoarseSolidCellsSerial();
  void deleteCoarseSolidCellsParallel();
  void performCutOff(MInt** offsets, MInt level_, MBool deleteMode = false);
  void keepOutsideBndryCellChildrenSerial(MInt* offsets, MInt level_);
  MInt nghborStencil(MInt, MInt);
  void createSolidCellLayer(MInt cellId, MInt solidLayer, MInt finalLayer_);
  void keepOutsideBndryCellChildrenParallel(MInt level_);
  void deleteOutsideCellsParallel(MInt level_);
  MInt getLastNonWindowCell(MInt no_consider, MInt last);
  void deleteCellReferences(MInt cell, MInt pos);
  void markInsideOutside(MInt** offsets, MInt level_);
  void markInsideOutside(MInt** offsets, MInt level_, MInt solver);
  void markSolverAffiliation(MInt level_);
  void excludeInsideOutside(MInt** offsets, MInt level_, MInt solver);
  void floodCells(std::stack<MInt>* fillStack, MChar marker);
  void floodCells(std::stack<MInt>* fillStack, MChar marker, MInt solver);
  void findChildLevelNeighbors(MInt** offsets, MInt level_);

  void reorderCellsHilbert();

  void parallelizeGrid();
  void updateHaloOffsets(MInt l, MInt noHalos, MInt* rfnCountHalosDom);
  void updateOffsets(MInt gridLevel);
  void updateInterRankNeighbors();
  void findHaloAndWindowCells(std::vector<std::vector<MInt>>& winCellIdsPerDomain,
                              std::vector<std::vector<MInt>>& haloCellIdsPerDomain);
  void findHaloAndWindowCellsKD(std::vector<std::vector<MInt>>& winCellIdsPerDomain,
                                std::vector<std::vector<MInt>>& haloCellIdsPerDomain);
  void updateGlobalIdsReferences();
  void collectHaloChildren(MInt parentId, std::vector<MInt>* cellIdsPerDomain);
  void collectWindowChildren(MInt parentId, std::vector<MInt>* cellIdsPerDomain);
  void reorderGlobalIdsDF();
  void traverseDFGlobalId(MInt parentId,
                          MLong* globalId_,
                          MLongScratchSpace& partitionCellList,
                          MFloatScratchSpace& workloadPerCell,
                          MFloat* currentWorkload,
                          MFloatScratchSpace& weight,
                          MFloatScratchSpace& workload,
                          MInt* j);
  void setCellWeights(MFloatScratchSpace& weight);
  void determineRankOffsets(MIntScratchSpace& offsets);

  void markSolverForRefinement(MInt level_, MInt solver);
  void concludeSolverRefinement(MInt level_);
  void markLocalSolverRefinement(MInt level_, MInt solver);
  void markPatchForSolverRefinement(MInt level_, MInt solver);
  void markBndForSolverRefinement(MInt level_, MInt solver);
  void markLocalBox(MInt level_, MInt patch, MInt solver);
  void markLocalRadius(MInt level_, MInt patch, MInt solver);
  void markLocalCylinder(MInt level_, MInt patch, MInt solver, MString patchStr);
  void markLocalCone(MInt level_, MInt patch, MInt solver);
  void markLocalRectangleAngled(MInt level_, MInt patch, MInt solver);
  void markLocalCartesianWedge(MInt level_, MInt patch, MInt solver);
  void markLocalFlatCone(MInt level_, MInt patch, MInt solver);
  void markLocalSlicedCone(MInt level_, MInt patch, MInt solver, MString patchStr);
  void markLocalHat(MInt level_, MInt patch, MInt solver);
  void markBndDistance(MInt level_, MInt solver);
  void propagateDistance(MInt level_, MInt distance, std::vector<MInt>& rfnBoundaryGroup, MInt solver);
  void propagationStep(MInt cellId, MInt rfnDistance, MInt finalDistance, MInt solver);
  MBool isInsideSlicedCone(const MFloat* const pointCoord,
                           const MFloat* const leftCoord,
                           const MFloat* const rightCoord,
                           const MFloat* const leftNormal,
                           const MFloat* const rightNormal,
                           const MFloat* const normalDifforigAB,
                           const MFloat leftR,
                           const MFloat rightR);
  MBool isInsideCylinder(const MFloat* const pointCoord,
                         const MFloat* const leftCoord,
                         const MFloat* const rightCoord,
                         const MFloat* const normalDiffAB,
                         const MFloat* const normalDiffBA,
                         const MFloat radius,
                         const MFloat innerRadius);
  MBool pointIsInside(MFloat* coordinates);
  MBool pointIsInsideSolver(MFloat* coordinates, MInt solver);
  MBool checkCellForCut(MInt id);

  void copyCell(MInt from, MInt to);
  void swapCells(MInt cellId1, MInt cellId2);

  void checkNeighborhood(MInt level_);
  void checkNeighborhoodDistance(MInt level_);
  void checkNeighborhoodIntegrity(MInt level_);
  void writeGridInformationPar();
  void writeGridInformation(MInt level_, MInt tag);

  // functions for parallel geometry
  void writeParallelGeometry();

  // functions used in the dynamic load balancing.
  void checkLoadBalance(MInt in_level);
  void dynamicLoadBalancing();
  void communicateInt(MIntScratchSpace& recvBuffer,
                      MIntScratchSpace& noCellsToReceive,
                      MIntScratchSpace& sendBuffer,
                      MIntScratchSpace& noCellsToSend);
  void communicateLong(MLongScratchSpace& recvBuffer,
                       MIntScratchSpace& noCellsToReceive,
                       MLongScratchSpace& sendBuffer,
                       MIntScratchSpace& noCellsToSend);
  void communicateDouble(MFloatScratchSpace& recvBuffer,
                         MIntScratchSpace& noCellsToReceive,
                         MFloatScratchSpace& sendBuffer,
                         MIntScratchSpace& noCellsToSend);
  void communicateIntToNeighbors(MIntScratchSpace& recvMem,
                                 MIntScratchSpace& noCellsToReceive,
                                 MIntScratchSpace& sendMem,
                                 MIntScratchSpace& noCellsToSend,
                                 MInt noVar);
  void communicateHaloGlobalIds(MInt level);

 private:
  const MPI_Comm m_mpiComm;
  MInt m_domainId;
  MInt m_noDomains;

  std::ostream outStream;

  Collector<GridgenCell<nDim>>* m_cells = nullptr;
  GridgenCell<nDim>* m_pCells = nullptr;
  MInt m_noCells;
  MLong m_partitionCellOffspringThreshold;
  MFloat m_partitionCellWorkloadThreshold;

  MInt m_initialRefinementLevelSerial;
  MInt m_minLevel;
  SolverRefinement* m_solverRefinement = nullptr;
  MInt m_maxUniformRefinementLevel;
  MInt m_maxRfnmntLvl;
  MInt m_weightMethod;
  MInt m_weightBndCells;
  MInt m_weightPatchCells;
  MBool m_weightSolverUniformLevel;
  MBool m_writeGridInformation;
  MBool m_checkGridLbValidity = true;

  MFloat m_reductionFactor;
  MString m_targetGridFileName = "";
  MString m_outputDir = "";
  MBool m_writeCoordinatesToGridFile = false;
  MInt m_keepOutsideBndryCellChildren;
  MInt* m_noSolidLayer = nullptr;

  // Global bounding box and information for multisolver grid
  MBool m_hasMultiSolverBoundingBox = false;
  std::vector<MFloat> m_multiSolverBoundingBox{};
  MInt m_multiSolverMinLevel = -1;
  MFloat m_multiSolverLengthLevel0 = -1.0;
  std::vector<MFloat> m_multiSolverCenterOfGravity{};

  MInt m_maxNoCells;
  MInt m_maxLevels;
  MString m_gridOutputFileName;
  MInt m_rfnCount;
  MInt m_rfnCountHalos;
  MInt* m_rfnCountHalosDom = nullptr;

  MInt m_noSolvers;
  MInt* m_noBndIdsPerSolver = nullptr;
  MBool** m_bndCutInfo = nullptr;

  MBool m_cutOff;

  MString m_parallelGeomFileName;

  Geometry<nDim>* m_STLgeometry;
  GeometryRoot* m_geometry = nullptr;
  MFloat* m_boundingBox = nullptr;
  MFloat* m_geometryExtents = nullptr;
  MFloat* m_centerOfGravity = nullptr;
  MInt m_decisiveDirection;
  MFloat* m_lengthOnLevel = nullptr;

  MInt** m_levelOffsets = nullptr;
  MInt m_maxNoChildren;
  MInt m_noNeighbors;

  MLong m_noTotalCells;
  MInt m_noPartitionCells;
  MLong m_noTotalPartitionCells;
  MInt m_noTotalHaloCells;
  std::vector<std::tuple<MInt, MLong, MFloat>> m_partitionCellList;
  MInt* m_noCellsPerDomain = nullptr;
  MInt* m_noPartitionCellsPerDomain = nullptr;
  MInt* m_noHaloCellsOnLevel = nullptr;
  MInt m_noNeighborDomains;
  MInt* m_neighborDomains = nullptr;
  MInt** m_haloCellOffsetsLevel = nullptr;
  MLong m_cellOffsetPar;

  MInt** m_haloCellOffsets = nullptr;
  std::map<MInt, MInt> m_cellIdLUT;


  // timers
  MInt m_t_comp_GG = -1;
  MInt m_t_readProperties = -1;
  MInt m_t_initMembers = -1;
  MInt m_t_initGeometry = -1;
  MInt m_t_createInitialGrid = -1;
  MInt m_t_parallelizeGrid = -1;
  MInt m_t_createStartGrid = -1;
  MInt m_t_createComputationalGrid = -1;
  MInt m_t_finalizeGrid = -1;
  MInt m_t_saveGrid = -1;
  MInt m_t_updateInterRankNeighbors = -1;

  // variables used by the dynamic load balancing.
  MInt m_noMissingParents;
  MBool m_hasBeenLoadBalanced;
};
#endif
