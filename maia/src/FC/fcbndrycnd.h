// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FCBNDRYCND_H
#define FCBNDRYCND_H

#include <map>
#include <vector>
#include "GEOM/geometryelement.h"
#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangridcellproperties.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/parallelio.h"
#include "compiler_config.h"
#include "fccellcollector.h"
#include "fcdescriptor.h"

template <MInt nDim>
class FcSolver;
template <class T>
class Collector;
template <MInt nDim>
class CartesianGrid;
template <MInt nDim>
class FcGridBndryCell;


template <MInt nDim>
class FcBndryCnd {
 public:
  template <MInt nDim_>
  friend class FcSolver;

  /// Types
  // Type for cell properties
  using Cell = GridCell;
  using SolverCell = FcCell;
  using Fd = FcDescriptor<nDim>;

 protected:
  using FcCellCollector = maia::fc::collector::FcCellCollector<nDim>;
  std::vector<FcGridBndryCell<nDim>> m_bndryCells;

  //! pointer to a member function data type
  FcSolver<nDim>* m_solver;
  FcGridBndryCell<nDim>* bndryCells;
  typedef void (FcBndryCnd::*BndryCndHandler)(MInt set);
  BndryCndHandler* bndryCndHandlerSystemMatrix;
  BndryCndHandler* bndryCndHandlerForce;

  MInt m_solverId;
  MInt m_noInternalCells;

  std::vector<MInt> m_bndryCndIds;     // Holds the different BC ids
  std::vector<MInt> m_bndryCndOffsets; // stores the starting positions

  std::vector<MInt> m_bndryCndSegIds;         // Holds the different segment ids
  std::vector<MInt> m_mapBndryCndSegId2Index; // maps global segmentId to local index (-1 if not available)

  MInt* m_mapBndryCndIdSegId{};

  MInt m_noSegments{}; // the number of all segments
  MInt* m_noBndryCellsPerSegment = nullptr;

  /// Read from property file
  // Basic
  MString m_bndryNormalMethod;
  MBool m_multiBC;

  // TIMERS
  MInt m_t_BCAll;

  MInt* m_subCellLayerDepth = nullptr;
  MFloat m_kFactor = 1e+12;

 public:
  FcBndryCnd(FcSolver<nDim>* solver);

  virtual ~FcBndryCnd();

  List<MInt>* m_bndryCellIds = nullptr;
  List<MInt>* m_sortedBndryCells = nullptr;
  MInt* m_boundarySurfaces = nullptr;

  void initBndryCnds();
  void setBndryCndHandler();
  void setBndryCndHandlerMb(MInt* bcTypeList);
  void createBoundaryCells();
  void sortBoundaryCells();
  void updateSystemMatrix();
  void updateForceVector();
  void updateTrianglePosition();
  void calcReactionForces();
  void writeOutModifiedBoundaries();
  void calculateCutPoints();
  void calcAvgFaceNormal(MInt index);
  void setBoundaryStlElements(std::vector<std::vector<MFloat>> cutPointList, MFloat* normal, MInt segId, MInt i);
  void writeOutSTLBoundaries(MInt index, MString prefix);
  void subCellIntegration(
      const MInt subCellLvl, MFloat* subCellParentCoord, const MInt childPos, const MInt pCellId, MFloat** Ke);
  void findCellsRequireSubCellIntegration();

  std::vector<std::vector<MFloat>> createTrianglesFromCutPoints(std::vector<std::vector<MFloat>> cutPointList,
                                                                MFloat* triangleNormal);
  std::vector<std::vector<MFloat>> createEdgesFromCutPoints(std::vector<std::vector<MFloat>> cutPointList);
  MFloat solveIntegrationOnTriangle(std::vector<MFloat> trianglePoints, MFloat* triangleNormal);
  MBool pointInTriangle(MFloat* A, MFloat* B, MFloat* C, std::vector<MFloat> P);
  void refineTriangle(MInt index);
  void getStlNodeList(const MInt cellId, std::vector<MInt>& nodeList);
  // Specific boundary conditions
  // Each boundary condition has additional properties,
  // they carry the number of the boundary condition
  // in their names (BC.0, BC.1, ...)

  // Dummy bc
  virtual void bc0(MInt index);

  // Fixation bc
  void bc8010(MInt index); // constraint fixation
  void bc8011(MInt index); // constraint fixation
  void bc8012(MInt index); // constraint fixation

  // Displacement bc
  void bc8020(MInt index); // constraint displacement

  // Surface tensions bc
  void bc8030(MInt index); // imposed loads

  // Surface tensions bc
  void bc8031(MInt index); // imposed loads
  void bc8032(MInt index); // imposed loads
  void bc8035(MInt index); // imposed loads

  // Surface tensions bc
  void bc8040(MInt index); // imposed loads
};

#endif
