// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDINTERPOLATION_H
#define STRUCTUREDINTERPOLATION_H

#include <vector>
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"
#include "globals.h"

template <MInt nDim>
class StructuredInterpolation {
 public:
  StructuredInterpolation(const MPI_Comm structuredCommunicator);
  StructuredInterpolation(MInt* noDonorCells, MFloat** donorCoordinates, MFloat** donorVariables,
                          const MPI_Comm structuredCommunicator);
  StructuredInterpolation(MInt* noDonorCellsDir, MFloat** donorCoordinates, const MPI_Comm structuredCommunicator);
  ~StructuredInterpolation();
  void interpolateAtPoint(MFloat* intPoint);
  void prepareInterpolationField(MInt* noReceiverCells, MFloat** receiverCoordinates);
  void interpolateField(MString, MFloat*);
  void loadDonorGrid();
  void loadDonorVariable(MString varName);
  void prepareInterpolation(MInt, MFloat**, MInt*);
  void interpolateVariables(MFloat**);
  MFloat getInterpolatedVariable(MInt, MInt);
  void prepareZonalInterpolation(MInt, MFloat**, MInt*, MBool);
  MFloat interpolateVariableZonal(MFloat*, MInt);
  MFloat getInterpolatedVariableZonal(MFloat*, MInt);

 protected:
  void buildDonorTree();
  inline void crossProduct(MFloat result[3], MFloat vec1[3], MFloat vec2[3]);
  inline MFloat scalarProduct(MFloat vec1[3], MFloat vec2[3]);
  inline MInt getCellIdfromCell(MInt origin, MInt incI, MInt incJ, MInt incK, MInt solverId);
  inline MBool approx(const MFloat&, const MFloat&, const MFloat);
  inline MInt getBlockId(MInt cellId);
  inline void trilinearInterpolation(MFloat*, MInt, MInt, MFloat*, MInt);
  inline void trilinearInterpolation(MFloat*, MInt, MFloat*, MInt);

  inline void computeInterpolationCoefficients(MFloat*, MInt);
  inline void transformPoint(MInt hexOrigin, MFloat intPoint[3], MFloat transformedPoint[3]);
  inline MInt findSurroundingHexahedron(MFloat intPoint[3], MInt centerCellId, MInt stencil);
  inline void nearestNeighbourInterpolation(MInt, MInt, MFloat*);
  inline void nearestNeighbourInterpolation(MInt, MFloat*);
  inline MInt cellIndex(MInt i, MInt j, MInt k, MInt solverId);
  inline MInt ic(MInt, MInt, MInt);
  void computeCellCentreCoordinates(MInt*, MFloatScratchSpace&, MInt, MInt);

  // index variables
  static const MInt xsd = 0;
  static const MInt ysd = 1;
  static const MInt zsd = 2;

  MInt** m_noDonorCellsDir = nullptr;
  MInt** m_noDonorPointsDir = nullptr;
  MInt* m_noDonorCells = nullptr;
  MInt* m_noDonorPoints = nullptr;
  MFloat** m_donorCoordinates = nullptr;
  MFloat** m_donorVariables = nullptr;
  MFloat* m_donorVar = nullptr;
  MInt m_noDonorDims;
  MInt m_totalNoDonorCells;
  MInt* m_donorBlockOffsets = nullptr;
  MInt m_noBlocks;

  const MPI_Comm m_StructuredComm;
  MInt m_domainId;
  MInt m_noDonorVariables;
  MBool m_donorIsCellCentered;
  MBool m_isFieldInterpolation;
  MFloat m_eps;

  MInt* m_donorOriginId = nullptr;                // needed if donorIds are saved for multiple interpolation
  MFloat** m_interpolationCoefficients = nullptr; // needed if coefficients are saved for multiple interpolation
  MFloat* m_donorDistance = nullptr;
  MFloat* m_globalDonorDistanceMin = nullptr;
  MBool m_hasInterpolationPartnerDomain;
  MInt* m_hasInterpolationPartnersZonal = nullptr;
  MInt* m_hasInterpolationPartnersZonalGlobal = nullptr;

  MFloat** m_transformedReceiverPoints = nullptr;
  MFloat** m_receiverVariables = nullptr;
  MBool* m_hasInterpolationPartners = nullptr;
  MInt m_noReceiverCells;
  MInt m_currentReceiverId;

  std::vector<Point<3>> m_donorPoints;
  KDtree<3>* m_donorTree;

  // pyramidPoints contains the ijk-combinations
  // for all possible tetraeders
  // inside a hexahedron, use together with function ic(tetraeder,side,dim)
  static constexpr MInt m_pyramidPoints[72] = {0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0,
                                               1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1,
                                               1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0};
};

#endif
