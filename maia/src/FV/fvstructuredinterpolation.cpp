// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredinterpolation.h"
#include "COMM/mpioverride.h"
#if defined(WITH_HDF5)
#include "IO/parallelio.h"
#endif
#include "fvstructuredsolver.h"

using namespace std;

template <MInt nDim>
constexpr MInt StructuredInterpolation<nDim>::m_pyramidPoints[72];

/** \brief Constructor for property given donor
 *
 * The constructor for when the grid and the Q file
 * are specified in the property file.
 *
 * \author Marian Albers, Nov 2015
 */
#if defined(WITH_HDF5)
template <MInt nDim>
StructuredInterpolation<nDim>::StructuredInterpolation(const MPI_Comm structuredCommunicator)
  : m_StructuredComm(structuredCommunicator),
    m_donorIsCellCentered(true),
    m_isFieldInterpolation(false),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  MPI_Comm_rank(m_StructuredComm, &m_domainId);

  loadDonorGrid();
  buildDonorTree();
}
#else
template <MInt nDim>
StructuredInterpolation<nDim>::StructuredInterpolation(const MPI_Comm structuredCommunicator)
  : m_StructuredComm(structuredCommunicator),
    m_donorIsCellCentered(true),
    m_isFieldInterpolation(false),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  MPI_Comm_rank(m_StructuredComm, &m_domainId);
  TERMM(1, "Activate HDF5 to use this function!");
}
#endif

/** \brief Constructor for manually given donor
 *
 * The constructor when the grid and the variable field
 * are handed over manually, donorCoordinates are the coordinates
 * at which the donorVariables are located, both have the amount of cells
 * given by noDonorCells
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
StructuredInterpolation<nDim>::StructuredInterpolation(MInt* noDonorCellsDir,
                                                       MFloat** donorCoordinates,
                                                       MFloat** donorVariables,
                                                       const MPI_Comm structuredCommunicator)
  : m_donorCoordinates(donorCoordinates),
    m_donorVariables(donorVariables),
    m_noBlocks(1),
    m_StructuredComm(structuredCommunicator),
    m_noDonorVariables(5),
    m_donorIsCellCentered(true),
    m_isFieldInterpolation(false),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  MPI_Comm_rank(m_StructuredComm, &m_domainId);
  mAlloc(m_noDonorCellsDir, m_noBlocks, nDim, AT_, 0, "m_noDonorCellsDir");
  mAlloc(m_donorBlockOffsets, m_noBlocks, AT_, 0, "m_donorBlockOffset");

  m_totalNoDonorCells = 1;
  for(MInt dim = 0; dim < nDim; dim++) {
    m_noDonorCellsDir[0][dim] = noDonorCellsDir[dim];
    m_totalNoDonorCells *= noDonorCellsDir[dim];
  }

  buildDonorTree();
}

template <MInt nDim>
StructuredInterpolation<nDim>::StructuredInterpolation(MInt* noDonorCellsDir,
                                                       MFloat** donorCoordinates,
                                                       const MPI_Comm structuredCommunicator)
  : m_donorCoordinates(donorCoordinates),
    m_noBlocks(1),
    m_StructuredComm(structuredCommunicator),
    m_noDonorVariables(5),
    m_donorIsCellCentered(true),
    m_isFieldInterpolation(false),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  MPI_Comm_rank(m_StructuredComm, &m_domainId);
  mAlloc(m_noDonorCellsDir, m_noBlocks, nDim, AT_, 0, "m_noDonorCellsDir");
  mAlloc(m_donorBlockOffsets, m_noBlocks, AT_, 0, "m_donorBlockOffset");

  m_totalNoDonorCells = 1;
  for(MInt dim = 0; dim < nDim; dim++) {
    m_noDonorCellsDir[0][dim] = noDonorCellsDir[dim];
    m_totalNoDonorCells *= noDonorCellsDir[dim];
  }

  buildDonorTree();
}


template <MInt nDim>
StructuredInterpolation<nDim>::~StructuredInterpolation() {
  if(m_isFieldInterpolation) {
    mDeallocate(m_donorOriginId);
    mDeallocate(m_transformedReceiverPoints);
    mDeallocate(m_hasInterpolationPartners);
    mDeallocate(m_donorVar);
    mDeallocate(m_hasInterpolationPartnersZonal);
    mDeallocate(m_hasInterpolationPartnersZonalGlobal);
    mDeallocate(m_globalDonorDistanceMin);
  } else {
    mDeallocate(m_donorOriginId);
    mDeallocate(m_interpolationCoefficients);
    mDeallocate(m_hasInterpolationPartners);
    mDeallocate(m_hasInterpolationPartnersZonal);
    mDeallocate(m_hasInterpolationPartnersZonalGlobal);
    mDeallocate(m_globalDonorDistanceMin);
  }

  delete m_donorTree;
}


/** \brief Builds a kd-tree
 * Creates a kd-tree from the predefined grid-data in
 * m_donorCoordinates
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
void StructuredInterpolation<nDim>::buildDonorTree() {
  m_log << "Building up kd-tree..." << endl;

  // cout << "totalNoCells: " << m_totalNoDonorCells << endl;

  // first create kd tree from whole grid
  for(MInt globalId = 0; globalId < m_totalNoDonorCells; globalId++) {
    Point<3> a(m_donorCoordinates[0][globalId], m_donorCoordinates[1][globalId], m_donorCoordinates[2][globalId],
               globalId);
    m_donorPoints.push_back(a);
  }

  m_log << "Created points for kd-tree" << endl;

  // build up the tree and fill it
  m_donorTree = new KDtree<3>(m_donorPoints);

  m_log << "Building up kd-tree... FINISHED!" << endl;
}

/** \brief interpolates variables at point
 * For a given 3D coordinate the method interpolates from
 * the donor grid
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
void StructuredInterpolation<nDim>::interpolateAtPoint(MFloat intPoint[3]) {
  // now go through own, fine cells and look for closest partner neighbour
  MFloat distance = 0;
  Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
  cout << "Finding nearest point..." << endl;
  // find point on the grid that is closest to intPoint
  MInt centerCellId = m_donorTree->nearest(pt, distance);
  cout << "Finding nearest point... FINISHED! Point"
       << " x: " << m_donorCoordinates[0][centerCellId] << " y: " << m_donorCoordinates[1][centerCellId]
       << " z: " << m_donorCoordinates[2][centerCellId] << endl;
  cout << "Finding surrounding hexahedron..." << endl;
  // now eight hexahedron could be candidates for a new home for intPoint,
  // all around centerCellId
  MInt hexahedronOriginId = findSurroundingHexahedron(intPoint, centerCellId, 1);
  MFloat interpolatedVariables[5];

  if(hexahedronOriginId != -1) {
    MFloat transformedPoint[3];
    // now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
    transformPoint(hexahedronOriginId, intPoint, transformedPoint);
    MInt currentBlockId = getBlockId(hexahedronOriginId);
    // interpolate variables at transformed coordinate
    trilinearInterpolation(transformedPoint, hexahedronOriginId, interpolatedVariables, currentBlockId);
  } else {
    // fallback to nearest neighbour interpolation
    nearestNeighbourInterpolation(centerCellId, interpolatedVariables);
  }
}

/** \brief interpolates a field
 * interpolates a given varName and
 * varF
 *
 * \author Marian Albers, Nov 2015
 */
#if defined(WITH_HDF5)
template <MInt nDim>
void StructuredInterpolation<nDim>::interpolateField(MString varName, MFloat* receiverVar) {
  // first load the current variable
  loadDonorVariable(varName);
  for(MInt cellId = 0; cellId < m_noReceiverCells; cellId++) {
    if(cellId % 50000 == 0 && m_domainId == 0) {
      cout << "Variable " << varName
           << " interpolation progress: " << (MInt)((MFloat)cellId / ((MFloat)m_noReceiverCells) * 100.0) << " percent"
           << endl;
    }

    if(m_hasInterpolationPartners[cellId]) {
      MFloat transformedPoint[3] = {m_transformedReceiverPoints[0][cellId], m_transformedReceiverPoints[1][cellId],
                                    m_transformedReceiverPoints[2][cellId]};

      // interpolate variables at transformed coordinate
      MInt currentBlockId = getBlockId(m_donorOriginId[cellId]);
      trilinearInterpolation(transformedPoint, m_donorOriginId[cellId], cellId, receiverVar, currentBlockId);
    } else {
      // fallback to nearest neighbour interpolation
      nearestNeighbourInterpolation(m_donorOriginId[cellId], cellId, receiverVar);
    }
  }
}
#else
template <MInt nDim>
void StructuredInterpolation<nDim>::interpolateField(MString varName, MFloat* receiverVar) {
  (void)varName;
  (void)receiverVar;
  TERMM(1, "Activate HDF5 to use this function!");
}
#endif

/** \brief Prepares interpolation for field
 * For a given 3D coordinate field the method
 * computes the transformed Points and stores them
 * to perform interpolations later
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
void StructuredInterpolation<nDim>::prepareInterpolationField(MInt* noReceiverCells, MFloat** receiverCoordinates) {
  m_isFieldInterpolation = true;
  m_noReceiverCells = noReceiverCells[0] * noReceiverCells[1] * noReceiverCells[2];
  MInt noTrilinear = 0;
  MInt noFallback = 0;
  mAlloc(m_donorOriginId, m_noReceiverCells, AT_, -1, "m_donorOriginId");
  mAlloc(m_transformedReceiverPoints, nDim, m_noReceiverCells, AT_, F0, "m_interpolationCoefficients");
  mAlloc(m_hasInterpolationPartners, m_noReceiverCells, AT_, "m_hasInterpolationPartners");
  mAlloc(m_donorVar, m_totalNoDonorCells, "m_donorVariables", F0, AT_);

  for(MInt cellId = 0; cellId < m_noReceiverCells; cellId++) {
    if(cellId % 50000 == 0) {
      cout << "Interpolation progress: " << (MInt)((MFloat)cellId / ((MFloat)m_noReceiverCells) * 100.0) << " percent"
           << endl;
    }
    m_currentReceiverId = cellId;
    MFloat intPoint[3] = {receiverCoordinates[0][cellId], receiverCoordinates[1][cellId],
                          receiverCoordinates[2][cellId]};
    // now go through own, fine cells and look for closest partner neighbour
    MFloat distance = 0;
    Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
    // find point on the grid that is closest to intPoint
    MInt centerCellId = m_donorTree->nearest(pt, distance);

    // now eight hexahedron could be candidates for a new home for intPoint,
    // all around centerCellId
    MInt hexahedronOriginId = -1;
    const MInt maxNghbrRadius = 4;
    for(MInt nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
      hexahedronOriginId = findSurroundingHexahedron(intPoint, centerCellId, nghbrRadius);
      if(hexahedronOriginId != -1) {
        break;
      }
    }

    if(hexahedronOriginId != -1) {
      MFloat transformedPoint[3] = {F0, F0, F0};
      // now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
      transformPoint(hexahedronOriginId, intPoint, transformedPoint);
      // interpolate variables at transformed coordinate

      for(MInt dim = 0; dim < nDim; dim++) {
        m_transformedReceiverPoints[dim][cellId] = transformedPoint[dim];
      }
      m_donorOriginId[cellId] = hexahedronOriginId;
      m_hasInterpolationPartners[cellId] = true;
      noTrilinear++;
    } else {
      // fallback to nearest neighbour interpolation
      m_hasInterpolationPartners[cellId] = false;
      m_donorOriginId[cellId] = centerCellId;
      noFallback++;
    }
  }

  MInt noLocal[3] = {noTrilinear, noFallback, m_noReceiverCells};
  MInt noGlobal[3] = {0, 0, 0};
  MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_StructuredComm, AT_, "noLocal", "noGlobal");
  if(m_domainId == 0) {
    cout << "Trilinear: " << noGlobal[0] << " (" << 100.0 * ((MFloat)noGlobal[0]) / ((MFloat)noGlobal[2]) << "%) "
         << "Fallback: " << noGlobal[1] << " (" << 100.0 * ((MFloat)noGlobal[1]) / ((MFloat)noGlobal[2]) << "%)"
         << endl;
  }
}

/** \brief Prepares interpolation neighbours and coefficients
 * For a given number of points the methods computes
 * the interpolation partners and coefficients for later use
 * use together with interpolateVariables()
 *
 * \author Marian Albers, Jan 2015
 */
template <MInt nDim>
void StructuredInterpolation<nDim>::prepareInterpolation(MInt noReceiverCells, MFloat** receiverCellCoordinates,
                                                         MInt* interpolationPartner) {
  m_isFieldInterpolation = false;
  MInt noTrilinear = 0;
  MInt noFallback = 0;
  m_noReceiverCells = noReceiverCells;

  if(noReceiverCells > 0) {
    mAlloc(m_donorOriginId, m_noReceiverCells, AT_, "m_donorOriginId");
    mAlloc(m_interpolationCoefficients, m_noReceiverCells, 8, AT_, "m_interpolationCoefficients");
    mAlloc(m_hasInterpolationPartners, m_noReceiverCells, AT_, "m_hasInterpolationPartners");

    for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      if(m_domainId == 0 && receiverId % 10000 == 0) {
        cout << "receiver no: " << receiverId << endl;
      }
      // now go through own, fine cells and look for closest partner neighbour
      MFloat intPoint[3] = {receiverCellCoordinates[0][receiverId], receiverCellCoordinates[1][receiverId],
                            receiverCellCoordinates[2][receiverId]};
      Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
      MFloat dist = 0;
      MInt closestCellId = m_donorTree->nearest(pt, dist);

      // now eight hexahedron could be candidates for a new home for intPoint,
      // all around closestCellId
      // MInt hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, 1);
      MInt hexahedronOriginId = -1;
      const MInt maxNghbrRadius = 8;
      for(MInt nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
        hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, nghbrRadius);
        if(hexahedronOriginId != -1) {
          break;
        }
      }

      if(hexahedronOriginId != -1) {
        m_hasInterpolationPartners[receiverId] = true;
        m_donorOriginId[receiverId] = hexahedronOriginId;
        MFloat transformedPoint[3];
        // now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
        transformPoint(hexahedronOriginId, intPoint, transformedPoint);
        // interpolate variables at transformed coordinate
        computeInterpolationCoefficients(transformedPoint, receiverId);
        noTrilinear++;
      } else {
        // fallback to nearest neighbour interpolation
        m_hasInterpolationPartners[receiverId] = false;
        m_donorOriginId[receiverId] = closestCellId;

        cout << "Fallback x: " << intPoint[0] << " y: " << intPoint[1] << " z: " << intPoint[2] << endl;
        noFallback++;
      }
      interpolationPartner[receiverId] = m_hasInterpolationPartners[receiverId];
    }
  }

  cout << "trilinar: " << noTrilinear << " fallback: " << noFallback << " noReceiverCells: " << noReceiverCells << endl;

  MInt noLocal[3] = {noTrilinear, noFallback, noReceiverCells};
  MInt noGlobal[3] = {0, 0, 0};
  MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_StructuredComm, AT_, "noLocal", "noGlobal");
  if(m_domainId == 0) {
    cout << "Trilinear: " << noGlobal[0] << " (" << 100.0 * ((MFloat)noGlobal[0]) / ((MFloat)noGlobal[2]) << "%)"
         << "Fallback: " << noGlobal[1] << " (" << 100.0 * ((MFloat)noGlobal[1]) / ((MFloat)noGlobal[2]) << "%)"
         << endl;
  }
}

/** \brief computes the interpolation coefficients
 * Interpolates all variables with precomputed interpolation coefficients
 *
 * \author Marian Albers, Jan 2015
 */

template <MInt nDim>
void StructuredInterpolation<nDim>::prepareZonalInterpolation(MInt noReceiverCells,
                                                              MFloat** receiverCellCoordinates,
                                                              MInt* interpolationPartner,
                                                              MBool hasInterpolationPartnerDomain) {
  m_isFieldInterpolation = false;
  MInt noTrilinear = 0;
  MInt noFallback = 0;
  m_noReceiverCells = noReceiverCells;
  m_hasInterpolationPartnerDomain = hasInterpolationPartnerDomain;

  mAlloc(m_donorOriginId, m_noReceiverCells, AT_, "m_donorOriginId");
  mAlloc(m_interpolationCoefficients, m_noReceiverCells, 8, AT_, "m_interpolationCoefficients");
  mAlloc(m_hasInterpolationPartnersZonal, m_noReceiverCells, AT_, 0, "m_hasInterpolationPartnersZonal");
  mAlloc(m_hasInterpolationPartnersZonalGlobal, m_noReceiverCells, AT_, 0, "m_hasInterpolationPartnersZonalGlobal");
  mAlloc(m_donorDistance, m_noReceiverCells, AT_, "m_donorDistance");
  mAlloc(m_globalDonorDistanceMin, m_noReceiverCells, AT_, F0, "m_globalDonorDistanceMin");

  if(m_hasInterpolationPartnerDomain == true) {
    for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      // if(m_domainId==0 && receiverId % 10000 == 0) {cout << "receiver no: " << receiverId << endl;}
      // now go through own, fine cells and look for closest partner neighbour
      MFloat intPoint[3] = {receiverCellCoordinates[0][receiverId], receiverCellCoordinates[1][receiverId],
                            receiverCellCoordinates[2][receiverId]};
      Point<3> pt(intPoint[0], intPoint[1], intPoint[2]);
      MFloat dist = 0;
      MInt closestCellId;
      closestCellId = m_donorTree->nearest(pt, dist);
      // now eight hexahedron could be candidates for a new home for intPoint,
      // all around closestCellId
      // MInt hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, 1);
      MInt hexahedronOriginId = -1;
      const MInt maxNghbrRadius = 4;
      for(MInt nghbrRadius = 1; nghbrRadius < maxNghbrRadius; nghbrRadius++) {
        hexahedronOriginId = findSurroundingHexahedron(intPoint, closestCellId, nghbrRadius);
        if(hexahedronOriginId != -1) {
          break;
        }
      }
      if(hexahedronOriginId != -1) {
        m_hasInterpolationPartnersZonal[receiverId] = true;
        m_donorOriginId[receiverId] = hexahedronOriginId;
        MFloat transformedPoint[3];
        // now hexahedronOriginId is the id of the hexahedron in which intPoint is immersed
        transformPoint(hexahedronOriginId, intPoint, transformedPoint);
        // interpolate variables at transformed coordinate
        computeInterpolationCoefficients(transformedPoint, receiverId);
        noTrilinear++;
      } else {
        // fallback to nearest neighbour interpolation
        m_hasInterpolationPartnersZonal[receiverId] = false;
        m_donorOriginId[receiverId] = closestCellId;
        m_donorDistance[receiverId] = dist;
        noFallback++;
      }
    }
  } else {
    for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      m_donorDistance[receiverId] =
          1000; // avoid to pick interpolation cells in receiverDomain using MPI_Allreduce(MIN) below
      m_hasInterpolationPartnersZonal[receiverId] = false;
    }
  }

  MPI_Allreduce(m_hasInterpolationPartnersZonal, m_hasInterpolationPartnersZonalGlobal, m_noReceiverCells, MPI_INT,
                MPI_SUM, m_StructuredComm, AT_, "m_hasInterpolationPartnersZonal",
                "m_hasInterpolationPartnersZonalGlobal");

  MPI_Allreduce(m_donorDistance, m_globalDonorDistanceMin, m_noReceiverCells, MPI_DOUBLE, MPI_MIN, m_StructuredComm,
                AT_, "m_donorDistance", "m_globalDonorDistanceMin");

  if(m_hasInterpolationPartnerDomain == true) {
    for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
      if(m_hasInterpolationPartnersZonalGlobal[receiverId] == 0
         && approx(m_donorDistance[receiverId], m_globalDonorDistanceMin[receiverId], m_eps)) {
        m_hasInterpolationPartnersZonal[receiverId] = true;
      }
    }
  }

  for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
    interpolationPartner[receiverId] = m_hasInterpolationPartnersZonal[receiverId];
  }


  MInt noLocal[3] = {noTrilinear, noFallback, noReceiverCells};
  MInt noGlobal[3] = {0, 0, 0};
  MPI_Allreduce(noLocal, noGlobal, 3, MPI_INT, MPI_SUM, m_StructuredComm, AT_, "noLocal", "noGlobal");
}


template <MInt nDim>
MFloat StructuredInterpolation<nDim>::interpolateVariableZonal(MFloat* donorVars, MInt cellIdBC) {
  MFloat interpolatedVariable = 0.0;

  if(m_hasInterpolationPartnersZonalGlobal[cellIdBC] > 0) {
    interpolatedVariable = getInterpolatedVariableZonal(donorVars, cellIdBC);
  } else if(m_hasInterpolationPartnersZonalGlobal[cellIdBC] == 0) {
    interpolatedVariable = donorVars[m_donorOriginId[cellIdBC]];
  }

  return interpolatedVariable;
}


template <MInt nDim>
void StructuredInterpolation<nDim>::interpolateVariables(MFloat** interpolatedVariables) {
  for(MInt receiverId = 0; receiverId < m_noReceiverCells; receiverId++) {
    if(m_hasInterpolationPartners[receiverId]) {
      for(MInt var = 0; var < m_noDonorVariables; var++) {
        interpolatedVariables[var][receiverId] = getInterpolatedVariable(receiverId, var);
      }
    } else {
      for(MInt var = 0; var < m_noDonorVariables; var++) {
        interpolatedVariables[var][receiverId] = m_donorVariables[var][m_donorOriginId[receiverId]];
      }
    }
  }
}

/** \brief Finds surrounding hexahedron for point
 * For a given cellId it finds out which of
 * the 8 eight surrounding hexahedrons contains
 * intPoint. If point is inside any of the 6
 * possible tetraeders inside the hexahedron it
 * is inside the hexahedron, the origin cellId
 * is returned then.
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
inline MInt StructuredInterpolation<nDim>::findSurroundingHexahedron(MFloat intPoint[3], MInt centerCellId,
                                                                     MInt stencil) {
  MFloat v1[3] = {F0, F0, F0};
  MFloat v2[3] = {F0, F0, F0};
  MFloat v3[3] = {F0, F0, F0};
  MFloat vp[3] = {F0, F0, F0};
  MFloat vn[3] = {F0, F0, F0};

  MBool isInside = false;

  for(MInt i = -stencil; i < stencil; i++) {
    for(MInt j = -stencil; j < stencil; j++) {
      for(MInt k = -stencil; k < stencil; k++) {
        // find out current blockId
        const MInt currentBlockId = getBlockId(centerCellId);
        const MInt currentOffset = m_donorBlockOffsets[currentBlockId];
        const MInt IJK = getCellIdfromCell(centerCellId, i, j, k, currentBlockId);

        // compute i,j,k from cellId, but use local cellIds
        MInt noLocalCellsDir[3] = {m_noDonorCellsDir[currentBlockId][0] - currentOffset,
                                   m_noDonorCellsDir[currentBlockId][1] - currentOffset,
                                   m_noDonorCellsDir[currentBlockId][2] - currentOffset};
        MInt trueI = (IJK % (noLocalCellsDir[2] * noLocalCellsDir[1])) % noLocalCellsDir[2];
        MInt trueJ = ((IJK - trueI) / noLocalCellsDir[2]) % noLocalCellsDir[1];
        MInt trueK = ((IJK - trueI) / noLocalCellsDir[2] - trueJ) / noLocalCellsDir[1];

        // check if inside regular grid
        if(trueI < 0 || trueI >= noLocalCellsDir[2] - 1 || trueJ < 0 || trueJ >= noLocalCellsDir[1] - 1 || trueK < 0
           || trueK >= noLocalCellsDir[0] - 1) {
          continue;
        }

        // now loop over all 6 possible tetragonal pyramids
        for(MInt tetra = 0; tetra < 6; tetra++) {
          // and now loop over all 4 sides of the pyramid
          isInside = true;
          for(MInt side = 0; side < 4; side++) {
            MInt ixp = (side + 1) % 4;
            MInt ixp2 = (ixp + 1) % 4;
            MInt ixp3 = (ixp2 + 1) % 4;

            // compute vectors
            for(MInt dim = 0; dim < nDim; dim++) {
              v1[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, side, 0), ic(tetra, side, 1),
                                                                  ic(tetra, side, 2), currentBlockId)]
                        - m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp, 0), ic(tetra, ixp, 1),
                                                                    ic(tetra, ixp, 2), currentBlockId)];
              v2[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp2, 0), ic(tetra, ixp2, 1),
                                                                  ic(tetra, ixp2, 2), currentBlockId)]
                        - m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp, 0), ic(tetra, ixp, 1),
                                                                    ic(tetra, ixp, 2), currentBlockId)];
              v3[dim] = m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp3, 0), ic(tetra, ixp3, 1),
                                                                  ic(tetra, ixp3, 2), currentBlockId)]
                        - m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp, 0), ic(tetra, ixp, 1),
                                                                    ic(tetra, ixp, 2), currentBlockId)];
              vp[dim] = intPoint[dim]
                        - m_donorCoordinates[dim][getCellIdfromCell(IJK, ic(tetra, ixp, 0), ic(tetra, ixp, 1),
                                                                    ic(tetra, ixp, 2), currentBlockId)];
            }

            // compute perpendicular vector vn
            crossProduct(vn, v1, v2);

            // check direction of normal vector and correct if necessary
            if(scalarProduct(v3, vn) > m_eps) {
              // invert direction of vn
              for(MInt dim = 0; dim < nDim; dim++) {
                vn[dim] = -vn[dim];
              }
            }

            // check scalar product, if vp*vn > 0 then point is outside of this pyramid
            if(scalarProduct(vp, vn) > m_eps) {
              // not inside pyramid
              isInside = false;
              break;
            }
          }

          // point is inside the current tetra and
          // thus inside the current hexahedron
          if(isInside) {
            return IJK;
          }
        }
      }
    }
  }

  // surprisingly the point isn't inside any of the 8 hexahedrons, return invalid cellId
  return -1;
}

/** \brief Transforms point from physical to computational space
 *
 * Transforms a given point from physical space
 * to computational space with help of the surrounding
 * eight points (which define a hexahedron). Iterative
 * solution with Newton to solve system of equations.
 * Benek (1987) - Chimera - A grid-embedding technique
 *
 * \author Marian Albers, Nov 23, 2015
 */
template <MInt nDim>
inline void StructuredInterpolation<nDim>::transformPoint(MInt hexOrigin, MFloat intPoint[3],
                                                          MFloat transformedPoint[3]) {
  MFloatScratchSpace a(3, 8, AT_, "a");
  const MInt currentBlockId = getBlockId(hexOrigin);

  for(MInt dim = 0; dim < nDim; dim++) {
    a(dim, 0) = m_donorCoordinates[dim][hexOrigin];
    a(dim, 1) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 1, 0, 0, currentBlockId)] - a(dim, 0);
    a(dim, 2) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 0, 1, 0, currentBlockId)] - a(dim, 0);
    a(dim, 3) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 0, 0, 1, currentBlockId)] - a(dim, 0);
    a(dim, 4) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 1, 1, 0, currentBlockId)] - a(dim, 0) - a(dim, 1)
                - a(dim, 2);
    a(dim, 5) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 1, 0, 1, currentBlockId)] - a(dim, 0) - a(dim, 1)
                - a(dim, 3);
    a(dim, 6) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 0, 1, 1, currentBlockId)] - a(dim, 0) - a(dim, 2)
                - a(dim, 3);
    a(dim, 7) = m_donorCoordinates[dim][getCellIdfromCell(hexOrigin, 1, 1, 1, currentBlockId)] - a(dim, 0) - a(dim, 1)
                - a(dim, 2) - a(dim, 3) - a(dim, 4) - a(dim, 5) - a(dim, 6) - a(dim, 7);
  }

  // set initial values 0.5
  MFloatScratchSpace dxezdxyz(3, 3, AT_, "dxezdxyz");
  MFloatScratchSpace rhs(3, AT_, "rhs");

  MFloat xi = F1B2;
  MFloat eta = F1B2;
  MFloat zeta = F1B2;

  MFloat xicor, etacor, zetacor;

  MInt noIterations = 50;
  for(MInt newton = 0; newton < noIterations; newton++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      dxezdxyz(dim, 0) = a(dim, 1) + a(dim, 4) * eta + a(dim, 5) * zeta + a(dim, 7) * eta * zeta;
      dxezdxyz(dim, 1) = a(dim, 2) + a(dim, 4) * xi + a(dim, 6) * zeta + a(dim, 7) * xi * zeta;
      dxezdxyz(dim, 2) = a(dim, 3) + a(dim, 5) * xi + a(dim, 6) * eta + a(dim, 7) * xi * eta;

      rhs[dim] = intPoint[dim]
                 - (a(dim, 0) + a(dim, 1) * xi + a(dim, 2) * eta + a(dim, 3) * zeta + a(dim, 4) * xi * eta
                    + a(dim, 5) * xi * zeta + a(dim, 6) * eta * zeta + a(dim, 7) * xi * eta * zeta);
    }

    MFloat jac = dxezdxyz(0, 0) * (dxezdxyz(1, 1) * dxezdxyz(2, 2) - dxezdxyz(1, 2) * dxezdxyz(2, 1))
                 + dxezdxyz(0, 1) * (dxezdxyz(1, 2) * dxezdxyz(2, 0) - dxezdxyz(1, 0) * dxezdxyz(2, 2))
                 + dxezdxyz(0, 2) * (dxezdxyz(1, 0) * dxezdxyz(2, 1) - dxezdxyz(1, 1) * dxezdxyz(2, 0));

    if(fabs(jac) < m_eps) {
      xicor = F0;
      etacor = F0;
      zetacor = F0;
    } else {
      xicor = (rhs[0] * (dxezdxyz(1, 1) * dxezdxyz(2, 2) - dxezdxyz(1, 2) * dxezdxyz(2, 1))
               + dxezdxyz(0, 1) * (dxezdxyz(1, 2) * rhs[2] - rhs[1] * dxezdxyz(2, 2))
               + dxezdxyz(0, 2) * (rhs[1] * dxezdxyz(2, 1) - dxezdxyz(1, 1) * rhs[2]))
              / jac;
      etacor = (dxezdxyz(0, 0) * (rhs[1] * dxezdxyz(2, 2) - dxezdxyz(1, 2) * rhs[2])
                + rhs[0] * (dxezdxyz(1, 2) * dxezdxyz(2, 0) - dxezdxyz(1, 0) * dxezdxyz(2, 2))
                + dxezdxyz(0, 2) * (dxezdxyz(1, 0) * rhs[2] - rhs[1] * dxezdxyz(2, 0)))
               / jac;
      zetacor = (dxezdxyz(0, 0) * (dxezdxyz(1, 1) * rhs[2] - rhs[1] * dxezdxyz(2, 1))
                 + dxezdxyz(0, 1) * (rhs[1] * dxezdxyz(2, 0) - dxezdxyz(1, 0) * rhs[2])
                 + rhs[0] * (dxezdxyz(1, 0) * dxezdxyz(2, 1) - dxezdxyz(1, 1) * dxezdxyz(2, 0)))
                / jac;
    }

    xi = xi + xicor;
    eta = eta + etacor;
    zeta = zeta + zetacor;

    MFloat sumcor = fabs(xicor) + fabs(etacor) + fabs(zetacor);

    if(sumcor < m_eps) {
      break;
    }
  }

  xi = mMax(mMin(xi, F1), F0);
  eta = mMax(mMin(eta, F1), F0);
  zeta = mMax(mMin(zeta, F1), F0);

  transformedPoint[0] = xi;
  transformedPoint[1] = eta;
  transformedPoint[2] = zeta;
}

/** \brief Trilinear Interpolation for field
 *
 * Performs a trilinear interpolation in computational
 * space for a given point dx. The 8 values are specified
 * at the edges of a uniform cube around point dx
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
inline void StructuredInterpolation<nDim>::trilinearInterpolation(MFloat dx[3], MInt hexOrigin, MInt receiverId,
                                                                  MFloat* receiverVar, MInt currentBlockId) {
  MFloatScratchSpace v(8, AT_, "v");

  v[0] = dx[0] * dx[1] * dx[2];
  v[1] = (F1 - dx[0]) * dx[1] * dx[2];
  v[2] = dx[0] * (F1 - dx[1]) * dx[2];
  v[3] = dx[0] * dx[1] * (F1 - dx[2]);
  v[4] = (F1 - dx[0]) * (F1 - dx[1]) * dx[2];
  v[5] = (F1 - dx[0]) * dx[1] * (F1 - dx[2]);
  v[6] = dx[0] * (F1 - dx[1]) * (F1 - dx[2]);
  v[7] = (F1 - dx[0]) * (F1 - dx[1]) * (F1 - dx[2]);

  // now the actual interpolation
  receiverVar[receiverId] = v[0] * m_donorVar[getCellIdfromCell(hexOrigin, 1, 1, 1, currentBlockId)]
                            + v[1] * m_donorVar[getCellIdfromCell(hexOrigin, 0, 1, 1, currentBlockId)]
                            + v[2] * m_donorVar[getCellIdfromCell(hexOrigin, 1, 0, 1, currentBlockId)]
                            + v[3] * m_donorVar[getCellIdfromCell(hexOrigin, 1, 1, 0, currentBlockId)]
                            + v[4] * m_donorVar[getCellIdfromCell(hexOrigin, 0, 0, 1, currentBlockId)]
                            + v[5] * m_donorVar[getCellIdfromCell(hexOrigin, 0, 1, 0, currentBlockId)]
                            + v[6] * m_donorVar[getCellIdfromCell(hexOrigin, 1, 0, 0, currentBlockId)]
                            + v[7] * m_donorVar[getCellIdfromCell(hexOrigin, 0, 0, 0, currentBlockId)];
}

/** \brief Trilinear interpolation for single point
 *
 * Performs a trilinear interpolation in computational
 * space for a given point dx. The 8 values are specified
 * at the edges of a uniform cube around point dx
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
inline void StructuredInterpolation<nDim>::trilinearInterpolation(MFloat dx[3],
                                                                  MInt hexOrigin,
                                                                  MFloat* receiverVariables,
                                                                  MInt currentBlockId) {
  MFloatScratchSpace v(8, AT_, "v");

  v[0] = dx[0] * dx[1] * dx[2];
  v[1] = (F1 - dx[0]) * dx[1] * dx[2];
  v[2] = dx[0] * (F1 - dx[1]) * dx[2];
  v[3] = dx[0] * dx[1] * (F1 - dx[2]);
  v[4] = (F1 - dx[0]) * (F1 - dx[1]) * dx[2];
  v[5] = (F1 - dx[0]) * dx[1] * (F1 - dx[2]);
  v[6] = dx[0] * (F1 - dx[1]) * (F1 - dx[2]);
  v[7] = (F1 - dx[0]) * (F1 - dx[1]) * (F1 - dx[2]);

  // now the actual interpolation
  for(MInt var = 0; var < m_noDonorVariables; var++) {
    receiverVariables[var] = v[0] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 1, 1, 1, currentBlockId)]
                             + v[1] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 0, 1, 1, currentBlockId)]
                             + v[2] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 1, 0, 1, currentBlockId)]
                             + v[3] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 1, 1, 0, currentBlockId)]
                             + v[4] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 0, 0, 1, currentBlockId)]
                             + v[5] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 0, 1, 0, currentBlockId)]
                             + v[6] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 1, 0, 0, currentBlockId)]
                             + v[7] * m_donorVariables[var][getCellIdfromCell(hexOrigin, 0, 0, 0, currentBlockId)];
  }
}


/** \brief Computes trilinear interpolation coefficients
 *
 * \author Marian Albers, Jan 2016
 */
template <MInt nDim>
inline void StructuredInterpolation<nDim>::computeInterpolationCoefficients(MFloat dx[3], MInt receiverId) {
  m_interpolationCoefficients[receiverId][0] = dx[0] * dx[1] * dx[2];
  m_interpolationCoefficients[receiverId][1] = (F1 - dx[0]) * dx[1] * dx[2];
  m_interpolationCoefficients[receiverId][2] = dx[0] * (F1 - dx[1]) * dx[2];
  m_interpolationCoefficients[receiverId][3] = dx[0] * dx[1] * (F1 - dx[2]);
  m_interpolationCoefficients[receiverId][4] = (F1 - dx[0]) * (F1 - dx[1]) * dx[2];
  m_interpolationCoefficients[receiverId][5] = (F1 - dx[0]) * dx[1] * (F1 - dx[2]);
  m_interpolationCoefficients[receiverId][6] = dx[0] * (F1 - dx[1]) * (F1 - dx[2]);
  m_interpolationCoefficients[receiverId][7] = (F1 - dx[0]) * (F1 - dx[1]) * (F1 - dx[2]);
}

/** \brief Trilinear interpolation for single point
 *
 * Performs trilinear interpolation from previously
 * computed interpolation coefficients
 * (see computeInterpolationCoefficients)
 *
 * \author Marian Albers, Jan 2016
 */
template <MInt nDim>
MFloat StructuredInterpolation<nDim>::getInterpolatedVariable(MInt receiverId, MInt var) {
  const MInt donorOrigin = m_donorOriginId[receiverId];
  const MInt currentBlockId = getBlockId(donorOrigin);
  const MFloat interpolatedVar = m_interpolationCoefficients[receiverId][0]
                                     * m_donorVariables[var][getCellIdfromCell(donorOrigin, 1, 1, 1, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][1]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 0, 1, 1, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][2]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 1, 0, 1, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][3]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 1, 1, 0, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][4]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 0, 0, 1, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][5]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 0, 1, 0, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][6]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 1, 0, 0, currentBlockId)]
                                 + m_interpolationCoefficients[receiverId][7]
                                       * m_donorVariables[var][getCellIdfromCell(donorOrigin, 0, 0, 0, currentBlockId)];

  return interpolatedVar;
}

template <MInt nDim>
MFloat StructuredInterpolation<nDim>::getInterpolatedVariableZonal(MFloat* donorVars, MInt receiverId) {
  const MInt donorOrigin = m_donorOriginId[receiverId];
  const MInt currentBlockId = getBlockId(donorOrigin);
  const MFloat interpolatedVar =
      m_interpolationCoefficients[receiverId][0] * donorVars[getCellIdfromCell(donorOrigin, 1, 1, 1, currentBlockId)]
      + m_interpolationCoefficients[receiverId][1] * donorVars[getCellIdfromCell(donorOrigin, 0, 1, 1, currentBlockId)]
      + m_interpolationCoefficients[receiverId][2] * donorVars[getCellIdfromCell(donorOrigin, 1, 0, 1, currentBlockId)]
      + m_interpolationCoefficients[receiverId][3] * donorVars[getCellIdfromCell(donorOrigin, 1, 1, 0, currentBlockId)]
      + m_interpolationCoefficients[receiverId][4] * donorVars[getCellIdfromCell(donorOrigin, 0, 0, 1, currentBlockId)]
      + m_interpolationCoefficients[receiverId][5] * donorVars[getCellIdfromCell(donorOrigin, 0, 1, 0, currentBlockId)]
      + m_interpolationCoefficients[receiverId][6] * donorVars[getCellIdfromCell(donorOrigin, 1, 0, 0, currentBlockId)]
      + m_interpolationCoefficients[receiverId][7] * donorVars[getCellIdfromCell(donorOrigin, 0, 0, 0, currentBlockId)];

  return interpolatedVar;
}


/**
 * \brief Loads a grid file
 *
 * Loads the grid file given in the property file.
 * No domain decompositioning, all domains read
 * whole grid file as it is necessary for neighbour
 * search.
 *
 * \author Marian Albers
 * \date Nov 2015
 */
#if defined(WITH_HDF5)
template <MInt nDim>
void StructuredInterpolation<nDim>::loadDonorGrid() {
  if(m_domainId == 0) {
    cout << "Reading in donor grid file..." << endl;
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section donorGridName
    <code>MInt donorGridName </code>\n
    default = <code> "" </code>\n \n
    Name of the donor grid file to interpolate from.\n
    Possible values are:\n
    <ul>
    <li>String containing file name</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRUCTURED</i>
  */
  MString donorGridName = Context::getSolverProperty<MString>("donorGridName", 0, AT_);

  MFloat translation[3] = {F0, F0, F0};
  MFloat scale[3] = {F1, F1, F1};

  /*! \property
    \page propertiesFVSTRCTRD
    \section donorTranslation
    <code>MInt translation </code>\n
    default = <code> 0.0, 0.0, 0.0 </code>\n \n
    Translation of the donor grid in 3 space dimensions.\n
    Possible values are:\n
    <ul>
    <li>Float</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRUCTURED</i>
  */
  if(Context::propertyExists("donorTranslation", 0)) {
    for(MInt dim = 0; dim < nDim; dim++) {
      translation[dim] = Context::getSolverProperty<MFloat>("donorTranslation", 0, AT_, &translation[dim], dim);
    }
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section donorScale
    <code>MInt m_donorScale </code>\n
    default = <code> 1.0, 1.0, 1.0 </code>\n \n
    Scaling of the donor grid.\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRUCTURED</i>
  */
  if(Context::propertyExists("donorScale", 0)) {
    for(MInt dim = 0; dim < nDim; dim++) {
      scale[dim] = Context::getSolverProperty<MFloat>("donorScale", 0, AT_, &scale[dim], dim);
    }
  }

  if(m_domainId == 0) {
    cout << "Translating donorGrid by deltaX: " << translation[0] << " deltaY: " << translation[1]
         << " deltaZ: " << translation[2] << endl;
    cout << "Scaling donorGrid by scaleX: " << scale[0] << " scaleY: " << scale[1] << " scaleZ: " << scale[2] << endl;
  }

  ParallelIoHdf5 pio(donorGridName, maia::parallel_io::PIO_READ, m_StructuredComm);
  // create the string to contain the datasetname in the file
  MInt blockId1 = 0;
  MString sBlockName = "/block";
  stringstream dummy1;
  dummy1 << blockId1 << "/";
  sBlockName += dummy1.str();
  m_noDonorDims = pio.getDatasetNoDims("x", sBlockName);

  if(m_domainId == 0) {
    cout << "Donor grid has " << m_noDonorDims << " dimensions" << endl;
  }

  m_noBlocks = -1;
  MInt noBlocksType = pio.getAttributeType("noBlocks", "");
  if(noBlocksType == 1) {
    pio.getAttribute(&m_noBlocks, "noBlocks", "");
  } else if(noBlocksType == 0) {
    MFloat noBlocksFloat = -1.0;
    pio.getAttribute<MFloat>(&noBlocksFloat, "noBlocks", "");
    m_noBlocks = (MInt)noBlocksFloat;
  }

  mAlloc(m_noDonorCellsDir, m_noBlocks, nDim, AT_, 0, "m_noDonorCellsDir");
  mAlloc(m_noDonorPointsDir, m_noBlocks, nDim, AT_, 0, "m_noDonorPointsDir");
  mAlloc(m_noDonorPoints, m_noBlocks, AT_, 1, "m_noDonorPoints");
  mAlloc(m_noDonorCells, m_noBlocks, AT_, 1, "m_noDonorCells");
  m_totalNoDonorCells = 0;
  for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
    MInt ijk_max[MAX_SPACE_DIMENSIONS] = {0};
    stringstream blockName;
    blockName << "/block" << blockId << "/";
    MString blockNameStr = blockName.str();
    std::vector<ParallelIo::size_type> tmp(nDim, 0);
    pio.getArraySize("x", blockNameStr, &tmp[0]);
    std::copy(tmp.begin(), tmp.end(), &ijk_max[0]);

    for(MInt j = 0; j < m_noDonorDims; j++) {
      if(m_donorIsCellCentered == true) {
        m_noDonorCellsDir[blockId][j] = ijk_max[j] - 1;
      } else {
        m_noDonorCellsDir[blockId][j] = ijk_max[j];
      }
      m_noDonorPointsDir[blockId][j] = ijk_max[j];
      // no of points in each solver
      m_noDonorPoints[blockId] *= m_noDonorPointsDir[blockId][j];
      // no of cells in each solver
      m_noDonorCells[blockId] *= m_noDonorCellsDir[blockId][j];
    }
    if(m_domainId == 0) {
      cout << "Donor solver " << blockId << " has " << m_noDonorCells[blockId] << " cells" << endl;
    }
    // total number of cells in all solvers together
    m_totalNoDonorCells += m_noDonorCells[blockId];
  }

  if(m_domainId == 0) {
    cout << "totalNoDonorCells: " << m_totalNoDonorCells << endl;
  }

  if(m_noDonorDims == 2) {
    m_totalNoDonorCells = 0;
    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      m_noDonorPointsDir[blockId][2] = m_noDonorPointsDir[blockId][1];
      m_noDonorPointsDir[blockId][1] = m_noDonorPointsDir[blockId][0];
      m_noDonorPointsDir[blockId][0] = 3;
      m_noDonorCellsDir[blockId][2] = m_noDonorCellsDir[blockId][1];
      m_noDonorCellsDir[blockId][1] = m_noDonorCellsDir[blockId][0];
      m_noDonorCellsDir[blockId][0] = 2;
      m_noDonorCells[blockId] =
          m_noDonorCellsDir[blockId][0] * m_noDonorCellsDir[blockId][1] * m_noDonorCellsDir[blockId][2];
      m_noDonorPoints[blockId] =
          m_noDonorPointsDir[blockId][0] * m_noDonorPointsDir[blockId][1] * m_noDonorPointsDir[blockId][2];
      m_totalNoDonorCells += m_noDonorCells[blockId];
    }
  }

  mAlloc(m_donorBlockOffsets, m_noBlocks, AT_, 0, "m_donorBlockOffset");
  if(m_noBlocks > 1) {
    for(MInt blockId = 1; blockId < m_noBlocks; blockId++) {
      m_donorBlockOffsets[blockId] = m_donorBlockOffsets[blockId - 1] + m_noDonorCells[blockId - 1];
    }
  }

  mAlloc(m_donorCoordinates, nDim, m_totalNoDonorCells, "m_donorGridCoordinates", AT_);
  for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
    MInt offset = m_donorBlockOffsets[blockId];
    stringstream blockName;
    blockName << "/block" << blockId << "/";
    MString blockNameStr = blockName.str();

    if(m_noDonorDims == 3) {
      MFloatScratchSpace gridCoordinates(m_noDonorPoints[blockId], AT_, "gridCoordinates");
      MString coordName[3] = {"x", "y", "z"};
      ParallelIo::size_type ioSize[3] = {m_noDonorPointsDir[blockId][0], m_noDonorPointsDir[blockId][1],
                                         m_noDonorPointsDir[blockId][2]};
      ParallelIo::size_type ioOffset[3] = {0, 0, 0};
      ParallelIo::size_type ioDims = m_noDonorDims;
      for(MInt dim = 0; dim < m_noDonorDims; dim++) {
        pio.readArray(&gridCoordinates[0], blockNameStr, coordName[dim], ioDims, ioOffset, ioSize);
        if(!m_donorIsCellCentered) {
          for(MInt cellId = 0; cellId < m_noDonorCells[blockId]; cellId++) {
            m_donorCoordinates[dim][offset + cellId] = gridCoordinates[cellId];
          }
        } else {
          computeCellCentreCoordinates(m_noDonorPointsDir[blockId], gridCoordinates, dim, blockId);
        }
        for(MInt cellId = 0; cellId < m_noDonorCells[blockId]; cellId++) {
          m_donorCoordinates[dim][offset + cellId] =
              translation[dim] + scale[dim] * m_donorCoordinates[dim][offset + cellId];
        }
      }
    } else {
      ParallelIo::size_type ioSize[2] = {m_noDonorPointsDir[blockId][1], m_noDonorPointsDir[blockId][2]};
      ParallelIo::size_type ioOffset[3] = {0, 0, 0};
      MFloatScratchSpace gridCoordinates(m_noDonorPoints[blockId], AT_, "gridCoordinates");

      for(MInt dim = 0; dim < nDim; dim++) {
        MString coordName[3] = {"x", "y", "z"};
        if(dim < 2) {
          pio.readArray(gridCoordinates.begin(), blockNameStr, coordName[dim], m_noDonorDims, ioOffset, ioSize);

          // duplicate the values in z-direction
          for(MInt i = 0; i < m_noDonorPointsDir[blockId][2]; i++) {
            for(MInt j = 0; j < m_noDonorPointsDir[blockId][1]; j++) {
              for(MInt k = 1; k < m_noDonorPointsDir[blockId][0]; k++) {
                MInt pointId2D = i + j * m_noDonorPointsDir[blockId][2];
                MInt pointId = i + (j + k * m_noDonorPointsDir[blockId][1]) * m_noDonorPointsDir[blockId][2];
                gridCoordinates(pointId) = gridCoordinates(pointId2D);
              }
            }
          }
        } else {
          MFloat minZ = F0;
          MFloat maxZ = F1;
          for(MInt i = 0; i < m_noDonorPointsDir[blockId][2]; i++) {
            for(MInt j = 0; j < m_noDonorPointsDir[blockId][1]; j++) {
              for(MInt k = 0; k < m_noDonorPointsDir[blockId][0]; k++) {
                MInt pointId = i + (j + k * m_noDonorPointsDir[blockId][1]) * m_noDonorPointsDir[blockId][2];
                gridCoordinates(pointId) = (maxZ - minZ) / (m_noDonorPointsDir[blockId][0] - 1) * k;
              }
            }
          }
        }

        computeCellCentreCoordinates(m_noDonorPointsDir[blockId], gridCoordinates, dim, blockId);
        for(MInt cellId = 0; cellId < m_noDonorCells[blockId]; cellId++) {
          m_donorCoordinates[dim][cellId] = translation[dim] + scale[dim] * m_donorCoordinates[dim][offset + cellId];
        }
      }
    }
  }

  MBool write2D3DGrid = false;
  if(write2D3DGrid) {
    const char* fileName = "rescaledGrid.hdf5";
    ParallelIoHdf5 rescaledGridFile(fileName, maia::parallel_io::PIO_REPLACE, m_StructuredComm);
    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      MInt offset = m_donorBlockOffsets[blockId];
      stringstream blockName;
      blockName << "/block" << blockId << "/";
      MString blockNameStr = blockName.str();

      ParallelIo::size_type ioSize[3] = {m_noDonorCellsDir[blockId][0], m_noDonorCellsDir[blockId][1],
                                         m_noDonorCellsDir[blockId][2]};
      rescaledGridFile.defineArray(maia::parallel_io::PIO_FLOAT, blockNameStr, "x", 3, ioSize);
      rescaledGridFile.defineArray(maia::parallel_io::PIO_FLOAT, blockNameStr, "y", 3, ioSize);
      rescaledGridFile.defineArray(maia::parallel_io::PIO_FLOAT, blockNameStr, "z", 3, ioSize);
      ParallelIo::size_type ioOffset[3] = {0, 0, 0};

      if(m_domainId == 0) {
        rescaledGridFile.writeArray(&m_donorCoordinates[0][offset + 0], blockNameStr, "x", nDim, ioOffset, ioSize);
        rescaledGridFile.writeArray(&m_donorCoordinates[1][offset + 0], blockNameStr, "y", nDim, ioOffset, ioSize);
        rescaledGridFile.writeArray(&m_donorCoordinates[2][offset + 0], blockNameStr, "z", nDim, ioOffset, ioSize);
      } else {
        ParallelIo::size_type ioEmptySize[3] = {0, 0, 0};
        MFloat empty = 0;
        rescaledGridFile.writeArray(&empty, blockNameStr, "x", nDim, ioOffset, ioEmptySize);
        rescaledGridFile.writeArray(&empty, blockNameStr, "y", nDim, ioOffset, ioEmptySize);
        rescaledGridFile.writeArray(&empty, blockNameStr, "z", nDim, ioOffset, ioEmptySize);
      }
    }
  }
  if(m_domainId == 0) {
    cout << "Reading in donor grid file... FINISHED!" << endl;
  }
}
#endif

/** \brief Loads a Q file
 *
 * loads a given variable (with varName) from donorFile
 *
 * \author Marian Albers, Nov 2015
 */
#if defined(WITH_HDF5)
template <MInt nDim>
void StructuredInterpolation<nDim>::loadDonorVariable(MString varName) {
  stringstream donorFileName;

  if(m_domainId == 0) {
    cout << "Reading in " << varName << " from donor file..." << endl;
  }
  // reading in a specified donor file from the properties file

  /*! \property
    \page propertiesFVSTRCTRD
    \section donorVars
    <code>MInt donorFile </code>\n
    default = <code> "" </code>\n \n
    Name of the donor var file to interpolate from.\n
    Possible values are:\n
    <ul>
    <li>String containing file name</li>
    </ul>
    Keywords: <i>INTERPOLATION, STRUCTURED</i>
  */
  MString donorFile = Context::getSolverProperty<MString>("donorVars", 0, AT_);

  // open the file
  ParallelIoHdf5 pio(donorFile, maia::parallel_io::PIO_READ, m_StructuredComm);

  for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
    stringstream blockName;
    blockName << "/block" << blockId << "/";
    MString blockNameStr = blockName.str();

    MBool fieldExists = pio.hasDataset(varName, blockNameStr);
    if(fieldExists) {
      if(m_domainId == 0) {
        cout << "Field " << varName << " exist, loading from donorFile. " << endl;
      }
      MInt offset = m_donorBlockOffsets[blockId];
      ParallelIo::size_type ioSize[3] = {0, 0, 0};
      ParallelIo::size_type ioOffset[3] = {0, 0, 0};
      if(m_noDonorDims == 3) {
        ioSize[0] = m_noDonorCellsDir[blockId][0];
        ioSize[1] = m_noDonorCellsDir[blockId][1];
        ioSize[2] = m_noDonorCellsDir[blockId][2];
      } else {
        ioSize[0] = m_noDonorCellsDir[blockId][1];
        ioSize[1] = m_noDonorCellsDir[blockId][2];
      }

      pio.readArray(&m_donorVar[offset], blockNameStr, varName, m_noDonorDims, ioOffset, ioSize);

      // enlarge 2D field to 3D
      if(m_noDonorDims == 2) {
        for(MInt k = 1; k < m_noDonorCellsDir[blockId][0]; k++) {
          for(MInt j = 0; j < m_noDonorCellsDir[blockId][1]; j++) {
            for(MInt i = 0; i < m_noDonorCellsDir[blockId][2]; i++) {
              const MInt cellId = cellIndex(i, j, k, blockId);
              const MInt cellId2D = cellIndex(i, j, 0, blockId);
              m_donorVar[cellId] = m_donorVar[cellId2D];
            }
          }
        }
      }
    } else {
      if(m_domainId == 0) {
        cout << "Field " << varName << " doesn't exist, setting zero. " << endl;
      }
      for(MInt k = 0; k < m_noDonorCellsDir[blockId][0]; k++) {
        for(MInt j = 0; j < m_noDonorCellsDir[blockId][1]; j++) {
          for(MInt i = 0; i < m_noDonorCellsDir[blockId][2]; i++) {
            const MInt cellId = cellIndex(i, j, k, blockId);
            m_donorVar[cellId] = F0;
          }
        }
      }
    }
  }

  if(m_domainId == 0) {
    cout << "Reading in " << varName << " from donor file... FINISHED!" << endl;
  }
}
#endif

/** \brief Computes cell centered coordinates
 *
 * Computes the cell-centered donorCoordinates from
 * node-centered grid
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
void StructuredInterpolation<nDim>::computeCellCentreCoordinates(MInt* noPoints,
                                                                 MFloatScratchSpace& coordinates,
                                                                 MInt dim,
                                                                 MInt blockId) {
  // function to compute the coordinates at cell centre
  // do it over I, J, K loop but change to one array

  for(MInt k = 0; k < m_noDonorCellsDir[blockId][0]; k++) {
    for(MInt j = 0; j < m_noDonorCellsDir[blockId][1]; j++) {
      for(MInt i = 0; i < m_noDonorCellsDir[blockId][2]; i++) {
        MInt pointId = i + (j + k * noPoints[1]) * noPoints[2];
        MInt IJK = pointId;
        MInt IP1JK = pointId + 1;
        MInt IJP1K = pointId + noPoints[2];
        MInt IP1JP1K = IJP1K + 1;
        MInt IJKP1 = pointId + noPoints[2] * noPoints[1];
        MInt IP1JKP1 = IJKP1 + 1;
        MInt IJP1KP1 = pointId + noPoints[2] + noPoints[2] * noPoints[1];
        MInt IP1JP1KP1 = IJP1KP1 + 1;
        MInt cellId = cellIndex(i, j, k, blockId);

        // average the coordinates for cell centre data
        m_donorCoordinates[dim][cellId] =
            F1B8
            * (coordinates[IJK] + coordinates[IP1JK] + coordinates[IJP1K] + coordinates[IP1JP1K] + coordinates[IJKP1]
               + coordinates[IP1JKP1] + coordinates[IJP1KP1] + coordinates[IP1JP1KP1]);
      }
    }
  }
}

/** \brief Nearest neighbour interpolation
 *
 * Fallback routine for when no trilinear interpolation
 * is possible. Takes values from nearest neighbour.
 *
 * \author Marian Albers, Nov 2015
 */
template <MInt nDim>
inline void StructuredInterpolation<nDim>::nearestNeighbourInterpolation(MInt donorId, MInt receiverId,
                                                                         MFloat* receiverVar) {
  receiverVar[receiverId] = m_donorVar[donorId];
}

template <MInt nDim>
inline void StructuredInterpolation<nDim>::nearestNeighbourInterpolation(MInt donorId, MFloat* receiverVariables) {
  for(MInt var = 0; var < m_noDonorVariables; var++) {
    receiverVariables[var] = m_donorVariables[var][donorId];
  }
}

template <MInt nDim>
inline void StructuredInterpolation<nDim>::crossProduct(MFloat result[3], MFloat vec1[3], MFloat vec2[3]) {
  result[xsd] = vec1[ysd] * vec2[zsd] - vec1[zsd] * vec2[ysd];
  result[ysd] = vec1[zsd] * vec2[xsd] - vec1[xsd] * vec2[zsd];
  result[zsd] = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
}

template <MInt nDim>
inline MFloat StructuredInterpolation<nDim>::scalarProduct(MFloat vec1[3], MFloat vec2[3]) {
  return vec1[xsd] * vec2[xsd] + vec1[ysd] * vec2[ysd] + vec1[zsd] * vec2[zsd];
}

template <MInt nDim>
inline MInt StructuredInterpolation<nDim>::getCellIdfromCell(MInt origin, MInt incI, MInt incJ, MInt incK,
                                                             MInt blockId) {
  return origin + incI + incJ * m_noDonorCellsDir[blockId][2]
         + incK * m_noDonorCellsDir[blockId][2] * m_noDonorCellsDir[blockId][1];
}

template <MInt nDim>
inline MInt StructuredInterpolation<nDim>::cellIndex(MInt i, MInt j, MInt k, MInt blockId) {
  const MInt offset = m_donorBlockOffsets[blockId];
  return offset + i + (j + k * m_noDonorCellsDir[blockId][1]) * m_noDonorCellsDir[blockId][2];
}

// returns the ijk-coordinate for one of the 4 side
// of one of the 6 tetraeders inside a hexahedron
template <MInt nDim>
inline MInt StructuredInterpolation<nDim>::ic(MInt tetra, MInt side, MInt dim) {
  return m_pyramidPoints[dim + side * 3 + tetra * 3 * 4];
}

template <MInt nDim>
inline MInt StructuredInterpolation<nDim>::getBlockId(MInt cellId) {
  MInt currentBlockId = 0;
  if(m_noBlocks > 1) {
    for(MInt blockId = 0; blockId < m_noBlocks; blockId++) {
      if(cellId < m_donorBlockOffsets[blockId + 1]) {
        currentBlockId = blockId;
        break;
      }
    }
  }

  return currentBlockId;
}

template <MInt nDim>
inline MBool StructuredInterpolation<nDim>::approx(const MFloat& a, const MFloat& b, const MFloat eps) {
  return abs(a - b) < eps;
}


// Explicit instantiations for 2D and 3D
template class StructuredInterpolation<2>;
template class StructuredInterpolation<3>;
