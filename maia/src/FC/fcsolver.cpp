// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fcsolver.h"

#include <cstring>
#include <functional>
#include <random>
#include <unordered_map>
#include <utility>
#include "COMM/mpioverride.h"
#include "GEOM/geometryelement.h"
#include "IO/parallelio.h"
#include "UTIL/maiamath.h"
#include "fcgridbndrycell.h"
#include "globals.h"

using namespace std;

template <MInt nDim>
FcSolver<nDim>::FcSolver(MInt id, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm)
  : maia::CartesianSolver<nDim, FcSolver<nDim>>(id, gridProxy_, comm), m_geometry(&geometry_) {
  TRACE();

  initTimer();
  RECORD_TIMER_START(m_t.solver);

  m_geometryIntersection = new GeometryIntersection<nDim>(&grid(), m_geometry);

  /*! \page propertyPage1
    \section isThermal
    <code>MBool FcSolver::m_isThermal</code>\n
    default = <code>false</code>\n\n
    Do thermal stresses exist.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_isThermal = false;
  m_isThermal = Context::getSolverProperty<MBool>("isThermal", m_solverId, AT_, &m_isThermal);
  m_cells.setThermal(m_isThermal);

  /*! \page propertyPage1
    \section polyDegree
    <code>MInt FcSolver::m_polyDeg</code>\n
    default = <code>0</code>\n\n
    Specifies the polynominal degree\n
    Keywords: <i>FINITE CELL</i>
  */
  m_polyDeg = 0;
  m_polyDeg = Context::getSolverProperty<MInt>("polyDeg", m_solverId, AT_, &m_polyDeg);
  m_cells.setMaxPRfnmnt(m_polyDeg);

  /*! \page propertyPage1
    \section testRun
    <code>MBool FcSolver::m_testRun</code>\n
    default = <code>false</code>\n\n
    Enables extra debug output.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_testRun = false;
  m_testRun = Context::getSolverProperty<MBool>("testRun", m_solverId, AT_, &m_testRun);

  if(m_testRun) {
    /*! \page propertyPage1
      \section analyticSolution
      <code>MFloat FcSolver::m_analyticSolution</code>\n
      default = <code>F1</code>\n\n
      Keywords: <i>FINITE CELL</i>
    */
    m_analyticSolution = F1;
    m_analyticSolution = Context::getSolverProperty<MFloat>("analyticSolution", m_solverId, AT_, &m_analyticSolution);

    /*! \page propertyPage1
      \section printEigenValues
      <code>MBool FcSolver::m_printEigenValues</code>\n
      default = <code>false</code>\n\n
      Keywords: <i>FINITE CELL</i>
    */
    m_printEigenValues = false;
    m_printEigenValues = Context::getSolverProperty<MBool>("printEigenValues", m_solverId, AT_, &m_printEigenValues);
  }

  /*! \page propertyPage1
    \section EModule
    <code>MFloat FcSolver::m_E</code>\n
    default = <code>100000.0</code>\n\n
    Defines the E-Module of the material.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_E = 100000.0;
  m_E = Context::getSolverProperty<MFloat>("EModule", m_solverId, AT_, &m_E);

  /*! \page propertyPage1
    \section epsBiCG
    <code>MFloat FcSolver::m_eps</code>\n
    default = <code>1e-12</code>\n\n
    Defines the limit for the iterative solver.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_eps = 1e-12;
  m_eps = Context::getSolverProperty<MFloat>("epsBiCG", m_solverId, AT_, &m_eps);

  /*! \page propertyPage1
    \section alpha
    <code>MFloat FcSolver::m_alpha</code>\n
    default = <code>1e-14</code>\n\n
    Defines the penalty factor for sub cell integration.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_alpha = 1e-14;
  m_alpha = Context::getSolverProperty<MFloat>("alpha", m_solverId, AT_, &m_alpha);

  /*! \page propertyPage1
    \section noIterations
    <code>MInt FcSolver::m_maxNoIterations</code>\n
    default = <code>10000</code>\n\n
    Specifies the maximum number of iterations.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_maxNoIterations = 10000;
  m_maxNoIterations = Context::getSolverProperty<MInt>("noIterations", m_solverId, AT_, &m_maxNoIterations);

  /*! \page propertyPage1
    \section noLoadSteps
    <code>MInt FcSolver::m_noLoadSteps</code>\n
    default = <code>10</code>\n\n
    Specifies the maximum number of sub steps for iterative calculation of displacements.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_noLoadSteps = 1;
  m_noLoadSteps = Context::getSolverProperty<MInt>("noLoadSteps", m_solverId, AT_, &m_noLoadSteps);

  /*! \page propertyPage1
    \section solveSoEIteratively
    <code>MInt FcSolver::m_solveSoEIteratively</code>\n
    default = <code>true</code>\n\n
    Specifies if the SoE is solved iteratively using BiCGStab or directly using Eigen SuperLU.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_solveSoEIteratively = true;
  m_solveSoEIteratively =
      Context::getSolverProperty<MBool>("solveSoEIteratively", m_solverId, AT_, &m_solveSoEIteratively);

  /*! \page propertyPage1
    \section fcInterpolationMethod
    <code>MString FcSolver::m_fcInterpolationMethod</code>\n
    default = <code>LAGRANGE_INTERP</code>\n\n
    Specifies interpolation method used.\n
    Keywords: <i>FINITE CELL</i>
  */
  m_fcInterpolationMethod = "LAGRANGE_INTERP";
  m_fcInterpolationMethod =
      Context::getSolverProperty<MString>("fcInterpolationMethod", m_solverId, AT_, &m_fcInterpolationMethod);

  // Instead of noCells() grid().tree().capacity() was used in m_cells.reset() before;
  // it was changed because there was no difference at that time and to avoid direct access of the gird()
  //  m_cells.reset(grid().noCells());
  m_cells.reset(grid().raw().treeb().capacity());
  m_cells.append(grid().noCells());

  // Update cell collector with information from grid
  updateCellCollectorFromGrid();

  setIntegrationWeights();

  // TODO: Is this necessary
  initJacobianMatrix();

  // Find also the diagonal neighbors
  grid().findEqualLevelNeighborsParDiagonal(false);

  // Instantiation of bndrycnd
  m_bndryCnd = new FcBndryCnd<nDim>(this);

  deactivateCells();

  // Allocate memory for the displacement vector
  allocateNodeVectors();
}

/** \brief Destructor
 *
 * \author Moritz Waldmann
 * \date 20.02.21
 *
 **/
template <MInt nDim>
FcSolver<nDim>::~FcSolver() {
  TRACE();
  mDeallocate(m_nodalDisplacements);
  mDeallocate(m_totalNodalDisplacements);
  mDeallocate(m_sysMatCompressed);
  mDeallocate(m_compressionIndexSys);

  mDeallocate(m_nodalLoadVector);
  mDeallocate(m_externalLoadVector);
  mDeallocate(m_internalLoadVector);
  mDeallocate(m_reactionForceVector);

  RECORD_TIMER_STOP(m_t.solver);

  averageTimer();
}

/** \brief Set halo cells, poissonRatio polynom degree and
 * the number of elements for each cell.
 * Furthermore, the lobatto weights and the lobatto points are set
 *
 * \author Moritz Waldmann
 * \date 20.02.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::updateCellCollectorFromGrid() {
  TRACE();
  // Store level information in cell collector for faster access:
  // halo cells are stored,
  // the poisson ratio is stored in case it differs from cell to cell,
  // the p-Refinement per cell is stored,
  // the number of points per cell is stored,
  for(MInt i = 0; i < a_noCells(); i++) {
    a_isHalo(i) = false;
    a_isWindow(i) = false;
    a_poissonRatio(i) = m_poissonRatio;
    a_pRfnmnt(i) = m_polyDeg;
    a_noNodes(i) = (a_pRfnmnt(i) + 2) * (a_pRfnmnt(i) + 2);
    a_isActive(i) = true;
    if(nDim == 3) a_noNodes(i) *= (a_pRfnmnt(i) + 2);

    for(MInt strain = 0; strain < m_noStrains; strain++) {
      a_elementStrains(i, strain) = F0;
    }

    for(MInt stress = 0; stress < m_noStrains; stress++) {
      a_elementStresses(i, stress) = F0;
    }

    switch(string2enum(m_fcInterpolationMethod)) {
      case EQUIDIST_LAGRANGE_INTERP: {
        for(MInt deg = 0; deg < a_pRfnmnt(i) + 2; deg++) {
          a_nodePosition(i, deg) = Fd::nodePosLobattoPoints(a_pRfnmnt(i), deg - 1);
        }
        break;
      }
      case LAGRANGE_INTERP: {
        for(MInt deg = 0; deg < a_pRfnmnt(i) + 2; deg++) {
          if(deg == 0) {
            a_nodePosition(i, deg) = -F1;
          } else if(deg < a_pRfnmnt(i) + 1) {
            a_nodePosition(i, deg) = Fd::nodePosLobattoPoints(a_pRfnmnt(i), deg - 1);
          } else {
            a_nodePosition(i, deg) = F1;
          }
        }
        break;
      }
      case EQUIDIST_LEGENDRE_INTERP: {
        for(MInt deg = 0; deg < a_pRfnmnt(i) + 2; deg++) {
          a_nodePosition(i, deg) = Fd::nodePosLobattoPoints(a_pRfnmnt(i), deg - 1);
        }
        break;
      }
      case LEGENDRE_INTERP: {
        for(MInt deg = 0; deg < a_pRfnmnt(i) + 2; deg++) {
          if(deg == 0) {
            a_nodePosition(i, deg) = -F1;
          } else if(deg < a_pRfnmnt(i) + 1) {
            a_nodePosition(i, deg) = Fd::nodePosLobattoPoints(a_pRfnmnt(i), deg - 1);
          } else {
            a_nodePosition(i, deg) = F1;
          }
        }
        break;
      }
      default: {
        mTerm(1, AT_, "ERROR: Interpolation method unknown in FC solver!!!");
      }
    }
  }

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt c = 0; c < noHaloCells(i); c++) {
      a_isHalo(haloCell(i, c)) = true;
    }
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt c = 0; c < noWindowCells(i); c++) {
      a_isWindow(windowCell(i, c)) = true;
    }
  }
}

/** \brief Set integration weights and integration nodes and store
 * them in arrays for later.
 *
 * \author Moritz Waldmann
 * \date 20.04.23
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::setIntegrationWeights() {
  TRACE();

  if(m_gaussPoints != nullptr) {
    mDeallocate(m_gaussPoints);
  }
  if(m_gaussPoints != nullptr) {
    mDeallocate(m_gaussWeights);
  }
  mAlloc(m_gaussPoints, m_polyDeg + 1, m_polyDeg + 2, "m_gaussPoints", F0, AT_);
  mAlloc(m_gaussWeights, m_polyDeg + 1, m_polyDeg + 2, "m_gaussWeights", F0, AT_);

  for(MInt p = 0; p < m_polyDeg + 1; p++) {
    Fd::gaussPoint(p + 1, m_gaussPoints[p], m_gaussWeights[p]);
    MInt halfLen = (m_polyDeg + 2 + 1) / 2;
    for(MInt n = 0; n < halfLen; n++) {
      MFloat saveP = m_gaussPoints[p][n];
      MFloat saveW = m_gaussWeights[p][n];
      m_gaussPoints[p][n] = m_gaussPoints[p][m_polyDeg + 1 - n];
      m_gaussPoints[p][m_polyDeg + 1 - n] = saveP;
      m_gaussWeights[p][n] = m_gaussWeights[p][m_polyDeg + 1 - n];
      m_gaussWeights[p][m_polyDeg + 1 - n] = saveW;
    }
  }
}
/** \brief Deactivates of all cells located outside of the geometry
 *
 * \author Moritz Waldmann
 * \date 20.02.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::deactivateCells() {
  // Store level information in cell collector for faster access:
  // halo cells are stored,
  // the poisson ratio is stored in case it differs from cell to cell,
  // the p-Refinement per cell is stored,
  // the number of points per cell is stored,
  // the lobatto points and weights are stored (are these points necessary??)
  for(MInt i = 0; i < a_noCells(); i++) {
    if(!c_isLeafCell(i)) continue;
    if(a_isBndryCell(i)) continue;

    MFloat point[nDim] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      point[d] = c_coordinate(i, d);
    }
    MBool outside = m_geometry->pointIsInside2(point);
    a_isActive(i) = !outside;
  }
}

/** \brief Calculate the inverse of the jacobian
 * matrix and store it in the cell collector. For square/cubic cells
 * the jacobian matrix is a diagonal matrix with J[i][i] = cellLength / 2
 * and J[i][j] = 0 for (i != j). Hence, the inverse of the jacobian
 * matrix is J^(-1)[i][i] = 2 / cellLength and J^(-1)[i][j] = 0 for (i != j).
 *
 * \author Moritz Waldmann
 * \date  20.02.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::initJacobianMatrix() {
  TRACE();

  m_log << "Initializing inverse Jacobian determinant... ";

  for(MInt i = 0; i < a_noCells(); i++) {
    a_invJacobian(i) = F2 / grid().cellLengthAtCell(i);
  }

  m_log << "done" << endl;
}

/** \brief Allocate the displacement vector for all points
 *
 * \author Moritz Waldmann
 * \date  20.02.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::allocateNodeVectors() {
  TRACE();

  createElementToNodeMapping();

  mAlloc(m_globalNodeIdOffsets, noDomains() + 1, "m_globalNodeIdOffsets", 0, AT_);

  // TODO: Should be implemented in a more efficient way.
  if(noDomains() > 1) {
    // In case of local p-Refinement, exchange the number of nodes per cell first
    this->exchangeData(&(a_noNodes(0)), 1);

    for(MInt d = 0; d < noDomains(); d++) {
      if(domainId() == d) {
        createElementToNodeMappingGlobal();
      }

      exchangeNodeIds(d);

      MPI_Allreduce(MPI_IN_PLACE, &m_totalNoNodesGlobal, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "m_totalNoNodesGlobal");
      m_globalNodeIdOffsets[d + 1] = m_totalNoNodesGlobal;
    }

    reorderLocalNodeIds();
  } else {
    createElementToNodeMappingGlobal();

    m_globalNodeIdOffsets[1] = m_totalNoNodesGlobal;
  }

  m_noInternalNodes = m_globalNodeIdOffsets[domainId() + 1] - m_globalNodeIdOffsets[domainId()];

  if(m_testRun) {
    consistencyCheck();

    if(noDomains() > 1) consistencyCheck2();
  }

  mAlloc(m_nodalDisplacements, m_totalNoNodes * nDim, "m_nodalDisplacements", F0, AT_);
  mAlloc(m_totalNodalDisplacements, m_totalNoNodes * nDim, "m_totalNodalDisplacements", F0, AT_);

  mAlloc(m_nodalLoadVector, m_totalNoNodes * nDim, "m_nodalLoadVector", F0, AT_);
  mAlloc(m_externalLoadVector, m_totalNoNodes * nDim, "m_externalLoadVector", F0, AT_);
  mAlloc(m_internalLoadVector, m_totalNoNodes * nDim, "m_internalLoadVector", F0, AT_);
  mAlloc(m_reactionForceVector, m_totalNoNodes * nDim, "m_reactionForceVector", F0, AT_);
  mAlloc(m_localToGlobalId, m_totalNoNodes, "m_localToGlobalId", -1, AT_);
  mAlloc(m_globalToLocalId, m_totalNoNodesGlobal, "m_globalToLocalId", -1, AT_);

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;
    for(MInt node = 0; node < a_noNodes(cell); node++) {
      MInt localNodeId = a_nodeIdsLocal(cell, node);
      MInt globalNodeId = a_nodeIdsGlobal(cell, node);
      if(m_localToGlobalId[localNodeId] == -1) {
        m_localToGlobalId[localNodeId] = globalNodeId;
      } else {
        if(m_localToGlobalId[localNodeId] != globalNodeId) {
          mTerm(1, AT_, "ERROR in local to global Mapping");
        }
      }
      if(m_globalToLocalId[globalNodeId] == -1) {
        m_globalToLocalId[globalNodeId] = localNodeId;
      } else {
        if(m_globalToLocalId[globalNodeId] != localNodeId) {
          mTerm(1, AT_, "ERROR in global to local Mapping");
        }
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::reorderLocalNodeIds() {
  TRACE();

  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = 2 * nDim;

  MInt counter = 0;
  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!c_isLeafCell(cell)) continue;
    if(!a_isActive(cell)) continue;
    for(MInt node = 0; node < a_noNodes(cell); node++) {
      if(a_nodeIdsGlobal(cell, node) >= m_globalNodeIdOffsets[domainId()]
         && a_nodeIdsGlobal(cell, node) < m_globalNodeIdOffsets[domainId() + 1]) {
        a_nodeIdsLocal(cell, node) = a_nodeIdsGlobal(cell, node) - m_globalNodeIdOffsets[domainId()];
      } else {
        a_nodeIdsLocal(cell, node) = -1;
      }
    }
  }

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!c_isLeafCell(cell)) continue;
    if(!a_isActive(cell)) continue;
    for(MInt node = 0; node < a_noNodes(cell); node++) {
      if(a_nodeIdsLocal(cell, node) == -1) {
        a_nodeIdsLocal(cell, node) =
            (m_globalNodeIdOffsets[domainId() + 1] - m_globalNodeIdOffsets[domainId()]) + counter;
        counter++;
        if(node < noVertices) {
          for(MInt dir = 0; dir < noVertices - 1; dir++) {
            if(c_neighborId(cell, Fd::nghbrCellOfVertex(node, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfVertex(node, dir));
            if(a_nodeIdsLocal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) > -1) {
              cerr << a_nodeIdsLocal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) << " OVERWRITTEN in REORDER" << endl;
            } else {
              a_nodeIdsLocal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) = a_nodeIdsLocal(cell, node);
            }
          }
        } else if(node < (noVertices + noEdges * a_pRfnmnt(cell))) {
          MInt pos = node - noVertices;
          MInt edge = 0;
          for(MInt e = 0; e < noEdges; e++) {
            if(e * a_pRfnmnt(cell) < pos) {
              edge++;
            } else {
              edge--;
              break;
            }
          }
          MInt point = pos - (edge * a_pRfnmnt(cell));
          pos = noVertices + (edge * a_pRfnmnt(cell));
          for(MInt dir = 0; dir < 3; dir++) {
            if(c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir));
            MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(cell);
            if(a_nodeIdsLocal(nghbrId, oppositePos + point) > -1) {
              cerr << a_nodeIdsLocal(nghbrId, oppositePos + point) << " OVERWRITTEN in REORDER" << endl;
            } else {
              a_nodeIdsLocal(nghbrId, oppositePos + point) = a_nodeIdsLocal(cell, pos + point);
            }
          }
        } else if(node < (noVertices + noEdges * a_pRfnmnt(cell) + noSurfaces * a_pRfnmnt(cell) * a_pRfnmnt(cell))) {
          MInt pos = node - noVertices + noEdges * a_pRfnmnt(cell);
          MInt surface = 0;
          for(MInt s = 0; s < noSurfaces; s++) {
            if(s * a_pRfnmnt(cell) * a_pRfnmnt(cell) < pos) {
              surface++;
            } else {
              surface--;
              break;
            }
          }
          MInt point = pos - (surface * a_pRfnmnt(cell) * a_pRfnmnt(cell));
          pos = noVertices + (noEdges * a_pRfnmnt(cell)) + (surface * a_pRfnmnt(cell) * a_pRfnmnt(cell));
          if(c_neighborId(cell, Fd::nghbrCellOfSurface(surface)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfSurface(surface));
          MInt oppositePos = noVertices + noEdges * a_pRfnmnt(cell)
                             + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
          if(a_nodeIdsLocal(nghbrId, oppositePos + point) > -1) {
            cerr << a_nodeIdsLocal(nghbrId, oppositePos + point) << " OVERWRITTEN in REORDER" << endl;
          } else {
            a_nodeIdsLocal(nghbrId, oppositePos + point) = a_nodeIdsLocal(cell, pos + point);
          }
        }
      }
    }
  }
}

/** \brief Send the global nodeIds from one domain to all neighbor domains
 *
 * \author Moritz Waldmann
 * \date  20.03.22
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::exchangeNodeIds(MInt currentDomain) {
  TRACE();

  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = 2 * nDim;

  // Create arrays for the data offsets and data count
  MIntScratchSpace windowDataPerRank(noNeighborDomains(), "windowDataPerRank", AT_);
  MIntScratchSpace haloDataPerRank(noNeighborDomains(), "haloDataPerRank", AT_);
  MIntScratchSpace windowDataOffsets(noNeighborDomains() + 1, "windowDataOffsets", AT_);
  MIntScratchSpace haloDataOffsets(noNeighborDomains() + 1, "haloDataOffsets", AT_);
  windowDataPerRank.fill(0);
  haloDataPerRank.fill(0);
  windowDataOffsets.fill(0);
  haloDataOffsets.fill(0);

  // Determine the data count and the data offsets
  MInt totalCounterWindow = 0;
  MInt totalCounterHalo = 0;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt dataCounterWindow = 0;
    MInt dataCounterHalo = 0;
    for(MInt w = 0; w < noWindowCells(n); w++) {
      dataCounterWindow += a_noNodes(windowCell(n, w));
    }
    for(MInt h = 0; h < noHaloCells(n); h++) {
      dataCounterHalo += a_noNodes(haloCell(n, h));
    }
    windowDataPerRank(n) = dataCounterWindow;
    haloDataPerRank(n) = dataCounterHalo;
    totalCounterWindow += dataCounterWindow;
    totalCounterHalo += dataCounterHalo;
    windowDataOffsets(n + 1) = totalCounterWindow;
    haloDataOffsets(n + 1) = totalCounterHalo;
  }

  // Allocate the buffers
  MInt* sendBuffers{};
  MInt* receiveBuffers{};
  MPI_Request* mpi_request{};
  mAlloc(sendBuffers, totalCounterWindow, "sendBuffers", 0, AT_);
  mAlloc(receiveBuffers, totalCounterHalo, "receiveBuffers", 0, AT_);
  mAlloc(mpi_request, noNeighborDomains(), "mpi_request", AT_);

  // Gather data
  if(domainId() == currentDomain) {
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MInt cnt = 0;
      for(MInt j = 0; j < noWindowCells(n); j++) {
        for(MInt node = 0; node < a_noNodes(windowCell(n, j)); node++) {
          sendBuffers[windowDataOffsets(n) + cnt] = a_nodeIdsGlobal(windowCell(n, j), node);
          cnt++;
        }
      }
    }

    // Send data
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MPI_Issend(&sendBuffers[windowDataOffsets(n)], windowDataPerRank(n), MPI_INT, neighborDomain(n), 0, mpiComm(),
                 &mpi_request[n], AT_, "sendBuffers[n]");
    }
  }

  // Receive data
  MPI_Status status;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    if(currentDomain != neighborDomain(n)) continue;
    MPI_Recv(&receiveBuffers[haloDataOffsets(n)], haloDataPerRank(n), MPI_INT, neighborDomain(n), 0, mpiComm(), &status,
             AT_, "receiveBuffers[n]");
  }
  if(domainId() == currentDomain) {
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MPI_Wait(&mpi_request[n], &status, AT_);
    }
  }

  // Scatter data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    if(currentDomain != neighborDomain(n)) continue;
    MInt cnt = 0;
    for(MInt j = 0; j < noHaloCells(n); j++) {
      for(MInt node = 0; node < a_noNodes(haloCell(n, j)); node++) {
        cnt++;
        if(receiveBuffers[haloDataOffsets(n) + cnt - 1] == -1) continue;
        if((a_nodeIdsGlobal(haloCell(n, j), node) > -1)) continue;
        a_nodeIdsGlobal(haloCell(n, j), node) = receiveBuffers[haloDataOffsets(n) + cnt - 1];
        // If the global nodeId is set in the halo, it should be set also in the neighboring cells sharing the node with
        // the halo cell
        if(node < noVertices) {
          for(MInt dir = 0; dir < noVertices - 1; dir++) {
            if(c_neighborId(haloCell(n, j), Fd::nghbrCellOfVertex(node, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(haloCell(n, j), Fd::nghbrCellOfVertex(node, dir));
            if(a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) > -1) {
              cerr << a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) << " OVERWRITTEN in EXCHANGE" << endl;
            } else {
              a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(node, dir)) = a_nodeIdsGlobal(haloCell(n, j), node);
            }
          }
        } else if(node < (noVertices + noEdges * a_pRfnmnt(haloCell(n, j)))) {
          MInt pos = node - noVertices;
          MInt edge = 0;
          for(MInt e = 0; e < noEdges; e++) {
            if(e * a_pRfnmnt(haloCell(n, j)) < pos) {
              edge++;
            } else {
              edge--;
              break;
            }
          }
          MInt point = pos - (edge * a_pRfnmnt(haloCell(n, j)));
          pos = noVertices + (edge * a_pRfnmnt(haloCell(n, j)));
          for(MInt dir = 0; dir < 3; dir++) {
            if(c_neighborId(haloCell(n, j), Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(haloCell(n, j), Fd::nghbrCellOfEdge(edge, dir));
            MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(haloCell(n, j));
            if(a_nodeIdsGlobal(nghbrId, oppositePos + point) > -1) {
              cerr << a_nodeIdsGlobal(nghbrId, oppositePos + point) << " OVERWRITTEN in EXCHANGE" << endl;
            } else {
              a_nodeIdsGlobal(nghbrId, oppositePos + point) = a_nodeIdsGlobal(haloCell(n, j), pos + point);
            }
          }
        } else if(node < (noVertices + noEdges * a_pRfnmnt(haloCell(n, j))
                          + noSurfaces * a_pRfnmnt(haloCell(n, j)) * a_pRfnmnt(haloCell(n, j)))) {
          MInt pos = node - noVertices + noEdges * a_pRfnmnt(haloCell(n, j));
          MInt surface = 0;
          for(MInt s = 0; s < noSurfaces; s++) {
            if(s * a_pRfnmnt(haloCell(n, j)) * a_pRfnmnt(haloCell(n, j)) < pos) {
              surface++;
            } else {
              surface--;
              break;
            }
          }
          MInt point = pos - (surface * a_pRfnmnt(haloCell(n, j)) * a_pRfnmnt(haloCell(n, j)));
          pos = noVertices + (noEdges * a_pRfnmnt(haloCell(n, j)))
                + (surface * a_pRfnmnt(haloCell(n, j)) * a_pRfnmnt(haloCell(n, j)));
          if(c_neighborId(haloCell(n, j), Fd::nghbrCellOfSurface(surface)) < 0) continue;
          MInt nghbrId = c_neighborId(haloCell(n, j), Fd::nghbrCellOfSurface(surface));
          MInt oppositePos = noVertices + noEdges * a_pRfnmnt(haloCell(n, j))
                             + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(haloCell(n, j)) * a_pRfnmnt(haloCell(n, j));
          if(a_nodeIdsGlobal(nghbrId, oppositePos + point) > -1) {
            cerr << a_nodeIdsGlobal(nghbrId, oppositePos + point) << " OVERWRITTEN in EXCHANGE" << endl;
          } else {
            a_nodeIdsGlobal(nghbrId, oppositePos + point) = a_nodeIdsGlobal(haloCell(n, j), pos + point);
          }
        }
      }
    }
  }

  mDeallocate(sendBuffers);
  mDeallocate(receiveBuffers);
  mDeallocate(mpi_request);
}

/** \brief Consistency check for the global nodeIds, checking whether the global nodeId
 * of one cell and the corresponding neighbor cells sharing the same node is identical.
 * If it is not, something went wrong!!!
 *
 * \author Moritz Waldmann
 * \date  10.05.22
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::consistencyCheck() {
  TRACE();

  // Set nodes:
  // Firstly, set the nodes vertex wise
  // Secondly, set the nodes edge wise
  // Thirdly, set the nodes surface wise
  // Finally, set the inner nodes which have no neighbor
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = 2 * nDim;

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!c_isLeafCell(cell)) continue;
    if(!a_isActive(cell)) continue;
    // Check the node ids for the vertex nodes
    for(MInt vertex = 0; vertex < noVertices; vertex++) {
      if(a_nodeIdsGlobal(cell, vertex) == -1) {
        MInt globalId = -1;
        MInt noDiffCandidates = 0;
        for(MInt dir = 0; dir < noVertices - 1; dir++) {
          if(c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir));
          if(a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(vertex, dir)) > -1) {
            if(globalId != a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(vertex, dir))) {
              noDiffCandidates++;
            }
          }
        }
        if(noDiffCandidates == 1) {
          a_nodeIdsGlobal(cell, vertex) = globalId;
        }
      } else {
        MInt id = a_nodeIdsGlobal(cell, vertex);
        for(MInt dir = 0; dir < noVertices - 1; dir++) {
          if(c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir));
          if(a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(vertex, dir)) != id) {
            stringstream errorMessage;
            errorMessage << "INCONSISTENCY IN GLOBAL NODE ID FOR GLOBAL CELL ID " << c_globalId(cell) << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        }
      }
    }
    // Check the node ids for the edge nodes
    for(MInt edge = 0; edge < noEdges; edge++) {
      MInt pos = noVertices + edge * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell); points++) {
        if(a_nodeIdsGlobal(cell, pos + points) == -1) {
          MInt globalId = -1;
          MInt noDiffCandidates = 0;
          for(MInt dir = 0; dir < 3; dir++) {
            if(c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir));
            MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(cell);
            if(a_nodeIdsGlobal(nghbrId, oppositePos + points) > -1) {
              if(globalId != a_nodeIdsGlobal(nghbrId, oppositePos + points)) {
                noDiffCandidates++;
              }
            }
          }
          if(noDiffCandidates == 1) {
            a_nodeIdsGlobal(cell, pos + points) = globalId;
          }
        } else {
          MInt id = a_nodeIdsGlobal(cell, pos + points);
          for(MInt dir = 0; dir < 3; dir++) {
            if(c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
            MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir));
            MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(cell);
            if(a_nodeIdsGlobal(nghbrId, oppositePos + points) != id) {
              stringstream errorMessage;
              errorMessage << "INCONSISTENCY IN GLOBAL NODE ID FOR GLOBAL CELL ID " << c_globalId(cell) << endl;
              mTerm(1, AT_, errorMessage.str());
            }
          }
        }
      }
    }
    // Check the node ids for the surfaces
    for(MInt surface = 0; surface < noSurfaces; surface++) {
      MInt pos = noVertices + noEdges * a_pRfnmnt(cell) + surface * a_pRfnmnt(cell) * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell) * a_pRfnmnt(cell); points++) {
        if(a_nodeIdsGlobal(cell, pos + points) == -1) {
          MInt globalId = -1;
          if(c_neighborId(cell, Fd::nghbrCellOfSurface(surface)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfSurface(surface));
          MInt oppositePos = noVertices + noEdges * a_pRfnmnt(cell)
                             + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
          if(globalId != a_nodeIdsGlobal(nghbrId, oppositePos + points)) {
            a_nodeIdsGlobal(cell, pos + points) = globalId;
          }
        } else {
          MInt id = a_nodeIdsGlobal(cell, pos + points);
          if(c_neighborId(cell, Fd::nghbrCellOfSurface(surface)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfSurface(surface));
          MInt oppositePos = noVertices + noEdges * a_pRfnmnt(cell)
                             + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
          if(a_nodeIdsGlobal(nghbrId, oppositePos + points) != id) {
            stringstream errorMessage;
            errorMessage << "INCONSISTENCY IN GLOBAL NODE ID FOR GLOBAL CELL ID " << c_globalId(cell) << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        }
      }
    }
  }
}

/** \brief Consistency check for the global nodeIds, checking whether the global nodeId
 * of a halo cell and the corresponding window cell on the neighbor domaine is identical.
 * If it is not, something went wrong!!!
 *
 * \author Moritz Waldmann
 * \date  10.05.22
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::consistencyCheck2() {
  TRACE();

  // Create arrays for the data offsets and data count
  MIntScratchSpace windowDataPerRank(noNeighborDomains(), "windowDataPerRank", AT_);
  MIntScratchSpace haloDataPerRank(noNeighborDomains(), "haloDataPerRank", AT_);
  MIntScratchSpace windowDataOffsets(noNeighborDomains() + 1, "windowDataOffsets", AT_);
  MIntScratchSpace haloDataOffsets(noNeighborDomains() + 1, "haloDataOffsets", AT_);
  windowDataPerRank.fill(0);
  haloDataPerRank.fill(0);
  windowDataOffsets.fill(0);
  haloDataOffsets.fill(0);

  // Determine the data count and the data offsets
  MInt totalCounterWindow = 0;
  MInt totalCounterHalo = 0;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt dataCounterWindow = 0;
    MInt dataCounterHalo = 0;
    for(MInt w = 0; w < noWindowCells(n); w++) {
      dataCounterWindow += a_noNodes(windowCell(n, w));
    }
    for(MInt h = 0; h < noHaloCells(n); h++) {
      dataCounterHalo += a_noNodes(haloCell(n, h));
    }
    windowDataPerRank(n) = dataCounterWindow;
    haloDataPerRank(n) = dataCounterHalo;
    totalCounterWindow += dataCounterWindow;
    totalCounterHalo += dataCounterHalo;
    windowDataOffsets(n + 1) = totalCounterWindow;
    haloDataOffsets(n + 1) = totalCounterHalo;
  }

  // Allocate the buffers
  MInt* sendBuffers{};
  MInt* receiveBuffers{};
  MPI_Request* mpi_request{};
  mAlloc(sendBuffers, totalCounterWindow, "sendBuffers", 0, AT_);
  mAlloc(receiveBuffers, totalCounterHalo, "receiveBuffers", 0, AT_);
  mAlloc(mpi_request, noNeighborDomains(), "mpi_request", AT_);

  // Gather data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt cnt = 0;
    for(MInt j = 0; j < noWindowCells(n); j++) {
      for(MInt node = 0; node < a_noNodes(windowCell(n, j)); node++) {
        sendBuffers[windowDataOffsets(n) + cnt] = a_nodeIdsGlobal(windowCell(n, j), node);
        cnt++;
      }
    }
  }

  // Send data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Issend(&sendBuffers[windowDataOffsets(n)], windowDataPerRank(n), MPI_INT, neighborDomain(n), 0, mpiComm(),
               &mpi_request[n], AT_, "sendBuffers[n]");
  }

  // Receive data
  MPI_Status status;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Recv(&receiveBuffers[haloDataOffsets(n)], haloDataPerRank(n), MPI_INT, neighborDomain(n), 0, mpiComm(), &status,
             AT_, "receiveBuffers[n]");
  }
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Wait(&mpi_request[n], &status, AT_);
  }

  // Scatter data and compare it to the data of the halo cell
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt cnt = 0;
    for(MInt j = 0; j < noHaloCells(n); j++) {
      for(MInt node = 0; node < a_noNodes(haloCell(n, j)); node++) {
        cnt++;
        if(a_nodeIdsGlobal(haloCell(n, j), node) != receiveBuffers[haloDataOffsets(n) + cnt - 1]) {
          stringstream errorMessage;
          errorMessage << "INCONSISTENCY IN GLOBAL NODE ID FOR GLOBAL CELL ID " << c_globalId(haloCell(n, j)) << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    }
  }

  mDeallocate(sendBuffers);
  mDeallocate(receiveBuffers);
  mDeallocate(mpi_request);
}


/** \brief Calculates the mapping from the element
 * to the node. The nodeId is global, i.e., it is unique
 * on the whole computational domain.
 *
 * \author Moritz Waldmann
 * \date 20.04.22
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::createElementToNodeMappingGlobal() {
  TRACE();

  // Set nodes:
  // Firstly, set the nodes vertex wise
  // Secondly, set the nodes edge wise
  // Thirdly, set the nodes surface wise
  // Finally, set the inner nodes which have no neighbor
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = 2 * nDim;

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(a_isHalo(cell)) continue;
    if(!c_isLeafCell(cell)) continue;
    if(!a_isActive(cell)) continue;
    // Set the node ids for the vertex nodes
    for(MInt vertex = 0; vertex < noVertices; vertex++) {
      if(a_nodeIdsGlobal(cell, vertex) == -1) {
        a_nodeIdsGlobal(cell, vertex) = m_totalNoNodesGlobal;
        m_totalNoNodesGlobal++;
        // Check the neighbors and set the node ids there
        for(MInt dir = 0; dir < noVertices - 1; dir++) {
          if(c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir));
          a_nodeIdsGlobal(nghbrId, Fd::vertexIdOfOppCell(vertex, dir)) = a_nodeIdsGlobal(cell, vertex);
        }
      }
    }
    // Set the node ids for the edge nodes
    for(MInt edge = 0; edge < noEdges; edge++) {
      MInt pos = noVertices + edge * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell); points++) {
        if(a_nodeIdsGlobal(cell, pos + points) == -1) {
          a_nodeIdsGlobal(cell, pos + points) = m_totalNoNodesGlobal;
          m_totalNoNodesGlobal++;
        }
      }
      // Check the neighbors and set the node ids there
      for(MInt dir = 0; dir < 3; dir++) {
        if(c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
        MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir));
        MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(cell);
        for(MInt points = 0; points < a_pRfnmnt(cell); points++) {
          a_nodeIdsGlobal(nghbrId, oppositePos + points) = a_nodeIdsGlobal(cell, pos + points);
        }
      }
    }
    // Set the node ids for the surfaces
    for(MInt surface = 0; surface < noSurfaces; surface++) {
      MInt pos = noVertices + noEdges * a_pRfnmnt(cell) + surface * a_pRfnmnt(cell) * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell) * a_pRfnmnt(cell); points++) {
        if(a_nodeIdsGlobal(cell, pos + points) == -1) {
          a_nodeIdsGlobal(cell, pos + points) = m_totalNoNodesGlobal;
          m_totalNoNodesGlobal++;
        }
      }
      // Check the neighbors and set the node ids there
      if(c_neighborId(cell, Fd::nghbrCellOfSurface(surface)) < 0) continue;
      MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfSurface(surface));
      MInt oppositePos =
          noVertices + noEdges * a_pRfnmnt(cell) + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell) * a_pRfnmnt(cell); points++) {
        a_nodeIdsGlobal(nghbrId, oppositePos + points) = a_nodeIdsGlobal(cell, pos + points);
      }
    }
    // Set the node ids for all inner points
    MInt noInnerNodes = a_pRfnmnt(cell) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
    for(MInt innerNodes = 0; innerNodes < noInnerNodes; innerNodes++) {
      MInt pos = noVertices + noEdges * a_pRfnmnt(cell) + noSurfaces * a_pRfnmnt(cell) * a_pRfnmnt(cell) + innerNodes;
      if(a_nodeIdsGlobal(cell, pos) == -1) {
        a_nodeIdsGlobal(cell, pos) = m_totalNoNodesGlobal;
        m_totalNoNodesGlobal++;
      }
    }
  }
}

/** \brief Calculates the mapping from the element
 * to the node. The nodeId is local, i.e., it is unique
 * on the domain of corresponding rank. A mapping from
 * the nodes to the elements is not necessary.
 *
 * \author Moritz Waldmann
 * \date 20.02.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::createElementToNodeMapping() {
  TRACE();

  m_totalNoNodes = 0;

  // reset the node ids stored for each cell
  // TODO: At the moment the allocation is done
  // statically, but it should be dynamic to safe memory
  for(MInt cell = 0; cell < a_noCells(); cell++) {
    for(MInt node = 0; node < a_noNodes(cell); node++) {
      a_nodeIdsLocal(cell, node) = -1;
      a_nodeIdsGlobal(cell, node) = -1;
    }
  }

  // Set nodes:
  // Firstly, set the nodes vertex wise
  // Secondly, set the nodes edge wise
  // Thirdly, set the nodes surface wise
  // Finally, set the inner nodes which have no neighbor
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = 2 * nDim;
  MInt innerCounter = 0;
  MInt surfaceCounter = 0;
  MInt edgeCounter = 0;
  MInt vertexCounter = 0;

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!c_isLeafCell(cell)) continue;
    if(!a_isActive(cell)) continue;
    // Set the node ids for the vertex nodes
    for(MInt vertex = 0; vertex < noVertices; vertex++) {
      if(a_nodeIdsLocal(cell, vertex) == -1) {
        a_nodeIdsLocal(cell, vertex) = m_totalNoNodes;
        m_totalNoNodes++;
        vertexCounter++;
        // Check the neighbors and set the node ids there
        for(MInt dir = 0; dir < noVertices - 1; dir++) {
          if(c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir)) < 0) continue;
          MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfVertex(vertex, dir));
          a_nodeIdsLocal(nghbrId, Fd::vertexIdOfOppCell(vertex, dir)) = a_nodeIdsLocal(cell, vertex);
        }
      }
    }
    // Set the node ids for the edge nodes
    for(MInt edge = 0; edge < noEdges; edge++) {
      MInt pos = noVertices + edge * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell); points++) {
        if(a_nodeIdsLocal(cell, pos + points) == -1) {
          a_nodeIdsLocal(cell, pos + points) = m_totalNoNodes;
          m_totalNoNodes++;
          edgeCounter++;
        }
      }
      // Check the neighbors and set the node ids there
      for(MInt dir = 0; dir < 3; dir++) {
        if(c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir)) < 0) continue;
        MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfEdge(edge, dir));
        MInt oppositePos = noVertices + Fd::edgeIdOfOppCell(edge, dir) * a_pRfnmnt(cell);
        for(MInt points = 0; points < a_pRfnmnt(cell); points++) {
          a_nodeIdsLocal(nghbrId, oppositePos + points) = a_nodeIdsLocal(cell, pos + points);
        }
      }
    }
    // Set the node ids for the surfaces
    for(MInt surface = 0; surface < noSurfaces; surface++) {
      MInt pos = noVertices + noEdges * a_pRfnmnt(cell) + surface * a_pRfnmnt(cell) * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell) * a_pRfnmnt(cell); points++) {
        if(a_nodeIdsLocal(cell, pos + points) == -1) {
          a_nodeIdsLocal(cell, pos + points) = m_totalNoNodes;
          m_totalNoNodes++;
          surfaceCounter++;
        }
      }
      // Check the neighbors and set the node ids there
      if(c_neighborId(cell, Fd::nghbrCellOfSurface(surface)) < 0) continue;
      MInt nghbrId = c_neighborId(cell, Fd::nghbrCellOfSurface(surface));
      MInt oppositePos =
          noVertices + noEdges * a_pRfnmnt(cell) + Fd::surfaceIdOfOppCell(surface) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
      for(MInt points = 0; points < a_pRfnmnt(cell) * a_pRfnmnt(cell); points++) {
        a_nodeIdsLocal(nghbrId, oppositePos + points) = a_nodeIdsLocal(cell, pos + points);
      }
    }
    // Set the node ids for all inner points
    MInt noInnerNodes = a_pRfnmnt(cell) * a_pRfnmnt(cell) * a_pRfnmnt(cell);
    for(MInt innerNodes = 0; innerNodes < noInnerNodes; innerNodes++) {
      MInt pos = noVertices + noEdges * a_pRfnmnt(cell) + noSurfaces * a_pRfnmnt(cell) * a_pRfnmnt(cell) + innerNodes;
      if(a_nodeIdsLocal(cell, pos) == -1) {
        a_nodeIdsLocal(cell, pos) = m_totalNoNodes;
        m_totalNoNodes++;
        innerCounter++;
      }
    }
  }

  cout << "Number of Nodes on cell vertices (no vertices per cell: " << noVertices << "): " << vertexCounter << endl;
  cout << "Number of Nodes on cell edges (no edges per cell: " << noEdges << ")     : " << edgeCounter << endl;
  cout << "Number of Nodes on cell surfaces (no surfaces per cell: " << noSurfaces << "): " << surfaceCounter << endl;
  cout << "Number of Nodes within the cells                          : " << innerCounter << endl;
  cout << "Number of Nodes total                                     : " << m_totalNoNodes << endl;
  cout << "Number of Cells total without halo cells (on domain " << domainId() << ")    : " << grid().noInternalCells()
       << endl;
  cout << "Number of Cells total with halo cells (on domain " << domainId() << ")       : " << a_noCells() << endl;
  cout << endl;
  m_log << "Number of Nodes on cell vertices (no vertices per cell: " << noVertices << "): " << vertexCounter << endl;
  m_log << "Number of Nodes on cell edges (no edges per cell: " << noEdges << ")      : " << edgeCounter << endl;
  m_log << "Number of Nodes on cell surfaces (no surfaces per cell: " << noSurfaces << "): " << surfaceCounter << endl;
  m_log << "Number of Nodes within the cells                          : " << innerCounter << endl;
  m_log << "Number of Nodes total                                     : " << m_totalNoNodes << endl;
  m_log << endl;
}


/** \brief Transformes the coordinate z from the reference element
 * to the global coordinate system.
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::interpolateGeometry(const MInt pCellId, const MFloat* z, MFloat* x) {
  TRACE();

  const MFloat cellHalfLength = c_cellLengthAtCell(pCellId) * F1B2;

  for(MInt d = 0; d < nDim; d++) {
    x[d] = c_coordinate(pCellId, d) + z[d] * cellHalfLength;
  }
}

/** \brief Transformes the coordinate z from the sub cell of a reference element
 * to the global coordinate system.
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::interpolateSubCellGeometry(const MInt pCellId, const MFloat* z, MFloat* x, const MInt subCellLevel,
                                                const MFloat* subCellCoord) {
  TRACE();

  const MFloat cellHalfLength = c_cellLengthAtCell(pCellId) * FFPOW2(subCellLevel + 1);

  for(MInt d = 0; d < nDim; d++) {
    x[d] = subCellCoord[d] + z[d] * cellHalfLength;
  }
}

/** \brief Transformes the coordinate x from the global coordinate system
 * to the local coordinate system.
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::transformToLocal(const MInt pCellId, const MFloat* x, MFloat* z) {
  TRACE();

  const MFloat F1BcellHalfLength = F1 / (c_cellLengthAtCell(pCellId) * F1B2);

  for(MInt d = 0; d < nDim; d++) {
    z[d] = (x[d] - c_coordinate(pCellId, d)) * F1BcellHalfLength;
  }
}


/** \brief Calculates the displacement of a cell by interpolation
 * That is, the product of the displacement interpolation matrix
 * and the displacement is summed up for each node.
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::calcElementDisplacements() {
  TRACE();

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    const MInt pCellId = cell;

    if(!a_isActive(pCellId)) {
      for(MInt d = 0; d < nDim; d++) {
        a_elementDispl(pCellId, d) = F0;
      }
      continue;
    }

    std::array<MFloat, nDim> z{};

    getDisplacementsAtPoint(pCellId, z.data(), m_totalNodalDisplacements, &a_elementDispl(pCellId, 0));
  }
}

template <MInt nDim>
void FcSolver<nDim>::getDisplacementsAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* result) {
  TRACE();

  // Calculate the lagrangian polynome for calculating the
  // element displacement by interpolation of the nodal displacement
  MFloatScratchSpace L_coef(nDim, a_noNodes(cellId) * nDim, AT_, "L_coef");

  L_coef.fill(F1);
  getDisplacementInterpolationMatrix(cellId, z, L_coef);

  if(m_testRun && m_polyDeg < 1) getDisplacementInterpolationMatrixDebug(cellId, z, L_coef);

  // Loop over all directions and all points to get the
  // displacement of the cell
  for(MInt dim = 0; dim < nDim; dim++) {
    MFloat displacement = F0;
    for(MInt node = 0; node < a_noNodes(cellId); node++) {
      for(MInt d = 0; d < nDim; d++) {
        MInt nodeId = a_nodeIdsLocal(cellId, node);
        displacement += L_coef(dim, node * nDim + d) * displ[nodeId * nDim + d];
      }
    }
    result[dim] = displacement;
  }
}

/** \brief Calculate the strain of a cell.
 * The matrix product of strain interpolation matrix and
 * displacement vector gives the average strain per cell.
 *
 * \author Moritz Waldmann
 * \date 15.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::calcElementStrains() {
  TRACE();

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;

    // Calculate the strain interpolation matrix
    MFloat z[nDim] = {F0};

    getStrainsAtPoint(cell, z, m_totalNodalDisplacements, &a_elementStrains(cell, 0));
  }
}

template <MInt nDim>
void FcSolver<nDim>::getStrainsAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* strains) {
  TRACE();

  MFloatScratchSpace Be(m_noStrains, nDim * a_noNodes(cellId), AT_, "Be");
  getStrainInterpolationMatrix(Be, cellId, z);

  // Calculate the strains from the matrix product Be * u
  // strain interpolation matrix times displacement
  for(MInt d1 = 0; d1 < m_noStrains; d1++) {
    MFloat strain = F0;
    for(MInt node = 0; node < a_noNodes(cellId); node++) {
      MInt nodeId = a_nodeIdsLocal(cellId, node);
      for(MInt d2 = 0; d2 < nDim; d2++) {
        strain += Be(d1, node * nDim + d2) * displ[nodeId * nDim + d2];
      }
    }
    strains[d1] = strain;
  }
}

/** \brief Calculate the stresses of a cell.
 * The stress tensor is defined as the matrix product of
 * the material tensor and the strains matrix.
 *
 * \author Moritz Waldmann
 * \date 15.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::calcElementStresses() {
  TRACE();

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;

    std::array<MFloat, nDim> z{};
    getStressesAtPoint(cell, z.data(), m_totalNodalDisplacements, &a_elementStresses(cell, 0));
  }
}

template <MInt nDim>
void FcSolver<nDim>::getStressesAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* stresses) {
  TRACE();

  ASSERT(string2enum(solverMethod()) == MAIA_LINEAR_ELASTIC, "This function is not valid for plastic simulations.");

  std::array<MFloat, m_noStrains> strains{};
  getStrainsAtPoint(cellId, z, displ, strains.data());

  MFloatScratchSpace C(m_noStrains, m_noStrains, AT_, "C");
  getMaterialTensor(C, cellId);
  for(MInt d1 = 0; d1 < m_noStrains; d1++) {
    MFloat stress = F0;
    for(MInt d2 = 0; d2 < m_noStrains; d2++) {
      const MFloat strain = strains[d2];
      stress += C(d2, d1) * strain;
    }
    stresses[d1] = stress;
  }
}

/** \brief Calculate the displacement interpolation matrix
 * by calculating the lagrangian polynom at each node.
 *
 * \author Moritz Waldmann
 * \date 01.03.21
 *
 * The algorithm works as follows:
 *
 *  - calculate the lagrange polynom for all possible points and for all
 *  cartesian directions in vector L
 *  - return the lagrange polynom at each node in each direction
 **/
template <MInt nDim>
void FcSolver<nDim>::getDisplacementInterpolationMatrix(const MInt pCellId, const MFloat* z,
                                                        MFloatScratchSpace& L_coef) {
  TRACE();

  const MInt p = a_pRfnmnt(pCellId);
  MFloatScratchSpace L(nDim, p + 2, AT_, "L");

  // initialize the arrays for the mapping
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = (nDim == 2) ? 4 : 6;

  // initialize the L matrix with 1.0
  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 0; i < p + 2; i++) {
      L(d, i) = F1;
    }
  }

  // calculate the lagrangian polynominal
  // see Lagrange polynominal
  // Maybe L can be stored once at the beginning
  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 0; i < p + 2; i++) {
      for(MInt j = 0; j < p + 2; j++) {
        if(i == j) continue;
        L(d, i) *= (z[d] - a_nodePosition(pCellId, j)) / (a_nodePosition(pCellId, i) - a_nodePosition(pCellId, j));
      }
    }
  }

  // Now that all polynominals are calculated, the polynominals
  // for the different directions need to be multiplied for each node
  // First the values at the vertices, then at the edges, than at the
  // surfaces and finally at the inner points are determined
  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt node = 0; node < a_noNodes(pCellId); node++) {
      // vertices
      if(node < noVertices) {
        for(MInt d = 0; d < nDim; d++) {
          if(d == dim) {
            L_coef(dim, node * nDim + d) = F1;
          } else {
            L_coef(dim, node * nDim + d) = F0;
          }
          for(MInt d2 = 0; d2 < nDim; d2++) {
            MInt pos = -1;
            if(Fd::vertexPosition(node, d2) < F0) {
              pos = 0;
            } else {
              pos = p + 1;
            }
            L_coef(dim, node * nDim + d) *= L(d2, pos);
          }
        }
        // edges
      } else if(node < noVertices + noEdges * p) {
        // Which edge do we consider right now?
        MInt edge = (MInt)floor((node - noVertices) / p);
        // Which point of the edge do we consider right now?
        MInt l = (node - noVertices) % p;
        for(MInt d = 0; d < nDim; d++) {
          if(d == dim) {
            L_coef(dim, node * nDim + d) = F1;
          } else {
            L_coef(dim, node * nDim + d) = F0;
          }
          for(MInt d2 = 0; d2 < nDim; d2++) {
            MInt pos = -1;
            if(Fd::edgePosition(edge, d2) < F0) {
              pos = 0;
            } else if(Fd::edgePosition(edge, d2) > F0) {
              pos = p + 1;
            } else {
              pos = l + 1;
            }
            L_coef(dim, node * nDim + d) *= L(d2, pos);
          }
        }
        // surfaces
      } else if(node < noVertices + noEdges * p + noSurfaces * p * p) {
        // Which surface do we consider right now?
        MInt surface = (MInt)floor((node - noVertices - noEdges * p) / (p * p));
        // Which point of the surface do we consider right now?
        MInt l[2] = {0};
        l[1] = (node - noVertices - noEdges * p) % (p * p);
        while(l[1] >= p) {
          l[0]++;
          l[1] -= p;
        }
        for(MInt d = 0; d < nDim; d++) {
          if(d == dim) {
            L_coef(dim, node * nDim + d) = F1;
          } else {
            L_coef(dim, node * nDim + d) = F0;
          }
          MInt counter = 0;
          for(MInt d2 = 0; d2 < nDim; d2++) {
            MInt pos = -1;
            if(Fd::surfacePosition(surface, d2) < F0) {
              pos = 0;
            } else if(Fd::surfacePosition(surface, d2) > F0) {
              pos = p + 1;
            } else {
              pos = l[counter] + 1;
              counter++;
            }
            L_coef(dim, node * nDim + d) *= L(d2, pos);
          }
        }
        // inner points
      } else {
        MInt l[3] = {0};
        // Which point of the inner points do we consider right now?
        l[2] = node - noVertices - noEdges * p - noSurfaces * p * p;
        while(l[2] >= (p * p)) {
          l[0]++;
          l[2] -= (p * p);
        }
        while(l[2] >= p) {
          l[1]++;
          l[2] -= p;
        }
        for(MInt d = 0; d < nDim; d++) {
          if(d == dim) {
            L_coef(dim, node * nDim + d) = F1;
          } else {
            L_coef(dim, node * nDim + d) = F0;
          }
          MInt counter = 0;
          for(MInt d2 = 0; d2 < nDim; d2++) {
            MInt pos = l[counter] + 1;
            counter++;
            L_coef(dim, node * nDim + d) *= L(d2, pos);
          }
        }
      }
    }
  }
}

/** \brief Calculate the derivative of the lagrangian polynom at each node
 *
 * \author Moritz Waldmann
 * \date 01.03.21
 *
 * The algorithm works as follows:
 *
 *  - calculate the lagrange polynom for all possible points and for all
 *  cartesian directions
 *  - calculate the lagrangian coefficient by multiplying the polynoms for
 *  all cartesian directions
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::getDerivativeOfDisplacementInterpolationMatrix(const MInt pCellId, const MFloat* z,
                                                                    MFloatScratchSpace& L_prime) {
  TRACE();

  const MInt p = a_pRfnmnt(pCellId);
  MFloatScratchSpace L(nDim, p + 2, AT_, "L");
  MFloatScratchSpace L_deriv(nDim, p + 2, AT_, "L");

  // initialize the arrays for the mapping
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = (nDim == 2) ? 4 : 6;

  // initialize the L matrix with 1.0 and L_deriv with 0.0
  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 0; i < p + 2; i++) {
      L(d, i) = F1;
      L_deriv(d, i) = F0;
    }
  }

  // calculate the lagrangian polynominal and its derivative
  // see Lagrange polynominal
  // Maybe L can be stored once at the beginning
  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 0; i < p + 2; i++) {
      for(MInt j = 0; j < p + 2; j++) {
        if(i == j) continue;
        L(d, i) *= (z[d] - a_nodePosition(pCellId, j)) / (a_nodePosition(pCellId, i) - a_nodePosition(pCellId, j));
        MFloat L_p = (F1 / (a_nodePosition(pCellId, i) - a_nodePosition(pCellId, j)));
        for(MInt l = 0; l < p + 2; l++) {
          if(l == i || l == j) continue;
          L_p *= (z[d] - a_nodePosition(pCellId, l)) / (a_nodePosition(pCellId, i) - a_nodePosition(pCellId, l));
        }
        L_deriv(d, i) += L_p;
      }
    }
  }

  // Now that all polynominals are calculated, the polynominals
  // for the different directions need to be multiplied for each node
  // First the values at the vertices, then at the edges, than at the
  // surfaces and finally at the inner points are determined
  for(MInt node = 0; node < a_noNodes(pCellId); node++) {
    // vertices
    if(node < noVertices) {
      for(MInt dir = 0; dir < nDim; dir++) {
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          // Vertices are located at the ends of the L array
          // so the position is either 0 or p + 1
          if(Fd::vertexPosition(node, d) < F0) {
            pos = 0;
          } else {
            pos = p + 1;
          }
          if(dir == d) {
            L_prime(node, dir) *= L_deriv(d, pos);
          } else {
            L_prime(node, dir) *= L(d, pos);
          }
        }
      }
      // edges
    } else if(node < noVertices + noEdges * p) {
      // Which edge do we consider right now?
      MInt edge = (MInt)floor((node - noVertices) / p);
      // Which point of the edge do we consider right now?
      MInt l = (node - noVertices) % p;
      for(MInt dir = 0; dir < nDim; dir++) {
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          if(Fd::edgePosition(edge, d) < F0) {
            pos = 0;
          } else if(Fd::edgePosition(edge, d) > F0) {
            pos = p + 1;
          } else {
            pos = l + 1;
          }
          if(dir == d) {
            L_prime(node, dir) *= L_deriv(d, pos);
          } else {
            L_prime(node, dir) *= L(d, pos);
          }
        }
      }
      // surfaces
    } else if(node < noVertices + noEdges * p + noSurfaces * p * p) {
      // Which surface do we consider right now?
      MInt surface = (MInt)floor((node - noVertices - noEdges * p) / (p * p));
      // Which point of the surface do we consider right now?
      MInt l[2] = {0};
      l[1] = (node - noVertices - noEdges * p) % (p * p);
      while(l[1] >= p) {
        l[0]++;
        l[1] -= p;
      }
      for(MInt dir = 0; dir < nDim; dir++) {
        MInt counter = 0;
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          if(Fd::surfacePosition(surface, d) < F0) {
            pos = 0;
          } else if(Fd::surfacePosition(surface, d) > F0) {
            pos = p + 1;
          } else {
            pos = l[counter] + 1;
            counter++;
          }
          if(dir == d) {
            L_prime(node, dir) *= L_deriv(d, pos);
          } else {
            L_prime(node, dir) *= L(d, pos);
          }
        }
      }
      // inner points
    } else {
      MInt l[3] = {0};
      // Which point of the inner points do we consider right now?
      l[2] = node - noVertices - noEdges * p - noSurfaces * p * p;
      while(l[2] >= (p * p)) {
        l[0]++;
        l[2] -= (p * p);
      }
      while(l[2] >= p) {
        l[1]++;
        l[2] -= p;
      }
      for(MInt dir = 0; dir < nDim; dir++) {
        MInt counter = 0;
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = l[counter] + 1;
          counter++;
          if(dir == d) {
            L_prime(node, dir) *= L_deriv(d, pos);
          } else {
            L_prime(node, dir) *= L(d, pos);
          }
        }
      }
    }
  }
}


/** \brief Calculate the strain interpolation matrix
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 * The following is done in this function:
 *
 * - Calculate the jacobian matrix and the derivative of the lagrangian polynomial
 * - The matrix product of both, the jacobian matrix and the lagrangian polynomial
 *   give the entries in the strain interpolation matrix
 **/
template <MInt nDim>
void FcSolver<nDim>::getStrainInterpolationMatrix(MFloatScratchSpace& Be, const MInt pCellId, const MFloat* z) {
  TRACE();

  MFloatScratchSpace L_prime(a_noNodes(pCellId), nDim, AT_, "L_prime");
  MFloatScratchSpace L_mul_Xinverse(a_noNodes(pCellId), nDim, AT_, "L_mul_Jac");

  // initialize L_prime with 1.0 and L_mul_Xinverse with 0.0
  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d = 0; d < nDim; d++) {
      L_prime(n, d) = F1;
      L_mul_Xinverse(n, d) = F0;
    }
  }
  getDerivativeOfDisplacementInterpolationMatrix(pCellId, z, L_prime);

  // Calculate the jacobian matrix
  MFloatScratchSpace X_inverse(nDim * nDim, AT_, "X_inverse");
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      if(d1 == d2) {
        X_inverse(d2 + d1 * nDim) = a_invJacobian(pCellId);
      } else {
        X_inverse(d2 + d1 * nDim) = F0;
      }
    }
  }
  if(m_testRun && m_polyDeg < 1) getJacobiMatrixDebug(pCellId, z, X_inverse);

  // Calculate the matrix product of L_prime and X_inverse
  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        L_mul_Xinverse(n, d1) += L_prime(n, d2) * X_inverse(d2 + d1 * nDim);
      }
    }
  }

  // Set the entries of the strain inverpolation matrix
  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        if(d1 == d2) {
          Be(d2, n * nDim + d1) = L_mul_Xinverse(n, d2);
        } else {
          Be(d2, n * nDim + d1) = F0;
        }
      }
      if(nDim == 2) {
        MInt pos = ((d1 + 1) < nDim) ? (d1 + 1) : (d1 - 1);
        Be(nDim, n * nDim + d1) = L_mul_Xinverse(n, pos);
      } else {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          MInt d3 = ((d2 + 1) < nDim) ? (d2 + 1) : 0;
          MInt d4 = ((d1 + 1) < nDim) ? (d1 + 1) : 0;
          if(d1 == d2) {
            Be(d2 + nDim, n * nDim + d1) = L_mul_Xinverse(n, d3);
          } else if(d2 == d4) {
            Be(d2 + nDim, n * nDim + d1) = F0;
          } else {
            Be(d2 + nDim, n * nDim + d1) = L_mul_Xinverse(n, d2);
          }
        }
      }
    }
  }
  if(m_testRun && m_polyDeg < 1) getStrainInterpolationMatrixDebug(pCellId, z, Be);
}

/** \brief Returns the material tensor (for now only valid for linear
 *   elastic material)
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::getMaterialTensor(MFloatScratchSpace& C, const MInt pCellId, const MInt node) {
  TRACE();

  switch(string2enum(solverMethod())) {
    case MAIA_LINEAR_ELASTIC: {
      const MFloat nu = a_poissonRatio(pCellId);
      const MFloat c1 = m_E * (F1 - nu) / ((F1 + nu) * (F1 - F2 * nu));
      const MFloat c2 = m_E * nu / ((F1 + nu) * (F1 - F2 * nu));
      const MFloat c3 = m_E / (F2 * (F1 + nu));
      if(m_first && m_testRun) {
        cout << "///////////////////////////////////////////////////////////////" << endl;
        cout << "///////////////////////////////////////////////////////////////" << endl;
        cout << endl;
        cout << "MATERIAL PROPERTIES: " << endl;
        cout << endl;
        cout << "E-Modul E = " << m_E << endl;
        cout << "Poisson Ratio nu = " << nu << endl;
        cout << "c1 = E * (F1 - nu) / ((F1 + nu) * (F1 - F2 * nu)) = " << c1 << endl;
        cout << "c2 = E * nu / ((F1 + nu) * (F1 - F2 * nu)) = " << c2 << endl;
        cout << "c3 = E / (F2 * (F1 + nu)) = " << c3 << endl;
        cout << endl;
        m_log << "///////////////////////////////////////////////////////////////" << endl;
        m_log << "///////////////////////////////////////////////////////////////" << endl;
        m_log << endl;
        m_log << "MATERIAL PROPERTIES: " << fixed << setprecision(5) << endl;
        m_log << endl;
        m_log << "E-Modul E = " << m_E << endl;
        m_log << "Poisson Ratio nu = " << nu << endl;
        m_log << "c1 = E * (F1 - nu) / ((F1 + nu) * (F1 - F2 * nu)) = " << c1 << endl;
        m_log << "c2 = E * nu / ((F1 + nu) * (F1 - F2 * nu)) = " << c2 << endl;
        m_log << "c3 = E / (F2 * (F1 + nu)) = " << c3 << endl;
        m_log << endl;
      }

      for(MInt dirA = 0; dirA < m_noStrains; dirA++) {
        for(MInt dirB = 0; dirB < m_noStrains; dirB++) {
          C(dirA, dirB) = F0;
        }
      }
      for(MInt dirA = 0; dirA < nDim; dirA++) {
        for(MInt dirB = 0; dirB < nDim; dirB++) {
          C(dirA, dirB) = (dirA == dirB) ? c1 : c2;
        }
      }
      for(MInt dirA = nDim; dirA < m_noStrains; dirA++) {
        for(MInt dirB = nDim; dirB < m_noStrains; dirB++) {
          C(dirA, dirB) = (dirA == dirB) ? c3 : F0;
        }
      }
      if(m_first && m_testRun) {
        cout << "Material tensor: " << endl;
        m_log << "Material tensor: " << endl;
        for(MInt dirA = 0; dirA < m_noStrains; dirA++) {
          for(MInt dirB = 0; dirB < m_noStrains; dirB++) {
            cout << C(dirA, dirB) << "\t\t\t";
            m_log << C(dirA, dirB) << "\t\t\t";
          }
          cout << endl;
          m_log << endl;
        }
        cout << endl;
        m_log << endl;
      }
      m_first = false;
      break;
    }
    case MAIA_NONLINEAR_BAR: {
      std::array<MFloat, nDim> coordA{};
      getCoordinatesOfNode(node, pCellId, coordA.data());
      const MInt nodeIdA = a_nodeIdsLocal(pCellId, node);

      for(MInt nodeB = node + 1; nodeB < a_noNodes(pCellId); nodeB++) {
        const MInt nodeIdB = a_nodeIdsLocal(pCellId, nodeB);
        const MFloat A = F1;
        const MFloat E = m_E;
        std::array<MFloat, nDim> coordB{};
        getCoordinatesOfNode(nodeB, pCellId, coordB.data());
        std::array<MFloat, nDim> diff{};
        for(MInt d = 0; d < nDim; d++) {
          diff[d] = coordB[d] - coordA[d];
        }
        const MFloat LSq = std::inner_product(&diff[0], &diff[nDim], &diff[0], .0);
        const MFloat L = sqrt(LSq);
        const MFloat pre = E * A / (L * L * L);
        std::array<MFloat, 2 * nDim> c{};
        std::array<MInt, 2 * nDim> gGlobal{};
        std::array<MInt, 2 * nDim> gLocal{};
        std::array<MFloat, 2 * nDim> sign{};
        for(MInt d = 0; d < nDim; d++) {
          gGlobal[d] = nodeIdA * nDim + d;
          gLocal[d] = node * nDim + d;
          c[d] = diff[d];
          sign[d] = F1;
        }
        for(MInt d = nDim; d < 2 * nDim; d++) {
          gGlobal[d] = nodeIdB * nDim + (d - nDim);
          gLocal[d] = nodeB * nDim + (d - nDim);
          c[d] = diff[d - nDim];
          sign[d] = -F1;
        }
        for(MInt n1 = 0; n1 < m_noStrains; n1++) {
          for(MInt n2 = 0; n2 < m_noStrains; n2++) {
            const MFloat disA = m_totalNodalDisplacements[gGlobal[n1]];
            const MFloat disB = m_totalNodalDisplacements[gGlobal[n2]];
            const MInt nA = gLocal[n1];
            const MInt nB = gLocal[n2];
            C(nA, nB) += pre * (F1 + F1B2 * (disA - disB)) * c[n1] * c[n2] * sign[n1] * sign[n2];
          }
        }
      }
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "Error!!! Solver method not implemented." << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }
}

/** \brief Calculates the stiffness matrix in one point, i.e., at a integration point
 *  For the calculation, the strain-interpolation matrix and the material tensor are
 *  required. Furthermore, the determinant of the Jacobi matrix is needed.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::getNodalStiffnessMatrix(const MInt pCellId, const MInt node, MFloatScratchSpace& Ke, MInt* nodePos,
                                             MFloat* z, MFloat determinant, const MFloat alpha,
                                             const MInt subCellLevel) {
  TRACE();

  for(MInt d1 = 0; d1 < nDim * a_noNodes(pCellId); d1++) {
    for(MInt d2 = 0; d2 < nDim * a_noNodes(pCellId); d2++) {
      Ke(d1, d2) = F0;
    }
  }

  if(solverMethod() == "MAIA_NONLINEAR_BAR") {
    getMaterialTensor(Ke, pCellId, node);
  } else {
    MFloatScratchSpace Be(m_noStrains, nDim * a_noNodes(pCellId), AT_, "Be");
    MFloatScratchSpace BeT(nDim * a_noNodes(pCellId), m_noStrains, AT_, "BeT");
    MFloatScratchSpace BeT_mul_C(nDim * a_noNodes(pCellId), m_noStrains, AT_, "BeT_mul_C");

    // Get the material tensor
    MFloatScratchSpace C(m_noStrains, m_noStrains, AT_, "C");
    getMaterialTensor(C, pCellId, node);

    // Get strain interpolation matrix and transpose it
    getStrainInterpolationMatrix(Be, pCellId, z);
    for(MInt d = 0; d < m_noStrains; d++) {
      for(MInt e = 0; e < a_noNodes(pCellId) * nDim; e++) {
        BeT(e, d) = Be(d, e);
      }
    }

    // Calculate the integration weigths
    MFloat weight = F1;
    for(MInt d = 0; d < nDim; d++) {
      weight *= m_gaussWeights[a_pRfnmnt(pCellId)][nodePos[d]];
    }

    if(m_testRun) {
      m_volume += determinant * weight * alpha;
    }
    // Multiply the integration weights and the determinant to the material tensor
    for(MInt d1 = 0; d1 < m_noStrains; d1++) {
      for(MInt d2 = 0; d2 < m_noStrains; d2++) {
        C(d1, d2) *= determinant * weight * alpha;
      }
    }

    if(m_testRun && subCellLevel < 1) {
      cout << "Strain-Interpolation-Matrix B^T_e at this point:" << endl;
      for(MInt e = 0; e < a_noNodes(pCellId) * nDim; e++) {
        for(MInt d = 0; d < m_noStrains; d++) {
          cout << BeT(e, d) << " ";
        }
        cout << endl;
      }
      cout << endl;
      cout << "Material tensor multiplied with the determinant " << determinant << " and the integration weights w "
           << weight << " and alpha " << alpha << ":" << endl;
      for(MInt d1 = 0; d1 < m_noStrains; d1++) {
        for(MInt d2 = 0; d2 < m_noStrains; d2++) {
          cout << C(d1, d2) << " ";
        }
        cout << endl;
      }
      cout << endl;
      cout << "Strain-Interpolation-Matrix B_e at this point:" << endl;
      for(MInt d = 0; d < m_noStrains; d++) {
        for(MInt e = 0; e < a_noNodes(pCellId) * nDim; e++) {
          cout << Be(d, e) << " ";
        }
        cout << endl;
      }
      cout << endl;
    }

    // Matrix Multiplication
    maia::math::multiplyMatrices(BeT, C, BeT_mul_C, nDim * a_noNodes(pCellId), m_noStrains, m_noStrains, m_noStrains);

    maia::math::multiplyMatrices(BeT_mul_C, Be, Ke, nDim * a_noNodes(pCellId), m_noStrains, m_noStrains,
                                 nDim * a_noNodes(pCellId));

    if(m_testRun && subCellLevel < 1) {
      cout << "Intermediate result, B^T_e * C * B_e * det(X_,z) * w:" << endl;
      for(MInt d = 0; d < a_noNodes(pCellId) * nDim; d++) {
        for(MInt e = 0; e < a_noNodes(pCellId) * nDim; e++) {
          cout << Ke(d, e) << " ";
        }
        cout << endl;
      }
      cout << endl;
      cout << "---------------------------------------" << endl;
      cout << endl;
    }
  }
}

/** \brief Calculates the element stiffness matrix by summing up the nodel stiffness matrices.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::getElementMatrix(MFloatScratchSpace& Ke, const MInt pCellId, const MFloat alpha,
                                      const MInt subCellLevel, const MFloat* subCellCoord) {
  TRACE();

  MInt nodePos[nDim] = {0};
  MFloat z[nDim] = {F0};
  MFloat determinant = F1;
  for(MInt d = 0; d < nDim; d++) {
    determinant *= a_jacobianMatrix(pCellId) * FFPOW2(subCellLevel);
  }

  if(subCellLevel < 1 && m_testRun) {
    cout << endl;
    cout << "The Element-Stiffness-Matrix is calculated " << endl;
    cout << "by adding together the resulting matrices of the product " << endl;
    cout << endl;
    cout << "B^T_e(z) * C(z) * B_e(z) * det(X_,z(z))" << endl;
    cout << endl;
    cout << "calculated at each of the following integration points z:" << endl;
    cout << endl;
    cout << "---------------------------------------" << endl;
    cout << endl;
  }

  for(MInt node = 0; node < a_noNodes(node); node++) {
    MFloatScratchSpace Mat_part(a_noNodes(pCellId) * nDim, a_noNodes(pCellId) * nDim, AT_, "Mat_part");

    const MInt nodesPerDir = a_pRfnmnt(pCellId) + 2;
    Fd::nodePosition(node, nodesPerDir, nodePos);

    for(MInt d = 0; d < nDim; d++) {
      z[d] = m_gaussPoints[a_pRfnmnt(pCellId)][nodePos[d]];
    }

    if(subCellLevel > 0) {
      std::array<MFloat, nDim> point{};
      interpolateSubCellGeometry(pCellId, z, point.data(), subCellLevel, subCellCoord);
      MBool outside = m_geometry->pointIsInside2(point.data());
      if(outside) {
        continue;
      } else {
        transformToLocal(pCellId, point.data(), z);
      }
    }

    if(m_testRun) {
      if(subCellLevel < 1) {
        cout << "Integration points used: ";
        for(MInt d = 0; d < nDim; d++) {
          if(subCellLevel < 1) {
            cout << z[d] << " ";
          }
        }
        cout << endl;
      }

      MFloatScratchSpace x(nDim, nDim, AT_, "x");
      for(MInt d1 = 0; d1 < nDim; d1++) {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          x(d1, d2) = F0;
        }
      }
      lagrangianPointsJacobi(pCellId, z, x);
      std::array<std::array<MFloat, nDim>, nDim> x_{};
      for(MInt d1 = 0; d1 < nDim; d1++) {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          x_[d1][d2] = x(d1, d2);
        }
      }
      MFloat determinantDebug = maia::math::determinant(x_);
      for(MInt d = 0; d < nDim; d++) {
        determinantDebug *= FFPOW2(subCellLevel);
      }
      if(!approx(determinantDebug, determinant, m_eps)) {
        stringstream errorMessage;
        errorMessage << "Error in calculating the determinant of the jaconbian matrix" << determinant << " "
                     << determinantDebug << endl;
        mTerm(1, AT_, errorMessage.str());
      }
      if(subCellLevel < 1) {
        cout << "determinant at this point: " << determinant << endl;
      }
    }

    getNodalStiffnessMatrix(pCellId, node, Mat_part, nodePos, z, determinant, alpha, subCellLevel);
    for(MInt m1 = 0; m1 < a_noNodes(pCellId) * nDim; m1++) {
      for(MInt m2 = 0; m2 < a_noNodes(pCellId) * nDim; m2++) {
        Ke(m1, m2) += Mat_part(m1, m2);
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::getInternalLoadsFromStresses() {
  TRACE();

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;
    MFloat z[nDim] = {F0};
    MInt nodePos[nDim] = {0};

    MFloat determinant = F1;
    for(MInt d = 0; d < nDim; d++) {
      determinant *= a_jacobianMatrix(cell);
    }

    for(MInt node = 0; node < a_noNodes(cell); node++) {
      const MInt nodesPerDir = a_pRfnmnt(cell) + 2;
      Fd::nodePosition(node, nodesPerDir, nodePos);

      for(MInt d = 0; d < nDim; d++) {
        z[d] = m_gaussPoints[a_pRfnmnt(cell)][nodePos[d]];
      }

      MFloat weight = F1;
      for(MInt d = 0; d < nDim; d++) {
        weight *= m_gaussWeights[a_pRfnmnt(cell)][nodePos[d]];
      }
      weight *= determinant;

      MFloatScratchSpace Be(m_noStrains, nDim * a_noNodes(cell), AT_, "Be");

      // Get strain interpolation matrix and transpose it
      getStrainInterpolationMatrix(Be, cell, z);

      std::array<MFloat, m_noStrains> stresses{};
      if(string2enum(solverMethod()) == MAIA_LINEAR_ELASTIC) {
        getStressesAtPoint(cell, z, m_totalNodalDisplacements, stresses.data());
      } else {
        for(MInt d = 0; d < m_noStrains; d++) {
          stresses[d] = a_nodalStresses(cell, node, d);
        }
      }
      // Loop over all directions and all integration points to get the
      // load at the nodes
      for(MInt n = 0; n < a_noNodes(cell); n++) {
        for(MInt dim = 0; dim < nDim; dim++) {
          MInt nodeId = a_nodeIdsLocal(cell, n);
          MFloat load = F0;
          for(MInt d = 0; d < m_noStrains; d++) {
            MFloat BeT = Be(d, n * nDim + dim);
            load += BeT * stresses[d];
          }
          m_internalLoadVector[nodeId * nDim + dim] += load * weight;
        }
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::getReactionForces() {
  TRACE();

  m_bndryCnd->calcReactionForces();
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::rhs() {
  TRACE();
}

/** \brief Applies von Neumann Boundary conditions, e.g.
 * body forces.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::rhsBnd() {
  TRACE();

  // Apply boundary condition
  m_bndryCnd->updateForceVector();
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::lhs() {
  TRACE();
}

/** \brief Applies Dirichlet Boundary conditions, e.g.
 * displacements.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::lhsBnd() {
  TRACE();

  m_globalBndryMatrix.clear();

  // Apply boundary condition
  m_bndryCnd->updateSystemMatrix();
}

/** \brief Initialize the solver and the boundary
 * conditions.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::initSolver() {
  TRACE();
}

/** \brief  Initialize all timer groups, timer, and sub timer required by the LB solver
 *  \author Miro Gondrum
 *  \date   18.02.2022
 */
template <MInt nDim>
void FcSolver<nDim>::initTimer() {
  TRACE();

  NEW_TIMER_GROUP(tg_solver, "FC Solver (solverId=" + std::to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_t.solver, "complete solver", tg_solver);

  NEW_SUB_TIMER_NOCREATE(m_t.initSolver, "init solver", m_t.solver);
  NEW_SUB_TIMER_NOCREATE(m_t.calcStiffMat, "calc Siffness Matrix per Level", m_t.initSolver);
  NEW_SUB_TIMER_NOCREATE(m_t.solutionStep, "SolutionStep", m_t.solver);
  NEW_SUB_TIMER_NOCREATE(m_t.subCellInt, "subcell Integration", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.forceBC, "boundary forces", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.displacementBC, "boundary displacement", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.solving, "Solution of System of Equations", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.exchange, "Exchange", m_t.solutionStep);
}

/** \brief  Average all LB timer over participating ranks
 *  \author Miro Gondrum
 *  \date   18.02.2022
 */
template <MInt nDim>
void FcSolver<nDim>::averageTimer() {
  TRACE();
  if(!grid().isActive()) return;
  // 0) map timer ids for safety
  std::vector<MInt> timerIds_;
  timerIds_.reserve(17);
  timerIds_.emplace_back(m_t.solver);
  timerIds_.emplace_back(m_t.initSolver);
  timerIds_.emplace_back(m_t.calcStiffMat);
  timerIds_.emplace_back(m_t.solutionStep);
  timerIds_.emplace_back(m_t.subCellInt);
  timerIds_.emplace_back(m_t.forceBC);
  timerIds_.emplace_back(m_t.displacementBC);
  timerIds_.emplace_back(m_t.solving);
  timerIds_.emplace_back(m_t.exchange);
  const MInt noTimers = timerIds_.size();
  // 1) fill buffer with local timer values
  std::vector<MFloat> timerValues_;
  timerValues_.reserve(noTimers);
  for(MInt i = 0; i < noTimers; i++) {
    timerValues_.emplace_back(RETURN_TIMER_TIME(timerIds_[i]));
  }
  // 2) collect values from all ranks
  MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_SUM, mpiComm(),
                AT_, "MPI_IN_PLACE", "timerValues_");
  // 3) perform averaging on timer and4) set new timer values
  const MInt noDomains_ = noDomains();
  for(MInt i = 0; i < noTimers; i++) {
    const MFloat meanValue = timerValues_[i] / noDomains_;
    SET_RECORD(timerIds_[i], meanValue);
  }
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::finalizeInitSolver() {
  TRACE();
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::preTimeStep() {
  TRACE();
}

/** \brief Solution step of the FC solver. The element stiffness matrix
 * is calculated for each cell. The element stiffness matrices are assembled,
 * boundary conditions are set and the system of equations is solved.
 *
 * \author Moritz Waldmann
 * \date 20.03.21
 *
 **/
template <MInt nDim>
MBool FcSolver<nDim>::solutionStep() {
  TRACE();

  if(m_testRun) {
    // For testing, fill the Displacement vector with the point coordinates
    // Thus, the cell strains must be unity the shear components vanish
    fillDisplacementVectorWithCoords();
    calcElementDisplacements();
    for(MInt i = 0; i < a_noCells(); i++) {
      if(!a_isActive(i)) continue;
      for(MInt d = 0; d < nDim; d++) {
        if(!approx(c_coordinate(i, d), a_elementDispl(i, d), 1e-14)) {
          stringstream errorMessage;
          errorMessage << "Error in interpolating the strains " << c_coordinate(i, d) << " " << a_elementDispl(i, d)
                       << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    }
    resetDisplacements();
  }

  switch(string2enum(solverMethod())) {
    case MAIA_LINEAR_ELASTIC: {
      // Reset vectors
      resetDisplacements(0);
      resetLoadVector(0);

      // Calculate stiffness matrix
      if(globalTimeStep == 1) {
        calculateStiffnessMatrix();
      }
      lhsBnd();
      assembleFinalMatrix();
      createCompressedMatrix();

      // Calculate Force vector
      rhsBnd();
      getInternalLoadsFromStresses();
      getReactionForces();
      updateLoadVector();

      // Solve system of equations
      solveSystemOfEquations();

      // Update total displacements and calculate strains and stresses
      m_bndryCnd->writeOutModifiedBoundaries();
      calcElementDisplacements();
      calcElementStrains();
      calcElementStresses();

      break;
    }
    case MAIA_NONLINEAR_BAR: {
      // Reset vectors and calculate external forces
      resetLoadVector(-1);
      rhsBnd();

      for(MInt step = 0; step < m_noLoadSteps; step++) {
        // Reset Displacements
        resetDisplacements(0);

        // Calculate force vector
        const MFloat lambda = ((MFloat)step + F1) / (MFloat)m_noLoadSteps;
        updateLoadVector(lambda);

        MInt k = 0;
        MFloat error = F1;

        while(error > 1e-8 && k < 5) {
          // Calculate stiffness matrix
          calculateStiffnessMatrix();
          lhsBnd();
          assembleFinalMatrix();
          createCompressedMatrix();

          // Solve system of equations
          solveSystemOfEquations();

          // Calculate internal forces and reaction forces
          resetLoadVector(1);

          // Calculate stiffness matrix
          calculateStiffnessMatrix();
          m_globalBndryMatrix.clear();
          assembleFinalMatrix();
          createCompressedMatrix();

          const MInt n = m_finalGlobalMatrix.size();
          const MInt m = (m_totalNoNodes * nDim);
          maia::math::multiplySparseMatrixVector(m_sysMatCompressed, m_compressionIndexSys, n, m, m_internalLoadVector,
                                                 m_totalNodalDisplacements);

          for(MInt i = 0; i < m_totalNoNodes; i++) {
            cout << "Displacements (node " << i << ")";
            for(MInt d = 0; d < nDim; d++) {
              cout << std::setprecision(8) << " " << m_totalNodalDisplacements[i * nDim + d];
            }
            cout << endl;
          }
          cout << endl;

          for(MInt i = 0; i < m_totalNoNodes; i++) {
            cout << "External forces (node " << i << ")";
            for(MInt d = 0; d < nDim; d++) {
              cout << std::setprecision(8) << " " << m_nodalLoadVector[i * nDim + d];
            }
            cout << endl;
          }
          cout << endl;

          for(MInt i = 0; i < m_totalNoNodes; i++) {
            cout << "Internal forces (node " << i << ")";
            for(MInt d = 0; d < nDim; d++) {
              cout << std::setprecision(8) << " " << m_internalLoadVector[i * nDim + d];
            }
            cout << endl;
          }
          cout << endl;

          getReactionForces();
          for(MInt i = 0; i < m_totalNoNodes; i++) {
            cout << "Reaction forces (node " << i << ")";
            for(MInt d = 0; d < nDim; d++) {
              cout << std::setprecision(8) << " " << m_reactionForceVector[i * nDim + d];
            }
            cout << endl;
          }
          cout << endl;

          updateLoadVector(lambda);

          const MFloat norm2Res = maia::math::norm(m_nodalLoadVector, m);
          const MFloat norm2Int = maia::math::norm(m_internalLoadVector, m);
          const MFloat norm2Reac = maia::math::norm(m_reactionForceVector, m);

          error = (fabs(norm2Int) > m_eps) ? norm2Res / (norm2Int) : norm2Res;
          if(m_testRun) {
            cout << fixed << setprecision(15) << endl;
            cout << "Error in step " << step << " and iteration " << k << " is: " << error << endl;
            cout << "Lambda in step " << step << " and iteration " << k << " is: " << lambda << endl;
            cout << "Residual in step " << step << " and iteration " << k << " is: " << norm2Res << endl;
            cout << "IntForce vector in step " << step << " and iteration " << k << " is: " << norm2Int << endl;
            cout << "ReacForce vector in step " << step << " and iteration " << k << " is: " << norm2Reac << endl;
            cout << fixed << setprecision(5) << endl;
          }
          k++;
        }
      }
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "Error!!! Solver method not implemented." << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  return true;
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::saveSolverSolution(MBool forceOutput, const MBool finalTimeStep) {
  TRACE();
  std::ignore = forceOutput;
  std::ignore = finalTimeStep;
}

// Not required at the moment
template <MInt nDim>
void FcSolver<nDim>::postTimeStep() {
  TRACE();
}

/**
 * \brief This function prepares a restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */
template <MInt nDim>
MBool FcSolver<nDim>::prepareRestart(MBool writeRestart, MBool& writeGridRestart) {
  TRACE();

  writeGridRestart = false;

  if(((globalTimeStep % m_restartInterval) == 0) || writeRestart) {
    writeRestart = true;

    if(m_adaptationSinceLastRestart) {
      writeGridRestart = true;
    }
  }

  return writeRestart;
}

/**
 * \brief This function writes restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */
template <MInt nDim>
void FcSolver<nDim>::writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString gridFileName,
                                      MInt* recalcIdTree) {
  TRACE();
  std::ignore = writeBackup;

  if(writeRestart) {
    stringstream fileName;
    fileName << outputDir() << "restart_" << globalTimeStep << ParallelIo::fileExt();

    if(m_recalcIds != nullptr) {
      // for multisolver recalc size needs to be raw grid size?
      for(MInt cellId = 0; cellId < maxNoGridCells(); cellId++) {
        m_recalcIds[cellId] = recalcIdTree[cellId];
      }

      MIntScratchSpace recalcIdsSolver(grid().tree().size(), AT_, "recalcIds");

      MInt recalcCounter = 0;
      for(MInt cell = 0; cell < grid().raw().m_noInternalCells; cell++) {
        // cerr <<"cell : " << cell << " recalc: " << m_recalcIds[cell] << endl;
        if(grid().raw().a_hasProperty(cell, Cell::IsHalo)) {
          continue;
        }
        MInt sId = grid().tree().grid2solver(m_recalcIds[cell]);
        if(sId > -1) {
          recalcIdsSolver[recalcCounter] = sId;
          recalcCounter++;
        }
      }
      ASSERT(recalcCounter == grid().noInternalCells(), "recalc ids size is wrong");

      saveRestartStrainsOnly(fileName.str().c_str(), gridFileName.c_str(), recalcIdsSolver.getPointer());
      cerr << "RecalcIds" << endl;
    } else {
      saveRestartStrainsOnly(fileName.str().c_str(), gridFileName.c_str(), m_recalcIds);
      cerr << "No RecalcIds" << endl;
    }
  }
}

/** \brief This function stores the flow information of the cells
 *        such as variables and attributes like u_velocity,density,etc.,
 *
 * \author Moritz Waldmann
 * \date 14.06.2021
 * In contrast to saveUvwOnly(const MChar* fileName) this function also stores the
 * information for restart in parallel.
 * The attribute 'name' of the variables are set to their according meaning.
 *
 * \param[in] fileName the name of the file to write to
 **/
template <MInt nDim>
void FcSolver<nDim>::saveRestartStrainsOnly(const MChar* fileName, const MChar* gridInputFileName, MInt* recalcIdTree) {
  TRACE();
  if(nDim != 2 && nDim != 3) {
    cerr << " In global function saveRestartStrainsOnly: wrong number of dimensions !" << endl;
    exit(0);
    return;
  }

  cerr << "Writing solution to file" << endl;
  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());
  parallelIo.defineScalar(PIO_INT, "noCells");

  // total #of cells in the grid == global cell index of the last cell in the last process:
  const MInt totalNoCells = domainOffset(noDomains());
  cerr << "totalNoCells " << totalNoCells << endl;

  //--specify helper functions
  auto defFloatArray = [&](const MString arrayName, const MString varName, const MInt length) {
    parallelIo.defineArray(PIO_FLOAT, arrayName, length);
    parallelIo.setAttribute(varName, "name", arrayName);
  };
  // function: define macroscopic variables
  auto defineMacroscopicVariables = [&](const MString name, const MString prefix) {
    // velocities
    MInt noIds = (m_testRun) ? IPOW2(nDim) + 1 : 0;
    const MInt noOutputVar = (nDim == 2) ? (6 + nDim + noIds) : (12 + nDim + noIds);
    MStringScratchSpace velNames(noOutputVar, AT_, "velNames");
    if(nDim == 2) {
      velNames[0] = "Epsilon_x";
      velNames[1] = "Epsilon_y";
      velNames[2] = "Epsilon_xy";
      velNames[3] = "Stress_x";
      velNames[4] = "Stress_y";
      velNames[5] = "Stress_xy";
      velNames[6] = "Displacement_x";
      velNames[7] = "Displacement_y";
      if(m_testRun) {
        velNames[5] = "GlobalNodeId0";
        velNames[6] = "GlobalNodeId1";
        velNames[7] = "GlobalNodeId2";
        velNames[8] = "GlobalNodeId3";
        velNames[9] = "DomainId";
      }
    } else {
      velNames[0] = "Epsilon_x";
      velNames[1] = "Epsilon_y";
      velNames[2] = "Epsilon_z";
      velNames[3] = "Epsilon_xy";
      velNames[4] = "Epsilon_yz";
      velNames[5] = "Epsilon_xz";
      velNames[6] = "Stress_x";
      velNames[7] = "Stress_y";
      velNames[8] = "Stress_z";
      velNames[9] = "Stress_xy";
      velNames[10] = "Stress_yz";
      velNames[11] = "Stress_xz";
      velNames[12] = "Displacement_x";
      velNames[13] = "Displacement_y";
      velNames[14] = "Displacement_z";
      if(m_testRun) {
        velNames[15] = "GlobalNodeId0";
        velNames[16] = "GlobalNodeId1";
        velNames[17] = "GlobalNodeId2";
        velNames[18] = "GlobalNodeId3";
        velNames[19] = "GlobalNodeId4";
        velNames[20] = "GlobalNodeId5";
        velNames[21] = "GlobalNodeId6";
        velNames[22] = "GlobalNodeId7";
        velNames[23] = "DomainId";
      }
    }

    for(MInt d = 0; d != noOutputVar; d++) {
      defFloatArray(name + std::to_string(d), prefix + velNames[d], totalNoCells);
    }
  };
  // function: write macroscopic variable
  auto writeMacroscopicStrain = [&](const MInt index) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_elementStrains(recalcIdTree != nullptr ? recalcIdTree[i] : i, index);
    }
    parallelIo.writeArray(tmp.getPointer(), "strain" + std::to_string(index));
  };
  // function: write macroscopic variable
  auto writeMacroscopicStress = [&](const MInt d, const MInt index) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_elementStresses(recalcIdTree != nullptr ? recalcIdTree[i] : i, d);
    }
    parallelIo.writeArray(tmp.getPointer(), "strain" + std::to_string(index));
  };
  // function: write macroscopic variable
  auto writeMacroscopicDisplacement = [&](const MInt d, const MInt index) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_elementDispl(recalcIdTree != nullptr ? recalcIdTree[i] : i, d);
    }
    parallelIo.writeArray(tmp.getPointer(), "strain" + std::to_string(index));
  };
  // function: write macroscopic variable
  auto writeGlobalNodeIds = [&](const MInt d, const MInt index) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_nodeIdsGlobal(recalcIdTree != nullptr ? recalcIdTree[i] : i, d);
    }
    parallelIo.writeArray(tmp.getPointer(), "strain" + std::to_string(index));
  };
  // function: write macroscopic variable
  auto writeDomainId = [&](const MInt d, const MInt index) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    std::ignore = d;
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = (MFloat)domainId();
    }
    parallelIo.writeArray(tmp.getPointer(), "strain" + std::to_string(index));
  };
  //--define arrays
  defineMacroscopicVariables("strain", "");

  //--define global attributes
  parallelIo.setAttribute(solverId(), "solverId");
  parallelIo.setAttribute(gridInputFileName, "gridFile", "");
  parallelIo.setAttribute(globalTimeStep, "globalTimeStep");
  parallelIo.setAttribute(totalNoCells, "noCells");

  //--write arrays
  // Set file offsets (first globalId and #of cells to be written by this process):
  const MPI_Offset firstGlobalId = domainOffset(domainId());
  const MPI_Offset localNoCells = noInternalCells();
  parallelIo.setOffset(localNoCells, firstGlobalId);

  cerr << "InternalCells " << localNoCells << endl;
  cerr << "Write out element strains" << endl;
  for(MInt index = 0; index < m_noStrains; index++) {
    writeMacroscopicStrain(index);
  }
  cerr << "Write out element stresses" << endl;
  for(MInt index = 0; index < m_noStrains; index++) {
    writeMacroscopicStress(index, index + m_noStrains);
  }
  cerr << "Write out element displacements" << endl;
  for(MInt d = 0; d < nDim; d++) {
    writeMacroscopicDisplacement(d, d + 2 * m_noStrains);
  }
  if(m_testRun) {
    cerr << "Write out global node Ids" << endl;
    for(MInt d = 0; d < IPOW2(nDim); d++) {
      writeGlobalNodeIds(d, d + nDim + 2 * m_noStrains);
    }
    writeDomainId(0, IPOW2(nDim) + nDim + 2 * m_noStrains);
  }
}

/**
 * \brief This function resets the grid-trigger after a restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */
template <MInt nDim>
void FcSolver<nDim>::reIntAfterRestart(MBool doneRestart) {
  TRACE();

  if(doneRestart) {
    m_adaptationSinceLastRestart = false;
  }
}

/**
 * \brief Sets the cell-weight for balancing and a restarting
 * \author Moritz Waldmann
 */
template <MInt nDim>
void FcSolver<nDim>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt offset = noCellsGrid * solverId();

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    const MInt id = gridCellId + offset;
    solverCellWeight[id] = F1; // 0.2;
  }
}

/**
 * \brief Computes the Assembled system matrix by summing the element stiffness matrices
 *
 * \author Moritz Waldmann
 * \date 10.01.22
 *
 * The indices of the global stiffness matrix are based on the globalId of each node.
 */
template <MInt nDim>
void FcSolver<nDim>::computeAssembledSystemMatrix(MFloatScratchSpace& Ke, const MInt pCellId) {
  TRACE();

  for(MInt m1 = 0; m1 < a_noNodes(pCellId); m1++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt m2 = 0; m2 < a_noNodes(pCellId); m2++) {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          // Calculate node ids (global und local)
          MInt nodeIdA = a_nodeIdsLocal(pCellId, m1); // nodeId of all nodes
          MInt nodeIdB = a_nodeIdsLocal(pCellId, m2); // nodeId of all nodes
          MInt rowPosCell = m1 * nDim + d1;           // position in Stiffness Matrix of Cell
          MInt colPosCell = m2 * nDim + d2;           // position in Stiffness Matrix of Cell
          MInt rowPosRank = nodeIdA * nDim + d1;      // position in assembled Stiffness Matrix
          MInt colPosRank = nodeIdB * nDim + d2;      // position in assembled Stiffness Matrix

          // add to global stiffness matrix in vector notation
          // store only non zero entries
          if(!approx(Ke(rowPosCell, colPosCell), F0, m_eps)) {
            // search in vector if the element already exists
            auto pos = std::make_pair(rowPosRank, colPosRank);
            auto iter = m_globalStiffnessMatrix.find(pos);
            // if the element exists add new value to the existing one
            if(iter != m_globalStiffnessMatrix.end()) {
              iter->second += Ke(rowPosCell, colPosCell);

              // if the element does not exist insert it
              // at the correct position in the vector
            } else {
              MFloat value = Ke(rowPosCell, colPosCell);
              m_globalStiffnessMatrix.insert({pos, value});
            }
          }
        }
      }
    }
  }
}

/**
 * \brief Computes the Assembled boundary matrix by summing the element boundary matrices
 *
 * \author Moritz Waldmann
 * \date 10.01.22
 *
 * The indices of the global boundary matrix are based on the globalId of each node.
 */
template <MInt nDim>
void FcSolver<nDim>::computeAssembledBndryMatrix(MFloatScratchSpace& Bndry, const MInt pCellId) {
  TRACE();

  for(MInt m1 = 0; m1 < a_noNodes(pCellId); m1++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt m2 = 0; m2 < a_noNodes(pCellId); m2++) {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          // Calculate node ids (global und local)
          MInt nodeIdA = a_nodeIdsLocal(pCellId, m1); // nodeId of all nodes
          MInt nodeIdB = a_nodeIdsLocal(pCellId, m2); // nodeId of all nodes
          MInt rowPosCell = m1 * nDim + d1;           // position in Bndry Matrix of Cell
          MInt colPosCell = m2 * nDim + d2;           // position in Bndry Matrix of Cell
          MInt rowPosRank = nodeIdA * nDim + d1;      // position in assembled Mass Matrix
          MInt colPosRank = nodeIdB * nDim + d2;      // position in assembled Mass Matrix

          // add to global Mass matrix in vector notation
          // store only non zero entries
          if(!approx(Bndry(rowPosCell, colPosCell), F0, m_eps)) {
            // search in vector if the element already exists
            auto pos = std::make_pair(rowPosRank, colPosRank);
            auto iter = m_globalBndryMatrix.find(pos);
            // if the element exists add new value to the existing one
            if(iter != m_globalBndryMatrix.end()) {
              if(m_addingPenalties) {
                iter->second += Bndry(rowPosCell, colPosCell);
              }

              // if the element does not exist insert it
              // at the correct position in the vector
            } else {
              MFloat value = Bndry(rowPosCell, colPosCell);
              m_globalBndryMatrix.insert({pos, value});
            }
          }
        }
      }
    }
  }
}

/**
 * \brief Computes the Assembled final matrix by summing all assambled matrices
 *
 * \author Moritz Waldmann
 * \date 10.01.22
 *
 * The indices of the global matrix are based on the globalId of each node.
 */
template <MInt nDim>
void FcSolver<nDim>::assembleFinalMatrix() {
  TRACE();

  m_finalGlobalMatrix.clear();

  auto iterMat = m_globalStiffnessMatrix.begin();
  while(iterMat != m_globalStiffnessMatrix.end()) {
    auto pos = iterMat->first;
    auto value = iterMat->second;
    m_finalGlobalMatrix.insert({pos, value});
    iterMat++;
  }

  iterMat = m_globalBndryMatrix.begin();
  while(iterMat != m_globalBndryMatrix.end()) {
    auto pos = iterMat->first;
    auto value = iterMat->second;
    auto iter = m_finalGlobalMatrix.find(pos);
    // if the element exists add new value to the existing one
    if(iter != m_finalGlobalMatrix.end()) {
      iter->second += value;

      // if the element does not exist insert it
      // at the correct position in the vector
    } else {
      m_finalGlobalMatrix.insert({pos, value});
    }
    iterMat++;
  }
}

template <MInt nDim>
void FcSolver<nDim>::calculateStiffnessMatrix() {
  TRACE();

  // Allocate vector for the global stiffness matrix, which stores only non-zero entries
  MFloat* eigenValues{};
  MFloat** Ke_final{};
  m_globalStiffnessMatrix.clear();

  MBool first = true;
  for(MInt cell = 0; cell < a_noCells(); cell++) {
    m_volume = 0.0;
    if(!a_isActive(cell)) continue;

    if(m_testRun) {
      cout << "///////////////////////////////////////////////////////////////" << endl;
      cout << "///////////////////////////////////////////////////////////////" << endl;
      cout << endl;
      cout << "CellId " << cell << endl;
      cout << "Cell center: ";
      for(MInt d = 0; d < nDim; d++) {
        cout << c_coordinate(cell, d) << " ";
      }
      cout << endl;
      cout << endl;
      cout << "Nodes of cell " << cell << endl;
      std::array<MInt, nDim> nodePos{};
      for(MInt p = 0; p < a_noNodes(cell); p++) {
        cout << "Global coordinates of Point " << p << " (localId = " << a_nodeIdsLocal(cell, p)
             << ", globalId = " << a_nodeIdsGlobal(cell, p) << "): ";
        for(MInt d = 0; d < nDim; d++) {
          cout << m_nodalDisplacements[a_nodeIdsLocal(cell, p) * nDim + d] << " ";
        }
        cout << endl;
        cout << "Local coordinates of Point " << p << " (localId = " << a_nodeIdsLocal(cell, p)
             << ", globalId = " << a_nodeIdsGlobal(cell, p) << "): ";
        for(MInt d = 0; d < nDim; d++) {
          cout << a_nodePosition(cell, nodePos[d]) << " (" << nodePos[d] << ") ";
        }
        cout << endl;
        cout << endl;
        if((nodePos[nDim - 1] + 1) < (a_pRfnmnt(cell) + 2)) {
          nodePos[nDim - 1]++;
        } else {
          nodePos[nDim - 1] = 0;
          if((nodePos[nDim - 2] + 1) < (a_pRfnmnt(cell) + 2)) {
            nodePos[nDim - 2]++;
          } else {
            nodePos[nDim - 2] = 0;
            nodePos[0]++;
          }
        }
      }
    }

    // Allocate element stiffness matrix
    if(first) {
      first = false;
      mAlloc(eigenValues, a_noNodes(cell) * nDim, "eigenValues", F0, AT_);
      mAlloc(Ke_final, a_noNodes(cell) * nDim, a_noNodes(cell) * nDim, "Ke_final", F0, AT_);
    }
    MFloatScratchSpace Ke(a_noNodes(cell) * nDim, a_noNodes(cell) * nDim, AT_, "Ke");

    // reset element stiffness matrix
    for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
      for(MInt m2 = 0; m2 < a_noNodes(cell) * nDim; m2++) {
        Ke(m1, m2) = F0;
        Ke_final[m1][m2] = F0;
      }
      eigenValues[m1] = F0;
    }

    if(m_testRun) {
      // For testing of the interpolation matrix use random points. The results of x1 and x2
      // should be equal
      std::array<MFloat, nDim> z{};
      for(MInt d = 0; d < nDim; d++) {
        z[d] = (static_cast<MFloat>(rand()) / static_cast<MFloat>(RAND_MAX)) * F2 - F1;
      }
      std::array<MFloat, nDim> x1{};
      std::array<MFloat, nDim> x2{};
      interpolateGeometryDebug(cell, z.data(), x1.data());
      interpolateGeometry(cell, z.data(), x2.data());
      for(MInt d = 0; d < nDim; d++) {
        if(!approx(x1[d], x2[d], 1e-10)) {
          stringstream errorMessage;
          errorMessage << "Error in interpolating the geometry " << x1[d] << " " << x2[d] << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    }

    // Calculate Element stiffness matrix without penalty factor
    const MFloat alpha = (a_maxSubCellLvl(cell) < 1) ? F1 : m_alpha;
    const MInt subCellLevel = 0;
    getElementMatrix(Ke, cell, alpha, subCellLevel);
    for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
      for(MInt m2 = 0; m2 < a_noNodes(cell) * nDim; m2++) {
        Ke_final[m1][m2] = Ke(m1, m2);
      }
    }

    std::array<MFloat, nDim> coords{};
    for(MInt d = 0; d < nDim; d++) {
      coords[d] = c_coordinate(cell, d);
    }
    for(MInt child = 0; child < IPOW2(nDim); child++) {
      const MInt childLvl = subCellLevel + 1;
      m_bndryCnd->subCellIntegration(childLvl, coords.data(), child, cell, Ke_final);
    }
    for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
      for(MInt m2 = 0; m2 < a_noNodes(cell) * nDim; m2++) {
        Ke(m1, m2) = Ke_final[m1][m2];
      }
    }

    if(m_testRun) {
      cout << "Element-Stiffness-Matrix of Cell :" << cell << endl;
      cout << "SubCellTree depth " << a_maxSubCellLvl(cell) << endl;
      cout << "Volume is " << m_volume << endl;
      for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
        for(MInt m2 = 0; m2 < a_noNodes(cell) * nDim; m2++) {
          cout << Ke_final[m1][m2] << " ";
        }
        cout << endl;
      }
      cout << endl;
      maia::math::calcEigenValues(Ke_final, eigenValues, a_noNodes(cell) * nDim);

      for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
        cout << m1 + 1 << ". Eigen value  = " << setprecision(12) << eigenValues[m1] << endl;
      }
      cout << endl;
      if(m_printEigenValues) {
        ofstream eigenOutput;
        eigenOutput.open("eigenValues.txt", ofstream::app);
        for(MInt m1 = 0; m1 < a_noNodes(cell) * nDim; m1++) {
          eigenOutput << m1 + 1 << ". Eigen value  = " << setprecision(12) << eigenValues[m1] << endl;
        }
        eigenOutput << endl;
      }
    }

    computeAssembledSystemMatrix(Ke, cell);
  }

  mDeallocate(Ke_final);
  mDeallocate(eigenValues);
}

template <MInt nDim>
void FcSolver<nDim>::updateDisplacements() {
  TRACE();

  for(MInt i = 0; i < m_totalNoNodes; i++) {
    for(MInt d = 0; d < nDim; d++) {
      m_totalNodalDisplacements[i * nDim + d] += m_nodalDisplacements[i * nDim + d];
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::resetDisplacements(const MInt mode) {
  TRACE();

  switch(mode) {
    case -1: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_totalNodalDisplacements[i * nDim + d] = F0;
          m_nodalDisplacements[i * nDim + d] = F0;
        }
      }
      break;
    }
    case 0: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_nodalDisplacements[i * nDim + d] = F0;
        }
      }
      break;
    }
    default: {
      TERMM(1, "Not implemented for this method (resetDisplacements).");
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::resetLoadVector(const MInt mode) {
  TRACE();

  switch(mode) {
    // DEfault: Reset all
    case -1: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_externalLoadVector[i * nDim + d] = F0;
          m_internalLoadVector[i * nDim + d] = F0;
          m_reactionForceVector[i * nDim + d] = F0;
          m_nodalLoadVector[i * nDim + d] = F0;
        }
      }
      break;
    }
    // Mode 0: resets all forces but the the total load vector
    case 0: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_externalLoadVector[i * nDim + d] = F0;
          m_internalLoadVector[i * nDim + d] = F0;
          m_reactionForceVector[i * nDim + d] = F0;
        }
      }
      break;
    }
    // Mode 1: resets only internal loads and reaction forces
    case 1: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_internalLoadVector[i * nDim + d] = F0;
          m_reactionForceVector[i * nDim + d] = F0;
        }
      }
      break;
    }
    // Mode 2 : resets only the external forces
    case 2: {
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        for(MInt d = 0; d < nDim; d++) {
          m_externalLoadVector[i * nDim + d] = F0;
        }
      }
      break;
    }
    default: {
      TERMM(1, "Not implemented for this method (resetLoadVector).");
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::updateLoadVector(const MFloat lambda) {
  TRACE();

  for(MInt i = 0; i < m_totalNoNodes; i++) {
    for(MInt d = 0; d < nDim; d++) {
      const MFloat extLoad = m_externalLoadVector[i * nDim + d] * lambda;
      const MFloat intLoad = m_internalLoadVector[i * nDim + d];
      const MFloat reacForce = m_reactionForceVector[i * nDim + d];
      m_nodalLoadVector[i * nDim + d] = extLoad - (intLoad - reacForce);
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::createCompressedMatrix() {
  TRACE();

  MInt length = m_finalGlobalMatrix.size();

  if(m_sysMatCompressed != nullptr) {
    cerr << "Deallocate m_sysMatCompressed" << endl;
    mDeallocate(m_sysMatCompressed);
  }
  if(m_compressionIndexSys != nullptr) {
    cerr << "Deallocate m_compressionIndexSys" << endl;
    mDeallocate(m_compressionIndexSys);
  }

  mAlloc(m_sysMatCompressed, length, "m_sysMatCompressed", F0, AT_);

  mAlloc(m_compressionIndexSys, length, 2, "m_compressionIndexSys", 0, AT_);

  MInt cnt = 0;
  auto iterMat = m_finalGlobalMatrix.begin();
  while(iterMat != m_finalGlobalMatrix.end()) {
    // TODO: Necessary??
    if(approx(iterMat->second, F0, m_eps)) {
      m_sysMatCompressed[cnt] = F0;
    } else {
      m_sysMatCompressed[cnt] = iterMat->second;
    }
    m_compressionIndexSys[cnt][0] = iterMat->first.first;
    m_compressionIndexSys[cnt][1] = iterMat->first.second;
    iterMat++;
    cnt++;
  }
}

/**
 * \brief Computes the Assembled system matrix by summing the element mass matrices
 *
 * \author Moritz Waldmann
 * \date 10.01.22
 *
 * Boundary conditions are applied first. Then, the system of equation is
 * solved using Eigen
 */
template <MInt nDim>
void FcSolver<nDim>::solveSystemOfEquations() {
  TRACE();

  if(m_testRun) {
    cout << "///////////////////////////////////////////////////////////////" << endl;
    cout << "///////////////////////////////////////////////////////////////" << endl;
    cout << endl;
    cout << "SOLUTION STEP STARTS:" << endl;
    cout << endl;
    cout << fixed << setprecision(12) << endl;
    cout << "FORCE VECTOR (sorted by the global point Id): " << endl;
    for(MInt dom = 0; dom < noDomains(); dom++) {
      MPI_Barrier(mpiComm(), AT_);
      if(domainId() != dom) continue;
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        cout << "Node Id " << i << " ";
        for(MInt d = 0; d < nDim; d++) {
          cout << m_nodalLoadVector[i * nDim + d] << " ";
        }
        cout << endl;
      }
    }
    cout << fixed << setprecision(5) << endl;
  }

  const MInt n = m_finalGlobalMatrix.size();
  MInt const m = (m_totalNoNodes * nDim);
  if(!m_solveSoEIteratively) {
    maia::math::solveSparseMatrix(m_sysMatCompressed, m_compressionIndexSys, n, m, m_nodalLoadVector,
                                  m_nodalDisplacements);
  } else {
    if(m_useEigen) {
      maia::math::solveSparseMatrixIterative(m_sysMatCompressed, m_compressionIndexSys, n, m, m_nodalLoadVector,
                                             m_nodalDisplacements);
    } else {
      solveMatrixIterativelyPreCon(m_sysMatCompressed, m_compressionIndexSys, n, m, m_nodalLoadVector,
                                   m_nodalDisplacements);
    }
  }

  // Update total displacements and calculate strains and stresses
  updateDisplacements();

  if(m_testRun) {
    MFloat sum = F0;
    cout << fixed << setprecision(12) << endl;
    cout << "DISPLACEMENT VECTOR (sorted by the global point Id): " << endl;
    for(MInt dom = 0; dom < noDomains(); dom++) {
      MPI_Barrier(mpiComm(), AT_);
      if(domainId() != dom) continue;
      for(MInt i = 0; i < m_totalNoNodes; i++) {
        cout << "Node Id " << i << " ";
        for(MInt d = 0; d < nDim; d++) {
          cout << m_totalNodalDisplacements[i * nDim + d] << " ";
        }
        cout << endl;
      }
    }
    cout << endl;

    for(MInt i = 0; i < m_totalNoNodes; i++) {
      MFloat multi = F0;
      for(MInt d = 0; d < nDim; d++) {
        multi += F1B2 * (m_totalNodalDisplacements[i * nDim + d] * m_nodalLoadVector[i * nDim + d]);
      }
      sum += multi;
    }
    if(noDomains() > 1) {
      MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "sum");
    }
    MFloat res = sqrt(fabs(m_analyticSolution + sum) / fabs(m_analyticSolution));
    cout << "SUM " << -sum << " " << res << endl;
    cout << fixed << setprecision(5) << endl;
  }
}

template <MInt nDim>
void FcSolver<nDim>::getCoordinatesOfNode(MInt node, MInt cell, MFloat* coordinates) {
  TRACE();

  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = (nDim == 2) ? 4 : 6;

  MFloat cellHalfLen = F1B2 * c_cellLengthAtCell(cell);

  if(node < noVertices) {
    for(MInt d = 0; d < nDim; d++) {
      MFloat coord = c_coordinate(cell, d);
      MFloat nodeCoord = coord + cellHalfLen * Fd::vertexPosition(node, d);
      coordinates[d] = nodeCoord;
    }
  } else if(node < noVertices + noEdges * a_pRfnmnt(cell)) {
    MInt edge = (MInt)floor((node - noVertices) / a_pRfnmnt(cell));
    MInt l = (node - noVertices) % a_pRfnmnt(cell);
    for(MInt d = 0; d < nDim; d++) {
      MFloat nodeCoord = F0;
      MFloat coord = c_coordinate(cell, d);
      if(approx(Fd::edgePosition(edge, d), F0, 0.0000001)) {
        nodeCoord = coord + cellHalfLen * a_nodePosition(cell, l + 1);
      } else {
        nodeCoord = coord + cellHalfLen * Fd::edgePosition(edge, d);
      }
      coordinates[d] = nodeCoord;
    }
  } else if(node < noVertices + noEdges * a_pRfnmnt(cell) + noSurfaces * a_pRfnmnt(cell) * a_pRfnmnt(cell)) {
    MInt surface = (MInt)floor((node - noVertices - noEdges * a_pRfnmnt(cell)) / (a_pRfnmnt(cell) * a_pRfnmnt(cell)));
    MInt l[2] = {0};
    l[1] = (node - noVertices - noEdges * a_pRfnmnt(cell)) % (a_pRfnmnt(cell) * a_pRfnmnt(cell));
    while(l[1] >= a_pRfnmnt(cell)) {
      l[0]++;
      l[1] -= a_pRfnmnt(cell);
    }
    MInt counter = 0;
    for(MInt d = 0; d < nDim; d++) {
      MFloat nodeCoord = F0;
      MFloat coord = c_coordinate(cell, d);
      if(approx(Fd::surfacePosition(surface, d), F0, 0.0000001)) {
        nodeCoord = coord + cellHalfLen * a_nodePosition(cell, l[counter] + 1);
        counter++;
      } else {
        nodeCoord = coord + cellHalfLen * Fd::surfacePosition(surface, d);
      }
      coordinates[d] = nodeCoord;
    }
  } else {
    MInt l[3] = {0};
    l[2] = node - noVertices - noEdges * a_pRfnmnt(cell) - noSurfaces * a_pRfnmnt(cell) * a_pRfnmnt(cell);
    while(l[2] >= (a_pRfnmnt(cell) * a_pRfnmnt(cell))) {
      l[0]++;
      l[2] -= (a_pRfnmnt(cell) * a_pRfnmnt(cell));
    }
    while(l[2] >= a_pRfnmnt(cell)) {
      l[1]++;
      l[2] -= a_pRfnmnt(cell);
    }
    MInt counter = 0;
    for(MInt d = 0; d < nDim; d++) {
      MFloat coord = c_coordinate(cell, d);
      MFloat nodeCoord = coord + cellHalfLen * a_nodePosition(cell, l[counter] + 1);
      counter++;
      coordinates[d] = nodeCoord;
    }
  }
}

// The reordered BiCGStab method for distributed memory computer systems
// International Conference on Computational Science, ICCS 2010
// Procedia Computer Science 1 (2012) 213–218
template <MInt nDim>
void FcSolver<nDim>::solveMatrixIteratively(MFloat* A_coeff, MInt** pos, const MInt n) {
  TRACE();

  // 0. Allocate vectors required for the BiCGStab
  MFloatScratchSpace r0(m_totalNoNodes * nDim, "r0", AT_);
  r0.fill(F0);
  MFloatScratchSpace r(m_totalNoNodes * nDim, "r", AT_);
  r.fill(F0);
  MFloatScratchSpace nu(m_totalNoNodes * nDim, "nu", AT_);
  nu.fill(F0);
  MFloatScratchSpace p(m_totalNoNodes * nDim, "p", AT_);
  p.fill(F0);
  MFloatScratchSpace s(m_totalNoNodes * nDim, "s", AT_);
  s.fill(F0);
  MFloatScratchSpace t(m_totalNoNodes * nDim, "t", AT_);
  t.fill(F0);
  MFloatScratchSpace b_calc(m_totalNoNodes * nDim, "b_calc", AT_);
  b_calc.fill(F0);

  MFloat rho = F0;
  MFloat alpha = F0;
  MFloat omega = F0;
  MFloat beta = F0;
  MFloat b = F0;

  // 1. Initial guess of x0 is 0
  // 2. Calculate the dotProduct of r0 and r0
  // 3. Set p0 and r equal to r0
  for(MInt i = 0; i < m_totalNoNodes * nDim; i++) {
    r0[i] = m_nodalLoadVector[i];
    r[i] = r0[i];
    p[i] = m_nodalLoadVector[i];
  }
  for(MInt i = 0; i < m_noInternalNodes * nDim; i++) {
    rho += (r0[i] * r0[i]);
    b += (m_nodalLoadVector[i] * m_nodalLoadVector[i]);
  }
  cout << setprecision(12) << endl;
  MPI_Allreduce(MPI_IN_PLACE, &rho, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "rho");
  MPI_Allreduce(MPI_IN_PLACE, &b, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "b");

  for(MInt j = 0; j < m_maxNoIterations; j++) {
    // 4. Recalculate the vector nu with nu = A * p_j
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      nu[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      nu[pos[k][0]] += A_coeff[k] * p[pos[k][1]];
    }

    // 5. Set alpha_j to the ratio of rho and the dotproduct of nu and r0
    MFloat dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (nu[k] * r0[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    alpha = rho / dotProduct;

    // 6. Determine the vector s with s_j = r_j - alpha_j * nu
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      s[k] = r[k] - alpha * nu[k];
    }

    // 7. Recalculate the vector t with t = A * s_j
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      t[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      t[pos[k][0]] += A_coeff[k] * s[pos[k][1]];
    }

    // 8. Set omega_j to the ratio of the dotproduct of t and s_j and the dotproduct of t and t
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (t[k] * s[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    omega = dotProduct;
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (t[k] * t[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    omega /= dotProduct;

    // 9. Calculate the new result vector x_(j+1) with
    //   x_(j+1) = x_j + alpha_j * p_j + omega_j * s_j and check the convergence
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      m_nodalDisplacements[k] += alpha * p[k] + omega * s[k];
    }

    // 9b. Check the convergence: c = ||b - A * x|| / ||b||
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      b_calc[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      b_calc[pos[k][0]] += A_coeff[k] * m_nodalDisplacements[pos[k][1]];
    }
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += ((m_nodalLoadVector[k] - b_calc[k]) * (m_nodalLoadVector[k] - b_calc[k]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    MFloat convergence = dotProduct / b;
    cout << "Convergence of iteration " << j << " is " << convergence << endl;
    if(convergence < m_eps) {
      cout << "REACHED CONVERGENCY" << endl;
      break;
    }

    // 10. Calculate the vector r_(j+1) for the next iteration with r_(j+1) = s_j - omega_j * t
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      r[k] = s[k] - omega * t[k];
    }

    // 11a. We need to store 1.0/rho_j for later
    beta = F1 / rho;

    // 11b. Set rho_(j+1) to the dotproduct of r_(j+1) and r_0
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (r[k] * r0[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    rho = dotProduct;

    // 12. Set beta_j to beta_j = rho_(j+1)/rho_j * alpha_j/omega_j
    beta *= rho * alpha / omega;

    // 13. Calculate p_(j+1) for the next iteration with p_(j+1) = r_(j+1) + beta * (p_j - omega * nu)
    for(MInt k = 0; k < m_totalNoNodes * nDim; k++) {
      p[k] = r[k] + beta * (p[k] - omega * nu[k]);
    }
  }
  cout << setprecision(5) << endl;
}

// The reordered BiCGStab method for distributed memory computer systems
// International Conference on Computational Science, ICCS 2010
// Procedia Computer Science 1 (2012) 213–218
template <MInt nDim>
void FcSolver<nDim>::solveMatrixIterativelyPreCon(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b,
                                                  MFloat* x) {
  TRACE();

  // 0. Allocate vectors required for the BiCGStab
  MFloatScratchSpace r0(m, "r0", AT_);
  r0.fill(F0);
  MFloatScratchSpace r(m, "r", AT_);
  r.fill(F0);
  MFloatScratchSpace nu(m, "nu", AT_);
  nu.fill(F0);
  MFloatScratchSpace nuHat(m, "nuHat", AT_);
  nuHat.fill(F0);
  MFloatScratchSpace p(m, "p", AT_);
  p.fill(F0);
  MFloatScratchSpace s(m, "s", AT_);
  s.fill(F0);
  MFloatScratchSpace t(m, "t", AT_);
  t.fill(F0);
  MFloatScratchSpace tHat(m, "tHat", AT_);
  tHat.fill(F0);
  MFloatScratchSpace z(m, "z", AT_);
  z.fill(F0);
  MFloatScratchSpace b_calc(m, "b_calc", AT_);
  b_calc.fill(F0);

  MFloatScratchSpace bestSolution(m, "bestSolution", AT_);
  bestSolution.fill(F0);

  MFloatScratchSpace M_inv(m, "M_inv", AT_);
  M_inv.fill(F0);

  MFloat rho = F0;
  MFloat alpha = F0;
  MFloat omega = F0;
  MFloat beta = F0;
  MFloat theta = F0;
  MFloat phi = F0;
  MFloat psi = F0;
  MFloat b_mag = F0;

  MFloat minConvergence = std::numeric_limits<MFloat>::max();
  MFloat currentConvergence = F0;
  MFloat previousConvergence = std::numeric_limits<MFloat>::max();
  MInt j = 0;

  exchangeVector(b);

  // 0. Initial guess of x0 is 0
  // 1. Calculate the dotProduct of r0 and r0
  // 2. Set r equal to r0
  for(MInt i = 0; i < m; i++) {
    r0[i] = b[i];
    r[i] = r0[i];
  }

  for(MInt i = 0; i < m_noInternalNodes * nDim; i++) {
    rho += (r0[i] * r0[i]);
    b_mag += (b[i] * b[i]);
  }

  // 3. M_inv is the inverse of diag(A) and we set z to M_inv * r0 and nuHat to z
  for(MInt k = 0; k < n; k++) {
    if(pos[k][0] == pos[k][1]) {
      M_inv[pos[k][0]] = F1 / A_coeff[k];
    }
  }
  exchangeVector(M_inv.data());

  for(MInt i = 0; i < m; i++) {
    z[i] = M_inv[i] * r0[i];
    nuHat[i] = z[i];
  }

  cout << setprecision(12) << endl;
  MPI_Allreduce(MPI_IN_PLACE, &rho, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "rho");
  MPI_Allreduce(MPI_IN_PLACE, &b_mag, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "b_mag");

  if(b_mag < m_eps) return;

  while(j < m_maxNoIterations
        || ((currentConvergence > minConvergence) && !approx(currentConvergence, previousConvergence, m_eps))) {
    // 4. Recalculate the vector nu with nu = A * nuHat
    for(MInt k = 0; k < m; k++) {
      nu[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      nu[pos[k][0]] += A_coeff[k] * nuHat[pos[k][1]];
    }
    exchangeVector(nu.data());

    // 5. Set alpha_j to the ratio of rho_j and the dotproduct of nu and r0
    MFloat dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (nu[k] * r0[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    alpha = rho / dotProduct;

    // 6. Determine the vector s with s = M_inv * nu
    for(MInt k = 0; k < m; k++) {
      s[k] = M_inv[k] * nu[k];
    }

    // 7. Recalculate the vector t with t = A * tHat and tHat = z - alpha_j * s
    for(MInt k = 0; k < m; k++) {
      t[k] = F0;
      tHat[k] = z[k] - alpha * s[k];
    }
    for(MInt k = 0; k < n; k++) {
      t[pos[k][0]] += (A_coeff[k] * tHat[pos[k][1]]);
    }
    exchangeVector(t.data());

    // 8. Calculate the new result vector x_(j+1) with
    //   x_(j+1) = x_j + alpha_j * nuHat and check the convergence
    for(MInt k = 0; k < m; k++) {
      x[k] += alpha * nuHat[k];
    }

    // 8b. Check the convergence: c = ||b - A * x|| / ||b||
    for(MInt k = 0; k < m; k++) {
      b_calc[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      b_calc[pos[k][0]] += (A_coeff[k] * x[pos[k][1]]);
    }
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += ((b[k] - b_calc[k]) * (b[k] - b_calc[k]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    MFloat convergence = dotProduct / b_mag;
    previousConvergence = currentConvergence;
    currentConvergence = convergence;
    if(convergence < minConvergence) {
      minConvergence = convergence;
      for(MInt k = 0; k < m; k++) {
        bestSolution[k] = x[k];
      }
    }
    if(convergence < m_eps) {
      break;
    }

    // 9. Calculate the vector r_(j+1) for the next iteration with r_(j+1) = r_j - alpha_j * nu
    for(MInt k = 0; k < m; k++) {
      r[k] = r[k] - alpha * nu[k];
    }

    // 10a. Set theta_j to the dotproduct of t and r_(j+1)
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (t[k] * r[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    theta = dotProduct;

    // 10b. Set phi_j to the dotproduct of t and t
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (t[k] * t[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    phi = dotProduct;

    // 10c. Set psi_j to the dotproduct of t and r0
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += (t[k] * r0[k]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    psi = dotProduct;

    // 11. Recalculate z with M_inv * t
    for(MInt k = 0; k < m; k++) {
      z[k] = M_inv[k] * t[k];
    }

    // 12. Set omega_j to the ratio of theta_j and phi_j
    omega = theta / phi;

    // 13. Calculate the new result vector x_(j+1) with
    //   x_(j+1) = x_(j+1) + omega_j * tHat and check the convergence
    for(MInt k = 0; k < m; k++) {
      x[k] += omega * tHat[k];
    }

    // 13b. Check the convergence: c = ||b - A * x|| / ||b||
    for(MInt k = 0; k < m; k++) {
      b_calc[k] = F0;
    }
    for(MInt k = 0; k < n; k++) {
      b_calc[pos[k][0]] += A_coeff[k] * x[pos[k][1]];
    }
    dotProduct = F0;
    for(MInt k = 0; k < m_noInternalNodes * nDim; k++) {
      dotProduct += ((b[k] - b_calc[k]) * (b[k] - b_calc[k]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &dotProduct, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dotProduct");
    convergence = dotProduct / b_mag;
    previousConvergence = currentConvergence;
    currentConvergence = convergence;
    if(convergence < minConvergence) {
      minConvergence = convergence;
      for(MInt k = 0; k < m; k++) {
        bestSolution[k] = x[k];
      }
    }
    if(convergence < m_eps) {
      break;
    }

    // 14. Recalculate r_(j+1) with r_(j+1) = r_(j+1) - omega * t
    for(MInt k = 0; k < m; k++) {
      r[k] = r[k] - omega * t[k];
    }

    // 15a. We need to store rho in beta
    beta = F1 / rho;

    // 15b. Set rho_(j+1) to -omega_j * psi_j
    rho = -omega * psi;

    // 14. Recalculate z with z = tHat - omega_j * z
    for(MInt k = 0; k < m; k++) {
      z[k] = tHat[k] - omega * z[k];
    }

    // 17. Set beta_j to beta_j = rho_(j+1)/rho_j * alpha_j/omega_j
    beta *= rho * alpha / omega;

    // 18. Recalculate nuHat with nuHat = z + beta_j * (nuHat - omega * s)
    for(MInt k = 0; k < m; k++) {
      nuHat[k] = z[k] + beta * (nuHat[k] - omega * s[k]);
    }
    j++;
    cout << j << " " << minConvergence << " " << currentConvergence << " " << previousConvergence << endl;
  }

  for(MInt k = 0; k < m; k++) {
    x[k] = bestSolution[k];
  }
  cout << "Solved SOE with " << j << " iterations and min convergence " << minConvergence << " and current convergence "
       << currentConvergence << endl;
  m_log << "Solved SOE with " << j << " iterations and min convergence " << minConvergence
        << " and current convergence " << currentConvergence << endl;
  cout << setprecision(5) << endl;
}

/** \brief Send the global nodeIds from one domain to all neighbor domains
 *
 * \author Moritz Waldmann
 * \date  20.03.22
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::exchangeVector(MFloat* vector) {
  TRACE();

  if(noNeighborDomains() < 1) return;

  // Create arrays for the data offsets and data count
  MIntScratchSpace windowDataPerRank(noNeighborDomains(), "windowDataPerRank", AT_);
  MIntScratchSpace haloDataPerRank(noNeighborDomains(), "haloDataPerRank", AT_);
  MIntScratchSpace windowDataOffsets(noNeighborDomains() + 1, "windowDataOffsets", AT_);
  MIntScratchSpace haloDataOffsets(noNeighborDomains() + 1, "haloDataOffsets", AT_);
  windowDataPerRank.fill(0);
  haloDataPerRank.fill(0);
  windowDataOffsets.fill(0);
  haloDataOffsets.fill(0);

  // Determine the data count and the data offsets
  MInt totalCounterWindow = 0;
  MInt totalCounterHalo = 0;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt dataCounterWindow = 0;
    MInt dataCounterHalo = 0;
    for(MInt w = 0; w < noWindowCells(n); w++) {
      if(!a_isActive(windowCell(n, w))) continue;
      for(MInt node = 0; node < a_noNodes(windowCell(n, w)); node++) {
        for(MInt dim = 0; dim < nDim; dim++) {
          dataCounterWindow++;
        }
      }
    }
    for(MInt h = 0; h < noHaloCells(n); h++) {
      if(!a_isActive(haloCell(n, h))) continue;
      for(MInt node = 0; node < a_noNodes(haloCell(n, h)); node++) {
        for(MInt dim = 0; dim < nDim; dim++) {
          dataCounterHalo++;
        }
      }
    }
    windowDataPerRank(n) = dataCounterWindow;
    haloDataPerRank(n) = dataCounterHalo;
    totalCounterWindow += dataCounterWindow;
    totalCounterHalo += dataCounterHalo;
    windowDataOffsets(n + 1) = totalCounterWindow;
    haloDataOffsets(n + 1) = totalCounterHalo;
  }

  // Allocate the buffers
  MFloat* sendBuffers{};
  MFloat* receiveBuffers{};
  MPI_Request* mpi_request{};
  mAlloc(sendBuffers, totalCounterWindow, "sendBuffers", F0, AT_);
  mAlloc(receiveBuffers, totalCounterHalo, "receiveBuffers", F0, AT_);
  mAlloc(mpi_request, noNeighborDomains(), "mpi_request", AT_);

  // Gather data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt cnt = 0;
    for(MInt j = 0; j < noWindowCells(n); j++) {
      if(!a_isActive(windowCell(n, j))) continue;
      for(MInt node = 0; node < a_noNodes(windowCell(n, j)); node++) {
        MInt nodeIdLocal = a_nodeIdsLocal(windowCell(n, j), node);
        for(MInt dim = 0; dim < nDim; dim++) {
          sendBuffers[windowDataOffsets(n) + cnt] = vector[nodeIdLocal * nDim + dim];
          cnt++;
        }
      }
    }
  }

  // Send data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Issend(&sendBuffers[windowDataOffsets(n)], windowDataPerRank(n), MPI_DOUBLE, neighborDomain(n), 0, mpiComm(),
               &mpi_request[n], AT_, "sendBuffers[n]");
  }

  // Receive data
  MPI_Status status;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Recv(&receiveBuffers[haloDataOffsets(n)], haloDataPerRank(n), MPI_DOUBLE, neighborDomain(n), 0, mpiComm(),
             &status, AT_, "receiveBuffers[n]");
  }
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Wait(&mpi_request[n], &status, AT_);
  }

  // Scatter data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt cnt = 0;
    for(MInt j = 0; j < noHaloCells(n); j++) {
      if(!a_isActive(haloCell(n, j))) continue;
      for(MInt node = 0; node < a_noNodes(haloCell(n, j)); node++) {
        MInt nodeIdLocal = a_nodeIdsLocal(haloCell(n, j), node);
        for(MInt dim = 0; dim < nDim; dim++) {
          vector[nodeIdLocal * nDim + dim] = receiveBuffers[haloDataOffsets(n) + cnt];
          cnt++;
        }
      }
    }
  }

  mDeallocate(sendBuffers);
  mDeallocate(receiveBuffers);
  mDeallocate(mpi_request);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// BELOW THIS LINE YOU CAN FIND ONLY DEBUG FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////

/** \brief Transformes the coordinate from the reference element to the global
 * coordinate system.
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 * The following is done in this function:
 *
 *  - calculates the lagrangian polynom for the cell
 *  - calculates the global coordinate of each lobatto point
 *  - calculates the global coordinate of point z of the reference
 *    element by summing the product of lagrangian polynome and
 *    global coordinates for each lobatte point
 *
 * The function should only be used for debug purpose (to check if the function
 * getDisplacementInterpolationMatrix works correctly.
 *
 **/
template <MInt nDim>
void FcSolver<nDim>::interpolateGeometryDebug(const MInt pCellId, const MFloat* z, MFloat* x) {
  TRACE();

  // p holds the number of lagrangian points in one direction without the corner points
  const MInt p = a_pRfnmnt(pCellId);

  // calculation of the lagrangian polynom
  MFloatScratchSpace L_coef(nDim, a_noNodes(pCellId) * nDim, AT_, "L_coef");
  L_coef.fill(F1);
  getDisplacementInterpolationMatrix(pCellId, z, L_coef);

  if(m_testRun && m_polyDeg < 1) getDisplacementInterpolationMatrixDebug(pCellId, z, L_coef);

  // Set the arrays containing the order of the surfaces, edges, and vertices
  // for the 2D or 3D case
  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = (nDim == 2) ? 4 : 6;

  // Calculate the global coordinate for each lagrangian point
  MFloatScratchSpace X(nDim, p + 2, AT_, "X");

  MFloat coord[nDim] = {F0};
  MFloat cellHalfLength = F1B2 * c_cellLengthAtCell(pCellId);

  for(MInt d = 0; d < nDim; d++)
    coord[d] = c_coordinate(pCellId, d);

  for(MInt d = 0; d < nDim; d++) {
    X(d, 0) = coord[d] - cellHalfLength;
    X(d, p + 1) = coord[d] + cellHalfLength;
  }

  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 1; i < p + 1; i++) {
      X(d, i) = coord[d] + Fd::nodePosLobattoPoints(p, i - 1) * cellHalfLength;
    }
  }

  // Loop over all points of the current cell
  // depending on the point number, the point is a
  // vertex, edge point, surface point, or inner point
  // The location of the point is calculated.
  // The product of the lagrangian polynom and the global
  // coordinate of the current point is added for each point
  for(MInt node = 0; node < a_noNodes(pCellId); node++) {
    if(node < noVertices) {
      for(MInt d = 0; d < nDim; d++) {
        MInt pos = -1;
        if(Fd::vertexPosition(node, d) < F0) {
          pos = 0;
        } else {
          pos = p + 1;
        }
        MFloat L = L_coef(d, node * nDim + d);
        x[d] += L * X(d, pos);
      }
    } else if(node < noVertices + noEdges * p) {
      MInt edge = (MInt)floor((node - noVertices) / p);
      MInt l = node - noVertices - edge * p;
      for(MInt d = 0; d < nDim; d++) {
        MInt pos = -1;
        if(Fd::edgePosition(edge, d) < F0) {
          pos = 0;
        } else if(Fd::edgePosition(edge, d) > F0) {
          pos = p + 1;
        } else {
          pos = l + 1;
        }
        MFloat L = L_coef(d, node * nDim + d);
        x[d] += L * X(d, pos);
      }
    } else if(node < noVertices + noEdges * p + noSurfaces * p * p) {
      MInt surface = (MInt)floor((node - noVertices - noEdges * p) / (p * p));
      MInt l[2] = {0};
      l[1] = (node - noVertices - noEdges * p) % (p * p);
      while(l[1] >= p) {
        l[0]++;
        l[1] -= p;
      }
      MInt counter = 0;
      for(MInt d = 0; d < nDim; d++) {
        MInt pos = -1;
        if(Fd::surfacePosition(surface, d) < F0) {
          pos = 0;
        } else if(Fd::surfacePosition(surface, d) > F0) {
          pos = p + 1;
        } else {
          pos = l[counter] + 1;
          counter++;
        }
        MFloat L = L_coef(d, node * nDim + d);
        x[d] += L * X(d, pos);
      }
    } else {
      MInt l[3] = {0};
      l[2] = node - noVertices - noEdges * p - noSurfaces * p * p;
      while(l[2] >= (p * p)) {
        l[0]++;
        l[2] -= (p * p);
      }
      while(l[2] >= p) {
        l[1]++;
        l[2] -= p;
      }
      MInt counter = 0;
      for(MInt d = 0; d < nDim; d++) {
        MInt pos = l[counter] + 1;
        counter++;
        MFloat L = L_coef(d, node * nDim + d);
        x[d] += L * X(d, pos);
      }
    }
  }
}

/** \brief fills the coordinates of each point in the displacement
 *  vector to test the calculation of the cell displacement for example
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 * The following is done in this function:
 *
 *
 * The function should only be used for debug purpose.
 **/
template <MInt nDim>
void FcSolver<nDim>::fillDisplacementVectorWithCoords() {
  TRACE();

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;
    for(MInt node = 0; node < a_noNodes(cell); node++) {
      MInt nodeId = a_nodeIdsLocal(cell, node);
      for(MInt d = 0; d < nDim; d++) {
        m_nodalDisplacements[nodeId * nDim + d] = 10000.0;
      }
    }
  }

  for(MInt cell = 0; cell < a_noCells(); cell++) {
    if(!a_isActive(cell)) continue;

    for(MInt node = 0; node < a_noNodes(cell); node++) {
      MFloat nodalCoord[nDim];
      MInt nodeId = a_nodeIdsLocal(cell, node);
      getCoordinatesOfNode(node, cell, nodalCoord);
      for(MInt d = 0; d < nDim; d++) {
        if(m_nodalDisplacements[nodeId * nDim + d] < 10000.0) {
          if(!approx(m_nodalDisplacements[nodeId * nDim + d], nodalCoord[d], m_eps)) {
            stringstream errorMessage;
            errorMessage << "ERROR in node coordinates" << m_nodalDisplacements[nodeId * nDim + d] << " "
                         << nodalCoord[d] << " " << c_coordinate(cell, d) << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        } else {
          m_nodalDisplacements[nodeId * nDim + d] = nodalCoord[d];
        }
      }
    }
  }
  updateDisplacements();
}

/** \brief Calculates the jacobian matrix
 *
 * \author Moritz Waldmann
 * \date 10.03.21
 *
 * Should be
 *
 *    0.5 * cellLength 0.0              0.0
 *    0.0              0.5 * cellLength 0.0
 *    0.0              0.0              0.5 * cellLength
 **/
template <MInt nDim>
void FcSolver<nDim>::lagrangianPointsJacobi(const MInt pCellId, const MFloat* z, MFloatScratchSpace& x_prime) {
  TRACE();

  const MInt p = a_pRfnmnt(pCellId);
  MFloatScratchSpace X(nDim, p + 2, AT_, "X");

  MFloatScratchSpace L_prime(a_noNodes(pCellId), nDim, AT_, "L_prime");
  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d = 0; d < nDim; d++) {
      L_prime(n, d) = F1;
    }
  }
  getDerivativeOfDisplacementInterpolationMatrix(pCellId, z, L_prime);

  MInt noVertices = IPOW2(nDim);
  MInt noEdges = (nDim == 2) ? 0 : 12;
  MInt noSurfaces = (nDim == 2) ? 4 : 6;

  MFloat coord[nDim] = {F0};
  MFloat cellHalfLength = F1B2 * c_cellLengthAtCell(pCellId);

  for(MInt d = 0; d < nDim; d++)
    coord[d] = c_coordinate(pCellId, d);

  for(MInt d = 0; d < nDim; d++) {
    X(d, 0) = coord[d] - cellHalfLength;
    X(d, p + 1) = coord[d] + cellHalfLength;
  }

  for(MInt d = 0; d < nDim; d++) {
    for(MInt i = 1; i < p + 1; i++) {
      X(d, i) = coord[d] + Fd::nodePosLobattoPoints(p, i - 1) * cellHalfLength;
    }
  }

  for(MInt node = 0; node < a_noNodes(pCellId); node++) {
    if(node < noVertices) {
      for(MInt dir = 0; dir < nDim; dir++) {
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          if(Fd::vertexPosition(node, d) < F0) {
            pos = 0;
          } else {
            pos = p + 1;
          }
          x_prime(d, dir) += L_prime(node, dir) * X(d, pos);
        }
      }
    } else if(node < noVertices + noEdges * p) {
      MInt edge = (MInt)floor((node - noVertices) / p);
      MInt l = node - noVertices - edge * p;
      for(MInt dir = 0; dir < nDim; dir++) {
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          if(Fd::edgePosition(edge, d) < F0) {
            pos = 0;
          } else if(Fd::edgePosition(edge, d) > F0) {
            pos = p + 1;
          } else {
            pos = l + 1;
          }
          x_prime(d, dir) += L_prime(node, dir) * X(d, pos);
        }
      }
    } else if(node < noVertices + noEdges * p + noSurfaces * p * p) {
      MInt surface = (MInt)floor((node - noVertices - noEdges * p) / (p * p));
      MInt l[2] = {0};
      l[1] = (node - noVertices - noEdges * p) % (p * p);
      while(l[1] >= p) {
        l[0]++;
        l[1] -= p;
      }
      for(MInt dir = 0; dir < nDim; dir++) {
        MInt counter = 0;
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = -1;
          if(Fd::surfacePosition(surface, d) < F0) {
            pos = 0;
          } else if(Fd::surfacePosition(surface, d) > F0) {
            pos = p + 1;
          } else {
            pos = l[counter] + 1;
            counter++;
          }
          x_prime(d, dir) += L_prime(node, dir) * X(d, pos);
        }
      }
    } else {
      MInt l[3] = {0};
      l[2] = node - noVertices - noEdges * p - noSurfaces * p * p;
      while(l[2] >= (p * p)) {
        l[0]++;
        l[2] -= (p * p);
      }
      while(l[2] >= p) {
        l[1]++;
        l[2] -= p;
      }
      for(MInt dir = 0; dir < nDim; dir++) {
        MInt counter = 0;
        for(MInt d = 0; d < nDim; d++) {
          MInt pos = l[counter] + 1;
          counter++;
          x_prime(d, dir) += L_prime(node, dir) * X(d, pos);
        }
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::getDisplacementInterpolationMatrixDebug(const MInt pCellId, const MFloat* z,
                                                             MFloatScratchSpace& L_coef_lagrange) {
  TRACE();

  MFloatScratchSpace L_coef(a_noNodes(pCellId), AT_, "L_coef");

  if(a_noNodes(pCellId) != IPOW2(nDim)) TERMM(1, "Not implemented for solver.");
  if(a_pRfnmnt(pCellId) != 0) TERMM(1, "Not implemented for solver.");

  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt i = 0; i < IPOW2(nDim); i++) {
      L_coef(i) = FFPOW2(nDim);
      for(MInt d = 0; d < nDim; d++) {
        L_coef(i) *= (F1 + Fd::vertexPosition(i, d) * z[d]);
      }
      for(MInt d = 0; d < nDim; d++) {
        MFloat L = L_coef_lagrange(dim, i * nDim + d);
        if(d == dim) {
          if(!approx(L_coef(i), L, 1e-14)) {
            stringstream errorMessage;
            errorMessage << "ERROR!! Displacement Interpolation for x " << i << " is: " << L_coef(i) << " and " << L
                         << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        } else {
          if(!approx(F0, L, 1e-14)) {
            stringstream errorMessage;
            errorMessage << "ERROR!! Displacement Interpolation for x " << i << " is: " << F0 << " and " << L << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        }
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::getJacobiMatrixDebug(const MInt pCellId, const MFloat* z, MFloatScratchSpace& jacobi_lagrange) {
  TRACE();

  MFloatScratchSpace jacobi(nDim * nDim, AT_, "jacobi");

  if(a_noNodes(pCellId) != IPOW2(nDim)) TERMM(1, "Not implemented for solver.");
  if(a_pRfnmnt(pCellId) != 0) TERMM(1, "Not implemented for solver.");

  MFloat halfLength = c_cellLengthAtCell(pCellId) * F1B2;
  for(MInt d1 = 0; d1 < nDim; d1++) {   // Derivative direction
    for(MInt d2 = 0; d2 < nDim; d2++) { // Coordinate direction
      jacobi(d1 * nDim + d2) = F0;
      for(MInt n = 0; n < IPOW2(nDim); n++) { // Number of vertices
        MFloat cCoord = c_coordinate(pCellId, d1) + Fd::vertexPosition(n, d1) * halfLength;
        MFloat jac = cCoord * FFPOW2(nDim);
        for(MInt d = 0; d < nDim; d++) {
          if(d == d2) {
            jac *= Fd::vertexPosition(n, d);
          } else {
            jac *= (F1 + Fd::vertexPosition(n, d) * z[d]);
          }
        }
        jacobi(d1 * nDim + d2) += jac;
      }
    }
  }
  for(MInt d1 = 0; d1 < nDim; d1++) {   // Derivative direction
    for(MInt d2 = 0; d2 < nDim; d2++) { // Coordinate direction
      if(!approx(jacobi(d1 * nDim + d2), F0, 1e-10)) {
        MFloat inv = F1 / jacobi(d1 * nDim + d2);
        if(!approx(inv, jacobi_lagrange(d1 * nDim + d2), 1e-10)) {
          stringstream errorMessage;
          errorMessage << "ERROR!! Jacobi for x" << d1 << " and z" << d2 << " is: " << jacobi(d1 * nDim + d2) << " and "
                       << F1 / jacobi_lagrange(d1 * nDim + d2) << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    }
  }
}

template <MInt nDim>
void FcSolver<nDim>::getStrainInterpolationMatrixDebug(const MInt pCellId, const MFloat* z,
                                                       MFloatScratchSpace& strain_lagrange) {
  TRACE();

  MFloatScratchSpace strain(m_noStrains, nDim * a_noNodes(pCellId), AT_, "strain");
  strain.fill(F0);

  MFloatScratchSpace L_prime(a_noNodes(pCellId), nDim, AT_, "L_prime");
  MFloatScratchSpace L_prime_mul_inv_jac(a_noNodes(pCellId), nDim, AT_, "L_prime_mul_inv_jac");

  MFloatScratchSpace jacobi(nDim * nDim, AT_, "jacobi");

  if(a_noNodes(pCellId) != IPOW2(nDim)) TERMM(1, "Not implemented for solver.");
  if(a_pRfnmnt(pCellId) != 0) TERMM(1, "Not implemented for solver.");

  MFloat halfLength = c_cellLengthAtCell(pCellId) * F1B2;
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      jacobi(d1 * nDim + d2) = F0;
      for(MInt n = 0; n < IPOW2(nDim); n++) {
        MFloat cCoord = c_coordinate(pCellId, d1) + Fd::vertexPosition(n, d1) * halfLength;
        MFloat jac = cCoord * FFPOW2(nDim);
        for(MInt d = 0; d < nDim; d++) {
          if(d == d2) {
            jac *= Fd::vertexPosition(n, d);
          } else {
            jac *= (F1 + Fd::vertexPosition(n, d) * z[d]);
          }
        }
        jacobi(d1 * nDim + d2) += jac;
      }
    }
  }
  maia::math::invert(jacobi.getPointer(), nDim, nDim);

  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      MFloat derivative = FFPOW2(nDim);
      for(MInt d2 = 0; d2 < nDim; d2++) {
        if(d2 == d1) {
          derivative *= Fd::vertexPosition(n, d2);
        } else {
          derivative *= (F1 + Fd::vertexPosition(n, d2) * z[d2]);
        }
      }
      L_prime(n, d1) = derivative;
    }
  }

  for(MInt n = 0; n < a_noNodes(pCellId); n++) {
    for(MInt d1 = 0; d1 < nDim; d1++) {
      L_prime_mul_inv_jac(n, d1) = F0;
      for(MInt d2 = 0; d2 < nDim; d2++) {
        L_prime_mul_inv_jac(n, d1) += (L_prime(n, d2) * jacobi(d2 * nDim + d1));
      }
    }
  }

  for(MInt s = 0; s < nDim; s++) {
    for(MInt n = 0; n < a_noNodes(pCellId); n++) {
      for(MInt d = 0; d < nDim; d++) {
        if(s == d) {
          strain(s, n * nDim + d) = L_prime_mul_inv_jac(n, d);
        }
      }
    }
  }
  for(MInt s = 0; s < m_noStrains - nDim; s++) {
    for(MInt n = 0; n < a_noNodes(pCellId); n++) {
      for(MInt d = 0; d < nDim; d++) {
        if(s == d) {
          if((d + 1) < nDim) {
            strain(s + nDim, n * nDim + d) = L_prime_mul_inv_jac(n, d + 1);
            strain(s + nDim, n * nDim + d + 1) = L_prime_mul_inv_jac(n, d);
          } else {
            strain(s + nDim, n * nDim + d) = L_prime_mul_inv_jac(n, 0);
            strain(s + nDim, n * nDim + 0) = L_prime_mul_inv_jac(n, d);
          }
        }
      }
    }
  }
  for(MInt s = 0; s < m_noStrains; s++) {
    for(MInt n = 0; n < a_noNodes(pCellId) * nDim; n++) {
      if(!approx(strain(s, n), strain_lagrange(s, n), m_eps)) {
        stringstream errorMessage;
        errorMessage << "ERROR!! Strain Interpolation for s " << s << " and x " << n << " is: " << strain(s, n)
                     << " and " << strain_lagrange(s, n) << endl;
        mTerm(1, AT_, errorMessage.str());
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
// ABOVE THIS LINE YOU CAN FIND ONLY DEBUG FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////


template class FcSolver<2>;
template class FcSolver<3>;
