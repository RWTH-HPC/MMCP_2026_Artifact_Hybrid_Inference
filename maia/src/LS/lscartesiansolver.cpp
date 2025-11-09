// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lscartesiansolver.h"
#include <algorithm>
#include <stack>
#include "COMM/mpioverride.h"
#include "COUPLER/lsfv.h"
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "globals.h"

using namespace std;
using namespace maia;
using namespace maia::ls;

#define MAX_NO_LS_BNDRY_CND 10
//#define LS_DEBUG

//----------------------------------------------------------------------------

template <MInt nDim>
LsCartesianSolver<nDim>::LsCartesianSolver(MInt solverId_, const MBool* propertiesGroups, GridProxy& gridProxy_,
                                           Geometry<nDim>& geometry_, const MPI_Comm comm)
  : maia::CartesianSolver<nDim, LsCartesianSolver<nDim>>(solverId_, gridProxy_, comm, true), m_geometry(&geometry_) {
  TRACE();

  // set default
  m_maxNoSets = 1;

  const MLong oldAllocatedBytes = allocatedBytes();

  // general levelSet activation!
  const MBool levelSet = propertiesGroups[LEVELSET];

  ASSERT(levelSet, "");

  if(grid().azimuthalPeriodicity() && nDim == 2) {
    mTerm(1, "Azimuthal periodicity is untested in 2D and it doesn't make sense in 2D!");
  }

  // leveset used for emerged moving boundaries
  m_levelSetLb = propertiesGroups[LEVELSET_LB];
  m_levelSetMb = propertiesGroups[LEVELSETMB] || m_levelSetLb;
  // levelset used for combustion
  m_combustion = propertiesGroups[COMBUSTION];
  m_LSSolver = propertiesGroups[LS_SOLVER];
  m_levelSetRans = propertiesGroups[LS_RANS];
  // leveset used for the tracking of free surfaces between to fluid phases
  m_freeSurface = (string2enum(solverMethod()) == MAIA_LEVELSET_SURFACE);

  m_time = F0;

  m_levelSetFv = false;
  if((levelSet || m_levelSetRans) && !m_levelSetMb && !m_combustion && !m_LSSolver && !m_freeSurface) {
    m_levelSetFv = true;
  }

  if(domainId() == 0) {
    cerr << "m_levelSetMb  : " << m_levelSetMb << endl;
    cerr << "m_combustion  : " << m_combustion << endl;
    cerr << "m_LSSolver     : " << m_LSSolver << endl;
    cerr << "m_levelSetRans: " << m_levelSetRans << endl;
    cerr << "m_levelSetFv  : " << m_levelSetFv << endl;
    cerr << "m_freeSurface: " << m_freeSurface << endl;
  }

  ASSERT(m_LSSolver || m_levelSetMb || m_combustion || m_levelSetRans || m_levelSetFv || m_freeSurface,
         "Unintentianal-ls-solver initialisation?!");

  readLevelSetProperties();
  allocateLevelSetMemory();

  if(isActive()) printAllocatedMemory(oldAllocatedBytes, "LsCartesianSolver", mpiComm());

  initializeTimers();
}


template <MInt nDim>
void LsCartesianSolver<nDim>::initSolver() {
  TRACE();

  ASSERT(m_LSSolver || m_combustion || m_levelSetMb || m_levelSetRans || m_levelSetFv || m_freeSurface,
         "No valid LS mode");

  // If solver inactive only mark all cells as halo cells
  // Inactive solvers can not have internal cells
  if(!isActive()) {
    this->setHaloCellsOnInactiveRanks();
    return;
  }

  grid().updateLeafCellExchange();

  // Set the restart time step correctly when using useNonSpecifiedRestartFile
  if(m_useNonSpecifiedRestartFile) {
    m_restartTimeStep = globalTimeStep;
  }

  // Init the exchange for azimuthal periodicity
  initAzimuthalExchange();

  if(m_levelSetMb && m_constructGField) return;

  if(m_restart) {
    restartLocalizedLevelSetCG();
  } else {
    initLocalizedLevelSetCG();
  }

  // NOTE: only possible after the cells are added to the collector!
  this->checkNoHaloLayers();

  computeGCellTimeStep();

  if(m_combustion) {
    determineSteadyFlameLength();
    return;
  }

  if(m_freeSurface) {
    return;
  }

  // if(m_initialRefinement) return;

  setGCellBndryProperty();

  if(!m_semiLagrange) {
    computeNormalVectors();
    computeCurvature();
  }

  finalizeLevelSetInitialization();

  buildMultipleLevelSet();

  checkHaloCells();

  if(!m_semiLagrange) {
    determinePropagationSpeed();
  }

  if(m_restart && m_GFieldFromSTLInitCheck) {
    writeRestartLevelSetFileCG(1, "restartLSGridCG_restart", "restartLSCG_restart");
  }
}


template <MInt nDim>
void LsCartesianSolver<nDim>::allocateLevelSetMemory() {
  TRACE();

  if(m_combustion) {
    m_maxFlameFrontPosition = (MFloat*)nullptr;
    m_minFlameFrontPosition = (MFloat*)nullptr;
    m_meanFlameFrontPosition = (MFloat*)nullptr;

    mAlloc(m_maxFlameFrontPosition, nDim, "m_maxFlameFrontPosition", -1000.0, AT_);
    mAlloc(m_minFlameFrontPosition, nDim, "m_minFlameFrontPosition", 1000.0, AT_);
    mAlloc(m_meanFlameFrontPosition, nDim, "m_meanFlameFrontPosition", F0, AT_);
  }

  m_cells.setMaxNoSets(m_maxNoSets);

  // set the ls-Collector type to reduce ls-Cell memory!
  if(m_combustion || m_freeSurface) {
    // combustion Type
    m_lsCollectorMode = 0;
  } else if(m_semiLagrange) {
    // semi-Lagrange Type with reinitialisation
    m_lsCollectorMode = 2;
    if(!m_guaranteeReinit && m_STLReinitMode == 2) {
      // semi-Lagrange without reinitialisation
      m_lsCollectorMode = 1;
    }
  } else {
    // just the ls-Solver or lsFv
    // TODO labels:LS,FV reduce memory further for applications which don't need all information!
    ASSERT(m_LSSolver || m_levelSetFv || m_levelSetRans, "");
    m_lsCollectorMode = 3;
  }

  m_cells.setLsCollectorType(m_lsCollectorMode);

  m_noBodiesToCompute = 0;

  if(m_semiLagrange) {
    if(m_LsRotate) {
      allocateRotatingLs();
    }

    mAlloc(m_semiLagrange_xShift_ref, m_noEmbeddedBodies * nDim, "m_semiLagrange_xShift_ref", F0, AT_);
  }

  m_cells.setGapClosing(m_closeGaps);
  m_cells.setMaxBodiesToCompute(m_noBodiesToCompute);
  m_cells.setReconstructOldG(m_reconstructOldG);
  m_cells.setRotatingLs(m_LsRotate);
  m_cells.setReinit(m_guaranteeReinit || m_STLReinitMode != 2);

  m_cells.reset(m_maxNoCells);

  if(m_LsRotate) {
    m_cells.fillContainingCell();
    if(!m_reconstructOldG) m_cells.fillContainingDomain();
  }

  // multiple level-set functions
  m_noGapCells = 0;
  m_noOldGapCells = 0;

  mAlloc(m_bandCells, m_maxNoSets, "m_bandCells", AT_);
  mAlloc(m_internalBandCells, m_maxNoSets, "m_internalBandCells", AT_);
  mAlloc(m_bandBndryCells, m_maxNoSets, "m_bandBndryCells", AT_);
  mAlloc(m_G0Cells, m_maxNoSets, "m_G0Cells", AT_);
  for(MInt i = 0; i < m_maxNoSets; i++) {
    m_bandCells[i].clear();
    m_internalBandCells[i].clear();
    m_bandBndryCells[i].clear();
    m_G0Cells[i].clear();
  }
  mAlloc(m_bandLayer, (m_gShadowWidth + 1) * m_maxNoSets, "m_bandLayer", 0, AT_);
  mAlloc(m_internalBandLayer, (m_gShadowWidth + 1) * m_maxNoSets, "m_internalBandLayer", 0, AT_);

  if(!m_semiLagrange) {
    mAlloc(m_gBndryCells, m_maxNoSets, "m_gBndryCells", AT_);
    for(MInt i = 0; i < m_maxNoSets; i++) {
      m_gBndryCells[i].clear();
    }
  }

  m_localMarksteinLength = nullptr;
  if(m_useLocalMarksteinLength) {
    mAlloc(m_localMarksteinLength, m_maxNoCells, "m_localMarksteinLength", F0, AT_);
    mTerm(1, AT_, "code should be debugged for multiple sets");
  }

  // memorz for reinitialisation
  if(!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2) {
    mAlloc(m_cellList, m_maxNoCells, "m_cellList", 0, AT_);
    mAlloc(m_phiRatioCells, m_maxNoCells, 2 * nDim * m_maxNoSets, "m_phiRatioCells", 0, AT_);
    mAlloc(m_correction, m_maxNoCells * m_maxNoSets, "m_correction", F0, AT_);
    mAlloc(m_d, m_maxNoCells * m_maxNoSets, "m_d", F0, AT_);
    mAlloc(m_phiRatio, m_maxNoCells, 2 * nDim * m_maxNoSets, "m_phiRatio", F0, AT_);
    mAlloc(m_signG, m_maxNoCells * m_maxNoSets, "m_signG", F0, AT_);
  }

  // Parallelization
  // allocate space for buffers
  if(noNeighborDomains() > 0) {
    if(isActive() && grid().noNeighborDomains() > 0) {
      // allocate space for send and receive buffers
      if(!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2) {
        mAlloc(m_gSendBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gSendBuffers", F0, AT_);
        mAlloc(m_gReceiveBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gReceiveBuffers", F0, AT_);
      }

      if(m_combustion) {
        mAlloc(m_intSendBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intSendBuffers", 0, AT_);
        mAlloc(m_intReceiveBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intReceiveBuffers", 0, AT_);
      }

      if(m_combustion || (!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2)) {
        mAlloc(mpi_request, grid().noNeighborDomains(), "mpi_request", AT_);
        mAlloc(mpi_recive, grid().noNeighborDomains(), "mpi_recive", AT_);
      }
    }
  }

  // initialize the flame surface area
  m_arcLength = F0;
}

/** \brief Discretizes the level set equation
 *
 * \author: Thomas Hoesgen
 * \date Februar 2020
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::allocateRotatingLs() {
  TRACE();

  mAlloc(m_bodyRadius, m_noEmbeddedBodies, "m_bodyRadius", F0, AT_);
  mAlloc(m_omega, m_noEmbeddedBodies * nDim, "m_omega", F0, AT_);
  mAlloc(m_bodyAngularVelocity, m_noEmbeddedBodies * nDim, "m_bodyAngularVelocity", F0, AT_);
  mAlloc(m_bodyAngularAcceleration, m_noEmbeddedBodies * nDim, "m_bodyAngularAcceleration", F0, AT_);
  mAlloc(m_semiLagrange_xRot_ref, m_noEmbeddedBodies * nDim, "m_semiLagrange_xRot_ref", F0, AT_);
  mAlloc(m_semiLagrange_xRot_STL, m_noEmbeddedBodies * nDim, "m_semiLagrange_xRot_STL", F0, AT_);
  if(!m_reconstructOldG) {
    mAlloc(m_initialGCell, m_maxNoCells, "m_initialGCell", 0, AT_);
    mAlloc(m_cellDomIds, 2 * m_maxNoCells, "m_cellDomIds", -1, AT_);
  }

  MFloat maRot = F0;
  MFloat rad = F0;
  MInt ind = -1;
  for(MInt b = 0; b < m_noEmbeddedBodies; b++) {
    rad = Context::getSolverProperty<MFloat>("bodyRadius", solverId(), AT_, b);
    m_bodyRadius[b] = rad;
    for(MInt i = 0; i < nDim; i++) {
      ind = b * nDim + i;
      maRot = Context::getSolverProperty<MFloat>("MaRot", solverId(), AT_, ind);
      if(abs(maRot) > F0) {
        m_bodiesToCompute.push_back(b);
        break;
      }
    }
  }

  m_noBodiesToCompute = m_bodiesToCompute.size();

  // Communication buffers only needed when solver is active
  if(!m_reconstructOldG && isActive()) {
    // Global Communication
    mAlloc(m_globalSndOffsets, grid().noDomains() + 1, "m_globalSndOffsets", 0, AT_);
    mAlloc(m_globalRcvOffsets, grid().noDomains() + 1, "m_globalRcvOffsets", 0, AT_);
  }
}

// ----------------------------------------------------------------------------------------


/** \brief Discretizes the level set equation
 *
 * \author: Daniel Hartmann
 * \date June 2007
 *
 *  The fifth-order UC scheme is automatically reduced to third and first-order when necessary
 *
 * \cleanup: Sohel Herff
 *  All discretization schemes except for UC5 have been removed since they are worse and/or not used.
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeLevelSetRHS() {
  TRACE();

  MInt nghbrL;
  MInt nghbrL2;
  MInt nghbrL3;
  MInt nghbrR;
  MInt nghbrR2;
  MInt nghbrR3;
  MInt id;
  MInt cellId;
  //--- end of initialization
  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    // reset the G right hand side
    for(MInt bandId = 0; bandId < a_noBandCells(set); bandId++) {
      cellId = a_bandCellId(bandId, set);
      a_levelSetRHS(cellId, set) = F0;
    }
    // UC5
    for(MInt bandId = 0; bandId < a_noBandCells(set); bandId++) {
      cellId = a_bandCellId(bandId, set);
      if(a_isHalo(cellId)) continue;
      for(MInt i = 0; i < nDim; i++) {
        if(a_extensionVelocityG(cellId, i, set) > F0) {
          nghbrL = c_neighborId(cellId, 2 * i);
          nghbrL2 = c_neighborId(nghbrL, 2 * i);
          nghbrL3 = c_neighborId(nghbrL2, 2 * i);
          nghbrR = c_neighborId(cellId, 2 * i + 1);
          nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);
          id = cellId;
          if((!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR, set)) /*||
                                                                           (a_isGBoundaryCellG( nghbrL2 ,  set ) ||
                                                                           a_isGBoundaryCellG( nghbrR ,  set ) )*/
          ) {
            // reduce to first order
            nghbrR = cellId;
            nghbrR2 = cellId;
            nghbrL2 = cellId;
            nghbrL3 = cellId;
          } else {
            if((!a_inBandG(nghbrL3, set) || !a_inBandG(nghbrR2, set)) /*||
                                                                              (a_isGBoundaryCellG( nghbrL3 ,  set ) ||
                                                                              a_isGBoundaryCellG( nghbrR2 ,  set ) )*/
            ) {
              // reduce to third order
              id = nghbrR;
              nghbrR = cellId;
              nghbrR2 = nghbrL2;
              nghbrL3 = nghbrL2;
            }
          }
          a_levelSetRHS(cellId, set) +=
              a_extensionVelocityG(cellId, i, set)
              * (-F1B30 * a_levelSetFunctionG(nghbrL3, set) + F1B4 * a_levelSetFunctionG(nghbrL2, set)
                 - a_levelSetFunctionG(nghbrL, set) + F1B3 * a_levelSetFunctionG(id, set)
                 + F1B2 * a_levelSetFunctionG(nghbrR, set) - F1B20 * a_levelSetFunctionG(nghbrR2, set));
        } else {
          nghbrL = c_neighborId(cellId, 2 * i);
          nghbrL2 = c_neighborId(nghbrL, 2 * i);
          nghbrR = c_neighborId(cellId, 2 * i + 1);
          nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);
          nghbrR3 = c_neighborId(nghbrR2, 2 * i + 1);
          id = cellId;
          if((!a_inBandG(nghbrR2, set) || !a_inBandG(nghbrL, set)) /*||
                                                                           (a_isGBoundaryCellG( nghbrR2 ,  set ) ||
                                                                           a_isGBoundaryCellG( nghbrL ,  set ) )*/
          ) {
            // reduce to first order
            nghbrL = cellId;
            nghbrL2 = cellId;
            nghbrR2 = cellId;
            nghbrR3 = cellId;
          } else {
            if((!a_inBandG(nghbrR3, set) || !a_inBandG(nghbrL2, set)) /*||
                                                                              (a_isGBoundaryCellG( nghbrR3 ,  set ) ||
                                                                              a_isGBoundaryCellG( nghbrL2 ,  set ) )*/
            ) {
              // reduce to third order
              id = nghbrL;
              nghbrL = cellId;
              nghbrL2 = nghbrR2;
              nghbrR3 = nghbrR2;
            }
          }
          a_levelSetRHS(cellId, set) +=
              a_extensionVelocityG(cellId, i, set)
              * (F1B30 * a_levelSetFunctionG(nghbrR3, set) - F1B4 * a_levelSetFunctionG(nghbrR2, set)
                 + a_levelSetFunctionG(nghbrR, set) - F1B3 * a_levelSetFunctionG(id, set)
                 - F1B2 * a_levelSetFunctionG(nghbrL, set) + F1B20 * a_levelSetFunctionG(nghbrL2, set));
        }
      }
      a_levelSetRHS(cellId, set) *= m_FgCellDistance;
    }
  }
}


// ----------------------------------------------------------------------------------------

/** \brief computes Level Set
 *
 * computes the level set RHS,
 * advances the G field via the gRungeKutta Method
 * and averages G from the fine grid to the coarser grid
 *
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::levelSetSolver() {
  TRACE();

  MBool timeStepCompleted = false;

  // exchange level set function on halo cells, cause they are needed for the Markstein length computation!

  computeLevelSetRHS();

  // advance the G field
  timeStepCompleted = gRungeKutta();

  // average G from the fine grid to the coarse grid
  if(timeStepCompleted) levelSetRestriction();

  return timeStepCompleted;
}

//-----------------------------------------------------------------------------


template <MInt nDim>
MBool LsCartesianSolver<nDim>::gRungeKutta() {
  TRACE();
  MFloat c;
  MFloat factor;
  //---end of initialization
  if(m_gRKMethod == 5) return true;

  // set old variables - changes due to multilevel
  for(MInt set = 0; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    if(m_gRKStep == 0)
      for(MInt id = 0; id < a_noBandCells(set); id++)
        a_oldLevelSetFunctionG(a_bandCellId(id, set), set) = a_levelSetFunctionG(a_bandCellId(id, set), set);
  }

  switch(m_gRKMethod) {
    case 0: {
      factor = m_gRKalpha[m_gRKStep] * timeStep();
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          if(a_isGBoundaryCellG(a_bandCellId(id, set), set)) continue;
          a_levelSetFunctionG(a_bandCellId(id, set), set) =
              a_oldLevelSetFunctionG(a_bandCellId(id, set), set) - factor * a_levelSetRHS(a_bandCellId(id, set), set);
        }
      }
      break;
    }
    case 1: {
      factor = m_gRKalpha[m_gRKStep] * timeStep();
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          if(a_isGBoundaryCellG(a_bandCellId(id, set), set)) continue;
          a_levelSetFunctionG(a_bandCellId(id, set), set) =
              a_levelSetFunctionG(a_bandCellId(id, set), set) * m_gRKalpha[m_gRKStep]
              + a_oldLevelSetFunctionG(a_bandCellId(id, set), set) * (F1 - m_gRKalpha[m_gRKStep])
              - a_levelSetRHS(a_bandCellId(id, set), set) * factor;
        }
      }
      break;
    }
    case 2: {
      MInt noIntegrationLayers = m_gBandWidth - 5;
      factor = m_gRKalpha[m_gRKStep] * timeStep();
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          if(a_isGBoundaryCellG(a_bandCellId(id, set), set)) continue;
          c = fabs(a_levelSetFunctionG(a_bandCellId(id, set), set)) * m_FgCellDistance - noIntegrationLayers;
          if(c > F0)
            a_levelSetRHS(a_bandCellId(id, set), set) = F0;
          else if(c > -3.0)
            a_levelSetRHS(a_bandCellId(id, set), set) *= (F2B27 * POW3(c) + F1B3 * POW2(c));
          a_levelSetFunctionG(a_bandCellId(id, set), set) =
              a_levelSetFunctionG(a_bandCellId(id, set), set) * m_gRKalpha[m_gRKStep]
              + a_oldLevelSetFunctionG(a_bandCellId(id, set), set) * (F1 - m_gRKalpha[m_gRKStep])
              - a_levelSetRHS(a_bandCellId(id, set), set) * factor;
        }
      }
      break;
    }
    case 3: {
      factor = m_gRKalpha[m_gRKStep] * timeStep();
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_bandLayer(0, set); id++) {
          if(a_isGBoundaryCellG(a_bandCellId(id, set), set)) continue;
          a_levelSetFunctionG(a_bandCellId(id, set), set) =
              a_levelSetFunctionG(a_bandCellId(id, set), set) * m_gRKalpha[m_gRKStep]
              + a_oldLevelSetFunctionG(a_bandCellId(id, set), set) * (F1 - m_gRKalpha[m_gRKStep])
              - a_levelSetRHS(a_bandCellId(id, set), set) * factor;
        }
      }
      reinitBand(m_startSet, m_noSets);
      break;
    }
      // the level-set equation is only solved at the front
    case 4: {
      factor = m_gRKalpha[m_gRKStep] * timeStep();
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_bandLayer(0, set); id++) {
          if(a_isGBoundaryCellG(a_bandCellId(id, set), set)) continue;
          a_levelSetFunctionG(a_bandCellId(id, set), set) =
              a_levelSetFunctionG(a_bandCellId(id, set), set) * m_gRKalpha[m_gRKStep]
              + a_oldLevelSetFunctionG(a_bandCellId(id, set), set) * (F1 - m_gRKalpha[m_gRKStep])
              - a_levelSetRHS(a_bandCellId(id, set), set) * factor;
        }
      }
      break;
    }
    case 5: {
      // property is initialized to this value, which does nothing
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown gRKMethod");
    }
  }
  // cerr << "m_gRKStep: " << m_gRKStep << endl;
  m_gRKStep++;

  if(m_gRKStep == m_nogRKSteps) {
    m_gRKStep = 0;
    if(m_LSSolver) {
      if(globalTimeStep > 0) {
        m_time += timeStep();
      }
    }
    return true;
  } else
    return false;
}


//-----------------------------------------------------------------------------


/** \brief Discretizes the level set equation - semi Lagrangian solution scheme
 *
 * \author: Claudia Guenther
 * \date 04/2012
 *
 *  implemented schemes: BACKWARDS_PAR, ROTATING_LS
 *
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::semiLagrangeTimeStep() {
  TRACE();

  switch(string2enum(m_levelSetDiscretizationScheme)) {
    case BACKWARDS_PAR: {
      MFloat xOld[3] = {F0, F0, F0};
      MInt containingCell = 0;
      MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};

      MInt position = 0;
      MFloat xCurrent[3] = {F0, F0, F0};
      MFloat xShift[3] = {F0, F0, F0};
      MFloat xShift_cur[3] = {F0, F0, F0};
#if defined LS_DEBUG || !defined NDEBUG
      MFloat shiftCheck[3] = {F0, F0, F0};
#endif

      //-------------------------

      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;

        const MFloat cellLength = c_cellLengthAtLevel(a_maxGCellLevel(set));
        const MFloat cellHalfLength = c_cellLengthAtLevel(a_maxGCellLevel(set) + 1);

        for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
          const MInt body = m_setToBodiesTable[set][b];
          computeBodyPropertiesForced(1, xCurrent, body, time() + timeStep(), true);

#if defined LS_DEBUG || !defined NDEBUG
          computeBodyPropertiesForced(1, shiftCheck, body, time());
          for(MInt d = 0; d < nDim; d++) {
            MFloat shift = xCurrent[d] - shiftCheck[d];
            if(shift > cellLength && !m_periodicMovement) {
              mTerm(1, AT_, "LevelSet-Shift is exceeding maximum shift-limit of a cellLength! Reduce timestep!");
            }
          }
#endif
          computeBodyPropertiesForced(1, xOld, body, 0.0);
          for(MInt d = 0; d < nDim; d++) {
            xShift[d] = xCurrent[d] - xOld[d];
            xShift_cur[d] = xShift[d] - m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body];
          }

          // 0. if required, shift reference field
          for(MInt d = 0; d < nDim; d++) {
            MFloat tmp = xShift_cur[d];
            if(abs(tmp) > 2 * cellLength && m_periodicMovement) {
              m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body] -= m_periodicDistance - cellLength;
              shiftOldLevelSetField(2 * d + 1, set, body);
              exchangeLeafDataLS<false>();
            } else {
              if(tmp + cellHalfLength > cellLength) {
                m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body] += cellLength;
                shiftOldLevelSetField(2 * d + 1, set, body);
                exchangeLeafDataLS<false>();
              } else if(tmp - cellHalfLength < -cellLength) {
                m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body] -= cellLength;
                shiftOldLevelSetField(2 * d, set, body);
                exchangeLeafDataLS<false>();
              }
            }
          }

          // interpolate current field from reference field
          for(MInt id = 0; id < a_noBandCells(set); id++) {
            const MInt cellId = a_bandCellId(id, set);
            if(!(a_bodyIdG(cellId, set) == body)) continue;
            // if cell is a boundary cell, process simle (0th order)
            MBool boundaryCell = false;
            for(MInt d = 0; d < m_noDirs; d++) {
              if(!a_hasNeighbor(cellId, d)) {
                a_levelSetFunctionG(cellId, set) = a_oldLevelSetFunctionG(cellId, set);
                boundaryCell = true;
                break;
              }
            }
            if(boundaryCell) continue;
            if(a_isHalo(cellId)) continue;
            // else, 1st order interpolation:
            // 1. find x_old
            for(MInt i = 0; i < nDim; i++) {
              xOld[i] =
                  c_coordinate(cellId, i) - (xShift[i] - m_semiLagrange_xShift_ref[i * m_noEmbeddedBodies + body]);
            }
            // 2. get containing cell
            containingCell = getContainingCell(cellId, xOld, set);
            // 3. set up interpolation stencil
            position = setUpLevelSetInterpolationStencil(containingCell, interpolationCells, xOld);
            // 4. interpolate level set
            MFloat phiNew = -99;
            if(position > -1) {
              phiNew = interpolateOldLevelSet(interpolationCells, xOld, set);
            } else {
              phiNew = a_oldLevelSetFunctionG(containingCell, set);
            }
            a_levelSetFunctionG(cellId, set) = phiNew;
          }
        }
      }
      break;
    }
    case ROTATING_LS: {
      MFloat xInitial[3];
      MFloat xCoord[3];
      MFloat xOld[3];
      MFloat xCurrent[3];
      MInt body;
      MInt set;
      MFloat phiNew = F0;
      MInt containingCell;
      MInt searchCell;
      MInt searchDomain;
      MInt cellId;
      MInt ind;

      MBool firstRun = false;
      std::vector<MInt> remCells;

      // Update the rotated angle
      for(MInt i = 0; i < m_noBodiesToCompute; i++) {
        body = m_bodiesToCompute[i];
        for(MInt d = 0; d < nDim; d++) {
          m_semiLagrange_xRot_ref[body * nDim + d] += m_omega[body * nDim + d] * timeStep();
          m_semiLagrange_xRot_STL[body * nDim + d] += m_omega[body * nDim + d] * timeStep();
        }

#ifndef NDEBUG
        if(domainId() == 0) {
          cerr << "TS: " << globalTimeStep << " B: " << body
               << " rot: " << m_semiLagrange_xRot_STL[body * nDim + 0] * 180.0 / PI << ","
               << m_semiLagrange_xRot_STL[body * nDim + 1] * 180.0 / PI << ","
               << m_semiLagrange_xRot_STL[body * nDim + 2] * 180.0 / PI << " / "
               << m_semiLagrange_xRot_ref[body * nDim + 0] * 180.0 / PI << ","
               << m_semiLagrange_xRot_ref[body * nDim + 1] * 180.0 / PI << ","
               << m_semiLagrange_xRot_ref[body * nDim + 2] * 180.0 / PI << " deg;"
               << " Omega: " << m_omega[body * nDim + 0] << " " << m_omega[body * nDim + 1] << " "
               << m_omega[body * nDim + 2] << endl;
        }
#endif
      }

      // Compute levelset
      if(m_reconstructOldG) {
        // Update halo information of containingCells list
        MInt noData = m_noBodiesToCompute;
        MLongScratchSpace tmp_data(a_noCells(), noData, AT_, "tmp_data");
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noWindowCells(i); j++) {
              cellId = windowCellId(i, j);
              if(a_containingCell(cellId, b) > -1) {
                tmp_data(cellId, b) = c_globalId(a_containingCell(cellId, b));
              } else {
                tmp_data(cellId, b) = -1;
              }
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
              cellId = grid().azimuthalWindowCell(i, j);
              if(a_containingCell(cellId, b) > -1) {
                tmp_data(cellId, b) = c_globalId(a_containingCell(cellId, b));
              } else {
                tmp_data(cellId, b) = -1;
              }
            }
          }
        }
        exchangeDataLS(&tmp_data(0, 0), noData);
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noHaloCells(i); j++) {
              cellId = haloCellId(i, j);
              if(tmp_data(cellId, b) > -1) {
                a_containingCell(cellId, b) = a_localId(tmp_data(cellId, b));
              } else {
                a_containingCell(cellId, b) = -1;
              }
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
              cellId = grid().azimuthalHaloCell(i, j);
              if(tmp_data(cellId, b) > -1) {
                a_containingCell(cellId, b) = a_localId(tmp_data(cellId, b));
              } else {
                a_containingCell(cellId, b) = -1;
              }
            }
          }
        }

        // Loop over all sets and compute levelset value and update containingCells list
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          body = m_bodiesToCompute[b];
          set = m_bodyToSetTable[body];

          // Compute position of bodycenter
          computeBodyPropertiesForced(1, xCurrent, body, time() + timeStep());
          computeBodyPropertiesForced(1, xOld, body, 0.0);

          for(MInt id = 0; id < a_noBandCells(set); id++) {
            cellId = a_bandCellId(id, set);
            if(!(a_bodyIdG(cellId, set) == body)) continue;
            if(a_isHalo(cellId)) continue;

            searchCell = a_containingCell(cellId, b);
            if(searchCell < 0) {
              for(MInt i = 0; i < nDim; i++)
                xCoord[i] = c_coordinate(cellId, i);

              getContainingCellFromNeighbor(b, cellId, xCoord, xOld);
              m_newCells.push_back(cellId);
            }
          }
          m_newCells.clear();

          // Loop over all bandCells in set
          for(MInt id = 0; id < a_noBandCells(set); id++) {
            cellId = a_bandCellId(id, set);
            if(!(a_bodyIdG(cellId, set) == body)) continue;
            if(a_isHalo(cellId)) continue;

            // Access information about which cell lies at xInitial. If no information is available take it from
            // neighborCell
            searchCell = a_containingCell(cellId, b);

            // Compute the reference position of the bandCell at t=0
            for(MInt i = 0; i < nDim; i++) {
              xCoord[i] = c_coordinate(cellId, i);
            }
            rotateLevelSet(1, xInitial, body, xCoord, xOld, &m_semiLagrange_xRot_ref[body * nDim]);
            for(MInt i = 0; i < nDim; i++) {
              xInitial[i] += xCurrent[i] - xOld[i];
            }

            if(grid().azimuthalPeriodicity()) {
              MBool shift = false;
              for(MInt d = 0; d < nDim; d++) {
                shift = shift || !approx(c_coordinate(searchCell, d), xInitial[d], F2 * c_cellLengthAtCell(searchCell));
              }
              if(shift) {
                MFloat coords[nDim];
                MFloat coordsCylSearch[nDim];
                MFloat coordsCyl[nDim];
                for(MInt d = 0; d < nDim; d++) {
                  coords[d] = c_coordinate(searchCell, d);
                }
                grid().raw().cartesianToCylindric(coords, coordsCylSearch);
                grid().raw().cartesianToCylindric(xInitial, coordsCyl);
                MInt side = grid().determineAzimuthalBoundarySide(xInitial);
                MInt fac = 0;
                if(side == -1) {
                  fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() - F1B2);
                } else if(side == 1) {
                  fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() + F1B2);
                } else {
                  mTerm(1, AT_, "Invalid side!");
                }
                grid().raw().rotateCartesianCoordinates(xInitial, fac * grid().azimuthalAngle());
              }
            }

            // Find containingCell
            containingCell = getContainingCell(searchCell, xInitial);
            if(containingCell != cellId) {
              m_rotatingReinitTrigger = 1;
            }

            // Compute new levelSet value and determine containingCell. Then update list.
            MInt dummyDomain = -1;
            processRotatingLevelSet(phiNew, containingCell, dummyDomain, xInitial, set);

            a_levelSetFunctionG(cellId, set) = phiNew;
            a_containingCell(cellId, b) = containingCell;
          }

          // Reset cells in shadow layer
          for(MInt i = 0; i < m_maxNoCells; i++) {
            if(!a_inBandG(i, set)) {
              a_containingCell(i, b) = -1;
            }
          }
        }

      } else {
        // Update halo information of containingCells list
        MInt noData = 2 * m_noBodiesToCompute;
        MIntScratchSpace tmp_data(a_noCells(), noData, AT_, "tmp_data");
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noWindowCells(i); j++) {
              cellId = windowCellId(i, j);
              if(a_containingCell(cellId, b) > -1) {
                tmp_data(cellId, b) = a_containingCell(cellId, b);
                tmp_data(cellId, m_noBodiesToCompute + b) = a_containingDomain(cellId, b);
              } else {
                tmp_data(cellId, b) = -1;
                tmp_data(cellId, m_noBodiesToCompute + b) = -1;
              }
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
              cellId = grid().azimuthalWindowCell(i, j);
              if(a_containingCell(cellId, b) > -1) {
                tmp_data(cellId, b) = a_containingCell(cellId, b);
                tmp_data(cellId, m_noBodiesToCompute + b) = a_containingDomain(cellId, b);
              } else {
                tmp_data(cellId, b) = -1;
                tmp_data(cellId, m_noBodiesToCompute + b) = -1;
              }
            }
          }
        }
        exchangeDataLS(&tmp_data(0, 0), noData);
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noHaloCells(i); j++) {
              cellId = haloCellId(i, j);
              if(tmp_data(cellId, b) > -1) {
                a_containingCell(cellId, b) = tmp_data(cellId, b);
                a_containingDomain(cellId, b) = tmp_data(cellId, m_noBodiesToCompute + b);
              } else {
                a_containingCell(cellId, b) = -1;
                a_containingDomain(cellId, b) = -1;
              }
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
              cellId = grid().azimuthalHaloCell(i, j);
              if(tmp_data(cellId, b) > -1) {
                a_containingCell(cellId, b) = tmp_data(cellId, b);
                a_containingDomain(cellId, b) = tmp_data(cellId, m_noBodiesToCompute + b);
              } else {
                a_containingCell(cellId, b) = -1;
                a_containingDomain(cellId, b) = -1;
              }
            }
          }
        }

        // Loop over all sets and compute levelset value and update containingCells list
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          body = m_bodiesToCompute[b];
          set = m_bodyToSetTable[body];

          // Compute position of bodycenter
          computeBodyPropertiesForced(1, xCurrent, body, time() + timeStep());
          computeBodyPropertiesForced(1, xOld, body, 0.0);

          // Allocate Buffers for global communication
          MIntScratchSpace noCellsToDom(grid().noDomains(), AT_, "noCellsToDom");
          noCellsToDom.fill(0);

          for(MInt id = 0; id < a_noBandCells(set); id++) {
            cellId = a_bandCellId(id, set);
            if(!(a_bodyIdG(cellId, set) == body)) continue;
            if(a_isHalo(cellId)) continue;

            searchCell = a_containingCell(cellId, b);
            searchDomain = a_containingDomain(cellId, b);
            if(searchCell < 0 || searchDomain < 0) {
              for(MInt i = 0; i < nDim; i++)
                xCoord[i] = c_coordinate(cellId, i);
              getContainingCellFromNeighbor(b, cellId, xCoord, xOld);
              searchDomain = a_containingDomain(cellId, b);
              m_newCells.push_back(cellId);
            }
            if(searchDomain > -1) {
              noCellsToDom[searchDomain]++;
            }
          }
          m_newCells.clear();

          prepareGlobalComm(&noCellsToDom[0]);
          MInt noCellsComm = mMax(m_globalSndOffsets[grid().noDomains()], m_globalRcvOffsets[grid().noDomains()]);

          MIntScratchSpace intSndBufSizeGlob(grid().noDomains(), AT_, "intSndBufSizeGlob");
          MIntScratchSpace intRcvBufSizeGlob(grid().noDomains(), AT_, "intRcvBufSizeGlob");
          MIntScratchSpace floatSndBufSizeGlob(grid().noDomains(), AT_, "floatSndBufSizeGlob");
          MIntScratchSpace floatRcvBufSizeGlob(grid().noDomains(), AT_, "floatRcvBufSizeGlob");
          MInt floatOffset = 3;
          MInt intOffset = 2;
          MIntScratchSpace intSndBufGlob(intOffset * noCellsComm, AT_, "intSndBufGlob");
          MIntScratchSpace intRcvBufGlob(intOffset * noCellsComm, AT_, "intRcvBufGlob");
          MFloatScratchSpace floatSndBufGlob(floatOffset * noCellsComm, AT_, "floatSndBufGlob");
          MFloatScratchSpace floatRcvBufGlob(floatOffset * noCellsComm, AT_, "floatRcvBufGlob");

          std::vector<std::vector<MInt>> cellIdsLoc;
          cellIdsLoc.resize(grid().noDomains());

          intSndBufGlob.fill(-1);
          intRcvBufGlob.fill(-1);
          floatSndBufGlob.fill(-F1);
          floatRcvBufGlob.fill(-F1);
          intSndBufSizeGlob.fill(0);
          intRcvBufSizeGlob.fill(0);
          floatSndBufSizeGlob.fill(0);
          floatRcvBufSizeGlob.fill(0);

          // Loop over all bandCells in set
          for(MInt id = 0; id < a_noBandCells(set); id++) {
            cellId = a_bandCellId(id, set);
            if(!(a_bodyIdG(cellId, set) == body)) continue;
            if(a_isHalo(cellId)) continue;

            // Access information about which cell lies at xInitial. If no information is available take it from
            // neighborCell
            searchCell = a_containingCell(cellId, b);
            searchDomain = a_containingDomain(cellId, b);

            // Skip outer layer of bandCells.
            /*
              if ( a_isGBoundaryCellG(cellId,set) ) {
              if ( a_levelSetFunctionG(cellId,set) > F0 ) {
              a_levelSetFunctionG(cellId,set) = m_outsideGValue;
              } else {
              a_levelSetFunctionG(cellId,set) = -m_outsideGValue;
              }
              m_containingCell[ b*m_maxNoCells + cellId ] = -1;
              m_containingDomain[ b*m_maxNoCells + cellId ] = -1;
              continue;
              }
            */
            // Compute the reference position of the bandCell at t=0
            for(MInt i = 0; i < nDim; i++) {
              xCoord[i] = c_coordinate(cellId, i);
            }
            rotateLevelSet(1, xInitial, body, xCoord, xOld, &m_semiLagrange_xRot_ref[body * nDim]);
            for(MInt i = 0; i < nDim; i++) {
              xInitial[i] += xCurrent[i] - xOld[i];
            }

            // Now, handle all cases where cell and containingCell are on the same domain (no communication neccessary).
            if(searchDomain == domainId()) {
              if(grid().azimuthalPeriodicity()) {
                MBool shift = false;
                for(MInt d = 0; d < nDim; d++) {
                  shift =
                      shift || !approx(c_coordinate(searchCell, d), xInitial[d], F2 * c_cellLengthAtCell(searchCell));
                }
                if(shift) {
                  MFloat coords[nDim];
                  MFloat coordsCylSearch[nDim];
                  MFloat coordsCyl[nDim];
                  for(MInt d = 0; d < nDim; d++) {
                    coords[d] = c_coordinate(searchCell, d);
                  }
                  grid().raw().cartesianToCylindric(coords, coordsCylSearch);
                  grid().raw().cartesianToCylindric(xInitial, coordsCyl);
                  MInt side = grid().determineAzimuthalBoundarySide(xInitial);
                  MInt fac = 0;
                  if(side == -1) {
                    fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() - F1B2);
                  } else if(side == 1) {
                    fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() + F1B2);
                  } else {
                    mTerm(1, AT_, "Invalid side!");
                  }
                  grid().raw().rotateCartesianCoordinates(xInitial, fac * grid().azimuthalAngle());
                }
              }

              // Find containingCell
              containingCell = getContainingCell(searchCell, xInitial);

              // This is a backup, if the list is incorrect after restarts
              if(firstRun && containingCell < 0) {
                remCells.push_back(cellId);
                continue;
              }
              // Compute new levelSet value and determine containingCell and containingDomain. Then update list.
              processRotatingLevelSet(phiNew, containingCell, searchDomain, xInitial, set);

              a_levelSetFunctionG(cellId, set) = phiNew;
              a_containingCell(cellId, b) = containingCell;
              a_containingDomain(cellId, b) = searchDomain;
            } else if(searchDomain < 0) {
              // If the list is broken after restart, use backup method. Otherwise, exit...
              if(firstRun) {
                remCells.push_back(cellId);
              } else {
                mTerm(1, AT_, "Containing Domain unkown in ROTATING_LS!");
              }
            } else {
              // If containingCell is on different domain. Put all neccessary information in send buffer to exchange
              // latter.
              cellIdsLoc[searchDomain].push_back(cellId);
              intSndBufGlob[m_globalSndOffsets[searchDomain] + intSndBufSizeGlob.p[searchDomain]] = searchCell;
              intSndBufSizeGlob.p[searchDomain]++;

              for(MInt dim = 0; dim < nDim; dim++) {
                floatSndBufGlob[floatOffset * m_globalSndOffsets[searchDomain] + floatSndBufSizeGlob.p[searchDomain]] =
                    xInitial[dim];
                floatSndBufSizeGlob.p[searchDomain]++;
              }
            }
          }


          // exchange
          exchangeBuffersGlobal(intSndBufGlob.getPointer(), intRcvBufGlob.getPointer(), intSndBufSizeGlob.getPointer(),
                                intRcvBufSizeGlob.getPointer(), m_globalSndOffsets, m_globalRcvOffsets, 3);
          exchangeBuffersGlobal(floatSndBufGlob.getPointer(), floatRcvBufGlob.getPointer(),
                                floatSndBufSizeGlob.getPointer(), floatRcvBufSizeGlob.getPointer(), m_globalSndOffsets,
                                m_globalRcvOffsets, 5, floatOffset);

          intSndBufGlob.fill(-1);
          floatSndBufGlob.fill(-F1);
          intSndBufSizeGlob.fill(0);
          floatSndBufSizeGlob.fill(0);

          // Now handle all cells for which the containingCell is on different domain
          for(MInt i = 0; i < grid().noDomains(); i++) {
            ind = m_globalRcvOffsets[i];
            for(MInt j = 0; j < intRcvBufSizeGlob(i); j++) {
              searchCell = intRcvBufGlob(ind + j);
              searchDomain = domainId();

              xInitial[0] = floatRcvBufGlob(ind * floatOffset + j * floatOffset + 0);
              xInitial[1] = floatRcvBufGlob(ind * floatOffset + j * floatOffset + 1);
              xInitial[2] = floatRcvBufGlob(ind * floatOffset + j * floatOffset + 2);

              if(grid().azimuthalPeriodicity()) {
                MBool shift = false;
                for(MInt d = 0; d < nDim; d++) {
                  shift =
                      shift || !approx(c_coordinate(searchCell, d), xInitial[d], F2 * c_cellLengthAtCell(searchCell));
                }
                if(shift) {
                  MFloat coords[nDim];
                  MFloat coordsCylSearch[nDim];
                  MFloat coordsCyl[nDim];
                  for(MInt d = 0; d < nDim; d++) {
                    coords[d] = c_coordinate(searchCell, d);
                  }
                  grid().raw().cartesianToCylindric(coords, coordsCylSearch);
                  grid().raw().cartesianToCylindric(xInitial, coordsCyl);
                  MInt side = grid().determineAzimuthalBoundarySide(xInitial);
                  MInt fac = 0;
                  if(side == -1) {
                    fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() - F1B2);
                  } else if(side == 1) {
                    fac = (MInt)((coordsCyl[1] - coordsCylSearch[1]) / grid().azimuthalAngle() + F1B2);
                  } else {
                    mTerm(1, AT_, "Invalid side!");
                  }
                  grid().raw().rotateCartesianCoordinates(xInitial, fac * grid().azimuthalAngle());
                }
              }

              // Find containingCell
              containingCell = getContainingCell(searchCell, xInitial);

              // Backup for restart
              if(firstRun && containingCell < 0) {
                intSndBufGlob(ind * intOffset + intSndBufSizeGlob.p[i]) = -2;
                intSndBufSizeGlob.p[i]++;
                intSndBufGlob(ind * intOffset + intSndBufSizeGlob.p[i]) = -2;
                intSndBufSizeGlob.p[i]++;
                floatSndBufGlob(ind * floatOffset + floatSndBufSizeGlob.p[i]) = F0;
                floatSndBufSizeGlob.p[i]++;
                continue;
              }

              // Compute levelSet value and new containingCell and containingDomain and put it Buffer for later exchange
              processRotatingLevelSet(phiNew, containingCell, searchDomain, xInitial, set);

              floatSndBufGlob(ind + floatSndBufSizeGlob.p[i]) = phiNew;
              floatSndBufSizeGlob.p[i]++;
              intSndBufGlob(ind * intOffset + intSndBufSizeGlob.p[i]) = containingCell;
              intSndBufSizeGlob.p[i]++;
              intSndBufGlob(ind * intOffset + intSndBufSizeGlob.p[i]) = searchDomain;
              intSndBufSizeGlob.p[i]++;
            }
          }

          intRcvBufGlob.fill(-1);
          floatRcvBufGlob.fill(-F1);
          intRcvBufSizeGlob.fill(0);
          floatRcvBufSizeGlob.fill(0);

          // exchange
          exchangeBuffersGlobal(intSndBufGlob.getPointer(), intRcvBufGlob.getPointer(), intSndBufSizeGlob.getPointer(),
                                intRcvBufSizeGlob.getPointer(), &m_globalRcvOffsets[0], &m_globalSndOffsets[0], 3,
                                intOffset);
          exchangeBuffersGlobal(floatSndBufGlob.getPointer(), floatRcvBufGlob.getPointer(),
                                floatSndBufSizeGlob.getPointer(), floatRcvBufSizeGlob.getPointer(),
                                &m_globalRcvOffsets[0], &m_globalSndOffsets[0], 5);

          // Retrive information from other domains and update levelSet and lists
          for(MInt i = 0; i < grid().noDomains(); i++) {
            ind = m_globalSndOffsets[i];
            for(MInt j = 0; j < floatRcvBufSizeGlob(i); j++) {
              cellId = cellIdsLoc[i][j];

              a_containingCell(cellId, b) = intRcvBufGlob(ind * intOffset + j * intOffset + 0);
              a_containingDomain(cellId, b) = intRcvBufGlob(ind * intOffset + j * intOffset + 1);

              // Backup solution at restart
              if(firstRun && a_containingCell(cellId, b) < -1) {
                remCells.push_back(cellId);
                continue;
              }

              a_levelSetFunctionG(cellId, set) = floatRcvBufGlob(ind + j);
            }
          }

          // Backup method for restart. Search each cell in remCells on each domain. Very slow!!!!!
          if(firstRun) {
            for(MInt dom = 0; dom < grid().noDomains(); dom++) {
              MInt cnt = 0;
              if(domainId() == dom) {
                cnt = remCells.size();
                if(cnt > 0) cerr << "D:" << domainId() << " Restart Backup LS!" << endl;
              }

              MPI_Bcast(&cnt, 1, MPI_INT, dom, globalMaiaCommWorld(), AT_, "cnt");

              for(MInt c = 0; c < cnt; c++) {
                containingCell = -1;
                searchDomain = -2;
                if(domainId() == dom) {
                  for(MInt i = 0; i < nDim; i++) {
                    xCoord[i] = c_coordinate(remCells[c], i);
                  }
                  rotateLevelSet(1, xInitial, body, xCoord, xOld, &m_semiLagrange_xRot_ref[body * nDim]);
                  for(MInt i = 0; i < nDim; i++) {
                    xInitial[i] += xCurrent[i] - xOld[i];
                  }
                }

                MPI_Bcast(&xInitial[0], 3, maia::type_traits<MFloat>::mpiType(), dom, globalMaiaCommWorld(), AT_,
                          "xInitial[0]");

                containingCell = getContainingCell(xInitial);

                if(containingCell > -1) {
                  // This should be save, since only internalCells are investigated in getContainingCell!
                  searchDomain = domainId();

                  processRotatingLevelSet(phiNew, containingCell, searchDomain, xInitial, set);
                }
                MPI_Allreduce(MPI_IN_PLACE, &searchDomain, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                              "searchDomain");

                if(searchDomain != dom) {
                  if(domainId() == searchDomain) {
                    MPI_Send(&containingCell, 1, MPI_INT, dom, 5, mpiComm(), AT_, "containingCell");
                    MPI_Send(&phiNew, 1, maia::type_traits<MFloat>::mpiType(), dom, 6, mpiComm(), AT_, "phiNew");
                  }
                  if(domainId() == dom) {
                    MPI_Recv(&containingCell, 1, MPI_INT, searchDomain, 5, mpiComm(), MPI_STATUS_IGNORE, AT_,
                             "containingCell");
                    MPI_Recv(&phiNew, 1, maia::type_traits<MFloat>::mpiType(), searchDomain, 6, mpiComm(),
                             MPI_STATUS_IGNORE, AT_, "phiNew");
                  }
                }

                if(domainId() == dom) {
                  cellId = remCells[c];
                  a_levelSetFunctionG(cellId, set) = phiNew;
                  a_containingCell(cellId, b) = containingCell;
                  a_containingDomain(cellId, b) = searchDomain;
                }
              }
            }
          }

          // Set firstRun to false
          firstRun = false;

          // Reset cells in shadow layer
          for(MInt i = 0; i < m_maxNoCells; i++) {
            if(!a_inBandG(i, set)) {
              a_containingCell(i, b) = -1;
              a_containingDomain(i, b) = -1;
            }
          }
        }
      }

      // Finaly, update the bodies angular velocities. Neccessary for boundry condition
      for(MInt i = 0; i < m_noBodiesToCompute; i++) {
        body = m_bodiesToCompute[i];
        rotateLevelSet(5, &m_bodyAngularVelocity[body * nDim], body, nullptr, nullptr,
                       &m_semiLagrange_xRot_STL[body * nDim]);
      }

      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "LsCartesianSolver::semiLagrangeTimeStep(): switch variable 'm_levelSetDiscretizationScheme' "
                      "with value "
                   << m_levelSetDiscretizationScheme << " not matching any case." << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  return true;
}
//-----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetConstrainedReinitialization(MInt methodId, MInt startSet, MInt endSet,
                                                                  MInt gapMode) {
  TRACE();

  NEW_TIMER_GROUP_STATIC(reInit, "Reinitialisation");
  NEW_TIMER_STATIC(t_reInit, "Total time - levelset Reinitialisation", reInit);
  NEW_SUB_TIMER_STATIC(t_c1, "preparation", t_reInit);
  NEW_SUB_TIMER_STATIC(t_c2, "solver", t_reInit);

  RECORD_TIMER_START(t_reInit);
  RECORD_TIMER_START(t_c1);

  // statistics only for start set! does not compute statistics for multiple level set functions!
  // #define REINITIALIZATION_STATISTICS

  MBool upwind;
  MInt cellId;
  MInt cellListSize;
  MIntScratchSpace nghbr(m_noDirs, AT_, "nghbr");
  MInt counter;
  MFloatScratchSpace dx(nDim, AT_, "dx");
  MFloat eps = 0.000001;
  MFloatScratchSpace res(m_noSets, AT_, "res");
  MFloat smoothingTerm = F0;
  MFloat sumOfD;
  MFloat sumOfPhi;
  MIntScratchSpace reinit(m_maxNoCells * m_noSets, AT_, "reinit");

#ifdef REINITIALIZATION_STATISTICS

  MFloat avgGradient, maxGradient, minGradient, meanGradient, temp;
  MInt nghbrId;
  MFloatScratchSpace x(a_noBandCells(startSet) * m_noDirs, nDim, AT_, "x");
  MFloat variance;
  MFloat factor;

#endif
  //--- end of initialization

  // compute statistics
#ifdef REINITIALIZATION_STATISTICS
  counter = 0;
  for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
    cellId = a_G0CellId(id, startSet);

    if(a_levelSetFunctionG(cellId, startSet) < F0) {
      for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
        if(a_hasNeighbor(cellId, dirId) > 0) {
          nghbrId = c_neighborId(cellId, dirId);
          if(a_isGZeroCell(nghbrId, startSet) && a_levelSetFunctionG(nghbrId, startSet) > F0) {
            factor = ABS(a_levelSetFunctionG(cellId, startSet))
                     / (ABS(a_levelSetFunctionG(nghbrId, startSet)) + ABS(a_levelSetFunctionG(cellId, startSet)));
            for(MInt i = 0; i < nDim; i++) {
              x(counter, i) = c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i));
            }
            counter++;
          }
        }
      }
    }
  }

  // compute the reinitialization sensor
  avgGradient = F0, maxGradient = F0, minGradient = 1000.0, meanGradient = 0;
  variance = F0;

  for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
    cellId = a_G0CellId(id, startSet);

    temp = F0;
    for(MInt i = 0; i < nDim; i++) {
      temp += POW2(a_levelSetFunctionSlope(cellId, i, startSet));
    }
    temp = ABS(sqrt(temp));

    avgGradient += ABS(temp - F1);
    meanGradient += temp;
    maxGradient = mMax(maxGradient, temp);
    minGradient = mMin(minGradient, temp);
    variance += POW2(temp);
  }
  avgGradient = avgGradient / (MFloat)a_noG0Cells(startSet);
  meanGradient = meanGradient / (MFloat)a_noG0Cells(startSet);
  variance = sqrt(variance) / (MFloat)a_noG0Cells(startSet);
  FILE* avg;
  avg = fopen("avgGradient", "a+");
  fprintf(avg, "%d", globalTimeStep);
  fprintf(avg, "  %f", avgGradient * 1000.0);
  fprintf(avg, "  %f", variance);
  fprintf(avg, "  %f", meanGradient);
  fprintf(avg, "  %f", maxGradient);
  fprintf(avg, "  %f", minGradient);
  fprintf(avg, "\n");
  fclose(avg);
#endif


  // reset the reinit flag
  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt cell = 0; cell < a_internalBandLayer(0, set); cell++)
      reinit[IDX_LSSET(cell, set)] = 1;
  }

  // determine the reinit property
  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
      cellId = a_internalBandCellId(id, set);
      for(MInt j = 0; j < m_noDirs; j++) {
        nghbr[j] = a_bandNghbrIdsG(cellId, j, set);
      }
      for(MInt i = 0; i < nDim; i++) {
        if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(nghbr[2 * i], set) > F0
           && a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(a_bandNghbrIdsG(nghbr[2 * i], 2 * i, set), set)
                  > F0) {
          continue;
        }
        if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(nghbr[2 * i + 1], set) > F0
           && a_levelSetFunctionG(cellId, set)
                      * a_levelSetFunctionG(a_bandNghbrIdsG(nghbr[2 * i + 1], 2 * i + 1, set), set)
                  > F0) {
          continue;
        }
        reinit[IDX_LSSET(id, set)] = 0;
        i = nDim;
      }
    }
  }

  // determine the signed distance of phi0 cells
  // -------------------------------------------
  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
      cellId = a_internalBandCellId(id, set);
      if(reinit[IDX_LSSET(id, set)] == 0) {
        m_d[IDX_LSSET(cellId, set)] = a_levelSetFunctionG(cellId, set);
        continue;
      }

      for(MInt j = 0; j < m_noDirs; j++) {
        nghbr[j] = a_bandNghbrIdsG(cellId, j, set);
      }

      for(MInt i = 0; i < nDim; i++) {
        upwind = false;
        dx[i] = F2 * m_gCellDistance;
        if(!a_isGZeroCell(nghbr[2 * i], set)) {
          nghbr[2 * i] = cellId;
          upwind = true;
          dx[i] -= m_gCellDistance;
        }
        if(!a_isGZeroCell(nghbr[2 * i + 1], set)) {
          nghbr[2 * i + 1] = cellId;
          upwind = true;
          dx[i] -= m_gCellDistance;
        }

        if(!upwind) {
          if(a_levelSetFunctionG(nghbr[2 * i], set) * a_levelSetFunctionG(nghbr[2 * i + 1], set) < F0) {
            if((a_levelSetFunctionG(nghbr[2 * i], set) - a_levelSetFunctionG(cellId, set))
                       * (a_levelSetFunctionG(nghbr[2 * i + 1], set) - a_levelSetFunctionG(cellId, set))
                   > F0
               || a_levelSetFunctionG(nghbr[2 * i], set)
                          * a_levelSetFunctionG(a_bandNghbrIdsG(nghbr[2 * i], 2 * i, set), set)
                      < F0
               || a_levelSetFunctionG(nghbr[2 * i + 1], set)
                          * a_levelSetFunctionG(a_bandNghbrIdsG(nghbr[2 * i + 1], 2 * i + 1, set), set)
                      < F0) {
              if(fabs(a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(nghbr[2 * i], set))
                 > fabs(a_levelSetFunctionG(nghbr[2 * i + 1], set) - a_levelSetFunctionG(cellId, set) + eps)) {
                nghbr[2 * i + 1] = cellId;
                dx[i] -= m_gCellDistance;
              }
              if(fabs(a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(nghbr[2 * i], set) + eps)
                 < fabs(a_levelSetFunctionG(nghbr[2 * i + 1], set) - a_levelSetFunctionG(cellId, set))) {
                nghbr[2 * i] = cellId;
                dx[i] -= m_gCellDistance;
              }
            }
          }
        }
      }

      // compute the sdf at the front
      m_d[IDX_LSSET(cellId, set)] = F0;
      for(MInt i = 0; i < nDim; i++) {
        dx[i] = mMax(eps, dx[i]);
        m_d[IDX_LSSET(cellId, set)] +=
            POW2((a_levelSetFunctionG(nghbr[2 * i + 1], set) - a_levelSetFunctionG(nghbr[2 * i], set)) / dx[i]);
      }
      m_d[IDX_LSSET(cellId, set)] = F1 / sqrt(m_d[IDX_LSSET(cellId, set)]);
      m_d[IDX_LSSET(cellId, set)] *= a_levelSetFunctionG(cellId, set);
    }
  }

  // parallel version: exchange d
  exchangeLs(m_d, 0, m_maxNoSets);

  // determine r for each cell and space direction (stored in m_phiRatio)
  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
      cellId = a_internalBandCellId(id, set);
      for(MInt j = 0; j < m_noDirs; j++) {
        m_phiRatioCells[cellId][IDX_LSSET(j, set)] = -1;
      }
      for(MInt j = 0; j < m_noDirs; j++) {
        nghbr[j] = a_bandNghbrIdsG(cellId, j, set);
        if(a_levelSetFunctionG(nghbr[j], set) * a_levelSetFunctionG(cellId, set) < F0) {
          m_phiRatioCells[cellId][IDX_LSSET(j, set)] = nghbr[j];
          m_phiRatio[cellId][IDX_LSSET(j, set)] = a_levelSetFunctionG(cellId, set) / a_levelSetFunctionG(nghbr[j], set);
        }
      }
    }
  }

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // reinitialization scheme CR-1
    // ----------------------------
    if(methodId == 1) {
      for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
        cellId = a_internalBandCellId(id, set);
        if((a_curvatureG(cellId, set) >= F0 && a_levelSetFunctionG(cellId, set) < F0)
           || (a_curvatureG(cellId, set) < F0 && a_levelSetFunctionG(cellId, set) > F0)) {
          m_correction[IDX_LSSET(cellId, set)] = F0;
          counter = 0;
          for(MInt i = 0; i < m_noDirs; i++) {
            if(m_phiRatioCells[cellId][IDX_LSSET(i, set)] != -1) {
              m_correction[IDX_LSSET(cellId, set)] += m_d[IDX_LSSET(m_phiRatioCells[cellId][IDX_LSSET(i, set)], set)]
                                                      * m_phiRatio[cellId][IDX_LSSET(i, set)];
              counter++;
            }
          }
          if(counter > 0)
            m_correction[IDX_LSSET(cellId, set)] = m_correction[IDX_LSSET(cellId, set)] / (MFloat)counter;
          else
            m_correction[IDX_LSSET(cellId, set)] = m_d[IDX_LSSET(cellId, set)];
        } else {
          m_correction[IDX_LSSET(cellId, set)] = m_d[IDX_LSSET(cellId, set)];
        }
      }
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        m_d[IDX_LSSET(cellId, set)] = m_correction[IDX_LSSET(cellId, set)];
      }
    }

    // reinitialization scheme CR-2
    // ----------------------------
    if(methodId == 2) {
      for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
        cellId = a_internalBandCellId(id, set);

        if((a_curvatureG(cellId, set) >= F0 && a_levelSetFunctionG(cellId, set) < F0)
           || (a_curvatureG(cellId, set) < F0 && a_levelSetFunctionG(cellId, set) > F0)) {
          counter = 0;
          sumOfD = 0;
          sumOfPhi = 0;
          for(MInt i = 0; i < m_noDirs; i++) {
            if(m_phiRatioCells[cellId][IDX_LSSET(i, set)] != -1) {
              sumOfD += m_d[IDX_LSSET(m_phiRatioCells[cellId][IDX_LSSET(i, set)], set)];
              sumOfPhi += a_levelSetFunctionG(m_phiRatioCells[cellId][IDX_LSSET(i, set)], set);
              counter++;
            }
          }
          if(counter > 0)
            m_correction[IDX_LSSET(cellId, set)] = a_levelSetFunctionG(cellId, set) * sumOfD / sumOfPhi;
          else
            m_correction[IDX_LSSET(cellId, set)] = m_d[IDX_LSSET(cellId, set)];
        } else {
          m_correction[IDX_LSSET(cellId, set)] = m_d[IDX_LSSET(cellId, set)];
        }
      }
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        m_d[IDX_LSSET(cellId, set)] = m_correction[IDX_LSSET(cellId, set)];
      }
    }

    // set phi_0 = d
    for(MInt id = 0; id < a_internalBandLayer(0, set); id++)
      if(reinit[IDX_LSSET(id, set)] != 0)
        a_levelSetFunctionG(a_internalBandCellId(id, set), set) +=
            m_omegaReinit
            * (m_d[IDX_LSSET(a_internalBandCellId(id, set), set)]
               - a_levelSetFunctionG(a_internalBandCellId(id, set), set));

    // compute sign(G) and a_hasPositiveSign(G)
    for(MInt id = 0; id < a_noInternalBandCells(set); id++) {
      cellId = a_internalBandCellId(id, set);
      m_signG[IDX_LSSET(cellId, set)] =
          a_levelSetFunctionG(cellId, set) / sqrt(POW2(a_levelSetFunctionG(cellId, set)) + POW2(smoothingTerm));
      if(a_levelSetFunctionG(cellId, set) > F0) {
        a_hasPositiveSign(cellId, set) = true;
      } else {
        a_hasPositiveSign(cellId, set) = false;
      }
    }

    exchangeLevelSet();

    // call the Eikonal solver
    // -----------------------
    // fill the cell list and levelSetFunction
    cellListSize = 0;
    for(MInt id = 0; id < a_noInternalBandCells(set); id++) {
      cellId = a_internalBandCellId(id, set);

      if(gapMode == 0) {
        if(!a_isGZeroCell(cellId, set) && a_nearGapG(cellId) > 0 && a_potentialGapCellClose(cellId)) {
          m_cellList[cellListSize] = cellId;
          cellListSize++;
        }
      } else {
        if(!a_isGZeroCell(cellId, set)) {
          m_cellList[cellListSize] = cellId;
          cellListSize++;
        }
      }
    }

    RECORD_TIMER_STOP(t_c1);
    RECORD_TIMER_START(t_c2);
    res[set] = firstOrderEikonalSolver(cellListSize, m_gReinitIterations, set);
    RECORD_TIMER_STOP(t_c2);
    RECORD_TIMER_START(t_c1);
#ifndef NDEBUG
    m_log << "Reinitialization finished at ts " << globalTimeStep << ": " << res[set] << endl;
#endif
  }

  // statistics
#ifdef REINITIALIZATION_STATISTICS
  counter = 0;
  for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
    cellId = a_G0CellId(id, startSet);
    if(a_levelSetFunctionG(cellId, startSet) < F0) {
      for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
        if(a_hasNeighbor(cellId, dirId) > 0) {
          nghbrId = c_neighborId(cellId, dirId);
          if(a_isGZeroCell(nghbrId, startSet) && a_levelSetFunctionG(nghbrId, startSet) > F0) {
            factor = ABS(a_levelSetFunctionG(cellId, startSet))
                     / (ABS(a_levelSetFunctionG(nghbrId, startSet)) + ABS(a_levelSetFunctionG(cellId, startSet)));
            for(MInt i = 0; i < nDim; i++) {
              x(counter, i) =
                  100
                  * (x(counter, i)
                     - (c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i))));
            }
            counter++;
          }
        }
      }
    }
  }
  MFloat deviation = F0;
  for(MInt it = 0; it < counter; it++) {
    factor = F0;
    for(MInt i = 0; i < nDim; i++)
      factor += POW2(x(it, i));
    deviation += sqrt(factor);
  }
  deviation = deviation * 1000.0 / (MFloat)counter;
  FILE* dev;
  dev = fopen("deviation", "a+");
  fprintf(dev, "%d", globalTimeStep);
  fprintf(dev, "  %f", deviation);
  fprintf(dev, "\n");
  fclose(dev);

  //
  computeNormalVectors();

  avgGradient = F0, maxGradient = F0, minGradient = 1000.0;
  variance = F0;

  MInt nghbrL, nghbrL2, nghbrR, nghbrR2;
  for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
    cellId = a_G0CellId(id, startSet);

    temp = F0;
    for(MInt i = 0; i < nDim; i++) {
      // compute the fourth-order gradient
      nghbrL = a_bandNghbrIdsG(cellId, 2 * i, startSet);
      nghbrL2 = a_bandNghbrIdsG(nghbrL, 2 * i, startSet);
      nghbrR = a_bandNghbrIdsG(cellId, 2 * i + 1, startSet);
      nghbrR2 = a_bandNghbrIdsG(nghbrR, 2 * i + 1, startSet);
      a_levelSetFunctionSlope(cellId, i, startSet) =
          m_FgCellDistance
          * (F2B3 * (a_levelSetFunctionG(nghbrR, startSet) - a_levelSetFunctionG(nghbrL, startSet))
             - F1B12 * (a_levelSetFunctionG(nghbrR2, startSet) - a_levelSetFunctionG(nghbrL2, startSet)));
      temp += POW2(a_levelSetFunctionSlope(cellId, i, startSet));
    }
    temp = ABS(sqrt(temp));

    avgGradient += ABS(temp - F1);
    meanGradient += temp;
    maxGradient = mMax(maxGradient, temp);
    minGradient = mMin(minGradient, temp);
    variance += POW2(temp);
  }
  avgGradient = avgGradient / (MFloat)a_noG0Cells(startSet);
  meanGradient = meanGradient / (MFloat)a_noG0Cells(startSet);
  variance = sqrt(variance) / (MFloat)a_noG0Cells(startSet);
  FILE* avg2;
  avg2 = fopen("avgGradientAfter", "a+");
  fprintf(avg2, "%d", globalTimeStep);
  fprintf(avg2, "  %f", avgGradient * 1000.0);
  fprintf(avg2, "  %f", variance);
  fprintf(avg2, "  %f", meanGradient);
  fprintf(avg2, "  %f", maxGradient);
  fprintf(avg2, "  %f", minGradient);
  fprintf(avg2, "\n");
  fclose(avg2);

#endif

  RECORD_TIMER_STOP(t_c1);
  RECORD_TIMER_STOP(t_reInit);
}
//-----------------------------------------------------------------------------

/** \brief level set high order constrained reinitialization
 *
 * \author Daniel Hartmann (probably), Stephan Schlimpert
 * \date unknown, June 2011
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetHighOrderConstrainedReinitialization(MInt methodId, MInt startSet, MInt endSet,
                                                                           MInt gapMode) {
  TRACE();

  MInt cellId;
  MInt cellListSize;
  MInt counter;
  MInt counter2;
  MFloat smoothingTerm = m_gCellDistance;
  MFloat sumOfPhi;


  MFloatScratchSpace res(m_noSets, AT_, "res");
  MInt noCellsToCorrect = 0;
  MIntScratchSpace cellsToCorrect(m_maxNoCells, AT_, "cellsToCorrect");
  MFloatScratchSpace factors(m_maxNoCells, AT_, "factors");

  // output
  MInt nghbrId;
  MFloat avgGradient;
  MFloat maxGradient;
  MFloat minGradient;
  MFloat meanGradient = F0;
  MFloat factor;
  MFloat temp;
  MFloat variance;
  MFloat** x;
  //--- end of initialization

  // statistics output for startSet only, no statistics for other sets given!
  if(!m_writeReinitializationStatistics) {
    x = (MFloat**)nullptr;
  } else {
    x = new MFloat*[a_noBandCells(startSet) * m_noDirs];
  }
  if(m_writeReinitializationStatistics) {
    for(MInt cell = 0; cell < a_noBandCells(startSet) * m_noDirs; cell++)
      x[cell] = new MFloat[nDim];
    // FOR PUBLICATION OUTPUT ONLY!!!
    // determine the original position of the interface
    counter = 0;
    for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
      cellId = a_G0CellId(id, startSet);

      if(a_levelSetFunctionG(cellId, startSet) < F0) {
        for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
          if(a_hasNeighbor(cellId, dirId) > 0) {
            nghbrId = c_neighborId(cellId, dirId);
            if(a_isGZeroCell(nghbrId, startSet) && a_levelSetFunctionG(nghbrId, startSet) > F0) {
              factor = ABS(a_levelSetFunctionG(cellId, startSet))
                       / (ABS(a_levelSetFunctionG(nghbrId, startSet)) + ABS(a_levelSetFunctionG(cellId, startSet)));
              for(MInt i = 0; i < nDim; i++) {
                x[counter][i] = c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i));
              }
              counter++;
            }
          }
        }
      }
    }
  }

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // compute sign(G) and a_hasPositiveSign(G)
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      m_signG[IDX_LSSET(cellId, set)] =
          a_levelSetFunctionG(cellId, set) / sqrt(POW2(a_levelSetFunctionG(cellId, set)) + POW2(smoothingTerm));
      if(a_levelSetFunctionG(cellId, set) > F0) {
        a_hasPositiveSign(cellId, set) = true;
      } else {
        a_hasPositiveSign(cellId, set) = false;
      }
    }


    // fill the cell list
    cellListSize = 0;
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      if(gapMode == 0) {
        if(!a_isGZeroCell(cellId, set) && a_nearGapG(cellId) > 0 && a_potentialGapCellClose(cellId)) {
          m_cellList[cellListSize] = cellId;
          cellListSize++;
        }
      } else {
        m_cellList[cellListSize] = cellId;
        cellListSize++;
      }
    }

    // determine r for each cell and space direction (stored in m_phiRatio)
    for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
      cellId = a_internalBandCellId(id, set);
      for(MInt j = 0; j < m_noDirs; j++)
        m_phiRatioCells[cellId][IDX_LSSET(j, set)] = -1;
      for(MInt j = 0; j < m_noDirs; j++) {
        if(a_levelSetFunctionG(a_bandNghbrIdsG(cellId, j, set), set) * a_levelSetFunctionG(cellId, set) < F0) {
          m_phiRatioCells[cellId][IDX_LSSET(j, set)] = a_bandNghbrIdsG(cellId, j, set);
          m_phiRatio[cellId][IDX_LSSET(j, set)] =
              a_levelSetFunctionG(cellId, set) / a_levelSetFunctionG(a_bandNghbrIdsG(cellId, j, set), set);
        }
      }
    }

    // compute the forcing terms and call the Eikonal solver

    // reinitialization scheme CR-1
    // ----------------------------
    if(methodId == 1) {
      counter2 = 0;
      noCellsToCorrect = 0;
      for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
        cellId = a_internalBandCellId(id, set);
        m_correction[IDX_LSSET(cellId, set)] = F0;
        counter = 0;
        for(MInt i = 0; i < m_noDirs; i++) {
          if(m_phiRatioCells[cellId][IDX_LSSET(i, set)] != -1) {
            factors.p[counter2++] =
                a_levelSetFunctionG(cellId, set) / a_levelSetFunctionG(m_phiRatioCells[cellId][IDX_LSSET(i, set)], set);
            counter++;
          }
        }
        if(counter > 0) cellsToCorrect.p[noCellsToCorrect++] = cellId;
      }

      res[set] = fifthOrderEikonalSolver(cellListSize, m_gReinitIterations, cellsToCorrect.getPointer(),
                                         noCellsToCorrect, factors.getPointer(), 1, set);
    }

    // reinitialization scheme CR-2
    // ----------------------------
    if(methodId == 2 || methodId == 4 || methodId == 6) {
      noCellsToCorrect = 0;
      for(MInt id = 0; id < a_internalBandLayer(0, set); id++) {
        cellId = a_internalBandCellId(id, set);
        counter = 0;
        sumOfPhi = 0;
        for(MInt i = 0; i < m_noDirs; i++) {
          if(m_phiRatioCells[cellId][IDX_LSSET(i, set)] != -1) {
            sumOfPhi += a_levelSetFunctionG(m_phiRatioCells[cellId][IDX_LSSET(i, set)], set);
            counter++;
          }
        }
        factors.p[noCellsToCorrect] = a_levelSetFunctionG(cellId, set) / sumOfPhi;
        if(counter > 0) cellsToCorrect.p[noCellsToCorrect++] = cellId;
      }
      if(methodId == 2)
        res[set] = fifthOrderEikonalSolver(cellListSize, m_gReinitIterations, cellsToCorrect.getPointer(),
                                           noCellsToCorrect, factors.getPointer(), 2, set);
      else if(methodId == 4)
        res[set] = fifthOrderEikonalSolver(cellListSize, m_gReinitIterations, cellsToCorrect.getPointer(),
                                           noCellsToCorrect, factors.getPointer(), 4, set);
      else if(methodId == 6)
        res[set] = fifthOrderEikonalSolver(cellListSize, m_gReinitIterations, cellsToCorrect.getPointer(),
                                           noCellsToCorrect, factors.getPointer(), 6, set);
    }
  }

  if(domainId() == 0) {
    // write out info
    FILE* datei = nullptr;
    stringstream reinitFile;
    reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
    datei = fopen((reinitFile.str()).c_str(), "a+");
    fprintf(datei, "  %-10.8f  reinitialized", res[startSet]);
    fclose(datei);
  }

  if(m_writeReinitializationStatistics) {
    // FOR PUBLICATION OUTPUT ONLY!!!
    counter = 0;
    for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
      cellId = a_G0CellId(id, startSet);
      if(a_levelSetFunctionG(cellId, startSet) < F0) {
        for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
          if(a_hasNeighbor(cellId, dirId) > 0) {
            nghbrId = c_neighborId(cellId, dirId);
            if(a_isGZeroCell(nghbrId, startSet) && a_levelSetFunctionG(nghbrId, startSet) > F0) {
              factor = ABS(a_levelSetFunctionG(cellId, startSet))
                       / (ABS(a_levelSetFunctionG(nghbrId, startSet)) + ABS(a_levelSetFunctionG(cellId, startSet)));
              for(MInt i = 0; i < nDim; i++) {
                x[counter][i] =
                    100
                    * (x[counter][i]
                       - (c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i))));
              }
              counter++;
            }
          }
        }
      }
    }
    MFloat deviation = F0;
    for(MInt it = 0; it < counter; it++) {
      factor = F0;
      for(MInt i = 0; i < nDim; i++)
        factor += POW2(x[it][i]);
      deviation += sqrt(factor);
    }
    deviation = deviation * 1000.0 / (MFloat)counter;
    FILE* dev = nullptr;
    dev = fopen("deviation", "a+");
    fprintf(dev, "%d", globalTimeStep);
    fprintf(dev, "  %f", deviation);
    fprintf(dev, "\n");
    fclose(dev);

    avgGradient = F0, maxGradient = F0, minGradient = 1000.0;
    variance = F0;

    MInt nghbrL;
    MInt nghbrL2;
    MInt nghbrR;
    MInt nghbrR2;
    for(MInt id = 0; id < a_noG0Cells(startSet); id++) {
      cellId = a_G0CellId(id, startSet);

      temp = F0;
      for(MInt i = 0; i < nDim; i++) {
        // compute the fourth-order gradient
        nghbrL = a_bandNghbrIdsG(cellId, 2 * i, startSet);
        nghbrL2 = a_bandNghbrIdsG(nghbrL, 2 * i, startSet);
        nghbrR = a_bandNghbrIdsG(cellId, 2 * i + 1, startSet);
        nghbrR2 = a_bandNghbrIdsG(nghbrR, 2 * i + 1, startSet);
        a_levelSetFunctionSlope(cellId, i, startSet) =
            m_FgCellDistance
            * (F2B3 * (a_levelSetFunctionG(nghbrR, startSet) - a_levelSetFunctionG(nghbrL, startSet))
               - F1B12 * (a_levelSetFunctionG(nghbrR2, startSet) - a_levelSetFunctionG(nghbrL2, startSet)));
        temp += POW2(a_levelSetFunctionSlope(cellId, i, startSet));
      }
      temp = ABS(sqrt(temp));

      avgGradient += ABS(temp - F1);
      meanGradient += temp;
      maxGradient = mMax(maxGradient, temp);
      minGradient = mMin(minGradient, temp);
      variance += POW2(temp);
    }
    avgGradient = avgGradient / (MFloat)a_noG0Cells(startSet);
    meanGradient = meanGradient / (MFloat)a_noG0Cells(startSet);
    variance = sqrt(variance) / (MFloat)a_noG0Cells(startSet);
    FILE* avg2;
    avg2 = fopen("avgGradientAfter", "a+");
    fprintf(avg2, "%d", globalTimeStep);
    fprintf(avg2, "  %f", avgGradient * 1000.0);
    fprintf(avg2, "  %f", variance);
    fprintf(avg2, "  %f", meanGradient);
    fprintf(avg2, "  %f", maxGradient);
    fprintf(avg2, "  %f", minGradient);
    fprintf(avg2, "\n");
    fclose(avg2);

    // free memory
    for(MInt cell = 0; cell < a_noBandCells(startSet) * m_noDirs; cell++)
      delete[] x[cell];
    delete[] x;
    x = nullptr;
  }
}


//-----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::maintainOuterBandLayers(MInt order, MInt startSet, MInt endSet) {
  TRACE();

  MInt cellId;
  MInt cellListSize;
  MFloat smoothingTerm = m_gCellDistance;
  MFloatScratchSpace res(m_noSets, AT_, "res");
  auto startBand = (MInt)(m_gBandWidth / 2);
  if(startBand > 6) startBand = 6;
  //--- end of initialization

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // compute sign(G) and a_hasPositiveSign(G)
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      m_signG[IDX_LSSET(cellId, set)] =
          a_levelSetFunctionG(cellId, set) / sqrt(POW2(a_levelSetFunctionG(cellId, set)) + POW2(smoothingTerm));
      if(a_levelSetFunctionG(cellId, set) > F0) {
        a_hasPositiveSign(cellId, set) = true;
      } else {
        a_hasPositiveSign(cellId, set) = false;
      }
    }

    // fill the cell list
    cellListSize = 0;
    for(MInt id = a_bandLayer(startBand, set); id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      m_cellList[cellListSize] = cellId;
      cellListSize++;
    }

    switch(order) {
      case 1:
        res[set] = firstOrderEikonalSolver(cellListSize, m_maintenanceIterations, set);
        break;
      case 5:
        res[set] =
            fifthOrderEikonalSolver(cellListSize, m_maintenanceIterations, (MInt*)nullptr, 0, (MFloat*)nullptr, 0, set);
        break;
      default: {
        stringstream errorMessage;
        errorMessage << "LsCartesianSolver::maintainOuterBandLayers(): switch variable 'order' with value " << order
                     << " not matching any case." << endl;
        mTerm(1, AT_, errorMessage.str());
      }
    }
  }
  // write out info
  FILE* datei;
  stringstream reinitFile;
  reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
  datei = fopen((reinitFile.str()).c_str(), "a+");
  fprintf(datei, "  %f  maintained", res[startSet]);
  fclose(datei);
}


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
/** \brief computes min, max and mean flame front position in a small region
 *
 * \author Daniel Hartmann, Stephan Schlimpert
 * \date unkown, February 2013
 *
 * last changes: parallelization
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::determineMinMaxMeanInterfacePosition() {
  TRACE();
  MInt counter;
  MInt cellId;
  MInt nghbrId;
  MFloat Fcounter;
  MFloat factor;
  MInt set = 0; // operates on zeroth level-set function only!
  MFloatScratchSpace x(a_noBandCells(set) * nDim, AT_, "x");

  //---

  // determine all the G0 points
  counter = 0;
  for(MInt id = 0; id < a_noG0Cells(set); id++) {
    cellId = a_G0CellId(id, set);
    if(a_isHalo(cellId)) continue;

    if(a_levelSetFunctionG(cellId, set) < F0) {
      for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
        if(a_hasNeighbor(cellId, dirId) > 0) {
          nghbrId = c_neighborId(cellId, dirId);
          if(a_isGZeroCell(nghbrId, set) && a_levelSetFunctionG(nghbrId, set) > F0) {
            factor = ABS(a_levelSetFunctionG(cellId, set))
                     / (ABS(a_levelSetFunctionG(nghbrId, set)) + ABS(a_levelSetFunctionG(cellId, set)));
            for(MInt i = 0; i < nDim; i++) {
              x.p[counter * nDim + i] =
                  c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i));
            }
            counter++;
          }
        }
      }
    }
  }

  // determine the min, max, and mean values
  Fcounter = F1 / (MFloat)counter;
  for(MInt i = 0; i < nDim; i++) {
    m_minFlameFrontPosition[i] = 10000.0;
    m_maxFlameFrontPosition[i] = -10000.0;
    m_meanFlameFrontPosition[i] = F0;
  }
  for(MInt p = 0; p < counter; p++) {
    for(MInt i = 0; i < nDim; i++) {
      m_minFlameFrontPosition[i] = mMin(m_minFlameFrontPosition[i], x.p[nDim * p + i]);
      m_maxFlameFrontPosition[i] = mMax(m_maxFlameFrontPosition[i], x.p[nDim * p + i]);
      m_meanFlameFrontPosition[i] += x.p[nDim * p + i] * Fcounter;
    }
  }
  // exchange min/max/mean flame front position
  for(MInt i = 0; i < nDim; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &m_minFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_minFlameFrontPosition[i]");
    MPI_Allreduce(MPI_IN_PLACE, &m_maxFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_maxFlameFrontPosition[i]");
    MPI_Allreduce(MPI_IN_PLACE, &m_meanFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_meanFlameFrontPosition[i]");
  }
}

//-----------------------------------------------------------------------------
/** \brief computes flame base angle of a steady flame surface
 *
 * \author Stephan Schlimpert
 * \date April 2012
 *
 * needed for determining flame front amplitude of a forced flame
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::determineSteadyFlameLength() {
  TRACE();

  if(m_steadyFlameLength > -F1) {
    m_steadyFlameAngle = -F1;
    m_log << "WARINING: steady flame front length should be defined correctly in your properties.cdl " << endl;
    m_steadyFlameAngle = atan(m_steadyFlameLength / m_realRadiusFlameTube) * 180. / PI;

    m_log << "steadyFlameLength is : " << m_steadyFlameLength << endl;
    m_log << "uncurved flame base angle in Grad: " << m_steadyFlameAngle << endl;
    m_log << "flame surface slope is: " << tan(m_steadyFlameAngle * PI / 180.0) << endl;
  }

  if(approx(m_steadyFlameLength, -F1, MFloatEps) && m_forcing) {
    MString errorMessage = "ERROR: steady flame length is not determined!!! should be defined for forced flames";
    mTerm(1, AT_, errorMessage);
  }
}

//-----------------------------------------------------------------------------
/** \brief computes min, max and mean flame front position in a small region
 *
 * \author Stephan Schlimpert
 * \date July 2011, February 2013
 *
 * last changes: parallelization
 *
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::determineMinMaxMeanRegionInterfacePosition(MFloat xRegN, MFloat xRegP, MFloat yRegN,
                                                                         MFloat yRegP, MInt set) {
  TRACE();

  MInt counter;
  MInt cellId;
  MInt nghbrId;
  MFloat Fcounter;
  MFloat factor;
  MFloatScratchSpace x(a_noBandCells(set) * nDim, AT_, "x");
  //---

  // determine all the G0 points
  counter = 0;
  for(MInt id = 0; id < a_noG0Cells(set); id++) {
    cellId = a_G0CellId(id, set);

    if(c_coordinate(cellId, 0) > xRegP) continue;
    if(c_coordinate(cellId, 0) < xRegN) continue;
    if(c_coordinate(cellId, 1) > yRegP) continue;
    if(c_coordinate(cellId, 1) < yRegN) continue;
    if(a_levelSetFunctionG(cellId, set) < F0) {
      for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
        if(a_hasNeighbor(cellId, dirId) > 0) {
          nghbrId = c_neighborId(cellId, dirId);
          if(a_isGZeroCell(nghbrId, set) && a_levelSetFunctionG(nghbrId, set) > F0) {
            factor = ABS(a_levelSetFunctionG(cellId, set))
                     / (ABS(a_levelSetFunctionG(nghbrId, set)) + ABS(a_levelSetFunctionG(cellId, set)));
            for(MInt i = 0; i < nDim; i++) {
              x.p[counter * nDim + i] =
                  c_coordinate(cellId, i) + factor * (c_coordinate(nghbrId, i) - c_coordinate(cellId, i));
            }
            counter++;
          }
        }
      }
    }
  }

  // determine the min, max, and mean values
  Fcounter = F1 / (MFloat)counter;
  for(MInt i = 0; i < nDim; i++) {
    m_minFlameFrontPosition[i] = 10000.0;
    m_maxFlameFrontPosition[i] = -10000.0;
    m_meanFlameFrontPosition[i] = F0;
  }
  for(MInt p = 0; p < counter; p++) {
    for(MInt i = 0; i < nDim; i++) {
      m_minFlameFrontPosition[i] = mMin(m_minFlameFrontPosition[i], x.p[nDim * p + i]);
      m_maxFlameFrontPosition[i] = mMax(m_maxFlameFrontPosition[i], x.p[nDim * p + i]);
      m_meanFlameFrontPosition[i] += x.p[nDim * p + i] * Fcounter;
    }
  }

  // exchange min/max/mean flame front position
  for(MInt i = 0; i < nDim; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &m_minFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_minFlameFrontPosition[i]");
    MPI_Allreduce(MPI_IN_PLACE, &m_maxFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_maxFlameFrontPosition[i]");
    MPI_Allreduce(MPI_IN_PLACE, &m_meanFlameFrontPosition[i], 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_meanFlameFrontPosition[i]");
  }
}


//-----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::reinitBand(MInt startSet, MInt endSet) {
  TRACE();

  MInt cellId;
  MInt cellListSize;
  MFloatScratchSpace res(m_noSets, AT_, "res");
  //---

  // based on Sussman, JCP 1994

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // compute sign(G) and a_hasPositiveSign(G)
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      m_signG[IDX_LSSET(cellId, set)] = a_levelSetFunctionG(cellId, set) / ABS(a_levelSetFunctionG(cellId, set));
      if(a_levelSetFunctionG(cellId, set) > F0) {
        a_hasPositiveSign(cellId, set) = true;
      } else {
        a_hasPositiveSign(cellId, set) = false;
      }
    }

    // call the Eikonal solver
    // -----------------------
    // fill the cell list and levelSetFunction
    cellListSize = 0;
    for(MInt id = a_bandLayer(0, set); id < a_noBandCells(set); id++) {
      m_cellList[cellListSize] = a_bandCellId(id, set);
      cellListSize++;
    }

    if(cellListSize > 0) {
      res[set] = firstOrderEikonalSolver(cellListSize, m_intermediateReinitIterations, set);
      m_log << "Reinitialization at ts " << globalTimeStep << ": " << res[set] << endl;
    } else {
      m_log << "Reinitialization skipped since no cell was found to reinitialize " << endl;
    }
  }
}


// ----------------------------------------------------------------------------------------


/** \brief Solves the equation  0 = sgn(phi) ( 1 - | nabla phi | ) on a uniform Cartesian grid
 *
 * \author Daniel Hartmann
 * \date April 2007
 *
 * Input:  q            - level set function
 *         m_signG      - smoothed sign of the level set function
 *         cellList     - list cells to reinitialize
 *         cellListSize - number of cells to reinitialize
 * Output: q            - reinitialized level set function
 *
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::firstOrderEikonalSolver(MInt cellListSize, MInt maxIterations, MInt set) {
  TRACE();


  MInt cellId;
  MInt iteration;
  MFloat a, b;
  MFloat G;
  MFloat dt = m_reinitCFL * m_gCellDistance;
  MFloatScratchSpace res(1, AT_, "res");
  MFloatScratchSpace globalRes(1, AT_, "globalRes");
  //---

  // Iteration loop of the Sussman method
  // ------------------------------------
  iteration = 0;
  res.p[0] = F0;

  if(maxIterations == 0) {
    return res.p[0];
  }

  while(iteration < maxIterations) {
    res.p[0] = F0;

    // discretization
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];
      a_levelSetRHS(cellId, set) = F0;
      for(MInt i = 0; i < nDim; i++) {
        a = (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(a_bandNghbrIdsG(cellId, 2 * i, set), set))
            * m_FgCellDistance;
        b = (a_levelSetFunctionG(a_bandNghbrIdsG(cellId, 2 * i + 1, set), set) - a_levelSetFunctionG(cellId, set))
            * m_FgCellDistance;
        a_levelSetRHS(cellId, set) +=
            F1B2
            * ((a_levelSetSign(cellId, set) + F1) * mMax(POW2(mMax(a, F0)), POW2(mMin(b, F0)))
               - (a_levelSetSign(cellId, set) - F1) * mMax(POW2(mMin(a, F0)), POW2(mMax(b, F0))));
      }
      a_levelSetRHS(cellId, set) = m_signG[IDX_LSSET(cellId, set)] * (F1 - sqrt(a_levelSetRHS(cellId, set)));
    }

    // solver
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];
      G = a_levelSetFunctionG(cellId, set) + a_levelSetRHS(cellId, set) * dt;
      if(G * a_levelSetFunctionG(cellId, set) < F0) G = a_levelSetFunctionG(cellId, set);
      a_levelSetFunctionG(cellId, set) = G;
      res.p[0] = mMax(res.p[0], ABS(a_levelSetRHS(cellId, set)));
    }
    iteration++;

    MFloat* q;
    q = (MFloat*)&a_levelSetFunctionG(0, 0);
    exchangeLs(q, set, 1);
    // exchange residual and data
    MPI_Allreduce(res.getPointer(), globalRes.getPointer(), 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "res.getPointer()",
                  "globalRes.getPointer()");

    res.p[0] = globalRes.p[0];

    if(res.p[0] < m_reinitConvergence) {
      return res.p[0];
    }
  }

  return res.p[0];
}


// ----------------------------------------------------------------------------------------

/** \brief Solves the equation  0 = sgn(phi) ( 1 - | nabla phi | ) on a uniform Cartesian grid
 *
 * \author Daniel Hartmann
 * \date April 2007
 *
 * Input:  q            - level set function
 *         m_signG      - smoothed sign of the level set function
 *         cellList     - list cells to reinitialize
 *         cellListSize - number of cells to reinitialize
 * Output: q            - reinitialized level set function
 *
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::secondOrderEikonalSolver(MFloat* q, const MInt* nghbrs, MInt cellListSize,
                                                         MInt maxIterations, MInt set) {
  TRACE();

  MInt cellId;
  MInt iteration;
  MFloat a, b;
  MFloat G;
  MFloat cminus, cplus;
  MFloat PHI[7];
  MFloat dt = m_reinitCFL * m_gCellDistance;
  MFloatScratchSpace res(1, AT_, "res");
  MFloatScratchSpace globalRes(1, AT_, "globalRes");
  //---

  // Iteration loop of the Sussman method
  // ------------------------------------
  iteration = 0;
  res.p[0] = F0;

  if(maxIterations == 0) {
    return res.p[0];
  }

  while(iteration < maxIterations) {
    res.p[0] = F0;

    // discretization
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];
      a_levelSetRHS(cellId, set) = F0;
      for(MInt i = 0; i < nDim; i++) {
        // compute the table of divided differences
        // entries 0-6: PHI(k-2,k-1), PHI(k-1,k), PHI(k+1,k), PHI(k+2,k+1)
        //              PHI(k-2,k), PHI(k-1,k+1), PHI(k,k+2)
        PHI[0] = (q[IDX_LSSET(nghbrs[IDX_LSSETDIR(cellId, 2 * i, set)], set)]
                  - q[IDX_LSSET(nghbrs[IDX_LSSETDIR(nghbrs[IDX_LSSETDIR(cellId, 2 * i, set)], 2 * i, set)], set)])
                 * m_FgCellDistance;
        PHI[1] = (q[IDX_LSSET(cellId, set)] - q[IDX_LSSET(nghbrs[IDX_LSSETDIR(cellId, 2 * i, set)], set)])
                 * m_FgCellDistance;
        PHI[2] = (q[IDX_LSSET(nghbrs[IDX_LSSETDIR(cellId, 2 * i + 1, set)], set)] - q[IDX_LSSET(cellId, set)])
                 * m_FgCellDistance;
        PHI[3] = (q[IDX_LSSET(nghbrs[IDX_LSSETDIR(nghbrs[IDX_LSSETDIR(cellId, 2 * i + 1, set)], 2 * i + 1, set)], set)]
                  - q[IDX_LSSET(cellId, set)])
                 * m_FgCellDistance;
        PHI[4] = (PHI[1] - PHI[0]) * F1B2 * m_FgCellDistance;
        PHI[5] = (PHI[2] - PHI[1]) * F1B2 * m_FgCellDistance;
        PHI[6] = (PHI[3] - PHI[2]) * F1B2 * m_FgCellDistance;
        // compute cminus and cplus using a minmod limiter
        if(PHI[4] * PHI[5] > F0) {
          if(ABS(PHI[4]) <= ABS(PHI[5]))
            cminus = PHI[4];
          else
            cminus = PHI[5];
        } else
          cminus = F0;
        if(PHI[5] * PHI[6] > F0) {
          if(ABS(PHI[5]) <= ABS(PHI[6]))
            cplus = PHI[5];
          else
            cplus = PHI[6];
        } else
          cplus = F0;
        a = PHI[1] + cminus * m_gCellDistance;
        b = PHI[2] - cplus * m_gCellDistance;
        a_levelSetRHS(cellId, set) +=
            F1B2
            * ((a_levelSetSign(cellId, set) + F1) * mMax(POW2(mMax(a, F0)), POW2(mMin(b, F0)))
               - (a_levelSetSign(cellId, set) - F1) * mMax(POW2(mMin(a, F0)), POW2(mMax(b, F0))));
      }
      a_levelSetRHS(cellId, set) = m_signG[IDX_LSSET(cellId, set)] * (F1 - sqrt(a_levelSetRHS(cellId, set)));
    }

    // solver
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];

      G = q[IDX_LSSET(cellId, set)] + a_levelSetRHS(cellId, set) * dt;
      q[IDX_LSSET(cellId, set)] = ABS(G) * a_levelSetSign(cellId, set);
      res.p[0] = mMax(res.p[0], ABS(a_levelSetRHS(cellId, set)));
    }
    iteration++;

    exchangeLs(q, set, 1);
    // exchange residual and data
    MPI_Allreduce(res.getPointer(), globalRes.getPointer(), 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "res.getPointer()",
                  "globalRes.getPointer()");
    res.p[0] = globalRes.p[0];

    if(res.p[0] < m_reinitConvergence) {
      return res.p[0];
    }
  }

  return res.p[0];
}


// ----------------------------------------------------------------------------------------

/** \brief Solves the equation  0 = sgn(phi) ( 1 - | nabla phi | ) on a uniform Cartesian grid
 *
 * \author Daniel Hartmann, Stephan Schlimpert, Claudia Guenther
 * \date April 2007, June 2011, Nov. 2012
 *
 * Input:  q            - level set function
 *         nghbrs       - list of valid neighbors of all cells in cellList
 *         cellList     - list cells to reinitialize
 *         cellListSize - number of cells to reinitialize
 *         maxIterations - maximum number of iterations
 *         crCells      - cells on which forcing for constrained reinitialization should be applied (only in CR-Modes)
 *         noCRCells    - number of cells in crCells (only in CR-Modes)
 *         factors      - factors needed for CR forcing (only in CR-Modes)
 *         crMode       - 0 if no constrained Reinitialization is applied, 1-3: constrained reinitialization modes
 *         set          - which level-set function should be reinitialized
 * Output: q            - reinitialized level set function
 *         res          - residual (maximum right hand side = max. deviation of abs(nabla(phi)) = 1
 *
 * change in June 2011: use of new scratch method and change of firstOrder array use
 * change in Nov 2012: unified all fifthOrderEikonal/CREikonal solvers in one method!
 * change in Jan 2013: method 6 included -> contains minIteration to guarantee fully reinitialization of band
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::fifthOrderEikonalSolver(MInt cellListSize, MInt maxIterations, MInt* crCells,
                                                        MInt noCRCells, MFloat* factors, MInt crMode, MInt set) {
  TRACE();

  MInt cellId;
  MInt iteration;
  MInt nghbrL, nghbrL2, nghbrL3, nghbrL4, nghbrL5, nghbrR, nghbrR2, nghbrR3, nghbrR4, nghbrR5;
  MFloat a, b;
  MFloat dPlusL3, dPlusL2, dPlusL, dPlus, dPlusR, dPlusR2;
  MFloat eps = 0.000001;
  MFloat IS0, IS1, IS2;
  MFloat alpha0, alpha1, alpha2, omega0, omega2;
  MFloat PsiMinus, PsiPlus;
  MFloat c, d;
  MBool converged = false; // only active for crMode == 4 (limited CR2)
  MFloat dt = m_reinitCFL * m_gCellDistance;
  MBoolScratchSpace firstOrder(cellListSize, AT_, "firstOrder");
  MFloat res = F0;
  MInt minIteration;
  if(crMode == 6)
    minIteration = m_gBandWidth / m_reinitCFL;
  else
    minIteration = 0;

  //---

  // Iteration loop of the Sussman method
  // ------------------------------------
  iteration = 0;

  if(!maxIterations) {
    return res;
  }

  // only method 6 resets level set function
  switch(crMode) {
    case 6: {
      // check all cells which have a stencil neighbor outside of the domain
      for(MInt cell = 0; cell < cellListSize; cell++) {
        firstOrder.p[cell] = false;
        cellId = m_cellList[cell];
        if(a_isHalo(cellId)) continue;

        for(MInt i = 0; i < nDim; i++) {
          nghbrL = a_bandNghbrIdsG(cellId, 2 * i, set);
          nghbrL2 = a_bandNghbrIdsG(nghbrL, 2 * i, set);
          nghbrL3 = a_bandNghbrIdsG(nghbrL2, 2 * i, set);
          nghbrL4 = a_bandNghbrIdsG(nghbrL3, 2 * i, set);
          nghbrL5 = a_bandNghbrIdsG(nghbrL4, 2 * i, set);
          nghbrR = a_bandNghbrIdsG(cellId, 2 * i + 1, set);
          nghbrR2 = a_bandNghbrIdsG(nghbrR, 2 * i + 1, set);
          nghbrR3 = a_bandNghbrIdsG(nghbrR2, 2 * i + 1, set);
          nghbrR4 = a_bandNghbrIdsG(nghbrR3, 2 * i + 1, set);
          nghbrR5 = a_bandNghbrIdsG(nghbrR4, 2 * i + 1, set);
          if(nghbrL5 == nghbrL4 || nghbrR5 == nghbrR4) {
            firstOrder.p[cell] = true;
            break;
          }
        }
      }
      break;
    }
    default: {
      // check all cells which have a stencil neighbor outside of the domain
      for(MInt cell = 0; cell < cellListSize; cell++) {
        firstOrder.p[cell] = false;
        cellId = m_cellList[cell];
        if(a_isHalo(cellId)) continue;
        for(MInt i = 0; i < nDim; i++) {
          nghbrL = a_bandNghbrIdsG(cellId, 2 * i, set);
          nghbrL2 = a_bandNghbrIdsG(nghbrL, 2 * i, set);
          nghbrL3 = a_bandNghbrIdsG(nghbrL2, 2 * i, set);
          nghbrL4 = a_bandNghbrIdsG(nghbrL3, 2 * i, set);
          nghbrL5 = a_bandNghbrIdsG(nghbrL4, 2 * i, set);
          nghbrR = a_bandNghbrIdsG(cellId, 2 * i + 1, set);
          nghbrR2 = a_bandNghbrIdsG(nghbrR, 2 * i + 1, set);
          nghbrR3 = a_bandNghbrIdsG(nghbrR2, 2 * i + 1, set);
          nghbrR4 = a_bandNghbrIdsG(nghbrR3, 2 * i + 1, set);
          nghbrR5 = a_bandNghbrIdsG(nghbrR4, 2 * i + 1, set);
          if(nghbrL5 == nghbrL4 || nghbrR5 == nghbrR4) {
            firstOrder.p[cell] = true;
            break;
          }
        }
      }
      break;
    }
  }

  while(iteration < maxIterations) {
    res = F0;

    // discretization
    // WENO-JP-5
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];
      if(a_isHalo(cellId)) continue;
      a_levelSetRHS(cellId, set) = F0;
      if(!firstOrder.p[cell]) {
        for(MInt i = 0; i < nDim; i++) {
          nghbrL = a_bandNghbrIdsG(cellId, 2 * i, set);
          nghbrL2 = a_bandNghbrIdsG(nghbrL, 2 * i, set);
          nghbrL3 = a_bandNghbrIdsG(nghbrL2, 2 * i, set);
          nghbrR = a_bandNghbrIdsG(cellId, 2 * i + 1, set);
          nghbrR2 = a_bandNghbrIdsG(nghbrR, 2 * i + 1, set);
          nghbrR3 = a_bandNghbrIdsG(nghbrR2, 2 * i + 1, set);
          dPlusL3 = a_levelSetFunctionG(nghbrL2, set) - a_levelSetFunctionG(nghbrL3, set);
          dPlusL2 = a_levelSetFunctionG(nghbrL, set) - a_levelSetFunctionG(nghbrL2, set);
          dPlusL = a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(nghbrL, set);
          dPlus = a_levelSetFunctionG(nghbrR, set) - a_levelSetFunctionG(cellId, set);
          dPlusR = a_levelSetFunctionG(nghbrR2, set) - a_levelSetFunctionG(nghbrR, set);
          dPlusR2 = a_levelSetFunctionG(nghbrR3, set) - a_levelSetFunctionG(nghbrR2, set);
          a = (dPlusL2 - dPlusL3) * m_FgCellDistance;
          b = (dPlusL - dPlusL2) * m_FgCellDistance;
          c = (dPlus - dPlusL) * m_FgCellDistance;
          d = (dPlusR - dPlus) * m_FgCellDistance;
          IS0 = 13.0 * POW2(a - b) + 3.0 * POW2(a - 3 * b);
          IS1 = 13.0 * POW2(b - c) + 3.0 * POW2(b + c);
          IS2 = 13.0 * POW2(c - d) + 3.0 * POW2(3 * c - d);
          alpha0 = F1 / POW2(eps + IS0);
          alpha1 = F6 / POW2(eps + IS1);
          alpha2 = F3 / POW2(eps + IS2);
          omega0 = alpha0 / (alpha0 + alpha1 + alpha2);
          omega2 = alpha2 / (alpha0 + alpha1 + alpha2);
          PsiMinus = F1B3 * omega0 * (a - 2 * b + c) + F1B6 * (omega2 - F1B2) * (b - 2 * c + d);
          a = (dPlusR2 - dPlusR) * m_FgCellDistance;
          b = (dPlusR - dPlus) * m_FgCellDistance;
          d = (dPlusL - dPlusL2) * m_FgCellDistance;
          IS0 = 13.0 * POW2(a - b) + 3.0 * POW2(a - 3 * b);
          IS1 = 13.0 * POW2(b - c) + 3.0 * POW2(b + c);
          IS2 = 13.0 * POW2(c - d) + 3.0 * POW2(3 * c - d);
          alpha0 = F1 / POW2(eps + IS0);
          alpha1 = F6 / POW2(eps + IS1);
          alpha2 = F3 / POW2(eps + IS2);
          omega0 = alpha0 / (alpha0 + alpha1 + alpha2);
          omega2 = alpha2 / (alpha0 + alpha1 + alpha2);
          PsiPlus = F1B3 * omega0 * (a - 2 * b + c) + F1B6 * (omega2 - F1B2) * (b - 2 * c + d);
          a = F1B12 * m_FgCellDistance * (-dPlusL2 + F7 * dPlusL + F7 * dPlus - dPlusR) - PsiMinus;
          b = F1B12 * m_FgCellDistance * (-dPlusL2 + F7 * dPlusL + F7 * dPlus - dPlusR) + PsiPlus;

          a_levelSetRHS(cellId, set) +=
              F1B2
              * ((a_levelSetSign(cellId, set) + F1) * mMax(POW2(mMax(a, F0)), POW2(mMin(b, F0)))
                 - (a_levelSetSign(cellId, set) - F1) * mMax(POW2(mMin(a, F0)), POW2(mMax(b, F0))));
        }
      } else {
        for(MInt i = 0; i < nDim; i++) {
          a = (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(a_bandNghbrIdsG(cellId, 2 * i, set), set))
              * m_FgCellDistance;
          b = (a_levelSetFunctionG(a_bandNghbrIdsG(cellId, 2 * i + 1, set), set) - a_levelSetFunctionG(cellId, set))
              * m_FgCellDistance;
          a_levelSetRHS(cellId, set) +=
              F1B2
              * ((a_levelSetSign(cellId, set) + F1) * mMax(POW2(mMax(a, F0)), POW2(mMin(b, F0)))
                 - (a_levelSetSign(cellId, set) - F1) * mMax(POW2(mMin(a, F0)), POW2(mMax(b, F0))));
        }
      }
    }

    // final computation of a_levelSetRHS without forcing -> distinguish between limited (crMode = 4 ) and not limited
    // schemes:
    if(crMode == 4) {
      converged = true;
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        if(a_isHalo(cellId)) continue;
        if(fabs(F1 - sqrt(a_levelSetRHS(cellId, set))) > 1.0) {
          converged = false;
          a_levelSetRHS(cellId, set) = m_signG[IDX_LSSET(cellId, set)] * (F1 - sqrt(a_levelSetRHS(cellId, set)));
        } else {
          a_levelSetRHS(cellId, set) = F0;
        }
      }
    } else {
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        if(a_isHalo(cellId)) continue;
        a_levelSetRHS(cellId, set) = m_signG[IDX_LSSET(cellId, set)] * (F1 - sqrt(a_levelSetRHS(cellId, set)));
      }
    }


    // forcing term - only applied if crMode > 0 (CR scheme applied)
    switch(crMode) {
      case 0: { // zero disables forcing
        break;
      }
      case 1: { // CR1 - see Diss/Papers Daniel Hartmann
        MInt cnt, overallCnt;
        MBool forcing;
        MFloat sum;
        overallCnt = 0;
        for(MInt k = 0; k < noCRCells; k++) {
          forcing = true;
          sum = F0;
          cnt = 0;
          if(a_isHalo(crCells[k])) continue;
          for(MInt i = 0; i < m_noDirs; i++) {
            if(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)] != -1) {
              if(a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set)
                     * a_levelSetFunctionG(crCells[k], set)
                 < F0) {
                sum += factors[overallCnt++] * a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set);
                cnt++;
              } else
                forcing = false;
            }
          }
          if(forcing) {
            a_levelSetRHS(crCells[k], set) +=
                m_relaxationFactor * m_FgCellDistance * (sum / (MFloat)cnt - a_levelSetFunctionG(crCells[k], set));
          }
        }
        break;
      }
      case 6:
      case 2: { // CR2 - see Diss/Papers Daniel Hartmann
        MBool forcing;
        MFloat sum;
        //        MInt cnt =0;
        for(MInt k = 0; k < noCRCells; k++) {
          forcing = true;
          sum = F0;
          if(a_isHalo(crCells[k])) continue;
          //          cnt=0;
          // check whether sign is changed
          if(a_levelSetSign(crCells[k], set) > 0 && a_levelSetFunctionG(crCells[k], set) < 0) {
            forcing = false;
          }
          if(a_levelSetSign(crCells[k], set) < 0 && a_levelSetFunctionG(crCells[k], set) > 0) {
            forcing = false;
          }
          for(MInt i = 0; i < m_noDirs; i++) {
            if(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)] != -1) {
              if(a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set)
                     * a_levelSetFunctionG(crCells[k], set)
                 < F0) {
                sum += a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set);
                //                cnt++;
              } else {
                forcing = false;
                break; // Stephan : to not force when in one direction no sign is changing Stephan: debug
              }
            }
          }
          if(forcing) {
            a_levelSetRHS(crCells[k], set) +=
                m_relaxationFactor * m_FgCellDistance * (factors[k] * sum - a_levelSetFunctionG(crCells[k], set));
          }
        }
        break;
      }
      case 3: { // CR3 - ? Conservative Level Set only!
        for(MInt k = 0; k < noCRCells; k++) {
          if(a_isHalo(crCells[k])) continue;
          a_levelSetRHS(crCells[k], set) +=
              m_relaxationFactor * m_FgCellDistance
              * (atanh(m_hypTanLSF[IDX_LSSET(crCells[k], set)] * F2 - F1) * m_gCellDistance
                 - a_levelSetFunctionG(crCells[k], set));
        }
        break;
      }
      case 4: { // limited CR2
        MBool forcing;
        MFloat sum;
        for(MInt k = 0; k < noCRCells; k++) {
          if(a_isHalo(crCells[k])) continue;
          if(abs(a_levelSetRHS(crCells[k], set)) < 0.000000000001) continue;
          sum = F0;
          forcing = true;
          for(MInt i = 0; i < m_noDirs; i++) {
            if(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)] != -1) {
              if(a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set)
                     * a_levelSetFunctionG(crCells[k], set)
                 < F0)
                sum += a_levelSetFunctionG(m_phiRatioCells[crCells[k]][IDX_LSSET(i, set)], set);
              else
                forcing = false;
            }
          }
          if(forcing)
            a_levelSetRHS(crCells[k], set) +=
                m_relaxationFactor * m_FgCellDistance * (factors[k] * sum - a_levelSetFunctionG(crCells[k], set));
        }
        break;
      }
      default: {
        mTerm(1, AT_, "Unknown crMode");
      }
    }

    // solver
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = m_cellList[cell];
      if(a_isHalo(cellId)) continue;

      a_levelSetFunctionG(cellId, set) += a_levelSetRHS(cellId, set) * dt;
      res = mMax(res, ABS(a_levelSetRHS(cellId, set)));
      // debug: if this line is used instead of the previous one, stephans old results can be estabished!
      //       res = mMax( res, ABS(a_levelSetRHS(cellId , set) *dt) );
    }

    iteration++;

    MFloat* q;
    q = (MFloat*)&a_levelSetFunctionG(0, 0);
    exchangeLs(q, set, 1);

    // exchange residual and data
    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "res");

    // Stephan: minIteration included (result inly different for crMode == 6 where fully reinitialization of band is
    // guaranteed)
    if((res < m_reinitConvergence && iteration > minIteration) || converged) {
      if(set == m_startSet) {
        if(domainId() == 0) {
          FILE* datei;
          stringstream reinitFile;
          reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
          datei = fopen((reinitFile.str()).c_str(), "a+");
          fprintf(datei, "  %d", iteration);
          fclose(datei);
        }
      }
      return res;
    }
  }

  if(set == m_startSet) {
    if(domainId() == 0) {
      FILE* datei;
      stringstream reinitFile;
      reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
      datei = fopen((reinitFile.str()).c_str(), "a+");
      fprintf(datei, "  %d", iteration);
      fclose(datei);
    }
  }
  // as during regridding the convergence criterium is decreased by a factor of 3, here we reset it to the original
  // convergence criteria
  m_reinitConvergence = m_reinitConvergenceReset;
  return res;
}


//----------------------------------------------------------------------------

/**
 * \fn void LsCartesianSolver::initializeGControlPoint()
 * \brief this function is used to initialize the control point.
 * \author Sitthikrit Leckpool, Sep. 2010
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::initializeGControlPoint() {
  TRACE();

  MFloat elapsedTime = F0;
  if(m_restart) {
    elapsedTime = time();
  }

  /*! \page propertiesLS
   \section levelSetCtrlPntMethod
   <code>MFloat LsCartesianSolver::m_GCtrlPntMethod </code>\n
   default = 2 \n \n
   Choose control point methods.\n
   Possible values are:
   <ul>
   <li>1 for single body.</li>
   <li>2 for multiple bodies.</li>
   </ul>
   Keywords: <i>LEVELSET</i>
  */
  m_GCtrlPntMethod = 2; // control point 2 as default
  m_GCtrlPntMethod = Context::getSolverProperty<MInt>("levelSetCtrlPntMethod", m_solverId, AT_, &m_GCtrlPntMethod);
  // set initial orientation and position
  // if 2D, set the 3rd dimension to zero and w[] to [0,0,1], in-plane xy rotation
  // Control point alway refers 3D coordinates
  MFloat u[3], v[3], w[3];
  // Example, Anchor the reference orgin at bounding box mid point, or any as preference.
  MFloat InitPos[3] = {F0, F0, F0};
  //   InitPos[0] = m_geometry->m_mbMidPnt[0];
  //   InitPos[1] = m_geometry->m_mbMidPnt[1];
  //   InitPos[2] = m_geometry->m_mbMidPnt[2];
  //   IF_CONSTEXPR(nDim == 2) InitPos[2] = 0.0;
  // Example, +30 degree xy plane orientation
  // Orientation about "current reference origin",
  // if orientation about geometry mid point, do "InitOrientation" before "InitPosition"
  // unit normalized vector is needed
  u[0] = 1.0;
  v[0] = 0.0;
  w[0] = 0.0;
  u[1] = 0.0;
  v[1] = 1.0;
  w[1] = 0.0;
  u[2] = 0.0;
  v[2] = 0.0;
  w[2] = 1.0;
  IF_CONSTEXPR(nDim == 2) {
    u[2] = v[2] = w[0] = w[1] = 0.0;
    w[2] = 1.0;
  }
  // initialize control points
  //
  //
  switch(m_GCtrlPntMethod) {
    case 1: {
      // not adjusted to multiple bodies...
      if(m_noBodyBndryCndIds > 1)
        mTerm(1, AT_,
              "You are using GCtrlPntMethod 1, which is not adjusted to multiple bodies, but you are trying to "
              "compute multiple bodies. Please check!");
      m_gCtrlPnt.CtrlPnt1_Initialize(m_geometry, nDim);
      m_gCtrlPnt.CtrlPnt1_InitOrientation(u, v, w);
      computeBodyPropertiesForced(1, InitPos, 0, elapsedTime);
      IF_CONSTEXPR(nDim == 2) InitPos[2] = 0.0;
      m_gCtrlPnt.CtrlPnt1_InitPosition(InitPos);
      break;
    }
    case 2: {
      m_gCtrlPnt.CtrlPnt2_Initialize(m_geometry, nDim);
      for(MInt body = 0; body < m_noEmbeddedBodies; body++) {
        const MInt bcId = m_bodyBndryCndIds[body];
        computeBodyPropertiesForced(1, InitPos, body, elapsedTime);
        IF_CONSTEXPR(nDim == 2) InitPos[2] = 0.0;
        m_gCtrlPnt.CtrlPnt2_InitOrientation(u, v, w, bcId);
        m_gCtrlPnt.CtrlPnt2_InitPosition(InitPos, bcId);
#if defined LS_DEBUG
        stringstream controlfile;
        controlfile << "Ctrl_Body_" << body << ".stl";
        m_gCtrlPnt.CtrlPnt2_CtrlPntToSTL((controlfile.str()).c_str(), bcId);
#endif
      }
      m_gCtrlPnt.CtrlPnt2_Update();
      break;
    }
    case 3:
      // TODO labels:LS,totest Validate
      // Attempt to create body-free version ~jv
      {
        m_gCtrlPnt.CtrlPnt1_Initialize(m_geometry, nDim);
        m_gCtrlPnt.CtrlPnt1_InitOrientation(u, v, w);
        IF_CONSTEXPR(nDim == 2) InitPos[2] = 0.0;
        m_gCtrlPnt.CtrlPnt1_InitPosition(InitPos);
        break;
      }
    default: {
      mTerm(1, AT_, "Unknown GCtrlPntMethod");
    }
  }
}

//----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::initializeIntegrationScheme() {
  TRACE();

  /*! \page propertiesLS
    \section nogRKSteps
    <code>MInt LsCartesianSolver::initializeIntegrationScheme::m_nogRKSteps </code>\n
    default = <code>0</code>\n \n
    Sets the number of generalized Runge Kutta steps. \n
    <ul>
    <li> non-negative integers </li>
    </ul>
    Keywords: <i>LEVELSET, INTEGRATION, RUNGE KUTTA </i>
  */
  m_nogRKSteps = Context::getSolverProperty<MInt>("nogRKSteps", m_solverId, AT_);

  // alpha
  m_gRKalpha = new MFloat[m_nogRKSteps];
  for(MInt i = 0; i < m_nogRKSteps; i++) {
    /*! \page propertiesLS
     \section grkalpha-step
     <code>MFloat LsCartesianSolver::m_gRKalpha </code>\n
     Sets the alpha coefficient for the Runge Kutta steps used in the integration scheme.
     Should contain m_nogRKSteps values.
     Keywords: <i>LEVELSET, RUNGEKUTTA</i>
    */
    m_gRKalpha[i] = Context::getSolverProperty<MFloat>("grkalpha-step", m_solverId, AT_, i);
  }

  // reset step
  m_gRKStep = 0;
  if(!m_restart && m_LSSolver) {
    m_time = 0.0;
    globalTimeStep = 0;
  }
}

//----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::initializeIntegrationScheme_semiLagrange() {
  TRACE();

  // reset step - not really required for semiLagrange, but do nevertheless for security reasons!
  m_gRKStep = 0;

  // for a restart the value will be over-written later
  for(MInt i = 0; i < nDim * m_noEmbeddedBodies; i++) {
    m_semiLagrange_xShift_ref[i] = F0;
  }

  initRotatingLS();
}


// ----------------------------------------------------------------------------------------


// sets the m_isBndryCell member of levelSetFrontCells
template <MInt nDim>
void LsCartesianSolver<nDim>::setGCellBndryProperty() {
  TRACE();

  // checks whether a neighbor exists in all directions
  // if not, the cell is a boundary cell
  // this is only valid for isotropic G grids!

  // reset m_isBndryCell property:

  for(MInt set = 0; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      MInt cellId = a_G0CellId(id, set);
      a_isBndryCellG(cellId) = false;
    }
  }

  // set m_isBndryCell property - cell is bndry cell if at least 1 level set function identifies it as such
  for(MInt set = 0; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      MInt cellId = a_G0CellId(id, set);
      for(MInt dir = 0; dir < 2 * nDim; dir++) {
        if(a_hasNeighbor(cellId, dir) == 0) {
          a_isBndryCellG(cellId) = true;
          break;
        }
      }
    }
  }
}


// ----------------------------------------------------------------------------------------


// sets the m_isBndryCell member of levelSetFrontCells
template <MInt nDim>
void LsCartesianSolver<nDim>::applyLevelSetBoundaryConditions() {
  TRACE();

  for(MInt set = m_startSet; set < m_noSets; set++) {
    switch(m_levelSetBoundaryCondition) {
      case 17516: {
        MFloat x, err = 0.05;
        MFloat deltaX = m_flameRadiusOffset;
        // two flame code
        if(m_twoFlames) {
          for(MInt id = 0; id < a_noGBndryCells(set); id++) {
            a_isGBoundaryCellG(a_gBndryCellId(id, set), set) = true;
            x = c_coordinate(id, 0);

            if((x > (-m_radiusFlameTube2 * 1.3 - err + m_xOffsetFlameTube))
               && (x < (m_radiusFlameTube2 * 1.3 + err + m_xOffsetFlameTube))) {
              a_levelSetFunctionG(a_gBndryCellId(id, set), set) =
                  ABS(c_coordinate(a_gBndryCellId(id, set), 0) - m_xOffsetFlameTube) - m_radiusFlameTube * 1.3;

            } else if((x > (-m_radiusFlameTube2 * 1.3 - err + m_xOffsetFlameTube2))
                      && (x < (m_radiusFlameTube2 * 1.3 + err + m_xOffsetFlameTube2))) {
              a_levelSetFunctionG(a_gBndryCellId(id, set), set) =
                  ABS(c_coordinate(a_gBndryCellId(id, set), 0) - m_xOffsetFlameTube2) - m_radiusFlameTube2 * 1.3;
            }
          }
          // one flame code
        } else {
          for(MInt id = 0; id < a_noGBndryCells(set); id++) {
            MInt cellId = a_gBndryCellId(id, set);
            a_isGBoundaryCellG(cellId, set) = true;
            a_levelSetFunctionG(cellId, set) =
                ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube) - (m_radiusFlameTube + deltaX);
            a_curvatureG(cellId, 0) = 0;
            for(MInt i = 0; i < nDim; i++)
              a_extensionVelocityG(cellId, i, set) = F0;
          }
        }
        break;
      }
      case 1751600: {
        //        MFloat x,z,err= 0.05;
        MFloat deltaX = m_flameRadiusOffset;
        // two flame code
        for(MInt id = 0; id < a_noGBndryCells(set); id++) {
          MInt cellId = a_gBndryCellId(id, set);
          a_isGBoundaryCellG(cellId, set) = true;
          //          x = c_coordinate( gCellId ,  0 );
          //          z = c_coordinate( gCellId ,  2 );

          MFloat minXG = ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube) - (m_jetHalfWidth + deltaX);
          MFloat minZG = ABS(c_coordinate(cellId, 2)) - (m_jetHalfLength + deltaX);
          MFloat maxG = mMax(minXG, minZG);
          a_levelSetFunctionG(cellId, set) =
              m_initialFlameHeight * sqrt(2.0 - 1.0) * maxG + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

          a_curvatureG(cellId, 0) = 0;

          for(MInt i = 0; i < nDim; i++)
            a_extensionVelocityG(cellId, i, set) = F0;
        }
        break;
      }

      case 5401000: {
        MFloat deltaX = 0.05;
        // two flame code
        for(MInt id = 0; id < a_noGBndryCells(set); id++) {
          MInt cellId = a_gBndryCellId(id, set);
          a_isGBoundaryCellG(cellId, set) = true;

          MFloat radius = sqrt(POW2(c_coordinate(cellId, 0)) + POW2(c_coordinate(cellId, 2)));

          MFloat maxG = ABS(radius) - (0.515 + deltaX);
          a_levelSetFunctionG(cellId, set) =
              m_initialFlameHeight * sqrt(2.0 - 1.0) * maxG + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

          a_curvatureG(cellId, 0) = 0;

          for(MInt i = 0; i < nDim; i++)
            a_extensionVelocityG(cellId, i, set) = F0;
        }
        break;
      }
      default: {
        break; // some testcases do not need boundary conditions
      }
    }
  }
}


// ----------------------------------------------------------------------------------------


/**
 * \fn void LsCartesianSolver::initializeGField()
 * \brief Initializes the solver values with the values of the undisturbed flow
  The values are given by the property file. The conservative and primitive
  variables are calculated and set with the given values as also the
  variables of the cells
*/
template <MInt nDim>
void LsCartesianSolver<nDim>::initializeGField() {
  TRACE();

  ASSERT(!m_GFieldInitFromSTL, "");

  MInt noCells = a_noCells();
  //---

  IF_CONSTEXPR(nDim == 2) {
    switch(m_initialCondition) {
      case 50431:
      case 504312: {
        MFloat radius;
        MFloat x, y;
        MFloat h = 2.5;
        MFloat w = 0.25;
        MFloat R0 = 1.5;
        MFloat dy = 2.5;
        if(m_initialCondition == 504312) {
          dy = -1.9;
        }

        for(MInt set = m_startSet; set < m_noSets; set++) {
          for(MInt cellId = 0; cellId < noCells; cellId++) {
            radius = POW2(c_coordinate(cellId, 0)) + POW2(c_coordinate(cellId, 1) - dy);
            radius = sqrt(radius);
            x = c_coordinate(cellId, 0);
            y = c_coordinate(cellId, 1) - dy;

            if(radius < R0) {
              if(y > h - R0) {
                if(ABS(x) > w) {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, sqrt(POW2(ABS(x) - w) + POW2(y - h + R0)));
                } else {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, ABS(y - h + R0));
                }
              } else {
                if(ABS(x) < w) {
                  a_levelSetFunctionG(cellId, set) = mMax(-ABS(ABS(x) - w), -ABS(y - h + R0));
                } else {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, ABS(x) - w);
                }
              }
            } else {
              if(ABS(x) > w || y > 0) {
                a_levelSetFunctionG(cellId, set) = R0 - radius;
              } else {
                a_levelSetFunctionG(cellId, set) = -sqrt(POW2(ABS(x) - w) + POW2(y + R0));
              }
            }
          }
        }

        break;
      }
      case 203: {
        // rotating disk problem
        MFloat radius;
        MFloat x, y;
        MFloat h = 2.5;
        MFloat w = 0.25;
        MFloat R0 = 1.5;
        MFloat dy = 2.5;
        for(MInt set = m_startSet; set < m_noSets; set++) {
          for(MInt cellId = 0; cellId < noCells; cellId++) {
            radius =
                POW2(c_coordinate(cellId, 0) - a_meanCoord(0)) + POW2(c_coordinate(cellId, 1) - a_meanCoord(1) - dy);
            radius = sqrt(radius);
            x = c_coordinate(cellId, 0);
            y = c_coordinate(cellId, 1) - dy;

            if(radius < R0) {
              if(y > h - R0) {
                if(ABS(x) > w) {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, sqrt(POW2(ABS(x) - w) + POW2(y - h + R0)));
                } else {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, ABS(y - h + R0));
                }
              } else {
                if(ABS(x) < w) {
                  a_levelSetFunctionG(cellId, set) = mMax(-ABS(ABS(x) - w), -ABS(y - h + R0));
                } else {
                  a_levelSetFunctionG(cellId, set) = mMin(R0 - radius, ABS(x) - w);
                }
              }
            } else {
              if(ABS(x) > w || y > 0) {
                a_levelSetFunctionG(cellId, set) = R0 - radius;
              } else {
                a_levelSetFunctionG(cellId, set) = -sqrt(POW2(ABS(x) - w) + POW2(y + R0));
              }
            }
          }
        }
        break;
      }

      case 450: {
        // cylinder moving boundary test

        /*! \page propertiesLS
          \section radius
          <code>MFloat R0</code>\n
          default:  0.5 \n \n
          Radius of moving cylinder in the the level-set test-case 450. \n
          Keywords: <i>LEVEL_SET, BODY</i>
        */

        MFloat R0 = F1B2;
        R0 = Context::getSolverProperty<MFloat>("radius", m_solverId, AT_, &R0);

        MFloat center[nDim];
        MInt set = 0;
        for(MInt i = 0; i < nDim; i++) {
          center[i] = F0;
        }
        if(Context::propertyExists("initialBodyCenter", m_solverId)) {
          if(Context::propertyLength("initialBodyCenter", m_solverId) != nDim) {
            mTerm(1, AT_, "Property 'initialBodyCenter' has invalid amount of entries!");
          }
          for(MInt i = 0; i < nDim; i++) {
            center[i] = Context::getSolverProperty<MFloat>("initialBodyCenter", m_solverId, AT_, i);
          }
        }
        MFloat radius;

        for(MInt cellId = 0; cellId < noCells; cellId++) {
          radius = POW2(c_coordinate(cellId, 0) - center[0]) + POW2(c_coordinate(cellId, 1) - center[1]);
          radius = sqrt(radius);
          a_levelSetFunctionG(cellId, set) = radius - R0;
        }
        break;
      }

      case 17516: {
        MFloat x, err = 0.1;
        MFloat deltaX = m_flameRadiusOffset;
        MInt set = 0;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          x = c_coordinate(cellId, 0);

          if(m_twoFlames) {
            err = 0.3;
            if((x > (-m_radiusFlameTube2 * 1.3 - err + m_xOffsetFlameTube))
               && (x < (m_radiusFlameTube2 * 1.3 + err + m_xOffsetFlameTube))) {
              a_levelSetFunctionG(cellId, set) =
                  sqrt(4.0 - 1.0) * (ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube) - m_radiusFlameTube * 1.3)
                  + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

            } else if((x > (-m_radiusFlameTube2 * 1.3 - err + m_xOffsetFlameTube2))
                      && (x < (m_radiusFlameTube2 * 1.3 + err + m_xOffsetFlameTube2))) {
              a_levelSetFunctionG(cellId, set) =
                  sqrt(4.0 - 1.0) * (ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube2) - m_radiusFlameTube2 * 1.3)
                  + c_coordinate(cellId, 1) - m_yOffsetFlameTube2;

              // level set function outside of the domain
            } else {
              a_levelSetFunctionG(cellId, set) = c_coordinate(cellId, 1) + 0.5;
            }

            // one flame initial level set function
          } else {
            if((x > -m_radiusFlameTube - deltaX - 0.5) && (x < m_radiusFlameTube + deltaX + 0.5)) {
              a_levelSetFunctionG(cellId, set) =
                  sqrt(4.0 - 1.0) * (ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube) - (m_radiusFlameTube + deltaX))
                  + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

              // level set function outside of the domain
            } else {
              a_levelSetFunctionG(cellId, set) = c_coordinate(cellId, 1) + 0.5;
            }
          }
        }

        for(MInt cellId = 0; cellId < noCells; cellId++) {
          for(MInt s = 1; s < m_noSets; s++)
            a_levelSetFunctionG(cellId, s) = a_levelSetFunctionG(cellId, 0);
        }
        break;
      }
      default: {
        break; // some testcases do not need initialization
      }
    }
  }
  else IF_CONSTEXPR(nDim == 3) {
    switch(m_initialCondition) {
      case 0: {
        // parallel inflow field with G sphere
        MFloat radius;
        MInt set = 0;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          radius = POW2(c_coordinate(cellId, 0)) + POW2(c_coordinate(cellId, 1)) + POW2(c_coordinate(cellId, 2));
          radius = sqrt(radius);

          a_levelSetFunctionG(cellId, set) = (3.0 - radius);
        }
        break;
      }

      case 450: {
        // cylinder moving boundary test
        MFloat R0 = F1B2;
        MInt set = 0;

        /*! \page propertiesLS
          \section radius
          <code>MFloat R0</code>\n
          default:  0.5 \n \n
          Radius of moving cylinder in the the level-set test-case 450. \n
          Keywords: <i>LEVEL_SET, BODY</i>
        */
        R0 = Context::getSolverProperty<MFloat>("radius", m_solverId, AT_, &R0);

        MFloat center[nDim];
        for(MInt i = 0; i < nDim; i++) {
          center[i] = F0;
        }
        if(Context::propertyExists("initialBodyCenter", m_solverId)) {
          if(Context::propertyLength("initialBodyCenter", m_solverId) != nDim) {
            mTerm(1, AT_, "Property 'initialBodyCenter' has invalid amount of entries!");
          }
          for(MInt i = 0; i < nDim; i++) {
            center[i] = Context::getSolverProperty<MFloat>("initialBodyCenter", m_solverId, AT_, i);
          }
        }
        MFloat radius;

        for(MInt cellId = 0; cellId < noCells; cellId++) {
          radius = F0;
          for(MInt i = 0; i < nDim; i++) {
            radius += POW2(c_coordinate(cellId, i) - center[i]);
          }
          radius = sqrt(radius);
          a_levelSetFunctionG(cellId, set) = radius - R0;
        }
        break;
      }

      case 1751600: {
        MFloat x, z;
        MFloat deltaX = m_flameRadiusOffset;
        MInt set = 0;
        MInt count = 0;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          x = c_coordinate(cellId, 0);
          z = c_coordinate(cellId, 2);

          if((x > -m_radiusFlameTube - deltaX - 0.5) && (x < m_radiusFlameTube + deltaX + 0.5)
             && (z > -m_jetHalfLength - 0.5) && (z < m_jetHalfLength + 0.5)) {
            MFloat minXG = ABS(c_coordinate(cellId, 0) - m_xOffsetFlameTube) - (m_jetHalfWidth + deltaX);
            MFloat minZG = ABS(c_coordinate(cellId, 2)) - (m_jetHalfLength + deltaX);
            MFloat maxG = mMax(minXG, minZG);
            a_levelSetFunctionG(cellId, set) =
                m_initialFlameHeight * sqrt(2.0 - 1.0) * maxG + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

            // level set function outside of the domain
          } else {
            a_levelSetFunctionG(cellId, set) = c_coordinate(cellId, 1) - 1.0;
          }
          if(approx(a_levelSetFunctionG(cellId, set), 0.0, MFloatEps)) count++;
        }

        for(MInt cellId = 0; cellId < noCells; cellId++) {
          for(MInt s = 1; s < m_noSets; s++)
            a_levelSetFunctionG(cellId, s) = a_levelSetFunctionG(cellId, 0);
        }
        break;
      }

      case 5401000: {
        MFloat x, z;
        MFloat deltaX = 0.05;
        MInt set = 0;
        // MInt count=0;
        MFloat radius;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          x = c_coordinate(cellId, 0);
          z = c_coordinate(cellId, 2);
          radius = sqrt(POW2(x) + POW2(z));

          if((radius < 0.5 + 0.3)) {
            MFloat maxG = ABS(radius) - (0.53 + deltaX);
            a_levelSetFunctionG(cellId, set) =
                m_initialFlameHeight * sqrt(2.0 - 1.0) * maxG + c_coordinate(cellId, 1) - m_yOffsetFlameTube;

            // level set function outside of the domain
          } else {
            a_levelSetFunctionG(cellId, set) = c_coordinate(cellId, 1) + 0.0;
          }
          // if(approx(a_levelSetFunctionG( cellId , set) , 0.0, MFloatEps))
          // count++;
        }

        for(MInt cellId = 0; cellId < noCells; cellId++) {
          for(MInt s = 1; s < m_noSets; s++)
            a_levelSetFunctionG(cellId, s) = a_levelSetFunctionG(cellId, 0);
        }
        break;
      }

      default: {
        break; // some testcases do not need initialization
      }
    }
  }

  // initialize all old cell variables
  if(m_semiLagrange) {
    for(MInt set = 0; set < m_noSets; set++) {
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        a_oldLevelSetFunctionG(cellId, set) = a_levelSetFunctionG(cellId, set);
      }
    }
  }
  m_log << "Level set function initialized on " << noCells << " cells" << endl;
}

//----------------------------------------------------------------------------


/**
 * \fn void LsCartesianSolver::constructGFieldFromSTL(MInt ConstructFlag)
 * \brief Used for initializing G field into the domain.
             ConstructFlag
             == 0: only reinitialization part will be called
                   (used cell-type depends on the initMode)
             == 1: only initialization part will be called
                   (which uses all cells in the domain)
             == 2: both initialization and reinitialization part will be called (currently unused!)
                   (however only cells in the band are used for the initialization,
                   only usefull after band-cells have been found!)
             == 3: only initialization part will be called
                   (which uses all band-cells in the domain)
                   only usefull after band-cells have been found!
 *
 * \author Sitthikrit Leckpool, Tim Wegmann
 * \date Sep. 2010, May 2018
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::constructGFieldFromSTL(MInt ConstructFlag) {
  TRACE();

  // check dimension
  IF_CONSTEXPR(nDim == 1) mTerm(1, AT_, "Initialize G Field from STL does not support 1D problems");
  //  MInt maxLevel = maxLevel();

  MInt initMode = 1; // (Recommended version!)
  /*! \page propertiesLS
     \section GFieldFromSTLInitMode
     default = 1 \n \n
     Choose method how to generate Levelset from STL.\n
     In the Initialization-part, possible values are:
     <ul>
     <li>0: only G0-cells are initialized by the stl </li>
     <li>1: G-cells specified by the ContstructFlag are initialized by the stl (for ContstructFlag=1: all G-cells, for
     ContstructFlag=2/3: all Band-cells,) </li> <li>2: G-cells specified by the ContstructFlag are initialized by a
     sphere (with radius 0.5 in the origin). Very special testcase, barely useful! </li> <li>3: only G0-cells are
     initialized by the stl </li>
     </ul>
  */
  initMode = Context::getSolverProperty<MInt>("GFieldFromSTLInitMode", m_solverId, AT_, &initMode);

  // define lambda for setting a_hasPositiveSign with MInts
  auto signToBool = [](const MInt& plusMinus) { return (plusMinus > 0) ? true : false; };

  // a) Initialization:
  //---------------------------------------------------------------------------------------------------
  if(ConstructFlag) {
    switch(ConstructFlag) {
      case 1:
        // Use all G-cells:
        {
          // 1) Reset a_isGZeroCell
          for(MInt set = m_startSet; set < m_noSets; set++) {
            for(MInt GcellId = 0; GcellId < a_noCells(); GcellId++) {
              a_isGZeroCell(GcellId, set) = false;
            }
          }

          // Initialization of the collected-levelset!
          for(MInt set = 0; set < m_startSet; set++) {
            for(MInt GcellId = 0; GcellId < a_noCells(); GcellId++) {
              a_levelSetFunctionG(GcellId, set) = -std::numeric_limits<MFloat>::infinity();
            }
          }

          // 2) Determine the Levelset-Sign
          for(MInt set = m_startSet; set < m_noSets; set++) {
            for(MInt GcellId = 0; GcellId < a_noCells(); GcellId++) {
              MFloat target[3] = {0, 0, 0};
              if(a_isHalo(GcellId)) continue;
              if(!c_isLeafCell(GcellId)) continue;
              // if(a_level(GcellId) < maxLevel ) continue;

              for(MInt Dir = 0; Dir < nDim; Dir++)
                target[Dir] = c_coordinate(GcellId, Dir);

              if(m_GCtrlPntMethod == 1 || m_GCtrlPntMethod == 3) m_gCtrlPnt.CtrlPnt1_quvw(target);
              a_hasPositiveSign(GcellId, set) = signToBool(determineLevelSetSignFromSTL(target, set));
              if(m_levelSetSign[set] < 0) {
                a_hasPositiveSign(GcellId, set) = signToBool(-a_levelSetSign(GcellId, set));
              }
              a_levelSetFunctionG(GcellId, set) = a_levelSetSign(GcellId, set);
            }
          }

          // 3) Update G0-Cells
          determineG0Cells();

          // 4) Determine the distance to the Levelset-Interface (levelSet-value)
          for(MInt set = m_startSet; set < m_noSets; set++) {
            MInt noCells = a_noG0Cells(set);
            if(initMode > 0 && initMode < 3) noCells = a_noCells();
            for(MInt count = 0; count < noCells; count++) {
              MInt cellId = -1;
              if(initMode == 0 || initMode == 3) {
                cellId = a_G0CellId(count, set);
              } else {
                cellId = count;
              }
              if(a_isHalo(cellId)) continue;
              if(!c_isLeafCell(cellId)) continue;
              //        if(a_level(cellId) < maxLevel ) continue;

              if(initMode == 2) {
                MFloat radius{};
                for(MInt n = 0; n < nDim; n++) {
                  radius += c_coordinate(cellId, n) * c_coordinate(cellId, n);
                }
                radius = sqrt(radius);

                MFloat lvs = radius - 0.5;
                a_levelSetFunctionG(cellId, set) = lvs;

              } else {
                MInt closestElement;
                IF_CONSTEXPR(nDim == 2) {
                  MFloat closestPoint[2] = {F0, F0};
                  MFloat refPoint[2] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1)};
                  a_levelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                }
                else IF_CONSTEXPR(nDim == 3) {
                  MFloat closestPoint[3] = {F0, F0, F0};
                  MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
                  a_levelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                }
              }
            }
          }
          break;
        }
      case 2:
      case 3:
      case 5:
        // Use only band-cells:
        {
          if(m_GCtrlPntMethod == 2 && ConstructFlag != 5) {
            // control point method 2 needs to update normal vector, bounding box and adt
            m_gCtrlPnt.CtrlPnt2_UpdateAllNormalVector();
            m_gCtrlPnt.CtrlPnt2_Update();
          }

          // Initialization of the collected-levelset!
          for(MInt set = 0; set < m_startSet; set++) {
            if(ConstructFlag == 5 && !m_geometryChange[set]) continue;
            for(MInt id = 0; id < a_noBandCells(set); id++) {
              const MInt cellId = a_bandCellId(id, set);
              a_levelSetFunctionG(cellId, set) = -std::numeric_limits<MFloat>::infinity();
            }
          }

          // 2) Determine the Levelset-Sign
          for(MInt set = m_startSet; set < m_noSets; set++) {
            if(ConstructFlag == 5 && !m_geometryChange[set]) continue;
            for(MInt id = 0; id < a_noBandCells(set); id++) {
              const MInt cellId = a_bandCellId(id, set);
              MFloat target[3] = {0, 0, 0};
              if(a_isHalo(cellId)) continue;
              if(!c_isLeafCell(cellId)) continue;
              // if(a_level(GcellId) < maxLevel ) continue;

              for(MInt Dir = 0; Dir < nDim; Dir++)
                target[Dir] = c_coordinate(cellId, Dir);

              if(m_GCtrlPntMethod == 1 || m_GCtrlPntMethod == 3) m_gCtrlPnt.CtrlPnt1_quvw(target);
              a_hasPositiveSign(cellId, set) = signToBool(determineLevelSetSignFromSTL(target, set));
              if(m_levelSetSign[set] < 0) {
                a_hasPositiveSign(cellId, set) = signToBool(-a_levelSetSign(cellId, set));
              }
              a_levelSetFunctionG(cellId, set) = a_levelSetSign(cellId, set);
            }
          }

          // 3) Update G0-Cells and Band-Cells
          for(MInt set = 0; set < m_noSets; set++) {
            std::vector<MInt>().swap(m_bandCells[set]);
          }
          determineG0Cells();
          determineBandCells();

          // 4) Determine the distance to the Levelset-Interface (levelSet-value)
          for(MInt set = m_startSet; set < m_noSets; set++) {
            if(ConstructFlag == 5 && !m_geometryChange[set]) continue;
            MInt noCells = a_noG0Cells(set);
            if(initMode > 0 && initMode < 3) noCells = a_noBandCells(set);

            for(MInt count = 0; count < noCells; count++) {
              MInt cellId = -1;
              if(initMode == 0 || initMode == 3) {
                cellId = a_G0CellId(count, set);
              } else {
                cellId = a_bandCellId(count, set);
              }
              if(a_isHalo(cellId)) continue;
              if(!c_isLeafCell(cellId)) continue;
              // if(a_level(GcellId) < maxLevel ) continue;

              if(initMode == 2) {
                MFloat radius{};
                for(MInt n = 0; n < nDim; n++) {
                  radius += c_coordinate(cellId, n) * c_coordinate(cellId, n);
                }
                radius = sqrt(radius);

                const MFloat lvs = radius - 0.5;
                a_levelSetFunctionG(cellId, set) = lvs;

              } else {
                MInt closestElement;
                IF_CONSTEXPR(nDim == 2) {
                  MFloat closestPoint[2] = {F0, F0};
                  MFloat refPoint[2] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1)};
                  // Reset newly added band cells
                  if(abs(a_levelSetFunctionG(cellId, set)) > 1.0) {
                    a_hasPositiveSign(cellId, set) = signToBool(determineLevelSetSignFromSTL(refPoint, set));
                    if(m_levelSetSign[set] < 0) {
                      a_hasPositiveSign(cellId, set) = signToBool(-a_levelSetSign(cellId, set));
                    }
                    a_levelSetFunctionG(cellId, set) = a_levelSetSign(cellId, set);
                  }
                  a_levelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                }
                else IF_CONSTEXPR(nDim == 3) {
                  MFloat closestPoint[3] = {F0, F0, F0};
                  MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
                  // Reset newly added band cells
                  if(abs(a_levelSetFunctionG(cellId, set)) > 1.0) {
                    a_hasPositiveSign(cellId, set) = signToBool(determineLevelSetSignFromSTL(refPoint, set));
                    if(m_levelSetSign[set] < 0) {
                      a_hasPositiveSign(cellId, set) = signToBool(-a_levelSetSign(cellId, set));
                    }
                    a_levelSetFunctionG(cellId, set) = a_levelSetSign(cellId, set);
                  }
                  a_levelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                }
              }
            }
          }
          break;
        }
      case 6:
        // Use all cells; Update oldLevelSet:
        {
          // Initialization of the collected-levelset!
          for(MInt set = 0; set < m_startSet; set++) {
            if(!m_geometryChange[set]) continue;
            for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
              a_oldLevelSetFunctionG(cellId, set) = -std::numeric_limits<MFloat>::infinity();
            }
          }

          // 2) Determine the Levelset-Sign
          for(MInt set = m_startSet; set < m_noSets; set++) {
            if(!m_geometryChange[set]) continue;
            for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
              MFloat target[3] = {0, 0, 0};
              if(a_isHalo(cellId)) continue;
              if(!c_isLeafCell(cellId)) continue;

              for(MInt Dir = 0; Dir < nDim; Dir++)
                target[Dir] = c_coordinate(cellId, Dir);

              if(m_GCtrlPntMethod == 1 || m_GCtrlPntMethod == 3) m_gCtrlPnt.CtrlPnt1_quvw(target);
              a_hasPositiveSign(cellId, set) = signToBool(determineLevelSetSignFromSTL(target, set));
              if(m_levelSetSign[set] < 0) {
                a_hasPositiveSign(cellId, set) = signToBool(-a_levelSetSign(cellId, set));
              }
              a_oldLevelSetFunctionG(cellId, set) = a_levelSetSign(cellId, set);
            }
          }

          // 4) Determine the distance to the Levelset-Interface (levelSet-value)
          for(MInt set = m_startSet; set < m_noSets; set++) {
            if(!m_geometryChange[set]) continue;
            for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
              if(a_isHalo(cellId)) continue;
              if(!c_isLeafCell(cellId)) continue;
              // if(a_level(GcellId) < maxLevel ) continue;

              if(initMode == 2) {
                MFloat radius{};
                for(MInt n = 0; n < nDim; n++) {
                  radius += c_coordinate(cellId, n) * c_coordinate(cellId, n);
                }
                radius = sqrt(radius);

                const MFloat lvs = radius - 0.5;
                a_oldLevelSetFunctionG(cellId, set) = lvs;

              } else {
                MInt closestElement;
                if(nDim == 2) {
                  MFloat closestPoint[2] = {F0, F0};
                  MFloat refPoint[2] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1)};
                  a_oldLevelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                } else if(nDim == 3) {
                  MFloat closestPoint[3] = {F0, F0, F0};
                  MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
                  a_oldLevelSetFunctionG(cellId, set) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
                }
              }
            }
          }
          break;
        }
      case 4:
        // on the run construction of newly refined cells for stationary sets
        {
          ASSERT(initMode == 1, "");
          ASSERT(m_STLReinitMode == 2, "");

          // Initialization of the levelset!
          for(auto it = m_refinedCells.begin(); it != m_refinedCells.end(); it++) {
            const MInt cellId = it->first;
            const MInt set = it->second;
            ASSERT(cellId > -1 && cellId < a_noCells(), "");
            if(set > 0) {
              ASSERT(!m_computeSet_backup[set] || m_maxLevelChange, "");
              ASSERT(set >= 0 && set < m_noSets, "");
              a_oldLevelSetFunctionG(cellId, set) = -std::numeric_limits<MFloat>::infinity();
            } else {
              for(MInt setI = m_startSet; setI < m_noSets; setI++) {
                if(!m_maxLevelChange && m_computeSet_backup[setI] && !m_virtualSurgery) continue;
                a_oldLevelSetFunctionG(cellId, setI) = -std::numeric_limits<MFloat>::infinity();
              }
            }
          }

          // 2) Determine the Levelset-Sign
          for(auto it = m_refinedCells.begin(); it != m_refinedCells.end(); it++) {
            const MInt cellId = it->first;
            const MInt set = it->second;
            if(a_isHalo(cellId)) continue;
            if(!c_isLeafCell(cellId)) continue;
            // if(a_level(GcellId) < maxLevel ) continue;
            MFloat target[3] = {0, 0, 0};
            for(MInt dir = 0; dir < nDim; dir++) {
              target[dir] = c_coordinate(cellId, dir);
            }
            ASSERT(m_GCtrlPntMethod == 2, "");
            if(set > 0) {
              a_hasPositiveSign(cellId, set) = signToBool(determineLevelSetSignFromSTL(target, set));
              if(m_levelSetSign[set] < 0) {
                a_hasPositiveSign(cellId, set) = signToBool(-a_levelSetSign(cellId, set));
              }
              a_oldLevelSetFunctionG(cellId, set) = a_levelSetSign(cellId, set);
            } else {
              for(MInt setI = m_startSet; setI < m_noSets; setI++) {
                if(!m_maxLevelChange && m_computeSet_backup[setI] && !m_virtualSurgery) continue;
                a_hasPositiveSign(cellId, setI) = signToBool(determineLevelSetSignFromSTL(target, setI));
                if(m_levelSetSign[setI] < 0) {
                  a_hasPositiveSign(cellId, setI) = signToBool(-a_levelSetSign(cellId, setI));
                }
                a_oldLevelSetFunctionG(cellId, setI) = a_levelSetSign(cellId, setI);
              }
            }
          }


          // 4) Determine the distance to the Levelset-Interface (levelSet-value)
          for(auto it = m_refinedCells.begin(); it != m_refinedCells.end(); it++) {
            const MInt cellId = it->first;
            const MInt set = it->second;
            if(a_isHalo(cellId)) continue;
            if(!c_isLeafCell(cellId)) continue;
            // if(a_level(GcellId) < maxLevel ) continue;
            if(set > 0) {
              MInt closestElement;
              IF_CONSTEXPR(nDim == 2) {
                MFloat closestPoint[2] = {F0, F0};
                MFloat refPoint[2] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1)};
                a_oldLevelSetFunctionG(cellId, set) *=
                    computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set, m_sphereRadiusLimit);
              }
              else IF_CONSTEXPR(nDim == 3) {
                MFloat closestPoint[3] = {F0, F0, F0};
                MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
                a_oldLevelSetFunctionG(cellId, set) *=
                    computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set, m_sphereRadiusLimit);
              }
            } else {
              for(MInt setI = m_startSet; setI < m_noSets; setI++) {
                if(!m_maxLevelChange && m_computeSet_backup[setI] && !m_virtualSurgery) continue;
                MInt closestElement;
                IF_CONSTEXPR(nDim == 2) {
                  MFloat closestPoint[2] = {F0, F0};
                  MFloat refPoint[2] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1)};
                  a_oldLevelSetFunctionG(cellId, setI) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, setI, m_sphereRadiusLimit);
                }
                else IF_CONSTEXPR(nDim == 3) {
                  MFloat closestPoint[3] = {F0, F0, F0};
                  MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
                  a_oldLevelSetFunctionG(cellId, setI) *=
                      computeDistanceFromSTL(refPoint, &closestElement, closestPoint, setI, m_sphereRadiusLimit);
                }
              }
            }
          }
          break;
        }
      default: {
        mTerm(1, AT_, "Unknown constructFlag");
      }
    }
  }

  // b) reinitialization
  //---------------------------------------------------------------------------------------------------
  if(ConstructFlag == 0 || ConstructFlag == 2) {
    if(m_STLReinitMode == 2) return;
    IF_CONSTEXPR(nDim == 2) {
      if(m_STLReinitMode != 0) {
        computeNormalVectors();
        computeCurvature();
        return;
      }
    }

    MInt tmp_gReinitIterations = m_gReinitIterations;
    m_gReinitIterations =
        Context::getSolverProperty<MInt>("gReinitIterationsForGFieldFromSTL", m_solverId, AT_, &m_gReinitIterations);

    /*! \page propertiesLS
    \section gReinitIterationsForGFieldFromSTL
    <code>MInt FvCartesianSolver::m_gReinitIterations </code>\n
    default = <code>no default value</code>\n \n
    This property triggers the maximum number of reinitialization steps to reach the convergence criterion
    same functionality as gReinitIterations in lssolver.cpp:
    m_gReinitIterations reads in the property before calling the function levelSetConstrainedReinitialization
    after the "Reinit", return the old value to m_gReinitIterations
    possible values are:
    <ul>
    <li> Any non-negative integer value</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, ITERATION, STEPS </i>
  */

    if(m_STLReinitMode == 0 || m_STLReinitMode == 3 || m_STLReinitMode == 1) {
      levelSetConstrainedReinitialization(2, m_startSet, m_noSets, 1);
    }
    computeNormalVectors();
    computeCurvature();
    if(m_STLReinitMode == 0 || m_STLReinitMode == 1) {
      levelSetConstrainedReinitialization(2, m_startSet, m_noSets, 1);
    } else if(m_STLReinitMode == 3) {
      levelSetHighOrderConstrainedReinitialization(2, m_startSet, m_noSets, 1);
    }

    m_gReinitIterations = tmp_gReinitIterations;

    if(m_STLReinitMode == 0 || m_STLReinitMode == 3) {
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        if(a_isHalo(cellId)) continue;
        if(!c_isLeafCell(cellId)) continue;
        if(a_level(cellId) != a_maxGCellLevel()) continue;
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          const MInt nghbrId = c_neighborId(cellId, dir);
          if(nghbrId == -1 || nghbrId == cellId) {
            MFloat refPoint[3] = {F0, F0, F0};
            for(MInt i = 0; i < nDim; i++) {
              refPoint[i] = c_coordinate(cellId, i);
            }
            MInt closestElement;
            MFloat closestPoint[3] = {F0, F0, F0};
            for(MInt set = m_startSet; set < m_noSets; set++) {
              MFloat sign = determineLevelSetSignFromSTL(refPoint, set);
              if(m_levelSetSign[set] < 0) sign *= -1.0;
              MFloat phi = sign * computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);

              a_levelSetFunctionG(cellId, set) = phi;
            }
          }
        }
      }
    } else if(m_STLReinitMode == 1) {
      //      However, with Reinitialization on all G-cell-levels.
      // Timw: inserted, just as above, to allow for a Correction in Mode1.
      //      This avoids errors in the G-Field initialization!
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        if(!c_isLeafCell(cellId)) continue;
        if(a_level(cellId) != a_maxGCellLevel()) continue;
        if(a_isHalo(cellId)) continue;
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          MInt nghbrId = c_neighborId(cellId, dir);
          if(nghbrId == -1 || nghbrId == cellId) {
            MFloat refPoint[3] = {F0, F0, F0};
            for(MInt i = 0; i < nDim; i++) {
              refPoint[i] = c_coordinate(cellId, i);
            }
            MInt closestElement;
            MFloat closestPoint[3] = {F0, F0, F0};
            for(MInt set = m_startSet; set < m_noSets; set++) {
              MFloat sign = determineLevelSetSignFromSTL(refPoint, set);
              if(m_levelSetSign[set] < 0) sign *= -1.0;
              MFloat phi = sign * computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
              a_levelSetFunctionG(cellId, set) = phi;
            }
          }
        }
      }
    }
  }

  exchangeLevelSet();
}

//---------------------------------------------------------------------------

template <MInt nDim>
void LsCartesianSolver<nDim>::getContainingCellFromNeighbor(MInt body, MInt cellId, MFloat* xCoord, MFloat* xOld) {
  TRACE();

  MInt nghbrId = -1;
  MInt containingCell = -1;
  MInt containingDomain = -1;

  MInt bodyId;
  bodyId = m_bodiesToCompute[body];

  std::vector<MInt> dir;
  dir.resize(2 * nDim);

  std::vector<MFloat> vel;
  vel.resize(nDim);
  MFloat maxVel = -1;
  MFloat minVel = std::numeric_limits<MFloat>::max();
  MInt dirMax = -1;
  MInt dirMin = -1;

  rotateLevelSet(2, &vel[0], bodyId, xCoord, xOld, &m_semiLagrange_xRot_STL[body * nDim]);
  for(MInt i = 0; i < nDim; i++) {
    if(abs(vel[i]) > maxVel) {
      dirMax = i;
      maxVel = abs(vel[i]);
    }
    if(abs(vel[i]) < minVel) {
      dirMin = i;
      minVel = abs(vel[i]);
    }
  }
  for(MInt i = 0; i < nDim; i++) {
    if(i == dirMax) {
      if(vel[i] >= F0) {
        dir[0] = i * 2;
        dir[2 * nDim - 1] = i * 2 + 1;
      } else {
        dir[0] = i * 2 + 1;
        dir[2 * nDim - 1] = i * 2;
      }
    } else if(i == dirMin) {
      if(vel[i] >= F0) {
        dir[nDim - 1] = i * 2;
        dir[nDim] = i * 2 + 1;
      } else {
        dir[nDim - 1] = i * 2 + 1;
        dir[nDim] = i * 2;
      }
    } else {
      if(vel[i] >= F0) {
        dir[1] = i * 2;
        dir[4] = i * 2 + 1;
      } else {
        dir[1] = i * 2 + 1;
        dir[4] = i * 2;
      }
    }
  }


  MInt offset = 0;
  MInt oDir[6] = {1, 2, 0, 2, 0, 1};

  for(MInt n = 0; n < 2; n++) {
    // check neighbors
    for(MInt i = 0; i < nDim; i++) {
      nghbrId = c_neighborId(cellId, dir[i + offset]);
      if(nghbrId == -1) continue;

      containingCell = a_containingCell(nghbrId, body);
      if(containingCell > -1 && std::find(m_newCells.begin(), m_newCells.end(), nghbrId) == m_newCells.end()) {
        a_containingCell(cellId, body) = containingCell;

        if(!m_reconstructOldG) {
          containingDomain = a_containingDomain(nghbrId, body);
          a_containingDomain(cellId, body) = containingDomain;
        }

        return;
      }
    }
    // check diagonal neighbors
    for(MInt i = 0; i < nDim; i++) {
      if(!a_hasNeighbor(cellId, dir[i + offset])) continue;
      nghbrId = c_neighborId(cellId, dir[i + offset]);
      for(MInt j = 0; j < (nDim - 1); j++) {
        if(!a_hasNeighbor(nghbrId, dir[oDir[i * 2 + j] + offset])) continue;
        nghbrId = c_neighborId(nghbrId, dir[oDir[i * 2 + j] + offset]);

        containingCell = a_containingCell(nghbrId, body);
        if(containingCell > -1 && std::find(m_newCells.begin(), m_newCells.end(), nghbrId) == m_newCells.end()) {
          a_containingCell(cellId, body) = containingCell;
          if(!m_reconstructOldG) {
            containingDomain = a_containingDomain(nghbrId, body);
            a_containingDomain(cellId, body) = containingDomain;
          }

          return;
        }
      }
    }
    offset = nDim;
  }
  mTerm(1, AT_, "Get containingCell from neighbor did not work!");
}


//----------------------------------------------------------------------------

/**
 * \fn MInt LsCartesianSolver::determineLevelSetSignFromSTL
 * \brief this function checks if the "target" coordinates is inside(return 1) or outside
 *         (return -1) STL and return level set sign.
 *
 * \author: Sitthikrit Leckpool, Sep. 2010
 */

template <MInt nDim>
MInt LsCartesianSolver<nDim>::determineLevelSetSignFromSTL(MFloat* target, MInt set) {
  TRACE();

  if(m_noBodyBndryCndIds > 0) {
    ASSERT(m_noBodiesInSet[set] <= m_noBodyBndryCndIds,
           to_string(m_noBodiesInSet[set]) + " " + to_string(m_noBodyBndryCndIds));
  }

  IF_CONSTEXPR(nDim == 2) {
    if(m_geometry->pointIsInsideMBElements(target, m_bodyBndryCndIds, m_setToBodiesTable[set], m_noBodiesInSet[set]))
      return -1;
  }
  else {
    if(m_geometry->pointIsInsideMBElements2(target, m_bodyBndryCndIds, m_setToBodiesTable[set], m_noBodiesInSet[set]))
      return -1;
  }
  return 1;
}


//----------------------------------------------------------------------------
/** \brief this function returns the exact distance from the "target" point to the STL.
 *
 * \author: Sitthikrit Leckpool, Sep. 2010
 */

template <MInt nDim>
MFloat LsCartesianSolver<nDim>::computeDistanceFromSTL(MFloat* target, MInt* closestElement, MFloat* closestPoint,
                                                       MInt set, MFloat sphereRadiusFactor) {
  TRACE();

  // class for checking if "target" can have normal distance to the STL

  GeometryElement<nDim>* mbelements = m_geometry->m_mbelements->a;

  MFloat q[3] = {F0, F0, F0};
  for(MInt dim = 0; dim < nDim; dim++)
    q[dim] = target[dim];

  /*********** find the set of nearest geometrical elements to the target *********/
  std::vector<MInt> nodeList;
  MInt noNodes_set = 0;
  *closestElement = -1;
  const MFloat cellLength = c_cellLengthAtLevel(a_maxGCellLevel(set));
  const MFloat cellHalfLength = c_cellLengthAtLevel(a_maxGCellLevel(set) + 1);
  // Initialize mindist
  MFloat mindist = -1.0;

  MFloat enlargeFactor = 1.01;

  const MFloat radiusLimitSphere = sphereRadiusFactor * m_gBandWidth;

  while(noNodes_set == 0 && enlargeFactor < radiusLimitSphere)
  // there is at least one element,otherwise extend the range by enlargeFactor
  {
    // Get all triangles in currentCell (target)
    const MFloat radius = cellHalfLength * enlargeFactor * sqrt(nDim);

    nodeList.clear();
    m_geometry->getSphereIntersectionMBElements(q, radius, nodeList);
    for(MInt node = 0; node < (signed)nodeList.size(); node++) {
      for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
        const MInt body = m_setToBodiesTable[set][b];
        const MInt bcId = m_bodyBndryCndIds[body];
        if(m_geometry->mbelements[nodeList[node]].m_bndCndId == bcId) {
          nodeList[noNodes_set] = nodeList[node];
          noNodes_set++;
        }
      }
    }

    // Enlarge the boundingbox
    enlargeFactor *= 1.5;
  }

  if(noNodes_set == 0) {
    mindist = F2 * m_gBandWidth * cellLength;
    return mindist;
  }

  // loop over the set of nearest geometrical elements
  for(MInt inode = 0; inode < noNodes_set; inode++) {
    MInt e = nodeList[inode];
    GeometryElement<nDim>* el = &mbelements[e];

    // check if the target can have normal distance to the current elements
    MInt haveNormal = 0;

    MFloat transformationMatrix[3][3] = {{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
    haveNormal = checkNormal().PointInsideTriangle(el, q, transformationMatrix);
    if(haveNormal) {
      // Calculate distance by normal
      IF_CONSTEXPR(nDim == 2) {
        MFloat n[3] = {0, 0, 0};
        MFloat quv[3] = {0, 0, 0};
        MFloat v[3] = {0, 0, 0}; // use only one vertex of the element
        MFloat vuv[3] = {0, 0, 0};
        for(MInt i = 0; i < nDim; i++) {
          v[i] = mbelements[e].m_vertices[0][i];
          n[i] = mbelements[e].m_normal[i];
        }

        checkNormal().rotation(q, quv, transformationMatrix);
        checkNormal().rotation(v, vuv, transformationMatrix);

        MFloat dist = fabs(quv[1] - vuv[1]);
        MFloat point[2] = {F0, F0};
        if(n[0] * (q[0] - v[0]) + n[1] * (q[1] - v[1]) > F0) {
          point[0] = q[0] - n[0] * dist;
          point[1] = q[1] - n[1] * dist;
        } else {
          point[0] = q[0] + n[0] * dist;
          point[1] = q[1] + n[1] * dist;
        }

        if(mindist < 0.0) {
          mindist = dist;
          *closestElement = e;
          closestPoint[0] = point[0];
          closestPoint[1] = point[1];
        } else {
          if(mindist > dist) {
            mindist = dist;
            *closestElement = e;
            closestPoint[0] = point[0];
            closestPoint[1] = point[1];
          }
        }
      }
      else {
        MFloat n[3] = {0, 0, 0};
        MFloat v[3] = {0, 0, 0}; // use only one vertex of the element
        for(MInt i = 0; i < nDim; i++) {
          n[i] = mbelements[e].m_normal[i];
          v[i] = mbelements[e].m_vertices[0][i];
        }

        MFloat d = -(n[0] * v[0] + n[1] * v[1] + n[2] * v[2]);
        MFloat dist = fabs(n[0] * q[0] + n[1] * q[1] + n[2] * q[2] + d);
        MFloat point[3] = {F0, F0, F0};
        if(n[0] * (q[0] - v[0]) + n[1] * (q[1] - v[1]) + n[2] * (q[2] - v[2]) > F0) {
          point[0] = q[0] - n[0] * dist;
          point[1] = q[1] - n[1] * dist;
          point[2] = q[2] - n[2] * dist;
        } else {
          point[0] = q[0] + n[0] * dist;
          point[1] = q[1] + n[1] * dist;
          point[2] = q[2] + n[2] * dist;
        }
        if(mindist < 0.0) {
          mindist = dist;
          *closestElement = e;
          closestPoint[0] = point[0];
          closestPoint[1] = point[1];
          closestPoint[2] = point[2];
        } else {
          if(mindist > dist) {
            mindist = dist;
            *closestElement = e;
            closestPoint[0] = point[0];
            closestPoint[1] = point[1];
            closestPoint[2] = point[2];
          }
        }
      }
    } else {
      MInt noVertices;
      IF_CONSTEXPR(nDim == 2) {
        noVertices = 2;
        // Calculate distance by points
        for(MInt i = 0; i < noVertices; i++) {
          MFloat v[3] = {0, 0, 0};
          for(MInt j = 0; j < nDim; j++) {
            v[j] = mbelements[e].m_vertices[i][j];
          }
          MFloat dist = sqrt(pow(q[0] - v[0], 2) + pow(q[1] - v[1], 2) + pow(q[2] - v[2], 2));
          if(mindist < 0.0) {
            mindist = dist;
            *closestElement = e;
            closestPoint[0] = v[0];
            closestPoint[1] = v[1];
          } else {
            if(mindist > dist) {
              mindist = dist;
              *closestElement = e;
              closestPoint[0] = v[0];
              closestPoint[1] = v[1];
            }
          }
        }
      }
      else {
        noVertices = 3;
        for(MInt i = 0; i < noVertices; i++) {
          MFloat v1[3] = {0, 0, 0}, v2[3] = {0, 0, 0};
          for(MInt j = 0; j < nDim; j++) {
            v1[j] = mbelements[e].m_vertices[i][j];
            v2[j] = mbelements[e].m_vertices[(i + 1) % noVertices][j];
          }
          // edge normal distance
          MFloat vp[3] = {v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]};
          MFloat vq[3] = {q[0] - v1[0], q[1] - v1[1], q[2] - v1[2]};
          // normalize by p
          MFloat norm_vp = sqrt(vp[0] * vp[0] + vp[1] * vp[1] + vp[2] * vp[2]);
          MFloat vp_hat[3] = {vp[0] / norm_vp, vp[1] / norm_vp, vp[2] / norm_vp};
          MFloat vq_hat[3] = {vq[0] / norm_vp, vq[1] / norm_vp, vq[2] / norm_vp};
          // scaling factor
          MFloat I = vp_hat[0] * vq_hat[0] + vp_hat[1] * vq_hat[1] + vp_hat[2] * vq_hat[2];
          MFloat dist = 0.0;
          MFloat point[3] = {F0, F0, F0};
          if((0.0 <= I) && (I <= 1.0)) { // edge normal distance
            MFloat x_e[3] = {v1[0] + I * vp[0], v1[1] + I * vp[1], v1[2] + I * vp[2]};
            dist = sqrt(POW2(q[0] - x_e[0]) + POW2(q[1] - x_e[1]) + POW2(q[2] - x_e[2]));
            point[0] = x_e[0];
            point[1] = x_e[1];
            point[2] = x_e[2];
          } else { // point-to-vertex distance
            dist = sqrt(vq[0] * vq[0] + vq[1] * vq[1] + vq[2] * vq[2]);
            point[0] = v1[0];
            point[1] = v1[1];
            point[2] = v1[2];
          }

          if(mindist < 0.0) {
            mindist = dist;
            *closestElement = e;
            closestPoint[0] = point[0];
            closestPoint[1] = point[1];
            closestPoint[2] = point[2];
          } else {
            if(mindist > dist) {
              mindist = dist;
              *closestElement = e;
              closestPoint[0] = point[0];
              closestPoint[1] = point[1];
              closestPoint[2] = point[2];
            }
          }
        }
      }
    }
  }
  return mindist;
}

//----------------------------------------------------------------------------
/**
 * \fn void LsCartesianSolver::spatiallyAdaptiveCorrectionFromSTL()
 * \brief this function does a correction based on the curvature of the geometry.
             The high curvature regions are corrected, while the low curvature regions
             are left untouched. The number of the corrected is controlled
             by the multiple level thresholds.
             remarks: this works only with control 2 (moving STL).
 * \author  Sitthikrit Leckpool, Sep. 2010
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::spatiallyAdaptiveCorrectionFromSTL() {
  TRACE();

  //----------- prepare correction ---------//
  computeNormalVectorsAtFront();
  computeCurvature();
  // update control points
  m_gCtrlPnt.CtrlPnt2_UpdateAllNormalVector();
  m_gCtrlPnt.CtrlPnt2_Update();

  //----------- compute gradient of curvature ---------------//
  // the gradient of curvature is computed by the central difference.
  MInt maxNoG0Cells = 0;
  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    maxNoG0Cells = mMax(maxNoG0Cells, a_noG0Cells(set));
  }
  MFloatScratchSpace gradCurvature(maxNoG0Cells * m_maxNoSets, AT_, "gradCurvature");

  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      MInt nghbrL, nghbrR, cellId;
      gradCurvature[IDX_LSSET(id, set)] = F0;
      cellId = a_G0CellId(id, set);
      for(MInt i = 0; i < nDim; i++) {
        MFloat gradComponent = F0;
        nghbrL = c_neighborId(cellId, 2 * i);
        nghbrR = c_neighborId(cellId, 2 * i + 1);
        gradComponent = (a_curvatureG(nghbrR, set) - a_curvatureG(nghbrL, set)) / (2.0 * m_gCellDistance);
        gradCurvature[IDX_LSSET(id, set)] += (gradComponent * gradComponent);
      }
      gradCurvature[IDX_LSSET(id, set)] = sqrt(gradCurvature[IDX_LSSET(id, set)]);
    }
  }
  //---------- mark cells ----------//

  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;

    /*! \page propertiesLS
  \section levelSetMaxCorrectionThresholdLevels
  <code>MInt LsCartesianSolver::spatiallyAdaptiveCorrectionFromSTL::maxThresholdLevels</code>\n
  default = <code>3</code>\n
  spatiallyAdaptiveCorrectionFromSTL does a correction
  based on the curvature of the geometry.
  The high curvature regions are corrected, while the low curvature regions
  are left untouched. The number of the corrected is controlled
  by the multiple level thresholds.
  remarks: this works only with control 2 (moving STL).\n
  Possible values are:
  <ul>
    <li>1 to ?</li>
  </ul>
  Keywords: <i>LEVEL SET, STL, CURVATURE, GEOMETRY</i>
*/
    MInt maxThresholdLevels = 3;
    maxThresholdLevels =
        Context::getSolverProperty<MInt>("levelSetMaxCorrectionThresholdLevels", m_solverId, AT_, &maxThresholdLevels);
    MIntScratchSpace marker(a_noG0Cells(set), AT_, "marker");
    list<MInt> m_unmarked;
    list<MInt>::iterator it_unmarked;
    // initialize marker
    for(MInt c = 0; c < a_noG0Cells(set); c++) {
      marker[c] = 0;
      m_unmarked.push_back(c);
    }

    // compute thresholds
    MInt curThresholdLevel = 1;
    while((curThresholdLevel <= maxThresholdLevels) && (m_unmarked.size() != 0)) {
      MFloat avgCurvature = 0.0;
      it_unmarked = m_unmarked.begin();
      // average gradient by RMS
      while(it_unmarked != m_unmarked.end()) {
        avgCurvature += (gradCurvature[IDX_LSSET(*it_unmarked, set)] * gradCurvature[IDX_LSSET(*it_unmarked, set)]);
        it_unmarked++;
      }
      avgCurvature /= (MFloat)m_unmarked.size();
      avgCurvature = sqrt(avgCurvature);
      // mark cells for reconstruction
      it_unmarked = m_unmarked.begin();
      while(it_unmarked != m_unmarked.end()) {
        if(gradCurvature[IDX_LSSET(*it_unmarked, set)] > avgCurvature) {
          marker[*it_unmarked] = curThresholdLevel;
          *it_unmarked = a_noG0Cells(set);
        }
        it_unmarked++;
      }
      // clean up marked cells from m_unmarked list
      m_unmarked.remove(a_noG0Cells(set));
      // increase threshold
      curThresholdLevel++;
    }
    //-------------- correct marked cells ------------------//
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      MInt cellId = a_G0CellId(id, set);
      m_d[IDX_LSSET(cellId, set)] = F0;
      m_correction[IDX_LSSET(cellId, set)] = F0;
      if(marker[id] > 0) {
        // Sitti Implementation - I (claudia) don't understand why to do this so complcated -> dangerous!
        // compute distance
        //           MFloat q[] = { 0.0, 0.0, 0.0 };
        //           // extract G0 point from the interface cell
        //           for( MInt dim=0; dim<nDim; dim++ )
        //             q[dim] = c_coordinate( cellId ,  dim )
        //             + a_levelSetFunctionG( cellId , set)
        //             * a_normalVectorG( cellId ,  dim , set );
        //         MInt closestElement;
        //       MFloat closestPoint[3] = {F0,F0,F0};
        //         MFloat dist_q = computeDistanceFromSTL( q, &closestElement, closestPoint, set );
        //         m_d[ IDX_LSSET(cellId, set) ] = dist_q;
        //         // determine sign
        //         MInt sign = determineLevelSetSignFromSTL( q, set );
        //         m_correction[ IDX_LSSET(cellId, set) ] = sign;
        //
        //         // correction
        //         a_levelSetFunctionG( cellId , set)  = a_levelSetFunctionG( cellId , set)  + m_levelSetSign[set] *
        //         (MFloat)sign*dist_q;

        // better use this:

        MInt closestElement;
        MFloat closestPoint[3] = {F0, F0, F0};
        MFloat refPoint[3] = {c_coordinate(cellId, 0), c_coordinate(cellId, 1), c_coordinate(cellId, 2)};
        MFloat dist_q = computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
        m_d[IDX_LSSET(cellId, set)] = dist_q;
        // determine sign
        MInt sign = determineLevelSetSignFromSTL(refPoint, set);
        m_correction[IDX_LSSET(cellId, set)] = sign;
        // correction
        a_levelSetFunctionG(cellId, set) = m_levelSetSign[set] * (MFloat)sign * dist_q;
      }
    }
  }
}


//----------------------------------------------------------------------------


// based on all level sets!
/** \brief regrid level set
 *
 * \author Daniel Hartmann
 * \date July 2007
 *
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::regridLevelSet() {
  TRACE();

  MBool regrid = false;

  // check whether the G grid needs to be regridded
  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      if(a_isHalo(a_bandCellId(id, set))) continue;
      if(a_regridTriggerG(a_bandCellId(id, set))) {
        regrid = true;
        m_log << a_bandCellId(id, set) << " with set " << set << " caused regridding " << endl;
        break;
      }
    }
    if(regrid) break;
  }

  return regrid;
}


// ----------------------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetReinitialization(MInt mode /*==1*/) {
  TRACE();
  //#define orderOfAccuracyTest
  MInt startSet;
  MInt endSet;
  MString reinitMethod = m_reinitMethod;

  switch(mode) {
    case 0:
      startSet = 0;
      endSet = 1;
      reinitMethod = m_gapReinitMethod;
      break;
    default:
      startSet = m_startSet;
      endSet = m_noSets;
      break;
  }

  if(m_semiLagrange) {
    ASSERT(reinitMethod == m_gapReinitMethod, "");
    ASSERT(m_guaranteeReinit, "");
  }


#ifndef orderOfAccuracyTest
  // compute the reinitialization sensor
  if(!mode || ((levelSetReinitializationTrigger() || m_guaranteeReinit) && !m_maintainOuterBandLayers)) {
    if(mode && globalTimeStep % m_reinitInterval != 0) return;

    // correction scheme
    if(m_GFieldInitFromSTL && m_GWithReConstruction && mode) {
      //       cerr << "correction" << endl;
      spatiallyAdaptiveCorrectionFromSTL();
    }

    //     cerr << "reinitialize..."<< endl;

    // reinitialize
    switch(string2enum(reinitMethod)) {
      case CR1: {
        levelSetConstrainedReinitialization(1, startSet, endSet, mode);
        break;
      }
      case CR2: {
        levelSetConstrainedReinitialization(2, startSet, endSet, mode);
        break;
      }
      case HCR1: {
        levelSetHighOrderConstrainedReinitialization(1, startSet, endSet, mode);
        break;
      }
      case HCR2: {
        levelSetHighOrderConstrainedReinitialization(2, startSet, endSet, mode);
        break;
      }
      case HCR2_LIMITED: {
        levelSetHighOrderConstrainedReinitialization(4, startSet, endSet, mode);
        break;
      }
      case HCR2_FULLREINIT: {
        levelSetHighOrderConstrainedReinitialization(6, startSet, endSet, mode);
        break;
      }
      case no:
      default: {
        break;
      }
    }
  } else {
    ASSERT(!m_semiLagrange, "");

    // maintain outer band cell layers
    switch(string2enum(reinitMethod)) {
      case CR1:
      case CR2: {
        maintainOuterBandLayers(1, startSet, endSet);
        break;
      }
      case HCR1:
      case HCR2:
      case HCR2_FULLREINIT: {
        maintainOuterBandLayers(5, startSet, endSet);
        break;
      }
      case no:
      default: {
        break;
      }
    }
  }
  if(domainId() == 0) {
    // write out info
    FILE* datei;
    stringstream reinitFile;
    reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
    datei = fopen((reinitFile.str()).c_str(), "a+");
    fprintf(datei, "\n");
    fclose(datei);
  }

#else
  MInt cellId;
  MFloat radius;

  switch(string2enum(reinitMethod)) {
    case CR1: {
      levelSetConstrainedReinitialization(1, startSet, endSet, mode);
      break;
    }
    case CR2: {
      levelSetConstrainedReinitialization(2, startSet, endSet, mode);
      break;
    }
    case HCR1: {
      levelSetHighOrderConstrainedReinitialization(1, startSet, endSet, mode);
      break;
    }
    case HCR2: {
      levelSetHighOrderConstrainedReinitialization(2, startSet, endSet, mode);
      break;
    }
    case HCR2_FULLREINIT: {
      levelSetHighOrderConstrainedReinitialization(6, startSet, endSet, mode);
      break;
    }
    case no:
    default: {
      break;
    }
  }

  // compute the error
  /// circle problem
  MFloat error;
  error = F0;
  MInt set = 0; // defined only for zero level set (m_noSets = 1)
  for(MInt id = 0; id < a_noG0Cells(set); id++) {
    cellId = a_G0CellId(id, set);
    error +=
        ABS(sqrt(POW2(c_coordinate(cellId, 0)) + POW2(c_coordinate(cellId, 1))) - 3 + a_levelSetFunctionG(cellId, set));
  }
  error /= a_noG0Cells(set);
  cerr << "Error after reinitialization: " << error << endl;
  FILE* datei;
  if(globalTimeStep % 1 == 0) {
    datei = fopen("ErrorInGamma", "a+");
    fprintf(datei, "%d", globalTimeStep);
    fprintf(datei, " %f", error);
    fprintf(datei, "\n");
    fclose(datei);
  }
  ///
#endif
}


//-----------------------------------------------------------------------------


// based on zeroth level-set function
template <MInt nDim>
MBool LsCartesianSolver<nDim>::levelSetReinitializationTrigger() {
  TRACE();

  MInt cellId;
  MFloat temp = F0, avgDeviation = F0, maxDeviation = F0;
  MInt globalNoG0Cells = F0;
  MInt set = 0;
  //---

  for(MInt id = 0; id < a_noG0Cells(set); id++) {
    cellId = a_G0CellId(id, set);
    if(a_isHalo(cellId)) continue;
    temp = F0;
    for(MInt i = 0; i < nDim; i++)
      temp += POW2(a_levelSetFunctionSlope(cellId, i, set));
    temp = POW2(sqrt(temp) - F1);

    avgDeviation += temp;
    maxDeviation = mMax(maxDeviation, sqrt(temp));
    globalNoG0Cells++;
  }

  MPI_Allreduce(MPI_IN_PLACE, &globalNoG0Cells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "globalNoG0Cells");
  MPI_Allreduce(MPI_IN_PLACE, &avgDeviation, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "avgDeviation");
  MPI_Allreduce(MPI_IN_PLACE, &maxDeviation, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "maxDeviation");
  avgDeviation = sqrt(avgDeviation / (MFloat)globalNoG0Cells);


  if(domainId() == 0) {
    // write out info
    FILE* datei;
    stringstream reinitFile;
    reinitFile << "Reinitialization_" << m_solverId << "_" << domainId();
    datei = fopen((reinitFile.str()).c_str(), "a+");
    fprintf(datei, " %d", globalTimeStep);
    fprintf(datei, "  %f", time());
    fprintf(datei, "  %f", globalTimeStep * timeStep());
    fprintf(datei, "  %f", avgDeviation);
    fprintf(datei, "  %f", maxDeviation);
    fprintf(datei, "  %f", m_reinitThreshold);
    fprintf(datei, "  %d", globalNoG0Cells);

    fclose(datei);
  }

  if(fabs(maxDeviation) > m_reinitThreshold && fabs(avgDeviation) > m_reinitThresholdAvg)
    return true;
  else
    return false;
}

// ----------------------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::setBandNewArrivals(MInt computingSet) {
  TRACE();

  MInt startSet = 0;
  MInt endSet = m_noSets;
  if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
    ASSERT(m_changedSet[computingSet], "");
  }

  //------------------


  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // may not be skipped as it is levelSet-value dependand

    // reset all cells at the tube boundary
    for(MInt id = 0; id < a_noBandBndryCells(set); id++) {
      const MInt cellId = a_bandBndryCellId(id, set);
      if(a_isHalo(cellId)) {
        continue;
      }
      if(a_level(cellId) > a_maxGCellLevel(0)) continue;

      for(MInt d = 0; d < m_noDirs; d++) {
        if(a_hasNeighbor(cellId, d) == 0) {
          continue;
        }
        if(!a_inBandG(c_neighborId(cellId, d), set)) {
          continue;
        }
        if(a_isGBoundaryCellG(c_neighborId(cellId, d), set)) {
          continue;
        }
        if(a_levelSetFunctionG(cellId, set) > F0) {
          a_levelSetFunctionG(cellId, set) = mMin(a_levelSetFunctionG(cellId, set),
                                                  a_levelSetFunctionG(c_neighborId(cellId, d), set) + m_gCellDistance);
        } else {
          a_levelSetFunctionG(cellId, set) = mMax(a_levelSetFunctionG(cellId, set),
                                                  a_levelSetFunctionG(c_neighborId(cellId, d), set) - m_gCellDistance);
        }
      }
    }
  }

  exchangeLevelSet();
}

// ----------------------------------------------------------------------------------------


// might be speeded up by more intelligent coding (for multiple level-set functions)
/** \brief updates the level set function on lower grid levels
 *
 * \author Daniel Hartmann, Stephan Schlimpert
 * \date unknown, September 2011
 *
 * last change: function optimized
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::updateLowerGridLevels(MInt computingSet) {
  TRACE();

  MInt startSet = 0;
  MInt endSet = m_noSets;
  // if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  // TODO labels:LS update testcases, so that the uncommented version above passes!
  //      this is due to restartVariables on lowerGridLevels not matching otherwise!
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
  }

  MBoolScratchSpace skip(a_noCells(), AT_, "skip");
  MIntScratchSpace tempCells(a_noCells(), AT_, "tempCells");

  //---


  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // this happens before, determineG0Cells, thus !m_changedSet can not be skipped!

    for(MInt cell = 0; cell < a_noCells(); cell++) {
      skip.p[cell] = false;
    }

    // average the level set function to parents of leaf cells
    MInt noTempCells = 0;
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      const MInt cellId = a_bandCellId(id, set);
      const MInt parentId = c_parentId(cellId);
      if(parentId > -1) {
        if(!skip.p[parentId]) {
          tempCells.p[noTempCells] = parentId;
          noTempCells++;
          a_levelSetFunctionG(parentId, set) = F0;
          MInt noChildren = c_noChildren(parentId);
          if(!(noChildren > 0)) {
            cerr << "parent: " << parentId << endl;
            cerr << "gcell: " << cellId << endl;
            mTerm(1, AT_, "this should not occur: parent g cell with m_noChildIdsG == 0...");
          }
          for(MInt child = 0; child < IPOW2(nDim); child++) {
            if(c_childId(parentId, child) < 0) continue;
            a_levelSetFunctionG(parentId, set) += a_levelSetFunctionG(c_childId(parentId, child), set);
          }
          a_levelSetFunctionG(parentId, set) /= (MFloat)noChildren;
          skip.p[parentId] = true;
        }
      }
    }

    // average the level set function all the way down
    // overwrite tempCells
    while(noTempCells > 0) {
      const MInt listSize = noTempCells;
      noTempCells = 0;
      for(MInt cell = 0; cell < listSize; cell++) {
        const MInt gCellId = tempCells.p[cell];
        const MInt parentId = c_parentId(gCellId);
        if(parentId > -1) {
          if(!skip.p[parentId]) {
            tempCells.p[noTempCells] = parentId;
            noTempCells++;
            a_levelSetFunctionG(parentId, set) = F0;
            MInt noChildren = c_noChildren(parentId);
            if(!(noChildren > 0)) {
              cerr << parentId << endl;
              cerr << "gcell: " << gCellId << endl;
              mTerm(1, AT_, "this should not occur: parent g cell with m_noChildIdsG == 0...");
            }
            for(MInt child = 0; child < IPOW2(nDim); child++) {
              if(c_childId(parentId, child) < 0) continue;
              a_levelSetFunctionG(parentId, set) += a_levelSetFunctionG(c_childId(parentId, child), set);
            }
            a_levelSetFunctionG(parentId, set) /= (MFloat)noChildren;
            skip.p[parentId] = true;
          }
        }
      }
    }
  }
}

template <MInt nDim>
void LsCartesianSolver<nDim>::updateAllLowerGridLevels(MInt computingSet) {
  TRACE();

  MInt startSet = 0;
  MInt endSet = m_noSets;
  // if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  // TODO labels:LS update testcases, so that the uncommented version above passes!
  //      this is due to restartVariables on lowerGridLevels not matching otherwise!
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
  }

  MBoolScratchSpace skip(a_noCells(), AT_, "skip");
  MIntScratchSpace tempCells(a_noCells(), AT_, "tempCells");

  //---


  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    // this happens before, determineG0Cells, thus !m_changedSet can not be skipped!

    for(MInt cell = 0; cell < a_noCells(); cell++) {
      skip.p[cell] = false;
    }

    // average the level set function to parents of leaf cells
    MInt noTempCells = 0;
    for(MInt gCellId = 0; gCellId < a_noCells(); gCellId++) {
      if(a_isHalo(gCellId)) continue;
      if(!c_isLeafCell(gCellId)) continue;
      MInt parentId = c_parentId(gCellId);
      if(parentId > -1) {
        if(!skip.p[parentId]) {
          tempCells.p[noTempCells] = parentId;
          noTempCells++;
          a_levelSetFunctionG(parentId, set) = F0;
          MInt noChildren = c_noChildren(parentId);
          if(!(noChildren > 0)) {
            cerr << "parent: " << parentId << endl;
            cerr << "gcell: " << gCellId << endl;
            mTerm(1, AT_, "this should not occur: parent g cell with m_noChildIdsG == 0...");
          }
          for(MInt child = 0; child < IPOW2(nDim); child++) {
            if(c_childId(parentId, child) < 0) continue;
            a_levelSetFunctionG(parentId, set) += a_levelSetFunctionG(c_childId(parentId, child), set);
          }
          a_levelSetFunctionG(parentId, set) /= (MFloat)noChildren;
          skip.p[parentId] = true;
        }
      }
    }

    // average the level set function all the way down
    // overwrite tempCells
    while(noTempCells > 0) {
      MInt listSize = noTempCells;
      noTempCells = 0;
      for(MInt cell = 0; cell < listSize; cell++) {
        MInt gCellId = tempCells.p[cell];
        MInt parentId = c_parentId(gCellId);
        if(parentId > -1) {
          if(!skip.p[parentId]) {
            tempCells.p[noTempCells] = parentId;
            noTempCells++;
            a_levelSetFunctionG(parentId, set) = F0;
            MInt noChildren = c_noChildren(parentId);
            if(!(noChildren > 0)) {
              cerr << parentId << endl;
              cerr << "gcell: " << gCellId << endl;
              mTerm(1, AT_, "this should not occur: parent g cell with m_noChildIdsG == 0...");
            }
            for(MInt child = 0; child < IPOW2(nDim); child++) {
              if(c_childId(parentId, child) < 0) continue;
              a_levelSetFunctionG(parentId, set) += a_levelSetFunctionG(c_childId(parentId, child), set);
            }
            a_levelSetFunctionG(parentId, set) /= (MFloat)noChildren;
            skip.p[parentId] = true;
          }
        }
      }
    }
  }
}


// ----------------------------------------------------------------------------------------


/** \brief sets properties:
 *      m_isGBndryCell[ set ]   - set
 *      m_isInBand       - set
 *      m_inShadowLayer   - set
 *      m_isInBand       - false
 *      m_inShadowLayer   - false
 *
 * computingSet modes:
 * -1 : (default) uses all sets and uses previous band-cells to search for G0-Cells
 *  0 : only determines the collected levelsets G0-Cells, but uses all cells for the search
 *
 * \authors Daniel Hartmann, Stephan Schlimpert, Claudia Guenther
 * \date July 2007, June 2011, June 2012
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::determineG0Cells(MInt computingSet) {
  TRACE();
  MInt startSet = 0;
  MInt endSet = m_noSets;
  if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
    ASSERT(m_changedSet[computingSet], "");
  }
  const MInt noSets = endSet - startSet;


  // if computingSet is explicitly set to zero, reset noBandcells
  // default is: computingSet = -1, zero is only used for operations on the collected level-set function
  if(computingSet == 0) m_bandCells[computingSet].clear();
  //---end of initialization

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    if(computingSet < 0 && m_levelSetMb) m_changedSet[set] = false;

    // reset g cut cell counter
    std::vector<MInt>().swap(m_G0Cells[set]);

    if(a_noBandCells(set) == 0) {
      // use all leaf-cells
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        a_inBandG(cellId, set) = false;
        a_isGBoundaryCellG(cellId, set) = false;
        a_isGZeroCell(cellId, set) = false;
      }
      for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
        if(!c_isLeafCell(cellId)) continue;
        if(approx(a_levelSetFunctionG(cellId, set), F0, MFloatEps) && !a_isHalo(cellId)) {
          a_inBandG(cellId, set) = true;
          a_isGBoundaryCellG(cellId, set) = true;
          m_G0Cells[set].push_back(cellId);
          a_isGZeroCell(cellId, set) = true;
          if(!a_wasGZeroCell(cellId, set)) m_changedSet[set] = true;

        } else {
          for(MInt d = 0; d < m_noDirs; d++) {
            MInt parentId = cellId;
            // consider level-jumps between neighbors:
            // important if working with very small geometries embedded in a huge domain
            // -> huge level differences
            MInt nghbrId = -1;
            if(!a_hasNeighbor(cellId, d)) {
              while(true) {
                parentId = c_parentId(parentId);
                if(parentId == -1) break;
                if(a_hasNeighbor(parentId, d)) {
                  if(c_noChildren(c_neighborId(parentId, d)) == 0) break;
                }
              }
              if(parentId == -1) continue;
              nghbrId = c_neighborId(parentId, d);
            } else {
              nghbrId = c_neighborId(cellId, d);
            }
            if(nghbrId < 0) continue;
            if(!c_isLeafCell(nghbrId)) continue;
            if(a_levelSetFunctionG(nghbrId, set) * a_levelSetFunctionG(cellId, set) < F0) {
              if(!a_isGZeroCell(cellId, set) && !a_isHalo(cellId)) {
                a_inBandG(cellId, set) = true;
                a_isGBoundaryCellG(cellId, set) = true;
                m_G0Cells[set].push_back(cellId);
                a_isGZeroCell(cellId, set) = true;
                if(!a_wasGZeroCell(cellId, set)) m_changedSet[set] = true;
              }
              if(parentId != cellId && !a_isGZeroCell(nghbrId, set) && !a_isHalo(nghbrId)) {
                a_inBandG(nghbrId, set) = true;
                a_isGBoundaryCellG(nghbrId, set) = true;
                m_G0Cells[set].push_back(nghbrId);
                a_isGZeroCell(nghbrId, set) = true;
                if(!a_wasGZeroCell(nghbrId, set)) m_changedSet[set] = true;
              }
            }
          }
        }
      }
    } else {
      // use only previous band-cells, which must be at the highest refinement level!
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        const MInt cellId = a_bandCellId(id, set);
        a_inBandG(cellId, set) = false;
        a_isGBoundaryCellG(cellId, set) = false;
        a_isGZeroCell(cellId, set) = false;
        if(a_isHalo(cellId)) continue;
        if(approx(a_levelSetFunctionG(cellId, set), F0, MFloatEps)) {
          a_inBandG(cellId, set) = true;
          a_isGBoundaryCellG(cellId, set) = true;
          m_G0Cells[set].push_back(cellId);
          a_isGZeroCell(cellId, set) = true;
          if(!a_wasGZeroCell(cellId, set)) m_changedSet[set] = true;
        } else {
          for(MInt d = 0; d < m_noDirs; d++) {
            if(!a_hasNeighbor(cellId, d)) continue;
            if(!c_isLeafCell(c_neighborId(cellId, d))) continue;
            if((a_levelSetFunctionG(c_neighborId(cellId, d), set) * a_levelSetFunctionG(cellId, set) < F0)) {
              a_inBandG(cellId, set) = true;
              a_isGBoundaryCellG(cellId, set) = true;
              m_G0Cells[set].push_back(cellId);
              a_isGZeroCell(cellId, set) = true;
              if(!a_wasGZeroCell(cellId, set)) m_changedSet[set] = true;
              break;
            }
          }
        }
      }
    }
    // determine the internal level set front cells
    a_internalBandLayer(0, set) = 0;
    std::vector<MInt>().swap(m_internalBandCells[set]);
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      if(!a_isHalo(a_G0CellId(id, set))) {
        a_internalBandLayer(0, set)++;
        m_internalBandCells[set].push_back(a_G0CellId(id, set));
      }
    }
  }

  // exchange G0 cells for all sets with other domains
  MBoolScratchSpace tmp_data(a_noCells(), noSets, AT_, "tmp_data");

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noWindowCells(i); j++) {
      MInt dataSet = 0;
      for(MInt set = startSet; set < endSet; set++) {
        tmp_data(windowCellId(i, j), dataSet) = a_isGZeroCell(windowCellId(i, j), set);
        dataSet++;
      }
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
      MInt dataSet = 0;
      for(MInt set = startSet; set < endSet; set++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        tmp_data(windowId, dataSet) = a_isGZeroCell(windowId, set);
        dataSet++;
      }
    }
  }

  exchangeDataLS(&tmp_data(0, 0), noSets);

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      const MInt haloCell = haloCellId(i, j);
      MInt dataSet = 0;
      for(MInt set = startSet; set < endSet; set++) {
        if(!a_isGZeroCell(haloCell, set)) {
          a_isGZeroCell(haloCell, set) = tmp_data(haloCellId(i, j), dataSet);
          if(a_isGZeroCell(haloCell, set)) {
            m_G0Cells[set].push_back(haloCell);
            a_inBandG(haloCell, set) = true;
            a_isGBoundaryCellG(haloCell, set) = true;
            if(!a_wasGZeroCell(haloCell, set)) m_changedSet[set] = true;
          }
        }
        dataSet++;
      }
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      const MInt haloCell = grid().azimuthalHaloCell(i, j);
      MInt dataSet = 0;
      for(MInt set = startSet; set < endSet; set++) {
        if(!a_isGZeroCell(haloCell, set)) {
          a_isGZeroCell(haloCell, set) = tmp_data(haloCell, dataSet);
          if(a_isGZeroCell(haloCell, set)) {
            m_G0Cells[set].push_back(haloCell);
            a_inBandG(haloCell, set) = true;
            a_isGBoundaryCellG(haloCell, set) = true;
            if(!a_wasGZeroCell(haloCell, set)) {
              m_changedSet[set] = true;
            }
          }
        }
        dataSet++;
      }
    }
  }

  // rebuild G0 cells for phi^0 as set union of G0 cells of phi^i
  if(m_buildCollectedLevelSetFunction && m_determineG0CellsMode && computingSet <= 0) {
    startSet = 1;
    MInt changeSet = 0;
    // reset G0 cells of phi^0
    for(MInt id = 0; id < a_noG0Cells(changeSet); id++) {
      MInt cellId = a_G0CellId(id, changeSet);
      a_inBandG(cellId, changeSet) = false;
      a_isGBoundaryCellG(cellId, changeSet) = false;
      a_isGZeroCell(cellId, changeSet) = false;
    }
    std::vector<MInt>().swap(m_G0Cells[changeSet]);
    for(MInt set = startSet; set < endSet; set++) {
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        MInt cellId = a_G0CellId(id, set);
        if(!a_isGZeroCell(cellId, changeSet)) {
          a_inBandG(cellId, changeSet) = true;
          a_isGBoundaryCellG(cellId, changeSet) = true;
          m_G0Cells[changeSet].push_back(cellId);
          a_isGZeroCell(cellId, changeSet) = true;
          if(!a_wasGZeroCell(cellId, changeSet)) m_changedSet[changeSet] = true;
        }
      }
    }
  }

  // NOTE: the changedSet information is currently only used to skip the identifyBodies call!
  // TODO labels:LS,totest check why the determineBandCells can currently not be skipped if the G0-cells stay the same
  if(m_levelSetMb) {
    for(MInt set = startSet; set < endSet; set++) {
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        if(a_wasGZeroCell(cellId, set) && !a_isGZeroCell(cellId, set)) {
          m_changedSet[set] = true;
        }
      }
    }

    MIntScratchSpace changedSet(m_noSets, AT_, "changedSet");
    for(MInt set = startSet; set < endSet; set++) {
      if(m_changedSet[set])
        changedSet[set] = 1;
      else
        changedSet[set] = -1;
    }

    // exchange m_changedSet-property
    MPI_Allreduce(MPI_IN_PLACE, &changedSet[0], m_noSets, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "changedSet[0]");

    for(MInt set = startSet; set < endSet; set++) {
      if(changedSet[set] > 0) {
        m_changedSet[set] = true;
      } else {
        m_changedSet[set] = false;
      }
      if(globalTimeStep <= 0) {
        m_changedSet[set] = true;
        m_changedSet[0] = true;
      }
    }

    if(computingSet >= 0) ASSERT(m_changedSet[computingSet], "");
  }
}


// ----------------------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::determineBandCells(MInt computingSet) {
  TRACE();
  MInt startSet = 0;
  MInt endSet = m_noSets;
  if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
    ASSERT(m_changedSet[computingSet], "");
  }

  MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

  MIntScratchSpace temp(a_noCells(), AT_, "temp");
  MIntScratchSpace lastLayer(a_noCells(), AT_, "lastLayer");

  //---

  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;

    // if( computingSet < 0 && !m_changedSet[set]  &&
    //  (m_levelSetMb)  && !m_LsRotate) continue;

    // reset band cell counters
    std::vector<MInt>().swap(m_bandBndryCells[set]);

    // put the level set G0-cells into the band
    a_bandLayer(0, set) = a_noG0Cells(set);
    MInt layerCount = 1;
    MInt cnt = 0;

    std::vector<MInt>().swap(m_bandCells[set]);
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      m_bandCells[set].push_back(a_G0CellId(id, set));
    }

    // only then extend the first layer
    // by adding all existing neighbors which are not Halo-Cells to the band and lastLayer
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      const MInt cellId = a_G0CellId(id, set);
      a_isGBoundaryCellG(cellId, set) = false;
      // activate all neighbors
      for(MInt dir = 0; dir < m_noDirs; dir++) {
        const MInt nghbrId = c_neighborId(cellId, dir);
        if(nghbrId == -1) continue;
        if(a_isHalo(nghbrId)) continue;
        if(a_inBandG(nghbrId, set)) continue;
        lastLayer.p[cnt++] = nghbrId;
        a_inBandG(nghbrId, set) = true;
        m_bandCells[set].push_back(nghbrId);
      }
    }

    // exchange band cells for the current set with other domains
    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noWindowCells(i); j++) {
        tmp_data(windowCellId(i, j)) = a_inBandG(windowCellId(i, j), set);
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        tmp_data(windowId) = a_inBandG(windowId, set);
      }
    }

    exchangeDataLS(&tmp_data(0), 1);

    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noHaloCells(i); j++) {
        const MInt halocell = haloCellId(i, j);
        if(!a_inBandG(halocell, set)) {
          a_inBandG(halocell, set) = tmp_data(halocell);
          if(a_inBandG(halocell, set)) {
            m_bandCells[set].push_back(halocell);
            lastLayer.p[cnt++] = halocell;
          }
        }
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        const MInt haloCell = grid().azimuthalHaloCell(i, j);
        if(!a_inBandG(haloCell, set)) {
          a_inBandG(haloCell, set) = tmp_data(haloCell);
          if(a_inBandG(haloCell, set)) {
            m_bandCells[set].push_back(haloCell);
            lastLayer.p[cnt++] = haloCell;
          }
        }
      }
    }

    a_bandLayer(layerCount, set) = a_noBandCells(set);

    // extend all other layers
    // no add all layers iterativly
    while(layerCount < m_gBandWidth) {
      MInt tempCnt = 0;
      for(MInt c = 0; c < cnt; c++) {
        a_isGBoundaryCellG(lastLayer.p[c], set) = false;
        // activate all neighbors
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          const MInt nghbrId = c_neighborId(lastLayer.p[c], dir);
          if(nghbrId == -1) continue;
          if(a_isHalo(nghbrId)) continue;
          if(a_inBandG(nghbrId, set)) continue;
          temp.p[tempCnt++] = nghbrId;
          a_inBandG(nghbrId, set) = true;
          m_bandCells[set].push_back(nghbrId);
        }
      }


      // exchange band cells for all sets with other domains
      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noWindowCells(i); j++) {
          tmp_data(windowCellId(i, j)) = a_inBandG(windowCellId(i, j), set);
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
          MInt windowId = grid().azimuthalWindowCell(i, j);
          tmp_data(windowId) = a_inBandG(windowId, set);
        }
      }

      exchangeDataLS(&tmp_data(0), 1);

      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noHaloCells(i); j++) {
          const MInt halocell = haloCellId(i, j);
          if(!a_inBandG(halocell, set)) {
            a_inBandG(halocell, set) = tmp_data(halocell);
            if(a_inBandG(halocell, set)) {
              m_bandCells[set].push_back(halocell);
              temp.p[tempCnt++] = halocell;
            }
          }
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
          const MInt haloCell = grid().azimuthalHaloCell(i, j);
          if(!a_inBandG(haloCell, set)) {
            a_inBandG(haloCell, set) = tmp_data(haloCell);
            if(a_inBandG(haloCell, set)) {
              m_bandCells[set].push_back(haloCell);
              temp.p[tempCnt++] = haloCell;
            }
          }
        }
      }

      layerCount++;
      a_bandLayer(layerCount, set) = a_noBandCells(set);
      cnt = tempCnt;
      for(MInt c = 0; c < tempCnt; c++) {
        lastLayer.p[c] = temp.p[c];
      }
    }


    // set the last-layer as a_isGBoundaryCellG!
    // determine the G boundary cells
    for(MInt c = 0; c < cnt; c++) {
      a_isGBoundaryCellG(lastLayer.p[c], set) = true;
      m_bandBndryCells[set].push_back(lastLayer.p[c]);
    }

    // set the internal band layers
    // only relevant for reinitialization!
    for(MInt layer = 1; layer < m_gBandWidth + 1; layer++) {
      a_internalBandLayer(layer, set) = a_internalBandLayer(layer - 1, set);
      for(MInt id = a_bandLayer(layer - 1, set); id < a_bandLayer(layer, set); id++) {
        if(!a_isHalo(a_bandCellId(id, set))) {
          a_internalBandLayer(layer, set)++;
          m_internalBandCells[set].push_back(a_bandCellId(id, set));
        }
      }
    }
  }
}


// ----------------------------------------------------------------------------------------


/** \brief updates the boundary cell list
 *
 * \author Daniel Hartmann, Claudia Guenther, Stephan Schlimpert
 * \date unknown, Jan 2014
 *
 * last changes: parallelization for laminar and turbulent slot flames
 *
 * \todo labels:LS most cases will not be suited for multiple level sets! check definition of bndryCell in gCell
 * context!parallel implementing?!
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::updateBndryCellList() {
  TRACE();

  for(MInt set = 0; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    switch(m_levelSetBoundaryCondition) {
      case 19940404: {
        MInt noGBndryCellsBackup = a_noGBndryCells(set);

        // inflow boundary
        MInt nghbrXL = -1, nghbrXR = -1, nghbrXL2 = -1, nghbrXR2 = -1, parent, cellId;
        MInt nghbrYL = -1, nghbrYR = -1, nghbrYL2 = -1, nghbrYR2 = -1;

        // reset g boundary information for level set domain boundary cells
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          a_isBndryCellG(a_bandCellId(id, set)) = false;
        }

        std::vector<MInt>().swap(m_gBndryCells[set]);

        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);
          if(a_isHalo(cellId)) continue;

          nghbrXL = c_neighborId(cellId, 0);
          nghbrXL2 = c_neighborId(nghbrXL, 0);
          nghbrXR = c_neighborId(cellId, 1);
          nghbrXR2 = c_neighborId(nghbrXR, 1);

          nghbrYL = c_neighborId(cellId, 2);
          nghbrYL2 = c_neighborId(nghbrYL, 2);
          nghbrYR = c_neighborId(cellId, 3);
          nghbrYR2 = c_neighborId(nghbrYR, 3);

          for(MInt i = 0; i < nDim; i++) {
            parent = c_parentId(cellId);

            if(a_hasNeighbor(cellId, 2 + i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
            }

            if(a_hasNeighbor(nghbrXL, 2 + i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrXL))
                if(a_inBandG(nghbrXL, set) && a_level(nghbrXL) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXL)) {
                  m_gBndryCells[set].push_back(nghbrXL);
                  a_isGBoundaryCellG(nghbrXL, set) = true;
                  a_isBndryCellG(nghbrXL) = true;
                }
            }

            if(a_hasNeighbor(nghbrXL2, 2 + i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrXL2))
                if(a_inBandG(nghbrXL2, set) && a_level(nghbrXL2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXL2)) {
                  m_gBndryCells[set].push_back(nghbrXL2);
                  a_isGBoundaryCellG(nghbrXL2, set) = true;
                  a_isBndryCellG(nghbrXL2) = true;
                }
            }

            if(a_hasNeighbor(nghbrXR, 2 + i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrXR))
                if(a_inBandG(nghbrXR, set) && a_level(nghbrXR) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXR)) {
                  m_gBndryCells[set].push_back(nghbrXR);
                  a_isGBoundaryCellG(nghbrXR, set) = true;
                  a_isBndryCellG(nghbrXR) = true;
                }
            }

            if(a_hasNeighbor(nghbrXR2, 2 + i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrXR2))
                if(a_inBandG(nghbrXR2, set) && a_level(nghbrXR2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXR2)) {
                  m_gBndryCells[set].push_back(nghbrXR2);
                  a_isGBoundaryCellG(nghbrXR2, set) = true;
                  a_isBndryCellG(nghbrXR2) = true;
                }
            }

            if(a_hasNeighbor(nghbrYL, i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrYL))
                if(a_inBandG(nghbrYL, set) && a_level(nghbrYL) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYL)) {
                  m_gBndryCells[set].push_back(nghbrYL);
                  a_isGBoundaryCellG(nghbrYL, set) = true;
                  a_isBndryCellG(nghbrYL) = true;
                }
            }

            if(a_hasNeighbor(nghbrYL2, i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrYL2))
                if(a_inBandG(nghbrYL2, set) && a_level(nghbrYL2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYL2)) {
                  m_gBndryCells[set].push_back(nghbrYL2);
                  a_isGBoundaryCellG(nghbrYL2, set) = true;
                  a_isBndryCellG(nghbrYL2) = true;
                }
            }

            if(a_hasNeighbor(nghbrYR, i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrYR))
                if(a_inBandG(nghbrYR, set) && a_level(nghbrYR) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYR)) {
                  m_gBndryCells[set].push_back(nghbrYR);
                  a_isGBoundaryCellG(nghbrYR, set) = true;
                  a_isBndryCellG(nghbrYR) = true;
                }
            }

            if(a_hasNeighbor(nghbrYR2, i) == 0) {
              if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                m_gBndryCells[set].push_back(cellId);
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
              }
              if(!a_isHalo(nghbrYR2))
                if(a_inBandG(nghbrYR2, set) && a_level(nghbrYR2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYR2)) {
                  m_gBndryCells[set].push_back(nghbrYR2);
                  a_isGBoundaryCellG(nghbrYR2, set) = true;
                  a_isBndryCellG(nghbrYR2) = true;
                }
            }

            if(parent > -1) {
              if(a_hasNeighbor(parent, 2 * i) == 0) {
                for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
                  if(c_childId(parent, ch) < 0) continue;
                  if(a_isHalo(c_childId(parent, ch))) continue;
                  if(a_isBndryCellG(c_childId(parent, ch))) continue;
                  if(!a_inBandG(c_childId(parent, ch), set)) continue;
                  if(a_level(c_childId(parent, ch)) != a_maxGCellLevel()) continue;
                  m_gBndryCells[set].push_back(c_childId(parent, ch));
                  a_isGBoundaryCellG(c_childId(parent, ch), set) = true;
                  a_isBndryCellG(c_childId(parent, ch)) = true;
                }
              }
            }

            if(parent > -1)
              if(a_hasNeighbor(parent, 2 * i + 1) == 0) {
                for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
                  if(c_childId(parent, ch) < 0) continue;
                  if(a_isHalo(c_childId(parent, ch))) continue;
                  if(a_isBndryCellG(c_childId(parent, ch))) continue;
                  if(!a_inBandG(c_childId(parent, ch), set)) continue;
                  if(a_level(c_childId(parent, ch)) != a_maxGCellLevel()) continue;
                  m_gBndryCells[set].push_back(c_childId(parent, ch));
                  a_isGBoundaryCellG(c_childId(parent, ch), set) = true;
                  a_isBndryCellG(c_childId(parent, ch)) = true;
                }
              }
          }
        }

        MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
        MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");
        // gather:
        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          sendBufferSize.p[i] = 0;
          for(MInt j = 0; j < noWindowCells(i); j++) {
            m_intSendBuffers[i][sendBufferSize.p[i]++] = a_isBndryCellG(windowCellId(i, j));
          }
        }
        if(grid().azimuthalPeriodicity()) {
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
              m_intSendBuffers[i][sendBufferSize.p[i]++] = a_isBndryCellG(grid().azimuthalWindowCell(i, j));
            }
          }
        }

        // exchange data -> send, receive
        exchangeIntBuffers(sendBufferSize.getPointer(), receiveBufferSize.getPointer(), 9, 1);

        // scatter:
        // update the window cell list with the halo cells

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
          // MInt windowBndryCnt = 0;
          // check, if received data size matches expected size:
          if(receiveBufferSize.p[i] != noHaloCells(i))
            mTerm(1, AT_, "this was not expected to happen: wrong number of window information...");
#endif
          for(MInt j = 0; j < noHaloCells(i); j++) {
            MInt haloCell = haloCellId(i, j);
            // check if window cell is a boundary cell, otherwise it'll be overwritten from the neighbor domains
            if(m_intReceiveBuffers[i][j]) { //! a_isBndryCellG( haloCell ) && !&a_isGBoundaryCellG( haloCell , 0))
              a_isBndryCellG(haloCell) = m_intReceiveBuffers[i][j];
              a_isGBoundaryCellG(haloCell, set) = m_intReceiveBuffers[i][j];
              if(a_isBndryCellG(haloCell)) {
                m_gBndryCells[set].push_back(haloCell);
              }
            }
          }
        }
        if(grid().azimuthalPeriodicity()) {
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            MInt offset = noHaloCells(grid().azimuthalNeighborDomain(i));
            for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
              MInt n = offset + j;
              MInt haloCell = grid().azimuthalHaloCell(i, j);
              if(m_intReceiveBuffers[i][n]) {
                a_isBndryCellG(haloCell) = m_intReceiveBuffers[i][n];
                a_isGBoundaryCellG(haloCell, set) = m_intReceiveBuffers[i][n];
                if(a_isBndryCellG(haloCell)) {
                  m_gBndryCells[set].push_back(haloCell);
                }
              }
            }
          }
        }
        if(noGBndryCellsBackup != a_noGBndryCells(set) || globalTimeStep == 0) {
          m_log << "Number of level set boundary cells changed: " << endl;
          m_log << a_noGBndryCells(set) << " G boundary cells created..." << endl;
        }
        break;
      }
      case 17516: {
        MInt noGBndryCellsBackup = a_noGBndryCells(set);
        MFloat deltaX = m_flameRadiusOffset;

        if(m_twoFlames) {
          std::vector<MInt>().swap(m_gBndryCells[set]);
          for(MInt id = 0; id < a_noBandCells(set); id++) {
            if(c_coordinate(a_bandCellId(id, set), 1) < F0 || c_coordinate(a_bandCellId(id, set), 1) > m_gCellDistance)
              continue;
            m_gBndryCells[set].push_back(a_bandCellId(id, set));
          }

        } else {
          // inflow boundary
          MInt nghbrXL = -1, nghbrXR = -1, nghbrXL2 = -1, nghbrXR2 = -1, parent, cellId;
          MInt nghbrYL = -1, nghbrYR = -1, nghbrYL2 = -1, nghbrYR2 = -1;

          // reset g boundary information for level set domain boundary cells
          for(MInt id = 0; id < a_noBandCells(set); id++) {
            a_isBndryCellG(a_bandCellId(id, set)) = false;
          }

          std::vector<MInt>().swap(m_gBndryCells[set]);

          for(MInt id = 0; id < a_noBandCells(set); id++) {
            cellId = a_bandCellId(id, set);
            if(a_isHalo(cellId)) continue;

            nghbrXL = c_neighborId(cellId, 0);
            nghbrXL2 = c_neighborId(nghbrXL, 0);
            nghbrXR = c_neighborId(cellId, 1);
            nghbrXR2 = c_neighborId(nghbrXR, 1);

            nghbrYL = c_neighborId(cellId, 2);
            nghbrYL2 = c_neighborId(nghbrYL, 2);
            nghbrYR = c_neighborId(cellId, 3);
            nghbrYR2 = c_neighborId(nghbrYR, 3);

            if(c_coordinate(cellId, 1) < (m_yOffsetFlameTube - 0.1)) continue;

            if(c_coordinate(cellId, 1) > (m_yOffsetFlameTube + 0.1)) continue;

            if(ABS(c_coordinate(cellId, 0)) > (m_radiusFlameTube + deltaX + 0.1)) continue;

            for(MInt i = 0; i < nDim; i++) {
              parent = c_parentId(cellId);

              if(a_hasNeighbor(cellId, 2 + i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
              }

              if(a_hasNeighbor(nghbrXL, 2 + i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrXL))
                  if(a_inBandG(nghbrXL, set) && a_level(nghbrXL) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXL)) {
                    //              if(a_inBandG(nghbrXL,  set ) && !a_isBndryCellG(nghbrXL))
                    m_gBndryCells[set].push_back(nghbrXL);
                    a_isGBoundaryCellG(nghbrXL, set) = true;
                    a_isBndryCellG(nghbrXL) = true;
                  }
              }

              if(a_hasNeighbor(nghbrXL2, 2 + i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrXL2))
                  if(a_inBandG(nghbrXL2, set) && a_level(nghbrXL2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXL2)) {
                    // if(a_inBandG(nghbrXL2,  set ) && !a_isBndryCellG(nghbrXL2))
                    m_gBndryCells[set].push_back(nghbrXL2);
                    a_isGBoundaryCellG(nghbrXL2, set) = true;
                    a_isBndryCellG(nghbrXL2) = true;
                  }
              }

              if(a_hasNeighbor(nghbrXR, 2 + i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrXR))
                  if(a_inBandG(nghbrXR, set) && a_level(nghbrXR) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXR)) {
                    // if(a_inBandG(nghbrXR,  set ) && !a_isBndryCellG(nghbrXR))
                    m_gBndryCells[set].push_back(nghbrXR);
                    a_isGBoundaryCellG(nghbrXR, set) = true;
                    a_isBndryCellG(nghbrXR) = true;
                  }
              }

              if(a_hasNeighbor(nghbrXR2, 2 + i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrXR2))
                  if(a_inBandG(nghbrXR2, set) && a_level(nghbrXR2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrXR2)) {
                    // if(a_inBandG(nghbrXR2,  set ) && !a_isBndryCellG(nghbrXR2))
                    m_gBndryCells[set].push_back(nghbrXR2);
                    a_isGBoundaryCellG(nghbrXR2, set) = true;
                    a_isBndryCellG(nghbrXR2) = true;
                  }
              }

              if(a_hasNeighbor(nghbrYL, i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrYL))
                  if(a_inBandG(nghbrYL, set) && a_level(nghbrYL) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYL)) {
                    // if(a_inBandG(nghbrYL,  set ) && !a_isBndryCellG(nghbrYL))
                    m_gBndryCells[set].push_back(nghbrYL);
                    a_isGBoundaryCellG(nghbrYL, set) = true;
                    a_isBndryCellG(nghbrYL) = true;
                  }
              }

              if(a_hasNeighbor(nghbrYL2, i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrYL2))
                  if(a_inBandG(nghbrYL2, set) && a_level(nghbrYL2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYL2)) {
                    // if(a_inBandG(nghbrYL2,  set ) && !a_isBndryCellG(nghbrYL2))
                    m_gBndryCells[set].push_back(nghbrYL2);
                    a_isGBoundaryCellG(nghbrYL2, set) = true;
                    a_isBndryCellG(nghbrYL2) = true;
                  }
              }

              if(a_hasNeighbor(nghbrYR, i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  // if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrYR))
                  if(a_inBandG(nghbrYR, set) && a_level(nghbrYR) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYR)) {
                    //              if(a_inBandG(nghbrYR,  set ) && !a_isBndryCellG(nghbrYR))
                    m_gBndryCells[set].push_back(nghbrYR);
                    a_isGBoundaryCellG(nghbrYR, set) = true;
                    a_isBndryCellG(nghbrYR) = true;
                  }
              }

              if(a_hasNeighbor(nghbrYR2, i) == 0) {
                if(a_level(cellId) == a_maxGCellLevel() && !a_isBndryCellG(cellId)) {
                  //              if(!a_isBndryCellG(cellId))
                  m_gBndryCells[set].push_back(cellId);
                  a_isGBoundaryCellG(cellId, set) = true;
                  a_isBndryCellG(cellId) = true;
                }
                if(!a_isHalo(nghbrYR2))
                  if(a_inBandG(nghbrYR2, set) && a_level(nghbrYR2) == a_maxGCellLevel() && !a_isBndryCellG(nghbrYR2)) {
                    //              if(a_inBandG(nghbrYR2,  set ) && !a_isBndryCellG(nghbrYR2))
                    m_gBndryCells[set].push_back(nghbrYR2);
                    a_isGBoundaryCellG(nghbrYR2, set) = true;
                    a_isBndryCellG(nghbrYR2) = true;
                  }
              }

              if(parent > -1) {
                if(a_hasNeighbor(parent, 2 * i) == 0) {
                  for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
                    if(c_childId(parent, ch) < 0) continue;
                    if(a_isHalo(c_childId(parent, ch))) continue;
                    if(a_isBndryCellG(c_childId(parent, ch))) continue;
                    if(!a_inBandG(c_childId(parent, ch), set)) continue;
                    if(a_level(c_childId(parent, ch)) != a_maxGCellLevel()) continue;
                    m_gBndryCells[set].push_back(c_childId(parent, ch));
                    a_isGBoundaryCellG(c_childId(parent, ch), set) = true;
                    a_isBndryCellG(c_childId(parent, ch)) = true;
                  }
                }
              }

              // parent = c_parentId(parent);
              //
              // parent = c_parentId( cellId );
              if(parent > -1)
                if(a_hasNeighbor(parent, 2 * i + 1) == 0) {
                  for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
                    if(c_childId(parent, ch) < 0) continue;
                    if(a_isHalo(c_childId(parent, ch))) continue;
                    if(a_isBndryCellG(c_childId(parent, ch))) continue;
                    if(!a_inBandG(c_childId(parent, ch), set)) continue;
                    if(a_level(c_childId(parent, ch)) != a_maxGCellLevel()) continue;
                    m_gBndryCells[set].push_back(c_childId(parent, ch));
                    a_isGBoundaryCellG(c_childId(parent, ch), set) = true;
                    a_isBndryCellG(c_childId(parent, ch)) = true;
                  }
                }
              // parent = c_parentId(parent);
              //
            }
          }
        }

        MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_isBndryCellG(windowCellId(i, j));
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            // check if window cell is a boundary cell, otherwise it'll be overwritten from the neighbor domains
            if(tmp_data(haloCell)) {
              //! a_isBndryCellG( haloCell ) && !&a_isGBoundaryCellG( haloCell , 0)){
              a_isBndryCellG(haloCell) = tmp_data(haloCell);
              a_isGBoundaryCellG(haloCell, set) = tmp_data(haloCell);
              if(a_isBndryCellG(haloCell)) {
                m_gBndryCells[set].push_back(haloCell);
              }
            }
          }
        }

        if(noGBndryCellsBackup != a_noGBndryCells(set) || globalTimeStep == 0) {
          m_log << "Number of level set boundary cells changed: " << endl;
          m_log << a_noGBndryCells(set) << " G boundary cells created..." << endl;
        }
        break;
      }

      case 1751600: {
        MInt noGBndryCellsBackup = a_noGBndryCells(set);
        MFloat deltaX = m_flameRadiusOffset;
        MInt nghbrId, nghbrIdT;
        MInt cnt;

        // reset g boundary information for level set domain boundary cells
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          a_isBndryCellG(a_bandCellId(id, set)) = false;
        }

        std::vector<MInt>().swap(m_gBndryCells[set]);

        for(MInt id = 0; id < a_noBandCells(set); id++) {
          MInt cellId = a_bandCellId(id, set);

          if(a_isHalo(cellId)) continue;

          if(a_level(cellId) != a_maxGCellLevel()) continue;


          IF_CONSTEXPR(nDim == 3) {
            if(c_coordinate(cellId, 1) < (m_yOffsetFlameTube - 0.05)) continue;

            if(c_coordinate(cellId, 1) > (m_yOffsetFlameTube + 0.1)) continue;

            if(ABS(c_coordinate(cellId, 0)) > (m_jetHalfWidth + 0.3)) continue;

            if(ABS(c_coordinate(cellId, 0)) < (m_jetHalfWidth - 0.1)
               && ABS(c_coordinate(cellId, 2)) < (m_jetHalfLength - 0.10))
              continue;

            if(ABS(c_coordinate(cellId, 2)) > (m_jetHalfLength + 0.3)) continue;
          }
          else {
            if(c_coordinate(cellId, 1) < (m_yOffsetFlameTube - 0.1)) continue;

            if(c_coordinate(cellId, 1) > (m_yOffsetFlameTube + 0.1)) continue;

            if(ABS(c_coordinate(cellId, 0)) > (m_radiusFlameTube + deltaX + 0.1)) continue;
          }
          for(MInt dirId = 0; dirId < nDim * 2; dirId++) {
            if(a_hasNeighbor(cellId, dirId) == 0 && a_level(cellId) == a_maxGCellLevel()) {
              a_isGBoundaryCellG(cellId, set) = true;
              a_isBndryCellG(cellId) = true;
              m_gBndryCells[set].push_back(cellId);
            }
            for(MInt dirIdS = 0; dirIdS < nDim * 2; dirIdS++) {
              if(dirId == dirIdS) continue;
              cnt = 0;
              nghbrId = cellId;
              while(cnt < 0) {
                MInt noNghbrIds = a_hasNeighbor(nghbrId, dirIdS);
                if(noNghbrIds == 0) break;

                nghbrId = c_neighborId(nghbrId, dirIdS);
                if(nghbrId == -1) break;

                if(a_isHalo(nghbrId)) break;
                if(!a_inBandG(nghbrId, set)) break;
                if(a_level(nghbrId) != a_maxGCellLevel()) break;
                if(!a_isBndryCellG(nghbrId)) {
                  a_isGBoundaryCellG(nghbrId, set) = true;
                  a_isBndryCellG(nghbrId) = true;
                  m_gBndryCells[set].push_back(nghbrId);
                }
                for(MInt dirIdT = 0; dirIdT < nDim * 2; dirIdT++) {
                  MInt noNghbrIdsT = a_hasNeighbor(nghbrId, dirIdT);
                  if(noNghbrIdsT == 0) continue;

                  nghbrIdT = c_neighborId(nghbrId, dirIdT);
                  if(nghbrIdT == -1) continue;

                  if(a_isHalo(nghbrIdT)) continue;

                  if(!a_inBandG(nghbrIdT, set)) continue;
                  if(!a_isBndryCellG(nghbrIdT)) {
                    a_isGBoundaryCellG(nghbrIdT, set) = true;
                    a_isBndryCellG(nghbrIdT) = true;
                    m_gBndryCells[set].push_back(nghbrIdT);
                  }
                }


                if(noDomains() > 1 && cnt == m_noHaloLayers + 1) {
                  stringstream errorMessage;
                  errorMessage << "cnt is equal to " << cnt << " exiting ....";
                  mTerm(1, AT_, errorMessage.str());
                }
                cnt++;
              }
            }
          }
        }

        MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_isBndryCellG(windowCellId(i, j));
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            if(tmp_data(haloCell)) {
              a_isBndryCellG(haloCell) = tmp_data(haloCell);
              a_isGBoundaryCellG(haloCell, set) = tmp_data(haloCell);
              if(a_isBndryCellG(haloCell)) {
                m_gBndryCells[set].push_back(haloCell);
              }
            }
          }
        }
        if(noGBndryCellsBackup != a_noGBndryCells(set) || globalTimeStep == 0) {
          m_log << "Number of level set boundary cells changed: " << endl;
          m_log << a_noGBndryCells(set) << " G boundary cells created..." << endl;
        }
        break;
      }

      // massive parallel
      case 5401000: {
        MInt noGBndryCellsBackup = a_noGBndryCells(set);
        // MFloat deltaX = m_flameRadiusOffset;
        MInt nghbrId, nghbrIdT;
        MInt cnt;
        MFloat radius;
        // reset g boundary information for level set domain boundary cells
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          a_isBndryCellG(a_bandCellId(id, set)) = false;
        }
        std::vector<MInt>().swap(m_gBndryCells[set]);
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          MInt cellId = a_bandCellId(id, set);

          radius = sqrt(POW2(c_coordinate(cellId, 0)) + POW2(c_coordinate(cellId, 2)));


          if(a_isHalo(cellId)) continue;
          if(a_level(cellId) != a_maxGCellLevel()) continue;
          IF_CONSTEXPR(nDim == 3) {
            if(c_coordinate(cellId, 1) < (m_yOffsetFlameTube - 0.05)) continue;
            if(c_coordinate(cellId, 1) > (m_yOffsetFlameTube + 0.1)) continue;
            if(radius > (0.55)) continue;
            if(radius < 0.48) continue;
          }
          else {
            // if(c_coordinate(cellId, 1)<(m_yOffsetFlameTube-0.1))
            // continue;
            // if(c_coordinate(cellId, 1)>(m_yOffsetFlameTube+0.1))
            // continue;
            // if(ABS(c_coordinate(cellId, 0))> (m_radiusFlameTube+deltaX+0.1))
            // continue;
          }
          for(MInt dirId = 0; dirId < nDim * 2; dirId++) {
            if(a_hasNeighbor(cellId, dirId) == 0 && a_level(cellId) == a_maxGCellLevel()) {
              if(!a_isBndryCellG(cellId)) {
                a_isGBoundaryCellG(cellId, set) = true;
                a_isBndryCellG(cellId) = true;
                m_gBndryCells[set].push_back(cellId);
              }

              // go in all other directions
              for(MInt dirIdS = 0; dirIdS < nDim * 2; dirIdS++) {
                if(dirId == dirIdS) continue;
                cnt = 0;
                nghbrId = cellId;
                while(cnt < 0) { // Stephan/Jerry
                  //              nghbrId= c_neighborId(nghbrId, dirIdS);
                  //              if(nghbrId==-1)
                  //                break;

                  MInt noNghbrIds = a_hasNeighbor(nghbrId, dirIdS);
                  if(noNghbrIds == 0) // use noNghbrs because of determineStructured...() which
                                      // sets the neighbor to the current cell if there is no
                                      // neighbor
                    break;

                  nghbrId = c_neighborId(nghbrId, dirIdS);
                  if(nghbrId == -1) // just to be sure
                    break;
                  // don't add halo cells to the list, the info will be send to the other domains
                  // and they update the window cells
                  if(a_isHalo(nghbrId)) break;
                  if(!a_inBandG(nghbrId, set)) break;
                  if(a_level(nghbrId) != a_maxGCellLevel()) break;
                  if(!a_isBndryCellG(nghbrId)) {
                    a_isGBoundaryCellG(nghbrId, set) = true;
                    a_isBndryCellG(nghbrId) = true;
                    m_gBndryCells[set].push_back(nghbrId);
                  }
                  // go in all other directions
                  for(MInt dirIdT = 0; dirIdT < nDim * 2; dirIdT++) {
                    //                if(dirIdT==dirId)
                    //    continue;
                    //  if(dirIdT==dirIdS)
                    //  continue;

                    MInt noNghbrIdsT = a_hasNeighbor(nghbrId, dirIdT);
                    if(noNghbrIdsT == 0) // use noNghbrs because of determineStructured...() which
                                         // sets the neighbor to the current cell if there is no
                                         // neighbor
                      continue;

                    nghbrIdT = c_neighborId(nghbrId, dirIdT);
                    if(nghbrIdT == -1) // just to be sure
                      continue;
                    if(a_isHalo(nghbrIdT)) continue;

                    if(!a_inBandG(nghbrIdT, set)) continue;
                    if(!a_isBndryCellG(nghbrIdT)) {
                      a_isGBoundaryCellG(nghbrIdT, set) = true;
                      a_isBndryCellG(nghbrIdT) = true;
                      m_gBndryCells[set].push_back(nghbrIdT);
                    }
                  }


                  if(noDomains() > 1 && cnt == m_noHaloLayers + 1) {
                    stringstream errorMessage;
                    errorMessage << "cnt is equal to " << cnt << " exiting ....";
                    mTerm(1, AT_, errorMessage.str());
                  }
                  cnt++;
                }
              }
            }
          }
        }
        // exchange halo boundary cells for all sets with other domains
        MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_isBndryCellG(windowCellId(i, j));
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            if(tmp_data(haloCell)) { //! a_isBndryCellG( haloCell ) &&
                                     //! !&a_isGBoundaryCellG( haloCell , 0)){
              a_isBndryCellG(haloCell) = tmp_data(haloCell);
              a_isGBoundaryCellG(haloCell, set) = tmp_data(haloCell);
              if(a_isBndryCellG(haloCell)) {
                m_gBndryCells[set].push_back(haloCell);
              }
            }
          }
        }
        if(noGBndryCellsBackup != a_noGBndryCells(set) || globalTimeStep == 0) {
          m_log << "Number of level set boundary cells changed: " << endl;
          m_log << a_noGBndryCells(set) << " G boundary cells created...xy" << endl;
        }
        break;
      }
      default: {
        break; // most testcases do not need this
      }
    }
  }
}


// ----------------------------------------------------------------------------------------


/** \brief resets level set cells outside the band!
 *         Update: reset only leaf-cells and limit the levelset value for lower-grid-levels!
 *
 * \author Daniel Hartmann, Stephan Schlimpert, Claudia Guenther, Tim Wegmann
 * \date unknown, September 2011, November 2012, October 2018
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::resetOutsideCells(MInt computingSet) {
  TRACE();

  if(m_levelSetRans) return;
  if(m_virtualSurgery) return;

  MInt startSet = 0;
  MInt endSet = m_noSets;
  if(m_buildCollectedLevelSetFunction && globalTimeStep > 0) startSet = 1;
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
    ASSERT(m_computeSet[computingSet], "");
    m_changedSet[computingSet] = true;
    ASSERT(m_changedSet[computingSet], "");
  }


  //---
  // set level set function outside of the band
  for(MInt set = startSet; set < endSet; set++) {
    if(!m_computeSet[set]) continue;
    if(!m_changedSet[set] && m_levelSetMb && a_maxGCellLevel(set) == maxRefinementLevel()) continue;

    for(MInt cell = 0; cell < a_noCells(); cell++) {
      if(a_inBandG(cell, set)) continue;
      if(a_level(cell) > a_maxGCellLevel(0)) continue;

      if(c_noChildren(cell) > 0 && m_levelSetMb) { // on lower levels:
        if(a_levelSetFunctionG(cell, set) > m_outsideGValue) {
          a_levelSetFunctionG(cell, set) = m_outsideGValue;
        } else if(a_levelSetFunctionG(cell, set) < -m_outsideGValue) {
          a_levelSetFunctionG(cell, set) = -m_outsideGValue;
        }
      } else { // leaf-cells which are not in the band!

        if(a_levelSetFunctionG(cell, set) > F0) {
          a_levelSetFunctionG(cell, set) = m_outsideGValue;
        } else {
          a_levelSetFunctionG(cell, set) = -m_outsideGValue;
        }

        if(!m_semiLagrange) {
          a_correctedBurningVelocity(cell, set) = F0;
          for(MInt i = 0; i < nDim; i++) {
            a_extensionVelocityG(cell, i, set) = F0;
            a_normalVectorG(cell, i, set) = F0;
          }
          a_curvatureG(cell, set) = F0;
        }
      }
    }
  }
}

template <MInt nDim>
void LsCartesianSolver<nDim>::resetOldOutsideCells() {
  TRACE();

  // set level set function outside of the band
  for(MInt set = 0; set < m_noSets; set++) {
    for(MInt cell = 0; cell < noInternalCells(); cell++) {
      if(a_oldLevelSetFunctionG(cell, set) > m_outsideGValue) {
        a_oldLevelSetFunctionG(cell, set) = m_outsideGValue;
      } else if(a_oldLevelSetFunctionG(cell, set) < -m_outsideGValue) {
        a_oldLevelSetFunctionG(cell, set) = -m_outsideGValue;
      }
    }
  }
}
// ----------------------------------------------------------------------------------------


/** \brief computes the curvature
 *
 * \author Daniel Hartmann, Stephan Schlimpert,
 * \date 2007, May 2011

 * computes the curvature with different order of discretization, the same order is also used
 * for the computation of the normal
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeCurvature(MInt computingSet) {
  TRACE();

  ASSERT(m_lsCollectorMode == 0 || m_lsCollectorMode == 2 || m_lsCollectorMode == 3, "");


  IF_CONSTEXPR(nDim == 2) {
    MInt cellId, nghbrL, nghbrL2, nghbrR, nghbrR2;
    // MBool boundary= false;
    MBool oneSidedR = false;
    MBool oneSidedL = false;
    MFloat FgCellDistanceL2R2 = m_FgCellDistance * F1B4, FgCellDistanceLR = m_FgCellDistance * F1B2;
    MFloat levelSetNegative = F0, levelSetPlus = F0;
    MFloat filterCoord = m_filterFlameTubeEdgesDistance + 0.0234;
    MFloat cR = -9999.9, cL = -9999.9; //,cR2=-9999.9,cL2=-9999.9;
    MBool filter = false;

    MInt startSet = 0;
    MInt endSet = m_noSets;
    if(computingSet >= 0) {
      startSet = computingSet;
      endSet = computingSet + 1;
    }
    //---

#ifdef standardStencil

    for(MInt set = startSet; set < endSet; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);

        if(a_isHalo(cellId)) continue;

        // compute the second derivatives
        for(MInt i = 0; i < nDim; i++) {
          grad2G[i] =
              (a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i + 1], set)
               + a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i], set) - F2 * a_levelSetFunctionG(cellId, set))
              * POW2(Fdx);
        }

        // compute the mixed derivatives

        // dummy:
        mixedDerivative = 0.01;
        mTerm(1, AT_, "ERROR: check code.. this part is assumed to not be used ");

        // compute the denominator
        denominator = F0;
        for(MInt i = 0; i < nDim; i++) {
          denominator += POW2(a_levelSetFunctionSlope(cellId, i, set));
        }
        denominator = F1 / POW2(denominator);
        denominator = sqrt(denominator * denominator * denominator);

        a_curvatureG(cellId, set) = -denominator
                                    * (grad2G[0] * POW2(a_levelSetFunctionSlope(cellId, 1, set))
                                       + grad2G[1] * POW2(a_levelSetFunctionSlope(cellId, 0, set))
                                       - F2 * a_levelSetFunctionSlope(cellId, 0, set)
                                             * a_levelSetFunctionSlope(cellId, 1, set) * mixedDerivative);
      }
    }
#else

    if(m_fourthOrderNormalCurvatureComputation) {
      for(MInt set = startSet; set < endSet; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);

          // reset curvature
          a_curvatureG(cellId, set) = F0;
          if(a_isHalo(cellId)) continue;

          for(MInt i = 0; i < nDim; i++) {
            oneSidedR = false;
            oneSidedL = false;
            filter = false;

            FgCellDistanceL2R2 = m_FgCellDistance * F1B4;
            FgCellDistanceLR = m_FgCellDistance * F1B2;

            nghbrL = c_neighborId(cellId, 2 * i);
            nghbrL2 = c_neighborId(nghbrL, 2 * i);

            nghbrR = c_neighborId(cellId, 2 * i + 1);
            nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);

            // reduce to second order
            if(!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR2, set)) {
              nghbrR2 = nghbrR;
              nghbrL2 = nghbrL;
              FgCellDistanceL2R2 = FgCellDistanceLR;
            }

            // reduce to second order on boundaries (like for the Landau simulations)
            if(a_hasNeighbor(nghbrL, 2 * i) == 0 || a_hasNeighbor(nghbrR, 2 * i + 1) == 0) {
              // reduce to second order
              nghbrR2 = nghbrR;
              nghbrL2 = nghbrL;
              FgCellDistanceL2R2 = FgCellDistanceLR;
            }

            // take onesided differences directly on the boundary
            // no neighbors in positive direction
            if(a_hasNeighbor(cellId, 2 * i + 1) == 0 || a_hasNeighbor(nghbrR, 2 * i + 1) == 0
               || a_hasNeighbor(nghbrR2, 2 * i + 1) == 0) {
              oneSidedL = true;

              a_curvatureG(cellId, set) = F0;
              /*
                - m_FgCellDistance *

                ( a * a_normalVectorG( cellId ,  i, set) ) +
                b * a_normalVectorG( nghbrL ,  i, set ) +
                c * a_normalVectorG( nghbrL2 ,  i, set ) );
              */
            }

            // no neighbors in negative direction
            if(a_hasNeighbor(cellId, 2 * i) == 0 || a_hasNeighbor(nghbrL, 2 * i) == 0
               || a_hasNeighbor(nghbrL2, 2 * i) == 0) {
              oneSidedR = true;

              a_curvatureG(cellId, set) = F0;
              /*
                m_FgCellDistance *

                ( a * a_normalVectorG( cellId ,  i, set ) +
                b * a_normalVectorG( nghbrR ,  i, set ) +
                c * a_normalVectorG( nghbrR2 ,  i, set ) );
              */
            }

            // filtering
            if(!oneSidedR && !oneSidedL && m_filterFlameTubeEdges) {
              if(c_coordinate(cellId, 1) < filterCoord) {
                // compute second order curvature values at surfaces (staggered points)

                cR = (a_normalVectorG(nghbrR, i, set) - a_normalVectorG(cellId, i, set)) * m_FgCellDistance;

                cL = (a_normalVectorG(cellId, i, set) - a_normalVectorG(nghbrL, i, set)) * m_FgCellDistance;

                // compute via second order interpolation the curvature back to the cell centers
                a_curvatureG(cellId, 0) += (F1B2 * (cR + cL));
                filter = true;
              }
            }

            if(!filter) {
              if(!oneSidedR && !oneSidedL) {
                a_curvatureG(cellId, set) +=

                    F4B3 * FgCellDistanceLR *

                        (a_normalVectorG(nghbrR, i, set) - a_normalVectorG(nghbrL, i, set))

                    - F1B3 * FgCellDistanceL2R2 *

                          (a_normalVectorG(nghbrR2, i, set) - a_normalVectorG(nghbrL2, i, set));
              }
            }
          }
          if(m_curvatureDamp) {
            levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
            levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

            a_curvatureG(cellId, set) = F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0))
                                        * a_curvatureG(cellId, set);
          }
        }
      }
    } else {
      for(MInt set = startSet; set < endSet; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);

          a_curvatureG(cellId, set) = F0;
          if(a_isHalo(cellId)) continue;

          if(!a_isGBoundaryCellG(cellId, set)) {
            // second order old code
            a_curvatureG(cellId, set) = F1B2 * m_FgCellDistance
                                        * (a_normalVectorG(a_bandNghbrIdsG(cellId, 1, set), 0, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 0, set), 0, set)
                                           + a_normalVectorG(a_bandNghbrIdsG(cellId, 3, set), 1, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 2, set), 1, set));

            if(m_curvatureDamp) {
              levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
              levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

              a_curvatureG(cellId, set) = F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0))
                                          * a_curvatureG(cellId, set);
            }
          }
        }
#endif
  }
}

// exchange curavture on halo cells, cause they are needed for the Markstein length computation!
exchangeDataLS(&a_curvatureG(0, 0), m_maxNoSets);
  }
  else {
    MInt cellId;
    MInt startSet = 0;
    MInt endSet = m_noSets;
    if(computingSet >= 0) {
      startSet = computingSet;
      endSet = computingSet + 1;
    }
    //---
    for(MInt set = startSet; set < endSet; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);
        if(a_isHalo(cellId)) continue;

        if(!a_isGBoundaryCellG(cellId, set)) {
#ifdef standardStencil
          // compute the second derivatives
          for(MInt i = 0; i < nDim; i++) {
            grad2G[i] =
                (a_levelSetFunctionG(c_neighborId(cellId, 2 * i + 1), set)
                 + a_levelSetFunctionG(c_neighborId(cellId, 2 * i), set) - F2 * a_levelSetFunctionG(cellId, set))
                * POW2(m_FgCellDistance);
          }

          // order: xy,xz,yz
          mixedDerivative[0] = F1B4 * POW2(m_FgCellDistance)
                               * (a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 0), 2 + 0), set)
                                  + a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 1), 2 + 1), set)
                                  - a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 0), 2 + 1), set)
                                  - a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 1), 2 + 0), set));
          mixedDerivative[1] = F1B4 * POW2(m_FgCellDistance)
                               * (a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 0), 4 + 0), set)
                                  + a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 1), 4 + 1), set)
                                  - a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 0), 4 + 1), set)
                                  - a_levelSetFunctionG(c_neighborId(c_neighborId(cellId, 1), 4 + 0), set));
          mixedDerivative[2] =
              F1B4 * POW2(m_FgCellDistance)
              * (a_levelSetFunctionG(c_neighborId(c_neighborId(c_neighborId(cellId, 2), 2 + 0), 4 + 0), set)
                 + a_levelSetFunctionG(c_neighborId(c_neighborId(c_neighborId(cellId, 2), 2 + 1), 4 + 1), set)
                 - a_levelSetFunctionG(c_neighborId(c_neighborId(c_neighborId(cellId, 2), 2 + 0), 4 + 1), set)
                 - a_levelSetFunctionG(c_neighborId(c_neighborId(c_neighborId(cellId, 2), 2 + 1), 4 + 0), set));

          // compute the denominator
          denominator = F0;
          for(MInt i = 0; i < nDim; i++) {
            denominator += POW2(a_levelSetFunctionSlope(cellId, i, set));
          }
          denominator = F1 / POW2(denominator);
          denominator = sqrt(denominator * denominator * denominator);

          a_curvatureG(cellId, set) =
              -denominator
              * (POW2(grad2G[0])
                     * (POW2(a_levelSetFunctionSlope(cellId, 1, set)) + POW2(a_levelSetFunctionSlope(cellId, 2, set)))
                 + POW2(grad2G[1])
                       * (POW2(a_levelSetFunctionSlope(cellId, 0, set)) + POW2(a_levelSetFunctionSlope(cellId, 2, set)))
                 + POW2(grad2G[2])
                       * (POW2(a_levelSetFunctionSlope(cellId, 0, set)) + POW2(a_levelSetFunctionSlope(cellId, 1, set)))
                 - F2
                       * (mixedDerivative[0] * a_levelSetFunctionSlope(cellId, 0, set)
                              * a_levelSetFunctionSlope(cellId, 1, set)
                          + mixedDerivative[1] * a_levelSetFunctionSlope(cellId, 0, set)
                                * a_levelSetFunctionSlope(cellId, 2, set)
                          + mixedDerivative[2] * a_levelSetFunctionSlope(cellId, 1, set)
                                * a_levelSetFunctionSlope(cellId, 2, set)));

#else

          a_curvatureG(cellId, set) = F1B2 * m_FgCellDistance
                                      * (a_normalVectorG(a_bandNghbrIdsG(cellId, 1, set), 0, set)
                                         - a_normalVectorG(a_bandNghbrIdsG(cellId, 0, set), 0, set)
                                         + a_normalVectorG(a_bandNghbrIdsG(cellId, 3, set), 1, set)
                                         - a_normalVectorG(a_bandNghbrIdsG(cellId, 2, set), 1, set)
                                         + a_normalVectorG(a_bandNghbrIdsG(cellId, 5, set), 2, set)
                                         - a_normalVectorG(a_bandNghbrIdsG(cellId, 4, set), 2, set));

          if(m_curvatureDamp) {
            MFloat levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
            MFloat levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

            a_curvatureG(cellId, set) *= F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          }
#endif
        }
      }
    }


    // exchange curvature
    exchangeDataLS(&a_curvatureG(0, 0), m_maxNoSets);
  }
}

/** \brief computes the curvature 2D version for periodic boundary conditions
 *
 * \author Stephan Schlimpert,
 * \date July 2011

 * computes the curvature with different order of discretization, the same order is also used
 * for the computation of the normal
 *
 * \todo labels:LS implement periodic boundaries for all discretizations!!! (only done for
 fourthOrderNormalCurvatureComputation2)
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeCurvaturePeriodic() {
  TRACE();
  IF_CONSTEXPR(nDim == 2) {
    MInt cellId, nghbrL, nghbrL2, nghbrR, nghbrR2;
    // MBool boundary= false;
    MBool oneSidedR = false;
    MBool oneSidedL = false;
    MFloat FgCellDistanceL2R2 = m_FgCellDistance * F1B4, FgCellDistanceLR = m_FgCellDistance * F1B2;
    MFloat levelSetNegative = F0, levelSetPlus = F0;
    //---


#ifdef standardStencil

    for(MInt set = 0; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);
        if(a_isHalo(cellId)) continue;

        // compute the second derivatives
        for(MInt i = 0; i < nDim; i++) {
          grad2G[i] =
              (a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i + 1], set)
               + a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i], set) - F2 * a_levelSetFunctionG(cellId, set))
              * POW2(Fdx);
        }

        // compute the mixed derivatives

        // dummy:
        mixedDerivative = 0.01;
        mTerm(1, AT_, "ERROR: check code.. this part is assumed to not be used ");

        // compute the denominator
        denominator = F0;
        for(MInt i = 0; i < nDim; i++) {
          denominator += POW2(a_levelSetFunctionSlope(cellId, i, set));
        }
        denominator = F1 / POW2(denominator);
        denominator = sqrt(denominator * denominator * denominator);

        a_curvatureG(cellId, set) = -denominator
                                    * (grad2G[0] * POW2(a_levelSetFunctionSlope(cellId, 1, set))
                                       + grad2G[1] * POW2(a_levelSetFunctionSlope(cellId, 0, set))
                                       - F2 * a_levelSetFunctionSlope(cellId, 0, set)
                                             * a_levelSetFunctionSlope(cellId, 1, set) * mixedDerivative);
      }
    }

#else

    if(m_fourthOrderNormalCurvatureComputation) {
      for(MInt set = 0; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);

          // reset curvature
          a_curvatureG(cellId, set) = F0;
          if(a_isHalo(cellId)) continue;

          for(MInt i = 0; i < nDim; i++) {
            oneSidedR = false;
            oneSidedL = false;

            FgCellDistanceL2R2 = m_FgCellDistance * F1B4;
            FgCellDistanceLR = m_FgCellDistance * F1B2;

            nghbrL = c_neighborId(cellId, 2 * i);
            nghbrL2 = c_neighborId(nghbrL, 2 * i);
            nghbrR = c_neighborId(cellId, 2 * i + 1);
            nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);

            if(a_hasNeighbor(nghbrL2, 2 * i) != 0 && a_hasNeighbor(nghbrR2, 2 * i + 1) != 0) {
              // reduce to second order on band boundaries
              if((!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR2, set))) {
                nghbrR2 = nghbrR;
                nghbrL2 = nghbrL;
                FgCellDistanceL2R2 = m_FgCellDistance * F1B2;
              }
            }

            if(!oneSidedR && !oneSidedL) {
              a_curvatureG(cellId, set) +=

                  F4B3 * FgCellDistanceLR *

                      (a_normalVectorG(nghbrR, i, set) - a_normalVectorG(nghbrL, i, set))

                  - F1B3 * FgCellDistanceL2R2 *

                        (a_normalVectorG(nghbrR2, i, set) - a_normalVectorG(nghbrL2, i, set));
            }
          }
          if(m_curvatureDamp) {
            levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
            levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

            a_curvatureG(cellId, set) *= F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          }
        }
      }
    } else {
      // second order code
      for(MInt set = 0; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);
          a_curvatureG(cellId, set) = F0;
          if(a_isHalo(cellId)) continue;

          if(!a_isGBoundaryCellG(cellId, set)) {
            a_curvatureG(cellId, set) = F1B2 * m_FgCellDistance
                                        * (a_normalVectorG(a_bandNghbrIdsG(cellId, 1, set), 0, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 0, set), 0, set)
                                           + a_normalVectorG(a_bandNghbrIdsG(cellId, 3, set), 1, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 2, set), 1, set));

            if(m_curvatureDamp) {
              levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
              levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

              a_curvatureG(cellId, set) *=
                  F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
            }
          }
        }
      }
    }
#endif
  }
  else {
    mTerm(1, AT_, "ERROR: check 3D implementation, not done ");
    MInt cellId, nghbrL, nghbrL2, nghbrR, nghbrR2;
    MBool oneSidedR = false;
    MBool oneSidedL = false;
    MFloat FgCellDistanceL2R2 = m_FgCellDistance * F1B4, FgCellDistanceLR = m_FgCellDistance * F1B2;
    MFloat levelSetNegative = F0, levelSetPlus = F0;
    //---
#ifdef standardStencil
    for(MInt set = 0; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);
        if(a_isHalo(cellId)) continue;
        // compute the second derivatives
        for(MInt i = 0; i < nDim; i++) {
          grad2G[i] =
              (a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i + 1], set)
               + a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i], set) - F2 * a_levelSetFunctionG(cellId, set))
              * POW2(Fdx);
        }
        // compute the mixed derivatives
        for(MInt i = 0; i < nDim; i++) {
          grad2G[i] =
              (a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i + 1], set)
               + a_levelSetFunctionG(m_cells[cellId].m_cells[2 * i], set) - F2 * a_levelSetFunctionG(cellId, set))
              * POW2(Fdx);
        }
        // compute the mixed derivatives
        // dummy:
        mixedDerivative = 0.01;
        mTerm(1, AT_, "ERROR: check code.. this part is assumed to not be used ");

        // compute the denominator
        denominator = F0;
        for(MInt i = 0; i < nDim; i++) {
          denominator += POW2(a_levelSetFunctionSlope(cellId, i, set));
        }
        denominator = F1 / POW2(denominator);
        denominator = sqrt(denominator * denominator * denominator);

        a_curvatureG(cellId, set) = -denominator
                                    * (grad2G[0] * POW2(a_levelSetFunctionSlope(cellId, 1, set))
                                       + grad2G[1] * POW2(a_levelSetFunctionSlope(cellId, 0, set))
                                       - F2 * a_levelSetFunctionSlope(cellId, 0, set)
                                             * a_levelSetFunctionSlope(cellId, 1, set) * mixedDerivative);
      }
    }
#else

    if(m_fourthOrderNormalCurvatureComputation) {
      for(MInt set = 0; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);
          // reset curvature
          a_curvatureG(cellId, set) = F0;
          if(a_isHalo(cellId)) continue;

          for(MInt i = 0; i < nDim; i++) {
            oneSidedR = false;
            oneSidedL = false;

            FgCellDistanceL2R2 = m_FgCellDistance * F1B4;
            FgCellDistanceLR = m_FgCellDistance * F1B2;

            nghbrL = c_neighborId(cellId, 2 * i);
            nghbrL2 = c_neighborId(nghbrL, 2 * i);
            nghbrR = c_neighborId(cellId, 2 * i + 1);
            nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);

            if(a_hasNeighbor(nghbrL2, 2 * i) != 0 && a_hasNeighbor(nghbrR2, 2 * i + 1) != 0) {
              // reduce to second order on band boundaries
              if((!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR2, set))) {
                nghbrR2 = nghbrR;
                nghbrL2 = nghbrL;
                FgCellDistanceL2R2 = m_FgCellDistance * F1B2;
              }
            }

            if(!oneSidedR && !oneSidedL) {
              a_curvatureG(cellId, set) +=

                  F4B3 * FgCellDistanceLR *

                      (a_normalVectorG(nghbrR, i, set) - a_normalVectorG(nghbrL, i, set))

                  - F1B3 * FgCellDistanceL2R2 *

                        (a_normalVectorG(nghbrR2, i, set) - a_normalVectorG(nghbrL2, i, set));
            }
          }
          if(m_curvatureDamp) {
            levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
            levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

            a_curvatureG(cellId, set) = F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0))
                                        * a_curvatureG(cellId, set);
          }
        }
      }
    } else {
      for(MInt set = 0; set < m_noSets; set++) {
        if(!m_computeSet[set]) continue;
        for(MInt id = 0; id < a_noBandCells(set); id++) {
          cellId = a_bandCellId(id, set);
          if(a_isHalo(cellId)) continue;

          if(!a_isGBoundaryCellG(cellId, set)) {
            // second order old code
            a_curvatureG(cellId, set) = F1B2 * m_FgCellDistance
                                        * (a_normalVectorG(a_bandNghbrIdsG(cellId, 1, set), 0, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 0, set), 0, set)
                                           + a_normalVectorG(a_bandNghbrIdsG(cellId, 3, set), 1, set)
                                           - a_normalVectorG(a_bandNghbrIdsG(cellId, 2, set), 1, set));

            if(m_curvatureDamp) {
              levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
              levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);

              a_curvatureG(cellId, set) = F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0))
                                          * a_curvatureG(cellId, set);
            }
          }
        }
      }
    }
#endif
  }
  //#endif
}

//#endif


//----------------------------------------------------------------------------


template <MInt nDim>
void LsCartesianSolver<nDim>::determinePropagationSpeed() {
  TRACE();

  ASSERT(m_lsCollectorMode == 0 || m_lsCollectorMode == 3, "");

  if(m_gRKStep == 0) {
    if(globalTimeStep % 100 == 0) {
      if(domainId() == 0) {
        cerr << "time step #" << globalTimeStep << " - time " << m_time << endl;
      }
    }

    switch(m_levelSetTestCase) {
      case 50431:
      case 504312: {
        for(MInt set = m_startSet; set < m_noSets; set++) {
          for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
            a_extensionVelocityG(cellId, 0, set) = PI / 31.4 * (-c_coordinate(cellId, 1));
            a_extensionVelocityG(cellId, 1, set) = PI / 31.4 * (c_coordinate(cellId, 0));
            if(m_levelSetTestCase == 504312) {
              a_extensionVelocityG(cellId, 0, set) *= -1;
              a_extensionVelocityG(cellId, 1, set) *= -1;
            }
            a_flameSpeedG(cellId, set) = F0;
          }
        }
        break;
      }

      case 203: {
        // Zalesak's problem
        for(MInt set = m_startSet; set < m_noSets; set++) {
          for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
            // flow velocities are not set for this levelSet-testcase!
            // fvSolverD().a_variable(cellId, 0) = fvSolverD().m_rhoInfinity * PI / 31.4 * (-c_coordinate(cellId, 1));
            // fvSolverD().a_variable(cellId, 1) = fvSolverD().m_rhoInfinity * PI / 31.4 * (c_coordinate(cellId, 0));
            a_extensionVelocityG(cellId, 0, set) = PI / 31.4 * (-c_coordinate(cellId, 1));
            a_extensionVelocityG(cellId, 1, set) = PI / 31.4 * (c_coordinate(cellId, 0));
            a_flameSpeedG(cellId, set) = F0;
          }
        }
        break;
      }
      case 3012: {
        /*************** Translation *************/
        MFloat Vtrans[3] = {F0, F0, F0};
        MFloatScratchSpace VtransBodies(m_noEmbeddedBodies, 3, AT_, "VtransBodies");
        for(MInt body = 0; body < m_noEmbeddedBodies; body++) {
          computeBodyPropertiesForced(2, Vtrans, body, time());
          for(MInt i = 0; i < 3; i++) {
            VtransBodies(body, i) = Vtrans[i];
          }
        }
        /***************************************************************************************************************/

        /************** Define Velocity over the G0 cells **************/
        for(MInt set = m_startSet; set < m_noSets; set++) {
          if(!m_computeSet[set]) {
            continue;
          }
          for(MInt id = 0; id < a_noG0Cells(set); id++) {
            MInt cellId = a_G0CellId(id, set);
            MInt body = a_bodyIdG(cellId, set);
            if(body > -1) {
              a_extensionVelocityG(cellId, 0, set) = VtransBodies(body, 0);
              a_extensionVelocityG(cellId, 1, set) = VtransBodies(body, 1);
              IF_CONSTEXPR(nDim == 3) { a_extensionVelocityG(cellId, 2, set) = VtransBodies(body, 2); }
            } else {
              a_extensionVelocityG(cellId, 0, set) = F0;
              a_extensionVelocityG(cellId, 1, set) = F0;
              IF_CONSTEXPR(nDim == 3) { a_extensionVelocityG(cellId, 2, set) = F0; }
            }
          }
        }


        /****************** Updating control points: need new time step **************/
        MFloat dx[3] = {F0, F0, F0};
        MFloat x1[3] = {F0, F0, F0};
        MFloat x0[3] = {F0, F0, F0};
        if(m_GCtrlPntMethod == 2 && m_initialCondition != 31) {
          // Movable STL
          MInt bcId = 0;
          for(MInt body = 0; body < m_noEmbeddedBodies; body++) {
            bcId = m_bodyBndryCndIds[body];
            computeBodyPropertiesForced(1, x1, body, time() + timeStep());
            computeBodyPropertiesForced(1, x0, body, time());
            for(MInt i = 0; i < 3; i++) {
              dx[i] = (x1[i] - x0[i]);
            }
            m_gCtrlPnt.CtrlPnt2_shiftSTL(bcId, dx);
          }
        }
        /********************* Update Level Set time ***************************/
        m_time += timeStep();
        break;
      }
      default:
        break;
    }
  }
}


// ----------------------------------------------------------------------------------------

/** \brief computes the normal in the level set band
 *
 * \author Daniel Hartmann, Stephan Schlimpert
 * \date 2007, May 2011
 *
 * computes the normal with different order of discretization, the same order is also used
 * for the computation of the curvature
 *
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::computeNormalVectors(MInt computingSet) {
  TRACE();

  ASSERT(m_lsCollectorMode == 0 || m_lsCollectorMode == 2 || m_lsCollectorMode == 3, m_lsCollectorMode);

  MFloat FgradG;
  MFloat epsilon = m_gCellDistance * 0.00000001;
  // const MFloat Fdx = F1B2 * m_FgCellDistance;
  MInt cellId, nghbrL, nghbrL2, nghbrR, nghbrR2;
  // MBool boundary= false;
  MBool oneSidedR = false;
  MBool oneSidedL = false;
  MFloat FgCellDistanceL2R2 = m_FgCellDistance * F1B4, FgCellDistanceLR = m_FgCellDistance * F1B2;
  MFloat filterCoord = m_filterFlameTubeEdgesDistance + 0.0234;
  MFloat a = -F3B2, b = F2, c = -F1B2;
  MBool filter = false;
  MFloat levelSlopeR = -9999.9, levelSlopeL = -9999.9;

  MInt startSet = 0;
  MInt endSet = m_noSets;
  if(computingSet >= 0) {
    startSet = computingSet;
    endSet = computingSet + 1;
  }

  exchangeLevelSet();

  //---
  if(m_fourthOrderNormalCurvatureComputation) {
    for(MInt set = startSet; set < endSet; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);

        for(MInt i = 0; i < nDim; i++) {
          oneSidedR = false;
          oneSidedL = false;
          filter = false;

          FgCellDistanceL2R2 = m_FgCellDistance * F1B4;
          FgCellDistanceLR = m_FgCellDistance * F1B2;

          nghbrL = c_neighborId(cellId, 2 * i);
          nghbrL2 = c_neighborId(nghbrL, 2 * i);
          nghbrR = c_neighborId(cellId, 2 * i + 1);
          nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);

          if((!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR2, set)
              || (a_isGBoundaryCellG(nghbrL2, set) || a_isGBoundaryCellG(nghbrR2, set)))) {
            // reduce to second order
            nghbrR2 = nghbrR;
            nghbrL2 = nghbrL;
            FgCellDistanceL2R2 = FgCellDistanceLR;
          }

          // reduce to second order on boundaries (like for the Landau simulations)
          if(a_hasNeighbor(nghbrL, 2 * i) == 0 || a_hasNeighbor(nghbrR, 2 * i + 1) == 0) {
            // reduce to second order
            nghbrR2 = nghbrR;
            nghbrL2 = nghbrL;
            FgCellDistanceL2R2 = FgCellDistanceLR;
          }
          /*
            if(a_hasNeighbor( cellId ,  2*i )==0 ||
            a_hasNeighbor( cellId ,  2*i+1 )==0
            ) {
            // reduce to second order
            nghbrR2 = nghbrR;
            nghbrL2 = nghbrL;
            FgCellDistanceL2R2 = FgCellDistanceLR;
            }
          */

          // take onesided differences directly on the boundary
          // no neighbors in positive direction
          if(a_hasNeighbor(cellId, 2 * i + 1) == 0) {
            oneSidedL = true;

            a_levelSetFunctionSlope(cellId, i, set) =

                -m_FgCellDistance *

                (a * a_levelSetFunctionG(cellId, set) + b * a_levelSetFunctionG(nghbrL, set)
                 + c * a_levelSetFunctionG(nghbrL2, set));
          }

          // no neighbors in negative direction
          if(a_hasNeighbor(cellId, 2 * i) == 0) {
            oneSidedR = true;

            a_levelSetFunctionSlope(cellId, i, set) =

                m_FgCellDistance *

                (a * a_levelSetFunctionG(cellId, set) + b * a_levelSetFunctionG(nghbrR, set)
                 + c * a_levelSetFunctionG(nghbrR2, set));
          }

          // filtering
          if(!oneSidedR && !oneSidedL && m_filterFlameTubeEdges) {
            if(c_coordinate(cellId, 1) < filterCoord) {
              // compute second order curvature values at surfaces (staggered points)

              levelSlopeR = (a_levelSetFunctionG(nghbrR, set) - a_levelSetFunctionG(cellId, set)) * m_FgCellDistance;

              levelSlopeL = (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(nghbrL, set)) * m_FgCellDistance;

              // compute via second order interpolation the curvature back to the cell centers
              a_levelSetFunctionSlope(cellId, i, set) = (F1B2 * (levelSlopeR + levelSlopeL));
              filter = true;
            }
          }

          if(!filter) {
            if(!oneSidedR && !oneSidedL) {
              a_levelSetFunctionSlope(cellId, i, set) =

                  F4B3 * FgCellDistanceLR *

                      (a_levelSetFunctionG(nghbrR, set) - a_levelSetFunctionG(nghbrL, set))

                  - F1B3 * FgCellDistanceL2R2 *

                        (a_levelSetFunctionG(nghbrR2, set) - a_levelSetFunctionG(nghbrL2, set));
            }
          }
          FgradG = POW2(a_levelSetFunctionSlope(cellId, 0, set));
          for(MInt k = 1; k < nDim; k++)
            FgradG += POW2(a_levelSetFunctionSlope(cellId, k, set));

          FgradG = F1 / mMax(epsilon, sqrt(FgradG));

          for(MInt k = 0; k < nDim; k++)
            a_normalVectorG(cellId, k, set) = -FgradG * a_levelSetFunctionSlope(cellId, k, set);
        }
      }
    }
  } else {
    for(MInt set = startSet; set < endSet; set++) {
      if(!m_computeSet[set]) continue;

      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);
        if(a_isHalo(cellId)) continue;
        // second order code
        /*
                for( MInt i=0; i<nDim; i++ )
                  a_levelSetFunctionSlope( cellId ,i, set) =
                    ( a_levelSetFunctionG( c_neighborId( cellId ,  2*i+1 ) , set)  -
                      a_levelSetFunctionG( c_neighborId( cellId ,  2*i ) , set)  ) * Fdx;
        */
        // Valgrind fix !! One-sided differences if one neighbor missing. It should be checked why this situation can
        // occur!!!!
        for(MInt i = 0; i < nDim; i++) {
          MInt n0 = c_neighborId(cellId, 2 * i);
          MInt n1 = c_neighborId(cellId, 2 * i + 1);
          if(n0 < 0) n0 = cellId;
          if(n1 < 0) n1 = cellId;
          a_levelSetFunctionSlope(cellId, i, set) = (a_levelSetFunctionG(n1, set) - a_levelSetFunctionG(n0, set))
                                                    / (c_coordinate(n1, i) - c_coordinate(n0, i));
        }

        FgradG = POW2(a_levelSetFunctionSlope(cellId, 0, set));
        for(MInt i = 1; i < nDim; i++)
          FgradG += POW2(a_levelSetFunctionSlope(cellId, i, set));
        FgradG = F1 / mMax(epsilon, sqrt(FgradG));
        for(MInt i = 0; i < nDim; i++)
          a_normalVectorG(cellId, i, set) = -FgradG * a_levelSetFunctionSlope(cellId, i, set);
      }
    }
  }

  // exchange normal vectors on halo cells, cause they are needed for the curvature computation!
  exchangeDataLS(&a_normalVectorG(0, 0, 0), m_maxNoSets * nDim);
}


// ----------------------------------------------------------------------------------------

/** \brief computes the normal in the level set band for periodic boundaries
 *
 * \author Stephan Schlimpert
 * \date July 2011
 *
 * computes the normal with different order of discretization, the same order is also used
 * for the computation of the curvature
 *
 * \todo labels:LS implementation of periodic boundaries for all discretizations !!! (only done for
 * fourthOrderNormalCurvatureComputation2)
 *
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::computeNormalVectorsPeriodic() {
  TRACE();

  MFloat FgradG;
  MFloat epsilon = m_gCellDistance * 0.00000001;
  const MFloat Fdx = F1B2 * m_FgCellDistance;
  MInt cellId, nghbrL, nghbrL2, nghbrR, nghbrR2;
  // MBool boundary= false;
  MBool oneSidedR = false;
  MBool oneSidedL = false;
  MFloat FgCellDistanceL2R2 = m_FgCellDistance * F1B4, FgCellDistanceLR = m_FgCellDistance * F1B2;

  //---

  if(m_fourthOrderNormalCurvatureComputation) {
    for(MInt set = 0; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);

        for(MInt i = 0; i < nDim; i++) {
          oneSidedR = false;
          oneSidedL = false;

          FgCellDistanceL2R2 = m_FgCellDistance * F1B4;
          FgCellDistanceLR = m_FgCellDistance * F1B2;

          nghbrL = c_neighborId(cellId, 2 * i);
          nghbrL2 = c_neighborId(nghbrL, 2 * i);
          nghbrR = c_neighborId(cellId, 2 * i + 1);
          nghbrR2 = c_neighborId(nghbrR, 2 * i + 1);

          if(a_hasNeighbor(nghbrL2, 2 * i) != 0 && a_hasNeighbor(nghbrR2, 2 * i + 1) != 0) {
            if((!a_inBandG(nghbrL2, set) || !a_inBandG(nghbrR2, set))) {
              // reduce to second order
              nghbrR2 = nghbrR;
              nghbrL2 = nghbrL;
              FgCellDistanceL2R2 = FgCellDistanceLR;
            }
          }

          if(!oneSidedR && !oneSidedL) {
            a_levelSetFunctionSlope(cellId, i, set) =

                F4B3 * FgCellDistanceLR *

                    (a_levelSetFunctionG(nghbrR, set) - a_levelSetFunctionG(nghbrL, set))

                - F1B3 * FgCellDistanceL2R2 *

                      (a_levelSetFunctionG(nghbrR2, set) - a_levelSetFunctionG(nghbrL2, set));
          }
        }
        FgradG = POW2(a_levelSetFunctionSlope(cellId, 0, set));
        for(MInt i = 1; i < nDim; i++)
          FgradG += POW2(a_levelSetFunctionSlope(cellId, i, set));

        FgradG = F1 / mMax(epsilon, sqrt(FgradG));

        for(MInt i = 0; i < nDim; i++)
          a_normalVectorG(cellId, i, set) = -FgradG * a_levelSetFunctionSlope(cellId, i, set);
      }
    }
  } else {
    for(MInt set = 0; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      for(MInt id = 0; id < a_noBandCells(set); id++) {
        cellId = a_bandCellId(id, set);

        // second order code

        for(MInt i = 0; i < nDim; i++)
          a_levelSetFunctionSlope(cellId, i, set) = (a_levelSetFunctionG(c_neighborId(cellId, 2 * i + 1), set)
                                                     - a_levelSetFunctionG(c_neighborId(cellId, 2 * i), set))
                                                    * Fdx;

        FgradG = POW2(a_levelSetFunctionSlope(cellId, 0, set));
        for(MInt i = 1; i < nDim; i++)
          FgradG += POW2(a_levelSetFunctionSlope(cellId, i, set));
        FgradG = F1 / mMax(epsilon, sqrt(FgradG));
        for(MInt i = 0; i < nDim; i++)
          a_normalVectorG(cellId, i, set) = -FgradG * a_levelSetFunctionSlope(cellId, i, set);
      }
    }
  }

  // exchange normal vectors on halo cells, cause they are needed for the curvature computation!
  exchangeDataLS(&a_normalVectorG(0, 0, 0), m_maxNoSets * nDim);
}


// ----------------------------------------------------------------------------------------


// author: Daniel Hartmann, 2007
template <MInt nDim>
void LsCartesianSolver<nDim>::computeNormalVectorsAtFront() {
  TRACE();

  MFloat FgradG;
  const MFloat Fdx = F1B2 * m_FgCellDistance;
  //---

  for(MInt set = 0; set < m_noSets; set++) {
    if(!m_computeSet[set]) continue;
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      for(MInt i = 0; i < nDim; i++)
        a_levelSetFunctionSlope(a_G0CellId(id, set), i, set) =
            (a_levelSetFunctionG(c_neighborId(a_G0CellId(id, set), 2 * i + 1), set)
             - a_levelSetFunctionG(c_neighborId(a_G0CellId(id, set), 2 * i), set))
            * Fdx;

      FgradG = POW2(a_levelSetFunctionSlope(a_G0CellId(id, set), 0, set));
      for(MInt i = 1; i < nDim; i++)
        FgradG += POW2(a_levelSetFunctionSlope(a_G0CellId(id, set), i, set));
      FgradG = F1 / sqrt(FgradG);
      for(MInt i = 0; i < nDim; i++)
        a_normalVectorG(a_G0CellId(id, set), i, set) = -FgradG * a_levelSetFunctionSlope(a_G0CellId(id, set), i, set);
    }
  }
}

//-----------------------------------------------------------------------------------------


/**
 *  \brief Extends the extensionVelocity from the G0 cells to the rest of the band cells
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date July 2020
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::extendVelocity(const MInt set) {
  // create the cell lists
  // cell list: internal cells without cut and bndry cells
  MInt cellListSize = 0;

  // reset scratch array
  MIntScratchSpace gWindowCell(m_maxNoCells, AT_, "gWindowCell");
  gWindowCell.fill(0);

  fill_n(m_cellList, m_maxNoCells, -1);

  for(MInt id = 0; id < a_noBandCells(set); id++) {
    const MInt cellId = a_bandCellId(id, set);
    gWindowCell[cellId] = 0;

    if(a_isGZeroCell(cellId, set)) {
      continue;
    }

    if(a_level(cellId) != a_maxGCellLevel()) {
      continue;
    }

    if(a_isHalo(cellId)) {
      continue;
    }

    MInt count = 0;
    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      if(a_hasNeighbor(cellId, dirId) && a_hasNeighbor(c_neighborId(cellId, dirId), dirId)
         && !a_isGBoundaryCellG(cellId, set) && !a_isGBoundaryCellG(c_neighborId(cellId, dirId), set)) {
        count++;
      }
    }

    if(count == m_noDirs) {
      m_cellList[cellListSize++] = cellId;
      if(a_isWindow(cellId)) {
        gWindowCell[cellId] = 1;
      }
    }
  }

  // exchange:
  // add halo cells to the list
  MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");
  // gather:
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferSize.p[i] = 0;
    for(MInt j = 0; j < noWindowCells(i); j++) {
      m_intSendBuffers[i][sendBufferSize.p[i]++] = gWindowCell[windowCellId(i, j)];
    }
  }
  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        m_intSendBuffers[i][sendBufferSize.p[i]++] = gWindowCell[grid().azimuthalWindowCell(i, j)];
      }
    }
  }

  // exchange data -> send, receive
  exchangeIntBuffers(sendBufferSize.getPointer(), receiveBufferSize.getPointer(), 11, 1);

  // scatter:
  // update the list and add halo cells which are needed for hyperbolic extension
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
    // check, if received data size matches expected size:
    if(receiveBufferSize.p[i] != noHaloCells(i))
      mTerm(1, AT_, "this was not expected to happen: wrong number of halo information...");
#endif
    for(MInt j = 0; j < noHaloCells(i); j++) {
      MInt haloCell = haloCellId(i, j);
      if(m_intReceiveBuffers[i][j] == 1) {
        m_cellList[cellListSize++] = haloCell;
      }
    }
    // cerr << "number of cells " << haloCnt << " added from neighbor domain " << grid().neighborDomain(i) << endl;
  }

  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      MInt offset = noHaloCells(grid().azimuthalNeighborDomain(i));
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        MInt n = offset + j;
        MInt haloCell = grid().azimuthalHaloCell(i, j);
        if(m_intReceiveBuffers[i][n] == 1) {
          m_cellList[cellListSize++] = haloCell;
        }
      }
    }
  }

  // TODO labels:LS maxNoCells vs a_noCells
  MFloatScratchSpace velocity(m_maxNoCells, AT_, "velocity");

  for(MInt n = 0; n < nDim; n++) {
    velocity.fill(0.0);

    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      const MInt cellId = a_G0CellId(id, set);
      if(a_isHalo(cellId)) {
        continue;
      }

      velocity[cellId] = a_extensionVelocityG(cellId, n, set);
    }

    for(MInt cell = 0; cell < cellListSize; cell++) {
      const MInt cellId = m_cellList[cell];
      velocity[cellId] = F0;
    }

    hyperbolicExtensionOpt(velocity.getPointer(), m_cellList, cellListSize, m_extVelConvergence, set);

    // TODO labels:LS Why? Seems redundant
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      const MInt cellId = a_G0CellId(id, set);
      a_extensionVelocityG(cellId, n, set) = velocity[cellId];
    }

    for(MInt cell = 0; cell < cellListSize; cell++) {
      const MInt cellId = m_cellList[cell];
      a_extensionVelocityG(cellId, n, set) = velocity[cellId];
    }
  }
}

// ----------------------------------------------------------------------------------------
/** \brief computes extension velocity for a flame
 *
 * \author Stephan Schlimpert
 * \date Februar 2012
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::computeExtensionVelocityGEQUPVMarksteinOpt(MFloat* FfluidDensity, MInt set) {
  TRACE();

  ASSERT(m_lsCollectorMode == 0, "");

  MInt cellId;
  MInt cellListSize;
  MFloatScratchSpace q(m_maxNoCells, AT_, "q");
  MIntScratchSpace gWindowCell(m_maxNoCells, AT_, "gWindowCell");
  MInt count = 0;
  MFloat FdampingDistanceFlameBase = F1 / m_dampingDistanceFlameBase;
  MFloat dampCoord = m_dampingDistanceFlameBase + m_yOffsetFlameTube;
  MFloat FextensionDampingDistanceFlameBase = F1 / m_dampingDistanceFlameBaseExtVel;
  MFloat extensionDampCoord = m_dampingDistanceFlameBaseExtVel + m_yOffsetFlameTube;
  MFloat levelSetPlus = F0, levelSetNegative = F0; //, curvatureFactor=F0;
  //---

  // create the cell lists
  // cell list: internal cells without cut and bndry cells
  // bc cells : boundary cells
  cellListSize = 0;

  // reset scratch array
  q.fill(F0);
  gWindowCell.fill(0);

  for(MInt c = 0; c < m_maxNoCells; c++) {
    m_cellList[c] = -1;
  }

  for(MInt id = 0; id < a_noBandCells(set); id++) {
    cellId = a_bandCellId(id, set);
    gWindowCell[cellId] = 0;
    if(a_isGZeroCell(cellId, set)) continue;

    count = 0;
    if(a_level(cellId) != a_maxGCellLevel()) continue;

    if(a_isHalo(cellId)) continue;

    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      if(a_hasNeighbor(cellId, dirId) && a_hasNeighbor(c_neighborId(cellId, dirId), dirId) &&
         //         a_hasNeighbor( cellId ,  opposite[dirId] ) &&
         // a_hasNeighbor( c_neighborId( cellId ,  opposite[dirId] )  ,  opposite[dirId] ) &&
         !a_isGBoundaryCellG(cellId, set) && !a_isGBoundaryCellG(c_neighborId(cellId, dirId), set)) {
        count++;
      }
    }
    if(count == m_noDirs) {
      m_cellList[cellListSize++] = cellId;
      // for debugging
      // a_stretchG( cellId , 0) = 1;
      if(a_isWindow(cellId)) gWindowCell[cellId] = 1;
    }
  }

  // exchange:
  // add halo cells to the list
  MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");
  // gather:
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferSize.p[i] = 0;
    for(MInt j = 0; j < noWindowCells(i); j++) {
      m_intSendBuffers[i][sendBufferSize.p[i]++] = gWindowCell[windowCellId(i, j)];
    }
  }
  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        m_intSendBuffers[i][sendBufferSize.p[i]++] = gWindowCell[grid().azimuthalWindowCell(i, j)];
      }
    }
  }

  // exchange data -> send, receive
  exchangeIntBuffers(sendBufferSize.getPointer(), receiveBufferSize.getPointer(), 11, 1);

  // scatter:
  // update the list and add halo cells which are needed for hyperbolic extension
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
    //      MInt haloCnt = 0;
    // check, if received data size matches expected size:
    if(receiveBufferSize.p[i] != noHaloCells(i))
      mTerm(1, AT_, "this was not expected to happen: wrong number of halo information...");
#endif
    for(MInt j = 0; j < noHaloCells(i); j++) {
      MInt haloCell = haloCellId(i, j);
      if(m_intReceiveBuffers[i][j] == 1) {
        m_cellList[cellListSize++] = haloCell;
        // for debugging
        // a_stretchG( haloCell , 0) = 1;
        //  haloCnt++;
      }
    }
    //      m_log << "number of cells " << haloCnt << " added from neighbor domain " << grid().neighborDomain(i) <<
    //      endl;
    // cerr << "number of cells " << haloCnt << " added from neighbor domain " << grid().neighborDomain(i) << endl;
  }

  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      MInt offset = noHaloCells(grid().azimuthalNeighborDomain(i));
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        MInt n = offset + j;
        MInt haloCell = grid().azimuthalHaloCell(i, j);
        if(m_intReceiveBuffers[i][n] == 1) {
          m_cellList[cellListSize++] = haloCell;
        }
      }
    }
  }


  // compute the extension velocity of G0 cut cells
  for(MInt id = 0; id < a_noG0Cells(set); id++) {
    cellId = a_G0CellId(id, set);

    a_correctedBurningVelocity(cellId, set) =
        a_flameSpeedG(cellId, set) * (F1 - a_curvatureG(cellId, set) * m_marksteinLength);

    for(MInt i = 0; i < nDim; i++) {
      a_extensionVelocityG(cellId, i, set) += a_normalVectorG(cellId, i, set) * m_rhoFlameTube * FfluidDensity[id]
                                              * a_correctedBurningVelocity(cellId, set);
    }
  }

  if(cellListSize < a_noG0Cells(set)) {
    cerr << "Error: number of cells in cell list is smaller than number of G0 cells" << cellListSize
         << " , noG0cells: " << a_noG0Cells(set) << endl;
  }

  IF_CONSTEXPR(nDim == 2) {
    // iterative extension
    for(MInt i = 0; i < nDim; i++) {
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue;

        q.p[cellId] = a_extensionVelocityG(cellId, i, set);
        if(c_coordinate(cellId, 1) < extensionDampCoord) {
          q.p[cellId] = a_extensionVelocityG(cellId, i, set) * FextensionDampingDistanceFlameBase
                        * (c_coordinate(cellId, 1) - m_yOffsetFlameTube);
        }
      }
      // needed otherwise solution for previous space direction is used
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        q.p[cellId] = F0;
      }

      hyperbolicExtensionOpt(q.getPointer(), m_cellList, cellListSize, m_extVelConvergence, set);

      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        a_extensionVelocityG(cellId, i, set) = q.p[cellId];
      }

      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        a_extensionVelocityG(cellId, i, set) = q.p[cellId];
      }
    }

    if(m_useCorrectedBurningVelocity) {
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue;
        q.p[cellId] = a_correctedBurningVelocity(cellId, set);
      }
      hyperbolicExtensionOpt(q.getPointer(), m_cellList, cellListSize, m_extVelConvergence, set);

      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        a_correctedBurningVelocity(cellId, set) = q.p[cellId];
      }
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        a_correctedBurningVelocity(cellId, set) = q.p[cellId];
      }
    }

    if(m_hyperbolicCurvature)
      hyperbolicExtensionOpt(&a_curvatureG(0, 0), m_cellList, cellListSize, m_extVelConvergence, set);

    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);

      if(c_coordinate(cellId, 1) < m_yOffsetFlameTube) {
        a_curvatureG(cellId, set) = F0;
        for(MInt i = 0; i < nDim; i++) {
          a_extensionVelocityG(cellId, i, set) = F0;
        }
        continue;
      }
      // linearly damping of curvature to flame base
      if(c_coordinate(cellId, 1) < dampCoord) {
        a_curvatureG(cellId, set) *= FdampingDistanceFlameBase * (c_coordinate(cellId, 1) - 0.0234);
        levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
        levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);
        a_curvatureG(cellId, set) *= F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
      }
      // linearly damping of extension velocity to flame base
      if(c_coordinate(cellId, 1) < extensionDampCoord) {
        for(MInt i = 0; i < nDim; i++) {
          a_extensionVelocityG(cellId, i, set) *=
              FextensionDampingDistanceFlameBase * (c_coordinate(cellId, 1) - m_yOffsetFlameTube);
        }
      }
    }
    // 3D
  }
  else {
    // iterative extension
    for(MInt i = 0; i < nDim; i++) {
      /*      MFloat tmp=F0;
      MFloat tmpM=2147483647;
      */
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue;


        q.p[cellId] = a_extensionVelocityG(cellId, i, set);

        if(c_coordinate(cellId, 1) < extensionDampCoord) {
          q.p[cellId] *= FextensionDampingDistanceFlameBase * (c_coordinate(cellId, 1) - m_yOffsetFlameTube);
        }
        if(c_coordinate(cellId, 1) < m_yOffsetFlameTube) {
          q.p[cellId] = F0;
        }
      }
      // needed otherwise solution for previous space direction is used
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        q.p[cellId] = F0;
      }

      hyperbolicExtensionOpt(q.getPointer(), m_cellList, cellListSize, m_extVelConvergence, set);

      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        a_extensionVelocityG(cellId, i, set) = q.p[cellId];
      }

      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        a_extensionVelocityG(cellId, i, set) = q.p[cellId];
      }
    }

    if(m_useCorrectedBurningVelocity) {
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue;
        q.p[cellId] = a_correctedBurningVelocity(cellId, set);
      }
      hyperbolicExtensionOpt(q.getPointer(), m_cellList, cellListSize, m_extVelConvergence, set);

      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        cellId = a_G0CellId(id, set);
        a_correctedBurningVelocity(cellId, set) = q.p[cellId];
      }
      for(MInt cell = 0; cell < cellListSize; cell++) {
        cellId = m_cellList[cell];
        a_correctedBurningVelocity(cellId, set) = q.p[cellId];
      }
    }

    if(m_hyperbolicCurvature)
      hyperbolicExtensionOpt(&a_curvatureG(0, 0), m_cellList, cellListSize, m_extVelConvergence, set);

    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);

      if(c_coordinate(cellId, 1) < m_yOffsetFlameTube) {
        a_curvatureG(cellId, set) = F0;
        for(MInt i = 0; i < nDim; i++) {
          a_extensionVelocityG(cellId, i, set) = F0;
        }
        continue;
      }
      // linearly damping of curvature to flame base
      if(c_coordinate(cellId, 1) < dampCoord) {
        a_curvatureG(cellId, set) *= FdampingDistanceFlameBase * (c_coordinate(cellId, 1) - m_yOffsetFlameTube);
        levelSetNegative = a_levelSetFunctionG(cellId, set) - (m_noReactionCells / m_curvatureDampFactor);
        levelSetPlus = a_levelSetFunctionG(cellId, set) + (m_noReactionCells / m_curvatureDampFactor);
        a_curvatureG(cellId, set) *= F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
      }
      // linearly damping of extension velocity to flame base
      if(c_coordinate(cellId, 1) < extensionDampCoord) {
        for(MInt i = 0; i < nDim; i++) {
          a_extensionVelocityG(cellId, i, set) *=
              FextensionDampingDistanceFlameBase * (c_coordinate(cellId, 1) - m_yOffsetFlameTube);
        }
      }
    }
  }

  DEBUG("LsCartesianSolver::computeExtensionVelocityGEQUPVMarksteinOpt return", MAIA_DEBUG_TRACE_OUT);
}

// ----------------------------------------------------------------------------------------


/** \brief hyperbolic extension of a quantity along the normals without applying outer bounday conditions
 *
 * \author Stephan Schlimpert
 * \date Februar 2012, January 2013
 *
 * This function extends a quantity along the normals of the G-field starting from an arbitrary contour
 *
 * last change: including minIteration to assure full extension of velocity field! + parallelization
 *
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::hyperbolicExtensionOpt(MFloat* q, MInt* cellList, MInt cellListSize,
                                                     MFloat convergenceCriterion, MInt set) {
  TRACE();
  m_log << "ERROR, WARNING, you are using hyperbolicExtensionOpt... look at the code and make sure the "
           "'optimization' below '//euler solver' is not commented since it is wrong"
        << endl;

  MInt cellId;
  MInt iteration;
  MInt nghbrL, nghbrR;
  MFloat dtEULER = m_extVelCFL * m_gCellDistance;
  MFloat eps = 0.000000001;
  MFloat qOld = F0;
  MInt minIteration;
  MFloatScratchSpace levelSetRHS(m_maxNoCells, AT_, "levelSetRHS");

  // compute minimum iterations to transport the velocity to the band cells
  minIteration = m_gBandWidth / m_extVelCFL;

  //---

  // initialize
  iteration = 0;
  MFloat res = 0;

  for(MInt cell = 0; cell < cellListSize; cell++) {
    cellId = cellList[cell];
    m_signG[IDX_LSSET(cellId, set)] =
        a_levelSetFunctionG(cellId, set) / mMax(eps, fabs(a_levelSetFunctionG(cellId, set)));
  }

  exchangeLs(q, set, 1);

  // Stephan: this condition assures a fully extension of the velocity with a certain convergence criterion,
  //          - minIteration: necessary because otherwise a wrong user defined extVelIteration number could lead to
  //          unwanted velocity oscillations
  //          - checking another time convergenceCriterion is necessary otherwise the minimum number of iterations could
  //          be reached but the convergence criterion is not reached!
  // TODO labels:LS,totest check wether extVelIterations could be removed?!
  // TODO labels:LS,toenhance unified method of reinitialization
  while(iteration < m_extVelIterations || iteration < minIteration || res < convergenceCriterion) {
    // apply outer boundary condition d(fExt_i) / d(x_i) = 0
    /*
      for( MInt cell = 0; cell < noBcCells;  cell++ ) {
      cellId = bcCellList[ cell ];
      q[ cellId ] = bc*q[ linkedInternalCell[cellId] ];
      }
    */
    res = F0;
    for(MInt cell = 0; cell < cellListSize; cell++) {
      cellId = cellList[cell];
      levelSetRHS[cellId] = F0;
      if(a_isHalo(cellId)) continue;
      for(MInt i = 0; i < nDim; i++) {
        nghbrL = c_neighborId(cellId, 2 * i);     //[ cellId*m_noDirs + 2*i ];
        nghbrR = c_neighborId(cellId, 2 * i + 1); // nghbrs[ cellId*m_noDirs + 2*i+1 ];
        levelSetRHS[cellId] += mMax(m_signG[IDX_LSSET(cellId, set)] * (-a_normalVectorG(cellId, i, set)), F0)
                                   * (q[cellId] - q[nghbrL]) * m_FgCellDistance
                               + mMin(m_signG[IDX_LSSET(cellId, set)] * (-a_normalVectorG(cellId, i, set)), F0)
                                     * (q[nghbrR] - q[cellId]) * m_FgCellDistance;
      }


      // Euler solver //this commented part may not be commented... this "optimization" by stephan introduces a
      // direction dependency of the solution. This has to be uncommented before using this function!!
      /*    for( MInt cell = 0; cell < cellListSize;  cell++ ) {
        cellId = cellList[ cell ];
        if( a_isHalo( cellId ) )
          continue;
      */
      qOld = q[cellId];
      q[cellId] = qOld - dtEULER * levelSetRHS[cellId];
      res = mMax(res, ABS(q[cellId] - qOld));
    }

    exchangeLs(q, set, 1);

    iteration++;

    MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "res");

    switch(m_levelSetBoundaryCondition) {
      default: {
        // check convergence, if minimum number of iteration is reached
        if(res < convergenceCriterion) { //&& iteration> minIteration){

          if(domainId() == 0) {
            FILE* datei;
            stringstream hyperbolicFile;
            hyperbolicFile << "hyperbolicExtension";
            if(m_noSets > 1) hyperbolicFile << "_s" << set;
            hyperbolicFile << "_" << domainId();
            datei = fopen((hyperbolicFile.str()).c_str(), "a+");
            fprintf(datei, "  %d", globalTimeStep);
            fprintf(datei, "  %-10.10f", res);
            fprintf(datei, "  %d", iteration);
            fprintf(datei, "\n");
            fclose(datei);
          }


          return iteration;
        }
      }
    }
  }

  if(domainId() == 0) {
    FILE* datei;
    stringstream hyperbolicFile;
    hyperbolicFile << "hyperbolicExtension";
    if(m_noSets > 1) hyperbolicFile << "_s" << set;
    hyperbolicFile << "_" << domainId();
    datei = fopen((hyperbolicFile.str()).c_str(), "a+");
    fprintf(datei, "  %d", globalTimeStep);
    fprintf(datei, "  %-10.10f", res);
    fprintf(datei, "  %d", iteration);
    fprintf(datei, "\n");
    fclose(datei);
  }

  return iteration;
}


//---------------------------------------------------------------------------


/** \brief computes the gcell time step
 *
 * \author Daniel Hartmann
 * \date July 2007
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeGCellTimeStep() {
  TRACE();
  if(m_combustion) {
  } else if(m_levelSetMb || m_levelSetFv) {
    // NOTE: the timeStep is copied from the fv-solver in transferTimeStep in the coupler!
    return;
  } else if(m_freeSurface) {
    // Note: the timeStep is set by the flow solver
    return;
  } else {
    ASSERT(m_timeStepMethod == 8, "");
    m_timeStep = m_cfl * m_gCellDistance;
  }
}


//-----------------------------------------------------------------------------


/** \brief Updates coarser G-grid levels
 *
 * This function updates the level set function, the normal vector, the extension velocity, and the burning velocity.
 * For efficiency it is assumed that each cell has 2^d children.
 * Furthermore, it is important that the lower level cells are completely covered by in-band children.
 *
 * \author Daniel Hartmann
 * \date 2008
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetRestriction() {
  TRACE();

  MInt cellId;
  MFloat factor;
  MFloat FnoChildren = F1 / FPOW2(nDim);
  // ---

  for(MInt set = 0; set < m_noSets; set++) {
    // reset all coarse grid cells
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      while(c_parentId(cellId) > -1) {
        cellId = c_parentId(cellId);
        a_levelSetFunctionG(cellId, set) = 0;
        if(!m_semiLagrange) {
          a_correctedBurningVelocity(cellId, set) = 0;
          a_curvatureG(cellId, set) = F0;
          for(MInt i = 0; i < nDim; i++) {
            a_extensionVelocityG(cellId, i, set) = 0;
            a_normalVectorG(cellId, i, set) = 0;
          }
        }
      }
    }

    // update all coarse grid cells
    for(MInt id = 0; id < a_noBandCells(set); id++) {
      cellId = a_bandCellId(id, set);
      factor = F1;
      while(c_parentId(cellId) > -1) {
        cellId = c_parentId(cellId);
        factor *= FnoChildren;
        a_levelSetFunctionG(cellId, set) += factor * a_levelSetFunctionG(a_bandCellId(id, set), set);
        if(!m_semiLagrange) {
          a_correctedBurningVelocity(cellId, set) += factor * a_correctedBurningVelocity(a_bandCellId(id, set), set);
          a_curvatureG(cellId, set) += factor * a_curvatureG(a_bandCellId(id, set), set);

          for(MInt i = 0; i < nDim; i++) {
            a_extensionVelocityG(cellId, i, set) += factor * a_extensionVelocityG(a_bandCellId(id, set), i, set);
            a_normalVectorG(cellId, i, set) += factor * a_normalVectorG(a_bandCellId(id, set), i, set);
          }
        }
      }
    }
  }
}


//-----------------------------------------------------------------------------


// computations are only based on zeroth level-set function!
/** \brief Computes the zero level set arc length (2D) using a first-order discrete delta function
 *
 * \authors Daniel Hartmann, Stephan Schlimpert
 * \date unknown, 16.06.2011, November 2011, Februar 2012, Februrar 2013
 *
 * last change: parallelization
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeZeroLevelSetArcLength() {
  TRACE();

  MFloat deltaFunction; //,deltaFunctionR;
  MFloat Dx0, Dxplus, Dxminus, Dy0, Dyplus, Dyminus, sigma, rhoS_L;

  const MFloat eps = F1 / FPOW10[10];
  MInt cellId;
  //  MFloat a=-F23B24,b=F7B8,c1=F1B8,d=-F1B24;
  // MBool boundary = false;
  MInt set = 0;

  if(!m_highOrderDeltaFunction) {
    m_arcLength = F0;
    m_massConsumption = F0;

    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      deltaFunction = F0;
      cellId = a_G0CellId(id, set);
      if(a_isHalo(cellId)) continue;

      if(a_level(cellId) != a_maxGCellLevel()) continue;
      // compute the discrete derivatives
      Dx0 = (a_levelSetFunctionG(c_neighborId(cellId, 1), set) - a_levelSetFunctionG(c_neighborId(cellId, 0), set))
            * m_FgCellDistance * F1B2;
      Dxplus =
          (a_levelSetFunctionG(c_neighborId(cellId, 1), set) - a_levelSetFunctionG(cellId, set)) * m_FgCellDistance;
      Dxminus =
          (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(c_neighborId(cellId, 0), set)) * m_FgCellDistance;
      Dy0 = (a_levelSetFunctionG(c_neighborId(cellId, 3), set) - a_levelSetFunctionG(c_neighborId(cellId, 2), set))
            * m_FgCellDistance * F1B2;
      Dyplus =
          (a_levelSetFunctionG(c_neighborId(cellId, 3), set) - a_levelSetFunctionG(cellId, set)) * m_FgCellDistance;
      Dyminus =
          (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(c_neighborId(cellId, 2), set)) * m_FgCellDistance;
      IF_CONSTEXPR(nDim == 3) {
        MFloat Dz0 = F0;
        MFloat Dzplus = F0;
        MFloat Dzminus = F0;
        Dz0 = (a_levelSetFunctionG(c_neighborId(cellId, 5), set) - a_levelSetFunctionG(c_neighborId(cellId, 4), set))
              * m_FgCellDistance * F1B2;
        Dzplus =
            (a_levelSetFunctionG(c_neighborId(cellId, 5), set) - a_levelSetFunctionG(cellId, set)) * m_FgCellDistance;
        Dzminus =
            (a_levelSetFunctionG(cellId, set) - a_levelSetFunctionG(c_neighborId(cellId, 4), set)) * m_FgCellDistance;
        sigma = sqrt(POW2(Dx0) + POW2(Dy0) + POW2(Dz0) + eps);

        // compute delta function (z contribution)
        if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 4), set) < F0)
          deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 4), set) * Dz0)
                           / (ABS(Dzminus) * ABS(sigma) * POW2(m_gCellDistance));
        if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 5), set) < F0)
          deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 5), set) * Dz0)
                           / (ABS(Dzplus) * ABS(sigma) * POW2(m_gCellDistance));
      }
      else {
        sigma = sqrt(POW2(Dx0) + POW2(Dy0) + eps);
      } /*
          if(a_levelSetFunctionG( cellId , set)  == F0 ||
          a_levelSetFunctionG( m_nghbrIdsG[ m_noDirs*cellId ] , set)  == F0 ||
          a_levelSetFunctionG( m_nghbrIdsG[ m_noDirs*cellId ] + 1, set)  == F0 ||
          a_levelSetFunctionG( m_nghbrIdsG[ m_noDirs*cellId ] + 2 , set)  == F0 ||
          a_levelSetFunctionG( m_nghbrIdsG[ m_noDirs*cellId ] + 3 , set)  == F0
          ){
          cerr << "arclength could be affected by levelsetfunction == 0" << endl;
          }*/
      // compute the delta function (x and y contributions)
      if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 0), set) < F0)
        deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 0), set) * Dx0)
                         / (ABS(Dxminus) * ABS(sigma) * POW2(m_gCellDistance));
      if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 1), set) < F0)
        deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 1), set) * Dx0)
                         / (ABS(Dxplus) * ABS(sigma) * POW2(m_gCellDistance));
      if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 2), set) < F0)
        deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 2), set) * Dy0)
                         / (ABS(Dyminus) * ABS(sigma) * POW2(m_gCellDistance));
      if(a_levelSetFunctionG(cellId, set) * a_levelSetFunctionG(c_neighborId(cellId, 3), set) < F0)
        deltaFunction += ABS(a_levelSetFunctionG(c_neighborId(cellId, 3), set) * Dy0)
                         / (ABS(Dyplus) * ABS(sigma) * POW2(m_gCellDistance));

      if(m_combustion && m_plenum) {
        IF_CONSTEXPR(nDim == 2) {
          m_arcLength += deltaFunction * POW2(m_gCellDistance);
          rhoS_L = m_rhoFlameTube * a_flameSpeedG(cellId, set) * (F1 - a_curvatureG(cellId, set) * m_marksteinLength);
          // rhoS_L = m_rhoFlameTube * a_correctedBurningVelocity(cellId, set);
          m_massConsumption += rhoS_L * deltaFunction * sigma * POW2(m_gCellDistance);
          // m_massConsumption += a_correctedBurningVelocity(cellId, set) * deltaFunction * sigma * POW2(
          // m_gCellDistance );
        }
        else {
          m_arcLength += deltaFunction * POW3(m_gCellDistance);
          rhoS_L = m_rhoFlameTube * a_flameSpeedG(cellId, set) * (F1 - a_curvatureG(cellId, set) * m_marksteinLength);
          // rhoS_L = m_rhoFlameTube * a_correctedBurningVelocity(cellId, set);
          m_massConsumption += rhoS_L * deltaFunction * sigma * POW3(m_gCellDistance);
          // m_massConsumption += a_correctedBurningVelocity(cellId, set) * deltaFunction * sigma * POW3(
          // m_gCellDistance );
        }
      } else {
        IF_CONSTEXPR(nDim == 2) {
          m_arcLength += deltaFunction * POW2(m_gCellDistance);
          rhoS_L = m_rhoInfinity * a_flameSpeedG(cellId, set) * (F1 - a_curvatureG(cellId, set) * m_marksteinLength);
          m_massConsumption += rhoS_L * deltaFunction * sigma * POW2(m_gCellDistance);
          // m_massConsumption += a_correctedBurningVelocity(cellId, set) * deltaFunction * sigma * POW2(
          // m_gCellDistance );
        }
        else {
          m_arcLength += deltaFunction * POW3(m_gCellDistance);
          rhoS_L = m_rhoInfinity * a_flameSpeedG(cellId, set) * (F1 - a_curvatureG(cellId, set) * m_marksteinLength);
          m_massConsumption += rhoS_L * deltaFunction * sigma * POW3(m_gCellDistance);
          // m_massConsumption += a_correctedBurningVelocity(cellId, set) * deltaFunction * sigma * POW3(
          // m_gCellDistance );
        }
      }
    }
  } else {
    /// \todo labels:LS ask stephan intended behaviour here
    // stringstream errorMessage;
    // errorMessage << "LsCartesianSolver::computeZeroLevelSetArcLength(): switch variable
    // 'm_highOrderDeltaFunction' with value " << m_highOrderDeltaFunction << " not matching any case." << endl;
    // mTerm(1, AT_, errorMessage.str());
  }

  MPI_Allreduce(MPI_IN_PLACE, &m_arcLength, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "m_arcLength");
  MPI_Allreduce(MPI_IN_PLACE, &m_massConsumption, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                "m_massConsumption");
}


//-----------------------------------------------------------------------------


/**\brief Specifies how the different bodies are assigned to different level-set functions
 *
 *   m_bodyToSetTable
 *   m_noBodiesInSet
 *   m_setToBodiesTable
 *   m_startSet
 *   m_noSets
 *   are initialized.
 *
 * \author Claudia Guenther
 * \date 29.02.2012
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setUpBodyToSetTable() {
  TRACE();

  // initialize Tables
  mAlloc(m_bodyToSetTable, m_noBodyBndryCndIds, "m_bodyToSetTable", 0, AT_);
  mAlloc(m_noBodiesInSet, m_maxNoSets, "m_noBodiesInSet", 0, AT_);
  mAlloc(m_setToBodiesTable, m_maxNoSets, m_noBodyBndryCndIds, "m_setToBodiesTable", 0, AT_);
  MIntScratchSpace bodiesinSet(m_noSets, AT_, "bodiesinSet");
  MInt sumbodies = 0;

  // 0: read and allocate properties
  m_startSet = (m_buildCollectedLevelSetFunction ? 1 : 0);
  m_noSets = m_noBodyBndryCndIds + m_startSet; // default
  m_noSets = Context::getSolverProperty<MInt>("nodifferentSets", m_solverId, AT_, &m_noSets);

  // 1: test if property is set
  if(Context::propertyExists("bodiesinSet", m_solverId)) {
    for(MInt i = 0; i < m_noSets; i++) {
      bodiesinSet[i] = Context::getSolverProperty<MInt>("bodiesinSet", m_solverId, AT_, &bodiesinSet[i], i);
      sumbodies = sumbodies + bodiesinSet[i];
    }

    if(sumbodies != m_noBodyBndryCndIds)
      mTerm(1, AT_,
            "The number of bodies devided onto the sets does not fit the total number of moving "
            "boundery-conditions for mode 11!");
    m_noSets = m_noSets + m_startSet; // increased for G0-level

    // 2: setup the levelset table
    MInt count = 0;
    for(MInt j = 0; j < m_noSets; j++) {
      for(MInt i = 0; i < bodiesinSet[j]; i++) {
        m_bodyToSetTable[count] = j + m_startSet;
        m_setToBodiesTable[m_startSet + j][i] = count;
        m_noBodiesInSet[m_startSet + j]++;
        count++;
      }
    }

    // 3: catch if the propery is not set
  } else {
    m_startSet = 1;
    m_noSets = m_noBodyBndryCndIds + m_startSet;
    for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
      m_bodyToSetTable[i] = i + m_startSet;
      m_setToBodiesTable[i + m_startSet][m_noBodiesInSet[i + m_startSet]] = i;
      m_noBodiesInSet[i + m_startSet]++;
    }
  }

  if(m_noSets > 0 && domainId() == 0 && m_noBodyBndryCndIds <= 10) {
    cerr << " noSets: " << m_noSets << endl;
    cerr << " m_bodyToSetTable: ";
    for(MInt i = 0; i < m_noBodyBndryCndIds; i++)
      cerr << " " << m_bodyToSetTable[i];
    cerr << endl;
    cerr << " m_noBodiesInSet: ";
    for(MInt i = 0; i < m_maxNoSets; i++)
      cerr << " " << m_noBodiesInSet[i];
    cerr << endl;
    cerr << " m_setToBodiesTable: ";
    for(MInt i = 0; i < m_maxNoSets; i++) {
      cerr << "s" << i << ": ";
      for(MInt j = 0; j < m_noBodiesInSet[i]; j++)
        cerr << " " << m_setToBodiesTable[i][j];
    }
    cerr << endl;
  }

  for(MInt i = 0; i < m_maxNoSets; i++) {
    if(m_noBodiesInSet[i] > m_noBodyBndryCndIds) mTerm(1, AT_, to_string(m_noBodiesInSet[i]));
    for(MInt j = 0; j < m_noBodyBndryCndIds; j++) {
      if(m_setToBodiesTable[i][j] > m_noEmbeddedBodies) mTerm(1, AT_, to_string(m_setToBodiesTable[i][j]));
    }
  }
  for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
    if(m_bodyToSetTable[i] > m_noSets) mTerm(1, AT_, to_string(m_bodyToSetTable[i]));
  }
}


//-----------------------------------------------------------------------------


/**
 * \fn void LsCartesianSolver::computeBodyPropertiesForced
 * \brief returns a specific property of the specifuec body
 *        used to provide a unique function for both level-set and moving boundary code
 *        return mode:
 *        1: body venter
 *        2: body velocity
 *        3: body acceleration
 *        4: body temperature
 *
 *
 * \author Claudia Guenther, Tim Wegmann
 * \date 03/2011
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::computeBodyPropertiesForced(MInt returnMode, MFloat* bodyData, MInt body,
                                                          MFloat elapsedTime, MBool printPosition) {
  TRACE();

  ASSERT(!m_combustion && !m_freeSurface, "");

  if(m_periodicMovement && globalTimeStep == m_restartTimeStep) {
    determinePeriodicDistance();
  }


  MFloat angle = F0;
  MBool& first = m_static_computeBodyProperties_first;
  MFloat(&amplitude)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_amplitude;
  MFloat(&freqFactor)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_freqFactor;
  MFloat(&initialBodyCenter)[s_maxNoEmbeddedBodies * 3] = m_static_computeBodyProperties_initialBodyCenter;
  MFloat& Strouhal = m_static_computeBodyProperties_Strouhal;
  MFloat(&mu)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_mu;
  MFloat(&mu2)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_mu2;
  MFloat(&liftStartAngle1)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftStartAngle1;
  MFloat(&liftEndAngle1)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftEndAngle1;
  MFloat(&liftStartAngle2)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftStartAngle2;
  MFloat(&liftEndAngle2)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftEndAngle2;
  MFloat(&circleStartAngle)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_circleStartAngle;
  MFloat(&normal)[s_maxNoEmbeddedBodies * 3] = m_static_computeBodyProperties_normal;
  MInt(&bodyToFunction)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_bodyToFunction;
  MFloat& omega = m_static_computeBodyProperties_omega;
  MFloat& rotAngle = m_static_computeBodyProperties_rotAngle;
  MFloat(&temperature)[s_maxNoEmbeddedBodies] = m_static_computeBodyProperties_temperature;


  if(first) {
    const MInt noEmbeddedBodies = m_noEmbeddedBodies;

    if(noEmbeddedBodies > s_maxNoEmbeddedBodies) {
      mTerm(1, AT_, "Error in computeBodyProperties: too many embedded Bodies!");
    }

    // 1: set default values:
    Strouhal = 0.2;
    for(MInt k = 0; k < s_maxNoEmbeddedBodies; k++) {
      amplitude[k] = 0.1;
      freqFactor[k] = 1.0;
      bodyToFunction[k] = 1;
      for(MInt i = 0; i < nDim; i++) {
        initialBodyCenter[k * nDim + i] = F0;
        normal[k * nDim + i] = F0;
      }
      normal[k * nDim + 0] = 1.0;
      liftStartAngle1[k] = F0;
      liftEndAngle1[k] = PI;
      liftStartAngle2[k] = 3.0 * PI;
      liftEndAngle2[k] = 4.0 * PI;
      circleStartAngle[k] = F0;
      temperature[k] = 1;
    }

    if(m_levelSetLb) {
      MFloat MaLb = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_);
      for(MInt k = 0; k < s_maxNoEmbeddedBodies; k++) {
        amplitude[k] = MaLb * LBCS;
      }
    }

    // 2: read Properties

    /*! \page propertiesLS
      \section amplitudes
      <code>MFloat* LsCartesianSolver::m_static_computeBodyProperties_amplitude</code>\n
      default = <code>none</code>\n \n
      Amplitude for body motion of embedded bodies. \n
      NOTE: also used in FV-MB solver for some special cases. \n \n
      Possible values are:
      <ul>
      <li>list of floating point numbers</li>
      </ul>
      Keywords: <i>LEVELSET, MOVING, BODY, BODY_MOTION</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      amplitude[i] = Context::getSolverProperty<MFloat>("amplitudes", m_solverId, AT_, &amplitude[i], i);
    }

    if(m_levelSetLb) {
      MFloat MaLb = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_);
      for(MInt k = 0; k < s_maxNoEmbeddedBodies; k++) {
        amplitude[k] *= MaLb / F1BCS;
      }
    }

    /*! \page propertiesLS
      \section freqFactors
      <code>MFloat* LsCartesianSolver::m_static_computeBodyProperties_freqFactor</code>\n
      default = <code>none</code>\n \n
      Set the frequency factors for prescribing body motion for all embedded bodies. \n
      NOTE: also used in FV-MB solver for some special cases. \n \n
      Possible values are:
      <ul>
        <li>list of positive floating point numbers</li>
      </ul>
      Keywords: <i>LEVELSET, MOVING, BODY, BODY_MOTION</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      freqFactor[i] = Context::getSolverProperty<MFloat>("freqFactors", m_solverId, AT_, &freqFactor[i], i);

    /*! \page propertiesLS
      \section bodyMovementFunctions
      <code>MInt* LsCartesianSolver::m_static_computeBodyProperties_bodyToFunction</code>\n
      default = <code>1</code>\n \n
      Prescribes the functions for the body movement. Check the switch case in
      computeBodyProperties() for what each case actually does. \n \n
      Possible values are:
      <ul>
        <li>integer from 1 to 7</li>
      </ul>
      Keywords: <i>LEVELSET, BODY, BODY_MOTION, MOVING</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      bodyToFunction[i] =
          Context::getSolverProperty<MInt>("bodyMovementFunctions", m_solverId, AT_, &bodyToFunction[i], i);
    }
    Strouhal = Context::getSolverProperty<MFloat>("Strouhal", m_solverId, AT_, &Strouhal);

    /*! \page propertiesLS
      \section initialBodyCenters
      <code>MFloat LsCartesianSolver::initialBodyCenter</code>\n
      For each body, nDim Float values.
      With this property, one can move an STL around to its initial position.
      In initialInsidePoints is the "real" center of the stl files.
      This property, is the movement from the "real" center to the initial
      center used during calculation.
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      for(MInt j = 0; j < nDim; j++) {
        initialBodyCenter[i * nDim + j] = Context::getSolverProperty<MFloat>(
            "initialBodyCenters", m_solverId, AT_, &initialBodyCenter[i * nDim + j], i * nDim + j);
      }
    }

    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      for(MInt j = 0; j < nDim; j++) {
        normal[i * nDim + j] = Context::getSolverProperty<MFloat>("bodyMotionNormals", m_solverId, AT_,
                                                                  &normal[i * nDim + j], i * nDim + j);
      }
    }
    /*! \page propertiesLS
      \section liftStartAngles1
      <code>MFloat LsCartesianSolver::liftStartAngle1</code>\n
      default = <code>0.0</code>\n \n
      For each body, sets the start angle of the translation.
      The translation is described by a:
      - bodyToFunction case 1: cosine function
      - bodyToFunction case 2: valve lift shifted quadratic sine function (first angle)
      <ul>
      <li>Between 0 and 4.0 * PI </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      liftStartAngle1[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles1", m_solverId, AT_, &liftStartAngle1[i], i);
    }

    /*! \page propertiesLS
      \section liftStartAngles2
      <code>MFloat LsCartesianSolver::liftStartAngle2</code>\n
      default = <code>3.0 * PI</code>\n \n
      For each body, sets the start angle of the translation.
      The translation is described by a:
      - bodyToFunction case 2: valve lift shifted quadratic sine function (second angle)
      <ul>
      <li>Between 0 and 4.0 * PI </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      liftStartAngle2[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles2", m_solverId, AT_, &liftStartAngle2[i], i);
    }

    /*! \page propertiesLS
      \section liftEndAngles1
      <code>MFloat LsCartesianSolver::liftEndAngle1</code>\n
      default = <code>PI</code>\n \n
      For each body, sets the end angle of the translation.
      <ul>
      <li>Between 0 and 4.0 * PI </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      liftEndAngle1[i] = Context::getSolverProperty<MFloat>("liftEndAngles1", m_solverId, AT_, &liftEndAngle1[i], i);
    }

    /*! \page propertiesLS
      \section liftEndAngles2
      <code>MFloat LsCartesianSolver::liftEndAngle2</code>\n
      default = <code>4.0 * PI</code>\n \n
      For each body, sets the end angle of the translation.
      <ul>
      <li>Between 0 and 4.0 * PI </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      liftEndAngle2[i] = Context::getSolverProperty<MFloat>("liftEndAngles2", m_solverId, AT_, &liftEndAngle2[i], i);
    }

    /*! \page propertiesLS
      \section circleStartAngles
      <code>MFloat LsCartesianSolver::circleStartAngle</code>\n
      default = <code>0.0</code>\n \n
      For each body, sets the start angle for the circular motion case.
      (i.e. bodyToFunction case 5).
      <ul>
      <li>Between 0 and 4.0 * PI </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      circleStartAngle[i] =
          Context::getSolverProperty<MFloat>("circleStartAngles", m_solverId, AT_, &circleStartAngle[i], i);
    }

    /*! \page propertiesLS
      \section rotAngle
      <code> MFloat LsPar::rotAngle </code>  \n
      default = <code>"0.0"</code>\n \n
      Used to rotate the primary direction of the body movement.
      information in the property-file is expected in degree!
      <ul>
      <li>Any positive Float</li>
      </ul>
      Keywords: <i>LEVELSET EMBEDED BOUNDARY, MOVEMENT FUNCTIONS </i>
    */

    rotAngle = 0.0;
    rotAngle = Context::getSolverProperty<MFloat>("rotAngle", m_solverId, AT_, &rotAngle);
    rotAngle *= -PI / 180;

    // 3: compute relevant values:
    const MFloat freq0 = Strouhal * m_referenceVelocity;
    const MFloat freq02 = Strouhal;
    for(MInt k = 0; k < noEmbeddedBodies; k++) {
      // when using mu  : has a dimension!
      // when using mu2 : dimensionless!
      mu[k] = freqFactor[k] * freq0 * F2 * PI;
      mu2[k] = freqFactor[k] * freq02 * F2 * PI;
    }

    if(m_levelSetLb) {
      MFloat MaLb = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_);
      const MFloat referenceLengthLb = Context::getSolverProperty<MFloat>("referenceLengthLB", m_solverId, AT_);
      const MFloat freq02Lb = Strouhal * MaLb * LBCS / referenceLengthLb;
      for(MInt k = 0; k < noEmbeddedBodies; k++) {
        // when using mu2 : dimensionless!
        mu2[k] = freqFactor[k] * freq02Lb;
      }
    }

    // if bodyMovementFunction is 6 or 7, adjust start and end angles:
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      if(bodyToFunction[i] == 6 || bodyToFunction[i] == 7) {
        liftStartAngle1[i] = liftStartAngle1[i] * PI;
        liftEndAngle1[i] = liftEndAngle1[i] * PI - liftStartAngle1[i];
      }
    }

    omega = freqFactor[body] * m_referenceVelocity * PI / (F2 * amplitude[body]);

    if(Context::propertyExists("bodyTemperature", m_solverId)) {
      for(MInt i = 0; i < noEmbeddedBodies; i++) {
        temperature[i] = Context::getSolverProperty<MFloat>("bodyTemperature", m_solverId, AT_, &temperature[i], i);
      }
    }

    // read cotroll-point input data sizes
    MInt noControllFiles = 0;
    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      if(bodyToFunction[i] == 14) {
        noControllFiles++;
      }
    }

    if(noControllFiles > 0) {
      mAlloc(m_forcedMotionInput, noEmbeddedBodies, "m_forcedMotionInput", AT_);
      for(MInt i = 0; i < noEmbeddedBodies; i++) {
        m_forcedMotionInput[i].clear();
      }


      // read cotroll-point input data file
      for(MInt i = 0; i < noEmbeddedBodies; i++) {
        if(bodyToFunction[i] == 14) {
          ifstream readFile;
          MFloat coord = NAN;
          MFloat time = NAN;

          stringstream filename;
          filename.clear();
          filename.str("");
          filename << "forcedMotion_" << i << ".txt";

          readFile.open(filename.str(), ios_base::in);

          if(!readFile) {
            mTerm(1, AT_, "Error reading forced motion input file!");
          }
          for(;;) { // read the forced motion file

            if((readFile.rdstate() & ifstream::eofbit) != 0) {
              break;
            }

            readFile >> time;
            readFile >> coord;

            m_forcedMotionInput[i].insert(make_pair(time, coord));
          }
          readFile.close();
        }
      }
    }


    first = false;
  }

  if((m_levelSetMb && !m_levelSetLb) || m_levelSetFv) {
    if(returnMode == 4) {
      bodyData[0] = temperature[body];
      return;
    }

    //--------------------------------

    switch(bodyToFunction[body]) {
      case 1: // cosine function
      {
        angle = mu[body] * elapsedTime - liftStartAngle1[body];
        switch(returnMode) {
          case 1: // return body center
            if(angle > liftEndAngle1[body]) angle = liftEndAngle1[body];
            if(angle > 0) {
              bodyData[0] = -amplitude[body] * cos(angle);
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = -amplitude[body];
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 2: // return body velocity
            if((angle > 0 && angle < liftEndAngle1[body])) {
              bodyData[0] = mu[body] * amplitude[body] * sin(angle);
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if((angle > 0 && angle < liftEndAngle1[body])) {
              bodyData[0] = mu[body] * mu[body] * amplitude[body] * cos(angle);
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 2: // valve lift shifted quadratic sine
      {
        angle = mu2[body] * elapsedTime;

        switch(returnMode) {
          case 1: // return body center
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 2: // return body velocity
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                            * normal[body * nDim + 0];
              bodyData[1] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                            * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                              * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 0];
              bodyData[1] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 2];
              }
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }

        break;
      }
      case 3: // piston movement
      {
        angle = mu2[body] * elapsedTime;

        switch(returnMode) {
          case 1: // return body center
            if(elapsedTime > F0) {
              bodyData[0] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 0];
              bodyData[1] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = -amplitude[body] * normal[body * nDim + 0];
              bodyData[1] = -amplitude[body] * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = -amplitude[body] * normal[body * nDim + 2]; }
            }
            break;
          case 2: // return body velocity
            if(elapsedTime > F0) {
              bodyData[0] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 0];
              bodyData[1] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(elapsedTime > F0) {
              bodyData[0] =
                  mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 0];
              bodyData[1] =
                  mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] =
                    mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }

        break;
      }
      case 4: // cosine function with normal
      {
        angle = mu2[body] * elapsedTime;

        switch(returnMode) {
          case 1: // return body center
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * cos(angle) * normal[body * nDim + 2]; }
            } else {
              bodyData[0] = amplitude[body] * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * normal[body * nDim + 2]; }
            }
            break;
          case 2: // return body velocity
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
               || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
              bodyData[0] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 5: // circular motion:
      {
        // initialBodyCenter: center of rotation
        // amplitude        : radius
        // only 2d rotation implemented until now!
        angle = mu2[body] * elapsedTime + circleStartAngle[body];

        switch(returnMode) {
          case 1: // return body center
            if(elapsedTime > F0) {
              bodyData[0] = amplitude[body] * cos(angle);
              bodyData[1] = amplitude[body] * sin(angle);
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = amplitude[body] * cos(circleStartAngle[body]);
              bodyData[1] = amplitude[body] * sin(circleStartAngle[body]);
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 2: // return body velocity
            if(elapsedTime > F0) {
              bodyData[0] = -amplitude[body] * sin(angle);
              bodyData[1] = amplitude[body] * cos(angle);
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(elapsedTime > F0) {
              bodyData[0] = -amplitude[body] * cos(angle);
              bodyData[1] = -amplitude[body] * sin(angle);
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 6: // simplified piston motion
      {
        angle = omega * elapsedTime - liftStartAngle1[body];
        if(angle > liftEndAngle1[body]) angle = F0;

        switch(returnMode) {
          case 1: // return body center
            if(elapsedTime > F0) {
              bodyData[0] = -amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = -amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = -amplitude[body] * cos(angle) * normal[body * nDim + 2]; }
            } else {
              bodyData[0] = -amplitude[body] * normal[body * nDim + 0];
              bodyData[1] = -amplitude[body] * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = -amplitude[body] * normal[body * nDim + 2]; }
            }
            break;
          case 2: // return body velocity
            if(elapsedTime > F0) {
              bodyData[0] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 2]; }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(elapsedTime > F0) {
              bodyData[0] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }


        break;
      }
      case 7: // valve lift shifted quadratic sine matching piston motion 6
      {
        angle = omega * elapsedTime - liftStartAngle1[body];
        if(angle > liftEndAngle1[body]) angle = F0;

        switch(returnMode) {
          case 1: // return body center
            if(angle > F0) {
              bodyData[0] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 2]; }
            } else {
              bodyData[0] = 0;
              bodyData[1] = 0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
            }
            break;
          case 2: // return body velocity
            if(angle > F0) {
              bodyData[0] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(angle > F0) {
              bodyData[0] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 0];
              bodyData[1] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 0];
              bodyData[1] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 2]; }
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }

        break;
      }
      case 8: // translational movement in normal-direction with constant Ma-number (Ma_trans)
      {
        // Ma_trans is directly given by the amplitude-factor for each body!
        // The body velocity is then based on the free-stream speed of sound.
        switch(returnMode) {
          case 1: // return body center
            bodyData[0] = amplitude[body] * m_referenceVelocity * elapsedTime * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * m_referenceVelocity * elapsedTime * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) {
              bodyData[2] = amplitude[body] * m_referenceVelocity * elapsedTime * normal[body * nDim + 2];
            }
            break;
          case 2: // return body velocity
            bodyData[0] = amplitude[body] * m_referenceVelocity * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * m_referenceVelocity * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * m_referenceVelocity * normal[body * nDim + 2]; }
            break;
          case 3: // return body acceleration
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;

          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 9: // sine function with normel (=> periodic motion around the initialBodyCenter!)
      {
        angle = mu2[body] * elapsedTime;

        switch(returnMode) {
          case 1: // return body center
            if(angle >= liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] = amplitude[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * sin(angle) * normal[body * nDim + 2]; }
            } else if(angle < liftStartAngle1[body]) {
              bodyData[0] = 0;
              bodyData[1] = 0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
            } else if(angle > liftEndAngle1[body]) {
              bodyData[0] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 2];
              }
            }
            break;
          case 2: // return body velocity
            if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 10: // piston motion equation:
      {
        // reference coordincate system is the cylindeer head!
        // l            : rod length
        // r            : crank radius
        // TDC          : distance at TDC from the cylinder head
        // normal       : motion normal
        // phi          : crank angle in radian

        const MFloat l = liftStartAngle1[body];
        const MFloat TDC = liftEndAngle1[body];
        const MFloat r = amplitude[body];

        const MFloat phi = crankAngle(elapsedTime, 1);

        switch(returnMode) {
          case 1: // return body center
            bodyData[0] =
                normal[body * nDim + 0] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));
            bodyData[1] =
                normal[body * nDim + 1] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));
            IF_CONSTEXPR(nDim == 3) {
              bodyData[2] =
                  normal[body * nDim + 2] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));
            }

            if(domainId() == 0 && printPosition) {
              const MFloat cad = crankAngle(time(), 0);
              cerr << "Crank-angle-degree : " << cad << endl;
              if((cad / 45) > 0.989 && (cad / 45) < 1.01) {
                cerr << "Physical-Time : " << time() << endl;
                cerr << "                " << time() * 0.002898783653689 << endl;
              }
              cerr << "Piston-position    : "
                   << (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))) << " in m "
                   << bodyData[1] * 0.075 << endl;
              cerr << "Piston-velocity    : "
                   << mu2[body]
                          * (r * sin(phi)
                             + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))))
                   << endl;
            }

            break;
          case 2: // return body velocity
            bodyData[0] =
                normal[body * nDim + 0] * mu2[body]
                * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));

            bodyData[1] =
                normal[body * nDim + 1] * mu2[body]
                * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
            IF_CONSTEXPR(nDim == 3) {
              bodyData[2] =
                  normal[body * nDim + 2] * mu2[body]
                  * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
            }
            break;
          case 3: // return body acceleration
            bodyData[0] =
                normal[body * nDim + 0] * mu2[body] * mu2[body]
                * (r * cos(phi)
                   + (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                   + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / POW3((sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
            bodyData[1] =
                normal[body * nDim + 1] * mu2[body] * mu2[body]
                * (r * cos(phi)
                   + (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                   + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / (POW3(sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
            IF_CONSTEXPR(nDim == 3) {
              bodyData[2] =
                  normal[body * nDim + 2] * mu2[body] * mu2[body]
                  * (r * cos(phi)
                     + (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                     + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / (POW3(sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
        }
        break;
      }
      case 11: // valve motion for single-valve engine:
      {
        // reference coordincate system is the cylindeer head!
        // L            : maximum valve lift (2*amplitude[body])
        // normal       : motion normal
        // phi          : crank angle in radian
        // cad          : crank angle degree
        // maxNoCycles  : maximum number of considered cycles


        const MFloat cad = crankAngle(elapsedTime, 0);

        MFloat IO = liftStartAngle1[body];
        MFloat IC = liftEndAngle1[body];
        MFloat EO = liftStartAngle2[body];
        MFloat EC = liftEndAngle2[body];

        // possible values are:
        // IO: -3          // CAD BTDC 3
        // IC: 180+47      // CAD ABDC 47
        // EO: 3*180-47;   // CAD BBDC 47
        // EC: 720+3;      // CAD ATDC 3

        switch(returnMode) {
          case 1: // return body center
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= 0 && cad < IC) {
              MFloat phi_i = (cad - IO) * PI / 180;
              MFloat delta_i = (IC - IO) * PI / 180;
              MFloat freq_i = 1 / (delta_i / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));
              }

            } else if(cad > IC && cad < EO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < EC) {
              MFloat phi_e = (cad - EO) * PI / 180;
              MFloat delta_e = (EC - EO) * PI / 180;
              MFloat freq_e = 1 / (delta_e / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
              }

            } else if(cad >= EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            if(domainId() == 0 && printPosition) {
              cerr << "Crank-angle-degree : " << cad << endl;
              cerr << "Valve-position     : " << bodyData[0] << " in mm" << bodyData[0] * 75 << endl;
            }

            break;
          case 2: // return body velocity
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= 0 && cad < IC) {
              MFloat phi_i = (cad - IO) * PI / 180;
              MFloat delta_i = (IC - IO) * PI / 180;
              MFloat freq_i = 1 / (delta_i / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] =
                  -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);
              bodyData[1] =
                  -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] =
                    -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);
              }

            } else if(cad > IC && cad < EO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < EC) {
              MFloat phi_e = (cad - EO) * PI / 180;
              MFloat delta_e = (EC - EO) * PI / 180;
              MFloat freq_e = 1 / (delta_e / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] =
                  -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
              bodyData[1] =
                  -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] =
                    -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
              }
            } else if(cad >= EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            if(domainId() == 0 && printPosition) {
              cerr << "Valve-velocity     : " << bodyData[0] << endl;
            }

            break;
          case 3: // return body acceleration
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= 0 && cad < IC) {
              MFloat phi_i = (cad - IO) * PI / 180;
              MFloat delta_i = (IC - IO) * PI / 180;
              MFloat freq_i = 1 / (delta_i / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                            * sin(freq_i * phi_i - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                            * sin(freq_i * phi_i - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                              * sin(freq_i * phi_i - phase);
              }

            } else if(cad > IC && cad < EO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < EC) {
              MFloat phi_e = (cad - EO) * PI / 180;
              MFloat delta_e = (EC - EO) * PI / 180;
              MFloat freq_e = 1 / (delta_e / 2 / PI);
              MFloat phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                            * sin(freq_e * phi_e - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                            * sin(freq_e * phi_e - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                              * sin(freq_e * phi_e - phase);
              }
            } else if(cad > EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
        }
        break;
      }
      case 12: // valve motion for inlet valve:
      {
        // reference coordincate system is the cylindeer head!
        // amplitude    : dimensionless maximum  valve lift
        // normal       : normalised motion normal
        // phi          : crank angle in radian
        // cad          : crank angle degree
        // maxNoCycles  : maximum number of considered cycles

        const MFloat cad = crankAngle(elapsedTime, 0);

        MFloat scaleToMeter = 0.075;

        // valve-timing in crank-angle-degree

        MFloat IO = liftStartAngle1[body];
        MFloat IC = liftEndAngle1[body];

        // possible values are:
        // IO: -24         // 24 CAD BTDC
        // IC: 232         // CAD ATDC

        MFloat phi = 0;
        MFloat delta = -1;
        MFloat freq = 0;
        MFloat phase = -PI / 2;

        switch(returnMode) {
          case 1: // return body center
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= IO && cad <= IC) {
              phi = (cad - IO) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));
              }

            } else if(cad >= 720 + IO) {
              phi = (cad - (720 + IO)) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));
              }

            } else if(cad >= IC && cad < 720 + IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            if(domainId() == 0 && printPosition) {
              cerr << "Crank-angle-degree           : " << cad << endl;
              cerr << "dimensionless-Inlet-Valve-position : " << amplitude[body] * (1 - sin(freq * phi - phase))
                   << endl;
              cerr << "dimensional-Inlet-Valve-position   : "
                   << amplitude[body] * (1 - sin(freq * phi - phase)) * scaleToMeter << " in m " << endl;
              cerr << "dimensionless-Inlet-Valve-velocity : "
                   << -amplitude[body] * mu2[body] * freq * cos(freq * phi - phase) << endl;
            }

            break;
          case 2: // return body velocity
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= IO && cad <= IC) {
              phi = (cad - IO) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              }

            } else if(cad >= 720 + IO) {
              phi = (cad - (720 + IO)) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              }

            } else if(cad >= IC && cad < 720 + IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            break;
          case 3: // return body acceleration
            if(cad < IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= IO && cad <= IC) {
              phi = (cad - IO) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                              * sin(freq * phi - phase);
              }

            } else if(cad >= 720 + IO) {
              phi = (cad - (720 + IO)) * PI / 180;
              delta = (IC - IO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                              * sin(freq * phi - phase);
              }

            } else if(cad >= IC && cad < 720 + IO) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
        }
        break;
      }
      case 13: // valve motion for outlet valve:
      {
        // reference coordincate system is the cylindeer head!
        // amplitude    : dimensionless maximum  valve lift
        // normal       : motion normal
        // phi          : crank angle in radian
        // cad          : crank angle degree
        // maxNoCycles  : maximum number of considered cycles
        const MFloat cad = crankAngle(elapsedTime, 0);

        // valve-timing in crank-angle-degree

        MFloat EO = liftStartAngle1[body];
        MFloat EC = liftEndAngle1[body];
        MFloat scaleToMeter = 0.075;

        // possible values are:
        // EO: 480;    // CAD ATDC
        // EC: 16;     // CAD ATDC

        MFloat phi = 0;
        MFloat delta = -1;
        MFloat freq = 0;
        MFloat phase = -PI / 2;

        switch(returnMode) {
          case 1: // return body center
            if(cad < EO && cad >= EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < 720) {
              phi = (cad - EO) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));
              }

            } else if(cad < EC) {
              phi = -(EC - cad) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));
              }

            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            if(domainId() == 0 && printPosition) {
              cerr << "Crank-angle-degree           : " << cad << endl;
              cerr << "dimensionless-Outlet-Valve-position : " << amplitude[body] * (1 - sin(freq * phi - phase))
                   << endl;
              cerr << "dimensional-Outlet-Valve-position   : "
                   << amplitude[body] * (1 - sin(freq * phi - phase)) * scaleToMeter << " in m " << endl;
              cerr << "dimensionless-Outlet-Valve-velocity : "
                   << -amplitude[body] * mu2[body] * freq * cos(freq * phi - phase) << endl;
            }

            break;
          case 2: // return body velocity
            if(cad < EO && cad >= EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < 720) {
              phi = (cad - EO) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              }


            } else if(cad < EC) {
              phi = -(EC - cad) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
              }
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }

            break;
          case 3: // return body acceleration
            if(cad < EO && cad >= EC) {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            } else if(cad >= EO && cad < 720) {
              phi = (cad - EO) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                              * sin(freq * phi - phase);
              }

            } else if(cad < EC) {
              phi = -(EC - cad) * PI / 180;
              delta = (720 + EC - EO) * PI / 180;
              freq = 1 / (delta / 2 / PI);
              phase = -PI / 2;

              bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                            * sin(freq * phi - phase);
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body]
                              * (1 + POW2(freq) * sin(freq * phi - phase));
              }
            } else {
              mTerm(1, AT_, "Unexpected crank-angle-degree!");
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
        }
        break;
      }
      case 14: // liniear interpolation between controll-points
      {
        // for defined motion along motion normal and position(time) input from file
        // using linear interpolation to get the position and
        // left-sided-stencil to get velocity and acceleration
        // note: the stencil was choosen intuitively and any other might work as well.
        auto upperIt = upper_bound(m_forcedMotionInput[body].begin(), m_forcedMotionInput[body].end(),
                                   make_pair(elapsedTime, -999999.9));

        auto lowerIt = next(upperIt, -1);

        // linear interpolation between controll-points
        const MFloat pos = (*lowerIt).second
                           + (elapsedTime - (*lowerIt).first) * ((*upperIt).second - (*lowerIt).second)
                                 / ((*upperIt).first - (*lowerIt).first);


        if(domainId() == 0 && returnMode == 1 && printPosition) {
          cerr << globalTimeStep << " body " << body << " time " << elapsedTime << " earlier " << (*lowerIt).first
               << " later " << (*upperIt).first << " earlier Pos " << (*lowerIt).second << " later Pos "
               << (*upperIt).second << " interp. Pos " << pos << endl;
        }


        switch(returnMode) {
          case 1: // return body center
          {
            bodyData[0] = normal[body * nDim + 0] * pos;
            bodyData[1] = normal[body * nDim + 1] * pos;
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = normal[body * nDim + 2] * pos; }
            break;
          }
          case 2: // return body velocity
          {
            // get previous position for left-sided time gradient
            MFloat prevTime = elapsedTime - timeStep();
            if(globalTimeStep <= 0) prevTime = 0;
            upperIt = upper_bound(m_forcedMotionInput[body].begin(), m_forcedMotionInput[body].end(),
                                  make_pair(prevTime, -999999.9));
            lowerIt = next(upperIt, -1);

            const MFloat prevPos = (*lowerIt).second
                                   + (prevTime - (*lowerIt).first) * ((*upperIt).second - (*lowerIt).second)
                                         / ((*upperIt).first - (*lowerIt).first);
            MFloat vel = (pos - prevPos) / timeStep();
            if(globalTimeStep <= 0) vel = 0;

            bodyData[0] = normal[body * nDim + 0] * vel;
            bodyData[1] = normal[body * nDim + 1] * vel;
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = normal[body * nDim + 2] * vel; }

            break;
          }
          case 3: // return body acceleration
          {
            // get previous and older position for left-sided time gradient
            MFloat prevTime = elapsedTime - timeStep();
            if(globalTimeStep <= 0) prevTime = 0;


            upperIt = upper_bound(m_forcedMotionInput[body].begin(), m_forcedMotionInput[body].end(),
                                  make_pair(prevTime, -999999.9));
            lowerIt = next(upperIt, -1);

            const MFloat prevPos = (*lowerIt).second
                                   + (prevTime - (*lowerIt).first) * ((*upperIt).second - (*lowerIt).second)
                                         / ((*upperIt).first - (*lowerIt).first);


            MFloat olderTime = elapsedTime - 2 * timeStep();
            if(globalTimeStep <= 0) olderTime = 0;
            upperIt = upper_bound(m_forcedMotionInput[body].begin(), m_forcedMotionInput[body].end(),
                                  make_pair(olderTime, -999999.9));
            lowerIt = next(upperIt, -1);


            const MFloat olderPos = (*lowerIt).second
                                    + (olderTime - (*lowerIt).first) * ((*upperIt).second - (*lowerIt).second)
                                          / ((*upperIt).first - (*lowerIt).first);
            MFloat a = (pos - 2 * prevPos + olderPos) / (timeStep() * timeStep());
            if(globalTimeStep <= 0) a = 0;

            bodyData[0] = normal[body * nDim + 0] * a;
            bodyData[1] = normal[body * nDim + 1] * a;
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = normal[body * nDim + 2] * a; }

            break;
          }
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      default:
        mTerm(1, AT_, "function type not implemented. Please check bodyMovementFunctions property!");
    }

    // add the initialBodyCenter to the bodyDato for the body-positioning:
    if(returnMode == 1) {
      for(MInt dir = 0; dir < nDim; dir++) {
        bodyData[dir] += initialBodyCenter[body * nDim + dir];
      }
    }

    // rotate the final result around the z-axis
    if(rotAngle > 0 || rotAngle < 0) {
      MFloat tmp0 = bodyData[0] * cos(rotAngle) + bodyData[1] * sin(rotAngle);
      MFloat tmp1 = bodyData[1] * cos(rotAngle) - bodyData[0] * sin(rotAngle);
      bodyData[0] = tmp0;
      bodyData[1] = tmp1;
      bodyData[2] = bodyData[2];
    }
    if(m_periodicMovement) {
      MFloat bodyTempData = bodyData[m_periodicDirection];
      bodyData[m_periodicDirection] =
          std::fmod(bodyTempData + initialBodyCenter[body * nDim + m_periodicDirection] + 0.5 * m_periodicDistance,
                    m_periodicDistance)
          - 0.5 * m_periodicDistance - initialBodyCenter[body * nDim + m_periodicDirection];
    }
  } else if(m_levelSetLb) {
    switch(bodyToFunction[body]) {
      case 1: {
        // translational movement in normal-direction with constant velocity
        angle = c_cellLengthAtLevel(maxLevel());
        switch(returnMode) {
          case 1: // return body center
            bodyData[0] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 2]; }
            break;
          case 2: // return body velocity
            bodyData[0] = amplitude[body] * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * normal[body * nDim + 2]; }
            break;
          case 3: // return body acceleration
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;

          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 2: // sine function with normel (=> periodic motion around the initialBodyCenter!)
      {
        angle = freqFactor[body] * mu2[body] * elapsedTime + liftStartAngle1[body];
        switch(returnMode) {
          case 1: // return body center
            if(angle >= liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] =
                  amplitude[body] * c_cellLengthAtLevel(maxLevel()) / mu2[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] =
                  amplitude[body] * c_cellLengthAtLevel(maxLevel()) / mu2[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = amplitude[body] * c_cellLengthAtLevel(maxLevel()) / mu2[body] * sin(angle)
                              * normal[body * nDim + 2];
              }
            } else if(angle < liftStartAngle1[body]) {
              bodyData[0] = 0;
              bodyData[1] = 0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
            } else if(angle > liftEndAngle1[body]) {
              bodyData[0] = 0;
              bodyData[1] = 0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
            }
            break;
          case 2: // return body velocity
            if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] = amplitude[body] * cos(angle) * normal[body * nDim + 0];
              bodyData[1] = amplitude[body] * cos(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) { bodyData[2] = amplitude[body] * cos(angle) * normal[body * nDim + 2]; }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
              bodyData[0] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 0];
              bodyData[1] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 1];
              IF_CONSTEXPR(nDim == 3) {
                bodyData[2] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 2];
              }
            } else {
              bodyData[0] = F0;
              bodyData[1] = F0;
              IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      case 3: {
        if(!m_virtualSurgery) {
          TERMM(1, "This mode cannot be used without the m_virtualSurgery flag!");
        }
        switch(returnMode) {
          case 1: // return body center
          {
            // Blending between two level-sets using linear interpolation.
            // Only the 0th level set changes. The other sets are static.
            // Define the interpolation coefficients
            MIntScratchSpace startTime(m_noInterpolationRegions, AT_, "startTime");
            MIntScratchSpace noTimeSteps(m_noInterpolationRegions, AT_, "noTimeSteps");
            MIntScratchSpace timeStep(m_noInterpolationRegions, AT_, "timeStep");
            MFloatScratchSpace alpha(m_noInterpolationRegions, AT_, "alpha");
            for(MInt r = 0; r < m_noInterpolationRegions; r++) {
              if(r >= m_approxNoInterpReg) {
                startTime(r) = m_interpStartTime[m_approxNoInterpReg - 1];
                noTimeSteps(r) = m_noInterpTimeSteps[m_approxNoInterpReg - 1];
              } else {
                startTime(r) = m_interpStartTime[r];
                noTimeSteps(r) = m_noInterpTimeSteps[r];
              }
              timeStep(r) = mMax(0, globalTimeStep - startTime(r));
            }

            for(MInt r = 0; r < m_noInterpolationRegions; r++) {
              timeStep(r) = mMin(timeStep(r), noTimeSteps(r));
              alpha(r) = (F1 * timeStep(r)) / (F1 * noTimeSteps(r));
            }

            if(domainId() == 0) {
              cout << "GlobalTimeStep = " << globalTimeStep << ": " << endl;
              for(MInt r = 0; r < m_noInterpolationRegions; r++) {
                cout << "For region " << r << ": startTime = " << startTime(r) << ", timeStep = " << timeStep(r)
                     << ", alpha = " << alpha(r) << endl;
              }
            }

            // Recalculate the collected level set on the before collected cells
            for(MInt i = 0; i < a_noCells(); i++) {
              MInt cellId = i;
              MInt regions = mMax(0, m_cellIsInDiffRegion[cellId]);
              MFloat phi1 = m_correctedDistances[i][0];
              MFloat phi2 = m_correctedDistances[i][1];
              a_levelSetFunctionG(cellId, 0) = phi1 * (F1 - alpha(regions)) + phi2 * alpha(regions);
            }
            bodyData[0] = F0;
            bodyData[1] = F0;
            if(nDim == 3) {
              bodyData[2] = F0;
            }
            break;
          }
          case 2: // return body velocity
            bodyData[0] = F0;
            bodyData[1] = F0;
            if(nDim == 3) {
              bodyData[2] = F0;
            }
            break;
          case 3: // return body acceleration
            bodyData[0] = F0;
            bodyData[1] = F0;
            if(nDim == 3) {
              bodyData[2] = F0;
            }
            break;
          default:
            bodyData[0] = F0;
            bodyData[1] = F0;
            if(nDim == 3) bodyData[2] = F0;
            break;
        }
        break;
      }
      default:
        mTerm(1, AT_, "function type not implemented. Please check bodyMovementFunctions property!");
    }

    // add the initialBodyCenter to the bodyDato for the body-positioning:
    if(returnMode == 1) {
      for(MInt dir = 0; dir < nDim; dir++) {
        bodyData[dir] += initialBodyCenter[body * nDim + dir];
      }
    }

    // rotate the final result around the z-axis
    if(rotAngle > 0 || rotAngle < 0) {
      MFloat tmp0 = bodyData[0] * cos(rotAngle) + bodyData[1] * sin(rotAngle);
      MFloat tmp1 = bodyData[1] * cos(rotAngle) - bodyData[0] * sin(rotAngle);
      bodyData[0] = tmp0;
      bodyData[1] = tmp1;
      bodyData[2] = bodyData[2];
    }
    return;
  }
}


//----------------------------------------------------------------------------


/**
 * \fn void LsCartesianSolver::identifyBodies(MInt mode)
 * \brief sets a_bodyIdG(gCells,set) for all sets exept for the collected levelset
 *        (this is done in buildCollectedLevelSet())
 *
 * mode 0: during initialisation, after a balance, restart, adaptation (default)
 * mode 1: faster version which shifts the bodyId based on the levelSetShift
 *
 * \author Claudia Guenther, Update Tim Wegmann
 * \date   04/2011, 01/2019
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::identifyBodies(MInt mode) {
  TRACE();

  MBool& first = m_static_identifyBodies_first;
  MFloat(&initialInsidePoints)[s_maxNoEmbeddedBodies * 3] = m_static_identifyBodies_initialInsidePoints;

  MFloat& shiftTime = m_static_identifyBodies_shiftTime;

  // a) initialise and read initialInsidePoints-property
  //   (used as a starting Point to identify the different bodies)
  if(first) {
    if(m_noBodyBndryCndIds > s_maxNoEmbeddedBodies) {
      mTerm(1, AT_, "Error in LsCartesianSolver<nDim>::identifyBodies. Too many embedded Bodies!");
    }

    // set default values:
    for(MInt k = 0; k < s_maxNoEmbeddedBodies; k++) {
      for(MInt i = 0; i < nDim; i++) {
        initialInsidePoints[k * nDim + i] = F0;
      }
    }
    // read Properties

    /*! \page propertiesLS
      \section initialInsidePoints
      <code>MFloat LsCartesianSolver::initialInsidePoints[s_maxNoEmbeddedBodies*3]</code>\n
      default = <code>0.0 </code>\n\n
      Defines a reference point that is located inside the body for each embedded body (in x-,y-,z-coordinates)
      at the beginning of the calculation.\n
      Keywords: <i>LEVEL_SET</i>\n
    */
    for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
      for(MInt j = 0; j < nDim; j++) {
        initialInsidePoints[i * nDim + j] = Context::getSolverProperty<MFloat>(
            "initialInsidePoints", m_solverId, AT_, &initialInsidePoints[i * nDim + j], i * nDim + j);
      }
    }

    // CAUTION: during initSolver m_phsicalTime has not yet been set in the fvSolver!
    //         This leeds to difficulties when initialising a restart!
    //         This an update of the levelSet-Movement to the initial-Inside-Points is done here
    //         and reverted again at the end of the function!

    if(m_restart) {
      ASSERT(globalTimeStep == m_restartTimeStep,
             "globalTimeStep = " + std::to_string(globalTimeStep)
                 + "; m_restartTimeStep = " + std::to_string(m_restartTimeStep));
      for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
        for(MInt j = 0; j < nDim; j++) {
          initialInsidePoints[i * nDim + j] =
              initialInsidePoints[i * nDim + j] + m_semiLagrange_xShift_ref[j * m_noEmbeddedBodies + i];
        }
      }
    }

    shiftTime = 0;
  }

  MFloat elapsedTime = time();

  // b) initialise the bodyId in set 0 with -1 (if this is the collected-set)
  if(m_startSet > 0) {
    for(MInt gCellId = 0; gCellId < a_noCells(); gCellId++) {
      a_bodyIdG(gCellId, 0) = -1;
    }
  }

  MIntScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

  if(first || mode == 0 || !m_semiLagrange) {
    stack<MInt> bodyStack;


    // c) set the bodyId of other sets
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet[set] && !m_maxLevelChange) continue;
      if(!m_changedSet[set] && m_levelSetMb) continue;

      // 1) initialise with -1
      for(MInt gCellId = 0; gCellId < a_noCells(); gCellId++) {
        a_bodyIdG(gCellId, set) = -1;
      }
      // 2) go over all bodies in the set
      for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
        const MInt body = m_setToBodiesTable[set][b];
        MFloat bodyCenter[3] = {F0, F0, F0};
        computeBodyPropertiesForced(1, bodyCenter, body, time());
        for(MInt i = 0; i < nDim; i++) {
          bodyCenter[i] += initialInsidePoints[body * nDim + i];
        }
        // Consider azimuthalPer
        if(m_LsRotate) {
          std::array<MFloat, nDim> rotCenter{};
          std::array<MFloat, nDim> bodyCenter_shift{};
          std::array<MFloat, nDim> rotAngle{};
          std::copy_n(&bodyCenter[0], nDim, &bodyCenter_shift[0]);
          std::copy_n(&m_semiLagrange_xRot_STL[body * nDim], nDim, &rotAngle[0]);
          for(MInt d = 0; d < nDim; d++) {
            rotAngle[d] *= -F1;
          }
          computeBodyPropertiesForced(1, rotCenter.data(), body, 0.0);
          rotateLevelSet(1, bodyCenter, body, bodyCenter_shift.data(), rotCenter.data(), &rotAngle[0]);
          if(grid().azimuthalPeriodicity()) {
            MFloat center = grid().azimuthalCenter();
            MFloat angle = grid().azimuthalAngle();
            grid().raw().cartesianToCylindric(bodyCenter, bodyCenter_shift.data());
            MInt fac = 0;
            if(bodyCenter_shift[1] > (center + F1B2 * angle)) {
              fac = (MInt)((bodyCenter_shift[1] - (center - F1B2 * angle)) / angle);
            } else if(bodyCenter_shift[1] < (center - F1B2 * angle)) {
              fac = (MInt)((bodyCenter_shift[1] - (center + F1B2 * angle)) / angle);
            }
            grid().raw().rotateCartesianCoordinates(bodyCenter, fac * angle);
          }
        }

        const MInt startCellId = getContainingCell(bodyCenter);

        if(startCellId > -1 && startCellId < a_noCells()) {
          bodyStack.push(startCellId);
        }
        MInt cellsAdded = 1;

        while(cellsAdded) {
          cellsAdded = 0;

          while(!bodyStack.empty()) {
            const MInt currentCell = bodyStack.top();
            bodyStack.pop();
            if(a_bodyIdG(currentCell, set) > -1) continue;
            if(c_noChildren(currentCell) > 0 && !a_isHalo(currentCell)) {
              for(MInt child = 0; child < c_noChildren(currentCell); child++) {
                const MInt childId = c_childId(currentCell, child);
                if(childId < 0) continue;
                if(a_isHalo(childId)) continue;
                if(a_level(childId) == a_maxGCellLevel(set) && !a_isGZeroCell(childId, set)
                   && a_levelSetFunctionG(childId, set) * ((MFloat)m_levelSetSign[set]) >= 0)
                  continue;
                bodyStack.push(childId);
              }
              continue;
            }
            a_bodyIdG(currentCell, set) = body;
            cellsAdded++;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(!a_hasNeighbor(currentCell, dir)) continue;
              const MInt nghbrId = c_neighborId(currentCell, dir);
              if(a_isHalo(nghbrId)) continue;
              if(a_bodyIdG(nghbrId, set) > -1) continue;
              if(!a_isGZeroCell(nghbrId, set) && a_levelSetFunctionG(nghbrId, set) * ((MFloat)m_levelSetSign[set]) >= 0)
                continue;
              bodyStack.push(nghbrId);
            }
          }


          MPI_Allreduce(MPI_IN_PLACE, &cellsAdded, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "cellsAdded");

          // exchange halo-information
          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noWindowCells(i); j++) {
              tmp_data(windowCellId(i, j)) = a_bodyIdG(windowCellId(i, j), set);
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
              MInt windowId = grid().azimuthalWindowCell(i, j);
              tmp_data(windowId) = a_bodyIdG(windowId, set);
            }
          }

          exchangeDataLS(&tmp_data(0), 1);

          for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
            for(MInt j = 0; j < noHaloCells(i); j++) {
              const MInt haloCell = haloCellId(i, j);
              if(a_bodyIdG(haloCell, set) < 0) {
                a_bodyIdG(haloCell, set) = tmp_data(haloCell);
                if(a_bodyIdG(haloCell, set) > -1) {
                  for(MInt dir = 0; dir < m_noDirs; dir++) {
                    if(!a_hasNeighbor(haloCell, dir)) continue;
                    MInt nghbrId = c_neighborId(haloCell, dir);
                    if(a_isHalo(nghbrId)) continue;
                    if(a_bodyIdG(nghbrId, set) > -1) continue;
                    if(!a_isGZeroCell(nghbrId, set)
                       && a_levelSetFunctionG(nghbrId, set) * ((MFloat)m_levelSetSign[set]) >= 0)
                      continue;
                    bodyStack.push(nghbrId);
                  }
                }
              }
            }
          }
          for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
            for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
              const MInt haloCell = grid().azimuthalHaloCell(i, j);
              if(a_bodyIdG(haloCell, set) < 0) {
                a_bodyIdG(haloCell, set) = tmp_data(haloCell);
                if(a_bodyIdG(haloCell, set) > -1) {
                  for(MInt dir = 0; dir < m_noDirs; dir++) {
                    if(!a_hasNeighbor(haloCell, dir)) {
                      continue;
                    }
                    MInt nghbrId = c_neighborId(haloCell, dir);
                    if(a_isHalo(nghbrId)) {
                      continue;
                    }
                    if(a_bodyIdG(nghbrId, set) > -1) {
                      continue;
                    }
                    if(!a_isGZeroCell(nghbrId, set)
                       && a_levelSetFunctionG(nghbrId, set) * ((MFloat)m_levelSetSign[set]) >= 0) {
                      continue;
                    }
                    bodyStack.push(nghbrId);
                  }
                }
              }
            }
          }
        }
      } // loop over all bodies in the set

      // d) extend bodyIds further for anticipated bndry cells
      //    this is done at the same time for all bodies in the same set!

      MInt layerCount = 1;
      MInt cnt = 0;
      MIntScratchSpace lastLayer(a_noCells(), AT_, "lastLayer");

      // 1) extend the first layer, by going over all G0-cells in the set
      //    and find their neighbors!
      //    set the bodyId and add the cell to the lastLayer-list
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        const MInt cellId = a_G0CellId(id, set);
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          const MInt nghbrId = c_neighborId(cellId, dir);
          if(nghbrId == -1) continue;
          if(a_isGZeroCell(nghbrId, set)) continue;
          if(maxLevel() > minLevel()) {
            ASSERT(a_bodyIdG(nghbrId, set) == -1 || a_bodyIdG(nghbrId, set) == a_bodyIdG(cellId, set),
                   "Overlapping bodies in the same set not possible: gather bodies "
                       + to_string(a_bodyIdG(nghbrId, set)) + " " + to_string(a_bodyIdG(cellId, set))
                       + "in different sets! ");
          }
          lastLayer.p[cnt++] = nghbrId;
          a_bodyIdG(nghbrId, set) = a_bodyIdG(cellId, set);
        }
      }

      // 2) add the halo-neighbor-cells
      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noWindowCells(i); j++) {
          tmp_data(windowCellId(i, j)) = a_bodyIdG(windowCellId(i, j), set);
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
          MInt windowId = grid().azimuthalWindowCell(i, j);
          tmp_data(windowId) = a_bodyIdG(windowId, set);
        }
      }

      exchangeDataLS(&tmp_data(0), 1);

      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noHaloCells(i); j++) {
          const MInt haloCell = haloCellId(i, j);
          if(a_bodyIdG(haloCell, set) < 0) {
            a_bodyIdG(haloCell, set) = tmp_data(haloCell);
            if(a_bodyIdG(haloCell, set) > -1) {
              lastLayer.p[cnt++] = haloCell;
            }
          }
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
          const MInt haloCell = grid().azimuthalHaloCell(i, j);
          if(a_bodyIdG(haloCell, set) < 0) {
            a_bodyIdG(haloCell, set) = tmp_data(haloCell);
            if(a_bodyIdG(haloCell, set) > -1) {
              lastLayer.p[cnt++] = haloCell;
            }
          }
        }
      }

      MIntScratchSpace temp(a_noCells(), AT_, "temp");

      // 3) extend all other layers
      while(layerCount < m_gBandWidth + 1) {
        MInt tempCnt = 0;
        for(MInt c = 0; c < cnt; c++) {
          if(a_bodyIdG(lastLayer.p[c], set) == -1) continue;
          for(MInt dir = 0; dir < m_noDirs; dir++) {
            MInt nghbrId = c_neighborId(lastLayer.p[c], dir);
            if(nghbrId == -1) continue;
            if(a_isGZeroCell(nghbrId, set)) continue;
            if(a_bodyIdG(nghbrId, set) > -1) continue;
            temp.p[tempCnt++] = nghbrId;
            a_bodyIdG(nghbrId, set) = a_bodyIdG(lastLayer.p[c], set);
          }
        }

        // set bodyId for all halo-cells
        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_bodyIdG(windowCellId(i, j), set);
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
            MInt windowId = grid().azimuthalWindowCell(i, j);
            tmp_data(windowId) = a_bodyIdG(windowId, set);
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            if(a_bodyIdG(haloCell, set) < 0) {
              a_bodyIdG(haloCell, set) = tmp_data(haloCell);
              if(a_bodyIdG(haloCell, set) > -1) {
                temp.p[tempCnt++] = haloCell;
              }
            }
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
            const MInt haloCell = grid().azimuthalHaloCell(i, j);
            if(a_bodyIdG(haloCell, set) < 0) {
              a_bodyIdG(haloCell, set) = tmp_data(haloCell);
              if(a_bodyIdG(haloCell, set) > -1) {
                temp.p[tempCnt++] = haloCell;
              }
            }
          }
        }

        layerCount++;
        cnt = tempCnt;
        for(MInt c = 0; c < tempCnt; c++)
          lastLayer.p[c] = temp.p[c];
      }
    }

    // revert body-movement-changes
    if(m_restart && globalTimeStep == m_restartTimeStep && first) {
      for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
        for(MInt j = 0; j < nDim; j++) {
          initialInsidePoints[i * nDim + j] =
              initialInsidePoints[i * nDim + j] - m_semiLagrange_xShift_ref[j * m_noEmbeddedBodies + i];
        }
      }
    }

    shiftTime = elapsedTime;

  } else if(mode == 1) {
    // NOTE: faster version
    // however, this leads to differences towards the prevoius version,
    // especially for gapClosing and different bodies in the same set, that are close to each other!

    if(domainId() == 0) cerr << "Using fast identifyBodies-Version" << endl;
    MBool anyShift = false;

    // c) set the bodyId of other sets
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      if(!m_changedSet[set] & m_levelSetMb) continue;

      const MFloat length = 0.5 * c_cellLengthAtLevel(a_maxGCellLevel()); // / sqrt(nDim);
      // 1) go over all bodies in the set
      for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
        MInt body = m_setToBodiesTable[set][b];
        MFloat xCurrent[3] = {F0, F0, F0};
        MFloat xOld[3] = {F0, F0, F0};
        MFloat xShift[3] = {F0, F0, F0};

        computeBodyPropertiesForced(1, xCurrent, body, time() + timeStep());
        computeBodyPropertiesForced(1, xOld, body, shiftTime);

        // 2) find the body-Movement since the latest shift and check if a shift is necessary
        for(MInt d = 0; d < nDim; d++) {
          xShift[d] = xCurrent[d] - xOld[d];
        }
        MBool needShift = false;
        for(MInt d = 0; d < nDim; d++) {
          if(xShift[d] > length) {
            needShift = true;
            anyShift = true;
            shiftTime = elapsedTime;
          }
        }

        // 3) shift the bodyId
        if(needShift) {
          if(domainId() == 0) {
            cerr << "Shifting fast identifyBodies-Version " << shiftTime << " " << xShift[0] << " " << xShift[1] << " "
                 << xCurrent[0] << " " << xOld[0] << " " << xCurrent[1] << " " << xOld[1] << endl;
          }
          for(MInt dir = 0; dir < nDim; dir++) {
            if(xShift[dir] > length) {
              for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
                if(!c_isLeafCell(cellId)) continue;
                if(a_bodyIdG(cellId, set) != body) continue;

                if(!a_hasNeighbor(cellId, dir)) continue;
                MInt nghbrId = c_neighborId(cellId, dir);
                if(!a_bodyIdG(nghbrId, set)) {
                  a_bodyIdG(nghbrId, set) = body;
                }
              }
            }
          }
        }
      }
    }

    if(anyShift) {
      exchangeDataLS(&(a_bodyIdG(0, 0)), m_maxNoSets);
    }

  } else {
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      if(!m_changedSet[set] & m_levelSetMb) continue;

      // 1) initialise with -1 outside of the body:
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        if(a_levelSetFunctionG(cellId, set) * ((MFloat)m_levelSetSign[set]) > 0 && !a_isGZeroCell(cellId, set)) {
          a_bodyIdG(cellId, set) = -1;
        } else {
          if(a_level(cellId) == a_maxGCellLevel(set)) {
            ASSERT(a_bodyIdG(cellId, set) > -1, "");
          }
        }
      }

      MInt tempCnt = 0;
      MInt layerCount = 1;
      MInt cnt = 0;
      MIntScratchSpace lastLayer(a_noCells(), AT_, "lastLayer");
      MIntScratchSpace temp(a_noCells(), AT_, "temp");

      // 2) check G0-Cells and set bodyId for neighbors of the g0Cells

      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        const MInt cellId = a_G0CellId(id, set);
#if defined LS_DEBUG || !defined NDEBUG
        ASSERT(a_bodyIdG(cellId, set) > -1,
               "G0-Cell should have a valid bodyId " << to_string(c_globalId(cellId)) + " " + to_string(set));
        MBool bodyFound = false;
        for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
          MInt body = m_setToBodiesTable[set][b];
          if(a_bodyIdG(cellId, set) == body) {
            bodyFound = true;
            break;
          }
        }
        ASSERT(bodyFound, "BodyId of the G0-Cell is not matching any of the bodies in the set!");
#endif
        temp.p[id] = cellId;
        tempCnt++;
      }

      for(MInt c = 0; c < tempCnt; c++) {
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          MInt nghbrId = c_neighborId(temp.p[c], dir);
          if(nghbrId == -1) continue;
          if(a_isGZeroCell(nghbrId, set)) continue;
          lastLayer.p[cnt++] = nghbrId;
          a_bodyIdG(nghbrId, set) = a_bodyIdG(temp.p[c], set);
        }
      }

      // 2) add the halo-neighbor-cells
      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noWindowCells(i); j++) {
          tmp_data(windowCellId(i, j)) = a_bodyIdG(windowCellId(i, j), set);
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
          MInt windowId = grid().azimuthalWindowCell(i, j);
          tmp_data(windowId) = a_bodyIdG(windowId, set);
        }
      }

      exchangeDataLS(&tmp_data(0), 1);

      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        for(MInt j = 0; j < noHaloCells(i); j++) {
          const MInt haloCell = haloCellId(i, j);
          if(a_bodyIdG(haloCell, set) < 0) {
            a_bodyIdG(haloCell, set) = tmp_data(haloCell);
            if(a_bodyIdG(haloCell, set) > -1) {
              lastLayer.p[cnt++] = haloCell;
            }
          }
        }
      }
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
          const MInt haloCell = grid().azimuthalHaloCell(i, j);
          if(a_bodyIdG(haloCell, set) < 0) {
            a_bodyIdG(haloCell, set) = tmp_data(haloCell);
            if(a_bodyIdG(haloCell, set) > -1) {
              lastLayer.p[cnt++] = haloCell;
            }
          }
        }
      }

      // 3) extend all other layers
      while(layerCount < m_gBandWidth + 1) {
        tempCnt = 0;
        for(MInt c = 0; c < cnt; c++) {
          // if(a_bodyIdG(  lastLayer.p[c] , set)  == -1) continue;
          ASSERT(a_bodyIdG(lastLayer.p[c], set) > -1, "");
          for(MInt dir = 0; dir < m_noDirs; dir++) {
            MInt nghbrId = c_neighborId(lastLayer.p[c], dir);
            if(nghbrId == -1) continue;
            if(a_isGZeroCell(nghbrId, set)) continue;
            if(a_bodyIdG(nghbrId, set) > -1) continue;
            temp.p[tempCnt++] = nghbrId;
            a_bodyIdG(nghbrId, set) = a_bodyIdG(lastLayer.p[c], set);
          }
        }

        // set bodyId for all halo-cells
        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_bodyIdG(windowCellId(i, j), set);
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
            MInt windowId = grid().azimuthalWindowCell(i, j);
            tmp_data(windowId) = a_bodyIdG(windowId, set);
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            if(a_bodyIdG(haloCell, set) < 0) {
              a_bodyIdG(haloCell, set) = tmp_data(haloCell);
              if(a_bodyIdG(haloCell, set) > -1) {
                temp.p[tempCnt++] = haloCell;
              }
            }
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
            const MInt haloCell = grid().azimuthalHaloCell(i, j);
            if(a_bodyIdG(haloCell, set) < 0) {
              a_bodyIdG(haloCell, set) = tmp_data(haloCell);
              if(a_bodyIdG(haloCell, set) > -1) {
                temp.p[tempCnt++] = haloCell;
              }
            }
          }
        }

        layerCount++;
        cnt = tempCnt;
        for(MInt c = 0; c < tempCnt; c++)
          lastLayer.p[c] = temp.p[c];
      }
    }
  }


  first = false;
}


//----------------------------------------------------------------------------

/**
 * \brief returns cellId which contains the point (-1 if outside fluid domain)
 *
 * if point is really located in the fluid domain, cellId is returned
 * otherwise, -1 is returned
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::getContainingCell(MFloat* point) {
  TRACE();

  MInt startCellId = -1;

  // find proper starting cell on lowest level
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    if(a_level(cellId) != minLevel()) continue;
    if(inCell(cellId, point)) {
      startCellId = cellId;
      while(c_noChildren(startCellId) > 0) {
        MInt position = 0;
        if(point[0] > c_coordinate(startCellId, 0)) position += 1;
        if(point[1] > c_coordinate(startCellId, 1)) position += 2;
        IF_CONSTEXPR(nDim == 3) {
          if(point[2] > c_coordinate(startCellId, 2)) position += 4;
        }
        // If cell has children but does not have the child inside
        // which the point is located, the cell itself is returned
        const MInt childId = c_childId(startCellId, position);
        if(childId > -1 && !a_isHalo(childId)) {
          startCellId = childId;
        } else {
          return startCellId;
        }
      }
      return startCellId;
    }
  }

  ASSERT(startCellId < a_noCells(), "");

  if(startCellId > -1) {
    ASSERT(!a_isHalo(startCellId), "");
  }

  return startCellId;
}


//----------------------------------------------------------------------------

/**
 * \brief returns cellId which contains the point (-1 if outside fluid domain)
 *
 * if point is really located in the fluid domain, cellId is returned
 * otherwise, -1 is returned
 *
 * only cells neighboring the given startCell are investigated
 *
 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::getContainingCell(MInt startCell, MFloat* point, MInt set) {
  TRACE();

  MInt nghbrId = -1;

  MInt nghbrId2, nghbrId3;
  MInt d1, d2, d3, e, w, g;
  std::vector<MInt> diagCells;
  std::vector<MInt> diag2Cells;
  std::vector<MInt> signs;
  std::map<MInt, std::vector<MInt>> dirCode;

  //--------

  ASSERT(a_level(startCell) == a_maxGCellLevel(set) || a_level(startCell) == maxLevel(),
         "searchCell not on maxLevel! " + to_string(startCell));

  if(startCell < 0) return -1;
  if(inCell(startCell, point)) return startCell;
  // search neighbors
  for(MInt nghbr = 0; nghbr < 2 * nDim; nghbr++) {
    nghbrId = c_neighborId(startCell, nghbr);
    if(nghbrId == -1) continue;

    if(inCell(nghbrId, point)) return nghbrId;
  }

  signs.resize(nDim);

  // search diagonal neighbors
  MInt dir[3] = {0, 2, 4};
  MInt oDir[6] = {2, 4, 0, 4, 0, 2};
  for(MInt i = 0; i < nDim; i++) {
    for(MInt j = 0; j < 2; j++) {
      d1 = dir[i] + j;
      nghbrId = c_neighborId(startCell, d1);
      if(nghbrId == -1) continue;
      e = d1 / 2;
      signs[e] = d1 + pow(-1, j);
      for(MInt k = 0; k < 2; k++) {
        for(MInt m = 0; m < 2; m++) {
          d2 = oDir[2 * e + k] + m;
          nghbrId2 = c_neighborId(nghbrId, d2);
          if(nghbrId2 == -1) continue;
          diagCells.push_back(nghbrId2);
          w = d2 / 2;
          signs[w] = d2 + pow(-1, m);
          for(MInt n = 0; n < 2; n++) {
            d3 = dir[nDim - e - w] + n;
            if(a_hasNeighbor(nghbrId2, d3)) {
              g = d3 / 2;
              signs[g] = d3 + pow(-1, n);
              nghbrId3 = c_neighborId(nghbrId2, d3);
              diag2Cells.push_back(nghbrId3);
              if(dirCode.find(nghbrId3) == dirCode.end()) dirCode.insert(pair<MInt, vector<MInt>>(nghbrId3, signs));
            }
          }
        }
      }
    }
  }

  std::sort(diagCells.begin(), diagCells.end());
  auto last1 = std::unique(diagCells.begin(), diagCells.end());
  diagCells.erase(last1, diagCells.end());
  for(std::vector<MInt>::iterator it = diagCells.begin(); it != diagCells.end(); it++) {
    if(inCell(*it, point)) return *it;
  }
  std::sort(diag2Cells.begin(), diag2Cells.end());
  auto last2 = std::unique(diag2Cells.begin(), diag2Cells.end());
  diag2Cells.erase(last2, diag2Cells.end());
  for(std::vector<MInt>::iterator it = diag2Cells.begin(); it != diag2Cells.end(); it++) {
    if(inCell(*it, point)) return *it;
  }

  // If containing cell is not found search second layer
  // cerr << "D:" << domainId() << " SecondLayer " << startCell << " " << c_globalId(startCell) << " " <<
  // c_coordinate(startCell,0) << " " << c_coordinate(startCell,1) << " " << c_coordinate(startCell,2) << " " <<
  // a_level(startCell) << endl;

  nghbrId = checkSecondLayerCells(diag2Cells, dirCode, point);
  if(nghbrId > -1) {
    if(a_level(nghbrId) != a_maxGCellLevel() || (!m_reconstructOldG && m_initialGCell[nghbrId] == 0)) {
      nghbrId = -1;
    }
  }

  if(nghbrId == -1) {
    cerr << "D:" << domainId() << " Entire domain " << startCell << "/" << c_globalId(startCell) << " "
         << c_coordinate(startCell, 0) << " " << c_coordinate(startCell, 1) << " " << c_coordinate(startCell, 2) << " "
         << a_level(startCell) << " Point " << point[0] << " " << point[1] << " " << point[2] << endl;

    nghbrId = getContainingCell(point);

    if(nghbrId > -1
       && (a_level(nghbrId) != a_maxGCellLevel() || (!m_reconstructOldG && m_initialGCell[nghbrId] == 0))) {
      nghbrId = -1;
    }
    if(nghbrId > -1) {
      cerr << "Found " << domainId() << " " << nghbrId << "/" << c_globalId(nghbrId) << " " << c_coordinate(nghbrId, 0)
           << " " << c_coordinate(nghbrId, 1) << " " << c_coordinate(nghbrId, 2) << " " << a_level(nghbrId) << endl;
      MFloat dist = sqrt(pow(c_coordinate(nghbrId, 0) - c_coordinate(startCell, 0), 2.0)
                         + pow(c_coordinate(nghbrId, 1) - c_coordinate(startCell, 1), 2.0)
                         + pow(c_coordinate(nghbrId, 2) - c_coordinate(startCell, 2), 2.0))
                    / c_cellLengthAtLevel(a_maxGCellLevel());
      cerr << "Point " << point[0] << " " << point[1] << " " << point[2] << " " << dist << endl;
    }

    if(nghbrId == -1) {
      nghbrId = getContainingCellHalo(point);
      if(nghbrId > -1) {
        cerr << "D:" << domainId() << " Cell found in halo cell" << nghbrId;
      } else {
        mTerm(1, AT_, "No containing cell found... exiting!");
      }
    }
  }

  return nghbrId;


  /*
  MInt noDiagNghbrs = 8;
  IF_CONSTEXPR(nDim == 3)
    noDiagNghbrs = 24;
  MInt diagNghbrs2D[9] = {0,2,0,3,1,2,1,3,0};
  MInt diagNghbrs3D[25] = {0,2,0,3,1,2,1,3,0,4,0,5,1,5,2,4,3,4,1,4,2,5,3,5,0};
  MInt* diagNghbrs = diagNghbrs2D;
  IF_CONSTEXPR(nDim == 3)
    diagNghbrs = diagNghbrs3D;

  // check diagonal neighbors
  for( MInt d = 0; d < noDiagNghbrs; d++ ) {
    if( !a_hasNeighbor( startCell, diagNghbrs[d] ) )
      continue;
    nghbrId = c_neighborId( startCell, diagNghbrs[d] );
    if( !a_hasNeighbor( nghbrId, diagNghbrs[d+1] ) )
      continue;
    nghbrId = c_neighborId( nghbrId, diagNghbrs[d+1] );

    searchCells.push_back(nghbrId);
    IF_CONSTEXPR(nDim == 3){
      if( !a_hasNeighbor( nghbrId, diagNghbrs[d+2] ) )
        continue;
      nghbrId = c_neighborId( nghbrId, diagNghbrs[d+2] );
    }
    if(inCell(nghbrId, point)) {
      return nghbrId;
    }

 }



  cerr << "D:" << domainId() << "Containing " << nghbrId;
  if ( nghbrId > -1 ) {
    cerr << " " << a_level(nghbrId) << " " << c_globalId(nghbrId);

    for ( MInt i=0; i<nDim; i++ ) {
      cerr << (c_coordinate(nghbrId,i)-c_coordinate(startCell,i))/c_cellLengthAtLevel(a_level(nghbrId)) << " ";
      tmp += pow( (c_coordinate(nghbrId,i)-c_coordinate(startCell,i)) , 2.0 );
    }
    tmp = sqrt(tmp)/c_cellLengthAtLevel(a_level(nghbrId));
    cerr << tmp << endl;
  }



  return nghbrId;
  */
}

//----------------------------------------------------------------------------


/**
 * \brief returns true if the point is located inside the cell, otherwise false
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::inCell(MInt cellId, MFloat* point) {
  TRACE();

  MFloat xmin, xmax;
  MFloat ymin, ymax;
  MFloat zmin, zmax;
  MFloat cellHalfLength = c_cellLengthAtLevel(a_level(cellId) + 1);

  xmin = c_coordinate(cellId, 0) - cellHalfLength;
  xmax = c_coordinate(cellId, 0) + cellHalfLength;
  ymin = c_coordinate(cellId, 1) - cellHalfLength;
  ymax = c_coordinate(cellId, 1) + cellHalfLength;

  IF_CONSTEXPR(nDim == 2) {
    if(point[0] <= xmax && point[0] >= xmin && point[1] <= ymax && point[1] >= ymin) {
      return true;
    } else {
      return false;
    }
  }
  else {
    zmin = c_coordinate(cellId, 2) - cellHalfLength;
    zmax = c_coordinate(cellId, 2) + cellHalfLength;

    if(point[0] <= xmax && point[0] >= xmin && point[1] <= ymax && point[1] >= ymin && point[2] <= zmax
       && point[2] >= zmin) {
      return true;
    } else {
      return false;
    }
  }
}


//----------------------------------------------------------------------------


/**
 * \brief sets up interpolation stencil for levelSet interpolation
 * mode = 0 -> trilinear interpolation
 * cellId: cell on the highes refinement level which contains the point

 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::setUpLevelSetInterpolationStencil(const MInt cellId, MInt* interpolationCells,
                                                                MFloat* point) {
  TRACE();

  std::function<MBool(const MInt, const MInt)> alwaysTrue = [&](const MInt cell, const MInt dim) {
    return static_cast<MBool>(grid().tree().hasNeighbor(cell, dim));
  };

  return this->setUpInterpolationStencil(cellId, interpolationCells, point, alwaysTrue, false);
}


//----------------------------------------------------------------------------

/**
 * \brief sets up interpolation stencil for levelSet interpolation
 * mode = 0 -> trilinear interpolation
 * cellId: cell on the highes refinement level which contains the point
 *
 * Important: cell must not be a boundary cell! stencil must exist, since here no further checks are performed!
 *
 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setUpLevelSetInterpolationStencil(MInt cellId, MInt* interpolationCells, MInt position) {
  interpolationCells[position] = cellId;

  IF_CONSTEXPR(nDim == 2) {
    MInt nghbrX, nghbrY;
    MInt posIncrementX, posIncrementY;
    MInt xNghbrDir, yNghbrDir;
    if(position % 2 == 0) {
      xNghbrDir = 1;
      nghbrX = c_neighborId(cellId, xNghbrDir);
      posIncrementX = 1;
    } else {
      xNghbrDir = 0;
      nghbrX = c_neighborId(cellId, xNghbrDir);
      posIncrementX = -1;
    }
    if(((MInt)(position / 2)) % 2 == 0) {
      yNghbrDir = 3;
      nghbrY = c_neighborId(cellId, yNghbrDir);
      posIncrementY = 2;
    } else {
      yNghbrDir = 2;
      nghbrY = c_neighborId(cellId, yNghbrDir);
      posIncrementY = -2;
    }
    interpolationCells[position + posIncrementX] = nghbrX;
    interpolationCells[position + posIncrementY] = nghbrY;
    interpolationCells[position + posIncrementX + posIncrementY] = c_neighborId(nghbrX, yNghbrDir);
  }
  else {
    MInt nghbrX, nghbrY, nghbrZ;
    MInt posIncrementX, posIncrementY, posIncrementZ;
    MInt xNghbrDir, yNghbrDir, zNghbrDir;
    if(position % 2 == 0) {
      xNghbrDir = 1;
      nghbrX = c_neighborId(cellId, xNghbrDir);
      posIncrementX = 1;
    } else {
      xNghbrDir = 0;
      nghbrX = c_neighborId(cellId, xNghbrDir);
      posIncrementX = -1;
    }
    if(((MInt)(position / 2)) % 2 == 0) {
      yNghbrDir = 3;
      nghbrY = c_neighborId(cellId, yNghbrDir);
      posIncrementY = 2;
    } else {
      yNghbrDir = 2;
      nghbrY = c_neighborId(cellId, yNghbrDir);
      posIncrementY = -2;
    }
    if((MInt)(position / 4) == 0) {
      zNghbrDir = 5;
      nghbrZ = c_neighborId(cellId, zNghbrDir);
      posIncrementZ = 4;
    } else {
      zNghbrDir = 4;
      nghbrZ = c_neighborId(cellId, zNghbrDir);
      posIncrementZ = -4;
    }
    interpolationCells[position + posIncrementX] = nghbrX;
    interpolationCells[position + posIncrementY] = nghbrY;
    interpolationCells[position + posIncrementZ] = nghbrZ;
    interpolationCells[position + posIncrementX + posIncrementZ] = c_neighborId(nghbrX, zNghbrDir);
    interpolationCells[position + posIncrementX + posIncrementY] = c_neighborId(nghbrX, yNghbrDir);
    interpolationCells[position + posIncrementY + posIncrementZ] = c_neighborId(nghbrY, zNghbrDir);
    interpolationCells[position + posIncrementX + posIncrementY + posIncrementZ] =
        c_neighborId(c_neighborId(nghbrX, yNghbrDir), zNghbrDir);
  }
}


//----------------------------------------------------------------------------

/**
 * \brief shifts the old level set field 1 cell in prescribed direction

 * \author Claudia Guenther
 * \date 04/2012
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::shiftOldLevelSetField(MInt dir, MInt set, MInt body) {
  TRACE();

  const MInt otherDir[6] = {1, 0, 3, 2, 5, 4};
  const MInt bcId = m_bodyBndryCndIds[body];

  MFloatScratchSpace tmp_field(a_noCells(), AT_, "tmp_field");
  MFloat dx[3] = {F0, F0, F0};

  // shift STL to reference oldGField-position!
  if(m_GCtrlPntMethod == 2) {
    for(MInt i = 0; i < nDim; i++) {
      dx[i] = m_semiLagrange_xShift_ref[i * m_noEmbeddedBodies + body];
    }
    m_gCtrlPnt.CtrlPnt2_shiftSTL(bcId, dx);
  }

  // backup of the original old levelset-field
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    tmp_field[cellId] = a_oldLevelSetFunctionG(cellId, set);
  }

  // shift the field for non-halo cells with valid neighbors
  // for invalid neighbor the distance is recomputed from the stl!
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    if(a_level(cellId) != a_maxGCellLevel(set)) continue;
    if(a_bodyIdG(cellId, set) != body) continue;
    ASSERT(!a_isHalo(cellId), "");
    const MInt nghbrId = c_neighborId(cellId, otherDir[dir]);
    if(nghbrId == -1 || nghbrId == cellId) {
      MFloat refPoint[3] = {F0, F0, F0};
      for(MInt i = 0; i < nDim; i++) {
        refPoint[i] = c_coordinate(cellId, i);
      }
      MInt closestElement;
      MFloat closestPoint[3] = {F0, F0, F0};
      MFloat sign = determineLevelSetSignFromSTL(refPoint, set);
      if(m_levelSetSign[set] < 0) sign *= -1.0;
      const MFloat phi = sign * computeDistanceFromSTL(refPoint, &closestElement, closestPoint, set);
      a_oldLevelSetFunctionG(cellId, set) = phi;
    } else {
      a_oldLevelSetFunctionG(cellId, set) = tmp_field[nghbrId];
    }
  }

  // shift STL back to original position!
  if(m_GCtrlPntMethod == 2) {
    for(MInt i = 0; i < nDim; i++) {
      dx[i] = -m_semiLagrange_xShift_ref[i * m_noEmbeddedBodies + body];
    }
    m_gCtrlPnt.CtrlPnt2_shiftSTL(bcId, dx);
  }
}


//----------------------------------------------------------------------------

/**
 * \brief interpolates old levelset function to a given coordinate
 * mode = 0 -> trilinear interpolation
      and the one in the fv-solver!
 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::interpolateOldLevelSet(MInt* interpolationCells, MFloat* point, MInt referenceSet) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt set) {
    return static_cast<MFloat>(a_oldLevelSetFunctionG(cellId, set));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(c_coordinate(cellId, id));
  };

  return this->template interpolateFieldData<true>(&interpolationCells[0], &point[0], referenceSet, scalarField,
                                                   coordinate);
}

/**
 * \brief interpolates old levelset function to a given coordinate
 * mode = 0 -> trilinear interpolation
 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::interpolateLevelSet(MInt* interpolationCells, MFloat* point, MInt referenceSet) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt set) {
    return static_cast<MFloat>(a_levelSetFunctionG(cellId, set));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(c_coordinate(cellId, id));
  };

  return this->template interpolateFieldData<true>(&interpolationCells[0], &point[0], referenceSet, scalarField,
                                                   coordinate);
}


//----------------------------------------------------------------------------


/**
 * \brief builds a collected level-set function in set 0
 *
 * mode = 1: (default) builds collected level-set function, and bodyId in set 0 for all cells
 * mode = 0: (during initialisation) builds collected level-set function in set 0 for all cells
 * mode = 2: (for adaptation and balance) just as mode 1, but without the forward step in gapWidth
 *
 * variables that are set:
 * a_levelSetFunction of set 0
 * a_bodyIdG(of set 0
 * m_secondBodytifies the second closest body of each cell in the collected level-set function
 * a_gapWidth  -> unsymmetric sum of the distances of the closest and second closest
                  bodies (truncate cells with negative Ls-values!)
 *
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::buildCollectedLevelSet(MInt mode) {
  TRACE();

  // distance to the  closest and second closest bodies of each cell
  MFloatScratchSpace distBody1(a_noCells(), AT_, "distBody1");
  MFloatScratchSpace distBody2(a_noCells(), AT_, "distBody2");

  // sum of the distances to the closest and second closest bodies (gap width if cell is located
  // inside a gap)
  MFloatScratchSpace sumBodiesDist(a_noCells(), AT_, "sumBodiesDist");

  if(m_noSets < 2) return;

  const MFloat eps1 = c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim);
  const MFloat eps2 = c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * m_gapDeltaMin;

  if(mode == 1 || mode == 2) {
    // set a_levelSetFunctionG, a_bodyIdG

    // move forward in time
    if(mode == 1 && globalTimeStep > 0 && m_closeGaps) {
      for(MInt region = 0; region < m_noGapRegions; region++) {
        m_minGapWidthDt1[region] = m_minGapWidth[region];
        m_minGapWidth[region] = std::numeric_limits<MFloat>::max();
      }
    }

    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      MFloat currentWidth = std::numeric_limits<MFloat>::max();

      // a) determine body0 and phi0
      a_levelSetFunctionG(cellId, 0) = a_levelSetFunctionG(cellId, 1);
      a_bodyIdG(cellId, 0) = a_bodyIdG(cellId, 1);
      a_secondBodyId(cellId) = -1;
      for(MInt set = 2; set < m_noSets; set++) {
        MFloat phi0 = a_levelSetFunctionG(cellId, 0);
        MFloat phi1 = a_levelSetFunctionG(cellId, set);
        MInt body0 = a_bodyIdG(cellId, 0);
        MInt body1 = a_bodyIdG(cellId, set);
        MInt body2 = a_secondBodyId(cellId);
        MInt body3 = a_bodyIdG(cellId, set);

        //
        if(phi0 >= F0 && phi1 >= F0) {
          if(phi1 < phi0) {
            body2 = body0;
            body0 = body1;
          } else {
            // body0 = body0;
            if(body2 == -1 && body3 > -1) body2 = body3;
          }
          phi0 = mMin(phi0, phi1);
        } else if(phi0 <= F0 && phi1 <= F0) {
          if(abs(phi1) > abs(phi0)) {
            body2 = body0;
            body0 = body1;
          } else {
            if(body2 == -1 && body3 > -1) body2 = body3;
          }
          phi0 = -sqrt(phi0 * phi0 + phi1 * phi1);
        } else if(phi0 * phi1 <= F0) {
          if(phi0 < F0) {
            // phi0 = phi0;
            // body0 = body0;
            if(body2 == -1 && body3 > -1) body2 = body3;
          } else {
            phi0 = phi1;
            body2 = body0;
            body0 = body1;
          }
          ASSERT(phi0 <= phi1, "");
        } else {
          cerr << "WHAT THE FUCK! " << phi0 << phi1 << endl;
        }

        a_levelSetFunctionG(cellId, 0) = phi0;
        a_bodyIdG(cellId, 0) = body0;
        a_secondBodyId(cellId) = body2;
      }

      if(a_secondBodyId(cellId) == a_bodyIdG(cellId, 0)) {
        a_secondBodyId(cellId) = -1;
      }
      ASSERT(a_secondBodyId(cellId) < m_noBodyBndryCndIds, to_string(cellId) + " " + to_string(a_secondBodyId(cellId)));
      ASSERT(
          a_secondBodyId(cellId) < 0
              || (m_bodyToSetTable[a_secondBodyId(cellId)] > -1 && m_bodyToSetTable[a_secondBodyId(cellId)] < m_noSets),
          "");

      if(m_closeGaps) {
        distBody1[cellId] = a_levelSetFunctionG(cellId, 0);

        // b) set distBody2
        if(a_secondBodyId(cellId) > -1) {
          distBody2[cellId] = a_levelSetFunctionG(cellId, m_bodyToSetTable[a_secondBodyId(cellId)]);
        } else {
          distBody2[cellId] = m_outsideGValue;
        }

        // c) set sumBodiesDist
        if(a_secondBodyId(cellId) > -1 && a_bodyIdG(cellId, 0) > -1 && fabs(distBody1[cellId]) < m_outsideGValue
           && fabs(distBody2[cellId]) < m_outsideGValue) {
          if(m_gapInitMethod == 0) {
            sumBodiesDist[cellId] = distBody1[cellId] + distBody2[cellId];
          } else {
            sumBodiesDist[cellId] = fabs(distBody1[cellId]) + fabs(distBody2[cellId]);
          }


        } else {
          // default-value
          if(m_gapInitMethod == 0) {
            sumBodiesDist[cellId] = -m_outsideGValue;
          } else {
            sumBodiesDist[cellId] = -2 * m_outsideGValue;
          }
        }

        // d) set gapWidth
        if(m_gapInitMethod == 0) {
          // This unsymmetrical gap-behaviour concentrates gapCells in the gap-Area!
          // but is the main cause of gap-cell loos eventhough the gap is shrinking!

          if(a_levelSetFunctionG(cellId, 0) > F0) {
            a_gapWidth(cellId) = sumBodiesDist[cellId];
          } else if(fabs(a_levelSetFunctionG(cellId, 0)) < eps1) {
            a_gapWidth(cellId) = sumBodiesDist[cellId];
          } else {
            a_gapWidth(cellId) = -m_outsideGValue;
          }
        } else {
          a_gapWidth(cellId) = sumBodiesDist[cellId];
        }

        // Identifies if current cell is inside a gap
        const MBool isInsideGap = (fabs(a_gapWidth(cellId)) <= eps2 && a_potentialGapCell(cellId) > 0) ? true : false;

        // e) set m_minGapWidth
        if(m_gapInitMethod > 0 && mode == 1) {
          // CAUTION: temporarily overwritting sumBodiesDist
          if(a_secondBodyId(cellId) > -1 && a_bodyIdG(cellId, 0) > -1 && fabs(distBody1[cellId]) < m_outsideGValue
             && fabs(distBody2[cellId]) < m_outsideGValue) {
            sumBodiesDist[cellId] = distBody1[cellId] + distBody2[cellId];
          } else {
            sumBodiesDist[cellId] = 2 * m_outsideGValue;
          }

          if(a_potentialGapCell(cellId) && sumBodiesDist[cellId] <= eps2 && isInsideGap) {
            MInt regionId = a_potentialGapCellClose(cellId) - 2;
            ASSERT(regionId >= 0 && regionId < m_noGapRegions, to_string(regionId));

            currentWidth = distBody1[cellId] + distBody2[cellId];
            /*
            if(distBody1[ cellId ] > F0 && distBody2[ cellId ] > F0 ) {
              currentWidth = distBody1[ cellId ] + distBody2[ cellId ];
            } else if (distBody1[ cellId ] * distBody2[ cellId ] < F0) {
            //sign change in levseSet-value, however, this doesn't mean that the bodies are overlapping!

              currentWidth = distBody1[ cellId ] + distBody2[ cellId ];
              //currentWidth = mMin(distBody1[ cellId ], distBody2[ cellId ]);
            } else if (distBody1[ cellId ] < F0 && distBody2[ cellId ] < F0 ) {
              currentWidth = distBody1[ cellId ] + distBody2[ cellId ];
            }
            */
            if(currentWidth < m_minGapWidth[regionId]) {
              m_minGapWidth[regionId] = currentWidth;
            }
          }
        }
      }

    } // loop over all cells

    if(m_closeGaps && m_gapInitMethod > 0 && mode == 1) {
      MPI_Allreduce(MPI_IN_PLACE, &m_minGapWidth[0], m_noGapRegions, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_,
                    "MPI_IN_PLACE", "m_minGapWidth[0]");

      for(MInt region = 0; region < m_noGapRegions; region++) {
#if defined LS_DEBUG || !defined NDEBUG
        if(domainId() == 0 && m_minGapWidth[region] < 10 * m_outsideGValue) {
          cerr << "-->Min-Gap-Width at " << globalTimeStep << " for region " << region << " is "
               << m_minGapWidth[region] << " . Limit for Closure is: " << eps2 << endl;
        }
#endif
        if(domainId() == 0 && m_minGapWidth[region] < 0) {
          cerr << "Caution: overlapping levelSet-bodies in region " << region << " !" << endl;
        }

        if(m_minGapWidth[region] < -eps2) {
          MFloat dif = -(m_minGapWidth[region] + eps2);
          if(domainId() == 0) cerr << "Caution: temporary increase in gapDeltaMin, from " << m_gapDeltaMin;
          m_gapDeltaMin = m_gapDeltaMin + 1.5 * dif;
          if(domainId() == 0) cerr << " to " << m_gapDeltaMin << endl;
          ASSERT(m_noGapRegions == 1, "ERROR: Increased gapDeltaMin not yet implemented for multiple-Gap-Regions!");
          // for this m_gapDeltaMin needs to be region-dependand!
          // However this also leads to a variable shift of the g0-field for the gap-closing!
          // This should best be done body-dependand!
        } else {
          m_gapDeltaMin = m_gapDeltaMinOrig;
        }
      }
    }


  } else if(mode == 0) { // set a_levelSetFunctionG in set 0
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      a_levelSetFunctionG(cellId, 0) = a_levelSetFunctionG(cellId, 1);
      for(MInt set = 2; set < m_noSets; set++) {
        MFloat phi0 = a_levelSetFunctionG(cellId, 0);
        MFloat phi1 = a_levelSetFunctionG(cellId, set);

        //
        if(phi0 >= F0 && phi1 >= F0) {
          phi0 = mMin(phi0, phi1);
        } else if(phi0 <= F0 && phi1 <= F0) {
          phi0 = -sqrt(phi0 * phi0 + phi1 * phi1);
        } else if(phi0 * phi1 <= F0) {
          if(phi0 < F0) {
            // phi0 = phi0;
          } else {
            phi0 = phi1;
          }
        }

        a_levelSetFunctionG(cellId, 0) = phi0;
      }
    }
  } else {
    mTerm(1, AT_, "Unknown mode in buildCollectedLevelSet()");
  }
}

template <MInt nDim>
void LsCartesianSolver<nDim>::initializeCollectedLevelSet(MInt mode) {
  TRACE();

  grid().findEqualLevelNeighborsParDiagonal(false);

  //#######################################################################################################
  // Step 0: Initialize arrays and variables
  if(m_cellIsInDiffRegion != nullptr) {
    cout << "Deallocate m_cellIsInDiffRegion" << endl;
    mDeallocate(m_cellIsInDiffRegion);
  }
  if(m_correctedDistances != nullptr) {
    cout << "Deallocate m_correctedDistance" << endl;
    mDeallocate(m_correctedDistances);
  }

  mAlloc(m_cellIsInDiffRegion, a_noCells(), "m_cellIsInDiffRegion", -2, AT_);
  mAlloc(m_correctedDistances, a_noCells(), (m_maxNoSets - m_startSet), "m_correctedDistances", F0, AT_);
  MFloatScratchSpace boundingBox(m_noInterpolationRegions * m_noDirs, AT_, "boundingBox");
  MFloatScratchSpace globalBoundingBox(2 * nDim, AT_, "globalBoundingBox");
  MFloatScratchSpace pointCoord(nDim, AT_, "pointCoord");
  MFloatScratchSpace allPointCoords(noDomains(), nDim, AT_, "allPointCoords");
  MFloat eps = 1e-16;

  // How many cells are lying between the both static level sets
  m_geometry->getBoundingBox(globalBoundingBox.getPointer());

  // Find cells in which the two level sets differ (only on leaf cell level).
  // Those areas are marked with -1. All other cells are initialize with -2.
  // The parent of a marked cell is also marked with the same value. If the
  // children of a cell have different values, i.e. -1 and -2, -1 is set for
  // the parent.
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(!c_isLeafCell(cellId)) continue;
    MFloat phi1 = a_levelSetFunctionG(cellId, 1);
    MFloat phi2 = a_levelSetFunctionG(cellId, 2);
    MFloat deltaPhi = phi1 - phi2;
    if((fabs(deltaPhi) < eps)) {
      continue;
    } else {
      m_cellIsInDiffRegion[cellId] = -1;
      MBool reachedMinLevel = false;
      MInt parentId = cellId;
      while(!reachedMinLevel) {
        if(c_parentId(parentId) > -1) {
          parentId = c_parentId(parentId);
          if(m_cellIsInDiffRegion[parentId] < -1) {
            m_cellIsInDiffRegion[parentId] = -1;
          } else {
            reachedMinLevel = true;
          }
        } else {
          reachedMinLevel = true;
          ASSERT(a_level(parentId) == minLevel(), "This should not happen!");
        }
      }
    }
  }

  //#######################################################################################################
  // Step 1: Find all regions with a level set difference automatically using the region growing method
  m_noInterpolationRegions = 0;
  MBool foundAllRegions = false;
  while(!foundAllRegions) {
    // Step 1a: Reset the point coordinates. A value outside the bounding box is used
    for(MInt d = 0; d < nDim; d++) {
      pointCoord(d) = globalBoundingBox(d + nDim) + 10.0 * c_cellLengthAtLevel(minLevel());
    }

    // Step 1b: If there is still a cell lying in the "diff" region, set the
    //         point coordinates to the coordinates of the cell
    //         This is only done for leaf cells and non-halos
    for(MInt c = 0; c < grid().noLocalPartitionCells(); c++) {
      MInt cellId = grid().localPartitionCellLocalIds(c);
      if(cellId < 0) continue;

      if(m_cellIsInDiffRegion[cellId] == -1) {
        for(MInt d = 0; d < nDim; d++) {
          pointCoord(d) = c_coordinate(cellId, d);
          cout << "Writing point for cell for domain " << domainId() << " and region " << m_noInterpolationRegions
               << ": " << pointCoord(d) << endl;
        }
        break;
      }
    }

    // Step 1c: The point coordinates are exchanged
    for(MInt d = 0; d < nDim; d++) {
      allPointCoords(domainId(), d) = pointCoord(d);
    }
    MPI_Allgather(MPI_IN_PLACE, nDim, MPI_DOUBLE, &allPointCoords[0], nDim, MPI_DOUBLE, mpiComm(), AT_, "MPI_IN_PLACE",
                  "allPointCoords");

    // Step 1d: Check for all domains if there is a point which is not outside
    //         of the bounding box (which is the initialized value)
    //         If there is a cell on one process, this process starts the region
    //         growing. If not there are no further regions and we can stop here
    foundAllRegions = true;
    for(MInt i = 0; i < noDomains(); i++) {
      MBool startPoint = true;
      for(MInt d = 0; d < nDim; d++) {
        if(allPointCoords(i, d) > globalBoundingBox(d + nDim)) startPoint = false;
      }
      if(startPoint) {
        for(MInt d = 0; d < nDim; d++) {
          pointCoord(d) = allPointCoords(i, d);
        }
        foundAllRegions = false;
        break;
      }
    }
    if(foundAllRegions) {
      break;
    }

    // Step 1e: Find a start Cell in a not yet concidered modified region
    MInt startCell = -1;
    for(MInt i = 0; i < a_noCells(); i++) {
      if(a_isHalo(i)) continue;
      if(c_parentId(i) > -1) continue;
      if(m_cellIsInDiffRegion[i] != -1) continue;

      MFloat halfLength = grid().cellLengthAtLevel(a_level(i)) * F1B2;
      MInt cnt = 0;

      for(MInt d = 0; d < nDim; d++) {
        if(abs(c_coordinate(i, d) - pointCoord(d)) <= halfLength) cnt++;
      }
      if(cnt == nDim) {
        startCell = i;
        cout << "START Cell Found on process " << domainId() << endl;
        ASSERT(a_level(startCell) == minLevel(), "This should not happen!!!");
        break;
      }
    }

    // Step 1f: Use the region growing to find all connected cells of the start cell
    //         Thus all cells belonging to the region can be found
    if(startCell > -1) {
      regionGrowing(startCell, m_noInterpolationRegions);
    }

    // Step 1g: In case of parallel simulations, the region growing has to be communicated
    if(noDomains() > 1) {
      // prepare communication
      MIntScratchSpace noSendWindowPerDomain(noNeighborDomains(), AT_, "noSendWindowPerDomain");
      MIntScratchSpace noReceiveHaloPerDomain(noNeighborDomains(), AT_, "noReceiveHaloPerDomain");
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        noSendWindowPerDomain[d] = noWindowCells(d);
        ASSERT(noSendWindowPerDomain[d] >= 0, "noSendWindowPerDomain[d] < 0");
      }
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        noReceiveHaloPerDomain[d] = noHaloCells(d);
        ASSERT(noReceiveHaloPerDomain[d] >= 0, "noReceiveHaloPerDomain[d] < 0");
      }

      MInt allSend = 0;
      MInt allReceive = 0;
      vector<MInt> offsetsSend;
      vector<MInt> offsetsReceive;
      for(MInt dom = 0; dom < noNeighborDomains(); dom++) {
        offsetsSend.push_back(allSend);
        offsetsReceive.push_back(allReceive);
        allSend += noSendWindowPerDomain[dom];
        allReceive += noReceiveHaloPerDomain[dom];
      }

      MInt noChanges = 1;
      MInt noIterations = 0;
      while(noChanges != 0) {
        noChanges = 0;

        MIntScratchSpace sndBufWin(allSend, AT_, "sndBufWin");
        MIntScratchSpace rcvBufHalo(allReceive, AT_, "rcvBufHalo");

        MPI_Request* mpi_request_;
        mAlloc(mpi_request_, noNeighborDomains(), "mpi_request_", AT_);

        for(MInt d = 0; d < noNeighborDomains(); d++) {
          if(noSendWindowPerDomain[d] > 0) {
            for(MInt c = 0; c < noSendWindowPerDomain[d]; c++) {
              sndBufWin[offsetsSend[d] + c] = m_cellIsInDiffRegion[windowCellId(d, c)];
            }
            MPI_Issend(&(sndBufWin[offsetsSend[d]]), noSendWindowPerDomain[d], MPI_INT, neighborDomain(d), 0, mpiComm(),
                       &mpi_request_[d], AT_, "(sndBufWin[offsetsSend[d]])");
          }
        }
        MPI_Status status_;
        for(MInt d = 0; d < noNeighborDomains(); d++) {
          if(noReceiveHaloPerDomain[d] > 0) {
            MPI_Recv(&(rcvBufHalo[offsetsReceive[d]]), noReceiveHaloPerDomain[d], MPI_INT, neighborDomain(d), 0,
                     mpiComm(), &status_, AT_, "(rcvBufHalo[offsetsReceive[d]])");
          }
        }

        for(MInt d = 0; d < noNeighborDomains(); d++) {
          if(noReceiveHaloPerDomain[d] > 0 && noSendWindowPerDomain[d] > 0) {
            MPI_Wait(&mpi_request_[d], &status_, AT_);
          }
        }

        for(MInt d = 0; d < noNeighborDomains(); d++) {
          if(noReceiveHaloPerDomain[d] > 0) {
            for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
              const MInt halo = haloCellId(d, c);
              if(c_parentId(halo) > -1) continue;
              ASSERT(a_level(halo) == minLevel(), "This should not happen!!!!");
              if((rcvBufHalo[c + offsetsReceive[d]] == m_noInterpolationRegions)
                 && (m_cellIsInDiffRegion[halo] == -1)) {
                regionGrowing(halo, m_noInterpolationRegions);
                noChanges++;
              }
            }
          }
        }
        // Are there still changes in the process, if yes, start all over again with the communication
        MPI_Allreduce(MPI_IN_PLACE, &noChanges, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "noChanges");
        mDeallocate(mpi_request_);
        noIterations++;
      }
      cout << "NO CHANGES " << noIterations << " " << m_noInterpolationRegions << endl;
    }
    m_noInterpolationRegions++;
  }

  cout << "m_noInterpolationRegions " << m_noInterpolationRegions << endl;
  // Find cells in which the two level sets differ (only on leaf cell level).
  // Those areas are marked with -1. All other cells are initialize with -2.
  // The parent of a marked cell is also marked with the same value. If the
  // children of a cell have different values, i.e. -1 and -2, -1 is set for
  // the parent.
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(c_parentId(cellId) > -1) continue;
    ASSERT(m_cellIsInDiffRegion[cellId] != -1, "It seems like some minLevel cells were skipped!");
    if(m_cellIsInDiffRegion[cellId] < 0) continue;
    for(MInt c = 0; c < c_noChildren(cellId); c++) {
      if(c_childId(cellId, c) < 0) continue;
      setChildRegions(c_childId(cellId, c), m_cellIsInDiffRegion[cellId]);
    }
  }

  if(noDomains() > 1) {
    // prepare communication
    MIntScratchSpace noSendWindowPerDomain(noNeighborDomains(), AT_, "noSendWindowPerDomain");
    MIntScratchSpace noReceiveHaloPerDomain(noNeighborDomains(), AT_, "noReceiveHaloPerDomain");
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      noSendWindowPerDomain[d] = noWindowCells(d);
      ASSERT(noSendWindowPerDomain[d] >= 0, "noSendWindowPerDomain[d] < 0");
    }
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      noReceiveHaloPerDomain[d] = noHaloCells(d);
      ASSERT(noReceiveHaloPerDomain[d] >= 0, "noReceiveHaloPerDomain[d] < 0");
    }

    MInt allSend = 0;
    MInt allReceive = 0;
    vector<MInt> offsetsSend;
    vector<MInt> offsetsReceive;
    for(MInt dom = 0; dom < noNeighborDomains(); dom++) {
      offsetsSend.push_back(allSend);
      offsetsReceive.push_back(allReceive);
      allSend += noSendWindowPerDomain[dom];
      allReceive += noReceiveHaloPerDomain[dom];
    }

    MIntScratchSpace sndBufWin(allSend, AT_, "sndBufWin");
    MIntScratchSpace rcvBufHalo(allReceive, AT_, "rcvBufHalo");

    MPI_Request* mpi_request_;
    mAlloc(mpi_request_, noNeighborDomains(), "mpi_request_", AT_);

    for(MInt d = 0; d < noNeighborDomains(); d++) {
      if(noSendWindowPerDomain[d] > 0) {
        for(MInt c = 0; c < noSendWindowPerDomain[d]; c++) {
          sndBufWin[offsetsSend[d] + c] = m_cellIsInDiffRegion[windowCellId(d, c)];
        }
        MPI_Issend(&(sndBufWin[offsetsSend[d]]), noSendWindowPerDomain[d], MPI_INT, neighborDomain(d), 0, mpiComm(),
                   &mpi_request_[d], AT_, "(sndBufWin[offsetsSend[d]])");
      }
    }
    MPI_Status status_;
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      if(noReceiveHaloPerDomain[d] > 0) {
        MPI_Recv(&(rcvBufHalo[offsetsReceive[d]]), noReceiveHaloPerDomain[d], MPI_INT, neighborDomain(d), 0, mpiComm(),
                 &status_, AT_, "(rcvBufHalo[offsetsReceive[d]])");
      }
    }

    for(MInt d = 0; d < noNeighborDomains(); d++) {
      if(noReceiveHaloPerDomain[d] > 0 && noSendWindowPerDomain[d] > 0) {
        MPI_Wait(&mpi_request_[d], &status_, AT_);
      }
    }

    for(MInt d = 0; d < noNeighborDomains(); d++) {
      if(noReceiveHaloPerDomain[d] > 0) {
        for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
          const MInt halo = haloCellId(d, c);
          m_cellIsInDiffRegion[halo] = rcvBufHalo[c + offsetsReceive[d]];
        }
      }
    }
  }
  //#######################################################################################################
  // Step 2: Initialize the array containing the coorectedDistances in case the band
  // width is not set high enough

  if(mode == 0) {
    m_refinedCells.clear();
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(!c_isLeafCell(cellId)) continue;
      if(!a_isHalo(cellId)) {
        if(m_cellIsInDiffRegion[cellId] >= 0) {
          m_refinedCells.insert(make_pair(cellId, 0));
        }
      }
    }
  }

  MInt noRefinedBandCells = (signed)m_refinedCells.size();
  MPI_Allreduce(MPI_IN_PLACE, &noRefinedBandCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                "noRefinedBandCells");

  // On the fly reconstruction of levelset band cells outside the G0 Set!
  cerr << "Reinitialisation of " << noRefinedBandCells << " cells " << endl;
  if(noRefinedBandCells > 0) {
    m_startSet = 1;
    // construct old-LevelSet data based on shifted geometry!
    constructGFieldFromSTL(4);

    exchangeDataLS(&(a_oldLevelSetFunctionG(0, 0)), m_maxNoSets);

    // copy old levelset to levelset for stationary bodies
    for(auto& m_refinedCell : m_refinedCells) {
      const MInt cellId = m_refinedCell.first;
      ASSERT(!a_isHalo(cellId), "");
      for(MInt setI = m_startSet; setI < m_noSets; setI++) {
        a_levelSetFunctionG(cellId, setI) = a_oldLevelSetFunctionG(cellId, setI);
      }
    }

    for(MInt lvl = maxLevel() - 1; lvl >= minLevel(); lvl--) {
      for(MInt cell = 0; cell < a_noCells(); cell++) {
        for(MInt set = 1; set < m_maxNoSets; set++) {
          if(a_level(cell) != lvl) continue;
          MFloat levelSetCoarse = F0;
          if(c_noChildren(cell) > 0) {
            for(MInt child = 0; child < IPOW2(nDim); child++) {
              if(c_childId(cell, child) < 0) continue;
              levelSetCoarse += a_levelSetFunctionG(c_childId(cell, child), set);
            }
            a_levelSetFunctionG(cell, set) = levelSetCoarse / (MFloat)c_noChildren(cell);
          }
        }
      }
    }

    exchangeLevelSet();

    m_refinedCells.clear();

    m_startSet = 0;
  }

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    MFloat phi1 = a_levelSetFunctionG(cellId, 1);
    MFloat phi2 = a_levelSetFunctionG(cellId, 2);
    m_correctedDistances[cellId][0] = phi1;
    m_correctedDistances[cellId][1] = phi2;
  }

  MFloat xCurrent[3] = {F0, F0, F0};
  for(MInt set = m_startSet; set < m_noSets; set++) {
    if(!m_computeSet_backup[set]) continue;
    for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
      const MInt body = m_setToBodiesTable[set][b];
      computeBodyPropertiesForced(1, xCurrent, body, time() + timeStep());
    }
  }
}

template <MInt nDim>
void LsCartesianSolver<nDim>::regionGrowing(MInt cellId, MInt region) {
  TRACE();

  if(m_cellIsInDiffRegion[cellId] == -1) {
    m_cellIsInDiffRegion[cellId] = region;

    for(MInt n = 0; n < IPOW3[nDim] - 1; n++) {
      if(grid().neighborList(cellId, n) < 0) continue;
      if(m_cellIsInDiffRegion[grid().neighborList(cellId, n)] != -1) continue;
      regionGrowing(grid().neighborList(cellId, n), region);
    }
  }
}

template <MInt nDim>
void LsCartesianSolver<nDim>::setChildRegions(MInt cellId, MInt region) {
  TRACE();

  if(m_cellIsInDiffRegion[cellId] == -1) {
    m_cellIsInDiffRegion[cellId] = region;

    for(MInt c = 0; c < c_noChildren(cellId); c++) {
      if(c_childId(cellId, c) < 0) continue;
      setChildRegions(c_childId(cellId, c), region);
    }
  }
}

//----------------------------------------------------------------------------


/**
 * \brief builds a collected level-set function in set 0 - no secondary variables are computed
 *
 * mode = 0: builds collected level-set function in band of set 0
 * mode = 1: builds collected level-set function in all cells of set 0
 *
 * near gap cells are skipped!
 *
 * variables that are set:
 * a_levelSetFunction of set 0
 *
 * Note: the commented lines below (self assignments) remain for documentation purposes.
 *
 * \author Claudia Guenther
 * \date 05/2011
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::reBuildCollectedLevelSet(MInt mode) {
  TRACE();

  if(m_noSets < 2) {
    return;
  }
  ASSERT(m_buildCollectedLevelSetFunction, "");

  MInt noCells = -1;
  if(mode)
    noCells = a_noCells();
  else
    noCells = a_noBandCells(0);

  for(MInt i = 0; i < noCells; i++) {
    MInt cellId = -1;
    if(mode) {
      cellId = i;
    } else {
      cellId = a_bandCellId(i, 0);
    }
    if(a_nearGapG(cellId) > 0 && a_potentialGapCellClose(cellId)) continue;
    a_levelSetFunctionG(cellId, 0) = a_levelSetFunctionG(cellId, 1);

    for(MInt set = 2; set < m_noSets; set++) {
      MFloat phi0 = a_levelSetFunctionG(cellId, 0);
      MFloat phi1 = a_levelSetFunctionG(cellId, set);

      if(phi0 >= F0 && phi1 >= F0) {
        phi0 = mMin(phi0, phi1);
      } else if(phi0 <= F0 && phi1 <= F0) {
        phi0 = -sqrt(phi0 * phi0 + phi1 * phi1);
      } else if(phi0 * phi1 <= F0) {
        if(phi0 < F0) {
          // phi0 = phi0;
        } else {
          phi0 = phi1;
        }
      }

      a_levelSetFunctionG(cellId, 0) = phi0;
    }
  }
}


//----------------------------------------------------------------------------

/**
 * \brief identifies cells located in non-resolvable gaps
 *
 * the following variables are set:
 * a_nearGapG in the lsSolver
 * a_wasGapCell and a_isGapCell in the fvSolver!
 *
 * \author Claudia Guenther, Update: Tim Wegmann
 * \date   04/2011, 01/2019
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::gapHandling() {
  TRACE();

  if(!m_closeGaps) return;

  MIntScratchSpace gapCells(a_noCells(), AT_, "gapCells");
  MIntScratchSpace lastLayer(a_noCells(), AT_, "lastLayer");
  MBoolScratchSpace isInsideGap(a_noCells(), AT_, "isInsideGap");
  MBoolScratchSpace isInsideGapTmp(a_noCells(), AT_, "isInsideGapTmp");
  stack<MInt> gapStack;

  // same as eps2 used in buildCollectedLevelSet()
  const MFloat eps = c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * m_gapDeltaMin;

  const MFloat eps2 = c_cellLengthAtLevel(a_maxGCellLevel(0)) * (5.0 / 4.0 * sqrt(nDim)) * m_gapDeltaMin;
  const MFloat eps3 = (m_gapInitMethod < 2)
                          ? c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * (1.5 + m_gapDeltaMin)
                          : c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * (2.5 + m_gapDeltaMin);

  MInt maxNoLayers = mMax(6, m_gBandWidth);

  // MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  // MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");

  MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

  //-------------------------

  // b) reset gap-statistics
  m_noOldGapCells = m_noGapCells;
  m_noGapCells = 0;
  MInt gapCellCounter = 0;
  MInt cellsAdded = 1;

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    isInsideGap(cellId) = false;
    // Identifies cells inside a gap
    if((fabs(a_gapWidth(cellId)) <= eps && a_potentialGapCell(cellId) > 0)
       && a_inBandG(cellId, 0)) { // only if cell is inside the G0-Band
      isInsideGap(cellId) = true;
    }
    isInsideGapTmp(cellId) = false;
    a_nearGapG(cellId) = false;
  }

  MBoolScratchSpace forceNoGaps(m_noGapRegions, AT_, "forceNoGaps");
  MIntScratchSpace surpressingGaps(m_noGapRegions, AT_, "surpressingGaps");

  for(MInt region = 0; region < m_noGapRegions; region++) {
    forceNoGaps[region] = false;
    surpressingGaps[region] = 0;
    if(m_forceNoGaps > 0) {
      // specify time or crank-angle for which gap-Cells are surpressed!
      // -> ensures a clean gapOpneing!
      const MFloat elapsedTime = time();
      const MFloat cad = crankAngle(elapsedTime, 0);
      if(((m_gapSign[region] < F0 && (cad > m_gapAngleOpen[region] || cad < m_gapAngleClose[region]))
          || (m_gapSign[region] > F0 && cad > m_gapAngleOpen[region] && cad < m_gapAngleClose[region]))
         && m_forceNoGaps == 1) {
        forceNoGaps[region] = true;
      } else if(m_forceNoGaps == 2
                && ((cad > m_gapAngleOpen[region] && cad < m_gapAngleClose[region])
                    || (cad > m_gapAngleOpen[region + 1] && cad < m_gapAngleClose[region + 1]))) {
        forceNoGaps[region] = true;
      }
    }
  }

  // c) add all isInsideGap-Cells to the gapStack and set a_nearGapG
  for(MInt id = 0; id < a_noBandCells(0); id++) {
    MInt cellId = a_bandCellId(id, 0);
    if(isInsideGap(cellId)) {
      const MInt region = a_potentialGapCellClose(cellId) - 2;
      if(m_forceNoGaps > 0 && forceNoGaps[region]) {
        surpressingGaps[region]++;
        continue;
      }
      gapStack.push(cellId);
      isInsideGapTmp(cellId) = true;
      m_noGapCells++;
      gapCells.p[gapCellCounter++] = cellId;
      a_nearGapG(cellId) = true;
    }
  }

  if(m_forceNoGaps > 0) {
#if defined LS_DEBUG || !defined NDEBUG

    MPI_Allreduce(MPI_IN_PLACE, &surpressingGaps, m_noGapRegions, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "surpressingGaps");
    if(domainId() == 0) {
      for(MInt region = 0; region < m_noGapRegions; region++) {
        if(surpressingGaps[region] > 0)
          cerr << "Surpressing " << surpressingGaps[region] << " gap-Cells in region " << region << endl;
      }
    }
#endif
  }

  // d) add all neighbors of all gapStack-Cells within the gapWidth to the gapCell-List!
  while(cellsAdded) {
    cellsAdded = 0;

    if(gapCellCounter) {
      while(!gapStack.empty()) {
        MInt gCellId = gapStack.top();
        gapStack.pop();
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          MInt nghbrId = c_neighborId(gCellId, dir);
          if(nghbrId == -1) continue;
          if(a_isHalo(nghbrId)) continue;
          if(!a_potentialGapCell(nghbrId)) continue;

          if(m_gapInitMethod == 0) {
            if(abs(a_gapWidth(nghbrId)) <= eps2 && a_gapWidth(nghbrId) >= F0 && !isInsideGapTmp(nghbrId)) {
              gapCells.p[gapCellCounter++] = nghbrId;
              isInsideGapTmp(nghbrId) = true;
              m_noGapCells++;
              a_nearGapG(nghbrId) = true;
              gapStack.push(nghbrId);
              cellsAdded++;
            }

          } else {
            const MInt region = a_potentialGapCellClose(nghbrId) - 2;
            if(m_forceNoGaps > 0 && forceNoGaps[region]) continue;


            if(a_gapWidth(nghbrId) >= F0 && !isInsideGapTmp(nghbrId) && a_gapWidth(nghbrId) <= eps3) {
              // Timw: limit eps3
              // before: no upper-limit!
              gapCells.p[gapCellCounter++] = nghbrId;
              isInsideGapTmp(nghbrId) = true;
              m_noGapCells++;
              a_nearGapG(nghbrId) = true;
              gapStack.push(nghbrId);
              cellsAdded++;
            }
          }
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &cellsAdded, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "cellsAdded");

    // add Halo-Cells to the gap-Stack and List
    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noWindowCells(i); j++) {
        tmp_data(windowCellId(i, j)) = isInsideGapTmp(windowCellId(i, j));
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        tmp_data(windowId) = isInsideGapTmp(windowId);
      }
    }

    exchangeDataLS(&tmp_data(0), 1);

    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noHaloCells(i); j++) {
        const MInt haloCell = haloCellId(i, j);
        if(!isInsideGapTmp(haloCell)) {
          isInsideGapTmp(haloCell) = tmp_data(haloCell);
          if(isInsideGapTmp(haloCell)) {
            m_noGapCells++;
            gapCells.p[gapCellCounter++] = haloCell;
            a_nearGapG(haloCell) = true;
            gapStack.push(haloCell);
          }
        }
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        const MInt haloCell = grid().azimuthalHaloCell(i, j);
        if(!isInsideGapTmp(haloCell)) {
          isInsideGapTmp(haloCell) = tmp_data(haloCell);
          if(isInsideGapTmp(haloCell)) {
            m_noGapCells++;
            gapCells.p[gapCellCounter++] = haloCell;
            a_nearGapG(haloCell) = true;
            gapStack.push(haloCell);
          }
        }
      }
    }
  }


  if(m_gapInitMethod == 0) {
    MPI_Allreduce(MPI_IN_PLACE, &gapCellCounter, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "gapCellCounter");

    // f) add additional Cells to the gap-Cell-List
    if(gapCellCounter) {
      gapCellCounter = 0;
      MInt layerCount = 2;

      // first added to the local gapCells-counter
      for(MInt id = 0; id < a_noBandCells(0); id++) {
        MInt cellId = a_bandCellId(id, 0);
        if(a_nearGapG(cellId)) {
          gapCells.p[gapCellCounter++] = cellId;
        }
      }

      while(layerCount < maxNoLayers) {
        MInt layerCells = 0;
        for(MInt gc = 0; gc < gapCellCounter; gc++) {
          MInt gCellId = gapCells.p[gc];
          for(MInt dir = 0; dir < m_noDirs; dir++) {
            if(!a_hasNeighbor(gCellId, dir)) continue;
            MInt nghbrId = c_neighborId(gCellId, dir);
            if(a_isHalo(nghbrId)) continue;
            if(!a_nearGapG(nghbrId)) {
              lastLayer.p[layerCells++] = nghbrId;
              a_nearGapG(nghbrId) = layerCount;
              m_noGapCells++;
            }
          }
        }

        // exchange a_nearGapG-information
        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noWindowCells(i); j++) {
            tmp_data(windowCellId(i, j)) = a_nearGapG(windowCellId(i, j));
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
            MInt windowId = grid().azimuthalWindowCell(i, j);
            tmp_data(windowId) = a_nearGapG(windowId);
          }
        }

        exchangeDataLS(&tmp_data(0), 1);

        for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
          for(MInt j = 0; j < noHaloCells(i); j++) {
            const MInt haloCell = haloCellId(i, j);
            if(!a_nearGapG(haloCell)) {
              a_nearGapG(haloCell) = tmp_data(haloCell);
              if(a_nearGapG(haloCell)) {
                lastLayer.p[layerCells++] = haloCell;
                m_noGapCells++;
              }
            }
          }
        }
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
            const MInt haloCell = grid().azimuthalHaloCell(i, j);
            if(!a_nearGapG(haloCell)) {
              a_nearGapG(haloCell) = tmp_data(haloCell);
              if(a_nearGapG(haloCell)) {
                lastLayer.p[layerCells++] = haloCell;
                m_noGapCells++;
              }
            }
          }
        }

        for(MInt lc = 0; lc < layerCells; lc++) {
          gapCells.p[lc] = lastLayer.p[lc];
        }
        gapCellCounter = layerCells;
        layerCount++;
      }
    }
  }
}


//----------------------------------------------------------------------------


/**
 * \brief Level-set value of the G0-level-set function is decreased by slitly more than
 *        half of the minimum resolvable gap width (set by m_gapDeltaMin)
 *        For all G0-band Cells and additionally
 *        for all band-Cells of the first set, which are not already band-cells of the first set
 *        (But why?)
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetGapCorrect() {
  TRACE();

  //  MFloat eps2 =  c_cellLengthAtLevel(a_maxGCellLevel(0)) * (1.0 + sqrt(nDim)/2.0);
  //  MFloat eps2 =  c_cellLengthAtLevel(a_maxGCellLevel(0)) * (1.0 + sqrt(2.0)/2.0);

  const MFloat eps2 = (m_gapInitMethod < 2)
                          ? c_cellLengthAtLevel(a_maxGCellLevel(0)) * (5.0 / 4.0 * sqrt(nDim)) * m_gapDeltaMin * 0.51
                          : c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * m_gapDeltaMin;

  if(gapCellsExist()) {
    for(MInt id = 0; id < a_noBandCells(0); id++) {
      MInt cellId = a_bandCellId(id, 0);
      a_levelSetFunctionG(cellId, 0) -= eps2;
    }
    for(MInt id = 0; id < a_noBandCells(1); id++) {
      MInt cellId = a_bandCellId(id, 1);
      if(a_inBandG(cellId, 0)) continue;
      a_levelSetFunctionG(cellId, 0) -= eps2;
    }
  }
}


//----------------------------------------------------------------------------


/**
 * \brief Returns if there are any gap cells present (globally in all domains)
 *
 *
 * \author Claudia Guenther, Update: Tim Wegmann
 * \date 04/2011, 12/2018
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::gapCellsExist() {
  TRACE();

  MInt noGapCells = m_noGapCells;
  MPI_Allreduce(MPI_IN_PLACE, &noGapCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "noGapCells");
  return noGapCells;
}


//----------------------------------------------------------------------------


/**
 * \brief Returns if there are any gap cells present on this rank/domain
 *
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::localGapCellsExist() {
  TRACE();

  return m_noGapCells;
}

//----------------------------------------------------------------------------


/**
 * \brief Level-set value of the G0-level-set function is increased by slitly more than
 *        half of the minimum resolvable gap width (set by m_gapDeltaMin)
 *        For all G0-band Cells and additionally
 *        for all band-Cells of the first set, which are not already band-cells of the first set
 *        (But why?)
 *
 * \author Claudia Guenther
 * \date 04/2011
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::levelSetGapRecorrect() {
  TRACE();

  // MFloat eps2 =  c_cellLengthAtLevel(a_maxGCellLevel(0)) * (1.0 + sqrt(nDim)/2.0);
  // MFloat eps2 =  c_cellLengthAtLevel(a_maxGCellLevel(0)) * (1.0 + sqrt(2.0)/2.0);

  const MFloat eps2 = (m_gapInitMethod < 2)
                          ? c_cellLengthAtLevel(a_maxGCellLevel(0)) * (5.0 / 4.0 * sqrt(nDim)) * m_gapDeltaMin * 0.51
                          : c_cellLengthAtLevel(a_maxGCellLevel(0)) * sqrt(nDim) * m_gapDeltaMin;


  if(gapCellsExist()) {
    for(MInt id = 0; id < a_noBandCells(0); id++) {
      MInt cellId = a_bandCellId(id, 0);
      a_levelSetFunctionG(cellId, 0) += eps2;
    }

    for(MInt id = 0; id < a_noBandCells(1); id++) {
      MInt cellId = a_bandCellId(id, 1);
      if(a_inBandG(cellId, 0)) continue;
      a_levelSetFunctionG(cellId, 0) += eps2;
    }
  }
}


//---------------------------------------------------------------------------

/** \brief finalizes the initialization of the level-set method, sets m_computeSet triggers and identifies
 * potentialGapCells
 *
 * \author Claudia Guenther
 * \date 04/2012
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::finalizeLevelSetInitialization() {
  TRACE();

  // set compute set identification (previously, all set to 1 for initialization, now set to real values)
  // Timw: now moved to a later position in prepareLevelSet
  //      as this function is currently called in the initial-refinement-part!
  if(!m_levelSetMb) {
    for(MInt set = 0; set < m_noSets; set++)
      m_computeSet[set] = m_computeSet_tmp[set];
  }

  if(m_buildCollectedLevelSetFunction) {
    setUpPotentialGapCells();

    // during initialisation, startSet needs to be zero
    // only afterwards it is increased for buildCollectedLevelSet!
    // ASSERT(m_startSet == 0, " " );
    // m_startSet = 1;

  } else {
    if(m_maxNoSets > 1 || m_GFieldInitFromSTL) {
      MInt body = -1;
      for(MInt set = 0; set < m_noSets; set++) {
        if(m_noBodiesInSet[set] > 0) body = m_setToBodiesTable[set][0];
        for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
          a_bodyIdG(cellId, set) = body;
        }
      }
    } else {
      for(MInt set = 0; set < m_noSets; set++) {
        for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
          a_bodyIdG(cellId, set) = 0;
        }
      }
    }
  }
}


//---------------------------------------------------------------------------


/**
 * \fn void LsCartesianSolver::setUpPotentialGapCells()
 * \brief Set up cells, that may be tagged as gap cells during the solver run!
 * Initialises the arrays, according to the specifiec properties:
 * - a_potentialGapCell
 * - a_potentialGapCellClose
 *
 * \author Claudia Guenther, Update Tim Wegmann
 * \date   06/2012, 12/2018
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setUpPotentialGapCells() {
  TRACE();

  if(!m_closeGaps) return;

  // 1) Allocate potential gap cell arrays
  MBool& first = m_static_setUpPotentialGapCells_first;
  MFloat(&normal)[s_maxNoEmbeddedBodies * 3] = m_static_setUpPotentialGapCells_normal;
  MFloat(&center)[s_maxNoEmbeddedBodies * 3] = m_static_setUpPotentialGapCells_center;
  MFloat(&radius)[s_maxNoEmbeddedBodies] = m_static_setUpPotentialGapCells_radius;
  MFloat(&height)[s_maxNoEmbeddedBodies] = m_static_setUpPotentialGapCells_height;
  MFloat(&normalClose)[s_maxNoEmbeddedBodies * 3] = m_static_setUpPotentialGapCells_normalClose;
  MFloat(&centerClose)[s_maxNoEmbeddedBodies * 3] = m_static_setUpPotentialGapCells_centerClose;
  MFloat(&radiusClose)[s_maxNoEmbeddedBodies] = m_static_setUpPotentialGapCells_radiusClose;
  MFloat(&heightClose)[s_maxNoEmbeddedBodies] = m_static_setUpPotentialGapCells_heightClose;
  MInt(&bodyClose)[s_maxNoEmbeddedBodies] = m_static_setUpPotentialGapCells_bodyClose;
  MInt& noGapRegionsClose = m_static_setUpPotentialGapCells_noGapRegionsClose;

  MFloat rVec[3], hVec[3];
  if(first) {
    if(m_noBodyBndryCndIds > s_maxNoEmbeddedBodies) {
      mTerm(1, AT_, "Error: Too many embedded Bodies!");
    }

    // 2) Initialize potential gap cell arrays
    for(MInt body = 0; body < s_maxNoEmbeddedBodies; body++) {
      for(MInt dir = 0; dir < nDim; dir++) {
        normal[body * nDim + dir] = F0;
        center[body * nDim + dir] = F0;
        normalClose[body * nDim + dir] = F0;
        centerClose[body * nDim + dir] = F0;
      }
      radius[body] = 1.0;
      height[body] = 1.0;
      radiusClose[body] = 1.0;
      heightClose[body] = 1.0;
      normal[body * nDim + 0] = 1.0;
      normalClose[body * nDim + 0] = 1.0;
      bodyClose[body] = -1;
    }

    // 3) Read potential gap cell properties

    /*! \page propertiesLS
    \section gapRegionNormals
    Default: 0, 0, 0 \n
    Set the normal vectors of the gap region.\n
    Possible values are:
    <ul>
      <li>Coordinates [noGapRegions * nDim]</li>
    </ul>
    Keywords: <i>LEVEL-SET, GAP REGION</i>
    */
    for(MInt i = 0; i < m_noGapRegions; i++) {
      for(MInt j = 0; j < nDim; j++) {
        normal[i * nDim + j] = Context::getSolverProperty<MFloat>("gapRegionNormals", m_solverId, AT_,
                                                                  &normal[i * nDim + j], i * nDim + j);
      }
    }

    /*! \page propertiesLS
    \section gapRegionCenters
    Default: 0, 0, 0 \n
    Set the center of the gap region.\n
    Possible values are:
    <ul>
      <li>Coordinates [m_noGapRegions * nDim]</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, APE</i>
    */
    for(MInt i = 0; i < m_noGapRegions; i++) {
      for(MInt j = 0; j < nDim; j++) {
        center[i * nDim + j] = Context::getSolverProperty<MFloat>("gapRegionCenters", m_solverId, AT_,
                                                                  &center[i * nDim + j], i * nDim + j);
      }
    }
    /*! \page propertiesLS
     \section gapRegionRadii
     <code>MFloat LsCartesianSolver::setUpPotentialGapCells::radius</code>\n
     default:  1.0 \n \n
     Set the radius for each gap region, where any cell within the radius
     around the gap center will be treated as a potencial gap cell.  \n
     Keywords: <i>FINITE_VOLUME, GAP CELL</i>
   */
    for(MInt i = 0; i < m_noGapRegions; i++)
      radius[i] = Context::getSolverProperty<MFloat>("gapRegionRadii", m_solverId, AT_, &radius[i], i);


    for(MInt i = 0; i < m_noGapRegions; i++)
      height[i] = Context::getSolverProperty<MFloat>("gapRegionHeights", m_solverId, AT_, &height[i], i);


    /*! \page propertiesLS
    \section noGapRegionsClose
    <code>MInt LsCartesianSolver::setUpPotentialGapCells::noGapRegionsClose</code> \n
    Default: 1\n
    Sets the number of gap regions, which will be closed by the solver.\n
    Possible values are:
    <ul>
      <li>Any positive integer.</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, APE, GAP CELL</i>
    */
    noGapRegionsClose = 0;
    noGapRegionsClose = Context::getSolverProperty<MInt>("noGapRegionsClose", m_solverId, AT_, &noGapRegionsClose);
    ASSERT(noGapRegionsClose <= m_noEmbeddedBodies, "");

    for(MInt i = 0; i < noGapRegionsClose; i++)
      for(MInt j = 0; j < nDim; j++)
        normalClose[i * nDim + j] = Context::getSolverProperty<MFloat>("gapRegionNormalsClose", m_solverId, AT_,
                                                                       &normalClose[i * nDim + j], i * nDim + j);

    for(MInt i = 0; i < noGapRegionsClose; i++)
      for(MInt j = 0; j < nDim; j++)
        centerClose[i * nDim + j] = Context::getSolverProperty<MFloat>("gapRegionCentersClose", m_solverId, AT_,
                                                                       &centerClose[i * nDim + j], i * nDim + j);


    /*! \page propertiesLS
      \section gapRegionRadiiClose
      <code>MFloat LsCartesianSolver::setUpPotentialGapCells::radiusClose</code>\n
      default:  1.0 \n \n
      Set the radius for each gap closure region, where any cell within the radius
      around the gap closure center will be treated as a potencial gap closure cell.  \n
      Keywords: <i>FINITE_VOLUME, GAP CELL</i>
    */
    for(MInt i = 0; i < noGapRegionsClose; i++)
      radiusClose[i] = Context::getSolverProperty<MFloat>("gapRegionRadiiClose", m_solverId, AT_, &radiusClose[i], i);


    /*! \page propertiesLS
      \section gapRegionHeightsClose
      <code>MFloat LsCartesianSolver::setUpPotentialGapCells::gapRegionHeightsClose </code> \n
      Default: none \n
      Sets the heights of gap regions,  which will be closed by the solver.\n
      Possible values are:
      <ul>
      <li>floating point values </li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, APE, GAP CELL</i>
    */
    for(MInt i = 0; i < noGapRegionsClose; i++)
      heightClose[i] = Context::getSolverProperty<MFloat>("gapRegionHeightsClose", m_solverId, AT_, &heightClose[i], i);


    /*! \page propertiesLS
    \section gapRegionBodyClose
    <code>MInt LsCartesianSolver::setUpPotentialGapCells::gapRegionBodyClose </code> \n
    Default: none \n
    Sets the gap regions on the body which will be closed by the solver.\n
    Possible values are:
    <ul>
      <li>integers </li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, APE, GAP CELL</i>
    */
    for(MInt i = 0; i < m_noGapRegions; i++) {
      bodyClose[i] = Context::getSolverProperty<MInt>("gapRegionBodyClose", m_solverId, AT_, &bodyClose[i], i);
      ASSERT(bodyClose[i] + 1 <= m_noEmbeddedBodies, "");
    }

    first = false;
  }

  ASSERT(m_closeGaps, "");

  // 4) Loop over all cells and fill the arrays!

  for(MInt gCellId = 0; gCellId < a_noCells(); gCellId++) {
    // 4.1) Initilised with zero
    a_potentialGapCell(gCellId) = 0;
    a_potentialGapCellClose(gCellId) = 0;

    if(a_level(gCellId) == a_maxGCellLevel(0)) {
      for(MInt region = 0; region < m_noGapRegions; region++) {
        MFloat hCur = F0;
        MFloat rCur = F0;
        for(MInt dir = 0; dir < nDim; dir++) {
          rVec[dir] = c_coordinate(gCellId, dir) - center[region * nDim + dir];
          hCur += rVec[dir] * normal[region * nDim + dir];
        }
        for(MInt dir = 0; dir < nDim; dir++) {
          hVec[dir] = hCur * normal[region * nDim + dir];
          rVec[dir] -= hVec[dir];
          rCur += rVec[dir] * rVec[dir];
        }
        rCur = sqrt(rCur);
        if(rCur < radius[region] && abs(hCur) < height[region]) {
          // 4.2 ) If a cell is within in the GapRegion:
          // set a_potentialGapCell to regionId +1
          MInt gapRegionId = region + 1;
          a_potentialGapCell(gCellId) = gapRegionId;
        }
      }

      for(MInt region = 0; region < noGapRegionsClose; region++) {
        MFloat hCur = F0;
        MFloat rCur = F0;
        for(MInt dir = 0; dir < nDim; dir++) {
          rVec[dir] = c_coordinate(gCellId, dir) - centerClose[region * nDim + dir];
          hCur += rVec[dir] * normalClose[region * nDim + dir];
        }
        for(MInt dir = 0; dir < nDim; dir++) {
          hVec[dir] = hCur * normalClose[region * nDim + dir];
          rVec[dir] -= hVec[dir];
          rCur += rVec[dir] * rVec[dir];
        }
        rCur = sqrt(rCur);
        if(rCur < radiusClose[region] && abs(hCur) < heightClose[region]) {
          if(region >= m_noGapRegions) {
            a_potentialGapCellClose(gCellId) = m_G0regionId;
            a_potentialGapCell(gCellId) = 0;
          } else {
            MInt closeRegionId = region + 1;
            if(bodyClose[region] > -1) {
              closeRegionId = bodyClose[region] + 1;
            }
            a_potentialGapCellClose(gCellId) = closeRegionId;
            a_potentialGapCell(gCellId) = closeRegionId + m_noGapRegions;
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
//----------------------EXCHANGE-FUNCTIONS: ---------------------------

template <MInt nDim>
void LsCartesianSolver<nDim>::exchangeIntBuffers(MInt* sendBufferSize, MInt* receiveBufferSize, MInt tag,
                                                 MInt datasize) {
  TRACE();

  // debugging-version:
  // additional check of the buffer-sizes!
#ifdef LS_DEBUG
  std::ignore = datasize;

  MIntScratchSpace sendData(grid().noNeighborDomains(), AT_, "sendData");
  MIntScratchSpace receiveData(grid().noNeighborDomains(), AT_, "receiveData");
  MPI_Status status;

  // send the size of the data set
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    // cerr << "im sending to: " << grid().neighborDomain(i) << " i " << i << endl;
    sendData.p[i] = sendBufferSize[i];
    MPI_Issend(&sendData.p[i], 1, MPI_INT, grid().neighborDomain(i), tag, mpiComm(), &mpi_request[i], AT_,
               "sendData.p[i]");
  }

  // receive the size of the data set
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    // cerr << "im receiving from: " << grid().neighborDomain(i) << endl;
    MPI_Recv(&receiveData.p[i], (MInt)1, MPI_INT, grid().neighborDomain(i), tag, mpiComm(), &status, AT_,
             "receiveData.p[i]");
    receiveBufferSize[i] = receiveData.p[i];
  }

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Wait(&mpi_request[i], &status, AT_);
  }

  // send
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    if(sendBufferSize[i] == 0) {
      continue;
    }
    MPI_Issend(m_intSendBuffers[i], sendBufferSize[i], MPI_INT, grid().neighborDomain(i), tag, mpiComm(),
               &mpi_request[i], AT_, "m_intSendBuffers[i]");
  }
  // receive
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    if(receiveBufferSize[i] == 0) {
      continue;
    }
    MPI_Recv(m_intReceiveBuffers[i], receiveBufferSize[i], MPI_INT, grid().neighborDomain(i), tag, mpiComm(), &status,
             AT_, "m_intReceiveBuffers[i]");
  }
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Wait(&mpi_request[i], &status, AT_);
  }

#else

  std::ignore = sendBufferSize;
  std::ignore = receiveBufferSize;
  std::ignore = tag;

  // 1. send
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSizeW = noWindowCells(i) * datasize;
    if(bufSizeW == 0) continue;
    MPI_Issend(m_intSendBuffers[i], bufSizeW, MPI_INT, grid().neighborDomain(i), 0, mpiComm(), &mpi_request[i], AT_,
               "m_intSendBuffers[i]");
  }

  // 2. receive
  MPI_Status status;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSizeH = noHaloCells(i) * datasize;
    if(bufSizeH == 0) continue;
    MPI_Recv(m_intReceiveBuffers[i], bufSizeH, MPI_INT, grid().neighborDomain(i), 0, mpiComm(), &status, AT_,
             "m_intReceiveBuffers[i]");
  }

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSizeW = noWindowCells(i) * datasize;
    const MInt bufSizeH = noHaloCells(i) * datasize;

    if(bufSizeW == 0 || bufSizeH == 0) {
      continue;
    }

    MPI_Wait(&mpi_request[i], MPI_STATUS_IGNORE, AT_);
  }

  // exchangeBuffer(noNghbrDomains, nghbrDomains, noHaloCells, noWindowCells, comm, haloBuffer, &windowBuffer[0], noDat
  // );

#endif
}

template <MInt nDim>
template <typename T>
void LsCartesianSolver<nDim>::exchangeBuffersGlobal(T* sendBuffer, T* receiveBuffer, MInt* sendBufferSize,
                                                    MInt* receiveBufferSize, MInt* sndOffs, MInt* rcvOffs, MInt tag,
                                                    MInt offset) {
  TRACE();

  MInt ind;
  const MPI_Datatype DTYPE = maia::type_traits<T>::mpiType();

  ScratchSpace<MPI_Request> mpiSndReqGlob(grid().noDomains(), AT_, "mpi_requestGlobal");
  mpiSndReqGlob.fill(MPI_REQUEST_NULL);

#ifdef LS_DEBUG
  MPI_Status status;
  MIntScratchSpace sendData(grid().noDomains(), AT_, "sendData");
  MIntScratchSpace receiveData(grid().noDomains(), AT_, "receiveData");

  // send the size of the data set
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    sendData.p[i] = sendBufferSize[i];
    MPI_Issend(&sendData.p[i], 1, MPI_INT, i, tag, mpiComm(), &mpiSndReqGlob[i], AT_, "sendData.p[i]");
  }

  // receive the size of the data set
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    MPI_Recv(&receiveData.p[i], 1, MPI_INT, i, tag, mpiComm(), &status, AT_, "receiveData.p[i]");
    receiveBufferSize[i] = receiveData.p[i];
  }
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    MPI_Wait(&mpiSndReqGlob[i], &status, AT_);
  }

  // send
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    if(sendBufferSize[i] == 0) continue;
    ind = sndOffs[i] * offset;
    MPI_Issend(&sendBuffer[ind], sendBufferSize[i], DTYPE, i, (tag + 1), mpiComm(), &mpiSndReqGlob[i], AT_,
               "sendBuffer[ind]");
  }
  // receive
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    if(receiveBufferSize[i] == 0) continue;
    ind = rcvOffs[i] * offset;
    MPI_Recv(&receiveBuffer[ind], receiveBufferSize[i], DTYPE, i, (tag + 1), mpiComm(), &status, AT_,
             "receiveBuffer[ind]");
  }
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    MPI_Wait(&mpiSndReqGlob[i], &status, AT_);
  }
#else
  ScratchSpace<MPI_Request> mpiRcvReqGlob(grid().noDomains(), AT_, "mpi_requestGlobal");
  mpiRcvReqGlob.fill(MPI_REQUEST_NULL);

  for(MInt i = 0; i < grid().noDomains(); i++) {
    receiveBufferSize[i] = offset * (rcvOffs[i + 1] - rcvOffs[i]);
  }

  // send
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    if(sendBufferSize[i] == 0) continue;
    ind = sndOffs[i] * offset;

    // if ( offset*(sndOffs[i+1]-sndOffs[i]) != sendBufferSize[i] ) {
    /*if ( domainId() == 3 ) {
    cerr << "D:" << domainId() << " Snd:" << sendBufferSize[i] << " to " << i << " / " << offset << " " << sndOffs[i+1]
    << " " << sndOffs[i] << " " << m_globalSndOffsets[grid().noDomains()] << endl;
    }*/

    MPI_Isend(&sendBuffer[ind], sendBufferSize[i], DTYPE, i, (tag + 1), mpiComm(), &mpiSndReqGlob[i], AT_,
              "sendBuffer[ind]");
  }

  // receive
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    if(receiveBufferSize[i] == 0) continue;
    ind = rcvOffs[i] * offset;

    // if ( offset*(rcvOffs[i+1]-rcvOffs[i]) != receiveBufferSize[i] ) {
    /*if ( domainId() == 2 ) {
    cerr << "D:" << domainId() << " Rcv:" << receiveBufferSize[i] << " from " << i << " / " << offset << " " <<
    rcvOffs[i+1] << " " << rcvOffs[i] << " " << m_globalRcvOffsets[grid().noDomains()] <<endl;
    }*/
    MPI_Irecv(&receiveBuffer[ind], receiveBufferSize[i], DTYPE, i, (tag + 1), mpiComm(), &mpiRcvReqGlob[i], AT_,
              "receiveBuffer[ind]");
  }
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    MPI_Wait(&mpiSndReqGlob[i], MPI_STATUS_IGNORE, AT_);
    MPI_Wait(&mpiRcvReqGlob[i], MPI_STATUS_IGNORE, AT_);
  }
#endif
}


//----------------------------------------------------------------------------

/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::exchangeAllLevelSetData() {
  TRACE();

  if(grid().noNeighborDomains() == 0 && grid().noAzimuthalNeighborDomains() == 0) return;

  exchangeDataLS(&(a_levelSetFunctionG(0, 0)), m_maxNoSets);
  if(!m_combustion) {
    exchangeDataLS(&(a_bodyIdG(0, 0)), m_maxNoSets);
  }

  if(m_semiLagrange) {
    exchangeDataLS(&(a_oldLevelSetFunctionG(0, 0)), m_maxNoSets);
  }
}

/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::exchangeLevelSet() {
  TRACE();

  if(grid().noNeighborDomains() == 0 && grid().noAzimuthalNeighborDomains() == 0) return;

  exchangeDataLS(&(a_levelSetFunctionG(0, 0)), m_maxNoSets);
}

/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::exchangeLs(MFloat* dataField, MInt firstset, MInt dataBlockSize) {
  TRACE();

  if(grid().noNeighborDomains() == 0 && grid().noAzimuthalNeighborDomains() == 0) return;

  MInt sendBufferCounter = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferCounter = 0;
    for(MInt j = 0; j < noWindowCells(i); j++) {
      memcpy((void*)&m_gSendBuffers[i][sendBufferCounter],
             (void*)&dataField[IDX_LSSET(windowCellId(i, j), firstset)],
             dataBlockSize * sizeof(MFloat));
      sendBufferCounter += dataBlockSize;
    }
  }

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSize = noWindowCells(i) * dataBlockSize;
    if(bufSize == 0) continue;
    MPI_Issend(m_gSendBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 0, mpiComm(), &mpi_request[i], AT_,
               "m_gSendBuffers[i]");
  }

  MPI_Status status;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSize = noHaloCells(i) * dataBlockSize;
    if(bufSize == 0) continue;
    MPI_Recv(m_gReceiveBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 0, mpiComm(), &status, AT_,
             "m_gReceiveBuffers[i]");
  }

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    const MInt bufSizeW = noWindowCells(i) * dataBlockSize;
    const MInt bufSizeH = noHaloCells(i) * dataBlockSize;
    if(bufSizeH == 0 || bufSizeW == 0) {
      continue;
    }
    MPI_Wait(&mpi_request[i], &status, AT_);
  }

  MInt receiveBufferCounter = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    receiveBufferCounter = 0;
    for(MInt j = 0; j < noHaloCells(i); j++) {
      memcpy((void*)&dataField[IDX_LSSET(haloCellId(i, j), firstset)],
             (void*)&m_gReceiveBuffers[i][receiveBufferCounter],
             dataBlockSize * sizeof(MFloat));
      receiveBufferCounter += dataBlockSize;
    }
  }

  if(grid().azimuthalPeriodicity()) {
    this->exchangeAzimuthalPer(&dataField[0], dataBlockSize, firstset);
  }
}
//----------------------------------------------------------------------------
/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::exchangeGapInfo() {
  TRACE();

  if(grid().noNeighborDomains() == 0 && grid().noAzimuthalNeighborDomains() == 0) return;

  if(!m_closeGaps) return;

  MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noWindowCells(i); j++) {
      tmp_data[windowCellId(i, j)] = a_nearGapG(windowCellId(i, j));
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
      MInt windowId = grid().azimuthalWindowCell(i, j);
      tmp_data(windowId) = a_nearGapG(windowId);
    }
  }

  exchangeDataLS(&tmp_data[0], 1);

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      a_nearGapG(haloCellId(i, j)) = tmp_data[haloCellId(i, j)];
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      const MInt haloCell = grid().azimuthalHaloCell(i, j);
      a_nearGapG(haloCell) = tmp_data[haloCell];
    }
  }
}

//----------------------------------------------------------------------------


/** \fn void LsCartesianSolver<nDim>::readLevelSetProperties()
 * \brief reads in the level set properties
 *
 * \author unknown, Stephan Schlimpert
 * \date unknown, September 2011, January 2013
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::readLevelSetProperties() {
  TRACE();

  MBool tmpFalse = false;
  MBool tmpTrue = true;

  m_referenceLength = F1;
  /*! \page propertiesLS
    \section referenceLength
    <code>MFloat LsCartesianSolver::m_referenceLength </code>\n
    default = <code>1.0</code>\n \n
    WARNING: Do NOT use any value different than 1.0 - The correct implementation of this is not checked, so it probably
    will not do what you think it does/should do. Don't use it unless you REALLY know what you are doing. Reference
    Length L - The length = 1.0 of the grid is scaled with L. Possible values are: <ul> <li>1.0 +- eps</li>
    </ul>
    Keywords: <i>LEVELSET, VARIABLES</i>
  */
  m_referenceLength = Context::getSolverProperty<MFloat>("referenceLength", m_solverId, AT_, &m_referenceLength);

  m_referenceVelocity = F1;
  /*! \page propertiesLS
    \section referenceLength
    <code>MFloat FvCartesianSolver::m_Pr </code>\n
    default = <code>sqrt(1/(1+0.5*(1.4-1)/2*Ma^2</code>\n \n
    Reference velocity used for non-dimensionalisation of some of the movement functions!
    NOTE: this is not consistently used, some versions use the velocity, other the speed of sound
          (both in the infinity state) as reference Velocity!
    Keywords: <i>LEVELSET, VARIABLES</i>
  */
  if(Context::propertyExists("referenceVelocity", m_solverId)) {
    m_referenceVelocity =
        Context::getSolverProperty<MFloat>("referenceVelocity", m_solverId, AT_, &m_referenceVelocity);
  } else {
    if(m_levelSetMb || m_levelSetFv) {
      // update referenceVelocity based on Ma
      // labels:LS,DOC NOTE: @Thomas this hack changes when a different initialCondition is used!
      //      best to specify a referenceVelocity instead!
      const MFloat Ma = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_);
      const MFloat gammaMinusOne = 1.4 - 1.0;
      const MFloat TInfinity = 1.0 / (1.0 + 0.5 * gammaMinusOne * POW2(Ma));
      m_referenceVelocity = sqrt(TInfinity);
    }
  }

  /*! \page propertiesLS
    \section virtualSurgery
    <code>MBool LsCartesianSolver::m_virtualSurgery </code>\n
    default = <code>false</code>\n \n
    Activates the virtual surgery environment of the level-set solver.
    When activated, the solver modifies the 0th level-set field by interpolating
    between the 1st and 2nd level-set field. The 1st and 2nd level-set field are static.
    The 0th level-set field must not be a collected level-set.
    Keywords: <i>LEVELSET, VARIABLES</i>
  */
  m_virtualSurgery = false;
  m_virtualSurgery = Context::getSolverProperty<MBool>("virtualSurgery", m_solverId, AT_, &m_virtualSurgery);

  /*! \page propertiesLS
    \section sphereRadiusLimit
    <code>MFloat LsCartesianSolver::m_sphereRadiusLimit </code>\n
    default = <code>5.0</code>\n \n
    The sphereRadiusLimit defines the maximum area in which the level-set is defined.
    This factor can be set, if the band should be small to save cells, but the level-set
    value outside of the band should still be calculated. Is only used in mode 4 of
    constructGFieldFromSTL.
    Keywords: <i>LEVELSET, VARIABLES</i>
  */
  m_sphereRadiusLimit = 5.0;
  m_sphereRadiusLimit = Context::getSolverProperty<MFloat>("sphereRadiusLimit", m_solverId, AT_, &m_sphereRadiusLimit);

  if(m_virtualSurgery) {
    /*! \page propertiesLS
      \section approxNoInterpRegions
      <code>MInt LsCartesianSolver::m_approxNoInterpRegions </code>\n
      default = <code></code>\n \n
      Number of Interpolations regions in virtual surgery environment. The number can only
      be approximated a priori.
      Keywords: <i>LEVELSET, VARIABLES</i>
    */
    m_approxNoInterpReg = Context::getSolverProperty<MInt>("approxNoInterpRegions", m_solverId, AT_);

    if(m_approxNoInterpReg <= 0) {
      stringstream errorMessage;
      errorMessage << "ERROR: Invalid number of interpolation regions. Set approxNoInterpRegions > 0!!!";
      mTerm(1, AT_, errorMessage.str());
    }
    /*! \page propertiesLS
      \section interpStartTimes
      <code>MInt LsCartesianSolver::m_interpStartTime </code>\n
      default = <code></code>\n \n
      This property defines the start time step of a virtual surgery for each interpolation
      region.
      Keywords: <i>LEVELSET, VARIABLES</i>
    */
    if(Context::propertyLength("interpStartTimes", m_solverId) == m_approxNoInterpReg) {
      mAlloc(m_interpStartTime, m_approxNoInterpReg, "m_interpStartTime", -1, AT_);
      for(MInt i = 0; i < m_approxNoInterpReg; i++) {
        m_interpStartTime[i] = Context::getSolverProperty<MInt>("interpStartTimes", m_solverId, AT_, i);
      }
    } else {
      stringstream errorMessage;
      errorMessage << "ERROR: Array length of interpStartTimes and approxNoInterpRegions do not match!!";
      mTerm(1, AT_, errorMessage.str());
    }

    /*! \page propertiesLS
      \section noInterpTimeSteps
      <code>MInt LsCartesianSolver::m_noInterpTimeSteps </code>\n
      default = <code></code>\n \n
      This property defines the number of time steps  used for the interpolation in a virtual
      surgery for each interpolation region.
      Keywords: <i>LEVELSET, VARIABLES</i>
    */
    if(Context::propertyLength("noInterpTimeSteps", m_solverId) == m_approxNoInterpReg) {
      mAlloc(m_noInterpTimeSteps, m_approxNoInterpReg, "m_noInterpTimeSteps", -1, AT_);
      for(MInt i = 0; i < m_approxNoInterpReg; i++) {
        m_noInterpTimeSteps[i] = Context::getSolverProperty<MInt>("noInterpTimeSteps", m_solverId, AT_, i);
      }
    } else {
      stringstream errorMessage;
      errorMessage << "ERROR: Array length of noInterpTimeSteps and approxNoInterpRegions do not match!!";
      mTerm(1, AT_, errorMessage.str());
    }
  }

  // property which determines whether the bodyId and old-LevelsetFuntion is relevant
  m_semiLagrange = (string2enum(solverMethod()) == MAIA_SEMI_LAGRANGE_LEVELSET
                    || string2enum(solverMethod()) == MAIA_RUNGE_KUTTA_MB_SEMI_LAGRANGE_LEVELSET
                    || string2enum(solverMethod()) == MAIA_SEMI_LAGRANGE_LEVELSET_LB);

  m_flameSpeed = 0.0;
  m_flameSpeed = Context::getSolverProperty<MFloat>("flameSpeed", m_solverId, AT_, &m_flameSpeed);

  m_jetHalfWidth = 0.5;
  m_radiusFlameTube = 0.5;
  m_yOffsetFlameTube = 0.04;
  m_xOffsetFlameTube = 0.0;
  m_marksteinLength = 0.0;
  m_marksteinLength = Context::getSolverProperty<MFloat>("marksteinLength", m_solverId, AT_, &m_marksteinLength);
  m_marksteinLengthPercentage = F1;
  m_marksteinLengthPercentage =
      Context::getSolverProperty<MFloat>("marksteinLengthPercentage", m_solverId, AT_, &m_marksteinLengthPercentage);
  m_marksteinLength *= m_marksteinLengthPercentage;

  m_xOffsetFlameTube = Context::getSolverProperty<MFloat>("xOffsetFlameTube", m_solverId, AT_, &m_xOffsetFlameTube);
  m_xOffsetFlameTube2 = Context::getSolverProperty<MFloat>("xOffsetFlameTube2", m_solverId, AT_, &m_xOffsetFlameTube);
  m_yOffsetFlameTube = Context::getSolverProperty<MFloat>("yOffsetFlameTube", m_solverId, AT_, &m_yOffsetFlameTube);
  m_yOffsetFlameTube2 = Context::getSolverProperty<MFloat>("yOffsetFlameTube2", m_solverId, AT_, &m_yOffsetFlameTube);
  m_radiusFlameTube = Context::getSolverProperty<MFloat>("radiusFlameTube", m_solverId, AT_, &m_radiusFlameTube);
  m_radiusFlameTube2 = Context::getSolverProperty<MFloat>("radiusFlameTube2", m_solverId, AT_, &m_radiusFlameTube);


  m_trackMovingBndry = true;
  m_trackMbStart = -1;
  m_trackMbEnd = numeric_limits<MInt>::max();

  m_trackMovingBndry = Context::getSolverProperty<MBool>("trackMovingBndry", m_solverId, AT_, &m_trackMovingBndry);
  m_trackMbStart = Context::getSolverProperty<MInt>("trackMbStart", m_solverId, AT_, &m_trackMbStart);
  m_trackMbEnd = Context::getSolverProperty<MInt>("trackMbEnd", m_solverId, AT_, &m_trackMbEnd);

  if(m_levelSetMb) {
    m_constructGField = true;
    m_constructGField = Context::getSolverProperty<MBool>("constructGField", m_solverId, AT_, &m_constructGField);
  }

  m_realRadiusFlameTube = 0.5;
  m_realRadiusFlameTube =
      Context::getSolverProperty<MFloat>("realRadiusFlameTube", m_solverId, AT_, &m_realRadiusFlameTube);
  /*! \page propertiesLS
      \section forcing
      <code>MInt LsCartesianSolver::m_forcing </code>\n
      default = <code>no default</code>\n \n
      Activates the excitation / forcing via sponge in 2D. \n \n
      Possible values are:
      <ul>
      <li> 0: inactive </li>
      <li> 1: active </li>
      </ul>
Keywords: <i> LEVELSET, FLAME, SPONGE</i>
  */
  m_forcing = false;
  m_forcing = Context::getSolverProperty<MBool>("forcing", m_solverId, AT_, &m_forcing);

  m_flameRadiusOffset = F0;
  m_flameRadiusOffset = Context::getSolverProperty<MFloat>("flameRadiusOffset", m_solverId, AT_, &m_flameRadiusOffset);

  m_twoFlames = false;
  m_twoFlames = Context::getSolverProperty<MBool>("twoFlames", m_solverId, AT_, &tmpFalse);

  m_initialFlameHeight = F1;
  m_initialFlameHeight =
      Context::getSolverProperty<MFloat>("initialFlameHeight", m_solverId, AT_, &m_initialFlameHeight);

  m_initialCondition = Context::getSolverProperty<MInt>("initialCondition", m_solverId, AT_);

  /*! \page propertiesLS
    \section jetHalfLength
    <code>MFloat FvCartesianSolver::m_jetHalfLength </code>\n
    default = <code>4.165</code>\n \n
    parameter used in generation of various jet inflows
    possible values are:
    <ul>
    <li> any positive float </li>
    </ul>
    Keywords: <i> jet  </i>
  */
  m_jetHalfLength = 4.165;
  m_jetHalfLength = Context::getSolverProperty<MFloat>("jetHalfLength", m_solverId, AT_, &m_jetHalfLength);

  m_filterFlameTubeEdgesDistance = -9999.9;
  m_filterFlameTubeEdgesDistance = Context::getSolverProperty<MFloat>("filterFlameTubeEdgesDistance", m_solverId, AT_,
                                                                      &m_filterFlameTubeEdgesDistance);

  m_filterFlameTubeEdges = Context::getSolverProperty<MBool>("filterFlameTubeEdges", m_solverId, AT_, &tmpFalse);

  m_noReactionCells = 0.026367201;
  m_noReactionCells = Context::getSolverProperty<MFloat>("noReactionCells", m_solverId, AT_, &m_noReactionCells);


  m_dampingDistanceFlameBase = 0.259;
  m_dampingDistanceFlameBase =
      Context::getSolverProperty<MFloat>("dampingDistanceFlameBase", m_solverId, AT_, &m_dampingDistanceFlameBase);

  m_dampingDistanceFlameBaseExtVel = 0.05;
  m_dampingDistanceFlameBaseExtVel = Context::getSolverProperty<MFloat>("dampingDistanceFlameBaseExtVel", m_solverId,
                                                                        AT_, &m_dampingDistanceFlameBaseExtVel);

  m_useCorrectedBurningVelocity = 0;
  m_useCorrectedBurningVelocity =
      Context::getSolverProperty<MBool>("useCorrectedBurningVelocity", m_solverId, AT_, &m_useCorrectedBurningVelocity);

  m_plenum = false;
  m_plenum = Context::getSolverProperty<MBool>("plenum", m_solverId, AT_, &tmpFalse);

  m_timeStepMethod = Context::getSolverProperty<MInt>("timeStepMethod", m_solverId, AT_);
  m_cfl = Context::getSolverProperty<MFloat>("cfl", m_solverId, AT_);


  m_steadyFlameLength = -F1;
  /*! \page propertiesLS
    \section steadyFlameLength
    <code>MFloat LsCartesianSolver::m_steadyFlameLength </code>\n
    no default \n \n
    This variable defines the steady flame front amplitude \n
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>COMBUSTION, STEADY, FLAME, FRONT, AMPLITUDE</i>
  */
  m_steadyFlameLength = Context::getSolverProperty<MFloat>("steadyFlameLength", m_solverId, AT_, &m_steadyFlameLength);

  m_maxNoCells = maxNoGridCells();

  m_maxNoSets = 1;
  /*! \page propertiesLS
    \section maxNoLevelSets
    <code>MInt LsCartesianSolver::m_maxNoSets </code>\n
    default = <code>1</code>\n \n
    This property fixes the maximum number of separate level set functions to be computed cells \n
    possible values are:
    <ul>
    <li>Positive integer values, should be O(1) -> more than 10 separate sets will become very slow and should be
    unnencessary</li>
    </ul>
    Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
  */
  m_maxNoSets = Context::getSolverProperty<MInt>("maxNoLevelSets", m_solverId, AT_, &m_maxNoSets);


  m_STLReinitMode = 2;
  /*! \page propertiesLS
     \section GFieldFromSTLReInitMode
     default = 2 \n \n
     Choose method how to Reinit Levelset from STL.\n
     In the Reinitialization part:
     </ul>
     <li>0: low order reinitialization </li>
     <li>1: low order reinitialization on all levels is included,
            to avoid initialization errors!  </li>
     <li>2: no reinitialisation at all (recommended version!) </li>
     <li>3: uses high order reinitialization </li>
     </ul>
     Keywords: <i>LEVELSET</i>
    */

  m_STLReinitMode = Context::getSolverProperty<MInt>("GFieldFromSTLReinitMode", m_solverId, AT_, &m_STLReinitMode);

  m_noSets = m_maxNoSets;
  m_startSet = 0;
  m_buildCollectedLevelSetFunction = false;
  m_determineG0CellsMode = 0;
  m_noBodyBndryCndIds = 1;
  m_closeGaps = false;

  if(m_levelSetMb) {
    m_noEmbeddedBodies = 1;
  }

  mAlloc(m_levelSetSign, m_maxNoSets, "m_levelSetSign", 1, AT_);
  mAlloc(m_computeSet, m_maxNoSets, "m_computeSet", true, AT_);
  mAlloc(m_computeSet_tmp, m_maxNoSets, "m_computeSet_tmp", true, AT_);
  mAlloc(m_computeSet_backup, m_maxNoSets, "m_computeSet_backup", true, AT_);
  mAlloc(m_changedSet, m_maxNoSets, "m_changedSet", true, AT_);

  if(this->m_adaptation) {
    if(!m_restart) {
      m_initialRefinement = true;
    }
  }

  /*! \page propertiesLS
  \section levelSetSign
  <code>MInt* LsCartesianSolver::m_levelSetSign</code>\n
  default = <code>1</code>\n \n
  This property specifies if the sign of the level-set function during initialization or reconstruction should be
  assumed to be as computed e.g. from the stl geometry (1) or inverted (-1).\n Possible values are: <ul> <li> 1 </li>
  <li> -1 </li>
  </ul>
  Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS, MOVING_BOUNDARY</i>
  */
  for(MInt i = 0; i < m_noSets; i++) {
    m_levelSetSign[i] = Context::getSolverProperty<MInt>("levelSetSign", m_solverId, AT_, &m_levelSetSign[i], i);
  }

  /*! \page propertiesLS
  \section computeSet
  <code>MBool* LsCartesianSolver::m_computeSet</code>\n
  default = <code>1</code>\n \n
  This property specifies if the level set function should be transported in time or if initialization at startup is
  sufficient.\n Possible values are: <ul> <li> 1 </li> <li> 0 </li>
  </ul>
  Keywords: <i>LEVELSET</i>
  */
  for(MInt i = 0; i < m_noSets; i++) {
    m_computeSet_tmp[i] = Context::getSolverProperty<MBool>("computeSet", m_solverId, AT_, &m_computeSet_tmp[i], i);
    m_computeSet_backup[i] = Context::getSolverProperty<MBool>("computeSet", m_solverId, AT_, &m_computeSet_tmp[i], i);
  }

  /*! \page propertiesLS
  \section GFieldInitFromSTL
  <code>MInt LsCartesianSolver::m_GFieldInitFromSTL </code>\n
  default = <code>0</code>\n \n
  This property triggers if the level-set function(s) are initialized from an .stl-geometry.
  <ul>
  <li>off - level-set functions are not initialized from .stl-geometries </li>
  <li>on - level-set functions are initialized from .stl-geometries </li>
  </ul>
  Keywords: <i>LEVELSET, LEVEL_SET_FROM_STL </i>
  */
  m_GFieldInitFromSTL = tmpFalse;
  m_GFieldInitFromSTL = Context::getSolverProperty<MBool>("GFieldInitFromSTL", m_solverId, AT_, &m_GFieldInitFromSTL);

  /*! \page propertiesLS
  \section GFieldFromSTLInitCheck
  <code>MBool LsCartesianSolver::m_GFieldFromSTLInitCheck </code>\n
  default = <code>0</code>\n \n
  This property triggers if an output is generated after the initialization. Careful: produces valgrind errors as only
  the level-set function is initialized, the other variables which are written out are not. <ul> <li>off - no check
  output is generated</li> <li>on - check output is generated after the initialization </li>
  </ul>
  Keywords: <i>LEVELSET, LEVEL_SET_FROM_STL </i>
  */
  m_GFieldFromSTLInitCheck = Context::getSolverProperty<MBool>("GFieldFromSTLInitCheck", m_solverId, AT_, &tmpFalse);

  if(m_GFieldInitFromSTL) {
    /*! \page propertiesLS
    \section levelSetWithSTLCorrection
    <code>MBool LsCartesianSolver::m_GWithReConstruction </code>\n
    default = <code>0</code>\n \n
    This property triggers if the level-set function(s) are corrected using the .stl-geometry
    <ul>
    <li>off - level-set stl correction is off </li>
    <li>on - level-set stl correction is on </li>
    </ul>
    Keywords: <i>LEVELSET, LEVEL_SET_FROM_STL </i>
    */
    m_GWithReConstruction = Context::getSolverProperty<MBool>("levelSetWithSTLCorrection", m_solverId, AT_, &tmpFalse);


    m_reconstructBand = Context::getSolverProperty<MInt>("reconstructBand", m_solverId, AT_, &m_reconstructBand);
  }

  if(m_maxNoSets > 1) {
    /*! \page propertiesLS
    \section buildCollectedLevelSetFunction
    <code>MInt LsCartesianSolver::m_buildCollectedLevelSetFunction </code>\n
    default = <code>0</code>\n \n
    This property triggers if the level-set functions from m_startSet to m_noSets are collected in set 0 \n
    possible values are:
    <ul>
    <li>off - 0th level-set function is computed normally and not collected </li>
    <li>on - 0th level-set function is collected based on the other level-set functions </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET </i>
    */
    m_buildCollectedLevelSetFunction = Context::getSolverProperty<MBool>("buildCollectedLevelSetFunction", m_solverId,
                                                                         AT_, &m_buildCollectedLevelSetFunction);
  }

  MInt movingBndryCndId = 3006;
  /*! \page propertiesLS
    \section movingBndryCndId
    <code>MInt* LsCartesianSolver::movingBndryCndId </code>\n
    default = <code>3006</code>\n \n
    This property specifies the standard moving boundary condition Id -> required for the identification of moving
    objects. If no other bndryCndIds are specified in \ref bodyBndryCndIds, this Id is assumed to be the Id of all
    moving bodies (in this case, just one body).\n Possible values are: <ul> <li>Valid boundary condition ids (see
    boundary condition ids in the code) -> not necessarily a moving bndry cnd Id, this is only required if you work with
    the MB part of the code! </li>
    </ul>
    Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
    */
  movingBndryCndId = Context::getSolverProperty<MInt>("movingBndryCndId", m_solverId, AT_, &movingBndryCndId);

  m_noEmbeddedBodies = 1;
  if(!Context::propertyExists("bodyBndryCndIds", m_solverId) && m_levelSetMb) {
    /*! \page propertiesLS
      \section noEmbeddedBodies
      <code>MInt* LsCartesianSolver::noEmbeddedBodies </code>\n
      default = <code>1</code>\n \n
      This property specifies the number of bodies released into the flow field.\n
      Possible values are:
      <ul>
      <li> 1 </li>
      </ul>
      Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
    */
    m_noEmbeddedBodies = Context::getSolverProperty<MInt>("noEmbeddedBodies", m_solverId, AT_, &m_noEmbeddedBodies);
    m_noBodyBndryCndIds = m_noEmbeddedBodies;

    if(m_noBodyBndryCndIds > 0) {
      mAlloc(m_bodyBndryCndIds, m_noBodyBndryCndIds, "m_bodyBndryCndIds", movingBndryCndId, AT_);
    }
  } else {
    /*! \page propertiesLS
      \section bodyBndryCndIds
      <code>MInt* LsCartesianSolver::m_bodyBndryCndIds</code>\n
      default = <code>movingBndryCndId</code>\n \n
      This property specifies the moving boundary condition Ids -> required for the identification of moving objects. If
      no bndryCndIds are specified here, the Id specified in \ref movingBndryCndId is assumed to be the Id of all moving
      bodies (in this case, just one body).\n Possible values are: <ul> <li> all valid Boundary Condition Ids  -> not
      necessarily a moving bndry cnd Id, this is only required if you work with the MB part of the code! \n \n
      </li>
      </ul>
      Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS, MOVING_BOUNDARY</i>
    */
    m_noBodyBndryCndIds = Context::propertyLength("bodyBndryCndIds", m_solverId);
    if(m_noBodyBndryCndIds > 0) {
      mAlloc(m_bodyBndryCndIds, m_noBodyBndryCndIds, "m_bodyBndryCndIds", 0, AT_);
    }
    for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
      m_bodyBndryCndIds[i] = Context::getSolverProperty<MInt>("bodyBndryCndIds", m_solverId, AT_, i);
    }
    if(m_noBodyBndryCndIds > 1) {
      m_noEmbeddedBodies = m_noBodyBndryCndIds;
    }
  }

  MBool useBodyToSetTable = Context::propertyExists("bodyBndryCndIds", m_solverId);
  if((m_GFieldInitFromSTL && useBodyToSetTable) || m_maxNoSets > 1) {
    setUpBodyToSetTable();
  }

  /*! \page propertiesLS
  \section highOrderDeltaFunction
  <code>MBool LsCartesianSolver::m_highOrderDeltaFunction </code>\n
  default = <code>0</code>\n \n
  test trigger for high order arclength calculcation which equals the flame surface area \n
  possible values are:
  <ul>
  <li>0 - off</li>
  <li>1 - on</li>
  </ul>
  Keywords: <i>LEVELSET, COMBUSTION, DELTA, FUNCTION, FLAME, SURFACE, AREA</i>
  */
  m_highOrderDeltaFunction = Context::getSolverProperty<MBool>("highOrderDeltaFunction", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
    \section fourthOrderNormalCurvatureComputation
    <code>MBool LsCartesianSolver::m_fourthOrderNormalCurvatureComputation </code>\n
    default = <code>0</code>\n \n
    test trigger for fourth order normal and curvature computation, should not be used for combustion \n
    possible values are:
    <ul>
    <li>0 - off</li>
    <li>1 - on</li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, FOURTH, ORDER, NORMAL, CURVATURE</i>
  */
  m_fourthOrderNormalCurvatureComputation =
      Context::getSolverProperty<MBool>("fourthOrderNormalCurvatureComputation", m_solverId, AT_, &tmpFalse);

  if(m_combustion && m_fourthOrderNormalCurvatureComputation) {
    stringstream errorMessage;
    errorMessage << "ERROR: fourthOrderNormalCurvatureCompuatation should not be used, pockets are not correctly "
                    "generated ... exiting";
    mTerm(1, AT_, errorMessage.str());
  }

  /*! \page propertiesLS
    \section curvatureDamp
    <code>MBool LsCartesianSolver::m_curvatureDamp </code>\n
    default = <code>0</code>\n \n
    temporarly trigger for damping the curvature and heat release\n
    possible values are:
    <ul>
    <li>0 - off</li>
    <li>1 - on</li>
    </ul>
    Keywords: <i>COMBUSTION, LEVELSET, DAMPING, CURVATURE</i>
  */
  m_curvatureDamp = Context::getSolverProperty<MBool>("curvatureDamp", m_solverId, AT_, &tmpFalse);

  m_curvatureDampFactor = F3;
  /*! \page propertiesLS
    \section curvatureDampFactor
    <code>MFloat LsCartesianSolver::curvatureDampFactor </code>\n
    default = <code>0</code>\n \n
    temporarly trigger for damping the curvature by a factor \n
    possible values are:
    <ul>
    <li> Non-negative floating point values of the order of 3.0 </li>
    </ul>
    Keywords: <i>COMBUSTION, LEVELSET, DAMPING, CURVATURE</i>
  */
  m_curvatureDampFactor =
      Context::getSolverProperty<MFloat>("curvatureDampFactor", m_solverId, AT_, &m_curvatureDampFactor);

  /*! \page propertiesLS
    \section sharpDamp
    <code>MBool LsCartesianSolver::m_sharpDamp </code>\n
    default = <code>0</code>\n \n
    temporarly trigger for sharp damping the curvature and heat release \n
    possible values are:
    <ul>
    <li>0 - off</li>
    <li>1 - on</li>
    </ul>
    Keywords: <i>COMBUSTION, LEVELSET, SHARP, DAMPING, CURVATURE</i>
  */
  m_sharpDamp = Context::getSolverProperty<MBool>("sharpDamp", m_solverId, AT_, &tmpFalse);


  /*! \page propertiesLS
    \section useLocalMarksteinLength
    <code>MBool LsCartesianSolver::m_useLocalMarksteinLength </code>\n
    default = <code>0</code>\n \n
    temporarly trigger for using the local markstein length when computing the neutral markstein length, (see testcases
    1990, 19901) \n possible values are: <ul> <li>0 - off</li> <li>1 - on</li>
    </ul>
    Keywords: <i>COMBUSTION, LEVELSET, LOCAL, MARKSTEIN, LENGTH</i>
  */
  m_useLocalMarksteinLength = Context::getSolverProperty<MBool>("useLocalMarksteinLength", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
    \section hyperbolicCurvature
    <code>MBool LsCartesianSolver::m_hyperbolicCurvature </code>\n
    default = <code>0</code>\n \n
    test trigger for the use of the hyperbolic extension of the curvature itselfs. \n
    possible values are:
    <ul>
    <li>0 - off</li>
    <li>1 - on</li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, CURVATURE, HYPERBOLIC, EXTENSION </i>
  */
  m_hyperbolicCurvature = Context::getSolverProperty<MBool>("hyperbolicCurvature", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
    \section gRKMethod
    <code>MInt LsCartesianSolver::m_gRKMethod </code>\n
    default = <code>no default value</code>\n \n
    This property sets the level set runge kutta method\n
    possible values are:
    <ul>
    <li>0 - </li>
    <li>1 - three step third order TVD Runge-kutta scheme</li>
    <li>2 - </li>
    <li>3 - </li>
    <li>4 - </li>
    <li>5 - no runge kutta is solved (fixes the level set function) </li>
    </ul>
    Keywords: <i>LEVELSET, RUNGE, KUTTA, SCHEME </i>
  */

  m_gRKMethod = Context::getSolverProperty<MInt>("gRKMethod", m_solverId, AT_);

  /*! \page propertiesLS
    \section levelSetDiscretizationScheme
    <code>MString LsCartesianSolver::m_levelSetDiscretizationScheme </code>\n
    default = <code>no default value</code>\n \n
    This property sets the level set discretization scheme of the level set equation.  \n
    possible values are:
    <ul>
    <li>US1 - first order upwind scheme</li>
    <li>UC3 - third order upwind scheme</li>
    <li>UC3_SB - third order upwind scheme</li>
    <li>UC5 - fifth order upwind scheme</li>
    <li>UC5_SB - fifth order upwind scheme</li>
    <li>WENO5 - fifth order upwind scheme</li>
    <li>WENO5_SB - fifth order upwind scheme</li>
    <li>UC11 - higher order upwind scheme</li>
    </ul>
    Keywords: <i>LEVELSET, DISCRETIZATION, METHOD, EQUATION, RHS </i>
  */

  m_levelSetDiscretizationScheme = Context::getSolverProperty<MString>("levelSetDiscretizationScheme", m_solverId, AT_);


  /*! \page propertiesLS
    \section LsRotate
    <code>MString LsCartesianSolver::m_LsRotate </code>\n
    default = <code> flase </code>\n \n
    This property needs to be true for rotating levelset
  */
  m_LsRotate = false;
  m_LsRotate = Context::getSolverProperty<MBool>("LsRotate", m_solverId, AT_, &m_LsRotate);

  if(m_LsRotate) {
    m_reconstructOldG = Context::getSolverProperty<MBool>("reconstructOldG", m_solverId, AT_, &m_reconstructOldG);
  } else {
    m_reconstructOldG = false;
  }

  /*! \page propertiesLS
    \section gBandWidth
    <code>MInt LsCartesianSolver::m_gBandWidth </code>\n
    default = <code>no default value</code>\n \n
    This property sets the level set band width in which the level set equation is solved. The outer band cells are
    reset to a constant level set function value G and are not taken into account for the level set computation. E.g.
    gBandWidth = 15 means 15 cells in each direction of the G=0 cells are put into the level set band  \n possible
    values are: <ul> <li> Any non-negative floating point smaller than the number of time steps \ref timeSteps </li>
    </ul>
    Keywords: <i>LEVELSET, BAND, WIDTH </i>
  */
  m_gBandWidth = Context::getSolverProperty<MInt>("gBandWidth", m_solverId, AT_);

  /*! \page propertiesLS
    \section gShadowWidth
    <code>MInt LsCartesianSolver::m_gShadowWidth </code>\n
    default = <code>no default value is set</code>\n \n
    This property sets the level set shadow width in which the level set cells are refined to the maximum G-Cell
    refinement level \ref maxGCellLevel.  \n possible values are: <ul> <li> Any non-negative floating point number which
    greater than the level set band width \ref m_gBandWith</li>
    </ul>
    Keywords: <i>LEVELSET, SHADOW, WIDTH, MAXIMUM, REFINEMENT, LEVEL </i>
  */

  m_gShadowWidth = Context::getSolverProperty<MInt>("gShadowWidth", m_solverId, AT_);

  /*! \page propertiesLS
    \section gShadowWidthRans
    <code>MInt LsCartesianSolver::m_gShadowWidthRans </code>\n
    default = <code>no default value is set</code>\n \n
    This property sets the level set shadow width for levelSetRans in which the level set cells are refined to the
    maximum G-Cell refinement level \ref maxGCellLevel.  \n possible values are: <ul> <li> Any non-negative floating
    point number which greater than the level set band width \ref m_gBandWith</li>
    </ul>
    Keywords: <i>LEVELSET, SHADOW, WIDTH, MAXIMUM, REFINEMENT, LEVEL </i>
  */
  if(m_levelSetRans) {
    m_gShadowWidthRans = Context::getSolverProperty<MInt>("gShadowWidthRans", m_solverId, AT_);
  }

  /*! \page propertiesLS
    \section gInnerBound
    <code>MInt LsCartesianSolver::m_gInnerBound </code>\n
    default = <code>no default value</code>\n \n
    .  \n
    possible values are:
    <ul>

    </ul>
    Keywords: <i>LEVELSET, METHOD </i>
  */
  m_gInnerBound = Context::getSolverProperty<MInt>("gInnerBound", m_solverId, AT_);

  ASSERT(m_gInnerBound >= 2, "Should be at least 2!");


  /*! \page propertiesLS
    \section maxGCellLevel
    <code>MInt LsCartesianSolver::m_maxGCellLevel </code>\n
    default = <code>no default value</code>\n \n
    This property sets the maximum level set refinement level of the level set grid. The grid width is calculated via
    \f$ m_gCellDistance = c_cellLengthAtLevel(m_maxGCellLevel) \f$. \n The maximum level set refinement lebel is limited
    by the maximum refinement level \ref maxRfnmtLevel  \n possible values are: <ul> <li> Any non-negative integer value
    smaller than the maximum refinement level </li>
    </ul>
    Keywords: <i>LEVELSET, MAXIMUM, GCELL, REFINEMENT, LEVEL </i>
  */
  if(Context::propertyLength("maxGCellLevel", m_solverId)) {
    if(Context::propertyLength("maxGCellLevel", m_solverId) == 1) {
      TERMM(1, "Warning: maxGCellLevel has only one set! Use maxRfnmntLvl instead!");
    }

    m_gCellLevelJump = true;

    ASSERT(m_semiLagrange && !m_LsRotate, "");
    // otherwise untested

    mAlloc(m_maxGCellLevel, m_maxNoSets, "m_maxGCellLevel", 0, AT_);

    ASSERT(Context::propertyLength("maxGCellLevel", m_solverId) == m_maxNoSets, "");

    for(MInt set = 0; set < m_maxNoSets; set++) {
      m_maxGCellLevel[set] = Context::getSolverProperty<MInt>("maxGCellLevel", m_solverId, AT_, set);
      ASSERT(m_maxGCellLevel[set] <= maxRefinementLevel() && m_maxGCellLevel[set] >= minLevel(), "");
    }

  } else {
    mAlloc(m_maxGCellLevel, 1, "m_maxGCellLevel", 0, AT_);
    m_maxGCellLevel[0] = Context::getSolverProperty<MInt>("maxRfnmntLvl", m_solverId, AT_);

    ASSERT(m_maxGCellLevel[0] == maxRefinementLevel(),
           "maxGCellLevel " + to_string(m_maxGCellLevel[0]) + " maxLevel " + to_string(maxRefinementLevel()));
  }

  m_maxLevelChange = false;
  if(globalTimeStep > 0 && m_restartFile && isActive() && maxLevel() != maxRefinementLevel()) {
    m_maxLevelChange = true;
    if(domainId() == 0) {
      cerr << "MaxLevel is " << maxLevel() << " and maxRefinementLevel is " << maxRefinementLevel()
           << " => temporarily reduction of m_maxGCellLevel to MaxLevel!" << endl;
    }
    m_maxGCellLevel[0] = maxLevel();
  }

  if(this->m_noSensors > 0) {
    if(globalTimeStep > 0 && m_restartFile && isActive() && this->m_adaptation
       && maxLevel() > this->m_maxSensorRefinementLevel[0]) {
      m_maxLevelChange = true;
      if(domainId() == 0) {
        cerr << "MaxLevel is " << maxLevel() << " and maxSensorLevel is " << this->m_maxSensorRefinementLevel[0]
             << " => temporarily increase of m_maxGCellLevel to MaxLevel!" << endl;
      }
      ASSERT(m_maxGCellLevel[0] == maxLevel(), "");
    }
  }

  m_gCellDistance = c_cellLengthAtLevel(m_maxGCellLevel[0]);
  m_FgCellDistance = F1 / m_gCellDistance;
  m_log << "smallest gCell length " << m_gCellDistance << endl;

  m_outsideGValue = (MFloat)(2 * m_gBandWidth * m_gCellDistance);
  m_log << "Outside G value: " << m_outsideGValue << endl;

  m_computeExtVel = 1;
  /*! \page propertiesLS
    \section computeExtVel
    <code>MInt LsCartesianSolver::m_computeExtVel </code>\n
    default = <code>1</code>\n \n
    This property sets the level set method to compute the corrected burning velocity and the extension velocity of the
    level set cells. \n possible values are: <ul> <li>1  - computes extension velocity (curvature effect controlled by
    the Markstein length), \ref marksteinLength</li> <li>10 - computes extension velocity (curvature effect controlled
    by the Markstein length), \ref marksteinLength</li> <li>2  - computes extension velocity for the combustion solver
    FVCombstnSolver (curvature effect controlled by the Markstein length), \ref marksteinLength</li> <li>5  -
    computes extension velocity for the G-Equation Progress Variable Model FVGequPvSolver (cuvrature effect
    controlled by the molecular diffusivity) </li> <li>6  - computes extension velocity for the G-Equation Progress
    Variable Model FVGequPvSolver (no curvature effect)</li> <li>60 - computes extension velocity for the G-Equation
    Progress Variable Model FVGequPvSolver (no curvature effect)</li> <li>7  - computes extension velocity for the
    G-Equation Progress Variable Model FVGequPvSolver (curvature effect controlled by the Markstein length), \ref
    marksteinLength</li> <li>70 - computes extension velocity for the G-Equation Progress Variable Model
    FVGequPvSolver (curvature effect controlled by the Markstein length + stretch effects), \ref
    marksteinLength</li> <li>71 - computes extension velocity for the G-Equation Progress Variable Model
    FVGequPvSolver (curvature effect controlled by the Markstein length for symmetric flames), \ref marksteinLength
    </li> <li>8  - computes extension velocity for the G-Equation Progress Variable Model FVGequPvSolver
    (computation of neutral Markstein length and growth rate of perturbed flame surfaces), \ref marksteinLength</li>
    <li>81 - computes extension velocity for the G-Equation Progress Variable Model FVGequPvSolver (stability test
    with use of neutral Markstein length (enforced version of case 8), neutral Markstein length is computed for every
    time step), \ref marksteinLength</li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, EXTENSION, VELOCITY, NEUTRAL, MARKSTEIN, LENGTH </i>
  */
  m_computeExtVel = Context::getSolverProperty<MInt>("computeExtVel", m_solverId, AT_, &m_computeExtVel);

  m_smoothExtVel = true;
  /*! \page propertiesLS
    \section smoothExtVel
    <code>MInt LsCartesianSolver::m_smoothExtVel </code>\n
    default = <code>1</code>\n \n
    This property triggers the smoothed computation of the corrected burning velocity in the band layer which is done by
    different computations, see \ref computeExtVel \n possible values are: <ul> <li>0 - computation of extension
    velocity only at the G zero front cells </li> <li>1 - computation of smoothed extension velocity for the band layer
    cells, see \ref computeExtVel </li>
    </ul>
    Keywords: <i>LEVELSET, SMOOTHED, EXTENSION, VELOCITY, BAND, LAYER, G0, ZERO, FRONT, CELLS </i>
  */
  m_smoothExtVel = Context::getSolverProperty<MBool>("smoothExtVel", m_solverId, AT_, &tmpTrue);

  m_extVelIterations = 10000;
  /*! \page propertiesLS
    \section extVelIterations
    <code>MInt LsCartesianSolver::m_extVelIterations </code>\n
    default = <code>10000</code>\n \n
    controls the maximum iterations to reach the convergence criteria of the hyperbolic extension, \ref extVelCFL, \ref
    extVelConvergence \n possible values are: <ul> <li>Non-negative integer values of the order of 1000 </li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, CURVATURE, HYPERBOLIC, EXTENSION </i>
  */
  m_extVelIterations = Context::getSolverProperty<MInt>("extVelIterations", m_solverId, AT_, &m_extVelIterations);

  /*! \page propertiesLS
    \section extVelConvergence
    <code>MInt LsCartesianSolver::m_extVelConvergence </code>\n
    no default valaue used\n \n
    controls the convergence of the hyperbolic extension, \ref extVelCFL, \ref extVelConvergence \n
    possible values are:
    <ul>
    <li>Non-negative floating point values of the order of 0.1</li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, COMBUSTION, HYPERBOLIC, EXTENSION </i>
  */
  m_extVelConvergence = Context::getSolverProperty<MFloat>("extVelConvergence", m_solverId, AT_);

  m_extVelCFL = 0.5;
  /*! \page propertiesLS
    \section extVelCFL
    <code>MFloat LsCartesianSolver::m_extVelCFL </code>\n
    default = <code>0.5</code>\n \n
    controls the CFL number of the pseudo time discretization of the hyperbolic extension, \ref extVelConvergence, \ref
    extVelIterations \n possible values are: <ul> <li>Non-negative floating point values of the order of 0.1</li>
    </ul>
    Keywords: <i>LEVELSET, COMBUSTION, COMBUSTION, HYPERBOLIC, EXTENSION </i>
  */
  m_extVelCFL = Context::getSolverProperty<MFloat>("extVelCFL", m_solverId, AT_, &m_extVelCFL);

  /*! \page propertiesLS
    \section reinitMethod
    <code>MString LsCartesianSolver::m_reinitMethod </code>\n
    default = <code>no default value</code>\n \n
    This property sets the reinitialization method, for more information see Dissertation D. Hartmann \n
    possible values are:
    <ul>
    <li>CR1 </li>
    <li>CR2 </li>
    <li>HCR1 </li>
    <li>HCR2 </li>
    <li>HCR2_LIMITED </li>
    <li>HCR2_FULLREINIT </li>
    <li>SUS5CR </li>
    <li>CR2PLUS </li>
    <li>RSU </li>
    <li>DL1 </li>
    <li>DL2 </li>
    <li>SUS_1 </li>
    <li>SUS_1PLUS </li>
    <li>SUS_2 </li>
    <li>SUS_WENO5 </li>
    <li>SUS_WENO5PLUS </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, METHOD </i>
  */
  m_reinitMethod = Context::getSolverProperty<MString>("reinitMethod", m_solverId, AT_);

  m_noGapRegions = 0;
  m_gapInitMethod = 2;
  m_forceNoGaps = 0;

  if(m_buildCollectedLevelSetFunction) {
    /*! \page propertiesLS
    \section gapReinitMethod
    <code>MString LsCartesianSolver::m_gapReinitMethod </code>\n
    default = <code>m_reinitMethod</code>\n \n
    This property sets the reinitialization method for the gap reinitialization for multiple level-set functions, for
    more information see Dissertation D. Hartmann \n possible values are: <ul> <li>CR1 </li> <li>CR2 </li> <li>HCR1
    </li> <li>HCR2 </li> <li>HCR2_LIMITED </li> <li>HCR2_FULLREINIT </li> <li>SUS5CR </li> <li>CR2PLUS </li> <li>RSU
    </li> <li>DL1 </li> <li>DL2 </li> <li>SUS_1 </li> <li>SUS_1PLUS </li> <li>SUS_2 </li> <li>SUS_WENO5 </li>
    <li>SUS_WENO5PLUS </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, METHOD </i>
  */
    m_gapReinitMethod = Context::getSolverProperty<MString>("gapReinitMethod", m_solverId, AT_, &m_reinitMethod);


    /*! \page propertiesLS
    \section gapReinitMethod
    <code>MBool LsCartesianSolver::m_closeGaps </code>\n
    default = <code>false</code>\n \n
    This property triggers if the reinitialization based closure procedure for narrow gaps between different bodies is
    called for the collected level-set function <ul> <li>0 </li> <li>1 </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET </i>
    */
    m_closeGaps = Context::getSolverProperty<MBool>("closeGaps", m_solverId, AT_, &m_closeGaps);

    if(!m_closeGaps) {
      /*! \page propertiesLS
        \section determineG0Mode
        <code>MBool LsCartesianSolver::m_determineG0CellsMode </code>\n
        default = <code>false</code>\n \n
        Triggers if G0 cells of combined level-set function are determined from phi^0 (mode=0) or if they are inherited
        from the individual level-set functions (mode=1) <ul> <li>0 </li> <li>1 </li>
        </ul>
        Keywords: <i>LEVELSET, MULTILEVELSET </i>
      */
      m_determineG0CellsMode =
          Context::getSolverProperty<MBool>("determineG0Mode", m_solverId, AT_, &m_determineG0CellsMode);
    } else { // m_closeGaps

      /*! \page propertiesLS
      \section noGapRegions
      <code>MFloat noGapRegions</code>\n
      Default: -1 \n
      Set how many regions can be gap regions.\n
      Possible values are:
      <ul>
        <li>Any positive integer number.</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, APE</i>
      */

      m_noGapRegions = Context::getSolverProperty<MInt>("noGapRegions", m_solverId, AT_, &m_noGapRegions);

      ASSERT(m_noGapRegions > 0, "");

      /*! \page propertiesLS
        \section G0regionId
        <code>MInt LsCartesianSolver::m_G0regionId </code>\n
        default = <code>-1</code>\n
        defines a region for which the combined levelSet is used to smooth the overlapping of
        different bodies. However no gapClosing where the combined levelSet is additionally modified
        is applied!
        NOTE: possibly a sponge needs to a applied in this region if any of the bodies is moving!
        <ul>
        <li> -1: not used </li> (1)
        <li> positiv initger: regionId for which this is used!</li>
        </ul>
        Keywords: <i> MOVING BOUNDARY, GAP, FINITE_VOLUME, LEVELSET</i>
      */
      if(Context::propertyExists("G0regionId", m_solverId)) {
        m_G0regionId = Context::getSolverProperty<MInt>("G0regionId", m_solverId, AT_, &m_G0regionId);
      }


      /*! \page propertiesLS
        \section gapInitMethod
        <code>MInt LsCartesianSolver::m_gapInitMethod </code>\n
        default = <code>2</code>\n
        Possible values are:
        <ul>
        <li> 0: all gap Cells are reseted and initialized (older version from Claudia) </li> (1)
        <li> 1: only new arising gap Cells are reseted! </li> (2)
        <li> 2: only new arising gap Cells are reseted and mew epsilons are used</li>
        </ul>
        Keywords: <i> MOVING BOUNDARY, GAP, FINITE_VOLUME</i>
      */

      m_gapInitMethod = Context::getSolverProperty<MInt>("gapInitMethod", 0, AT_, &m_gapInitMethod);

      /*! \page propertiesLS
        \section forceNoGaps
        <code>MInt LsCartesianSolver::m_forceNoGaps </code>\n
        default = <code>0</code>\n
        triggers the surpression of gap cells
        <ul>
        <li> 0: do not surpress gap cells </li> (1)
        <li> 1: surpress gap cells for multi-valve engine </li> (2)
        <li> 2: surpress gap cells for single-valve engine</li>
        </ul>
        Keywords: <i> MOVING BOUNDARY, GAP, LS</i>
      */
      m_forceNoGaps = Context::getSolverProperty<MInt>("forceNoGaps", m_solverId, AT_, &m_forceNoGaps);

      if(m_forceNoGaps == 1) { // surpress gap cells for multi-valve engine

        mAlloc(m_gapAngleClose, m_noGapRegions, "m_gapAngleClose", F0, AT_);
        mAlloc(m_gapAngleOpen, m_noGapRegions, "m_gapAngleOpen", F0, AT_);
        mAlloc(m_gapSign, m_noGapRegions, "m_gapSign", F1, AT_);

        for(MInt i = 0; i < m_noGapRegions; i++) {
          m_gapAngleClose[i] =
              Context::getSolverProperty<MFloat>("gapAngleClose", m_solverId, AT_, &m_gapAngleClose[i], i);
          m_gapAngleOpen[i] =
              Context::getSolverProperty<MFloat>("gapAngleOpen", m_solverId, AT_, &m_gapAngleClose[i], i);
          if(m_gapAngleClose[i] < 0.0) {
            m_gapAngleClose[i] = 720 + m_gapAngleClose[i];
          }
          if(m_gapAngleOpen[i] < 0.0) {
            m_gapAngleOpen[i] = 720 + m_gapAngleOpen[i];
          }
          if(m_gapAngleOpen[i] > m_gapAngleClose[i]) {
            m_gapSign[i] = -1;
          }
        }
      } else if(m_forceNoGaps == 2) { // surpress gap cells for single-valve engine

        ASSERT(m_noGapRegions == 1, "Mode only intended for 1 gap-region!");

        mAlloc(m_gapAngleClose, m_noGapRegions + 1, "m_gapAngleClose", F0, AT_);
        mAlloc(m_gapAngleOpen, m_noGapRegions + 1, "m_gapAngleOpen", F0, AT_);
        mAlloc(m_gapSign, m_noGapRegions + 1, "m_gapSign", F1, AT_);

        ASSERT(m_noGapRegions + 1 == Context::propertyLength("gapAngleOpen", m_solverId), "");

        for(MInt i = 0; i < m_noGapRegions + 1; i++) {
          m_gapAngleClose[i] =
              Context::getSolverProperty<MFloat>("gapAngleClose", m_solverId, AT_, &m_gapAngleClose[i], i);
          m_gapAngleOpen[i] =
              Context::getSolverProperty<MFloat>("gapAngleOpen", m_solverId, AT_, &m_gapAngleClose[i], i);
          if(m_gapAngleClose[i] < 0.0) {
            m_gapAngleClose[i] = 720 + m_gapAngleClose[i];
          }
          if(m_gapAngleOpen[i] < 0.0) {
            m_gapAngleOpen[i] = 720 + m_gapAngleOpen[i];
          }
          if(m_gapAngleOpen[i] > m_gapAngleClose[i]) {
            m_gapSign[i] = -1;
          }
        }
      } else {
        mAlloc(m_gapAngleClose, m_noGapRegions, "m_gapAngleClose", F0, AT_);
        mAlloc(m_gapAngleOpen, m_noGapRegions, "m_gapAngleOpen", F0, AT_);
        mAlloc(m_gapSign, m_noGapRegions, "m_gapSign", F1, AT_);
      }
    }

    m_gapDeltaMin = F1;
    m_gapDeltaMinOrig = F1;
    /*! \page propertiesLS
      \section deltaMin
      <code>MBool LsCartesianSolver::m_gapDeltaMin </code>\n
      default = <code>1</code>\n \n
      <ul>
      <li>0 </li>
      <li>1 </li>
      </ul>
      Keywords: <i>LEVELSET, MULTILEVELSET </i>
    */
    if(Context::propertyExists("deltaMin", m_solverId)) {
      m_gapDeltaMin = Context::getSolverProperty<MFloat>("deltaMin", m_solverId, AT_, &m_gapDeltaMin);
    }
    m_gapDeltaMinOrig = m_gapDeltaMin;
  }

  /*! \page propertiesLS
    \section gReinitIterations
    <code>MInt LsCartesianSolver::m_gReinitIterations </code>\n
    default = <code>no default value</code>\n \n
    This property triggers the maximum number of reinitialization steps to reach the convergence criterion \ref
    reinitConvergence \n possible values are: <ul> <li> Any non-negative integer value</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, ITERATION, STEPS </i>
  */
  m_gReinitIterations = Context::getSolverProperty<MInt>("gReinitIterations", m_solverId, AT_);

  m_minReinitializationSteps = m_gReinitIterations;
  /*! \page propertiesLS
    \section minReinitializationSteps
    <code>MInt LsCartesianSolver::m_minReinitializationSteps </code>\n
    default = <code> m_gReinitIterations </code>\n \n
    This property triggers the minimum number of reinitialization steps to reach the convergence criterion \ref
    reinitConvergence \n possible values are: <ul> <li> Any non-negative integer value</i>
    </ul>
    Keywords: <i>LEVELSET, MINIMUM REINITIALIZATION, ITERATION, STEPS </i>
  */
  m_minReinitializationSteps =
      Context::getSolverProperty<MInt>("minReinitializationSteps", m_solverId, AT_, &m_minReinitializationSteps);

  m_maintenanceIterations = m_gReinitIterations;
  /*! \page propertiesLS
    \section maintenanceIterations
    <code>MInt LsCartesianSolver::m_maintenanceIterations </code>\n
    default = <code>gReinitIterations</code>\n \n
    This property triggers the maximum number of iterations for maintaining the outer band layers \n
    possible values are:
    <ul>
    <li> Any non-negative integer value</li>
    </ul>
    Keywords: <i>LEVELSET, MAINTENANCE, ITERATION, STEPS </i>
  */
  m_maintenanceIterations =
      Context::getSolverProperty<MInt>("maintenanceIterations", m_solverId, AT_, &m_maintenanceIterations);

  /*! \page propertiesLS
    \section guaranteeReinit
    <code>MBool LsCartesianSolver::m_guaranteeReinit </code>\n
    default = <code>0</code>\n \n
    this property guarantees the reinitialization at every time step \n
    possible values are:
    <ul>
    <li>0 - off</li>
    <li>1 - on</li>
    </ul>
    Keywords: <i>LEVELSET, GUARANTEE, REINITIALIZATION </i>
  */
  m_guaranteeReinit = Context::getSolverProperty<MBool>("guaranteeReinit", m_solverId, AT_, &tmpFalse);

  if(!m_guaranteeReinit && m_closeGaps) {
    ASSERT(m_gapDeltaMin < 1.00000001, "Reinitialisation is required if deltaMin > 1 is used!");
  }

  /*! \page propertiesLS
    \section reinitCFL
    <code>MInt LsCartesianSolver::m_reinitCFL </code>\n
    default = <code>1</code>\n \n
    This property sets the reinitialization cfl number which is used to calculate the pseudo time step for the
    reinitialization  \f$ \Delta t_{pseudo} = CFL_{reinit} * \Delta x \f$, with  \f$ \Delta x \f$ being the smallest G
    cell width, see \ref maxGCellLevel  \n possible values are: <ul> <li> Any non-negative floating point number smaller
    than 1</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, CFL, PSEUDO, TIME, STEP </i>
  */
  m_reinitCFL = Context::getSolverProperty<MFloat>("reinitCFL", m_solverId, AT_);

  m_intermediateReinitIterations = 1;
  /*! \page propertiesLS
      \section intermediateReinitIterations
      <code>MInt LsCartesianSolver::m_intermediateReinitIterations </code>\n
      default = <code>1</code>\n \n
      This property triggers the maximum number of iterations for reinitializing the level set function at every runge
     kutta step. Only used for the runge kutta method 3, \ref gRKMethod. \n possible values are: <ul> <li> Any
     non-negative integer value</li>
      </ul>
      Keywords: <i>LEVELSET, REINITIALIZATION, ITERATION, STEPS, RUNGE, KUTTA </i>
    */
  m_intermediateReinitIterations = Context::getSolverProperty<MInt>("intermediateReinitIterations", m_solverId, AT_,
                                                                    &m_intermediateReinitIterations);

  /*! \page propertiesLS
    \section reinitConvergence
    <code>MFloat LsCartesianSolver::m_reinitConvergence </code>\n
    default = <code>no default value</code>\n \n
    This property sets the convergence criterion of the reinitialization. The criterion should be of the order of
    10^{-3} - 10^{-10}, see also \ref gReinitIterations \n possible values are: <ul> <li> Any non-negative floating
    point number smaller than 1</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, CONVERGENCE, CRITERION, ITERATION, STEPS </i>
  */
  m_reinitConvergence = Context::getSolverProperty<MFloat>("reinitConvergence", m_solverId, AT_);
  m_reinitConvergenceReset = m_reinitConvergence;

  m_reinitThreshold = F0;
  /*! \page propertiesLS
    \section reinitThreshold
    <code>MFloat LsCartesianSolver::m_reinitThreshold </code>\n
    default = <code>0</code>\n \n
    This property sets the reinitialization threshold. The important parameter for the reinitialization is the
    maintenance of the property \f$ \left| \nabla G \right| = 1\f$. \n The deviation to this property is the so called
    reinitialization threshold calculated for all level set front cells (G0 cells) via  \n \f$ \Delta_{thresh} = \left(
    \left| \nabla G \right| - 1 \right)^{2} \f$. If the maximum threshold and the averaged threshold, see \ref
    reinitThresholdAvg, \n is reached the level set will be reinitialized otherwise maintained. \n possible values are:
    <ul>
    <li> Any non-negative floating point number </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, THRESHOLD, G0, CELLS </i>
  */
  m_reinitThreshold = Context::getSolverProperty<MFloat>("reinitThreshold", m_solverId, AT_, &m_reinitThreshold);

  m_reinitThresholdAvg = F0;
  /*! \page propertiesLS
    \section reinitThresholdAvg
    <code>MFloat LsCartesianSolver::m_reinitThresholdAvg </code>\n
    default = <code>0</code>\n \n
    This property sets the averaged reinitialization threshold. The important parameter for the reinitialization is the
    maintenance of the property \f$ \left| \nabla G \right| = 1\f$. \n The deviation to this property is the so called
    reinitialization threshold calculated for all level set front cells (G0 cells), see \ref reinitThreshold, and the
    averaged threshold  \n \f$ \overline{\Delta_{threshAvg}} =  \overline{\left( \left| \nabla G \right| - 1
    \right)^{2}} \f$ of the G0 cells. \n If both criterions are reached the level set will be reinitialized otherwise
    maintained. \n possible values are: <ul> <li> Any non-negative floating point number </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, THRESHOLD, G0, CELLS </i>
  */
  m_reinitThresholdAvg =
      Context::getSolverProperty<MFloat>("reinitThresholdAvg", m_solverId, AT_, &m_reinitThresholdAvg);

  m_omegaReinit = F1;
  /*! \page propertiesLS
    \section omegaReinit
    <code>MInt LsCartesianSolver::m_omegaReinit </code>\n
    default = <code>1</code>\n \n
    This property is used for the constrained reinitialization methods Cr1, CR2, CR2PLUS and RSU, see \ref reinitMethod.
    \n possible values are: <ul> <li> Any non-negative floating point number </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, OMEGA </i>
  */
  m_omegaReinit = Context::getSolverProperty<MFloat>("omegaReinit", m_solverId, AT_, &m_omegaReinit);

  m_relaxationFactor = 0.5;
  /*! \page propertiesLS
    \section relaxationFactor
    <code>MInt LsCartesianSolver::m_relaxationFactor </code>\n
    default = <code>0.5</code>\n \n
    This property forces the reinitialization.\n
    possible values are:
    <ul>
    <li> Any non-negative floating point number smaller than one</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, RELAXATION, FACTOR, CELLS, TO, CORRECT </i>
  */
  m_relaxationFactor = Context::getSolverProperty<MFloat>("relaxationFactor", m_solverId, AT_, &m_relaxationFactor);

  m_levelSetTestCase = 0;
  /*! \page propertiesLS
    \section levelSetTestCase
    <code>MInt LsCartesianSolver::m_levelSetTestCase </code>\n
    default = <code>0</code>\n \n
    This property triggers the computation of the propagation speed for a few test cases and triggers the calculation of
    a few level set ouput informations, see file levelSetData2 in LsCartesianSolver::writeLevelSet(). \n possible
    values are: <ul> <li>0 - no test case chosen</li> <li>1 - ?</li> <li>200 - diagonal convection</li> <li>201 - merge
    problem</li> <li>202 - oscillating circle problem, sine wave around a cirlce</li> <li>203 - Zalesak's problem </li>
    <li>204 - Enright's test</li>
    <li>205 - Harmonic forcing problem (Preetham07)</li>
    <li>206 - simple expansion/compression?</li>
    <li>207 - periodic array of 5 vortices in a 3x1 domain</li>
    <li>208 - rotation</li>
    <li>1001 - ?</li>
    <li>1002 - ?</li>
    </ul>
    3D test cases:
    <ul>
    <li> 300 - convection </li>
    <li> 301,30100 - merge/expand </li>
    <li> 302 - merge/shrink </li>
    <li> 303 - Oscillating circle problem </li>
    <li> 305 - Enright test </li>
    <li> 306 - Cold vortex-flame interaction </li>
    <li> 366 - Oscillating sphere </li>
    <li> 302 - Pancake sphere </li>
    <li> 3005 - Enright's test by STL </li>
    <li> 3006 - Zalesak Disk by STL, refered to Dupont's and Back compensation </li>
    <li> 3009 - General rigid body motion template </li>
    </ul>
    Keywords: <i>LEVELSET, CASES, PROPAGATION, SPEED, OUTPUT </i>
  */
  m_levelSetTestCase = Context::getSolverProperty<MInt>("levelSetTestCase", m_solverId, AT_, &m_levelSetTestCase);

  m_levelSetBoundaryCondition = m_levelSetTestCase;
  /*! \page propertiesLS
    \section levelSetBoundaryCondition
    <code>MInt LsCartesianSolver::m_levelSetBoundaryCondition </code>\n
    default = <code>m_levelSetTestCase</code>\n \n
    This property triggers the use of special boundary conditions for each specific level set case, see also \ref
    levelSetTestCase \n possible values are: <ul> <li>1 - 2D bunsen flame </li> <li>1949 - ? </li> <li>1957 - ? </li>
    <li>1966,1967,1968 - ? </li>
    <li>17511,175110 - forced response </li>
    <li>17512 - forced response bunsen flame  </li>
    <li>17513,17514,17515,17516,17518 - forced response bunsen flame -> fixes the level set function on the flame tube
    edges, \ref radiusFlameTube, \ref xOffsetFlameTube(2) for one or two flames, \ref twoFlames </li> <li>17517 - forced
    response symmetric bunsen flame -> fixes the level set function on the flame tube edges, \ref radiusFlameTube, \ref
    xOffsetFlameTube </li>

    </ul>
    Keywords: <i>LEVELSET, SMOOTHED, EXTENSION, VELOCITY, BAND, LAYER, G0, ZERO, FRONT, CELLS </i>
  */
  m_levelSetBoundaryCondition =
      Context::getSolverProperty<MInt>("levelSetBoundaryCondition", m_solverId, AT_, &m_levelSetBoundaryCondition);

  m_levelSetBC = "SYMMETRIC";
  /*! \page propertiesLS
   \section levelSetBC
   <code>MString FvCartesianSolverPar::m_levelSetBC </code>\n
   default = <code>"SYMMETRIC"</code>\n \n
   This property triggers the use of simplifications in the GField-calculations for special testcases and setups,
   possible values are: <ul> <li>SYMMETRIC </li> <li>PERIODIC </li> <li>NONE </li>
   </ul>
   If the Levelset-Field is neither symmetric nor periodic, the NONE specification has to be selected!
   Keywords: <i>LEVELSET BAND, LAYER, G0, ZERO, FRONT, CELLS </i>
 */

  m_levelSetBC = Context::getSolverProperty<MString>("levelSetBC", m_solverId, AT_, &m_levelSetBC);

  // only relevant for m_levelSetBoundaryCondition = 1751600!
  m_noHaloLayers = 1;
  /*! \page propertiesLS
    \section noHaloLayers
    <code>MInt LsCartesianSolver::m_noHaloLayers </code>\n
    default = <code>m_noHaloLayers</code>\n \n
    This property sets the number of halo layers for the level set grid. \n
    possible values are:
    <ul>
    <li> Any non-negative integer number greater than one for parallel computation</li>
    </ul>
    Keywords: <i>LEVELSET, NUMBER, HALO, LAYERS, PARALLEL, COMPUTATION </i>
  */
  m_noHaloLayers = Context::getSolverProperty<MInt>("noHaloLayers", m_solverId, AT_, &m_noHaloLayers);

  if(m_fourthOrderNormalCurvatureComputation) {
    ASSERT(m_noHaloLayers > 1, "Not enough halo layer for fourth order curvature computation");
  }

  m_reinitInterval = 1;
  /*! \page propertiesLS
    \section reinitInterval
    <code>MInt LsCartesianSolver::m_reinitInterval </code>\n
    default = <code>1</code>\n \n
    This property triggers the reinitialization interval, but it doesn't ensures the reinitialization, this is depend on
    wether the thresholds, \ref reinitThreshold, \ref reinitThresholdAvg \n or could be ensured by the property \ref
    guaranteeReinit. \n possible values are: <ul> <li> Any non-negative integer number smaller than the maximum number
    of time steps \ref timeSteps</li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, INTERVAL </i>
  */
  m_reinitInterval = Context::getSolverProperty<MInt>("reinitInterval", m_solverId, AT_, &m_reinitInterval);

  /*! \page propertiesLS
    \section maintainOuterBandLayers
    <code>MBool LsCartesianSolver::m_maintainOuterBandLayers </code>\n
    default = <code>0</code>\n \n
    This property guarantees that the outer band layers are maintained, see also \ref reinitInterval and \ref
    guaranteeReinit. \n possible values are: <ul> <li>0 - off</li> <li>1 - on</li>
    </ul>
    Keywords: <i>LEVELSET, MAINTENANCE, GUARANTEE, OUTER, BAND, LAYERS </i>
  */
  m_maintainOuterBandLayers = Context::getSolverProperty<MBool>("maintainOuterBandLayers", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
    \section writeReinitializationStatistics
    <code>MInt LsCartesianSolver::m_writeReinitializationStatistics </code>\n
    default = <code>0</code>\n \n
    This property triggers the reinitialization statistics. If the statistics are turned on you get the ASCII files
    deviation and avgGradientAfter containing the statistics for the G0 cells (level set front cells): \n deviation:
    <ul>
    <li> \ref globalTimeStep</li>
    <li> deviation (?)</li>
    </ul>

    avgGradientAfter:
    <ul>
    <li> \ref globalTimeStep</li>
    <li> averaged deviation (here called: averaged gradient), see \ref reinitThresholdAvg </li>
    <li> variance  \f$ \sigma =  \overline{\Delta_{threshAvg}}^{2} \f$ </li>
    <li> mean level set function gradient  \f$  \overline{\left| \nabla G \right| }\f$ </li>
    <li> maximum level set function gradient  \f$  \left| \nabla G \right|_{MAX}\f$ </li>
    <li> minimum level set function gradient  \f$  \left| \nabla G \right|_{MIN}\f$ </li>
    </ul>

    The property is only meaningful for the renitialization methods, \ref reinitMethod,
    HCR1,HCR2,HCR2_LIMITED,HCR2_FULLREINIT. The statistics should be included also for the other reinitialization
    methods. \n possible values are: <ul> <li>0 - off </li> <li>1 - on </li>
    </ul>
    Keywords: <i>LEVELSET, REINITIALIZATION, STATISTICS, DEVIATION, AVGGRADIENTAFTER, G0, CELLS, GRADIENT </i>
  */
  m_writeReinitializationStatistics =
      Context::getSolverProperty<MBool>("writeReinitializationStatistics", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
    \section interpolateFlowFieldToFlameFront
    <code>MBool LsCartesianSolver::m_interpolateFlowFieldToFlameFront </code>\n
    default = <code>0</code>\n \n
    This property triggers the transfer of the velocity v from the flow to the G-grid and interpolates the flow field to
    the flame front, see \ref computeExtVel (cases 5..81) \n possible values are: <ul> <li>off - transfer is done but
    without interpolation (possible for all test cases) </li> <li>on - transfer is done with interpolations (probably
    not possible for all test cases) </li>
    </ul>
    Keywords: <i>LEVELSET, TRANSFER, FLOW, FIELD, INTERPOLATION </i>
  */
  m_interpolateFlowFieldToFlameFront =
      Context::getSolverProperty<MBool>("interpolateFlowFieldToFlameFront", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllLevelSetFunctions
 <code>MBool LsCartesianSolver::m_writeOutAllLevelSetFunctions </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set functions are written out (.vtk-file), or just the zeroth level-set function \n
 possible values are:
 <ul>
 <li>off - output of the zeroth level-set function only </li>
 <li>on - output of all level-set functions (file may become large!) </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllLevelSetFunctions =
      Context::getSolverProperty<MBool>("writeOutAllLevelSetFunctions", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllExtensionVelocities
 <code>MBool LsCartesianSolver::m_writeOutAllExtensionVelocities </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set extension velocities are written out (.vtk-file), or just respect to the zeroth
 level-set function \n possible values are: <ul> <li>off - output of the extension velocities with respecot to the
 zeroth level-set function only </li> <li>on - output of all level-set extension velocities (file may become large!)
 </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllExtensionVelocities =
      Context::getSolverProperty<MBool>("writeOutAllExtensionVelocities", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllCurvatures
 <code>MBool LsCartesianSolver::m_writeOutAllCurvatures </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set function curvatures are written out (.vtk-file), or just the zeroth level-set
 function curvature \n possible values are: <ul> <li>off - output of the zeroth level-set function curvature only </li>
 <li>on - output of all level-set function curvatures (file may become large!) </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllCurvatures = Context::getSolverProperty<MBool>("writeOutAllCurvatures", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllCorrectedBurningVelocity
 <code>MBool LsCartesianSolver::m_writeOutAllCorrectedBurningVelocity </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set function corrected burning velocities are written out (.vtk-file), or just the
 zeroth level-set function corrected burning velocity\n possible values are: <ul> <li>off - output of the zeroth
 level-set function corrected burning velocity only </li> <li>on - output of all level-set function corrected burning
 velocities (file may become large!) </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllCorrectedBurningVelocity =
      Context::getSolverProperty<MBool>("writeOutAllCorrectedBurningVelocity", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllFlameSpeeds
 <code>MBool LsCartesianSolver::m_writeOutAllFlameSpeeds </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set function flame speeds are written out (.vtk-file), or just the zeroth level-set
 function flame speed \n possible values are: <ul> <li>off - output of the zeroth level-set function flame speed only
 </li> <li>on - output of all level-set function flame speeds (file may become large!) </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllFlameSpeeds = Context::getSolverProperty<MBool>("writeOutAllFlameSpeeds", m_solverId, AT_, &tmpFalse);

  /*! \page propertiesLS
 \section writeOutAllNormalVectors
 <code>MBool LsCartesianSolver::m_writeOutAllNormalVectors </code>\n
 default = <code>0</code>\n \n
 This property triggers if all level-set function normal vectors are written out (.vtk-file) \n
 possible values are:
 <ul>
 <li>off - output of the normal vectors </li>
 <li>on - output of the normal vectors (file may become large!) </li>
 </ul>
 Keywords: <i>LEVELSET, MULTILEVELSET </i>
 */
  m_writeOutAllNormalVectors =
      Context::getSolverProperty<MBool>("writeOutAllNormalVectors", m_solverId, AT_, &tmpFalse);

  if(string2enum(m_reinitMethod) == HCR2_FULLREINIT) {
    // output info min iteration
    MInt minIteration = m_gBandWidth / m_reinitCFL;
    m_log << "reinitMethod: HCR2_FULLREINIT -> full reinitialization guaranteed after " << minIteration << " iterations"
          << endl;
  }

  if(m_combustion) {
    // compute minimum iterations to transport the velocity to the band cells
    MInt minIteration = m_gBandWidth / m_extVelCFL;
    m_log << "number of minimum iterations necessary to extend velocity field on band " << minIteration << endl;
    m_log << "number of iterations set to extend velocity field on band " << m_extVelIterations << endl;
  }

  /*! \page propertiesLS
    \section engineSetup
    <code>MBool LsCartesianSolver::m_engineSetup </code>\n
    Triggers specific stuff for TINA or other engine applications
    default = <code>false</code>\n \n
    Keywords: <i>LEVELSET, ADAPTATION</i>
  */
  m_engineSetup = Context::getSolverProperty<MBool>("engineSetup", m_solverId, AT_, &m_engineSetup);

  if(m_LsRotate) {
    mAlloc(m_outerBandWidth, maxRefinementLevel() + 1, "m_outerBandWidth", 0, AT_);
    m_outerBandWidth[maxRefinementLevel()] = m_gShadowWidth;
    for(MInt lvl = maxRefinementLevel() - 1; lvl >= 0; lvl--) {
      MInt dif = (maxRefinementLevel() - lvl);
      m_outerBandWidth[lvl] = (m_outerBandWidth[lvl + 1] / 2) + 1 + dif;
    }
  } else if(m_combustion) {
    mAlloc(m_outerBandWidth, maxRefinementLevel() + 1, "m_outerBandWidth", 0, AT_);
    m_outerBandWidth[maxRefinementLevel()] = m_gShadowWidth;
    for(MInt lvl = maxRefinementLevel() - 1; lvl >= 0; lvl--) {
      MInt factor = IPOW2((maxRefinementLevel() - lvl));
      MInt currentWidth = m_gShadowWidth / factor;
      if(currentWidth * factor < m_gShadowWidth) {
        currentWidth++;
      }
      m_outerBandWidth[lvl] = currentWidth;
    }
  } else if(m_engineSetup) {
    mAlloc(m_outerBandWidth, maxRefinementLevel() + 1, "m_outerBandWidth", 0, AT_);
    m_outerBandWidth[maxRefinementLevel()] = m_gShadowWidth + 2;
    for(MInt lvl = maxRefinementLevel() - 1; lvl >= 0; lvl--) {
      MInt val = m_outerBandWidth[lvl + 1] / 2;
      if(val * 2 < m_outerBandWidth[lvl + 1]) {
        val++;
      }
      m_outerBandWidth[lvl] = val;
    }
  } else if(m_levelSetRans) {
    mAlloc(m_outerBandWidth, maxRefinementLevel(), "m_outerBandWidth", 0, AT_);
    m_outerBandWidth[maxRefinementLevel() - 1] = m_gShadowWidthRans;
    for(MInt i = maxRefinementLevel() - 2; i >= 0; i--) {
      m_outerBandWidth[i] = (m_outerBandWidth[i + 1] / 2) + 1;
    }
  } else {
    mAlloc(m_outerBandWidth, maxRefinementLevel(), "m_outerBandWidth", 0, AT_);
    m_outerBandWidth[maxRefinementLevel() - 1] = m_gShadowWidth;
    for(MInt i = maxRefinementLevel() - 2; i >= 0; i--) {
      m_outerBandWidth[i] = (m_outerBandWidth[i + 1] / 2) + 1;
    }
  }

  m_refineDiagonals = true;
  /*! \page propertiesLS
    \section refineDiagonals
    <code>MBool LsCartesianSolver::m_refineDiagonals </code>\n
    default = <code>true</code>\n \n
    Determines whether the diagonal cells for the levelSet band should be refined as well!
    <ul>
    <li>1.0 +- eps</li>
    </ul>
    Keywords: <i>LEVELSET, ADAPTATION</i>
  */
  if(m_LsRotate) m_refineDiagonals = false;
  m_refineDiagonals = Context::getSolverProperty<MBool>("refineDiagonals", m_solverId, AT_, &m_refineDiagonals);

  m_bodyIdOutput = m_semiLagrange;
  /*! \page propertiesLS
    \section bodyIdOutput
    <code>MBool LsCartesianSolver::m_bodyIdOutput </code>\n
    default = <code>true</code>\n \n
    Determines whether the additionale bodyIdOutput output should be written!
    Default is true for semiLagrange, but the output is not necessary for a restart and
    should be set to false for large applications, to reduce the amount of data!
    Keywords: <i>LEVELSET, SEMI-LAGRANGE</i>
  */
  m_bodyIdOutput = Context::getSolverProperty<MBool>("bodyIdOutput", m_solverId, AT_, &m_bodyIdOutput);


  if(m_noInitGFieldFromSTLBndCndIds + m_noBodyBndryCndIds > 0) {
    mAlloc(m_initGFieldFromSTLBndCndIds, m_noInitGFieldFromSTLBndCndIds + m_noBodyBndryCndIds,
           "m_initGFieldFromSTLBndCndIds", -1, AT_);
  }
  MInt defaultValue = -1;
  for(MInt i = 0; i < m_noBodyBndryCndIds; i++) {
    m_initGFieldFromSTLBndCndIds[i] = m_bodyBndryCndIds[i];
  }
  m_noInitGFieldFromSTLBndCndIds += m_noBodyBndryCndIds;
  for(MInt i = m_noBodyBndryCndIds; i < m_noInitGFieldFromSTLBndCndIds; i++) {
    m_initGFieldFromSTLBndCndIds[i] = Context::getSolverProperty("GFieldInitFromSTLBndCndIds", m_solverId, AT_,
                                                                 &defaultValue, i - m_noBodyBndryCndIds);
  }


  /*! \page propertiesLS
    \section LsGeometryChange
    <code>MBool LsCartesianSolver::m_geometryChange </code>\n
    default = <code>false</code>\n \n
    Used to trigger an adaptation right after the restart, which can be used to change
    the geometry within the ls-BandWidth. Useful for setup changes.
    Keywords: <i>LEVELSET, SEMI-LAGRANGE, GEOMETRY</i>
  */

  if(Context::propertyExists("LsGeometryChange", m_solverId)) {
    m_forceAdaptation = true;

    mAlloc(m_geometryChange, m_noSets, "geometryChange", false, AT_);

    MBool anychange = false;
    for(MInt set = 0; set < m_noSets; set++) {
      m_geometryChange[set] =
          Context::getSolverProperty<MBool>("LsGeometryChange", m_solverId, AT_, &m_geometryChange[set], set);
      if(m_geometryChange[set]) {
        anychange = true;
        if(domainId() == 0) {
          cerr << "Changing geometry of set " << set << " upon restart!" << endl;
        }
      }
    }

    if(!anychange) {
      mDeallocate(m_geometryChange);
      ASSERT(m_geometryChange == nullptr, "");
    }
  }

  m_periodicMovement = false;
  if(Context::propertyExists("periodicMovement", m_solverId)) {
    m_periodicMovement = Context::getSolverProperty<MBool>("periodicMovement", m_solverId, AT_, &m_periodicMovement);
    m_periodicDirection = Context::getSolverProperty<MInt>("periodicDirection", m_solverId, AT_, &m_periodicDirection);
  }
  /*! \page propertiesLS
   \section weightBaseCell
   <code>MFloat LS::m_weightBaseCell</code>\n
   default = <code>1.0</code>\n \n
   Weight applied for any level-set cell during static weight computation for
   domain decomposition during balance, good value could be 0.01
   Keywords: <i> LS, WEIGHTING, BALANCE</i>
 */
  m_weightBaseCell = 1.0;
  m_weightBaseCell = Context::getSolverProperty<MFloat>("weightBaseCell", solverId(), AT_, &m_weightBaseCell);

  /*! \page propertiesLS
    \section weightLeafCell
    <code>MFloat LS::m_weightLeafCell</code>\n
    default = <code>1.0</code>\n \n
    Weight applied for any level-set leaf-cell during static weight computation
    for domain decomposition during balance, good value could be 0.05
    Keywords: <i> LS, WEIGHTING, BALANCE</i>
  */
  m_weightLeafCell = 1.0;
  m_weightLeafCell = Context::getSolverProperty<MFloat>("weightLeafCell", solverId(), AT_, &m_weightLeafCell);

  /*! \page propertiesLS
    \section weightBandCell
    <code>MFloat LS::m_weightBandCell</code>\n
    default = <code>1.0</code>\n \n
    Weight applied for any level-set leaf-cell band cell during static weight computation
    for domain decomposition during balance, good value could be 0.1
    Keywords: <i> LS, WEIGHTING, BALANCE</i>
  */
  m_weightBandCell = 1.0;
  m_weightBandCell = Context::getSolverProperty<MFloat>("weightBandCell", solverId(), AT_, &m_weightBandCell);

  /*! \page propertiesLS
    \section weightMulitSolverFactor
    <code>MFloat LS::m_weightMulitSolverFactor</code>\n
    default = <code>1.0</code>\n \n
    Mutli-solver weight factor applied to all level-set cell weights for static weight
    computation for domain decomposition during balance. 1.0 for single solver application
    , otherwise setup dependent.
    Keywords: <i> LS, WEIGHTING, BALANCE</i>
  */
  m_weightMulitSolverFactor = 1.0;
  m_weightMulitSolverFactor =
      Context::getSolverProperty<MFloat>("weightMulitSolverFactor", solverId(), AT_, &m_weightMulitSolverFactor);

  /*! \page propertiesLS
    \section limitWeights
    <code>MBool LS::m_limitWeights</code>\n
    default = <code>false</code>\n \n
    Limit weight of level-set cells by a factor of the largest weight,
    to ensure a more even distribution of solver memory across ranks.
    Keywords: <i> LS, WEIGHTING, BALANCE</i>
  */
  m_limitWeights = false;
  m_limitWeights = Context::getSolverProperty<MBool>("limitDLBWeights", solverId(), AT_, &m_limitWeights);

  /*! \page propertiesLS
    \section weightLevelSet
    <code>MBool LS::m_weightLevelSet</code>\n
    default = <code>true</code>\n \n
    Triggers if level-set cells should be considered as DLB weights.
    Can be deactivated if the level-set is not updated each time-step
    can no overhead exists i.e., for static boundaries.
    Keywords: <i> LS, WEIGHTING, BALANCE</i>
  */
  m_weightLevelSet = true;
  m_weightLevelSet = Context::getSolverProperty<MBool>("weightLevelSet", solverId(), AT_, &m_weightLevelSet);
}

/// Stuff to be done before the timestep.
/// \author Sven Berger, Update: Tim Wegmann
/// \tparam nDim
template <MInt nDim>
void LsCartesianSolver<nDim>::preTimeStep() {
  TRACE();

  if(m_combustion || m_freeSurface) {
    exchangeAllLevelSetData();

    setGCellBndryProperty();

    updateBndryCellList();

    if(string2enum(m_levelSetBC) == PERIODIC) {
      computeNormalVectorsPeriodic();

      computeCurvaturePeriodic();
    } else {
      computeNormalVectors();

      computeCurvature();
    }
    // constructExtensionVelocity(); //this is now done in the coupler
  } else if(!m_levelSetMb) {
    setGCellBndryProperty();

    if(!m_semiLagrange) {
      computeNormalVectors();

      computeCurvature();

      determinePropagationSpeed();
    }

  } else if(m_levelSetMb) {
    if(m_constructGField != 0) return;

    checkHaloCells();

    static MBool firstRun = true;

    if(firstRun) {
      for(MInt set = 0; set < m_noSets; set++) {
        m_computeSet[set] = m_computeSet_tmp[set];
        m_changedSet[set] = true;
      }
      firstRun = false;
    }

    for(MInt set = 0; set < m_noSets; set++) {
      ASSERT(m_computeSet[set] == m_computeSet_backup[set], "ERROR in m_computeSet");
    }

    m_time += timeStep();
  }
}

/**
 * \brief  handle multiple levelsets
 *         build the collected levelset
 *         mode = 1: (default) called after each timeStep
 *         mode = 2: (for adaptation and balance) calling buildCollectedLevelSet with mode 2
 *                   without the forward step in gapWidth
 * \author Tim Wegmann
 */
///
/// \author Sven Berger
/// \tparam nDim
template <MInt nDim>
void LsCartesianSolver<nDim>::buildMultipleLevelSet(const MInt mode) {
  TRACE();

  if(!m_buildCollectedLevelSetFunction) return;

  NEW_TIMER_GROUP_STATIC(buildMultipleLevelSet, "MultipleLevelset");
  NEW_TIMER_STATIC(t_multiLevelSet, "Total time - multiple-levelset", buildMultipleLevelSet);
  NEW_SUB_TIMER_STATIC(t_c1, "identifyBodies", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c2, "collectedLevelSet", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c3, "gapHandling", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c4, "gapCells", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c5, "gapReInit", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c6, "Tube0", t_multiLevelSet);
  NEW_SUB_TIMER_STATIC(t_c7, "curv+normal", t_multiLevelSet);

  RECORD_TIMER_START(t_multiLevelSet);


  if(mode == 2) {
    // initialisation version
    identifyBodies(0);
  } else {
    // TODO labels:LS,TIMERS update default to faster shifting version after further testing!
    RECORD_TIMER_START(t_c1);
    identifyBodies();
    RECORD_TIMER_STOP(t_c1);
  }

  RECORD_TIMER_START(t_c2);
  // generares the coorect collected levelSet!
  // the bodies need to be identified before!
  buildCollectedLevelSet(mode);

  buildLevelSetTube(0);
  RECORD_TIMER_STOP(t_c2);

  RECORD_TIMER_START(t_c3);
  gapHandling();
  RECORD_TIMER_STOP(t_c3);

  RECORD_TIMER_START(t_c4);
  // if gaps exist, solve this
  if(m_closeGaps && gapCellsExist()) {
    // levelSet value reduction for G0-Band cells
    levelSetGapCorrect();

    RECORD_TIMER_START(t_c6);
    buildLevelSetTube(0);
    RECORD_TIMER_STOP(t_c6);

    setBandNewArrivals(0);

    RECORD_TIMER_START(t_c7);
    if(m_guaranteeReinit) {
      computeNormalVectors(0);
      computeCurvature(0);
    }
    RECORD_TIMER_STOP(t_c7);

    RECORD_TIMER_START(t_c5);
    if(m_guaranteeReinit) {
      levelSetReinitialization(0);
    }
    RECORD_TIMER_STOP(t_c5);

    // levelSet value increase for G0-Band cells
    levelSetGapRecorrect();

    reBuildCollectedLevelSet(1);

    RECORD_TIMER_START(t_c6);
    buildLevelSetTube(0);
    RECORD_TIMER_STOP(t_c6);

    setBandNewArrivals(0);
  }

  RECORD_TIMER_STOP(t_c4);

  exchangeAllLevelSetData();

  if(!m_semiLagrange) {
    computeNormalVectors(0);
    computeCurvature(0);
    exchangeAllLevelSetData();
  }

  RECORD_TIMER_STOP(t_multiLevelSet);
}


/// Solution of the levelset solver
/// \author Sven Berger
/// \tparam nDim
template <MInt nDim>
MBool LsCartesianSolver<nDim>::_levelSetSolutionStep() {
  TRACE();

  if(m_freeSurface) {
    updateBndryCellList();
    return levelSetSolver();

  } else if(m_combustion) {
    updateBndryCellList();
    applyLevelSetBoundaryConditions();
    return levelSetSolver();

  } else if(!m_levelSetMb && m_semiLagrange) {
    applyLevelSetBoundaryConditions();
    return semiLagrangeTimeStep();

  } else if(!m_levelSetMb) {
    applyLevelSetBoundaryConditions();
    return levelSetSolver();

  } else if(m_levelSetMb) {
    ASSERT(!m_constructGField, "");
    ASSERT(m_semiLagrange, "");
    return semiLagrangeTimeStep();
  }

  return true;
}


/// Rebuilds the grid and reinitializes (this is the 2nd verion used by the combustion code)
/// \author Sven Berger
/// \tparam nDim

template <MInt nDim>
MBool LsCartesianSolver<nDim>::finalizeLevelSet_(const MInt t_levelSet, const MInt t_output) {
  TRACE();

  // Necessary to avoid 'unused-parameter' warning when timers are disabled
  static_cast<void>(t_levelSet);
  static_cast<void>(t_output);

  NEW_SUB_TIMER_STATIC(t_levelSetGGrid, "GGrid", t_levelSet);
  NEW_SUB_TIMER_STATIC(t_levelSetFastGGrid, "FastGGrid", t_levelSet);
  NEW_SUB_TIMER_STATIC(t_levelSetCR, "CR", t_levelSet);

  RECORD_TIMER_START(t_levelSetGGrid);
  buildLevelSetTube();
  setBandNewArrivals();
  RECORD_TIMER_STOP(t_levelSetGGrid);

  // reinitializes the level set function
  RECORD_TIMER_START(t_levelSetCR);
  levelSetReinitialization();
  RECORD_TIMER_STOP(t_levelSetCR);

  RECORD_TIMER_START(t_levelSetFastGGrid);
  fastBuildLevelSetTubeCG();
  RECORD_TIMER_STOP(t_levelSetFastGGrid);

  if(m_combustion) {
    determineMinMaxMeanInterfacePosition();
  }

  return true;
}
/// Solution of the levelset solver
/// \author Sven Berger
/// \tparam nDim
template <MInt nDim>
MBool LsCartesianSolver<nDim>::solutionStep() {
  TRACE();

  if(m_constructGField) return true;

  if(!m_trackMovingBndry || globalTimeStep < m_trackMbStart || globalTimeStep > m_trackMbEnd) {
    return true;
  }

  RECORD_TIMER_START(m_timers[Timers::TimeInt]);
  m_firstSolutionExchange = true;
  MBool timeStepCompleted = _levelSetSolutionStep();
  exchangeLeafDataLS<true>();
  RECORD_TIMER_STOP(m_timers[Timers::TimeInt]);

  return timeStepCompleted;
}


/// Rebuilds the grid and reinitializes
/// \author Sven Berger
/// \tparam nDim
template <MInt nDim>
void LsCartesianSolver<nDim>::finalizeLevelSet() {
  TRACE();

  RECORD_TIMER_START(m_timers[Timers::Finalize]);

  if(m_combustion || m_freeSurface) {
    buildLevelSetTube();
    setBandNewArrivals();

    // rebuilds the tube and reinitializes the level set function
    levelSetReinitialization();

  } else if(m_semiLagrange) {
    if(m_constructGField) {
      return;
    }

    testCellsCG();

    // for the individual sets
    RECORD_TIMER_START(m_timers[Timers::BuildTube]);
    buildLevelSetTube();
    RECORD_TIMER_STOP(m_timers[Timers::BuildTube]);

    RECORD_TIMER_START(m_timers[Timers::SetBand]);
    setBandNewArrivals();
    RECORD_TIMER_STOP(m_timers[Timers::SetBand]);

    // for the collected levelset
    RECORD_TIMER_START(m_timers[Timers::BuildMultiple]);
    buildMultipleLevelSet();
    RECORD_TIMER_STOP(m_timers[Timers::BuildMultiple]);

  } else {
    testCellsCG();

    buildLevelSetTube();

    setBandNewArrivals();

    computeNormalVectors();

    computeCurvature();

    levelSetReinitialization();

    buildMultipleLevelSet();
  }

  RECORD_TIMER_STOP(m_timers[Timers::Finalize]);
}

/**
 * \
 * \author Jannik Borgelt
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::determinePeriodicDistance() {
  TRACE();

  // Determine local minimum ans maximum
  MFloat minPos = std::numeric_limits<MFloat>::max();
  MFloat maxPos = std::numeric_limits<MFloat>::lowest();
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    MInt level = a_level(cellId);
    MFloat cellHalfLength = 0.5 * c_cellLengthAtLevel(level);
    if(c_coordinate(cellId, m_periodicDirection) - cellHalfLength < minPos) {
      minPos = c_coordinate(cellId, m_periodicDirection) - cellHalfLength;
    }
    if(c_coordinate(cellId, m_periodicDirection) + cellHalfLength > maxPos) {
      maxPos = c_coordinate(cellId, m_periodicDirection) + cellHalfLength;
    }
  }

  // Exchange global minimum and maximum
  MPI_Allreduce(MPI_IN_PLACE, &maxPos, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "maxPos", "1");
  MPI_Allreduce(MPI_IN_PLACE, &minPos, 1, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "minPos", "1");

  m_periodicDistance = maxPos - minPos;
}

/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::initLocalizedLevelSetCG() {
  TRACE();
  MInt noCells;

  m_log << "Initializing local levelset CG" << endl;
  if(domainId() == 0) {
    cerr << "Initializing local levelset CG" << endl;
  }

  m_time = 0.0;

  // initialise minGapWidth and minGapWidthDt1#
  if(m_closeGaps) {
    m_minGapWidth.clear();
    m_minGapWidthDt1.clear();
    for(MInt region = 0; region < m_noGapRegions; region++) {
      m_minGapWidth.push_back(std::numeric_limits<MFloat>::max());
      m_minGapWidthDt1.push_back(std::numeric_limits<MFloat>::max());
    }
  }

  if(m_GFieldInitFromSTL) initializeGControlPoint();

  for(MInt set = 0; set < m_noSets; set++)
    std::vector<MInt>().swap(m_bandCells[set]);

  if(m_semiLagrange) {
    initializeIntegrationScheme_semiLagrange();
  } else {
    initializeIntegrationScheme();
  }

  if(!m_initFromRestartFile) {
    // Iteration 1: (on minLevel-/boundingBoxLevel-basis)
    //------------------------------------------------
    // a) create G-Base-Grid
    // if( domainId()==0 ) { cerr << "G grid initialization: Iteration-1" << endl; }

    createBaseGgridCG();

    // c) initialize the G-field on minLevel-/boundingBoxLevel-basis
    if(m_GFieldInitFromSTL) {
      constructGFieldFromSTL(1);
    } else {
      initializeGField();
    }

    // d) determine G-cell-properties
    if(m_buildCollectedLevelSetFunction) {
      buildCollectedLevelSet(0);
    }

    determineG0Cells();

    // Iteration 4:
    //------------------------------------------------
    // if( domainId()==0 ) { cerr << "G grid initialization: Iteration-4" << endl; }


    // a) refine or coarsen g0 cells!
    m_log << "createGgridCG is starting" << endl;
    createGgridCG();
    m_log << "createGgridCG is done" << endl;

    // b) initialize again
    // if( m_GFieldInitFromSTL ) constructGFieldFromSTL(1);
    if(m_GFieldInitFromSTL) {
      determineG0Cells();
      determineBandCells();
      m_log << "Initialize G Field from stl data...";
      constructGFieldFromSTL(3);
      m_log << "ok" << endl;
      for(MInt set = 0; set < m_noSets; set++) {
        std::vector<MInt>().swap(m_bandCells[set]);
      }
    } else {
      initializeGField();
    }

    // c) determine G-cell-properties
    if(m_buildCollectedLevelSetFunction) {
      buildCollectedLevelSet(0);
    }

    m_log << "Determine G cells in Gamma...";
    determineG0Cells();
    m_log << "ok" << endl;

    m_log << "Determine the narrow computation band...";
    determineBandCells();
    m_log << "ok" << endl;

    m_log << "Update G boundary cell list...";
    updateBndryCellList();
    m_log << "ok" << endl;

    m_log << "smallest gCell length " << m_gCellDistance << endl;

    m_log << "Reset outside cells...";
    resetOutsideCells();
    m_log << "ok" << endl;

    if(domainId() == 0) {
      cerr << "G grid Reinitialization " << endl;
    }

    // moving boundary
    if(m_GFieldInitFromSTL) {
      constructGFieldFromSTL(0);
    }

    // collected level-set function
    if(m_buildCollectedLevelSetFunction) {
      buildCollectedLevelSet(0);
    }

    // claudia: check if this is really necessary here...
    if(m_GFieldInitFromSTL) {
      resetOutsideCells();
    }

    if(m_levelSetRans) {
      constructGFieldFromSTL(1);
      updateAllLowerGridLevels();
      // initialization output for checking purpose
      if(m_GFieldFromSTLInitCheck) {
        writeRestartLevelSetFileCG(1, "restartLSGridCG_init", "restartLSCG_init");
      }
      return;
    }

// usefull-debug-output:
#ifdef LS_DEBUG
    MInt globalTimeStep_temp = globalTimeStep;
    globalTimeStep = 0;
    writeRestartLevelSetFileCG(1, "restartLSGridCG_debug_0", "restartLSCG_debug_0");
    globalTimeStep = globalTimeStep_temp;
#endif

    // if(m_initialRefinement) return;

    // initialize fields...
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      for(MInt set = 0; set < m_maxNoSets; set++) {
        if(!m_semiLagrange) {
          for(MInt d = 0; d < nDim; d++) {
            a_extensionVelocityG(cellId, d, set) = F0;
          }
        } else {
          a_oldLevelSetFunctionG(cellId, set) = a_levelSetFunctionG(cellId, set);
        }
      }
      if(m_maxNoSets > 1 && m_closeGaps) {
        a_potentialGapCell(cellId) = 0;
        a_potentialGapCellClose(cellId) = 0;
      }
    }


    // initialization output for checking purpose
    if(m_GFieldFromSTLInitCheck) {
      writeRestartLevelSetFileCG(1, "restartLSGridCG_init", "restartLSCG_init");
    }

    exchangeAllLevelSetData();

    if(m_virtualSurgery) {
      initializeCollectedLevelSet(0);
    }

    m_log << "*** G grid initialization finished ***" << endl;

  } else { // initFromRestartFile

    if(domainId() == 0) {
      cerr << " By loading localized level-set file: ";
    }
    stringstream levelSetFileName;
    levelSetFileName << restartDir() << "restartLSCG_init";
    if(!m_useNonSpecifiedRestartFile) levelSetFileName << "_0";
    levelSetFileName << ParallelIo::fileExt();


    if(domainId() == 0) {
      cerr << (levelSetFileName.str()).c_str() << endl;
    }

    noCells = loadLevelSetGridFlowVarsParCG((levelSetFileName.str()).c_str());

    ASSERT(a_noCells() == noCells, " size of collector is not noCells: " << a_noCells() << " noCells: " << noCells);

    generateListOfGExchangeCellsCG();

    for(MInt set = 0; set < m_noSets; set++) {
      std::vector<MInt>().swap(m_bandCells[set]);
    }

    determineG0Cells();

    createGgridCG();

    determineG0Cells();

    // set g_properties[0,1]
    determineBandCells();

    setGCellBndryProperty();

    if(m_levelSetRans) {
      updateAllLowerGridLevels();
    }

    updateBndryCellList();

    m_log << "Reset outside cells...";
    resetOutsideCells();
    m_log << "ok" << endl;

    // initialize fields...
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(m_maxNoSets > 1 && m_closeGaps) {
        a_potentialGapCell(cellId) = 0;
        a_potentialGapCellClose(cellId) = 0;
      }
    }

    if(m_combustion) {
      determineMinMaxMeanInterfacePosition();
      computeZeroLevelSetArcLength();
    }

    if(m_initialRefinement) return;

    if(m_GFieldFromSTLInitCheck) {
      writeRestartLevelSetFileCG(1, "restartLSGridCG_reinit", "restartLSCG_reinit");
    }

    exchangeAllLevelSetData();
  }
}


/**
 * \brief
 * \author Sohel Herff, Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::createBaseGgridCG() {
  TRACE();

  MBool& firstRun = m_static_createBaseGgrid_firstRun;

  // assure that the function is only called once,
  // otherwise cells are added multiple-times!
  if(!firstRun) return;
  firstRun = false;

  m_log << "createBaseGgridCG" << endl;

  if(minLevel() == 0) {
    mTerm(1, AT_, "minLevel() is 0");
  }

  // establish a base g cell grid
  for(MInt fc = 0; fc < grid().tree().size(); fc++) {
    a_appendCollector();
    m_cells.erase(a_noCells() - 1);
  }

  generateListOfGExchangeCellsCG();
}


/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
// this function does some sanity assert checks
template <MInt nDim>
void LsCartesianSolver<nDim>::testCellsCG() {
  TRACE();

  if(!g_multiSolverGrid) {
    for(MInt gridCellId = 0; gridCellId < grid().raw().treeb().size(); gridCellId++) {
      ASSERT(grid().tree().solver2grid(gridCellId) == gridCellId, "");
      ASSERT(grid().tree().grid2solver(gridCellId) == gridCellId, "");
    }
  }

  for(MInt gc = 0; gc < a_noCells(); gc++) {
    ASSERT(grid().tree().grid2solver(grid().tree().solver2grid(gc)) == gc,
           "proxy-grid map is incorrect: gc, grid().tree().solver2grid(gc), "
           "grid().tree().grid2solver(grid().tree().solver2grid(gc)): "
               << gc << " , " << grid().tree().solver2grid(gc) << " , "
               << grid().tree().grid2solver(grid().tree().solver2grid(gc)));
  }
}

/**
 * \brief  set a_regridTriggerG
 *
 * default: initRegrid = false
 *
 * \author Sohel Herff & Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::createGgridCG(MBool initRegrid) {
  TRACE();

  m_log << "Reset lsSolver-regrid trigger at timestep: " << globalTimeStep << " init regrid: " << initRegrid << endl;

  vector<MInt> lastLayer;
  vector<MInt> tmp;

  // MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  // MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");

  MBoolScratchSpace inShadowLayer(a_noCells(), AT_, "inShadowLayer");
  MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

  MBool newVersion = false;
  MInt startSet = 0;
  if(newVersion) startSet = m_startSet;

  // 1) reset regrid property for all cells
  //    set inShadowLayer property for all G0-cells
  for(MInt gc = 0; gc < a_noCells(); gc++) {
    MBool isGZero = false;
    for(MInt set = startSet; set < m_noSets; set++) {
      if(/*newVersion &&*/ !m_computeSet_backup[set]) continue;
      if(a_isGZeroCell(gc, set)) {
        isGZero = true;
        break;
      }
    }
    inShadowLayer[gc] = isGZero;
    a_regridTriggerG(gc) = false;
  }


  for(MInt set = startSet; set < m_noSets; set++) { // m_startSet
    if(newVersion && !m_computeSet_backup[set]) continue;
    for(MInt id = 0; id < a_noG0Cells(set); id++) {
      const MInt currentCellId = a_G0CellId(id, set);
      MBool alreadyAdded = false;
      // loop over prevoius sets and check if the cell was already a G0 cell there
      // to avaoid double entries!
      if(set > 0) {
        for(MInt s = set - 1; s >= startSet; s--) { // m_startSet
          if(newVersion && !m_computeSet_backup[s]) continue;
          if(a_isGZeroCell(currentCellId, s)) {
            alreadyAdded = true;
            break;
          }
        }
      }
      if(!alreadyAdded) {
        lastLayer.emplace_back(currentCellId);
      }
    }
  }

  // 4. cover shadow band locally in the domain
  MInt itCount = 0;
  while(itCount < m_gShadowWidth) {
    // ii. build the next layer
    tmp.clear();
    for(MInt c = 0; c < (signed)lastLayer.size(); c++) {
      const MInt cellId = lastLayer[c];
      for(MInt d = 0; d < m_noDirs; d++) {
        if(a_hasNeighbor(cellId, d) == 0) continue;
        const MInt nghbrId = c_neighborId(cellId, d);
        if(!inShadowLayer[nghbrId]) {
          tmp.emplace_back(nghbrId);
          if(itCount > m_gShadowWidth - m_gInnerBound) {
            a_regridTriggerG(nghbrId) = true;
          }
        }
        inShadowLayer[nghbrId] = true;
      }
    }

    // exchange next layer with neighboring domains
    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noWindowCells(i); j++) {
        tmp_data(windowCellId(i, j)) = inShadowLayer(windowCellId(i, j));
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        tmp_data(windowId) = inShadowLayer(windowId);
      }
    }

    exchangeDataLS(&tmp_data(0), 1);

    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
      // check, if received data size matches expected size:
      if(!(receiveBufferSize.p[i] == noHaloCells(i))) {
        stringstream errorMessage{};
        errorMessage << "this was not expected to happen: wrong number of halo information, buf="
                     << receiveBufferSize.p[i] << "noGHaloCells=" << noHaloCells(i) << ", shadow layer";
        m_log << errorMessage.str() << endl;
        cerr << errorMessage.str() << endl;
      }
#endif
      for(MInt j = 0; j < noHaloCells(i); j++) {
        const MInt haloCell = haloCellId(i, j);
        if(!inShadowLayer[haloCell]) {
          inShadowLayer[haloCell] = tmp_data(haloCell);
          if(inShadowLayer[haloCell]) {
            tmp.emplace_back(haloCell);
            if(itCount > m_gShadowWidth - m_gInnerBound) {
              a_regridTriggerG(haloCell) = true;
            }
          }
        }
      }
    }
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        const MInt haloCell = grid().azimuthalHaloCell(i, j);
        if(!inShadowLayer[haloCell]) {
          inShadowLayer[haloCell] = tmp_data(haloCell);
          if(inShadowLayer[haloCell]) {
            tmp.emplace_back(haloCell);
            if(itCount > m_gShadowWidth - m_gInnerBound) {
              a_regridTriggerG(haloCell) = true;
            }
          }
        }
      }
    }

    // iii. swap data
    //      continue to the next layer
    lastLayer.clear();
    for(MInt c = 0; c < (signed)tmp.size(); c++) {
      lastLayer.emplace_back(tmp[c]);
    }
    itCount++;
  }

  m_log << itCount << " layers built - ";

  if(initRegrid) return;

  // 4. update G on all coarse grid cells
  // assign the value of the first child cell to the parentId that exists
  // TODO labels:LS delete this update of lower-gridLevels!
  //      once updateLowerGridLevels is called in prepareRestart
  //      otherwise this changes the testcases on lower grid levels
  for(MInt level = maxRefinementLevel() - 1; level > minLevel() - 1; level--) {
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_level(cellId) == level) {
        for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
          if(c_childId(cellId, ch) < 0) continue;
          for(MInt set = 0; set < m_noSets; set++) {
            if(level >= a_maxGCellLevel(set)) continue;
            a_levelSetFunctionG(cellId, set) = a_levelSetFunctionG(c_childId(cellId, ch), set);
            if(m_semiLagrange) {
              a_oldLevelSetFunctionG(cellId, set) = a_oldLevelSetFunctionG(c_childId(cellId, ch), set);
            }
          }
          break;
        }
      }
    }
  }

  // set bodyId for multiple sets but without collected Levelset!
  if(!m_buildCollectedLevelSetFunction && !m_combustion && !m_freeSurface) {
    if(m_maxNoSets > 1 || m_GFieldInitFromSTL) {
      MInt body = -1;
      for(MInt set = 0; set < m_noSets; set++) {
        if(m_noBodiesInSet[set] > 0) {
          body = m_setToBodiesTable[set][0];
        }
        for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
          a_bodyIdG(cellId, set) = body;
        }
      }
    } else {
      for(MInt set = 0; set < m_noSets; set++) {
        for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
          a_bodyIdG(cellId, set) = 0;
        }
      }
    }
  }
}


/**
 * \brief
 * \author Sohel Herff, Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::generateListOfGExchangeCellsCG() {
  TRACE();

  // Initialize Halo- and Window-Cells
  for(MInt gCell = 0; gCell < a_noCells(); gCell++) {
    a_isHalo(gCell) = false;
    a_isWindow(gCell) = false;
  }

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      a_isHalo(haloCellId(i, j)) = true;
    }
    for(MInt j = 0; j < noWindowCells(i); j++) {
      a_isWindow(windowCellId(i, j)) = true;
    }
  }

  // Mark azimuthal periodic exchange halos
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      a_isHalo(grid().azimuthalHaloCell(i, j)) = true;
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalUnmappedHaloCells(); i++) {
    a_isHalo(grid().azimuthalUnmappedHaloCell(i)) = true;
  }
}

/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 *
 *
 * mode = -1 : (default) check for regridding and use all sets
 * mode =  0 : ignore regridding and only use the zero-levelset!
 *
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::buildLevelSetTube(MInt mode) {
  TRACE();

  // timer-for mode 1:
  NEW_TIMER_GROUP_STATIC(buildLevelSet, "Level Set - build tube");
  NEW_TIMER_STATIC(t_levelSetTube, "Total time - build tube", buildLevelSet);
  NEW_SUB_TIMER_STATIC(t_a1, "updateLowerGridLevels", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a2, "determineG0Cells", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a3, "determineBandCells", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a4, "regridCheck", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a5, "createGGrid", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a6, "group", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a7, "updateBnd", t_levelSetTube);
  NEW_SUB_TIMER_STATIC(t_a8, "resetOutside", t_levelSetTube);

  // time for mode 0:
  NEW_TIMER_GROUP_STATIC(buildLevelSetZero, "Zero Level Set - build tube");
  NEW_TIMER_STATIC(t_levelSetTubeZero, "Total time - build tube zero", buildLevelSetZero);
  NEW_SUB_TIMER_STATIC(t_b1, "determineG0Cells", t_levelSetTubeZero);
  NEW_SUB_TIMER_STATIC(t_b2, "determineBandCells", t_levelSetTubeZero);
  NEW_SUB_TIMER_STATIC(t_b3, "resetOutside", t_levelSetTubeZero);

  if(mode == 0) {
    RECORD_TIMER_START(t_levelSetTubeZero);
  } else {
    RECORD_TIMER_START(t_levelSetTube);
  }


  // if in the complete global domain, there are no G0 cells present, currently no LVS is computed...
  MInt status = 1;
  if(mode != 0) {
    status = mMin(1, a_noG0Cells(0));
    for(MInt set = 1; set < m_noSets; set++) {
      status = mMin(status, a_noG0Cells(set));
    }
  }

  if(mode != 0) {
    RECORD_TIMER_START(t_a1);
  }

  // TODO labels:LS uncomment the below and all lines with
  // avoid updateLowerGridLevels
  // and update testcases accordingly! The leafCell computation doesn't need the lower levels!
  // but testcases fail otherwise!
  // if(!m_levelSetMb) {
  updateLowerGridLevels(mode);
  //}
  if(mode != 0) {
    RECORD_TIMER_STOP(t_a1);
  }

  if(mode == 0) {
    RECORD_TIMER_START(t_b1);
  } else {
    RECORD_TIMER_START(t_a2);
  }
  determineG0Cells(mode);
  if(mode == 0) {
    RECORD_TIMER_STOP(t_b1);
  } else {
    RECORD_TIMER_STOP(t_a2);
  }

  if(mode == 0) {
    RECORD_TIMER_START(t_b2);
  } else {
    RECORD_TIMER_START(t_a3);
  }
  determineBandCells(mode);
  if(mode == 0) {
    RECORD_TIMER_STOP(t_b2);
  } else {
    RECORD_TIMER_STOP(t_a3);
  }

  MInt regrid = 0;
  // check if the level set status changed (i,e, an interface appeared or disappeared)
  // only check, if not called for gap closing and other operations on the collected level-set function [mode 0]
  if(mode != 0) {
    RECORD_TIMER_START(t_a4);
    regrid = (MInt)regridLevelSet();

    MInt newStatus = mMin(1, a_noG0Cells(0));
    for(MInt set = 1; set < m_noSets; set++) {
      newStatus = mMin(newStatus, a_noG0Cells(set));
    }

    MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "status");
    MPI_Allreduce(MPI_IN_PLACE, &newStatus, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "newStatus");
    MPI_Allreduce(MPI_IN_PLACE, &regrid, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "regrid");

    if(m_reconstructOldG) {
      MPI_Allreduce(MPI_IN_PLACE, &m_rotatingReinitTrigger, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "m_rotatingReinitTrigger");
      if(m_rotatingReinitTrigger) {
        reconstructOldGField();
      }
    }

    if(regrid != 0) m_log << "body-movement caused regridding" << endl;

    if(newStatus != status) {
      regrid = 1;
      m_log << "changed number of G0-cells from " << status << " to " << newStatus << " caused regridding" << endl;
    }

    RECORD_TIMER_STOP(t_a4);
  }

  if(regrid != 0) {
    m_forceAdaptation = true;
    // This should force a global mesh-adaptation though the grid-controller at the next time-Step!

    if(domainId() == 0) {
      cerr << "LS-Solver is forcing a mesh-adaptation at time step " << globalTimeStep << " ...";
    }

    m_log << "LS-Solver is forcing a mesh-adaptation at time step " << globalTimeStep << endl;
    m_reinitConvergence = m_reinitConvergenceReset; // this prevents a double reduction of the reinit criterion which
                                                    // can appearently happen in the unified run loop
    m_reinitConvergence = m_reinitConvergence / F3;

    m_log << "reinitConvergence is temporarily decreased by a factor of 3 " << m_reinitConvergence << endl;

    if(!m_levelSetMb || (!this->m_adaptation && m_levelSetMb)) {
      determineG0Cells();

      // reset m_computeSet, since if cells are deleted, G0Cells etc. have to be updated for all sets!
      for(MInt set = 0; set < m_noSets; set++) {
        m_computeSet_tmp[set] = m_computeSet[set];
        m_computeSet[set] = true;
        m_changedSet[set] = true;
      }

      for(MInt set = 0; set < m_noSets; set++)
        std::vector<MInt>().swap(m_bandCells[set]);

      RECORD_TIMER_START(t_a5);
      createGgridCG();
      RECORD_TIMER_STOP(t_a5);

      RECORD_TIMER_START(t_a6);
      determineG0Cells();
      determineBandCells();

      // use updateLowerGridLevels instead of levelSetRestriction, it is more accurate!
      if(m_combustion) {
        levelSetRestriction();
      } else {
        updateLowerGridLevels(mode);
      }
      RECORD_TIMER_STOP(t_a6);

      RECORD_TIMER_START(t_a7);
      setGCellBndryProperty();
      updateBndryCellList();
      RECORD_TIMER_STOP(t_a7);

      // reset m_computeSet, since if cells are deleted, G0Cells etc. had to be updated for all sets!
      for(MInt set = 0; set < m_noSets; set++) {
        m_computeSet[set] = m_computeSet_tmp[set];
      }

      if(domainId() == 0) cerr << " finished." << endl;
    }
  }

  if(!m_levelSetMb) {
    // TODO labels:LS,TIMERS fix timers for different modes
    /* RECORD_TIMER_START(t_a7); */
    updateBndryCellList();
    /* RECORD_TIMER_STOP(t_a7); */
  }

  if(mode == 0) {
    RECORD_TIMER_START(t_b3);
  } else {
    RECORD_TIMER_START(t_a8);
  }
  resetOutsideCells(mode);
  if(mode == 0) {
    RECORD_TIMER_STOP(t_b3);
  } else {
    RECORD_TIMER_STOP(t_a8);
  }

  if(mode == 0) {
    RECORD_TIMER_STOP(t_levelSetTubeZero);
  } else {
    RECORD_TIMER_STOP(t_levelSetTube);
  }
}


/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::fastBuildLevelSetTubeCG() {
  TRACE();


  NEW_TIMER_GROUP(fastBuildLevelSet, "Level Set - build tube");
  NEW_TIMER(t_levelSet, "Total time - build tube", fastBuildLevelSet);
  NEW_SUB_TIMER(t_updateL, "updateLowerGridLevels", t_levelSet);
  NEW_SUB_TIMER(t_detBand, "determineBandCells", t_levelSet);
  NEW_SUB_TIMER(t_detG0, "determineG0Cells", t_levelSet);


  RECORD_TIMER_START(t_levelSet);

  RECORD_TIMER_START(t_updateL);
  updateLowerGridLevels();
  RECORD_TIMER_STOP(t_updateL);

  RECORD_TIMER_START(t_detBand);
  determineG0Cells();
  RECORD_TIMER_STOP(t_detBand);

  RECORD_TIMER_START(t_detG0);
  determineBandCells();
  RECORD_TIMER_STOP(t_detG0);

  updateBndryCellList(); // Stephan

  RECORD_TIMER_STOP(t_levelSet);
}

/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::restartLocalizedLevelSetCG() {
  TRACE();
  stringstream fileName;
  MInt noCells;
  //---

  m_log << " Restarting localized level-set... " << endl;
  if(domainId() == 0) {
    cerr << " Restarting localized level-set... " << endl;
  }

  // initialise minGapWidth and minGapWidthDt1
  if(m_closeGaps) {
    m_minGapWidth.clear();
    m_minGapWidthDt1.clear();
    for(MInt region = 0; region < m_noGapRegions; region++) {
      m_minGapWidth.push_back(std::numeric_limits<MFloat>::max());
      m_minGapWidthDt1.push_back(std::numeric_limits<MFloat>::max());
    }
  }

  if(m_GFieldInitFromSTL) initializeGControlPoint();

  for(MInt set = 0; set < m_noSets; set++)
    std::vector<MInt>().swap(m_bandCells[set]);

  if(m_semiLagrange) {
    initializeIntegrationScheme_semiLagrange();
  } else {
    initializeIntegrationScheme();
  }

  // TODO labels:LS,toremove run_unified remove this, globalTimeStep should not be set by the solver
  if(!m_levelSetMb) globalTimeStep = m_restartTimeStep;

  stringstream levelSetFileName;
  if(!g_multiSolverGrid) {
    levelSetFileName << restartDir() << "restartLSCG";
  } else {
    levelSetFileName << restartDir() << "restartLSCG_" << m_solverId;
  }
  if(!m_useNonSpecifiedRestartFile) levelSetFileName << "_" << m_restartTimeStep;
  levelSetFileName << ParallelIo::fileExt();
  if(domainId() == 0) {
    cerr << " Loading " << (levelSetFileName.str()).c_str() << endl;
  }

  noCells = loadLevelSetGridFlowVarsParCG((levelSetFileName.str()).c_str());
  ASSERT(a_noCells() == noCells && noCells == grid().tree().size(),
         " size of collector is not noCells: " << a_noCells() << " noCells: " << noCells);

  generateListOfGExchangeCellsCG();

  testCellsCG();

  // set g0 cells on the base grid
  determineG0Cells();

  for(MInt set = 0; set < m_noSets; set++) {
    std::vector<MInt>().swap(m_bandCells[set]);
    // very important trigger: determineG0Cells will reset inBand, inGamma,... for all cells
  }

  determineG0Cells();

  determineBandCells();

  exchangeLevelSet();

  setGCellBndryProperty();

  updateBndryCellList();

  resetOutsideCells();

  exchangeAllLevelSetData();

  if(m_virtualSurgery && isActive()) {
    initializeCollectedLevelSet(2);
  }

  if(m_levelSetRans) {
    updateAllLowerGridLevels();
  }

  // initialize fields...
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(m_maxNoSets > 1 && m_closeGaps) {
      a_potentialGapCell(cellId) = 0;
      a_potentialGapCellClose(cellId) = 0;
    }
  }

  if(m_combustion) {
    determineMinMaxMeanInterfacePosition();
    computeZeroLevelSetArcLength();
  }
}


/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::loadLevelSetGridFlowVarsParCG(const MChar* fileName) {
  TRACE();


  IF_CONSTEXPR(nDim != 2 && nDim != 3) {
    cerr << " In global function loadGridFlowVarsPar: wrong number of dimensions !" << endl;
    mTerm(1, AT_);
    return 0;
  }

  using namespace maia::parallel_io;
  ParallelIo data(fileName, PIO_READ, mpiComm());

  if(m_LSSolver) {
    // set level set time
    data.getAttribute(&m_time, "levelSetTime");
  }

  MInt noCells = 0;
  MInt noGridCell = grid().tree().size();
  MInt noInternalGCells = 0;

  for(MInt i = 0; i < noGridCell; i++) {
    noCells++;
    if(grid().tree().hasProperty(i, Cell::IsHalo)) {
      continue;
    }
    ++noInternalGCells;
  }

  MFloatScratchSpace tmpFloat(noCells, AT_, "tmpFloat");
  MLongScratchSpace tmpInt(noCells, AT_, "tmpInt");
  tmpFloat.fill(-2.0);
  tmpInt.fill(-2);

  vector<MString> variableNames = data.getDatasetNames(1);
  MString dataVarName;
  MInt dataVarId;

  // This should be the same for all the variables
  ParallelIo::size_type dimLen = noInternalGCells;
  ParallelIo::size_type start = grid().domainOffset(domainId()) - grid().raw().m_32BitOffset;
  // set offset for all read operations
  data.setOffset(dimLen, start);

  data.getAttribute(&dataVarName, "name", variableNames[0]);
  ASSERT(dataVarName == "G" || dataVarName == "G_0", "ERROR retsart-file, wrong Variable-name ordering");
  data.readArray(&tmpFloat[0], variableNames[0]);

  dataVarName.clear();
  dataVarId = -1;
  for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
    data.getAttribute(&dataVarName, "name", variableNames[i]);
    if("regrid" == dataVarName) {
      dataVarId = i;
      break;
    }
  }

  data.readArray(&tmpInt[0], variableNames[dataVarId]);

  // exchange G0-values to get correct values on halo cells!
  // otherwise halo cells will not be detected as g cells
  exchangeDataLS(&tmpFloat[0], 1);
  exchangeDataLS(&tmpInt[0], 1);

  for(MInt i = 0; i < noGridCell; ++i) {
    a_appendCollector();

    a_regridTriggerG(i) = tmpInt[i];
    a_levelSetFunctionG(i, 0) = tmpFloat[i];
  }

  if(m_maxNoSets > 1) {
    for(MInt set = 1; set < m_noSets; set++) {
      dataVarName.clear();
      dataVarId = -1;
      string tmps = "G_" + to_string(set);
      for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
        data.getAttribute(&dataVarName, "name", variableNames[i]);
        if(tmps == dataVarName) {
          dataVarId = i;
          break;
        }
      }
      if(dataVarId == -1)
        mTerm(1, AT_, "ERROR: More levelsets specified in the properties file than existing in the restart-file!");

      data.readArray(&tmpFloat[0], variableNames[dataVarId]);
      for(MInt i = 0; i < noGridCell; ++i) {
        a_levelSetFunctionG(i, set) = tmpFloat[i];
      }
    }
  }

  if(m_semiLagrange && m_maxNoSets == 1) {
    dataVarName.clear();
    dataVarId = -1;
    for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
      data.getAttribute(&dataVarName, "name", variableNames[i]);
      if("oldG" == dataVarName) {
        dataVarId = i;
        break;
      }
    }
    data.readArray(&tmpFloat[0], variableNames[dataVarId]);
    for(MInt i = 0; i < noGridCell; ++i) {
      a_oldLevelSetFunctionG(i, 0) = tmpFloat[i];
    }
  } else if(m_semiLagrange && m_maxNoSets > 1) {
    for(MInt set = 0; set < m_noSets; set++) {
      string tmps = "oldG_" + to_string(set);
      dataVarName.clear();
      dataVarId = -1;
      for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
        data.getAttribute(&dataVarName, "name", variableNames[i]);
        if(tmps == dataVarName) {
          dataVarId = i;
          break;
        }
      }
      data.readArray(&tmpFloat[0], variableNames[dataVarId]);
      for(MInt i = 0; i < noGridCell; ++i) {
        a_oldLevelSetFunctionG(i, set) = tmpFloat[i];
      }
    }
  }


  if(m_semiLagrange) {
    if(m_LsRotate) {
      if(!m_reconstructOldG) {
        MInt body;
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          body = m_bodiesToCompute[b];

          string tmps = "containingCell_" + to_string(body);
          dataVarName.clear();
          dataVarId = -1;
          for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
            data.getAttribute(&dataVarName, "name", variableNames[i]);
            if(tmps == dataVarName) {
              dataVarId = i;
              break;
            }
          }
          data.readArray(&tmpInt[0], variableNames[dataVarId]);

          for(MInt i = 0; i < noGridCell; ++i) {
            a_containingCell(i, b) = tmpInt[i];
          }
          globalToLocalIdsContainingCells();
        }

        dataVarName.clear();
        dataVarId = -1;
        for(MInt i = 0; i < (signed)variableNames.size(); ++i) {
          data.getAttribute(&dataVarName, "name", variableNames[i]);
          if("initialGCell" == dataVarName) {
            dataVarId = i;
            break;
          }
        }
        data.readArray(&tmpInt[0], variableNames[dataVarId]);
        for(MInt i = 0; i < noGridCell; ++i) {
          m_initialGCell[i] = tmpInt[i];
        }
      }

      for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
        for(MInt j = 0; j < nDim; j++) {
          MInt id = i * nDim + j;
          string tmps = "SL_xRot_" + to_string(id);
          data.readScalar(&m_semiLagrange_xRot_ref[id], tmps.c_str());
        }
      }

      if(m_reconstructOldG) {
        for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
          for(MInt j = 0; j < nDim; j++) {
            MInt id = i * nDim + j;
            string tmps = "SL_xRot_STL" + to_string(id);
            data.readScalar(&m_semiLagrange_xRot_STL[id], tmps.c_str());
          }
        }
      } else {
        std::copy_n(&m_semiLagrange_xRot_ref[0], nDim * m_noEmbeddedBodies, &m_semiLagrange_xRot_STL[0]);
      }
    }

    for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
      for(MInt j = 0; j < nDim; j++) {
        MInt id = i * nDim + j;
        string tmps = "SL_xShift_" + to_string(id);
        data.readScalar(&m_semiLagrange_xShift_ref[id], tmps.c_str());
      }
    }
  }

  if(m_noGapRegions > 0) {
    for(MInt i = 0; i < m_noGapRegions; i++) {
      string tmps = "deltaMin_" + to_string(i);
      data.readScalar(&m_minGapWidth[i], tmps.c_str());
    }
  }


  // generate the window and halo cells of the level set
  generateListOfGExchangeCellsCG();
  // exchange to get correct values for halo cells

  exchangeAllLevelSetData();

  return a_noCells();
}

/**
 * \brief write solver restart file and matching grid file!
 * NOTE: for debugging purposes only!
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::writeRestartLevelSetFileCG(MBool writeRestart,
                                                         const MString& levelSetFileNameGGrid,
                                                         const MString& levelSetFileNameGSol) {
  if(m_restartInterval == -1) return;
  if(((globalTimeStep % m_restartInterval) == 0 && globalTimeStep > m_restartTimeStep) || writeRestart) {
    if(domainId() == 0) {
      cerr << endl;
      cerr << "Writing levelset CG restart file at time step " << globalTimeStep << " ..." << endl;
    }
    writeRestart = true;
  }
  if(!writeRestart) return;

  std::ignore = levelSetFileNameGGrid;
  stringstream gridFile;
  stringstream gridFilePath;
  gridFile << "grid_" << globalTimeStep << ".Netcdf";
  gridFilePath << outputDir() << "grid_" << globalTimeStep << ".Netcdf";
  MIntScratchSpace recalcIds(grid().raw().treeb().size(), AT_, "recalcIds");

  grid().raw().saveGrid((gridFilePath.str()).c_str(), recalcIds.begin());

  ASSERT(m_currentFileName.empty(), "");
  m_currentFileName = levelSetFileNameGSol;

  writeRestartFile(true, false, (gridFile.str()).c_str(), recalcIds.begin());

  m_currentFileName.clear();
}

/**
 * \brief
 * \author Sohel Herff
 */
template <MInt nDim>
MBool LsCartesianSolver<nDim>::levelSetAdaptationTrigger() {
  TRACE();
  MInt mode = -1;

  MIntScratchSpace globalRegrid(1, AT_, "globalRegrid");
  MIntScratchSpace regrid(1, AT_, "regrid");
  MIntScratchSpace globalStatus(1, AT_, "globalStatus");
  MIntScratchSpace status(1, AT_, "status");
  MIntScratchSpace newStatus(1, AT_, "newStatus");

  // if in the complete global domain, there are no G0 cells present, currently no LVS is computed...
  status.p[0] = mMin(1, a_noG0Cells(0));
  for(MInt set = 1; set < m_noSets; set++) {
    status.p[0] = mMin(status.p[0], a_noG0Cells(set));
  }
  MPI_Allreduce(status.getPointer(), globalStatus.getPointer(), 1, MPI_INT, MPI_MAX, mpiComm(), AT_,
                "status.getPointer()", "globalStatus.getPointer()");
  status.p[0] = globalStatus.p[0];

  // required here, since determineG0Cells may rely on lower grid levels...
  updateLowerGridLevels(mode);
  determineG0Cells(mode);
  determineBandCells(mode);
  updateBndryCellList();
  regrid.p[0] = 0;
  // check if the level set status changed (i,e, an interface appeared or disappeared)
  // only check, if not called for gap closing and other operations on the collected level-set function [mode 0]
  regrid.p[0] = (MInt)regridLevelSet();
  newStatus.p[0] = mMin(1, a_noG0Cells(0));
  for(MInt set = 1; set < m_noSets; set++) {
    newStatus.p[0] = mMin(newStatus.p[0], a_noG0Cells(set));
  }

  MPI_Allreduce(newStatus.getPointer(), globalStatus.getPointer(), 1, MPI_INT, MPI_MAX, mpiComm(), AT_,
                "newStatus.getPointer()", "globalStatus.getPointer()");
  newStatus.p[0] = globalStatus.p[0];

  MPI_Allreduce(regrid.getPointer(), globalRegrid.getPointer(), 1, MPI_INT, MPI_MAX, mpiComm(), AT_,
                "regrid.getPointer()", "globalRegrid.getPointer()");
  regrid.p[0] = globalRegrid.p[0];

  if(newStatus.p[0] != status.p[0]) regrid.p[0] = true;


  return regrid.p[0];
}

/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::refineCell(const MInt gridCellId) {
  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
    const MInt childId = grid().raw().treeb().child(gridCellId, child);

    if(childId == -1) continue;

    // @ansgar_pls_adapt2 test this for ls solver!
    //    --> this sort of skip at least causes one assert in compactCells to fail (see testcase
    //        LS/2D_twoLsBlock_slottedDisk_different_rotations_single with partitionCellMaxNoOffspring=1
    // Skip if cell is a partition level ancestor and its child was not newly created
    if(!grid().raw().a_hasProperty(childId, Cell::WasNewlyCreated)
       && grid().raw().a_hasProperty(gridCellId, Cell::IsPartLvlAncestor)) {
      continue;
    }

    if(!g_multiSolverGrid) ASSERT(grid().raw().a_hasProperty(childId, Cell::WasNewlyCreated), "");

    // If solver is inactive all cells musst be halo cells!
    if(!isActive()) ASSERT(grid().raw().a_isHalo(childId), "");

    if(grid().azimuthalPeriodicity()) {
      // Solver flag will be set in proxy
      // The cartesianGrid does not know the geometry, therefore it needs to be checked
      // that the child actually is inside the geometry
      if(grid().checkOutsideGeometry(childId) == 1) {
        grid().raw().setSolver(childId, solverId(), false);
        continue;
      }
    }

    // If child exists in grid but is not located inside solver geometry
    if(!grid().solverFlag(childId, solverId())) continue;

    const MInt solverChildId = this->createCellId(childId);

    if(!g_multiSolverGrid) ASSERT(solverChildId == childId, "");

    // avoid faulty initialisation of childs
    for(MInt set = m_startSet; set < m_noSets; set++) {
      a_inBandG(solverChildId, set) = a_inBandG(solverCellId, set);
    }

    for(MInt set = 0; set < m_noSets; set++) {
      a_levelSetFunctionG(solverChildId, set) = a_levelSetFunctionG(solverCellId, set);
      if(m_semiLagrange) {
        a_oldLevelSetFunctionG(solverChildId, set) = a_oldLevelSetFunctionG(solverCellId, set);
      }
      if(m_levelSetMb) a_bodyIdG(solverChildId, set) = a_bodyIdG(solverCellId, set);
    }

    if(!a_isHalo(solverCellId) && globalTimeStep > 0) {
      if(m_virtualSurgery) {
        if((a_levelSetFunctionG(solverCellId, 1) * a_levelSetFunctionG(solverCellId, 2)) < 0) {
          m_refinedCells.insert(make_pair(solverCellId, 0));
        }
      } else {
        MInt cellAdded = 0;
        for(MInt set = m_startSet; set < m_noSets; set++) {
          if(a_inBandG(solverCellId, set)) {
            if(!m_computeSet_backup[set] || m_maxLevelChange) {
              cellAdded++;
              if(cellAdded > 1) {
                auto it0 = m_refinedCells.find(solverChildId);
                ASSERT(it0 != m_refinedCells.end(), "");
                m_refinedCells.erase(it0);
                m_refinedCells.insert(make_pair(solverChildId, 0));
              } else {
                m_refinedCells.insert(make_pair(solverChildId, set));
              }
            }

            if(m_computeSet_backup[set]) {
              // tested only for m_maxLevelChange, in this case the levelSet is interpolated
              // and the old-levelSet is reconstructed to maintain mass-conservation!
              MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
              MFloat cellPos[3] = {c_coordinate(solverChildId, 0), c_coordinate(solverChildId, 1),
                                   c_coordinate(solverChildId, nDim - 1)};
              const MInt position = setUpLevelSetInterpolationStencil(solverCellId, interpolationCells, cellPos);
              if(position > -1) {
                if(!m_maxLevelChange) {
                  a_oldLevelSetFunctionG(solverChildId, set) = interpolateOldLevelSet(interpolationCells, cellPos, set);
                }
                a_levelSetFunctionG(solverChildId, set) = interpolateLevelSet(interpolationCells, cellPos, set);
              }
            }
          }
        }
      }
    }

    if(m_LsRotate) {
      for(MInt b = 0; b < m_noBodiesToCompute; b++) {
        a_containingCell(solverChildId, b) = -1;
        if(!m_reconstructOldG) {
          a_containingDomain(solverChildId, b) = -1;
          m_initialGCell[solverChildId] = 0;
        }
      }
    }
  }
}

/**
 * \brief
 * \author Sohel Herff (m_levelSetMb part by: Tim Wegmann)
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::removeChilds(const MInt gridCellId) {
  // If solver is inactive cell musst never be a internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  for(MInt set = 0; set < m_noSets; set++) {
    a_inBandG(solverCellId, set) = false;
  }


  ASSERT(solverCellId > -1 && solverCellId < m_cells.size(), "solverCellId is: " << solverCellId);
  if(!g_multiSolverGrid) ASSERT(solverCellId == gridCellId, "");

  for(MInt set = 0; set < m_noSets; set++) {
    a_levelSetFunctionG(solverCellId, set) = F0;
    if(m_semiLagrange) {
      a_oldLevelSetFunctionG(solverCellId, set) = F0;
    }
  }
  MInt noChildren = 0;

  for(MInt c = 0; c < grid().m_maxNoChilds; c++) {
    MInt childId = c_childId(solverCellId, c);
    if(childId < 0) continue;
    noChildren++;
    for(MInt set = 0; set < m_noSets; set++) {
      a_levelSetFunctionG(solverCellId, set) += a_levelSetFunctionG(childId, set);
      if(m_semiLagrange) {
        a_oldLevelSetFunctionG(solverCellId, set) += a_oldLevelSetFunctionG(childId, set);
      }
    }

    if(m_LsRotate) {
      for(MInt b = 0; b < m_noBodiesToCompute; b++) {
        a_containingCell(childId, b) = -1;
        if(!m_reconstructOldG) {
          a_containingDomain(childId, b) = -1;
          if(!a_isHalo(childId) && m_initialGCell[childId] == 1) {
            MInt parentId = c_parentId(childId);
            cerr << c_globalId(childId) << "," << c_globalId(parentId) << "," << a_level(childId) << ","
                 << a_isHalo(childId) << "," << a_isHalo(parentId) << endl;
            mTerm(1, AT_, "Initial G cell is deleted!");
          }
        }
      }
    }
    for(MInt set = m_startSet; set < m_noSets; set++) {
      a_inBandG(solverCellId, set) = a_inBandG(solverCellId, set) || a_inBandG(childId, set);
    }

    this->removeCellId(childId);
  }
  for(MInt set = 0; set < m_noSets; set++) {
    a_levelSetFunctionG(solverCellId, set) /= noChildren;
    if(m_semiLagrange) {
      a_oldLevelSetFunctionG(solverCellId, set) /= noChildren;
    }
  }

  if(m_maxLevelChange && globalTimeStep > 0 && !a_isHalo(solverCellId)) {
    MInt cellAdded = 0;
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(a_inBandG(solverCellId, set)) {
        cellAdded++;
        if(cellAdded > 1) {
          auto it0 = m_refinedCells.find(solverCellId);
          ASSERT(it0 != m_refinedCells.end(), "");
          m_refinedCells.erase(it0);
          m_refinedCells.insert(make_pair(solverCellId, 0));
        } else {
          m_refinedCells.insert(make_pair(solverCellId, set));
        }
      }
    }
  }

  if(!g_multiSolverGrid) {
    ASSERT((grid().raw().treeb().size() - m_cells.size()) <= grid().m_maxNoChilds, "");
  }
}

/**
 * \brief
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::removeCell(const MInt gridCellId) {
  // If solver is inactive cell musst never be a internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(gridCellId > -1 && gridCellId < grid().raw().treeb().size() && solverCellId > -1
             && solverCellId < m_cells.size() && grid().tree().solver2grid(solverCellId) == gridCellId,
         "");

  this->removeCellId(solverCellId);
}

// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
template <MInt nDim>
void LsCartesianSolver<nDim>::resizeGridMap() {
  grid().resizeGridMap(m_cells.size());
}

/**
 * \brief
 * \author Sohel Herff, Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::swapCells(const MInt cellId0, const MInt cellId1) {
  const MInt size = m_cells.size();
  m_cells.append();
  m_cells.erase(size);
  m_cells.copy(cellId0, size);
  m_cells.copy(cellId1, cellId0);
  m_cells.copy(size, cellId1);
  m_cells.erase(size);
  m_cells.size(size);

  if(m_LsRotate) {
    // cellId1 is moved ahead
    if(m_reconstructOldG) {
      // Property is already been swapped
      if(m_swapIds.find(cellId1) == m_swapIds.end()) m_swapIds[cellId1] = cellId0;
    } else {
      if(m_initialGCell[cellId1] == 1) m_swapIds[cellId1] = cellId0;
      std::swap(m_initialGCell[cellId0], m_initialGCell[cellId1]);
      m_initialGCell[cellId1] = 0;
    }
  }

  if(!m_refinedCells.empty()) {
    auto it0 = m_refinedCells.find(cellId0);
    auto it1 = m_refinedCells.find(cellId1);
    if(it0 != m_refinedCells.end() && it1 != m_refinedCells.end()) {
      std::swap(it0->second, it1->second);
    } else if(it0 != m_refinedCells.end()) {
      MInt set = it0->second;
      m_refinedCells.erase(it0);
      m_refinedCells.insert(make_pair(cellId1, set));
    } else if(it1 != m_refinedCells.end()) {
      MInt set = it1->second;
      m_refinedCells.erase(it1);
      m_refinedCells.insert(make_pair(cellId0, set));
    }
  }


  if(!m_oldG0Cells.empty()) {
    auto it0 = m_oldG0Cells.find(cellId0);
    auto it1 = m_oldG0Cells.find(cellId1);
    if(it0 != m_oldG0Cells.end() && it1 != m_oldG0Cells.end()) {
      std::swap(it0->second, it1->second);
    } else if(it0 != m_oldG0Cells.end()) {
      const MInt set = it0->second;
      m_oldG0Cells.erase(it0);
      m_oldG0Cells.insert(make_pair(cellId1, set));
    } else if(it1 != m_oldG0Cells.end()) {
      const MInt set = it1->second;
      m_oldG0Cells.erase(it1);
      m_oldG0Cells.insert(make_pair(cellId0, set));
    }
  }
}


/**
 * \brief
 * \author Lennart Schneiders
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::swapProxy(const MInt cellId0, const MInt cellId1) {
  grid().swapGridIds(cellId0, cellId1);
}


/**
 * \brief sets the cell-weight for balancing and a restarting
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt offset = noCellsGrid * solverId();

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    const MInt id = gridCellId + offset;
    solverCellWeight[id] = 0.0;
    if(a_isHalo(cellId)) continue;
    solverCellWeight[id] = m_weightBaseCell * m_weightMulitSolverFactor;
    if(c_noChildren(cellId) > 0) continue;
    solverCellWeight[id] = m_weightLeafCell * m_weightMulitSolverFactor;
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet[set]) continue;
      if(a_inBandG(cellId, set)) {
        solverCellWeight[id] = m_weightBandCell * m_weightMulitSolverFactor;
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------

/**
 * \brief This function prepares a restart that is handled by the grid-controller!
 * \author Tim Wegmann
 */

template <MInt nDim>
MBool LsCartesianSolver<nDim>::prepareRestart(MBool writeRestart, MBool& writeGridRestart) {
  TRACE();

  writeGridRestart = false;

  if(((globalTimeStep % m_restartInterval) == 0) || writeRestart) {
    writeRestart = true;

    if(m_adaptationSinceLastRestart) {
      writeGridRestart = true;
    }

    // update updateLowerGridLevels before a restart!
    // avoid updateLowerGridLevels
    // if(m_levelSetMb && isActive()) {
    //  updateLowerGridLevels();
    //}
  }

  return writeRestart;
}

/**
 * \brief This function resets the grid-trigger after a restart that is handled by the grid-controller!
 * \author Tim Wegmann
 */

template <MInt nDim>
void LsCartesianSolver<nDim>::reIntAfterRestart(MBool doneRestart) {
  TRACE();

  if(doneRestart) {
    m_adaptationSinceLastRestart = false;
  }
}

/**
 * \brief This function writes restart that is handled by the grid-controller!
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::writeRestartFile(const MBool writeRestart, const MBool writeBackup,
                                               const MString gridFileName, MInt* recalcIdTree) {
  TRACE();

  m_currentGridFileName = gridFileName;

  if(writeRestart) {
    saveRestartFile(writeBackup, &recalcIdTree[0]);
  }
}

/**
 * \brief This function writes restart that is handled by the grid-controller!
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::saveRestartFile(const MBool writeBackup, MInt* recalcIds) {
  TRACE();

  if(domainId() == 0) {
    cerr << "Writing levelset restart file for solver " << m_solverId << " at time step " << globalTimeStep << " ...";
  }

  MBool debugOutput = false;
#ifdef LS_DEBUG
  debugOutput = true;
#endif


#if defined LS_DEBUG
  for(MInt body = 0; body < m_noEmbeddedBodies; body++) {
    const MInt bcId = m_bodyBndryCndIds[body];
    stringstream controlfile;
    controlfile << "Ctrl_Body_" << body << "_" << globalTimeStep << ".stl";
    m_gCtrlPnt.CtrlPnt2_CtrlPntToSTL((controlfile.str()).c_str(), bcId);
  }
#endif

  MInt noCells;
  MInt noInternalCellIds;
  std::vector<MInt> recalcIdsSolver(0);
  std::vector<MInt> reOrderedCells(0);
  this->calcRecalcCellIdsSolver(recalcIds, noCells, noInternalCellIds, recalcIdsSolver, reOrderedCells);

  MInt noIdVars = 1;
  if(m_LsRotate && !m_reconstructOldG) noIdVars += 1;
  MInt noDbVars = m_maxNoSets;
  if(m_semiLagrange) noDbVars += m_maxNoSets;
  if(m_semiLagrange && m_bodyIdOutput) noDbVars += m_maxNoSets;
  if(m_LsRotate && !m_reconstructOldG) noDbVars += m_noBodiesToCompute;
  if(debugOutput) {
    noDbVars = noDbVars + 4;
  }
  const MInt noIdParams = 0;
  MInt noDbParams = m_semiLagrange ? nDim * m_noEmbeddedBodies : 0;
  noDbParams = noDbParams + (MInt)m_noGapRegions;
  if(m_LsRotate) noDbParams += (2 * nDim * m_noBodiesToCompute);
  MIntScratchSpace idVariables(noCells * noIdVars, AT_, "idVariables");
  MFloatScratchSpace dbVariables(noCells * noDbVars, AT_, "dbVariables");
  MIntScratchSpace idParameters(noIdParams, AT_, "idParameters");
  MFloatScratchSpace dbParameters(noDbParams, AT_, "dbParameters");
  vector<MString> dbVariablesName;
  vector<MString> idVariablesName;
  vector<MString> dbParametersName;
  vector<MString> idParametersName;
  vector<MString> name;

  MFloatScratchSpace levelSet(noCells, AT_, "levelSet");
  MIntScratchSpace regridL(noCells, AT_, "regridL");
  MInt tmpSize = debugOutput ? 1 : 0;
  MFloatScratchSpace tmpW(noCells * tmpSize, AT_, "tmpw");
  regridL.fill(3);

  if(m_levelSetRans) {
    updateLowerGridLevels();
  }

  if(m_maxNoSets == 1) {
    name.clear();
    name.push_back("G");
    if(grid().newMinLevel() < 0) {
      for(MInt cell = 0; cell < noCells; cell++) {
        levelSet[cell] = a_levelSetFunctionG(cell, 0);
      }
    } else {
      for(MInt cell = 0; cell < noCells; cell++) {
        levelSet[cell] = a_levelSetFunctionG(reOrderedCells[cell], 0);
      }
    }
    this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  } else {
    for(MInt set = 0; set < m_noSets; set++) {
      string tmps = "G_" + to_string(set);
      name.clear();
      name.push_back(tmps);
      if(grid().newMinLevel() < 0) {
        for(MInt cell = 0; cell < noCells; cell++) {
          levelSet[cell] = a_levelSetFunctionG(cell, set);
        }
      } else {
        for(MInt cell = 0; cell < noCells; cell++) {
          levelSet[cell] = a_levelSetFunctionG(reOrderedCells[cell], set);
        }
      }
      this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
    }
  }

  if(m_semiLagrange) {
    if(m_maxNoSets == 1) {
      name.clear();
      name.push_back("oldG");
      if(grid().newMinLevel() < 0) {
        for(MInt cell = 0; cell < noCells; cell++) {
          levelSet[cell] = a_oldLevelSetFunctionG(cell, 0);
        }
      } else {
        for(MInt cell = 0; cell < noCells; cell++) {
          levelSet[cell] = a_oldLevelSetFunctionG(reOrderedCells[cell], 0);
        }
      }
      this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
    } else {
      for(MInt set = 0; set < m_noSets; set++) {
        string tmps = "oldG_" + to_string(set);
        name.clear();
        name.push_back(tmps);
        if(grid().newMinLevel() < 0) {
          for(MInt cell = 0; cell < noCells; cell++) {
            levelSet[cell] = a_oldLevelSetFunctionG(cell, set);
          }
        } else {
          for(MInt cell = 0; cell < noCells; cell++) {
            levelSet[cell] = a_oldLevelSetFunctionG(reOrderedCells[cell], set);
          }
        }
        this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
      }
    }
    if(m_bodyIdOutput) {
      for(MInt i = 0; i < m_noSets; i++) {
        string tmps = "bodyId_" + to_string(i);
        name.clear();
        name.push_back(tmps);
        if(grid().newMinLevel() < 0) {
          for(MInt cell = 0; cell < noCells; cell++) {
            levelSet[cell] = a_bodyIdG(cell, i);
          }
        } else {
          for(MInt cell = 0; cell < noCells; cell++) {
            levelSet[cell] = a_bodyIdG(reOrderedCells[cell], i);
          }
        }
        this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
      }
    }
    if(m_LsRotate) {
      if(!m_reconstructOldG) {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          MInt body = m_bodiesToCompute[b];

          string tmps = "containingCell_" + to_string(body);
          name.clear();
          name.push_back(tmps);

          localToGlobalIdsContainingCells();
          for(MInt cell = 0; cell < noCells; cell++) {
            levelSet[cell] = a_containingCell(cell, b);
          }
          globalToLocalIdsContainingCells();

          this->collectVariables(levelSet.begin(), dbVariables, name, dbVariablesName, 1, noCells);
        }

        name.clear();
        name.push_back("initialGCell");
        for(MInt cell = 0; cell < noCells; cell++) {
          regridL[cell] = m_initialGCell[cell];
        }
        this->collectVariables(regridL.begin(), idVariables, name, idVariablesName, 1, noCells);
      }

      for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
        for(MInt j = 0; j < nDim; j++) {
          MInt id = i * nDim + j;
          string tmps = "SL_xRot_" + to_string(id);
          this->collectParameters(m_semiLagrange_xRot_ref[id], dbParameters, tmps.c_str(), dbParametersName);
        }
      }

      for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
        for(MInt j = 0; j < nDim; j++) {
          MInt id = i * nDim + j;
          string tmps = "SL_xRot_STL" + to_string(id);
          this->collectParameters(m_semiLagrange_xRot_STL[id], dbParameters, tmps.c_str(), dbParametersName);
        }
      }
    }


    for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
      for(MInt j = 0; j < nDim; j++) {
        const MInt id = i * nDim + j;
        const string tmps = "SL_xShift_" + to_string(id);
        this->collectParameters(m_semiLagrange_xShift_ref[id], dbParameters, tmps.c_str(), dbParametersName);
      }
    }
  }


  if(m_noGapRegions > 0) {
    for(MInt i = 0; i < m_noGapRegions; i++) {
      string tmps = "deltaMin_" + to_string(i);
      this->collectParameters(m_minGapWidthDt1[i], dbParameters, tmps.c_str(), dbParametersName);
    }
  }
  name.clear();
  name.push_back("regrid");
  if(grid().newMinLevel() < 0) {
    for(MInt cell = 0; cell < noCells; cell++) {
      regridL[cell] = a_regridTriggerG(cell);
    }
  } else {
    for(MInt cell = 0; cell < noCells; cell++) {
      regridL[cell] = a_regridTriggerG(reOrderedCells[cell]);
    }
  }
  this->collectVariables(regridL.begin(), idVariables, name, idVariablesName, 1, noCells);

  if(debugOutput) {
    if(grid().newMinLevel() > 0) mTerm(1, AT_, "Not implemented yet!");
    name.clear();
    name.push_back("Window");
    for(MInt i = 0; i < noCells; i++) {
      if(a_isWindow(i)) {
        tmpW[i] = -1;
      } else {
        tmpW[i] = domainId();
      }
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

    name.clear();
    name.push_back("globalId");
    for(MInt i = 0; i < noCells; i++) {
      tmpW[i] = c_globalId(i);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
    name.clear();
    name.push_back("BandCells");
    for(MInt i = 0; i < noCells; i++) {
      tmpW[i] = -2;
      for(MInt set = 0; set < m_noSets; set++) {
        if(a_isGBoundaryCellG(i, set)) tmpW[i] = -1;
        if(a_inBandG(i, set)) tmpW[i] = set;
      }
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
    name.clear();
    name.emplace_back("G0-Cells");
    for(MInt i = 0; i < noCells; i++) {
      tmpW[i] = -1;
      for(MInt set = 0; set < m_noSets; set++) {
        for(MInt id = 0; id < a_noG0Cells(set); id++) {
          MInt cellId = a_G0CellId(id, set);
          tmpW[cellId] = set;
        }
      }
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

    /*  name.clear();
        name.push_back("gapWidth");
        for (MInt i = 0; i < noCells; i++){
        tmpW[i] = a_gapWidth(i);
        }
        this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1,noCells);
        name.clear();
        name.push_back("secondBodyId");
        for (MInt i = 0; i < noCells; i++){
        tmpW[i] = a_secondBodyId( i );
        }
        this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1,noCells);
    */
  }

  stringstream levelSetFileName;
  levelSetFileName.clear();
  levelSetFileName.str("");
  if(m_currentFileName.empty()) {
    if(!g_multiSolverGrid) {
      levelSetFileName << outputDir() << "restartLSCG";
    } else {
      levelSetFileName << outputDir() << "restartLSCG_" << m_solverId;
    }

    if(!m_useNonSpecifiedRestartFile) {
      levelSetFileName << "_" << globalTimeStep;
    }

    levelSetFileName << ParallelIo::fileExt();

  } else {
    levelSetFileName << outputDir() << m_currentFileName;
    if(!m_useNonSpecifiedRestartFile) levelSetFileName << "_" << globalTimeStep;
    levelSetFileName << ParallelIo::fileExt();
  }

  MFloat time = -1;
  if(m_LSSolver) time = m_time;

  MInt* pointerRecalcIds = (recalcIds == nullptr) ? nullptr : recalcIdsSolver.data();
  if(writeBackup) {
    stringstream levelSetBackupFileName;
    levelSetBackupFileName.clear();
    levelSetBackupFileName.str("");
    levelSetBackupFileName << outputDir() << "restartLSCGBackup_" << globalTimeStep;
    levelSetBackupFileName << ParallelIo::fileExt();
    if(domainId() == 0) cerr << "Writing level set (backup) for the ls-solver... ";

    this->saveGridFlowVars((levelSetBackupFileName.str()).c_str(), m_currentGridFileName.c_str(), noCells,
                           noInternalCellIds, dbVariables, dbVariablesName, 0, idVariables, idVariablesName, 0,
                           dbParameters, dbParametersName, idParameters, idParametersName, pointerRecalcIds, time);
  }

  this->saveGridFlowVars((levelSetFileName.str()).c_str(), m_currentGridFileName.c_str(), noCells, noInternalCellIds,
                         dbVariables, dbVariablesName, 0, idVariables, idVariablesName, 0, dbParameters,
                         dbParametersName, idParameters, idParametersName, pointerRecalcIds, time);

  if(domainId() == 0) cerr << "ok" << endl;
}

// --------------------------------------------------------------------------------------


/*
 * \brief finalize levelSet solver for rotating levelSet
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::finalizeInitSolver() {
  TRACE();
  if(m_combustion && globalTimeStep < 1) {
    initSolver();
  }

  // Nothing to be done if solver is not active
  if(!isActive()) return;

  if(m_levelSetMb) {
    if(m_constructGField) return;

    checkHaloCells();

    for(MInt set = 0; set < m_noSets; set++) {
      m_computeSet[set] = m_computeSet_tmp[set];
      m_changedSet[set] = true;
    }

    for(MInt set = 0; set < m_noSets; set++) {
      ASSERT(m_computeSet[set] == m_computeSet_backup[set], "ERROR in m_computeSet");
    }

    if(m_trackMovingBndry != 0 && globalTimeStep >= m_trackMbStart && globalTimeStep < m_trackMbEnd) {
      // set a_wasGZeroCell to a_isGZeroCell before the timeStep!
      for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
        for(MInt set = 0; set < m_noSets; set++) {
          a_wasGZeroCell(cellId, set) = a_isGZeroCell(cellId, set);
        }
      }

      if(!m_semiLagrange) {
        setGCellBndryProperty();

        computeNormalVectors();

        computeCurvature();

        determinePropagationSpeed();
      }


      testCellsCG();

      // if buildCollectedLevelSet, this should only be working on the individual-sets!
      buildLevelSetTube();
      setBandNewArrivals();

      if(!m_semiLagrange) {
        computeNormalVectors();
        computeCurvature();
        levelSetReinitialization();
      }

      // wokring on the collected levelSet!
      buildMultipleLevelSet();
    }

    if(m_LsRotate) {
      MInt ind;
      for(MInt i = 0; i < m_noEmbeddedBodies; i++) {
        for(MInt j = 0; j < nDim; j++) {
          ind = i * nDim + j;
          m_bodyAngularAcceleration[ind] = F0;
        }
        rotateLevelSet(5, &m_bodyAngularVelocity[i * nDim], i, nullptr, nullptr, &m_semiLagrange_xRot_STL[i * nDim]);
      }
      updateContainingGCells(1);
      copyWindowToHaloIds();
    }
  }
}

/*
 * \brief: Initialize lists for rotating levelset
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::initRotatingLS() {
  TRACE();

  if(!m_LsRotate) return;

  MFloat maRot;
  MInt ind;
  for(MInt b = 0; b < m_noEmbeddedBodies; b++) {
    for(MInt i = 0; i < nDim; i++) {
      ind = b * nDim + i;
      maRot = Context::getSolverProperty<MFloat>("MaRot", solverId(), AT_, ind);
      m_omega[b * nDim + i] = m_referenceLength / m_bodyRadius[b] * m_referenceVelocity * maRot;
    }
  }
  if(m_restart) {
    if(m_reconstructOldG) {
      m_rotatingReinitTrigger = 1;
      resetContainingGCells();
    }
    return;
  }
  if(m_initialRefinement) return;

  resetContainingGCells();
}

/*
 * \brief: resets list for rotating levelset
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::resetContainingGCells() {
  TRACE();

  for(MInt cellId = 0; cellId < m_maxNoCells; cellId++) {
    if(a_isHalo(cellId)) continue;

    for(MInt b = 0; b < m_noBodiesToCompute; b++) {
      a_containingCell(cellId, b) = cellId;
    }
  }

  if(m_reconstructOldG) return;

  for(MInt cellId = 0; cellId < m_maxNoCells; cellId++) {
    if(a_isHalo(cellId)) continue;
    m_initialGCell[cellId] = 0;
    for(MInt b = 0; b < m_noBodiesToCompute; b++) {
      a_containingCell(cellId, b) = cellId;
      a_containingDomain(cellId, b) = domainId();
    }
  }


  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    if(a_level(cellId) == a_maxGCellLevel()) {
      m_initialGCell[cellId] = 1;
      MInt parentId = cellId;
      for(MInt level = a_level(cellId); level > minLevel(); level--) {
        parentId = c_parentId(parentId);
        m_initialGCell[parentId] = 1;
      }
    }
  }
}

/*
 * \brief: update list for rotating levelset
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::updateContainingGCells(MInt mode) {
  TRACE();

  MInt body;
  MInt set;

  if(mode == 1) {
    for(MInt cellId = 0; cellId < m_maxNoCells; cellId++) {
      for(MInt b = 0; b < m_noBodiesToCompute; b++) {
        body = m_bodiesToCompute[b];
        set = m_bodyToSetTable[body];
        if(!a_inBandG(cellId, set)) {
          a_containingCell(cellId, b) = -1;
          if(!m_reconstructOldG) {
            a_containingDomain(cellId, b) = -1;
          }
        }
      }
    }
  } else if(mode == 0) {
    if(m_reconstructOldG) {
      for(MInt cellId = 0; cellId < m_maxNoCells; cellId++) {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          body = m_bodiesToCompute[b];
          set = m_bodyToSetTable[body];
          if(a_inBandG(cellId, set)) {
            MInt contCell = a_containingCell(cellId, b);
            if(contCell > -1) {
              for(MInt level = grid().maxUniformRefinementLevel(); level < grid().maxRefinementLevel(); level++) {
                if(m_swapIds.find(contCell) != m_swapIds.end()) {
                  contCell = m_swapIds[contCell];
                }
              }
              a_containingCell(cellId, b) = contCell;
            }
          } else {
            a_containingCell(cellId, b) = -1;
          }
        }
      }
    } else {
      MInt domId;

      // send the size of the data set

      MIntScratchSpace noBandCells(grid().noDomains(), AT_, "noBandCells");
      noBandCells.fill(0);
      for(MInt b = 0; b < m_noBodiesToCompute; b++) {
        body = m_bodiesToCompute[b];
        set = m_bodyToSetTable[body];
        noBandCells[domainId()] += a_noBandCells(set);
      }
      MPI_Allreduce(MPI_IN_PLACE, &noBandCells[0], grid().noDomains(), MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                    "noBandCells[0]");

      MInt noCellsComm = m_swapIds.size();
      MIntScratchSpace noCellsToDom(grid().noDomains(), AT_, "noCellsToDom");
      noCellsToDom.fill(0);

      for(MInt i = 0; i < grid().noDomains(); i++) {
        if(noBandCells[i] > 0) {
          noCellsToDom[i] = noCellsComm;
        }
      }

      prepareGlobalComm(&noCellsToDom[0]);

      noCellsComm = mMax(m_globalSndOffsets[grid().noDomains()], m_globalRcvOffsets[grid().noDomains()]);


      MIntScratchSpace sndData(2 * noCellsComm, AT_, "sndData");
      MIntScratchSpace sndDataSize(grid().noDomains(), AT_, "sndDataSize");
      sndData.fill(-1);
      sndDataSize.fill(0);
      MIntScratchSpace rcvData(2 * noCellsComm, AT_, "rcvData");
      MIntScratchSpace rcvDataSize(grid().noDomains(), AT_, "rcvDataSize");
      rcvData.fill(-1);
      rcvDataSize.fill(0);

      std::vector<std::map<MInt, MInt>> swapIdsGlobal;
      swapIdsGlobal.resize(grid().noDomains());

      // send data
      for(std::map<MInt, MInt>::iterator it = m_swapIds.begin(); it != m_swapIds.end(); it++) {
        for(MInt d = 0; d < grid().noDomains(); d++) {
          if(domainId() == d) continue;
          if(noBandCells[d] <= 0) continue;
          sndData[2 * m_globalSndOffsets[d] + sndDataSize.p[d]] = it->first;
          sndDataSize.p[d]++;
          sndData[2 * m_globalSndOffsets[d] + sndDataSize.p[d]] = it->second;
          sndDataSize.p[d]++;
        }
      }

      exchangeBuffersGlobal(sndData.getPointer(), rcvData.getPointer(), sndDataSize.getPointer(),
                            rcvDataSize.getPointer(), &m_globalSndOffsets[0], &m_globalRcvOffsets[0], 9, 2);
      for(MInt i = 0; i < grid().noDomains(); i++) {
        MInt ind = 2 * m_globalRcvOffsets[i];
        for(MInt j = 0; j < rcvDataSize(i); j += 2) {
          swapIdsGlobal[i].insert(make_pair(rcvData[ind + j], rcvData[ind + j + 1]));
        }
      }

      for(std::map<MInt, MInt>::iterator it = m_swapIds.begin(); it != m_swapIds.end(); it++) {
        swapIdsGlobal[domainId()].insert(make_pair(it->first, it->second));
      }

      for(MInt cellId = 0; cellId < m_maxNoCells; cellId++) {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          body = m_bodiesToCompute[b];
          set = m_bodyToSetTable[body];
          if(a_inBandG(cellId, set)) {
            MInt contCell = a_containingCell(cellId, b);
            if(contCell > -1) {
              domId = a_containingDomain(cellId, b);
              if(swapIdsGlobal[domId].find(contCell) != swapIdsGlobal[domId].end()) {
                a_containingCell(cellId, b) = swapIdsGlobal[domId][contCell];
              }
            }
          } else {
            a_containingCell(cellId, b) = -1;
            a_containingDomain(cellId, b) = -1;
          }
        }
      }
    }
  }

  m_swapIds.clear();
}

/*
 * \brief: Updates a list which contains the respective windowCellId of each haloCell
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::copyWindowToHaloIds() {
  TRACE();

  MInt haloId, windowId;

  if(m_reconstructOldG) return;

  for(MInt i = 0; i < 2 * m_maxNoCells; i++) {
    m_cellDomIds[i] = -1;
  }

  MIntScratchSpace tmp_data(a_noCells() * 3, AT_, "tmp_data");
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noWindowCells(i); j++) {
      windowId = windowCellId(i, j);
      tmp_data[windowId * 3] = windowId;
      tmp_data[windowId * 3 + 1] = domainId();
      tmp_data[windowId * 3 + 2] = m_initialGCell[windowId];
    }
  }
  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        windowId = grid().azimuthalWindowCell(i, j);
        tmp_data[windowId * 3] = windowId;
        tmp_data[windowId * 3 + 1] = domainId();
        tmp_data[windowId * 3 + 2] = m_initialGCell[windowId];
      }
    }
  }

  // exchange data -> send, receive
  exchangeDataLS(&tmp_data[0], 3);

  // scatter:
  // update the halo and window cell lists
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      haloId = haloCellId(i, j);
      m_cellDomIds[haloId * 2 + 0] = tmp_data[haloId * 3];
      m_cellDomIds[haloId * 2 + 1] = tmp_data[haloId * 3 + 1];
      m_initialGCell[haloId] = tmp_data[haloId * 3 + 2];
    }
  }
  if(grid().azimuthalPeriodicity()) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
        haloId = grid().azimuthalHaloCell(i, j);
        m_cellDomIds[haloId * 2 + 0] = tmp_data[haloId * 3];
        m_cellDomIds[haloId * 2 + 1] = tmp_data[haloId * 3 + 1];
        m_initialGCell[haloId] = tmp_data[haloId * 3 + 2];
      }
    }
  }
}


/**
 * \brief   Check that the important properies and values in the the halo-Cells are set correctly!
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::checkHaloCells() {
  TRACE();

#if defined LS_DEBUG || !defined NDEBUG

  const MFloat eps0 = 1e-8;

  MInt noChecks = 3 * m_noSets;
  if(m_closeGaps) noChecks += 2;

  MFloatScratchSpace cellCheck(a_noCells(), noChecks, AT_, "cellCheck");
  cellCheck.fill(std::numeric_limits<MFloat>::max());

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    for(MInt set = 0; set < m_noSets; set++) {
      cellCheck(cellId, set * 3) = (MFloat)a_inBandG(cellId, set);
      cellCheck(cellId, set * 3 + 1) = a_levelSetFunctionG(cellId, set);
      if(!m_combustion) {
        cellCheck(cellId, set * 3 + 2) = (MFloat)a_bodyIdG(cellId, set);
      }
    }
    if(m_closeGaps) {
      cellCheck(cellId, 3 * m_noSets) = (MFloat)a_secondBodyId(cellId);
      cellCheck(cellId, 3 * m_noSets + 1) = (MFloat)a_nearGapG(cellId);
    }
  }

  exchangeData(&cellCheck(0), noChecks);

  for(MInt cellId = noInternalCells(); cellId < a_noCells(); cellId++) {
    // changes to avoid updateLowerGridLevels each TS
    if(!c_isLeafCell(cellId)) continue;

    // Since azimuthal periodic halos are not exact cartesian matches, these
    // checks are not valid for them
    if(grid().azimuthalPeriodicity() && grid().isPeriodic(cellId)) continue;

    for(MInt set = 0; set < m_noSets; set++) {
      ASSERT((MInt)cellCheck(cellId, set * 3) == a_inBandG(cellId, set), " " + to_string(set));
      ASSERT(fabs(cellCheck(cellId, set * 3 + 1) - a_levelSetFunctionG(cellId, set)) < eps0,
             to_string(a_levelSetFunctionG(cellId, set)) + " " + to_string(cellCheck(cellId, set * 3 + 1)));
      if(!m_combustion) {
        ASSERT((MInt)cellCheck(cellId, set * 3 + 2) == a_bodyIdG(cellId, set),
               to_string(cellCheck(cellId, set * 3 + 2)) + " " + to_string(a_bodyIdG(cellId, set)));
      }
    }
    if(m_closeGaps) {
      // ASSERT((MInt)cellCheck(cellId,3 * m_noSets) == a_secondBodyId( cellId ), to_string(cellCheck(cellId,3 *
      // m_noSets)) + " " + to_string(a_secondBodyId( cellId )));
      ASSERT((MInt)cellCheck(cellId, 3 * m_noSets + 1) == a_nearGapG(cellId),
             to_string(cellCheck(cellId, 3 * m_noSets + 1)) + " " + to_string(a_nearGapG(cellId)));
    }
  }

#endif
}


/**
 * \brief determines the value of 'data' in the given cell by recusively volumetric averaging among all its offsprings
 * \author Lennart Schneiders
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::reduceData(const MInt cellId, MFloat* data, const MInt dataBlockSize) {
  MFloat vol = cellVolumeAtCell(cellId);
  if(c_noChildren(cellId) > 0) {
    vol = F0;
    for(MInt d = 0; d < dataBlockSize; d++) {
      data[dataBlockSize * cellId + d] = F0;
    }
    for(MInt child = 0; child < IPOW2(nDim); child++) {
      MInt childId = c_childId(cellId, child);
      if(childId < 0) continue;
      if(c_noChildren(childId) == 0) continue;
      MFloat volc = reduceData(childId, data, dataBlockSize);
      for(MInt d = 0; d < dataBlockSize; d++) {
        data[dataBlockSize * cellId + d] += volc * data[dataBlockSize * childId + d];
      }
      vol += volc;
    }
    for(MInt d = 0; d < dataBlockSize; d++) {
      data[dataBlockSize * cellId + d] /= mMax(1e-14, vol);
    }
  }
  return vol;
}

// --------------------------------------------------------------------------------------


/**
 * \brief   sets the interfaceCell-array
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setInterfaceList(MIntScratchSpace& inList) {
  TRACE();

  inList.fill(0);

  const MInt startSet = m_reconstructBand > 0 ? 0 : m_startSet;
  const MInt endSet = m_noSets;


  // based on the levelset-values
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    if(inList[cellId] != 0) continue;

    for(MInt set = startSet; set < endSet; set++) {
      if(m_reconstructBand > 0 && !m_computeSet_backup[set]) {
        continue;
      }

      MBool addParents = false;
      if(approx(a_levelSetFunctionG(cellId, set), F0, MFloatEps)) {
        // TOBETESTED:
        // if(m_adaptationLevel < a_maxGCellLevel(set)) {
        if(a_level(cellId) < a_maxGCellLevel(set) && a_level(cellId) < this->m_maxSensorRefinementLevel[0]) {
          inList[cellId] = 1;
        }
        addParents = true;
      } else {
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          if(a_hasNeighbor(cellId, dir) > 0) {
            const MInt nghbrId = c_neighborId(cellId, dir);
            if((a_levelSetFunctionG(nghbrId, set) * a_levelSetFunctionG(cellId, set) < F0)) {
              // if(m_adaptationLevel < a_maxGCellLevel(set)) {
              if(a_level(cellId) < a_maxGCellLevel(set) && a_level(cellId) < this->m_maxSensorRefinementLevel[0]) {
                inList[cellId] = 1;
              }
              addParents = true;
              break;
            }
          }
        }
      }

      if(addParents) {
        MInt parentId = c_parentId(cellId);
        while(parentId > -1 && parentId < a_noCells()) {
          if(a_level(parentId) < this->m_maxSensorRefinementLevel[0]) {
            inList[parentId] = 1;
          }
          parentId = c_parentId(parentId);
        }
      }
    }
  }

  // Exchange the listCount on all Domains
#if defined LS_DEBUG || !defined NDEBUG
  if(m_reconstructBand > 0) {
    ASSERT(m_buildCollectedLevelSetFunction, "");
  }

  MInt listCount = 0;
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(inList[cellId] > 0) {
      ++listCount;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &listCount, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "listCount");

  if(listCount == 0 && !m_maxLevelChange && this->m_maxSensorRefinementLevel[0] > minLevel()) {
    mTerm(1, AT_, "No Cells found for refinement!");
  }
#endif

  // the levelset-solver has only one layer of haloCells, the exchange is necessary!
  exchangeDataLS(&inList[0]);
}

/**
 * \brief actually doing some pre-balance-stuff
 * \author Tim Wegmann, Sohel Herff
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::resetSolver() {
  /* Do something */
}

/**
 * \brief reset the solver during balancing
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::resetSolverFull() {
  for(MInt set = 0; set < m_noSets; set++)
    std::vector<MInt>().swap(m_bandCells[set]);
}

/**
 * \brief
 * \author Jannik Borgelt
 */
/// \brief Return data size to be communicated during DLB for a grid cell and given data id
template <MInt nDim>
MInt LsCartesianSolver<nDim>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  // Inactive ranks do not have any data to communicate
  if(!isActive()) {
    return 0;
  }

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0 || cellId >= noInternalCells()) {
    return 0;
  }

  MInt dataSize = 0;

  switch(dataId) {
    case 0: {
      dataSize = m_noSets;
      break;
    }
    case 1: {
      dataSize = 1;
      break;
    }
    case 2: {
      dataSize = (m_semiLagrange) ? m_noSets : 1;
      break;
    }
    case 3: {
      if(m_semiLagrange) {
        if(!m_reconstructOldG && m_LsRotate) {
          dataSize = 1;
        } else
          TERMM(1, "Unknown data id for !m_reconstructOldG && m_LsRotate.");
      } else {
        dataSize = nDim;
      }
      break;
    }
    case 4: {
      if(m_semiLagrange) {
        if(!m_reconstructOldG && m_LsRotate) {
          dataSize = 1;
        } else
          TERMM(1, "Unknown data id for !m_reconstructOldG && m_LsRotate.");
      } else {
        dataSize = nDim;
      }
      break;
    }
    default: {
      TERMM(1, "Unknown data id.");
      break;
    }
  }
  return dataSize;
}

/**
 * \brief
 * \author Jannik Borgelt
 */
/// \brief Store the solver data for a given data id ordered in the given buffer for DLB
template <MInt nDim>
void LsCartesianSolver<nDim>::getCellDataDlb(const MInt dataId, const MInt oldNoCells,
                                             const MInt* const bufferIdToCellId, MFloat* const data) {
  TRACE();

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    MInt dataSize = cellDataSizeDlb(dataId, gridCellId);

    switch(dataId) {
      case 0: {
        std::copy_n(&a_levelSetFunctionG(cellId, 0), dataSize, &data[localBufferId * dataSize]);
        break;
      }
      case 2: {
        if(m_semiLagrange) {
          std::copy_n(&a_oldLevelSetFunctionG(cellId, 0), dataSize, &data[localBufferId * dataSize]);
        } else {
          std::copy_n(&a_normalVectorG(cellId, 0, 0), dataSize, &data[localBufferId * dataSize]);
        }
        break;
      }
      case 3: {
        std::copy_n(&a_curvatureG(cellId, 0), dataSize, &data[localBufferId * dataSize]);
        break;
      }
      case 4: {
        std::copy_n(&a_levelSetFunctionSlope(cellId, 0, 0), dataSize, &data[localBufferId * dataSize]);
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
    localBufferId++;
  }
}

/**
 * \brief
 * \author Jannik Borgelt
 */
/// \brief Store the solver data for a given data id ordered in the given buffer for DLB
template <MInt nDim>
void LsCartesianSolver<nDim>::getCellDataDlb(const MInt dataId, const MInt oldNoCells,
                                             const MInt* const bufferIdToCellId, MInt* const data) {
  TRACE();

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    switch(dataId) {
      case 1: {
        data[localBufferId] = a_regridTriggerG(cellId);
        break;
      }
      case 3: {
        std::copy_n(&m_initialGCell[cellId], 1, &data[localBufferId]);
        break;
      }
      case 4: {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          data[localBufferId * m_noBodiesToCompute + b] = a_containingCell(cellId, b);
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
    localBufferId++;
  }
}

/**
 * \brief
 * \author Jannik Borgelt
 */
/// \brief Set the solver cell data after DLB
template <MInt nDim>
void LsCartesianSolver<nDim>::setCellDataDlb(const MInt dataId, const MFloat* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set the variables if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    switch(dataId) {
      case 0: {
        std::copy_n(data, noInternalCells() * m_noSets, &a_levelSetFunctionG(0, 0));
        break;
      }
      case 2: {
        if(m_semiLagrange) {
          std::copy_n(data, noInternalCells() * m_noSets, &a_oldLevelSetFunctionG(0, 0));
        } else {
          std::copy_n(data, noInternalCells(), &a_curvatureG(0, 0));
        }
        break;
      }
      case 3: {
        std::copy_n(data, noInternalCells() * nDim, &a_normalVectorG(0, 0, 0));
        break;
      }
      case 4: {
        std::copy_n(data, noInternalCells() * nDim, &a_levelSetFunctionSlope(0, 0, 0));
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
    }
  }
}

/**
 * \brief
 * \author Jannik Borgelt
 */
/// \brief Set the solver cell data after DLB
template <MInt nDim>
void LsCartesianSolver<nDim>::setCellDataDlb(const MInt dataId, const MInt* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set the variables if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    switch(dataId) {
      case 1: {
        for(MInt i = 0; i < noInternalCells(); i++) {
          a_regridTriggerG(i) = (MBool)data[i];
        }
        break;
      }
      case 3: {
        std::copy_n(data, noInternalCells(), &m_initialGCell[0]);
        break;
      }
      case 4: {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          for(MInt i = 0; i < noInternalCells(); i++) {
            a_containingCell(i, b) = data[b * noInternalCells() + i];
          }
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
    }
  }
}

/**
 * \brief
 * \author Jannik Borgelt
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::balancePre() {
  TRACE();

  // Set reinitialization stage
  m_loadBalancingReinitStage = 0;

  // Update the grid proxy for this solver
  grid().update();

  if(!grid().isActive()) {
    // Reset parallelization information if solver is not active
    updateDomainInfo(-1, -1, MPI_COMM_NULL, AT_);
  } else {
    // Set new domain info for solver
    updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);
  }

  // Reset cell, surface, boundary cell data and deallocate halo/window cell arrays
  resetSolverFull();

  // Reset all cells
  m_cells.clear();

  // Return if solver is not active
  if(!grid().isActive()) {
    return;
  }

  grid().updateLeafCellExchange();

  // check for empry cell collector
  ASSERT(m_cells.size() == 0, "");

  // Resize cell collector to internal cells
  m_cells.append(grid().noInternalCells());

  // Check that global ids are sorted
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    if(grid().domainOffset(domainId()) + (MLong)cellId != c_globalId(cellId)) {
      TERMM(1, "Global id mismatch.");
    }
    m_cells.erase(cellId);
  }
}


/**
 * \brief
 * \author Jannik Borgelt
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::balancePost() {
  TRACE();

  m_loadBalancingReinitStage = 1;

  // If m_buildCollectedLevelSetFunction = false, identifyBodies() is never called.
  // Setting bodyId to -1, however, leads to errors in cut-cell generation.
  MInt defaultBodyId = (m_buildCollectedLevelSetFunction ? -1 : 0);
  for(MInt i = 0; i < noInternalCells(); i++) {
    for(MInt j = 0; j < m_noSets; j++) {
      if(m_semiLagrange) {
        a_bodyIdG(i, j) = defaultBodyId;
      }
    }
  }

  // append the halo-cells and erase
  m_cells.append(grid().tree().size() - m_cells.size());
  for(MInt cellId = grid().noInternalCells(); cellId < a_noCells(); cellId++) {
    m_cells.erase(cellId);
  }

  testCellsCG();

  // Nothing to do if solver is not active
  if(!grid().isActive()) {
    return;
  }

  // reallocate exchange-storages:
  if(grid().noDomains() > 1) {
    if(m_combustion) {
      mDeallocate(m_intSendBuffers);
      mDeallocate(m_intReceiveBuffers);

      mAlloc(m_intSendBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intSendBuffers", 0, AT_);
      mAlloc(m_intReceiveBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intReceiveBuffers", 0, AT_);
    }

    if(!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2) {
      mDeallocate(m_gSendBuffers);
      mDeallocate(m_gReceiveBuffers);

      mAlloc(m_gSendBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gSendBuffers", F0, AT_);
      mAlloc(m_gReceiveBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gReceiveBuffers", F0, AT_);
    }

    if(m_combustion || (!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2)) {
      mDeallocate(mpi_request);
      mDeallocate(mpi_recive);

      mAlloc(mpi_request, grid().noNeighborDomains(), "mpi_request", AT_);
      mAlloc(mpi_recive, grid().noNeighborDomains(), "mpi_recive", AT_);
    }
  }

  generateListOfGExchangeCellsCG();
  this->checkNoHaloLayers();

  // initAzimuthalExchange
  initAzimuthalExchange();

  exchangeAllLevelSetData();

  m_adaptationSinceLastRestart = true;
  m_loadBalancingReinitStage = 2;
}

/**
 * \brief
 * \author Sohel Herff, Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::balance(const MInt* const noCellsToReceiveByDomain,
                                      const MInt* const noCellsToSendByDomain, const MInt* const sortedCellId,
                                      const MInt oldNoCells) {
  TRACE();

  NEW_TIMER_GROUP(t_initTimer, "balance solver");
  NEW_TIMER(t_timertotal, "balance solver", t_initTimer);
  NEW_SUB_TIMER(t_variables, "variables", t_timertotal);
  NEW_SUB_TIMER(t_communicator, "communicator", t_timertotal);

  RECORD_TIMER_START(t_timertotal);
  RECORD_TIMER_START(t_variables);

  // save solver-data in grid-format (necessary for multi-solver-balancing!)
  // while this might not be the fastes way to ensure a multi-solver balance,
  // it is ginuelly simple and clear to understand what is happening!
  // This needs to happen before the proxy-update!! (otherwise the solver2grid might change)
  MFloatScratchSpace lsValuesBalance(oldNoCells, m_noSets, FUN_, "lsValuesBalance");
  MIntScratchSpace regridBalance(oldNoCells, FUN_, "regridBalance");
  MFloatScratchSpace curvatureBalance((!m_semiLagrange) * oldNoCells, FUN_, "curvatureBalance");
  MFloatScratchSpace normalVectorsBalance((!m_semiLagrange) * oldNoCells, nDim, FUN_, "normalVectorsBalance");
  MFloatScratchSpace slopeBalance((!m_semiLagrange) * oldNoCells, nDim, FUN_, "slopeBalance");

  MFloatScratchSpace oldLsValuesBalance(m_semiLagrange * oldNoCells, m_semiLagrange * m_noSets, FUN_,
                                        "oldLsValuesBalance");

  MLongScratchSpace containingCellsBalance((m_LsRotate && !m_reconstructOldG) * oldNoCells,
                                           (m_LsRotate && !m_reconstructOldG) * m_noBodiesToCompute, FUN_,
                                           "containingCellsBalance");
  MLongScratchSpace initialGCellBalance((m_LsRotate && !m_reconstructOldG) * oldNoCells, FUN_, "initialGCellBalance");

  regridBalance.fill(-2);
  lsValuesBalance.fill(-900);
  curvatureBalance.fill(-900);
  normalVectorsBalance.fill(-900);
  slopeBalance.fill(-900);


  oldLsValuesBalance.fill(-900);
  containingCellsBalance.fill(-900);
  initialGCellBalance.fill(-900);

  for(MInt cellId = 0; cellId < grid().tree().size(); cellId++) {
    MInt gridCellId = grid().tree().solver2grid(cellId);
    for(MInt set = 0; set < m_noSets; set++) {
      lsValuesBalance(gridCellId, set) = a_levelSetFunctionG(cellId, set);
    }
    if(a_regridTriggerG(cellId)) regridBalance(gridCellId) = 1;
  }

  if(m_semiLagrange) {
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      MInt gridCellId = grid().tree().solver2grid(cellId);
      for(MInt set = 0; set < m_noSets; set++) {
        oldLsValuesBalance(gridCellId, set) = a_oldLevelSetFunctionG(cellId, set);
      }
      if(!m_reconstructOldG && m_LsRotate) {
        for(MInt b = 0; b < m_noBodiesToCompute; b++) {
          containingCellsBalance(gridCellId, b) = a_containingCell(cellId, b);
        }
        initialGCellBalance(gridCellId) = m_initialGCell[cellId];
      }
    }

  } else {
    for(MInt cellId = 0; cellId < grid().tree().size(); cellId++) {
      MInt gridCellId = grid().tree().solver2grid(cellId);
      for(MInt dim = 0; dim < nDim; dim++) {
        normalVectorsBalance(gridCellId, dim) = a_normalVectorG(cellId, dim, 0);
        slopeBalance(gridCellId, dim) = a_levelSetFunctionSlope(cellId, dim, 0);
      }
      curvatureBalance(gridCellId) = a_curvatureG(cellId, 0);
    }
  }


  // update of the proxy
  grid().update();

  // Just reset parallelization information if solver is not active
  // NOTE: inactive ranks not supported here since the data sizes etc are not determined for the
  // solver but for the whole grid -> use balancePre/Post instead!
  if(!isActive()) {
    updateDomainInfo(-1, -1, MPI_COMM_NULL, AT_);
    TERMM(1, "fixme: inactive ranks not supported in balance(); implement balancePre/Post!");
    return;
  }

  // Set new domain info for solver
  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);

  // This is only working of the same domains as in the grid are used!
  ASSERT(domainId() == grid().raw().domainId(), "");
  ASSERT(noDomains() == grid().raw().noDomains(), "");

  // data-to be saved during the balancing:
  //- a_levelSetFunctionG
  //- regridTrigger
  // if semiLagrange:
  //- a_bodyIdG
  //- oldLevelSetFunction

  MFloatScratchSpace levelSet(noCellsToReceiveByDomain[noDomains()], m_noSets, FUN_, "levelSet");
  MIntScratchSpace regridTrigger(noCellsToReceiveByDomain[noDomains()], FUN_, "regridTrigger");

  MFloatScratchSpace curvature((!m_semiLagrange) * noCellsToReceiveByDomain[noDomains()], FUN_, "curvature");
  MFloatScratchSpace normalVectors((!m_semiLagrange) * noCellsToReceiveByDomain[noDomains()], nDim, FUN_,
                                   "normalVectors");
  MFloatScratchSpace slopes((!m_semiLagrange) * noCellsToReceiveByDomain[noDomains()], nDim, FUN_, "slopes");

  MFloatScratchSpace oldLevelSet(m_semiLagrange * noCellsToReceiveByDomain[noDomains()], m_semiLagrange * m_noSets,
                                 FUN_, "oldLevelSet");

  MLongScratchSpace containingCells((m_LsRotate && !m_reconstructOldG) * noCellsToReceiveByDomain[noDomains()],
                                    (m_LsRotate && !m_reconstructOldG) * m_noBodiesToCompute, FUN_, "containingCells");
  MLongScratchSpace initialGCell((m_LsRotate && !m_reconstructOldG) * noCellsToReceiveByDomain[noDomains()], FUN_,
                                 "initialCell");

  levelSet.fill(-1.0);
  regridTrigger.fill(-2);
  slopes.fill(-1.0);
  curvature.fill(-1.0);
  normalVectors.fill(-1.0);


  oldLevelSet.fill(-1.0);
  containingCells.fill(-1);
  initialGCell.fill(-1);

  maia::mpi::communicateData(&lsValuesBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, m_noSets, &levelSet[0]);

  maia::mpi::communicateData(&regridBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &regridTrigger[0]);

  if(m_semiLagrange) {
    maia::mpi::communicateData(&oldLsValuesBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                               noCellsToSendByDomain, noCellsToReceiveByDomain, m_noSets, &oldLevelSet[0]);

    if(!m_reconstructOldG && m_LsRotate) {
      maia::mpi::communicateData(&containingCellsBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(),
                                 mpiComm(), noCellsToSendByDomain, noCellsToReceiveByDomain, m_noBodiesToCompute,
                                 &containingCells[0]);
      maia::mpi::communicateData(&initialGCellBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                                 noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &initialGCell[0]);
    }

  } else {
    maia::mpi::communicateData(&normalVectorsBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                               noCellsToSendByDomain, noCellsToReceiveByDomain, nDim, &normalVectors[0]);
    maia::mpi::communicateData(&slopeBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                               noCellsToSendByDomain, noCellsToReceiveByDomain, nDim, &slopes[0]);
    maia::mpi::communicateData(&curvatureBalance[0], oldNoCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                               noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &curvature[0]);
  }

  // reset all data!
  resetSolverFull();

  // Reset all cells
  m_cells.clear();

  // Iterate over all received grid cells
  for(MInt gridCellId = 0; gridCellId < noCellsToReceiveByDomain[noDomains()]; gridCellId++) {
    // Determine cell id and append to collector
    const MInt cellId = m_cells.size();

    if(!grid().raw().treeb().solver(gridCellId, m_solverId)) continue;

    ASSERT(grid().tree().solver2grid(cellId) == gridCellId, "");

    if(grid().domainOffset(domainId()) + (MLong)cellId != c_globalId(cellId)) {
      mTerm(1, AT_, "Global id mismatch.");
    }
    m_cells.append();
    m_cells.erase(cellId);

    // Set solver data
    for(MInt j = 0; j < m_noSets; j++) {
      a_levelSetFunctionG(cellId, j) = levelSet(gridCellId, j);

      if(m_semiLagrange) {
        a_bodyIdG(cellId, j) = -1;
        a_oldLevelSetFunctionG(cellId, j) = oldLevelSet(gridCellId, j);
      }
    }

    if(!m_reconstructOldG && m_LsRotate) {
      for(MInt b = 0; b < m_noBodiesToCompute; b++) {
        a_containingCell(cellId, b) = containingCells(gridCellId, b);
      }
      m_initialGCell[cellId] = initialGCell[gridCellId];
    }

    if(!m_semiLagrange) {
      for(MInt dim = 0; dim < nDim; dim++) {
        a_normalVectorG(cellId, dim, 0) = normalVectors(gridCellId, dim);
        a_levelSetFunctionSlope(cellId, dim, 0) = slopes(gridCellId, dim);
      }
      a_curvatureG(cellId, 0) = curvature(gridCellId);
    }

    if(regridTrigger[gridCellId] > 0) {
      a_regridTriggerG(cellId) = true;
    }
  }

  // append the halo-cells:
  m_cells.append(grid().tree().size() - m_cells.size());
  for(MInt cellId = noInternalCells(); cellId < a_noCells(); cellId++) {
    m_cells.erase(cellId);
  }

  testCellsCG();

  RECORD_TIMER_STOP(t_variables);
  RECORD_TIMER_START(t_communicator);

  // reallocate exchange-storages:
  if(grid().noDomains() > 1) {
    if(m_combustion) {
      mDeallocate(m_intSendBuffers);
      mDeallocate(m_intReceiveBuffers);

      mAlloc(m_intSendBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intSendBuffers", 0, AT_);
      mAlloc(m_intReceiveBuffers, grid().noNeighborDomains(), m_maxNoCells, "m_intReceiveBuffers", 0, AT_);
    }

    if(!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2) {
      mDeallocate(m_gSendBuffers);
      mDeallocate(m_gReceiveBuffers);

      mAlloc(m_gSendBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gSendBuffers", F0, AT_);
      mAlloc(m_gReceiveBuffers, grid().noNeighborDomains(), m_maxNoSets * m_maxNoCells, "m_gReceiveBuffers", F0, AT_);
    }

    if(m_combustion || (!m_semiLagrange || m_guaranteeReinit || m_STLReinitMode != 2)) {
      mDeallocate(mpi_request);
      mDeallocate(mpi_recive);

      mAlloc(mpi_request, grid().noNeighborDomains(), "mpi_request", AT_);
      mAlloc(mpi_recive, grid().noNeighborDomains(), "mpi_recive", AT_);
    }
  }

  generateListOfGExchangeCellsCG();

  this->checkNoHaloLayers();

  startLoadTimer(AT_);

  // Initialize the azimuthal periodic exchange
  initAzimuthalExchange();

  exchangeAllLevelSetData();

  m_adaptationSinceLastRestart = true;

  RECORD_TIMER_STOP(t_communicator);
  RECORD_TIMER_STOP(t_timertotal);
  DISPLAY_TIMER(t_timertotal);
}

/**
 * \brief re-inits the solver after the balance
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::finalizeBalance() {
  if(!grid().isActive()) {
    return;
  }

  reInitSolver(false);

  // Update local ids of window cell on halo rank or reconstructOldG
  if(m_LsRotate) {
    if(m_reconstructOldG) {
      reconstructOldGField();
    } else {
      globalToLocalIdsContainingCells();
      copyWindowToHaloIds();
    }
  }
}


/**
 * \brief rotates a coordinate
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::rotateLevelSet(MInt returnMode, MFloat* cellData, MInt body, const MFloat* xCoord,
                                             const MFloat* xCenter, const MFloat* angle) {
  TRACE();

  IF_CONSTEXPR(nDim == 2) mTerm(1, AT_, "rotateLevelSet needs to be updated to 2D!");

  switch(returnMode) {
    case 1: // current coordinate
      cellData[0] = xCenter[0] + (cos(angle[1]) * cos(angle[0])) * (xCoord[0] - xCenter[0])
                    + (cos(angle[1]) * sin(angle[0])) * (xCoord[1] - xCenter[1])
                    + (-sin(angle[1])) * (xCoord[2] - xCenter[2]);

      cellData[1] =
          xCenter[1]
          + (sin(angle[2]) * sin(angle[1]) * cos(angle[0]) - cos(angle[2]) * sin(angle[0])) * (xCoord[0] - xCenter[0])
          + (sin(angle[2]) * sin(angle[1]) * sin(angle[0]) + cos(angle[2]) * cos(angle[0])) * (xCoord[1] - xCenter[1])
          + (sin(angle[2]) * cos(angle[1])) * (xCoord[2] - xCenter[2]);

      cellData[2] =
          xCenter[2]
          + (cos(angle[2]) * sin(angle[1]) * cos(angle[0]) + sin(angle[2]) * sin(angle[0])) * (xCoord[0] - xCenter[0])
          + (cos(angle[2]) * sin(angle[1]) * sin(angle[0]) - sin(angle[2]) * cos(angle[0])) * (xCoord[1] - xCenter[1])
          + (cos(angle[2]) * cos(angle[1])) * (xCoord[2] - xCenter[2]);

      break;

    case 2: // current velocity

      cellData[0] =
          ((m_omega[body * nDim + 1] * cos(angle[0]) + m_omega[body * nDim + 2] * sin(angle[0]) * cos(angle[1]))
               * (xCoord[2] - xCenter[2])
           - (m_omega[body * nDim + 0] - m_omega[body * nDim + 2] * sin(angle[1])) * (xCoord[1] - xCenter[1]));
      cellData[1] =
          ((m_omega[body * nDim + 0] - m_omega[body * nDim + 2] * sin(angle[1])) * (xCoord[0] - xCenter[0])
           - (-m_omega[body * nDim + 1] * sin(angle[0]) + m_omega[body * nDim + 2] * cos(angle[0]) * cos(angle[1]))
                 * (xCoord[2] - xCenter[2]));
      cellData[2] =
          ((-m_omega[body * nDim + 1] * sin(angle[0]) + m_omega[body * nDim + 2] * cos(angle[0]) * cos(angle[1]))
               * (xCoord[1] - xCenter[1])
           - (m_omega[body * nDim + 1] * cos(angle[0]) + m_omega[body * nDim + 2] * sin(angle[0]) * cos(angle[1]))
                 * (xCoord[0] - xCenter[0]));

      break;

    case 3: // current acceleration

      MFloat omeg[3], omegRad[3];

      rotateLevelSet(5, omeg, body, nullptr, nullptr, &m_semiLagrange_xRot_STL[body * nDim]);
      rotateLevelSet(2, omegRad, body, xCoord, xCenter, &m_semiLagrange_xRot_STL[body * nDim]);

      cellData[0] = omeg[1] * omegRad[2] - omeg[2] * omegRad[1];
      cellData[1] = omeg[2] * omegRad[0] - omeg[0] * omegRad[2];
      cellData[2] = omeg[0] * omegRad[1] - omeg[1] * omegRad[0];

      break;

    case 4: // current angle

      cellData[0] = angle[0];
      cellData[1] = angle[1];
      cellData[2] = angle[2];

      break;

    case 5: // current angular velocity

      cellData[0] =
          (-m_omega[body * nDim + 1] * sin(angle[0]) + m_omega[body * nDim + 2] * cos(angle[0]) * cos(angle[1]));
      cellData[1] =
          (m_omega[body * nDim + 1] * cos(angle[0]) + m_omega[body * nDim + 2] * sin(angle[0]) * cos(angle[1]));
      cellData[2] = (m_omega[body * nDim + 0] - m_omega[body * nDim + 2] * sin(angle[1]));

      break;

    case 6: // current angular acceleration

      cellData[0] = F0;
      cellData[1] = F0;
      cellData[2] = F0;

      break;

    default:
      cellData[0] = F0;
      cellData[1] = F0;
      cellData[2] = F0;
  }
}

/**
 * \brief searches halo cells for point
 *
 * \author Thomas Hoesgen
 * \date 01/2019
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::getContainingCellHalo(MFloat* point) {
  TRACE();

  MInt cellId = -1;

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      cellId = haloCellId(i, j);
      if(a_level(cellId) != a_maxGCellLevel() || (!m_reconstructOldG && m_initialGCell[cellId] == 0)) continue;

      if(inCell(cellId, point)) {
        cerr << "Halo " << cellId << " " << c_coordinate(cellId, 0) << " " << c_coordinate(cellId, 1) << " "
             << c_coordinate(cellId, 2) << endl;

        return cellId;
      }
    }
  }

  return -1;
}


/**
 * \brief computes levelset value and processes containingCell
 *
 * \author Thomas Hoesgen
 * \date 01/2019
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::processRotatingLevelSet(MFloat& phi, MInt& cellId, MInt& domId, MFloat* point, MInt set) {
  TRACE();

  if(a_level(cellId) != a_maxGCellLevel() || (!m_reconstructOldG && m_initialGCell[cellId] == 0)) {
    mTerm(1, AT_, "ContainingCell is not initialGCell or on maxLevel!");
  } else {
    const MInt magic_number = 8; // pow(2, nDim);
    std::array<MInt, magic_number> interpolationCells = {0, 0, 0, 0, 0, 0, 0, 0};
    MInt position = 0;
    // Set up interpolation stencil
    position = setUpLevelSetInterpolationStencil(cellId, interpolationCells.data(), point);
    // Check if all interpolationCells are initialGCells
    if(position > -1) {
      for(MInt i = 0; i < 8; i++) {
        if(a_level(interpolationCells[i]) != a_maxGCellLevel()
           || (!m_reconstructOldG && m_initialGCell[interpolationCells[i]] == 0)) {
          mTerm(1, AT_, "interpolationCell is not initialGCell!");
          position = -1;
          break;
        }
      }
    }
    // Interpolate level set
    if(position > -1) {
      phi = interpolateOldLevelSet(interpolationCells.data(), point, set);
    } else {
      phi = a_oldLevelSetFunctionG(cellId, set);
    }

    if(a_isHalo(cellId)) {
      if(!m_reconstructOldG) {
        domId = m_cellDomIds[cellId * 2 + 1];
        cellId = m_cellDomIds[cellId * 2 + 0];
      }
    }
  }
}

/**
 * \brief  checks if a child lies outSide of the domain!
 *         necessary for refinement at the bndry!
 * \author Tim Wegmann
 */

template <MInt nDim>
MInt LsCartesianSolver<nDim>::cellOutside(const MFloat* coords, const MInt level, const MInt gridCellId) {
  // TODO labels:LS,toremove remove and update fv-combustion testcases accordingly!
  if(m_combustion) return -1;

  if(m_engineSetup) {
    return -1;
  }

  if(m_virtualSurgery) {
    return -1;
  }

  std::ignore = gridCellId;

  static constexpr MInt cornerIndices[8][3] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
                                               {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
  MFloat corner[3] = {0, 0, 0};
  MBool outside = true;
  MFloat cellHalfLength = F1B2 * c_cellLengthAtLevel(level);

  for(MInt i = 0; i < m_noCorners; i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      corner[dim] = coords[dim] + cornerIndices[i][dim] * cellHalfLength;
    }
    IF_CONSTEXPR(nDim == 2) {
      if(!m_geometry->pointIsInside(corner)) outside = false;
      // pointIsInside == true if Point is outside fluid domain
    }
    else {
      if(!m_geometry->pointIsInside2(corner)) outside = false;
      // pointIsInside == true if Point is outside fluid domain
    }
  }

  // Why? pointIsInside checks for the outer geometry
  // if(m_levelSetSign[0] < 0) {
  // outside = !outside;
  //}

  return outside;
}

/**
 * \brief  Return the number of Ls load types.
 *
 * Type 1: number of band-cells
 * Type 2: number of G0-cells
 *
 * \author Tim Wegmann
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::noLoadTypes() const {
  TRACE();
  // band-Cells and G0-Cells
  const MInt noLsLoadTypes = m_weightLevelSet ? 3 : 1;
  return noLsLoadTypes;
}


/// \brief Return the default weights for all load quantities
template <MInt nDim>
void LsCartesianSolver<nDim>::getDefaultWeights(MFloat* weights, std::vector<MString>& names) const {
  TRACE();

  // TODO labels:LS set sensible default values
  weights[0] = 0.1;
  names[0] = "ls_cell";
  MInt count = 1;

  if(m_weightLevelSet) {
    weights[1] = 0.1;
    names[1] = "ls_band_cell";
    count++;

    weights[2] = 0.5;
    names[2] = "ls_g0_cell";
    count++;
  }

  if(noLoadTypes() != count) {
    TERMM(1, "Count does not match noLoadTypes.");
  }
}


/**
 * \brief  Return the cumulative load quantities on this domain.
 *
 * \param[out] loadQuantities Storage for load quantities.
 *
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::getLoadQuantities(MInt* const loadQuantities) const {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // reset
  for(MInt type = 0; type < noLoadTypes(); type++) {
    loadQuantities[type] = 0;
  }


  loadQuantities[0] = noInternalCells();

  MInt noBandCells = 0;
  MInt noG0Cells = 0;
  if(m_weightLevelSet) {
    for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
      MBool inAnyBand = false;
      for(MInt set = (MInt)m_buildCollectedLevelSetFunction; set < m_noSets; set++) {
        if(a_inBandG(cellId, set) && !inAnyBand) {
          noBandCells++;
          inAnyBand = true;
        }
        if(a_isGZeroCell(cellId, set)) {
          noG0Cells++;
          break;
        }
      }
    }
  }

  loadQuantities[1] = noBandCells;
  loadQuantities[2] = noG0Cells;
}


/**
 * \brief  Return the load of a single cell (given computational weights).
 *
 * \param[in] cellId Requested grid cell id.
 * \param[in] weights Computational weights for different simulation components.
 * \return Cell load.
 *
 * \author Tim Wegmann
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::getCellLoad(const MInt gridCellId, const MFloat* const weights) const {
  TRACE();
  ASSERT(isActive(), "solver is not active");

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0) {
    return 0;
  }

  if(cellId < 0 || cellId >= grid().noInternalCells()) {
    TERMM(1, "The given cell id is invalid.");
  }

  // Default cell load
  MFloat cellLoad = weights[0];

  if(m_weightLevelSet) {
    for(MInt set = (MInt)m_buildCollectedLevelSetFunction; set < m_noSets; set++) {
      if(a_inBandG(cellId, set)) {
        cellLoad = weights[1];
        if(a_isGZeroCell(cellId, set)) {
          cellLoad = weights[2];
          return cellLoad;
        }
      }
    }
  }
  return cellLoad;
}

/**
 * \brief  Limit weight of base cell
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::limitWeights(MFloat* weights) {
  if(m_limitWeights) {
    weights[0] = mMax(weights[0], 0.01 * mMax(weights[1], weights[2]));
  }
}

/**
 * \brief  Get solver timings
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings,
                                               const MBool NotUsed(allTimings)) {
  TRACE();

  const MString namePrefix = "b" + std::to_string(solverId()) + "_";

  const MFloat load = returnLoadRecord();
  const MFloat idle = returnIdleRecord();

  solverTimings.emplace_back(namePrefix + "loadLsCartesianSolver", load);
  solverTimings.emplace_back(namePrefix + "idleLsCartesianSolver", idle);

#ifdef MAIA_TIMER_FUNCTION
  solverTimings.emplace_back(namePrefix + "timeIntegration", RETURN_TIMER_TIME(m_timers[Timers::TimeInt]));
  solverTimings.emplace_back(namePrefix + "firstEx", RETURN_TIMER_TIME(m_timers[Timers::FirstEx]));
  solverTimings.emplace_back(namePrefix + "postTS", RETURN_TIMER_TIME(m_timers[Timers::PostTime]));
  solverTimings.emplace_back(namePrefix + "finalize", RETURN_TIMER_TIME(m_timers[Timers::Finalize]));
  solverTimings.emplace_back(namePrefix + "buildTube", RETURN_TIMER_TIME(m_timers[Timers::BuildTube]));
  solverTimings.emplace_back(namePrefix + "setBand", RETURN_TIMER_TIME(m_timers[Timers::SetBand]));
  solverTimings.emplace_back(namePrefix + "multiple", RETURN_TIMER_TIME(m_timers[Timers::BuildMultiple]));
#endif
}


/**
 * \brief  Return decomposition information, i.e. number of local elements,...
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) {
  TRACE();

  MInt noBandCells = 0;
  MInt noG0Cells = 0;
  MInt noCells = noInternalCells();
  if(m_weightLevelSet) {
    for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
      MBool inAnyBand = false;
      for(MInt set = (MInt)m_buildCollectedLevelSetFunction; set < m_noSets; set++) {
        if(a_inBandG(cellId, set) && !inAnyBand) {
          noBandCells++;
          inAnyBand = true;
        }
        if(a_isGZeroCell(cellId, set)) {
          noG0Cells++;
          break;
        }
      }
    }
  }

  const MString namePrefix = "b" + std::to_string(solverId()) + "_";
  domainInfo.emplace_back(namePrefix + "noLsCells", noCells);
  domainInfo.emplace_back(namePrefix + "noLsBandCells", noBandCells);
  domainInfo.emplace_back(namePrefix + "noLsG0Cells", noG0Cells);
}

/**
 * \brief returns the cellId of the cell containing point in the second layer neighbor cells.
 *        Is to be called from getContainingCell
 * \author Thomas Hoesgen
 * \date 03/2019
 */
template <MInt nDim>
MInt LsCartesianSolver<nDim>::checkSecondLayerCells(std::vector<MInt>& diag2Cells,
                                                    std::map<MInt, std::vector<MInt>>& dirCode, MFloat* point) {
  TRACE();

  MInt nghbrId = -1;

  std::vector<MInt> secondLayer;
  MInt d[3];
  MInt nghbrId2, nghbrId3;
  MInt d1, d2, d3, e, w; //,g;
  MInt dir[3] = {0, 2, 4};
  MInt oDir[6] = {2, 4, 0, 4, 0, 2};

  for(MInt& diag2Cell : diag2Cells) {
    d[0] = dirCode[diag2Cell][0];
    d[1] = dirCode[diag2Cell][1];
    d[2] = dirCode[diag2Cell][2];

    // add neighbors
    for(MInt nghbr = 0; nghbr < 2 * nDim; nghbr++) {
      if(nghbr == d[nghbr / 2]) continue;
      nghbrId = c_neighborId(diag2Cell, nghbr);
      if(nghbrId == -1) continue;

      secondLayer.push_back(nghbrId);
    }

    // add diagonal neighbors
    for(MInt i = 0; i < nDim; i++) {
      for(MInt j = 0; j < 2; j++) {
        d1 = dir[i] + j;
        nghbrId = c_neighborId(diag2Cell, d1);
        if(nghbrId == -1) continue;
        e = d1 / 2;
        for(MInt k = 0; k < 2; k++) {
          for(MInt m = 0; m < 2; m++) {
            d2 = oDir[2 * e + k] + m;
            nghbrId2 = c_neighborId(nghbrId, d2);
            if(nghbrId2 == -1) continue;
            w = d2 / 2;
            // if ( !(d1 == d[i]) && !(d2 == d[w]) )
            secondLayer.push_back(nghbrId2);

            for(MInt n = 0; n < 2; n++) {
              d3 = dir[nDim - e - w] + n;
              if(a_hasNeighbor(nghbrId2, d3)) {
                nghbrId3 = c_neighborId(nghbrId2, d3);

                // g = d3 / 2;
                // if ( !(d1 == d[i]) && !(d2 == d[w]) && !(d3 == d[g] ) )
                secondLayer.push_back(nghbrId3);
              }
            }
          }
        }
      }
    }
  }

  std::sort(secondLayer.begin(), secondLayer.end());
  auto last = std::unique(secondLayer.begin(), secondLayer.end());
  secondLayer.erase(last, secondLayer.end());
  for(MInt& it : secondLayer) {
    if(inCell(it, point)) return it;
  }

  return -1;
}


/**
 * \brief Determines offsets for global communication to get minimum memory requirements
 * \author Thomas Hoesgen
 * \date 05/2019
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::prepareGlobalComm(MInt* noCellsToDom) {
  TRACE();

  ScratchSpace<MPI_Request> requestGlobal(grid().noDomains(), AT_, "requestGlobal");
  MPI_Status status;
  MInt tag = 3;

  requestGlobal.fill(MPI_REQUEST_NULL);

  m_globalSndOffsets[0] = 0;
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) {
      m_globalSndOffsets[i + 1] = m_globalSndOffsets[i];
    } else {
      m_globalSndOffsets[i + 1] = m_globalSndOffsets[i] + noCellsToDom[i];
    }
  }


  MIntScratchSpace noCellsToSend(grid().noDomains(), AT_, "sendData");
  MIntScratchSpace noCellsToReceive(grid().noDomains(), AT_, "receiveData");
  noCellsToSend.fill(0);
  noCellsToReceive.fill(0);

  m_globalRcvOffsets[0] = 0;
  // send the size of the data set
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    noCellsToSend.p[i] = noCellsToDom[i];
    MPI_Issend(&noCellsToSend.p[i], 1, MPI_INT, i, tag, mpiComm(), &requestGlobal[i], AT_, "noCellsToSend.p[i]");
  }

  // receive the size of the data set
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) {
      m_globalRcvOffsets[i + 1] = m_globalRcvOffsets[i];
    } else {
      MPI_Recv(&noCellsToReceive.p[i], 1, MPI_INT, i, tag, mpiComm(), &status, AT_, "noCellsToReceive.p[i]");
      m_globalRcvOffsets[i + 1] = m_globalRcvOffsets[i] + noCellsToReceive.p[i];
    }
  }
  for(MInt i = 0; i < grid().noDomains(); i++) {
    if(domainId() == i) continue;
    MPI_Wait(&requestGlobal[i], &status, AT_);
  }
}

/** \page sensorsLS
 *
 *  \section interface Interface
 *
 *  This sensors ensures a band of refined cells around the level set.<br>
 *  Property: <code>INTERFACE</code>
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::sensorInterface(std::vector<std::vector<MFloat>>& sensors,
                                              std::vector<std::bitset<64>>& sensorCellFlag,
                                              std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) {
  m_log << "   - Sensor preparation for the interface sensor" << endl;

  MIntScratchSpace inList(a_noCells(), AT_, "inList");
  setInterfaceList(inList);

  for(MInt level = minLevel(); level < this->m_maxSensorRefinementLevel[sen]; level++) {
    this->markSurrndCells(inList, m_outerBandWidth[level], level, m_refineDiagonals);
  }

  ASSERT(a_noCells() == grid().tree().size(), "");

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    if(inList(cellId) == 0) {
      // if(a_oldLevelSetFunctionG(cellId , 0)  > - m_outsideGValue  &&
      //   a_oldLevelSetFunctionG(cellId , 0)  < m_outsideGValue ) continue;
      // keep the are of the old-levelSet-function refined!
      if(a_level(cellId) == minLevel()) continue;
      if(inList(c_parentId(cellId))) continue;
      if(c_noChildren(cellId) == 0) {
        const MInt gridCellId = grid().tree().solver2grid(cellId);
        sensors[sensorOffset + sen][gridCellId] = -1.0;
        sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      }
      if(!m_reconstructOldG && m_LsRotate && globalTimeStep > 0) {
        if(m_initialGCell[cellId] == 1) {
          const MInt gridCellId = grid().tree().solver2grid(cellId);
          sensors[sensorOffset + sen][gridCellId] = 1.0;
          sensorCellFlag[gridCellId][sensorOffset + sen] = true;
        }
      }
    } else {
      ASSERT(inList(cellId) > 0, "");
      const MInt gridCellId = grid().tree().solver2grid(cellId);
      if(a_level(cellId) < this->m_maxSensorRefinementLevel[sen]) { // refine cell
        if(c_noChildren(cellId) > 0) continue;
        if(a_isHalo(cellId)) continue;
        sensors[sensorOffset + sen][gridCellId] = 1.0;
        sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      }
    }
  }

  // find additional small-geometries:

  if(globalTimeStep < 0 && !m_combustion) {
    //&& m_adaptationLevel < this->m_maxSensorRefinementLevel[sen]
    //(minLevel() + 2)
    /*&& !m_reconstructBand*/
    //&& !m_combustion) {
    if(m_engineSetup && m_adaptationLevel < this->m_maxSensorRefinementLevel[sen]) {
      for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
        for(MInt set = 0; set < m_noSets; set++) {
          if(!m_computeSet_backup[set]) continue;
          if(a_inBandG(cellId, set)) {
            const MInt gridCellId = grid().tree().solver2grid(cellId);
            sensors[sensorOffset + sen][gridCellId] = 1;
            sensorCellFlag[gridCellId][sensorOffset + sen] = 1;
          }
        }
      }
    } else if(m_adaptationLevel == minLevel() && m_reconstructBand < 1) {
      for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
        for(MInt set = 0; set < m_noSets; set++) {
          if(a_inBandG(cellId, set)) {
            const MInt gridCellId = grid().tree().solver2grid(cellId);
            sensors[sensorOffset + sen][gridCellId] = 1;
            sensorCellFlag[gridCellId][sensorOffset + sen] = 1;
          }
        }
      }
    }
  }


  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
}

/**
 * \brief set the sensors for a single adaptation step
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::setSensors(std::vector<std::vector<MFloat>>& sensors,
                                         std::vector<MFloat>& sensorWeight,
                                         std::vector<std::bitset<64>>& sensorCellFlag,
                                         std::vector<MInt>& sensorSolverId) {
  TRACE();

  if(this->m_noSensors < 1) {
    return;
  }

  // If solver is inactive the sensor arrays still need to be set to optain the
  // correct offsets
  if(!isActive()) {
    const MInt sensorOffset = (signed)sensors.size();
    ASSERT(sensorOffset == 0 || grid().raw().treeb().noSolvers() > 1, "");
    sensors.resize(sensorOffset + this->m_noSensors, vector<MFloat>(grid().raw().m_noInternalCells, F0));
    sensorWeight.resize(sensorOffset + this->m_noSensors, -1);
    sensorCellFlag.resize(grid().raw().m_noInternalCells, sensorOffset + this->m_noSensors);
    sensorSolverId.resize(sensorOffset + this->m_noSensors, solverId());
    ASSERT(sensorOffset + this->m_noSensors < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

    for(MInt sen = 0; sen < this->m_noSensors; sen++) {
      sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
    }

    return;
  }

  // necessary here, as the sensors are also set on lower Grid-levels!
  // for avoid updateLowerGridLevels
  // updateLowerGridLevels();

  // the sensors are added at the end of the previous sensor-vectors!
  const auto sensorOffset = (signed)sensors.size();
  ASSERT(sensorOffset == 0 || grid().raw().treeb().noSolvers() > 1, "");
  sensors.resize(sensorOffset + this->m_noSensors, vector<MFloat>(grid().raw().m_noInternalCells, F0));
  sensorWeight.resize(sensorOffset + this->m_noSensors, -1);
  sensorCellFlag.resize(grid().raw().m_noInternalCells, sensorOffset + this->m_noSensors);
  sensorSolverId.resize(sensorOffset + this->m_noSensors, solverId());
  ASSERT(sensorOffset + this->m_noSensors < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

  // Sohels-sensor-method:
  if(m_combustion) {
    MIntScratchSpace regrid(1, AT_, "regrid");

    // FIXME labels:LS,toenhance
    // more efficient solution can be found!!:
    // levelSetAdaptationTrigger is also run in the application loop at the moment to figure out if the grid controller
    // adaptation has to be forced for the ls-solver it is also included here in case that a different solver wants to
    // adapt since in that case this function has to do certain things even if levelSetAdaptationTrigger() is false
    regrid.p[0] = levelSetAdaptationTrigger();
    regrid.p[0] = 1;
    // FIXME labels:LS
    // only using this for initial refinement right now
    // should be removed once the proper initialAdaptation routines are used instead of this one
    if(globalTimeStep == 0 || globalTimeStep == -1) {
      regrid.p[0] = 1;
    }

    // Use all sets for the mesh-adaptation!
    for(MInt set = 0; set < m_noSets; set++) {
      m_computeSet_tmp[set] = m_computeSet[set];
      m_computeSet[set] = true;
    }

    // FIXME labels:LS
    // to remove, its replaced by levelSetAdaptationTrigger
    // From this part on a level-Set based adaptation will be triggered
    m_forceAdaptation = true;

    if(regrid.p[0] == 0) {
      return;
    }

    MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
    MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");

    MInt tmpCount = 0;
    MIntScratchSpace tmp(m_maxNoCells, AT_, "tmp");
    MInt lastLayerCount = 0;
    MIntScratchSpace lastLayer(m_maxNoCells, AT_, "lastLayer");
    MBoolScratchSpace isAdded(m_maxNoCells, AT_, "isAdded");
    MInt listCount = 0;
    MIntScratchSpace list(m_maxNoCells, AT_, "list");
    MBoolScratchSpace coarseningFlag(grid().raw().treeb().size(), AT_, "coarseningFlag");

    for(MInt i = 0; i < a_noCells(); i++) {
      a_isBndryCellG(i) = false;
    }

    for(MInt c = 0; c < m_maxNoCells; c++) {
      isAdded.p[c] = false;
    }

    for(MInt c = 0; c < grid().raw().treeb().size(); c++) {
      coarseningFlag.p[c] = true;
    }

    // if this is not done the g0 cells are sometimes not determined correctly after adaptation
    determineG0Cells(0);
    // first layer consists of g0 cells
    MInt endSet = m_noSets;
    if(m_buildCollectedLevelSetFunction) endSet = 1;
    for(MInt set = 0; set < endSet; set++) {
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        MInt cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue; // dont add halo cells here
        lastLayer.p[lastLayerCount] = cellId;
        list.p[lastLayerCount] = cellId;
        isAdded.p[cellId] = true;
        lastLayerCount++;
        listCount++;
      }
    }

    ASSERT(lastLayerCount == listCount, "");

    // now add the correct halo cells to the list
    // gather:
    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      sendBufferSize.p[i] = 0;
      for(MInt j = 0; j < noWindowCells(i); j++) {
        m_intSendBuffers[i][sendBufferSize.p[i]++] = isAdded(windowCellId(i, j));
      }
    }
    if(grid().azimuthalPeriodicity()) {
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
          m_intSendBuffers[i][sendBufferSize.p[i]++] = isAdded(grid().azimuthalWindowCell(i, j));
        }
      }
    }
    // exchange data -> send, receive
    exchangeIntBuffers(sendBufferSize.getPointer(), receiveBufferSize.getPointer(), 4, 1);
    // scatter:
    // update the halo information
    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
      // check, if received data size matches expected size:
      if(!(receiveBufferSize.p[i] == noHaloCells(i))) {
        cerr << "this was not expected to happen: wrong number of halo information, buf=" << receiveBufferSize.p[i]
             << "noGHaloCells=" << noHaloCells(i) << ", has been added" << endl;
        ;
      }
#endif
      for(MInt j = 0; j < noHaloCells(i); j++) {
        MInt gc = haloCellId(i, j);
        if(m_intReceiveBuffers[i][j]) {
          lastLayer.p[lastLayerCount] = gc;
          lastLayerCount++;
          isAdded.p[gc] = true;
        }
      }
    }
    if(grid().azimuthalPeriodicity()) {
      for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
        MInt offset = noHaloCells(grid().azimuthalNeighborDomain(i));
        for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
          MInt n = offset + j;
          MInt haloCell = grid().azimuthalHaloCell(i, j);
          if(m_intReceiveBuffers[i][n]) {
            lastLayer.p[lastLayerCount] = haloCell;
            lastLayerCount++;
            isAdded.p[haloCell] = true;
          }
        }
      }
    }


    //////////////////////////////
    MInt level = minLevel();

    // calculate the number of layers of the current level
    MInt factor = IPOW2((a_maxGCellLevel(0) - level));
    MInt currentWidth = m_gShadowWidth / factor;
    if(currentWidth * factor < m_gShadowWidth) { // round up the currentWidth value
      currentWidth++;
    }
    if(currentWidth == 1) {
      cerr << "WARNING: you only have 1 layer of cells for your level set band on the coarsest level right now... you "
              "may want to rethink your level set grid resolution"
           << endl;
    }
    MInt currentLayer = 0;

    tmpCount = 0;
    // add g0 cells to the list
    for(MInt i = 0; i < lastLayerCount; i++) {
      if(!isAdded.p[lastLayer.p[i]]) {
        list.p[listCount] = lastLayer.p[i];
        listCount++;
        isAdded.p[lastLayer.p[i]] = true;
      }
    }
    currentLayer = 0;
    while(currentLayer < currentWidth) {
      // find next layer of cells
      for(MInt i = 0; i < lastLayerCount; i++) {
        MInt gc = lastLayer.p[i];

        // if gc is on a higher level than the current level, take the parent until you find a parent with the correct
        // level
        while(a_level(gc) > level) {
          gc = c_parentId(gc);
        }

        for(MInt d = 0; d < m_noDirs; d++) {
          if(a_hasNeighbor(gc, d)) {
            if(!isAdded.p[c_neighborId(gc, d)]) {
              tmp.p[tmpCount] = c_neighborId(gc, d);
              isAdded.p[tmp.p[tmpCount]] = true;
              tmpCount++;
            }
          }
        }
      }

      // check if halo cells must be coarsened
      // gather:
      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
        sendBufferSize.p[i] = 0;
        for(MInt j = 0; j < noWindowCells(i); j++) {
          m_intSendBuffers[i][sendBufferSize.p[i]++] = isAdded(windowCellId(i, j));
        }
      }
      if(grid().azimuthalPeriodicity()) {
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
            m_intSendBuffers[i][sendBufferSize.p[i]++] = isAdded[grid().azimuthalWindowCell(i, j)];
          }
        }
      }
      // exchange data -> send, receive
      exchangeIntBuffers(sendBufferSize.getPointer(), receiveBufferSize.getPointer(), 4, 1);
      // scatter:
      // update the halo information
      for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
#ifdef LS_DEBUG
        // check, if received data size matches expected size:
        if(!(receiveBufferSize.p[i] == noHaloCells(i))) {
          m_log << "this was not expected to happen: wrong number of halo information, buf=" << receiveBufferSize.p[i]
                << "noGHaloCells=" << noHaloCells(i) << ", has been added" << endl;
          ;
        }
#endif
        for(MInt j = 0; j < noHaloCells(i); j++) {
          MInt gc = haloCellId(i, j);
          if(m_intReceiveBuffers[i][j]) {
            tmp.p[tmpCount] = gc;
            isAdded.p[tmp.p[tmpCount]] = true;
            tmpCount++;
          }
        }
      }
      if(grid().azimuthalPeriodicity()) {
        for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
          MInt offset = noHaloCells(grid().azimuthalNeighborDomain(i));
          for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
            MInt n = offset + j;
            MInt haloCell = grid().azimuthalHaloCell(i, j);
            if(m_intReceiveBuffers[i][n]) {
              tmp.p[tmpCount] = haloCell;
              isAdded.p[tmp.p[tmpCount]] = true;
              tmpCount++;
            }
          }
        }
      }

      // add new layer to list
      for(MInt i = 0; i < tmpCount; i++) {
        list.p[listCount] = tmp.p[i];
        listCount++;
        // and replace the last layer list with the tmp list
        lastLayer.p[i] = tmp.p[i];
      }
      lastLayerCount = tmpCount;
      tmpCount = 0;
      currentLayer++;
    }

    // add all children to the list
    for(MInt i = 0; i < listCount; i++) {
      MInt gc = list.p[i];
      if(c_noChildren(gc) > 0) {
        for(MInt child = 0; child < IPOW2(nDim); child++) {
          if(c_childId(gc, child) == -1) {
            continue;
          }
          MInt cId = c_childId(gc, child);
          if(!isAdded.p[cId]) {
            isAdded.p[cId] = true;
            list.p[listCount] = cId;
            listCount++;
          }
        }
      }
    }

    // set sensor for all cells that have been added to the list
    for(MInt c = 0; c < listCount; c++) {
      MInt gridCellId = grid().tree().solver2grid(list.p[c]);
      coarseningFlag.p[gridCellId] = false;
      if(a_level(list.p[c]) < a_maxGCellLevel(0)) { // refine cell
        if(c_noChildren(list.p[c]) > 0) continue;
        if(a_isHalo(list.p[c])) continue;
        sensors[sensorOffset][gridCellId] = 1.0;
        sensorCellFlag[gridCellId][sensorOffset] = 1;
      }
    }

    // now coarsen all cells that are not tagged and that have no children
    // only do this for the grid proxy though
    for(MInt i = 0; i < grid().tree().size(); i++) {
      if(!coarseningFlag.p[grid().tree().solver2grid(i)]) {
        continue;
      }
      if(grid().tree().level(i) == minLevel()) {
        continue;
      }
      if(grid().tree().noChildren(i) == 0) {
        MInt solverCellId = i;
        if(grid().tree().hasProperty(solverCellId, Cell::IsHalo)) {
          continue;
        }
        MInt gridCellId = grid().tree().solver2grid(solverCellId); // use correct grid id for the sensor
        sensors[sensorOffset][gridCellId] = -1.0;
        sensorCellFlag[gridCellId][sensorOffset] = 1;
      }
    }

  } else { // using Tims-sensor-method:

    // only set sensors if the adaptation was globally triggered in buildLevelSetTubeCG!
    if(this->m_adapts && ((m_forceAdaptation && globalTimeStep > 0) || (globalTimeStep < 0))) {
      if(domainId() == 0) {
        cerr << "Setting " << this->m_noSensors << " sensor(s) for the ls-solver adaptation!" << endl;
      }

      for(MInt sen = 0; sen < this->m_noSensors; sen++) {
        (this->*(this->m_sensorFnPtr[sen]))(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen);
      }
    }

    /*
    // debug-output:
    for(MInt i=0;i<grid().raw().treeb().size();i++){
      grid().raw().treeb().weight(i)=1;
    }
    m_bodyIdOutput = true;
    MInt globalTimeStep_temp = globalTimeStep;
    globalTimeStep = m_adaptationLevel;
    // for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    //  const MInt gridCellId = grid().tree().solver2grid(cellId);
    //  a_bodyIdG(cellId, 0) = sensors[sensorOffset][gridCellId];
    //}
    writeRestartLevelSetFileCG(1, "restartLSGridCG_sensors" , "restartLSCG_sensors");
    globalTimeStep = globalTimeStep_temp;
    */

    ASSERT(m_freeIndices.empty(), "");
    m_freeIndices.clear();

    m_refinedCells.clear();
  }
}

/**
 * \brief  reinit the solver after a single adaptation step
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::postAdaptation() {
  TRACE();

  this->compactCells();
  // URL
  if(!g_multiSolverGrid) {
    for(MInt gridCellId = 0; gridCellId < grid().raw().treeb().size(); gridCellId++) {
      ASSERT(grid().tree().solver2grid(gridCellId) == gridCellId, "");
      ASSERT(grid().tree().grid2solver(gridCellId) == gridCellId, "");
    }
  }
  //}

  grid().updateOther();

  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);

  // Nothing further to be done if solver inactive
  if(!isActive()) return;

  grid().updateLeafCellExchange();

  ASSERT(!m_constructGField, "");

  generateListOfGExchangeCellsCG();

  this->checkNoHaloLayers();

  // Initialize the azimuthal periodic exchange
  initAzimuthalExchange();

  exchangeAllLevelSetData();
  exchangeGapInfo();

  // previous part from initialRefinement!
  if(globalTimeStep < 0 && m_geometryChange == nullptr) {
    if(m_combustion) {
      initSolver();
    } else {
      for(MInt set = 0; set < m_noSets; set++) {
        std::vector<MInt>().swap(m_bandCells[set]);
      }

      if(m_GFieldInitFromSTL) constructGFieldFromSTL(1);

      if(m_buildCollectedLevelSetFunction) buildCollectedLevelSet(0);

      determineG0Cells();
      createGgridCG();

      if(m_GFieldInitFromSTL) {
        determineG0Cells();
        determineBandCells();
        m_log << "Initialize G Field from stl data...";
        constructGFieldFromSTL(3);
        m_log << "ok" << endl;
        for(MInt set = 0; set < m_noSets; set++) {
          std::vector<MInt>().swap(m_bandCells[set]);
        }
      }

      if(m_buildCollectedLevelSetFunction) buildCollectedLevelSet(0);
      determineG0Cells();
      determineBandCells();
      updateBndryCellList();
      resetOutsideCells();
      if(m_GFieldInitFromSTL) constructGFieldFromSTL(0);
      if(m_buildCollectedLevelSetFunction) buildCollectedLevelSet(0);
      if(m_GFieldInitFromSTL) resetOutsideCells();

      if(m_LsRotate) resetContainingGCells();
      exchangeAllLevelSetData();
      checkHaloCells();
    }
  } else if(globalTimeStep < 0 && m_geometryChange != nullptr) {
    for(MInt set = 0; set < m_noSets; set++) {
      std::vector<MInt>().swap(m_bandCells[set]);
    }

    determineG0Cells();
    determineBandCells();
    m_log << "Initialize G Field from mew stl data...";
    constructGFieldFromSTL(5);
    m_log << "ok" << endl;

    // NOTE: two loops ensure that new bandCells are also already initialised in the
    //      first adaptationLoop. This increases the load but ensures the current
    //      construction of the new levelSet Field!
    //      thus only necessary if this is only called once!
    if(maxRefinementLevel() - maxUniformRefinementLevel() < 2) {
      for(MInt set = 0; set < m_noSets; set++) {
        std::vector<MInt>().swap(m_bandCells[set]);
      }

      determineG0Cells();
      determineBandCells();
      updateBndryCellList();
      resetOutsideCells();

      determineG0Cells();
      determineBandCells();
      m_log << "Initialize G Field from mew stl data...";
      constructGFieldFromSTL(5);
      m_log << "ok" << endl;
    }

    if(m_buildCollectedLevelSetFunction) buildCollectedLevelSet(2);

    reInitSolver(true);

    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      for(MInt set = m_startSet; set < m_noSets; set++) {
        if(!m_geometryChange[set]) continue;
        a_oldLevelSetFunctionG(cellId, set) = a_levelSetFunctionG(cellId, set);
      }
    }

    if(m_LsRotate) resetContainingGCells();

    // exchange required afterwards
    exchangeDataLS(&(a_oldLevelSetFunctionG(0, 0)), m_maxNoSets);

  } else { // globalTimeStep > -1
    if(m_combustion) {
      if(globalTimeStep == 0) {
        initSolver();
      }
    } else {
      MInt noRefinedBandCells = (signed)m_refinedCells.size();
      MPI_Allreduce(MPI_IN_PLACE, &noRefinedBandCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                    "noRefinedBandCells");

      // On the fly reconstruction of levelset band cells outside the G0 Set!
      if(noRefinedBandCells > 0 && !m_virtualSurgery) {
        if(domainId() == 0 && m_reconstructBand > 1) {
          cerr << "Reinitialisation of " << noRefinedBandCells << " cells " << endl;
        }
        if(domainId() == 0 && m_reconstructBand == 1) {
          cerr << "Interpolation of " << noRefinedBandCells << " cells " << endl;
        }

        m_newRefinedBand = true;

        ASSERT(m_semiLagrange, "");
        ASSERT(m_reconstructBand > 0,
               "CAUTION: refining a band cell after the initial adaptation without reconstructBand!");

        ASSERT(m_GFieldInitFromSTL, "");

        // construct old-LevelSet data based on shifted geometry!
        if(m_reconstructBand > 1) {
          constructGFieldFromSTL(4);
          exchangeDataLS(&(a_oldLevelSetFunctionG(0, 0)), m_maxNoSets);

          if(m_maxLevelChange) {
            resetOldOutsideCells();
          }

          // copy old levelset to levelset for stationary bodies
          for(auto& m_refinedCell : m_refinedCells) {
            const MInt cellId = m_refinedCell.first;
            ASSERT(!a_isHalo(cellId), "");
            const MInt set = m_refinedCell.second;
            if(set > 0) {
              if(!m_computeSet_backup[set]) {
                a_levelSetFunctionG(cellId, set) = a_oldLevelSetFunctionG(cellId, set);
              }
            } else {
              for(MInt setI = m_startSet; setI < m_noSets; setI++) {
                if(m_computeSet_backup[setI]) {
                  continue;
                }
                a_levelSetFunctionG(cellId, setI) = a_oldLevelSetFunctionG(cellId, setI);
              }
            }
          }
        } else {
          // interpolate level set
          for(auto& m_refinedCell : m_refinedCells) {
            const MInt cellId = m_refinedCell.first;
            const MInt set = m_refinedCell.second;
            if(set == 0) {
              continue;
            }
            std::array<MInt, 8> interpolationCells = {0, 0, 0, 0, 0, 0, 0, 0};
            std::array<MFloat, nDim> cellPos = {};
            for(MInt i = 0; i < nDim; i++) {
              cellPos[i] = c_coordinate(cellId, i);
            }
            const MInt parent = c_parentId(cellId);
            const MInt position = setUpLevelSetInterpolationStencil(parent, interpolationCells.data(), cellPos.data());
            if(position > -1) {
              a_levelSetFunctionG(cellId, set) = interpolateLevelSet(interpolationCells.data(), cellPos.data(), set);
              a_oldLevelSetFunctionG(cellId, set) =
                  interpolateOldLevelSet(interpolationCells.data(), cellPos.data(), set);
            }
          }
        }
        exchangeLevelSet();

        reInitSolver(m_forceAdaptation);

        if(grid().azimuthalPeriodicity()) {
          this->exchangeAzimuthalPer(&(a_oldLevelSetFunctionG(0, 0)), m_maxNoSets);
        }
      }
    }
  }

  // debug-output:
  /*
  for(MInt i=0;i<grid().raw().treeb().size();i++){
    grid().raw().treeb().weight(i)=1;
  }
  MInt globalTimeStep_temp = globalTimeStep;
  globalTimeStep = m_adaptationLevel;
  writeRestartLevelSetFileCG(1, "restartLSGridCG_init" , "restartLSCG_init");
  globalTimeStep = globalTimeStep_temp;
  */
  m_adaptationLevel++;
}


/**
 * \brief  reinit the solver after the full adaptation loop!
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::finalizeAdaptation() {
  TRACE();

  // Set properties and return if solver is inactive
  if(!isActive()) {
    m_forceAdaptation = false;
    m_adaptationSinceLastRestart = false;
    return;
  }

  if(m_combustion || m_freeSurface) {
    finalizeLevelSet_(0, 0);
    m_forceAdaptation = false;
    m_adaptationSinceLastRestart = true;
    return;
  }

  if(globalTimeStep > 0) { // part from reinitAfterAdaptation:


    reInitSolver(m_forceAdaptation);

    if(m_LsRotate) {
      updateContainingGCells();
      copyWindowToHaloIds();
    }

    if(m_maxLevelChange) {
      if(domainId() == 0) {
        cerr << "Ls-MaxLevel after adaptation: " << maxLevel() << endl;
      }
    }

  } else if(m_geometryChange == nullptr) { // part from postInitialRefinementStep:

    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      for(MInt set = 0; set < m_maxNoSets; set++) {
        if(!m_semiLagrange) {
          for(MInt d = 0; d < nDim; d++) {
            a_extensionVelocityG(cellId, d, set) = F0;
          }
        } else {
          a_oldLevelSetFunctionG(cellId, set) = a_levelSetFunctionG(cellId, set);
        }
      }
      if(m_maxNoSets > 1 && m_closeGaps) {
        a_potentialGapCell(cellId) = 0;
        a_potentialGapCellClose(cellId) = 0;
      }
    }

    setGCellBndryProperty();

    if(!m_semiLagrange) {
      computeNormalVectors();
      computeCurvature();
    }

    finalizeLevelSetInitialization();

    buildMultipleLevelSet(2);

    if(m_GFieldFromSTLInitCheck) {
      // FIXME labels:LS hack to initialise weights,
      // necessary for writing a restart-file out of the solver, after an adaptation
      // move into the meshAdaptation to initialice new cells with weight 1/0
      for(MInt i = 0; i < grid().raw().treeb().size(); i++) {
        grid().raw().treeb().weight(i) = 1;
      }

      MInt globalTimeStep_temp = globalTimeStep;
      globalTimeStep = 0;
      if(!m_initFromRestartFile) {
        writeRestartLevelSetFileCG(1, "restartLSGridCG_init", "restartLSCG_init");
      } else {
        writeRestartLevelSetFileCG(1, "restartLSGridCG_reinit", "restartLSCG_reinit");
      }
      globalTimeStep = globalTimeStep_temp;
    }

    if(m_LSSolver) initSolver();

  } else { // globalTimeStep = -1 && m_geometryChange != nullptr

    finalizeLevelSetInitialization();

    if(!m_oldG0Cells.empty()) {
      MInt startSet = 0;
      if(m_buildCollectedLevelSetFunction) startSet = 1;

      // check new G0-cells and delete entries which were were also a G0 cell in the old geometry!
      for(MInt set = startSet; set < m_noSets; set++) {
        if(!m_geometryChange[set]) continue;
        for(MInt id = 0; id < a_noG0Cells(set); id++) {
          const MInt cellId = a_G0CellId(id, set);
          if(a_isHalo(cellId)) continue;
          auto it0 = m_oldG0Cells.find(cellId);
          if(it0 != m_oldG0Cells.end()) {
            const MInt oldSet = it0->second;
            if(oldSet == set) {
              m_oldG0Cells.erase(it0);
            }
          }
        }
      }
    }
  }

  if((m_maxLevelChange || m_geometryChange != nullptr) && m_trackMovingBndry) {
    // shift STL back!
    ASSERT(m_GCtrlPntMethod == 2, "Only working for controlPoint-Method 2!");
    MFloat dx[3] = {F0, F0, F0};
    MFloat xCurrent[3] = {F0, F0, F0};
    MFloat xOld[3] = {F0, F0, F0};
    const MFloat usedTime = m_geometryChange != nullptr ? time() : time() + timeStep();
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet_backup[set]) continue;
      for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
        const MInt body = m_setToBodiesTable[set][b];
        const MInt bcId = m_bodyBndryCndIds[body];
        computeBodyPropertiesForced(1, xCurrent, body, usedTime);
        computeBodyPropertiesForced(1, xOld, body, 0.0);

        for(MInt d = 0; d < nDim; d++) {
          // dx[d] = -(xCurrent[d] - xOld[d]);
          dx[d] = -m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body];
        }
        m_gCtrlPnt.CtrlPnt2_shiftSTL(bcId, dx);
      }
    }
  }

  // reset properties back to default after the sucessful adaptation
  if(m_maxLevelChange && (maxLevel() == maxRefinementLevel() || maxLevel() == this->m_maxSensorRefinementLevel[0])) {
    m_maxLevelChange = false;
  }
  if(m_geometryChange != nullptr) {
    mDeallocate(m_geometryChange);
  }

  exchangeAllLevelSetData();

  if(m_virtualSurgery) {
    if(globalTimeStep < 0) {
      initializeCollectedLevelSet(0);
    } else {
      initializeCollectedLevelSet(1);
    }
  }

  // ASSERT sanity checks
  testCellsCG();

  checkHaloCells();

  m_forceAdaptation = false;
  m_adaptationSinceLastRestart = true;
}

/**
 * \brief  prepare solver adaptation
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::prepareAdaptation() {
  TRACE();

  ASSERT(m_freeIndices.empty(), "");

  m_adaptationLevel = maxUniformRefinementLevel();
  m_newRefinedBand = false;
  m_oldG0Cells.clear();

  if(!isActive()) return;

  if(this->m_noSensors <= 0) {
    return;
  }

  if(maxRefinementLevel() > maxLevel() && globalTimeStep > 0) {
    m_forceAdaptation = true;
    if(domainId() == 0) {
      cerr << "Max Refinement Level has been increased at restart. " << endl;
      cerr << "Initialising interpolation at adaptation!" << endl;
    }

    // recalculate outside-G-value for max-Level-change!
    if(m_maxLevelChange) {
      if(!Context::propertyExists("maxGCellLevel", m_solverId)) {
        m_maxGCellLevel[0] = Context::getSolverProperty<MInt>("maxRfnmntLvl", m_solverId, AT_);
      } else {
        m_maxGCellLevel[0] = Context::getSolverProperty<MInt>("maxGCellLevel", m_solverId, AT_, 0);
      }

      m_gCellDistance = c_cellLengthAtLevel(m_maxGCellLevel[0]);
      m_FgCellDistance = F1 / m_gCellDistance;
      m_log << "smallest gCell length " << m_gCellDistance << endl;

      m_outsideGValue = (MFloat)(2 * m_gBandWidth * m_gCellDistance);
      m_log << "Outside G value: " << m_outsideGValue << endl;
    }
  }

  if(this->m_maxSensorRefinementLevel[0] < maxLevel()) {
    m_forceAdaptation = true;
    if(domainId() == 0) {
      cerr << "Max Refinement Level has been decreased at restart." << endl;
      cerr << "Initialising interpolation at adaptation!" << endl;
    }
    if(Context::propertyLength("maxRfnmntLvl", m_solverId) == 1) {
      m_maxGCellLevel[0] = this->m_maxSensorRefinementLevel[0];
    } else {
      for(MInt set = 0; set < m_maxNoSets; set++) {
        m_maxGCellLevel[set] = this->m_maxSensorRefinementLevel[0];
      }
    }

    m_gCellDistance = c_cellLengthAtLevel(m_maxGCellLevel[0]);
    m_FgCellDistance = F1 / m_gCellDistance;
    m_log << "smallest gCell length " << m_gCellDistance << endl;

    m_outsideGValue = (MFloat)(2 * m_gBandWidth * m_gCellDistance);
    m_log << "Outside G value: " << m_outsideGValue << endl;
  }

  if(grid().newMinLevel() > maxUniformRefinementLevel()) m_forceAdaptation = true;


  if(globalTimeStep < 0 && m_geometryChange != nullptr) {
    // save old G0-Cells, so that a sponging can be applied to the fv-solution
    // if they are no longer a G0-cell for the new geometry!

    MInt startSet = 0;
    if(m_buildCollectedLevelSetFunction) startSet = 1;

    for(MInt set = startSet; set < m_noSets; set++) {
      if(!m_geometryChange[set]) continue;
      for(MInt id = 0; id < a_noG0Cells(set); id++) {
        const MInt cellId = a_G0CellId(id, set);
        if(a_isHalo(cellId)) continue;
        ASSERT(a_isGZeroCell(cellId, set), "");
        m_oldG0Cells.insert(make_pair(cellId, set));
      }
    }
  }

  // the levelset of a moving STL will be reconstructed during the adaptation!
  // shift all STLs to the correct place
  if((m_geometryChange != nullptr || m_maxLevelChange) && m_trackMovingBndry) {
    ASSERT(m_GCtrlPntMethod == 2, "Only working for controlPoint-Method 2!");
    MFloat dx[3] = {F0, F0, F0};
    MFloat xCurrent[3] = {F0, F0, F0};
    MFloat xOld[3] = {F0, F0, F0};
    const MFloat usedTime = m_geometryChange != nullptr ? time() : time() + timeStep();
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_computeSet_backup[set]) continue;
      for(MInt b = 0; b < m_noBodiesInSet[set]; b++) {
        const MInt body = m_setToBodiesTable[set][b];
        const MInt bcId = m_bodyBndryCndIds[body];
        computeBodyPropertiesForced(1, xCurrent, body, usedTime);
        computeBodyPropertiesForced(1, xOld, body, 0.0);

        for(MInt d = 0; d < nDim; d++) {
          // dx[d] = xCurrent[d] - xOld[d];
          dx[d] = m_semiLagrange_xShift_ref[d * m_noEmbeddedBodies + body];
        }
        m_gCtrlPnt.CtrlPnt2_shiftSTL(bcId, dx);
      }
    }
  }
}

/**
 * \brief Wrapper for finalizeLevelSet()
 * \author Thomas Hoesgen
 * \date 11/2019
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::postTimeStep() {
  TRACE();

  if(!m_trackMovingBndry || globalTimeStep < m_trackMbStart || globalTimeStep > m_trackMbEnd) {
    return;
  }

  RECORD_TIMER_START(m_timers[Timers::PostTime]);

  if(m_combustion) {
    m_forceAdaptation = levelSetAdaptationTrigger();
    if(!m_forceAdaptation) {
      finalizeLevelSet_(0, 0);
    }

  } else if(m_freeSurface) {
    finalizeLevelSet_(0, 0);

  } else {
    finalizeLevelSet();
  }

  RECORD_TIMER_STOP(m_timers[Timers::PostTime]);
}


/**
 * \brief  reinit levelset solver after adaptation or balance
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::reInitSolver(const MBool regrid) {
  TRACE();

  // reset
  for(MInt set = 0; set < m_noSets; set++) {
    std::vector<MInt>().swap(m_bandCells[set]);
    m_computeSet[set] = true;
    m_changedSet[set] = true;
  }

  // recreate a_regridTriggerG
  if(regrid) {
    determineG0Cells();
    if(m_buildCollectedLevelSetFunction) determineG0Cells(0);

    createGgridCG(); // requires all G0-cells
  } else {
    determineG0Cells();
  }

  determineBandCells();

  if(!m_semiLagrange) {
    setGCellBndryProperty();
    updateBndryCellList();
  }

  resetOutsideCells();

  setUpPotentialGapCells();
  buildMultipleLevelSet(2);
  exchangeGapInfo();

  updateLowerGridLevels();

  exchangeAllLevelSetData();

  checkHaloCells();

  // reset m_computeSet, to the original values before the adaptation
  if(globalTimeStep > -1) {
    for(MInt set = 0; set < m_noSets; set++) {
      m_computeSet[set] = m_computeSet_backup[set];
    }
  }
}

/** \brief help-function for engine calculations which returns the crank-angle
 *         for a given time
 *         mode = 0: return CAD in range of (0-720)
 *         mode = 1: return accumulated crankAnge in radian
 *
 *  \author Tim Wegmann
 */
template <MInt nDim>
MFloat LsCartesianSolver<nDim>::crankAngle(const MFloat elapsedTime, const MInt mode) {
  TRACE();

  MFloat& Strouhal = m_static_crankAngle_Strouhal;
  MFloat& initialCad = m_static_crankAngle_initialCad;

  if(initialCad < 0) {
    Strouhal = Context::getSolverProperty<MFloat>("Strouhal", m_solverId, AT_, &Strouhal);
    initialCad = 0;
    initialCad = Context::getSolverProperty<MFloat>("initialCrankAngle", m_solverId, AT_, &initialCad);
  }

  return maia::math::crankAngle(elapsedTime, Strouhal, initialCad, mode);
}

/** \brief
 *  \Init azimuthal ls exchange
 *  \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::initAzimuthalExchange() {
  TRACE();

  if(grid().noAzimuthalNeighborDomains() == 0) return;

  MUint sndSize = maia::mpi::getBufferSize(grid().azimuthalHaloCells());
  MFloatScratchSpace haloBuff(sndSize * nDim, AT_, "haloBuff");
  haloBuff.fill(0);
  MUint rcvSize = maia::mpi::getBufferSize(grid().azimuthalWindowCells());
  MFloatScratchSpace windowBuff(rcvSize * nDim, AT_, "windowBuff");
  windowBuff.fill(0);

  this->m_azimuthalCartRecCoord.clear();
  this->m_azimuthalCartRecCoord.reserve(rcvSize * nDim);

  // Get azimuthal image coordinate
  MFloat angle = grid().azimuthalAngle();
  MInt sndCnt = 0;
  MFloat haloCoords[nDim];
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      MInt haloId = grid().azimuthalHaloCell(i, j);
      ASSERT(grid().raw().a_hasProperty(grid().tree().solver2grid(haloId), Cell::IsPeriodic),
             "azimuthal halo is not isPeriodic!");

      for(MInt d = 0; d < nDim; d++) {
        haloCoords[d] = c_coordinate(haloId, d);
      }

      MInt side = grid().determineAzimuthalBoundarySide(haloCoords);

      grid().rotateCartesianCoordinates(haloCoords, (side * angle));

      for(MInt d = 0; d < nDim; d++) {
        haloBuff[sndCnt++] = haloCoords[d];
      }
    }
  }

  maia::mpi::exchangeBuffer(grid().azimuthalNeighborDomains(), grid().azimuthalWindowCells(),
                            grid().azimuthalHaloCells(), mpiComm(), haloBuff.getPointer(), windowBuff.getPointer(),
                            nDim);

  MInt rcvCnt = 0;
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
      for(MInt d = 0; d < nDim; d++) {
        this->m_azimuthalCartRecCoord[rcvCnt] = windowBuff[rcvCnt];
        rcvCnt++;
      }
    }
  }
}

/** \brief rotates stl for rotating ls
 *  \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::rotateSTL(MInt direction, MInt body, MFloat* center) {
  TRACE();

  MFloat phiRot[nDim];
  const MInt bcId = m_bodyBndryCndIds[body];

  for(MInt n = 0; n < nDim; n++) {
    phiRot[n] = -1 * direction * m_semiLagrange_xRot_STL[body * nDim + n];
  }

  m_gCtrlPnt.CtrlPnt2_rotateSTL(bcId, phiRot, center);
}

/** \brief rotates stl for rotating ls
 *  \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::rotateSTL(MInt direction) {
  TRACE();

  MFloat center[nDim];

  for(MInt b = 0; b < m_noBodiesToCompute; b++) {
    const MInt body = m_bodiesToCompute[b];

    computeBodyPropertiesForced(1, center, body, time());

    rotateSTL(direction, body, center);
  }
}

/** \brief returns the globalIds of the containingCells for rotating ls with reconstructOldG=false
 *  \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::localToGlobalIdsContainingCells() {
  TRACE();

  for(MInt b = 0; b < m_noBodiesToCompute; b++) {
    MInt body = m_bodiesToCompute[b];
    MInt set = m_bodyToSetTable[body];

    MIntScratchSpace noCellsToDom(grid().noDomains(), AT_, "noCellsToDom");
    noCellsToDom.fill(0);

    for(MInt id = 0; id < a_noBandCells(set); id++) {
      MInt cellId = a_bandCellId(id, set);
      if(!(a_bodyIdG(cellId, set) == body)) continue;
      if(grid().tree().hasProperty(cellId, Cell::IsHalo)) continue;
      if(a_level(cellId) != a_maxGCellLevel()) continue;

      MInt domId = a_containingDomain(cellId, b);
      if(domId > -1) {
        noCellsToDom[domId]++;
      }
    }

    startLoadTimer(AT_);

    prepareGlobalComm(&noCellsToDom[0]);

    MInt noCellsComm = mMax(m_globalSndOffsets[grid().noDomains()], m_globalRcvOffsets[grid().noDomains()]);

    MLongScratchSpace sndBufGlob(noCellsComm, AT_, "sndBufGlob");
    MLongScratchSpace rcvBufGlob(noCellsComm, AT_, "rcvBufGlob");
    MIntScratchSpace sndBufSizeGlob(grid().noDomains(), AT_, "sndBufSizeGlob");
    MIntScratchSpace rcvBufSizeGlob(grid().noDomains(), AT_, "rcvBufSizeGlob");

    sndBufGlob.fill(-1);
    rcvBufGlob.fill(-1);
    sndBufSizeGlob.fill(0);
    rcvBufSizeGlob.fill(0);

    std::vector<std::vector<MInt>> cellIdsLoc;
    cellIdsLoc.resize(grid().noDomains());

    for(MInt cell = 0; cell < a_noCells(); cell++) {
      if(a_isHalo(cell)) continue;
      ASSERT(a_domainId(c_globalId(cell)) == domainId(), "globalId conversion is broken!");
      MInt cellId = a_containingCell(cell, b);
      if(cellId > -1 && a_inBandG(cell, set) && a_bodyIdG(cell, set) == body && a_level(cell) == a_maxGCellLevel()) {
        MInt domId = a_containingDomain(cell, b);
        if(domainId() == domId) {
          ASSERT(a_localId(c_globalId(cellId)) == cellId, "GlobalId conversion is broken!");
          a_containingCell(cell, b) = c_globalId(cellId);
        } else {
          cellIdsLoc[domId].push_back(cell);
          sndBufGlob[m_globalSndOffsets[domId] + sndBufSizeGlob.p[domId]] = cellId;
          sndBufSizeGlob.p[domId]++;
        }
      } else {
        a_containingCell(cell, b) = -1;
      }
    }

    exchangeBuffersGlobal(sndBufGlob.getPointer(), rcvBufGlob.getPointer(), sndBufSizeGlob.getPointer(),
                          rcvBufSizeGlob.getPointer(), &m_globalSndOffsets[0], &m_globalRcvOffsets[0], 3);

    sndBufGlob.fill(-1);
    sndBufSizeGlob.fill(0);
    for(MInt i = 0; i < grid().noDomains(); i++) {
      MInt ind = m_globalRcvOffsets[i];
      for(MInt j = 0; j < rcvBufSizeGlob(i); j++) {
        MInt cellId = rcvBufGlob(ind + j);
        ASSERT(a_localId(c_globalId(cellId)) == cellId, "GlobalId conversion is broken!");
        sndBufGlob(ind + sndBufSizeGlob.p[i]) = c_globalId(cellId);
        sndBufSizeGlob.p[i]++;
      }
    }
    rcvBufGlob.fill(-1);
    rcvBufSizeGlob.fill(0);

    exchangeBuffersGlobal(sndBufGlob.getPointer(), rcvBufGlob.getPointer(), sndBufSizeGlob.getPointer(),
                          rcvBufSizeGlob.getPointer(), &m_globalRcvOffsets[0], &m_globalSndOffsets[0], 3);

    for(MInt i = 0; i < grid().noDomains(); i++) {
      MInt ind = m_globalSndOffsets[i];
      for(MInt j = 0; j < rcvBufSizeGlob(i); j++) {
        MInt cellId = cellIdsLoc[i][j];
        a_containingCell(cellId, b) = rcvBufGlob(ind + j);
      }
    }
  }
}

/** \brief stores the globalIds of the containingCells for rotating ls with reconstructOldG=false
 *  \author Thomas Hoesgen
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::globalToLocalIdsContainingCells() {
  TRACE();

  MIntScratchSpace noCellsToDom(grid().noDomains(), AT_, "noCellsToDom");
  noCellsToDom.fill(0);

  for(MInt b = 0; b < m_noBodiesToCompute; b++) {
    for(MInt i = 0; i < a_noCells(); ++i) {
      if(grid().tree().hasProperty(i, Cell::IsHalo)) continue;
      if(a_level(i) != a_maxGCellLevel()) continue;

      MLong globalId = a_containingCell(i, b);
      if(globalId > -1) {
        MInt domId = a_domainId(globalId);
        noCellsToDom[domId]++;
      }
    }

    prepareGlobalComm(&noCellsToDom[0]);

    MInt noCellsComm = mMax(m_globalSndOffsets[grid().noDomains()], m_globalRcvOffsets[grid().noDomains()]);

    MLongScratchSpace sndBufGlob(noCellsComm, AT_, "sndBufGlob");
    MLongScratchSpace rcvBufGlob(noCellsComm, AT_, "rcvBufGlob");
    MIntScratchSpace sndBufSizeGlob(grid().noDomains(), AT_, "sndBufSizeGlob");
    MIntScratchSpace rcvBufSizeGlob(grid().noDomains(), AT_, "rcvBufSizeGlob");
    sndBufGlob.fill(-1);
    rcvBufGlob.fill(-1);
    sndBufSizeGlob.fill(0);
    rcvBufSizeGlob.fill(0);

    std::vector<std::vector<MInt>> cellIdsLoc;
    cellIdsLoc.resize(grid().noDomains());

    for(MInt i = 0; i < a_noCells(); ++i) {
      if(grid().tree().hasProperty(i, Cell::IsHalo)) continue;

      ASSERT(a_domainId(c_globalId(i)) == domainId(), "globalId conversion is broken!");

      MLong globalId = a_containingCell(i, b);
      if(globalId > -1 && a_level(i) == a_maxGCellLevel()) {
        MInt domId = a_domainId(globalId);
        if(domainId() == domId) {
          ASSERT(c_globalId(a_localId(globalId)) == globalId, "GlobalId conversion is broken!");
          a_containingCell(i, b) = a_localId(globalId);
          a_containingDomain(i, b) = domId;
        } else {
          cellIdsLoc[domId].push_back(i);
          sndBufGlob[m_globalSndOffsets[domId] + sndBufSizeGlob.p[domId]] = globalId;
          sndBufSizeGlob.p[domId]++;
        }
      } else {
        a_containingCell(i, b) = -1;
        a_containingDomain(i, b) = -1;
      }
    }

    exchangeBuffersGlobal(sndBufGlob.getPointer(), rcvBufGlob.getPointer(), sndBufSizeGlob.getPointer(),
                          rcvBufSizeGlob.getPointer(), &m_globalSndOffsets[0], &m_globalRcvOffsets[0], 3);
    sndBufGlob.fill(-1);
    sndBufSizeGlob.fill(0);

    for(MInt i = 0; i < grid().noDomains(); i++) {
      MInt ind = m_globalRcvOffsets[i];
      for(MInt j = 0; j < rcvBufSizeGlob(i); j++) {
        MLong globalId = rcvBufGlob(ind + j);
        ASSERT(c_globalId(a_localId(globalId)) == globalId, "GlobalId conversion is broken!");
        sndBufGlob(ind + sndBufSizeGlob.p[i]) = a_localId(globalId);
        sndBufSizeGlob.p[i]++;
      }
    }
    rcvBufGlob.fill(-1);
    rcvBufSizeGlob.fill(0);

    exchangeBuffersGlobal(sndBufGlob.getPointer(), rcvBufGlob.getPointer(), sndBufSizeGlob.getPointer(),
                          rcvBufSizeGlob.getPointer(), &m_globalRcvOffsets[0], &m_globalSndOffsets[0], 3);

    for(MInt i = 0; i < grid().noDomains(); i++) {
      MInt ind = m_globalSndOffsets[i];
      for(MInt j = 0; j < rcvBufSizeGlob(i); j++) {
        MInt containingId = rcvBufGlob(ind + j);
        MInt cellId = cellIdsLoc[i][j];

        a_containingCell(cellId, b) = containingId;
        a_containingDomain(cellId, b) = i;
      }
    }
  }
}


/// \brief Change local into global ids.
///
/// \author Thomas Hoesgen
template <MInt nDim>
void LsCartesianSolver<nDim>::localToGlobalIds() {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  if(!m_reconstructOldG && m_LsRotate) {
    updateContainingGCells();
    localToGlobalIdsContainingCells();
  }
}

/// Reconstructs oldG for rotatigLS with reconstructOldG = true
///
/// \author Thomas Hoesgen
template <MInt nDim>
void LsCartesianSolver<nDim>::reconstructOldGField() {
  TRACE();

  if(domainId() == 0) cerr << "Reconstruct oldG Field...";

  mAlloc(m_geometryChange, m_noSets, "geometryChange", false, AT_);
  for(MInt i = 0; i < m_noBodiesToCompute; i++) {
    MInt body = m_bodiesToCompute[i];
    MInt set = m_bodyToSetTable[body];
    m_geometryChange[set] = true;
  }
  rotateSTL(1);
  constructGFieldFromSTL(6);
  /*
  constructGFieldFromSTL(5);
  if(m_buildCollectedLevelSetFunction) buildCollectedLevelSet(0);
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    for(MInt set = m_startSet; set < m_noSets; set++) {
      if(!m_geometryChange[set]) continue;
      a_oldLevelSetFunctionG(cellId, set) = a_levelSetFunctionG(cellId, set);
    }
  }
  */
  rotateSTL(-1);
  mDeallocate(m_geometryChange);
  // Reset the rotated angle
  for(MInt i = 0; i < m_noBodiesToCompute; i++) {
    MInt body = m_bodiesToCompute[i];
    for(MInt d = 0; d < nDim; d++) {
      m_semiLagrange_xRot_ref[body * nDim + d] = F0;
    }
  }
  resetContainingGCells();

  m_rotatingReinitTrigger = 0;
  if(domainId() == 0) cerr << " done" << endl;
}

/**
 * \brief Cartesian leaf-level exchange data with additional timer
 * \author Tim Wegmann
 *
 */
template <MInt nDim>
template <MBool currentLevelSet>
void LsCartesianSolver<nDim>::exchangeLeafDataLS() {
  TRACE();

  if(m_firstSolutionExchange) {
    RECORD_TIMER_START(m_timers[Timers::FirstEx]);
  }

  std::function<MFloat&(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt set) -> MFloat& {
    IF_CONSTEXPR(currentLevelSet) { return a_levelSetFunctionG(cellId, set); }
    else {
      return a_oldLevelSetFunctionG(cellId, set);
    }
  };

  this->exchangeLeafData(scalarField, m_maxNoSets);

  if(grid().azimuthalPeriodicity()) {
    IF_CONSTEXPR(currentLevelSet) { this->exchangeAzimuthalPer(&a_levelSetFunctionG(0, 0), m_maxNoSets); }
    else {
      this->exchangeAzimuthalPer(&a_oldLevelSetFunctionG(0, 0), m_maxNoSets);
    }
  }

  if(m_firstSolutionExchange) {
    RECORD_TIMER_STOP(m_timers[Timers::FirstEx]);
    m_firstSolutionExchange = false;
  }
}


/**
 * \brief exchange data on all levels with additional timer
 * \author Tim Wegmann
 *
 */
template <MInt nDim>
template <typename T>
void LsCartesianSolver<nDim>::exchangeDataLS(T* data, const MInt dataSize) {
  TRACE();

  if(m_firstSolutionExchange) {
    RECORD_TIMER_START(m_timers[Timers::FirstEx]);
  }

  exchangeData(data, dataSize);

  if(grid().azimuthalPeriodicity()) {
    this->exchangeAzimuthalPer(data, dataSize);
  }

  if(m_firstSolutionExchange) {
    RECORD_TIMER_STOP(m_timers[Timers::FirstEx]);
    m_firstSolutionExchange = false;
  }
}


/**
 * \brief Initializes the communication timers
 * \author Tim Wegmann
 */
template <MInt nDim>
void LsCartesianSolver<nDim>::initializeTimers() {
  TRACE();

  std::fill(m_timers.begin(), m_timers.end(), -1);

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "LS (solverId = " + to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::Solver], "total object lifetime", m_timerGroup);
  RECORD_TIMER_START(m_timers[Timers::Solver]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TimeInt], "Time-integration", m_timers[Timers::Solver]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FirstEx], "firstLSExchange", m_timers[Timers::TimeInt]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PostTime], "PostTimeStep", m_timers[Timers::Solver]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Finalize], "finalizeLS", m_timers[Timers::PostTime]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BuildTube], "buildLSTube", m_timers[Timers::Finalize]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SetBand], "setLSBand", m_timers[Timers::Finalize]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BuildMultiple], "buildMultipleLS", m_timers[Timers::Finalize]);
}

template class LsCartesianSolver<2>;
template class LsCartesianSolver<3>;
