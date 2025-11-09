// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lpt.h"
#include <chrono>
#include <cmath>
#include <memory>
#include <random>
#include <thread>
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "lpt_props.h"
#include "lptcollision.h"
#include "lptellipsoiddistance.h"

using namespace std;
using namespace maia;
using namespace maia::lpt;

/** Constructor for class LPT
 *
 *  \brief Constructor for LPT particles
 *
 *  \author Rudie Kunnen, Sven Berger
 *  \date   November 2010
 */
template <MInt nDim>
LPT<nDim>::LPT(const MInt solverId_, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm)
  : maia::CartesianSolver<nDim, LPT<nDim>>(solverId_, gridProxy_, comm, true), m_geometry(&geometry_) {
  TRACE();

  // create state files
  m_material = std::make_unique<MaterialState<nDim>>(solverId());

  // set wether the LPT-solver is used non-dimensionally
  // this needs to be known before the initialisation of the solver!
  m_nonDimensional = false;
  m_nonDimensional = Context::getSolverProperty<MBool>("nonDimensionaliseLPT", solverId(), AT_, &m_nonDimensional);

  // To set particle release Time for particles in isoTropic Turbulence
  // Release Timestep is determined by a Skewness of 0.44 of fluid velocity
  // skips Solution step until desired skewness is reached
  m_skewnessTimeStep = Context::getSolverProperty<MInt>("skewnessTimeStep", solverId(), AT_, &m_skewnessTimeStep);

  // creates the LPT Solver timers
  initializeTimers();
}

template <MInt nDim>
void LPT<nDim>::initSolver() {

  m_time = 0.0;

  // 1) read all properties
  readModelProps();

  // deleted in new version
  // readAdditionalProperties();

  m_cells.setLptCollectorCoupling(m_massCoupling, m_momentumCoupling, m_heatCoupling);

  if(m_evaporation) {
    m_cells.setLptCollectorNoSpecies(1);
  }

  if(m_ellipsoids) {
    m_cells.setLptCollectorSlopes();
  }

  // load properties related to intial generation of particles
  if(m_spawnParticles) {
    readSpawnProps();
  }

  if(m_momentumCoupling || m_heatCoupling) {
    readMomentumCouplingProps();
  }
  if(m_ellipsoids) {
    readEllipsoidProps();
  }

  initParticleVector();

  if(m_activeSecondaryBUp || m_activePrimaryBUp) {
    m_sprayModel = std::unique_ptr<SprayModel<nDim>>(new SprayModel<nDim>());
    m_sprayModel->init(this);
    m_log << "Spray model has been activated." << std::endl;
  }

  // If solver inactive only add dummy cells and mark all cells as halo cells
  // Inactive solvers can not have internal cells
  if(!isActive()) {
    m_cells.clear();
    m_cells.reset(grid().maxNoCells());
    this->setHaloCellsOnInactiveRanks();
    initBndryCells();
    return;
  }

  initCollector();

  initBndryCells();

  initGridProperties();

  findSpawnCellId();

  updateExchangeCells();
  grid().updateLeafCellExchange();
  initMPI();

  grid().findEqualLevelNeighborsParDiagonal(false);

  initialCondition();

  if(a_timeStepComputationInterval() > 0 && globalTimeStep % a_timeStepComputationInterval() == 0) {
    forceTimeStep(computeTimeStep());
  }

  countParticlesInCells();

  // LPT-statistics
  MInt noCells = a_noCells();
  MInt minCells = noCells;
  MInt maxCells = noCells;
  MInt noPart = a_noParticles();
  MInt minParts = noPart;
  MInt maxParts = noPart;

  MPI_Allreduce(MPI_IN_PLACE, &noCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noCells");
  MPI_Allreduce(MPI_IN_PLACE, &noPart, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noPart");
  MPI_Allreduce(MPI_IN_PLACE, &maxCells, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxCells");
  MPI_Allreduce(MPI_IN_PLACE, &minCells, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minCells");
  MPI_Allreduce(MPI_IN_PLACE, &maxParts, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxParts");
  MPI_Allreduce(MPI_IN_PLACE, &minParts, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minParts");

  if(domainId() == 0) {
    cerr << "LPT cell-statistics : " << minCells << " - " << maxCells << " avg " << noCells / noDomains() << endl;
    cerr << "LPT particle-statistics : " << minParts << " - " << maxParts << " avg " << noPart / noDomains() << endl;
  }
}

/**
 *  \brief init LPT cell collector
 *
 *  \author Tim Wegmann
 *  \date   October 2020
 */
template <MInt nDim>
void LPT<nDim>::initCollector() {
  TRACE();

  m_cells.clear();
  m_cells.reset(grid().maxNoCells());
  ASSERT(grid().tree().size() < grid().maxNoCells(), "Increase collector size!");
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    m_cells.append();
    a_isHalo(cellId) = false;
    a_isValidCell(cellId) = true;
    a_isWindow(cellId) = false;
    a_regridTrigger(cellId) = false;
    a_noParticlesInCell(cellId) = 0;
    a_noEllipsoidsInCell(cellId) = 0;
    a_bndryCellId(cellId) = -1;
    for(MInt n = 0; n < nDim; n++) {
      a_fluidVelocity(cellId, n) = 0.0;
    }
    a_fluidDensity(cellId) = 0.0;
    a_fluidPressure(cellId) = 0.0;
    a_fluidTemperature(cellId) = 0.0;
    if(m_evaporation) {
      a_fluidSpecies(cellId) = 0.0;
    }

    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(cellId, i) = 0;
      }
      a_workFlux(cellId) = 0;
    }
    if(m_heatCoupling) {
      a_heatFlux(cellId) = 0;
    }
    if(m_massCoupling) {
      a_massFlux(cellId) = 0;
    }
    if(m_ellipsoids) {
      for(MInt varId = 0; varId < nDim; varId++) {
        for(MInt dir = 0; dir < nDim; dir++) {
          a_velocitySlope(cellId, varId, dir) = F0;
        }
      }
    }
  }
  for(MInt cellId = grid().noInternalCells(); cellId < grid().tree().size(); cellId++) {
    m_cells.append();
    a_isHalo(cellId) = true;
    a_isWindow(cellId) = false;
    a_isValidCell(cellId) = true;
    a_regridTrigger(cellId) = false;
    a_noParticlesInCell(cellId) = 0;
    a_noEllipsoidsInCell(cellId) = 0;
    a_bndryCellId(cellId) = -1;
    for(MInt n = 0; n < nDim; n++) {
      a_fluidVelocity(cellId, n) = 0.0;
    }
    a_fluidDensity(cellId) = 0.0;
    a_fluidPressure(cellId) = 0.0;
    a_fluidTemperature(cellId) = 0.0;
    if(m_evaporation) {
      a_fluidSpecies(cellId) = 0.0;
    }
    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(cellId, i) = 0;
      }
      a_workFlux(cellId) = 0;
    }
    if(m_heatCoupling) {
      a_heatFlux(cellId) = 0;
    }
    if(m_massCoupling) {
      a_massFlux(cellId) = 0;
    }
    if(m_ellipsoids) {
      for(MInt varId = 0; varId < nDim; varId++) {
        for(MInt dir = 0; dir < nDim; dir++) {
          a_velocitySlope(cellId, varId, dir) = F0;
        }
      }
    }
  }
}

/**
 *  \fn void LPT::initParticleVector()
 *
 *  \brief init LPT spherical and ellipsoids particle vectors
 *         (i.e. set length and read non-dimensionalisation Parameters!)
 *  \author Tim Wegmann
 *  \date   August 2021
 */
template <MInt nDim>
void LPT<nDim>::initParticleVector() {
  TRACE();

  // if(isActive()) {
  m_partList.reserve(static_cast<MUlong>(m_maxNoParticles));
  if(m_ellipsoids) m_partListEllipsoid.reserve(static_cast<MUlong>(m_maxNoParticles));
  // }

  // set common/base static values

  /*! \page propertyPageLPT
      \section particleInterpolationOrder
      <code>MInt LPTBase::s_interpolationOrder</code>\n
      default = <code>2</code> \n \n
      Choose velocity-interpolation order: \n
      <ul>
      <li><code>2</code> Second-order</li>
      <li><code>3</code> Third-order</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  MInt intOrder = 2;
  intOrder = Context::getSolverProperty<MInt>("particleInterpolationOrder", solverId(), AT_, &intOrder);
  m_partList[0].s_interpolationOrder = intOrder;

  MInt intMethod = 1;
  intMethod = Context::getSolverProperty<MInt>("particleInterpolationMethod", solverId(), AT_, &intMethod);
  m_partList[0].s_interpolationMethod = intMethod;

  MFloat distFactorImp = 0.5;
  distFactorImp =
      Context::getSolverProperty<MFloat>("particleInterpolationDistFactor", solverId(), AT_, &distFactorImp);
  m_partList[0].s_distFactorImp = distFactorImp;

  m_partList[0].s_backPtr = this;

  // set ellipsoids static values
  if(m_ellipsoids) {
    m_partListEllipsoid[0].s_interpolationOrder = intOrder;
    m_partListEllipsoid[0].s_interpolationMethod = intMethod;
    m_partListEllipsoid[0].s_backPtr = this;
    m_partListEllipsoid[0].s_distFactorImp = distFactorImp;
  }


  // set particle static values:

  if(m_nonDimensional) {
    const MFloat Re = Context::getSolverProperty<MFloat>("Re", solverId(), AT_);
    m_partList[0].s_Re = Re;
    if(m_ellipsoids) m_partListEllipsoid[0].s_Re = Re;

    m_FrMag = 0;
    for(MInt i = 0; i < nDim; i++) {
      const MFloat Frm = Context::getSolverProperty<MFloat>("Frm", solverId(), AT_, i);
      m_partList[0].s_Frm[i] = Frm;
      if(m_ellipsoids) m_partListEllipsoid[0].s_Frm[i] = Frm;
      m_FrMag += POW2(Frm);
    }
    m_FrMag = sqrt(m_FrMag);

    if(m_heatCoupling || m_evaporation) {
      m_partList[0].s_Pr = Context::getSolverProperty<MFloat>("Pr", solverId(), AT_);
      m_partList[0].s_Sc = Context::getSolverProperty<MFloat>("Sc", solverId(), AT_);
    }
    if(m_activeSecondaryBUp) m_partList[0].s_We = Context::getSolverProperty<MFloat>("We", solverId(), AT_);

  } else {
    m_partList[0].s_Re = F1;
    if(m_ellipsoids) m_partListEllipsoid[0].s_Re = F1;
    MFloat defaultGravity = 0.0;
    MFloat gravity = 0;
    for(MInt i = 0; i < nDim; i++) {
      gravity = Context::getSolverProperty<MFloat>("particleGravity", solverId(), AT_, &defaultGravity, i);
      m_partList[0].s_Frm[i] = gravity;
      if(m_ellipsoids) m_partListEllipsoid[0].s_Frm[i] = gravity;
      m_FrMag += POW2(gravity);
    }
    m_FrMag = sqrt(m_FrMag);
  }

  /*! \page propertyPageLPT
      \section s_lengthFactor
      <code>MFloat LPTSpherical::s_lengthFactor</code>\n
      default = <code>1.0</code> \n \n
      Length factor conversion between the grid reference length as defined by the STL
      (aka STL-length) to the reference length used in the LPT solver!
      For numerical reasons and sake of clarity it is advised to use the
      reference-diameter of the particles (d_p0) as a reference length in the LPT-solver!
      L_ref / L_grid \n
      Keywords: <i>PARTICLE, NON-DIMENSIONALISATION,REFERENCE-LENGTH</i>
   */
  MFloat lengthFactor = F1;
  lengthFactor = Context::getSolverProperty<MFloat>("lengthFactor", solverId(), AT_, &lengthFactor);
  m_partList[0].s_lengthFactor = lengthFactor;
  if(m_ellipsoids) m_partListEllipsoid[0].s_lengthFactor = lengthFactor;
}

/**
 *  \brief init LPT boundary-cells used for the wall-collision
 *
 *  \author Tim Wegmann
 *  \date   October 2020
 */
template <MInt nDim>
void LPT<nDim>::initBndryCells() {
  TRACE();

  // if(!m_wallCollisions) return;

  // init data structure
  constexpr MInt maxNoSurfaces = 1;
  mAlloc(m_bndryCells, m_maxNoBndryCells, nDim, 0, 0, maxNoSurfaces, 0, "m_bndryCells", AT_);

  // fill on own data
  m_noOuterBndryCells = 0;
  m_bndryCells->resetSize(0);
}

/**
 *  \fn void LPT::initMPI()
 *
 *  \brief Init MPI buffers and communication datastructures.
 *  \author Sven Berger, Tim Wegmann
 *  \date   August 2015
 */
template <MInt nDim>
void LPT<nDim>::initMPI(const MBool fullReinit) {
  TRACE();

  if(noDomains() < 1) {
    return;
  }

  const MInt noNghbrDomains = grid().noNeighborDomains();

  if(fullReinit) {
    m_queueToSend.resize(noNghbrDomains);
    m_sendSize.resize(noNghbrDomains);
    m_recvSize.resize(noNghbrDomains);
    m_sendBuffer.resize(noNghbrDomains);
    m_recvBuffer.resize(noNghbrDomains);
    m_intSendBuffer.resize(noNghbrDomains);
    m_intRecvBuffer.resize(noNghbrDomains);
    m_mpi_reqSendFloat.resize(noNghbrDomains);
    m_mpi_reqRecvFloat.resize(noNghbrDomains);
    m_mpi_reqSendSize.resize(noNghbrDomains);
    m_mpi_reqSendInt.resize(noNghbrDomains);
    m_mpi_reqRecvInt.resize(noNghbrDomains);
    m_mpi_statusProbe.resize(noNghbrDomains);

    // if particles and ellipsoids can be present, the larger buffer size is chosen
    const MInt bufferSize = !m_ellipsoids ? bufSize<LPTSpherical<nDim>>()
                                          : mMax(bufSize<LPTSpherical<nDim>>(), bufSize<LPTEllipsoidal<nDim>>());
    const MInt intBufferSize = !m_ellipsoids
                                   ? intBufSize<LPTSpherical<nDim>>()
                                   : mMax(intBufSize<LPTSpherical<nDim>>(), intBufSize<LPTEllipsoidal<nDim>>());

    for(MInt i = 0; i < noNghbrDomains; i++) {
      m_sendSize[i] = 0;
      m_recvSize[i] = 0;
      m_sendBuffer[i] = make_unique<MFloat[]>(bufferSize);
      m_recvBuffer[i] = make_unique<MFloat[]>(bufferSize);
      m_intSendBuffer[i] = make_unique<MInt[]>(intBufferSize);
      m_intRecvBuffer[i] = make_unique<MInt[]>(intBufferSize);
    }

    // initialize the buffer
    for(MInt i = 0; i < noNghbrDomains; i++) {
      for(MInt j = 0; j < bufferSize; j++) {
        m_sendBuffer[i][j] = 0.0;
        m_recvBuffer[i][j] = 0.0;
      }
      for(MInt j = 0; j < intBufferSize; j++) {
        m_intSendBuffer[i][j] = 0;
        m_intRecvBuffer[i][j] = 0;
      }
    }


    if(m_nonBlockingComm) {
      m_noSourceTerms = 0;
      if(m_momentumCoupling) m_noSourceTerms++;
      if(m_momentumCoupling) m_noSourceTerms = m_noSourceTerms + nDim + 1;
      if(m_heatCoupling) m_noSourceTerms++;

      m_sourceSend.resize(noNghbrDomains);
      m_sourceRecv.resize(noNghbrDomains);

      if(m_sourceSendRequest != nullptr) mDeallocate(m_sourceSendRequest);
      mAlloc(m_sourceSendRequest, noNghbrDomains, "sourceSendReq", AT_);
      if(m_sourceRecvRequest != nullptr) mDeallocate(m_sourceRecvRequest);
      mAlloc(m_sourceRecvRequest, noNghbrDomains, "sourceRecvReq", AT_);

      if(m_flowSendRequest != nullptr) {
        mDeallocate(m_flowSendRequest);
      }
      mAlloc(m_flowSendRequest, noNghbrDomains, "flowSendReq", AT_);
      if(m_flowRecvRequest != nullptr) {
        mDeallocate(m_flowRecvRequest);
      }
      mAlloc(m_flowRecvRequest, noNghbrDomains, "flowRecvReq", AT_);

      if(m_checkAdaptationSendRequest != nullptr) {
        mDeallocate(m_checkAdaptationSendRequest);
      }
      mAlloc(m_checkAdaptationSendRequest, noDomains(), "checkAdapSendReq", AT_);
      if(m_checkAdaptationRecvRequest != nullptr) {
        mDeallocate(m_checkAdaptationRecvRequest);
      }
      mAlloc(m_checkAdaptationRecvRequest, noDomains(), "checkAdapRecvReq", AT_);

      m_checkAdaptationSend.resize(noDomains());
      m_checkAdaptationRecv.resize(noDomains());

      if(m_ellipsoids) {
        if(m_slopesSendRequest != nullptr) {
          mDeallocate(m_slopesSendRequest);
        }
        mAlloc(m_slopesSendRequest, noNghbrDomains, "slopesSendReq", AT_);
        if(m_slopesRecvRequest != nullptr) {
          mDeallocate(m_slopesRecvRequest);
        }
        mAlloc(m_slopesRecvRequest, noNghbrDomains, "slopesRecvReq", AT_);
      }
    }
  }

  if(m_nonBlockingComm) {
    MInt maxNoHaloCells = 0;
    MInt maxNoWindowCells = 0;

    ScratchSpace<MInt> windowCellsCnt(noNeighborDomains(), AT_, "noWindowCells");
    ScratchSpace<MInt> haloCellsCnt(noNeighborDomains(), AT_, "noHaloCells");

    for(MInt i = 0; i < noNghbrDomains; i++) {
      maxNoHaloCells = mMax(maxNoHaloCells, grid().noLeafHaloCells(i));
      maxNoWindowCells = mMax(maxNoWindowCells, grid().noLeafWindowCells(i));
      m_sourceSendRequest[i] = MPI_REQUEST_NULL;
      m_sourceRecvRequest[i] = MPI_REQUEST_NULL;
      m_flowSendRequest[i] = MPI_REQUEST_NULL;
      m_flowRecvRequest[i] = MPI_REQUEST_NULL;
      if(m_ellipsoids) {
        m_slopesSendRequest[i] = MPI_REQUEST_NULL;
        m_slopesRecvRequest[i] = MPI_REQUEST_NULL;
      }
      windowCellsCnt[i] = grid().noLeafWindowCells(i);
      haloCellsCnt[i] = grid().noLeafHaloCells(i);
    }

    for(MInt d = 0; d < noDomains(); d++) {
      m_checkAdaptationRecvRequest[d] = MPI_REQUEST_NULL;
      m_checkAdaptationSendRequest[d] = MPI_REQUEST_NULL;
    }

    if(m_flowRecv != nullptr) {
      mDeallocate(m_flowRecv);
    }
    mAlloc(m_flowRecv, noNeighborDomains(), &haloCellsCnt[0], PV.noVars() + m_evaporation, "m_flowRecv", AT_);
    if(m_flowSend != nullptr) {
      mDeallocate(m_flowSend);
    }
    mAlloc(m_flowSend, noNeighborDomains(), &windowCellsCnt[0], PV.noVars() + m_evaporation, "m_flowSend", AT_);


    for(MInt i = 0; i < noNghbrDomains; i++) {
      m_sourceSend[i] = make_unique<MFloat[]>(maxNoHaloCells * m_noSourceTerms);
      m_sourceRecv[i] = make_unique<MFloat[]>(maxNoWindowCells * m_noSourceTerms);
    }


    for(MInt i = 0; i < noNghbrDomains; i++) {
      for(MInt j = 0; j < grid().noLeafHaloCells(i) * m_noSourceTerms; j++) {
        m_sourceSend[i][j] = 0.0;
      }
      for(MInt j = 0; j < grid().noLeafWindowCells(i) * m_noSourceTerms; j++) {
        m_sourceRecv[i][j] = 0.0;
      }
    }
    for(MInt i = 0; i < noDomains(); i++) {
      m_checkAdaptationSend[i] = 0;
      m_checkAdaptationRecv[i] = 0;
    }

    if(m_ellipsoids) {
      if(m_slopesRecv != nullptr) {
        mDeallocate(m_slopesRecv);
      }
      mAlloc(m_slopesRecv, noNeighborDomains(), &haloCellsCnt[0], nDim * nDim, "m_slopesRecv", AT_);
      if(m_slopesSend != nullptr) {
        mDeallocate(m_slopesSend);
      }
      mAlloc(m_slopesSend, noNeighborDomains(), &windowCellsCnt[0], nDim * nDim, "m_slopesSend", AT_);
    }
  }

  if(m_respawn && fullReinit) {
    struct {
      MInt val;
      MInt rank;
    } local{}, global{};
    local.val = (MInt)m_respawnCells.size();
    local.rank = domainId();
    global.val = 0;
    global.rank = -1;
    MPI_Allreduce(&local, &global, 1, MPI_2INT, MPI_MAXLOC, mpiComm(), AT_, "local", "global");
    m_respawnDomain = global.rank;
    MInt* particleRespawnCellsPerDomain = nullptr;
    if(domainId() == m_respawnDomain) {
      particleRespawnCellsPerDomain = new MInt[noDomains()];
    }
    MPI_Gather(&local.val, 1, MPI_INT, particleRespawnCellsPerDomain, 1, MPI_INT, m_respawnDomain, mpiComm(), AT_,
               "local.val", "particleRespawnCellsPerDomain");


    if(domainId() == m_respawnDomain) {
      MInt globalNoRespawnCells = 0;
      m_noRespawnDomains = 0;
      for(MInt i = 0; i < noDomains(); ++i) {
        globalNoRespawnCells += particleRespawnCellsPerDomain[i];
        if(particleRespawnCellsPerDomain[i] > 0) {
          m_noRespawnDomains++;
        }
      }
      if(globalNoRespawnCells == 0) {
        stringstream errorMessage;
        errorMessage << "globally no RespawnCells found despite the set properties particleRespawn " << m_respawn
                     << " and general Respawn (particleRespawnType==0) or massiveParallel ";
        mTerm(1, AT_, errorMessage.str());
      }
      m_respawnGlobalDomainOffsets = new MInt[m_noRespawnDomains + 1];
      m_respawnGlobalDomainOffsets[0] = 0;
      m_respawnDomainRanks = new MInt[m_noRespawnDomains];
      MInt counter = 0;
      for(MInt i = 0; i < noDomains(); ++i) {
        if(particleRespawnCellsPerDomain[i] > 0) {
          m_respawnDomainRanks[counter] = i;
          m_respawnGlobalDomainOffsets[counter + 1] =
              m_respawnGlobalDomainOffsets[counter] + particleRespawnCellsPerDomain[i];
          ++counter;
        }
      }

      delete[] particleRespawnCellsPerDomain;
    } else {
      m_noRespawnDomains = 1;
    }
  }
}

/**
 * \brief calculate Grid Properties for periodic boundaries, adapted from rigedbodies
 * \author Oliver Groll, Johannes Grafen
 * \date March 2022
 */
template <MInt nDim>
void LPT<nDim>::initGridProperties() {
  computeBoundingBox(m_globBbox.begin());
  MInt periodic = 0;
  for(MInt n = 0; n < nDim; n++) {
    m_globDomainLength[n] = m_globBbox[n + nDim] - m_globBbox[n];
    periodic += grid().raw().periodicCartesianDir(n);
  }

  if(periodic) {
    m_periodicBC = true;
  }
}

/**
 * \fn void LPT::initialCondition()
 *
 * \brief Initialize particles
 * \author Christoph Siewert, Sven Berger, Tim Wegmann
 * \date   June 2015
 */
template <MInt nDim>
void LPT<nDim>::initialCondition() {
  TRACE();

  if(m_restart) {
    const MInt noParticles = loadParticleRestartFile();
    m_log << noParticles << " particles loaded from restart-file." << endl;

    // particle
    calculateTerminalVelocities();

    checkParticles();

  } else {
    MInt noParticles = 0;

    // particle initialization methods
    //  0 = general init (default)
    //  1 = centered in all inflow cells
    //  2,3 = provide 'part.txt' file with diameter, density ratio and coordinates
    //      (the corresponding cellIds are determined, velocity is zero)
    //  4 = provide 'part.txt' file with diameter, density ratio, coordinates
    //      and velocity (the corresponding cellIds are determined)
    //  5 = no initialisation (only primary spray injection)


    switch(m_initializationMethod) {
      // general particle initialisation
      case 0: {
        MFloat* dia = nullptr;
        MFloat* particleDiameters = nullptr;
        MInt sumOfRates = 0;

        // ellipsoids
        MFloat particleEllipsoidInitNoPartperCell = F0;
        vector<pair<MFloat, MFloat>> ellipsoidsInitVector;
        MInt noEllipsoidCat = 0;

        /*! \page propertyPageLPT
        \section particleInitNoPartperCell
        <code>MFloat LPT::particleInitNoPartperCell</code>\n
        default = <code>0.0</code> \n \n
        number of particles to be initialized per cell. \n
        It can be a double value not only a MInt, but if it is double,
        the number of initialized Particle depends on the number of processors
        in massive parallel concept \n
        Keywords: <i>PARTICLE</i>
        */
        MFloat particleInitNoPartperCell = 0.0;
        particleInitNoPartperCell = Context::getSolverProperty<MFloat>("particleInitNoPartperCell", solverId(), AT_,
                                                                       &particleInitNoPartperCell);

        if(particleInitNoPartperCell > 0) {
          m_log << "LPT.cpp:initialCondition(): Initializing " << particleInitNoPartperCell << " particles per cell."
                << endl;

          /*! \page propertyPageLPT
              \section particleDiameters
              <code>MFloat particleDiameters</code>\n
              default = <code>(0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007,
             0.00008, 0.00009, 0.00010)</code> \n \n diameters of particles. \n Keywords:
             <i>PARTICLE</i>
           */
          array<MFloat, 10> particleDiametersDefaults = {0.00001, 0.00002, 0.00003, 0.00004, 0.00005,
                                                         0.00006, 0.00007, 0.00008, 0.00009, 0.00010};
          MInt noDifferentPart = 10;
          if(Context::propertyExists("particleDiameters", solverId())) {
            noDifferentPart = Context::propertyLength("particleDiameters", solverId());
            particleDiameters = new MFloat[noDifferentPart];
            for(MInt i = 0; i < noDifferentPart; i++) {
              particleDiameters[i] = Context::getSolverProperty<MFloat>("particleDiameters", solverId(), AT_, i);
            }
          } else {
            particleDiameters = new MFloat[noDifferentPart];
            copy(&particleDiametersDefaults[0], &particleDiametersDefaults[9], particleDiameters);
          }
          /*! \page propertyPageLPT
              \section particleInitRatesOfDiameters
              <code>MInt(*) particleInitQuotasOfDiameters</code>\n
              default = <code>1</code> \n \n
              quota of particles with associated diameter. \n
              Keywords: <i>PARTICLE</i>
           */
          MIntScratchSpace particleInitQuotasOfDiameters(noDifferentPart, AT_, "particleInitQuotasOfDiameters");
          MInt count = 0;
          if(Context::propertyExists("particleInitRatesOfDiameters", solverId())) {
            count = Context::propertyLength("particleInitRatesOfDiameters", solverId());
            if(count == noDifferentPart) {
              for(MInt i = 0; i < count; i++) {
                particleInitQuotasOfDiameters[i] =
                    Context::getSolverProperty<MInt>("particleInitRatesOfDiameters", solverId(), AT_, i);
              }
              // const MInt* partQuotas =
              // *(Context::getSolverProperty<MInt>("particleInitRatesOfDiameters",
              //                                                    solverId(), AT_));
              // const MInt* partQuotas = nullptr;
              // copy(partQuotas, &partQuotas[noDifferentPart], &particleInitQuotasOfDiameters[0]);
            } else {
              stringstream errorMessage;
              errorMessage << " LPT::initialize() with m_particleInitializationMethod"
                           << " = -1 : Try to initialize particles with " << noDifferentPart
                           << " different diameters (property particleDiameters), but found " << count
                           << " different rates (property particleInitRatesOfDiameters) " << endl;
              mTerm(1, AT_, errorMessage.str());
            }
          } // otherwise define diameters to be equally likely
          else {
            for(MInt i = 0; i < noDifferentPart; ++i) {
              particleInitQuotasOfDiameters[i] = 1;
            }
          }

          // to randomly distribute the diameters among particles declare array
          // dia with the number of diameters being dependent on the previously
          // declared quotas
          sumOfRates =
              accumulate(&particleInitQuotasOfDiameters[0], &particleInitQuotasOfDiameters[0] + noDifferentPart, 0);
          dia = new MFloat[sumOfRates];
          MInt counter = 0;
          for(MInt i = 0; i < noDifferentPart; ++i) {
            for(MInt j = 0; j < particleInitQuotasOfDiameters[i]; ++j) {
              dia[counter] = particleDiameters[i];
              counter++;
            }
          }
          shuffle(dia, &dia[sumOfRates], randomInit());

          calculateTerminalVelocities();
        }

        if(m_ellipsoids) {
          MFloat* ellipsoidSemiMinorAxis = nullptr;
          MFloat* ellipsoidAspectRatios = nullptr;
          const MInt noEllipsoidCatDefault = 6;
          array<MFloat, noEllipsoidCatDefault> ellipsoidSemiMinorAxisDefaults = {
              0.000031748, 0.000039850263, 0.00002773445, 0.000034668, 0.0000259842, 0.000031498026};
          array<MFloat, noEllipsoidCatDefault> ellipsoidsAspectRatioDefaults = {2.0, 2.0, 3.0, 3.0, 4.0, 4.0};

          particleEllipsoidInitNoPartperCell = F0;
          particleEllipsoidInitNoPartperCell = Context::getSolverProperty<MFloat>(
              "particleEllipsoidInitNoPartperCell", solverId(), AT_, &particleEllipsoidInitNoPartperCell);

          if(particleEllipsoidInitNoPartperCell > 0) {
            m_log << "LPT.cpp:initialCondition(): Initializing " << particleEllipsoidInitNoPartperCell
                  << " ellipsoids per cell." << endl;

            if(Context::propertyExists("particleEllipsoidSemiMinorAxis", solverId())) {
              noEllipsoidCat = Context::propertyLength("particleEllipsoidSemiMinorAxis", solverId());
              if(noEllipsoidCat != Context::propertyLength("particleEllipsoidAspectRatio", solverId())) {
                mTerm(1, AT_,
                      "Error the properties particleEllipsoidSemiMinorAxis and particleEllipsoidAspectRatio "
                      " do not have the same count!");
              }
              ellipsoidSemiMinorAxis = new MFloat[noEllipsoidCat];
              ellipsoidAspectRatios = new MFloat[noEllipsoidCat];
              for(MInt i = 0; i < noEllipsoidCat; i++) {
                ellipsoidSemiMinorAxis[i] =
                    Context::getSolverProperty<MFloat>("particleEllipsoidSemiMinorAxis", solverId(), AT_, i);
                ellipsoidAspectRatios[i] =
                    Context::getSolverProperty<MFloat>("particleEllipsoidAspectRatio", solverId(), AT_, i);
              }
            } else {
              noEllipsoidCat = noEllipsoidCatDefault;
              ellipsoidSemiMinorAxis = new MFloat[noEllipsoidCat];
              ellipsoidAspectRatios = new MFloat[noEllipsoidCat];
              copy(&ellipsoidSemiMinorAxisDefaults[0], &ellipsoidSemiMinorAxisDefaults[noEllipsoidCat - 1],
                   ellipsoidSemiMinorAxis);
              copy(&ellipsoidsAspectRatioDefaults[0], &ellipsoidsAspectRatioDefaults[noEllipsoidCat - 1],
                   ellipsoidAspectRatios);
            }

            MIntScratchSpace ellipsoidInitQuotasPerCat(noEllipsoidCat, AT_, "particleInitQuotasOfDiameters");
            MInt count = 0;
            if(Context::propertyExists("particleEllipsoidInitRatesPerCat", solverId())) {
              count = Context::propertyLength("particleEllipsoidInitRatesPerCat", solverId());
              if(count == noEllipsoidCat) {
                for(MInt i = 0; i < count; i++) {
                  ellipsoidInitQuotasPerCat[i] =
                      Context::getSolverProperty<MInt>("particleEllipsoidInitRatesPerCat", solverId(), AT_, i);
                }
              } else {
                stringstream errorMessage;
                errorMessage << " LPT::initialize() with m_particleInitializationMethod"
                             << " = -1 : Try to initialize particles with " << noEllipsoidCat
                             << " different ellipsoidal categories, but found " << count
                             << " different rates (property particleEllipsoidInitRatesPerCat) " << endl;
                mTerm(1, AT_, errorMessage.str());
              }
            } else { // otherwise define cats to be equally likely
              for(MInt i = 0; i < noEllipsoidCat; ++i) {
                ellipsoidInitQuotasPerCat[i] = 1;
              }
            }

            for(MInt i = 0; i < noEllipsoidCat; ++i) {
              for(MInt j = 0; j < ellipsoidInitQuotasPerCat[i]; ++j) {
                ellipsoidsInitVector.emplace_back(ellipsoidSemiMinorAxis[i], ellipsoidAspectRatios[i]);
              }
            }
            shuffle(ellipsoidsInitVector.begin(), ellipsoidsInitVector.end(), randomInit());
          }
          delete[] ellipsoidSemiMinorAxis;
          delete[] ellipsoidAspectRatios;
        } // end of ellipsoids properties reading

        if(particleEllipsoidInitNoPartperCell > 0 || particleInitNoPartperCell > 0) {
          /*! \page propertyPageLPT
              \section particleInitOffsets
              <code>MFloat(*) particleInitOffsets</code>\n
              default = <code>bounding box == no Offsets</code> \n \n
              provide offsets to restrict particle initialization \n
              Keywords: <i>PARTICLE</i>
           */
          MFloatScratchSpace particleInitOffsets(nDim * 2, AT_, "particleInitOffsets");
          m_geometry->getBoundingBox(&particleInitOffsets[0]);
          for(MInt i = 0; i < nDim * 2; ++i) {
            particleInitOffsets[i] =
                Context::getSolverProperty<MFloat>("particleInitOffsets", solverId(), AT_, &particleInitOffsets[i], i);
          }

          MInt teller = 0;
          MInt sayer = 0;
          if(m_ellipsoids) sayer = domainId() % noEllipsoidCat;
          MInt noParticleCells = 0;
          MBool randomLocWthCell = true;
          /*! \page propertyPageLPT
            \section particleRandLocWthCell
            <code>MBool particleRandomizeLocWthCell</code>\n
            default = true \n \n
            Randomize the location of the added particles within each cell. \n
            Keywords: <i>PARTICLE</i>
         */
          randomLocWthCell =
              Context::getSolverProperty<MBool>("particleRandLocWthCell", solverId(), AT_, &randomLocWthCell);

          /*! \page propertyPageLPT
               \section particleInitVelocity
               <code>MFloat[dim] particleInitVelocity</code>\n
               default = <code>0.0, 0.0, 0.0</code> \n \n
               Velocity of the particles to be initialized.\n
               Keywords: <i>PARTICLE</i>
               */
          array<MFloat, nDim> particleInitVelocity{};
          for(MInt i = 0; i < nDim; ++i) {
            particleInitVelocity.at(i) = Context::getSolverProperty<MFloat>("particleInitVelocity", solverId(), AT_,
                                                                            &particleInitVelocity.at(i), i);
          }

          // Note: particleInitOffsets is multiplied by referenceLength
          for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
            if(!c_isLeafCell(cellId)) continue;
            // don't use bndryCells!
            MBool hasNoNeighbor = false;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(!c_hasNeighbor(cellId, m_revDir[dir])
                 && (c_parentId(cellId) == -1 || !c_hasNeighbor(c_parentId(cellId), m_revDir[dir]))) {
                hasNoNeighbor = true;
                break;
              }
            }
            if(hasNoNeighbor) continue;
            // check if particle is in bounding box
            MBool inside = true;
            for(MInt direction = 0; direction < nDim; ++direction) {
              if(c_coordinate(cellId, direction) < particleInitOffsets[direction]
                 || c_coordinate(cellId, direction) > particleInitOffsets[nDim + direction]) {
                inside = false;
              }
            }
            if(!inside) continue;
            if(inside) {
              ++noParticleCells;
              while((m_partList.size() / MFloat(noParticleCells)) < particleInitNoPartperCell) {
                addParticle(cellId, dia[teller], m_material->densityRatio(), static_cast<MInt>(randomLocWthCell));

                m_partList.back().m_velocity = particleInitVelocity;

                // loops through the diameters
                teller = (teller + 1) % sumOfRates;
                if(teller == 0) {
                  shuffle(dia, &dia[sumOfRates], randomInit());
                }
                noParticles++;
              }

              if(m_ellipsoids) {
                while((m_partListEllipsoid.size() / MFloat(noParticleCells)) < particleEllipsoidInitNoPartperCell) {
                  addEllipsoid(cellId, ellipsoidsInitVector[sayer].first, ellipsoidsInitVector[sayer].second,
                               m_material->densityRatio(), static_cast<MInt>(randomLocWthCell));

                  m_partListEllipsoid.back().m_velocity = particleInitVelocity;

                  sayer = (sayer + 1) % ellipsoidsInitVector.size(); // loops through the cats
                  if(sayer == 0) {
                    shuffle(ellipsoidsInitVector.begin(), ellipsoidsInitVector.end(), randomInit());
                  }
                  ++noParticles;
                }
              }
            }
          }
          cerr << domainId() << ": Particles initialized: " << noParticles << endl;
          exchangeOffset(noParticles);
        }

        delete[] particleDiameters;
        break;
      }
      case 1: {
        MFloat diam = 0.1;
        MInt inflowDir = 1;
        /*! \page propertyPageLPT
          \section inflowDir
          <code>MInt particleInflowDir</code>\n
          default = <code>1.0</code> \n \n
          Defines the inflow direction of particles initialized at the inflow! <i>PARTICLE</i>
        */
        inflowDir = Context::getSolverProperty<MInt>("particleInflowDir", solverId(), AT_, &inflowDir);
        MInt revDir = m_revDir[inflowDir];
        for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
          if(!c_isLeafCell(cellId)) continue;
          if(!c_hasNeighbor(cellId, revDir)
             && (c_parentId(cellId) == -1 || !c_hasNeighbor(c_parentId(cellId), revDir))) {
            addParticle(cellId, diam, m_material->densityRatio(), 0);

            // correct initial velocity!
            // The position is already at the cartesian cell center!
            for(MInt t = 0; t < nDim; t++) {
              m_partList.back().m_oldVel.at(t) = 0;
              m_partList.back().m_velocity.at(t) = 0;
            }
            noParticles++;
          }
        }
        cerr << "domainId " << domainId() << " created " << noParticles << endl;

        if(m_ellipsoids) mTerm(1, AT_, "Init mode 1 not yet implemented for ellipsoidal particles!");

        exchangeOffset(noParticles);
        break;
      }

      case 2:
      case 3: // particles without defined velocity
      case 4: {
        ifstream readFile;
        array<MFloat, nDim> thisCoord{};
        array<MFloat, nDim> thisVel{};
        MFloat thisDiameter = NAN;
        MFloat density = NAN;

        readFile.open("part.txt", ios_base::in);
        if(!readFile) {
          mTerm(1, AT_, "Error reading particle input file 'part.txt'.");
        }

        for(;;) { // read the particle parameters from file:
          //   3D     2D
          //  (1)    (1)     diameter
          //  (2)    (2)     density
          // (3-5)  (3-4)    x, y(, z) coordinates
          // (6-8)  (5-6)    u, v(, w) velocities
          readFile >> thisDiameter;
          if((readFile.rdstate() & ifstream::eofbit) != 0) {
            break;
          }
          readFile >> density;
          for(MInt i = 0; i < nDim; i++) {
            readFile >> thisCoord.at(i);
          }

          if(m_initializationMethod == 4) {
            for(MInt i = 0; i < nDim; i++) {
              readFile >> thisVel.at(i);
            }
          }

          // find cellId of new particle
          MInt thisCellId = grid().findContainingLeafCell(&thisCoord[0]);

          // it is possible that there is more than one valid location to eliminate the possibility
          // of creating duplicate particles on multiple domains check via mpi
          if(noDomains() > 1) {
            MInt spawnDomain = domainId();
            MInt minDomain = -1;
            if(thisCellId == -1) {
              spawnDomain = std::numeric_limits<MInt>::max();
            }

            // Only create particle on the domain with the lowest Id
            MPI_Allreduce(&spawnDomain, &minDomain, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "spawnDomain", "minDomain");

            if(minDomain != domainId()) {
              continue;
            }
          }

          // particle is not inside this domain
          if(thisCellId == -1) {
            if(noDomains() == 1) {
              m_log << "Particle initialisation file contains particles outside of the geometry." << endl;
            }
            continue;
          }

          addParticle(thisCellId, thisDiameter, density);
          for(MInt i = 0; i < nDim; i++) {
            m_partList.back().m_position.at(i) = thisCoord.at(i);
          }
          for(MInt i = 0; i < nDim; i++) {
            m_partList.back().m_velocity.at(i) = thisVel.at(i);
          }
          for(MInt i = 0; i < nDim; i++) {
            m_partList.back().m_accel.at(i) = 0.0;
          }
          for(MInt i = 0; i < nDim; i++) {
            m_partList.back().m_oldPos.at(i) = thisCoord.at(i);
          }
          for(MInt i = 0; i < nDim; i++) {
            m_partList.back().m_oldVel.at(i) = thisVel.at(i);
          }
          noParticles++;
        }
        readFile.close();

        m_log << noParticles << " particles from file introduced in this solver." << endl;

        if(m_ellipsoids) {
          MInt noEllipsoids = 0;
          ifstream readFileE;
          array<MFloat, nDim> thisCoordE{};
          array<MFloat, nDim> thisVelE{};
          array<MFloat, nDim> thisOrientE{};
          MFloat thisAspectRatio = NAN;
          MFloat thisSemiMinorAxis = NAN;
          MFloat densityE = NAN;

          readFile.open("ellipsoids.txt", ios_base::in);
          if(!readFile) {
            mTerm(1, AT_, "Error reading ellipsoidal particle input file 'ellipsoids.txt'.");
          }

          for(;;) { // read the ellipsoidal particle parameters from file:
            //   3D      2D
            //   (1)     (1)     aspect ratio
            //   (2)     (2)     semi minor axis
            //   (3)     (3)     density ratio
            //  (4-6)   (4-5)    x, y (, z) coordinates
            //  (7-9)   (6-7)    x, y (, z) orientation
            // (10-12)  (8-9)    u, v (, w) velocities
            readFile >> thisAspectRatio;
            readFile >> thisSemiMinorAxis;
            if((readFile.rdstate() & ifstream::eofbit) != 0) {
              break;
            }
            readFile >> densityE;
            for(MInt i = 0; i < nDim; i++) {
              readFile >> thisCoordE.at(i);
            }

            if(m_initializationMethod >= 3) {
              for(MInt i = 0; i < nDim; i++) {
                readFile >> thisOrientE.at(i);
                thisOrientE.at(i) = thisOrientE.at(i) * PI / 180.0;
              }
            }
            if(m_initializationMethod == 4) {
              for(MInt i = 0; i < nDim; i++) {
                readFile >> thisVelE.at(i);
              }
            }

            // find cellId of new particle
            MInt thisCellId = grid().findContainingLeafCell(&thisCoordE[0]);

            // it is possible that there is more than one valid location to eliminate the possibility
            // of creating duplicate particles on multiple domains check via mpi
            if(noDomains() > 1) {
              MInt spawnDomain = domainId();
              MInt minDomain = -1;
              if(thisCellId == -1) {
                spawnDomain = std::numeric_limits<MInt>::max();
              }

              // Only create particle on the domain with the lowest Id
              MPI_Allreduce(&spawnDomain, &minDomain, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "spawnDomain", "minDomain");

              if(minDomain != domainId()) {
                continue;
              }
            }

            // particle is not inside this domain
            if(thisCellId == -1) {
              if(noDomains() == 1) {
                m_log << "Ellipsoidal particle initialisation file contains particles outside of the geometry." << endl;
              }
              continue;
            }

            addEllipsoid(thisCellId, thisSemiMinorAxis, thisAspectRatio, densityE);
            for(MInt i = 0; i < nDim; i++) {
              m_partListEllipsoid.back().m_position.at(i) = thisCoordE.at(i);
              m_partListEllipsoid.back().m_oldPos.at(i) = thisCoordE.at(i);
              m_partListEllipsoid.back().m_velocity.at(i) = thisVelE.at(i);
              m_partListEllipsoid.back().m_oldVel.at(i) = thisVelE.at(i);
              m_partListEllipsoid.back().m_accel.at(i) = F0;
              m_partListEllipsoid.back().m_oldAccel.at(i) = F0;
              m_partListEllipsoid.back().m_angularVel.at(i) = F0;
              m_partListEllipsoid.back().m_oldAngularVel.at(i) = F0;
              m_partListEllipsoid.back().m_angularAccel.at(i) = F0;
              m_partListEllipsoid.back().m_oldAngularAccel.at(i) = F0;
            }
            if(m_initializationMethod >= 3) {
              maia::math::rotation2quaternion(thisOrientE.begin(), m_partListEllipsoid.back().m_quaternion.begin());
            } else {
              // Rotation to identity rotation matrix
              m_partListEllipsoid.back().m_quaternion.at(0) = F1; // qw
              m_partListEllipsoid.back().m_quaternion.at(1) = F0; // qx
              m_partListEllipsoid.back().m_quaternion.at(2) = F0; // qy
              m_partListEllipsoid.back().m_quaternion.at(3) = F0; // qz
            }
            for(MInt i = 0; i < 4; i++) {
              m_partListEllipsoid.back().m_oldQuaternion.at(i) = m_partListEllipsoid.back().m_quaternion.at(i);
            }

            noEllipsoids++;
            noParticles++;
          }
          readFile.close();

          m_log << noEllipsoids << " ellipsoidal particles from file introduced in this solver." << endl;
        }
        exchangeOffset(noParticles);
        break;
      }

      case 5: {
        ASSERT(m_activePrimaryBUp, "");
        break;
      }
      default:
        mTerm(1, AT_,
              "LPT::initialize() switch variable 'm_particleInitializationMethod' not "
              "matching any case");
    }
  }
}


/*  \brief runified run loop function call
 *
 *  \date October-2020
 *  \author Sven Berger
 */
template <MInt nDim>
void LPT<nDim>::finalizeInitSolver() {
  TRACE();

  if(m_spawnParticles && !m_restart) {
    m_particleResiduum = 1.0;
  }

  if(!isActive()) return;

  if(domainId() == 0) {
    cerr << "Finding cartesian neighbors for LPT-solver!" << endl;
  }
  grid().findCartesianNghbIds();


  if(domainId() == 0) {
    initSummary();
  }

  if(!m_restart) {
    if(domainId() == 0) {
      cerr << "Writing initial particle file" << endl;
    }
    writePartData();
  }

  setRegridTrigger();
}


/**
 * \brief  prepare solver adaptation
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::prepareAdaptation() {
  TRACE();

  ASSERT(m_freeIndices.empty(), "");
  m_adaptationLevel = maxUniformRefinementLevel();
  m_forceAdaptation = false;

  if(!isActive()) {
    return;
  }

  receiveFlowField();
  receiveVelocitySlopes();
  waitForSendReqs();

  if(globalTimeStep > 0) {
    countParticlesInCells();
  }

  // add refinement around spawnCellId!
  // NOTE: a_noParticlesInCell will be corrected after the adaptation!
  if(m_spawnCellId > -1) {
    a_noParticlesInCell(m_spawnCellId) += 20;
    if(m_ellipsoids) a_noEllipsoidsInCell(m_spawnCellId) += 20;
  }

  this->exchangeData(&a_noParticlesInCell(0), 1);
  if(m_ellipsoids) this->exchangeData(&a_noEllipsoidsInCell(0), 1);


  // add fluxes to window-cells
  if(!m_nonBlockingComm) {
    exchangeSourceTerms();
  } else {
    sendSourceTerms();
    receiveSourceTerms();
    waitForSendReqs();
  }
  // reset fluxes on halo-cells
  resetSourceTerms(noInternalCells());

  // set bndry-Cell Id for adaptation on leaf-cells and all parants
  // and exchange information, to set this on halo-cells and
  // partition level-shift ancestors in a different rank
  // as required for interface sensor!!
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_bndryCellId(cellId) = -1;
  }
  for(MInt id = 0; id < m_bndryCells->size(); id++) {
    const MInt cellId = m_bndryCells->a[id].m_cellId;
    a_bndryCellId(cellId) = 1;
    MInt parentId = c_parentId(cellId);
    while(parentId > -1 && parentId < a_noCells()) {
      a_bndryCellId(parentId) = 1;
      parentId = c_parentId(parentId);
    }
  }
  this->exchangeData(&a_bndryCellId(0), 1);

  // remove all bndryCells
  if(m_wallCollisions) {
    m_bndryCells->resetSize(0);
    m_noOuterBndryCells = 0;
  }
}

/**
 * \brief set the sensors for a single adaptation step
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::setSensors(std::vector<std::vector<MFloat>>& sensors,
                           std::vector<MFloat>& sensorWeight,
                           std::vector<std::bitset<64>>& sensorCellFlag,
                           std::vector<MInt>& sensorSolverId) {
  TRACE();

  MInt noSensors = this->m_noInitialSensors;
  if(globalTimeStep > 0) noSensors = this->m_noSensors;

  // If solver is inactive the sensor arrays still need to be set to optain the
  // correct offsets
  const MInt sensorOffset = (signed)sensors.size();
  ASSERT(sensorOffset == 0 || grid().raw().treeb().noSolvers() > 1, "");
  sensors.resize(sensorOffset + noSensors, vector<MFloat>(grid().raw().m_noInternalCells, F0));
  sensorWeight.resize(sensorOffset + noSensors, -1);
  sensorCellFlag.resize(grid().raw().m_noInternalCells, sensorOffset + noSensors);
  sensorSolverId.resize(sensorOffset + noSensors, solverId());
  ASSERT(sensorOffset + noSensors < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

  if(domainId() == 0) {
    cerr << "Setting " << noSensors << " sensors for LPT-Solver adaptation." << endl;
  }

  for(MInt sen = 0; sen < noSensors; sen++) {
    sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
  }

  ASSERT(m_freeIndices.empty(), "");
  m_freeIndices.clear();

  if(!isActive()) return;

  for(MInt sen = 0; sen < noSensors; sen++) {
    (this->*(this->m_sensorFnPtr[sen]))(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen);
  }

  /*
  if(m_spawnCellId > -1) {
    const MInt gridCellId = grid().tree().solver2grid(m_spawnCellId);
    const MInt partent = c_parentId(m_spawnCellId);
    const MInt parentGridCellId = partent > -1 ? grid().tree().solver2grid(partent) : -1;
    if(parentGridCellId > -1 && sensors[sensorOffset][parentGridCellId] < 1) {
      if(sensors[sensorOffset][gridCellId] < 1) {
        cerr << "Injector is not refined!" << endl;
      }
    }
  }
  */
}


/** \page sensorsLPT
 *
 *  \section part Particles
 *
 *  This sensor represents the number of particles.<br>
 *  Property: <code>PARTICLES</code>
 */
/**
 * \brief set the sensors around the particle positions
 *       to get the markedCells
 * \author Tim Wegmann, Sven Berger
 */
template <MInt nDim>
void LPT<nDim>::sensorParticle(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                               std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) {
  static constexpr MFloat particleLimit = 1;

  std::function<MFloat(const MInt)> noPart = [&](const MInt cellId) {
    return static_cast<MFloat>(a_noParticlesInCell(cellId)) + static_cast<MFloat>(a_noEllipsoidsInCell(cellId));
  };

  this->sensorLimit(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen, noPart, particleLimit, &m_bandWidth[0],
                    true);
}


/**
 * \brief set the sensors around the boundary positions
 *       to get the markedCells
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::sensorInterface(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                                std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) {
  MIntScratchSpace bandWidth(maxRefinementLevel() + 1, AT_, "bandWidth");
  MInt range1 = 2;
  MInt range2 = 1;
  bandWidth[maxRefinementLevel() - 1] = range1;
  for(MInt i = maxRefinementLevel() - 2; i >= 0; i--) {
    bandWidth[i] = (bandWidth[i + 1] / 2) + 1 + range2;
  }

  static constexpr MFloat limit = 1;
  std::function<MFloat(const MInt)> isBndry = [&](const MInt cellId) {
    return static_cast<MFloat>(a_bndryCellId(cellId));
  };

  this->sensorLimit(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen, isBndry, limit, &bandWidth[0], true,
                    false);
}


/*  \brief post Adaptation call
 *
 *  \date October-2020
 *  \author Sven Berger, updates Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::postAdaptation() {
  TRACE();

  this->compactCells();

  m_freeIndices.clear();

  grid().updateOther();

  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);
  this->checkNoHaloLayers();

  if(!isActive()) return;

  grid().updateLeafCellExchange();
  initMPI(true);

  if(globalTimeStep < 0) {
    findSpawnCellId();
    if(m_spawnCellId > -1) {
      a_noParticlesInCell(m_spawnCellId) += 1;
      if(m_ellipsoids) a_noEllipsoidsInCell(m_spawnCellId) += 1;
    }
  }

  this->exchangeData(&a_noParticlesInCell(0), 1);
  this->exchangeData(&a_bndryCellId(0), 1);
  if(m_ellipsoids) this->exchangeData(&a_noEllipsoidsInCell(0), 1);

  resetSourceTerms(noInternalCells());

  m_adaptationLevel++;
}


/**
 * \brief  reinit the solver after the full adaptation loop!
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::finalizeAdaptation() {
  TRACE();

  if(!isActive()) return;

  if(domainId() == 0) {
    cerr << "Finding cartesian neighbors for LPT-solver!" << endl;
  }
  grid().findCartesianNghbIds();

  updateExchangeCells();

  findSpawnCellId();

  m_cellToNghbrHood.clear();
  m_cellToNghbrHoodInterpolation.clear();

  // update particle cellIds
  // this is probably still faster than having to swap all particle-cellIds during
  // multiple adaptation steps!
  for(MInt i = 0; i < a_noParticles(); i++) {
    m_partList[i].m_neighborList.clear();
    const MInt cellId = grid().findContainingLeafCell(&m_partList[i].m_position[0]);
    ASSERT(cellId > -1, "");
    m_partList[i].m_cellId = cellId;
    m_partList[i].m_oldCellId = cellId;
    ASSERT(m_partList[i].m_cellId > 0, "ERROR: Illegal cellId after adaptation!");
    ASSERT(c_isLeafCell(cellId) && c_level(cellId) == maxRefinementLevel(), "");
  }
  for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
    m_partListEllipsoid[i].m_neighborList.clear();
    const MInt cellId = grid().findContainingLeafCell(&m_partListEllipsoid[i].m_position[0]);
    ASSERT(cellId > -1, "");
    m_partListEllipsoid[i].m_cellId = cellId;
    m_partListEllipsoid[i].m_oldCellId = cellId;
    ASSERT(m_partListEllipsoid[i].m_cellId > 0, "ERROR: Illegal cellId after adaptation!");
    ASSERT(c_isLeafCell(cellId) && c_level(cellId) == maxRefinementLevel(), "");
  }

  countParticlesInCells();
  setRegridTrigger();

  updateFluidFraction();

  checkParticles();
  checkCells();

  // unset bndry-Cell id after adaptation
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_bndryCellId(cellId) = -1;
  }

  // LPT-statistics
  MInt noCells = a_noCells();
  MInt minCells = noCells;
  MInt maxCells = noCells;
  MInt noPart = a_noParticles();
  MInt minParts = noPart;
  MInt maxParts = noPart;
  MInt noEllips = a_noEllipsoidalParticles();
  MInt minEllips = noEllips;
  MInt maxEllips = noEllips;

  MPI_Allreduce(MPI_IN_PLACE, &noCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noCells");
  MPI_Allreduce(MPI_IN_PLACE, &noPart, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noPart");
  MPI_Allreduce(MPI_IN_PLACE, &maxCells, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxCells");
  MPI_Allreduce(MPI_IN_PLACE, &minCells, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minCells");
  MPI_Allreduce(MPI_IN_PLACE, &maxParts, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxParts");
  MPI_Allreduce(MPI_IN_PLACE, &minParts, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minParts");
  if(m_ellipsoids) {
    MPI_Allreduce(MPI_IN_PLACE, &noEllips, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noEllips");
    MPI_Allreduce(MPI_IN_PLACE, &maxEllips, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxEllips");
    MPI_Allreduce(MPI_IN_PLACE, &minEllips, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minEllips");
  }

  if(domainId() == 0) {
    cerr << "LPT cell-statistics : " << minCells << " - " << maxCells << " avg " << noCells / noDomains() << endl;
    cerr << "LPT particle-statistics : " << minParts << " - " << maxParts << " avg " << noPart / noDomains() << endl;
    if(m_ellipsoids)
      cerr << "LPT ellipsoid-statistics : " << minEllips << " - " << maxEllips << " avg " << noEllips / noDomains()
           << endl;
  }
}

template <MInt nDim>
void LPT<nDim>::resizeGridMap() {
  grid().resizeGridMap(m_cells.size());
}

/**
 * \brief
 * \author Lennart Schneiders
 */
template <MInt nDim>
void LPT<nDim>::swapProxy(const MInt cellId0, const MInt cellId1) {
  grid().swapGridIds(cellId0, cellId1);
}


/**
 * \brief refining a cartesian cell
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::refineCell(const MInt gridCellId) {
  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(solverCellId > -1 && solverCellId < m_cells.size(), "solverCellId is: " << solverCellId);

  MInt noNewChildren = 0;

  for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
    const MInt gridChildId = grid().raw().treeb().child(gridCellId, child);
    if(gridChildId == -1) continue;

    // Skip if cell is a partition level ancestor and its child was not newly created
    if(!grid().raw().a_hasProperty(gridChildId, Cell::WasNewlyCreated)
       && grid().raw().a_hasProperty(gridCellId, Cell::IsPartLvlAncestor)) {
      continue;
    }

    // If solver is inactive all cells musst be halo cells!
    if(!isActive()) ASSERT(grid().raw().a_isHalo(gridChildId), "");
    // If child exists in grid but is not located inside solver geometry
    if(!grid().solverFlag(gridChildId, solverId())) continue;

    const MInt solverChildId = this->createCellId(gridChildId);
    noNewChildren++;

    for(MInt n = 0; n < nDim; n++) {
      a_fluidVelocity(solverChildId, n) = 0.0;
    }
    a_fluidDensity(solverChildId) = 0.0;
    a_fluidPressure(solverChildId) = 0.0;
    a_fluidTemperature(solverChildId) = 0.0;

    if(m_evaporation) {
      a_fluidSpecies(solverChildId) = 0.0;
    }
    a_isValidCell(solverChildId) = true;

    a_bndryCellId(solverChildId) = a_bndryCellId(solverCellId);
    if(m_ellipsoids) {
      for(MInt varId = 0; varId < nDim; varId++) {
        for(MInt dir = 0; dir < nDim; dir++) {
          a_velocitySlope(solverChildId, varId, dir) = F0;
        }
      }
    }
  }

  if(noNewChildren == 0) return;

  const MInt noParticles = a_noParticlesInCell(solverCellId);
  const MInt noEllipsoids = a_noEllipsoidsInCell(solverCellId);

  //NOTE: setting an estimate of the number particles in the cell
  //      is necessary as this is used to set the sensors for further refinement!
  //      the true value will be updated in finalizeAdaptation!
  const MInt noParticlesEach = noParticles / noNewChildren;
  const MInt noEllipsoidsEach = noEllipsoids / noNewChildren;

  for(MInt child = 0; child < c_noChildren(solverCellId) -1; child++) {
    const MInt childId = c_childId(solverCellId,child);
    if(childId < 0) continue;
    a_noParticlesInCell(childId) = noParticles > 0 ? noParticlesEach : 0;
    a_noEllipsoidsInCell(childId) = noEllipsoids > 0 ? noEllipsoidsEach : 0;

    // average of solver parent
    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(childId, i) = a_momentumFlux(solverCellId, i) / noNewChildren;
      }
      a_workFlux(childId) = a_workFlux(solverCellId) / noNewChildren;
    }
    if(m_heatCoupling) {
      a_heatFlux(childId) = a_heatFlux(solverCellId) / noNewChildren;
    }
    if(m_massCoupling) {
      a_massFlux(childId) = a_massFlux(solverCellId) / noNewChildren;
    }
  }
  const MInt remainingParticles = noParticles - ((c_noChildren(solverCellId) - 1) * noParticlesEach);
  const MInt remainingEllipsoids = noEllipsoids - ((c_noChildren(solverCellId) - 1) * noEllipsoidsEach);

  const MInt remainingChild = c_childId(solverCellId, c_noChildren(solverCellId) - 1);
  if(remainingChild > 0) {
    a_noParticlesInCell(remainingChild) = noParticles > 0 ? remainingParticles : 0;
    a_noEllipsoidsInCell(remainingChild) = noEllipsoids > 0 ? remainingEllipsoids : 0;

    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(remainingChild, i) = a_momentumFlux(solverCellId, i) / noNewChildren;
      }
      a_workFlux(remainingChild) = a_workFlux(solverCellId) / noNewChildren;
    }
    if(m_heatCoupling) {
      a_heatFlux(remainingChild) = a_heatFlux(solverCellId) / noNewChildren;
    }
    if(m_massCoupling) {
      a_massFlux(remainingChild) = a_massFlux(solverCellId) / noNewChildren;
    }
  }

  // reset values on new non-leaf cell
  if(noNewChildren > 0 ){
    a_noParticlesInCell(solverCellId) = 0;
    a_noEllipsoidsInCell(solverCellId) = 0;

    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(solverCellId, i) = 0;
      }
      a_workFlux(solverCellId) = 0;
    }
    if(m_heatCoupling) {
      a_heatFlux(solverCellId) = 0;
    }
    if(m_massCoupling) {
      a_massFlux(solverCellId) = 0;
    }
  }
}

/**
 * \brief coarseining a cartesian cell
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::removeChilds(const MInt gridCellId) {
  // If solver is inactive cell musst never be a internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(solverCellId > -1 && solverCellId < m_cells.size(), "solverCellId is: " << solverCellId);

  // reset values of new leaf cell
  a_noParticlesInCell(solverCellId) = 0;
  a_noEllipsoidsInCell(solverCellId) = 0;
  if(m_momentumCoupling) {
    for(MInt i = 0; i < nDim; i++) {
      a_momentumFlux(solverCellId, i) = 0;
    }
    a_workFlux(solverCellId) = 0;
  }
  if(m_heatCoupling) {
    a_heatFlux(solverCellId) = 0;
  }
  if(m_massCoupling) {
    a_massFlux(solverCellId) = 0;
  }

  for(MInt c = 0; c < grid().m_maxNoChilds; c++) {
    const MInt childId = c_childId(solverCellId, c);
    if(childId < 0) continue;

    // sum values of childs which will be removed
    a_noParticlesInCell(solverCellId) += a_noParticlesInCell(childId);
    a_noEllipsoidsInCell(solverCellId) += a_noEllipsoidsInCell(childId);
    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(solverCellId, i) += a_momentumFlux(childId, i);
      }
      a_workFlux(solverCellId) += a_workFlux(childId);
    }
    if(m_heatCoupling) {
      a_heatFlux(solverCellId) += a_heatFlux(childId);
    }
    if(m_massCoupling) {
      a_massFlux(solverCellId) += a_massFlux(childId);
    }

    this->removeCellId(childId);
  }
}

/**
 * \brief
 */
template <MInt nDim>
void LPT<nDim>::removeCell(const MInt gridCellId) {
  // If solver is inactive cell musst never be a internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(gridCellId > -1 && gridCellId < grid().raw().treeb().size() && solverCellId > -1
             && solverCellId < m_cells.size() && grid().tree().solver2grid(solverCellId) == gridCellId,
         "");
  ASSERT(c_noChildren(solverCellId) == 0, "");
  this->removeCellId(solverCellId);
}


/**
 * \brief  checks if a child lies outSide of the domain!
 *         necessary for refinement at the bndry!
 * \author Tim Wegmann
 */

template <MInt nDim>
MInt LPT<nDim>::cellOutside(const MFloat* coords, const MInt level, const MInt gridCellId) {
  // currently all cells with the refineFlag are keep during adaptation!
  return -1;

  std::ignore = gridCellId;
  std::ignore = coords[0];
  std::ignore = level;

  /*
  static constexpr MInt cornerIndices[8][3] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
                                               {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
  static constexpr const MInt noCorners = (nDim == 2) ? 4 : 8;

  MFloat corner[3] = {0, 0, 0};
  MBool outside = true;
  MFloat cellHalfLength = F1B2 * c_cellLengthAtLevel(level);


  for(MInt i = 0; i < noCorners; i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      corner[dim] = coords[dim] + cornerIndices[i][dim] * cellHalfLength;
    }
    IF_CONSTEXPR(nDim == 2) {
      if(!m_geometry->pointIsInside(corner)) outside = false;
      // pointIsInside == true if Point is outside fluid domain
    } else {
      if(!m_geometry->pointIsInside2(corner)) outside = false;
      // pointIsInside == true if Point is outside fluid domain
    }
  }

  return outside;
  */
}


/**
 * \brief swap Cell information when the grid collector is restructured
 *        during adaptation cleanup
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::swapCells(const MInt cellId0, const MInt cellId1) {
  std::swap(m_cells.properties(cellId1), m_cells.properties(cellId0));

  std::swap(m_cells.volumeFraction(cellId1), m_cells.volumeFraction(cellId0));
  std::swap(m_cells.noParticles(cellId1), m_cells.noParticles(cellId0));
  std::swap(m_cells.noEllipsoids(cellId1), m_cells.noEllipsoids(cellId0));

  for(MInt n = 0; n < PV.noVars(); n++) {
    std::swap(m_cells.fluidVariable(cellId1, n), m_cells.fluidVariable(cellId0, n));
  }

  if(m_evaporation) {
    std::swap(m_cells.fluidSpecies(cellId1), m_cells.fluidSpecies(cellId0));
  }

  if(m_massCoupling) {
    std::swap(m_cells.massFlux(cellId1), m_cells.massFlux(cellId0));
  }
  if(m_heatCoupling) {
    std::swap(m_cells.heatFlux(cellId1), m_cells.heatFlux(cellId0));
  }
  if(m_momentumCoupling) {
    for(MInt i = 0; i < nDim; i++) {
      std::swap(m_cells.momentumFlux(cellId1, i), m_cells.momentumFlux(cellId0, i));
    }
    std::swap(m_cells.workFlux(cellId1), m_cells.workFlux(cellId0));
  }
  std::swap(m_cells.bndryCellId(cellId1), m_cells.bndryCellId(cellId0));

  // only ellipsoids require velocity slopes
  if(m_ellipsoids) {
    for(MInt varId = 0; varId < nDim; varId++) {
      for(MInt dim = 0; dim < nDim; dim++) {
        std::swap(m_cells.velocitySlope(cellId1, varId, dim), m_cells.velocitySlope(cellId0, varId, dim));
      }
    }
  }
}


/**
 * \brief set cell weight for partitioning
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  ASSERT(isActive(), "solver is not active");

  countParticlesInCells();

  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt offset = noCellsGrid * solverId();

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    const MInt id = gridCellId + offset;
    solverCellWeight[id] = m_weightBaseCell * m_weightMulitSolverFactor;
    if(!c_isLeafCell(cellId)) continue;
    solverCellWeight[id] = m_weightLeafCell * m_weightMulitSolverFactor;
    if(a_noParticlesInCell(cellId) == 0 && a_noEllipsoidsInCell(cellId) == 0) continue;
    solverCellWeight[id] =
        m_weightMulitSolverFactor
        * (m_weightParticleCell + (a_noParticlesInCell(cellId) + a_noEllipsoidsInCell(cellId)) * m_weightParticle);
  }

  // additionally increase weight at spawncellId
  if(m_engineSetup && m_sprayModel != nullptr && m_spawnCellId > -1) {
    const MInt gridCellId = grid().tree().solver2grid(m_spawnCellId);
    const MInt id = gridCellId + offset;
    solverCellWeight[id] += m_weightSpawnCell * m_weightMulitSolverFactor;
  }
}

/**
 * \brief  Get solver timings
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings,
                                 const MBool NotUsed(allTimings)) {
  TRACE();

  const MString namePrefix = "b" + std::to_string(solverId()) + "_";

  const MFloat load = returnLoadRecord();
  const MFloat idle = returnIdleRecord();

  solverTimings.emplace_back(namePrefix + "loadLptSolver", load);
  solverTimings.emplace_back(namePrefix + "idleLptSolver", idle);

#ifdef MAIA_TIMER_FUNCTION
  solverTimings.emplace_back(namePrefix + "timeIntegration", RETURN_TIMER_TIME(m_timers[Timers::TimeInt]));
  solverTimings.emplace_back(namePrefix + "motion", RETURN_TIMER_TIME(m_timers[Timers::Motion]));
  solverTimings.emplace_back(namePrefix + "energy", RETURN_TIMER_TIME(m_timers[Timers::Energy]));
  solverTimings.emplace_back(namePrefix + "injection", RETURN_TIMER_TIME(m_timers[Timers::Injection]));
  solverTimings.emplace_back(namePrefix + "secBreakUp", RETURN_TIMER_TIME(m_timers[Timers::Breakup]));
  solverTimings.emplace_back(namePrefix + "wall", RETURN_TIMER_TIME(m_timers[Timers::Wall]));
  solverTimings.emplace_back(namePrefix + "source", RETURN_TIMER_TIME(m_timers[Timers::SourceTerms]));
  solverTimings.emplace_back(namePrefix + "regridExc", RETURN_TIMER_TIME(m_timers[Timers::Exchange2]));
  solverTimings.emplace_back(namePrefix + "exchange", RETURN_TIMER_TIME(m_timers[Timers::Exchange]));
#endif
}

/// \brief Return the default weights for all load quantities
template <MInt nDim>
void LPT<nDim>::getDefaultWeights(MFloat* weights, std::vector<MString>& names) const {
  TRACE();

  weights[0] = 0.1;
  names[0] = "lpt_cell";
  MInt count = 1;

  weights[1] = 0.9;
  names[1] = "lpt_particleCell";
  count++;

  weights[2] = 0.1;
  names[2] = "lpt_particle";
  count++;

  if(m_weightSourceCells) {
    weights[3] = 0.01;
    names[3] = "lpt_source";
    count++;
  }

  if(noLoadTypes() != count) {
    TERMM(1, "Count does not match noLoadTypes.");
  }
}

/**
 * \brief  Limit weight of base cell
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::limitWeights(MFloat* weights) {
  if(m_limitWeights) {
    weights[0] = mMax(weights[0], 0.01 * weights[1]);
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
void LPT<nDim>::getLoadQuantities(MInt* const loadQuantities) const {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // reset
  for(MInt type = 0; type < noLoadTypes(); type++) {
    loadQuantities[type] = 0;
  }

  MInt noLeafCells = 0;
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    if(c_isLeafCell(cellId) && a_isValidCell(cellId)) {
      noLeafCells++;
    }
  }

  loadQuantities[0] = noLeafCells;

  MInt noCellsWithParticles = 0;
  MInt noCellsWithSources = 0;

  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    if(a_noParticlesInCell(cellId) == 0 && a_noEllipsoidsInCell(cellId) == 0) {
      if(m_massCoupling && a_massFlux(cellId) < 0.0) {
        noCellsWithSources++;
      }
      continue;
    }
    noCellsWithParticles++;
  }
  loadQuantities[1] = noCellsWithParticles;
  loadQuantities[2] = a_noSphericalParticles() + a_noEllipsoidalParticles();


  if(m_weightSourceCells) {
    loadQuantities[3] = noCellsWithSources;
  }
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
MFloat LPT<nDim>::getCellLoad(const MInt gridCellId, const MFloat* const weights) const {
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
  MFloat cellLoad = 0.0;
  if(c_isLeafCell(cellId) && a_isValidCell(cellId)) {
    cellLoad = weights[0];
  }

  if(a_noParticlesInCell(cellId) > 0 || a_noEllipsoidsInCell(cellId) > 0) {
    cellLoad = weights[1];
    cellLoad += a_noParticlesInCell(cellId) * weights[2];
  } else if(m_weightSourceCells) {
    if(m_massCoupling && a_massFlux(cellId) < 0.0) {
      cellLoad = weights[3];
    }
  }

  if(m_sprayModel != nullptr && cellId == m_spawnCellId) {
    cellLoad += 5 * weights[1];
  }


  return cellLoad;
}


/**
 * \brief  Return decomposition information, i.e. number of local elements,...
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) {
  TRACE();

  const MString namePrefix = "b" + std::to_string(solverId()) + "_";

  const MInt noInternalCells = isActive() ? grid().noInternalCells() : 0;

  domainInfo.emplace_back(namePrefix + "noLptCells", noInternalCells);

  // Number of Ls-Band-Cells
  MInt noCellsWithParticles = 0;
  MInt noLeafCells = 0;
  MInt noParticles = 0;

  if(isActive()) {
    noParticles = a_noParticles();
    for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
      if(a_noParticlesInCell(cellId) > 0 || a_noEllipsoidsInCell(cellId) > 0) {
        noCellsWithParticles++;
      }
      if(c_isLeafCell(cellId) && a_isValidCell(cellId)) {
        noLeafCells++;
      }
    }
  }

  domainInfo.emplace_back(namePrefix + "noLptLeafCells", noLeafCells);
  domainInfo.emplace_back(namePrefix + "noLptParticleCells", noCellsWithParticles);
  domainInfo.emplace_back(namePrefix + "noLptParticles", noParticles);
}


/**
 *  \brief prepare the LPT timeStep
 *  \date November-20
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::preTimeStep() {
  TRACE();

  RECORD_TIMER_START(m_timers[Timers::PreTime]);

  // Exchange flow field wich has been set by one or more flow solver via the coupler(s)
  // exchangeData(&a_fluidVariable(0, 0), nDim+2);

  // m_partCellMapCurrent = false;

  // advanceTimeStep
  m_time += m_timeStep;
  ASSERT(m_timeStep > MFloatEps, "TS not set!");


  m_timeStepOld = m_timeStep;
  m_timeStepUpdated = false;
  m_sumEvapMass = 0;
  m_noSendParticles = 0;

  if(m_skipLPT) {
    return;
  }

  resetSourceTerms(0);

  // reset creation time of all particles
  for(MInt i = 0; i < a_noParticles(); i++) {
    m_partList[i].m_creationTime = 0.0;
  }
  for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
    m_partListEllipsoid[i].m_creationTime = 0.0;
  }

  RECORD_TIMER_STOP(m_timers[Timers::PreTime]);
}

/**
 *  \brief reset all source terms with cellId larger than the offset
 *  \date November-20
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::resetSourceTerms(const MInt offset) {
  TRACE();

  for(MInt cellId = offset; cellId < a_noCells(); cellId++) {
    if(m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        a_momentumFlux(cellId, i) = 0;
      }
      a_workFlux(cellId) = 0;
    }
    if(m_heatCoupling) {
      a_heatFlux(cellId) = 0;
    }
    if(m_massCoupling) {
      a_massFlux(cellId) = 0;
    }
  }
}

/**
 *  \brief Time stepping of particles
 *         NOTE: the solutionStep is split into sub-steps, to allow for an interleafed computation
 *               when coupled with a different solver!
 *  \date December-09
 *  \author Rudie Kunnen, Sven Berger, Tim Wegmann
 */
template <MInt nDim>
MBool LPT<nDim>::solutionStep() {
  TRACE();

  if((m_partList.empty() || m_partListEllipsoid.empty()) && !m_spawnParticles && !m_activePrimaryBUp
     && globalTimeStep == 0 && noDomains() == 1) {
    return true;
  }

  // To set particle release Time for particles in isoTropic Turbulence
  // Release Timestep is determined by a Skewness of 0.44 of fluid velocity
  if(m_skewnessTimeStep != 0 && globalTimeStep < m_skewnessTimeStep) return true;

  RECORD_TIMER_START(m_timers[Timers::TimeInt]);

  // apply arbitrary LPT load to test&validate interleaved coupling
  // remaining for documentation and testing purposes
  if(m_sleepLPT > 0.0) {
    cerr0 << "Sleeping LPT-solver for " << m_sleepLPT << " milli-seconds!" << endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(m_sleepLPT));
  }

  if(m_skipLPT) {
    cerr0 << "Skipping LPT execution!" << endl;
    return true;
  }


  MBool completed = false;
  switch(m_noSolutionSteps) {
    case 1: {
      completed = oneStageSolutionStep();
      break;
    }
    case 5: {
      completed = fiveStageSolutionStep();
      break;
    }
    default: {
      mTerm(1, AT_, "SolutionStep method not implemented in LPT!");
      break;
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::TimeInt]);
  return completed;
}


/**
 *  \brief Splitting the solutionStep into 5 separate parts to match
 *         the 5-stage Runge-Kutta time-stepping for inter-leafed computation!
 *  \date August 21
 *  \author Tim Wegmann
 */
template <MInt nDim>
MBool LPT<nDim>::fiveStageSolutionStep() {
  NEW_TIMER_GROUP_STATIC(t_part, "particle timestep");
  NEW_TIMER_STATIC(partTS, "particle timestep", t_part);
  NEW_SUB_TIMER_STATIC(spawn, "spawn", partTS);
  NEW_SUB_TIMER_STATIC(motion, "motion", partTS);
  NEW_SUB_TIMER_STATIC(transfer, "transfer", partTS);
  NEW_SUB_TIMER_STATIC(evap, "evaporation", partTS);
  NEW_SUB_TIMER_STATIC(couplingTerms, "couplingTerms", partTS);
  NEW_SUB_TIMER_STATIC(collision1, "wall-collision", partTS);
  NEW_SUB_TIMER_STATIC(primBU, "primaryBreakUp", partTS);
  NEW_SUB_TIMER_STATIC(secBU, "secondaryBreakUp", partTS);
  NEW_SUB_TIMER_STATIC(respawn, "particle-respawn", partTS);
  NEW_SUB_TIMER_STATIC(deleteP, "Remove-particles", partTS);

  RECORD_TIMER_START(partTS);

  m_solutionStep++;

  switch(m_solutionStep) {
    case 1: {
      // blocking injection step
      RECORD_TIMER_START(primBU);
      // some debug checks
      checkParticles();
      checkCells();

      if(m_activePrimaryBUp) {
        sprayInjection();
      }
      RECORD_TIMER_STOP(primBU);

      RECORD_TIMER_STOP(partTS);
      return false;
      break;
    }
    case 2: {
      // receive flow field data from previous TS
      receiveFlowField();
      // if ellipsoids are activated also receive velocity slopes from previous TS
      receiveVelocitySlopes();

      if(m_nonBlockingComm && m_activePrimaryBUp) {
        RECORD_TIMER_START(transfer);
        receiveParticles<true>();
        m_primaryExchange = false;
        RECORD_TIMER_STOP(transfer);
      }

      // fully non-blocking movement of particles
      RECORD_TIMER_START(motion);
      motionEquation(0);
      RECORD_TIMER_STOP(motion);

      RECORD_TIMER_START(collision1);
      wallCollision();
      RECORD_TIMER_STOP(collision1);

      if(m_nonBlockingComm) {
        RECORD_TIMER_START(transfer);
        exchangeParticles(false, 0);
        RECORD_TIMER_STOP(transfer);
      }

      RECORD_TIMER_STOP(partTS);
      return false;
      break;
    }
    case 3: {
      if(!m_nonBlockingComm) {
        // blocking exchange step
        RECORD_TIMER_START(transfer);
        exchangeParticles(false, 0);
        RECORD_TIMER_STOP(transfer);
      } else {
        // non-blocking receive
        RECORD_TIMER_START(transfer);
        receiveParticles<false>();
        RECORD_TIMER_STOP(transfer);
      }

      // non-blocking eveporation
      RECORD_TIMER_START(evap);
      evaporation(0);
      RECORD_TIMER_STOP(evap);

      if(m_nonBlockingComm) {
        RECORD_TIMER_START(couplingTerms);
        coupling(0);
        RECORD_TIMER_STOP(couplingTerms);
        sendSourceTerms();

        if(this->m_adaptation && this->m_noSensors > 0) {
          checkRegridTrigger();
        }
      }

      RECORD_TIMER_STOP(partTS);

      return false;
      break;
    }
    case 4: {
      if(!m_nonBlockingComm) {
        // non-blocking computation of source-terms
        RECORD_TIMER_START(couplingTerms);
        coupling(0);
        RECORD_TIMER_STOP(couplingTerms);
      } else {
        // wait and receive source terms
        RECORD_TIMER_START(transfer);
        receiveSourceTerms();
        RECORD_TIMER_STOP(transfer);
      }

      RECORD_TIMER_STOP(partTS);
      return false;
      break;
    }
    case 5: {
      // non-blocking step unless using outdated particle spawn/respawn or particle-collision
      // only particle data-structure related
      RECORD_TIMER_START(secBU);
      RECORD_TIMER_START(m_timers[Timers::Breakup]);
      if(m_activeSecondaryBUp) {
        m_sprayModel->secondaryBreakUp(m_timeStep);
      }
      RECORD_TIMER_STOP(m_timers[Timers::Breakup]);
      RECORD_TIMER_STOP(secBU);

      // check for collided particles to be copied to other ranks
      RECORD_TIMER_START(transfer);
      if(m_collisions < 5 && m_collisions > 0) {
        exchangeParticles(true, 0);
      }
      RECORD_TIMER_STOP(transfer);

      RECORD_TIMER_START(spawn);
      if(m_spawnParticles) {
        spawnTimeStep();
      }
      RECORD_TIMER_STOP(spawn);

      RECORD_TIMER_START(deleteP);
      removeInvalidParticles(false);
      RECORD_TIMER_STOP(deleteP);

      RECORD_TIMER_START(respawn);
      if(m_respawn) {
        particleRespawn();
      }
      RECORD_TIMER_STOP(respawn);

      RECORD_TIMER_START(deleteP);
      updateFluidFraction();

      if(a_timeStepComputationInterval() > 0 && globalTimeStep % a_timeStepComputationInterval() == 0) {
        forceTimeStep(computeTimeStep());
      }

      // gather post-processing data
      if(m_sprayModel != nullptr) {
        const MInt count = m_sprayModel->m_injDataSize;
        vector<MFloat> tmp(count + m_addInjDataCnt);
        for(MInt i = 0; i < count; i++) {
          tmp[i] = m_sprayModel->m_injData[i];
        }
        tmp[count] = m_sumEvapMass;
        tmp[count + 1] = m_sprayModel->m_noRTsecBreakUp;
        tmp[count + 2] = m_sprayModel->m_noKHsecBreakUp;
        tmp[count + 3] = m_noSendParticles;

        m_injData.insert(make_pair(globalTimeStep, tmp));
      }

      RECORD_TIMER_STOP(deleteP);
      RECORD_TIMER_STOP(partTS);

      return true;
      break;
    }
    default: {
      mTerm(1, AT_, "Invalid solution Step!");
      break;
    }
  }
}

/**
 *  \brief single solution step for pure LPT computation or without interleafed coupling!
 *  \date August 21
 *  \author Tim Wegmann
 */
template <MInt nDim>
MBool LPT<nDim>::oneStageSolutionStep() {
  NEW_TIMER_GROUP_STATIC(t_part, "particle timestep");
  NEW_TIMER_STATIC(partTS, "particle timestep", t_part);
  NEW_SUB_TIMER_STATIC(spawn, "spawn", partTS);
  NEW_SUB_TIMER_STATIC(motion, "motion", partTS);
  NEW_SUB_TIMER_STATIC(transfer, "transfer", partTS);
  NEW_SUB_TIMER_STATIC(evap, "evaporation", partTS);
  NEW_SUB_TIMER_STATIC(couplingTerms, "couplingTerms", partTS);
  NEW_SUB_TIMER_STATIC(collision1, "wall-collision", partTS);
  NEW_SUB_TIMER_STATIC(primBU, "primaryBreakUp", partTS);
  NEW_SUB_TIMER_STATIC(secBU, "secondaryBreakUp", partTS);
  NEW_SUB_TIMER_STATIC(respawn, "particle-respawn", partTS);
  NEW_SUB_TIMER_STATIC(deleteP, "Remove-particles", partTS);

  RECORD_TIMER_START(partTS);

  ASSERT(!m_nonBlockingComm, "");

  // some debug checks
  checkParticles();
  checkCells();

  RECORD_TIMER_START(primBU);
  if(m_activePrimaryBUp) {
    sprayInjection();
  }
  RECORD_TIMER_STOP(primBU);

  // receive flow field data from previous TS
  receiveFlowField();
  // if ellipsoids are activated also receive velocity slopes from previous TS
  receiveVelocitySlopes();

  // move particles
  RECORD_TIMER_START(motion);
  motionEquation(0);
  RECORD_TIMER_STOP(motion);

  RECORD_TIMER_START(collision1);
  wallCollision();
  RECORD_TIMER_STOP(collision1);

  RECORD_TIMER_START(transfer);
  exchangeParticles(false, 0);
  RECORD_TIMER_STOP(transfer);

  RECORD_TIMER_START(evap);
  evaporation(0);
  RECORD_TIMER_STOP(evap);

  RECORD_TIMER_START(couplingTerms);
  coupling(0);
  RECORD_TIMER_STOP(couplingTerms);

  RECORD_TIMER_START(secBU);
  RECORD_TIMER_START(m_timers[Timers::Breakup]);
  if(m_activeSecondaryBUp) {
    m_sprayModel->secondaryBreakUp(m_timeStep);
  }
  RECORD_TIMER_STOP(m_timers[Timers::Breakup]);
  RECORD_TIMER_STOP(secBU);

  // check for collided particles to be copied to other ranks
  RECORD_TIMER_START(transfer);
  if(m_collisions < 5 && m_collisions > 0) {
    exchangeParticles(true, 0);
  }
  RECORD_TIMER_STOP(transfer);

  RECORD_TIMER_START(spawn);
  if(m_spawnParticles) {
    spawnTimeStep();
  }
  RECORD_TIMER_STOP(spawn);

  RECORD_TIMER_START(deleteP);
  removeInvalidParticles(false);
  RECORD_TIMER_STOP(deleteP);

  // reduceParticles(); // TODO labels:LPT activate and make settable

  // advanceParticles(0);

  RECORD_TIMER_START(respawn);
  if(m_respawn) {
    particleRespawn();
  }
  RECORD_TIMER_STOP(respawn);

  updateFluidFraction();

  if(this->m_adaptation && this->m_noSensors > 0) {
    checkRegridTrigger();
  }

  if(a_timeStepComputationInterval() > 0 && globalTimeStep % a_timeStepComputationInterval() == 0) {
    forceTimeStep(computeTimeStep());
  }

  // gather post-processing data
  if(m_sprayModel != nullptr) {
    const MInt count = m_sprayModel->m_injDataSize;
    vector<MFloat> tmp(count + m_addInjDataCnt);
    for(MInt i = 0; i < count; i++) {
      tmp[i] = m_sprayModel->m_injData[i];
    }
    tmp[count] = m_sumEvapMass;
    tmp[count + 1] = m_sprayModel->m_noRTsecBreakUp;
    tmp[count + 2] = m_sprayModel->m_noKHsecBreakUp;
    tmp[count + 3] = m_noSendParticles;

    m_injData.insert(make_pair(globalTimeStep, tmp));
  }

  RECORD_TIMER_STOP(partTS);
  return true;
}

/**
 *  \brief after the LPT time-Step!
 *
 *  \date November-20
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::postTimeStep() {
  TRACE();

  NEW_TIMER_GROUP_STATIC(t_post, "LPT postTimeStep");
  NEW_TIMER_STATIC(postTS, "postTimeStep", t_post);

  RECORD_TIMER_START(postTS);

  if(m_skipLPT) {
    return;
  }
  if(m_sleepLPT > 0.0) {
    std::this_thread::sleep_for(std::chrono::milliseconds(m_sleepLPT));
  }

  // some debug checks
  checkParticles();
  checkCells();

  m_solutionStep = 0;


  removeInvalidParticles(true);
  advanceParticles(0);

  RECORD_TIMER_STOP(postTS);
}

/**
 *  \brief before writing a restart file
 *
 *  \date November-20
 *  \author Tim Wegmann
 */
template <MInt nDim>
MBool LPT<nDim>::prepareRestart(MBool force, MBool& writeGridRestart) {
  MBool writeRestart = false;

  if((((globalTimeStep % this->m_restartInterval) == 0 && (globalTimeStep > m_restartTimeStep)) || force)) {
    if(isActive()) {
      countParticlesInCells();
      if(!m_nonBlockingComm) {
        exchangeSourceTerms();
      } else {
        sendSourceTerms();
        receiveSourceTerms();
        receiveFlowField();
        receiveVelocitySlopes();
        waitForSendReqs();
      }
    }
    writeRestart = true;
  }

  if(m_restart && !m_restartFile) writeGridRestart = true;

  return (force || writeRestart);
}

/**
 *  \brief writing a restart File as regularly called from the grid-controller!
 *
 *  \date January-2021
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::writeRestartFile(const MBool writeRestart, const MBool, const MString gridFileName,
                                 MInt* recalcIdTree) {
  if(writeRestart) {
    checkParticles();

    writeCellSolutionFile(gridFileName, &recalcIdTree[0]);

    // write the particle based restart-file
    writeParticleRestartFile();
  }
}

/**
 *  \brief writing a restart File
 *  NOTE: for debugging purposes only!
 *  \date January-2021
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::saveDebugRestartFile() {
  writeParticleRestartFile();

  MBool writeCellFile = true;
  if(grid().hasInactiveRanks()) writeCellFile = false;

  // write an additional cell based solution file (useful for debuging!)
  if(writeCellFile) {
    stringstream gridFile;
    stringstream gridFilePath;
    gridFile << "grid_LPTDebug_" << globalTimeStep << ".Netcdf";
    gridFilePath << outputDir() << "grid_LPTDebug_" << globalTimeStep << ".Netcdf";
    MIntScratchSpace recalcIds(grid().raw().treeb().size(), AT_, "recalcIds");

    grid().raw().saveGrid((gridFilePath.str()).c_str(), recalcIds.begin());

    writeCellSolutionFile((gridFile.str()).c_str(), &recalcIds[0]);
  }
}

/**
 *  \brief spawn new particles and compute their timeStep
 *  NOTE: this is an older functionality which is only used in testcases...
 *  \date August-2021
 *  \author ???, update Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::spawnTimeStep() {
  // current number of particles
  const MInt oldNum = a_noParticles();

  ASSERT(!m_nonBlockingComm, "");

  spawnParticles();

  exchangeParticles(false, oldNum, false);

  motionEquation(oldNum);

  exchangeParticles(false, oldNum);

  evaporation(oldNum);

  coupling(oldNum);
}

template <MInt nDim>
void LPT<nDim>::sprayInjection() {
  RECORD_TIMER_START(m_timers[Timers::Injection]);
  const MInt oldNum = a_noParticles();

  m_sprayModel->injection(m_timeStep);

  // check for particles to be copied to other solvers
  // this is necessary for particles during injection
  // for which a leaf-halo cell could be found
  // where the injector dimension lies within a 2 minLevel cells
  // meaning that at least lower level halo cells can be found and an exchange is possible
  // instead of a broadcast!
  MBool allowNonLeaf = !m_sprayModel->m_broadcastInjected;
  m_primaryExchange = true;
  exchangeParticles(false, oldNum, allowNonLeaf);
  if(!m_nonBlockingComm) {
    m_primaryExchange = false;
  }
  RECORD_TIMER_STOP(m_timers[Timers::Injection]);
}


template <MInt nDim>
void LPT<nDim>::removeInvalidParticles(const MBool alsoFullyEvaporated) {
  // check list for inactive particles and remove these
  // inactive particles are generally particles which:
  // a) have been communicated to a different rank
  // b) have been fully evaporated

  if(!m_respawn) {
    for(MInt i = a_noParticles() - 1; i >= 0; i--) {
      if(!alsoFullyEvaporated && m_partList[i].fullyEvaporated()) continue;
      if(m_partList[i].isInvalid()) {
        if(!m_partList[i].fullyEvaporated() && !m_partList[i].wasSend()) {
          if(m_partList[i].m_cellId < 0 || m_partList[i].m_cellId >= a_noCells() || !a_isHalo(m_partList[i].m_cellId)) {
            cerr << "Removing at " << m_partList[i].m_position[0] << " " << m_partList[i].m_position[1] << " "
                 << m_partList[i].m_position[nDim - 1] << " " << m_partList[i].m_cellId << " "
                 << m_partList[i].m_oldCellId << " " << m_partList[i].hadWallColl() << " " << m_partList[i].m_oldPos[0]
                 << " " << m_partList[i].m_oldPos[1] << " " << m_partList[i].m_oldPos[nDim - 1] << " "
                 << m_partList[i].m_partId << endl;
            if(m_partList[i].m_cellId > -1 && m_partList[i].m_cellId < a_noCells()) {
              cerr << "Cell-stats " << a_isHalo(m_partList[i].m_cellId) << " " << a_isBndryCell(m_partList[i].m_cellId)
                   << " " << a_isValidCell(m_partList[i].m_cellId) << endl;
              if(m_partList[i].m_oldCellId > -1 && m_partList[i].m_oldCellId < a_noCells()) {
                cerr << "Old-cell-stats " << a_isBndryCell(m_partList[i].m_oldCellId) << " "
                     << a_isValidCell(m_partList[i].m_oldCellId) << endl;
              }
              m_partList[i].checkCellChange(nullptr, false);
              cerr << "After-update: " << m_partList[i].m_cellId << " " << m_partList[i].isInvalid() << " "
                   << a_isValidCell(m_partList[i].m_cellId) << endl;
            }
          }
        }
        m_partList.erase(m_partList.begin() + i);
      }
    }
    if(m_ellipsoids) {
      for(MInt i = a_noEllipsoidalParticles() - 1; i >= 0; i--) {
        if(m_partListEllipsoid[i].isInvalid()) {
          if(!m_partListEllipsoid[i].fullyEvaporated() && !m_partListEllipsoid[i].wasSend()) {
            if(m_partListEllipsoid[i].m_cellId < 0 || m_partListEllipsoid[i].m_cellId >= a_noCells()
               || !a_isHalo(m_partListEllipsoid[i].m_cellId)) {
              cerr << "Removing at " << m_partListEllipsoid[i].m_position[0] << " "
                   << m_partListEllipsoid[i].m_position[1] << " " << m_partListEllipsoid[i].m_position[nDim - 1] << " "
                   << m_partListEllipsoid[i].m_cellId << " " << m_partListEllipsoid[i].m_oldCellId << " "
                   << m_partListEllipsoid[i].hadWallColl() << " " << m_partListEllipsoid[i].m_oldPos[0] << " "
                   << m_partListEllipsoid[i].m_oldPos[1] << " " << m_partListEllipsoid[i].m_oldPos[nDim - 1] << " "
                   << m_partListEllipsoid[i].m_partId << endl;
              if(m_partListEllipsoid[i].m_cellId > -1 && m_partListEllipsoid[i].m_cellId < a_noCells()) {
                cerr << "Cell-stats " << a_isHalo(m_partListEllipsoid[i].m_cellId) << " "
                     << a_isBndryCell(m_partListEllipsoid[i].m_cellId) << " "
                     << a_isValidCell(m_partListEllipsoid[i].m_cellId) << endl;
                if(m_partListEllipsoid[i].m_oldCellId > -1 && m_partListEllipsoid[i].m_oldCellId < a_noCells()) {
                  cerr << "Old-cell-stats " << a_isBndryCell(m_partListEllipsoid[i].m_oldCellId) << " "
                       << a_isValidCell(m_partListEllipsoid[i].m_oldCellId) << endl;
                }
                m_partListEllipsoid[i].checkCellChange(nullptr, false);
                cerr << "After-update: " << m_partListEllipsoid[i].m_cellId << " " << m_partListEllipsoid[i].isInvalid()
                     << " " << a_isValidCell(m_partListEllipsoid[i].m_cellId) << endl;
              }
            }
          }
          m_partListEllipsoid.erase(m_partListEllipsoid.begin() + i);
        }
      }
    }
  }
}


/**
 *  \brief set accessor for a_noParticlesInCell
 *         which is used for the particle-sensor during adaptation
 *  \date Aug. 2021
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::countParticlesInCells() {
  TRACE();

  // reset
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_noParticlesInCell(cellId) = 0;
    a_noEllipsoidsInCell(cellId) = 0;
  }

  // fast, local count
  for(MInt i = 0; i < a_noParticles(); i++) {
    const MInt cellId = m_partList[i].m_cellId;
    if(m_partList[i].isInvalid()) {
      mTerm(1, AT_, "Counting invalid particle!");
    }
    ASSERT(cellId > -1, "");
    ASSERT(!a_isHalo(cellId), "Particle in halo-cell!");
    ASSERT(c_isLeafCell(cellId), "Particle in non-leaf cell!");
    a_noParticlesInCell(cellId)++;
  }
  // Ellipsoidal particles
  for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
    const MInt cellId = m_partListEllipsoid[i].m_cellId;
    if(m_partListEllipsoid[i].isInvalid()) {
      mTerm(1, AT_, "Counting invalid ellipsoid!");
    }
    ASSERT(cellId > -1, "");
    ASSERT(!a_isHalo(cellId), "Ellipsoid in halo-cell!");
    if(globalTimeStep > 0) ASSERT(c_isLeafCell(cellId), "Ellipsoid in non-leaf cell!");
    a_noEllipsoidsInCell(cellId)++;
  }
}

/**
 *  \brief set cell-particle mapping
 *         which is used to merge particle in reduceParticles!
 *  \author Sven Berger
 */
template <MInt nDim>
void LPT<nDim>::perCellStats() {
  TRACE();

  // sort particles after id to make results consistent and faster searching!
  // own_sort(m_partList, sort_particleAfterPartIds<LPTSpherical<nDim>>::compare);

  // update mapping from cellId to particles
  m_cellToPartMap.clear();
  // generate mapping cell to particles
  for(MInt i = 0; i < a_noParticles(); i++) {
    m_cellToPartMap.insert({m_partList[i].m_cellId, &m_partList[i]});
    ASSERT(!a_isHalo(m_partList[i].m_cellId), "Particle in halo-cell!");
  }
  if(m_ellipsoids) {
    m_cellToEllipsMap.clear();
    for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
      m_cellToEllipsMap.insert({m_partListEllipsoid[i].m_cellId, &m_partListEllipsoid[i]});
      ASSERT(!a_isHalo(m_partListEllipsoid[i].m_cellId), "Ellipsoid in halo-cell!");
    }
  }
}

/// Reduce the number of active particles by using some artificial cell-based criterias
/// \tparam nDim
template <MInt nDim>
void LPT<nDim>::reduceParticles() {
  countParticlesInCells();
  perCellStats();

  vector<MInt> cellIds;

  // determine critical cells by checking conditions
  static constexpr MInt cellPartLimit = 150;
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    if(a_noParticlesInCell(cellId) < cellPartLimit) {
      continue;
    }
    cellIds.emplace_back(cellId);
  }

  // TODO labels:LPT make settable
  MInt totalMerged = 0;
  // 1. remove all particles with low relative velocity
  const MFloat dvLimit = 0.5;
  // 2. remove all particles with a small diameter
  const MFloat diamLimit = 1E-5;
  // 3. combine the rest by merging with the next valid parcel
  for(auto const& cellId : cellIds) {
    // iterator over all particles in this cell
    auto particlesInCell = m_cellToPartMap.equal_range(cellId);
    LPTSpherical<nDim>* validPartToMerge = nullptr;

    MInt filterCount = 0;
    MInt mergedParticles = 0;
    for(auto i = particlesInCell.first; i != particlesInCell.second; i++) {
      auto& particle = i->second;
      const MFloat magDv =
          sqrt(POW2(particle->m_accel.at(0)) + POW2(particle->m_accel.at(1)) + POW2(particle->m_accel.at(2)))
          * m_timeStep;
      const MFloat partD = particle->m_diameter;
      if(magDv < dvLimit && partD < diamLimit) {
        if(validPartToMerge == nullptr) {
          validPartToMerge = particle;
        } else {
          mergeParticles(validPartToMerge, particle);

          mergedParticles++;
          validPartToMerge = nullptr;
        }
        filterCount++;
      }
    }


    if(filterCount > 0) {
      totalMerged += mergedParticles;
    }
  }

  // sort by diameter smallest last (otherwise deletion takes forever)
  own_sort(m_partList, sortDesc_particleAfterDiameter<LPTSpherical<nDim>>::compare);

  // erase particles with 0 mass
  MInt noDeleted = 0;
  for(MInt i = (MInt)m_partList.size() - 1; i >= 0; i--) {
    if(smallParticle(m_partList[i])) {
#ifdef LPT_DEBUG
      cerr << "deleted small particle " << i << " at timestep " << m_timeStep << endl;
#endif
      ++noDeleted;
      m_partList.erase(m_partList.begin() + i);
    }
    //    else {
    //      break;
    //    }
  }

  // recover original sort
  own_sort(m_partList, sort_particleAfterPartIds<LPTSpherical<nDim>>::compare);

  if(noDeleted > 0) {
    cerr << domainId() << ": total particles merged " << totalMerged << " deleted " << noDeleted << endl;
  }
}

/// Merge both particles to A
/// \tparam nDim
/// \param partA
/// \param partB
template <MInt nDim>
void LPT<nDim>::mergeParticles(LPTSpherical<nDim>* partA, LPTSpherical<nDim>* partB) {
  if(smallParticle(*partA) || smallParticle(*partB)) {
    TERM(-1);
  }

  // mass average value
  auto massAvg = [](const MFloat massA, MFloat* const valueA, const MFloat massB, const MFloat* const valueB) {
    for(MInt dim = 0; dim < nDim; dim++) {
      valueA[dim] = (massA * valueA[dim] + massB * valueB[dim]) / (massA + massB);
    }
  };
  // calculate new direction of parcel while keeping kinetic energy constant
  auto massAvg2 = [](const MFloat massA, MFloat* const valueA, const MFloat massB, const MFloat* const valueB) {
    // calculate unit direction vector
    array<MFloat, nDim> dir{};
    MFloat resultingV =
        pow((massA * POW2(math::norm(valueA, nDim)) + massB * POW2(math::norm(valueB, nDim))) / (massA + massB), 0.5);
    for(MInt dim = 0; dim < nDim; ++dim) {
      dir.at(dim) = massA * valueA[dim] + massB * valueB[dim];
    }
    math::normalize(dir);

    for(MInt dim = 0; dim < nDim; dim++) {
      valueA[dim] = dir.at(dim) * resultingV;
    }
  };

  const MInt noA = partA->m_noParticles;
  const MInt noB = partB->m_noParticles;

  const MFloat massA = partA->sphericalMass() * noA;
  const MFloat massB = partB->sphericalMass() * noB;

  massAvg(massA, &partA->m_oldPos[0], massB, &partB->m_oldPos[0]);
  massAvg(massA, &partA->m_position[0], massB, &partB->m_position[0]);
  massAvg(massA, &partA->m_accel[0], massB, &partB->m_accel[0]);
  massAvg(massA, &partA->m_oldAccel[0], massB, &partB->m_oldAccel[0]);
  massAvg2(massA, &partA->m_oldVel[0], massB, &partB->m_oldVel[0]);

  //  cerr << "EKIN before " << massA * POW2(math::norm(&partA->m_velocity[0], 3)) + massB * POW2(math::norm
  //  (&partB->m_velocity[0], 3)) << endl;

  //  auto temp = partA->m_velocity;
  //  cerr << "asddas " << temp[0] << endl;
  //  massAvg(massA, &temp[0], massB, &partB->m_velocity[0]);
  //  cerr << "before " << temp[0] << endl;
  massAvg2(massA, &partA->m_velocity[0], massB, &partB->m_velocity[0]);
  //  cerr << "after " << partA->m_velocity[0] << endl;
  //  cerr << "EKIN after " << (massA+massB) * POW2(math::norm(&partA->m_velocity[0], 3)) << endl;
  //  cerr << "EKIN after2 " << (massA+massB) * POW2(math::norm(&temp[0], 3)) << endl;
  //  TERM(-1);
  partA->m_creationTime = (massA * partA->m_creationTime + massB * partB->m_creationTime) / (massA + massB);
  partA->m_temperature = (massA * partA->m_temperature + massB * partB->m_temperature) / (massA + massB);
  partA->m_breakUpTime = (massA * partA->m_breakUpTime + massB * partB->m_breakUpTime) / (massA + massB);
  partA->m_shedDiam = (massA * partA->m_shedDiam + massB * partB->m_shedDiam) / (massA + massB);
  partA->m_dM = (massA * partA->m_dM + massB * partB->m_dM) / (massA + massB);

  // calculate new diameter based on mass -> which does not keep SMD constant....
  //  partA->m_diameter = pow((partA->m_noParticles * POW3(partA->m_diameter) + partB->m_noParticles * POW3
  //      (partB->m_diameter))/(partA->m_noParticles + partB->m_noParticles), 1.0/3.0);

  // to keep SMD constant merge to SMD value diameter
  partA->m_diameter = (noA * POW3(partA->m_diameter) + noB * POW3(partB->m_diameter))
                      / (noA * POW2(partA->m_diameter) + noB * POW2(partB->m_diameter));
  partB->m_diameter = 0; // mark for deletion

  const MFloat newMassPerDrop = partA->sphericalMass();
  partA->m_noParticles = (massA + massB) / newMassPerDrop;
}


template <MInt nDim>
void LPT<nDim>::motionEquation(const MInt offset) {
  RECORD_TIMER_START(m_timers[Timers::Motion]);
  ASSERT(grid().checkNghbrIds(), "");

#ifdef _OPENMP
#pragma omp parallel default(none) shared(offset)
#endif
  {
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(MInt i = offset; i < a_noParticles(); i++) {
      if(!m_partList[i].isInvalid()) {
        m_partList[i].motionEquation();
      }
    }
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(MInt i = offset; i < a_noEllipsoidalParticles(); i++) {
      if(!m_partListEllipsoid[i].isInvalid()) {
        m_partListEllipsoid[i].motionEquation();
      }
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::Motion]);
}

template <MInt nDim>
void LPT<nDim>::evaporation(const MInt offset) {
  if(!m_heatCoupling && !m_evaporation) return;

  RECORD_TIMER_START(m_timers[Timers::Energy]);
  ASSERT(grid().checkNghbrIds(), "");

#ifdef _OPENMP
#pragma omp parallel default(none) shared(offset)
#endif
  {
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(MInt i = offset; i < a_noParticles(); i++) {
      if(!m_partList[i].isInvalid()) {
        m_partList[i].energyEquation();
      }
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::Energy]);
}

template <MInt nDim>
void LPT<nDim>::coupling(MInt offset) {
  MBool afterBalance = false;
  if(offset > 0) {
    RECORD_TIMER_START(m_timers[Timers::SourceTerms]);
  } else {
    offset = 0;
    afterBalance = true;
  }

  ASSERT(grid().checkNghbrIds(), "");

  for(MInt i = offset; i < a_noParticles(); i++) {
    if(!m_partList[i].isInvalid() || m_partList[i].fullyEvaporated()) {
      m_partList[i].coupling();
    }
  }
  for(MInt i = offset; i < a_noEllipsoidalParticles(); i++) {
    if(!m_partListEllipsoid[i].isInvalid()) {
      m_partListEllipsoid[i].coupling();
    }
  }
  if(!afterBalance) {
    RECORD_TIMER_STOP(m_timers[Timers::SourceTerms]);
  }
}

template <MInt nDim>
void LPT<nDim>::advanceParticles(const MInt offset) {
  ASSERT(grid().checkNghbrIds(), "");
#ifdef _OPENMP
#pragma omp parallel default(none) shared(offset)
#endif
  {
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(MInt i = offset; i < a_noParticles(); i++) {
      if(!m_partList[i].isInvalid()) {
        m_partList[i].advanceParticle();
      }
    }
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for(MInt i = offset; i < a_noEllipsoidalParticles(); i++) {
      if(!m_partListEllipsoid[i].isInvalid()) {
        m_partListEllipsoid[i].advanceParticle();
      }
    }
  }
}

/** \brief collision step with LPT wall
 *  \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::wallCollision() {
  if(!m_wallCollisions) return;

  RECORD_TIMER_START(m_timers[Timers::Wall]);

  // loop over all particles and call particleCollision step if necessary
  for(MInt i = 0; i < a_noParticles(); i++) {
    const MInt cellId = m_partList[i].m_cellId;
    const MInt oldCellId = m_partList[i].m_oldCellId;

    if(cellId < 0) continue;

    // if is or was in a boundary cell
    // if the cell leaves the cartesin grid the current cellId may be -1,
    // but in this case the oldCellId should be a bndryCell!
    if(a_isBndryCell(cellId) || a_isBndryCell(oldCellId)) {
      m_partList[i].particleWallCollision();
      m_partList[i].wallParticleCollision();
    }
  }

  // for ellipsoids
  for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
    const MInt cellId = m_partListEllipsoid[i].m_cellId;
    const MInt oldCellId = m_partListEllipsoid[i].m_oldCellId;
    if(cellId < 0) continue;
    if(a_isBndryCell(cellId) || a_isBndryCell(oldCellId)) {
      m_partListEllipsoid[i].particleWallCollision();
      m_partListEllipsoid[i].wallParticleCollision();
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Wall]);
}

/** void LPT::particleRespawn()
 *
 *  \brief respawn of Particles
 *
 *  version: December-09, split into seperate function from timestep() Dez-2013
 *  @author Rudie Kunnen, Christoph Siewert
 *
 */
template <MInt nDim>
void LPT<nDim>::particleRespawn() {
  partType thisPart{};
  list<partType> partRespawnList; // queue<partType> partQueue;

  {
    auto deleteOffset = partition(m_partList.begin(), m_partList.end(), activeParticle<nDim>);
    auto offset = deleteOffset;
    while(offset != m_partList.end()) {
      if(offset->toBeRespawn()) {
        thisPart.partId = offset->m_partId;
        thisPart.cellId = offset->m_cellId;
        thisPart.diam = offset->m_diameter;
        thisPart.densRatio = offset->m_densityRatio;
        partRespawnList.push_back(thisPart); // partQueue.push(thisPart);
      }
      ++offset;
    }
    m_partList.erase(deleteOffset, m_partList.end());
  }

  if(noDomains() > 1) {
    MInt const respawnSize = 3;
    MInt const respawnSendSize = respawnSize + 3;
    if(domainId() == m_respawnDomain) {
      MIntScratchSpace arrayNoOfParts(noDomains(), AT_, "arrayNoOfParts");
      arrayNoOfParts.p[m_respawnDomain] = (MInt)partRespawnList.size();
      MPI_Gather(MPI_IN_PLACE, 1, MPI_INT, &arrayNoOfParts.p[0], 1, MPI_INT, m_respawnDomain, mpiComm(), AT_, "INPLACE",
                 "arrayNoOfParts");

      partRespawnList.sort(sort_respawnParticleAfterCellIds<partType>::compare);

      MInt globalNoRespawnPart = 0;
      for(MInt i = 0; i < noDomains(); ++i) {
        globalNoRespawnPart += arrayNoOfParts.p[i];
        arrayNoOfParts.p[i] *= respawnSize;
      }
      MIntScratchSpace displs((noDomains() + 1), AT_, "displs");
      displs.p[0] = 0;
      for(MInt i = 1; i <= noDomains(); ++i) {
        displs.p[i] = displs.p[i - 1] + arrayNoOfParts.p[i - 1];
      }

      MFloatScratchSpace recvbuf(max(displs.p[noDomains()], 1), AT_, "recvbuf");
      MInt counter = displs.p[m_respawnDomain];
      while(!partRespawnList.empty()) {
        thisPart = partRespawnList.front();
        recvbuf.p[counter++] = MFloat(thisPart.partId);
        recvbuf.p[counter++] = thisPart.diam;
        recvbuf.p[counter++] = thisPart.densRatio;
        partRespawnList.pop_front(); // partQueue.pop();
      }
      MPI_Gatherv(MPI_IN_PLACE, arrayNoOfParts.p[m_respawnDomain], MPI_DOUBLE, &recvbuf.p[0], &arrayNoOfParts.p[0],
                  &displs.p[0], MPI_DOUBLE, m_respawnDomain, mpiComm(), AT_, "INPLACE", "m_respawnDomain");

      // generate random things on top of the id, dia, densRation. 1 which cell, 2 random values for
      // coordinates
      MIntScratchSpace bufferSizes(m_noRespawnDomains, AT_, "bufferSizes");
      for(MInt domain = 0; domain < m_noRespawnDomains; ++domain) {
        bufferSizes.p[domain] = 0;
      }
      MIntScratchSpace cellIds(globalNoRespawnPart, AT_, "cellIds");
      MFloatScratchSpace randCoord1(globalNoRespawnPart, AT_, "randCoord1");
      MFloatScratchSpace randCoord2(globalNoRespawnPart, AT_, "randCoord2");
      uniform_int_distribution<MInt> dist_cells(0, m_respawnGlobalDomainOffsets[m_noRespawnDomains] - 1);
      uniform_real_distribution<MFloat> dist_coord(0, 1);
      for(MInt i = 0; i < globalNoRespawnPart; ++i) {
        cellIds.p[i] = dist_cells(randomRespawn());
        randCoord1.p[i] = dist_coord(randomRespawn()) - 0.5;
        randCoord2.p[i] = dist_coord(randomRespawn()) - 0.5;

        for(MInt domain = 1; domain <= m_noRespawnDomains; ++domain) {
          if(cellIds.p[i] < m_respawnGlobalDomainOffsets[domain]) {
            bufferSizes.p[domain - 1] += 1;
            break;
          }
        }
      }
      MInt ownRankInRespawnDomains = -1;
      ScratchSpace<MPI_Request> receiveRequest(m_noRespawnDomains, AT_, "receiveRequest");
      for(MInt i = 0; i < m_noRespawnDomains; ++i) {
        receiveRequest[i] = MPI_REQUEST_NULL;
        if(m_respawnDomainRanks[i] == domainId()) {
          ownRankInRespawnDomains = i;
          continue;
        }

        MPI_Isend(&bufferSizes.p[i], 1, MPI_INT, m_respawnDomainRanks[i], LPT_MPI_PARTICLE_RESPAWN_TAG, mpiComm(),
                  &receiveRequest[i], AT_, "bufferSizes");
      }
      MIntScratchSpace countsPerDomain(m_noRespawnDomains, AT_, "countsPerDomain");
      countsPerDomain.p[0] = 0;
      for(MInt domain = 1; domain < m_noRespawnDomains; ++domain) {
        countsPerDomain.p[domain] = countsPerDomain.p[domain - 1];
        if((domain - 1) != ownRankInRespawnDomains) {
          countsPerDomain.p[domain] += bufferSizes.p[domain - 1] * respawnSendSize;
        }
      }
      MFloatScratchSpace sendbuf(respawnSendSize * (globalNoRespawnPart - bufferSizes.p[ownRankInRespawnDomains]), AT_,
                                 "sendbuf");
      for(MInt i = 0; i < globalNoRespawnPart; ++i) {
        MInt thisDomain = -1;
        for(MInt domain = 1; domain <= m_noRespawnDomains; ++domain) {
          if(cellIds.p[i] < m_respawnGlobalDomainOffsets[domain]) {
            thisDomain = domain - 1;
            break;
          }
        }
        if(thisDomain == ownRankInRespawnDomains) {
          continue;
        }
        MInt address = countsPerDomain.p[thisDomain];
        for(MInt inner = 0; inner < respawnSize; ++inner) {
          sendbuf.p[address + inner] = recvbuf.p[i * respawnSize + inner];
        }
        sendbuf.p[address + respawnSize + 0] = MFloat(cellIds.p[i] - m_respawnGlobalDomainOffsets[thisDomain]);
        sendbuf.p[address + respawnSize + 1] = randCoord1.p[i];
        sendbuf.p[address + respawnSize + 2] = randCoord2.p[i];

        countsPerDomain.p[thisDomain] += respawnSendSize;
      }
      // reset countsPerDomain
      countsPerDomain.p[0] = 0;
      for(MInt domain = 1; domain < m_noRespawnDomains; ++domain) {
        countsPerDomain.p[domain] = countsPerDomain.p[domain - 1];
        if((domain - 1) != ownRankInRespawnDomains) {
          countsPerDomain.p[domain] += bufferSizes.p[domain - 1] * respawnSendSize;
        }
      }
      MPI_Waitall(m_noRespawnDomains, &receiveRequest[0], MPI_STATUSES_IGNORE, AT_);
      for(MInt domain = 0; domain < m_noRespawnDomains; ++domain) {
        if((domain == ownRankInRespawnDomains) || (bufferSizes.p[domain] == 0)) {
          receiveRequest[domain] = MPI_REQUEST_NULL;
        } else {
          MInt address = countsPerDomain.p[domain];
          MPI_Isend(&sendbuf.p[address], bufferSizes.p[domain] * respawnSendSize, MPI_DOUBLE,
                    m_respawnDomainRanks[domain], LPT_MPI_PARTICLE_RESPAWN_TAG, mpiComm(), &receiveRequest[domain], AT_,
                    "sendbuf");
        }
      }

      // auspacken
      if(bufferSizes.p[ownRankInRespawnDomains] != 0) {
        for(MInt i = 0; i < globalNoRespawnPart; ++i) {
          if((cellIds.p[i] < m_respawnGlobalDomainOffsets[ownRankInRespawnDomains + 1])
             && (cellIds.p[i] >= m_respawnGlobalDomainOffsets[ownRankInRespawnDomains])) {
            LPTSpherical<nDim> thisParticle;
            MFloat x[3];
            MFloat v[3];

            counter = i * respawnSize;
            MInt id = cellIds.p[i] - m_respawnGlobalDomainOffsets[ownRankInRespawnDomains];
            // obtain random positions within this cell
            MFloat cellLength = c_cellLengthAtCell(m_respawnCells.at((MUlong)id));

            x[m_respawn - 1] = m_respawnPlane;

            x[m_respawn % 3] =
                cellLength * randCoord1.p[i] + c_coordinate(m_respawnCells.at((MUlong)id), m_respawn % 3);
            x[(m_respawn + 1) % 3] =
                cellLength * randCoord2.p[i] + c_coordinate(m_respawnCells.at((MUlong)id), (m_respawn + 1) % 3);

            thisParticle.m_cellId = m_respawnCells.at((MUlong)id);
            thisParticle.m_partId = MInt(recvbuf.p[counter++]);
            thisParticle.m_diameter = recvbuf.p[counter++];
            thisParticle.m_densityRatio = recvbuf.p[counter];
            thisParticle.updateProperties();
            thisParticle.firstStep() = true;
            interpolateVariablesLS<0, nDim>(m_respawnCells[id], x, v);
            for(MInt j = 0; j < nDim; j++) {
              thisParticle.m_position[j] = x[j];
              thisParticle.m_oldPos[j] = x[j];
              thisParticle.m_velocity[j] =
                  v[j] + m_terminalVelocity[thisParticle.m_diameter] * m_partList[0].s_Frm[j] / m_FrMag;
              thisParticle.m_oldVel[j] =
                  v[j] + m_terminalVelocity[thisParticle.m_diameter] * m_partList[0].s_Frm[j] / m_FrMag;
              thisParticle.m_accel[j] = F0;
              thisParticle.m_oldFluidVel[j] = a_fluidVelocity(thisParticle.m_cellId, j);
            }

            thisParticle.m_oldFluidDensity = a_fluidDensity(thisParticle.m_cellId);

            m_partList.push_back(thisParticle);
          }
        }
      }

      MPI_Waitall(m_noRespawnDomains, &receiveRequest[0], MPI_STATUSES_IGNORE, AT_);
    } else {
      MInt noPartsToRespawn = 0;
      MPI_Request receiveRequest = MPI_REQUEST_NULL;
      if(!m_respawnCells.empty()) {
        MPI_Irecv(&noPartsToRespawn, 1, MPI_INT, m_respawnDomain, LPT_MPI_PARTICLE_RESPAWN_TAG, mpiComm(),
                  &receiveRequest, AT_, "noPartsToRespawn");
      }

      auto localNoOfParts = (MInt)partRespawnList.size();
      MPI_Gather(&localNoOfParts, 1, MPI_INT, nullptr, 1, MPI_INT, m_respawnDomain, mpiComm(), AT_, "localNoOfParts",
                 "nullptr");
      MFloatScratchSpace sendbuf(max(localNoOfParts, 1) * respawnSize, AT_, "sendbuf");
      MInt counter = 0;
      while(!partRespawnList.empty()) {
        thisPart = partRespawnList.front();
        sendbuf.p[counter++] = MFloat(thisPart.partId);
        sendbuf.p[counter++] = thisPart.diam;
        sendbuf.p[counter++] = thisPart.densRatio;
        partRespawnList.pop_front(); // partQueue.pop();
      }

      MPI_Gatherv(&sendbuf.p[0], localNoOfParts * respawnSize, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE,
                  m_respawnDomain, mpiComm(), AT_, "sendbuf", "nullptr");

      if(!m_respawnCells.empty()) {
        MPI_Wait(&receiveRequest, MPI_STATUS_IGNORE, AT_);

        if(noPartsToRespawn != 0) {
          MFloatScratchSpace recvbuf(noPartsToRespawn * respawnSendSize, AT_, "recvbuf");
          MPI_Recv(&recvbuf.p[0], noPartsToRespawn * respawnSendSize, MPI_DOUBLE, m_respawnDomain,
                   LPT_MPI_PARTICLE_RESPAWN_TAG, mpiComm(), MPI_STATUS_IGNORE, AT_, "recvbuf");

          for(MInt part = 0; part < noPartsToRespawn; ++part) {
            LPTSpherical<nDim> thisParticle;
            MInt id;
            MInt address;
            MFloat x[3];
            MFloat v[3];

            address = part * respawnSendSize;
            id = MInt(recvbuf.p[address + respawnSize + 0]);
            MFloat cellLength = c_cellLengthAtCell(m_respawnCells.at((MUlong)id));

            x[m_respawn - 1] = m_respawnPlane;

            x[m_respawn % 3] = cellLength * recvbuf.p[address + respawnSize + 1]
                               + c_coordinate(m_respawnCells.at((MUlong)id), m_respawn % 3);
            x[(m_respawn + 1) % 3] = cellLength * recvbuf.p[address + respawnSize + 2]
                                     + c_coordinate(m_respawnCells.at((MUlong)id), (m_respawn + 1) % 3);
            thisParticle.m_cellId = m_respawnCells.at((MUlong)id);
            thisParticle.m_partId = MInt(recvbuf.p[address + 0]);
            thisParticle.m_diameter = recvbuf.p[address + 1];
            thisParticle.m_densityRatio = recvbuf.p[address + 2];
            thisParticle.updateProperties();
            thisParticle.firstStep() = true;
            interpolateVariablesLS<0, nDim>(m_respawnCells[id], x, v);
            for(MInt j = 0; j < nDim; j++) {
              thisParticle.m_position[j] = x[j];
              thisParticle.m_oldPos[j] = x[j];
              thisParticle.m_velocity[j] =
                  v[j] + m_terminalVelocity[thisParticle.m_diameter] * m_partList[0].s_Frm[j] / m_FrMag;
              thisParticle.m_oldVel[j] =
                  v[j] + m_terminalVelocity[thisParticle.m_diameter] * m_partList[0].s_Frm[j] / m_FrMag;
              thisParticle.m_accel[j] = F0;
              thisParticle.m_oldFluidVel[j] = a_fluidVelocity(thisParticle.m_cellId, j);
            }
            thisParticle.m_oldFluidDensity = a_fluidDensity(thisParticle.m_cellId);
            m_partList.push_back(thisParticle);
          }
        }
      }
    }
  }
  if(m_ellipsoids) {
    mTerm(-1, AT_, "particleRespawn not implemented for ellipsoidal particle");
  }
}

/**
 * \fn void LPT::calculateTerminalVelocities()
 *
 *  \brief calculates the terminal velocity of particles used for initialization
 *
 *  version: Dez-2013
 *  @author Christoph Siewert
 *
 */
template <MInt nDim>
void LPT<nDim>::calculateTerminalVelocities() {
  /*! \page propertyPageLPT
      \section particleDiameters
      <code>MFloat particleDiameter</code>\n
      default = <code>(0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.00006, 0.00007, 0.00008,
     0.00009, 0.00010)</code> \n \n diameters of particles. \n Keywords: <i>PARTICLE</i>
   */
  MFloat particleDiametersDefaults[10] = {0.00001, 0.00002, 0.00003, 0.00004, 0.00005,
                                          0.00006, 0.00007, 0.00008, 0.00009, 0.00010};
  MInt noDiffParticles = 10;
  //  MFloat* helpPointer = nullptr;
  if(Context::propertyExists("particleDiameters", solverId())) {
    noDiffParticles = Context::propertyLength("particleDiameters", solverId());
  } // else {
    //  helpPointer = &particleDiametersDefaults[0];
    //}
  MFloatScratchSpace particleDiameters(noDiffParticles, AT_, "particleDiameters");
  for(MInt i = 0; i < noDiffParticles; i++) {
    particleDiameters[i] =
        Context::getSolverProperty<MFloat>("particleDiameters", solverId(), AT_, &particleDiametersDefaults[i], i);
  }

  m_terminalVelocity.clear();

  if((m_FrMag > 1.0e-12) || (m_dragModelType == 0)) {
    for(MInt i = 0; i < noDiffParticles; ++i) {
      const MFloat tmp = particleDiameters[i];
      m_terminalVelocity.insert(make_pair(tmp, 0.0));
    }
  } else {
    MInt const maxIter = 1000;
    MFloat epsilon = 1.0e-8, min = F0, max = F1;
    MFloat value = F0;
    for(MInt i = 0; i < noDiffParticles; ++i) {
      MInt side = 0;
      MFloat fr = NAN;
      MFloat fs = 0;
      MFloat ft = 0;

      MInt n = 0;
      for(n = 1; n <= maxIter; n++) {
        value = (fs * max - ft * min) / (fs - ft);
        fr = 0;
        if(fabs(fr) < epsilon * value) {
          break;
        }


        if(fr * ft > 0) {
          max = value;
          ft = fr;
          if(side == -1) {
            fs /= 2;
          }
          side = -1;
        } else {
          min = value;
          fs = fr;
          if(side == +1) {
            ft /= 2;
          }
          side = +1;
        }
      }
      m_terminalVelocity.insert(pair<MFloat, MFloat>(particleDiameters[i], value));
    }
  }
}

/**
 *  \brief Calls exchange for the different particle types
 *  \author  Laurent Andre
 */
template <MInt nDim>
void LPT<nDim>::exchangeParticles(const MBool collOnly, const MInt offset, const MBool allowNonLeaf) {
  if(m_ellipsoids) {
    exchangeParticles_<LPTEllipsoidal<nDim>>(collOnly, offset, allowNonLeaf);
  } else {
    exchangeParticles_<LPTSpherical<nDim>>(collOnly, offset, allowNonLeaf);
  }
}

/**
 *  \brief exchange particles to neighboring domains,
 *         can be either non-blocking or blocking communication!
 *  @author Rudie Kunnen, Sven Berger, Tim Wegmann, Johannes Grafen
 */
template <MInt nDim>
template <class LPTParticle>
void LPT<nDim>::exchangeParticles_(const MBool collOnly, const MInt offset, const MBool allowNonLeaf) {
  TRACE();

  if(noDomains() == 1) return;

  MBool foundCell = false;
  RECORD_TIMER_START(m_timers[Timers::Exchange]);

  vector<LPTParticle>& partList = a_particleList<LPTParticle>();
  MInt noParticles = a_noParticles<LPTParticle>();

  // check for transfer
  for(MInt id = offset; id < noParticles; id++) {
    auto& part = partList[id];
    if((collOnly && !part.hasCollided()) || // if coll only skip all non collided particles
       (part.wasSend())) {                  // exclude ghost particles
      continue;
    }

    foundCell = true;
    if(!collOnly && m_collisions > 0 && part.isWindow()) {
      // particle must be copied to the halo cell(s) of neighboring solver(s)
      // this is only necessary when collisions are activated otherwise this just leads to
      // overhead...
      foundCell = pushToQueue<LPTParticle>(m_pointsToHalo, id);
    } else if(!m_periodicBC && part.reqSend()) {
      // particle must be copied to a window cell of a neighboring solver
      // reduce the search to only halo-cells!

      if(a_isHalo(part.m_cellId)) {
        foundCell = pushToQueue<LPTParticle>(m_pointsToWindow, id);
        if(allowNonLeaf && !foundCell) {
          // if allowing for non-leaf cell transfer, check for all halo cells
          for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
            for(MInt i = 0; i < noHaloCells(d); i++) {
              if(part.m_cellId == grid().haloCell(d, i)) {
                sendQueueType<nDim> thisSend{};
                thisSend.partPos = id;
                thisSend.toCellId = i;
                m_queueToSend[d].push(thisSend);
                foundCell = true;
              }
            }
          }
        }
      }
    } else if(m_periodicBC && part.reqSend()) {
      std::array<MFloat, nDim> tmpShift{F0};
      if(grid().isPeriodic(part.m_cellId)) {
        for(MInt n = 0; n < nDim; n++) {
          if(grid().raw().periodicCartesianDir(n)) {
            MFloat cellCoordinate = c_coordinate(part.m_cellId, n);
            MInt sign1 = maia::math::sgn(m_globBbox[n] - cellCoordinate);
            MInt sign2 = maia::math::sgn(m_globBbox[n + nDim] - cellCoordinate);
            MFloat sign = F0;
            if(sign1 == 1 && sign2 == 1) {
              sign = F1;
            } else if(sign1 == -1 && sign2 == -1) {
              sign = -F1;
            } else {
              sign = F0;
            }
            tmpShift[n] = m_globDomainLength[n] * sign;
          }
        }
        // new position of particle
        for(MInt n = 0; n < nDim; n++) {
          partList[id].m_position[n] += tmpShift[n];
        }
      }
      // find cell where particle is supposed to be send to
      if(a_isHalo(part.m_cellId)) {
        foundCell = pushToQueue<LPTParticle>(m_pointsToWindow, id);
        if(allowNonLeaf && !foundCell) {
          // if allowing for non-leaf cell transfer, check for all halo cells
          for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
            for(MInt i = 0; i < noHaloCells(d); i++) {
              if(part.m_cellId == grid().haloCell(d, i)) {
                sendQueueType<nDim> thisSend{};
                thisSend.partPos = id;
                thisSend.toCellId = i;
                m_queueToSend[d].push(thisSend);
                foundCell = true;
              }
            }
          }
        }
      }
    }

    if(!foundCell) {
      cerr << " globalTimeStep " << globalTimeStep << endl;
      cerr << part.m_cellId << endl;
      cerr << a_isHalo(part.m_cellId) << endl;
      m_log << "TRANSFER status &1 or &4: no corresponding cell found! Original cell " << part.m_cellId
            << " properties " << endl;
      TERMM(-1, "????");
    }
  }

  // copy particles to the neighboring domains
  // blocking comminucation version
  if(!m_nonBlockingComm || collOnly) {
    sendAndReceiveParticles<LPTParticle>(allowNonLeaf);
  } else {
    if(allowNonLeaf) {
      sendParticles<true, LPTParticle>();
      m_nonBlockingStage = 2;
    } else {
      sendParticles<false, LPTParticle>();
      m_nonBlockingStage = 1;
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::Exchange]);
}

/**
 *  \brief Particle transfer between solvers
 *  The particles to be sent are in the queueToSend array
 *  MBool mayExist is true for particle transfer after collisions, since
 *  in that case the particle may already be found on this rank
 *  \author Rudie Kunnen, Sven Berger, Tim Wegmann
 */
template <MInt nDim>
template <class LPTParticle>
void LPT<nDim>::sendAndReceiveParticles(const MBool allowNonLeaf) {
  TRACE();

  MPI_Status status{};

  // 1) exchange the number of particles to be sent:
  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    m_sendSize[d] = m_queueToSend[d].size();
    if(m_sendSize[d] > m_exchangeBufferSize) {
      mTerm(1, AT_, "Increase particle exchange buffer size!");
    }
    MPI_Isend(&m_sendSize[d], 1, MPI_INT, grid().neighborDomain(d), mpiTag("PARTICLE_COUNT"), mpiComm(),
              &m_mpi_reqSendSize[d], AT_, "noToSend");
  }

  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    MPI_Recv(&m_recvSize[d], 1, MPI_INT, grid().neighborDomain(d), mpiTag("PARTICLE_COUNT"), mpiComm(), &status, AT_,
             "m_recvSize");
  }

  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    MPI_Wait(&m_mpi_reqSendSize[d], &status, AT_);
  }

  // 2) prepare/write the two send buffers (MInt, MFloat)
  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    if(m_queueToSend[d].empty()) {
      // skip nothing to send
      m_mpi_reqSendFloat[d] = MPI_REQUEST_NULL;
      m_mpi_reqSendInt[d] = MPI_REQUEST_NULL;
      continue;
    }

    vector<LPTParticle*> particlesToSend;
    vector<MInt> cellIds;
    vector<LPTParticle>& partList = a_particleList<LPTParticle>();

    while(!m_queueToSend[d].empty()) {
      sendQueueType<nDim> thisSend = m_queueToSend[d].front();
      m_queueToSend[d].pop();
      MInt partPos = thisSend.partPos;

      particlesToSend.push_back(&partList[partPos]);
      cellIds.push_back(thisSend.toCellId);
    }

    MInt noParticles = particlesToSend.size();

    packParticles(particlesToSend, m_intSendBuffer[d].get(), m_sendBuffer[d].get(), cellIds);

    for(auto* part : particlesToSend) {
      if(part->reqSend()) {
        part->wasSend() = true;
        part->reqSend() = false;
        part->isInvalid() = true;
      }
    }

    // 3) send the buffers
    MPI_Isend(&m_sendBuffer[d][0], noParticles * elemPerP<LPTParticle>(), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("PARTICLE_FLOAT"), mpiComm(), &m_mpi_reqSendFloat[d], AT_, "m_sendBuffer");

    MPI_Isend(&m_intSendBuffer[d][0], noParticles * intElemPerP<LPTParticle>(), MPI_INT, grid().neighborDomain(d),
              mpiTag("PARTICLE_INT"), mpiComm(), &m_mpi_reqSendInt[d], AT_, "m_intSendBuffer");
  }

  // 4) receive the buffers
  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    if(m_recvSize[d] > 0) {
      MPI_Recv(&m_recvBuffer[d][0], mMin(bufSize<LPTParticle>(), m_recvSize[d] * elemPerP<LPTParticle>()), MPI_DOUBLE,
               grid().neighborDomain(d), mpiTag("PARTICLE_FLOAT"), mpiComm(), &status, AT_, "m_recvBuffer");
      MPI_Recv(&m_intRecvBuffer[d][0], mMin(intBufSize<LPTParticle>(), m_recvSize[d] * intElemPerP<LPTParticle>()),
               MPI_INT, grid().neighborDomain(d), mpiTag("PARTICLE_INT"), mpiComm(), &status, AT_, "m_intRecvBuffer");
    }
  }

  // 5) wait for completion of communication
  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    MPI_Wait(&m_mpi_reqSendFloat[d], &status, AT_);
    MPI_Wait(&m_mpi_reqSendInt[d], &status, AT_);
  }

  // 6) interpret the receive buffers and add incomming particles
  for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
    // nothing to receive skip
    if(m_recvSize[d] == 0) continue;

    if(!m_periodicBC) {
      unpackParticles<LPTParticle, false>(m_recvSize[d], &m_intRecvBuffer[d][0], &m_recvBuffer[d][0], d, allowNonLeaf);
    } else {
      unpackParticles<LPTParticle, true>(m_recvSize[d], &m_intRecvBuffer[d][0], &m_recvBuffer[d][0], d, allowNonLeaf);
    }
  }
}


/**
 *  \brief send particles non-blocking
 *  \author Tim Wegmann
 */
template <MInt nDim>
template <MBool allNeighbors, class LPTParticle>
void LPT<nDim>::sendParticles() {
  TRACE();

  if(m_nonBlockingStage != 0) {
    mTerm(1, AT_, "Incorrect non-blocking stage!");
  }

  RECORD_TIMER_START(m_timers[Timers::Exchange1]);

  // wait on all previously send before send new
  if(m_openParticleInjSend) {
    MPI_Waitall(noNeighborDomains(), &m_mpi_reqSendFloat[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(noNeighborDomains(), &m_mpi_reqSendInt[0], MPI_STATUSES_IGNORE, AT_);
    m_openParticleInjSend = false;
  }

  if(m_openParticleSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_mpi_reqSendFloat[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_mpi_reqSendInt[0], MPI_STATUSES_IGNORE, AT_);
    m_openParticleSend = false;
  }

  const MInt noNeighborsSend = allNeighbors ? grid().noNeighborDomains() : grid().noLeafSendNeighborDomains();
  const MInt noNeighborsRecv = allNeighbors ? grid().noNeighborDomains() : grid().noLeafRecvNeighborDomains();

  // 3) prepare/write the two send buffers (MInt, MFloat)
  for(MInt n = 0; n < noNeighborsSend; n++) {
    MInt d = n;
    if(!allNeighbors) {
      d = grid().leafSendNeighborDomain(n);
    }
    if(m_queueToSend[d].empty()) {
      // skip nothing to send
      m_mpi_reqSendFloat[n] = MPI_REQUEST_NULL;
      m_intSendBuffer[d][0] = -intElemPerP<LPTParticle>();
      MPI_Isend(&m_intSendBuffer[d][0], 0, MPI_INT, grid().neighborDomain(d), mpiTag("PARTICLE_INT"), mpiComm(),
                &m_mpi_reqSendInt[n], AT_, "m_intSendBuffer");
      continue;
    }

    vector<LPTParticle*> particlesToSend;
    RECORD_TIMER_START(m_timers[Timers::Exchange3]);
    vector<MInt> cellIds;
    vector<LPTParticle>& partList = a_particleList<LPTParticle>();

    while(!m_queueToSend[d].empty()) {
      sendQueueType<nDim> thisSend = m_queueToSend[d].front();
      m_queueToSend[d].pop();
      MInt partPos = thisSend.partPos;

      particlesToSend.push_back(&partList[partPos]);
      cellIds.push_back(thisSend.toCellId);
    }

    MInt noParticles = particlesToSend.size();
    packParticles(particlesToSend, &m_intSendBuffer[d][0], &m_sendBuffer[d][0], cellIds);

    for(auto& part : particlesToSend) {
      if(part->reqSend()) {
        part->wasSend() = true;
        part->reqSend() = false;
        part->isInvalid() = true;
      }
    }
    RECORD_TIMER_STOP(m_timers[Timers::Exchange3]);

    // 3) send the buffers
    MPI_Isend(&m_sendBuffer[d][0], noParticles * elemPerP<LPTParticle>(), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("PARTICLE_FLOAT"), mpiComm(), &m_mpi_reqSendFloat[n], AT_, "m_sendBuffer");

    MPI_Isend(&m_intSendBuffer[d][0], noParticles * intElemPerP<LPTParticle>(), MPI_INT, grid().neighborDomain(d),
              mpiTag("PARTICLE_INT"), mpiComm(), &m_mpi_reqSendInt[n], AT_, "m_intSendBuffer");
  }

  if(allNeighbors) {
    m_openParticleInjSend = true;
  } else {
    m_openParticleSend = true;
  }

  // 1) Probe and get count of the receiving int-buffer
  //   to find the number received particles for each rank
  //   from all ranks which have already called the send and receive those buffers
  const MBool loadWasRunning = this->isLoadTimerRunning();
  if(loadWasRunning) {
    this->stopLoadTimer(AT_);
    this->startIdleTimer(AT_);
    this->disableDlbTimers();
  }

  for(MInt n = 0; n < noNeighborsRecv; n++) {
    MInt d = n;
    if(!allNeighbors) {
      d = grid().leafRecvNeighborDomain(n);
    }
    m_recvSize[d] = -99;
    MInt flag{};
    MPI_Iprobe(grid().neighborDomain(d), mpiTag("PARTICLE_INT"), mpiComm(), &flag, &m_mpi_statusProbe[n],
               "statusProbe");
    if(flag) {
      MPI_Get_count(&m_mpi_statusProbe[n], MPI_INT, &m_recvSize[d], "m_intRecvBuffer");

      m_recvSize[d] = m_recvSize[d] / intElemPerP<LPTParticle>();

      MPI_Irecv(&m_intRecvBuffer[d][0], m_recvSize[d] * intElemPerP<LPTParticle>(), MPI_INT, grid().neighborDomain(d),
                mpiTag("PARTICLE_INT"), mpiComm(), &m_mpi_reqRecvInt[n], AT_, "m_intRecvBuffer");

      if(m_recvSize[d] > 0) {
        MPI_Irecv(&m_recvBuffer[d][0], m_recvSize[d] * elemPerP<LPTParticle>(), MPI_DOUBLE, grid().neighborDomain(d),
                  mpiTag("PARTICLE_FLOAT"), mpiComm(), &m_mpi_reqRecvFloat[n], AT_, "m_recvBuffer");
      } else {
        m_mpi_reqRecvFloat[n] = MPI_REQUEST_NULL;
      }
    }
  }

  if(loadWasRunning) {
    this->reEnableDlbTimers();
    this->stopIdleTimer(AT_);
    this->startLoadTimer(AT_);
  }
  RECORD_TIMER_STOP(m_timers[Timers::Exchange1]);
}


/**
 *  \brief  Calls particle receive function for different particle types
 *  \author Laurent Andre
 */
template <MInt nDim>
template <MBool allNeighbors>
void LPT<nDim>::receiveParticles() {
  if(m_ellipsoids) {
    receiveParticles_<allNeighbors, LPTEllipsoidal<nDim>>();
  } else {
    receiveParticles_<allNeighbors, LPTSpherical<nDim>>();
  }
}


/**
 *  \brief receive particles that have been sent non-blocking before
 *  \author Tim Wegmann
 */
template <MInt nDim>
template <MBool allNeighbors, class LPTParticle>
void LPT<nDim>::receiveParticles_() {
  TRACE();

  if(m_nonBlockingStage < 1) {
    mTerm(1, AT_, "Incorrect non-blocking stage!");
  }
  RECORD_TIMER_START(m_timers[Timers::Exchange]);
  RECORD_TIMER_START(m_timers[Timers::Exchange1]);

  RECORD_TIMER_START(m_timers[Timers::Exchange4]);

  // 1) Probe and get count of the receiving int-buffer
  //   to find the number received particles for each rank
  //   from all ranks, which could not be received at the send-stage

  const MBool loadWasRunning = this->isLoadTimerRunning();
  if(loadWasRunning) {
    this->stopLoadTimer(AT_);
    this->startIdleTimer(AT_);
    this->disableDlbTimers();
  }

  const MInt noNeighborsRecv = allNeighbors ? grid().noNeighborDomains() : grid().noLeafRecvNeighborDomains();

  MBool allReceived = false;
  while(!allReceived) {
    // check if all have been received
    allReceived = true;
    for(MInt n = 0; n < noNeighborsRecv; n++) {
      MInt d = n;
      if(!allNeighbors) {
        d = grid().leafRecvNeighborDomain(n);
      }
      if(m_recvSize[d] == -99) {
        allReceived = false;
        break;
      }
    }

    // loop aver all domains and check for receive again
    for(MInt n = 0; n < noNeighborsRecv; n++) {
      MInt d = n;
      if(!allNeighbors) {
        d = grid().leafRecvNeighborDomain(n);
      }
      MInt flag{};
      MPI_Iprobe(grid().neighborDomain(d), mpiTag("PARTICLE_INT"), mpiComm(), &flag, &m_mpi_statusProbe[n],
                 "statusProbe");
      if(flag) {
        MPI_Get_count(&m_mpi_statusProbe[n], MPI_INT, &m_recvSize[d], "m_intRecvBuffer");
        m_recvSize[d] = m_recvSize[d] / intElemPerP<LPTParticle>();
        MPI_Irecv(&m_intRecvBuffer[d][0], m_recvSize[d] * intElemPerP<LPTParticle>(), MPI_INT, grid().neighborDomain(d),
                  mpiTag("PARTICLE_INT"), mpiComm(), &m_mpi_reqRecvInt[n], AT_, "m_intRecvBuffer");

        if(m_recvSize[d] > 0) {
          MPI_Irecv(&m_recvBuffer[d][0], m_recvSize[d] * elemPerP<LPTParticle>(), MPI_DOUBLE, grid().neighborDomain(d),
                    mpiTag("PARTICLE_FLOAT"), mpiComm(), &m_mpi_reqRecvFloat[n], AT_, "m_recvBuffer");
        } else {
          m_mpi_reqRecvFloat[n] = MPI_REQUEST_NULL;
        }
      }
    }
  }

  if(loadWasRunning) {
    this->reEnableDlbTimers();
    this->stopIdleTimer(AT_);
    this->startLoadTimer(AT_);
  }

  RECORD_TIMER_STOP(m_timers[Timers::Exchange4]);

  // 1) wait until I have received all
  MPI_Waitall(noNeighborsRecv, &m_mpi_reqRecvFloat[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborsRecv, &m_mpi_reqRecvInt[0], MPI_STATUSES_IGNORE, AT_);


  RECORD_TIMER_START(m_timers[Timers::Exchange5]);

  MBool allowNonLeaf = false;
  if(m_nonBlockingStage == 2) allowNonLeaf = true;

  // 2) interpret the receive buffers and add incomming particles
  for(MInt n = 0; n < noNeighborsRecv; n++) {
    MInt d = n;
    if(!allNeighbors) {
      d = grid().leafRecvNeighborDomain(n);
    }
    // nothing to receive skip
    if(m_recvSize[d] == 0) continue;
    if(!m_periodicBC) {
      unpackParticles<LPTParticle, false>(m_recvSize[d], &m_intRecvBuffer[d][0], &m_recvBuffer[d][0], d, allowNonLeaf);
    } else {
      unpackParticles<LPTParticle, true>(m_recvSize[d], &m_intRecvBuffer[d][0], &m_recvBuffer[d][0], d, allowNonLeaf);
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Exchange5]);

  RECORD_TIMER_STOP(m_timers[Timers::Exchange1]);
  m_nonBlockingStage = 0;
  RECORD_TIMER_STOP(m_timers[Timers::Exchange]);
}


/// \brief In case of multisolver, communicate particle id offset to ensure
/// consistent particle ids across the solvers
///
/// \author Sven Berger
/// \date   July 2015
/// \param[in] noParticles Number of particles created in this domain.
template <MInt nDim>
void LPT<nDim>::exchangeOffset(MInt noParticles) {
  if(noDomains() > 1) {
    MIntScratchSpace offsets(noDomains() + 1, AT_, "offsets");
    MPI_Allgather(&noParticles, 1, MPI_INT, &offsets[0], 1, MPI_INT, mpiComm(), AT_, "noParticles", "offsets");

    // determine the particle offset for this solver
    const MInt totalOffset = accumulate(&offsets[0], &offsets[domainId()], 0);

    // determine the maximal partId
    MLong maxPartId = accumulate(&offsets[0], &offsets[noDomains()], 0);

    // output maxPartId
    if(domainId() == 0) {
      m_log << "Domain 0 has " << noParticles << " particles initialized, out of globally " << maxPartId << " particles"
            << endl;
    }

    // now the particles can be consistently numbered throughout the
    // entire domain
    for(auto& part : m_partList) {
      part.m_partId += totalOffset;
    }
    for(auto& ellip : m_partListEllipsoid) {
      ellip.m_partId += totalOffset;
    }
  }
}

/**
 *  \brief write particle snapshot to Netcdf file (using parallel output)
 *
 *  @author Jerry Grimmen, Apr-10
 */
template <MInt nDim>
void LPT<nDim>::writePartData() {
  TRACE();

  // sort particles after paricle-id!
  own_sort(m_partList, sort_particleAfterPartIds<LPTBase<nDim>>::compare);

  queue<partListIteratorConst<nDim>> ncmpiPartQueue;
  MInt ncmpiPartQueueSize = 0;

  auto i1 = m_partList.begin();
  while(i1 != m_partList.end()) {
    if(!(*i1).isInvalid() && ((*i1).m_position[0] > m_xCutOff)) {
      ++ncmpiPartQueueSize;
      ncmpiPartQueue.push(i1);
    }
    i1++;
  }

  // Get Variable partCount[noDomains()]
  MIntScratchSpace ncmpiPartCount(noDomains(), AT_, "ncmpiPartCount");
  MPI_Allgather(&ncmpiPartQueueSize, 1, MPI_INT, &ncmpiPartCount[0], 1, MPI_INT, mpiComm(), AT_, "ncmpiPartQueueSize",
                "ncmpiPartCount");

  // Calculate ncmpiPartCountMax
  ParallelIo::size_type ncmpiPartCountMax = 0;
  for(MInt i = 0; i < noDomains(); ++i) {
    ncmpiPartCountMax += ncmpiPartCount[i];
  }

  // If there are 0 particle in the whole domain, don't write particle data.
  if(ncmpiPartCountMax != 0) {
    // stringstream ncmpistream;
    const MString ncmpiFileName =
        outputDir() + "partData_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

    using namespace maia::parallel_io;
    ParallelIo parallelIo(ncmpiFileName, PIO_REPLACE, mpiComm());

    // Define Attribute timestep (double).
    parallelIo.setAttribute(globalTimeStep, "globalTimestep");
    parallelIo.setAttribute(globalTimeStep, "particleTimestep");
    if(globalTimeStep == 0) {
      parallelIo.setAttribute(0.0, "time");
    } else {
      parallelIo.setAttribute(m_time, "time");
    }
    // Define Dimension partCount [noDomains()] & Define Variable partCount
    // (int) [partCount].
    parallelIo.defineArray(PIO_INT, "partCount", noDomains());

    // Define Dimension partId [ncmpiPartCountMax] & Define Variable partId
    // (int) [partId].
    parallelIo.defineArray(PIO_LONG, "partId", ncmpiPartCountMax);

    // Define Dimension partPos [3 * ncmpiPartCountMax] & Define Variable
    // partPos (double) [partPos].
    parallelIo.defineArray(PIO_FLOAT, "partPos", 3 * ncmpiPartCountMax);

    // Define Dimension partVel [3 * ncmpiPartCountMax] & Define Variable
    // partVel (double) [partVel].
    parallelIo.defineArray(PIO_FLOAT, "partVel", 3 * ncmpiPartCountMax);

    // Define Dimension partDia [ncmpiPartCountMax] & Define Variable partDia
    // (double) [partDia].
    parallelIo.defineArray(PIO_FLOAT, "partDia", ncmpiPartCountMax);

    if(m_activePrimaryBUp || m_activeSecondaryBUp || m_wallCollisions) {
      parallelIo.defineArray(PIO_INT, "partParceledNo", ncmpiPartCountMax);
    }

    if(m_heatCoupling || m_evaporation) {
      parallelIo.defineArray(PIO_FLOAT, "partTemp", ncmpiPartCountMax);
    }

    if(m_domainIdOutput) {
      parallelIo.defineArray(PIO_INT, "domainId", ncmpiPartCountMax);
    }

    parallelIo.setOffset(1, domainId());
    parallelIo.writeArray(&ncmpiPartQueueSize, "partCount");

    // Create arrays to hold the particle datas.
    ParallelIo::size_type ncmpiStart = 0;
    ParallelIo::size_type ncmpiCount = ncmpiPartCount[domainId()];

    if(ncmpiPartCount[domainId()] > 0) {
      for(MInt i = 0; i < domainId(); i++) {
        ncmpiStart += ncmpiPartCount[i];
      }
    }

    ASSERT(ncmpiStart + ncmpiPartCount[domainId()] <= ncmpiPartCountMax,
           "ERROR: Invalid number of particles " + to_string(ncmpiStart + ncmpiPartCount[domainId()]) + " > "
               + to_string(ncmpiPartCountMax));

    ASSERT(ncmpiStart + ncmpiPartCount[domainId()] <= ncmpiPartCountMax,
           "ERROR: Invalid number of particles " + to_string(ncmpiStart + ncmpiPartCount[domainId()]) + " > "
               + to_string(ncmpiPartCountMax));

    MLongScratchSpace ncmpiPartId(ncmpiPartCount[domainId()], AT_, "ncmpiPartId");
    MIntScratchSpace ncmpiPartParceledNo(ncmpiPartCount[domainId()], AT_, "ncmpiPartParceledNo");
    MFloatScratchSpace ncmpiDiameter(ncmpiPartCount[domainId()], AT_, "ncmpiDiameter");
    MFloatScratchSpace ncmpiTemp(ncmpiPartCount[domainId()], AT_, "ncmpiTemp");
    MFloatScratchSpace ncmpiPartCoords(3 * ncmpiPartCount[domainId()], AT_, "ncmpiPartCoords");
    //    MFloatScratchSpace ncmpiPartCoordsM(3 * ncmpiPartCount[domainId()], AT_,
    //    "ncmpiPartCoordsM");
    MFloatScratchSpace ncmpiPartVel(3 * ncmpiPartCount[domainId()], AT_, "ncmpiPartVel");
    MInt ncmpiCountId = 0;

    // Put all particle data in arrays.
    while(!ncmpiPartQueue.empty()) {
      ncmpiPartId[ncmpiCountId] = (*(ncmpiPartQueue.front())).m_partId;
      ncmpiPartParceledNo[ncmpiCountId] = (*(ncmpiPartQueue.front())).m_noParticles;
      ncmpiDiameter[ncmpiCountId] = (*(ncmpiPartQueue.front())).m_diameter;
      ncmpiTemp[ncmpiCountId] = (*(ncmpiPartQueue.front())).m_temperature;
      ncmpiPartCoords[(3 * ncmpiCountId)] = (*(ncmpiPartQueue.front())).m_position[0];
      ncmpiPartCoords[(3 * ncmpiCountId) + 1] = (*(ncmpiPartQueue.front())).m_position[1];
      ncmpiPartCoords[(3 * ncmpiCountId) + 2] = (*(ncmpiPartQueue.front())).m_position[2];
      ncmpiPartVel[(3 * ncmpiCountId)] = (*(ncmpiPartQueue.front())).m_velocity[0];
      ncmpiPartVel[(3 * ncmpiCountId) + 1] = (*(ncmpiPartQueue.front())).m_velocity[1];
      ncmpiPartVel[(3 * ncmpiCountId) + 2] = (*(ncmpiPartQueue.front())).m_velocity[2];

      ncmpiPartQueue.pop();
      ++ncmpiCountId;
    }

    // Put the arrays in the Netcdf file.
    parallelIo.setOffset(ncmpiCount, ncmpiStart);
    parallelIo.writeArray(ncmpiPartId.begin(), "partId");
    parallelIo.writeArray(ncmpiDiameter.begin(), "partDia");

    if(m_activePrimaryBUp || m_activeSecondaryBUp || m_wallCollisions) {
      parallelIo.writeArray(ncmpiPartParceledNo.begin(), "partParceledNo");
    }

    if(m_heatCoupling || m_evaporation) {
      parallelIo.writeArray(ncmpiTemp.begin(), "partTemp");
    }
    if(m_domainIdOutput) {
      MFloatScratchSpace ncmpiDomain(ncmpiPartCount[domainId()], AT_, "ncmpiDomain");
      ncmpiDomain.fill(domainId());
      parallelIo.writeArray(ncmpiDomain.begin(), "domainId");
    }

    ncmpiCount *= 3;
    ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;
    parallelIo.setOffset(ncmpiCount, ncmpi3Start);
    parallelIo.writeArray(ncmpiPartCoords.begin(), "partPos");
    parallelIo.writeArray(ncmpiPartVel.begin(), "partVel");

    if(domainId() == 0) {
      cerr << "Writing particle solution file at time step " << globalTimeStep << endl;
    }

  } else {
    if(domainId() == 0) {
      m_log << "WARNING: Write called but no particles present!" << endl;
      cerr << "WARNING: Write called but no particles present!" << endl;
    }
  }
  if(m_ellipsoids) {
    // Make ncmpiPartQueue, Queue of "Real" Particle and ncmpiPartQueueSize (MInt),
    // size of the queue.
    queue<ellipsListIterator<nDim>> ncmpiEllipsoidQueue;
    ncmpiPartQueueSize = 0;

    auto i2 = m_partListEllipsoid.begin();
    while(i2 != m_partListEllipsoid.end()) {
      if(!(*i2).isInvalid() && ((*i2).m_position[0] > m_xCutOff)) {
        ++ncmpiPartQueueSize;
        ncmpiEllipsoidQueue.push(i2);
      }
      i2++;
    }

    // Get Variable partCount[noDomains()]
    MIntScratchSpace ncmpiEllipsoidCount(noDomains(), AT_, "ncmpiEllipsoidCount");
    MPI_Allgather(&ncmpiPartQueueSize, 1, MPI_INT, &ncmpiEllipsoidCount[0], 1, MPI_INT, mpiComm(), AT_,
                  "ncmpiPartQueueSize", "ncmpiEllipsoidCount");

    // Calculate ncmpiPartCountMax
    ncmpiPartCountMax = 0;
    for(MInt i = 0; i < noDomains(); ++i) {
      ncmpiPartCountMax += ncmpiEllipsoidCount[i];
    }

    // If there are 0 particle in the whole domain, don't write particle data.
    if(ncmpiPartCountMax != 0) {
      const MString ncmpiFileName =
          outputDir() + "partEllipsoid_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

      using namespace maia::parallel_io;
      ParallelIo parallelIo(ncmpiFileName, PIO_REPLACE, mpiComm());

      // Define Attribute timestep (double).
      parallelIo.setAttribute(globalTimeStep + 1, "globalTimestep");
      parallelIo.setAttribute(globalTimeStep, "particleTimestep");
      parallelIo.setAttribute(m_time, "time");

      parallelIo.defineArray(PIO_INT, "partCount", noDomains());
      parallelIo.defineArray(PIO_LONG, "partId", ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partSemiMinorAxis", ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partAspectRatio", ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partDia", ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partDensityRatio", ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partPos", 3 * ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partVel", 3 * ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partAngVel", 3 * ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partMajorAxis", 3 * ncmpiPartCountMax);
      parallelIo.defineArray(PIO_FLOAT, "partQuat", 4 * ncmpiPartCountMax);
      // EndDef (Go into data mode).

      ParallelIo::size_type ncmpiStart = domainId();
      ParallelIo::size_type ncmpiCount = 1;

      parallelIo.setOffset(ncmpiCount, ncmpiStart);

      // Put Variable partCount = ncmpiPartQueueSize at Position MPI_Rank
      parallelIo.writeArray(&ncmpiPartQueueSize, "partCount");

      // Create arrays to hold the particle data.
      ncmpiStart = 0;

      for(MInt i = 0; i < domainId(); ++i) {
        ncmpiStart += ncmpiEllipsoidCount[i];
      }
      ncmpiCount = ncmpiEllipsoidCount[domainId()];
      if(ncmpiStart >= ncmpiPartCountMax) {
        if(ncmpiCount == 0) {
          ncmpiStart = 0;
        } else {
          mTerm(1, AT_, "Error in m_ellipsoids ncmpiStart >= ncmpiPartCountMax but ncmpiCount != 0");
        }
      }

      MLongScratchSpace ncmpiPartId(ncmpiPartCountMax, AT_, "ncmpiEllipsId");
      MFloatScratchSpace ncmpiSemiMinorAxis(ncmpiPartCountMax, AT_, "ncmpiSemiMinorAxis");
      MFloatScratchSpace ncmpiAspectRatio(ncmpiPartCountMax, AT_, "ncmpiAspectRatio");
      MFloatScratchSpace ncmpiEqDia(ncmpiPartCountMax, AT_, "ncmpiEquivalentDiameter");
      MFloatScratchSpace ncmpiDensityRatio(ncmpiPartCountMax, AT_, "ncmpiDensityRatio");
      MFloatScratchSpace ncmpiPartCoords(3 * ncmpiPartCountMax, AT_, "ncmpiCoords");
      MFloatScratchSpace ncmpiPartVel(3 * ncmpiPartCountMax, AT_, "ncmpiEllipsVel");
      MFloatScratchSpace ncmpiPartAngVel(3 * ncmpiPartCountMax, AT_, "ncmpiEllipsAngVel");
      MFloatScratchSpace ncmpiPartMajorAxis(3 * ncmpiPartCountMax, AT_, "ncmpiEllipsMajorAxis");
      MFloatScratchSpace ncmpiPartQuat(4 * ncmpiPartCountMax, AT_, "ncmpiEllipsQuat");
      MInt ncmpiCountId = 0;

      // Put all particle data in arrays.
      while(!ncmpiEllipsoidQueue.empty()) {
        LPTEllipsoidal<nDim>& particle = *ncmpiEllipsoidQueue.front();
        ncmpiPartId[ncmpiCountId] = particle.m_partId;
        ncmpiSemiMinorAxis[ncmpiCountId] = particle.m_semiMinorAxis;
        ncmpiAspectRatio[ncmpiCountId] = particle.m_aspectRatio;
        ncmpiEqDia[ncmpiCountId] = particle.equivalentDiameter();
        ncmpiDensityRatio[ncmpiCountId] = particle.m_densityRatio;
        std::array<MFloat, nDim> majorAxis{};
        particle.calculateMajorAxisOrientation(majorAxis.begin());
        for(MInt n = 0; n < nDim; n++) {
          ncmpiPartCoords[(3 * ncmpiCountId) + n] = particle.m_position[n];
          ncmpiPartVel[(3 * ncmpiCountId) + n] = particle.m_velocity[n];
          ncmpiPartAngVel[(3 * ncmpiCountId) + n] = particle.m_angularVel[n];
          ncmpiPartMajorAxis[(3 * ncmpiCountId) + n] = majorAxis[n];
        }
        for(MInt n = 0; n < 4; n++) {
          ncmpiPartQuat[(4 * ncmpiCountId) + n] = particle.m_quaternion[n];
        }

        ncmpiEllipsoidQueue.pop();
        ++ncmpiCountId;
      }

      // Put the arrays in the Netcdf file.
      parallelIo.setOffset(ncmpiCount, ncmpiStart);
      parallelIo.writeArray(&ncmpiPartId[0], "partId");
      parallelIo.writeArray(&ncmpiSemiMinorAxis[0], "partSemiMinorAxis");
      parallelIo.writeArray(&ncmpiAspectRatio[0], "partAspectRatio");
      parallelIo.writeArray(&ncmpiEqDia[0], "partDia");
      parallelIo.writeArray(&ncmpiDensityRatio[0], "partDensityRatio");

      ncmpiCount *= 3;
      ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;
      parallelIo.setOffset(ncmpiCount, ncmpi3Start);
      parallelIo.writeArray(&ncmpiPartCoords[0], "partPos");
      parallelIo.writeArray(&ncmpiPartVel[0], "partVel");
      parallelIo.writeArray(&ncmpiPartAngVel[0], "partAngVel");
      parallelIo.writeArray(&ncmpiPartMajorAxis[0], "partMajorAxis");

      ParallelIo::size_type ncmpi4Start = 4 * ncmpiStart;
      ncmpiCount /= 3;
      ncmpiCount *= 4;
      parallelIo.setOffset(ncmpiCount, ncmpi4Start);
      parallelIo.writeArray(&ncmpiPartQuat[0], "partQuat");

      // Free Memory
      while(!ncmpiEllipsoidQueue.empty()) {
        ncmpiEllipsoidQueue.pop();
      }
    }
  }
}

/**
 * \brief  Write particle restart file
 *
 * @author Jerry Grimmen, Sven Berger
 */
template <MInt nDim>
void LPT<nDim>::writeParticleRestartFile() {
  TRACE();

  // if(grid().hasInactiveRanks()) {
  grid().updatePartitionCellOffsets();
  //}

  m_log << "Writing LPT restart file..." << endl;

  // Write a collective Netcdf restart file.
  // First, create a new collective Netcdf file. (Leaves file open and in define mode).
  MInt noParticles = 0;
  for(MInt id = 0; id < a_noParticles(); id++) {
    if(m_partList[id].isInvalid()) {
      cerr << "Invalid particle when writing restart-file!" << endl;
      continue;
    }
    noParticles++;
  }
  MIntScratchSpace ncmpiPartCount(noDomains(), AT_, "ncmpiPartCount");
  const MInt noLocalPartitionCells =
      grid().localPartitionCellOffsetsRestart(1) - grid().localPartitionCellOffsetsRestart(0);
  const MInt noGlobalPartitionCells = grid().localPartitionCellOffsetsRestart(2);

  MIntScratchSpace noParticlesPerLocalPartitionCell(noLocalPartitionCells, AT_, "noParticlesPerLocalPartitionCell");
  noParticlesPerLocalPartitionCell.fill(0);

  MPI_Allgather(&noParticles, 1, MPI_INT, &ncmpiPartCount[0], 1, MPI_INT, mpiComm(), AT_, "noParticles",
                "ncmpiPartCount");

  // Calculate global number of particles
  ParallelIo::size_type globalNoParticles = 0;
  for(MInt i = 0; i < noDomains(); ++i) {
    globalNoParticles += ncmpiPartCount[i];
  }

  if(domainId() == 0) {
    cerr << "Write Particle Restart at Time step " << globalTimeStep << endl;
    cerr << "for number of particles: " << globalNoParticles << " and " << noGlobalPartitionCells << " partition cells!"
         << endl;
  }

  m_log << "Time step " << globalTimeStep << " -- this proc. has " << noParticles << " particles of "
        << globalNoParticles << endl;


  ParallelIo::size_type ncmpiStart = grid().localPartitionCellOffsetsRestart(0);
  ParallelIo::size_type ncmpiCount = noLocalPartitionCells;

  auto sortByGId = [&](const LPTBase<nDim>& i, const LPTBase<nDim>& j) {
    const MLong globalId1 = c_globalId(i.m_cellId);
    const MLong globalId2 = c_globalId(j.m_cellId);
    if(globalId1 != globalId2) {
      return (c_globalId(i.m_cellId) < c_globalId(j.m_cellId));
    } else {
      return (i.m_partId < j.m_partId);
    }
  };

  if(a_noParticles() > 0) {
    own_sort(m_partList, sortByGId);

    MInt localPartitionCellCounter = 0;
    for(auto i1 = m_partList.begin(); i1 != m_partList.end(); i1++) {
      MLong particleGlobalCellId = c_globalId(i1->m_cellId);
      if(i1->isInvalid()) continue;
      if(a_isHalo(i1->m_cellId)) {
        cerr << "Particle in halo-cell!" << endl;
        cerr << i1->hadWallColl() << " " << i1->isInvalid() << " " << i1->reqSend() << " " << i1->m_partId << endl;
        MInt cellId2 = grid().findContainingLeafCell(&i1->m_position[0]);
        MInt cellId3 = grid().findContainingLeafCell(&i1->m_position[0], i1->m_cellId, true);
        cerr << i1->m_cellId << " " << cellId2 << " " << cellId3 << endl;
      }
      if((localPartitionCellCounter + 1 < noLocalPartitionCells)) {
        ASSERT(grid().localPartitionCellGlobalIdsRestart(localPartitionCellCounter + 1) > -1, "");
        while(particleGlobalCellId >= grid().localPartitionCellGlobalIdsRestart(localPartitionCellCounter + 1)) {
          localPartitionCellCounter++;
          if(localPartitionCellCounter + 1 >= noLocalPartitionCells) {
            break;
          }
        }
      }
      noParticlesPerLocalPartitionCell[localPartitionCellCounter] += 1;
    }
  }

  MString ncmpiFileName =
      outputDir() + "restartPart_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

  using namespace maia::parallel_io;
  ParallelIo parallelIo(ncmpiFileName, PIO_REPLACE, mpiComm());

  parallelIo.setAttribute(globalTimeStep, "timestep");
  parallelIo.setAttribute(globalTimeStep, "particleTimestep");
  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    MInt globalInjStep = m_sprayModel->m_injStep;
    MPI_Allreduce(MPI_IN_PLACE, &globalInjStep, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "globalInjStep");
    parallelIo.setAttribute(globalInjStep, "injStep");

    MFloat globalTimeSOI = m_sprayModel->timeSinceSOI();
    MPI_Allreduce(MPI_IN_PLACE, &globalTimeSOI, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "INPLACE", "globalTimeSOI");
    parallelIo.setAttribute(globalTimeSOI, "timeSinceSOI");
  }
  parallelIo.setAttribute(m_time, "time");

  // count the number of times random numbers were drawn for predictable restart
  MInt count = m_PRNGSpawnCount;
  MFloat particleResiduum = 0;

  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    count = m_sprayModel->m_PRNGPrimBreakUpCount;
    particleResiduum = std::numeric_limits<MFloat>::max();
  }

  if(m_spawnParticles) {
    particleResiduum = std::numeric_limits<MFloat>::max();
  }

  if(domainId() != m_spawnDomainId) {
    ASSERT(count == 0, "");
  } else {
#ifndef NDEBUG
    // re-compute the number of times the PRNG has been used!
    mt19937_64 PRNG;
    MInt j = 0;
    if(m_spawnParticles) {
      PRNG.seed(m_spawnSeed);
    } else if(m_sprayModel) {
      PRNG.seed(m_sprayModel->m_spraySeed);
    }
    const MInt maxCount = 1000000000;
    for(MInt i = 0; i < maxCount; i++) {
      if(!m_activePrimaryBUp && !m_activeSecondaryBUp && PRNG == m_PRNGSpawn) {
        break;
      } else if((m_activePrimaryBUp || m_activeSecondaryBUp) && PRNG == m_sprayModel->m_PRNGPrimBreakUp) {
        break;
      }
      j++;
      PRNG();
    }
    if(j < maxCount) {
      ASSERT(m_PRNGSpawnCount == j
                 || ((m_activePrimaryBUp || m_activeSecondaryBUp) && j == m_sprayModel->m_PRNGPrimBreakUpCount),
             "");
    } else {
      cerr << "PRNG-state could not be checked" << endl;
    }
#endif
  }
  if(m_spawnCellId > -1) {
    particleResiduum = m_particleResiduum;
  }

  MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "count");
  MPI_Allreduce(MPI_IN_PLACE, &particleResiduum, 1, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "INPLACE", "particleResiduum");

  parallelIo.setAttribute(count, "spawnCount");
  parallelIo.setAttribute(particleResiduum, "particleResiduum");

  parallelIo.defineArray(PIO_INT, "partCount", noGlobalPartitionCells);

  if(globalNoParticles > 0) {
    parallelIo.defineArray(PIO_LONG, "partId", globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "partDia", globalNoParticles);
    parallelIo.defineArray(PIO_INT, "partStatus", globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "partPos", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "partVel", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "partAccel", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "oldPos", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "oldVel", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "oldAccel", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_INT, "partParceledNo", globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "creationTime", globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "oldFluidVelocity", 3 * globalNoParticles);
    parallelIo.defineArray(PIO_FLOAT, "oldFluidDensity", globalNoParticles);
    if(m_activePrimaryBUp || m_activeSecondaryBUp) {
      parallelIo.defineArray(PIO_FLOAT, "breakUpTime", globalNoParticles);
    }
    if(m_heatCoupling || m_evaporation) {
      parallelIo.defineArray(PIO_FLOAT, "partTemp", globalNoParticles);
      parallelIo.defineArray(PIO_FLOAT, "partDM", globalNoParticles);
    }
    if(m_activeSecondaryBUp) {
      parallelIo.defineArray(PIO_FLOAT, "shedD", globalNoParticles);
    }
  }

  // Put Variable partCount = ncmpiPartQueueSize at Position MPI_Rank
  parallelIo.setOffset(ncmpiCount, ncmpiStart);
  parallelIo.writeArray(&noParticlesPerLocalPartitionCell[0], "partCount");

  // Write particle data.
  ncmpiStart = 0;
  for(MInt i = 0; i < domainId(); ++i) {
    ncmpiStart += ncmpiPartCount[i];
  }
  ncmpiCount = ncmpiPartCount[domainId()];
  if(ncmpiStart >= globalNoParticles) {
    if(ncmpiCount == 0) {
      ncmpiStart = 0;
    } else {
      mTerm(1, AT_, "Error ncmpiStart >= globalNoParticles but ncmpiCount != 0");
    }
  }
  ASSERT(ncmpiCount == a_noParticles(), "");

  if(globalNoParticles > 0) {
    MLongScratchSpace ncmpiPartId(ncmpiCount, AT_, "ncmpiPartId");
    MFloatScratchSpace ncmpiDiameter(ncmpiCount, AT_, "ncmpiDiameter");
    MFloatScratchSpace ncmpiDensityRatio(ncmpiCount, AT_, "ncmpiDensityRatio");
    MIntScratchSpace ncmpiPartStatus(ncmpiCount, AT_, "ncmpiPartStatus");
    MFloatScratchSpace ncmpiPartCoords(3 * ncmpiCount, AT_, "PartCoords");
    MFloatScratchSpace ncmpiPartVel(3 * ncmpiCount, AT_, "ncmpiPartVel");
    MFloatScratchSpace ncmpiPartAccel(3 * ncmpiCount, AT_, "ncmpiPartAccel");
    MFloatScratchSpace ncmpiOldCoords(3 * ncmpiCount, AT_, "ncmpiOldCoords");
    MFloatScratchSpace ncmpiOldVel(3 * ncmpiCount, AT_, "ncmpiOldVel");
    MFloatScratchSpace ncmpiOldAccel(3 * ncmpiCount, AT_, "ncmpiOldAccel");
    MIntScratchSpace ncmpiPartParceledNo(ncmpiCount, AT_, "ncmpiPartParceledNo");
    MFloatScratchSpace ncmpiTemp(ncmpiCount, AT_, "ncmpiTemp");
    MFloatScratchSpace ncmpiDM(ncmpiCount, AT_, "ncmpiDM");
    MFloatScratchSpace ncmpiCT(ncmpiCount, AT_, "ncmpiCT");
    MFloatScratchSpace ncmpiBT(ncmpiCount, AT_, "ncmpiBT");
    MFloatScratchSpace ncmpiOldFluidDensity(ncmpiCount, AT_, "ncmpiOldFluidDensity");
    MFloatScratchSpace ncmpiOldFluidVelocity(3 * ncmpiCount, AT_, "ncmpiOldFluidVelocity");
    MFloatScratchSpace ncmpishedD(ncmpiCount, AT_, "ncmpishedD");

    MInt id = 0;

    for(MInt i = 0; i < a_noParticles(); i++) {
      if(m_partList[i].isInvalid()) continue;
      ncmpiPartId[id] = m_partList[i].m_partId;
      ncmpiDiameter[id] = m_partList[i].m_diameter;
      ncmpiDensityRatio[id] = m_partList[i].m_densityRatio;
      ncmpiPartStatus[id] = (m_partList[i].firstStep() ? 1 : 0);
      ncmpiPartCoords[(3 * id)] = m_partList[i].m_position[0];
      ncmpiPartCoords[(3 * id) + 1] = m_partList[i].m_position[1];
      ncmpiPartCoords[(3 * id) + 2] = m_partList[i].m_position[2];
      ncmpiPartVel[(3 * id)] = m_partList[i].m_velocity[0];
      ncmpiPartVel[(3 * id) + 1] = m_partList[i].m_velocity[1];
      ncmpiPartVel[(3 * id) + 2] = m_partList[i].m_velocity[2];
      ncmpiPartAccel[(3 * id)] = m_partList[i].m_accel[0];
      ncmpiPartAccel[(3 * id) + 1] = m_partList[i].m_accel[1];
      ncmpiPartAccel[(3 * id) + 2] = m_partList[i].m_accel[2];
      ncmpiOldCoords[(3 * id)] = m_partList[i].m_oldPos[0];
      ncmpiOldCoords[(3 * id) + 1] = m_partList[i].m_oldPos[1];
      ncmpiOldCoords[(3 * id) + 2] = m_partList[i].m_oldPos[2];
      ncmpiOldVel[(3 * id)] = m_partList[i].m_oldVel[0];
      ncmpiOldVel[(3 * id) + 1] = m_partList[i].m_oldVel[1];
      ncmpiOldVel[(3 * id) + 2] = m_partList[i].m_oldVel[2];
      ncmpiOldAccel[(3 * id)] = m_partList[i].m_oldAccel[0];
      ncmpiOldAccel[(3 * id) + 1] = m_partList[i].m_oldAccel[1];
      ncmpiOldAccel[(3 * id) + 2] = m_partList[i].m_oldAccel[2];
      ncmpiPartParceledNo[id] = m_partList[i].m_noParticles;
      ncmpiTemp[id] = m_partList[i].m_temperature;
      ncmpiDM[id] = m_partList[i].m_dM;
      ncmpiCT[id] = 0.0; // m_partList[i].m_creationTime;
      ncmpiBT[id] = m_partList[i].m_breakUpTime;
      if(m_activeSecondaryBUp) {
        ncmpishedD[id] = m_partList[i].m_shedDiam;
      }

      ncmpiOldFluidDensity[id] = m_partList[i].m_oldFluidDensity;
      ncmpiOldFluidVelocity[(3 * id)] = m_partList[i].m_oldFluidVel[0];
      ncmpiOldFluidVelocity[(3 * id) + 1] = m_partList[i].m_oldFluidVel[1];
      ncmpiOldFluidVelocity[(3 * id) + 2] = m_partList[i].m_oldFluidVel[2];

      id = id + 1;
    }

    parallelIo.setOffset(ncmpiCount, ncmpiStart);

    parallelIo.writeArray(ncmpiPartId.begin(), "partId");
    parallelIo.writeArray(ncmpiDiameter.begin(), "partDia");
    parallelIo.writeArray(ncmpiPartStatus.begin(), "partStatus");
    parallelIo.writeArray(ncmpiPartParceledNo.begin(), "partParceledNo");
    parallelIo.writeArray(ncmpiCT.begin(), "creationTime");

    if(m_activePrimaryBUp || m_activeSecondaryBUp) {
      parallelIo.writeArray(ncmpiBT.begin(), "breakUpTime");
    }

    if(m_heatCoupling || m_evaporation) {
      parallelIo.writeArray(ncmpiTemp.begin(), "partTemp");
      parallelIo.writeArray(ncmpiDM.begin(), "partDM");
    }

    if(m_activeSecondaryBUp) {
      parallelIo.writeArray(ncmpishedD.begin(), "shedD");
    }

    parallelIo.writeArray(ncmpiOldFluidDensity.begin(), "oldFluidDensity");

    ncmpiCount *= 3;
    ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;
    parallelIo.setOffset(ncmpiCount, ncmpi3Start);
    parallelIo.writeArray(ncmpiPartCoords.begin(), "partPos");
    parallelIo.writeArray(ncmpiPartVel.begin(), "partVel");
    parallelIo.writeArray(ncmpiPartAccel.begin(), "partAccel");
    parallelIo.writeArray(ncmpiOldCoords.begin(), "oldPos");
    parallelIo.writeArray(ncmpiOldVel.begin(), "oldVel");
    parallelIo.writeArray(ncmpiOldAccel.begin(), "oldAccel");
    parallelIo.writeArray(ncmpiOldFluidVelocity.begin(), "oldFluidVelocity");
  }
  // recover original sorting
  own_sort(m_partList, sort_particleAfterPartIds<LPTBase<nDim>>::compare);

  if(m_ellipsoids) {
    MInt noEllipsoids = 0;
    for(MInt id = 0; id < a_noEllipsoidalParticles(); id++) {
      if(m_partListEllipsoid[id].isInvalid()) {
        cerr << "Invalid particle when writing restart-file!" << endl;
        continue;
      }
      noEllipsoids++;
    }

    MIntScratchSpace noEllipsoidsPerLocalPartitionCell(noLocalPartitionCells, AT_, "noEllipsoidsPerLocalPartitionCell");
    noEllipsoidsPerLocalPartitionCell.fill(0);
    MIntScratchSpace ncmpiEllipsoidCount(noDomains(), AT_, "ncmpiEllipsoidCount");
    MPI_Allgather(&noEllipsoids, 1, MPI_INT, &ncmpiEllipsoidCount[0], 1, MPI_INT, mpiComm(), AT_, "noEllipsoids",
                  "ncmpiEllipsoidCount");

    // Calculate global no ellipsoids
    ParallelIo::size_type globalNoEllipsoids = 0;
    for(MInt i = 0; i < noDomains(); ++i) {
      globalNoEllipsoids += ncmpiEllipsoidCount[i];
    }

    m_log << "Time step " << globalTimeStep << " -- this proc. has " << noEllipsoids << " ellipsoids of "
          << globalNoEllipsoids << endl;
    if(domainId() == 0) {
      cerr << "Write Ellipsoid Restart at Time step " << globalTimeStep << endl;
      cerr << "for number of ellipsoids: " << globalNoEllipsoids << " and " << noGlobalPartitionCells
           << " partition cells!" << endl;
    }

    ncmpiStart = grid().localPartitionCellOffsetsRestart(0);
    ncmpiCount = noLocalPartitionCells;

    // If there are 0 particle in the whole domain, don't write particle data.
    if(a_noEllipsoidalParticles() > 0) {
      own_sort(m_partListEllipsoid, sortByGId);

      MInt localPartitionCellCounter = 0;
      for(auto i1 = m_partListEllipsoid.begin(); i1 != m_partListEllipsoid.end(); i1++) {
        MLong particleGlobalCellId = c_globalId(i1->m_cellId);
        if(i1->isInvalid()) continue;
        if(a_isHalo(i1->m_cellId)) {
          cerr << "Particle in halo-cell!" << endl;
          cerr << i1->hadWallColl() << " " << i1->isInvalid() << " " << i1->reqSend() << " " << i1->m_partId << endl;
          MInt cellId2 = grid().findContainingLeafCell(&i1->m_position[0]);
          MInt cellId3 = grid().findContainingLeafCell(&i1->m_position[0], i1->m_cellId, true);
          cerr << i1->m_cellId << " " << cellId2 << " " << cellId3 << endl;
        }
        if((localPartitionCellCounter + 1 < noLocalPartitionCells)) {
          ASSERT(grid().localPartitionCellGlobalIdsRestart(localPartitionCellCounter + 1) > -1, "");
          while(particleGlobalCellId >= grid().localPartitionCellGlobalIdsRestart(localPartitionCellCounter + 1)) {
            localPartitionCellCounter++;
            if(localPartitionCellCounter + 1 >= noLocalPartitionCells) {
              break;
            }
          }
        }
        noEllipsoidsPerLocalPartitionCell[localPartitionCellCounter] += 1;
      }
    }

    ncmpiFileName =
        outputDir() + "restartPartEllipsoid_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

    ParallelIo parallelIoE(ncmpiFileName, PIO_REPLACE, mpiComm());

    // Define Attribute timestep (double).
    parallelIoE.setAttribute(globalTimeStep, "timestep");
    parallelIoE.setAttribute(globalTimeStep, "particleTimestep");
    parallelIoE.defineArray(PIO_INT, "partCount", noGlobalPartitionCells);


    if(globalNoEllipsoids > 0) {
      parallelIoE.defineArray(PIO_LONG, "partId", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partSemiMinorAxis", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partDia", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partAspectRatio", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partDensityRatio", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_INT, "partStatus", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partPos", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partVel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partAccel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partAngVel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partAngAccel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldPos", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldVel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldAccel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldAngVel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldAngAccel", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partQuat", 4 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "partOldQuat", 4 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "creationTime", globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "oldFluidVelocity", 3 * globalNoEllipsoids);
      parallelIoE.defineArray(PIO_FLOAT, "oldFluidDensity", globalNoEllipsoids);
      if(m_heatCoupling) parallelIoE.defineArray(PIO_FLOAT, "partTemp", globalNoEllipsoids);
      // EndDef (Go into data mode).

      // Put Variable partCount = ncmpiPartQueueSize at Position MPI_Rank
      parallelIoE.setOffset(ncmpiCount, ncmpiStart);
      parallelIoE.writeArray(&noEllipsoidsPerLocalPartitionCell[0], "partCount");

      // Write particle data.
      ncmpiStart = 0;
      for(MInt i = 0; i < domainId(); ++i) {
        ncmpiStart += ncmpiEllipsoidCount[i];
      }
      ncmpiCount = ncmpiEllipsoidCount[domainId()];
      if(ncmpiStart >= globalNoEllipsoids) {
        if(ncmpiCount == 0) {
          ncmpiStart = 0;
        } else {
          mTerm(1, AT_, "m_ellipsoids ncmpiStart >= globalNoEllipsoids but ncmpiCount != 0");
        }
      }
      ASSERT(ncmpiCount == a_noEllipsoidalParticles(), "");

      MLongScratchSpace ncmpiPartId(ncmpiCount, AT_, "ncmpiPartIdE");
      MFloatScratchSpace ncmpiSemiMinorAxis(ncmpiCount, AT_, "ncmpiSemiMinorAxis");
      MFloatScratchSpace ncmpiEquivalentDia(ncmpiCount, AT_, "ncmpiEquivalentDia");
      MFloatScratchSpace ncmpiAspectRatio(ncmpiCount, AT_, "ncmpiAspectRatio");
      MFloatScratchSpace ncmpiDensityRatio(ncmpiCount, AT_, "ncmpiDensityRatio");
      MIntScratchSpace ncmpiPartStatus(ncmpiCount, AT_, "ncmpiPartStatus");
      MFloatScratchSpace ncmpiPartCoords(3 * ncmpiCount, AT_, "ncmpiPartCoords");
      MFloatScratchSpace ncmpiPartVel(3 * ncmpiCount, AT_, "ncmpiPartVel");
      MFloatScratchSpace ncmpiPartAccel(3 * ncmpiCount, AT_, "ncmpiPartAccel");
      MFloatScratchSpace ncmpiPartAngVel(3 * ncmpiCount, AT_, "ncmpiPartAngVel");
      MFloatScratchSpace ncmpiPartAngAccel(3 * ncmpiCount, AT_, "ncmpiPartAngAccel");
      MFloatScratchSpace ncmpiPartOldCoords(3 * ncmpiCount, AT_, "ncmpiPartOldCoords");
      MFloatScratchSpace ncmpiPartOldVel(3 * ncmpiCount, AT_, "ncmpiPartOldVel");
      MFloatScratchSpace ncmpiPartOldAccel(3 * ncmpiCount, AT_, "ncmpiPartOldAccel");
      MFloatScratchSpace ncmpiPartOldAngVel(3 * ncmpiCount, AT_, "ncmpiPartOldAngVel");
      MFloatScratchSpace ncmpiPartOldAngAccel(3 * ncmpiCount, AT_, "ncmpiPartOldAngAccel");
      MFloatScratchSpace ncmpiPartQuat(4 * ncmpiCount, AT_, "ncmpiPartQuaternion");
      MFloatScratchSpace ncmpiPartOldQuat(4 * ncmpiCount, AT_, "ncmpiPartOldQuaternion");
      MFloatScratchSpace ncmpiTemp(ncmpiCount, AT_, "ncmpiTemp");
      MFloatScratchSpace ncmpiCT(ncmpiCount, AT_, "ncmpiCT");
      MFloatScratchSpace ncmpiOldFluidDensity(ncmpiCount, AT_, "ncmpiOldFluidDensity");
      MFloatScratchSpace ncmpiOldFluidVelocity(3 * ncmpiCount, AT_, "ncmpiOldFluidVelocity");

      for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
        if(m_partListEllipsoid[i].isInvalid()) continue;
        ncmpiPartId[i] = m_partListEllipsoid[i].m_partId;
        ncmpiSemiMinorAxis[i] = m_partListEllipsoid[i].m_semiMinorAxis;
        ncmpiEquivalentDia[i] = m_partListEllipsoid[i].equivalentDiameter();
        ncmpiAspectRatio[i] = m_partListEllipsoid[i].m_aspectRatio;
        ncmpiDensityRatio[i] = m_partListEllipsoid[i].m_densityRatio;
        ncmpiPartStatus[i] = (m_partListEllipsoid[i].firstStep() ? 1 : 0);
        for(MInt n = 0; n < nDim; n++) {
          ncmpiPartCoords[(3 * i) + n] = m_partListEllipsoid[i].m_position[n];
          ncmpiPartVel[(3 * i) + n] = m_partListEllipsoid[i].m_velocity[n];
          ncmpiPartAccel[(3 * i) + n] = m_partListEllipsoid[i].m_accel[n];
          ncmpiPartAngVel[(3 * i) + n] = m_partListEllipsoid[i].m_angularVel[n];
          ncmpiPartAngAccel[(3 * i) + n] = m_partListEllipsoid[i].m_angularAccel[n];
          ncmpiPartOldCoords[(3 * i) + n] = m_partListEllipsoid[i].m_oldPos[n];
          ncmpiPartOldVel[(3 * i) + n] = m_partListEllipsoid[i].m_oldVel[n];
          ncmpiPartOldAccel[(3 * i) + n] = m_partListEllipsoid[i].m_oldAccel[n];
          ncmpiPartOldAngVel[(3 * i) + n] = m_partListEllipsoid[i].m_oldAngularVel[n];
          ncmpiPartOldAngAccel[(3 * i) + n] = m_partListEllipsoid[i].m_oldAngularAccel[n];
          ncmpiOldFluidVelocity[(3 * i) + n] = m_partListEllipsoid[i].m_oldFluidVel[n];
        }
        for(MInt n = 0; n < 4; n++) {
          ncmpiPartQuat[(4 * i) + n] = m_partListEllipsoid[i].m_quaternion[n];
          ncmpiPartOldQuat[(4 * i) + n] = m_partListEllipsoid[i].m_oldQuaternion[n];
        }
        ncmpiOldFluidDensity[i] = m_partListEllipsoid[i].m_oldFluidDensity;
        ncmpiTemp[i] = m_partListEllipsoid[i].m_temperature;
        ncmpiCT[i] = m_partListEllipsoid[i].m_creationTime;
      }

      parallelIoE.setOffset(ncmpiCount, ncmpiStart);
      parallelIoE.writeArray(ncmpiPartId.begin(), "partId");
      parallelIoE.writeArray(ncmpiSemiMinorAxis.begin(), "partSemiMinorAxis");
      parallelIoE.writeArray(ncmpiEquivalentDia.begin(), "partDia");
      parallelIoE.writeArray(ncmpiAspectRatio.begin(), "partAspectRatio");
      parallelIoE.writeArray(ncmpiDensityRatio.begin(), "partDensityRatio");
      parallelIoE.writeArray(ncmpiPartStatus.begin(), "partStatus");
      parallelIoE.writeArray(ncmpiCT.begin(), "creationTime");
      parallelIoE.writeArray(ncmpiOldFluidDensity.begin(), "oldFluidDensity");
      if(m_heatCoupling) parallelIoE.writeArray(ncmpiTemp.begin(), "partTemp");

      ncmpiCount *= 3;
      ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;
      parallelIoE.setOffset(ncmpiCount, ncmpi3Start);
      parallelIoE.writeArray(ncmpiPartCoords.begin(), "partPos");
      parallelIoE.writeArray(ncmpiPartVel.begin(), "partVel");
      parallelIoE.writeArray(ncmpiPartAccel.begin(), "partAccel");
      parallelIoE.writeArray(ncmpiPartAngVel.begin(), "partAngVel");
      parallelIoE.writeArray(ncmpiPartAngAccel.begin(), "partAngAccel");
      parallelIoE.writeArray(ncmpiPartOldCoords.begin(), "partOldPos");
      parallelIoE.writeArray(ncmpiPartOldVel.begin(), "partOldVel");
      parallelIoE.writeArray(ncmpiPartOldAccel.begin(), "partOldAccel");
      parallelIoE.writeArray(ncmpiPartOldAngVel.begin(), "partOldAngVel");
      parallelIoE.writeArray(ncmpiPartOldAngAccel.begin(), "partOldAngAccel");

      ParallelIo::size_type ncmpi4Start = 4 * ncmpiStart;
      ncmpiCount /= 3;
      ncmpiCount *= 4;
      parallelIoE.setOffset(ncmpiCount, ncmpi4Start);
      parallelIoE.writeArray(ncmpiPartQuat.begin(), "partQuat");
      parallelIoE.writeArray(ncmpiPartOldQuat.begin(), "partOldQuat");
    }
    own_sort(m_partListEllipsoid, sort_particleAfterPartIds<LPTBase<nDim>>::compare);
  }
}

/**
 * \brief  Read particle restart file
 *         and create new particle instances accordingly
 *
 * @author Jerry Grimmen, Apr-10
 */
template <MInt nDim>
MInt LPT<nDim>::loadParticleRestartFile() {
  TRACE();

  // if(grid().hasInactiveRanks()) {
  grid().updatePartitionCellOffsets();
  //}

  using namespace maia::parallel_io;

  // Open collective Netcdf file.
  MString particleRestartFilename =
      restartDir() + "restartPart_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

  // 1. Loading spherical particles
  {
    // check if file exists
    struct stat buffer {};
    if(stat(particleRestartFilename.c_str(), &buffer) != 0) {
      if(domainId() == 0) {
        cerr << "WARNING: starting with restart but no particle file present!" << endl;
        m_log << "WARNING: starting with restart but no particle file present!" << endl;
      }
      m_restartFile = false;
      return 0;
    }

    // check that the timeStep in the restart-file matches the current globalTimeStep!
    MFloat loadedTimeStep = 0;
    ParallelIo parallelIo(particleRestartFilename, PIO_READ, mpiComm());

    parallelIo.getAttribute(&loadedTimeStep, "timestep");
    if(MInt(loadedTimeStep) != globalTimeStep) {
      stringstream errorMessage;
      errorMessage << "Error! restartTimeStep = " << globalTimeStep << " differs from timestep in restartPart"
                   << ParallelIo::fileExt() << ", which is " << loadedTimeStep << endl;
      mTerm(1, AT_, errorMessage.str());
    }

    parallelIo.getAttribute(&m_time, "time");

    if(m_activePrimaryBUp || m_activeSecondaryBUp) {
      parallelIo.getAttribute(&m_sprayModel->m_injStep, "injStep");
      parallelIo.getAttribute(&m_sprayModel->timeSinceSOI(), "timeSinceSOI");
    }


    // Get the number of particle each domain has.
    const MInt noGlobalPartitionCells = grid().localPartitionCellOffsetsRestart(2);
    const MInt localPartitionCellBeginning = grid().localPartitionCellOffsetsRestart(0);
    const MInt localPartitionCellEnd = grid().localPartitionCellOffsetsRestart(1);

    MIntScratchSpace noParticlesPerPartitionCell(noGlobalPartitionCells, AT_, "noParticlesPerPartitionCell");

    if(domainId() == 0) {
      cerr << "noGlobalPartitionCells " << noGlobalPartitionCells << endl;
    }

    // load number of particle in partition cell for all partition cells!
    parallelIo.setOffset(noGlobalPartitionCells, 0);
    parallelIo.readArray(&noParticlesPerPartitionCell[0], "partCount");

    // total number of particles
    MInt numberOfParts = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfParts)
#endif
    for(MInt i = 0; i < noGlobalPartitionCells; ++i) {
      numberOfParts += noParticlesPerPartitionCell[i];
    }
    if(domainId() == 0) {
      cerr << "Global number of particles loaded from restart file " << numberOfParts << endl;
    }
    m_log << "Global number of particles loaded from restart file " << numberOfParts << endl;


    m_PRNGSpawn.seed(m_spawnSeed);
    m_particleResiduum = 0;
    if(domainId() == m_spawnDomainId) {
      parallelIo.getAttribute(&m_PRNGSpawnCount, "spawnCount");
      m_PRNGSpawn.discard(m_PRNGSpawnCount);

      m_log << "PRNG state after restart " << randomSpawn(0) << endl;

      parallelIo.getAttribute(&m_particleResiduum, "particleResiduum");
    }

    ParallelIo::size_type beginPartId = 0;
    beginPartId = 0;
    for(MInt i = 0; i < localPartitionCellBeginning; ++i) {
      beginPartId += noParticlesPerPartitionCell[i];
    }

    MInt localNoPart = noParticlesPerPartitionCell[localPartitionCellBeginning];
    for(MInt i = (localPartitionCellBeginning + 1); i < localPartitionCellEnd; ++i) {
      localNoPart += noParticlesPerPartitionCell[i];
    }

    MInt checkPartCount = 0;
    MPI_Allreduce(&localNoPart, &checkPartCount, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "localNoPart");

    if(domainId() == 0 && checkPartCount != numberOfParts) {
      // numberOfParts is only set on the root-rank!
      TERM(-1);
    }

    ASSERT(localNoPart >= 0, "ERROR: Invalid number of particles to be loaded during restart");

    if(localNoPart == 0) {
      beginPartId = 0;
    }
    // if the last domain has nothing to load, beginPartId is
    // one index after the dimension of the array

    // Preparing arrays and getting particle datas.
    if(numberOfParts > 0) {
      MLongScratchSpace ncmpiPartId(localNoPart, AT_, "ncmpiPartId");
      MFloatScratchSpace ncmpiDiameter(localNoPart, AT_, "ncmpiDiameter");
      MFloatScratchSpace ncmpiDensityRatio(localNoPart, AT_, "ncmpiDensityRatio");
      MIntScratchSpace ncmpiPartStatus(localNoPart, AT_, "ncmpiPartStatus");
      MFloatScratchSpace ncmpiPartCoords(3 * localNoPart, AT_, "PartCoords");
      MFloatScratchSpace ncmpiPartVel(3 * localNoPart, AT_, "ncmpiPartVel");
      MFloatScratchSpace ncmpiPartAccel(3 * localNoPart, AT_, "ncmpiPartAccel");
      MFloatScratchSpace ncmpiOldCoords(3 * localNoPart, AT_, "ncmpiOldCoords");
      MFloatScratchSpace ncmpiOldVel(3 * localNoPart, AT_, "ncmpiOldVel");
      MFloatScratchSpace ncmpiOldAccel(3 * localNoPart, AT_, "ncmpiOldAccel");
      MIntScratchSpace ncmpiPartParceledNo(localNoPart, AT_, "ncmpiPartParceledNo");
      MFloatScratchSpace ncmpiTemp(localNoPart, AT_, "ncmpiTemp");
      MFloatScratchSpace ncmpiDM(localNoPart, AT_, "ncmpiDM");
      MFloatScratchSpace ncmpiCT(localNoPart, AT_, "ncmpiCT");
      MFloatScratchSpace ncmpiBT(localNoPart, AT_, "ncmpiBT");
      MFloatScratchSpace ncmpiOldFluidDensity(localNoPart, AT_, "ncmpiOldFluidDensity");
      MFloatScratchSpace ncmpiOldFluidVelocity(3 * localNoPart, AT_, "ncmpiOldFluidVelocity");
      MFloatScratchSpace ncmpiShedD(localNoPart, AT_, "ncmpiShedD");

      parallelIo.setOffset(localNoPart, beginPartId);
      parallelIo.readArray(ncmpiPartId.begin(), "partId");
      parallelIo.readArray(ncmpiDiameter.begin(), "partDia");
      // parallelIo.readArray(ncmpiDensityRatio.begin(), "partPpPf");
      parallelIo.readArray(ncmpiPartStatus.begin(), "partStatus");
      parallelIo.readArray(ncmpiPartParceledNo.begin(), "partParceledNo");
      parallelIo.readArray(ncmpiCT.begin(), "creationTime");

      if(m_activePrimaryBUp || m_activeSecondaryBUp) {
        parallelIo.readArray(ncmpiBT.begin(), "breakUpTime");
      }


      if(m_activeSecondaryBUp) {
        parallelIo.readArray(ncmpiShedD.begin(), "shedD");
      }

      if(m_heatCoupling || m_evaporation) {
        parallelIo.readArray(ncmpiTemp.begin(), "partTemp");
        parallelIo.readArray(ncmpiDM.begin(), "partDM");
      }
      parallelIo.readArray(ncmpiOldFluidDensity.begin(), "oldFluidDensity");

      ParallelIo::size_type ncmpi3Start = 3 * beginPartId;
      ParallelIo::size_type ncmpi3Count = 3 * localNoPart;
      parallelIo.setOffset(ncmpi3Count, ncmpi3Start);
      parallelIo.readArray(ncmpiPartCoords.begin(), "partPos");
      parallelIo.readArray(ncmpiPartVel.begin(), "partVel");
      parallelIo.readArray(ncmpiPartAccel.begin(), "partAccel");
      parallelIo.readArray(ncmpiOldFluidVelocity.begin(), "oldFluidVelocity");

      LPTSpherical<nDim> thisParticle;

      for(MLong i = 0; i < localNoPart; ++i) {
        thisParticle.m_position[0] = ncmpiPartCoords[(3 * i)];
        thisParticle.m_position[1] = ncmpiPartCoords[(3 * i) + 1];
        thisParticle.m_position[2] = ncmpiPartCoords[(3 * i) + 2];
        // find matching cellId
        const MInt cellId = grid().findContainingLeafCell(&thisParticle.m_position[0]);
        if(cellId == -1) continue;
        if(a_isHalo(cellId)) {
          cerr << "Particle in halo-cell!" << endl;
          continue;
        }

        thisParticle.m_cellId = cellId;
        ASSERT(thisParticle.m_cellId >= 0, "Invalid cellId! " + to_string(thisParticle.m_cellId));
        ASSERT(c_isLeafCell(thisParticle.m_cellId), "No leaf cell... " + std::to_string(thisParticle.m_cellId));

        thisParticle.m_partId = ncmpiPartId[i];
        thisParticle.m_noParticles = ncmpiPartParceledNo[i];
        thisParticle.m_temperature = ncmpiTemp[i];
        thisParticle.m_dM = ncmpiDM[i];
        thisParticle.m_creationTime = ncmpiCT[i];
        thisParticle.m_breakUpTime = ncmpiBT[i];
        thisParticle.m_diameter = ncmpiDiameter[i];
        thisParticle.m_densityRatio = ncmpiDensityRatio[i];
        thisParticle.firstStep() = (ncmpiPartStatus[i] > 0);
        thisParticle.m_velocity[0] = ncmpiPartVel[(3 * i)];
        thisParticle.m_velocity[1] = ncmpiPartVel[(3 * i) + 1];
        thisParticle.m_velocity[2] = ncmpiPartVel[(3 * i) + 2];
        thisParticle.m_accel[0] = ncmpiPartAccel[(3 * i)];
        thisParticle.m_accel[1] = ncmpiPartAccel[(3 * i) + 1];
        thisParticle.m_accel[2] = ncmpiPartAccel[(3 * i) + 2];

        thisParticle.m_oldCellId = cellId;
        for(MInt j = 0; j < nDim; j++) {
          thisParticle.m_oldPos[j] = thisParticle.m_position[j];
          thisParticle.m_oldVel[j] = thisParticle.m_velocity[j];
          thisParticle.m_oldAccel[j] = thisParticle.m_accel[j];
        }

        thisParticle.m_oldFluidDensity = ncmpiOldFluidDensity[i];
        for(MInt j = 0; j < nDim; j++) {
          thisParticle.m_oldFluidVel[j] = ncmpiOldFluidVelocity[(3 * i) + j];
        }

        thisParticle.updateProperties();
        if(m_activeSecondaryBUp) {
          thisParticle.m_shedDiam = ncmpiShedD[i];
        } else {
          thisParticle.m_shedDiam = thisParticle.m_diameter;
        }

        m_partList.push_back(thisParticle);
      }

      MLong sumPart = 0;
      MLong noPart = a_noParticles();
      MPI_Allreduce(&noPart, &sumPart, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_, "INPLACE", "sumpart");

      // check if everything is loaded
      if(domainId() == 0) {
        cerr << "Loaded " << numberOfParts << " particles of " << sumPart << endl;

        if(numberOfParts != sumPart) {
          TERMM(-1, "invalid number of particles loaded");
        }
      }
    }
  }

  if(m_ellipsoids) {
    particleRestartFilename =
        restartDir() + "restartPartEllipsoid_" + getIdentifier() + to_string(globalTimeStep) + ParallelIo::fileExt();

    ParallelIo parallelIo2(particleRestartFilename, PIO_READ, mpiComm());

    MInt loadedTimeStep2 = 0;
    parallelIo2.getAttribute(&loadedTimeStep2, "timestep");
    if(loadedTimeStep2 != globalTimeStep) {
      stringstream errorMessage;
      errorMessage << endl
                   << "Error! restartTimeStep = " << globalTimeStep
                   << " differs from timestep in "
                      "restartPartEllipsoid"
                   << ParallelIo::fileExt() << ", which is " << loadedTimeStep2 << endl;
      mTerm(1, AT_, errorMessage.str());
    }

    // Get the number of particle each domain has.
    const MInt noGlobalPartitionCells = grid().localPartitionCellOffsetsRestart(2);
    const MInt localPartitionCellBeginning = grid().localPartitionCellOffsetsRestart(0);
    const MInt localPartitionCellEnd = grid().localPartitionCellOffsetsRestart(1);

    MIntScratchSpace noEllipsoidsPerPartitionCell(noGlobalPartitionCells, AT_, "noParticlesPerPartitionCell");

    if(domainId() == 0) {
      cerr << "noGlobalPartitionCells " << noGlobalPartitionCells << endl;
    }

    // load number of particle in partition cell for all partition cells!
    parallelIo2.setOffset(noGlobalPartitionCells, 0);
    parallelIo2.readArray(&noEllipsoidsPerPartitionCell[0], "partCount");

    // total number of ellipsoids
    MInt numberOfEllips = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numberOfEllips)
#endif
    for(MInt i = 0; i < noGlobalPartitionCells; ++i) {
      numberOfEllips += noEllipsoidsPerPartitionCell[i];
    }

    if(domainId() == 0) {
      cerr << "Global number of ellipsoids loaded from restart file " << numberOfEllips << endl;
    }
    m_log << "Global number of ellipsoids loaded from restart file " << numberOfEllips << endl;

    ParallelIo::size_type beginPartId = 0;
    beginPartId = 0;
    for(MInt i = 0; i < localPartitionCellBeginning; ++i) {
      beginPartId += noEllipsoidsPerPartitionCell[i];
    }

    MInt localNoPart = noEllipsoidsPerPartitionCell[localPartitionCellBeginning];
    for(MInt i = (localPartitionCellBeginning + 1); i < localPartitionCellEnd; ++i) {
      localNoPart += noEllipsoidsPerPartitionCell[i];
    }

    MInt checkPartCount = 0;
    MPI_Allreduce(&localNoPart, &checkPartCount, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "localNoPart");

    if(domainId() == 0 && checkPartCount != numberOfEllips) {
      // numberOfEllps is only set on the root-rank!
      TERM(-1);
    }

    ASSERT(localNoPart >= 0, "ERROR: Invalid number of ellipsoids to be loaded during restart");

    if(localNoPart == 0) {
      beginPartId = 0;
    }

    // Preparing arrays and getting particle datas.
    if(numberOfEllips > 0) {
      MLongScratchSpace ncmpiPartId(localNoPart, AT_, "ncmpiPartId");
      MFloatScratchSpace ncmpiSemiMinorAxis(localNoPart, AT_, "ncmpiSemiMinorAxis");
      MFloatScratchSpace ncmpiAspectRatio(localNoPart, AT_, "ncmpiAspectRatio");
      MFloatScratchSpace ncmpiDensityRatio(localNoPart, AT_, "ncmpiDensityRatio");
      MIntScratchSpace ncmpiPartStatus(localNoPart, AT_, "ncmpiPartStatus");
      MFloatScratchSpace ncmpiPartCoords(3 * localNoPart, AT_, "ncmpiPartCoords");
      MFloatScratchSpace ncmpiPartVel(3 * localNoPart, AT_, "ncmpiPartVel");
      MFloatScratchSpace ncmpiPartAccel(3 * localNoPart, AT_, "ncmpiPartAccel");
      MFloatScratchSpace ncmpiPartAngVel(3 * localNoPart, AT_, "ncmpiPartVel");
      MFloatScratchSpace ncmpiPartAngAccel(3 * localNoPart, AT_, "ncmpiPartAngAccel");
      MFloatScratchSpace ncmpiPartOldCoords(3 * localNoPart, AT_, "ncmpiPartOldCoords");
      MFloatScratchSpace ncmpiPartOldVel(3 * localNoPart, AT_, "ncmpiPartOldVel");
      MFloatScratchSpace ncmpiPartOldAccel(3 * localNoPart, AT_, "ncmpiPartOldAccel");
      MFloatScratchSpace ncmpiPartOldAngVel(3 * localNoPart, AT_, "ncmpiPartOldAngVel");
      MFloatScratchSpace ncmpiPartOldAngAccel(3 * localNoPart, AT_, "ncmpiPartOldAngAccel");
      MFloatScratchSpace ncmpiPartQuat(4 * localNoPart, AT_, "ncmpiPartQuaternions");
      MFloatScratchSpace ncmpiPartOldQuat(4 * localNoPart, AT_, "ncmpiPartOldQuaternions");
      MFloatScratchSpace ncmpiTemp(localNoPart, AT_, "ncmpiTemp");
      MFloatScratchSpace ncmpiCT(localNoPart, AT_, "ncmpiCT");
      MFloatScratchSpace ncmpiOldFluidDensity(localNoPart, AT_, "ncmpiOldFluidDensity");
      MFloatScratchSpace ncmpiOldFluidVelocity(3 * localNoPart, AT_, "ncmpiOldFluidVelocity");

      parallelIo2.setOffset(localNoPart, beginPartId);
      parallelIo2.readArray(ncmpiPartId.begin(), "partId");
      parallelIo2.readArray(ncmpiSemiMinorAxis.begin(), "partSemiMinorAxis");
      parallelIo2.readArray(ncmpiAspectRatio.begin(), "partAspectRatio");
      parallelIo2.readArray(ncmpiDensityRatio.begin(), "partDensityRatio");
      parallelIo2.readArray(ncmpiPartStatus.begin(), "partStatus");
      parallelIo2.readArray(ncmpiCT.begin(), "creationTime");
      if(m_heatCoupling) parallelIo2.readArray(ncmpiTemp.begin(), "partTemp");
      parallelIo2.readArray(ncmpiOldFluidDensity.begin(), "oldFluidDensity");

      ParallelIo::size_type localNoPart3 = 3 * localNoPart;
      ParallelIo::size_type beginPartId3 = 3 * beginPartId;
      parallelIo2.setOffset(localNoPart3, beginPartId3);
      parallelIo2.readArray(ncmpiPartCoords.begin(), "partPos");
      parallelIo2.readArray(ncmpiPartVel.begin(), "partVel");
      parallelIo2.readArray(ncmpiPartAccel.begin(), "partAccel");
      parallelIo2.readArray(ncmpiPartAngVel.begin(), "partAngVel");
      parallelIo2.readArray(ncmpiPartAngAccel.begin(), "partAngAccel");
      parallelIo2.readArray(ncmpiPartOldCoords.begin(), "partOldPos");
      parallelIo2.readArray(ncmpiPartOldVel.begin(), "partOldVel");
      parallelIo2.readArray(ncmpiPartOldAccel.begin(), "partOldAccel");
      parallelIo2.readArray(ncmpiPartOldAngVel.begin(), "partOldAngVel");
      parallelIo2.readArray(ncmpiPartOldAngAccel.begin(), "partOldAngAccel");

      ParallelIo::size_type localNoPart4 = 4 * localNoPart;
      ParallelIo::size_type beginPartId4 = 4 * beginPartId;
      parallelIo2.setOffset(localNoPart4, beginPartId4);
      parallelIo2.readArray(ncmpiPartQuat.begin(), "partQuat");
      parallelIo2.readArray(ncmpiPartOldQuat.begin(), "partOldQuat");

      LPTEllipsoidal<nDim> thisParticleEllipsoid;

      // If any particle, write particle to the particle list m_partList.
      for(MLong i = 0; i < localNoPart; ++i) {
        for(MInt n = 0; n < nDim; n++)
          thisParticleEllipsoid.m_position[n] = ncmpiPartCoords[(3 * i) + n];
        // find matching cellId
        const MInt cellId = grid().findContainingLeafCell(&thisParticleEllipsoid.m_position[0]);
        if(cellId == -1) continue;
        if(a_isHalo(cellId)) {
          cerr << "Particle in halo-cell!" << endl;
          continue;
        }

        thisParticleEllipsoid.m_cellId = cellId;
        ASSERT(thisParticleEllipsoid.m_cellId >= 0, "Invalid cellId! " + to_string(thisParticleEllipsoid.m_cellId));
        ASSERT(c_isLeafCell(thisParticleEllipsoid.m_cellId),
               "No leaf cell... " + std::to_string(thisParticleEllipsoid.m_cellId));
        thisParticleEllipsoid.m_partId = ncmpiPartId[i];
        thisParticleEllipsoid.m_semiMinorAxis = ncmpiSemiMinorAxis[i];
        thisParticleEllipsoid.m_aspectRatio = ncmpiAspectRatio[i];
        thisParticleEllipsoid.initEllipsoialProperties();
        thisParticleEllipsoid.m_densityRatio = ncmpiDensityRatio[i];
        thisParticleEllipsoid.m_temperature = ncmpiTemp[i];
        thisParticleEllipsoid.m_creationTime = ncmpiCT[i];
        thisParticleEllipsoid.firstStep() = (ncmpiPartStatus[i] > 0);
        thisParticleEllipsoid.m_oldCellId = cellId;
        thisParticleEllipsoid.m_oldFluidDensity = ncmpiOldFluidDensity[i];

        for(MInt n = 0; n < nDim; n++) {
          // velocites and accelerations
          thisParticleEllipsoid.m_velocity[n] = ncmpiPartVel[(3 * i) + n];
          thisParticleEllipsoid.m_accel[n] = ncmpiPartAccel[(3 * i) + n];
          thisParticleEllipsoid.m_angularVel[n] = ncmpiPartAngVel[(3 * i) + n];
          thisParticleEllipsoid.m_angularAccel[n] = ncmpiPartAngAccel[(3 * i) + n];
          // old position, velocities and accelerations
          thisParticleEllipsoid.m_oldPos[n] = ncmpiPartOldCoords[(3 * i) + n];
          thisParticleEllipsoid.m_oldVel[n] = ncmpiPartOldVel[(3 * i) + n];
          thisParticleEllipsoid.m_oldAccel[n] = ncmpiPartOldAccel[(3 * i) + n];
          thisParticleEllipsoid.m_oldAngularVel[n] = ncmpiPartOldAngVel[(3 * i) + n];
          thisParticleEllipsoid.m_oldAngularAccel[n] = ncmpiPartOldAngAccel[(3 * i) + n];
          // copy fluid velocity
          thisParticleEllipsoid.m_oldFluidVel[n] = ncmpiOldFluidVelocity[(3 * i) + n];
        }
        // store current and old quaternions
        for(MInt n = 0; n < 4; n++) {
          thisParticleEllipsoid.m_quaternion[n] = ncmpiPartQuat[(4 * i) + n];
          thisParticleEllipsoid.m_oldQuaternion[n] = ncmpiPartOldQuat[(4 * i) + n];
        }
        thisParticleEllipsoid.updateProperties();

        m_partListEllipsoid.push_back(thisParticleEllipsoid);
      }
      MLong sumEllips = 0;
      MLong noEllips = a_noEllipsoidalParticles();
      MPI_Allreduce(&noEllips, &sumEllips, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_, "INPLACE",
                    "sumpart");

      // check if everything is loaded
      if(domainId() == 0) {
        cerr << "Loaded " << numberOfEllips << " particles of " << sumEllips << endl;

        if(numberOfEllips != sumEllips) {
          TERMM(-1, "invalid number of ellipsoids loaded");
        }
      }
    }
  }

  MInt theSize = a_noParticles();
  own_sort(m_partList, sort_particleAfterPartIds<LPTBase<nDim>>::compare);
  if(m_ellipsoids) {
    theSize += a_noEllipsoidalParticles();
    own_sort(m_partListEllipsoid, sort_particleAfterPartIds<LPTBase<nDim>>::compare);
  }

  m_addedParticle = 0;
  MLong noPart = a_noParticles();
  MPI_Allreduce(MPI_IN_PLACE, &noPart, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_, "INPLACE", "noPart");

  if(m_ellipsoids) {
    MLong noEllipsoids = a_noEllipsoidalParticles();
    MPI_Allreduce(MPI_IN_PLACE, &noEllipsoids, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_, "INPLACE",
                  "noPart");
    noPart += noEllipsoids;
  }

  if(domainId() == m_spawnDomainId) {
    m_addedParticle = noPart;
  }

  return theSize;
}

/**
 * \brief  Write lpt cell based solution file
 *         currently ony used for debugging!
 *
 * @author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::writeCellSolutionFile(const MString& gridFileName, MInt* recalcIds) {
  TRACE();

  MInt noCells;
  MInt noInternalCellIds;
  std::vector<MInt> recalcIdsSolver(0);
  std::vector<MInt> reOrderedCells(0);
  this->calcRecalcCellIdsSolver(recalcIds, noCells, noInternalCellIds, recalcIdsSolver, reOrderedCells);
  MInt* pointerRecalcIds = (recalcIds == nullptr) ? nullptr : recalcIdsSolver.data();

  if(grid().newMinLevel() > 0) {
    cerr0 << "Skipping LPT solution file when min-level changes are applied!" << endl;
    return;
  }

  MBool debugOutput = false;
#ifdef LPT_DEBUG
  debugOutput = true;
#endif

  const MInt noIdParams = 0;
  const MInt noDbParams = 0;
  const MInt noIdVars = 1 + m_ellipsoids;
  // VolumeFraction + flowVariables + noSpecies + massCoupling + heatCoupling + momentumCoupling
  MInt noDbVars = 1 + PV.noVars() + m_evaporation + m_massCoupling + m_heatCoupling + (nDim + 1) * m_momentumCoupling;
  if(m_ellipsoids) {
    noDbVars += (nDim * nDim); // if ellipsoids are used also write out the velocity slopes
  }
  if(debugOutput) {
    if(this->m_adaptation && this->m_noSensors > 0) {
      noDbVars = noDbVars + 1;
    }
    if(m_wallCollisions) {
      noDbVars = noDbVars + 1;
    }
  }

  MIntScratchSpace idVariables(noCells * noIdVars, AT_, "idVariables");
  MFloatScratchSpace dbVariables(noCells * noDbVars, AT_, "dbVariables");
  MIntScratchSpace idParameters(noIdParams, AT_, "idParameters");
  MFloatScratchSpace dbParameters(noDbParams, AT_, "dbParameters");
  vector<MString> dbVariablesName;
  vector<MString> idVariablesName;
  vector<MString> dbParametersName;
  vector<MString> idParametersName;
  vector<MString> name;

  // TODO labels:LPT,IO @Julian, avoid using buffer
  MIntScratchSpace tmp(noCells, AT_, "tmp");
  MFloatScratchSpace tmpW(noCells, AT_, "tmpw");

  // gather solver data:
  name.clear();
  name.push_back("noParticlesInCell");
  for(MInt cell = 0; cell < noCells; cell++) {
    tmp[cell] = a_noParticlesInCell(cell);
  }
  this->collectVariables(tmp.begin(), idVariables, name, idVariablesName, 1, noCells);

  if(m_ellipsoids) {
    name.clear();
    name.push_back("noEllipsoidsInCell");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmp[cell] = a_noEllipsoidsInCell(cell);
    }
    this->collectVariables(tmp.begin(), idVariables, name, idVariablesName, 1, noCells);
  }

  if(m_massCoupling) {
    name.clear();
    name.push_back("massFlux");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_massFlux(cell);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  if(m_momentumCoupling) {
    for(MInt i = 0; i < nDim; i++) {
      name.clear();
      string tmps = "momentumFlux" + to_string(i);
      name.push_back(tmps);
      for(MInt cell = 0; cell < noCells; cell++) {
        tmpW[cell] = a_momentumFlux(cell, i);
      }
      this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
    }
    name.clear();
    name.push_back("workFlux");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_workFlux(cell);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  if(m_heatCoupling) {
    name.clear();
    name.push_back("heatFlux");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_heatFlux(cell);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  name.clear();
  name.push_back("volumeFraction");
  for(MInt cell = 0; cell < noCells; cell++) {
    tmpW[cell] = a_volumeFraction(cell);
  }
  this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

  for(MInt i = 0; i < nDim; i++) {
    name.clear();
    string tmps = "fluidVelocity_" + to_string(i);
    name.push_back(tmps);
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_fluidVelocity(cell, i);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  name.clear();
  name.push_back("fluidDensity");
  for(MInt cell = 0; cell < noCells; cell++) {
    tmpW[cell] = a_fluidDensity(cell);
  }
  this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

  name.clear();
  name.push_back("fluidPressure");
  for(MInt cell = 0; cell < noCells; cell++) {
    tmpW[cell] = a_fluidPressure(cell);
  }
  this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

  name.clear();
  name.push_back("fluidTemperture");
  for(MInt cell = 0; cell < noCells; cell++) {
    tmpW[cell] = a_fluidTemperature(cell);
  }
  this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);

  if(m_evaporation) {
    name.clear();
    name.push_back("fluidSpecies");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_fluidSpecies(cell);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  if(m_ellipsoids) {
    for(MInt i = 0; i < nDim; i++) {
      for(MInt j = 0; j < nDim; j++) {
        name.clear();
        string tmps = "fluidVelocitySlopes_" + to_string(i) + "_" + to_string(j);
        name.push_back(tmps);
        for(MInt cell = 0; cell < noCells; cell++) {
          tmpW[cell] = a_velocitySlope(cell, i, j);
        }
        this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
      }
    }
  }

  if(debugOutput && this->m_adaptation && this->m_noSensors > 0) {
    name.clear();
    name.push_back("regridTrigger");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_regridTrigger(cell);
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  if(debugOutput && m_wallCollisions) {
    name.clear();
    name.push_back("boundaryCellId");
    for(MInt cell = 0; cell < noCells; cell++) {
      tmpW[cell] = a_bndryCellId(cell);
      if(!a_isValidCell(cell)) {
        tmpW[cell] = -2;
      }
    }
    this->collectVariables(tmpW.begin(), dbVariables, name, dbVariablesName, 1, noCells);
  }

  // build file name
  stringstream fileName;
  fileName.clear();
  fileName.str("");
  fileName << outputDir() << "solutionLPT_" << getIdentifier(true) << globalTimeStep << ParallelIo::fileExt();

  this->saveGridFlowVars((fileName.str()).c_str(), gridFileName.c_str(), noCells, noInternalCellIds, dbVariables,
                         dbVariablesName, 0, idVariables, idVariablesName, 0, dbParameters, dbParametersName,
                         idParameters, idParametersName, pointerRecalcIds, -1);
}

/// \brief Add new particle
///
/// \author Sven Berger
/// \date   June 2015
/// \param[in] cellId CellId of the cell in which the particle is to be created.
/// \param[in] diameter Diameter of the added particle.
/// \param[in] particleDensityRatio The density ratio of the added particle.
/// \param[in] random Add random displacement to particle position.
///  \param[in] addMode 0: as initial-condition
///                     1: as particle spawn
///                     2: as spray primaryBreakUp
///                     3: as spray secondaryBreakUp
template <MInt nDim>
MInt LPT<nDim>::addParticle(const MInt cellId, const MFloat diameter, const MFloat particleDensityRatio,
                            const MInt random, const MInt addMode, const MFloat* velocity, const MFloat* position,
                            const MInt parcelSize) {
  LPTSpherical<nDim> thisParticle;
  const MFloat epsilon = MFloatEps;
  const MFloat cellLength = 0.5 * c_cellLengthAtCell(cellId) - epsilon;

  thisParticle.m_cellId = cellId;
  thisParticle.m_oldCellId = cellId;
  thisParticle.m_diameter = diameter;
  thisParticle.m_shedDiam = diameter;
  thisParticle.m_temperature = m_material->T();
  thisParticle.m_noParticles = parcelSize;
  thisParticle.m_dM = 0;

  if(addMode == 0) {
    for(MInt j = 0; j < nDim; j++) {
      // add particle around the cell-coordinate
      // include a small offset w.r.t.(with respect to) the cell coordinates, so that the
      // interpolation stencil is unambiguously defined later on
      thisParticle.m_position[j] = thisParticle.m_oldPos[j] =
          c_coordinate(cellId, j) + epsilon - random * (-cellLength + 2 * cellLength * rand() / (RAND_MAX + 1.0));
      thisParticle.m_velocity[j] = thisParticle.m_oldVel[j] = a_fluidVelocity(cellId, j);

      thisParticle.m_accel[j] = 0.0;
    }
  } else if(addMode == 1 || addMode == 2) {
    // for spray-injection the position and velocity of new particles is specified
    for(MInt j = 0; j < nDim; j++) {
      thisParticle.m_position[j] = position[j];
      thisParticle.m_oldPos[j] = c_coordinate(cellId, j);
      thisParticle.m_velocity[j] = velocity[j];
      thisParticle.m_oldVel[j] = velocity[j];
      thisParticle.m_accel[j] = 0.0;
    }
  } else if(addMode == 3) {
    for(MInt j = 0; j < nDim; j++) {
      thisParticle.m_position[j] = position[j] + MFloatEps;
      thisParticle.m_oldPos[j] = c_coordinate(cellId, j);
      thisParticle.m_velocity[j] = velocity[j];
      thisParticle.m_oldVel[j] = velocity[j];
      thisParticle.m_accel[j] = 0.0;
    }
  } else {
    mTerm(1, AT_, "Unknown particle add-Mode!");
  }

  thisParticle.m_oldFluidDensity = a_fluidDensity(thisParticle.m_cellId);
  for(MInt i = 0; i < nDim; i++) {
    thisParticle.m_oldFluidVel[i] = a_fluidVelocity(thisParticle.m_cellId, i);
  }
  thisParticle.m_densityRatio = particleDensityRatio;

  thisParticle.updateProperties();
  thisParticle.firstStep() = true;

  // update the cellId and compute createTime!
  MFloat time = 0;
  if(addMode == 1) {
    thisParticle.checkCellChange(m_spawnCoord);
  } else if(addMode == 2) {
    ASSERT(m_activePrimaryBUp, "");
    thisParticle.checkCellChange(m_spawnCoord, !m_sprayModel->m_broadcastInjected);
    uniform_real_distribution<MFloat> randMov(0, m_timeStep); // - MFloatEps);
    time = randMov(m_sprayModel->randomPrimBreakUp(1));
    if(time > m_timeStep) {
      time = m_timeStep;
    }
  } else if(addMode == 3) {
    ASSERT(m_activeSecondaryBUp, "");

    thisParticle.m_fluidDensity = a_fluidDensity(thisParticle.m_cellId);
    for(MInt i = 0; i < nDim; i++) {
      thisParticle.m_fluidVel[i] = a_fluidVelocity(thisParticle.m_cellId, i);
    }
  }


  thisParticle.m_creationTime = time;
  thisParticle.m_breakUpTime = -time;

  if(thisParticle.m_creationTime < 0 || thisParticle.m_creationTime > m_timeStep) {
    mTerm(1, AT_, "Invalid creation Time");
  }

#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // create a unique particle id
    thisParticle.m_partId = static_cast<MLong>(domainId() * 1E9 + m_addedParticle);
    m_addedParticle++;
    m_partList.push_back(thisParticle);
  }

  return a_noParticles() - 1;
}

/**
 * \brief Add new ellipsoidal particle
 * \author Laurent Andre
 * \date September 2022
 *
 * \param cellId
 * \param semiMinorAxis
 * \param aspectRatio
 * \param particleDensityRatio
 * \param random
 *
 * \returns
 */
template <MInt nDim>
MInt LPT<nDim>::addEllipsoid(const MInt cellId, const MFloat semiMinorAxis, const MFloat aspectRatio,
                             const MFloat particleDensityRatio, const MInt random, const MInt addMode,
                             const MFloat* velocity, const MFloat* position) {
  LPTEllipsoidal<nDim> thisParticle;
  const MFloat epsilon = MFloatEps;
  const MFloat cellLength = 0.5 * c_cellLengthAtCell(cellId) - epsilon;

  thisParticle.m_cellId = cellId;
  thisParticle.m_oldCellId = cellId;
  thisParticle.m_temperature = m_material->T();
  thisParticle.m_semiMinorAxis = semiMinorAxis;
  thisParticle.m_aspectRatio = aspectRatio;
  thisParticle.initEllipsoialProperties();

  if(addMode == 0) {
    for(MInt j = 0; j < nDim; j++) {
      // add particle around the cell-coordinate
      // include a small offset w.r.t.(with respect to) the cell coordinates, so that the
      // interpolation stencil is unambiguously defined later on
      thisParticle.m_position[j] = thisParticle.m_oldPos[j] =
          c_coordinate(cellId, j) + epsilon - random * (-cellLength + 2 * cellLength * rand() / (RAND_MAX + 1.0));
      thisParticle.m_velocity[j] = thisParticle.m_oldVel[j] = a_fluidVelocity(cellId, j);
      thisParticle.m_angularVel[j] = F0;
      thisParticle.m_accel[j] = F0;
      thisParticle.m_angularAccel[j] = F0;
    }
  } else if(addMode == 1 || addMode == 2) {
    // for spray-injection the position and velocity of new particles is specified
    for(MInt j = 0; j < nDim; j++) {
      thisParticle.m_position[j] = position[j];
      thisParticle.m_oldPos[j] = c_coordinate(cellId, j);
      thisParticle.m_velocity[j] = velocity[j];
      thisParticle.m_oldVel[j] = velocity[j];
      thisParticle.m_angularVel[j] = F0;
      thisParticle.m_accel[j] = F0;
      thisParticle.m_angularAccel[j] = F0;
    }
  } else if(addMode == 3) {
    for(MInt j = 0; j < nDim; j++) {
      thisParticle.m_position[j] = position[j] + MFloatEps;
      thisParticle.m_oldPos[j] = c_coordinate(cellId, j);
      thisParticle.m_velocity[j] = velocity[j];
      thisParticle.m_oldVel[j] = velocity[j];
      thisParticle.m_angularVel[j] = F0;
      thisParticle.m_accel[j] = F0;
      thisParticle.m_angularAccel[j] = F0;
    }
  } else {
    mTerm(1, AT_, "Unknown particle add-Mode!");
  }

  if(m_ellipsoidRandomOrientation != 0) {
    mt19937_64 randomOrientationGen(m_ellipsoidRandomOrientationSeed);
    uniform_real_distribution<MFloat> dist(F0, F1);
    MFloat alpha = dist(randomOrientationGen);
    MFloat beta = dist(randomOrientationGen);
    MFloat gamma = dist(randomOrientationGen);
    thisParticle.m_quaternion[0] = sqrt(1 - alpha) * sin(2 * PI * beta);
    thisParticle.m_quaternion[1] = sqrt(1 - alpha) * cos(2 * PI * beta);
    thisParticle.m_quaternion[2] = sqrt(alpha) * sin(2 * PI * gamma);
    thisParticle.m_quaternion[3] = sqrt(alpha) * cos(2 * PI * gamma);
    MFloat normFactor = F0;
    for(MFloat k : thisParticle.m_quaternion) {
      normFactor += k * k;
    }
    normFactor = sqrt(normFactor);
    for(MInt i = 0; i < 4; i++) {
      thisParticle.m_quaternion[i] /= normFactor;
      thisParticle.m_oldQuaternion[i] = thisParticle.m_quaternion[i];
    }
  }

  thisParticle.m_oldFluidDensity = a_fluidDensity(thisParticle.m_cellId);
  for(MInt i = 0; i < nDim; i++) {
    thisParticle.m_oldFluidVel[i] = a_fluidVelocity(thisParticle.m_cellId, i);
  }
  thisParticle.m_densityRatio = particleDensityRatio;

  thisParticle.updateProperties();
  thisParticle.firstStep() = true;

  // update the cellId and compute createTime!
  MFloat time = 0;
  if(addMode == 1) {
    thisParticle.checkCellChange(m_spawnCoord);
  } else if(addMode == 2) {
    ASSERT(m_activePrimaryBUp, "");
    thisParticle.checkCellChange(m_spawnCoord, !m_sprayModel->m_broadcastInjected);
    uniform_real_distribution<MFloat> randMov(0, m_timeStep); // - MFloatEps);
    time = randMov(m_sprayModel->randomPrimBreakUp(1));
    if(time > m_timeStep) {
      time = m_timeStep;
    }
  } else if(addMode == 3) {
    ASSERT(m_activeSecondaryBUp, "");

    thisParticle.m_fluidDensity = a_fluidDensity(thisParticle.m_cellId);
    for(MInt i = 0; i < nDim; i++) {
      thisParticle.m_fluidVel[i] = a_fluidVelocity(thisParticle.m_cellId, i);
    }
  }

  thisParticle.m_creationTime = time;

  if(thisParticle.m_creationTime < 0 || thisParticle.m_creationTime > m_timeStep) {
    mTerm(1, AT_, "Invalid creation Time");
  }

#ifdef _OPENMP
#pragma omp critical
#endif
  {
    // create a unique particle id
    thisParticle.m_partId = static_cast<MLong>(domainId() * 1E9 + m_addedParticle);
    m_addedParticle++;
    m_partListEllipsoid.push_back(thisParticle);
  }

  return a_noParticles() + a_noEllipsoidalParticles() - 1;
}


/// \brief Spawn new particles (used for testing)
///
/// \author Sven Berger
/// \date   June 2015
template <MInt nDim>
void LPT<nDim>::spawnParticles() {
  // particles are spawned on another domain
  if(m_spawnCellId < 0) {
    if(noDomains() > 1) {
      // receive injected particles instead
      recvInjected();
    }
    m_particleResiduum = 0.0;
    return;
  }
  ASSERT(domainId() == m_spawnDomainId, "");

  const MInt prevNo = m_partList.empty() ? 0 : m_partList.size() - 1;
  const auto noNewParticles = (MInt)(m_spawnParticlesCount * m_timeStep + m_particleResiduum);
  m_particleResiduum += m_spawnParticlesCount * m_timeStep - noNewParticles;

  if(noNewParticles >= 1) {
    MFloatScratchSpace spawnParticlesInitVelo(nDim, AT_, "spawnParticlesInitVelo");

    for(MInt i = 0; i < noNewParticles; i++) {
      m_PRNGSpawnCount +=
          randomVectorInCone(&spawnParticlesInitVelo[0], &m_spawnDir[0], m_spawnVelocity, m_spawnParticlesConeAngle,
                             m_spawnEmittDist, randomSpawn(0), m_spawnDistSigmaCoeff, 0);

      addParticle(m_spawnCellId, m_spawnDiameter, m_material->densityRatio(), 0, 1, &spawnParticlesInitVelo[0],
                  &m_spawnCoord[0], 1);
    }


    broadcastInjected(prevNo);

#ifdef LPT_DEBUG
    cerr << "added " << noNewParticles << " new particles with diameter " << m_spawnDiameter << " to cell "
         << m_spawnCellId << endl;
#endif
  }
}


/// \brief Print information at the start of simulation
///
/// \author Sven Berger
/// \date   January 2018
template <MInt nDim>
void LPT<nDim>::initSummary() {
  cerr << "//////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cerr << "/// Particle INITIALISATION INFORMATION                                                        ///" << endl;
  cerr << "///--------------------------------------------------------------------------------------------///" << endl;
  cerr << "spawning is active                  : " << boolalpha << m_spawnParticles << endl;
  cerr << "spray model is active               : " << boolalpha << (MBool)(m_activePrimaryBUp || m_activeSecondaryBUp)
       << endl;
  cerr << "momentum coupling is active         : " << boolalpha << m_momentumCoupling << endl;
  if(m_momentumCoupling || m_heatCoupling) {
    cerr << "  * redistribution is active        : " << boolalpha << m_couplingRedist << endl;
    cerr << "    - redistribution area is set to : " << m_noRedistLayer << endl;
  }
  cerr << "drag modelling is active            : " << boolalpha << (m_dragModelType > 0) << endl;
  if(m_dragModelType > 0) {
    cerr << "  * drag model used is              : " << m_dragModelType << endl;
  }
  if(m_ellipsoids) {
    cerr << "lift modelling is active            : " << boolalpha << (m_liftModelType > 0) << endl;
    if(m_liftModelType > 0) {
      cerr << "  * lift model used is              : " << m_liftModelType << endl;
    }
    cerr << "torque modelling is active          : " << boolalpha << (m_torqueModelType > 0) << endl;
    if(m_torqueModelType > 0) {
      cerr << "  * torque model used is            : " << m_torqueModelType << endl;
    }
  }
  cerr << "//////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
}

/// \brief receive particles generated during injection
//         if the injection spreads across non-neighbor domains on the min-Level!
/// \author Sven Berger
/// \date   September 2018
template <MInt nDim>
void LPT<nDim>::recvInjected() {
  TRACE();

  MInt partReceive = 0;

  // wait for number
  MPI_Bcast(&partReceive, 1, MPI_INT, m_spawnDomainId, mpiComm(), AT_, "part");

#ifndef NDEBUG
  MInt noParticlesB = -static_cast<MInt>(m_partList.size());
  // Note: cast to signed type first to prevent overflow
#endif

  if(partReceive > 0) {
    // receive
    MPI_Bcast(&m_intRecvBuffer[0][0], partReceive * intElemPerP<LPTSpherical<nDim>>(), MPI_INT, m_spawnDomainId,
              mpiComm(), AT_, "buffer");
    MPI_Bcast(&m_recvBuffer[0][0], partReceive * elemPerP<LPTSpherical<nDim>>(), MPI_DOUBLE, m_spawnDomainId, mpiComm(),
              AT_, "buffer");

    // unpack and check whether they are on this domain
    unpackParticles<LPTSpherical<nDim>, true>(partReceive, m_intRecvBuffer[0].get(), m_recvBuffer[0].get());
  }

#ifndef NDEBUG
  noParticlesB += m_partList.size();
  MPI_Allreduce(MPI_IN_PLACE, &noParticlesB, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noParticlesB");
#endif
}

/// \brief Pack particles for sending
///
/// \author Sven Berger
/// \date   September 2018
/// \param particlesToSend particles that need to be send
/// \param intBuffer integer buffer to use
/// \param floatBuffer float buffer to use
template <MInt nDim>
void LPT<nDim>::packParticles(vector<LPTSpherical<nDim>*>& particlesToSend, MInt* intBuffer, MFloat* floatBuffer,
                              vector<MInt> cellIds) {
  MInt z = 0;
  if(intBuffer != nullptr) {
    for(MInt id = 0; id < (MInt)particlesToSend.size(); id++) {
      // splitting long into two ints
      intBuffer[z++] = (MInt)(particlesToSend[id]->m_partId >> 32);
      intBuffer[z++] = MInt(particlesToSend[id]->m_partId);
      intBuffer[z++] = cellIds[id];
      intBuffer[z++] = particlesToSend[id]->m_noParticles;
      intBuffer[z++] = ((particlesToSend[id]->firstStep() ? 1 : 0) + (particlesToSend[id]->reqSend() ? 2 : 0));
      intBuffer[z++] = (MInt)(particlesToSend[id]->hadWallColl());
      if(id == 0 && z != intElemPerP<LPTSpherical<nDim>>()) {
        mTerm(1, AT_, "Buffersize missmatching!");
      }
    }
  }

  if(floatBuffer != nullptr) {
    z = 0;
    for(MInt id = 0; id < (MInt)particlesToSend.size(); id++) {
      floatBuffer[z++] = particlesToSend[id]->m_diameter;
      floatBuffer[z++] = particlesToSend[id]->m_densityRatio;
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_position.at(j);
      }
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_velocity.at(j);
      }
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_accel.at(j);
      }
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_oldPos.at(j);
      }
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_oldVel.at(j);
      }
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_oldAccel.at(j);
      }
      floatBuffer[z++] = particlesToSend[id]->m_creationTime;
      floatBuffer[z++] = particlesToSend[id]->m_breakUpTime;
      floatBuffer[z++] = particlesToSend[id]->m_shedDiam;
      floatBuffer[z++] = particlesToSend[id]->m_temperature;
      floatBuffer[z++] = particlesToSend[id]->m_dM;
      floatBuffer[z++] = particlesToSend[id]->m_heatFlux;
      floatBuffer[z++] = particlesToSend[id]->m_fluidVelMag;

      floatBuffer[z++] = particlesToSend[id]->m_oldFluidDensity;
      for(MInt j = 0; j < nDim; j++) {
        floatBuffer[z++] = particlesToSend[id]->m_oldFluidVel[j];
      }
      if(id == 0 && z != elemPerP<LPTSpherical<nDim>>()) {
        mTerm(1, AT_, "Buffersize missmatching!");
      }
      if(intBuffer == nullptr) {
        MInt partId1 = (MInt)(particlesToSend[id]->m_partId >> 32);
        MInt partId2 = MInt(particlesToSend[id]->m_partId);
        floatBuffer[z++] = (MFloat)partId1;
        floatBuffer[z++] = (MFloat)partId2;
        floatBuffer[z++] = particlesToSend[id]->m_noParticles;
      }
    }
  }

  if(intBuffer != nullptr && floatBuffer != nullptr) {
    m_noSendParticles++;
  }
}


/**
 * \brief Pack ellipsoids to be sent
 * \author Laurent Andre
 * \date September 2022
 */
template <MInt nDim>
void LPT<nDim>::packParticles(vector<LPTEllipsoidal<nDim>*>& particlesToSend, MInt* intBuffer, MFloat* floatBuffer,
                              vector<MInt> cellIds) {
  if(!m_ellipsoids) {
    mTerm(1, AT_, "Only allowed for ellipsoids!");
  }

  MInt z = 0;
  if(intBuffer != nullptr) {
    for(MInt id = 0; id < (MInt)particlesToSend.size(); id++) {
      // splitting long into two ints
      intBuffer[z++] = (MInt)(particlesToSend[id]->m_partId >> 32);
      intBuffer[z++] = MInt(particlesToSend[id]->m_partId);
      intBuffer[z++] = cellIds[id];
      intBuffer[z++] = ((particlesToSend[id]->firstStep() ? 1 : 0) + (particlesToSend[id]->reqSend() ? 2 : 0));
      intBuffer[z++] = (MInt)(particlesToSend[id]->hadWallColl());
      if(id == 0 && z != intElemPerP<LPTEllipsoidal<nDim>>()) {
        mTerm(1, AT_, "Buffersize missmatching!");
      }
    }
  }

  if(floatBuffer != nullptr) {
    z = 0;
    for(MInt id = 0; id < (MInt)particlesToSend.size(); id++) {
      floatBuffer[z++] = particlesToSend[id]->m_semiMinorAxis;
      floatBuffer[z++] = particlesToSend[id]->m_aspectRatio;
      floatBuffer[z++] = particlesToSend[id]->m_densityRatio;
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_position.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_velocity.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_accel.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_angularVel.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_angularAccel.at(j);
      for(MInt j = 0; j < 4; j++)
        floatBuffer[z++] = particlesToSend[id]->m_quaternion.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldPos.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldVel.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldAccel.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldAngularVel.at(j);
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldAngularAccel.at(j);
      for(MInt j = 0; j < 4; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldQuaternion.at(j);
      floatBuffer[z++] = particlesToSend[id]->m_creationTime;
      floatBuffer[z++] = particlesToSend[id]->m_temperature;
      floatBuffer[z++] = particlesToSend[id]->m_heatFlux;
      floatBuffer[z++] = particlesToSend[id]->m_fluidVelMag;
      floatBuffer[z++] = particlesToSend[id]->m_oldFluidDensity;
      for(MInt j = 0; j < nDim; j++)
        floatBuffer[z++] = particlesToSend[id]->m_oldFluidVel[j];

      if(id == 0 && z != elemPerP<LPTEllipsoidal<nDim>>()) {
        mTerm(1, AT_, "Buffersize missmatching!");
      }
      if(intBuffer == nullptr) {
        MInt partId1 = (MInt)(particlesToSend[id]->m_partId >> 32);
        MInt partId2 = MInt(particlesToSend[id]->m_partId);
        floatBuffer[z++] = (MFloat)partId1;
        floatBuffer[z++] = (MFloat)partId2;
      }
    }
  }
}


/// \brief Unpack particles
///
/// \author Sven Berger
/// \date   September 2018
/// \tparam t_search search for valid cellId, otherwise use known cellId
template <MInt nDim>
template <class LPTParticle, MBool t_search>
void LPT<nDim>::unpackParticles(const MInt num, const MInt* intBuffer, const MFloat* floatBuffer, const MInt domain,
                                const MBool allowNonLeaf) {
  MInt z = 0;
  MInt h = 0;
  MInt i = -1;
  for(MInt id = 0; id < num; id++) {
    LPTParticle thisParticle;
    unpackParticle(thisParticle, intBuffer, z, floatBuffer, h, i, id);

    thisParticle.initProperties();
    thisParticle.firstStep() = ((i % 2) > 0);

    // search for new cellId
    if(t_search) {
      const MInt cellId = grid().findContainingLeafCell(&thisParticle.m_position[0]);

      if(cellId >= 0 && !a_isHalo(cellId)) {
        thisParticle.m_cellId = cellId;
      } else {
        continue;
      }
    } else { // use cellId from halo/window cell collector
      if(i < 2) {
        // before the particle-particle collision step,
        // when window particles are communicated back to the halo part for collision
        thisParticle.m_cellId = haloCellId(domain, thisParticle.m_cellId);
        thisParticle.wasSend() = true;
        thisParticle.toBeDeleted() = false;
        thisParticle.isInvalid() = true;
        ASSERT(m_collisions > 0, "");
      } else {
        // NOTE: at this point the cellId represents the entry to the matching
        //      halo/windowcellId in the corresponding containers!
        if(thisParticle.m_cellId < 0) {
          mTerm(1, AT_, "ERROR: Particle has no valid cell!");
        }
        // regular exchange from halo to window => thus already knowing the matching cellId
        // and avoiding the search for the matching cellId!
        thisParticle.m_cellId = windowCellId(domain, thisParticle.m_cellId);
        // update window/halo cell properties only
        thisParticle.updateProperties(false);

        // loop down and find leaf-windowcell for particles
        // which where on a non-leaf halo-cell
        if(allowNonLeaf && !c_isLeafCell(thisParticle.m_cellId)) {
          thisParticle.m_cellId = grid().findContainingLeafCell(&thisParticle.m_position[0], thisParticle.m_cellId);
          thisParticle.m_oldCellId = thisParticle.m_cellId;

          if(!thisParticle.firstStep()) {
            mTerm(1, AT_, "Non-leaf comm. of not injected cell requested!");
          }
        }
        if(!c_isLeafCell(thisParticle.m_cellId)) {
          mTerm(1, AT_, "ERROR: Particle has no valid leaf-cell!");
        }
      }
    }

    if(thisParticle.m_cellId < 0) {
      mTerm(1, AT_, "ERROR: Particle has no valid cell!");
    }

    if(thisParticle.m_oldCellId < 0) {
      thisParticle.m_oldCellId = grid().findContainingLeafCell(&thisParticle.m_oldPos[0], thisParticle.m_cellId);
    }

    std::vector<LPTParticle>& particleList = a_particleList<LPTParticle>();
    particleList.push_back(thisParticle);
  }
}


/**
 * \brief Unpack spherical particle from buffers
 */
template <MInt nDim>
void LPT<nDim>::unpackParticle(LPTSpherical<nDim>& thisParticle, const MInt* intBuffer, MInt& z,
                               const MFloat* floatBuffer, MInt& h, MInt& step, MInt id) {
  if(intBuffer != nullptr) {
    MInt a = intBuffer[z++];
    MInt b = intBuffer[z++];
    MLong partId = ((MLong)a << 32 | ((MLong)b & 0xFFFFFFFFL));
    thisParticle.m_partId = partId;
    thisParticle.m_cellId = intBuffer[z++];
    thisParticle.m_noParticles = intBuffer[z++];
    step = intBuffer[z++];
    thisParticle.hadWallColl() = (MBool)intBuffer[z++];
    if(id == 0 && z != intElemPerP<LPTSpherical<nDim>>()) {
      mTerm(1, AT_, "Buffersize missmatching!");
    }
  }

  thisParticle.m_diameter = floatBuffer[h++];
  thisParticle.m_densityRatio = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_position.at(l) = floatBuffer[h++];
  }
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_velocity.at(l) = floatBuffer[h++];
  }
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_accel.at(l) = floatBuffer[h++];
  }
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_oldPos.at(l) = floatBuffer[h++];
  }
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_oldVel.at(l) = floatBuffer[h++];
  }
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_oldAccel.at(l) = floatBuffer[h++];
  }
  thisParticle.m_creationTime = floatBuffer[h++];
  thisParticle.m_breakUpTime = floatBuffer[h++];
  thisParticle.m_shedDiam = floatBuffer[h++];
  thisParticle.m_temperature = floatBuffer[h++];
  thisParticle.m_dM = floatBuffer[h++];
  thisParticle.m_heatFlux = floatBuffer[h++];
  thisParticle.m_fluidVelMag = floatBuffer[h++];

  thisParticle.m_oldFluidDensity = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++) {
    thisParticle.m_oldFluidVel[l] = floatBuffer[h++];
  }
  if(id == 0 && h != elemPerP<LPTSpherical<nDim>>()) {
    mTerm(1, AT_, "Buffersize missmatching!");
  }

  if(intBuffer == nullptr) {
    MInt a = (MInt)floatBuffer[h++];
    MInt b = (MInt)floatBuffer[h++];
    MLong partId = ((MLong)a << 32 | ((MLong)b & 0xFFFFFFFFL));
    thisParticle.m_partId = partId;
    thisParticle.m_noParticles = (MInt)floatBuffer[h++];
  }
}


/**
 * \brief Unpack ellipsoidal particle and fill buffers
 * \author Laurent Andre
 * \date Februar 2022
 */
template <MInt nDim>
void LPT<nDim>::unpackParticle(LPTEllipsoidal<nDim>& thisParticle, const MInt* intBuffer, MInt& z,
                               const MFloat* floatBuffer, MInt& h, MInt& step, MInt id) {
  if(intBuffer != nullptr) {
    MInt a = intBuffer[z++];
    MInt b = intBuffer[z++];
    MLong partId = ((MLong)a << 32 | ((MLong)b & 0xFFFFFFFFL));
    thisParticle.m_partId = partId;
    thisParticle.m_cellId = intBuffer[z++];
    step = intBuffer[z++];
    thisParticle.hadWallColl() = (MBool)intBuffer[z++];
    if(id == 0 && z != intElemPerP<LPTEllipsoidal<nDim>>()) {
      mTerm(1, AT_, "Buffersize missmatching!");
    }
  }
  thisParticle.m_semiMinorAxis = floatBuffer[h++];
  thisParticle.m_aspectRatio = floatBuffer[h++];
  thisParticle.initEllipsoialProperties();
  thisParticle.m_densityRatio = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_position.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_velocity.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_accel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_angularVel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_angularAccel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < 4; l++)
    thisParticle.m_quaternion.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldPos.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldVel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldAccel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldAngularVel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldAngularAccel.at(l) = floatBuffer[h++];
  for(MInt l = 0; l < 4; l++)
    thisParticle.m_oldQuaternion.at(l) = floatBuffer[h++];
  thisParticle.m_creationTime = floatBuffer[h++];
  thisParticle.m_temperature = floatBuffer[h++];
  thisParticle.m_heatFlux = floatBuffer[h++];
  thisParticle.m_fluidVelMag = floatBuffer[h++];
  thisParticle.m_oldFluidDensity = floatBuffer[h++];
  for(MInt l = 0; l < nDim; l++)
    thisParticle.m_oldFluidVel[l] = floatBuffer[h++];

  if(id == 0 && h != elemPerP<LPTEllipsoidal<nDim>>()) {
    mTerm(1, AT_, "Buffersize missmatching!");
  }

  if(intBuffer == nullptr) {
    MInt a = (MInt)floatBuffer[h++];
    MInt b = (MInt)floatBuffer[h++];
    MLong partId = ((MLong)a << 32 | ((MLong)b & 0xFFFFFFFFL));
    thisParticle.m_partId = partId;
  }
}


/// \brief Send newly created particles
///
/// \author Sven Berger
/// \date   October 2018
/// \param[in] prevNumPart Previous number of particles.
template <MInt nDim>
void LPT<nDim>::broadcastInjected(const MUint prevNumPart) {
  TRACE();

  if(noDomains() > 1) {
    vector<LPTSpherical<nDim>*> particlesToSend;
    vector<MInt> cellIds;
    vector<LPTSpherical<nDim>>& partList = a_particleList<LPTSpherical<nDim>>();
    // collect particles which are marked as inactive
    for(MInt i = prevNumPart; i < a_noParticles(); i++) {
      if(partList[i].toBeDeleted() || partList[i].toBeRespawn()) {
        particlesToSend.push_back(&partList[i]);
        cellIds.push_back(-1);
#ifdef LPT_DEBUG
        cerr << domainId() << ": " << partList[i].m_partId << " is to be send to all domains during injection." << endl;
        << endl;
#endif
      }
    }

    MInt numPartSend = static_cast<MInt>(particlesToSend.size());

    // send number of particles
    MPI_Bcast(&numPartSend, 1, MPI_INT, m_spawnDomainId, mpiComm(), AT_, "intSendBuffer");

#ifdef LPT_DEBUG
    cerr << domainId() << ": "
         << " is sending " << numPartSend << " particles" << endl;
#endif

    ASSERT(m_exchangeBufferSize
               > max(numPartSend * intElemPerP<LPTSpherical<nDim>>(), numPartSend * elemPerP<LPTSpherical<nDim>>()),
           "ERROR: Exchange buffer to small");

    if(numPartSend > 0) {
      // pack particles
      packParticles(particlesToSend, m_intSendBuffer[0].get(), m_sendBuffer[0].get(), cellIds);

      // send particles
      MPI_Bcast(&m_intSendBuffer[0][0], numPartSend * intElemPerP<LPTSpherical<nDim>>(), type_traits<MInt>::mpiType(),
                m_spawnDomainId, mpiComm(), AT_, "intSendBuffer");

      MPI_Bcast(&m_sendBuffer[0][0], numPartSend * elemPerP<LPTSpherical<nDim>>(), type_traits<MFloat>::mpiType(),
                m_spawnDomainId, mpiComm(), AT_, "intSendBuffer");
    }

    for(auto part : particlesToSend) {
      part->wasSend() = true;
      part->reqSend() = false;
      part->isInvalid() = true;
    }


#ifndef NDEBUG
    MInt particlesReceived = 0;
    MPI_Allreduce(MPI_IN_PLACE, &particlesReceived, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE",
                  "particlesReceived");
    if(particlesReceived != numPartSend) {
      cerr << "total particles received: " << particlesReceived << endl;
      TERMM(-1, "ERROR: Inconsistent number of particles");
    }
#endif
  }
}

template <MInt nDim>
void LPT<nDim>::updateFluidFraction() {

#ifdef _OPENMP
#pragma omp single
#endif
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_volumeFraction(cellId) = 0.0;
  }

#ifdef _OPENMP
#pragma omp single
#endif
  for(MInt i = 0; i < a_noParticles(); i++) {
    const MInt cellId = m_partList[i].m_cellId;
    // particle Volume Fraction = particleVolume * parcelSize / volume of the cartesian cell!
    a_volumeFraction(cellId) +=
        4.0 / 3.0 * PI * POW3(0.5 * m_partList[i].m_diameter) * m_partList[i].m_noParticles / c_cellVolume(cellId);
  }
  if(m_ellipsoids) {
#ifdef _OPENMP
#pragma omp single
#endif
    for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
      const MInt cellId = m_partListEllipsoid[i].m_cellId;
      // particle Volume Fraction = particleVolume / volume of the cartesian cell!
      a_volumeFraction(cellId) += m_partListEllipsoid[i].particleVolume() / c_cellVolume(cellId);
    }
  }

  // correct cell-volume for bndrycells
  for(MInt bndryCellId = 0; bndryCellId < m_bndryCells->size(); bndryCellId++) {
    const MInt lptCellId = m_bndryCells->a[bndryCellId].m_cellId;
    a_volumeFraction(lptCellId) *= c_cellVolume(lptCellId) / m_bndryCells->a[bndryCellId].m_volume;
  }

#if !defined NDEBUG
  /*
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_volumeFraction(cellId) > 0.9) {
      cerr << "Particle volume in cell is almost larger than the cartesian cell! " << a_volumeFraction(cellId)
           << " This is not possible, particles should have left the cell before! " << endl;
    }
  }
  */
#endif
}

/*! \brief Set exchange properties in LPT cell collector (a_isHalo/a_isWindow)
 *         and set m_pointsToHalo and m_pointsToWindow lists
 * \author Tim Wegmann
 * \date November 2020
 */
template <MInt nDim>
void LPT<nDim>::updateExchangeCells() {
  TRACE();

  const MInt noNghbrDomains = grid().noNeighborDomains();

  // set lpt cell collector property
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_isHalo(cellId) = false;
    a_isWindow(cellId) = false;
  }

  m_pointsToHalo.clear();
  m_pointsToHalo.resize(noNghbrDomains);

  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noWindowCells(i); j++) {
      const MInt windowCell = windowCellId(i, j);
      a_isWindow(windowCell) = true;
      if(!c_isLeafCell(windowCell)) continue;
      // set pointsToHalo for leaf-Cell window-Cells.
      // this is used in the pushToQueue-function in particleTransfer!
      // Meaning that the particle is only transferred if the corresponding cellId
      // is found in this list!
      const MInt windowCellOffset = j;
      m_pointsToHalo[i].insert(make_pair(windowCell, windowCellOffset));
    }
  }

  m_pointsToWindow.clear();
  m_pointsToWindow.resize(noNghbrDomains);

  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noHaloCells(i); j++) {
      const MInt haloCell = haloCellId(i, j);
      a_isHalo(haloCell) = true;
      if(!c_isLeafCell(haloCell)) continue;
      const MInt HaloCellOffset = j;
      // set pointsToHalo for leaf-Cell halo-Cells.
      // this is used in the pushToQueue-function in particleTransfer!
      // Meaning that the particle is only transferred if the corresponding cellId
      // is found in this list!
      // NOTE: the search is now only done for matching cells(in this case halo-cells!)
      m_pointsToWindow[i].insert(make_pair(haloCell, HaloCellOffset));
    }
  }
}

/*! \brief add an LPT bndryCell to the collector and fill the data
 * \author Tim Wegmann
 * \date November 2020
 */
template <MInt nDim>
void LPT<nDim>::addBndryCell(const MInt cellId, MFloat* bndryData, MFloat* surfaceN, MFloatScratchSpace& surfaceC,
                             MFloat* surfaceV) {
  TRACE();

  ASSERT(a_bndryCellId(cellId) < 0, "ERROR, cell already has a bndryCell!");

  // currently the LPT-bndryCell consist of the following data-structure:
  //- matching cellId
  //- coordinate of the cut-Cell
  //  (note, that in the FV-solver, the coordinate is only the shift from the un-cut center!)
  //- volume of the cut-Cell
  //- a single bndry-Surface with the following data:
  //  - surface normal
  //  - surface velocity (only unqual zero for moving bndries!)
  //  - surface centroid coordinate
  //  - coordinate of 3 points on the cell-edges of the surface

  const MInt bndryCellId = m_bndryCells->size();
  m_bndryCells->append();
  a_bndryCellId(cellId) = bndryCellId;

  m_bndryCells->a[bndryCellId].m_cellId = cellId;
  m_bndryCells->a[bndryCellId].m_volume = bndryData[0];
  for(MInt n = 0; n < nDim; n++) {
    m_bndryCells->a[bndryCellId].m_coordinates[n] = c_coordinate(cellId, n) + bndryData[1 + n];
  }
  const MInt noSurf = 1;
  m_bndryCells->a[bndryCellId].m_noSrfcs = noSurf;
  for(MInt s = 0; s < noSurf; s++) {
    for(MInt n = 0; n < nDim; n++) {
      m_bndryCells->a[bndryCellId].m_srfcs[s]->m_normal[n] = surfaceN[noSurf * s + n];
      m_bndryCells->a[bndryCellId].m_srfcs[s]->m_velocity[n] = surfaceV[noSurf * s + n];
      m_bndryCells->a[bndryCellId].m_srfcs[s]->m_coordinates[n] = surfaceC(noSurf * s + n, 0);
      for(MInt i = 0; i < nDim; i++) {
        m_bndryCells->a[bndryCellId].m_srfcs[s]->m_planeCoordinates[i][n] = surfaceC(noSurf * s + n, i + 1);
      }
    }
  }
}


template <MInt nDim>
void LPT<nDim>::writeCollData() {
  m_collisionModel->writeCollData();
}

/// \brief some debug check for partcles
/// \author Tim Wegmann
/// \date   Aug 2021
template <MInt nDim>
void LPT<nDim>::checkParticles() {
  TRACE();

  if(globalTimeStep < 0) return;

    // check that particle starts on a valid location
#if defined LPT_DEBUG || !defined NDEBUG

  // calculate weights of the redistribution
  static const MInt noOfRedistLayers = m_couplingRedist ? ((m_noRedistLayer > 0) ? m_noRedistLayer : 2) : 0;

  if(!((noOfRedistLayers == 0 && !m_couplingRedist) || m_couplingRedist)) {
    mTerm(1, AT_, "Invalid configuration noOfRedistLayers > 0, but coupling redistribution is not active!");
  }

  for(MInt part = 0; part < a_noParticles(); part++) {
    if(m_partList[part].isInvalid()) {
      if(!m_partList[part].fullyEvaporated()) continue;
      if(m_partList[part].m_noParticles < 1) {
        cerr << "Fully evap. has neg. parcel count! " << m_partList[part].m_noParticles << " " << m_partList[part].m_dM
             << endl;
      }
      if(m_partList[part].m_dM > 0.000001) {
        cerr << "Fully evap. has large mass-transfer " << m_partList[part].m_dM << " " << m_partList[part].m_partId
             << endl;
      }
      continue;
    }
    const MInt cellId = m_partList[part].m_cellId;
    if(cellId < 0 || cellId > noInternalCells()) {
      cerr << "Invalid cellId " << cellId << " " << m_partList[part].hasCollided() << " " << a_isHalo(cellId) << " "
           << c_isLeafCell(cellId) << " " << m_partList[part].isWindow() << " " << m_partList[part].reqSend() << " "
           << m_partList[part].wasSend() << " " << m_partList[part].hadWallColl() << " " << part << " "
           << m_partList[part].m_partId << endl;
      mTerm(1, AT_, "Invalid cellId!");
    }
    if(m_partList[part].isInvalid()) {
      mTerm(1, AT_, "Invalid particle at TS beginning!");
    }
    if(a_isHalo(cellId) || !c_isLeafCell(cellId) || c_noChildren(cellId) > 0) {
      mTerm(1, AT_, "Particle in halo or non-leaf cell!");
    }
    // check particle coordinate:
    const MFloat halfCellLength = c_cellLengthAtLevel(c_level(cellId) + 1);
    std::array<MFloat, 3> origCellC = {};
    for(MInt i = 0; i < nDim; i++) {
      origCellC[i] = c_coordinate(cellId, i);
    }
    // before particle motion, the particle should be located in the correct cell
    if(fabs(origCellC[0] - m_partList[part].m_position[0]) > halfCellLength
       || fabs(origCellC[1] - m_partList[part].m_position[1]) > halfCellLength
       || fabs(origCellC[2] - m_partList[part].m_position[2]) > halfCellLength) {
      cerr << "ERROR: particle " << m_partList[part].m_partId << " is in on incorrect cell "
           << sqrt(POW2(origCellC[0] - m_partList[part].m_position[0])
                   + POW2(origCellC[1] - m_partList[part].m_position[1])
                   + POW2(origCellC[2] - m_partList[part].m_position[2]))
           << " > " << halfCellLength << " and particle location: " << m_partList[part].m_position[0] << " "
           << m_partList[part].m_position[1] << " " << m_partList[part].m_position[2] << " "
           << m_partList[part].hadWallColl() << " " << a_level(cellId) << " " << origCellC[0] << " " << origCellC[1]
           << " " << origCellC[2] << endl;
      // mTerm(1, AT_, "Invalid cell for particle position!" );
      m_partList[part].isInvalid() = true;
      continue;
    }

    if(m_partList[part].m_diameter < m_sizeLimit || isnan(m_partList[part].m_diameter)) {
      cerr << "Invalid particle diameter : " << m_partList[part].m_diameter << " " << m_sizeLimit << " "
           << m_partList[part].firstStep() << " " << m_partList[part].hadWallColl() << " "
           << m_partList[part].m_breakUpTime << " " << m_partList[part].m_partId << endl;
      // mTerm(1, AT_, "Invalid particle diameter!");
    }

    if(m_partList[part].m_noParticles < 0) {
      cerr << m_partList[part].m_noParticles << endl;
      mTerm(1, AT_, "Invalid parcel count!");
    }

    if(m_partList[part].m_oldCellId < 0 && !m_periodicBC) {
      mTerm(1, AT_, "Invalid old CellId!");
    }

    if(m_partList[part].m_oldFluidDensity < MFloatEps) {
      cerr << a_fluidDensity(cellId) << " " << a_isHalo(cellId) << " " << a_level(cellId) << " " << part << " "
           << a_noParticles() << " " << m_partList[part].m_oldFluidDensity << endl;
      mTerm(1, AT_, "Invalid old fluid density!");
    }

    for(MInt i = 0; i < nDim; i++) {
      if(isnan(m_partList[part].m_oldFluidVel[i])) {
        mTerm(1, AT_, "Invalid old fluid velocity!");
      }
    }

    if(m_partList[part].m_temperature < 0) {
      mTerm(1, AT_, "Invalid particle temperature!");
    }

    if(m_partList[part].m_shedDiam < 0) {
      cerr << "Negative shed-diameter " << m_partList[part].m_partId << endl;
    }
    if(m_partList[part].m_noParticles < 1) {
      cerr << "Negative parcel-count " << m_partList[part].m_partId << endl;
    }
    if(m_partList[part].m_dM > 0.000001) {
      cerr << "Large mass-transfer " << m_partList[part].m_dM << " " << m_partList[part].m_partId << endl;
    }
  }

  // ellipsoids
  for(MInt part = 0; part < a_noEllipsoidalParticles(); part++) {
    const MInt cellId = m_partListEllipsoid[part].m_cellId;

    if(cellId < 0 || cellId > noInternalCells()) {
      cerr << "Invalid cellId " << cellId << " " << m_partListEllipsoid[part].hasCollided() << " " << a_isHalo(cellId)
           << " " << c_isLeafCell(cellId) << " " << m_partListEllipsoid[part].isWindow() << " "
           << m_partListEllipsoid[part].reqSend() << " " << m_partListEllipsoid[part].wasSend() << " "
           << m_partListEllipsoid[part].hadWallColl() << " " << part << " " << m_partListEllipsoid[part].m_partId
           << endl;
      mTerm(1, AT_, "Invalid cellId!");
    }
    if(m_partListEllipsoid[part].isInvalid()) {
      mTerm(1, AT_, "Invalid ellipsoid at TS beginning!");
    }


    if(a_isHalo(cellId) || !c_isLeafCell(cellId) || c_noChildren(cellId) > 0) {
      mTerm(1, AT_, "Ellipsoid in halo or non-leaf cell!");
    }

    // check particle coordinate:
    const MFloat halfCellLength = c_cellLengthAtLevel(c_level(cellId) + 1);

    std::array<MFloat, 3> origCellC = {};
    for(MInt i = 0; i < nDim; i++) {
      origCellC[i] = c_coordinate(cellId, i);
    }

    // before particle motion, the particle should be located in the correct cell
    if(fabs(origCellC[0] - m_partListEllipsoid[part].m_position[0]) > halfCellLength
       || fabs(origCellC[1] - m_partListEllipsoid[part].m_position[1]) > halfCellLength
       || fabs(origCellC[2] - m_partListEllipsoid[part].m_position[2]) > halfCellLength) {
      cerr << "ERROR: ellipsoid " << m_partListEllipsoid[part].m_partId << " is in on incorrect cell "
           << sqrt(POW2(origCellC[0] - m_partListEllipsoid[part].m_position[0])
                   + POW2(origCellC[1] - m_partListEllipsoid[part].m_position[1])
                   + POW2(origCellC[2] - m_partListEllipsoid[part].m_position[2]))
           << " > " << halfCellLength << " and ellipsoid location: " << m_partListEllipsoid[part].m_position[0] << " "
           << m_partListEllipsoid[part].m_position[1] << " " << m_partListEllipsoid[part].m_position[2] << " "
           << m_partListEllipsoid[part].hadWallColl() << " " << a_level(cellId) << " " << origCellC[0] << " "
           << origCellC[1] << " " << origCellC[2] << endl;
      // mTerm(1, AT_, "Invalid cell for particle position!" );
      m_partListEllipsoid[part].isInvalid() = true;
      continue;
    }

    if(m_partListEllipsoid[part].equivalentDiameter() < m_sizeLimit
       || isnan(m_partListEllipsoid[part].equivalentDiameter())) {
      cerr << "Invalid ellipsoid diameter : " << m_partListEllipsoid[part].equivalentDiameter() << " " << m_sizeLimit
           << " " << m_partListEllipsoid[part].firstStep() << " " << m_partListEllipsoid[part].hadWallColl() << " "
           << m_partListEllipsoid[part].m_partId << endl;
      // mTerm(1, AT_, "Invalid particle diameter!");
    }

    if(m_partListEllipsoid[part].m_oldCellId < 0 && !m_periodicBC) {
      mTerm(1, AT_, "Invalid old CellId!");
    }

    if(m_partListEllipsoid[part].m_oldFluidDensity < 0) {
      cerr << a_fluidDensity(cellId) << " " << a_isHalo(cellId) << " " << a_level(cellId) << "" << part << " "
           << a_noEllipsoidalParticles() << endl;
      mTerm(1, AT_, "Invalid old fluid density!");
    }

    for(MInt i = 0; i < nDim; i++) {
      if(isnan(m_partListEllipsoid[part].m_oldFluidVel[i])) {
        mTerm(1, AT_, "Invalid old fluid velocity!");
      }
    }

    if(m_partListEllipsoid[part].m_temperature < 0) {
      mTerm(1, AT_, "Invalid particle temperature!");
    }
  }
#endif
}


/// \brief some debug check for cells
/// \author Tim Wegmann
/// \date   Aug 2021
template <MInt nDim>
void LPT<nDim>::checkCells() {
  TRACE();

// check that particle starts on a valid location
#ifdef LPT_DEBUG
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(!c_isLeafCell(cellId)) continue;
    if(a_fluidDensity(cellId) < MFloatEps) {
      mTerm(1, AT_, "Invalid density in cell!");
    }

    if(isnan(a_massFlux(cellId))) {
      mTerm(1, AT_, "Invalid mass-flux in cell!");
    }
  }

#endif
}

/**
 * \brief actually doing some pre-balance-stuff
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::resetSolver() {
  TRACE();

  if(!isActive()) {
    m_injData.clear();
  } else {
    receiveFlowField();
    receiveVelocitySlopes();
    waitForSendReqs();
  }

  // set correct size of m_injData on inactive ranks!
  MInt vectorSize = 0;
  if(m_activePrimaryBUp && grid().hasInactiveRanks()) {
    if(!m_injData.empty()) {
      vectorSize = m_injData.size();
    }
    MPI_Allreduce(MPI_IN_PLACE, &vectorSize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, AT_, "INPLACE", "vectorSize");

    if(vectorSize > 0 && !isActive()) {
      m_injData.clear();
      const MInt count = m_sprayModel->m_injDataSize + m_addInjDataCnt;
      MInt dummyTimeStep = -1;
      for(MInt i = 0; i < vectorSize; i++) {
        vector<MFloat> tmp(count);
        for(MInt j = 0; j < count; j++) {
          tmp[j] = 0.0;
        }
        m_injData.insert(make_pair(dummyTimeStep, tmp));
        dummyTimeStep--;
      }
    } else if(isActive() && domainId() == m_spawnDomainId) {
      if(vectorSize != (signed)m_injData.size()) {
        cerr << "Injector size differs! " << vectorSize << " " << m_injData.size() << endl;
      }
    }
  }

  if(!isActive()) {
    return;
  }

  cerr0 << "LPT-time step is : " << m_timeStep << endl;

  removeInvalidParticles(true);

  checkParticles(); // check particles before balance!

  countParticlesInCells(); // used to determine the data size necessary for sending

  // remove all bndryCells
  if(m_wallCollisions) {
    m_bndryCells->resetSize(0);
    m_noOuterBndryCells = 0;
  }

  // NOTE:latest point to write a debug file before the balance:
  /*
  MInt backup = globalTimeStep;
  globalTimeStep = -1;
  saveDebugRestartFile();
  globalTimeStep = backup;
  */

  perCellStats(); // cell-particle mapping necessary to pack particles for a given cell

  MPI_Allreduce(MPI_IN_PLACE, &m_particleResiduum, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "INPLACE",
                "m_particleResiduum");

  MPI_Allreduce(MPI_IN_PLACE, &m_PRNGSpawnCount, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "m_PRNGSpawnCount");

  MPI_Allreduce(MPI_IN_PLACE, &m_spawnDomainId, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "m_spawnDomainId");

  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    MInt count = m_sprayModel->m_PRNGPrimBreakUpCount;
    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "m_PRNGPrimBreakUpCount");
    m_sprayModel->m_PRNGPrimBreakUpCount = count;

    MFloat timeSinceSOI = 0;
    if(m_spawnCellId > -1) {
      timeSinceSOI = m_sprayModel->timeSinceSOI();
    }
    MPI_Allreduce(MPI_IN_PLACE, &timeSinceSOI, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "INPLACE", "timeSinceSOI");
    m_sprayModel->timeSinceSOI() = timeSinceSOI;
  }

  // exchange the injection-data
  if(m_activePrimaryBUp) {
    MInt vectorSizeBU = vectorSize;
    vectorSize = m_injData.size();
    const MInt count = m_sprayModel->m_injDataSize + m_addInjDataCnt;

    MPI_Allreduce(MPI_IN_PLACE, &vectorSize, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "vectorSize");
    const MInt totalSize = vectorSize * (count + 1);

    if(vectorSize > 0) {
      MFloatScratchSpace injData(vectorSize * (count + 1), AT_, "injData");
      injData.fill(0.0);
      if(vectorSize != vectorSizeBU && grid().hasInactiveRanks()) {
        cerr << domainId() << "Vector size differs " << vectorSizeBU << " " << vectorSize << endl;
      }

      if(domainId() == m_spawnDomainId) {
        ASSERT((signed)m_injData.size() == vectorSize, "");
        MInt i = 0;
        for(auto it = m_injData.begin(); it != m_injData.end(); it++) {
          const MFloat timeStep = (MFloat)it->first;
          injData(i++) = timeStep;
          for(MInt j = 0; j < count; j++) {
            const MFloat data = (it->second)[j];
            injData(i++) = data;
          }
        }
      } else {
        MInt i = 0;
        for(auto it = m_injData.begin(); it != m_injData.end(); it++) {
          i = i + m_sprayModel->m_injDataSize + 1;
          const MFloat sumEvap = (it->second)[m_sprayModel->m_injDataSize];
          const MInt numRT = (it->second)[m_sprayModel->m_injDataSize + 1];
          const MInt numKH = (it->second)[m_sprayModel->m_injDataSize + 2];
          const MInt numSend = (it->second)[m_sprayModel->m_injDataSize + 3];
          injData(i++) = sumEvap;
          injData(i++) = numRT;
          injData(i++) = numKH;
          injData(i++) = numSend;
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &injData[0], totalSize, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "INPLACE", "injData");

      // reset injection data with summed data
      m_injData.clear();

      for(MInt i = 0; i < totalSize; i = i + (count + 1)) {
        const MInt timeStep = (MInt)injData(i);
        vector<MFloat> tmp(count);
        for(MInt j = 0; j < count; j++) {
          tmp[j] = injData[i + j + 1];
        }
        m_injData.insert(make_pair(timeStep, tmp));
      }
    }
  }

#ifdef LPT_DEBUG
  MInt noParticles = a_noParticles();

  MPI_Allreduce(MPI_IN_PLACE, &noParticles, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noParticles");

  if(domainId() == 0) {
    cerr << noParticles << " #particles before balancing!" << endl;
  }

  if(m_ellipsoids) {
    MInt noEllipsoids = a_noEllipsoidalParticles();
    MPI_Allreduce(MPI_IN_PLACE, &noEllipsoids, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noEllipsoids");
    cerr0 << noEllipsoids << " #ellipsoids before balancing!" << endl;
  }

#endif
}

/**
 * \brief reset the solver during balancing
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::resetSolverFull() {
  m_queueToSend.clear();
  m_sendBuffer.clear();
  m_recvBuffer.clear();
  m_intSendBuffer.clear();
  m_intRecvBuffer.clear();
  m_mpi_reqSendFloat.clear();
  m_mpi_reqSendInt.clear();
  m_mpi_reqSendSize.clear();
  m_mpi_statusProbe.clear();

  m_partList.clear();
  m_partListEllipsoid.clear();
  m_cells.clear();
  if(m_wallCollisions) {
    m_bndryCells->resetSize(0);
  }

  m_cellToNghbrHood.clear();
  m_cellToNghbrHoodInterpolation.clear();

  m_pointsToHalo.clear();
  m_pointsToWindow.clear();

  m_cellToPartMap.clear();
  m_cellToEllipsMap.clear();
}

/**
 * \brief actually doing same pre-partitioneing stuff during balcne
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::cancelMpiRequests() {
  if(isActive()) {
    countParticlesInCells();
  }
}


/**
 * \brief reset the solver during balancing
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::balancePre() {
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

  // Reset cell, boundary cell data and deallocate halo/window cell arrays
  resetSolverFull();

  // check for empry cell and particle collector
  ASSERT(m_cells.size() == 0, "");
  ASSERT(a_noParticles() == 0, "");
  ASSERT(a_noEllipsoidalParticles() == 0, "");

  // Return if solver is not active
  if(!grid().isActive()) {
    m_cells.reset(grid().maxNoCells());
    this->setHaloCellsOnInactiveRanks();
    return;
  }

  initCollector();

  // Check that global ids are sorted
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    if(grid().domainOffset(domainId()) + (MLong)cellId != c_globalId(cellId)) {
      TERMM(1, "Global id mismatch.");
    }
  }

  updateExchangeCells();
  grid().updateLeafCellExchange();
  initMPI();


  // this can only be done, after the halo property has been set above
  this->checkNoHaloLayers();
}

/**
 * \brief Reinitialize solver after setting solution data in DLB.
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::balancePost() {
  TRACE();

  m_loadBalancingReinitStage = 1;

  ASSERT(m_cells.size() == c_noCells() && c_noCells() == a_noCells(), "");
}

/**
 * \brief Return data size to be communicated during DLB for a grid cell and given data id
 * \author Tim Wegmann
 */
template <MInt nDim>
MInt LPT<nDim>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  // Inactive ranks do not have any data to communicate
  if(!isActive()) {
    return 0;
  }

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0 || cellId >= noInternalCells()) {
    return 0;
  }

  // only noParticlesInCell and particle information is balanced
  if(dataId == 0) return 2; // 2 because we spherical and elllipsoidal particles

  if(a_noParticlesInCell(cellId) == 0 && a_noEllipsoidsInCell(cellId) == 0) return 0;

  MInt dataSize = 0;

  switch(dataId) {
    case 1:
      dataSize = (elemPerP<LPTSpherical<nDim>>() + 3) * a_noParticlesInCell(cellId)
                 + (elemPerP<LPTEllipsoidal<nDim>>() + 2) * a_noEllipsoidsInCell(cellId);
      break;
    case 2: // spheres particleId(2) and noParcels, ellipsoids only particleId(2)
      dataSize = 3 * a_noParticlesInCell(cellId) + 2 * a_noEllipsoidsInCell(cellId);
      break;
    default:
      TERMM(1, "Unknown data id.");
      break;
  }

  return dataSize;
}

/**
 * \brief Store the solver data for a given data id ordered in the given buffer for DLB
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                               MFloat* const data) {
  TRACE();

  if(dataId != 1) mTerm(1, AT_, "Only Type 1 should be FLOAT");

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    if(a_noParticlesInCell(cellId) > 0) {
      auto particlesInCell = m_cellToPartMap.equal_range(cellId);
      // loop over all particles in the cell
      vector<LPTSpherical<nDim>*> particlesToSend;
      vector<MInt> cellIds;
      for(auto it = particlesInCell.first; it != particlesInCell.second; it++) {
        auto& particle = it->second;
        particlesToSend.push_back(particle);
      }
      if((MInt)particlesToSend.size() != a_noParticlesInCell(cellId)) {
        cerr << particlesToSend.size() << " " << a_noParticlesInCell(cellId) << " miss-count!" << endl;
        mTerm(1, AT_, "Particle count not matching!");
      }
      packParticles(particlesToSend, nullptr, &data[localBufferId], cellIds);
      localBufferId += ((elemPerP<LPTSpherical<nDim>>() + 3) * particlesToSend.size());
    }

    if(a_noEllipsoidsInCell(cellId) > 0) {
      auto ellipsoidsInCell = m_cellToEllipsMap.equal_range(cellId);
      vector<LPTEllipsoidal<nDim>*> ellipsoidsToSend;
      vector<MInt> cellIds;
      for(auto it = ellipsoidsInCell.first; it != ellipsoidsInCell.second; it++) {
        auto& particle = it->second;
        ellipsoidsToSend.push_back(particle);
      }
      if((MInt)ellipsoidsToSend.size() != a_noEllipsoidsInCell(cellId)) {
        cerr << ellipsoidsToSend.size() << " " << a_noEllipsoidsInCell(cellId) << " miss-count!" << endl;
        mTerm(1, AT_, "Ellipsoids count not matching!");
      }
      packParticles(ellipsoidsToSend, nullptr, &data[localBufferId], cellIds);
      localBufferId += ((elemPerP<LPTEllipsoidal<nDim>>() + 2) * ellipsoidsToSend.size());
    }
  }
}

template <MInt nDim>
void LPT<nDim>::getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                               MInt* const data) {
  TRACE();

  if(dataId != 2 && dataId != 0) mTerm(1, AT_, "Type 0 or 2 should be INT");

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    if(dataId == 0) {
      data[localBufferId++] = a_noParticlesInCell(cellId);
      data[localBufferId++] = a_noEllipsoidsInCell(cellId);
      continue;
    }

    if(a_noParticlesInCell(cellId) > 0) {
      auto particlesInCell = m_cellToPartMap.equal_range(cellId);
      // loop over all particles in the cell
      vector<LPTSpherical<nDim>*> particlesToSend;
      for(auto it = particlesInCell.first; it != particlesInCell.second; it++) {
        auto& particle = it->second;
        data[localBufferId] = (MInt)(particle->m_partId >> 32);
        data[localBufferId + 1] = MInt(particle->m_partId);
        data[localBufferId + 2] = particle->m_noParticles;
        // ASSERT(data[localBufferId + 2] > 0 , "");
        if(data[localBufferId + 2] < 0) {
          cerr << "Setting strange exchange buffer!" << endl;
        }
        localBufferId += 3;
      }
    }

    if(a_noEllipsoidsInCell(cellId) > 0) {
      auto ellipsoidsInCell = m_cellToEllipsMap.equal_range(cellId);
      vector<LPTEllipsoidal<nDim>*> ellipsoidsToSend;
      for(auto it = ellipsoidsInCell.first; it != ellipsoidsInCell.second; it++) {
        auto& particle = it->second;
        data[localBufferId] = (MInt)(particle->m_partId >> 32);
        data[localBufferId + 1] = MInt(particle->m_partId);
        localBufferId += 2;
      }
    }
  }
}

/**
 * \brief Set the solver cell data after DLB
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::setCellDataDlb(const MInt dataId, const MFloat* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  if(dataId != 1) mTerm(1, AT_, "Only Type 1 should be FLOAT");

  // Set the no-particles in cell if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    // gather the number of particles this rank will have after the balance
    MInt noParticlesOnRank = 0;
    MInt noEllipsoidsOnRank = 0;
    for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
      noParticlesOnRank += a_noParticlesInCell(cellId);
      noEllipsoidsOnRank += a_noEllipsoidsInCell(cellId);
    }

    unpackParticles<LPTSpherical<nDim>, true>(noParticlesOnRank, nullptr, &data[0]);
    if(a_noParticles() != noParticlesOnRank) {
      cerr << domainId() << " has incorrect count " << a_noSphericalParticles() << " " << noParticlesOnRank << endl;
      mTerm(1, AT_, "Particle gone missing!");
    }

    if(m_ellipsoids) {
      unpackParticles<LPTEllipsoidal<nDim>, true>(noEllipsoidsOnRank, nullptr,
                                                  &data[(elemPerP<LPTSpherical<nDim>>() + 3) * noParticlesOnRank]);
      if(a_noEllipsoidalParticles() != noEllipsoidsOnRank) {
        cerr << domainId() << " has incorrect count " << a_noParticles<LPTEllipsoidal<nDim>>() << " "
             << noEllipsoidsOnRank << endl;
        mTerm(1, AT_, "Ellipsoid gone missing!");
      }
    }
  }
}

template <MInt nDim>
void LPT<nDim>::setCellDataDlb(const MInt dataId, const MInt* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  if(dataId != 2 && dataId != 0) mTerm(1, AT_, "Type 0 or 2 should be INT");

  if(m_loadBalancingReinitStage == 0) {
    if(dataId == 0) {
      MInt counter = 0;
      for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
        a_noParticlesInCell(cellId) = data[counter++];
        a_noEllipsoidsInCell(cellId) = data[counter++];
      }
    } else if(dataId == 2) {
      /*
      MInt i = 0;
      for(MInt part = 0; part < a_noParticles(); part++){
        MInt a = data[i];
        MInt b = data[i + 1 ];
        MLong partId = ((MLong)a << 32 | ((MLong)b & 0xFFFFFFFFL));
        //m_partList[part].m_partId = partId;
        //m_partList[part].m_noParticles = data[i + 2];
        //ASSERT(m_partList[part].m_noParticles > 0 , "");
        if(data[i + 2] < 0 || data[i + 2] > 10000) {
          cerr << "Strange parcel-count!" << data[i + 2] << " " << partId << endl;
        }
        i = i + 3;
      }
      */
    }
  }
}


/**
 * \brief Get global solver variables that need to be communicated for DLB
 *        (same on each domain)
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::getGlobalSolverVars(vector<MFloat>& globalFloatVars, vector<MInt>& globalIntVars) {
  TRACE();

  globalIntVars.push_back(m_PRNGSpawnCount);

  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    globalIntVars.push_back(m_sprayModel->m_PRNGPrimBreakUpCount);
  }

  globalFloatVars.push_back(m_particleResiduum);
  globalFloatVars.push_back(m_time);
  globalFloatVars.push_back(m_timeStep);

  if(m_activePrimaryBUp) {
    globalFloatVars.push_back(m_sprayModel->timeSinceSOI());
  }

  if(m_activePrimaryBUp && !m_injData.empty()) {
    for(auto it = m_injData.begin(); it != m_injData.end(); it++) {
      const MFloat timeStep = (MFloat)it->first;
      globalFloatVars.push_back(timeStep);
      for(MInt i = 0; i < m_sprayModel->m_injDataSize + m_addInjDataCnt; i++) {
        const MFloat data = (it->second)[i];
        globalFloatVars.push_back(data);
      }
    }
    globalIntVars.push_back(m_injData.size());
  } else if(m_activePrimaryBUp) {
    globalIntVars.push_back(0);
  }
}

/**
 * \brief Set global solver variables for DLB (same on each domain)
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::setGlobalSolverVars(vector<MFloat>& globalFloatVars, vector<MInt>& globalIntVars) {
  TRACE();

  m_PRNGSpawnCount = globalIntVars[0];

  m_particleResiduum = globalFloatVars[0];
  m_time = globalFloatVars[1];
  m_timeStep = globalFloatVars[2];

  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    m_sprayModel->m_PRNGPrimBreakUpCount = globalIntVars[1];
  }
  if(m_activePrimaryBUp) {
    m_sprayModel->timeSinceSOI() = globalFloatVars[3];
  }


  if(m_activePrimaryBUp) {
    const MInt vectorSize = globalIntVars[2];
    m_injData.clear();

    const MInt count = m_sprayModel->m_injDataSize + m_addInjDataCnt;
    for(MInt i = 0; i < vectorSize * (count + 1); i = i + 1 + count) {
      const MInt timeStep = (MInt)globalFloatVars[i + m_addInjDataCnt];
      vector<MFloat> tmp(count);
      for(MInt j = 0; j < count; j++) {
        tmp[j] = globalFloatVars[i + j + 1 + m_addInjDataCnt];
      }
      m_injData.insert(make_pair(timeStep, tmp));
    }
  }
}


/**
 * \brief re-inits the solver after the balance
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::finalizeBalance() {
  TRACE();

  if(!isActive()) return;

  grid().findCartesianNghbIds();

  findSpawnCellId();

  // reset random-number generator
  if(m_activePrimaryBUp || m_activeSecondaryBUp) {
    m_sprayModel->m_PRNGPrimBreakUp.seed(m_sprayModel->m_spraySeed);
  }

  if(domainId() != m_spawnDomainId) {
    m_particleResiduum = 0;
    if(m_activePrimaryBUp || m_activeSecondaryBUp) {
      m_sprayModel->m_PRNGPrimBreakUpCount = 0;
    }
  } else {
    m_sprayModel->m_PRNGPrimBreakUp.discard(m_sprayModel->m_PRNGPrimBreakUpCount);
  }

  m_cellToNghbrHood.clear();
  m_cellToNghbrHoodInterpolation.clear();

  // removeInvalidParticles(true);

  checkParticles();
  checkCells();

  updateFluidFraction();

  setRegridTrigger();

  // reset injection data on all ranks which do not have the injection-cell
  if(m_spawnCellId == -1 && m_sprayModel != nullptr) {
    const MInt vectorSize = m_injData.size();
    m_injData.clear();
    const MInt count = m_sprayModel->m_injDataSize + m_addInjDataCnt;
    vector<MFloat> tmp(count);
    for(MInt i = 0; i < count; i++) {
      tmp[i] = 0.0;
    }
    for(MInt i = 0; i < vectorSize; i++) {
      MInt stepId = globalTimeStep - vectorSize + i + 1;
      m_injData.insert(make_pair(stepId, tmp));
    }
  }

#ifdef LPT_DEBUG
  MInt noParticles = a_noParticles();

  MPI_Allreduce(MPI_IN_PLACE, &noParticles, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noParticles");

  if(domainId() == 0) {
    cerr << noParticles << " #particles after balancing!" << endl;
  }

  if(m_ellipsoids) {
    MInt noEllipsoids = a_noEllipsoidalParticles();
    MPI_Allreduce(MPI_IN_PLACE, &noEllipsoids, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noEllipsoids");
    cerr0 << noEllipsoids << " #ellipsoids after balancing!" << endl;
  }
#endif

  /*
  MInt backup = globalTimeStep;
  globalTimeStep = -2;
  saveDebugRestartFile();
  globalTimeStep = backup;
  */

  if(globalTimeStep < 0) {
    // add refinement around spawnCellId!
    // NOTE: a_noParticlesInCell will be corrected after the adaptation!
    if(m_spawnCellId > -1) {
      a_noParticlesInCell(m_spawnCellId) += 1;
      if(m_ellipsoids) a_noEllipsoidsInCell(m_spawnCellId) += 1;
    }
  } else {
    // get new neighbor list and distribution weights
    for(MInt i = 0; i < a_noParticles(); i++) {
      m_partList[i].resetWeights();
    }
    for(MInt i = 0; i < a_noEllipsoidalParticles(); i++) {
      m_partListEllipsoid[i].resetWeights();
    }
    // re-compute source terms
    resetSourceTerms(0);
    m_sumEvapMass = 0;
    if(!m_skipLPT) {
      coupling(-1);
      m_sumEvapMass = 0;
      if(m_nonBlockingComm) {
        sendSourceTerms();
        waitForSendReqs();
        receiveSourceTerms();
      } else {
        exchangeSourceTerms();
      }
    }
  }

  // LPT-statistics
  MInt noCells = a_noCells();
  MInt minCells = noCells;
  MInt maxCells = noCells;
  MInt noPart = a_noParticles();
  MInt minParts = noPart;
  MInt maxParts = noPart;

  MPI_Allreduce(MPI_IN_PLACE, &noCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noCells");
  MPI_Allreduce(MPI_IN_PLACE, &noPart, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "noPart");
  MPI_Allreduce(MPI_IN_PLACE, &maxCells, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxCells");
  MPI_Allreduce(MPI_IN_PLACE, &minCells, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minCells");
  MPI_Allreduce(MPI_IN_PLACE, &maxParts, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "INPLACE", "maxParts");
  MPI_Allreduce(MPI_IN_PLACE, &minParts, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "INPLACE", "minParts");

  if(domainId() == 0) {
    cerr << "LPT cell-statistics : " << minCells << " - " << maxCells << " avg " << noCells / noDomains() << endl;
    cerr << "LPT particle-statistics : " << minParts << " - " << maxParts << " avg " << noPart / noDomains() << endl;
  }
}

/**
 * \brief  set a_regridTriggerG around all particles to allow for an adaptation check
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::setRegridTrigger() {
  TRACE();

  if(!this->m_adaptation || this->m_noSensors == 0) return;

  m_log << "Reset LPT-regrid trigger at timestep: " << globalTimeStep << endl;

  MBoolScratchSpace inShadowLayer(a_noCells(), AT_, "inShadowLayer");
  MBoolScratchSpace tmp_data(a_noCells(), AT_, "tmp_data");

  vector<MInt> lastLayer;
  vector<MInt> tmp;

  // reset property and mark cells with particles as first layer
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_regridTrigger(cellId) = false;
    if(a_noParticlesInCell(cellId) > 0 || a_noEllipsoidsInCell(cellId) > 0 || cellId == m_spawnCellId) {
      inShadowLayer(cellId) = true;
      lastLayer.emplace_back(cellId);
    }
  }

  // 4. cover shadow band locally in the domain
  MInt layerCount = 0;
  while(layerCount < m_bandWidth[maxRefinementLevel() - 1]) {
    // ii. build the next layer
    tmp.clear();
    for(MInt c = 0; c < (signed)lastLayer.size(); c++) {
      const MInt cellId = lastLayer[c];
      for(MInt d = 0; d < m_noDirs; d++) {
        if(c_hasNeighbor(cellId, d) == 0) {
          // additionally check for regrid when rebuilding
          if(c_parentId(cellId) > -1 && c_hasNeighbor(c_parentId(cellId), d) && !a_isHalo(cellId)) {
            if(a_isValidCell(cellId)) {
              m_forceAdaptation = true;
            }
          }
          continue;
        }
        const MInt nghbrId = c_neighborId(cellId, d);
        if(!inShadowLayer[nghbrId]) {
          tmp.emplace_back(nghbrId);
          if(layerCount > m_bandWidth[maxRefinementLevel() - 1] - m_innerBound) {
            a_regridTrigger(nghbrId) = true;
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

    this->exchangeData(&tmp_data(0), 1);

    for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
      for(MInt j = 0; j < noHaloCells(i); j++) {
        const MInt haloCell = haloCellId(i, j);
        if(!inShadowLayer[haloCell]) {
          inShadowLayer[haloCell] = tmp_data(haloCell);
          if(inShadowLayer[haloCell]) {
            tmp.emplace_back(haloCell);
            if(layerCount > m_bandWidth[maxRefinementLevel() - 1] - m_innerBound) {
              a_regridTrigger(haloCell) = true;
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
    layerCount++;
  }


  MPI_Allreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, mpiComm(), AT_, "MPI_IN_PLACE", "status");

  if(domainId() == 0 && m_forceAdaptation) {
    cerr << "LPT-Solver is forcing a mesh-adaptation at time step " << globalTimeStep << endl;
  }
}


/**
 * \brief  check if any of the cells with the regrid Trigger also have particles inside
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::checkRegridTrigger() {
  if(!this->m_adaptation || this->m_noSensors < 1) return;

  for(MInt id = 0; id < a_noParticles(); id++) {
    if(m_partList[id].isInvalid()) continue;
    const MInt cellId = m_partList[id].m_cellId;
    if(a_isHalo(cellId)) continue;
    if(a_regridTrigger(cellId)) {
      m_forceAdaptation = true;
      break;
    }
  }
  for(MInt id = 0; id < a_noEllipsoidalParticles(); id++) {
    if(m_partListEllipsoid[id].isInvalid()) continue;
    const MInt cellId = m_partListEllipsoid[id].m_cellId;
    if(a_isHalo(cellId)) continue;
    if(a_regridTrigger(cellId)) {
      m_forceAdaptation = true;
      break;
    }
  }

  if(!m_nonBlockingComm) {
    // moved exchange for forceAdaptation() call!
    MPI_Allreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, mpiComm(), AT_, "MPI_IN_PLACE", "status");

    if(domainId() == 0 && m_forceAdaptation) {
      cerr << "LPT-Solver is forcing a mesh-adaptation at time step " << globalTimeStep << endl;
    }
  } else {
    MPI_Iallreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, mpiComm(),
                   &m_checkAdaptationRecvRequest[0], AT_, "MPI_IN_PLACE", "status");
    m_openRegridSend = true;
  }
}

template <MInt nDim>
void LPT<nDim>::recvRegridTrigger() {
  if(!m_nonBlockingComm) {
    return;
  }
  if(!isActive()) {
    return;
  }

  // wait for all send and receive
  if(m_openRegridSend) {
    RECORD_TIMER_START(m_timers[Timers::Exchange]);
    RECORD_TIMER_START(m_timers[Timers::Exchange2]);
    MPI_Wait(&m_checkAdaptationRecvRequest[0], MPI_STATUS_IGNORE, AT_);
    m_openRegridSend = false;
    RECORD_TIMER_STOP(m_timers[Timers::Exchange2]);
    RECORD_TIMER_STOP(m_timers[Timers::Exchange]);
  }
}

/**
 * \brief  computes the non-dimensional LPT time-step on the assumption
 *         that a particle may only travel
 *         for a complete cell length with cfl = 1 at current velocity!
 * \author Tim Wegmann
 */
template <MInt nDim>
MFloat LPT<nDim>::computeTimeStep() {
  TRACE();

  m_timeStepUpdated = true;

  MFloat timeStep = std::numeric_limits<MFloat>::max();

  for(MInt id = 0; id < a_noParticles(); id++) {
    if(m_partList[id].isInvalid()) continue;
    const MFloat dx = c_cellLengthAtCell(m_partList[id].m_cellId) / m_partList[id].s_lengthFactor;
    for(MInt d = 0; d < nDim; ++d) {
      const MFloat u_d = fabs(m_partList[id].m_velocity[d]);
      MFloat dt_dir = 0;
      if(abs(m_cfl + 1.0) < MFloatEps && abs(m_Ma + 1.0) > MFloatEps) {
        dt_dir = m_Ma / SQRT3 * dx / u_d;
      } else if(abs(m_cfl + 1.0) > MFloatEps)
        dt_dir = m_cfl * dx / u_d;
      timeStep = mMin(timeStep, dt_dir);
    }
  }

  // Ellipsoids
  for(MInt id = 0; id < a_noEllipsoidalParticles(); id++) {
    if(m_partListEllipsoid[id].isInvalid()) continue;
    const MFloat dx = c_cellLengthAtCell(m_partListEllipsoid[id].m_cellId) / m_partListEllipsoid[id].s_lengthFactor;
    for(MInt d = 0; d < nDim; ++d) {
      const MFloat u_d = fabs(m_partListEllipsoid[id].m_velocity[d]);
      MFloat dt_dir = 0;
      if(abs(m_cfl + 1.0) < MFloatEps && abs(m_Ma + 1.0) > MFloatEps) {
        dt_dir = m_Ma / SQRT3 * dx / u_d;
      } else if(abs(m_cfl + 1.0) > MFloatEps)
        dt_dir = m_cfl * dx / u_d;
      timeStep = mMin(timeStep, dt_dir);
    }
  }

  // consider injection speed
  if(m_activePrimaryBUp && m_spawnCellId > -1) {
    const MFloat dx = c_cellLengthAtCell(maxRefinementLevel()) / m_partList[0].s_lengthFactor;
    const MFloat u_d = fabs(m_sprayModel->injectionSpeed());
    MFloat dt_dir = m_cfl * dx / u_d;
    timeStep = mMin(timeStep, dt_dir);
  }

  // consider spawn velocity
  if(m_spawnParticles && m_spawnCellId > -1) {
    const MFloat dx = c_cellLengthAtCell(m_spawnCellId) / m_partList[0].s_lengthFactor;
    const MFloat u_d = fabs(m_spawnVelocity);
    MFloat dt_dir = m_cfl * dx / u_d;
    timeStep = mMin(timeStep, dt_dir);
  }

  // read timeStep from property file
  // i.e. used to resolve evaporation time-span
  MBool useFixed = false;
  if(a_noParticles() > 0) useFixed = true;
  if(a_noEllipsoidalParticles() > 0) useFixed = true;
  if(m_spawnCellId > -1 && m_sprayModel != nullptr && m_sprayModel->soonInjection(m_time + timeStep * 100)) {
    // check if injection will occure in the next time-steps
    useFixed = true;
  }

  if(useFixed && Context::propertyExists("fixedTimeStep", solverId())) {
    timeStep = Context::getSolverProperty<MFloat>("fixedTimeStep", m_solverId, AT_);
    if(domainId() == m_spawnDomainId) {
      cerr << "Applying LPT fixed timeStep of " << timeStep << endl;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &timeStep, 1, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "timeStep");

  if(domainId() == 0) {
    cerr << "LPT-timestep could be " << timeStep << " is " << m_timeStep << endl;
  }

  return timeStep;
}

/**
 * \brief  exchanges the source terms on the LPT grid, this is necessary
 *         a) before writing the cell solution file for consistent/rank-independant data
 *         b) before the adaptation, so that all relevant data is stored on window-cells
 *         c) before transferring the source terms to the coupled solver, so that all relevant data
 *            is stored on window cells
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::exchangeSourceTerms() {
  if(m_nonBlockingComm) return;
  if(noDomains() > 1) {
    if(m_massCoupling) {
      maia::mpi::reverseExchangeAddData(grid().neighborDomains(), grid().haloCells(), grid().windowCells(), mpiComm(),
                                        &a_massFlux(0), 1);
    }
    if(m_momentumCoupling) {
      maia::mpi::reverseExchangeAddData(grid().neighborDomains(), grid().haloCells(), grid().windowCells(), mpiComm(),
                                        &a_momentumFlux(0, 0), nDim);
      maia::mpi::reverseExchangeAddData(grid().neighborDomains(), grid().haloCells(), grid().windowCells(), mpiComm(),
                                        &a_workFlux(0), 1);
    }
    if(m_heatCoupling) {
      maia::mpi::reverseExchangeAddData(grid().neighborDomains(), grid().haloCells(), grid().windowCells(), mpiComm(),
                                        &a_heatFlux(0), 1);
    }
  }
  // after the values on halo-cells have been added to the matching window cells, the halo-cell values
  // must be resettet
  resetSourceTerms(noInternalCells());
}

/**
 * \brief  exchanges the source terms on the LPT grid, non-Blocking version from above
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::sendSourceTerms() {
  if(!m_nonBlockingComm) {
    mTerm(1, AT_, "Calling non-Blocking exchange for blocking setup!");
  }

  // start receiving source terms
  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    MPI_Irecv(&m_sourceRecv[d][0], m_noSourceTerms * grid().noLeafWindowCells(d), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("SOURCE_TERMS"), mpiComm(), &m_sourceRecvRequest[n], AT_, "m_sourceRecv");
  }

  // wait until all sends have been received before sending new
  if(m_openSourceSend) {
    MPI_Waitall(grid().noLeafRecvNeighborDomains(), &m_sourceSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openSourceSend = false;
  }

  // send source terms
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafHaloCells(d); j++) {
      const MInt cellId = grid().leafHaloCell(d, j);
      if(m_massCoupling) {
        m_sourceSend[d][buffId++] = a_massFlux(cellId);
      }
      if(m_momentumCoupling) {
        for(MInt i = 0; i < nDim; i++) {
          m_sourceSend[d][buffId++] = a_momentumFlux(cellId, i);
        }
        m_sourceSend[d][buffId++] = a_workFlux(cellId);
      }
      if(m_heatCoupling) {
        m_sourceSend[d][buffId++] = a_heatFlux(cellId);
      }
    }
    MPI_Isend(&m_sourceSend[d][0], m_noSourceTerms * grid().noLeafHaloCells(d), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("SOURCE_TERMS"), mpiComm(), &m_sourceSendRequest[n], AT_, "m_sourceSend");
  }
  m_openSourceSend = true;

  // reset sources on halo cells to zero after send
  resetSourceTerms(noInternalCells());
}

/**
 * \brief  exchanges the source terms on the LPT grid, non-Blocking version from above
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::receiveSourceTerms() {
  if(!m_nonBlockingComm) return;

  RECORD_TIMER_START(m_timers[Timers::Exchange]);

  // wait for all receives before using buffers
  MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_sourceRecvRequest[0], MPI_STATUSES_IGNORE, AT_);

  // add received source terms
  RECORD_TIMER_START(m_timers[Timers::Exchange5]);

  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafWindowCells(d); j++) {
      const MInt cellId = grid().leafWindowCell(d, j);
      if(m_massCoupling) {
        a_massFlux(cellId) += m_sourceRecv[d][buffId++];
      }
      if(m_momentumCoupling) {
        for(MInt i = 0; i < nDim; i++) {
          a_momentumFlux(cellId, i) += m_sourceRecv[d][buffId++];
        }
        a_workFlux(cellId) += m_sourceRecv[d][buffId++];
      }
      if(m_heatCoupling) {
        a_heatFlux(cellId) += m_sourceRecv[d][buffId++];
      }
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Exchange5]);
  RECORD_TIMER_STOP(m_timers[Timers::Exchange]);
}

/**
 * \brief  exchanges the flow Field data on the LPT grid, non-Blocking version
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::sendFlowField() {
  if(!m_nonBlockingComm) {
    return;
  }
  m_receiveFlowField = true;

  // start receiving flow field data
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MPI_Irecv(&m_flowRecv[d][0], (PV.noVars() + m_evaporation) * grid().noLeafHaloCells(d), MPI_DOUBLE,
              grid().neighborDomain(d), mpiTag("FLOW_FIELD"), mpiComm(), &m_flowRecvRequest[n], AT_, "m_flowRecv");
  }

  // wait until all sends have been processed before sending new
  if(m_openFlowSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_flowSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openFlowSend = false;
  }


  // send flow field data
  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafWindowCells(d); j++) {
      const MInt cellId = grid().leafWindowCell(d, j);
      std::copy_n(&a_fluidVariable(cellId, 0), PV.noVars(), &m_flowSend[d][buffId]);
      buffId += PV.noVars();
      if(m_evaporation) {
        std::copy_n(&a_fluidSpecies(cellId), 1, &m_flowSend[d][buffId]);
        buffId++;
      }
    }
    MPI_Isend(&m_flowSend[d][0], (PV.noVars() + m_evaporation) * grid().noLeafWindowCells(d), MPI_DOUBLE,
              grid().neighborDomain(d), mpiTag("FLOW_FIELD"), mpiComm(), &m_flowSendRequest[n], AT_, "m_flowSend");
  }

  m_openFlowSend = true;
}

/**
 * \brief  exchanges the flow Field data on the LPT grid, non-Blocking version
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::receiveFlowField() {
  if(!m_nonBlockingComm) return;
  if(!m_receiveFlowField) return;

  RECORD_TIMER_START(m_timers[Timers::Exchange]);

  MPI_Waitall(grid().noLeafRecvNeighborDomains(), &m_flowRecvRequest[0], MPI_STATUSES_IGNORE, AT_);

  // set data on halo-cells
  RECORD_TIMER_START(m_timers[Timers::Exchange5]);
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafHaloCells(d); j++) {
      const MInt cellId = grid().leafHaloCell(d, j);
      std::copy_n(&m_flowRecv[d][buffId], PV.noVars(), &a_fluidVariable(cellId, 0));
      buffId += PV.noVars();
      if(m_evaporation) {
        std::copy_n(&m_flowRecv[d][buffId], 1, &a_fluidSpecies(cellId));
        buffId++;
      }
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Exchange5]);
  RECORD_TIMER_STOP(m_timers[Timers::Exchange]);
  m_receiveFlowField = false;
}

/**
 * \brief  wait on all mpi send calls, prior to balance, adaptation or IO
 * \author Tim Wegmann
 */
template <MInt nDim>
void LPT<nDim>::waitForSendReqs() {
  if(!m_nonBlockingComm) {
    return;
  }

  if(m_openParticleInjSend) {
    MPI_Waitall(noNeighborDomains(), &m_mpi_reqSendFloat[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(noNeighborDomains(), &m_mpi_reqSendInt[0], MPI_STATUSES_IGNORE, AT_);
    m_openParticleInjSend = false;
  }

  if(m_openParticleSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_mpi_reqSendFloat[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_mpi_reqSendInt[0], MPI_STATUSES_IGNORE, AT_);
    m_openParticleSend = false;
  }

  if(m_openSourceSend) {
    MPI_Waitall(grid().noLeafRecvNeighborDomains(), &m_sourceSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openSourceSend = false;
  }

  if(m_openFlowSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_flowSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openFlowSend = false;
  }

  if(m_openRegridSend) {
    MPI_Wait(&m_checkAdaptationRecvRequest[0], MPI_STATUS_IGNORE, AT_);
    m_openRegridSend = false;
  }

  if(m_openSlopesSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_slopesSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openSlopesSend = false;
  }
}

/**
 * \brief  exchanges the velocity slopes on the LPT grid, non-Blocking version
 * \author Laurent Andre
 */
template <MInt nDim>
void LPT<nDim>::sendVelocitySlopes() {
  if(!m_ellipsoids) return;
  if(!m_nonBlockingComm) return;

  m_receiveVelocitySlopes = true;

  // start receiving velocity slopes
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MPI_Irecv(&m_slopesRecv[d][0], (nDim * nDim) * grid().noLeafHaloCells(d), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("VELOCITY_SLOPES"), mpiComm(), &m_slopesRecvRequest[n], AT_, "m_slopesRecv");
  }

  // wait until all sends have been processed before sending new
  if(m_openSlopesSend) {
    MPI_Waitall(grid().noLeafSendNeighborDomains(), &m_slopesSendRequest[0], MPI_STATUSES_IGNORE, AT_);
    m_openSlopesSend = false;
  }

  // send velocity slopes
  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafWindowCells(d); j++) {
      const MInt cellId = grid().leafWindowCell(d, j);
      std::copy_n(&a_velocitySlope(cellId, 0, 0), nDim * nDim, &m_slopesSend[d][buffId]);
      buffId += nDim * nDim;
    }
    MPI_Isend(&m_slopesSend[d][0], (nDim * nDim) * grid().noLeafWindowCells(d), MPI_DOUBLE, grid().neighborDomain(d),
              mpiTag("VELOCITY_SLOPES"), mpiComm(), &m_slopesSendRequest[n], AT_, "m_slopesSend");
  }

  m_openSlopesSend = true;
}

/**
 * \brief  exchanges the velocity slopes on the LPT grid, non-Blocking version
 * \author Laurent Andre
 */
template <MInt nDim>
void LPT<nDim>::receiveVelocitySlopes() {
  if(!m_ellipsoids) return; // only if ellipsoids are activated
  if(!m_nonBlockingComm) return;
  if(!m_receiveVelocitySlopes) return;

  RECORD_TIMER_START(m_timers[Timers::Exchange]);

  MPI_Waitall(grid().noLeafRecvNeighborDomains(), &m_slopesRecvRequest[0], MPI_STATUSES_IGNORE, AT_);

  // set data on halo-cells
  RECORD_TIMER_START(m_timers[Timers::Exchange5]);
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MInt buffId = 0;
    for(MInt j = 0; j < grid().noLeafHaloCells(d); j++) {
      const MInt cellId = grid().leafHaloCell(d, j);
      std::copy_n(&m_slopesRecv[d][buffId], nDim * nDim, &a_velocitySlope(cellId, 0, 0));
      buffId += nDim * nDim;
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Exchange5]);
  RECORD_TIMER_STOP(m_timers[Timers::Exchange]);

  m_receiveVelocitySlopes = false;
}

/**
 * \brief save a full solutionOutput at timeSteps closesed to a specified crankAngle Interval
 *        save a solution slice at timeSteps closesed to a specified slice Interval
 * \author Tim Wegmann
 *
 */
template <MInt nDim>
void LPT<nDim>::crankAngleSolutionOutput() {
  TRACE();

  if(m_sprayModel == nullptr) return;

  static const MInt cadStart = m_sprayModel->m_injectionCrankAngle;
  static const MInt cadEnd = m_sprayModel->m_injectionCrankAngle + 90;
  static const MInt cadInterval = 1;
  const MInt cadMaxIter = (cadEnd - cadStart) / cadInterval;

  const MFloat cad = maia::math::crankAngle(m_time, m_sprayModel->m_Strouhal, m_sprayModel->m_initialCad, 0);
  const MFloat cad_prev =
      maia::math::crankAngle(m_time - m_timeStep, m_sprayModel->m_Strouhal, m_sprayModel->m_initialCad, 0);
  const MFloat cad_next =
      maia::math::crankAngle(m_time + m_timeStep, m_sprayModel->m_Strouhal, m_sprayModel->m_initialCad, 0);

  if(cad > cadEnd) return;
  if(cad < (cadStart - cadInterval)) return;

  // iterate for full solution output
  for(MInt i = 0; i < cadMaxIter; i++) {
    const MFloat cadTarget = cadStart + i * cadInterval;

    if(cad < cadTarget && cad_next > cadTarget) {
      if(fabs(cad - cadTarget) <= fabs(cad_next - cadTarget)) {
        if(domainId() == 0) {
          cerr << "Saving particle output at crankAngle " << cad << endl;
        }
        writePartData();
        break;
      }
    } else if(cad > cadTarget && cad_prev < cadTarget) {
      if(abs(cad - cadTarget) < fabs(cad_prev - cadTarget)) {
        if(domainId() == 0) {
          cerr << "Saving particle output at crankAngle " << cad << endl;
        }
        writePartData();
        break;
      }
    }
  }
}

/**
 * \brief Initializes the communication timers
 * \author Tim Wegmann
 *
 */
template <MInt nDim>
void LPT<nDim>::initializeTimers() {
  TRACE();

  std::fill(m_timers.begin(), m_timers.end(), -1);

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "LPT (solverId = " + to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::SolverType], "total object lifetime", m_timerGroup);
  RECORD_TIMER_START(m_timers[Timers::SolverType]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TimeInt], "Time-integration", m_timers[Timers::SolverType]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PreTime], "PreTimeStep", m_timers[Timers::SolverType]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Motion], "Motion", m_timers[Timers::TimeInt]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Energy], "Energy", m_timers[Timers::TimeInt]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Wall], "Wall", m_timers[Timers::TimeInt]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SourceTerms], "Sources", m_timers[Timers::TimeInt]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Injection], "injection", m_timers[Timers::TimeInt]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Breakup], "breakup", m_timers[Timers::TimeInt]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange], "exchange", m_timers[Timers::SolverType]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange1], "mpi", m_timers[Timers::Exchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange2], "prepare", m_timers[Timers::Exchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange3], "pack", m_timers[Timers::Exchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange4], "receive", m_timers[Timers::Exchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Exchange5], "unpack", m_timers[Timers::Exchange]);
}

/**
 * \brief reduce data to lower level
 * \author Tim Wegmann
 *
 */
template <MInt nDim>
MFloat LPT<nDim>::reduceData(const MInt cellId, MFloat* data, const MInt dataBlockSize, const MBool average) {
  MFloat vol = a_bndryCellId(cellId) < 0 ? c_cellVolume(cellId) : m_bndryCells->a[a_bndryCellId(cellId)].m_volume;

  if(c_noChildren(cellId) > 0) {
    vol = F0;
    for(MInt d = 0; d < dataBlockSize; d++) {
      data[dataBlockSize * cellId + d] = F0;
    }
    for(MInt child = 0; child < c_noChildren(cellId); child++) {
      MInt childId = c_childId(cellId, child);
      if(childId < 0) {
        continue;
      }
      if(c_noChildren(childId) == 0 && !a_isValidCell(childId)) {
        continue;
      }
      MFloat volc = reduceData(childId, data, dataBlockSize);
      for(MInt d = 0; d < dataBlockSize; d++) {
        data[dataBlockSize * cellId + d] += volc * data[dataBlockSize * childId + d];
      }
      vol += volc;
    }
    if(average) {
      for(MInt d = 0; d < dataBlockSize; d++) {
        data[dataBlockSize * cellId + d] /= mMax(1e-14, vol);
      }
    }
  }
  return vol;
}

/*! \brief calculates interpolated variables (in the range a, b) for a given position in a given cell
 *
 * \author Sven Berger, update Laurent Andre
 * \date January 2016, August 2022
 *
 * \param[in] cellId id of the cell used for interpolation
 * \param[in] position interpolation point
 * \param[out] interpolatedVar array containing the interpolated variables
 *
 */
template <MInt nDim>
template <MInt a, MInt b, MBool interpolateVelocitySlopes>
void LPT<nDim>::interpolateVariablesLS(const MInt cellId, const MFloat* const position, MFloat* const interpolatedVar) {
  ASSERT((b - a) > 0, "ERROR: Difference between b and a needs to be positive!");
  ASSERT(cellId >= 0, "ERROR: Invalid cellId!");

  const MFloat originX = c_coordinate(cellId, 0);
  const MFloat originY = c_coordinate(cellId, 1);
  const MFloat originZ = c_coordinate(cellId, 2);

  ///
  /// 1. Find all direct neighboring cells
  ///
  std::vector<MInt> nghbrList;
#ifdef _OPENMP
#pragma omp critical
#endif
  {
    if(m_cellToNghbrHoodInterpolation.count(cellId) > 0) {
      nghbrList = m_cellToNghbrHoodInterpolation[cellId];
      // std::cout << "New use existing hood with " << nghbrList.size() << " for " << cellId << std::endl;
    } else {
      grid().findDirectNghbrs(cellId, nghbrList);
      // std::cout << "New find hood with " << nghbrList.size() << std::endl;
      if(nghbrList.size() < 12) {
        // std::cout << "New find better hood with " << nghbrList.size() << std::endl;
        grid().findNeighborHood(cellId, 2, nghbrList);
#ifndef NDEBUG
        if(nghbrList.size() < 14) {
          std::cout << "WARNING: Only found " << nghbrList.size() << " neighbors during interpolation!" << std::endl;
        }
#endif
      }
      m_cellToNghbrHoodInterpolation.emplace(make_pair(cellId, nghbrList));
    }
  }

  const MInt pseudoVarCount = (interpolateVelocitySlopes ? nDim * nDim + b - a : b - a);
  std::vector<std::array<MFloat, pseudoVarCount>> pseudoCellVar{};
  std::vector<std::array<MFloat, nDim>> pseudoCellPos{};

  auto nghbrIt = nghbrList.begin();
  while(nghbrIt != nghbrList.end()) {
    const MInt nghbrId = *nghbrIt;
    const MInt bndryCellId = a_bndryCellId(nghbrId);
    if(bndryCellId > -1) {
      const auto& bndryCell = m_bndryCells->a[bndryCellId];
      for(MInt srfcId = 0; srfcId < bndryCell.m_noSrfcs; srfcId++) {
        auto& srfc = *bndryCell.m_srfcs[srfcId];
        pseudoCellPos.emplace_back();
        auto& pos = pseudoCellPos.back();
        for(MInt n = 0; n < nDim; n++) {
          pos[n] = srfc.m_coordinates[n];
        }
        // TODO labels:LPT Distinguish between different boundary condtions
        pseudoCellVar.emplace_back();
        auto& vars = pseudoCellVar.back();
        for(MInt i = a; i < b; i++) {
          vars[i] = a_fluidVariable(nghbrId, i);
        }
        if(interpolateVelocitySlopes) {
          for(MInt i = 0; i < nDim; i++) {
            for(MInt j = 0; j < nDim; j++) {
              vars[b + nDim * i + j] = a_velocitySlope(nghbrId, i, j);
            }
          }
        }
      }
      nghbrIt = nghbrList.erase(nghbrIt);
    } else {
      ++nghbrIt;
    }
  }

  const MUint noInterpolationPoints = nghbrList.size() + pseudoCellPos.size();

  ///
  /// 2a. if only a low number of neighbors is found use distance weighted averaging
  ///
  if(noInterpolationPoints < 10) {
    std::vector<MFloat> dist(noInterpolationPoints, 0);
    // Internal cells
    for(MUint nId = 0; nId < nghbrList.size(); nId++) {
      const MInt nghbrId = nghbrList[nId];
      dist[nId] = 0;
      for(MInt i = 0; i < nDim; i++) {
        dist[nId] += POW2(c_coordinate(nghbrId, i) - position[i]);
      }
    }
    // Boundary surfaces
    for(MUint pId = 0; pId < pseudoCellPos.size(); pId++) {
      IF_CONSTEXPR(nDim == 2) {
        dist[pId + nghbrList.size()] =
            POW2(pseudoCellPos[pId][0] - position[0]) + POW2(pseudoCellPos[pId][1] - position[1]);
      }
      else {
        dist[pId + nghbrList.size()] = POW2(pseudoCellPos[pId][0] - position[0])
                                       + POW2(pseudoCellPos[pId][1] - position[1])
                                       + POW2(pseudoCellPos[pId][2] - position[2]);
      }
    }

    MFloat sum = std::accumulate(dist.begin(), dist.end(), 0.0);

    for(MInt k = a; k < b; k++) {
      MFloat var = 0;
      interpolatedVar[k] = 0.0;

      // Internal cells
      for(MUint nId = 0; nId < nghbrList.size(); nId++) {
        const MInt nghbrId = nghbrList[nId];
        var = a_fluidVariable(nghbrId, k);
        interpolatedVar[k] += dist[nId] / sum * var;
      }
      // Boundary surfaces
      for(MUint pId = 0; pId < pseudoCellPos.size(); pId++) {
        var = pseudoCellVar[pId][k - a];
        interpolatedVar[k] += dist[pId + nghbrList.size()] / sum * var;
      }
    }

    // interpolate velocity slopes
    if(interpolateVelocitySlopes) {
      for(MInt i = 0; i < nDim; i++) {
        for(MInt j = 0; j < nDim; j++) {
          MFloat var = 0;
          interpolatedVar[b + nDim * i + j - a] = 0.0;
          // Internal cells
          for(MUint nId = 0; nId < nghbrList.size(); nId++) {
            const MInt nghbrId = nghbrList[nId];
            var = a_velocitySlope(nghbrId, i, j);
            interpolatedVar[b + nDim * i + j - a] += dist[nId] / sum * var;
          }
          // Boundary surfaces
          for(MUint pId = 0; pId < pseudoCellPos.size(); pId++) {
            var = pseudoCellVar[pId][nDim * i + j + b];
            interpolatedVar[b + nDim * i + j - a] += dist[pId + nghbrList.size()] / sum * var;
          }
        }
      }
    }

#ifndef NDEBUG
    //    cout << "WARNING: Not enough neighbors found to use LS using weighted averaging instead!" << endl;
    m_log << "WARNING: Not enough neighbors found to use LS using weighted averaging instead!" << std::endl;
#endif
  } else {
    ///
    /// 2b. build the least square matrix LSMatrix...
    ///
    MFloatTensor LSMatrix(2 * nDim + pow(2, nDim - 1), 2 * nDim + pow(2, nDim - 1));
    MFloatTensor rhs(2 * nDim + pow(2, nDim - 1));
    //    MFloatTensor weights(2 * nDim + pow(2, nDim - 1));
    MFloatTensor matInv(2 * nDim + pow(2, nDim - 1), 2 * nDim + pow(2, nDim - 1));

    LSMatrix.set(0.0);
    //    weights.set(1.0);
    matInv.set(0.0);

    MFloat sumX = 0;
    MFloat sumY = 0;
    MFloat sumZ = 0;
    MFloat sumXY = 0;
    MFloat sumXZ = 0;
    MFloat sumYZ = 0;
    MFloat sumX2 = 0;
    MFloat sumY2 = 0;
    MFloat sumZ2 = 0;
    MFloat sumXYZ = 0;
    MFloat sumX2Y = 0;
    MFloat sumX2Z = 0;
    MFloat sumXY2 = 0;
    MFloat sumY2Z = 0;
    MFloat sumXZ2 = 0;
    MFloat sumYZ2 = 0;
    MFloat sumX3 = 0;
    MFloat sumY3 = 0;
    MFloat sumZ3 = 0;
    MFloat sumX3Y = 0;
    MFloat sumX3Z = 0;
    MFloat sumXY3 = 0;
    MFloat sumY3Z = 0;
    MFloat sumXZ3 = 0;
    MFloat sumYZ3 = 0;
    MFloat sumX2Y2 = 0;
    MFloat sumX2Z2 = 0;
    MFloat sumX2YZ = 0;
    MFloat sumY2Z2 = 0;
    MFloat sumXY2Z = 0;
    MFloat sumXYZ2 = 0;
    MFloat sumX4 = 0;
    MFloat sumY4 = 0;
    MFloat sumZ4 = 0;


    // 2.1 Calculate intermediate variables
    for(MUint nId = 1; nId < noInterpolationPoints; nId++) {
      MFloat x = 0;
      MFloat y = 0;
      MFloat z = 0;
      if(nId >= nghbrList.size()) {
        const MUint pId = nId - nghbrList.size();
        x = pseudoCellPos[pId][0] - originX;
        y = pseudoCellPos[pId][1] - originY;
        IF_CONSTEXPR(nDim == 3) { z = pseudoCellPos[pId][2] - originZ; }
      } else {
        const MInt nghbrId = nghbrList[nId];
        x = c_coordinate(nghbrId, 0) - originX;
        y = c_coordinate(nghbrId, 1) - originY;
        IF_CONSTEXPR(nDim == 3) { z = c_coordinate(nghbrId, 2) - originZ; }
      }

      sumX += x;
      sumY += y;
      sumZ += z;
      sumXY += x * y;
      sumXZ += x * z;
      sumYZ += y * z;
      sumX2 += x * x;
      sumY2 += y * y;
      sumZ2 += z * z;
      sumXYZ += x * y * z;
      sumX2Y += x * x * y;
      sumX2Z += x * x * z;
      sumXY2 += x * y * y;
      sumY2Z += y * y * z;
      sumXZ2 += z * z * x;
      sumYZ2 += z * z * y;
      sumX3 += x * x * x;
      sumY3 += y * y * y;
      sumZ3 += z * z * z;
      sumX3Y += x * x * x * y;
      sumX3Z += x * x * x * z;
      sumXY3 += y * y * y * x;
      sumY3Z += y * y * y * z;
      sumXZ3 += z * z * z * x;
      sumYZ3 += z * z * z * y;
      sumX2Y2 += x * x * y * y;
      sumX2Z2 += x * x * z * z;
      sumX2YZ += x * x * y * z;
      sumY2Z2 += y * y * z * z;
      sumXY2Z += x * y * y * z;
      sumXYZ2 += x * y * z * z;
      sumX4 += x * x * x * x;
      sumY4 += y * y * y * y;
      sumZ4 += z * z * z * z;
    }

    // 2.2 Use intermediate variables to build LSMatrix
    IF_CONSTEXPR(nDim == 2) {
      // diagonal
      LSMatrix(0, 0) = sumX4;
      LSMatrix(1, 1) = sumY4;
      LSMatrix(2, 2) = sumX2Y2;
      LSMatrix(3, 3) = sumX2;
      LSMatrix(4, 4) = sumY2;
      LSMatrix(5, 5) = noInterpolationPoints;

      // first column/row
      LSMatrix(0, 1) = LSMatrix(1, 0) = sumX2Y2;
      LSMatrix(0, 2) = LSMatrix(2, 0) = sumX3Y;
      LSMatrix(0, 3) = LSMatrix(3, 0) = sumX3;
      LSMatrix(0, 4) = LSMatrix(4, 0) = sumX2Y;
      LSMatrix(0, 5) = LSMatrix(5, 0) = sumX2;

      // second column/row
      LSMatrix(1, 2) = LSMatrix(2, 1) = sumXY3;
      LSMatrix(1, 3) = LSMatrix(3, 1) = sumXY2;
      LSMatrix(1, 4) = LSMatrix(4, 1) = sumY3;
      LSMatrix(1, 5) = LSMatrix(5, 1) = sumY2;

      // third column/row
      LSMatrix(2, 3) = LSMatrix(3, 2) = sumX2Y;
      LSMatrix(2, 4) = LSMatrix(4, 2) = sumXY2;
      LSMatrix(2, 5) = LSMatrix(5, 2) = sumXY;

      // fourth column/row
      LSMatrix(3, 4) = LSMatrix(4, 3) = sumXY;
      LSMatrix(3, 5) = LSMatrix(5, 3) = sumX;

      // fifth coulmn/row
      LSMatrix(4, 5) = LSMatrix(5, 4) = sumY;
    }
    else {
      // diagonal
      LSMatrix(0, 0) = sumX4;
      LSMatrix(1, 1) = sumY4;
      LSMatrix(2, 2) = sumZ4;
      LSMatrix(3, 3) = sumX2Y2;
      LSMatrix(4, 4) = sumX2Z2;
      LSMatrix(5, 5) = sumY2Z2;
      LSMatrix(6, 6) = sumX2;
      LSMatrix(7, 7) = sumY2;
      LSMatrix(8, 8) = sumZ2;
      LSMatrix(9, 9) = noInterpolationPoints;

      // first column/row
      LSMatrix(0, 1) = LSMatrix(1, 0) = sumX2Y2;
      LSMatrix(0, 2) = LSMatrix(2, 0) = sumX2Z2;
      LSMatrix(0, 3) = LSMatrix(3, 0) = sumX3Y;
      LSMatrix(0, 4) = LSMatrix(4, 0) = sumX3Z;
      LSMatrix(0, 5) = LSMatrix(5, 0) = sumX2YZ;
      LSMatrix(0, 6) = LSMatrix(6, 0) = sumX3;
      LSMatrix(0, 7) = LSMatrix(7, 0) = sumX2Y;
      LSMatrix(0, 8) = LSMatrix(8, 0) = sumX2Z;
      LSMatrix(0, 9) = LSMatrix(9, 0) = sumX2;

      // second column/row
      LSMatrix(1, 2) = LSMatrix(2, 1) = sumY2Z2;
      LSMatrix(1, 3) = LSMatrix(3, 1) = sumXY3;
      LSMatrix(1, 4) = LSMatrix(4, 1) = sumXY2Z;
      LSMatrix(1, 5) = LSMatrix(5, 1) = sumY3Z;
      LSMatrix(1, 6) = LSMatrix(6, 1) = sumXY2;
      LSMatrix(1, 7) = LSMatrix(7, 1) = sumY3;
      LSMatrix(1, 8) = LSMatrix(8, 1) = sumY2Z;
      LSMatrix(1, 9) = LSMatrix(9, 1) = sumY2;

      // third column/row
      LSMatrix(2, 3) = LSMatrix(3, 2) = sumXYZ2;
      LSMatrix(2, 4) = LSMatrix(4, 2) = sumXZ3;
      LSMatrix(2, 5) = LSMatrix(5, 2) = sumYZ3;
      LSMatrix(2, 6) = LSMatrix(6, 2) = sumXZ2;
      LSMatrix(2, 7) = LSMatrix(7, 2) = sumYZ2;
      LSMatrix(2, 8) = LSMatrix(8, 2) = sumZ3;
      LSMatrix(2, 9) = LSMatrix(9, 2) = sumZ2;

      // fourth column, fourth line
      LSMatrix(3, 4) = LSMatrix(4, 3) = sumX2YZ;
      LSMatrix(3, 5) = LSMatrix(5, 3) = sumXY2Z;
      LSMatrix(3, 6) = LSMatrix(6, 3) = sumX2Y;
      LSMatrix(3, 7) = LSMatrix(7, 3) = sumXY2;
      LSMatrix(3, 8) = LSMatrix(8, 3) = sumXYZ;
      LSMatrix(3, 9) = LSMatrix(9, 3) = sumXY;

      // fifth coulmn, fifth line
      LSMatrix(4, 5) = LSMatrix(5, 4) = sumXYZ2;
      LSMatrix(4, 6) = LSMatrix(6, 4) = sumX2Z;
      LSMatrix(4, 7) = LSMatrix(7, 4) = sumXYZ;
      LSMatrix(4, 8) = LSMatrix(8, 4) = sumXZ2;
      LSMatrix(4, 9) = LSMatrix(9, 4) = sumXZ;

      // sixth column, sixth line
      LSMatrix(5, 6) = LSMatrix(6, 5) = sumXYZ;
      LSMatrix(5, 7) = LSMatrix(7, 5) = sumY2Z;
      LSMatrix(5, 8) = LSMatrix(8, 5) = sumYZ2;
      LSMatrix(5, 9) = LSMatrix(9, 5) = sumYZ;

      // seventh column, seventh line
      LSMatrix(6, 7) = LSMatrix(7, 6) = sumXY;
      LSMatrix(6, 8) = LSMatrix(8, 6) = sumXZ;
      LSMatrix(6, 9) = LSMatrix(9, 6) = sumX;

      // eigth column, eigth line
      LSMatrix(7, 8) = LSMatrix(8, 7) = sumYZ;
      LSMatrix(7, 9) = LSMatrix(9, 7) = sumY;

      // nineth column, nineth line
      LSMatrix(8, 9) = LSMatrix(9, 8) = sumZ;
    }

    // 2.3 scale solution since the condition number gets really bad for large discrepancies in cell size
    const MFloat normalizationFactor = FPOW2(2 * c_level(cellId)) / c_cellLengthAtLevel(0);
    for(MInt i = 0; i < 2 * nDim + pow(2, nDim - 1); i++) {
      for(MInt j = 0; j < 2 * nDim + pow(2, nDim - 1); j++) {
        LSMatrix(i, j) *= normalizationFactor;
      }
    }

    // 2.4 determine the pseudoinverse using SVD
    maia::math::invert(LSMatrix, matInv, 2 * nDim + pow(2, nDim - 1), 2 * nDim + pow(2, nDim - 1));

    ///
    /// 3. Build rhs and solve the LS problem
    ///
    const MInt totalNumberOfVars = interpolateVelocitySlopes ? b + nDim * nDim : b;
    for(MInt k = a; k < totalNumberOfVars; k++) {
      MFloat sumVar = 0;
      if(k < b) {
        sumVar = a_fluidVariable(cellId, k);
      } else {
        sumVar = a_velocitySlope(cellId, std::floor((k - b) / 3), (k - b) % 3);
      }
      MFloat sumVarX = 0;
      MFloat sumVarY = 0;
      MFloat sumVarZ = 0;
      MFloat sumVarXY = 0;
      MFloat sumVarXZ = 0;
      MFloat sumVarYZ = 0;
      MFloat sumVarX2 = 0;
      MFloat sumVarY2 = 0;
      MFloat sumVarZ2 = 0;

      rhs.set(0.0);

      // 3.1 Calculate intermediate variables
      for(MUint nId = 1; nId < noInterpolationPoints; nId++) {
        MFloat x = 0;
        MFloat y = 0;
        MFloat z = 0;
        MFloat var = 0;
        if(nId >= nghbrList.size()) {
          const MUint pId = nId - nghbrList.size();
          x = pseudoCellPos[pId][0] - originX;
          y = pseudoCellPos[pId][1] - originY;
          var = pseudoCellVar[pId][k - a];
          IF_CONSTEXPR(nDim == 3) { z = pseudoCellPos[pId][2] - originZ; }
        } else {
          const MInt nghbrId = nghbrList[nId];
          x = c_coordinate(nghbrId, 0) - originX;
          y = c_coordinate(nghbrId, 1) - originY;
          if(k < b) {
            var = a_fluidVariable(nghbrId, k);
          } else {
            var = a_velocitySlope(nghbrId, std::floor((k - b) / 3), (k - b) % 3);
          }
          IF_CONSTEXPR(nDim == 3) { z = c_coordinate(nghbrId, 2) - originZ; }
        }

        sumVar += var;
        sumVarX += var * x;
        sumVarY += var * y;
        sumVarZ += var * z;
        sumVarXY += var * x * y;
        sumVarXZ += var * x * z;
        sumVarYZ += var * y * z;
        sumVarX2 += var * x * x;
        sumVarY2 += var * y * y;
        sumVarZ2 += var * z * z;
      }

      // 3.2 assign intermediate variables to rhs
      IF_CONSTEXPR(nDim == 2) {
        rhs[0] = sumVarX2;
        rhs[1] = sumVarY2;
        rhs[2] = sumVarXY;
        rhs[3] = sumVarX;
        rhs[4] = sumVarY;
        rhs[5] = sumVar;
      }
      else {
        rhs(0) = sumVarX2;
        rhs(1) = sumVarY2;
        rhs(2) = sumVarZ2;
        rhs(3) = sumVarXY;
        rhs(4) = sumVarXZ;
        rhs(5) = sumVarYZ;
        rhs(6) = sumVarX;
        rhs(7) = sumVarY;
        rhs(8) = sumVarZ;
        rhs(9) = sumVar;
      }

      // 3.3 scale rhs
      for(MInt i = 0; i < 2 * nDim + pow(2, nDim - 1); i++) {
        rhs(i) *= normalizationFactor;
      }

      MFloatTensor coeff(2 * nDim + pow(2, nDim - 1));
      coeff.set(F0);

      // 3.4 determine coefficients of the aim function polynomial
      for(MInt i = 0; i < 2 * nDim + pow(2, nDim - 1); i++) {
        for(MInt j = 0; j < 2 * nDim + pow(2, nDim - 1); j++) {
          coeff(i) += matInv(i, j) * rhs(j);
        }
      }

      ///
      /// 4. determine solution
      ///
      MFloat positionX = position[0] - originX;
      MFloat positionY = position[1] - originY;
      MFloat positionZ = 0;
      IF_CONSTEXPR(nDim == 3) { positionZ = position[2] - originZ; }
      MFloat result = coeff((MInt)(2 * nDim + pow(2, nDim - 1) - 1));
      IF_CONSTEXPR(nDim == 2) {
        result += coeff(0) * positionX * positionX;
        result += coeff(1) * positionY * positionY;
        result += coeff(2) * positionX * positionY;
        result += coeff(3) * positionX;
        result += coeff(4) * positionY;
      }
      else {
        result += coeff(0) * positionX * positionX;
        result += coeff(1) * positionY * positionY;
        result += coeff(2) * positionZ * positionZ;
        result += coeff(3) * positionX * positionY;
        result += coeff(4) * positionX * positionZ;
        result += coeff(5) * positionY * positionZ;
        result += coeff(6) * positionX;
        result += coeff(7) * positionY;
        result += coeff(8) * positionZ;
      }

      interpolatedVar[k - a] = result;
    }
  }

  // TRACE_OUT();
}

template void LPT<2>::interpolateVariablesLS<0, 2, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 2, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 3, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 3, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 5, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 5, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 6, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<2>::interpolateVariablesLS<0, 6, false>(const MInt, const MFloat* const, MFloat* const);

template void LPT<3>::interpolateVariablesLS<0, 3, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 3, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 4, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 4, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 6, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 6, false>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 7, true>(const MInt, const MFloat* const, MFloat* const);
template void LPT<3>::interpolateVariablesLS<0, 7, false>(const MInt, const MFloat* const, MFloat* const);


template class LPT<3>;
