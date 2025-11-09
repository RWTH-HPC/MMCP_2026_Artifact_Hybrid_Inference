// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "rigidbodies.h"
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "IO/parallelio.h"
#include "UTIL/maiamath.h"
#include "typetraits.h"


#include <algorithm>
#include <functional>
#include <iomanip>

/**
 * \brief C'tor for the RigidBodies solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param solverId_ Unique solver id
 * \param gridProxy_ Grid proxy representing the solvers grid information
 * \param geometry_ Holding all relevant geometry information for this solver
 * \param comm_ MPI communicator between all ranks responsible for this solver
 */

template <MInt nDim>
RigidBodies<nDim>::RigidBodies(const MInt solverId_, GridProxy& gridProxy_, Geom& geometry_, MPI_Comm comm_)
  : maia::CartesianSolver<nDim, RigidBodies<nDim>>(solverId_, gridProxy_, comm_), m_geometry(&geometry_) {
  TRACE();

  initTimers();
  readProperties();

#ifdef RB_DEBUG
  m_debugFileStream.open("d" + std::to_string(domainId()) + ".txt");
#endif

  initBodyData();

  if(!domainId() && !m_restart) {
    initBodyLog();
  }

  // create body communication mappings
  initGridProperties();
  updateInfoDiameter(true);
  updateBodyDomainConnections(true);
  exchangeEdgeBodyData();
// for indirect neighbors (large bodies, to be improved)
// exchangeNeighborConnectionInfo();
#ifdef RB_DEBUG
  printBodyDomainConnections();
#endif
  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

/**
 * \brief D'tor for the RigidBodies solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
RigidBodies<nDim>::~RigidBodies() {
  TRACE();
  // Close stream
  if(!domainId()) {
    m_anglog.close();
  }

#ifdef RB_DEBUG
  m_debugFileStream.close();
#endif
  // Deallocate memory

  // stop timer
  RECORD_TIMER_STOP(m_timers[Timers::Class]);

  averageTimer();
  printScalingVariables();
}

template <MInt nDim>
void RigidBodies<nDim>::initTimers() {
  TRACE();

  // Create timer group and coupler timer, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timers[Timers::TimerGroup], "RigidBodies");
  NEW_TIMER_NOCREATE(m_timers[Timers::Class], "Total object lifetime", m_timers[Timers::TimerGroup]);
  RECORD_TIMER_START(m_timers[Timers::Class]);

  // Create and start constructor timer
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Constructor], "Constructor", m_timers[Timers::Class]);
  RECORD_TIMER_START(m_timers[Timers::Constructor]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SolutionStep], "Solution Step", m_timers[Timers::Class]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Predict], "Predict", m_timers[Timers::SolutionStep]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Correct], "Correct", m_timers[Timers::SolutionStep]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Communication], "Communication", m_timers[Timers::SolutionStep]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::UpdateCommunication], "Update Communication", m_timers[Timers::SolutionStep]);
}

template <MInt nDim>
void RigidBodies<nDim>::averageTimer() {
  TRACE();
  if(!grid().isActive()) return;

  m_timerType = "max";
  m_timerType = Context::getSolverProperty<MString>("timerType", m_solverId, AT_, &m_timerType);

  // 0) map timer ids for safety

  const MInt noTimers = m_timers.size();

  // 1) fill buffer with local timer values
  std::vector<MFloat> timerValues_;
  timerValues_.reserve(noTimers);
  for(MInt i = 0; i < noTimers; i++) {
    timerValues_.emplace_back(RETURN_TIMER_TIME(m_timers[i]));
  }

  // 2) collect values from all ranks
  if(m_timerType == "average") {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_SUM, mpiComm(),
                  AT_, "MPI_IN_PLACE", "timerValues_");
  } else {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_MAX, mpiComm(),
                  AT_, "MPI_IN_PLACE", "timerValues_");
  }

  // 3) perform averaging on timer and4) set new timer values
  if(m_timerType == "average") {
    const MInt noDomains_ = noDomains();
    for(MInt i = 0; i < noTimers; i++) {
      const MFloat meanValue = timerValues_[i] / noDomains_;
      SET_RECORD(m_timers[i], meanValue);
    }
  } else {
    for(MInt i = 0; i < noTimers; i++) {
      SET_RECORD(m_timers[i], timerValues_[i]);
    }
  }
}

/**
 * \brief Reading all necessary properties for the RigidBodies solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::readProperties() {
  TRACE();

  m_outputDir = "./";
  m_outputDir = Context::getSolverProperty<MString>("outputDir", m_solverId, AT_, &m_outputDir);

  m_printKineticEnergy = false;
  Context::getSolverProperty<MBool>("printKineticEnergy", m_solverId, AT_, &m_printKineticEnergy);

  m_intervalBodyLog = 1;
  m_intervalBodyLog = Context::getSolverProperty<MInt>("intervalBodyLog", m_solverId, AT_, &m_intervalBodyLog);

  // restart related Properties
  m_restart = false;
  m_restart = Context::getSolverProperty<MBool>("restartFile", m_solverId, AT_, &m_restart);

  m_restartTimeStep = 0;
  m_restartTimeStep = Context::getSolverProperty<MInt>("restartTimeStep", m_solverId, AT_, &m_restartTimeStep);

  m_restartDir = Context::getSolverProperty<MString>("restartDir", m_solverId, AT_, &m_outputDir);

  m_bodyCenterInitMethod = "POINT";
  m_bodyCenterInitMethod =
      Context::getSolverProperty<MString>("bodyCenterInitMethod", m_solverId, AT_, &m_bodyCenterInitMethod);

  if(!m_restart) {
    if(m_bodyCenterInitMethod == "BOX_SEED") {
      boxSeedInitialize();
    } else if(m_bodyCenterInitMethod == "POINT") {
      m_noEmbeddedBodies = Context::getSolverProperty<MInt>("noEmbeddedBodies", m_solverId, AT_, &m_noEmbeddedBodies);

      const MInt num_initialBodyCenters = Context::propertyLength("initialBodyCenters", m_solverId);
      for(MInt i = 0; i < num_initialBodyCenters; i++) {
        m_initialBodyCenters.emplace_back(Context::getSolverProperty<MFloat>("initialBodyCenters", m_solverId, AT_, i));
      }
    }
  }

  m_maxNoBodies = m_noEmbeddedBodies;
  m_maxNoBodies = Context::getSolverProperty<MInt>("maxNoBodies", m_solverId, AT_, &m_maxNoBodies);

  m_forcedMotion = Context::getSolverProperty<MBool>("forcedMotion", m_solverId, AT_, &m_forcedMotion);

  m_bodyType = 1;
  m_bodyType = Context::getSolverProperty<MInt>("bodyTypeMb", m_solverId, AT_, &m_bodyType);

  const MInt num_initBodyRadius = Context::propertyLength("bodyRadius", m_solverId);
  for(MInt i = 0; i < num_initBodyRadius; i++) {
    m_initBodyRadius.emplace_back(Context::getSolverProperty<MFloat>("bodyRadius", m_solverId, AT_, i));
  }

  const MInt num_initBodyRadii = Context::propertyLength("bodyRadii", m_solverId);
  for(MInt i = 0; i < num_initBodyRadii; i++) {
    m_initBodyRadii.emplace_back(Context::getSolverProperty<MFloat>("bodyRadii", m_solverId, AT_, i));
  }

  m_bodyHeight = 1.0;
  m_bodyHeight = Context::getSolverProperty<MFloat>("bodyHeight", m_solverId, AT_, &m_bodyHeight);

  m_uRef = 1.0;
  m_uRef = Context::getSolverProperty<MFloat>("uRef", m_solverId, AT_, &m_uRef);

  m_logBodyData = false;
  m_logBodyData = Context::getSolverProperty<MBool>("logBodyData", m_solverId, AT_, &m_logBodyData);

  // Array-like properties
  const MInt num_initBodyDensityRatios = Context::propertyLength("bodyDensityRatio", m_solverId);
  for(MInt i = 0; i < num_initBodyDensityRatios; i++) {
    m_initBodyDensityRatios.emplace_back(Context::getSolverProperty<MFloat>("bodyDensityRatio", m_solverId, AT_, i));
  }

  if(Context::propertyExists("gravity", m_solverId)) {
    for(MInt n = 0; n < nDim; n++) {
      gravity[n] = Context::getSolverProperty<MFloat>("gravity", m_solverId, AT_, n);
      gravity[n] /= POW2(m_uRef);
    }
  }

  // TODO labels:DOC,toenhance Change to per-body value
  m_uniformBodyTemperature = F1;
  m_uniformBodyTemperature =
      Context::getSolverProperty<MFloat>("uniformBodyTemperature", m_solverId, AT_, &m_uniformBodyTemperature);

  // allow translation
  m_translation = true;
  m_translation = Context::getSolverProperty<MBool>("translation", m_solverId, AT_, &m_translation);

  // allow rotation
  m_rotation = true;
  m_rotation = Context::getSolverProperty<MBool>("integrateRotation", m_solverId, AT_, &m_rotation);

  // allow rotation only around z-axis
  m_rotXYaxis = true;
  m_rotXYaxis = Context::getSolverProperty<MBool>("rotXYaxis", m_solverId, AT_, &m_rotXYaxis);
}

/**
 * \brief Allocating and initializing body data necessary for the solver run
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::initBodyData() {
  TRACE();

  if(m_restart) {
    loadBodiesSizeAndPosition();
  }

  findLocalBodies();

  m_bodies.reset(m_maxNoBodies);

  // Create Mapping
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt globDomain = neighborDomain(i);
    m_remoteDomainLocalBodies[globDomain] = std::vector<MInt>();
    m_homeDomainRemoteBodies[globDomain] = std::vector<MInt>();
  }

  m_bodies.append(m_noLocalBodies);

  if(m_noLocalBodies == 0 && !m_restart) return;

  // Init bodyRadius and bodyRadii via bodyRadius
  if(!m_initBodyRadius.empty()) {
    // If only one radius is given, use for all bodies
    if(m_initBodyRadius.size() == 1) {
      const MFloat defaultRadius = m_initBodyRadius[0];
      std::fill_n(&a_bodyRadius(0), m_noLocalBodies, defaultRadius);
      std::fill_n(&a_bodyRadii(0, 0), m_noLocalBodies * nDim, defaultRadius);
    } else {
      for(MUint bodyId = 0; bodyId < m_initBodyRadius.size(); bodyId++) {
        a_bodyRadius(bodyId) = m_initBodyRadius[getGlobalBodyId(bodyId)];
        std::fill_n(&a_bodyRadii(bodyId, 0), nDim, m_initBodyRadius[getGlobalBodyId(bodyId)]);
      }
    }
  }

  // Init bodyRadii via bodyRadii
  if(!m_initBodyRadii.empty()) {
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      // If only one set of radii is given, use for all bodies
      if(m_initBodyRadii.size() == nDim) {
        for(MInt n = 0; n < nDim; n++) {
          a_bodyRadii(bodyId, n) = m_initBodyRadii[n];
        }
      } else {
        for(MInt n = 0; n < nDim; n++) {
          a_bodyRadii(bodyId, n) = m_initBodyRadii[getGlobalBodyId(bodyId) * nDim + n];
        }
      }
    }
  }

  for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
    for(MInt n = 0; n < nRot; n++) {
      a_angularVelocityT1(bodyId, n) = 0.0;
      a_angularAccelerationT1(bodyId, n) = 0.0;
    }
  }

  if(!m_forcedMotion) {
    if(m_initBodyDensityRatios.size() > 1) {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        a_bodyDensityRatio(bodyId) = m_initBodyDensityRatios[getGlobalBodyId(bodyId)];
      }
    } else if(m_initBodyDensityRatios.size() == 1) {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        a_bodyDensityRatio(bodyId) = m_initBodyDensityRatios[0];
      }
    } else {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        a_bodyDensityRatio(bodyId) = 1.0;
      }
    }

    // Init rotational data
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      for(MInt n = 0; n < nRot; n++) {
        a_bodyInertia(bodyId, n) = 1.0;
        a_angularVelocityT1B2(bodyId, n) = 0.0;
        a_angularVelocityBodyT1(bodyId, n) = 0.0;
        a_angularVelocityBodyT1B2(bodyId, n) = 0.0;
        a_angularAccelerationBody(bodyId, n) = 0.0;
        a_torqueT1(bodyId, n) = 0.0;
      }
    }

    // Init bodyQuaternion and moment of inertia
    IF_CONSTEXPR(nDim == 3) {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        for(MInt n = 0; n < nQuat - 1; n++) {
          a_bodyQuaternionT1B2(bodyId, n) = 0;
          a_bodyQuaternionT1(bodyId, n) = 0;
        }
        a_bodyQuaternionT1B2(bodyId, nQuat - 1) = 1.0;
        a_bodyQuaternionT1(bodyId, nQuat - 1) = 1.0;
        a_bodyInertia(bodyId, 0) =
            (POW2(a_bodyRadii(bodyId, 1)) + POW2(a_bodyRadii(bodyId, 2))) / 5.0 * a_bodyMass(bodyId);
        a_bodyInertia(bodyId, 1) =
            (POW2(a_bodyRadii(bodyId, 0)) + POW2(a_bodyRadii(bodyId, 2))) / 5.0 * a_bodyMass(bodyId);
        a_bodyInertia(bodyId, 2) =
            (POW2(a_bodyRadii(bodyId, 0)) + POW2(a_bodyRadii(bodyId, 1))) / 5.0 * a_bodyMass(bodyId);
      }
    }
  }

  if(m_restart) {
    loadBodyRestartFile();
    // Init remaining values
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      a_bodyHeatFlux(bodyId) = 0.0;
      const MInt tmpGlobalId = getGlobalBodyId(bodyId);

      a_bodyTemperature(bodyId) = m_bResFile.bodyTemperature(tmpGlobalId);
      for(MInt n = 0; n < nDim; n++) {
        a_bodyCenter(bodyId, n) = m_bResFile.bodyCenter(tmpGlobalId, n);
        a_bodyVelocity(bodyId, n) = m_bResFile.bodyVelocity(tmpGlobalId, n);
        a_bodyAcceleration(bodyId, n) = m_bResFile.bodyAcceleration(tmpGlobalId, n);
      }
      for(MInt r = 0; r < nRot; r++) {
        a_angularVelocityBodyT1B2(bodyId, r) = m_bResFile.angularVelocityBodyT1B2(tmpGlobalId, r);
        a_angularAccelerationBody(bodyId, r) = m_bResFile.angularAccelerationBody(tmpGlobalId, r);
      }
      for(MInt q = 0; q < nQuat; q++) {
        a_bodyQuaternionT1B2(bodyId, q) = m_bResFile.bodyQuaternionT1B2(tmpGlobalId, q);
        a_bodyQuaternionT1(bodyId, q) = m_bResFile.bodyQuaternionT1(tmpGlobalId, q);
      }
    }
  } else {
    // Set init body values
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      for(MInt n = 0; n < nDim; n++) {
        a_bodyCenter(bodyId, n) = m_initialBodyCenters[getGlobalBodyId(bodyId) * nDim + n];
        a_bodyVelocity(bodyId, n) = 0.0;
        a_bodyAcceleration(bodyId, n) = 0.0;
        a_bodyForce(bodyId, n) = 0.0;
      }
      a_bodyTemperature(bodyId) = 1.0;
      a_bodyHeatFlux(bodyId) = 0.0;
    }
  }

  // Transfer initialValues to *Dt1 fields.
  if(!m_forcedMotion) {
    advanceBodies();
  }
}

/**
 * \brief Initalize log-file
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::initBodyLog() {
  TRACE();

  // delete output file if it already exists
  MBool body_exists = fileExists("bodies.log");
  MBool angBody_exists = fileExists("angvel.log");

  if(body_exists) {
    std::remove("bodies.log");
    m_log << "replaced body log";
  }

  if(angBody_exists) {
    std::remove("angvel.log");
    m_log << "replaced angular body log";
  }

  m_log.open("bodies.log", std::ios::app);
  m_anglog.open("angvel.log", std::ios::app);

  if(!m_logBodyData) {
    if(!domainId()) {
      std::stringstream logEntry;
      std::stringstream logEntryAng;

      logEntry << ">>> enable body log by using: 'logBodyData = true' "
               << "\n";
      logEntryAng << ">>> enable body log by using: 'logBodyData = true' "
                  << "\n";

      m_log << logEntry.str();
      m_log.close();

      m_anglog << logEntryAng.str();
      m_anglog.close();
    }
    return;
  }

  MString delimiter = " ";
  m_logVars.insert(m_logVars.end(), {"t", "body", "x", "y", "z", "u", "v", "w", "c_x", "c_y", "c_z"});

  m_logAngVars.insert(m_logAngVars.end(),
                      {"t", "body", "ang_vel_x", "ang_vel_y", "ang_vel_z", "tourque_x", "tourque_y", "tourque_z"});

  IF_CONSTEXPR(nDim == 2) {
    std::vector<MString> removeList{"z", "w", "c_z"};
    std::vector<MString> removeListAng{"ang_vel_z", "tourque_z"};
    for(const MString& r : removeList) {
      m_logVars.erase(std::remove(m_logVars.begin(), m_logVars.end(), r), m_logVars.end());
    }

    for(const MString& r : removeListAng) {
      m_logAngVars.erase(std::remove(m_logAngVars.begin(), m_logAngVars.end(), r), m_logAngVars.end());
    }
  }

  for(const MString& var : m_logVars) {
    m_log << var << delimiter;
  }

  for(const MString& var : m_logAngVars) {
    m_anglog << var << delimiter;
  }

  m_log << "\n";
  m_anglog << "\n";

  m_log.close();
  m_anglog.close();
}

template <MInt nDim>
void RigidBodies<nDim>::printScalingVariables() {
  TRACE();

  std::map<MString, MInt> scalingVars;
  scalingVars["localBodies"] = m_noLocalBodies;
  scalingVars["connectedBodies"] = noConnectingBodies();
  scalingVars["neighborDomains"] = noNeighborDomains();

  std::vector<MInt> recvBuffer;
  recvBuffer.resize(noDomains());
  for(const auto& [name, value] : scalingVars) {
    MPI_Gather(&value, 1, maia::type_traits<MInt>::mpiType(), recvBuffer.data(), 1, maia::type_traits<MInt>::mpiType(),
               0, mpiComm(), AT_, "send", "recv");

    if(!domainId()) {
      // average
      MFloat avgValue = 0;
      avgValue = std::accumulate(recvBuffer.begin(), recvBuffer.end(), avgValue);
      avgValue /= noDomains();
      std::cout << "Average " + name + " per Domain are: " << avgValue << std::endl;

      // maximum
      const MInt maxValue = *std::max_element(recvBuffer.begin(), recvBuffer.end());
      std::cout << "Max " + name + " per Domain are: " << maxValue << std::endl;

      // minimum
      const MInt minValue = *std::min_element(recvBuffer.begin(), recvBuffer.end());
      std::cout << "Min " + name + " per Domain are: " << minValue << std::endl;
    }
  }
}


template <MInt nDim>
void RigidBodies<nDim>::totalKineticEnergy() {
  TRACE();

  if(!(globalTimeStep % m_solutionInterval == 0) || !m_printKineticEnergy) return;
  std::vector<MFloat> recvBuffer;
  recvBuffer.resize(noDomains());

  MFloat summedKineticEnergy = 0.0;
  for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
    summedKineticEnergy += F1B2 * a_bodyDensityRatio(bodyId) * a_volume(bodyId) * POW2(a_bodyVelocityMag(bodyId));
  }

  MPI_Gather(&summedKineticEnergy, 1, maia::type_traits<MFloat>::mpiType(), recvBuffer.data(), 1,
             maia::type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_, "send", "recv");

  MFloat totalKineticEnergy = 0.0;
  if(!domainId()) {
    totalKineticEnergy = std::accumulate(recvBuffer.begin(), recvBuffer.end(), totalKineticEnergy);
    std::cout << "Total Kinetic Energy of Bodies is: " << totalKineticEnergy << std::endl;
  }
}


template <MInt nDim>
void RigidBodies<nDim>::boxSeedInitialize() {
  TRACE();

  std::array<MFloat, nDim> bMin{};
  std::array<MFloat, nDim> bMax{};

  for(MInt n = 0; n < nDim; n++) {
    bMin[n] = Context::getSolverProperty<MFloat>("seedBoxMin", m_solverId, AT_, n);
    bMax[n] = Context::getSolverProperty<MFloat>("seedBoxMax", m_solverId, AT_, n);
  }

  MInt bodiesPerEdge = Context::getSolverProperty<MInt>("bodiesPerEdge", m_solverId, AT_, &bodiesPerEdge);
  m_noEmbeddedBodies = pow(bodiesPerEdge, nDim);

  std::array<std::vector<MFloat>, nDim> bodyGrid;
  for(MInt n = 0; n < nDim; n++) {
    bodyGrid[n] = maia::math::linSpace(bMin[n], bMax[n], bodiesPerEdge);
  }

  m_initialBodyCenters.resize(MLong(m_noEmbeddedBodies * nDim), F0);
  for(MInt i = 0; i < bodiesPerEdge; i++) {
    MFloat x = bodyGrid[0][i];
    for(MInt j = 0; j < bodiesPerEdge; j++) {
      MFloat y = bodyGrid[1][j];
      for(MInt k = 0; k < bodiesPerEdge; k++) {
        MFloat z = bodyGrid[2][k];
        m_initialBodyCenters[i * nDim * pow(bodiesPerEdge, 2) + j * nDim * bodiesPerEdge + k * nDim] = x;
        m_initialBodyCenters[i * nDim * pow(bodiesPerEdge, 2) + j * nDim * bodiesPerEdge + k * nDim + 1] = y;
        m_initialBodyCenters[i * nDim * pow(bodiesPerEdge, 2) + j * nDim * bodiesPerEdge + k * nDim + 2] = z;
      }
    }
  }
}

/**
 * \brief Resetting the 2-Step predictor-corrector cycle
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::preTimeStep() {
  m_newTimeStep = true;
}

/**
 * \brief
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \returns True if the predictor-corrector cycle has finished
 */
template <MInt nDim>
MBool RigidBodies<nDim>::solutionStep() {
  if(m_newTimeStep) {
#ifdef RB_DEBUG
    m_debugFileStream
        << "------------TIMESTEP: " << globalTimeStep
        << " -----------------------------------------------------------------------------------------------"
        << std::endl;
#endif

    RECORD_TIMER_START(m_timers[Timers::SolutionStep]);

    RECORD_TIMER_START(m_timers[Timers::Predict]);
    computeBodies();
    RECORD_TIMER_STOP(m_timers[Timers::Predict]);

    RECORD_TIMER_START(m_timers[Timers::Communication]);
    exchangeKinematicData();
    RECORD_TIMER_STOP(m_timers[Timers::Communication]);

    m_newTimeStep = false;

    RECORD_TIMER_STOP(m_timers[Timers::SolutionStep]);

    return false;
  } else {
    RECORD_TIMER_START(m_timers[Timers::SolutionStep]);

    RECORD_TIMER_START(m_timers[Timers::Communication]);
    exchangeFsiData();
    RECORD_TIMER_STOP(m_timers[Timers::Communication]);

    RECORD_TIMER_START(m_timers[Timers::Correct]);
    correctBodies();
    RECORD_TIMER_STOP(m_timers[Timers::Correct]);

    RECORD_TIMER_START(m_timers[Timers::Communication]);
    exchangeKinematicData();
    RECORD_TIMER_STOP(m_timers[Timers::Communication]);

    RECORD_TIMER_START(m_timers[Timers::UpdateCommunication]);
    updateConnections();
    RECORD_TIMER_STOP(m_timers[Timers::UpdateCommunication]);

    advanceBodies();
    RECORD_TIMER_STOP(m_timers[Timers::SolutionStep]);

    return true;
  }
}

/**
 * \brief Calculates the new state for all bodies
 *
 * This is done via prediction or by using an analytic motion equation (forcedMotion)
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::computeBodies() {
  TRACE();

  // dimensionless RB time
  const MFloat time = globalTimeStep * m_timestep;

  if(m_forcedMotion) {
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      computeBodyPropertiesForced(1, &m_bodies.bodyCenter(bodyId, 0), bodyId, time);
      computeBodyPropertiesForced(2, &m_bodies.bodyVelocity(bodyId, 0), bodyId, time);
      computeBodyPropertiesForced(3, &m_bodies.bodyAcceleration(bodyId, 0), bodyId, time);
    }
  } else {
    predictorStep();
  }
}

/**
 * \brief exchange of predictor & corrector step output
 *
 * local(edge) bodies' kinematic data is send to remote domains +
 * remote bodies' kinematic data is received from home domains
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeKinematicData() {
  TRACE();

  MInt noRemoteDomains = m_remoteDomainLocalBodies.size();
  MInt noHomeDomains = m_homeDomainRemoteBodies.size();

  std::vector<MPI_Request> mpi_send_req(noRemoteDomains, MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noHomeDomains, MPI_REQUEST_NULL);

  MInt bufSizePerBody = nDim * 3 + nRot * 2 + nQuat;
  // allocate receive buffer per homeDomain
  std::vector<std::vector<MFloat>> rcvBufferAll(noHomeDomains);
  MInt idx = 0;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    rcvBufferAll[idx].resize(bodies.size() * bufSizePerBody);
    idx++;
  }

  // receive predicted body state
  idx = 0;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    MPI_Irecv(rcvBufferAll[idx].data(), rcvBufferAll[idx].size(), maia::type_traits<MFloat>::mpiType(), homeDomain, 2,
              mpiComm(), &mpi_recv_req[idx], AT_, "recv predict");
    idx++;
  }

  idx = 0;
  MInt paraCount = 0;
  MInt bodyCount = 0;
  std::vector<std::vector<MFloat>> sendBufferAll(noRemoteDomains);
  // pack send buffer and send predicted body states to all remote domains
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    // pack
    auto& sendBuffer = sendBufferAll[idx];
    sendBuffer.resize(bufSizePerBody * bodies.size());
    bodyCount = 0;
    for(const auto& body : bodies) {
      paraCount = 0;
      for(MInt n = 0; n < nDim; n++) {
        sendBuffer[bodyCount + n] = a_bodyCenter(body, n);
        // apply shift if periodic BCs are active
        if(m_periodicBC) {
          sendBuffer[bodyCount + n] += m_periodicShift[remoteDomain][body][n];
        }
        sendBuffer[bodyCount + nDim + n] = a_bodyVelocity(body, n);
        sendBuffer[bodyCount + 2 * nDim + n] = a_bodyAcceleration(body, n);
      }
      paraCount += 3 * nDim;
      for(MInt r = 0; r < nRot; r++) {
        sendBuffer[bodyCount + paraCount + r] = a_angularVelocity(body, r);
        sendBuffer[bodyCount + paraCount + nRot + r] = a_angularVelocityBody(body, r);
      }
      paraCount += 2 * nRot;
      for(MInt q = 0; q < nQuat; q++) {
        sendBuffer[bodyCount + paraCount + q] = a_bodyQuaternionT1(body, q);
      }
      bodyCount += bufSizePerBody;
    }
    // send
    MPI_Isend(sendBuffer.data(), sendBuffer.size(), maia::type_traits<MFloat>::mpiType(), remoteDomain, 2, mpiComm(),
              &mpi_send_req[idx], AT_, "send predict");
    idx++;
  }

  // Wait for MPI exchange to complete
  for(MInt i = 0; i < noRemoteDomains; i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
  }

  for(MInt i = 0; i < noHomeDomains; i++) {
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }

  // unpack buffer -> same as packing
  idx = 0;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    bodyCount = 0;
    for(const auto& body : bodies) {
      paraCount = 0;
      for(MInt n = 0; n < nDim; n++) {
        a_bodyCenter(body, n) = rcvBufferAll[idx][bodyCount + n];
        a_bodyVelocity(body, n) = rcvBufferAll[idx][bodyCount + nDim + n];
        a_bodyAcceleration(body, n) = rcvBufferAll[idx][bodyCount + nDim * 2 + n];
      }
      paraCount += 3 * nDim;
      for(MInt r = 0; r < nRot; r++) {
        a_angularVelocity(body, r) = rcvBufferAll[idx][bodyCount + paraCount + r];
        a_angularVelocityBody(body, r) = rcvBufferAll[idx][bodyCount + paraCount + nRot + r];
      }
      paraCount += 2 * nRot;
      for(MInt q = 0; q < nQuat; q++) {
        a_bodyQuaternionT1(body, q) = rcvBufferAll[idx][bodyCount + paraCount + q];
      }
      bodyCount += bufSizePerBody;
    }
    idx++;
  }

  // Exchange kinematic Data for self-mapped Bodies
  auto aDBCopy = m_associatedDummyBodies;
  for(const auto& [body, dummyId] : aDBCopy) {
    // if a remote body has a associated dummy body, but is not intersecting with a halo cell anymore -> delete dummy
    // body
    const MInt intersectingCellId =
        grid().raw().intersectingWithHaloCells(&a_bodyCenter(body, 0), m_infoDiameter / 2.0, 0, true);
    if(intersectingCellId == -1) {
      deleteDummyBody(body);
      continue;
    }
    for(MInt n = 0; n < nDim; n++) {
      // periodic shift has to be updated?
      m_periodicShift[domainId()][body][n] = calculatePeriodicShift(intersectingCellId, n);
      a_bodyCenter(dummyId, n) = a_bodyCenter(body, n) + m_periodicShift[domainId()][body][n];
      a_bodyVelocity(dummyId, n) = a_bodyVelocity(body, n);
      a_bodyAcceleration(dummyId, n) = a_bodyAcceleration(body, n);
    }
    for(MInt r = 0; r < nRot; r++) {
      a_angularVelocity(dummyId, r) = a_angularVelocity(body, r);
      a_angularVelocityBody(dummyId, r) = a_angularVelocityBody(body, r);
    }
    for(MInt q = 0; q < nQuat; q++) {
      a_bodyQuaternionT1(dummyId, q) = a_bodyQuaternionT1(body, q);
    }
  }
}


/**
 * \brief exchange of summed forces and torques
 *
 * local(edge) bodies' parts of body forces are received from remote domains +
 * remote bodies' parts of body forces is send to their home domains
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeFsiData() {
  TRACE();

  // exchange local partial sum of the remote body force to the bodies home domain and add it
  std::vector<MPI_Request> mpi_send_req(noHomeDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noRemoteDomains(), MPI_REQUEST_NULL);

  MInt noRemoteDomains = m_remoteDomainLocalBodies.size();
  MInt noHomeDomains = m_homeDomainRemoteBodies.size();

  MInt bufSizePerBody = nDim * 2;
  // allocate receive buffer per remote domain and create tmp sum buffer for each body
  std::vector<std::vector<MFloat>> rcvBufferAll(noRemoteDomains);
  MInt idx = 0;
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    rcvBufferAll[idx].resize(bodies.size() * bufSizePerBody);
    idx++;
  }

  // receive predicted body state
  idx = 0;
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    MPI_Irecv(rcvBufferAll[idx].data(), rcvBufferAll[idx].size(), maia::type_traits<MFloat>::mpiType(), remoteDomain, 2,
              mpiComm(), &mpi_recv_req[idx], AT_, "recv predict");
    idx++;
  }

  idx = 0;
  MInt bodyCount = 0;
  std::vector<std::vector<MFloat>> sendBufferAll(noHomeDomains);
  // pack send buffer and send predicted body states to home domains
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    // pack
    auto& sendBuffer = sendBufferAll[idx];
    sendBuffer.resize(bufSizePerBody * bodies.size());
    bodyCount = 0;
    for(const auto& body : bodies) {
      for(MInt n = 0; n < nDim; n++) {
        // when sending the FSI-Data to the homeDomain, we not only need the forces of the remoteBody, but also the
        // forces of its associated dummy body, in case the remoteBody is self-mapped
        if(hasAssociatedDummyBody(body)) {
          a_bodyForce(body, n) += a_bodyForce(m_associatedDummyBodies[body], n);
          a_bodyTorque(body, n) += a_bodyTorque(m_associatedDummyBodies[body], n);
        }
        sendBuffer[bodyCount + n] = a_bodyForce(body)[n];
        sendBuffer[bodyCount + nDim + n] = a_bodyTorque(body)[n];
      }
      bodyCount += bufSizePerBody;
    }
    // send
    MPI_Isend(sendBuffer.data(), sendBuffer.size(), maia::type_traits<MFloat>::mpiType(), homeDomain, 2, mpiComm(),
              &mpi_send_req[idx], AT_, "send predict");
    idx++;
  }

  // Wait for MPI exchange to complete
  for(MInt i = 0; i < noHomeDomains; i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
  }

  for(MInt i = 0; i < noRemoteDomains; i++) {
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }

  // add up required parts of remotely calculated surface force + torque
  idx = 0;
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    bodyCount = 0;
    for(const auto& body : bodies) {
      for(MInt n = 0; n < nDim; n++) {
        a_bodyForce(body)[n] += rcvBufferAll[idx][bodyCount + n];
        a_bodyTorque(body)[n] += rcvBufferAll[idx][bodyCount + nDim + n];
      }
      bodyCount += bufSizePerBody;
    }
    idx++;
  }

  // Exchange FSI data for self-mapped Bodies
  for(const auto& [body, dummyId] : m_associatedDummyBodies) {
    for(MInt n = 0; n < nDim; n++) {
      a_bodyForce(body, n) += a_bodyForce(dummyId, n);
      a_bodyTorque(body, n) += a_bodyTorque(dummyId, n);
    }
  }
}

template <MInt nDim>
void RigidBodies<nDim>::updateConnections() {
  TRACE();
  findTransferBodies();
  exchangeTransferBodies();
  checkDummyBodiesForSwap();
  updateInfoDiameter();
  updateBodyDomainConnections(false);
  exchangeEdgeBodyData();
// exchangeNeighborConnectionInfo();
#ifdef RB_DEBUG
  printBodyDomainConnections();
#endif
  totalKineticEnergy();
}

/**
 * \brief Correcting the state of all bodies using all available external forces/fluxes
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::correctBodies() {
  TRACE();

  if(!m_forcedMotion) {
    collideBodies();
    correctorStep();
  }

  logBodyData();
}

/**
 * \brief Writing a simple log file containing the current state for all bodies
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::logBodyData() {
  TRACE();
  if(!m_logBodyData || !(globalTimeStep % m_intervalBodyLog == 0)) return;
  // mpi_gatherv body data to write file with domain 0
  // bodyId + pos + vel + acc + ang_vel + torque
  const MInt bufSizePerBody = 1 + nDim * 5;
  std::vector<MFloat> sendBuffer(MLong(bufSizePerBody * m_noLocalBodies));
  for(MInt i = 0; i < m_noLocalBodies; i++) {
    MInt bodyId = i;
    sendBuffer[MLong(i * bufSizePerBody)] = (MFloat)bodyId;
    for(MInt n = 0; n < nDim; n++) {
      sendBuffer[i * bufSizePerBody + 1 + n] = a_bodyCenter(bodyId, n);
      sendBuffer[i * bufSizePerBody + 1 + 1 * nDim + n] = a_bodyVelocity(bodyId, n);
      sendBuffer[i * bufSizePerBody + 1 + 2 * nDim + n] = a_bodyAcceleration(bodyId, n);
      sendBuffer[i * bufSizePerBody + 1 + 3 * nDim + n] = a_angularVelocity(bodyId, n);
      sendBuffer[i * bufSizePerBody + 1 + 4 * nDim + n] = a_bodyTorque(bodyId, n);
    }
  }

  std::vector<MInt> noToRecv(noDomains());
  const MUint noToSend = sendBuffer.size();
  exchangeBufferLengthsAllToRoot(noToRecv, noToSend);

  MInt count = 0;
  std::vector<MInt> displacements(noDomains());
  for(MInt n = 0; n < noDomains(); n++) {
    // noToRecv[n] *= bufSizePerBody;
    displacements[n] = count;
    count += noToRecv[n];
  }

  // create receive buffer
  std::vector<MFloat> rcvBuffer(MLong(m_noEmbeddedBodies * bufSizePerBody));
  MPI_Gatherv(sendBuffer.data(), sendBuffer.size(), maia::type_traits<MFloat>::mpiType(), rcvBuffer.data(),
              noToRecv.data(), displacements.data(), maia::type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_, "send",
              "recv");

  if(domainId() != 0) return;

  MString delimiter = " ";
  std::stringstream logEntry;
  std::stringstream logEntryAng;

  for(MInt b = 0; b < m_noEmbeddedBodies; b++) {
    MInt bodyId = (MInt)rcvBuffer[MLong(b * bufSizePerBody)];

    logEntry << globalTimeStep << delimiter;
    logEntry << bodyId << delimiter;

    logEntryAng << globalTimeStep << delimiter;
    logEntryAng << bodyId << delimiter;

    std::array<MFloat, nDim> center{};
    std::array<MFloat, nDim> velocity{};
    std::array<MFloat, nDim> acceleration{};
    std::array<MFloat, nDim> ang_vel{};
    std::array<MFloat, nDim> torque{};

    for(MInt n = 0; n < nDim; n++) {
      center[n] = rcvBuffer[b * bufSizePerBody + 1 + n];
      velocity[n] = rcvBuffer[b * bufSizePerBody + 1 + nDim + n];
      acceleration[n] = rcvBuffer[b * bufSizePerBody + 1 + 2 * nDim + n];
      ang_vel[n] = rcvBuffer[b * bufSizePerBody + 1 + 3 * nDim + n];
      torque[n] = rcvBuffer[b * bufSizePerBody + 1 + 4 * nDim + n];
    }

    // bodies log
    for(MInt n = 0; n < nDim; n++) {
      logEntry << std::scientific << center[n] << delimiter;
    }

    for(MInt n = 0; n < nDim; n++) {
      logEntry << std::scientific << velocity[n] << delimiter;
    }

    for(MInt n = 0; n < nDim; n++) {
      logEntry << std::scientific << acceleration[n] << delimiter;
    }

    // angvel log
    for(MInt n = 0; n < nDim; n++) {
      logEntryAng << std::scientific << ang_vel[n] << delimiter;
    }

    for(MInt n = 0; n < nDim; n++) {
      logEntryAng << std::scientific << torque[n] << delimiter;
    }

    logEntry << "\n";
    logEntryAng << "\n";
  }

  m_log.open("bodies.log", std::ios::app);
  m_log << logEntry.str();
  m_log.close();

  m_anglog.open("angvel.log", std::ios::app);
  m_anglog << logEntryAng.str();
  m_anglog.close();
}

/**
 * \brief Copy the body data for time t to t-1 to prepare the next timestep
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::advanceBodies() {
  TRACE();

  m_bodies.advanceBodies();
}

/**
 * \brief Peform a prediction for all bodies
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::predictorStep() {
  TRACE();

  const MFloat dt = m_timestep;

  if(m_translation) {
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      // Translational values
      for(MInt n = 0; n < nDim; n++) {
        // Determine acceleration
        a_bodyAcceleration(bodyId, n) = a_bodyAccelerationOld(bodyId, n);

        // Determine velocities
        a_bodyVelocity(bodyId, n) = a_bodyVelocityOld(bodyId, n)
                                    + 0.5 * dt * (a_bodyAcceleration(bodyId, n) + a_bodyAccelerationOld(bodyId, n));

        // Determine position
        a_bodyCenter(bodyId, n) =
            a_bodyCenterOld(bodyId, n) + dt * a_bodyVelocityOld(bodyId, n)
            + 0.25 * POW2(dt) * (a_bodyAcceleration(bodyId, n) + a_bodyAccelerationOld(bodyId, n));
      }

      if(m_rotation) {
        predictRotation(bodyId);
      }

      // Scalar values
      // Determine temperature
      a_bodyTemperature(bodyId) =
          a_bodyTemperatureOld(bodyId) + dt * a_bodyHeatFlux(bodyId) / (m_bodyHeatCapacity * a_bodyMass(bodyId));
    }
  } else {
    if(m_rotation) {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        predictRotation(bodyId);
      }
    }
  }
}

template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsRemoteToHome(std::vector<MInt>& noToRecv, std::vector<MInt>& noToSend) {
  TRACE();

  std::vector<MPI_Request> mpi_send_req(noHomeDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noRemoteDomains(), MPI_REQUEST_NULL);

  // recv from every neighbor
  MInt count = 0;
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    MPI_Irecv(&noToRecv[count], 1, MPI_INT, remoteDomain, 2, mpiComm(), &mpi_recv_req[count], AT_, "recv from remote");
    count++;
  }

  // send
  count = 0;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    MPI_Isend(&noToSend[count], 1, MPI_INT, homeDomain, 2, mpiComm(), &mpi_send_req[count], AT_, "send to home");
    count++;
  }

  // wait for communication to finish
  for(MInt n = 0; n < noHomeDomains(); n++) {
    MPI_Wait(&mpi_send_req[n], MPI_STATUS_IGNORE, AT_);
  }

  for(MInt n = 0; n < noRemoteDomains(); n++) {
    MPI_Wait(&mpi_recv_req[n], MPI_STATUS_IGNORE, AT_);
  }
}

/**
 * \brief exchange of Buffer legnths for further exchanges
 *
 * serves as an mpi probe function
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsNeighbor(std::vector<MInt>& noToRecv,
                                                      std::vector<std::vector<MInt>>& bodyList) {
  TRACE();

  // exchange no of edge bodies to neighboring domains
  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // recv from every neighbor
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Irecv(&noToRecv[n], 1, MPI_INT, neighborDomain(n), 2, mpiComm(), &mpi_recv_req[n], AT_, "recv predict");
  }

  // create send buffer
  std::vector<MUint> sendBuffer(noNeighborDomains(), 0);
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    sendBuffer[d] = bodyList[d].size();
  }

  // send to every neighbor
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Isend(&sendBuffer[i], 1, MPI_INT, neighborDomain(i), 2, mpiComm(), &mpi_send_req[i], AT_, "send predict");
  }

  // Wait for MPI exchange to complete
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }
}

/**
 * \brief exchange of Buffer legnths for further exchanges
 *
 * serves as an mpi probe function
 *
 * \param[in] noToRecv count of bodies to receive for a neighbor domain
 * \param[in] bodyList map like remoteDomainLocalBodies
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsNeighbor(std::vector<MInt>& noToRecv,
                                                      std::map<MInt, std::vector<MInt>>& bodyList) {
  TRACE();

  // exchange no of edge bodies to neighboring domains
  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // recv from every neighbor
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Irecv(&noToRecv[n], 1, MPI_INT, neighborDomain(n), 2, mpiComm(), &mpi_recv_req[n], AT_, "recv predict");
  }

  // create send buffer
  std::vector<MUint> sendBuffer(noNeighborDomains(), 0);
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    MInt globDomain = neighborDomain(d);
    sendBuffer[d] = bodyList[globDomain].size();
  }

  // send to every neighbor
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Isend(&sendBuffer[i], 1, MPI_INT, neighborDomain(i), 2, mpiComm(), &mpi_send_req[i], AT_, "send predict");
  }

  // Wait for MPI exchange to complete
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }
}

/**
 * \brief exchange of Buffer lengths all to root
 *
 * serves as an mpi probe function
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsAllToRoot(std::vector<MInt>& noToRecv, const MInt noToSend) {
  TRACE();

  MPI_Request mpi_send_req = MPI_REQUEST_NULL;
  std::vector<MPI_Request> mpi_recv_req(noDomains(), MPI_REQUEST_NULL);

  // recv from every domain and wait
  if(!domainId()) {
    for(MInt n = 0; n < noDomains(); n++) {
      MPI_Irecv(&noToRecv[n], 1, MPI_INT, n, 2, mpiComm(), &mpi_recv_req[n], AT_, "recv buffer length");
    }
  }

  // send to root
  MPI_Isend(&noToSend, 1, MPI_INT, 0, 2, mpiComm(), &mpi_send_req, AT_, "send buffer length");

  // Wait for MPI exchange to complete
  MPI_Waitall(noDomains(), mpi_recv_req.data(), MPI_STATUS_IGNORE, AT_);

  // Wait for MPI exchange to complete
  MPI_Wait(&mpi_send_req, MPI_STATUS_IGNORE, AT_);
}

/**
 * \brief exchange of Buffer lengths all to all
 *
 * serves as an mpi probe function
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsAllToAll(std::vector<MInt>& noToRecv, const MInt noToSend) {
  TRACE();
  // exchange no of edge bodies to neighboring domains
  std::vector<MPI_Request> mpi_send_req(noDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noDomains(), MPI_REQUEST_NULL);

  // recv from every domain and wait
  for(MInt n = 0; n < noDomains(); n++) {
    MPI_Irecv(&noToRecv[n], 1, MPI_INT, n, 2, mpiComm(), &mpi_recv_req[n], AT_, "recv buffer length");
  }

  // send to root
  for(MInt n = 0; n < noDomains(); n++) {
    MPI_Isend(&noToSend, 1, MPI_INT, n, 2, mpiComm(), &mpi_send_req[n], AT_, "send buffer length");
  }
  // Wait for MPI exchange to complete
  MPI_Waitall(noDomains(), mpi_send_req.data(), MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noDomains(), mpi_recv_req.data(), MPI_STATUSES_IGNORE, AT_);
}

/**
 * \brief exchange of Buffer lengths all to all in direct neighborhood
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBufferLengthsNeighbors(std::vector<MInt>& noToRecv, const MInt noToSend) {
  TRACE();
  // exchange no of edge bodies to neighboring domains
  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // recv from every domain and wait
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Irecv(&noToRecv[n], 1, MPI_INT, neighborDomain(n), 2, mpiComm(), &mpi_recv_req[n], AT_, "recv buffer length");
  }

  // send to root
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Isend(&noToSend, 1, MPI_INT, neighborDomain(n), 2, mpiComm(), &mpi_send_req[n], AT_, "send buffer length");
  }
  // Wait for MPI exchange to complete
  MPI_Waitall(noNeighborDomains(), mpi_send_req.data(), MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborDomains(), mpi_recv_req.data(), MPI_STATUSES_IGNORE, AT_);
}


/**
 *  \brief for new restartFile format
 *
 *
 *
 */
template <MInt nDim>
MLong RigidBodies<nDim>::getDomainOffset() {
  TRACE();

  std::vector<MInt> noLocalBodiesPerPartition(noDomains());
  exchangeBufferLengthsAllToAll(noLocalBodiesPerPartition, m_noLocalBodies);

  MLong offset{};
  for(MInt i = 0; i < domainId(); i++) {
    offset += noLocalBodiesPerPartition[i];
  }

  return offset;
}


/**
 * \brief exchange of Body Variables all to all
 *
 * every domain has every current data for every body relevant to perform a solution step
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeBodyVariablesAllToAll() {
  TRACE();

  m_bResFile.reset(noEmbeddedBodies());
  m_bResFile.append(noEmbeddedBodies());

  // variables
  std::vector<MFloat*> transVarsRcv{&m_bResFile.bodyCenter(0, 0), &m_bResFile.bodyVelocity(0, 0),
                                    &m_bResFile.bodyAcceleration(0, 0)};

  std::vector<MFloat*> rotVarsRcv{&m_bResFile.angularVelocityBodyT1B2(0, 0), &m_bResFile.angularAccelerationBody(0, 0)};

  std::vector<MFloat*> quatVarsRcv{&m_bResFile.bodyQuaternionT1B2(0, 0), &m_bResFile.bodyQuaternionT1(0, 0)};

  // bodyId + temperature + trans + rot + quat
  const MInt noTransVars = 3;
  const MInt noRotVars = 2;
  const MInt noQuatVars = 2;
  const MInt bufSizePerBody = 1 + 1 + noTransVars * nDim + noRotVars * nRot + noQuatVars * nQuat;
  MInt varSizeCount = 0;
  MInt bodySizeCount = 0;
  // create send buffer
  std::vector<MFloat> sendBuffer(MLong(m_noLocalBodies * bufSizePerBody));
  const MUint noToSend = sendBuffer.size();

  // fill send buffer, if bodies are in collector
  if(m_bodies.size() > 0) {
    // send variables
    std::vector<MFloat*> transVars{&a_bodyCenter(0, 0), &a_bodyVelocity(0, 0), &a_bodyAcceleration(0, 0)};
    // MInt noTransVars = transVars.size();

    std::vector<MFloat*> rotVars{&a_angularVelocityBodyT1B2(0, 0), &a_angularAccelerationBody(0, 0)};
    // MInt noRotVars = rotVars.size();

    std::vector<MFloat*> quatVars{&a_bodyQuaternionT1B2(0, 0), &a_bodyQuaternionT1(0, 0)};
    // MInt noQuatVars = quatVars.size();

    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      sendBuffer[bodySizeCount] = (MFloat)getGlobalBodyId(bodyId);
      varSizeCount += 1;

      // temperature
      sendBuffer[bodySizeCount + varSizeCount] = a_bodyTemperature(bodyId);
      varSizeCount += 1;

      // trans vars
      for(MInt v = 0; v < noTransVars; v++) {
        for(MInt n = 0; n < nDim; n++) {
          sendBuffer[bodySizeCount + varSizeCount + n] = transVars[v][nDim * bodyId + n];
        }
        varSizeCount += nDim;
      }
      //  rot vars
      for(MInt v = 0; v < noRotVars; v++) {
        for(MInt r = 0; r < nRot; r++) {
          sendBuffer[bodySizeCount + varSizeCount + r] = rotVars[v][nRot * bodyId + r];
        }
        varSizeCount += nRot;
      }

      // quat vars
      for(MInt v = 0; v < noQuatVars; v++) {
        for(MInt q = 0; q < nQuat; q++) {
          sendBuffer[bodySizeCount + varSizeCount + q] = quatVars[v][nQuat * bodyId + q];
        }
        varSizeCount += nQuat;
      }
      bodySizeCount += bufSizePerBody;
      varSizeCount = 0;
    }
  }

  // allocate receive buffer
  std::vector<MInt> noToRecv(noDomains());
  exchangeBufferLengthsAllToAll(noToRecv, noToSend);
  std::vector<MFloat> rcvBuffer(noEmbeddedBodies() * bufSizePerBody);

  // displacements
  MInt count = 0;
  std::vector<MInt> displacements(noDomains());
  for(MInt n = 0; n < noDomains(); n++) {
    displacements[n] = count;
    count += noToRecv[n];
  }

  // exchange data
  MPI_Allgatherv(sendBuffer.data(), sendBuffer.size(), maia::type_traits<MFloat>::mpiType(), rcvBuffer.data(),
                 noToRecv.data(), displacements.data(), maia::type_traits<MFloat>::mpiType(), mpiComm(), AT_, "send",
                 "recv");

  // unpack receive buffer
  varSizeCount = 0;
  bodySizeCount = 0;
  for(MInt i = 0; i < noEmbeddedBodies(); i++) {
    const MInt bodyId = (MInt)rcvBuffer[bodySizeCount];
    varSizeCount += 1;

    // temperature
    m_bResFile.bodyTemperature(bodyId) = rcvBuffer[bodySizeCount + varSizeCount];
    varSizeCount += 1;

    // trans vars
    for(MInt v = 0; v < noTransVars; v++) {
      for(MInt n = 0; n < nDim; n++) {
        transVarsRcv[v][nDim * bodyId + n] = rcvBuffer[bodySizeCount + varSizeCount + n];
      }
      varSizeCount += nDim;
    }
    //  rot vars
    for(MInt v = 0; v < noRotVars; v++) {
      for(MInt r = 0; r < nRot; r++) {
        rotVarsRcv[v][nRot * bodyId + r] = rcvBuffer[bodySizeCount + varSizeCount + r];
      }
      varSizeCount += nRot;
    }

    // quat vars
    for(MInt v = 0; v < noQuatVars; v++) {
      for(MInt q = 0; q < nQuat; q++) {
        quatVarsRcv[v][nQuat * bodyId + q] = rcvBuffer[bodySizeCount + varSizeCount + q];
      }
      varSizeCount += nQuat;
    }
    bodySizeCount += bufSizePerBody;
    varSizeCount = 0;
  }
}

/**
 * \brief searches in all initial Bodies for local bodies
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::findLocalBodies() {
  TRACE();

  // loop over all bodies initial bodies
  for(MInt globalBodyId = 0; globalBodyId < noEmbeddedBodies(); globalBodyId++) {
    std::array<MFloat, nDim> tmpCenter{};
    for(MInt n = 0; n < nDim; n++) {
      if(m_restart) {
        tmpCenter[n] = m_bResFile.bodyCenter(globalBodyId, n);
      } else {
        tmpCenter[n] = m_initialBodyCenters[globalBodyId * nDim + n];
      }
    }
    MBool insideLocalBbox = grid().raw().pointInLocalBoundingBox(tmpCenter.data());
    if(insideLocalBbox) {
      MInt localPCellId = grid().raw().findContainingPartitionCell(tmpCenter.data(), -1, nullptr);
      if(localPCellId != -1) {
        m_globalBodyIds.emplace_back(globalBodyId);
        m_noLocalBodies++;
      }
    }
  }
}

template <MInt nDim>
void RigidBodies<nDim>::initGridProperties() {
  TRACE();
  // Counting periodic directions
  MInt periodic = 0;
  for(MInt n = 0; n < nDim; n++) {
    m_globDomainLength[n] = globBBox(nDim + n) - globBBox(n);
    periodic += grid().raw().periodicCartesianDir(n);
  }

  if(periodic) {
    m_periodicBC = true;

    constexpr MFloat eps = 1e-9;
    MBool domainIsCube = true;
    for(MInt n = 0; n < nDim; n++) {
      domainIsCube = domainIsCube && (m_globDomainLength[0] - m_globDomainLength[n]) < eps;
    }

    if(!domainIsCube) {
      std::cerr << "\033[0;31m Warning:\033[0m Domain is not a cube - periodic BCs only work on cube domains"
                << std::endl;
    }
  }

  updateMaxLevel(maxLevel());
}

template <MInt nDim>
void RigidBodies<nDim>::updateInfoDiameter(const MBool initCall) {
  TRACE();
  std::vector<MInt> l_Bodies;
  MFloat maxBodyExtension = F0;

  if(initCall) {
    for(MInt i = 0; i < m_noLocalBodies; i++) {
      l_Bodies.push_back(i);
    }
  } else {
    for(MInt bodyId = 0; bodyId < noCollectorBodies(); bodyId++) {
      l_Bodies.push_back(bodyId);
    }
  }

  if(m_initBodyRadii.empty()) {
    for(const MInt& bodyId : l_Bodies) {
      maxBodyExtension = std::max(maxBodyExtension, 2 * a_bodyRadius(bodyId));
    }
  } else {
    for(const MInt& bodyId : l_Bodies) {
      for(MInt n = 0; n < nDim; n++) {
        maxBodyExtension = std::max(maxBodyExtension, 2 * a_bodyRadii(bodyId, n));
      }
    }
  }
  // Body extension + buffer
  constexpr const MInt noCellsClearance = 4.0;
  m_infoDiameter = maxBodyExtension + noCellsClearance * c_cellLengthAtMaxLevel();
}

/**
 * \brief creates list with all edge bodies that are no longer contained within local domain
 *
 * + removes them from local & edge body list & body mapping
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::findTransferBodies() {
  TRACE();
  m_transferBodies.clear();
  m_transferBodies.resize(noNeighborDomains());

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    /* 1) if local bodies center is contained in a neighbor domains halo cell
    ** move local body -> transfer body
    */
    for(const MInt& bodyId : m_remoteDomainLocalBodies[neighborDomain(i)]) {
      MBool outside = -1 != grid().raw().findContainingHaloCell(&a_bodyCenter(bodyId, 0), -1, i, true, nullptr);
      if(outside) {
        // Already get new index of body to be transfered in collector in transferBodies
        m_transferBodies[i].push_back(bodyId);
#ifdef RB_DEBUG
        m_debugFileStream << "Found transfer body with globalId " << getGlobalBodyId(bodyId)
                          << " will be sent to domain " << i << std::endl;
#endif
      }
    }
  }
}

/**
 * \brief exchanges transfer bodies
 *
 * sends them to new home domain +
 * receives them from old home domain
 *
 * puts them into local bodies list (comperable state as in find local bodies)
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
// TODO, just loop over old remote Domains of body for sending
template <MInt nDim>
void RigidBodies<nDim>::exchangeTransferBodies() {
  TRACE();
  // transfer to new home domain
  // exchange no of transfer bodies to neighboring domains
  std::vector<MInt> noToRecv(noNeighborDomains());
  MInt totalBufferLength = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    totalBufferLength += m_transferBodies[i].size();
  }

  exchangeBufferLengthsNeighbors(noToRecv, totalBufferLength);

  // send new home domain Id to every other home domain
  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // create receive buffer
  std::vector<std::vector<MFloat>> rcvBufferAll(noNeighborDomains());
  MInt bufSizePerBody = 2;

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    rcvBufferAll[n].resize(MLong(noToRecv[n] * bufSizePerBody));
  }

  // non-blocking receive
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    if(!noToRecv[i]) {
      continue;
    }

    MPI_Irecv(rcvBufferAll[i].data(), rcvBufferAll[i].size(), maia::type_traits<MFloat>::mpiType(), neighborDomain(i),
              2, mpiComm(), &mpi_recv_req[i], AT_, "recv predict");
  }

  // create sendBuffer that is send to every domain
  MInt idx = 0;
  std::vector<MFloat> sendBufferAll(MLong(totalBufferLength * bufSizePerBody));
  std::list<std::pair<MInt, MInt>> localBodyRemoveList;
  for(MUint i = 0; i < m_transferBodies.size(); i++) {
    MUint domain = i;
    MUint newHomeDomain = neighborDomain(domain);
    for(const auto& body : m_transferBodies[i]) {
      sendBufferAll[MLong(idx * bufSizePerBody)] = (MFloat)getGlobalBodyId(body);
      sendBufferAll[idx * bufSizePerBody + 1] = (MFloat)newHomeDomain; // <- send new globalDomainId
      idx++;

      auto it = std::find_if(localBodyRemoveList.cbegin(), localBodyRemoveList.cend(),
                             [&body](const std::pair<MInt, MInt>& localB) { return localB.second == body; });
      if(it == localBodyRemoveList.cend()) {
        localBodyRemoveList.push_front(std::make_pair(newHomeDomain, body));
      }
    }
  }

  // non-blocking send
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    if(sendBufferAll.empty()) {
      continue;
    }

    MPI_Isend(sendBufferAll.data(), sendBufferAll.size(), maia::type_traits<MFloat>::mpiType(), neighborDomain(i), 2,
              mpiComm(), &mpi_send_req[i], AT_, "send predict");
  }

  // wait for communication to finish
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
  }

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }

  // transferBodies change from local to remote status -> swap position in collector
  // only executed on sending domain
  // update mappings and collector accordingly

  // Due to swapping, localIds in localBodyRemoveList might lose validity
  // This mapping starts as unity and is filled if a localBody becomes a swap target
  std::map<MInt, MInt> tmpMapping;
  for(const auto& [newDomainId, body] : localBodyRemoveList) {
    tmpMapping[body] = body;
  }
  for(const auto& [newDomainId, body] : localBodyRemoveList) {
    const MInt tmpCollectorId = tmpMapping[body];
    const MInt remoteBodyId = localToRemoteBody(tmpCollectorId, body, newDomainId);

    // Fill mapping if one of the remaining localBodies was swapped
    for(auto& [newDomainId_, body_] : localBodyRemoveList) {
      if(body_ == remoteBodyId) {
        tmpMapping[body_] = tmpCollectorId;
      }
    }
  }

  // unpack buffer
  // 1) Find all bodies, that are transfered to the current domain (remoteBodiesToLocal)
  // 2) Find all bodies, that changed there homeDomain in the same timestep, and update mapping accordingly, if relevant
  std::vector<std::pair<MInt, MInt>> remoteBodiesToLocal;
  std::vector<std::array<MInt, 2>> remoteBodiesNeighDomain;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt b = 0; b < noToRecv[n]; b++) {
      MInt globalBodyId = MInt(rcvBufferAll[n][MLong(b * bufSizePerBody)]);
      MInt newHomeDomain = MInt(rcvBufferAll[n][MLong(b * bufSizePerBody + 1)]);
      const MInt localId = getLocalBodyId(globalBodyId);

#ifdef RB_DEBUG
      m_debugFileStream << "received: " << globalBodyId << " with localId " << localId << " and new homeDomain "
                        << newHomeDomain << "from " << neighborDomain(n) << std::endl;
#endif

      if(newHomeDomain == domainId()) {
        if(localId == -1) {
          mTerm(1, AT_,
                "ERROR: Body " + std::to_string(globalBodyId) + " has not been on domain " + std::to_string(domainId())
                    + " before !");
        }
        MInt oldHomeDomainId = neighborDomain(n);
        remoteBodiesToLocal.emplace_back(oldHomeDomainId, globalBodyId);
      } else if(localId != -1) {
        remoteBodiesNeighDomain.push_back({newHomeDomain, globalBodyId});
      }
    }
  }

  // 3) make received bodies local
  // If body was previously remote on domain, it gets now local -> change position of body in collector of
  // receiving domain (append to local bodies)
  for(MUint i = 0; i < remoteBodiesToLocal.size(); i++) {
    if(m_noRemoteBodies <= 0) {
      mTerm(1, AT_, "Domain: " + std::to_string(domainId()) + " has no remote bodies that can change to local status!");
    }
    MInt oldHomeDomainId = remoteBodiesToLocal[i].first;
    MInt globalId = remoteBodiesToLocal[i].second;
    MInt localId = getLocalBodyId(globalId);

    const MInt newIdx = m_noLocalBodies;

#ifdef RB_DEBUG
    m_debugFileStream << "B localId: " << localId << " globalId " << globalId << " newIdx " << newIdx
                      << " globalIdNewIdx: " << getGlobalBodyId(newIdx) << " oldHomeDomain " << oldHomeDomainId
                      << std::endl;
#endif

    if(noConnectingBodies() > 1) {
      m_bodies.swap(localId, newIdx);
      iter_swap(m_globalBodyIds.begin() + localId, m_globalBodyIds.begin() + newIdx);
    }
    m_noLocalBodies++;
    m_noRemoteBodies--;

    // we just append to local bodies and therefore do not need to update this mapping
    m_remoteDomainLocalBodies[oldHomeDomainId].push_back(newIdx);

    // update localIds in mapping
    std::map<MInt, std::vector<MInt>> oldMapping(m_homeDomainRemoteBodies);
    for(auto& [domain, bodies] : oldMapping) {
      // first delete the new local Body from Remote Body mapping -> body is now on current domain
      for(const MInt& body : bodies) {
        if(body == localId) {
#ifdef RB_DEBUG
          m_debugFileStream << "delete body " << body << " domain " << domain << std::endl;
#endif
          std::vector<MInt>* c = &m_homeDomainRemoteBodies[domain];
          c->erase(std::remove(c->begin(), c->end(), localId), c->end());
        }
      }
      // when making a local body remote, all the ids of the remote bodies, that have a lower id than the new local body
      // have to be modified. This is done in a second loop to prevent accidental deletion of a body
      // if you want to remove local id 1, local id 0 is raised to local id 1
      for(MInt& body : m_homeDomainRemoteBodies[domain]) {
        if(body < localId) {
          body++;
        }
      }
    }
  }

#ifdef RB_DEBUG
  printBodyDomainConnections();
#endif

  // 4) now after all localBodies and the corresponding remoteDomains are updated
  // we potentielly need to update the homeDomain of our remote Bodies
  // search in mapping to find body and update new HomeDomain
  MBool found = false;
  for(const auto& arr : remoteBodiesNeighDomain) {
    const MInt newHomeDomain = arr[0];
    const MInt globalId = arr[1];

    const MInt NewLocalId = getLocalBodyId(globalId);

    for(const auto& [domain, bodies] : m_homeDomainRemoteBodies) {
      if(domain == domainId()) continue;
      for(const auto& bodyId : bodies) {
        if(getGlobalBodyId(bodyId) == globalId) {
          found = true;
          // a remote Body was transfered to a new homedomain, update:
          if(domain != newHomeDomain) {
            std::vector<MInt>* c = &m_homeDomainRemoteBodies[domain];
            c->erase(std::remove(c->begin(), c->end(), NewLocalId), c->end());

            m_homeDomainRemoteBodies[newHomeDomain].push_back(NewLocalId);
            break;
          }
#ifdef RB_DEBUG
          m_debugFileStream
              << "something went wrong, domain should be diff from newHomeDomain in remoteBodiesNeighDomain! "
              << "Domain is " << domain << " newHomeDomain is " << newHomeDomain << std::endl;
#endif
        }
      }
      if(found) break;
    }
  }
#ifdef RB_DEBUG
  m_debugFileStream << "end of exchange: " << std::endl;
  printBodyDomainConnections();
#endif
}


/***
 * \brief checks if body center leaves domain and the corresponding dummybody center enters therefore the domain
 *
 * \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 */
template <MInt nDim>
void RigidBodies<nDim>::checkDummyBodiesForSwap() {
  TRACE();

  if(m_associatedDummyBodies.empty()) return;

  for(const auto& [bodyId, dummyId] : m_associatedDummyBodies) {
    MBool outside = -1 != grid().raw().findContainingHaloCell(&a_bodyCenter(bodyId, 0), -1, domainId(), true, nullptr);
    if(outside) {
      // swap body and its associated dummybody
      m_bodies.swap(bodyId, dummyId);
#ifdef RB_DEBUG
      m_debugFileStream << "body " << bodyId << " and associated dummyBody " << dummyId << " were swapped!"
                        << std::endl;
#endif
    }
  }
}


/**
 * \brief updates local/ edge body domain connections
 *
 * updates all local lists + mappings
 *
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>s
 */
template <MInt nDim>
void RigidBodies<nDim>::updateBodyDomainConnections(MBool /*initCall*/) {
  TRACE();

#ifdef RB_DEBUG
  m_debugFileStream << "entering updateBodyDomainConnections on domain " << domainId() << std::endl;
  printBodyDomainConnections();
#endif

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    /* 1) if local body are intersecting with Halo Partition Cells
    ** copy local body -> edge body
    */
    MInt globDomain = neighborDomain(i);
    for(MInt bodyId = 0; bodyId < noConnectingBodies(); bodyId++) {
      MBool alreadyInside = std::binary_search(m_remoteDomainLocalBodies[globDomain].begin(),
                                               m_remoteDomainLocalBodies[globDomain].end(), bodyId);

      MBool alreadySelfMapped = m_periodicBC && hasAssociatedDummyBody(bodyId);

#ifdef RB_DEBUG
      m_debugFileStream << "D" << globDomain << " lBody " << bodyId << " already in mapping " << alreadyInside
                        << " selfmapped " << alreadySelfMapped << std::endl;
      for(auto&& b : m_remoteDomainLocalBodies[globDomain]) {
        m_debugFileStream << "[" << globDomain << "," << b << "]" << std::endl;
      }
      m_debugFileStream << std::endl;
#endif
      if(alreadyInside || alreadySelfMapped) {
        continue;
      }

      MInt intersectingCellId =
          grid().raw().intersectingWithHaloCells(&a_bodyCenter(bodyId, 0), m_infoDiameter / 2.0, i, true);
      if(intersectingCellId != -1) {
        if(domainId() == globDomain && m_periodicBC
           && grid().raw().a_hasProperty(intersectingCellId, Cell::IsPeriodic)) {
          // Hande Self-Mapping Case: A domain does not have to send something to itself, so the dummyBody does not
          // appear in the domain - body mappings, but a dummy has still to be added to collector, as the body has now
          // to different positions of the bodycenter, in the periodic case
          createDummyBody(bodyId);
        } else if(isLocalBody(bodyId)) {
          // Regular handling of edgeBodies, that have a remoteDomain
          std::vector<MInt>* bodyIds = &(m_remoteDomainLocalBodies[globDomain]);
          if(!std::binary_search(bodyIds->begin(), bodyIds->end(), bodyId)) {
            // Insert at the right position to keep this mapping sorted by the local body id
            bodyIds->insert(std::upper_bound(bodyIds->begin(), bodyIds->end(), bodyId), bodyId);
          }
        }
      }
    }


    /* 2) if edge body is no longer intersecting halo partition cells
    ** remove edge body -> []
    */
    std::vector<MInt> removeList;
    for(const MInt& bodyId : m_remoteDomainLocalBodies[globDomain]) {
      // checking only partition halo cells
      MInt intersectingCellId =
          grid().raw().intersectingWithHaloCells(&a_bodyCenter(bodyId, 0), m_infoDiameter / 2.0, i, true);
      if(intersectingCellId == -1) {
        // body is no longer intersecting haloPartitionCell (body has left domain ->
        // body that was previously marked as transferBody and now needs to be removed from edgeBodies of homeDomain)
        removeList.push_back(bodyId);
      }
    }

#ifdef RB_DEBUG
    if(!removeList.empty()) {
      MString removeListS = "[";
      for(const MInt& body : removeList) {
        removeListS += std::to_string(body) + ", ";
      }
      removeListS += "] ";
      m_debugFileStream << "RemoveList for domain " << domainId() << " " << removeListS
                        << " neighborDomain: " << neighborDomain(i) << std::endl;
    }
#endif

    for(const auto& body : removeList) {
      // body is not remote on other domain, remove it from mapping and edgeBodies
      std::vector<MInt>* bodies = &m_remoteDomainLocalBodies[neighborDomain(i)];
      bodies->erase(std::remove(bodies->begin(), bodies->end(), body), bodies->end());
    }
  }

  /* 3) if edge body intersects a periodic halo partition cell
  ** calculate periodic shift
  */
  if(m_periodicBC) {
    m_periodicShift.clear();
    // determine periodicShift for bodies that are remote on Neighbors
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MInt globDomain = neighborDomain(i);

      if(globDomain == domainId()) continue;

      for(const auto& bodyId : m_remoteDomainLocalBodies[globDomain]) {
        std::array<MFloat, nDim> center{};
        for(MInt n = 0; n < nDim; n++) {
          center[n] = a_bodyCenter(bodyId, n);
        }
        const MInt intersectingCellId =
            grid().raw().intersectingWithHaloCells(center.data(), m_infoDiameter / 2.0, i, true);
        std::array<MFloat, nDim> tmpShift{};
        for(MInt n = 0; n < nDim; n++) {
          tmpShift[n] = calculatePeriodicShift(intersectingCellId, n);
        }
        const MInt globalDomain = neighborDomain(i);
        m_periodicShift[globalDomain][bodyId] = tmpShift;
      }
    }

    // Now handle Self-Mapped Bodies
    for(const auto& [localBodyId, dummyBodyId] : m_associatedDummyBodies) {
      std::array<MFloat, nDim> center{};
      for(MInt n = 0; n < nDim; n++) {
        center[n] = a_bodyCenter(localBodyId, n);
      }
      const MInt intersectingCellId =
          grid().raw().intersectingWithHaloCells(center.data(), m_infoDiameter / 2.0, 0, true);
      std::array<MFloat, nDim> tmpShift{};
      for(MInt n = 0; n < nDim; n++) {
        tmpShift[n] = calculatePeriodicShift(intersectingCellId, n);
      }
      m_periodicShift[domainId()][localBodyId] = tmpShift;
    }

#ifdef RB_DEBUG
    MString periodicShiftS{};
    for(const auto& [domain, bodies] : m_periodicShift) {
      periodicShiftS += "on domain " + std::to_string(domain) + " periodic Shifts are: \n";
      for(const auto& [bodyId, shift] : bodies) {
        periodicShiftS += "  bodyId: " + std::to_string(bodyId) + " [";
        for(MInt n = 0; n < nDim; n++) {
          periodicShiftS += std::to_string(shift[n]) + " ";
        }
        periodicShiftS += "] \n";
      }
    }
    m_debugFileStream << periodicShiftS;
#endif
  }
}

/**
 * \brief appends new DummyBody to collector
 *
 * \author Johannes Grafen <johannes.grafen@rwth-aachen.de>s
 */
template <MInt nDim>
void RigidBodies<nDim>::createDummyBody(MInt bodyId) {
  TRACE();

  MInt newBodyId = noCollectorBodies();
  m_noDummyBodies++;
  m_associatedDummyBodies[bodyId] = newBodyId;
  m_globalBodyIds.push_back(getGlobalBodyId(bodyId));
  m_bodies.append();
#ifdef RB_DEBUG
  m_debugFileStream << "Self-mapping case, adding dummy body to collector at " << newBodyId
                    << " associated with bodyId " << bodyId << std::endl;
#endif
  // Initialize the dummyBuddy in collector
  a_bodyRadius(newBodyId) = a_bodyRadius(bodyId);
  a_bodyDensityRatio(newBodyId) = a_bodyDensityRatio(bodyId);
  a_bodyTemperature(newBodyId) = a_bodyTemperature(bodyId);
  a_bodyHeatFlux(newBodyId) = a_bodyHeatFlux(bodyId);
  for(MInt n = 0; n < nDim; n++) {
    a_bodyRadii(newBodyId, n) = a_bodyRadii(bodyId, n);
    a_bodyInertia(newBodyId, n) = a_bodyInertia(bodyId, n);
  }
}

/**
 * \brief periodic shift for position for further data exchange
 *
 * shift will always point inwards
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *
 * \param direction[in][3]
 * \param intersectingCellId[in]
 *
 * \returns shift for given direction depending on the intersecting halo cell id
 */
template <MInt nDim>
MFloat RigidBodies<nDim>::calculatePeriodicShift(MInt intersectingCellId, MInt direction) {
  TRACE();
  MFloat shift = F0;
  const MInt n = direction;
  const MBool periodicDir = grid().raw().periodicCartesianDir(n);
  const MBool periodicCell = grid().raw().a_hasProperty(intersectingCellId, Cell::IsPeriodic);
  if(m_periodicBC && periodicDir && periodicCell) {
    if(grid().raw().periodicCartesianDir(n)) {
      const MFloat cellCoordinate = grid().raw().a_coordinate(intersectingCellId, n);
      const MInt sign1 = maia::math::sgn(globBBox(n) - cellCoordinate);
      const MInt sign2 = maia::math::sgn(globBBox(n + nDim) - cellCoordinate);
      MFloat sign = F0;
      if(sign1 == 1 && sign2 == 1) {
        sign = 1.0;
      } else if(sign1 == -1 && sign2 == -1) {
        sign = -1.0;
      }
      shift = m_globDomainLength[n] * sign;
    }
  }
  return shift;
}

/**
 * \brief exchanges relevant body connections with new neighbor domains
 *
 * and updates remaining (remote) body lists & mappings
 *
 * \author Oliver Groll <oliver.groll@rwth-aachen.de>
 *         Johannes Grafen <johannes.grafen@rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::exchangeEdgeBodyData() {
  TRACE();
  // exchange no of edge bodies to neighboring domains
  std::vector<MInt> noToRecv(noNeighborDomains());
  exchangeBufferLengthsNeighbor(noToRecv, m_remoteDomainLocalBodies);

  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // create receive buffer
  std::vector<std::vector<MFloat>> rcvBufferAll(noNeighborDomains());
  MInt bufSizePerBody = 1 + nDim * 1 + nRot * 2 + nQuat * 2 + 3;

  // if not shape, different body Radii have to be send
  if(a_bodyType() != Shapes::Sphere) {
    bufSizePerBody += nDim;
  }

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    rcvBufferAll[n].resize(noToRecv[n] * bufSizePerBody);
  }

  // non-blocking receive
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    if(!noToRecv[i]) {
      continue;
    }

    MPI_Irecv(rcvBufferAll[i].data(), rcvBufferAll[i].size(), maia::type_traits<MFloat>::mpiType(), neighborDomain(i),
              2, mpiComm(), &mpi_recv_req[i], AT_, "recv predict");
  }

  // create send buffer
  MInt idx = 0;
  MInt paraCount = 0;
  std::vector<std::vector<MFloat>> sendBufferAll(noNeighborDomains());
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt globDomain = neighborDomain(i);
    if(domainId() == globDomain) continue;

    auto& tmpSendBuffer = sendBufferAll[i];
    tmpSendBuffer.resize(bufSizePerBody * m_remoteDomainLocalBodies[globDomain].size());
    idx = 0;
    for(const auto& body : m_remoteDomainLocalBodies[globDomain]) {
      paraCount = 0;
      // send corresponding globalId of body
      tmpSendBuffer[idx * bufSizePerBody] = (MFloat)getGlobalBodyId(body);
      paraCount++;

      for(MInt n = 0; n < nDim; n++) {
        tmpSendBuffer[idx * bufSizePerBody + paraCount + n] = a_bodyCenter(body, n);
      }
      paraCount += nDim;

      // necessary to exchange quaternions?
      for(MInt r = 0; r < nRot; r++) {
        tmpSendBuffer[idx * bufSizePerBody + paraCount + r] = a_angularAccelerationBody(body, r);
        tmpSendBuffer[idx * bufSizePerBody + paraCount + nRot + r] = a_angularVelocityBodyT1B2(body, r);
      }
      paraCount += 2 * nRot;

      for(MInt q = 0; q < nQuat; q++) {
        tmpSendBuffer[idx * bufSizePerBody + paraCount + q] = a_bodyQuaternionT1(body, q);
        tmpSendBuffer[idx * bufSizePerBody + paraCount + nQuat + q] = a_bodyQuaternionT1B2(body, q);
      }
      paraCount += 2 * nQuat;

      tmpSendBuffer[idx * bufSizePerBody + paraCount] = a_bodyRadius(body);
      paraCount++;
      if(a_bodyType() != Shapes::Sphere) {
        for(MInt n = 0; n < nDim; n++) {
          tmpSendBuffer[idx * bufSizePerBody + paraCount + n] = a_bodyRadii(body, n);
        }
        paraCount += nDim;
      }

      tmpSendBuffer[idx * bufSizePerBody + paraCount] = a_bodyDensityRatio(body);
      paraCount++;
      tmpSendBuffer[idx * bufSizePerBody + paraCount] = a_bodyTemperature(body);
      paraCount++;
      idx++;
    }
  }

  // non-blocking send
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt globDomain = neighborDomain(i);
    if(m_remoteDomainLocalBodies[globDomain].empty()) {
      continue;
    }

    MPI_Isend(sendBufferAll[i].data(), sendBufferAll[i].size(), maia::type_traits<MFloat>::mpiType(), neighborDomain(i),
              2, mpiComm(), &mpi_send_req[i], AT_, "send predict");
  }

  // wait for communication to finish
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt globDomain = neighborDomain(i);
    if(!m_remoteDomainLocalBodies[globDomain].empty()) {
      MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
    }

    if(noToRecv[i]) {
      MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
    }
  }

  // All bodies, that are marked as remote on domain are added to remove list
  // it is checked if data for the body that was already remote is received
  // if not it is not relevant anymore and can be deleted from the mapping
  std::list<std::pair<MInt, MInt>> remoteBodyRemoveList;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    for(const MInt rBody : bodies) {
      remoteBodyRemoveList.push_back(std::make_pair(homeDomain, rBody));
    }
  }

#ifdef RB_DEBUG
  MString removeListS = "[";
  for(const auto& body : remoteBodyRemoveList) {
    removeListS += "[" + std::to_string(body.first) + "," + std::to_string(body.second) + "]";
  }
  removeListS += "] ";
  m_debugFileStream << "Before: RemoveList for domain " << domainId() << " " << removeListS << std::endl;
  printBodyDomainConnections();
#endif

  // unpack buffer
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt globalDomain = neighborDomain(n);
    MInt mapIdx = 0;
    for(MInt b = 0; b < noToRecv[n]; b++) {
      paraCount = 0;
      MInt remoteBodyId = rcvBufferAll[n][b * bufSizePerBody];
      paraCount++;

#ifdef RB_DEBUG
      m_debugFileStream << "domain: " << domainId() << " receives remoteBodyId: " << remoteBodyId << " from "
                        << globalDomain << " with locId " << getLocalBodyId(remoteBodyId) << std::endl;
#endif

      MInt localBodyId = getLocalBodyId(remoteBodyId);
      if(localBodyId != -1) {
        // received data for localBodyId, body is still relevant and therefore has to be eliminated again from list
        std::list<std::pair<MInt, MInt>>* a = &remoteBodyRemoveList;
        a->erase(std::remove_if(a->begin(), a->end(),
                                [&localBodyId](std::pair<MInt, MInt> x) { return x.second == localBodyId; }),
                 a->end());

        for(const auto& [domain, bodies] : m_homeDomainRemoteBodies) {
          if(bodies.empty()) continue;

          for(const auto& bodyId : bodies) {
            // check if body was previously in m_transferBodies -> homeDomain has changed and
            // mapping needs to be updated accordingly
            if(domain != globalDomain && localBodyId == bodyId) {
              // Body was transfered and has a new homeDomain
              std::vector<MInt>* c = &m_homeDomainRemoteBodies[domain];
              c->erase(std::remove(c->begin(), c->end(), localBodyId), c->end());
              m_homeDomainRemoteBodies[globalDomain].push_back(localBodyId);
              break;
            }
          }
        }

        // Ensure the bodies in the mapping are in the received order
        // This is important for all communication buffers to be consistent
        std::vector<MInt>& remoteBodies = m_homeDomainRemoteBodies[globalDomain];
        auto actualPos = std::find(remoteBodies.begin(), remoteBodies.end(), localBodyId);
        auto wantedPos = remoteBodies.begin() + mapIdx;
        std::iter_swap(actualPos, wantedPos);
        mapIdx++;
        m_debugFileStream << "SWAP SWAP";
        for(auto&& r : remoteBodies) {
          m_debugFileStream << r << ",";
        }
        m_debugFileStream << std::endl;

      } else {
        // we want to maintain the structure in the collector (localBodies | remoteBodies | dummyBodies)
        // therefore all relevant Mappings have to be updated accordingly
        // if there are no dummyBodies in collector we can simply append to collector
        if(m_associatedDummyBodies.empty()) {
          m_bodies.append();
          m_globalBodyIds.push_back(remoteBodyId);
          m_noRemoteBodies++;
          m_homeDomainRemoteBodies[globalDomain].push_back(m_globalBodyIds.size() - 1);
          mapIdx++;
        } else {
          // if there are already dummy Bodies in the collector
          // a new Remote Body is inserted in the collector in front of all dummy bodies
          MInt newIdx = noConnectingBodies();
          m_bodies.insert(newIdx);
          m_globalBodyIds.insert(m_globalBodyIds.begin() + newIdx, remoteBodyId);
          m_noRemoteBodies++;
          m_homeDomainRemoteBodies[globalDomain].push_back(newIdx);
          mapIdx++;
          for(auto& [bodyId, dummyId] : m_associatedDummyBodies) {
            dummyId++;
          }
        }
      }

      const MInt newLocalBodyId = getLocalBodyId(remoteBodyId);

      for(MInt d = 0; d < nDim; d++) {
        a_bodyCenter(newLocalBodyId, d) = rcvBufferAll[n][b * bufSizePerBody + paraCount + d];
      }
      paraCount += nDim;

      for(MInt r = 0; r < nRot; r++) {
        a_angularAccelerationBody(newLocalBodyId, r) = rcvBufferAll[n][b * bufSizePerBody + paraCount + r];
        a_angularVelocityBodyT1B2(newLocalBodyId, r) = rcvBufferAll[n][b * bufSizePerBody + paraCount + nRot + r];
      }
      paraCount += 2 * nRot;

      for(MInt q = 0; q < nQuat; q++) {
        a_bodyQuaternionT1(newLocalBodyId, q) = rcvBufferAll[n][b * bufSizePerBody + paraCount + q];
        a_bodyQuaternionT1B2(newLocalBodyId, q) = rcvBufferAll[n][b * bufSizePerBody + paraCount + nQuat + q];
      }
      paraCount += 2 * nQuat;

      a_bodyRadius(newLocalBodyId) = rcvBufferAll[n][b * bufSizePerBody + paraCount];
      paraCount++;
      if(a_bodyType() != Shapes::Sphere) {
        for(MInt d = 0; d < nDim; d++) {
          a_bodyRadii(newLocalBodyId, d) = rcvBufferAll[n][b * bufSizePerBody + paraCount + d];
        }
        paraCount += nDim;
      } else {
        std::fill_n(&a_bodyRadii(newLocalBodyId, 0), nDim, a_bodyRadius(newLocalBodyId));
      }

      a_bodyDensityRatio(newLocalBodyId) = rcvBufferAll[n][b * bufSizePerBody + paraCount];
      paraCount++;
      a_bodyTemperature(newLocalBodyId) = rcvBufferAll[n][b * bufSizePerBody + paraCount];
      paraCount++;
      a_bodyHeatFlux(newLocalBodyId) = 0.0;

      // if domain receives a new remote Body, that has to be self-mapped on remoteDomain, create dummyBody
      const MInt remoteBodyIntersectsWithDomainId =
          grid().raw().intersectingWithHaloCells(&a_bodyCenter(newLocalBodyId, 0), m_infoDiameter / 2.0, 0, true);
      if(m_periodicBC && !hasAssociatedDummyBody(newLocalBodyId) && remoteBodyIntersectsWithDomainId != -1
         && grid().raw().a_hasProperty(remoteBodyIntersectsWithDomainId, Cell::IsPeriodic)) {
        createDummyBody(newLocalBodyId);

        std::array<MFloat, nDim> center{};
        for(MInt d = 0; d < nDim; d++) {
          center[n] = a_bodyCenter(newLocalBodyId, n);
        }

        for(MInt d = 0; d < nDim; d++) {
          m_periodicShift[domainId()][newLocalBodyId][d] = calculatePeriodicShift(remoteBodyIntersectsWithDomainId, d);
        }
        // initialize dummyBody kinematic data with the data of the remote Body
        initRemoteDummyBodyData(newLocalBodyId);
      }
    }
  }

#ifdef RB_DEBUG
  if(!remoteBodyRemoveList.empty()) {
    removeListS = "[";
    for(const auto& body : remoteBodyRemoveList) {
      removeListS += "[" + std::to_string(body.first) + "," + std::to_string(body.second) + "]";
    }
    removeListS += "] ";
    m_debugFileStream << "After: RemoveList for domain " << domainId() << " " << removeListS << std::endl;
  }
#endif

  deleteIrrelevantBodies(remoteBodyRemoveList);
}

/**
 * \brief  initalize kinematic data of dummyBody that is associated with a remoteBody
 *
 * \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 * \param[in]  bodyId , localBodyId of the remoteBody
 */
template <MInt nDim>
void RigidBodies<nDim>::initRemoteDummyBodyData(const MInt bodyId) {
  TRACE();

  const MInt dummyId = m_associatedDummyBodies[bodyId];
  if(isRemoteBody(bodyId)) {
    // Exchange kinematic Data for remote self-mapped Bodies
    for(MInt n = 0; n < nDim; n++) {
      a_bodyCenter(dummyId, n) = a_bodyCenter(bodyId, n) + m_periodicShift[domainId()][bodyId][n];
      a_bodyVelocity(dummyId, n) = a_bodyVelocity(bodyId, n);
      a_bodyAcceleration(dummyId, n) = a_bodyAcceleration(bodyId, n);
      a_bodyForce(dummyId, n) = -a_bodyForce(bodyId, n);
      a_bodyTorque(dummyId, n) = a_bodyTorque(bodyId, n);
    }
    for(MInt r = 0; r < nRot; r++) {
      a_angularVelocity(dummyId, r) = a_angularVelocity(bodyId, r);
      a_angularVelocityBody(dummyId, r) = a_angularVelocityBody(bodyId, r);
    }
    for(MInt q = 0; q < nQuat; q++) {
      a_bodyQuaternionT1(dummyId, q) = a_bodyQuaternionT1(bodyId, q);
    }
  } else {
#ifdef RB_DEBUG
    m_debugFileStream << "Something went wrong, body with localId: " << bodyId
                      << " is not a remoteBody -> dummyBodyData cannot be initialized!" << std::endl;
#endif
    return;
  }
}

/** \brief  this function transforms a localBody to a remoteBody on the current domain,
 *          necessary if a body crosses a partition boundary to another partition
 *
 *  \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 *  \param[in] oldBodyId is the Id the body had when it was a localBody
 *  \param[in] domain the new homeDomain of the body
 */
template <MInt nDim>
MInt RigidBodies<nDim>::localToRemoteBody(const MInt collectorId, const MInt oldBodyId, const MInt domain) {
  TRACE();

  if(m_noLocalBodies <= 0) {
    mTerm(1, AT_, "Domain: " + std::to_string(domainId()) + " has no local bodies that can change to remote status!");
  }

  const MInt newIdx = m_noLocalBodies - 1;

  if(hasAssociatedDummyBody(oldBodyId)) {
    deleteDummyBody(oldBodyId);
  }

  // If more than one localBody in collector -> do swap, else not necessary to swap, just mark local body as remote
  if(noConnectingBodies() > 0) {
    m_bodies.swap(collectorId, newIdx);
    iter_swap(m_globalBodyIds.begin() + collectorId, m_globalBodyIds.begin() + newIdx);
  }

  m_noLocalBodies--;
  m_noRemoteBodies++;

  /* Delete from local Bodies, now a remote Body
  / Not every local Body has a corresponding remoteDomain, therefore, we need to check each body, already
  / contained in the mapping. A local body can potentielly be remote to multiple domains
  / -> iterate over all domains in mapping
  */
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt globDomain = neighborDomain(n);
    std::vector<MInt>* bodies = &m_remoteDomainLocalBodies[globDomain];

    // Check if the swapped body is in the mapping
    auto swappedIt = std::find(bodies->begin(), bodies->end(), newIdx);
    if(swappedIt != bodies->end()) {
      // Swapped body is in the mapping, delete its old entry and keep the other for the swapped Body
      bodies->erase(std::remove(bodies->begin(), bodies->end(), newIdx), bodies->end());
    } else {
      // Swapped body is NOT in the mapping, just delete the entry of the old local body
      bodies->erase(std::remove(bodies->begin(), bodies->end(), collectorId), bodies->end());
    }

    if(hasAssociatedDummyBody(newIdx)) {
      m_associatedDummyBodies.erase(newIdx);
      m_associatedDummyBodies[newIdx] = collectorId;
    }
  }

  m_homeDomainRemoteBodies[domain].push_back(newIdx);

  return newIdx;
}


/** \brief Delete remote bodies from collector that are not relevant anymore
 *
 *  \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 *  \param[in] removeList list containing pair of homeDomain and localId of remoteBody
 */
template <MInt nDim>
void RigidBodies<nDim>::deleteIrrelevantBodies(std::list<std::pair<MInt, MInt>>& removeList) {
  TRACE();

  // Sort in descending order
  removeList.sort([](const std::pair<MInt, MInt>& a, const std::pair<MInt, MInt>& b) { return a.second > b.second; });

  std::vector<MInt> removeListGlobIds;
  for(const auto& [homeDomain, body] : removeList) {
    removeListGlobIds.push_back(getGlobalBodyId(body));
  }

  MInt collectorShiftIdx = 0;
  for(const auto& [homeDomain, body] : removeList) {
    m_bodyWasDeleted = true;
#ifdef RB_DEBUG
    m_debugFileStream << "Deleting remotebody " << body << " on domain " << domainId() << std::endl;
#endif

    m_globalBodyIds.erase(
        std::remove(m_globalBodyIds.begin(), m_globalBodyIds.end(), removeListGlobIds[collectorShiftIdx]),
        m_globalBodyIds.end());

    m_bodies.removeAndShift(body);
    m_noRemoteBodies--;

    m_homeDomainRemoteBodies[homeDomain].erase(
        std::remove(m_homeDomainRemoteBodies[homeDomain].begin(), m_homeDomainRemoteBodies[homeDomain].end(), body),
        m_homeDomainRemoteBodies[homeDomain].end());

    collectorShiftIdx++;
    // In case a remoteBody has also an associated dummy Body
    auto aDBCopy = m_associatedDummyBodies;
    for(auto& [bodyId, dummyId] : aDBCopy) {
      if(bodyId == body) {
        deleteDummyBody(bodyId, collectorShiftIdx);
      }
    }
  }

#ifdef RB_DEBUG
  m_debugFileStream << "after deletion: " << std::endl;
  printBodyDomainConnections();
#endif

  // Adjust localIds in mapping accordingly
  for(std::list<std::pair<MInt, MInt>>::const_iterator it = removeList.cbegin(); it != removeList.cend(); it++) {
    for(auto& [domain, bodies] : m_homeDomainRemoteBodies) {
      // adjust localIds
      for(MInt& body : bodies) {
#ifdef RB_DEBUG
        m_debugFileStream << "body " << body << " it " << it->second << std::endl;
#endif
        if(body > it->second) {
          body--;
        }
      }
    }
    for(auto& [bodyId, dummyId] : m_associatedDummyBodies) {
#ifdef RB_DEBUG
      m_debugFileStream << "body " << bodyId << " with Dummybody " << dummyId << " it " << it->second << std::endl;
#endif
      if(dummyId > it->second) {
        dummyId--;
      }
    }
  }

#ifdef RB_DEBUG
  m_debugFileStream << "after deletion and update of mappings: " << std::endl;
  printBodyDomainConnections();
#endif
}

/** \brief Delete dummy bodies from collector that are not relevant anymore
 *
 *  \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 *  \param[in] bodyId corresponds to the localId of the body of which the associated dummy body has to be deleted
 *  \param[in] collectorShift default=0, if this function is called mutiple times to delete multiple bodies
 *             we generate a certain Shift inside the collector
 */
template <MInt nDim>
void RigidBodies<nDim>::deleteDummyBody(const MInt bodyId, const MInt collectorShift) {
  TRACE();

  const MInt dummyBodyId = m_associatedDummyBodies[bodyId] - collectorShift;

#ifdef RB_DEBUG
  m_debugFileStream << "Deleting Dummy Body " << m_associatedDummyBodies[bodyId] << " with corresponding bodyId "
                    << bodyId << std::endl;
#endif

  m_associatedDummyBodies.erase(bodyId);
  m_globalBodyIds.erase(std::remove(m_globalBodyIds.begin(), m_globalBodyIds.end(), dummyBodyId),
                        m_globalBodyIds.end());

  m_bodies.removeAndShift(dummyBodyId);
  m_noDummyBodies--;

  m_bodyWasDeleted = true;

#ifdef RB_DEBUG
  printBodyDomainConnections();
#endif
}

template <MInt nDim>
void RigidBodies<nDim>::exchangeNeighborConnectionInfo() {
  TRACE();
  // for every remote body:
  // check every neighbordomain halo cells for intersection
  // contains globalIds
  std::vector<std::vector<MInt>> bodyIdsForRemoteDomain;
  std::vector<std::vector<MInt>> homeDomainIdsForRemoteDomain;

  // contains localIds
  std::vector<std::vector<MInt>> bodyIdsForHomeDomain;
  std::vector<std::vector<MInt>> remoteDomainIdsForHomeDomain;
  std::vector<std::vector<std::array<MFloat, nDim>>> shiftForHomeDomain;

  // body Ids + homeDomain Ids to send to remote domains
  bodyIdsForRemoteDomain.clear();
  bodyIdsForRemoteDomain.resize(noNeighborDomains());

  homeDomainIdsForRemoteDomain.clear();
  homeDomainIdsForRemoteDomain.resize(noNeighborDomains());

  // body Ids + remoteDomain Ids + shift to send to home domains
  bodyIdsForHomeDomain.clear();
  bodyIdsForHomeDomain.resize(noHomeDomains());

  remoteDomainIdsForHomeDomain.clear();
  remoteDomainIdsForHomeDomain.resize(noHomeDomains());

  shiftForHomeDomain.clear();
  shiftForHomeDomain.resize(noHomeDomains());

  MInt count = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    count = 0;
    for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
      for(const auto& body : bodies) {
        // check neighbor domain halo cells
        std::array<MFloat, nDim> center{};
        for(MInt n = 0; n < nDim; n++) {
          center[n] = a_bodyCenter(body, n);
        }
        const MInt intersectingCellId =
            grid().raw().intersectingWithHaloCells(center.data(), m_infoDiameter / 2.0, i, true);
        if(intersectingCellId != -1) {
          // data for new remote domain
          bodyIdsForRemoteDomain[i].push_back(getGlobalBodyId(body));
          homeDomainIdsForRemoteDomain[i].push_back(homeDomain);

          // data for home domain
          const MInt remoteDomain = neighborDomain(i);
          // dont send body to its own home domain
          if(homeDomain == remoteDomain) {
            continue;
          }
          bodyIdsForHomeDomain[count].push_back(getGlobalBodyId(body));
          remoteDomainIdsForHomeDomain[count].push_back(remoteDomain);
          // shift
          std::array<MFloat, nDim> shiftVector{};
          for(MInt n = 0; n < nDim; n++) {
            const MFloat shift = calculatePeriodicShift(intersectingCellId, n);
            shiftVector[n] = shift;
          }
          shiftForHomeDomain[count].push_back(shiftVector);
        }
      }
      count++;
    }
  }

  std::map<MInt, std::vector<MInt>> tmpRemoteMap;

  exchangeNeighborNeighborRemote(bodyIdsForRemoteDomain, homeDomainIdsForRemoteDomain, tmpRemoteMap);
  exchangeNeighborNeighborHome(bodyIdsForHomeDomain, remoteDomainIdsForHomeDomain, shiftForHomeDomain);

  // transfer tmp maps to member mappings
  for(auto [homeDomainId, remoteBodyIds] : tmpRemoteMap) {
    for(auto remoteBodyId : remoteBodyIds) {
      if(homeDomainId == domainId()) continue;
      const std::vector<MInt>* bodyIds = &(m_homeDomainRemoteBodies[homeDomainId]);
      // check if remote body is not in mapping yet
      if(!std::binary_search(bodyIds->begin(), bodyIds->end(), remoteBodyId)) {
        std::cout << "In Timestep " << globalTimeStep << "indirect Neighborstuff" << std::endl;
        // PseudoLocal, as just remote body for domain, but body is contained in collector
        // const MInt newPseudoLocalBodyId = getLocalBodyId(remoteBodyId);
        m_homeDomainRemoteBodies[homeDomainId].push_back(remoteBodyId);
      }
    }
  }
}

template <MInt nDim>
void RigidBodies<nDim>::exchangeNeighborNeighborRemote(std::vector<std::vector<MInt>>& bodyIdsForRemoteDomain,
                                                       std::vector<std::vector<MInt>>& homeDomainIdsForRemoteDomain,
                                                       std::map<MInt, std::vector<MInt>>& tmpRemoteMap) {
  TRACE();
  std::vector<MPI_Request> mpi_send_req(noNeighborDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noNeighborDomains(), MPI_REQUEST_NULL);

  // body Ids, home domain
  const MInt bufSizePerBody = 1 + 1;
  std::vector<MInt> noToSend(noNeighborDomains());
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    noToSend[n] = bodyIdsForRemoteDomain[n].size();
  }
  // pack send buffer
  std::vector<std::vector<MInt>> sendBufferAll(noNeighborDomains());
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    sendBufferAll[n].resize(noToSend[n] * bufSizePerBody);
    for(MInt b = 0; b < noToSend[n]; b++) {
      sendBufferAll[n][b * bufSizePerBody] = bodyIdsForRemoteDomain[n][b];
      sendBufferAll[n][b * bufSizePerBody + 1] = homeDomainIdsForRemoteDomain[n][b];
    }
  }
  // exchange buffer lengths
  std::vector<MInt> noToRecv(noNeighborDomains());
  exchangeBufferLengthsNeighbor(noToRecv, bodyIdsForRemoteDomain);

  // create receive buffer
  std::vector<std::vector<MInt>> recvBufferAll(noNeighborDomains());
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    recvBufferAll[n].resize(noToRecv[n] * bufSizePerBody);
  }

  // recv
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Irecv(recvBufferAll[n].data(), recvBufferAll[n].size(), maia::type_traits<MInt>::mpiType(), neighborDomain(n),
              2, mpiComm(), &mpi_recv_req[n], AT_, "recv neighbor-neighbor exchange");
  }

  // send
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MPI_Isend(sendBufferAll[n].data(), sendBufferAll[n].size(), maia::type_traits<MInt>::mpiType(), neighborDomain(n),
              2, mpiComm(), &mpi_send_req[n], AT_, "send neighbor-neighbor exchange");
  }

  // wait for communication to finish
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
    MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
  }

  // unpack recv Buffer and fill remoteBodyHomeDomains
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt b = 0; b < noToRecv[n]; b++) {
      const MInt remoteBodyId = getLocalBodyId(recvBufferAll[n][b * bufSizePerBody]);
      const MInt homeDomain = recvBufferAll[n][b * bufSizePerBody + 1];
      const std::vector<MInt>* bodyIds = &(tmpRemoteMap[homeDomain]);
      // check if remote body is not in mappin yet
      if(!std::binary_search(bodyIds->begin(), bodyIds->end(), remoteBodyId)) {
        tmpRemoteMap[homeDomain].push_back(remoteBodyId);
      }
    }
  }
}

template <MInt nDim>
void RigidBodies<nDim>::exchangeNeighborNeighborHome(
    std::vector<std::vector<MInt>>& bodyIdsForHomeDomain,
    std::vector<std::vector<MInt>>& remoteDomainIdsForHomeDomain,
    std::vector<std::vector<std::array<MFloat, nDim>>>& shiftForHomeDomain) {
  TRACE();
  std::vector<MPI_Request> mpi_send_req(noHomeDomains(), MPI_REQUEST_NULL);
  std::vector<MPI_Request> mpi_recv_req(noRemoteDomains(), MPI_REQUEST_NULL);

  // exchange buffer lengths to home domains
  std::vector<MInt> noToRecv(noRemoteDomains());

  // body Id + remote domain id + shift
  const MInt bufSizePerBody = 1 + 1 + nDim;
  std::vector<MInt> noToSend(noHomeDomains());
  for(MInt n = 0; n < noHomeDomains(); n++) {
    noToSend[n] = bodyIdsForHomeDomain[n].size();
  }

  exchangeBufferLengthsRemoteToHome(noToRecv, noToSend);

  // create receive buffer
  std::vector<std::vector<MFloat>> recvBufferAll(noRemoteDomains());
  for(MInt n = 0; n < noRemoteDomains(); n++) {
    recvBufferAll[n].resize(noToRecv[n] * bufSizePerBody);
  }

  // pack send buffer
  std::vector<std::vector<MFloat>> sendBufferAll(noHomeDomains());
  for(MInt n = 0; n < noHomeDomains(); n++) {
    sendBufferAll[n].resize(noToSend[n] * bufSizePerBody);
    for(MInt b = 0; b < noToSend[n]; b++) {
      sendBufferAll[n][b * bufSizePerBody] = bodyIdsForHomeDomain[n][b];
      sendBufferAll[n][b * bufSizePerBody + 1] = remoteDomainIdsForHomeDomain[n][b];
      for(MInt i = 0; i < nDim; i++) {
        sendBufferAll[n][b * bufSizePerBody + 2 + i] = shiftForHomeDomain[n][b][i];
      }
    }
  }

  // recv
  MInt count = 0;
  for(const auto& [remoteDomain, bodies] : m_remoteDomainLocalBodies) {
    if(noToRecv[count]) {
      MPI_Irecv(recvBufferAll[count].data(), recvBufferAll[count].size(), maia::type_traits<MFloat>::mpiType(),
                remoteDomain, 2, mpiComm(), &mpi_recv_req[count], AT_, "recv neighbor-neighbor exchange");
    }
    count++;
  }

  // send
  count = 0;
  for(const auto& [homeDomain, bodies] : m_homeDomainRemoteBodies) {
    if(noToSend[count]) {
      MPI_Isend(sendBufferAll[count].data(), sendBufferAll[count].size(), maia::type_traits<MFloat>::mpiType(),
                homeDomain, 2, mpiComm(), &mpi_send_req[count], AT_, "send neighbor-neighbor exchange");
    }
    count++;
  }

  // wait for communication to finish
  for(MInt i = 0; i < noHomeDomains(); i++) {
    if(noToSend[i]) {
      MPI_Wait(&mpi_send_req[i], MPI_STATUS_IGNORE, AT_);
    }
  }

  for(MInt i = 0; i < noRemoteDomains(); i++) {
    if(noToRecv[i]) {
      MPI_Wait(&mpi_recv_req[i], MPI_STATUS_IGNORE, AT_);
    }
  }

  // unpack buffer
  const MInt currNoRemoteDomains = noRemoteDomains();
  for(MInt n = 0; n < currNoRemoteDomains; n++) {
    for(MInt b = 0; b < noToRecv[n]; b++) {
      const MInt localBodyId = getLocalBodyId((MInt)recvBufferAll[n][b * bufSizePerBody]);
      const MInt remoteDomainId = (MInt)recvBufferAll[n][b * bufSizePerBody + 1];

      std::array<MFloat, nDim> shift{};
      for(MInt i = 0; i < nDim; i++) {
        shift[i] = recvBufferAll[n][b * bufSizePerBody + 2 + i];
      }

      if(m_remoteDomainLocalBodies.find(remoteDomainId) != m_remoteDomainLocalBodies.end()) {
        std::vector<MInt>* bodyIds = &(m_remoteDomainLocalBodies[remoteDomainId]);
        if(std::binary_search(bodyIds->begin(), bodyIds->end(), localBodyId)) {
          continue;
        }
      }
      m_remoteDomainLocalBodies[remoteDomainId].push_back(localBodyId);
      m_periodicShift[remoteDomainId][localBodyId] = shift;
    }
  }
}

/**
 * \brief debugging: print all necessary mappings for each partition in specific file
 *
 * \author Johannes Grafen <johannes.grafen@rwth-aachen.de>
 *
 */
template <MInt nDim>
void RigidBodies<nDim>::printBodyDomainConnections() {
  TRACE();

  std::vector<MString> sBodies;
  MInt idx = 0;

  sBodies.emplace_back("");
  if(m_noLocalBodies != 0) {
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      sBodies[idx] += std::to_string(bodyId) + " ";
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(m_noRemoteBodies > 0) {
    for(MInt bodyId = m_noLocalBodies; bodyId < noConnectingBodies(); bodyId++) {
      sBodies[idx] += std::to_string(bodyId) + " ";
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(m_noDummyBodies > 0) {
    for(MInt bodyId = noConnectingBodies(); bodyId < noCollectorBodies(); bodyId++) {
      sBodies[idx] += std::to_string(bodyId) + " ";
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(!m_transferBodies.empty()) {
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      if(!m_transferBodies[i].empty()) {
        MString sTmp;
        for(MInt body : m_transferBodies[i]) {
          sTmp += std::to_string(body) + " ";
        }
        sBodies[idx] += std::to_string(neighborDomain(i)) + ": [" + sTmp + " ] | ";
      }
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(!m_globalBodyIds.empty()) {
    for(MUint localBodyId = 0; localBodyId < m_globalBodyIds.size(); localBodyId++) {
      sBodies[idx] +=
          "| localId: " + std::to_string(localBodyId) + " globalId: " + std::to_string(m_globalBodyIds[localBodyId]);
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(!m_remoteDomainLocalBodies.empty()) {
    for(const auto& [domain, bodies] : m_remoteDomainLocalBodies) {
      for(const auto& body : bodies) {
        sBodies[idx] += std::string(" ") + "[" + std::to_string(domain) + ", " + std::to_string(body) + "]";
      }
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(!m_homeDomainRemoteBodies.empty()) {
    for(const auto& [domain, bodies] : m_homeDomainRemoteBodies) {
      for(const auto& body : bodies) {
        sBodies[idx] += std::string(" ") + "[" + std::to_string(domain) + ", " + std::to_string(body) + "]";
      }
    }
  }
  idx++;

  sBodies.emplace_back("");
  if(!m_associatedDummyBodies.empty()) {
    for(const auto& [body, dummy] : m_associatedDummyBodies) {
      sBodies[idx] += std::string(" ") + "[" + std::to_string(body) + ", " + std::to_string(dummy) + "]";
    }
  }
  idx++;

  m_debugFileStream << "--------------BODY INFO for domain: " << domainId() << "---------------------" << std::endl
                    << " contains bodyIds [ " << sBodies[0] << " ]" << std::endl
                    << " contains remoteBodyIds [ " << sBodies[1] << " ]" << std::endl
                    << " contains dummyBodyIds [ " << sBodies[2] << " ]" << std::endl
                    << " has transferBodies for neighbourdomain " << sBodies[3] << std::endl
                    << " has following pairings: " << sBodies[4] << std::endl
                    << " remote domain - local body pairings: " << sBodies[5] << std::endl
                    << " home domain - remote body pairings: " << sBodies[6] << std::endl
                    << " body - associated dummy pairings: " << sBodies[7] << std::endl
                    << " collector has capacity: " << m_bodies.capacity() << " and size " << m_bodies.size()
                    << " MaxNoBodies is " << m_maxNoBodies << std::endl;
}

/**
 * \brief Perform the prediction for the rotational state of a single body in 3D
 *
 * For reference, see L.J.H. Seelen et al.,
 * Improved bodyQuaternion-based integration scheme for rigid body motion
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyId Body id for which the prediction is performed
 */
template <>
inline void RigidBodies<3>::predictRotation(const MInt bodyId) {
  const MFloat dt = m_timestep;

  MFloat angularVelocityBodyT3B4[nRot]{};
  // Advance angular velocity from time n+1/2 to time n+3/4
  for(MInt n = 0; n < nRot; n++) {
    angularVelocityBodyT3B4[n] =
        a_angularVelocityBodyT1B2(bodyId, n) + 0.25 * dt * a_angularAccelerationBody(bodyId, n);
  }
  // Transform angular velocity to world frame
  MFloat angularVelocityT3B4[nRot]{};
  transformToWorldFrame(&a_bodyQuaternionT1B2(bodyId, 0), angularVelocityBodyT3B4, angularVelocityT3B4);

  // Advance bodyQuaternion from time n+1/2 to time n+1
  advanceQuaternion(angularVelocityT3B4, dt / 2, &a_bodyQuaternionT1B2(bodyId, 0), &a_bodyQuaternionT1(bodyId, 0));

  // Advance angular velocity from time n+1/2 to time n+3/4
  for(MInt n = 0; n < nRot; n++) {
    a_angularVelocityBodyT1(bodyId, n) =
        a_angularVelocityBodyT1B2(bodyId, n) + 0.5 * dt * a_angularAccelerationBody(bodyId, n);
  }
  transformToWorldFrame(&a_bodyQuaternionT1(bodyId, 0), &a_angularVelocityBodyT1(bodyId, 0),
                        &a_angularVelocityT1(bodyId, 0));
}

/**
 * \brief Empty function for the 2D rotation prediction
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param MInt Body id for wich the prediction is performed
 */
template <>
inline void RigidBodies<2>::predictRotation(const MInt /*bodyId*/) {
  std::cerr << "predictRotation: No implementation for 2D!" << std::endl;
}

/**
 * \brief Correct the state of all bodies using external forces/fluxes/torques
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::correctorStep() {
  TRACE();

  const MFloat dt = m_timestep;
  if(m_translation) {
    for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
      // Translational values
      for(MInt n = 0; n < nDim; n++) {
        // Correct acceleration
        a_bodyAcceleration(bodyId, n) =
            a_bodyForce(bodyId, n) / a_bodyMass(bodyId)
            + gravity[n] * (1.0 - 1.0 / a_bodyDensityRatio(bodyId)); // Subtract bouyancy force

        // Correct velocities
        a_bodyVelocity(bodyId, n) = a_bodyVelocityOld(bodyId, n)
                                    + 0.5 * dt * (a_bodyAcceleration(bodyId, n) + a_bodyAccelerationOld(bodyId, n));

        // Correct position
        a_bodyCenter(bodyId, n) =
            a_bodyCenterOld(bodyId, n) + dt * a_bodyVelocityOld(bodyId, n)
            + 0.25 * POW2(dt) * (a_bodyAcceleration(bodyId, n) + a_bodyAccelerationOld(bodyId, n));
      }

      if(m_rotation) {
        correctRotation(bodyId);
      }

      // Scalar values
      // Correct temperature
      a_bodyTemperature(bodyId) = F1B2
                                  * (a_bodyTemperature(bodyId) + a_bodyTemperatureOld(bodyId)
                                     + dt * a_bodyHeatFlux(bodyId) / (m_bodyHeatCapacity * a_bodyMass(bodyId)));
    }

  } else {
    if(m_rotation) {
      for(MInt bodyId = 0; bodyId < m_noLocalBodies; bodyId++) {
        correctRotation(bodyId);
      }
    }
  }
}

/**
 * \brief Perform the correction for the rotational state for a single body in 3D
 *
 * For reference, see L.J.H. Seelen et al.,
 * Improved bodyQuaternion-based integration scheme for rigid body motion
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyId Body id for which the correction is performed
 */
template <>
inline void RigidBodies<3>::correctRotation(const MInt bodyId) {
  const MFloat dt = m_timestep;

  MFloat torqueBodyT1[nRot];
  // Transform torque to the body frame
  transformToBodyFrame(&a_bodyQuaternionT1(bodyId, 0), &a_torqueT1(bodyId, 0), torqueBodyT1);

  // Calculate angular momentum in the body frame
  MFloat angularMomentumBody[nRot]{};
  for(MInt n = 0; n < nRot; n++) {
    angularMomentumBody[n] = a_bodyInertia(bodyId, n) * a_angularVelocityBodyT1(bodyId, n);
  }

  // Calculate convective part of the angular acceleration
  MFloat convectiveTerm[nRot]{};
  maia::math::cross(&a_angularVelocityBodyT1(bodyId, 0), &angularMomentumBody[0], &convectiveTerm[0]);

  // Calculate angular acceleration
  MFloat angularAccelerationBodyT1[nRot]{};
  for(MInt n = 0; n < nRot; n++) {
    angularAccelerationBodyT1[n] = 1 / a_bodyInertia(bodyId, n) * (torqueBodyT1[n] - convectiveTerm[n]);
  }
  transformToWorldFrame(&a_bodyQuaternionT1(bodyId, 0), angularAccelerationBodyT1, &a_angularAccelerationT1(bodyId, 0));

  // Advance angular velocity from time t+1/2 to time t+3/4
  MFloat angularVelocityBodyT3B2[nRot]{};
  for(MInt n = 0; n < nRot; n++) {
    angularVelocityBodyT3B2[n] = a_angularVelocityBodyT1B2(bodyId, n) + dt * angularAccelerationBodyT1[n];
  }

  // Advance bodyQuaternion from time t+1/2 to time t+3/2
  MFloat bodyQuaternionT3B2[nQuat]{};
  advanceQuaternion(&a_angularVelocityT1(bodyId, 0), dt, &a_bodyQuaternionT1B2(bodyId, 0), bodyQuaternionT3B2);

  // Transform angular velocity to the world frame
  MFloat angularVelocityT3B2[nRot]{};
  transformToWorldFrame(bodyQuaternionT3B2, angularVelocityBodyT3B2, angularVelocityT3B2);

  // Advance time step
  constexpr MLong size_rot = nRot * sizeof(MFloat);
  constexpr MLong size_quat = nQuat * sizeof(MFloat);

  memcpy(&a_angularAccelerationBody(bodyId, 0), angularAccelerationBodyT1, size_rot);
  memcpy(&a_angularVelocityBodyT1B2(bodyId, 0), angularVelocityBodyT3B2, size_rot);
  memcpy(&a_angularVelocityT1B2(bodyId, 0), angularVelocityT3B2, size_rot);
  memcpy(&a_bodyQuaternionT1B2(bodyId, 0), bodyQuaternionT3B2, size_quat);
}

/**
 * \brief Empty function for the 2D rotation correction
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param MInt Body id for wich the correction is performed
 */
template <>
inline void RigidBodies<2>::correctRotation(const MInt /*bodyId*/) {
  std::cerr << "correctRotation: No implementation for 2D!" << std::endl;
}

/**
 * \brief Transform a vector form the body frame to the world frame
 *
 * For reference, see L.J.H. Seelen et al.,
 * Improved bodyQuaternion-based integration scheme for rigid body motion
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyQuaternion[in]   bodyQuaternion describing the rotational state of the body frame
 * \param vectorBody[in]   Vector defined in the body frame
 * \param vectorWorld[out] Vector defined in the world frame
 */
template <MInt nDim>
inline void RigidBodies<nDim>::transformToWorldFrame(const MFloat* const bodyQuaternion,
                                                     const MFloat* const vectorBody,
                                                     MFloat* const vectorWorld) {
  const MFloat uv = std::inner_product(bodyQuaternion, &bodyQuaternion[3], vectorBody, 0.0);
  const MFloat uu = std::inner_product(bodyQuaternion, &bodyQuaternion[3], bodyQuaternion, 0.0);
  const MFloat ss = bodyQuaternion[3] * bodyQuaternion[3];
  MFloat ucv[nRot]{};
  IF_CONSTEXPR(nDim == 3) { maia::math::cross(bodyQuaternion, vectorBody, &ucv[0]); }
  else {
    return;
  }
  for(MInt n = 0; n < nRot; n++) {
    vectorWorld[n] = 2.0 * uv * bodyQuaternion[n] + (ss - uu) * vectorBody[n] + 2.0 * bodyQuaternion[3] * -(ucv[n]);
  }
}

/**
 * \brief Transform a vector form the world frame to the body frame
 *
 * For reference, see L.J.H. Seelen et al.,
 * Improved bodyQuaternion-based integration scheme for rigid body motion
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyQuaternion[in]   bodyQuaternion describing the rotational state of the body frame
 * \param vectorWorld[in]   Vector defined in the world frame
 * \param vectorBody[out] Vector defined in the body frame
 */
template <MInt nDim>
inline void RigidBodies<nDim>::transformToBodyFrame(const MFloat* const bodyQuaternion,
                                                    const MFloat* const vectorWorld,
                                                    MFloat* const vectorBody) {
  const MFloat uv = std::inner_product(bodyQuaternion, &bodyQuaternion[3], vectorWorld, 0.0);
  const MFloat uu = std::inner_product(bodyQuaternion, &bodyQuaternion[3], bodyQuaternion, 0.0);
  const MFloat ss = bodyQuaternion[3] * bodyQuaternion[3];
  MFloat ucv[nRot]{};
  IF_CONSTEXPR(nDim == 3) { maia::math::cross(bodyQuaternion, vectorWorld, &ucv[0]); }
  else {
    return;
  }
  for(MInt n = 0; n < nRot; n++) {
    vectorBody[n] = 2.0 * uv * bodyQuaternion[n] + (ss - uu) * vectorWorld[n] + 2.0 * bodyQuaternion[3] * ucv[n];
  }
}

/**
 * \brief Transform a vector form the world frame to the body frame
 *
 * This version modifies the vector in place.
 *
 * For reference, see L.J.H. Seelen et al.,
 * Improved bodyQuaternion-based integration scheme for rigid body motion
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyQuaternion[in]     bodyQuaternion describing the rotational state of the body frame
 * \param vectorBody[in,out] Vector to be transformed from the world to body frame
 */
template <MInt nDim>
inline void RigidBodies<nDim>::transformToBodyFrame(const MFloat* const bodyQuaternion, MFloat* const vec) {
  const MFloat vectorWorld[3]{vec[0], vec[1], vec[2]};
  transformToBodyFrame(bodyQuaternion, vectorWorld, vec);
}

/**
 * \brief Advance the rotational state by one timestep for a given angular velocity
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  angularVelocity Angular velocity
 * \param[in]  dt              Time step
 * \param[in]  qIn             bodyQuaternion before the rotation
 * \param[out] qOut            bodyQuaternion after the rotation
 */
template <MInt nDim>
inline void RigidBodies<nDim>::advanceQuaternion(const MFloat* const angularVelocity,
                                                 const MFloat dt,
                                                 const MFloat* const qIn,
                                                 MFloat* const qOut) {
  const MFloat angularVelMag = sqrt(std::inner_product(angularVelocity, &angularVelocity[nRot], angularVelocity, 0.0));

  MFloat qIncrement[4]{0.0, 0.0, 0.0, 1.0};

  if(angularVelMag > 0.0) {
    const MFloat angle = 0.5 * angularVelMag * dt;
    const MFloat tmp = sin(angle) / angularVelMag;
    qIncrement[0] = tmp * angularVelocity[0];
    qIncrement[1] = tmp * angularVelocity[1];
    qIncrement[2] = tmp * angularVelocity[2];
    qIncrement[3] = cos(angle);
  }

  maia::math::quatMult(qIncrement, qIn, qOut);
}

/**
 * \brief Transform the quaternion to Euler angles
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  quaternion quaternion describing the rotational state
 * \param[out] angles     Euler angles describing the rotational state
 */
template <MInt nDim>
void RigidBodies<nDim>::quaternionToEuler(const MFloat* const quaternion, MFloat* const angles) {
  TRACE();
  const MFloat& x = quaternion[0];
  const MFloat& y = quaternion[1];
  const MFloat& z = quaternion[2];
  const MFloat& w = quaternion[3];

  MFloat& roll = angles[0];
  MFloat& pitch = angles[1];
  MFloat& yaw = angles[2];

  // roll (x-axis rotation)
  MFloat sinr_cosp = 2 * (w * x + y * z);
  MFloat cosr_cosp = 1 - 2 * (x * x + y * y);
  roll = std::atan2(sinr_cosp, cosr_cosp);

  // pitch (y-axis rotation)
  MFloat sinp = 2 * (w * y - z * x);
  if(std::abs(sinp) >= 1) {
    pitch = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
  } else {
    pitch = std::asin(sinp);
  }

  // yaw (z-axis rotation)
  MFloat siny_cosp = 2 * (w * z + x * y);
  MFloat cosy_cosp = 1 - 2 * (y * y + z * z);
  yaw = std::atan2(siny_cosp, cosy_cosp);
}

/**
 * \brief Analytic motion equations for forced motions
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  returnMode Define whether location(0), velocity(1) or acceleration(2) is returned
 * \param[out] bodyData Body data for a single body
 * \param[in]  bodyId Body id for which is motion is calculated
 * \param[in]  time Current time
 */
template <MInt nDim>
void RigidBodies<nDim>::computeBodyPropertiesForced(MInt returnMode, MFloat* bodyData, MInt bodyId, MFloat time) {
  TRACE();

  MFloat elapsedTime = time;
  MFloat angle = F0;
  MBool& first = m_static_computeBodyProperties_first;
  MFloat*& amplitude = m_static_computeBodyProperties_amplitude;
  MFloat*& freqFactor = m_static_computeBodyProperties_freqFactor;
  MFloat*& initialBodyCenter = m_static_computeBodyProperties_initialBodyCenter;
  MFloat*& Strouhal = m_static_computeBodyProperties_Strouhal;
  MFloat*& mu2 = m_static_computeBodyProperties_mu2;
  MFloat*& liftStartAngle1 = m_static_computeBodyProperties_liftStartAngle1;
  MFloat*& liftEndAngle1 = m_static_computeBodyProperties_liftEndAngle1;
  MFloat*& liftStartAngle2 = m_static_computeBodyProperties_liftStartAngle2;
  MFloat*& liftEndAngle2 = m_static_computeBodyProperties_liftEndAngle2;
  MFloat*& circleStartAngle = m_static_computeBodyProperties_circleStartAngle;
  MFloat*& normal = m_static_computeBodyProperties_normal;
  MInt*& bodyToFunction = m_static_computeBodyProperties_bodyToFunction;
  MFloat*& rotAngle = m_static_computeBodyProperties_rotAngle;


  // Could be moved to readProperties ... ~jv
  if(first) {
    // Allocate memory for body movement function data
    mAlloc(m_static_computeBodyProperties_amplitude, noEmbeddedBodies(), "m_static_computeBodyProperties_amplitude",
           AT_);
    m_static_computeBodyProperties_amplitude[0] = F0;
    mAlloc(m_static_computeBodyProperties_freqFactor, noEmbeddedBodies(), "m_static_computeBodyProperties_freqFactor",
           AT_);
    mAlloc(m_static_computeBodyProperties_initialBodyCenter,
           noEmbeddedBodies() * nDim,
           "m_static_computeBodyProperties_initialBodyCenter",
           AT_);
    mAlloc(m_static_computeBodyProperties_Strouhal, noEmbeddedBodies(), "m_static_computeBodyProperties_Strouhal", AT_);
    mAlloc(m_static_computeBodyProperties_mu2, noEmbeddedBodies(), "m_static_computeBodyProperties_mu2", AT_);
    mAlloc(m_static_computeBodyProperties_liftStartAngle1,
           noEmbeddedBodies(),
           "m_static_computeBodyProperties_liftStartAngle1",
           AT_);
    mAlloc(m_static_computeBodyProperties_liftEndAngle1, noEmbeddedBodies(),
           "m_static_computeBodyProperties_liftEndAngle1", AT_);
    mAlloc(m_static_computeBodyProperties_liftStartAngle2,
           noEmbeddedBodies(),
           "m_static_computeBodyProperties_liftStartAngle2",
           AT_);
    mAlloc(m_static_computeBodyProperties_liftEndAngle2, noEmbeddedBodies(),
           "m_static_computeBodyProperties_liftEndAngle2", AT_);
    mAlloc(m_static_computeBodyProperties_circleStartAngle,
           noEmbeddedBodies(),
           "m_static_computeBodyProperties_circleStartAngle",
           AT_);
    mAlloc(m_static_computeBodyProperties_normal, noEmbeddedBodies() * nDim, "m_static_computeBodyProperties_normal",
           AT_);
    mAlloc(m_static_computeBodyProperties_bodyToFunction,
           noEmbeddedBodies(),
           "m_static_computeBodyProperties_bodyToFunction",
           AT_);
    mAlloc(m_static_computeBodyProperties_rotAngle, noEmbeddedBodies(), "m_static_computeBodyProperties_rotAngle", AT_);

    // 1: set default values:
    for(MInt k = 0; k < noEmbeddedBodies(); k++) {
      Strouhal[k] = 0.2;
      amplitude[k] = 1.0;
      freqFactor[k] = 1.0;
      bodyToFunction[k] = 0;
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
    }

    // 2: read Properties

    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      amplitude[i] = Context::getSolverProperty<MFloat>("amplitudes", m_solverId, AT_, &amplitude[i], i);
    }
    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      freqFactor[i] = Context::getSolverProperty<MFloat>("freqFactors", m_solverId, AT_, &freqFactor[i], i);
    }

    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      bodyToFunction[i] =
          Context::getSolverProperty<MInt>("bodyMovementFunctions", m_solverId, AT_, &bodyToFunction[i], i);
    }

    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      Strouhal[i] = Context::getSolverProperty<MFloat>("Strouhal", m_solverId, AT_, &Strouhal[i], i);
    }

    /*! \page propertyPage1
\section initialBodyCenters
<code>MFloat LsCartesianSolver::initialBodyCenter</code>\n
      For each body, nDim Float values.
      With this property, one can move an STL around to its initial position.
      In initialInsidePoints is the "real" center of the stl files.
      This property, is the movement from the "real" center to the initial center used during
calculation.

Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      for(MInt j = 0; j < nDim; j++) {
        initialBodyCenter[i * nDim + j] = Context::getSolverProperty<MFloat>(
            "initialBodyCenters", m_solverId, AT_, &initialBodyCenter[i * nDim + j], i * nDim + j);
      }
    }

    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      for(MInt j = 0; j < nDim; j++) {
        normal[i * nDim + j] = Context::getSolverProperty<MFloat>("bodyMotionNormals", m_solverId, AT_,
                                                                  &normal[i * nDim + j], i * nDim + j);
      }
    }

    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      liftStartAngle1[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles1", m_solverId, AT_, &liftStartAngle1[i], i);
    }

    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      liftStartAngle2[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles2", m_solverId, AT_, &liftStartAngle2[i], i);
    }

    /*! \page propertyPage1
    \section liftEndAngles1
    <code>MFloat LsCartesianSolver::liftEndAngle1</code>\n
    default = <code>PI</code>\n \n
          For each body, sets the end angle of the translation.
    <ul>
    <li>Between 0 and 4.0 * PI </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      liftEndAngle1[i] = Context::getSolverProperty<MFloat>("liftEndAngles1", m_solverId, AT_, &liftEndAngle1[i], i);
    }

    /*! \page propertyPage1
    \section liftEndAngles2
    <code>MFloat LsCartesianSolver::liftEndAngle2</code>\n
    default = <code>4.0 * PI</code>\n \n
          For each body, sets the end angle of the translation.
    <ul>
    <li>Between 0 and 4.0 * PI </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      liftEndAngle2[i] = Context::getSolverProperty<MFloat>("liftEndAngles2", m_solverId, AT_, &liftEndAngle2[i], i);
    }

    /*! \page propertyPage1
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
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      circleStartAngle[i] =
          Context::getSolverProperty<MFloat>("circleStartAngles", m_solverId, AT_, &circleStartAngle[i], i);
    }

    /*! \page propertyPage1
    \section rotAngle
    <code> MFloat LsPar::rotAngle </code>  \n
    default = <code>"0.0"</code>\n \n
    Used for specific bodyToFunction settings in the Ls-Solver to rotate the primary direction of the
    body movement. <ul> <li>Any positive Float</li>
    </ul>
    Keywords: <i>LEVELSET EMBEDED BOUNDARY, MOVEMENT FUNCTIONS </i>
   */

    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      rotAngle[i] = Context::getSolverProperty<MFloat>("rotAngle", m_solverId, AT_, &rotAngle[i], i);
      rotAngle[i] *= -PI / 180;
    }

    // 3: compute relevant values:
    for(MInt k = 0; k < noEmbeddedBodies(); k++) {
      // when using mu  : has a dimension!
      // when using mu2 : dimensionless!
      mu2[k] = freqFactor[k] * Strouhal[k] * F2 * PI;
    }

    // if bodyMovementFunction is 6 or 7, adjust start and end angles:
    for(MInt i = 0; i < noEmbeddedBodies(); i++) {
      if(bodyToFunction[i] == 6 || bodyToFunction[i] == 7 || bodyToFunction[i] == 9) {
        liftStartAngle1[i] = liftStartAngle1[i] * 2 * PI;
        liftEndAngle1[i] = liftEndAngle1[i] * 2 * PI;
      }
    }

    first = false;
  }


  //--------------------------------

  switch(bodyToFunction[bodyId]) {
    case 0:
      // Dummy movement function
      for(MInt n = 0; n < nDim; n++) {
        bodyData[n] = F0;
      }
      break;

    case 8:
      // translational movement in normal-direction with constant velocity
      switch(returnMode) {
        case 1: // return body center
          for(MInt n = 0; n < nDim; n++) {
            bodyData[n] = amplitude[bodyId] * elapsedTime * normal[bodyId * nDim + n];
          }
          break;
        case 2: // return body velocity
          for(MInt n = 0; n < nDim; n++) {
            bodyData[n] = amplitude[bodyId] * normal[bodyId * nDim + n];
          }
          break;
        case 3: // return body acceleration
          for(MInt n = 0; n < nDim; n++) {
            bodyData[n] = F0;
          }
          break;

        default:
          for(MInt n = 0; n < nDim; n++) {
            bodyData[n] = F0;
          }
          break;
      }
      break;

    case 9:                              // sine function with normal (=> periodic motion around the initialBodyCenter!)
      angle = mu2[bodyId] * elapsedTime; // F1B2;
      switch(returnMode) {
        case 1: // return body center
          if(angle >= liftStartAngle1[bodyId] && angle <= liftEndAngle1[bodyId]) {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = amplitude[bodyId] / mu2[bodyId] * sin(angle) * normal[bodyId * nDim + n];
            }
          } else if(angle < liftStartAngle1[bodyId]) {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = 0;
            }
          } else if(angle > liftEndAngle1[bodyId]) {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = amplitude[bodyId] / mu2[bodyId] * sin(liftEndAngle1[bodyId]) * normal[bodyId * nDim + n];
            }
          }
          break;
        case 2: // return body velocity
          if(angle > liftStartAngle1[bodyId] && angle <= liftEndAngle1[bodyId]) {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = amplitude[bodyId] * cos(angle) * normal[bodyId * nDim + n];
            }
          } else {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = 0;
            }
          }
          break;
        case 3: // return body acceleration
          if(angle > liftStartAngle1[bodyId] && angle <= liftEndAngle1[bodyId]) {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = -amplitude[bodyId] * mu2[bodyId] * sin(angle) * normal[bodyId * nDim + n];
            }
          } else {
            for(MInt n = 0; n < nDim; n++) {
              bodyData[n] = 0;
            }
          }
          break;
        default:
          for(MInt n = 0; n < nDim; n++) {
            bodyData[n] = F0;
          }
          break;
      }
      break;

    default:
      mTerm(1, AT_, "function type not implemented. Please check bodyMovementFunctions property!");
  }

  // add the initialBodyCenter to the bodyData for the body position:
  if(returnMode == 1) {
    for(MInt dir = 0; dir < nDim; dir++) {
      bodyData[dir] += initialBodyCenter[bodyId * nDim + dir];
    }
  }

  // rotate the final result around the z-axis
  if(rotAngle[bodyId] > 0 || rotAngle[bodyId] < 0) {
    MFloat tmp0 = bodyData[0] * cos(rotAngle[bodyId]) + bodyData[1] * sin(rotAngle[bodyId]);
    MFloat tmp1 = bodyData[1] * cos(rotAngle[bodyId]) - bodyData[0] * sin(rotAngle[bodyId]);
    bodyData[0] = tmp0;
    bodyData[1] = tmp1;
    bodyData[2] = bodyData[2];
  }
}


/**
 * \brief Calculates the closest distance between a given point and a sphere
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the sphere
 */
template <MInt nDim>
MFloat RigidBodies<nDim>::getDistanceSphere(const MFloat* const coordinates, const MInt bodyId) {
  TRACE();
  MFloat dist = .0;
  for(MInt n = 0; n < nDim; n++) {
    dist += POW2(coordinates[n] - a_bodyCenter(bodyId, n));
  }

  dist = sqrt(dist) - a_bodyRadius(bodyId);

  return dist;
}

/**
 * \brief Calculates the closest distance between a given point and a piston
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the piston
 */
template <MInt nDim>
MFloat RigidBodies<nDim>::getDistancePiston(const MFloat* const coordinates, const MInt bodyId) {
  TRACE();
  const MFloat x = coordinates[0] - a_bodyCenter(bodyId, 0);
  const MFloat y = coordinates[1] - a_bodyCenter(bodyId, 1);
  const MFloat z = coordinates[2] - a_bodyCenter(bodyId, 2);

  MFloat dist_cylinder = sqrt(POW2(x) + POW2(z)) - a_bodyRadius(bodyId);
  MFloat dist_caps = -y - m_bodyHeight / 2;

  return std::max(dist_cylinder, dist_caps);
}

/**
 * \brief Calculates the closest distance between a given point and a cube
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the cube
 */
template <MInt nDim>
MFloat RigidBodies<nDim>::getDistanceCube(const MFloat* const coordinates, const MInt bodyId) {
  TRACE();
  const MFloat x = coordinates[0] - a_bodyCenter(bodyId, 0);
  const MFloat y = coordinates[1] - a_bodyCenter(bodyId, 1);

  const MFloat lim_x = std::max(x, -x);
  const MFloat lim_y = std::max(y, -y);

  MFloat dist = std::max(lim_x, lim_y);

  IF_CONSTEXPR(nDim == 3) {
    const MFloat z = coordinates[2] - a_bodyCenter(bodyId, 2);
    const MFloat lim_z = std::max(z, -z);
    dist = std::max(dist, lim_z);
  }

  return dist - a_bodyRadius(bodyId);
}

/**
 * \brief Calculates the closest distance between a given point and a ellipsoid
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the ellipsoid
 */
template <MInt nDim>
MFloat RigidBodies<nDim>::getDistanceEllipsoid(const MFloat* const coordinates, const MInt bodyId) {
  TRACE();
  MFloat relPos[nDim]{};
  for(MInt n = 0; n < nDim; n++) {
    relPos[n] = coordinates[n] - a_bodyCenter(bodyId, n);
  }

  const MFloat* const bodyQuaternion = &a_bodyQuaternionT1(bodyId, 0);
  transformToBodyFrame(bodyQuaternion, relPos);

  MFloat dist = findClosestPointEllipsoid(&relPos[0], bodyId);

  return dist;
}

/**
 * \brief Calculates the closest distance between a given point and a ellipsoid 2D
 *
 * Here, the point is given in the body frame of reference
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the ellipsoid
 */
template <>
MFloat RigidBodies<2>::findClosestPointEllipsoid(const MFloat* const /*relPos*/, const MInt /*bodyId*/) {
  std::cerr << "Find Ellipse closest point 2D needs to be implemented! " << std::endl;
  return 0.0;
}

/**
 * \brief Calculates the closest distance between a given point and a ellipsoid in 3D
 *
 * Here, the point is given in the body frame of reference
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param coordinates[in] Point from which the distance is calculated
 * \param bodyId[in]      Body id to which the distance is calculated
 *
 * \returns Closest distance to the ellipsoid
 */
template <>
MFloat RigidBodies<3>::findClosestPointEllipsoid(const MFloat* const relPos, const MInt bodyId) {
  TRACE();
  std::ignore = bodyId;
  MFloat e[3] = {a_bodyRadii(bodyId, 0), a_bodyRadii(bodyId, 1), a_bodyRadii(bodyId, 2)};
  MFloat y[3] = {relPos[0], relPos[1], relPos[2]};
  MFloat x[3]{};

  // Determine reflections for y to the first octant.
  bool reflect[3];
  int i;
  int j;
  for(i = 0; i < 3; ++i) {
    reflect[i] = (y[i] < .0);
  }

  // Determine the axis order for decreasing extents.
  int permute[3];
  if(e[0] < e[1]) {
    if(e[2] < e[0]) {
      permute[0] = 1;
      permute[1] = 0;
      permute[2] = 2;
    } else if(e[2] < e[1]) {
      permute[0] = 1;
      permute[1] = 2;
      permute[2] = 0;
    } else {
      permute[0] = 2;
      permute[1] = 1;
      permute[2] = 0;
    }
  } else {
    if(e[2] < e[1]) {
      permute[0] = 0;
      permute[1] = 1;
      permute[2] = 2;
    } else if(e[2] < e[0]) {
      permute[0] = 0;
      permute[1] = 2;
      permute[2] = 1;
    } else {
      permute[0] = 2;
      permute[1] = 0;
      permute[2] = 1;
    }
  }

  std::array<MInt, 3> invpermute{};
  for(i = 0; i < 3; ++i) {
    invpermute[permute[i]] = i;
  }

  std::array<MFloat, 3> locE{};
  std::array<MFloat, 3> locY{};
  for(i = 0; i < 3; ++i) {
    j = permute[i];
    locE[i] = e[j];
    locY[i] = y[j];
    if(reflect[j]) {
      locY[i] = -locY[i];
    }
  }

  std::array<MFloat, 3> locX{};
  for(i = 0; i < 3; ++i) {
    ASSERT(!(locY[i] < F0), "");
    locY[i] = mMax(MFloatEps, locY[i]); // guarantee non-zero entries
  }
  MFloat distance = distancePointEllipsoid(locE, locY, locX);

  // Restore the axis order and reflections.
  for(i = 0; i < 3; ++i) {
    j = invpermute[i];
    if(reflect[i]) {
      locX[j] = -locX[j];
    }
    x[i] = locX[j];
  }

  std::ignore = x;

  if((POW2(y[0] / e[0]) + POW2(y[1] / e[1]) + POW2(y[2] / e[2])) < F1) {
    distance = -distance;
  }

  return distance;
}

/**
 * \brief Calculates the closest distance between a given point and a ellipsoid in 3D
 *
 * Here, the point is given in the body frame of reference _and_ the axis have been sorted by decreasing extension
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param e[3][in]  Extension of the ellipsoid in each direction
 * \param y[3][in]  Point from which the distance is calculated
 * \param x[3][out] Closest point on the ellipsoids surface
 *
 * \returns Closest distance to the ellipsoid
 */
template <MInt nDim>
inline MFloat RigidBodies<nDim>::distancePointEllipsoid(const std::array<MFloat, 3> e, const std::array<MFloat, 3> y,
                                                        std::array<MFloat, 3> x) {
  constexpr MFloat eps0 = 1e-12;
  constexpr MFloat eps1 = 1e-6;
  const MFloat dist0 = c_cellLengthAtMaxLevel();
  const MFloat dist1 = 10.0 * dist0;
  MFloat distance;
  // Bisect to compute the root of F(t) for t >= -e2*e2.
  MFloat esqr[3] = {e[0] * e[0], e[1] * e[1], e[2] * e[2]};
  MFloat ey[3] = {e[0] * y[0], e[1] * y[1], e[2] * y[2]};
  MFloat r[3];
  MFloat t0 = -esqr[2] + ey[2];
  MFloat t1 = -esqr[2] + sqrt(ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]);
  MFloat t = t0;
  const int imax = 2 * std::numeric_limits<MFloat>::max_exponent;
  MFloat eps = eps1;
  if(fabs(t1 - t0) < eps) t1 = t0 + F2 * eps;
  MInt i = 0;
  distance = 1.0;
  while(fabs(t1 - t0) > eps && i < imax) {
    while(fabs(t1 - t0) > eps && i < imax) {
      i++;
      t = F1B2 * (t0 + t1);
      r[0] = ey[0] / (t + esqr[0]);
      r[1] = ey[1] / (t + esqr[1]);
      r[2] = ey[2] / (t + esqr[2]);
      MFloat f = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] - F1;
      // f==0 also means convergence, intentionally?
      if(f >= F0) {
        t0 = t;
      }
      if(f <= F0) {
        t1 = t;
      }
    }
    x[0] = esqr[0] * y[0] / (t + esqr[0]);
    x[1] = esqr[1] * y[1] / (t + esqr[1]);
    x[2] = esqr[2] * y[2] / (t + esqr[2]);
    distance = sqrt(POW2(x[0] - y[0]) + POW2(x[1] - y[1]) + POW2(x[2] - y[2]));
    // adjust eps: small close to surface, larger further away from surface
    eps = eps0 + mMin(F1, mMax(F0, (distance - dist0) / (dist1 - dist0))) * (eps1 - eps0);
  }

  if(fabs(t0 - t1) > eps) {
    std::cerr << std::setprecision(16) << "ellipsoid dist not converged! " << i << " " << t1 - t0 << std::endl;
  }

  return distance;
}

/**
 * \brief Write restart file for the RigidBody solver
 *
 * Checks if the restart file should be saved in the current time step
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] writeRestart True if a restart file should be written
 * \param[in] writeBackup True if a backup restart file should be written (e.g. for the initial conditions)
 * \param[in] gridFileName Name of the corresponding grid file name
 * \param[in] recalcIds Use id ordering give in grid file
 */
template <MInt nDim>
void RigidBodies<nDim>::writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString /**/,
                                         MInt* /*recaldIdTree*/) {
  TRACE();
  if(writeRestart || globalTimeStep % restartInterval() == 0) {
    saveBodyRestartFile(writeBackup);
  }
}

/**
 * \brief Write restart file for the RigidBody solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param MBool backup True if backup restart file should be written (e.g. for the initial conditions)
 */
template <MInt nDim>
void RigidBodies<nDim>::saveBodyRestartFile(const MBool) {
  TRACE();

  using namespace maia::parallel_io;

  const MLong noLocal = m_noLocalBodies;
  const MLong noGlobal = noEmbeddedBodies();

  // Find order in which the data will be written
  // Sorting by containing partition cell ensures consistency between runs with different partitioning
  std::vector<MInt> localPartitionCellOfBody(noLocal, -1);
  for(MUint i = 0; i < noLocal; i++) {
    localPartitionCellOfBody[i] = grid().raw().findContainingPartitionCell(&a_bodyCenter(i, 0), -1, nullptr);
  }
  std::vector<int> indices(noLocal);
  // Init unity index mapping
  std::iota(indices.begin(), indices.end(), 0);
  // Sort indices based on the local partition cell index
  std::sort(indices.begin(), indices.end(),
            [&](int a, int b) -> bool { return localPartitionCellOfBody[a] < localPartitionCellOfBody[b]; });

  if(domainId() == 0) {
    std::cerr << "writing body-restart-data ... (rigid bodies) ";
  }

  static constexpr MLong DOF_SCALAR = 1;
  static constexpr MLong DOF_TRANS = nDim;
  static constexpr MLong DOF_ROT = nRot;
  static constexpr MLong DOF_QUAT = nDim + 1;

  const MLong offset = getDomainOffset();

  std::cout << "Domain " << domainId() << " has " << noLocal << " and offset " << offset << std::endl;

  std::stringstream fileNameStream;
  fileNameStream << outputDir() << "restartBodyData_" << globalTimeStep;
  fileNameStream << ParallelIo::fileExt();

  const MString fileName = fileNameStream.str();

  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());

  // Creating file header.
  parallelIo.defineScalar(PIO_INT, "noBodies");

  // Define scalars
  parallelIo.defineArray(PIO_FLOAT, "bodyTemperature", noGlobal * DOF_SCALAR);

  // Define vectors
  parallelIo.defineArray(PIO_FLOAT, "bodyCenter", noGlobal * DOF_TRANS);
  parallelIo.defineArray(PIO_FLOAT, "bodyVelocity", noGlobal * DOF_TRANS);
  parallelIo.defineArray(PIO_FLOAT, "bodyAcceleration", noGlobal * DOF_TRANS);

  // Define rotational vectors
  parallelIo.defineArray(PIO_FLOAT, "angularVelocityBodyT1B2", noGlobal * DOF_ROT);
  parallelIo.defineArray(PIO_FLOAT, "angularAccelerationBody", noGlobal * DOF_ROT);

  // Define rotational bodyQuaternion
  parallelIo.defineArray(PIO_FLOAT, "bodyQuaternionT1B2", noGlobal * DOF_QUAT);
  parallelIo.defineArray(PIO_FLOAT, "bodyQuaternionT1", noGlobal * DOF_QUAT);

  // Write data
  parallelIo.writeScalar(noEmbeddedBodies(), "noBodies");

  // Scalar
  parallelIo.setOffset(noLocal * DOF_SCALAR, offset * DOF_SCALAR);
  {
    MFloatScratchSpace bufTemperature(noLocal * DOF_SCALAR, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      bufTemperature[i] = a_bodyTemperature(indices[i]);
    }
    parallelIo.writeArray(bufTemperature.data(), "bodyTemperature");
  }

  // Translation
  parallelIo.setOffset(noLocal * DOF_TRANS, offset * DOF_TRANS);
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_TRANS, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_TRANS; n++) {
        writeBuffer[i * DOF_TRANS + n] = a_bodyCenter(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "bodyCenter");
  }
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_TRANS, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_TRANS; n++) {
        writeBuffer[i * DOF_TRANS + n] = a_bodyVelocity(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "bodyVelocity");
  }
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_TRANS, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_TRANS; n++) {
        writeBuffer[i * DOF_TRANS + n] = a_bodyAcceleration(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "bodyAcceleration");
  }

  // Rotation (vector)
  parallelIo.setOffset(noLocal * DOF_ROT, offset * DOF_ROT);
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_ROT, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_ROT; n++) {
        writeBuffer[i * DOF_ROT + n] = a_angularVelocityBodyT1B2(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "angularVelocityBodyT1B2");
  }
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_ROT, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_ROT; n++) {
        writeBuffer[i * DOF_ROT + n] = a_angularAccelerationBody(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "angularAccelerationBody");
  }

  // Rotation (quaternions)
  parallelIo.setOffset(noLocal * DOF_QUAT, offset * DOF_QUAT);
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_QUAT, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_QUAT; n++) {
        writeBuffer[i * DOF_QUAT + n] = a_bodyQuaternionT1B2(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "bodyQuaternionT1B2");
  }
  {
    MFloatScratchSpace writeBuffer(noLocal * DOF_QUAT, AT_, "writeBuffer");
    for(MInt i = 0; i < noLocal; i++) {
      for(MInt n = 0; n < DOF_QUAT; n++) {
        writeBuffer[i * DOF_QUAT + n] = a_bodyQuaternionT1(indices[i], n);
      }
    }
    parallelIo.writeArray(writeBuffer.data(), "bodyQuaternionT1");
  }

  if(domainId() == 0) {
    std::cerr << "ok" << std::endl;
  }
}

/**
 * \brief loadind no of embedded bodies and it's position from restart file
 *
 * \author Grafen <johannes.grafen@rwth-aachen.de>
 *
 */
template <MInt nDim>
void RigidBodies<nDim>::loadBodiesSizeAndPosition() {
  TRACE();

  using namespace maia::parallel_io;

  std::stringstream fileNameStream;
  fileNameStream.clear();

  fileNameStream << m_restartDir << "restartBodyData_" << m_restartTimeStep;
  fileNameStream << ParallelIo::fileExt();

  const MString fileName = fileNameStream.str();

  if(domainId() == 0) {
    std::cerr << "loading noBodies and positions " << fileNameStream.str() << " at " << globalTimeStep << "...";
  }

  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  parallelIo.setOffset(1, 0);
  parallelIo.readScalar(&m_noEmbeddedBodies, "noBodies");
  m_maxNoBodies = m_noEmbeddedBodies;

  m_bResFile.reset(noEmbeddedBodies());
  m_bResFile.append(noEmbeddedBodies());
  parallelIo.setOffset(noEmbeddedBodies() * nDim, 0);
  parallelIo.readArray(&m_bResFile.bodyCenter(0, 0), "bodyCenter");

  if(domainId() == 0) {
    std::cerr << "reading noBodies and initial position from restart-file done" << std::endl;
  }
}

/**
 * \brief Loads restart file for the RigidBodies solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::loadBodyRestartFile() {
  TRACE();

  using namespace maia::parallel_io;

  std::stringstream fileNameStream;
  fileNameStream.clear();

  fileNameStream << m_restartDir << "restartBodyData_" << m_restartTimeStep;
  fileNameStream << ParallelIo::fileExt();

  const MString fileName = fileNameStream.str();

  const MLong DOF = noEmbeddedBodies();
  const MLong DOF_TRANS = DOF * nDim;
  const MLong DOF_ROT = DOF * nDim;
  const MLong DOF_QUAT = DOF * (nDim + 1);

  if(domainId() == 0) {
    std::cerr << "loading body restart file " << fileNameStream.str() << " at " << globalTimeStep << "...";
  }

  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  parallelIo.setOffset(DOF, 0);
  parallelIo.readArray(&m_bResFile.bodyTemperature(0), "bodyTemperature");

  parallelIo.setOffset(DOF_TRANS, 0);
  parallelIo.readArray(&m_bResFile.bodyVelocity(0, 0), "bodyVelocity");
  parallelIo.readArray(&m_bResFile.bodyAcceleration(0, 0), "bodyAcceleration");

  parallelIo.setOffset(DOF_ROT, 0);
  parallelIo.readArray(&m_bResFile.angularVelocityBodyT1B2(0, 0), "angularVelocityBodyT1B2");
  parallelIo.readArray(&m_bResFile.angularAccelerationBody(0, 0), "angularAccelerationBody");

  parallelIo.setOffset(DOF_QUAT, 0);
  parallelIo.readArray(&m_bResFile.bodyQuaternionT1B2(0, 0), "bodyQuaternionT1B2");
  parallelIo.readArray(&m_bResFile.bodyQuaternionT1(0, 0), "bodyQuaternionT1");

  if(domainId() == 0) {
    std::cerr << "reading bodyRestartFile done" << std::endl;
  }
}

/**
 * \brief Calculate collisions forces for all bodies
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void RigidBodies<nDim>::collideBodies() {
  TRACE();
  for(MInt bodyA = 0; bodyA < noConnectingBodies(); bodyA++) {
    // collide local with local Bodies
    for(MInt bodyB = 0; bodyB < noConnectingBodies(); bodyB++) {
      if(bodyA == bodyB) {
        continue;
      }
      collideSpheres(bodyA, bodyB);
    }
  }
}

/**
 * \brief Calculate collision force between two spheres
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param bodyA[in] Id of the first body
 * \param bodyB[in] Id of the second body
 */
template <MInt nDim>
void RigidBodies<nDim>::collideSpheres(const MInt bodyA, const MInt bodyB) {
  TRACE();
  // TODO labels:toenhance support bodies with different radii
  const MFloat S = 2.0 * c_cellLengthAtMaxLevel();
  const MFloat D = a_bodyRadius(bodyA) + a_bodyRadius(bodyB);
  // const MFloat termv = m_uRef;

  MFloat delta = 0.0;
  for(MInt n = 0; n < nDim; n++) {
    delta += POW2(a_bodyCenter(bodyA, n) - a_bodyCenter(bodyB, n));
  }
  delta = sqrt(delta);

  const MFloat dist = delta - D;

  if(dist > S) {
    return;
  }

  if(dist < F0) {
    std::cerr << "\033[0;31m Warning:\033[0m potential overlap for bodies" << bodyA << " " << bodyB << std::endl;
  }

  const MFloat M = F1B2 * (a_bodyDensityRatio(bodyA) * a_volume(bodyA) + a_bodyDensityRatio(bodyB) * a_volume(bodyB));
  const MFloat C0 = 8.0 * M * POW2(2.0) / c_cellLengthAtMaxLevel();

#ifdef RB_DEBUG
  std::array<MFloat, nDim> debug{};
#endif
  for(MInt i = 0; i < nDim; i++) {
    const MFloat df = C0 * POW2(mMax(F0, -(dist - S) / S)) * (a_bodyCenter(bodyA, i) - a_bodyCenter(bodyB, i)) / dist;
#ifdef RB_DEBUG
    debug[i] = df;
#endif
    a_bodyForce(bodyA, i) += df;
  }
#ifdef RB_DEBUG
  if(dist > F0) {
    std::cerr << "Contact force applied " << debug[0] << " " << debug[1] << " " << dist << " " << S << std::endl;
  }
#endif
}

// Explicit instatiation
template class RigidBodies<2>;
template class RigidBodies<3>;
