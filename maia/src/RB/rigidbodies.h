// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _H_MAIA_RIGIDBODIES_
#define _H_MAIA_RIGIDBODIES_

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/alloc.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/maiamath.h"
#include "cartesiansolver.h"
#include "rigidbodiescollector.h"
#include "typetraits.h"

#include <algorithm>
#include <fstream>
#include <iostream>

// Activate this pragma to get detailed debug output
//#define RB_DEBUG

template <MInt nDim>
class RigidBodies : public maia::CartesianSolver<nDim, RigidBodies<nDim>> {
  // Cartesian solver
 public:
  using CartesianSolver = typename maia::CartesianSolver<nDim, RigidBodies>;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  using Geom = Geometry<nDim>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;

  using Status = maia::rb::collector::Status;

  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::grid;
  using CartesianSolver::haloCellId;
  using CartesianSolver::localPartitionCellOffsets;
  using CartesianSolver::m_adaptation;
  using CartesianSolver::m_bandWidth;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_innerBandWidth;
  using CartesianSolver::m_Ma;
  using CartesianSolver::m_outerBandWidth;
  using CartesianSolver::m_Re;
  using CartesianSolver::m_recalcIds;
  using CartesianSolver::m_residualInterval;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solutionOffset;
  using CartesianSolver::m_solutionTimeSteps;
  using CartesianSolver::m_solverId;
  using CartesianSolver::maxLevel;
  using CartesianSolver::maxNoGridCells;
  using CartesianSolver::maxRefinementLevel;
  using CartesianSolver::maxUniformRefinementLevel;
  using CartesianSolver::minLevel;
  using CartesianSolver::mpiComm;
  using CartesianSolver::neighborDomain;
  using CartesianSolver::noDomains;
  using CartesianSolver::noHaloCells;
  using CartesianSolver::noNeighborDomains;
  using CartesianSolver::noWindowCells;
  using CartesianSolver::outputDir;
  using CartesianSolver::restartDir;
  using CartesianSolver::restartInterval;
  using CartesianSolver::solverId;
  using CartesianSolver::solverMethod;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCellId;

  Geom* m_geometry;

  /// Access the solver's geometry
  Geom& geometry() const { return *m_geometry; }

  void preTimeStep() override;
  MBool solutionStep() override;
  void postTimeStep(){};

  void initSolver() override{};
  void finalizeInitSolver() override{};
  void saveSolverSolution(const MBool, const MBool){};
  void cleanUp(){};

  MFloat time() const override { return 0.0; };
  MInt noInternalCells() const override { return 0; };

  MBool prepareRestart(MBool writeRestart, MBool& /*writeRestartGrid*/) override { return writeRestart; }
  void reIntAfterRestart(MBool /*doneRestart*/) override{};
  void writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString /**/,
                        MInt* /*recaldIdTree*/) override;

  MInt a_hasNeighbor(const MInt, const MInt) const { mTerm(1, AT_, "Not implemented for this solver!"); }

  MPI_Comm mpiComm() const { return globalMaiaCommWorld(); }

  // Adaptation related functions wich need to be overriden
  // Most of them are meanignless for the RB solver itself
  void prepareAdaptation(std::vector<std::vector<MFloat>>& NotUsed(sensors),
                         std::vector<MFloat>& NotUsed(sensorWeight),
                         std::vector<std::bitset<64>>& NotUsed(sensorCellFlag),
                         std::vector<MInt>& NotUsed(sensorSolverId)) override{};
  void prepareAdaptation() override{};
  void setSensors(std::vector<std::vector<MFloat>>& /*sensors*/, std::vector<MFloat>& /*sensorWeight*/,
                  std::vector<std::bitset<64>>& /*sensorCellFlag*/, std::vector<MInt>& /*sensorSolverId*/) override{};
  // TODO dummy functions that are required to allow adaptation for another solver in a multisolver
  // computation without changing the grid for the DG solver
  void resizeGridMap() {
    // Note: resize with current tree size since the number of cells should not change
    grid().resizeGridMap(grid().tree().size());
  }
  void postAdaptation() override{};
  void finalizeAdaptation() override{};
  void removeChilds(const MInt) override{};

  void reinitAfterAdaptation() override{};
  void refineCell(const MInt) override{};
  void swapCells(const MInt, const MInt) override{};
  void swapProxy(const MInt cellId0, const MInt cellId1) { grid().swapGridIds(cellId0, cellId1); }
  void setCellWeights(MFloat* /*solverCellWeight*/) override{};

 private:
  static constexpr MInt nRot = (nDim == 3) ? 3 : 1;
  static constexpr MInt nQuat = (nDim == 3) ? 4 : 0;

  // TODO labels:toenhance Shift types and rm dummy
  struct Shapes {
    enum {
      Dummy,
      Sphere,
      Piston,
      Cube,
      Ellipsoid,

      _count
    };
  };

  // Timestep control
  MBool m_newTimeStep = true;


  // Datastrucutre: collector (localBodies | remoteBodies | dummyBodies)
  maia::rb::collector::RigidBodyCollector<nDim> m_bodies;

  // additional collector to maintain structure of deprecated restartFiles, to be removed
  maia::rb::collector::RigidBodyCollector<nDim> m_bResFile;

  MFloat m_bodyHeight{};

  // TODO labels:toenhance Use proper reference values incstead, use Prandtl ...
  const MFloat m_bodyHeatCapacity = 1.0;
  MFloat m_uniformBodyTemperature{};

  // Should be removed ~jv
  MString m_outputDir;

  // Properties
  MInt m_maxNoBodies{};
  MInt m_noEmbeddedBodies = 1;
  MString m_bodyCenterInitMethod{};
  MInt m_bodyType = 1;
  MInt m_restartTimeStep = 0;
  MBool m_forcedMotion = false;
  MBool m_translation = true;
  MBool m_rotation = true;
  MBool m_rotXYaxis = true;
  MBool m_restart = false;
  std::array<MFloat, nDim> gravity{};
  MString m_restartDir;
  std::vector<MString> m_logVars;
  std::vector<MString> m_logAngVars;

  // Reference values
  MFloat m_uRef{};

  // Time integration stuff
  MFloat m_timestep = 0;
  MInt m_timesteps = 0;

  std::vector<MFloat> m_initialBodyCenters;
  std::vector<MFloat> m_initBodyDensityRatios;
  // For Spheres, i.e. one radius per body
  std::vector<MFloat> m_initBodyRadius;
  // For Ellipsoids, i.e. nDim radii per body
  std::vector<MFloat> m_initBodyRadii;

  // ... computeBodyProperties, not sure why named as static
  MBool m_static_computeBodyProperties_first = true;
  MFloat* m_static_computeBodyProperties_amplitude = nullptr;
  MFloat* m_static_computeBodyProperties_freqFactor = nullptr;
  MFloat* m_static_computeBodyProperties_initialBodyCenter = nullptr;
  MFloat* m_static_computeBodyProperties_Strouhal = nullptr;
  MFloat* m_static_computeBodyProperties_mu = nullptr;
  MFloat* m_static_computeBodyProperties_mu2 = nullptr;
  MFloat* m_static_computeBodyProperties_liftStartAngle1 = nullptr;
  MFloat* m_static_computeBodyProperties_liftEndAngle1 = nullptr;
  MFloat* m_static_computeBodyProperties_liftStartAngle2 = nullptr;
  MFloat* m_static_computeBodyProperties_liftEndAngle2 = nullptr;
  MFloat* m_static_computeBodyProperties_circleStartAngle = nullptr;
  MFloat* m_static_computeBodyProperties_normal = nullptr;
  MInt* m_static_computeBodyProperties_bodyToFunction = nullptr;
  MFloat* m_static_computeBodyProperties_omega = nullptr;
  MFloat* m_static_computeBodyProperties_rotAngle = nullptr;

  void predictorStep();
  void correctorStep();
  void advanceBodies();

  inline void predictRotation(const MInt bodyId);
  inline void correctRotation(const MInt bodyId);

  inline void transformToWorldFrame(const MFloat* const quaternion, const MFloat* const vectorBody,
                                    MFloat* const vectorWorld);

  inline void transformToBodyFrame(const MFloat* const quaternion, const MFloat* const vectorWorld,
                                   MFloat* const vectorBody);

  inline void transformToBodyFrame(const MFloat* const quaternion, MFloat* const vecInOut);

  inline void advanceQuaternion(const MFloat* const angularVelocity, const MFloat dt, const MFloat* const qIn,
                                MFloat* const qOut);

  inline void quaternionToEuler(const MFloat* const quaternion, MFloat* const angles);

  void initTimers();
  void averageTimer();
  void readProperties();
  void initBodyData();

  struct Timers {
    enum {
      TimerGroup,
      Class,

      Constructor,
      SolutionStep,

      Predict,
      Correct,
      Communication,
      UpdateCommunication,

      _count
    };
  };

  std::array<MInt, Timers::_count> m_timers;
  MString m_timerType;

  MFloat findClosestPointEllipsoid(const MFloat* const relPos, const MInt bodyId);
  inline MFloat distancePointEllipsoid(const std::array<MFloat, 3> e, const std::array<MFloat, 3> y,
                                       std::array<MFloat, 3> x);

  void collideBodies();
  void collideSpheres(const MInt bodyA, const MInt bodyB);

  MString outputDir() const { return m_outputDir; }

  MFloat c_cellLengthAtMaxLevel() const { return grid().cellLengthAtLevel(currMaxLevel()); }

  MFloat currMaxLevel() const { return m_currMaxLevel; }

  MFloat globBBox(MInt n) { return grid().raw().globalBoundingBox()[n]; }

  // communication relations:
  // number of bodies that are local (center of mass is located in partition) for specific partition, contains
  // localBodyIds
  MInt m_noLocalBodies = 0;
  // number of bodies that are partially on partition, but it's center is located on neighbouring domain, contains
  // remoteBodyIds
  MInt m_noRemoteBodies = 0;
  // number of dummy bodies, that are added to collector, if self-mapping occurs (due to periodic boundaries,
  // the body would have two different body centers)
  MInt m_noDummyBodies = 0;

  // vector of indirect neighbouring domains of a local body
  std::vector<MInt> m_indirectHomeDomains;
  // all bodies that are marked to be transfered
  std::vector<std::vector<MInt>> m_transferBodies;
  // in the case of self-mapping contains the actual bodyId and the bodyId of the corresponding "dummy" body
  std::map<MInt, MInt> m_associatedDummyBodies;
  // maps every global remote domain id to local body ids, which are edgebodies in other domains
  std::map<MInt, std::vector<MInt>> m_remoteDomainLocalBodies;
  // maps for every global remote Domain id to remote (global) body ids
  std::map<MInt, std::vector<MInt>> m_homeDomainRemoteBodies;
  // map from domain to body to periodic shift
  std::map<MInt, std::map<MInt, std::array<MFloat, nDim>>> m_periodicShift;
  // vector of globalBodyIds
  std::vector<MInt> m_globalBodyIds;

  MBool m_periodicBC = false;
  std::array<MFloat, nDim> m_globDomainLength;
  // Transforms a local body to a remote body
  // collectorId and oldBodyId might differ because of multiple swaps
  MInt localToRemoteBody(const MInt collectorId, const MInt oldBodyId, const MInt domain);
  void deleteIrrelevantBodies(std::list<std::pair<MInt, MInt>>& removeList);
  void deleteDummyBody(const MInt bodyId, const MInt collectorShift = 0);

 public:
  RigidBodies(const MInt solverId, GridProxy& gridProxy, Geom& geometry, MPI_Comm comm);
  ~RigidBodies();

  std::array<MFloat, 2> m_distFac{};
  MFloat m_infoDiameter{};
  MInt m_maxSensorLevel{};
  MInt m_currMaxLevel{};
  MBool m_printKineticEnergy = false;
  MBool m_logBodyData = false;
  MBool m_bodyWasDeleted = false;
  MInt m_intervalBodyLog = 1;
  std::ofstream m_log;
  std::ofstream m_anglog;
  std::ofstream m_debugFileStream;

  void initBodyLog();
  void logBodyData();
  void printScalingVariables();
  void totalKineticEnergy();

  void computeBodies();
  void correctBodies();
  void checkDummyBodiesForSwap();
  void updateConnections();
  void exchangeKinematicData();
  void exchangeFsiData();
  void boxSeedInitialize();
  void linSpace(MFloat start, MFloat end, MInt num, std::vector<MFloat>& linspaced);
  void computeBodyPropertiesForced(MInt returnMode, MFloat* bodyData, MInt bodyId, MFloat time);

  void initGridProperties();
  void updateInfoDiameter(const MBool initCall = false);
  void updateMaxLevel(MInt newMaxLevel) { m_currMaxLevel = newMaxLevel; }
  void findLocalBodies();
  void findTransferBodies();
  void updateBodyDomainConnections(MBool initCall);
  void createDummyBody(MInt bodyId);

  // exchange methods
  void exchangeEdgeBodyData();
  void initRemoteDummyBodyData(const MInt bodyId);
  void exchangeNeighborConnectionInfo();
  void exchangeTransferBodies();
  void exchangeNeighborNeighborRemote(std::vector<std::vector<MInt>>& bodyIdsForRemoteDomain,
                                      std::vector<std::vector<MInt>>& homeDomainIdsForRemoteDomain,
                                      std::map<MInt, std::vector<MInt>>& tmpRemoteMap);
  void exchangeNeighborNeighborHome(std::vector<std::vector<MInt>>& bodyIdsForHomeDomain,
                                    std::vector<std::vector<MInt>>& remoteDomainIdsForHomeDomain,
                                    std::vector<std::vector<std::array<MFloat, nDim>>>& shiftForHomeDomain);
  MFloat calculatePeriodicShift(MInt intersectingCellId, MInt direction);

  void exchangeBufferLengthsNeighbor(std::vector<MInt>& noToRecv, std::vector<std::vector<MInt>>& bodyList);
  void exchangeBufferLengthsNeighbor(std::vector<MInt>& noToRecv, std::map<MInt, std::vector<MInt>>& bodyList);

  void exchangeBufferLengthsAllToRoot(std::vector<MInt>& noToRecv, const MInt noToSend);
  void exchangeBufferLengthsAllToAll(std::vector<MInt>& noToRecv, const MInt noToSend);
  void exchangeBufferLengthsNeighbors(std::vector<MInt>& noToRecv, const MInt noToSend);
  void exchangeBufferLengthsRemoteToHome(std::vector<MInt>& noToRecv, std::vector<MInt>& noToSend);
  void exchangeBodyVariablesAllToAll();
  MLong getDomainOffset();

  // debug
  void printBodyDomainConnections();

  /**
   * \brief Calculates the closest distance between a given point and body
   *
   * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
   *
   * \param coordinates[in] Point from which the distance is calculated
   * \param bodyId[in]      Body id to which the distance is calculated
   *
   * \returns Closest distance to the body
   */
  template <MInt bodyType>
  MFloat getDistance(const MFloat* const coordinates, const MInt bodyId) {
    ASSERT(bodyId > -1 && bodyId < noCollectorBodies(), "BodyId " << bodyId);
    MFloat dist = F1;
    if(bodyType == Shapes::Sphere) {
      dist = getDistanceSphere(coordinates, bodyId);
    } else if(bodyType == Shapes::Piston) {
      dist = getDistancePiston(coordinates, bodyId);
    } else if(bodyType == Shapes::Cube) {
      dist = getDistanceCube(coordinates, bodyId);
    } else if(bodyType == Shapes::Ellipsoid) {
      dist = getDistanceEllipsoid(coordinates, bodyId);
    } else {
      mTerm(1, AT_, "BodyType not implemented!");
    }

    return dist;
  }

  MFloat getDistanceSphere(const MFloat* coordinates, const MInt bodyId);
  MFloat getDistancePiston(const MFloat* coordinates, const MInt bodyId);
  MFloat getDistanceCube(const MFloat* coordinates, const MInt bodyId);
  MFloat getDistanceEllipsoid(const MFloat* coordinates, const MInt bodyId);

  MInt size() const { return noCollectorBodies(); }
  MInt getLocalBodyId(MInt globalBodyId) {
    std::vector<MInt>::iterator itr = std::find(m_globalBodyIds.begin(), m_globalBodyIds.end(), globalBodyId);
    if(itr != m_globalBodyIds.cend()) {
      return std::distance(m_globalBodyIds.begin(), itr);
    }
    return -1;
  }
  inline MInt getGlobalBodyId(MInt localBodyId) { return m_globalBodyIds.at(localBodyId); }

  inline MBool isLocalBody(MInt bodyId) { return bodyId < m_noLocalBodies; }

  inline MBool isRemoteBody(MInt bodyId) { return bodyId >= m_noLocalBodies; }

  inline MInt noEmbeddedBodies() const { return m_noEmbeddedBodies; }

  inline MInt noConnectingBodies() const { return m_noLocalBodies + m_noRemoteBodies; }

  inline MInt noCollectorBodies() const { return noConnectingBodies() + m_noDummyBodies; }

  inline MInt noHomeDomains() const { return m_homeDomainRemoteBodies.size(); }

  inline MInt noRemoteDomains() const { return m_remoteDomainLocalBodies.size(); }

  MBool hasAssociatedDummyBody(const MInt bodyId) const {
    if(m_associatedDummyBodies.empty()) {
      return false;
    }
    return m_associatedDummyBodies.find(bodyId) != m_associatedDummyBodies.end();
  }

  void addForce(const MInt bodyId, const std::array<MFloat, nDim>& force) {
    for(MInt n = 0; n < nDim; n++) {
      m_bodies.bodyForce(bodyId, n) += force[n];
    }
  }

  void addTorque(const MInt bodyId, const std::array<MFloat, nRot>& torque) {
    if(!m_rotXYaxis && nRot == 3) {
      a_bodyTorque(bodyId, 0) = 0.0;
      a_bodyTorque(bodyId, 1) = 0.0;
      a_bodyTorque(bodyId, 2) += torque[2];
    } else {
      for(MInt n = 0; n < nRot; n++) {
        a_bodyTorque(bodyId, n) += torque[n];
      }
    }
  }

  void resetForce() {
    // Reset bodyForces
    std::fill_n(a_bodyForce(0), noCollectorBodies() * nDim, 0.0);
  }

  void resetTorque() {
    // Reset bodyTorques
    std::fill_n(a_bodyTorque(0), noCollectorBodies() * nRot, 0.0);
  }

  MFloat a_bodyMass(const MInt bodyId) { return a_bodyDensityRatio(bodyId) * a_volume(bodyId); }

  MFloat a_volume(const MInt bodyId) {
    IF_CONSTEXPR(nDim == 2) { std::cout << "Revise volume for 2D!" << std::endl; }

    MFloat volume = 1.0;
    if(a_bodyType() == Shapes::Sphere || a_bodyType() == Shapes::Ellipsoid) {
      volume = PI * F4B3 * a_bodyRadii(bodyId, 0) * a_bodyRadii(bodyId, 1) * a_bodyRadii(bodyId, 2);
    } else {
      std::cerr << "Implement volume formula!" << std::endl;
    }
    return volume;
  }

  MFloat& a_bodyRadii(const MInt bodyId, const MInt n) { return m_bodies.bodyRadii(bodyId, n); }
  MFloat& a_bodyRadius(const MInt bodyId) { return m_bodies.bodyRadius(bodyId); }

  void getVelocity(const MInt bodyId, std::array<MFloat, nDim>& velocity) {
    for(MInt n = 0; n < nDim; n++) {
      velocity[n] = m_bodies.bodyVelocity(bodyId, n);
    }
  }

  void incrementHeatFlux(const MInt bodyId, MFloat heatflux) { m_bodies.bodyHeatFlux(bodyId) += heatflux; }

  MFloat getHeatFlux(const MInt bodyId) { return m_bodies.bodyHeatFlux(bodyId); }

  MFloat* a_bodyForce(const MInt bodyId) { return &m_bodies.bodyForce(bodyId, 0); }
  MFloat& a_bodyForce(const MInt bodyId, const MInt dim) { return m_bodies.bodyForce(bodyId, dim); }
  MFloat a_bodyForce(const MInt bodyId, const MInt dim) const { return m_bodies.bodyForce(bodyId, dim); }

  MFloat& a_bodyCenter(const MInt bodyId, const MInt dim) { return m_bodies.bodyCenter(bodyId, dim); }
  MFloat a_bodyCenter(const MInt bodyId, const MInt dim) const { return m_bodies.bodyCenter(bodyId, dim); }

  MFloat& a_bodyCenterOld(const MInt bodyId, const MInt dim) { return m_bodies.bodyCenterOld(bodyId, dim); }
  MFloat a_bodyCenterOld(const MInt bodyId, const MInt dim) const { return m_bodies.bodyCenterOld(bodyId, dim); }

  MFloat& a_bodyVelocity(const MInt bodyId, const MInt dim) { return m_bodies.bodyVelocity(bodyId, dim); }
  MFloat a_bodyVelocity(const MInt bodyId, const MInt dim) const { return m_bodies.bodyVelocity(bodyId, dim); }

  MFloat& a_bodyVelocityOld(const MInt bodyId, const MInt dim) { return m_bodies.bodyVelocityOld(bodyId, dim); }
  MFloat a_bodyVelocityOld(const MInt bodyId, const MInt dim) const { return m_bodies.bodyVelocityOld(bodyId, dim); }

  MFloat& a_bodyAcceleration(const MInt bodyId, const MInt dim) { return m_bodies.bodyAcceleration(bodyId, dim); }
  MFloat a_bodyAcceleration(const MInt bodyId, const MInt dim) const { return m_bodies.bodyAcceleration(bodyId, dim); }

  MFloat& a_bodyAccelerationOld(const MInt bodyId, const MInt dim) { return m_bodies.bodyAccelerationOld(bodyId, dim); }
  MFloat a_bodyAccelerationOld(const MInt bodyId, const MInt dim) const {
    return m_bodies.bodyAccelerationOld(bodyId, dim);
  }

  MFloat& a_bodyTemperature(const MInt bodyId) { return m_bodies.bodyTemperature(bodyId); }
  MFloat a_bodyTemperature(const MInt bodyId) const { return m_bodies.bodyTemperature(bodyId); }

  MFloat& a_bodyTemperatureOld(const MInt bodyId) { return m_bodies.bodyTemperatureOld(bodyId); }
  MFloat a_bodyTemperatureOld(const MInt bodyId) const { return m_bodies.bodyTemperatureOld(bodyId); }

  MFloat& a_bodyHeatFlux(const MInt bodyId) { return m_bodies.bodyHeatFlux(bodyId); }
  MFloat a_bodyHeatFlux(const MInt bodyId) const { return m_bodies.bodyHeatFlux(bodyId); }

  MFloat& a_bodyDensityRatio(const MInt bodyId) { return m_bodies.bodyDensityRatio(bodyId); }
  MFloat a_bodyDensityRatio(const MInt bodyId) const { return m_bodies.bodyDensityRatio(bodyId); }

  MFloat& a_bodyInertia(const MInt bodyId, const MInt dim) { return m_bodies.bodyInertia(bodyId, dim); }
  MFloat a_bodyInertia(const MInt bodyId, const MInt dim) const { return m_bodies.bodyInertia(bodyId, dim); }

  MFloat& a_bodyQuaternionT1(const MInt bodyId, const MInt dim) { return m_bodies.bodyQuaternionT1(bodyId, dim); }
  MFloat a_bodyQuaternionT1(const MInt bodyId, const MInt dim) const { return m_bodies.bodyQuaternionT1(bodyId, dim); }

  MFloat& a_bodyQuaternionT1B2(const MInt bodyId, const MInt dim) { return m_bodies.bodyQuaternionT1B2(bodyId, dim); }
  MFloat a_bodyQuaternionT1B2(const MInt bodyId, const MInt dim) const {
    return m_bodies.bodyQuaternionT1B2(bodyId, dim);
  }

  MFloat& a_angularVelocityT1(const MInt bodyId, const MInt dim) { return m_bodies.angularVelocityT1(bodyId, dim); }
  MFloat a_angularVelocityT1(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularVelocityT1(bodyId, dim);
  }

  MFloat& a_angularVelocityT1B2(const MInt bodyId, const MInt dim) { return m_bodies.angularVelocityT1B2(bodyId, dim); }
  MFloat a_angularVelocityT1B2(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularVelocityT1B2(bodyId, dim);
  }

  MFloat& a_angularVelocityBodyT1(const MInt bodyId, const MInt dim) {
    return m_bodies.angularVelocityBodyT1(bodyId, dim);
  }
  MFloat a_angularVelocityBodyT1(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularVelocityBodyT1(bodyId, dim);
  }

  MFloat& a_angularVelocityBodyT1B2(const MInt bodyId, const MInt dim) {
    return m_bodies.angularVelocityBodyT1B2(bodyId, dim);
  }
  MFloat a_angularVelocityBodyT1B2(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularVelocityBodyT1B2(bodyId, dim);
  }

  MFloat& a_angularAccelerationT1(const MInt bodyId, const MInt dim) {
    return m_bodies.angularAccelerationT1(bodyId, dim);
  }
  MFloat a_angularAccelerationT1(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularAccelerationT1(bodyId, dim);
  }

  MFloat& a_angularAccelerationBody(const MInt bodyId, const MInt dim) {
    return m_bodies.angularAccelerationBody(bodyId, dim);
  }
  MFloat a_angularAccelerationBody(const MInt bodyId, const MInt dim) const {
    return m_bodies.angularAccelerationBody(bodyId, dim);
  }

  MFloat& a_torqueT1(const MInt bodyId, const MInt dim) { return m_bodies.torqueT1(bodyId, dim); }
  MFloat a_torqueT1(const MInt bodyId, const MInt dim) const { return m_bodies.torqueT1(bodyId, dim); }

  Status& a_status(const MInt bodyId) { return m_bodies.status(bodyId); }
  Status a_status(const MInt bodyId) const { return m_bodies.status(bodyId); }

  MFloat* a_bodyTorque(const MInt bodyId) { return &m_bodies.torqueT1(bodyId, 0); }
  MFloat& a_bodyTorque(const MInt bodyId, const MInt dim) { return a_torqueT1(bodyId, dim); }
  MFloat a_bodyTorque(const MInt bodyId, const MInt dim) const { return a_torqueT1(bodyId, dim); }

  MFloat& a_angularVelocity(const MInt bodyId, const MInt dim) { return a_angularVelocityT1(bodyId, dim); }
  MFloat a_angularVelocity(const MInt bodyId, const MInt dim) const { return a_angularVelocityT1(bodyId, dim); }

  MFloat& a_angularVelocityBody(const MInt bodyId, const MInt dim) { return a_angularVelocityBodyT1(bodyId, dim); }
  MFloat a_angularVelocityBody(const MInt bodyId, const MInt dim) const { return a_angularVelocityBodyT1(bodyId, dim); }

  void getAngularVelocity(const MInt bodyId, std::array<MFloat, nRot>& angularVelocity) {
    for(MInt n = 0; n < nRot; n++) {
      angularVelocity[n] = a_angularVelocityT1(bodyId, n);
    }
  }

  MFloat a_bodyVelocityMag(const MInt bodyId) {
    MFloat velMag = 0.0;
    for(MInt n = 0; n < nDim; n++) {
      velMag += POW2(a_bodyVelocity(bodyId, n));
    }
    return sqrt(velMag);
  }

  MFloat* getBodyProperties(MInt returnMode, MInt bodyId);

  MInt a_bodyType() const { return m_bodyType; }

  void setTimestep(const MFloat timestep) { m_timestep = timestep; }

  void saveBodyRestartFile(const MBool backup);

  void loadBodyRestartFile();

  void loadBodiesSizeAndPosition();
  // dummy function for compatibility
  MInt a_noCells() const { return 0; }
};

#endif
