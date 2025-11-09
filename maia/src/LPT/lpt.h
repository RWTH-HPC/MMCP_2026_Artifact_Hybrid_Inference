// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef _LPT_H
#define _LPT_H

#include <algorithm>
#include <limits>
#include <list>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include "GRID/cartesiangrid.h"
#include "cartesiansolver.h"
#include "lptbase.h"
#include "lptbndrycell.h"
#include "lptcellcollector.h"
#include "lptcollision.h"
#include "lptellipsoidal.h"
#include "lptlib.h"
#include "lptspherical.h"
#include "lptspray.h"

//#define LPT_DEBUG

template <MInt nDim>
class SprayModel;
class CartesianSolver;

template <MInt nDim>
class MaterialState;

template <MInt nDim>
class ParticleCollision;

template <MInt nDim>
class LPTEllipsoidal;

template <MInt nDim>
using particleSendQueue = std::vector<std::queue<maia::lpt::sendQueueType<nDim>>>;

namespace maia {
namespace lpt {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    SolverType,
    PreTime,
    TimeInt,

    Motion,
    Energy,
    Wall,
    SourceTerms,

    Injection,
    Breakup,

    Exchange,
    Exchange1,
    Exchange2,
    Exchange3,
    Exchange4,
    Exchange5,

    // Special enum value used to initialize timer array
    _count
  };
};

} // namespace lpt
} // namespace maia


template <MInt nDim>
class LPT : public maia::CartesianSolver<nDim, LPT<nDim>> {
  template <MInt nDim_>
  friend class LPTBase;

  template <MInt nDim_>
  friend class LPTSpherical;

  template <MInt nDim_>
  friend class LPTEllipsoidal;

  template <MInt nDim_>
  friend class SprayModel;

  template <MInt nDim_>
  friend class CouplingParticle;

  template <MInt nDim_, class SysEqn>
  friend class CouplerFvParticle;

  template <MInt nDim_, MInt nDist, class SysEqn>
  friend class LbLpt;

  template <MInt nDim_, class SysEqn>
  friend class PostProcessingFvLPT;

  template <MInt nDim_>
  friend class PostProcessingLbLPT;

  template <MInt nDim_>
  friend class PostProcessingLPT;

  template <MInt nDim_>
  friend class ParticleCollision;

  template <MInt nDim_>
  friend class MaterialState;

 public:
  using CartesianSolver = typename maia::CartesianSolver<nDim, LPT>;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  using Geom = Geometry<nDim>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;
  using LptCellCollector = maia::lpt::collector::LptCells<nDim>;

  // used CartesianSolver
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::getIdentifier;
  using CartesianSolver::grid;
  using CartesianSolver::haloCellId;
  using CartesianSolver::isActive;
  using CartesianSolver::localPartitionCellOffsets;
  using CartesianSolver::m_bandWidth;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_innerBandWidth;
  using CartesianSolver::m_Ma;
  using CartesianSolver::m_noDirs;
  using CartesianSolver::m_outerBandWidth;
  using CartesianSolver::m_Re;
  using CartesianSolver::m_recalcIds;
  using CartesianSolver::m_residualInterval;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_revDir;
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
  using CartesianSolver::returnIdleRecord;
  using CartesianSolver::returnLoadRecord;
  using CartesianSolver::solverId;
  using CartesianSolver::solverMethod;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCellId;

  using CartesianSolver::exchangeData;

  Geom* m_geometry;

  LptCellCollector m_cells;

  struct PrimitiveVariables;
  PrimitiveVariables PV;

  Collector<LPTBndryCell<nDim>>* m_bndryCells = nullptr;

  /// Access the solver's geometry
  Geom& geometry() const { return *m_geometry; }

  /// Solver constructor
  LPT(const MInt solverId, GridProxy& gridProxy_, Geom& geometry_, const MPI_Comm comm);

  virtual ~LPT() {
  };

  LPT(const LPT&) = delete;            // Prevent copy-construction
  LPT& operator=(const LPT&) = delete; // Prevent assignment

  // partcile collision class
  std::unique_ptr<ParticleCollision<nDim>> m_collisionModel;

  //------------------- run-loop and grid-controller function calls -------------------------------

  void preTimeStep() override;
  MBool solutionStep() override;
  void postTimeStep() override;

  void initSolver() override;
  void finalizeInitSolver() override;
  void saveSolverSolution(const MBool /*unused*/, const MBool /*unused*/) override{};
  void cleanUp() override{};


  // adaptation related functions
  void prepareAdaptation() override;
  void setSensors(std::vector<std::vector<MFloat>>&, std::vector<MFloat>&, std::vector<std::bitset<64>>&,
                  std::vector<MInt>&) override;
  void swapProxy(MInt, MInt) override;
  void resizeGridMap() override;
  void postAdaptation() override;
  void finalizeAdaptation() override;
  void removeChilds(const MInt) override;
  void removeCell(const MInt) override;
  void refineCell(const MInt) override;
  void swapCells(const MInt, const MInt) override;
  MInt cellOutside(const MFloat*, const MInt, const MInt) override;
  void sensorParticle(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                      std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorInterface(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;

  /// \brief Returns the LPT adaptation forcing
  MBool forceAdaptation() override {
    this->startLoadTimer(AT_);

    recvRegridTrigger();

    // If LPTSolver is inactive on one of the ranks we need this allreduce!
    if(grid().hasInactiveRanks()) {
      MPI_Allreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, grid().raw().mpiComm(), AT_,
                    "MPI_IN_PLACE", "m_forceAdaptation");
    }

    this->stopLoadTimer(AT_);

    return m_forceAdaptation;
  }

  // restart related functions
  MBool prepareRestart(MBool force, MBool&) override;
  void writeRestartFile(const MBool, const MBool, const MString, MInt*) override;
  void writeParticleRestartFile();
  void writeCellSolutionFile(const MString&, MInt*);
  void reIntAfterRestart(MBool /*unused*/) override{};
  void saveDebugRestartFile();
  void writeRestartFile(MBool /*unused*/) override{};

  // partition related functions
  void setCellWeights(MFloat*) override;
  void getSolverTimings(std::vector<std::pair<std::string, MFloat>>&, const MBool allTimings) override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const override;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;
  void limitWeights(MFloat*) override;
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override;
  MInt noSolverTimers(const MBool allTimings) override {
#ifdef MAIA_TIMER_FUNCTION
    TERMM_IF_COND(!allTimings, "FIXME: reduced timings mode not yet supported by LPT.");
    static const MInt noAdditionTimers = 9;
    return 2 + noAdditionTimers;
#else
    return 2;
#endif
  }
  MInt noLoadTypes() const override { return 3 + m_weightSourceCells; };

  // balance related functions
  void resetSolverFull();
  void resetSolver() override;
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override;
  void cancelMpiRequests() override;

  MInt noCellDataDlb() const override {
    if(grid().isActive()) {
      return 3;
      // 0: cell int data
      // 1: particle float data
      // 2: particle int data
    } else {
      return 0;
    }
  }
  MInt cellDataTypeDlb(const MInt dataId) const override {
    if(dataId == 1) {
      return MFLOAT;
    } else if(dataId == 0 || dataId == 2) {
      return MINT;
    } else {
      TERMM(1, "solverCelldataType: invalid data id");
    }
  };
  MInt cellDataSizeDlb(const MInt dataId, const MInt cellId) override;
  /// Return solver data for DLB
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MInt* const data) override;
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override;

  /// Set solver data for DLB
  void setCellDataDlb(const MInt dataId, const MInt* const data) override;
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override;

  /// Get/set global solver variables during DLB
  void getGlobalSolverVars(std::vector<MFloat>& globalFloatVars, std::vector<MInt>& globalIdVars) override;
  void setGlobalSolverVars(std::vector<MFloat>& globalFloatVars, std::vector<MInt>& globalIdVars) override;

  MBool hasSplitBalancing() const override { return true; }

  MFloat computeTimeStep();
  MFloat timeStep() { return m_timeStep; }
  void forceTimeStep(const MFloat timeStep) { m_timeStep = timeStep; }
  MInt a_timeStepComputationInterval() { return m_timeStepComputationInterval; }

  MInt a_hasNeighbor(const MInt cellId, const MInt dir) const { return grid().tree().hasNeighbor(cellId, dir); }

  MInt mpiTag(const MString exchangeType) {
    switch((LPTmpiTag)string2enum(exchangeType)) {
      case PARTICLE_COUNT: {
        return 0;
      }
      case PARTICLE_FLOAT: {
        if(!m_primaryExchange) {
          return 1;
        } else {
          return 5;
        }
      }
      case PARTICLE_INT: {
        if(!m_primaryExchange) {
          return 2;
        } else {
          return 6;
        }
      }
      case SOURCE_TERMS: {
        return 3;
      }
      case FLOW_FIELD: {
        return 7;
      }
      case CHECK_ADAP: {
        return 8;
      }
      case VELOCITY_SLOPES: {
        return 14;
      }
      default: {
        mTerm(1, AT_, "Unknown mpiTag");
        return -1;
      }
    }
  }

  MFloat reduceData(const MInt cellId, MFloat* data, const MInt dataBlockSize = 1, const MBool average = true);

  //---------------------- accessors to LPT cell- or particle collector ---------------------------

  template <class LPTParticle = LPTSpherical<nDim>>
  MInt a_noParticles() {
    IF_CONSTEXPR(std::is_same_v<LPTParticle, LPTSpherical<nDim>>) { return m_partList.size(); }
    IF_CONSTEXPR(std::is_same_v<LPTParticle, LPTEllipsoidal<nDim>>) { return m_partListEllipsoid.size(); }
    mTerm(-1, AT_, "Unknown particle type");
  }
  MInt a_noSphericalParticles() const { return m_partList.size(); }
  MInt a_noEllipsoidalParticles() const { return m_partListEllipsoid.size(); }
  MInt a_noCells() const { return m_cells.size(); }

  template <class LPTParticle>
  std::vector<LPTParticle>& a_particleList() {
    IF_CONSTEXPR(std::is_same_v<LPTParticle, LPTSpherical<nDim>>) { return m_partList; }
    IF_CONSTEXPR(std::is_same_v<LPTParticle, LPTEllipsoidal<nDim>>) { return m_partListEllipsoid; }
    mTerm(-1, AT_, "Unknown particle type");
  }

  MFloat& a_volumeFraction(const MInt cellId) { return m_cells.volumeFraction(cellId); }
  MFloat a_volumeFraction(const MInt cellId) const { return m_cells.volumeFraction(cellId); }

  MInt& a_noParticlesInCell(const MInt cellId) { return m_cells.noParticles(cellId); }
  MInt a_noParticlesInCell(const MInt cellId) const { return m_cells.noParticles(cellId); }

  MInt& a_noEllipsoidsInCell(const MInt cellId) { return m_cells.noEllipsoids(cellId); }
  MInt a_noEllipsoidsInCell(const MInt cellId) const { return m_cells.noEllipsoids(cellId); }

  MInt& a_bndryCellId(const MInt cellId) { return m_cells.bndryCellId(cellId); }
  MInt a_bndryCellId(const MInt cellId) const { return m_cells.bndryCellId(cellId); }

  MFloat& a_fluidVariable(const MInt cellId, const MInt var) { return m_cells.fluidVariable(cellId, var); }
  MFloat a_fluidVariable(const MInt cellId, const MInt var) const { return m_cells.fluidVariable(cellId, var); }

  MFloat& a_fluidVelocity(const MInt cellId, const MInt dim) { return m_cells.fluidVariable(cellId, PV.VV[dim]); }
  MFloat a_fluidVelocity(const MInt cellId, const MInt dim) const { return m_cells.fluidVariable(cellId, PV.VV[dim]); }

  MFloat& a_fluidDensity(const MInt cellId) { return m_cells.fluidVariable(cellId, PV.RHO); }
  MFloat a_fluidDensity(const MInt cellId) const { return m_cells.fluidVariable(cellId, PV.RHO); }

  MFloat& a_fluidPressure(const MInt cellId) { return m_cells.fluidVariable(cellId, PV.P); }
  MFloat a_fluidPressure(const MInt cellId) const { return m_cells.fluidVariable(cellId, PV.P); }

  MFloat& a_fluidTemperature(const MInt cellId) { return m_cells.fluidVariable(cellId, PV.T); }
  MFloat a_fluidTemperature(const MInt cellId) const { return m_cells.fluidVariable(cellId, PV.T); }

  MFloat& a_fluidSpecies(const MInt cellId) { return m_cells.fluidSpecies(cellId); }
  MFloat a_fluidSpecies(const MInt cellId) const { return m_cells.fluidSpecies(cellId); }

  MFloat& a_heatFlux(const MInt cellId) { return m_cells.heatFlux(cellId); }
  MFloat a_heatFlux(const MInt cellId) const { return m_cells.heatFlux(cellId); }

  MFloat& a_massFlux(const MInt cellId) { return m_cells.massFlux(cellId); }
  MFloat a_massFlux(const MInt cellId) const { return m_cells.massFlux(cellId); }

  MFloat& a_momentumFlux(const MInt cellId, const MInt dim) { return m_cells.momentumFlux(cellId, dim); }
  MFloat a_momentumFlux(const MInt cellId, const MInt dim) const { return m_cells.momentumFlux(cellId, dim); }

  MFloat& a_workFlux(const MInt cellId) { return m_cells.workFlux(cellId); }
  MFloat a_workFlux(const MInt cellId) const { return m_cells.workFlux(cellId); }

  MFloat& a_velocitySlope(const MInt cellId, const MInt varId, const MInt dir) {
    return m_cells.velocitySlope(cellId, varId, dir);
  }
  MFloat a_velocitySlope(const MInt cellId, const MInt varId, const MInt dir) const {
    return m_cells.velocitySlope(cellId, varId, dir);
  }

  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  MBool a_isHalo(const MInt cellId) const { return m_cells.hasProperty(cellId, LptCell::IsHalo); }
  maia::lpt::cell::BitsetType::reference a_isHalo(const MInt cellId) {
    return m_cells.hasProperty(cellId, LptCell::IsHalo);
  }

  MBool a_isWindow(const MInt cellId) const { return m_cells.hasProperty(cellId, LptCell::IsWindow); }
  maia::lpt::cell::BitsetType::reference a_isWindow(const MInt cellId) {
    return m_cells.hasProperty(cellId, LptCell::IsWindow);
  }

  MBool a_isValidCell(const MInt cellId) const { return m_cells.hasProperty(cellId, LptCell::IsValid); }
  maia::lpt::cell::BitsetType::reference a_isValidCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, LptCell::IsValid);
  }

  MBool a_regridTrigger(const MInt cellId) const { return m_cells.hasProperty(cellId, LptCell::RegridTrigger); }
  maia::lpt::cell::BitsetType::reference a_regridTrigger(const MInt cellId) {
    return m_cells.hasProperty(cellId, LptCell::RegridTrigger);
  }

  MBool a_isBndryCell(const MInt cellId) const override { return a_bndryCellId(cellId) > -1; }

  const MFloat& a_coordinate(const MInt cellId, const MInt dim) const {
    if(a_isBndryCell(cellId)) {
      return m_bndryCells->a[a_bndryCellId(cellId)].m_coordinates[dim];
    }
    return grid().tree().coordinate(cellId, dim);
  }


  //--------------------------LPT internal functions and data ----------------------------------------

  template <class P, MInt a, MInt b>
  void interpolateAndCalcWeights(P& particle, const MInt cellId, const MFloat* position, const MFloat* interpolatedVar,
                                 std::vector<MFloat>& weight);

  template <MInt a, MInt b, MBool interpolateVelocitySlopes>
  void interpolateVariablesLS(const MInt cellId, const MFloat* const position, MFloat* const interpolatedVar);
  template <MInt a, MInt b>
  void interpolateVariablesLS(const MInt cellId, const MFloat* const position, MFloat* const interpolatedVar) {
    interpolateVariablesLS<a, b, false>(cellId, position, interpolatedVar);
  }

  void updateFluidFraction();
  MFloat injectionDuration() const { return m_sprayModel != nullptr ? m_sprayModel->injectionDuration() : m_time; }

  MFloat timeSinceSOI() const { return m_sprayModel != nullptr ? m_sprayModel->timeSinceSOI() : m_time; }

  MInt injectorCellId() const {
    if(m_sprayModel != nullptr) return m_spawnCellId;
    if(domainId() == 0) {
      return 1;
    } else {
      return -1;
    }
  }

  // TODO labels:LPT,PP move to postData as soon as its available!
  std::map<MInt, std::vector<MFloat>> m_injData;


  MInt globalNoParticles() {
    MInt sumPart = 0;
    MInt noPart = m_partList.size();
    MPI_Allreduce(&noPart, &sumPart, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "sumpart");
    return sumPart;
  }

  MInt globalNoEllipsoids() {
    MInt sumPart = 0;
    MInt noPart = m_partListEllipsoid.size();
    MPI_Allreduce(&noPart, &sumPart, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "INPLACE", "sumpart");
    return sumPart;
  }

  template <class LPTParticle>
  MBool pushToQueue(std::vector<std::map<MInt, MInt>>& pointsTo, const MInt partPos) {
    std::vector<LPTParticle>& partList = a_particleList<LPTParticle>();
    maia::lpt::sendQueueType<nDim> thisSend{};
    for(MInt d = 0; d < grid().noNeighborDomains(); d++) {
      auto pos = pointsTo[d].find(partList[partPos].m_cellId);
      // particle is inside a halo cell so send it...
      if(pos != pointsTo[d].end()) {
        thisSend.partPos = partPos;
        thisSend.toCellId = pos->second;
        m_queueToSend[d].push(thisSend);
        return true;
      }
    }
    return false;
  }

 protected:
  // Types
  std::vector<LPTSpherical<nDim>> m_partList;
  std::vector<LPTEllipsoidal<nDim>> m_partListEllipsoid;

  std::map<MInt, std::vector<MInt>> m_cellToNghbrHood;
  std::map<MInt, std::vector<MInt>> m_cellToNghbrHoodInterpolation;

  using Timers = maia::lpt::Timers_;
  // Timers
  // Timer group which holds all solver-wide timers
  MInt m_timerGroup = -1;
  // Stores all solver-wide timers
  std::array<MInt, Timers::_count> m_timers{};

  MInt m_sleepLPT = -1;
  MBool m_skipLPT = false;

  MBool m_domainIdOutput = false;

  MFloat m_weightBaseCell = 0.0;
  MFloat m_weightLeafCell = 0.05;
  MFloat m_weightParticleCell = 0.1;
  MFloat m_weightParticle = 0.1;
  MFloat m_weightSpawnCell = 5.0;
  MFloat m_weightMulitSolverFactor = 1;

  MBool m_limitWeights = false;
  MBool m_weightSourceCells = false;

  MFloat m_time = 0.0;

  // randon number for particle spawn at spanDomainId
  std::mt19937_64 m_PRNGSpawn;
  MInt m_PRNGSpawnCount = 0;
  // random number for particle respawn at respawnDomain
  std::mt19937_64 m_PRNGRespawn;
  // randon number for particle initialisation on all ranks
  // the value is not needed for a restart!
  std::mt19937_64 m_PRNGInit;

  particleSendQueue<nDim> m_queueToSend;

  std::vector<std::map<MInt, MInt>> m_pointsToHalo;
  std::vector<std::map<MInt, MInt>> m_pointsToWindow;

  std::map<MFloat, MFloat> m_terminalVelocity;

  MFloat m_timeStep{};
  MFloat m_timeStepOld{};
  MInt m_timeStepComputationInterval{};
  MBool m_timeStepUpdated = true;

  // parallel:
  MBool m_primaryExchange = false;
  MInt m_noSourceTerms = 0;
  MPI_Request* m_sourceSendRequest = nullptr;
  MPI_Request* m_sourceRecvRequest = nullptr;
  std::vector<std::unique_ptr<MFloat[]>> m_sourceRecv{};
  std::vector<std::unique_ptr<MFloat[]>> m_sourceSend{};

  MPI_Request* m_checkAdaptationRecvRequest = nullptr;
  MPI_Request* m_checkAdaptationSendRequest = nullptr;
  std::vector<MInt> m_checkAdaptationRecv{};
  std::vector<MInt> m_checkAdaptationSend{};

  MPI_Request* m_flowSendRequest = nullptr;
  MPI_Request* m_flowRecvRequest = nullptr;
  MFloat** m_flowRecv = nullptr;
  MFloat** m_flowSend = nullptr;

  MPI_Request* m_slopesSendRequest = nullptr;
  MPI_Request* m_slopesRecvRequest = nullptr;
  MFloat** m_slopesRecv = nullptr;
  MFloat** m_slopesSend = nullptr;

  std::vector<MInt> m_sendSize{};
  std::vector<MInt> m_recvSize{};

  std::vector<MPI_Request> m_mpi_reqSendFloat{};
  std::vector<MPI_Request> m_mpi_reqSendInt{};
  std::vector<MPI_Request> m_mpi_reqSendSize{};
  std::vector<MPI_Request> m_mpi_reqRecvFloat{};
  std::vector<MPI_Request> m_mpi_reqRecvInt{};
  std::vector<MPI_Status> m_mpi_statusProbe{};
  std::vector<std::unique_ptr<MFloat[]>> m_sendBuffer{};
  std::vector<std::unique_ptr<MFloat[]>> m_recvBuffer{};
  std::vector<std::unique_ptr<MInt[]>> m_intSendBuffer{};
  std::vector<std::unique_ptr<MInt[]>> m_intRecvBuffer{};
  MBool m_nonBlockingComm = false;
  MInt m_nonBlockingStage = -1;
  MBool m_receiveFlowField = false;
  MBool m_receiveVelocitySlopes = false;

  MInt m_noSendParticles = 0;

  MPI_Status m_mpi_statusInt{};
  MPI_Status m_mpi_statusFloat{};

  MBool m_openParticleInjSend = false;
  MBool m_openParticleSend = false;

  MBool m_openSourceSend = false;
  MBool m_openFlowSend = false;
  MBool m_openRegridSend = false;
  MBool m_openSlopesSend = false;

  MBool m_couplingRedist = false;
  MInt m_noRedistLayer = 1;

  MBool m_momentumCoupling = false;
  MBool m_heatCoupling = false;
  MBool m_massCoupling = false;
  MBool m_evaporation = false;

  MFloat m_sumEvapMass = 0;

  MFloat m_particleResiduum = 0.0;

  void calculateTerminalVelocities();

  MFloat calculateAverageRep();

  MFloat calculateVRMS();

  MFloat calculateTotalKineticEnergy();

  void particleRespawn();

  void spawnParticles();

  MInt addParticle(const MInt cellId, const MFloat diameter, const MFloat particleDensityRatio, const MInt random = 1,
                   const MInt addMode = 0, const MFloat* velocity = nullptr, const MFloat* position = nullptr,
                   const MInt parcelSize = 1);

  MInt addEllipsoid(const MInt cellId, const MFloat semiMinorAxis, const MFloat aspectRatio,
                    const MFloat particleDensityRatio, const MInt random = 1, const MInt addMode = 0,
                    const MFloat* velocity = nullptr, const MFloat* position = nullptr);

  MInt loadParticleRestartFile();

  void writePartData();
  void crankAngleSolutionOutput();

  // Particle sizes for exchange
  template <class LPTParticle>
  MInt elemPerP() {
    return LPTParticle::s_floatElements + LPTBase<nDim>::s_floatElements;
  }
  template <class LPTParticle>
  MInt intElemPerP() {
    return LPTParticle::s_intElements + LPTBase<nDim>::s_intElements;
  }
  template <class LPTParticle>
  MInt bufSize() {
    return m_exchangeBufferSize * elemPerP<LPTParticle>() + 1;
  }
  template <class LPTParticle>
  MInt intBufSize() {
    return m_exchangeBufferSize * intElemPerP<LPTParticle>() + 1;
  }

  void exchangeParticles(const MBool collOnly, const MInt offset, const MBool allowNonLeaf = false);
  template <class LPTParticle>
  void exchangeParticles_(const MBool collOnly, const MInt offset, const MBool allowNonLeaf = false);
  template <class LPTParticle>
  void sendAndReceiveParticles(const MBool allowNonLeaf);
  template <MBool allNeighbors, class LPTParticle>
  void sendParticles();
  template <MBool allNeighbors>
  void receiveParticles();
  template <MBool allNeighbors, class LPTParticle>
  void receiveParticles_();
  void exchangeOffset(MInt noParticles);
  void exchangeSourceTerms();
  void sendSourceTerms();
  void receiveSourceTerms();
  void sendFlowField();
  void receiveFlowField();
  void waitForSendReqs();
  void sendVelocitySlopes();
  void receiveVelocitySlopes();

  void checkParticles();
  void checkCells();
  void resetSourceTerms(const MInt);

  template <class T, class Pred>
  void own_sort(T& c, Pred& p) {
    using tag = typename std::iterator_traits<typename T::iterator>::iterator_category;
    sort_impl_(c, p, tag());
  }

  template <class T, class Pred>
  void sort_impl_(T& c, Pred& p, std::random_access_iterator_tag /*unused*/) {
    std::sort(c.begin(), c.end(), p);
  }

  template <class T, class Pred>
  void own_remove_if(T& c, Pred& p) {
    using tag = typename std::iterator_traits<typename T::iterator>::iterator_category;
    remove_if_impl_(c, p, tag());
  }

  template <class T, class Pred>
  void remove_if_impl_(T& c, Pred& p, std::random_access_iterator_tag /*unused*/) {
    c.erase(std::remove_if(c.begin(), c.end(), p), c.end());
  }

  MBool smallParticle(const LPTSpherical<nDim>& particle) { return (particle.m_diameter < m_sizeLimit); }

  inline std::mt19937_64& randomSpawn(const MInt calls) {
    ASSERT(domainId() == m_spawnDomainId, "");
    m_PRNGSpawnCount += calls;
    return m_PRNGSpawn;
  }
  inline std::mt19937_64& randomRespawn() {
    ASSERT(domainId() == m_respawnDomain, "");
    return m_PRNGRespawn;
  }
  inline std::mt19937_64& randomInit() {
    ASSERT(!m_restart, "");
    return m_PRNGInit;
  }

 private:
  MInt m_skewnessTimeStep = -1;

  MBool m_periodicBC = false;
  MBool m_engineSetup = false;
  std::array<MFloat, nDim> m_globDomainLength{};
  std::array<MFloat, nDim * 2> m_globBbox{0};

  MFloat m_sizeLimit{};
  MLong m_spawnSeed{};
  MBool m_activePrimaryBUp{};
  MBool m_activeSecondaryBUp{};

  MInt m_maxNoParticles{};

  MFloat m_xCutOff = -1000.0;

  MBool m_restart = false;

  MBool m_ellipsoids = false;
  MInt m_ellipsoidRandomOrientation = 1;
  MLong m_ellipsoidRandomOrientationSeed{};

  MInt m_initializationMethod = 0;
  MFloat m_FrMag = 0.0;
  MInt m_dragModelType = 1;
  MInt m_liftModelType = 0;
  MInt m_torqueModelType = 0;
  MInt m_motionEquationType = 0;

  // respawn properties
  MBool m_respawn = false;
  MFloat m_respawnPlane{};
  MInt m_respawnDomain{};
  MInt m_noRespawnDomains{};
  MInt* m_respawnDomainRanks{};
  MInt* m_respawnGlobalDomainOffsets{};
  std::vector<MInt> m_respawnCells;

  MBool m_wallCollisions = true;
  MInt m_maxNoBndryCells = 0;
  MInt m_collisions = 0;
  MInt m_noOuterBndryCells = 0;

  MInt m_exchangeBufferSize = 1000;

  MFloat m_cfl = 1.0;

  std::unique_ptr<SprayModel<nDim>> m_sprayModel;
  static constexpr MInt m_addInjDataCnt = 4;

  // particle spawn properties:
  MBool m_spawnParticles = false;
  MInt m_spawnCellId = -1;
  MInt m_spawnDomainId = -1;
  MFloat m_spawnParticlesConeAngle = 0.0;
  MFloat m_spawnDir[nDim]{};
  MFloat m_spawnDiameter = 0.00001;
  MInt m_spawnParticlesCount = 1;
  MFloat m_spawnCoord[nDim]{};
  MFloat m_spawnVelocity = 0.0;
  MInt m_spawnEmittDist = string2enum("PART_EMITT_DIST_UNIFORM");
  MFloat m_spawnDistSigmaCoeff = 1.0;
  MInt m_sprayAngleModel = 0;

  std::unique_ptr<MaterialState<nDim>> m_material;

  MLong m_addedParticle = 0;

  std::multimap<MInt, LPTSpherical<nDim>*> m_cellToPartMap;
  std::multimap<MInt, LPTEllipsoidal<nDim>*> m_cellToEllipsMap;

  MInt m_adaptationLevel = 0;

  // DLB reinitialization stage
  MInt m_loadBalancingReinitStage = -1;

  MInt m_shadowWidth{};
  MInt m_innerBound = 2;
  MBool m_forceAdaptation = false;

  // particle functions (placed in lpt_props.h)
  void readModelProps();
  void readSpawnProps();
  void readEllipsoidProps();
  void readMomentumCouplingProps();
  void findSpawnCellId();

  MBool fiveStageSolutionStep();
  MBool oneStageSolutionStep();

  void initMPI(const MBool = true);
  void initialCondition();
  void initSummary();
  void initCollector();
  void initBndryCells();
  void initParticleVector();
  void initGridProperties();
  void initParticleLog();

  void writeParticleLog();

  void computeBoundingBox(MFloat* globalBbox) {
    // if has cutoff, find solver specific bouding box, else use bounding box of raw grid
    if(Context::propertyExists("cutOffCoordinates", m_solverId)
       || Context::propertyExists("cutOffDirections", m_solverId)) {
      MFloat bboxLocal[nDim * 2];
      for(MInt i = 0; i < nDim; i++) {
        bboxLocal[i] = std::numeric_limits<MFloat>::max();
        bboxLocal[nDim + i] = std::numeric_limits<MFloat>::lowest();
      }
      // Find local max and min values
      for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
        for(MInt i = 0; i < nDim; i++) {
          bboxLocal[i] =
              mMin(bboxLocal[i],
                   grid().tree().coordinate(cellId, i) - F1B2 * grid().cellLengthAtLevel(grid().tree().level(cellId)));
          bboxLocal[nDim + i] =
              mMax(bboxLocal[nDim + i],
                   grid().tree().coordinate(cellId, i) + F1B2 * grid().cellLengthAtLevel(grid().tree().level(cellId)));
        }
      }
      for(MInt i = 0; i < 2 * nDim; i++)
        globalBbox[i] = bboxLocal[i];
      // reduce over all domains to find the global min and max values
      MPI_Allreduce(MPI_IN_PLACE, &globalBbox[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                    "m_bbox[0]");
      MPI_Allreduce(MPI_IN_PLACE, &globalBbox[nDim], nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "m_bbox[nDim]");
    } else {
      for(MInt i = 0; i < 2 * nDim; i++) {
        globalBbox[i] = grid().raw().globalBoundingBox()[i];
      }
    }
  }

  void addBndryCell(const MInt, MFloat*, MFloat*, MFloatScratchSpace&, MFloat*);

  // void addBndryCell(const MInt cellId, std::array<MFloat, nDim + 1>& nbndryData, std::array<MFloat, nDim>& surfaceN,
  //                  std::array<std::array<MFloat, nDim>, nDim + 1>& surfaceC, std::array<MFloat, nDim>& surfaceV);


  void packParticles(std::vector<LPTSpherical<nDim>*>& particlesToSend, MInt* intBuffer, MFloat* floatBuffer,
                     std::vector<MInt>);
  void packParticles(std::vector<LPTEllipsoidal<nDim>*>& particlesToSend, MInt* intBuffer, MFloat* floatBuffer,
                     std::vector<MInt>);

  template <class LPTParticle, MBool t_search>
  void unpackParticles(const MInt num, const MInt* intBuffer, const MFloat* floatBuffer, const MInt domainId = -1,
                       const MBool allowNonLeaf = false);
  void unpackParticle(LPTSpherical<nDim>& thisParticle, const MInt* intBuffer, MInt& z, const MFloat* floatBuffer,
                      MInt& h, MInt& step, MInt id);
  void unpackParticle(LPTEllipsoidal<nDim>& thisParticle, const MInt* intBuffer, MInt& z, const MFloat* floatBuffer,
                      MInt& h, MInt& step, MInt id);

  void broadcastInjected(const MUint prevNumPart);
  void recvInjected();

  void sendInjected(const MUint prevNumPart);

  void updateExchangeCells();

  void motionEquation(const MInt offset);
  void motionEquationEllipsoids(const MInt offset);
  void evaporation(const MInt offset);
  void coupling(MInt offset);
  void advanceParticles(const MInt offset);
  void perCellStats();
  void countParticlesInCells();
  void reduceParticles();
  void mergeParticles(LPTSpherical<nDim>* partA, LPTSpherical<nDim>* partB);
  void spawnTimeStep();
  void sprayInjection();
  // void collision();
  void removeInvalidParticles(const MBool);

  void wallCollision();
  void particleWallCollisionStep(const MInt);
  void writeCollData();

  void setRegridTrigger();
  void checkRegridTrigger();
  void recvRegridTrigger();

  void initializeTimers();

 public:
  MFloat m_sutherlandConstant = -1;
  MFloat m_sutherlandPlusOne = -1;

  MFloat time() const override { return m_time; }
  MInt noVariables() const override {
    // noParticlesInCell
    return 1;
  }

  MBool m_nonDimensional = false;

  MInt m_noSolutionSteps = 1;
  MInt m_solutionStep = 0;
  MInt m_restartTimeStep{};

  //-------------- Pass-Through accesors to the grid proxy ------------------------------------------

  MInt c_noCells() const { return grid().tree().size(); }

  MInt c_childId(const MInt cellId, const MInt pos) const { return grid().tree().child(cellId, pos); }

  MInt c_noChildren(const MInt cellId) const { return grid().tree().noChildren(cellId); }

  MInt c_parentId(const MInt cellId) const { return grid().tree().parent(cellId); }

  MFloat c_cellLengthAtLevel(const MInt level) const { return grid().cellLengthAtLevel(level); }

  MFloat c_cellLengthAtCell(const MInt cellId) const { return c_cellLengthAtLevel(c_level(cellId)); }

  MInt c_level(const MInt cellId) const { return grid().tree().level(cellId); }

  MInt a_level(const MInt cellId) const { return c_level(cellId); }

  MFloat c_cellVolume(const MInt cellId) const { return grid().cellVolumeAtLevel(grid().tree().level(cellId)); }

  MFloat c_coordinate(const MInt cellId, const MInt dir) const { return grid().tree().coordinate(cellId, dir); }

  MBool c_isLeafCell(const MInt cellId) const { return grid().tree().isLeafCell(cellId); }

  MLong c_globalId(const MInt cellId) const { return grid().tree().globalId(cellId); }

  // override accessors as already defined in the solver!
  MInt noInternalCells() const override { return grid().noInternalCells(); }
  MInt domainId() const override { return grid().domainId(); }
  MInt noDomains() const override { return grid().noDomains(); }

  MInt c_hasNeighbor(const MInt cellId, const MInt dir) const { return grid().tree().hasNeighbor(cellId, dir); }

  MInt c_neighborId(const MInt cellId, const MInt dir) const { return grid().tree().neighbor(cellId, dir); }

  inline MBool a_validCell(const MInt cellId) {
    if(cellId < 0 || cellId > a_noCells()) {
      return false;
    }
    return true;
  }


  //-------------- empty functions for post-processing ------------------------------------------

  void loadSampleVariables(MInt timeStep) {
    std::cerr << "loadSampleVariables DgCartesianSolver " << timeStep << std::endl;
  };
  void getSampleVariables(MInt cellId, const MFloat*& vars) {
    std::cerr << "getSampleVariables LPT " << cellId << " " << vars << std::endl;
  };
  void getSampleVariables(const MInt cellId, std::vector<MFloat>& /*vars*/) {
    std::cerr << "getSampleVariables LPT " << cellId << std::endl;
  };
  void calculateHeatRelease() { std::cerr << "calculateHeatRelease DgCartesianSolver " << std::endl; }
  void getHeatRelease(MFloat*& heatRelease) {
    std::cerr << "getHeatRelease DgCartesianSolver " << heatRelease << std::endl;
  }
};

/// \brief Static indices for accessing primitive variables
/// in nDim spatial dimensions
template <MInt nDim>
struct LPT<nDim>::PrimitiveVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, 3> VV = {0, 1, 2};
  static constexpr MInt RHO = nDim;
  static constexpr MInt P = nDim + 1;
  static constexpr MInt T = nDim + 2;

  static constexpr MInt m_noVars = nDim + 3;
  constexpr MInt noVars() const { return m_noVars; };
};

#endif
