// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDCONTROLLER_H_
#define GRIDCONTROLLER_H_

#include <functional>
#include "COMM/mpioverride.h"
#include "COUPLER/coupling.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "cartesiangrid.h"
#include "cartesiangridproxy.h"
#include "partition.h"
#include "solver.h"

#include "LB/lbsolverdxqy.h"

namespace maia {
namespace grid {

// Create struct for easy timer identification
struct Timers {
  // Enum to store timer "names"
  enum {
    Controller,
    DLB,
    Trigger,
    Partition,
    Prepare,
    BalanceGrid,
    BalanceSolvers,
    BalanceCouplers,
    IO,

    Adaptation,
    MeshAdaptation,

    // Special enum value used to initialize timer array
    _count
  };
};

// Timers for each solver
struct SolverTimers {
  // Enum to store timer "names"
  enum {
    Total,
    DLB,
    // CancelMPI,
    Reset,
    // CommGlobalVars,

    CalcDataSizes,
    CalcDataSizesMpiBlocking,
    CalcDataSizesMpi,

    // LocalToGlobalIds
    Balance,

    Redistribute,
    RedistributeMpiBlocking,
    RedistributeMpi,

    BalancePre,

    SetData,
    SetDataMpiBlocking,
    SetDataMpi,

    // GlobalToLocalIds
    BalancePost,
    DlbOther,

    Adaptation,

    // Special enum value used to initialize timer array
    _count
  };
};


/// Grid controller manages adaptation and load balancing in multi-solver environment.
template <MInt nDim>
class Controller {
 public:
  using Grid = CartesianGrid<nDim>;
  using Cell = typename CartesianGrid<nDim>::Cell;
  using GridProxy = typename maia::grid::Proxy<nDim>;

  // Constructor/destructor
  Controller(Grid* grid_, std::vector<std::unique_ptr<Solver>>* solvers,
             std::vector<std::unique_ptr<Coupling>>* couplers);
  ~Controller() {
    if(m_timersInitialized) { // Timers are not created if there is no Cartesian grid
      RECORD_TIMER_STOP(m_timers[Timers::Controller]);
    }
  }

  // Dynamic load balancing
  MBool balance(const MBool force = false, const MBool finalTimeStep = false, const MBool adaptation = false);
  MBool adaptation(const MBool force = false);
  void writeRestartFile(const MBool, const MBool);
  void printDomainStatistics(const MString& status = "");
  void savePartitionFile() {
    if(!m_balance || !gridb().wasBalancedAtLeastOnce()) {
      return;
    }
    gridb().savePartitionFile();
  }

  // Update the partition workloads in the grid file using given weights
  MBool updateGridPartitionWorkloads();

  // Cast of adaptation interval needed for LB
  void castAdaptationIntervalToMultipleOfCoarsestTimeStep(MInt maxLevel, MInt maxUniformRefinementLevel);

  void logTimerStatistics(const MString&);
  MBool isDlbTimeStep();
  MBool isAdaptationTimeStep();

 private:
  // Grid accessors
  // Accessor to grid is named "gridb" to indicate that this is the base grid and not a
  // solver-specific grid
  Grid& gridb() { return *m_grid; }
  const Grid& gridb() const { return *m_grid; }
  // constexpr const GridProxy& gridProxy() const { return m_gridProxy; }
  // GridProxy& gridProxy() { return m_gridProxy; }

  MInt domainId() { return gridb().domainId(); }
  MInt noDomains() { return gridb().noDomains(); }

  // Solver accessors
  Solver& solver(const MInt solverId) { return *m_solvers->at(solverId); }
  const Solver& solver(const MInt solverId) const { return (const Solver&)m_solvers->at(solverId); }
  MInt noSolvers() const { return m_solvers->size(); }
  MInt solverLocalRootDomain(Solver* const solver);

  // Coupler accessors
  Coupling& coupler(const MInt couplerId) { return *m_couplers->at(couplerId); }
  const Coupling& coupler(const MInt couplerId) const { return *m_couplers->at(couplerId); }
  MInt noCouplers() const { return m_couplers->size(); }

  void initDlbProperties();
  void initTimers();

  MBool needLoadBalancing(const MFloat localRunTime, MFloat* const loads, MFloat& imbalance);

  void partition(MLong* partitionCellOffsets, MLong* globalIdOffsets, const MBool onlyPartitionOffsets);
  MBool loadBalancingPartition(const MFloat* loads, const MFloat imbalance, MLong* const partitionCellOffsets,
                               MLong* const globalIdOffsets);

  void updateWeightsAndWorkloads(const std::vector<MFloat>& weights, const MBool restoreDefaultWeights);
  void getSpecifiedSolverWeights(std::vector<MFloat>& weights);

  void accumulateCellWeights();

  void storeTimings();
  void storeLoadsAndWeights(const MFloat* const loads, const MInt noLoadTypes, const MInt* const loadQuantities,
                            const MFloat domainWeight, const MFloat* const weights);

  void loadBalancingCalcNewGlobalOffsets(const MLong* const oldPartitionCellOffsets,
                                         const MLong* const newPartitionCellOffsets,
                                         MLong* const globalOffsets);

  void loadBalancingCalcNoCellsToSend(const MLong* const offsets,
                                      MInt* const noCellsToSendByDomain,
                                      MInt* const noCellsToReceiveByDomain,
                                      MInt* const sortedCellId,
                                      MInt* const bufferIdToCellId);

  void computeWeights(const MFloat* loads, const MFloat domainWeight, std::vector<MFloat>& weights);
  MInt globalNoLoadTypes();
  void getLoadQuantities(MInt* const loadQuantities);

  void estimateParameters(MInt m, MInt n, const MFloat* const A, const MFloat* const b, MFloat* const x);

  void communicateGlobalSolverVars(Solver* const solver);

  void determineDataSizesDlb(const MInt solverId, const MInt mode, const MInt* const noCellsToSend,
                             const MInt* const bufferIdToCellId, std::vector<std::vector<MInt>>& sendSizeVector,
                             std::vector<std::vector<MInt>>& recvSizeVector);

  void redistributeDataDlb(const MInt id, const MInt mode, std::vector<MInt>& sendSizeVector,
                           std::vector<MInt>& recvSizeVector, const MInt* const bufferIdToCellId, const MInt noCells,
                           std::vector<MInt*>& intDataRecv, std::vector<MLong*>& longDataRecv,
                           std::vector<MFloat*>& floatDataRecv, std::vector<MInt>& dataTypes);

  void setDataDlb(const MInt solverId, const MInt mode, std::vector<MInt*>& intDataRecv,
                  std::vector<MLong*>& longDataRecv, std::vector<MFloat*>& floatDataRecv, std::vector<MInt>& dataTypes,
                  const MBool freeMemory);

  void resetAllTimer();

  // Data members
  Grid* m_grid;
  // GridProxy& m_gridProxy;
  const std::vector<std::unique_ptr<Solver>>* m_solvers;
  const std::vector<std::function<void(const MInt)>> m_refineCellSolver;
  const std::vector<std::function<void(const MInt)>> m_removeChildsSolver;
  const std::vector<std::function<void(const MInt, const MInt)>> m_swapProxySolver;
  const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>> m_cellOutside;
  const std::vector<std::function<void()>> m_resizeGridMapSolver;
  const std::vector<std::function<void(const MInt)>> m_removeCellSolver;

  const std::vector<std::unique_ptr<Coupling>>* m_couplers;

  // Mesh adaptation
  MBool m_adaptation = false;
  MInt m_adaptationInterval = 0;
  MInt m_adaptationStart = -1;
  MInt m_adaptationStop = -1;
  MInt m_lastAdaptationTimeStep = 0;
  MInt m_noAdaptations = 0;
  MInt m_noMaxAdaptations = -1;

  // Sensors
  MInt m_noSensors{};

  // Dynamic load balancing
  MBool m_balance = false;
  // enforce balance even for low imbalance
  MInt m_forceBalance = 0;
  // Performance output switch if balance is not enabled
  MBool m_performanceOutput = false;
  // Debug mode
  MBool m_debugBalance = false;

  // Load balancing mode
  MInt m_loadBalancingMode = -1;
  // Load balancing timer mode
  MInt m_loadBalancingTimerMode = 0;

  // Main properties for dynamic load balancing
  MInt m_loadBalancingInterval = -1;        // Load balancing interval
  MInt m_loadBalancingOffset = 0;           // Load balancing offset
  MInt m_loadBalancingStartTimeStep = 0;    // Load balancing start timestep
  MInt m_loadBalancingStopTimeStep = -1;    // Load balancing stop timestep
  MInt m_loadBalancingTimerStartOffset = 0; // Offset for DLB timings
  MBool m_forceLoadBalancing = false;       // Switch to force DLB
  MBool m_testDynamicLoadBalancing = false; // Switch for DLB testcases
  // Switch for DLB testcases with partition level shift and testing of partition cell update
  MBool m_testUpdatePartitionCells = false;
  // Thresholds for testing/forcing partition cell updates
  MInt m_testUpdatePartCellsOffspringThreshold = 10;
  MFloat m_testUpdatePartCellsWorkloadThreshold = 100.0;
  // starting with large value, to trigger first balance
  MInt m_nAdaptationsSinceBalance = 9999999;

  MBool m_syncTimeSteps = false;
  MBool m_syncTimerSteps = true;

  // Number of performed DLB steps
  MInt m_dlbStep = 0;
  // Last time step at which DLB was performed
  MInt m_lastLoadBalancingTimeStep = 0;

  // Last time step at which timers were reset
  MInt m_dlbLastResetTimeStep = 0;

  // Imbalance threshold to trigger DLB (default 5%)
  MFloat m_dlbImbalanceThreshold = 0.05;
  // Performance variation threshold to skip a DLB step if the timings on a domain are not reliable
  MFloat m_maxPerformanceVarThreshold = 0.15;
  // DLB partition method to use
  MInt m_dlbPartitionMethod = DLB_PARTITION_DEFAULT;
  // Update the partition cells before load balancing, i.e., introduce new or remove existing
  // partition cells on higher levels depending on the size/weight of the local subtree (i.e.
  // increase/decrease partition level shifts)
  MBool m_dlbUpdatePartitionCells = false;
  // Indicates if a grid file with updated partition cells should be written with the next restart
  MBool m_saveGridNewPartitionCells = false;

  // Additional parameters for DLB_PARTITION_SHIFT_OFFSETS method
  MBool m_dlbSmoothGlobalShifts = true;
  // Number of intermediate DLB steps for which the offsets are only shifted locally
  MInt m_dlbNoLocalShifts = 0;
  // Number of final DLB steps for which the offsets are only shifted locally (after reverting to
  // the best found configuration until then)
  MInt m_dlbNoFinalLocalShifts = 0;
  // Limit for the maximum relative workload on a single domain
  MFloat m_dlbMaxWorkloadLimit = 1.5;

  // DLB timer records at previous timestep
  MFloat m_dlbPreviousLocalRunTime = 0.0;
  MFloat m_dlbPreviousLocalIdleTime = 0.0;

  // Timings of time steps since last DLB step (for DLB imbalance evaluation)
  std::vector<MFloat> m_dlbTimings;

  // Switch to output timings
  MBool m_outputDlbTimings = false;
  // Storage for timings of all timesteps (for performance evaluations)
  std::vector<MInt> m_dlbTimeStepsAll;
  std::vector<MFloat> m_dlbRunTimeAll;
  std::vector<MFloat> m_dlbIdleTimeAll;

  // Previous load distribution for output as a histogram, plus additional information
  std::vector<MInt> m_previousLoadBins{};
  std::array<MFloat, 5> m_previousLoadInfo{};
  MInt m_previousLoadInfoStep = -1;

  // Domain weights for DLB (i.e. processing capacities)
  std::vector<MFloat> m_domainWeights;
  // Switch to enable use of domain weight ratios during shifting of offsets
  MBool m_useDomainFactor = false;

  // Store last direction in which each offset was shifted
  std::vector<MInt> m_lastOffsetShiftDirection{};

  // DLB: store best time per time step and the corresponding 'optimal' partition-cell
  // offset (revert to this partitioning after a finite number of load balancing
  // steps to get the best performance in a non-adaptive simulation)
  MFloat m_timePerStepTotal = std::numeric_limits<MFloat>::max();
  MLong m_optPartitionCellOffsetTotal = -1;
  MFloat m_imbalance = std::numeric_limits<MFloat>::max();

  // balance after loadBalancingOffset timeSteps after any adaptation
  MBool m_balanceAfterAdaptation = false;
  MInt m_balanceAdaptationInterval = 1;

  // use DLB weights at restart
  MBool m_dlbRestartWeights = false;
  MFloat* m_dlbLastWeights = nullptr;
  MFloat* m_dlbStaticWeights = nullptr;
  MInt m_dlbStaticWeightMode = -1;

  // Timers
  // Status of timer initialization
  MBool m_timersInitialized = false;
  // Timer group which holds all controller-wide timers
  MInt m_timerGroup = -1;
  // Stores all controller-wide timers
  std::array<MInt, Timers::_count> m_timers{};
  // Timer groups for solver (and coupler) timers
  std::vector<MInt> m_solverTimerGroups{};
  // Stores all solver (and coupler) timers
  std::vector<std::array<MInt, SolverTimers::_count>> m_solverTimers{};

  // Create temporary function pointer vectors to create const class members (lennart)
  std::vector<std::function<void(const MInt)>> refineCellVec() {
    std::vector<std::function<void(const MInt)>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::refineCell, &solver(i), std::placeholders::_1);
    }
    return vec;
  }
  std::vector<std::function<void(const MInt)>> removeChildsVec() {
    std::vector<std::function<void(const MInt)>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::removeChilds, &solver(i), std::placeholders::_1);
    }
    return vec;
  }
  std::vector<std::function<void(const MInt, const MInt)>> swapProxyVec() {
    std::vector<std::function<void(const MInt, const MInt)>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::swapProxy, &solver(i), std::placeholders::_1, std::placeholders::_2);
    }
    return vec;
  }
  std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>> cellOutsideVec() {
    std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::cellOutside, &solver(i), std::placeholders::_1, std::placeholders::_2,
                         std::placeholders::_3);
    }
    return vec;
  }
  std::vector<std::function<void()>> resizeGridMapVec() {
    std::vector<std::function<void()>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::resizeGridMap, &solver(i));
    }
    return vec;
  }
  std::vector<std::function<void(const MInt)>> removeCellVec() {
    std::vector<std::function<void(const MInt)>> vec(noSolvers());
    for(MInt i = 0; i < noSolvers(); i++) {
      vec[i] = std::bind(&Solver::removeCell, &solver(i), std::placeholders::_1);
    }
    return vec;
  }

  // Restart
  MString m_currentGridFileName;
  MBool m_useNonSpecifiedRestartFile = false;
  MString m_outputDir;
  MInt* m_recalcIds;
};


template <MInt nDim>
Controller<nDim>::Controller(Grid* grid_, std::vector<std::unique_ptr<Solver>>* solvers,
                             std::vector<std::unique_ptr<Coupling>>* couplers)
  : m_grid(grid_),
    m_solvers(solvers),
    m_refineCellSolver(refineCellVec()),
    m_removeChildsSolver(removeChildsVec()),
    m_swapProxySolver(swapProxyVec()),
    m_cellOutside(cellOutsideVec()),
    m_resizeGridMapSolver(resizeGridMapVec()),
    m_removeCellSolver(removeCellVec()),
    m_couplers(couplers) {
  if(m_grid == nullptr) {
    return;
  }

  // Read properties for adaptive mesh refinement

  /*! \property
  \page propertiesAMR
  \section adaptation
  <code>MBool m_adaptation </code>\n
  default = <code>false</code>\n
  Triggers adaptive mesh refinement
  Possible values are:
  <ul>
  <li> true </li>
  <li> false </li>
  </ul>
  Keywords: <i> Mesh Adaptation, Cartesian Grid</i>
  */
  m_adaptation = false;
  m_adaptation = Context::getBasicProperty<MBool>("adaptation", AT_, &m_adaptation);

  /*! \property
  \page propertiesAMR
  \section adaptationInterval
  <code>MInt m_adaptationInterval </code>\n
  default = <code>0</code>\n
  Number of timesteps between mesh adaptations.
  Possible values are:
  <ul>
  <li> integer >=0 </li>
  </ul>
  Keywords: <i> Mesh Adaptation, Cartesian Grid</i>
  */
  m_adaptationInterval = 0;
  m_adaptationInterval = Context::getBasicProperty<MInt>("adaptationInterval", AT_, &m_adaptationInterval);

  /*! \property
  \page propertiesAMR
  \section adaptationStart
  <code>MInt m_adaptationStart </code>\n
  default = <code>0</code>\n
  First possible time step with adaptation,
  i.e. adaptation is skipped before this time step!
  Possible values are:
  <ul>
  <li> integer >=0 </li>
  </ul>
  Keywords: <i> Mesh Adaptation, Cartesian Grid</i>
  */
  m_adaptationStart = 0;
  m_adaptationStart = Context::getBasicProperty<MInt>("adaptationStart", AT_, &m_adaptationStart);

  /*! \property
  \page propertiesAMR
  \section adaptationStop
  <code>MInt m_adaptationStop </code>\n
  default = <code>max. Int</code>\n
  Last possible time step with adaptation,
  i.e. adaptation is skipped after this time step!
  Possible values are:
  <ul>
  <li> integer >=0 </li>
  </ul>
  Keywords: <i> Mesh Adaptation, Cartesian Grid</i>
  */
  m_adaptationStop = std::numeric_limits<MInt>::max();
  m_adaptationStop = Context::getBasicProperty<MInt>("adaptationStop", AT_, &m_adaptationStop);

  /*! \property
  \page propertiesAMR
  \section noMaxAdaptations
  <code>MInt m_noMaxAdaptations </code>\n
  default = <code>max. Int</code>\n
  Maximum number of adaptation calls,
  i.e. adaptation is skipped after this number of adaptation calls!
  Possible values are:
  <ul>
  <li> integer >=0 </li>
  </ul>
  Keywords: <i> Mesh Adaptation, Cartesian Grid</i>
  */
  m_noMaxAdaptations = std::numeric_limits<MInt>::max();
  m_noMaxAdaptations = Context::getBasicProperty<MInt>("noMaxAdaptations", AT_, &m_noMaxAdaptations);

  m_lastAdaptationTimeStep = g_restartTimeStep;
  m_dlbLastResetTimeStep = g_restartTimeStep;

  // Read properties for write-restart-file:
  const MInt maxNoCell = gridb().maxNoCells();
  m_useNonSpecifiedRestartFile =
      Context::getBasicProperty<MBool>("useNonSpecifiedRestartFile", AT_, &m_useNonSpecifiedRestartFile);

  m_currentGridFileName = gridb().m_gridInputFileName;
  m_outputDir = Context::getBasicProperty<MString>("outputDir", AT_);

  m_recalcIds = (MInt*)nullptr;
  mAlloc(m_recalcIds, maxNoCell, "m_recalcIds", -1, AT_);
  for(MInt i = 0; i < maxNoCell; i++) {
    m_recalcIds[i] = i;
  }

  initDlbProperties();

  initTimers();

  m_syncTimeSteps = Context::getBasicProperty<MBool>("syncTimeSteps", AT_, &m_syncTimeSteps);
}


/// Read Dynamic Load Balancing properties
template <MInt nDim>
void Controller<nDim>::initDlbProperties() {
  m_loadBalancingInterval = 0;
  if(Context::propertyExists("onlineRestartInterval", -1)) {
    if(domainId() == 0) {
      std::cerr << "Property 'onlineRestartInterval' is deprecated, please rename it to "
                   "'loadBalancingInterval'."
                << std::endl;
    }
    m_loadBalancingInterval = Context::getBasicProperty<MInt>("onlineRestartInterval", AT_, &m_loadBalancingInterval);
  }

  /*! \property
    \page propertiesDLB
    \section loadBalancingInterval
    <code>MInt CartesianGrid::m_loadBalancingInterval</code>\n
    default = <code>0</code> (disabled)\n \n
    Time step interval in which dynamic load balancing is performed.\n
    possible values are:
    <ul>
    <li>any positive integer</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingInterval = Context::getBasicProperty<MInt>("loadBalancingInterval", AT_, &m_loadBalancingInterval);

  /*! \property
    \page propertiesDLB
    \section balance
    <code>MBool Controller::m_balance</code>\n
    default = <code>(loadBalancingInterval > 0)</code>\n \n
    Enables dynamic load balancing.
    Is activated automatically when a positive load balancing interval is specified.\n
    <ul>
    <li>0: disabled</li>
    <li>1: enabled</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_balance = (m_loadBalancingInterval > 0);
  m_balance = Context::getBasicProperty<MBool>("balance", AT_, &m_balance);

  m_forceBalance = 0;
  // 0: force never
  // 1: force once
  // 2: force always
  m_forceBalance = Context::getBasicProperty<MInt>("forceBalance", AT_, &m_forceBalance);

  /*! \property
    \page propertiesDLB
  \section balanceAfterAdaptation
  <code>MFloat GridController::m_balanceAfterAdaptation </code>\n
  default = <code>m_adaptation</code>\n \n
  Performance a balance after the balanceOffset of timeSteps after each adaptation.
  Meaning that after each adaptation after loadBalancingTimerStartOffset, the time are resetted and
  then after the loadBalancingOffset after the adaptation a balance is performed!
  <ul> <li>true/false</li> </ul>
  Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, TRIGGER, GLOBAL</i>
  */
  m_balanceAfterAdaptation = (m_adaptation && m_balance);
  m_balanceAfterAdaptation = Context::getBasicProperty<MBool>("balanceAfterAdaptation", AT_, &m_balanceAfterAdaptation);

  /*! \property
    \page propertiesDLB
  \section balanceAfterAdaptation
  <code>MInt MAIAGridController::m_balanceAdaptationInterval </code>\n
  default = <code>1</code>\n \n
  Interval after how many adaptations a balance is triggered, default is 1 meaning after every
  Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, TRIGGER, GLOBAL</i>
  */
  m_balanceAdaptationInterval = 1;
  m_balanceAdaptationInterval =
      Context::getBasicProperty<MInt>("balanceAfterAdaptationInterval", AT_, &m_balanceAdaptationInterval);

  m_loadBalancingOffset = 0;
  m_loadBalancingOffset = Context::getBasicProperty<MInt>("loadBalancingOffset", AT_, &m_loadBalancingOffset);

  const MBool balanceOnlyAfterAdapt = (m_balanceAfterAdaptation && m_loadBalancingInterval <= 0);
  if(balanceOnlyAfterAdapt) {
    m_loadBalancingInterval = std::numeric_limits<MInt>::max();
  }
  m_balance = m_balance || m_balanceAfterAdaptation;

  // Return if not enabled
  if(!m_balance) {
    // Note: always display imbalance even though balancing is not enabled
    m_performanceOutput = true;
    m_performanceOutput = Context::getBasicProperty<MBool>("performanceOutput", AT_, &m_performanceOutput);
    m_loadBalancingMode = 1;
    // Determine interval for performance output as power of 10
    MInt interval = std::max((MInt)pow(10.0, std::floor(std::log10((MFloat)g_timeSteps / 10.0))), 1);
    if(g_timeSteps / interval > 20.0) { // Increase interval to limit total number of evaluations
      interval = interval * 5;          // Now interval is in {10,50,100,500,1000,...}
    }

    m_loadBalancingInterval = std::max(interval, 10);
    // Return if balance and performanceOutput are both disabled
    if(!m_performanceOutput) {
      m_balance = false;
      m_log << "Performance output: disabled." << std::endl;
      return;
    } else {
      m_balance = true; // Set s.t. balance() does not immediately return
      m_log << "Performance output: enabled - every " << m_loadBalancingInterval << " time steps." << std::endl;
    }
    m_log << "Dynamic load balancing disabled." << std::endl;
    return;
  }

  /*! \property
    \page propertiesDLB
    \section loadBalancingMode
    <code>MInt Controller::m_loadBalancingMode</code>\n
    default = <code>0</code>\n \n
    Dynamic load balancing mode.
    <ul>
    <li>0: continuous, i.e. periodically, using fixed weights specified in the solver</li>
    <li>1: finite number of DLB steps possible, e.g. for static workload distributions; estimation
    of computational weights based on runtime measurements</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingMode = 0;
  m_loadBalancingMode = Context::getBasicProperty<MInt>("loadBalancingMode", AT_, &m_loadBalancingMode);

  /*! \property
    \page propertiesDLB
    \section loadBalancingTimerMode
    <code>MInt Controller::m_loadBalancingTimerMode</code>\n
    default = <code>0</code>\n \n
    Load balancing timer mode for runtime measurements.\n
    Measure compute/idle times based on:
    <ul>
    <li>0: wall time using MPI_Wtime()</li>
    <li>1: process cpu time using clock_gettime(...)</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingTimerMode = 0;
  m_loadBalancingTimerMode = Context::getBasicProperty<MInt>("loadBalancingTimerMode", AT_, &m_loadBalancingTimerMode);

  /*! \property
    \page propertiesDLB
    \section loadBalancingStartTimeStep
    <code>MInt CartesianGrid::m_loadBalancingStartTimeStep</code>\n
    default = <code>0</code>\n \n
    Start time step for dynamic load balancing.\n
    possible values are:
    <ul>
    <li>any integer \>=0</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingStartTimeStep = 0;
  m_loadBalancingStartTimeStep =
      Context::getBasicProperty<MInt>("loadBalancingStartTimeStep", AT_, &m_loadBalancingStartTimeStep);

  /*! \property
    \page propertiesDLB
    \section loadBalancingStopTimeStep
    <code>MInt CartesianGrid::m_loadBalancingStopTimeStep</code>\n
    default = <code>-1</code>\n \n
    Stop time step for dynamic load balancing.\n
    possible values are:
    <ul>
    <li>any positive integer</li>
    <li>-1 : no stop time step</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingStopTimeStep = std::numeric_limits<MInt>::max();
  m_loadBalancingStopTimeStep =
      Context::getBasicProperty<MInt>("loadBalancingStopTimeStep", AT_, &m_loadBalancingStopTimeStep);

  /*! \property
    \page propertiesDLB
    \section loadBalancingTimerStartOffset
    <code>MInt CartesianGrid::m_loadBalancingTimerStartOffset</code>\n
    default = <code>20\% of load balancing interval</code>\n \n
    Number of time steps after each DLB step which are not considered for
    determining the performance.\n
    possible values are:
    <ul>
    <li>any integer \>=0 and \<=m_loadBalancingInterval</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_loadBalancingTimerStartOffset = std::floor(0.2 * m_loadBalancingInterval);
  // Set default timer start offset if balance is used only after adaptation
  // TODO labels:DLB case with balance interval + balance after adapt?
  if(balanceOnlyAfterAdapt) {
    if(m_loadBalancingOffset == 0) {
      m_loadBalancingTimerStartOffset = std::floor(0.2 * m_adaptationInterval);
    } else {
      m_loadBalancingTimerStartOffset = std::floor(0.2 * m_loadBalancingOffset);
    }
  }
  m_loadBalancingTimerStartOffset =
      Context::getBasicProperty<MInt>("loadBalancingTimerStartOffset", AT_, &m_loadBalancingTimerStartOffset);

  // TODO labels:DLB add checks for timerStartOffset for adaptation/+balancingoffset
  if(m_loadBalancingTimerStartOffset >= m_loadBalancingInterval) {
    TERMM(1, "DLB timerStartOffset = " + std::to_string(m_loadBalancingTimerStartOffset)
                 + " must be smaller than dlb-interval = " + std::to_string(m_loadBalancingInterval) + ".");
  }

  /*! \property
    \page propertiesDLB
    \section forceLoadBalancing
    <code>MBool CartesianGrid::m_forceLoadBalancing</code>\n
    default = <code>false</code>\n \n
    Force dynamic load balancing for testing (e.g. even if the partitioning
    did not change)\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_forceLoadBalancing = Context::getBasicProperty<MBool>("forceLoadBalancing", AT_, &m_forceLoadBalancing);

  /*! \property
    \page propertiesDLB
    \section testDynamicLoadBalancing
    <code>MBool CartesianGrid::m_testDynamicLoadBalancing</code>\n
    default = <code>false</code>\n \n
    Testcase switch for dynamic load balancing.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_testDynamicLoadBalancing = false;
  m_testDynamicLoadBalancing =
      Context::getBasicProperty<MBool>("testDynamicLoadBalancing", AT_, &m_testDynamicLoadBalancing);

  /*! \property
    \page propertiesDLB
    \section testUpdatePartitionCells
    <code>MBool GridController::m_testUpdatePartitionCells</code>\n
    default = <code>false</code>\n \n
    Testcase switch for updating partition cells before dynamic load balancing.\n
    Keywords: <i>PARALLEL, PARTITION, DYNAMIC, LOAD, BALANCING, GLOBAL, TESTING</i>
  */
  m_testUpdatePartitionCells = false;
  m_testUpdatePartitionCells =
      Context::getBasicProperty<MBool>("testUpdatePartitionCells", AT_, &m_testUpdatePartitionCells);

  if(m_testUpdatePartitionCells) {
    /*! \property
    \page propertiesDLB
      \section testUpdatePartitionCellsOffspringThreshold
      <code>MInt GridController::m_testUpdatePartitionCellsOffspringThreshold</code>\n
      default = <code>10</code>\n \n
      Offspring threshold for testing the update partition cells before dynamic load balancing.\n
      Keywords: <i>PARALLEL, PARTITION, DYNAMIC, LOAD, BALANCING, GLOBAL, TESTING</i>
    */
    m_testUpdatePartCellsOffspringThreshold = 10;
    m_testUpdatePartCellsOffspringThreshold = Context::getBasicProperty<MInt>(
        "testUpdatePartCellsOffspringThreshold", AT_, &m_testUpdatePartCellsOffspringThreshold);

    /*! \property
    \page propertiesDLB
      \section testUpdatePartitionCellsWorkloadThreshold
      <code>MFloat GridController::m_testUpdatePartitionCellsWorkloadThreshold</code>\n
      default = <code>100.0</code>\n \n
      Workload threshold for testing the update partition cells before dynamic load balancing.\n
      Keywords: <i>PARALLEL, PARTITION, DYNAMIC, LOAD, BALANCING, GLOBAL, TESTING</i>
    */
    m_testUpdatePartCellsWorkloadThreshold = 100.0;
    m_testUpdatePartCellsWorkloadThreshold = Context::getBasicProperty<MFloat>(
        "testUpdatePartCellsWorkloadThreshold", AT_, &m_testUpdatePartCellsWorkloadThreshold);

    m_log << "Testing partition cell update before DLB: offspringThreshold = "
          << m_testUpdatePartCellsOffspringThreshold
          << "; workloadThreshold = " << m_testUpdatePartCellsWorkloadThreshold << std::endl;
  }

  /*! \property
    \page propertiesDLB
    \section debugBalance
    <code>MBool CartesianGrid::m_debugBalance</code>\n
    default = <code>false</code>\n \n
    Enable additional debug output for dynamic load balancing.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_debugBalance = false;
  m_debugBalance = Context::getBasicProperty<MBool>("debugBalance", AT_, &m_debugBalance);

  if(m_testDynamicLoadBalancing) {
    // Force DLB for testcases and enable additional debug output
    m_forceLoadBalancing = true;
    m_debugBalance = true;
  }

  /*! \property
    \page propertiesDLB
    \section outputDlbTimings
    <code>MBool CartesianGrid::m_outputDlbTimings</code>\n
    default = <code>false</code>\n \n
    Switch to output dynamic load balancing timings.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_outputDlbTimings = false;
  m_outputDlbTimings = Context::getBasicProperty<MBool>("outputDlbTimings", AT_, &m_outputDlbTimings);

  m_log << "Dynamic load balancing activated (mode = " << m_loadBalancingMode
        << ", interval = " << m_loadBalancingInterval << ", startTimeStep = " << m_loadBalancingStartTimeStep
        << ", stopTimeStep = " << m_loadBalancingStopTimeStep
        << ", timerStartOffset = " << m_loadBalancingTimerStartOffset << ", force = " << m_forceLoadBalancing
        << ", test = " << m_testDynamicLoadBalancing << ")" << std::endl;

  /*! \property
    \page propertiesDLB
    \section dlbPartitionMethod
    <code>MInt CartesianGrid::m_dlbPartitionMethod</code>\n
    default = <code>DLB_PARTITION_DEFAULT</code>\n \n
    Dynamic load balancing partition method.\n
    possible values are:
    <ul>
    <li>DLB_PARTITION_DEFAULT : default method (combination of weighting and offset shifting
      methods)</li>
    <li>DLB_PARTITION_WEIGHT: weighting method</li>
    <li>DLB_PARTITION_SHIFT_OFFSETS: offset shifting method</li>
    <li>DLB_PARTITION_TEST: method for testing DLB reinitialization</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  if(Context::propertyExists("dlbPartitionMethod", -1)) {
    const MString dlbPartitionMethod = Context::getBasicProperty<MString>("dlbPartitionMethod", AT_);
    m_dlbPartitionMethod = string2enum(dlbPartitionMethod);
  }

  // Use weighting method for testing DLB
  if(m_testDynamicLoadBalancing) {
    m_dlbPartitionMethod = DLB_PARTITION_WEIGHT;
  }

  /*! \property
    \page propertiesDLB
    \section dlbUpdatePartitionCells
    <code>MBool GridController::m_dlbUpdatePartitionCells</code>\n
    default = <code>false</code>\n \n
    Switch for updating/refiltering of partition cells before dynamic load balancing.\n
    Keywords: <i>PARALLEL, PARTITION, DYNAMIC, LOAD, BALANCING, GLOBAL</i>
  */
  m_dlbUpdatePartitionCells = false;
  m_dlbUpdatePartitionCells =
      Context::getBasicProperty<MBool>("dlbUpdatePartitionCells", AT_, &m_dlbUpdatePartitionCells);

  m_useDomainFactor = false;
  m_useDomainFactor = Context::getBasicProperty<MBool>("dlbPartitionDomainFactor", AT_, &m_useDomainFactor);

  m_dlbSmoothGlobalShifts = true;
  m_dlbSmoothGlobalShifts = Context::getBasicProperty<MBool>("dlbSmoothGlobalShifts", AT_, &m_dlbSmoothGlobalShifts);

  m_dlbNoLocalShifts = 0;
  m_dlbNoLocalShifts = Context::getBasicProperty<MInt>("dlbNoLocalShifts", AT_, &m_dlbNoLocalShifts);

  m_dlbNoFinalLocalShifts = 0;
  m_dlbNoFinalLocalShifts = Context::getBasicProperty<MInt>("dlbNoFinalLocalShifts", AT_, &m_dlbNoFinalLocalShifts);

  // Note: use 0 or negative to disable max workload limitation
  m_dlbMaxWorkloadLimit = 1.5;
  m_dlbMaxWorkloadLimit = Context::getBasicProperty<MFloat>("dlbMaxWorkloadLimit", AT_, &m_dlbMaxWorkloadLimit);
  if(m_dlbMaxWorkloadLimit > 0.0 && m_dlbMaxWorkloadLimit < 1.0) {
    TERMM(1, "DLB: maximum workload limit needs to be > 1.0, is " + std::to_string(m_dlbMaxWorkloadLimit)
                 + "; or set to <= 0.0 to disable this feature.");
  }

  /*! \property
    \page propertiesDLB
    \section dlbImbalanceThreshold
    <code>MFloat GridController::m_dlbImbalanceThreshold</code>\n
    default = <code>0.05</code>\n \n
    Imbalance percentage threshold value for triggering dynamic load balancing.\n
    The imbalance percentage is defined as: I_% = (t_max - t_avg)/t_max * N/(N-1), with t_max and
    t_avg the maximum and the average run time among all N domains. It quantifies the 'badness' of
    the imbalance and gives the "percentage of resources available for parallelism that is wastet"
    (Reference: DeRose 2007, Detecting Application Load Imbalance on High End Massively Parallel
    Systems)\n
    possible values
    are: <ul> <li>any positive value \<1.0</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, TRIGGER, GLOBAL</i>
  */
  m_dlbImbalanceThreshold = 0.05;
  m_dlbImbalanceThreshold = Context::getBasicProperty<MFloat>("dlbImbalanceThreshold", AT_, &m_dlbImbalanceThreshold);

  /*! \property
    \page propertiesDLB
    \section maxPerformanceVarThreshold
    <code>MFloat GridController::m_maxPerformanceVarThreshold </code>\n
    default = <code>0.15</code>\n \n
    Performance variation threshold to skip a DLB step if the performance on a domain deviates over
    time and the timings cannot be reliably used for the estimation of computation weights etc.\n
    possible values are:
    <ul> <li>any positive value</li> </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, TRIGGER, GLOBAL</i>
  */
  m_maxPerformanceVarThreshold = 0.15;
  m_maxPerformanceVarThreshold =
      Context::getBasicProperty<MFloat>("maxPerformanceVarThreshold", AT_, &m_maxPerformanceVarThreshold);

  m_lastLoadBalancingTimeStep = m_loadBalancingStartTimeStep;

  if(m_loadBalancingMode == 1) {
    m_log << "Dynamic load balancing: using partition method #" << m_dlbPartitionMethod
          << "; imbalance percentage threshold: " << m_dlbImbalanceThreshold * 100.0
          << "%; maxPerfVarThreshold: " << m_maxPerformanceVarThreshold << std::endl;
    m_log << "DLB settings: domainFactor=" << m_useDomainFactor << "; noLocalShifts=" << m_dlbNoLocalShifts
          << "; noFinalLocalShifts=" << m_dlbNoFinalLocalShifts << "; maxWorkloadLimit=" << m_dlbMaxWorkloadLimit
          << "; smooth shifts=" << m_dlbSmoothGlobalShifts << std::endl;
  }

  m_dlbRestartWeights = false;
  if(m_loadBalancingMode == 0 && m_dlbPartitionMethod == DLB_PARTITION_WEIGHT) {
    m_dlbRestartWeights = true;
  }
  m_dlbRestartWeights = Context::getBasicProperty<MBool>("dlbRestartWeights", AT_, &m_dlbRestartWeights);

  m_dlbStaticWeights = nullptr;
  //-1: never apply static weights
  // 0: always apply static weights
  //>0: apply static weights after n-adaptations
  m_dlbStaticWeightMode = -1;
  if(Context::propertyExists("dlbStaticWeights", 0)) {
    m_dlbStaticWeightMode = 0;
    if(Context::propertyLength("dlbStaticWeights", 0) != globalNoLoadTypes()) {
      cerr0 << "Dlb static Load size does not match!" << globalNoLoadTypes() << " "
            << Context::propertyLength("dlbStaticWeights", 0) << std::endl;
    } else {
      mAlloc(m_dlbStaticWeights, globalNoLoadTypes(), "m_dlbStaticWeight", AT_);
      for(MInt i = 0; i < globalNoLoadTypes(); i++) {
        m_dlbStaticWeights[i] = Context::getBasicProperty<MFloat>("dlbStaticWeights", AT_, i);
      }
    }
    m_dlbStaticWeightMode = Context::getBasicProperty<MInt>("dlbStaticWeightMode", AT_, &m_dlbStaticWeightMode);
  }
}


/// Initialize timers
template <MInt nDim>
void Controller<nDim>::initTimers() {
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "GridController (noSolvers = " + std::to_string(noSolvers()) + ")");
  m_timers.fill(-1);
  NEW_TIMER_NOCREATE(m_timers[Timers::Controller], "total object lifetime", m_timerGroup);
  RECORD_TIMER_START(m_timers[Timers::Controller]);

  // Balancing timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::DLB], "DLB", m_timers[Timers::Controller]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Trigger], "Trigger", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Partition], "Partition", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Prepare], "Prepare", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BalanceGrid], "Balance grid", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BalanceSolvers], "Balance solvers", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BalanceCouplers], "Balance couplers", m_timers[Timers::DLB]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::IO], "I/O", m_timers[Timers::DLB]);

  // Adaptation timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Adaptation], "Adaptation", m_timers[Timers::Controller]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MeshAdaptation], "Mesh Adaptation", m_timers[Timers::Adaptation]);

  const MInt noSolversAndCouplers = noSolvers() + noCouplers();
  m_solverTimerGroups.resize(noSolversAndCouplers);
  m_solverTimers.resize(noSolversAndCouplers);

  for(MInt b = 0; b < noSolversAndCouplers; b++) {
    m_solverTimerGroups[b] = -1;
    m_solverTimers[b].fill(-1);

    const MBool isSolver = (b < noSolvers());
    const MString groupName = (isSolver) ? "solverId = " + std::to_string(solver(b).solverId())
                                         : "couplerId = " + std::to_string(coupler(b - noSolvers()).couplerId());

    NEW_TIMER_GROUP_NOCREATE(m_solverTimerGroups[b], "GridController: " + groupName);
    NEW_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::Total], "total", m_solverTimerGroups[b]);
    const MInt solverTimer = m_solverTimers[b][SolverTimers::Total];

    // Balancing timers
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::DLB], "Dynamic load balancing", solverTimer);
    const MInt dlbTimer = m_solverTimers[b][SolverTimers::DLB];

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::Reset], "Reset solver", dlbTimer);
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::CalcDataSizes], "Determine data sizes", dlbTimer);
    const MInt dataSizeTimer = m_solverTimers[b][SolverTimers::CalcDataSizes];

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::CalcDataSizesMpiBlocking], "initial blocking MPI",
                           dataSizeTimer);
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::CalcDataSizesMpi], "remaining MPI", dataSizeTimer);

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::Balance], "Balance", dlbTimer);
    const MInt balanceTimer = m_solverTimers[b][SolverTimers::Balance];

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::Redistribute], "Redistribute data", balanceTimer);
    const MInt redistTimer = m_solverTimers[b][SolverTimers::Redistribute];

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::RedistributeMpiBlocking], "initial blocking MPI",
                           redistTimer);
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::RedistributeMpi], "remaining MPI", redistTimer);

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::BalancePre], "Balance pre", balanceTimer);

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::SetData], "Set data", balanceTimer);
    const MInt setDataTimer = m_solverTimers[b][SolverTimers::SetData];

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::SetDataMpiBlocking], "initial blocking MPI", setDataTimer);
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::SetDataMpi], "remaining MPI", setDataTimer);

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::BalancePost], "Balance post", balanceTimer);

    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::DlbOther], "Other", dlbTimer);

    // Adaptation timers
    NEW_SUB_TIMER_NOCREATE(m_solverTimers[b][SolverTimers::Adaptation], "Adaptation", solverTimer);
  }
}


/// Check if load balancing is necessary and if yes, call appropriate grid & solver methods.
///
/// Return true if actual load balancing took place.
template <MInt nDim>
MBool Controller<nDim>::balance(const MBool force, const MBool finalTimeStep, const MBool adaptation) {
  // Return if not enabled or in serial
  if(!m_balance || noDomains() == 1) {
    return false;
  }

  // Reset balance status
  gridb().m_wasBalanced = false;

  std::vector<std::pair<MFloat, MString>> durations{};
  auto logDuration = [&durations](const MFloat time, const MString comment, const MBool fromTimer = false) {
    const MFloat duration = (fromTimer) ? time : wallTime() - time;
    durations.push_back(std::make_pair(duration, comment));
  };

  RECORD_TIMER_START(m_timers[Timers::DLB]);

  const MFloat dlbStartTime = wallTime();
  const MInt oldAllocatedBytes = allocatedBytes();

  const MInt noDlbTimers = maia::dlb::g_dlbTimerController.noDlbTimers();
  if(noDlbTimers == 0 && m_loadBalancingMode == 1 && !m_testDynamicLoadBalancing) {
    std::cerr << "There are no DLB timers, but loadBalancingMode is 1 and testing is off; "
                 "switching to loadBalancingMode = 0."
              << std::endl;
    m_loadBalancingMode = 0;
  }

  // Accumulate timer records of all dlb timers
  MFloat localRunTime = 0.0;
  MFloat localIdleTime = 0.0;

  for(MInt i = 0; i < noDlbTimers; i++) {
    const MFloat loadRecord = maia::dlb::g_dlbTimerController.returnLoadRecord(i, m_loadBalancingTimerMode);
    const MFloat idleRecord = maia::dlb::g_dlbTimerController.returnIdleRecord(i, m_loadBalancingTimerMode);

    if(m_loadBalancingMode == 1 && !m_testDynamicLoadBalancing && (loadRecord < 0.0 || idleRecord < 0.0)) {
      TERMM(1, "Load/Idle record for dlb timer #" + std::to_string(i) + " is less than zero on global domain #"
                   + std::to_string(domainId()) + ": " + std::to_string(loadRecord) + ", "
                   + std::to_string(idleRecord));
    }

    localRunTime += loadRecord;
    localIdleTime += idleRecord;
  }

  const MInt timeStep = globalTimeStep;
  if(!adaptation) { // Skip storing timings for call from adaptation
    // Runtime and idle-time of last step
    const MFloat localRunTimeLastStep = localRunTime - m_dlbPreviousLocalRunTime;
    const MFloat localIdleTimeLastStep = localIdleTime - m_dlbPreviousLocalIdleTime;

    // Store current timings for next evaluation
    m_dlbPreviousLocalRunTime = localRunTime;
    m_dlbPreviousLocalIdleTime = localIdleTime;

    // Store timing for imbalance evaluation
    m_dlbTimings.push_back(localRunTimeLastStep);

    // Store runtime and idle time for performance evaluations
    if(m_outputDlbTimings) {
      m_dlbTimeStepsAll.push_back(timeStep);
      m_dlbRunTimeAll.push_back(localRunTimeLastStep);
      m_dlbIdleTimeAll.push_back(localIdleTimeLastStep);

      // Write timings to disk at final time step
      if(finalTimeStep) {
        storeTimings();
      }
    }
  }

  // Check if this is a timer reset time step
  MBool resetTimeStep = (m_lastLoadBalancingTimeStep + m_loadBalancingTimerStartOffset == timeStep
                         || timeStep == m_loadBalancingStartTimeStep + m_loadBalancingTimerStartOffset);

  // TODO labels:DLB @ansgar FIXME
  if(m_balanceAfterAdaptation) {
    // resetTimeStep = ((timeStep - m_loadBalancingTimerStartOffset) % m_loadBalancingInterval == 0
    resetTimeStep = resetTimeStep || (timeStep - m_lastAdaptationTimeStep == m_loadBalancingTimerStartOffset);
    // Reset timer at adaptation time step (if this is no DLB step, see below)
    resetTimeStep = resetTimeStep || (timeStep - m_lastAdaptationTimeStep == 0);
  }

  MBool dlbTimeStep = ((timeStep - m_loadBalancingOffset) % m_loadBalancingInterval == 0);
  if(m_balanceAfterAdaptation) {
    dlbTimeStep = dlbTimeStep
                  || (((timeStep - m_lastAdaptationTimeStep) == m_loadBalancingOffset)
                      && (m_nAdaptationsSinceBalance >= m_balanceAdaptationInterval || m_dlbStep < 2)
                      && (m_loadBalancingStopTimeStep < 0 || globalTimeStep < m_loadBalancingStopTimeStep));
  }

  if(resetTimeStep && !dlbTimeStep) {
    // reset all timer
    resetAllTimer();

    // ... and return
    RECORD_TIMER_STOP(m_timers[Timers::DLB]);
    return false;
  }

  // Quit early if balancing is not forced and interval does not match
  if(!force && !finalTimeStep && (m_loadBalancingInterval <= 0 || !dlbTimeStep)) {
    RECORD_TIMER_STOP(m_timers[Timers::DLB]);
    return false;
  }

  // const MBool atDlbInterval = (timeStep % m_loadBalancingInterval == 0);
  const MBool initialBalance = (timeStep == -1 && force);
  const MBool initialAdaptationOnly = (timeStep == -1 && m_performanceOutput);
  const MBool afterLastDlbStep = (timeStep > m_loadBalancingStopTimeStep && m_loadBalancingStopTimeStep != -1);
  const MBool beforeFirstDlbStep = (timeStep <= m_loadBalancingStartTimeStep);
  // Performance output step without load balancing
  const MBool performanceOutput =
      (m_performanceOutput || finalTimeStep || afterLastDlbStep || beforeFirstDlbStep) && !initialBalance;

  // Return if this is not a DLB time step and not the initial balancing
  if(!dlbTimeStep && (!initialBalance || initialAdaptationOnly)) {
    // Continue if forced or at final time step to determine imbalance before returning
    if((!force && !finalTimeStep) || initialAdaptationOnly) {
      RECORD_TIMER_STOP(m_timers[Timers::DLB]);
      return false;
    }
  }

  if(!performanceOutput) {
    m_log << "Dynamic load balancing at time step " << timeStep << std::endl;
    cerr0 << "=== Dynamic load balancing at time step " << timeStep << std::endl;
  }

  logTimerStatistics("before balance");

  // Determine if this is the last DLB step
  const MBool lastDlbStep =
      (timeStep + m_loadBalancingInterval > m_loadBalancingStopTimeStep && m_loadBalancingStopTimeStep != -1)
      && !m_balanceAfterAdaptation;
  // Determine if this is a DLB step at which the partitioning should be reverted to the best so far
  // to improve it by using local shifts only
  const MInt isDlbRevertStep =
      (m_loadBalancingMode == 1 && !lastDlbStep && m_loadBalancingStopTimeStep != -1 && m_dlbNoFinalLocalShifts > 0
       && timeStep == m_loadBalancingStopTimeStep - m_loadBalancingInterval * (m_dlbNoFinalLocalShifts + 1));


  // Determine the number of timesteps that are timed
  // MInt noTimeSteps = timeStep - m_lastLoadBalancingTimeStep - m_loadBalancingTimerStartOffset;
  MInt noTimeSteps = timeStep - m_dlbLastResetTimeStep;
  if(timeStep < 0) noTimeSteps = 0;

  if(m_loadBalancingMode == 1) {
    ASSERT(noTimeSteps > 0, "ERROR: noTimeSteps = " + std::to_string(noTimeSteps));
  }

  if(noTimeSteps != (MInt)m_dlbTimings.size()) {
    m_log << "DLB: Number of timings does not match: " + std::to_string(m_dlbTimings.size()) + " "
                 + std::to_string(noTimeSteps) + " (beforeFirstStep = " + std::to_string(beforeFirstDlbStep)
                 + "; afterLastDlbStep = " + std::to_string(afterLastDlbStep) + ")"
          << std::endl;
    std::cerr << "WARNING: number of timings mismatch! " << noTimeSteps << " " << m_dlbTimings.size() << std::endl;
    noTimeSteps = m_dlbTimings.size();
  }


  std::vector<MFloat> timings(m_dlbTimings.begin(), m_dlbTimings.end());
  // Sort timings for evaluation
  std::sort(timings.begin(), timings.end());
  const MInt noSamples = timings.size();
  if(noSamples < 1 && !initialBalance) {
    TERMM(1, "ERROR: number of samples is < 1");
  }

  const MBool truncatedMean = !m_balanceAfterAdaptation;
  // Compute 25% truncated mean for computing loads for static grids
  const MFloat trim = 0.25;
  const MInt lowerBound = std::floor(trim * noSamples);
  const MInt upperBound = std::max(MInt(std::floor((1 - trim) * noSamples)), std::min(noSamples, 1));
  const MInt noSamplesTruncated = upperBound - lowerBound;

  if(m_loadBalancingMode == 1 && noSamplesTruncated < 1 && !initialBalance) {
    std::cerr << "ERROR no samples truncated " << noSamplesTruncated << std::endl;
    /* return 0; */
  }

  const MFloat truncatedMeanRunTime =
      std::accumulate(&timings[lowerBound], &timings[upperBound], 0.0) / noSamplesTruncated;
  const MFloat meanRunTime = std::accumulate(timings.begin(), timings.end(), 0.0) / noSamples;
  const MFloat meanRunTimeTrigger = (truncatedMean) ? truncatedMeanRunTime : meanRunTime;

  /* const MFloat maxLocalRunTime = *std::max_element(timings.begin(), timings.end()); */
  /* const MFloat minLocalRunTime = *std::min_element(timings.begin(), timings.end()); */

  // Compute average time per step (all timings)
  MFloat timePerStep = (localRunTime + localIdleTime) / noTimeSteps;
  MFloat maxRunTime = localRunTime / noTimeSteps;
  const MFloat localTimePerStep = timePerStep;
  // Communicate and take maximum value -> assure same value on all domains
  MPI_Allreduce(MPI_IN_PLACE, &timePerStep, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "timePerStep");
  MPI_Allreduce(MPI_IN_PLACE, &maxRunTime, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxRunTime");

  MFloat maxIdleTime = localIdleTime / noTimeSteps;
  MFloat minIdleTime = localIdleTime / noTimeSteps;
  MPI_Allreduce(MPI_IN_PLACE, &maxIdleTime, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxIdleTime");
  MPI_Allreduce(MPI_IN_PLACE, &minIdleTime, 1, MPI_DOUBLE, MPI_MIN, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "minIdleTime");

  const MFloat minIdleTimeRel = minIdleTime / timePerStep;

  m_log << timeStep << " * Average time per step: " << timePerStep << " ; local: " << localTimePerStep
        << ", idle/comp = " << localIdleTime / localRunTime
        << ", idle/timePerStep = " << localIdleTime / (noTimeSteps * timePerStep) << std::endl;
  m_log << timeStep << " * Relative idle time: max = " << maxIdleTime / timePerStep << ", min "
        << minIdleTime / timePerStep << std::endl;
  m_log << timeStep << " * maxRunTime " << maxRunTime << std::endl;

  logDuration(dlbStartTime, "DLB preparation/timers");

  const MFloat imbalanceStartTime = wallTime();
  RECORD_TIMER_START(m_timers[Timers::Trigger]);
  MFloat imbalance = -1.0;
  ScratchSpace<MFloat> loads(noDomains(), FUN_, "loads");
  MBool loadBalance = initialBalance;
  if(!loadBalance) { // Skip this for initial balance since there are no timings
    loadBalance = needLoadBalancing(meanRunTimeTrigger, &loads[0], imbalance) || force || (m_forceBalance > 0);
  }
  if(m_forceBalance == 1) {
    m_forceBalance = 0;
  }
  RECORD_TIMER_STOP(m_timers[Timers::Trigger]);
  logDuration(imbalanceStartTime, "Imbalance evaluation");

  // Determine variance of timings to detect performance variations
  // TODO labels:CONTROLLER,DLB exclude ranks with low load?
  MFloat localRunTimeVariance = 0.0;
  for(MInt i = 0; i < noSamples; i++) {
    localRunTimeVariance += POW2(timings[i] - meanRunTime);
  }
  localRunTimeVariance /= noSamples;
  const MFloat localRunTimeStdev = std::sqrt(localRunTimeVariance);

  const MFloat performanceVariation = localRunTimeStdev / meanRunTime;
  MFloat perfVarMax[2];
  perfVarMax[0] = performanceVariation;
  // Note: use only ranks above a specific load
  perfVarMax[1] = (loads[domainId()] > 0.85) ? performanceVariation : 0.0;

  MPI_Allreduce(MPI_IN_PLACE, &perfVarMax, 2, type_traits<MFloat>::mpiType(), MPI_MAX, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "perfVarMax");
  const MFloat maxPerformanceVariation = perfVarMax[1];

  const MBool performanceVariationCheck = (maxPerformanceVariation > m_maxPerformanceVarThreshold && !m_adaptation);

  m_log << timeStep << " * DLB: Performance variation: load>0.85 " << maxPerformanceVariation << "; all "
        << perfVarMax[0] << std::endl;

  m_log << timeStep << " * Average timings: avgRunTime = " << localRunTime / noTimeSteps
        << ", truncated = " << truncatedMean << " " << truncatedMeanRunTime << ", diff "
        << meanRunTime - truncatedMeanRunTime << ", noTimeSteps = " << noTimeSteps << std::endl;

  m_log << timeStep << " * timePerStep = " << timePerStep << "; imbalance = " << imbalance
        << "; maxRunTime = " << maxRunTime << "; minIdleTimeRel = " << minIdleTimeRel << std::endl;

  if(m_loadBalancingMode == 1 && performanceVariationCheck && !lastDlbStep && !performanceOutput
     && !m_testDynamicLoadBalancing) {
    /* std::cerr << domainId() << " performance variation " << minLocalRunTime << " " */
    /*           << maxLocalRunTime << " " << localRunTimeVariance << " " << localRunTimeStdev */
    /*           << " " << localRunTimeStdev/meanRunTime << std::endl; */

    std::stringstream perfMessage;
    perfMessage << "DLB: Performance variation " << maxPerformanceVariation << ", skip DLB step." << std::endl;
    m_log << perfMessage.str();
    if(domainId() == 0) {
      std::cout << perfMessage.str();
    }

    resetAllTimer();

    m_lastLoadBalancingTimeStep = timeStep;
    m_nAdaptationsSinceBalance = 0;
    RECORD_TIMER_STOP(m_timers[Timers::DLB]);
    return false;
  }

  if(m_loadBalancingMode == 0) {
    m_lastLoadBalancingTimeStep = timeStep;
    m_nAdaptationsSinceBalance = 0;
  }

  // Check if the current time per timestep is the best one so far. If so, store the new best
  // time/timestep and the corresponding local domain offset. (Note: includes idle times, also the
  // imbalance should not be too large, accept also better imbalance with slightly increased time
  // per step as the latter might be due e.g. network latencies and slower communication)
  MBool newBestTimePerStep = false;
  if(m_loadBalancingMode == 1 && !performanceOutput && !performanceVariationCheck && !m_testDynamicLoadBalancing
     && ((timePerStep < m_timePerStepTotal && imbalance < 1.05 * m_imbalance)
         || (imbalance < m_imbalance && timePerStep < 1.05 * m_timePerStepTotal))) {
    if(m_optPartitionCellOffsetTotal == -1) {
      m_timePerStepTotal = -1.0;
      m_imbalance = -1.0;
    }
    m_log << timeStep << " * Storing new best domain partitioning: timePerStep = " << timePerStep << " ("
          << m_timePerStepTotal << "); imbalance = " << imbalance << " (" << m_imbalance << "), maxRunTime "
          << maxRunTime << " minIdleTimeRel " << minIdleTimeRel << std::endl;
    newBestTimePerStep = true;
    m_timePerStepTotal = timePerStep;
    m_imbalance = imbalance;
    m_optPartitionCellOffsetTotal = gridb().m_localPartitionCellOffsets[0];

    RECORD_TIMER_START(m_timers[Timers::IO]);
    // Store partition file for restarting with the same domain offsets on the
    // same number of domains
    // This will overwrite the partition_n[noDomains].[ext] file every time a
    // better partitioning is found.
    std::stringstream partitionFileName;
    partitionFileName << "partition_n" << noDomains();
    gridb().savePartitionFile(partitionFileName.str(), m_optPartitionCellOffsetTotal);
    RECORD_TIMER_STOP(m_timers[Timers::IO]);
  } else if(m_testDynamicLoadBalancing && finalTimeStep) {
    // Always write partition file at final time step for testing
    std::stringstream partitionFileName;
    partitionFileName << "partition_n" << noDomains();
    gridb().savePartitionFile(partitionFileName.str(), gridb().m_localPartitionCellOffsets[0]);
  }

  // Return if nothing needs to be done
  if(!loadBalance && !m_forceLoadBalancing && (isDlbRevertStep == 0)) {
    m_log << " * no load imbalance detected at timestep " << timeStep << "!" << std::endl;
    // Continue if this is the last DLB step and the time per step is worse than
    // the best one, else return here
    if(!(lastDlbStep && !newBestTimePerStep)) {
      RECORD_TIMER_STOP(m_timers[Timers::DLB]);
      return false;
    }
  } else {
    if(!performanceOutput) {
      m_log << " * load imbalance detected at timestep " << timeStep << "!" << std::endl;
    }
  }

  // Return here for final time step/or performance output after determining imbalance etc.
  if(performanceOutput) {
    // clear timings such that only the next x timesteps are used for the following performance
    // output
    m_log << "Performance evaluation: clear timings at timestep " << timeStep << "!" << std::endl;

    resetAllTimer();

    RECORD_TIMER_STOP(m_timers[Timers::DLB]);
    return false;
  }

  // Determine imbalance statistics *BEFORE* online restart
  printDomainStatistics("before load balancing");

  // reset newMinLeven for balancing, otherwise cells with zero noOffsprings occour!
  MInt backupLevel = -1;
  if(gridb().m_newMinLevel > 0) {
    backupLevel = gridb().m_newMinLevel;
    gridb().m_newMinLevel = -1;
  }

  gridb().computeGlobalIds();
  gridb().storeMinLevelCells();

  // Cancel any open MPI (receive) requests already opened to allow for interleaved execution, since
  // these might lead to conflicting messages
  for(MInt i = 0; i < noSolvers(); i++) {
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DlbOther]);
    solver(i).cancelMpiRequests();
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DlbOther]);

    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
  }

  MBool partitionCellChange = false;
  const MFloat updatePartCellsStartTime = wallTime();
  // If enabled: Check if the partition cells need to be updated/changed, i.e., decrease/increase
  // the partition level shift
  //
  // Note: in case of a change in the partition cells the balancing needs to be performed (even if
  // the determined new partitioning is the same as after the partition cell update)
  // TODO labels:CONTROLLER,DLB last DLB or revert step!
  if(m_dlbUpdatePartitionCells && !lastDlbStep) {
    MInt offspringThreshold = gridb().m_partitionCellOffspringThreshold;
    MFloat workloadThreshold = gridb().m_partitionCellWorkloadThreshold;

    // For testing: set different thresholds to force partition cell changes
    if(m_testUpdatePartitionCells && m_dlbStep % 2 == 0) {
      offspringThreshold = m_testUpdatePartCellsOffspringThreshold;
      workloadThreshold = m_testUpdatePartCellsWorkloadThreshold;
    }

    // Update partition cells
    partitionCellChange = gridb().updatePartitionCells(offspringThreshold, workloadThreshold);

    if(partitionCellChange) {
      // Grid with new partition cells will be written with next restart files
      m_saveGridNewPartitionCells = true;

      // Reset best found partitioning since not useful anymore with changed partition cells
      m_optPartitionCellOffsetTotal = -1;
    }
  }
  logDuration(updatePartCellsStartTime, "Update partition cells");

  // Storage for new partitioning
  MLongScratchSpace partitionCellOffsets(noDomains() + 1, AT_, "partitionCellOffsets");
  MLongScratchSpace globalIdOffsets(noDomains() + 1, AT_, "globalIdOffsets");

  const MFloat partitionStartTime = wallTime();
  // Determine new partitioning based on load balancing mode
  if(m_loadBalancingMode == 0 || timeStep == -1) { // Use mode 0 always for initial balance
    RECORD_TIMER_START(m_timers[Timers::Partition]);
    accumulateCellWeights();

    // Compute new domain decomposition
    partition(&partitionCellOffsets[0], &globalIdOffsets[0], false);

    // Increase DLB step
    m_dlbStep++;

    RECORD_TIMER_STOP(m_timers[Timers::Partition]);
  } else {
    // Static mode (i.e. fixed number of DLB steps)
    RECORD_TIMER_START(m_timers[Timers::Partition]);

    MBool newPartition = false;
    // Check if this is the last DLB step
    if((!lastDlbStep || m_testDynamicLoadBalancing) && ((isDlbRevertStep == 0) || newBestTimePerStep)) {
      // Determine new partitioning, i.e. new partition cell offsets and domain offsets
      newPartition = loadBalancingPartition(&loads[0], imbalance, &partitionCellOffsets[0], &globalIdOffsets[0]);
    } else {
      // Check if the current configuration is not the best one (wrt. time per timestep)
      if(((!newBestTimePerStep && !m_adaptation) || (isDlbRevertStep != 0)) && m_optPartitionCellOffsetTotal != -1) {
        m_log << timeStep << " * DLB reverting to best configuration." << std::endl;
        // Last DLB step: since there is still some load imbalance, revert to the best partitioning
        // found

        // Gather current partition cell offsets
        MLongScratchSpace localPartitionCellOffsets(noDomains() + 1, AT_, "localPartitionCellOffsets");
        MPI_Allgather(&gridb().m_localPartitionCellOffsets[0], 1, type_traits<MLong>::mpiType(),
                      &localPartitionCellOffsets[0], 1, type_traits<MLong>::mpiType(), gridb().mpiComm(), AT_,
                      "gridb().m_localPartitionCellOffsets[0]", "localPartitionCellOffsets[0]");
        localPartitionCellOffsets[noDomains()] = gridb().m_noPartitionCellsGlobal;

        // Gather 'optimal' partition cell offsets
        MPI_Allgather(&m_optPartitionCellOffsetTotal, 1, type_traits<MLong>::mpiType(), &partitionCellOffsets[0], 1,
                      type_traits<MLong>::mpiType(), gridb().mpiComm(), AT_, "m_optPartitionCellOffsetTotal",
                      "partitionCellOffsets[0]");
        partitionCellOffsets[noDomains()] = gridb().m_noPartitionCellsGlobal;

        // Check if this is a new partitioning, i.e. different from the old one
        newPartition =
            (std::mismatch(partitionCellOffsets.begin(), partitionCellOffsets.end(), &localPartitionCellOffsets[0]))
                .first
            != partitionCellOffsets.end();

        loadBalancingCalcNewGlobalOffsets(&localPartitionCellOffsets[0], &partitionCellOffsets[0], &globalIdOffsets[0]);
      } else {
        // Current partitioning is accepted
        // Note: this will not work if m_forceLoadBalancing is activated, since
        // the number of cells to send, etc. are not determined.
        newPartition = false;
      }
    }

    RECORD_TIMER_STOP(m_timers[Timers::Partition]);

    // Return here if partitioning did not change
    if(!newPartition && !partitionCellChange && !m_forceLoadBalancing) {
      m_log << "Dynamic load balancing: load imbalance detected but partition did not change at "
               "timestep "
            << timeStep << "!" << std::endl;
      RECORD_TIMER_STOP(m_timers[Timers::DLB]);
      return false;
    }
  }
  logDuration(partitionStartTime, "New partitioning");

  RECORD_TIMER_START(m_timers[Timers::Prepare]);
  const MFloat prepareStartTime = wallTime();

  // Disable all running DLB timers and store status
  ScratchSpace<MBool> dlbTimersStatus(std::max(noDlbTimers, 1), AT_, "dlbTimersStatus");
  maia::dlb::g_dlbTimerController.disableAllDlbTimers(&dlbTimersStatus[0]);

  // Reset solvers before load balancing
  for(MInt i = 0; i < noSolvers(); i++) {
    const MFloat resetStartTime = wallTime();
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Reset]);
    solver(i).resetSolver();
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Reset]);


    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DlbOther]);
    // Communicate global solver variables from current local rank 0 to all domains
    communicateGlobalSolverVars(&solver(i));
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DlbOther]);

    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
    logDuration(resetStartTime, "  Reset solver #" + std::to_string(i));
  }

  const MFloat prepareGridStartTime = wallTime();
  gridb().deletePeriodicConnection(false);

  std::vector<MInt>().swap(gridb().m_minLevelCells);

  gridb().localToGlobalIds();

  const MInt oldNoCells = gridb().treeb().size();

  MIntScratchSpace noCellsToReceiveByDomain(noDomains() + 1, AT_, "noCellsToReceiveByDomain");
  MIntScratchSpace noCellsToSendByDomain(noDomains() + 1, AT_, "noCellsToSendByDomain");
  // Mapping: gridCellId -> sort index in buffer
  MIntScratchSpace sortedCellId(oldNoCells, AT_, "sortedCellId");
  // Reverse mapping: sort index in buffer -> gridCellId
  MIntScratchSpace bufferIdToCellId(oldNoCells, AT_, "bufferIdToCellId");
  bufferIdToCellId.fill(-1);

  loadBalancingCalcNoCellsToSend(&globalIdOffsets[0], &noCellsToSendByDomain[0], &noCellsToReceiveByDomain[0],
                                 &sortedCellId[0], &bufferIdToCellId[0]);
  logDuration(prepareGridStartTime, "  Prepare grid");

  std::vector<std::vector<MInt>> dataSendSize{};
  std::vector<std::vector<MInt>> dataRecvSize{};

  // Determine communication data sizes for each solver
  for(MInt i = 0; i < noSolvers(); i++) {
    const MFloat dataStartTime = wallTime();
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::CalcDataSizes]);
    determineDataSizesDlb(i, 0, &noCellsToSendByDomain[0], &bufferIdToCellId[0], dataSendSize, dataRecvSize);
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::CalcDataSizes]);

    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DlbOther]);
    // Change local to global ids here, since local ids previously needed to determine data sizes
    solver(i).localToGlobalIds();
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DlbOther]);

    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
    logDuration(dataStartTime, "  Data sizes solver #" + std::to_string(i));
  }


  std::vector<std::vector<MInt>> dataSendSizeCoupler{};
  std::vector<std::vector<MInt>> dataRecvSizeCoupler{};

  // Determine communication data sizes for each coupler
  for(MInt i = 0; i < noCouplers(); i++) {
    const MFloat dataStartTime = wallTime();
    const MInt timerId = noSolvers() + i;
    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::CalcDataSizes]);
    determineDataSizesDlb(i, 1, &noCellsToSendByDomain[0], &bufferIdToCellId[0], dataSendSizeCoupler,
                          dataRecvSizeCoupler);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::CalcDataSizes]);

    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Total]);
    logDuration(dataStartTime, "  Data sizes coupler #" + std::to_string(i));
  }


  RECORD_TIMER_STOP(m_timers[Timers::Prepare]);
  logDuration(prepareStartTime, "Prepare balancing");

  const MFloat balanceGridStartTime = wallTime();
  RECORD_TIMER_START(m_timers[Timers::BalanceGrid]);
  // Rebalance grid
  gridb().balance(&noCellsToReceiveByDomain[0], &noCellsToSendByDomain[0], &sortedCellId[0], &partitionCellOffsets[0],
                  &globalIdOffsets[0]);
  RECORD_TIMER_STOP(m_timers[Timers::BalanceGrid]);
  logDuration(balanceGridStartTime, "Balance grid");

  RECORD_TIMER_START(m_timers[Timers::BalanceSolvers]);
  const MFloat balanceSolversStartTime = wallTime();
  // Rebalance solvers
  for(MInt i = 0; i < noSolvers(); i++) {
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Balance]);
    const MFloat balanceSolverStartTime = wallTime();

    if(!solver(i).hasSplitBalancing()) {
      solver(i).balance(&noCellsToReceiveByDomain[0], &noCellsToSendByDomain[0], &sortedCellId[0], oldNoCells);
    } else {
      // Variables to store pointers to allocated memory
      std::vector<MInt*> intDataRecv{};
      std::vector<MLong*> longDataRecv{};
      std::vector<MFloat*> floatDataRecv{};
      std::vector<MInt> dataTypes{};

      // TODO labels:CONTROLLER,DLB,TIMERS create shorthand for: timer start, function call, timer stop, logDuration?
      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Redistribute]);
      // Communicate all solver data
      redistributeDataDlb(i, 0, dataSendSize[i], dataRecvSize[i], &bufferIdToCellId[0], oldNoCells, intDataRecv,
                          longDataRecv, floatDataRecv, dataTypes);
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Redistribute]);
      logDuration(RETURN_TIMER(m_solverTimers[i][SolverTimers::Redistribute]),
                  "Redistribute solver #" + std::to_string(i), true);

      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::BalancePre]);
      // Perform reinitialization steps prior to setting solver data
      // i.e. apply balance to the proxy!
      solver(i).balancePre();
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::BalancePre]);
      logDuration(RETURN_TIMER(m_solverTimers[i][SolverTimers::BalancePre]), "BalancePre solver #" + std::to_string(i),
                  true);

      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::SetData]);
      setDataDlb(i, 0, intDataRecv, longDataRecv, floatDataRecv, dataTypes, false);
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::SetData]);
      logDuration(RETURN_TIMER(m_solverTimers[i][SolverTimers::SetData]), "SetData1 solver #" + std::to_string(i),
                  true);

      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::DlbOther]);
      // Change global to local ids in the solver
      solver(i).globalToLocalIds();
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DlbOther]);

      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::BalancePost]);
      // Perform reinitialization steps after setting solver data the first time
      solver(i).balancePost();
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::BalancePost]);
      logDuration(RETURN_TIMER(m_solverTimers[i][SolverTimers::BalancePost]),
                  "BalancePost solver #" + std::to_string(i), true);

      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::SetData]);
      // Set solver data again (if required by solver) and deallocate buffers
      setDataDlb(i, 0, intDataRecv, longDataRecv, floatDataRecv, dataTypes, true);
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::SetData]);
      logDuration(RETURN_TIMER(m_solverTimers[i][SolverTimers::SetData]), "SetData2 solver #" + std::to_string(i),
                  true);

      intDataRecv.clear();
      longDataRecv.clear();
      floatDataRecv.clear();
    }

    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Balance]);
    logDuration(balanceSolverStartTime, "Balance solver #" + std::to_string(i));

    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
  }
  RECORD_TIMER_STOP(m_timers[Timers::BalanceSolvers]);
  logDuration(balanceSolversStartTime, "Balance solvers");

  RECORD_TIMER_START(m_timers[Timers::BalanceCouplers]);
  const MFloat balanceCouplersStartTime = wallTime();
  // Rebalance couplers (if required)
  for(MInt i = 0; i < noCouplers(); i++) {
    const MInt timerId = noSolvers() + i;
    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Total]);
    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::DLB]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Balance]);

    // Variables to store pointers to allocated memory
    std::vector<MInt*> intDataRecv{};
    std::vector<MLong*> longDataRecv{};
    std::vector<MFloat*> floatDataRecv{};
    std::vector<MInt> dataTypes{};

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Redistribute]);
    // Communicate all coupler data
    redistributeDataDlb(i, 1, dataSendSizeCoupler[i], dataRecvSizeCoupler[i], &bufferIdToCellId[0], oldNoCells,
                        intDataRecv, longDataRecv, floatDataRecv, dataTypes);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Redistribute]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::BalancePre]);
    // Perform reinitialization steps prior to setting solver data
    coupler(i).balancePre();
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::BalancePre]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::SetData]);
    setDataDlb(i, 1, intDataRecv, longDataRecv, floatDataRecv, dataTypes, false);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::SetData]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::BalancePost]);
    // Perform reinitialization steps after setting solver data the first time
    coupler(i).balancePost();
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::BalancePost]);

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::SetData]);
    // Set coupler data again (if required by coupler) and deallocate buffers
    setDataDlb(i, 1, intDataRecv, longDataRecv, floatDataRecv, dataTypes, true);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::SetData]);

    intDataRecv.clear();
    longDataRecv.clear();
    floatDataRecv.clear();

    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Balance]);

    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::DLB]);
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Total]);
  }
  RECORD_TIMER_STOP(m_timers[Timers::BalanceCouplers]);
  logDuration(balanceCouplersStartTime, "Balance couplers");

  const MFloat finalizeBalanceStartTime = wallTime();
  // finalize balance for solvers and couplers
  for(MInt i = 0; i < noSolvers(); i++) {
    const MFloat finalizeSolverStartTime = wallTime();
    solver(i).finalizeBalance();
    for(MInt j = 0; j < noCouplers(); j++) {
      const MFloat finalizeCouplerStartTime = wallTime();
      coupler(j).finalizeBalance(i);
      logDuration(finalizeCouplerStartTime,
                  "Finalize balance #" + std::to_string(i) + " coupler #" + std::to_string(j));
    }
    logDuration(finalizeSolverStartTime, "Finalize balance solver #" + std::to_string(i));
  }
  logDuration(finalizeBalanceStartTime, "Finalize balance total");

  // Determine imbalance statistics *AFTER* online restart
  printDomainStatistics("after load balancing");

  // Enable all previously active DLB timers again
  maia::dlb::g_dlbTimerController.enableAllDlbTimers(&dlbTimersStatus[0]);

  // Reset DLB timer records
  maia::dlb::g_dlbTimerController.resetRecords();

  m_dlbPreviousLocalRunTime = 0.0;
  m_dlbPreviousLocalIdleTime = 0.0;

  // ... reset timings ...
  resetAllTimer();

  m_lastLoadBalancingTimeStep = timeStep;
  m_nAdaptationsSinceBalance = 0;

  // Write timings to file
  if(m_outputDlbTimings) {
    storeTimings();
  }

  logDuration(dlbStartTime, "Balance total");
  logDurations(durations, "DLB", gridb().mpiComm(), globalDomainId(), globalNoDomains());

  const MFloat dlbTimeTotal = wallTime() - dlbStartTime;
  std::stringstream dlbMessage;
  dlbMessage << "=== Dynamic load balancing performed at timestep " << timeStep << "! Duration: " << dlbTimeTotal
             << " s" << std::endl;
  m_log << dlbMessage.str();
  cerr0 << dlbMessage.str();

  printAllocatedMemory(oldAllocatedBytes, "gridcontroller::balance()", gridb().mpiComm());

  // restore backup value for next restartFile!
  if(backupLevel > 0) {
    gridb().m_newMinLevel = backupLevel;
  }

  RECORD_TIMER_STOP(m_timers[Timers::DLB]);
  return true;
}


/// Compute new domain decomposition
template <MInt nDim>
void Controller<nDim>::partition(MLong* partitionCellOffsets,
                                 MLong* globalIdOffsets,
                                 const MBool onlyPartitionOffsets) {
  TRACE();
  if(gridb().m_maxPartitionLevelShift > 0) {
    cerr0 << "WARNING: Load balancing with partition level shifts might not fully work!" << std::endl;
  }

  // Note: noPartitionCells = 0 allowed after updatePartitionCells()!
  const MLong noPartitionCells = gridb().m_noPartitionCells;

  // Update noOffsprings and workLoad for all cells and calculate accumulated workload
  MFloat totalWorkload = 0.0;
  MFloatScratchSpace partitionCellsWorkload(std::max(noPartitionCells, 1L), AT_, "partitionCellsWorkload");
  MLongScratchSpace partitionCellsGlobalId(std::max(noPartitionCells, 1L), AT_, "partitionCellsGlobalId");

  gridb().calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), gridb().treeb().size());

  for(MUint i = 0; i < noPartitionCells; i++) {
    const MLong globalCellId = gridb().m_localPartitionCellGlobalIds[i]; // sorted by global id
    const MInt cellId = gridb().globalIdToLocalId(globalCellId, true);

    partitionCellsWorkload(i) = gridb().a_workload(cellId);
    partitionCellsGlobalId(i) = gridb().a_globalId(cellId);
    ASSERT(globalCellId == partitionCellsGlobalId(i),
           "global cell id mismatch! " + std::to_string(globalCellId) + " " + std::to_string(partitionCellsGlobalId(i))
               + " " + std::to_string(cellId));

    ASSERT(i == 0 || partitionCellsGlobalId(i) > partitionCellsGlobalId(i - 1),
           "Partition cells not sorted by global id: " + std::to_string(partitionCellsGlobalId(i))
               + " <= " + std::to_string(partitionCellsGlobalId(i - 1)));
    totalWorkload += gridb().a_workload(cellId);
  }

  MPI_Allreduce(MPI_IN_PLACE, &totalWorkload, 1, MPI_DOUBLE, MPI_SUM, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "totalWorkload");

  MBool calcGlobalOffsets = false;

  if(gridb().m_partitionParallelSplit) {
    maia::grid::partitionParallelSplit(&partitionCellsWorkload[0], noPartitionCells,
                                       gridb().m_localPartitionCellOffsets[0], static_cast<MLong>(noDomains()),
                                       static_cast<MLong>(domainId()), gridb().mpiComm(), &partitionCellOffsets[0]);
    partitionCellOffsets[noDomains()] = gridb().m_noPartitionCellsGlobal;

    if(!onlyPartitionOffsets) {
      // Determine the globalIdOffsets later (required for load balancing mode 0)
      calcGlobalOffsets = true;
    }
  } else {
    gridb().partitionParallel(gridb().m_noPartitionCells, gridb().m_localPartitionCellOffsets[0],
                              &partitionCellsWorkload[0], &partitionCellsGlobalId[0], totalWorkload,
                              partitionCellOffsets, globalIdOffsets, true);

    if(!onlyPartitionOffsets && gridb().m_maxPartitionLevelShift > 0) {
      // Redetermine the global offsets below in case of a partition level shift since
      // partitionParallel() does not include the correction of the offsets as done in
      // correctDomainOffsetsAtPartitionLevelShifts()
      calcGlobalOffsets = true;
    }
  }

  // Determine the globalIdOffsets
  if(calcGlobalOffsets) {
    MLongScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
    MLongScratchSpace localPartitionCellOffsets(noDomains() + 1, AT_, "localPartitionCellOffsets");
    // Determine the number of partition cells on all domains and the partition cell offsets
    gridb().determineNoPartitionCellsAndOffsets(&localPartitionCellCounts[0], &localPartitionCellOffsets[0]);

    // Calculate new global domain offsets from the new partition cell offsets
    loadBalancingCalcNewGlobalOffsets(&localPartitionCellOffsets[0], &partitionCellOffsets[0], &globalIdOffsets[0]);
  }

  if(domainId() == 0) std::cerr << std::endl;
}

/// Print statistics regarding the cell distribution among domains.
template <MInt nDim>
void Controller<nDim>::printDomainStatistics(const MString& status) {
  TRACE();

  MIntScratchSpace minNoSolverCells(noSolvers(), FUN_, "minNoSolverCells");
  MIntScratchSpace maxNoSolverCells(noSolvers(), FUN_, "maxNoSolverCells");
  MIntScratchSpace avgNoSolverCells(noSolvers(), FUN_, "avgNoSolverCells");
  MLongScratchSpace globalNoSolverCells(noSolvers(), FUN_, "globalNoSolverCells");

  MLong noGridLeafCells = 0;
  minNoSolverCells.fill(0);
  maxNoSolverCells.fill(0);
  avgNoSolverCells.fill(0);
  globalNoSolverCells.fill(0);

  for(MInt cellId = 0; cellId < gridb().treeb().size(); cellId++) {
    if(!gridb().a_hasProperty(cellId, Cell::IsHalo) && gridb().a_noChildren(cellId) == 0) {
      noGridLeafCells++;
    }
  }

  MLong globalGridNoLeafCells = noGridLeafCells;
  MLong globalGridNoCells = (MLong)gridb().noInternalCells();
  MInt maxGridNoCells = gridb().noInternalCells();
  MInt minGridNoCells = gridb().noInternalCells();

  for(MInt i = 0; i < noSolvers(); i++) {
    minNoSolverCells(i) = solver(i).noInternalCells();
    maxNoSolverCells(i) = solver(i).noInternalCells();
    globalNoSolverCells(i) = solver(i).noInternalCells();
  }

  MPI_Allreduce(MPI_IN_PLACE, &globalGridNoLeafCells, 1, MPI_LONG, MPI_SUM, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "globalGridNoLeafCells");
  MPI_Allreduce(MPI_IN_PLACE, &globalGridNoCells, 1, MPI_LONG, MPI_SUM, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "globalGridNoCells");
  MPI_Allreduce(MPI_IN_PLACE, &maxGridNoCells, 1, MPI_INT, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxGridNoCells");
  MPI_Allreduce(MPI_IN_PLACE, &minGridNoCells, 1, MPI_INT, MPI_MIN, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "minGridNoCells");
  MPI_Allreduce(MPI_IN_PLACE, &maxNoSolverCells[0], noSolvers(), MPI_INT, MPI_MAX, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "maxNoSolverCells[0]");
  MPI_Allreduce(MPI_IN_PLACE, &minNoSolverCells[0], noSolvers(), MPI_INT, MPI_MIN, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "minNoSolverCells[0]");
  MPI_Allreduce(MPI_IN_PLACE, &globalNoSolverCells[0], noSolvers(), MPI_LONG, MPI_SUM, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "globalNoSolverCells[0]");

  for(MInt i = 0; i < noSolvers(); i++) {
    avgNoSolverCells(i) = (MInt)(((MFloat)globalNoSolverCells(i)) / ((MFloat)noDomains()));
  }

  const MInt avgGridNoLeafCells = (MInt)(((MFloat)globalGridNoLeafCells) / ((MFloat)noDomains()));
  const MInt avgGridNoCells = (MInt)(((MFloat)globalGridNoCells) / ((MFloat)noDomains()));
  MInt devNoGridCells = std::abs(gridb().noInternalCells() - avgGridNoCells);
  MInt devNoGridLeafCells = std::abs(noGridLeafCells - avgGridNoLeafCells);


  MPI_Allreduce(MPI_IN_PLACE, &devNoGridCells, 1, MPI_INT, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "devNoGridCells");
  MPI_Allreduce(MPI_IN_PLACE, &devNoGridLeafCells, 1, MPI_INT, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "devNoGridLeafCells");
  if(domainId() == 0) {
    std::cerr << "Domain statistics ";
    if(!status.empty()) std::cerr << status << " ";
    std::cerr << "at global time step " << globalTimeStep << ": avg no cells=" << avgGridNoCells
              << ", min=" << minGridNoCells << ", max=" << maxGridNoCells << ", max deviation=" << devNoGridCells
              << ", max deviation leaf=" << devNoGridLeafCells << ", total=" << globalGridNoCells << std::endl;
    for(MInt i = 0; i < noSolvers(); i++) {
      std::cerr << "Solver-statistics for solver " << i << " min= " << minNoSolverCells(i)
                << " max= " << maxNoSolverCells(i) << " total= " << globalNoSolverCells(i)
                << " avg= " << avgNoSolverCells(i) << std::endl;
    }
  }
}


/**
 * \brief performs mesh adaptation
 * \author Lennart Schneiders
 */
template <MInt nDim>
MBool Controller<nDim>::adaptation(const MBool force) {
  if(m_grid == nullptr) {
    return false;
  }

  // Reset adaptation status
  gridb().m_wasAdapted = false;

  if(!m_adaptation
     || (!force
         && ((m_adaptationInterval <= 0 || globalTimeStep % m_adaptationInterval != 0)
             || globalTimeStep < m_adaptationStart))) {
    return false;
  }

  RECORD_TIMER_START(m_timers[Timers::Adaptation]);

  cerr0 << "=== Initialising solver-adaptation at time-step " << globalTimeStep << std::endl;
  // TODO labels:TIMERS,ADAPTATION adaptation timers for new split version -> initTimers()
  logTimerStatistics("before adaptation");

  // Disable all DLB timers
  const MInt noDlbTimers = maia::dlb::g_dlbTimerController.noDlbTimers();
  ScratchSpace<MBool> dlbTimersStatus(std::max(noDlbTimers, 1), AT_, "dlbTimersStatus");
  maia::dlb::g_dlbTimerController.disableAllDlbTimers(&dlbTimersStatus[0]);

  // TODO labels:CONTROLLER remove if all solvers have splitAdaptation
  MBool splitAdaptation = true;
  for(MInt i = 0; i < noSolvers(); i++) {
    if(!solver(i).m_splitAdaptation) {
      std::cerr << "Update Adaptation to splitAdaptation! solver " << i << std::endl;
      splitAdaptation = false;
    }
  }

  // Cancel any open MPI requests in the solver since these might conflict with communication
  // during adaptation
  for(MInt i = 0; i < noSolvers(); i++) {
    solver(i).cancelMpiRequests();
  }

  if(splitAdaptation) { // new split-Adaptation version

    // prepare Adaptations
    for(MInt i = 0; i < noSolvers(); i++) {
      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Adaptation]);
      solver(i).prepareAdaptation();
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Adaptation]);
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
    }

    for(MInt i = 0; i < noCouplers(); i++) {
      const MInt timerId = noSolvers() + i;
      RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Total]);
      RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Adaptation]);
      coupler(i).prepareAdaptation();
      RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Adaptation]);
      RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Total]);
    }

    // loop over all available levels (and a specified additional number of steps)!
    // prefered over while-loop!
    MInt addedAdaptationSteps = 0;
    addedAdaptationSteps = Context::getBasicProperty<MInt>("addedAdaptationSteps", AT_, &addedAdaptationSteps);
    for(MInt level = gridb().maxUniformRefinementLevel(); level < gridb().maxRefinementLevel() + addedAdaptationSteps;
        level++) {
      const MLong oldLocalCnt = gridb().noInternalCells();
      const MLong oldCnt = gridb().noCellsGlobal();
      std::vector<std::vector<MFloat>> sensors;
      std::vector<std::bitset<CartesianGrid<nDim>::m_maxNoSensors>> sensorCellFlag(gridb().noInternalCells());
      std::vector<MFloat> sensorWeight;
      std::vector<MInt> sensorSolverId;

      // set sensors
      for(MInt i = 0; i < noSolvers(); i++) {
        solver(i).setSensors(sensors, sensorWeight, sensorCellFlag, sensorSolverId);
      }

      ASSERT(sensors.size() < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitsetsize!");

      // saveSensorData
      // Determine of at least one solver wants to save its sensor data
      MBool saveSensorData = false;
      for(MInt i = 0; i < noSolvers(); i++) {
        saveSensorData = saveSensorData || solver(i).m_saveSensorData;
      }
      if(saveSensorData) {
        // Write grid file
        std::stringstream gridFileName;
        gridFileName << "sensorDataGrid_" << std::to_string(level) << "_" << globalTimeStep << ParallelIo::fileExt();
        gridb().saveGrid((m_outputDir + gridFileName.str()).c_str(), m_recalcIds);
        // Write sensor data
        for(MInt i = 0; i < noSolvers(); i++) {
          if(solver(i).m_saveSensorData) {
            solver(i).saveSensorData(sensors, level, gridFileName.str(), m_recalcIds);
          }
        }
      }

      RECORD_TIMER_START(m_timers[Timers::MeshAdaptation]);
      gridb().meshAdaptation(sensors, sensorWeight, sensorCellFlag, sensorSolverId, m_refineCellSolver,
                             m_removeChildsSolver, m_removeCellSolver, m_swapProxySolver, m_cellOutside,
                             m_resizeGridMapSolver);
      RECORD_TIMER_STOP(m_timers[Timers::MeshAdaptation]);

      // compacts cells, and emidiate after adaptation
      for(MInt i = 0; i < noSolvers(); i++) {
        RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
        RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Adaptation]);
        solver(i).postAdaptation();
        RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Adaptation]);
        RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
      }

      // call couplers only after all solvers have rebild their grid-structure!
      for(MInt j = 0; j < noCouplers(); j++) {
        const MInt timerId = noSolvers() + j;
        RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Total]);
        RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::Adaptation]);
        coupler(j).postAdaptation();
        RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Adaptation]);
        RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::Total]);
      }

      m_log << "Mesh adaptation: " << gridb().noCellsGlobal() - oldCnt << " cells generated"
            << " (before: " << oldCnt << ", now: " << gridb().noCellsGlobal() << ")." << std::endl;

      MBool skipLoop = false;
      for(MInt i = 0; i < noSolvers(); i++) {
        if(solver(i).m_singleAdaptation) {
          skipLoop = true;
        }
      }

      if(skipLoop && globalTimeStep > 0) break;

      if(globalTimeStep < 0 && m_loadBalancingInterval > 0) {
        MInt balanceTrigger = 0;
        MInt deltaCells = gridb().noInternalCells() - oldLocalCnt;
        if(gridb().noInternalCells() + deltaCells * gridb().m_maxNoChilds > gridb().maxNoCells()) {
          balanceTrigger = 1;
        }
        MPI_Allreduce(MPI_IN_PLACE, &balanceTrigger, 1, MPI_INT, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                      "balanceTrigger");
        if(balanceTrigger || level == (gridb().maxRefinementLevel() - 1)) {
          cerr0 << "Performing intermediate load balancing step." << std::endl;
#ifdef MAIA_GRID_SANITY_CHECKS
          gridb().gridSanityChecks();
          gridb().checkWindowHaloConsistency(true);
#endif
          // preserve loadBalancing mode and
          // switch to weight based mode for the balance during initial adaptation!
          const MInt backUp = m_loadBalancingMode;
          m_loadBalancingMode = 0;
          balance(true, false, true);
          m_loadBalancingMode = backUp;
        }
      }
    }

    // reinit solver after adaptation
    for(MInt i = 0; i < noSolvers(); i++) {
      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Total]);
      RECORD_TIMER_START(m_solverTimers[i][SolverTimers::Adaptation]);
      solver(i).finalizeAdaptation();
      for(MInt j = 0; j < noCouplers(); j++) {
        coupler(j).finalizeAdaptation(i);
      }
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Adaptation]);
      RECORD_TIMER_STOP(m_solverTimers[i][SolverTimers::Total]);
    }

  } else { // TODO labels:CONTROLLER delete below
    const MLong oldCnt = gridb().noCellsGlobal();
    std::vector<std::vector<MFloat>> sensors;
    std::vector<std::bitset<CartesianGrid<nDim>::m_maxNoSensors>> sensorCellFlag(gridb().noInternalCells());
    std::vector<MFloat> sensorWeight;
    std::vector<MInt> sensorSolverId;

    for(MInt i = 0; i < noSolvers(); i++) {
      std::cerr << "preparing  adaptation for solver " << i << std::endl;
      solver(i).prepareAdaptation(sensors, sensorWeight, sensorCellFlag, sensorSolverId);
    }
    ASSERT(sensors.size() < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

    gridb().meshAdaptation(sensors, sensorWeight, sensorCellFlag, sensorSolverId, m_refineCellSolver,
                           m_removeChildsSolver, m_removeCellSolver, m_swapProxySolver, m_cellOutside,
                           m_resizeGridMapSolver);

    for(MInt i = 0; i < noSolvers(); i++) {
      // in multisolver case: translate tree ids back to solver-local ids
      solver(i).reinitAfterAdaptation();
    }

    m_log << "Mesh adaptation: " << gridb().noCellsGlobal() - oldCnt << " cells generated"
          << " (before: " << oldCnt << ", now: " << gridb().noCellsGlobal() << ")." << std::endl;
  }

  // Enable all previously active DLB timers again
  maia::dlb::g_dlbTimerController.enableAllDlbTimers(&dlbTimersStatus[0]);

  // reset all timers after an adaptation:
  // Tim version
  // resetAllTimer();

  // set the last adaptation timeStep
  m_lastAdaptationTimeStep = globalTimeStep;
  m_nAdaptationsSinceBalance = m_nAdaptationsSinceBalance + 1;

  printDomainStatistics("after adaptation");

  cerr0 << "=== Finished adaptation" << std::endl;
  m_log << "Finished adaptation at time-step " << globalTimeStep << std::endl;
  RECORD_TIMER_STOP(m_timers[Timers::Adaptation]);

  return true;
}

/**
 * \brief This function handels the setCellWeights-functionality of all solvers and
 *        writes the accumulated-cellweight of all solvers into the grid.
 * \author Tim Wegmann
 */
template <MInt nDim>
void Controller<nDim>::accumulateCellWeights() {
  TRACE();

  const MInt noGridCells = gridb().treeb().size();

  // Reset cell weights
  for(MInt gridCellId = 0; gridCellId < noGridCells; gridCellId++) {
    gridb().treeb().weight(gridCellId) = 0.0;
  }

  if(!m_dlbRestartWeights || m_dlbLastWeights == nullptr) {
    // use static cell weights

    // TODO labels:CONTROLLER why not use an array of size noGridCells, set to zero
    // and add cell weights of each solver?
    const MInt sizeTimesSolver = noGridCells * noSolvers();
    MFloatScratchSpace solverweight(sizeTimesSolver, FUN_, "solverweight");
    // might be necessary for dimensionless weight-comparison between different-solvers,
    // to devide the weight of ech solver-cell by its maximum-possible-value!
    // MFloatScratchSpace solvermaxweights(noSolvers(),FUN_, "solvermaxweights");
    solverweight.fill(F0);

    for(MInt i = 0; i < noSolvers(); i++) {
      if(solver(i).isActive()) { // Only for solvers active on this ranks
        solver(i).setCellWeights(&solverweight[0]);
      }
    }

    // Accumulate weights of all solvers
    for(MInt i = 0; i < noSolvers(); i++) {
      for(MInt gridCellId = 0; gridCellId < noGridCells; gridCellId++) {
        const MInt id = gridCellId + noGridCells * i;
        gridb().treeb().weight(gridCellId) += solverweight[id];
      }
    }


  } else {
    // use weights from last DLB
    MInt weightOffset = 0;
    for(MInt solverId = 0; solverId < noSolvers(); solverId++) {
      if(solver(solverId).isActive()) {
        for(MInt cellId = 0; cellId < noGridCells; cellId++) {
          if(gridb().a_hasProperty(cellId, Cell::IsHalo)) {
            continue;
          }
          gridb().a_weight(cellId) += solver(solverId).getCellLoad(cellId, &m_dlbLastWeights[weightOffset]);
        }
      }
      weightOffset += solver(solverId).noLoadTypes();
    }
  }
}


/// \brief Determine the cell weights using the given weighting factors for the different solver load
///        quantities and update the partition cell workloads.
template <MInt nDim>
void Controller<nDim>::updateWeightsAndWorkloads(const std::vector<MFloat>& weights,
                                                 const MBool restoreDefaultWeights) {
  TRACE();

  const MInt noGridCells = gridb().treeb().size();

  // Restore default cell weights: each cell has weight 1.
  // else reset all weights to 0
  const MFloat initWeight = (restoreDefaultWeights) ? 1.0 : 0.0;
  for(MInt gridCellId = 0; gridCellId < noGridCells; gridCellId++) {
    gridb().a_weight(gridCellId) = initWeight;
  }

  if(!restoreDefaultWeights) {
    // Update the cell weights
    MInt weightOffset = 0;
    for(MInt solverId = 0; solverId < noSolvers(); solverId++) {
      if(solver(solverId).isActive()) { // Only add cell loads for active solvers
        for(MInt cellId = 0; cellId < noGridCells; cellId++) {
          if(gridb().a_hasProperty(cellId, Cell::IsHalo)) continue;
          gridb().a_weight(cellId) += solver(solverId).getCellLoad(cellId, &weights[weightOffset]);
        }
      }
      weightOffset += solver(solverId).noLoadTypes();
    }
  }

  // Calculate the partition-cell workloads
  gridb().calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), noGridCells);
}


/// \brief Update the partition cell workloads in the grid file using user specified or default
///        weights for the solver load quantities.
template <MInt nDim>
MBool Controller<nDim>::updateGridPartitionWorkloads() {
  TRACE();

  // Check if update of partition workloads is enabled
  MBool updateGridPartitionWorkloads = false;
  updateGridPartitionWorkloads =
      Context::getBasicProperty<MBool>("updateGridPartitionWorkloads", AT_, &updateGridPartitionWorkloads);
  if(!updateGridPartitionWorkloads) {
    return false;
  }

  m_log << "Updating partition cell workloads and save them to grid... " << std::endl;

  MBool restore = false;
  restore = Context::getBasicProperty<MBool>("restoreDefaultWorkloads", AT_, &restore);

  // Number of load quantities (should be the same for every rank)
  const MInt noLoadTypes = globalNoLoadTypes();
  std::vector<MFloat> weights{};

  if(noLoadTypes < 1 || restore) {
    restore = true;
    m_log << "Restoring default weights" << std::endl;
  } else {
    m_log << "Using specified solver weights" << std::endl;
    // Get the weights for all solvers
    getSpecifiedSolverWeights(weights);
  }

  // Use the determined weights to update the partition workloads (or restore default weighting)
  updateWeightsAndWorkloads(weights, restore);
  // Save new partition cell workloads to grid file
  gridb().savePartitionCellWorkloadsGridFile();
  m_log << "done" << std::endl;

  if(domainId() == 0) {
    std::cout << "Updated cell weights. Restart solver with 'updateGridPartitionWorkloads = false'" << std::endl;
  }
  return true;
}

/// \brief Return the specified (or default) solver weights for all solvers
template <MInt nDim>
void Controller<nDim>::getSpecifiedSolverWeights(std::vector<MFloat>& weights) {
  TRACE();

  // Determine the weights for all solvers
  const MInt noLoadTypes = globalNoLoadTypes();
  weights.resize(noLoadTypes);
  std::fill(weights.begin(), weights.end(), 0.0);

  MInt offset = 0;
  for(MInt i = 0; i < noSolvers(); i++) {
    const MInt solverCount = solver(i).noLoadTypes();
    std::vector<MString> names(solverCount);

    // Get default weights (+names) for this solver
    solver(i).getDefaultWeights(&weights[offset], names);

    // Check for user specified solver weights
    const MString propName = "solverWeights_" + std::to_string(i);
    if(Context::propertyExists(propName)) {
      if(Context::propertyLength(propName) != solverCount) {
        TERMM(1, "wrong length of property '" + propName + "', should be of length " + std::to_string(solverCount));
      }
      for(MInt j = 0; j < solverCount; j++) {
        weights[offset + j] = Context::getBasicProperty<MFloat>(propName, AT_, j);
      }
    }

    for(MInt j = 0; j < solverCount; j++) {
      m_log << "Solver #" << i << ", weight #" << j << ": " << names[j] << ", " << weights[offset + j] << std::endl;
    }

    offset += solverCount;
  }
}


template <MInt nDim>
void Controller<nDim>::writeRestartFile(const MBool forceRestart, const MBool writeBackup) {
  if(m_grid == nullptr) {
    return;
  }
  MBoolScratchSpace writeRestart(noSolvers(), FUN_, "writeRestart");
  MBoolScratchSpace writeGridRestart(noSolvers(), FUN_, "writeGridRestart");

  // disable load timers:
  std::vector<MInt> dlbTimersStatus(noSolvers());

  for(MInt i = 0; i < noSolvers(); i++) {
    dlbTimersStatus[i] = solver(i).dlbTimersEnabled();
    if(dlbTimersStatus[i] != 0) {
      solver(i).disableDlbTimers(); // Disable DLB timers if active
    }
  }

  for(MInt i = 0; i < noSolvers(); i++) {
    // prepareRestart() and reInitAfterRestart must be called for all ranks!
    writeRestart[i] = solver(i).prepareRestart(forceRestart, writeGridRestart[i]);
  }

  MBool allsolversrestart = false;
  MBool gridrestart = false;

  for(MInt i = 0; i < noSolvers(); i++) {
    if(writeRestart[i]) allsolversrestart = true;
    if(writeRestart[i] && writeGridRestart[i]) gridrestart = true;
  }

  gridb().deletePeriodicConnection(true);

  if(gridrestart || (allsolversrestart && m_saveGridNewPartitionCells)) {
    accumulateCellWeights();

    {
      std::stringstream s;
      s << "restartGrid";
      if(writeBackup) s << "Backup";
      if(!m_useNonSpecifiedRestartFile || writeBackup) s << "_00" << globalTimeStep;
      // TODO labels:CONTROLLER,IO,totest change to the version below, fix testcases
      // if (!m_useNonSpecifiedRestartFile || writeBackup) {
      //  s << "_" << setw(8) << setfill('0') << globalTimeStep;
      //}
      s << ParallelIo::fileExt();
      m_currentGridFileName = s.str();

      // Cancel any open MPI requests in the solver since these might conflict with communication
      // necessary when writing the grid file
      for(MInt i = 0; i < noSolvers(); i++) {
        solver(i).cancelMpiRequests();
      }

      cerr0 << "Saving adapted grid file for all solvers " << s.str() << "... ";

      gridb().saveGrid((m_outputDir + m_currentGridFileName).c_str(), m_recalcIds);

      cerr0 << "ok." << std::endl;

      m_saveGridNewPartitionCells = false;

    }
  }

  for(MInt c = 0; c < noCouplers(); c++) {
    coupler(c).writeRestartFile(allsolversrestart);
  }

  if(allsolversrestart) {
    cerr0 << "Writing restart files at time step " << globalTimeStep << " :" << std::endl;

    for(MInt i = 0; i < noSolvers(); i++) {
      if(solver(i).isActive()) {
        solver(i).writeRestartFile(writeRestart[i], writeBackup, m_currentGridFileName, m_recalcIds);
      }
    }
    // savePartitionFile();
  }

  for(MInt i = 0; i < noSolvers(); i++) {
    solver(i).reIntAfterRestart(writeRestart[i]);
  }

  gridb().restorePeriodicConnection();

  // reanable all previously running timers:
  for(MInt i = 0; i < noSolvers(); i++) {
    if(dlbTimersStatus[i] != 0) {
      solver(i).enableDlbTimers(); // Reenable previously active DLB timers
    }
  }
}


/// \brief Return global number of load types of all solvers
template <MInt nDim>
MInt Controller<nDim>::globalNoLoadTypes() {
  MInt noLoadTypes = 0;
  for(MInt i = 0; i < noSolvers(); i++) {
    noLoadTypes += solver(i).noLoadTypes();
  }
  return noLoadTypes;
}


/// \brief Return load quantities of all solvers on this domain
template <MInt nDim>
void Controller<nDim>::getLoadQuantities(MInt* const loadQuantities) {
  // Fill with zeros since inactive solvers will not write anything
  std::fill_n(&loadQuantities[0], globalNoLoadTypes(), 0);

  MInt offset = 0;
  for(MInt i = 0; i < noSolvers(); i++) {
    solver(i).getLoadQuantities(&loadQuantities[offset]);
    offset += solver(i).noLoadTypes();
  }
}


/// \brief Return if dynamic load balancing is needed and compute domain loads.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] localRunTime The local run time, i.e. computation time.
/// \param[out] load Computed domain loads.
/// \param[out] imbalance Computed imbalance.
template <MInt nDim>
MBool Controller<nDim>::needLoadBalancing(const MFloat localRunTime, MFloat* const loads, MFloat& imbalance) {
  TRACE();

  // Nothing to do in serial
  if(noDomains() == 1) {
    return false;
  }

  MFloatScratchSpace localRunTimes(noDomains(), FUN_, "localRunTimes");
  MFloatScratchSpace localIdleTimes(noDomains(), FUN_, "localIdleTimes");

  // Gather the local runtimes on all ranks
  MPI_Allgather(&localRunTime, 1, type_traits<MFloat>::mpiType(), &localRunTimes[0], 1, type_traits<MFloat>::mpiType(),
                gridb().mpiComm(), AT_, "localRunTime", "localRunTimes[0]");

  // Calculate average and maximum runtime of all ranks
  const MFloat averageRunTime = std::accumulate(localRunTimes.begin(), localRunTimes.end(), 0.0) / noDomains();
  const MFloat maxRunTime = *std::max_element(localRunTimes.begin(), localRunTimes.end());

  // Imbalance percentage: I_% = (t_max - t_avg)/t_max * N/(N-1)
  // (DeRose 2007, Detecting Application Load Imbalance on High End Massively Parallel Systems)
  imbalance = (maxRunTime - averageRunTime) / maxRunTime * noDomains() / (noDomains() - 1.0);
  const MBool loadBalance = (imbalance > m_dlbImbalanceThreshold);

  // Compute domain loads as fraction of local run time to average run time
  for(MInt i = 0; i < noDomains(); i++) {
    const MFloat load = localRunTimes[i] / averageRunTime;
    loads[i] = load;
  }

  const MFloat maxLoad = *std::max_element(&loads[0], &loads[0] + noDomains());
  const MFloat minLoad = *std::min_element(&loads[0], &loads[0] + noDomains());

  char imb[10];
  std::sprintf(imb, "%.2f", imbalance * 100.0);
  std::stringstream message;
  message << globalTimeStep << " * Imbalance percentage: " << imb << "%, t_avg = " << averageRunTime
          << ", t_max = " << maxRunTime << std::endl;
  message << globalTimeStep << " * Loads: max = " << maxLoad << ", min = " << minLoad << std::endl;
  m_log << message.str();
  cerr0 << message.str();

  // TODO labels:CONTROLLER,DLB generalize histogram output for other purposes?
  const MFloat binWidth = 0.05;
  const MInt noBins = std::ceil(maxLoad / binWidth) + 1;
  std::vector<MInt> loadBins(noBins);
  std::fill(loadBins.begin(), loadBins.end(), 0);
  for(MInt i = 0; i < noDomains(); i++) {
    const MInt binId = std::floor(loads[i] / binWidth);
    loadBins[binId] += 1;
  }

  const MInt sumBins = std::accumulate(loadBins.begin(), loadBins.end(), 0);
  if(sumBins != noDomains()) {
    m_log << "ERROR in load bins, count " << sumBins << " != " << noDomains() << std::endl;
  }

  const MInt noBinsPrev = m_previousLoadBins.size();
  const MInt maxNoBins = std::max(noBins, noBinsPrev);
  const MInt maxBinCountCurr = *std::max_element(loadBins.begin(), loadBins.end());
  const MInt maxBinCountPrev =
      (noBinsPrev > 0) ? *std::max_element(m_previousLoadBins.begin(), m_previousLoadBins.end()) : 0;
  const MInt maxBinCount = std::max(maxBinCountCurr, maxBinCountPrev);

  const MInt firstBin =
      (std::find_if(loadBins.begin(), loadBins.end(), [](MInt i) { return i > 0; }) - loadBins.begin()) - 1;
  const MInt firstBinPrev =
      (std::find_if(m_previousLoadBins.begin(), m_previousLoadBins.end(), [](MInt i) { return i > 0; })
       - m_previousLoadBins.begin())
      - 1;

  const MInt minFirstBin = (noBinsPrev > 0) ? std::max(std::min(firstBin, firstBinPrev), 0) : std::max(firstBin, 0);

  // Set maximum line length
  const MInt maxLineLength = 256;
  const MString dashes(maxLineLength, '-');
  // Create sprintf buffer and string stream
  MChar b[maxLineLength];
  std::stringstream hist;

  snprintf(b, maxLineLength, " |%.126s|\n", dashes.c_str());
  MString separatorTop(b);

  snprintf(b, maxLineLength, " |%.18s|%.42s|%.42s|%.21s|\n", dashes.c_str(), dashes.c_str(), dashes.c_str(),
           dashes.c_str());
  MString separator(b);

  // Create header
  snprintf(b, maxLineLength,
           " | Load distribution at global timestep %-8d%*s - imbalance %5.2f%% - t_avg %5.3e - t_max %5.3e - loads "
           "[%5.3f,%5.3f] |\n",
           globalTimeStep, 2, " ", imbalance * 100.0, averageRunTime, maxRunTime, minLoad, maxLoad);
  hist << "\n" << separatorTop << b;

  if(noBinsPrev > 0) {
    snprintf(b, maxLineLength,
             " | Load distribution at global timestep %-8d%*s - imbalance %5.2f%% - t_avg %5.3e - t_max %5.3e - loads "
             "[%5.3f,%5.3f] |\n",
             m_previousLoadInfoStep, 2, " ", m_previousLoadInfo[0], m_previousLoadInfo[1], m_previousLoadInfo[2],
             m_previousLoadInfo[3], m_previousLoadInfo[4]);
    hist << separatorTop << b;
  }

  snprintf(b, maxLineLength, " |     load bin     | %-32s%8d | %-32s%8d |   curr/prev count   |\n",
           "current distribution - timestep", globalTimeStep, "previous distribution - timestep",
           m_previousLoadInfoStep);
  hist << separator << b << separator;

  MInt checksum = 0, checksumPrev = 0;
  // Output horizontal histogram
  for(MInt i = maxNoBins - 1; i >= minFirstBin; i--) {
    const MInt maxBarWidth = 40;
    const MInt count = (i < noBins) ? loadBins[i] : 0;
    const MInt countPrev = (i < noBinsPrev) ? m_previousLoadBins[i] : 0;
    // Display at least one char if the count is > 0
    const MInt width = (count > 0) ? std::max((MInt)(maxBarWidth * (MFloat)count / (MFloat)maxBinCount), 1) : 0;
    const MInt widthPrev =
        (countPrev > 0) ? std::max((MInt)(maxBarWidth * (MFloat)countPrev / (MFloat)maxBinCount), 1) : 0;
    const MString bar(width, '#');
    const MString barPrev(widthPrev, '@');

    if(approx(i * binWidth, 1.0, 0.0001)) {
      // Add horizontal separator marks between overloaded/underloaded ranks
      snprintf(b, maxLineLength, " | [%6.3f, %-6.3f) |_%-40s_|_%-40s_| %8d | %8d |\n", i * binWidth, (i + 1) * binWidth,
               bar.c_str(), barPrev.c_str(), count, countPrev);
    } else {
      snprintf(b, maxLineLength, " | [%6.3f, %-6.3f) | %-40s | %-40s | %8d | %8d |\n", i * binWidth, (i + 1) * binWidth,
               bar.c_str(), barPrev.c_str(), count, countPrev);
    }
    hist << b;

    // Make sure all ranks appear in the histogram
    checksum += count;
    checksumPrev += countPrev;
  }
  hist << separator;

  if(checksum != noDomains() || (checksumPrev != noDomains() && noBinsPrev > 0)) {
    m_log << "ERROR in load histogram" << std::endl;
    cerr0 << "ERROR in load histogram" << std::endl;
  } else {
    m_log << hist.str() << std::endl;
    cerr0 << hist.str() << std::endl;
  }

  // Store current information for next output
  m_previousLoadBins = loadBins;
  m_previousLoadInfoStep = globalTimeStep;
  m_previousLoadInfo[0] = imbalance * 100.0;
  m_previousLoadInfo[1] = averageRunTime;
  m_previousLoadInfo[2] = maxRunTime;
  m_previousLoadInfo[3] = minLoad;
  m_previousLoadInfo[4] = maxLoad;

  return loadBalance;
}


/// \brief Calculate new global domain offsets given the current and new global partition cell
///        offsets.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] oldPartitionCellOffsets The current (old) partition cell offsets.
/// \param[in] newPartitionCellOffsets The new partition cell offsets.
/// \param[out] globalOffsets Pointer to storage of size noDomains()+1 that will hold the new global
///             domain offsets on exit.
template <MInt nDim>
void Controller<nDim>::loadBalancingCalcNewGlobalOffsets(const MLong* const oldPartitionCellOffsets,
                                                         const MLong* const newPartitionCellOffsets,
                                                         MLong* const globalOffsets) {
  TRACE();

  const MInt oldNoPartitionCells = gridb().m_noPartitionCells;
  MLongScratchSpace partitionCellsGlobalId(oldNoPartitionCells, AT_, "partitionCellsGlobalId");

  for(MInt i = 0; i < oldNoPartitionCells; i++) {
    const MLong globalCellId = gridb().m_localPartitionCellGlobalIds[i]; // sorted by global id

    partitionCellsGlobalId(i) = globalCellId;
    ASSERT(i == 0 || partitionCellsGlobalId(i) > partitionCellsGlobalId(i - 1),
           "Partition cell global ids not in ascending order.");
  }

  // Fill with zeros
  std::fill_n(&globalOffsets[0], noDomains() + 1, 0);

  // Calculate global offsets for the given new partition cell offsets
  for(MInt i = 1; i < noDomains(); i++) {
    const MLong partitionCellId = newPartitionCellOffsets[i];
    // Determine if this partition cell is currently on this domain
    const MBool hasPartitionCell = (oldPartitionCellOffsets[domainId()] <= partitionCellId
                                    && partitionCellId < oldPartitionCellOffsets[domainId() + 1]);

    // Set global domain offset if the partition cell is present, else it remains zero
    if(hasPartitionCell) {
      // Determine local partition cell index and the global cell id
      const MLong partitionCellLocalId = partitionCellId - oldPartitionCellOffsets[domainId()];
      globalOffsets[i] = partitionCellsGlobalId[partitionCellLocalId];

      // Local cell id of the partition cell and its parent
      MInt currentId = gridb().m_localPartitionCellLocalIds[partitionCellLocalId];
      MInt parentId = gridb().a_parentId(currentId);

      MInt shift = 0;
      // Partition level shift: if the partition cell is not on the min-level check the global id of
      // its parent, if 'globalId(parent) == globalId(partitionCell) - 1' the parent is a partition
      // level ancestor that belongs to this domain, i.e., it has no offspring cells on the previous
      // domain. Continue to loop up the parent cells and check this condition to find the shift
      // that yields the correct domain offset.
      while(gridb().a_level(currentId) != gridb().minLevel()
            && gridb().a_globalId(parentId) == gridb().a_globalId(currentId) - 1) {
        shift++;
        currentId = parentId;
        parentId = gridb().a_parentId(currentId);
      }

      TERMM_IF_COND(shift > 0 && gridb().m_maxPartitionLevelShift == 0,
                    "Error: domain offset has a shift but the maximum partition level shift is 0.");

      // Correct global domain offset
      globalOffsets[i] -= shift;
    }
  }

  // Set only once since distributed via allreduce with sum to all domains
  if(domainId() == 0) {
    globalOffsets[noDomains()] = gridb().domainOffset(noDomains());
  }

  // Combine the new domain offset information of all ranks, i.e. each new
  // global domain offset is determined by the domain where the corresponding
  // partition cell is present, other domains will add a value of zero to this offset
  MPI_Allreduce(MPI_IN_PLACE, &globalOffsets[0], noDomains() + 1, type_traits<MLong>::mpiType(), MPI_SUM,
                gridb().mpiComm(), AT_, "MPI_IN_PLACE", "globalOffsets[0]");
}


/// \brief Based on new domain offsets calculate the number of cells to
///        send/receive to/from each domain.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] offsets New domain offsets.
/// \param[out] noCellsToSendByDomains Holds the number of cells to send to each domain.
/// \param[out] noCellsToRecvByDomain Holds the number of cells to receive from each domain.
/// \param[out] sortedCellId Buffer index for each cell, such that data is sorted by global id.
/// \param[out] bufferIdToCellId Reversed map of sortedCellId.
template <MInt nDim>
void Controller<nDim>::loadBalancingCalcNoCellsToSend(const MLong* const offsets,
                                                      MInt* const noCellsToSendByDomain,
                                                      MInt* const noCellsToReceiveByDomain,
                                                      MInt* const sortedCellId,
                                                      MInt* const bufferIdToCellId) {
  TRACE();

  const MInt noCells = gridb().treeb().size();

  MIntScratchSpace targetDomainsByCell(noCells, AT_, "targetDomainsByCell");
  std::fill_n(&targetDomainsByCell[0], noCells, -1);
  std::fill_n(&noCellsToSendByDomain[0], noDomains() + 1, 0);

  // Loop over all cells and determine to which domain they belong under the new domain
  // decomposition
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    // Skip halo cells
    if(gridb().a_hasProperty(cellId, Cell::IsHalo)) {
      continue;
    }

    const MLong globalCellId = gridb().a_globalId(cellId);
    // Search for global cell id in domain offsets array
    auto lowerBound = std::lower_bound(&offsets[0], &offsets[0] + noDomains(), globalCellId);
    const MInt dist = std::distance(&offsets[0], lowerBound);
    // Check if this cell is a domain offset (i.e. in the offsets list)
    const MBool isDomainOffset = (*lowerBound == globalCellId);
    // Determine target domain id
    const MInt globalDomain = (isDomainOffset) ? dist : dist - 1;

    // Tag cell with target domain
    targetDomainsByCell[cellId] = globalDomain;
    // Increase number of cells to send to this domain
    noCellsToSendByDomain[globalDomain]++;
  }

  // Use all-to-all communication to determine how many cells to send/receive to/from each domain
  MPI_Alltoall(&noCellsToSendByDomain[0], 1, type_traits<MInt>::mpiType(), &noCellsToReceiveByDomain[0], 1,
               type_traits<MInt>::mpiType(), gridb().mpiComm(), AT_, "noCellsToSendByDomain[0]",
               "noCellsToReceiveByDomain[0]");

  // Determine total number of cells to send/recv (store as last entry)
  noCellsToSendByDomain[noDomains()] =
      std::accumulate(&noCellsToSendByDomain[0], &noCellsToSendByDomain[0] + noDomains(), 0);
  noCellsToReceiveByDomain[noDomains()] =
      std::accumulate(&noCellsToReceiveByDomain[0], &noCellsToReceiveByDomain[0] + noDomains(), 0);

  if(noCellsToSendByDomain[noDomains()] < 1) {
    TERMM(1, std::to_string(domainId()) + " noCellsToSend = " + std::to_string(noCellsToSendByDomain[noDomains()]));
  }
  if(noCellsToReceiveByDomain[noDomains()] < 1) {
    TERMM(1, std::to_string(domainId()) + " noCellsToRecv = " + std::to_string(noCellsToReceiveByDomain[noDomains()]));
  }

  // Sort cells depending on the domain to send to
  std::fill_n(&sortedCellId[0], noCells, -1);

  // @ansgar TODO labels:CONTROLLER,DLB this might be very inefficient for large number of cores
  MInt currentBufferId = 0;
  for(MInt dom = 0; dom < noDomains(); ++dom) {
    // Map from global cell id to local cell id for cells that will belong to the current domain
    std::map<MLong, MInt> cellMap;
    for(MInt cellId = 0; cellId < noCells; ++cellId) {
      // Halo cells are automatically skipped as they have tag[] of -1.
      // Halo cells will have a sortedCellId[] of -1.
      if(targetDomainsByCell[cellId] == dom) {
        cellMap[gridb().a_globalId(cellId)] = cellId;
      }
    }

    // Set buffer location for all cells in global id order -> received data is stored in
    // global id order and does not need to be resorted
    for(auto&& cell : cellMap) {
      const MInt cellId = cell.second;
      sortedCellId[cellId] = currentBufferId;
      bufferIdToCellId[currentBufferId] = cellId; // map from buffer id to cell id
      ++currentBufferId;
    }
  }
}


/// Communicate the solver variables which should be the same on all ranks from the solver local rank
/// 0 to all other domains. This is needed if during load balancing the ranks of a solver change,
/// i.e. ranks becoming active/inactive due to partitioning changes
template <MInt nDim>
void Controller<nDim>::communicateGlobalSolverVars(Solver* const solver) {
  TRACE();

  // Determine solver local root domain
  const MInt localRootGlobalDomain = solverLocalRootDomain(solver);

  std::vector<MInt> globalIdVars(0);
  std::vector<MFloat> globalFloatVars(0);

  // Request global solver variables
  solver->getGlobalSolverVars(globalFloatVars, globalIdVars);

  const MInt noIdVars = globalIdVars.size();
  const MInt noFloatVars = globalFloatVars.size();

  // Distribute all global variables from local rank 0 to all domains
  MPI_Bcast(&globalIdVars[0], noIdVars, type_traits<MInt>::mpiType(), localRootGlobalDomain, gridb().mpiComm(), AT_,
            "globalIdVars[0]");
  MPI_Bcast(&globalFloatVars[0], noFloatVars, type_traits<MFloat>::mpiType(), localRootGlobalDomain, gridb().mpiComm(),
            AT_, "globalFloatVars[0]");

  // Set variables in solver
  solver->setGlobalSolverVars(globalFloatVars, globalIdVars);
}


/// Determine the data sizes for a solver/coupler that need to be communicated during load balancing
template <MInt nDim>
void Controller<nDim>::determineDataSizesDlb(const MInt id, const MInt mode, const MInt* const noCellsToSend,
                                             const MInt* const bufferIdToCellId,
                                             std::vector<std::vector<MInt>>& sendSizeVector,
                                             std::vector<std::vector<MInt>>& recvSizeVector) {
  TRACE();

  const MBool isSolver = (mode == 0); // mode 0: solver; mode 1: coupler
  const MInt timerId = (isSolver) ? id : noSolvers() + id;
  // Get the number of quantities that need to be communicated
  MInt dataCount = (isSolver) ? solver(id).noCellDataDlb() : coupler(id).noCellDataDlb();

  RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::CalcDataSizesMpiBlocking]);
  // Note: since inactive ranks might not return the correct count, take the max among all domains
  MPI_Allreduce(MPI_IN_PLACE, &dataCount, 1, type_traits<MInt>::mpiType(), MPI_MAX, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "dataCount");
  RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::CalcDataSizesMpiBlocking]);

  std::vector<MInt> dataSendSize(dataCount * (globalNoDomains() + 1));
  std::vector<MInt> dataRecvSize(dataCount * (globalNoDomains() + 1));

  // Determine data size to send/receive to/from each domain for all quantities
  // last entry holds the sum
  for(MInt dataId = 0; dataId < dataCount; dataId++) {
    const MInt offset = dataId * (globalNoDomains() + 1);
    MInt bufferId = 0;

    // Loop over all domains and determine the data send size for this domain
    for(MInt domain = 0; domain < globalNoDomains(); domain++) {
      MInt domainDataSize = 0;

      // Loop over all grid cells that are send to this domain
      for(MInt i = 0; i < noCellsToSend[domain]; i++) {
        // Inquire the data size for this cell
        // TODO labels:DLB @ansgar inquire all data sizes at once?
        const MInt cellId = bufferIdToCellId[bufferId];
        ASSERT(cellId > -1, "");
        const MInt dataSize =
            (isSolver) ? solver(id).cellDataSizeDlb(dataId, cellId) : coupler(id).cellDataSizeDlb(dataId, cellId);
        domainDataSize += dataSize;
        bufferId++;
      }

      // Store total data size to send to this domain for this quantity
      dataSendSize[offset + domain] = domainDataSize;
    }

    // TODO labels:DLB @ansgar perform just one alltoall for all data ids?
    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::CalcDataSizesMpi]);
    // Exchange data send size with all other domains and store as receive size
    MPI_Alltoall(&dataSendSize[offset], 1, type_traits<MInt>::mpiType(), &dataRecvSize[offset], 1,
                 type_traits<MInt>::mpiType(), gridb().mpiComm(), AT_, "dataSendSize[offset]", "dataRecvSize[offset]");
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::CalcDataSizesMpi]);

    // Compute the total data send size and store as last entry
    dataSendSize[offset + globalNoDomains()] =
        std::accumulate(&dataSendSize[offset], &dataSendSize[offset + globalNoDomains()], 0);

    // Compute the total data receive size and store as last entry
    dataRecvSize[offset + globalNoDomains()] =
        std::accumulate(&dataRecvSize[offset], &dataRecvSize[offset + globalNoDomains()], 0);
  }

  sendSizeVector.push_back(dataSendSize);
  recvSizeVector.push_back(dataRecvSize);
}


/// Communicate all solver data for load balancing according to the send/recv sizes
template <MInt nDim>
void Controller<nDim>::redistributeDataDlb(const MInt id, const MInt mode, std::vector<MInt>& dataSendSize,
                                           std::vector<MInt>& dataRecvSize, const MInt* const bufferIdToCellId,
                                           const MInt oldNoCells, std::vector<MInt*>& intDataRecv,
                                           std::vector<MLong*>& longDataRecv, std::vector<MFloat*>& floatDataRecv,
                                           std::vector<MInt>& dataTypes) {
  TRACE();

  const MBool isSolver = (mode == 0); // mode 0: solver; mode 1: coupler
  const MInt timerId = (isSolver) ? id : noSolvers() + id;
  // Solver pointer for easier access
  Solver* const solverP = (isSolver) ? &solver(id) : nullptr;
  Coupling* const couplerP = (isSolver) ? nullptr : &coupler(id);

  RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::RedistributeMpiBlocking]);

  // Determine solver local root domain
  // TODO labels:CONTROLLER,DLB @ansgar_coupler check this!
  const MInt localRootGlobalDomain = (isSolver) ? solverLocalRootDomain(solverP) : 0;

  // Get the number of quantities that need to be communicated
  MInt dataCount = (isSolver) ? solverP->noCellDataDlb() : couplerP->noCellDataDlb();

  // Note: since inactive ranks might not return the correct count, take the max among all domains
  MPI_Allreduce(MPI_IN_PLACE, &dataCount, 1, type_traits<MInt>::mpiType(), MPI_MAX, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "dataCount");
  RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::RedistributeMpiBlocking]);

  MInt intDataCount = 0;
  MInt longDataCount = 0;
  MInt floatDataCount = 0;

  // Get solver/solver data and exchange
  for(MInt dataId = 0; dataId < dataCount; dataId++) {
    const MInt offset = dataId * (globalNoDomains() + 1);
    // Inquire data type (communicate from local rank 0 to avoid errors with inactive ranks) and
    // total data send size for this quantity
    // MInt dataType = (solverP->domainId() == 0) ? solverP->cellDataTypeDlb(dataId) : -1;
    MInt dataType = -1;
    if(isSolver) {
      dataType = (solverP->domainId() == 0) ? solverP->cellDataTypeDlb(dataId) : -1;
    } else {
      // TODO labels:CONTROLLER,DLB @ansgar_coupler
      /* dataType = (couplerP->domainId() == 0) ? couplerP->cellDataTypeDlb(dataId) : -1; */
      dataType = couplerP->cellDataTypeDlb(dataId);
    }

    RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);
    MPI_Bcast(&dataType, 1, type_traits<MInt>::mpiType(), localRootGlobalDomain, gridb().mpiComm(), AT_, "dataType");
    RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);

    const MInt dataSize = dataSendSize[offset + globalNoDomains()];
    const MInt recvSize = dataRecvSize[offset + globalNoDomains()];

    // Check data type (MInt, MLong or MFloat for now)
    switch(dataType) {
      case MINT: {
        dataTypes.push_back(MINT);
        // Allocate (persistent) data receive buffer and store pointer
        MInt* newIntData = nullptr;
        mAlloc(newIntData, std::max(recvSize, 1), "newIntData", AT_);
        intDataRecv.push_back(newIntData);

        // Data send buffer, temporary since not needed anymore after exchange
        ScratchSpace<MInt> intDataSend(std::max(dataSize, 1), AT_, "intDataSend");

        if(dataSize > 0) {
          // Get the solver/solver data and store in send buffer (sorted according to buffer mapping)
          if(isSolver) {
            solverP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &intDataSend[0]);
          } else {
            couplerP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &intDataSend[0]);
          }
        }

        RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);
        // Exchange
        maia::mpi::exchangeData(&intDataSend[0], globalDomainId(), globalNoDomains(), gridb().mpiComm(), 1,
                                &dataSendSize[offset], &dataRecvSize[offset], &intDataRecv[intDataCount][0]);
        RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);

        intDataCount++;
        break;
      }
      case MLONG: {
        dataTypes.push_back(MLONG);
        // Allocate (persistent) data receive buffer and store pointer
        MLong* newLongData = nullptr;
        mAlloc(newLongData, std::max(recvSize, 1), "newLongData", AT_);
        longDataRecv.push_back(newLongData);

        // Data send buffer, temporary since not needed anymore after exchange
        ScratchSpace<MLong> longDataSend(std::max(dataSize, 1), AT_, "longDataSend");

        if(dataSize > 0) {
          // Get the solver/solver data and store in send buffer (sorted according to buffer mapping)
          if(isSolver) {
            solverP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &longDataSend[0]);
          } else {
            couplerP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &longDataSend[0]);
          }
        }

        RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);
        // Exchange
        maia::mpi::exchangeData(&longDataSend[0], globalDomainId(), globalNoDomains(), gridb().mpiComm(), 1,
                                &dataSendSize[offset], &dataRecvSize[offset], &longDataRecv[longDataCount][0]);
        RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);

        longDataCount++;
        break;
      }
      case MFLOAT: {
        dataTypes.push_back(MFLOAT);
        // Allocate (persistent) data receive buffer and store pointer
        MFloat* newFloatData = nullptr;
        mAlloc(newFloatData, std::max(recvSize, 1), "newFloatData", AT_);
        floatDataRecv.push_back(newFloatData);

        // Data send buffer, temporary since not needed anymore after exchange
        ScratchSpace<MFloat> floatDataSend(std::max(dataSize, 1), AT_, "floatDataSend");

        if(dataSize > 0) {
          // Get the solver/solver data and store in send buffer
          if(isSolver) {
            solverP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &floatDataSend[0]);
          } else {
            couplerP->getCellDataDlb(dataId, oldNoCells, &bufferIdToCellId[0], &floatDataSend[0]);
          }
        }

        RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);
        // Exchange
        maia::mpi::exchangeData(&floatDataSend[0], globalDomainId(), globalNoDomains(), gridb().mpiComm(), 1,
                                &dataSendSize[offset], &dataRecvSize[offset], &floatDataRecv[floatDataCount][0]);
        RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::RedistributeMpi]);

        floatDataCount++;
        break;
      }
      default:
        TERMM(1, "Unknown data type: " + std::to_string(dataType) + ".");
        break;
    }
  }
}


/// Set the solver/coupler data for load balancing
template <MInt nDim>
void Controller<nDim>::setDataDlb(const MInt id, const MInt mode, std::vector<MInt*>& intDataRecv,
                                  std::vector<MLong*>& longDataRecv, std::vector<MFloat*>& floatDataRecv,
                                  std::vector<MInt>& dataTypes, const MBool freeMemory) {
  TRACE();

  const MBool isSolver = (mode == 0); // mode 0: solver; mode 1: coupler
  const MInt timerId = (isSolver) ? id : noSolvers() + id;
  // Pointer for easier access
  Solver* const solverP = (isSolver) ? &solver(id) : nullptr;
  Coupling* const couplerP = (isSolver) ? nullptr : &coupler(id);

  // Get the number of quantities that need to be communicated
  MInt dataCount = (isSolver) ? solverP->noCellDataDlb() : couplerP->noCellDataDlb();

  RECORD_TIMER_START(m_solverTimers[timerId][SolverTimers::SetDataMpiBlocking]);
  // Note: since inactive ranks might return the correct count, take the max among all domains
  MPI_Allreduce(MPI_IN_PLACE, &dataCount, 1, type_traits<MInt>::mpiType(), MPI_MAX, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "dataCount");
  RECORD_TIMER_STOP(m_solverTimers[timerId][SolverTimers::SetDataMpiBlocking]);

  MInt intDataCount = 0;
  MInt longDataCount = 0;
  MInt floatDataCount = 0;

  // Set solver/solver data
  for(MInt dataId = 0; dataId < dataCount; dataId++) {
    const MInt dataType = dataTypes[dataId];

    switch(dataType) {
      case MINT: {
        if(isSolver) {
          solverP->setCellDataDlb(dataId, &intDataRecv[intDataCount][0]);
        } else {
          couplerP->setCellDataDlb(dataId, &intDataRecv[intDataCount][0]);
        }
        if(freeMemory) {
          mDeallocate(intDataRecv[intDataCount]);
        }
        intDataCount++;
        break;
      }
      case MLONG: {
        if(isSolver) {
          solverP->setCellDataDlb(dataId, &longDataRecv[longDataCount][0]);
        } else {
          couplerP->setCellDataDlb(dataId, &longDataRecv[longDataCount][0]);
        }
        if(freeMemory) {
          mDeallocate(longDataRecv[longDataCount]);
        }
        longDataCount++;
        break;
      }
      case MFLOAT: {
        if(isSolver) {
          solverP->setCellDataDlb(dataId, &floatDataRecv[floatDataCount][0]);
        } else {
          couplerP->setCellDataDlb(dataId, &floatDataRecv[floatDataCount][0]);
        }
        if(freeMemory) {
          mDeallocate(floatDataRecv[floatDataCount]);
        }
        floatDataCount++;
        break;
      }
      default:
        TERMM(1, "Unknown data type.");
        break;
    }
  }
}


/// \brief Store timings of all timesteps for all domains for performance
///        evaluations (i.e. runtime and idle time).
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim>
void Controller<nDim>::storeTimings() {
  TRACE();
  using namespace maia::parallel_io;

  RECORD_TIMER_START(m_timers[Timers::IO]);

  // Number of timings
  const MInt noTimings = m_dlbTimeStepsAll.size();
  const MInt noValues = 2; // runTime and idleTime

  if(noTimings == 0) {
    RECORD_TIMER_STOP(m_timers[Timers::IO]);
    return;
  }

  // Timings output file name
  std::stringstream fileName;
  fileName << m_outputDir << "timings_" << noDomains() << "_" << globalTimeStep << ParallelIo::fileExt();

  m_log << "Write timings to file: " << fileName.str() << std::endl;

  ParallelIo file(fileName.str(), PIO_REPLACE, gridb().mpiComm());

  file.defineArray(PIO_INT, "timeStep", noTimings);

  // Dimensions of timings: noDomains x noTimings x noValues
  ParallelIo::size_type dimSizes[] = {noDomains(), noTimings, noValues};

  file.defineArray(PIO_FLOAT, "timings", 3, &dimSizes[0]);

  file.setAttribute("domain index", "dim_0", "timings");
  file.setAttribute("time step index", "dim_1", "timings");
  file.setAttribute("timings index", "dim_2", "timings");

  file.setAttribute("run time", "var_0", "timings");
  file.setAttribute("idle time", "var_1", "timings");

  // Assemble timesteps and timings
  MIntScratchSpace timeSteps(noTimings, AT_, "timeSteps");
  MFloatScratchSpace data(noTimings, noValues, AT_, "data");
  for(MInt i = 0; i < noTimings; i++) {
    timeSteps[i] = m_dlbTimeStepsAll[i];
    data(i, 0) = m_dlbRunTimeAll[i];
    data(i, 1) = m_dlbIdleTimeAll[i];
  }

  // root writes timestep data
  if(domainId() == 0) {
    file.setOffset(noTimings, 0);
  } else {
    file.setOffset(0, 0);
  }
  file.writeArray(&timeSteps[0], "timeStep");

  // Write timings of all domains
  file.setOffset(1, domainId(), 3);
  file.writeArray(&data[0], "timings");

  RECORD_TIMER_STOP(m_timers[Timers::IO]);
}


/// Store domain loads and additional information to file
template <MInt nDim>
void Controller<nDim>::storeLoadsAndWeights(const MFloat* const loads, const MInt noLoadTypes,
                                            const MInt* const loadQuantities, const MFloat domainWeight,
                                            const MFloat* const weights) {
  TRACE();
  using namespace maia::parallel_io;

  RECORD_TIMER_START(m_timers[Timers::IO]);

  if(domainId() == 0) {
    std::cerr << "Store loads and weights" << std::endl;
  }

  std::stringstream fileName;
  fileName << m_outputDir << "loads_" << noDomains() << "_" << std::setw(8) << std::setfill('0') << globalTimeStep;
  fileName << ParallelIo::fileExt();

  ParallelIo file(fileName.str(), PIO_REPLACE, gridb().mpiComm());

  file.defineArray(PIO_FLOAT, "weights", noLoadTypes);
  file.defineArray(PIO_FLOAT, "loads", noDomains());
  file.defineArray(PIO_FLOAT, "domainWeights", noDomains());

  // Dimensions of load quantities: noDomains x noLoadTypes
  ParallelIo::size_type dimSizes[] = {noDomains(), noLoadTypes};
  file.defineArray(PIO_INT, "loadQuantities", 2, &dimSizes[0]);

  // root writes estimated weights
  if(domainId() == 0) {
    file.setOffset(noLoadTypes, 0);
  } else {
    file.setOffset(0, 0);
  }
  file.writeArray(&weights[0], "weights");

  file.setOffset(1, domainId());
  file.writeArray(&loads[domainId()], "loads");
  file.writeArray(&domainWeight, "domainWeights");

  // Write load quantities of all domains
  file.setOffset(1, domainId(), 2);
  file.writeArray(&loadQuantities[0], "loadQuantities");

  RECORD_TIMER_STOP(m_timers[Timers::IO]);
}

/// Determine the global domain id of the solver local root domain
template <MInt nDim>
MInt Controller<nDim>::solverLocalRootDomain(Solver* const solver) {
  MInt localRootGlobalDomain = (solver->domainId() == 0) ? globalDomainId() : 0;
  MPI_Allreduce(MPI_IN_PLACE, &localRootGlobalDomain, 1, type_traits<MInt>::mpiType(), MPI_SUM, gridb().mpiComm(), AT_,
                "MPI_IN_PLACE", "localRootGlobalDomain");
  return localRootGlobalDomain;
}

/* brief:   Cast adaptation interval to multiple to of timestep on coarsest level
   details: Needed for mesh adaptation with LB in case of Dupuis initialization of
            newly created cells.
   author:  Philipp Brokof, Moritz Waldmann
*/
template <MInt nDim>
void Controller<nDim>::castAdaptationIntervalToMultipleOfCoarsestTimeStep(MInt maxLevel,
                                                                          MInt maxUniformRefinementLevel) {
  this->m_adaptationInterval = this->m_adaptationInterval / IPOW2(maxLevel - maxUniformRefinementLevel)
                               * IPOW2(maxLevel - maxUniformRefinementLevel);
  this->m_adaptationStart = this->m_adaptationStart / IPOW2(maxLevel - maxUniformRefinementLevel)
                            * IPOW2(maxLevel - maxUniformRefinementLevel);
  std::cout << "Set adaptationStart to: " << this->m_adaptationStart
            << " and adaptationInterval to: " << this->m_adaptationInterval << "\n";
}

template <MInt nDim>
void Controller<nDim>::resetAllTimer() {
  // Reset timer records
  maia::dlb::g_dlbTimerController.resetRecords();

  m_dlbPreviousLocalRunTime = 0.0;
  m_dlbPreviousLocalIdleTime = 0.0;

  // ... clear timings ...
  std::vector<MFloat>().swap(m_dlbTimings);

  m_dlbLastResetTimeStep = globalTimeStep;

  m_log << "Resetting DLB timers at timestep " << globalTimeStep << std::endl;
  cerr0 << "Resetting DLB timers at timestep " << globalTimeStep << std::endl;
}

template <MInt nDim>
void Controller<nDim>::logTimerStatistics(const MString& status) {
  m_log << "DLB: Timer statistics ";
  if(!status.empty()) m_log << status << " ";
  m_log << "at global time step " << globalTimeStep << std::endl;

  // Accumulate timer records of all dlb timers
  MFloat localRunTime = 0.0;
  MFloat localIdleTime = 0.0;
  const MInt noDlbTimers = maia::dlb::g_dlbTimerController.noDlbTimers();
  for(MInt i = 0; i < noDlbTimers; i++) {
    const MFloat loadRecord = maia::dlb::g_dlbTimerController.returnLoadRecord(i, m_loadBalancingTimerMode);
    const MFloat idleRecord = maia::dlb::g_dlbTimerController.returnIdleRecord(i, m_loadBalancingTimerMode);

    if(m_loadBalancingMode == 1 && !m_testDynamicLoadBalancing && (loadRecord < 0.0 || idleRecord < 0.0)) {
      TERMM(1, "Load/Idle record for dlb timer #" + std::to_string(i) + " is less than zero on global domain #"
                   + std::to_string(domainId()));
    }

    localRunTime += loadRecord;
    localIdleTime += idleRecord;
  }

  const MInt noTimeSteps = globalTimeStep - m_dlbLastResetTimeStep;

  MFloat timePerStep = (localRunTime + localIdleTime) / noTimeSteps;
  MFloat maxRunTime = localRunTime / noTimeSteps;
  const MFloat localTimePerStep = timePerStep;

  // Communicate and take maximum value -> assure same value on all domains
  MPI_Allreduce(MPI_IN_PLACE, &timePerStep, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "timePerStep");
  MPI_Allreduce(MPI_IN_PLACE, &maxRunTime, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxRunTime");

  MFloat maxIdleTime = localIdleTime / noTimeSteps;
  MFloat minIdleTime = localIdleTime / noTimeSteps;
  MPI_Allreduce(MPI_IN_PLACE, &maxIdleTime, 1, MPI_DOUBLE, MPI_MAX, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxIdleTime");
  MPI_Allreduce(MPI_IN_PLACE, &minIdleTime, 1, MPI_DOUBLE, MPI_MIN, gridb().mpiComm(), AT_, "MPI_IN_PLACE",
                "minIdleTime");


  m_log << globalTimeStep << " Average time per step: " << timePerStep << " ; local: " << localTimePerStep
        << ", idle/comp = " << localIdleTime / localRunTime
        << ", idle/timePerStep = " << localIdleTime / (noTimeSteps * timePerStep) << std::endl;
  m_log << globalTimeStep << " Relative idle time: max = " << maxIdleTime / timePerStep << ", min "
        << minIdleTime / timePerStep << std::endl;
  m_log << globalTimeStep << " maxRunTime " << maxRunTime << std::endl;
}

template <MInt nDim>
MBool Controller<nDim>::isDlbTimeStep() {
  TRACE();
  MBool dlbTimeStep = ((globalTimeStep - m_loadBalancingOffset) % m_loadBalancingInterval == 0)
                      && (m_loadBalancingStopTimeStep < 0 || globalTimeStep < m_loadBalancingStopTimeStep);
  if(m_balanceAfterAdaptation) {
    dlbTimeStep = dlbTimeStep
                  || (((globalTimeStep - m_lastAdaptationTimeStep) == m_loadBalancingOffset)
                      && (m_nAdaptationsSinceBalance >= m_balanceAdaptationInterval || m_dlbStep < 2)
                      && (m_loadBalancingStopTimeStep < 0 || globalTimeStep < m_loadBalancingStopTimeStep));

    if((m_syncTimerSteps && ((globalTimeStep - m_lastAdaptationTimeStep) < m_loadBalancingOffset)) || m_syncTimeSteps) {
      cerr0 << "Applying Barrier for correct timer!" << std::endl;
      if(m_syncTimeSteps && !((globalTimeStep - m_lastAdaptationTimeStep) < m_loadBalancingOffset)) {
        solver(0).startLoadTimer(AT_);
      }
      MPI_Barrier(gridb().mpiComm(), AT_);
      if(m_syncTimeSteps && !((globalTimeStep - m_lastAdaptationTimeStep) < m_loadBalancingOffset)) {
        solver(0).stopLoadTimer(AT_);
      }
    }
  }

  return dlbTimeStep;
}
template <MInt nDim>
MBool Controller<nDim>::isAdaptationTimeStep() {
  TRACE();

  return (m_adaptationInterval > 0 && globalTimeStep % m_adaptationInterval == 0 && globalTimeStep > m_adaptationStart);
}


} // namespace grid
} // namespace maia


#endif // ifndef GRIDCONTROLLER_H_
