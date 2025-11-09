// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SOLVER_H
#define SOLVER_H

#include <bitset>
#include <limits>
#include <set>
#include "GEOM/geometry.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/dlbtimer.h"
#include "UTIL/timer.h"
#include "application.h"

class Coupling;
class Application;

/** \brief Parent class of all solvers
 * This class is the base for all solvers.
 * I.e. all solver class (e.g. structured or cartesian)
 * are derived from this class.
 */
class Solver {
  // Declared friend to allow for functions which should only be used in Application
  friend Application;
  template <MInt nDim, class ppType>
  friend class PostProcessing;

 private:
  // The path to the testcase directory
  MString m_testcaseDir;

  // The path to the output directory
  MString m_outputDir;

  // The path to the restart directory
  MString m_restartDir;

  // The MPI communicator to be used by this solver
  MPI_Comm m_mpiComm;

  // The domain id of this solver on the current MPI communicator
  MInt m_domainId;

  // The number of domains of this solver on the current MPI communicator
  MInt m_noDomains;

  // The current solver status in the execution recipe (solver enabled/disabled in current step)
  MBool m_solverStatus = true;

  MString m_solverMethod;
  MString m_solverType;

  const MInt m_noDim;

  MInt m_dlbTimerId = -1;

  static std::map<MInt, MString> m_aliases;

 protected:
  /// the Reynolds number
  MFloat m_Re{};

  /// the Mach number
  MFloat m_Ma{};

  /// The number of timesteps before writing the next solution file
  MInt m_solutionInterval;
  MInt m_solutionOffset{};
  std::set<MInt> m_solutionTimeSteps;

  /// The number of timesteps before writing the next restart file
  MInt m_restartInterval;
  MInt m_restartTimeStep;
  MInt m_restartOffset;
  MString m_solutionOutput;
  MBool m_useNonSpecifiedRestartFile = false;
  MBool m_initFromRestartFile;

  /// The number of timesteps before writing the next residual
  MInt m_residualInterval;

  /// a unique solver identifier
  const MInt m_solverId;

  MFloat* m_outerBandWidth = nullptr;
  MFloat* m_innerBandWidth = nullptr;
  MInt* m_bandWidth = nullptr;

  // TODO: remove one, they are identically!
  MBool m_restart = false;
  MBool m_restartFile = false;

 public:
  std::set<MInt> m_freeIndices;

  MBool m_singleAdaptation = false;

  // TODO labels:toremove delete once all solvers and applications have splitAdaptation
  MBool m_splitAdaptation = true; // TODO labels:toremove remove

  // Save sensor data to file
  MBool m_saveSensorData = false;

  // protected doesnot allow instances of non-derived classes
 protected:
  Solver(const MInt solverId, const MPI_Comm comm, const MBool isActive = true);

 public:
  MString getIdentifier(const MBool useSolverId = false, const MString preString = "", const MString postString = "_");

  virtual ~Solver() = default;

  /// Return the number of internal cells within this solver
  virtual MInt noInternalCells() const = 0;

  /// Return the time
  virtual MFloat time() const = 0;

  /// Return the number of variables
  virtual MInt noVariables() const { TERMM(-1, "Not implemented for this solver"); }

  /// Return the dimensionalization parameters of this solver
  virtual void getDimensionalizationParams(std::vector<std::pair<MFloat, MString>>& /*dimParams*/) const {
    TERMM(-1, "getDimensionalizationParams()");
  };

  /// Set new domain information
  void updateDomainInfo(const MInt domainId, const MInt noDomains, const MPI_Comm mpiComm, const MString& loc) {
    m_domainId = domainId;
    m_noDomains = noDomains;
    m_mpiComm = mpiComm;

    // DEBUG output
    m_log << "Solver #" << solverId() << " updateDomainInfo: domainId=" << domainId << ", noDomains=" << noDomains
          << " from " << loc << std::endl;
    //#ifndef NDEBUG
    //    std::cerr << "Solver #" << solverId() << " updateDomainInfo: domainId=" << domainId
    //              << ", noDomains=" << noDomains << " from " << loc << std::endl;
    //#endif
  }

  // relevant for Postprocessing
  // virtual MInt& a_level(const MInt) { TERMM(1, "Not implemented for this solver"); };
  virtual MFloat& a_slope(const MInt, MInt const, const MInt) { TERMM(1, "Not implemented for this solver"); };
  virtual MBool a_isBndryCell(const MInt) const { TERMM(1, "Not implemented for this solver"); };
  virtual MFloat& a_FcellVolume(MInt) { TERMM(1, "Not implemented for this solver"); };
  virtual MInt getCurrentTimeStep() const {
    TERMM(1, "Not implemented for this solver");
    return -1;
  }
  virtual void accessSampleVariables(MInt, MFloat*&) { TERMM(1, "Not implemented for this solver"); };
  virtual void getSampleVariableNames(std::vector<MString>& NotUsed(varNames)) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual MBool a_isBndryGhostCell(MInt /*cellId*/) const { return false; }
  // only implemented in LBSolver
  virtual void saveCoarseSolution() { TERMM(1, "Not implemented for this solver"); }; // only implemented in LBSolver

  virtual void getSolverSamplingProperties(std::vector<MInt>& NotUsed(samplingVarIds),
                                           std::vector<MInt>& NotUsed(noSamplingVars),
                                           std::vector<std::vector<MString>>& NotUsed(samplingVarNames),
                                           const MString NotUsed(featureName) = "") {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void initSolverSamplingVariables(const std::vector<MInt>& NotUsed(varIds),
                                           const std::vector<MInt>& NotUsed(noSamplingVars)) {
    TERMM(1, "Not implemented for this solver");
  };

  virtual void calcSamplingVariables(const std::vector<MInt>& NotUsed(varIds), const MBool NotUsed(exchange)) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void calcSamplingVarAtPoint(const MFloat* NotUsed(point), const MInt NotUsed(id),
                                      const MInt NotUsed(sampleVarId), MFloat* NotUsed(state),
                                      const MBool NotUsed(interpolate) = false) {
    TERMM(1, "Not implemented for this solver");
  };


  /// Perform load balancing
  virtual void balance(const MInt* const NotUsed(noCellsToReceiveByDomain),
                       const MInt* const NotUsed(noCellsToSendByDomain),
                       const MInt* const NotUsed(targetDomainsByCell),
                       const MInt NotUsed(oldNoCells)) {
    TERMM(1, "Not implemented for this solver");
  }

  /// Return if load balancing for solver is split into multiple methods or implemented in balance()
  virtual MBool hasSplitBalancing() const { return false; }
  virtual void balancePre() { TERMM(1, "Not implemented for this solver"); };
  virtual void balancePost() { TERMM(1, "Not implemented for this solver"); };

  virtual void finalizeBalance() { TERMM(1, "Not implemented for this solver"); };

  /// Reset the solver/solver for load balancing
  virtual void resetSolver() { TERMM(1, "Not implemented for this solver " + std::to_string(m_solverId)); };

  /// Cancel open mpi (receive) requests in the solver (e.g. due to interleaved execution)
  virtual void cancelMpiRequests(){};

  /// Set cell weights
  virtual void setCellWeights(MFloat* /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  virtual MInt noLoadTypes() const { return -1; };
  virtual void getDefaultWeights(MFloat* NotUsed(weights), std::vector<MString>& NotUsed(names)) const {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void getLoadQuantities(MInt* const NotUsed(loadQuantities)) const {};
  virtual MFloat getCellLoad(const MInt NotUsed(cellId), const MFloat* const NotUsed(weights)) const { return -1.0; };
  virtual void limitWeights(MFloat* NotUsed(weights)){};

  virtual void localToGlobalIds(/*const MInt* const NotUsed(offsets)*/){
      // todo labels:toenhance fix this by writing a stub in offending solver with proper warning message or reason why
      // impl is not
      // needed
      // std::cerr << "ERROR: calling localToGlobalIds() within solver without implementation!"<<std::endl;
      // TERMM(1, "Not implemented for this solver");
  };
  virtual void globalToLocalIds(/*const MInt* const NotUsed(offsets)*/){
      // todo labels:toenhance fix this by writing a stub in offending solver with proper warning message or reason why
      // impl is not
      // needed
      // std::cerr << "ERROR: calling globalToLocalIds() within solver without implementation!"<<std::endl;
      //    TERMM(1, "Not implemented for this solver");
  };

  /// Methods to inquire solver data information
  virtual MInt noCellDataDlb() const { return 0; };
  virtual MInt cellDataTypeDlb(const MInt NotUsed(dataId)) const { return -1; };
  virtual MInt cellDataSizeDlb(const MInt NotUsed(dataId), const MInt NotUsed(cellId)) { return -1; };

  // Note: create function for every datatype, else overloaded function with template parameter is
  // not called in derived solver class when using a Solver pointer.
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MInt* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MLong* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MInt* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MLong* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented for solver.");
  }

  // Methods to get/set global solver variables that should be the same for all ranks
  virtual void getGlobalSolverVars(std::vector<MFloat>& NotUsed(globalFloatVars),
                                   std::vector<MInt>& NotUsed(globalIntVars)){
      // todo labels:toenhance fix this by writing a stub in offending solver with proper warning message or reason why
      // impl is not
      // needed
      // std::cerr << "ERROR: calling getGlobalSolverVars() within solver without implementation!"<<std::endl;
      //    TERMM(1, "Not implemented for this solver");
  };
  virtual void setGlobalSolverVars(std::vector<MFloat>& NotUsed(globalFloatVars),
                                   std::vector<MInt>& NotUsed(globalIdVars)){
      // todo labels:toenhance fix this by writing a stub in offending solver with proper warning message or reason why
      // impl is not
      // needed
      // std::cerr << "ERROR: calling setGlobalSolverVars() within solver without implementation!"<<std::endl;
      //    TERMM(1, "Not implemented for this solver");
  };

  // Enable timers for dynamic load balancing, disabled by default since not supported by all
  // solvers/run loops
  inline void enableDlbTimers() { maia::dlb::g_dlbTimerController.enableDlbTimers(m_dlbTimerId); }

  // re-enable timers, and keep them running
  inline void reEnableDlbTimers() { maia::dlb::g_dlbTimerController.reEnableDlbTimer(m_dlbTimerId); }

  // Temporarily disable timers, e.g. during adaptation
  inline void disableDlbTimers() { maia::dlb::g_dlbTimerController.disableDlbTimers(m_dlbTimerId); }

  inline MBool dlbTimersEnabled() { return maia::dlb::g_dlbTimerController.dlbTimersEnabled(m_dlbTimerId); }

  inline void startLoadTimer(const MString name) { maia::dlb::g_dlbTimerController.startLoadTimer(m_dlbTimerId, name); }

  inline void stopLoadTimer(const MString& name) { maia::dlb::g_dlbTimerController.stopLoadTimer(m_dlbTimerId, name); }

  inline void stopIdleTimer(const MString& name) { maia::dlb::g_dlbTimerController.stopIdleTimer(m_dlbTimerId, name); }

  inline void startIdleTimer(const MString& name) {
    maia::dlb::g_dlbTimerController.startIdleTimer(m_dlbTimerId, name);
  }

  inline MBool isLoadTimerRunning() { return maia::dlb::g_dlbTimerController.isLoadTimerRunning(m_dlbTimerId); };

  virtual MInt noSolverTimers(const MBool NotUsed(allTimings)) { return maia::dlb::g_dlbTimerController.noSubTimers(); }

  virtual void getSolverTimings(std::vector<std::pair<MString, MFloat>>& NotUsed(solverTimings),
                                const MBool NotUsed(allTimings)) {
    TERMM(1, "not implemented in base class");
  }

  virtual void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& NotUsed(domainInfo)) {
    TERMM(1, "not implemented for solver");
  };

  void setDlbTimer(const MInt timerId) {
    if(m_dlbTimerId != -1) {
      TERMM(1, "m_dlbTimerId already set");
    }
    m_dlbTimerId = timerId;
  }

  /// Prepare the solver/solver for adaptation and collect refinement sensors
  /// TODO labels:toremove remove once all solvers use the split adaptation
  virtual void prepareAdaptation(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<MFloat>& /*unused*/,
                                 std::vector<std::bitset<64>>& /*unused*/, std::vector<MInt>& /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };

  /// Reinit the solver/solver after adaptation
  /// TODO labels:toremove remove once all solvers use the split adaptation
  virtual void reinitAfterAdaptation() { TERMM(1, "Not implemented for this solver"); };

  /// prepare adaptation for split adaptation before the adaptation loop
  virtual void prepareAdaptation() { TERMM(1, "Not implemented for this solver"); };

  /// set solver sensors for split adaptation within the adaptation loop
  virtual void setSensors(std::vector<std::vector<MFloat>>&, std::vector<MFloat>&, std::vector<std::bitset<64>>&,
                          std::vector<MInt>&) {
    TERMM(1, "Not implemented for this solver");
  };

  virtual void saveSensorData(const std::vector<std::vector<MFloat>>& /*sensors*/, const MInt& /*level*/,
                              const MString& /*gridFileName*/, const MInt* const /*recalcIds*/){};

  /// post adaptation for split adaptation within the adaptation loop
  virtual void postAdaptation() { TERMM(1, "Not implemented for this solver"); };

  /// finalize adaptation for split sadptation after the adaptation loop
  virtual void finalizeAdaptation() { TERMM(1, "Not implemented for this solver"); };

  /// Refine the given cell
  virtual void refineCell(const MInt /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Coarsen the given cell
  virtual void removeChilds(const MInt /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Remove the given cell
  virtual void removeCell(const MInt /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Swap the given cells
  virtual void swapCells(const MInt /*unused*/, const MInt /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Swap the given cells
  virtual void swapProxy(const MInt /*unused*/, const MInt /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Check whether cell is outside the fluid domain
  virtual MInt cellOutside(const MFloat* /*unused*/, const MInt /*unused*/, const MInt /*unused*/) { return -1; };

  /// Swap the given cells
  virtual void resizeGridMap() { TERMM(1, "Not implemented for this solver"); };

  /// Prepare the solvers for a grid-restart
  virtual MBool prepareRestart(MBool /*unused*/, MBool& /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  // Reinit the solvers after a restart
  virtual void reIntAfterRestart(MBool /*unused*/) { TERMM(1, "Not implemented for this solver"); };

  /// Return the MPI communicator used by this solver
  MPI_Comm mpiComm() const { return m_mpiComm; }

  /// Return the domainId (rank)
  virtual MInt domainId() const { return m_domainId; }

  /// Return the total number of domains (total number of ranks in current MPI
  /// communicator)
  virtual MInt noDomains() const { return m_noDomains; }

  /// Return if the solver is active on this rank; needs to be implemented in derived solver since
  /// access to gridproxy not possible here
  virtual MBool isActive() const {
    TERMM(1, "Not implemented for solver!");
    return false;
  }

  /// Set the solver status to indicate if the solver is currently active in the execution recipe
  /// Note: might be required for couplers to check which solvers currently require coupling
  void setSolverStatus(const MBool status) { m_solverStatus = status; }

  /// Get the solver status indicating if the solver is currently active in the execution recipe
  MBool getSolverStatus() { return m_solverStatus; }

  /// Return the testcase directory
  MString testcaseDir() const { return m_testcaseDir; }

  /// Return the directory for output files
  MString outputDir() const { return m_outputDir; }

  /// Return the directory for restart files
  MString restartDir() const { return m_restartDir; }

  /// Return the solverMethod of this solver
  MString solverMethod() const { return m_solverMethod; }

  /// Return the solverType of this solver
  MString solverType() const { return m_solverType; }

  /// Return the restart interval of this solver
  MInt restartInterval() const { return m_restartInterval; }

  /// Return the restart interval of this solver
  MInt restartTimeStep() const { return m_restartTimeStep; }

  /// Return the solverId
  MInt solverId() const { return m_solverId; };

  MBool restartFile() { return m_restartFile; };

  MInt readSolverSamplingVarNames(std::vector<MString>& varNames, const MString featureName = "") const;

  /// Return if the restart time step can be determined from the restart file (for
  /// useNonSpecifiedRestartFile = true) see determineRestartTimeStep()
  virtual MBool hasRestartTimeStep() const { return false; }

  virtual MBool forceAdaptation() { return false; };
  // If false is commented in, compilation with PGI will fail!
  virtual void preTimeStep() = 0;

  virtual void postTimeStep() = 0;

  virtual void initSolver() = 0;

  virtual void finalizeInitSolver() = 0;

  virtual void saveSolverSolution(const MBool NotUsed(forceOutput) = false,
                                  const MBool NotUsed(finalTimeStep) = false) = 0;

  virtual void cleanUp() = 0;

  // If false is commented in, compilation with PGI will fail!
  virtual MBool solutionStep() { TERMM(-1, "Empty body of Solver::solutionStep()"); };
  virtual void preSolutionStep(MInt /*mode*/) { TERMM(-1, "Empty body of Solver::preSolutionStep()"); };

  virtual MBool postSolutionStep() { TERMM(-1, "Empty body of Solver::postSolutionStep()"); };


  virtual MBool solverConverged() { return false; };

  virtual void getInterpolatedVariables(MInt /*cellId*/, const MFloat* /*position*/, MFloat* /*vars*/) {
    TERMM(-1, "Empty body of Solver::getInterpolatedVariables()");
  };

  virtual void loadRestartFile() { TERMM(-1, "Empty body of Solver::loadRestartFile()"); };

  virtual MInt determineRestartTimeStep() const {
    TERMM(-1, "Empty body of Solver::determineRestartTimeStep()");
    return -1;
  };

  virtual void writeRestartFile(MBool /*unused*/) { TERMM(-1, "Empty body of Solver::writeRestartFile()"); };

  virtual void writeRestartFile(const MBool /*unused*/, const MBool /*unused*/, const MString /*unused*/,
                                MInt* /*unused*/) {
    TERMM(-1, "Empty body of Solver::writeRestartFile()");
  };

  virtual void setTimeStep() { TERMM(-1, "Empty body of Solver::setTimeStep()"); };

  virtual void implicitTimeStep() { TERMM(-1, "Empty body of Solver::implicitTimeStep()"); };

  virtual void prepareNextTimeStep() { return; };

 protected:
  // Methods for solver timers needed for dynamic load balancing
  MFloat returnLoadRecord() const { return maia::dlb::g_dlbTimerController.returnLoadRecord(m_dlbTimerId); }
  MFloat returnIdleRecord() const { return maia::dlb::g_dlbTimerController.returnIdleRecord(m_dlbTimerId); }

 private:
  void initAdaptation();
};

#endif // SOLVER_H
