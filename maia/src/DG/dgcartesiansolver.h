// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSOLVER_H_
#define DGSOLVER_H_

#include <array>
#include <memory>
#include <vector>
#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangrid.h"
#include "POST/postprocessing.h"
#include "POST/samplingdata.h"
#include "cartesiansolver.h"
#include "dgcartesianboundarycondition.h"
#include "dgcartesianelementcollector.h"
#include "dgcartesiangridmap.h"
#include "dgcartesianhelementcollector.h"
#include "dgcartesianinterpolation.h"
#include "dgcartesianslices.h"
#include "dgcartesiansponge.h"
#include "dgcartesiansurfacecollector.h"
#include "dgcartesiantimers.h"
#include "filter.h"
#include "sbpcartesianinterpolation.h"

template <MInt nDim, class SysEqn>
class DgBoundaryCondition;
template <MInt nDim, class SysEqn>
class DgBoundaryConditionFactory;
template <MInt nDim, class SysEqn>
class CouplingDg;
template <MInt nDim, class CouplingX>
class CouplingDgApe;
template <MInt nDim, class SysEqn>
class DgCcAcousticPerturb;
template <MInt nDim, class ppType>
class PostProcessing;

template <MInt nDim, class SysEqn>
class DgCartesianSolver : public maia::CartesianSolver<nDim, DgCartesianSolver<nDim, SysEqn>> {
  friend class DgBoundaryCondition<nDim, SysEqn>;
  template <MInt nDim_, class CouplingX>
  friend class CouplingDgApe;
  friend class DgCcAcousticPerturb<nDim, SysEqn>;
  friend class CouplingDg<nDim, SysEqn>;
  friend class DgSlices<nDim, SysEqn>;
  friend class maia::CartesianSolver<nDim, DgCartesianSolver<nDim, SysEqn>>;
  template <MInt nDim_, class ppType>
  friend class PostProcessing;

  // Types
 private:
  using ElementCollector = maia::dg::collector::ElementCollector<nDim, SysEqn>;
  using HElementCollector = maia::dg::collector::HElementCollector<nDim, SysEqn>;
  using SurfaceCollector = maia::dg::collector::SurfaceCollector<nDim, SysEqn>;
  using Timers = maia::dg::Timers_;
  using BC = typename DgBoundaryConditionFactory<nDim, SysEqn>::ReturnType;

 public:
  using Geom = Geometry<nDim>;
  using CartesianSolver = typename maia::CartesianSolver<nDim, DgCartesianSolver>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;
  using CartesianSolver::disableDlbTimers;
  using CartesianSolver::domainId;
  using CartesianSolver::enableDlbTimers;
  using CartesianSolver::getIdentifier;
  using CartesianSolver::grid;
  using CartesianSolver::isActive;
  using CartesianSolver::m_restart;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_restartTimeStep;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solverId;
  using CartesianSolver::m_useNonSpecifiedRestartFile;
  using CartesianSolver::mpiComm;
  using CartesianSolver::noDomains;
  using CartesianSolver::outputDir;
  using CartesianSolver::readSolverSamplingVarNames;
  using CartesianSolver::restartDir;
  using CartesianSolver::returnIdleRecord;
  using CartesianSolver::returnLoadRecord;
  using CartesianSolver::solverId;
  using CartesianSolver::startLoadTimer;
  using CartesianSolver::stopLoadTimer;
  using CartesianSolver::updateDomainInfo;

  // Type for cell properties
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  // Methods
 public:
  DgCartesianSolver(const MInt solverId, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm);
  ~DgCartesianSolver() override;
  void run();
  MInt getCurrentTimeStep() const override { return m_timeStep; }
  /// Force time step externally
  void forceTimeStep(const MFloat dt) {
    m_externalDt = dt;
    if(m_restart) {
      // m_dt is only updated via calcTimeStep(), not necessarily at a restart
      m_dt = m_externalDt;
    }
  }

  // STUBS required after method removal
  MFloat time() const override { return m_time; }
  MInt noVariables() const override { return SysEqn::noVars(); }

  void loadSampleVariables(MInt timeStep) {
    std::cerr << "loadSampleVariables DgCartesianSolver " << timeStep << std::endl;
  };
  void getSampleVariables(MInt cellId, const MFloat*& vars) {
    std::cerr << "getSampleVariables DgCartesianSolver " << cellId << " " << vars << std::endl;
  };
  void getSampleVariables(const MInt cellId, std::vector<MFloat>& /*vars*/) {
    std::cerr << "getSampleVariables DgCartesianSolver " << cellId << std::endl;
  };
  void calculateHeatRelease() { std::cerr << "calculateHeatRelease DgCartesianSolver " << std::endl; }
  void getHeatRelease(MFloat*& heatRelease) {
    std::cerr << "getHeatRelease DgCartesianSolver " << heatRelease << std::endl;
  }

  /// Access the solver's geometry
  constexpr const Geom& geometry() const { return m_geometry; }

 private:
  Geom& m_geometry; ///< Reference to geometry object

  /// Access the solver's geometry (non-const version)
  Geom& geometry() { return m_geometry; }

 public:
  // These methods are necessary for the postprocessing solver to work
  const MFloat& a_coordinate(const MInt cellId, const MInt dim) const { return grid().tree().coordinate(cellId, dim); }
  MInt a_noCells() const { TERMM(1, "Make this function return something meaningful in future!"); }
  MInt a_hasNeighbor(const MInt cellId, const MInt dir) const { return grid().tree().hasNeighbor(cellId, dir); }
  MInt a_level(const MInt) const { TERMM(1, "Make this function return something meaningful in future!"); }

 private:
  // Constructor methods
  void initDgCartesianSolver();
  void initTimers();
  void setNumericalProperties();
  void setInputOutputProperties();
  void allocateAndInitSolverMemory();

  // Initialization methods for the solver
 public:
  void initSolver() override;
  void finalizeInitSolver() override;

 private:
  void initGridMap();
  void loadGridMap(const MString& gridMapFileName);
  void checkGridMap(const MString& donorGridFileName);
  void initElements();
  MBool needElementForCell(const MInt cellId);
  MInt calculateNeededNoSurfaces();
  MBool pointIsInside(const MFloat* const coordinates);
  void createElement(const MInt cellId);
  void initPolynomialDegree();
  void initInterpolation();
  void initNodeCoordinates();
  void initJacobian();
  void checkCellProperties();
  std::array<MInt, 2 * nDim> getBoundaryConditionIds(const MInt cellId);
  void initSurfaces();
  MInt initBoundarySurface(const MInt elementId, const MInt dir);
  MInt initInnerSurface(const MInt elementId, const MInt dir);
  MInt initMpiSurface(const MInt elementId, const MInt dir);
  MBool hasSurface(const MInt elementId, const MInt dir);
  MInt createSurface(const MInt elementId, const MInt dir, MInt nghbrId = -1);
  void calcSurfaceNodeCoordinates(const MInt surfaceId);
  void calcElementNodeCoordinates(const MInt elementId);
  void initMpiExchange();
  void updateNodeVariables();
  void extendNodeVariables();
  void extendNodeVariablesSingleDirection(const MInt extendDir, const MFloat extendOffset, const MFloat extendLimit);
  void updateNodeVariablesSingleElement(const MInt elementId,
                                        const MInt extendDir,
                                        MIntTensor hForwardSurfaceId,
                                        MIntTensor hReverseSurfaceId,
                                        std::vector<MBool>& cellHasUpdatedMeans,
                                        std::vector<MInt>& SurfaceWantsToMPISend,
                                        MBool& noMoreUpdates);
  void initSimulationStatistics();
  void initSolverObjects();
  void initData();
  void initMainLoop();

  // Methods for handling elements and cells
 private:
  MInt getLevelByElementId(const MInt elementId);
  //  MInt getLevelOfNeighborElement(const MInt elementId, const MInt dir);
  MInt getLevelOfDirectNeighborCell(const MInt cellId, const MInt dir);
  MBool hasNeighborCell(const MInt elementId, const MInt dir);
  //  MBool hasNeighborElement(const MInt elementId, const MInt dir);
  MBool isMpiSurface(const MInt id) const;
  void updateWeightsAndWorkloads();
  void savePartitionCellWorkloads();
  void writeEOC();
  MInt getHElementId(const MInt elementId) const;

  // Finalization methods for the solver
 public:
  void cleanUp() override;

 private:
  void finalizeMpiExchange();

  // Initialization methods for the solution
  virtual void initialCondition();
  void outputInitSummary();

  // Input/output
  void saveSolutionFile();
  void saveSolutionFile(const MString& suffix);
  virtual void saveRestartFile();
  void loadRestartFile() override;
  void saveNodalData(const MString& fileNameBase,
                     const MInt noVars,
                     const std::vector<MString>& varNames,
                     const MFloat* const data) const;
  void saveNodeVarsFile();
  void loadNodeVarsFile();
  MString getRestartFileName(const MInt timeStep, const MInt useNonSpecifiedRestartFile);

  // Physical properties and system of equations
  // Method to create new boundary condition
  BC make_bc(MInt bcId);

  // Main DG discretization methods
  void calcDgTimeDerivative(const MFloat t, const MFloat tStage);
  void prolongToSurfaces();
  void prolongToSurfaces(const MInt begin, const MInt end);
  void applyForwardProjection();
  void calcVolumeIntegral();
  template <class F>
  void calcVolumeIntegral(const MInt noElements, ElementCollector& elem, F& fluxFct);
  void calcBoundarySurfaceFlux(MFloat t);
  void calcInnerSurfaceFlux();
  void calcMpiSurfaceFlux();
  template <class F>
  void calcRegularSurfaceFlux(const MInt begin, const MInt end, SurfaceCollector& surf, F& riemannFct);
  void calcSurfaceIntegral();
  void calcSurfaceIntegral(const MInt begin, const MInt end, ElementCollector& elem, SurfaceCollector& surf,
                           HElementCollector& helem, const MInt noHElements);
  void applySurfaceIntegral(MFloat* rhs, const MInt polyDeg, const MInt noNodes1D, const MInt srfcId, const MInt side,
                            const MFloat* flux, SurfaceCollector& surf);
  void applyJacobian();
  void applyJacobian(const MInt noElements, ElementCollector& elem);
  void calcSourceTerms(MFloat t);
  template <class F>
  void calcSourceTerms(MFloat t, const MInt noElements, ElementCollector& elem, F& sourceFct);
  void applyExternalSourceTerms(const MFloat time);

  // Parallelization
  MBool isMpiRoot() const;
  MBool hasMpiExchange() const;
  void startMpiSurfaceExchange();
  void finishMpiSurfaceExchange();
  void cancelMpiRequests() override;

  // Time integration
 public:
  // methods for NEW unified run loop
  void preTimeStep() override;
  MBool solutionStep() override;
  void postTimeStep();
  // DG-Solver is converged if the final time is reached
  MBool solverConverged() { return m_finalTimeStep; };

  void saveSolverSolution(const MBool forceOutput, const MBool finalTimeStep) override;
  void writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString gridFileName,
                        MInt* recalcIdTree) override;
  MBool prepareRestart(MBool, MBool&) override;
  virtual void reIntAfterRestart(MBool){};

  // TODO labels:DG dummy functions that are required to allow adaptation for another solver in a multisolver
  // computation without changing the grid for the DG solver
  void resizeGridMap() {
    // Note: resize with current tree size since the number of cells should not change
    grid().resizeGridMap(grid().tree().size());
  }
  void swapProxy(const MInt cellId0, const MInt cellId1) { grid().swapGridIds(cellId0, cellId1); }

  void prepareAdaptation() override{};
  void setSensors(std::vector<std::vector<MFloat>>& NotUsed(sensors),
                  std::vector<MFloat>& NotUsed(sensorWeight),
                  std::vector<std::bitset<64>>& NotUsed(sensorCellFlag),
                  std::vector<MInt>& NotUsed(sensorSolverId)) override{};
  void postAdaptation() override{};
  void finalizeAdaptation() override {
    // TODO labels:DG,GRID check this, previously in reinitAfterAdaptation()
    grid().updateOther();
    updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);
  };


  MBool step(const MFloat externalDt = -std::numeric_limits<MFloat>::infinity(), const MBool substep = false);

  MFloat calcTimeStep();

  /// * Functions for postprocessing/sampling *

  /// Return cell id belonging to an element id/index
  MInt getCellIdByIndex(const MInt index) { return m_elements.cellId(index); }

  /// Return the element id containing a given point
  MInt getIdAtPoint(const MFloat* point, const MBool globalUnique = false) {
    return getElementIdAtPoint(point, globalUnique);
  }

  void calcSamplingVarAtPoint(const MFloat* point, const MInt elementId, const MInt sampleVarId, MFloat* state,
                              const MBool interpolate) override;
  void getSolverSamplingProperties(std::vector<MInt>& samplingVarIds, std::vector<MInt>& noSamplingVars,
                                   std::vector<std::vector<MString>>& samplingVarNames,
                                   const MString featureName = "") override;

  // Sampling functions which have nothing to do for the DG solver at the moment
  virtual void initSolverSamplingVariables(const std::vector<MInt>& NotUsed(varIds),
                                           const std::vector<MInt>& NotUsed(noSamplingVars)){};
  virtual void initInterpolationForCell(const MInt NotUsed(cellId)){};
  virtual void calcSamplingVariables(const std::vector<MInt>& NotUsed(varIds), const MBool NotUsed(exchange)){};

 private:
  void timeStepRk(const MFloat t, const MFloat dt, const MInt substep = -1);
  void subTimeStepRk(const MFloat dt, const MInt stage, const MInt totalSize, const MFloat* const rhs,
                     MFloat* const variables, MFloat* const timeIntStorage);
  void analyzeTimeStep(MFloat t, MFloat runTimeRelative, MFloat runTimeTotal, MInt timeStep, MFloat dt);
  void calcErrorNorms(const MFloat t, std::vector<MFloat>& L2Error, std::vector<MFloat>& LInfError,
                      std::vector<MFloat>& L2ErrLocal, std::vector<MFloat>& LInfErrLocal);
  virtual void resetRHS();
  void resetExternalSources();
  void resetBuffer(const MInt totalSize, MFloat* const buffer);

  // Post-Processing
  //  MInt getCellIdAtPoint(const MFloat* point);
  MInt getElementIdAtPoint(const MFloat* point, MBool globalUnique = false);
  MBool calcStateAtPoint(const MFloat* point, MFloat* state);
  void calcStateAtPoint(const MFloat* point, const MInt elementId, MFloat* state);
  void calcStateAtPoint(const MFloat* point, const MInt elementId, const MInt noVars, const MFloat* const u,
                        MFloat* state);
  //  MBool calcNodeVarsAtPoint(const MFloat* point, MFloat* nodeVars);
  //  void
  //  calcNodeVarsAtPoint(const MFloat* const point, const MInt elementId,
  //  MFloat* const nodeVars);


 public:
  MInt noInternalCells() const override { return grid().noInternalCells(); }

 private:
  // Sponge
  DgSponge<nDim, SysEqn>& sponge() { return m_sponge; }
  MBool useSponge() const { return m_useSponge; }

  // H/p-refinement
  void initMortarProjections();
  void initPrefinement();
  template <MBool forward, MInt noVars>
  void calcMortarProjectionH(const MInt srfcId, const MInt dir, MFloat* source, MFloat* destination,
                             ElementCollector& elem, SurfaceCollector& surf);
  template <MBool forward, MInt noVars>
  void calcMortarProjectionP(const MInt srfcId, const MInt dir, MFloat* source, MFloat* destination,
                             ElementCollector& elem, SurfaceCollector& surf);
  template <MBool forward, MInt noVars>
  void calcMortarProjection(const MInt srfcId, const MInt dir, MFloat* source, MFloat* destination,
                            ElementCollector& elem, SurfaceCollector& surf);
  template <MBool forward>
  MFloat* mortarP(const MInt sourcePolyDeg, const MInt targetPolyDeg);
  template <MBool forward>
  MFloat* mortarH(const MInt polyDeg, const MInt position);
  void exchangeMpiSurfacePolyDeg();
  void adaptiveRefinement(const MInt timeStep);
  void calcErrorEstimate(std::vector<MFloat>& errorEstimate);
  void interpolateElement(const MInt elementId, const MInt newPolyDeg);
  MBool hasAdaptivePref() const;
  MBool needHElementForCell(const MInt cellId);
  void createHElement(const MInt cellId);
  void initHElements();
  MInt createHMPISurfaces(const MInt elementId, const MInt dir);
  MBool hasPref() const;
  MBool isAdaptationTimeStep(const MInt timeStep) const;


 public:
  // Dynamic load balancing
  void resetSolver() override;
  void setCellWeights(MFloat* solverCellWeight) override;

  // Partitioning
  MInt noLoadTypes() const override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;

  void balance(const MInt* const NotUsed(noCellsToReceiveByDomain),
               const MInt* const NotUsed(noCellsToSendByDomain),
               const MInt* const NotUsed(targetDomainsByCell),
               const MInt NotUsed(oldNoCells)) override {
    TERMM(1, "Use split balancing methods for DG solver instead of balance().");
  };
  // DG solver requires the balancing to be split into separate methods
  MBool hasSplitBalancing() const override { return 1; }
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override{};

  void localToGlobalIds() override;
  void globalToLocalIds() override;

  /// Methods to inquire solver data information
  MInt noCellDataDlb() const override {
    if(grid().isActive()) {
      return noDgCartesianSolverCellData() + noBcSolverCellData();
    } else {
      return 0;
    }
  };
  MInt cellDataTypeDlb(const MInt dataId) const override;
  MInt cellDataSizeDlb(const MInt dataId, const MInt cellId) override;
  template <typename dataType>
  void getCellDataDlb_(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                       dataType* const data);
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MInt* const data) override {
    getCellDataDlb_(dataId, oldNoCells, bufferIdToCellId, data);
  };
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override {
    getCellDataDlb_(dataId, oldNoCells, bufferIdToCellId, data);
  };
  template <typename dataType>
  void setCellDataDlb_(const MInt dataId, const dataType* const data);
  void setCellDataDlb(const MInt dataId, const MInt* const data) override { setCellDataDlb_(dataId, data); };
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override { setCellDataDlb_(dataId, data); };

  void getGlobalSolverVars(std::vector<MFloat>& NotUsed(globalFloatVars),
                           std::vector<MInt>& NotUsed(globalIntVars)) override;
  void setGlobalSolverVars(std::vector<MFloat>& NotUsed(globalFloatVars),
                           std::vector<MInt>& NotUsed(globalIdVars)) override;

  MInt noSolverTimers(const MBool allTimings) override {
#ifdef MAIA_TIMER_FUNCTION
    if(allTimings) {
      return 28;
    } else {
      return 5;
    }
#else
    return 2;
#endif
  }

  void getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings, const MBool allTimings) override;
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override;

 private:
  MBool m_weightDgRbcElements = false;

  MInt noDgCartesianSolverCellData() const { return CellData::count; };
  MInt noBcSolverCellData() const {
    MInt count = 0;
    for(auto&& bc : m_boundaryConditions) {
      count += bc->noCellDataDlb();
    }
    return count;
  };

  // DLB reinitialization stage
  MInt m_loadBalancingReinitStage = -1;

  struct CellData;

  static const std::array<MInt, CellData::count> s_cellDataTypeDlb;

  // // void checkRhsForNan(MString tag);

  // // // Member variables
  // // //////////////////////////////////////////////////////////////////////////////
  // // // RULES FOR DEFAULT MEMBER INITIALIZATION
  //
  // If the member is an aggregate and cannot be aggregate-initialized,
  // value-initialize it (e.g. 'm_var()') in the member initializer list of the
  // constructor(s).
  // Else if there is a typical, sensible default value/constructor, use it.
  // Else use a value for initialization that, if unchanged, leads to an early
  // and deterministic abort so it can be detected and fixed immediately.
  //
  // REMINDER: Never leave a variable uninitialized!
  //////////////////////////////////////////////////////////////////////////////

  // Flow control
  // Current time step
  MInt m_timeStep = -1;
  // Total number of time steps
  MInt m_timeSteps = -1;
  // Identify first time step
  MBool m_firstTimeStep = false;
  // Specify number of time steps between recalculation of m_dt
  MInt m_calcTimeStepInterval = -1;
  // Current Runge Kutta stage
  MInt m_rkStage = -1;
  // Number of Runge Kutta stages
  MInt m_noRkStages = -1;
  // Time Integration Coefficient A
  std::vector<MFloat> m_timeIntegrationCoefficientsA{};
  // Time Integration Coefficient B
  std::vector<MFloat> m_timeIntegrationCoefficientsB{};
  // Time Integration Coefficient C
  std::vector<MFloat> m_timeIntegrationCoefficientsC{};
  // Identify if there are open MPI recv/send surface exchange requests
  MBool m_mpiRecvRequestsOpen = false;
  MBool m_mpiSendRequestsOpen = false;

  // Numerical properties
  // SBP Mode
  MBool m_sbpMode = false;
  // SBP Operator
  MString m_sbpOperator = "";
  // Initial polynomial degree for all elements
  MInt m_initPolyDeg = -1;
  // Minimum polynomial degree for all elements
  MInt m_minPolyDeg = -1;
  // Maximum polynomial degree for all elements
  MInt m_maxPolyDeg = -1;
  // Initial number of nodes 1D for all elements
  MInt m_initNoNodes1D = -1;
  // Minimum number of nodes 1D for all elements
  MInt m_minNoNodes1D = -1;
  // Maximum number of nodes 1D for all elements
  MInt m_maxNoNodes1D = -1;
  // Integration method for the DG method
  MInt m_dgIntegrationMethod = -1;
  // Time integration scheme for the DG method
  MInt m_dgTimeIntegrationScheme = -1;
  // Polynomial type for the DG method
  MInt m_dgPolynomialType = -1;
  // Physical time at which to start calculations
  MFloat m_startTime = 0.0;
  // Physical time at which to end calculations
  MFloat m_finalTime = 0.0;
  // Indicates if the final time step is reached
  MBool m_finalTimeStep = false;
  // Non-dimensional Courant-Friedrichs-Levy number
  MFloat m_cfl = 1.0;

  // Interpolation
  // Interpolation information (access by polynomial degree and number of nodes)
  std::vector<std::vector<DgInterpolation>> m_interpolation{};
  // Global domain volume (needed for error norm calculations)
  MFloat m_globalVolume = -1.0;
  // Local domain volume (needed for error norm calculations)
  MFloat m_localVolume = -1.0;
  // Specify whether to calculate error norms
  MBool m_calcErrorNorms = true;
  // Specify number of time steps between performance output
  /* MInt m_aliveInterval = -1; */
  // Specify number of time steps between two analyses
  MInt m_analysisInterval = -1;
  // Number of nodes (1D) to use for error analysis
  MInt m_noAnalysisNodes = -1;
  // Polynomial degree for error analysis
  MInt m_polyDegAnalysis = -1;
  // Number of significant digits in error analysis output
  MInt m_noErrorDigits = -1;
  // Interpolation information for error analysis
  DgInterpolation m_interpAnalysis{};
  // Volume weights for error analysis
  MFloatTensor m_wVolumeAnalysis{};
  // Vandermonde matrices for error analysis (access by polynomial degree
  // and number of nodes of original solution)
  std::vector<std::vector<MFloatTensor>> m_vdmAnalysis{};

  // Parallelization
  // Number of domains with which surface data needs to be exchanged
  MInt m_noExchangeNghbrDomains = -1;
  // Domain ids of neighbors with which surface data needs to be exchanged
  std::vector<MInt> m_exchangeNghbrDomains{};
  // Surface ids of MPI surfaces for each neighbor domain
  std::vector<std::vector<MInt>> m_mpiSurfaces{};
  // Surface ids of MPI surfaces for each neighbor domain for receiving (needed for periodic
  // boundary conditions)
  std::vector<std::vector<MInt>> m_mpiRecvSurfaces{};
#ifdef DG_USE_MPI_BUFFERS
  // Send buffers for each neighbor domain
  std::vector<std::vector<MFloat>> m_sendBuffers{};
  // Receive buffers for each neighbor domain
  std::vector<std::vector<MFloat>> m_recvBuffers{};
#endif
#ifdef DG_USE_MPI_DERIVED_TYPES
  // MPI datatypes for sending (access by neighbor domain id)
  std::vector<MPI_Datatype> m_sendTypes{};
  // MPI datatypes for receiving (access by neighbor domain id)
  std::vector<MPI_Datatype> m_recvTypes{};
#endif
  // Send requests for each neighbor domain
  std::vector<MPI_Request> m_sendRequests{};
  // Receive requests for each neighbor domain
  std::vector<MPI_Request> m_recvRequests{};

  // If true, always save solution file after last time step
  MBool m_alwaysSaveFinalSolution = true;
  // If true, always save restart file after last time step
  MBool m_alwaysSaveFinalRestart = true;
  // If true, save node variables to solution/restart files
  MBool m_saveNodeVariablesToSolutionFile = false;
  // Save time derivative
  MBool m_writeTimeDerivative = false;
  // Sponge eta will be included in data output
  MBool m_writeSpongeEta = false;
  // How many minutes before the job ends a restart file should be written
  MInt m_noMinutesEndAutoSave = 0;
  // How often (timesteps) to check if an end autosave file should be written
  MInt m_endAutoSaveCheckInterval = 0;
  // Update cell weights and min-cell workloads and save them to the grid file.
  MBool m_updateCellWeights = false;
  // Factor by which the DG contributions to the cell weights are scaled.
  MFloat m_weightPerNode = 1.0;
  // Additional weight per DG element independent of the polynomial degree.
  MFloat m_weightPerElement = 0.0;
  // Weight per cell
  MFloat m_weightPerCell = 0.0;
  // Restore default cell weights in grid file.
  MBool m_restoreDefaultWeights = false;
  // Save intial solution to solution file
  MBool m_writeInitialSolution = true;
  // Save node vars to node vars file
  MBool m_writeNodeVarsFile = true;

  // Auxiliary class for collecting and saving point data
  PointData<nDim, DgCartesianSolver<nDim, SysEqn>> m_pointData{*this};
  // Auxiliary class for collecting and saving surface data
  SurfaceData<nDim, DgCartesianSolver<nDim, SysEqn>> m_surfaceData{*this};
  // Auxiliary class for collecting and saving volume data
  VolumeData<nDim, DgCartesianSolver<nDim, SysEqn>> m_volumeData{*this};
  // Auxiliary class for creating 2D slices from 3D simulations
  DgSlices<nDim, SysEqn> m_slice{*this};

  // Grid map
  // Offsets to donor grid
  std::vector<maia::dg::GridMapOffset> m_gridMapOffsets{};
  // Maximum number of offsets on all domains
  MInt m_maxNoGridMapOffsets = -1;

  // Sponge
  DgSponge<nDim, SysEqn> m_sponge{-1, MPI_COMM_NULL};
  MBool m_useSponge = false;

  // H/P-refinement
  // If true, pRefinement is utilized based on additional properties
  MInt m_pref = -1;
  // Adaptive refinement setup
  MInt m_adaptivePref = -1;
  // Error threshold
  MFloat m_adaptiveThreshold = -1;
  // Interval for adaptive refinement
  MInt m_adaptiveInterval = -1;
  // Storage for mortar projections
  std::vector<MFloat*> m_projectionMatrixPointersH{};
  std::vector<MFloat> m_projectionMatricesH{};
  std::vector<MFloat*> m_projectionMatrixPointersP{};
  std::vector<MFloat> m_projectionMatricesP{};
  // Collector for h-refined elements
  HElementCollector m_helements;

  // Storage for the coordinates and polyDeg of p-refinement patches
  std::vector<std::array<MFloat, 2 * nDim>> m_prefPatchesCoords{};
  std::vector<MFloat> m_prefPatchesPolyDeg{};
  std::vector<MString> m_prefPatchesOperators{};
  std::vector<MInt> m_prefPatchesNoNodes1D{};

  // Geometry Intersection
  GeometryIntersection<nDim>* m_geometryIntersection = nullptr;

  // Physical properties and system of equations
  // System of equation-specific methods and variables
  SysEqn m_sysEqn;
  // Boundary condition factory to translate bcIds into objects
  DgBoundaryConditionFactory<nDim, SysEqn> m_boundaryConditionFactory{*this};
  // Container with boundary conditions
  std::vector<BC> m_boundaryConditions{};
  // Enable/disable support for cut-off boundaries
  MBool m_useCutOffBoundaries = false;
  // Store cut-off boundaries
  std::array<MInt, 2 * nDim> m_cutOffBoundaryConditionIds{};
  // Count number of surfaces with cut-off boundary conditions
  std::array<MInt, 2 * nDim> m_noCutOffBoundarySurfaces{};

  // Ramp up the external source term in time
  MBool m_useSourceRampUp = false;
  // Time at which the external source term is fully ramped up
  MFloat m_sourceRampUpTime = 0.0;

  // Memory management
  // Maximum number of surfaces that can be created
  MInt m_maxNoSurfaces = -1;
  // Total number of values (DOF * noVars * noInternalCells)
  MInt m_internalDataSize = -1;
  // Total number of nodes
  MInt m_noTotalNodesXD = -1;

  // Element variables
  // Collector of elements
  ElementCollector m_elements;

  // Surface variables
  // Collector of surfaces
  SurfaceCollector m_surfaces;
  // Number of boundary surfaces
  MInt m_noBoundarySurfaces = -1;
  // Number of inner surfaces
  MInt m_noInnerSurfaces = -1;
  // Number of MPI surfaces
  MInt m_noMpiSurfaces = -1;
  // Offset for (i.e. id of first) boundary surfaces
  //  MInt m_boundarySurfacesOffset = -1;
  // Offset for (i.e. id of first) inner surfaces
  MInt m_innerSurfacesOffset = -1;
  // Offset for (i.e. id of first) MPI surfaces
  MInt m_mpiSurfacesOffset = -1;

  // Time derivative calculation
  // Current simulation time
  MFloat m_time = 0.0;
  // Current time step size
  MFloat m_dt = -1.0;
  // Current external time step size set from outside the solver
  MFloat m_externalDt = -1.0;

  // Simulation statistics
  // Minimum polynomial degree on this domain
  MInt m_statLocalMinPolyDeg = -1;
  // Maximum polynomial degree on this domain
  MInt m_statLocalMaxPolyDeg = -1;
  // Minimum refinement level on this domain
  MInt m_statLocalMinLevel = -1;
  // Maximum refinement level on this domain
  MInt m_statLocalMaxLevel = -1;
  // Number of all cells on this domain
  MInt m_statLocalNoCells = -1;
  // Number of internal cells on this domain
  MInt m_statLocalNoInternalCells = -1;
  // Number of halo cells on this domain
  MInt m_statLocalNoHaloCells = -1;
  // Collector size for cells on this domain
  MInt m_statLocalMaxNoCells = -1;
  // Number of all elements on this domain
  MInt m_statLocalNoElements = -1;
  // Number of all helements on this domain
  MInt m_statLocalNoHElements = -1;
  // Number of all surfaces on this domain
  MInt m_statLocalNoSurfaces = -1;
  // Number of boundary surfaces on this domain
  MInt m_statLocalNoBoundarySurfaces = -1;
  // Number of inner surfaces on this domain
  MInt m_statLocalNoInnerSurfaces = -1;
  // Number of MPI surfaces on this domain
  MInt m_statLocalNoMpiSurfaces = -1;
  // Collector size for surfaces on this domain
  MInt m_statLocalMaxNoSurfaces = -1;
  // Number of active cells on this domain
  MInt m_statLocalNoActiveCells = -1;
  // Number of active degrees of freedom on this domain
  MInt m_statLocalNoActiveDOFs = -1;
  // Local number of active degrees of freedom for each polynomial degree
  std::vector<MInt> m_statLocalNoActiveDOFsPolyDeg{};
  // Sum of polynomial degrees for all elements on this domain
  MInt m_statLocalPolyDegSum = -1;
  // Sum of h-refined surfaces
  MInt m_statLocalNoHrefSurfs = -1;
  // Sum of p-refined surfaces
  MInt m_statLocalNoPrefSurfs = -1;
  // Minimum polynomial degree on all domains
  MInt m_statGlobalMinPolyDeg = -1;
  // Maximum polynomial degree on all domains
  MInt m_statGlobalMaxPolyDeg = -1;
  // Minimum refinement level on all domains
  MInt m_statGlobalMinLevel = -1;
  // Maximum refinement level on all domains
  MInt m_statGlobalMaxLevel = -1;
  // Number of all cells on all domains
  MLong m_statGlobalNoCells = -1;
  // Number of internal cells on all domains
  MLong m_statGlobalNoInternalCells = -1;
  // Number of halo cells on all domains
  MLong m_statGlobalNoHaloCells = -1;
  // Collector size for cells on all domains
  MLong m_statGlobalMaxNoCells = -1;
  // Number of all elements on all domains
  MLong m_statGlobalNoElements = -1;
  // Number of all surfaces on all domains
  MLong m_statGlobalNoSurfaces = -1;
  // Number of boundary surfaces on all domains
  MLong m_statGlobalNoBoundarySurfaces = -1;
  // Number of inner surfaces on all domains
  MLong m_statGlobalNoInnerSurfaces = -1;
  // Number of MPI surfaces on all domains
  MLong m_statGlobalNoMpiSurfaces = -1;
  // Collector size for surfaces on all domains
  MLong m_statGlobalMaxNoSurfaces = -1;
  // Number of active cells on all domains
  MLong m_statGlobalNoActiveCells = -1;
  // Number of active degrees of freedom on all domains
  MLong m_statGlobalNoActiveDOFs = -1;
  // Sum of polynomial degrees for all elements on this domain
  MInt m_statGlobalPolyDegSum = -1;
  // Sum of h-refined surfaces
  MLong m_statGlobalNoHrefSurfs = -1;
  // Sum of p-refined surfaces
  MLong m_statGlobalNoPrefSurfs = -1;
  // Sum of number of helements
  MLong m_statGlobalNoHElements = -1;
  // Sum of number of surfaces with cut-off boundary condition (max 6 for 3D)
  std::array<MLong, 6> m_statGlobalNoCutOffBoundarySurfaces{};

  // Timers
  // Timer group which holds all solver-wide timers
  MInt m_timerGroup = -1;
  // Stores all solver-wide timers
  std::array<MInt, Timers::_count> m_timers{};

  // Initialization status
  MBool m_isInitTimers = false;
  MBool m_isInitSolver = false;
  MBool m_isInitData = false;
  MBool m_isInitMainLoop = false;
  MBool m_isInitMpiExchange = false;

  // Variables related to main loop execution
  MInt m_noAnalyzeTimeSteps = -1;
  MFloat m_analyzeTimeStart = 0.0;
  MFloat m_loopTimeStart = 0.0;
  std::time_t m_endAutoSaveTime = -1;
  MBool m_endAutoSaveWritten = false;
  MFloat m_outputTime = 0.0;
};

// Ensure that only one type of MPI communication is enabled. For more
// information, please refer to config.h.
#if defined(DG_USE_MPI_BUFFERS) && defined(DG_USE_MPI_DERIVED_TYPES)
#error Both methods of doing MPI surface exchanges are enabled. Pick one.
#endif
#if !defined(DG_USE_MPI_BUFFERS) && !defined(DG_USE_MPI_DERIVED_TYPES)
#error No method of doing MPI surface exchanges is enabled. Pick one.
#endif

// Ensure that if DG_USE_MPI_DERIVED_TYPES and DG_USE_MPI_WAITSOME are
// not enabled at the same time
#if defined(DG_USE_MPI_DERIVED_TYPES) && defined(DG_USE_MPI_WAITSOME)
#error It makes no sense to use MPI_Waitsome with the derived data types.
#endif


// Struct to differentiate cell (element) data
template <MInt nDim, class SysEqn>
struct DgCartesianSolver<nDim, SysEqn>::CellData {
  static constexpr const MInt count = 5;

  static constexpr const MInt ELEM_CELL_ID = 0;
  static constexpr const MInt ELEM_POLY_DEG = 1;
  static constexpr const MInt ELEM_NO_NODES_1D = 2;
  static constexpr const MInt ELEM_VARIABLES = 3;
  static constexpr const MInt ELEM_NODE_VARS = 4;
};


// Data types of cell data
template <MInt nDim, class SysEqn>
const std::array<MInt, DgCartesianSolver<nDim, SysEqn>::CellData::count>
    DgCartesianSolver<nDim, SysEqn>::s_cellDataTypeDlb = {{MINT, MINT, MINT, MFLOAT, MFLOAT}};


/**
 * \brief Calculate the volume integral for all elements and update
 *        m_rightHandSide.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] noElements Number of elements.
 * \param[in] elem Pointer to elements.
 * \param[in] fluxFct Object which provides the calcFlux() function used to
 *                    compute the physical flux.
 */
template <MInt nDim, class SysEqn>
template <class F>
void DgCartesianSolver<nDim, SysEqn>::calcVolumeIntegral(const MInt noElements, ElementCollector& elem, F& fluxFct) {
  TRACE();

  const MInt noVars = SysEqn::noVars();
  const MInt* const polyDegs = &elem.polyDeg(0);
  const MInt* const noNodes = &elem.noNodes1D(0);

  // Create temporary storage for flux values
  const MInt maxNoNodesXD = elem.maxNoNodesXD();
  std::vector<MFloat> flux(maxNoNodesXD * noVars * nDim);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(flux)
#endif
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt polyDeg = polyDegs[elementId];
    const MInt noNodes1D = noNodes[elementId];
    const MFloat* const dhat = &m_interpolation[polyDeg][noNodes1D].m_Dhat[0];
    MFloat* const ut = &elem.rightHandSide(elementId);
    MInt index = 0;

    // Calculate flux
    fluxFct.calcFlux(&elem.nodeVars(elementId), &elem.variables(elementId), noNodes1D, &flux[0]);

    IF_CONSTEXPR(nDim == 2) {
      // Copy flux to time derivative
      MFloatTensor f(&flux[0], noNodes1D, noNodes1D, nDim, noVars);
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt n = 0; n < noVars; n++) {
            // auto df = 0.0;
            for(MInt l = 0; l < noNodes1D; l++) {
              ut[index] += dhat[i * noNodes1D + l] * f(l, j, 0, n) + dhat[j * noNodes1D + l] * f(i, l, 1, n);
            }
            index++;
          }
        }
      }
    }
    else IF_CONSTEXPR(nDim == 3) {
      // Copy flux to time derivative
      MFloatTensor f(&flux[0], noNodes1D, noNodes1D, noNodes1D, nDim, noVars);

      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D; k++) {
            for(MInt n = 0; n < noVars; n++) {
              for(MInt l = 0; l < noNodes1D; l++) {
                ut[index] += dhat[i * noNodes1D + l] * f(l, j, k, 0, n) + dhat[j * noNodes1D + l] * f(i, l, k, 1, n)
                             + dhat[k * noNodes1D + l] * f(i, j, l, 2, n);
              }
              index++;
            }
          }
        }
      }
    }
  }
}


/**
 * \brief Calculate the numerical flux for a regular (i.e. inner or MPI)
 *        surface.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] begin Index of the first surface to consider.
 * \param[in] end Index+1 of the last surface to consider.
 * \param[in] surf Pointer to surfaces.
 * \param[in] riemannFct Object which provides the calcRiemann() function used
 *                       to compute the numerical flux.
 */
template <MInt nDim, class SysEqn>
template <class F>
void DgCartesianSolver<nDim, SysEqn>::calcRegularSurfaceFlux(const MInt begin,
                                                             const MInt end,
                                                             SurfaceCollector& surf,
                                                             F& riemannFct) {
  TRACE();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  // Loop over all surfaces, calculate the numerical flux
  for(MInt srfcId = begin; srfcId < end; srfcId++) {
    MFloat* flux = &surf.flux(srfcId);
    MFloat* nodeVarsL = &surf.nodeVars(srfcId, 0);
    MFloat* nodeVarsR = &surf.nodeVars(srfcId, 1);
    MFloat* stateL = &surf.variables(srfcId, 0);
    MFloat* stateR = &surf.variables(srfcId, 1);
    const MInt noNodes1D = surf.noNodes1D(srfcId);
    const MInt dirId = surf.orientation(srfcId);

    riemannFct.calcRiemann(nodeVarsL, nodeVarsR, stateL, stateR, noNodes1D, dirId, flux);
  }
}


/**
 * \brief Calculates the source terms for each node and adds them to the time
 *        derivative of the conservative variables.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-23
 *
 * \param[in] t Current time.
 * \param[in] noElements Number of elements.
 * \param[in] elem Pointer to elements.
 * \param[in] sourceFct Object which provides the calcSource() function used to
 *                      calculate the source terms.
 */
template <MInt nDim, class SysEqn>
template <class F>
void DgCartesianSolver<nDim, SysEqn>::calcSourceTerms(MFloat t,
                                                      const MInt noElements,
                                                      ElementCollector& elem,
                                                      F& sourceFct) {
  TRACE();

  const MInt* const noNodes1D = &elem.noNodes1D(0);

  const MInt maxDataBlockSize = ipow(m_maxNoNodes1D, nDim) * SysEqn::noVars();
  std::vector<MFloat> sources(maxDataBlockSize);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(sources)
#endif
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt noNodesXD = elem.noNodesXD(elementId);
    const MInt dataBlockSize = noNodesXD * SysEqn::noVars();
    sourceFct.calcSource(&elem.nodeVars(elementId), &elem.variables(elementId), noNodes1D[elementId], t,
                         &elem.nodeCoordinates(elementId), &sources[0]);

    MFloat* const rhs = &elem.rightHandSide(elementId);
    for(MInt dataId = 0; dataId < dataBlockSize; dataId++) {
      rhs[dataId] += sources[dataId];
    }
  }
}


/// \brief Get solver cell (element) data for load balancing.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// This method returns the requested data for all elements, e.g. the polynomial
/// degree or the variables of all elements.
///
/// \param[in] dataId Requested data id.
/// \param[in] oldNoCells Current (old) number of cells before load balancing.
/// \param[in] bufferIdToCellId Mapping from buffer location to corresponding cell id.
/// \param[out] data Pointer to storage for requested data.
template <MInt nDim, class SysEqn>
template <typename DataType>
void DgCartesianSolver<nDim, SysEqn>::getCellDataDlb_(const MInt dataId,
                                                      const MInt oldNoCells,
                                                      const MInt* const bufferIdToCellId,
                                                      DataType* const data) {
  TRACE();

  // Check for unsorted cells, not supported yet, data needs to be resorted into buffers
  MInt prevCellId = -1;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt mapping = bufferIdToCellId[i];
    if(mapping != -1) {
      if(mapping <= prevCellId) {
        TERMM(1, "Error: assembling data buffers for unsorted cells not supported yet in DG solver.");
      }
      prevCellId = mapping;
    }
  }

  const MInt noElements = m_elements.size();

  if(dataId > -1 && dataId < noDgCartesianSolverCellData()) {
    // DG solver cell data
    switch(dataId) {
      case CellData::ELEM_CELL_ID:
        std::copy_n(&m_elements.cellId(0), noElements, data);
        break;
      case CellData::ELEM_POLY_DEG:
        std::copy_n(&m_elements.polyDeg(0), noElements, data);
        break;
      case CellData::ELEM_NO_NODES_1D:
        std::copy_n(&m_elements.noNodes1D(0), noElements, data);
        break;
      case CellData::ELEM_VARIABLES: {
        MInt dataOffset = 0;
        // Store all element variables in the data buffer
        for(MInt elementId = 0; elementId < noElements; elementId++) {
          const MInt noNodesXD = m_elements.noNodesXD(elementId);
          const MInt dataSize = noNodesXD * SysEqn::noVars();
          std::copy_n(&m_elements.variables(elementId), dataSize, &data[dataOffset]);
          dataOffset += dataSize;
        }
        break;
      }
      case CellData::ELEM_NODE_VARS: {
        MInt dataOffset = 0;
        // Store all element node variables in the data buffer
        for(MInt elementId = 0; elementId < noElements; elementId++) {
          const MInt noNodesXD = m_elements.noNodesXD(elementId);
          const MInt dataSize = noNodesXD * SysEqn::noNodeVars();
          std::copy_n(&m_elements.nodeVars(elementId), dataSize, &data[dataOffset]);
          dataOffset += dataSize;
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id (" + std::to_string(dataId) + ").");
        break;
    }
  } else if(dataId >= noDgCartesianSolverCellData() && dataId < noCellDataDlb()) {
    // Boundary condition cell data
    MInt offset = noDgCartesianSolverCellData();
    for(auto&& bc : m_boundaryConditions) {
      const MInt bcNoCellData = bc->noCellDataDlb();
      if(dataId >= offset && dataId < offset + bcNoCellData) {
        bc->getCellDataDlb(dataId - offset, data);
        break;
      }
      offset += bcNoCellData;
    }
  } else {
    // TODO labels:DG support exchange of CC mean vars data -> needs to be performed via gridcontroller!
    TERMM(1, "Invalid dataId (" + std::to_string(dataId) + ").");
  }
}


/// \brief Set solver cell (element) data.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// This method sets the given data for all elements, e.g. the polynomial
/// degree or the variables of all elements.
///
/// \param[in] dataId Requested data id.
/// \param[in] data Pointer to storage of data.
template <MInt nDim, class SysEqn>
template <typename DataType>
void DgCartesianSolver<nDim, SysEqn>::setCellDataDlb_(const MInt dataId, const DataType* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!grid().isActive()) {
    return;
  }

  const MInt noElements = m_elements.size();

  if(dataId > -1 && dataId < noDgCartesianSolverCellData()) {
    // Skip if this is the wrong reinitialization stage
    if(m_loadBalancingReinitStage != 0) {
      return;
    }

    // DG solver cell data
    switch(dataId) {
      case CellData::ELEM_CELL_ID:
        std::copy_n(data, noElements, &m_elements.cellId(0));
        break;
      case CellData::ELEM_POLY_DEG:
        std::copy_n(data, noElements, &m_elements.polyDeg(0));
        break;
      case CellData::ELEM_NO_NODES_1D:
        std::copy_n(data, noElements, &m_elements.noNodes1D(0));
        break;
      case CellData::ELEM_VARIABLES: {
        MInt dataOffset = 0;
        for(MInt elementId = 0; elementId < noElements; elementId++) {
          // NOTE: assumes polynomial degree is already set!
          const MInt noNodesXD = m_elements.noNodesXD(elementId);
          const MInt dataSize = noNodesXD * SysEqn::noVars();
          std::copy_n(&data[dataOffset], dataSize, &m_elements.variables(elementId));
          dataOffset += dataSize;
        }
        break;
      }
      case CellData::ELEM_NODE_VARS: {
        MInt dataOffset = 0;
        for(MInt elementId = 0; elementId < noElements; elementId++) {
          // NOTE: assumes polynomial degree is already set!
          const MInt noNodesXD = m_elements.noNodesXD(elementId);
          const MInt dataSize = noNodesXD * SysEqn::noNodeVars();
          std::copy_n(&data[dataOffset], dataSize, &m_elements.nodeVars(elementId));
          dataOffset += dataSize;
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
  } else if(dataId >= noDgCartesianSolverCellData() && dataId < noCellDataDlb()) {
    // Note: this should happen after boundary conditions are created, skip if
    // this is the wrong reinitialization stage
    if(m_loadBalancingReinitStage != 2) {
      return;
    }

    // Boundary condition cell data
    MInt offset = noDgCartesianSolverCellData();
    for(auto&& bc : m_boundaryConditions) {
      const MInt bcNoCellData = bc->noCellDataDlb();
      if(dataId >= offset && dataId < offset + bcNoCellData) {
        bc->setCellDataDlb(dataId - offset, data);
        break;
      }
      offset += bcNoCellData;
    }
  } else {
    // TODO labels:DG support exchange of CC mean vars data -> needs to be performed via gridcontroller!
    // Note: when setting the data for the first time the boundary conditions are not created yet
    // and the total cell data count does not match
    if(m_loadBalancingReinitStage == 2) {
      TERMM(1, "Invalid dataId: " + std::to_string(dataId));
    }
  }
}

#endif // DGSOLVER_H_
