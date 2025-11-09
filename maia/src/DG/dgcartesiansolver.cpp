// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgcartesiansolver.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include "COMM/mpioverride.h"
#include "GEOM/geometryelement.h"
#include "IO/logtable.h"
#include "IO/parallelio.h"
#include "UTIL/hilbert.h"
#include "UTIL/maiamath.h"
#include "dgcartesianmortar.h"
#include "globals.h"
#include "sbpcartesianmortar.h"
#include "typetraits.h"

// Enable/disable compiler-specific diagnostics
#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic push
// #pragma GCC diagnostic error "-Wold-style-cast"
#endif

using namespace std;
using namespace maia; // Omit once this is all done within 'maia' namespace

/// \brief Constructor of the DG solver reads properties and allocates solver
///       resources (as far as they are known at instantiation time).
///
/// \author Michael Schlottke
/// \date   October 2012
///
/// The constructor does NOT prepare the solver for starting a simulation,
/// however. You need to call initSolverObjects() to set up all necessary
/// methods/data structures when you want to run a simulation.
///
///
template <MInt nDim, class SysEqn>
DgCartesianSolver<nDim, SysEqn>::DgCartesianSolver(const MInt solverId_, GridProxy& gridProxy_,
                                                   Geometry<nDim>& geometry_, const MPI_Comm comm)
  : maia::CartesianSolver<nDim, DgCartesianSolver<nDim, SysEqn>>(solverId_, gridProxy_, comm, true),
    m_geometry(geometry_),
    m_sponge(solverId_, comm),
    m_sysEqn(solverId_) {
  TRACE();

  // Store currently used memory
  const MLong previouslyAllocated = allocatedBytes();

  // Initialize all solver-wide timers
  initTimers();

  // Initialize basic properties
  initDgCartesianSolver();

  // Input/output
  setInputOutputProperties();

  // Numerical method
  setNumericalProperties();

  // Allocate general DG solver memory
  allocateAndInitSolverMemory();

  // Initialize boundary condition factory
  dg::bc::factory::Init<nDim, SysEqn>::init(m_boundaryConditionFactory);

  if(gridProxy_.isActive()) {
    // Print information on used memory
    printAllocatedMemory(previouslyAllocated, "DgCartesianSolver (solverId = " + to_string(m_solverId) + ")",
                         mpiComm());
  }

  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}


///
/// \brief Frees resources that were allocated in the constructor and
/// that are not automatically released when the member variables are
/// destroyed.
///
/// \author Michael Schlottke
/// \date   October 2012
template <MInt nDim, class SysEqn>
DgCartesianSolver<nDim, SysEqn>::~DgCartesianSolver() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Destructor]);
  RECORD_TIMER_STOP(m_timers[Timers::Destructor]);
  RECORD_TIMER_STOP(m_timers[Timers::SolverType]);
}


/**
 * \brief Initialize all solver-wide timers and start the solver timer.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-24
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initTimers() {
  TRACE();

  // Invalidate timer ids
  std::fill(m_timers.begin(), m_timers.end(), -1);

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "DgCartesianSolver (solverId = " + to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::SolverType], "total object lifetime", m_timerGroup);
  RECORD_TIMER_START(m_timers[Timers::SolverType]);

  // Create & start constructor timer
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Constructor], "Constructor", m_timers[Timers::SolverType]);
  RECORD_TIMER_START(m_timers[Timers::Constructor]);

  // Create regular solver-wide timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Run], "run", m_timers[Timers::SolverType]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::RunInit], "run initialization", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitSolverObjects], "initialize solver infrastructure",
                         m_timers[Timers::RunInit]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitMortarProjection], "initMortarProjection",
                         m_timers[Timers::InitSolverObjects]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitData], "initialize solver data", m_timers[Timers::RunInit]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitialCondition], "initialCondition", m_timers[Timers::InitData]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitMainLoop], "initialize main loop", m_timers[Timers::RunInit]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MainLoop], "main loop", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcTimeStep], "calcTimeStep", m_timers[Timers::MainLoop]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ResetExternalSources], "resetExternalSources", m_timers[Timers::MainLoop]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::RungeKuttaStep], "timeStepRk", m_timers[Timers::MainLoop]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TimeDeriv], "calcDgTimeDerivative", m_timers[Timers::RungeKuttaStep]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ResetRHS], "resetRHS", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Prolong], "prolong to surfaces", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ForwardProjection], "forward projection", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SurfExchange], "MPI surface exchange", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SurfExchangeComm], "communication", m_timers[Timers::SurfExchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SECommSend], "send", m_timers[Timers::SurfExchangeComm]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SECommRecv], "recv", m_timers[Timers::SurfExchangeComm]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SurfExchangeCopy], "copy operations", m_timers[Timers::SurfExchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SECopySend], "send", m_timers[Timers::SurfExchangeCopy]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SECopyRecv], "recv", m_timers[Timers::SurfExchangeCopy]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SurfExchangeWait], "waiting", m_timers[Timers::SurfExchange]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SEWaitSend], "send", m_timers[Timers::SurfExchangeWait]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SEWaitRecv], "recv", m_timers[Timers::SurfExchangeWait]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::VolInt], "calcVolumeIntegral", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Flux], "flux calculation", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FluxBndry], "calcBoundarySurfaceFlux", m_timers[Timers::Flux]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FluxInner], "calcInnerSurfaceFlux", m_timers[Timers::Flux]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FluxMPI], "calcMpiSurfaceFlux", m_timers[Timers::Flux]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SurfInt], "surface integrals", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Jacobian], "applyJacobian", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Sources], "calcSourceTerms", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ExternalSources], "applyExternalSourceTerms", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Sponge], "calcSpongeTerms", m_timers[Timers::TimeDeriv]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TimeInt], "time integration", m_timers[Timers::RungeKuttaStep]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MainLoopIO], "I/O", m_timers[Timers::MainLoop]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Analysis], "solution analysis", m_timers[Timers::MainLoop]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CleanUp], "cleanUp", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Destructor], "Destructor", m_timers[Timers::SolverType]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::AdaptiveRefinement], "adaptive refinement", m_timers[Timers::MainLoop]);

  // Create accumulated timers that monitor selected subsystems
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Accumulated], "selected accumulated timers", m_timers[Timers::SolverType]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::IO], "IO", m_timers[Timers::Accumulated]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SaveSolutionFile], "saveSolutionFile", m_timers[Timers::IO]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SaveRestartFile], "saveRestartFile", m_timers[Timers::IO]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::AnalyzeTimeStep], "analyzeTimeStep", m_timers[Timers::Accumulated]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MPI], "MPI", m_timers[Timers::Accumulated]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MPIComm], "communication", m_timers[Timers::MPI]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MPICopy], "copy operations", m_timers[Timers::MPI]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::MPIWait], "waiting", m_timers[Timers::MPI]);

  // Set status to initialized
  m_isInitTimers = true;
}


/** \fn void DgCartesianSolver<nDim, SysEqn>::initDgCartesianSolver()
 * \brief Initializes basic properties of the solver.
 *
 * \author Michael Schlottke
 * \date October 2012
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initDgCartesianSolver() {
  TRACE();

  /*! \property
    \page propertyPageDG DG
    \section maxNoSurfaces
    <code>MInt DgCartesianSolver::m_maxNoSurfaces </code>\n
    default = <code>none</code>\n \n
    Set the maximum number of surfaces that may be created. Default value is -1,
    in which case the number of surfaces is determined automatically on solver
    startup.\n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>SURFACES</i>
  */

  m_maxNoSurfaces = -1;
  m_maxNoSurfaces = Context::getSolverProperty<MInt>("maxNoSurfaces", m_solverId, AT_, &m_maxNoSurfaces);

  /*! \property
    \page propertyPageDG DG
    \section timeSteps
    <code>MInt DgCartesianSolver::m_timeSteps </code>\n
    default = <code>none</code>\n \n
    Set the maximum number of time steps to run the simulation for.\n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>TIME STEPS, FLOW CONTROL</i>
  */

  m_timeSteps = Context::getSolverProperty<MInt>("timeSteps", m_solverId, AT_);
}


/**
 * \details Reads properties associated with input/output.
 *
 * \author Michael Schlottke
 * \date   October 2012
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::setInputOutputProperties() {
  TRACE();

  /*! \property
    \page propertyPageDG DG
    \section alwaysSaveFinalSolution
    <code>MBool DgCartesianSolver::m_alwaysSaveFinalSolution</code>\n
    default = <code>true if m_solutionInterval > 1, otherwise false</code>\n\n
    If enabled, a solution file is saved after the last time step, no matter
    which value is set for m_solutionInterval.\n
    Possible values are:
    <ul>
      <li>0 - file is not always saved</li>
      <li>1 - file is always saved</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O, SOLUTION</i>
  */
  m_alwaysSaveFinalSolution = m_solutionInterval > 0;
  m_alwaysSaveFinalSolution =
      Context::getSolverProperty<MBool>("alwaysSaveFinalSolution", m_solverId, AT_, &m_alwaysSaveFinalSolution);

  /*! \property
    \page propertyPageDG DG
    \section alwaysSaveFinalRestart
    <code>MBool DgCartesianSolver::m_alwaysSaveFinalRestart</code>\n
    default = <code>true if m_restartInterval > 1, otherwise false</code>\n\n
    If enabled, a restart file is saved after the last time step, no matter
    which value is set for m_restartInterval.\n
    Possible values are:
    <ul>
      <li>0 - file is not always saved</li>
      <li>1 - file is always saved</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O, RESTART</i>
  */
  m_alwaysSaveFinalRestart = m_restartInterval > 0;
  m_alwaysSaveFinalRestart =
      Context::getSolverProperty<MBool>("alwaysSaveFinalRestart", m_solverId, AT_, &m_alwaysSaveFinalRestart);

  /*! \property
    \page propertyPageDG DG
    \section saveNodeVariablesToSolutionFile
    <code>MBool DgCartesianSolver::m_saveNodeVariablesToSolutionFile</code>\n
    default = <code>false</code>\n\n
    If enabled, the node variables will be saved to the solution files as well.
    They are always saved to the restart files.\n
    Possible values are:
    <ul>
      <li>0 - node variables are not stored in solution file</li>
      <li>1 - node variables are written to solution file</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O, SOLUTION, NODE_VARIABLES</i>
  */
  m_saveNodeVariablesToSolutionFile = false;
  m_saveNodeVariablesToSolutionFile = Context::getSolverProperty<MBool>("saveNodeVariablesToSolutionFile", m_solverId,
                                                                        AT_, &m_saveNodeVariablesToSolutionFile);

  /*! \property
    \page propertyPageDG DG
    \section analysisInterval
    <code>MInt DgCartesianSolver::m_analysisInterval</code>\n
    default = <code>10</code>\n \n
    Set the number of time steps after which an analysis step is performed. If
    the analysis interval is set to 0, no analyses are performed, except
    before the first time step and after the last time step.\n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, ANALYSIS, ERROR_NORM</i>
  */
  m_analysisInterval = 0;
  // m_analysisInterval
  m_analysisInterval = Context::getSolverProperty<MInt>("analysisInterval", m_solverId, AT_, &m_analysisInterval);
  if(m_analysisInterval < 0) {
    TERMM(1, "Analysis interval must be >= 0 (is: " + to_string(m_analysisInterval) + ").");
  }

  // Initialize I/O properties of helper class for point data gathering & saving
  m_pointData.setInputOutputProperties();
  // ..., for surface data
  m_surfaceData.setInputOutputProperties();
  // ... and for volume data
  m_volumeData.setInputOutputProperties();

  // Initialize I/O properties of dg slices
  m_slice.setProperties();

  /*! \property
    \page propertyPageDG DG
    \section aliveInterval
    <code>MInt DgCartesianSolver::m_aliveInterval</code>\n
    default = <code>AnalysisInterval / 10</code>\n \n
    Set the number of time steps after which performance information is given
    to the user. If the alive interval is set to 0, no alive signals are given.
    \n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN</i>
  */
  /* m_aliveInterval = m_analysisInterval / 10; */
  /* m_aliveInterval = */
  /*     Context::getSolverProperty<MInt>("aliveInterval", m_solverId, AT_, &m_aliveInterval); */
  /* if (m_aliveInterval < 0) { */
  /*   TERMM(1, "Alive interval must be >= 0 (is: " + to_string(m_analysisInterval) + ")."); */
  /* } */

  /*! \property
      \page propertyPageDG DG
      \section writeSpongeEta
      <code>MInt DgCartesianSolver::m_writeSpongeEta</code>\n
      default = false \n\n
      If not 0 spongeEta can be viewed with paraview.\n\n
      Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  m_writeSpongeEta = false;
  m_writeSpongeEta = Context::getSolverProperty<MBool>("writeSpongeEta", m_solverId, AT_, &m_writeSpongeEta);

  /*! \property
    \page propertyPageDG DG
    \section writeTimeDerivative
    <code>MInt DgCartesianSolver::m_writeTimeDerivative</code>\n
    default = 0 \n\n
    If 0 time derivative will not be saved. \n\n
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */

  m_writeTimeDerivative = 0;
  m_writeTimeDerivative =
      Context::getSolverProperty<MBool>("writeTimeDerivative", m_solverId, AT_, &m_writeTimeDerivative);

  /*!\property
     \page propertyPageDG DG
     \section noMinutesEndAutoSave
     <code> MInt DgCartesianSolver::m_noMinutesEndAutoSave</code>\
     Sets how many minutes before the job ends a restart file
     should be written (DEFAULT: 10).\n
     Only works if the environment variable MAIA_JOB_END_TIME contains
     the end time of the job in Unix time.
     Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
   */
  m_noMinutesEndAutoSave = 10;
  m_noMinutesEndAutoSave =
      Context::getSolverProperty<MInt>("noMinutesEndAutoSave", m_solverId, AT_, &m_noMinutesEndAutoSave);

  /*! \property
   *  \page propertyPageDG DG
   *  \section endAutoSaveCheckInterval
   *  <code> MInt DgCartesianSolver::m_endAutoSaveCheckInterval</code>\
   *  Sets how often to check if the jobs is about end.
   *  The interval is measured in timesteps.
   *  Keywords <i>DISCONTINUOUS_GALERKIN, I/0</i>
   */
  m_endAutoSaveCheckInterval = 50;
  m_endAutoSaveCheckInterval =
      Context::getSolverProperty<MInt>("endAutoSaveCheckInterval", m_solverId, AT_, &m_endAutoSaveCheckInterval);

  /*! \property
   *  \page propertyPageDG DG
   *  \section updateCellWeights
   *  <code>MInt DgCartesianSolver::m_updateCellWeights</code>\n
   *  default = <code>0</code>\n
   *  Activate computation of new cell weights and min-cell workloads and update the grid file. In
   *  the next solver run this will affect the grid partitioning in order to obtain a more balanced
   *  workload.
   *  <ul>
   *  <li><code>0</code> deactivated</li>
   *  <li><code>1</code> activated</li>
   *  </ul>
   *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
   */
  /* m_updateCellWeights = false; */
  /* m_updateCellWeights = Context::getSolverProperty<MBool>("updateCellWeights", */
  /*                                   m_solverId, AT_, &m_updateCellWeights); */

  /*! \property
   *  \page propertyPageDG DG
   *  \section weightDgRbcElements
   *  <code>MInt DgCartesianSolver::m_weightDgRbcElements</code>\n
   *  default = <code>0</code>\n
   *  Compute and add weight for RBC elements when using load balancing.
   *  <ul>
   *  <li><code>0</code> deactivated</li>
   *  <li><code>1</code> activated</li>
   *  </ul>
   *  Keywords: <i>DISCONTINUOUS_GALERKIN, DYNAMIC LOAD BALANCING</i>
   */
  m_weightDgRbcElements = false;
  m_weightDgRbcElements =
      Context::getSolverProperty<MBool>("weightDgRbcElements", m_solverId, AT_, &m_weightDgRbcElements);

  /* if (m_updateCellWeights) { */
  {
    /*! \property
     *  \page propertyPageDG DG
     *  \section weightPerNode
     *  <code>MFloat DgCartesianSolver::m_weightPerNode</code>\n
     *  default = <code>1.0</code>\n
     *  Factor by which the DG contributions to the cell weights are scaled if updateCellWeights is
     *  activated.
     *  <ul>
     *  <li><code>weight &ge; 0.0</code> deactivated</li>
     *  </ul>
     *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
     */
    m_weightPerNode = 1.0;
    m_weightPerNode = Context::getSolverProperty<MFloat>("weightPerNode", m_solverId, AT_, &m_weightPerNode);

    /*! \property
     *  \page propertyPageDG DG
     *  \section weightPerElement
     *  <code> DgCartesianSolver::m_weightPerElement</code>\n
     *  default = <code>0.0</code>\n
     *  Additional weight per DG element which is independent of the polynomial degree (see
     *  weightPerNode and updateCellWeights).
     *  <ul>
     *  <li><code>weight &ge; 0.0</code> deactivated</li>
     *  </ul>
     *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
     */
    m_weightPerElement = 0.0;
    m_weightPerElement = Context::getSolverProperty<MFloat>("weightPerElement", m_solverId, AT_, &m_weightPerElement);

    /*! \property
     *  \page propertyPageDG DG
     *  \section restoreDefaultWeights
     *  <code>MInt DgCartesianSolver::m_restoreDefaultWeights</code>\n
     *  default = <code>0</code>\n
     *  Determines if the default cell weights are restored, i.e. each cell is assigned the weight
     *  1.0 (updateCellWeights needs to be activated).
     *  <ul>
     *  <li><code>0</code> deactivated</li>
     *  <li><code>1</code> activated</li>
     *  </ul>
     *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
     */
    /* m_restoreDefaultWeights = false; */
    /* m_restoreDefaultWeights = Context::getSolverProperty<MBool>("restoreDefaultWeights", */
    /*                                 m_solverId, AT_, &m_restoreDefaultWeights); */
  }

  /*! \property
   *  \page propertyPageDG DG
   *  \section writeInitialSolution
   *  <code>MInt DgCartesianSolver::m_writeInitialSolution</code>\n
   *  default = <code>1</code>\n
   *  Should the initial solution be written to a solution file?
   *  <ul>
   *  <li><code>0</code> deactivated</li>
   *  <li><code>1</code> activated</li>
   *  </ul>
   *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
   */
  m_writeInitialSolution = true;
  m_writeInitialSolution =
      Context::getSolverProperty<MBool>("writeInitialSolution", m_solverId, AT_, &m_writeInitialSolution);

  /*! \property
   *  \page propertyPageDG DG
   *  \section writeNodeVarsFile
   *  <code>MInt DgCartesianSolver::m_writeNodeVarsFile</code>\n
   *  default = <code>1</code>\n
   *  Write the node vars to a node Vars file?
   *  <ul>
   *  <li><code>0</code> deactivated</li>
   *  <li><code>1</code> activated</li>
   *  </ul>
   *  Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
   */
  m_writeNodeVarsFile = true;
  m_writeNodeVarsFile = Context::getSolverProperty<MBool>("writeNodeVarsFile", m_solverId, AT_, &m_writeNodeVarsFile);
}


/** \fn void DgCartesianSolver<nDim, SysEqn>::setNumericalProperties()
 * \brief Reads properties associated with the numerical method.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::setNumericalProperties() {
  TRACE();

  /*! \property
      \page propertyPageDG DG
      \section startTime
      <code>MInt DgCartesianSolver::m_startTime </code>\n
      default = <code>no default</code>\n \n
      Starting time of the calculation. \n \n
      Possible values are:
      <ul>
      <li> Positive floating point values. </li>
      </ul>
      Keywords: <i> GENERAL, TIME_INTEGRATION </i>
  */
  m_startTime = Context::getSolverProperty<MFloat>("startTime", m_solverId, AT_);

  /*! \property
      \page propertyPageDG DG
      \section finalTime
      <code>MInt DgCartesianSolver::m_finalTime </code>\n
      default = <code>no default</code>\n \n
      Final time of the calculation. \n \n
      Possible values are:
      <ul>
      <li> Positive floating point values. </li>
      </ul>
      Keywords: <i> GENERAL, TIME_INTEGRATION </i>
  */
  m_finalTime = Context::getSolverProperty<MFloat>("finalTime", m_solverId, AT_);

  if(m_startTime > m_finalTime || m_startTime < F0) {
    cerr << "Illegal setup of the time integration parameters. Aborting..." << endl;
    m_log << "Illegal setup of the time integration parameters. Aborting..." << endl;
    mTerm(1, AT_, "Illegal setup of the time integration parameters.");
  }

  /*! \property
      \page propertyPageDG DG
      \section calcTimeStepInterval
      <code>MInt DgCartesianSolver::m_calcTimeStepInterval </code>\n
      default = <code>0</code>\n \n
      Set the number of time steps after which size of timeStep is
      recalculated. A zero-value implies that the time step is only
      calculated once.\n
      Possible values are:
      <ul>
      <li> integer >= 0 </li>
      </ul>
      Keywords: <i> GENERAL, TIME_INTEGRATION </i>
  */
  m_calcTimeStepInterval = 0;
  m_calcTimeStepInterval =
      Context::getSolverProperty<MInt>("calcTimeStepInterval", m_solverId, AT_, &m_calcTimeStepInterval);
  if(m_calcTimeStepInterval < 0) {
    TERMM(1, "Recalculation interval for time step must be >= 0 (is: " + to_string(m_calcTimeStepInterval) + ").");
  }

  /*! \property
      \page propertyPageDG DG
      \section sbpMode
      <code>MBool DgCartesianSolver::m_sbpMode </code>\n
      default = <code>0</code>\n \n
      Exchanges DG operators with SBP operators.\n\n
      Keywords: <i> SBP </i>
  */
  m_sbpMode = false;
  m_sbpMode = Context::getSolverProperty<MBool>("sbpMode", m_solverId, AT_, &m_sbpMode);

  /*! \property
      \page propertyPageDG DG
      \section sbpOperator
      <code>MString DgCartesianSolver::m_sbpOperator </code>\n
      default = <code>""</code>\n \n
      Defines which SBP Operator is to be used.
      Keywords: <i> SBP </i>
  */
  m_sbpOperator = "";
  m_sbpOperator = Context::getSolverProperty<MString>("sbpOperator", m_solverId, AT_, &m_sbpOperator);

  /*! \property
      \page propertyPageDG DG
      \section dgIntegrationMethod
      <code>MInt DgCartesianSolver::m_dgIntegrationMethod </code>\n
      default = <code>DG_INTEGRATE_GAUSS</code>\n \n
      Specifies the (spatial) integration method to be used in the DG solver.\n\n
      Possible values are:
      <ul>
      <li> DG_INTEGRATE_GAUSS </li>
      <li> DG_INTEGRATE_GAUSS_LOBATTO </li>
      </ul>
      Keywords: <i> DISCONTINUOUS_GALERKIN, SPATIAL_INTEGRATION </i>
  */

  MString dgIntegrationMethod = m_sbpMode ? "DG_INTEGRATE_GAUSS_LOBATTO" : "DG_INTEGRATE_GAUSS";
  m_dgIntegrationMethod =
      string2enum(Context::getSolverProperty<MString>("dgIntegrationMethod", m_solverId, AT_, &dgIntegrationMethod));

  /*!  \property
       \page propertyPageDG DG
       \section dgTimeIntegrationScheme
       <code>MInt DgCartesianSolver::m_dgTimeIntegrationScheme </code>\n
       default = <code>DG_TIMEINTEGRATION_NIEGEMANN_4_14</code>\n \n
       Specifies the (temporal) integration method to be used in the DG solver.\n\n
       Possible values are:
       <ul>
       <li> DG_TIMEINTEGRATION_CARPENTER_4_5 </li>
       <li> DG_TIMEINTEGRATION_TOULORGEC_4_8 </li>
       <li> DG_TIMEINTEGRATION_NIEGEMANN_4_14 </li>
       <li> DG_TIMEINTEGRATION_NIEGEMANN_4_13 </li>
       <li> DG_TIMEINTEGRATION_TOULORGEC_3_7 </li>
       <li> DG_TIMEINTEGRATION_TOULORGEF_4_8 </li>
       </ul>
       Keywords: <i> DISCONTINUOUS_GALERKIN, TEMPORAL_INTEGRATION </i>
   */
  MString dgTimeIntegrationScheme = "DG_TIMEINTEGRATION_CARPENTER_4_5";
  m_dgTimeIntegrationScheme = string2enum(
      Context::getSolverProperty<MString>("dgTimeIntegrationScheme", m_solverId, AT_, &dgTimeIntegrationScheme));

  // Define the number of Rk-steps according to the scheme being used and it's coeffiecients
  m_noRkStages = -1;
  switch(m_dgTimeIntegrationScheme) {
    case DG_TIMEINTEGRATION_CARPENTER_4_5: {
      default:
        m_noRkStages = 5;

        m_timeIntegrationCoefficientsA.resize(m_noRkStages);
        m_timeIntegrationCoefficientsB.resize(m_noRkStages);
        m_timeIntegrationCoefficientsC.resize(m_noRkStages);

        m_timeIntegrationCoefficientsA = {{F0, 567301805773.0 / 1357537059087.0, 2404267990393.0 / 2016746695238.0,
                                           3550918686646.0 / 2091501179385.0, 1275806237668.0 / 842570457699.0}};
        m_timeIntegrationCoefficientsB = {{1432997174477.0 / 9575080441755.0, 5161836677717.0 / 13612068292357.0,
                                           1720146321549.0 / 2090206949498.0, 3134564353537.0 / 4481467310338.0,
                                           2277821191437.0 / 14882151754819.0}};
        m_timeIntegrationCoefficientsC = {{F0, 1432997174477.0 / 9575080441755.0, 2526269341429.0 / 6820363962896.0,
                                           2006345519317.0 / 3224310063776.0, 2802321613138.0 / 2924317926251.0}};
        break;
    }

    case DG_TIMEINTEGRATION_TOULORGEC_4_8: {
      m_noRkStages = 8;

      m_timeIntegrationCoefficientsA.resize(m_noRkStages);
      m_timeIntegrationCoefficientsB.resize(m_noRkStages);
      m_timeIntegrationCoefficientsC.resize(m_noRkStages);

      m_timeIntegrationCoefficientsA = {{F0, 0.7212962482279240, 0.0107733657161298, 0.5162584698930970,
                                         1.7301002866322010, 5.2001293044030760, -0.7837058945416420,
                                         0.5445836094332190}};
      m_timeIntegrationCoefficientsB = {{0.2165936736758085, 0.1773950826411583, 0.0180253861162329, 0.0847347637254149,
                                         0.8129106974622483, 1.9034160304227600, 0.1314841743399048,
                                         0.2082583170674149}};
      m_timeIntegrationCoefficientsC = {{F0, 0.2165936736758085, 0.2660343487538170, 0.2840056122522720,
                                         0.3251266843788570, 0.4555149599187530, 0.7713219317101170,
                                         0.9199028964538660}};
      break;
    }

    case DG_TIMEINTEGRATION_NIEGEMANN_4_14: {
      m_noRkStages = 14;

      m_timeIntegrationCoefficientsA.resize(m_noRkStages);
      m_timeIntegrationCoefficientsB.resize(m_noRkStages);
      m_timeIntegrationCoefficientsC.resize(m_noRkStages);

      m_timeIntegrationCoefficientsA = {{F0, 0.718801210867241, 0.778533117342157, 0.005328279665404, 0.855297993402928,
                                         3.95641382457746, 1.57805753805874, 2.08370945525741, 0.748333418276161,
                                         0.703286110656336, -0.001391709611768, 0.093207536963746, 0.951420047087595,
                                         7.11515716939226}};
      m_timeIntegrationCoefficientsB = {{0.036776245431967, 0.313629660755396, 0.153184869186903, 0.003009708681818,
                                         0.332629379064611, 0.244025140535086, 0.371887923959228, 0.620412622158244,
                                         0.152404317302874, 0.076089492741927, 0.007760421404098, 0.002464728475538,
                                         0.078034834004939, 5.50597772702696}};
      m_timeIntegrationCoefficientsC = {{F0, 0.036776245431967, 0.124968526272503, 0.244617770227770, 0.247614953107042,
                                         0.296931112038247, 0.397814964580264, 0.527085458944033, 0.698126999417570,
                                         0.819089083535213, 0.852705988709862, 0.860471181746282, 0.862706037696998,
                                         0.873421312760098}};
      break;
    }

    case DG_TIMEINTEGRATION_NIEGEMANN_4_13: {
      m_noRkStages = 13;

      m_timeIntegrationCoefficientsA.resize(m_noRkStages);
      m_timeIntegrationCoefficientsB.resize(m_noRkStages);
      m_timeIntegrationCoefficientsC.resize(m_noRkStages);

      m_timeIntegrationCoefficientsA = {{F0, 0.6160178650170565, 0.4449487060774118, 1.0952033345276178,
                                         1.2256030785959187, 0.2740182222332805, 0.0411952089052647, 0.1797084899153560,
                                         1.1771530652064288, 0.4078831463120878, 0.8295636426191777, 4.7895970584252288,
                                         0.6606671432964504}};
      m_timeIntegrationCoefficientsB = {{0.0271990297818803, 0.1772488819905108, 0.0378528418949694, 0.6086431830142991,
                                         0.2154313974316100, 0.2066152563885843, 0.0415864076069797, 0.0219891884310925,
                                         0.9893081222650993, 0.0063199019859826, 0.3749640721105318, 1.6080235151003195,
                                         0.0961209123818189}};
      m_timeIntegrationCoefficientsC = {{F0, 0.0271990297818803, 0.0952594339119365, 0.1266450286591127,
                                         0.1825883045699772, 0.3737511439063931, 0.5301279418422206, 0.5704177433952291,
                                         0.5885784947099155, 0.6160769826246714, 0.6223252334314046, 0.6897593128753419,
                                         0.9126827615920843}};
      break;
    }

    case DG_TIMEINTEGRATION_TOULORGEC_3_7: {
      m_noRkStages = 7;

      m_timeIntegrationCoefficientsA.resize(m_noRkStages);
      m_timeIntegrationCoefficientsB.resize(m_noRkStages);
      m_timeIntegrationCoefficientsC.resize(m_noRkStages);

      m_timeIntegrationCoefficientsA = {{F0, 0.808316387498383, 1.503407858773331, 1.053064525050744, 1.463149119280508,
                                         0.659288128108783, 1.667891931891068}};
      m_timeIntegrationCoefficientsB = {{0.0119705267309784, 0.8886897793820711, 0.4578382089261419, 0.5790045253338471,
                                         0.3160214638138484, 0.2483525368264122, 0.0677123095940884}};
      m_timeIntegrationCoefficientsC = {{F0, 0.0119705267309784, 0.1823177940361990, 0.5082168062551849,
                                         0.6532031220148590, 0.8534401385678250, 0.9980466084623790}};
      break;
    }

    case DG_TIMEINTEGRATION_TOULORGEF_4_8: {
      m_noRkStages = 8;

      m_timeIntegrationCoefficientsA.resize(m_noRkStages);
      m_timeIntegrationCoefficientsB.resize(m_noRkStages);
      m_timeIntegrationCoefficientsC.resize(m_noRkStages);

      m_timeIntegrationCoefficientsA = {{F0, 0.5534431294501569, -0.0106598757020349, 0.5515812888932000,
                                         1.8857903775587410, 5.7012957427932640, -2.1139039656647930,
                                         0.5339578826675280}};
      m_timeIntegrationCoefficientsB = {{0.0803793688273695, 0.5388497458569843, 0.0197497440903196, 0.0991184129733997,
                                         0.7466920411064123, 1.6795842456188940, 0.2433728067008188,
                                         0.1422730459001373}};
      m_timeIntegrationCoefficientsC = {{F0, 0.0803793688273695, 0.3210064250338430, 0.3408501826604660,
                                         0.3850364824285470, 0.5040052477534100, 0.6578977561168540,
                                         0.9484087623348481}};
      break;
    }
  }

  // Sanity check to ensure a valid time integration scheme has been selected
  if(m_noRkStages == -1) {
    TERMM(1, "Invalid property value for time integration scheme");
  }

  MString dgPolynomialType = "DG_POLY_LEGENDRE";
  m_dgPolynomialType =
      string2enum(Context::getSolverProperty<MString>("dgPolynomialType", m_solverId, AT_, &dgPolynomialType));

  /*! \property
    \page propertyPageDG DG
    \section initPolyDeg
    <code>MInt DgCartesianSolver::m_initPolyDeg </code>\n
    default = <code>no default value</code>\n \n
    Set the initial polynomial order in the elements.\n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, POLYNOMIAL_DEGREE</i>
  */
  m_initPolyDeg = -1;
  m_initPolyDeg = Context::getSolverProperty<MInt>("initPolyDeg", m_solverId, AT_, &m_initPolyDeg);

  if(m_initPolyDeg < 0) {
    mTerm(1, AT_, "Negative polynomial order set in properties.");
  }


  /*! \property
    \page propertyPageDG DG
    \section minPolyDeg
    <code>MInt DgCartesianSolver::m_minPolyDeg </code>\n
    default = <code>m_initPolyDeg</code>\n \n
    Set the minimum polynomial degree for p-refinement.
    If the value is not set or does not make sense with respect to
    the initial polynomial degree, it is automatically reset.\n
    Possible values are:
    <ul>
      <li>any integer >= 0 and <= m_maxPolyDeg</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, POLYNOMIAL_DEGREE</i>
  */
  m_minPolyDeg = numeric_limits<MInt>::max();
  m_minPolyDeg = Context::getSolverProperty<MInt>("minPolyDeg", m_solverId, AT_, &m_minPolyDeg);

  /*! \property
    \page propertyPageDG DG
    \section maxPolyDeg
    <code>MInt DgCartesianSolver::m_maxPolyDeg </code>\n
    default = <code>m_initPolyDeg</code>\n \n
    Set the maximum polynomial degree for p-refinement.
    If the value is not set or does not make sense with respect to
    the initial polynomial degree, it is automatically reset.\n
    Possible values are:
    <ul>
      <li>any integer >= 0 and >= m_minPolyDeg</li>
    </dgsolver.cpp:843:    m_noAnalysisNodes = 2 * (m_maxPolyDeg + 1);
    ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, POLYNOMIAL_DEGREE</i>
  */
  m_maxPolyDeg = -1;
  m_maxPolyDeg = Context::getSolverProperty<MInt>("maxPolyDeg", m_solverId, AT_, &m_maxPolyDeg);

  // Make sure that all values for the polynomimal order are consistent
  // (for Gauss-Lobatto, the minimum polynomial degree is 1)
  if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
    m_initPolyDeg = max(1, m_initPolyDeg);
    m_maxPolyDeg = max(1, m_maxPolyDeg);
    m_minPolyDeg = max(1, m_minPolyDeg);
  } else {
    m_initPolyDeg = max(0, m_initPolyDeg);
    m_maxPolyDeg = max(0, m_maxPolyDeg);
    m_minPolyDeg = max(0, m_minPolyDeg);
  }

  m_minPolyDeg = min(m_minPolyDeg, m_initPolyDeg);
  m_maxPolyDeg = max(m_maxPolyDeg, m_initPolyDeg);

  // Determine range in noNodes1D which is needed for SBP mode
  if(m_sbpMode) {
    m_initNoNodes1D = -1;
    m_initNoNodes1D = Context::getSolverProperty<MInt>("initNoNodes", m_solverId, AT_, &m_initNoNodes1D);

    m_minNoNodes1D = Context::getSolverProperty<MInt>("minNoNodes", m_solverId, AT_, &m_initNoNodes1D);

    m_maxNoNodes1D = Context::getSolverProperty<MInt>("maxNoNodes", m_solverId, AT_, &m_initNoNodes1D);

    // Make sure that all values for the number of nodes  are consistent
    m_initNoNodes1D = max(2, m_initNoNodes1D);
    m_maxNoNodes1D = max(2, m_maxNoNodes1D);
    m_minNoNodes1D = max(2, m_minNoNodes1D);

    m_minNoNodes1D = min(m_minNoNodes1D, m_initNoNodes1D);
    m_maxNoNodes1D = max(m_maxNoNodes1D, m_initNoNodes1D);

  } else { // DG
    m_initNoNodes1D = m_initPolyDeg + 1;
    m_minNoNodes1D = m_minPolyDeg + 1;
    m_maxNoNodes1D = m_maxPolyDeg + 1;
  }

  /*! \property
    \page propertyPageDG DG
    \section calcErrorNorms
    <code>MInt DgCartesianSolver::m_calcErrorNorms </code>\n
    default = <code>1</code>\n \n
    Specify whether to calculate error norms\n
    Possible values are:
    <ul>
      <li>0 - no error norm calculation</li>
      <li>1 - calculate error norms</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, POLYNOMIAL_DEGREE, ERROR_NORM</i>
  */
  m_calcErrorNorms = 1;
  m_calcErrorNorms = Context::getSolverProperty<MBool>("calcErrorNorms", m_solverId, AT_, &m_calcErrorNorms);

  if(m_calcErrorNorms) {
    /*! \property
      \page propertyPageDG DG
      \section noAnalysisNodes
      <code>MInt DgCartesianSolver::m_noAnalysisNodes </code>\n
      default = <code>2 * (m_maxPolyDeg+1)</code>\n \n
      Set the number of nodes (1D) that should be used for the error norm
      calculation.\n
      Possible values are:
      <ul>
        <li>any integer >= 2 * (m_maxPolyDeg+1)</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, POLYNOMIAL_DEGREE, ERROR_NORM</i>
    */
    m_noAnalysisNodes = 2 * (m_maxPolyDeg + 1);
    m_noAnalysisNodes = Context::getSolverProperty<MInt>("noAnalysisNodes", m_solverId, AT_, &m_noAnalysisNodes);
    if(m_noAnalysisNodes < 2 * m_maxPolyDeg) {
      const MString warning =
          "WARNING: Polynomial degree for error norms should be at least twice the maximum polynomial degree";
      m_log << warning << std::endl;
      cerr0 << warning << std::endl;
    }
    m_polyDegAnalysis = m_noAnalysisNodes - 1;

    if(m_sbpMode) {
      m_noAnalysisNodes = m_initNoNodes1D;
      m_polyDegAnalysis = m_initPolyDeg;
    }

    /*! \property
      \page propertyPageDG DG
      \section noErrorDigits
      <code>MInt DgCartesianSolver::m_noErrorDigits </code>\n
      default = <code>6</code>\n \n
      Set the number of significant digits in the error analysis output.\n
      Possible values are:
      <ul>
        <li>any integer >= 1</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, ERROR_NORM, SIGNIFICANT_DIGITS</i>
    */
    m_noErrorDigits = 6;
    m_noErrorDigits = Context::getSolverProperty<MInt>("noErrorDigits", m_solverId, AT_, &m_noErrorDigits);
  }
  /*! \property
    \page propertyPageDG DG
    \section useSponge
    <code>MBool DgCartesianSolver::m_useSponge</code>\n
    If this property is set to false, sponge is deactivated and no ouput for
    spongeEta is generated.\n
    Keywords: <i>DISCONTINUOUS_GALERKIN</i>
  */
  m_useSponge = false;
  m_useSponge = Context::getSolverProperty<MBool>("useSponge", m_solverId, AT_, &m_useSponge);

  MInt count = 0;
  MString nameCoords = "prefCoordinates_" + to_string(count);
  MString namePolyDeg = "prefPolyDeg_" + to_string(count);
  MString nameOperator = "prefOperator_" + to_string(count);
  MString nameNoNodes1D = "prefNoNodes_" + to_string(count);

  // Loop over all p-refinement patches that are defined in the properties
  while(Context::propertyExists(nameCoords, m_solverId) && Context::propertyExists(namePolyDeg, m_solverId)) {
    // Get polynomial degree and check whether it fullfills the conditions
    const MInt polyDeg = Context::getSolverProperty<MInt>(namePolyDeg, m_solverId, AT_);

    if(polyDeg < m_minPolyDeg || polyDeg > m_maxPolyDeg) {
      TERMM(1, "p-refinement: Polynomial degree of patch " + to_string(count)
                   + " doesn't fit into range defined by minimum and maximum "
                     "polynomial degree");
    }

    // Get the coordinates and check whether they fullfill the conditions
    // Check if enough coordinates are defined
    if(Context::propertyLength(nameCoords, m_solverId) != 2 * nDim) {
      TERMM(1, "p-refinement: wrong number of coordinates for patch " + to_string(count));
    }

    array<MFloat, 2 * nDim> coords;
    // Read and store patch coordinates and patch polynomial degree
    for(MInt i = 0; i < 2 * nDim; i++) {
      coords[i] = Context::getSolverProperty<MFloat>(nameCoords, m_solverId, AT_, i);
    }

    MString prefOperator = "";
    prefOperator = Context::getSolverProperty<MString>(nameOperator, m_solverId, AT_, &prefOperator);

    MInt prefNoNodes1D = polyDeg + 1;
    prefNoNodes1D = Context::getSolverProperty<MInt>(nameNoNodes1D, m_solverId, AT_, &prefNoNodes1D);

    m_prefPatchesCoords.push_back(coords);
    m_prefPatchesPolyDeg.push_back(polyDeg);
    m_prefPatchesOperators.push_back(prefOperator);
    m_prefPatchesNoNodes1D.push_back(prefNoNodes1D);

    // Increase counter and generate new names
    count++;

    nameCoords = "prefCoordinates_" + to_string(count);
    namePolyDeg = "prefPolyDeg_" + to_string(count);
    nameOperator = "prefOperator_" + to_string(count);
    nameNoNodes1D = "prefNoNodes_" + to_string(count);
  }

  /*! \property
    \page propertyPageDG DG
    \section pref
    <code>MInt DgCartesianSolver::m_pref</code>\n
    Variable to inicate whether p-Refinement is utilized or not (DEFAULT: 1)\n
    If set to 1, the additional info need to be set in the variables
    prefCoordinates_i and prefPolyDeg_i, with i indicating the
    index of the p-refinement patches
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  m_pref = 1;
  m_pref = Context::getSolverProperty<MInt>("pref", m_solverId, AT_, &m_pref);

  /*! \property
    \page propertyPageDG DG
    \section adaptivePref
    <code>MInt DgCartesianSolver::m_adaptivePref</code>\n
    Variable to choose adaptive refinement method (DEFAULT: -1)\n
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  MString adaptiveMethod = "DG_ADAPTIVE_NONE";
  m_adaptivePref = string2enum(Context::getSolverProperty<MString>("adaptivePref", m_solverId, AT_, &adaptiveMethod));

  if(hasAdaptivePref()) {
    if(g_multiSolverGrid) {
      TERMM(1, "Multi-solver does not work with adaptive p-refinement yet");
    }
    /*! \property
      \page propertyPageDG DG
      \section adaptiveInterval
      <code>MInt DgCartesianSolver::m_dgAdaptiveInterval</code>\n
      Variable to choose adaptive refinement interval (DEFAULT: -1)\n
      Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
    */
    m_adaptiveInterval = Context::getSolverProperty<MInt>("adaptiveInterval", m_solverId, AT_, &m_adaptiveInterval);

    /*! \property
    \page propertyPageDG DG
    \section adaptiveThreshold
    <code>MInt DgCartesianSolver::m_dgAdaptiveThreshold</code>\n
    Variable to choose adaptive refinement threshold (DEFAULT: -1)\n
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
    m_adaptiveThreshold =
        Context::getSolverProperty<MFloat>("adaptiveThreshold", m_solverId, AT_, &m_adaptiveThreshold);
    if(m_adaptiveInterval < 0 || m_adaptiveThreshold < 0) {
      stringstream errorMessage;
      errorMessage << "Error: Adaptive setting variables cannot be set to < 0 "
                      "(adaptiveInterval/adaptiveThreshold)"
                   << endl;
      TERMM(1, errorMessage.str());
    }
  }

  // Read in optional cut-off boundary information
  m_useCutOffBoundaries = false;
  fill(m_cutOffBoundaryConditionIds.begin(), m_cutOffBoundaryConditionIds.end(), -1);
  // Only do anything if there are cutOffBoundaryConditionIds in the properties
  if(Context::propertyExists("cutOffBoundaryConditionIds", m_solverId)) {
    // If the property exists, the default is to use the cut-off BCs
    /*! \property
      \page propertyPageDG DG
      \section useCutOffBoundaries
      <code>MInt DgCartesianSolver::m_useCutOffBoundaries </code>\n
      default = <code>true</code>\n \n
      Allow to disable cut-off boundary conditions even if they are specified.\n
      Possible values are:
      <ul>
        <li>true</li>
        <li>false</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, BOUNDARY_CONDITIONS, CUT-OFF</i>
    */
    m_useCutOffBoundaries = true;
    m_useCutOffBoundaries =
        Context::getSolverProperty<MBool>("useCutOffBoundaries", m_solverId, AT_, &m_useCutOffBoundaries);

    // If it was not explicitly disabled by the user, read in cut-off BC ids
    if(m_useCutOffBoundaries) {
      for(MInt i = 0; i < 2 * nDim; i++) {
        /*! \property
          \page propertyPageDG DG
          \section cutOffBoundaryConditionIds
          <code>MInt DgCartesianSolver::m_cutOffBoundaryConditionIds</code>\n
          default = -1</code>\n \n
          Specify boundary conditions for each spatial direction in case there
          is no geometry intersecting some cells.\n
          Possible values are:
          <ul>
            <li>a valid boundary condition id</li>
          </ul>
          Keywords: <i>DISCONTINUOUS_GALERKIN, BOUNDARY_CONDITIONS, CUT-OFF</i>
        */
        // m_cutOffBoundaryConditionIds[i]
        m_cutOffBoundaryConditionIds[i] = Context::getSolverProperty<MInt>("cutOffBoundaryConditionIds", m_solverId,
                                                                           AT_, &m_cutOffBoundaryConditionIds[i], i);
      }
    }
  }

  /*! \property
    \page propertyPageDG DG
    \section useSourceRampUp
    <code>MBool DgCartesianSolver::m_useSourceRampUp </code>\n
    default = <code>false</code>\n \n
    Select whether the external source terms should be ramped up in time. If enabled, the
    ramp-up time property (sourceRampUpTime) has to be specified as well.\n
    Keywords: <i>DG, COUPLING, SOURCE_TERM, FILTER, RAMP_UP</i>
  */
  m_useSourceRampUp = false;
  m_useSourceRampUp = Context::getSolverProperty<MBool>("useSourceRampUp", solverId(), AT_, &m_useSourceRampUp);

  if(m_useSourceRampUp) {
    /*! \property
      \page propertyPageDG DG
      \section sourceRampUpTime
      <code>MFloat DgCartesianSolver::m_sourceRampUpTime</code>\n
      default = <code>none</code>\n \n
      Set the source term ramp-up time, i.e. the duration after which the external source terms are
      fully ramped up.\n
      Keywords: <i>DG, COUPLING, SOURCE_TERM, FILTER, RAMP_UP</i>
    */
    m_sourceRampUpTime = Context::getSolverProperty<MFloat>("sourceRampUpTime", solverId(), AT_);
  }
}


/**
 * \brief Allocates main memory resources.
 *
 * \author Michael Schlottke
 * \date   October 2012
 *
 * Some of the container sizes depend directly on the number of leaf cells,
 * others are controlled by a property.
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::allocateAndInitSolverMemory() {
  TRACE();

  // Count leaf cells to determine number of elements
  MInt noElements = 0;
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    // Check if element needs to be created
    if(needElementForCell(cellId)) {
      noElements++;
    }
  }

  // Determine and set max. number of surfaces if not specified in properties
  if(m_maxNoSurfaces == -1) {
    // TODO labels:DG in some cases with partition level shift the calculated number was too small!
    m_maxNoSurfaces = calculateNeededNoSurfaces() + 100;
    m_log << "Max. no surfaces was not specified in the properties file (or "
             "set to -1). It was therefore calculated automatically for the "
             "current domain: "
          << m_maxNoSurfaces << endl;
  }

  if(isActive()) {
    MInt maxNoElements = noElements;
    MPI_Allreduce(MPI_IN_PLACE, &maxNoElements, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_,
                  "MPI_IN_PLACE", "maxNoElements");
    MInt minNoElements = noElements;
    MPI_Allreduce(MPI_IN_PLACE, &minNoElements, 1, type_traits<MInt>::mpiType(), MPI_MIN, mpiComm(), AT_,
                  "MPI_IN_PLACE", "minNoElements");
    MInt maxNoSurfaces = m_maxNoSurfaces;
    MPI_Allreduce(MPI_IN_PLACE, &maxNoSurfaces, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_,
                  "MPI_IN_PLACE", "maxNoSurfaces");

    stringstream message;
    message << "Solver #" << solverId() << " - maximum number of DG elements among ranks: " << maxNoElements
            << std::endl;
    message << "Solver #" << solverId() << " - minimum number of DG elements among ranks: " << minNoElements
            << std::endl;
    message << "Solver #" << solverId() << " - maximum number of DG surfaces among ranks: " << maxNoSurfaces
            << std::endl;
    m_log << message.str();
    cerr0 << message.str();
  }

  // Elements (allocated with maximum polynomial degree to support p-ref.)
  m_elements.maxPolyDeg(m_maxPolyDeg);
  m_elements.maxNoNodes1D(m_maxNoNodes1D);
  m_elements.noNodeVars(SysEqn::noNodeVars());
  m_elements.reset(noElements);

  // Surfaces (allocated with maximum polynomial degree to support p-ref.)
  m_surfaces.maxPolyDeg(m_maxPolyDeg);
  m_surfaces.maxNoNodes1D(m_maxNoNodes1D);
  m_surfaces.noNodeVars(SysEqn::noNodeVars());
  m_surfaces.reset(m_maxNoSurfaces);

  // Count coarse elements that smaller elements and need an h-element
  MInt noHElements = 0;
  for(MInt cellId = 0; cellId < grid().noCells(); cellId++) {
    if(needHElementForCell(cellId)) {
      noHElements++;
    }
  }

  // H-elements
  m_helements.reset(noHElements);
}

template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::writeEOC() {
  vector<MFloat> L2Error(m_sysEqn.noVars());
  vector<MFloat> LInfError(m_sysEqn.noVars());
  vector<MFloat> L2ErrLocal(m_sysEqn.noVars());
  vector<MFloat> LInfErrLocal(m_sysEqn.noVars());
  calcErrorNorms(m_finalTime, L2Error, LInfError, L2ErrLocal, LInfErrLocal);

  cout << "WRITE EOC" << endl;
  // Write output for EOC analysis
  if(m_calcErrorNorms && isMpiRoot()) {
  }
}

/**
 * \brief Main method to run this solver.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-26
 *
 * After creation of a solver instance, this should be the only method that needs
 * to be called from outside.
 *
 * In the future, when we have "true" multi solver support, this needs to be
 * rewritten properly, but for now this is where the magic takes place.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::run() {
  TRACE();

  RECORD_TIMER_START(m_timers[Timers::Run]);
  RECORD_TIMER_START(m_timers[Timers::RunInit]);

  // Set polynomial degree, init interpolation, create surfaces etc.
  initSolverObjects();

  // Apply initial conditions or load restart file
  initData();

  // Perform remaining steps needed before main loop execution
  initMainLoop();

  RECORD_TIMER_STOP(m_timers[Timers::RunInit]);

  // Run main loop
  MBool finalTimeStep = false;
  while(!finalTimeStep) {
    finalTimeStep = step();
  }

  // Clean up after the simulation is done
  cleanUp();

  RECORD_TIMER_STOP(m_timers[Timers::Run]);
}

/// \brief Initialize solver data, i.e. set initial conditions or load restart
/// file.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-19
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initData() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitData]);

  // Abort if solver not initialized
  if(!m_isInitSolver) {
    TERMM(1, "Solver was not initialized.");
  }

  if(m_restart) {
    m_log << "Restart is enabled." << endl;

    // Load restart file
    loadRestartFile();

    // Load node variable file
    if(SysEqn::noNodeVars() > 0 && !SysEqn::hasTimeDependentNodeVars()) {
      loadNodeVarsFile();
    }
  } else {
    // Apply the initial conditions
    initialCondition();

    // Reset global time step and simulation time
    m_timeStep = 0;
    m_time = m_startTime;
    m_firstTimeStep = true;
  }

  // Reset current Runge Kutta stage
  m_rkStage = 0;

  // Update all node variables, if they are used
  if(SysEqn::noNodeVars() > 0) {
    updateNodeVariables();
    // TODO labels:DG,toenhance not required at a restart when loading nodevars, but could be used to overwrite saved
    // node vars/change extension
    extendNodeVariables();
  }

  // Set status to initialized
  m_isInitData = true;

  RECORD_TIMER_STOP(m_timers[Timers::InitData]);
}


/// \brief Perform all operations that prepare the execution of the main loop.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-19
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initMainLoop() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitMainLoop]);

  // Abort if data not initialized
  if(!m_isInitData) {
    TERMM(1, "Solver data was not initialized.");
  }

  if(!m_restart) {
    // Save initial state before entering main loop
    if(m_writeInitialSolution) {
      saveSolutionFile();
    }

    // Save node variables for restarting if they are constant in time
    if(SysEqn::noNodeVars() > 0 && !SysEqn::hasTimeDependentNodeVars()) {
      if(m_writeNodeVarsFile) {
        saveNodeVarsFile();
      }
    }

    // Save initial sampling data before entering main loop
    m_pointData.save(false);
    m_surfaceData.save(false);
    m_volumeData.save(false);

    // Save inital slice data before entering main loop
    m_slice.save(false);
  }

  // Init end auto save
  m_endAutoSaveTime = -1;
  m_endAutoSaveWritten = false;
  const char* envJobEndTime = getenv("MAIA_JOB_END_TIME");

  if(envJobEndTime) {
    // Check if env variable only contains digits
    const MString jobEndTime(envJobEndTime);
    MBool onlyDigits = true;
    for(auto&& character : jobEndTime) {
      if(!isdigit(character)) {
        m_log << "Warning: the environment variable MAIA_JOB_END_TIME "
                 "included non-digit characters.\n"
              << "Saving a restart file before the compute job ends is NOT "
                 "active."
              << endl;
        onlyDigits = false;
        break;
      }
    }
    // Only activate if MAIA_JOB_END_TIME is a proper integer
    if(onlyDigits) {
      m_endAutoSaveTime = stoi(jobEndTime);
      m_endAutoSaveTime -= m_noMinutesEndAutoSave * 60;
      m_log << "Activated automatic restart file writing at " << ctime(&m_endAutoSaveTime)
            << " (in Unix time: " << m_endAutoSaveTime << ")." << endl;
    }
  }

  // Output initialization summary
  outputInitSummary();

  // Init main loop
  analyzeTimeStep(m_time, F0, F0, m_timeStep, F0);

  // Init timing for run time statistics
  m_loopTimeStart = wallTime();
  m_analyzeTimeStart = wallTime();
  m_outputTime = 0.0;
  m_noAnalyzeTimeSteps = 0;

  // Set status to initialized
  m_isInitMainLoop = true;

  RECORD_TIMER_STOP(m_timers[Timers::InitMainLoop]);

  // If multiphysics-optimized parallelization is used, prolong and start
  // tranmitting
  if(g_splitMpiComm) {
    RECORD_TIMER_START(m_timers[Timers::MainLoop]);
    RECORD_TIMER_START(m_timers[Timers::RungeKuttaStep]);
    RECORD_TIMER_START(m_timers[Timers::TimeDeriv]);
    prolongToSurfaces();
    applyForwardProjection();
    startMpiSurfaceExchange();
    RECORD_TIMER_STOP(m_timers[Timers::TimeDeriv]);
    RECORD_TIMER_STOP(m_timers[Timers::RungeKuttaStep]);
    RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
  }
}


/// \brief Advance the solution by one time step.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-19
///
/// \param[in] externalDt External time step size to use for current time step.
/// \param[in] substep Indicator if only a single Runge-Kutta substep/stage
///                    should be performed.
///
/// \return True if this is the last time step, i.e. if the main loop should be
///         terminated.
///
/// NOTE: can be removed once the transition to the unified run loop is finished
// TODO labels:DG,totest,toremove @ansgar_unified check if everything has been transfered to unified run loop, then
// remove
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::step(const MFloat externalDt, const MBool substep) {
  TRACE();
  TERMM(1, "deprecated");

  // Abort if main loop not initialized
  if(!m_isInitMainLoop) {
    TERMM(1, "Main loop was not initialized.");
  }

  RECORD_TIMER_START(m_timers[Timers::MainLoop]);

  startLoadTimer(AT_);

  // Check if this is the first Runge-Kutta stage, i.e. start of new timestep
  if(m_rkStage == 0) {
    // Increment time step
    m_timeStep++;

    // Calculate time step if it is not already set
    if(externalDt < 0.0) {
      // Calculate time step at first time step or for correct modulo
      if(m_firstTimeStep == true || (m_calcTimeStepInterval > 0 && m_timeStep % m_calcTimeStepInterval == 0)) {
        m_dt = calcTimeStep();
        m_firstTimeStep = false;
      }
    } else {
      m_dt = externalDt;
    }

    // Perform adaptive refinement if needed
    if(isAdaptationTimeStep(m_timeStep)) {
      // Cancel open MPI (receive) requests
      cancelMpiRequests();
      adaptiveRefinement(m_timeStep);
    }
  }

  // Determine if this is the last time step and reset time step if necessary
  MBool finalTimeStep = false;
  if(m_finalTime - m_time - m_dt < 1.0E-10) {
    finalTimeStep = true;
    m_dt = m_finalTime - m_time;
  } else if(m_timeStep == m_timeSteps) {
    finalTimeStep = true;
  } else if(m_timeStep == std::max(m_restartTimeStep, 0) + m_timeSteps) {
    // TODO labels:DG @ansgar_unified for coupled combustion case, step() should not be used anymore by other
    // run loops!
    finalTimeStep = true;
  }

  if(!substep) {
    // Run Runge-Kutta step
    timeStepRk(m_time, m_dt);
    // Reset current Runge-Kutta stage
    m_rkStage = 0;
  } else {
    // Perform just a single Runge-Kutta stage
    timeStepRk(m_time, m_dt, m_rkStage);
    // Update current Runge-Kutta stage
    m_rkStage = (m_rkStage + 1) % m_noRkStages;
  }

  stopLoadTimer(AT_);

  // Check if time step is completed and return if not
  if(m_rkStage != 0) {
    RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
    return finalTimeStep;
  }

  disableDlbTimers();

  // Time step completed
  // Update physical time and counters
  m_time += m_dt;
  m_noAnalyzeTimeSteps++;

  // Write out solution in intervals
  if((m_solutionInterval > 0 && m_timeStep % m_solutionInterval == 0) || (finalTimeStep && m_alwaysSaveFinalSolution)) {
    RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);

    const MFloat tOutputStart = wallTime();
    saveSolutionFile();
    const MFloat tOutputEnd = wallTime();

    // Update output time so that inner loop does not consider it
    m_outputTime += tOutputEnd - tOutputStart;

    RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);
  }

  // Write out restart file in intervals
  if((m_restartInterval > 0 && m_timeStep % m_restartInterval == 0) || (finalTimeStep && m_alwaysSaveFinalRestart)) {
    RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);

    const MFloat tOutputStart = wallTime();
    saveRestartFile();
    const MFloat tOutputEnd = wallTime();

    // Update output time so that inner loop does not consider it
    m_outputTime += tOutputEnd - tOutputStart;

    RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);
  }

  // Check regularly if job is about to end
  if(m_endAutoSaveTime != -1 && !m_endAutoSaveWritten && m_timeStep % m_endAutoSaveCheckInterval == 0) {
    // synchronize ranks (prevent deadlocks)
    MBool writeEndAutoSave = false;
    if(isMpiRoot()
       && m_endAutoSaveTime
              <= chrono::duration_cast<chrono::seconds>(chrono::system_clock::now().time_since_epoch()).count()) {
      writeEndAutoSave = true;
    }
    MPI_Bcast(&writeEndAutoSave, 1, MPI_C_BOOL, 0, mpiComm(), AT_, "writeEndAutoSave");
    // Write out restart file if job is about to end
    if(writeEndAutoSave) {
      RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);

      const MFloat tOutputStart = wallTime();
      saveRestartFile();
      time_t now =
          std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      m_log << "Finished end auto save at " << std::ctime(&now) << " (in Unix time: " << now << ")." << endl;
      time_t jobEndTime = m_endAutoSaveTime + m_noMinutesEndAutoSave * 60;
      m_log << "Job will end at " << std::ctime(&jobEndTime) << " (in Unix time: " << jobEndTime << ")." << endl;

      m_endAutoSaveWritten = true;

      const MFloat tOutputEnd = wallTime();

      // Update output time so that inner loop does not consider it
      m_outputTime += tOutputEnd - tOutputStart;

      RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);
    }
  }

  RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);
  // Produce sampling data output
  m_pointData.save(finalTimeStep);
  m_surfaceData.save(finalTimeStep);
  m_volumeData.save(finalTimeStep);
  // Produce slice data output
  m_slice.save(finalTimeStep);
  RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);

  // Analyze error norms and print info to user
  if((m_analysisInterval > 0 && m_timeStep % m_analysisInterval == 0) || (finalTimeStep)) {
    RECORD_TIMER_START(m_timers[Timers::Analysis]);

    const MFloat timeDiff = wallTime() - m_analyzeTimeStart - m_outputTime;
    const MFloat runTimeRelative = timeDiff * noDomains() / m_noAnalyzeTimeSteps / m_statGlobalNoActiveDOFs;
    const MFloat runTimeTotal = wallTime() - m_loopTimeStart;

    analyzeTimeStep(m_time, runTimeRelative, runTimeTotal, m_timeStep, m_dt);
    if(finalTimeStep && isMpiRoot()) {
      cout << endl;
      cout << "----------------------------------------"
           << "----------------------------------------" << endl;
      cout << " MAIA finished.   Final time: " << m_time << "   Time steps: " << m_timeStep << endl;
      cout << "----------------------------------------"
           << "----------------------------------------" << endl;
    }

    // Reset time + counters
    m_analyzeTimeStart = wallTime();
    m_outputTime = 0.0;
    m_noAnalyzeTimeSteps = 0;

    RECORD_TIMER_STOP(m_timers[Timers::Analysis]);
  }
  /* } else if (m_aliveInterval > 0 && m_timeStep % m_aliveInterval == 0 && isMpiRoot()) { */
  /*   // Calculate total run time for printing */
  /*   const MFloat runTimeTotal = wallTime() - m_loopTimeStart; */

  /*   // Print out processing information to user */
  /*   printf("#t/s: %8d | dt: %.4e | Sim. time: %.4e | Run time: %.4e s\n", m_timeStep, m_dt,
   * m_time, */
  /*          runTimeTotal); */
  /* } */

  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);

  enableDlbTimers();

  return finalTimeStep;
}


/// \brief Perform pre-time-step operations, e.g. set new dt if required
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::preTimeStep() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::MainLoop]);

  // Abort if main loop not initialized
  if(!m_isInitMainLoop) {
    TERMM(1, "Main loop was not initialized.");
  }

  // Check that this is the first Runge-Kutta stage, i.e. start of new timestep
  if(m_rkStage != 0) {
    TERMM(1,
          "preTimeStep should only be called at the start of a new timestep: rkStage = " + std::to_string(m_rkStage));
  }

  // Increment time step
  m_timeStep++;
  ASSERT(m_timeStep == globalTimeStep, "Error: time step inconsistent.");

  // Calculate time step (or use external dt) at first time step or for correct modulo
  if(m_firstTimeStep == true || (m_calcTimeStepInterval > 0 && m_timeStep % m_calcTimeStepInterval == 0)) {
    m_dt = calcTimeStep();
    m_firstTimeStep = false;
  }

  // Perform adaptive refinement if needed
  if(isAdaptationTimeStep(m_timeStep)) {
    // Cancel open MPI (receive) requests
    cancelMpiRequests();

    adaptiveRefinement(m_timeStep);
  }

  // Determine if this is the last time step and reset time step size if necessary
  m_finalTimeStep = false;
  if(m_finalTime - m_time - m_dt < 1.0E-10) {
    m_finalTimeStep = true;
    m_dt = m_finalTime - m_time;
  } else if(g_multiSolverGrid && m_timeStep == std::max(m_restartTimeStep, 0) + m_timeSteps) {
    m_finalTimeStep = true;
  } else if(!g_multiSolverGrid && m_timeStep == m_timeSteps) {
    // TODO labels:DG change default time step behaviour of DG solver?
    m_finalTimeStep = true;
  }

  // Reset external (coupling) sources at beginning of a new time step
  resetExternalSources();

  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
}


/// \brief Perform one Runge-Kutta step/stage
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::solutionStep() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::MainLoop]);

  // Perform just a single Runge-Kutta step/stage
  timeStepRk(m_time, m_dt, m_rkStage);
  // Update current Runge-Kutta stage
  m_rkStage = (m_rkStage + 1) % m_noRkStages;

  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
  // Return if time step is completed
  return (m_rkStage == 0);
}


/// \brief Perform post-time-step operations, e.g. advance time, error analysis
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::postTimeStep() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::MainLoop]);

  // Time step completed
  // Update physical time and counters
  m_time += m_dt;
  m_noAnalyzeTimeSteps++;

  // Analyze error norms and print info to user
  if((m_analysisInterval > 0 && m_timeStep % m_analysisInterval == 0) || (m_finalTimeStep)) {
    RECORD_TIMER_START(m_timers[Timers::Analysis]);

    const MFloat timeDiff = wallTime() - m_analyzeTimeStart - m_outputTime;
    const MFloat runTimeRelative = timeDiff * noDomains() / m_noAnalyzeTimeSteps / m_statGlobalNoActiveDOFs;
    const MFloat runTimeTotal = wallTime() - m_loopTimeStart;

    analyzeTimeStep(m_time, runTimeRelative, runTimeTotal, m_timeStep, m_dt);
    if(m_finalTimeStep && isMpiRoot()) {
      cout << endl;
      cout << "----------------------------------------"
           << "----------------------------------------" << endl;
      cout << " MAIA finished.   Final time: " << m_time << "   Time steps: " << m_timeStep << endl;
      cout << "----------------------------------------"
           << "----------------------------------------" << endl;
    }

    // Reset time + counters
    m_analyzeTimeStart = wallTime();
    m_outputTime = 0.0;
    m_noAnalyzeTimeSteps = 0;

    RECORD_TIMER_STOP(m_timers[Timers::Analysis]);
  }
  // TODO labels:DG,toremove moved to run_unified() but no output of sim-time/dt
  /* } else if (m_aliveInterval > 0 && m_timeStep % m_aliveInterval == 0 && isMpiRoot()) { */
  /*   // Calculate total run time for printing */
  /*   const MFloat runTimeTotal = wallTime() - m_loopTimeStart; */

  /*   // Print out processing information to user */
  /*   printf("#t/s: %8d | dt: %.4e | Sim. time: %.4e | Run time: %.4e s\n", m_timeStep, m_dt,
   * m_time, */
  /*          runTimeTotal); */
  /* } */

  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
}


/// \brief Return true if the solver wants to write a restart file
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::prepareRestart(MBool writeRestart, MBool& NotUsed(writeGridRestart)) {
  TRACE();

  if(!isActive()) return false;

  if((m_restartInterval > 0 && m_timeStep % m_restartInterval == 0) || (m_finalTimeStep && m_alwaysSaveFinalRestart)) {
    return true;
  }

  return writeRestart;
}


/// \brief Write a restart file
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::writeRestartFile(const MBool writeRestart,
                                                       const MBool NotUsed(writeBackup),
                                                       const MString NotUsed(gridFileName),
                                                       MInt* NotUsed(recalcIdTree)) {
  TRACE();
  CHECK_TIMERS_IO();
  RECORD_TIMER_START(m_timers[Timers::MainLoop]);

  // Write out restart file
  if(writeRestart) {
    RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);

    const MFloat tOutputStart = wallTime();
    saveRestartFile();
    const MFloat tOutputEnd = wallTime();

    // Update output time so that inner loop does not consider it
    m_outputTime += tOutputEnd - tOutputStart;

    RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);
  }

  // TODO labels:DG integrate this in the unified run loop/gridcontroller for all solvers!
  /*   // Check regularly if job is about to end */
  /*   if (m_endAutoSaveTime != -1 && !m_endAutoSaveWritten */
  /*       && m_timeStep % m_endAutoSaveCheckInterval == 0) { */
  /*     // synchronize ranks (prevent deadlocks) */
  /*     MBool writeEndAutoSave = false; */
  /*     if (isMpiRoot() */
  /*         && m_endAutoSaveTime */
  /*                 <= chrono::duration_cast<chrono::seconds>( */
  /*                       chrono::system_clock::now().time_since_epoch()) */
  /*                       .count()) { */
  /*       writeEndAutoSave = true; */
  /*     } */
  /*     MPI_Bcast(&writeEndAutoSave, 1, MPI_C_BOOL, 0, mpiComm(), AT_, "writeEndAutoSave" ); */
  /*     // Write out restart file if job is about to end */
  /*     if (writeEndAutoSave) { */
  /*       RECORD_TIMER_START(m_timers[Timers::MainLoopIO]); */

  /*       const MFloat tOutputStart = wallTime(); */
  /*       saveRestartFile(); */
  /*       time_t now = std::chrono::duration_cast<std::chrono::seconds>( */
  /*                         std::chrono::system_clock::now().time_since_epoch()) */
  /*                         .count(); */
  /*       m_log << "Finished end auto save at " << std::ctime(&now) */
  /*               << " (in Unix time: " << now << ")." << endl; */
  /*       time_t jobEndTime = m_endAutoSaveTime + m_noMinutesEndAutoSave * 60; */
  /*       m_log << "Job will end at " << std::ctime(&jobEndTime) */
  /*               << " (in Unix time: " << jobEndTime << ")." << endl; */

  /*       m_endAutoSaveWritten = true; */

  /*       const MFloat tOutputEnd = wallTime(); */

  /*       // Update output time so that inner loop does not consider it */
  /*       m_outputTime += tOutputEnd - tOutputStart; */

  /*       RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]); */
  /*     } */
  /*   } */
  /* } */

  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
}


/// \brief Save the solver solution, i.e. write solution files and sampling/slice output
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveSolverSolution(const MBool forceOutput, const MBool finalTimeStep) {
  TRACE();
  CHECK_TIMERS_IO();
  RECORD_TIMER_START(m_timers[Timers::MainLoop]);
  RECORD_TIMER_START(m_timers[Timers::MainLoopIO]);

  // Write out solution in intervals, at the final time step, or if forced
  if((m_solutionInterval > 0 && m_timeStep % m_solutionInterval == 0)
     || ((m_finalTimeStep || finalTimeStep) && m_alwaysSaveFinalSolution) || forceOutput) {
    const MFloat tOutputStart = wallTime();
    saveSolutionFile();
    const MFloat tOutputEnd = wallTime();

    // Update output time so that inner loop does not consider it
    m_outputTime += tOutputEnd - tOutputStart;
  }

  // Produce sampling data output
  m_pointData.save(finalTimeStep);
  m_surfaceData.save(finalTimeStep);
  m_volumeData.save(finalTimeStep);
  // Produce slice data output
  m_slice.save(finalTimeStep);

  RECORD_TIMER_STOP(m_timers[Timers::MainLoopIO]);
  RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
}


/// \brief Initialize the solver
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initSolver() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Run]);

  // Nothing to be done if solver is not active
  if(!isActive()) {
    // TODO labels:DG inactive ranks need to call createGridSlice as well
    TERMM_IF_COND(m_slice.enabled(), "Fixme: DG slices with inactive ranks not supported.");
    return;
  }

  RECORD_TIMER_START(m_timers[Timers::RunInit]);

  // Set polynomial degree, init interpolation, create surfaces etc.
  initSolverObjects();

  // Apply initial conditions or load restart file
  initData();

  // Note: moved to finalizeInitSolver, coupler need to be initialized first
  // Perform remaining steps needed before main loop execution
  // initMainLoop();

  RECORD_TIMER_STOP(m_timers[Timers::RunInit]);
}


/// \brief Finalization of the solver initialization
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::finalizeInitSolver() {
  TRACE();

  // Nothing to be done if solver is not active
  if(!isActive()) return;

  RECORD_TIMER_START(m_timers[Timers::RunInit]);
  // Perform remaining steps needed before main loop execution
  initMainLoop();
  RECORD_TIMER_STOP(m_timers[Timers::RunInit]);
}


/** \brief Initializes the solver. This must be called before using any of the
 *         discretization methods, and should be followed (after the simulation
 *         run) by a call to cleanUp().
 *
 * \author Michael Schlottke
 * \date   October 2012
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initSolverObjects() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitSolverObjects]);

  // NOTE: if anything is added/changed here, make the same changes in balancePre()/balancePost()

  // Init grid-to-grid map
  initGridMap();

  // Init elements
  initElements();

  // Init h-elements
  initHElements();

  // Set the polynomial degree in the elements
  initPolynomialDegree();

  // Calculate all necessary interpolation values
  initInterpolation();

  // Calculate the coordinates of the integration nodes in the elements
  initNodeCoordinates();

  // Calculate the inverse Jacobian for the mapping to the reference element
  initJacobian();

  // Check cell properties
  checkCellProperties();

  // Create all surfaces (except between internal cells and halo cells)
  initSurfaces();

  // Initialize all necessary routines & data arrays for MPI communication
  if(hasMpiExchange()) {
    initMpiExchange();
  }

  // Check whether sponge is activated
  if(useSponge()) {
    m_log << "Initializing sponge... ";
    sponge().init(m_maxPolyDeg, &grid(), &m_elements, &m_surfaces, &m_boundaryConditions, &m_sysEqn, mpiComm());
    m_log << "done" << endl;
  }

  // For all halo cells send a list of polynomial degrees to neighboring domain
  // (p-refinement)
  if(hasMpiExchange() && hasPref()) {
    exchangeMpiSurfacePolyDeg();
  }

  // Calculate all necessary mortar projections
  initMortarProjections();

  // Load points at which the states should be written out
  m_pointData.init();
  // Load surface points on which the states should be written out
  m_surfaceData.init();
  // Load volume points on which the states should be written out
  m_volumeData.init();

  // Load coordinates of slice intercept
  m_slice.init();

  // Initialize statistics about the simulation itself and write it to m_log
  initSimulationStatistics();

  // Set status to initialized
  m_isInitSolver = true;

  RECORD_TIMER_STOP(m_timers[Timers::InitSolverObjects]);
}


/// \brief Determine grid-to-grid mapping.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-11-29
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initGridMap() {
  TRACE();

  /*! \property
    \page propertyPageDG DG
    \section donorGridFileName
    <code>MInt DgCartesianSolver::donorGridFileName </code>\n
    default = <code>n/a</code>\n \n
    Name of donor grid (on which sources terms have been computed).
    Possible values are:
    <ul>
      <li>any valid file name</li>
    </ul>
    Keywords: <i>COUPLING, GRIDMAP, PARALLEL, FILENAME</i>
   */

  // Return fast if no donor grid is specified
  if(!Context::propertyExists("donorGridFileName", m_solverId)) {
    return;
  }

  /*! \property
    \page propertyPageDG DG
    \section loadGridMap
    <code>MInt DgCartesianSolver::loadGridMap </code>\n
    default = <code>0</code>\n \n
    Grid map will be loaded if specified to 1.
    Possible values are:
    <ul>
      <li>1 : to load grid-map</li>
    </ul>
    Keywords: <i>COUPLING, GRIDMAP, PARALLEL</i>
   */

  // Return if no grid map should be loaded
  if(!(Context::getSolverProperty<MBool>("loadGridMap", m_solverId, AT_))) {
    return;
  }

  /*! \property
    \page propertyPageDG DG
    \section gridMapFileName
    <code>MString DgCartesianSolver::gridMapFileName </code>\n
    default = <code>n/a</code>\n \n
    Name of the grid-map (mapping btn. donor and recipient) file that is generated.\n
    Possible values are:
    <ul>
      <li>any valid file name</li>
    </ul>
    Keywords: <i>COUPLING, GRIDMAP, PARALLEL, FILENAME</i>
   */

  // Get grid map file name if specified
  MString gridMapFileName = outputDir() + "gridmap" + ParallelIo::fileExt();
  gridMapFileName = Context::getSolverProperty<MString>("gridMapFileName", m_solverId, AT_, &gridMapFileName);
  // Load grid map file (grid map was already generated in CartesianGrid)
  loadGridMap(gridMapFileName);

  /*! \property
    \page propertyPageDG DG
    \section checkGridMap
    <code>MBool check</code>\n
    default = <code>1</code>\n \n
    Enable/disable sanity checks for the grid-map.\n
    Possible values are:
    <ul>
    <li>0: disabled</li>
    <li>1: enabled</li>
    </ul>
    Keywords: <i>COUPLING, GRIDMAP, PARALLEL, FILENAME</i>
  */
  // Check grid map if desired
  MBool check = true;
  check = Context::getSolverProperty<MBool>("checkGridMap", m_solverId, AT_, &check);
  if(check) {
    // Get name of donor grid file
    const MString donorGridFileName = Context::getSolverProperty<MString>("donorGridFileName", m_solverId, AT_);

    checkGridMap(donorGridFileName);
  }
}


/// \brief Load previously created grid map.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-17
///
/// \param[in] gridMapFileName File name of grid map.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::loadGridMap(const MString& gridMapFileName) {
  TRACE();

  m_log << "Loading grip map from " << gridMapFileName << "..." << endl;

  // Read grid map
  const MInt noCells = grid().noInternalCells();
  MIntScratchSpace gridMap(noCells, AT_, "gridMap");
  MIntScratchSpace gridMapNoOffspring(noCells, AT_, "gridMapNoOffspring");
  using namespace parallel_io;
  ParallelIo gridMapFile(gridMapFileName, PIO_READ, mpiComm());
  gridMapFile.setOffset(noCells, grid().domainOffset(domainId()));
  gridMapFile.readArray(&gridMap[0], "cellIds");
  gridMapFile.readArray(&gridMapNoOffspring[0], "noOffspring");

  // Temporary storage for the number of offspring cells mapped to target cells
  MIntScratchSpace noOffspring(noCells, AT_, "noOffspring");

  // Clear previous grid map offsets
  m_gridMapOffsets.clear();

  // Determine grid map offsets
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    // Skip if mapped donor cell id is -1, i.e., if there is no donor cell
    if(gridMap[cellId] == -1) {
      continue;
    }

    // Store first local target cell id and first mapped donor cell global id
    const MInt firstTargetCellId = cellId;
    const MInt firstDonorGlobalId = gridMap[cellId];

    // Store the number of donor offspring cells for the first considered target
    // cell
    noOffspring[0] = gridMapNoOffspring[cellId];

    // Find consecutive mapped ids (account for child cells on the donor grid),
    // i.e., consecutive cells on the target grid that have one or more donor
    // cells
    MInt noDonorCells = 1 + gridMapNoOffspring[cellId];
    while(cellId + 1 < noCells) {
      if(gridMap[cellId + 1] == firstDonorGlobalId + noDonorCells) {
        // If next mapped target cell id is contiguous, increase donor size and
        // proceed with next
        cellId++;
        noDonorCells += 1 + gridMapNoOffspring[cellId];
        noOffspring[cellId - firstTargetCellId] = gridMapNoOffspring[cellId];
      } else if(gridMap[cellId + 1] == firstDonorGlobalId + noDonorCells - 1) {
        // One-to-many mapping: the next target cell is also mapped to the
        // previous donor cell.
        cellId++;
        // Set number of offspring to -1 to indicate one-to-many mapping
        noOffspring[cellId - firstTargetCellId] = -1;
      } else {
        // Otherwise exit while loop
        break;
      }
    }

    // Determine number of contiguous mapped target cells
    const MInt noTargetCells = cellId - firstTargetCellId + 1;

    // Create & store grid map offset
    m_gridMapOffsets.push_back({firstTargetCellId, noTargetCells, firstDonorGlobalId, noDonorCells,
                                std::vector<MInt>(noOffspring.begin(), noOffspring.end())});
  }

  m_log << "Found " << m_gridMapOffsets.size() << " contiguous mapped region(s):" << endl;
  for(size_t i = 0; i < m_gridMapOffsets.size(); i++) {
    const auto& o = m_gridMapOffsets[i];
    m_log << "- offset: " << i << "; firstTargetCellId: " << o.m_firstTargetCellId
          << "; globalId: " << grid().tree().globalId(o.m_firstTargetCellId) << "; noDonorCells: " << o.m_noDonorCells
          << "; noTargetCells: " << o.m_noTargetCells << "; firstDonorGlobalId: " << o.m_firstDonorGlobalId << endl;
  }

  // Determine maximum amount of offsets across all domains
  m_maxNoGridMapOffsets = m_gridMapOffsets.size();
  MPI_Allreduce(MPI_IN_PLACE, &m_maxNoGridMapOffsets, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_,
                "MPI_IN_PLACE", "m_maxNoGridMapOffsets");
  m_log << "Maximum number of grid map offsets on all domains: " << m_maxNoGridMapOffsets << endl;
}


/// \brief Perform some sanity checks on loaded grid map.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-11-04
///
/// \param[in] donorGridFileName Name of donor grid file.
///
/// Note: tests should be sorted by increasing test time, i.e. simpler tests
/// should go first.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::checkGridMap(const MString& donorGridFileName) {
  TRACE();

  //------------------------------
  //------------------------------
  //------------------------------
  if(domainId() == 0) {
    cerr << "Warning: DgCartesianSolver::checkGridMap deactivated! (fix coordinate comparison)" << endl;
  }
  return;
  //------------------------------
  //------------------------------
  //------------------------------

  m_log << "Checking grip map for donor grid " << donorGridFileName << "...";
  const MFloat startTime = wallTime();

  // Determine the local number of donor cells and the maximum number of donor
  // cells across all grid map offsets on this domain
  MInt noDonorCellsLocal = 0;
  MInt gridMapMaxNoDonorCells = 0;
  for(auto&& offset : m_gridMapOffsets) {
    noDonorCellsLocal += offset.m_noDonorCells;
    gridMapMaxNoDonorCells = max(gridMapMaxNoDonorCells, offset.m_noDonorCells);
  }


  //////////////////////////////////////////////////////////////////////////////
  // Check 1: total number of donor cells matches number of cells in grid file
  //////////////////////////////////////////////////////////////////////////////
  using namespace maia::parallel_io;
  ParallelIo file(donorGridFileName, PIO_READ, mpiComm());
  MInt noDonorCellsFile;
  // file.readScalar(&noDonorCellsFile, "noCells");
  file.getAttribute(&noDonorCellsFile, "noCells");
  MInt noDonorCellsGlobal;
  MPI_Allreduce(&noDonorCellsLocal, &noDonorCellsGlobal, 1, type_traits<MInt>::mpiType(), MPI_SUM, mpiComm(), AT_,
                "noDonorCellsLocal", "noDonorCellsGlobal");
  if(noDonorCellsFile != noDonorCellsGlobal) {
    TERMM(1, "Number of mapped donor cells does not match number of cells in "
             "donor grid file: mapped: "
                 + to_string(noDonorCellsGlobal) + "; grid: " + to_string(noDonorCellsFile));
  }


  //////////////////////////////////////////////////////////////////////////////
  // Check 2: coordinates of donor cells match those of target cells
  //////////////////////////////////////////////////////////////////////////////
  // Define epsilon according to CartesianGrid constructor
  const MFloat eps = 1.0 / FPOW2(30) * grid().lengthLevel0();

  // Successively read coordinates of donor grid and check them against target
  // grid coordinates
  MFloatScratchSpace coordinates(max(gridMapMaxNoDonorCells, 1), AT_, "coordinates");
  const array<MString, 3> dirs = {{"x", "y", "z"}};
  for(MInt offsetId = 0; offsetId < m_maxNoGridMapOffsets; offsetId++) {
    // Store whether this is a valid iteration or if this is just done to make
    // sure that reading from the grid file in parallel works correctly
    const MBool valid = (static_cast<size_t>(offsetId) < m_gridMapOffsets.size());

    // Set offset dependent on whether this is a valid iteration
    const MInt noDonorCells = valid ? m_gridMapOffsets[offsetId].m_noDonorCells : 0;
    const MInt firstDonorGlobalId = valid ? m_gridMapOffsets[offsetId].m_firstDonorGlobalId : 0;
    file.setOffset(noDonorCells, firstDonorGlobalId);

    // Read and compare coordinates
    for(MInt dir = 0; dir < nDim; dir++) {
      // Read from file
      file.readArray(&coordinates[0], "coordinates_" + to_string(dir));

      // Skip comparison if not a valid iteration
      if(!valid) {
        continue;
      }

      const MInt noTargetCells = m_gridMapOffsets[offsetId].m_noTargetCells;
      MInt donorCellId = 0;
      // Initialize donorLength to bogus value s.t. it always triggers an abort
      // if this value is (erroneously) used before a valid length is set
      MFloat donorLength = -numeric_limits<MFloat>::infinity();

      // Compare coordinates
      for(MInt i = 0; i < noTargetCells; i++) {
        const MInt cellId = m_gridMapOffsets[offsetId].m_firstTargetCellId + i;
        const MFloat targetCoordinate = grid().tree().coordinate(cellId, dir);

        // Check if this element belongs to a one-to-many mapping
        if(m_gridMapOffsets[offsetId].m_noOffspring[i] == -1) {
          const MFloat lastDonorCoordinate = coordinates[donorCellId - 1];
          // Check if target cell is contained in the donor cell, which is the
          // same as the last considered and matching donor cell, for which the
          // cell length  is stored in 'donorLength'.
          if(targetCoordinate < lastDonorCoordinate - 0.5 * donorLength
             || targetCoordinate > lastDonorCoordinate + 0.5 * donorLength) {
            TERMM(1, "Target cell of one-to-many mapping not contained in "
                     "donor cell.");
          }
          // Continue with next target cell
          continue;
        }

        // Compare coordinates of matching target and donor cell
        const MFloat donorCoordinate = coordinates[donorCellId];
        if(!approx(targetCoordinate, donorCoordinate, eps)) {
          const MInt globalId = grid().tree().globalId(cellId);
          TERMM(1, dirs[dir] + "-coordinates for target cell " + to_string(cellId) + " (globalId: "
                       + to_string(globalId) + ") and donor cell " + to_string(firstDonorGlobalId + donorCellId)
                       + " do not match: target: " + to_string(targetCoordinate)
                       + "; donor: " + to_string(donorCoordinate));
        }
        donorCellId++;

        const MFloat targetLength = grid().cellLengthAtCell(cellId);
        // Store cell length of matching donor cell for one-to-many check
        donorLength = targetLength;

        // Check if donor offspring cells are contained in target cell
        for(MInt offspringId = 0; offspringId < m_gridMapOffsets[offsetId].m_noOffspring[i]; offspringId++) {
          const MFloat offspringCoordinate = coordinates[donorCellId + offspringId];
          if(offspringCoordinate < targetCoordinate - 0.5 * targetLength
             || offspringCoordinate > targetCoordinate + 0.5 * targetLength) {
            TERMM(1, "Donor offspring not contained in target cell.");
          }
        }
        // Increase by the number of donor child cells for this target cell
        donorCellId += m_gridMapOffsets[offsetId].m_noOffspring[i];
      }
    }
  }

  // Show duration to allow estimate about how costly this check was
  m_log << "OK (duration: " << (wallTime() - startTime) << " s)" << endl;
}


/// \brief Initialize all elements by iterating over all cells and creating an
///        element for each internal leaf cell.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-13
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initElements() {
  TRACE();

  m_log << "Initializing elements... ";
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    // Check if element needs to be created
    if(needElementForCell(cellId)) {
      createElement(cellId);
    }
  }
  m_log << "done" << endl;
}

/// \brief Initialize all helements by iterating over all cells and creating an
///        element for each coarse cell in a h-refined grid.
///
/// \author Sven Berger
/// \date March 2015
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initHElements() {
  TRACE();

  m_log << "Initializing h-refined elements... ";

  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    // Check if h-element needs to be created
    if(needHElementForCell(cellId)) {
      createHElement(cellId);
    }
  }

  m_log << "done" << endl;
}


/// \brief Return true if element is needed for cell, false otherwise.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-21
///
/// \param[in] cellId Id of cell that should be checked.
///
/// \return True if element is needed, false otherwise.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::needElementForCell(const MInt cellId) {
  // TRACE();

  MBool needElement = true;

  // Do not create element if...
  // ... cell has child cells or
  if(grid().tree().hasChildren(cellId)
     // ... cell center is outside the geometry
     || !pointIsInside(&grid().tree().coordinate(cellId, 0))) {
    needElement = false;
  }
  return needElement;
}


/// \brief Return true if h-element is needed for cell, false otherwise.
///
/// \author Sven Berger
/// \date March 2015
///
/// \param[in] cellId Id of cell that should be checked.
///
/// \return True if h-element is needed, false otherwise.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::needHElementForCell(const MInt cellId) {
  // TRACE();

  // Only leaf cells have elements
  if(grid().tree().hasChildren(cellId)) {
    return false;
  }

  // Check neighbors in every direction. If neighbor has children, then neighbor
  // is h-refined and coarse element requires an h-element
  for(MInt dir = 0; dir < 2 * nDim; dir++) {
    if(grid().tree().hasNeighbor(cellId, dir)) {
      const MInt neighborId = grid().tree().neighbor(cellId, dir);
      if(grid().tree().hasChildren(neighborId)) {
        return true;
      }
    }
  }

  return false;
}


/// \brief Determine the number of surfaces that need to be created on this
/// domain.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-11-27
///
/// \return The number of needed DG surfaces.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::calculateNeededNoSurfaces() {
  TRACE();

  // Loop over all internal cells
  MInt noSurfaces = 0;
  for(MInt cellId = 0; cellId < grid().noInternalCells(); cellId++) {
    // Skip cell if it is not a leaf cell (i.e. if it has child cells)
    if(grid().tree().hasChildren(cellId)) {
      continue;
    }

    // Loop over all directions to check if surfaces need to be created
    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      if(!grid().tree().hasNeighbor(cellId, dir)) {
        // If it does not have a neighbor, cell is boundary cell or neighbor is
        // coarse, thus one surface is needed
        noSurfaces++;
      } else {
        // Otherwise check if surface needs to be created
        const MInt neighborId = grid().tree().neighbor(cellId, dir);
        const MBool neighborIsHalo = (neighborId >= grid().noInternalCells());

        if(grid().tree().hasChildren(neighborId)) {
          // Neighbor is refined
          if(neighborIsHalo) {
            // Only create surface towards refined cells if neighbor is halo
            // cell (in this case 2 ^ (nDim - 1) surfaces are needed). This
            // slightly overestimates the actual number of needed surfaces in
            // case no all refined cells exist
            noSurfaces += IPOW2(nDim - 1);
          }
        } else {
          // Neighbor is at same level
          if(dir % 2 == 1 || neighborIsHalo) {
            // Count surfaces only in positive direction or if neighbor is halo
            // cell
            noSurfaces++;
          }
        }
      }
    }
  }

  return noSurfaces;
}


/// \brief Check if point is inside the computational domain.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-04-23
///
/// \param[in] coordinates The coordinates of the point to check.
///
/// \return If the point is inside the computational domain (i.e. not inside
///         "non-fluid" areas), return true. Return false otherwise.
///
/// FIXME labels:DG This has been shamelessly copied from GridgenPar. This should be
/// implemented in the grid or geometry classes and used from there. However,
/// currently CartesianGrid has a mess of pointIsInside methods and it is not
/// clear which one to use.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::pointIsInside(const MFloat* const coordinates) {
  // TRACE();

  // Get bounding box of grid
  const MFloat* const box = geometry().boundingBox();

  // Shoot rays in each direction
  MBool isInside = true;
  for(MInt dir = 0; dir < nDim; dir++) {
    // Cast ray from the point in question in the -ve 'dir'-direction
    array<MFloat, 2 * nDim> target;
    for(MInt i = 0; i < nDim; i++) {
      target[i] = coordinates[i];
      target[i + nDim] = coordinates[i];
    }
    // This line ensures that the second point is well outside the bounding box
    target[dir] = box[dir] - (target[dir] - box[dir]);

    // Get list of intersecting geometry elements
    std::vector<MInt> nodeList;
    geometry().getLineIntersectionElements(&target[0], nodeList);

    // Point is inside iff the number of cuts is uneven
    if(nodeList.size() % 2 == 0) {
      isInside = false;
      break;
    }
  }

  return isInside;
}


/// \brief Create element for cell with id `cellId`.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-21
///
/// \param[in] cellId Cell id of cell for which element should be created.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::createElement(const MInt cellId) {
  // TRACE();

  // Determine element id and create element in collector
  const MInt elementId = m_elements.size();
  m_elements.append();

  // Set cell id for forward references to cells
  m_elements.cellId(elementId) = cellId;

  // Reset surface ids
  for(MInt dir = 0; dir < 2 * nDim; dir++) {
    m_elements.surfaceIds(elementId, dir) = -1;
  }
}

/// \brief Create h-element for cell with id `cellId`.
///
/// \author Sven Berger
/// \date March 2015
///
/// \param[in] cellId Cell id of cell for which h-element should be created.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::createHElement(const MInt cellId) {
  // TRACE();

  // Determine h-element id and create h-element in collector
  const MInt hElementId = m_helements.size();
  m_helements.append();

  // Set elementId so that it can be used in loops over the h-element collector
  m_helements.elementId(hElementId) = m_elements.getElementByCellId(cellId);

  // Reset h-refined surface ids
  for(MInt dir = 0; dir < 2 * nDim; dir++) {
    for(MInt pos = 0; pos < 2 * (nDim - 1); pos++) {
      m_helements.hrefSurfaceIds(hElementId, dir, pos) = -1;
    }
  }
}


/**
 * \brief Calculate and set initial polynomial degree and number of nodes in all elements
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-03
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initPolynomialDegree() {
  TRACE();

  m_log << "Initializing polynomial degree and number of nodes in the elements... ";

  if(!hasPref()) {
    fill_n(&m_elements.polyDeg(0), m_elements.size(), m_initPolyDeg);
    fill_n(&m_elements.noNodes1D(0), m_elements.size(), m_initNoNodes1D);
  } else {
    // If p-refinement is used, select method based on property settings
    initPrefinement();
  }

  // Once the polynomial degree is set, calculate the total data size and number
  // of nodes on this domain.
  // The internal data size counts the total number of MFloats allocated for
  // the conservative variables
  m_internalDataSize = m_elements.size() * m_elements.maxNoNodesXD() * SysEqn::noVars();

  // The total number of nodes counts the number of nodes that are in use, i.e.
  // not counting unused but allocated nodes for elements with a polynomial
  // degree less than the maximum
  m_noTotalNodesXD = 0;
  for(MInt i = 0; i < m_elements.size(); i++) {
    m_noTotalNodesXD += m_elements.noNodesXD(i);
  }

  m_log << "done" << endl;
}


/// \brief Set polynomial degree for static p-refinement case.
///
/// \author Sven Berger, Bjoern Peeters
/// \date   November 2014, November 2015
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initPrefinement() {
  TRACE();

  m_log << "Static p-refinement is enabled." << endl;

  MIntScratchSpace changed(m_elements.size(), AT_, "changed");
  fill(changed.begin(), changed.end(), 0);
  const MInt noElements = m_elements.size();

  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);

    // Static p-refinement is done through patches that specify a region and a
    // polynomial degree. Elements that are in two or more patches get the
    // polynomial degree of the last matching patch.

    // Initialize polynomial degree with initPolyDeg, which is used if no
    // matching patch is found
    MInt elemPolyDeg = m_initPolyDeg;
    MInt elemNoNodes1D = m_initNoNodes1D;

    // Iterate over all defined patches to check if they contain the current
    // cell
    for(std::vector<MFloat>::size_type patchId = 0; patchId < m_prefPatchesPolyDeg.size(); patchId++) {
      MBool inPatch = true;

      for(MInt i = 0; i < nDim; i++) {
        const MFloat cellCoord = grid().tree().coordinate(cellId, i);

        if(cellCoord < m_prefPatchesCoords[patchId][i]) {
          inPatch = false;
        }
        if(cellCoord > m_prefPatchesCoords[patchId][i + nDim]) {
          inPatch = false;
        }
      }

      // Set the new polynomial degree if the element lies inside the patch
      if(inPatch) {
        elemPolyDeg = m_prefPatchesPolyDeg[patchId];
        elemNoNodes1D = m_prefPatchesNoNodes1D[patchId];

        changed[elementId]++;
      }
    }
    m_elements.polyDeg(elementId) = elemPolyDeg;
    m_elements.noNodes1D(elementId) = elemNoNodes1D;
  }

  m_log << "p-refinement"
        << " changed the polynomial degree "
        << "on " << count_if(changed.begin(), changed.end(), [](MInt i) { return i > 0; }) << " elements in total for "
        << accumulate(changed.begin(), changed.end(), 0) << " times." << endl;
}


/**
 * \brief Calculates necessary coefficients for interpolation and stores them
 *        once for the whole solver.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initInterpolation() {
  TRACE();

  m_log << "Initializing interpolation data... ";

  // Interpolation
  m_interpolation.clear();
  m_interpolation.resize(m_maxPolyDeg + 1);
  for(MInt polyDeg = m_minPolyDeg; polyDeg <= m_maxPolyDeg; polyDeg++) {
    m_interpolation[polyDeg].clear();
    m_interpolation[polyDeg].resize(m_maxNoNodes1D + 1);
  }

  // Convert integers to enums
  auto polyType = static_cast<DgPolynomialType>(m_dgPolynomialType);
  auto intMethod = static_cast<DgIntegrationMethod>(m_dgIntegrationMethod);

  // Create DgInterpolation object for each possible occurring polynomial
  // degree AND number of nodes (trivial for DG, necessary for SBP)

  for(MInt i = m_minPolyDeg; i <= m_maxPolyDeg; i++) {
    MInt prefIndex = -1;
    for(std::vector<MFloat>::size_type j = 0; j < m_prefPatchesPolyDeg.size(); j++) {
      if(i == (MInt)m_prefPatchesPolyDeg[j]) {
        prefIndex = j;
        break;
      }
    }

    MString refOperator = m_sbpOperator;
    MInt refNoNodes1D = m_sbpMode ? m_initNoNodes1D : (i + 1);
    if(prefIndex != -1) {
      refOperator = m_prefPatchesOperators[prefIndex];
      refNoNodes1D = m_prefPatchesNoNodes1D[prefIndex];
    }
    m_interpolation[i][refNoNodes1D].init(i, polyType, refNoNodes1D, intMethod, m_sbpMode, refOperator);
  }

  // Local and global volume
  // Calculate the local and global volume for use in error normalizations
  m_globalVolume = F0;
  m_localVolume = F0;
  for(MInt elementId = 0; elementId < m_elements.size(); elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    m_localVolume += pow(grid().cellLengthAtCell(cellId), nDim);
  }

  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::MPI]);
  RECORD_TIMER_START(m_timers[Timers::MPIComm]);
  MPI_Allreduce(&m_localVolume, &m_globalVolume, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "m_localVolume",
                "m_globalVolume");
  RECORD_TIMER_STOP(m_timers[Timers::MPIComm]);
  RECORD_TIMER_STOP(m_timers[Timers::MPI]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);

  // Init interpolation for error analysis
  if(m_calcErrorNorms) {
    // Init interpolation object for error analysis
    m_interpAnalysis.init(m_polyDegAnalysis, polyType, m_noAnalysisNodes, intMethod, m_sbpMode, m_sbpOperator);

    // Create volume weights for error analysis
    const MInt noNodesAnalysis1D = m_noAnalysisNodes;
    const MInt noNodesAnalysis1D3 = (nDim == 3) ? noNodesAnalysis1D : 1;
    m_wVolumeAnalysis.resize(noNodesAnalysis1D, noNodesAnalysis1D, noNodesAnalysis1D3);
    const MFloatVector& wInt = m_interpAnalysis.m_wInt;
    for(MInt i = 0; i < noNodesAnalysis1D; i++) {
      for(MInt j = 0; j < noNodesAnalysis1D; j++) {
        for(MInt k = 0; k < noNodesAnalysis1D3; k++) {
          m_wVolumeAnalysis(i, j, k) = wInt[i] * wInt[j] * (nDim == 3 ? wInt[k] : F1);
        }
      }
    }

    // Create analysis Vandermonde matrices to interpolate the solution from the
    // calculation nodes to the analysis nodes
    m_vdmAnalysis.clear();
    m_vdmAnalysis.resize(m_maxPolyDeg + 1);
    for(MInt polyDeg = m_minPolyDeg; polyDeg <= m_maxPolyDeg; polyDeg++) {
      m_vdmAnalysis[polyDeg].clear();
      m_vdmAnalysis[polyDeg].resize(m_maxNoNodes1D + 1);
    }


    // Create matrices
    for(MInt polyDeg = m_minPolyDeg; polyDeg <= m_maxPolyDeg; polyDeg++) {
      if(m_sbpMode) {
        for(MInt noNodes1D = m_minNoNodes1D; noNodes1D <= m_maxNoNodes1D; noNodes1D++) {
          m_vdmAnalysis[polyDeg][noNodes1D].resize(m_noAnalysisNodes, noNodes1D);
          const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];
          ASSERT(noNodes1D == m_noAnalysisNodes,
                 "For now SBP Error Analysis only supports direct evaluation at regular nodes "
                 "without interpolation...(noNodes1D = "
                     + to_string(noNodes1D) + ", noAnalysisNodes = " + to_string(m_noAnalysisNodes) + ")");

          dg::interpolation::calcLinearInterpolationMatrix(noNodes1D,
                                                           &interp.m_nodes[0],
                                                           m_noAnalysisNodes,
                                                           &m_interpAnalysis.m_nodes[0],
                                                           &m_vdmAnalysis[polyDeg][noNodes1D][0]);
        }
      } else {
        const MInt noNodes1D = polyDeg + 1;
        m_vdmAnalysis[polyDeg][noNodes1D].resize(m_noAnalysisNodes, noNodes1D);
        const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];
        dg::interpolation::calcPolynomialInterpolationMatrix(noNodes1D,
                                                             &interp.m_nodes[0],
                                                             m_noAnalysisNodes,
                                                             &m_interpAnalysis.m_nodes[0],
                                                             &interp.m_wBary[0],
                                                             &m_vdmAnalysis[polyDeg][noNodes1D][0]);
      }
    }
  }

  m_log << "done" << endl;
}

/**
 * \brief Calculates the coordinates of the integration nodes within each
 *        element.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-28
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initNodeCoordinates() {
  TRACE();

  m_log << "Initializing node coordinates... ";

  const MInt noElements = m_elements.size();

  // Iterate over all elements and calculate the integration node coordinates
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    calcElementNodeCoordinates(elementId);
  }

  m_log << "done" << endl;
}


/**
 * \brief Calculates the coordinates of the integration nodes within the
 *        element.
 *
 * \author Sven Berger
 * \date Februar 2015
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcElementNodeCoordinates(const MInt elementId) {
  // Create a shallow tensor for the node coordinates
  const MInt cellId = m_elements.cellId(elementId);
  const MInt lvl = grid().tree().level(cellId);
  const MInt polyDeg = m_elements.polyDeg(elementId);
  const MInt noNodes1D = m_elements.noNodes1D(elementId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  MFloatTensor nodeCoordinates(&m_elements.nodeCoordinates(elementId), noNodes1D, noNodes1D, noNodes1D3, nDim);

  // Iterate over all nodes and set the coordinates
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        // Node coordinates = cell center
        //     + 1/2 cell length * normalized ([-1,1]) node coordinate
        nodeCoordinates(i, j, k, 0) =
            grid().tree().coordinate(cellId, 0)
            + F1B2 * grid().cellLengthAtLevel(lvl) * m_interpolation[polyDeg][noNodes1D].m_nodes[i];
        nodeCoordinates(i, j, k, 1) =
            grid().tree().coordinate(cellId, 1)
            + F1B2 * grid().cellLengthAtLevel(lvl) * m_interpolation[polyDeg][noNodes1D].m_nodes[j];
        IF_CONSTEXPR(nDim == 3) {
          nodeCoordinates(i, j, k, 2) =
              grid().tree().coordinate(cellId, 2)
              + F1B2 * grid().cellLengthAtLevel(lvl) * m_interpolation[polyDeg][noNodes1D].m_nodes[k];
        }
      }
    }
  }
}


/**
 * \brief Calculates the inverse Jacobian for each element.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-23
 *
 *  The Jacobian describes the mapping of the arbitarily-sized (but
 *  squared or cubic) elements to the reference element
 *  [-1,1] x [-1,1] x [-1,1]. For regular elements it is always
 *  (element length)/2. However, since we only ever need the inverse Jacobian,
 *  the inverse is saved, i.e. 2/(element length), to make the code more
 *  efficient.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initJacobian() {
  TRACE();

  m_log << "Initializing inverse Jacobian determinant... ";

  const MInt noElements = m_elements.size();
  MFloat* invJacobians = &m_elements.invJacobian(0);

  for(MInt i = 0; i < noElements; i++) {
    const MInt cellId = m_elements.cellId(i);
    invJacobians[i] = F2 / grid().cellLengthAtCell(cellId);
  }

  m_log << "done" << endl;
}


/**
 * \brief Check all relevant bit properties in the cells.
 *
 * \author Michael Schlottke-Lakemper
 * \date   November 2017
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::checkCellProperties() {
  TRACE();

  m_log << "Initializing cell properties... ";

  const MInt noCells = grid().noCells();

  // Make sure that all cells that are not internal cells are marked as halo, as this is an
  // implicit assumption used throughout the solver
  for(MInt c = grid().noInternalCells(); c < noCells; c++) {
    if(!grid().tree().hasProperty(c, Cell::IsHalo)) {
      TERMM(1, "Cell " + to_string(c) + " should be marked as a halo cell");
    }
  }

  m_log << "done" << endl;
}

/// \brief Create for all elements and directions surfaces if necessary.
///
/// \author Sven Berger
/// \date   Februar 2014
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initSurfaces() {
  TRACE();

  const MInt noElements = m_elements.size();
  const MInt noDirs = 2 * nDim;
  const MInt* const cellIds = &m_elements.cellId(0);

  //////////////////////////////////////////////////////////////////////////////
  // 1) Create boundary surfaces
  //////////////////////////////////////////////////////////////////////////////
  m_log << "Creating boundary surfaces... ";
  //  m_boundarySurfacesOffset = m_surfaces.size();
  m_noBoundarySurfaces = 0;

  // Identify all interface cells, i.e., cells that are intersected by at least one geometry
  // element
  MBoolScratchSpace isInterface(grid().tree().size(), AT_, "isInterface");
  this->identifyBoundaryCells(&isInterface[0]);

  // Delete m_geometryIntersection if already allocated
  if(m_geometryIntersection) {
    delete m_geometryIntersection;
  }
  // Initialize geometry intersection object
  m_geometryIntersection = new GeometryIntersection<nDim>(&grid(), &geometry());

  // Loop over all elements and all directions and check if boundary
  // surfaces have to be created and which boundary condition ids are present
  MIntScratchSpace bcIds(noElements, noDirs, AT_, "bcIds");
  fill(bcIds.begin(), bcIds.end(), -1);
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    // Skip elements that are not at an interface
    const MInt cellId = cellIds[elementId];
    if(!isInterface[cellId]) {
      continue;
    }

    // Get boundary condition ids and save them to scratch
    array<MInt, 2 * nDim> ids = getBoundaryConditionIds(cellId);

    // Check for neighboring elements and prevent creating a boundary surface
    for(MInt dir = 0; dir < noDirs; dir++) {
      const MInt nId = grid().tree().neighbor(cellId, dir);
      if(ids[dir] > -1 && nId > -1 && needElementForCell(nId)) {
        ids[dir] = -1;
      }
    }

    copy(ids.begin(), ids.end(), &bcIds(elementId, 0));
  }

  // Find additional boundary conditions for boundary elements that are not
  // intersected by geometry but that do not have a neighbor. This is mainly
  // necessary to support concave domains and situations where the intersected
  // cell does not have an element because the cell center is outside of the
  // domain.
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = cellIds[elementId];
    for(MInt dir = 0; dir < noDirs; dir++) {
      // Skip if element already has boundary condition for this direction
      if(bcIds(elementId, dir) > -1) {
        continue;
      }

      // Skip if no neighbor cell exists
      if(!grid().tree().hasNeighbor(cellId, dir)) {
        continue;
      }

      // Skip if neighbor is not an interface cell
      const MInt neighborCellId = grid().tree().neighbor(cellId, dir);
      if(!isInterface[neighborCellId]) {
        continue;
      }

      // Skip if neighbor element exists
      if(needElementForCell(neighborCellId)) {
        continue;
      }

      // Get boundary condition ids for neighbor
      auto ids = getBoundaryConditionIds(neighborCellId);

      // Check if neighbor cell has boundary condition in the direction of
      // current element and if yes, save it for current element
      const MInt oppositeDir = 2 * (dir / 2) + 1 - (dir % 2);
      if(ids[oppositeDir] > -1) {
        bcIds(elementId, dir) = ids[oppositeDir];
      }
    }
  }

  // Find cut-off boundaries, i.e., cells that do not have a neighbor while at
  // the same time not being intersected by geometry
  if(m_useCutOffBoundaries) {
    // Reset cut-off boundary counter
    fill(m_noCutOffBoundarySurfaces.begin(), m_noCutOffBoundarySurfaces.end(), 0);

    // Find and mark cut-off boundaries
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt cellId = cellIds[elementId];
      for(MInt dir = 0; dir < noDirs; dir++) {
        // Skip if no cut-off boundary condition is specified for this direction
        if(m_cutOffBoundaryConditionIds[dir] < 0) {
          continue;
        }

        // Skip if element already has boundary condition for this direction
        if(bcIds(elementId, dir) > -1) {
          continue;
        }

        // Skip if neighbor cell exists
        if(grid().tree().hasNeighbor(cellId, dir)) {
          continue;
        }

        // Skip if neighbor of parent cell exists (will be handled by normal
        // h-refinement)
        if(grid().tree().parent(cellId) > -1 && grid().tree().hasNeighbor(grid().tree().parent(cellId), dir)) {
          continue;
        }

        // Set boundary condition id to cut-off boundary condition id
        bcIds(elementId, dir) = m_cutOffBoundaryConditionIds[dir];

        // Count number of cut-off surfaces
        m_noCutOffBoundarySurfaces[dir]++;
      }
    }
  }

  // Obtain list of unique boundary condition ids
  set<MInt> localUniqueBcIds(bcIds.begin(), bcIds.end());

  // Remove -1 (represents faces without a boundary condition id)
  localUniqueBcIds.erase(-1);

  // Exchange all local unique boundary conditions ids and determine all global
  // unique boundary condition ids. Each rank will create all boundary
  // conditions such that a global communication is possible for e.g.
  // reading/writing boundary condition restart data.

  // Local number of unique boundary condition ids
  MInt noBcIds = localUniqueBcIds.size();

  // Counts and offsets for Allgatherv
  ScratchSpace<MInt> countsBcIds(noDomains(), FUN_, "countsBcIds");
  ScratchSpace<MInt> offsetsBcIds(noDomains(), FUN_, "offsetsBcIds");

  // Gather the number of unique boundary condition ids on each domain
  MPI_Allgather(&noBcIds, 1, MPI_INT, &countsBcIds[0], 1, MPI_INT, mpiComm(), AT_, "noBcIds", "countsBcIds[0]");

  // Compute offsets
  offsetsBcIds[0] = 0;
  for(MInt dId = 1; dId < noDomains(); dId++) {
    offsetsBcIds[dId] = offsetsBcIds[dId - 1] + countsBcIds[dId - 1];
  }

  // Send buffer
  // Avoid dereferencing a zero length array
  ScratchSpace<MInt> localBcIds(max(noBcIds, 1), AT_, "localBcIds");
  std::copy(localUniqueBcIds.begin(), localUniqueBcIds.end(), &localBcIds[0]);

  // Receive buffer
  const MInt totalNoBcIds = offsetsBcIds[noDomains() - 1] + countsBcIds[noDomains() - 1];
  ScratchSpace<MInt> globalBcIds(totalNoBcIds, FUN_, "globalBcIds");

  // Gather and distribute all boundary condition ids from all domains
  MPI_Allgatherv(&localBcIds[0], noBcIds, MPI_INT, &globalBcIds[0], &countsBcIds[0], &offsetsBcIds[0], MPI_INT,
                 mpiComm(), AT_, "localBcIds[0]", "globalBcIds[0]");

  // Obtain list of global unique boundary condition ids
  set<MInt> uniqueBcIds(globalBcIds.begin(), globalBcIds.end());

  // Loop over all boundary conditions and create surfaces, then the boundary
  // condition object (so that the result is sorted by bcId)
  std::vector<MInt> boundaryConditionIds;
  std::vector<std::pair<MInt, MInt>> bcIntervals;
  for(const auto& boundaryConditionId : uniqueBcIds) {
    const MInt begin = m_noBoundarySurfaces;

    for(MInt elementId = 0; elementId < noElements; elementId++) {
      for(MInt dir = 0; dir < noDirs; dir++) {
        // Skip direction if boundary condition id does not match
        if(boundaryConditionId != bcIds(elementId, dir)) {
          continue;
        }
        if(hasSurface(elementId, dir)) {
          continue;
        }

        // Skip if boundary surface creation if periodic connectivity is given in this direction
        // Determine periodic direction pDir(+x,+y,+z) based on dir(+-x,+-y,+-z)
        // --> pDir = (MInt) dir/2
        if(grid().periodicCartesianDir(dir / 2) != 0) {
          // Add additional check to avoid unintentional errors (i.e., users must set a periodic
          // direction *and* use bcId=0 on all relevant boundaries)
          if(boundaryConditionId != 0) {
            TERMM(1,
                  "If you use periodic BCs, you must set the BC id of corresponding boundaries "
                  "to zero.");
          }
          continue;
        }

        initBoundarySurface(elementId, dir);
      }
    }
    const MInt end = m_noBoundarySurfaces;

    // Store boundary condition id and the corresponding surface interval
    // Boundary conditions are created & initialized after all remaining
    // surfaces are created
    boundaryConditionIds.push_back(boundaryConditionId);
    bcIntervals.emplace_back(begin, end);
  }

  m_log << "done" << endl;

  //////////////////////////////////////////////////////////////////////////////
  // 2) Create inner surfaces
  //////////////////////////////////////////////////////////////////////////////
  m_log << "creating inner surfaces... ";

  m_innerSurfacesOffset = m_surfaces.size();
  m_noInnerSurfaces = 0;

  for(MInt elementId = 0; elementId < noElements; elementId++) {
    for(MInt dir = 0; dir < noDirs; dir++) {
      if(hasSurface(elementId, dir)) {
        continue;
      }
      initInnerSurface(elementId, dir);
    }
  }

  m_log << "done" << endl;

  //////////////////////////////////////////////////////////////////////////////
  // 3) Create MPI surfaces
  //////////////////////////////////////////////////////////////////////////////
  m_log << "creating MPI surfaces... ";

  m_mpiSurfacesOffset = m_surfaces.size();
  m_noMpiSurfaces = 0;

  if(hasMpiExchange()) {
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      for(MInt dir = 0; dir < noDirs; dir++) {
        const MInt cellId = m_elements.cellId(elementId);
        const MInt nghbrId = grid().tree().neighbor(cellId, dir);
        MBool partLvlAncestorNghbr = false;

        // In case of partition level shifts there can be elements at level jumps that have both
        // internal and halo cells as neighbors on the higher level. Check if the neighbor cell is a
        // candidate for such a case
        if(nghbrId > -1 && grid().tree().hasProperty(nghbrId, Cell::IsPartLvlAncestor)
           && grid().tree().hasChildren(nghbrId)) {
          // Check all childs of the neighbor for a cell and continue process with initMpiSurface()
          // below to create any missing h-mpi surfaces (and skip internal childs of the neighbor)
          for(MInt child = 0; child < IPOW2(nDim); child++) {
            const MInt childId = grid().tree().child(nghbrId, child);
            if(childId > -1 && grid().tree().hasProperty(childId, Cell::IsHalo)) {
              partLvlAncestorNghbr = true;
            }
          }
        }

        if(hasSurface(elementId, dir) && !partLvlAncestorNghbr) {
          continue;
        }
        initMpiSurface(elementId, dir);
      }
    }
  }

  m_log << "done" << endl;

  // Reset previous boundary condition objects
  m_boundaryConditions.clear();

  // Create & initialize boundary condition
  m_log << "creating and initializing boundary conditions... ";
  for(MUint bcId = 0; bcId < boundaryConditionIds.size(); bcId++) {
    m_boundaryConditions.emplace_back(make_bc(boundaryConditionIds[bcId]));

    const MInt begin = bcIntervals[bcId].first;
    const MInt end = bcIntervals[bcId].second;
    m_boundaryConditions[bcId]->init(begin, end);
  }
  m_log << "done" << endl;

  // Sanity check: all elements should have surfaces in each direction
  m_log << "Sanity check: all elements must have surfaces in each "
           "direction... ";
  MInt missing = 0;
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    for(MInt dir = 0; dir < noDirs; dir++) {
      if(!hasSurface(elementId, dir)) {
        // Write surface information to log
        m_log << "\nmissing surface for element " << elementId << " in direction " << dir << " (coordinates: ";
        const MInt cellId = m_elements.cellId(elementId);
        for(MInt i = 0; i < nDim; i++) {
          m_log << grid().tree().coordinate(cellId, i) << " ";
        }
        m_log << ")" << std::endl;
        m_log.flush();

        // Count number of missing surfaces for abort message
        missing++;
      }
    }
  }

  // Abort if there are missing surfaces
  if(missing) {
    TERMM(1, "There are " + to_string(missing)
                 + " surfaces missing (check m_log for more details). Did "
                   "you maybe forget to specify cut-off boundary conditions?");
  } else {
    m_log << "done" << endl;
  }

  // Ensure that there are no boundary surfaces if all directions are periodic
  if(grid().periodicCartesianDir(0) && grid().periodicCartesianDir(1) && (nDim == 2 || grid().periodicCartesianDir(2))
     && m_noBoundarySurfaces > 0) {
    TERMM(1, "There are boundary surfaces although every direction has a periodic boundaries. "
             "Please check what is wrong here.\n");
  }
}


/// \brief Check if a surface exists for the element in the given direction.
///
/// \author Sven Berger
/// \date   November 2014
///
/// \param[in] elementId Id of element to check.
/// \param[in] dir Direction element-> surface (0..5 -> -x,+x,-y,+y,-z,+z)
///
/// \return True if surface exists, false otherwise.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::hasSurface(const MInt elementId, const MInt dir) {
  TRACE();
  return m_elements.surfaceIds(elementId, dir) > -1;
}


/// \brief Check if the surface of the element in the given direction
/// is a boundary surface. If it is init, a boundary surface.
///
/// \author Sven Berger
/// \date   Februar 2014
///
/// \param[in] elementId Id of element to check.
/// \param[in] dir Direction element->surface (0..5 -> -x,+x,-y,+y,-z,+z)
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::initBoundarySurface(const MInt elementId, const MInt dir) {
  TRACE();

  const MInt surfaceId = createSurface(elementId, dir);
  m_noBoundarySurfaces++;

  return surfaceId;
}


/// \brief Check if the surface of the element in the given
/// direction is an inner surface without h/p-refinement.
/// If it is, init an inner surface.
///
/// \author Sven Berger
/// \date   Februar 2014
///
/// \param[in] elementId Id of element to check.
/// \param[in] dir Direction element->surface (0..5 -> -x,+x,-y,+y,-z,+z)
///
/// If the neighboring cell/element is h-refined with respect to the current
/// (coarse) cell/element, no surfaces are created. Surfaces are always created
/// from the refined to the coarses elements, and there is one surface per
/// refined element (i.e. in 2D there will be at most 2 surfaces for each
/// fine-coarse interface, and in 3D there will be at most 5 surfaces). This is
/// different form initMpiSurface, where surfaces are created from the point of
/// view of the coarse element as well.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::initInnerSurface(const MInt elementId, const MInt dir) {
  TRACE();

  MInt cellId = m_elements.cellId(elementId);

  // No surface if there is no neighbor cell on current and parent level
  if(!grid().tree().hasNeighbor(cellId, dir)) {
    const MInt parentId = grid().tree().parent(cellId);
    if(parentId > -1 && grid().tree().hasNeighbor(parentId, dir)) {
      // If parent cell exists and it has neighbor in correct direction, use
      // parent to determine neighbor cell
      cellId = parentId;
    } else {
      // Otherwise there is no neighbor in the specified direction, so quit
      // without creating a surface
      return -1;
    }
  }

  const MInt nghbrId = grid().tree().neighbor(cellId, dir);

  // If child exists skip since surfaces are created on highest level
  if(grid().tree().hasChildren(nghbrId)) {
    return -1;
  }

  // No surface if neighbor is a halo cell
  if(grid().tree().hasProperty(nghbrId, Cell::IsHalo)) {
    return -1;
  }

  // No surface if there is no neighbor element
  const MInt nghbrElementId = m_elements.getElementByCellId(nghbrId);
  if(nghbrElementId == -1) {
    return -1;
  }

  // Create surface
  const MInt surfaceId = createSurface(elementId, dir);
  m_noInnerSurfaces++;

  return surfaceId;
}


/// \brief Check if the surface of the element in the given direction
/// is a MPI surface. If it is, init a MPI surface.
///
/// \author Sven Berger
/// \date   Februar 2014
///
/// \param[in] elementId Id of element to check.
/// \param[in] dir Direction element->surface (0..5 -> -x,+x,-y,+y,-z,+z)
///
/// If the corresponding halo cell is refined in comparison to the current cell,
/// *two* surfaces are created. This is different to initInnterSurface, where
/// surfaces are always created from the fine elements, never from the coarse
/// element.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::initMpiSurface(const MInt elementId, const MInt dir) {
  TRACE();

  MInt cellId = m_elements.cellId(elementId);

  // No surface if there is no neighbor cell on current and parent level
  if(!grid().tree().hasNeighbor(cellId, dir)) {
    const MInt parentId = grid().tree().parent(cellId);
    if(parentId > -1 && grid().tree().hasNeighbor(parentId, dir)) {
      // If parent cell exists and it has neighbor in correct direction, use
      // parent to determine neighbor cell
      cellId = parentId;
    } else {
      // Otherwise there is no neighbor in the specified direction, so quit
      // without creating a surface
      return -1;
    }
  }

  const MInt nghbrId = grid().tree().neighbor(cellId, dir);

  // Surface only with halo cells as neighbor OR
  // (in case of partition level shifts) if the neighbor has children (createHMPISurfaces will
  // handle/check such cases)
  if(!grid().tree().hasProperty(nghbrId, Cell::IsHalo) && !grid().tree().hasChildren(nghbrId)) {
    return -1;
  }

  MInt surfaceId = -1;

  // If child exists create multiple MPI surfaces
  if(grid().tree().hasChildren(nghbrId)) {
    // Create multiple surfaces
    surfaceId = createHMPISurfaces(elementId, dir);
  } else {
    // Create surface
    surfaceId = createSurface(elementId, dir);
    m_noMpiSurfaces++;
  }

  return surfaceId;
}


/// \brief Determine if element is boundary is cut by geometry elements and
///        return corresponding boundary condition ids.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014
///
/// \param[in] elementId Id of element for which to return boundary condition
///                      ids.
///
/// \return An array of length 2*nDim with the boundary condition id for each
///         element face (or -1 if the element face is not at a physical
///         boundary).
template <MInt nDim, class SysEqn>
array<MInt, 2 * nDim> DgCartesianSolver<nDim, SysEqn>::getBoundaryConditionIds(const MInt cellId) {
  TRACE();

  // Set target region to be the minimum and maximum cell coordinates, i.e. of
  // the 4 (2D) or 6 (3D) vertices
  array<MFloat, 2 * nDim> targetRegion;
  for(MInt dir = 0; dir < nDim; dir++) {
    targetRegion[dir] = grid().tree().coordinate(cellId, dir) - 0.5 * grid().cellLengthAtCell(cellId);
    targetRegion[dir + nDim] = grid().tree().coordinate(cellId, dir) + 0.5 * grid().cellLengthAtCell(cellId);
  }

  // Get cell coordinates and half length
  const MFloat cellHalfLength = 0.5 * grid().cellLengthAtCell(cellId);
  array<MFloat, nDim> coordinates;
  for(MInt i = 0; i < nDim; i++) {
    coordinates[i] = grid().tree().coordinate(cellId, i);
  }

  // Check if there are elements in the cell
  const MBool hasElement = needElementForCell(cellId);

  // Get list of intersecting elements
  std::vector<MInt> nodeList;
  geometry().getIntersectionElements(&targetRegion[0], nodeList, cellHalfLength, &coordinates[0]);

  array<MInt, 2 * nDim> bcIds;
  bcIds.fill(-1);

  // Process each geometric element ('node')
  static constexpr MFloat standardBasis[MAX_SPACE_DIMENSIONS][MAX_SPACE_DIMENSIONS] = {
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  for(MInt n = 0; n < (signed)nodeList.size(); n++) {
    // Obtain element normal
    const GeometryElement<nDim>& elem = geometry().elements[nodeList[n]];
    // Obtain geometry vertices
    array<MFloat, nDim * nDim> geometryPoints;
    elem.getVertices(&geometryPoints[0]);

    // Obtain normal
    array<MFloat, nDim> normal;
    elem.calcNormal(&geometryPoints[0], &normal[0]);

    // Check orientation of element normal (axis aligned)
    MInt orientation = -1;
    for(MInt d = 0; d < nDim; d++) {
      const MFloat dot = inner_product(&normal[0], &normal[nDim], standardBasis[d], 0.0);
      if(approx(std::fabs(dot), 1.0, MFloatEps)) {
        orientation = d;
        break;
      }
    }

    // Set boundary id
    if(orientation != -1) {
      // ... for axis-aligned cases

      // Find closest cell face (i.e. determine +ve or -ve placement w.r.t. to the
      // cell center)
      array<MFloat, nDim> centroid;
      elem.calcCentroid(&geometryPoints[0], &centroid[0]);
      const MFloat pos = grid().tree().coordinate(cellId, orientation) + 0.5 * grid().cellLengthAtCell(cellId);
      const MFloat neg = grid().tree().coordinate(cellId, orientation) - 0.5 * grid().cellLengthAtCell(cellId);
      const MInt direction = (std::fabs(neg - centroid[orientation]) < std::fabs(pos - centroid[orientation]))
                                 ? 2 * orientation
                                 : 2 * orientation + 1;

      // Abort if boundary condition id for this direction was already set to a
      // different value
      if(bcIds[direction] > -1 && bcIds[direction] != elem.m_bndCndId) {
        stringstream ss;
        ss << "Bad cell: " << cellId << ". Cell has two different boundary condition ids (" << bcIds[direction] << ", "
           << elem.m_bndCndId << ") in direction " << direction << ". This feature is not yet implemented." << endl;
        TERMM(1, ss.str());
      }

      // Set boundary id
      bcIds[direction] = elem.m_bndCndId;

    } else {
      // ... for non axis-aligned cases

      {
        // TODO labels:DG,GEOM,toenhance Needs to be adapted for cases where the outer boundary is non axis-aligned
        //
        //  /DESCRIPTION/
        //
        //  Currently, the excerpt of code below only works for the cases where the
        //  INNER boundary only is non axis-aligned, e.g., an airfoil (2D), or a wing (3D)
        //  inside the domain. The outer boundary HAS to be axis aligned, e.g., has to be
        //  a square (2D) or cube (3D). This happens due to the following condition:
        //
        //    direction = (normal[i] < 0) ? 2 * i : 2 * i + 1;
        //
        //  That condition is only capable of assigning the correct direction for INNER
        //  boundaries with the normal vector pointing outwards the geometry, i.e., pointing
        //  inside the domain (for outer boundaries this just need to be flipped in order to
        //  work). However, we currently can't find a proper way to to differentiate if a geometry
        //  element belongs to the inner or outer boundary.
        //
        //  /POSSIBLE FIX/
        //
        //  A simplistic way of dealing with this issue can be by creating different
        //  boundary condition IDs for inner and outer boundaries. Then, the code below
        //  can be modified to check if the current geometry element has a m_bndCndId
        //  that corresponds to a inner or outer boundary and assign the direction correctly.
        //  (Remember to also change the geometry.toml files accordingly)
        //
        for(MInt i = 0; i < nDim; i++) {
          if(!approx(normal[i], 0.0, MFloatEps)) {
            MInt direction = (normal[i] < 0) ? 2 * i : 2 * i + 1;
            // If elements exist in the cell, use the opposite normal direction of geometry face
            if(hasElement) {
              direction = 2 * (direction / 2) + 1 - (direction % 2);
            }

            // Abort if h-ref is used along a non axis-aligned boundary
            // First check if a neighbor on same level can theoretically exists -> if not a coarser nghbr might exists
            array<MFloat, nDim> nghbrCellCoord;
            for(MInt d = 0; d < nDim; ++d) {
              nghbrCellCoord[d] = grid().tree().coordinate(cellId, d);
            }
            nghbrCellCoord[i] += (2 * (direction % 2) - 1) * grid().cellLengthAtCell(cellId);
            const MBool nghbrIsInside = geometry().pointIsInside2(&nghbrCellCoord[0]);
            if((!nghbrIsInside && grid().tree().level(cellId) != getLevelOfDirectNeighborCell(cellId, direction))
               || (grid().tree().level(cellId) < getLevelOfDirectNeighborCell(cellId, direction))) {
              stringstream ss;
              ss << "Bad cell: " << cellId << ". There is h-refinement being used along a non axis-aligned surface in "
                 << "direction " << direction << ". This feature is not yet implemented." << endl;
              TERMM(1, ss.str());
            }

            // Abort if boundary condition id for this direction was already set to a
            // different value
            if(bcIds[direction] > -1 && bcIds[direction] != elem.m_bndCndId) {
              stringstream ss;
              ss << "Bad cell: " << cellId << ". Cell has two different boundary condition ids (" << bcIds[direction]
                 << ", " << elem.m_bndCndId << ") in direction " << direction
                 << ". This feature is not yet implemented." << endl;
              TERMM(1, ss.str());
            }

            // Set boundary id
            bcIds[direction] = elem.m_bndCndId;
          }
        }
      }

      // Candidates for geometry intersection
      std::vector<CutCandidate<nDim>> cutCandidates(1);
      // Assign cellId as the only candidate
      cutCandidates[0].cellId = cellId;

      m_geometryIntersection->computeCutPointsFromSTL(cutCandidates);

      const MInt noCutPoints = cutCandidates[0].noCutPoints;

      // Set target points according to the number of cut points and dimensions (2D or 3D)
      vector<MFloat> targetPoints(noCutPoints * nDim);
      for(MInt cutPoint = 0; cutPoint < noCutPoints; cutPoint++) {
        for(MInt i = 0; i < nDim; i++) {
          targetPoints[(cutPoint * nDim) + i] = cutCandidates[0].cutPoints[cutPoint][i];
        }
      }

      // Recalculate the normal for cutPoints
      elem.calcNormal(&targetPoints[0], &normal[0]);

      // Check the new orientation for non-axis aligned cases
      orientation = 0;
      MFloat maxNormal = fabs(normal[0]);
      for(MInt d = 0; d < nDim; d++) {
        if(fabs(normal[d]) > maxNormal) {
          orientation = d;
          maxNormal = fabs(normal[d]);
        }
      }

      {
        // Set boundary ids for faces that are not intersected by geometry
        // Find closest cell face (i.e. determine +ve or -ve placement w.r.t. to the cell center)
        array<MFloat, nDim> centroid;
        elem.calcCentroid(&targetPoints[0], &centroid[0]);
        const MFloat pos = grid().tree().coordinate(cellId, orientation) + 0.5 * grid().cellLengthAtCell(cellId);
        const MFloat neg = grid().tree().coordinate(cellId, orientation) - 0.5 * grid().cellLengthAtCell(cellId);
        MInt direction = (std::fabs(neg - centroid[orientation]) < std::fabs(pos - centroid[orientation]))
                             ? 2 * orientation
                             : 2 * orientation + 1;

        // Needed when a geometry vertice is inside the cell and the wrong direction is assigned.
        // Mirror direction assignment if:
        // 1. No elements exist in the cell AND no neighbor exists in the specified direction OR
        // 2. Element exists in the cell AND neighbor exists in the speficied direction
        if((!hasElement && !grid().tree().hasNeighbor(cellId, direction))
           || (hasElement && grid().tree().hasNeighbor(cellId, direction))) {
          direction = 2 * (direction / 2) + 1 - (direction % 2);
        }
        // Abort if h-ref is used along a non axis-aligned boundary
        // First check if a neighbor on same level can theoretically exists -> if not a coarser nghbr might exists
        array<MFloat, nDim> nghbrCellCoord;
        for(MInt d = 0; d < nDim; ++d) {
          nghbrCellCoord[d] = grid().tree().coordinate(cellId, d);
        }
        nghbrCellCoord[orientation] += (2 * (direction % 2) - 1) * grid().cellLengthAtCell(cellId);
        const MBool nghbrIsInside = geometry().pointIsInside2(&nghbrCellCoord[0]);
        if((!nghbrIsInside && grid().tree().level(cellId) != getLevelOfDirectNeighborCell(cellId, direction))
           || (grid().tree().level(cellId) < getLevelOfDirectNeighborCell(cellId, direction))) {
          stringstream ss;
          ss << "Bad cell: " << cellId << ". There is h-refinement being used along a non axis-aligned surface "
             << "in direction " << direction << ". This feature is not yet implemented." << endl;
          TERMM(1, ss.str());
        }
        // Abort if boundary condition id for this direction was already set to a
        // different value
        if(bcIds[direction] > -1 && bcIds[direction] != elem.m_bndCndId) {
          stringstream ss;
          ss << "Bad cell: " << cellId << ". Cell has two different boundary condition ids (" << bcIds[direction]
             << ", " << elem.m_bndCndId << ") in direction " << direction << ". This feature is not yet implemented."
             << endl;
          TERMM(1, ss.str());
        }

        // Set boundary id
        bcIds[direction] = elem.m_bndCndId;
      }
    }
  }

  return bcIds;
}


/*
 * \brief Creates multiple surfaces between an element and its refined neighbor
 * elements for an inter-domain (MPI) boundary.
 *
 * \author Sven Berger
 * \date   April 2015
 *
 * \details This method assumes that elementId is always a valid element, i.e.
 * that it is at an MPI boundary and that it is the coarser element.
 *
 * \param[in] elementId The element id of the coarse element.
 *
 * \param[in] dir Direction id of neighbor elements relative to the main element
 *                (i.e. 0 = -x, 1 = +x, 2 = -y etc.)
 *
 * \return The surface id of the first newly created surface.
 */
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::createHMPISurfaces(const MInt elementId, const MInt dir) {
  // TRACE();

  // Determine neighbor element
  const MInt cellId = m_elements.cellId(elementId);

  // Get neighbor cells from child level
  vector<MInt> nghbrIds;
  const MInt parentNghbrId = grid().tree().neighbor(cellId, dir);

  // Arrays give child id of child level neighbors for a neighbor cell in a
  // given direction
  static constexpr MInt dirNghbrs[6][4] = {
      {1, 3, 5, 7}, // -x direction
      {0, 2, 4, 6}, // +x direction
      {2, 3, 6, 7}, // -y direction
      {0, 1, 4, 5}, // +y direction
      {4, 5, 6, 7}, // -z direction
      {0, 1, 2, 3}  // +z direction
  };

  MInt noNonHaloChilds = 0;
  // Store every existing child level neighbor element
  for(MInt i = 0; i < 2 * (nDim - 1); i++) {
    if(grid().tree().hasChild(parentNghbrId, dirNghbrs[dir][i])) {
      // Check if child is an internal cell (partition level shift), skip
      const MInt childId = grid().tree().child(parentNghbrId, dirNghbrs[dir][i]);
      if(!grid().tree().hasProperty(childId, Cell::IsHalo)) {
        ASSERT(grid().tree().hasProperty(parentNghbrId, Cell::IsPartLvlAncestor),
               "this case is only valid with a neighboring partition level ancestor cell");
        noNonHaloChilds++;
        continue;
      }
      nghbrIds.push_back(grid().tree().child(parentNghbrId, dirNghbrs[dir][i]));
    }
  }

  // Partition level shift: only neighboring non-halo childs, nothing to do
  if(noNonHaloChilds == 2 * (nDim - 1)) {
    return -1;
  }

  if(nghbrIds.empty()) {
    stringstream errorMessage;
    errorMessage << "Error: child cells of neighboring halo cell not found. "
                    "Try setting 'newCreateWindowCellsF = 0'."
                 << endl;
    TERMM(1, errorMessage.str());
  }

  // Create a surface for each child level neighbor element
  MInt firstSurfaceId = -1;
  for(const auto& nghbrId : nghbrIds) {
    // First, create normal surface
    const MInt srfcId = createSurface(elementId, dir, nghbrId);

    // Store first created surface for the return value
    if(firstSurfaceId == -1) {
      firstSurfaceId = srfcId;
    }

    // Then, add surface id of new surface to corresponding h-element
    for(MInt hElementId = 0; hElementId < m_helements.size(); hElementId++) {
      // Find h-element of the current element
      if(m_helements.elementId(hElementId) == elementId) {
        // Add surface id at first unused location
        for(MInt hSurfId = 0; hSurfId < 2 * (nDim - 1); hSurfId++) {
          if(m_helements.hrefSurfaceIds(hElementId, dir, hSurfId) == -1) {
            m_helements.hrefSurfaceIds(hElementId, dir, hSurfId) = srfcId;
            break;
          }
        }
        break;
      }
    }

    // Set h-refinement specific variables
    m_surfaces.fineCellId(srfcId) = nghbrId;

    // Increase counter
    m_noMpiSurfaces++;
  }

  return firstSurfaceId;
}


/*
 * \brief Creates a surface between two neighboring elements.
 *
 * \author Michael Schlottke
 * \date   October 2012
 *
 * \details This method assumes that elementId is always a valid element.
 *          nghbrId might be a non-existing element (i.e. -1).
 *
 * \param[in] elementId The element id of the "main" adjacent element. This must
 *                      be a valid element!
 * \param[in] dir Direction id of neighbor element relative to the main element
 *                (i.e. 0 = -x, 1 = +x, 2 = -y etc.)
 * \param[in] nghbrId CellId of the neighbor cell (needs to be set *only* for
 *                    halo cells)
 *
 * \return The surface id of the newly created surface.
 */
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::createSurface(const MInt elementId, const MInt dir, MInt nghbrId) {
  // TRACE();

  // Determine neighbor element
  const MInt cellId = m_elements.cellId(elementId);
  const MInt parentId = grid().tree().parent(cellId);

  // Get neighbor from current or parent cell (or leave id as -1 if no neighbor
  // found)
  if(nghbrId == -1) {
    nghbrId = grid().tree().neighbor(cellId, dir);
    if(nghbrId == -1 && parentId > -1 && grid().tree().hasNeighbor(parentId, dir)) {
      nghbrId = grid().tree().neighbor(parentId, dir);
    }
  }
  const MInt nghbrElementId = m_elements.getElementByCellId(nghbrId);

  // Determine new surface id and create the surface in collector
  const MInt srfcId = m_surfaces.size();
  m_surfaces.append();

  // Calculate opposite direction
  const MInt oppositeDir = 2 * (dir / 2) + 1 - (dir % 2);

  // Calculate and set the surface orientation (0: x, 1: y, 2: z)
  const MInt orientation = dir / 2;
  m_surfaces.orientation(srfcId) = orientation;

  // Calculate and set the neighboring element ids (element in -ve direction is
  // first)
  const MInt nghbrSideId = dir % 2;
  const MInt elementSideId = (nghbrSideId + 1) % 2;
  m_surfaces.nghbrElementIds(srfcId, nghbrSideId) = nghbrElementId;
  m_surfaces.nghbrElementIds(srfcId, elementSideId) = elementId;

  // Set surfaceId for element
  m_elements.surfaceIds(elementId, dir) = srfcId;

  // Set the surfaceId for the neighbor element if it is not already set
  if(nghbrElementId > -1 && m_elements.surfaceIds(nghbrElementId, oppositeDir) == -1) {
    m_elements.surfaceIds(nghbrElementId, oppositeDir) = srfcId;
  }

  // Get levels of element and neighbor element (if existing)
  const MInt level = getLevelByElementId(elementId);
  MInt nghbrLvl = -1;
  if(nghbrElementId > -1) {
    nghbrLvl = getLevelByElementId(nghbrElementId);
  }

  // Reset fine cell id
  m_surfaces.fineCellId(srfcId) = -1;

  // Set h-elements if a neighbor with a lower level exists
  if(nghbrLvl != -1 && level > nghbrLvl) {
    // Iterate over h-elements
    for(MInt hElementId = 0; hElementId < m_helements.size(); hElementId++) {
      // Find h-element for neighbor element
      if(m_helements.elementId(hElementId) == nghbrElementId) {
        m_surfaces.fineCellId(srfcId) = cellId;
        // Add surface id at first unused location
        for(MInt hSurfId = 0; hSurfId < 2 * (nDim - 1); hSurfId++) {
          if(m_helements.hrefSurfaceIds(hElementId, oppositeDir, hSurfId) == -1) {
            m_helements.hrefSurfaceIds(hElementId, oppositeDir, hSurfId) = srfcId;
            break;
          }
        }
        break;
      }
    }
  }

  // Set polynomial degree to maximum of both elements' polynomial degrees
  // Set number of nodes 1D to maximum of both elements' number of nodes 1D
  if(nghbrElementId > -1) {
    m_surfaces.polyDeg(srfcId) = max(m_elements.polyDeg(elementId), m_elements.polyDeg(nghbrElementId));
    m_surfaces.noNodes1D(srfcId) = max(m_elements.noNodes1D(elementId), m_elements.noNodes1D(nghbrElementId));
  } else {
    m_surfaces.polyDeg(srfcId) = m_elements.polyDeg(elementId);
    m_surfaces.noNodes1D(srfcId) = m_elements.noNodes1D(elementId);
  }


  // Reset level in case there is a neighbor cell
  if(nghbrId > -1) {
    nghbrLvl = grid().tree().level(nghbrId);
  }

  // Compute surface coordinates by using the higher level cell
  if(level > nghbrLvl) {
    const MFloat cellLength = grid().cellLengthAtCell(cellId);
    m_surfaces.coords(srfcId, orientation) =
        grid().tree().coordinate(cellId, orientation) + (static_cast<MFloat>(nghbrSideId) - F1B2) * cellLength;
    for(MInt i = 0; i < nDim; i++) {
      if(i != orientation) {
        m_surfaces.coords(srfcId, i) = grid().tree().coordinate(cellId, i);
      }
    }
  } else {
    const MFloat cellLength = grid().cellLengthAtCell(nghbrId);
    m_surfaces.coords(srfcId, orientation) =
        grid().tree().coordinate(nghbrId, orientation) + (static_cast<MFloat>(1 - nghbrSideId) - F1B2) * cellLength;
    for(MInt i = 0; i < nDim; i++) {
      if(i != orientation) {
        m_surfaces.coords(srfcId, i) = grid().tree().coordinate(nghbrId, i);
      }
    }
  }

  // Set internal side id
  const MInt elementIdL = m_surfaces.nghbrElementIds(srfcId, 0);
  const MInt elementIdR = m_surfaces.nghbrElementIds(srfcId, 1);
  // If a valid nghbrCell exists only one side -> mark this side
  // otherwise keep -1
  m_surfaces.internalSideId(srfcId) = -1;
  if(elementIdL == -1) {
    m_surfaces.internalSideId(srfcId) = 1;
  }
  if(elementIdR == -1) {
    m_surfaces.internalSideId(srfcId) = 0;
  }

  // Calculate node coordinates
  calcSurfaceNodeCoordinates(srfcId);

  // Determine & set global surface id
  const MLong globalCellId = grid().tree().globalId(m_elements.cellId(elementId));
  if(nghbrId == -1 || nghbrLvl == -1) {
    // If no neighbor cell exists, take the global id from the cell
    m_surfaces.globalId(srfcId) = 2 * nDim * globalCellId + dir;
  } else {
    // If neighbor cell exists, take the higher level cell or if both have the
    // same level take the cell with the lower global cell id
    const MLong globalNghbrId = grid().tree().globalId(nghbrId);
    if(level > nghbrLvl || (globalCellId < globalNghbrId && level == nghbrLvl)) {
      m_surfaces.globalId(srfcId) = 2 * nDim * globalCellId + dir;
    } else {
      m_surfaces.globalId(srfcId) = 2 * nDim * globalNghbrId + oppositeDir;
    }
  }

  return srfcId;
}


/*
 * \brief Calc coordinates of the nodes on the surface.
 *
 * \author Sven Berger
 * \date   Februar 2015
 *
 * \param[in] srfcId The surface id of the surface for which the node
 *  coordinates are calculated.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcSurfaceNodeCoordinates(const MInt srfcId) {
  const MInt internalSide = (m_surfaces.internalSideId(srfcId) <= 0) ? 0 : 1;
  const MInt elementId = m_surfaces.nghbrElementIds(srfcId, internalSide);
  const MInt cellId = m_elements.cellId(elementId);

  // Compute surface node coordinates
  const MInt polyDeg = m_surfaces.polyDeg(srfcId);
  const MInt noNodes1D = m_surfaces.noNodes1D(srfcId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloat length = grid().cellLengthAtCell(cellId);
  MFloatTensor nodeCoordinates(&m_surfaces.nodeCoords(srfcId), noNodes1D, noNodes1D3, nDim);

  // Iterate over all nodes and set the coordinates
  // Node coordinates = surface center
  //   + (1/2 element length * normalized ([-1,1]) node coordinate)
  switch(m_surfaces.orientation(srfcId)) {
    case 0: // x-direction
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          nodeCoordinates(j, k, 0) = m_surfaces.coords(srfcId, 0);
          nodeCoordinates(j, k, 1) =
              m_surfaces.coords(srfcId, 1) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[j];
          IF_CONSTEXPR(nDim == 3) {
            nodeCoordinates(j, k, 2) =
                m_surfaces.coords(srfcId, 2) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[k];
          }
        }
      }
      break;

    case 1: // y-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          nodeCoordinates(i, k, 0) =
              m_surfaces.coords(srfcId, 0) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[i];
          nodeCoordinates(i, k, 1) = m_surfaces.coords(srfcId, 1);
          IF_CONSTEXPR(nDim == 3) {
            nodeCoordinates(i, k, 2) =
                m_surfaces.coords(srfcId, 2) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[k];
          }
        }
      }
      break;

    case 2: // z-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D3; j++) {
          nodeCoordinates(i, j, 0) =
              m_surfaces.coords(srfcId, 0) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[i];
          nodeCoordinates(i, j, 1) =
              m_surfaces.coords(srfcId, 1) + F1B2 * length * m_interpolation[polyDeg][noNodes1D].m_nodes[j];
          IF_CONSTEXPR(nDim == 3) { nodeCoordinates(i, j, 2) = m_surfaces.coords(srfcId, 2); }
        }
      }
      break;

    default:
      mTerm(1, AT_, "Bad orientation");
  }
}


/* Determines if MPI exchange has to take place, i.e. if computation is in parallel or a periodic
 * BC is specified If there is another condition than periodicity where an MPI exchange shall
 * occur, please modify
 */
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::hasMpiExchange() const {
  TRACE();

  // If InitMpiExchange has not been performed yet: Initialize it, if the grid has neighbor
  // domains. Else if InitMpiExchange has been performed: m_noExchangeNghbrDomains describes
  // number of neighbor domains, with which actually data has to be exchanged. Look also in
  // initMpiExchange().
  const MBool neighbors = m_isInitMpiExchange ? m_noExchangeNghbrDomains > 0 : grid().noNeighborDomains() > 0;

  // Return true if executed on more than one domain *or* if there exist neighbor domains (the
  // latter may occur in case of periodicity)
  return (noDomains() > 1 || neighbors);
}


/*
 * \brief Initialize data arrays that are needed for MPI communication, i.e.
 *        buffers and persistent MPI requests.
 *
 * \author Michael Schlottke
 * \date   July 2013
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initMpiExchange() {
  TRACE();

  m_log << "Initialize MPI exchange... ";

  // Create map domain ids to exchange neighbor domain ids
  map<MInt, MInt> exchangeNghbrDomainMap;

  // Reset number of exchange neighbor domains (i.e. domains with which actual
  // data needs to be exchanged, as opposed to neighbor domains, with which just
  // halo/window cells are shared).
  m_noExchangeNghbrDomains = 0;

  // Initialize resize containers to max. number of exchange domains
  m_exchangeNghbrDomains.resize(grid().noNeighborDomains());
  m_mpiSurfaces.resize(grid().noNeighborDomains());

  // Fill MPI surfaces container with surfaces that need to be exchanged with
  // the respective domain
  const MInt begin = m_mpiSurfacesOffset;
  const MInt end = m_mpiSurfacesOffset + m_noMpiSurfaces;
  for(MInt srfcId = begin; srfcId < end; srfcId++) {
    const MInt internalSideId = m_surfaces.internalSideId(srfcId);
    const MInt elementId = m_surfaces.nghbrElementIds(srfcId, internalSideId);
    const MInt internalCellId = m_elements.cellId(elementId);
    const MInt nghbrCellDir = 2 * m_surfaces.orientation(srfcId) + 1 - internalSideId;

    // Get halo cell
    MInt haloCellId = grid().tree().neighbor(internalCellId, nghbrCellDir);
    if(haloCellId == -1) {
      // If halo cell does not exist, get parent-level cell
      haloCellId = grid().tree().neighbor(grid().tree().parent(internalCellId), nghbrCellDir);
    }

    // Get child-level halo cells if halo cell has children
    if(grid().tree().hasChildren(haloCellId)) {
      haloCellId = m_surfaces.fineCellId(srfcId);
    }

    // Determine domain id
    MInt exchangeNghbrDomainId = -1;
    for(MInt i = 0; i < noDomains() + 1; i++) {
      if(grid().domainOffset(i) > grid().tree().globalId(haloCellId)) {
        exchangeNghbrDomainId = i - 1;
        break;
      }
    }
    ASSERT(exchangeNghbrDomainId > -1, "Could not find domain id!");

    // Add new exchange neighbor domain if not yet existing
    if(exchangeNghbrDomainMap.count(exchangeNghbrDomainId) == 0) {
      m_exchangeNghbrDomains[m_noExchangeNghbrDomains] = exchangeNghbrDomainId;
      exchangeNghbrDomainMap[exchangeNghbrDomainId] = m_noExchangeNghbrDomains;
      m_noExchangeNghbrDomains++;
    }

    // Add surface to container
    m_mpiSurfaces[exchangeNghbrDomainMap[exchangeNghbrDomainId]].push_back(srfcId);
  }

  // Reset size of containers to actual size
  m_exchangeNghbrDomains.resize(m_noExchangeNghbrDomains);
  m_mpiSurfaces.resize(m_noExchangeNghbrDomains);

  // Sort the MPI surfaces container according to the global surface id
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    sort(m_mpiSurfaces[i].begin(), m_mpiSurfaces[i].end(),
         [this](const MInt a, const MInt b) { return (m_surfaces.globalId(a) < m_surfaces.globalId(b)); });
  }

  // Build a container for receiving periodic interfaces
  m_mpiRecvSurfaces.resize(m_noExchangeNghbrDomains);
  m_mpiRecvSurfaces = m_mpiSurfaces;

  // Modification for periodic boundary conditions. This is required since in case of periodic
  // boundary conditions with periodicity on the same MPI rank, two surfaces exist with the same
  // global surface id.
  if(grid().periodicCartesianDir(0) || grid().periodicCartesianDir(1)
     || (nDim == 3 && grid().periodicCartesianDir(2))) {
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      // If a global surface id exists twice within one domain, change its order for receiving
      for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size() - 1; j++) {
        if(m_surfaces.globalId(m_mpiSurfaces[i][j]) == m_surfaces.globalId(m_mpiSurfaces[i][j + 1])) {
          swap(m_mpiRecvSurfaces[i][j], m_mpiRecvSurfaces[i][j + 1]);
        }
      }
    }
  }


  // TODO labels:DG,toremove Remove this debugging output once multi-solver simulations are properly tested
  // TODO labels:DG add test to check if mpi surfaces match on exchange domains!
  // stringstream ss;
  // ss << solverId() << " " << domainId() << " " << m_noExchangeNghbrDomains << " exchange neighb:";
  // for (MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
  //  ss << " " << m_exchangeNghbrDomains[i];
  //}
  // ss << endl;
  // for (MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
  //  ss << solverId() << " " << domainId() << " " << m_mpiSurfaces[i].size()
  //     << " MPI surfaces for domain " << m_exchangeNghbrDomains[i] << ":";
  //  for (MInt j = 0; j < static_cast<MInt>(m_mpiSurfaces[i].size()); j++) {
  //    ss << " " << m_surfaces.globalId(m_mpiSurfaces[i][j]);
  //  }
  //  ss << endl;
  //}
  // ss << endl;
  // cout << ss.str() << endl;
  // m_log << ss.str() << endl;

  // IMPORTANT:
  // If anything in the following loops is changed, please check if
  // updateNodeVariables() needs to be changed as well, since it contains more
  // or less the same code.

#ifdef DG_USE_MPI_BUFFERS
  // Resize send/receive buffers
  m_sendBuffers.resize(m_noExchangeNghbrDomains);
  m_recvBuffers.resize(m_noExchangeNghbrDomains);
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MInt size = 0;
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      // Use max. noNodes to account for p-refined surfaces
      const MInt noNodes1D = m_maxNoNodes1D;
      const MInt noNodesXD = ipow(noNodes1D, nDim - 1);
      const MInt dataBlockSize = noNodesXD * SysEqn::noVars();
      size += dataBlockSize;
    }

    m_sendBuffers[i].resize(size);
    m_recvBuffers[i].resize(size);
  }
#endif
#ifdef DG_USE_MPI_DERIVED_TYPES
#error Does not work anymore - need to fix for p-refinement!
  // Create and commit MPI derived datatypes for send/recv operations
  m_sendTypes.resize(m_noExchangeNghbrDomains);
  m_recvTypes.resize(m_noExchangeNghbrDomains);
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    // Note: This method uses non-MAIA types for interfacing with an external
    //       library
    // Count specifies the number of solvers in the derived data type
    const MInt count = m_mpiSurfaces[i].size();

    // If count is zero, skip this domain (since there is no data to exchange)
    // and send the data type to a null type. This is used later to determine
    // that no data needs to be exchanged.
    if(count == 0) {
      m_sendTypes[i] = MPI_DATATYPE_NULL;
      m_recvTypes[i] = MPI_DATATYPE_NULL;
      continue;
    }

    vector<MInt> lengths(count);
    vector<MPI_Aint> displacementsSend(count);
    vector<MPI_Aint> displacementsRecv(count);

    // Fill lengths and displacements vectors
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      const MInt srfcId = m_mpiSurfaces[i][j];
      // Set solver length
      const MInt noNodes1D = m_surfaces.noNodes1D(srfcId);
      const MInt noNodesXD = ipow(noNodes1D, nDim - 1);
      const MInt dataBlockSize = noNodesXD * SysEqn::noVars();
      lengths[j] = dataBlockSize;

      // Set displacements
      const MInt internalSideId = m_surfaces.internalSideId(srfcId);
      MPI_Get_address(m_surfaces.variables(srfcId, internalSideId), &displacementsSend[j], AT_);
      MPI_Get_address(m_surfaces.variables(srfcId, 1 - internalSideId), &displacementsRecv[j], AT_);
    }

    // Create MPI data types (hindexed since we're using absolute addresses for
    // the displacement)
    MPI_Type_create_hindexed(count, &lengths[0], &displacementsSend[0], type_traits<MFloat>::mpiType(), &m_sendTypes[i],
                             AT_);
    MPI_Type_create_hindexed(count, &lengths[0], &displacementsRecv[0], type_traits<MFloat>::mpiType(), &m_recvTypes[i],
                             AT_);

    // Commit data types
    MPI_Type_commit(&m_sendTypes[i], AT_);
    MPI_Type_commit(&m_recvTypes[i], AT_);
  }
#endif

  // Create persistent requests
  m_sendRequests.resize(m_noExchangeNghbrDomains);
  m_recvRequests.resize(m_noExchangeNghbrDomains);
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
#ifdef DG_USE_MPI_BUFFERS
    MPI_Send_init(&m_sendBuffers[i][0], m_sendBuffers[i].size(), type_traits<MFloat>::mpiType(),
                  m_exchangeNghbrDomains[i], domainId(), mpiComm(), &m_sendRequests[i], AT_, "m_sendBuffers[i][0]");
    MPI_Recv_init(&m_recvBuffers[i][0], m_recvBuffers[i].size(), type_traits<MFloat>::mpiType(),
                  m_exchangeNghbrDomains[i], m_exchangeNghbrDomains[i], mpiComm(), &m_recvRequests[i], AT_,
                  "m_recvBuffers[i][0]");
#endif
#ifdef DG_USE_MPI_DERIVED_TYPES
    MPI_Send_init(MPI_BOTTOM, 1, m_sendTypes[i], m_exchangeNghbrDomains[i], domainId(), mpiComm(), &m_sendRequests[i],
                  AT_, "MPI_BOTTOM");
    MPI_Recv_init(MPI_BOTTOM, 1, m_recvTypes[i], m_exchangeNghbrDomains[i], m_exchangeNghbrDomains[i], mpiComm(),
                  &m_recvRequests[i], AT_, "MPI_BOTTOM");
#endif
  }

  m_isInitMpiExchange = true;
  m_log << "done" << endl;
}


/*
 * \brief Calculate and - if necessary - exchange all statistical information
 *        about the simulation itself (e.g. total number of active cells,
 *        elements, surfaces etc.)
 *
 * \author Michael Schlottke
 * \date   July 2013
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initSimulationStatistics() {
  TRACE();

  m_log << "Initializing simulation statistics... ";

  // Set or determine local statistics

  // Minimum & maximum grid refinement level
  m_statLocalMinLevel = grid().minLevel();
  m_statLocalMaxLevel = grid().maxLevel();

  // Total cell counts
  m_statLocalNoCells = grid().noCells();
  m_statLocalNoInternalCells = grid().noInternalCells();
  m_statLocalNoHaloCells = grid().noCells() - grid().noInternalCells();
  m_statLocalMaxNoCells = grid().raw().treeb().capacity();

  // Total element counts
  m_statLocalNoElements = m_elements.size();
  m_statLocalNoHElements = m_helements.size();

  // Total surface counts
  m_statLocalNoSurfaces = m_surfaces.size();
  m_statLocalNoBoundarySurfaces = m_noBoundarySurfaces;
  m_statLocalNoInnerSurfaces = m_noInnerSurfaces;
  m_statLocalNoMpiSurfaces = m_noMpiSurfaces;
  m_statLocalMaxNoSurfaces = m_maxNoSurfaces;

  // Active cell/DOF count
  m_statLocalNoActiveCells = m_elements.size();
  m_statLocalNoActiveDOFs =
      accumulate(&m_elements.noNodes1D(0), &m_elements.noNodes1D(0) + m_elements.size(), 0,
                 [](const MInt result, const MInt noNodes1D) { return result + ipow(noNodes1D, nDim); });

  // Count active DOF per polynomial degree (if there are multiple polyDegs)
  const MInt noPolyDegs = m_maxPolyDeg - m_minPolyDeg + 1;
  if(noPolyDegs > 1) {
    m_statLocalNoActiveDOFsPolyDeg.assign(noPolyDegs, 0);
    const MInt noElements = m_elements.size();
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt polyDeg = m_elements.polyDeg(elementId);
      const MInt noDOFs = m_elements.noNodesXD(elementId);
      m_statLocalNoActiveDOFsPolyDeg[polyDeg - m_minPolyDeg] += noDOFs;
    }
  }


  // Minimum & maximum polynomial degrees
  m_statLocalMinPolyDeg = *min_element(&m_elements.polyDeg(0), &m_elements.polyDeg(0) + m_elements.size());
  m_statLocalMaxPolyDeg = *max_element(&m_elements.polyDeg(0), &m_elements.polyDeg(0) + m_elements.size());

  // Form sum of polynomial degress for average calculation
  m_statLocalPolyDegSum = accumulate(&m_elements.polyDeg(0), &m_elements.polyDeg(0) + m_elements.size(), 0);

  // Number of h- and p-refined surfaces
  m_statLocalNoHrefSurfs = count_if(&m_surfaces.fineCellId(0), &m_surfaces.fineCellId(0) + m_surfaces.size(),
                                    [](const MInt id) { return id != -1; });
  m_statLocalNoPrefSurfs = 0;
  // Only counted for inner surfaces
  // TODO labels:DG Also count for MPI surfaces
  for(MInt srfcId = m_innerSurfacesOffset; srfcId < m_innerSurfacesOffset + m_noInnerSurfaces; srfcId++) {
    const MInt elementIdL = m_surfaces.nghbrElementIds(srfcId, 0);
    const MInt elementIdR = m_surfaces.nghbrElementIds(srfcId, 1);
    if(m_elements.polyDeg(elementIdL) != m_elements.polyDeg(elementIdR)) {
      m_statLocalNoPrefSurfs++;
    }
  }

  // Determine global statistics

  // Pack buffer
  MLong buffer[26];
  // Min
  buffer[0] = m_statLocalMinLevel;
  buffer[1] = m_statLocalMinPolyDeg;
  // Max
  buffer[2] = m_statLocalMaxLevel;
  buffer[3] = m_statLocalMaxPolyDeg;
  // Sum
  buffer[4] = m_statLocalNoCells;
  buffer[5] = m_statLocalNoInternalCells;
  buffer[6] = m_statLocalNoHaloCells;
  buffer[7] = m_statLocalMaxNoCells;
  buffer[8] = m_statLocalNoElements;
  buffer[9] = m_statLocalNoSurfaces;
  buffer[10] = m_statLocalNoBoundarySurfaces;
  buffer[11] = m_statLocalNoInnerSurfaces;
  buffer[12] = m_statLocalNoMpiSurfaces;
  buffer[13] = m_statLocalMaxNoSurfaces;
  buffer[14] = m_statLocalNoActiveCells;
  buffer[15] = m_statLocalNoActiveDOFs;
  buffer[16] = m_statLocalPolyDegSum;
  buffer[17] = m_statLocalNoHrefSurfs;
  buffer[18] = m_statLocalNoPrefSurfs;
  buffer[19] = m_statLocalNoHElements;
  buffer[20] = m_noCutOffBoundarySurfaces[0];
  buffer[21] = m_noCutOffBoundarySurfaces[1];
  buffer[22] = m_noCutOffBoundarySurfaces[2];
  buffer[23] = m_noCutOffBoundarySurfaces[3];
  IF_CONSTEXPR(nDim == 3) {
    buffer[24] = m_noCutOffBoundarySurfaces[4];
    buffer[25] = m_noCutOffBoundarySurfaces[5];
  }
  else {
    buffer[24] = 0;
    buffer[25] = 0;
  }

  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::MPI]);
  RECORD_TIMER_START(m_timers[Timers::MPIComm]);

  // Exchange statistics
  // Operation: min
  MPI_Allreduce(MPI_IN_PLACE, &buffer[0], 2, MPI_INT, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "buffer[0]");
  // Operation: max
  MPI_Allreduce(MPI_IN_PLACE, &buffer[2], 2, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "buffer[2]");
  // Operation: sum
  MPI_Allreduce(MPI_IN_PLACE, &buffer[4], 22, maia::type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_,
                "MPI_IN_PLACE", "buffer[4]");

  RECORD_TIMER_STOP(m_timers[Timers::MPIComm]);
  RECORD_TIMER_STOP(m_timers[Timers::MPI]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);

  // Unpack buffer
  m_statGlobalMinLevel = buffer[0];
  m_statGlobalMinPolyDeg = buffer[1];
  m_statGlobalMaxLevel = buffer[2];
  m_statGlobalMaxPolyDeg = buffer[3];
  m_statGlobalNoCells = buffer[4];
  m_statGlobalNoInternalCells = buffer[5];
  m_statGlobalNoHaloCells = buffer[6];
  m_statGlobalMaxNoCells = buffer[7];
  m_statGlobalNoElements = buffer[8];
  m_statGlobalNoSurfaces = buffer[9];
  m_statGlobalNoBoundarySurfaces = buffer[10];
  m_statGlobalNoInnerSurfaces = buffer[11];
  m_statGlobalNoMpiSurfaces = buffer[12];
  m_statGlobalMaxNoSurfaces = buffer[13];
  m_statGlobalNoActiveCells = buffer[14];
  m_statGlobalNoActiveDOFs = buffer[15];
  m_statGlobalPolyDegSum = buffer[16];
  m_statGlobalNoHrefSurfs = buffer[17];
  m_statGlobalNoPrefSurfs = buffer[18];
  m_statGlobalNoHElements = buffer[19];
  m_statGlobalNoCutOffBoundarySurfaces[0] = buffer[20];
  m_statGlobalNoCutOffBoundarySurfaces[1] = buffer[21];
  m_statGlobalNoCutOffBoundarySurfaces[2] = buffer[22];
  m_statGlobalNoCutOffBoundarySurfaces[3] = buffer[23];
  m_statGlobalNoCutOffBoundarySurfaces[4] = buffer[24];
  m_statGlobalNoCutOffBoundarySurfaces[5] = buffer[25];

  m_log << "done" << endl;
}


/**
 * \brief Clean up after a simulation run.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-31
 *
 * Clean up only stuff that cannot (or maybe is not) properly handled
 * in the destructor of the solver. So e.g. do not clean up every single
 * std::vector that was used as a member variable, but free file resources,
 * free MPI resources, destroy manually allocated memory etc.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::cleanUp() {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    RECORD_TIMER_STOP(m_timers[Timers::Run]);
    return;
  }

  RECORD_TIMER_START(m_timers[Timers::CleanUp]);

  // Abort if solver not initialized
  if(!m_isInitSolver) {
    TERMM(1, "Solver was not initialized.");
  }

  // Free persistent MPI requests
  if(hasMpiExchange()) {
    finalizeMpiExchange();
  }

  RECORD_TIMER_STOP(m_timers[Timers::CleanUp]);
  RECORD_TIMER_STOP(m_timers[Timers::Run]);
}


/// \brief Cancel open MPI (receive) requests
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2018-10-17
///
/// Cancel opened receive request that do not have a matching send initiated yet.
/// Note: canceling send requests might cause MPI errors depending on the MPI implementation.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::cancelMpiRequests() {
  TRACE();

  // Complete already started communication in split-MPI mode
  if(g_splitMpiComm && m_mpiSendRequestsOpen) {
    RECORD_TIMER_START(m_timers[Timers::MainLoop]);
    RECORD_TIMER_START(m_timers[Timers::RungeKuttaStep]);
    RECORD_TIMER_START(m_timers[Timers::TimeDeriv]);
    finishMpiSurfaceExchange();
    RECORD_TIMER_STOP(m_timers[Timers::TimeDeriv]);
    RECORD_TIMER_STOP(m_timers[Timers::RungeKuttaStep]);
    RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);

    // Wait for send requests to finish
    MPI_Waitall(m_noExchangeNghbrDomains, &m_sendRequests[0], MPI_STATUSES_IGNORE, AT_);
    m_mpiSendRequestsOpen = false;
  }

  // Cancel opened receive requests
  if(m_mpiRecvRequestsOpen) {
    std::vector<MBool> waitForCancel(m_noExchangeNghbrDomains, false);
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      if(m_recvRequests[i] != MPI_REQUEST_NULL) {
        MPI_Cancel(&m_recvRequests[i], AT_);
        waitForCancel[i] = true;
      }
    }
    // Wait for all requests until they are canceled
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      if(waitForCancel[i]) {
        MPI_Wait(&m_recvRequests[i], MPI_STATUS_IGNORE, AT_);
      }
    }
    m_mpiRecvRequestsOpen = false;
  }
}


/*
 * \brief Free resources that were allocated in initMpiExchange().
 *
 * \author Michael Schlottke
 * \date   July 2013
 *
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::finalizeMpiExchange() {
  TRACE();

  if(m_noExchangeNghbrDomains > 0) {
    // Wait for open send requests to finish
    if(m_mpiSendRequestsOpen) {
      MPI_Waitall(m_noExchangeNghbrDomains, &m_sendRequests[0], MPI_STATUSES_IGNORE, AT_);
    }

    // Cancel and free requests if they are non-null
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      if(m_sendRequests[i] != MPI_REQUEST_NULL) {
        // Note: canceling send requests might cause MPI errors depending on the MPI
        // implementation https://github.com/mpi-forum/mpi-forum-historic/issues/479
        /* MPI_Cancel(&m_sendRequests[i], AT_); */
        MPI_Request_free(&m_sendRequests[i], AT_);
      }

      if(m_recvRequests[i] != MPI_REQUEST_NULL) {
        // Cancel opened receive request that do not have a matching send initiated yet
        if(m_mpiRecvRequestsOpen) {
          MPI_Cancel(&m_recvRequests[i], AT_);
        }
        MPI_Request_free(&m_recvRequests[i], AT_);
      }
    }
  }

  m_mpiSendRequestsOpen = false;
  m_mpiRecvRequestsOpen = false;

#ifdef DG_USE_MPI_DERIVED_TYPES
  // Free derived data types
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MPI_Type_free(&m_sendTypes[i], AT_);
    MPI_Type_free(&m_recvTypes[i], AT_);
  }
#endif
}


/**
 * \brief Set the initial condition in all elements
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-03
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initialCondition() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitialCondition]);

  m_log << "Applying initial conditions... ";

  const MFloat time = m_startTime;

  // First apply initial conditions from coupling
  // Note: changed with unified run loop, if you need to overwrite an initial condition set by the
  // coupling initialCondition() needs to be called again from the coupler!

  // Iterate over all elements
  const MInt noElements = m_elements.size();
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt noNodes1D = m_elements.noNodes1D(elementId);
    const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
    // Set node vars to minimum 1 to avoid errors in Tensor
    // TODO labels:DG,totest Check if this is really sensible
    MFloatTensor nodeVars(&m_elements.nodeVars(elementId), noNodes1D, noNodes1D, noNodes1D3,
                          max(SysEqn::noNodeVars(), 1));
    MFloatTensor u(&m_elements.variables(elementId), noNodes1D, noNodes1D, noNodes1D3, SysEqn::noVars());
    MFloatTensor x(&m_elements.nodeCoordinates(elementId), noNodes1D, noNodes1D, noNodes1D3, nDim);
    // Loop over all nodes and apply initial condition at all integration points
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          // Apply initial condition
          m_sysEqn.calcInitialCondition(time, &x(i, j, k, 0), &nodeVars(i, j, k, 0), &u(i, j, k, 0));
        }
      }
    }
  }

  m_log << "done" << endl;

  RECORD_TIMER_STOP(m_timers[Timers::InitialCondition]);
}


/// \brief Update all node variables at the surfaces.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-09-24
///
/// This method executes the following three steps:
/// - prolong all node variables from the elements to the surfaces
/// - apply the forward projection algorithm to the node variables
/// - exchange the node variable information for MPI surfaces
///
/// This method only needs to be called once after a change to the node
/// variables.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::updateNodeVariables() {
  TRACE();

  m_log << "Update node variable data on surfaces... ";

  const MInt* const polyDegs = &m_elements.polyDeg(0);
  const MInt* const noNodes = &m_elements.noNodes1D(0);
  const MInt* const surfaceIds = &m_elements.surfaceIds(0, 0);
  const MInt noElements = m_elements.size();

  //////////////////////////////////////////////////////////////////////////////
  // Prolong to surfaces
  //////////////////////////////////////////////////////////////////////////////

  // IMPORTANT:
  // If anything in this part of the method is changed, please check if
  // prolongToSurfaces() needs to be changed as well, since it contains more
  // or less the same code.
  using namespace dg::interpolation;

  if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Loop over all elements in given range
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt surfaceIdOffset = elementId * 2 * nDim;
      const MInt polyDeg = polyDegs[elementId];
      const MInt noNodes1D = noNodes[elementId];
      const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];

      // Extrapolate the solution to each surface on the faces
      for(MInt dir = 0; dir < 2 * nDim; dir++) {
        const MInt srfcId = surfaceIds[surfaceIdOffset + dir];
        const MInt side = 1 - dir % 2;

        MFloat* src = &m_elements.nodeVars(elementId);
        MFloat* dest = &m_surfaces.nodeVars(srfcId, side);
        prolongToFaceGauss<nDim, SysEqn::noNodeVars()>(src, dir, noNodes1D, &interp.m_LFace[0][0],
                                                       &interp.m_LFace[1][0], dest);
      }
    }
  } else if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Loop over all elements in given range
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt surfaceIdOffset = elementId * 2 * nDim;
      const MInt noNodes1D = noNodes[elementId];

      // Extrapolate the solution to each surface on the faces
      for(MInt dir = 0; dir < 2 * nDim; dir++) {
        const MInt srfcId = surfaceIds[surfaceIdOffset + dir];
        const MInt side = 1 - dir % 2;

        MFloat* src = &m_elements.nodeVars(elementId);
        MFloat* dest = &m_surfaces.nodeVars(srfcId, side);
        prolongToFaceGaussLobatto<nDim, SysEqn::noNodeVars()>(src, dir, noNodes1D, dest);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Forward projection for hp-refined surfaces
  //////////////////////////////////////////////////////////////////////////////

  // IMPORTANT:
  // If anything in this part of the method is changed, please check if
  // applyForwardProjection() needs to be changed as well, since it contains
  // more  or less the same code.

  const MInt noVars = SysEqn::noNodeVars();
  const MInt maxNoNodes1D = m_maxNoNodes1D;
  const MInt maxNoNodes1D3 = (nDim == 3) ? maxNoNodes1D : 1;
  const MInt noHElements = m_helements.size();
  const MInt noDirs = 2 * nDim;
  const MInt noSurfs = 2 * (nDim - 1);

  // Copy prolong results to h-refined surfaces (the prolong step stored its
  // results only in the first refined surface)
  if(noHElements > 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt hElementId = 0; hElementId < noHElements; hElementId++) {
      const MInt coarseElementId = m_helements.elementId(hElementId);
      for(MInt dir = 0; dir < noDirs; dir++) {
        const MInt coarseSrfcId = m_elements.surfaceIds(coarseElementId, dir);
        for(MInt pos = 0; pos < noSurfs; pos++) {
          const MInt hSrfcId = m_helements.hrefSurfaceIds(hElementId, dir, pos);
          const MInt side = 1 - dir % 2;

          // Skip if this surface does not exist
          if(hSrfcId == -1) {
            continue;
          }

          // Copy values of prolong step to all remaining refined surfaces
          if(coarseSrfcId != hSrfcId) {
            const MInt noNodes1D = noNodes[coarseElementId];
            const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
            const MInt noNodesXD = noNodes1D * noNodes1D3;

            // Copy nodeVars from prolong step to surface
            copy_n(&m_surfaces.nodeVars(coarseSrfcId, side),
                   noNodesXD * SysEqn::noNodeVars(),
                   &m_surfaces.nodeVars(hSrfcId, side));
          }
        }
      }
    }
  }

  // Apply mortar projection
  MFloatTensor projected(maxNoNodes1D, maxNoNodes1D3, noVars);

  // Apply projection to h- (and possibly p-)refined surfaces
  if(noHElements > 0) {
#ifdef _OPENMP
#pragma omp parallel for firstprivate(projected)
#endif
    for(MInt hElementId = 0; hElementId < noHElements; hElementId++) {
      for(MInt dir = 0; dir < noDirs; dir++) {
        for(MInt pos = 0; pos < noSurfs; pos++) {
          const MInt hSrfcId = m_helements.hrefSurfaceIds(hElementId, dir, pos);
          const MInt side = 1 - dir % 2;

          // Skip if this surface does not exist
          if(hSrfcId == -1) {
            continue;
          }

          const MInt surfaceNoNodes1D = m_surfaces.noNodes1D(hSrfcId);
          const MInt surfaceNoNodes1D3 = (nDim == 3) ? surfaceNoNodes1D : 1;
          const MInt size = surfaceNoNodes1D * surfaceNoNodes1D3 * noVars;

          // Project the coarse side to the surface
          calcMortarProjection<dg::mortar::forward, SysEqn::noNodeVars()>(
              hSrfcId, dir, &m_surfaces.nodeVars(hSrfcId, side), &projected[0], m_elements, m_surfaces);

          // Copy results of projection back to surface
          copy_n(&projected[0], size, &m_surfaces.nodeVars(hSrfcId, side));
        }
      }
    }
  }

  // Apply projection to pure p-refined surfaces
#ifdef _OPENMP
#pragma omp parallel for firstprivate(projected)
#endif
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt surfaceIdOffset = elementId * 2 * nDim;

    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      // Define auxiliary variables for better readability
      const MInt srfcId = surfaceIds[surfaceIdOffset + dir];

      // Skip boundary surfaces as they are *always* conforming
      if(srfcId < m_innerSurfacesOffset) {
        continue;
      }

      const MInt surfacePolyDeg = m_surfaces.polyDeg(srfcId);
      const MInt elementPolyDeg = m_elements.polyDeg(elementId);
      const MInt side = 1 - dir % 2;
      const MInt surfaceNoNodes1D = m_surfaces.noNodes1D(srfcId);
      const MInt surfaceNoNodes1D3 = (nDim == 3) ? surfaceNoNodes1D : 1;
      const MInt size = surfaceNoNodes1D * surfaceNoNodes1D3 * noVars;

      // Skip h-refined surfaces since they have been already projected
      if(m_surfaces.fineCellId(srfcId) != -1 && m_surfaces.fineCellId(srfcId) != m_elements.cellId(elementId)) {
        continue;
      }

      // Calculate forward projection for lower polyDeg elements
      if(surfacePolyDeg > elementPolyDeg) {
        // Calculate forward projection
        calcMortarProjection<dg::mortar::forward, SysEqn::noNodeVars()>(srfcId, dir, &m_surfaces.nodeVars(srfcId, side),
                                                                        &projected[0], m_elements, m_surfaces);

        // Copy results of projection back to surface
        copy_n(&projected[0], size, &m_surfaces.nodeVars(srfcId, side));
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // MPI exchange
  //////////////////////////////////////////////////////////////////////////////

  // Exit here if there is nothing to communicate
  if(!hasMpiExchange()) {
    return;
  }

  // IMPORTANT:
  // If anything in this part of the method is changed, please check if
  // initMpiExchange()/startMpiSurfaceExchange()/finishMpiSurfaceExchange()
  // needs to be changed as well, since they contain more or less the same code.

  // Set up buffers and requests
  vector<vector<MFloat>> sendBuffers(m_noExchangeNghbrDomains);
  vector<vector<MFloat>> recvBuffers(m_noExchangeNghbrDomains);
  const MInt dataBlockSize = ipow(m_maxNoNodes1D, nDim - 1) * SysEqn::noNodeVars();
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    const MInt size = m_mpiSurfaces[i].size() * dataBlockSize;
    sendBuffers[i].resize(size);
    recvBuffers[i].resize(size);
  }
  ScratchSpace<MPI_Request> sendRequests(m_noExchangeNghbrDomains, AT_, "sendRequests");
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  ScratchSpace<MPI_Request> recvRequests(m_noExchangeNghbrDomains, AT_, "recvRequests");
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);

  // Start receiving
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MPI_Irecv(&recvBuffers[i][0], recvBuffers[i].size(), type_traits<MFloat>::mpiType(), m_exchangeNghbrDomains[i],
              m_exchangeNghbrDomains[i], mpiComm(), &recvRequests[i], AT_, "recvBuffers[i][0]");
  }

  // Fill send buffers
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MInt size = 0;
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      const MInt srfcId = m_mpiSurfaces[i][j];
      const MInt sideId = m_surfaces.internalSideId(srfcId);

      // Copy nodeVars data
      const MFloat* data = &m_surfaces.nodeVars(srfcId, sideId);
      copy_n(data, dataBlockSize, &sendBuffers[i][size]);
      size += dataBlockSize;
    }
    ASSERT(size == static_cast<MInt>(sendBuffers[i].size()), "Data size does not match buffer size.");
  }

  // Start sending
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MPI_Isend(&sendBuffers[i][0], sendBuffers[i].size(), type_traits<MFloat>::mpiType(), m_exchangeNghbrDomains[i],
              domainId(), mpiComm(), &sendRequests[i], AT_, "sendBuffers[i][0]");
  }

  // Finish receiving
  MPI_Waitall(m_noExchangeNghbrDomains, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Unpack receive buffers
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MInt size = 0;
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      const MInt srfcId = m_mpiSurfaces[i][j];
      const MInt sideId = 1 - m_surfaces.internalSideId(srfcId);

      // Copy nodeVars data
      MFloat* data = &m_surfaces.nodeVars(srfcId, sideId);
      std::copy_n(&recvBuffers[i][size], dataBlockSize, data);
      size += dataBlockSize;
    }
    ASSERT(size == static_cast<MInt>(recvBuffers[i].size()), "Data size does not match buffer size.");
  }

  // Finish sending
  MPI_Waitall(m_noExchangeNghbrDomains, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  m_log << "done" << endl;
}


/// \brief Extends nodeVars from given planes to given directions.
///
/// \author Patrick Antony (patrick) <patrick@aia.rwth-aachen.de>
/// \date 2019-09-18
///
/// This method checks if the meanExtension feature has been
/// activated properly in the properties file and triggers the extension
/// for each specified direction/plane-pair.
///
/// This method must be called once after updateNodeVars and will
/// call updateNodeVars itself again.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::extendNodeVariables() {
  // check if this feature is activated properly in the .toml file
  MBool nodeVarsExtension = false;
  if(Context::propertyExists("nodeVarsExtension", m_solverId)) {
    nodeVarsExtension = Context::getSolverProperty<MBool>("nodeVarsExtension", m_solverId, AT_, &nodeVarsExtension);
  } else {
    const MBool extensionSetImplicit = Context::propertyExists("meanExtendDirections", m_solverId)
                                       || Context::propertyExists("meanExtendOffsets", m_solverId)
                                       || Context::propertyExists("meanExtendLimits", m_solverId);
    nodeVarsExtension = extensionSetImplicit;
  }

  // return if neither flag(explicit) nor parameters(implicit) are set
  if(!nodeVarsExtension) {
    return;
  }

  // read in directions and offsets (which specify planes)
  const MInt noExtendDirs = Context::propertyLength("meanExtendDirections", m_solverId);
  const MInt noExtendOffsets = Context::propertyLength("meanExtendOffsets", m_solverId);
  const MInt noExtendLimits = Context::propertyLength("meanExtendLimits", m_solverId);
  if(noExtendOffsets != noExtendDirs || noExtendLimits != noExtendDirs) {
    TERMM(1, "ERROR: # of specified directions (" + to_string(noExtendDirs) + "), # of offsets ("
                 + to_string(noExtendOffsets) + "), and # of limits (" + to_string(noExtendLimits) + ") do not match!");
  }

  std::vector<MInt> meanExtendDirections(noExtendDirs);
  std::vector<MFloat> meanExtendOffsets(noExtendDirs);
  std::vector<MFloat> meanExtendLimits(noExtendDirs);

  for(MInt i = 0; i < noExtendDirs; i++) {
    meanExtendDirections[i] = Context::getSolverProperty<MInt>("meanExtendDirections", m_solverId, AT_, i);
    meanExtendOffsets[i] = Context::getSolverProperty<MFloat>("meanExtendOffsets", m_solverId, AT_, i);
    meanExtendLimits[i] = Context::getSolverProperty<MFloat>("meanExtendLimits", m_solverId, AT_, i);
  }

  // trigger extension for each direction/offset-pair
  for(MInt i = 0; i < noExtendDirs; i++) {
    if(meanExtendDirections[i] >= 0 && meanExtendDirections[i] < 2 * nDim) {
      extendNodeVariablesSingleDirection(meanExtendDirections[i], meanExtendOffsets[i], meanExtendLimits[i]);
    }
  }
}


/// \brief Set nodeVars upstream of given plane to values of elements in plane.
///
/// \author Patrick Antony (patrick) <patrick@aia.rwth-aachen.de>
/// \date 2019-04-29
///
/// This method executes the following steps until it hits the domain border:
/// - extend nodeVars values in given direction orthogonal to given plane
/// --update nodeVars in elements and surfaces along the direction
/// -if a domain boundary is updated, updated surfaces are communicated and
///  the updated values are propagated on the next domain until no domain sends
///  updated values or updates its own elements.
/// afterwards updateNodeVars is called to update nodeVars on surfaces
/// orthogonal to the extend direction.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::extendNodeVariablesSingleDirection(const MInt extendDir,
                                                                         const MFloat extendOffset,
                                                                         const MFloat extendLimit) {
  m_log << "Propagating mean values in direction " << extendDir << " from offset " << extendOffset << "..."
        << std::endl;
  // set up and pre-compute grid and extension variables
  const MInt noElements = m_elements.size();
  const MInt* const noNodes1D = &m_elements.noNodes1D(0);
  const MInt extendSide = extendDir % 2;
  const MInt extendDimension = (extendDir - extendSide) / 2;
  const MInt donorDir = 2 * extendDimension + (1 - extendSide);
  const MInt noHElements = m_helements.size();
  const MInt noSurfs = 2 * (nDim - 1);

  // precompute map from surfaceId to surfaceId of h-refined child surfaces
  // for the surface, which is shared between the coarse element and one of the fine elements,
  // this gives for every position the surface of the child surface of that position. for child
  // surfaces, it gives number of surfaces+1, because they should not acces anything surface
  // related. for unrefined surfaces, all entries are -1
  // surfIdForward has the recv side's surfaces, because forward mortar projection is used there
  // to
  // get from coarse(1 element) to fine (multiple child elements)
  // surfIdReverse has the donor side's surfaces, because reverse mortar projection is used there
  // to
  // get from fine (multiple child elements) to coarse(1 element)
  MIntTensor hForwardSurfaceId(m_surfaces.size(), noSurfs);
  MIntTensor hReverseSurfaceId(m_surfaces.size(), noSurfs);
  hForwardSurfaceId.set(-1);
  hReverseSurfaceId.set(-1);
  for(MInt hElemId = 0; hElemId < noHElements; hElemId++) {
    const MInt elemId = m_helements.elementId(hElemId);
    // get global surface id for this surface
    const MInt surfIdForward = m_elements.surfaceIds(elemId, extendDir);
    const MInt surfIdReverse = m_elements.surfaceIds(elemId, donorDir);
    for(MInt pos = 0; pos < noSurfs; pos++) {
      hForwardSurfaceId(surfIdForward, pos) = m_helements.hrefSurfaceIds(hElemId, extendDir, pos);
      const MInt hReverseSurfaceIdFine = m_helements.hrefSurfaceIds(hElemId, donorDir, pos);
      hReverseSurfaceId(surfIdReverse, pos) = hReverseSurfaceIdFine;
      // this prevents access on a valid surface, but marks children surfs as refined
      if(hReverseSurfaceId(surfIdReverse, pos) != surfIdReverse && hReverseSurfaceIdFine > 0) {
        hReverseSurfaceId(hReverseSurfaceIdFine, 0) = m_surfaces.size() + 1;
      }
    }
  }
  // this is needed for hrefinement on mpisurfaces; for the lack of booltensor, these are ints
  MIntTensor surfaceHasRecvd(m_surfaces.size());
  surfaceHasRecvd.set(0.0);
  // mark those elements as updated, which are on the plane specified by the offset
  std::vector<MBool> elementUpdated(noElements, false);
  std::vector<MBool> elementOutsideLimit(noElements, false);
  std::vector<MInt> mpiSurfaceSendSize(m_surfaces.size(), 0);
  for(MInt elemId = 0; elemId < noElements; elemId++) {
    // assume here: all elements have either only old nodeVars or only new nodeVars
    const MFloat pSideCoordinate =
        m_surfaces.coords(m_elements.surfaceIds(elemId, extendDimension * 2 + 1), extendDimension);
    const MFloat mSideCoordinate =
        m_surfaces.coords(m_elements.surfaceIds(elemId, extendDimension * 2), extendDimension);
    elementUpdated[elemId] = (pSideCoordinate >= extendOffset && mSideCoordinate <= extendOffset);
    // mark elements that are outside the extend limit coordinates
    if((!extendSide && pSideCoordinate <= extendLimit && mSideCoordinate <= extendLimit)
       || (extendSide && pSideCoordinate >= extendLimit && mSideCoordinate >= extendLimit)) {
      elementOutsideLimit[elemId] = true;
    }
    // mark mpi surfaces with initially updated NodeVars as ready to send data to neighbor domain
    if(elementUpdated[elemId]) {
      const MInt sId = m_elements.surfaceIds(elemId, extendDir);
      if(sId >= m_mpiSurfacesOffset) {
        const MInt buffSize = ipow(noNodes1D[elemId], nDim - 1) * SysEqn::noNodeVars();
        mpiSurfaceSendSize[sId] = buffSize;
        if(hForwardSurfaceId(sId, 0) > -1) {
          for(MInt pos = 0; pos < noSurfs; pos++) {
            const MInt hElemId = m_surfaces.nghbrElementIds(sId, extendSide);
            mpiSurfaceSendSize[m_elements.surfaceIds(hElemId, extendDir)] = buffSize;
          }
        }
      }
    }
  }

  for(MInt extendIterator = 0; extendIterator < 4 * noDomains() + 1; extendIterator++) {
    if(extendIterator == 4 * noDomains()) {
      // terminate because strict upper bound for number if iterations is reached
      for(MInt elemId = 0; elemId < noElements; elemId++) {
        if(elementUpdated[elemId] == 0) {
          m_log << "r" << domainId() << " e" << elemId << " not updated" << std::endl;
        }
      }
      TERMM(1, "Maximum amount of extension iterations reached.");
    }
    for(MInt elementCounter = 0; elementCounter < noElements; elementCounter++) {
      MBool noMoreUpdates = true;
      for(MInt elemId = 0; elemId < noElements; elemId++) {
        // update neighbor element in extend direction, skipping elements that
        // have no neighbor in extend direction to update
        const MInt srfcId = m_elements.surfaceIds(elemId, extendDir);
        if(srfcId < m_innerSurfacesOffset) {
          continue;
        }
        if(srfcId < m_mpiSurfacesOffset) {
          const MInt updateeElementId = m_surfaces.nghbrElementIds(srfcId, extendSide);
          if(elementUpdated[elemId] && !(elementUpdated[updateeElementId])
             && !(elementOutsideLimit[updateeElementId])) {
            updateNodeVariablesSingleElement(updateeElementId, extendDir, hForwardSurfaceId, hReverseSurfaceId,
                                             elementUpdated, mpiSurfaceSendSize, noMoreUpdates);
            noMoreUpdates = false;
          }
        }
      }
      if(noMoreUpdates) {
        break;
      }
    }

    // MPI communication of surfaces updated in this iteration:
    MInt extendDomainId;
    MPI_Comm_rank(mpiComm(), &extendDomainId);
    // exchange information about whether any domain has updated nodeVars to communicate
    MInt thisRankIsDone = 1;
    for(MInt a = 0; a < m_surfaces.size(); a++) {
      if(mpiSurfaceSendSize[a] > 0) {
        thisRankIsDone = 0;
      }
    }
    MInt allRanksAreDone;
    MPI_Allreduce(&thisRankIsDone, &allRanksAreDone, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "thisRankIsDone",
                  "allRanksAreDone");
    if(allRanksAreDone) {
      break;
    }

    // exchange information on how much data to recv from each neighbor domain
    // allocate send and recv metadata buffer and fill send metadata buffer
    std::vector<std::vector<MInt>> nghbrSendSizes(m_noExchangeNghbrDomains);
    std::vector<std::vector<MInt>> nghbrRecvSizes(m_noExchangeNghbrDomains);

    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      nghbrSendSizes[i].resize(m_mpiSurfaces[i].size());
      nghbrRecvSizes[i].resize(m_mpiSurfaces[i].size());
      for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
        nghbrSendSizes[i][j] = mpiSurfaceSendSize[m_mpiSurfaces[i][j]];
        mpiSurfaceSendSize[m_mpiSurfaces[i][j]] = 0;
      }
    }

    // Exchange metadata buffers.
    std::vector<MPI_Request> metaDataRecvRequests(m_noExchangeNghbrDomains);
    std::vector<MPI_Request> metaDataSendRequests(m_noExchangeNghbrDomains);
    for(MInt nghbrIndex = 0; nghbrIndex < m_noExchangeNghbrDomains; nghbrIndex++) {
      MPI_Irecv(&nghbrRecvSizes[nghbrIndex][0], nghbrRecvSizes[nghbrIndex].size(), type_traits<MInt>::mpiType(),
                m_exchangeNghbrDomains[nghbrIndex], 0, mpiComm(), &metaDataRecvRequests[nghbrIndex], AT_,
                "nghbrRecvSizes[nghbrIndex][0]");
      MPI_Isend(&nghbrSendSizes[nghbrIndex][0], nghbrSendSizes[nghbrIndex].size(), type_traits<MInt>::mpiType(),
                m_exchangeNghbrDomains[nghbrIndex], 0, mpiComm(), &metaDataSendRequests[nghbrIndex], AT_,
                "nghbrSendSizes[nghbrIndex][0]");
    }
    MPI_Waitall(m_noExchangeNghbrDomains, &metaDataRecvRequests[0], MPI_STATUS_IGNORE, AT_);

    // Exchange actual nodevars to neighbordomains.
    // Set up buffers.
    vector<vector<MFloat>> nodeVarsSendBuffers(m_noExchangeNghbrDomains);
    vector<vector<MFloat>> nodeVarsRecvBuffers(m_noExchangeNghbrDomains);
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      const MInt sendSize = accumulate(nghbrSendSizes[i].begin(), nghbrSendSizes[i].end(), 0);
      const MInt recvSize = accumulate(nghbrRecvSizes[i].begin(), nghbrRecvSizes[i].end(), 0);
      nodeVarsSendBuffers[i].resize(sendSize);
      nodeVarsRecvBuffers[i].resize(recvSize);

      // Fill sendbuffer.
      MInt size = 0;
      for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
        if(nghbrSendSizes[i][j] < 1) {
          continue;
        }
        const MInt srfcId = m_mpiSurfaces[i][j];
        // The side of the surface from which the nodevars are copied should
        // not matter, as both sides have to contain updated nodeVars in an
        // appropriate basis.
        const MFloat* data = &m_surfaces.nodeVars(srfcId, 1 - extendSide);
        copy(data, data + nghbrSendSizes[i][j], &nodeVarsSendBuffers[i][size]);
        size += nghbrSendSizes[i][j];
      }
      ASSERT(size == static_cast<MInt>(nodeVarsSendBuffers[i].size()), "Data size does not match buffer size.");
    }
    std::vector<MPI_Request> nodeVarsRecvRequests(m_noExchangeNghbrDomains);
    std::vector<MPI_Request> nodeVarsSendRequests(m_noExchangeNghbrDomains);

    // Mpi-exchange buffer.
    for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
      if(nodeVarsRecvBuffers[i].size() > 0) {
        MPI_Irecv(&nodeVarsRecvBuffers[i][0], nodeVarsRecvBuffers[i].size(), type_traits<MFloat>::mpiType(),
                  m_exchangeNghbrDomains[i], m_exchangeNghbrDomains[i], mpiComm(), &nodeVarsRecvRequests[i], AT_,
                  "nodeVarsRecvBuffers[i][0]");
      } else {
        nodeVarsRecvRequests[i] = MPI_REQUEST_NULL;
      }
      if(nodeVarsSendBuffers[i].size() > 0) {
        MPI_Isend(&nodeVarsSendBuffers[i][0], nodeVarsSendBuffers[i].size(), type_traits<MFloat>::mpiType(),
                  m_exchangeNghbrDomains[i], domainId(), mpiComm(), &nodeVarsSendRequests[i], AT_,
                  "nodeVarsSendBuffers[i][0]");
      } else {
        nodeVarsSendRequests[i] = MPI_REQUEST_NULL;
      }
    }
    // Since the send buffer will not get modified after this point, this wait
    // can be moved to a later point, it is posted here for readability only.
    MPI_Waitall(m_noExchangeNghbrDomains, &nodeVarsSendRequests[0], MPI_STATUSES_IGNORE, AT_);

    // This wait can not be posted later than here, since the recvbuffer is
    // read after this!
    MPI_Waitall(m_noExchangeNghbrDomains, &nodeVarsRecvRequests[0], MPI_STATUSES_IGNORE, AT_);
    // Unpack recvbuffers into surfaces and trigger updates on element(s) at surfaces.
    for(MInt nghbrIndex = 0; nghbrIndex < m_noExchangeNghbrDomains; nghbrIndex++) {
      // Unpack recvbuffer.
      MInt size = 0;
      for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[nghbrIndex].size(); j++) {
        if(nghbrRecvSizes[nghbrIndex][j] < 1) {
          continue;
        }
        const MInt srfcId = m_mpiSurfaces[nghbrIndex][j];
        surfaceHasRecvd[srfcId] = 1.0;
        // Copy nodeVars data.
        copy_n(&nodeVarsRecvBuffers[nghbrIndex][size],
               nghbrRecvSizes[nghbrIndex][j],
               &m_surfaces.nodeVars(srfcId, 1 - extendSide));
        copy_n(&nodeVarsRecvBuffers[nghbrIndex][size],
               nghbrRecvSizes[nghbrIndex][j],
               &m_surfaces.nodeVars(srfcId, extendSide));
        // Set up update.
        const MInt recvElementId = m_surfaces.nghbrElementIds(srfcId, extendSide);
        // Coarse cells with fine neighbors on a different domain
        //("reverse h-Refined") update only if all child surfaces have received
        // updated nodeVars.
        MBool srfcIsReverseHrefTemp = false;
        for(MInt pos = 0; pos < noSurfs; pos++) {
          if(hReverseSurfaceId(srfcId, pos) > -1) {
            srfcIsReverseHrefTemp = true;
          }
        }
        const MBool srfcIsReverseHref = srfcIsReverseHrefTemp;
        MBool allSurfsReady = true;
        if(srfcIsReverseHref) {
          const MInt mainSurfId = m_elements.surfaceIds(recvElementId, donorDir);
          for(MInt pos = 0; pos < noSurfs; pos++) {
            if(surfaceHasRecvd[hReverseSurfaceId(mainSurfId, pos)] < 1) {
              allSurfsReady = false;
            }
          }
        }
        if(allSurfsReady) {
          MBool dummy; // this is not used here
          const MInt updateeElementId = m_surfaces.nghbrElementIds(srfcId, extendSide);
          // skip elements outside the extend limit coordinates
          if(!elementOutsideLimit[updateeElementId]) {
            updateNodeVariablesSingleElement(recvElementId, extendDir, hForwardSurfaceId, hReverseSurfaceId,
                                             elementUpdated, mpiSurfaceSendSize, dummy);
            elementUpdated[recvElementId] = true;
          }
        }
        size += nghbrRecvSizes[nghbrIndex][j];
      }
      ASSERT(size == static_cast<MInt>(nodeVarsRecvBuffers[nghbrIndex].size()),
             "Data size does not match buffer size.");
    }
    // The sendSurfaceBufferSize might be  modified in the next iteration,
    // so we have to wait for this request no later than here.
    MPI_Waitall(m_noExchangeNghbrDomains, &metaDataSendRequests[0], MPI_STATUS_IGNORE, AT_);
  }
  m_log << "done" << std::endl;
  updateNodeVariables();
}


/// \brief Sets nodeVars of an element to values on the surface opposite to extendDir.
///
/// \author Patrick Antony (patrick) <patrick@aia.rwth-aachen.de>
/// \date 2019-04-29
///
/// This method sets the nodeVars of an element to values on the surface opposite to extendDir.
/// It then updates the nodeVars on the surface in extendDir of the element,
/// and notes the surfaces buffer size if it is on a domain boundary.
/// If the element is h-refined, it also calls itself for the child elements.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::updateNodeVariablesSingleElement(const MInt elementId,
                                                                       const MInt extendDir,
                                                                       MIntTensor hForwardSurfaceId,
                                                                       MIntTensor hReverseSurfaceId,
                                                                       std::vector<MBool>& elementUpdated,
                                                                       std::vector<MInt>& mpiSurfaceSendSize,
                                                                       MBool& noMoreUpdates) {
  // Precompute variables.
  const MInt* const noNodes1D = &m_elements.noNodes1D(0);
  const MInt* const surfaceNoNodes1D = &m_surfaces.noNodes1D(0);
  const MInt extendSide = extendDir % 2;
  const MInt extendDimension = (extendDir - extendSide) / 2;
  const MInt donorDir = 2 * extendDimension + (1 - extendSide);
  const MInt noVars = max(SysEqn::noNodeVars(), 1);
  const MInt noSurfs = 2 * (nDim - 1);
  const MInt maxNoNodes1D = noNodes1D[elementId];
  const MInt maxNoNodes1D3 = (nDim == 3) ? maxNoNodes1D : 1;
  MInt extendDomainId;
  MPI_Comm_rank(mpiComm(), &extendDomainId);
  // Set up Ids and data tensors.
  const MInt donorSrfcId = m_elements.surfaceIds(elementId, donorDir);
  const MInt recvSrfcId = m_elements.surfaceIds(elementId, extendDir);
  MFloatTensor nodeVars(&m_elements.nodeVars(elementId), maxNoNodes1D, maxNoNodes1D, maxNoNodes1D3, noVars);
  const MInt noNodesSecondSurfDim = extendDimension == 2 ? maxNoNodes1D : maxNoNodes1D3;
  MFloatTensor elemNodeVarsSlice(maxNoNodes1D, noNodesSecondSurfDim, noVars);
  /* const MFloatTensor donorNodeVars(&m_surfaces.nodeVars(donorSrfcId, 1 - extendSide), */
  /*                                    maxNoNodes1D, noNodesSecondSurfDim, noVars); */
  MFloatTensor recvSurfaceNodeVarsMSide(&m_surfaces.nodeVars(recvSrfcId, extendSide), maxNoNodes1D,
                                        noNodesSecondSurfDim, noVars);
  MFloatTensor recvSurfaceNodeVarsPSide(&m_surfaces.nodeVars(recvSrfcId, 1 - extendSide), maxNoNodes1D,
                                        noNodesSecondSurfDim, noVars);

  // Update receiver cell and its receiver side surface.
  const MBool pRefinedDonorSide = m_surfaces.polyDeg(donorSrfcId) > m_elements.polyDeg(elementId);
  const MBool pRefinedRecvSide = m_surfaces.polyDeg(recvSrfcId) > m_elements.polyDeg(elementId);
  const MBool srfcIsForwardHref = hForwardSurfaceId(donorSrfcId, 0) > -1;

  // Get default node variables for those variables that are not extended
  MFloatScratchSpace defaultNodeVars(SysEqn::noNodeVars(), AT_, "defaultNodeVars");
  m_sysEqn.getDefaultNodeVars(defaultNodeVars.getPointer());
  // Store if variables needs to be extended
  MBoolScratchSpace extendNodeVar(SysEqn::noNodeVars(), AT_, "extendNodeVars");
  for(MInt iVars = 0; iVars < noVars; iVars++) {
    extendNodeVar[iVars] = m_sysEqn.extendNodeVar(iVars);
  }

  // This abuses the fact that the coarse elements surface and the main position's fine surface
  // share the srfcID
  MBool srfcIsReverseHrefTemp = false;
  for(MInt pos = 0; pos < noSurfs; pos++) {
    if(hReverseSurfaceId(donorSrfcId, pos) > -1) {
      srfcIsReverseHrefTemp = true;
    }
  }
  const MBool srfcIsReverseHref = srfcIsReverseHrefTemp;
  MBool srfcIsMainPosTemp = false;
  for(MInt pos = 0; pos < noSurfs; pos++) {
    if(hReverseSurfaceId(donorSrfcId, pos) == donorSrfcId) {
      srfcIsMainPosTemp = true;
    }
  }
  const MBool srfcIsMainPos = srfcIsMainPosTemp;
  if(srfcIsForwardHref) { // One donor element, multiple recv/updatee elements
    // Call update for children, then proceed to update self.
    for(MInt pos = 1; pos < noSurfs; pos++) {
      const MInt hSrfcId = hForwardSurfaceId(donorSrfcId, pos);
      const MInt hElemId = m_surfaces.nghbrElementIds(hSrfcId, extendSide);
      updateNodeVariablesSingleElement(hElemId, extendDir, hForwardSurfaceId, hReverseSurfaceId, elementUpdated,
                                       mpiSurfaceSendSize, noMoreUpdates);
    }
  }
  // h(p) is excluded here, because h already includes p refinement(if necessary)
  if(!srfcIsReverseHref) { // one donor/one recv
    // project/copy from donorVars(fine surface) to element(coarse).
    if(pRefinedDonorSide) { // p projection required
      const MInt maxNoNodes1DFine = surfaceNoNodes1D[donorSrfcId];
      const MInt maxNoNodes1D3Fine = nDim == 3 ? maxNoNodes1DFine : 1;
      const MInt noNodesSecondSurfDimFine = extendDimension == 2 ? maxNoNodes1DFine : maxNoNodes1D3Fine;
      // P-ref mortar projection projects onto same size buffers (high degree),
      // even though lower p'nomial degree has fewer dofs.
      MFloatTensor projectedRight(maxNoNodes1DFine, noNodesSecondSurfDimFine, noVars);
      calcMortarProjection<dg::mortar::reverse, SysEqn::noNodeVars()>(donorSrfcId, donorDir,
                                                                      &m_surfaces.nodeVars(donorSrfcId, 1 - extendSide),
                                                                      &projectedRight[0], m_elements, m_surfaces);
      for(MInt a = 0; a < maxNoNodes1D; ++a) {
        for(MInt b = 0; b < maxNoNodes1D3; ++b) {
          for(MInt iVars = 0; iVars < noVars; iVars++) {
            elemNodeVarsSlice(a, b, iVars) = projectedRight(a, b, iVars);
          }
        }
      }
    } else { // no p projection required
      copy_n(&m_surfaces.nodeVars(donorSrfcId, 1 - extendSide),
             maxNoNodes1D * noNodesSecondSurfDim * noVars,
             &elemNodeVarsSlice[0]);
    }
    // copy nodeVars to element
    for(MInt a = 0; a < maxNoNodes1D; ++a) {
      for(MInt b = 0; b < maxNoNodes1D; ++b) {
        for(MInt c = 0; c < maxNoNodes1D3; ++c) {
          // determine which iterators are need for the surface-like elemNodeVarsSlice
          const MInt iNodes2D = extendDimension == 0 ? b : a;
          const MInt iNodes3D = extendDimension == 2 ? b : c;
          for(MInt iVars = 0; iVars < noVars; iVars++) {
            nodeVars(a, b, c, iVars) =
                (extendNodeVar[iVars]) ? elemNodeVarsSlice(iNodes2D, iNodes3D, iVars) : defaultNodeVars(iVars);
          }
        }
      }
    }
    // Mark element as updated.
    elementUpdated[elementId] = true;
    noMoreUpdates = false;
    // Copy from elem(coarse) and map to recvSurf(fine)
    if(pRefinedRecvSide) {
      // Set up projection to next surface.
      // Skip h-refined surfaces, as they will be handled separately.
      if(hForwardSurfaceId(recvSrfcId, 0) <= 0) {
        const MInt maxNoNodes1DSurf = surfaceNoNodes1D[recvSrfcId];
        const MInt maxNoNodes1D3Surf = nDim == 3 ? maxNoNodes1DSurf : 1;
        const MInt noNodesSecondRefSurfDim = extendDimension == 2 ? maxNoNodes1DSurf : maxNoNodes1D3Surf;
        MFloatTensor projected(maxNoNodes1DSurf, noNodesSecondRefSurfDim, noVars);
        // Project data to next surface.
        calcMortarProjection<dg::mortar::forward, SysEqn::noNodeVars()>(
            recvSrfcId, extendDir, &m_surfaces.nodeVars(donorSrfcId, 1 - extendSide), &projected[0], m_elements,
            m_surfaces);
        // Update next surface(both sides).
        copy_n(&projected[0],
               maxNoNodes1DSurf * noNodesSecondRefSurfDim * noVars,
               &m_surfaces.nodeVars(recvSrfcId, extendSide));
        copy_n(&projected[0],
               maxNoNodes1DSurf * noNodesSecondRefSurfDim * noVars,
               &m_surfaces.nodeVars(recvSrfcId, 1 - extendSide));
        // Mark recv surface for mpi exchange if on domain border.
        if(recvSrfcId >= m_mpiSurfacesOffset) {
          mpiSurfaceSendSize[recvSrfcId] = maxNoNodes1DSurf * maxNoNodes1D3Surf * SysEqn::noNodeVars();
        }
      }
    } else { // no projection needed, just copy to both sides of recv surface
      copy_n(&elemNodeVarsSlice[0], maxNoNodes1D * noNodesSecondSurfDim * noVars, &recvSurfaceNodeVarsMSide[0]);
      copy_n(&elemNodeVarsSlice[0], maxNoNodes1D * noNodesSecondSurfDim * noVars, &recvSurfaceNodeVarsPSide[0]);
      // Mark recv surface for mpi exchange if on domain border.
      if(recvSrfcId >= m_mpiSurfacesOffset) {
        mpiSurfaceSendSize[recvSrfcId] = maxNoNodes1D * maxNoNodes1D3 * SysEqn::noNodeVars();
      }
    }
  }
  if(srfcIsReverseHref && srfcIsMainPos) { // multiple donors, one receive element
    // Check if all four pos have updated means.
    MBool allDonorsHaveUpdatedMeans = true;
    for(MInt pos = 0; pos < noSurfs; pos++) {
      const MInt hElemId = m_surfaces.nghbrElementIds(hReverseSurfaceId(donorSrfcId, pos), 1 - extendSide);
      // This case happens if href happens across domain boundaries and is handled outside of this
      // function.
      if(hElemId < 0) {
        allDonorsHaveUpdatedMeans = true;
        continue;
      }
      if(elementUpdated[hElemId] == false) {
        allDonorsHaveUpdatedMeans = false;
        elementUpdated[elementId] = false;
        noMoreUpdates = false;
      }
    }
    // Only when all donor surfaces are updated, the recv element assembles its nodevars.
    if(allDonorsHaveUpdatedMeans) {
      const MInt maxNoNodes1DSurf = surfaceNoNodes1D[donorSrfcId];
      const MInt maxNoNodes1D3Surf = nDim == 3 ? maxNoNodes1DSurf : 1;
      const MInt noNodesSecondRefSurfDim = extendDimension == 2 ? maxNoNodes1DSurf : maxNoNodes1D3Surf;
      const MInt coarseNextSrfcId = m_elements.surfaceIds(elementId, extendDir);
      // mark refined cell as updated with hElemId
      elementUpdated[elementId] = true;
      // mark recv surface for mpi exchange if on domain border
      if(coarseNextSrfcId >= m_mpiSurfacesOffset) {
        mpiSurfaceSendSize[coarseNextSrfcId] = maxNoNodes1D * maxNoNodes1D3 * SysEqn::noNodeVars();
      }
      noMoreUpdates = false;
      MFloatTensor projectedSum(maxNoNodes1DSurf, noNodesSecondRefSurfDim, noVars);
      projectedSum.set(0.0);
      for(MInt pos = 0; pos < noSurfs; pos++) {
        const MInt hSrfId = hReverseSurfaceId(donorSrfcId, pos);
        MFloatTensor projected(maxNoNodes1DSurf, noNodesSecondRefSurfDim, noVars);
        projected.set(0.0);
        calcMortarProjection<dg::mortar::reverse, SysEqn::noNodeVars()>(
            hSrfId, donorDir, &m_surfaces.nodeVars(hSrfId, 1 - extendSide), &projected[0], m_elements, m_surfaces);
        for(MInt i = 0; i < maxNoNodes1D; i++) {
          for(MInt j = 0; j < noNodesSecondSurfDim; j++) {
            for(MInt k = 0; k < noVars; k++) {
              projectedSum(i, j, k) += projected(i, j, k);
            }
          }
        }
      }
      // Copy projected values to element.
      // This loop structure is suboptimal regarding memory acces, but easier to read.
      // If better performance is desired here, consider using separate loops for each possible
      // direction.
      MFloatTensor hNodeVarsLeft(&m_elements.nodeVars(elementId), maxNoNodes1D, maxNoNodes1D, maxNoNodes1D3, noVars);
      MFloatTensor hNodeVarsLeftSurf1(&m_surfaces.nodeVars(coarseNextSrfcId, extendSide), maxNoNodes1D,
                                      noNodesSecondSurfDim, noVars);
      MFloatTensor hNodeVarsLeftSurf2(&m_surfaces.nodeVars(coarseNextSrfcId, 1 - extendSide), maxNoNodes1D,
                                      noNodesSecondSurfDim, noVars);
      for(MInt a = 0; a < maxNoNodes1D; ++a) {
        for(MInt b = 0; b < maxNoNodes1D; ++b) {
          for(MInt c = 0; c < maxNoNodes1D3; ++c) {
            const MInt iNodes2D = extendDimension == 0 ? b : a;
            const MInt iNodes3D = extendDimension == 2 ? b : c;
            for(MInt iVars = 0; iVars < noVars; iVars++) {
              if(extendNodeVar[iVars]) {
                hNodeVarsLeft(a, b, c, iVars) = projectedSum(iNodes2D, iNodes3D, iVars);
                elemNodeVarsSlice(iNodes2D, iNodes3D, iVars) = projectedSum(iNodes2D, iNodes3D, iVars);
                hNodeVarsLeftSurf1(iNodes2D, iNodes3D, iVars) = projectedSum(iNodes2D, iNodes3D, iVars);
                hNodeVarsLeftSurf2(iNodes2D, iNodes3D, iVars) = projectedSum(iNodes2D, iNodes3D, iVars);
              } else {
                hNodeVarsLeft(a, b, c, iVars) = defaultNodeVars(iVars);
                elemNodeVarsSlice(iNodes2D, iNodes3D, iVars) = defaultNodeVars(iVars);
                hNodeVarsLeftSurf1(iNodes2D, iNodes3D, iVars) = defaultNodeVars(iVars);
                hNodeVarsLeftSurf2(iNodes2D, iNodes3D, iVars) = defaultNodeVars(iVars);
              }
            }
          }
        }
      }
      elementUpdated[elementId] = true;
      noMoreUpdates = false;
    }
  }
  // If the receiver surface is h(p)-refined(because the element after the receive element is),
  //  copy the updated nodeVars to all child surfaces and project.
  //  This guarantees that the next element can use these as is.
  if(hForwardSurfaceId(recvSrfcId, 0) > 0) {
    const MInt maxNoNodes1DFine = surfaceNoNodes1D[recvSrfcId];
    const MInt maxNoNodes1D3Fine = nDim == 3 ? maxNoNodes1DFine : 1;
    const MInt noNodesSecondFineSurfDim = extendDimension == 2 ? maxNoNodes1DFine : maxNoNodes1D3Fine;
    for(MInt pos = 0; pos < noSurfs; pos++) {
      const MInt hSrfId = hForwardSurfaceId(recvSrfcId, pos);
      copy_n(&elemNodeVarsSlice(0, 0, 0),
             maxNoNodes1D * noNodesSecondSurfDim * noVars,
             &m_surfaces.nodeVars(hSrfId, 1 - extendSide));
    }
    for(MInt pos = 0; pos < noSurfs; pos++) {
      const MInt hSrfId = hForwardSurfaceId(recvSrfcId, pos);
      MFloatTensor projected(maxNoNodes1DFine, noNodesSecondFineSurfDim, noVars);
      projected.set(0.0);
      calcMortarProjection<dg::mortar::forward, SysEqn::noNodeVars()>(
          hSrfId, extendDir, &m_surfaces.nodeVars(hSrfId, 1 - extendSide), &projected[0], m_elements, m_surfaces);
      MFloatTensor nodeVarsNextSrfPSide(&m_surfaces.nodeVars(hSrfId, extendSide), maxNoNodes1DFine,
                                        noNodesSecondFineSurfDim, noVars);
      MFloatTensor nodeVarsNextSrfMSide(&m_surfaces.nodeVars(hSrfId, 1 - extendSide), maxNoNodes1DFine,
                                        noNodesSecondFineSurfDim, noVars);
      for(MInt iNodes2D = 0; iNodes2D < maxNoNodes1DFine; iNodes2D++) {
        for(MInt iNodes3D = 0; iNodes3D < noNodesSecondFineSurfDim; iNodes3D++) {
          for(MInt iVars = 0; iVars < noVars; iVars++) {
            nodeVarsNextSrfPSide(iNodes2D, iNodes3D, iVars) = projected(iNodes2D, iNodes3D, iVars);
            nodeVarsNextSrfMSide(iNodes2D, iNodes3D, iVars) = projected(iNodes2D, iNodes3D, iVars);
          }
        }
      }
      // Mark next left surface of child element for mpi exchange if on domain border.
      if(recvSrfcId >= m_mpiSurfacesOffset) {
        mpiSurfaceSendSize[hSrfId] = maxNoNodes1DFine * maxNoNodes1D3Fine * noVars;
      }
    }
  }
}


/**
 * \brief Print initialization summary to user and m_log
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-29
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::outputInitSummary() {
  TRACE();

  // Get proper names from numeric properties
  const MString polynomialType = "Legendre";
  const MString integrationMethod = (m_dgIntegrationMethod == 0) ? "Gauss" : "Gauss-Lobatto";

  using namespace maia::logtable;

  // INIT SUMMARY FRAME
  Frame summary(" SOLVER " + std::to_string(solverId()) + " INITIALIZATION SUMMARY AT TIME STEP "
                + std::to_string(m_timeStep));

  // PROBLEM SUMMARY GROUP
  Group& problem = summary.addGroup(" PROBLEM SUMMARY");
  problem.addData("System of equations", m_sysEqn.sysEqnName());
  problem.addData("Number of dimensions", nDim);
  Data& vars = problem.addData("Number of variables", m_sysEqn.noVars());
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    if(i == 0) {
      vars.addData("Conservative variable name(s)", m_sysEqn.consVarNames(i));
    } else {
      vars.addData("", m_sysEqn.consVarNames(i));
    }
  }
  problem.addData("Restart", m_restart);
  if(m_restart) {
    problem.addData("Initial condition", getRestartFileName(m_timeStep, m_useNonSpecifiedRestartFile));
  } else {
    problem.addData("Initial condition", m_sysEqn.m_initialCondition);
  }
  problem.addData("Start time (non-dimensionalized)", m_time);
  problem.addData("Final time (non-dimensionalized)", m_finalTime);
  problem.addData("CFL", m_sysEqn.cfl());
  problem.addData("Recalculation interval for time step", m_calcTimeStepInterval);
  problem.addData("Time step", m_timeStep);
  problem.addData("Maximum number of time steps", m_timeSteps);

  // DISCRETIZATION SUMMARY GROUP
  Group& discret = summary.addGroup("DISCRETIZATION SUMMARY");
  discret.addData("Initial polynomial degree", m_initPolyDeg);
  discret.addData("Minimum polynomial degree (limit)", m_minPolyDeg);
  discret.addData("Maximum polynomial degree (limit)", m_maxPolyDeg);
  discret.addData("Polynomial type", polynomialType);
  discret.addData("Integration method", integrationMethod);
  MString timeIntegrationScheme;
  switch(m_dgTimeIntegrationScheme) {
    case DG_TIMEINTEGRATION_CARPENTER_4_5: {
      default:
        timeIntegrationScheme = "Carpenter 4/5";
        break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_4_8: {
      timeIntegrationScheme = "Toulorge C 4/8";
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_14: {
      timeIntegrationScheme = "Niegemann 4/14";
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_13: {
      timeIntegrationScheme = "Niegemann 4/13";
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_3_7: {
      timeIntegrationScheme = "Toulorge C 3/7";
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEF_4_8: {
      timeIntegrationScheme = "Toulorge F 4/8";
      break;
    }
  }
  discret.addData("Time integration scheme", timeIntegrationScheme);

  // PARALLELIZATION SUMMARY
  Group& parallel = summary.addGroup("PARALLELIZATION SUMMARY");
  parallel.addData("Domain id", domainId());
  parallel.addData("Number of neighbor domains", grid().noNeighborDomains());
  parallel.addData("Number of exchange neighbor domains", m_noExchangeNghbrDomains);
  parallel.addData("Total number of domains", noDomains());

  // GRID SUMMARY LOCAL
  Group& gridLocal = summary.addGroup("GRID SUMMARY (LOCAL)");
  gridLocal.addData("Minimum used polynomial degree", m_statLocalMinPolyDeg);
  gridLocal.addData("Maximum used polynomial degree", m_statLocalMaxPolyDeg);
  gridLocal.addData("Average used polynomial degree", 1.0 * m_statLocalPolyDegSum / m_statLocalNoActiveCells);
  gridLocal.addBlank();
  gridLocal.addData("Minimum grid level", m_statLocalMinLevel);
  gridLocal.addData("Maximum grid level", m_statLocalMaxLevel);
  gridLocal.addBlank();
  Data& noCellsLocal = gridLocal.addData("Number of cells", m_statLocalNoCells);
  noCellsLocal.addData("internal cells", m_statLocalNoInternalCells);
  noCellsLocal.addData("halo cells", m_statLocalNoHaloCells);
  gridLocal.addData("Number of active cells", m_statLocalNoActiveCells);
  gridLocal.addData("Number of active DOFs", m_statLocalNoActiveDOFs);
  // Number of active DOFs per polynomial degree
  const MInt noPolyDegs = m_maxPolyDeg - m_minPolyDeg + 1;
  if(noPolyDegs > 1) {
    for(MInt i = 0; i < noPolyDegs; i++) {
      MString name = "Number of active DOFs polyDeg=" + std::to_string(m_minPolyDeg + i);
      gridLocal.addData(name, m_statLocalNoActiveDOFsPolyDeg[i]);
    }
  }
  gridLocal.addData("Max. number of cells (collector size)", m_statLocalMaxNoCells);
  gridLocal.addData("Memory utilization", 100.0 * m_statLocalNoCells / m_statLocalMaxNoCells);
  gridLocal.addBlank();
  gridLocal.addData("Number of elements", m_statLocalNoElements);
  gridLocal.addData("Number of helements", m_statLocalNoHElements);
  gridLocal.addBlank();
  Data& surfaces = gridLocal.addData("Number of surfaces", m_statLocalNoSurfaces);
  surfaces.addData("boundary surfaces", m_statLocalNoBoundarySurfaces);
  surfaces.addData("inner surfaces", m_statLocalNoInnerSurfaces);
  surfaces.addData("MPI surfaces", m_statLocalNoMpiSurfaces);
  gridLocal.addData("Max. number of surfaces (collector size)", m_statLocalMaxNoSurfaces);
  gridLocal.addData("Memory utilization", 100.0 * m_statLocalNoSurfaces / m_statLocalMaxNoSurfaces);

  // GRID SUMMARY (GLOBAL)
  Group& gridGlobal = summary.addGroup("GRID SUMMARY (GLOBAL)");
  gridGlobal.addData("Minimum used polynomial degree", m_statGlobalMinPolyDeg);
  gridGlobal.addData("Maximum used polynomial degree", m_statGlobalMaxPolyDeg);
  gridGlobal.addData("Average used polynomial degree", 1.0 * m_statGlobalPolyDegSum / m_statGlobalNoActiveCells);
  gridGlobal.addBlank();
  gridGlobal.addData("Minimum grid level", m_statGlobalMinLevel);
  gridGlobal.addData("Maximum grid level", m_statGlobalMaxLevel);
  gridGlobal.addBlank();
  Data& noCellsGlbl = gridGlobal.addData("Number of cells", m_statGlobalNoCells);
  noCellsGlbl.addData("internal cells", m_statGlobalNoInternalCells);
  noCellsGlbl.addData("halo cells", m_statGlobalNoHaloCells);
  gridGlobal.addData("Number of active cells", m_statGlobalNoActiveCells);
  gridGlobal.addData("Number of active DOFs", m_statGlobalNoActiveDOFs);
  gridGlobal.addData("Max. number of cells (collector size)", m_statGlobalMaxNoCells);
  gridGlobal.addData("Memory utilization", 100.0 * m_statGlobalNoCells / m_statGlobalMaxNoCells);
  gridGlobal.addBlank();
  gridGlobal.addData("Number of elements", m_statGlobalNoElements);
  gridGlobal.addData("Number of helements", m_statGlobalNoHElements);
  gridGlobal.addBlank();
  Data& surfacesGlobal = gridGlobal.addData("Number of surfaces", m_statGlobalNoSurfaces);
  surfacesGlobal.addData("boundary surfaces", m_statGlobalNoBoundarySurfaces);
  surfacesGlobal.addData("inner surfaces", m_statGlobalNoInnerSurfaces);
  surfacesGlobal.addData("MPI surfaces", m_statGlobalNoMpiSurfaces);
  surfacesGlobal.addData("Surfaces marked for p-ref", m_statGlobalNoPrefSurfs);
  surfacesGlobal.addData("Surfaces marked for h-ref", m_statGlobalNoHrefSurfs);
  gridGlobal.addData("Max. number of surfaces (collector size)", m_statGlobalMaxNoSurfaces);
  gridGlobal.addData("Memory utilization", 100.0 * m_statGlobalNoSurfaces / m_statGlobalMaxNoSurfaces);

  // BOUNDARY CONDITION SUMMARY
  Group& boundary = summary.addGroup("BOUNDARY CONDITIONS");
  for(auto&& bc : m_boundaryConditions) {
    boundary.addData(bc->name(), bc->id());
  }

  // CUT-OFF BOUNDARY SUMMARY
  Group& cutOffBoundary = summary.addGroup("CUT-OFF BOUNDARY CONDITIONS");
  cutOffBoundary.addData("Enabled?", m_useCutOffBoundaries);
  if(m_useCutOffBoundaries) {
    std::array<MString, 6> dirs = {{"-x", "+x", "-y", "+y", "-z", "+z"}};
    for(MInt i = 0; i < 2 * nDim; i++) {
      const MString bc = dirs[i] + " (bcId: " + to_string(m_cutOffBoundaryConditionIds[i]) + ")";
      cutOffBoundary.addData(bc, m_statGlobalNoCutOffBoundarySurfaces[i]);
    }
  }

  // SBP SUMMARY
  Group& sbp = summary.addGroup("SBP MODE");
  sbp.addData("Enabled?", m_sbpMode);
  sbp.addData("Operator", m_sbpOperator);

  std::string s = summary.buildString();

  // Print output to terminal only on root domain, print output to m_log on
  // all domains
  if(isMpiRoot()) {
    cout << s << std::endl;
  }
  m_log << s << endl;
}


/**
 * \brief Saves all available data to disk.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-22
 *
 * Calls saveSolutionFile(const MString& suffix) with the current
 * m_timeStep as the suffix.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveSolutionFile() {
  TRACE();

  stringstream ss;
  ss << setw(8) << setfill('0') << m_timeStep;
  saveSolutionFile(ss.str());
}


/**
 * \brief Saves all available data to disk.
 *
 * \author Michael Schlottke, Sven Berger
 * \date 2013-04-19
 *
 * \param[in] suffix The suffix that should be appended to the generic output
 *                   name (e.g. the current time step).
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveSolutionFile(const MString& suffix) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::IO]);
  RECORD_TIMER_START(m_timers[Timers::SaveSolutionFile]);

  m_log << "Saving solution file ... ";

  // Create ParallelIo instance
  using namespace parallel_io;
  const MString fileName =
      outputDir() + "solution_" + getIdentifier(g_multiSolverGrid, "b") + suffix + ParallelIo::fileExt();
  const MInt noElements = m_elements.size();
  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());

  // If adaptive refinement is active set to correct adapted grid file name
  if(hasAdaptivePref()) {
    parallelIo.setAttribute(suffix + grid().gridInputFileName(), "gridFile");
  } else {
    parallelIo.setAttribute(grid().gridInputFileName(), "gridFile");
  }

  // for g_multiSolverGrids we need the solverId in the file
  parallelIo.setAttribute(solverId(), "solverId");

  // Determine data offset for each element in output buffer
  MIntScratchSpace elementOffset(noElements, AT_, "elementOffset");
  MIntScratchSpace elementNodes(noElements, AT_, "elementNodes");
  for(MInt cellId = 0, offset = 0, elementId = 0; elementId < noElements; cellId++) {
    // If cellId is not the one of the current element (i.e. the current cell
    // does not have an element), add the default offset (based on initPolyDeg)
    if(cellId != m_elements.cellId(elementId)) {
      offset += ipow(m_initNoNodes1D, nDim);
      continue;
    }

    // Otherwise use current offset for the current element
    elementOffset[elementId] = offset;

    // Increase offset based on element polynomial degree
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    offset += noNodesXD;

    // Set number of nodes to use later and proceed with next element
    elementNodes[elementId] = noNodesXD;
    elementId++;
  }

  // Determine local data array size
  // array size = #nodes of elements + (#cells - #elements) * default cell size
  // (cells without elements use polyDeg = initPolyDeg)
  const MInt localNoNodes = m_noTotalNodesXD + (grid().noInternalCells() - noElements) * ipow(m_initNoNodes1D, nDim);

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes;
  ParallelIo::calcOffset(localNoNodes, &nodesOffset, &globalNoNodes, mpiComm());

  // Grid file name and solver type
  parallelIo.setAttribute("DG", "solverType");
  if(m_sbpMode) {
    parallelIo.setAttribute((MInt)m_sbpMode, "sbpMode");
  }
  parallelIo.setAttribute(m_timeStep, "timeStep");
  parallelIo.setAttribute(m_time, "time");

  // Define arrays in file
  // Polyomial degree
  parallelIo.defineArray(PIO_UCHAR, "polyDegs", grid().domainOffset(noDomains()));
  // Number of nodes in SBP mode
  if(m_sbpMode) {
    parallelIo.defineArray(PIO_UCHAR, "noNodes1D", grid().domainOffset(noDomains()));
  }


  // Get information about integration method and polynomial type
  MString dgIntegrationMethod = m_sbpMode ? "DG_INTEGRATE_GAUSS_LOBATTO" : "DG_INTEGRATE_GAUSS";
  dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", m_solverId, AT_, &dgIntegrationMethod);
  MString dgPolynomialType = "DG_POLY_LEGENDRE";
  dgPolynomialType = Context::getSolverProperty<MString>("dgPolynomialType", m_solverId, AT_, &dgPolynomialType);

  // Add integration method & polynomial type to grid file as attributes
  parallelIo.setAttribute(dgIntegrationMethod, "dgIntegrationMethod", "polyDegs");
  parallelIo.setAttribute(dgPolynomialType, "dgPolynomialType", "polyDegs");

  // Set counter to get correct variable names
  MInt varId = 0;

  // Solution
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    const MString name = "variables" + to_string(varId++);
    parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
    parallelIo.setAttribute(m_sysEqn.consVarNames(i), "name", name);
  }

  // Node variables
  if(m_saveNodeVariablesToSolutionFile) {
    for(MInt i = 0; i < SysEqn::noNodeVars(); i++) {
      const MString name = "variables" + to_string(varId++);
      parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
      parallelIo.setAttribute(m_sysEqn.nodeVarNames(i), "name", name);
    }
  }

  // Time derivative
  if(m_writeTimeDerivative) {
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      const MString name = "variables" + to_string(varId++);
      parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
      parallelIo.setAttribute(m_sysEqn.consVarNames(i) + "_tDeriv", "name", name);
    }
  }

  if(useSponge() && m_writeSpongeEta > 0) {
    const MString name = "variables" + to_string(varId);
    parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
    parallelIo.setAttribute("spongeEta", "name", name);
  }

  // Write data to disk
  // Create scratch space with polynomial degree and write it to file
  parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
  MUcharScratchSpace polyDegs(grid().noInternalCells(), FUN_, "polyDegs");
  fill_n(&polyDegs[0], grid().noInternalCells(), m_initPolyDeg);
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    polyDegs.p[cellId] = m_elements.polyDeg(elementId);
  }
  parallelIo.writeArray(polyDegs.begin(), "polyDegs");

  // Create scratch space with noNodes and write it to file
  if(m_sbpMode) {
    parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
    MUcharScratchSpace noNodes1D(grid().noInternalCells(), FUN_, "noNodes1D");
    fill_n(&noNodes1D[0], grid().noInternalCells(), m_initNoNodes1D);
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt cellId = m_elements.cellId(elementId);
      noNodes1D.p[cellId] = m_elements.noNodes1D(elementId);
    }
    parallelIo.writeArray(noNodes1D.begin(), "noNodes1D");
  }

  varId = 0;
  parallelIo.setOffset(localNoNodes, nodesOffset);

  // Solution
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    fill(buffer.begin(), buffer.end(), 0.0);
    for(MInt e = 0; e < noElements; e++) {
      MFloat* const b = &buffer[elementOffset[e]];
      const MFloat* const v = &m_elements.variables(e) + i;
      for(MInt n = 0; n < elementNodes[e]; n++) {
        b[n] = v[n * m_sysEqn.noVars()];
      }
    }
    const MString name = "variables" + to_string(varId++);
    parallelIo.writeArray(&buffer[0], name);
  }

  // Node variables
  if(m_saveNodeVariablesToSolutionFile) {
    for(MInt i = 0; i < SysEqn::noNodeVars(); i++) {
      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
      fill(buffer.begin(), buffer.end(), 0.0);
      for(MInt e = 0; e < noElements; e++) {
        MFloat* const b = &buffer[elementOffset[e]];
        const MFloat* const v = &m_elements.nodeVars(e) + i;
        for(MInt n = 0; n < elementNodes[e]; n++) {
          b[n] = v[n * SysEqn::noNodeVars()];
        }
      }
      const MString name = "variables" + to_string(varId++);
      parallelIo.writeArray(&buffer[0], name);
    }
  }

  // Time derivative
  if(m_writeTimeDerivative) {
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
      fill(buffer.begin(), buffer.end(), 0.0);
      for(MInt e = 0; e < noElements; e++) {
        MFloat* const b = &buffer[elementOffset[e]];
        const MFloat* const r = &m_elements.rightHandSide(e) + i;
        for(MInt n = 0; n < elementNodes[e]; n++) {
          b[n] = r[n * m_sysEqn.noVars()];
        }
      }
      const MString name = "variables" + to_string(varId++);
      parallelIo.writeArray(&buffer[0], name);
    }
  }

  // SpongeEta
  if(useSponge() && m_writeSpongeEta > 0) {
    const MInt noSpongeElements = sponge().noSpongeElements();
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    fill(buffer.begin(), buffer.end(), 0.0);
    for(MInt e = 0; e < noSpongeElements; e++) {
      MInt elementId = sponge().elementId(e);
      MFloat* const b = &buffer[elementOffset[elementId]];
      for(MInt n = 0; n < elementNodes[elementId]; n++) {
        b[n] = sponge().spongeEta(e, n);
      }
    }
    const MString name = "variables" + to_string(varId);
    parallelIo.writeArray(&buffer[0], name);
  }
  m_log << "done" << endl;

  RECORD_TIMER_STOP(m_timers[Timers::SaveSolutionFile]);
  RECORD_TIMER_STOP(m_timers[Timers::IO]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
}


/**
 * \brief Saves a file to disk with all information that is necessary to restart
 *        the calculations from here.
 *
 * \author Michael Schlottke, Sven Berger
 * \date 2012-12-22
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveRestartFile() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::IO]);
  RECORD_TIMER_START(m_timers[Timers::SaveRestartFile]);

  // Determine file name and grid file name
  const MString fileName = getRestartFileName(m_timeStep, m_useNonSpecifiedRestartFile);

  stringstream ss;
  ss << "restart_" << getIdentifier(g_multiSolverGrid, "b") << setw(8) << setfill('0') << m_timeStep;
  const MString suffix = ss.str();

  // Determine data offset for each element in output buffer
  const MInt noElements = m_elements.size();
  MIntScratchSpace elementOffset(noElements, AT_, "elementOffset");
  MIntScratchSpace elementNodes(noElements, AT_, "elementNodes");
  for(MInt cellId = 0, offset = 0, elementId = 0; elementId < noElements; cellId++) {
    // If cellId is not the one of the current element (i.e. the current cell
    // does not have an element), add the default offset (based on initPolyDeg)
    if(cellId != m_elements.cellId(elementId)) {
      offset += ipow(m_initNoNodes1D, nDim);
      continue;
    }

    // Otherwise use current offset for the current element
    elementOffset[elementId] = offset;

    // Increase offset based on element polynomial degree
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    offset += noNodesXD;

    // Set number of nodes to use later and proceed with next element
    elementNodes[elementId] = noNodesXD;
    elementId++;
  }

  // Determine local data array size
  // array size = #nodes of elements + (#cells - #elements) * default cell size
  // (cells without elements use polyDeg = initPolyDeg)
  const MInt localNoNodes = m_noTotalNodesXD + (grid().noInternalCells() - noElements) * ipow(m_initNoNodes1D, nDim);

  // Create data out object to save data to disk
  using namespace parallel_io;
  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());

  // Set attributes
  // If adaptive refinement is active set to correct adapted grid file name
  if(hasAdaptivePref()) {
    parallelIo.setAttribute(suffix + grid().gridInputFileName(), "gridFile");
  } else {
    parallelIo.setAttribute(grid().gridInputFileName(), "gridFile");
  }
  parallelIo.setAttribute("DG", "solverType");
  if(m_sbpMode) {
    parallelIo.setAttribute((MInt)m_sbpMode, "sbpMode");
  }
  parallelIo.setAttribute(solverId(), "solverId");

  // Set time variables
  parallelIo.defineScalar(PIO_INT, "timeStep");
  parallelIo.defineScalar(PIO_FLOAT, "time");
  parallelIo.defineScalar(PIO_FLOAT, "dt");

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes;
  ParallelIo::calcOffset(localNoNodes, &nodesOffset, &globalNoNodes, mpiComm());

  // Define arrays in file
  // Polynomial degree
  parallelIo.defineArray(PIO_UCHAR, "polyDegs", grid().domainOffset(noDomains()));
  // Number of nodes in SBP mode
  if(m_sbpMode) {
    parallelIo.defineArray(PIO_UCHAR, "noNodes1D", grid().domainOffset(noDomains()));
  }

  // Get information about integration method and polynomial type
  MString dgIntegrationMethod = m_sbpMode ? "DG_INTEGRATE_GAUSS_LOBATTO" : "DG_INTEGRATE_GAUSS";
  dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", m_solverId, AT_, &dgIntegrationMethod);
  MString dgPolynomialType = "DG_POLY_LEGENDRE";
  dgPolynomialType = Context::getSolverProperty<MString>("dgPolynomialType", m_solverId, AT_, &dgPolynomialType);
  // Add integration method & polynomial type to restart file as attributes
  parallelIo.setAttribute(dgIntegrationMethod, "dgIntegrationMethod", "polyDegs");
  parallelIo.setAttribute(dgPolynomialType, "dgPolynomialType", "polyDegs");

  // Set counter to get correct variable names
  MInt varId = 0;

  // Solution
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    const MString name = "variables" + to_string(varId++);
    parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
    parallelIo.setAttribute(m_sysEqn.consVarNames(i), "name", name);
  }

  // Node variables
  if(SysEqn::hasTimeDependentNodeVars()) {
    for(MInt i = 0; i < SysEqn::noNodeVars(); i++) {
      const MString name = "variables" + to_string(varId++);
      parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
      parallelIo.setAttribute(m_sysEqn.nodeVarNames(i), "name", name);
    }
  }

  // Note: kept here to have a template for writing a data array that is only defined for all
  // elements but should be filled and written on a cell basis for visualization
  // Coupling variables
  // if (hasCoupling()) {
  //  for (auto& cc : m_couplingConditions) {
  //    for (MInt i = 0; i < cc->noRestartVars(); i++) {
  //      const MString name = "variables" + to_string(varId++);
  //      parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
  //      parallelIo.setAttribute(cc->restartVarName(i), "name", name);
  //      parallelIo.setAttribute("couplingCondition", "type", name);
  //      parallelIo.setAttribute(cc->name(), "ccName", name);
  //      parallelIo.setAttribute(cc->id(), "ccId", name);
  //    }
  //  }
  //}

  // Boundary conditions
  MInt bcId = 0;
  const MInt noBcIds = m_boundaryConditions.size();
  std::vector<ParallelIo::size_type> bcNodesOffset(noBcIds, 0), bcLocalNoNodes(noBcIds, 0), bcGlobalNoNodes(noBcIds, 0);

  for(auto&& bc : m_boundaryConditions) {
    const MInt bcNoVars = bc->noRestartVars();
    if(bcNoVars > 0) {
      bcLocalNoNodes[bcId] = bc->getLocalNoNodes();
      ParallelIo::calcOffset(bcLocalNoNodes[bcId], &bcNodesOffset[bcId], &bcGlobalNoNodes[bcId], mpiComm());
      for(MInt i = 0; i < bcNoVars; i++) {
        const MString name = "variables" + to_string(varId++);
        parallelIo.defineArray(PIO_FLOAT, name, bcGlobalNoNodes[bcId]);
        parallelIo.setAttribute(bc->restartVarName(i), "name", name);
      }
    }
    bcId++;
  }

  // Write data to disk

  // Write time variables to file
  parallelIo.writeScalar(m_timeStep, "timeStep");
  parallelIo.writeScalar(m_time, "time");
  parallelIo.writeScalar(m_dt, "dt");

  // Create scratch space with polynomial degree and write it to file
  parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
  MUcharScratchSpace polyDegs(grid().noInternalCells(), FUN_, "polyDegs");
  fill_n(&polyDegs[0], grid().noInternalCells(), m_initPolyDeg);
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    polyDegs.p[cellId] = m_elements.polyDeg(elementId);
  }
  parallelIo.writeArray(polyDegs.begin(), "polyDegs");

  // Create scratch space with noNodes and write it to file
  if(m_sbpMode) {
    parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
    MUcharScratchSpace noNodes1D(grid().noInternalCells(), FUN_, "noNodes1D");
    fill_n(&noNodes1D[0], grid().noInternalCells(), m_initNoNodes1D);
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt cellId = m_elements.cellId(elementId);
      noNodes1D.p[cellId] = m_elements.noNodes1D(elementId);
    }
    parallelIo.writeArray(noNodes1D.begin(), "noNodes1D");
  }

  varId = 0;
  parallelIo.setOffset(localNoNodes, nodesOffset);

  // Solution
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    // Solution
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    fill(buffer.begin(), buffer.end(), 0.0);
    for(MInt e = 0; e < noElements; e++) {
      MFloat* const b = &buffer[elementOffset[e]];
      const MFloat* const v = &m_elements.variables(e) + i;
      for(MInt n = 0; n < elementNodes[e]; n++) {
        b[n] = v[n * m_sysEqn.noVars()];
      }
    }
    const MString name = "variables" + to_string(varId++);
    parallelIo.writeArray(&buffer[0], name);
  }

  // Node variables
  if(SysEqn::hasTimeDependentNodeVars()) {
    for(MInt i = 0; i < SysEqn::noNodeVars(); i++) {
      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
      fill(buffer.begin(), buffer.end(), 0.0);
      for(MInt e = 0; e < noElements; e++) {
        MFloat* const b = &buffer[elementOffset[e]];
        const MFloat* const v = &m_elements.nodeVars(e) + i;
        for(MInt n = 0; n < elementNodes[e]; n++) {
          b[n] = v[n * SysEqn::noNodeVars()];
        }
      }
      const MString name = "variables" + to_string(varId++);
      parallelIo.writeArray(&buffer[0], name);
    }
  }

  // Note: keep this, see comment above at definition of coupling variables
  // Coupling variables
  // if (hasCoupling()) {
  //  for (auto& cc : m_couplingConditions) {
  //    for (MInt i = 0; i < cc->noRestartVars(); i++) {
  //      // Load restart variable from coupling class
  //      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
  //      cc->getRestartVariable(i, &buffer[0]);
  //
  //      // Re-arrange data in buffer to add gaps for each cell without element
  //      MInt offset = m_noTotalNodesXD;
  //      for (MInt e = noElements - 1; e >= 0; e--) {
  //        // Update offset
  //        const MInt noNodesXD = m_elements.noNodesXD(e);
  //        offset -= noNodesXD;
  //
  //        // Create pointers to data
  //        MFloat* const source = &buffer[offset];
  //        MFloat* const destination = &buffer[elementOffset[e]];
  //
  //        // No copy operation needed if pointers are equal
  //        if (destination == source) {
  //          continue;
  //        }
  //
  //        // Copy (move) data
  //        copy_backward(source, source + noNodesXD, destination + noNodesXD);
  //
  //        // Fill up original storage location with zeros
  //        // The min(...) expression makes sure that
  //        // - if the original and the final storage location overlap, only the
  //        //   values that are not in the final storage location are reset
  //        // - if the original and the final storage location DO NOT overlap,
  //        // - the entire original region is reset with zeros
  //        fill_n(source, min(static_cast<MLong>(noNodesXD), destination - source), 0.0);
  //      }
  //
  //      // Write data from buffer to file
  //      const MString name = "variables" + to_string(varId++);
  //      parallelIo.writeArray(&buffer[0], name);
  //    }
  //  }
  //}

  // Boundary conditions
  bcId = 0;
  for(auto&& bc : m_boundaryConditions) {
    parallelIo.setOffset(bcLocalNoNodes[bcId], bcNodesOffset[bcId]);
    for(MInt i = 0; i < bc->noRestartVars(); i++) {
      // Avoid dereferencing a zero length array
      MFloatScratchSpace buffer(max(bcLocalNoNodes[bcId], 1l), AT_, "buffer");
      fill(buffer.begin(), buffer.end(), 0.0);

      // Get boundary condition restart variable
      bc->getRestartVariable(i, &buffer[0]);
      const MString name = "variables" + to_string(varId++);
      parallelIo.writeArray(&buffer[0], name);
    }
    bcId++;
  }

  RECORD_TIMER_STOP(m_timers[Timers::SaveRestartFile]);
  RECORD_TIMER_STOP(m_timers[Timers::IO]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
}


/**
 * \brief Load restart file with all information that is necessary to restart
 *        the calculations from here.
 *
 * \author Michael Schlottke, Sven Berger
 * \date 2013-09-23
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::loadRestartFile() {
  TRACE();

  m_log << "Loading restart file... ";

  // Determine file name and grid file name
  const MString fileName = getRestartFileName(m_restartTimeStep, m_useNonSpecifiedRestartFile);

  // Create big data object to load data from disk
  using namespace parallel_io;
  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  // Get attributes
  MString gridFileName;
  parallelIo.getAttribute(&gridFileName, "gridFile");

  // Get time variables
  parallelIo.readScalar(&m_timeStep, "timeStep");
  parallelIo.readScalar(&m_time, "time");
  parallelIo.readScalar(&m_dt, "dt");

  // Load polynomial degrees
  MIntScratchSpace polyDegs(grid().noInternalCells(), AT_, "buffer");
  if(!parallelIo.hasDataset("polyDegs", 1)) {
    TERMM(1, "ERROR: restart file is not valid (missing polynomial degree data).");
  }
  parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
  parallelIo.readArray(&polyDegs[0], "polyDegs");

  // update elements
  const MInt noElements = m_elements.size();
  m_noTotalNodesXD = 0;
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    const MInt polyDeg = polyDegs[cellId];
    const MInt noNodes1D = polyDeg + 1;
    const MInt noNodesXD = ipow(polyDeg + 1, nDim);
    m_elements.polyDeg(elementId) = polyDeg;
    m_elements.noNodes1D(elementId) = noNodes1D;
    m_noTotalNodesXD += noNodesXD;
  }

  // Determine data offset for each element in output buffer
  MIntScratchSpace elementOffset(noElements, AT_, "elementOffset");
  for(MInt cellId = 0, offset = 0, elementId = 0; elementId < noElements; cellId++) {
    // If cellId is not the one of the current element (i.e. the current cell
    // does not have an element), add the default offset (based on initPolyDeg)
    if(cellId != m_elements.cellId(elementId)) {
      offset += ipow(m_initNoNodes1D, nDim);
      continue;
    }

    // Otherwise use current offset for the current element
    elementOffset[elementId] = offset;

    // Increase offset based on element polynomial degree
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    offset += noNodesXD;

    // Set number of nodes to use later and proceed with next element
    elementId++;
  }

  // Set data offsets and array sizes
  // Determine local data array size
  // array size = #nodes of elements + (#cells - #elements) * default cell size
  // (cells without elements use polyDeg = initPolyDeg)
  const MInt localNoNodes = m_noTotalNodesXD + (grid().noInternalCells() - noElements) * ipow(m_initNoNodes1D, nDim);

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes;
  ParallelIo::calcOffset(localNoNodes, &nodesOffset, &globalNoNodes, mpiComm());

  // Load data from disk
  parallelIo.setOffset(localNoNodes, nodesOffset);
  MInt varId = 0;

  // Solution
  for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    const MString name = "variables" + to_string(varId++);
    parallelIo.readArray(&buffer[0], name);
    for(MInt e = 0; e < noElements; e++) {
      const MFloat* const b = &buffer[elementOffset[e]];
      MFloat* const v = &m_elements.variables(e) + i;
      const MInt noNodesXD = m_elements.noNodesXD(e);
      for(MInt n = 0; n < noNodesXD; n++) {
        v[n * m_sysEqn.noVars()] = b[n];
      }
    }
  }

  // Node variables
  if(SysEqn::hasTimeDependentNodeVars()) {
    for(MInt i = 0; i < m_sysEqn.noNodeVars(); i++) {
      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
      const MString name = "variables" + to_string(varId++);
      parallelIo.readArray(&buffer[0], name);
      for(MInt e = 0; e < noElements; e++) {
        const MFloat* const b = &buffer[elementOffset[e]];
        MFloat* const v = &m_elements.nodeVars(e) + i;
        const MInt noNodesXD = m_elements.noNodesXD(e);
        for(MInt n = 0; n < noNodesXD; n++) {
          v[n * m_sysEqn.noNodeVars()] = b[n];
        }
      }
    }
  }

  // Note: keep this, might be useful as a reference sometime
  // Coupling variables
  // if (hasCoupling()) {
  //  for (auto& cc : m_couplingConditions) {
  //    for (MInt i = 0; i < cc->noRestartVars(); i++) {
  //      // Load data into local buffer
  //      MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
  //      const MString name = "variables" + to_string(varId++);
  //      parallelIo.readArray(&buffer[0], name);
  //
  //      // Re-arrange data in buffer to make it "gap-less", i.e. for each
  //      // element the data is adjacent to the neighboring elements
  //      MInt offset = 0;
  //      for (MInt e = 0; e < noElements; e++) {
  //        // Create pointers to data
  //        const MFloat* const b = &buffer[elementOffset[e]];
  //        MFloat* const v = &buffer[offset];
  //
  //        // Update offset
  //        const MInt noNodesXD = m_elements.noNodesXD(e);
  //        offset += noNodesXD;
  //
  //        // No copy operation needed if pointers are equal
  //        if (b == v) {
  //          continue;
  //        }
  //
  //        // Copy (move) data
  //        copy(b, b + noNodesXD, v);
  //      }
  //
  //      // Store restart variable in coupling class
  //      cc->setRestartVariable(i, &buffer[0]);
  //    }
  //  }
  //}

  // Boundary conditions
  for(auto&& bc : m_boundaryConditions) {
    const MInt bcNoVars = bc->noRestartVars();
    if(bcNoVars > 0) {
      const MInt bcLocalNoNodes = bc->getLocalNoNodes();
      ParallelIo::size_type bcNodesOffset, bcGlobalNoNodes;
      ParallelIo::calcOffset(bcLocalNoNodes, &bcNodesOffset, &bcGlobalNoNodes, mpiComm());
      parallelIo.setOffset(bcLocalNoNodes, bcNodesOffset);
      for(MInt i = 0; i < bcNoVars; i++) {
        // Avoid dereferencing a zero length array
        MFloatScratchSpace buffer(max(bcLocalNoNodes, 1), AT_, "buffer");
        const MString name = "variables" + to_string(varId++);
        parallelIo.readArray(&buffer[0], name);

        // Store restart variable in boundary condition
        bc->setRestartVariable(i, &buffer[0]);
      }
    }
  }

  m_log << "done" << endl;
}


/// \brief Save nodal data to file.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-26
///
/// \param[in] fileNameBase Base name of file (i.e., without directory or file
///                         extension) to write noda data to.
/// \param[in] noVars Number of variables to save to file.
/// \param[in] varNames List of variable names.
/// \param[in] data Nodal data to save.
///
/// This function allows to save nodal data to a file. The data needs to be
/// stored with a data solver size per element of the (maximum number of nodes) x
/// (number of variables), i.e, it has to have the same structure as the
/// m_variables field in m_elements. If noVars > 1, it is expected that the
/// variables are stored in an array-of-structs layout, i.e., var0 var1 var2
/// var0 var1 var2 var0 ...
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveNodalData(const MString& fileNameBase,
                                                    const MInt noVars,
                                                    const vector<MString>& varNames,
                                                    const MFloat* const data) const {
  TRACE();
  CHECK_TIMERS_IO();

  // Sanity check
  if(noVars < 1) {
    TERMM(1, "noVars must by >= 1");
  }

  // Determine data offset for each element in output buffer
  const MInt noElements = m_elements.size();
  MIntScratchSpace elementOffset(noElements, AT_, "elementOffset");
  MIntScratchSpace elementNodes(noElements, AT_, "elementNodes");
  for(MInt cellId = 0, offset = 0, elementId = 0; elementId < noElements; cellId++) {
    // If cellId is not the one of the current element (i.e. the current cell
    // does not have an element), add the default offset (based on initPolyDeg)
    if(cellId != m_elements.cellId(elementId)) {
      offset += ipow(m_initNoNodes1D, nDim);
      continue;
    }

    // Otherwise use current offset for the current element
    elementOffset[elementId] = offset;

    // Increase offset based on element polynomial degree
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    offset += noNodesXD;

    // Set number of nodes to use later and proceed with next element
    elementNodes[elementId] = noNodesXD;
    elementId++;
  }

  // Determine local data array size
  // array size = #nodes of elements + (#cells - #elements) * default cell size
  // (cells without elements use polyDeg = initPolyDeg)
  const MInt localNoNodes = m_noTotalNodesXD + (grid().noInternalCells() - noElements) * ipow(m_initNoNodes1D, nDim);

  // Create data out object to save data to disk
  const MString fileName = outputDir() + fileNameBase + ParallelIo::fileExt();
  using namespace parallel_io;
  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());

  // Set attributes
  parallelIo.setAttribute(grid().gridInputFileName(), "gridFile");
  parallelIo.setAttribute("DG", "solverType");
  if(m_sbpMode) {
    parallelIo.setAttribute((MInt)m_sbpMode, "sbpMode");
  }
  // @ansgar TODO labels:DG,toenhance change this in the future to always write the solver id, this will change a lot
  // of reference files, e.g., the nodevars
  if(grid().raw().treeb().noSolvers() > 1) {
    parallelIo.setAttribute(solverId(), "solverId");
  }

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes;
  ParallelIo::calcOffset(localNoNodes, &nodesOffset, &globalNoNodes, mpiComm());

  // Define arrays in file
  // Polyomial degree
  parallelIo.defineArray(PIO_UCHAR, "polyDegs", grid().domainOffset(noDomains()));
  // Number of nodes in SBP mode
  if(m_sbpMode) {
    parallelIo.defineArray(PIO_UCHAR, "noNodes1D", grid().domainOffset(noDomains()));
  }

  // Get information about integration method and polynomial type
  MString dgIntegrationMethod = m_sbpMode ? "DG_INTEGRATE_GAUSS_LOBATTO" : "DG_INTEGRATE_GAUS";
  dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", m_solverId, AT_, &dgIntegrationMethod);
  MString dgPolynomialType = "DG_POLY_LEGENDRE";
  dgPolynomialType = Context::getSolverProperty<MString>("dgPolynomialType", m_solverId, AT_, &dgPolynomialType);

  // Add integration method & polynomial type to grid file as attributes
  parallelIo.setAttribute(dgIntegrationMethod, "dgIntegrationMethod", "polyDegs");
  parallelIo.setAttribute(dgPolynomialType, "dgPolynomialType", "polyDegs");

  // Add arrays to file
  for(MInt i = 0; i < noVars; i++) {
    const MString name = "variables" + to_string(i);
    parallelIo.defineArray(PIO_FLOAT, name, globalNoNodes);
    parallelIo.setAttribute(varNames[i], "name", name);
  }

  const MInt maxNoNodesXD = ipow(m_maxNoNodes1D, nDim);
  const MInt elementDataSize = noVars * maxNoNodesXD;

  // Write data to disk
  // Create scratch space with polynomial degree and write it to file
  parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
  MUcharScratchSpace polyDegs(grid().noInternalCells(), FUN_, "polyDegs");
  fill_n(&polyDegs[0], grid().noInternalCells(), m_initPolyDeg);
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    polyDegs.p[cellId] = m_elements.polyDeg(elementId);
  }
  parallelIo.writeArray(polyDegs.begin(), "polyDegs");

  // Create scratch space with noNodes and write it to file
  if(m_sbpMode) {
    parallelIo.setOffset(grid().noInternalCells(), grid().domainOffset(domainId()));
    MUcharScratchSpace noNodes1D(grid().noInternalCells(), FUN_, "noNodes1D");
    fill_n(&noNodes1D[0], grid().noInternalCells(), m_initNoNodes1D);
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt cellId = m_elements.cellId(elementId);
      noNodes1D.p[cellId] = m_elements.noNodes1D(elementId);
    }
    parallelIo.writeArray(noNodes1D.begin(), "noNodes1D");
  }

  parallelIo.setOffset(localNoNodes, nodesOffset);
  for(MInt i = 0; i < noVars; i++) {
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    fill(buffer.begin(), buffer.end(), 0.0);
    for(MInt e = 0; e < noElements; e++) {
      MFloat* const b = &buffer[elementOffset[e]];
      const MInt dataOffset = e * elementDataSize + i;
      const MFloat* const v = &data[dataOffset];
      for(MInt n = 0; n < elementNodes[e]; n++) {
        b[n] = v[n * noVars];
      }
    }
    const MString name = "variables" + to_string(i);
    parallelIo.writeArray(&buffer[0], name);
  }
}


/// \brief Save node variables to file.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-26
///
/// This is mainly used as a restart file if the node variables are constant in
/// time.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::saveNodeVarsFile() {
  TRACE();
  CHECK_TIMERS_IO();

  m_log << "Saving file with node variable data for restarting... ";

  // Add solver number in case of multisolver execution
  stringstream ss;
  ss << "nodevars" << getIdentifier(g_multiSolverGrid, "_b", "");

  // Assemble node variable names
  vector<MString> nodeVarNames(SysEqn::noNodeVars());
  for(MInt nodeVar = 0; nodeVar < SysEqn::noNodeVars(); nodeVar++) {
    nodeVarNames[nodeVar] = m_sysEqn.nodeVarNames(nodeVar);
  }

  // Save node variables
  saveNodalData(ss.str(), SysEqn::noNodeVars(), nodeVarNames, &m_elements.nodeVars(0));

  m_log << "done" << endl;
}


/// \brief Load node variables from file.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-26
///
/// This is mainly used to restart from a previously stored file if the node
/// variables are constant in time.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::loadNodeVarsFile() {
  TRACE();
  CHECK_TIMERS_IO();

  m_log << "Loading nodevars file... ";

  // Add solver number in case of multisolver execution
  stringstream ss;
  ss << restartDir() << "nodevars" << getIdentifier(g_multiSolverGrid, "_b", "") << ParallelIo::fileExt();

  // Create big data object to load data from disk
  using namespace parallel_io;
  ParallelIo parallelIo(ss.str(), PIO_READ, mpiComm());

  // Determine data offset for each element in output buffer
  const MInt noElements = m_elements.size();
  MIntScratchSpace elementOffset(noElements, AT_, "elementOffset");
  for(MInt cellId = 0, offset = 0, elementId = 0; elementId < noElements; cellId++) {
    // If cellId is not the one of the current element (i.e. the current cell
    // does not have an element), add the default offset (based on initPolyDeg)
    if(cellId != m_elements.cellId(elementId)) {
      offset += ipow(m_initNoNodes1D, nDim);
      continue;
    }

    // Otherwise use current offset for the current element
    elementOffset[elementId] = offset;

    // Increase offset based on element polynomial degree
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    offset += noNodesXD;

    // Set number of nodes to use later and proceed with next element
    elementId++;
  }

  // Set data offsets and array sizes
  // Determine local data array size
  // array size = #nodes of elements + (#cells - #elements) * default cell size
  // (cells without elements use polyDeg = initPolyDeg)
  const MInt localNoNodes = m_noTotalNodesXD + (grid().noInternalCells() - noElements) * ipow(m_initNoNodes1D, nDim);

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes;
  ParallelIo::calcOffset(localNoNodes, &nodesOffset, &globalNoNodes, mpiComm());

  // Load data from disk
  parallelIo.setOffset(localNoNodes, nodesOffset);
  for(MInt i = 0; i < m_sysEqn.noNodeVars(); i++) {
    MFloatScratchSpace buffer(localNoNodes, AT_, "buffer");
    const MString name = "variables" + to_string(i);
    parallelIo.readArray(&buffer[0], name);
    for(MInt e = 0; e < noElements; e++) {
      const MFloat* const b = &buffer[elementOffset[e]];
      MFloat* const v = &m_elements.nodeVars(e) + i;
      const MInt noNodesXD = m_elements.noNodesXD(e);
      for(MInt n = 0; n < noNodesXD; n++) {
        v[n * m_sysEqn.noNodeVars()] = b[n];
      }
    }
  }

  m_log << "done" << endl;
}


/**
 * \brief Return name of a restart file for the given time step.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-09-22
 *
 * \param[in] timeStep Time step to calculate the name for.
 * \param[in] useNonSpecifiedRestartFile If true, return restart file name
 *                                       without time step code.
 * \return Name of the restart file.
 */
template <MInt nDim, class SysEqn>
MString DgCartesianSolver<nDim, SysEqn>::getRestartFileName(const MInt timeStep,
                                                            const MInt useNonSpecifiedRestartFile) {
  TRACE();

  stringstream ss;
  if(useNonSpecifiedRestartFile) {
    ss << restartDir() << "restart_" << getIdentifier(g_multiSolverGrid, "b", "") << ParallelIo::fileExt();
  } else {
    ss << restartDir() << "restart_" << getIdentifier(g_multiSolverGrid, "b") << setw(8) << setfill('0') << timeStep
       << ParallelIo::fileExt();
  }
  const MString fileName = ss.str();

  return fileName;
}


template <MInt nDim, class SysEqn>
typename DgCartesianSolver<nDim, SysEqn>::BC DgCartesianSolver<nDim, SysEqn>::make_bc(MInt bcId) {
  TRACE();

  BC bc = m_boundaryConditionFactory.create(bcId);

  return bc;
}

/// \brief Main routine to calculate the discontinuous Galerkin time derivative.
///        After this method was called, m_rightHandSide contains the time
///        derivative of the conservative variables at the integration points.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2013-04-18
///
/// \param[in] t Physical time at which the time derivative is calculated.
/// \param[in] tStage Pseudo time of the Runge-Kutta stage (needed for
///                   time-dependent source terms).
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcDgTimeDerivative(const MFloat NotUsed(t), const MFloat tStage) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::TimeDeriv]);

  // Reset the time derivative to zero
  resetRHS();

  // Calculate solution on all surfaces and start exchange of surface data
  // In split MPI mode most of the time this was already done at the end of
  // timeStepRk(). If there are no open MPI requests it needs to be done here.
  if(!g_splitMpiComm || !m_mpiRecvRequestsOpen) {
    // Extrapolate the solution vector U from the element volume to all surfaces
    prolongToSurfaces();

    // Apply forward mortar projection to h- and/or p-refined surfaces
    applyForwardProjection();

    // Start exchanging surface data
    startMpiSurfaceExchange();
  }

  // Calculate the volume integral and update dU/dt
  calcVolumeIntegral();

  // Calculate the fluxes on the internal surfaces
  calcInnerSurfaceFlux();

  // Exclude communication from domain load and add to idle time
  // Finish exchanging surface data
  finishMpiSurfaceExchange();

  // Calculate the fluxes on the boundary surfaces
  // Note: this was previously called after calcVolumeIntegral(). It was moved
  // in order to make the radiation boundary conditions work in parallel without
  // any additional MPI communication.
  calcBoundarySurfaceFlux(tStage);

  // Calculate the fluxes on the MPI surfaces
  calcMpiSurfaceFlux();

  // Calculate the integrals for all surfaces and update dU/dt
  calcSurfaceIntegral();

  // Apply Jacobian to dU/dt
  applyJacobian();

  // Calculate source terms and add them to dU/dt
  calcSourceTerms(tStage);

  // Add external (coupling) source terms
  applyExternalSourceTerms(tStage);

  // Calculate sponge terms and add them to dU/dt
  if(useSponge()) {
    RECORD_TIMER_START(m_timers[Timers::Sponge]);
    sponge().calcSourceTerms();
    RECORD_TIMER_STOP(m_timers[Timers::Sponge]);
  }

  RECORD_TIMER_STOP(m_timers[Timers::TimeDeriv]);
}

/// \brief Extrapolate the solution from inside the elements to the surfaces.
///
/// \author Sven Berger
/// \date   November 2014
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::prolongToSurfaces() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Prolong]);

  prolongToSurfaces(0, m_elements.size());

  RECORD_TIMER_STOP(m_timers[Timers::Prolong]);
}


/// \brief Extrapolate the solution in the given range of elements
///        from the elements to the surfaces
///
/// \author Sven Berger
/// \date   November 2014
///
/// \param[in] begin Index of the first element to consider.
/// \param[in] end Index+1 of the last element to consider.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::prolongToSurfaces(const MInt begin, const MInt end) {
  TRACE();

  const MInt* polyDegs = &m_elements.polyDeg(0);
  const MInt* noNodes = &m_elements.noNodes1D(0);
  const MInt* surfaceIds = &m_elements.surfaceIds(0, 0);

  using namespace dg::interpolation;

  // IMPORTANT:
  // If anything in this method is changed, please check if
  // updateNodeVariables() needs to be changed as well, since it contains more
  // or less the same code.

  if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Loop over all elements in given range
    for(MInt elementId = begin; elementId < end; elementId++) {
      const MInt surfaceIdOffset = elementId * 2 * nDim;
      const MInt polyDeg = polyDegs[elementId];
      const MInt noNodes1D = noNodes[elementId];
      const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];

      // Extrapolate the solution to each surface on the faces
      for(MInt dir = 0; dir < 2 * nDim; dir++) {
        const MInt srfcId = surfaceIds[surfaceIdOffset + dir];
        const MInt side = 1 - dir % 2;

        MFloat* src = &m_elements.variables(elementId);
        MFloat* dest = &m_surfaces.variables(srfcId, side);

        prolongToFaceGauss<nDim, SysEqn::noVars()>(src, dir, noNodes1D, &interp.m_LFace[0][0], &interp.m_LFace[1][0],
                                                   dest);
      }
    }
  } else if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Loop over all elements in given range
    for(MInt elementId = begin; elementId < end; elementId++) {
      const MInt surfaceIdOffset = elementId * 2 * nDim;
      const MInt noNodes1D = noNodes[elementId];

      // Extrapolate the solution to each surface on the faces
      for(MInt dir = 0; dir < 2 * nDim; dir++) {
        const MInt srfcId = surfaceIds[surfaceIdOffset + dir];
        const MInt side = 1 - dir % 2;

        MFloat* src = &m_elements.variables(elementId);
        MFloat* dest = &m_surfaces.variables(srfcId, side);

        prolongToFaceGaussLobatto<nDim, SysEqn::noVars()>(src, dir, noNodes1D, dest);
      }
    }
  }
}


/// Apply forward mortar projections (hp-refinement).
///
/// \author Sven Berger, Michael Schlottke
/// \date 2015-03-02
///
/// This method does two things: first, all missing results from the prolong
/// step are copied to the remaining surfaces for h-refined faces. Then, both h-
/// and p-refinement mortar projections are first applied to h-refined surfaces
/// and then p-refinement mortar projections are applied to purely p-refined
/// surfaces.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::applyForwardProjection() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::ForwardProjection]);

  // IMPORTANT:
  // If anything in this method is changed, please check if
  // updateNodeVariables() needs to be changed as well, since it contains more
  // or less the same code.

  const MInt* surfaceIds = &m_elements.surfaceIds(0, 0);
  const MInt noVars = SysEqn::noVars();
  const MInt maxNoNodes1D = m_maxNoNodes1D;
  const MInt maxNoNodes1D3 = (nDim == 3) ? maxNoNodes1D : 1;
  const MInt noHElements = m_helements.size();
  const MInt noDirs = 2 * nDim;
  const MInt noSurfs = 2 * (nDim - 1);

  //////////////////////////////////////////////////////////////////////////////
  // Copy prolong results to h-refined surfaces (the prolong step stored its
  // results only in the first refined surface)
  //////////////////////////////////////////////////////////////////////////////
  if(noHElements > 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt hElementId = 0; hElementId < noHElements; hElementId++) {
      const MInt coarseElementId = m_helements.elementId(hElementId);
      for(MInt dir = 0; dir < noDirs; dir++) {
        const MInt coarseSrfcId = m_elements.surfaceIds(coarseElementId, dir);
        for(MInt pos = 0; pos < noSurfs; pos++) {
          const MInt hSrfcId = m_helements.hrefSurfaceIds(hElementId, dir, pos);
          const MInt side = 1 - dir % 2;

          // Skip if this surface does not exist
          if(hSrfcId == -1) {
            continue;
          }

          // Copy values of prolong step to all remaining refined surfaces
          if(coarseSrfcId != hSrfcId) {
            const MInt noNodes1D = m_elements.noNodes1D(coarseElementId);
            const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
            const MInt size = noNodes1D * noNodes1D3 * SysEqn::noVars();

            // Copy values from prolong step to surface
            copy_n(&m_surfaces.variables(coarseSrfcId, side), size, &m_surfaces.variables(hSrfcId, side));
          }
        }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // Apply mortar projection
  //////////////////////////////////////////////////////////////////////////////
  const MInt noElements = m_elements.size();
  MFloatTensor projected(maxNoNodes1D, maxNoNodes1D3, noVars);

  // Apply projection to h- (and possibly p-)refined surfaces
  if(noHElements > 0) {
#ifdef _OPENMP
#pragma omp parallel for firstprivate(projected)
#endif
    for(MInt hElementId = 0; hElementId < noHElements; hElementId++) {
      for(MInt dir = 0; dir < noDirs; dir++) {
        for(MInt pos = 0; pos < noSurfs; pos++) {
          const MInt hSrfcId = m_helements.hrefSurfaceIds(hElementId, dir, pos);
          const MInt side = 1 - dir % 2;

          // Skip if this surface does not exist
          if(hSrfcId == -1) {
            continue;
          }

          const MInt surfaceNoNodes1D = m_surfaces.noNodes1D(hSrfcId);
          const MInt surfaceNoNodes1D3 = (nDim == 3) ? surfaceNoNodes1D : 1;
          const MInt size = surfaceNoNodes1D * surfaceNoNodes1D3 * noVars;

          // Project the coarse side to the surface
          calcMortarProjection<dg::mortar::forward, SysEqn::noVars()>(
              hSrfcId, dir, &m_surfaces.variables(hSrfcId, side), &projected[0], m_elements, m_surfaces);

          // Copy results of projection back to surface
          copy_n(&projected[0], size, &m_surfaces.variables(hSrfcId, side));
        }
      }
    }
  }

  // Apply projection to pure p-refined surfaces
#ifdef _OPENMP
#pragma omp parallel for firstprivate(projected)
#endif
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt surfaceIdOffset = elementId * 2 * nDim;

    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      // Define auxiliary variables for better readability
      const MInt srfcId = surfaceIds[surfaceIdOffset + dir];

      // Skip boundary surfaces as they are *always* conforming
      if(srfcId < m_innerSurfacesOffset) {
        continue;
      }

      const MInt surfacePolyDeg = m_surfaces.polyDeg(srfcId);
      const MInt elementPolyDeg = m_elements.polyDeg(elementId);
      const MInt side = 1 - dir % 2;
      const MInt surfaceNoNodes1D = m_surfaces.noNodes1D(srfcId);
      const MInt surfaceNoNodes1D3 = (nDim == 3) ? surfaceNoNodes1D : 1;
      const MInt size = surfaceNoNodes1D * surfaceNoNodes1D3 * noVars;

      // Skip h-refined surfaces since they have been already projected
      if(m_surfaces.fineCellId(srfcId) != -1 && m_surfaces.fineCellId(srfcId) != m_elements.cellId(elementId)) {
        continue;
      }

      // Calculate forward projection for lower polyDeg elements
      if(surfacePolyDeg > elementPolyDeg) {
        // Calculate forward projection
        calcMortarProjection<dg::mortar::forward, SysEqn::noVars()>(srfcId, dir, &m_surfaces.variables(srfcId, side),
                                                                    &projected[0], m_elements, m_surfaces);

        // Copy results of projection back to surface
        copy_n(&projected[0], size, &m_surfaces.variables(srfcId, side));
      }
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::ForwardProjection]);
}


/// \brief Calculate the volume integral for all elements and update
///        m_rightHandSide.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcVolumeIntegral() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::VolInt]);

  calcVolumeIntegral(m_elements.size(), m_elements, m_sysEqn);

  RECORD_TIMER_STOP(m_timers[Timers::VolInt]);
}


/**
 * \brief Calculate the numerical flux on the boundary surfaces and update
 *        m_flux.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] t Time at which the boundary conditions should be set (necessary
 *              for time-dependent b.c.'s).
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcBoundarySurfaceFlux(MFloat t) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Flux]);
  RECORD_TIMER_START(m_timers[Timers::FluxBndry]);

  for(auto&& bc : m_boundaryConditions) {
    bc->apply(t);
  }

  RECORD_TIMER_STOP(m_timers[Timers::FluxBndry]);
  RECORD_TIMER_STOP(m_timers[Timers::Flux]);
}


/**
 * \brief Calculate the numerical flux on the internal surfaces and update
 *        m_flux.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-24
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcInnerSurfaceFlux() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Flux]);
  RECORD_TIMER_START(m_timers[Timers::FluxInner]);

  const MInt begin = m_innerSurfacesOffset;
  const MInt end = m_innerSurfacesOffset + m_noInnerSurfaces;
  calcRegularSurfaceFlux(begin, end, m_surfaces, m_sysEqn);

  RECORD_TIMER_STOP(m_timers[Timers::FluxInner]);
  RECORD_TIMER_STOP(m_timers[Timers::Flux]);
}


/**
 * \brief Calculate the numerical flux on the MPI surfaces and update
 *        m_flux.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-24
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcMpiSurfaceFlux() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Flux]);
  RECORD_TIMER_START(m_timers[Timers::FluxMPI]);

  const MInt begin = m_mpiSurfacesOffset;
  const MInt end = m_mpiSurfacesOffset + m_noMpiSurfaces;
  calcRegularSurfaceFlux(begin, end, m_surfaces, m_sysEqn);

  RECORD_TIMER_STOP(m_timers[Timers::FluxMPI]);
  RECORD_TIMER_STOP(m_timers[Timers::Flux]);
}


/// \brief Calculate the surface integral for all faces of element and  update
///        dU/dt.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-02-07
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcSurfaceIntegral() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::SurfInt]);

  calcSurfaceIntegral(0, m_elements.size(), m_elements, m_surfaces, m_helements, m_helements.size());

  RECORD_TIMER_STOP(m_timers[Timers::SurfInt]);
}


/// \brief Calculate the surface integral for all faces of each element and
///        update dU/dt.
///
/// \author Michael Schlottke, Sven Berger
/// \date 2014-02-07
/// \param[in] begin Index of the first element to consider.
/// \param[in] end Index+1 of the last element to consider.
/// \param[in] elem Pointer to elements.
/// \param[in] surf Pointer to surfaces.
/// \param[in] helem Pointer to h-elements.
/// \param[in] noHElements Number of h-elements.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcSurfaceIntegral(const MInt begin,
                                                          const MInt end,
                                                          ElementCollector& elem,
                                                          SurfaceCollector& surf,
                                                          HElementCollector& helem,
                                                          const MInt noHElements) {
  TRACE();

  const MInt noVars = SysEqn::noVars();
  const MInt* surfaceIds = &elem.surfaceIds(0, 0);

  /////////////////////////////////////////////////////////////////////////////
  /// 1. Calculate surface integral for all surfaces except for the coarse
  /// element side of h-refined surfaces.
  /////////////////////////////////////////////////////////////////////////////

#ifdef _OPENMP
#pragma omp parallel for
#endif
  // Loop over all elements and calculate surface integrals
  for(MInt elementId = begin; elementId < end; elementId++) {
    const MInt surfaceIdOffset = elementId * 2 * nDim;

    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      // Calculate auxiliary variables for better readability
      const MInt srfcId = surfaceIds[surfaceIdOffset + dir];
      const MInt srfcPolyDeg = surf.polyDeg(srfcId);
      const MInt srfcNodes1D = surf.noNodes1D(srfcId);
      const MInt srfcNodes1D3 = (nDim == 3) ? srfcNodes1D : 1;
      const MInt polyDeg = elem.polyDeg(elementId);
      const MInt noNodes1D = elem.noNodes1D(elementId);
      MFloatTensor f(&surf.flux(srfcId), srfcNodes1D, srfcNodes1D3, noVars);

      const MInt side = 1 - dir % 2;

      // Skip coarse surfaces they need to be handled separately
      if(surf.fineCellId(srfcId) != -1 && surf.fineCellId(srfcId) != elem.cellId(elementId)) {
        continue;
      }

      // Pure p-refinement case (srfcPolyDeg is never small than polyDeg)
      if(srfcPolyDeg > polyDeg) {
        MFloatTensor fp(srfcNodes1D, srfcNodes1D3, noVars);

        // Calculate reverse projection
        calcMortarProjection<dg::mortar::reverse, SysEqn::noVars()>(srfcId, dir, &f[0], &fp[0], elem, surf);

        // Use projected flux for integration
        applySurfaceIntegral(&elem.rightHandSide(elementId), polyDeg, noNodes1D, srfcId, side, &fp[0], surf);

      } else {
        // Use original flux for integration
        applySurfaceIntegral(&elem.rightHandSide(elementId), polyDeg, noNodes1D, srfcId, side, &f[0], surf);
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// 2. Calculate surface integral for h-refined surfaces.
  /////////////////////////////////////////////////////////////////////////////
  if(noHElements > 0) {
    const MInt noDirs = 2 * nDim;
    const MInt noSurfs = 2 * (nDim - 1);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt hElementId = 0; hElementId < noHElements; hElementId++) {
      const MInt coarseElementId = helem.elementId(hElementId);
      const MInt polyDeg = elem.polyDeg(coarseElementId);
      const MInt noNodes1D = elem.noNodes1D(coarseElementId);

      for(MInt dir = 0; dir < noDirs; dir++) {
        const MInt coarseSrfcId = elem.surfaceIds(coarseElementId, dir);

        for(MInt pos = 0; pos < noSurfs; pos++) {
          const MInt hSrfcId = helem.hrefSurfaceIds(hElementId, dir, pos);
          const MInt side = 1 - dir % 2;

          // Skip if this surface does not exist
          if(hSrfcId == -1) {
            continue;
          }

          const MInt srfcNodes1D = surf.noNodes1D(hSrfcId);
          const MInt srfcNodes1D3 = (nDim == 3) ? srfcNodes1D : 1;
          MFloatTensor f(&surf.flux(hSrfcId), srfcNodes1D, srfcNodes1D3, noVars);
          MFloatTensor fp(srfcNodes1D, srfcNodes1D3, noVars);

          // Calculate reverse projection
          calcMortarProjection<dg::mortar::reverse, SysEqn::noVars()>(hSrfcId, dir, &f[0], &fp[0], elem, surf);

          // Use projected flux for integration
          applySurfaceIntegral(&elem.rightHandSide(coarseElementId), polyDeg, noNodes1D, coarseSrfcId, side, &fp[0],
                               surf);
        }
      }
    }
  }
}


/// \brief Calculate the surface integral for a face of an element and
///        update dU/dt.
///
/// \author Michael Schlottke, Sven Berger
/// \date March 2015
/// \param[in] rhs Time derivatives of element variables where the surface
///                integral is added to.
/// \parma[in] polyDeg Polynomial degree of element.
/// \param[in] srfcId Surface id from which the flux is to be used.
/// \param[in] side Position of element with respect to the surface. A value of
///                 zero means the element is in the -ve coordinate direction
///                 with respect to the surface, a value of one means the
///                 element is in the +ve direction.
/// \param[in] flux The surface flux.
/// \param[in] surf Pointer to surfaces.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::applySurfaceIntegral(MFloat* rhs, const MInt polyDeg, const MInt noNodes1D,
                                                           const MInt srfcId, const MInt side, const MFloat* flux,
                                                           SurfaceCollector& surf) {
  // TRACE();
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MInt noVars = SysEqn::noVars();
  const MInt srfcNodes1D = surf.noNodes1D(srfcId);
  const MInt srfcNodes1D3 = (nDim == 3) ? srfcNodes1D : 1;

  const MFloatTensor f(const_cast<MFloat*>(flux), srfcNodes1D, srfcNodes1D3, noVars);

  // Calculate surface integral
  MFloatTensor ut(rhs, noNodes1D, noNodes1D, noNodes1D3, noVars);
  const MInt dirId = surf.orientation(srfcId);
  const MInt index = (side == 0) ? (noNodes1D - 1) : 0;
  const MInt sign = (side == 0) ? 1 : -1;
  const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];


  if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS) {
    // Use different loops depending on the surface orientation
    switch(dirId) {
      case 0:
        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D; j++) {
            for(MInt k = 0; k < noNodes1D3; k++) {
              for(MInt n = 0; n < noVars; n++) {
                ut(i, j, k, n) += sign * f(j, k, n) * interp.m_LhatFace[1 - side][i];
              }
            }
          }
        }
        break;

      case 1:
        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D; j++) {
            for(MInt k = 0; k < noNodes1D3; k++) {
              for(MInt n = 0; n < noVars; n++) {
                ut(i, j, k, n) += sign * f(i, k, n) * interp.m_LhatFace[1 - side][j];
              }
            }
          }
        }
        break;

      case 2:
        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D; j++) {
            for(MInt k = 0; k < noNodes1D3; k++) {
              for(MInt n = 0; n < noVars; n++) {
                ut(i, j, k, n) += sign * f(i, j, n) * interp.m_LhatFace[1 - side][k];
              }
            }
          }
        }
        break;
      default:
        mTerm(1, AT_, "Bad direction id");
    }
  } else if(m_dgIntegrationMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
    // Use different loops depending on the surface orientation
    switch(dirId) {
      case 0:
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt n = 0; n < noVars; n++) {
              ut(index, j, k, n) += sign * f(j, k, n) * interp.m_LhatFace[1 - side][index];
            }
          }
        }
        break;

      case 1:
        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt n = 0; n < noVars; n++) {
              ut(i, index, k, n) += sign * f(i, k, n) * interp.m_LhatFace[1 - side][index];
            }
          }
        }
        break;

      case 2:
        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D3; j++) {
            for(MInt n = 0; n < noVars; n++) {
              ut(i, j, index, n) += sign * f(i, j, n) * interp.m_LhatFace[1 - side][index];
            }
          }
        }
        break;

      default:
        mTerm(1, AT_, "Bad direction id");
    }
  }
}


/// \brief Adds the negative of the inverse Jacobian to the time derivative.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::applyJacobian() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Jacobian]);

  applyJacobian(m_elements.size(), m_elements);

  RECORD_TIMER_STOP(m_timers[Timers::Jacobian]);
}


/**
 * \brief Adds the negative of the inverse Jacobian to the time derivative.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-23
 *
 * \param[in] noElements Number of elements.
 * \param[in] elem Pointer to elements.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::applyJacobian(const MInt noElements, ElementCollector& elem) {
  TRACE();

  MFloat* const invJacobians = &elem.invJacobian(0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt dataBlockSize = elem.noNodesXD(elementId) * SysEqn::noVars();
    const MFloat invJacobian = invJacobians[elementId];

    MFloat* const rhs = &elem.rightHandSide(elementId);
    for(MInt dataId = 0; dataId < dataBlockSize; dataId++) {
      rhs[dataId] *= -invJacobian;
    }
  }
}


/// \brief Calculates the source terms for each node and adds them to the time
///        derivative of the conservative variables.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcSourceTerms(MFloat t) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Sources]);

  calcSourceTerms(t, m_elements.size(), m_elements, m_sysEqn);

  RECORD_TIMER_STOP(m_timers[Timers::Sources]);
}


/// \brief Add the external coupling source terms to the right hand side
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::applyExternalSourceTerms(const MFloat time) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::ExternalSources]);

  const MInt* noNodes1D = &m_elements.noNodes1D(0);

  // Source term ramp up factor (linear in time)
  const MFloat rampUpFactor =
      (m_useSourceRampUp) ? maia::filter::slope::linear(m_startTime, m_startTime + m_sourceRampUpTime, time) : 1.0;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt elementId = 0; elementId < m_elements.size(); elementId++) {
    const MInt dataBlockSize = ipow(noNodes1D[elementId], nDim) * SysEqn::noVars();
    MFloat* const rhs = &m_elements.rightHandSide(elementId);
    const MFloat* const sources = &m_elements.externalSource(elementId);

    for(MInt dataId = 0; dataId < dataBlockSize; dataId++) {
      rhs[dataId] += rampUpFactor * sources[dataId];
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::ExternalSources]);
}


/// \brief Return true if this is the root rank of the solver MPI communicator.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-01-13
///
/// \return True if the rank is zero in the solver MPI communicator.
template <MInt nDim, class SysEqn>
inline MBool DgCartesianSolver<nDim, SysEqn>::isMpiRoot() const {
  return (domainId() == 0);
}


/**
 * \brief Start sending the window-side data and start receiving the halo-side
 *        data for all MPI surfaces.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-24
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::startMpiSurfaceExchange() {
  TRACE();

  // IMPORTANT:
  // If anything in this method is changed, please check if
  // updateNodeVariables() needs to be changed as well, since it contains more
  // or less the same code.

  RECORD_TIMER_START(m_timers[Timers::SurfExchange]);
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::MPI]);

  // Start receiving
  RECORD_TIMER_START(m_timers[Timers::MPIComm]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeComm]);
  RECORD_TIMER_START(m_timers[Timers::SECommRecv]);
  if(!m_mpiRecvRequestsOpen && m_noExchangeNghbrDomains > 0) {
    MPI_Startall(m_noExchangeNghbrDomains, &m_recvRequests[0], AT_);
    m_mpiRecvRequestsOpen = true;
  }
  RECORD_TIMER_STOP(m_timers[Timers::SECommRecv]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeComm]);
  RECORD_TIMER_STOP(m_timers[Timers::MPIComm]);

  // Finish previous sending
  RECORD_TIMER_START(m_timers[Timers::MPIWait]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_START(m_timers[Timers::SEWaitSend]);
  if(m_mpiSendRequestsOpen && m_noExchangeNghbrDomains > 0) {
    MPI_Waitall(m_noExchangeNghbrDomains, &m_sendRequests[0], MPI_STATUSES_IGNORE, AT_);
    m_mpiSendRequestsOpen = false;
  }
  RECORD_TIMER_STOP(m_timers[Timers::SEWaitSend]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_STOP(m_timers[Timers::MPIWait]);

  // Pack send buffers
#ifdef DG_USE_MPI_BUFFERS
  RECORD_TIMER_START(m_timers[Timers::MPICopy]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeCopy]);
  RECORD_TIMER_START(m_timers[Timers::SECopySend]);

  // Use max. polynomial degree to account for p-refined surfaces
  const MInt dataBlockSize = ipow(m_maxNoNodes1D, nDim - 1) * SysEqn::noVars();

  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MInt size = 0;
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      const MInt srfcId = m_mpiSurfaces[i][j];
      const MInt sideId = m_surfaces.internalSideId(srfcId);

      // Copy state data
      const MFloat* data = &m_surfaces.variables(srfcId, sideId);
      copy(data, data + dataBlockSize, &m_sendBuffers[i][size]);
      size += dataBlockSize;
    }
    ASSERT(size == static_cast<MInt>(m_sendBuffers[i].size()), "Data size does not match buffer size.");
  }
  RECORD_TIMER_STOP(m_timers[Timers::SECopySend]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeCopy]);
  RECORD_TIMER_STOP(m_timers[Timers::MPICopy]);
#endif

  // Start sending
  RECORD_TIMER_START(m_timers[Timers::MPIComm]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeComm]);
  RECORD_TIMER_START(m_timers[Timers::SECommSend]);
  if(m_noExchangeNghbrDomains > 0) {
    MPI_Startall(m_noExchangeNghbrDomains, &m_sendRequests[0], AT_);
  }
  RECORD_TIMER_STOP(m_timers[Timers::SECommSend]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeComm]);
  RECORD_TIMER_STOP(m_timers[Timers::MPIComm]);

  m_mpiSendRequestsOpen = true;

  RECORD_TIMER_STOP(m_timers[Timers::MPI]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchange]);
}


/**
 * \brief Finish sending the window-side data and finish receiving the halo-side
 *        data for all MPI surfaces.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-07-24
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::finishMpiSurfaceExchange() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::SurfExchange]);
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::MPI]);

  // IMPORTANT:
  // If anything in this method is changed, please check if
  // updateNodeVariables() needs to be changed as well, since it contains more
  // or less the same code.

  if(noDomains() > 1 && (!m_mpiRecvRequestsOpen || !m_mpiSendRequestsOpen)) {
    TERMM(1, "MPI requests are not open: receive=" + std::to_string(m_mpiRecvRequestsOpen)
                 + ", send=" + std::to_string(m_mpiSendRequestsOpen));
  }

#ifdef DG_USE_MPI_WAITSOME
  TERMM(1, "This code is untested and will not work with p-refinement.");

  // Init status array and counter
  vector<MInt> indices(m_noExchangeNghbrDomains);
  MInt noFinished;

  // Start with receiving data and unpack the buffers as they are filled
  // with incoming data
  RECORD_TIMER_START(m_timers[Timers::MPIWait]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_START(m_timers[Timers::SEWaitRecv]);
  MPI_Waitsome(m_noExchangeNghbrDomains, &m_recvRequests[0], &noFinished, &indices[0], MPI_STATUSES_IGNORE, AT_);
  while(noFinished != MPI_UNDEFINED) {
    for(MInt i = 0; i < noFinished; i++) {
      const MInt index = indices[i];
      MInt size = 0;
      for(size_t j = 0; j < m_mpiSurfaces[index].size(); j++) {
        const MInt srfcId = m_mpiSurfaces[index][j];
        const MInt noNodes1D = m_surfaces.noNodes1D(srfcId);
        const MInt noNodesXD = ipow(noNodes1D, nDim - 1);
        const MInt sideId = 1 - m_surfaces.internalSideId(srfcId);

        MFloat* data = m_surfaces.variables(srfcId, sideId);
        const MInt dataBlockSize = noNodesXD * SysEqn::noVars();
        copy_n(&m_recvBuffers[index][size], dataBlockSize, data);
        size += dataBlockSize;
      }
      ASSERT(size == static_cast<MInt>(m_recvBuffers[index].size()), "Data size does not match buffer size.");
    }

    // Wait for next batch
    MPI_Waitsome(m_noExchangeNghbrDomains, &m_recvRequests[0], &noFinished, &indices[0], MPI_STATUSES_IGNORE, AT_);
  }
  RECORD_TIMER_STOP(m_timers[Timers::SEWaitRecv]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_STOP(m_timers[Timers::MPIWait]);
#else

  // Finish receiving
  RECORD_TIMER_START(m_timers[Timers::MPIWait]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_START(m_timers[Timers::SEWaitRecv]);
  if(m_noExchangeNghbrDomains > 0) {
    MPI_Waitall(m_noExchangeNghbrDomains, &m_recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  }
  RECORD_TIMER_STOP(m_timers[Timers::SEWaitRecv]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeWait]);
  RECORD_TIMER_STOP(m_timers[Timers::MPIWait]);

  // Unpack receive buffers
#ifdef DG_USE_MPI_BUFFERS
  RECORD_TIMER_START(m_timers[Timers::MPICopy]);
  RECORD_TIMER_START(m_timers[Timers::SurfExchangeCopy]);
  RECORD_TIMER_START(m_timers[Timers::SECopyRecv]);

  // Use max. polynomial degree to account for p-refined surfaces
  const MInt dataBlockSize = ipow(m_maxNoNodes1D, nDim - 1) * SysEqn::noVars();

  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MInt size = 0;
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      // Here the surface id from m_mpiRecvSurfaces is used to handle the case of communication
      // of a domain with itself due to periodically connected boundaries
      const MInt srfcId = m_mpiRecvSurfaces[i][j];


      // Use max. polynomial degree to account for p-refined surfaces
      const MInt sideId = 1 - m_surfaces.internalSideId(srfcId);
      MFloat* data = &m_surfaces.variables(srfcId, sideId);
      copy_n(&m_recvBuffers[i][size], dataBlockSize, data);
      size += dataBlockSize;
    }
    ASSERT(size == static_cast<MInt>(m_recvBuffers[i].size()), "Data size does not match buffer size.");
  }
  RECORD_TIMER_STOP(m_timers[Timers::SECopyRecv]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchangeCopy]);
  RECORD_TIMER_STOP(m_timers[Timers::MPICopy]);
#endif // DG_USE_MPI_BUFFERS
#endif // DG_USE_MPI_WAITSOME

  m_mpiRecvRequestsOpen = false;

  // Start new receive requests
  if(m_noExchangeNghbrDomains > 0) {
    MPI_Startall(m_noExchangeNghbrDomains, &m_recvRequests[0], AT_);
    m_mpiRecvRequestsOpen = true;
  }

  RECORD_TIMER_STOP(m_timers[Timers::MPI]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
  RECORD_TIMER_STOP(m_timers[Timers::SurfExchange]);
}


/// \brief Calculate the next time step.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2013-04-26
///
/// \return The new time step size.
template <MInt nDim, class SysEqn>
MFloat DgCartesianSolver<nDim, SysEqn>::calcTimeStep() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcTimeStep]);

  MFloat minDt = numeric_limits<MFloat>::infinity();

  // First check if external coupling overrides the time step calculation
  const MBool forcedTimeStep = (m_externalDt > 0.0);

  // Time step calculation (if forced to check that the external dt is small enough)
  {
    // Calculate time step from system of equations class
    const MInt noElements = m_elements.size();
    const MInt* noNodes1D = &m_elements.noNodes1D(0);
    const MFloat* invJacobians = &m_elements.invJacobian(0);

#ifdef _OPENMP
#pragma omp parallel for reduction(min : minDt)
#endif
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      MFloat ts = m_sysEqn.getTimeStep(&m_elements.nodeVars(elementId), &m_elements.variables(elementId),
                                       noNodes1D[elementId], invJacobians[elementId], m_sbpMode);
      minDt = min(minDt, ts);
    }
  }

  if(forcedTimeStep) {
    if(m_externalDt > minDt) {
      cerr << " WARNING: external dt larger than computed time step size on global domain " << globalDomainId() << "("
           << m_externalDt << " > " << minDt << ")" << std::endl;
    }
    // Calculate global minimum computed time step
    MPI_Allreduce(MPI_IN_PLACE, &minDt, 1, type_traits<MFloat>::mpiType(), MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                  "minDt");
    m_log << "DG-Solver #" << solverId() << " using external time step size: " << m_externalDt
          << " (computed global minimum: " << minDt << ")" << std::endl;
    minDt = m_externalDt;
  }

  // Check for NaN in time step
  if(std::isnan(minDt)) {
    TERMM(1,
          "Time step is NaN at time step " + to_string(m_timeStep) + " on global domain "
              + to_string(globalDomainId()));
  }

  // Check negative time step
  if(minDt < 0.0) {
    TERMM(1,
          "Time step is less than zero at time step " + to_string(m_timeStep) + " on global domain "
              + to_string(globalDomainId()));
  }

  // Check ridiculuously small time step that indicates *most probably* an error
  // Note: the threshold here was arbitrarily chosen, but seemed "low enough"
  if(minDt < 1.0e-100) {
    TERMM(1,
          "Time step is less than 1.0e-100 at time step " + to_string(m_timeStep) + " on global domain "
              + to_string(globalDomainId()));
  }

  // Calculate global minimum time step
  MPI_Allreduce(MPI_IN_PLACE, &minDt, 1, type_traits<MFloat>::mpiType(), MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                "minDt");

  RECORD_TIMER_STOP(m_timers[Timers::CalcTimeStep]);

  return minDt;
}


/**
 * \brief Time integration using the five-stage fourth-order low-storage
 *        Runge-Kutta scheme as described in the book by
 *        Hesthaven and Warburton (2008), p. 64, Eqn. 3.5.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-26
 *
 * \param[in] t Current physical time (needed for boundary conditions etc.).
 * \param[in] dt Current time step size.
 * \param[in] substep Perform only the specified Runge-Kutta substep/stage.
 *
 * The 'substep' option allows to perform only a single Runge-Kutta substep,
 * which is necessary for interleaving the Runge-Kutta stages in a coupled
 * simulation to avoid idle time. By default the complete timestep is performed.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::timeStepRk(const MFloat t, const MFloat dt, const MInt substep) {
  TRACE();

  RECORD_TIMER_START(m_timers[Timers::RungeKuttaStep]);

  // Save total data length
  const MInt totalSize = m_internalDataSize;

  // General Runge-Kutta time step calculation
  for(MInt s = 0; s < m_noRkStages; s++) {
    // Perform only the given single Runge-Kutta substep if specified
    if(substep != -1 && s != substep) {
      continue;
    }

    calcDgTimeDerivative(t, t + m_timeIntegrationCoefficientsC[s] * dt);

    RECORD_TIMER_START(m_timers[Timers::TimeInt]);
    subTimeStepRk(dt, s, totalSize, &m_elements.rightHandSide(0), &m_elements.variables(0),
                  &m_elements.timeIntStorage(0));
    RECORD_TIMER_STOP(m_timers[Timers::TimeInt]);


    // Determine if this is the final Runge-Kutta stage before an adaptation
    // time step (error estimates are calculated with the current solution on
    // all surfaces, so prevent overwriting the surface states if split MPI
    // communication is active)
    const MBool lastRkStage = (s == (m_noRkStages - 1));
    const MBool adaptationTimeStep = (lastRkStage && isAdaptationTimeStep(m_timeStep + 1));

    // If multiphysics-optimized parallelization is used, prolong and start
    // tranmitting
    if(g_splitMpiComm && !adaptationTimeStep) {
      RECORD_TIMER_START(m_timers[Timers::TimeDeriv]);
      prolongToSurfaces();
      applyForwardProjection();
      startMpiSurfaceExchange();
      RECORD_TIMER_STOP(m_timers[Timers::TimeDeriv]);
    }
  }


  RECORD_TIMER_STOP(m_timers[Timers::RungeKuttaStep]);
}


/// \brief Perform one Runge-Kutta substep on the given elements.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// \param[in] dt Current time step size.
/// \param[in] stage Runge-Kutta stage to perform.
/// \param[in] totalSize Total data length.
/// \param[in] rhs Pointer to right hand side, i.e. time derivative.
/// \param[in] variables Pointer to variables.
/// \param[in] timeIntStorage Pointer to time integration storage.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::subTimeStepRk(const MFloat dt,
                                                    const MInt stage,
                                                    const MInt totalSize,
                                                    const MFloat* const rhs,
                                                    MFloat* const variables,
                                                    MFloat* const timeIntStorage) {
  TRACE();

  // Get pointers for a more concise code
  MFloat* const p = variables;
  MFloat* const k = timeIntStorage;

  // Calculate auxiliary variables to save on operations
  MFloat bDt = -1.0;
  bDt = m_timeIntegrationCoefficientsB[stage] * dt;

  // Stage 0
  if(stage == 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt i = 0; i < totalSize; i++) {
      k[i] = rhs[i];
      p[i] += k[i] * bDt;
    }
  } else {
    // Stage 1-End
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt i = 0; i < totalSize; i++) {
      k[i] = rhs[i] - k[i] * m_timeIntegrationCoefficientsA[stage];
      p[i] += k[i] * bDt;
    }
  }
}


/**
 * \brief Calculates and prints the L^2 and L^inf errors
 *        (using calcErrorNorms()) for the current time.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-05-01
 *
 * \param[in] t Time at which to analyze the solution.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::analyzeTimeStep(MFloat t, MFloat runTimeRelative, MFloat runTimeTotal,
                                                      MInt timeStep, MFloat dt) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::AnalyzeTimeStep]);
  // Calculate error norms only if enabled in properties
  vector<MFloat> L2Error(m_sysEqn.noVars());
  vector<MFloat> LInfError(m_sysEqn.noVars());
  vector<MFloat> L2ErrLocal(m_sysEqn.noVars());
  vector<MFloat> LInfErrLocal(m_sysEqn.noVars());
  if(m_calcErrorNorms) {
    calcErrorNorms(t, L2Error, LInfError, L2ErrLocal, LInfErrLocal);
    // If any of the error measures are NaN, write a solution file before the
    // simulation is aborted
    if(any_of(L2Error.begin(), L2Error.end(), [](const MFloat e) { return std::isnan(e); })
       || any_of(LInfError.begin(), LInfError.end(), [](const MFloat e) { return std::isnan(e); })) {
      saveSolutionFile("nan_" + to_string(timeStep));
    }

    // Check error norms for NaN and abort if found
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      if(std::isnan(L2ErrLocal[i])) {
        TERMM(1, "L^2 error for variable '" + m_sysEqn.consVarNames(i) + "' is NaN at time step " + to_string(timeStep)
                     + " on global domain " + to_string(globalDomainId()));
      }
      if(std::isnan(LInfErrLocal[i])) {
        TERMM(1, "L^inf error for variable '" + m_sysEqn.consVarNames(i) + "' is NaN at time step "
                     + to_string(timeStep) + " on global domain " + to_string(globalDomainId()));
      }
    }
  }

  const MInt precision = m_noErrorDigits;
  const MInt fieldwith = precision + 6;

  if(isMpiRoot()) {
    stringstream log;

    log << endl;
    log << "----------------------------------------"
        << "----------------------------------------" << endl;
    log << " Solver " << solverId() << " running '" << m_sysEqn.sysEqnName() << "' with N = " << m_initPolyDeg
        << " and maxLevel = " << grid().maxLevel() << endl;
    log << "----------------------------------------"
        << "----------------------------------------" << endl;
    log << " No. timesteps: " << timeStep << endl;
    log << " dt:            " << scientific << dt << endl;
    log << " Run time:      " << runTimeTotal << " s" << endl;
    log << " Time/DOF/step: " << runTimeRelative << " s" << endl;

    if(m_calcErrorNorms) {
      streamsize ss = log.precision();
      log << setprecision(precision);
      log << " Variable:    ";
      for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
        log << "  " << left << setw(fieldwith) << setfill(' ') << m_sysEqn.consVarNames(i);
      }
      log << right << endl;
      log << " L^2 error:   ";
      for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
        log << "  " << L2Error[i];
      }
      log << endl;
      log << " L^inf error: ";
      for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
        log << "  " << LInfError[i];
      }
      log << endl;
      log << setprecision(ss);
    }

    log << "----------------------------------------"
        << "----------------------------------------" << endl;
    log << " Simulation time: " << t << endl;
    log << "----------------------------------------"
        << "----------------------------------------" << endl;

    m_log << log.str() << endl;
    cout << log.str() << endl;

    // Determine if this is the last time step and reset time step if necessary
    MBool finalTimeStep = false;
    if(m_finalTime - m_time - m_dt < 1.0E-10) {
      finalTimeStep = true;
    } else if(m_timeStep == m_timeSteps) {
      finalTimeStep = true;
    }

    if(m_calcErrorNorms && isMpiRoot() && finalTimeStep) {
      std::ofstream eocStream;
      eocStream.open("EOCout.csv");
      eocStream << m_statGlobalNoActiveDOFs << "," << L2Error[0] << "," << LInfError[0];
      eocStream.close();

      std::ofstream perfStream;
      perfStream.open("TimerOut.csv");
      perfStream << m_statGlobalNoActiveCells << "," << m_timeStep << "," << runTimeTotal << "," << runTimeRelative;
      perfStream.close();
    }
  }

  // write local errors
  if(m_calcErrorNorms) {
    stringstream logLocal;
    logLocal << setprecision(precision) << scientific;
    logLocal << endl << " Variable:         ";
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      logLocal << "  " << left << setw(fieldwith) << setfill(' ') << m_sysEqn.consVarNames(i);
    }
    logLocal << right << endl;
    logLocal << " Local L^2 error:  ";
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      logLocal << "  " << L2ErrLocal[i];
    }
    logLocal << endl;
    logLocal << " Local L^inf error:";
    for(MInt i = 0; i < m_sysEqn.noVars(); i++) {
      logLocal << "  " << LInfErrLocal[i];
    }
    logLocal << endl;

    m_log << logLocal.str() << endl;
  }

  RECORD_TIMER_STOP(m_timers[Timers::AnalyzeTimeStep]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
}


/**
 * \brief Calculate the L^2 and L^infinity error norms at a given time. The
 *        error is calculated in comparison to the initial condition.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-05-01
 *
 * \param[in] t Physical time at which to calculate the error norms.
 * \param[out] L2Error Storage for the L^2 error (one value/variable).
 * \param[out] LInfError Storage for the L^infinity error (one value/variable).
 * \param[out] L2ErrLocal Storage for the local L^2 error (one value/variable).
 * \param[out] LInfErrLocal Storage for the local L^infinity error (one
 *             value/variable).
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcErrorNorms(const MFloat t,
                                                     vector<MFloat>& L2Error,
                                                     vector<MFloat>& LInfError,
                                                     vector<MFloat>& L2ErrLocal,
                                                     vector<MFloat>& LInfErrLocal) {
  TRACE();
  const MInt noVars = m_sysEqn.noVars();
  L2Error.assign(noVars, F0);
  LInfError.assign(noVars, F0);
  L2ErrLocal.assign(noVars, F0);
  LInfErrLocal.assign(noVars, F0);

  // Allocate space for temporary values
  MFloatScratchSpace uExact(noVars, FUN_, "uExact");
  const MInt noNodesAnalysis1D = m_noAnalysisNodes;
  const MInt noNodesAnalysis1D3 = (nDim == 3) ? noNodesAnalysis1D : 1;
  MFloatTensor u(noNodesAnalysis1D, noNodesAnalysis1D, noNodesAnalysis1D3, noVars);
  MFloatTensor x(noNodesAnalysis1D, noNodesAnalysis1D, noNodesAnalysis1D3, nDim);
  // Set node vars to minimum 1 to avoid errors in Tensor
  // TODO labels:DG Check if this is really sensible
  MFloatTensor nodeVars(noNodesAnalysis1D, noNodesAnalysis1D, noNodesAnalysis1D3, max(SysEqn::noNodeVars(), 1));

  // Loop over internal elements and calculate local errors
  const MInt noElements = m_elements.size();
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    // Interpolate solution from regular nodes to analysis nodes
    const MInt polyDegElement = m_elements.polyDeg(elementId);
    const MInt noNodesElement = m_elements.noNodes1D(elementId);
    auto vdm = m_vdmAnalysis[polyDegElement][noNodesElement];
    dg::interpolation::interpolateNodes<nDim>(&m_elements.variables(elementId), &vdm[0], noNodesElement,
                                              m_noAnalysisNodes, noVars, &u[0]);

    // Interpolation node locations from regular nodes to analysis nodes
    dg::interpolation::interpolateNodes<nDim>(&m_elements.nodeCoordinates(elementId), &vdm[0], noNodesElement,
                                              m_noAnalysisNodes, nDim, &x[0]);

    // Calculate error norms
    const MFloat jacobian = pow(F1 / m_elements.invJacobian(elementId), nDim);
    for(MInt i = 0; i < noNodesAnalysis1D; i++) {
      for(MInt j = 0; j < noNodesAnalysis1D; j++) {
        for(MInt k = 0; k < noNodesAnalysis1D3; k++) {
          // TODO labels:DG Check if we need to initialize coupling quantities here as
          //      well (instead of just interpolating)
          m_sysEqn.calcInitialCondition(t, &x(i, j, k, 0), &nodeVars(i, j, k, 0), &uExact[0]);
          for(MInt v = 0; v < noVars; v++) {
            const MFloat diff = uExact[v] - u(i, j, k, v);
            L2ErrLocal[v] += pow(diff, F2) * m_wVolumeAnalysis(i, j, k) * jacobian;
            LInfErrLocal[v] = max(LInfErrLocal[v], fabs(diff));
          }
        }
      }
    }
  }

  // Calculate global errors
  RECORD_TIMER_START(m_timers[Timers::MPI]);
  RECORD_TIMER_START(m_timers[Timers::MPIComm]);
  MPI_Allreduce(&L2ErrLocal[0], &L2Error[0], noVars, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "L2ErrLocal[0]",
                "L2Error[0]");
  MPI_Allreduce(&LInfErrLocal[0], &LInfError[0], noVars, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "LInfErrLocal[0]",
                "LInfError[0]");
  RECORD_TIMER_STOP(m_timers[Timers::MPIComm]);
  RECORD_TIMER_STOP(m_timers[Timers::MPI]);

  // Finalize L^2 error calculation
  for(MInt v = 0; v < noVars; v++) {
    L2ErrLocal[v] = sqrt(L2ErrLocal[v] / m_localVolume);
    L2Error[v] = sqrt(L2Error[v] / m_globalVolume);
  }
}


/// \brief Reset the time derivative of the conservative variables to zero.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::resetRHS() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::ResetRHS]);

  resetBuffer(m_internalDataSize, &m_elements.rightHandSide(0));

  RECORD_TIMER_STOP(m_timers[Timers::ResetRHS]);
}


template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::resetExternalSources() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::ResetExternalSources]);

  resetBuffer(m_internalDataSize, &m_elements.externalSource(0));

  RECORD_TIMER_STOP(m_timers[Timers::ResetExternalSources]);
}


/**
 * \brief Reset the given buffer to zero.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-23
 *
 * \param[in] totalSize Total data size.
 * \param[in] buffer Pointer to data buffer to be reset.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::resetBuffer(const MInt totalSize, MFloat* const buffer) {
  TRACE();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt i = 0; i < totalSize; i++) {
    buffer[i] = 0.0;
  }
}


/**
 *  \brief returns the elementId of a element containing a given point
 *
 *  \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>
 *  \date 24.01.2014 (update 01.08.17)
 *
 *  \param[in] point Evaluation point.
 *  \param[in] globalUnique Specify if a point should be globally unique (i.e. only appears once).
 *  \param[out] returns elementId.
 *
 *  If the parameter "globalUnique" is set to true it is assured that every point is globally
 *  associated to only one element. This allows to avoid problems for points located on a domain
 *  boundary that should only appear on a single domain, e.g. when writing point data. By default
 *  globalUnique is set to false.
 */
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::getElementIdAtPoint(const MFloat* point, MBool globalUnique) {
  TRACE();

  MInt foundElementId = -1;
  const MInt noElements = m_elements.size();

  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    const MFloat cellLength = grid().cellLengthAtCell(cellId);
    MBool pointInCell = true;
    MBool takeCell = true;

    for(MInt i = 0; i < nDim; i++) {
      if(!(fabs(grid().tree().coordinate(cellId, i) - point[i]) <= cellLength * 0.5)) {
        pointInCell = false;
      }
    }

    // Set globalUnique to true in order to avoid that points located on domain boundaries are
    // found multiple times (i.e. on different domains). This ensures that each point is globally
    // unique
    if(pointInCell && globalUnique) {
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        const MFloat distance = grid().tree().coordinate(cellId, dimId) - point[dimId];
        // Check if the point is located on a cell edge
        if(approx(fabs(distance), cellLength * 0.5, MFloatEps)) {
          // Check the relative position of the point with respect to the cell center
          if(distance > 0.0) {
            // Point is on surface in negative coordinate direction, thus neighborDir is either
            // equal to 0, 2, or 4
            const MInt neighborDir = 2 * dimId;
            const MInt surfaceId = m_elements.surfaceIds(elementId, neighborDir);
            // Check for a MPI surface, if this is the case the unique point belongs to the
            // neighboring element on another domain
            if(isMpiSurface(surfaceId)) {
              takeCell = false;
            }
          }
        }
      }
    }

    if(pointInCell && takeCell) {
      foundElementId = elementId;
      break;
    }
  }

  return foundElementId;
}


/**
 *  \brief returns the cellId of a cell containing a given point
 *
 *  \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>
 *  \date 24.01.2014
 *
 *  \param[in] point Evaluation point.
 *  \param[out] returns cellId.
 *
 */
// template <MInt nDim, class SysEqn>
// MInt DgCartesianSolver<nDim, SysEqn>::getCellIdAtPoint(const MFloat* point) {
//  TRACE();
//
//  const MInt elementId = getElementIdAtPoint(point);
//
//  return (elementId == -1) ? -1 : m_elements.cellId(elementId);
//}

/**
 *  \brief Calculates the state vector at a given Point.
 *
 *  \author Fabian Klemp <f.klemp@aia.rwth-aachen.de>
 *  \date 2016-07-22
 *
 *  \param[in] point Evaluation point.
 *  \param[out] state State vector at a given point.
 *  \details state should have the dimension SysEqn::noVars();
 */
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::calcStateAtPoint(const MFloat* const point, MFloat* const state) {
  TRACE();
  const MInt elementId = getElementIdAtPoint(point);

  if(elementId == -1) {
    return false;
  }

  calcStateAtPoint(point, elementId, state);

  return true;
}


/**
 *  \brief Calculates the state vector at a given Point.
 *
 *  \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>
 *  \date 24.01.2014
 *
 *  \param[in] point Evaluation point.
 *  \param[in] elementId Element ID of the evaluation point.
 *  \param[in] noVars Number of variables of given element data field.
 *  \param[in] u Element data field to evaluate.
 *  \param[out] state State vector at a given point.
 *  \details state should have the dimension 'noVars'.
 */
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcStateAtPoint(const MFloat* const point, const MInt elementId,
                                                       const MInt noVars, const MFloat* const u, MFloat* const state) {
  TRACE();

  const MInt cellId = m_elements.cellId(elementId);
  const MInt polyDeg = m_elements.polyDeg(elementId);
  const MInt noNodes1D = m_elements.noNodes1D(elementId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];

  const MFloat cellLength = grid().cellLengthAtCell(cellId);

  // Vector saves the point with coordinates on the interval [17,1]
  MFloat pointOnUnitInt[nDim];

  const MFloatTensor U(const_cast<MFloat*>(u), noNodes1D, noNodes1D, noNodes1D3, noVars);

  // Because Lagrangefunctions are only defined on the interval [-1,1]
  // the coordinates are transformed to the interval.
  for(MInt n = 0; n < nDim; n++) {
    pointOnUnitInt[n] = F2 * (point[n] - grid().tree().coordinate(cellId, n)) / cellLength;
  }

  if(m_sbpMode) {
    IF_CONSTEXPR(nDim == 2) {
      dg::interpolation::calcBilinearInterpolation(pointOnUnitInt, &interp.m_nodes[0], noNodes1D, noVars, u, state);
    }
    else {
      dg::interpolation::calcTrilinearInterpolation(pointOnUnitInt, &interp.m_nodes[0], noNodes1D, noVars, u, state);
    }
  } else {
    /*
     * Transform the coordinates of the Point to the intervall [-1,1]:
     *  P(x,y,z) --> P(xi, eta, mu)
     *
     * To calculate the state Q at a Point P we use the interpolation formula:
     *  Q(x,y,z) = Q(x(xi), y(eta), z(mu))
     *           = \Sum_{i,j,k}^{N} q_{ijk} * l_i(xi) * l_j(eta) * l_k(mu)
     **/

    MFloatVector lagrangePoly[nDim];
    for(MInt n = 0; n < nDim; n++) {
      lagrangePoly[n].resize(noNodes1D);
      fill_n(&lagrangePoly[n][0], noNodes1D, F0);
    }

    for(MInt n = 0; n < nDim; n++) {
      dg::interpolation::calcLagrangeInterpolatingPolynomials(pointOnUnitInt[n], polyDeg, &interp.m_nodes[0],
                                                              &interp.m_wBary[0], &lagrangePoly[n][0]);
    }

    fill_n(state, noVars, 0.0);

    IF_CONSTEXPR(nDim == 2) {
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt v = 0; v < noVars; v++) {
              state[v] += U(i, j, k, v) * lagrangePoly[0][i] * lagrangePoly[1][j];
            }
          }
        }
      }
    }
    else { // nDim == 3
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt v = 0; v < noVars; v++) {
              state[v] += U(i, j, k, v) * lagrangePoly[0][i] * lagrangePoly[1][j] * lagrangePoly[2][k];
            }
          }
        }
      }
    }
  }
}

///// \brief Calculate the node variables at a given point
// template <MInt nDim, class SysEqn>
// MBool DgCartesianSolver<nDim, SysEqn>::calcNodeVarsAtPoint(const MFloat* const point,
//                                                      MFloat* const nodeVars) {
//  TRACE();
//  const MInt elementId = getElementIdAtPoint(point);
//
//  if (elementId == -1) {
//    return false;
//  }
//
//  calcStateAtPoint(point, elementId, SysEqn::noNodeVars(), &m_elements.nodeVars(elementId),
//                   nodeVars);
//
//  return true;
//}


/// Calculate the node variables at a given point in an element
// template <MInt nDim, class SysEqn>
// void DgCartesianSolver<nDim, SysEqn>::calcNodeVarsAtPoint(const MFloat* const point,
//                                                   const MInt elementId,
//                                                   MFloat* const nodeVars) {
//  calcStateAtPoint(point, elementId, SysEqn::noNodeVars(), &m_elements.nodeVars(elementId),
//                   nodeVars);
//}


/// Calculate the state variables at a given point in an element
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcStateAtPoint(const MFloat* const point,
                                                       const MInt elementId,
                                                       MFloat* const state) {
  calcStateAtPoint(point, elementId, SysEqn::noVars(), &m_elements.variables(elementId), state);
}


/// \brief Calculate the state vector at a given point and for the specified sampling variable
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcSamplingVarAtPoint(const MFloat* const point, const MInt elementId,
                                                             const MInt sampleVarId, MFloat* const state,
                                                             const MBool NotUsed(interpolate)) {
  switch(sampleVarId) {
    case DG_VARS: {
      // Evaluate solution state
      calcStateAtPoint(point, elementId, state);
      break;
    }
    case DG_NODEVARS: {
      // Evaluate node variables
      calcStateAtPoint(point, elementId, SysEqn::noNodeVars(), &m_elements.nodeVars(elementId), state);
      break;
    }
    case DG_SOURCETERMS: {
      // Evaluate external source terms
      calcStateAtPoint(point, elementId, SysEqn::noVars(), &m_elements.externalSource(elementId), state);
      break;
    }
    default: {
      TERMM(1, "sampling variable not supported");
    }
  }
}


/// \brief Calculate all necessary mortar projection matrices
///
/// \author Sven Berger
/// \date   Januar 2014
///
/// Note: If the way how memory is allocated for the matrices is changed here,
///     there is a good chance that they have to be changed in mortarH() and/or
///     mortarP() as well.
/// Note: Projections have been taken from:
///     Kopriva, D. A.; Woodruff, S. L. & Hussaini, M.:
///     Computation of electromagnetic scattering with a non-conforming
///     discontinuous spectral element method,
///     Int. Journal for Numerical Methods in Engineering, 2002, 53, 105-222.
/// Note: Further information can be found in:
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::initMortarProjections() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitMortarProjection]);

  m_log << "Initializing mortar projections... ";

  using namespace dg::interpolation;

  // Init h-refinement mortar projection matrices
  // All four matrices for a particular polynomial degree are stored
  // coonsecutively. A second data structure with pointers stores the location
  // for each matrix.

  m_projectionMatricesH.clear();
  m_projectionMatrixPointersH.clear();
  // Allocate storage for all h-refinement projection matrices
  MInt sizeH = 0;
  for(MInt noNodes1D = m_minNoNodes1D; noNodes1D <= m_maxNoNodes1D; noNodes1D++) {
    sizeH += 4 * noNodes1D * noNodes1D;
  }
  m_projectionMatricesH.resize(sizeH);

  // Store pointers to matrices
  m_projectionMatrixPointersH.resize(4 * (m_maxNoNodes1D + 1));
  MInt offset = 0;
  for(MInt noNodes1D = m_minNoNodes1D; noNodes1D <= m_maxNoNodes1D; noNodes1D++) {
    for(MInt p = 0; p < 4; p++) {
      m_projectionMatrixPointersH[4 * noNodes1D + p] = &m_projectionMatricesH[offset];
      offset += noNodes1D * noNodes1D;
    }
  }

  // Calculate matrices
  if(m_sbpMode) {
    for(MInt noNodes1D = m_minNoNodes1D; noNodes1D <= m_maxNoNodes1D; noNodes1D++) {
      using namespace sbp::mortar;
      // Calculate fwd,bwd x lwr,upr projections matrices (all at onces)
      calcMortarProjectionMatrixHSBP(m_initPolyDeg, m_sbpOperator, mortarH<dg::mortar::forward>(noNodes1D, lower),
                                     mortarH<dg::mortar::forward>(noNodes1D, upper),
                                     mortarH<dg::mortar::reverse>(noNodes1D, lower),
                                     mortarH<dg::mortar::reverse>(noNodes1D, upper));
    }
  } else {
    for(MInt polyDeg = m_minPolyDeg; polyDeg <= m_maxPolyDeg; polyDeg++) {
      const MInt noNodes1D = polyDeg + 1;
      const DgInterpolation& interp = m_interpolation[polyDeg][noNodes1D];
      const MFloat* const nodes = &interp.m_nodes[0];
      const MFloat* const wBary = &interp.m_wBary[0];
      using namespace dg::mortar;

      // Calculate first forward projection matrix
      calcMortarProjectionMatrixHForward(polyDeg, nodes, wBary, dg::mortar::lower,
                                         mortarH<dg::mortar::forward>(noNodes1D, dg::mortar::lower));

      // Calculate second forward projection matrix
      calcMortarProjectionMatrixHForward(polyDeg, nodes, wBary, dg::mortar::upper,
                                         mortarH<dg::mortar::forward>(noNodes1D, dg::mortar::upper));

      // Calculate first reverse projection matrix
      calcMortarProjectionMatrixHReverse(polyDeg, nodes, wBary, dg::mortar::lower,
                                         mortarH<dg::mortar::reverse>(noNodes1D, dg::mortar::lower));

      // Calculate second reverse projection matrix
      calcMortarProjectionMatrixHReverse(polyDeg, nodes, wBary, dg::mortar::upper,
                                         mortarH<dg::mortar::reverse>(noNodes1D, dg::mortar::upper));
    }
  }

  // Init p-refinement mortar projection matrices
  // The forward/reverse matrices for a particular polyDegHi/polyDegLo
  // combination are stored consecutively. A second data structure stores
  // pointers to the beginning of each matrix.

  m_projectionMatricesP.clear();
  m_projectionMatrixPointersP.clear();
  // Allocate storage for all p-refinement projection matrices
  MInt sizeP = 0;
  for(MInt polyDegHi = m_minPolyDeg; polyDegHi <= m_maxPolyDeg; polyDegHi++) {
    for(MInt polyDegLo = m_minPolyDeg; polyDegLo < polyDegHi; polyDegLo++) {
      sizeP += 2 * (polyDegHi + 1) * (polyDegLo + 1);
    }
  }
  m_projectionMatricesP.resize(sizeP);

  // Store pointers to matrices
  m_projectionMatrixPointersP.resize(m_maxPolyDeg * (m_maxPolyDeg + 1));
  sizeP = 0;
  for(MInt polyDegHi = m_minPolyDeg; polyDegHi <= m_maxPolyDeg; polyDegHi++) {
    for(MInt polyDegLo = m_minPolyDeg; polyDegLo < polyDegHi; polyDegLo++) {
      for(MInt p = 0; p < 2; p++) {
        const MInt idx = 2 * (polyDegHi * (polyDegHi - 1) / 2 + polyDegLo) + p;
        m_projectionMatrixPointersP[idx] = &m_projectionMatricesP[sizeP];
        sizeP += (polyDegHi + 1) * (polyDegLo + 1);
      }
    }
  }

  // Calculate matrices
  if(m_sbpMode) {
    // Read SBP Projection Matrices
    for(MInt polyDegHi = m_minPolyDeg; polyDegHi <= m_maxPolyDeg; polyDegHi++) {
      for(MInt polyDegLo = m_minPolyDeg; polyDegLo < polyDegHi; polyDegLo++) {
        auto operatorLo = m_sbpOperator;
        auto operatorHi = m_sbpOperator;

        for(MUint j = 0; j < m_prefPatchesPolyDeg.size(); j++) {
          if(polyDegLo == (MInt)m_prefPatchesPolyDeg[j]) {
            operatorLo = m_prefPatchesOperators[j];
          }
          if(polyDegHi == (MInt)m_prefPatchesPolyDeg[j]) {
            operatorHi = m_prefPatchesOperators[j];
          }
        }

        sbp::mortar::calcMortarProjectionMatrixPSBP(operatorLo, operatorHi, polyDegLo, polyDegHi,
                                                    mortarP<dg::mortar::forward>(polyDegLo, polyDegHi),
                                                    mortarP<dg::mortar::reverse>(polyDegLo, polyDegHi));
      }
    }

  } else {
    // Calculate DG Projection Matrices by using Lagrange interpolation
    for(MInt polyDegHi = m_minPolyDeg; polyDegHi <= m_maxPolyDeg; polyDegHi++) {
      const MInt noNodes1DHi = polyDegHi + 1;
      const DgInterpolation& interpHi = m_interpolation[polyDegHi][noNodes1DHi];
      for(MInt polyDegLo = m_minPolyDeg; polyDegLo < polyDegHi; polyDegLo++) {
        const MInt noNodes1DLo = polyDegLo + 1;
        const DgInterpolation& interpLo = m_interpolation[polyDegLo][noNodes1DLo];

        const MFloat* const nodesHi = &interpHi.m_nodes[0];
        const MFloat* const wBaryHi = &interpHi.m_wBary[0];
        const MFloat* const nodesLo = &interpLo.m_nodes[0];
        const MFloat* const wBaryLo = &interpLo.m_wBary[0];

        // Calculate forward projection matrix from lower to higher polynomial
        // degree
        dg::mortar::calcMortarProjectionMatrixP(polyDegLo, nodesLo, wBaryLo, polyDegHi, nodesHi,
                                                mortarP<dg::mortar::forward>(polyDegLo, polyDegHi));

        // Calculate reverse projection matrix from higher to lower polynomial
        // degree
        dg::mortar::calcMortarProjectionMatrixP(polyDegHi, nodesHi, wBaryHi, polyDegLo, nodesLo,
                                                mortarP<dg::mortar::reverse>(polyDegLo, polyDegHi));
      }
    }
  }
  m_log << "done" << endl;

  RECORD_TIMER_STOP(m_timers[Timers::InitMortarProjection]);
}


/// Calculate the forward/reverse mortar projection.
///
/// \author Sven Berger
/// \date November 2014
///
/// \tparam forward True calculates the forward, false the reverse projection.
/// \param[in] srfcId Surface ID of the surface to be projected.
/// \param[in] dir Direction of the projection.
/// \param[in] source Pointer to the values to be projected.
/// \param[out] destination Pointer to the destination of the result.
/// \param[in] elem Pointer to elements.
/// \param[in] surf Pointer to surfaces.
template <MInt nDim, class SysEqn>
template <MBool forward, MInt noVars>
void DgCartesianSolver<nDim, SysEqn>::calcMortarProjection(const MInt srfcId,
                                                           const MInt dir,
                                                           MFloat* source,
                                                           MFloat* destination,
                                                           ElementCollector& elem,
                                                           SurfaceCollector& surf) {
  calcMortarProjectionH<forward, noVars>(srfcId, dir, source, destination, elem, surf);
  calcMortarProjectionP<forward, noVars>(srfcId, dir, source, destination, elem, surf);
}


/// \brief Calculate the p-refinement mortar projection.
///
/// \author Sven Berger
/// \date March 2014
///
/// \param[in] srfcId Surface ID of the surface to be projected.
/// \param[in] dir Direction of the projection.
/// \param[in] source Pointer to the values to be projected.
/// \param[out] destination Pointer to the destination of the result.
/// \param[in] elem Pointer to elements.
/// \param[in] surf Pointer to surfaces.
///
/// See also Chapter 3 in
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
template <MInt nDim, class SysEqn>
template <MBool forward, MInt noVars>
void DgCartesianSolver<nDim, SysEqn>::calcMortarProjectionP(const MInt srfcId,
                                                            const MInt dir,
                                                            MFloat* source,
                                                            MFloat* destination,
                                                            ElementCollector& elem,
                                                            SurfaceCollector& surf) {
  // TRACE();

  // Skip if this is not a p-refined surface
  const MInt side = 1 - dir % 2;
  const MInt elementId = surf.nghbrElementIds(srfcId, side);
  const MInt elementPolyDeg = elem.polyDeg(elementId);
  const MInt surfacePolyDeg = surf.polyDeg(srfcId);
  if(elementPolyDeg == surfacePolyDeg) {
    return;
  }

  // Calculate auxiliary variables for better readability
  const MInt surfaceNodes1D = surf.noNodes1D(srfcId);
  const MInt surfaceNodes1D3 = (nDim == 3) ? surfaceNodes1D : 1;
  const MInt elementNodes1D = elem.noNodes1D(elementId);
  const MInt elementNodes1D3 = (nDim == 3) ? elementNodes1D : 1;

  MFloatTensor src;
  MFloatTensor dest;

  if(forward) {
    src = MFloatTensor(source, elementNodes1D, elementNodes1D3, noVars);
    dest = MFloatTensor(destination, surfaceNodes1D, surfaceNodes1D3, noVars);
  } else {
    src = MFloatTensor(source, surfaceNodes1D, surfaceNodes1D3, noVars);
    dest = MFloatTensor(destination, surfaceNodes1D, surfaceNodes1D3, noVars);
  }

  // Reset destination matrix to 0
  dest.set(0.0);

  // Determine matrix sizes
  const MInt size0 = forward ? surfaceNodes1D : elementNodes1D;
  const MInt size1 = forward ? elementNodes1D : surfaceNodes1D;

  // Select projection matrix
  MFloatMatrix p(mortarP<forward>(elementPolyDeg, surfacePolyDeg), size0, size1);


  // Apply projection
  IF_CONSTEXPR(nDim == 2) {
    // 2D version
    for(MInt i = 0; i < size0; i++) {
      for(MInt j = 0; j < size1; j++) {
        for(MInt var = 0; var < noVars; var++) {
          dest(i, 0, var) += src(j, 0, var) * p(i, j);
        }
      }
    }
  }
  else {
    // 3D version
    for(MInt i = 0; i < size0; i++) {
      for(MInt s = 0; s < size0; s++) {
        for(MInt j = 0; j < size1; j++) {
          for(MInt k = 0; k < size1; k++) {
            for(MInt var = 0; var < noVars; var++) {
              dest(i, s, var) += src(j, k, var) * p(i, j) * p(s, k);
            }
          }
        }
      }
    }
  }
}


/// \brief Calculate the h-refinement mortar projection.
///
/// \author Sven Berger
/// \date March 2014
///
/// \param[in] srfcId Surface ID of the surface to be projected.
/// \param[in] dir Direction of the projection.
/// \param[in] source Pointer to the values to be projected.
/// \param[out] destination Pointer to the destination of the result.
/// \param[in] elem Pointer to elements.
/// \param[in] surf Pointer to surfaces.
///
/// See also Chapter 3 in
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
template <MInt nDim, class SysEqn>
template <MBool forward, MInt noVars>
void DgCartesianSolver<nDim, SysEqn>::calcMortarProjectionH(const MInt srfcId,
                                                            const MInt dir,
                                                            MFloat* source,
                                                            MFloat* destination,
                                                            ElementCollector& elem,
                                                            SurfaceCollector& surf) {
  // TRACE(); //causes a huge performance hit

  // Skip if this is not an h-refined surface
  const MInt side = 1 - dir % 2;
  const MInt elementId = surf.nghbrElementIds(srfcId, side);
  if(surf.fineCellId(srfcId) == -1 || surf.fineCellId(srfcId) == elem.cellId(elementId)) {
    return;
  }

  // Determine polynomial degree(s) and number of nodes
  const MInt elementPolyDeg = elem.polyDeg(elementId);
  const MInt surfacePolyDeg = surf.polyDeg(srfcId);
  const MInt elementNoNodes1D = elem.noNodes1D(elementId);
  const MInt surfaceNoNodes1D = surf.noNodes1D(srfcId);
  const MInt noNodes1D = forward ? elementNoNodes1D : surfaceNoNodes1D;
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;

  // Select source/destination storage
  MFloatTensor src(source, noNodes1D, noNodes1D3, noVars);
  MFloatTensor dest(destination, noNodes1D, noNodes1D3, noVars);

  // Reset destination matrix to 0
  dest.set(0.0);

  // Select correct projection matrix for h-refinement
  const MInt cellIdR = elem.cellId(elementId);
  const MInt cellIdL = surf.fineCellId(srfcId);
  const MInt orientation = surf.orientation(srfcId);

  // Select projection matrix
  const MInt dirA = (orientation == 0) ? 1 : 0;
  const MInt positionA = (grid().tree().coordinate(cellIdR, dirA) < grid().tree().coordinate(cellIdL, dirA))
                             ? dg::mortar::upper
                             : dg::mortar::lower;

  const MFloatMatrix pA(mortarH<forward>(noNodes1D, positionA), noNodes1D, noNodes1D);

  // Apply projection
  IF_CONSTEXPR(nDim == 2) {
    // 2D version
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt var = 0; var < noVars; var++) {
          dest(i, 0, var) += src(j, 0, var) * pA(i, j);
        }
      }
    }
  }
  else {
    // 3D version
    // Select additional projection matrix
    const MInt dirB = (orientation == 2) ? 1 : 2;
    const MInt positionB = (grid().tree().coordinate(cellIdL, dirB) < grid().tree().coordinate(cellIdR, dirB))
                               ? dg::mortar::lower
                               : dg::mortar::upper;

    const MFloatMatrix pB(mortarH<forward>(noNodes1D, positionB), noNodes1D, noNodes1D);

    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D; k++) {
          for(MInt l = 0; l < noNodes1D; l++) {
            for(MInt var = 0; var < noVars; var++) {
              dest(i, j, var) += src(k, l, var) * pA(i, k) * pB(j, l);
            }
          }
        }
      }
    }
  }

  // Copy results to source if p-projection is necessary
  // TODO labels:DG Move this p-refinement code out of a h-refinement method
  if(elementPolyDeg != surfacePolyDeg) {
    const MInt size = src.dim0() * src.dim1() * src.dim2();
    copy_n(&dest[0], size, &src[0]);
  }
}


/// Returns the forward/reverse projection matrix for p-refinement.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-03
///
/// \tparam forward True returns the forward, false the reverse projection
///                 matrix.
/// \param[in] sourcePolyDeg Polynomial degree of projection source.
/// \param[in] targetPolyDeg Polynomial degree of projection target.
///
/// \return Pointer to projection matrix.
///
/// Note: the source polynomial degree must always be lower than the target
///       polynomial degree.
template <MInt nDim, class SysEqn>
template <MBool forward>
MFloat* DgCartesianSolver<nDim, SysEqn>::mortarP(const MInt sourcePolyDeg, const MInt targetPolyDeg) {
  // Sanity check
  ASSERT(sourcePolyDeg < targetPolyDeg,
         "source polynomial degree (= " + to_string(sourcePolyDeg)
             + ") must be lower than target polynomial degree (= " + to_string(targetPolyDeg) + ")");

  // Calculate matrix index from polyomial degrees
  const MInt index = 2 * (targetPolyDeg * (targetPolyDeg - 1) / 2 + sourcePolyDeg) + (1 - forward);

  // Return pointer to first element in projection matrix
  return m_projectionMatrixPointersP[index];
}


/// Returns the forward/reverse projection matrix for h-refinement.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-03
///
/// \tparam forward True returns the forward, false the reverse projection
///                 matrix.
/// \param[in] polyDeg Polynomial degree of projection.
/// \param[in] position Relative position of smaller element (0 means lower
///                     element, 1 means upper element; for more info see
///                     comment on calcMortarProjectionMatrix...H{A,B})
///
/// \return Pointer to projection matrix.
template <MInt nDim, class SysEqn>
template <MBool forward>
MFloat* DgCartesianSolver<nDim, SysEqn>::mortarH(const MInt noNodes1D, const MInt position) {
  // Sanity check
  ASSERT(position == 0 || position == 1, "position must be `0` or `1` (=" + to_string(position) + ")");

  // Calculate matrix index relative from polynomial degree and position
  const MInt index = 4 * noNodes1D + 2 * (1 - forward) + position;

  // Return pointer to first element in projection matrix
  return m_projectionMatrixPointersH[index];
}


/// Exchange polynomial degrees for MPI surfaces (necessary for adaptive
/// p-refinement)
///
/// \author Sven Berger
/// \date   Februar 2015
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::exchangeMpiSurfacePolyDeg() {
  TRACE();

  m_log << "Exchanging polynomial degrees of MPI surfaces... ";

  // Create communication buffers
  vector<vector<MInt>> sendBuffer(m_noExchangeNghbrDomains);
  vector<vector<MInt>> recvBuffer(m_noExchangeNghbrDomains);
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    sendBuffer[i].resize(m_mpiSurfaces[i].size());
    recvBuffer[i].resize(m_mpiSurfaces[i].size());
  }

  // Exchange data with each neighbor exchange domain
  vector<MPI_Request> recvRequests(m_noExchangeNghbrDomains, MPI_REQUEST_NULL);
  vector<MPI_Request> sendRequests(m_noExchangeNghbrDomains, MPI_REQUEST_NULL);
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    // Copy data from surfaces to send buffer
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      sendBuffer[i][j] = m_surfaces.polyDeg(m_mpiSurfaces[i][j]);
    }

    // Begin MPI exchange
    MPI_Irecv(&recvBuffer[i][0], m_mpiSurfaces[i].size(), maia::type_traits<MInt>::mpiType(), m_exchangeNghbrDomains[i],
              m_exchangeNghbrDomains[i], mpiComm(), &recvRequests[i], AT_, "recvBuffer[i][0]");
    MPI_Isend(&sendBuffer[i][0], m_mpiSurfaces[i].size(), maia::type_traits<MInt>::mpiType(), m_exchangeNghbrDomains[i],
              domainId(), mpiComm(), &sendRequests[i], AT_, "sendBuffer[i][0]");
  }

  // Wait for MPI exchange to complete
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    MPI_Wait(&sendRequests[i], MPI_STATUS_IGNORE, AT_);
    MPI_Wait(&recvRequests[i], MPI_STATUS_IGNORE, AT_);
  }

  // Copy data from receive buffer to surfaces
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    for(vector<MInt>::size_type j = 0; j < m_mpiSurfaces[i].size(); j++) {
      const MInt srfcId = m_mpiSurfaces[i][j];
      const MInt currentPolyDeg = m_surfaces.polyDeg(srfcId);
      const MInt receivedPolyDeg = recvBuffer[i][j];

      // If new polynomial degree is lower than current polynomial degree,
      // change surface polynomial degree and recalculate node coordinates
      if(currentPolyDeg < receivedPolyDeg) {
        m_surfaces.polyDeg(srfcId) = receivedPolyDeg;
        m_surfaces.noNodes1D(srfcId) = receivedPolyDeg + 1;
        // Recompute surface node coordinates
        calcSurfaceNodeCoordinates(srfcId);
      }
    }
  }

  m_log << "done" << endl;
}


/// \brief Apply adaptive refinement (right now: only p-refinement)
///
/// \author Sven Berger
/// \date   Februar 2015
///
/// \param[in] timeStep Timestep at which adaptive refinement has been called
///
/// Note:
///     Further information can be found in:
///     Chapter 4 of "Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, Sven Berger, RWTH Aachen University, 2014".
///
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::adaptiveRefinement(const MInt timeStep) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::AdaptiveRefinement]);

  // Empty set to hold surface ids of refined elements
  set<MInt> adaptedSrfcs;

  // Vector to hold error estimate
  const MInt noElements = m_elements.size();
  vector<MFloat> errorEstimate(noElements, 0.0);

  // Update the error estimate of the elements
  calcErrorEstimate(errorEstimate);

  // Determine maximum error on local domain
  MFloat localErrorEstimate = *max_element(errorEstimate.begin(), errorEstimate.end());
  MFloat buffer = -1.0;

  // exchange local max. error estimate to obtain global max. error
  if(noDomains() > 1) {
    MPI_Allreduce(&localErrorEstimate, &buffer, 1, maia::type_traits<MFloat>::mpiType(), MPI_MAX, mpiComm(), AT_,
                  "localErrorEstimate", "buffer");
  }

  const MFloat maxErrorEstimate = (buffer > localErrorEstimate) ? buffer : localErrorEstimate;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  // Determine which elments need to be refined
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt cellId = m_elements.cellId(elementId);
    const MFloat error = errorEstimate[elementId];

    // Initialize adapted polynomial degree to current value
    MInt adaptedPolyDeg = m_elements.polyDeg(elementId);

    ////////////////////////////////////////////////////////////
    /// Adaptive refinement method 1
    /// Refine cells for x>0 +1 in each adaptive refinement step
    ////////////////////////////////////////////////////////////
    if(m_adaptivePref == DG_ADAPTIVE_TEST) {
      if(grid().tree().coordinate(cellId, 0) > 0) {
        adaptedPolyDeg += 1;
      }
    }

    ////////////////////////////////////////////////////////////
    /// Adaptive refinement method 2
    /// Sets the error relative to the maximum error
    ///
    /// Properties:
    /// The property adaptiveThreshold determines the error threshold
    /// relative to the max. error. E.g.:
    ///
    ///   *smaller value -> more elements are set to a higher polynomial degree
    ///
    ///   *larger value -> a smaller number of elements is set to each
    ///                    higher polynomial degree
    ///
    /// Note: Further information can be found in:
    /// Chapter 4 of "Implementation and validation of an adaptive
    /// hp-refinement method for the discontinuous Galerkin spectral
    /// element method. Master thesis, Sven Berger,
    /// RWTH Aachen University, 2014".
    ////////////////////////////////////////////////////////////
    else if(m_adaptivePref == DG_ADAPTIVE_GRADIENT) {
      // Determine relative error threshold
      MFloat lvlThreshold = maxErrorEstimate * m_adaptiveThreshold;

      // Iterate over all possible polynomial degrees (e.g. between minPolyDeg
      // and maxPolyDeg)
      for(MInt refLvl = m_maxPolyDeg; refLvl > m_minPolyDeg; refLvl--) {
        // If error is below current threshold, set new polynomialDegree
        if(error > lvlThreshold || refLvl == m_minPolyDeg) {
          adaptedPolyDeg = refLvl;
          break;
        }
        // Decrease threshold
        lvlThreshold *= m_adaptiveThreshold;
      }
    } else {
      TERMM(1, "Invalid adaptive refinement case.");
    }

    // Skip elements for which no refinement is necessary
    if(adaptedPolyDeg == m_elements.polyDeg(elementId)) {
      continue;
    }

    // Do not allow elements to be refined over m_maxPolyDeg
    if(adaptedPolyDeg > m_maxPolyDeg) {
      m_log << "WARNING: Element polynomial degree would exceed maximum "
               "polynomial degree."
            << endl;
      cout << "WARNING: Element polynomial degree would exceed maximum "
              "polynomial degree."
           << endl;
      adaptedPolyDeg = m_maxPolyDeg;
    }

    // Do not allow elements to be coarsed lower than m_minPolyDeg
    if(adaptedPolyDeg < m_minPolyDeg) {
      m_log << "WARNING: Element polynomial degree would be coarsed below "
               "minimum polynomial degree."
            << endl;
      cout << "WARNING: Element polynomial degree would be coarsed below "
              "minimum polynomial degree."
           << endl;
      adaptedPolyDeg = m_minPolyDeg;
    }

    // Interpolate element to degree specified by adaptedPolyDeg
    interpolateElement(elementId, adaptedPolyDeg);
  }

  const MInt begin = 0;
  const MInt end = m_surfaces.size();

  // Adapt all surfaces to the new polynomial degree
  MInt noAdaptedSrfcs = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : noAdaptedSrfcs)
#endif
  for(MInt srfcId = begin; srfcId < end; srfcId++) {
    const MInt srfcPolyDeg = m_surfaces.polyDeg(srfcId);
    const MInt internalSide = m_surfaces.internalSideId(srfcId);
    const MInt elementIdL = m_surfaces.nghbrElementIds(srfcId, (internalSide == -1) ? 0 : internalSide);
    const MInt elementIdR = m_surfaces.nghbrElementIds(srfcId, (internalSide == -1) ? 1 : internalSide);
    const MInt polyDegL = m_elements.polyDeg(elementIdL);
    const MInt polyDegR = m_elements.polyDeg(elementIdR);
    const MInt noNodes1DL = m_elements.noNodes1D(elementIdL);
    const MInt noNodes1DR = m_elements.noNodes1D(elementIdR);

    // If the maximum polynomial degree of the adjacent elements has changed,
    // update polynomial degree and re-calculate surface node coordinates
    if(max(polyDegL, polyDegR) != srfcPolyDeg) {
      m_surfaces.polyDeg(srfcId) = max(polyDegL, polyDegR);
      m_surfaces.noNodes1D(srfcId) = max(noNodes1DL, noNodes1DR);
      calcSurfaceNodeCoordinates(srfcId);
      noAdaptedSrfcs++;
    }
  }

  // Update MPI surface polynomial degree
  if(hasMpiExchange()) {
    exchangeMpiSurfacePolyDeg();
  }

  // Recalculate values that depend on element polynomial degree
  m_noTotalNodesXD = 0;
  for(MInt i = 0; i < m_elements.size(); i++) {
    const MInt noNodesXD = m_elements.noNodesXD(i);
    m_noTotalNodesXD += noNodesXD;
  }

  cout << "Grid has been refined at timestep " << timeStep << " (adapted surfaces: " << noAdaptedSrfcs << " )" << endl;
  m_log << "Grid has been refined at timestep " << timeStep << " (adapted surfaces: " << noAdaptedSrfcs << " )" << endl;

  RECORD_TIMER_STOP(m_timers[Timers::AdaptiveRefinement]);
}


/// \brief Calculate error estimate for adaptive hp-refinement.
///
/// \author Sven Berger
/// \date   Februar 2015
///
/// \param[out] errorEstimate Reference to a vector of size noElements to store
///                           the calculated error.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::calcErrorEstimate(vector<MFloat>& errorEstimate) {
  TRACE();

  if(m_adaptivePref == DG_ADAPTIVE_GRADIENT) {
    // Calculate the gradient over all surfaces to determine the maximum
    // gradient for each element
    const MInt begin = m_innerSurfacesOffset;
    const MInt end = m_surfaces.size();
    const MInt noVars = SysEqn::noVars();
    const MInt noElements = m_elements.size();

    MFloatScratchSpace error(m_surfaces.size(), FUN_, "Error estimate");
    fill_n(&error[0], m_surfaces.size(), 0.0);

// Loop over all surfaces, calculate error estimates, and update elements
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt srfcId = begin; srfcId < end; srfcId++) {
      const MInt noNodes1D = m_surfaces.noNodes1D(srfcId);
      const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
      const MFloat* stateL = &m_surfaces.variables(srfcId, 0);
      const MFloat* stateR = &m_surfaces.variables(srfcId, 1);


      const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noNodes1D3, noVars);
      const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noNodes1D3, noVars);

      // Calculate error estimate as difference between left and right state
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D3; j++) {
          for(MInt var = 0; var < noVars; var++) {
            error[srfcId] += fabs(uL(i, j, var) - uR(i, j, var));
          }
        }
      }

      // Normalize by number of nodes on surface
      error[srfcId] = error[srfcId] / (noNodes1D * noNodes1D3);
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // sum element error contributions
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      const MInt noSrfcs = 2 * nDim;
      for(MInt i = 0; i < noSrfcs; i++) {
        const MInt srfcId = m_elements.surfaceIds(elementId, i);
        errorEstimate[elementId] += error[srfcId];
      }
    }
  }
}


/// \brief Interpolate an element to a different polynomial degree
///
/// \author Sven Berger
/// \date   Februar 2015
///
///
/// \param[in] elementId ElementId of the element to be interpolated to a new
///                      polynomial degree.
/// \param[in] adaptedPolyDeg New polynomial degree of the element.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::interpolateElement(const MInt elementId, const MInt adaptedPolyDeg) {
  TRACE();

  // Interpolate element data to new integration nodes at new polynomial degree
  // in temporary storage, then copy result to destination location
  // (do not use Scratch space to remain thread-safe!)
  const MInt adaptedNodes1D = adaptedPolyDeg + 1;
  const MInt adaptedNodes1D3 = (nDim == 3) ? adaptedNodes1D : 1;
  const MInt noVars = SysEqn::noVars();
  const MInt elemPolyDeg = m_elements.polyDeg(elementId);
  const MInt elemNodes1D = m_elements.noNodes1D(elementId);
  const MInt size = adaptedNodes1D * adaptedNodes1D * adaptedNodes1D3 * noVars;
  vector<MFloat> buffer(size);
  using namespace dg::interpolation;

  // Variables
  if(elemPolyDeg < adaptedPolyDeg) {
    interpolateNodes<nDim>(&m_elements.variables(elementId), mortarP<dg::mortar::forward>(elemPolyDeg, adaptedPolyDeg),
                           elemNodes1D, adaptedNodes1D, noVars, &buffer[0]);
    copy(buffer.begin(), buffer.end(), &m_elements.variables(elementId));

    // Time Integration Storage
    interpolateNodes<nDim>(&m_elements.timeIntStorage(elementId),
                           mortarP<dg::mortar::forward>(elemPolyDeg, adaptedPolyDeg), elemNodes1D, adaptedNodes1D,
                           noVars, &buffer[0]);
    copy(buffer.begin(), buffer.end(), &m_elements.timeIntStorage(elementId));
  } else {
    interpolateNodes<nDim>(&m_elements.variables(elementId), mortarP<dg::mortar::reverse>(adaptedPolyDeg, elemPolyDeg),
                           elemNodes1D, adaptedNodes1D, noVars, &buffer[0]);
    copy(buffer.begin(), buffer.end(), &m_elements.variables(elementId));

    // Time Integration Storage
    interpolateNodes<nDim>(&m_elements.timeIntStorage(elementId),
                           mortarP<dg::mortar::reverse>(adaptedPolyDeg, elemPolyDeg), elemNodes1D, adaptedNodes1D,
                           noVars, &buffer[0]);
    copy(buffer.begin(), buffer.end(), &m_elements.timeIntStorage(elementId));
  }

  // Update polynomial degree
  m_elements.polyDeg(elementId) = adaptedPolyDeg;
  m_elements.noNodes1D(elementId) = adaptedNodes1D;

  // Recalculate node coordinates of the element
  calcElementNodeCoordinates(elementId);
}


/// Check if an element has a neighbor cell in the given direction.
///
/// \author Sven Berger
/// \date   Februar 2015
///
/// \param[in] elementId Element for which neighbor level is requested.
/// \param[in] dir Direction of neighbor element (-x,+x,-y,... = 0,1,2,...).
///
/// \return If an cell exists in given direction.
///
/// This method is smart and does not just return false if there is no
/// same-level neighbor, but also checks if there are refined (smaller) or
/// coarsened (larger) neighbor element(s).
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::hasNeighborCell(const MInt elementId, const MInt dir) {
  // TRACE();

  const MInt cellId = m_elements.cellId(elementId);
  MInt nghbrCellId = grid().tree().neighbor(cellId, dir);

  // If neigbor cell does not exist on current level, check parent level
  if(nghbrCellId < 0 && grid().tree().hasParent(cellId)) {
    nghbrCellId = grid().tree().neighbor(grid().tree().parent(cellId), dir);
  }

  return nghbrCellId >= 0;
}


/// Get h-refinement (cell) level of element.
///
/// \author Sven Berger
/// \date   Februar 2015
///
/// \param[in] elementId Id of the element for which level is requested.
/// \return The level of the corresponding cell.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::getLevelByElementId(const MInt elementId) {
  // TRACE();
  ASSERT(elementId >= 0, "Invalid elementId");

  const MInt cellId = m_elements.cellId(elementId);
  const MInt level = grid().tree().level(cellId);

  return level;
}


/// Get h-refinement (cell) level of neighbor element.
///
/// \author Sven Berger
/// \date   Februar 2015
///
/// \param[in] elementId Element for which neighbor level is requested.
/// \param[in] dir Direction of neighbor element (-x,+x,-y,... = 0,1,2,...).
///
/// \return Level of the neighbor element or -1 if neighbor does not exist.
// template <MInt nDim, class SysEqn>
// MInt DgCartesianSolver<nDim, SysEqn>::getLevelOfNeighborElement(const MInt elementId,
// const MInt dir) {
//  // TRACE();
//  ASSERT(elementId >= 0, "Invalid elementId");
//
//  const MInt cellId = m_elements.cellId(elementId);
//
//  return getLevelOfDirectNeighborCell(cellId, dir);
//}


/// Get h-refinement (cell) level of neighbor cell.
///
/// \author Sven Berger, Rodrigo Miguez
/// \date   Februar 2015, October 2019
///
/// \param[in] cellId for which neighbor level is requested.
/// \param[in] dir Direction of neighbor cell (-x,+x,-y,... = 0,1,2,...).
///
/// \return Level of the neighbor element or -1 if neighbor does not exist.
///
/// This method is smart and does not just return -1 if there is no same-level
/// neighbor, but also checks if there are refined (smaller) or coarsened
/// (larger) neighbor element(s).
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::getLevelOfDirectNeighborCell(const MInt cellId, const MInt dir) {
  // TRACE();

  const MInt nghbrCellId = grid().tree().neighbor(cellId, dir);
  const MInt cellLvl = grid().tree().level(cellId);

  // Check for finer nghbr. The finer nghbr must be directly adjacent to current cell.
  if(nghbrCellId > 0 && grid().tree().hasChildren(nghbrCellId)) {
    TERMM(1, "The following is not yet tested! If it works, delete this TERMM!");
    for(MInt child = 0; child < IPOW2(nDim); child++) {
      if(!(childCodePro[dir] & (1 << child))) continue;
      if(grid().tree().child(nghbrCellId, child) > -1) return cellLvl + 1;
    }
  }

  // If neigbor cell does not exist on current level, check parent level
  if(nghbrCellId < 0 && grid().tree().hasParent(cellId)) {
    if(grid().tree().neighbor(grid().tree().parent(cellId), dir)) {
      return cellLvl - 1;
    }
  }

  return nghbrCellId < 0 ? -1 : cellLvl;
}


/// \brief Return true if adaptive hp-refinement is activated.
///
/// \author Sven Berger
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::hasAdaptivePref() const {
  TRACE();

  return (m_adaptivePref > 0 && hasPref());
}

/// \brief Return true if p-refinement is set
///
/// \author Bjoern Peeters
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::hasPref() const {
  TRACE();

  return (m_pref == 1 && m_minPolyDeg != m_maxPolyDeg);
}


/// \brief Return if the given timestep is an adaptation timestep.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-11-10
///
/// \param[in] timeStep Time step to check if it is an adaptation time step.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::isAdaptationTimeStep(const MInt timeStep) const {
  TRACE();

  // Check if adaptive p-refinement is enabled: The adaptation is done once
  // after the first time step (i.e. before the second time step), and then
  // before each adaptation interval
  return (hasAdaptivePref() && (timeStep == 2 || timeStep % m_adaptiveInterval == 0));
}


/// \brief Return true if a surface is a MPI surface
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-01
///
/// \param[in] id Surface id to check.
template <MInt nDim, class SysEqn>
MBool DgCartesianSolver<nDim, SysEqn>::isMpiSurface(const MInt id) const {
  TRACE();

  const MInt begin = m_mpiSurfacesOffset;
  const MInt end = m_mpiSurfacesOffset + m_noMpiSurfaces;
  MBool isMpiSrfc = false;

  if(begin <= id && id < end) {
    isMpiSrfc = true;
  }

  return isMpiSrfc;
}


/// \brief Return h-element id of a given element (if it exists).
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-03-06
///
/// \param[in] elementId Element id to check.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::getHElementId(const MInt elementId) const {
  TRACE();

  MInt hElementId = -1;

  const MInt noHElements = m_helements.size();
  for(MInt hId = 0; hId < noHElements; hId++) {
    if(m_helements.elementId(hId) == elementId) {
      hElementId = hId;
      break;
    }
  }

  return hElementId;
}


/// \brief Reset solver (for load balancing)
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::resetSolver() {
  TRACE();

  finalizeMpiExchange();

  // Set to -1 such that the new size is calculated in allocateAndInitSolverMemory
  m_maxNoSurfaces = -1;

  // Cleanup communication memory
  for(MInt i = 0; i < m_noExchangeNghbrDomains; i++) {
    m_mpiSurfaces[i].clear();
    m_sendBuffers[i].clear();
    m_recvBuffers[i].clear();
  }
  m_mpiSurfaces.clear();
  m_sendBuffers.clear();
  m_recvBuffers.clear();
  m_exchangeNghbrDomains.clear();
  m_sendRequests.clear();
  m_recvRequests.clear();
  m_noExchangeNghbrDomains = 0;

  m_isInitSolver = false;
}


/// \brief Return the number of DG load types.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::noLoadTypes() const {
  TRACE();

  MInt totalNoDgLT = 0;
  // Add a load type for each polynomial degree (even if there might be no element with this
  // polyomial degree)
  totalNoDgLT += m_maxPolyDeg - m_minPolyDeg + 1;

  // Add load type for boundary condition elements (assumes only one polynomial degree)
  if(m_weightDgRbcElements) {
    totalNoDgLT += 1;
  }

  return totalNoDgLT;
}


/// \brief Return the default weights for all load quantities
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getDefaultWeights(MFloat* weights, std::vector<MString>& names) const {
  TRACE();

  const MInt noPolyDegs = m_maxPolyDeg - m_minPolyDeg + 1;
  for(MInt i = 0; i < noPolyDegs; i++) {
    weights[i] = 0.2;
    names[i] = "dg_node_p" + std::to_string(m_minPolyDeg + i);
  }

  if(m_weightDgRbcElements) {
    weights[noPolyDegs] = 0.5;
    names[noPolyDegs] = "dg_node_rbc";
  }
}


/// \brief Return the cumulative load quantities on this domain.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[out] loadQuantities Storage for load quantities.
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getLoadQuantities(MInt* const loadQuantities) const {
  TRACE();

  // Reset
  std::fill_n(&loadQuantities[0], noLoadTypes(), 0);

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  const MInt noPolyDegs = m_maxPolyDeg - m_minPolyDeg + 1;

  if(noPolyDegs > 1) {
    // Total DOFs of different polynomial degrees
    for(MInt i = 0; i < noPolyDegs; i++) {
      loadQuantities[i] = m_statLocalNoActiveDOFsPolyDeg[i];
    }
  } else {
    // Only one polynomial degree
    loadQuantities[0] = m_statLocalNoActiveDOFs;
  }

  // Add number of boundary condition element nodes
  if(m_weightDgRbcElements) {
    loadQuantities[noPolyDegs] = 0;
    for(auto&& bc : m_boundaryConditions) {
      loadQuantities[noPolyDegs] += bc->getLocalNoNodes();
    }
  }
}


/// \brief Return the load of a single cell (given computational weights).
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] cellId Requested grid cell id.
/// \param[in] weights Computational weights for different simulation components.
/// \return Cell load.
template <MInt nDim, class SysEqn>
MFloat DgCartesianSolver<nDim, SysEqn>::getCellLoad(const MInt gridCellId, const MFloat* const weights) const {
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

  const MInt elementId = m_elements.getElementByCellId(cellId);

  // Default cell load
  MFloat cellLoad = 0.0;

  // Add load if there is an element
  if(elementId != -1) {
    const MInt polyDeg = m_elements.polyDeg(elementId);
    const MInt noNodesXD = m_elements.noNodesXD(elementId);
    // Weight the number of nodes and add a constant weight for each element
    // Note: m_weightPerElement should only be != 0 when getCellLoad is called from setCellWeights!
    cellLoad = m_weightPerElement + noNodesXD * weights[polyDeg - m_minPolyDeg];

    // Add weight for boundary condition element
    if(m_weightDgRbcElements) {
      for(auto&& bc : m_boundaryConditions) {
        if(bc->hasBcElement(elementId)) {
          cellLoad += noNodesXD * weights[m_maxPolyDeg - m_minPolyDeg + 1];
        }
      }
    }
  }

  return cellLoad;
}


/// Set cell weights with constant weighting factor weightPerNode
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  ASSERT(isActive(), "solver is not active");

  const MInt noCells = grid().noInternalCells();
  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt noWeights = noLoadTypes();
  const MInt offset = solverId() * noCellsGrid;

  std::vector<MFloat> weights(noWeights);
  std::fill(weights.begin(), weights.end(), m_weightPerNode);

  for(MInt cellId = 0; cellId < noCells; cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    solverCellWeight[offset + gridCellId] = getCellLoad(gridCellId, &weights[0]);
  }
}


/// \brief Change local into global ids.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::localToGlobalIds() {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  const MInt domainOffset = grid().domainOffset(domainId());
  const MInt noElements = m_elements.size();
  // Store the global cell id in all elements
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt globalCellId = m_elements.cellId(elementId) + domainOffset;
    m_elements.cellId(elementId) = globalCellId;
  }
}


/// \brief Change global into local ids.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::globalToLocalIds() {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  const MInt domainOffset = grid().domainOffset(domainId());
  const MInt noElements = m_elements.size();
  // Change global cell ids to local cell ids
  for(MInt elementId = 0; elementId < noElements; elementId++) {
    const MInt localCellId = m_elements.cellId(elementId) - domainOffset;
    m_elements.cellId(elementId) = localCellId;
  }
}


/// \brief Get data type of cell data.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] dataId Data id of requested cell data.
/// \return Data type of cell data.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::cellDataTypeDlb(const MInt dataId) const {
  TRACE();

  if(!isActive()) {
    TERMM(1, "Error: cellDataTypeDlb() might give wrong results on inactive ranks.");
    return -1;
  }

  MInt dataType = -1;
  if(dataId > -1 && dataId < noDgCartesianSolverCellData()) {
    // DG solver cell data
    dataType = s_cellDataTypeDlb[dataId];
  } else if(dataId >= noDgCartesianSolverCellData() && dataId < noCellDataDlb()) {
    // Boundary condition cell data
    MInt offset = noDgCartesianSolverCellData();
    for(auto&& bc : m_boundaryConditions) {
      const MInt bcNoCellData = bc->noCellDataDlb();
      if(dataId >= offset && dataId < offset + bcNoCellData) {
        dataType = bc->cellDataTypeDlb(dataId - offset);
        break;
      }
      offset += bcNoCellData;
    }
  } else {
    TERMM(1, "The requested dataId is not valid: " + to_string(dataId) + " (" + to_string(noDgCartesianSolverCellData())
                 + ", " + to_string(noCellDataDlb()) + ")");
  }

  return dataType;
}


/// \brief Return data size of cell data.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] dataId Data id of requested cell data.
/// \param[in] gridCellId Requested grid cell id.
/// \return Data size of requested cell data for given cell.
template <MInt nDim, class SysEqn>
MInt DgCartesianSolver<nDim, SysEqn>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  TRACE();

  // Inactive ranks do not have any data to communicate
  if(!isActive()) {
    return 0;
  }

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0) {
    return 0;
  }

  MInt dataSize = 0;

  if(dataId > -1 && dataId < noDgCartesianSolverCellData()) {
    // DG solver cell data
    const MInt elementId = m_elements.getElementByCellId(cellId);

    if(elementId != -1) {
      const MInt noNodesXD = m_elements.noNodesXD(elementId);
      switch(dataId) {
        case CellData::ELEM_CELL_ID:
        case CellData::ELEM_POLY_DEG:
          dataSize = 1;
          break;
        case CellData::ELEM_NO_NODES_1D:
          dataSize = 1;
          break;
        case CellData::ELEM_VARIABLES:
          dataSize = noNodesXD * SysEqn::noVars();
          break;
        case CellData::ELEM_NODE_VARS:
          dataSize = noNodesXD * SysEqn::noNodeVars();
          break;
        default:
          TERMM(1, "Unknown data id. (" + to_string(dataId) + ")");
          break;
      }
    }
  } else if(dataId >= noDgCartesianSolverCellData() && dataId < noCellDataDlb()) {
    // Boundary condition cell data
    MInt offset = noDgCartesianSolverCellData();
    for(auto&& bc : m_boundaryConditions) {
      const MInt bcNoCellData = bc->noCellDataDlb();
      if(dataId >= offset && dataId < offset + bcNoCellData) {
        dataSize = bc->cellDataSizeDlb(dataId - offset, cellId);
        break;
      }
      offset += bcNoCellData;
    }
  } else {
    TERMM(1, "The requested dataId is not valid.");
  }

  return dataSize;
}


/// \brief Reinitialize solver prior to setting solution data in DLB.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::balancePre() {
  TRACE();

  // Set reinitialization stage
  m_loadBalancingReinitStage = 0;

  // Store currently used memory
  const MLong previouslyAllocated = allocatedBytes();

  // Set to prevent initialization (e.g. due to coupling) that overwrites redistributed data
  m_restart = true;

  // Update the grid proxy for this solver
  grid().update();

  // Just reset parallelization information if solver is not active and cleanup containers
  if(!isActive()) {
    updateDomainInfo(-1, -1, MPI_COMM_NULL, AT_);
    m_elements.reset(0);
    m_helements.reset(0);
    m_surfaces.reset(0);
    return;
  }

  // Set new domain info for solver
  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);

  allocateAndInitSolverMemory();

  // Print information on used memory
  printAllocatedMemory(previouslyAllocated, "DgCartesianSolver (solverId = " + to_string(m_solverId) + ")", mpiComm());

  RECORD_TIMER_START(m_timers[Timers::RunInit]);
  RECORD_TIMER_START(m_timers[Timers::InitSolverObjects]);
  // Initialization from initSolver()
  initGridMap();
  initElements();
  initHElements();
  // Note: Polynomial degree is communicated and set from gridcontroller
  RECORD_TIMER_STOP(m_timers[Timers::InitSolverObjects]);
  RECORD_TIMER_STOP(m_timers[Timers::RunInit]);
}


/// \brief Reinitialize solver after setting solution data.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::balancePost() {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set reinitialization stage
  m_loadBalancingReinitStage = 1;

  RECORD_TIMER_START(m_timers[Timers::RunInit]);
  RECORD_TIMER_START(m_timers[Timers::InitSolverObjects]);

  initInterpolation();

  // Set new internal data size
  m_internalDataSize = m_elements.size() * pow(m_maxNoNodes1D, nDim) * SysEqn::noVars();

  // Calculate new total number of nodes
  m_noTotalNodesXD = 0;
  for(MInt i = 0; i < m_elements.size(); i++) {
    const MInt noNodesXD = m_elements.noNodesXD(i);
    m_noTotalNodesXD += noNodesXD;
  }

  initNodeCoordinates();

  initJacobian();

  checkCellProperties();

  // Prevent the RBC from overwriting its solution
  m_restart = true;

  initSurfaces();

  if(noDomains() > 1) {
    initMpiExchange();
  }

  if(useSponge()) {
    m_log << "Reinitializing sponge... ";
    sponge().init(m_maxPolyDeg, &grid(), &m_elements, &m_surfaces, &m_boundaryConditions, &m_sysEqn, mpiComm());
    m_log << "done" << endl;
  }

  if(noDomains() > 1 && hasPref()) {
    exchangeMpiSurfacePolyDeg();
  }

  initMortarProjections();

  // @ansgar TODO labels:DG,PP wont work, no cleanup, save data prior to DLB? also the point ordering in the
  // output file might change
  m_pointData.init();
  m_surfaceData.init();
  m_volumeData.init();

  m_slice.init();

  if(m_pointData.enabled() || m_surfaceData.enabled()) {
    TERMM(1, "DLB for point/surface/volumeData not supported yet");
  }

  initSimulationStatistics();

  // From initData()
  // Reset current Runge Kutta stage
  m_rkStage = 0;

  // Update all node variables, if they are used
  if(SysEqn::noNodeVars() > 0) {
    updateNodeVariables();
    extendNodeVariables(); // @ansgar TODO labels:DG,totest,toremove check if needed!
  }

  m_isInitSolver = true;
  m_isInitData = true;
  // TODO labels:DG endAutoSaveTime
  m_isInitMainLoop = true;

  // outputInitSummary();
  RECORD_TIMER_STOP(m_timers[Timers::InitSolverObjects]);
  RECORD_TIMER_STOP(m_timers[Timers::RunInit]);

  // Set reinitialization stage
  m_loadBalancingReinitStage = 2;
}


/// Get global solver variables that need to be communicated for DLB (same on each domain)
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getGlobalSolverVars(std::vector<MFloat>& globalFloatVars,
                                                          std::vector<MInt>& globalIntVars) {
  TRACE();

  globalFloatVars.push_back(m_time);
  globalFloatVars.push_back(m_dt);

  globalIntVars.push_back(m_timeStep);
  globalIntVars.push_back(m_firstTimeStep);
  globalIntVars.push_back(m_noAnalyzeTimeSteps);
}


/// Set global solver variables for DLB (same on each domain)
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::setGlobalSolverVars(std::vector<MFloat>& globalFloatVars,
                                                          std::vector<MInt>& globalIntVars) {
  TRACE();

  m_time = globalFloatVars[0];
  m_dt = globalFloatVars[1];

  m_timeStep = globalIntVars[0];
  m_firstTimeStep = globalIntVars[1];
  m_noAnalyzeTimeSteps = globalIntVars[2];
}


/// Get solver timings
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings,
                                                       const MBool allTimings) {
  TRACE();
  const MString namePrefix = "b" + std::to_string(solverId()) + "_";

  const MFloat load = returnLoadRecord();
  const MFloat idle = returnIdleRecord();

  solverTimings.emplace_back(namePrefix + "loadDgCartesianSolver", load);
  solverTimings.emplace_back(namePrefix + "idleDgCartesianSolver", idle);

#ifdef MAIA_TIMER_FUNCTION
  solverTimings.emplace_back(namePrefix + "timeStepRk", RETURN_TIMER_TIME(m_timers[Timers::RungeKuttaStep]));

  if(allTimings) {
    // Full set of timings
    solverTimings.emplace_back(namePrefix + "calcDgTimeDerivative", RETURN_TIMER_TIME(m_timers[Timers::TimeDeriv]));

    solverTimings.emplace_back(namePrefix + "resetRHS", RETURN_TIMER_TIME(m_timers[Timers::ResetRHS]));
    solverTimings.emplace_back(namePrefix + "prolong_to_surfaces", RETURN_TIMER_TIME(m_timers[Timers::Prolong]));
    solverTimings.emplace_back(namePrefix + "forward_projection",
                               RETURN_TIMER_TIME(m_timers[Timers::ForwardProjection]));

    solverTimings.emplace_back(namePrefix + "MPI_surface_exchange", RETURN_TIMER_TIME(m_timers[Timers::SurfExchange]));
    solverTimings.emplace_back(namePrefix + "communication", RETURN_TIMER_TIME(m_timers[Timers::SurfExchangeComm]));
    solverTimings.emplace_back(namePrefix + "copy_operations", RETURN_TIMER_TIME(m_timers[Timers::SurfExchangeCopy]));
    solverTimings.emplace_back(namePrefix + "waiting", RETURN_TIMER_TIME(m_timers[Timers::SurfExchangeWait]));
    solverTimings.emplace_back(namePrefix + "waiting_send", RETURN_TIMER_TIME(m_timers[Timers::SEWaitSend]));
    solverTimings.emplace_back(namePrefix + "waiting_recv", RETURN_TIMER_TIME(m_timers[Timers::SEWaitRecv]));

    solverTimings.emplace_back(namePrefix + "calcVolumeIntegral", RETURN_TIMER_TIME(m_timers[Timers::VolInt]));

    solverTimings.emplace_back(namePrefix + "flux_calculation", RETURN_TIMER_TIME(m_timers[Timers::Flux]));
    solverTimings.emplace_back(namePrefix + "calcBoundarySurfaceFlux", RETURN_TIMER_TIME(m_timers[Timers::FluxBndry]));
    solverTimings.emplace_back(namePrefix + "calcInnerSurfaceFlux", RETURN_TIMER_TIME(m_timers[Timers::FluxInner]));
    solverTimings.emplace_back(namePrefix + "calcMpiSurfaceFlux", RETURN_TIMER_TIME(m_timers[Timers::FluxMPI]));

    solverTimings.emplace_back(namePrefix + "surface_integrals", RETURN_TIMER_TIME(m_timers[Timers::SurfInt]));
    solverTimings.emplace_back(namePrefix + "applyJacobian", RETURN_TIMER_TIME(m_timers[Timers::Jacobian]));
    solverTimings.emplace_back(namePrefix + "calcSourceTerms", RETURN_TIMER_TIME(m_timers[Timers::Sources]));
    solverTimings.emplace_back(namePrefix + "resetExtSources",
                               RETURN_TIMER_TIME(m_timers[Timers::ResetExternalSources]));
    solverTimings.emplace_back(namePrefix + "applyExternalSources",
                               RETURN_TIMER_TIME(m_timers[Timers::ExternalSources]));
    solverTimings.emplace_back(namePrefix + "calcSpongeTerms", RETURN_TIMER_TIME(m_timers[Timers::Sponge]));

    solverTimings.emplace_back(namePrefix + "time_integration", RETURN_TIMER_TIME(m_timers[Timers::TimeInt]));

    solverTimings.emplace_back(namePrefix + "IO", RETURN_TIMER_TIME(m_timers[Timers::MainLoopIO]));
    solverTimings.emplace_back(namePrefix + "solution_analysis", RETURN_TIMER_TIME(m_timers[Timers::Analysis]));
    solverTimings.emplace_back(namePrefix + "adaptive_refinement",
                               RETURN_TIMER_TIME(m_timers[Timers::AdaptiveRefinement]));
  } else {
    // Reduced/essential set of timings
    // MPI exchange
    solverTimings.emplace_back(namePrefix + "MPI_surface_exchange", RETURN_TIMER_TIME(m_timers[Timers::SurfExchange]));
    // Boundary conditions
    solverTimings.emplace_back(namePrefix + "calcBoundarySurfaceFlux", RETURN_TIMER_TIME(m_timers[Timers::FluxBndry]));
  }
#endif
}


/// Return decomposition information, i.e. number of local elements,...
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getDomainDecompositionInformation(
    std::vector<std::pair<MString, MInt>>& domainInfo) {
  TRACE();

  const MString namePrefix = "b" + std::to_string(solverId()) + "_";

  // Number of DG elements
  const MInt noElements = m_elements.size();
  domainInfo.emplace_back(namePrefix + "noDgElements", noElements);

  // Number of additional boundary condition elements
  MInt noBcElements = 0;
  for(auto&& bc : m_boundaryConditions) {
    for(MInt elementId = 0; elementId < noElements; elementId++) {
      if(bc->hasBcElement(elementId)) {
        noBcElements++;
      }
    }
  }
  domainInfo.emplace_back(namePrefix + "noDgBcElements", noBcElements);
}


/// \brief Return sampling properties for the DG solver
template <MInt nDim, class SysEqn>
void DgCartesianSolver<nDim, SysEqn>::getSolverSamplingProperties(std::vector<MInt>& samplingVarIds,
                                                                  std::vector<MInt>& noSamplingVars,
                                                                  std::vector<std::vector<MString>>& samplingVarNames,
                                                                  const MString featureName) {
  TRACE();

  // Read sampling variable names (for a specific feature)
  std::vector<MString> varNamesList;
  MInt noSampleVars = readSolverSamplingVarNames(varNamesList, featureName);

  // Set default sampling variables if none specified
  if(noSampleVars == 0) {
    varNamesList.push_back("DG_VARS");
    noSampleVars = 1;
  }

  for(MInt i = 0; i < noSampleVars; i++) {
    const MInt samplingVar = string2enum(varNamesList[i]);
    std::vector<MString> varNames;

    auto samplingVarIt = std::find(samplingVarIds.begin(), samplingVarIds.end(), samplingVar);
    if(samplingVarIt != samplingVarIds.end()) {
      TERMM(1, "Sampling variable '" + varNamesList[i] + "' already specified.");
    }

    switch(samplingVar) {
      case DG_VARS: {
        const MInt noVars = SysEqn::noVars();

        samplingVarIds.push_back(DG_VARS);
        noSamplingVars.push_back(noVars);

        varNames.resize(noVars);
        for(MInt varId = 0; varId < noVars; varId++) {
          // TODO labels:DG needs to be fixed if conservative and primitive variables are different
          varNames[varId] = SysEqn::consVarNames(varId);
        }

        samplingVarNames.push_back(varNames);
        break;
      }
      case DG_NODEVARS: {
        const MInt noVars = SysEqn::noNodeVars();

        samplingVarIds.push_back(DG_NODEVARS);
        noSamplingVars.push_back(noVars);

        varNames.resize(noVars);
        for(MInt varId = 0; varId < noVars; varId++) {
          varNames[varId] = SysEqn::nodeVarNames(varId);
        }
        samplingVarNames.push_back(varNames);
        break;
      }
      case DG_SOURCETERMS: {
        const MInt noVars = SysEqn::noVars();

        samplingVarIds.push_back(DG_SOURCETERMS);
        noSamplingVars.push_back(noVars);

        varNames.resize(noVars);
        for(MInt varId = 0; varId < noVars; varId++) {
          varNames[varId] = "source_" + SysEqn::consVarNames(varId);
        }
        samplingVarNames.push_back(varNames);
        break;
      }
      default: {
        TERMM(1, "Unknown sampling variable: " + varNamesList[i]);
        break;
      }
    }
  }
}


// Enable/disable compiler-specific diagnostics
#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic pop
#endif
