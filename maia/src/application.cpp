// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "application.h"

#include <chrono>

#if not defined(MAIA_DISABLE_DG)
#include "COUPLER/dgccacousticperturbsourcefiles.h"
#include "DG/dgcartesiansolver.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "DG/dgcartesiansyseqnlinearscalaradv.h"
#endif

#include "ACA/acasolver.h"
#include "COMM/mpioverride.h"
#include "COUPLER/couplerfvmultilevel.h"
#include "COUPLER/couplerlbfv.h"
#include "COUPLER/couplerlbfveemultiphase.h"
#include "COUPLER/couplerlblb.h"
#include "COUPLER/fvcartesianinterpolation.h"
#include "COUPLER/fvmbzonal.h"
#include "COUPLER/fvparticle.h"
#include "COUPLER/fvzonalrtv.h"
#include "COUPLER/fvzonalstg.h"
#include "COUPLER/lbdgape.h"
#include "COUPLER/lblpt.h"
#include "COUPLER/lbrb.h"
#include "COUPLER/lsfv.h"
#include "COUPLER/lsfvcombustion.h"
#include "COUPLER/lsfvmb.h"
#include "COUPLER/lslb.h"
#include "COUPLER/lslbsurface.h"
#include "FC/fcsolver.h"
#include "FV/fvcartesianapesolver2d.h"
#include "FV/fvcartesianapesolver3d.h"
#include "FV/fvcartesiansolverxd.h"
#include "FV/fvcartesiansyseqndetchem.h"
#include "FV/fvcartesiansyseqneegas.h"
#include "FV/fvcartesiansyseqnns.h"
#include "FV/fvcartesiansyseqnrans.h"
#ifndef MAIA_DISABLE_STRUCTURED
#include "FV/fvstructuredsolver2d.h"
#include "FV/fvstructuredsolver3d.h"
#endif
#include "GRID/cartesiangridcontroller.h"
#include "GRID/cartesiangridgenpar.h"
#include "IO/parallelio.h"
#include "LB/lbsolverdxqy.h"
#include "LB/lbsolverfactory.h"
#include "LPT/lpt.h"
#include "LS/lscartesiansolver.h"
#include "LS/lscartesiansolverfactory.h"
#include "POST/postdata.h"
#include "POST/postprocessing.h"
#include "POST/postprocessingcontroller.h"
#include "POST/postprocessingdg.h"
#include "POST/postprocessingfv.h"
#include "POST/postprocessingfvlpt.h"
#include "POST/postprocessinglb.h"
#include "POST/postprocessinglblpt.h"
#include "POST/postprocessinglpt.h"
#include "RB/rigidbodies.h"
#include "executionrecipe.h"

using namespace std;
using namespace maia;

// TODO labels:toenhance move these to solvertraits.h to allow usage in other places (e.g. in the coupling
// conditions)
namespace {

// Auxiliary type traits to support selecting the correct FV-APE solver class
template <MInt nDim>
struct FvApeSolverXD {};
template <>
struct FvApeSolverXD<2> {
  using type = FvApeSolver2D;
};
template <>
struct FvApeSolverXD<3> {
  using type = FvApeSolver3D;
};

// Auxiliary type traits to support selecting the correct structured solver class
#ifndef MAIA_DISABLE_STRUCTURED
template <MInt nDim>
struct FvStructuredSolverXD {};
template <>
struct FvStructuredSolverXD<2> {
  using type = FvStructuredSolver2D;
};
template <>
struct FvStructuredSolverXD<3> {
  using type = FvStructuredSolver3D;
};
#endif

// Auxiliary type traits to support selecting the correct LS solver class
template <MInt nDim>
struct LsCartesianSolverXD {};
template <>
struct LsCartesianSolverXD<2> {
  using type = LsCartesianSolver<2>;
};
template <>
struct LsCartesianSolverXD<3> {
  using type = LsCartesianSolver<3>;
};

template <MInt nDim, class SysEqn>
struct LsFvCombustionXD {};
template <>
struct LsFvCombustionXD<2, FvSysEqnNS<2>> {
  using type = LsFvCombustion<2, FvSysEqnNS<2>>;
};
template <>
struct LsFvCombustionXD<3, FvSysEqnNS<3>> {
  using type = LsFvCombustion<3, FvSysEqnNS<3>>;
};
template <>
struct LsFvCombustionXD<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>> {
  using type = LsFvCombustion<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
};
template <>
struct LsFvCombustionXD<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>> {
  using type = LsFvCombustion<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
};

template <MInt nDim, class SysEqn>
struct LsFvMbXD {};
template <>
struct LsFvMbXD<2, FvSysEqnNS<2>> {
  using type = LsFvMb<2, FvSysEqnNS<2>>;
};
template <>
struct LsFvMbXD<3, FvSysEqnNS<3>> {
  using type = LsFvMb<3, FvSysEqnNS<3>>;
};
template <>
struct LsFvMbXD<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>> {
  using type = LsFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
};
template <>
struct LsFvMbXD<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>> {
  using type = LsFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
};

template <MInt nDim, class SysEqn>
struct CouplingLsFvXD {};
template <>
struct CouplingLsFvXD<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>> {
  using type = CouplingLsFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
};
template <>
struct CouplingLsFvXD<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>> {
  using type = CouplingLsFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
};
template <>
struct CouplingLsFvXD<2, FvSysEqnNS<2>> {
  using type = CouplingLsFv<2, FvSysEqnNS<2>>;
};
template <>
struct CouplingLsFvXD<3, FvSysEqnNS<3>> {
  using type = CouplingLsFv<3, FvSysEqnNS<3>>;
};
} // namespace

//--------------------------------------------------------------------------

/** \fn  Application::Application
 * \brief Gets initial data from the property file and creates the solvers and methods
 *
 * This constructor creates a list of solvers and corresponding methods.
 * The list entries are created according to the definitions in the
 * property file.
 * While the solver type is already specified here, the method type
 * will be specified in the Methods constructor (see also
 * methods.h).
 * Basic application properties are also read here.
 */
Application::Application() {
  TRACE();

  /*! \page propertiesGlobal
      \section multiSolverGrid
      <code>MInt g_multiSolverGrid</code> \n
      default = <code>0</code> \n \n
      Check if a true multi-solver simulation with one grid should be performed. This flag should be
      temporary and may be used for some multi-solver specific shenanigans.\n
      Possible values are:
      <ul>
      <li><code>0</code> (off)</li>
      <li><code>1</code> (on)</li>
      </ul>
      Keywords: <i>GENERAL, GLOBAL</i>
  */
  g_multiSolverGrid = false;
  g_multiSolverGrid = Context::getBasicProperty<MBool>("multiSolverGrid", AT_, &g_multiSolverGrid);

  /*! \page propertiesGlobal
      \section flowSolver
      <code>MBool flowSolver</code> \n
      default = <code>0</code> \n \n
      Check if flow solver is activated in property file \n
      Possible values are:
      <ul>
      <li><code>0</code> (off)</li>
      <li><code>1</code> (on)</li>
      </ul>
      Keywords: <i>GENERAL, GLOBAL</i>
  */
  const MBool flowSolverDefault = false;
  const auto flowSolver = Context::getBasicProperty<MBool>("flowSolver", AT_, &flowSolverDefault);

  globalTimeStep = 1;

  /*! \page propertiesGlobal
    \section dualTimeStepping
    <code>MBool m_dualTimeStepping</code>\n
    default = <code>false</code>\n\n
    Activates the dual-time stepping\n
    Possible values are:\n
    <ul>\n
    <li>0, 1</li>\n
    </ul>\n
    Keywords: <i>GENERAL</i>\n
  */
  m_dualTimeStepping = false;
  m_dualTimeStepping = Context::getBasicProperty<MBool>("dualTimeStepping", AT_, &m_dualTimeStepping);

  /*! \page propertiesGlobal
    \section gridGenerator
    <code>MBool gridGenerator</code> \n
    default = <code>0</code> \n \n
    Check if grid generator is activated in property file \n
    Possible values are:
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>
    Keywords: <i>GENERAL, GLOBAL</i>
  */
  const MBool gridGenDefault = false;
  const auto gridGeneration = Context::getBasicProperty<MBool>("gridGenerator", AT_, &gridGenDefault);

  // Check that either grid generation or the flow solver is activated
  if(gridGeneration && flowSolver) {
    TERMM(1, "Error: grid generation and flow solver activated.");
  }
  if(!gridGeneration && !flowSolver) {
    TERMM(1, "Error: grid generation and flow solver both deactivated.");
  }

  g_dynamicLoadBalancing = false;
  /*! \page propertiesDLB
    \section dynamicLoadBalancing
    <code>MBool Application::g_dynamicLoadBalancing </code>\n
    default = <code>0</code>\n \n
    Trigger the use of the dynamic load balancing.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, DYNAMIC, LOAD, BALANCING, TRIGGER</i>
  */
  g_dynamicLoadBalancing = Context::getBasicProperty<MBool>("dynamicLoadBalancing", AT_, &g_dynamicLoadBalancing);

  g_splitMpiComm = false;
  /*! \page propertiesGlobal
    \section splitMpiComm
    <code>MInt g_splitMpiComm </code>\n
    default = <code>0</code>\n \n
    Trigger the use of split MPI communication.\n
    Possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, MASSIVE, MPI, TRIGGER</i>
  */
  g_splitMpiComm = Context::getBasicProperty<MBool>("splitMpiComm", AT_, &g_splitMpiComm);
  if(g_splitMpiComm) {
    m_log << "Split MPI communication activated." << std::endl;
  }

  /*! \page propertiesHPC
    \section writeSolverTimings
    <code>MInt Application::m_writeSolverTimings</code>\n
    default = <code>0</code>\n \n
    Trigger the collection and output of solver timings and domain decomposition information.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, TIMINGS, PERFORMANCE, HPC, TRIGGER</i>
  */
  m_writeSolverTimings = false;
  m_writeSolverTimings = Context::getBasicProperty<MBool>("writeSolverTimings", AT_, &m_writeSolverTimings);

  if(m_writeSolverTimings) {
    /*! \page propertiesHPC
      \section solverTimingsWriteInterval
      <code>MInt Application::m_solverTimingsWriteInterval</code>\n
      default = <code>-1</code>\n \n
      Set the write interval for the solver timings output (enabled by writeSolverTimings). If
      negative the timings are only saved once at the final time step.\n
      Keywords: <i>PARALLEL, TIMINGS, PERFORMANCE</i>
    */
    m_solverTimingsWriteInterval = -1;
    m_solverTimingsWriteInterval =
        Context::getBasicProperty<MInt>("solverTimingsWriteInterval", AT_, &m_solverTimingsWriteInterval);

    /*! \page propertiesHPC
      \section solverTimingsSampleInterval
      <code>MInt Application::m_solverTimingsSampleInterval</code>\n
      default = <code>1</code>\n \n
      Set the sampling interval of the solver timings (enabled by writeSolverTimings).\n
      Keywords: <i>PARALLEL, TIMINGS, PERFORMANCE</i>
    */
    m_solverTimingsSampleInterval = 1;
    m_solverTimingsSampleInterval =
        Context::getBasicProperty<MInt>("solverTimingsSampleInterval", AT_, &m_solverTimingsSampleInterval);
    if(m_solverTimingsSampleInterval <= 0) {
      TERMM(1, "solverTimingsSampleInterval <= 0");
    }

    /*! \page propertiesHPC
      \section writeAllSolverTimings
      <code>MInt Application::m_writeAllSolverTimings</code>\n
      default = <code>true</code>\n \n
      Trigger the output of ALL solver timings and domain decomposition information.\n
      If disabled only a reduced subset of all timings is collected. This can be useful for large scale simulations to
      keep the timings files smaller and avoid storing most probably useless timing data. Keywords: <i>PARALLEL,
      TIMINGS, PERFORMANCE, TRIGGER</i>
    */
    m_writeAllSolverTimings = true;
    m_writeAllSolverTimings = Context::getBasicProperty<MBool>("writeAllSolverTimings", AT_, &m_writeAllSolverTimings);

    m_log << "Output of solver timings enabled: sampleInterval = " << m_solverTimingsSampleInterval
          << ", writeInterval = " << m_solverTimingsWriteInterval << ", allTimings = " << m_writeAllSolverTimings
          << endl;
  }

#ifdef DISABLE_OUTPUT
  const MString message = "NOTE: skipping solution and restart IO for performance testing!";
  m_log << message << std::endl;
  cerr0 << message << std::endl;
#endif


  // Get number of spatial dimensions
  const MInt nDim = read_nDim();

  m_noSolvers = 1;
  /*! \page propertiesGlobal
    \section noSolvers
    <code>MInt m_noSolvers </code> \n
    default = <code>1</code> \n \n
    Sets the number of solvers.\n
    Possible values are:
    <ul>
    <li>positive integers</li>
    </ul>
    Keywords: <i>MULTISOLVER</i>
  */
  m_noSolvers = Context::getBasicProperty<MInt>("noSolvers", AT_, &m_noSolvers);

  m_initialAdaptation = true;
  /*! \page propertiesAMR
    \section initialAdaptation
    <code>MBool m_initialAdaptation </code> \n
    default = <code>true</code> \n \n
    Activates initial adaptation. \n
    Possible values are:
    <ul>
    <li>false, true</li>
    </ul>
    Keywords: <i>MULTISOLVER, ADAPTATION</i>
  */
  m_initialAdaptation = Context::getBasicProperty<MBool>("initialAdaptation", AT_, &m_initialAdaptation);

  DEBUG("Application::() noSolvers " << m_noSolvers, MAIA_DEBUG_LEVEL1);
  if(globalDomainId() == 0) {
    cerr << "m_noSolvers: " << m_noSolvers << endl;
  }

  if(gridGeneration) {
    if(nDim == 2) {
      GridgenPar<2>(globalMaiaCommWorld(), m_noSolvers);
    } else {
      GridgenPar<3>(globalMaiaCommWorld(), m_noSolvers);
    }
    return;
  }

  // Running with true multi-solver support but just one solver does not make sense (for now)
  if(g_multiSolverGrid && m_noSolvers == 1) {
    mTerm(1, AT_, "Enabling multi-solver without having more than one solver does not make sense");
  }

  MBool addSolverToGrid = false;
  addSolverToGrid = Context::getBasicProperty<MBool>("addSolverToGrid", AT_, &addSolverToGrid);
  // multiSolverGrid needs to be enabled for 2 or more solvers (except for the case when the second solver (post-data)
  // is added during startup of the simulation)
  if(m_noSolvers > 1 && !g_multiSolverGrid && !(m_noSolvers == 2 && addSolverToGrid)) {
    TERMM(1, "Number of solvers > 1, but multiSolverGrid not enabled!");
  }

  g_timeSteps = Context::getBasicProperty<MInt>("timeSteps", AT_);
  DEBUG("Application::() g_timeSteps " << g_timeSteps, MAIA_DEBUG_LEVEL1);

  /*! \page propertiesGlobal
       \section restartBackupInterval
       <code>MInt FvCartesianSolverXD::m_restartBackupInterval</code>\n
       default = <code>25000</code>\n
       Defines the interval in which restart backup files are created\n
       Possible values are:
       <ul>
       <li>positive integers</li>
       </ul>
       Keywords: <i>FINITE_VOLUME, RESTART, I/O</i>
   */
  m_restartBackupInterval = 25000;
  m_restartBackupInterval = Context::getBasicProperty<MInt>("restartBackupInterval", AT_, &m_restartBackupInterval);

  /*! \page propertiesGlobal
    \section maxIterations
    <code>MInt Application::m_maxIterations</code>\n
    default = <code>1</code>\n \n
    Sets the maximal number of iterations\n
    Possible values are:
    <ul>
    <li>any integer n greater equal 1</li>
    </ul>
    Keywords: <i>APPLICATION</i>
  */
  m_maxIterations = 1;
  m_maxIterations = Context::getBasicProperty<MInt>("maxIterations", AT_, &m_maxIterations);

  m_solvers.resize(m_noSolvers);

  /*! \page propertiesGlobal
    \section noCouplers
    <code>MInt Application::m_noCouplers</code>\n
    default = <code>0</code>\n \n
    Sets the number of couplers\n
    Possible values are:
    <ul>
    <li>any integer n greater equal 1</li>
    </ul>
    Keywords: <i>APPLICATION</i>
  */
  m_noCouplers = 0;
  if(m_noSolvers > 1) {
    m_noCouplers = Context::getBasicProperty<MInt>("noCouplers", AT_, &m_noCouplers);
  }
  m_couplers.resize(m_noCouplers);

  // Post-Processing
  m_postProcessing = false;
  m_noPostProcessing = 0;

#ifndef DISABLE_OUTPUT
  m_postProcessing = Context::getBasicProperty<MBool>("postProcessing", AT_, &m_postProcessing);
  m_noPostProcessing = Context::getBasicProperty<MInt>("noPostProcessing", AT_, &m_noPostProcessing);
  if(m_postProcessing && m_noPostProcessing == 0) {
    TERMM(1, "postProcessing is true, but the number of postprocessing solvers is 0!");
  }
#endif

  /*! \page propertiesHPC
    \section displayMemoryStatistics
    <code>MBool Application::m_displayMemoryStatistics</code>\n
    default = <code>true</code>\n \n
    Controls if memory statistics determined and displayed for the user.\n
    Possible values are:
    <ul>
    <li>false, true</li>
    </ul>
    Keywords: <i>APPLICATION, MEMORY, HPC</i>
  */
  m_displayMemoryStatistics = true;
  m_displayMemoryStatistics =
      Context::getBasicProperty<MBool>("displayMemoryStatistics", AT_, &m_displayMemoryStatistics);

  /*! \page postProcessing
    \section ppAfterTS
    <code>MBool Application::m_ppAfterTS</code>\n
    default = <code>false</code>\n \n
    Controls the position of the post-processing in the run-loop.
    For interleaved, non-blocking coupling it can be meaningful to
    perform the in-solution post-processing before the post-time-step
    and post-couple\n
    Possible values are:
    <ul>
    <li>false, true</li>
    </ul>
    Keywords: <i>APPLICATION, HPC</i>
  */
  m_ppAfterTS = false;
  m_ppAfterTS = Context::getBasicProperty<MBool>("postProcessingAfterTS", AT_, &m_ppAfterTS);
}


Application::~Application() {
  for(auto& coupler : m_couplers) {
    coupler.reset();
  }
  m_couplers.clear();

  for(auto& solver : m_solvers) {
    solver.reset();
  }
  m_solvers.clear();
}

/// Initialize the collection of solver/coupler timings for performance evaluations
void Application::initTimings() {
  TRACE();

  // Return if not enabled
  if(!m_writeSolverTimings) {
    return;
  }

  const MInt noSolversAndCouplers = m_noSolvers + m_noCouplers;
  const MBool allTimings = m_writeAllSolverTimings;
  // Determine total number of timers and domain decomposition information of each solver/coupler
  m_noGlobalSolverTimers = 0;

  for(MInt i = 0; i < noSolversAndCouplers; i++) {
    const MInt noTimers = (i < m_noSolvers) ? m_solvers[i]->noSolverTimers(allTimings)
                                            : m_couplers[i - m_noSolvers]->noCouplingTimers(allTimings);
    m_noGlobalSolverTimers += noTimers;
  }

  MInt maxNoGlobalSolverTimers = -1;
  MPI_Allreduce(&m_noGlobalSolverTimers, &maxNoGlobalSolverTimers, 1, maia::type_traits<MInt>::mpiType(), MPI_MAX,
                globalMaiaCommWorld(), AT_, "m_noGlobalSolverTimers", "maxNoGlobalSolverTimers");
  if(maxNoGlobalSolverTimers != m_noGlobalSolverTimers) {
    TERMM(1, "Error: number of global solver timings does not match on all domains.");
  }

  m_solverTimings.clear();
  m_solverTimings.resize(m_noGlobalSolverTimers);

  m_solverTimingsPrevTime.clear();
  m_solverTimingsPrevTime.resize(m_noGlobalSolverTimers);

  m_solverTimingsTimeStep.clear();
  m_solverTimingsTimeStep.reserve(m_maxNoSolverTimings);

  // Reserve memory for timings (should be enough to avoid reallocation during solver run)
  for(MInt timerId = 0; timerId < m_noGlobalSolverTimers; timerId++) {
    m_solverTimings[timerId].clear();
    m_solverTimings[timerId].reserve(m_maxNoSolverTimings);

    m_solverTimingsPrevTime[timerId] = 0.0;
  }

  for(MInt i = 0; i < noSolversAndCouplers; i++) {
    const MInt noTimers = (i < m_noSolvers) ? m_solvers[i]->noSolverTimers(allTimings)
                                            : m_couplers[i - m_noSolvers]->noCouplingTimers(allTimings);
    // Get solver timer names
    std::vector<std::pair<MString, MFloat>> timings{};
    (i < m_noSolvers) ? m_solvers[i]->getSolverTimings(timings, allTimings)
                      : m_couplers[i - m_noSolvers]->getCouplingTimings(timings, allTimings);
    ASSERT((MInt)timings.size() == noTimers, "number of timings mismatch #" + std::to_string(i) + ": "
                                                 + std::to_string(timings.size()) + " != " + std::to_string(noTimers));

    for(MInt j = 0; j < noTimers; j++) {
      m_solverTimingsNames.push_back(timings[j].first);
    }
  }

  m_log << "Initialized collection of solver timings: " << m_noGlobalSolverTimers << std::endl;
}


/** \brief The main loop of maia.
 *
 *
 * \author Thomas Hoesgen
 * \date 10/2019
 */
template <MInt nDim>
void Application::run() {
  TRACE();
  cerr0 << endl << "=== MAIA RUN LOOP ===" << endl << endl;
  const MPI_Comm comm = globalMaiaCommWorld();
  const MFloat runTimeStart = wallTime();

  auto logDuration = [&](const MFloat timeStart, const MString comment) {
    logDuration_(timeStart, "RUN", comment, comm, globalDomainId(), globalNoDomains());
  };

  // Create some timers
  NEW_TIMER_GROUP(tg_totalMainLoop, "total main loop");
  NEW_TIMER(t_totalMainLoop, "total main loop", tg_totalMainLoop);
  NEW_SUB_TIMER(t_adaptation, "adaptation", t_totalMainLoop);
  NEW_SUB_TIMER(t_preTimeStep, "preTimeStep", t_totalMainLoop);
  NEW_SUB_TIMER(t_preCouple, "preCouple", t_totalMainLoop);
  NEW_SUB_TIMER(t_timeStep, "timeStep", t_totalMainLoop);
  NEW_SUB_TIMER(t_postTimeStep, "postTimeStep", t_totalMainLoop);
  NEW_SUB_TIMER(t_postCouple, "postCouple", t_totalMainLoop);
  NEW_SUB_TIMER(t_implicitTimeStep, "Implicit time-step", t_totalMainLoop);

#ifndef DISABLE_OUTPUT
  NEW_SUB_TIMER(t_solutionOutput, "solutionOutput", t_totalMainLoop);
#endif

  NEW_SUB_TIMER(t_balance, "balance", t_totalMainLoop);

  RECORD_TIMER_START(t_totalMainLoop);

  for(MInt i = 0; i < m_noSolvers; i++) {
    // Create separate output directory for this solver if it does not exist yet
    if(globalDomainId() == 0) {
      const MString outputDir = Context::getSolverProperty<MString>("outputDir", i, AT_);
      createDir(outputDir);
    }
  }

  const MFloat geometryTimeStart = wallTime();
  // Create Geometries (for solvers using cartesian grid)
  vector<typename GeometryXD<nDim>::type*> geometries;
  for(MInt solverId = 0; solverId < m_noSolvers; solverId++) {
    const MString solverType = Context::getSolverProperty<MString>("solvertype", solverId, AT_);
    if(gridType(solverType) == MAIA_GRID_CARTESIAN) {
      geometries.push_back(new typename GeometryXD<nDim>::type(solverId, comm));
    } else {
      geometries.push_back(nullptr);
    }
  }
  logDuration(geometryTimeStart, "Create geometries");

  // Enable auto-save feature:
  MInt noMinutesEndAutoSave = 5;
  noMinutesEndAutoSave = Context::getBasicProperty<MInt>("noMinutesEndAutoSave", AT_, &noMinutesEndAutoSave);

  MInt endAutoSaveTime = -1;
  const char* envJobEndTime = getenv("MAIA_JOB_END_TIME");
  if(envJobEndTime) {
    const MString jobEndTime(envJobEndTime);
    MBool onlyDigits = true;
    for(auto&& character : jobEndTime) {
      if(!isdigit(character)) {
        m_log << "Warning: the environment variable MAIA_JOB_END_TIME includes "
              << "non-digit characters.\n"
              << "Saving a restart file before the compute job ends is NOT active." << endl;
        onlyDigits = false;
        break;
      }
    }
    if(onlyDigits) {
      endAutoSaveTime = stoi(jobEndTime);
      endAutoSaveTime -= noMinutesEndAutoSave * 60;
      time_t time_point = endAutoSaveTime;
      m_log << "Activated automatic restart file writing at " << ctime(&time_point)
            << " (in Unix time: " << endAutoSaveTime << ")." << endl;
    }
  }

  // Seed RNG
  /*! \page propertiesHPC
    \section seedRNGWithTime
    <code>MBool MAIAApplication::run()::seedRNGWithTime</code>\n
    default = <code>false</code>\n \n
    Determines if the random number generator is seeded with time.
    <ul>
      <li>true</li>
      <li>false</li>
    </ul>
    Keywords: <i>RANDOM, RNG, SEED</i>
  */
  const MBool seedRNGWithTime = false;
  if(Context::getBasicProperty<MBool>("seedRNGWithTime", AT_, &seedRNGWithTime)) {
    if(Context::propertyExists("RNGSeed")) mTerm(1, "Property RNGSeed and seedRNGWithTime are exclusive!");
    srand(time(nullptr));
  } else {
    /*! \page propertiesHPC
      \section RNGSeed
      <code>MInt MAIAApplication::run()::seed</code>\n
      default = <code>0</code>\n \n
      Sets the seed for the random number generator.
      <ul>
        <li>any positive integer</li>
      </ul>
      Keywords: <i>RANDOM, RNG, SEED</i>
    */
    MInt seed = 0;
    seed = Context::getBasicProperty<MInt>("RNGSeed", AT_, &seed);
    if(seed < 0) mTerm(1, "RNGSeed has to be a positive integer.");
    srand(static_cast<MUint>(seed));
  }

  // Create Grid
  cerr0 << "=== Create grid..." << endl;
  const MFloat gridTimeStart = wallTime();

  CartesianGrid<nDim>* grid = nullptr;
#ifndef MAIA_DISABLE_STRUCTURED
  StructuredGrid<nDim>* gridStructured = nullptr;
#endif
  MBool structuredGridExist = false;
  MBool cartGridExist = false;

  for(MInt i = 0; i < m_noSolvers; i++) {
    const MString solverType = Context::getSolverProperty<MString>("solvertype", i, AT_);
    const MInt solverGridType = gridType(solverType);

    if(solverGridType == MAIA_GRID_STRUCTURED) {
      if(!structuredGridExist) {
        // Create Grid
        cerr0 << "=== Create structured grid..." << endl;
#ifndef MAIA_DISABLE_STRUCTURED
        gridStructured = new StructuredGrid<nDim>(i, comm);
#endif
        cerr0 << "=== done." << endl;
        structuredGridExist = true;
      }
    } else if(solverGridType == MAIA_GRID_CARTESIAN) {
      if(!cartGridExist) {
        // Create Grid
        cerr0 << "=== Create Cartesian grid..." << endl;
        // TODO: Determination of maxNoCells is already done at environment level
        MInt maxNoCells = Context::getBasicProperty<MInt>("maxNoCells", AT_);
        if(Context::propertyExists("noDomains")) {
          const MInt testNoDomains = Context::getBasicProperty<MInt>("noDomains", AT_, 0);
          if(globalNoDomains() < testNoDomains) {
            // Here, the number of maxNoCells is scaled. This is useful if a test
            // case is specified to run with a certain number of ranks
            // ('noDomains' property in run.toml) but needs to be run on a lower
            // number of mpiranks, p.e. by running maia on an accelerator.
            maxNoCells *= testNoDomains / (MFloat)globalNoDomains();
            cerr0 << "noDomain > number of used mpi ranks! Therefore, increasing maxNoCells: " << maxNoCells
                  << std::endl;
          }
        }
        grid = new CartesianGrid<nDim>(maxNoCells, (geometries[0])->boundingBox(), comm);
        cerr0 << "=== done." << endl;
        cartGridExist = true;
      }
    } else if(solverGridType == MAIA_GRID_NONE) {
      continue;
    } else {
      TERMM(1, "Invalid grid type: " + std::to_string(solverGridType));
    }
  }
  logDuration(gridTimeStart, "Create grid");

  // Read in the propertiesGroups, which were before set in initMehtods()
  MBool* propertiesGroups = readPropertiesGroups();
  std::vector<MBool> isActive(m_noSolvers, false);

  // Create Solvers
  const MFloat createSolversTimeStart = wallTime();
  cerr0 << "=== Create solvers..." << endl;

  for(MInt i = 0; i < m_noSolvers; i++) {
    const MFloat createSolverTimeStart = wallTime();
    MBool active = false;
    const MString solverType = Context::getSolverProperty<MString>("solvertype", i, AT_);
    const MInt solverGridType = gridType(solverType);

    if(solverGridType == MAIA_GRID_STRUCTURED) {
      ASSERT(structuredGridExist, "");
      cerr0 << "=== Create structured solver..." << endl;
#ifndef MAIA_DISABLE_STRUCTURED
      createSolver<nDim>(gridStructured, i, propertiesGroups, comm);
#endif
      cerr0 << "=== done." << endl;
    } else if(solverGridType == MAIA_GRID_CARTESIAN) {
      ASSERT(cartGridExist, "");
      cerr0 << "=== Create " << solverType << "..." << endl;
      createSolver<nDim>(grid, i, geometries[i], propertiesGroups, active);
      cerr0 << "=== done." << endl;
    } else if(solverGridType == MAIA_GRID_NONE) {
      cerr0 << "=== Create gridless solver " << solverType << "..." << endl;
      createSolver<nDim>(i, comm);
      cerr0 << "=== done." << endl;
    } else {
      TERMM(1, "Invalid grid type: " + std::to_string(solverGridType));
    }

    isActive[i] = active;
    logDuration(createSolverTimeStart, "Create solver #" + std::to_string(i));
  }
  logDuration(createSolversTimeStart, "Create solvers");

  // No longer needed at that point
  delete[] propertiesGroups;

  // Create solver couplers
  const MFloat couplersTimeStart = wallTime();
  cerr0 << "=== Create couplers..." << endl;
  for(MInt i = 0; i < m_noCouplers; i++) {
    createCoupler<nDim>(i);
  }
  cerr0 << "=== done." << endl;
  logDuration(couplersTimeStart, "Create couplers");

  globalTimeStep = 0;
  MInt restartTimeStep = -1;

  MBool restartFile = false;
  restartFile = Context::getBasicProperty<MBool>("restartFile", AT_, &restartFile);

  // Set the correct global time step (at a restart)
  globalTimeStep = getInitialTimeStep(restartFile, comm);
  restartTimeStep = globalTimeStep;
  g_restartTimeStep = globalTimeStep;

  // Create grid controller
  const MFloat controllerTimeStart = wallTime();
  cerr0 << "=== Create grid controller..." << endl;
  typename ::maia::grid::Controller<nDim> controller(grid, &m_solvers, &m_couplers);
  cerr0 << "=== done." << endl;
  logDuration(controllerTimeStart, "Create grid controller");

  // Create PostProcessing controller
  PostProcessingController<nDim>* ppController = nullptr;
  std::vector<PostProcessingInterface*> pp;

  // Init solvers - solverId in propertiesFile determines solver order
  const MFloat initSolversTimeStart = wallTime();
  cerr0 << "=== Init solvers..." << endl;
  for(MInt i = 0; i < m_noSolvers; i++) {
    const MFloat initSolverTimeStart = wallTime();
    m_solvers[i]->initSolver();
    logDuration(initSolverTimeStart, "Init solver #" + std::to_string(i));

    if(m_displayMemoryStatistics) {
      writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "solver init #" + std::to_string(i));
    }
  }
  logDuration(initSolversTimeStart, "Init solvers");
  cerr0 << "=== done." << endl;


  // Init couplers - couplerId in propertiesFile determines coupler order
  const MFloat initCouplersTimeStart = wallTime();
  cerr0 << "=== Init couplers..." << endl;
  for(auto&& coupler : m_couplers) {
    coupler->init();
  }
  logDuration(initCouplersTimeStart, "Init couplers");
  if(m_displayMemoryStatistics) {
    writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "coupler init");
  }

  // Switch to fully disable/ignore all DLB timers
  MBool ignoreDlbTimers = false;
  ignoreDlbTimers = Context::getBasicProperty<MBool>("ignoreDlbTimers", AT_, &ignoreDlbTimers);
  // Prevent disabled DLB timers to be used for any performance output since its useless in this case
  if(ignoreDlbTimers) {
    const MBool performanceOutput = Context::getBasicProperty<MBool>("performanceOutput", AT_);
    TERMM_IF_COND(performanceOutput,
                  "ERROR: ignoreDlbTimers=true should only be used together with performanceOutput=false");
  }

  // Create DLB timers for solvers and couplers
  maia::dlb::g_dlbTimerController.createDlbTimers(m_noSolvers + m_noCouplers, ignoreDlbTimers);
  m_log << "Created " << m_noSolvers + m_noCouplers << " DLB timers (ignore = " << ignoreDlbTimers << ")" << std::endl;

  // Assign each solver/coupler its DLB timer id
  for(MInt i = 0; i < m_noSolvers; i++) {
    m_solvers[i]->setDlbTimer(i);
  }
  for(MInt i = 0; i < m_noCouplers; i++) {
    m_couplers[i]->setDlbTimer(m_noSolvers + i);
  }

  cerr0 << "=== done." << endl;

  MBool forceInitialAdaptation = false;
  for(MInt solverId = 0; solverId < m_noSolvers; solverId++) {
    if(m_solvers[solverId]->forceAdaptation()) {
      forceInitialAdaptation = true;
      cerr0 << "=== Forcing Initial adaptation: " << endl;
      break;
    }
  }

  // Initial adaptation - can be turned off by setting initialAdaptation = false
  if(!restartFile || forceInitialAdaptation) {
    cerr0 << "=== Initial adaptation..." << endl;
    const MInt gtbak = globalTimeStep;
    globalTimeStep = -1;
    RECORD_TIMER_START(t_adaptation);
    controller.adaptation(m_initialAdaptation);
    RECORD_TIMER_STOP(t_adaptation);
    globalTimeStep = gtbak;

    MBool outputInitialAdaptation = false;
    outputInitialAdaptation =
        Context::getBasicProperty<MBool>("outputInitialAdaptation", AT_, &outputInitialAdaptation);
    if(!restartFile && outputInitialAdaptation) {
      controller.writeRestartFile(true, false);
    }
    cerr0 << "=== done." << endl;
  }

  // Finalize solvers
  const MFloat finalizeInitTimeStart = wallTime();
  cerr0 << "=== Finalize initialization of solvers and couplers..." << endl;

  for(MInt i = 0; i < m_noSolvers; i++) {
    const MFloat finalizeInitSolverTimeStart = wallTime();
    m_solvers[i]->finalizeInitSolver();
    if(geometries[i] != nullptr) {
      geometries[i]->logStatistics();
    }
    // Some solvers need coupler information before they can be finalized!
    // NOTE: If couplers were to know solverIds of the solvers they couple, a call of each coupler after each finalize
    // could be prevented!
    for(auto&& coupler : m_couplers) {
      coupler->finalizeSubCoupleInit(i);
    }
    logDuration(finalizeInitSolverTimeStart, "Finalize initialization solver #" + std::to_string(i));
  }
  logDuration(finalizeInitTimeStart, "Finalize initialization");

  // Finalize couplers
  const MFloat finalizeCouplerInitTimeStart = wallTime();
  for(auto&& coupler : m_couplers) {
    coupler->finalizeCouplerInit();
  }
  cerr0 << "=== done." << endl;


  logDuration(finalizeCouplerInitTimeStart, "Finalize coupler init");

  if(m_postProcessing) {
    pp.clear();
    cerr0 << "=== Init postprocessing..." << endl;
    for(MInt ppId = 0; ppId < (MInt)m_noPostProcessing; ppId++) {
      pp.push_back(createPostProcessing<nDim>(ppId));
    }

    ppController = new PostProcessingController<nDim>(pp);

    ppController->init();
    ppController->preSolve();

    cerr0 << "=== done." << endl;
  }

  MBool outputInitialCondition = false;
  outputInitialCondition = Context::getBasicProperty<MBool>("outputInitialCondition", AT_, &outputInitialCondition);
  MBool writePostprocessingPre = false;
  writePostprocessingPre = Context::getBasicProperty<MBool>("writePostprocessingPre", AT_, &writePostprocessingPre);
  const MBool initialBackup = globalTimeStep > 0 && globalTimeStep % m_restartBackupInterval == 0;
  if((outputInitialCondition && (!restartFile || forceInitialAdaptation)) || writePostprocessingPre) {
    controller.writeRestartFile(true, initialBackup);
  }

  // Update partition cell workloads if enabled and save to the grid file, cleanup, then quit.
  if(controller.updateGridPartitionWorkloads()) {
    cleanUp();
    return;
  }

  // Initialize timings for dynamic load balancing
  initTimings();

  // Create execution recipe in executionrecipe.h
  ExecutionRecipe* recipe = createRecipe();

  // Memory statistics for each process
  if(m_displayMemoryStatistics) {
    writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "After initialization");
  }

  // Alive output during simulation
  /*! \page propertiesHPC
    \section aliveInterval
    <code>MInt Application::run_unified()::aliveInterval</code>\n
    default = <code>0 (disabled)</code>\n \n
    Set the number of time steps after which performance information is given
    to the user. If the alive interval is set to 0, no alive signals are given.\n
    Possible values are:
    <ul>
      <li>any integer >= 0</li>
    </ul>
    Keywords: <i>GENERAL, PERFORMANCE, HPC, OUTPUT</i>
  */
  MInt aliveInterval = 0;
  aliveInterval = Context::getBasicProperty<MInt>("aliveInterval", AT_, &aliveInterval);
  if(aliveInterval < 0) {
    TERMM(1, "Alive interval must be >= 0 (is: " + to_string(aliveInterval) + ").");
  }
  const MFloat loopTimeStart = wallTime();
  MFloat lastRunTime = 0.0;
  MInt lastGlobalTimeStep = globalTimeStep;

  /// Time-step loop (loop over all time steps):
  MBool finalTimeStep = (globalTimeStep) >= (g_timeSteps + restartTimeStep);

  MBool debugOutputAfterBalance = false;
  debugOutputAfterBalance = Context::getBasicProperty<MBool>("debugOutputAfterBalance", AT_, &debugOutputAfterBalance);

  // set finale timeStep if auto-save feature is enabled!
  if(endAutoSaveTime != -1
     && endAutoSaveTime
            <= chrono::duration_cast<chrono::seconds>(chrono::system_clock::now().time_since_epoch()).count()) {
    finalTimeStep = true;
  }

  logDuration(runTimeStart, "Full initialization");

  // Enable DLB timers
  maia::dlb::g_dlbTimerController.enableAllDlbTimers();

  while(globalTimeStep < (g_timeSteps + restartTimeStep)) {
    globalTimeStep++;
    finalTimeStep = (globalTimeStep) >= (g_timeSteps + restartTimeStep);

    // Advance time step loop
    MBool advanceTimeStep = false;
    while(!advanceTimeStep) {
      // Check whether a solver adaptation is allowed during a specific iteration step
      if(recipe->callAdaptation()) {

        MBool force = false;
        for(MInt solverId = 0; solverId < m_noSolvers; solverId++) {
          if(m_solvers[solverId]->forceAdaptation()) {
            force = true;
            m_log << "Solver " << solverId << " is forcing a mesh-adaptation at time-step " << globalTimeStep << endl;
            break;
          }
        }

        maia::dlb::g_dlbTimerController.disableAllDlbTimers();

        if(force || controller.isAdaptationTimeStep()) {
          storeTimingsAndSolverInformation(true);
        }

        RECORD_TIMER_START(t_adaptation);
        const MBool wasAdapted = controller.adaptation(force);
        RECORD_TIMER_STOP(t_adaptation);

        if(wasAdapted && m_displayMemoryStatistics) {
          writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "After adaptation");
        }

        maia::dlb::g_dlbTimerController.enableAllDlbTimers();
      }

      // Call preTimeStep() in each solver via the executionrecipe
      RECORD_TIMER_START(t_preTimeStep);
      recipe->preTimeStep();
      RECORD_TIMER_STOP(t_preTimeStep);

      // Call preCouple() in each coupler via the executionrecipe
      RECORD_TIMER_START(t_preCouple);
      recipe->preCouple();
      RECORD_TIMER_STOP(t_preCouple);

      // Advance the time step according to the specified recipe
      RECORD_TIMER_START(t_timeStep);
      recipe->timeStep();
      RECORD_TIMER_STOP(t_timeStep);

      // Post-Processing InSolve before postTS and postCouple
      if(m_postProcessing && m_ppAfterTS) {
        ppController->setStep(recipe->a_step());
        ppController->inSolve(finalTimeStep);
      }

      // Call postTimeStep() in each solver via the executionrecipe
      RECORD_TIMER_START(t_postTimeStep);
      recipe->postTimeStep();
      RECORD_TIMER_STOP(t_postTimeStep);

      // Call postCouple() in each coupler via the executionrecipe
      RECORD_TIMER_START(t_postCouple);
      recipe->postCouple();
      RECORD_TIMER_STOP(t_postCouple);

      // Post-Processing InSolve after postTS and postCoupl
      if(m_postProcessing && !m_ppAfterTS) {
        ppController->setStep(recipe->a_step());
        ppController->inSolve(finalTimeStep);
      }

      // Determine wich solvers are active during each iterarion step
      advanceTimeStep = recipe->updateCallOrder();
      // URL
    }

    // Dual time stepping - probably not working!
    if(m_dualTimeStepping) {
      RECORD_TIMER_START(t_implicitTimeStep);
      for(MInt i = 0; i < m_noSolvers; i++) {
        m_solvers[i]->implicitTimeStep();
      }
      RECORD_TIMER_STOP(t_implicitTimeStep);
      globalTimeStep++;
    }

    // Check if each solvers solution is converged. If so, write a restartFile and exit time step loop.
    MBool completed = true;
    for(MInt i = 0; i < m_noSolvers; i++) {
      // Note: solverConverged() is a prototype that only exists in solver.h and always returns false
      // ansgar: for the DG solver, solverConverged() returns true if the final solution time
      // T_end is reached
      if(m_solvers[i]->isActive()) {
        // TODO labels:totest if there are inactive ranks for a solver these do not have the converged info!
        completed &= m_solvers[i]->solverConverged();
      }
    }

    // TODO labels:noissue this is just to ensure that in a multisolver computation with inactive ranks the
    // converged info does not lead to ranks exiting the run loop while others continue!
    if(m_noSolvers > 1) {
#ifndef NDEBUG
      if(completed) {
        std::cerr << "Warning: one or multiple solvers in a multisolver computation returned a "
                     "'completed' status, however this information might not be consistent among "
                     "ranks if there are inactive ranks!"
                  << std::endl;
      }
#endif

      completed = false;
    }

    finalTimeStep = finalTimeStep || completed;

    // set finale timeStep if auto-save feature is enabled!
    if(endAutoSaveTime != -1
       && endAutoSaveTime
              <= chrono::duration_cast<chrono::seconds>(chrono::system_clock::now().time_since_epoch()).count()) {
      finalTimeStep = true;

      time_t now =
          std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

      m_log << "Finished end auto save at " << std::ctime(&now) << " (in Unix time: " << now << ")." << endl;
      time_t remaining_time = endAutoSaveTime - now;
      m_log << "Job will end at" << std::ctime(&remaining_time) << " (in Unix time: " << remaining_time << ")." << endl;
    }

#ifndef DISABLE_OUTPUT
    const MBool forceOutput = false;
    RECORD_TIMER_START(t_solutionOutput);

    // Write the restart files via the gridController
    const MBool backup = (globalTimeStep % m_restartBackupInterval == 0) && (globalTimeStep > 0);
    controller.writeRestartFile(finalTimeStep, backup);
    //    extRestartFile(finalTimeStep);

    // Solver specific solution output
    // TODO labels:GRID save also adapted grid if there is no restart file written at the same time and update the
    // grid file name
    for(MInt i = 0; i < m_noSolvers; i++) {
      if(m_solvers[i]->isActive()) {
        m_solvers[i]->disableDlbTimers();
        m_solvers[i]->saveSolverSolution(forceOutput, finalTimeStep);
        m_solvers[i]->enableDlbTimers();
      }
    }

    // Post-Processing InSolve write solution files
    if(m_postProcessing) {
      ppController->ppSolution();
    }

    RECORD_TIMER_STOP(t_solutionOutput);
#endif

    // Save timings for dynamic load balancing for this time step
    maia::dlb::g_dlbTimerController.disableAllDlbTimers();
    collectTimingsAndSolverInformation(finalTimeStep);
    maia::dlb::g_dlbTimerController.enableAllDlbTimers();

    // If all solvers are convered the time step loop is terminated
    // if (completed) { break; } // Note: moved to end of loop to allow DLB imbalance evaluation

    // Prepare balance and next time step
    // E.g. fvMbSolver needs to update his old variables with the new variables after it wrote its restartFile but
    // before load balancing can be performed
    for(MInt i = 0; i < m_noSolvers; i++) {
      m_solvers[i]->prepareNextTimeStep();

      // log statistics to check the balance/distribution
      // if(i == 0 && globalTimeStep % 50 == 0) {
      //  m_solvers[i]->disableDlbTimers();
      //  controller.logTimerStatistics("regular at timeStep");
      //  m_solvers[i]->enableDlbTimers();
      //}
    }

    // Load balancing
    maia::dlb::g_dlbTimerController.disableAllDlbTimers();
    RECORD_TIMER_START(t_balance);
    if(controller.isDlbTimeStep()) {
      storeTimingsAndSolverInformation(true);
    }

    const MBool wasBalanced = controller.balance(false, finalTimeStep, false);
    RECORD_TIMER_STOP(t_balance);

    if(wasBalanced && m_displayMemoryStatistics) {
      writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "After load balancing");
    }

    maia::dlb::g_dlbTimerController.enableAllDlbTimers();

    // Alive output during simulation
    if(aliveInterval > 0 && globalTimeStep % aliveInterval == 0 && globalDomainId() == 0) {
      // Calculate total run time for printing
      const MFloat runTimeTotal = wallTime() - loopTimeStart;
      const MFloat runTimePerStep =
          (runTimeTotal - lastRunTime) / static_cast<MFloat>(globalTimeStep - lastGlobalTimeStep);
      lastRunTime = runTimeTotal;
      lastGlobalTimeStep = globalTimeStep;

      // Print out processing information to user
      printf("#t/s: %8d | Run time: %.4e s | Time/step: %.4e s\n", globalTimeStep, runTimeTotal, runTimePerStep);
    }

    if(debugOutputAfterBalance && wasBalanced) {
      if(globalDomainId() == 0) {
        printf("------------Write debug restart after balance------------");
      }
      controller.writeRestartFile(true, false);
    }

    // If all solvers are converged the time step loop is terminated
    if(completed || finalTimeStep) {
      break;
    }
  } // end globalTimeStep loop

  // Disable DLB timers for remaining steps
  maia::dlb::g_dlbTimerController.disableAllDlbTimers();

  writeMemoryStatistics(comm, globalNoDomains(), globalDomainId(), AT_, "After run loop");

  // If grid was balanced save final partition file for restarting
  controller.savePartitionFile();

  // Post-Processing PostSolve
  // @Jannik: this must happen before the cleanUp, where data-structure might be destroyed
  //          and destructors called!
  if(m_postProcessing) {
    ppController->postSolve();
  }

  cleanUp();

  RECORD_TIMER_STOP(t_totalMainLoop);
}


/// \brief Return the initial global time step
MInt Application::getInitialTimeStep(const MBool restartFile, const MPI_Comm comm) {
  TRACE();
  if(!restartFile) {
    return 0;
  } else {
    MBool useNonSpecifiedRestartFile = false;
    useNonSpecifiedRestartFile =
        Context::getBasicProperty<MBool>("useNonSpecifiedRestartFile", AT_, &useNonSpecifiedRestartFile);

    if(!useNonSpecifiedRestartFile) {
      // Use given restart time step from properties
      return Context::getBasicProperty<MInt>("restartTimeStep", AT_);
    } else {
      // Restart time step not specified in properties, use the first suitable solver on rank 0 to
      // load time step from restart file and broadcast result
      MInt timeStep = -1;
      // TODO labels:totest,noissue this will not work if rank 0 does not have a solver with restart time step
      // information!
      if(globalDomainId() == 0) {
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(m_solvers[solver]->hasRestartTimeStep()) {
            timeStep = m_solvers[solver]->determineRestartTimeStep();
            break;
          }
        }
      }
      MPI_Bcast(&timeStep, 1, maia::type_traits<MInt>::mpiType(), 0, comm, AT_, "timeStep");
      m_log << "Determined global time step from restart file: " << timeStep << std::endl;
      return timeStep;
    }
  }
}


/** \brief Reads in the properties groups
 * \author Thomas Hoesgen
 * \date 10/2019
 */
MBool* Application::readPropertiesGroups() {
  TRACE();

  MBool tmpFalse = false;

  auto* propertiesGroups = new MBool[PROPERTIESGROUPS_COUNT];

  fill(propertiesGroups, propertiesGroups + PROPERTIESGROUPS_COUNT, tmpFalse);

  propertiesGroups[LEVELSETMB] = Context::getBasicProperty<MBool>("levelSetMb", AT_, &tmpFalse);
  propertiesGroups[LEVELSET] = Context::getBasicProperty<MBool>("levelSet", AT_, &tmpFalse);
  propertiesGroups[LB] = Context::getBasicProperty<MBool>("lb", AT_, &tmpFalse);
  propertiesGroups[LS_SOLVER] = Context::getBasicProperty<MBool>("lsSolver", AT_, &tmpFalse);
  propertiesGroups[COMBUSTION] = Context::getBasicProperty<MBool>("combustion", AT_, &tmpFalse);
  propertiesGroups[DG] = Context::getBasicProperty<MBool>("DG", AT_, &tmpFalse);
  propertiesGroups[LS_RANS] = Context::getBasicProperty<MBool>("levelSetRans", AT_, &tmpFalse);

  MString couplerTypeDefault = "";
  MString couplerType = Context::getBasicProperty<MString>("couplerType_0", AT_, &couplerTypeDefault);
  propertiesGroups[LEVELSET_LB] = couplerType == "COUPLER_LS_LB";

  if(globalDomainId() == 0) {
    cerr << "FV-MB is " << ((propertiesGroups[LEVELSETMB]) ? "on" : "off") << endl;
    cerr << "FV-LS is " << ((propertiesGroups[LEVELSET]) ? "on" : "off") << endl;
    cerr << "LB is " << ((propertiesGroups[LB]) ? "on" : "off") << endl;
    cerr << "LB LS is " << ((propertiesGroups[LEVELSET_LB]) ? "on" : "off") << endl;
    cerr << "LS SOLVER is " << ((propertiesGroups[LS_SOLVER]) ? "on" : "off") << endl;
    cerr << "COMBUSTION is " << ((propertiesGroups[COMBUSTION]) ? "on" : "off") << endl;
    cerr << "DG is " << ((propertiesGroups[DG]) ? "on" : "off") << endl;
    cerr << "LS-RANS is " << ((propertiesGroups[LS_RANS]) ? "on" : "off") << endl;
  }

  return propertiesGroups;
}


/// \brief Return the type of grid for a given solver type
MInt Application::gridType(const MString solverType) {
  switch(string2enum(solverType)) {
    case MAIA_STRUCTURED: {
      return MAIA_GRID_STRUCTURED;
    }
    case MAIA_ACOUSTIC_ANALOGY: {
      return MAIA_GRID_NONE; // No grid required for acoustic extrapolation (e.g. FWH)
    }
    default: { // Default Cartesian grid
      return MAIA_GRID_CARTESIAN;
    }
  }
}


/** \brief This function handels the creation of new solvers
 * \author Thomas Hoesgen
 * \date 10/2019
 */
template <MInt nDim>
void Application::createSolver(CartesianGrid<nDim>* grid, const MInt solverId, Geometry<nDim>* geometry,
                               MBool* propertiesGroups, MBool& isActive) {
  TRACE();


  // Determine solver type from property file
  const MString solverType = Context::getSolverProperty<MString>("solvertype", solverId, AT_);

  grid::Proxy<nDim>* gridProxy = nullptr;
  gridProxy = new grid::Proxy<nDim>(solverId, *grid, *geometry);
  isActive = gridProxy->isActive();

  // Now create the specified solver
  switch(string2enum(solverType)) {
    case MAIA_LEVELSET_SOLVER: {
      m_solvers[solverId] = LsCartesianSolverFactory<nDim>::create(solverId, propertiesGroups, *gridProxy, *geometry,
                                                                   gridProxy->mpiComm());
      break;
    }
    case MAIA_FINITE_VOLUME: {
      MInt noSpecies = 0;
      noSpecies = Context::getSolverProperty<MInt>("noSpecies", solverId, AT_, &noSpecies);
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solverId)) {
        fvSystemEquations = string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solverId, AT_));
      }
      MString ransMethod = "";
      if(Context::propertyExists("ransMethod", solverId)) {
        ransMethod = Context::getSolverProperty<MString>("ransMethod", solverId, AT_);
      }
      // cerr << ransMethod << endl;
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          switch(string2enum(ransMethod)) {
            case RANS_SA:
            case RANS_SA_DV: {
              m_solvers[solverId] =
                  make_unique<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
                      solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
              break;
            }
            case RANS_FS: {
              m_solvers[solverId] =
                  make_unique<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>>(
                      solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
              break;
            }
            case RANS_KOMEGA: {
              m_solvers[solverId] =
                  make_unique<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_KOMEGA>>>>(
                      solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
              break;
            }
            default: {
              TERMM(1, "Unknown RANS model for solver " + to_string(solverId));
            }
          }
          break;
        }
        case FV_SYSEQN_EEGAS: {
          IF_CONSTEXPR(nDim == 2)
          mTerm(1, AT_, "FVSysEqnEEGas is not tested for 2D at all and probably does not make much sense!");
          else m_solvers[solverId] = make_unique<FvCartesianSolverXD<nDim, FvSysEqnEEGas<nDim>>>(
              solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        case FV_SYSEQN_NS: {
          m_solvers[solverId] = make_unique<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>>(
              solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        case FV_SYSEQN_DETCHEM: {
          m_solvers[solverId] = make_unique<FvCartesianSolverXD<nDim, FvSysEqnDetChem<nDim>>>(
              solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case MAIA_FV_APE: {
      // TODO labels:FV switch to FV-solver with fvSystemEquations=APE?
      m_solvers[solverId] = make_unique<typename FvApeSolverXD<nDim>::type>(solverId, 0, propertiesGroups, *gridProxy,
                                                                            *geometry, gridProxy->mpiComm());
      break;
    }
    case MAIA_FV_MB: {
      MInt noSpecies = 0;
      noSpecies = Context::getSolverProperty<MInt>("noSpecies", solverId, AT_, &noSpecies);
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solverId)) {
        fvSystemEquations = string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solverId, AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          m_solvers[solverId] =
              make_unique<FvMbCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
                  solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        case FV_SYSEQN_NS: {
          m_solvers[solverId] = make_unique<FvMbCartesianSolverXD<nDim, FvSysEqnNS<nDim>>>(
              solverId, noSpecies, propertiesGroups, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case MAIA_DISCONTINUOUS_GALERKIN: {
      const MInt dgSystemEquations =
          string2enum(Context::getSolverProperty<MString>("dgSystemEquations", solverId, AT_));
      switch(dgSystemEquations) {
        case DG_SYSEQN_ACOUSTICPERTURB: {
          m_solvers[solverId] = make_unique<DgCartesianSolver<nDim, DgSysEqnAcousticPerturb<nDim>>>(
              solverId, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        case DG_SYSEQN_LINEARSCALARADV: {
          m_solvers[solverId] = make_unique<DgCartesianSolver<nDim, DgSysEqnLinearScalarAdv<nDim>>>(
              solverId, *gridProxy, *geometry, gridProxy->mpiComm());
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case MAIA_LATTICE_BOLTZMANN: {
      m_solvers[solverId] =
          maia::lb::LbSolverFactory<nDim>::create(solverId, *gridProxy, *geometry, gridProxy->mpiComm());
      break;
    }
    case MAIA_FINITE_CELL: {
      m_solvers[solverId] = make_unique<FcSolver<nDim>>(solverId, *gridProxy, *geometry, gridProxy->mpiComm());
      break;
    }
    case MAIA_PARTICLE: {
      IF_CONSTEXPR(nDim == 3) {
        m_solvers[solverId] = make_unique<LPT<nDim>>(solverId, *gridProxy, *geometry, gridProxy->mpiComm());
      }
      else {
        TERMM(-1, "LPT not currently supported in 2D");
      }
      break;
    }
    case MAIA_RIGID_BODIES: {
      m_solvers[solverId] = make_unique<RigidBodies<nDim>>(solverId, *gridProxy, *geometry, gridProxy->mpiComm());
      break;
    }
    case MAIA_POST_DATA: {
      m_solvers[solverId] = make_unique<PostData<nDim>>(solverId, *gridProxy, *geometry, gridProxy->mpiComm());
      if(m_postDataSolverId != -1) TERMM(1, "Currently not more than 1 PostData possible!");
      m_postDataSolverId = solverId;
      break;
    }
    default: {
      TERMM(1, "Unknown solver type '" + solverType + "', exiting ... ");
    }
  }

  m_log << "Created solver #" << std::to_string(solverId) << ": " << solverType
        << m_solvers[solverId]->getIdentifier(false, " with alias: ", "") << std::endl;
}

/** \brief This function handles the creation of new structured solvers
 * \author Julian Stemmermann
 * \date 05/2020
 */
#ifndef MAIA_DISABLE_STRUCTURED
template <MInt nDim>
void Application::createSolver(StructuredGrid<nDim>* grid, const MInt solverId, MBool* propertiesGroups,
                               const MPI_Comm comm) {
  TRACE();


  // Determine solver type from property file
  const MString solverType = Context::getSolverProperty<MString>("solvertype", solverId, AT_);

  // Now create the specified solver
  switch(string2enum(solverType)) {
    case MAIA_STRUCTURED: {
      m_solvers[solverId] =
          make_unique<typename FvStructuredSolverXD<nDim>::type>(solverId, grid, propertiesGroups, comm);
      break;
    }
    default: {
      TERMM(1, "Unknown solver type '" + solverType + "', exiting ... ");
    }
  }

  m_log << "Created solver #" << std::to_string(solverId) << ": " << solverType
        << m_solvers[solverId]->getIdentifier(false, " with alias: ", "") << std::endl;
}
#endif


/// \brief Handles the creation of gridless solvers
template <MInt nDim>
void Application::createSolver(const MInt solverId, const MPI_Comm comm) {
  TRACE();

  // Determine solver type from property file
  const MString solverType = Context::getSolverProperty<MString>("solvertype", solverId, AT_);

  switch(string2enum(solverType)) {
    case MAIA_ACOUSTIC_ANALOGY: {
      m_solvers[solverId] = make_unique<AcaSolver<nDim>>(solverId, comm);
      break;
    }
    default: {
      TERMM(1, "Unknown solver type '" + solverType + "', exiting ... ");
    }
  }

  m_log << "Created solver #" << std::to_string(solverId) << ": " << solverType
        << m_solvers[solverId]->getIdentifier(false, " with alias: ", "") << std::endl;
}


/** \brief This function handels the creation of new couplers
 * \author Thomas Hoesgen
 * \date 10/2019
 */
template <MInt nDim>
void Application::createCoupler(MInt couplerId) {
  TRACE();

  const MString couplerType = Context::getBasicProperty<MString>("couplerType_" + std::to_string(couplerId), AT_);

  cerr0 << "=== Create coupler #" << couplerId << " " << couplerType << std::endl;

  // Determine which solvers the coupler is supposed to couple
  // TODO labels:COUPLER allow more than 2 solvers to be coupled, e.g. for multilevel
  MInt solverToCoupleLength = Context::propertyLength("solversToCouple_" + std::to_string(couplerId));
  std::vector<MInt> solversToCouple;
  for(MInt solver = 0; solver < solverToCoupleLength; solver++) {
    solversToCouple.push_back(-1);
    if(Context::propertyExists("solversToCouple_" + std::to_string(couplerId))) {
      solversToCouple[solver] =
          Context::getBasicProperty<MInt>("solversToCouple_" + std::to_string(couplerId), AT_, nullptr, solver);
      if(solversToCouple[solver] < 0 || solversToCouple[solver] >= m_noSolvers) {
        TERMM(1, "Invalid solver id: " + std::to_string(solversToCouple[solver]));
      }
    } else {
      TERMM(1, "solversToCouple_" + std::to_string(couplerId) + " has to be specified!");
    }
  }

  // Create the coupler
  //
  switch(string2enum(couplerType)) {
    case COUPLER_LB_DG_APE: {
      const MInt nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[0], AT_);
      switch(nDist) {
        case 19: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          const auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          const auto dgSolver =
              static_cast<DgCartesianSolver<3, DgSysEqnAcousticPerturb<3>>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbDgApe<3, 19, SysEqnLb>>(couplerId, lbSolver, dgSolver);
          break;
        }
        case 27: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          const auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          const auto dgSolver =
              static_cast<DgCartesianSolver<3, DgSysEqnAcousticPerturb<3>>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbDgApe<3, 27, SysEqnLb>>(couplerId, lbSolver, dgSolver);
          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions! Only working for q19 and q27, yet.");
        }
      }
      break;
    }
    case COUPLER_LS_FV_MB: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          m_couplers[couplerId] =
              make_unique<typename LsFvMbXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>::type>(
                  couplerId,
                  static_cast<typename LsCartesianSolverXD<nDim>::type*>(m_solvers[solversToCouple[0]].get()),
                  static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                      m_solvers[solversToCouple[1]].get()));
          break;
        }
        case FV_SYSEQN_NS:
        default: {
          m_couplers[couplerId] = make_unique<typename LsFvMbXD<nDim, FvSysEqnNS<nDim>>::type>(
              couplerId,
              static_cast<typename LsCartesianSolverXD<nDim>::type*>(m_solvers[solversToCouple[0]].get()),
              static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
      }

      break;
    }
    case COUPLER_FV_ZONAL_RTV: {
      MString ransMethod = "";
      if(Context::propertyExists("ransMethod", solversToCouple[0])) {
        ransMethod = Context::getSolverProperty<MString>("ransMethod", solversToCouple[0], AT_);
      }
      switch(string2enum(ransMethod)) {
        case RANS_SA:
        case RANS_SA_DV: {
          m_couplers[couplerId] = make_unique<FvZonalRTV<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
              couplerId,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                  m_solvers[solversToCouple[0]].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
        case RANS_FS: {
          m_couplers[couplerId] = make_unique<FvZonalRTV<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>>(
              couplerId,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>*>(
                  m_solvers[solversToCouple[0]].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
        default: {
          TERMM(1, "Unknown RANS model for solver " + to_string(solversToCouple[0]));
        }
      }
      break;
    }
    case COUPLER_FV_ZONAL_STG: {
      MString ransMethod = "";
      if(Context::propertyExists("ransMethod", solversToCouple[0])) {
        ransMethod = Context::getSolverProperty<MString>("ransMethod", solversToCouple[0], AT_);
      }
      switch(string2enum(ransMethod)) {
        case RANS_SA:
        case RANS_SA_DV: {
          m_couplers[couplerId] = make_unique<FvZonalSTG<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
              couplerId,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                  m_solvers[solversToCouple[0]].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
        case RANS_FS: {
          m_couplers[couplerId] = make_unique<FvZonalSTG<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>>(
              couplerId,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>*>(
                  m_solvers[solversToCouple[0]].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
        default: {
          TERMM(1, "Unknown RANS model for solver " + to_string(solversToCouple[0]));
        }
      }
      break;
    }
    case COUPLER_FV_MB_ZONAL: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          m_couplers[couplerId] =
              make_unique<CouplerFvMbZonal<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
                  couplerId,
                  static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                      m_solvers[solversToCouple[0]].get()),
                  static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                      m_solvers[solversToCouple[1]].get()));
          break;
        }
        case FV_SYSEQN_NS:
        default: {
          m_couplers[couplerId] = make_unique<CouplerFvMbZonal<nDim, FvSysEqnNS<nDim>>>(
              couplerId,
              static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[0]].get()),
              static_cast<FvMbCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
      }
      break;
    }

    case COUPLER_CARTESIAN_INTERPOLATION: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          m_couplers[couplerId] =
              make_unique<FvCartesianInterpolation<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>,
                                                   FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
                  couplerId,
                  static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                      m_solvers[solversToCouple[0]].get()),
                  static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                      m_solvers[solversToCouple[1]].get()));
          break;
        }
        case FV_SYSEQN_NS:
        default: {
          m_couplers[couplerId] = make_unique<
              FvCartesianInterpolation<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>, FvSysEqnNS<nDim>>>(
              couplerId,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                  m_solvers[solversToCouple[0]].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[1]].get()));
          break;
        }
      }
      break;
    }

    case COUPLER_FV_MULTILEVEL: {
      // Note: couple all solvers (FV) in a multilevel computation
      auto fvSolvers = std::vector<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>{};
      for(MInt s = 0; s < (MInt)solversToCouple.size(); s++) {
        fvSolvers.push_back(
            static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[s]].get()));
      }
      m_couplers[couplerId] = make_unique<CouplerFvMultilevel<nDim, FvSysEqnNS<nDim>>>(couplerId, fvSolvers);

      break;
    }

    case COUPLER_FV_MULTILEVEL_INTERPOLATION: {
      // Multilevel interpolation: use coarse grid restart and interpolate it onto the finer grid(s)/solver(s)
      auto fvSolvers = std::vector<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>{};
      for(MInt s = 0; s < (MInt)solversToCouple.size(); s++) {
        fvSolvers.push_back(
            static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[s]].get()));
      }
      m_couplers[couplerId] =
          make_unique<CouplerFvMultilevelInterpolation<nDim, FvSysEqnNS<nDim>>>(couplerId, fvSolvers);
      break;
    }

    case COUPLER_FV_DG_APE: {
      m_couplers[couplerId] = make_unique<DgCcAcousticPerturb<nDim, FvSysEqnNS<nDim>>>(
          couplerId,
          static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[0]].get()),
          static_cast<DgCartesianSolver<nDim, DgSysEqnAcousticPerturb<nDim>>*>(m_solvers[solversToCouple[1]].get()));
      break;
    }
    case COUPLER_LS_LB: {
      MInt nDist = 27;
      nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[1], AT_, &nDist);

      switch(nDist) {
        case 9: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<2>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<2, 9>;
          auto lbSolver = static_cast<LbSolverDxQy<2, 9, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLb<2, 9, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        case 19: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<3>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLb<3, 19, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        case 27: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<3>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLb<3, 27, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions!");
        }
      }

      break;
    }

    case COUPLER_LB_RB: {
      MInt nDist = 27;
      nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[0], AT_, &nDist);

      switch(nDist) {
        case 9: {
          using SysEqnLb = lb::LbSysEqnIncompressible<2, 9>;
          auto lbSolver = static_cast<LbSolverDxQy<2, 9, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          auto rigidBodies = static_cast<RigidBodies<2>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbRb<2, 9, SysEqnLb>>(couplerId, lbSolver, rigidBodies);

          break;
        }
        case 19: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          auto rigidBodies = static_cast<RigidBodies<3>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbRb<3, 19, SysEqnLb>>(couplerId, lbSolver, rigidBodies);

          break;
        }
        case 27: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          auto rigidBodies = static_cast<RigidBodies<3>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbRb<3, 27, SysEqnLb>>(couplerId, lbSolver, rigidBodies);

          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions!");
        }
      }
      break;
    }

    case COUPLER_LB_LPT: {
      MInt nDist = 27;
      nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[1], AT_, &nDist);

      switch(nDist) {
        case 9: {
          IF_CONSTEXPR(nDim == 3) { mTerm(1, "LPT not supported in 2D!"); }
          break;
        }
        case 19: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          auto particleSolver = static_cast<LPT<3>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbLpt<3, 19, SysEqnLb>>(couplerId, particleSolver, lbSolver);

          break;
        }
        case 27: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
          auto particleSolver = static_cast<LPT<3>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LbLpt<3, 27, SysEqnLb>>(couplerId, particleSolver, lbSolver);

          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions!");
        }
      }
      break;
    }

    case COUPLER_LS_LB_SURFACE: {
      MInt nDist = 27;
      nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[1], AT_, &nDist);

      switch(nDist) {
        case 9: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<2>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<2, 9>;
          auto lbSolver = static_cast<LbSolverDxQy<2, 9, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLbSurface<2, 9, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        case 19: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<3>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLbSurface<3, 19, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        case 27: {
          auto lsSolver = static_cast<typename LsCartesianSolverXD<3>::type*>(m_solvers[solversToCouple[0]].get());
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[1]].get());
          m_couplers[couplerId] = make_unique<LsLbSurface<3, 27, SysEqnLb>>(couplerId, lsSolver, lbSolver);

          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions!");
        }
      }

      break;
    }

    case COUPLER_LS_FV: {
      auto lsSolver = static_cast<typename LsCartesianSolverXD<nDim>::type*>(m_solvers[solversToCouple[0]].get());
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          auto fvSolvers = (FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*)
                               m_solvers[solversToCouple[1]]
                                   .get();
          m_couplers[couplerId] =
              make_unique<typename CouplingLsFvXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>::type>(
                  couplerId, lsSolver, fvSolvers);
          break;
        }
        default: {
          auto fvSolvers = (FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*)m_solvers[solversToCouple[1]].get();
          m_couplers[couplerId] =
              make_unique<typename CouplingLsFvXD<nDim, FvSysEqnNS<nDim>>::type>(couplerId, lsSolver, fvSolvers);
        }
      }

      break;
    }

    case COUPLER_LS_FV_COMBUSTION: {
      auto lsSolverId = solversToCouple[0];
      auto fvSolverId = solversToCouple[1];

      m_couplers[couplerId] = make_unique<typename LsFvCombustionXD<nDim, FvSysEqnNS<nDim>>::type>(
          couplerId,
          static_cast<typename LsCartesianSolverXD<nDim>::type*>(m_solvers[lsSolverId].get()),
          static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[fvSolverId].get()));
      break;
    }

    case COUPLER_LB_FV_EE_MULTIPHASE: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      if(nDim != 3) mTerm(1, AT_, "LB-FV-Euler-Euler-Multiphase only implemented for nDim = 3");

      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          mTerm(1, AT_, "RANS not supported in combination with LB-FV-EE-Multiphase!");
          break;
        }
        case FV_SYSEQN_EEGAS: {
          MInt nDist = 27;
          nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[0], AT_, &nDist);

          switch(nDist) {
            case 9: {
              mTerm(1, "nDist = 9 not supported with EE Multiphase!");

              break;
            }
            case 19: {
              using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
              auto lbSolver = static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
              auto fvSolver =
                  static_cast<FvCartesianSolverXD<3, FvSysEqnEEGas<3>>*>(m_solvers[solversToCouple[1]].get());

              m_couplers[couplerId] = make_unique<CouplerLbFvEEMultiphase<3, 19, SysEqnLb, FvSysEqnEEGas<3>>>(
                  couplerId, lbSolver, fvSolver);

              break;
            }
            case 27: {
              using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
              auto lbSolver = static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[0]].get());
              auto fvSolver =
                  static_cast<FvCartesianSolverXD<3, FvSysEqnEEGas<3>>*>(m_solvers[solversToCouple[1]].get());

              m_couplers[couplerId] = make_unique<CouplerLbFvEEMultiphase<3, 27, SysEqnLb, FvSysEqnEEGas<3>>>(
                  couplerId, lbSolver, fvSolver);

              break;
            }
            default: {
              mTerm(1, "Unknown number of distributions!");
            }
          }
          break;
        }
        case FV_SYSEQN_NS:
        default: {
          mTerm(1, AT_, "This SysEqn is not supported in combination with LB-FV-EE-Multiphase!");
          break;
        }
      }
      break;
    }

    case COUPLER_LB_LB: {
      MInt nDist = 27;
      nDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[0], AT_, &nDist);
      {
        const MInt tempNDist = Context::getSolverProperty<MInt>("noDistributions", solversToCouple[1], AT_, &nDist);
        if(nDist != tempNDist) mTerm(1, AT_, "Can't couple LB solvers with different noDistributions!");
      }
      switch(nDist) {
        case 9: {
          using SysEqnLb = lb::LbSysEqnIncompressible<2, 9>;
          std::vector<LbSolverDxQy<2, 9, SysEqnLb>*> lbSolvers;
          lbSolvers.push_back(static_cast<LbSolverDxQy<2, 9, SysEqnLb>*>(m_solvers[solversToCouple[0]].get()));
          lbSolvers.push_back(static_cast<LbSolverDxQy<2, 9, SysEqnLb>*>(m_solvers[solversToCouple[1]].get()));
          m_couplers[couplerId] = make_unique<CouplerLbLb<2, 9, SysEqnLb>>(couplerId, lbSolvers);
          break;
        }
        case 19: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 19>;
          std::vector<LbSolverDxQy<3, 19, SysEqnLb>*> lbSolvers;
          lbSolvers.push_back(static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[0]].get()));
          lbSolvers.push_back(static_cast<LbSolverDxQy<3, 19, SysEqnLb>*>(m_solvers[solversToCouple[1]].get()));
          m_couplers[couplerId] = make_unique<CouplerLbLb<3, 19, SysEqnLb>>(couplerId, lbSolvers);
          break;
        }
        case 27: {
          using SysEqnLb = lb::LbSysEqnIncompressible<3, 27>;
          std::vector<LbSolverDxQy<3, 27, SysEqnLb>*> lbSolvers;
          lbSolvers.push_back(static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[0]].get()));
          lbSolvers.push_back(static_cast<LbSolverDxQy<3, 27, SysEqnLb>*>(m_solvers[solversToCouple[1]].get()));
          m_couplers[couplerId] = make_unique<CouplerLbLb<3, 27, SysEqnLb>>(couplerId, lbSolvers);
          break;
        }
        default: {
          mTerm(1, "Unknown number of distributions!");
        }
      }
      break;
    }

    case COUPLER_FV_PARTICLE: {
      auto particleSolver = static_cast<LPT<nDim>*>(m_solvers[solversToCouple[1]].get());
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", solversToCouple[1])) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", solversToCouple[1], AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          auto fvSolver = static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
              m_solvers[solversToCouple[0]].get());
          if(nDim == 3) {
            m_couplers[couplerId] =
                make_unique<CouplerFvParticle<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>>(
                    couplerId, particleSolver, fvSolver);
          } else {
            mTerm(1, AT_, "FV Particle not supported in 2D!");
          }
          break;
        }
        case FV_SYSEQN_NS:
        default: {
          auto fvSolver =
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[solversToCouple[0]].get());
          if(nDim == 3) {
            m_couplers[couplerId] =
                make_unique<CouplerFvParticle<nDim, FvSysEqnNS<nDim>>>(couplerId, particleSolver, fvSolver);
          } else {
            mTerm(1, AT_, "FV Particle not supported in 2D!");
          }

          break;
        }
      }
      break;
    }

    default: {
      TERMM(1, "Unknown coupler type '" + couplerType + "', exiting ... ");
    }
  }

  m_log << "Created coupler #" << std::to_string(couplerId) << ": " << couplerType << " for solvers";
  for(MInt solver = 0; solver < 2; solver++) {
    m_log << " " << solversToCouple[solver];
  }
  m_log << std::endl;
}


/** \brief This function handels the creation of the execution recipe
 * \author Thomas Hoesgen
 * \date 10/2019
 */
ExecutionRecipe* Application::createRecipe() {
  TRACE();

  ExecutionRecipe* recipe = nullptr;

  const MString recipeType = Context::getBasicProperty<MString>("executionRecipe", AT_);

  switch(string2enum(recipeType)) {
    case RECIPE_BASE: {
      recipe = new ExecutionRecipe(&m_solvers, &m_couplers);
      break;
    }
    case RECIPE_INTRASTEP: {
      recipe = new ExecutionRecipeIntraStepCoupling(&m_solvers, &m_couplers);
      break;
    }
    case RECIPE_ITERATION: {
      recipe = new ExecutionRecipeSolutionIteration(&m_solvers, &m_couplers);
      break;
    }
    default: {
      TERMM(1, "Unknown execution recipe '" + recipeType + "', exiting ... ");
    }
  }

  m_log << "Created execution recipe: " << recipeType << " for " << m_noSolvers << " solvers" << std::endl;

  return recipe;
}

/** \brief This function handels the creation of the Postprocessing classes
 * \author Jannik Borgelt
 * \date 07/2020
 */
template <MInt nDim>
PostProcessingInterface* Application::createPostProcessing(MInt postprocessingId) {
  TRACE();

  PostProcessingInterface* pp = nullptr;

  const MString postprocessingType =
      Context::getBasicProperty<MString>("postProcessingType_" + std::to_string(postprocessingId), AT_);

  // FIXME labels:PP,DOC read with _$postprocessingId and allow for multiple postProcessingSolverIds
  //       to specify the second/third solver Ids in coupled cases!
  const MInt postprocessingSolverId =
      Context::getBasicProperty<MInt>("postProcessingSolverIds", AT_, nullptr, postprocessingId);

  PostData<nDim>* postDataPointer = nullptr;
  // postprocessing for special output that does not need a PostData Solver
  if(m_postDataSolverId == -1) {
    if(globalDomainId() == 0)
      cerr << "\033[0;31m#### WARNING:\033[0m You are using a postprocessing routine without a PostData Solver!"
           << endl;
    m_postDataSolverId = postprocessingSolverId;
  } else {
    postDataPointer = static_cast<PostData<nDim>*>(m_solvers[m_postDataSolverId].get());
  }

  // Create the PostProcessing
  //
  switch(string2enum(postprocessingType)) {
    case POSTPROCESSING_FV: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", postprocessingSolverId)) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", postprocessingSolverId, AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          pp = new PostProcessingFv<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>(
              postprocessingId,
              postDataPointer,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                  m_solvers[postprocessingSolverId].get()));
          break;
        }
        case FV_SYSEQN_EEGAS: {
          IF_CONSTEXPR(nDim == 2)
          mTerm(1, AT_, "FVSysEqnEEGas is not tested for 2D at all and probably does not make much sense!");
          else pp = new PostProcessingFv<nDim, FvSysEqnEEGas<nDim>>(
              postprocessingId,
              postDataPointer,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnEEGas<nDim>>*>(m_solvers[postprocessingSolverId].get()));
          break;
        }
        case FV_SYSEQN_NS: {
          pp = new PostProcessingFv<nDim, FvSysEqnNS<nDim>>(
              postprocessingId,
              postDataPointer,
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[postprocessingSolverId].get()));
          break;
        }
        case FV_SYSEQN_DETCHEM: {
          pp = new PostProcessingFv<nDim, FvSysEqnDetChem<nDim>>(
              postprocessingId,
              static_cast<PostData<nDim>*>(m_solvers[m_postDataSolverId].get()),
              static_cast<FvCartesianSolverXD<nDim, FvSysEqnDetChem<nDim>>*>(m_solvers[postprocessingSolverId].get()));
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case POSTPROCESSING_DG: {
      const MInt dgSystemEquations =
          string2enum(Context::getSolverProperty<MString>("dgSystemEquations", postprocessingSolverId, AT_));
      switch(dgSystemEquations) {
        case DG_SYSEQN_ACOUSTICPERTURB: {
          pp = new PostProcessingDg<nDim, DgSysEqnAcousticPerturb<nDim>>(
              postprocessingId,
              postDataPointer,
              static_cast<DgCartesianSolver<nDim, DgSysEqnAcousticPerturb<nDim>>*>(
                  m_solvers[postprocessingSolverId].get()));
          break;
        }
        case DG_SYSEQN_LINEARSCALARADV: {
          pp = new PostProcessingDg<nDim, DgSysEqnLinearScalarAdv<nDim>>(
              postprocessingId,
              postDataPointer,
              static_cast<DgCartesianSolver<nDim, DgSysEqnLinearScalarAdv<nDim>>*>(
                  m_solvers[postprocessingSolverId].get()));
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case POSTPROCESSING_LB: {
      pp = new PostProcessingLb<nDim>(postprocessingId, postDataPointer,
                                      static_cast<LbSolver<nDim>*>(m_solvers[postprocessingSolverId].get()));
      break;
    }
    case POSTPROCESSING_FVLPT: {
      MInt fvSystemEquations = string2enum("FV_SYSEQN_NS");
      if(Context::propertyExists("fvSystemEquations", postprocessingSolverId)) {
        fvSystemEquations =
            string2enum(Context::getSolverProperty<MString>("fvSystemEquations", postprocessingSolverId, AT_));
      }
      switch(fvSystemEquations) {
        case FV_SYSEQN_RANS: {
          if(nDim == 3) {
            pp = new PostProcessingFvLPT<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>(
                postprocessingId,
                postDataPointer,
                static_cast<FvCartesianSolverXD<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>*>(
                    m_solvers[postprocessingSolverId].get()),
                static_cast<LPT<nDim>*>(m_solvers[postprocessingSolverId + 1].get()));
          } else {
            mTerm(1, AT_, "LPT not supported in 2D!");
          }
          break;
        }
        case FV_SYSEQN_NS: {
          if(nDim == 3) {
            pp = new PostProcessingFvLPT<nDim, FvSysEqnNS<nDim>>(
                postprocessingId,
                postDataPointer,
                static_cast<FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>*>(m_solvers[postprocessingSolverId].get()),
                static_cast<LPT<nDim>*>(m_solvers[postprocessingSolverId + 1].get()));
          } else {
            mTerm(1, AT_, "LPT not supported in 2D!");
          }
          break;
        }
        default:
          TERMM(1, "Unsupported system of equations.");
      }
      break;
    }
    case POSTPROCESSING_LBLPT: {
      if(nDim == 3) {
        pp = new PostProcessingLbLPT<nDim>(postprocessingId,
                                           postDataPointer,
                                           static_cast<LbSolver<nDim>*>(m_solvers[postprocessingSolverId].get()),
                                           static_cast<LPT<nDim>*>(m_solvers[postprocessingSolverId + 1].get()));
      } else {
        mTerm(1, AT_, "LPT not supported in 2D!");
      }
      break;
    }
    default: {
      TERMM(1, "Unknown postprocessing type '" + postprocessingType + "', exiting ... ");
    }
  }

  m_log << "Created postprocessing #" << std::to_string(postprocessingId) << ": " << postprocessingType << " for solver"
        << " " << postprocessingSolverId << std::endl;

  return pp;
}

/// Collect the timings of all solvers and domain decomposition information
void Application::collectTimingsAndSolverInformation(const MBool finalTimeStep) {
  // Return if not enabled
  if(!m_writeSolverTimings) {
    return;
  }

  const MBool isSamplingStep = (globalTimeStep % m_solverTimingsSampleInterval == 0);
  const MInt noSolversAndCouplers = m_noSolvers + m_noCouplers;
  const MBool allTimings = m_writeAllSolverTimings;

  if(isSamplingStep) {
    // Store time step
    m_solverTimingsTimeStep.push_back(globalTimeStep);

    // TODO labels:TIMERS store domainInfo for each sampling time step/for each changed domain decomposition?
    // Determine domain decomposition information of each solver
    m_domainInfo.clear();
    for(MInt solverId = 0; solverId < m_noSolvers; solverId++) {
      m_solvers[solverId]->getDomainDecompositionInformation(m_domainInfo);
    }
    for(MInt i = 0; i < m_noCouplers; i++) {
      m_couplers[i]->getDomainDecompositionInformation(m_domainInfo);
    }
  }

  // TODO labels:TIMERS if a sample interval N > 1 is given, just sample every N-th step (current status ) or
  // compute an average? Note: for computing an average the number of timings might not correspond
  // to the sample-interval
  MInt timerIndex = 0;
  for(MInt j = 0; j < noSolversAndCouplers; j++) {
    const MBool isSolver = (j < m_noSolvers);
    const MInt noTimers = (isSolver) ? m_solvers[j]->noSolverTimers(allTimings)
                                     : m_couplers[j - m_noSolvers]->noCouplingTimers(allTimings);

    // Get timing records of solver/coupler
    std::vector<std::pair<MString, MFloat>> timings{};
    (isSolver) ? m_solvers[j]->getSolverTimings(timings, allTimings)
               : m_couplers[j - m_noSolvers]->getCouplingTimings(timings, allTimings);

    if((MInt)timings.size() != noTimers) {
      TERMM(1, "Wrong number of solver timings returned by getSolverTimings(): is " + std::to_string(timings.size())
                   + ", should be " + std::to_string(noTimers));
    }

    for(MInt i = 0; i < noTimers; i++) {
      const MFloat timing = timings[i].second;
      const MString name = timings[i].first;
      // Compute time difference
      MFloat diff = timing - m_solverTimingsPrevTime[timerIndex];

      if(name != m_solverTimingsNames[timerIndex]) {
        TERMM(1, "Timer name does not match.");
      }

      // This happens if the DLB timers were reset from the gridcontroller
      if(diff < 0.0) {
        diff = timing;
      }

      // Only store timing if it should be written to the timings file
      if(isSamplingStep) {
        m_solverTimings[timerIndex].push_back(diff);
      }

      m_solverTimingsPrevTime[timerIndex] = timing;
      timerIndex++;
    }
  }

  // Store timings if this is the final time step or a timings write time step
  storeTimingsAndSolverInformation(finalTimeStep);
}


/// Store timings for all solvers and domain decomposition information on all domains
void Application::storeTimingsAndSolverInformation(const MBool finalTimeStep) {
  TRACE();
  using namespace maia::parallel_io;

  // Return if not enabled
  if(!m_writeSolverTimings) {
    return;
  }

  // Check if timings should be written at this step
  const MBool writeStep = (m_solverTimingsWriteInterval > 0)
                              ? (globalTimeStep % m_solverTimingsWriteInterval == 0 || finalTimeStep)
                              : finalTimeStep;
  if(!writeStep) {
    return;
  }

  const MInt noValues = m_solverTimings.size();
  if(noValues == 0) {
    return;
  }
  const MInt noTimings = m_solverTimings[0].size();
  if(noTimings == 0) {
    return;
  }

  // Timings output file name
  stringstream fileName;
  fileName << "solverTimings_n" << globalNoDomains() << "_" << globalTimeStep << ParallelIo::fileExt();

  // check for existing file
  if(globalDomainId() == 0) {
    ifstream readFile;
    readFile.open(fileName.str());
    if(readFile) {
      stringstream fileNameNew;
      fileNameNew << "solverTimings_n" << globalNoDomains() << "_" << globalTimeStep << "_bu" << ParallelIo::fileExt();
      std::rename(fileName.str().c_str(), fileNameNew.str().c_str());
    }
  }

  ParallelIo file(fileName.str(), PIO_REPLACE, globalMaiaCommWorld());
  file.defineArray(PIO_INT, "timeStep", noTimings);

  const MInt noInfo = m_domainInfo.size();

  // TODO labels:TIMERS store domainInfo for each changed configuration when using DLB?
  // Domain information
  ParallelIo::size_type dimSizesInfo[] = {globalNoDomains(), noInfo};
  file.defineArray(PIO_INT, "domainInfo", 2, &dimSizesInfo[0]);
  file.setAttribute("domain index", "dim_0", "domainInfo");
  file.setAttribute("information index", "dim_1", "domainInfo");

  for(MInt i = 0; i < noInfo; i++) {
    std::stringstream var;
    var << "var_" << i;
    file.setAttribute(m_domainInfo[i].first, var.str(), "domainInfo");
  }

  // Dimensions of timings: noDomains x noTimings x noValues
  ParallelIo::size_type dimSizesTimings[] = {globalNoDomains(), noTimings, noValues};

  file.defineArray(PIO_FLOAT, "timings", 3, &dimSizesTimings[0]);

  file.setAttribute("domain index", "dim_0", "timings");
  file.setAttribute("time step index", "dim_1", "timings");
  file.setAttribute("timings index", "dim_2", "timings");

  for(MInt i = 0; i < noValues; i++) {
    std::stringstream var;
    var << "var_" << i;
    file.setAttribute(m_solverTimingsNames[i], var.str(), "timings");
  }

  // Assemble domain information
  MFloatScratchSpace info(noInfo, AT_, "data");
  for(MInt i = 0; i < noInfo; i++) {
    info[i] = m_domainInfo[i].second;
  }

  // Assemble timings
  MFloatScratchSpace data(noTimings, noValues, AT_, "data");
  for(MInt i = 0; i < noTimings; i++) {
    for(MInt j = 0; j < noValues; j++) {
      data(i, j) = m_solverTimings[j][i];
    }
  }

  // root writes timestep data
  if(globalDomainId() == 0) {
    file.setOffset(noTimings, 0);
  } else {
    file.setOffset(0, 0);
  }
  file.writeArray(&m_solverTimingsTimeStep[0], "timeStep");

  // Write domain information
  file.setOffset(1, globalDomainId(), 2);
  file.writeArray(&info[0], "domainInfo");

  // Write timings of all domains
  file.setOffset(1, globalDomainId(), 3);
  file.writeArray(&data[0], "timings");


  // Clear collected timings
  for(MInt timerId = 0; timerId < m_noGlobalSolverTimers; timerId++) {
    m_solverTimings[timerId].clear();
    m_solverTimingsTimeStep.clear();
  }
}

/// call cleanup functions in all solvers and couplers before returning from the application
void Application::cleanUp() {
  TRACE();

  // Clean up
  for(MInt i = 0; i < m_noSolvers; i++) {
    m_solvers[i]->cleanUp();
  }
  for(auto&& coupler : m_couplers) {
    coupler->cleanUp();
  }
}

template void Application::run<2>();
template void Application::run<3>();
