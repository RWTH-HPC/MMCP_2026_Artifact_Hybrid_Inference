// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "solver.h"
#include "IO/context.h"

using namespace std;

std::map<MInt, MString> Solver::m_aliases;

//---------------------------------------------------------------------------
//
/*! \fn Solver constructor
 * \brief Data for the solvers will be collected from the property file
 *        Data: nDim, Re, Ma, CFL, angle, gamma
 */
Solver::Solver(const MInt solverId, const MPI_Comm comm, const MBool isActive)
  : m_mpiComm(comm), m_noDim(read_nDim()), m_solverId(solverId) {
  TRACE();

  if(isActive) {
    // Determine domainId and number of domains if this solver is active
    MPI_Comm_rank(mpiComm(), &m_domainId);
    MPI_Comm_size(mpiComm(), &m_noDomains);
  } else {
    m_domainId = -1;
    m_noDomains = -1;
  }

  // Find aliases, default is solverId
  if(Context::solverPropertyExists("solverAlias", solverId)) {
    MString solverAlias = std::to_string(m_solverId);
    solverAlias = Context::getSolverProperty<MString>("solverAlias", m_solverId, AT_, &solverAlias);
    m_aliases[m_solverId] = solverAlias;
  }

  // Get output and restart directories
  MString testcaseDir = "./";
  m_testcaseDir = Context::getSolverProperty<MString>("testcaseDir", m_solverId, AT_, &testcaseDir);
  m_outputDir = Context::getSolverProperty<MString>("outputDir", m_solverId, AT_, &testcaseDir);
  m_restartDir = Context::getSolverProperty<MString>("restartDir", m_solverId, AT_, &m_outputDir);
  m_solutionOutput =
      Context::getSolverProperty<MString>("solutionOutput", m_solverId, AT_, &m_outputDir); // Naming a la FV
  // Context::getSolverProperty<MString>("solutionDir", m_solverId, AT_, &m_outputDir); // TODO: unify?
  m_outputDir = testcaseDir + m_outputDir;
  m_restartDir = testcaseDir + m_restartDir;
  m_solutionOutput = testcaseDir + m_solutionOutput;

  /*! \page propertiesGlobal
    \section restartFile
    <code>MInt Solver::m_restartFile</code>\n
    default = <code>false</code>\n\n
    This property determines the restart.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL</i>
  */
  m_restartFile = false;
  m_restartFile = Context::getSolverProperty<MBool>("restartFile", m_solverId, AT_, &m_restartFile);

  /*! \page propertiesGlobal
   * \section restartFile
  <code>MInt Solver::m_restart</code>\n
  default = <code>false</code>\n\n
  This property determines the restart.
  <ul>
  <li><code>false</code> begin from initial condition</li>
  <li><code>true</code> restart from file</li>
  </ul>\n
  Keywords: <i>GENERAL, GLOBAL</i>
  */
  m_restart = false;
  m_restart = Context::getSolverProperty<MBool>("restartFile", m_solverId, AT_, &m_restart);

  /*! \page propertiesGlobal
    \section residualInterval
    <code>MInt Solver::m_residualInterval</code>\n
    default = <code></code>\n\n
    Controls the interval in which the Residual should be\n
    written out to the Residual file.\n
    <ul>
    <li>Any integer value > 0 </li>
    </ul>\n
    Keywords: <i>FV, SOLVER, RESIDUAL</i>
  */
  m_residualInterval = Context::getSolverProperty<MInt>("residualInterval", m_solverId, AT_);

  /*! \page propertiesGlobal
    \section restartInterval
    <code>MInt Solver::m_restartInterval</code>\n
    default = <code></code>\n\n
    Controls the interval in which a restart file should be written\n
    <ul>
    <li>Any integer value > 0 </li>
    </ul>\n
    Keywords: <i>FV, SOLVER, RESTART</i>
  */
  m_restartInterval = -1;
  m_restartInterval = Context::getSolverProperty<MInt>("restartInterval", m_solverId, AT_, &m_restartInterval);

  m_restartOffset = 0;
  /*! \page propertiesGlobal
    \section restartOffset
    <code>MInt Solver::m_restartOffset/</code>\n
    default = <code>0</code>\n\n
    Sets an offset for writing restart files. Restart files are only written after the offset timestep.
    <ul>
    <li>See enum <code>SolverType</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, SOLVER, SOLVER_TYPE, SOLVERS</i>
  */
  m_restartOffset = Context::getSolverProperty<MInt>("restartOffset", m_solverId, AT_, &m_restartOffset);

  m_useNonSpecifiedRestartFile = false;
  m_useNonSpecifiedRestartFile =
      Context::getSolverProperty<MBool>("useNonSpecifiedRestartFile", m_solverId, AT_, &m_useNonSpecifiedRestartFile);

  /*! \page propertiesGlobal
    \section initFromRestartFile
    <code>MBool LbSolver::m_initFromRestartFile</code>\n
    default = <code>false</code>\n\n
    This property defines if a it should be initialized from a restart file.
    <ul>
    <li><code>true</code> (off)</li>
    <li><code>false</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE BOLTZMANN</i>
    */
  m_initFromRestartFile = false;
  m_initFromRestartFile =
      Context::getSolverProperty<MBool>("initFromRestartFile", m_solverId, AT_, &m_initFromRestartFile);

  if(!m_restart && !m_initFromRestartFile) {
    m_restartTimeStep = 0;
  } else {
    if(!m_useNonSpecifiedRestartFile) {
      m_restartTimeStep = Context::getSolverProperty<MInt>("restartTimeStep", m_solverId, AT_);
    } else {
      m_restartTimeStep = 0;
    }
  }

  m_solutionInterval = Context::getSolverProperty<MInt>("solutionInterval", m_solverId, AT_);

  /*! \page propertiesGlobal
    \section solvertype
    <code>MInt Solver::m_solverType</code>\n
    default = <code>no default set</code>\n\n
    Determines which discretization method is used.
    <ul>
      <li>See enum <code>SolverType</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, SOLVER, SOLVER_TYPE, SOLVERS</i>
  */
  m_solverType = Context::getSolverProperty<MString>("solvertype", m_solverId, AT_);

  MString solverMethod = "defaultMethod";
  m_solverMethod = Context::getSolverProperty<MString>("solverMethod", m_solverId, AT_, &solverMethod);
}

/// Read sampling variables names, store in vector and return the number of sampling variables
MInt Solver::readSolverSamplingVarNames(std::vector<MString>& varNames, const MString featureName) const {
  MInt noVars = 0;
  MString propName = "samplingVariables"; // default property
  if(featureName != ""
     && Context::propertyExists("samplingVariables" + featureName, m_solverId)) { // feature specific property
    propName += featureName;
  }

  if(Context::propertyExists(propName, m_solverId)) {
    // Number of sampling variables
    noVars = Context::propertyLength(propName, solverId());

    varNames.clear();
    for(MInt i = 0; i < noVars; i++) {
      const MString samplingVarName = Context::getSolverProperty<MString>(propName, solverId(), AT_, i);
      varNames.push_back(samplingVarName);
    }
  }
  return noVars;
}

MString Solver::getIdentifier(const MBool useSolverId, const MString preString, const MString postString) {
  if(Context::solverPropertyExists("solverAlias", m_solverId)) {
    return preString + m_aliases.at(m_solverId) + postString;
  } else if(useSolverId) {
    return preString + std::to_string(m_solverId) + postString;
  } else {
    return "";
  }
}
