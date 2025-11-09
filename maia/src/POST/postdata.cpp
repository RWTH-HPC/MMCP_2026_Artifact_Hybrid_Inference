// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postdata.h"
#include <algorithm>
#include <stack>
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "globals.h"

using namespace std;

// TODO labels:PP these variable names only work for FV. This means they probably should not be stored here.
template <>
const std::vector<std::vector<MString>> PostData<3>::m_averageVariableName = {
    {"um", "vm", "wm", "rhom", "pm", "cm"},           // primitive
    {"u'", "v'", "w'", "u'v'", "v'w'", "w'u'", "p'"}, // square
    {"u''", "v''", "w''"},                            // skewness
    {"u'''", "v'''", "w'''"},                         // kurtosis
    {"hm", "c'", "h'"},                               // statisticCombustionAnalysis
    {"vort_x", "vort_y", "vort_z"},                   // averageVorticities
    {"c0", "dc0_dx", "dc0_dy", "dc0_dz"},             // averageSpeedOfSound
    {"wxv_x", "wxv_y", "wxv_z"},                      // lamb0
    {"gradm_u_x", "gradm_v_y", "gradm_w_z"},          // du
    {"gradm_rho_x", "gradm_rho_y", "gradm_rho_z"},    // drho
    {"gradm_p_x", "gradm_p_y", "gradm_p_z"},          // dp
    {"gradm_u_x", "gradm_u_y", "gradm_u_z", "gradm_v_x", "gradm_v_y", "gradm_v_z", "gradm_w_x", "gradm_w_y",
     "gradm_w_z"},                                       // gradu
    {"ugradu_x", "ugradu_y", "ugradu_z"},                // ugradu
    {"ugradrho_x", "ugradrho_y", "ugradrho_z"},          // ugradrho
    {"mgrad_p_rho_x", "mgrad_p_rho_y", "mgrad_p_rho_z"}, // gradrhop
    {"rhodivu_x", "rhodivu_y", "rhodivu_z"},             // rhodivu
    {"correlation_var"}};

template <>
const std::vector<std::vector<MString>> PostData<2>::m_averageVariableName = {
    {"um", "vm", "rhom", "pm", "cm"},                     // primitive
    {"u'", "v'", "u'v'", "p'"},                           // square
    {"u''", "v''"},                                       // skewness
    {"u'''", "v'''"},                                     // kurtosis
    {"hm", "c'", "h'"},                                   // statisticCombustionAnalysis
    {"vort_z"},                                           // averageVorticities
    {"c0", "dc0_dx", "dc0_dy"},                           // averageSpeedOfSound
    {"wxv_x", "wxv_y"},                                   // lamb0
    {"gradm_u_x", "gradm_v_y"},                           // du
    {"gradm_rho_x", "gradm_rho_y"},                       // drho
    {"gradm_p_x", "gradm_p_y"},                           // dp
    {"gradm_u_x", "gradm_u_y", "gradm_v_x", "gradm_v_y"}, // gradu
    {"ugradu_x", "ugradu_y"},                             // ugradu
    {"ugradrho_x", "ugradrho_y"},                         // ugradrho
    {"mgrad_p_rho_x", "mgrad_p_rho_y"},                   // gradrhop
    {"rhodivu_x", "rhodivu_y"},                           // rhodivu
    {"correlation_var"}};

template <MInt nDim>
PostData<nDim>::PostData(MInt solverId_, GridProxy& gridProxy_, Geom& geometry_, const MPI_Comm comm)

  : maia::CartesianSolver<nDim, PostData<nDim>>(solverId_, gridProxy_, comm, true),
    m_geometry(&geometry_),
    m_time(0.0) {
  TRACE();

  readProperties();

  for(MInt i = 0; i < (MInt)m_propertyName.size(); i++) {
    m_variableOffset.push_back(std::make_pair(std::numeric_limits<MInt>::max(), -std::numeric_limits<MInt>::max()));
  }
}

template <MInt nDim>
void PostData<nDim>::readProperties() {
  TRACE();

  /*! \page propertyPage1
    \section outputFormat
    <code>MString PostData::m_outputFormat </code>\n
    default = <code> "NETCDF" </code>\n
    Defines the output file format of the post data solver\n
    Possible values are:
    <ul>
    <li>NETCDF</li>
    </ul>
    Keywords: <i>POST DATA, FILE FORMAT, NETCDF, VTU</i>
  */
  m_outputFormat = "NETCDF";
  m_outputFormat = Context::getSolverProperty<MString>("outputFormat", m_solverId, AT_, &m_outputFormat);

  /*! \page propertyPage1
    \section solutionOffset
    <code>MInt PostSolver::m_solutionOffset</code>\n
    default = <code>0</code>\n \n
    which time step to start writing out solution
    Possible values are:
    <ul>
    <li> Int </li>
    </ul>
    Keywords: <i>output</i>
  */
  m_solutionOffset = 0;
  m_solutionOffset = Context::getSolverProperty<MInt>("solutionOffset", m_solverId, AT_, &m_solutionOffset);

  m_noVariables = 0;
  m_noVariables = Context::getSolverProperty<MInt>("noVariables", m_solverId, AT_, &m_noVariables);
}

template <MInt nDim>
void PostData<nDim>::copyGridProperties() {
  TRACE();

  // d) Set Properties for valid-grid-cells
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt c = 0; c < noHaloCells(i); c++) {
      a_isHalo(haloCellId(i, c)) = true;
    }
    for(MInt j = 0; j < noWindowCells(i); j++) {
      a_isWindow(windowCellId(i, j)) = true;
    }
  }

  // Inactive can only have halo cells, however noNeighborDomains=-1
  if(!isActive()) {
    for(MInt c = 0; c < a_noCells(); c++) {
      a_isHalo(c) = true;
    }
  }
}


template <MInt nDim>
void PostData<nDim>::initSolver() {
  TRACE();

  if(!isActive()) {
    m_cells.clear();
    m_cells.reset(grid().maxNoCells());
    this->setHaloCellsOnInactiveRanks();
    return;
  }

  setAndAllocateSolverData(true);
}


template <MInt nDim>
void PostData<nDim>::setAndAllocateSolverData(const MBool fullReset) {
  TRACE();

  const MInt maxNoCells = maxNoGridCells();

  mAlloc(m_recalcIds, maxNoCells, "m_recalcIds", -1, AT_);
  for(MInt i = 0; i < maxNoCells; i++) {
    m_recalcIds[i] = i;
  }

  m_cells.setNoVariables(noVariables());

  m_cells.reset(grid().raw().treeb().capacity());
  m_cells.append(c_noCells());

  copyGridProperties();

  // NOTE: setting default here
  setVariableNames();

  if(fullReset) {
    for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
      for(MInt v = 0; v < noVariables(); v++) {
        a_variable(cellId, v) = 0;
      }
    }
  }
}


template <MInt nDim>
void PostData<nDim>::setVariableNames() {
  TRACE();

  m_variablesName.clear();

  for(MInt i = 0; i < noVariables(); i++) {
    MString s = getVariableName(i);
    auto it = std::find(m_variablesName.begin(), m_variablesName.end(), s);
    if(it == m_variablesName.end()) {
      m_variablesName.push_back(getVariableName(i));
    }
  }
}


template <MInt nDim>
MString PostData<nDim>::getVariableName(MInt offset) {
  MString s = "var" + std::to_string(offset);
  if(offset >= getPropertyVariableOffset("primitive").first && offset < getPropertyVariableOffset("primitive").second) {
    MInt index = offset - getPropertyVariableOffset("primitive").first;
    s = m_averageVariableName[0][index];
  }
  if(offset >= getPropertyVariableOffset("square").first && offset < getPropertyVariableOffset("square").second) {
    MInt index = offset - getPropertyVariableOffset("square").first;
    s = m_averageVariableName[1][index];
  }
  if(offset >= getPropertyVariableOffset("skewness").first && offset < getPropertyVariableOffset("skewness").second) {
    MInt index = offset - getPropertyVariableOffset("skewness").first;
    s = m_averageVariableName[2][index];
  }
  if(offset >= getPropertyVariableOffset("kurtosis").first && offset < getPropertyVariableOffset("kurtosis").second) {
    MInt index = offset - getPropertyVariableOffset("kurtosis").first;
    s = m_averageVariableName[3][index];
  }
  if(offset >= getPropertyVariableOffset("statisticCombustionAnalysis").first
     && offset < getPropertyVariableOffset("statisticCombustionAnalysis").second) {
    MInt index = offset - getPropertyVariableOffset("statisticCombustionAnalysis").first;
    s = m_averageVariableName[4][index];
  }
  if(offset >= getPropertyVariableOffset("averageVorticity").first
     && offset < getPropertyVariableOffset("averageVorticity").second) {
    MInt index = offset - getPropertyVariableOffset("averageVorticity").first;
    s = m_averageVariableName[5][index];
  }
  if(offset >= getPropertyVariableOffset("averageSpeedOfSound").first
     && offset < getPropertyVariableOffset("averageSpeedOfSound").second) {
    MInt index = offset - getPropertyVariableOffset("averageSpeedOfSound").first;
    s = m_averageVariableName[6][index];
  }
  if(offset >= getPropertyVariableOffset("lamb0").first && offset < getPropertyVariableOffset("lamb0").second) {
    MInt index = offset - getPropertyVariableOffset("lamb0").first;
    s = m_averageVariableName[7][index];
  }
  if(offset >= getPropertyVariableOffset("du").first && offset < getPropertyVariableOffset("du").second) {
    MInt index = offset - getPropertyVariableOffset("du").first;
    s = m_averageVariableName[8][index];
  }
  if(offset >= getPropertyVariableOffset("drho").first && offset < getPropertyVariableOffset("drho").second) {
    MInt index = offset - getPropertyVariableOffset("drho").first;
    s = m_averageVariableName[9][index];
  }
  if(offset >= getPropertyVariableOffset("dp").first && offset < getPropertyVariableOffset("dp").second) {
    MInt index = offset - getPropertyVariableOffset("dp").first;
    s = m_averageVariableName[10][index];
  }
  if(offset >= getPropertyVariableOffset("gradu").first && offset < getPropertyVariableOffset("gradu").second) {
    MInt index = offset - getPropertyVariableOffset("gradu").first;
    s = m_averageVariableName[11][index];
  }
  if(offset >= getPropertyVariableOffset("ugradu").first && offset < getPropertyVariableOffset("ugradu").second) {
    MInt index = offset - getPropertyVariableOffset("ugradu").first;
    s = m_averageVariableName[12][index];
  }
  if(offset >= getPropertyVariableOffset("ugradrho").first && offset < getPropertyVariableOffset("ugradrho").second) {
    MInt index = offset - getPropertyVariableOffset("ugradrho").first;
    s = m_averageVariableName[13][index];
  }
  if(offset >= getPropertyVariableOffset("gradprho").first && offset < getPropertyVariableOffset("gradprho").second) {
    MInt index = offset - getPropertyVariableOffset("gradprho").first;
    s = m_averageVariableName[14][index];
  }
  if(offset >= getPropertyVariableOffset("rhodivu").first && offset < getPropertyVariableOffset("rhodivu").second) {
    MInt index = offset - getPropertyVariableOffset("rhodivu").first;
    s = m_averageVariableName[15][index];
  }
  if(offset >= getPropertyVariableOffset("correlation").first
     && offset < getPropertyVariableOffset("correlation").second) {
    MInt index = offset - getPropertyVariableOffset("correlation").first;
    s = m_averageVariableName[16][0 /*index*/] + to_string(index);
  }
  return s;
}


template <MInt nDim>
void PostData<nDim>::finalizeInitSolver() {
  TRACE();

  if(!isActive()) return;

  if(m_restart) {
    // TODO labels:PP allow to skip loading restart, e.g. if we want to start the averaging after adding the pp-solver
    // at some point and just need everything set to zero instead?
    loadRestartFile();
  }
}


template <MInt nDim>
void PostData<nDim>::writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString gridFileName,
                                      MInt* recalcIdTree) {
  TRACE();

  // post-data has no cell-data
  if(noVariables() == 0) return;

  m_currentGridFileName = gridFileName;

  if(m_recalcIds != nullptr) {
    for(MInt cellId = 0; cellId < maxNoGridCells(); cellId++) {
      m_recalcIds[cellId] = recalcIdTree[cellId];
    }
  }

  if(writeRestart) {
    saveRestartFile(writeBackup);
  }
}


template <MInt nDim>
MBool PostData<nDim>::prepareRestart(MBool writeRestart, MBool& writeGridRestart) {
  TRACE();

  writeGridRestart = false;

  // write intermediate restart file
  if(m_restartInterval == -1 && !writeRestart) {
    writeRestart = false;
  }

  MInt relativeTimeStep = globalTimeStep - m_restartOffset;
  MInt relativeRestartTimeStep = m_restartTimeStep - m_restartOffset;

  if(((relativeTimeStep % m_restartInterval) == 0 && relativeTimeStep > relativeRestartTimeStep) || writeRestart) {
    writeRestart = true;
  }

  if(m_forceWriteRestart) {
    writeRestart = true;
    if(m_adaptationSinceLastRestart) {
      writeGridRestart = true;
    }
  }

  return writeRestart;
}

/**
 * \brief This function resets the grid-trigger after a restart that is handled by the grid-controller!
 * \author Tim Wegmann
 */

template <MInt nDim>
void PostData<nDim>::reIntAfterRestart(MBool doneRestart) {
  TRACE();

  if(doneRestart) {
    m_adaptationSinceLastRestart = false;
  }

  m_forceWriteRestart = false;
}

template <MInt nDim>
void PostData<nDim>::saveRestartFile(const MBool writeBackup) {
  TRACE();

  ASSERT(noVariables() > 0, "No variables set for PostData restart!");
  ASSERT(!m_variablesName.empty(), "");
  ASSERT((MInt)m_variablesName.size() == noVariables(), "");

  // TODO labels:IO use same restart file name format for all solvers, e.g. restart_s[solverId]_[timeStep].[ext] with
  // zero padded formatted timeStep
  stringstream fileName;
  fileName << outputDir() << "restartFile_" << solverId() << "_" << globalTimeStep << ParallelIo::fileExt();

  saveDataFile(writeBackup, fileName.str(), noVariables(), m_variablesName, &a_variable(0, 0));
}

template <MInt nDim>
void PostData<nDim>::saveDataFile(const MBool writeBackup, const MString fileName, const MInt noVars,
                                  std::vector<MString>& variablesName, MFloat* variables) {
  TRACE();

  ASSERT(noVars > 0, "No variables set for PostData restart!");
  ASSERT((MInt)variablesName.size() == noVars, "");

  if(domainId() == 0) {
    cerr << "Writing post data file '" << fileName << "' for solver " << solverId() << " at time step "
         << globalTimeStep << " ...";
  }

  vector<MInt> reOrderedCells;
  MInt countInternal = 0;
  if(grid().newMinLevel() > 0) {
    if(domainId() == 0) {
      cerr << "Increasing minLevel for solver " << m_solverId;
    }
    this->reOrderCellIds(reOrderedCells);
    for(MUint i = 0; i < reOrderedCells.size(); i++) {
      if(a_isHalo(reOrderedCells[i])) continue;
      countInternal++;
    }
    m_recalcIds = nullptr;
  }
  const MInt noCells = grid().newMinLevel() < 0 ? noInternalCells() : reOrderedCells.size();
  const MInt noInternalCellIds = grid().newMinLevel() < 0 ? noInternalCells() : countInternal;

  MFloatScratchSpace dbVariables(noCells * noVars, AT_, "dbVariables");
  MIntScratchSpace idVariables(noCells * 2, AT_, "idVariables");
  MFloatScratchSpace dbParameters(0, AT_, "dbParameters");
  MIntScratchSpace idParameters(1, AT_, "idParameters");
  vector<MString> dbVariablesName;
  vector<MString> idVariablesName;
  vector<MString> dbParametersName;
  vector<MString> idParametersName;

  this->collectVariables(variables, dbVariables, variablesName, dbVariablesName, noVars, noCells);

  this->collectParameters(globalTimeStep, idParameters, "globalTimeStep", idParametersName);

  MIntScratchSpace recalcIdsSolver(grid().tree().size(), AT_, "recalcIds");
  if(m_recalcIds != nullptr) {
    MInt cellId = 0;
    for(MInt gridcell = 0; gridcell < grid().raw().m_noInternalCells; gridcell++) {
      if(grid().raw().a_hasProperty(gridcell, Cell::IsHalo)) continue;
      MInt l_solverId = grid().tree().grid2solver(m_recalcIds[gridcell]);
      if(l_solverId > -1) {
        recalcIdsSolver[cellId] = l_solverId;
        cellId++;
        ASSERT(grid().solverFlag(m_recalcIds[gridcell], m_solverId), "");
      }
    }
    ASSERT(cellId == grid().noInternalCells(), "recalc ids size is wrong");
  }

  MString gridFile = (m_currentGridFileName == "") ? grid().gridInputFileName() : m_currentGridFileName;

  saveGridFlowVars(fileName.c_str(), gridFile.c_str(), noCells, noInternalCellIds, dbVariables, dbVariablesName, 0,
                   idVariables, idVariablesName, 0, dbParameters, dbParametersName, idParameters, idParametersName,
                   recalcIdsSolver.begin(), -1);

  // TODO labels:PP this is still specific to restart files
  if(writeBackup) {
    stringstream backupFileName;
    backupFileName.clear();
    backupFileName.str("");
    backupFileName << outputDir() << "restartFileBackup_" << solverId() << "_" << globalTimeStep;
    backupFileName << ParallelIo::fileExt();
    if(domainId() == 0) cerr << "Writing post data (backup) for solver " << solverId() << "... ";

    saveGridFlowVars((backupFileName.str()).c_str(), gridFile.c_str(), noCells, noInternalCellIds, dbVariables,
                     dbVariablesName, 0, idVariables, idVariablesName, 0, dbParameters, dbParametersName, idParameters,
                     idParametersName, recalcIdsSolver.begin(), -1);
  }


  if(domainId() == 0) cerr << "ok" << endl;
}


template <MInt nDim>
void PostData<nDim>::loadRestartFile() {
  TRACE();

  stringstream fileName;
  fileName << restartDir() << "restartFile_" << solverId() << "_" << globalTimeStep << ParallelIo::fileExt();

  if(domainId() == 0) {
    cerr << "loading post data for solver " << solverId() << " at time step " << globalTimeStep << " ...";
  }

  vector<MString> name;
  for(MInt v = 0; v < noVariables(); v++) {
    name.push_back(m_variablesName[v]);
  }

  m_log << "loading postprocessing variables ... ";

  loadGridFlowVars((fileName.str()).c_str(), noVariables(), name);

  m_log << "ok" << endl;
}


/** \brief load a file for averaging
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * checks for attribute "isMeanFile", if found all variables are loaded (e.g. means and statistical moments).
 *
 * \param[in] fileName name of data file to load
 *
 **/
template <MInt nDim>
void PostData<nDim>::loadMeanFile(const MString fileName) {
  TRACE();

  if(!isActive()) return;

  ParallelIo parallelIo(fileName, maia::parallel_io::PIO_READ, mpiComm());

  if(parallelIo.hasAttribute(MString("isMeanFile"))) {
    // This should be the same for all the variables
    ParallelIo::size_type dimLen = noInternalCells();
    ParallelIo::size_type start = domainOffset(domainId()) - grid().bitOffset();

    // set offset for all read operations
    parallelIo.setOffset(dimLen, start);

    vector<MString> var_names = parallelIo.getDatasetNames(1);
    m_fileNoVars = var_names.size();
    const MInt noCells = noInternalCells();

    // Note: allocate new storage for all variables from the file
    mAlloc(m_averagedVars, noCells, m_fileNoVars, "m_averagedVars", AT_);
    ScratchSpace<MFloat> tmpVars(noCells, AT_, "tmpVars");

    m_fileVarNames.clear();
    for(MInt varId = 0; varId < m_fileNoVars; varId++) {
      // Check for variable 'name' attribute
      if(parallelIo.hasAttribute("name", var_names[varId])) {
        MString varName = "";
        parallelIo.getAttribute(&varName, "name", var_names[varId]);
        m_fileVarNames.push_back(varName);
      } else {
        TERMM(1, "Variable has no attribute 'name'.");
      }

      parallelIo.readArray(tmpVars.begin(), var_names[varId]);
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        m_averagedVars[cellId][varId] = tmpVars[cellId];
      }
    }
    TERMM_IF_COND(m_fileVarNames.size() > 0 && (MInt)m_fileVarNames.size() != m_fileNoVars,
                  "Error: not every variable in '" + MString(fileName) + "' has the 'name' attribute.");

    m_isMeanFile = true;
  } else {
    TERMM(1, "FIXME try to load regular data file");
    m_isMeanFile = false;
    // loadGridFlowVars(...);
  }
}


/**
 * \brief This function reads the parallel Netcdf cartesian grid cell based solution/restart file
 *        currently used in postData!
 * \author Jannik Borgelt
 */
template <MInt nDim>
void PostData<nDim>::loadGridFlowVars(const MChar* fileName, MInt noVariables, vector<MString> names) {
  TRACE();

  // File loading.
  stringstream variables;

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  // This should be the same for all the variables
  ParallelIo::size_type dimLen = noInternalCells();
  ParallelIo::size_type start = domainOffset(domainId()) - grid().bitOffset();

  // set offset for all read operations
  parallelIo.setOffset(dimLen, start);

  for(MInt vId = 0; vId < noVariables; vId++) {
    MString varName = names[vId]; // varNames[vId];//

    // Load our variables
    MFloatScratchSpace tmpVar((MInt)dimLen, AT_, "tmpVar");

    parallelIo.readArray(tmpVar.getPointer(), "variables" + to_string(vId));
    for(MInt i = 0; i < (MInt)dimLen; ++i) {
      a_variable(i, vId) = tmpVar.p[i];
    }
  }
}

template <MInt nDim>
void PostData<nDim>::saveSolverSolution(const MBool, const MBool) {
  TRACE();

  if(!isActive()) return;

  MInt offset = m_solutionOffset;

  // solution output
  // ---------------
  if((((globalTimeStep - offset) % m_solutionInterval) == 0 && globalTimeStep >= offset)
     || (m_solutionTimeSteps.count(globalTimeStep) > 0)) {
    if(domainId() == 0)
      cerr << "solverId: " << m_solverId << ", Writing " << m_outputFormat << " output at time step " << globalTimeStep
           << ", time " << m_time << " ... ";

    m_log << "Writing " << m_outputFormat << " output at time step " << globalTimeStep << ", time " << m_time
          << " for postData ... ";

    MFloatScratchSpace dbVariables(a_noCells() * noVariables(), AT_, "dbVariables");
    MIntScratchSpace idVariables(a_noCells() * 2, AT_, "idVariables");
    MFloatScratchSpace dbParameters(4, AT_, "dbParameters");
    MIntScratchSpace idParameters(2, AT_, "idParameters");
    vector<MString> dbVariablesName;
    vector<MString> idVariablesName;
    vector<MString> dbParametersName;
    vector<MString> idParametersName;
    vector<MString> name;

    stringstream fileName;

    fileName << m_solutionOutput << "ppQOUT_" << globalTimeStep;

    for(MInt v = 0; v < noVariables(); v++) {
      name.push_back(m_variablesName[v]);
    }

    this->collectVariables(&a_variable(0, 0), dbVariables, name, dbVariablesName, noVariables(), a_noCells());

    this->collectParameters(globalTimeStep, idParameters, "globalTimeStep", idParametersName);

    switch(string2enum(m_outputFormat)) {
      case NETCDF: {
        fileName << ParallelIo::fileExt();
        // TODO: currently not working for adaptation, check save restartFile
        //      and get solver recalcIds!
        //
        saveGridFlowVars((fileName.str()).c_str(), grid().gridInputFileName().c_str(), a_noCells(), noInternalCells(),
                         dbVariables, dbVariablesName, 0, idVariables, idVariablesName, 0, dbParameters,
                         dbParametersName, idParameters, idParametersName, m_recalcIds, -1);
        break;
      }
      case VTK:
      case VTU:
      case VTP: {
        break;
      }
      default: {
        stringstream errorMessage;
        errorMessage << "PostData::saveSolverSolution(): switch variable 'm_outputFormat' with value " << m_outputFormat
                     << " not matching any case." << endl;
        mTerm(1, AT_, errorMessage.str());
      }
    }

    if(domainId() == 0) cerr << "ok" << endl;
  }
}

/**
 * \brief  Get solver timings
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings, const MBool) {
  TRACE();

  const MString namePrefix = "s" + std::to_string(solverId()) + "_";

  const MFloat load = returnLoadRecord();
  const MFloat idle = returnIdleRecord();

  solverTimings.emplace_back(namePrefix + "loadPostData", load);
  solverTimings.emplace_back(namePrefix + "idlePostData", idle);
}

/**
 * \brief  Return the number of Ls load types.
 *
 * Type 1: number of cells
 * Type 2: ??? number cells used for probing
 * Type 3: ??? number cells used for averaging
 *
 * \author Tim Wegmann
 */
template <MInt nDim>
MInt PostData<nDim>::noLoadTypes() const {
  TRACE();

  const MInt noLsLoadTypes = 1;

  return noLsLoadTypes;
}

/**
 * \brief  Return the default weights for all load quantities
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::getDefaultWeights(MFloat* weights, std::vector<MString>& names) const {
  TRACE();

  // TODO set sensible default values
  weights[0] = 1.0;
  names[0] = "pp_cell";
  MInt count = 1;

  // TODO: add other load types!

  if(noLoadTypes() != count) {
    TERMM(1, "Count does not match noLoadTypes.");
  }
}

/**
 * \brief  Return decomposition information, i.e. number of local elements,...
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) {
  TRACE();

  const MString namePrefix = "s" + std::to_string(solverId()) + "_";

  // Number of Post-Data-Cells
  const MInt noCells = c_noCells();

  domainInfo.emplace_back(namePrefix + "noPostCells", noCells);

  // additionally add cells for probing
}

/**
 * \brief sets the cell-weight for balancing and a restarting
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt offset = noCellsGrid * solverId();

  for(MInt cellId = 0; cellId < c_noCells(); cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    const MInt id = gridCellId + offset;
    // TODO: find good value for postData!
    solverCellWeight[id] = 0.2;
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
void PostData<nDim>::getLoadQuantities(MInt* const loadQuantities) const {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // reset
  for(MInt type = 0; type < noLoadTypes(); type++) {
    loadQuantities[type] = 0;
  }

  loadQuantities[0] = c_noCells();

  // TODO: add other load types!
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
MFloat PostData<nDim>::getCellLoad(const MInt gridCellId, const MFloat* const weights) const {
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

  cellLoad = weights[0];

  // TODO: add other load types!

  return cellLoad;
}

/**
 * \brief  prepare solver adaptation
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::prepareAdaptation() {
  TRACE();

  ASSERT(m_freeIndices.empty(), "");

  m_adaptationLevel = maxUniformRefinementLevel();

  if(!isActive()) return;
}

/**
 * \brief set the sensors for a single adaptation step
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::setSensors(std::vector<std::vector<MFloat>>& sensors,
                                std::vector<MFloat>& sensorWeight,
                                std::vector<std::bitset<64>>& sensorCellFlag,
                                std::vector<MInt>& sensorSolverId) {
  TRACE();

  const MInt sensorOffset = (signed)sensors.size();
  ASSERT(sensorOffset == 0 || grid().raw().treeb().noSolvers() > 1, "");
  sensors.resize(sensorOffset + this->m_noSensors, vector<MFloat>(grid().raw().m_noInternalCells, F0));
  sensorWeight.resize(sensorOffset + this->m_noSensors, -1);
  sensorCellFlag.resize(grid().raw().m_noInternalCells, sensorOffset + this->m_noSensors);
  sensorSolverId.resize(sensorOffset + this->m_noSensors, solverId());
  ASSERT(sensorOffset + this->m_noSensors < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

  for(MInt sen = 0; sen < this->m_noSensors; sen++) {
    sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
  }

  ASSERT(m_freeIndices.empty(), "");
  m_freeIndices.clear();

  if(!isActive()) return;
  if(!this->m_adapts) return;

  if(domainId() == 0) {
    cerr << "Setting " << this->m_noSensors << " sensors for post-Data adaptation." << endl;
  }

  for(MInt sen = 0; sen < this->m_noSensors; sen++) {
    (this->*(this->m_sensorFnPtr[sen]))(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen);
  }
}


/**
 * \brief  reinit the solver after a single adaptation step
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::postAdaptation() {
  TRACE();

  if(isActive()) {
    this->compactCells();
    m_freeIndices.clear();
  } else {
    grid().updateGridMap();
  }

  grid().updateOther();

  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);
  this->checkNoHaloLayers();

  // Nothing further to be done if solver inactive
  if(!isActive()) return;

  copyGridProperties();

  m_adaptationLevel++;
}

/**
 * \brief  reinit the solver after the full adaptation loop!
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::finalizeAdaptation() {
  TRACE();

  m_forceAdaptation = false;
  m_adaptationSinceLastRestart = true;
}

/**
 * \brief
 * \author Tim Wegmann, Jannik Borgelt
 */
template <MInt nDim>
void PostData<nDim>::refineCell(const MInt gridCellId) {
  TRACE();

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
    const MInt childId = grid().raw().treeb().child(gridCellId, child);
    if(childId == -1) continue;

    if(!grid().raw().a_hasProperty(childId, Cell::WasNewlyCreated)
       && grid().raw().a_hasProperty(gridCellId, Cell::IsPartLvlAncestor)) {
      continue;
    }

    if(!g_multiSolverGrid) ASSERT(grid().raw().a_hasProperty(childId, Cell::WasNewlyCreated), "");

    // If solver is inactive all cells musst be halo cells!
    if(!isActive()) ASSERT(grid().raw().a_isHalo(childId), "");
    // If child exists in grid but is not located inside solver geometry
    if(!grid().solverFlag(childId, solverId())) continue;

    const MInt solverChildId = this->createCellId(childId);

    if(!g_multiSolverGrid) ASSERT(solverChildId == childId, "");

    for(MInt v = 0; v < noVariables(); v++) {
      a_variable(solverChildId, v) = a_variable(solverCellId, v);
    }

    // check-child-values
#if !defined NDEBUG
    for(MInt v = 0; v < noVariables(); v++) {
      if(std::isnan(a_variable(solverChildId, v))) {
        cerr << "Invalid-value in refined-Cell! "
             << " " << solverChildId << " in rank " << domainId() << endl;
      }
    }
#endif
  }
}


template <MInt nDim>
void PostData<nDim>::removeChilds(const MInt gridCellId) {
  TRACE();
  // If solver is inactive cell must never be an internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(solverCellId > -1 && solverCellId < m_cells.size(), "solverCellId is: " << solverCellId);

  if(!g_multiSolverGrid) ASSERT(solverCellId == gridCellId, "");

  for(MInt c = 0; c < grid().m_maxNoChilds; c++) {
    MInt childId = c_childId(solverCellId, c);
    if(childId < 0) continue;
    this->removeCellId(childId);
  }

  if(!g_multiSolverGrid) {
    ASSERT((grid().raw().treeb().size() - m_cells.size()) <= grid().m_maxNoChilds, "");
  }
}


/**
 * \brief
 */
template <MInt nDim>
void PostData<nDim>::removeCell(const MInt gridCellId) {
  TRACE();
  // If solver is inactive cell musst never be a internal cell
  if(!isActive()) {
    ASSERT(grid().raw().a_isHalo(gridCellId), "");
  }

  const MInt solverCellId = grid().tree().grid2solver(gridCellId);

  ASSERT(gridCellId > -1 && gridCellId < grid().raw().treeb().size() && solverCellId > -1
             && solverCellId < m_cells.size() && grid().tree().solver2grid(solverCellId) == gridCellId,
         "");

  this->removeCellId(solverCellId);
}


template <MInt nDim>
void PostData<nDim>::resizeGridMap() {
  grid().resizeGridMap(m_cells.size());
}


/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::swapCells(const MInt cellId0, const MInt cellId1) {
  TRACE();

  const MInt size = m_cells.size();
  m_cells.append();
  m_cells.erase(size);
  m_cells.copy(cellId0, size);
  m_cells.copy(cellId1, cellId0);
  m_cells.copy(size, cellId1);
  m_cells.erase(size);
  m_cells.size(size);
}


/**
 * \brief
 * \author Tim Wegmann
 */
template <MInt nDim>
void PostData<nDim>::swapProxy(const MInt cellId0, const MInt cellId1) {
  grid().swapGridIds(cellId0, cellId1);
}


/**
 * \brief  checks if a child lies outSide of the domain!
 *         necessary for refinement at the bndry!
 * \author Tim Wegmann
 */

template <MInt nDim>
MInt PostData<nDim>::cellOutside(const MFloat* coords, const MInt level, const MInt gridCellId) {
  return -1;

  std::ignore = coords;
  std::ignore = level;
  std::ignore = gridCellId;
}


/// \brief Reset the solver prior to load balancing
/// \author Thomas Hoesgen
template <MInt nDim>
void PostData<nDim>::resetSolver() {
  TRACE();
  mDeallocate(m_recalcIds);
}

/// \brief Reinitialize solver for DLB prior to setting solution data.
///
/// \author Thomas Hoesgen
template <MInt nDim>
void PostData<nDim>::balancePre() {
  TRACE();

  // Note: every function that allocates persistent memory should first deallocate that memory, for
  // this to work initialize all pointers with nullptr in the class definition

  // Set reinitialization stage
  m_loadBalancingReinitStage = 0;

  // Store currently used memory
  /* const MLong previouslyAllocated = allocatedBytes(); */

  // Update the grid proxy for this solver
  grid().update();

  if(!grid().isActive()) {
    // Reset parallelization information if solver is not active
    updateDomainInfo(-1, -1, MPI_COMM_NULL, AT_);
  } else {
    // Set new domain info for solver
    updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);
  }

  // Return if solver is not active
  if(!grid().isActive()) {
    m_cells.reset(grid().maxNoCells());
    this->setHaloCellsOnInactiveRanks();
    return;
  }

  // Resize cell collector to internal cells
  m_cells.reset(grid().raw().treeb().capacity());

  setAndAllocateSolverData(false);
  this->checkNoHaloLayers();
}

/// \brief Reinitialize solver for DLB after to setting solution data.
///
/// \author Thomas Hoesgen
template <MInt nDim>
void PostData<nDim>::balancePost() {
  TRACE();

  m_loadBalancingReinitStage = 1;
  m_adaptationSinceLastRestart = true;

  // Nothing to do if solver is not active
  if(!grid().isActive()) return;

  m_loadBalancingReinitStage = 2;
}

/// \brief Reinitialize solver after all data structures have been recreated
/// \author Thomas Hoesgen
template <MInt nDim>
void PostData<nDim>::finalizeBalance() {
  TRACE();

  // Nothing to do if solver is not active
  if(!grid().isActive()) return;

  m_loadBalancingReinitStage = -1;
}

/// \brief Return data size to be communicated during DLB for a grid cell and given data id
template <MInt nDim>
MInt PostData<nDim>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  // Inactive ranks do not have any data to communicate
  if(!isActive()) {
    return 0;
  }

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0 || cellId >= noInternalCells()) {
    return 0;
  }

  MInt dataSize = 0;

  switch(dataId) {
    case 0: // variables
      dataSize = noVariables();
      break;
    default:
      TERMM(1, "Unknown data id.");
      break;
  }

  return dataSize;
}

template class PostData<2>;
template class PostData<3>;
