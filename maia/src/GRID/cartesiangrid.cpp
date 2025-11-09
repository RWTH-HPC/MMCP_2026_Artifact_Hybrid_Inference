// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangrid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stack>
#include <sys/stat.h>
#include <unordered_map>
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "UTIL/hilbert.h"
#include "UTIL/maiamath.h"
#include "cartesiangridio.h"
#include "globals.h"
#include "partition.h"
#include "typetraits.h"

using namespace std;
using namespace maia;


template <MInt nDim>
CartesianGrid<nDim>::CartesianGrid(MInt maxCells, const MFloat* const bBox, const MPI_Comm comm,
                                   const MString& fileName)
  : m_maxNoCells(maxCells), m_mpiComm(comm), m_paraViewPlugin((MBool)!fileName.empty()), m_maxNoNghbrs(m_noDirs) {
  TRACE();

  if(paraViewPlugin()) {
    cerr << "Initialising MAIA for the ParaView-Plugin" << endl;
  }


  // Initialize bounding box (only makes sense if multi-solver is disabled)
  if(!g_multiSolverGrid || m_paraViewPlugin) {
    copy_n(bBox, 2 * nDim, &m_boundingBox[0]);
    // Remove the following line once transition to new grid format is complete
    // (see also setupWindowHaloCellConnectivity() -> moved to gridio::loadGrid())
    copy_n(bBox, 2 * nDim, &m_boundingBoxBackup[0]);
  }

  // Determine domainId and number of domains
  MPI_Comm_rank(mpiComm(), &m_domainId);
  MPI_Comm_size(mpiComm(), &m_noDomains);

  const MLong oldAllocatedBytes = allocatedBytes();
  NEW_TIMER_GROUP(gridTimer, "Cartesian Grid");
  NEW_TIMER(timertotal, "Total", gridTimer);
  NEW_SUB_TIMER(timerAllocProp, "allocProperties", timertotal);
  NEW_SUB_TIMER(timerLoadGridPar, "loadGrid", timertotal);
  NEW_SUB_TIMER(timerCreateMB, "window/halo generation", timertotal);

  RECORD_TIMER_START(timertotal);
  RECORD_TIMER_START(timerAllocProp);

  if(!m_paraViewPlugin) {
    m_partitionCellOffspringThreshold = 50000;
    /*! \page propertyPage1
      \section partitionCellMaxNoOffspring
      <code>MInt CartesianGrid::m_partitionCellOffspringThreshold </code>\n
      default = <code>50000</code>\n \n
      controls the filtering of cells in massive parallel file writing\n
      possible values are:
      <ul>
      <li>Non-negative int value of the order less than number of cells on each domains</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, PARALLEL, OUTPUT, FILTER, MAX, MINCELL, SIZE</i>
    */
    if(Context::propertyExists("minCellMaxSize")) {
      mTerm(1, AT_, "Error: Property minCellMaxSize is deprecated. Please rename to partitionCellOffspringThreshold.");
    } else if(Context::propertyExists("partitionCellMaxNoOffspring")) {
      if(domainId() == 0)
        cerr << "Warning: Property partitionCellMaxNoOffspring is deprecated. Please rename to "
                "partitionCellOffspringThreshold."
             << endl;
      m_partitionCellOffspringThreshold =
          Context::getBasicProperty<MInt>("partitionCellMaxNoOffspring", AT_, &m_partitionCellOffspringThreshold);
    } else {
      m_partitionCellOffspringThreshold =
          Context::getBasicProperty<MInt>("partitionCellOffspringThreshold", AT_, &m_partitionCellOffspringThreshold);
    }

    m_partitionCellWorkloadThreshold = 50000.;
    /*! \page propertyPage1
      \section partitionCellWorkloadThreshold
      <code>MFloat CartesianGrid::m_partitionCellWorkloadThreshold </code>\n
      default = <code>50000</code>\n \n
      controls the filtering of cells in massive parallel file writing\n
      possible values are:
      <ul>
      <li>Non-negative int value of the order less than number of cells on each domains</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, PARALLEL, OUTPUT, FILTER, MAX, MINCELL, SIZE</i>
    */
    m_partitionCellWorkloadThreshold =
        Context::getBasicProperty<MFloat>("partitionCellWorkloadThreshold", AT_, &m_partitionCellWorkloadThreshold);

    m_maxUniformRefinementLevel = -1;
    if(Context::propertyExists("maxGeometricalRfnLvl")) {
      m_maxUniformRefinementLevel = Context::getBasicProperty<MInt>("maxGeometricalRfnLvl", AT_);
      if(domainId() == 0)
        cerr << "Warning: Property maxGeometricalRfnLvl is deprecated and should be provided as "
                "maxUniformRefinementLevel by the grid file."
             << endl;
    }

    m_maxRfnmntLvl = Context::getBasicProperty<MInt>("maxRfnmntLvl", AT_);
    if(m_maxRfnmntLvl > 63) {
      mTerm(1, AT_,
            "Error: m_maxRfnmntLvl > 31 not supported. Long bitshift from FPOW2() only provides 63 correct entries.");
    }

    m_allowInterfaceRefinement = false;
    m_allowInterfaceRefinement =
        Context::getBasicProperty<MBool>("allowInterfaceRefinement", AT_, &m_allowInterfaceRefinement);

    m_cutOff = false;
    m_cutOff = (MBool)Context::getBasicProperty<MInt>("cutOff", AT_, &m_cutOff);

    if(Context::propertyExists("cutOffDirections")) {
      m_cutOff = true;
    }
    /*! \page propertyPage1
      \section periodicCartesianDir
      <code>MInt CartesianGrid::m_newMinLevel </code>\n
      default = <code>0</code>\n \n
      Specify a new minLevel to which the grid will be raised when writing the new restartGrid file!
      <ul>
      <li>Any integer between the old minLevel and the m_maxUniformRefinementLevel!</li>
      </ul>
      Keywords: <i>GRID </i>
    */
    m_newMinLevel = -1;
    m_newMinLevel = Context::getBasicProperty<MInt>("newMinLevel", AT_, &m_newMinLevel);

    m_noPeriodicCartesianDirs = 0;
    m_periodicCartesianDir.fill(0.0);

    /*! \page propertyPage1
      \section periodicCartesianDir
      <code>MFloat CartesianGrid::m_periodicCartesianDir </code>\n
      default = <code>0</code>\n \n
      Space direction in which the grid should be periodic.\n
      possible values are:
      <ul>
      <li>Array with three true/false entries for each Cartesian direction</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, PERIODIC</i>
    */
    if(Context::propertyExists("periodicCartesianDir")) {
      for(MInt dir = 0; dir < nDim; dir++) {
        m_periodicCartesianDir[dir] =
            Context::getBasicProperty<MInt>("periodicCartesianDir", AT_, &m_periodicCartesianDir[dir], dir);
      }
    } else if(Context::propertyExists("periodicDir")) {
      for(MInt dir = 0; dir < nDim; dir++) {
        m_periodicCartesianDir[dir] =
            Context::getBasicProperty<MInt>("periodicDir", AT_, &m_periodicCartesianDir[dir], dir);
      }
    }
    // Read properties requiered for azimuthal periodicy
    m_azimuthalPer = false;
    m_azimuthalPer = Context::getBasicProperty<MBool>("azimuthalPer", AT_, &m_azimuthalPer);
    if(m_azimuthalPer) {
      MInt cnt = 0;
      for(MInt dir = 0; dir < nDim; dir++) {
        m_azimuthalAngle = PI / 180.0 * Context::getBasicProperty<MFloat>("azimuthalAngle", AT_, &m_azimuthalAngle);
        m_azimuthalPerCenter[dir] =
            Context::getBasicProperty<MFloat>("azimuthalPerCenter", AT_, &m_azimuthalPerCenter[dir], dir);
        if(m_periodicCartesianDir[dir] == 0) {
          ASSERT(nDim == 3, "2D grid but periodicCartesianDir is 0!");
          m_azimuthalAxialDir = dir;
        } else {
          ASSERT(cnt < 2, "Azimuthal periodicity is used, but more than 2 peridodic directions!");
          m_azimuthalPeriodicDir[cnt] = dir;
          cnt++;
        }
      }
      ASSERT(cnt == 2, "Azimuthal periodicity is used, but less than 2 peridodic directions!");
    }
    for(MInt dir = 0; dir < nDim; dir++) {
      if(m_periodicCartesianDir[dir]) m_noPeriodicCartesianDirs++;
    }
    m_log << "Number of periodic Cartesian directions is " << m_noPeriodicCartesianDirs << " ("
          << m_periodicCartesianDir[0] << "," << m_periodicCartesianDir[1] << "," << m_periodicCartesianDir[2] << ")"
          << endl;

    // Get output and restart directories
    MString testcaseDir = "./";
    testcaseDir = Context::getBasicProperty<MString>("testcaseDir", AT_, &testcaseDir);

    /*! \page propertyPage1
      \section outputDir
      <code>MInt CartesianGrid::m_outputDir </code>\n
      default = <code></code>\n \n
      Sets the name of the output directory.
      possible values are:
      <ul>
      <li>Strings</li>
      </ul>
      Keywords: <i>GRID, OUTPUT</i>
    */
    m_outputDir = Context::getBasicProperty<MString>("outputDir", AT_);

    m_restart = false;
    m_restart = Context::getBasicProperty<MBool>("restartFile", AT_, &m_restart);

    MBool pp = false;
    pp = Context::getBasicProperty<MBool>("postProcessing", AT_, &pp);

    /*! \page propertyPage1
      \section restartDir
      <code>MInt CartesianGrid::m_outputDir </code>\n
      default = <code>m_outputDir</code>\n \n
      Sets the name of the directory that is used to restart.
      possible values are:
      <ul>
      <li>Strings</li>
      </ul>
      Keywords: <i>GRID, RESTART</i>
    */

    m_restartDir = Context::getBasicProperty<MString>("restartDir", AT_, &m_outputDir);
    if(!m_restart && !pp) {
      m_restartDir = testcaseDir + m_outputDir;
    } else {
      m_restartDir = testcaseDir + m_restartDir;
    }
    m_outputDir = testcaseDir + m_outputDir;

    /*! \page propertyPage1
      \section noHaloLayers
      <code>MInt CartesianGrid::m_noHaloLayers </code>\n
      default = <code>2</code>\n \n
      number of halo layers
      possible values are:
      <ul>
      <li>any non negative integer value</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, FLOW, SOLVER, PARALLEL, HALO, LAYER</i>
    */
    m_noHaloLayers = 2;
    m_noHaloLayers = Context::getBasicProperty<MInt>("noHaloLayers", AT_, &m_noHaloLayers);


    m_loadGridPartition = false;
    /*! \page propertyPage1
      \section loadGridPartition
      <code>MBool CartesianGrid::m_loadGridPartition </code>\n
      default = <code>0</code>\n \n
      enables/disables that a specific grid partitioning is loaded (old concept with two separate
      grids)\n
      possible values are:
      <ul>
      <li>0 (off), 1 (on)</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, FLOW, SOLVER, PARALLEL, PARTITIONING</i>
    */
    m_loadGridPartition = Context::getBasicProperty<MBool>("loadGridPartition", AT_, &m_loadGridPartition);

    /*! \page propertyPage1
      \section loadPartition
      <code>MBool CartesianGrid::m_loadPartition </code>\n
      default = <code>false</code>\n \n
      enables/disables that a specific grid partitioning is loaded
      Keywords: <i>GRID, PARALLEL, PARTITIONING</i>
    */
    m_loadPartition = false;
    m_loadPartition = Context::getBasicProperty<MBool>("loadPartition", AT_, &m_loadPartition);

    /*! \page propertyPage1
      \section partitionParallelSplit
      <code>MBool CartesianGrid::m_partitionParallelSplit</code>\n
      default = <code>false</code>\n \n
      enables/disables the partitioning with a parallel splitting method
      possible values are:
      <ul>
      <li>false (off), true (on)</li>
      </ul>
      Keywords: <i>GRID, PARALLEL, PARTITIONING</i>
    */
    m_partitionParallelSplit = false;
    m_partitionParallelSplit =
        Context::getBasicProperty<MBool>("partitionParallelSplit", AT_, &m_partitionParallelSplit);

    /*! \page propertyPage1
      \section maxNoSurfaces
      <code> MInt CartesianGrid::offset </code>\n
      default = <code>0</code>\n \n
      Set the offset for globalIds to 2,147,483,647 (32Bit). This is necessary to facilitate testcases for global Ids
      exceeding 2,147,483,647 (32Bit). Otherwise testcases would be unreasonably large.
      Possible values are:
      <ul>
        <li>0 - start counting globalIds with 0.</li>
        <li>1 - start counting globalIds with 2,147,483,647</li>
      </ul>
      Keywords: <i>32Bit</i>
    */
    // To test as many routines as possible using the offset, it has to be switched between 0 and
    // std::numeric_limits<MInt>::max(). The place, where the offset has to be set to 0 depends on whether it is a
    // restartFile using a restartGrid, where all globalIds already start with an offset.
    MBool offset = false;
    offset = Context::getBasicProperty<MBool>("bitOffset", AT_, &offset);
    if(offset) m_32BitOffset = std::numeric_limits<MInt>::max();


    // Philipp Brokof, Moritz Waldmann:
    // Initialize status flag for additional LB-grid checks
    /*! \page propertyPage1
      \section lebmGridChecks
      <code>MBool CartesianGrid::m_lbGridChecks </code>\n
      default = <code>false</code>\n \n
      switch for using extra grid checks for LB.\n
      Keywords: <i>GRID, LB, ADAPTATION</i>
    */
    m_lbGridChecks = false;
    if(Context::propertyExists("lbGridChecks")) {
      m_lbGridChecks = Context::getBasicProperty<MBool>("lbGridChecks", AT_);
    }

    /*! \page propertyPage1
      \section coarseRation
      <code>MFloat CartesianGrid::coarseRatio </code>\n
      default = <code>0.2</code>\n \n
      adjusts the factor between refinement and coarsening threshold.\n
      (See BA Philipp Brokof, Implementation of a sensor based adaptive
      mesh refinement ..., variable K_c)\n
      Keywords: <i>GRID, LB, ADAPTATION</i>
    */
    m_coarseRatio = 0.2;
    m_coarseRatio = Context::getBasicProperty<MFloat>("coarseRatio", AT_, &m_coarseRatio);

    /*! \page propertyPage1
      \section coarseRation
      <code>MFloat CartesianGrid::allowCoarsening </code>\n
      default = <code>true</code>\n \n
      Allow coarsening of cells, during adaptation, otherwise cells are only added!)\n
      Keywords: <i>GRID, ADAPTATION</i>
    */
    m_allowCoarsening = true;
    m_allowCoarsening = Context::getBasicProperty<MBool>("allowCoarsening", AT_, &m_allowCoarsening);

    /*! \page propertyPage1
      \section partition cells
      <code>MBool MAIACartesianGrid::updatePartitionCellsOnRestart </code>\n
      default = <code>true</code>\n \n
      Sets whether partition cells shall be updated when writing a restart file.
      Default triggers a new partition when restarting, when setting to false,
      the partition cells are not re-computed when writing the restart file.
      Should be set to FALSE, for LPT applications, as the particle offset
      is based on the current partitioncells and is not recomputed.
      Keywords: <i>GRID, PARTITION</i>
    */

    m_updatePartitionCellsOnRestart = true;
    m_updatePartitionCellsOnRestart =
        Context::getBasicProperty<MBool>("updatePartitionCellsOnRestart", AT_, &m_updatePartitionCellsOnRestart);


    m_lowMemAdaptation = true;
    m_lowMemAdaptation = Context::getBasicProperty<MBool>("lowMemAdaptation", AT_, &m_lowMemAdaptation);
    if(m_lowMemAdaptation) {
      m_log << "Adaptation mode: low memory" << endl;
    } else {
      m_log << "Adaptation mode: default" << endl;
    }

  } else {
    // PARAVIEW PLUGIN
    m_loadGridPartition = 0;
    m_partitionCellOffspringThreshold = 50000;
    m_partitionCellWorkloadThreshold = 50000.;
    m_noPeriodicCartesianDirs = 0;
    m_periodicCartesianDir.fill(0.0);
    m_noHaloLayers = 1;
    m_lbGridChecks = false;
    m_coarseRatio = 0.2;
  }

  createPaths();

  RECORD_TIMER_STOP(timerAllocProp);

  m_log << "Initializing collector | "
        << "#maxCells: " << maxCells << " | #dimensions: " << nDim << endl;

  // Initialize tree
  m_tree.reset(maxCells);
  MInt noSolvers = 1;
  // Set number of solvers
  if(g_multiSolverGrid && !m_paraViewPlugin) {
    noSolvers = Context::getBasicProperty<MInt>("noSolvers", AT_);
    m_log << "Number of solver: " << noSolvers << std::endl;
  } else if(m_paraViewPlugin) {
    // in case of a paraview Plugin we need to reset the number of solvers
    // since no property file is read from the plugin
    m_gridInputFileName = fileName;
    ParallelIo grid(m_gridInputFileName, maia::parallel_io::PIO_READ, mpiComm());
    if(grid.hasAttribute("noSolvers")) {
      grid.getAttribute(&noSolvers, "noSolvers");
    }
  }

  m_maxLevel = 0;

  m_addSolverToGrid = false;
  if(!m_paraViewPlugin) {
    m_addSolverToGrid = Context::getBasicProperty<MBool>("addSolverToGrid", AT_, &m_addSolverToGrid);
    m_referenceSolver = Context::getBasicProperty<MInt>("referenceSolver", AT_, &m_referenceSolver);

    if(domainId() == 0 && m_addSolverToGrid) {
      cerr << "Adding additional solver to the grid file based on solver " << m_referenceSolver << endl;
    }

    if(!g_multiSolverGrid && m_addSolverToGrid) {
      noSolvers++;
    }
  }

  if(domainId() == 0) {
    cerr << "Grid consisting of " << noSolvers << " solvers." << endl;
  }

  m_tree.setNoSolvers(noSolvers);

  // assemble grid input file name
  // note: requires noSolvers information in the grid-tree
  if(!m_paraViewPlugin) {
    m_gridInputFileName = "";
    setGridInputFilename();
  }

  m_log << "Collector initialized" << endl;

  // --- Properties related to solver specific halo generation --- //
  /*! \page propertyPage1
    \section haloMode
    <code>haloMode </code>\n
    default = <code>0</code>\n \n
    Possible values are:
    <ul>
      <li><code>0</code> Just takes all window/halo cells created by raw grid without accounting for
                         varying # of haloLayers (deprecated)</li>
      <li><code>1</code> Creates solver specific #halos</li>
      <li><code>2</code> Same as <code>1</code>, but addiditionally checks connectivity on leaf level
                         and deletes unecessary cells on all lower levels(recommended)</li>
    </ul>
    Keywords: <i> GRID, WINDOW/HALO </i>
  */
  m_haloMode = 0;
  if(!m_paraViewPlugin) {
    m_haloMode = Context::getBasicProperty<MInt>("haloMode", AT_, &m_haloMode);
    TERMM_IF_NOT_COND(m_haloMode >= 0 && m_haloMode <= 2, "Check your input, dude!");
    if(m_haloMode == 0 && domainId() == 0)
      cerr << "\033[0;31m#### WARNING:\033[0m haloMode==0 is deprecated!!!" << endl;
  } else {
    m_haloMode = 0;
  }

  // read Solver specific number of halo layer
  mAlloc(m_noSolverHaloLayers, treeb().noSolvers(), "m_noSolverHaloLayers", AT_);
  for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
    // initialise with default
    m_noSolverHaloLayers[solver] = m_noHaloLayers;
    if(!m_paraViewPlugin && Context::propertyExists("noHaloLayers", solver)) {
      m_noSolverHaloLayers[solver] = Context::getSolverProperty<MInt>("noHaloLayers", solver, AT_);
      ASSERT(m_noSolverHaloLayers[solver] < WINDOWLAYER_MAX, "");
      m_noHaloLayers = std::max(m_noHaloLayers, m_noSolverHaloLayers[solver]);
      ASSERT(m_noSolverHaloLayers[solver] <= m_noHaloLayers, "");
      m_log << " -- solverId " << solver << " -- no Solver Halo-Layers : " << m_noSolverHaloLayers[solver] << endl;
    }
  }

  if(!m_paraViewPlugin && Context::propertyExists("checkRefinementHoles")) {
    mAlloc(m_checkRefinementHoles, treeb().noSolvers(), "m_checkRefinementHoles", AT_);
    /*! \page propertyPage1
     \section checkRefinementHoles
     <code>MBool* CartesianGrid::m_checkRefinementHoles </code>\n
     default = <code>nullptr</code>\n \n
     Set bool for each solver to check for refinement holes and islands
     If activated less level-Jumps are created by the solution-adaptive mesh-refinement!
     Keywords: <i>GRID, ADAPTATION</i>
   */
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      m_checkRefinementHoles[solver] = Context::getSolverProperty<MBool>("checkRefinementHoles", solver, AT_);
    }

    mAlloc(m_diagSmoothing, treeb().noSolvers(), "m_diagSmoothing", AT_);
    /*! \page propertyPage1
     \section checkRefinementHoles
     <code>MBool* CartesianGrid::m_diagSmoothing </code>\n
     default = <code>nullptr</code>\n \n
     Set bool for each solver to check for refinement holes and islands
     If activated less level-Jumps are created by the solution-adaptive mesh-refinement!
     Keywords: <i>GRID, ADAPTATION</i>
   */
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      m_diagSmoothing[solver] = false;
      m_diagSmoothing[solver] =
          Context::getSolverProperty<MBool>("diagSmoothing", solver, AT_, &m_diagSmoothing[solver]);
    }
  }

  if(!m_paraViewPlugin && Context::propertyExists("identicalSolvers")) {
    m_noIdenticalSolvers = Context::propertyLength("identicalSolvers") / 2;
    ASSERT(m_noIdenticalSolvers <= treeb().noSolvers(), "");
    mAlloc(m_identicalSolvers, m_noIdenticalSolvers * 2, "m_identicalSolvers", AT_);
    mAlloc(m_identicalSolverMaxLvl, m_noIdenticalSolvers, "m_identicalSolverMaxLvl", AT_);
    mAlloc(m_identicalSolverLvlJumps, m_noIdenticalSolvers, "m_identicalSolverLvlJumps", AT_);
    for(MInt s = 0; s < m_noIdenticalSolvers * 2; s++) {
      m_identicalSolvers[s] = Context::getBasicProperty<MInt>("identicalSolvers", AT_, s);
    }
    for(MInt pair = 0; pair < m_noIdenticalSolvers; pair++) {
      const MInt slaveId = m_identicalSolvers[pair * 2];
      const MInt masterId = m_identicalSolvers[pair * 2 + 1];
      MInt slaveMaxLvl = Context::getSolverProperty<MInt>("maxRfnmntLvl", slaveId, AT_);
      if(Context::propertyExists("identicalSolverMaxLevel")) {
        slaveMaxLvl = Context::getSolverProperty<MInt>("identicalSolverMaxLevel", slaveId, AT_);
      }
      const MInt masterMaxLvl = Context::getSolverProperty<MInt>("maxRfnmntLvl", masterId, AT_);
      m_identicalSolverMaxLvl[pair] = slaveMaxLvl;
      m_identicalSolverLvlJumps[pair] = masterMaxLvl - slaveMaxLvl;
      cerr0 << "Slave solverId " << slaveId << " has maxRefinementLevel " << slaveMaxLvl << " and "
            << m_identicalSolverLvlJumps[pair] << " refinement-Jumps towards the master!" << endl;
    }
  }

  m_log << "/***** STARTING FLOW SOLVER ******/" << endl;
  loadGridFile(timerLoadGridPar, timerCreateMB);

  if(!m_paraViewPlugin) {
    // After grid was loaded, initialize grid mapping if necessary
    if(Context::propertyExists("donorGridFileName")) {
      initGridMap();
    }

    if(m_addSolverToGrid) {
      if(Context::propertyExists("solverBoundingBox")) {
        if(domainId() == 0) {
          cerr << "Applying solver bounding box to solver " << treeb().noSolvers() - 1 << endl;
        }
        std::array<MFloat, nDim * 2> boxCoordinates;
        for(MInt i = 0; i < nDim * 2; i++) {
          boxCoordinates[i] = Context::getBasicProperty<MFloat>("solverBoundingBox", AT_, i);
        }
        MInt solverMaxLvl = maxRefinementLevel();
        solverMaxLvl = Context::getBasicProperty<MInt>("solverMaxLevel", AT_, &solverMaxLvl);
        // loop top down and remove solver property if:
        //- the cell is outside the box and
        //- none of the cells children belong to the solver!
        const MInt solverId = treeb().noSolvers() - 1;
        for(MInt lvl = maxRefinementLevel(); lvl >= minLevel(); lvl--) {
          for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
            if(a_level(cellId) != lvl) continue;
            if(!treeb().solver(cellId, treeb().noSolvers() - 1)) continue;
            if(a_hasChildren(cellId, treeb().noSolvers() - 1)) continue;
            for(MInt i = 0; i < nDim; i++) {
              if(a_coordinate(cellId, i) < boxCoordinates[i] || a_coordinate(cellId, i) > boxCoordinates[nDim + i]) {
                treeb().solver(cellId, treeb().noSolvers() - 1) = false;
              }
            }
            if(a_level(cellId) > solverMaxLvl) {
              treeb().solver(cellId, solverId) = false;
            }
          }

          if(Context::propertyExists("addedSolverMaxLevel")) {
            MInt addedSolverMaxLevel = Context::getBasicProperty<MInt>("addedSolverMaxLevel", AT_);
            for(MInt cellId = 0; cellId < treeb().size(); cellId++) {
              if(!treeb().solver(cellId, treeb().noSolvers() - 1)) continue;
              if(a_level(cellId) > addedSolverMaxLevel) {
                treeb().solver(cellId, treeb().noSolvers() - 1) = false;
              }
            }
            MInt* recalcIds = nullptr;
            mAlloc(recalcIds, maxNoCells(), "recalcIds", -1, AT_);

            ASSERT(m_tree.noSolvers() > 1, "");
            g_multiSolverGrid = true;

            saveGrid((m_outputDir + m_gridInputFileName).c_str(), recalcIds);
          }
        }
        // exchange solver bitset
        // WH_old
        if(m_haloMode > 0) {
          exchangeSolverBitset(&m_tree.solverBits(0));
        } else {
          maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
                                    m_tree.size());
        }
      }
      MInt* recalcIds = nullptr;
      mAlloc(recalcIds, maxNoCells(), "recalcIds", -1, AT_);
      saveGrid((m_outputDir + m_gridInputFileName + "+addSolver").c_str(), recalcIds);
    }

    // Compute and store local bounding box
    computeLocalBoundingBox(&m_localBoundingBox[0]);
  }

  if(noDomains() > 1) {
    m_log << " -- minLevel for current domain: " << m_minLevel << endl;
    m_log << " -- maxLevel for current domain: " << m_maxLevel << endl;
  } else {
    m_log << " -- minLevel = " << m_minLevel << endl;
    m_log << " -- maxLevel = " << m_maxLevel << endl;
  }

  m_log << " [ Grid file has " << m_tree.size() << " cells ] " << endl;

  RECORD_TIMER_STOP(timertotal);
  DISPLAY_TIMER(timertotal);

  // Store grid cell volume
  // maxCellVolumeLevel is calculated up to level 32 because at some parts of the code
  // the gridCellVolume up to level 32 is needed
  MInt maxCellVolumeLevel = (m_maxLevel < 32) ? 32 : m_maxLevel;
  mAlloc(m_gridCellVolume, maxCellVolumeLevel, "m_gridCellVolume", AT_);
  for(MInt level = 0; level < maxCellVolumeLevel; level++) {
    m_gridCellVolume[level] = F1;
    for(MInt spaceId = 0; spaceId < nDim; spaceId++) {
      m_gridCellVolume[level] *= cellLengthAtLevel(level);
    }
  }

  // check that all minLevels are refined

  if(!m_paraViewPlugin && m_newMinLevel > 0) {
    // check if the grid is already sufficiently refiend!
    MInt noUnrefinedMinCells = 0;
    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_level(cellId) < m_newMinLevel) {
        if(a_noChildren(cellId) < 1) noUnrefinedMinCells++;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &noUnrefinedMinCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "noUnrefinedMinCells");

    if(noUnrefinedMinCells > 0) {
      if(domainId() == 0) {
        cerr << "Mincells not yet sufficiently refined for minLevel increase from " << m_minLevel << " to "
             << m_newMinLevel << "!" << endl;
      }
    } else if(noDomains() == 1) {
      // rename the file-name from which the restart was made
      std::rename((m_restartDir + m_gridInputFileName).c_str(),
                  (m_restartDir + m_gridInputFileName + "backup").c_str());

      m_targetGridMinLevel = m_newMinLevel;
      m_updatedPartitionCells = false;
      MInt* recalcIds = nullptr;
      mAlloc(recalcIds, maxNoCells(), "recalcIds", -1, AT_);
      if(m_maxUniformRefinementLevel < m_newMinLevel) {
        m_maxUniformRefinementLevel = m_newMinLevel;
      }

      saveGrid((m_outputDir + m_gridInputFileName).c_str(), recalcIds);

      cerr0 << "Successfully written new grid file! Exiting MAIA." << endl;
      mTerm(0, AT_, "Finished writing new grid restartFile!");
    }
  }

  if(!m_paraViewPlugin) {
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      MInt inactiveMinCells = 0;
      for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
        if(a_level(cellId) != minLevel()) continue;
        if(a_isHalo(cellId)) continue;
        if(!a_solver(cellId, solver)) {
          inactiveMinCells++;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &inactiveMinCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                    "inactiveMinCells");

      if(inactiveMinCells > 0) {
        cerr0 << "Solver " << solver << " has " << inactiveMinCells << " inactive mincells!" << endl;
      }
    }
  }
  printAllocatedMemory(oldAllocatedBytes, "CartesianGrid", mpiComm());
}


/** Destructor
 *
 */
template <MInt nDim>
CartesianGrid<nDim>::~CartesianGrid() {
  TRACE();

  mDeallocate(m_paths1D);
  mDeallocate(m_paths2D);
  mDeallocate(m_paths3D);
  mDeallocate(m_neighborList);
}


//-----------------------------------------------------------------------------


/** \brief Read grid file name from properties or restart files
    \author Lennart Schneiders
  */
template <MInt nDim>
void CartesianGrid<nDim>::setGridInputFilename(MBool forceRestart) {
  TRACE();

  stringstream s;
  MString name;
  MBool tmpFalse = false;

  const MBool restart = forceRestart || m_restart;

  const MBool useNonSpecifiedRestartFile =
      restart ? Context::getBasicProperty<MBool>("useNonSpecifiedRestartFile", AT_, &tmpFalse) : true;

  const MInt restartTimeStep = useNonSpecifiedRestartFile ? 0 : Context::getBasicProperty<MInt>("restartTimeStep", AT_);
  MBool multilevel = false;
  m_zonal = false;
  multilevel = Context::getBasicProperty<MBool>("multilevel", AT_, &multilevel);
  m_zonal = Context::getBasicProperty<MBool>("zonal", AT_, &m_zonal);

  // NOTE: gridSolverId is the reference solverId used by the grid to read the gridFileName from the
  //      correcponding restartVariables file!
  const MInt gridSolverId = 0;

  MString solverType = Context::getSolverProperty<MString>("solvertype", gridSolverId, AT_);


  if(solverType == "MAIA_LS_SOLVER" || solverType == "MAIA_MULTI_LS") {
    // the ls solver needs special treatment right now since cartesian grid is expected but it can also have its own
    // grid here the cartesian grid file name is required... ls solver right now has to have a cartesian grid file name
    // 'grid.Netcdf' this will be adjusted in the near future
    if(domainId() == 0) {
      cerr << "cartesian grid file name has to be 'grid.Netcdf'" << endl;
    }
    m_gridInputFileName = "grid.Netcdf";
    MBool m_adaptation = Context::getBasicProperty<MBool>("adaptation", AT_, &m_adaptation);
    if(m_adaptation && m_restart) {
      MInt restartTime = Context::getBasicProperty<MInt>("restartTimeStep", AT_);

      stringstream reinitFile;
      reinitFile << "restartGrid_00" << restartTime << ".Netcdf";

      m_gridInputFileName = reinitFile.str();
      cerr << "grid file name changed to : " << m_gridInputFileName << endl;
    }

  } else {
    if(restart) {
      stringstream varFileName;
      if(solverType == "MAIA_LATTICE_BOLTZMANN") {
        // Get restart filename for LB
        varFileName << m_restartDir;
        varFileName << "restart_" << restartTimeStep << ParallelIo::fileExt();
      } else if(solverType == "MAIA_DISCONTINUOUS_GALERKIN") {
        // Get restart filename for DG
        if(useNonSpecifiedRestartFile) {
          varFileName << m_restartDir << "restart" << ParallelIo::fileExt();
        } else if(g_multiSolverGrid) {
          varFileName << m_restartDir << "restart_b" << gridSolverId << "_" << setw(8) << setfill('0')
                      << restartTimeStep << ParallelIo::fileExt();
        } else {
          varFileName << m_restartDir << "restart_" << setw(8) << setfill('0') << restartTimeStep
                      << ParallelIo::fileExt();
        }
      } else if(solverType == "MAIA_DG_MULTISOLVER") {
        // Get restart filename for DG
        if(useNonSpecifiedRestartFile) {
          varFileName << m_restartDir << "restart_b0" << ParallelIo::fileExt();
        } else {
          varFileName << m_restartDir << "restart_b0_" << setw(8) << setfill('0') << restartTimeStep
                      << ParallelIo::fileExt();
        }
      } else if(solverType == "MAIA_LEVELSET_SOLVER") {
        varFileName << m_restartDir << "restartLSCG";
        if(treeb().noSolvers() > 1) varFileName << "_" << gridSolverId;
        if(!useNonSpecifiedRestartFile) varFileName << "_" << restartTimeStep;
        varFileName << ParallelIo::fileExt();
      } else {
        // Get restart filename for FV
        if(!multilevel) {
          if(m_zonal) {
            varFileName << m_restartDir << "restartVariables" << gridSolverId;
          } else {
            varFileName << m_restartDir << "restartVariables";
          }
          if(!useNonSpecifiedRestartFile) varFileName << "_" << restartTimeStep;
          varFileName << ParallelIo::fileExt();
        } else {
          // For multilevel, use solver-specific file name for restart files
          const MInt maxLength = 256;
          array<MChar, maxLength> buffer{};
          snprintf(buffer.data(), maxLength, "restart_b00_t%08d", restartTimeStep);
          varFileName << m_restartDir << MString(buffer.data()) << ParallelIo::fileExt();
        }
      }
      const MString fileName = varFileName.str();
      m_log << "Trying to read grid input file name from restart file " << fileName << endl;
      if(fileExists(fileName.c_str())) {
        ParallelIo parallelIo(fileName.c_str(), maia::parallel_io::PIO_READ, mpiComm());
        parallelIo.getAttribute(&name, "gridFile");
        s << name;
      } else {
        MString tmp = Context::getBasicProperty<MString>("gridInputFileName", AT_);
        s << tmp;
      }
    } else {
      MString tmp = Context::getBasicProperty<MString>("gridInputFileName", AT_);
      s << tmp;
    }
    m_gridInputFileName = s.str();
    struct stat buffer {};
    if(restart && stat((m_restartDir + m_gridInputFileName).c_str(), &buffer) != 0) {
      mTerm(1, AT_, "Grid input file " + m_restartDir + m_gridInputFileName + " does not exist.");
    }
    m_log << "Cartesian Grid input file name is " << m_gridInputFileName
          << " (full path: " << m_restartDir + m_gridInputFileName << ")." << endl;
  }
}
//-----------------------------------------------------------------------------


/** \brief Load grid from disk and setup communication
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
void CartesianGrid<nDim>::loadGridFile(const MInt timerLoadGridPar, const MInt timerCreateMB) {
  TRACE();

  // 1. read grid from disk
  RECORD_TIMER_START(timerLoadGridPar);
  loadGrid(m_restartDir + m_gridInputFileName);
  RECORD_TIMER_STOP(timerLoadGridPar);

  // ---
  // At this point, parentIds, childIds, and nghbrIds are local ids with respect to the collector

  // 2. create window and halo cells
  if(noDomains() > 1 || m_noPeriodicCartesianDirs > 0) {
    RECORD_TIMER_START(timerCreateMB);
    // also calls storeMinLevelCells and computeLeafLevel!
    setupWindowHaloCellConnectivity();
    RECORD_TIMER_STOP(timerCreateMB);
  } else {
    storeMinLevelCells();
    computeLeafLevel();
  }

#ifdef MAIA_GRID_SANITY_CHECKS
  gridSanityChecks();
  checkWindowHaloConsistency(true);
#endif
}


//-----------------------------------------------------------------------------


/** \brief Find domain id containing a given global cell id
 *
 * \author ansgar
 * \date Nov 2021
 *
 * \param[in] globalId global cell id
 * \return domain id containing globalId
 *
 * Note: this function utilized a linear search before, which was quite inefficient for large cases on O(100000) cores
 * and e.g. slowed down the window/halo cell generation significantly
 */
template <MInt nDim>
MInt CartesianGrid<nDim>::findNeighborDomainId(const MLong globalId, const MInt noDomains, const MLong* domainOffsets) {
  if(globalId < domainOffsets[0] || globalId > domainOffsets[noDomains]) {
    TERMM(1, "Invalid global id: " + std::to_string(globalId)
                 + "; domainOffset[0]=" + std::to_string(domainOffsets[noDomains])
                 + "; domainOffset[noDomains]=" + std::to_string(domainOffsets[noDomains]));
    return -1;
  }

  // Search for global cell id in (the sorted) domain offsets array in logarithmic time
  auto lowerBound = std::lower_bound(&domainOffsets[0], &domainOffsets[0] + noDomains, globalId);
  const MInt dist = std::distance(&domainOffsets[0], lowerBound);
  // Check if this cell is a domain offset (i.e. in the offsets list)
  const MBool isDomainOffset = (*lowerBound == globalId);
  // Determine neighbor domain id
  const MInt domain = (isDomainOffset) ? dist : dist - 1;
  return domain;
}


//-----------------------------------------------------------------------------


/** \brief Find index in array [first,last) with matching entry 'val'
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
template <class ITERATOR, typename U>
MInt CartesianGrid<nDim>::findIndex(ITERATOR first, ITERATOR last, const U& val) {
  if(distance(first, last) <= 0) return -1;
  auto idx = (MInt)distance(first, find(first, last, val));
  if(idx == (MInt)distance(first, last)) idx = -1;
  return idx;
}


//-----------------------------------------------------------------------------


/** \brief Find neighbor domain index or create new if not existing
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
template <typename TA, typename TB, typename TC>
MInt CartesianGrid<nDim>::setNeighborDomainIndex(const MInt ndom, vector<TA>& vecA, vector<TB>& vecB,
                                                 vector<TC>& vecC) {
  ASSERT(ndom > -1, "");
  if(m_nghbrDomainIndex[ndom] < 0) {
    vecA.emplace_back();
    vecB.emplace_back();
    vecC.emplace_back();
    m_nghbrDomainIndex[ndom] = (signed)m_nghbrDomains.size();
    m_nghbrDomains.push_back(ndom);
  }
  return m_nghbrDomainIndex[ndom];
}

/** \brief Find neighbor domain index or create new if not existing
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
template <typename TA, typename TB>
MInt CartesianGrid<nDim>::setNeighborDomainIndex(const MInt ndom, vector<TA>& vecA, vector<TB>& vecB) {
  ASSERT(ndom > -1, "");
  if(m_nghbrDomainIndex[ndom] < 0) {
    vecA.emplace_back();
    vecB.emplace_back();
    m_nghbrDomainIndex[ndom] = (signed)m_nghbrDomains.size();
    m_nghbrDomains.push_back(ndom);
  }
  return m_nghbrDomainIndex[ndom];
}

/** \brief Find azimuthal neighbor domain index or create new if not existing
    \author Thomas Hoesgen
    \date March 2021
  */
template <MInt nDim>
MInt CartesianGrid<nDim>::setAzimuthalNeighborDomainIndex(const MInt ndom, vector<vector<MInt>>& vecA,
                                                          vector<vector<MInt>>& vecB) {
  ASSERT(ndom > -1, "");
  if(m_azimuthalNghbrDomainIndex[ndom] < 0) {
    vecA.emplace_back();
    vecB.emplace_back();
    m_azimuthalHigherLevelConnectivity.emplace_back();
    m_azimuthalNghbrDomainIndex[ndom] = (signed)m_azimuthalNghbrDomains.size();
    m_azimuthalNghbrDomains.push_back(ndom);
  }
  return m_azimuthalNghbrDomainIndex[ndom];
}

/** \brief Find azimuthal neighbor domain index or create new if not existing
    \author Thomas Hoesgen
    \date March 2021
  */
template <MInt nDim>
MInt CartesianGrid<nDim>::setAzimuthalNeighborDomainIndex(const MInt ndom, vector<vector<MInt>>& vecA,
                                                          vector<vector<MLong>>& vecB) {
  ASSERT(ndom > -1, "");
  if(m_azimuthalNghbrDomainIndex[ndom] < 0) {
    vecA.emplace_back();
    vecB.emplace_back();
    m_azimuthalHigherLevelConnectivity.emplace_back();
    m_azimuthalNghbrDomainIndex[ndom] = (signed)m_azimuthalNghbrDomains.size();
    m_azimuthalNghbrDomains.push_back(ndom);
  }
  return m_azimuthalNghbrDomainIndex[ndom];
}

//-----------------------------------------------------------------------------


/** \brief Create window and halo cell connectivity between domains
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
void CartesianGrid<nDim>::setupWindowHaloCellConnectivity() {
  TRACE();

  auto logDuration = [this](const MFloat timeStart, const MString comment) {
    logDuration_(timeStart, "WINDOW/HALO", comment, mpiComm(), domainId(), noDomains());
  };
  const MFloat winHaloTimeStart = wallTime();

  std::vector<MInt>().swap(m_nghbrDomains);
  std::vector<MInt>().swap(m_nghbrDomainIndex);
  std::vector<std::vector<MInt>>().swap(m_haloCells);
  std::vector<std::vector<MInt>>().swap(m_windowCells);
  std::vector<std::unordered_map<MInt, M32X4bit<true>>>().swap(m_windowLayer_);

  cerr0 << "Establish window/halo connectivity..." << endl;

  /* if( Context::propertyExists("targetGridFileName") ) { */
  /*   TERMM(1, "deprecated"); */
  /* MString targetGridFileName = Context::getBasicProperty<MString>("targetGridFileName",
   * AT_, &targetGridFileName); */

  std::vector<MInt>().swap(m_nghbrDomains);
  m_nghbrDomainIndex.resize(noDomains());
  std::fill_n(m_nghbrDomainIndex.begin(), noDomains(), -1);

  if(m_azimuthalPer) {
    m_azimuthalNghbrDomains.clear();
    m_azimuthalNghbrDomainIndex.resize(noDomains());
    std::fill_n(m_azimuthalNghbrDomainIndex.begin(), noDomains(), -1);

    m_azimuthalHaloCells.clear();
    m_azimuthalUnmappedHaloCells.clear();
    m_azimuthalWindowCells.clear();
  }

  // 1. create window-halo-connectivity on m_minLevel
  const MFloat minExchangeTimeStart = wallTime();
  createMinLevelExchangeCells();
  logDuration(minExchangeTimeStart, "Create min level exchange cells");

  // 2. iteratively create window-halo-connectivity on m_minLevel+1...m_maxLevel
  const MFloat higherExchangeTimeStart = wallTime();
  if(m_maxLevel > m_minLevel) {
    // WH_old
    if(m_haloMode > 0)
      createHigherLevelExchangeCells();
    else
      createHigherLevelExchangeCells_old();
  }
  logDuration(higherExchangeTimeStart, "Create higher level exchange cells");

#ifdef MAIA_GRID_SANITY_CHECKS
  checkWindowHaloConsistency();
#endif

  // exchange cell properties
  exchangeProperties();

  // create communication data
  updateHaloCellCollectors();

  // WH_old
  // Exchange solver info
  if(m_haloMode > 0)
    exchangeSolverBitset(&m_tree.solverBits(0));
  else
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
                              m_tree.size());

  if(m_azimuthalPer) {
    if(noAzimuthalNeighborDomains() > 0) {
      maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                                &m_tree.solverBits(0), m_tree.size());
    }
    // Set unmapped halos as true for each solver
    // Solver bit is updated in proxy
    for(MInt i = 0; i < noAzimuthalUnmappedHaloCells(); i++) {
      MInt cellId = azimuthalUnmappedHaloCell(i);
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        m_tree.solver(cellId, solver) = true;
      }
    }
    // To ensure that azimuthal halo cells are fully refined (have all children)
    // or are not refined at all, solver Bits need to be checked!
    correctAzimuthalSolverBits();
  }

  storeMinLevelCells();

  computeLeafLevel();

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      a_hasProperty(m_haloCells[i][j], Cell::IsHalo) = true;
      a_hasProperty(m_haloCells[i][j], Cell::IsWindow) = false;
    }
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      a_hasProperty(m_windowCells[i][j], Cell::IsHalo) = false;
      a_hasProperty(m_windowCells[i][j], Cell::IsWindow) = true;
    }
  }
  if(m_azimuthalPer) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsHalo) = true;
        a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsWindow) = false;
      }
      for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
        a_hasProperty(m_azimuthalWindowCells[i][j], Cell::IsHalo) = false;
        a_hasProperty(m_azimuthalWindowCells[i][j], Cell::IsWindow) = true;
      }
    }
    for(MInt j = 0; j < noAzimuthalUnmappedHaloCells(); j++) {
      a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsHalo) = true;
    }
  }

  createGlobalToLocalIdMapping();

  MFloat cellCnt[4] = {F0, F0, F0, F0};
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    cellCnt[0] += (MFloat)m_haloCells[i].size();
    cellCnt[1] += (MFloat)m_windowCells[i].size();
  }
  cellCnt[2] = (MFloat)m_tree.size();
  cellCnt[3] = (MFloat)noNeighborDomains();
  m_log << "Grid local halo/window cell ratio: " << cellCnt[0] << " " << cellCnt[1] << endl;
  m_log << "Grid local halo/window cell ratio: " << cellCnt[0] << "; relative " << cellCnt[0] / cellCnt[2] << "; "
        << cellCnt[1] << "; relative " << cellCnt[1] / cellCnt[2] << endl;
#ifndef NDEBUG
  MFloat maxHaloCellRatio = cellCnt[0] / cellCnt[2];
  MPI_Allreduce(MPI_IN_PLACE, &maxHaloCellRatio, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxHaloCellRatio");
  m_log << "Grid maximum halo cell ratio: " << maxHaloCellRatio << endl;

  MFloat maxWindowCellRatio = cellCnt[1] / cellCnt[2];
  MPI_Allreduce(MPI_IN_PLACE, &maxWindowCellRatio, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxWindowCellRatio");
  m_log << "Grid maximum window cell ratio: " << maxWindowCellRatio << endl;

  MPI_Allreduce(MPI_IN_PLACE, cellCnt, 4, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "cellCnt");
  cellCnt[0] /= cellCnt[2];
  cellCnt[1] /= cellCnt[2];
  cellCnt[3] /= (MFloat)noDomains();
  m_log << "Grid global halo/window cell ratio: " << cellCnt[0] << " " << cellCnt[1] << endl;
  MInt maxNoNeighborDomains = noNeighborDomains();
  MPI_Allreduce(MPI_IN_PLACE, &maxNoNeighborDomains, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxNoNeighborDomains");
  m_log << "Grid maximum number of neighbor domains: " << maxNoNeighborDomains << ", average: " << cellCnt[3] << endl;
#endif

  cerr0 << "done." << endl;
  logDuration(winHaloTimeStart, "Window/halo total");
}


//-------------------------------------------------------------------------------------------


/** \brief Create window-halo-connectivity based on hilbertIndex matching for regular and periodic exchange cells
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
void CartesianGrid<nDim>::createMinLevelExchangeCells() {
  TRACE();

  if(domainId() == 0) cerr << "  * create minLevel exchange cells...";

  const MInt oldNoCells = m_tree.size();
  ScratchSpace<MInt> noHalos(noDomains(), AT_, "noHalos");
  ScratchSpace<MInt> noWindows(noDomains(), AT_, "noWindows");
  ScratchSpace<MLong> hilbertOffsets(noDomains() + 1, AT_, "hilbertOffsets");
  vector<vector<MInt>> halos;
  vector<vector<MLong>> hilbertIds;
  vector<vector<MLong>> recvHilbertIds;
  MInt noMinLevelCells = 0;
  MFloat bbox[2 * nDim];
  for(MInt i = 0; i < nDim; i++) {
    bbox[i] = numeric_limits<MFloat>::max();
    bbox[nDim + i] = numeric_limits<MFloat>::lowest();
  }
  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    for(MInt i = 0; i < nDim; i++) {
      bbox[i] = mMin(bbox[i], a_coordinate(cellId, i) - F1B2 * cellLengthAtLevel(a_level(cellId)));
      bbox[nDim + i] = mMax(bbox[nDim + i], a_coordinate(cellId, i) + F1B2 * cellLengthAtLevel(a_level(cellId)));
    }
  }
  // Note: use epsilon based on lenght0 (was F1B2 * m_lengthLevel0 +- 1e-12 before), required for
  // multisolver grids with multisolver bounding box and shifted center of gravity
  const MFloat halfLength0 = 0.5 * ((F1 + F1 / FPOW2(30)) * m_lengthLevel0);
  for(MInt i = 0; i < nDim; i++) {
    ASSERT(bbox[i] > m_centerOfGravity[i] - halfLength0 && bbox[nDim + i] < m_centerOfGravity[i] + halfLength0,
           "bbox " + std::to_string(bbox[i]) + "," + std::to_string(bbox[nDim + i]) + ", center = "
               + std::to_string(m_centerOfGravity[i]) + ", halfLenght0 = " + std::to_string(halfLength0));
  }
  MPI_Allreduce(MPI_IN_PLACE, &bbox[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[0]");
  MPI_Allreduce(MPI_IN_PLACE, &bbox[nDim], nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[nDim]");

  if(((MLong)pow(2.0, m_targetGridMinLevel * nDim)) > std::numeric_limits<MLong>::max()) {
    mTerm(1, AT_, "Error: minLevel cell size is bigger than std::numeric_limits<MInt>::max()");
  }

  // Count the min-level cells on this domain (including partition level ancestors)
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    a_isToDelete(cellId) = false;
    if(a_level(cellId) == m_minLevel) {
      noMinLevelCells++;
    }
  }
  ASSERT(noMinLevelCells > 0, "Error: no min-level cell (including partition level ancestors) found.");

  unordered_multimap<MLong, MInt> hilbertToLocal;
  hilbertToLocal.reserve(noMinLevelCells + noMinLevelCells / 5);
  noHalos.fill(0);
  noWindows.fill(0);

  MLong minHilbertIndex = std::numeric_limits<MLong>::max();
  MLong prevHilbertIndex = std::numeric_limits<MLong>::min();
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_level(cellId) != m_minLevel) continue;
    const MLong hilbertId = hilbertIndexGeneric(&a_coordinate(cellId, 0));
    hilbertToLocal.insert(make_pair(hilbertId, cellId));
    if(a_hasProperty(cellId, Cell::IsHalo)) {
      MInt ndom = findNeighborDomainId(a_globalId(cellId));
      ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor), "");
      ASSERT(ndom > -1, "");
      MInt idx = setNeighborDomainIndex(ndom, halos, hilbertIds);
      ASSERT(idx > -1, "");
      halos[idx].push_back(cellId);
      hilbertIds[idx].push_back(hilbertId);
      noHalos[ndom]++;
    }
    ASSERT((cellId < m_noInternalCells) != a_hasProperty(cellId, Cell::IsHalo), "");
    ASSERT(a_hasProperty(cellId, Cell::IsHalo)
               == (a_hasProperty(cellId, Cell::IsPartLvlAncestor) && (cellId >= m_noInternalCells)),
           "");
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    minHilbertIndex = mMin(minHilbertIndex, hilbertId);
    ASSERT(hilbertId > prevHilbertIndex,
           "Min level cells not sorted by Hilbert id: " + to_string(hilbertId) + " > " + to_string(prevHilbertIndex));
    prevHilbertIndex = hilbertId;
  }
  MPI_Allgather(&minHilbertIndex, 1, MPI_LONG, &hilbertOffsets[0], 1, MPI_LONG, mpiComm(), AT_, "minHilbertIndex",
                "hilbertOffsets[0]");
  hilbertOffsets[0] = 0;
  hilbertOffsets[noDomains()] = std::numeric_limits<MLong>::max();

  // some domains may not have any min level cells at all (partition level shift)
  for(MInt cpu = noDomains() - 1; cpu >= 0; cpu--) {
    if(hilbertOffsets[cpu] == std::numeric_limits<MLong>::max()) {
      hilbertOffsets[cpu] = hilbertOffsets[cpu + 1];
    }
    ASSERT(hilbertOffsets[cpu] <= hilbertOffsets[cpu + 1], "");
  }

  // const MInt maxNoAdjacentCells = ICUBE[2*m_noHaloLayers+1];
  const MInt maxNoAdjacentCells = ipow(2 * m_noHaloLayers + 1, 3);
  ScratchSpace<MInt> cellList(maxNoAdjacentCells, AT_, "cellList");
  MInt ifcDirs[m_noDirs];

  vector<pair<MInt, MInt>> interfaceCells;
  for(MInt cellId = 0; cellId < oldNoCells; cellId++) {
    if(a_level(cellId) > m_minLevel) continue;
    for(MInt dir = 0; dir < m_noDirs; dir++) {
      if(a_hasNeighbor(cellId, dir) > 0) continue;
      interfaceCells.emplace_back(cellId, dir);
    }
  }

  // find direct and diagonal neighboring cells
  MUint cnt = 0;
  while(cnt < interfaceCells.size()) {
    fill(&ifcDirs[0], &ifcDirs[0] + m_noDirs, 0);
    const MInt cellId = interfaceCells[cnt].first;
    while(cnt < interfaceCells.size() && cellId == interfaceCells[cnt].first) {
      ifcDirs[interfaceCells[cnt].second] = 1;
      cnt++;
    }
    MInt cellCnt = 0;
    cellList[cellCnt++] = cellId;
    for(MInt dim = 0; dim < nDim; dim++) {
      const MInt cellCnt0 = cellCnt;
      for(MInt ori = 0; ori < 2; ori++) {
        const MInt dir = 2 * dim + ori;
        if(!ifcDirs[dir]) continue;
        for(MInt c = 0; c < cellCnt0; c++) {
          MInt nextId = cellList[c];
          for(MInt layer = 0; layer < m_noHaloLayers; layer++) {
            if(a_hasNeighbor(nextId, dir) > 0)
              nextId = a_neighborId(nextId, dir);
            else
              nextId = createAdjacentHaloCell(nextId, dir, &hilbertOffsets[0], hilbertToLocal, bbox, &noHalos[0], halos,
                                              hilbertIds);
            if(nextId < 0) break;
            cellList[cellCnt++] = nextId;
          }
        }
      }
    }
  }

  MPI_Alltoall(&noHalos[0], 1, MPI_INT, &noWindows[0], 1, MPI_INT, mpiComm(), AT_, "noHalos[0]", "noWindows[0]");

  for(MInt i = 0; i < noDomains(); i++) {
    if(noHalos[i] == 0 && noWindows[i] > 0) {
      ASSERT(m_nghbrDomainIndex[i] < 0, "");
      setNeighborDomainIndex(i, halos, hilbertIds);
    }
  }

  ScratchSpace<MPI_Request> sendReq(noNeighborDomains(), AT_, "sendReq");
  sendReq.fill(MPI_REQUEST_NULL);
  recvHilbertIds.resize(noNeighborDomains());
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    ASSERT(m_nghbrDomains[i] > -1 && m_nghbrDomains[i] < noDomains(), "");
    recvHilbertIds[i].resize(noWindows[m_nghbrDomains[i]]);
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt k = 0; k < noHalos[m_nghbrDomains[i]]; k++) {
      ASSERT(hilbertIds[i][k] >= hilbertOffsets[m_nghbrDomains[i]]
                 && hilbertIds[i][k] < hilbertOffsets[m_nghbrDomains[i] + 1],
             to_string(hilbertIds[i][k]) + " " + to_string(hilbertOffsets[m_nghbrDomains[i]]) + " "
                 + to_string(hilbertOffsets[m_nghbrDomains[i] + 1]));
    }
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Issend(&hilbertIds[i][0], noHalos[m_nghbrDomains[i]], MPI_LONG, m_nghbrDomains[i], 34, mpiComm(), &sendReq[i],
               AT_, "hilbertIds[i][0]");
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Recv(&recvHilbertIds[i][0], noWindows[m_nghbrDomains[i]], MPI_LONG, m_nghbrDomains[i], 34, mpiComm(),
             MPI_STATUS_IGNORE, AT_, "recvHilbertIds[i][0]");
  }
  if(noNeighborDomains() > 0) MPI_Waitall(noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt k = 0; k < noWindows[m_nghbrDomains[i]]; k++) {
      ASSERT(recvHilbertIds[i][k] >= hilbertOffsets[domainId()]
                 && recvHilbertIds[i][k] < hilbertOffsets[domainId() + 1],
             "");
    }
  }

  vector<vector<MLong>> sendGlobalIds;
  vector<vector<MLong>> recvGlobalIds;
  sendGlobalIds.resize(noNeighborDomains());
  recvGlobalIds.resize(noNeighborDomains());
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    ASSERT(m_nghbrDomains[i] > -1 && m_nghbrDomains[i] < noDomains(), "");
    sendGlobalIds[i].resize(noWindows[m_nghbrDomains[i]]);
    recvGlobalIds[i].resize(noHalos[m_nghbrDomains[i]]);
  }

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    m_windowCells.push_back(vector<MInt>());
    for(MInt k = 0; k < noWindows[m_nghbrDomains[i]]; k++) {
      MInt windowId = -1;
      auto range = hilbertToLocal.equal_range(recvHilbertIds[i][k]);
      for(auto it = range.first; it != range.second; ++it) {
        if(!a_hasProperty(it->second, Cell::IsPeriodic) && !a_hasProperty(it->second, Cell::IsHalo)) {
          windowId = it->second;
        }
      }
      ASSERT(windowId > -2 && windowId < m_noInternalCells, to_string(windowId) + " " + to_string(m_noInternalCells));
      if(windowId > -1) {
        ASSERT(a_level(windowId) == m_minLevel, "");
        m_windowCells[i].push_back(windowId);
        sendGlobalIds[i][k] = a_globalId(windowId);
      } else {
        sendGlobalIds[i][k] = -1;
      }
    }
  }

  // send back globalId when the queried cell exists, -1 if not existing
  sendReq.fill(MPI_REQUEST_NULL);
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Issend(&sendGlobalIds[i][0], noWindows[m_nghbrDomains[i]], MPI_LONG, m_nghbrDomains[i], 35, mpiComm(),
               &sendReq[i], AT_, "sendGlobalIds[i][0]");
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MPI_Recv(&recvGlobalIds[i][0], noHalos[m_nghbrDomains[i]], MPI_LONG, m_nghbrDomains[i], 35, mpiComm(),
             MPI_STATUS_IGNORE, AT_, "recvGlobalIds[i][0]");
  }
  if(noNeighborDomains() > 0) MPI_Waitall(noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);

  {
    ScratchSpace<MInt> posMapping(m_tree.size(), AT_, "posMapping");
    posMapping.fill(-1);
    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      posMapping[cellId] = cellId;
    }

    ScratchSpace<MInt> delFlag(m_tree.size(), AT_, "delFlag");
    delFlag.fill(0);
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt k = 0; k < noHalos[m_nghbrDomains[i]]; k++) {
        MInt cellId = halos[i][k];
        MLong globalId = recvGlobalIds[i][k];
        a_globalId(cellId) = globalId;
        if(globalId < 0) {
          ASSERT(cellId >= oldNoCells, "Attempting to delete parent gap halo cell.");
          delFlag[cellId] = 1;
          posMapping[cellId] = -1;
        }
      }
    }

    for(MInt cellId = oldNoCells; cellId < m_tree.size(); cellId++) {
      if(delFlag[cellId]) {
        MInt otherCellId = deleteCell(cellId);
        if(otherCellId > -1) {
          ASSERT(delFlag[otherCellId] || posMapping[otherCellId] == otherCellId, "");
          if(!delFlag[otherCellId]) posMapping[otherCellId] = cellId;
          delFlag[cellId] = delFlag[otherCellId];
          delFlag[otherCellId] = 0;
          cellId--;
        }
      }
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      m_haloCells.push_back(vector<MInt>());
      for(MInt k = 0; k < noHalos[m_nghbrDomains[i]]; k++) {
        MLong globalId = recvGlobalIds[i][k];
        if(globalId > -1) {
          ASSERT(posMapping[halos[i][k]] > -1, "");
          MInt cellId = posMapping[halos[i][k]];
          ASSERT(cellId > -1 && cellId < m_tree.size(), "");
          ASSERT(a_globalId(cellId) == globalId, "");
          m_haloCells[i].push_back(cellId);
        }
      }
    }
  }

  for(MInt i = 0; i < nDim; i++) {
    if(m_periodicCartesianDir[i] > 0) {
      m_periodicCartesianLength[i] = bbox[nDim + i] - bbox[i];
    } else {
#ifdef MAIA_PGI_COMPILER
      m_periodicCartesianLength[i] = numeric_limits<MFloat>::quiet_NaN();
#else
      m_periodicCartesianLength[i] = numeric_limits<MFloat>::signaling_NaN();
#endif
    }
  }

  IF_CONSTEXPR(nDim == 3) {
    if(m_azimuthalPer) {
      // First the azimuthal bounding box is determined
      // It is used determine where azimuthal halos cells are required
      for(MInt i = 0; i < nDim; i++) {
        MFloat dirLength = bbox[i + nDim] - bbox[i];
        MInt nBins = (MInt)(dirLength / cellLengthAtLevel(m_minLevel));
        m_azimuthalBbox.minCoord[i] = bbox[i];
        m_azimuthalBbox.boxLength[i] = dirLength;
        m_azimuthalBbox.nBins[i] = nBins;
      }

      m_azimuthalBbox.init(m_azimuthalAngle, m_azimuthalPeriodicDir[0], m_azimuthalPeriodicDir[1], m_noHaloLayers);
      MInt perSize1 = m_azimuthalBbox.nBins[m_azimuthalBbox.perDir1];
      MInt perSize2 = m_azimuthalBbox.nBins[m_azimuthalBbox.perDir2];
      MFloatScratchSpace minPer1(perSize1, AT_, "minPer1");
      minPer1.fill(std::numeric_limits<MFloat>::max());
      MFloatScratchSpace maxPer1(perSize1, AT_, "maxPer1");
      maxPer1.fill(-std::numeric_limits<MFloat>::max());
      MFloatScratchSpace minPer2(perSize2, AT_, "minPer2");
      minPer2.fill(std::numeric_limits<MFloat>::max());
      MFloatScratchSpace maxPer2(perSize2, AT_, "maxPer2");
      maxPer2.fill(-std::numeric_limits<MFloat>::max());

      MFloatScratchSpace minPer3(perSize1 * perSize2, AT_, "minPer3");
      minPer3.fill(std::numeric_limits<MFloat>::max());
      MFloatScratchSpace maxPer3(perSize1 * perSize2, AT_, "maxPer3");
      maxPer3.fill(-std::numeric_limits<MFloat>::max());

      MFloat coordsCyl[3];
      for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
        if(a_level(cellId) != m_minLevel) continue;

        MInt ind1 = m_azimuthalBbox.index(&a_coordinate(cellId, 0), m_azimuthalBbox.perDir1);
        MInt ind2 = m_azimuthalBbox.index(&a_coordinate(cellId, 0), m_azimuthalBbox.perDir2);
        cartesianToCylindric(&a_coordinate(cellId, 0), coordsCyl);
        minPer1[ind1] = mMin(coordsCyl[1], minPer1[ind1]);
        maxPer1[ind1] = mMax(coordsCyl[1], maxPer1[ind1]);
        minPer2[ind2] = mMin(coordsCyl[1], minPer2[ind2]);
        maxPer2[ind2] = mMax(coordsCyl[1], maxPer2[ind2]);

        minPer3[ind1 * perSize2 + ind2] = mMin(coordsCyl[1], minPer3[ind1 * perSize2 + ind2]);
        maxPer3[ind1 * perSize2 + ind2] = mMax(coordsCyl[1], maxPer3[ind1 * perSize2 + ind2]);
      }
      MPI_Allreduce(MPI_IN_PLACE, &minPer1[0], perSize1, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                    "minPer1[0]");
      MPI_Allreduce(MPI_IN_PLACE, &maxPer1[1], perSize1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "maxPer1[0]");
      MPI_Allreduce(MPI_IN_PLACE, &minPer2[0], perSize2, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                    "minPer2[0]");
      MPI_Allreduce(MPI_IN_PLACE, &maxPer2[1], perSize2, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "maxPer2[0]");

      MPI_Allreduce(MPI_IN_PLACE, &minPer3[0], perSize1 * perSize2, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                    "minPer3[0]");
      MPI_Allreduce(MPI_IN_PLACE, &maxPer3[1], perSize1 * perSize2, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                    "maxPer3[0]");

      MFloat minPhi = minPer1[0];
      MFloat maxPhi = maxPer1[0];
      for(MInt p = 0; p < perSize1; p++) {
        minPhi = mMin(minPhi, minPer1[p]);
        maxPhi = mMax(maxPhi, maxPer1[p]);
      }
      for(MInt p = 0; p < perSize2; p++) {
        minPhi = mMin(minPhi, minPer2[p]);
        maxPhi = mMax(maxPhi, maxPer2[p]);
      }
      MFloat center = F1B2 * (minPhi + maxPhi);

      std::copy_n(&minPer3[0], perSize1 * perSize2, &m_azimuthalBbox.minPer3[0]);
      std::copy_n(&maxPer3[0], perSize1 * perSize2, &m_azimuthalBbox.maxPer3[0]);
      m_azimuthalBbox.center = center;

      // Now create the azimuthal periodic halo cells
      noHalos.fill(0);
      noWindows.fill(0);
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        halos[i].clear();
        hilbertIds[i].clear();
      }
      hilbertToLocal.clear();
      for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
        if(a_level(cellId) != m_minLevel) continue;
        const MLong hilbertId = hilbertIndexGeneric(&a_coordinate(cellId, 0));
        hilbertToLocal.insert(make_pair(hilbertId, cellId));
      }

      // find direct and diagonal neighboring cells
      cnt = 0;
      ScratchSpace<MInt> cellLayer(maxNoAdjacentCells, AT_, "cellLayer");
      while(cnt < interfaceCells.size()) {
        fill(&ifcDirs[0], &ifcDirs[0] + m_noDirs, 0);
        const MInt cellId = interfaceCells[cnt].first;
        while(cnt < interfaceCells.size() && cellId == interfaceCells[cnt].first) {
          ifcDirs[interfaceCells[cnt].second] = 1;
          cnt++;
        }
        MInt cellCnt = 0;
        cellList[cellCnt] = cellId;
        cellLayer[cellCnt++] = 0;
        for(MInt dim = 0; dim < nDim; dim++) {
          const MInt cellCnt0 = cellCnt;
          for(MInt ori = 0; ori < 2; ori++) {
            const MInt dir = 2 * dim + ori;
            if(!ifcDirs[dir]) continue;
            for(MInt c = 0; c < cellCnt0; c++) {
              MInt nextId = cellList[c];
              MInt startLayer = cellLayer[c];
              for(MInt layer = startLayer; layer < m_noHaloLayers; layer++) {
                if(a_hasNeighbor(nextId, dir) > 0)
                  nextId = a_neighborId(nextId, dir);
                else
                  nextId = createAzimuthalHaloCell(nextId, dir, &hilbertOffsets[0], hilbertToLocal, bbox, &noHalos[0],
                                                   halos, hilbertIds);
                if(nextId < 0) break;
                cellList[cellCnt] = nextId;
                cellLayer[cellCnt++] = layer;
              }
            }
          }
        }
      }

      MPI_Alltoall(&noHalos[0], 1, MPI_INT, &noWindows[0], 1, MPI_INT, mpiComm(), AT_, "noHalos[0]", "noWindows[0]");

      for(MInt i = 0; i < noDomains(); i++) {
        if(noHalos[i] == 0 && noWindows[i] > 0) {
          setAzimuthalNeighborDomainIndex(i, halos, hilbertIds);
        }
      }

      ScratchSpace<MPI_Request> sendReqAzi(noAzimuthalNeighborDomains(), AT_, "sendReqAzi");
      sendReqAzi.fill(MPI_REQUEST_NULL);
      recvHilbertIds.resize(noAzimuthalNeighborDomains());
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        ASSERT(m_azimuthalNghbrDomains[i] > -1 && m_azimuthalNghbrDomains[i] < noDomains(), "");
        recvHilbertIds[i].resize(noWindows[m_azimuthalNghbrDomains[i]]);
      }
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt k = 0; k < noHalos[m_azimuthalNghbrDomains[i]]; k++) {
          ASSERT(hilbertIds[i][k] >= hilbertOffsets[m_azimuthalNghbrDomains[i]]
                     && hilbertIds[i][k] < hilbertOffsets[m_azimuthalNghbrDomains[i] + 1],
                 to_string(hilbertIds[i][k]) + " " + to_string(hilbertOffsets[m_azimuthalNghbrDomains[i]]) + " "
                     + to_string(hilbertOffsets[m_azimuthalNghbrDomains[i] + 1]));
        }
      }
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        MPI_Issend(&hilbertIds[i][0], noHalos[m_azimuthalNghbrDomains[i]], MPI_LONG, m_azimuthalNghbrDomains[i], 34,
                   mpiComm(), &sendReqAzi[i], AT_, "hilbertIds[i][0]");
      }
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        MPI_Recv(&recvHilbertIds[i][0], noWindows[m_azimuthalNghbrDomains[i]], MPI_LONG, m_azimuthalNghbrDomains[i], 34,
                 mpiComm(), MPI_STATUS_IGNORE, AT_, "recvHilbertIds[i][0]");
      }
      if(noAzimuthalNeighborDomains() > 0)
        MPI_Waitall(noAzimuthalNeighborDomains(), &sendReqAzi[0], MPI_STATUSES_IGNORE, AT_);
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt k = 0; k < noWindows[m_azimuthalNghbrDomains[i]]; k++) {
          ASSERT(recvHilbertIds[i][k] >= hilbertOffsets[domainId()]
                     && recvHilbertIds[i][k] < hilbertOffsets[domainId() + 1],
                 "");
        }
      }

      sendGlobalIds.resize(noAzimuthalNeighborDomains());
      recvGlobalIds.resize(noAzimuthalNeighborDomains());
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        ASSERT(m_azimuthalNghbrDomains[i] > -1 && m_azimuthalNghbrDomains[i] < noDomains(), "");
        sendGlobalIds[i].resize(noWindows[m_azimuthalNghbrDomains[i]]);
        recvGlobalIds[i].resize(noHalos[m_azimuthalNghbrDomains[i]]);
      }

      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        m_azimuthalWindowCells.push_back(vector<MInt>());
        for(MInt k = 0; k < noWindows[m_azimuthalNghbrDomains[i]]; k++) {
          MInt windowId = -1;
          auto range = hilbertToLocal.equal_range(recvHilbertIds[i][k]);
          for(auto it = range.first; it != range.second; ++it) {
            if(!a_hasProperty(it->second, Cell::IsPeriodic) && !a_hasProperty(it->second, Cell::IsHalo)) {
              windowId = it->second;
            }
          }
          ASSERT(windowId > -2 && windowId < m_noInternalCells,
                 to_string(windowId) + " " + to_string(m_noInternalCells));
          if(windowId > -1) {
            ASSERT(a_level(windowId) == m_minLevel, "");
            m_azimuthalWindowCells[i].push_back(windowId);
            sendGlobalIds[i][k] = a_globalId(windowId);
          } else {
            sendGlobalIds[i][k] = -1;
          }
        }
      }

      // send back globalId when the queried cell exists, -1 if not existing
      sendReqAzi.fill(MPI_REQUEST_NULL);
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        MPI_Issend(&sendGlobalIds[i][0], noWindows[m_azimuthalNghbrDomains[i]], MPI_LONG, m_azimuthalNghbrDomains[i],
                   35, mpiComm(), &sendReqAzi[i], AT_, "sendGlobalIds[i][0]");
      }
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        MPI_Recv(&recvGlobalIds[i][0], noHalos[m_azimuthalNghbrDomains[i]], MPI_LONG, m_azimuthalNghbrDomains[i], 35,
                 mpiComm(), MPI_STATUS_IGNORE, AT_, "recvGlobalIds[i][0]");
      }
      if(noAzimuthalNeighborDomains() > 0)
        MPI_Waitall(noAzimuthalNeighborDomains(), &sendReqAzi[0], MPI_STATUSES_IGNORE, AT_);

      {
        for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
          m_azimuthalHaloCells.push_back(vector<MInt>());
          for(MInt k = 0; k < noHalos[m_azimuthalNghbrDomains[i]]; k++) {
            MInt cellId = halos[i][k];
            MLong globalId = recvGlobalIds[i][k];
            if(globalId > -1) {
              ASSERT(cellId > -1 && cellId < m_tree.size(), "");
              a_globalId(cellId) = globalId;
              m_azimuthalHaloCells[i].push_back(cellId);
            } else { // Still keep it! Cell might be necessary. Solver flag will be adjusted in proxy
              ASSERT(cellId > -1 && cellId < m_tree.size(), "");
              a_globalId(cellId) = -1;
              m_azimuthalUnmappedHaloCells.push_back(cellId);
              m_azimuthalUnmappedHaloDomains.push_back(m_azimuthalNghbrDomains[i]);
            }
          }
        }
      }

      // Determine the gridBndryCells. These cells are the outer cells of the grid
      // Tagging them is necessary for the refinement of the azimuthal periodic halo cells
      std::vector<MInt>().swap(m_gridBndryCells);
      MBool adaptation = false;
      adaptation = Context::getBasicProperty<MBool>("adaptation", AT_, &adaptation);
      if(adaptation) {
        if(domainId() == 0)
          cerr << endl << "Azimuthal periodic cells are not refined at grid bndry if adaption is on!" << endl;
      } else {
        if(domainId() == 0) cerr << "Set gridBndryCells on minLevel." << endl;
        MIntScratchSpace gridBndryCells(m_tree.size(), AT_, "gridBndryCells");
        gridBndryCells.fill(-1);
        MBool isBndry = false;
        for(MInt cellId = 0; cellId < oldNoCells; cellId++) {
          if(a_level(cellId) != m_minLevel) continue;
          isBndry = false;
          for(MInt dir = 0; dir < m_noDirs; dir++) {
            if(a_hasNeighbor(cellId, dir) > 0) continue;
            if(!m_periodicCartesianDir[dir / 2]) {
              MFloat coords[nDim];
              for(MInt d = 0; d < nDim; d++) {
                coords[d] = a_coordinate(cellId, d);
              }
              coords[dir / 2] += ((dir % 2) == 0 ? -F1 : F1) * cellLengthAtLevel(m_minLevel);
              if(coords[dir / 2] < bbox[dir / 2] || coords[dir / 2] > bbox[nDim + dir / 2]) continue;
            }
            isBndry = true;
          }
          if(isBndry) {
            if(!a_isHalo(cellId)) m_gridBndryCells.push_back(cellId);
            gridBndryCells[cellId] = 1;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir) > 0) {
                MInt nghbrId = a_neighborId(cellId, dir);
                if(!a_isHalo(nghbrId)) m_gridBndryCells.push_back(nghbrId);
                gridBndryCells[nghbrId] = 1;
              }
            }
          }
        }
        maia::mpi::exchangeData(m_nghbrDomains, m_windowCells, m_haloCells, mpiComm(), gridBndryCells.data());
        for(MInt i = 0; i < noNeighborDomains(); i++) {
          for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
            if(gridBndryCells[m_windowCells[i][j]] == 1) {
              m_gridBndryCells.push_back(m_windowCells[i][j]);
            }
          }
        }
      }
    }
  }

  // It's very polite to clear the data, but we don't need it
  /*  halos.clear();
    hilbertIds.clear();
    recvHilbertIds.clear();
    hilbertToLocal.clear();*/

  // Set window/halo flags
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      ASSERT(!a_hasProperty(m_haloCells[i][j], Cell::IsWindow), "halo cell is marked as window");
      a_hasProperty(m_haloCells[i][j], Cell::IsHalo) = true;
    }
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      ASSERT(!a_isHalo(m_windowCells[i][j]), "window cell is marked as halo");
      a_hasProperty(m_windowCells[i][j], Cell::IsWindow) = true;
    }
  }
  if(m_azimuthalPer) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        ASSERT(!a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsWindow), "azimuthal halo cell is marked as window");
        a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsHalo) = true;
      }
      for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
        ASSERT(!a_isHalo(m_azimuthalWindowCells[i][j]), "azimuthal window cell is marked as halo");
        a_hasProperty(m_azimuthalWindowCells[i][j], Cell::IsWindow) = true;
      }
    }
    for(MInt j = 0; j < noAzimuthalUnmappedHaloCells(); j++) {
      ASSERT(!a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsWindow),
             "azimuthal halo cell is marked as window");
      a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsHalo) = true;
    }
  }

  if(m_haloMode > 0) {
    tagActiveWindowsAtMinLevel(m_minLevel);
#ifndef NDEBUG
    checkWindowLayer("tagActiveWindowsAtMinLevel completed: ");
#endif
    // Exchange solverBits also in case of 1 solver, because solverBits are used to disable halo cells, which
    // are created, but in tagActiveWindowsAtMinLevel seen to be superfluous
    // if(treeb().noSolvers() > 1) {
    exchangeSolverBitset(&m_tree.solverBits(0));
    //}

    if(m_haloMode == 2) {
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        std::vector<MInt> haloCells;
        haloCells.reserve(m_haloCells[i].size());
        for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
          const MInt haloCellId = m_haloCells[i][j];
          if(m_tree.solverBits(haloCellId).any())
            haloCells.push_back(haloCellId);
          else {
            removeCell<false>(haloCellId);
          }
        }
        m_haloCells[i] = haloCells;

        //
        std::vector<MInt> windowCells;
        windowCells.reserve(m_windowCells[i].size());
        std::set<MInt> erased;
        for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
          const MInt windowCellId = m_windowCells[i][j];
          ASSERT(m_windowLayer_[i].find(windowCellId) != m_windowLayer_[i].end()
                     || erased.find(windowCellId) != erased.end(),
                 "");
          MBool validWindow = false;
          for(MInt solverId = 0; solverId < treeb().noSolvers(); ++solverId) {
            if(isSolverWindowCell(i, windowCellId, solverId)) {
              windowCells.push_back(windowCellId);
              validWindow = true;
              break;
            }
          }
          if(!validWindow) {
            const MBool ok = m_windowLayer_[i].erase(windowCellId);
            ASSERT(ok || erased.find(windowCellId) != erased.end(), "");
            erased.insert(windowCellId);
            MBool isWindow = false;
            for(const auto& w : m_windowLayer_) {
              if(w.find(windowCellId) != w.end()) {
                isWindow = true;
                break;
              }
            }
            if(!isWindow) a_hasProperty(windowCellId, Cell::IsWindow) = false;
          }
        }
        m_windowCells[i] = windowCells;
      }

      compactCells();
    }

#ifndef NDEBUG
    // for debugging, check consistency before new level
    checkWindowHaloConsistency(false, "createMinLevelExchangeCells completed: ");
#endif
  }


  cerr0 << " done." << endl;
}


//-------------------------------------------------------------------------------------------

/** \brief Iteratively create window and halo cells on levels m_minLevel+1...m_maxLevel,
 *         starting from m_minLevel exchange cells
    \author Lennart Schneiders
    \date October 2017
  */
// clang-format off
template <MInt nDim>
void CartesianGrid<nDim>::createHigherLevelExchangeCells(
    const MInt onlyLevel, const MBool duringMeshAdaptation,
    const std::vector<std::function<void(const MInt)>>& refineCellSolver, const std::vector<std::function<void(const MInt)>>& removeCellSolver, const MBool forceLeafLvlCorrection) {
  TRACE();

  auto logDuration = [this] (const MFloat timeStart, const MString comment) {
    logDuration_(timeStart, "WINDOW/HALO HIGHER LEVEL", comment, mpiComm(), domainId(), noDomains());
  };

  ScratchSpace<MInt> plaMap(m_maxNoCells, AT_, "plaMap");
  vector<MLong> refineChildIds;
  vector<MLong> recvChildIds;
  vector<vector<MInt>> haloMap;
  vector<vector<MInt>> windMap;
  vector<MInt> refineIds;

  if(domainId() == 0 && onlyLevel < 0) cerr << "  * create higher level exchange cells...";

  plaMap.fill(-1);
  if(m_maxPartitionLevelShift > 0) {
    for(MUint i = 0; i < m_partitionLevelAncestorIds.size(); i++) {
      plaMap[m_partitionLevelAncestorIds[i]] = i;
    }
  }

  if(m_azimuthalPer && onlyLevel <= m_minLevel) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      m_azimuthalHigherLevelConnectivity[i].clear();
    }
  }

  // given window-halo cell connectivity on minLevel established in setupWindowHaloCellConnectivity(),
  // iteratively create connectivity on ascending levels
  for(MInt level = m_minLevel; level < m_maxLevel; level++) {
    if(onlyLevel > -1 && onlyLevel != level) continue;

    const MFloat levelTimeStart = wallTime();
#ifndef NDEBUG
    // for debugging, check consistency before new level
    checkWindowHaloConsistency(false, "createHigherLevelExchangeCells level=" + to_string(level) + ": ");
#endif

    refineIds.clear();
    MInt haloCnt = 0;
    MInt windowCnt = 0;
    haloMap.resize(noNeighborDomains());
    windMap.resize(noNeighborDomains());

    // Sets to hold halo/window cells for each neighbor domain to allow fast searching if cells are
    // already present
    vector<std::set<MInt>> haloCellSet(noNeighborDomains());
    vector<std::set<MInt>> windCellSet(noNeighborDomains());

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      haloMap[i].resize(m_haloCells[i].size());
      windMap[i].resize(m_windowCells[i].size());
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        haloMap[i][j] = haloCnt++;
      }
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        windMap[i][j] = windowCnt++;
      }

      haloCellSet[i].insert(m_haloCells[i].begin(), m_haloCells[i].end());
      windCellSet[i].insert(m_windowCells[i].begin(), m_windowCells[i].end());
    }

    recvChildIds.resize(haloCnt * m_maxNoChilds);
    fill(recvChildIds.begin(), recvChildIds.end(), -1);
    refineChildIds.clear();
    refineChildIds.resize(windowCnt * m_maxNoChilds);
    fill(refineChildIds.begin(), refineChildIds.end(), -1);

#ifndef NDEBUG
    for (MInt cellId = 0; cellId < m_tree.size(); ++cellId) {
      if (!a_isToDelete(cellId)) {
        if (a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPeriodic)) {
          const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          const MInt maxNghbrDomId =
            findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
          if (domainId()>=minNghbrDomId && domainId()<=maxNghbrDomId) {
            TERMM_IF_NOT_COND(a_hasChildren(cellId), "d=" + to_string(domainId()) + " nghbrs="
                + to_string(minNghbrDomId) + "|" + to_string(maxNghbrDomId) + " isHalo=" + to_string(a_isHalo(cellId)));
          }
        }
      }
    }

    TERMM_IF_COND(bitOffset() && m_maxPartitionLevelShift>0, "Not a good idea to use bitOffset==true & partLvlShift!");
    if (!bitOffset()) {
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        const MInt cellId = m_windowCells[i][j];

        if (!a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
          const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          // I don't know why but some cells in testcase LS/3D_minLevelIncrease have a_noOffsprings(cellId)==0
          const MInt maxNghbrDomId =
            findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)std::max(1,a_noOffsprings(cellId))-1));
          TERMM_IF_NOT_COND(minNghbrDomId==maxNghbrDomId && minNghbrDomId==domainId(),
            "Window cell is not partLvlAncestor, but seems to have children on different domains!");
        }
      }
    }
    }
#endif

    /// --- PartLvlShift (*) --- ///
    // In the following we provide all window cells, whose parent cell is partLvlAncestor and has children on more
    // then one domain, with the information they need and store it in the variable 'windowTemp'

    // In periodic case we might need to add same window cell twice, therefore multimap. Save for all
    // window children of partLvlAncestor parents with children on more then one domain the connectivity to all halos.
    std::multimap<std::pair<MInt/*nghbrDomainIdx*/,MInt/*cellId*/>, M32X4bit<true>> windowsTemp;
    if (m_maxPartitionLevelShift>0) {
      // 1) Refine all partLvlAncestor cells with level+1 children on more than one domain
      std::set<MInt> partLvlAncSet;
      for(MUint i = 0; i < m_partitionLevelAncestorIds.size(); i++) {
        const MInt cellId = m_partitionLevelAncestorIds[i];
        if (a_level(cellId)!=level) continue;

        ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor),
               "not a partition level ancestor: cellId" + std::to_string(cellId) + " halo"
                   + std::to_string(a_hasProperty(cellId, Cell::IsHalo)) + " domain" + std::to_string(domainId()));
        ASSERT(a_hasProperty(cellId, Cell::IsHalo) || findNeighborDomainId(a_globalId(cellId)) == domainId(), "");
        ASSERT(a_noChildren(cellId)>0, "PartLvlAncestor cell must have children!");
        ASSERT((signed)i==plaMap[cellId], "");

        // Guess: during meshAdaptation partLvlAncestor cells may already have halo children
        // Note: if cell is halo cell, we may have the situation that neither of the children on level+1 is on current domain
        std::set<MInt> nghbrDomIds;
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          if(m_partitionLevelAncestorChildIds[i * m_maxNoChilds + child] > -1) {
            nghbrDomIds.insert(findNeighborDomainId(m_partitionLevelAncestorChildIds[i * m_maxNoChilds + child]));
          }
        }
        const MBool multikultiCell = nghbrDomIds.size()>1;

        const MInt noChildren = a_noChildren(cellId);
        if(multikultiCell) {
          partLvlAncSet.insert(cellId);

          MLong* childIds = &m_partitionLevelAncestorChildIds[i * m_maxNoChilds];
          if(accumulate(childIds, childIds + m_maxNoChilds, static_cast<MLong>(-1),
                        [](const MLong& a, const MLong& b) { return std::max(a, b); })
             > -1) {
            refineCell(cellId, childIds, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
            // Note: If function is called during meshAdaptation all children should already exist,
            //       i.e.a_noChildren(cellId)==noChildren
            if (a_noChildren(cellId)>noChildren)
              refineIds.push_back(cellId); //TODO_SS labels:GRID,totest is it correct? is it necessary?

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;

              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              ASSERT(a_globalId(childId) > -1, "");
              const MInt ndom = findNeighborDomainId(a_globalId(childId));
              ASSERT(ndom > -1, "");

              if (a_isHalo(childId)) {
                // There can be partLvlHalo childs and non-partLvlHalo childs (the latter are newly craeted ones)
                ASSERT(ndom!=domainId(), "");
              } else {
                ASSERT(ndom==domainId(), "");
              }
              //TODO_SS labels:GRID,totest check the following
              a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
            }
            // The following check fails during meshAdaptation
            //TERMM_IF_NOT_COND(a_noChildren(cellId)-noChildren==noHaloChilds, "Expected non-partLvlAncestor halo child!");
            // Debug output
/*            stringstream ss;
            ss << " PartLvlAncestor Refinement: d=" << domainId() << " cellId=" << cellId << " (halo=" << a_isHalo(cellId)
              << "): Created " << noHaloChilds << " halos on level=" << level+1;
            for (MInt child = 0; child < m_maxNoChilds; child++) {
              ss << "\n\t1) childId=" << childInfos[3*child];
              if (childInfos[3*child+1]==domainId())
                ss << " isWindow";
              else
                ss << " isHalo nghbrDom=" << childInfos[3*child+1] << " isPartLvlHalo=" << childInfos[3*child+2];
              ///cout << ss.str() << endl;
            }*/

          }
        }
      }

      // 2)
      // Determine all pairs of neighbor domains to connect:
      // For simplicity, if a window cell has children on level+1, which reside on different domains, then the
      // window/halo connectivity of the window cell is propagated to all the children, meaning that if domains
      // 2,5,11,23 have current cell as halo cell, they also have all children as halo cells.
      // Determine window cells on current domain with children on level+1, which reside on different domain (*).
      // Determine all neighbors, which have current window cell as halo cell and inform domain (*) about these neighbors.
      vector<MInt> sendDomIdx;
      vector<M32X4bit<>::type> sendData;
      MBoolScratchSpace cellFlag(m_tree.size(), AT_, "cellFlag");
      fill(cellFlag.begin(), cellFlag.end(), false);
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
          const MInt cellId = m_windowCells[i][j];
          ASSERT(m_windowLayer_[i].find(cellId)!=m_windowLayer_[i].end(), "");

          if (partLvlAncSet.find(cellId)!=partLvlAncSet.end()) {
            ASSERT(a_level(cellId) == level && a_noChildren(cellId) > 0, "");
            // Save all domains which own at least one child of current cell
            std::set<MInt> nghbrs;
            // Boolean to verify that current window cell has at least one window child
            MBool windowChild = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt cnt = windMap[i][j];

              // The purpose of the following is to tag childs which don't reside on current domain because of pls
              const MInt pId = plaMap[cellId];
              ASSERT(pId>-1, "");
              if(m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child] > -1) {
                refineChildIds[m_maxNoChilds * cnt + child] = m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child];
                const MInt domainChild = findNeighborDomainId(m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child]);
                if (domainChild != domainId()) {
                  nghbrs.insert(domainChild);
                } else {
                  windowChild = true;
                  const MInt childId = a_childId(cellId, child);
                  ASSERT(childId>-1, " ");
                  // In periodic case duplicate window cells possible
                  TERMM_IF_COND(windowsTemp.find(std::make_pair(m_nghbrDomains[i],childId))!=windowsTemp.end()
                             && m_noPeriodicCartesianDirs==0, " ");
                  M32X4bit<true> windowLayer = m_windowLayer_[i].find(childId)!=m_windowLayer_[i].end()
                    ? m_windowLayer_[i].at(childId) :  m_windowLayer_[i].at(cellId);//M32X4bit<true>{};
                  windowsTemp.insert({std::make_pair(m_nghbrDomains[i],childId), windowLayer});
                }
              }
            }
            ASSERT(windowChild, "At least one child must reside on current domain!");
            ASSERT(nghbrs.size()>0, "At least one child must reside on a different domain!");

            for (const auto ndom : nghbrs) {
              const MInt idx = findIndex(m_nghbrDomains.begin(), m_nghbrDomains.end(), ndom);
              ASSERT(idx > -1, "");
              ASSERT(m_nghbrDomains[idx] == ndom, "");
              // In periodic case a neighbor domain can have the cell twice, meaning that current window cell is sending
              // the data to two halo cells of the same domain;
              // In periodic case we can also have the situation that a child belongs to a neighbor domain,
              // and that neighbor domain also has that cell as a periodic halo cell.
              if (ndom != m_nghbrDomains[i] || m_noPeriodicCartesianDirs>0) {
                sendDomIdx.push_back(idx);
                sendData.push_back(m_nghbrDomains[i]);
                sendData.push_back(a_globalId(cellId));
                M32X4bit<true> windowLayer = m_windowLayer_[i].at(cellId);
                sendData.push_back(windowLayer.data());
              }
              if (!cellFlag[cellId]) {
                sendDomIdx.push_back(idx);
                sendData.push_back(domainId());
                sendData.push_back(a_globalId(cellId));
                M32X4bit<true> windowLayer = m_windowLayer_[i].at(cellId);
                sendData.push_back(windowLayer.data());
              }
            }
            cellFlag[cellId] = true;
          }
        }
      }

      // exchange redirected communication info
      vector<MInt> recvOffs;
      vector<M32X4bit<>::type> recvData;
      maia::mpi::exchangeScattered(m_nghbrDomains, sendDomIdx, sendData, mpiComm(), recvOffs, recvData, 3);

      // create window cells if this domain is affected by redirected communication links
      const MInt noNgbrDomainsBak = noNeighborDomains();
      for(MInt i = 0; i < noNgbrDomainsBak; i++) {
        set<MInt> check4PeriodicHalos;
        for(MInt j = recvOffs[i]; j < recvOffs[i + 1]; j++) {
          auto ndom = static_cast<MInt>(recvData[3 * j]);
          const MLong globalCellId = recvData[3 * j + 1];
          if (ndom==domainId()) { //might occur in periodic case
            // If in periodic case a window cell appears multiple times, skip it one time
            ASSERT(m_noPeriodicCartesianDirs>0, "");
            if (check4PeriodicHalos.find(globalCellId)==check4PeriodicHalos.end()) {
              check4PeriodicHalos.insert(globalCellId);
              continue;
            }
          }
          auto windowLayer = recvData[3 * j + 2];
          ASSERT(ndom>-1, "");
          ASSERT(m_globalToLocalId.find(globalCellId)!=m_globalToLocalId.end(), "");
          const MInt cellId = m_globalToLocalId[globalCellId];
          ASSERT(cellId > -1, "");
          ASSERT(a_globalId(cellId) == globalCellId, "");
          ASSERT(findNeighborDomainId(globalCellId) != domainId(), "");
//          const MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells, m_windowLayer_);
//          TERMM_IF_NOT_COND(idx > -1, "");

          // New neighbor domain, create new cell sets
//          if(idx == static_cast<MInt>(windCellSet.size())) {
//            windCellSet.emplace_back();
//            haloCellSet.emplace_back();
//          }

          //const MInt pId = plaMap[cellId];
          //TERMM_IF_NOT_COND(pId>-1, "");

          MBool foundWindow = false;
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            const MInt childId = a_childId(cellId, child);
            if(childId < 0 || a_isHalo(childId)) continue;
            foundWindow = true;

            // In periodic case duplicate window cells possible
            ASSERT(windowsTemp.find(std::make_pair(ndom,childId))==windowsTemp.end()
                       || m_noPeriodicCartesianDirs>0, " ");

//            m_windowCells[idx].push_back(childId);
//            windCellSet[idx].insert(childId);
            windowsTemp.insert({std::make_pair(ndom,childId), M32X4bit<true>(windowLayer)});
          }
          ASSERT(foundWindow, "");
        }
      }
    }
    /// --- PartLvlShift Ends --- ///


    // since setNeighborDomainIndex is called in tagActiveWindows2_ we might get new neighbors, so save old state
    std::vector<MInt> nghbrDomains(m_nghbrDomains);

    tagActiveWindows2_(refineChildIds, level);

    // We might have added a neighbor domain
    haloCellSet.resize(noNeighborDomains());
    windCellSet.resize(noNeighborDomains());

// set globalIds of each windowCells's child
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        const MInt cellId = m_windowCells[i][j];

        if(a_level(cellId) == level && a_noChildren(cellId) > 0) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            const MInt cnt = windMap[i][j];

            if(refineChildIds[m_maxNoChilds * cnt + child] > 0) {
              const MInt childId = a_childId(cellId, child);
              ASSERT(childId > -1, "");
              ASSERT((a_globalId(static_cast<MInt>(childId)) >= m_domainOffsets[domainId()]
                      && a_globalId(static_cast<MInt>(childId)) < m_domainOffsets[domainId() + 1])
                         || a_hasProperty(cellId, Cell::IsPartLvlAncestor),
                     to_string(a_globalId(static_cast<MInt>(childId))) + " "
                         + to_string(m_domainOffsets[domainId()]) + " " + to_string(m_domainOffsets[domainId() + 1]));
              refineChildIds[m_maxNoChilds * cnt + child] = a_globalId(static_cast<MInt>(childId));
            }
          }
        }
      }
    }

    // exchange childIds
    maia::mpi::exchangeBuffer(nghbrDomains,//m_nghbrDomains,
                              m_haloCells,
                              m_windowCells,
                              mpiComm(),
                              refineChildIds.data(),
                              recvChildIds.data(),
                              m_maxNoChilds);

#ifndef NDEBUG
    if(m_maxPartitionLevelShift > 0) {
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
          MInt cellId = m_haloCells[i][j];
          if(a_level(cellId) == level) {
            MInt pId = plaMap[cellId];
            if(m_maxPartitionLevelShift > 0 && pId > -1) {
              MInt cnt = haloMap[i][j];
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                const MLong childId = m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child];
                ASSERT(childId == recvChildIds[m_maxNoChilds * cnt + child],
                       std::to_string(childId) + " != " + std::to_string(recvChildIds[m_maxNoChilds * cnt + child])
                           + "; c" + std::to_string(cellId) + "; g" + std::to_string(a_globalId(cellId)) + "; p"
                           + std::to_string(pId));
              }
            }
          }
        }
      }
    }
#endif


    // This is the continuation of 'PartLvlShift (*)'
    if (m_maxPartitionLevelShift>0) {
      for (const auto& item : windowsTemp) {
        // ndom -> domain to which current window cell needs to send the data to
        const MInt ndom = item.first.first;
        // childId -> cellId of the current window cell
        const MInt childId = item.first.second;
        ASSERT(a_hasProperty(a_parentId(childId), Cell::IsPartLvlAncestor), "");
        const M32X4bit<true> windowLayer = item.second;
        const MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells, m_windowLayer_);
        ASSERT(idx > -1, "");
//        TERMM_IF_NOT_COND(a_isHalo(a_parentId(childId)) && a_hasProperty(a_parentId(childId), Cell::IsPartLvlAncestor), "");
        ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell), "");

        // New neighbor domain, create new cell sets
        if(idx == static_cast<MInt>(windCellSet.size())) {
          windCellSet.emplace_back();
          haloCellSet.emplace_back();
        }

        if (m_noPeriodicCartesianDirs==0) {
          if(windCellSet[idx].find(childId) != windCellSet[idx].end())
            continue;
        } else {
          // Note: we can also have 4 halo cells mapped to same window cell (see DG/3D_cube_linearscalaradv_linear_periodic)
          //TERMM_IF_COND(windowsTemp.count(std::make_pair(ndom,childId))>2 || std::count(m_windowCells[idx].begin(), m_windowCells[idx].end(), childId)>2, "");
          if ((signed)windowsTemp.count(std::make_pair(ndom,childId))<=std::count(m_windowCells[idx].begin(), m_windowCells[idx].end(), childId)) //performance wise this is a shit
            continue;
        }

        m_windowCells[idx].push_back(childId);
        windCellSet[idx].insert(childId);

        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if (m_tree.solver(childId, solver) && !isSolverWindowCell(idx,childId,solver) && windowLayer.get(solver)<=m_noSolverHaloLayers[solver]) {
            const MInt w = std::max(windowLayer.get(solver), m_noSolverHaloLayers[solver]);
            m_windowLayer_[idx][childId].set(solver, w/*windowLayer.get(solver)*//*m_noSolverHaloLayers[solver]*/);
          }
        }

        //TODO_SS labels:GRID,toenhance find a better way
        if (m_windowLayer_[idx].find(childId)==m_windowLayer_[idx].end()) {
          for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
//            if (m_tree.solver(childId, solver)) {
              //TERMM_IF_NOT_COND(m_tree.solver(a_parentId(childId), solver), "");
              // We can have cases that a partLvlAnc cell belongs to a solver, but one of its children
              // does not. So we can not just propagate the information from parent cell to children cell.
              // So we invalidate those window cells by setting layer to m_noSolverHaloLayers[solver]+1.
              const MInt w = m_tree.solver(childId, solver) ?
                std::max(windowLayer.get(solver), m_noSolverHaloLayers[solver]):
                m_noSolverHaloLayers[solver]+1;
              m_windowLayer_[idx][childId].set(solver, w/*windowLayer.get(solver)*//*m_noSolverHaloLayers[solver]*//*m_noHaloLayers+1*/);
//            }
          }
        }
      }
    }


    // update regular window cells between this domain and m_nghbrDomains[i] to include the level+1 cells
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      if (m_windowCells[i].size()==0) continue;
      vector<MInt> oldWindowVec(m_windowCells[i]);
      std::set<MInt> oldWindCellSet(windCellSet[i]);

      std::vector<MInt>().swap(m_windowCells[i]);
      windCellSet[i].clear();

      for(MInt j = 0; j < (signed)oldWindowVec.size(); j++) {
        MInt cellId = oldWindowVec[j];

        ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet[i].find(cellId) == windCellSet[i].end(),
               "duplicate window cell: " + std::to_string(cellId) + " " + std::to_string(a_globalId(cellId)) + " "
                   + std::to_string(j));

        m_windowCells[i].push_back(cellId);
        windCellSet[i].insert(cellId);

        ASSERT(cellId < m_tree.size(), to_string(level));

        if(a_level(cellId) == level && a_noChildren(cellId) > 0) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt cnt = windMap[i][j];

            const MInt childId = a_childId(cellId, child);
            if(childId < 0) continue;

            // Skip child if it is already part of the old window cell vector and is handled by the
            // outer loop
            if(oldWindCellSet.find(childId) != oldWindCellSet.end()) {
              ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                     "not a partition level ancestor or partition cell");
              continue;
            }

            const MLong globalChildId = refineChildIds[m_maxNoChilds * cnt + child];
            if(globalChildId > -1) {
              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              // ASSERT( childId > -1 && childId < m_noInternalCells, to_string(childId) + " " +
              // to_string(m_noInternalCells) );
              MInt ndom = findNeighborDomainId(globalChildId);
              if(ndom == domainId()) {
                ASSERT(globalChildId == a_globalId(childId), "");

                ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet[i].find(childId) == windCellSet[i].end(),
                       "duplicate window cell: " + std::to_string(childId));

                m_windowCells[i].push_back(childId);
//                if (m_windowLayer_[i].find(childId)==m_windowLayer_[i].end()) {
                  // I think the following fails because of testingFix_partitionLevelShift
                  //TERMM_IF_NOT_COND(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell), "");
                  // Note: it is dangerous to set windowLayer=m_noHaloLayers+1, since in case
                  //       !isSolverWindowCell(i,childId,solver) for all solvers, we might get
                  //       a halo cell which is completly disabled
//                  for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
//                    if (m_tree.solver(childId, solver) && !isSolverWindowCell(i,childId,solver)) {
//                      m_windowLayer_[i][childId].set(solver, m_noHaloLayers+1/*1*/);
//                    }
//                  }
#ifndef NDEBUG
                  MBool validWindow = false;
                  for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
                    if (isSolverWindowCell(i,childId,solver)) {
                      validWindow = true;
                      break;
                    }
                  }
                  TERMM_IF_NOT_COND(validWindow, "");
#endif
//                }
                windCellSet[i].insert(childId);
              } 
            }
          }
        }
      }
    }


    // create level+1 halo cells existing on other domains
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      if (m_haloCells[i].size()==0) continue;
      vector<MInt> oldHaloVec(m_haloCells[i]);
      std::set<MInt> oldHaloCellSet(haloCellSet[i]);

      std::vector<MInt>().swap(m_haloCells[i]);
      haloCellSet[i].clear();

      for(MInt j = 0; j < (signed)oldHaloVec.size(); j++) {
        MInt cellId = oldHaloVec[j];

        ASSERT(haloCellSet[i].find(cellId) == haloCellSet[i].end(),
               "duplicate halo cell: " + std::to_string(cellId) + " " + std::to_string(a_globalId(cellId)) + " "
                   + std::to_string(j));

        m_haloCells[i].push_back(cellId);
        haloCellSet[i].insert(cellId);

        if(a_level(cellId) == level) {
          MInt cnt = haloMap[i][j];

          if(accumulate(&recvChildIds[m_maxNoChilds * cnt], &recvChildIds[m_maxNoChilds * cnt] + m_maxNoChilds,
                        static_cast<MLong>(-1), [](const MLong& a, const MLong& b) { return std::max(a, b); })
             > -1) {
            refineCell(cellId, &recvChildIds[m_maxNoChilds * cnt], a_hasProperty(cellId, Cell::IsPartLvlAncestor));
            refineIds.push_back(cellId);

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;

              // Skip halo child already present in case of partition level shifts
              if(oldHaloCellSet.find(childId) != oldHaloCellSet.end()) {
                ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                       "not a partition level ancestor or partition cell");
                continue;
              }

              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              ASSERT(a_globalId(childId) > -1, "");
              MInt ndom = findNeighborDomainId(a_globalId(childId));
              ASSERT(ndom > -1, "");
              if(recvChildIds[m_maxNoChilds * cnt + child] > -1)
                ASSERT(a_globalId(childId) == recvChildIds[m_maxNoChilds * cnt + child], "");
              if(recvChildIds[m_maxNoChilds * cnt + child] < 0) ASSERT(ndom == domainId(), "");
              if(ndom == m_nghbrDomains[i]) {
                if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                  // In case of a partition level shift skip if the child is already a halo cell
                  if(haloCellSet[i].find(childId) != haloCellSet[i].end()) {
                    continue;
                  }
                }

                m_haloCells[i].push_back(childId);
                haloCellSet[i].insert(childId);
              } else if(m_maxPartitionLevelShift > 0) {
                //if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && a_hasProperty(cellId, Cell::IsPeriodic))
                //  mTerm(1, AT_, "check code here!!!");
#ifndef NDEBUG
                if (ndom==domainId() && a_hasProperty(cellId, Cell::IsPeriodic))
                  TERMM_IF_NOT_COND(a_hasProperty(cellId, Cell::IsPartLvlAncestor), "");
#endif
                // Here also the case ndom==domainId(), where child halo lies on current domain because of periodicty
                // is dealt with
                if(ndom != domainId() || a_hasProperty(cellId, Cell::IsPeriodic)) {
                  const MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells, m_windowLayer_);
                  ASSERT(idx > -1, "");

                  // New neighbor domain, create new cell sets
                  if(idx == static_cast<MInt>(windCellSet.size())) {
                    windCellSet.emplace_back();
                    haloCellSet.emplace_back();
                  }

                  if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                    // In case of a partition level shift skip if the child is already a halo cell
                    if(haloCellSet[idx].find(childId) != haloCellSet[idx].end()) {
                      continue;
                    }
                  }

                  m_haloCells[idx].push_back(childId);
                  haloCellSet[idx].insert(childId);
                }
              }
              a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
            }
            if(a_noChildren(cellId) == 0) mTerm(1, AT_, "all children gone");
          }
        }
      }
    }

    // also create level+1 halo cells for local partitionLevelAncestors with children on different domains
    if(m_maxPartitionLevelShift > 0) {
      for(MUint i = 0; i < m_partitionLevelAncestorIds.size(); i++) {
        const MInt cellId = m_partitionLevelAncestorIds[i];

        ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor),
               "not a partition level ancestor: cellId" + std::to_string(cellId) + " halo"
                   + std::to_string(a_hasProperty(cellId, Cell::IsHalo)) + " domain" + std::to_string(domainId()));

        if(a_hasProperty(cellId, Cell::IsHalo)) continue;
        ASSERT(findNeighborDomainId(a_globalId(cellId)) == domainId(), "");
        if(a_level(cellId) == level) {
          MLong* childIds = &m_partitionLevelAncestorChildIds[i * m_maxNoChilds];
          if(accumulate(childIds, childIds + m_maxNoChilds, static_cast<MLong>(-1),
                        [](const MLong& a, const MLong& b) { return std::max(a, b); })
             > -1) {
            refineCell(cellId, childIds, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
            refineIds.push_back(cellId);

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;
              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              a_globalId(childId) = childIds[child];
              MInt ndom = findNeighborDomainId(a_globalId(childId));
              ASSERT(ndom > -1, "");
              if(ndom != domainId()) {
                MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells, m_windowLayer_);
                ASSERT(idx > -1, "");

                // New neighbor domain, create new cell sets
                if(idx == static_cast<MInt>(windCellSet.size())) {
                  windCellSet.emplace_back();
                  haloCellSet.emplace_back();
                }

                if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                  // In case of a partition level shift skip if the child is already a halo cell
                  if(haloCellSet[idx].find(childId) != haloCellSet[idx].end()) {
                    continue;
                  }
                }

                m_haloCells[idx].push_back(childId);
                haloCellSet[idx].insert(childId);
              }
              a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
            }
          }
          if(a_noChildren(cellId) == 0) mTerm(1, AT_, "all children gone");
        }
      }
    }

    // Refinement of the azimuthal periodic halo cells
    if(m_azimuthalPer) {
      refineChildIds.clear();
      recvChildIds.clear();

      // In this function cells are tagged for refinement
      tagAzimuthalHigherLevelExchangeCells(refineChildIds, recvChildIds, level);

      vector<vector<MLong>> shiftWindows;
      shiftWindows.resize(noDomains());
      MIntScratchSpace noShiftWindows(noNeighborDomains(), AT_, "noShiftWindows");
      noShiftWindows.fill(0);

      MInt cnt = 0;
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldWindowVec(m_azimuthalWindowCells[i]);
	vector<MInt> oldHigherLevelCon(m_azimuthalHigherLevelConnectivity[i]);

        m_azimuthalWindowCells[i].clear();
	m_azimuthalHigherLevelConnectivity[i].clear();

        for(MInt j = 0; j < (signed)oldWindowVec.size(); j++) {
          MInt cellId = oldWindowVec[j];

	  if(find(oldHigherLevelCon.begin(), oldHigherLevelCon.end(), j) != oldHigherLevelCon.end()) {
	    m_azimuthalHigherLevelConnectivity[i].push_back(m_azimuthalWindowCells[i].size());
	  }	  
          m_azimuthalWindowCells[i].push_back(cellId);

          ASSERT(cellId < m_tree.size(), to_string(level));

          if(a_level(cellId) == level) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MLong globalChildId = refineChildIds[m_maxNoChilds * cnt + child];
              // If the correspondig halo cell will be refiened, add child to azimuthal window cell vector
              if(globalChildId > -1) {
                MInt ndom = findNeighborDomainId(globalChildId);
                // Since azimuthal halo cells and azimuthal window cells are not complete copies
                // in the sense of the cartesian grid, the corresponding internal cells for the children of
                // an azimuthal halo cell might not lie in the internal cell which is mapped to the parent
                // halo cell. Therefore, it might be a child of an azimuthal halo cell is connected to a
                // differnt mpi rank. This is handled by the shiftWindows.
                if(ndom == domainId()) {
                  const MInt childId = globalIdToLocalId(globalChildId, true);
		  if(level == m_minLevel && a_level(childId) <= m_minLevel) {
		    m_azimuthalHigherLevelConnectivity[i].push_back(m_azimuthalWindowCells[i].size());
		  }
                  m_azimuthalWindowCells[i].push_back(childId);
                } else {
                  ASSERT(m_nghbrDomainIndex[ndom] > -1, "");
                  shiftWindows[m_nghbrDomainIndex[ndom]].push_back(globalChildId);
                  shiftWindows[m_nghbrDomainIndex[ndom]].push_back(azimuthalNeighborDomain(i));
                  noShiftWindows[m_nghbrDomainIndex[ndom]]++;
                }
              }
            }
          }
          cnt++;
        }
      }
      
      vector<vector<MInt>> shiftHalos;
      shiftHalos.resize(noDomains());
      vector<vector<MInt>> shiftHaloDoms;
      shiftHaloDoms.resize(noDomains());
      MIntScratchSpace noShiftHalos(noDomains(), AT_, "noShiftHalos");
      noShiftHalos.fill(0);
      cnt = 0;
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldHaloVec(m_azimuthalHaloCells[i]);

        m_azimuthalHaloCells[i].clear();

        for(MInt j = 0; j < (signed)oldHaloVec.size(); j++) {
          MInt cellId = oldHaloVec[j];
        
          m_azimuthalHaloCells[i].push_back(cellId);

          if(a_level(cellId) == level) {
            if(accumulate(&recvChildIds[m_maxNoChilds * cnt], &recvChildIds[m_maxNoChilds * cnt] + m_maxNoChilds,
                          static_cast<MLong>(-2), [](const MLong& a, const MLong& b) { return std::max(a, b); })
               > -2) {

              // Refine the azimuthal halo cell
              refineCell(cellId, nullptr, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
              refineIds.push_back(cellId);

              for(MInt child = 0; child < m_maxNoChilds; child++) {
                const MInt childId = a_childId(cellId, child);
                ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
                a_hasProperty(childId, Cell::IsHalo) = a_hasProperty(cellId, Cell::IsHalo);
                a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
                a_globalId(childId) = recvChildIds[m_maxNoChilds * cnt + child];
                for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
                  treeb().solver(childId, solver) = treeb().solver(cellId, solver);
                }

                // If the globalId of the child == -1, it has no corresponding interal window cell.
                // Since it still might be required is becomes an unmapped halo cell
                if(a_globalId(childId) == -1) {
                  m_azimuthalUnmappedHaloCells.push_back(childId);
                  m_azimuthalUnmappedHaloDomains.push_back(azimuthalNeighborDomain(i));
                  continue;
                }

                ASSERT(a_globalId(childId) > -1,"");
                MInt ndom = findNeighborDomainId(a_globalId(childId));

                // Since azimuthal halo cells and azimuthal window cells are not complete copies
                // in the sense of the cartesian grid, the corresponding internal cells for the children of
                // an azimuthal halo cell might not lie in the internal cell which is mapped to the parent
                // halo cell. Therefore, it might be a child of an azimuthal halo cell is connected to a
                // differnt mpi rank. This is handled by the shiftHalos.
                if(ndom == azimuthalNeighborDomain(i)) {
                  m_azimuthalHaloCells[i].push_back(childId);
                } else {
                  shiftHalos[ndom].push_back(childId);
                  shiftHaloDoms[ndom].push_back(azimuthalNeighborDomain(i));
                  noShiftHalos[ndom]++;
                }
              }
            }
          }
          cnt++;
        }
      }

      // Finally, the shiftWindows and shiftHalos are mapped.
      vector<vector<MLong>> newWindows;
      newWindows.resize(noNeighborDomains());
      MIntScratchSpace noNewWindows(noNeighborDomains(), AT_, "noNewWindows");
      noNewWindows.fill(0);
      vector<MInt> noValsToReceive;
      noValsToReceive.assign(noNeighborDomains(), 1);
      maia::mpi::exchangeBuffer(m_nghbrDomains, noValsToReceive, mpiComm(), noShiftWindows.data(), noNewWindows.getPointer());

    
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        newWindows[i].resize(2*noNewWindows[i]);
      }

      ScratchSpace<MPI_Request> sendReq(noNeighborDomains(), AT_, "sendReq");
      sendReq.fill(MPI_REQUEST_NULL);

      const MInt magic_mpi_tag = 12;
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noShiftWindows[i] == 0) { continue;
}
        MPI_Issend(&shiftWindows[i][0], 2 * noShiftWindows[i], type_traits<MLong>::mpiType(), neighborDomain(i), magic_mpi_tag, mpiComm(), &sendReq[i], AT_, "shiftWindows[i][0]");
      }
      
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noNewWindows[i] == 0) { continue;
}
        MPI_Recv(&newWindows[i][0], 2 * noNewWindows[i], type_traits<MLong>::mpiType(), neighborDomain(i), magic_mpi_tag, mpiComm(), MPI_STATUS_IGNORE, AT_, "newWindows[i][0]");
      }

      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noShiftWindows[i] == 0) { continue;
}
        MPI_Wait(&sendReq[i], MPI_STATUSES_IGNORE, AT_);
      }

      // Sort windows - optimize as in updateHaloCellCollectors?
      for(MInt i = 0; i < noDomains(); i++) {
	shiftWindows[i].clear();
	for(MInt d = 0; d < noDomains(); d++) {
	  MInt nghbrDomId = m_nghbrDomainIndex[d];
	  if(nghbrDomId == -1) { continue;
}
	  for(MInt j = 0; j < noNewWindows[nghbrDomId]; j++ ) {
	    MInt haloDomain = newWindows[nghbrDomId][j*2 + 1];
            if(haloDomain == i) { shiftWindows[i].push_back(newWindows[nghbrDomId][j*2]);
}
	  }
	}
      }

      // New windows
      for(MInt i = 0; i < noDomains(); i++) {
        for(MUint j = 0; j < shiftWindows[i].size(); j++ ) {
          const MLong globalChildId = shiftWindows[i][j];
          ASSERT(findNeighborDomainId(globalChildId) == domainId(),"Wrong domain!");

          const MInt childId = globalIdToLocalId(globalChildId, true);
          MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);

	  if(level == m_minLevel && a_level(childId) <= m_minLevel) {
	    m_azimuthalHigherLevelConnectivity[idx].push_back(m_azimuthalWindowCells[idx].size());
	  }
          m_azimuthalWindowCells[idx].push_back(childId);
        }
      }
      
      // Sort halos
      for(MInt i = 0; i < noDomains(); i++) {
        vector<MInt> shiftHalosBck(shiftHalos[i]);
        shiftHalos[i].clear();
        for(MInt d = 0; d < noDomains(); d++) {
          for(MInt j = 0; j < noShiftHalos[i]; j++ ) {
            MInt nghbrDomain = shiftHaloDoms[i][j];

            if(nghbrDomain == d) { shiftHalos[i].push_back(shiftHalosBck[j]);
}
          }
        }
      }

      // Create new halos
      for(MInt i = 0; i < noDomains(); i++) {
        for(MInt j = 0; j < noShiftHalos[i]; j++ ) {
          const MInt childId = shiftHalos[i][j];
          MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);
          
          m_azimuthalHaloCells[idx].push_back(childId);
        }
      }

      // Refine unmapped halos
      if(!refineCellSolver.empty()) {
	if(domainId() == 0) { cerr << "Attention: Unmapped halo are not refined during adaptation!" << endl;
}
      } else {
	shiftWindows.clear();
	recvChildIds.clear();
	vector<MInt> recvChildDomainIds;
	tagAzimuthalUnmappedHaloCells(shiftWindows, recvChildIds, recvChildDomainIds, level);

	// Newly mapped windows
        // If a child of an unmapped azimuthal halo cell has an adequate internal cell
        // This child becomes an regular azimuthal halo cell. Thus, also the internal window cell
        // is added to the azimuthal window cell vector
	for(MInt i = 0; i < noDomains(); i++) {
	  for(MUint j = 0; j < shiftWindows[i].size(); j++ ) {
	    const MLong globalChildId = shiftWindows[i][j];
	    if(globalChildId < 0) { continue;
}

	    ASSERT(findNeighborDomainId(globalChildId) == domainId(),"Wrong domain! " + to_string(globalChildId) + " " + to_string(findNeighborDomainId(globalChildId)) + " " + to_string(domainId()));

	    MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);

	    const MInt childId = globalIdToLocalId(globalChildId, true);
	    m_azimuthalWindowCells[idx].push_back(childId);
	  }
	}

	// Newly mapped halos
        // If a child of an unmapped azimuthal halo cell has an adequate internal cell
        // This child becomes an regular azimuthal halo cell. Thus, it is added to the azimuthal halo cell vector
        // Otherwise the child becomes an unmapped azimuthal halo cell
	MInt noUnmappedCells = noAzimuthalUnmappedHaloCells();
	for(MInt i = 0; i < noUnmappedCells; i++ ) {
	  MInt cellId = azimuthalUnmappedHaloCell(i);
	  if(a_level(cellId) == level) {
	    if(accumulate(&recvChildIds[m_maxNoChilds * i], &recvChildIds[m_maxNoChilds * i] + m_maxNoChilds,
                          static_cast<MLong>(-2), [](const MLong& a, const MLong& b) { return std::max(a, b); })
               > -2) {

	      refineCell(cellId, nullptr, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
	      refineIds.push_back(cellId);

	      for(MInt child = 0; child < m_maxNoChilds; child++) {
		const MInt childId = a_childId(cellId, child);

		ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
		a_hasProperty(childId, Cell::IsHalo) = a_hasProperty(cellId, Cell::IsHalo);
		a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
		a_globalId(childId) = recvChildIds[m_maxNoChilds * i + child];
		for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
		  treeb().solver(childId, solver) = treeb().solver(cellId, solver);
		}

		MInt dom = recvChildDomainIds[m_maxNoChilds * i + child];

		if(a_globalId(childId) == -1) {
		  m_azimuthalUnmappedHaloCells.push_back(childId);
		  m_azimuthalUnmappedHaloDomains.push_back(dom);
		  continue;
		}
              
		ASSERT(findNeighborDomainId(a_globalId(childId)) == dom, "Wrong domain! " + to_string(a_globalId(childId)) + " " + to_string(findNeighborDomainId(a_globalId(childId))) + " " + to_string(dom) + " " + to_string(domainId()) + " " + to_string(cellId));

		MInt idx = setAzimuthalNeighborDomainIndex(dom, m_azimuthalHaloCells, m_azimuthalWindowCells);

		m_azimuthalHaloCells[idx].push_back(childId);
	      }
	    }
	  }
        }
      }

      // Update grid bndry cells
      MBool adaptation = false;
      adaptation = Context::getBasicProperty<MBool>("adaptation", AT_, &adaptation);
      if(!adaptation) {
        std::array<MFloat, 2*nDim> bbox;
	for(MInt i = 0; i < nDim; i++) {
	  bbox[i] = numeric_limits<MFloat>::max();
	  bbox[nDim + i] = numeric_limits<MFloat>::lowest();
	}
	for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
	  for(MInt i = 0; i < nDim; i++) {
	    bbox[i] = mMin(bbox[i], a_coordinate(cellId, i) - F1B2 * cellLengthAtLevel(a_level(cellId)));
	    bbox[nDim + i] = mMax(bbox[nDim + i], a_coordinate(cellId, i) + F1B2 * cellLengthAtLevel(a_level(cellId)));
	  }
	}
	MPI_Allreduce(MPI_IN_PLACE, &bbox[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[0]");
	MPI_Allreduce(MPI_IN_PLACE, &bbox[nDim], nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[nDim]");
	MBool isBndry = false;
	for(MInt i = 0; i < m_noInternalCells; i++) {
	  MInt cellId =i;  
	  if(a_level(cellId) != level+1) { continue;
}
	  isBndry = false;
	  for(MInt dir = 0; dir < m_noDirs; dir++) {
	    if(a_hasNeighbor(cellId, dir) > 0) { continue;
}
	    if(a_hasParent(cellId) && a_hasNeighbor(a_parentId(cellId), dir) > 0) {
	      MInt nghbrParentId = a_neighborId(a_parentId(cellId),dir);
	      if(a_noChildren(nghbrParentId) == 0) {
		continue;
	      }
	    }
	    if(!m_periodicCartesianDir[dir/2]) {
              std::array<MFloat, nDim> coords;
	      for(MInt d = 0; d < nDim; d++) {
		coords[d] = a_coordinate(cellId,d);
	      }
	      coords[dir / 2] += ((dir % 2) == 0 ? -F1 : F1) * cellLengthAtLevel(m_minLevel);
	      if(coords[dir / 2] < bbox[dir / 2] || coords[dir / 2] > bbox[nDim + dir / 2]) {
		continue;
}
	    }
	    isBndry = true;
	  }
	  if(isBndry) {
	    if(!a_isHalo(cellId)) { m_gridBndryCells.push_back(cellId);
}
	    for(MInt dir = 0; dir < m_noDirs; dir++) {
	      if(a_hasNeighbor(cellId, dir) > 0) {
		MInt nghbrId = a_neighborId(cellId, dir);
		if(!a_isHalo(nghbrId)) {
		  m_gridBndryCells.push_back(nghbrId);
		}
	      }
	    }
	  }
	}
      }
    }

    // sort halo and window cells by globalId to get matching connectivity
    // this may not even be required, check whether ordering is implicity matching for
    // a complex case with multiple partition level shifts
    if(m_maxPartitionLevelShift > 0) {
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        sort(m_haloCells[i].begin(), m_haloCells[i].end(),
             [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
        sort(m_windowCells[i].begin(), m_windowCells[i].end(),
             [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
      }
    }

    // Exchange solverBits also in case of 1 solver, because solverBits are used to disable halo cells, which
    // are created, but in tagActiveWindowsAtMinLevel seen to be superfluous; if we can ensure that in case of
    // noSolvers==0, all cells in m_windowCells have a windowLayer<=m_noSolverHaloLayers, we might recomment the
    // following if-statement
    //if(treeb().noSolvers() > 1) { // Exchange solver info
      exchangeSolverBitset(&m_tree.solverBits(0));
    //}
    if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
      maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(), &m_tree.solverBits(0),
				m_tree.size());

      // To ensure that azimuthal halo cells are fully refined (have all children)
      // or are not refined at all, solver Bits need to be checked!
      correctAzimuthalSolverBits();
    }

    // trigger refineCell() for halo cells on the solvers
    if(!refineCellSolver.empty()) {
      // initially empty during solver startup, no solver refinement needed then.
      for(auto& cellId : refineIds) {
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) {
            refineCellSolver[solver](cellId); // call refineCell() function in each solver
          }
        }
      }
      refineIds.clear();
    }

    // When we arrive at the highest level, we can optionally delete uneccesary cells
    if (m_haloMode==2 && (forceLeafLvlCorrection || (onlyLevel==-1 && level==m_maxLevel-1))) {
      tagActiveWindowsOnLeafLvl3(level+1, duringMeshAdaptation, removeCellSolver);
    }

    // --- Exchange Cell::IsPartLvlAncestor & noOffsprings --- //
    // exchange some window cell data
    ScratchSpace<MInt> bprops(m_tree.size(), 2, AT_, "bprops");
    bprops.fill(-1);
    // Note: changed loop to iterate over whole tree, internal and halo cells not sorted!
    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      if(a_isHalo(cellId) || a_isToDelete(cellId)) {
        continue;
      }
      bprops(cellId, 0) = (MInt)a_hasProperty(cellId, Cell::IsPartLvlAncestor);
      bprops(cellId, 1) = a_noOffsprings(cellId);

      ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
             std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                 + std::to_string(bprops(cellId, 1)));
    }

    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), bprops.getPointer(), 2);

    // Note: loop only over the current haloCells since only for these data is exchanged, this does
    // not work properly if in case of a partition level shift some halo cell is not part of
    // m_haloCells yet such that its property IsPartLvlAncestor is reset!
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        const MInt cellId = m_haloCells[i][j];
        ASSERT(!a_isToDelete(cellId), "");
        ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
               std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                   + std::to_string(bprops(cellId, 1)));
        a_hasProperty(cellId, Cell::IsPartLvlAncestor) = (MBool)bprops(cellId, 0);
        a_noOffsprings(cellId) = bprops(cellId, 1);
      }
    }
    // --- Exchange Cell::IsPartLvlAncestor & noOffsprings --- //

    // Set window/halo flags
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      // Halo flag are already set in refineCell
#ifndef NDEBUG
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        ASSERT(!a_hasProperty(m_haloCells[i][j], Cell::IsWindow), "halo cell is marked as window");
        ASSERT(a_hasProperty(m_haloCells[i][j], Cell::IsHalo), "halo cell flag not set!");
        // a_hasProperty(m_haloCells[i][j], Cell::IsHalo) = true;
        // a_hasProperty(m_haloCells[i][j], Cell::IsWindow) = false;
      }
#endif
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        ASSERT(!a_isHalo(m_windowCells[i][j]), "window cell is marked as halo");
        // a_hasProperty( m_windowCells[i][j], Cell::IsHalo ) = false;
        a_hasProperty(m_windowCells[i][j], Cell::IsWindow) = true;
      }
    }
    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
          ASSERT(!a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsWindow), "azimuthal halo cell is marked as window");
          a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsHalo) = true;
        }
        for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
          ASSERT(!a_isHalo(m_azimuthalWindowCells[i][j]), "azimuthal window cell is marked as halo");
          a_hasProperty(m_azimuthalWindowCells[i][j], Cell::IsWindow) = true;
        }
      }
      for(MInt j = 0; j < noAzimuthalUnmappedHaloCells(); j++ ) {
        ASSERT(!a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsWindow), "azimuthal halo cell is marked as window");
      a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsHalo) = true;
      }
    }

#ifndef NDEBUG
    // Sanity check, that parent of a halo cell is also a halo cell
    for (MInt i = 0; i < noNeighborDomains(); i++) {
      // We need to reinit set, because in tagActiveWindowsOnLeafLvl3 the ordering might change 
      haloCellSet[i].clear();
      std::copy(m_haloCells[i].begin(), m_haloCells[i].end(), std::inserter(haloCellSet[i], haloCellSet[i].begin()));
      for (MInt j = 0; j < (signed)m_haloCells[i].size(); ++j) {
        MInt parentId = m_haloCells[i][j];
        TERMM_IF_NOT_COND(a_isHalo(parentId), "");
        char isSolverHalo = 0;
        for (MInt solver = 0; solver < treeb().noSolvers(); ++solver) {
          if (m_tree.solver(parentId,solver)) {
            isSolverHalo = isSolverHalo | (1 << solver);
          }
        }
        while(a_parentId(parentId)>-1) {
          parentId = a_parentId(parentId);
          TERMM_IF_NOT_COND(a_isHalo(parentId) || a_hasProperty(parentId, Cell::IsPartLvlAncestor), "");
          TERMM_IF_NOT_COND(haloCellSet[i].find(parentId)!=haloCellSet[i].end()
              || a_hasProperty(parentId, Cell::IsPartLvlAncestor), "");
          for (MInt solver = 0; solver < treeb().noSolvers(); ++solver) {
            if (isSolverHalo & (1<<solver)) {
              TERMM_IF_NOT_COND(m_tree.solver(parentId,solver), "");
            }
            if (m_tree.solver(parentId,solver)) {
              isSolverHalo = isSolverHalo | (1 << solver);
            }
          }
        }
      }
    }
    // for debugging, check consistency before new level
    checkWindowHaloConsistency(false, "createHigherLevelExchangeCells completed on level="+to_string(level)+": ");
#endif

    logDuration(levelTimeStart, "level #"+std::to_string(level)+" total");
  }

#ifndef NDEBUG
  if (onlyLevel==-1 || forceLeafLvlCorrection)
    checkWindowLayer("createHigherLevelExchangeCells completed: ");
#endif

  if(domainId() == 0 && onlyLevel < 0) cerr << endl;
}
//clang-format on

//WH_OLD
/** \brief Iteratively create window and halo cells on levels m_minLevel+1...m_maxLevel,
 *         starting from m_minLevel exchange cells
    \author Lennart Schneiders
    \date October 2017
  */

template <MInt nDim>
void CartesianGrid<nDim>::createHigherLevelExchangeCells_old(
    const MInt onlyLevel, const std::vector<std::function<void(const MInt)>>& refineCellSolver) {
  TRACE();

  ScratchSpace<MInt> plaMap(m_maxNoCells, AT_, "plaMap");
  vector<MLong> refineChildIds;
  vector<MLong> recvChildIds;
  vector<vector<MInt>> haloMap;
  vector<vector<MInt>> windMap;
  vector<MInt> refineIds;

  if(domainId() == 0 && onlyLevel < 0) cerr << "  * create higher level exchange cells...";

  plaMap.fill(-1);
  if(m_maxPartitionLevelShift > 0) {
    for(MUint i = 0; i < m_partitionLevelAncestorIds.size(); i++) {
      plaMap[m_partitionLevelAncestorIds[i]] = i;
    }
  }

  if(m_azimuthalPer && onlyLevel <= m_minLevel) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      m_azimuthalHigherLevelConnectivity[i].clear();
    }
  }

  // given window-halo cell connectivity on minLevel established in setupWindowHaloCellConnectivity(),
  // iteratively create connectivity on ascending levels
  for(MInt level = m_minLevel; level < m_maxLevel; level++) {
    if(onlyLevel > -1 && onlyLevel != level) continue;

    // for debugging, check consistency before new level
    // checkWindowHaloConsistency();

    refineIds.clear();
    MInt haloCnt = 0;
    MInt windowCnt = 0;
    haloMap.resize(noNeighborDomains());
    windMap.resize(noNeighborDomains());

    // Sets to hold halo/window cells for each neighbor domain to allow fast searching if cells are
    // already present
    vector<std::set<MInt>> haloCellSet(noNeighborDomains());
    vector<std::set<MInt>> windCellSet(noNeighborDomains());

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      haloMap[i].resize(m_haloCells[i].size());
      windMap[i].resize(m_windowCells[i].size());
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        haloMap[i][j] = haloCnt++;
      }
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        windMap[i][j] = windowCnt++;
      }

      haloCellSet[i].insert(m_haloCells[i].begin(), m_haloCells[i].end());
      windCellSet[i].insert(m_windowCells[i].begin(), m_windowCells[i].end());
    }

    recvChildIds.resize(haloCnt * m_maxNoChilds);
    fill(recvChildIds.begin(), recvChildIds.end(), -1);
    refineChildIds.clear();
    tagActiveWindows(refineChildIds, level);


// set globalIds of each windowCells's child
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        MInt cellId = m_windowCells[i][j];

        if(a_level(cellId) == level && a_noChildren(cellId) > 0) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt cnt = windMap[i][j];

            if(refineChildIds[m_maxNoChilds * cnt + child] > 0) {
              ASSERT(a_childId(cellId, child) > -1, "");
              ASSERT((a_globalId(static_cast<MInt>(a_childId(cellId, child))) >= m_domainOffsets[domainId()]
                      && a_globalId(static_cast<MInt>(a_childId(cellId, child))) < m_domainOffsets[domainId() + 1])
                         || a_hasProperty(cellId, Cell::IsPartLvlAncestor),
                     to_string(a_globalId(static_cast<MInt>(a_childId(cellId, child)))) + " "
                         + to_string(m_domainOffsets[domainId()]) + " " + to_string(m_domainOffsets[domainId() + 1]));
              refineChildIds[m_maxNoChilds * cnt + child] = a_globalId(static_cast<MInt>(a_childId(cellId, child)));
            }

            MInt pId = plaMap[cellId];
            if(m_maxPartitionLevelShift > 0 && pId > -1) {
              if(m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child] > -1) {
                ASSERT(refineChildIds[m_maxNoChilds * cnt + child] < 0
                           || refineChildIds[m_maxNoChilds * cnt + child]
                                  == m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child],
                       std::to_string(refineChildIds[m_maxNoChilds * cnt + child]) + " plac"
                           + std::to_string(m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child]) + "; c"
                           + std::to_string(cellId) + "; pId" + std::to_string(pId) + "; cnt" + std::to_string(cnt));

                refineChildIds[m_maxNoChilds * cnt + child] =
                    m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child];
              }
            }
          }
        }
      }
    }

    // exchange childIds
    maia::mpi::exchangeBuffer(m_nghbrDomains,
                              m_haloCells,
                              m_windowCells,
                              mpiComm(),
                              refineChildIds.data(),
                              recvChildIds.data(),
                              m_maxNoChilds);

#ifndef NDEBUG
    if(m_maxPartitionLevelShift > 0) {
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
          MInt cellId = m_haloCells[i][j];
          if(a_level(cellId) == level) {
            MInt pId = plaMap[cellId];
            if(m_maxPartitionLevelShift > 0 && pId > -1) {
              MInt cnt = haloMap[i][j];
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                const MLong childId = m_partitionLevelAncestorChildIds[pId * m_maxNoChilds + child];
                ASSERT(childId == recvChildIds[m_maxNoChilds * cnt + child],
                       std::to_string(childId) + " != " + std::to_string(recvChildIds[m_maxNoChilds * cnt + child])
                           + "; c" + std::to_string(cellId) + "; g" + std::to_string(a_globalId(cellId)) + "; p"
                           + std::to_string(pId));
              }
            }
          }
        }
      }
    }
#endif


    // initiate communication between two other domains, if one of my childs lies on a different domain (due to
    // partitionLevelShift)
    if(m_maxPartitionLevelShift > 0) {
      // determine all pairs of neighbor domains to connect
      vector<MInt> sendDomIdx;
      vector<MLong> sendData;
      // Note: this assumed that internal and halo cells are separated which is not true if called from meshAdaptation!
      // MBoolScratchSpace cellFlag(m_noInternalCells, AT_, "cellFlag");
      MBoolScratchSpace cellFlag(m_tree.size(), AT_, "cellFlag");
      fill(cellFlag.begin(), cellFlag.end(), false);

      for(MInt i = 0; i < noNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
          const MInt cellId = m_windowCells[i][j];
          ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "window cell marked as halo");

          if(a_level(cellId) == level && a_noChildren(cellId) > 0) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MLong childId = refineChildIds[m_maxNoChilds * windMap[i][j] + child];
              if(childId < 0) continue;
              MInt ndom = findNeighborDomainId(childId);
              ASSERT(ndom > -1 && ndom < noDomains(), "");
              if(ndom != domainId()) {
                ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor), "");
                MInt idx = findIndex(m_nghbrDomains.begin(), m_nghbrDomains.end(), ndom);
                ASSERT(idx > -1, "");
                ASSERT(m_nghbrDomains[idx] == ndom, "");
                if(ndom != m_nghbrDomains[i]) {
                  sendDomIdx.push_back(idx);
                  sendData.push_back(m_nghbrDomains[i]);
                  sendData.push_back(childId);
                }
                if(!cellFlag[cellId]) {
                  sendDomIdx.push_back(idx);
                  sendData.push_back(domainId());
                  sendData.push_back(childId);
                }
              }
            }
          }
          cellFlag[cellId] = true;
        }
      }

      // exchange redirected communication info
      vector<MInt> recvOffs;
      vector<MLong> recvData;
      maia::mpi::exchangeScattered(m_nghbrDomains, sendDomIdx, sendData, mpiComm(), recvOffs, recvData, 2);

      // create window cells if this domain is affected by redirected communication links
      const MInt noNgbrDomainsBak = noNeighborDomains();
      for(MInt i = 0; i < noNgbrDomainsBak; i++) {
        for(MInt j = recvOffs[i]; j < recvOffs[i + 1]; j++) {
          auto ndom = static_cast<MInt>(recvData[2 * j]);
          MLong globalChildId = recvData[2 * j + 1];
          if(ndom > -1) {
            MInt childId = m_globalToLocalId[globalChildId];
            ASSERT(childId > -1, "");
            ASSERT(a_globalId(childId) == globalChildId, "");
            ASSERT(findNeighborDomainId(globalChildId) == domainId(), "");
            MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells);
            ASSERT(idx > -1, "");

            // New neighbor domain, create new cell sets
            if(idx == static_cast<MInt>(windCellSet.size())) {
              windCellSet.emplace_back();
              haloCellSet.emplace_back();
            }

            if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
              // In case of a partition level shift skip if the child is already a window cell
              if(windCellSet[idx].find(childId) != windCellSet[idx].end()) {
                continue;
              }
            }

            ASSERT(windCellSet[idx].find(childId) == windCellSet[idx].end(),
                   "duplicate window cell: " + std::to_string(childId));

            m_windowCells[idx].push_back(childId);
            windCellSet[idx].insert(childId);
          }
        }
      }

    } // ( m_maxPartitionLevelShift > 0 )


    // update regular window cells between this domain and m_nghbrDomains[i] to include the level+1 cells
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      vector<MInt> oldWindowVec(m_windowCells[i]);
      std::set<MInt> oldWindCellSet(windCellSet[i]);

      std::vector<MInt>().swap(m_windowCells[i]);
      windCellSet[i].clear();

      for(MInt j = 0; j < (signed)oldWindowVec.size(); j++) {
        MInt cellId = oldWindowVec[j];

        // TODO labels:GRID fix duplicate window cells in periodic cases!
        ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet[i].find(cellId) == windCellSet[i].end(),
               "duplicate window cell: " + std::to_string(cellId) + " " + std::to_string(a_globalId(cellId)) + " "
                   + std::to_string(j));

        m_windowCells[i].push_back(cellId);
        windCellSet[i].insert(cellId);

        ASSERT(cellId < m_tree.size(), to_string(level));

        if(a_level(cellId) == level && a_noChildren(cellId) > 0) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt cnt = windMap[i][j];

            const MInt childId = a_childId(cellId, child);
            if(childId < 0) continue;

            // Skip child if it is already part of the old window cell vector and is handled by the
            // outer loop
            if(oldWindCellSet.find(childId) != oldWindCellSet.end()) {
              ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                     "not a partition level ancestor or partition cell");
              continue;
            }

            const MLong globalChildId = refineChildIds[m_maxNoChilds * cnt + child];
            if(globalChildId > -1) {
              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              // ASSERT( childId > -1 && childId < m_noInternalCells, to_string(childId) + " " +
              // to_string(m_noInternalCells) );
              MInt ndom = findNeighborDomainId(globalChildId);
              if(ndom == domainId()) {
                ASSERT(globalChildId == a_globalId(childId), "");

                // TODO labels:GRID fix duplicate window cells in periodic cases!
                ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet[i].find(childId) == windCellSet[i].end(),
                       "duplicate window cell: " + std::to_string(childId));

                m_windowCells[i].push_back(childId);
                windCellSet[i].insert(childId);
              }
            }
          }
        }
      }
    }


    // create level+1 halo cells existing on other domains
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      vector<MInt> oldHaloVec(m_haloCells[i]);
      std::set<MInt> oldHaloCellSet(haloCellSet[i]);

      std::vector<MInt>().swap(m_haloCells[i]);
      haloCellSet[i].clear();

      for(MInt j = 0; j < (signed)oldHaloVec.size(); j++) {
        MInt cellId = oldHaloVec[j];

        ASSERT(haloCellSet[i].find(cellId) == haloCellSet[i].end(),
               "duplicate halo cell: " + std::to_string(cellId) + " " + std::to_string(a_globalId(cellId)) + " "
                   + std::to_string(j));

        m_haloCells[i].push_back(cellId);
        haloCellSet[i].insert(cellId);

        if(a_level(cellId) == level) {
          MInt cnt = haloMap[i][j];

          if(accumulate(&recvChildIds[m_maxNoChilds * cnt], &recvChildIds[m_maxNoChilds * cnt] + m_maxNoChilds,
                        static_cast<MLong>(-1), [](const MLong& a, const MLong& b) { return std::max(a, b); })
             > -1) {
            // NOTE: this leads to different noChilds on window/halo-Cells!
            // TODO labels:GRID,totest
            // check if halo cell, that is to be refined, is on the second periodic layer
            // if so, do not refine this cell
            // if(a_hasProperty(cellId,Cell::IsPeriodic)){

            // }

            refineCell(cellId, &recvChildIds[m_maxNoChilds * cnt], a_hasProperty(cellId, Cell::IsPartLvlAncestor));
            refineIds.push_back(cellId);

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;

              // Skip halo child already present in case of partition level shifts
              if(oldHaloCellSet.find(childId) != oldHaloCellSet.end()) {
                ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                       "not a partition level ancestor or partition cell");
                continue;
              }

              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              ASSERT(a_globalId(childId) > -1, "");
              MInt ndom = findNeighborDomainId(a_globalId(childId));
              ASSERT(ndom > -1, "");
              if(recvChildIds[m_maxNoChilds * cnt + child] > -1)
                ASSERT(a_globalId(childId) == recvChildIds[m_maxNoChilds * cnt + child], "");
              if(recvChildIds[m_maxNoChilds * cnt + child] < 0) ASSERT(ndom == domainId(), "");
              if(ndom == m_nghbrDomains[i]) {
                if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                  // In case of a partition level shift skip if the child is already a halo cell
                  if(haloCellSet[i].find(childId) != haloCellSet[i].end()) {
                    continue;
                  }
                }

                m_haloCells[i].push_back(childId);
                haloCellSet[i].insert(childId);
              } else if(m_maxPartitionLevelShift > 0) {
                if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && a_hasProperty(cellId, Cell::IsPeriodic))
                  mTerm(1, AT_, "check code here!!!");
                if(ndom != domainId()) {
                  MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells);
                  ASSERT(idx > -1, "");

                  // New neighbor domain, create new cell sets
                  if(idx == static_cast<MInt>(windCellSet.size())) {
                    windCellSet.emplace_back();
                    haloCellSet.emplace_back();
                  }

                  if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                    // In case of a partition level shift skip if the child is already a halo cell
                    if(haloCellSet[idx].find(childId) != haloCellSet[idx].end()) {
                      continue;
                    }
                  }

                  m_haloCells[idx].push_back(childId);
                  haloCellSet[idx].insert(childId);
                }
              }
              a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
            }
            if(a_noChildren(cellId) == 0) mTerm(1, AT_, "all children gone");
          }
        }
      }
    }

    // also create level+1 halo cells for local partitionLevelAncestors with children on different domains
    if(m_maxPartitionLevelShift > 0) {
      for(MUint i = 0; i < m_partitionLevelAncestorIds.size(); i++) {
        const MInt cellId = m_partitionLevelAncestorIds[i];

        ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor),
               "not a partition level ancestor: cellId" + std::to_string(cellId) + " halo"
                   + std::to_string(a_hasProperty(cellId, Cell::IsHalo)) + " domain" + std::to_string(domainId()));

        if(a_hasProperty(cellId, Cell::IsHalo)) continue;
        ASSERT(findNeighborDomainId(a_globalId(cellId)) == domainId(), "");
        if(a_level(cellId) == level) {
          MLong* childIds = &m_partitionLevelAncestorChildIds[i * m_maxNoChilds];
          if(accumulate(childIds, childIds + m_maxNoChilds, static_cast<MLong>(-1),
                        [](const MLong& a, const MLong& b) { return std::max(a, b); })
             > -1) {
            refineCell(cellId, childIds, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
            refineIds.push_back(cellId);

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;
              ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
              a_globalId(childId) = childIds[child];
              MInt ndom = findNeighborDomainId(a_globalId(childId));
              ASSERT(ndom > -1, "");
              if(ndom != domainId()) {
                MInt idx = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells);
                ASSERT(idx > -1, "");

                // New neighbor domain, create new cell sets
                if(idx == static_cast<MInt>(windCellSet.size())) {
                  windCellSet.emplace_back();
                  haloCellSet.emplace_back();
                }

                if(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell)) {
                  // In case of a partition level shift skip if the child is already a halo cell
                  if(haloCellSet[idx].find(childId) != haloCellSet[idx].end()) {
                    continue;
                  }
                }

                m_haloCells[idx].push_back(childId);
                haloCellSet[idx].insert(childId);
              }
              a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
            }
          }
          if(a_noChildren(cellId) == 0) mTerm(1, AT_, "all children gone");
        }
      }
    }

    // Refinement of the azimuthal periodic halo cells
    if(m_azimuthalPer) {
      refineChildIds.clear();
      recvChildIds.clear();

      // In this function cells are tagged for refinement
      tagAzimuthalHigherLevelExchangeCells(refineChildIds, recvChildIds, level);

      vector<vector<MLong>> shiftWindows;
      shiftWindows.resize(noDomains());
      MIntScratchSpace noShiftWindows(noNeighborDomains(), AT_, "noShiftWindows");
      noShiftWindows.fill(0);

      MInt cnt = 0;
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldWindowVec(m_azimuthalWindowCells[i]);
	vector<MInt> oldHigherLevelCon(m_azimuthalHigherLevelConnectivity[i]);

        m_azimuthalWindowCells[i].clear();
	m_azimuthalHigherLevelConnectivity[i].clear();

        for(MInt j = 0; j < (signed)oldWindowVec.size(); j++) {
          MInt cellId = oldWindowVec[j];

	  if(find(oldHigherLevelCon.begin(), oldHigherLevelCon.end(), j) != oldHigherLevelCon.end()) {
	    m_azimuthalHigherLevelConnectivity[i].push_back(m_azimuthalWindowCells[i].size());
	  }	  
          m_azimuthalWindowCells[i].push_back(cellId);

          ASSERT(cellId < m_tree.size(), to_string(level));

          if(a_level(cellId) == level) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MLong globalChildId = refineChildIds[m_maxNoChilds * cnt + child];
              // If the correspondig halo cell will be refiened, add child to azimuthal window cell vector
              if(globalChildId > -1) {
                MInt ndom = findNeighborDomainId(globalChildId);
                // Since azimuthal halo cells and azimuthal window cells are not complete copies
                // in the sense of the cartesian grid, the corresponding internal cells for the children of
                // an azimuthal halo cell might not lie in the internal cell which is mapped to the parent
                // halo cell. Therefore, it might be a child of an azimuthal halo cell is connected to a
                // differnt mpi rank. This is handled by the shiftWindows.
                if(ndom == domainId()) {
                  const MInt childId = globalIdToLocalId(globalChildId, true);
		  if(level == m_minLevel && a_level(childId) <= m_minLevel) {
		    m_azimuthalHigherLevelConnectivity[i].push_back(m_azimuthalWindowCells[i].size());
		  }
                  m_azimuthalWindowCells[i].push_back(childId);
                } else {
                  ASSERT(m_nghbrDomainIndex[ndom] > -1, "");
                  shiftWindows[m_nghbrDomainIndex[ndom]].push_back(globalChildId);
                  shiftWindows[m_nghbrDomainIndex[ndom]].push_back(azimuthalNeighborDomain(i));
                  noShiftWindows[m_nghbrDomainIndex[ndom]]++;
                }
              }
            }
          }
          cnt++;
        }
      }
      
      vector<vector<MInt>> shiftHalos;
      shiftHalos.resize(noDomains());
      vector<vector<MInt>> shiftHaloDoms;
      shiftHaloDoms.resize(noDomains());
      MIntScratchSpace noShiftHalos(noDomains(), AT_, "noShiftHalos");
      noShiftHalos.fill(0);
      cnt = 0;
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldHaloVec(m_azimuthalHaloCells[i]);

        m_azimuthalHaloCells[i].clear();

        for(MInt j = 0; j < (signed)oldHaloVec.size(); j++) {
          MInt cellId = oldHaloVec[j];
        
          m_azimuthalHaloCells[i].push_back(cellId);

          if(a_level(cellId) == level) {
            if(accumulate(&recvChildIds[m_maxNoChilds * cnt], &recvChildIds[m_maxNoChilds * cnt] + m_maxNoChilds,
                          static_cast<MLong>(-2), [](const MLong& a, const MLong& b) { return std::max(a, b); })
               > -2) {

              // Refine the azimuthal halo cell
              refineCell(cellId, nullptr, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
              refineIds.push_back(cellId);

              for(MInt child = 0; child < m_maxNoChilds; child++) {
                const MInt childId = a_childId(cellId, child);
                ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
                a_hasProperty(childId, Cell::IsHalo) = a_hasProperty(cellId, Cell::IsHalo);
                a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
                a_globalId(childId) = recvChildIds[m_maxNoChilds * cnt + child];
                for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
                  treeb().solver(childId, solver) = treeb().solver(cellId, solver);
                }

                // If the globalId of the child == -1, it has no corresponding interal window cell.
                // Since it still might be required is becomes an unmapped halo cell
                if(a_globalId(childId) == -1) {
                  m_azimuthalUnmappedHaloCells.push_back(childId);
                  m_azimuthalUnmappedHaloDomains.push_back(azimuthalNeighborDomain(i));
                  continue;
                }

                ASSERT(a_globalId(childId) > -1,"");
                MInt ndom = findNeighborDomainId(a_globalId(childId));

                // Since azimuthal halo cells and azimuthal window cells are not complete copies
                // in the sense of the cartesian grid, the corresponding internal cells for the children of
                // an azimuthal halo cell might not lie in the internal cell which is mapped to the parent
                // halo cell. Therefore, it might be a child of an azimuthal halo cell is connected to a
                // differnt mpi rank. This is handled by the shiftHalos.
                if(ndom == azimuthalNeighborDomain(i)) {
                  m_azimuthalHaloCells[i].push_back(childId);
                } else {
                  shiftHalos[ndom].push_back(childId);
                  shiftHaloDoms[ndom].push_back(azimuthalNeighborDomain(i));
                  noShiftHalos[ndom]++;
                }
              }
            }
          }
          cnt++;
        }
      }

      // Finally, the shiftWindows and shiftHalos are mapped.
      vector<vector<MLong>> newWindows;
      newWindows.resize(noNeighborDomains());
      MIntScratchSpace noNewWindows(noNeighborDomains(), AT_, "noNewWindows");
      noNewWindows.fill(0);
      vector<MInt> noValsToReceive;
      noValsToReceive.assign(noNeighborDomains(), 1);
      maia::mpi::exchangeBuffer(m_nghbrDomains, noValsToReceive, mpiComm(), noShiftWindows.data(), noNewWindows.getPointer());

    
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        newWindows[i].resize(2*noNewWindows[i]);
      }

      ScratchSpace<MPI_Request> sendReq(noNeighborDomains(), AT_, "sendReq");
      sendReq.fill(MPI_REQUEST_NULL);

      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noShiftWindows[i] == 0) continue;
        MPI_Issend(&shiftWindows[i][0], 2 * noShiftWindows[i], type_traits<MLong>::mpiType(), neighborDomain(i), 12, mpiComm(), &sendReq[i], AT_, "shiftWindows[i][0]");
      }
      
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noNewWindows[i] == 0) continue;
        MPI_Recv(&newWindows[i][0], 2 * noNewWindows[i], type_traits<MLong>::mpiType(), neighborDomain(i), 12, mpiComm(), MPI_STATUS_IGNORE, AT_, "newWindows[i][0]");
      }

      for(MInt i = 0; i < noNeighborDomains(); i++) {
        if(noShiftWindows[i] == 0) continue;
        MPI_Wait(&sendReq[i], MPI_STATUSES_IGNORE, AT_);
      }

      // Sort windows - optimize as in updateHaloCellCollectors?
      for(MInt i = 0; i < noDomains(); i++) {
	shiftWindows[i].clear();
	for(MInt d = 0; d < noDomains(); d++) {
	  MInt nghbrDomId = m_nghbrDomainIndex[d];
	  if(nghbrDomId == -1) continue;
	  for(MInt j = 0; j < noNewWindows[nghbrDomId]; j++ ) {
	    MInt haloDomain = newWindows[nghbrDomId][j*2 + 1];
            if(haloDomain == i) shiftWindows[i].push_back(newWindows[nghbrDomId][j*2]);
	  }
	}
      }

      // New windows
      for(MInt i = 0; i < noDomains(); i++) {
        for(MUint j = 0; j < shiftWindows[i].size(); j++ ) {
          const MLong globalChildId = shiftWindows[i][j];
          ASSERT(findNeighborDomainId(globalChildId) == domainId(),"Wrong domain!");

          const MInt childId = globalIdToLocalId(globalChildId, true);
          MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);

	  if(level == m_minLevel && a_level(childId) <= m_minLevel) {
	    m_azimuthalHigherLevelConnectivity[idx].push_back(m_azimuthalWindowCells[idx].size());
	  }
          m_azimuthalWindowCells[idx].push_back(childId);
        }
      }
      
      // Sort halos
      for(MInt i = 0; i < noDomains(); i++) {
        vector<MInt> shiftHalosBck(shiftHalos[i]);
        shiftHalos[i].clear();
        for(MInt d = 0; d < noDomains(); d++) {
          for(MInt j = 0; j < noShiftHalos[i]; j++ ) {
            MInt nghbrDomain = shiftHaloDoms[i][j];

            if(nghbrDomain == d) shiftHalos[i].push_back(shiftHalosBck[j]);
          }
        }
      }

      // Create new halos
      for(MInt i = 0; i < noDomains(); i++) {
        for(MInt j = 0; j < noShiftHalos[i]; j++ ) {
          const MInt childId = shiftHalos[i][j];
          MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);
          
          m_azimuthalHaloCells[idx].push_back(childId);
        }
      }

      // Refine unmapped halos
      if(!refineCellSolver.empty()) {
	if(domainId() == 0) cerr << "Attention: Unmapped halo are not refined during adaptation!" << endl;
      } else {
	shiftWindows.clear();
	recvChildIds.clear();
	vector<MInt> recvChildDomainIds;
	tagAzimuthalUnmappedHaloCells(shiftWindows, recvChildIds, recvChildDomainIds, level);

	// Newly mapped windows
        // If a child of an unmapped azimuthal halo cell has an adequate internal cell
        // This child becomes an regular azimuthal halo cell. Thus, also the internal window cell
        // is added to the azimuthal window cell vector
	for(MInt i = 0; i < noDomains(); i++) {
	  for(MUint j = 0; j < shiftWindows[i].size(); j++ ) {
	    const MLong globalChildId = shiftWindows[i][j];
	    if(globalChildId < 0) continue;

	    ASSERT(findNeighborDomainId(globalChildId) == domainId(),"Wrong domain! " + to_string(globalChildId) + " " + to_string(findNeighborDomainId(globalChildId)) + " " + to_string(domainId()));

	    MInt idx = setAzimuthalNeighborDomainIndex(i, m_azimuthalHaloCells, m_azimuthalWindowCells);

	    const MInt childId = globalIdToLocalId(globalChildId, true);
	    m_azimuthalWindowCells[idx].push_back(childId);
	  }
	}

	// Newly mapped halos
        // If a child of an unmapped azimuthal halo cell has an adequate internal cell
        // This child becomes an regular azimuthal halo cell. Thus, it is added to the azimuthal halo cell vector
        // Otherwise the child becomes an unmapped azimuthal halo cell
	MInt noUnmappedCells = noAzimuthalUnmappedHaloCells();
	for(MInt i = 0; i < noUnmappedCells; i++ ) {
	  MInt cellId = azimuthalUnmappedHaloCell(i);
	  if(a_level(cellId) == level) {
	    if(accumulate(&recvChildIds[m_maxNoChilds * i], &recvChildIds[m_maxNoChilds * i] + m_maxNoChilds,
                          static_cast<MLong>(-2), [](const MLong& a, const MLong& b) { return std::max(a, b); })
               > -2) {

	      refineCell(cellId, nullptr, a_hasProperty(cellId, Cell::IsPartLvlAncestor));
	      refineIds.push_back(cellId);

	      for(MInt child = 0; child < m_maxNoChilds; child++) {
		const MInt childId = a_childId(cellId, child);

		ASSERT(childId > -1 && childId < m_tree.size(), to_string(childId) + " " + to_string(m_tree.size()));
		a_hasProperty(childId, Cell::IsHalo) = a_hasProperty(cellId, Cell::IsHalo);
		a_hasProperty(childId, Cell::IsPeriodic) = a_hasProperty(cellId, Cell::IsPeriodic);
		a_globalId(childId) = recvChildIds[m_maxNoChilds * i + child];
		for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
		  treeb().solver(childId, solver) = treeb().solver(cellId, solver);
		}

		MInt dom = recvChildDomainIds[m_maxNoChilds * i + child];

              if(a_globalId(childId) == -1) {
                m_azimuthalUnmappedHaloCells.push_back(childId);
                m_azimuthalUnmappedHaloDomains.push_back(dom);
                continue;
              }
              
              ASSERT(findNeighborDomainId(a_globalId(childId)) == dom, "Wrong domain! " + to_string(a_globalId(childId)) + " " + to_string(findNeighborDomainId(a_globalId(childId))) + " " + to_string(dom) + " " + to_string(domainId()) + " " + to_string(cellId));

              MInt idx = setAzimuthalNeighborDomainIndex(dom, m_azimuthalHaloCells, m_azimuthalWindowCells);

              m_azimuthalHaloCells[idx].push_back(childId);
            }
          }
        }
        }
      }

      // Update grid bndry cells
      MBool adaptation = false;
      adaptation = Context::getBasicProperty<MBool>("adaptation", AT_, &adaptation);
      if(!adaptation) {
        std::array<MFloat, 2 * nDim> bbox;
	for(MInt i = 0; i < nDim; i++) {
	  bbox[i] = numeric_limits<MFloat>::max();
	  bbox[nDim + i] = numeric_limits<MFloat>::lowest();
	}
	for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
	  for(MInt i = 0; i < nDim; i++) {
	    bbox[i] = mMin(bbox[i], a_coordinate(cellId, i) - F1B2 * cellLengthAtLevel(a_level(cellId)));
	    bbox[nDim + i] = mMax(bbox[nDim + i], a_coordinate(cellId, i) + F1B2 * cellLengthAtLevel(a_level(cellId)));
	  }
	}
	MPI_Allreduce(MPI_IN_PLACE, &bbox[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[0]");
	MPI_Allreduce(MPI_IN_PLACE, &bbox[nDim], nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "bbox[nDim]");
	MBool isBndry = false;
	for(MInt i = 0; i < m_noInternalCells; i++) {
	  MInt cellId =i;  
	  if(a_level(cellId) != level+1) { continue;
}
	  isBndry = false;
	  for(MInt dir = 0; dir < m_noDirs; dir++) {
	    if(a_hasNeighbor(cellId, dir) > 0) { continue;
}
	    if(a_hasParent(cellId) && a_hasNeighbor(a_parentId(cellId), dir) > 0) {
	      MInt nghbrParentId = a_neighborId(a_parentId(cellId),dir);
	      if(a_noChildren(nghbrParentId) == 0) {
		continue;
	      }
	    }
	    if(!m_periodicCartesianDir[dir/2]) {
              std::array<MFloat, nDim> coords;
	      for(MInt d = 0; d < nDim; d++) {
		coords[d] = a_coordinate(cellId,d);
	      }
	      coords[dir / 2] += ((dir % 2) == 0 ? -F1 : F1) * cellLengthAtLevel(m_minLevel);
	      if(coords[dir / 2] < bbox[dir / 2] || coords[dir / 2] > bbox[nDim + dir / 2]) {
		continue;
}
	    }
	    isBndry = true;
	  }
	  if(isBndry) {
	    if(!a_isHalo(cellId)) { m_gridBndryCells.push_back(cellId);
}
	    for(MInt dir = 0; dir < m_noDirs; dir++) {
	      if(a_hasNeighbor(cellId, dir) > 0) {
		MInt nghbrId = a_neighborId(cellId, dir);
		if(!a_isHalo(nghbrId)) {
		  m_gridBndryCells.push_back(nghbrId);
		}
	      }
	    }
	  }
	}
      }
    }

    // sort halo and window cells by globalId to get matching connectivity
    // this may not even be required, check whether ordering is implicity matching for
    // a complex case with multiple partition level shifts
    if(m_maxPartitionLevelShift > 0) {
      for(MInt i = 0; i < noNeighborDomains(); i++) {
        sort(m_haloCells[i].begin(), m_haloCells[i].end(),
             [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
        sort(m_windowCells[i].begin(), m_windowCells[i].end(),
             [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
      }
    }


    if(treeb().noSolvers() > 1) { // Exchange solver info
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
                                m_tree.size());
      if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
        maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(), &m_tree.solverBits(0),
                                  m_tree.size());

        // To ensure that azimuthal halo cells are fully refined (have all children)
        // or are not refined at all, solver Bits need to be checked!
	correctAzimuthalSolverBits();
      }
    }

    // trigger refineCell() for halo cells on the solvers
    if(!refineCellSolver.empty()) {
      // initially empty during solver startup, no solver refinement needed then.
      for(auto& cellId : refineIds) {
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) {
            refineCellSolver[solver](cellId); // call refineCell() function in each solver
          }
        }
      }
      refineIds.clear();
    }

    // Set window/halo flags
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        ASSERT(!a_hasProperty(m_haloCells[i][j], Cell::IsWindow), "halo cell is marked as window");
        a_hasProperty(m_haloCells[i][j], Cell::IsHalo) = true;
        // a_hasProperty(m_haloCells[i][j], Cell::IsWindow) = false;
      }
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        ASSERT(!a_isHalo(m_windowCells[i][j]), "window cell is marked as halo");
        // a_hasProperty( m_windowCells[i][j], Cell::IsHalo ) = false;
        a_hasProperty(m_windowCells[i][j], Cell::IsWindow) = true;
      }
    }
    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
          ASSERT(!a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsWindow), "azimuthal halo cell is marked as window");
          a_hasProperty(m_azimuthalHaloCells[i][j], Cell::IsHalo) = true;
        }
        for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
          ASSERT(!a_isHalo(m_azimuthalWindowCells[i][j]), "azimuthal window cell is marked as halo");
          a_hasProperty(m_azimuthalWindowCells[i][j], Cell::IsWindow) = true;
        }
      }
      for(MInt j = 0; j < noAzimuthalUnmappedHaloCells(); j++ ) {
        ASSERT(!a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsWindow), "azimuthal halo cell is marked as window");
      a_hasProperty(m_azimuthalUnmappedHaloCells[j], Cell::IsHalo) = true;
      }
    }
  }

  if(domainId() == 0 && onlyLevel < 0) cerr << endl;
  
}

//-------------------------------------------------------------------------------------------

/* \param[in] data data field, whre each field of the bitset corresponds to a solver
 * \param[in] defaultVal
 *
 * The concept of the following is to exchange 'data' from window cells to halo cells. 'data' is a bitset
 * with elements for each solver. If for the i-th solver a cell is not a window cell, the value 'defaultVal'
 * is exchanged instead of data[cellId][i]. One use case are the solverBits.
 */
template <MInt nDim>
template <std::size_t N>
void CartesianGrid<nDim>::exchangeSolverBitset(std::bitset<N>* const data, const MBool defaultVal) {
  TRACE();

  //maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
  //                              m_tree.size());
  //  return;

  // There might be the case that a halo cell resides inside the domain of a solver, but since that solver
  // does not require that many layers we disable that cell for that solver. Sounds crazy, but trust me.
  static_assert(N <= 64, "conversion to ulong not appropriate, change to ullong!");
  ScratchSpace<MPI_Request> recvRequests(max(1, noNeighborDomains()), AT_, "recvRequests");
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  MInt receiveCount = 0;
  for (const auto& vecHalo : m_haloCells) receiveCount+=vecHalo.size();
  ScratchSpace<MUlong> haloBuffer(max(1, receiveCount), AT_, "haloBuffer");
  for (MInt i = 0, offset = 0; i < noNeighborDomains(); i++) {
    const MInt noHaloCells = m_haloCells[i].size();
    if (noHaloCells<1) continue;

    MPI_Irecv(&haloBuffer[offset], noHaloCells, type_traits<MUlong>::mpiType(), m_nghbrDomains[i],
        m_nghbrDomains[i], mpiComm(), &recvRequests[i], AT_, "haloBuffer[offset]");

    offset += noHaloCells;
  }


  ScratchSpace<MPI_Request> sendRequests(max(1, noNeighborDomains()), AT_, "sendRequests");
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  MInt sendCount = 0;
  for (const auto& vecWindow : m_windowCells) sendCount+=vecWindow.size();
  ScratchSpace<MUlong> tmp_data(sendCount, AT_, "tmp_data");
  for (MInt i = 0, idx = 0; i < noNeighborDomains(); i++) {
    const MInt offset = idx;
    const MInt noWindowCells = m_windowCells[i].size();
    if (noWindowCells<1) continue;
    for (const auto cellId : m_windowCells[i]) {
      const auto backup = data[cellId].to_ulong();
      ASSERT(m_windowLayer_[i].find(cellId)!=m_windowLayer_[i].end(), "You don't know what you are doing!");
      for (MInt solverId = 0; solverId < treeb().noSolvers(); ++solverId) {
        if (!isSolverWindowCell(i, cellId, solverId)) {
          data[cellId][solverId] = defaultVal;
        }
      }
      tmp_data[idx++] = data[cellId].to_ulong();
      data[cellId] = std::bitset<N>(backup);
    }
    ASSERT(idx-offset==noWindowCells, "");

    MPI_Isend(&tmp_data[offset], noWindowCells, type_traits<MUlong>::mpiType(), m_nghbrDomains[i], domainId(),
        mpiComm(), &sendRequests[i], AT_, "tmp_data");
  }

  // Finish MPI communication
  MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  for (MInt i = 0, idx = 0; i < noNeighborDomains(); i++) {
    for (const auto cellId : m_haloCells[i])
      data[cellId] = std::bitset<N>(haloBuffer[idx++]);
  }
}


/** \brief Create a new halo cell (candidate) as neighbor of cellId in direction dir and find matching neighbor domain
   from hilbertIndex \author Lennart Schneiders \date October 2017
  */
template <MInt nDim>
MInt CartesianGrid<nDim>::createAdjacentHaloCell(const MInt cellId, const MInt dir, const MLong* hilbertOffsets,
                                                     unordered_multimap<MLong, MInt>& hilbertToLocal,
                                                     const MFloat* bbox, MInt* const noHalos,
                                                     vector<vector<MInt>>& halos, vector<vector<MLong>>& hilbertIds) {
  static constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};
  const MFloat cellLength = cellLengthAtCell(cellId);

  MFloat coords[3];
  MFloat dcoords[3];
  MFloat dummyCoords[3];
  for(MInt i = 0; i < nDim; i++) {
    coords[i] = a_coordinate(cellId, i);
  }
  coords[dir / 2] += ((dir % 2) == 0 ? -F1 : F1) * cellLength;

  for(MInt i = 0; i < nDim; i++) {
    if(!m_periodicCartesianDir[dir / 2] && (coords[dir / 2] < bbox[dir / 2] || coords[dir / 2] > bbox[nDim + dir / 2]))
      return -1;
  }

  MBool insd = true;
  for(MInt i = 0; i < nDim; i++) {
    if(coords[i] < bbox[i] || coords[i] > bbox[nDim + i]) insd = false;
  }
  MBool isPeriodic = !insd;

  if(m_azimuthalPer && isPeriodic) return -1;

  MFloat dx[3];
  for(MInt i = 0; i < nDim; i++) {
    dx[i] = bbox[nDim + i] - bbox[i];
  }
  for(MInt i = 0; i < nDim; i++) {
    dummyCoords[i] = coords[i];
  }
  if(isPeriodic) {
    for(MInt i = 0; i < nDim; i++) {
      if(coords[i] < bbox[i])
        dummyCoords[i] += dx[i];
      else if(coords[i] > bbox[nDim + i])
        dummyCoords[i] -= dx[i];
    }
  }

  for(MInt i = 0; i < nDim; i++) {
    if(dummyCoords[i] < bbox[i] || dummyCoords[i] > bbox[nDim + i]) {
      mTerm(1, AT_, "Halo cell coords outside domain");
    }
  }

  const MLong hilbertIndex = hilbertIndexGeneric(&dummyCoords[0]);
  MInt ndom = -1;
  while(hilbertIndex >= hilbertOffsets[ndom + 1]) {
    ndom++;
    if(ndom == noDomains() - 1) break;
  }
  ASSERT(ndom > -1, "");
  ASSERT(hilbertIndex >= hilbertOffsets[ndom] && hilbertIndex < hilbertOffsets[ndom + 1],
         to_string(hilbertIndex) + " " + to_string(hilbertOffsets[ndom]) + " " + to_string(hilbertOffsets[ndom + 1]));

  const MInt haloCellId = m_tree.size();
  m_tree.append();

  for(MInt i = 0; i < nDim; i++) {
    a_level(haloCellId) = a_level(cellId);
  }

  for(MInt i = 0; i < nDim; i++) {
    a_coordinate(haloCellId, i) = coords[i];
  }

  hilbertToLocal.insert(make_pair(hilbertIndex, haloCellId));

  for(MInt i = 0; i < m_noDirs; i++) {
    a_neighborId(haloCellId, i) = -1;
  }
  for(MInt i = 0; i < m_maxNoChilds; i++) {
    a_childId(haloCellId, i) = -1;
  }
  a_parentId(haloCellId) = -1;

  a_neighborId(cellId, dir) = haloCellId;
  a_neighborId(haloCellId, revDir[dir]) = cellId;
  for(MInt otherDir = 0; otherDir < m_noDirs; otherDir++) {
    if(otherDir / 2 == dir / 2) continue;
    if(a_hasNeighbor(cellId, otherDir) == 0) continue;
    MInt nghbrId = a_neighborId(cellId, otherDir);
    if(a_hasNeighbor(nghbrId, dir) > 0) {
      nghbrId = a_neighborId(nghbrId, dir);
      a_neighborId(nghbrId, revDir[otherDir]) = haloCellId;
      a_neighborId(haloCellId, otherDir) = nghbrId;
    }
  }
  for(MInt ndir = 0; ndir < m_noDirs; ndir++) {
    if(a_hasNeighbor(haloCellId, ndir) > 0) continue;
    for(MInt i = 0; i < nDim; i++) {
      dummyCoords[i] = a_coordinate(haloCellId, i);
    }
    dummyCoords[ndir / 2] += ((ndir % 2) == 0 ? -F1 : F1) * cellLength;
    for(MInt i = 0; i < nDim; i++) {
      dcoords[i] = dummyCoords[i];
    }
    for(MInt i = 0; i < nDim; i++) {
      if(dummyCoords[i] < bbox[i])
        dummyCoords[i] += dx[i];
      else if(dummyCoords[i] > bbox[nDim + i])
        dummyCoords[i] -= dx[i];
    }
    const MLong nghbrHilbertId = hilbertIndexGeneric(&dummyCoords[0]);
    MInt nghbrId = -1;
    auto range = hilbertToLocal.equal_range(nghbrHilbertId);
    for(auto it = range.first; it != range.second; ++it) {
      MFloat dist = F0;
      for(MInt i = 0; i < nDim; i++) {
        dist += POW2(a_coordinate(it->second, i) - dcoords[i]);
      }
      dist = sqrt(dist);
      if(dist < 0.001 * cellLength) {
        if(nghbrId > -1) cerr << "duplicate " << hilbertToLocal.count(nghbrHilbertId) << endl;
        nghbrId = it->second;
      }
    }
    if(nghbrId < 0) continue;

#ifndef NDEBUG
    for(MInt i = 0; i < nDim; i++) {
      ASSERT(fabs(a_coordinate(nghbrId, i) - a_coordinate(haloCellId, i)) < cellLength * 1.0001, "");
    }
#endif
    a_neighborId(nghbrId, revDir[ndir]) = haloCellId;
    a_neighborId(haloCellId, ndir) = nghbrId;
  }

  a_resetProperties(haloCellId);
  a_hasProperty(haloCellId, Cell::IsHalo) = true;
  a_hasProperty(haloCellId, Cell::IsPeriodic) = isPeriodic;

  MInt idx = setNeighborDomainIndex(ndom, halos, hilbertIds);
  ASSERT(idx > -1, "");
  halos[idx].push_back(haloCellId);
  hilbertIds[idx].push_back(hilbertIndex);
  noHalos[ndom]++;

  return haloCellId;
}


// -------------------------------------------------------------------------------------------


/** \brief Fill m_nghbrDomains, m_haloCells and m_windowCells with data from temporary buffers
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
void CartesianGrid<nDim>::updateHaloCellCollectors() {
  TRACE();
  const auto noNghbrDomains0 = (signed)m_nghbrDomains.size();
  const auto noAzimuthalNghbrDomains0 = (signed)m_azimuthalNghbrDomains.size();
  if(noNghbrDomains0 == 0 && noAzimuthalNghbrDomains0 == 0) {
    // serial computations can have neighbor domains, for example, when using
    // periodic boundaries, and that requires halo cells.
    //
    // That is, if there aren't any neighbor domains, halo cells are not required.
    return;
  }

  ScratchSpace<MInt> domMap(noNghbrDomains0, AT_, "domMap");

  vector<pair<MInt, MInt>> tmpDoms;
  tmpDoms.reserve(noNghbrDomains0);

  // sort neighbor domain ids to avoid send/recv deadlocks
  for(MInt i = 0; i < noNghbrDomains0; i++) {
    if(m_haloCells[i].size() > 0 || m_windowCells[i].size() > 0) {
      tmpDoms.push_back(make_pair(m_nghbrDomains[i], i));
    }
  }
  sort(tmpDoms.begin(), tmpDoms.end()); // sort by ascending neighbor domain id

  std::vector<MInt>().swap(m_nghbrDomains);
  ASSERT(m_nghbrDomainIndex.size() >= static_cast<size_t>(noDomains()),
         "m_nghbrDomainIndex size is " + std::to_string(m_nghbrDomainIndex.size())
             + " which is smaller than noDomains() = " + std::to_string(noDomains()));
  std::fill_n(m_nghbrDomainIndex.begin(), noDomains(), -1);
  for(MInt i = 0; i < (signed)tmpDoms.size(); i++) {
    domMap[i] = tmpDoms[i].second;
    m_nghbrDomainIndex[tmpDoms[i].first] = m_nghbrDomains.size();
    m_nghbrDomains.push_back(tmpDoms[i].first);
    if(m_nghbrDomains[i] == domainId() && m_noPeriodicCartesianDirs == 0) {
      TERMM(1, "Not supposed to happen: connection to self without periodicity.");
    }
    if(m_nghbrDomains[i] == domainId()) m_log << domainId() << ": periodic connection to self." << endl;
  }


  // sort window and halo cells by globalId
  const MBool sortHaloWindowCells = false;
  if(sortHaloWindowCells) {
    for(MInt i = 0; i < noNghbrDomains0; i++) {
      sort(m_haloCells[i].begin(), m_haloCells[i].end(),
           [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
      sort(m_windowCells[i].begin(), m_windowCells[i].end(),
           [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });
    }
  }


  // write window/halo cells
  vector<std::vector<MInt>> haloCellBak(m_haloCells);
  vector<std::vector<MInt>> windowCellBak(m_windowCells);
  m_haloCells.resize(noNeighborDomains());
  m_windowCells.resize(noNeighborDomains());
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    m_haloCells[i].resize(haloCellBak[domMap[i]].size());
    m_windowCells[i].resize(windowCellBak[domMap[i]].size());
    ASSERT(m_haloCells[i].size() >= haloCellBak[domMap[i]].size(), "");
    ASSERT(m_windowCells[i].size() >= windowCellBak[domMap[i]].size(), "");
    copy(haloCellBak[domMap[i]].begin(), haloCellBak[domMap[i]].end(), m_haloCells[i].begin());
    copy(windowCellBak[domMap[i]].begin(), windowCellBak[domMap[i]].end(), m_windowCells[i].begin());
  }

  // TODO_SS labels:GRID,toenhance use move operation
  if (m_windowLayer_.size()>0 && m_haloMode>0) { //WH_old
    ASSERT(noNeighborDomains()<=(signed)m_windowLayer_.size(), "");
    std::vector<std::unordered_map<MInt, M32X4bit<true>>> windowLayerBak(m_windowLayer_);
    m_windowLayer_.resize(noNeighborDomains());
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      m_windowLayer_[i] = windowLayerBak[domMap[i]];
    }
  }


  if(m_azimuthalPer) {
    ScratchSpace<MInt> domMapAzi(noAzimuthalNghbrDomains0, AT_, "domMapAzi");
    // sort neighbor domain ids to avoid send/recv deadlocks
    tmpDoms.clear();
    tmpDoms.reserve(noAzimuthalNghbrDomains0);
    for(MInt i = 0; i < noAzimuthalNghbrDomains0; i++) {
      if(m_azimuthalHaloCells[i].size() > 0 || m_azimuthalWindowCells[i].size() > 0) {
        tmpDoms.push_back(make_pair(m_azimuthalNghbrDomains[i], i));
      }
    }
    sort(tmpDoms.begin(), tmpDoms.end()); // sort by ascending neighbor domain id

    m_azimuthalNghbrDomains.clear();
    ASSERT(m_azimuthalNghbrDomainIndex.size() >= static_cast<size_t>(noDomains()),
           "m_azimuthalNghbrDomainIndex size is " + std::to_string(m_azimuthalNghbrDomainIndex.size())
           + " which is smaller than noDomains() = " + std::to_string(noDomains()));
    std::fill_n(m_azimuthalNghbrDomainIndex.begin(), noDomains(), -1);
    for(MInt i = 0; i < (signed)tmpDoms.size(); i++) {
      domMapAzi[i] = tmpDoms[i].second;
      m_azimuthalNghbrDomainIndex[tmpDoms[i].first] = m_azimuthalNghbrDomains.size();
      m_azimuthalNghbrDomains.push_back(tmpDoms[i].first);
      if(m_azimuthalNghbrDomains[i] == domainId()) m_log << domainId() << ": periodic connection to self." << endl;
    }

    // sort azimuthal window and halo cells by globalId
    if(sortHaloWindowCells) {
      vector<MInt> posMapAzi;
      for(MInt i = 0; i < noAzimuthalNghbrDomains0; i++) {
        sort(m_azimuthalHaloCells[i].begin(), m_azimuthalHaloCells[i].end(),
             [this](const MInt& a, const MInt& b) { return a_globalId(a) < a_globalId(b); });

	vector<pair<MInt, MInt>> azimuthalWindowCells;
	for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
	  azimuthalWindowCells.push_back(make_pair(m_azimuthalWindowCells[i][j], j));
	}
	m_azimuthalWindowCells[i].clear();
        sort(azimuthalWindowCells.begin(), azimuthalWindowCells.end(),
             [&](const auto& a, const auto& b) { return a_globalId(a.first) < a_globalId(b.first); });
	m_azimuthalWindowCells[i].resize(azimuthalWindowCells.size());
	posMapAzi.resize(azimuthalWindowCells.size());
	for(MInt j = 0; j < (signed)azimuthalWindowCells.size(); j++) {
	  m_azimuthalWindowCells[i][j] = azimuthalWindowCells[j].first;
	  posMapAzi[azimuthalWindowCells[j].second] = j;
	}
	vector<MInt> higherLevelConnectivityBak(m_azimuthalHigherLevelConnectivity[i]);
	m_azimuthalHigherLevelConnectivity[i].resize(higherLevelConnectivityBak.size());
	for(MInt c = 0; c < (signed)higherLevelConnectivityBak.size(); c++) {
	  m_azimuthalHigherLevelConnectivity[i].push_back(posMapAzi[higherLevelConnectivityBak[c]]);
	}
	azimuthalWindowCells.clear();
	higherLevelConnectivityBak.clear();
      }
    }

    // write window/halo cells
    haloCellBak.clear();
    windowCellBak.clear();
    haloCellBak = m_azimuthalHaloCells;
    windowCellBak = m_azimuthalWindowCells;
    vector<vector<MInt>> higherLevelConnectivityBak2(m_azimuthalHigherLevelConnectivity);
    m_azimuthalHaloCells.resize(noAzimuthalNeighborDomains());
    m_azimuthalWindowCells.resize(noAzimuthalNeighborDomains());
    m_azimuthalHigherLevelConnectivity.resize(noAzimuthalNeighborDomains());
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      m_azimuthalHaloCells[i].resize(haloCellBak[domMapAzi[i]].size());
      m_azimuthalWindowCells[i].resize(windowCellBak[domMapAzi[i]].size());
      m_azimuthalHigherLevelConnectivity[i].resize(higherLevelConnectivityBak2[domMapAzi[i]].size());
      ASSERT(m_azimuthalHaloCells[i].size() >= haloCellBak[domMapAzi[i]].size(), "");
      ASSERT(m_azimuthalWindowCells[i].size() >= windowCellBak[domMapAzi[i]].size(), "");
      copy(haloCellBak[domMapAzi[i]].begin(), haloCellBak[domMapAzi[i]].end(), m_azimuthalHaloCells[i].begin());
      copy(windowCellBak[domMapAzi[i]].begin(), windowCellBak[domMapAzi[i]].end(), m_azimuthalWindowCells[i].begin());
      copy(higherLevelConnectivityBak2[domMapAzi[i]].begin(), higherLevelConnectivityBak2[domMapAzi[i]].end(), m_azimuthalHigherLevelConnectivity[i].begin());
    }
  }
}


// -------------------------------------------------------------------------------------------


/** \brief Flag m_noHaloLayers layers of halo cells adjacent to internal cells on this domain
    \author Lennart Schneiders
    \date October 2017
  */
template <MInt nDim>
void CartesianGrid<nDim>::tagActiveWindows(vector<MLong>& refineChildIds, MInt level) {
  constexpr MInt maxNghbrDoms = 1024;
  MIntScratchSpace nghbrList(27 * m_maxNoChilds, AT_, "nghbrList");
  vector<std::bitset<maxNghbrDoms>> nghrDomFlag(m_tree.size());
  vector<std::bitset<maxNghbrDoms>> nghrDomFlagBak;
  vector<std::bitset<maxNghbrDoms>> nghrDomFlagBak2;
  set<MInt> exchangeCellList;
  if(noNeighborDomains() > maxNghbrDoms) {
    mTerm(1, AT_, "Too many neighbor domains(" + std::to_string(noNeighborDomains()) + "). Increase bitset size.");
  }
  // temporary solution for partitionLevelShifts, might produce a lot of halo cells, more clever way to flag active
  // windows should be found
  // TODO labels:GRID @ansgar_pls_adapt fix this!
  constexpr MBool testingFix_partitionLevelShift = true; // false;

  // TODO labels:GRID @timw_multiSolverHalos remove old version!
  constexpr MBool multiSolverHaloLayer = true; // should be true

  MInt windowCnt = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++)
    windowCnt += m_windowCells[i].size();
  refineChildIds.resize(windowCnt * m_maxNoChilds);
  fill(refineChildIds.begin(), refineChildIds.end(), -1);


  // identify window cells on the next higher level, or coarser leaf window cells if not existing
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      const MInt cellId = m_windowCells[i][j];

      //      if ( a_level( cellId ) != level ) continue; //however, such cells should not exist, yet
      if(a_noChildren(cellId) > 0) {
        if(a_level(cellId) == level) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt childId = a_childId(cellId, child);
            if(childId < 0) continue;

            // fix for unconsistent window tags with double window cells (domain interface and periodic)
            // const MBool doubleWindow = exchangeCellList.find(cellId) != exchangeCellList.end();
            // if(!doubleWindow){
            exchangeCellList.insert(childId);
            //}
          }
        }
      } else {
        exchangeCellList.insert(cellId);
      }
    }
  }

  // add halo cells
  // NOTE: the loop below or similar assumptions that
  //       internal/halo cells are sorted are wrong during adaptation
  // for (MInt cellId = m_noInternalCells; cellId < m_tree.size(); cellId++) {
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_isHalo(cellId)) {
      exchangeCellList.insert(cellId);
    }
  }

  // exchangeCellList conaints all window and halo-Cells,
  // and additionally all children of all windowCells

  if(m_maxPartitionLevelShift > 0) {
    // exchange some window cell data
    ScratchSpace<MInt> bprops(m_tree.size(), 2, AT_, "bprops");
    bprops.fill(-1);
    // Note: changed loop to iterate over whole tree, internal and halo cells not sorted!
    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      if(a_isHalo(cellId) || a_isToDelete(cellId)) {
        continue;
      }
      bprops(cellId, 0) = (MInt)a_hasProperty(cellId, Cell::IsPartLvlAncestor);
      bprops(cellId, 1) = a_noOffsprings(cellId);

      ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
             std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                 + std::to_string(bprops(cellId, 1)));
    }
    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), bprops.getPointer(), 2);

    // Note: loop only over the current haloCells since only for these data is exchanged, this does
    // not work properly if in case of a partition level shift some halo cell is not part of
    // m_haloCells yet such that its property IsPartLvlAncestor is reset!
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        const MInt cellId = m_haloCells[i][j];
        ASSERT(!a_isToDelete(cellId), "");
        ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
               std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                   + std::to_string(bprops(cellId, 1)));
        a_hasProperty(cellId, Cell::IsPartLvlAncestor) = (MBool)bprops(cellId, 0);
        a_noOffsprings(cellId) = bprops(cellId, 1);
      }
    }
  }

  for(MInt i = m_tree.size(); i--;) {
    for(MInt j = 0; j < maxNghbrDoms; j++) {
      ASSERT(!nghrDomFlag[i][j], "");
    }
  }

  // mark all halo cells as nghbrdomain
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      const MInt cellId = m_haloCells[i][j];

      nghrDomFlag[cellId][i] = true;

      if(m_maxPartitionLevelShift > 0) {
        if(a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
          MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          MInt maxNghbrDomId =
              findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)));
          for(MInt k = 0; k < noNeighborDomains(); k++) {
            if(m_nghbrDomains[k] >= minNghbrDomId && m_nghbrDomains[k] <= maxNghbrDomId) {
              nghrDomFlag[cellId][k] = true;
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                MInt childId = a_childId(cellId, child);
                if(childId > -1) {
                  nghrDomFlag[childId][k] = true;
                }
              }
            }
          }
        }
      }
    }
  }

  map<MInt, MInt> solverCellsAdded;
  solverCellsAdded.clear();

  // extend halo cells by the defined number of halo layers
  // NOTE: m_noHaloLayers corresponds to the number of halo-Layers of the grid
  //      for multiSolver applications some solvers may have a lower number of halo-Layers!
  for(MInt layer = 0; layer < m_noHaloLayers; layer++) {
    nghrDomFlagBak.assign(nghrDomFlag.begin(), nghrDomFlag.end());

    for(MInt cellId : exchangeCellList) {
      if(a_level(cellId) > (level + 1)) continue;
      if(a_level(cellId) < (level + 1) && a_noChildren(cellId) > 0) continue;
      const MInt counter = getAdjacentGridCells(cellId, nghbrList, level);
      // check diffs when called instead!!!!
      for(MInt n = 0; n < counter; n++) {
        MInt nghbrId = nghbrList[n];
        if(nghbrId < 0) continue;

        // if the neighbor doesn't belong to all solvers
        // an additional layer needs to be added to the layer is below the noHaloLayers
        // for the solver, which does not have the neighbor!
        // to ensure that this each solver has at least the number of halo layers required!
        if(g_multiSolverGrid && multiSolverHaloLayer && !m_paraViewPlugin) {
          for(MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
            if(!treeb().solver(nghbrId, solverId) && layer < m_noSolverHaloLayers[solverId] /*&& layer > 0*/) {
              solverCellsAdded.insert(make_pair(nghbrId, solverId));
            }
          }
        }

        for(MInt i = 0; i < noNeighborDomains(); i++) {
          if(nghrDomFlagBak[nghbrId][i]) {
            nghrDomFlag[cellId][i] = true;
          }
        }

        // fix for unconsistent window tags with double window cells (domain interface and periodic)
        // MBool tagWindow = false;
        // for ( MInt i = 0; i < noNeighborDomains(); i++ ) {

        //   if(nghrDomFlagBak[ nghbrId ][ i ]){
        //     //nghrDomFlag[ cellId ][ i ] = true;
        //     tagWindow = true;
        //     break;
        //   }
        // }
        // if(tagWindow){
        //   for ( MInt i = 0; i < noNeighborDomains(); i++ ) {
        //     nghrDomFlag[ cellId ][ i ] = true;
        //   }
        // }


        if(m_maxPartitionLevelShift > 0) {
          MInt rootId = a_parentId(cellId) > -1 ? a_parentId(cellId) : cellId;
          if(a_hasProperty(rootId, Cell::IsPartLvlAncestor)) {
            MInt minNghbrDomId = findNeighborDomainId(a_globalId(rootId));
            MInt maxNghbrDomId =
                findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(rootId) + (MLong)a_noOffsprings(rootId)));
            for(MInt k = 0; k < noNeighborDomains(); k++) {
              if(m_nghbrDomains[k] >= minNghbrDomId && m_nghbrDomains[k] <= maxNghbrDomId) {
                nghrDomFlag[cellId][k] = true;
              }
            }
          }
        }
      }
    }

    // now extend one additional layer around cells which might be missing the layer
    if(!solverCellsAdded.empty()) {
      nghrDomFlagBak2.assign(nghrDomFlag.begin(), nghrDomFlag.end());

      for(map<MInt, MInt>::iterator it = solverCellsAdded.begin(); it != solverCellsAdded.end(); it++) {
        const MInt cellId = it->first;
        const MInt solverId = it->second;

        // get additional cells
        // NOTE: diagonal neighbors are necessary as well, however checkNoHaloLayers
        //      only checks for direct neighbors!
        const MInt counter = getAdjacentGridCells(cellId, nghbrList, level, true);

        MBool anyNeighbor = false;
        const MInt size = noNeighborDomains();
        MBoolScratchSpace addedCell(size, AT_, "addedCell");
        fill(addedCell.begin(), addedCell.end(), false);
        for(MInt n = 0; n < counter; n++) {
          const MInt nghbrId = nghbrList[n];
          if(nghbrId < 0) continue;
          MInt parentId = a_parentId(cellId);
          if(parentId > -1) {
            if(!treeb().solver(parentId, solverId)) continue;
          }
          // labels:GRID testcase hack
          if(m_zonal) {
            MBool checkParentChilds = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId > -1) {
                if(!treeb().solver(childId, solverId)) {
                  checkParentChilds = true;
                  break;
                }
              }
            }
            if(checkParentChilds) continue;
          }
          if(treeb().solver(nghbrId, solverId)) anyNeighbor = true;
          // if( !treeb().solver( nghbrId, solverId ) &&
          //    a_level(nghbrId) == a_level(cellId)) continue;
          for(MInt i = 0; i < noNeighborDomains(); i++) {
            if(nghrDomFlagBak2[nghbrId][i]) {
              if(!nghrDomFlag[cellId][i]) {
                addedCell[i] = true;
              }
              nghrDomFlag[cellId][i] = true;
            }
          }
        }
        // remove the cell from the list again if:
        // no neighbors for the cell and intended solver could be found!
        // Meaning that the solver has a lower max-level or different boundingBox!
        if(!anyNeighbor) {
          // for ( MInt n = 0; n < counter; n++ ) {
          // const MInt nghbrId = nghbrList[ n ];
          for(MInt i = 0; i < noNeighborDomains(); i++) {
            if(addedCell[i]) {
              nghrDomFlag[cellId][i] = false;
            }
          }
          //}
        }
      }
      solverCellsAdded.clear();
      nghrDomFlagBak2.clear();
    }
  }

  // not openmp compatible
  MInt cnt = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      MInt cellId = m_windowCells[i][j];
      if(a_level(cellId) == level) {
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          MInt childId = a_childId(cellId, child);
          if(childId > -1) {
            if(nghrDomFlag[childId][i]) {
              refineChildIds[m_maxNoChilds * cnt + child] = 1;
            } else if(testingFix_partitionLevelShift && a_level(cellId) < m_minLevel + m_maxPartitionLevelShift) {
              refineChildIds[m_maxNoChilds * cnt + child] = 1;
            }
          }
        }
      }
      cnt++;
    }
  }
}


/** \brief
 * NOTE: Actually the most clever way of tagging is to go from highest level to lower levels, but that
 *       means quite some work to do!
    \date June 2021
  */
// clang-format off
template <MInt nDim>
void CartesianGrid<nDim>::tagActiveWindows2_(vector<MLong>& refineChildIds, const MInt level) {

  MIntScratchSpace nghbrList(27 * m_maxNoChilds, AT_, "nghbrList");

  // We rely on the fact, that exchangeCellList is sorted by the tuple's first element, then by its 2nd element etc.
  std::set<std::tuple<MInt/*layer*/,MInt/*cellId*/, MInt/*nghbrIdx*/,MInt/*solverId*/>> exchangeCellList;
  // partLvlShift
  auto newNghbrDomainsIndex = m_nghbrDomainIndex;
  auto newNghbrDomains = m_nghbrDomains;

#ifndef NDEBUG
  std::set<MInt> allWindows;
#endif

  // The algorithm for m_maxPartitionLevelShift>0 requires a_hasProperty(cellId, Cell::IsPartLvlAncestor) &
  // a_noOffsprings(cellId) of the halo cells to be set. This should be the case at this place.

  // Get halo candidates
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      const MInt cellId = m_haloCells[i][j];
      // If this is called during meshAdaptation and grid has partLvlShift, then we might have halo cells for which
      // a_level(cellId)>level
      ASSERT(a_level(cellId)<=level || m_maxPartitionLevelShift>0, "");
      if (a_level(cellId)>level) continue;

      // Note: currently we add all halos as candidates. A more sophisticated choice would be to add only those
      //       who are adajcent to cells, which are not halos of same neighbor (the neighbor can still be a halo)

      //TODO_SS labels:GRID what about the case a_noChildren(cellId)<IPOW2(nDim)
      if (a_level(cellId)==level || (a_noChildren(cellId)==0 && a_level(cellId)<level)) {
        for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
          // For the following to work we need to exchange the solverBits before tagActiveWindowsAtMinLevel is
          // called (see comments in tagActiveWindowsAtMinLevel)
          // if (treeb().solver(cellId, solverId) || treeb().noSolvers()==1)
          if (treeb().solver(cellId, solverId)) {
            exchangeCellList.insert({std::make_tuple(0, cellId, i, solverId)});
          }
        }
      }

      // Before we do the following we need to exchange the property IsPartLvlAncesotr and a_noOffsping (see above)
      if (a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
        if (a_level(cellId)==level || (a_noChildren(cellId)==0 && a_level(cellId)<level)) {

          // TODO_SS labels:GRID what about the case a_noChildren(cellId)<IPOW2(nDim)
          const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          const MInt maxNghbrDomId =
            findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));

          for (MInt nghbrDom = minNghbrDomId; nghbrDom <= maxNghbrDomId; ++nghbrDom) {
            MInt k = newNghbrDomainsIndex[nghbrDom];
            if (nghbrDom==domainId()) continue; //TODO_SS labels:GRID not sure if this is correct
            //if (k<0 && nghbrDom==domainId()) continue; // periodic???
            if (k<0) {
              // We have a new neighbor
              k = newNghbrDomains.size();
              newNghbrDomainsIndex[nghbrDom] = k;
              newNghbrDomains.push_back(nghbrDom);
            }
            //if (i==k) continue;
//            if(m_nghbrDomains[k] >= minNghbrDomId && m_nghbrDomains[k] <= maxNghbrDomId) { }
            for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
              // For the following to work we need to exchange the solverBits before tagActiveWindowsAtMinLevel is
              // called (see comments in tagActiveWindowsAtMinLevel)
              // if (treeb().solver(cellId, solverId) || treeb().noSolvers()==1)
              if (treeb().solver(cellId, solverId)) {
                exchangeCellList.insert({std::make_tuple(0, cellId, k, solverId)});

                for(MInt child = 0; child < m_maxNoChilds; child++) {
                  const MInt childId = a_childId(cellId, child);
                  if(childId > -1) {

                    //TODO_SS labels:GRID exchange noOffspring of halo cells (check if it is already done, especially those halo
                    //         cells which are not yet in m_haloCells, but are created because of partLvlShift)
                    const MInt minNghbrDomId2 = findNeighborDomainId(a_globalId(childId));
                    const MInt maxNghbrDomId2 =
                      findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(childId) + (MLong)a_noOffsprings(childId)-1));
                    MInt lay = (newNghbrDomains[k] >= minNghbrDomId2 && newNghbrDomains[k] <= maxNghbrDomId2) ? 0 : 1;
                    // TODO_SS labels:GRID since halo data are not exchanged yet, safety approach
                    if (a_isHalo(childId)) lay = 0;

                    exchangeCellList.insert({std::make_tuple(lay, childId, k, solverId)});
                    if (!a_isHalo(childId)) {
                      if (treeb().solver(childId, solverId)) {
                        //TODO_SS labels:GRID,totest following assert fails in periodic case when nghbrDom==domainId() --> check
                        //TERMM_IF_COND(lay==0 && !a_hasProperty(childId, Cell::IsPartLvlAncestor), to_string(minNghbrDomId2) + " " + to_string(maxNghbrDomId2));
                        //TODO_SS labels:GRID can we ensure that the newly found domain will also find current domain?
                        m_windowLayer_[setNeighborDomainIndex(nghbrDom, m_haloCells, m_windowCells, m_windowLayer_)][childId].set(solverId, lay);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

//  vector<std::set<MInt>> windCellSet(noNeighborDomains());
//  for(MInt i = 0; i < noNeighborDomains(); i++) {
//    windCellSet[i].insert(m_windowCells[i].begin(), m_windowCells[i].end());
//  }

  // Get window candidates on lower level
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      const MInt cellId = m_windowCells[i][j];
#ifndef NDEBUG
      allWindows.insert(cellId);
      TERMM_IF_COND(m_windowLayer_[i].find(cellId)==m_windowLayer_[i].end(), std::to_string(m_maxPartitionLevelShift));
      TERMM_IF_COND(m_windowLayer_[i].at(cellId).all(), "");
#endif

      if (a_level(cellId)==level) {
        for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
          const MInt windowLayer = m_windowLayer_[i].at(cellId).get(solverId);
          // TODO_SS labels:GRID what about the case that a_noChildren<IPOW2(nDim)
          if (!a_hasChildren(cellId, solverId) && windowLayer<m_noSolverHaloLayers[solverId]) {
            // labels:GRID,DLB With DLB the following check fails: Bug?!
//            if (windowLayer<1 && !a_hasProperty(cellId, Cell::IsPartLvlAncestor))
//              TERMM(1, "");
            exchangeCellList.insert({std::make_tuple(windowLayer, cellId, i, solverId)});
          }
        }
      }

      // 1) window cell (a_level(cellId)==level) with IPOW2(nDim) number of window children
      // 2) window cell (a_level(cellId)==level) with <IPOW2(nDim) number of window children
      if (a_hasProperty(cellId, Cell::IsPartLvlAncestor) && a_level(cellId)==level) {
        ASSERT(a_noChildren(cellId)>0, "");
        ASSERT(!a_isHalo(cellId), "");

        MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
        MInt maxNghbrDomId =
          findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
#ifndef NDEBUG
        // All neighbor domains, who have children above this current window cell, should already be contained in
        // m_nghbrDomains
        for (MInt nD = minNghbrDomId; nD<=maxNghbrDomId; ++nD) {
          if (nD==domainId()) continue;
          MBool found = false;
          for (MInt ii = 0; ii < noNeighborDomains(); ii++) {
            if (m_nghbrDomains[ii]== nD) {
              found = true;
              break;
            }
          }
          TERMM_IF_NOT_COND(found, to_string(a_isHalo(cellId)) + " " + to_string(a_noOffsprings(cellId)));
        }
#endif
        if(m_nghbrDomains[i] >= minNghbrDomId && m_nghbrDomains[i] <= maxNghbrDomId) {

          //TODO_SS labels:GRID do we need to put cellId into exchangeCellList in case a_noChildren(cellId)<IPOW2(nDim)
          for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
            if (treeb().solver(cellId, solverId)) {
              exchangeCellList.insert({std::make_tuple(0, cellId, i, solverId)});
            }
          }

//          if (a_noChildren(cellId)==IPOW2(nDim)) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              const MInt childId = a_childId(cellId, child);
              // Note: Due to cutOff we might have less then IPOW[nDim] children
              if (childId<0) continue;
              //TODO_SS labels:GRID,totest check the following assert: I guess during meshAdaptation is might fail
              //TERMM_IF_NOT_COND(!a_isHalo(childId), "");

              minNghbrDomId = findNeighborDomainId(a_globalId(childId));
              maxNghbrDomId =
                findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(childId) + (MLong)a_noOffsprings(childId)-1));
  //            for(MInt k = 0; k < noNeighborDomains(); k++) { }
                //if (domainId()==m_nghbrDomains[k]) continue; // check if this makes sense
                //if(m_nghbrDomains[k] >= minNghbrDomId && m_nghbrDomains[k] <= maxNghbrDomId) { }
              MInt lay = (m_nghbrDomains[i] >= minNghbrDomId && m_nghbrDomains[i] <= maxNghbrDomId) ? 0 : 1;
              // TODO_SS labels:GRID since halo data are not exchanged yet, safety approach
              if (a_isHalo(childId)) lay = 0;

              for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
                //TODO_SS labels:GRID,totest maybe we could directly check solverBits of children, but first check
                // if solverBits of halo children are already set
                if (treeb().solver(cellId, solverId)) {
                  exchangeCellList.insert({std::make_tuple(lay, childId, i, solverId)});
                  if (!a_isHalo(childId)) {
                    if (treeb().solver(childId, solverId)) {
//                      TERMM_IF_COND(lay==0 && !a_hasProperty(childId, Cell::IsPartLvlAncestor), to_string(domainId()) + " " + to_string(m_nghbrDomains[i]) + " " + to_string(minNghbrDomId) + " " + to_string(maxNghbrDomId));
                      m_windowLayer_[i][childId].set(solverId, lay);
                    }
                  }
                }
              }
            }
          }
          }
    }
  }


  //
  std::vector<std::unordered_map<MInt, M32X4bit<true>>> haloLayer_(newNghbrDomains.size());
  for(MInt layer = 0; layer < m_noHaloLayers; layer++) {

    const auto it_begin = exchangeCellList.lower_bound(std::make_tuple(layer,std::numeric_limits<MInt>::min(),
          std::numeric_limits<MInt>::min(),std::numeric_limits<MInt>::min()));
    const auto it_end = exchangeCellList.upper_bound(std::make_tuple(layer,std::numeric_limits<MInt>::max(),
          std::numeric_limits<MInt>::max(),std::numeric_limits<MInt>::max()));
    const MInt noCandidates = std::distance(it_begin, it_end);

    MInt cnt = 0;
    for (auto it = it_begin; it!=it_end;) {
      const MInt cellId = get<1>(*it);
      TERMM_IF_NOT_COND(cellId > -1 && cellId < m_tree.size(), "");
      const MBool isHaloCellId = a_isHalo(cellId);

      // 1) If cellId is on level 'level+1' it must be a window cell, and we are only looking for neighbors on same level
      // 2) If cellId is window cell, but not on 'level+1' we are only interested in neighbors on 'level+1', since we
      //    should already have tagged windows on 'level' (and those windows on 'level' without children are put into
      //    exchangeCellList above)
      // 3) If cellId is halo cell, it must be on level<='level': we allow both searching on lower and upper level
#ifndef NDEBUG
      TERMM_IF_COND(a_level(cellId)==level+1 && allWindows.find(a_parentId(cellId))==allWindows.end()
          && !a_hasProperty(a_parentId(cellId), Cell::IsPartLvlAncestor), "");

      // The following check fails when this function is called during meshAdaptation
      //TERMM_IF_COND(a_isHalo(cellId) && a_level(cellId)>level && !a_hasProperty(cellId, Cell::IsPartLvlAncestor)/*&& m_maxPartitionLevelShift==0*/, "");
#endif

      const MInt counter = a_level(cellId)==level+1 ? getAdjacentGridCells5<false,false>(cellId, nghbrList.data())
                                     : (isHaloCellId ? getAdjacentGridCells5<true,true>(cellId, nghbrList.data())
                                         : getAdjacentGridCells5<true,false>(cellId, nghbrList.data(), level+1));

      // It can happen that we reached noCandidates, but by coincidence the next cellId is the same, but on layer+1
      //while (cellId==it->first.second && cnt<noCandidates) {
      while (cnt<noCandidates && cellId==get<1>(*it)) {
        ASSERT(get<0>(*it)==layer, "Send a bug report to the C++ vendor!");
        const MInt nghbrDomIdx = get<2>(*it);
        const MInt solverId = get<3>(*it);

        for(MInt n = 0; n < counter; n++) {
          const MInt nghbrId = nghbrList[n];
          ASSERT(nghbrId > -1 && nghbrId < m_tree.size(), to_string(nghbrId) + "|" + to_string(m_tree.size()));
          ASSERT(a_level(nghbrId)<=level+1, to_string(a_level(nghbrId)) + "|" + to_string(level+1));

          // At this place the solver tags of halo cells should already be set, i.e. solverBits are already exchanged.
          if (!treeb().solver(nghbrId, solverId)) continue;

          // nghbrId is either
          // 1) a halo cell
          // 2) itself a window cell with no parent or its parent is a window cell or a window cell on level+1,
          //    which not yet in allWindows
          const MBool isHalo = a_isHalo(nghbrId);
#ifndef NDEBUG
          if (m_maxPartitionLevelShift==0) {
            TERMM_IF_COND(isHalo && (allWindows.find(a_parentId(nghbrId))!=allWindows.end()
                  || allWindows.find(nghbrId)!=allWindows.end()), "");
          }
#endif

          //TODO_SS labels:GRID,totest check the last condition, if testcases fail and if it speeds up
          if (!isHalo && a_level(nghbrId)<level+1) continue;

          // If parent is window, it must be in isSolverWindowCell(nghbrDomIdx, a_parentId(nghbrId), solverId)
          if (a_level(nghbrId)==level+1 && !a_isHalo(a_parentId(nghbrId)) && !isSolverWindowCell(nghbrDomIdx, a_parentId(nghbrId), solverId))
            continue;

          // Sanity check that cell is in m_windowCells, in case it's not a halo
#ifndef NDEBUG
          TERMM_IF_COND(!isHalo
              && allWindows.find(a_parentId(nghbrId))==allWindows.end()
              && allWindows.find(nghbrId)!=allWindows.end()
              && !a_hasProperty(a_parentId(nghbrId), Cell::IsPartLvlAncestor), "");

          // Check in case of no partition level shift
          if (m_maxPartitionLevelShift==0) {
            TERMM_IF_COND(a_level(nghbrId)==level+1 && allWindows.find(a_parentId(nghbrId))==allWindows.end(), "");
          } else {
            /* Actually I would expect for a_level(nghbrId)==level+1 the following possible situations:
             * 1) allWindows.find(a_parentId(nghbrId))!=allWindows.end() && !a_hasProperty(a_parentId(nghbrId), Cell::IsPartLvlAncestor)
             * 2) allWindows.find(a_parentId(nghbrId))!=allWindows.end() && a_hasProperty(a_parentId(nghbrId), Cell::IsPartLvlAncestor)
             * 3) allWindows.find(a_parentId(nghbrId))==allWindows.end() && a_hasProperty(a_parentId(nghbrId), Cell::IsPartLvlAncestor) && a_isHalo(a_parentId(nghbrId))
             * But since we are tagging currently way too many cells, we can have:
             * allWindows.find(a_parentId(nghbrId))==allWindows.end() && !a_hasProperty(a_parentId(nghbrId), Cell::IsPartLvlAncestor)
            */
          }
#endif

          /* Currently nghbrId could be also a halo cell, e.g. if domain c has halo cells on domain a, we have to go
            * over the halos of domain b
            *
            *  -------    -------
            *        |    |
            *    a   | b  |  c
            *        |    |
            */
            MInt nghbrDomIdx2 = nghbrDomIdx;
            if (m_maxPartitionLevelShift>0) {
              if (!isHalo) {
                const MInt ndom = newNghbrDomains[nghbrDomIdx];
                //TODO_SS labels:GRID can we ensure that the newly found domain will also find current domain?
                nghbrDomIdx2 = setNeighborDomainIndex(ndom, m_haloCells, m_windowCells, m_windowLayer_);
              }
            }
            const MInt windowLayer = isHalo ? haloLayer_[nghbrDomIdx][nghbrId].get(solverId) 
                                      : m_windowLayer_[nghbrDomIdx2][nghbrId].get(solverId);

            if (windowLayer>layer+1) {

              isHalo ? haloLayer_[nghbrDomIdx][nghbrId].set(solverId, layer+1) 
                      : m_windowLayer_[nghbrDomIdx2].at(nghbrId).set(solverId, layer+1);
              ASSERT((isHalo && haloLayer_[nghbrDomIdx][nghbrId].get(solverId)==layer+1) 
                  || (!isHalo && m_windowLayer_[nghbrDomIdx2].at(nghbrId).get(solverId)==layer+1), "");

              if (layer+1<m_noSolverHaloLayers[solverId])
                exchangeCellList.insert({std::make_tuple(layer+1, nghbrId, nghbrDomIdx, solverId)});
            }
          //}
        }
        ++cnt;
        ++it;
      }
      if (cnt==noCandidates) break;
    }
  } //loop over layers


  MInt cnt = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      MInt cellId = m_windowCells[i][j];
      if(a_level(cellId) == level) {
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          const MInt childId = a_childId(cellId, child);
          if(childId > -1) {
            if (m_windowLayer_[i].find(childId)!=m_windowLayer_[i].end())  {
              ASSERT(!m_windowLayer_[i].at(childId).all(), "");
              refineChildIds[m_maxNoChilds * cnt + child] = 1;
            }
          }
        }
      }
      cnt++;
    }
  }

  // DEBUG output
}



/** \brief Don't touch this function, because you don't know what you are doing
 * NOTE: Actually the most clever way of tagging is to go from highest level to lower levels (if we want
 *       to have the predefined number of halo layers on leaf level), but that means quite some work to do!
 * \date June 2021
 */
template <MInt nDim>
void CartesianGrid<nDim>::tagActiveWindowsOnLeafLvl3(const MInt maxLevel, const MBool duringMeshAdaptation,
                              const std::vector<std::function<void(const MInt)>>& removeCellSolver) {

  MIntScratchSpace nghbrList(27 * m_maxNoChilds, AT_, "nghbrList");

  // We don't need to exchange a_hasProperty(cellId, Cell::IsPartLvlAncestor) and a_noOffsprings(cellId)
  // here, because the halo cells on level=maxLevel are never partLvlAnc and have always
  // a_noOffsprings(cellId)==1.

#ifndef NDEBUG
  set<MInt> halos;
  // Get halo candidates
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      const MInt cellId = m_haloCells[i][j];
      halos.insert(cellId);
    }
  }
#endif

  std::unordered_map<MInt, M32X4bit<true>> haloLayers;

  // Note: when this is called during meshAdaptation, then changes will only occur above m_maxUniformRefinementLevel
  for (MInt level = m_minLevel; level <= maxLevel; level++) {
    std::set<std::tuple<MInt/*halolayer*/,MInt/*cellId*/,MInt/*solverId*/>> exchangeCellList;

  // Find first layer window cells
  for (MInt i = 0; i < noNeighborDomains(); i++) {
    for (auto& item : m_windowLayer_[i]) {
      const MInt cellId = item.first;

      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        // TODO_SS labels:GRID what about a_hasChildren(cellId, solverId)<IPOW2(nDim)
        if (a_level(cellId)==level || (a_level(cellId)==level-1 && !a_hasChildren(cellId, solverId))) {
          if (item.second.get(solverId)<=1) {
            ASSERT(m_tree.solver(cellId, solverId), "Cell was identifed to be a valid windowCell for"
                " solver=" + to_string(solverId) + ", but does not belong to that solver");
            exchangeCellList.insert({std::make_tuple(0, cellId, solverId)});
          }
        }
      }
    }

    // Find halo cells on level-1 with no children
    for (const auto& item : haloLayers) {
      const MInt cellId = item.first;
      if (a_level(cellId)==level-1) {
        for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
          const MInt haloLayer = item.second.get(solverId);
          if (!a_hasChildren(cellId, solverId) && haloLayer<m_noSolverHaloLayers[solverId]) {
            exchangeCellList.insert({std::make_tuple(haloLayer, cellId, solverId)});
            // TODO_SS labels:GRID what about a_hasChildren(cellId, solverId)<IPOW2(nDim)
          } /*else if (a_noChildren(cellId)<IPOW2(nDim) && haloLayer<=m_noSolverHaloLayers[solverId]) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;
              if (treeb().solver(childId, solverId)) {
                haloLayers[childId].set(solverId, haloLayer);
                if (haloLayer<m_noSolverHaloLayers[solverId]) {
                  exchangeCellList.insert({std::make_pair(haloLayer, childId), solverId});
                }
              }
            }
          }*/
        }
      }
    }

    // In the following halo cells with children on current domain are also identified as 1st layer 'window' cells
    if (m_maxPartitionLevelShift>0) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        const MInt cellId = m_haloCells[i][j];
        if (a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
          const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          const MInt maxNghbrDomId =
            findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
          if (domainId()<minNghbrDomId || domainId()>maxNghbrDomId) continue;
          for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
            if (m_tree.solver(cellId, solverId) && a_hasChildren(cellId, solverId)) {
              // Guess: The folowing assert fails e.g. if partitionCellMaxNoOffspring==1 and one solver is at least one
              //        level coarser then the other one.
//              TERMM_IF_NOT_COND(a_hasChildren(cellId, solverId), "");
              //TODO_SS labels:GRID what about noChildren<IPOW2(nDim)
              if (a_level(cellId)==level || (a_level(cellId)==level-1 && !a_hasChildren(cellId, solverId))) {
                exchangeCellList.insert({std::make_tuple(0, cellId, solverId)});
                // It may happen that one domain has no window cell at this level at all so we add this cell to
                // haloLayers with layer==0
                haloLayers[cellId].set(solverId,0);
              }
            }
          }
        }
      }
    }
  }


  // In the following no distinction is made regarding the neighbor. The current domain is extended by the
  // required number of halo layers, and halo cells required by current domain are inserted in haloLayers.
  for(MInt layer = 0; layer < m_noHaloLayers; layer++) {

    const auto it_begin = exchangeCellList.lower_bound(std::make_tuple(layer,std::numeric_limits<MInt>::min(),
          std::numeric_limits<MInt>::min()));
    const auto it_end = exchangeCellList.upper_bound(std::make_tuple(layer,std::numeric_limits<MInt>::max(),
          std::numeric_limits<MInt>::min()));
    const MInt noCandidates = std::distance(it_begin, it_end);

    MInt cnt = 0;
    for (auto it = it_begin; it!=it_end;) {
      const MInt cellId = get<1>(*it);
      ASSERT(cellId > -1 && cellId < m_tree.size(), "");

      const MInt counter = a_level(cellId)==level ? getAdjacentGridCells5<false,false>(cellId, nghbrList.data())
                                         : getAdjacentGridCells5<true,false>(cellId, nghbrList.data(), level);
      // TEST1
//      const MInt counter = a_level(cellId)==level ? getAdjacentGridCells5<false,true>(cellId, nghbrList.data())
//                                         : getAdjacentGridCells5<true,false>(cellId, nghbrList.data());


      while (cnt<noCandidates && cellId==get<1>(*it)) {
        ASSERT(get<0>(*it)==layer, "Send a bug report to the C++ vendor!");
        const MInt solverId = get<2>(*it);

        //TODO_SS labels:GRID the case hat one of the neighbors has a_noChildren<IPOW2(nDim)
        for(MInt n = 0; n < counter; n++) {
          const MInt nghbrId = nghbrList[n];
          ASSERT(nghbrId > -1 && nghbrId < m_tree.size(), "");

          //TEST1
/*          TERMM_IF_NOT_COND(a_level(nghbrId)==level || a_level(nghbrId)==level-1, "");
          if (!a_isHalo(nghbrId) || !treeb().solver(nghbrId, solverId))
            continue;
          if (a_level(nghbrId)==level-1 && !a_hasChildren(nghbrId,solverId)) continue;
          if (a_level(nghbrId)==level-1) {
            TERMM_IF_NOT_COND(a_noChildren(nghbrId)<IPOW2(nDim) || a_noChildren(nghbrId)==0, "");
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(nghbrId, child);
              if(childId < 0) continue;

              TERMM_IF_NOT_COND(treeb().solver(childId,solverId) && a_isHalo(childId), "");
              const MInt haloLayer = haloLayers[childId].get(solverId);

              if (haloLayer>layer+1) {

                haloLayers[childId].set(solverId, layer+1);

                if (layer+1<m_noSolverHaloLayers[solverId])
                  exchangeCellList.insert({std::make_pair(layer+1, childId), solverId});
              }
            }
            continue;
          }*/

#ifndef NDEBUG
          // Sanity checks
          TERMM_IF_NOT_COND(a_level(nghbrId)==level, "");
          // It might happen that a cell is disabled, but the isHalo flag is still set, that' why the a_isHalo test is not reliable
          TERMM_IF_NOT_COND(((halos.find(nghbrId)!=halos.end())==(findNeighborDomainId(a_globalId(nghbrId))!=domainId()))
                            || a_hasProperty(nghbrId, Cell::IsPeriodic) || !treeb().solver(nghbrId, solverId), "");
#endif

          // Skip if cell is not halo
          if (!a_isHalo(nghbrId) || !treeb().solver(nghbrId, solverId)/* || a_hasProperty(nghbrId, Cell::IsPartLvlAncestor)*/)
            continue;

          const MInt haloLayer = haloLayers[nghbrId].get(solverId);

          if (haloLayer>layer+1) {

            haloLayers[nghbrId].set(solverId, layer+1);
            ASSERT(haloLayers[nghbrId].get(solverId)==layer+1, "");

            if (layer+1<m_noSolverHaloLayers[solverId])
              exchangeCellList.insert({std::make_tuple(layer+1, nghbrId, solverId)});
          }
        }
        ++cnt;
        ++it;
      }
      if (cnt==noCandidates)
        break;
    }
  } //loop over layers
  }


  // Traverse through the tree to mark all required halos on all levels
  std::unordered_map<MInt,M8X1bit<false>> validHalos;
  for (auto it = haloLayers.cbegin(); it!=haloLayers.cend(); ++it) {
    const MInt cellId = it->first;
    ASSERT(validHalos.find(cellId)==validHalos.end(), "Duplicate halo cells!");
    auto& flag = validHalos[cellId];
    for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
      if (it->second.get(solverId)<=m_noSolverHaloLayers[solverId]) {
        ASSERT(m_tree.solver(cellId,solverId), "");
        flag.set(solverId,1);
      }
    }
  }

//#ifndef NDEBUG
  // Sanity check that all parents of valid halo cells are also marked as valid
  for (const auto& item : validHalos) {
    MInt parentId = item.first;
    while (a_parentId(parentId)>-1) {
      parentId = a_parentId(parentId);
      const auto it = validHalos.find(parentId);

      // Case 1: noPartLvlShift, then parent must be halo, but it can reside on current domain in periodic case
      // Case 2: partLvlShift, then parent can be window
      TERMM_IF_NOT_COND(m_maxPartitionLevelShift>0 || (a_isHalo(parentId) && it!=validHalos.end()), "");
      TERMM_IF_NOT_COND(m_maxPartitionLevelShift>0 || a_hasProperty(parentId, Cell::IsPeriodic) || findNeighborDomainId(a_globalId(parentId))!=domainId(), "");
      TERMM_IF_NOT_COND(it!=validHalos.end() || (a_hasProperty(parentId, Cell::IsPartLvlAncestor) && !a_isHalo(parentId)), "");

      if (it==validHalos.end()) break;

      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        TERMM_IF_COND(item.second.get(solverId) && !it->second.get(solverId), " solverId=" + to_string(solverId) +
            ": Parent of valid halo cell is no halo!");
      }
    }
  }
//#endif

  // ---Disable all halos not required & communicate with respective windows --- //

  // Receive information from halo cells, whether the respective window cell is required
  ScratchSpace<MPI_Request> recvRequests(max(1, noNeighborDomains()), AT_, "recvRequests");
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  MInt receiveCount = 0;
  for (const auto& vecWindow : m_windowCells) receiveCount+=vecWindow.size();
  ScratchSpace<M8X1bit<>::type> windowBuffer(max(1, receiveCount), AT_, "windowBuffer");
  for (MInt i = 0, offset = 0; i < noNeighborDomains(); i++) {
    const MInt noWindowCells = m_windowCells[i].size();
    if (noWindowCells<1) continue;

    MPI_Irecv(&windowBuffer[offset], noWindowCells, type_traits<M8X1bit<>::type>::mpiType(), m_nghbrDomains[i],
        m_nghbrDomains[i], mpiComm(), &recvRequests[i], AT_, "windowBuffer[offset]");

    offset += noWindowCells;
  }

  // The meshAdaptation forces us to keep children of pla's
//  const MBool hasPartLvlShift = (m_maxPartitionLevelShift > 0);
  const MInt tmpScratchSize = m_tree.size();//(hasPartLvlShift) ? m_tree.size() : 1;
  ScratchSpace<MBool> isHaloPartLvlAncestor(tmpScratchSize, AT_, "isHaloPartLvlAncestor");
  isHaloPartLvlAncestor.fill(false);

  // Determine relevant halo partition level ancestor window/halo cells that need to be preserved
  if(m_maxPartitionLevelShift > 0) {
    for(auto& i : m_partitionLevelAncestorIds) {
      ASSERT(a_hasProperty(i, Cell::IsPartLvlAncestor),
                        "cell is not a partition level ancestor: " + std::to_string(i));
      // Mark halo partition level ancestors (which are part of the missing subtree of the grid)
      if(a_hasProperty(i, Cell::IsHalo)) {
        isHaloPartLvlAncestor[i] = true;
      }
      // Mark all halo childs of partition level ancestors (required for exchange of e.g. global
      // ids!)
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(i, child);
        if(childId > -1 && a_isHalo(childId)) {
          isHaloPartLvlAncestor[childId] = true;
        }
      }
    }
  }
  // Should we keep all isHaloPartLvlAncestor cells, because they are kept in meshAdaptation?!

  struct comp {
    MBool operator()(const std::pair<MInt, MInt>& m1, const std::pair<MInt, MInt>& m2) const {
      if(m1.first < m2.first) {
        return true;
      } else if (m1.first==m2.first) {
        if(m1.second < m2.second) {
          return true;
        } else {
          return false;
        }
      } else {
        return false;
      }
    }
  };

  const MInt threesholdLvl = !duringMeshAdaptation ? 0/*m_minLevel*/ : m_maxUniformRefinementLevel;
  ScratchSpace<MPI_Request> sendRequests(max(1, noNeighborDomains()), AT_, "sendRequests");
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  MInt sendCount = 0;
  for (const auto& vecHalo : m_haloCells) sendCount+=vecHalo.size();
  ScratchSpace<M8X1bit<>::type> tmp_data(sendCount, AT_, "tmp_data");
  std::map<std::pair<MInt/*level*/,MInt/*cellId*/>,M8X1bit<false>,comp> toDelete;
  for (MInt i = 0, idx = 0; i < noNeighborDomains(); i++) {
    if (m_maxPartitionLevelShift==0) toDelete.clear();
    const MInt offset = idx;
    const MInt noHaloCells = m_haloCells[i].size();
    if (noHaloCells<1) continue;
    std::vector<MInt> haloCells;
    haloCells.reserve(noHaloCells);

    // Temporary fix: In case two halos are mapped to same window cell, keep them (occurs in periodic case)
    ScratchSpace<MBool> isDuplicateHalo(m_noPeriodicCartesianDirs>0 ? m_tree.size() : 1, AT_, "isDuplicateHalo");
    isDuplicateHalo.fill(false);
    if (m_noPeriodicCartesianDirs>0) {
      std::map<MLong,MInt> haloGlobalIds;
      for (const auto cellId : m_haloCells[i]) {
        const MLong globalId = a_globalId(cellId);
        if (haloGlobalIds.find(globalId)!=haloGlobalIds.end()) {
          isDuplicateHalo[haloGlobalIds[globalId]] = true;
          isDuplicateHalo[cellId] = true;
        } else
          haloGlobalIds.insert({std::make_pair(globalId, cellId)});
      }
    }

    for (const auto cellId : m_haloCells[i]) {
      const MBool isValidHalo = validHalos.find(cellId)!=validHalos.end();

#ifndef NDEBUG
      // Check that partLvlAnc halo cell, with children on current domain is identified as validHalo
      const MInt minNghbrDomId1 = findNeighborDomainId(a_globalId(cellId));
      const MInt maxNghbrDomId1 =
        findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
      if (domainId()>=minNghbrDomId1 && domainId()<=maxNghbrDomId1 && !a_hasProperty(cellId, Cell::IsPeriodic))
        TERMM_IF_NOT_COND(!a_hasProperty(cellId, Cell::IsPartLvlAncestor) || isValidHalo, "Cell is partLvlAncestor, but not found as valid halo!");
#endif

      // TODO_SS labels:GRID,toenhance do it more efficiently
      // Don't delete anything on minLevel
      const M8X1bit<> solver = (a_level(cellId)<=threesholdLvl
                             || isHaloPartLvlAncestor[cellId]
                             || (m_noPeriodicCartesianDirs>0 && isDuplicateHalo[cellId]))
        ? M8X1bit<false>{M8X1bit<true>{}.data()} : (isValidHalo ? validHalos[cellId] : M8X1bit<false>{});

      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        if (m_tree.solver(cellId,solverId) && !solver.get(solverId)) {
          auto key = std::make_pair(-a_level(cellId),cellId);
          toDelete[key].set(solverId,true);
        }
        if (a_level(cellId)>threesholdLvl && !isHaloPartLvlAncestor[cellId] && m_tree.solver(cellId,solverId))
          m_tree.solver(cellId,solverId) = solver.get(solverId);
      }
      tmp_data[idx++] = solver.data();
      if (m_tree.solverBits(cellId).any() || isHaloPartLvlAncestor[cellId]) {
        haloCells.push_back(cellId);
      } else {
//        a_hasProperty(cellId, Cell::IsHalo) = false;
//        a_isToDelete(cellId)=true;
//        m_freeIndices.insert(cellId);
//        removeCell<false>(cellId);
//        ASSERT(toDelete.find(std::make_pair(-a_level(cellId),cellId))!=toDelete.end(), "");
        //TODO_SS labels:GRID,toenhance Check if it might make sense, to set  M8X1bit<false>(M8X1bit<true>.data()), to force
        //         the solvers to delete the cells, if they exist in the solvers.
        // We need the following because, we might have partLvlAnc cell belonging to a solver,
        // but one of its children does not. So if we just propagate the informations from parent cell
        // we might have created halo cells, which are not used by any solver. We need to delete them.
        toDelete.insert({std::make_pair(-a_level(cellId),cellId), M8X1bit<false>()});
      }
    }
    m_haloCells[i] = haloCells;
    if (m_maxPartitionLevelShift==0) {
      // In case of partLvlShift we need to first loop over all neigbhors, because a cell marked for deletion
      // might still have children, which belong to a different neighbor and are therefore not deleted yet
      for (const auto& item : toDelete) {
        const MInt cellId = item.first.second;
        for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
          if (!removeCellSolver.empty() && item.second.get(solverId)) {
            removeCellSolver[solverId](cellId);
          }
        }
        if (m_tree.solverBits(cellId).none()) {
          ASSERT(!isHaloPartLvlAncestor[cellId], "");
          removeCell<false>(cellId);

        }
      }
    }

    ASSERT(idx-offset==noHaloCells, "");

    MPI_Isend(&tmp_data[offset], noHaloCells, type_traits<M8X1bit<>::type>::mpiType(), m_nghbrDomains[i], domainId(),
        mpiComm(), &sendRequests[i], AT_, "tmp_data");
  }
  if (m_maxPartitionLevelShift>0) {
    for (const auto& item : toDelete) {
      const MInt cellId = item.first.second;
      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        if (!removeCellSolver.empty() && item.second.get(solverId)) {
          removeCellSolver[solverId](cellId);
        }
      }
      if (m_tree.solverBits(cellId).none()) {
        removeCell<false>(cellId);
      }
    }
  }

  // Finish MPI communication
  MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  for (MInt i = 0, idx = 0; i < noNeighborDomains(); i++) {
    // In periodic cases some window cells may appear twice; we have ensured above that those cells are always kept
    std::vector<MInt> windowCells;
    windowCells.reserve(m_windowCells[i].size());
    for (const auto cellId : m_windowCells[i]) {
      ASSERT(m_windowLayer_[i].find(cellId)!=m_windowLayer_[i].end(), "");
      M8X1bit<> solver(windowBuffer[idx++]);
      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        if (!solver.get(solverId) && isSolverWindowCell(i,cellId,solverId))
          m_windowLayer_[i].at(cellId).set(solverId, m_noHaloLayers+1);
      }
      MBool validWindow = false;
      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        if (isSolverWindowCell(i,cellId,solverId)) {
          windowCells.push_back(cellId);
          validWindow=true;
          break;
        }
      }
      ASSERT(solver.any()==validWindow, "");
      // If a window cell was already checked and not deleted, we are not allowed to delete it now
      if (!validWindow) {
        //TODO_SS labels:GRID periodic cases
        const MBool ok = m_windowLayer_[i].erase(cellId);
        ASSERT(ok, "");
        MBool isWindow = false;
        for (const auto& w : m_windowLayer_) {
          if (w.find(cellId)!=w.end()) {
            isWindow = true;
            break;
          }
        }
        if (!isWindow)
          a_hasProperty(cellId, Cell::IsWindow) = false;
      }
    }
    m_windowCells[i] = windowCells;
  }

  if  (!duringMeshAdaptation)
    compactCells();

  // DEBUG output
}


/** \brief
 * Actually the most clever way of tagging is to go from highest level to lower levels, but that
 * means quite some work to do!
    \date June 2021
  */
template <MInt nDim>
void CartesianGrid<nDim>::tagActiveWindowsAtMinLevel(const MInt level) {
  TRACE();

  MIntScratchSpace nghbrList(27 * m_maxNoChilds, AT_, "nghbrList");

  // We rely on the fact, that exchangeCellList is sorted by the tuple's first element, then by its 2nd element etc.
  std::set<std::tuple<MInt/*layer*/,MInt/*cellId*/, MInt/*nghbrIdx*/,MInt/*solverId*/>> exchangeCellList;
  std::vector<std::unordered_map<MInt, M32X4bit<true>>> haloLayer_(noNeighborDomains());

  ASSERT(m_windowLayer_.empty(), "");
  m_windowLayer_.resize(noNeighborDomains());

  if(m_maxPartitionLevelShift > 0) {
    // exchange some window cell data
    ScratchSpace<MInt> bprops(m_tree.size(), 2, AT_, "bprops");
    bprops.fill(-1);
    // Note: changed loop to iterate over whole tree, internal and halo cells not sorted!
    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      if(a_isHalo(cellId) || a_isToDelete(cellId)) {
        continue;
      }
      bprops(cellId, 0) = (MInt)a_hasProperty(cellId, Cell::IsPartLvlAncestor);
      bprops(cellId, 1) = a_noOffsprings(cellId);

      ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
             std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                 + std::to_string(bprops(cellId, 1)));
    }
    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), bprops.getPointer(), 2);

    // Note: loop only over the current haloCells since only for these data is exchanged, this does
    // not work properly if in case of a partition level shift some halo cell is not part of
    // m_haloCells yet such that its property IsPartLvlAncestor is reset!
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        const MInt cellId = m_haloCells[i][j];
        ASSERT(!a_isToDelete(cellId), "");
        ASSERT(bprops(cellId, 0) > -1 && bprops(cellId, 1) > 0,
               std::to_string(cellId) + " " + std::to_string(bprops(cellId, 0)) + " "
                   + std::to_string(bprops(cellId, 1)));
        a_hasProperty(cellId, Cell::IsPartLvlAncestor) = (MBool)bprops(cellId, 0);
        a_noOffsprings(cellId) = bprops(cellId, 1);
      }
    }
  }

  // Get halo candidates
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      const MInt cellId = m_haloCells[i][j];
      ASSERT(a_level(cellId)==level, "We expect all halo cells at this point to be on minLevel");

      // Note: currently we add all halos as candidates. A more sophisticated choice would be to add only those
      //       who are adajcent to cells, which are not halos of same neighbor (the neighbor can still be a halo)
      for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
        // For the following to work we need to first exchange the solverBits; Check if it speeds up the whole stuff.
        // if (treeb().solver(cellId, solverId) || treeb().noSolvers()==1)
        exchangeCellList.insert({std::make_tuple(0, cellId, i, solverId)});
      }

      if(m_maxPartitionLevelShift > 0) {
        if(a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
          const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
          const MInt maxNghbrDomId =
              findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
          for(MInt k = 0; k < noNeighborDomains(); k++) {
            if(m_nghbrDomains[k] >= minNghbrDomId && m_nghbrDomains[k] <= maxNghbrDomId) {
              // Note: currently we add all halos as candidates. A more sophisticated choice would be to add only those
              //       who are adajcent to cells, which are not halos of same neighbor (the neighbor can still be a halo)
              for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
                // For the following to work we need to first exchange the solverBits; Check if it speeds up the whole stuff.
                // if (treeb().solver(cellId, solverId) || treeb().noSolvers()==1)
                exchangeCellList.insert({std::make_tuple(0, cellId, k, solverId)});
              }
            }
          }
        }
      }

    }
  }

  vector<std::set<MInt>> windCellSet(noNeighborDomains());
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    windCellSet[i].insert(m_windowCells[i].begin(), m_windowCells[i].end());
  }
  //
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      const MInt cellId = m_windowCells[i][j];

      if (a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
        const MInt minNghbrDomId = findNeighborDomainId(a_globalId(cellId));
        const MInt maxNghbrDomId =
          findNeighborDomainId(mMin(m_noCellsGlobal - 1, a_globalId(cellId) + (MLong)a_noOffsprings(cellId)-1));
#ifndef NDEBUG
        for (MInt nD = minNghbrDomId; nD<=maxNghbrDomId; ++nD) {
          if (nD==domainId()) continue;
          MBool found = false;
          for (MInt ii = 0; ii < noNeighborDomains(); ii++) {
            if (m_nghbrDomains[ii]== nD) {
              found = true;
              break;
            }
          }
          TERMM_IF_NOT_COND(found, to_string(a_isHalo(cellId)) + " " + to_string(a_noOffsprings(cellId)));
        }
#endif
        for (MInt solverId = 0; solverId < treeb().noSolvers(); solverId++) {
          if (treeb().solver(cellId, solverId)) {
//            for(MInt k = 0; k < noNeighborDomains(); k++) {
              //if (domainId()==m_nghbrDomains[k]) continue; // check if this makes sense
              if(m_nghbrDomains[i/*k*/] >= minNghbrDomId && m_nghbrDomains[i/*k*/] <= maxNghbrDomId) {
                TERMM_IF_NOT_COND(windCellSet[i].find(cellId)!=windCellSet[i].end(), "");
                m_windowLayer_[i/*k*/][cellId].set(solverId, 0);
                exchangeCellList.insert({std::make_tuple(0, cellId, i/*k*/, solverId)});
              }
//            }
          }
        }
      }
    }
  }



  //
  for(MInt layer = 0; layer < m_noHaloLayers; layer++) {

    const auto it_begin = exchangeCellList.lower_bound(std::make_tuple(layer,std::numeric_limits<MInt>::min(),
          std::numeric_limits<MInt>::min(), std::numeric_limits<MInt>::min()));
    const auto it_end = exchangeCellList.upper_bound(std::make_tuple(layer+1,std::numeric_limits<MInt>::min(),
          std::numeric_limits<MInt>::min(), std::numeric_limits<MInt>::min()));
    const MInt noCandidates = std::distance(it_begin, it_end);

    MInt cnt = 0;
    for (auto it = it_begin; it!=it_end;) {
      ASSERT(get<0>(*it)==layer, "Send a bug report to the C++ vendor!");
      const MInt cellId = get<1>(*it);

      // Only look for nghbrs on same level
      const MInt counter = getAdjacentGridCells5<false,false>(cellId, nghbrList.data());

      while (cnt<noCandidates && cellId==get<1>(*it)) {
        ASSERT(get<0>(*it)==layer, "Send a bug report to the C++ vendor!");
        const MInt nghbrDomIdx = get<2>(*it);
        ASSERT((signed)m_windowLayer_.size()>nghbrDomIdx, "");
        const MInt solverId = get<3>(*it);

        for(MInt n = 0; n < counter; n++) {
          const MInt nghbrId = nghbrList[n];
          ASSERT(a_level(nghbrId)==level, "Nghbr expected to be on same level");

          // if 'treeb().solver(nghbrId, solverId)' is run before the bitsets are exhchanged, it will always result
          // in false for halo cells.
          // nghbrId must be either haloCell or window cell (at least for m_maxPartitionLevelShift==0).
          // Since solverBits haven't yet been exchanged, we cannot check solver affiliation of halo cells
          const MBool isHalo = a_isHalo(nghbrId);
          if (!isHalo && !treeb().solver(nghbrId, solverId)) continue;

          // Sanity check that cell is in m_windowCells (with partLvlShift this check doesn't work always)
          if (!isHalo && windCellSet[nghbrDomIdx].find(nghbrId)==windCellSet[nghbrDomIdx].end() && m_maxPartitionLevelShift>0) continue;
#ifndef NDEBUG
          TERMM_IF_COND(!isHalo && windCellSet[nghbrDomIdx].find(nghbrId)==windCellSet[nghbrDomIdx].end(), "This must be a window cell!");
#endif

          MInt windowLayer = isHalo ? haloLayer_[nghbrDomIdx][nghbrId].get(solverId) 
                                    : m_windowLayer_[nghbrDomIdx][nghbrId].get(solverId);
          if (windowLayer>layer+1) {

            isHalo ? haloLayer_[nghbrDomIdx][nghbrId].set(solverId, layer+1) 
                   : m_windowLayer_[nghbrDomIdx].at(nghbrId).set(solverId, layer+1);
            ASSERT((isHalo && haloLayer_[nghbrDomIdx][nghbrId].get(solverId)==layer+1) 
                || (!isHalo && m_windowLayer_[nghbrDomIdx].at(nghbrId).get(solverId)==layer+1), "");

            if (layer+1<m_noSolverHaloLayers[solverId])
              exchangeCellList.insert({std::make_tuple(layer+1, nghbrId, nghbrDomIdx, solverId)});
          }
        }
        ++cnt;
        ++it;
      }
      if (cnt==noCandidates) break;
    }
  } //loop over layers

  // All cells which are in m_windowCells but still not in m_windowLayer_ are appended to the latter (I don't know if
  // it's necessary or clever to do it)
  for (MInt i = 0; i < noNeighborDomains(); ++i) {
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      const MInt cellId = m_windowCells[i][j];
      ASSERT(!a_isHalo(cellId), "");
      if (m_windowLayer_[i].find(cellId)==m_windowLayer_[i].end()) {
        for (MInt solverId = 0; solverId < treeb().noSolvers(); ++solverId)
          if (m_tree.solver(cellId, solverId)) {
            // My guess: actually it should be sufficient to have m_noHaloLayers+1, but since the respective halos are
            // already created, we need to accept all cells in m_windowCells, which are not yet in m_windowLayer_,
            // to have matching buffer sizes in MPI communications
            m_windowLayer_[i][cellId].set(solverId, a_hasProperty(cellId, Cell::IsPartLvlAncestor) ? 1 : (m_haloMode==2) ? m_noHaloLayers+1 : m_noSolverHaloLayers[solverId]);
          }
      }
    }
  }

  // DEBUG output
}
// clang-format on
//-------------------------------------------------------------------------------------------


/**
 * \brief Retrieves all direct and diagonal neighboring cells of the given cell on the child level if available
 * \author Lennart Schneiders
 */
template <MInt nDim>
MInt CartesianGrid<nDim>::getAdjacentGridCells(MInt cellId, MIntScratchSpace& adjacentCells, MInt level,
                                               MBool diagonalNeighbors) {
  set<MInt> nghbrs;
  for(MInt dir0 = 0; dir0 < m_noDirs; dir0++) {
    MInt nghbrId0 = -1;
    if(a_hasNeighbor(cellId, dir0) > 0)
      nghbrId0 = a_neighborId(cellId, dir0);
    else if(a_parentId(cellId) > -1) {
      if(a_hasNeighbor(a_parentId(cellId), dir0) > 0) {
        nghbrId0 = a_neighborId(a_parentId(cellId), dir0);
      }
    }

    if(nghbrId0 < 0) continue;
    if(a_noChildren(nghbrId0) > 0 && a_level(nghbrId0) <= level) {
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        if(!childCode[dir0][child]) continue;
        if(a_childId(nghbrId0, child) > -1) nghbrs.insert(a_childId(nghbrId0, child));
      }
    } else {
      nghbrs.insert(nghbrId0);
    }

    if(diagonalNeighbors) {
      for(MInt dir1 = 0; dir1 < m_noDirs; dir1++) {
        if((dir1 / 2) == (dir0 / 2)) continue;
        MInt nghbrId1 = -1;
        if(a_hasNeighbor(nghbrId0, dir1) > 0)
          nghbrId1 = a_neighborId(nghbrId0, dir1);
        else if(a_parentId(nghbrId0) > -1)
          if(a_hasNeighbor(a_parentId(nghbrId0), dir1) > 0) nghbrId1 = a_neighborId(a_parentId(nghbrId0), dir1);
        if(nghbrId1 < 0) continue;
        if(a_noChildren(nghbrId1) > 0 && a_level(nghbrId1) <= level) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            if(!childCode[dir0][child] || !childCode[dir1][child]) continue;
            if(a_childId(nghbrId1, child) > -1) nghbrs.insert(a_childId(nghbrId1, child));
          }
        } else
          nghbrs.insert(nghbrId1);
        IF_CONSTEXPR(nDim == 3) {
          for(MInt dir2 = 0; dir2 < m_noDirs; dir2++) {
            if(((dir2 / 2) == (dir0 / 2)) || ((dir2 / 2) == (dir1 / 2))) continue;
            MInt nghbrId2 = -1;
            if(a_hasNeighbor(nghbrId1, dir2) > 0)
              nghbrId2 = a_neighborId(nghbrId1, dir2);
            else if(a_parentId(nghbrId1) > -1)
              if(a_hasNeighbor(a_parentId(nghbrId1), dir2) > 0) nghbrId2 = a_neighborId(a_parentId(nghbrId1), dir2);
            if(nghbrId2 < 0) continue;
            if(a_noChildren(nghbrId2) > 0 && a_level(nghbrId2) <= level) {
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                if(!childCode[dir0][child] || !childCode[dir1][child] || !childCode[dir2][child]) continue;
                if(a_childId(nghbrId2, child) > -1) nghbrs.insert(a_childId(nghbrId2, child));
              }
            } else
              nghbrs.insert(nghbrId2);
          }
        }
      }
    }
  }


  MInt cnt = 0;
  for(MInt nghbr : nghbrs) {
    ASSERT(a_level(nghbr) <= level + 1, "");
    ASSERT(cnt < (signed)adjacentCells.size(), "");
    adjacentCells[cnt] = nghbr;
    cnt++;
  }
  return cnt;
}

// Object oriented version of getAdjacent... which is fast as a bullet
// same level <false, false>
// one lower or/and one upper allowed
template <MInt nDim>
template <MBool finer, MBool coarser>
inline MInt CartesianGrid<nDim>::getAdjacentGridCells5(const MInt cellId, MInt* const adjacentCells, const MInt level,
                                                       const MBool diagonalNeighbors) {
  static constexpr const MInt dimConverter[3][2] = {{1, 2}, {0, 2}, {0, 1}};
  set<MInt> nghbrs;

  if(diagonalNeighbors) {
    for(MInt dim = 0; dim < nDim; ++dim) {
      if(nDim == 2)
        getAdjacentGridCells1d5<finer, coarser>(nghbrs, cellId, dim, dimConverter[dim][0]);
      else if(nDim == 3)
        getAdjacentGridCells1d5<finer, coarser>(nghbrs, cellId, dim, dimConverter[dim][0], dimConverter[dim][1]);
    }
  } else {
    for(MInt dim = 0; dim < nDim; ++dim) {
      getAdjacentGridCells1d5<finer, coarser>(nghbrs, cellId, dim);
    }
  }

  ASSERT(nghbrs.size() <= 27 * m_maxNoChilds, "");

  if(!finer && !coarser) {
    // same level
    const MInt noNghbrs = nghbrs.size();
    std::copy(nghbrs.begin(), nghbrs.end(), adjacentCells);
    return noNghbrs;
  } else if(level == -1) {
    // all levels
    const MInt noNghbrs = nghbrs.size();
    std::copy(nghbrs.begin(), nghbrs.end(), adjacentCells);
    return noNghbrs;
  } else {
    MInt cnt = 0;
    for(auto cellId_ : nghbrs) {
      if(a_level(cellId_) == level) {
        adjacentCells[cnt++] = cellId_;
      }
    }
    return cnt;
  }
}

/**
 * \brief
 * \author me
 */
template <MInt nDim>
template <MBool finer, MBool coarser>
inline void CartesianGrid<nDim>::getAdjacentGridCells1d5(set<MInt>& nghbrs, const MInt cellId, const MInt dim,
                                                         const MInt dimNext, const MInt dimNextNext,
                                                         const uint_fast8_t childCodes) {
  TRACE();

  for(MInt dir = 2 * dim; dir < 2 * dim + 2; ++dir) {
    auto childCodes_ = childCodes & childCodePro[dir];
    if(a_hasNeighbor(cellId, dir) > 0) {
      const MInt nghbrId0 = a_neighborId(cellId, dir);
      if(finer && a_noChildren(nghbrId0) > 0) {
        // The neighbor might have children, but which are not adjacent
        MBool childFound = false;
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          if(/*!(childCodes_ & (1<<child))*/ !childCode[dir][child]) continue;
          if(a_childId(nghbrId0, child) > -1) {
            const MInt childId = a_childId(nghbrId0, child);
            // TODO_SS labels:GRID Can this cell have children and what to do in that case?!
            // if (nghbrs.find(childId)!=nghbrs.end()) continue;
            nghbrs.insert(childId);
            childFound = true;
            if(dimNext > -1)
              getAdjacentGridCells1d5<false, true>(nghbrs, childId, dimNext, dimNextNext, -1, childCodes_);
            if(dimNextNext > -1)
              getAdjacentGridCells1d5<false, true>(nghbrs, childId, dimNextNext, dimNext, -1, childCodes_);
          }
        }
        if(!childFound) { // TODO_SS labels:GRID,totest Check if this is required
          // if (nghbrs.find(nghbrId0)!=nghbrs.end()) continue;
          nghbrs.insert(nghbrId0);
          if(dimNext > -1)
            getAdjacentGridCells1d5<finer, coarser>(nghbrs, nghbrId0, dimNext, dimNextNext, -1, childCodes_);
          if(dimNextNext > -1)
            getAdjacentGridCells1d5<finer, coarser>(nghbrs, nghbrId0, dimNextNext, dimNext, -1, childCodes_);
        }
      } else {
        // if (nghbrs.find(nghbrId0)!=nghbrs.end()) continue;
        nghbrs.insert(nghbrId0);
        if(dimNext > -1)
          getAdjacentGridCells1d5<finer, coarser>(nghbrs, nghbrId0, dimNext, dimNextNext, -1, childCodes_);
        if(dimNextNext > -1)
          getAdjacentGridCells1d5<finer, coarser>(nghbrs, nghbrId0, dimNextNext, dimNext, -1, childCodes_);
      }
    } else if(coarser && a_parentId(cellId) > -1) {
      if(a_hasNeighbor(a_parentId(cellId), dir) > 0) {
        const MInt nghbrId0 = a_neighborId(a_parentId(cellId), dir);
        // TODO_SS labels:GRID Can this cell have children and what to do in that case?!
        // if (nghbrs.find(nghbrId0)!=nghbrs.end()) continue;
        nghbrs.insert(nghbrId0);
        // TODO_SS labels:GRID,totest check the following
        // TERMM_IF_COND(a_hasChildren(nghbrId0), "");
        if(dimNext > -1) getAdjacentGridCells1d5<true, false>(nghbrs, nghbrId0, dimNext, dimNextNext, -1, childCodes_);
        if(dimNextNext > -1)
          getAdjacentGridCells1d5<true, false>(nghbrs, nghbrId0, dimNextNext, dimNext, -1, childCodes_);
      }
    }
  }
}

//-------------------------------------------------------------------------


/** Refine the cells with the given ids.
 *
 * \author Lennart Schneiders
 * \date 12.12.2012
 *
 * This function only creates the geometrical data for new cells,
 * the initialization of flow variables must be dealt with by the
 * corresponding solver classes.
 *
 */
template <MInt nDim>
void CartesianGrid<nDim>::refineCell(
    const MInt cellId, const MLong* const refineChildIds, const MBool mayHaveChildren,
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver,
    const bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()> refineFlag) {
  static constexpr MInt noInternalConnections = (nDim == 2) ? 4 : 12;
  static constexpr MInt connectionDirs[12] = {1, 1, 3, 3, 1, 1, 3, 3, 5, 5, 5, 5};
  static constexpr MInt childs0[12] = {0, 2, 0, 1, 4, 6, 4, 5, 0, 1, 2, 3};
  static constexpr MInt childs1[12] = {1, 3, 2, 3, 5, 7, 6, 7, 4, 5, 6, 7};
  static constexpr MInt dirStencil[3][8] = {
      {0, 1, 0, 1, 0, 1, 0, 1}, {2, 2, 3, 3, 2, 2, 3, 3}, {4, 4, 4, 4, 5, 5, 5, 5}};
  static constexpr MInt sideIds[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {0, 0, 1, 1, 0, 0, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}};
  static constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};
  static constexpr MInt otherSide[2] = {1, 0};
  static constexpr MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                               {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};

  ASSERT(m_maxRfnmntLvl >= m_maxLevel, "");
  if(a_level(cellId) >= m_maxRfnmntLvl) return;

  const MInt childLevel = a_level(cellId) + 1;
  const MFloat childCellLength = cellLengthAtLevel(childLevel);

  MInt noOutsideChilds = 0;
  for(MInt c = 0; c < m_maxNoChilds; c++) {
    MInt isOutside[maia::grid::tree::Tree<nDim>::maxNoSolvers()];
    std::fill_n(isOutside, maia::grid::tree::Tree<nDim>::maxNoSolvers(), -1);

    if(!mayHaveChildren && a_childId(cellId, c) > -1) {
      // mTerm(1, AT_, "Not supposed to happen." );
    }

    if(a_childId(cellId, c) > -1) continue;

    if(refineChildIds != nullptr && refineChildIds[c] < 0) {
      a_childId(cellId, c) = -1;
      continue;
    }

    MFloat coords[nDim];
    for(MInt i = 0; i < nDim; i++) {
      coords[i] = a_coordinate(cellId, i) + F1B2 * signStencil[c][i] * childCellLength;
    }
    MInt allOutside = 0;
    if(!cellOutsideSolver.empty()) {
      allOutside = 1;
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(refineFlag[solver]) {
          isOutside[solver] = cellOutsideSolver[solver](coords, childLevel, cellId);
        } else {
          isOutside[solver] = 1;
        }
        allOutside = mMin(allOutside, isOutside[solver]); // call cellOutsideSolver() function in each solver
      }
    }
    if(allOutside > 0) {
      noOutsideChilds++;
      a_childId(cellId, c) = -1;
      continue;
    }

    MInt childId;
    if(m_freeIndices.size() > 0) {
      auto it = m_freeIndices.begin();
      childId = *(it);
      m_freeIndices.erase(it);
      ASSERT(childId > -1 && childId < m_tree.size(), "");
    } else {
      childId = m_tree.size();
      m_tree.append();
    }

    for(MInt i = 0; i < nDim; i++) {
      a_coordinate(childId, i) = coords[i];
    }

    if(refineChildIds != nullptr) {
      a_globalId(childId) = refineChildIds[c];
    } else {
      a_globalId(childId) = -1;
    }

    treeb().resetSolver(childId);
    MBool any = false;
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(refineFlag[solver] && isOutside[solver] < 1) treeb().solver(childId, solver) = true;
      any = any || m_tree.solver(childId, solver);
    }
    ASSERT(any, "");

    a_level(childId) = childLevel;

    a_childId(cellId, c) = childId;
    a_parentId(childId) = cellId;
    a_isToDelete(childId) = false;
    a_resetProperties(childId);
    a_hasProperty(childId, Cell::IsHalo) = (refineChildIds != nullptr);
    // a_hasProperty( cellId, Cell::IsHalo);
    //    a_hasProperty( childId, Cell::IsPeriodic) = a_hasProperty( cellId, Cell::IsPeriodic);

    for(MInt k = 0; k < m_maxNoChilds; k++)
      a_childId(childId, k) = -1;
    for(MInt d = 0; d < 2 * nDim; d++) {
      a_neighborId(childId, d) = -1;
    }

    if(treeb().noSolvers() == 1) {
      treeb().solver(childId, 0) = true; // remove if multisolver info set correctly!!!
    }

    treeb().resetIsLeafCell(childId);

    MInt nodeId[3];
    MInt revNodeId[3];
    for(MInt i = 0; i < nDim; i++) {
      nodeId[i] = sideIds[i][c] * IPOW2(i);
      revNodeId[i] = otherSide[sideIds[i][c]] * IPOW2(i);
    }
    for(MInt d = 0; d < nDim; d++) {
      const MInt dir = dirStencil[d][c];
      if(a_hasNeighbor(cellId, dir) == 0) {
        a_neighborId(childId, dir) = -1;
      } else if(a_noChildren(a_neighborId(cellId, dir)) == 0) {
        a_neighborId(childId, dir) = -1;
      } else {
        const MInt childNode = c - nodeId[d] + revNodeId[d];
        const MInt nghbrId = a_childId(a_neighborId(cellId, dir), childNode);
        if(nghbrId < 0) {
          a_neighborId(childId, dir) = -1;
          continue;
        }
        a_neighborId(childId, dir) = nghbrId;
        a_neighborId(nghbrId, revDir[dir]) = childId;
      }
    }
    a_hasProperty(childId, Cell::WasNewlyCreated) = true;

    // Set a default weight and the number of offspring
    a_weight(childId) = 1.0;
    a_noOffsprings(childId) = 1;
  }

  if(noOutsideChilds == m_maxNoChilds) {
    // all childs are outside!
    // NOTE: this occours when there is a difference in the geometry used in the grid-generation
    //      geometryroot.cpp-isPointInsideNode and
    //      geometry.cpp-pointIsInside of the geometry during the solver run!
    // FIX: call again but this time without the solver-Outside check,
    //     meaning that all cells will be refined instead!
    // NOTE: How is that a fix?
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>> cellOutsideEmpty;
    refineCell(cellId, nullptr, false, cellOutsideEmpty, refineFlag);
    return;
  }


  ASSERT(noOutsideChilds < m_maxNoChilds,
         to_string(noOutsideChilds) + "Check your geometry: all childs are outside or should not be refined! ");
  ASSERT(a_noChildren(cellId) > 0, "");
  treeb().resetIsLeafCell(cellId);

  for(MInt c = 0; c < noInternalConnections; c++) {
    const MInt dir = connectionDirs[c];
    const MInt child0 = a_childId(cellId, childs0[c]);
    const MInt child1 = a_childId(cellId, childs1[c]);
    if(child0 > -1 && child1 > -1) {
      a_neighborId(child0, dir) = child1;
      a_neighborId(child1, revDir[dir]) = child0;
    }
  }

  a_hasProperty(cellId, Cell::WasRefined) = true;
}


// --------------------------------------------------------------------------------------


/**
 * \brief removes the children of the given cell
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::removeChilds(const MInt cellId) {
  TRACE();

  if(a_noChildren(cellId) <= 0) {
    mTerm(1, AT_, "Unexpected situation 1 in CartesianGrid::removeChilds(MInt childId)");
  }
  for(MInt i = 0; i < m_maxNoChilds; i++) {
    const MInt childId = a_childId(cellId, i);

    if(childId > -1) {
      removeCell<true>(childId);
      a_childId(cellId, i) = -1;
    }

    treeb().resetIsLeafCell(cellId);
  }

  a_hasProperty(cellId, Cell::WasCoarsened) = true;
}


/**
 * \brief removes the children of the given cell
 * \author Lennart Schneiders
 */
template <MInt nDim>
template <MBool removeChilds_>
void CartesianGrid<nDim>::removeCell(const MInt cellId) {
  static constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};

  ASSERT(cellId > -1 && cellId < treeb().size(), "");
  ASSERT(a_noChildren(cellId) == 0, "No children expected for child. ");
  ASSERT(!a_hasProperty(cellId, Cell::IsPartitionCell) || a_hasProperty(cellId, Cell::IsHalo),
         "Error: cannot remove non-halo partition cell");

  for(MInt dir = 0; dir < m_noDirs; dir++) {
    if(a_hasNeighbor(cellId, dir) > 0) {
      MInt nghbrId = a_neighborId(cellId, dir);
      if((nghbrId < 0) || (nghbrId >= treeb().size())) {
        mTerm(1, AT_,
              "Unexpected situation 2 in CartesianGrid::removeChildIds() >> " + to_string(cellId) + " " + to_string(dir)
                  + " " + to_string(nghbrId) + " " + to_string(treeb().size()));
      }
      a_neighborId(nghbrId, revDir[dir]) = -1;
      a_neighborId(cellId, dir) = -1;
    } else {
      a_neighborId(cellId, dir) = -1;
    }
  }
  const MInt parentId = a_parentId(cellId);
  a_parentId(cellId) = -1;
  a_globalId(cellId) = -1;
  treeb().resetSolver(cellId);

  a_resetProperties(cellId);

  a_isToDelete(cellId) = true;

  a_level(cellId) = -1;
  treeb().resetIsLeafCell(cellId);

  if(cellId == (treeb().size() - 1)) {
    treeb().size(treeb().size() - 1);
  } else {
    m_freeIndices.insert(cellId);
  }

  // return early if called from inside removeChilds
  if(removeChilds_ || parentId == -1) return;

  ASSERT(parentId < treeb().size(), "");

  //
  for(MInt i = 0; i < m_maxNoChilds; i++) {
    if(a_childId(parentId, i) == cellId) {
      a_childId(parentId, i) = -1;
      break;
    }
  }

  //
  treeb().resetIsLeafCell(parentId);
  for(MInt solverId = 0; solverId < m_tree.noSolvers(); solverId++) {
    if(m_tree.solver(parentId, solverId)) {
      a_isLeafCell(parentId, solverId) =
          !a_hasChildren(parentId, solverId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor);
    }
  }
}

//----------------------------------------------------------------------------------


/** \brief Necessary for extractPointIdsFromGrid function
 *
 * This function creates all possible paths from a cell to its
 * adjacent neighbors (in a uniform grid) and saves them.
 * Of course the possible paths could also be hardcoded but I'm
 * too lazy, so feel free to take the output of this function
 * and to hardcode it into a fixed array.
 *
 */
template <MInt nDim>
void CartesianGrid<nDim>::createPaths() {
  TRACE();

  MInt offset[6] = {2, 2, 4, 4, 6, 6};
  MInt offset2[4] = {2, 2, 0, 0};
  MInt allNeighbors[(nDim == 2) ? 8 : 12];
  MInt secondLevelNeighbors[8];

  IF_CONSTEXPR(nDim == 2) {
    MInt allNeighbors2D[8] = {0, 1, 2, 3, 0, 1, 2, 3};
    for(MInt i = 0; i < 8; i++) {
      allNeighbors[i] = allNeighbors2D[i];
    }
    // Holds the corresponding neighbor directions for the binary codes according to D2Q9.
    // 9 means in 2D that the corresponding code doesn't exist.
    MInt neighborCode2D[11] = {9, 0, 1, 9, 2, 6, 5, 9, 3, 7, 4};
    for(MInt i = 0; i < 11; i++) {
      m_neighborCode[i] = neighborCode2D[i];
    }
  }
  else {
    MInt allNeighbors3D[12] = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
    for(MInt i = 0; i < 12; i++) {
      allNeighbors[i] = allNeighbors3D[i];
    }
    // Holds the corresponding neighbor directions for the binary codes according to D3Q27.
    // 27 means in 3D that the corresponding code doesn't exist.
    MInt neighborCode3D[43] = {27, 0,  1,  27, 2,  6,  8,  27, 3,  7,  9, 27, 27, 27, 27, 27, 4,  10, 12, 27, 14, 18,
                               22, 27, 16, 20, 24, 27, 27, 27, 27, 27, 5, 11, 13, 27, 15, 19, 23, 27, 17, 21, 25};
    for(MInt i = 0; i < 43; i++) {
      m_neighborCode[i] = neighborCode3D[i];
    }
  }
  m_counter1D = 0;
  m_counter2D = 0;
  m_counter3D = 0;

  // Determine counters
  // All 1D neighbors
  for(MInt i = 0; i < nDim * 2; i++) {
    m_counter1D++;
    // Remaining possible neighbors
    for(MInt j = 0; j < (nDim * 2 - 2); j++) {
      m_counter2D++;
      // Remaining possible neighbors
      for(MInt k = 0; k < nDim * 2 - 2 - 2; k++) {
        m_counter3D++;
      }
    }
  }

  // Allocate memory
  mDeallocate(m_paths1D);
  mAlloc(m_paths1D, m_counter1D, "m_paths1D", AT_);
  mDeallocate(m_paths2D);
  mAlloc(m_paths2D, m_counter2D, 2, "m_paths2D", AT_);
  IF_CONSTEXPR(nDim == 3) {
    mDeallocate(m_paths3D);
    mAlloc(m_paths3D, m_counter3D, 3, "m_paths3D", AT_);
  }

  m_counter1D = 0;
  m_counter2D = 0;
  m_counter3D = 0;

  // save the paths...
  for(MInt i = 0; i < nDim * 2; i++) {
    // Save 1D paths
    m_paths1D[m_counter1D] = allNeighbors[i];
    m_counter1D++;
    // Remaining possible neighbors
    for(MInt j = 0; j < (nDim * 2 - 2); j++) {
      // Fill secondeLevelNeighbors twice!
      secondLevelNeighbors[j] = allNeighbors[offset[i] + j];
      secondLevelNeighbors[j + nDim * 2 - 2] = allNeighbors[offset[i] + j];
    }
    for(MInt j = 0; j < (nDim * 2 - 2); j++) {
      // Save 2D paths
      m_paths2D[m_counter2D][0] = allNeighbors[i];
      m_paths2D[m_counter2D][1] = secondLevelNeighbors[j];
      m_counter2D++;
      // Remaining possible neighbors
      IF_CONSTEXPR(nDim == 3) {
        MInt thirdLevelNeighbors[2];
        for(MInt k = 0; k < nDim * 2 - 2 - 2; k++) {
          // Fill thirdLevelNeighbors
          thirdLevelNeighbors[k] = secondLevelNeighbors[offset2[j] + k];
          // Save 3D paths
          m_paths3D[m_counter3D][0] = allNeighbors[i];
          m_paths3D[m_counter3D][1] = secondLevelNeighbors[j];
          m_paths3D[m_counter3D][2] = thirdLevelNeighbors[k];

          m_counter3D++;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

/**
 * Wrapper function to choose between two different adaptation modes based on m_lowMemAdaptation
 * True: uses a new adaptation procedure in which information of partition level ancestors is exchanged
 * serially for every neighbor domain. This avoids vast allocations of memory but is potentially slower
 * False: uses the old adpatation procedure in which information of partition level ancestors in exchanged
 * for all domains at once. This can require large amounts of memory for high numbers of neighbor domains
 **/
template <MInt nDim>
void CartesianGrid<nDim>::meshAdaptation(
    vector<vector<MFloat>>& sensors, vector<MFloat>& sensorWeight, vector<bitset<64>>& sensorCellFlag,
    vector<MInt>& sensorSolverId, const std::vector<std::function<void(const MInt)>>& refineCellSolver,
    const std::vector<std::function<void(const MInt)>>& removeChildsSolver,
    const std::vector<std::function<void(const MInt)>>& removeCellSolver,
    const std::vector<std::function<void(const MInt, const MInt)>>& swapCellsSolver,
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver,
    const std::vector<std::function<void()>>& resizeGridMapSolver) {
  TRACE();

  if(m_lowMemAdaptation) {
    this->meshAdaptationLowMem(sensors, sensorWeight, sensorCellFlag, sensorSolverId, refineCellSolver,
                               removeChildsSolver, removeCellSolver, swapCellsSolver, cellOutsideSolver,
                               resizeGridMapSolver);
  } else {
    this->meshAdaptationDefault(sensors, sensorWeight, sensorCellFlag, sensorSolverId, refineCellSolver,
                                removeChildsSolver, removeCellSolver, swapCellsSolver, cellOutsideSolver,
                                resizeGridMapSolver);
  }
}


/**
 * \brief Performs mesh adaptation and continuously updates the window/halo cells
 * \author Lennart Schneiders
 *
 * \param[in] sensors Holds the mesh refinement sensor(s) for each cell
 * \param[in] sensorWeight Weighting of this sensor in relation to other sensors if greater than zero, otherwise
 * indicating a true/false-type sensor \param[in] sensorCellFlag Indicator whether a sensor was set for this particular
 * cell \param[out] oldCellId Mapping from new grid cells to their original cellId, or -1 if they did not exist
 * previously \param[out] oldNoCells Number of cells in the map
 */
template <MInt nDim>
void CartesianGrid<nDim>::meshAdaptationLowMem(
    vector<vector<MFloat>>& sensors, vector<MFloat>& sensorWeight, vector<bitset<64>>& sensorCellFlag,
    vector<MInt>& sensorSolverId, const std::vector<std::function<void(const MInt)>>& refineCellSolver,
    const std::vector<std::function<void(const MInt)>>& removeChildsSolver,
    const std::vector<std::function<void(const MInt)>>& removeCellSolver,
    const std::vector<std::function<void(const MInt, const MInt)>>& swapCellsSolver,
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver,
    const std::vector<std::function<void()>>& resizeGridMapSolver) {
  TRACE();

  //----------------------------
  // 0. init
  const MInt noCells = treeb().size();
  const MInt maxNoCells = m_maxNoCells;
  const MInt maxLevel = mMin(m_maxRfnmntLvl, m_maxLevel + 1);

  const MInt noSensors = (signed)sensorWeight.size();
  ASSERT(sensorWeight.size() == sensors.size(), "");
  ASSERT(sensorWeight.size() == sensorSolverId.size(), "");
  ScratchSpace<MFloat> sensorCnt(noSensors, 2, AT_, "sensorCnt");
  ScratchSpace<MFloat> sensorThresholds(noSensors, 2, AT_, "sensorThresholds");
  ScratchSpace<bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()>> refineFlag(maxNoCells, AT_, "refineFlag");
  ScratchSpace<bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()>> coarseFlag(maxNoCells, AT_, "coarseFlag");
  ASSERT(m_freeIndices.empty(), "");
  m_freeIndices.clear();
  sensorCnt.fill(F0);
  for(MInt c = 0; c < maxNoCells; c++) {
    refineFlag[c].reset();
    coarseFlag[c].reset();
  }
  // these flags indicate which cells have changed for the subsequent reinitialization by the solvers
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    a_hasProperty(cellId, Cell::IsToDelete) = false;
    a_hasProperty(cellId, Cell::WasNewlyCreated) = false;
    a_hasProperty(cellId, Cell::WasCoarsened) = false;
    a_hasProperty(cellId, Cell::WasRefined) = false;

    // Reset the window cell flag,
    // is set for all preserved window cells when added to m_windowCells
    a_hasProperty(cellId, Cell::IsWindow) = false;
  }

  //----------------------------
  // 1. exchange sensors and compute RMS values
  // TODO labels:GRID move these parts to the cartesian solver and do this sort of smoothing
  //       on the solver grid in setSensors, after all sensors have been set there!!!!
  //       the cartesiangrid should only get one refine/coarse flags for each cell from the solver!
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    if(a_isToDelete(cellId)) continue;

    // Note: sensors for partition level ancestors might not be set correctly!
    if(a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
      continue;
    }

    for(MInt s = 0; s < noSensors; s++) {
      if(sensorCellFlag[cellId][s]) {
        ASSERT(sensorWeight[s] < F0 || std::fabs(sensors[s][cellId]) < MFloatMaxSqrt,
               "invalid sensor value: " + std::to_string(sensors[s][cellId]));
        sensorCnt(s, 0) += POW2(sensors[s][cellId]);
        sensorCnt(s, 1) += F1;
      }
    }
  }

  if(noSensors > 0) {
    MPI_Allreduce(MPI_IN_PLACE, &sensorCnt(0), 2 * noSensors, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "sensorCnt(0)");
  }

  // TODO labels:GRID,ADAPTATION introduce sensor type to replace this weight hack
  for(MInt s = 0; s < noSensors; s++) {
    if(sensorWeight[s] < F0) { // true/false sensor
      sensorThresholds(s, 0) = -1e-12;
      sensorThresholds(s, 1) = 1e-12;
    } else { // continuous sensor
      sensorThresholds(s, 1) = sensorWeight[s] * sqrt(mMax(1e-14, (sensorCnt(s, 0) / sensorCnt(s, 1))));
      sensorThresholds(s, 0) = m_coarseRatio * sensorThresholds(s, 1);
    }
  }

  //---------------------------(1)

  const MBool hasPartLvlShift = (m_maxPartitionLevelShift > 0);
  const MInt tmpScratchSize = (hasPartLvlShift) ? m_tree.size() : 1;

  ScratchSpace<MBool> isHaloPartLvlAncestor(tmpScratchSize, AT_, "isHaloPartLvlAncestor");

  isHaloPartLvlAncestor.fill(false);

  // Determine relevant halo partition level ancestor window/halo cells that need to be preserved
  std::vector<std::vector<MInt>> oldWindowCells;
  std::vector<MInt> oldNoWindowCells;
  MInt totalNoWindowCells = 0;

  for(MInt d = 0; d < noNeighborDomains(); d++) {
    oldWindowCells.push_back(m_windowCells[d]);
    oldNoWindowCells.push_back(noWindowCells(d));
    totalNoWindowCells += noWindowCells(d);
  }

  MInt recvSize = (m_maxPartitionLevelShift > 0) ? totalNoWindowCells : 1;
  ScratchSpace<MBool> recvIsHaloPartLvlAncestor(recvSize, AT_, "recvIsHaloPartLvlAncestor");

  if(m_maxPartitionLevelShift > 0) {
    for(auto& i : m_partitionLevelAncestorIds) {
      TERMM_IF_NOT_COND(a_hasProperty(i, Cell::IsPartLvlAncestor),
                        "cell is not a partition level ancestor: " + std::to_string(i));
      // Mark halo partition level ancestors (which are part of the missing subtree of the grid)
      if(a_hasProperty(i, Cell::IsHalo)) {
        isHaloPartLvlAncestor[i] = true;
      }
      // Mark all halo childs of partition level ancestors (required for exchange of e.g. global
      // ids!)
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(i, child);
        if(childId > -1 && a_isHalo(childId)) {
          isHaloPartLvlAncestor[childId] = true;
        }
      }
    }

    // Reverse exchange markers from halo cells to corresponding window cells
    if(noNeighborDomains() > 0) {
#ifndef PVPLUGIN
      maia::mpi::reverseExchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &isHaloPartLvlAncestor[0],
                                     &recvIsHaloPartLvlAncestor[0]);
#else
      TERMM(1, "not working for compilation of paraview plugin");
#endif
    }
  }

  //----------------------------
  // 2. tag cells for coarsening / refinement

  // init parent cells with coarse flag true and set false if any child intervenes
  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "");
    ASSERT(!a_isToDelete(cellId), "");
    if(a_level(cellId) > m_maxUniformRefinementLevel) {
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(a_hasChildren(cellId, solver)) continue;

        const MInt parentId = a_parentId(cellId, solver);
        ASSERT(parentId > -1, "");
        // Partition level shift, parent is a halo cell and cannot be coarsened
        if(a_isHalo(parentId) || a_hasProperty(parentId, Cell::IsPartLvlAncestor)) continue;

        for(MInt s = 0; s < noSensors; s++) {
          if(sensorSolverId[s] != solver) continue;
          if(sensorCellFlag[cellId][s] || sensorCellFlag[parentId][s]) {
            // only coarsen cells if they are above the newMinLevel!
            if(a_level(cellId) > m_newMinLevel && m_allowCoarsening) {
              coarseFlag[a_parentId(cellId)][solver] = true;
            }
          }
        }
      }
    }
  }

  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "");
    ASSERT(!a_isToDelete(cellId), "");
    const MInt level = a_level(cellId);
    //---
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_tree.solver(cellId, solver)) continue;
      if(a_hasChildren(cellId, solver)) continue;
      const MInt parentId = a_parentId(cellId, solver);
      // Partition level shift, if parent is a halo it cannot be coarsened and e.g.
      // sensors/coarseFlag[parentId] should not be accessed!
      const MBool partLvlAncestorParent = (parentId > -1) ? a_hasProperty(parentId, Cell::IsPartLvlAncestor) : false;

      // cell should be marked for refinement
      MBool refine = false;
      // cell should be marked for coarsened
      MBool coarsen = true;
      // cell should be marked for refinement regardless
      // whether a different sensor marks this cell for coarsening
      MBool forceRefine = false;
      for(MInt s = 0; s < noSensors; s++) {
        if(sensorSolverId[s] != solver) continue;
        if(sensorCellFlag[cellId][s]) {
          coarsen = coarsen && (sensors[s][cellId] < sensorThresholds(s, 0));
        }
        if(level > m_maxUniformRefinementLevel && !partLvlAncestorParent && sensorCellFlag[parentId][s]) {
          coarsen = coarsen && (sensors[s][parentId] < sensorThresholds(s, 0));
        }
        // this sensor would refine the cell
        const MBool refineSensor = sensors[s][cellId] > sensorThresholds(s, 1) && sensorCellFlag[cellId][s];
        refine = refine || refineSensor;

        // a true/false sensor marks the cell for refinement
        if(refineSensor && sensorWeight[s] < F0) {
          forceRefine = true;
        }
      }
      if(level < m_maxRfnmntLvl) {
        refineFlag[cellId][solver] = refine;
      }

      if(level > m_maxUniformRefinementLevel && !partLvlAncestorParent) {
        coarseFlag[parentId][solver] = coarseFlag[parentId][solver] && coarsen;
      }
      // force a refinement of a cell:
      // for all minLevel cells which should be raised to the newMinLevel
      // for the cells forced by any true/false sensor!
      if((forceRefine && level < m_maxRfnmntLvl) || a_level(cellId) < m_newMinLevel) {
        refineFlag[cellId][solver] = true;
        coarseFlag[cellId][solver] = false;
        if(level > m_maxUniformRefinementLevel) {
          coarseFlag[parentId][solver] = false;
        }
      }
    }
  }

  // 3. Exchange since consistency checks need neighbor values
  if(noNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), coarseFlag.getPointer(),
                              m_tree.size());
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), refineFlag.getPointer(),
                              m_tree.size());
  }

  //----------------------------
  // 4. consistency check for internal cells
  if(m_checkRefinementHoles != nullptr) {
    // necessary for solution adaptive refinement (i.e. at shocks)
    const MInt adaptationHoleLimit = nDim;
    MBool gridUpdated = true;

    MInt loopCount = 0;
    while(gridUpdated && loopCount < 25) {
      MInt noChanges = 0;
      loopCount++;
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_checkRefinementHoles[solver]) continue;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          if(a_isToDelete(cellId)) continue;
          if(a_hasProperty(cellId, Cell::IsHalo)) continue;
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) continue;

          // check, if current cell is a hole in the refinement => refine cell aswell
          if(!refineFlag[cellId][solver]
             || (a_parentId(cellId, solver) > -1 && coarseFlag[a_parentId(cellId, solver)][solver])) {
            MInt checkRefinedNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                const MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(refineFlag[nghbrId][solver] || a_hasChildren(nghbrId, solver)) {
                  checkRefinedNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkRefinedNeighbors > adaptationHoleLimit) {
              refineFlag[cellId][solver] = true;
              coarseFlag[cellId][solver] = false;
              const MInt parentId = a_parentId(cellId, solver);
              if(parentId > -1) coarseFlag[parentId][solver] = false;
              noChanges += 1;
            }
          }

          // check, if current cell coarsening would create a hole in the refinement => do not coarsen cell
          const MInt parentId = a_parentId(cellId, solver);
          if(parentId < 0) continue;
          if(coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            MInt checkRefinedNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(parentId, dir, solver) > 0) {
                const MInt parNghbrId = a_neighborId(parentId, dir, solver);
                if(!coarseFlag[parNghbrId][solver]
                   && (refineFlag[parNghbrId][solver] || a_hasChildren(parNghbrId, solver))) {
                  checkRefinedNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkRefinedNeighbors > adaptationHoleLimit) {
              coarseFlag[parentId][solver] = false;
              noChanges += 1;
            }
          }

          // Check if the coarsening of neighboring cells would create a refined cell island => coarsen cell aswell
          if(!coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            MInt checkCoarseNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(parentId, dir, solver) > 0) {
                const MInt parNghbrId = a_neighborId(parentId, dir, solver);
                if(coarseFlag[parNghbrId][solver]
                   || (!refineFlag[parNghbrId][solver] && !a_hasChildren(parNghbrId, solver))) {
                  checkCoarseNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkCoarseNeighbors > currentHoleLimit) {
              coarseFlag[parentId][solver] = true;
              noChanges += 1;
            }
          }

          // check if the current cell refinement creates a refined cell island => do not refine cell
          if(refineFlag[cellId][solver]) {
            MInt checkRefinedNeighbors = 0;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                const MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(refineFlag[nghbrId][solver] || a_hasChildren(nghbrId, solver)) {
                  checkRefinedNeighbors += 1;
                }
              }
            }
            if(checkRefinedNeighbors == 0) {
              refineFlag[cellId][solver] = false;
              noChanges += 1;
            }
          }
        }
      }
      if(noChanges == 0) {
        gridUpdated = false;
      }
    }
  }

  // smooth diagonal level jumps,
  // even if already existing in the base grid!
  if(m_diagSmoothing != nullptr) {
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_diagSmoothing[solver]) continue;
      for(MInt lvl = maxLevel - 1; lvl >= minLevel(); lvl--) {
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          if(a_level(cellId) != lvl) continue;
          if(!a_isLeafCell(cellId)) continue;
          if(!m_tree.solver(cellId, solver)) continue;
          if(refineFlag[cellId][solver]) continue;
          for(MInt i = 0; i < m_noDirs; i++) {
            const MInt nghbrId = a_neighborId(cellId, i, solver);
            if(nghbrId < 0) continue;
            if(a_isLeafCell(nghbrId)) continue;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(nghbrId, child, solver);
              if(childId < 0) continue;
              if(!a_isLeafCell(childId)) {
                refineFlag[cellId][solver] = true;
                coarseFlag[cellId][solver] = false;
                MInt parentId = a_parentId(cellId, solver);
                if(parentId > -1) {
                  coarseFlag[cellId][solver] = false;
                }
                break;
              }
            }
            for(MInt j = 0; j < m_noDirs; j++) {
              if(j / 2 == i / 2) continue;
              const MInt diagNghbrId = a_neighborId(nghbrId, j, solver);
              if(diagNghbrId < 0) continue;
              if(a_isLeafCell(diagNghbrId)) continue;
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                MInt childId = a_childId(diagNghbrId, child, solver);
                if(childId < 0) continue;
                if(!a_isLeafCell(childId)) {
                  refineFlag[cellId][solver] = true;
                  coarseFlag[cellId][solver] = false;
                  MInt parentId = a_parentId(cellId, solver);
                  if(parentId > -1) {
                    coarseFlag[cellId][solver] = false;
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Additional consistency checks fixing contradictory flags
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    if(a_isToDelete(cellId)) continue;

    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_tree.solver(cellId, solver)) continue;
      if(a_hasChildren(cellId, solver)) continue;
      const MInt parentId = a_parentId(cellId, solver);
      if(refineFlag[cellId][solver] && coarseFlag[cellId][solver]) {
        refineFlag[cellId][solver] = false;
        coarseFlag[cellId][solver] = false;
      }

      // TODO labels:GRID,ADAPTATION,totest check for multisolver!!!!!
      if(refineFlag[cellId][solver]) {
        if(parentId > -1) {
          if(coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            coarseFlag[parentId][solver] = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(parentId, child, solver);
              if(childId > -1) {
                refineFlag[childId][solver] = false;
              }
            }
          }
        }
      }
    }
  }
  //---------------------------(4)

  // ensure the same refinement between two different solvers in the combined/overlapping region!
  if(m_noIdenticalSolvers > 0) {
    for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
      for(MInt pair = 0; pair < m_noIdenticalSolvers; pair++) {
        const MInt slaveId = m_identicalSolvers[pair * 2];
        const MInt masterId = m_identicalSolvers[pair * 2 + 1];
        if(treeb().solver(cellId, slaveId) && treeb().solver(cellId, masterId)) {
          if(!coarseFlag[cellId][masterId]) {
            coarseFlag[cellId][slaveId] = coarseFlag[cellId][masterId];
          }
          if(m_identicalSolverLvlJumps[pair] == 0) {
            refineFlag[cellId][slaveId] = refineFlag[cellId][masterId];
          } else if(m_identicalSolverLvlJumps[pair] == 1) {
            if(a_level(cellId) < m_identicalSolverMaxLvl[pair] && a_hasChildren(cellId, masterId)) {
              if(!a_hasChildren(cellId, slaveId)) {
                refineFlag[cellId][slaveId] = true;
              }
              coarseFlag[cellId][slaveId] = false;
            }
          } else {
            if(a_level(cellId) < m_identicalSolverMaxLvl[pair] && a_hasChildren(cellId, masterId)
               && a_childId(cellId, 0) > -1 && treeb().solver(a_childId(cellId, 0), masterId)
               && a_hasChildren(a_childId(cellId, 0), masterId)) {
              if(!a_hasChildren(cellId, slaveId)) {
                refineFlag[cellId][slaveId] = true;
              }
              coarseFlag[cellId][slaveId] = false;
            }
          }
        }
      }
    }
  }


  //----------------------------
  // 5. Coarsening step: loop from maxLevel to minLevel and exchange coarseFlags to retrieve
  //    consistent refinement at domain interfaces; update window/halo collectors
  for(MInt level = m_maxRfnmntLvl - 1; level >= m_maxUniformRefinementLevel; level--) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      // TODO labels:GRID,ADAPTATION,totest check for multisolver setting.. neighborIds solver aware!!
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(!a_hasChildren(cellId, solver)) continue;

        // Partition level ancestor cell cannot be coarsened (exept if its a partition cell)
        if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
          coarseFlag[cellId][solver] = false;
        }

        if(coarseFlag[cellId][solver]) {
          if(!coarsenCheck(cellId, solver)) {
            coarseFlag[cellId][solver] = false;
          }
        }
      }
    }

    if(noNeighborDomains() > 0) {
      //      exchangeSolverBitset(coarseFlag.data());
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), coarseFlag.data(),
                                m_tree.size());
    }

    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      // TODO labels:GRID,ADAPTATION,totest check for multisolver setting.. neighborIds solver aware!!
      if(level < m_maxRfnmntLvl - 1) {
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) continue;
          if(!refineFlag[cellId][solver] && !coarseFlag[cellId][solver]) {
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(m_azimuthalPer && a_hasProperty(nghbrId, Cell::IsPeriodic)) continue;
                if(a_hasChildren(nghbrId, solver) > 0) {
                  for(MInt child = 0; child < m_maxNoChilds; child++) {
                    if(!childCode[dir][child]) continue;
                    MInt childId = a_childId(nghbrId, child, solver);
                    if(childId < 0) continue;
                    if(refineFlag[childId][solver]) {
                      refineFlag[cellId][solver] = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if(noNeighborDomains() > 0) {
      //      exchangeSolverBitset(refineFlag.data());
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), refineFlag.data(),
                                m_tree.size());
    }

    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(m_azimuthalPer && a_hasProperty(cellId, Cell::IsPeriodic)) continue;

      // Partition level ancestor cell cannot be coarsened (exept if its a partition cell)
      if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
        continue;
      }

      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(coarseFlag[cellId][solver]) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt childId = a_childId(cellId, child, solver);
            if(childId > -1) {
              coarseFlag[childId][solver] = false;
              refineFlag[childId][solver] = false;
            }
          }
          if(a_hasChildren(cellId, solver)) {
            removeChildsSolver[solver](cellId);
            // call removeChilds() function in each solver before child links are deleted
            refineFlag[cellId][solver] = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId > -1) {
                treeb().solver(childId, solver) = false;
              }
            }
          }
        }
      }
      // check if the child can also be removed from the cartesian/combined grid
      // if the child is no-longer in any of the solvers, the grid child will be removed!
      if(a_noChildren(cellId) > 0) {
        MBool rmChilds = true;
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(a_hasChildren(cellId, solver)) {
            rmChilds = false;
          }
        }
        if(rmChilds) {
          removeChilds(cellId);
        }
      }
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      std::set<MInt> deleted;
      // WH_old
      if(m_haloMode > 0) {
        for(auto it_ = m_windowLayer_[i].begin(); it_ != m_windowLayer_[i].end();) {
          if(a_isToDelete(it_->first)) {
            it_ = m_windowLayer_[i].erase(it_);
          } else {
            // We might have deleted a cell from a solver
            for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
              if(it_->second.get(solver) < WINDOWLAYER_MAX /*<=m_noSolverHaloLayers[solver]*/
                 && !m_tree.solver(it_->first, solver))
                it_->second.set(solver, WINDOWLAYER_MAX);
              if(it_->second.get(solver) > m_noSolverHaloLayers[solver]) it_->second.set(solver, WINDOWLAYER_MAX);
            }
            if(it_->second.all()) {
              deleted.insert(it_->first);
              it_ = m_windowLayer_[i].erase(it_);
            } else
              ++it_;
          }
        }
      }

      // 1. coarsen window cells
#ifndef NDEBUG
      std::set<MInt> windCellSet;
      std::set<MInt> haloCellSet;
#endif

      m_windowCells[i].clear();

      auto it = oldWindowCells[i].begin();

      for(it = oldWindowCells[i].begin(); it < oldWindowCells[i].end(); it++) {
        const MInt cellId = *it;
        if(!a_isToDelete(cellId) && (deleted.find(cellId) == deleted.end() || m_haloMode == 0)) { // WH_old
          ASSERT(cellId < treeb().size(), "");

#ifndef NDEBUG
          // TODO labels:GRID,ADAPTATION fix duplicate window cells in periodic cases!
          ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet.find(cellId) == windCellSet.end(),
                 "duplicate window cell: " + std::to_string(cellId));
          windCellSet.insert(cellId);
#endif

          m_windowCells[i].push_back(cellId);
        }
      }

      // 2. coarsen halo cells
      vector<MInt> oldHaloVec(m_haloCells[i]);
      m_haloCells[i].clear();
      for(it = oldHaloVec.begin(); it < oldHaloVec.end(); it++) {
        const MInt cellId = *it;
        if(!a_isToDelete(cellId)) {
          ASSERT(cellId < treeb().size(), "");

#ifndef NDEBUG
          ASSERT(haloCellSet.find(cellId) == haloCellSet.end(), "duplicate halo cell: " + std::to_string(cellId));
          haloCellSet.insert(cellId);
#endif

          if(m_tree.solverBits(cellId).any() || m_haloMode == 0) m_haloCells[i].push_back(cellId);
        }
      }
      oldHaloVec.clear();
    }
  }
  //---------------------------(5)


  //----------------------------
  // 6. Remove all halo cells, for levels > maxUniformRefinementLevel
  //    since the number of halo layers is inconsistent at this point
  //    Note: preserve halo partition level ancestors of the missing grid subtree in case of a
  //    partition level shift

  {
    for(MInt level = maxLevel - 1; level >= m_maxUniformRefinementLevel; level--) {
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        if(a_isToDelete(cellId)) continue;
        if(!a_hasProperty(cellId, Cell::IsHalo)) continue;
        if(a_level(cellId) == level) {
          if(a_noChildren(cellId) > 0) {
            // Preserve partition level ancestor halos (and their childs) if they are part of the
            // missing subtree of the grid in case of a partition level shift
            if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
              if(std::find(m_partitionLevelAncestorIds.begin(), m_partitionLevelAncestorIds.end(), cellId)
                 != m_partitionLevelAncestorIds.end()) {
                continue;
              }
            }

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;
              refineFlag[childId].reset();
              coarseFlag[childId].reset();
            }
            for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
              if(!m_tree.solver(cellId, solver)) continue;
              if(a_hasChildren(cellId, solver)) {
                removeChildsSolver[solver](cellId);
              }
            }
            removeChilds(cellId);
            coarseFlag[cellId].reset();
            refineFlag[cellId].reset();
          }
        }
      }
    }

    ScratchSpace<MBool> isWindowPartLvlAncestor_(tmpScratchSize, AT_, "isWindowPartLvLAncestor");
    MInt cnt = 0;
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      isWindowPartLvlAncestor_.fill(false);
      // now take data from recvIsHaloPartLvlAncestor and store it for this neighborDomain
      if(hasPartLvlShift) {
        for(MInt j = 0; j < oldNoWindowCells[i]; j++) {
          const MInt cellId = oldWindowCells[i][j];
          isWindowPartLvlAncestor_(cellId) = recvIsHaloPartLvlAncestor[cnt];
          cnt++;
        }
      }

      m_windowCells[i].clear();

      auto it = oldWindowCells[i].begin();
      for(it = oldWindowCells[i].begin(); it < oldWindowCells[i].end(); it++) {
        const MInt cellId = *it;

        // Add partition level ancestor window cells if required
        const MBool addWindow = (hasPartLvlShift) ? isWindowPartLvlAncestor_(cellId) : false;

        if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addWindow) {
          m_windowCells[i].push_back(cellId);
          // Set window cell flag for preserved window cells
          a_hasProperty(cellId, Cell::IsWindow) = true;
        }
      }

      // WH_old
      if(m_haloMode > 0) {
        for(auto it_ = m_windowLayer_[i].begin(); it_ != m_windowLayer_[i].end();) {
          const MInt cellId = it_->first;

          // Add partition level ancestor window cells if required
          const MBool addWindow = (hasPartLvlShift) ? isWindowPartLvlAncestor_(cellId) : false;

          if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addWindow)
            ++it_;
          else
            it_ = m_windowLayer_[i].erase(it_);
        }
      }

      vector<MInt> oldHaloVec(m_haloCells[i]);
      m_haloCells[i].clear();
      for(it = oldHaloVec.begin(); it < oldHaloVec.end(); it++) {
        const MInt cellId = *it;

        // Add partition level ancestor halo cells if required
        const MBool addHalo = (hasPartLvlShift) ? isHaloPartLvlAncestor[cellId] : false;

        if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addHalo) {
          ASSERT(cellId > -1 && cellId < treeb().size(), to_string(cellId) + " " + to_string(treeb().size()));
          m_haloCells[i].push_back(cellId);
        }
      }
      oldHaloVec.clear();
    }

    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldAzimuthalWindowVec(m_azimuthalWindowCells[i]);

        m_azimuthalWindowCells[i].clear();
        auto it = oldAzimuthalWindowVec.begin();
        for(MInt j = 0; j < (signed)oldAzimuthalWindowVec.size(); j++) {
          // Remove higher level azimuthal window cells from collector
          // If an azimuthal halo cell on a higher level is mapped to an
          // azimuthal window cell on a lower level. This azimuthal window cell
          // musst also be removed from the collector
          const MInt cellId = oldAzimuthalWindowVec[j];
          MBool keepWindow = true;
          if(find(m_azimuthalHigherLevelConnectivity[i].begin(), m_azimuthalHigherLevelConnectivity[i].end(), j)
             != m_azimuthalHigherLevelConnectivity[i].end()) {
            keepWindow = false;
          }
          if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel && keepWindow) {
            m_azimuthalWindowCells[i].push_back(cellId);
            a_hasProperty(cellId, Cell::IsWindow) = true;
          }
        }
        oldAzimuthalWindowVec.clear();

        vector<MInt> oldAzimuthalHaloVec(m_azimuthalHaloCells[i]);
        m_azimuthalHaloCells[i].clear();
        for(it = oldAzimuthalHaloVec.begin(); it < oldAzimuthalHaloVec.end(); it++) {
          const MInt cellId = *it;
          if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
            ASSERT(cellId > -1 && cellId < treeb().size(), "");
            m_azimuthalHaloCells[i].push_back(cellId);
          }
        }
        oldAzimuthalHaloVec.clear();
      }

      // Remove higher level unmapped azimuthal halo cells
      vector<MInt> oldUnmappedHaloVec(m_azimuthalUnmappedHaloCells);
      vector<MInt> oldUnmappedHaloDomVec(m_azimuthalUnmappedHaloDomains);
      m_azimuthalUnmappedHaloCells.clear();
      m_azimuthalUnmappedHaloDomains.clear();
      for(MUint i = 0; i < oldUnmappedHaloVec.size(); i++) {
        MInt cellId = oldUnmappedHaloVec[i];
        if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
          ASSERT(cellId > -1 && cellId < treeb().size(), "");
          m_azimuthalUnmappedHaloCells.push_back(cellId);
          m_azimuthalUnmappedHaloDomains.push_back(oldUnmappedHaloDomVec[i]);
        }
      }
      oldUnmappedHaloVec.clear();
      oldUnmappedHaloDomVec.clear();

      // Update gridBndryCells
      vector<MInt> oldGridBndryVec(m_gridBndryCells);
      m_gridBndryCells.clear();
      for(MUint i = 0; i < oldGridBndryVec.size(); i++) {
        MInt cellId = oldGridBndryVec[i];
        if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
          m_gridBndryCells.push_back(cellId);
        }
      }
      oldGridBndryVec.clear();
    }
  }

  //---------------------------(6)
  //----------------------------
  // 7. Refinement step: loop from minLevel to maxLevel and exchange refineFlags to retrieve
  //    a consistent refinement at domain interfaces; update window/halo connectivity at each level
  for(MInt level = m_maxUniformRefinementLevel; level < maxLevel; level++) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(a_hasChildren(cellId, solver)) continue;
        if(refineFlag[cellId][solver]) {
          if(!refineCheck(cellId, solver, cellOutsideSolver)) {
            // comment out for stl interface refinement, causes weird remaining cells
            refineFlag[cellId][solver] = 0;
          }
        }
      }

      if(refineFlag[cellId].any()) {
        refineCell(cellId, nullptr, false, cellOutsideSolver, refineFlag[cellId]);

        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!a_hasChildren(cellId, solver) && false) { // sohel false
            // this is somehow destroying multisolver cases
            // if one solver already has a cell the other solver doesnt get it even if it needs it..
            // i think the solver flag is not set correctly in the refineCell function of cartesiangrid
            refineFlag[cellId][solver] = false;
          }
          if(refineFlag[cellId][solver]) {
            // Set solver flag for child cells
            setChildSolverFlag(cellId, solver, cellOutsideSolver[solver]);
            refineCellSolver[solver](cellId); // call refineCell() function in each solver
          }
        }
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          MInt childId = a_childId(cellId, child);
          if(childId > -1 && a_hasProperty(childId, Cell::WasNewlyCreated)) {
            // Timw: only reset the child-flags if the cell was neawly created
            // otherwise this stops a refinement of the childs of a different solver!
            coarseFlag[childId].reset();
            refineFlag[childId].reset();
          }
        }
      }
    }

    setLevel();

    computeGlobalIds();

    // WH_old
    if(m_haloMode > 0)
      // noOffspring required in createHigherLevelExchangeCells
      calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), m_tree.size());

    // Update child ids of partition level ancestors (stored as global ids)
    m_partitionLevelAncestorChildIds.clear();
    const MInt noPartLvlAncestorIds = m_partitionLevelAncestorIds.size();
    for(MInt i = 0; i < noPartLvlAncestorIds; i++) {
      const MInt cellId = m_partitionLevelAncestorIds[i];

      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(cellId, child);
        const MLong childGlobalId = (childId > -1) ? a_globalId(childId) : -1;
        m_partitionLevelAncestorChildIds.push_back(childGlobalId);
      }
    }

    if(noDomains() > 1 || m_noPeriodicCartesianDirs > 0) {
      updateHaloCellCollectors();

      // WH_old
      if(m_haloMode > 0) {
        createHigherLevelExchangeCells(level, true, refineCellSolver, removeCellSolver, level == maxLevel - 1);
        // createHigherLevelExchangeCells(level, swapCellsSolver, refineCellSolver);
      } else
        createHigherLevelExchangeCells_old(level, refineCellSolver);
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        MInt cellId = m_haloCells[i][j];
        if(a_level(cellId) != level) continue;
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          const MInt childId = a_childId(cellId, child);
          if(childId < 0) continue;

          // Partition level shift: do not reset flags if the child of a halo cell is an internal
          // cell and already existed!
          if(!a_isHalo(childId)) {
            ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                   "child is not a partition level ancestor or a partition cell");
            a_hasProperty(cellId, Cell::WasRefined) = true;
          } else {
            coarseFlag[childId].reset();
            refineFlag[childId].reset();
            a_hasProperty(childId, Cell::WasNewlyCreated) = true;
            a_hasProperty(cellId, Cell::WasRefined) = true;
          }
        }
      }
    }

    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
          MInt cellId = m_azimuthalHaloCells[i][j];
          if(a_level(cellId) != level) continue;
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            const MInt childId = a_childId(cellId, child);
            if(childId < 0) continue;
            if(!a_isHalo(childId)) {
              mTerm(1, AT_, "");
            } else {
              coarseFlag[childId].reset();
              refineFlag[childId].reset();
              a_hasProperty(childId, Cell::WasNewlyCreated) = true;
              a_hasProperty(cellId, Cell::WasRefined) = true;
            }
          }
        }
      }
    }
  }

  //---------------------------(7)

  //----------------------------
  // 8. finalize; synchonize window/halo cell order
  setLevel();
  computeGlobalIds();
  for(auto& rgm : resizeGridMapSolver) {
    rgm(); // make sure grid2solver is large enough
  }
  compactCells(swapCellsSolver);
  // Update after cells have been shifted (min-level cells, globalToLocalId, partitionCellIds, etc)
  computeGlobalIds();
  for(auto& rgm : resizeGridMapSolver) {
    rgm(); // update grid2solver after cell swapping
  }
  if(noDomains() > 1 || m_noPeriodicCartesianDirs > 0) {
    updateHaloCellCollectors();
  }
  computeLeafLevel();
  //---------------------------(8)


  //----------------------------
  // 9. finalize
  if(noNeighborDomains() > 0 || noAzimuthalNeighborDomains() > 0) {
    exchangeProperties();
  }

  // Update local bounding box
  computeLocalBoundingBox(&m_localBoundingBox[0]);
  //---------------------------(9)

  if(treeb().noSolvers() > 1 && noNeighborDomains() > 0) { // Exchange solver info
    // WH_old
    if(m_haloMode > 0) {
      exchangeSolverBitset(&m_tree.solverBits(0));
    } else
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
                                m_tree.size());
  }

  if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                              &m_tree.solverBits(0), m_tree.size());
    // To ensure that azimuthal halo cells are fully refined (have all children)
    // or are not refined at all, solver Bits need to be checked!
    correctAzimuthalSolverBits();
  }


  // TODO labels:GRID,ADAPTATION update noOffspring? not updated during adaptation --> noOffspring required by
  // createHigherLevelExchangeCells, so added above
  // calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), m_tree.size());

  /*
    for ( MInt solver = 0; solver < treeb().noSolvers(); solver++ ){
      MInt debugHalo=0;
      MInt debugWindow=0;
      for ( MInt i = 0; i < noNeighborDomains(); i++ ) {
        for ( MInt j = 0; j < (signed)m_haloCells[i].size(); j++ ) {
          MInt id=m_haloCells[i][j];
          if(m_tree.solver( id, solver )){
            debugHalo++;
          }
        }
        for ( MInt j = 0; j < (signed)m_windowCells[i].size(); j++ ) {
          MInt id=m_windowCells[i][j];
          if(m_tree.solver( id, solver )){
            debugWindow++;
          }
        }

        cerr << " neighbor: " << i << "--------------" << endl;
        cerr <<"count of halos for solver: " << solver << " : " << debugHalo << endl;
        cerr <<"count of windows for solver: " << solver << " : " << debugWindow << endl;
        cerr <<"total size halo: " << m_haloCells[i].size() << " window: " << m_windowCells[i].size() << endl;
      }
    }*/

  // Update the partition cells when writing the next restart file (in case this is not previously
  // done via updatePartitionCells() during balancing)
  // TODO labels:GRID,DLB,ADAPTATION partition files written at final time step might not be usable since partition
  // cells are updated during new grid output
  m_updatedPartitionCells = false;

#ifdef MAIA_GRID_SANITY_CHECKS
  gridSanityChecks();
  checkWindowHaloConsistency(true);
  if(m_haloMode > 0) // WH_old
    checkWindowLayer("meshAdaptation completed: ");
#endif

  // Set adaptation status
  m_wasAdapted = true;

  //  cout << "dom: " << domainId() << " noNeighborDomains: " << noNeighborDomains() << " windowCell: " <<
  //  windowCell_pre_size << " noWindowCells: " << noWindowCells_pre_size<< " isWindowPartLvlAncestor: " <<
  //  isWindowPartLvlAncestor_size << endl;
}


/**
 * \brief Performs mesh adaptation and continuously updates the window/halo cells
 * \author Lennart Schneiders
 *
 * \param[in] sensors Holds the mesh refinement sensor(s) for each cell
 * \param[in] sensorWeight Weighting of this sensor in relation to other sensors if greater than zero, otherwise
 * indicating a true/false-type sensor \param[in] sensorCellFlag Indicator whether a sensor was set for this particular
 * cell \param[out] oldCellId Mapping from new grid cells to their original cellId, or -1 if they did not exist
 * previously \param[out] oldNoCells Number of cells in the map
 */
template <MInt nDim>
void CartesianGrid<nDim>::meshAdaptationDefault(
    vector<vector<MFloat>>& sensors, vector<MFloat>& sensorWeight, vector<bitset<64>>& sensorCellFlag,
    vector<MInt>& sensorSolverId, const std::vector<std::function<void(const MInt)>>& refineCellSolver,
    const std::vector<std::function<void(const MInt)>>& removeChildsSolver,
    const std::vector<std::function<void(const MInt)>>& removeCellSolver,
    const std::vector<std::function<void(const MInt, const MInt)>>& swapCellsSolver,
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver,
    const std::vector<std::function<void()>>& resizeGridMapSolver) {
  TRACE();

  //----------------------------
  // 0. init
  const MInt noCells = treeb().size();
  const MInt maxNoCells = m_maxNoCells;
  const MInt maxLevel = mMin(m_maxRfnmntLvl, m_maxLevel + 1);

  const MInt noSensors = (signed)sensorWeight.size();
  ASSERT(sensorWeight.size() == sensors.size(), "");
  ASSERT(sensorWeight.size() == sensorSolverId.size(), "");
  ScratchSpace<MFloat> sensorCnt(noSensors, 2, AT_, "sensorCnt");
  ScratchSpace<MFloat> sensorThresholds(noSensors, 2, AT_, "sensorThresholds");
  ScratchSpace<bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()>> refineFlag(maxNoCells, AT_, "refineFlag");
  ScratchSpace<bitset<maia::grid::tree::Tree<nDim>::maxNoSolvers()>> coarseFlag(maxNoCells, AT_, "coarseFlag");
  ASSERT(m_freeIndices.empty(), "");
  std::set<MInt>().swap(m_freeIndices);
  sensorCnt.fill(F0);
  for(MInt c = 0; c < maxNoCells; c++) {
    refineFlag[c].reset();
    coarseFlag[c].reset();
  }
  // these flags indicate which cells have changed for the subsequent reinitialization by the solvers
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    a_hasProperty(cellId, Cell::IsToDelete) = false;
    a_hasProperty(cellId, Cell::WasNewlyCreated) = false;
    a_hasProperty(cellId, Cell::WasCoarsened) = false;
    a_hasProperty(cellId, Cell::WasRefined) = false;

    // Reset the window cell flag,
    // is set for all preserved window cells when added to m_windowCells
    a_hasProperty(cellId, Cell::IsWindow) = false;
  }

  //----------------------------
  // 1. exchange sensors and compute RMS values
  // TODO labels:GRID move these parts to the cartesian solver and do this sort of smoothing
  //       on the solver grid in setSensors, after all sensors have been set there!!!!
  //       the cartesiangrid should only get one refine/coarse flags for each cell from the solver!
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    if(a_isToDelete(cellId)) continue;

    // Note: sensors for partition level ancestors might not be set correctly!
    if(a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
      continue;
    }

    for(MInt s = 0; s < noSensors; s++) {
      if(sensorCellFlag[cellId][s]) {
        ASSERT(sensorWeight[s] < F0 || std::fabs(sensors[s][cellId]) < MFloatMaxSqrt,
               "invalid sensor value: " + std::to_string(sensors[s][cellId]));
        sensorCnt(s, 0) += POW2(sensors[s][cellId]);
        sensorCnt(s, 1) += F1;
      }
    }
  }

  if(noSensors > 0) {
    MPI_Allreduce(MPI_IN_PLACE, &sensorCnt(0), 2 * noSensors, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "sensorCnt(0)");
  }

  // TODO labels:GRID,ADAPTATION introduce sensor type to replace this weight hack
  for(MInt s = 0; s < noSensors; s++) {
    if(sensorWeight[s] < F0) { // true/false sensor
      sensorThresholds(s, 0) = -1e-12;
      sensorThresholds(s, 1) = 1e-12;
    } else { // continuous sensor
      sensorThresholds(s, 1) = sensorWeight[s] * sqrt(mMax(1e-14, (sensorCnt(s, 0) / sensorCnt(s, 1))));
      sensorThresholds(s, 0) = m_coarseRatio * sensorThresholds(s, 1);
    }
  }

  //---------------------------(1)

  const MBool hasPartLvlShift = (m_maxPartitionLevelShift > 0);
  const MInt tmpScratchSize = (hasPartLvlShift) ? m_tree.size() : 1;

  ScratchSpace<MBool> isHaloPartLvlAncestor(tmpScratchSize, AT_, "isHaloPartLvlAncestor");
  // TODO labels:GRID,ADAPTATION make this work using less memory, e.g., use one bit for each neighbor?
  ScratchSpace<MBool> isWindowPartLvlAncestor(tmpScratchSize, noNeighborDomains(), AT_, "isWindowPartLvlAncestor");

  isHaloPartLvlAncestor.fill(false);
  isWindowPartLvlAncestor.fill(false);

  // Determine relevant halo partition level ancestor window/halo cells that need to be preserved
  if(m_maxPartitionLevelShift > 0) {
    MInt totalNoWindowCells = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      totalNoWindowCells += noWindowCells(d);
    }
    ScratchSpace<MBool> recvIsHaloPartLvlAncestor(std::max(totalNoWindowCells, 1), AT_, "recvIsHaloPartLvlAncestor");

    for(auto& i : m_partitionLevelAncestorIds) {
      TERMM_IF_NOT_COND(a_hasProperty(i, Cell::IsPartLvlAncestor),
                        "cell is not a partition level ancestor: " + std::to_string(i));
      // Mark halo partition level ancestors (which are part of the missing subtree of the grid)
      if(a_hasProperty(i, Cell::IsHalo)) {
        isHaloPartLvlAncestor[i] = true;
      }
      // Mark all halo childs of partition level ancestors (required for exchange of e.g. global
      // ids!)
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(i, child);
        if(childId > -1 && a_isHalo(childId)) {
          isHaloPartLvlAncestor[childId] = true;
        }
      }
    }

    // Reverse exchange markers from halo cells to corresponding window cells
    if(noNeighborDomains() > 0) {
#ifndef PVPLUGIN
      maia::mpi::reverseExchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &isHaloPartLvlAncestor[0],
                                     &recvIsHaloPartLvlAncestor[0]);
#else
      TERMM(1, "not working for compilation of paraview plugin");
#endif
    }

    MInt cnt = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      for(MInt j = 0; j < noWindowCells(d); j++) {
        const MInt cellId = windowCell(d, j);
        // Store marker for window cell separately for each neighbor domain!
        isWindowPartLvlAncestor(cellId, d) = recvIsHaloPartLvlAncestor[cnt];
        cnt++;
      }
    }
  }

  //----------------------------
  // 2. tag cells for coarsening / refinement

  // init parent cells with coarse flag true and set false if any child intervenes
  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "");
    ASSERT(!a_isToDelete(cellId), "");
    if(a_level(cellId) > m_maxUniformRefinementLevel) {
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(a_hasChildren(cellId, solver)) continue;

        const MInt parentId = a_parentId(cellId, solver);
        ASSERT(parentId > -1, "");
        // Partition level shift, parent is a halo cell and cannot be coarsened
        if(a_isHalo(parentId) || a_hasProperty(parentId, Cell::IsPartLvlAncestor)) continue;

        for(MInt s = 0; s < noSensors; s++) {
          if(sensorSolverId[s] != solver) continue;
          if(sensorCellFlag[cellId][s] || sensorCellFlag[parentId][s]) {
            // only coarsen cells if they are above the newMinLevel!
            if(a_level(cellId) > m_newMinLevel && m_allowCoarsening) {
              coarseFlag[a_parentId(cellId)][solver] = true;
            }
          }
        }
      }
    }
  }

  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "");
    ASSERT(!a_isToDelete(cellId), "");
    const MInt level = a_level(cellId);
    //---
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_tree.solver(cellId, solver)) continue;
      if(a_hasChildren(cellId, solver)) continue;
      const MInt parentId = a_parentId(cellId, solver);
      // Partition level shift, if parent is a halo it cannot be coarsened and e.g.
      // sensors/coarseFlag[parentId] should not be accessed!
      const MBool partLvlAncestorParent = (parentId > -1) ? a_hasProperty(parentId, Cell::IsPartLvlAncestor) : false;

      // cell should be marked for refinement
      MBool refine = false;
      // cell should be marked for coarsened
      MBool coarsen = true;
      // cell should be marked for refinement regardless
      // whether a different sensor marks this cell for coarsening
      MBool forceRefine = false;
      for(MInt s = 0; s < noSensors; s++) {
        if(sensorSolverId[s] != solver) continue;
        if(sensorCellFlag[cellId][s]) {
          coarsen = coarsen && (sensors[s][cellId] < sensorThresholds(s, 0));
        }
        if(level > m_maxUniformRefinementLevel && !partLvlAncestorParent && sensorCellFlag[parentId][s]) {
          coarsen = coarsen && (sensors[s][parentId] < sensorThresholds(s, 0));
        }
        // this sensor would refine the cell
        const MBool refineSensor = sensors[s][cellId] > sensorThresholds(s, 1) && sensorCellFlag[cellId][s];
        refine = refine || refineSensor;

        // a true/false sensor marks the cell for refinement
        if(refineSensor && sensorWeight[s] < F0) {
          forceRefine = true;
        }
      }
      if(level < m_maxRfnmntLvl) {
        refineFlag[cellId][solver] = refine;
      }

      if(level > m_maxUniformRefinementLevel && !partLvlAncestorParent) {
        coarseFlag[parentId][solver] = coarseFlag[parentId][solver] && coarsen;
      }
      // force a refinement of a cell:
      // for all minLevel cells which should be raised to the newMinLevel
      // for the cells forced by any true/false sensor!
      if((forceRefine && level < m_maxRfnmntLvl) || a_level(cellId) < m_newMinLevel) {
        refineFlag[cellId][solver] = true;
        coarseFlag[cellId][solver] = false;
        if(level > m_maxUniformRefinementLevel) {
          coarseFlag[parentId][solver] = false;
        }
      }
    }
  }

  // 3. Exchange since consistency checks need neighbor values
  if(noNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), coarseFlag.getPointer(),
                              m_tree.size());
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), refineFlag.getPointer(),
                              m_tree.size());
  }

  //----------------------------
  // 4. consistency check for internal cells
  if(m_checkRefinementHoles != nullptr) {
    // necessary for solution adaptive refinement (i.e. at shocks)
    const MInt adaptationHoleLimit = nDim;
    MBool gridUpdated = true;

    MInt loopCount = 0;
    while(gridUpdated && loopCount < 25) {
      MInt noChanges = 0;
      loopCount++;
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_checkRefinementHoles[solver]) continue;
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          if(a_isToDelete(cellId)) continue;
          if(a_hasProperty(cellId, Cell::IsHalo)) continue;
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) continue;

          // check, if current cell is a hole in the refinement => refine cell aswell
          if(!refineFlag[cellId][solver]
             || (a_parentId(cellId, solver) > -1 && coarseFlag[a_parentId(cellId, solver)][solver])) {
            MInt checkRefinedNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                const MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(refineFlag[nghbrId][solver] || a_hasChildren(nghbrId, solver)) {
                  checkRefinedNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkRefinedNeighbors > adaptationHoleLimit) {
              refineFlag[cellId][solver] = true;
              coarseFlag[cellId][solver] = false;
              const MInt parentId = a_parentId(cellId, solver);
              if(parentId > -1) coarseFlag[parentId][solver] = false;
              noChanges += 1;
            }
          }

          // check, if current cell coarsening would create a hole in the refinement => do not coarsen cell
          const MInt parentId = a_parentId(cellId, solver);
          if(parentId < 0) continue;
          if(coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            MInt checkRefinedNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(parentId, dir, solver) > 0) {
                const MInt parNghbrId = a_neighborId(parentId, dir, solver);
                if(!coarseFlag[parNghbrId][solver]
                   && (refineFlag[parNghbrId][solver] || a_hasChildren(parNghbrId, solver))) {
                  checkRefinedNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkRefinedNeighbors > adaptationHoleLimit) {
              coarseFlag[parentId][solver] = false;
              noChanges += 1;
            }
          }

          // Check if the coarsening of neighboring cells would create a refined cell island => coarsen cell aswell
          if(!coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            MInt checkCoarseNeighbors = 0;
            MInt currentHoleLimit = adaptationHoleLimit;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(parentId, dir, solver) > 0) {
                const MInt parNghbrId = a_neighborId(parentId, dir, solver);
                if(coarseFlag[parNghbrId][solver]
                   || (!refineFlag[parNghbrId][solver] && !a_hasChildren(parNghbrId, solver))) {
                  checkCoarseNeighbors += 1;
                }
              } else {
                // Lower limit at boundaries
                currentHoleLimit--;
              }
            }
            if(checkCoarseNeighbors > currentHoleLimit) {
              coarseFlag[parentId][solver] = true;
              noChanges += 1;
            }
          }

          // check if the current cell refinement creates a refined cell island => do not refine cell
          if(refineFlag[cellId][solver]) {
            MInt checkRefinedNeighbors = 0;
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                const MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(refineFlag[nghbrId][solver] || a_hasChildren(nghbrId, solver)) {
                  checkRefinedNeighbors += 1;
                }
              }
            }
            if(checkRefinedNeighbors == 0) {
              refineFlag[cellId][solver] = false;
              noChanges += 1;
            }
          }
        }
      }
      if(noChanges == 0) {
        gridUpdated = false;
      }
    }
  }

  // smooth diagonal level jumps,
  // even if already existing in the base grid!
  if(m_diagSmoothing != nullptr) {
    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_diagSmoothing[solver]) continue;
      for(MInt lvl = maxLevel - 1; lvl >= minLevel(); lvl--) {
        for(MInt cellId = 0; cellId < noCells; cellId++) {
          if(a_level(cellId) != lvl) continue;
          if(!a_isLeafCell(cellId)) continue;
          if(!m_tree.solver(cellId, solver)) continue;
          if(refineFlag[cellId][solver]) continue;
          for(MInt i = 0; i < m_noDirs; i++) {
            const MInt nghbrId = a_neighborId(cellId, i, solver);
            if(nghbrId < 0) continue;
            if(a_isLeafCell(nghbrId)) continue;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(nghbrId, child, solver);
              if(childId < 0) continue;
              if(!a_isLeafCell(childId)) {
                refineFlag[cellId][solver] = true;
                coarseFlag[cellId][solver] = false;
                MInt parentId = a_parentId(cellId, solver);
                if(parentId > -1) {
                  coarseFlag[cellId][solver] = false;
                }
                break;
              }
            }
            for(MInt j = 0; j < m_noDirs; j++) {
              if(j / 2 == i / 2) continue;
              const MInt diagNghbrId = a_neighborId(nghbrId, j, solver);
              if(diagNghbrId < 0) continue;
              if(a_isLeafCell(diagNghbrId)) continue;
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                MInt childId = a_childId(diagNghbrId, child, solver);
                if(childId < 0) continue;
                if(!a_isLeafCell(childId)) {
                  refineFlag[cellId][solver] = true;
                  coarseFlag[cellId][solver] = false;
                  MInt parentId = a_parentId(cellId, solver);
                  if(parentId > -1) {
                    coarseFlag[cellId][solver] = false;
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // Additional consistency checks fixing contradictory flags
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    if(a_isToDelete(cellId)) continue;

    for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
      if(!m_tree.solver(cellId, solver)) continue;
      if(a_hasChildren(cellId, solver)) continue;
      const MInt parentId = a_parentId(cellId, solver);
      if(refineFlag[cellId][solver] && coarseFlag[cellId][solver]) {
        refineFlag[cellId][solver] = false;
        coarseFlag[cellId][solver] = false;
      }

      // TODO labels:GRID,ADAPTATION,totest check for multisolver!!!!!
      if(refineFlag[cellId][solver]) {
        if(parentId > -1) {
          if(coarseFlag[parentId][solver]) {
            ASSERT(!a_isHalo(parentId) && !a_hasProperty(parentId, Cell::IsPartLvlAncestor),
                   "Partition level ancestor parent cell has coarseFlag set.");
            coarseFlag[parentId][solver] = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(parentId, child, solver);
              if(childId > -1) {
                refineFlag[childId][solver] = false;
              }
            }
          }
        }
      }
    }
  }
  //---------------------------(4)

  // ensure the same refinement between two different solvers in the combined/overlapping region!
  if(m_noIdenticalSolvers > 0) {
    for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
      for(MInt pair = 0; pair < m_noIdenticalSolvers; pair++) {
        const MInt slaveId = m_identicalSolvers[pair * 2];
        const MInt masterId = m_identicalSolvers[pair * 2 + 1];
        if(treeb().solver(cellId, slaveId) && treeb().solver(cellId, masterId)) {
          if(!coarseFlag[cellId][masterId]) {
            coarseFlag[cellId][slaveId] = coarseFlag[cellId][masterId];
          }
          if(m_identicalSolverLvlJumps[pair] == 0) {
            refineFlag[cellId][slaveId] = refineFlag[cellId][masterId];
          } else if(m_identicalSolverLvlJumps[pair] == 1) {
            if(a_level(cellId) < m_identicalSolverMaxLvl[pair] && a_hasChildren(cellId, masterId)) {
              if(!a_hasChildren(cellId, slaveId)) {
                refineFlag[cellId][slaveId] = true;
              }
              coarseFlag[cellId][slaveId] = false;
            }
          } else {
            if(a_level(cellId) < m_identicalSolverMaxLvl[pair] && a_hasChildren(cellId, masterId)
               && a_childId(cellId, 0) > -1 && treeb().solver(a_childId(cellId, 0), masterId)
               && a_hasChildren(a_childId(cellId, 0), masterId)) {
              if(!a_hasChildren(cellId, slaveId)) {
                refineFlag[cellId][slaveId] = true;
              }
              coarseFlag[cellId][slaveId] = false;
            }
          }
        }
      }
    }
  }


  //----------------------------
  // 5. Coarsening step: loop from maxLevel to minLevel and exchange coarseFlags to retrieve
  //    consistent refinement at domain interfaces; update window/halo collectors
  for(MInt level = m_maxRfnmntLvl - 1; level >= m_maxUniformRefinementLevel; level--) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      // TODO labels:GRID,ADAPTATION,totest check for multisolver setting.. neighborIds solver aware!!
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(!a_hasChildren(cellId, solver)) continue;

        // Partition level ancestor cell cannot be coarsened (exept if its a partition cell)
        if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
          coarseFlag[cellId][solver] = false;
        }

        if(coarseFlag[cellId][solver]) {
          if(!coarsenCheck(cellId, solver)) {
            coarseFlag[cellId][solver] = false;
          }
        }
      }
    }

    if(noNeighborDomains() > 0) {
      //      exchangeSolverBitset(coarseFlag.data());
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), coarseFlag.data(),
                                m_tree.size());
    }

    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      // TODO labels:GRID,ADAPTATION,totest check for multisolver setting.. neighborIds solver aware!!
      if(level < m_maxRfnmntLvl - 1) {
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!m_tree.solver(cellId, solver)) continue;
          if(a_hasChildren(cellId, solver)) continue;
          if(!refineFlag[cellId][solver] && !coarseFlag[cellId][solver]) {
            for(MInt dir = 0; dir < m_noDirs; dir++) {
              if(a_hasNeighbor(cellId, dir, solver) > 0) {
                MInt nghbrId = a_neighborId(cellId, dir, solver);
                if(m_azimuthalPer && a_hasProperty(nghbrId, Cell::IsPeriodic)) continue;
                if(a_hasChildren(nghbrId, solver) > 0) {
                  for(MInt child = 0; child < m_maxNoChilds; child++) {
                    if(!childCode[dir][child]) continue;
                    MInt childId = a_childId(nghbrId, child, solver);
                    if(childId < 0) continue;
                    if(refineFlag[childId][solver]) {
                      refineFlag[cellId][solver] = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if(noNeighborDomains() > 0) {
      //      exchangeSolverBitset(refineFlag.data());
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), refineFlag.data(),
                                m_tree.size());
    }

    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(m_azimuthalPer && a_hasProperty(cellId, Cell::IsPeriodic)) continue;

      // Partition level ancestor cell cannot be coarsened (exept if its a partition cell)
      if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
        continue;
      }

      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(coarseFlag[cellId][solver]) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt childId = a_childId(cellId, child, solver);
            if(childId > -1) {
              coarseFlag[childId][solver] = false;
              refineFlag[childId][solver] = false;
            }
          }
          if(a_hasChildren(cellId, solver)) {
            removeChildsSolver[solver](cellId);
            // call removeChilds() function in each solver before child links are deleted
            refineFlag[cellId][solver] = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId > -1) {
                treeb().solver(childId, solver) = false;
              }
            }
          }
        }
      }
      // check if the child can also be removed from the cartesian/combined grid
      // if the child is no-longer in any of the solvers, the grid child will be removed!
      if(a_noChildren(cellId) > 0) {
        MBool rmChilds = true;
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(a_hasChildren(cellId, solver)) {
            rmChilds = false;
          }
        }
        if(rmChilds) {
          removeChilds(cellId);
        }
      }
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      std::set<MInt> deleted;
      // WH_old
      if(m_haloMode > 0) {
        for(auto it_ = m_windowLayer_[i].begin(); it_ != m_windowLayer_[i].end();) {
          if(a_isToDelete(it_->first)) {
            it_ = m_windowLayer_[i].erase(it_);
          } else {
            // We might have deleted a cell from a solver
            for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
              if(it_->second.get(solver) < WINDOWLAYER_MAX /*<=m_noSolverHaloLayers[solver]*/
                 && !m_tree.solver(it_->first, solver))
                it_->second.set(solver, WINDOWLAYER_MAX);
              if(it_->second.get(solver) > m_noSolverHaloLayers[solver]) it_->second.set(solver, WINDOWLAYER_MAX);
            }
            if(it_->second.all()) {
              deleted.insert(it_->first);
              it_ = m_windowLayer_[i].erase(it_);
            } else
              ++it_;
          }
        }
      }

      // 1. coarsen window cells
      vector<MInt> oldWindowVec(m_windowCells[i]);
#ifndef NDEBUG
      std::set<MInt> windCellSet;
      std::set<MInt> haloCellSet;
#endif

      std::vector<MInt>().swap(m_windowCells[i]);
      auto it = oldWindowVec.begin();
      for(it = oldWindowVec.begin(); it < oldWindowVec.end(); it++) {
        const MInt cellId = *it;
        if(!a_isToDelete(cellId) && (deleted.find(cellId) == deleted.end() || m_haloMode == 0)) { // WH_old
          ASSERT(cellId < treeb().size(), "");

#ifndef NDEBUG
          // TODO labels:GRID,ADAPTATION fix duplicate window cells in periodic cases!
          ASSERT(m_noPeriodicCartesianDirs > 0 || windCellSet.find(cellId) == windCellSet.end(),
                 "duplicate window cell: " + std::to_string(cellId));
          windCellSet.insert(cellId);
#endif

          m_windowCells[i].push_back(cellId);
        }
      }
      oldWindowVec.clear();

      // 2. coarsen halo cells
      vector<MInt> oldHaloVec(m_haloCells[i]);
      std::vector<MInt>().swap(m_haloCells[i]);
      for(it = oldHaloVec.begin(); it < oldHaloVec.end(); it++) {
        const MInt cellId = *it;
        if(!a_isToDelete(cellId)) {
          ASSERT(cellId < treeb().size(), "");

#ifndef NDEBUG
          ASSERT(haloCellSet.find(cellId) == haloCellSet.end(), "duplicate halo cell: " + std::to_string(cellId));
          haloCellSet.insert(cellId);
#endif

          if(m_tree.solverBits(cellId).any() || m_haloMode == 0) m_haloCells[i].push_back(cellId);
        }
      }
      oldHaloVec.clear();
    }
  }
  //---------------------------(5)


  //----------------------------
  // 6. Remove all halo cells, for levels > maxUniformRefinementLevel
  //    since the number of halo layers is inconsistent at this point
  //    Note: preserve halo partition level ancestors of the missing grid subtree in case of a
  //    partition level shift

  {
    for(MInt level = maxLevel - 1; level >= m_maxUniformRefinementLevel; level--) {
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        if(a_isToDelete(cellId)) continue;
        if(!a_hasProperty(cellId, Cell::IsHalo)) continue;
        if(a_level(cellId) == level) {
          if(a_noChildren(cellId) > 0) {
            // Preserve partition level ancestor halos (and their childs) if they are part of the
            // missing subtree of the grid in case of a partition level shift
            if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsPartitionCell)) {
              if(std::find(m_partitionLevelAncestorIds.begin(), m_partitionLevelAncestorIds.end(), cellId)
                 != m_partitionLevelAncestorIds.end()) {
                continue;
              }
            }

            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId < 0) continue;
              refineFlag[childId].reset();
              coarseFlag[childId].reset();
            }
            for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
              if(!m_tree.solver(cellId, solver)) continue;
              if(a_hasChildren(cellId, solver)) {
                removeChildsSolver[solver](cellId);
              }
            }
            removeChilds(cellId);
            coarseFlag[cellId].reset();
            refineFlag[cellId].reset();
          }
        }
      }
    }
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      vector<MInt> oldWindowVec(m_windowCells[i]);
      std::vector<MInt>().swap(m_windowCells[i]);
      auto it = oldWindowVec.begin();
      for(it = oldWindowVec.begin(); it < oldWindowVec.end(); it++) {
        const MInt cellId = *it;

        // Add partition level ancestor window cells if required
        const MBool addWindow = (hasPartLvlShift) ? isWindowPartLvlAncestor(cellId, i) : false;

        if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addWindow) {
          m_windowCells[i].push_back(cellId);
          // Set window cell flag for preserved window cells
          a_hasProperty(cellId, Cell::IsWindow) = true;
        }
      }
      oldWindowVec.clear();

      // WH_old
      if(m_haloMode > 0) {
        for(auto it_ = m_windowLayer_[i].begin(); it_ != m_windowLayer_[i].end();) {
          const MInt cellId = it_->first;

          // Add partition level ancestor window cells if required
          const MBool addWindow = (hasPartLvlShift) ? isWindowPartLvlAncestor(cellId, i) : false;

          if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addWindow)
            ++it_;
          else
            it_ = m_windowLayer_[i].erase(it_);
        }
      }

      vector<MInt> oldHaloVec(m_haloCells[i]);
      std::vector<MInt>().swap(m_haloCells[i]);
      for(it = oldHaloVec.begin(); it < oldHaloVec.end(); it++) {
        const MInt cellId = *it;

        // Add partition level ancestor halo cells if required
        const MBool addHalo = (hasPartLvlShift) ? isHaloPartLvlAncestor[cellId] : false;

        if((a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) || addHalo) {
          ASSERT(cellId > -1 && cellId < treeb().size(), to_string(cellId) + " " + to_string(treeb().size()));
          m_haloCells[i].push_back(cellId);
        }
      }
      oldHaloVec.clear();
    }

    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        vector<MInt> oldAzimuthalWindowVec(m_azimuthalWindowCells[i]);

        m_azimuthalWindowCells[i].clear();
        auto it = oldAzimuthalWindowVec.begin();
        for(MInt j = 0; j < (signed)oldAzimuthalWindowVec.size(); j++) {
          // Remove higher level azimuthal window cells from collector
          // If an azimuthal halo cell on a higher level is mapped to an
          // azimuthal window cell on a lower level. This azimuthal window cell
          // musst also be removed from the collector
          const MInt cellId = oldAzimuthalWindowVec[j];
          MBool keepWindow = true;
          if(find(m_azimuthalHigherLevelConnectivity[i].begin(), m_azimuthalHigherLevelConnectivity[i].end(), j)
             != m_azimuthalHigherLevelConnectivity[i].end()) {
            keepWindow = false;
          }
          if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel && keepWindow) {
            m_azimuthalWindowCells[i].push_back(cellId);
            a_hasProperty(cellId, Cell::IsWindow) = true;
          }
        }
        oldAzimuthalWindowVec.clear();

        vector<MInt> oldAzimuthalHaloVec(m_azimuthalHaloCells[i]);
        m_azimuthalHaloCells[i].clear();
        for(it = oldAzimuthalHaloVec.begin(); it < oldAzimuthalHaloVec.end(); it++) {
          const MInt cellId = *it;
          if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
            ASSERT(cellId > -1 && cellId < treeb().size(), "");
            m_azimuthalHaloCells[i].push_back(cellId);
          }
        }
        oldAzimuthalHaloVec.clear();
      }

      // Remove higher level unmapped azimuthal halo cells
      vector<MInt> oldUnmappedHaloVec(m_azimuthalUnmappedHaloCells);
      vector<MInt> oldUnmappedHaloDomVec(m_azimuthalUnmappedHaloDomains);
      m_azimuthalUnmappedHaloCells.clear();
      m_azimuthalUnmappedHaloDomains.clear();
      for(MUint i = 0; i < oldUnmappedHaloVec.size(); i++) {
        MInt cellId = oldUnmappedHaloVec[i];
        if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
          ASSERT(cellId > -1 && cellId < treeb().size(), "");
          m_azimuthalUnmappedHaloCells.push_back(cellId);
          m_azimuthalUnmappedHaloDomains.push_back(oldUnmappedHaloDomVec[i]);
        }
      }
      oldUnmappedHaloVec.clear();
      oldUnmappedHaloDomVec.clear();

      // Update gridBndryCells
      vector<MInt> oldGridBndryVec(m_gridBndryCells);
      m_gridBndryCells.clear();
      for(MUint i = 0; i < oldGridBndryVec.size(); i++) {
        MInt cellId = oldGridBndryVec[i];
        if(a_level(cellId) > -1 && a_level(cellId) <= m_maxUniformRefinementLevel) {
          m_gridBndryCells.push_back(cellId);
        }
      }
      oldGridBndryVec.clear();
    }
  }

  //---------------------------(6)
  //----------------------------
  // 7. Refinement step: loop from minLevel to maxLevel and exchange refineFlags to retrieve
  //    a consistent refinement at domain interfaces; update window/halo connectivity at each level
  for(MInt level = m_maxUniformRefinementLevel; level < maxLevel; level++) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(a_level(cellId) != level) continue;
      if(a_isToDelete(cellId)) continue;
      if(a_hasProperty(cellId, Cell::IsHalo)) continue;

      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(!m_tree.solver(cellId, solver)) continue;
        if(a_hasChildren(cellId, solver)) continue;
        if(refineFlag[cellId][solver]) {
          if(!refineCheck(cellId, solver, cellOutsideSolver)) {
            refineFlag[cellId][solver] = 0;
          }
        }
      }

      if(refineFlag[cellId].any()) {
        refineCell(cellId, nullptr, false, cellOutsideSolver, refineFlag[cellId]);

        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(!a_hasChildren(cellId, solver) && false) { // sohel false
            // this is somehow destroying multisolver cases
            // if one solver already has a cell the other solver doesnt get it even if it needs it..
            // i think the solver flag is not set correctly in the refineCell function of cartesiangrid
            refineFlag[cellId][solver] = false;
          }
          if(refineFlag[cellId][solver]) {
            // Set solver flag for child cells
            setChildSolverFlag(cellId, solver, cellOutsideSolver[solver]);
            refineCellSolver[solver](cellId); // call refineCell() function in each solver
          }
        }
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          MInt childId = a_childId(cellId, child);
          if(childId > -1 && a_hasProperty(childId, Cell::WasNewlyCreated)) {
            // Timw: only reset the child-flags if the cell was neawly created
            // otherwise this stops a refinement of the childs of a different solver!
            coarseFlag[childId].reset();
            refineFlag[childId].reset();
          }
        }
      }
    }

    setLevel();

    computeGlobalIds();

    // WH_old
    if(m_haloMode > 0)
      // noOffspring required in createHigherLevelExchangeCells
      calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), m_tree.size());

    // Update child ids of partition level ancestors (stored as global ids)
    m_partitionLevelAncestorChildIds.clear();
    const MInt noPartLvlAncestorIds = m_partitionLevelAncestorIds.size();
    for(MInt i = 0; i < noPartLvlAncestorIds; i++) {
      const MInt cellId = m_partitionLevelAncestorIds[i];

      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(cellId, child);
        const MLong childGlobalId = (childId > -1) ? a_globalId(childId) : -1;
        m_partitionLevelAncestorChildIds.push_back(childGlobalId);
      }
    }

    if(noDomains() > 1 || m_noPeriodicCartesianDirs > 0) {
      updateHaloCellCollectors();

      // WH_old
      if(m_haloMode > 0) {
        createHigherLevelExchangeCells(level, true, refineCellSolver, removeCellSolver, level == maxLevel - 1);
        // createHigherLevelExchangeCells(level, swapCellsSolver, refineCellSolver);
      } else
        createHigherLevelExchangeCells_old(level, refineCellSolver);
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        MInt cellId = m_haloCells[i][j];
        if(a_level(cellId) != level) continue;
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          const MInt childId = a_childId(cellId, child);
          if(childId < 0) continue;

          // Partition level shift: do not reset flags if the child of a halo cell is an internal
          // cell and already existed!
          if(!a_isHalo(childId)) {
            ASSERT(a_hasProperty(childId, Cell::IsPartLvlAncestor) || a_hasProperty(childId, Cell::IsPartitionCell),
                   "child is not a partition level ancestor or a partition cell");
            a_hasProperty(cellId, Cell::WasRefined) = true;
          } else {
            coarseFlag[childId].reset();
            refineFlag[childId].reset();
            a_hasProperty(childId, Cell::WasNewlyCreated) = true;
            a_hasProperty(cellId, Cell::WasRefined) = true;
          }
        }
      }
    }

    if(m_azimuthalPer) {
      for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
        for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
          MInt cellId = m_azimuthalHaloCells[i][j];
          if(a_level(cellId) != level) continue;
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            const MInt childId = a_childId(cellId, child);
            if(childId < 0) continue;
            if(!a_isHalo(childId)) {
              mTerm(1, AT_, "");
            } else {
              coarseFlag[childId].reset();
              refineFlag[childId].reset();
              a_hasProperty(childId, Cell::WasNewlyCreated) = true;
              a_hasProperty(cellId, Cell::WasRefined) = true;
            }
          }
        }
      }
    }
  }

  //---------------------------(7)

  //----------------------------
  // 8. finalize; synchonize window/halo cell order
  setLevel();
  computeGlobalIds();
  for(auto& rgm : resizeGridMapSolver) {
    rgm(); // make sure grid2solver is large enough
  }
  compactCells(swapCellsSolver);
  // Update after cells have been shifted (min-level cells, globalToLocalId, partitionCellIds, etc)
  computeGlobalIds();
  for(auto& rgm : resizeGridMapSolver) {
    rgm(); // update grid2solver after cell swapping
  }
  if(noDomains() > 1 || m_noPeriodicCartesianDirs > 0) {
    updateHaloCellCollectors();
  }
  computeLeafLevel();
  //---------------------------(8)


  //----------------------------
  // 9. finalize
  if(noNeighborDomains() > 0 || noAzimuthalNeighborDomains() > 0) {
    exchangeProperties();
  }

  // Update local bounding box
  computeLocalBoundingBox(&m_localBoundingBox[0]);
  //---------------------------(9)

  if(treeb().noSolvers() > 1 && noNeighborDomains() > 0) { // Exchange solver info
    // WH_old
    if(m_haloMode > 0) {
      exchangeSolverBitset(&m_tree.solverBits(0));
    } else
      maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.solverBits(0),
                                m_tree.size());
  }

  if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                              &m_tree.solverBits(0), m_tree.size());
    // To ensure that azimuthal halo cells are fully refined (have all children)
    // or are not refined at all, solver Bits need to be checked!
    correctAzimuthalSolverBits();
  }


  // TODO labels:GRID,ADAPTATION update noOffspring? not updated during adaptation --> noOffspring required by
  // createHigherLevelExchangeCells, so added above
  // calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), m_tree.size());

  /*
    for ( MInt solver = 0; solver < treeb().noSolvers(); solver++ ){
      MInt debugHalo=0;
      MInt debugWindow=0;
      for ( MInt i = 0; i < noNeighborDomains(); i++ ) {
        for ( MInt j = 0; j < (signed)m_haloCells[i].size(); j++ ) {
          MInt id=m_haloCells[i][j];
          if(m_tree.solver( id, solver )){
            debugHalo++;
          }
        }
        for ( MInt j = 0; j < (signed)m_windowCells[i].size(); j++ ) {
          MInt id=m_windowCells[i][j];
          if(m_tree.solver( id, solver )){
            debugWindow++;
          }
        }

        cerr << " neighbor: " << i << "--------------" << endl;
        cerr <<"count of halos for solver: " << solver << " : " << debugHalo << endl;
        cerr <<"count of windows for solver: " << solver << " : " << debugWindow << endl;
        cerr <<"total size halo: " << m_haloCells[i].size() << " window: " << m_windowCells[i].size() << endl;
      }
    }*/

  // Update the partition cells when writing the next restart file (in case this is not previously
  // done via updatePartitionCells() during balancing)
  // TODO labels:GRID,DLB,ADAPTATION partition files written at final time step might not be usable since partition
  // cells are updated during new grid output
  m_updatedPartitionCells = false;

#ifdef MAIA_GRID_SANITY_CHECKS
  gridSanityChecks();
  checkWindowHaloConsistency(true);
  if(m_haloMode > 0) // WH_old
    checkWindowLayer("meshAdaptation completed: ");
#endif

  // Set adaptation status
  m_wasAdapted = true;
}


// --------------------------------------------------------------------------------------


/**
 * \brief Exchange properties of window/halo cells
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::exchangeProperties() {
  ScratchSpace<MInt> isPeriodic(m_tree.size() - m_noInternalCells, AT_, "isPeriodic");
  for(MInt cellId = m_noInternalCells; cellId < m_tree.size(); cellId++) {
    ASSERT(a_hasProperty(cellId, Cell::IsHalo), "not a halo cell");
    isPeriodic[cellId - m_noInternalCells] = a_hasProperty(cellId, Cell::IsPeriodic);
  }
  maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.properties(0),
                            m_tree.size());
  if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                              &m_tree.properties(0), m_tree.size());
  }
  for(MInt cellId = m_noInternalCells; cellId < m_tree.size(); cellId++) {
    a_hasProperty(cellId, Cell::IsHalo) = true;
    a_hasProperty(cellId, Cell::IsWindow) = false;
    a_hasProperty(cellId, Cell::IsPeriodic) = isPeriodic[cellId - m_noInternalCells];
  }
}


// --------------------------------------------------------------------------------------


/**
 * \brief Removes all holes in the cell collector and moves halo cells to the back of the collector
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::compactCells(
    const std::vector<std::function<void(const MInt, const MInt)>>& swapCellsSolver) {
  MIntScratchSpace oldCellId(m_maxNoCells, AT_, "oldCellId");
  MIntScratchSpace isToDelete(m_maxNoCells, AT_, "isToDelete");
  oldCellId.fill(-1);
  isToDelete.fill(0);
  for(auto& i : m_freeIndices) {
    isToDelete(i) = 1;
  }
  std::set<MInt>().swap(m_freeIndices);

  // 1. determine number of cells and internal cells
  MInt noCells = 0;
  m_noInternalCells = 0;
  oldCellId.fill(-1);
  for(MInt cellId = 0; cellId < treeb().size(); cellId++) {
    if(isToDelete(cellId)) continue;
    oldCellId(cellId) = cellId;
    if(!a_hasProperty(cellId, Cell::IsHalo)) m_noInternalCells++;
    noCells++;
  }

  // 2. remove holes created by previously deleted cells and move halo cells to the back
  MInt otherId = treeb().size() - 1;
  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    if(isToDelete(cellId) || a_hasProperty(cellId, Cell::IsHalo)) {
      while(isToDelete(otherId) || a_hasProperty(otherId, Cell::IsHalo)) {
        otherId--;
      }
      ASSERT(cellId < otherId, "");
      swapCells(cellId, otherId);
      for(auto& swp : swapCellsSolver) {
        swp(cellId, otherId); // call swapCells() function in each solver
      }
      std::swap(oldCellId(cellId), oldCellId(otherId));
      std::swap(isToDelete(cellId), isToDelete(otherId));
      ASSERT(!a_hasProperty(cellId, Cell::IsHalo), "");
      ASSERT(isToDelete(otherId) || a_hasProperty(otherId, Cell::IsHalo), "");
    }
    ASSERT(!a_hasProperty(cellId, Cell::IsHalo) && !isToDelete(cellId), "");
    ASSERT(a_level(cellId) > -1, "");
  }


  // 3. remove holes in the range of halo cells
  otherId = treeb().size() - 1;
  for(MInt cellId = m_noInternalCells; cellId < noCells; cellId++) {
    if(isToDelete(cellId)) {
      while(isToDelete(otherId)) {
        otherId--;
      }
      ASSERT(cellId < otherId, "");
      ASSERT(otherId >= noCells, "");
      ASSERT(a_hasProperty(otherId, Cell::IsHalo) && !isToDelete(otherId), "");
      swapCells(cellId, otherId);
      for(auto& swp : swapCellsSolver) {
        swp(cellId, otherId); // call swapCells() function in each solver
      }
      std::swap(oldCellId(cellId), oldCellId(otherId));
      std::swap(isToDelete(cellId), isToDelete(otherId));
    }
    ASSERT(a_hasProperty(cellId, Cell::IsHalo) && !isToDelete(cellId), "");
    ASSERT(a_level(cellId) > -1, "");
  }

  // 4. update window/halo cell ids
  treeb().size(noCells);
  MIntScratchSpace newCellId(m_maxNoCells, AT_, "newCellId");
  newCellId.fill(-1);
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(oldCellId(cellId) < 0) continue;
    newCellId(oldCellId(cellId)) = cellId;
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      MInt cellId = m_haloCells[i][j];
      if(newCellId(cellId) > -1) {
        m_haloCells[i][j] = newCellId(cellId);
      }
    }
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      MInt cellId = m_windowCells[i][j];
      if(newCellId(cellId) > -1) {
        m_windowCells[i][j] = newCellId(cellId);
      }
    }

    // TODO_SS labels:GRID,toenhance maybe think of something more efficient
    if(m_haloMode > 0) { // WH_old
      std::unordered_map<MInt, M32X4bit<true>> windowLayerBak_(m_windowLayer_[i]);
      m_windowLayer_[i].clear();
      for(auto m : windowLayerBak_) {
        if(newCellId(m.first) > -1) {
          m_windowLayer_[i].insert({newCellId(m.first), m.second});
        }
      }
    }
  }

  if(m_azimuthalPer) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        MInt cellId = m_azimuthalHaloCells[i][j];
        if(newCellId(cellId) > -1) {
          m_azimuthalHaloCells[i][j] = newCellId(cellId);
        }
      }
      for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
        MInt cellId = m_azimuthalWindowCells[i][j];
        if(newCellId(cellId) > -1) {
          m_azimuthalWindowCells[i][j] = newCellId(cellId);
        }
      }
    }
    for(MInt c = 0; c < noAzimuthalUnmappedHaloCells(); c++) {
      MInt cellId = azimuthalUnmappedHaloCell(c);
      if(newCellId(cellId) > -1) {
        m_azimuthalUnmappedHaloCells[c] = newCellId(cellId);
      }
    }
  }

  // Update partition level ancestor data
  {
    std::vector<MInt> oldPartLevelAncestorIds(m_partitionLevelAncestorIds);
    const MInt noPartLvlAncestorIds = oldPartLevelAncestorIds.size();

    // m_partitionLevelAncestorChildIds.clear();
    // Note: the variables below seem to be not required at the moment
    std::vector<MInt>().swap(m_partitionLevelAncestorIds);
    std::vector<MInt>().swap(m_partitionLevelAncestorNghbrDomains);
    std::vector<std::vector<MInt>>().swap(m_partitionLevelAncestorHaloCells);
    std::vector<std::vector<MInt>>().swap(m_partitionLevelAncestorWindowCells);
    m_noHaloPartitionLevelAncestors = 0;

    // Update partition level ancestor ids and their child (global) ids
    for(MInt i = 0; i < noPartLvlAncestorIds; i++) {
      const MInt oldId = oldPartLevelAncestorIds[i];
      const MInt newId = newCellId(oldId);
      m_partitionLevelAncestorIds.push_back(newId);

      // TODO labels:GRID Since we store the global ids of the childs, there is no need to update them.
      //      If we would do so, we need to keep in mind that some children may reside on
      //      neighbor domains, so we need to communicate.
      //      for(MInt child = 0; child < m_maxNoChilds; child++) {
      //        const MInt childId = a_childId(newId, child);
      //        const MLong childGlobalId = (childId > -1) ? a_globalId(childId) : -1;
      //        m_partitionLevelAncestorChildIds.push_back(childGlobalId);
      //      }
    }
  }
}


// --------------------------------------------------------------------------------------


/**
 * \brief swap two cells; the window/halo cell arrays have to be update elsewhere for performance reasons
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::swapCells(const MInt cellId0, const MInt cellId1) {
  static constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};
  if(cellId1 == cellId0) return;
  ASSERT(cellId1 > -1 && cellId0 > -1 && cellId1 < treeb().size() && cellId0 < treeb().size(),
         "Invalid cell range " << cellId1 << " " << cellId0 << " " << treeb().size());

  // TODO labels:GRID use swap() of grid tree class instead?
  // swap parent-child links, the sequence is important in case cellId0 and cellId1 are in parent-child relation!
  MInt parentId0 = a_parentId(cellId0);
  MInt parentId1 = a_parentId(cellId1);
  MInt child0 = -1;
  MInt child1 = -1;
  if(parentId0 > -1) {
    for(MInt k = 0; k < m_maxNoChilds; k++) {
      if(a_childId(parentId0, k) == cellId0) {
        child0 = k;
        break;
      }
    }
  }
  if(parentId1 > -1) {
    for(MInt k = 0; k < m_maxNoChilds; k++) {
      if(a_childId(parentId1, k) == cellId1) {
        child1 = k;
        break;
      }
    }
  }
  for(MInt k = 0; k < m_maxNoChilds; k++) {
    if(a_childId(cellId0, k) > -1) {
      ASSERT(a_parentId(a_childId(cellId0, k)) == cellId0, "Wrong parent!");
      a_parentId(a_childId(cellId0, k)) = cellId1;
    }
  }
  for(MInt k = 0; k < m_maxNoChilds; k++) {
    if(a_childId(cellId1, k) > -1) {
      ASSERT(a_parentId(a_childId(cellId1, k)) == cellId1, "Wrong parent!");
      a_parentId(a_childId(cellId1, k)) = cellId0;
    }
  }
  if(child0 > -1) a_childId(parentId0, child0) = cellId1;
  if(child1 > -1) a_childId(parentId1, child1) = cellId0;

  for(MInt k = 0; k < m_maxNoChilds; k++) {
    std::swap(a_childId(cellId1, k), a_childId(cellId0, k));
  }
  std::swap(a_parentId(cellId0), a_parentId(cellId1));

  // swap forward/reverse neighbor links (in opposite directions, otherwise not working when cellId0 and cellId1 are
  // neighbors)
  for(MInt dir = 0; dir < m_noDirs; dir++) {
    MInt nghbrId0 = a_neighborId(cellId0, dir);
    MInt nghbrId1 = a_neighborId(cellId1, revDir[dir]);
    if(nghbrId0 > -1) {
      ASSERT(a_neighborId(nghbrId0, revDir[dir]) == cellId0, "Reverse link dead");
      a_neighborId(nghbrId0, revDir[dir]) = cellId1;
    }
    if(nghbrId1 > -1) {
      ASSERT(a_neighborId(nghbrId1, dir) == cellId1, "Reverse link dead");
      a_neighborId(nghbrId1, dir) = cellId0;
    }
  }
  for(MInt dir = 0; dir < m_noDirs; dir++) {
    std::swap(a_neighborId(cellId1, dir), a_neighborId(cellId0, dir));
  }


  // swap other data
  std::swap(a_level(cellId1), a_level(cellId0));
  std::swap(m_tree.properties(cellId0), m_tree.properties(cellId1));
  std::swap(m_tree.solverBits(cellId0), m_tree.solverBits(cellId1));
  std::swap(a_globalId(cellId1), a_globalId(cellId0));
  std::swap(m_tree.leafCellBits(cellId1), m_tree.leafCellBits(cellId0));
  std::swap(a_weight(cellId1), a_weight(cellId0));
  std::swap(a_workload(cellId1), a_workload(cellId0));
  std::swap(a_noOffsprings(cellId1), a_noOffsprings(cellId0));
  for(MInt i = 0; i < nDim; i++) {
    std::swap(a_coordinate(cellId1, i), a_coordinate(cellId0, i));
  }
}


// --------------------------------------------------------------------------------------

/**
 * \brief checks if the given cell's children may be removed
 * \author Lennart Schneiders, Andreas Lintermann
 * \date 23.04.2019
 *
 * The check is performed by:
 *   1. running over all children and
 *     1.1 checking if any of the existing children has children (do not coarsen)
 *     1.2 (switch) diagonal level jumps exist (do not coarsen)
 *   2. additionally checking if the current cell has a missing neighbor. In this case and if it has children, it has to
 *be a boundary cell and it's children need to be marked as non-coarsenable cells (past function call).
 *
 * \param[in] cellId the cell id to check
 * \param[in] solver the solver id
 *
 **/
template <MInt nDim>
MBool CartesianGrid<nDim>::coarsenCheck(const MInt cellId, const MInt solver) {
  // Enables grid consistency checks(coarsening) for lattice boltzmann methods
  MBool extendCheckTo2ndNeighbor = m_lbGridChecks;
  MBool restrictDiagonalLevelJumps = m_lbGridChecks;
  MBool lbBndCheck = m_lbGridChecks;

  ASSERT(m_tree.solver(cellId, solver), "");

  // Check that coarse neighbor has all children
  // --> checks if coarse neighbor contains a boundary
  // --> checks that no interface parents with missing children are created
  // Use enhanced computations of directions in above check (See refine check)
  // for checking diagonal neighbors (corners) in 3D, too.
  if(lbBndCheck) {
    for(MInt dir = 0; dir < m_noDirs; dir++) {
      if(!a_hasNeighbor(cellId, dir, solver)) /*continue*/
        return false;                         // Avoid coarsening of boundary cells
      MInt cartesianNeighbor = a_neighborId(cellId, dir, solver);
      if(a_noChildren(cartesianNeighbor) != m_maxNoChilds && a_noChildren(cartesianNeighbor) != 0) return false;
      MInt d0 = dir / 2;
      MInt d1 = (d0 + 1) % nDim;
      for(MInt p = 0; p < 2; p++) {
        MInt dir1 = 2 * d1 + p;
        if(!a_hasNeighbor(cartesianNeighbor, dir1, solver)) continue;
        MInt diagonalNeighbor = a_neighborId(cartesianNeighbor, dir1, solver);
        if(a_noChildren(diagonalNeighbor) != m_maxNoChilds && a_noChildren(diagonalNeighbor) != 0) return false;
        IF_CONSTEXPR(nDim == 3) {
          MInt d2 = (d1 + 1) % nDim;
          for(MInt q = 0; q < 2; q++) {
            MInt dir2 = 2 * d2 + q;
            if(!a_hasNeighbor(diagonalNeighbor, dir2, solver)) continue;
            MInt cornerNeighbor = a_neighborId(diagonalNeighbor, dir2, solver);
            if(a_noChildren(cornerNeighbor) != m_maxNoChilds && a_noChildren(cornerNeighbor) != 0) return false;
          }
        }
      }
    }
  }


  // 1. run over all children
  for(MInt child = 0; child < m_maxNoChilds; child++) {
    MInt childId = a_childId(cellId, child, solver);

    if(childId < 0) continue;
    if(a_hasChildren(childId, solver)) return false;
    for(MInt i = 0; i < m_noDirs; i++) {
      if(childCode[i][child]) continue;
      if(a_hasNeighbor(childId, i, solver) > 0) {
        // 1.1 check if any of the existing children has children (do not coarsen)
        if(a_hasChildren(a_neighborId(childId, i, solver), solver) > 0) {
          return false;
        }

        // Extend to second Cartesian neighbor of child to prevent double interpolation in LB case
        if(extendCheckTo2ndNeighbor) {
          if(a_hasNeighbor(a_neighborId(childId, i, solver), i, solver)) {
            if(a_hasChildren(a_neighborId(a_neighborId(childId, i, solver), i, solver)) > 0) {
              return false;
            }
          }
        }

        // Change compiler flag to if statement
        // #ifdef RESTRICT_DIAGONAL_LEVEL_JUMPS
        // 1.3. switch to restrict diagonal level jumps
        if(restrictDiagonalLevelJumps) {
          for(MInt j = 0; j < m_noDirs; j++) {
            if(j / 2 == i / 2) continue;
            if(a_hasNeighbor(a_neighborId(childId, i, solver), j, solver) > 0) {
              if(a_hasChildren(a_neighborId(a_neighborId(childId, i, solver), j, solver), solver)) {
                return false;
              }

              // Extend to second diagonal neighbor of child to prevent double interpolation in LB case
              if(extendCheckTo2ndNeighbor) {
                MInt diagNeighbor = a_neighborId(a_neighborId(childId, i, solver), j, solver);
                if(a_hasNeighbor(diagNeighbor, i, solver) > 0) {
                  // Test check for children along path ??
                  if(a_hasChildren(a_neighborId(diagNeighbor, i, solver), solver)) return false;

                  if(a_hasNeighbor(a_neighborId(diagNeighbor, i, solver), j, solver) > 0) {
                    if(a_hasChildren(a_neighborId(a_neighborId(diagNeighbor, i, solver), j, solver), solver)) {
                      return false;
                    }
                  }
                }
              }

              IF_CONSTEXPR(nDim == 3) {
                for(MInt k = 0; k < m_noDirs; k++) {
                  // if ( k/2 == i/2 || k/2 == i/2) continue; // TODO labels:GRID old version, was this a typo
                  // and meant something else?
                  if(k / 2 == i / 2) continue;
                  if(a_hasNeighbor(a_neighborId(a_neighborId(childId, i, solver), j, solver), k, solver) > 0) {
                    if(a_hasChildren(a_neighborId(a_neighborId(a_neighborId(childId, i, solver), j, solver), k, solver),
                                     solver)
                       > 0) {
                      return false;
                    }

                    // Extend to second corner neighbor of child to prevent double interplation in LB case
                    if(extendCheckTo2ndNeighbor) {
                      MInt cornerNeighbor =
                          a_neighborId(a_neighborId(a_neighborId(childId, i, solver), j, solver), k, solver);
                      if(a_hasNeighbor(cornerNeighbor, i, solver) > 0) {
                        // Test check for children along path ??
                        if(a_hasChildren(a_neighborId(cornerNeighbor, i, solver), solver)) return false;

                        if(a_hasNeighbor(a_neighborId(cornerNeighbor, i, solver), j, solver) > 0) {
                          // Test check for children along path ??
                          if(a_hasChildren(a_neighborId(a_neighborId(cornerNeighbor, i, solver), j, solver), solver))
                            return false;

                          if(a_hasNeighbor(a_neighborId(a_neighborId(cornerNeighbor, i, solver), j, solver), k, solver)
                             > 0) {
                            if(a_hasChildren(a_neighborId(
                                   a_neighborId(a_neighborId(cornerNeighbor, i, solver), j, solver), k, solver)))
                              return false;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        //#endif
      }
    }
  }

  // Not needed for lb adaptation
  if(!m_lbGridChecks) {
    // 2. additionally checks if the current cell has a missing neighbor. In this case and if it has children, it has to
    // be a boundary cell and
    //    it's children need to be marked as non-coarsenable cells (this is past this function horizon).
    if(!m_allowInterfaceRefinement)
      for(MInt i = 0; i < m_noDirs; i++)
        if(!a_hasNeighbor(cellId, i, solver)) return false;
  }
  return true;
}


// --------------------------------------------------------------------------------------


/**
 * \brief checks if the given cell may be refined
 * \author Lennart Schneiders
 */
template <MInt nDim>
MBool CartesianGrid<nDim>::refineCheck(
    const MInt cellId, const MInt solver,
    const std::vector<std::function<MInt(const MFloat*, const MInt, const MInt)>>& cellOutsideSolver) {
  // Enable grid checks for LB (refinement)
  MBool restrictDiagonalLevelJumps = m_lbGridChecks;
  MBool extendCheckTo2ndNeighbor = m_lbGridChecks;

  ASSERT(m_tree.solver(cellId, solver), "");
  for(MInt dir = 0; dir < m_noDirs; dir++) {
    if(a_hasNeighbor(cellId, dir, solver)) {
      if(restrictDiagonalLevelJumps) {
        MInt d0 = dir / 2;
        MInt d1 = (d0 + 1) % nDim;
        for(MInt p = 0; p < 2; p++) {
          MInt dir1 = 2 * d1 + p;
          if(a_hasNeighbor(a_neighborId(cellId, dir, solver), dir1, solver) == 0) {
            return false;
          }

          // Ensure two diagonal neigbors to avoid double interpolation in LB
          // For readability introduce names of cells
          MInt firstDiagNghbor;
          if(extendCheckTo2ndNeighbor) {
            firstDiagNghbor = a_neighborId(a_neighborId(cellId, dir, solver), dir1, solver);
            if(a_hasNeighbor(firstDiagNghbor, dir, solver)) {
              if(a_hasNeighbor(a_neighborId(firstDiagNghbor, dir, solver), dir1, solver) == 0) {
                return false;
              }
            } else {
              return false;
            }
          }

          // From last merge: else if. But this crashes 3D case. Thus change to if
          IF_CONSTEXPR(nDim == 3) {
            MInt d2 = (d1 + 1) % nDim;
            for(MInt q = 0; q < 2; q++) {
              MInt dir2 = 2 * d2 + q;
              if(a_hasNeighbor(a_neighborId(a_neighborId(cellId, dir, solver), dir1, solver), dir2, solver) == 0) {
                return false;
              }

              // Ensure two diagonal neigbors (corner) to avoid double interpolation in LB
              // The step-wise procedure is necessary to ensure "filling" between neighbors
              MInt firstCornerNeighbor =
                  a_neighborId(a_neighborId(a_neighborId(cellId, dir, solver), dir1, solver), dir2, solver);
              if(extendCheckTo2ndNeighbor) {
                if(a_hasNeighbor(firstCornerNeighbor, dir, solver)) {
                  if(a_hasNeighbor(a_neighborId(firstCornerNeighbor, dir, solver), dir1, solver)) {
                    if(a_hasNeighbor(a_neighborId(a_neighborId(firstCornerNeighbor, dir, solver), dir1, solver), dir2,
                                     solver)
                       == 0) {
                      return false;
                    }
                  } else {
                    return false;
                  }
                } else {
                  return false;
                }
              }
            }
          }
        }
      }
    } else if(m_allowInterfaceRefinement && hasCutOff()) {
      // continue;
      // check if the cell is a cutOff cell
      // meaning that the largest parent doesn't have a neighbor in that direction!

      MInt parentId = a_parentId(cellId, solver);
      if(parentId > -1 && !a_hasNeighbor(parentId, dir, solver)) {
        continue;
      } else if(parentId > -1) {
        return false;
      }

    } else if(m_allowInterfaceRefinement) {
      MFloat coords[nDim];
      for(MInt k = 0; k < nDim; k++) {
        coords[k] = a_coordinate(cellId, k);
      }
      coords[dir / 2] += ((dir % 2 == 0) ? -F1 : F1) * cellLengthAtCell(cellId);
      MInt isOutside = cellOutsideSolver[solver](coords, a_level(cellId), cellId);
      if(isOutside == 0) {
        return false;
      } else if(isOutside == -1) {
        MInt parentId = a_parentId(cellId, solver);
        while(parentId > -1) {
          if(a_hasNeighbor(parentId, dir, solver)) {
            return false;
          }
          parentId = a_parentId(parentId, solver);
        }
      }
    } else {
      return false;
    }
  }
  return true;
}

//-------------------------------------------------------------------------------------


/** \brief Set minimum und maximum cell levels
 *
 *
 * This function must be called after a dummy cell is appended
 * otherwise the loop must be changed to size().
 */
template <MInt nDim>
void CartesianGrid<nDim>::setLevel() {
  TRACE();
  m_maxLevel = -1;

  for(MInt i = 0; i < m_tree.size(); i++) { // m_noRegularCells
    if(a_isToDelete(i)) continue;
    if(a_hasProperty(i, Cell::IsHalo)) continue;
    m_maxLevel = (m_maxLevel < a_level(i)) ? a_level(i) : m_maxLevel;
  }

  MPI_Allreduce(MPI_IN_PLACE, &m_maxLevel, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "m_maxLevel");

  if(m_minLevel < 0 || m_maxLevel < m_minLevel || m_maxLevel > 31 || m_minLevel > m_maxLevel) {
    mTerm(1, AT_, "Inconsistent min/max levels: " + to_string(m_minLevel) + "/" + to_string(m_maxLevel));
  }
}


/** \brief Deletes a cell (without collector fragmentation) and returns new collector size
 *
 * \author ?, adjusted by Claudia Guenther, Nov 2012
 */
template <MInt nDim>
MInt CartesianGrid<nDim>::deleteCell(const MInt id) {
  if(m_tree.size() == 0) {
    mTerm(1, AT_, " Error in CartesianGrid::deleteCell(), collector is empty.");
  }
  if(id >= m_tree.size()) {
    mTerm(1, AT_, " Error in CartesianGrid::deleteCell(), cell out of range.");
  }

  a_resetProperties(id);

  // If the cell to delete is the last cell, return -1, otherwise return the previously largest id
  const MInt lastId = m_tree.size() - 1;
  const MInt returnValue = (lastId <= id) ? -1 : lastId;

  // Erase cell, move last cell to gap, shrink tree
  m_tree.removeAndFill(id);

  return returnValue;
}


//-------------------------------------------------------------------------------------

/** \brief Return the hilbert index for the given coordinates in 2D/3D
 *
 *  This function normalizes the cell coordinates and calls
 *  the hilbertIndex function.
 *
 *  **IMPORTANT**: The `baseLevel` is the normal refinement level - 1 (e.g. if
 *                 your refinement level is 3, base level is 2).
 *  NOTE: this **IMPORTANT** information seems to be outdated!
 */
template <MInt nDim>
MLong CartesianGrid<nDim>::hilbertIndexGeneric(const MFloat* const coords,
                                               const MFloat* const center,
                                               const MFloat length0,
                                               const MInt baseLevel) const {
  // TRACE();
  ASSERT(length0 > 0.0, "length needs to be > 0.");
  // Relate to unit cube
  MFloat x[nDim];
  for(MInt i = 0; i < nDim; i++) {
    x[i] = (coords[i] - center[i] + 1e-12 + F1B2 * length0) / length0;
    ASSERT(x[i] > F0 && x[i] < F1,
           "Normalized coordinate outside of range (0,1): " + std::to_string(x[i])
               + "; length0 = " + std::to_string(length0) + "; baseLevel = " + std::to_string(baseLevel)
               + "; coord = " + std::to_string(coords[i]) + "; center = " + std::to_string(center[i]));
  }
  const MLong hilbertId = maia::grid::hilbert::index<nDim>(&x[0], (MLong)baseLevel);
  ASSERT(hilbertId > -1, "Invalid hilbert id");
  return hilbertId;
}


template <MInt nDim>
void CartesianGrid<nDim>::saveGrid(const MChar* fileName, const std::vector<std::vector<MInt>>& haloCells,
                                   const std::vector<std::vector<MInt>>& windowCells,
                                   const std::vector<std::vector<MInt>>& azimuthalHaloCells,
                                   const std::vector<MInt>& azimuthalUnmappedHaloCells,
                                   const std::vector<std::vector<MInt>>& azimuthalWindowCells, MInt* recalcIdTree) {
  TRACE();
  maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), haloCells,
                                            windowCells, azimuthalHaloCells, azimuthalUnmappedHaloCells,
                                            azimuthalWindowCells, recalcIdTree);
}

template <MInt nDim>
void CartesianGrid<nDim>::saveGrid(const MChar* fileName, MInt* recalcIdTree) {
  TRACE();

  if(m_newMinLevel < 0) {
    maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), m_haloCells,
                                              m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells,
                                              m_azimuthalWindowCells, recalcIdTree);

  } else { // if possible Increase MinLevel when writing the restart-File!
    // if the m_targetGridMinLevel has already been updated, the new grid restartFile
    // is forced without a backup and further check!

    if(m_newMinLevel == m_targetGridMinLevel && noDomains() == 1) {
      maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), m_haloCells,
                                                m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells,
                                                m_azimuthalWindowCells, recalcIdTree);
    } else if(noDomains() > 1) {
      MInt backup = m_newMinLevel;
      m_newMinLevel = -1;
      maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), m_haloCells,
                                                m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells,
                                                m_azimuthalWindowCells, recalcIdTree);
      m_newMinLevel = backup;

    } else {
      // check that all minLevels are refined
      MInt noUnrefinedMinCells = 0;
      for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
        if(a_isHalo(cellId)) continue;
        if(a_level(cellId) < m_newMinLevel) {
          if(a_noChildren(cellId) < 1) noUnrefinedMinCells++;
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &noUnrefinedMinCells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                    "noUnrefinedMinCells");

      // if the minLevel cells are not sufficiently refiened a reguar output is wirtten
      if(noUnrefinedMinCells > 0) {
        MInt backup = m_newMinLevel;
        if(domainId() == 0) {
          cerr << "Mincells not yet sufficiently refined for minLevel increase from " << m_minLevel << " to "
               << m_newMinLevel << "!" << endl;
        }
        m_newMinLevel = -1;
        maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), m_haloCells,
                                                  m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells,
                                                  m_azimuthalWindowCells, recalcIdTree);
        m_newMinLevel = backup;

      } else { // writing the restartGrid with increased minLevel!

        // write the original grid as backup
        MInt backupMinLevel = m_newMinLevel;
        m_newMinLevel = -1;
        m_targetGridMinLevel = m_minLevel;
        std::stringstream s;
        s << "restartGrid_backup_" << globalTimeStep << ParallelIo::fileExt();
        MString newGridFileName = s.str();

        maia::grid::IO<CartesianGrid<nDim>>::save(
            *this, (m_outputDir + newGridFileName).c_str(), static_cast<Collector<void>*>(nullptr), m_haloCells,
            m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells, m_azimuthalWindowCells, recalcIdTree);

        m_newMinLevel = backupMinLevel;
        m_targetGridMinLevel = m_newMinLevel;

        if(domainId() == 0) {
          cerr << "Increasing minLevel from " << m_minLevel << " to " << m_newMinLevel << endl;
        }

        MInt backupUniform = m_maxUniformRefinementLevel;
        if(m_maxUniformRefinementLevel < m_newMinLevel) {
          if(domainId() == 0) {
            cerr << "Increasing maxUniformRefinementLevel from " << m_maxUniformRefinementLevel << " to "
                 << m_newMinLevel << endl;
          }
          m_maxUniformRefinementLevel = m_newMinLevel;
        }
        // writing the new grid
        maia::grid::IO<CartesianGrid<nDim>>::save(*this, fileName, static_cast<Collector<void>*>(nullptr), m_haloCells,
                                                  m_windowCells, m_azimuthalHaloCells, m_azimuthalUnmappedHaloCells,
                                                  m_azimuthalWindowCells, recalcIdTree);

        m_maxUniformRefinementLevel = backupUniform;
      }
    }
  }
}

//-----------------------------------------------------------------------


/** \brief Load a grid file writen with saveGridDomainPar
    \author Jerry Grimmen
    \date 06.2012
  */

template <MInt nDim>
void CartesianGrid<nDim>::loadGrid(const MString& fileName) {
  TRACE();
  maia::grid::IO<CartesianGrid<nDim>>::load(*this, fileName);
}


//-----------------------------------------------------------------------


/** \brief Global to local id mapping
    \author Lennart Schneiders
  */
template <MInt nDim>
inline MInt CartesianGrid<nDim>::globalIdToLocalId(const MLong& globalId, const MBool termIfNotExisting) {
  ASSERT(m_globalToLocalId.size() > 0, "Error: m_globalToLocalId is empty.");
  auto it = m_globalToLocalId.find(globalId);
  if(it != m_globalToLocalId.end()) {
    return it->second;
  } else {
    TERMM_IF_COND(termIfNotExisting, "Mapping from global to local id not found.");
    return -1;
  }
}


//-----------------------------------------------------------------------


/** \brief Reset cell to default values
    \author Lennart Schneiders
  */
template <MInt nDim>
inline void CartesianGrid<nDim>::resetCell(const MInt& cellId) {
  a_resetProperties(cellId);
  m_tree.resetSolver(cellId);
  m_tree.resetIsLeafCell(cellId);
  a_parentId(cellId) = -1;
  a_globalId(cellId) = -1;
  a_level(cellId) = -1;
  a_noOffsprings(cellId) = 0;
  a_weight(cellId) = F1;
  a_workload(cellId) = F0;
  fill(&a_childId(cellId, 0), &a_childId(cellId, 0) + m_maxNoChilds, -1);
  fill(&a_neighborId(cellId, 0), &a_neighborId(cellId, 0) + m_maxNoNghbrs, -1);
}


// -----------------------------------------------------------------------------


/// \brief Create grid map and partition file if requested by user.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-07
///
/// If donorGridFileName is set and createGridMap is true, create grid map. If
/// grid map already exists, regenerate unless forceGridMapGeneration is false
/// (default: true). If saveDonorGridPartition is true (default: false), also
/// create and store donor grid partition by calling saveDonorGridPartition.
template <MInt nDim>
void CartesianGrid<nDim>::initGridMap() {
  TRACE();

  // Return fast if no donor grid is specified
  if(!Context::propertyExists("donorGridFileName")) {
    return;
  }

  /*! \page propertyPage1
  \section createGridMap
  default = None \n \n
  Create grid map. \n
  Possible values are:
  <ul>
  <li>0 Deactivated </li>
  <li>1 Activated </li>
  </ul>
  Keywords: <i>ALL</i>
  */
  // Return if no grid map should be created
  if(!Context::getBasicProperty<MBool>("createGridMap", AT_)) {
    return;
  }

  // Get name of donor grid file
  const MString donorGridFileName = Context::getBasicProperty<MString>("donorGridFileName", AT_);

  // Get grid map file name if specified
  MString gridMapFileName = m_outputDir + "gridmap.Netcdf";
  gridMapFileName = Context::getBasicProperty<MString>("gridMapFileName", AT_, &gridMapFileName);

  /*! \page propertyPage1
  \section forceGridMapGeneration
  default = 1 \n \n
  Forces a regeneration of grid map, even if a grid map already exists. \n
  Possible values are:
  <ul>
  <li>0 Deactivated </li>
  <li>1 Activated </li>
  </ul>
  Keywords: <i>ALL</i>
  */
  // If grid map should not be reused, re-generate it
  MBool regenerate = true;
  regenerate = Context::getBasicProperty<MBool>("forceGridMapGeneration", AT_, &regenerate);
  if(regenerate || !ParallelIo::fileExists(gridMapFileName, mpiComm())) {
    createGridMap(donorGridFileName, gridMapFileName);
  }

  /*! \page propertyPage1
     \section saveDonorGridPartition
     <code>MBool CartesianGrid::savePartition</code>\n
     default = <code>false</code>\n \n
     Triggers if the donor grid partition is created and stored.
    <ul>
    <li>0 (false-deactivated)</li>
    <li>1 (true-activated)</li>
    </ul>
     Keywords: <i>OUTPUT, </i>
  */
  // Create donor grid partition file if desired
  MBool savePartition = false;
  savePartition = Context::getBasicProperty<MBool>("saveDonorGridPartition", AT_, &savePartition);
  if(savePartition) {
    // Get partition file name if specified
    MString gridPartitionFileName = Context::getBasicProperty<MString>("gridPartitionFileName", AT_);

    // Load relevant information from grid map and store to grid partition file
    saveDonorGridPartition(gridMapFileName, gridPartitionFileName);
  }
}


/// \brief Create file that contains mapping from donor to target cell.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-16
///
/// \param[in] donorGridFileName Name of the original (donor) grid file.
/// \param[in] gridMapFileName Name of the file where grid map should be stored.
///
/// The grid map file contains an array 'cellIds' with a value for each cell of
/// the loaded grid file. The value is either -1 if no corresponding cell exists
/// in the donor grid, or the corresponding cell id (global id). Additionally,
/// there is an array 'noOffspring' that contains the number of mapped donor
/// cell offspring for each target cell. This value is 0 if either a target cell
/// has zero or only one donor cell. Otherwise, for a target cell without
/// children it contains the number of offspring of the donor cell, i.e. the
/// number of all smaller cells contained in this cell on the donor grid that
/// are mapped to this target cell.
///
/// Note: Currently only one-to-one and multiple-to-one mappings are allowed,
///       multiple-to-one mappings are not supported yet. That is, the donor
///       cells must be smaller or at the same level as the target cells, while
///       coarser donor cells are not allowed yet.
template <MInt nDim>
void CartesianGrid<nDim>::createGridMap(const MString& donorGridFileName, const MString& gridMapFileName) {
  TRACE();

  m_log << "Creating grid map from donor grid " << donorGridFileName << " to target grid " << m_gridInputFileName
        << " and writing the results to " << gridMapFileName << "..." << endl;

  //////////////////////////////////////////////////////////////////////////////
  // Step 0: sanity checks
  //////////////////////////////////////////////////////////////////////////////

  // Open donor grid for reading
  using namespace parallel_io;
  ParallelIo donorGrid(donorGridFileName, PIO_READ, mpiComm());

  // Check extents
  array<MFloat, nDim> donorCenter;
  // donorGrid.setOffset(nDim, 0);
  // donorGrid.readArray(&donorCenter[0], "centerOfGravity");
  donorGrid.getAttribute(&donorCenter[0], "centerOfGravity", nDim);
  MFloat donorLengthLevel0 = NAN;
  // donorGrid.readScalar(&donorLengthLevel0, "lengthLevel0");
  donorGrid.getAttribute(&donorLengthLevel0, "lengthLevel0");
  const array<MString, 3> dirs = {{"x", "y", "z"}};
  for(MInt dir = 0; dir < nDim; dir++) {
    const array<MFloat, 2> donorExtent = {
        {donorCenter[dir] - 0.5 * donorLengthLevel0, donorCenter[dir] + 0.5 * donorLengthLevel0}};
    const array<MFloat, 2> targetExtent = {
        {m_centerOfGravity[dir] - 0.5 * lengthLevel0(), m_centerOfGravity[dir] + 0.5 * lengthLevel0()}};
    if(donorExtent[0] < targetExtent[0]) {
      TERMM(1, "Donor grid extents exceed target grid in negative " + dirs[dir] + "-direction");
    }
    if(donorExtent[1] > targetExtent[1]) {
      TERMM(1, "Donor grid extents exceed target grid in positive " + dirs[dir] + "-direction");
    }
  }

  // Define epsilon for floating point comparisons (argh!)... the original
  // definition is taken from the constructor but it must be clear to anyone
  // that this is a less-than-optimal solution
  const MFloat eps = 1.0 / FPOW2(30) * m_lengthLevel0;

  // Check cell length at min level
  MInt donorMinLevel = 0;
  // donorGrid.readScalar(&donorMinLevel, "minLevel");
  donorGrid.getAttribute(&donorMinLevel, "minLevel");
  const MFloat donorMinLevelLength = donorLengthLevel0 * FFPOW2(donorMinLevel);
  const MFloat targetMinLevelLength = cellLengthAtLevel(minLevel());
  if(fabs(donorMinLevelLength - targetMinLevelLength) > eps) {
    TERMM(1,
          "Length of min level cells do not match between donor and target "
          "grid: donor: "
              + to_string(donorMinLevelLength) + "; target: " + to_string(targetMinLevelLength));
  }

  // Check if grid centers are displaced by an integer multiple of the min
  // level length
  for(MInt dir = 0; dir < nDim; dir++) {
    const MFloat displacement = fabs(donorCenter[dir] - m_centerOfGravity[dir]);
    const MFloat quotient = displacement / targetMinLevelLength;
    if(!isApproxInt(quotient, eps)) {
      TERMM(1, "The grid centers are displaced in the " + dirs[dir]
                   + "-direction by a non-integer multiple of the length of a "
                     "partition cell: "
                   + to_string(quotient));
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Step 1: partition donor grid
  //////////////////////////////////////////////////////////////////////////////

  // Read partition cells and partition donor grid
  // const MInt noPartitionCells = donorGrid.getArraySize("partitionCellsId");
  MInt noPartitionCells = 0;
  donorGrid.getAttribute(&noPartitionCells, "noPartitionCells");

  // Read data only on MPI root by setting the offsets accordingly
  donorGrid.setOffset(isMpiRoot() ? noPartitionCells : 0, 0);

  // Ensure that there are no partition level shifts as this is not yet supported
  {
    // MIntScratchSpace partitionCellsLvlDiff(noPartitionCells, AT_, "partitionCellsLvlDiff");
    // donorGrid.readArray(&partitionCellsLvlDiff[0], "partitionCellsLvlDiff");
    /// const MBool hasDiff = any_of(partitionCellsLvlDiff.cbegin(),
    //                               partitionCellsLvlDiff.cend(),
    //                               [](MInt diff) { return diff != 0; });
    MInt hasDiff = 0;
    donorGrid.getAttribute(&hasDiff, "maxPartitionLevelShift");
    if(hasDiff != 0) {
      TERMM(1, "partition level shifts not supported but level difference found");
    }
  }

  // Determine offsets
  MIntScratchSpace partitionCellsId(noPartitionCells, AT_, "partitionCellsId");
  MIntScratchSpace globalIdOffsets(noDomains(), AT_, "globalIdOffsets");
  {
    // Determine partition cell offsets
    MIntScratchSpace offsets(noDomains() + 1, AT_, "offsets");
    MFloatScratchSpace partitionCellsWorkLoad(noPartitionCells, AT_, "partitionCellsWorkLoad");
    // donorGrid.readArray(&partitionCellsWorkLoad[0], "partitionCellsWorkLoad");
    donorGrid.readArray(&partitionCellsWorkLoad[0], "partitionCellsWorkload");
    if(isMpiRoot() && noDomains() > 1) {
      grid::optimalPartitioningSerial(&partitionCellsWorkLoad[0], noPartitionCells, noDomains(), &offsets[0]);
    }

    // Determine global id offsets
    // MIntScratchSpace partitionCellsNoOffsprings(noPartitionCells, AT_,
    //                                       "partitionCellsNoOffsprings");
    // donorGrid.readArray(&partitionCellsNoOffsprings[0], "partitionCellsNoOffsprings");
    // donorGrid.readArray(&partitionCellsId[0], "partitionCellsId");
    donorGrid.readArray(&partitionCellsId[0], "partitionCellsGlobalId");
    if(isMpiRoot() && noDomains() > 1) {
      grid::partitionCellToGlobalOffsets(&offsets[0], &partitionCellsId[0], noDomains(), &globalIdOffsets[0]);
    }
  }

  // Distribute global id offsets
  MPI_Bcast(&globalIdOffsets[0], noDomains(), type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "globalIdOffsets[0]");

  // Distribute partition cells ids
  MPI_Bcast(&partitionCellsId[0], noPartitionCells, type_traits<MInt>::mpiType(), 0, mpiComm(), AT_,
            "partitionCellsId[0]");


  //////////////////////////////////////////////////////////////////////////////
  // Step 2: calculate Hilbert indices for donor grid and current cells
  //////////////////////////////////////////////////////////////////////////////

  // Determine Hilbert index for all partition cells for which coordinates were read
  MIntScratchSpace hilbertIds(noPartitionCells, AT_, "hilbertIds");
  const MInt noCells = m_noInternalCells;
  MIntScratchSpace localHilbertIds(noCells, AT_, "localHilbertIds");
  {
    // Determine offset and length for reading coordinate information
    // const MInt globalCount    = donorGrid.getArraySize("parentId");
    MInt globalCount = 0;
    donorGrid.getAttribute(&globalCount, "noCells");
    const MInt globalIdOffset = globalIdOffsets[domainId()];
    const MInt localCount = (domainId() == noDomains() - 1)
                                ? globalCount - globalIdOffsets[domainId()]
                                : globalIdOffsets[domainId() + 1] - globalIdOffsets[domainId()];

    // Read coordinates
    MFloatScratchSpace coordinates(nDim * localCount, AT_, "coordinates");
    /*donorGrid.setOffset(localCount, globalIdOffset);
    for (MInt i = 0; i < nDim; i++) {
      const MString name = "coordinates_" + to_string(i);
      donorGrid.readArray(&coordinates[i], name, nDim);
    }*/
    MLongScratchSpace minLevelCellsTreeId(noPartitionCells, AT_, "minLevelCellsTreeId");
    {
      donorGrid.setOffset(noPartitionCells, 0);
      donorGrid.readArray(minLevelCellsTreeId.data(), "minLevelCellsTreeId");
    }

    // Determine Hilbert index
    fill(hilbertIds.begin(), hilbertIds.end(), -1);
    for(MInt i = 0; i < noPartitionCells; i++) {
      // Skip partition cells that are outside the range
      const MInt globalId = partitionCellsId[i];
      if(globalId < globalIdOffset || globalId >= globalIdOffset + localCount) {
        continue;
      }

      // Determine coordinates relative to unit cube
      array<MFloat, nDim> x;
      const MInt localCellId = globalId - globalIdOffset;
      maia::grid::hilbert::treeIdToCoordinates<nDim>(&coordinates[localCellId * nDim], minLevelCellsTreeId[i],
                                                     (MLong)m_minLevel, &donorCenter[0], donorLengthLevel0);
      for(MInt j = 0; j < nDim; j++) {
        x[j] = (coordinates[localCellId * nDim + j] - m_centerOfGravity[j] + 0.5 * lengthLevel0()) / lengthLevel0();
      }

      // Calculate Hilbert index
      hilbertIds[i] = maia::grid::hilbert::index<nDim>(&x[0], minLevel());
    }

    // Determine Hilbert index for all cells on current domain
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      // Skip cells that are not partition cells
      if(a_level(cellId) != minLevel()) {
        continue;
      }

      // Determine coordinates relative to unit cube
      array<MFloat, nDim> x;
      for(MInt j = 0; j < nDim; j++) {
        x[j] = (a_coordinate(cellId, j) - m_centerOfGravity[j] + 0.5 * lengthLevel0()) / lengthLevel0();
      }

      // Calculate Hilbert index
      localHilbertIds[cellId] = maia::grid::hilbert::index<nDim>(&x[0], minLevel());
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // Step 3: exchange Hilbert ids of donor grid with all domains (parallel only)
  //////////////////////////////////////////////////////////////////////////////

  // Determine how much data is to exchange and how to receive it
  {
    MIntScratchSpace dataCount(noDomains(), AT_, "dataCount");
    const MInt noCellsFound = count_if(hilbertIds.begin(), hilbertIds.end(), [](const MInt a) { return a != -1; });
    dataCount[domainId()] = noCellsFound;
    MPI_Allgather(MPI_IN_PLACE, 1, maia::type_traits<MInt>::mpiType(), &dataCount[0], 1,
                  maia::type_traits<MInt>::mpiType(), mpiComm(), AT_, "MPI_IN_PLACE", "dataCount[0]");
    MIntScratchSpace displacements(noDomains(), AT_, "displacements");
    displacements[0] = 0;
    for(MInt i = 1; i < noDomains(); i++) {
      displacements[i] = displacements[i - 1] + dataCount[i - 1];
    }
    MIntScratchSpace sendBuffer(noCellsFound, AT_, "sendBuffer");

    // Exchange partition cell ids
    for(MInt count = 0, i = 0; i < noPartitionCells; i++) {
      if(hilbertIds[i] == -1) {
        continue;
      }
      sendBuffer[count++] = partitionCellsId[i];
    }
    MPI_Allgatherv(&sendBuffer[0], noCellsFound, maia::type_traits<MInt>::mpiType(), &partitionCellsId[0],
                   &dataCount[0], &displacements[0], maia::type_traits<MInt>::mpiType(), mpiComm(), AT_,
                   "sendBuffer[0]", "partitionCellsId[0]");

    // Exchange Hilbert ids
    for(MInt count = 0, i = 0; i < noPartitionCells; i++) {
      if(hilbertIds[i] == -1) {
        continue;
      }
      sendBuffer[count++] = hilbertIds[i];
    }
    MPI_Allgatherv(&sendBuffer[0], noCellsFound, maia::type_traits<MInt>::mpiType(), &hilbertIds[0], &dataCount[0],
                   &displacements[0], maia::type_traits<MInt>::mpiType(), mpiComm(), AT_, "sendBuffer[0]",
                   "hilbertIds[0]");
  }


  //////////////////////////////////////////////////////////////////////////////
  // Step 4: determine grid map and write to file
  //////////////////////////////////////////////////////////////////////////////

  // Sort partition cell ids (and Hilbert ids) by Hilbert id
  {
    ScratchSpace<array<MInt, 2>> s(noPartitionCells, AT_, "s");
    for(MInt i = 0; i < noPartitionCells; i++) {
      s[i][0] = hilbertIds[i];
      s[i][1] = partitionCellsId[i];
    }
    sort(s.begin(), s.end(), [](const array<MInt, 2>& a, const array<MInt, 2>& b) { return a[0] < b[0]; });
    for(MInt i = 0; i < noPartitionCells; i++) {
      hilbertIds[i] = s[i][0];
      partitionCellsId[i] = s[i][1];
    }
  }

  MIntScratchSpace gridMap(noCells, AT_, "gridMap");
  fill(gridMap.begin(), gridMap.end(), -1);

  MInt firstMappedPartitionCellId = std::numeric_limits<MInt>::max();
  MInt noMappedPartitionCells = 0;

  // Determine matching partition cells, the first mapped partition cell id and the number
  // of mapped partition cells on this domain
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    // Skip cells that are not partition cells
    if(a_level(cellId) != minLevel()) {
      continue;
    }

    // Find corresponding Hilbert id
    const MInt hilbertId = localHilbertIds[cellId];
    MInt* const lower = lower_bound(hilbertIds.data(), (hilbertIds.data() + hilbertIds.size()), hilbertId);

    // If id was not found, no mapped cell exists; continue with next cell
    if(lower == (hilbertIds.data() + hilbertIds.size()) || *lower != hilbertId) {
      continue;
    }

    // Store corresponding partition cell id to grid map
    const MInt partitionCellId = distance(hilbertIds.data(), lower);
    gridMap[cellId] = partitionCellsId[partitionCellId];

    // Determine first partition cell id
    firstMappedPartitionCellId = min(firstMappedPartitionCellId, partitionCellId);

    // Increase number of mapped min-cells
    noMappedPartitionCells++;
  }

  // First mapped cell id on this domain
  const MInt firstMappedCellId = (noMappedPartitionCells > 0) ? partitionCellsId[firstMappedPartitionCellId] : 0;

  // Determine the total number of donor cells on this domain
  MInt noDonorCells = 0;
  {
    // Read partitionCellsNoOffsprings from donor grid
    MIntScratchSpace partitionCellsNoOffsprings(max(noMappedPartitionCells, 1), AT_, "partitionCellsNoOffsprings");
    MIntScratchSpace partitionCellsId2(noPartitionCells, AT_, "partitionCellsId2");
    // const MInt partitionCellOffset
    //    = (noMappedPartitionCells > 0) ? firstMappedPartitionCellId : 0;
    // donorGrid.setOffset(noMappedPartitionCells, partitionCellOffset);
    // donorGrid.readArray(&partitionCellsNoOffsprings[0], "partitionCellsNoOffsprings");
    donorGrid.setOffset(noPartitionCells, 0);
    donorGrid.readArray(&partitionCellsId2[0], "partitionCellsGlobalId");
    MInt globalCount;
    donorGrid.getAttribute(&globalCount, "noCells");
    for(MInt i = 0; i < noMappedPartitionCells; i++) {
      MInt id = firstMappedPartitionCellId + i;
      MInt nextId = (id == noPartitionCells - 1) ? globalCount : partitionCellsId2[id + 1];
      partitionCellsNoOffsprings[i] = nextId - partitionCellsId2[id];
      ASSERT(partitionCellsNoOffsprings[i] > 0, "");
    }

    // Sum up the number of offsprings of all partition cells on this domain
    for(MInt i = 0; i < noMappedPartitionCells; i++) {
      noDonorCells += partitionCellsNoOffsprings[i];
    }
  }

  // Calculate the donor grid domain offset
  ParallelIo::size_type offset, totalCount;
  ParallelIo::calcOffset(noDonorCells, &offset, &totalCount, mpiComm());

  // Read the noChildIds array of the donor grid
  MIntScratchSpace noChildIds(max(noDonorCells, 1), AT_, "noChildIds");
  donorGrid.setOffset(noDonorCells, offset);
  // donorGrid.readArray(&noChildIds[0], "noChildIds");
  vector<MUchar> cellInfo(noDonorCells);
  donorGrid.readArray(cellInfo.data(), "cellInfo");
  for(MInt i = 0; i < noDonorCells; i++) {
    const MUint childCnt = static_cast<MUint>(cellInfo[i]) & 15u;
    noChildIds[i] = (MInt)childCnt;
  }

  // Storage for the number of donor cell offspring for each target cell
  MIntScratchSpace gridMapNoOffspring(noCells, AT_, "gridMapNoOffspring");
  fill(gridMapNoOffspring.begin(), gridMapNoOffspring.end(), 0);

  // Auxiliary method: Return the number of offspring, i.e. child cells,
  // child-child cells etc., for a given cell based on a list of number of
  // children for each cell, ordered depth-first. For a cell without children, 0
  // is returned.
  auto getNoOffspring = [](const MInt cellId, const MInt* const noChilds) {
    // Initialize with the number of child cells of considered cell
    MInt noOffspring = noChilds[cellId];

    // Loop until all childs of childs of ... are taken into account
    // If there are only leaf cells left 'noOffspring' wont increase any more
    // since the number of children is zero from there on and the loop will end.
    MInt offspringId = 1;
    while(offspringId <= noOffspring) {
      noOffspring += noChilds[cellId + offspringId];
      offspringId++;
    }
    return noOffspring;
  };

  // Determine all matching offspring cells of the mapped partition cells
  MInt noOneToMultipleMappings = 0;
  for(MInt cellId = 0; cellId < noCells;) {
    // Check if this is a partition cell with mapped, corresponding donor cell and
    // continue with the next cell if not. Non-partition cells will never be
    // considered as the else-condition below will take care of them.
    MInt donorGlobalCellId = gridMap[cellId];
    if(donorGlobalCellId == -1) {
      cellId++;
      continue;
    }

    // Store child cells in grid map
    if(a_noChildren(cellId) == 0) {
      // No child cells on target grid, store the number of offspring on the
      // donor grid (0 if the donor cell has no children)
      const MInt donorLocalCellId = donorGlobalCellId - firstMappedCellId;
      gridMapNoOffspring[cellId] = getNoOffspring(donorLocalCellId, &noChildIds[0]);

      // Continue with next cell
      cellId++;
    } else {
      // partition cell on target grid has children
      // Loop as long as the current cell is a descendant of the considered min
      // cell
      const MInt partitionCellLevel = a_level(cellId);
      // 'cellId' belongs to the partition cell, check the next cell in the first
      // iteration
      MInt currentCellId = cellId + 1;

      while(currentCellId < noCells && a_level(currentCellId) > partitionCellLevel) {
        MInt donorLocalCellId = donorGlobalCellId - firstMappedCellId;
        // No child cells on the donor grid
        if(noChildIds[donorLocalCellId] == 0) {
          // Store mapping for matching cell
          gridMap[cellId] = donorGlobalCellId;
          const MInt parentLevel = a_level(cellId);

          // Continue with next cell, this will either be a child of the last
          // cell, i.e. there is a one-to-many mapping, or it has the same or a
          // lower level, i.e. the last cell belongs to a one-to-one mapping.
          cellId++;
          // Loop over all offspring of the current target cell if there are any
          while(cellId < noCells && a_level(cellId) > parentLevel) {
            gridMap[cellId] = donorGlobalCellId; // The same for all offspring
            // Use -1 to indicate 'no additional donor cell' for this cell
            gridMapNoOffspring[cellId] = -1;
            noOneToMultipleMappings++;

            // Continue with next cell
            cellId++;
          }
          // Move on to the next donor cell (one-to-one or one-to-many finished)
          donorGlobalCellId++;
        } else { // Donor cell has children
          if(a_noChildren(cellId) == 0) {
            // Target cell has no children, i.e. many-to-one mapping
            gridMap[cellId] = donorGlobalCellId;

            // Store number of offspring on the donor grid
            donorLocalCellId = donorGlobalCellId - firstMappedCellId;
            const MInt noOffspring = getNoOffspring(donorLocalCellId, &noChildIds[0]);
            gridMapNoOffspring[cellId] = noOffspring;

            // Increase the global donor cell id by the number of mapped cells
            donorGlobalCellId += noOffspring + 1;
          } else { // Both cells have children
            // Store mapping and move to next donor cell (and next target cell)
            gridMap[cellId] = donorGlobalCellId;
            donorGlobalCellId++;
          }
          // Continue with next cell
          cellId++;
        }
        // Set current cell id which to check in while() loop if it is still a
        // descendant of the considered partition cell
        currentCellId = cellId;
      }
    }
  }

  // Open grid map file and write to it
  ParallelIo gridMapFile(gridMapFileName, PIO_REPLACE, mpiComm());
  gridMapFile.setAttribute(donorGridFileName, "donorGridFileName");
  gridMapFile.setAttribute(m_gridInputFileName, "gridFile");
  gridMapFile.defineScalar(PIO_INT, "noCells");
  gridMapFile.defineArray(PIO_INT, "cellIds", m_domainOffsets[noDomains()]);
  gridMapFile.defineArray(PIO_INT, "noOffspring", m_domainOffsets[noDomains()]);
  gridMapFile.setOffset(noCells, m_domainOffsets[domainId()]);
  gridMapFile.writeScalar(m_domainOffsets[noDomains()], "noCells");
  gridMapFile.writeArray(&gridMap[0], "cellIds");
  gridMapFile.writeArray(&gridMapNoOffspring[0], "noOffspring");

  m_log << "Created & saved grid map file." << endl;
  if(noOneToMultipleMappings > 0) {
    m_log << "WARNING: There were " << noOneToMultipleMappings
          << " cells on the target grid that are finer than the "
             "corresponding donor cell."
          << endl;
  }
}


/// \brief Create and store donor grid partition for volume-coupled simulations.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-07
///
/// \param[in] gridMapFileName Name of the grid map that was already generated.
/// \param[in] gridPartitionFileName Name of the partition file to use.
///
/// Store the grid partition information for a given donor grid and for the
/// current number of MPI ranks. Works by finding the global cell id of the
/// first mapped cell in the grid map file.
template <MInt nDim>
void CartesianGrid<nDim>::saveDonorGridPartition(const MString& gridMapFileName, const MString& gridPartitionFileName) {
  TRACE();

  using namespace maia::parallel_io;

  // Open grid map file and read mapping information
  ParallelIo gridMapFile(gridMapFileName, PIO_READ, mpiComm());
  MString donorGridFileName;
  gridMapFile.getAttribute(&donorGridFileName, "donorGridFileName");
  gridMapFile.setOffset(m_noInternalCells, m_domainOffsets[domainId()]);
  MIntScratchSpace gridMap(m_noInternalCells, AT_, "gridMap");
  gridMapFile.readArray(&gridMap[0], "cellIds");

  m_log << "Saving grid partition for donor grid " << donorGridFileName << " and writing the results to "
        << gridPartitionFileName << "..." << endl;

  // Determine global cell id of first mapped cell on this domain
  auto c = find_if(gridMap.begin(), gridMap.end(), [](const MInt a) { return a != -1; });
  const MInt cellId = (c == end(gridMap)) ? -1 : *c;

  // Open grid partition file and write to it
  ParallelIo file(gridPartitionFileName, PIO_REPLACE, mpiComm());
  file.setAttribute(m_gridInputFileName, "targetGridFileName");
  file.setAttribute(donorGridFileName, "gridFile");
  file.defineArray(PIO_INT, "firstCellIds", noDomains());
  file.setOffset(1, domainId());
  file.writeArray(&cellId, "firstCellIds");

  m_log << "Saved grid partition file." << endl;
}


/// \brief Return partition cell offsets based on grid partition file.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-07
///
/// \param[in] gridPartitionFileName Name of grid partition file to use.
/// \param[in] partitionCellsId Point to list of partition cells ids.
/// \param[in] noPartitionCells Number of partition cells.
/// \param[out] offsets Resulting offsets.
///
/// This method is a drop-in replacement for maia::grid::partition(...). It is
/// intended for coupled multi-solver simulations, where both solvers need to
/// exchange volume data and thus corresponding partition cells and their subtrees
/// must end up on the same MPI rank.
template <MInt nDim>
void CartesianGrid<nDim>::loadDonorGridPartition(const MLong* const partitionCellsId, const MInt noPartitionCells) {
  TRACE();

  MLongScratchSpace partitionCellOffsets(noDomains() + 1, AT_, "partitionCellOffsets");

  if(domainId() == 0) {
    MString gridPartitionFileName = m_outputDir + "partition.Netcdf";
    gridPartitionFileName = Context::getBasicProperty<MString>("gridPartitionFileName", AT_, &gridPartitionFileName);

    // Open grid partition file and read out first cell ids for each domain
    using namespace maia::parallel_io;
    ParallelIo file(gridPartitionFileName, PIO_READ, MPI_COMM_SELF);
    const MInt noTargetDomains = file.getArraySize("firstCellIds");
    file.setOffset(noTargetDomains, 0);
    MIntScratchSpace firstCellIds(noTargetDomains, AT_, "firstCellIds");
    file.readArray(&firstCellIds[0], "firstCellIds");

    // Sanity check: make sure that first cell ids are monotonically increasing
    MInt previous = -1;
    for(auto&& first : firstCellIds) {
      // Skip if first id is -1 as it means there are no cells on that domain
      if(first == -1) {
        continue;
      }

      // Abort if cell id is less than previous one
      if(first < previous) {
        TERMM(1, "Cell ids not monotonically increasing: " + to_string(first) + " < " + to_string(previous)
                     + ". Did you forget to set 'targetGridFileName' when "
                       "generating the donor grid file?");
      }

      // Store id for next iteration
      previous = first;
    }

    // Determine partition cell offsets. For each domain, find the position of the first
    // cell id in the list of mincell ids. This position is the desired offset.
    // First, store pointers for convenience
    const MLong* first = partitionCellsId;
    const MLong* const last = first + noPartitionCells;
    MInt donorDomainId = 0;
    for(MInt d = 0; d < noTargetDomains; d++) {
      // Skip target domain if it does not receive donor cells
      const MInt firstCellId = firstCellIds[d];
      if(firstCellId == -1) {
        continue;
      }

      // Search list of partition cells for first cell id and store pointer
      const auto bound = lower_bound(first, last, firstCellId);
      if(distance(bound, last) == 0) {
        // If value was not found, somewhere (here or in the generation of the
        // grid partition file) an error must have occurred
        TERMM(1, "partition cell not found, this should not happen. Go fix your code!");
      } else {
        // If value was found, store the distance between the beginning of the
        // array and the found value as the offset that is searched for
        partitionCellOffsets[donorDomainId] = distance(partitionCellsId, bound);

        // Reduce search space for increased efficiency
        first = bound;

        // Increment domain id
        donorDomainId++;
      }
    }

    partitionCellOffsets[noDomains()] = noPartitionCells;
    m_log << "Cartesian grid loaded grid partition from " << gridPartitionFileName << endl;

    m_domainOffsets[0] = 0;
    for(MInt i = 1; i < noDomains(); i++) {
      m_domainOffsets[i] = partitionCellsId[partitionCellOffsets[i]];
    }
    m_domainOffsets[noDomains()] = m_noCellsGlobal;

    // Sanity check: donorDomainId must be equal to number of domains after last
    // iteration
    if(donorDomainId != noDomains()) {
      TERMM(1, "mismatch between number of donor domains in partition file and "
               "actual number of "
               "donor domains");
    }
  }
  MPI_Bcast(partitionCellOffsets.data(), noDomains() + 1, MPI_LONG, 0, mpiComm(), AT_, "partitionCellOffsets.data()");
  MPI_Bcast(m_domainOffsets, noDomains() + 1, MPI_LONG, 0, mpiComm(), AT_, "m_domainOffsets");
  m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];
  m_localPartitionCellOffsets[0] =
      partitionCellOffsets[domainId()]; // begin of the local partitionCells in the gobal partitionCell array
  m_localPartitionCellOffsets[1] = partitionCellOffsets[domainId() + 1]; // end index of the gobal partitionCell array
  m_localPartitionCellOffsets[2] =
      m_noPartitionCellsGlobal; // end of the local partitionCells in the gobal partitionCell array
  m_noInternalCells = m_domainOffsets[domainId() + 1] - m_domainOffsets[domainId()];
}


/// \brief Create the mapping from global to local cell ids
template <MInt nDim>
void CartesianGrid<nDim>::createGlobalToLocalIdMapping() {
  TRACE();

  std::map<MLong, MInt>().swap(m_globalToLocalId);
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    // multiple periodic cells with identical globalId may appear
    // if neccessary, m_globalToLocalId should be changed to a multimap
    if(a_hasProperty(cellId, Cell::IsPeriodic)) continue;

    if(a_hasProperty(cellId, Cell::IsToDelete)) continue;

    if(a_globalId(cellId) > -1) {
      ASSERT(m_globalToLocalId.count(a_globalId(cellId)) == 0,
             "Global id already in map: " + std::to_string(a_globalId(cellId)) + " " + std::to_string(cellId) + " "
                 + std::to_string(m_globalToLocalId[a_globalId(cellId)]));
      m_globalToLocalId[a_globalId(cellId)] = cellId;
    } else if(!a_isToDelete(cellId)) {
      TERMM(1,
            "Error: invalid global id for cell #" + std::to_string(cellId) + ": " + std::to_string(a_globalId(cellId)));
    }
  }
}


// -----------------------------------------------------------------------------
/** \brief Change global child/neighbor/parent cell ids into local cell ids. Requires that
 *         m_globalToLocalId contains the current mapping from global to local cell ids (see
 *         createGlobalToLocalIdMapping()).
 *
 * \author Stephan Schlimpert
 * \date 9.12.2012
 */
template <MInt nDim>
void CartesianGrid<nDim>::changeGlobalToLocalIds() {
  TRACE();

  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    //  local childIds
    for(MInt childId = 0; childId < m_maxNoChilds; childId++) {
      if(a_childId(cellId, childId) == -1) {
        continue;
      }

      if(m_globalToLocalId.count(a_childId(cellId, childId)) > 0) {
        a_childId(cellId, childId) = m_globalToLocalId[a_childId(cellId, childId)];
      } else {
        a_childId(cellId, childId) = -1;
      }
    }
    // this is commented out, because the noChildIds of the non leaf level halo cells on the second layer would be
    // corrected here to zero such that certain parts of the code will use different neighbor and surface informations
    // especially at the cut boundaries. This leads to serial-parallel inconsistencies because m_noChildIds is used as a
    // trigger for some algorithms
    /*
    a_noChildren( cellId ) = 0;
    for (MInt childId = 0; childId < m_maxNoChilds; childId++) {
      if (a_childId( cellId ,  childId ) > -1) {
        a_noChildren( cellId )++;
      }
    }
    */
    // local neighbor Ids
    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      if(a_neighborId(cellId, dirId) > -1) {
        if(m_globalToLocalId.count(a_neighborId(cellId, dirId)) > 0) {
          a_neighborId(cellId, dirId) = m_globalToLocalId[a_neighborId(cellId, dirId)];
        } else {
          a_neighborId(cellId, dirId) = -1;
        }
      }
    }

    // local parentId
    if(a_parentId(cellId) == -1) {
      a_parentId(cellId) = -1;
    } else {
      if(m_globalToLocalId.count(a_parentId(cellId)) > 0) {
        a_parentId(cellId) = m_globalToLocalId[a_parentId(cellId)];
      } else {
        a_parentId(cellId) = -1;
      }
    }
  }
}


// -----------------------------------------------------------------------------


/// Convert parent ids, neighbor ids, and child ids from local to global cell ids
template <MInt nDim>
void CartesianGrid<nDim>::localToGlobalIds() {
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_parentId(cellId) > -1) {
      a_parentId(cellId) = a_globalId(static_cast<MInt>(a_parentId(cellId)));
    }

    for(MInt i = 0; i < m_noDirs; i++) {
      if(a_hasNeighbor(cellId, i)) {
        a_neighborId(cellId, i) = a_globalId(static_cast<MInt>(a_neighborId(cellId, i)));
      }
    }

    for(MInt i = 0; i < IPOW2(nDim); i++) {
      if(a_childId(cellId, i) > -1) {
        a_childId(cellId, i) = a_globalId(static_cast<MInt>(a_childId(cellId, i)));
      }
    }
  }
}


/**
 * \brief Update number of internal cells, recalculate domain offsets, recalculate global ids for all
 *  internal cells and update global ids of all halo cells.
 *
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::computeGlobalIds() {
  // Update list of min-level cell ids (required later)
  storeMinLevelCells();

  // Update number of internal cells
  m_noInternalCells = 0;
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    if(a_isToDelete(cellId)) continue;
    m_noInternalCells++;
  }

  // Exchange number of internal cells per domain
  const MInt noCells = m_noInternalCells;
  MIntScratchSpace noCellsPerDomain(noDomains(), AT_, "noCellsPerDomain");
  MPI_Allgather(&noCells, 1, type_traits<MInt>::mpiType(), &noCellsPerDomain[0], 1, type_traits<MInt>::mpiType(),
                mpiComm(), AT_, "noCells", "noCellsPerDomain[0]");
  ASSERT(noCellsPerDomain[domainId()] == noCells, "Local number of cells does not match.");

  if(m_domainOffsets == nullptr) {
    mAlloc(m_domainOffsets, noDomains() + 1, "m_domainOffsets", 0L, AT_);
  }

  // Re-calculate domain offsets
  m_domainOffsets[0] = m_32BitOffset;
  for(MInt d = 0; d < noDomains(); d++) {
    m_domainOffsets[d + 1] = m_domainOffsets[d] + static_cast<MLong>(noCellsPerDomain[d]);
  }

  m_noCellsGlobal = m_domainOffsets[noDomains()] - m_32BitOffset;

  // Note: requirement is that min-level cells are stored, done above!
  const MInt firstMinLvlCell = m_minLevelCells[0];
  MInt localCnt = 0;

  // Check the first min-level cell: if it is a partition-level ancestor and a halo it is the root
  // cell of the first local partition cell. The subtree starting from that min-level halo cell with
  // all locally relevant partition-level ancestor cells is contained in the local tree,
  // descendStoreGlobalId will skip any halo cell and will assign the first local cell the global id
  // corresponding to the domain offset and then continue with all offspring cells.
  if(a_hasProperty(firstMinLvlCell, Cell::IsPartLvlAncestor) && a_hasProperty(firstMinLvlCell, Cell::IsHalo)) {
    descendStoreGlobalId(firstMinLvlCell, localCnt);
  }

  // Re-calculate global ids for all internal cells starting from the min-level cells
  // Iteration starts at 1 if first min-lvl cell is a halo partLvlAncestor (see above)
  for(MUint i = std::max(0, std::min(1, localCnt)); i < m_minLevelCells.size(); i++) {
    const MInt cellId = m_minLevelCells[i];

    // skip halos, the only relevant halo min-level cell in case of a partition level shift is
    // handled above
    if(a_hasProperty(cellId, Cell::IsHalo)) continue;
    descendStoreGlobalId(cellId, localCnt);
  }

  // Update global ids for all halo cells
  if(noNeighborDomains() > 0) {
    ScratchSpace<MLong> globalId(m_tree.size(), AT_, "globalId");
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
        MInt cellId = m_windowCells[i][j];
        globalId(cellId) = a_globalId(cellId);
        ASSERT(globalId(cellId) >= m_domainOffsets[domainId()] && globalId(cellId) < m_domainOffsets[domainId() + 1],
               "");
      }
    }
    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &globalId[0], 1);

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
        MInt cellId = m_haloCells[i][j];
        a_globalId(cellId) = globalId(cellId);
      }
    }
  }

  if(m_azimuthalPer && noAzimuthalNeighborDomains() > 0) {
    ScratchSpace<MLong> globalId(m_tree.size(), AT_, "globalId");
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalWindowCells[i].size(); j++) {
        MInt cellId = m_azimuthalWindowCells[i][j];
        globalId(cellId) = a_globalId(cellId);
        ASSERT(globalId(cellId) >= m_domainOffsets[domainId()] && globalId(cellId) < m_domainOffsets[domainId() + 1],
               "");
      }
    }
    maia::mpi::exchangeData(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                            &globalId[0], 1);
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        MInt cellId = m_azimuthalHaloCells[i][j];
        a_globalId(cellId) = globalId(cellId);
      }
    }
  }

  // Update the global to local id mapping
  createGlobalToLocalIdMapping();

  // Note: when updating the partition cells it can be possible that a rank currently does not have
  // a partition cell anymore
  TERMM_IF_NOT_COND(m_noPartitionCells > 0 || m_updatingPartitionCells,
                    "Error: number of partition cells needs to be at least 1.");

  updatePartitionCellInformation();
}


/// \brief Update the partition cell local/global ids and the partition cell global offsets
template <MInt nDim>
void CartesianGrid<nDim>::updatePartitionCellInformation() {
  // update partition cell global ids
  MInt noLocalPartitionCells = 0;
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_isToDelete(cellId)) {
      continue;
    }
    if(a_hasProperty(cellId, Cell::IsHalo)) {
      continue;
    }

    if(a_hasProperty(cellId, Cell::IsPartitionCell)) {
      m_localPartitionCellGlobalIds[noLocalPartitionCells] = a_globalId(cellId);
      noLocalPartitionCells++;
    }
  }

  if(noLocalPartitionCells == 0 && !m_updatingPartitionCells) {
    TERMM(1, "noLocalPartitionCells == 0");
    return;
  }

  if(noLocalPartitionCells != m_noPartitionCells) {
    TERMM(-1, "Mismatch in partitionCells " + to_string(m_noPartitionCells) + "/" + to_string(noLocalPartitionCells));
  }

  // sort by globalIds
  sort(m_localPartitionCellGlobalIds, m_localPartitionCellGlobalIds + noLocalPartitionCells);

  // determine offset within GLOBAL!!! partition cells
  MIntScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
  MPI_Allgather(&noLocalPartitionCells, 1, MPI_INT, &localPartitionCellCounts[0], 1, MPI_INT, mpiComm(), AT_,
                "noLocalPartitionCells", "localPartitionCellCounts[0]");

  // sum up all previous domains partition cells
  MInt offset = 0;
  for(MInt dId = 0; dId < domainId(); dId++) {
    offset += localPartitionCellCounts[dId];
  }

  // Set local partition cell offset, the next offset and the number of global partition cells
  m_localPartitionCellOffsets[0] = offset;
  m_localPartitionCellOffsets[1] = m_localPartitionCellOffsets[0] + noLocalPartitionCells;
  m_localPartitionCellOffsets[2] = m_noPartitionCellsGlobal;

  // Create list of partition cell local ids (also ordered by global id)
  for(MInt i = 0; i < m_noPartitionCells; i++) {
    m_localPartitionCellLocalIds[i] = globalIdToLocalId(m_localPartitionCellGlobalIds[i]);
  }
}


/**
 * \brief Recursively descend subtree to reset global id for internal cells
 *
 * Note: This algorithm only makes sense if called *in-order* for all local min cells /// or ///
 *       in case of a partition level shift starting with the min-level halo partition level
 *       ancestor whos offspring is the first local partition cell and the called *in-order* for all
 *       local min-level cells
 *
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::descendStoreGlobalId(const MInt cellId, MInt& localCnt) {
  // Update global id (unless it is a halo cell)
  if(!a_hasProperty(cellId, Cell::IsHalo)) {
    ASSERT(localCnt < m_noInternalCells, "");
    a_globalId(cellId) = m_domainOffsets[domainId()] + localCnt++;
    ASSERT(a_globalId(cellId) >= m_domainOffsets[domainId()] && a_globalId(cellId) < m_domainOffsets[domainId() + 1],
           "Error: global id outside of global id range for this domain.");
  }

  // Descend tree to all children
  for(MInt child = 0; child < ipow(2, nDim); child++) {
    if(a_childId(cellId, child) < 0) {
      continue;
    }
    descendStoreGlobalId(a_childId(cellId, child), localCnt);
  }
}


/// Recursively update noOffsprings and workLoad for all descendants of `cellId`
/// NOTE: this will not work properly for partition level ancestors in case of partition level
/// shifts, use calculateNoOffspringsAndWorkload() instead
template <MInt nDim>
void CartesianGrid<nDim>::descendNoOffsprings(const MLong cellId, const MLong offset) {
  ASSERT(!a_hasProperty(cellId, Cell::IsPartLvlAncestor), "not supposed to be called for a partition level ancestor");
  a_noOffsprings(cellId) = 1;
  a_workload(cellId) = a_weight(cellId);
  if(a_noChildren(cellId) > 0) {
    for(MInt child = 0; child < ipow(2, nDim); child++) {
      if(a_childId(cellId, child) < 0) continue;
      MLong childId = a_childId(cellId, child) - offset;
      descendNoOffsprings(childId, offset);
      a_noOffsprings(cellId) += a_noOffsprings(childId);
      a_workload(cellId) += a_workload(childId);
    }
  }
}


/**
 * \brief Store cell ids of all min-level cells.
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::storeMinLevelCells(const MBool updateMinlevel) {
  // Create map from hilbert id to all internal min cells
  m_minLevelCells.clear();
  map<MLong, MInt> minLevelCells;
  const MInt minLevel = updateMinlevel ? m_newMinLevel : m_minLevel;
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) continue; // min-level halo partition level ancestor handled below
    if(a_isToDelete(cellId)) continue;
    if(a_level(cellId) == minLevel) {
      const MLong hilbertId =
          updateMinlevel ? generateHilbertIndex(cellId, m_newMinLevel) : generateHilbertIndex(cellId);
      TERMM_IF_COND(minLevelCells.count(hilbertId), "Error: duplicate hilbertId.");
      minLevelCells.insert(pair<MLong, MInt>(hilbertId, cellId));
    }
  }

  MBool found = false;
  MInt minLvlHaloPartLvlAncestor = -1;
  for(auto& i : m_partitionLevelAncestorIds) {
    TERMM_IF_NOT_COND(a_hasProperty(i, Cell::IsPartLvlAncestor),
                      "cell is not a partition level ancestor: " + std::to_string(i));
    TERMM_IF_COND(a_isToDelete(i), "Error: partition level ancestor marked for deletion.");
    // If there is a min-level halo partition level ancestor this is the min-level cell for the
    // first partition cell on this domain, it will be stored in m_minLevelCells as first entry
    // (lowest hilbert id)!
    if(a_level(i) == minLevel && a_hasProperty(i, Cell::IsHalo)) {
      TERMM_IF_COND(found, "there should only be a single min-level halo partition level ancestor "
                           "in the list of partitionLevelAncestorIds! Other halo cells with "
                           "these properties shouldnt be in that list!");
      const MInt hilbertId = updateMinlevel ? generateHilbertIndex(i, m_newMinLevel) : generateHilbertIndex(i);
      TERMM_IF_COND(minLevelCells.count(hilbertId), "duplicate hilbertId.");
      minLevelCells.insert(pair<MInt, MInt>(hilbertId, i));
      minLvlHaloPartLvlAncestor = i;
      found = true;

      // DEBUG
      /* m_log << "found a min-level halo partition level ancestor: " << i << " " << hilbertId */
      /*         << std::endl; */
    }
  }

  for(auto& minLevelCell : minLevelCells) {
    m_minLevelCells.push_back(minLevelCell.second);
  }

  // add halo min level cells to the back
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    // Skip the min-level halo partition level ancestor which was already added as a first entry to
    // m_minLevelCells above
    if(found && cellId == minLvlHaloPartLvlAncestor) continue;

    if(a_level(cellId) == minLevel && a_hasProperty(cellId, Cell::IsHalo)) {
      TERMM_IF_COND(a_isToDelete(cellId), "Error: min-level halo marked for deletion.");
      m_minLevelCells.push_back(cellId);
    }
  }
}


/**
 * \brief Compute distance to leaf level
 * \author Lennart Schneiders
 */
template <MInt nDim>
void CartesianGrid<nDim>::computeLeafLevel() {
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    m_tree.resetIsLeafCell(cellId);
  }
  for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
    for(MInt solverId = 0; solverId < m_tree.noSolvers(); solverId++) {
      if(m_tree.solver(cellId, solverId)) {
        a_isLeafCell(cellId, solverId) =
            !a_hasChildren(cellId, solverId) && !a_hasProperty(cellId, Cell::IsPartLvlAncestor);
      }
    }
  }

  if(noNeighborDomains() > 0) {
    maia::mpi::exchangeBitset(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), &m_tree.leafCellBits(0),
                              m_tree.size());
  }

  if(m_azimuthalPer) {
    if(noAzimuthalNeighborDomains() > 0) {
      maia::mpi::exchangeBitset(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                                &m_tree.leafCellBits(0), m_tree.size());
    }

    // Since connection between two different cell levels is possible in
    // azimuthal periodicity. Azimuthal halos need to be adjusted.
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        MInt cellId = azimuthalHaloCell(i, j);
        for(MInt solverId = 0; solverId < m_tree.noSolvers(); solverId++) {
          if(m_tree.solver(cellId, solverId)) {
            a_isLeafCell(cellId, solverId) =
                !a_hasChildren(cellId, solverId) && !a_hasProperty(cellId, Cell::IsPartLvlAncestor);
          }
        }
      }
    }
    for(MInt i = 0; i < noAzimuthalUnmappedHaloCells(); i++) {
      MInt cellId = azimuthalUnmappedHaloCell(i);
      for(MInt solverId = 0; solverId < m_tree.noSolvers(); solverId++) {
        if(m_tree.solver(cellId, solverId)) {
          a_isLeafCell(cellId, solverId) =
              !a_hasChildren(cellId, solverId) && !a_hasProperty(cellId, Cell::IsPartLvlAncestor);
        }
      }
    }
  }
}


/**
 * \brief Balance the grid according to the given cell send/recv counts.
 * \author Jerry Grimmen, Lennart Schneiders, Ansgar Niemoeller
 */
template <MInt nDim>
void CartesianGrid<nDim>::balance(const MInt* const noCellsToReceiveByDomain,
                                  const MInt* const noCellsToSendByDomain,
                                  const MInt* const sortedCellId,
                                  const MLong* const partitionCellOffsets,
                                  const MLong* const globalIdOffsets) {
  TRACE();

  const MInt noCells = m_tree.size();
  const MInt receiveSize = noCellsToReceiveByDomain[noDomains()];

  // Note: cellInfo needs to be communicated since it cannot be fully reconstructed from the
  // communicated data but it is required to communicate the partition level ancestor data and
  // setup the missing subtrees of the grid
  ScratchSpace<MUchar> cellInfo(receiveSize, AT_, "cellInfo");
  {
    // Temporary buffers
    ScratchSpace<MUint> cellInfoSend(noCells, AT_, "cellInfoSend");
    ScratchSpace<MUint> cellInfoRecv(receiveSize, AT_, "cellInfoRecv");

    // TODO labels:GRID,DLB since the cell info needs to be communicated maybe some of the other grid information is
    // not required anymore if the tree is rebuild with the cell info?
    for(MInt i = 0; i < noCells; i++) {
      // TODO labels:GRID,DLB make cellInfo computation a function that can be reused?
      const MUint noChilds = (MUint)a_noChildren(i);
      const MUint isMinLevel = (a_level(i) == m_minLevel);
      const MInt parent = m_globalToLocalId[a_parentId(i)];
      MUint position = 0;
      if(parent > -1) {
        for(MUint j = 0; j < (unsigned)m_maxNoChilds; j++) {
          if(a_childId(parent, j) == a_globalId(i)) {
            position = j;
          }
        }
      }
      const MUint tmpBit = noChilds | (position << 4) | (isMinLevel << 7);
      cellInfoSend[i] = tmpBit;
    }

    // Redistribute cell info
    maia::mpi::communicateData(&cellInfoSend[0], noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                               noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &cellInfoRecv[0]);

    for(MInt i = 0; i < receiveSize; i++) {
      // Convert cell info
      cellInfo[i] = static_cast<MUchar>(cellInfoRecv[i]);
    }
  }


  // ScratchSpaces to hold cell datas
  MLongScratchSpace parentId(receiveSize, FUN_, "parentId");
  MLongScratchSpace childIds(receiveSize, IPOW2(nDim), FUN_, "childIds");
  MLongScratchSpace nghbrIds(receiveSize, m_noDirs, FUN_, "nghbrIds");
  MLongScratchSpace globalId(receiveSize, FUN_, "globalId");
  MIntScratchSpace level(receiveSize, FUN_, "level");

  MFloatScratchSpace coordinates(receiveSize, nDim, FUN_, "coordinates");
  MFloatScratchSpace weight(receiveSize, FUN_, "weight");

  ScratchSpace<maia::grid::tree::SolverBitsetType> solver(receiveSize, FUN_, "solver");
  ScratchSpace<maia::grid::tree::PropertyBitsetType> properties(receiveSize, AT_, "properties");

  // Reset ScratchSpaces
  parentId.fill(-1);
  childIds.fill(-1);
  nghbrIds.fill(-1);
  globalId.fill(-1);
  level.fill(-1);

  coordinates.fill(-1.0);
  weight.fill(1.0);

  solver.fill(maia::grid::tree::SolverBitsetType());
  properties.fill(maia::grid::tree::PropertyBitsetType());

  // Communicate Cell Data
  maia::mpi::communicateData(&a_parentId(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &parentId[0]);
  maia::mpi::communicateData(&a_childId(0, 0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, IPOW2(nDim), &childIds[0]);
  maia::mpi::communicateData(&a_neighborId(0, 0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, m_noDirs, &nghbrIds[0]);
  maia::mpi::communicateData(&a_level(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &level[0]);
  maia::mpi::communicateData(&a_globalId(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &globalId[0]);
  maia::mpi::communicateData(&a_coordinate(0, 0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, nDim, &coordinates[0]);
  maia::mpi::communicateData(&a_weight(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                             noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &weight[0]);
  maia::mpi::communicateBitsetData(&m_tree.solverBits(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                                   noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &solver[0]);

  maia::mpi::communicateBitsetData(&m_tree.properties(0), noCells, sortedCellId, noDomains(), domainId(), mpiComm(),
                                   noCellsToSendByDomain, noCellsToReceiveByDomain, 1, &properties[0]);

  std::map<MLong, MInt>().swap(m_globalToLocalId);

  // Reset tree
  m_tree.clear();

  // Storage for min-level cell information
  vector<MInt> localMinLevelCells(0);
  vector<MLong> minLevelCellsTreeId(0);
  ScratchSpace<MInt> localMinLevelId(receiveSize, AT_, "localMinLevelId");
  localMinLevelId.fill(-1);

  vector<MLong> partitionCells(0);

  // Iterate over all received cells and add each cell individually (cells are sorted by global id).
  for(MInt i = 0; i < receiveSize; i++) {
    // Determine cell id and append to collector
    const MInt cellId = m_tree.size();
    if(globalIdOffsets[domainId()] + (MLong)cellId != globalId(i)) {
      TERMM(1,
            "Global id mismatch." + to_string(m_domainOffsets[domainId()] + (MLong)cellId) + " "
                + to_string(globalId(i)));
    }
    m_tree.append();

    // Set tree information
    a_parentId(cellId) = parentId(i);
    for(MInt j = 0; j < IPOW2(nDim); j++) {
      a_childId(cellId, j) = childIds(i, j);
    }
    for(MInt j = 0; j < m_noDirs; j++) {
      a_neighborId(cellId, j) = nghbrIds(i, j);
    }
    a_globalId(cellId) = globalId(i);
    a_level(cellId) = level(i);
    for(MInt j = 0; j < nDim; j++) {
      a_coordinate(cellId, j) = coordinates(i, j);
    }
    a_weight(cellId) = weight(i);
    m_tree.solverBits(cellId) = solver(i);

    // Reset offsprings/workload and rebuild list of min cells
    a_noOffsprings(cellId) = 0;
    a_workload(cellId) = F0;

    // Store min-level cell information
    if(a_level(cellId) == m_minLevel) {
      localMinLevelId[i] = (MInt)localMinLevelCells.size();
      localMinLevelCells.push_back(i);

      MLong treeId = -1;
      maia::grid::hilbert::coordinatesToTreeId<nDim>(treeId, &a_coordinate(cellId, 0), (MLong)m_targetGridMinLevel,
                                                     m_targetGridCenterOfGravity, m_targetGridLengthLevel0);
      minLevelCellsTreeId.push_back(treeId);
    }

    // Set partition cell and partition level ancestor properties
    const maia::grid::tree::PropertyBitsetType cellProp = properties(i);
    const MBool isPartitionCell = cellProp[maia::grid::cell::p(Cell::IsPartitionCell)];
    const MBool isPartLvlAncestor = cellProp[maia::grid::cell::p(Cell::IsPartLvlAncestor)];
    m_tree.resetProperties(cellId);
    a_hasProperty(cellId, Cell::IsPartitionCell) = isPartitionCell;
    a_hasProperty(cellId, Cell::IsPartLvlAncestor) = isPartLvlAncestor;

    // Store global ids of partition cells
    if(isPartitionCell) {
      partitionCells.push_back(globalId(i));
    }
  }
  m_noInternalCells = m_tree.size();

  m_log << "Partition level shift: level diff of first cell is " << a_level(0) - m_minLevel << std::endl;
#ifndef NDEBUG
  std::cerr << domainId() << " Partition level shift: level diff of first cell is " << a_level(0) - m_minLevel
            << std::endl;
#endif

  std::vector<MInt>().swap(m_minLevelCells);

  // Set new global domain offsets
  std::copy_n(&globalIdOffsets[0], noDomains() + 1, &m_domainOffsets[0]);

  m_localPartitionCellOffsets[0] = partitionCellOffsets[domainId()];
  m_localPartitionCellOffsets[1] = partitionCellOffsets[domainId() + 1];
  m_localPartitionCellOffsets[2] = m_noPartitionCellsGlobal;
  m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];

  const MLong noCellsCheck = globalIdOffsets[domainId() + 1] - globalIdOffsets[domainId()];
  if(noCellsCheck != m_noInternalCells) {
    TERMM(1, "Wrong number of internal cells: " + std::to_string(m_noInternalCells)
                 + " != " + std::to_string(noCellsCheck));
  }

  TERMM_IF_COND(m_noPartitionCells <= 0, "Cannot allocate array with " + to_string(m_noPartitionCells) + " elements.");

  // Reallocate partition cell arrays with new size
  mDeallocate(m_localPartitionCellGlobalIds);
  mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
         AT_);
  mDeallocate(m_localPartitionCellLocalIds);
  mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, AT_);

  const MInt maxPartLvl = m_minLevel + m_maxPartitionLevelShift;
  // Store partition cell global ids and determine number of offsprings for all cells starting at
  // the partition level (required for communicatePartitionLevelAncestorData())
  for(MInt i = 0; i < m_noPartitionCells; i++) {
    auto cellId = (MInt)(partitionCells[i] - m_domainOffsets[domainId()]);
    const MInt partLevel = a_level(cellId);
    TERMM_IF_NOT_COND(a_hasProperty(cellId, Cell::IsPartitionCell), "Partition cell flag not set.");
    TERMM_IF_COND(partitionCells[i] != a_globalId(cellId), "Error: partition level cell global id mismatch.");
    TERMM_IF_COND(partLevel < m_minLevel || partLevel > maxPartLvl, "Invalid level for partition cell.");

    m_localPartitionCellGlobalIds[i] = partitionCells[i];
    // TODO labels:DLB,GRID @ansgar does not work properly for PLS!
    descendNoOffsprings(cellId, m_domainOffsets[domainId()]);
  }

  // Determine global min and max levels
  setLevel();

  // Rebuild global-to-local id mapping and convert global ids back to local ids
  {
    createGlobalToLocalIdMapping();

    const MInt noMinLevelCells = localMinLevelCells.size();
    if(noMinLevelCells == 0) {
      localMinLevelCells.push_back(-1);
      minLevelCellsTreeId.push_back(-1);
    }

    // Communicate the partition level ancestor data and create the missing parts of the tree from
    // the min-level up to the partition level
    maia::grid::IO<CartesianGrid<nDim>>::communicatePartitionLevelAncestorData(
        *this, noMinLevelCells, &localMinLevelCells[0], &localMinLevelId[0], &minLevelCellsTreeId[0], &cellInfo[0]);

    for(MInt i = m_noInternalCells; i < m_tree.size(); ++i) {
      ASSERT(a_hasProperty(i, Cell::IsPartLvlAncestor), "Cell is not a partition level ancestor.");
      a_hasProperty(i, Cell::IsHalo) = true;
    }

    // Set neighbor ids
    maia::grid::IO<CartesianGrid<nDim>>::propagateNeighborsFromMinLevelToMaxLevel(*this);

    // Update the global to local id mapping and change global ids in the tree into local ids
    createGlobalToLocalIdMapping();
    changeGlobalToLocalIds();

    // Rebuild window/halo cell information
    setupWindowHaloCellConnectivity();
  }

  setLevel();

  // TODO labels:GRID,DLB @ansgar_pls is this still required here? after windowHalo?
  computeGlobalIds();

  // Update local bounding box
  computeLocalBoundingBox(&m_localBoundingBox[0]);

  // Set balance status
  m_wasBalancedAtLeastOnce = true;
  m_wasBalanced = true;

#ifdef MAIA_GRID_SANITY_CHECKS
  gridSanityChecks();
#endif
}

// --------------------------------------------------------------------------------------

/**
 * \brief Delete the Connection of periodic-Halo-Cells at the bounding box.
 * \author Tim Wegmann
 */
template <MInt nDim>
void CartesianGrid<nDim>::deletePeriodicConnection(const MBool saveBackup) {
  TRACE();

  m_neighborBackup.clear();

  // For all halo cells that are periodic, reset globalId to -1
  // For all non-halo neighbors of periodic halo cells, reset neighbor id to -1
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      MInt cellId = m_haloCells[i][j];
      if(!a_hasProperty(cellId, Cell::IsPeriodic)) continue;
      if(!saveBackup) a_globalId(cellId) = -1;
      for(MInt dir = 0; dir < m_noDirs; dir++) {
        if(a_hasNeighbor(cellId, dir) == 0) continue;
        MInt nghbrId = a_neighborId(cellId, dir);
        if(!a_hasProperty(nghbrId, Cell::IsHalo)) {
          if(saveBackup) m_neighborBackup.push_back(make_tuple(cellId, nghbrId, dir));
          a_neighborId(cellId, dir) = -1;
          a_neighborId(nghbrId, m_revDir[dir]) = -1;
        }
      }
    }
  }

  if(m_azimuthalPer) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        MInt cellId = m_azimuthalHaloCells[i][j];
        ASSERT(a_hasProperty(cellId, Cell::IsPeriodic), "");
        if(!saveBackup) a_globalId(cellId) = -1;
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          if(a_hasNeighbor(cellId, dir) == 0) continue;
          MInt nghbrId = a_neighborId(cellId, dir);
          if(!a_hasProperty(nghbrId, Cell::IsHalo)) {
            if(saveBackup) m_neighborBackup.push_back(make_tuple(cellId, nghbrId, dir));
            a_neighborId(cellId, dir) = -1;
            a_neighborId(nghbrId, m_revDir[dir]) = -1;
          }
        }
      }
    }
  }
}

// --------------------------------------------------------------------------------------


/**
 * \brief Delete the Connection of periodic-Halo-Cells at the bounding box.
 * \author Tim Wegmann
 */
template <MInt nDim>
void CartesianGrid<nDim>::restorePeriodicConnection() {
  TRACE();

  for(MUint i = 0; i < m_neighborBackup.size(); i++) {
    MInt cellId = get<0>(m_neighborBackup[i]);
    MInt nghbrId = get<1>(m_neighborBackup[i]);
    MInt dir = get<2>(m_neighborBackup[i]);
    a_neighborId(cellId, dir) = nghbrId;
    a_neighborId(nghbrId, m_revDir[dir]) = cellId;
  }
  m_neighborBackup.clear();
}


// --------------------------------------------------------------------------------------

/** \brief Caluclate the number of offsprings and the workload for each cell
 *
 * 1. Assign each cell with the max possible level (as those can't have childs)
 *    the value 1 as offSprings and the value m_weight as workload.
 * 2. Assign each cell with one level up (the parents of the cells of the previous step)
 *    the value 1 + value of each child if it has any (offSprings) and the value
 *    m_weight + workload of each child if it has any (workload).
 * 3. Repeat step 2 for each level until ones reachs the minimum level.
 *
 * \tparam[in] T Cell type
 *
 **/
template <MInt nDim>
template <typename CELLTYPE>
void CartesianGrid<nDim>::calculateNoOffspringsAndWorkload(Collector<CELLTYPE>* input_cells, MInt input_noCells) {
  TRACE();
  vector<vector<MInt>> partitionLevelAncestorHaloCells;
  vector<vector<MInt>> partitionLevelAncestorWindowCells;

  if(m_maxPartitionLevelShift > 0) {
    determinePartLvlAncestorHaloWindowCells(partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells);
  } else {
    partitionLevelAncestorHaloCells.resize(noNeighborDomains());
    partitionLevelAncestorWindowCells.resize(noNeighborDomains());
  }

  maia::grid::IO<CartesianGrid<nDim>>::calculateNoOffspringsAndWorkload(
      *this, input_cells, input_noCells, m_haloCells, m_windowCells, partitionLevelAncestorWindowCells,
      partitionLevelAncestorHaloCells);
}


template <MInt nDim>
void CartesianGrid<nDim>::partitionParallel(const MInt tmpCount, const MLong tmpOffset,
                                            const MFloat* const partitionCellsWorkload,
                                            const MLong* const partitionCellsGlobalId, const MFloat totalWorkload,
                                            MLong* partitionCellOffsets, MLong* globalIdOffsets,
                                            MBool computeOnlyPartition) {
  maia::grid::IO<CartesianGrid<nDim>>::partitionParallel(*this, tmpCount, tmpOffset, partitionCellsWorkload,
                                                         partitionCellsGlobalId, totalWorkload, partitionCellOffsets,
                                                         globalIdOffsets, computeOnlyPartition);
}


/// \brief Find the partition cell containing a given coordinate (if any). If a valid solverId is
///        given it is also checked that the found partition cell belongs to this solver.
template <MInt nDim>
template <MBool t_correct>
MInt CartesianGrid<nDim>::findContainingPartitionCell(const MFloat* const coord, const MInt solverId,
                                                      function<MFloat*(MInt, MFloat* const)> correctCellCoord) {
  ASSERT(solverId >= -1 && solverId < treeb().noSolvers(), "Invalid solver id " + to_string(solverId));

  // array<MFloat, nDim> tmp;
  const MInt noLocalPartitionCells = m_localPartitionCellOffsets[1] - m_localPartitionCellOffsets[0];

  // loop over partition cells
  for(MInt i = 0; i < noLocalPartitionCells; i++) {
    const MInt curId = m_localPartitionCellLocalIds[i];

    MBool isInCell = pointWthCell<t_correct, false>(coord, curId, correctCellCoord);

    if(isInCell) {
      if(solverId < 0) {
        return curId;
      } else {
        // If a solverId is given check if partition cell belongs to solver
        if(a_solver(curId, solverId)) {
          return curId;
        } else {
          return -1;
        }
      }
    }
  }
  return -1;
}

/* J. Arvo
from "Graphics Gems", Academic Press, 1990
*/
template <MInt nDim>
MBool CartesianGrid<nDim>::boxSphereIntersection(const MFloat* bMin, const MFloat* bMax,
                                                 const MFloat* const sphereCenter, MFloat const radius) {
  MFloat dmin = 0;
  for(MInt n = 0; n < nDim; n++) {
    if(sphereCenter[n] < bMin[n]) {
      dmin += pow((sphereCenter[n] - bMin[n]), 2);
    } else if(sphereCenter[n] > bMax[n]) {
      dmin += pow((sphereCenter[n] - bMax[n]), 2);
    }
  }
  return (dmin <= pow(radius, 2));
}

/* J. Arvo
from "Graphics Gems", Academic Press, 1990
*/
template <MInt nDim>
MBool CartesianGrid<nDim>::hollowBoxSphereIntersection(const MFloat* bMin, const MFloat* bMax,
                                                       const MFloat* const sphereCenter, MFloat const radius) {
  MFloat dmin = 0;
  MBool face = false;
  for(MInt i = 0; i < nDim; i++) {
    if(sphereCenter[i] < bMin[i]) {
      face = true;
      dmin += pow(sphereCenter[i] - bMin[i], 2);
    } else if(sphereCenter[i] > bMax[i]) {
      face = true;
      dmin += pow(sphereCenter[i] - bMax[i], 2);
    } else if(sphereCenter[i] - bMin[i] <= radius) {
      face = true;
    } else if(bMax[i] - sphereCenter[i] <= radius) {
      face = true;
    }
  }
  return (face && (dmin <= pow(radius, 2)));
}

template <MInt nDim>
MInt CartesianGrid<nDim>::intersectingWithPartitioncells(MFloat* const center, MFloat const radius) {
  const MInt noLocalPartitionCells = m_localPartitionCellOffsets[1] - m_localPartitionCellOffsets[0];
  MFloat bMin[nDim];
  MFloat bMax[nDim];

  MBool intersectsACell = false;
  MInt intersectingCellId = -1;
  for(MInt i = 0; i < noLocalPartitionCells; i++) {
    const MInt curId = m_localPartitionCellLocalIds[i];
    const MFloat _halfCellLength = halfCellLength(curId);
    const MFloat* cCoord = m_tree.coordinate(curId);

    for(MInt n = 0; n < nDim; n++) {
      bMin[n] = cCoord[n] - _halfCellLength;
      bMax[n] = cCoord[n] + _halfCellLength;
    }

    intersectsACell = boxSphereIntersection(bMin, bMax, center, radius);

    if(intersectsACell) {
      intersectingCellId = curId;
    }
  }
  return intersectingCellId;
}

template <MInt nDim>
MInt CartesianGrid<nDim>::intersectingWithHaloCells(MFloat* const center, MFloat const radius, MInt domainId,
                                                    MBool onlyPartition) {
  MBool intersectsACell = false;
  MInt intersectingCellId = -1;
  std::vector<MFloat> bMin(nDim);
  std::vector<MFloat> bMax(nDim);

  for(MInt i = 0; i < noHaloCells(domainId); i++) {
    const MInt curId = haloCell(domainId, i);
    // skip all non-partition cells
    if(onlyPartition && !a_hasProperty(curId, Cell::IsPartitionCell)) {
      continue;
    }

    const MFloat _halfCellLength = halfCellLength(curId);
    const MFloat* cCoord = m_tree.coordinate(curId);

    for(MInt n = 0; n < nDim; n++) {
      // std::cout << " GRID POS: " << cCoord[n];
      bMin[n] = cCoord[n] - _halfCellLength;
      bMax[n] = cCoord[n] + _halfCellLength;
    }

    intersectsACell = boxSphereIntersection(bMin.data(), bMax.data(), center, radius);
    if(intersectsACell) {
      intersectingCellId = curId;
      break;
    }
  }
  return intersectingCellId;
}

template <MInt nDim>
MBool CartesianGrid<nDim>::intersectingWithLocalBoundingBox(MFloat* const center, MFloat const radius) {
  const MFloat* bMin = &m_localBoundingBox[0];
  const MFloat* bMax = &m_localBoundingBox[nDim];

  MBool intersecting = boxSphereIntersection(bMin, bMax, center, radius);
  return intersecting;
}

template <MInt nDim>
MInt CartesianGrid<nDim>::ancestorPartitionCellId(const MInt cellId) const {
  TERMM(1, "function not used anywhere, check this code!");

  if(a_hasProperty(cellId, Cell::IsPartitionCell)) {
    return cellId;
  }

  MInt parentId = a_parentId(cellId);
  while(!a_hasProperty(parentId, Cell::IsPartitionCell)) {
    parentId = a_parentId(parentId);
  }
  return parentId;
}

template <MInt nDim>
MInt CartesianGrid<nDim>::findContainingHaloPartitionCell(const MFloat* const coord, const MInt solverId) {
  // loop over all halo cells
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt i = 0; i < noHaloCells(n); i++) {
      const MInt curId = haloCell(n, i);

      if(!a_hasProperty(curId, Cell::IsPartitionCell)) {
        continue;
      }

      MBool isInCell = true;

      const MFloat* cellCoords = m_tree.coordinate(curId);
      const MFloat _halfCellLength = halfCellLength(curId);

      for(MInt dimId = 0; dimId < nDim; dimId++) {
        if(coord[dimId] < cellCoords[dimId] - _halfCellLength or coord[dimId] >= cellCoords[dimId] + _halfCellLength) {
          isInCell = false;
          break;
        }
      }

      if(isInCell) {
        if(solverId < 0) {
          return curId;
        } else {
          // If a solverId is given check if partition cell belongs to solver
          if(a_solver(curId, solverId)) {
            return curId;
          } else {
            return -1;
          }
        }
      }
    }
  }
  return -1;
}

template <MInt nDim>
template <MBool t_correct>
MInt CartesianGrid<nDim>::findContainingHaloCell(const MFloat* const coord, const MInt solverId, MInt domainId,
                                                 MBool onlyPartitionCells,
                                                 function<MFloat*(MInt, MFloat* const)> correctCellCoord) {
  ASSERT(solverId >= -1 && solverId < treeb().noSolvers(), "Invalid solver id " + to_string(solverId));

  MInt n = domainId;
  for(MInt i = 0; i < noHaloCells(n); i++) {
    const MInt curId = haloCell(n, i);

    if(onlyPartitionCells && !a_hasProperty(curId, Cell::IsPartitionCell)) {
      continue;
    }

    MBool isInCell = pointWthCell<t_correct, false>(coord, curId, correctCellCoord);

    if(isInCell) {
      if(solverId < 0) {
        return curId;
      } else {
        // If a solverId is given check if partition cell belongs to solver
        if(a_solver(curId, solverId)) {
          return curId;
        } else {
          return -1;
        }
      }
    }
  }
  return -1;
}


/** \brief returns the local id of the leaf cell (global or for solver) containing some coordinates p
 *
 * \author Thomas Schilden
 * \date 10.05.2016
 *
 * The algorithm does the following:
 *
 *  1. Check if the point is inside the local bounding box, else return early
 *  2. Loops over the partition cells to find the partition cell that contains p
 *  3. loops down the children to the leaf cell
 *  definition of inside: all coordinates fulfill: x_c,min <= x_p < x_c,max
 *
 * \param[in]  coordinates of the point
 * \param[out] local cell Id, -1 if not in domain
 *
 **/
template <MInt nDim>
template <MBool t_correct>
MInt CartesianGrid<nDim>::findContainingLeafCell(
    const MFloat* coord, function<MFloat*(MInt, MFloat* const)> correctCellCoord, const MInt solverId) {
  TRACE();

  // Return early if point is not inside the local bounding box
  if(!pointInLocalBoundingBox(coord)) {
    return -1;
  }

  ASSERT(solverId >= -1 && solverId < treeb().noSolvers(), "Invalid solver id " + to_string(solverId));
  // Check if the lookup should be global or on a solver basis using the given solverId
  const MBool global = (solverId < 0);

  MInt cellId = findContainingPartitionCell<t_correct>(coord, solverId, correctCellCoord);

  // return -1 if not in local partition cells
  if(cellId == -1) {
    return -1;
  }

  // loop down to leaf cell (global or for solver) containing coord
  while((global) ? a_hasChildren(cellId) : a_hasChildren(cellId, solverId)) {
    MInt childId = 0;
    for(; childId < IPOW2(nDim); childId++) {
      const MInt childCellId = (global) ? a_childId(cellId, childId) : a_childId(cellId, childId, solverId);
      if(childCellId < 0) continue;
      MBool isInCell = pointWthCell<t_correct, false>(coord, childCellId, correctCellCoord);
      if(isInCell) {
        cellId = childCellId;
        break;
      }
    }
    if(childId == IPOW2(nDim)) {
      mTerm(1, AT_, "No matching child during loop down!");
    }
  }

  if(global) {
    ASSERT(m_tree.isLeafCell(cellId), "No leaf cell... " + std::to_string(cellId));
  } else {
    ASSERT(!a_hasChildren(cellId, solverId), "No leaf cell... " + std::to_string(cellId));
  }
  return cellId;
}


/** \brief returns the local id of the leaf cell (global or for solver) containing some coordinates p
 *         based on a given startId
 *         NOTE: also cells on the same rank, but far away from the startId are considered!
 *
 * \author Sven Berger, Tim Wegmann
 * \date March 2019
 *
 * The algorithm does the following:
 *
 *  1. Determine direction of the coordinate relative to the center of the start cell
 *  2. Find neighbor cell including diagonals in that direction
 *
 * \param[in]  coordinates of the point
 * \param[in]  cellId from which to start searching
 * \param[in]  allowNonLeafHalo option to return the non-leaf halo cell which contains the coordiante
 *                              instead of returning -1 as not found! Thus the neighbor-Domain which
 *                              contains the leaf cell is known and can be used for communication!
 * \param[out] local cell Id, -1 if not found
 *
 **/
template <MInt nDim>
template <MBool t_correct>
MInt CartesianGrid<nDim>::findContainingLeafCell(const MFloat* const coord, const MInt startId,
                                                 function<MFloat*(MInt, MFloat* const)> correctCellCoord,
                                                 const MInt solverId, const MBool allowNonLeafHalo) {
  TRACE();

  // Check if the lookup should be global or on a solver basis using the given solverId
  const MBool global = (solverId < 0);

  MInt cellId = startId;
  array<MFloat, nDim> tmp;

  auto findExistingNghbr = [&](const MInt _cellId, const MInt _dir, const MInt _solverId) {
    MInt curId = _cellId;
    if(_solverId < 0) {
      while(curId > -1 && a_neighborId(_cellId, _dir) < 0) {
        curId = a_parentId(_cellId);
      }
    } else {
      while(curId > -1 && a_neighborId(curId, _dir, _solverId) < 0) {
        curId = a_parentId(curId, _solverId);
      }
    }
    if(curId > -1 && solverId > -1) {
      return a_neighborId(curId, _dir, _solverId);
    } else if(curId > -1 && solverId < 0) {
      return a_neighborId(_cellId, _dir);
    } else {
      return static_cast<MLong>(-1);
    }
  };

  set<MInt> checkedCells;
  // iterate overall neighbors till found
  while(!pointWthCell<t_correct, true>(coord, cellId, correctCellCoord)) {
    // current cellCoords
    const MFloat* tmpCoord = m_tree.coordinate(cellId);
    const MFloat* cellCoords = t_correct ? (correctCellCoord)(cellId, &tmp[0]) : tmpCoord;
    const MFloat halfCellLength = cellLengthAtLevel(a_level(cellId) + 1);

    checkedCells.insert(cellId);

    MInt direction = -1;

    // determine which cell side is crossed
    for(MInt i = 0; i < nDim; i++) {
      const MFloat difference = coord[i] - cellCoords[i];
      const MFloat absDifference = fabs(difference);
      if(absDifference > halfCellLength) {
        const MInt temp = 2 * i + ((difference > 0.0) ? 1 : 0);
        const MLong neighborId = findExistingNghbr(cellId, temp, solverId);
        if(neighborId > -1) {
          direction = temp;
          cellId = neighborId;
          break;
        } else if(a_parentId(cellId) > -1
                  && pointWthCell<t_correct, true>(coord, a_parentId(cellId), correctCellCoord)) {
          // point is within the parent cell but not within any of the cells-neighbors
          // this means, that the matching child where the point should be is outside the STL geometry
          return -1;
        }
      }
    }

    if(direction == -1) {
      // no valid cell in any direction of coords
      // meaning that the cell is to far (in all directions) from the startId
      // start a general search in the entire domain instead
      return findContainingLeafCell<t_correct>(coord, correctCellCoord, solverId);
    }

    // If tmpCoord is located on the cell surface of two cells, this loop would jump between these cells until the end
    // of time.
    if(checkedCells.find(cellId) != checkedCells.end()) break;
  }

  // loop down to leaf cell (global or for solver) containing coord
  while((global) ? a_hasChildren(cellId) : a_hasChildren(cellId, solverId)) {
    MInt childId = 0;
    for(; childId < IPOW2(nDim); childId++) {
      const MInt childCellId = (global) ? a_childId(cellId, childId) : a_childId(cellId, childId, solverId);
      if(childCellId < 0) {
        continue;
      }
      MInt cnt = 0;
      const MFloat* cellCoords = m_tree.coordinate(childCellId);
      // TODO FIXME labels:GRID
      //      const MFloat* cellCoords = t_correct ? (*correctCellCoord)(curId, &tmp[0]): tmpCoord;
      const MFloat constHalfCellLength = halfCellLength(childCellId);

      for(MInt dimId = 0; dimId < nDim; dimId++) {
        MFloat dist = coord[dimId] - cellCoords[dimId];
        if(coord[dimId] >= cellCoords[dimId] - constHalfCellLength
           && coord[dimId] < cellCoords[dimId] + constHalfCellLength) {
          cnt++;
        } else {
          const MFloat increaseOrderOfMagnitude = 10.0; // MFloatEps is to small.
          if(approx(fabs(dist), constHalfCellLength, increaseOrderOfMagnitude * MFloatEps)) {
            MInt dir = (dist > F0) ? 1 : 0;
            if(a_neighborId(childCellId, dimId + dir) > -1) {
              if(a_parentId(a_neighborId(childCellId, dimId + dir)) != a_parentId(childCellId)) {
                // Point on surface of parent cell
                cnt++;
              } else {
                if(dist > F0) {
                  // Point in center of parent cell. Choose cell with smaller cellCoords[dimId].
                  cnt++;
                }
              }
            } else {
              // Point musst be on surface of parent cell. Else neighbor would exist.
              cnt++;
            }
          }
        }
      }

      if(cnt == nDim) {
        cellId = childCellId;
        break;
      }
    }
    if(childId == IPOW2(nDim)) {
      if(!a_isHalo(cellId)) {
        mTerm(1, AT_, "No matching child during loop down!");
      } else {
        if(!allowNonLeafHalo) {
          cout << "dom: " << domainId() << " no matching child during loop down" << endl;
          return -1;
        } else {
          return cellId;
        }
      }
    }
  }

  if(global) {
    ASSERT(m_tree.isLeafCell(cellId), "No leaf cell... " + std::to_string(cellId));
  } else {
    ASSERT(!a_hasChildren(cellId, solverId), "No leaf cell... " + std::to_string(cellId));
  }

  return cellId;
}


/** \brief propagates a given distance from a given list of cells on equal level neighbours
 *
 * \author Thomas Schilden
 * \date 18.10.2016
 *
 * The algorithm does the following:
 *
 *  1. propagates locally from the local cells in the list in parameter
 *  2. exchanges the distance and repeats the propagation on the halo cells
 *
 * \param[in]  array of cells to propagate from
 * \param[in]  number of cells to propagate from
 * \param[in]  array to store the distance
 * \param[in]  final distance
 *
 **/
template <MInt nDim>
void CartesianGrid<nDim>::propagateDistance(std::vector<MInt>& list, MIntScratchSpace& distMem, MInt dist) {
  TRACE();

  m_log << "running distance propagation for distance " << dist << endl;

  for(const auto gridCellId : list)
    propagationStep(gridCellId, 0, distMem.data(), dist);

  if(noNeighborDomains() > 0) {
    MIntScratchSpace recvMemOffset(noNeighborDomains() + 1, AT_, "recvMemOffset");
    recvMemOffset.fill(0);
    MInt allRecv = 0;
    for(MInt dom = 0; dom < noNeighborDomains(); dom++) {
      allRecv += (signed)m_haloCells[dom].size();
      recvMemOffset[dom + 1] = allRecv;
    }
    MIntScratchSpace recvMem(allRecv, AT_, "recvMem");
    MInt checkHalos = 1;
    MInt nmbrComm = 0;
    m_log << "checking halos" << endl;
    while(checkHalos) {
      checkHalos = 0;
      // exchange distance
      exchangeNotInPlace(distMem.data(), recvMem.data());
      // update halo cells and propagate
      for(MInt dom = 0; dom < noNeighborDomains(); dom++) {
        for(MInt cell = 0; cell < (signed)m_haloCells[dom].size(); cell++) {
          MInt windowDist = recvMem[recvMemOffset[dom] + cell];
          MInt haloDist = distMem[m_haloCells[dom][cell]];
          if(windowDist < haloDist) {
            propagationStep(m_haloCells[dom][cell], windowDist, distMem.data(), dist);
            checkHalos = 1;
          }
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &checkHalos, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "checkHalos");
      if(checkHalos) {
        nmbrComm++;
        m_log << "checking halos again " << nmbrComm << endl;
      }
    }
  }
}

/** \brief starts the propagation on a cell and continues calling the function for
 *         the neighbours
 *
 * \author Thomas Schilden
 * \date 18.10.2016
 *
 * \param[in]  cellId
 * \param[in]  current distance
 * \param[in]  array to store the distance
 * \param[in]  final distance
 *
 **/
template <MInt nDim>
void CartesianGrid<nDim>::propagationStep(MInt cellId, MInt dist, MInt* distMem, MInt endDist) {
  TRACE();
  if(dist < distMem[cellId]) {
    distMem[cellId] = dist;
    if(dist < endDist) {
      for(MInt i = 0; i < 2 * nDim; i++) {
        const MInt nghbrId = a_neighborId(cellId, i);
        if(nghbrId > -1) propagationStep(nghbrId, dist + 1, distMem, endDist);
      }
    }
  }
}

/** \brief communicates a variable from windows to halos but does not overwrite halo
 *         values. Window values are stored in a buffer.
 *
 * \author Thomas Schilden
 * \date 18.10.2016
 *
 * \param[in]  variable to exchange
 * \param[in]  buffer memory
 *
 **/
template <MInt nDim>
template <typename DATATYPE>
void CartesianGrid<nDim>::exchangeNotInPlace(DATATYPE* exchangeVar, DATATYPE* recvMem) {
  TRACE();
  maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), exchangeVar, recvMem);
}

template <MInt nDim>
template <typename DATATYPE>
void CartesianGrid<nDim>::generalExchange(MInt noVars, const MInt* vars, DATATYPE** exchangeVar, MInt noDomSend,
                                          MInt* domSend, const MInt* noCellsSendPerDom, const MInt* cellIdsSend,
                                          MInt noDomRecv, MInt* domRecv, const MInt* noCellsRecvPerDom,
                                          const MInt* cellIdsRecv) {
  TRACE();

  // recv info
  MIntScratchSpace recvMemOffset(noDomRecv + 1, AT_, "recvMemOffset");
  recvMemOffset.fill(0);
  MInt allRecv = 0;
  for(MInt dom = 0; dom < noDomRecv; dom++) {
    allRecv += noCellsRecvPerDom[dom];
    recvMemOffset[dom + 1] = allRecv;
  }
  ScratchSpace<DATATYPE> recvMem(allRecv, AT_, "recvMem");

  // send info
  MIntScratchSpace sendMemOffset(noDomSend + 1, AT_, "sendMemOffset");
  sendMemOffset.fill(0);
  MInt allSend = 0;
  for(MInt dom = 0; dom < noDomSend; dom++) {
    allSend += noCellsSendPerDom[dom];
    sendMemOffset[dom + 1] = allSend;
  }
  ScratchSpace<DATATYPE> sendMem(allSend, AT_, "sendMem");

  DATATYPE* ptr = &(exchangeVar[0][0]);

  for(MInt i = 0; i < noVars; i++) {
    MInt var = vars[i];

    // fill the send mem
    for(MInt j = 0; j < allSend; j++)
      sendMem[j] = ptr[cellIdsSend[j] * noVars + var];

    // communicate
    ScratchSpace<MPI_Request> mpi_request_(noDomSend, AT_, "mpi_request_");

    for(MInt dom = 0; dom < noDomSend; dom++)
      MPI_Issend(&(sendMem[sendMemOffset[dom]]), noCellsSendPerDom[dom] * sizeof(DATATYPE), MPI_CHAR, domSend[dom], 0,
                 mpiComm(), &mpi_request_[dom], AT_, "(sendMem[sendMemOffset[dom]])");

    MPI_Status status_;
    for(MInt dom = 0; dom < noDomRecv; dom++)
      MPI_Recv(&(recvMem[recvMemOffset[dom]]), noCellsRecvPerDom[dom] * sizeof(DATATYPE), MPI_CHAR, domRecv[dom], 0,
               mpiComm(), &status_, AT_, "(recvMem[recvMemOffset[dom]])");

    for(MInt dom = 0; dom < noDomSend; dom++)
      MPI_Wait(&mpi_request_[dom], &status_, AT_);

    // distribute the recv mem
    for(MInt j = 0; j < allRecv; j++)
      ptr[cellIdsRecv[j] * noVars + var] = recvMem[j];
  }
}


/** \brief Overload method of createGridSlice to provide usage with or without
 *         return of cellIds which are in the slice
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 01.06.2016
 **/
template <MInt nDim>
void CartesianGrid<nDim>::createGridSlice(const MString& axis, const MFloat intercept, const MString& fileName) {
  TRACE();

  createGridSlice(axis, intercept, fileName, -1, nullptr, nullptr, nullptr, nullptr);
}


/** \brief Create a valid 2D grid which represents a slice from a 3D grid.
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 06.01.2016
 *
 *  \param[in] axis Axis to which the slice should be orthogonal ("x", "y", or "z").
 *  \param[in] intercept Absolute value of the axis position for the slice.
 *  \param[in] fileName File name (including output directory path) to which the slice grid should be written.
 *  \param[out] noSliceCellIds Stores the number of cells in the slice from this domain. If nullptr
 *                             is given, no data is stored.
 *  \param[out] sliceCellIds Stores a list of all cells that are in the slice. Thus, a pointer to
 *                           storage of size >= m_noInternalCells should be provided. If nullptr is
 *                           given, no data is stored.
 *  \param[out] noSliceHilbertIds Stores the number of hilbertIds in the slice from this domain. If
 *                                nullptr is given, no data is stored.
 *  \param[out] sliceHilbertIds Stores a list of all hilbertIds of the slice. Thus, a pointer of the
 *                              storage of size >= m_noInternalCells * 3 should be provided. If
 *                              nullptr is given, no data is stored.
 *
 *  This method is structured in 4 parts:
 *  Step 1: Identify cells in slice and sort them by their sliceHilbertId. First the data is saved
 *          a collector and after identification sorted. Then the data is distrubted into single
 *          arrays.
 *  Step 2: Create MPI communicator only for domains with slice cells
 *  Step 3: Gather all data for slice grid file. First the hilbertId information is exchanged
 *          globally. Then the global offsets for the local data can be determined on every domain.
 *          It is important to understand, that the number of hilbertIds is different on every
 *          domain and that they can be distributed in in any order, since a slice of a 3D hilbert
 *          curve is performed. After that the neighbor domains are determined to exchange cellIds
 *          to determine the cell mapping to slice globalId completly. At the end the data arrays
 *          are determined.
 *  Step 4: Write grid file with slice data. Since the data is sorted by slice HilbertId, the data
 *          writing is done in data chunks with the same hilbertId. If a domain has less hilbertIds
 *          than others, then dummy calls a performed.
 **/
template <MInt nDim>
void CartesianGrid<nDim>::createGridSlice(const MString& axis,
                                          const MFloat intercept,
                                          const MString& fileName,
                                          const MInt solverId,
                                          MInt* const noSliceCellIds,
                                          MInt* const sliceCellIds,
                                          MInt* const noSliceHilbertIds,
                                          MInt* const sliceHilbertInfo,
                                          MInt* const noSliceContHilbertIds,
                                          MInt* const sliceContiguousHilbertInfo) {
  TRACE();

  // Check if maximum and minimum cell levels are equal
  if(m_minLevel == m_maxLevel) {
    if(domainId() == 0) {
      // TODO FIXME labels:GRID what is the reason for this?
      std::cerr << "Warning: CartesianGrid::createGridSlice minLevel and maxLevel grid are equal! "
                << "(hang-up might occur)" << std::endl;
    }
  }

  // Dimension check
  IF_CONSTEXPR(nDim != 3) { TERMM(1, "Can not create a 2D slice from a 2D grid."); }

  m_log << "Creating 2D slice from grid... ";

  // New number of dimension
  const MInt nDimSlice = 2;

  // Check for slice orientation
  MInt axisNum;
  // Reference to relevant coordinates
  array<MInt, nDimSlice> coordArray{};
  if(axis == "x") {
    axisNum = 0;
    coordArray = {{1, 2}};
  } else if(axis == "y") {
    axisNum = 1;
    coordArray = {{0, 2}};
  } else if(axis == "z") {
    axisNum = 2;
    coordArray = {{0, 1}};
  } else {
    TERMM(1, "Unknown axis! Check your property file.");
  }

  MBool isAtCenterOfGravity = false;
  // Check whether intercept is located at the centroid of the specified axis
  if(approx(intercept, m_targetGridCenterOfGravity[axisNum], MFloatEps)) {
    isAtCenterOfGravity = true;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 1: identify cells in slice
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Get center of gravity (needed for hilbertId in slice)
  std::array<MFloat, nDimSlice> centerOfGravity{};
  centerOfGravity[0] = m_centerOfGravity[coordArray[0]];
  centerOfGravity[1] = m_centerOfGravity[coordArray[1]];

  std::array<MFloat, nDimSlice> targetGridCenterOfGravity{};
  targetGridCenterOfGravity[0] = m_targetGridCenterOfGravity[coordArray[0]];
  targetGridCenterOfGravity[1] = m_targetGridCenterOfGravity[coordArray[1]];

  // Contains all cell ids relevant to slice
  MInt noSlicePartitionCells = 0;
  MInt noSliceMinLevelCells = 0;
  MInt noSliceCells = 0;
  MLong noLeafCells = 0;
  // SliceCellCollector for all local sliceCells with localCellId, hilbertId, isPartitionCell and
  // isMinCell information. This collector is useful for sorting this data in groups
  // After sorting the data is separated into single arrays
  ScratchSpace<array<MInt, 4>> sliceCellCollector(m_noInternalCells, AT_, "sliceCellCollector");
  for(MInt i = 0; i < m_noInternalCells; i++) {
    sliceCellCollector[i].fill(0);
  }
  // Determine partition cells of slice
  const MInt noLocalPartitionCells = m_localPartitionCellOffsets[1] - m_localPartitionCellOffsets[0];

  // Scratch to access partitionCells in following for loop
  MIntScratchSpace localPartitionCellGlobalIds(noLocalPartitionCells + 1, AT_, "localPartitionCellGlobalIds");
  copy_n(&m_localPartitionCellGlobalIds[0], noLocalPartitionCells, &localPartitionCellGlobalIds[0]);
  // This partitionCell would be located on the next domain, slice cell search on this
  // domains is perfmored until this global cellId
  localPartitionCellGlobalIds[noLocalPartitionCells] = m_noInternalCells + m_domainOffsets[domainId()];

  // Count how many cells belong to a hilbertId (for this solver) on this domain. Needed for mapping and file writing
  std::map<MInt, MInt> hilbertIdCount;
  std::map<MInt, MInt> hilbertIdCountSolver;

  for(MInt i = 0, cellId = 0; i < noLocalPartitionCells; i++) {
    // get local partition cellId
    const MInt partitionCellId = localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
    MInt isPartitionCell = 0;
    // Skip partition cell if not intersected by slice
    MFloat eps = 0.0;
    if(isAtCenterOfGravity) {
      eps = MFloatEps;
    }
    if(fabs((intercept + eps) - a_coordinate(partitionCellId, axisNum)) < halfCellLength(partitionCellId)) {
      isPartitionCell = 1;
      noSlicePartitionCells++;
    } else {
      // Skip all offsprings
      cellId = partitionCellId + a_noOffsprings(partitionCellId);
    }

    // Search for all cells between two partitionCells and collect information
    const MInt nextPartitionCell = localPartitionCellGlobalIds[i + 1] - m_domainOffsets[domainId()];
    while(cellId < nextPartitionCell) {
      // Skip cell if not intersected by slice
      if(!(fabs((intercept + eps) - a_coordinate(cellId, axisNum)) < halfCellLength(cellId))) {
        cellId++;
        continue;
      }

      // save partitionCell indicator
      if(cellId == partitionCellId) {
        sliceCellCollector[noSliceCells][2] = isPartitionCell;
        isPartitionCell = 0;
      }

      // Count hilbertIds of slice in 2D grid
      const MFloat* const sliceCoord = &a_coordinate(cellId, 0);
      // Coordinates mapped to unit cube
      array<MFloat, 2> normedSliceCoord{};
      normedSliceCoord[0] =
          (sliceCoord[coordArray[0]] - targetGridCenterOfGravity[0] + 1e-12 + m_targetGridLengthLevel0 * 0.5)
          / m_targetGridLengthLevel0;
      normedSliceCoord[1] =
          (sliceCoord[coordArray[1]] - targetGridCenterOfGravity[1] + 1e-12 + m_targetGridLengthLevel0 * 0.5)
          / m_targetGridLengthLevel0;
      const MLong sliceHilbertId =
          maia::grid::hilbert::index<2>(&normedSliceCoord[0], static_cast<MLong>(m_targetGridMinLevel));
      // count cells per slice hilbertId
      hilbertIdCount[sliceHilbertId]++;
      sliceCellCollector[noSliceCells][1] = sliceHilbertId;

      // Count solver-cells per slice hilbertId
      if(solverId > -1 && a_solver(cellId, solverId)) hilbertIdCountSolver[sliceHilbertId]++;

      // Check if current cell is a minLevelCell
      if(a_level(cellId) == m_minLevel) {
        sliceCellCollector[noSliceCells][3] = 1;
        noSliceMinLevelCells++;
      }

      // Count number of leaf cells
      if(a_isLeafCell(cellId)) {
        noLeafCells++;
      }
      // Add cell to list of slice cells
      sliceCellCollector[noSliceCells][0] = cellId;
      noSliceCells++;

      // next cell id
      cellId++;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 1.1: sort slice cells by slice hilbertÍd
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Allocate memory for slice  partition cells
  MIntScratchSpace slicePartitionCells(noSlicePartitionCells + 1, AT_, "slicePartitionCells");
  // Allocate memory for slice minLevel cells
  MIntScratchSpace sliceMinLevelCells(noSliceMinLevelCells, AT_, "sliceMinLevelCells");
  // Allocate memory for slice cells
  MIntScratchSpace cellIdsInSlice(noSliceCells, AT_, "cellIdsInSlice");
  // save hilbert id keys for easy access and sort slice cell arrays
  MIntScratchSpace hIdKeys(hilbertIdCount.size(), AT_, "hIdKeys");
  // Count and save hilbertIds for partition and minLevelCells
  std::map<MInt, MInt> hilbertIdMinLevelCellsCount;
  std::map<MInt, MInt> hilbertIdPartitionCellsCount;

  if(noSliceCells > 0) {
    // create hilbert key array for easy access
    MInt noHilbertIds = 0;
    for(auto const& ent1 : hilbertIdCount) {
      hIdKeys[noHilbertIds] = ent1.first;
      noHilbertIds++;
    }

    // sort slice cell info by slice hilbertId
    std::stable_sort(sliceCellCollector.begin(), sliceCellCollector.begin() + noSliceCells,
                     [](const array<MInt, 4>& a, const array<MInt, 4>& b) { return a[1] < b[1]; });

    // Split slice cell info into seperate arrays
    for(MInt i = 0, j = 0, k = 0; i < noSliceCells; i++) {
      cellIdsInSlice[i] = sliceCellCollector[i][0];
      // fill slicePartitionCells
      if(sliceCellCollector[i][2] == 1) {
        slicePartitionCells[j] = sliceCellCollector[i][0];
        hilbertIdPartitionCellsCount[sliceCellCollector[i][1]]++;
        j++;
      }
      // fill sliceMinLevelCells
      if(sliceCellCollector[i][3] == 1) {
        sliceMinLevelCells[k] = sliceCellCollector[i][0];
        hilbertIdMinLevelCellsCount[sliceCellCollector[i][1]]++;
        k++;
      }
    }
    // Save total number of slice partitionCells (need for noOffsprings search)
    slicePartitionCells[noSlicePartitionCells] = cellIdsInSlice[noSliceCells - 1] + 1;
  }

  if(solverId < 0) {
    // return noSliceCells and cellIds
    if(noSliceCellIds != nullptr) {
      *noSliceCellIds = noSliceCells;
    }
    if(sliceCellIds != nullptr && noSliceCells > 0) {
      std::copy_n(&cellIdsInSlice[0], noSliceCells, sliceCellIds);
    }
  } else {
    // Find number of cells and cell ids for given solverId
    MIntScratchSpace cellIdsInSliceSolver(noSliceCells, AT_, "cellIdsInSliceSolver");
    MInt noSliceCellsSolver = 0;
    for(MInt i = 0; i < noSliceCells; i++) {
      if(a_solver(cellIdsInSlice[i], solverId)) {
        cellIdsInSliceSolver[noSliceCellsSolver] = cellIdsInSlice[i];
        noSliceCellsSolver++;
      }
    }

    // return number of cells for solver and cell ids
    if(noSliceCellIds != nullptr) {
      *noSliceCellIds = noSliceCellsSolver;
    }
    if(sliceCellIds != nullptr && noSliceCellsSolver > 0) {
      std::copy_n(&cellIdsInSliceSolver[0], noSliceCellsSolver, sliceCellIds);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 2: create MPI communicator only for domains with slice cells
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Create a new MPI Communicator for all slice domains
  MPI_Comm mpiCommSlice = MPI_COMM_NULL;
  MIntScratchSpace sliceRanks(globalNoDomains(), AT_, "sliceRanks");

  // Get Rank if points are relevant
  MInt rank = -1;
  if(noSliceCells > 0) {
    MPI_Comm_rank(mpiComm(), &rank);
  }
  // Combine all ranks to get relevant ranks
  MPI_Allgather(&rank, 1, type_traits<MInt>::mpiType(), &sliceRanks[0], 1, type_traits<MInt>::mpiType(), mpiComm(), AT_,
                "rank", "sliceRanks[0]");

  // Check for relevant ranks and save them to create the new communicator
  const MInt noRelDomains = count_if(sliceRanks.begin(), sliceRanks.end(), [](const MInt a) { return a != -1; });
  MIntScratchSpace relDomains(noRelDomains, AT_, "relDomains");
  // A domainMap is needed to refer from sliceDomain to global domains
  map<MInt, MInt> domainMap;
  MInt position = 0;
  for(auto&& slice : sliceRanks) {
    if(slice != -1) {
      relDomains[position] = slice;
      domainMap[slice] = position;
      position++;
    }
  }
  // Create new point data mpi group
  MPI_Group globalGroup, localGroup;
  MPI_Comm_group(mpiComm(), &globalGroup, AT_, "globalGroup");
  MPI_Group_incl(globalGroup, noRelDomains, &relDomains[0], &localGroup, AT_);

  // Create new communicator and clean up
  MPI_Comm_create(mpiComm(), localGroup, &mpiCommSlice, AT_, "mpiCommSlice");

  MPI_Group_free(&globalGroup, AT_);
  MPI_Group_free(&localGroup, AT_);

  // Leave function if not relevant to slice
  if(noSliceCells < 1) {
    return;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3: gather all data for slice grid file
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // This part determines the new global cells order by slice hilbertId. Information about global
  // order must be avialable on all domains that the local order can be determined on every domain.
  // The new global order determines the offset for file writing and this is computed in this step.

  // Determine total global number of slice cells for each domain
  MInt noSliceDomains = -1;
  MPI_Comm_size(mpiCommSlice, &noSliceDomains);
  // Get domainId in slice mpi comm to access domainOffset
  MInt sliceDomain = -1;
  MPI_Comm_rank(mpiCommSlice, &sliceDomain);

  // Determine total no of slice hilbertIds
  MIntScratchSpace noLocalHilbertIds(noSliceDomains, AT_, "noLocalHilbertIds");
  noLocalHilbertIds[sliceDomain] = hilbertIdCount.size();

  // Gather the number of unique number of slice hilbert ids on each domain
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, &noLocalHilbertIds[0], 1, MPI_INT, mpiCommSlice, AT_, "MPI_IN_PLACE",
                "noLocalHilbertIds[0]");

  // Compute offsets for recieveing data
  MIntScratchSpace offsetsRecvData(noSliceDomains, AT_, "offsetsRecvData");
  offsetsRecvData[0] = 0;
  for(MInt i = 1; i < noSliceDomains; i++) {
    offsetsRecvData[i] = offsetsRecvData[i - 1] + noLocalHilbertIds[i - 1] * 6;
  }
  MInt noTotalRecvData = offsetsRecvData[noSliceDomains - 1] + noLocalHilbertIds[noSliceDomains - 1] * 6;

  // Create exchange data (every domain gets info about hilbert ids and offsets for cellInfo,
  // partitionCells and minLevelCells). From this offsets for file writing are determined for every
  // variable indepently.
  MIntScratchSpace sendHilbertData(noLocalHilbertIds[sliceDomain] * 6, AT_, "sendHilbertData");
  for(MUlong i = 0; i < hIdKeys.size(); i++) {
    sendHilbertData[i * 6] = hIdKeys[i];
    sendHilbertData[i * 6 + 1] = domainId();
    sendHilbertData[i * 6 + 2] = hilbertIdCount[hIdKeys[i]];
    sendHilbertData[i * 6 + 3] = hilbertIdPartitionCellsCount[hIdKeys[i]];
    sendHilbertData[i * 6 + 4] = hilbertIdMinLevelCellsCount[hIdKeys[i]];
    sendHilbertData[i * 6 + 5] = hilbertIdCountSolver[hIdKeys[i]];
  }

  // Exchange number of slice hilbert ids
  MIntScratchSpace recvHilbertData(noTotalRecvData, AT_, "recvHilbertData");
  MInt noSendData = sendHilbertData.size();
  MIntScratchSpace noRecvData(noSliceDomains, AT_, "noRecvData");
  for(MInt i = 0; i < noSliceDomains; i++) {
    noRecvData[i] = noLocalHilbertIds[i] * 6;
  }
  MPI_Allgatherv(&sendHilbertData[0], noSendData, MPI_INT, &recvHilbertData[0], &noRecvData[0], &offsetsRecvData[0],
                 MPI_INT, mpiCommSlice, AT_, "sendHilbertData[0]", "recvHilbertData[0]");

  // Create sliceGlobalHilbertInfo
  ScratchSpace<array<MInt, 6>> sliceGlobalHilbertInfo(noTotalRecvData / 6, AT_, "sliceGlobalHilbertInfo");
  // Store globally recieved data in sliceGlobalHilbertInfo
  for(MInt i = 0; i < noTotalRecvData / 6; i++) {
    sliceGlobalHilbertInfo[i][0] = recvHilbertData[i * 6];
    sliceGlobalHilbertInfo[i][1] = recvHilbertData[i * 6 + 1];
    sliceGlobalHilbertInfo[i][2] = recvHilbertData[i * 6 + 2];
    sliceGlobalHilbertInfo[i][3] = recvHilbertData[i * 6 + 3];
    sliceGlobalHilbertInfo[i][4] = recvHilbertData[i * 6 + 4];
    sliceGlobalHilbertInfo[i][5] = recvHilbertData[i * 6 + 5];
  }

  // Sort sliceGlobalHilbertInfo by domainId...
  std::stable_sort(sliceGlobalHilbertInfo.begin(), sliceGlobalHilbertInfo.end(),
                   [](const array<MInt, 6>& a, const array<MInt, 6>& b) { return a[1] < b[1]; });
  // ... and then sort stable by hilbertId -> list ordered by hilbertId and then nested for each hilbertId by domainId
  std::stable_sort(sliceGlobalHilbertInfo.begin(), sliceGlobalHilbertInfo.end(),
                   [](const array<MInt, 6>& a, const array<MInt, 6>& b) { return a[0] < b[0]; });

  // Determine hilbertId offsets (cells in slice file will be sorted by hilbertId)
  // For I/O it is necessary to call writeData() on all domains equally. If a slice has less then
  // maxNoHilbertIds then its array is filled with zero. Then a cell of writeData() will happen, but
  // without writing data.
  MInt maxNoHilbertIds = *std::max_element(noLocalHilbertIds.begin(), noLocalHilbertIds.end());
  MIntScratchSpace hilbertDomainOffset(maxNoHilbertIds, AT_, "hilbertDomainOffset");
  hilbertDomainOffset.fill(0);
  MIntScratchSpace hilbertPartitionDomainOffset(maxNoHilbertIds, AT_, "hilbertDomainOffset");
  hilbertPartitionDomainOffset.fill(0);
  MIntScratchSpace hilbertMinLevelDomainOffset(maxNoHilbertIds, AT_, "hilbertDomainOffset");
  hilbertMinLevelDomainOffset.fill(0);

  // Offsets for solver cells
  MIntScratchSpace hilbertDomainOffsetSolver(maxNoHilbertIds, AT_, "hilbertDomainOffsetSolver");
  hilbertDomainOffsetSolver.fill(0);

  // Calculate offsets by hilbert ids for output arrays
  // The data in sliceGlobalHilbertInfo is sorted by their hilbertId on the first level. On the
  // second level by domainId. Example:
  //                             [hilbertId, domaindId, noCells, noPartitionCells noMinLevelCells]
  // sliceGlobalHilbertInfo[0] = [0, 0, 20, 1, 1]
  // sliceGlobalHilbertInfo[1] = [0, 2, 10, 1, 0]
  // sliceGlobalHilbertInfo[2] = [1, 1, 6, 0 ,0]
  // sliceGlobalHilbertInfo[3] = [1, 2, 8, 1 ,1]
  //
  // To determine the global offset for the local data (for each local hilbertId), count the number
  // of cells of all hilbertIds below the current one. If this hilbertId is distributed on
  // multiple domains, then count also the next number of nodes till the local domainid is reached
  // in sliceGlobalHilbertInfo. In the example above the offset for hilbertId 1 on domain 1 would
  // be 30.
  for(MInt i = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
    MInt offset = 0, offsetPart = 0, offsetMin = 0, offsetSolver = 0, j = 0;
    // Count until current hilbertId is reached
    while(sliceGlobalHilbertInfo[j][0] < hIdKeys[i]) {
      offset += sliceGlobalHilbertInfo[j][2];
      offsetPart += sliceGlobalHilbertInfo[j][3];
      offsetMin += sliceGlobalHilbertInfo[j][4];
      offsetSolver += sliceGlobalHilbertInfo[j][5];
      j++;
    }
    // Count until this domain is reached
    while(sliceGlobalHilbertInfo[j][1] < domainId()) {
      offset += sliceGlobalHilbertInfo[j][2];
      offsetPart += sliceGlobalHilbertInfo[j][3];
      offsetMin += sliceGlobalHilbertInfo[j][4];
      offsetSolver += sliceGlobalHilbertInfo[j][5];
      j++;
    }
    hilbertDomainOffset[i] = offset;
    hilbertPartitionDomainOffset[i] = offsetPart;
    hilbertMinLevelDomainOffset[i] = offsetMin;

    if(solverId > -1) {
      hilbertDomainOffsetSolver[i] = offsetSolver;
      TERMM_IF_COND(!g_multiSolverGrid && offsetSolver != offset, "Error: fixme");
    }
  }

  // Contains offset for every domain to determine global slice cellId
  MIntScratchSpace domainOffset(noSliceDomains + 1, AT_, "domainOffset");
  // Collect number of slice cells for each domain
  MPI_Allgather(&noSliceCells, 1, type_traits<MInt>::mpiType(), &domainOffset[0], 1, type_traits<MInt>::mpiType(),
                mpiCommSlice, AT_, "noSliceCells", "domainOffset[0]");

  // Calculate offsets from received noSliceCells
  for(MInt i = 0, offset = 0, tmp = 0; i < (noSliceDomains + 1); i++) {
    tmp = domainOffset[i];
    domainOffset[i] = offset;
    offset += tmp;
  }

  MInt noLocalHilbertIdsSolver = 0;

  if(solverId < 0) {
    // return noHilbertIds
    if(noSliceHilbertIds != nullptr) {
      *noSliceHilbertIds = noLocalHilbertIds[sliceDomain];
    }
    // return HilbertInfo (how many cells per HilbertId and offset for file writing)
    if(sliceHilbertInfo != nullptr) {
      MIntScratchSpace hilbertInfo(noLocalHilbertIds[sliceDomain] * 3, AT_, "hilbertInfo");
      for(MInt i = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
        hilbertInfo[i * 3] = hIdKeys[i];
        hilbertInfo[i * 3 + 1] = hilbertIdCount[hIdKeys[i]];
        hilbertInfo[i * 3 + 2] = hilbertDomainOffset[i];
      }
      std::copy_n(&hilbertInfo[0], noLocalHilbertIds[sliceDomain] * 3, sliceHilbertInfo);
    }
  } else {
    TERMM_IF_COND(noSliceHilbertIds == nullptr || sliceHilbertInfo == nullptr,
                  "Error: solverId given, pointers need to be != nullptr.");

    MIntScratchSpace hilbertInfoSolver(noLocalHilbertIds[sliceDomain] * 3, AT_, "hilbertInfo");
    MInt noHIdsSolver = 0;
    for(MInt i = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
      if(hilbertIdCountSolver[hIdKeys[i]] > 0) {
        hilbertInfoSolver[noHIdsSolver * 3] = hIdKeys[i];
        hilbertInfoSolver[noHIdsSolver * 3 + 1] = hilbertIdCountSolver[hIdKeys[i]];
        hilbertInfoSolver[noHIdsSolver * 3 + 2] = hilbertDomainOffsetSolver[i];
        noHIdsSolver++;
      }
    }

    noLocalHilbertIdsSolver = noHIdsSolver;
    // return noHilbertIds for solver
    *noSliceHilbertIds = noHIdsSolver;
    // return HilbertInfo (how many solver-cells per HilbertId and offset for file writing)
    std::copy_n(&hilbertInfoSolver[0], noHIdsSolver * 3, sliceHilbertInfo);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3.1: determine neighbor domains and compute partition data
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // The neighbor domains are needed to create correct mapping of neighbor cells which are on other
  // domains

  // Reference to relevant neighbor cells
  array<MInt, nDimSlice * 2> nghbrArray{};
  if(axis == "x") {
    nghbrArray = {{2, 3, 4, 5}};
  } else if(axis == "y") {
    nghbrArray = {{0, 1, 4, 5}};
  } else if(axis == "z") {
    nghbrArray = {{0, 1, 2, 3}};
  }

  // Collector for related domains of slice to create global cellId map
  std::set<MInt> relatedDomains;

  // Check for neighbor cell of all slice cells on other domain
  for(MInt i = 0; i < noSliceMinLevelCells; i++) { // TODO labels:GRID this is not exactly what the comment above states
    // Only check for neighbors in slice directions
    for(MInt j = 0; j < (nDimSlice * 2); j++) {
      // Skip if no neighbor in current direction
      if(!a_hasNeighbor(sliceMinLevelCells[i], nghbrArray[j])) {
        continue;
      }
      // Save related domain
      MInt nId = a_neighborId(sliceMinLevelCells[i], nghbrArray[j]);
      const MInt nghbrDomainId = findNeighborDomainId(a_globalId(nId));
      if(domainId() != nghbrDomainId) {
        relatedDomains.insert(nghbrDomainId);
      }
    }
  }

  // Determine partition level shift in slice
  MInt partitionLevelShift = 0;
  for(MInt i = 0; i < noSlicePartitionCells; ++i) {
    const MInt levelDiff = a_level(slicePartitionCells[i]) - m_minLevel;
    partitionLevelShift = mMax(levelDiff, partitionLevelShift);
  }

  // Determine partition level ancestors (parent cells of partition cells)
  set<MLong> partitionLevelAncestorIds;
  for(MInt i = 0; i < noSliceCells; i++) {
    const MInt cellId = cellIdsInSlice[i];
    if(a_hasProperty(cellId, Cell::IsPartLvlAncestor) && !a_hasProperty(cellId, Cell::IsHalo)) {
      partitionLevelAncestorIds.insert(cellId);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3.2: Map (and exchange) globalIds to new globalIds in slice
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Create idMap to map from globalId to global slice cellId for the local cells
  map<MInt, MInt> idMap;
  // CellId -1 should be invariant
  idMap[-1] = -1;
  // Map global cellIds of locally founded slice cells to global slice cellIds
  for(MInt i = 0, j = 0, h = 0; i < noSliceCells; i++, j++) {
    // set next hilbertr offset and reset counter
    if(j == hilbertIdCount[hIdKeys[h]]) {
      h++;
      j = 0;
    }
    idMap[cellIdsInSlice[i] + m_domainOffsets[domainId()]] = j + hilbertDomainOffset[h];
  }

  // If slice has no related domains, skip slice data exchange step
  if(!relatedDomains.empty()) {
    // Get unique domainIds form set
    std::vector<MInt> sliceRelatedDomains(relatedDomains.begin(), relatedDomains.end());

    // All cellIds from other domains will be saved in a single vector.
    // Determine offset for this vector for each
    MInt noRecvCells = 0;
    MIntScratchSpace relDomainOffsets(sliceRelatedDomains.size() + 1, AT_, "relDomainOffsets");
    for(MUlong i = 0; i < sliceRelatedDomains.size(); i++) {
      relDomainOffsets[i] = noRecvCells;
      noRecvCells +=
          domainOffset[domainMap[sliceRelatedDomains[i]] + 1] - domainOffset[domainMap[sliceRelatedDomains[i]]];
    }
    // Save total number of cellIds which will be received
    relDomainOffsets[sliceRelatedDomains.size()] = noRecvCells;

    // Scratch to save all cellIds from other domains
    MIntScratchSpace recvIds(noRecvCells, AT_, "recvIds");
    recvIds.fill(-1);

    // Save mpi requests
    ScratchSpace<MPI_Request> sendRequests(sliceRelatedDomains.size(), AT_, "sendRequests");
    fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
    ScratchSpace<MPI_Request> recvRequests(sliceRelatedDomains.size(), AT_, "recvRequests");
    fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);

    // Start receiving
    for(MUlong i = 0; i < sliceRelatedDomains.size(); i++) {
      MInt d = sliceRelatedDomains[i];
      MPI_Irecv(&recvIds[relDomainOffsets[i]], domainOffset[domainMap[d] + 1] - domainOffset[domainMap[d]],
                type_traits<MInt>::mpiType(), domainMap[d], domainMap[d], mpiCommSlice, &recvRequests[i], AT_,
                "recvIds[relDomainOffsets[i]]");
    }

    // Start sending
    for(MUlong i = 0; i < sliceRelatedDomains.size(); i++) {
      MPI_Isend(&cellIdsInSlice[0], noSliceCells, type_traits<MInt>::mpiType(), domainMap[sliceRelatedDomains[i]],
                domainMap[domainId()], mpiCommSlice, &sendRequests[i], AT_, "cellIdsInSlice[0]");
    }

    // Finish receiving
    MPI_Waitall(sliceRelatedDomains.size(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);

    // Finish sending
    MPI_Waitall(sliceRelatedDomains.size(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

    // Map relevant cellIds for each related domain
    for(MUlong j = 0; j < sliceRelatedDomains.size(); j++) {
      MInt d = sliceRelatedDomains[j];
      std::vector<MInt> hId(noLocalHilbertIds[domainMap[d]]);
      std::vector<MInt> hCount(noLocalHilbertIds[domainMap[d]]);
      std::vector<MInt> hOffsets(noLocalHilbertIds[domainMap[d]]);

      // Calculate hilbertDomainOffset for related domain
      for(MUlong k = 0, i = 0; i < hId.size(); k++) {
        if(sliceGlobalHilbertInfo[k][1] == d) {
          hId[i] = sliceGlobalHilbertInfo[k][0];
          hCount[i] = sliceGlobalHilbertInfo[k][2];

          MInt m = 0;
          MInt offset = 0;
          while(sliceGlobalHilbertInfo[m][0] < hId[i]) {
            offset += sliceGlobalHilbertInfo[m][2];
            m++;
          }
          while(sliceGlobalHilbertInfo[m][1] < d) {
            offset += sliceGlobalHilbertInfo[m][2];
            m++;
          }
          hOffsets[i] = offset;
          i++;
        }
      }

      // Map global cellIds of related cellIds to global slice cellIds
      for(MInt i = 0, n = 0, h = 0; i < (relDomainOffsets[j + 1] - relDomainOffsets[j]); i++, n++) {
        if(n == hCount[h]) {
          h++;
          n = 0;
        }
        idMap[recvIds[i + relDomainOffsets[j]] + m_domainOffsets[d]] = n + hOffsets[h];
      }
    }
  }
  // Check if only cells where found, but no partitionCells (e.g. parent is located on this domain,
  //  but the first partitionCell which belongs to the slice is on next domain
  MBool onlySliceCells = false;
  if(noSlicePartitionCells < 1) {
    onlySliceCells = true;
  }

  // Create articifal mapped id for partitioncell on next domain, and update slicePartitionCells.
  // Needed for computation of workload.
  if(noSlicePartitionCells > 0) {
    std::map<MInt, MInt>::iterator it;
    MInt lastKeyId = hIdKeys.size() - 1;
    MInt val = hilbertIdCount[hIdKeys[lastKeyId]] + hilbertDomainOffset[lastKeyId];

    it = idMap.find(slicePartitionCells[noSlicePartitionCells] + m_domainOffsets[domainId()]);

    auto result = std::find_if(idMap.begin(), idMap.end(), [val](const auto& mo) { return mo.second == val; });

    if(it == idMap.end() && result == idMap.end()) {
      idMap[slicePartitionCells[noSlicePartitionCells] + m_domainOffsets[domainId()]] =
          hilbertIdCount[hIdKeys[lastKeyId]] + hilbertDomainOffset[lastKeyId];
    } else if(result != idMap.end()) {
      MInt foundKey = result->first;
      slicePartitionCells[noSlicePartitionCells] = foundKey - m_domainOffsets[domainId()];
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3.3: determine partition cell information and workload
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // New number of childs in a refinement step
  const MUint noChildInRef = pow(2, nDimSlice);

  // Create Variables to hold slice grid information
  MIntScratchSpace partitionCellsId(noSlicePartitionCells, AT_, "partitionCellsId");
  MFloatScratchSpace partitionCellsWorkload(noSlicePartitionCells, AT_, "partitionCellsWorkLoad");
  MInt maxNoOffsprings = 0;
  MFloat maxWorkload = 0.0;
  MFloat totalWorkload = 0.0;

  for(MInt i = 0, c = 0, partitionCellNoOffsprings = 0; i < noSlicePartitionCells; i++) {
    partitionCellsId[i] = idMap[slicePartitionCells[i] + m_domainOffsets[domainId()]];
    // Count number of offsprings for partitionCells (start with 1)
    partitionCellNoOffsprings = 1;
    // Search all sliceCells between 2 slice partitionCells for offsprings. The id mapping is used
    // here since the mapped ids are sorted by the slice hilbertId
    while(idMap[a_globalId(cellIdsInSlice[c])] < idMap[slicePartitionCells[i + 1] + m_domainOffsets[domainId()]]) {
      // Only count cells which are on a higher level than current partitionCell
      if(a_level(cellIdsInSlice[c]) > a_level(slicePartitionCells[i])) {
        partitionCellNoOffsprings++;
      }
      c++;
      if(c == noSliceCells) {
        break;
      }
    }

    // partitionCellsWorkload
    partitionCellsWorkload[i] = partitionCellNoOffsprings;
    // totalWorkload
    totalWorkload += partitionCellsWorkload[i];
    // maxNoOfssprings
    maxNoOffsprings = mMax(partitionCellNoOffsprings, maxNoOffsprings);
    // maxWorkload
    maxWorkload = mMax(partitionCellsWorkload[i], maxWorkload);
  }

  // Get total number of partition level ancestors
  MLong totalNoPartitionLevelAncestors = 0;
  MLong localPartitionLevelAncestorCount = (signed)partitionLevelAncestorIds.size();
  MPI_Allreduce(&localPartitionLevelAncestorCount, &totalNoPartitionLevelAncestors, 1, MPI_LONG, MPI_SUM, mpiCommSlice,
                AT_, "localPartitionLevelAncestorCount", "totalNoPartitionLevelAncestors");
  // max partition level shift
  MInt maxPartitionLevelShift = -1;
  MPI_Allreduce(&partitionLevelShift, &maxPartitionLevelShift, 1, MPI_INT, MPI_MAX, mpiCommSlice, AT_,
                "partitionLevelShift", "maxPartitionLevelShift");

  // Get max number of offsprings and calculate maxCPU
  MPI_Allreduce(MPI_IN_PLACE, &maxNoOffsprings, 1, MPI_INT, MPI_MAX, mpiCommSlice, AT_, "MPI_IN_PLACE",
                "maxNoOffsprings");
  // Get total workload by sum
  MPI_Allreduce(MPI_IN_PLACE, &totalWorkload, 1, MPI_DOUBLE, MPI_SUM, mpiCommSlice, AT_, "MPI_IN_PLACE",
                "totalWorkload");
  MPI_Allreduce(MPI_IN_PLACE, &maxWorkload, 1, MPI_DOUBLE, MPI_MAX, mpiCommSlice, AT_, "MPI_IN_PLACE", "maxWorkload");
  MFloat maxNoCPUs = totalWorkload / maxWorkload;


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3.4: determine min level cell information
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Mapping of child positions from 3D to 2D
  array<MInt, 8> childArray{};
  if(axis == "x") {
    childArray = {{0, 0, 1, 1, 2, 2, 3, 3}};
  } else if(axis == "y") {
    childArray = {{0, 1, 0, 1, 2, 3, 2, 3}};
  } else if(axis == "z") {
    childArray = {{0, 1, 2, 3, 0, 1, 2, 3}};
  }

  // Get boundingBox and decisiveDirection for the slice
  const array<MFloat, nDimSlice* 2> boundingBox = {{m_boundingBox[coordArray[0]], m_boundingBox[coordArray[1]],
                                                    m_boundingBox[coordArray[0] + 3],
                                                    m_boundingBox[coordArray[1] + 3]}};
  const array<MFloat, nDimSlice* 2> targetGridBoundingBox = {
      {m_targetGridBoundingBox[coordArray[0]], m_targetGridBoundingBox[coordArray[1]],
       m_targetGridBoundingBox[coordArray[0] + 3], m_targetGridBoundingBox[coordArray[1] + 3]}};
  array<MFloat, nDimSlice * 2> geometryExtents{};
  MInt decisiveDirection = 0;
  for(MInt dir = 0; dir < nDimSlice; dir++) {
    geometryExtents[dir] = targetGridBoundingBox[dir + nDimSlice] - targetGridBoundingBox[dir];
    decisiveDirection = geometryExtents[dir] > geometryExtents[decisiveDirection] ? dir : decisiveDirection;
  }

  // minLevelCellsTreeId
  MLongScratchSpace sliceMinLevelCellsTreeId(noSliceMinLevelCells, AT_, "sliceMinLevelCellsTreeId");
  for(MInt i = 0; i < noSliceMinLevelCells; ++i) {
    std::array<MFloat, nDimSlice> coord = {
        {a_coordinate(sliceMinLevelCells[i], coordArray[0]), a_coordinate(sliceMinLevelCells[i], coordArray[1])}};
    maia::grid::hilbert::coordinatesToTreeId<nDimSlice>(sliceMinLevelCellsTreeId[i], &coord[0],
                                                        (MLong)m_targetGridMinLevel, &targetGridCenterOfGravity[0],
                                                        m_targetGridLengthLevel0);
  }

  // minLevelCellsNghbrIds
  MLongScratchSpace sliceMinLevelCellsNghbrIds(noSliceMinLevelCells, 2 * nDimSlice, AT_, "sliceMinLevelCellsNghbrIds");
  for(MInt i = 0; i < noSliceMinLevelCells; ++i) {
    for(MInt j = 0; j < (nDimSlice * 2); j++) {
      MInt nghbrId = a_neighborId(sliceMinLevelCells[i], nghbrArray[j]);
      if(nghbrId == -1) {
        sliceMinLevelCellsNghbrIds(i, j) = -1;
      } else {
        sliceMinLevelCellsNghbrIds(i, j) = idMap[a_globalId(nghbrId)];
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3.5: determine cell information
  //////////////////////////////////////////////////////////////////////////////////////////////////

  MInt maxLevel = -1;
  ScratchSpace<MUchar> sliceCellInfo(noSliceCells, AT_, "sliceCellInfo");
  // cell info
  for(MInt i = 0; i < noSliceCells; ++i) {
    MLong cellId = cellIdsInSlice[i];
    MUint noChilds = 0;
    MUint curChildPos = 0;
    // Counting how many childs of current Cell belong to the slice
    for(MInt j = 0; j < pow(2, nDim) && curChildPos < noChildInRef; j++) {
      if(a_childId(cellIdsInSlice[i], j) != -1) {
        // Check to ensure that child cellId is also in slice
        if(fabs(intercept - a_coordinate(a_childId(cellIdsInSlice[i], j), axisNum))
           > halfCellLength(a_childId(cellIdsInSlice[i], j))) {
          continue;
        }
        noChilds += 1;
        curChildPos++;
      }
    }
    // Position of current cell in the child array of its parent
    MUint pos = 0;
    if(a_parentId(cellId) > -1) {
      MInt parentId = a_parentId(cellId);
      for(MUint j = 0; j < (unsigned)m_maxNoChilds; j++) {
        if(a_childId(parentId, j) == cellId) {
          pos = childArray[j];
        }
      }
    }
    // MinLevelCell indicator
    MUint isMinLvl = 0;
    if(a_level(cellId) == m_minLevel) {
      isMinLvl = (MUint)1;
    }
    MUint tmpBit = noChilds | (pos << 4) | (isMinLvl << 7);
    sliceCellInfo[i] = static_cast<MUchar>(tmpBit);

    // Determine max level in slice
    if(a_level(cellId) > maxLevel) {
      maxLevel = a_level(cellId);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 4: write grid file with slice data
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Create slice grid flle
  using namespace maia::parallel_io;

  // Create File
  const MString gridFileName = fileName + ParallelIo::fileExt();
  ParallelIo file(gridFileName, PIO_REPLACE, mpiCommSlice);

  // Determine offsets and total number of cells
  ParallelIo::size_type sliceCellsOffset, noTotalCells;
  ParallelIo::calcOffset(noSliceCells, &sliceCellsOffset, &noTotalCells, mpiCommSlice);
  // If only none partitionCells are found, set noSlicePartitionCells to 0, to write no data in
  // partitionCell array
  if(onlySliceCells) {
    noSlicePartitionCells = 0;
  }
  // Determine offsets and total number of partition cells
  ParallelIo::size_type partitionOffset, noPartitionCells;
  ParallelIo::calcOffset(noSlicePartitionCells, &partitionOffset, &noPartitionCells, mpiCommSlice);

  ParallelIo::size_type minLevelCellsOffset, noTotalMinLevelCells;
  ParallelIo::calcOffset(noSliceMinLevelCells, &minLevelCellsOffset, &noTotalMinLevelCells, mpiCommSlice);

  MFloat avgWorkload = totalWorkload / ((MFloat)noPartitionCells);
  MFloat avgOffspring = ((MFloat)noTotalCells) / ((MFloat)noPartitionCells);

  // Define arrays
  file.defineArray(PIO_LONG, "partitionCellsGlobalId", noPartitionCells);
  file.defineArray(PIO_FLOAT, "partitionCellsWorkload", noPartitionCells);
  file.defineArray(PIO_LONG, "minLevelCellsTreeId", noTotalMinLevelCells);
  file.defineArray(PIO_LONG, "minLevelCellsNghbrIds", 2 * nDimSlice * noTotalMinLevelCells);
  file.defineArray(PIO_UCHAR, "cellInfo", noTotalCells);

  const MInt noSolvers = m_tree.noSolvers();
  const MBool writeSolver = (noSolvers > 1 || g_multiSolverGrid);
  ScratchSpace<MUchar> solverBits(max(noSliceCells, 1), AT_, "solverBits");
  if(writeSolver) {
    file.defineArray(PIO_UCHAR, "solver", noTotalCells);

    for(MInt i = 0; i < noSliceCells; ++i) {
      const MLong cellId = cellIdsInSlice[i];
      MUint tmpBit = 0;
      for(MInt solver = 0; solver < noSolvers; solver++) {
        if(m_tree.solver(cellId, solver)) {
          tmpBit |= (1 << solver);
        }
      }
      solverBits[i] = static_cast<MUchar>(tmpBit);
    }
  }

  // Set attributes
  MInt tstep = 0;
  MPI_Allreduce(MPI_IN_PLACE, &maxLevel, 1, MPI_INT, MPI_MAX, mpiCommSlice, AT_, "MPI_IN_PLACE", "maxLevel");
  // Get total number of leaf cells
  MPI_Allreduce(MPI_IN_PLACE, &noLeafCells, 1, MPI_LONG, MPI_SUM, mpiCommSlice, AT_, "MPI_IN_PLACE", "noLeafCells");

  // Count number of cells per solver and set as attribute
  if(g_multiSolverGrid) {
    for(MInt b = 0; b < noSolvers; b++) {
      MLong solverCount = 0;
      for(MInt i = 0; i < noSliceCells; i++) {
        const MLong cellId = cellIdsInSlice[i];
        if(m_tree.solver(cellId, b) == true) {
          solverCount++;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &solverCount, 1, MPI_LONG, MPI_SUM, mpiCommSlice, AT_, "MPI_IN_PLACE", "solverCount");
      file.setAttributes(&solverCount, "noCells_" + std::to_string(b), 1);
    }
  }

  // Set all attributes
  file.setAttributes(&nDimSlice, "nDim", 1);
  file.setAttributes(&noSolvers, "noSolvers", 1);
  file.setAttributes(&tstep, "globalTimeStep", 1);
  file.setAttributes(&noTotalCells, "noCells", 1);
  file.setAttributes(&noLeafCells, "noLeafCells", 1);
  file.setAttributes(&noTotalMinLevelCells, "noMinLevelCells", 1);
  file.setAttributes(&noPartitionCells, "noPartitionCells", 1);
  file.setAttributes(&totalNoPartitionLevelAncestors, "noPartitionLevelAncestors", 1);
  file.setAttributes(&m_minLevel, "minLevel", 1);
  file.setAttributes(&maxLevel, "maxLevel", 1);
  file.setAttributes(&m_maxUniformRefinementLevel, "maxUniformRefinementLevel", 1);
  file.setAttributes(&maxPartitionLevelShift, "maxPartitionLevelShift", 1);
  file.setAttributes(&m_lengthLevel0, "lengthLevel0", 1);
  file.setAttributes(&centerOfGravity[0], "centerOfGravity", nDimSlice);
  file.setAttributes(&boundingBox[0], "boundingBox", 2 * nDimSlice);

  // Add additional multisolver information if grid is reordered by a different Hilbert curve
  if(m_hasMultiSolverBoundingBox) {
    file.setAttributes(&m_targetGridLengthLevel0, "multiSolverLengthLevel0", 1);
    file.setAttributes(&m_targetGridMinLevel, "multiSolverMinLevel", 1);
    file.setAttributes(&targetGridCenterOfGravity[0], "multiSolverCenterOfGravity", nDimSlice);
    file.setAttributes(&targetGridBoundingBox[0], "multiSolverBoundingBox", 2 * nDimSlice);
  }

  file.setAttributes(&m_reductionFactor, "reductionFactor", 1);
  file.setAttributes(&decisiveDirection, "decisiveDirection", 1);
  file.setAttributes(&totalWorkload, "totalWorkload", 1);
  file.setAttributes(&maxWorkload, "partitionCellMaxWorkload", 1);
  file.setAttributes(&avgWorkload, "partitionCellAverageWorkload", 1);
  file.setAttributes(&maxNoOffsprings, "partitionCellMaxNoOffspring", 1);
  file.setAttributes(&avgOffspring, "partitionCellAverageNoOffspring", 1);
  file.setAttributes(&m_partitionCellWorkloadThreshold, "partitionCellWorkloadThreshold", 1);
  file.setAttributes(&m_partitionCellOffspringThreshold, "partitionCellOffspringThreshold", 1);
  file.setAttributes(&maxNoCPUs, "maxNoBalancedCPUs", 1);

  // Save slice axis and intercept to identify the grid file
  file.setAttribute(axis, "sliceAxis");
  file.setAttribute(intercept, "sliceIntercept");

  MBool optimizedSliceIo = true;
  optimizedSliceIo = Context::getBasicProperty<MBool>("optimizedSliceIo", AT_, &optimizedSliceIo);

  // Write data in chunks of contiguous hilbert-ids/min-cells
  if(optimizedSliceIo) {
    std::vector<MInt> contHilbertIdsPartitionCounts{};
    std::vector<MInt> contHilbertIdsPartitionOffset{};

    std::vector<MInt> contHilbertIdCount{};
    std::vector<MInt> contHilbertIdOffset{};

    std::vector<MInt> contHilbertIdMinCellCount{};
    std::vector<MInt> contHilbertIdMinCellOffset{};

    { // Min-cells
      MInt contHIdMinCellCount = -1;
      std::vector<MInt> minCellCount{};
      std::vector<MInt> minCellOffset{};
      // Store min-cell counts/offsets
      for(MInt i = 0, localMinOffset = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
        const MInt count = hilbertIdMinLevelCellsCount[hIdKeys[i]];
        const MInt offset = hilbertMinLevelDomainOffset[i];
        if(localMinOffset < noSliceMinLevelCells) {
          if(count > 0) {
            minCellCount.push_back(count);
            minCellOffset.push_back(offset);
          }
        }
        localMinOffset += count;
      }

      contHIdMinCellCount = (noSliceMinLevelCells > 0) ? minCellCount[0] : -1;
      if(noSliceMinLevelCells == 1) { // Only one min-cell
        contHilbertIdMinCellCount.push_back(contHIdMinCellCount);
        contHilbertIdMinCellOffset.push_back(minCellOffset[0]);
      }

      // Find contiguous min-cells and store counts/offsets for writing
      for(MInt h = 1, i = 0; h < noSliceMinLevelCells; h++) {
        if(minCellCount[h - 1] + minCellOffset[h - 1] == minCellOffset[h]) {
          contHIdMinCellCount += minCellCount[h];
        } else {
          contHilbertIdMinCellCount.push_back(contHIdMinCellCount);
          contHilbertIdMinCellOffset.push_back(minCellOffset[i]);

          i = h;
          contHIdMinCellCount = minCellCount[h];
        }

        // Last index
        if(h == noSliceMinLevelCells - 1) {
          contHilbertIdMinCellCount.push_back(contHIdMinCellCount);
          contHilbertIdMinCellOffset.push_back(minCellOffset[i]);
        }
      }
    }

    { // Cells and partition cells
      MInt contHIdPartitionCellCount = hilbertIdPartitionCellsCount[hIdKeys[0]];
      MInt contHIdCount = hilbertIdCount[hIdKeys[0]];

      if(noLocalHilbertIds[sliceDomain] == 1) { // Only one hilbert id
        contHilbertIdsPartitionCounts.push_back(contHIdPartitionCellCount);
        contHilbertIdsPartitionOffset.push_back(hilbertPartitionDomainOffset[0]);

        contHilbertIdCount.push_back(contHIdCount);
        contHilbertIdOffset.push_back(hilbertDomainOffset[0]);
      }

      // Find contiguous cells/partition-cells and store counts/offsets for writing
      for(MInt h = 1, i = 0; h < noLocalHilbertIds[sliceDomain]; h++) {
        // Check if current hilbertId is contiguous; if so increase the current cell counts
        if(hilbertIdPartitionCellsCount[hIdKeys[h - 1]] + hilbertPartitionDomainOffset[h - 1]
           == hilbertPartitionDomainOffset[h]) {
          TERMM_IF_NOT_COND(hilbertIdCount[hIdKeys[h - 1]] + hilbertDomainOffset[h - 1] == hilbertDomainOffset[h],
                            "Error: cells not contiguous");
          contHIdPartitionCellCount += hilbertIdPartitionCellsCount[hIdKeys[h]];
          contHIdCount += hilbertIdCount[hIdKeys[h]];
        } else {
          // Not contiguous: store cell counts/offsets and start over with current cell
          contHilbertIdsPartitionCounts.push_back(contHIdPartitionCellCount);
          contHilbertIdsPartitionOffset.push_back(hilbertPartitionDomainOffset[i]);

          contHilbertIdCount.push_back(contHIdCount);
          contHilbertIdOffset.push_back(hilbertDomainOffset[i]);

          i = h;
          contHIdPartitionCellCount = hilbertIdPartitionCellsCount[hIdKeys[h]];
          contHIdCount = hilbertIdCount[hIdKeys[h]];
        }

        // Last hilbert id, store final counts/offsets
        if(h == noLocalHilbertIds[sliceDomain] - 1) {
          contHilbertIdsPartitionCounts.push_back(contHIdPartitionCellCount);
          contHilbertIdsPartitionOffset.push_back(hilbertPartitionDomainOffset[i]);

          contHilbertIdCount.push_back(contHIdCount);
          contHilbertIdOffset.push_back(hilbertDomainOffset[i]);
        }
      }
    }
    const MInt noContHilbertIds = contHilbertIdsPartitionCounts.size();

    if(solverId > -1 && noSliceContHilbertIds != nullptr) {
      TERMM_IF_COND(sliceContiguousHilbertInfo == nullptr,
                    "Error: sliceContiguousHilbertInfo is a nullptr but noSliceContHilbertIds is not.");

      std::vector<MInt> hIdCountSolver;
      std::vector<MInt> hDomainOffsetSolver;
      for(MInt i = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
        if(hilbertIdCountSolver[hIdKeys[i]] > 0) {
          hIdCountSolver.push_back(hilbertIdCountSolver[hIdKeys[i]]);
          hDomainOffsetSolver.push_back(hilbertDomainOffsetSolver[i]);
        }
      }

      std::vector<MInt> contHilbertIdCountSolver{};
      std::vector<MInt> contHilbertIdOffsetSolver{};

      // Note: can be empty if domain has no cells of this solver
      MInt contHIdCountSolver = (hIdCountSolver.size() > 0) ? hIdCountSolver[0] : 0;

      if(noLocalHilbertIdsSolver == 1) { // Only one hilbert id
        contHilbertIdCountSolver.push_back(contHIdCountSolver);
        contHilbertIdOffsetSolver.push_back(hDomainOffsetSolver[0]);
      }

      // Find contiguous cells and store counts/offsets for writing
      for(MInt h = 1, i = 0; h < noLocalHilbertIdsSolver; h++) {
        // Check if current hilbertId is contiguous; if so increase the current cell counts
        if(hIdCountSolver[h - 1] + hDomainOffsetSolver[h - 1] == hDomainOffsetSolver[h]) {
          contHIdCountSolver += hIdCountSolver[h];
        } else {
          // Not contiguous: store cell counts/offsets and start over with current cell
          contHilbertIdCountSolver.push_back(contHIdCountSolver);
          contHilbertIdOffsetSolver.push_back(hDomainOffsetSolver[i]);

          i = h;
          contHIdCountSolver = hIdCountSolver[h];
        }

        // Last hilbert id, store final count/offset
        if(h == noLocalHilbertIdsSolver - 1) {
          contHilbertIdCountSolver.push_back(contHIdCountSolver);
          contHilbertIdOffsetSolver.push_back(hDomainOffsetSolver[i]);
        }
      }

      const MInt noContHilbertIdsSolver = contHilbertIdCountSolver.size();
      *noSliceContHilbertIds = noContHilbertIdsSolver;

      if(noContHilbertIdsSolver > 0) {
        MIntScratchSpace contHilbertInfo(noContHilbertIdsSolver * 3, AT_, "contHilbertInfo");
        for(MInt i = 0; i < noContHilbertIdsSolver; i++) {
          contHilbertInfo[i * 3] = -1; // Note: cell id not required at the moment
          contHilbertInfo[i * 3 + 1] = contHilbertIdCountSolver[i];
          contHilbertInfo[i * 3 + 2] = contHilbertIdOffsetSolver[i];
        }
        std::copy_n(&contHilbertInfo[0], noContHilbertIdsSolver * 3, sliceContiguousHilbertInfo);
      }
    } else if(solverId < 0) {
      // return number of contiguous hilbert id chunks
      if(noSliceContHilbertIds != nullptr) {
        *noSliceContHilbertIds = noContHilbertIds;
      }
      // return contiguous HilbertInfo (how many cells per chunks of HilbertIds and offset for file writing)
      if(sliceContiguousHilbertInfo != nullptr) {
        MIntScratchSpace contHilbertInfo(noContHilbertIds * 3, AT_, "contHilbertInfo");
        for(MInt i = 0; i < noContHilbertIds; i++) {
          contHilbertInfo[i * 3] = -1; // Note: cell id not required at the moment
          contHilbertInfo[i * 3 + 1] = contHilbertIdCount[i];
          contHilbertInfo[i * 3 + 2] = contHilbertIdOffset[i];
        }
        std::copy_n(&contHilbertInfo[0], noContHilbertIds * 3, sliceContiguousHilbertInfo);
      }
    }

    const MFloat writeTimeStart = wallTime();
    for(MInt i = 0, localOffset = 0, localPartOffset = 0; i < hilbertDomainOffset.size0(); i++) {
      // Write data for every contiguous range of hilbertIds on this domain
      if(i < noContHilbertIds) {
        // Write partitionCells
        file.setOffset(contHilbertIdsPartitionCounts[i], contHilbertIdsPartitionOffset[i]);
        file.writeArray(&partitionCellsId[0] + localPartOffset, "partitionCellsGlobalId");
        file.writeArray(&partitionCellsWorkload[0] + localPartOffset, "partitionCellsWorkload");
        localPartOffset += contHilbertIdsPartitionCounts[i];

        // Write cellInfo
        file.setOffset(contHilbertIdCount[i], contHilbertIdOffset[i]);
        file.writeArray(&sliceCellInfo[0] + localOffset, "cellInfo");
        // Write solver bits
        if(writeSolver) {
          file.writeArray(&solverBits[0] + localOffset, "solver");
        }
        localOffset += contHilbertIdCount[i];
      } else {
        // dummy calls of writing function, if this domain has finished writing data, but other domains not
        file.setOffset(0, 0);
        file.writeArray(&partitionCellsId[0], "partitionCellsGlobalId");
        file.writeArray(&partitionCellsWorkload[0], "partitionCellsWorkload");

        file.writeArray(&sliceCellInfo[0], "cellInfo");
        if(writeSolver) {
          file.writeArray(&solverBits[0], "solver");
        }
      }
    }

    const MInt noContMinCells = contHilbertIdMinCellCount.size();
    for(MInt i = 0, localMinOffset = 0; i < hilbertDomainOffset.size0(); i++) {
      // Write data for every contiguous range of min-cells on this domain
      if(i < noContMinCells) {
        // Write min level cells tree id array
        file.setOffset(contHilbertIdMinCellCount[i], contHilbertIdMinCellOffset[i]);
        file.writeArray(&sliceMinLevelCellsTreeId[0] + localMinOffset, "minLevelCellsTreeId");
        // Write min level cells neighbor ids array
        file.setOffset(contHilbertIdMinCellCount[i] * 2 * nDimSlice, contHilbertIdMinCellOffset[i] * 2 * nDimSlice);
        file.writeArray(&sliceMinLevelCellsNghbrIds[0] + localMinOffset * 2 * nDimSlice, "minLevelCellsNghbrIds");
        localMinOffset += contHilbertIdMinCellCount[i];
      } else {
        // Dummy call if all minCell data is written or no minCell data is on this domain
        file.setOffset(0, 0);
        file.writeArray(&partitionCellsId[0], "minLevelCellsTreeId");
        file.writeArray(&partitionCellsId[0], "minLevelCellsNghbrIds");
      }
    }

    const MFloat writeTimeTotal = wallTime() - writeTimeStart;
    if(sliceDomain == 0) std::cerr << "Slice grid " << gridFileName << " write time: " << writeTimeTotal << std::endl;
  } else { // !optimizedSliceIo
    // Note: this is quite inefficient on a large number of cores or just for many min/partition cells. The optimized
    // version above should be used instead, however, for debugging/checking the output the original slow version below
    // should be kept.
    const MFloat writeTimeStart = wallTime();

    // Make sure anyone requesing the contiguous info uses the optimized version were the info is available
    if(noSliceContHilbertIds != nullptr || sliceContiguousHilbertInfo != nullptr) {
      TERMM(1, "Error: using non-optimized slice IO but contiguous Hilbert info requested by passing non-nullptr as "
               "arguments.");
    }

    // Write data in chunks which have the same hilbertId. The number of calls of the writing function
    // is equal on all domains
    for(MInt i = 0, localOffset = 0, localPartOffset = 0, localMinOffset = 0; i < hilbertDomainOffset.size0(); i++) {
      if(i < noLocalHilbertIds[sliceDomain]) {
        // Write data for every hilbertId on this domain

        // Write partitionCells
        file.setOffset(hilbertIdPartitionCellsCount[hIdKeys[i]], hilbertPartitionDomainOffset[i]);
        file.writeArray(&partitionCellsId[0] + localPartOffset, "partitionCellsGlobalId");
        file.writeArray(&partitionCellsWorkload[0] + localPartOffset, "partitionCellsWorkload");
        localPartOffset += hilbertIdPartitionCellsCount[hIdKeys[i]];

        // NOTE: cellInfo was written at last before, however this seems to occassionally result in the cellInfo array
        // not being fully written to the file (observed for the case when each slice-rank writes only a single cellInfo
        // entry)
        // Write cellInfo
        file.setOffset(hilbertIdCount[hIdKeys[i]], hilbertDomainOffset[i]);
        file.writeArray(&sliceCellInfo[0] + localOffset, "cellInfo");

        // Write solver bits
        if(writeSolver) {
          file.writeArray(&solverBits[0] + localOffset, "solver");
        }

        // Write only minCell arrays if some are found
        if(localMinOffset < noSliceMinLevelCells) {
          // Write min level cells tree id array
          file.setOffset(hilbertIdMinLevelCellsCount[hIdKeys[i]], hilbertMinLevelDomainOffset[i]);
          file.writeArray(&sliceMinLevelCellsTreeId[0] + localMinOffset, "minLevelCellsTreeId");
          // Write min level cells neighbor ids array
          file.setOffset(hilbertIdMinLevelCellsCount[hIdKeys[i]] * 2 * nDimSlice,
                         hilbertMinLevelDomainOffset[i] * 2 * nDimSlice);
          file.writeArray(&sliceMinLevelCellsNghbrIds[0] + localMinOffset * 2 * nDimSlice, "minLevelCellsNghbrIds");
        } else {
          // Dummy call if all minCell data is written or no minCell data is on this domain
          file.setOffset(0, 0);
          file.writeArray(&partitionCellsId[0], "minLevelCellsTreeId");
          file.writeArray(&partitionCellsId[0], "minLevelCellsNghbrIds");
        }
        localMinOffset += hilbertIdMinLevelCellsCount[hIdKeys[i]];

        localOffset += hilbertIdCount[hIdKeys[i]];
      } else {
        // dummy calls of writing function, if this domain has finished writing data, but other
        // domains not
        file.setOffset(0, 0);
        file.writeArray(&partitionCellsId[0], "partitionCellsGlobalId");
        file.writeArray(&partitionCellsWorkload[0], "partitionCellsWorkload");

        file.writeArray(&partitionCellsId[0], "minLevelCellsTreeId");
        file.writeArray(&partitionCellsId[0], "minLevelCellsNghbrIds");

        file.writeArray(&sliceCellInfo[0], "cellInfo");
        if(writeSolver) {
          file.writeArray(&solverBits[0], "solver");
        }
      }
    }
    const MFloat writeTimeTotal = wallTime() - writeTimeStart;
    if(sliceDomain == 0) std::cerr << "Slice grid " << gridFileName << " write time: " << writeTimeTotal << std::endl;
  }

  // Close the grid file explicitely here and finish writing before the MPI communicator is freed
  file.close();

  // Free the now unneeded MPI communicator
  MPI_Comm_free(&mpiCommSlice, AT_, "mpiCommSlice");

  m_log << "done!" << endl;
}


/// \brief Identify all leaf cells of one solver which are mapped to a given grid cell.
///
/// \author ansgar
///
/// \param[in] donorSolverId Donor solver id.
/// \param[in] gridCellId Grid cell to create mapping for.
/// \param[out] mappedLeafCells Vector of all mapped leaf cells of the donor solver.
template <MInt nDim>
void CartesianGrid<nDim>::createLeafCellMapping(const MInt donorSolverId, const MInt gridCellId,
                                                std::vector<MInt>& mappedLeafCells, MBool allChilds) {
  TRACE();

  mappedLeafCells.clear();

  std::stack<MInt> cellStack;

  // Check if this cell is used by the donor solver
  if(a_solver(gridCellId, donorSolverId)) {
    // Initialize cell stack to process
    cellStack.push(gridCellId);

    // Iterate until cellstack is empty
    // else add all children to cell stack
    while(!cellStack.empty()) {
      // Get current cell to check and remove from stack
      const MInt cellId = cellStack.top();
      cellStack.pop();

      // Check if this is a leaf cell of the donor solver
      if(a_isLeafCell(cellId, donorSolverId)) {
        // Add to mapped cells
        mappedLeafCells.push_back(cellId);
      } else {
        // Not a leaf cell, add all child cells to the cell stack and continue
        for(MInt childId = 0; childId < m_maxNoChilds; childId++) {
          if(a_hasChild(cellId, childId, donorSolverId)) {
            cellStack.push(a_childId(cellId, childId, donorSolverId));
            if(allChilds) {
              mappedLeafCells.push_back(a_childId(cellId, childId, donorSolverId));
            }
          }
        }
      }
    }
  } else {
    // Check parent cells until cell belonging to donor solver is found or there is no parent cell
    // anymore
    MBool foundMappedParent = false;
    MInt parentId = gridCellId;
    while(!foundMappedParent && a_hasParent(parentId)) {
      parentId = a_parentId(parentId);
      if(a_solver(parentId, donorSolverId) && a_isLeafCell(parentId, donorSolverId)) {
        foundMappedParent = true;
      }
    }

    if(foundMappedParent) {
      mappedLeafCells.push_back(parentId);
    }
  }

  // Sort mapped leaf cell ids
  std::sort(mappedLeafCells.begin(), mappedLeafCells.end());
}


/// \brief Determine the number of partition cells on each domain and the corresponding offsets.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[out] noPartitionCells Pointer to storage of size noDomains() that will hold the number of
///             partition cells on each domain.
/// \param[out] partitionCellOffsets Pointer to storage of size noDomains()+1 that will hold the
///             partition cell offsets for all domains and the total number of partition cells as
///             the last entry.
template <MInt nDim>
void CartesianGrid<nDim>::determineNoPartitionCellsAndOffsets(MLong* const noPartitionCells,
                                                              MLong* const partitionCellOffsets) {
  TRACE();

  const MLong localNoPartitionCells = m_noPartitionCells;
  // Gather local number of partition cells on all domains
  MPI_Allgather(&localNoPartitionCells, 1, type_traits<MLong>::mpiType(), &noPartitionCells[0], 1,
                type_traits<MLong>::mpiType(), mpiComm(), AT_, "localNoPartitionCells", "noPartitionCells[0]");

  // Partition cell offset on first domain is zero
  partitionCellOffsets[0] = 0;
  // Determine partition cell offsets for all domains (except first), last entry contains total
  // number of partition cells
  for(MInt i = 1; i < noDomains() + 1; i++) {
    partitionCellOffsets[i] = partitionCellOffsets[i - 1] + noPartitionCells[i - 1];
  }

  // Check partition cell count
  if(partitionCellOffsets[noDomains()] != m_noPartitionCellsGlobal) {
    TERMM(1, "determineNoPartitionCellsAndOffsets(): Partition cell count does not match: "
                 + std::to_string(partitionCellOffsets[noDomains()])
                 + " != " + std::to_string(m_noPartitionCellsGlobal));
  }
}


/// \brief Save given partitioning to file.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] partitionFileNameBase Output base file name.
/// \param[in] partitionCellOffset Local partition-cell offset.
template <MInt nDim>
void CartesianGrid<nDim>::savePartitionFile(const MString& partitionFileNameBase, const MLong partitionCellOffset) {
  TRACE();

  using namespace maia::parallel_io;

  stringstream partitionFileName;
  partitionFileName << m_outputDir << partitionFileNameBase << ParallelIo::fileExt();

  // Open grid partition file and write to it
  ParallelIo file(partitionFileName.str(), PIO_REPLACE, mpiComm());

  const MLong noPartitionCells = m_noPartitionCellsGlobal;
  file.setAttributes(&noPartitionCells, "noPartitionCells", 1);

  file.defineArray(PIO_LONG, "partitionCellOffset", noDomains());
  file.setOffset(1, domainId());
  file.writeArray(&partitionCellOffset, "partitionCellOffset");

  m_log << "Saved partition to file '" << partitionFileName.str() << "'." << endl;
}


/// \brief Save current grid partitioning to file.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim>
void CartesianGrid<nDim>::savePartitionFile() {
  TRACE();

  stringstream fileName;
  fileName << "partition_n" << noDomains() << "_" << globalTimeStep;

  savePartitionFile(fileName.str(), m_localPartitionCellOffsets[0]);
}


/// \brief Load a grid partitioning from a file.
template <MInt nDim>
void CartesianGrid<nDim>::loadPartitionFile(const MString& partitionFileName, MLong* partitionCellOffsets) {
  TRACE();

  using namespace maia::parallel_io;

  if(domainId() == 0) {
    // Open grid partition file
    ParallelIo file(partitionFileName, PIO_READ, MPI_COMM_SELF);

    MLong noPartitionCells = -1;
    file.getAttribute(&noPartitionCells, "noPartitionCells");
    TERMM_IF_NOT_COND(noPartitionCells == m_noPartitionCellsGlobal, "global number of partition cell mismatch");

    if(file.getArraySize("partitionCellOffset", 0) != noDomains()) {
      TERMM(1, "array size does not match number of domains");
    }
    file.setOffset(noDomains(), 0);
    file.readArray(partitionCellOffsets, "partitionCellOffset");
  }

  // Distribute partition cells ids
  MPI_Bcast(&partitionCellOffsets[0], noDomains(), type_traits<MLong>::mpiType(), 0, mpiComm(), AT_,
            "partitionCellOffsets");

  m_log << "Loaded partition from file '" << partitionFileName << "'." << endl;
}


/// \brief Update the partition cell workloads in the grid file
template <MInt nDim>
void CartesianGrid<nDim>::savePartitionCellWorkloadsGridFile() {
  TRACE();
  using namespace parallel_io;

  ParallelIo parallelIo(m_restartDir + gridInputFileName(), PIO_APPEND, mpiComm());

  // Set offset and total size for partition-cell workloads
  const MInt noLocalPartitionCells = m_noPartitionCells;
  parallelIo.setOffset(noLocalPartitionCells, m_localPartitionCellOffsets[0]);

  // Write new partition-cell workloads to grid file
  ScratchSpace<MFloat> partitionCellsWorkload(noLocalPartitionCells, AT_, "partitionCellsWorkload");
  MFloat totalWorkload = 0.0;
  for(MInt i = 0; i < noLocalPartitionCells; i++) {
    const MInt partitionCellId = m_localPartitionCellLocalIds[i];
    partitionCellsWorkload[i] = a_workload(partitionCellId);
    totalWorkload += a_workload(partitionCellId);
  }

  MPI_Allreduce(MPI_IN_PLACE, &totalWorkload, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "totalWorkload");
  // TODO labels:GRID,IO @ansgar set other attributes as well, e.g. partitionCellMaxWorkload, ...
  parallelIo.setAttributes(&totalWorkload, "totalWorkload", 1);
  parallelIo.writeArray(&partitionCellsWorkload[0], "partitionCellsWorkload");
}


template <MInt nDim>
template <MBool t_correct, MBool insideLimit>
MBool CartesianGrid<nDim>::pointWthCell(const MFloat* const coord, const MInt cellId,
                                        std::function<MFloat*(MInt, MFloat* const)> correctCellCoord) const {
  array<MFloat, nDim> tmp;
  const MFloat* tmpCoord = m_tree.coordinate(cellId);
  const MFloat* cellCenter = t_correct ? (correctCellCoord)(cellId, &tmp[0]) : tmpCoord;
  const MFloat _halfCellLength = halfCellLength(cellId);
  MBool isInCell = true;

  if(insideLimit) {
    IF_CONSTEXPR(nDim == 3) {
      return fabs(cellCenter[0] - coord[0]) <= _halfCellLength && fabs(cellCenter[1] - coord[1]) <= _halfCellLength
             && fabs(cellCenter[2] - coord[2]) <= _halfCellLength;
    }
    return fabs(cellCenter[0] - coord[0]) <= _halfCellLength && fabs(cellCenter[1] - coord[1]) <= _halfCellLength;
  } else {
    for(MInt dimId = 0; dimId < nDim; dimId++) {
      if(coord[dimId] < cellCenter[dimId] - _halfCellLength || coord[dimId] >= cellCenter[dimId] + _halfCellLength) {
        isInCell = false;
        break;
      }
    }
  }

  return isInCell;
}

/// Compute the local bounding box of all internal cells
template <MInt nDim>
void CartesianGrid<nDim>::computeLocalBoundingBox(MFloat* const bbox) {
  TRACE();

  for(MInt i = 0; i < nDim; i++) {
    bbox[i] = numeric_limits<MFloat>::max();
    bbox[nDim + i] = numeric_limits<MFloat>::lowest();
  }

  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_hasProperty(cellId, Cell::IsHalo)) {
      continue;
    }

    const MFloat halfCellLength = F1B2 * cellLengthAtLevel(a_level(cellId));
    for(MInt i = 0; i < nDim; i++) {
      bbox[i] = mMin(bbox[i], a_coordinate(cellId, i) - halfCellLength);
      bbox[nDim + i] = mMax(bbox[nDim + i], a_coordinate(cellId, i) + halfCellLength);
    }
  }
  // See createMinLevelExchangeCells()
  const MFloat halfLength0 = 0.5 * ((F1 + F1 / FPOW2(30)) * m_lengthLevel0);
  for(MInt i = 0; i < nDim; i++) {
    ASSERT(bbox[i] > m_centerOfGravity[i] - halfLength0 && bbox[nDim + i] < m_centerOfGravity[i] + halfLength0,
           "local bounding box error");
  }

  m_localBoundingBoxInit = true;
}

/// Checks if given point is inside the local bounding box
template <MInt nDim>
MBool CartesianGrid<nDim>::pointInLocalBoundingBox(const MFloat* coord) {
  // Not initialized, fast check not possible
  if(!m_localBoundingBoxInit) {
    return true;
  }

  const MFloat eps = 1e-12;
  for(MInt i = 0; i < nDim; i++) {
    if(coord[i] < m_localBoundingBox[i] - eps || coord[i] > m_localBoundingBox[nDim + i] + eps) {
      return false;
    }
  }
  return true;
}


/// \brief Perform grid sanity checks
template <MInt nDim>
void CartesianGrid<nDim>::gridSanityChecks() {
  m_log << "Performing grid sanity checks... ";

  const MInt noCells = m_tree.size();

  // Check that number internal cells and the interval between domain offsets is consistent
  TERMM_IF_NOT_COND(m_noInternalCells == m_domainOffsets[domainId() + 1] - m_domainOffsets[domainId()],
                    "Error: number of internal cells does not match with interval between domain offsets.");

  // Partition cell checks
  {
    ScratchSpace<MBool> isPartitionCell(noCells, AT_, "isPartitionCell");
    isPartitionCell.fill(false);

    MLong lastGlobalId = -1;
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      const MInt localId = m_localPartitionCellLocalIds[i];
      const MLong globalId = m_localPartitionCellGlobalIds[i];

      // Check that partition cell property is set
      TERMM_IF_NOT_COND(a_hasProperty(localId, Cell::IsPartitionCell), "Error: Cell is not a partition cell.");
      // Check the partition cell global id
      TERMM_IF_NOT_COND(globalId == a_globalId(localId), "Error: partition cell global id mismatch.");
      // Check that partition cells are sorted by global ids
      TERMM_IF_NOT_COND(globalId > lastGlobalId, "Error: partition cells not sorted by global id.");

      isPartitionCell[localId] = true;
      lastGlobalId = a_globalId(localId);
    }

    // Make sure that all other internal cells dont have the partition cell property set
    for(MInt i = 0; i < noCells; i++) {
      if(isPartitionCell[i] || a_isToDelete(i) || a_isHalo(i)) continue;
      TERMM_IF_COND(a_hasProperty(i, Cell::IsPartitionCell),
                    "Error: cell marked as partition cell but not in partition cell list " + std::to_string(i));
    }

    TERMM_IF_NOT_COND(m_localPartitionCellOffsets[2] == m_noPartitionCellsGlobal,
                      "Error: global number of partition cells mismatch.");

    MLongScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
    MLongScratchSpace localPartitionCellOffsets(noDomains() + 1, AT_, "localPartitionCellOffsets");
    // Determine the number of partition cells on all domains and the partition cell offsets
    determineNoPartitionCellsAndOffsets(&localPartitionCellCounts[0], &localPartitionCellOffsets[0]);

    TERMM_IF_NOT_COND(m_localPartitionCellOffsets[0] == localPartitionCellOffsets[domainId()],
                      "Error: partition cell offset mismatch.");
    TERMM_IF_NOT_COND(m_localPartitionCellOffsets[1] == localPartitionCellOffsets[domainId() + 1],
                      "Error: next partition cell offset mismatch.");
  }

  // Check the global to local id mapping
  for(MInt i = 0; i < noCells; i++) {
    if(a_isToDelete(i) || a_hasProperty(i, Cell::IsPeriodic)) continue;
    TERMM_IF_NOT_COND(m_globalToLocalId[a_globalId(i)] == i,
                      "Error in global to local id mapping: l" + to_string(i) + " g" + to_string(a_globalId(i)) + " g2l"
                          + to_string(m_globalToLocalId[a_globalId(i)]));
  }

  // Check partition level ancestors
  {
    ScratchSpace<MBool> isPartLvlAncestor(noCells, AT_, "isPartLvlAncestor");
    isPartLvlAncestor.fill(false);

    MInt noLocalPartLvlAncestors = 0;
    MInt count = 0;
    for(auto& i : m_partitionLevelAncestorIds) {
      TERMM_IF_NOT_COND(a_hasProperty(i, Cell::IsPartLvlAncestor), "Error: cell is not a partition level ancestor.");
      TERMM_IF_COND(a_isToDelete(i), "Error: partition level ancestor marked for deletion.");
      isPartLvlAncestor[i] = true;

      // Check global ids of childs
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        const MInt childId = a_childId(i, child);
        const MLong childGlobalId = (childId > -1) ? a_globalId(childId) : -1;
        const MInt index = count * m_maxNoChilds + child;
        TERMM_IF_NOT_COND(m_partitionLevelAncestorChildIds[index] == childGlobalId,
                          "Error: partition level ancestor child id mismatch. "
                              + std::to_string(m_partitionLevelAncestorChildIds[index])
                              + " != " + std::to_string(childGlobalId) + "; " + std::to_string(i) + "; "
                              + std::to_string(child));
      }
      count++;

      if(!a_isHalo(i)) {
        noLocalPartLvlAncestors++;
      }
    }

    // Check the number of internal partition level ancestor cells
    TERMM_IF_NOT_COND(noLocalPartLvlAncestors == m_noPartitionLevelAncestors,
                      "Error: number of local partition level ancestors mismatch.");

    // All other internal cells shouldnt be marked
    for(MInt i = 0; i < noCells; i++) {
      if(isPartLvlAncestor[i] || a_isToDelete(i) || a_isHalo(i)) continue;
      TERMM_IF_COND(a_hasProperty(i, Cell::IsPartLvlAncestor),
                    "Error: cell marked as partition level ancestor but is not in corresponding list.");
    }

    // TODO labels:GRID check other partition level ancestor info
  }

  // Check min-level cells
  {
    ScratchSpace<MBool> isMinLevelCell(noCells, AT_, "isMinLevelCell");
    isMinLevelCell.fill(false);

    for(auto& minLevelCell : m_minLevelCells) {
      TERMM_IF_NOT_COND(a_level(minLevelCell) == m_minLevel, "Error: level of min-level cell mismatch.");
      isMinLevelCell[minLevelCell] = true;
    }

    for(MInt i = 0; i < noCells; i++) {
      if(isMinLevelCell[i] || a_isToDelete(i)) continue;
      TERMM_IF_COND(a_level(i) <= m_minLevel, "Error: non-min-level cell has level <= minLevel");
    }
  }

  m_log << "done" << std::endl;
}


/// \brief Checks consistency of window/halo cells
///
/// TODO labels:GRID copied/adapted from gridproxy, generalize to avoid duplicate code?
template <MInt nDim>
void CartesianGrid<nDim>::checkWindowHaloConsistency(const MBool fullCheck, const MString a) {
  TRACE();
  if(globalNoDomains() == 1) {
    return;
  }

  m_log << "checkWindowHaloConsistency... ";

  // WH_old
  if(m_haloMode > 0) {
    ScratchSpace<MBool> isWindowCell(m_tree.size(), AT_, "isWindowCell");
    ScratchSpace<MBool> isHaloCell(m_tree.size(), AT_, "isHaloCell");
    isWindowCell.fill(false);
    isHaloCell.fill(false);

    ScratchSpace<MPI_Request> recvRequests(std::max(1, noNeighborDomains()), AT_, "recvRequests");
    ScratchSpace<MPI_Request> sendRequests(std::max(1, noNeighborDomains()), AT_, "sendRequests");

    for(MInt solver = 0; solver < treeb().noSolvers(); ++solver) {
      //////////////////////////////////////////////////////////////////////////////////////////////////
      // Check #1: number of window/halo cells match
      //////////////////////////////////////////////////////////////////////////////////////////////////
      // Start receiving number of window cells from each neighbor domain
      fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
      MIntScratchSpace noWindowCellsRecv(std::max(1, noNeighborDomains()), AT_, "noWindowCellsRecv");
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        noWindowCellsRecv[d] = -1;
        MPI_Irecv(&noWindowCellsRecv[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), neighborDomain(d),
                  mpiComm(), &recvRequests[d], AT_, "noWindowCellsRecv[d]");
      }

      // Start sending number of window cells to each neighbor domain
      fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
      MIntScratchSpace noWindowCellsSend(std::max(1, noNeighborDomains()), AT_, "noWindowCellsSend");
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        // //TODO_SS labels:GRID,COMM currently without duplicate window cells for periodic ...
        noWindowCellsSend[d] = noSolverWindowCells(d, solver); // noWindowCells(d);
        MPI_Isend(&noWindowCellsSend[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), domainId(), mpiComm(),
                  &sendRequests[d], AT_, "noWindowCellsSend[d]");
      }

      // Finish MPI communication
      MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
      MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

      // Check if received number of window cells matches the local number of halo cells
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        if(noWindowCellsRecv[d] != noSolverHaloCells(d, solver, true) /*noHaloCells(d)*/) {
          TERMM(1, "Cartesian Grid : Number of window cells from domain " + to_string(neighborDomain(d))
                       + " does not match local number of halo cells; window: " + to_string(noWindowCellsRecv[d])
                       + " ,halo: " + to_string(noSolverHaloCells(d, solver, true) /*noHaloCells(d)*/) + " " + a
                       + " d=" + to_string(domainId()) + " solver=" + to_string(solver));
        }
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Check #2: grid global ids of window/halo cells match
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Start receiving window cell global ids from each neighbor domain
    fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
    // const MInt totalNoWindowCellsRecv = accumulate(noWindowCellsRecv.begin(), noWindowCellsRecv.end(), 0);
    MInt totalNoWindowCellsRecv = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++)
      totalNoWindowCellsRecv += noHaloCells(d);

    MIntScratchSpace windowCellsRecv(max(1, totalNoWindowCellsRecv), AT_, "windowCellsRecv");
    fill(windowCellsRecv.begin(), windowCellsRecv.end(), -1);
    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      if(noHaloCells(d) > 0) {
        MPI_Irecv(&windowCellsRecv[offset], noHaloCells(d), type_traits<MInt>::mpiType(), neighborDomain(d),
                  neighborDomain(d), mpiComm(), &recvRequests[d], AT_, "windowCellsRecv[offset]");
      }
      offset += noHaloCells(d);
    }

    // Start sending window cell global ids to each neighbor domain
    fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
    // const MInt totalNoWindowCellsSend = accumulate(noWindowCellsSend.begin(), noWindowCellsSend.end(), 0);
    MInt totalNoWindowCellsSend = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++)
      totalNoWindowCellsSend += noWindowCells(d);
    MIntScratchSpace windowCellsSend(max(1, totalNoWindowCellsSend), AT_, "windowCellsSend");

    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noWindowCells(d) /*noWindowCellsSend[d]*/; c++) {
        //      TERMM_IF_NOT_COND(a_hasProperty(windowCell(d, c), Cell::IsWindow), a+"not a window cell");
        windowCellsSend[offset + c] = a_globalId(windowCell(d, c));
        isWindowCell(windowCell(d, c)) = true;
      }
      if(noWindowCells(d) > 0) {
        MPI_Isend(&windowCellsSend[offset], noWindowCells(d), type_traits<MInt>::mpiType(), neighborDomain(d),
                  domainId(), mpiComm(), &sendRequests[d], AT_, "windowCellsSend[offset]");
      }
      offset += noWindowCells(d);
    }

    // Finish MPI communication
    MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

    // Check if received window cell global ids match the local halo cell global ids
    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noHaloCells(d); c++) {
        const MInt cellId = haloCell(d, c);
        TERMM_IF_NOT_COND(a_isHalo(cellId), "not a halo cell");
        isHaloCell(haloCell(d, c)) = true;
        const MInt globalId = a_globalId(cellId);
        // If halo cell has periodic flag, its global id should be -1
        if(windowCellsRecv[offset + c] != globalId && !(a_hasProperty(cellId, Cell::IsPeriodic) && globalId == -1)) {
          TERMM(1, a + "Global id of window cell " + to_string(c) + " from domain " + to_string(neighborDomain(d))
                       + " does not match local halo cell gobal id (" + to_string(windowCellsRecv[offset + c]) + " vs. "
                       + to_string(a_globalId(haloCell(d, c))) + ")");
        }
      }
      offset += noHaloCells(d);
    }

    // So if-statements below do not fail
    if(m_azimuthalPer) {
      for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
        for(MInt c = 0; c < noAzimuthalWindowCells(d); c++) {
          isWindowCell(azimuthalWindowCell(d, c)) = true;
        }
        for(MInt c = 0; c < noAzimuthalHaloCells(d); c++) {
          isHaloCell(azimuthalHaloCell(d, c)) = true;
        }
      }
      for(MInt c = 0; c < noAzimuthalUnmappedHaloCells(); c++) {
        isHaloCell(azimuthalUnmappedHaloCell(c)) = true;
      }
    }

    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      // check that only cells in m_windowCells are marked as window!
      if(!isWindowCell(cellId)) {
        TERMM_IF_NOT_COND(!a_hasProperty(cellId, Cell::IsWindow),
                          "cell is marked as window but not in m_windowCells: " + std::to_string(cellId) + a);
      }
      // check that only cells in m_haloCells are marked as halo!
      if(!isHaloCell(cellId) && !a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
        TERMM_IF_NOT_COND(!a_isHalo(cellId),
                          "cell is marked as halo but not in m_haloCells: " + std::to_string(cellId));
      }
    }
  } else {
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Check #1: number of window/halo cells match
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Start receiving number of window cells from each neighbor domain
    ScratchSpace<MPI_Request> recvRequests(std::max(1, noNeighborDomains()), AT_, "recvRequests");
    MIntScratchSpace noWindowCellsRecv(std::max(1, noNeighborDomains()), AT_, "noWindowCellsRecv");
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      noWindowCellsRecv[d] = -1;
      MPI_Irecv(&noWindowCellsRecv[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), neighborDomain(d), mpiComm(),
                &recvRequests[d], AT_, "noWindowCellsRecv[d]");
    }

    // Start sending number of window cells to each neighbor domain
    ScratchSpace<MPI_Request> sendRequests(std::max(1, noNeighborDomains()), AT_, "sendRequests");
    MIntScratchSpace noWindowCellsSend(std::max(1, noNeighborDomains()), AT_, "noWindowCellsSend");
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      noWindowCellsSend[d] = noWindowCells(d);
      MPI_Isend(&noWindowCellsSend[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), domainId(), mpiComm(),
                &sendRequests[d], AT_, "noWindowCellsSend[d]");
    }

    // Finish MPI communication
    MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

    // Check if received number of window cells matches the local number of halo cells
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      if(noWindowCellsRecv[d] != noHaloCells(d)) {
        TERMM(1, "Cartesian Grid : Number of window cells from domain " + to_string(neighborDomain(d))
                     + " does not match local number of halo cells; window: " + to_string(noWindowCellsRecv[d])
                     + " ,halo: " + to_string(noHaloCells(d)));
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Check #2: grid global ids of window/halo cells match
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // Start receiving window cell global ids from each neighbor domain
    fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
    const MInt totalNoWindowCellsRecv = accumulate(noWindowCellsRecv.begin(), noWindowCellsRecv.end(), 0);

    MLongScratchSpace windowCellsRecv(max(1, totalNoWindowCellsRecv), AT_, "windowCellsRecv");
    fill(windowCellsRecv.begin(), windowCellsRecv.end(), -1);
    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      if(noHaloCells(d) > 0) {
        MPI_Irecv(&windowCellsRecv[offset], noHaloCells(d), type_traits<MLong>::mpiType(), neighborDomain(d),
                  neighborDomain(d), mpiComm(), &recvRequests[d], AT_, "windowCellsRecv[offset]");
      }
      offset += noHaloCells(d);
    }

    // Start sending window cell global ids to each neighbor domain
    fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
    const MInt totalNoWindowCellsSend = accumulate(noWindowCellsSend.begin(), noWindowCellsSend.end(), 0);
    MLongScratchSpace windowCellsSend(max(1, totalNoWindowCellsSend), AT_, "windowCellsSend");

    ScratchSpace<MBool> isWindowCell(m_tree.size(), AT_, "isWindowCell");
    ScratchSpace<MBool> isHaloCell(m_tree.size(), AT_, "isHaloCell");
    isWindowCell.fill(false);
    isHaloCell.fill(false);

    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noWindowCellsSend[d]; c++) {
        TERMM_IF_NOT_COND(a_hasProperty(windowCell(d, c), Cell::IsWindow), "not a window cell");
        windowCellsSend[offset + c] = a_globalId(windowCell(d, c));
        isWindowCell(windowCell(d, c)) = true;
      }
      if(noWindowCells(d) > 0) {
        MPI_Isend(&windowCellsSend[offset], noWindowCells(d), type_traits<MLong>::mpiType(), neighborDomain(d),
                  domainId(), mpiComm(), &sendRequests[d], AT_, "windowCellsSend[offset]");
      }
      offset += noWindowCells(d);
    }

    // Finish MPI communication
    MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
    MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

    // Check if received window cell global ids match the local halo cell global ids
    for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noHaloCells(d); c++) {
        const MInt cellId = haloCell(d, c);
        TERMM_IF_NOT_COND(a_isHalo(cellId), "not a halo cell");
        isHaloCell(haloCell(d, c)) = true;
        const MLong globalId = a_globalId(cellId);
        // If halo cell has periodic flag, its global id should be -1
        if(windowCellsRecv[offset + c] != globalId && !(a_hasProperty(cellId, Cell::IsPeriodic) && globalId == -1)) {
          TERMM(1, "Global id of window cell " + to_string(c) + " from domain " + to_string(neighborDomain(d))
                       + " does not match local halo cell gobal id (" + to_string(windowCellsRecv[offset + c]) + " vs. "
                       + to_string(a_globalId(haloCell(d, c))) + ")");
        }
      }
      offset += noHaloCells(d);
    }

    // So if-statements below do not fail
    if(m_azimuthalPer) {
      for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
        for(MInt c = 0; c < noAzimuthalWindowCells(d); c++) {
          isWindowCell(azimuthalWindowCell(d, c)) = true;
        }
        for(MInt c = 0; c < noAzimuthalHaloCells(d); c++) {
          isHaloCell(azimuthalHaloCell(d, c)) = true;
        }
      }
      for(MInt c = 0; c < noAzimuthalUnmappedHaloCells(); c++) {
        isHaloCell(azimuthalUnmappedHaloCell(c)) = true;
      }
    }

    for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
      // check that only cells in m_windowCells are marked as window!
      if(!isWindowCell(cellId)) {
        TERMM_IF_NOT_COND(!a_hasProperty(cellId, Cell::IsWindow),
                          "cell is marked as window but not in m_windowCells: " + std::to_string(cellId));
      }
      // check that only cells in m_haloCells are marked as halo!
      if(!isHaloCell(cellId) && !a_hasProperty(cellId, Cell::IsPartLvlAncestor)) {
        TERMM_IF_NOT_COND(!a_isHalo(cellId),
                          "cell is marked as halo but not in m_haloCells: " + std::to_string(cellId));
      }
    }
  }

  // Check parent/child/neighbor global ids (if != -1)
  // Note: not working if the global ids are not computed/updated yet (e.g. during adaptation)!
  // Note: assumes that internal and halo cells are separated
  if(fullCheck && noNeighborDomains() > 0) {
    ScratchSpace<MLong> parentData(m_tree.size(), AT_, "parentData");
    ScratchSpace<MLong> childData(m_tree.size(), m_maxNoChilds, AT_, "childData");
    parentData.fill(-1);
    childData.fill(-1);

    for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
      TERMM_IF_COND(a_isHalo(cellId), "cell is a halo" + a);
      const MInt parentId = a_parentId(cellId);
      parentData[cellId] = (parentId > -1) ? a_globalId(parentId) : -1;

      for(MInt j = 0; j < m_maxNoChilds; j++) {
        const MInt childId = a_childId(cellId, j);
        childData(cellId, j) = (childId > -1) ? a_globalId(childId) : -1;
      }
    }

    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), parentData.getPointer(), 1);

    maia::mpi::exchangeData(m_nghbrDomains, m_haloCells, m_windowCells, mpiComm(), childData.getPointer(),
                            m_maxNoChilds);

    for(MInt d = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noHaloCells(d); c++) {
        const MInt cellId = haloCell(d, c);
        TERMM_IF_NOT_COND(a_isHalo(cellId), "not a halo cell");
        const MInt parentId = a_parentId(cellId);
        if(parentId > -1) {
          TERMM_IF_NOT_COND(a_globalId(parentId) == parentData[cellId],
                            "Halo cell parent global id mismatch: cell" + std::to_string(cellId) + " parent"
                                + std::to_string(parentId) + "; " + std::to_string(a_globalId(parentId))
                                + " != " + std::to_string(parentData[cellId]));
        }

        for(MInt j = 0; j < m_maxNoChilds; j++) {
          const MInt childId = a_childId(cellId, j);
          if(childId > -1) {
            TERMM_IF_NOT_COND(a_globalId(childId) == childData(cellId, j),
                              "Halo cell child global id mismatch: cell" + std::to_string(cellId) + " child"
                                  + std::to_string(childId) + "; " + std::to_string(a_globalId(childId))
                                  + " != " + std::to_string(childData(cellId, j)));
          }
        }
      }
    }
  }

  m_log << "done" << std::endl;

  if(m_azimuthalPer) {
    checkAzimuthalWindowHaloConsistency();
  }
}


/// \brief Write the cell data of the local tree to some file for debugging purposes
template <MInt nDim>
void CartesianGrid<nDim>::dumpCellData(const MString name) {
  ofstream logfile;
  logfile.open(name + "_" + std::to_string(domainId()));

  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    logfile << cellId << " g" << a_globalId(cellId) << " l" << a_level(cellId) << " p" << a_parentId(cellId);
    for(MUint j = 0; j < (unsigned)m_maxNoChilds; j++) {
      logfile << " c" << a_childId(cellId, j);
    }
    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      logfile << " n" << a_neighborId(cellId, dirId);
    }
    logfile << " " << a_propertiesToString(cellId);
    logfile << " w" << a_weight(cellId);
    logfile << " o" << a_noOffsprings(cellId) << std::endl;
  }
  logfile << std::endl;

  for(MUint i = 0; i < m_minLevelCells.size(); i++) {
    logfile << "m" << i << " " << m_minLevelCells[i] << std::endl;
  }
  logfile << std::endl;

  for(MInt i = 0; i < m_noPartitionCells; i++) {
    logfile << "p" << i << " g" << m_localPartitionCellGlobalIds[i] << " l" << m_localPartitionCellLocalIds[i]
            << " level" << a_level(m_localPartitionCellLocalIds[i]) << std::endl;
  }
  logfile << std::endl;

  MInt count = 0;
  for(auto& id : m_globalToLocalId) {
    logfile << "g2l" << count++ << " g" << id.first << " l" << id.second << std::endl;
  }
  logfile << std::endl;

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    logfile << "window " << i << " " << m_nghbrDomains[i] << " " << m_windowCells[i].size() << std::endl;
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      logfile << "window " << m_nghbrDomains[i] << " " << j << " w" << m_windowCells[i][j] << " g"
              << a_globalId(m_windowCells[i][j]) << std::endl;
    }
    logfile << std::endl;

    logfile << "halo " << i << " " << m_nghbrDomains[i] << " " << m_haloCells[i].size() << std::endl;
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); j++) {
      logfile << "halo " << m_nghbrDomains[i] << " " << j << " h" << m_haloCells[i][j]

              << " g" << a_globalId(m_haloCells[i][j]) << std::endl;
    }
    logfile << std::endl;
  }
}


/// \brief Determine new partition cells (i.e. in/decrease the partition level shifts) and change
///        the grid accordingly such that these partition cells are used for load balancing
///
/// Note: should only be called from gridcontroller::balance() which restores a consistent grid with
/// e.g. all partition level ancestor information! Also updating the partition cells might lead to
/// domains which temporarily dont have any partition cell anymore!
template <MInt nDim>
MBool CartesianGrid<nDim>::updatePartitionCells(MInt offspringThreshold, MFloat workloadThreshold) {
  TRACE();
  m_log << "Update partition cells: offspringThreshold=" << offspringThreshold
        << "; workloadThreshold=" << workloadThreshold << std::endl;

  MBool partitionCellChange = false;

  // update noOffspring and workloads
  calculateNoOffspringsAndWorkload(static_cast<Collector<void>*>(nullptr), m_tree.size());

  // Determine new partition cells (similar to gridio::saveGrid())
  // 1. Start with collecting min-level cells
  std::vector<MInt> tmpMinLevelCells;
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_level(cellId) == m_minLevel && !a_hasProperty(cellId, Cell::IsPeriodic)) {
      ASSERT(a_parentId(cellId) == -1, "");
      tmpMinLevelCells.push_back(cellId);
    }
  }

  std::vector<MInt> partitionCellsFiltered;
  while(true) { // Repeat partition cell filtering until enough partition cells are present
    partitionCellsFiltered.clear();
    std::vector<MInt> tmpPartitionCells(tmpMinLevelCells);

    // Iterate over all cells to consider
    while(!tmpPartitionCells.empty()) {
      const MInt cellId = tmpPartitionCells.front();
      tmpPartitionCells.erase(tmpPartitionCells.begin());

      // Multisolver: make sure that children of solver leaf-cells cannot become partition cells!
      MBool isSolverLeafCell = false;
      for(MInt b = 0; b < treeb().noSolvers(); b++) {
        if(a_isLeafCell(cellId, b)) {
          isSolverLeafCell = true;
          break;
        }
      }

      // Check if cell has too many offspring or too high workload (and is not a leaf cell)
      if((a_noOffsprings(cellId) > offspringThreshold || a_workload(cellId) > workloadThreshold) && !isSolverLeafCell) {
        MInt cnt = 0;
        // Add child cells to list of cells to consider as partition cells
        for(MInt j = 0; j < IPOW2(nDim); ++j) {
          if(a_childId(cellId, j) > -1) {
            tmpPartitionCells.insert(tmpPartitionCells.begin() + cnt, a_childId(cellId, j));
            cnt++;
          }
        }
      } else if(!a_isHalo(cellId)) {
        // Add as partition cell
        partitionCellsFiltered.push_back(cellId);
      }
    }

    MLong globalNoNewPartitionCells = partitionCellsFiltered.size();
    MPI_Allreduce(MPI_IN_PLACE, &globalNoNewPartitionCells, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_,
                  "MPI_IN_PLACE", "globalNoNewPartitionCells");
    cerr0 << "Found " << globalNoNewPartitionCells << " new partition cells." << std::endl;
    m_log << "Found " << globalNoNewPartitionCells << " new partition cells." << std::endl;

    // Make sure that there are enough partition cells (at least one per domain)
    // TODO labels:GRID add a minimum allowed number of partition cells on average per domain
    if(globalNoNewPartitionCells < noDomains()) {
      offspringThreshold /= 2;
      workloadThreshold *= 0.5;
      cerr0 << "Error: too few partition cells; Repeating filtering of partition cells with "
               "lowered thresholds..."
            << std::endl;
    } else {
      break;
    }
  }

  const MInt oldNoPartCells = m_noPartitionCells;
  const MInt newNoPartCells = partitionCellsFiltered.size();

  sort(partitionCellsFiltered.begin(), partitionCellsFiltered.end());

  if(oldNoPartCells != newNoPartCells) {
    // Number of partition cells changed
    partitionCellChange = true;
  } else {
    // Number of partition cells did not change, check the partition cells for any change
    for(MInt i = 0; i < newNoPartCells; i++) {
      const MInt newPartCellId = partitionCellsFiltered[i];
      if(!a_hasProperty(newPartCellId, Cell::IsPartitionCell)) {
        partitionCellChange = true;
        break;
      }
    }
  }

  MInt globalPartitionCellChange = partitionCellChange;
  MPI_Allreduce(MPI_IN_PLACE, &globalPartitionCellChange, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_,
                "MPI_IN_PLACE", "globalPartitionCellChange");

  // If there is no change of any partition cell globally, nothing to do
  if(!globalPartitionCellChange) {
    cerr0 << "Partition cells did not change." << std::endl;
    m_log << "Partition cells did not change." << std::endl;
    return false;
  }

  // Set 'updating partition cells'-flag (allow case with 0 partition cells in eg. computeGlobalIds)
  m_updatingPartitionCells = true;

  // Store if the current first min-level cell is a halo partition level ancestor (i.e. is the
  // ancestor of the first local partition cell)
  const MInt firstMinLvlCell = m_minLevelCells[0];
  MInt minLvlHaloPartLvlAncestor = -1;
  if(a_hasProperty(firstMinLvlCell, Cell::IsPartLvlAncestor) && a_hasProperty(firstMinLvlCell, Cell::IsHalo)) {
    minLvlHaloPartLvlAncestor = firstMinLvlCell;
  }

  // Reset cell flags
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    a_hasProperty(cellId, Cell::IsPartitionCell) = false;
    a_hasProperty(cellId, Cell::IsPartLvlAncestor) = false;
  }

  // Set new partition cell flags
  MInt maxLevelDiff = 0;
  for(MInt i = 0; i < newNoPartCells; i++) {
    const MInt newPartCellId = partitionCellsFiltered[i];
    a_hasProperty(newPartCellId, Cell::IsPartitionCell) = true;

    maxLevelDiff = std::max(maxLevelDiff, a_level(newPartCellId) - m_minLevel);
  }

  m_noPartitionCells = newNoPartCells;

  m_noPartitionCellsGlobal = m_noPartitionCells;
  MPI_Allreduce(MPI_IN_PLACE, &m_noPartitionCellsGlobal, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_,
                "MPI_IN_PLACE", "m_noPartitionCellsGlobal");

  MPI_Allreduce(MPI_IN_PLACE, &maxLevelDiff, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxLevelDiff");
  cerr0 << "New maximum partition level shift: " << maxLevelDiff << std::endl;
  m_log << "New maximum partition level shift: " << maxLevelDiff << std::endl;
  m_maxPartitionLevelShift = maxLevelDiff;

  // Reallocate partition cell arrays with new size
  // Note: data is set in computeGlobalIds, should not be required before
  mDeallocate(m_localPartitionCellGlobalIds);
  mAlloc(m_localPartitionCellGlobalIds, std::max(m_noPartitionCells, 1), "m_localPartitionCellGlobalIds",
         static_cast<MLong>(-1), AT_);
  mDeallocate(m_localPartitionCellLocalIds);
  mAlloc(m_localPartitionCellLocalIds, std::max(m_noPartitionCells, 1), "m_localPartitionCellLocalIds", -1, AT_);

  MIntScratchSpace isPartitionLevelAncestor(m_noInternalCells, AT_, "isPartitionLevelAncestor");
  isPartitionLevelAncestor.fill(0);

  MInt noPartitionLevelAncestorsLocal = 0;
  // Determine partition level ancestor cells and set the corresponding cell property
  for(MInt i = 0; i < newNoPartCells; i++) {
    const MInt cellId = partitionCellsFiltered[i];
    MInt parentId = a_parentId(cellId);

    // Loop up the parent cells of all partition cells
    while(parentId > -1) {
      a_hasProperty(parentId, Cell::IsPartLvlAncestor) = true;

      if(!a_isHalo(parentId)) {
        noPartitionLevelAncestorsLocal++;
      }

      for(MInt b = 0; b < treeb().noSolvers(); b++) {
        TERMM_IF_COND(a_isLeafCell(parentId, b), "Error: partition level ancestor is the leaf cell of a solver.");
      }

      parentId = a_parentId(parentId);
    }
  }

  m_noPartitionLevelAncestors = noPartitionLevelAncestorsLocal;
  // Determine global number of partition level ancestor cells
  m_noPartitionLevelAncestorsGlobal = m_noPartitionLevelAncestors;
  MPI_Allreduce(MPI_IN_PLACE, &m_noPartitionLevelAncestorsGlobal, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(),
                AT_, "MPI_IN_PLACE", "noPartitionLevelAncestorsGlobal");

  std::vector<MInt>().swap(m_partitionLevelAncestorIds);
  std::vector<MLong>().swap(m_partitionLevelAncestorChildIds);

  // Note: only the min-level halo partition level ancestor needs to be stored since it is required
  // for computeGlobalIds()
  if(minLvlHaloPartLvlAncestor > -1) {
    TERMM_IF_NOT_COND(a_isHalo(minLvlHaloPartLvlAncestor), "Error: not a halo cell.");
    a_hasProperty(minLvlHaloPartLvlAncestor, Cell::IsPartLvlAncestor) = true;
    m_partitionLevelAncestorIds.push_back(minLvlHaloPartLvlAncestor);
  }

  // Update partition cell local/global ids and partition cell offsets
  updatePartitionCellInformation();

  // exchange properties
  exchangeProperties();

  m_updatingPartitionCells = false;
  // Set status, checked when writing a grid restart file to avoid a refiltering of partition cells
  m_updatedPartitionCells = true;

  return true;
}


/// \brief Store partition level ancestor window/halo cells
template <MInt nDim>
void CartesianGrid<nDim>::determinePartLvlAncestorHaloWindowCells(
    std::vector<std::vector<MInt>>& partLvlAncestorHaloCells,
    std::vector<std::vector<MInt>>& partLvlAncestorWindowCells) {
  TRACE();

  partLvlAncestorHaloCells.clear();
  partLvlAncestorWindowCells.clear();

  partLvlAncestorHaloCells.resize(noNeighborDomains());
  partLvlAncestorWindowCells.resize(noNeighborDomains());

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt j = 0; j < (signed)m_haloCells[i].size(); ++j) {
      if(a_hasProperty(m_haloCells[i][j], Cell::IsPartLvlAncestor)) {
        partLvlAncestorHaloCells[i].push_back(m_haloCells[i][j]);
      }
    }
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); ++j) {
      if(a_hasProperty(m_windowCells[i][j], Cell::IsPartLvlAncestor)) {
        partLvlAncestorWindowCells[i].push_back(m_windowCells[i][j]);
      }
    }
  }
}

// Sets the solver flag for cells which have been newly created in the grid
template <MInt nDim>
void CartesianGrid<nDim>::setChildSolverFlag(
    const MInt cellId, const MInt solver,
    const std::function<MInt(const MFloat*, const MInt, const MInt)>& cellOutsideSolver) {
  TRACE();

  static constexpr MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                               {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};

  const MInt childLevel = a_level(cellId) + 1;
  const MFloat childCellLength = cellLengthAtLevel(childLevel);

  for(MInt c = 0; c < m_maxNoChilds; c++) {
    MInt childId = a_childId(cellId, c);

    if(childId < 0) continue;
    if(treeb().solver(childId, solver)) continue;

    MFloat coords[nDim];
    for(MInt i = 0; i < nDim; i++) {
      coords[i] = a_coordinate(cellId, i) + F1B2 * signStencil[c][i] * childCellLength;
    }
    MInt isOutside = cellOutsideSolver(coords, childLevel, cellId);

    if(isOutside < 1) {
      treeb().solver(childId, solver) = true;
    }
  }
}

/** \brief Transform cartesian cell coordinate to cylindrcal coordinats using the
    / using the azimuthal periodic center
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::cartesianToCylindric(const MFloat* coords, MFloat* coordsCyl) {
  TRACE();

  MFloat radius = F0;

  for(MInt d = 0; d < nDim; d++) {
    if(m_periodicCartesianDir[d] == 0) {
      coordsCyl[nDim - 1] = coords[d];
    } else {
      radius += pow((coords[d] - m_azimuthalPerCenter[d]), 2);
    }
  }
  radius = sqrt(radius);

  MFloat fac = F0;
  if(coords[m_azimuthalPeriodicDir[0]] >= F0) {
    fac = 1.0;
  } else {
    fac = -1.0;
  }
  MFloat phi = fac * acos(coords[m_azimuthalPeriodicDir[1]] / radius);

  coordsCyl[0] = radius;
  coordsCyl[1] = phi;
}

/** \brief Rotate caresian coordinates by angle. Azimuthal periodic center is used.
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::rotateCartesianCoordinates(MFloat* coords, MFloat angle) {
  TRACE();

  MFloat tmpCoords[nDim];

  for(MInt d = 0; d < nDim; d++) {
    tmpCoords[d] = coords[d] - m_azimuthalPerCenter[d];
  }

  coords[m_azimuthalPeriodicDir[0]] =
      tmpCoords[m_azimuthalPeriodicDir[0]] * cos(angle) - tmpCoords[m_azimuthalPeriodicDir[1]] * sin(angle);
  coords[m_azimuthalPeriodicDir[1]] =
      tmpCoords[m_azimuthalPeriodicDir[0]] * sin(angle) + tmpCoords[m_azimuthalPeriodicDir[1]] * cos(angle);

  for(MInt d = 0; d < nDim; d++) {
    coords[d] += m_azimuthalPerCenter[d];
  }
}


/** \brief Tag azimuthal halo cells for refinement
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::tagAzimuthalHigherLevelExchangeCells(vector<MLong>& refineChildIds,
                                                               vector<MLong>& recvChildIds,
                                                               MInt level) {
  TRACE();

  if(noAzimuthalNeighborDomains() == 0) return;

  static constexpr MFloat cornerStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                                 {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
  MIntScratchSpace nghbrList(200, AT_, "nghbrList");
  MFloat childCoords[3] = {F0, F0, F0};
  MFloat coordsCyl[3] = {F0, F0, F0};
  MFloat dxCyl = m_azimuthalAngle;

  MInt noCells = treeb().size();
  MBoolScratchSpace gridBndryCells(noCells, AT_, "gridBndryCells");
  gridBndryCells.fill(false);
  for(auto it = m_gridBndryCells.begin(); it != m_gridBndryCells.end(); it++) {
    MInt cellId = *it;
    gridBndryCells[cellId] = true;
  }

  MInt windowCnt = 0;
  MInt haloCnt = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    windowCnt += m_azimuthalWindowCells[i].size();
    haloCnt += m_azimuthalHaloCells[i].size();
  }

  refineChildIds.resize(windowCnt * m_maxNoChilds);
  recvChildIds.resize(haloCnt * m_maxNoChilds);

  MIntScratchSpace refinedWindows(windowCnt, AT_, "refinedWindows");
  refinedWindows.fill(0);
  MIntScratchSpace refineCandidates(haloCnt, AT_, "refineCandidates");
  refineCandidates.fill(0);

  // Check which azimuthal window allows possible halo refinement
  // Only halo cells for which the corresponding internal cell is refined
  // are allowed to be refined
  MInt cnt = 0;
  MIntScratchSpace windowRefineCnt(noAzimuthalNeighborDomains(), AT_, "windowRefineCnt");
  windowRefineCnt.fill(0);
  MInt windowRefineCntTotal = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalWindowCells(i); j++) {
      MInt cellId = m_azimuthalWindowCells[i][j];

      MBool refine = false;
      if(a_level(cellId) == level) {
        if(a_noChildren(cellId) > 0) {
          refine = true;
        }
      }

      if(refine) {
        refinedWindows[cnt] = 1;
        windowRefineCnt[i]++;
        windowRefineCntTotal++;
      }
      cnt++;
    }
  }

  MIntScratchSpace windowMap(windowRefineCntTotal, AT_, "windowMap");
  cnt = 0;
  windowRefineCntTotal = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalWindowCells(i); j++) {
      if(refinedWindows[cnt] > 0) {
        MInt cellId = m_azimuthalWindowCells[i][j];
        windowMap[windowRefineCntTotal] = cellId;
        windowRefineCntTotal++;
      }
      cnt++;
    }
  }

  maia::mpi::exchangeBuffer(m_azimuthalNghbrDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                            refinedWindows.getPointer(), refineCandidates.getPointer());

  // Compute hilbertIds of refineCandidates childs
  cnt = 0;
  MInt haloRefineCntTotal = 0;
  MIntScratchSpace haloRefineCnt(noAzimuthalNeighborDomains(), AT_, "haloRefineCnt");
  haloRefineCnt.fill(0);
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalHaloCells(i); j++) {
      if(refineCandidates[cnt] > 0) {
        haloRefineCnt[i]++;
        haloRefineCntTotal++;
      }
      cnt++;
    }
  }

  MIntScratchSpace haloMap(haloRefineCntTotal, AT_, "haloMap");
  MIntScratchSpace refineHalo(noCells, AT_, "refineHalo");
  cnt = 0;
  haloRefineCntTotal = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalHaloCells(i); j++) {
      MInt cellId = m_azimuthalHaloCells[i][j];
      refineHalo[cellId] = refineCandidates[cnt];
      if(refineCandidates[cnt] > 0) {
        haloMap[haloRefineCntTotal] = cellId;
        haloRefineCntTotal++;
      }
      cnt++;
    }
  }

  MIntScratchSpace layerId(noCells, AT_, "layerId");
  layerId.fill(-1);
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_isHalo(cellId)) {
      continue;
    }
    if(a_level(cellId) != level) {
      continue;
    }
    for(MInt d = 0; d < m_noDirs; d++) {
      MInt nghbrId = a_neighborId(cellId, d);
      if(nghbrId <= -1) {
        continue;
      }
      if(!a_isHalo(nghbrId)) {
        continue;
      }
      layerId[nghbrId] = 1;

      MInt addLayer = true;
      if(!a_hasProperty(nghbrId, Cell::IsPeriodic)) {
        if(a_noChildren(nghbrId) != 0) {
          addLayer = false;
        }
      } else {
        if(refineHalo[nghbrId] > 0) {
          addLayer = false;
        }
      }

      for(MInt d2 = 0; d2 < m_noDirs; d2++) {
        if(!addLayer && d == d2) {
          continue;
        }
        MInt nghbrId2 = a_neighborId(nghbrId, d2);
        if(nghbrId2 <= -1) {
          continue;
        }
        if(a_isHalo(nghbrId2) && layerId[nghbrId2] <= -1) {
          layerId[nghbrId2] = 2;
        }
        for(MInt d3 = 0; d3 < m_noDirs; d3++) {
          if(d3 == d || d3 == d2) {
            continue;
          }
          MInt nghbrId3 = a_neighborId(nghbrId2, d3);
          if(nghbrId3 <= -1) {
            continue;
          }
          if(a_isHalo(nghbrId3) && layerId[nghbrId3] <= -1) {
            layerId[nghbrId3] = 3;
          }
        }
      }
    }
  }

  MLongScratchSpace haloChildIds(2 * haloRefineCntTotal * m_maxNoChilds, AT_, "haloChildIds");
  cnt = 0;
  ASSERT(m_noHaloLayers == 2, ""); // If m_noHaloLayers != 2 everything will break!
  for(MInt i = 0; i < haloRefineCntTotal; i++) {
    MInt cellId = haloMap[i];

    MBool refine = false;

    // Only refine halo cells which has a refined internal cell in its vicinity.
    if(a_level(cellId) == level) {
      if(layerId[cellId] > 0) {
        refine = true;
      }
    }

    if(!refine) {
      fill(&haloChildIds[i * 2 * m_maxNoChilds], &haloChildIds[i * 2 * m_maxNoChilds] + (2 * m_maxNoChilds), -2);
      continue;
    }

    MInt parentId = cellId;
    for(MInt lvl = level; lvl > m_minLevel; lvl--) {
      parentId = a_parentId(parentId);
    }

    cartesianToCylindric(&a_coordinate(parentId, 0), coordsCyl);
    MInt side = m_azimuthalBbox.azimuthalSide(coordsCyl[1]);

    // If cell is supposed to be refined. Calculate hilbertIds of the rotated coordinate of possible childs
    MFloat hLength = F1B2 * cellLengthAtLevel(a_level(cellId));
    for(MInt child = 0; child < m_maxNoChilds; child++) {
      for(MInt d = 0; d < nDim; d++) {
        childCoords[d] = a_coordinate(cellId, d) + cornerStencil[child][d] * F1B2 * hLength;
      }
      rotateCartesianCoordinates(childCoords, side * dxCyl);

      MLong hilbertId =
          hilbertIndexGeneric(&childCoords[0], m_targetGridCenterOfGravity, m_targetGridLengthLevel0, level + 1);
      haloChildIds[i * 2 * m_maxNoChilds + 2 * child] = hilbertId;
      hilbertId = hilbertIndexGeneric(&childCoords[0], m_targetGridCenterOfGravity, m_targetGridLengthLevel0, level);
      haloChildIds[i * 2 * m_maxNoChilds + 2 * child + 1] = hilbertId;
    }
  }

  ScratchSpace<MPI_Request> mpi_send_req(noAzimuthalNeighborDomains(), AT_, "mpi_send_req");
  ScratchSpace<MPI_Request> mpi_recv_req(noAzimuthalNeighborDomains(), AT_, "mpi_recv_req");
  mpi_send_req.fill(MPI_REQUEST_NULL);
  mpi_recv_req.fill(MPI_REQUEST_NULL);

  MLongScratchSpace windowChildIds(2 * windowRefineCntTotal * m_maxNoChilds, AT_, "windowChildIds");
  cnt = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    if(windowRefineCnt[i] > 0) {
      MInt bufSize = 2 * m_maxNoChilds * windowRefineCnt[i];
      MPI_Irecv(&windowChildIds[cnt], bufSize, MPI_LONG, m_azimuthalNghbrDomains[i], 2, mpiComm(), &mpi_recv_req[i],
                AT_, "windowChildIds[cnt]");
      cnt += bufSize;
    }
  }

  cnt = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    if(haloRefineCnt[i] > 0) {
      MInt bufSize = 2 * m_maxNoChilds * haloRefineCnt[i];
      MPI_Isend(&haloChildIds[cnt], bufSize, MPI_LONG, m_azimuthalNghbrDomains[i], 2, mpiComm(), &mpi_send_req[i], AT_,
                "haloChildIds[cnt]");
      cnt += bufSize;
    }
  }

  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    if(windowRefineCnt[i] > 0) MPI_Wait(&mpi_recv_req[i], MPI_STATUSES_IGNORE, AT_);
    if(haloRefineCnt[i] > 0) MPI_Wait(&mpi_send_req[i], MPI_STATUSES_IGNORE, AT_);
  }

  // Now an internal cell is searched with the same hilbertId as the rotate halo child coordinate
  for(MInt i = 0; i < windowRefineCntTotal; i++) {
    MInt cellId = windowMap[i];

    MBool gridBndry = gridBndryCells[cellId];
    MBool noRefine = false;

    for(MInt child = 0; child < m_maxNoChilds; child++) {
      MLong haloHilbertId_L1 = windowChildIds[i * 2 * m_maxNoChilds + 2 * child];
      MLong haloHilbertId_L0 = windowChildIds[i * 2 * m_maxNoChilds + 2 * child + 1];
      windowChildIds[i * 2 * m_maxNoChilds + 2 * child] = -2;
      windowChildIds[i * 2 * m_maxNoChilds + 2 * child + 1] = -2;

      if(haloHilbertId_L1 < -1) {
        windowChildIds[i * m_maxNoChilds + child] = -2;
        continue;
      }

      MBool found = false;
      // First check in the childs of the parent window cell
      for(MInt child2 = 0; child2 < m_maxNoChilds; child2++) {
        MInt childId = a_childId(cellId, child2);
        if(childId < 0) continue;

        MLong hilbertId = hilbertIndexGeneric(&a_coordinate(childId, 0), m_targetGridCenterOfGravity,
                                              m_targetGridLengthLevel0, level + 1);

        if(hilbertId == haloHilbertId_L1) {
          windowChildIds[i * m_maxNoChilds + child] = a_globalId(childId);

          found = true;
          break;
        }
      }
      if(found) continue;

      // Now the neighboring cells of the internal cell are searched
      // They can either be on the child level or the current level
      // An azimuthal periodic connection between to cells, which are not on the same grid level
      // is possible!
      const MInt counter = getAdjacentGridCells(cellId, nghbrList, level);
      for(MInt n = 0; n < counter; n++) {
        MInt nghbrId = nghbrList[n];
        if(nghbrId < 0) continue;
        if(a_level(nghbrId) == level + 1) {
          MLong hilbertId = hilbertIndexGeneric(&a_coordinate(nghbrId, 0), m_targetGridCenterOfGravity,
                                                m_targetGridLengthLevel0, level + 1);
          if(hilbertId == haloHilbertId_L1) {
            windowChildIds[i * m_maxNoChilds + child] = a_globalId(nghbrId);

            found = true;
            break;
          }
        } else if(a_level(nghbrId) == level) {
          MLong hilbertId = hilbertIndexGeneric(&a_coordinate(nghbrId, 0), m_targetGridCenterOfGravity,
                                                m_targetGridLengthLevel0, level);

          if(hilbertId == haloHilbertId_L0) {
            windowChildIds[i * m_maxNoChilds + child] = a_globalId(nghbrId);
            found = true;
            break;
          }
        }
      }
      // If no adequate internal cellis found for on of the possible halo childs.
      // the halo cell will not be refined at all, alas it is a gridBndryCell, meaing
      // the rotated coordinate of the possible halo cell child could lie outside the
      // cartesian grid. In this case this cells becomes an unmapped halo cell and requires
      // special treatment
      if(!found) {
        if(gridBndry) {
          windowChildIds[i * m_maxNoChilds + child] = -1;
        } else {
          noRefine = true;
          break;
        }
      }
    }
    if(noRefine) {
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        windowChildIds[i * m_maxNoChilds + child] = -2;
      }
    }
  }

  cnt = 0;
  windowRefineCntTotal = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalWindowCells(i); j++) {
      if(refinedWindows[cnt] > 0) {
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          refineChildIds[cnt * m_maxNoChilds + child] = windowChildIds[windowRefineCntTotal * m_maxNoChilds + child];
        }
        windowRefineCntTotal++;
      } else {
        auto it = refineChildIds.begin() + (cnt * m_maxNoChilds);
        fill(it, it + m_maxNoChilds, -1);
      }
      cnt++;
    }
  }

  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    mpi_send_req[i] = MPI_REQUEST_NULL;
    mpi_recv_req[i] = MPI_REQUEST_NULL;
  }


  cnt = 0;
  haloChildIds.fill(-2);
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    if(haloRefineCnt[i] > 0) {
      MInt bufSize = m_maxNoChilds * haloRefineCnt[i];
      MPI_Irecv(&haloChildIds[cnt], bufSize, MPI_DOUBLE, m_azimuthalNghbrDomains[i], 2, mpiComm(), &mpi_recv_req[i],
                AT_, "haloChildIds[cnt]");
      cnt += bufSize;
    }
  }

  cnt = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    if(windowRefineCnt[i] > 0) {
      MInt bufSize = m_maxNoChilds * windowRefineCnt[i];
      MPI_Isend(&windowChildIds[cnt], bufSize, MPI_DOUBLE, m_azimuthalNghbrDomains[i], 2, mpiComm(), &mpi_send_req[i],
                AT_, "windowChildIds[cnt]");
      cnt += bufSize;
    }
  }

  MPI_Waitall(noAzimuthalNeighborDomains(), &mpi_recv_req[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noAzimuthalNeighborDomains(), &mpi_send_req[0], MPI_STATUSES_IGNORE, AT_);

  // GlobalIds of internal cells are send to halo rank.
  cnt = 0;
  haloRefineCntTotal = 0;
  for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < noAzimuthalHaloCells(i); j++) {
      if(refineCandidates[cnt] > 0) {
        for(MInt child = 0; child < m_maxNoChilds; child++) {
          recvChildIds[cnt * m_maxNoChilds + child] = haloChildIds[haloRefineCntTotal * m_maxNoChilds + child];
        }
        haloRefineCntTotal++;
      } else {
        auto it = recvChildIds.begin() + (cnt * m_maxNoChilds);
        fill(it, it + m_maxNoChilds, -2);
      }
      cnt++;
    }
  }
}


/** \brief Create a new azimuthal halo cell as neighbor of cellId in direction dir and find matching neighbor domain
   from hilbertIndex
   \author Thomas Hoesgen
  */
template <MInt nDim>
MInt CartesianGrid<nDim>::createAzimuthalHaloCell(const MInt cellId, const MInt dir, const MLong* hilbertOffsets,
                                                  unordered_multimap<MLong, MInt>& hilbertToLocal, const MFloat* bbox,
                                                  MInt* const noHalos, vector<vector<MInt>>& halos,
                                                  vector<vector<MLong>>& hilbertIds) {
  static constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};
  const MFloat cellLength = cellLengthAtCell(cellId);

  MFloat coords[3];
  MFloat coordsCyl[3];
  MFloat dcoords[3];
  MFloat dummyCoords[3];
  MFloat dxCyl = F0;
  for(MInt i = 0; i < nDim; i++) {
    coords[i] = a_coordinate(cellId, i);
  }
  coords[dir / 2] += ((dir % 2) == 0 ? -F1 : F1) * cellLength;

  if(!m_periodicCartesianDir[dir / 2] && (coords[dir / 2] < bbox[dir / 2] || coords[dir / 2] > bbox[nDim + dir / 2])) {
    return -1;
  }

  MBool isPeriodic = false;

  for(MInt i = 0; i < nDim; i++) {
    dummyCoords[i] = coords[i];
  }

  cartesianToCylindric(coords, coordsCyl);
  dxCyl = m_azimuthalAngle;

  if(coordsCyl[1] < m_azimuthalBbox.azimuthalBoundary(coords, -1) - MFloatEps) {
    isPeriodic = true;
    rotateCartesianCoordinates(dummyCoords, -dxCyl);
  } else if(coordsCyl[1] > m_azimuthalBbox.azimuthalBoundary(coords, 1) + MFloatEps) {
    isPeriodic = true;
    rotateCartesianCoordinates(dummyCoords, dxCyl);
  } else if(approx(dxCyl, 0.5 * PI, pow(10, -12))) {
    if(m_azimuthalBbox.azimuthalCenter(/*coords*/) > F0 && m_azimuthalBbox.azimuthalCenter(/*coords*/) < 0.5 * PI) {
      if(dir / 2 == 1 && coords[1] < bbox[1]) {
        isPeriodic = true;
        rotateCartesianCoordinates(dummyCoords, -dxCyl);
      }
      if(dir / 2 == 2 && coords[2] < bbox[2]) {
        isPeriodic = true;
        rotateCartesianCoordinates(dummyCoords, dxCyl);
      }
    } else {
      mTerm(1, AT_, "Do some implementing!");
    }
  }

  if(!isPeriodic) return -1;

  const MInt hilbertIndex = hilbertIndexGeneric(&dummyCoords[0]);

  MInt ndom = -1;
  while(hilbertIndex >= hilbertOffsets[ndom + 1]) {
    ndom++;
    if(ndom == noDomains() - 1) break;
  }
  ASSERT(ndom > -1, "");
  ASSERT(hilbertIndex >= hilbertOffsets[ndom] && hilbertIndex < hilbertOffsets[ndom + 1],
         to_string(hilbertIndex) + " " + to_string(hilbertOffsets[ndom]) + " " + to_string(hilbertOffsets[ndom + 1]));

  const MInt haloCellId = m_tree.size();
  m_tree.append();

  for(MInt i = 0; i < nDim; i++) {
    a_level(haloCellId) = a_level(cellId);
  }

  for(MInt i = 0; i < nDim; i++) {
    a_coordinate(haloCellId, i) = coords[i];
  }

  hilbertToLocal.insert(make_pair(hilbertIndex, haloCellId));

  for(MInt i = 0; i < m_noDirs; i++) {
    a_neighborId(haloCellId, i) = -1;
  }
  for(MInt i = 0; i < m_maxNoChilds; i++) {
    a_childId(haloCellId, i) = -1;
  }
  a_parentId(haloCellId) = -1;

  a_neighborId(cellId, dir) = haloCellId;
  a_neighborId(haloCellId, revDir[dir]) = cellId;
  for(MInt otherDir = 0; otherDir < m_noDirs; otherDir++) {
    if(otherDir / 2 == dir / 2) continue;
    if(a_hasNeighbor(cellId, otherDir) == 0) continue;
    MInt nghbrId = a_neighborId(cellId, otherDir);
    if(a_hasNeighbor(nghbrId, dir) > 0) {
      nghbrId = a_neighborId(nghbrId, dir);
      a_neighborId(nghbrId, revDir[otherDir]) = haloCellId;
      a_neighborId(haloCellId, otherDir) = nghbrId;
    }
  }

  for(MInt ndir = 0; ndir < m_noDirs; ndir++) {
    if(a_hasNeighbor(haloCellId, ndir) > 0) continue;
    for(MInt i = 0; i < nDim; i++) {
      dummyCoords[i] = a_coordinate(haloCellId, i);
    }
    dummyCoords[ndir / 2] += ((ndir % 2) == 0 ? -F1 : F1) * cellLength;
    for(MInt i = 0; i < nDim; i++) {
      dcoords[i] = dummyCoords[i];
    }
    MLong nghbrHilbertId = hilbertIndexGeneric(&dummyCoords[0]);
    MInt nghbrId = -1;
    auto range = hilbertToLocal.equal_range(nghbrHilbertId);
    for(auto it = range.first; it != range.second; ++it) {
      MFloat dist = F0;
      for(MInt i = 0; i < nDim; i++) {
        dist += POW2(a_coordinate(it->second, i) - dcoords[i]);
      }
      dist = sqrt(dist);
      if(dist < 0.001 * cellLength) {
        if(nghbrId > -1) cerr << "duplicate " << hilbertToLocal.count(nghbrHilbertId) << endl;
        nghbrId = it->second;
      }
    }
    if(nghbrId > -1) {
      a_neighborId(nghbrId, revDir[ndir]) = haloCellId;
      a_neighborId(haloCellId, ndir) = nghbrId;
#ifndef NDEBUG
      for(MInt i = 0; i < nDim; i++) {
        ASSERT(fabs(a_coordinate(nghbrId, i) - a_coordinate(haloCellId, i)) < cellLength * 1.0001, "");
      }
#endif
      continue;
    }

    for(MInt i = 0; i < nDim; i++) {
      cartesianToCylindric(dummyCoords, coordsCyl);
      dxCyl = m_azimuthalAngle;
      if(coordsCyl[1] < m_azimuthalBbox.azimuthalBoundary(coords, -1)) {
        rotateCartesianCoordinates(dummyCoords, -dxCyl);
      } else if(coordsCyl[1] > m_azimuthalBbox.azimuthalBoundary(coords, 1)) {
        rotateCartesianCoordinates(dummyCoords, dxCyl);
      }
    }
    nghbrHilbertId = hilbertIndexGeneric(&dummyCoords[0]);
    nghbrId = -1;
    range = hilbertToLocal.equal_range(nghbrHilbertId);
    for(auto it = range.first; it != range.second; ++it) {
      MFloat dist = F0;
      for(MInt i = 0; i < nDim; i++) {
        dist += POW2(a_coordinate(it->second, i) - dcoords[i]);
      }
      dist = sqrt(dist);
      if(dist < 0.001 * cellLength) {
        if(nghbrId > -1) cerr << "duplicate " << hilbertToLocal.count(nghbrHilbertId) << endl;
        nghbrId = it->second;
      }
    }
    if(nghbrId > -1) {
      a_neighborId(nghbrId, revDir[ndir]) = haloCellId;
      a_neighborId(haloCellId, ndir) = nghbrId;
#ifndef NDEBUG
      for(MInt i = 0; i < nDim; i++) {
        ASSERT(fabs(a_coordinate(nghbrId, i) - a_coordinate(haloCellId, i)) < cellLength * 1.0001, "");
      }
#endif
    }
  }

  a_resetProperties(haloCellId);
  a_hasProperty(haloCellId, Cell::IsHalo) = true;
  a_hasProperty(haloCellId, Cell::IsPeriodic) = isPeriodic;

  MInt idx = setAzimuthalNeighborDomainIndex(ndom, halos, hilbertIds);
  ASSERT(idx > -1, "");
  halos[idx].push_back(haloCellId);
  hilbertIds[idx].push_back(hilbertIndex);
  noHalos[ndom]++;

  return haloCellId;
}

/** \brief Tag azimuthal halo cells for refinement
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::tagAzimuthalUnmappedHaloCells(vector<vector<MLong>>& shiftWindows,
                                                        vector<MLong>& haloChildIds,
                                                        vector<MInt>& haloDomainIds,
                                                        MInt level) {
  TRACE();

  static constexpr MFloat cornerStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                                 {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
  MFloat childCoords[3] = {F0, F0, F0};
  MFloat coordsCyl[3] = {F0, F0, F0};
  MFloat dxCyl = m_azimuthalAngle;

  MInt cellCnt = 0;
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_level(cellId) != level + 1) continue;
    if(a_isHalo(cellId)) continue;
    cellCnt++;
  }

  unordered_multimap<MLong, MInt> hilbertToLocal;
  hilbertToLocal.reserve(cellCnt * 1.5);
  ScratchSpace<MLong> hilbertOffsets(noDomains() + 1, AT_, "hilbertOffsets");
  MLong minHilbertIndex = std::numeric_limits<MLong>::max();
  for(MInt cellId = 0; cellId < m_tree.size(); cellId++) {
    if(a_isHalo(cellId)) continue;
    MLong hilbertId = hilbertIndexGeneric(&a_coordinate(cellId, 0), m_targetGridCenterOfGravity,
                                          m_targetGridLengthLevel0, m_minLevel);
    minHilbertIndex = mMin(minHilbertIndex, hilbertId);
    if(a_level(cellId) != level + 1) continue;
    hilbertId =
        hilbertIndexGeneric(&a_coordinate(cellId, 0), m_targetGridCenterOfGravity, m_targetGridLengthLevel0, level + 1);
    hilbertToLocal.insert(make_pair(hilbertId, cellId));
  }
  MPI_Allgather(&minHilbertIndex, 1, MPI_LONG, &hilbertOffsets[0], 1, MPI_LONG, mpiComm(), AT_, "minHilbertIndex",
                "hilbertOffsets[0]");
  hilbertOffsets[0] = 0;
  hilbertOffsets[noDomains()] = std::numeric_limits<MLong>::max();

  if(m_noHaloLayers != 2) mTerm(1, AT_, "Only working with noHaloLayers=2");

  MIntScratchSpace noSndHilberts(noDomains(), AT_, "noSndHilberts");
  noSndHilberts.fill(0);

  shiftWindows.resize(noDomains());

  haloDomainIds.resize(m_maxNoChilds * noAzimuthalUnmappedHaloCells());
  fill(haloDomainIds.begin(), haloDomainIds.end(), -1);
  haloChildIds.resize(m_maxNoChilds * noAzimuthalUnmappedHaloCells());
  fill(haloChildIds.begin(), haloChildIds.end(), -1);

  MIntScratchSpace candidateMap(noAzimuthalUnmappedHaloCells(), AT_, "candidateMap");
  candidateMap.fill(-1);
  MLongScratchSpace candidateHilbertIds(m_maxNoChilds * noAzimuthalUnmappedHaloCells(), AT_, "candidateHilbertIds");
  candidateHilbertIds.fill(-1);
  MInt noRefineCandidates = 0;

  // Almost Always refine unmapped halo cells
  MIntScratchSpace nghbrList(200, AT_, "nghbrList");
  for(MInt i = 0; i < noAzimuthalUnmappedHaloCells(); i++) {
    MInt cellId = azimuthalUnmappedHaloCell(i);
    if(a_level(cellId) != level) continue;

    MBool refine = false;
    const MInt counter = getAdjacentGridCells(cellId, nghbrList, (level - 1), true);
    for(MInt n = 0; n < counter; n++) {
      MInt nghbrId = nghbrList[n];
      if(nghbrId < 0) continue;
      if(a_isHalo(nghbrId)) continue;
      if(a_noChildren(nghbrId) > 0) {
        refine = true;
      }
    }

    // First the hilbertIds of the child cells are calculated.
    // The hilbertIds are communicated to the mpi rank which has the
    // respective hilbertIdRange.
    // Then, it is checked if an internal cell can be found with the same hilbertId
    // to create an azimuthal window halo mapping.
    if(refine) {
      candidateMap[noRefineCandidates] = i;

      MInt parentId = cellId;
      for(MInt lvl = level; lvl > m_minLevel; lvl--) {
        parentId = a_parentId(parentId);
      }

      // It is more robust to determine the side with the parent on the minLevel
      cartesianToCylindric(&a_coordinate(parentId, 0), coordsCyl);
      MInt side = m_azimuthalBbox.azimuthalSide(coordsCyl[1]);

      MFloat hLength = F1B2 * cellLengthAtLevel(a_level(cellId));
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        for(MInt d = 0; d < nDim; d++) {
          childCoords[d] = a_coordinate(cellId, d) + cornerStencil[child][d] * F1B2 * hLength;
        }
        rotateCartesianCoordinates(childCoords, side * dxCyl);

        MLong hilbertId =
            hilbertIndexGeneric(childCoords, m_targetGridCenterOfGravity, m_targetGridLengthLevel0, m_minLevel);

        MInt dom = -1;
        for(MInt d = 0; d < noDomains(); d++) {
          if(hilbertOffsets[d + 1] > hilbertId) {
            dom = d;
            break;
          }
        }

        hilbertId = hilbertIndexGeneric(childCoords, m_targetGridCenterOfGravity, m_targetGridLengthLevel0, level + 1);

        haloDomainIds[i * m_maxNoChilds + child] = dom;
        candidateHilbertIds[noRefineCandidates * m_maxNoChilds + child] = hilbertId;
        noSndHilberts[dom]++;
      }
      noRefineCandidates++;
    } else {
      fill(&haloChildIds[i * m_maxNoChilds], &haloChildIds[i * m_maxNoChilds] + m_maxNoChilds, -2);
    }
  }


  // Now the hilbertIds are communicated
  MIntScratchSpace sndOffsets(noDomains() + 1, AT_, "sndOffsets");
  sndOffsets[0] = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    sndOffsets[d + 1] = sndOffsets[d] + noSndHilberts[d];
    noSndHilberts[d] = 0;
  }

  MLongScratchSpace sndHilbertIds(m_maxNoChilds * noRefineCandidates, AT_, "sndHilbertIds");
  for(MInt i = 0; i < noRefineCandidates; i++) {
    MInt candidateId = candidateMap[i];
    for(MInt child = 0; child < m_maxNoChilds; child++) {
      MInt dom = haloDomainIds[candidateId * m_maxNoChilds + child];
      sndHilbertIds[sndOffsets[dom] + noSndHilberts[dom]] = candidateHilbertIds[i * m_maxNoChilds + child];
      noSndHilberts[dom]++;
    }
  }


  MIntScratchSpace noRcvHilberts(noDomains(), AT_, "noRcvHilberts");
  MPI_Alltoall(&noSndHilberts[0], 1, MPI_INT, &noRcvHilberts[0], 1, MPI_INT, mpiComm(), AT_, "noSndHilberts[0]",
               "noRcvHilberts[0]");


  MInt rcvHilbertsTotal = 0;
  MIntScratchSpace rcvOffsets(noDomains() + 1, AT_, "rcvOffsets");
  rcvOffsets[0] = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    rcvOffsets[d + 1] = rcvOffsets[d] + noRcvHilberts[d];
    rcvHilbertsTotal += noRcvHilberts[d];
  }
  MLongScratchSpace rcvHilbertIds(rcvHilbertsTotal, AT_, "rcvHilbertIds");


  ScratchSpace<MPI_Request> sendReq(noDomains(), AT_, "sendReq");
  sendReq.fill(MPI_REQUEST_NULL);

  for(MInt i = 0; i < noDomains(); i++) {
    if(noSndHilberts[i] == 0) continue;
    MPI_Issend(&sndHilbertIds[sndOffsets[i]], noSndHilberts[i], type_traits<MLong>::mpiType(), i, 24, mpiComm(),
               &sendReq[i], AT_, "sndHilbertIds[sndOffsets[i]");
  }

  for(MInt i = 0; i < noDomains(); i++) {
    if(noRcvHilberts[i] == 0) continue;
    MPI_Recv(&rcvHilbertIds[rcvOffsets[i]], noRcvHilberts[i], type_traits<MLong>::mpiType(), i, 24, mpiComm(),
             MPI_STATUS_IGNORE, AT_, "rcvHilbertIds[rcvOffsets[i]]");
  }

  for(MInt i = 0; i < noDomains(); i++) {
    if(noSndHilberts[i] == 0) continue;
    MPI_Wait(&sendReq[i], MPI_STATUSES_IGNORE, AT_);
  }


  // Check if internal cell can be found which has same hilbertId
  for(MInt i = 0; i < noDomains(); i++) {
    for(MInt j = 0; j < noRcvHilberts[i]; j++) {
      MLong hilbertId = rcvHilbertIds[rcvOffsets[i] + j];

      MInt windowId = -1;
      auto range = hilbertToLocal.equal_range(hilbertId);
      for(auto it = range.first; it != range.second; ++it) {
        if(!a_hasProperty(it->second, Cell::IsPeriodic) && !a_hasProperty(it->second, Cell::IsHalo)) {
          windowId = it->second;
        }
      }
      ASSERT(windowId > -2, to_string(windowId) + " " + to_string(m_noInternalCells));

      if(windowId > -1) {
        rcvHilbertIds[rcvOffsets[i] + j] = a_globalId(windowId);
        shiftWindows[i].push_back(rcvHilbertIds[rcvOffsets[i] + j]);
      } else {
        rcvHilbertIds[rcvOffsets[i] + j] = -1;
      }
    }
  }

  // Communicated globalIds of found internal cells back to halo rank
  sendReq.fill(MPI_REQUEST_NULL);

  for(MInt i = 0; i < noDomains(); i++) {
    if(noRcvHilberts[i] == 0) continue;
    MPI_Issend(&rcvHilbertIds[rcvOffsets[i]], noRcvHilberts[i], type_traits<MLong>::mpiType(), i, 12, mpiComm(),
               &sendReq[i], AT_, "rcvHilbertIds[rcvOffsets[i]]");
  }

  for(MInt i = 0; i < noDomains(); i++) {
    if(noSndHilberts[i] == 0) continue;
    MPI_Recv(&sndHilbertIds[sndOffsets[i]], noSndHilberts[i], type_traits<MLong>::mpiType(), i, 12, mpiComm(),
             MPI_STATUS_IGNORE, AT_, "sndHilbertIds[sndOffsets[i]");
  }

  for(MInt i = 0; i < noDomains(); i++) {
    if(noRcvHilberts[i] == 0) continue;
    MPI_Wait(&sendReq[i], MPI_STATUSES_IGNORE, AT_);
  }

  for(MInt d = 0; d < noDomains(); d++) {
    noSndHilberts[d] = 0;
  }

  for(MInt i = 0; i < noRefineCandidates; i++) {
    MInt candidateId = candidateMap[i];
    for(MInt child = 0; child < m_maxNoChilds; child++) {
      MInt dom = haloDomainIds[candidateId * m_maxNoChilds + child];
      haloChildIds[candidateId * m_maxNoChilds + child] = sndHilbertIds[sndOffsets[dom] + noSndHilberts[dom]++];
    }
  }
}

/** \brief Checks consistency of aimuthal window/halo cells
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::checkAzimuthalWindowHaloConsistency() {
  TRACE();
  if(globalNoDomains() == 1) {
    return;
  }

  m_log << "checkWindowHaloConsistency azimuthal... ";

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #1: number of azimuthal window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving number of window cells from each neighbor domain
  ScratchSpace<MPI_Request> recvRequests(std::max(1, noAzimuthalNeighborDomains()), AT_, "recvRequests");
  MIntScratchSpace noWindowCellsRecv(std::max(1, noAzimuthalNeighborDomains()), AT_, "noWindowCellsRecv");
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    noWindowCellsRecv[d] = -1;
    MPI_Irecv(&noWindowCellsRecv[d], 1, type_traits<MInt>::mpiType(), azimuthalNeighborDomain(d),
              azimuthalNeighborDomain(d), mpiComm(), &recvRequests[d], AT_, "noWindowCellsRecv[d]");
  }

  // Start sending number of window cells to each neighbor domain
  ScratchSpace<MPI_Request> sendRequests(std::max(1, noAzimuthalNeighborDomains()), AT_, "sendRequests");
  MIntScratchSpace noWindowCellsSend(std::max(1, noAzimuthalNeighborDomains()), AT_, "noWindowCellsSend");
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    noWindowCellsSend[d] = noAzimuthalWindowCells(d);
    MPI_Isend(&noWindowCellsSend[d], 1, type_traits<MInt>::mpiType(), azimuthalNeighborDomain(d), domainId(), mpiComm(),
              &sendRequests[d], AT_, "noWindowCellsSend[d]");
  }

  // Finish MPI communication
  MPI_Waitall(noAzimuthalNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noAzimuthalNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received number of window cells matches the local number of halo cells
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    if(noWindowCellsRecv[d] != noAzimuthalHaloCells(d)) {
      TERMM(1, "Cartesian Grid : Number of azimuthal window cells from domain " + to_string(azimuthalNeighborDomain(d))
                   + " does not match local number of azimuthal halo cells; window: " + to_string(noWindowCellsRecv[d])
                   + " ,halo: " + to_string(noAzimuthalHaloCells(d)));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #2: grid global ids of window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving window cell global ids from each neighbor domain
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsRecv = accumulate(noWindowCellsRecv.begin(), noWindowCellsRecv.end(), 0);

  MIntScratchSpace windowCellsRecv(max(1, totalNoWindowCellsRecv), AT_, "windowCellsRecv");
  fill(windowCellsRecv.begin(), windowCellsRecv.end(), -1);
  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    if(noAzimuthalHaloCells(d) > 0) {
      MPI_Irecv(&windowCellsRecv[offset], noAzimuthalHaloCells(d), type_traits<MInt>::mpiType(),
                azimuthalNeighborDomain(d), azimuthalNeighborDomain(d), mpiComm(), &recvRequests[d], AT_,
                "windowCellsRecv[offset]");
    }
    offset += noAzimuthalHaloCells(d);
  }

  // Start sending window cell global ids to each neighbor domain
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsSend = accumulate(noWindowCellsSend.begin(), noWindowCellsSend.end(), 0);
  MIntScratchSpace windowCellsSend(max(1, totalNoWindowCellsSend), AT_, "windowCellsSend");

  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    for(MInt c = 0; c < noWindowCellsSend[d]; c++) {
      TERMM_IF_NOT_COND(a_hasProperty(azimuthalWindowCell(d, c), Cell::IsWindow), "not a window cell");
      windowCellsSend[offset + c] = a_globalId(azimuthalWindowCell(d, c));
    }
    if(noAzimuthalWindowCells(d) > 0) {
      MPI_Isend(&windowCellsSend[offset], noAzimuthalWindowCells(d), type_traits<MInt>::mpiType(),
                azimuthalNeighborDomain(d), domainId(), mpiComm(), &sendRequests[d], AT_, "windowCellsSend[offset]");
    }
    offset += noAzimuthalWindowCells(d);
  }

  // Finish MPI communication
  MPI_Waitall(noAzimuthalNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noAzimuthalNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received window cell global ids match the local halo cell global ids
  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    for(MInt c = 0; c < noAzimuthalHaloCells(d); c++) {
      const MInt cellId = azimuthalHaloCell(d, c);
      TERMM_IF_NOT_COND(a_isHalo(cellId), "not a halo cell");
      const MInt globalId = a_globalId(cellId);
      // If halo cell has periodic flag, its global id should be -1
      if(windowCellsRecv[offset + c] != globalId && !(a_hasProperty(cellId, Cell::IsPeriodic) && globalId == -1)) {
        TERMM(1, "Global id of window cell " + to_string(c) + " from domain " + to_string(azimuthalNeighborDomain(d))
                     + " does not match local halo cell gobal id (" + to_string(windowCellsRecv[offset + c]) + " vs. "
                     + to_string(a_globalId(azimuthalHaloCell(d, c))) + ")");
      }
    }
    offset += noAzimuthalHaloCells(d);
  }

  m_log << "done" << std::endl;
}

/** \brief Corrects solver bits on azimuthal halo cells
           Ensure that azimuthal halo cells have all children
           or are not refined at all. Further ensure that
           parent cells have solverBit = true
    \author Thomas Hoesgen
    \date November 2020
*/
template <MInt nDim>
void CartesianGrid<nDim>::correctAzimuthalSolverBits() {
  TRACE();

  MInt maxLvl = mMin(m_maxRfnmntLvl, m_maxLevel);
  for(MInt level = m_maxUniformRefinementLevel; level < maxLvl; level++) {
    for(MInt i = 0; i < noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)m_azimuthalHaloCells[i].size(); j++) {
        MInt cellId = m_azimuthalHaloCells[i][j];
        if(a_level(cellId) != level) continue;
        for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
          if(m_tree.solver(cellId, solver)) {
            MBool noRefine = false;
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId == -1) continue;
              if(!m_tree.solver(childId, solver)) {
                noRefine = true;
                break;
              }
            }
            if(noRefine) {
              for(MInt child = 0; child < m_maxNoChilds; child++) {
                MInt childId = a_childId(cellId, child);
                if(childId == -1) continue;
                m_tree.solver(childId, solver) = false;
              }
            }
          } else {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId == -1) continue;
              m_tree.solver(childId, solver) = false;
            }
          }
        }
      }
    }


    for(MInt i = 0; i < (signed)m_azimuthalUnmappedHaloCells.size(); i++) {
      MInt cellId = m_azimuthalUnmappedHaloCells[i];
      if(a_level(cellId) != level) continue;
      for(MInt solver = 0; solver < treeb().noSolvers(); solver++) {
        if(m_tree.solver(cellId, solver)) {
          MBool noRefine = false;
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt childId = a_childId(cellId, child);
            if(childId == -1) continue;
            if(!m_tree.solver(childId, solver)) {
              noRefine = true;
              break;
            }
          }
          if(noRefine) {
            for(MInt child = 0; child < m_maxNoChilds; child++) {
              MInt childId = a_childId(cellId, child);
              if(childId == -1) continue;
              m_tree.solver(childId, solver) = false;
            }
          }
        } else {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            MInt childId = a_childId(cellId, child);
            if(childId == -1) continue;
            m_tree.solver(childId, solver) = false;
          }
        }
      }
    }
  }
}


/** \brief Checks variable m_windowLayer_
    \date 21th century
  */
template <MInt nDim>
void CartesianGrid<nDim>::checkWindowLayer(const MString text) {
  TRACE();

  TERMM_IF_COND((signed)m_windowLayer_.size() < noNeighborDomains(), text);
  // Ensure that all cells in m_windowLayer_ have proper layer
  for(MInt d = 0; d < noNeighborDomains(); ++d) {
    for(const auto item : m_windowLayer_[d]) {
      const MInt cellId = item.first;
      TERMM_IF_NOT_COND(cellId > -1 && cellId < m_tree.size(), "");
      TERMM_IF_COND(item.second.all(), "This entry is unnecessary!");

      for(MInt solver = 0; solver < treeb().noSolvers(); ++solver) {
        TERMM_IF_COND(item.second.get(solver) <= m_noSolverHaloLayers[solver] /*M32X4bit<>::MAX*/
                          && !m_tree.solver(cellId, solver),
                      text);
      }
    }
  }

  // Check that all cells in m_windowCells are also in m_windowLayer_
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    std::set<MInt> allWs(m_windowCells[i].begin(), m_windowCells[i].end());
    for(MInt j = 0; j < (signed)m_windowCells[i].size(); j++) {
      const MInt cellId = m_windowCells[i][j];
      if(m_windowLayer_[i].find(cellId) == m_windowLayer_[i].end()) {
        stringstream ss;
        ss << " domainId=" << domainId() << ": i=" << i << " j=" << j << " cellId=" << cellId
           << " globalId=" << a_globalId(cellId) << " nghbrDom=" << m_nghbrDomains[i] << endl;
        TERMM(1, text + ss.str());
      }
    }

    // Check that all cells in m_windowLayers_ are also in m_windowCells
    std::set<std::pair<MInt, MInt>> del;
    for(auto m : m_windowLayer_[i]) {
      if(allWs.find(m.first) == allWs.end()) {
        // Note: If following assert fails, one could try to uncomment the next line
        // del.insert(make_pair(i,m.first)); continue;
        std::stringstream ss;
        ss << text << " d=" << domainId() << " nghbrDom=" << m_nghbrDomains[i] << "(" << i
           << ") minLevel=" << m_minLevel << " cell_level=" << a_level(m.first) << " isHalo=" << a_isHalo(m.first)
           << " cellId=" << m.first << " cellLayer=" << m.second.get(0) << "|" << m.second.get(1)
           << " parentId=" << a_parentId(m.first) << " parentIsHalo=" << a_isHalo(std::max(0, (int)a_parentId(m.first)))
           << " " << isSolverWindowCell(i, m.first, 0);
        TERMM(1, ss.str());
      }
    }
    for(auto d : del) {
      m_windowLayer_[d.first].erase(m_windowLayer_[d.first].find(d.second));
    }

    if(m_noPeriodicCartesianDirs == 0) TERMM_IF_NOT_COND(m_windowCells[i].size() == m_windowLayer_[i].size(), text);
  }
}
