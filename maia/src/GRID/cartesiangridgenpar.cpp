// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangridgenpar.h"

#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "UTIL/functions.h"
#include "UTIL/hilbert.h"

using namespace std;
using namespace maia;

#include <vector>
#include "GEOM/geometryanalytic.h"
#include "GEOM/geometryroot.h"
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"

// prevent 64bit to 32bit integer narrowing conversion
#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wconversion"
#endif

/** \brief generates a Cartesian grid in parallel
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * The complete algorithm consitst of the following steps:
 *
 *   1. Read properties
 *   2. Init members
 *   3. Init geometry
 *   4. Create initial grid
 *     4.1 Creates an initial cube around the geometry.
 *     4.2 Refines the initial cube to the initialGeometricalRfnLvl.
 *         In each iteration the following is performed:
 *       4.2.1 update the level offsets:
 *           this gets rid of the delettion of lower levels by using the following structure.
 *           Based on the m_initialGeometricalRfnLevel the first initial cube is located at the
 *           beginning or at the end of the collector:
 *             - m_initialGeometricalRfnLevel is even: locate at the beginning
 *             - m_initialGeometricalRfnLevel is odd: locate at the end
 *           The next following level is always written to the other end of the collector, i.e.,
 *           if the first cube is located at the end of the collector, the next higher level starts
 *           at the beginning of the collector. This method iterates, cells are written to the front
 *           back of the collector continuously. This way it is not necessary to delete the lower levels
 *           since they are just overwritten by one of the upcoming levels. This way, at the end of
 *           the refinement the initial level is located at the beginning of the collector. Before refining
 *           and writing to a level at the end of the collector the offsets are determined by estimating the
 *           possible number of cells on the upcoming level.
 *       4.2.2 check if there is enough memory
 *       4.x.x each solver marks its cells for refinement
 *       4.2.3 refine the grid one further level
 *       4.2.4 update the neighborhood and solver information on the newly created cells
 *       4.2.5 cells outside the geometry are deleted
 *       4.2.6 cells outside the individual solvers are unmarked
 *     4.3 Remove the link to the parents on the initial level
 *     4.4 Update the number of generated cells and print information
 *     4.5 Reorder cells after Hilbert id
 *     4.6 If the code runs parallel, parallelize the grid.
 *  5. Create start grid (if required), can be found in createComputationalMultisolverGrid
 *     5.1. Run over all levels that are still to be refined and do the following:
 *       5.1.1 Update the m_levelOffsets of the next level. This generates a new
 *             offset range for the new level to create under the assumption
 *             that all cells on the current level are refined.
 *       5.1.2 Update the halo cell offsets. This call the updateHaloCellOffsets(...)
 *             function.
 *       5.1.3 Check memory availability.
 *       5.1.4 Refine the normal cells and the halo cells. This calls refineGrid(...)
 *             similar to the serial case. To refine the halo cells the new halo cell
 *             offsets are provided.
 *       5.1.5 Update the number of cells on this domain, since it has changed.
 *       5.1.6 Find the neighbors for all newly generated normal and halo cells.
 *             This calls findChildLevelNeighbors(...) and provides the offsets for
 *             the normal cells and the halo cells of the current level.
 *       5.1.7 Delete all normal outside cells. This has still to be done for the
 *             halo cells.
 *     5.2 Again update the number of cells.
 *  6. Create computational multisolver grid (if required)
 *     6.1 This calls either the method for patches
 *         or the method for local boundary refinement
 *       6.1a.2 Run over all levels which still need to be refined. Depending on the properties
 *              for this level (is it a box patch or a sphere patch?) the grid is refined
 *         6.1a.2.1 Check what type of patch should be used and mark the according cells. This calls
 *                  either refineLocalBox(...), which marks the cells according to the defined
 *                  coordinates in a box-style or refineLocalRadius(...), which marks the cells
 *                  according to the coordinates of the center of a sphere and a radius. If we
 *                  run in parallel mode, also do this for the halos. The functions refineLocalBox(...)
 *                  and refineLocalRadius(...) already update the offsets for the normal cells
 *                  and the halos.
 *         6.1a.2.2 Check memory availability.
 *         6.1a.2.3 Refine the cells that have previously been marked. This calls the function
 *                  refineGridPatch(...) which only refines cells that have previously been
 *                  marked. If we run in parallel mode, do this also for the halo cells.
 *         6.1b.2.4 Update the number of cells.
 *         6.1b.2.5 Find the neighbors for the new cells. Calls findChildLevelNeighbors(...).
 *                  If we run in parallel mode, this is also done for the halo cells.
 *         6.1b.2.6 Delete all outside cells, do this in serial or in parallel. This calls
 *                  either deleteOutsideCellsSerial(...) or deleteOutsideCellsParallel(...)
 *       6.1a.1 Initialize the distances for all upcoming levels. The distances are based
 *              on the value provided by the property localMinBoundaryThreshold and by
 *              the property smoothDistance. The first defines the distance in cell units
 *              on the current level, which should be refined. The latter property defines
 *              how many cells should be inbetween the levels.
 *         6.3b.2 Run over all levels that are still to be refined and do the following:
 *         6.3b.2.1 Do the boundary propagation. This calls boundaryPropagation(...)
 *                  with the current level and the final distance to use.
 *         6.3b.2.2 Dry run for the normal cells: This finds out how many cells have been
 *                  marked in boundaryPropagation(...) and recalulcates the new offsets
 *                  for the new level. Mark the according cells which are inside the given
 *                  distance: calls markBndDistance(...)
 *         6.2b.2.3 Check memory availability.
 *         6.2b.2.4 Update the distances for the new level. A distance on a given level
 *                  is twice the distance on a higher level.
 *         6.2b.2.5 Refine the cells that have previously been marked. This calls the function
 *                  refineGridPatch(...) which only refines cells that have previously been
 *                  marked. If we run in parallel mode, do this also for the halo cells.
 *         6.2b.2.6 Update the number of cells.
 *         6.2b.2.7 Find the neighbors for the new cells. Calls findChildLevelNeighbors(...).
 *                  If we run in parallel mode, this is also done for the halo cells.
 *         6.2b.2.8 Delete all outside cells, do this in serial or in parallel. This calls
 *                  either deleteOutsideCellsSerial(...) or deleteOutsideCellsParallel(...)
 *     6.3 Check the refinement validity for LB.
 *     6.4 We are finished, run the final setup. This calls updateInterRankNeighbors() to
 *         exchange the neighborhood.
 *  7. Write grid to file
 *  8. Cleaning up
 *
 * \param[in] dimensions the dimensionality of the problem
 *
 **/
template <MInt nDim>
GridgenPar<nDim>::GridgenPar(const MPI_Comm comm, const MInt noSolvers)
  : m_mpiComm(comm), outStream(nullptr), m_noSolvers(noSolvers) {
  TRACE();

  MPI_Comm_rank(mpiComm(), &m_domainId);
  MPI_Comm_size(mpiComm(), &m_noDomains);

  if(globalDomainId() == 0) {
    outStream.rdbuf(std::cout.rdbuf());
  }

  m_log << "Parallel grid generator started on process " << globalDomainId() << endl;
  std::cout << "Parallel grid generator started on process " << globalDomainId() << endl;

  // 0. Initialize timers
  initTimers();
  RECORD_TIMER_START(m_t_comp_GG);

  // 1. Read properties
  readProperties();

  // 2. Init members
  initMembers();

  // 3. Init geometry
  initGeometry();

  writeMemoryStatistics(mpiComm(), noDomains(), domainId(), AT_, "Gridgen after init");

  // 3 - 4
  gridAlignCutOff();

  // 4. Create initial grid
  createInitialGrid();

  if(noDomains() > 1) parallelizeGrid();

  // 5. and 6.
  createComputationalMultisolverGrid();

  // 6.4 we are finished, run the final setup
  finalizeGrid();

  // for debugging
  if(m_writeGridInformation) writeGridInformationPar();

  // 7. Write grid to file
  saveGrid();

  writeMemoryStatistics(mpiComm(), noDomains(), domainId(), AT_, "Gridgen finalize");

  // 8. Cleaning up
  m_log << "  (8) Cleaning up" << endl;
  outStream << "  (8) Cleaning up" << endl;

  RECORD_TIMER_STOP(m_t_comp_GG);
}

template <MInt nDim>
void GridgenPar<nDim>::initTimers() {
  TRACE();

  NEW_TIMER_GROUP(tgrp_GG, "Parallel grid generation");
  NEW_TIMER(t_comp_GG, "complete grid generation", tgrp_GG);
  NEW_SUB_TIMER(t_readProperties, "(1) read properties", t_comp_GG);
  NEW_SUB_TIMER(t_initMembers, "(2) init members", t_comp_GG);
  NEW_SUB_TIMER(t_initGeometry, "(3) init geometry", t_comp_GG);
  NEW_SUB_TIMER(t_createInitialGrid, "(4) create initial grid", t_comp_GG);
  NEW_SUB_TIMER(t_parallelizeGrid, "(5) distribute grid among processes", t_comp_GG);
  NEW_SUB_TIMER(t_createStartGrid, "(6) create start grid", t_comp_GG);
  NEW_SUB_TIMER(t_createComputationalGrid, "(7) create computational grid", t_comp_GG);
  NEW_SUB_TIMER(t_finalizeGrid, "(8) prepare grid for IO", t_comp_GG);
  NEW_SUB_TIMER(t_saveGrid, "(9) parallel grid I/O", t_comp_GG);

  m_t_comp_GG = t_comp_GG;
  m_t_readProperties = t_readProperties;
  m_t_initMembers = t_initMembers;
  m_t_initGeometry = t_initGeometry;
  m_t_createInitialGrid = t_createInitialGrid;
  m_t_parallelizeGrid = t_parallelizeGrid;
  m_t_createStartGrid = t_createStartGrid;
  m_t_createComputationalGrid = t_createComputationalGrid;
  m_t_finalizeGrid = t_finalizeGrid;
  m_t_saveGrid = t_saveGrid;
}

/**
 * \fn void GridgenPar<nDim>::readProperties()
 * \brief reads necessary properties from the property file
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::readProperties() {
  TRACE();

  RECORD_TIMER_START(m_t_readProperties);

  m_log << "  (1) reading properties" << endl;
  outStream << "  (1) reading properties" << endl;

  mAlloc(m_solverRefinement, m_noSolvers, "m_solverRefinement", AT_);
  /*! \page propertyPage1
    \section initialGeometricalRfnLvl
    <code>MInt GridGenPar::m_initialRefinementLevelSerial </code>\n
    default = <code>n/a</code>\n \n
    Specifies the coarsest refinement level that is retained and where ordering by the Hilbert
    curve takes place (deprecated).\n
    Possible values are:
    <ul>
      <li>Non-negative int value</li>
    </ul>
    Keywords: <i>GRID GENERATOR, PARALLEL, REFINEMENT</i>
   */

  /*! \page propertyPage1
    \section minLevel
    <code>MInt GridGenPar::m_initialRefinementLevelSerial </code>\n
    default = <code>n/a</code>\n \n
    Specifies the coarsest refinement level that is retained and where ordering by the Hilbert
    curve takes place.\n
    Possible values are:
    <ul>
      <li>Non-negative int value</li>
    </ul>
    Keywords: <i>GRID GENERATOR, PARALLEL, REFINEMENT</i>
   */
  if(Context::propertyExists("initialGeometricalRfnLvl", 0)) {
    // mTerm(1, AT_, "Error: Property initialGeometricalRfnLvl is deprecated. Please rename to minLevel.");
    cerr << "Warning: Property initialGeometricalRfnLvl is deprecated. Please rename to minLevel." << endl;
    m_initialRefinementLevelSerial = Context::getBasicProperty<MInt>("initialGeometricalRfnLvl", AT_);
  } else {
    m_initialRefinementLevelSerial = Context::getBasicProperty<MInt>("minLevel", AT_);
  }

  m_minLevel = m_initialRefinementLevelSerial;
  // currently, these equal. However, it should be possible to have m_initialRefinementLevelSerial < m_minLevel, because
  // for big grids serial grid generation will run out of memory!

  // to minimize
  m_maxUniformRefinementLevel = numeric_limits<MInt>::max();
  // to maximize
  m_maxRfnmntLvl = numeric_limits<MInt>::min();
  // "maximize"
  m_cutOff = false;
  m_weightPatchCells = 0;
  m_weightBndCells = 0;
  m_weightSolverUniformLevel = 0;
  m_weightSolverUniformLevel =
      Context::getBasicProperty<MBool>("weightSolverUniformLevel", AT_, &m_weightSolverUniformLevel);

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    readSolverProperties(solver);
  }

  /*! \page propertyPage1
    \section weightMethod
    <code>MInt GridGenPar::m_weightMethod </code>\n
    default = <code>0</code>\n \n
    Specifies the weight method that should be used to partition the grid.\n
    Possible values are:
    <ul>
      <li>0: all cells have the same weight (= 1.0)</li>
      <li>N: different weight methods are supported (see `setCellWeights()` method)</li>
      <li>3: set specific weight for each level (using 'setWeightValue' from max to min level)</li>
    </ul>
    Keywords: <i>GRID GENERATOR, PARALLEL, CELL WEIGHTS</i>
   */
  m_weightMethod = 0;
  m_weightMethod = Context::getBasicProperty<MInt>("weightMethod", AT_, &m_weightMethod);

  m_writeGridInformation = 0;
  m_writeGridInformation = Context::getBasicProperty<MBool>("writeGridInformation", AT_, &m_writeGridInformation);

  /*! \page propertyPage1
    \section reductionFactor
    <code>MFloat GridGenPar::m_reductionFactor </code>\n
    default = <code>1.0</code>\n \n
    Factor by which the length of the level-0 cell is multiplied. This can be
    used to obtain cells of a different cell length, independent of the largest
    extent of the bounding box.\n
    Possible values are:
    <ul>
      <li>1.0 <= reductionFactor < 2.0</li>
    </ul>
    Keywords: <i>GRID GENERATOR, PARALLEL, GEOMETRY, LENGTH</i>
   */
  m_reductionFactor = 1.0;
  m_reductionFactor = Context::getBasicProperty<MFloat>("reductionFactor", AT_, &m_reductionFactor);

  /*! \page propertyPage1
    \section gridOutputFileName
    <code>MString GridGenPar::m_gridOutputFileName </code>\n
    default = <code>n/a</code>\n \n
    Name of the grid file that is generated.\n
    Possible values are:
    <ul>
      <li>any valid file name</li>
    </ul>
    Keywords: <i>GRID GENERATOR, PARALLEL, FILENAME</i>
   */
  m_gridOutputFileName = Context::getBasicProperty<MString>("gridOutputFileName", AT_);

  /*! \page propertyPage1
    \section writeCoordinatesToGridFile
    <code>MString GridGenPar::m_writeCoordinatesToGridFile</code>\n
    default = <code>false</code>\n \n
    Write all cell coordinates to the grid file, e.g. for external processing without the grid reader.\n
    Keywords: <i>GRID, GENERATOR, OUTPUT, COORDINATES</i>
   */
  m_writeCoordinatesToGridFile = false;
  m_writeCoordinatesToGridFile =
      Context::getBasicProperty<MBool>("writeCoordinatesToGridFile", AT_, &m_writeCoordinatesToGridFile);

  stringstream errorMsg;
  if(m_reductionFactor < 1.0 || m_reductionFactor >= 2.0) {
    errorMsg << "ERROR: reductionFactor must be between 1.0 and 2.0, but is " << m_reductionFactor << endl;
    m_log << errorMsg.str();
    mTerm(1, AT_, errorMsg.str());
  }

  m_maxNoCells = Context::getBasicProperty<MInt>("maxNoCells", AT_, &m_maxNoCells);

  m_partitionCellOffspringThreshold = 50000;
  /*! \page propertyPage1
    \section partitionCellOffspringThreshold
    <code>MInt GridGenPar::m_partitionCellOffspringThreshold </code>\n
    default = <code>50000</code>\n \n
    controls the filtering of cells in massive parallel file writing\n
    possible values are:
    <ul>
    <li>Non-negative int value of the order less than number of cells on each domains</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, OUTPUT, FILTER, MAX, SIZE</i>
   */
  if(Context::propertyExists("minCellMaxSize", 0)) {
    mTerm(1, AT_, "Error: Property minCellMaxSize is deprecated. Please rename to partitionCellOffspringThreshold.");
  } else if(Context::propertyExists("partitionCellMaxNoOffspring", 0)) {
    cerr << "Warning: Property partitionCellMaxNoOffspring is deprecated. Please rename to "
            "partitionCellOffspringThreshold."
         << endl;
    MInt tmpPartitionCellMaxSize = (MInt)m_partitionCellOffspringThreshold;
    m_partitionCellOffspringThreshold =
        (MLong)(Context::getBasicProperty<MInt>("partitionCellMaxNoOffspring", AT_, &tmpPartitionCellMaxSize));
  } else {
    MInt tmpPartitionCellMaxSize = (MInt)m_partitionCellOffspringThreshold;
    m_partitionCellOffspringThreshold =
        (MLong)(Context::getBasicProperty<MInt>("partitionCellOffspringThreshold", AT_, &tmpPartitionCellMaxSize));
  }

  m_partitionCellWorkloadThreshold = 50000.;
  /*! \page propertyPage1
    \section partitionCellWorkloadThreshold
    <code>MFloat GridGenPar::m_partitionCellWorkloadThreshold </code>\n
    default = <code>50000</code>\n \n
    controls the filtering of cells in massive parallel file writing\n
    possible values are:
    <ul>
    <li>Non-negative int value of the order less than number of cells on each domains</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, OUTPUT, FILTER, MAX, SIZE</i>
   */
  m_partitionCellWorkloadThreshold =
      Context::getBasicProperty<MFloat>("partitionCellWorkloadThreshold", AT_, &m_partitionCellWorkloadThreshold);

  /*! \page propertyPage1
    \section targetGridFileName
    <code>MString GridGenPar::m_targetGridFileName </code>\n
    default = <code>""</code>\n \n
    If non-empty, use target grid center of gravity and lengthLevel0 when
    ordering cells by Hilbert id.\n
    possible values are:
    <ul>
    <li>Full path to target grid file.</li>
    <li>"" (empty string)</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, COUPLING</i>
   */
  m_targetGridFileName = "";
  m_targetGridFileName = Context::getBasicProperty<MString>("targetGridFileName", AT_, &m_targetGridFileName);

  if(Context::propertyExists("multiSolverBoundingBox")) {
    m_hasMultiSolverBoundingBox = true;
    /*! \page propertyPage1
       \section multiSolverBoundingBox
       <code>std::vector<MFloat> CartesianGrid::m_multiSolverBoundingBox</code>\n
       default = <code>none</code>\n \n
       Defines the global bounding box by two points for a (multisolver) grid. To be used in
       combination with 'multiSolverMinLevel'. Specifying both parameters allows to create a
       corresponding single solver grid that is ordered by the same Hilbert curve such that single
       solver restart/mean-vars files can be used in a coupled multisolver simulation.\n\n
       Keywords: <i>GRID, HILBERT, MULTISOLVER, COUPLING</i>
    */
    if(Context::propertyLength("multiSolverBoundingBox") != 2 * nDim) {
      TERMM(1, "multiSolverBoundingBox: wrong number of coordinates for bouding box "
                   + to_string(Context::propertyLength("multiSolverBoundingBox")));
    }

    m_multiSolverBoundingBox.resize(2 * nDim);
    for(MInt i = 0; i < 2 * nDim; i++) {
      m_multiSolverBoundingBox[i] = Context::getBasicProperty<MFloat>("multiSolverBoundingBox", AT_, i);
    }

    if(Context::propertyExists("multiSolverMinLevel")) {
      /*! \page propertyPage1
         \section multiSolverMinLevel
         <code>MInt CartesianGrid::m_multiSolverMinLevel</code>\n
         default = <code>none</code>\n \n
         Defines the global min-level for a (multisolver) grid. See 'multiSolverBoundingBox'.
         Keywords: <i>GRID, HILBERT, MULTISOLVER, COUPLING</i>
      */
      m_multiSolverMinLevel = Context::getBasicProperty<MInt>("multiSolverMinLevel", AT_);

      if(m_multiSolverMinLevel < 0) {
        TERMM(1, "ERROR: invalid multiSolverMinLevel " + std::to_string(m_multiSolverMinLevel));
      }
    } else {
      m_multiSolverMinLevel = m_minLevel;
    }

    if(m_targetGridFileName != "") {
      TERMM(1, "ERROR: multiSolverBoundingBox and targetGridFileName specified.");
    }
  }

  // Get output directory
  MString testcaseDir = "./";
  /*! \page propertyPage1
    \section testcaseDir
    <code>MString testcaseDir </code>\n
    default = <code>"./"</code>\n \n
    Main directory where MAIA runs.\n
    Possible values are:
    <ul>
      <li>any valid directory name</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, DIRECTORY</i>
   */
  testcaseDir = Context::getBasicProperty<MString>("testcaseDir", AT_, &testcaseDir);

  /*! \page propertyPage1
    \section outputDir
    <code>MString GridgenPar::m_outputDir </code>\n
    default = <code>n/a</code>\n \n
    Directory where the generated grid is placed (relative to `testcaseDir`).\n
    Possible values are:
    <ul>
      <li>any valid directory name</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, DIRECTORY</i>
   */
  m_outputDir = testcaseDir + Context::getBasicProperty<MString>("outputDir", AT_);

  struct stat info;
  // check if dir is valid
  if(stat(m_outputDir.c_str(), &info) != 0) {
    std::stringstream ss;
    ss << "Output directory is invalid: ";
    ss << m_outputDir.c_str();
    mTerm(1, AT_, ss.str());
  }

  /*! \page propertyPage1
    \section checkGridLbValidity
    <code>MString GridgenPar::m_checkGridLbValidity</code>\n
    default = <code>true</code>\n \n
    Enable/disable the grid LB validity check in the grid generator to detect errorneous cells for LB computations.\n
    Keywords: <i>GRID, GENERATOR, PARALLEL, VALIDITY, CHECK, LB,</i>
   */
  m_checkGridLbValidity = true;
  m_checkGridLbValidity = Context::getBasicProperty<MBool>("checkGridLbValidity", AT_, &m_checkGridLbValidity);

  /*! \page propertyPage1
    \section ggp_keepOutsideBndryCellChildren
    <code>MInt GridgenPar::m_keepOutsideBndryCellChildren</code>\n
    default:  0 \n \n
      <ul>
     <li>0 (delets all outside boundary children)</li>
     <li>1 (keep outside boundary children, if the outside-cell's parent has a boundary-cell-neighbor without
    children)</li> <li>2 (keep outside boundary children, if the outside-cell's parent is a boundary-cell)</li>
     </ul>
    Triggers whether and which outside boundary children are keept or deleted. \n
    Keywords: <i>FINITE_VOLUME, GRID GENERATION, BOUNDARY</i>
  */
  m_keepOutsideBndryCellChildren = 0;
  m_keepOutsideBndryCellChildren =
      Context::getBasicProperty<MInt>("ggp_keepOutsideBndryCellChildren", AT_, &m_keepOutsideBndryCellChildren);

  if(m_keepOutsideBndryCellChildren == 3) {
    if(m_reductionFactor <= 1.0) {
      mTerm(1, AT_, "ERROR: m_reductionFactor must be > 1 to use keepOutsideBndryCellChildren in mode 3");
    }
    if(noDomains() > 1) {
      mTerm(1, AT_, "ERROR: keepOutsideBndryCellChildren in mode 3 is not implemented for parallel usage yet!!!");
    }
    mAlloc(m_noSolidLayer, m_noSolvers, "noSolidLayer", AT_);
    for(MInt s = 0; s < m_noSolvers; s++) {
      m_noSolidLayer[s] = 1;
      m_noSolidLayer[s] = Context::getSolverProperty<MInt>("noSolidLayer", s, AT_, &m_noSolidLayer[s]);
      if(domainId() == 0) {
        cerr << "Applying " << m_noSolidLayer[s] << " solid layers for solver " << s << endl;
      }
    }
  }

  RECORD_TIMER_STOP(m_t_readProperties);
}

// how does propertyExists work with solver values?
template <MInt nDim>
void GridgenPar<nDim>::readSolverProperties(MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  // the highest level of definite uniform refinement
  if(Context::propertyExists("maxGeometricalRfnLvl", solver)) {
    cerr << "Warning: Property maxGeometricalRfnLvl is deprecated. Please rename to maxUniformRefinementLevel." << endl;
    bp->maxUniformRefinementLevel = Context::getSolverProperty<MInt>("maxGeometricalRfnLvl", solver, AT_);
  } else {
    bp->maxUniformRefinementLevel = Context::getSolverProperty<MInt>("maxUniformRefinementLevel", solver, AT_);
  }

  MInt maxRequestedRfnmntLvl = bp->maxUniformRefinementLevel;

  m_maxUniformRefinementLevel = mMin(bp->maxUniformRefinementLevel, m_maxUniformRefinementLevel);

  // the highest level of refinement
  bp->maxRfnmntLvl = Context::getSolverProperty<MInt>("maxRfnmntLvl", solver, AT_);

  // highest level of boundary refinement
  bp->maxBoundaryRfnLvl = Context::getSolverProperty<MInt>("maxBoundaryRfnLvl", solver, AT_);

  // sanity check
  stringstream errorMsg;
  if(m_minLevel > bp->maxUniformRefinementLevel) {
    errorMsg << "ERROR: minLevel is larger than maxUniformRefinementLevel of solver " << solver << endl;
  }
  if(m_minLevel > bp->maxRfnmntLvl) {
    errorMsg << "ERROR: minLevel is larger than maxRfnmntLvl" << endl;
  }
  if(bp->maxUniformRefinementLevel > bp->maxRfnmntLvl) {
    errorMsg << "ERROR: maxUniformRefinementLevel is larger than maxRfnmntLvl" << endl;
  }
  if(bp->maxBoundaryRfnLvl > bp->maxRfnmntLvl) {
    errorMsg << "ERROR: maxBoundaryRfnLvl is larger than maxRfnmntLvl" << endl;
  }

  if(errorMsg.str() != "") {
    m_log << errorMsg.str();
    mTerm(1, AT_, errorMsg.str());
  }

  // check local refinement properties
  // 0: no local refinement
  // 1: patch refinement
  // 2: boundary refinement
  // 3: patch and boundary refinement
  /*! \page propertyPage1
    \section localRfnMethods
    <code>MInt GridgenPar::localRfnMethods</code>\n
    default = <code>0</code>\n \n
    Sets the local refinement method \n
     <ul>
     <li>0 (no local refinement)</li>
     <li>1 (patch refinement)</li>
     <li>2 (boundary refinement)</li>
     <li>3 (patch and boundary refinement)</li>
     </ul>
    Keywords: <i>FINITE_VOLUME, GRID GENERATION</i>
  */

  MInt localRfnMethod = 0;
  localRfnMethod = Context::getSolverProperty<MInt>("localRfnMethod", solver, AT_, &localRfnMethod);

  if((localRfnMethod == 1 || localRfnMethod == 3) && (bp->maxRfnmntLvl > bp->maxUniformRefinementLevel)) {
    /*! \page propertyPage1
      \section localRfnLvlMethods
      <code>MInt GridgenPar::localRfnLevelMethods</code>\n
      default = <code>0</code>\n \n
      Sets the local refinement level method \n
       <ul>
       <li></li>
       </ul>
      Keywords: <i>FINITE_VOLUME, GRID GENERATION</i>
    */
    MString localRfnLevelMethods = Context::getSolverProperty<MString>("localRfnLvlMethods", solver, AT_);

    /*! \page propertyPage1
      \section weightPatchCells
      <code>MInt GridGenPar::m_weightPatchCells </code>\n
      default = <code>0</code>\n \n
      Controls the initial static load balancing of patch cells.
      It improves the initial cell distribution.
      <ul>
      <li>0 - Inactive</li>
      <li>1 - Active</li>
      </ul>
      Keywords: <i>GRID, GENERATOR, PARALLEL, LOAD BALANCING</i>
     */
    m_weightPatchCells = mMax(m_weightPatchCells,
                              Context::getSolverProperty<MInt>("weightPatchCells", solver, AT_, &m_weightPatchCells));

    // sanity check
    if(localRfnLevelMethods.front() == '-' || localRfnLevelMethods.back() == '-')
      mTerm(1, AT_, "ERROR: localRfnLvlMethods begins or ends with hyphen!");

    if(localRfnLevelMethods.find(" ") != MString::npos) mTerm(1, AT_, "ERROR: localRfnLvlMethods contains space!");

    MString::size_type prev_pos = 0;
    MString::size_type nextMarker = 0;

    // save substrings marked with '-' as the patch-refinement methods
    while(nextMarker != MString::npos) {
      if(prev_pos > 0) {
        prev_pos++;
      }
      nextMarker = localRfnLevelMethods.find("-", prev_pos);
      MString substring(localRfnLevelMethods.substr(prev_pos, nextMarker - prev_pos));
      bp->localRfnLevelMethods.push_back(substring);
      prev_pos = nextMarker;
    }

    MBool reducedNoPatches = false;
    // limit to bp->maxRfnmntLvl
    if(bp->noLocalPatchRfnLvls() + bp->maxUniformRefinementLevel > bp->maxRfnmntLvl) {
      cerr << "WARNING: highest patch level exceeds maxRfnmntLvl in solver " << solver
           << ", higher patches will be omitted!" << endl;
      bp->localRfnLevelMethods.resize(bp->maxRfnmntLvl - bp->maxUniformRefinementLevel);
      reducedNoPatches = true;
    }

    // how many patches do we have
    MInt numPatches = 0;
    for(MInt l = 0; l < bp->noLocalPatchRfnLvls(); l++)
      numPatches += bp->noPatchesPerLevel(l);

    // allocate vector containing the offset to each rfnLvl in the patchProperties matrix
    mAlloc(bp->localRfnLevelPropertiesOffset, bp->noLocalPatchRfnLvls() + 1, "localRfnLevelPropertiesOffset", 0, AT_);
    // allocate vector containing the number of properties for each patch
    mAlloc(bp->noLocalRfnPatchProperties, numPatches, "noLocalRfnPatchProperties", 0, AT_);

    // calculate number of properties required for each patchRfnLvl and set offsets
    bp->localRfnLevelPropertiesOffset[0] = 0;
    MInt count_req = 0;
    MInt s = 0;
    for(MInt i = 0; i < bp->noLocalPatchRfnLvls(); i++) {
      const MString lvlStr = bp->localRfnLevelMethods[i];
      for(MString::size_type j = 0; j < lvlStr.size(); j++) {
        const MString patchStr = lvlStr.substr(j, 1);
        if(patchStr == "B") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim;
        } else if(patchStr == "R") {
          bp->noLocalRfnPatchProperties[s] = nDim + 1;
        } else if(patchStr == "C") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 1;
        } else if(patchStr == "T") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 2;
        } else if(patchStr == "O") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 2;
        } else if(patchStr == "H") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 3;
        } else if(patchStr == "S") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 2;
        } else if(patchStr == "W") {
          bp->noLocalRfnPatchProperties[s] = nDim + 5;
        } else if(patchStr == "F") {
        } else if(patchStr == "A") {
          bp->noLocalRfnPatchProperties[s] = 2 * nDim + 2;
        } else if(patchStr == "N") {
          bp->noLocalRfnPatchProperties[s] = 4 * nDim + 2;
        } else {
          TERMM(1, "Unknown patch type: '" + patchStr + "'");
        }
        count_req += bp->noLocalRfnPatchProperties[s];
        s++;
      }
      bp->localRfnLevelPropertiesOffset[i + 1] = s;
    }

    // allocate matrix containing the properties for each patch
    mAlloc(bp->localRfnPatchProperties, numPatches, bp->noLocalRfnPatchProperties, "localRfnPatchProperties", AT_);


    /*! \page propertyPage1
     \section localRfnLevelProperties
     <code>MInt gridgenpar::property</code>\n
     default = <code>no default value</code>\n \n
     List of floating numbers to define the parameters of the chosen localRfnLvlMethods.
     For example for BR-B one would need to give the following list of values:
     corner_point1 + corner_point2 + center_point + radius + corner_point3 + corner_point4
     Possible values are:
     <ul>
       <li>Coordinates of a box patch (xmin, ymin, zmin, xmax, ymax, zmax).</li>
     </ul>
     Keywords: <i>Grid generation</i>
   */
    MInt count_prov = Context::propertyLength("localRfnLevelProperties", solver);

    // sanity check
    if(count_prov < count_req || (count_prov != count_req && !reducedNoPatches)) {
      mTerm(1, AT_,
            "ERROR: number of localRfnLevelProperties does not match the requested value! " + std::to_string(count_prov)
                + " != " + std::to_string(count_req));
    }

    // load properties of patch refinement methods
    MInt lvl = 0;
    MInt j = 0;
    for(MInt i = 0; i < count_req; i++) {
      bp->localRfnPatchProperties[lvl][j] =
          Context::getSolverProperty<MFloat>("localRfnLevelProperties", solver, AT_, i);

      j++;
      if(j == bp->noLocalRfnPatchProperties[lvl]) {
        j = 0;
        lvl++;
      }
    }

    maxRequestedRfnmntLvl = mMax(bp->maxUniformRefinementLevel + bp->noLocalPatchRfnLvls(), maxRequestedRfnmntLvl);
  } else {
    if(localRfnMethod == 1 || localRfnMethod == 3) {
      cerr << "WARNING: local patch refinement requested (localRfnMethod = 1 or 3),"
           << " but maxRfnmntLvl is equal to maxUniformRefinementLevel";
    }
  }

  bp->noLocalBndRfnLvls = 0;
  if((localRfnMethod == 2 || localRfnMethod == 3) && (bp->maxBoundaryRfnLvl > bp->maxUniformRefinementLevel)) {
    /*! \page propertyPage1
      \section weightBndCells
      <code>MInt gridgenpar::m_weightBndCells</code>\n
      default = <code>0</code>\n \n
      If true, compute cell weights and the sum of the weights in ParallelizeGrid()
      Possible values are:
      <ul>
        <li>0 or 1</li>
      </ul>
      Keywords: <i>Grid generation</i>
    */
    m_weightBndCells =
        mMax(m_weightBndCells, Context::getSolverProperty<MInt>("weightBndCells", solver, AT_, &m_weightBndCells));

    bp->noLocalBndRfnLvls = bp->maxBoundaryRfnLvl - bp->maxUniformRefinementLevel;
    /*! \page propertyPage1
      \section localBndRfnMethod
      <code>MInt gridgenpar::m_localBndRfnMethod</code>\n
      default = <code>0</code>\n \n
      the method to determine distance of all levels in markComputationalGridBndRfn
      Possible values are:
      <ul>
      <li> 0: refinment distance = level_0*localMinBoundaryThreshold. smoothDistance based on level_0 </li>
      <li> 1: refinment distance = level_0/localMinBoundaryThreshold. smoothDistance based on current lvl </li>
      <li> 2: localBndRfnDistance represents refinement distance. smoothDistance based on current level </li>
      </ul>
      Keywords: <i>Grid generation</i>
    */
    // how do i want to prescribe the layers of refinement
    bp->localBndRfnMethod = 0;
    bp->localBndRfnMethod =
        Context::getSolverProperty<MInt>("localBndRfnMethod", solver, AT_, &(bp->localBndRfnMethod));


    bp->noLocalRfnBoundaryIds = Context::propertyLength("localRfnBoundaryIds", solver);

    if(!(bp->noLocalRfnBoundaryIds > 0))
      mTerm(1, AT_, "ERROR: localRfnBoundaryIds was requested but is not in property file");

    mAlloc(bp->localRfnBoundaryIds, bp->noLocalRfnBoundaryIds, "locaRfnBoundaryIds", 0, AT_);

    for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++) {
      bp->localRfnBoundaryIds[i] = Context::getSolverProperty<MInt>("localRfnBoundaryIds", solver, AT_, i);
    }

    /*! \page propertyPage1
     \section localMinBoundaryThreshold
     <code> MInt GridGenPar::m_localMinBoundaryThreshold </code>  \n
     default = <code>""</code>\n \n
     Distance in cell units on the current level, which should be further refined.
     <ul>
     <li>Any positive integer</li>
     </ul>
     Keywords: <i>GRID, GENERATOR, BOUNDARY, REFINEMENT,</i>
    */

    // if there are less entries as in localRfnBoundaryIds, then last value of
    // localMinBoundaryThreshold is default for the missing entries
    MInt noLocalMinBoundaryThreshold = Context::propertyLength("localMinBoundaryThreshold", solver);

    if((noLocalMinBoundaryThreshold > bp->noLocalRfnBoundaryIds) || !(noLocalMinBoundaryThreshold > 0))
      mTerm(1, AT_, "ERROR: localMinBoundaryThreshold has wrong length");

    mAlloc(bp->localMinBoundaryThreshold, bp->noLocalRfnBoundaryIds, "localMinBoundaryThreshold",
           Context::getSolverProperty<MInt>("localMinBoundaryThreshold", solver, AT_, noLocalMinBoundaryThreshold - 1),
           AT_);

    for(MInt i = 0; i < noLocalMinBoundaryThreshold; i++) {
      bp->localMinBoundaryThreshold[i] = Context::getSolverProperty<MInt>("localMinBoundaryThreshold", solver, AT_, i);
    }


    if(bp->localBndRfnMethod == 2) {
      /*! \page propertyPage1
       \section localBndRfnDistance
       default = <code>""</code>\n \n
       Distance in STL units to refine first layer around boundary.
       <ul>
       <li>Array of positive floating point numbers. Last one is taken as default value if entries are missing.</li>
       </ul>
       Keywords: <i>GRID, GENERATOR, BOUNDARY, REFINEMENT,</i>
      */
      if(Context::propertyExists("localBndRfnDistance", solver)) {
        const MInt noLocalBndRfnDistance = Context::propertyLength("localBndRfnDistance", solver);
        mAlloc(bp->localBndRfnDistance, bp->noLocalRfnBoundaryIds, "localBndRfnDistance",
               Context::getSolverProperty<MFloat>("localBndRfnDistance", solver, AT_, noLocalBndRfnDistance - 1), AT_);
        for(MInt i = 0; i < noLocalBndRfnDistance; i++) {
          bp->localBndRfnDistance[i] = Context::getSolverProperty<MFloat>("localBndRfnDistance", solver, AT_, i);
        }
      } else {
        mTerm(1, AT_, "ERROR: localBndRfnMethod == 2 requires localBndRfnDistance property");
      }
      /*! \page propertyPage1
       \section localBndRfnMinLvlDiff
       default = <code>""</code>\n \n
       Min target level to refine down from boundary level. Useful to avoid
       further refinement than p.e. provided by already existing patch.
       <ul>
       <li>Positive integer number.</li>
       </ul>
       Keywords: <i>GRID, GENERATOR, BOUNDARY, REFINEMENT,</i>
      */
      if(Context::propertyExists("localBndRfnMinLvlDiff", solver)) {
        const MInt noLocalBndRfnMinLvlDiff = Context::propertyLength("localBndRfnMinLvlDiff", solver);
        mAlloc(bp->localBndRfnMinLvlDiff, bp->noLocalRfnBoundaryIds, "localBndRfnMinLvlDiff",
               Context::getSolverProperty<MInt>("localBndRfnMinLvlDiff", solver, AT_, noLocalBndRfnMinLvlDiff - 1),
               AT_);
        for(MInt i = 0; i < noLocalBndRfnMinLvlDiff; i++) {
          bp->localBndRfnMinLvlDiff[i] = Context::getSolverProperty<MInt>("localBndRfnMinLvlDiff", solver, AT_, i);
        }
      } else {
        // if property is not given, use default value to refine down to min level
        mAlloc(bp->localBndRfnMinLvlDiff, bp->noLocalRfnBoundaryIds, "localBndRfnMinLvlDiff", bp->noLocalBndRfnLvls,
               AT_);
      }
    }


    /*! \page propertyPage1
      \section smoothDistance
      <code>MInt GridGenPar::property </code>\n
      default = <code>1</code>\n \n
      The property defines how many cells should be inbetween the refinement levels.
      <ul>
      <li>Integer values larged than 1.</li>
      </ul>
      Keywords: <i>GRID, CUTOFF</i>
    */
    // if smoothDistance is not in property file, then '1' is set for all ids
    // if there are less entries as in localRfnBoundaryIds, then last value of
    // smoothDistance is default for the missing entries
    if(Context::propertyExists("smoothDistance", solver)) {
      MInt noSmoothDistance = Context::propertyLength("smoothDistance", solver);

      if(noSmoothDistance <= bp->noLocalRfnBoundaryIds && noSmoothDistance > 0) {
        mAlloc(bp->smoothDistance, bp->noLocalRfnBoundaryIds, "smoothDistance",
               Context::getSolverProperty<MInt>("smoothDistance", solver, AT_, noSmoothDistance - 1), AT_);

        for(MInt i = 0; i < noSmoothDistance; i++) {
          bp->smoothDistance[i] = Context::getSolverProperty<MInt>("smoothDistance", solver, AT_, i);
        }
      } else {
        mTerm(1, AT_, "ERROR: smoothDistance is wrongly defined in property file");
      }
    } else {
      mAlloc(bp->smoothDistance, bp->noLocalRfnBoundaryIds, "smoothDistance", 1, AT_);
    }


    /*! \page propertyPage1
      \section localRfnLvlDiff
      <code>MInt GridGenPar::m_localRfnLvlDiff </code>\n
      default = <code>[0,0,...]</code>\n \n
      Level difference between maxBoundaryRfnLvl and the level at the
      corresponding entry in localRfnBoundaryIds. Last value holds for
      non-specified level differences.
      Keywords: <i>GRID, BOUNDARY REFINEMENT</i>
    */

    // if localRfnLvlDiff is not in property file, then '0' is set for all ids
    // if there are less entries as in localRfnBoundaryIds, then last value of
    // smoothDistance is default for the missing entries
    if(Context::propertyExists("localRfnLvlDiff", solver)) {
      MInt noLocalRfnLvlDiff = Context::propertyLength("localRfnLvlDiff", solver);
      if(noLocalRfnLvlDiff <= bp->noLocalRfnBoundaryIds && noLocalRfnLvlDiff > 0) {
        mAlloc(bp->localRfnLvlDiff, bp->noLocalRfnBoundaryIds, "localRfnLvlDiff",
               Context::getSolverProperty<MInt>("localRfnLvlDiff", solver, AT_, noLocalRfnLvlDiff - 1), AT_);

        for(MInt i = 0; i < noLocalRfnLvlDiff; i++)
          bp->localRfnLvlDiff[i] = Context::getSolverProperty<MInt>("localRfnLvlDiff", solver, AT_, i);
      } else {
        mTerm(1, AT_, "ERROR: localRfnLvlDiff is wrongly defined in property file");
      }
    } else {
      mAlloc(bp->localRfnLvlDiff, bp->noLocalRfnBoundaryIds, "localRfnLvlDiff", 0, AT_);
    }

    MInt minLocalRfnLvlDiff = numeric_limits<MInt>::max();
    for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++) {
      minLocalRfnLvlDiff = mMin(minLocalRfnLvlDiff, bp->localRfnLvlDiff[i]);
    }

    if(minLocalRfnLvlDiff < 0) mTerm(1, AT_, "ERROR: localRfnLvlDiff has negative values");

    if(bp->maxRfnmntLvl - minLocalRfnLvlDiff <= bp->maxUniformRefinementLevel) {
      cerr << "WARNING: localRfnLvlDiff of solver " << solver << " leads to no further refinement"
           << "         boundary refinement is disabled" << endl;
      bp->noLocalBndRfnLvls = 0;
    }

    maxRequestedRfnmntLvl = mMax(bp->maxBoundaryRfnLvl - minLocalRfnLvlDiff, maxRequestedRfnmntLvl);
  } else {
    if(localRfnMethod == 2 || localRfnMethod == 3) {
      cerr << "WARNING: local boundary refinement requested (localRfnMethod = 2 or 3),"
           << " but maxBoundaryRfnLvl is equal to maxUniformRefinementLevel";
    }
  }

  m_maxRfnmntLvl = mMax(m_maxRfnmntLvl, maxRequestedRfnmntLvl);

  /*! \page propertyPage1
  \section cutOff
  <code>MBool GridGenPar::m_cutOff </code>\n
  default = <code>0</code>\n \n
  defines if a cutoff is applied to the grid and if so, at which levels
  <ul>
  <li> 0 - no cutOff  </li>
  <li> 1 - cutOff at initialLevel </li>
  <li> 2 - cutOff at levels <= initialLevel </li>
  <li> 3 - cutOff at levels >= initialLevel</li>
  <li> 34- cutOff at all levels </li>
  </ul>
  Keywords: <i>GRID, CUTOFF</i>
 */
  bp->cutOff = 0;
  bp->cutOff = Context::getSolverProperty<MInt>("cutOff", solver, AT_, &(bp->cutOff));

  if(bp->cutOff) {
    m_cutOff = true;

    /*! \page propertyPage1
      \section cutOffMethod
      <code>MInt GridGenPar::readSolverProperties::cutOffMethod </code>\n
      default = <code>""</code>\n \n
      Specifies the cutoff method to be used.\n
      <ul>
      <li> B - box cut-off  </li>
      <li> iB- inverse box cut-off (keep everything outside box)  </li>
      <li> P - plane cut-off  </li>
      <li> C - cylindrical slice cut-off (azimuthal Periodicity) </li>
      </ul>
      Keywords: <i>GRID, CUTOFF</i>
    */
    MString cutOffMethod = "";
    cutOffMethod = Context::getSolverProperty<MString>("cutOffMethod", solver, AT_, &cutOffMethod);

    if(cutOffMethod == "") {
      errorMsg << "ERROR: cutOffMethod is not defined in property file for solver " << solver << endl;
      m_log << errorMsg.str();
      mTerm(1, AT_, errorMsg.str());
    }

    MString::size_type prev_pos = 0, pos = 0;

    while((pos = cutOffMethod.find("-", pos)) != MString::npos) {
      MString substring(cutOffMethod.substr(prev_pos, pos - prev_pos));
      bp->cutOffMethods.push_back(substring);
      prev_pos = ++pos;
    }

    bp->cutOffMethods.push_back(cutOffMethod.substr(prev_pos, pos - prev_pos));

    MIntScratchSpace nmbrValuesPerMethod((signed)bp->cutOffMethods.size(), AT_, "nmbrValuesPerMethod");
    MInt* tmpPointer = nmbrValuesPerMethod.begin();

    MInt count_req = 0;
    for(MInt i = 0; i < (signed)bp->cutOffMethods.size(); i++) {
      if(bp->cutOffMethods[i] == "P") {
        nmbrValuesPerMethod[i] = 2 * nDim;
      } else if(bp->cutOffMethods[i] == "B") {
        nmbrValuesPerMethod[i] = 2 * nDim;
      } else if(bp->cutOffMethods[i] == "iB") {
        nmbrValuesPerMethod[i] = 2 * nDim;
      } else if(bp->cutOffMethods[i] == "C") {
        nmbrValuesPerMethod[i] = 4;
      }
      count_req += nmbrValuesPerMethod[i];
    }

    // read properties
    /*! \page propertyPage1
      \section cutOffCoordinates
      <code>MFloat GRIDGENPAR::cutOffCoordinates</code>\n
      no default \n \n
      Coordinates for the cut off (box,plane) \n
      Keywords: <i>GRID</i>
    */
    MInt count_prov = Context::propertyLength("cutOffCoordinates", solver);
    if(count_req == count_prov) {
      mAlloc(bp->cutOffCoordinates, (signed)bp->cutOffMethods.size(), tmpPointer, "cutOffCoordinates", AT_);

      MInt pos2 = 0;
      for(MInt i = 0; i < (signed)bp->cutOffMethods.size(); i++) {
        for(MInt j = 0; j < nmbrValuesPerMethod[i]; j++, pos2++) {
          bp->cutOffCoordinates[i][j] = Context::getSolverProperty<MFloat>("cutOffCoordinates", solver, AT_, pos2);
        }
      }
    } else {
      errorMsg << "ERROR: cutOffCoordinates are not defined or not properly defined in property file for solver "
               << solver << endl;
      m_log << errorMsg.str();
      mTerm(1, AT_, errorMsg.str());
    }


    MIntScratchSpace nmbrLayersPerMethod((signed)bp->cutOffMethods.size(), AT_, "nmbrLayersPerMethod");
    tmpPointer = nmbrLayersPerMethod.begin();

    count_req = 0;
    for(MInt i = 0; i < (signed)bp->cutOffMethods.size(); i++) {
      if(bp->cutOffMethods[i] == "P") {
        nmbrLayersPerMethod[i] = 1;
      } else if(bp->cutOffMethods[i] == "B") {
        nmbrLayersPerMethod[i] = 2 * nDim;
      } else if(bp->cutOffMethods[i] == "iB") {
        nmbrLayersPerMethod[i] = 2 * nDim;
      } else if(bp->cutOffMethods[i] == "C") {
        nmbrLayersPerMethod[i] = 1;
      }
      count_req += nmbrLayersPerMethod[i];
    }

    /*! \page propertyPage1
      \section cutOffCoordinates
      <code>MFloat GRIDGENPAR::cutOffNmbrLayers</code>\n
      default: 0 \n \n
      Allow for additional layers around the cutOff. These additional cells are remaining in the grid!
      This is very useful for all zonal methods!
      \n
      Keywords: <i>GRID, CUTOFF, ZONAL</i>
    */
    mAlloc(bp->cutOffNmbrLayers, (signed)bp->cutOffMethods.size(), tmpPointer, "cutOffNmbrLayers", AT_);
    if(Context::propertyExists("cutOffNmbrLayers", solver)) {
      count_prov = Context::propertyLength("cutOffNmbrLayers", solver);
      if(count_prov == count_req) {
        MInt pos2 = 0;
        for(MInt i = 0; i < (signed)bp->cutOffMethods.size(); i++) {
          for(MInt j = 0; j < nmbrLayersPerMethod[i]; j++, pos2++) {
            bp->cutOffNmbrLayers[i][j] = Context::getSolverProperty<MInt>("cutOffNmbrLayers", solver, AT_, pos2);
          }
        }
      } else {
        errorMsg << "ERROR: cutOffNmbrLayers are not properly defined in property file for solver " << solver << endl;
        m_log << errorMsg.str();
        mTerm(1, AT_, errorMsg.str());
      }
    } else {
      for(MInt i = 0; i < (signed)bp->cutOffMethods.size(); i++)
        for(MInt j = 0; j < nmbrLayersPerMethod[i]; j++)
          bp->cutOffNmbrLayers[i][j] = 0;
    }
  }
}

/** \brief initializes the member variables
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::initMembers() {
  TRACE();

  RECORD_TIMER_START(m_t_initMembers);

  m_log << "  (2) initializing member variables" << endl;
  outStream << "  (2) initializing member variables" << endl;

  mAlloc(m_boundingBox, 2 * nDim, "m_boundingBox", AT_);
  mAlloc(m_geometryExtents, nDim, "m_geometryExtents", AT_);
  mAlloc(m_centerOfGravity, nDim, "m_centerOfGravity", AT_);

  mAlloc(m_cells, m_maxNoCells, nDim, 0, "m_cells", AT_);
  m_pCells = m_cells->a;

  m_noCells = -1;
  m_noTotalCells = -1;
  m_noPartitionCells = -1;
  m_noTotalPartitionCells = -1;
  m_noTotalHaloCells = -1;
  m_decisiveDirection = 0;
  m_noNeighborDomains = 0;
  m_maxNoChildren = IPOW2(nDim);
  m_maxLevels = 31;
  m_noNeighbors = 2 * nDim;
  m_cellOffsetPar = 0;

  // the offsets of the levels
  mAlloc(m_levelOffsets, m_maxRfnmntLvl + 1, 2, "m_levelOffsets", -1, AT_);
  mAlloc(m_haloCellOffsetsLevel, m_maxRfnmntLvl + 1, 2, "m_haloCellOffsetsLevel", -1, AT_);
  mAlloc(m_noHaloCellsOnLevel, m_maxRfnmntLvl + 1, "m_noHaloCellsOnLevel", 0, AT_);
  mAlloc(m_noCellsPerDomain, noDomains(), "m_noCellsPerDomain", 0, AT_);

  // the lengths on the different levels
  mAlloc(m_lengthOnLevel, m_maxLevels, "m_lengthOnLevel", 0.0, AT_);

  if(m_minLevel % 2 == 0) {
    m_levelOffsets[0][0] = 0;
    m_levelOffsets[0][1] = 1;
    m_levelOffsets[1][0] = m_maxNoCells - m_maxNoChildren;
    m_levelOffsets[1][1] = m_maxNoCells;
  } else {
    m_levelOffsets[0][0] = m_maxNoCells - 1;
    m_levelOffsets[0][1] = m_maxNoCells;
    m_levelOffsets[1][0] = 0;
    m_levelOffsets[1][1] = m_maxNoChildren;
  }

  // dynamic load balancing
  m_noMissingParents = 0;
  m_hasBeenLoadBalanced = false;

  RECORD_TIMER_STOP(m_t_initMembers);
}

namespace {

// Auxiliary type traits to support selecting the correct geometry class
template <MInt nDim>
struct GeometryXD {};
template <>
struct GeometryXD<2> {
  using type = Geometry2D;
};
template <>
struct GeometryXD<3> {
  using type = Geometry3D;
};

} // namespace

/** \brief initializes the geometry
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * - Creates a new Geometry-object depending on the dimensionality of the problem.
 * - Prints out further information about the geometry
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::initGeometry() {
  TRACE();

  RECORD_TIMER_START(m_t_initGeometry);

  m_log << "  (3) initializing geometry" << endl;
  outStream << "  (3) initializing geometry" << endl;

  m_geometry = new GeometryRoot(m_noSolvers, nDim, mpiComm());

  mAlloc(m_noBndIdsPerSolver, m_geometry->noNodes(), AT_, 0, "noBndIds");
  for(MInt solver = 0; solver < m_geometry->noNodes(); solver++) {
    m_noBndIdsPerSolver[solver] = m_geometry->noSegmentsOfNode(solver);
  }
  mAlloc(m_bndCutInfo, m_geometry->noNodes(), m_noBndIdsPerSolver, AT_, "bndCutInfo");

  m_geometry->boundingBox(m_boundingBox);


  // Use the given multisolver bounding box as bounding box for the grid generation (and not only for
  // sorting by a different Hilbert curve if the given min-levels are different)
  if(m_hasMultiSolverBoundingBox && m_multiSolverMinLevel == m_minLevel) {
    m_log << "Using multisolver bounding box information from property file" << std::endl;
    for(MInt i = 0; i < nDim; i++) {
      TERMM_IF_COND(m_boundingBox[i] < m_multiSolverBoundingBox[i],
                    "Multisolver bounding box error (dim: " + std::to_string(i) + "): "
                        + std::to_string(m_boundingBox[i]) + " < " + std::to_string(m_multiSolverBoundingBox[i]));
      TERMM_IF_COND(m_boundingBox[nDim + i] > m_multiSolverBoundingBox[nDim + i],
                    "Multisolver bounding box error (dim: " + std::to_string(i)
                        + "): " + std::to_string(m_boundingBox[nDim + i]) + " > "
                        + std::to_string(m_multiSolverBoundingBox[nDim + i]));
      m_boundingBox[i] = m_multiSolverBoundingBox[i];
      m_boundingBox[nDim + i] = m_multiSolverBoundingBox[nDim + i];
    }
  }


  for(MInt dir = 0; dir < nDim; dir++) {
    m_geometryExtents[dir] = m_boundingBox[dir + nDim] - m_boundingBox[dir];
    m_decisiveDirection = m_geometryExtents[dir] > m_geometryExtents[m_decisiveDirection] ? dir : m_decisiveDirection;
    m_centerOfGravity[dir] = m_boundingBox[dir] + 0.5 * (m_boundingBox[dir + nDim] - m_boundingBox[dir]);
  }

  m_lengthOnLevel[0] = (F1 + F1 / FPOW2(30)) * m_reductionFactor * m_geometryExtents[m_decisiveDirection];
  for(MInt l = 1; l < m_maxLevels; l++)
    m_lengthOnLevel[l] = m_lengthOnLevel[l - 1] * 0.5;

  m_log << "     + center of gravity: ";
  for(MInt dir = 0; dir < nDim; dir++)
    m_log << m_centerOfGravity[dir] << " ";
  m_log << "\n";
  m_log << "     + decisive direction: " << m_decisiveDirection << "\n";
  m_log << "     + geometry extents: ";
  for(MInt dir = 0; dir < nDim; dir++)
    m_log << m_geometryExtents[dir] << " ";
  m_log << "\n";
  m_log << "     + bounding box: ";
  for(MInt dir = 0; dir < m_noNeighbors; dir++)
    m_log << m_boundingBox[dir] << " ";
  m_log << endl;

  outStream << "     + center of gravity: ";
  for(MInt dir = 0; dir < nDim; dir++)
    outStream << m_centerOfGravity[dir] << " ";
  outStream << "\n";
  outStream << "     + decisive direction: " << m_decisiveDirection << "\n";
  outStream << "     + geometry extents: ";
  for(MInt dir = 0; dir < nDim; dir++)
    outStream << m_geometryExtents[dir] << " ";
  outStream << "\n";
  outStream << "     + bounding box: ";
  for(MInt dir = 0; dir < m_noNeighbors; dir++)
    outStream << m_boundingBox[dir] << " ";
  outStream << endl;

  RECORD_TIMER_STOP(m_t_initGeometry);
}


/** \brief aligns the cutOffCoordinates with the grid
 *
 * \author Thomas schilden
 * \date 20.02.2016
 *
 *  aligns the cutOffCoordinates of boxes with the grid on the
 *  m_minLevel
 *  this is necessary for cutOffs at levels higher than
 *  m_minLevel
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::gridAlignCutOff() {
  MFloat cellLength = m_lengthOnLevel[m_minLevel];
  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    SolverRefinement* bp = m_solverRefinement + solver;
    for(MInt c = 0; c < (signed)bp->cutOffMethods.size(); c++) {
      if(!(bp->cutOffMethods[c] == "B" || bp->cutOffMethods[c] == "iB")) continue;
      for(MInt dim = 0; dim < nDim; dim++) {
        for(MInt dir = 0; dir < 2; dir++) {
          MFloat relCoord = bp->cutOffCoordinates[c][dim + dir * nDim] - m_centerOfGravity[dim];
          MInt n = -1;
          if(bp->cutOffMethods[c] == "B") {
            n = (MInt)floor(relCoord / cellLength + 0.5);
          } else if(bp->cutOffMethods[c] == "iB") {
            n = (MInt)ceil(relCoord / cellLength + 0.5);
          }
          bp->cutOffCoordinates[c][dim + dir * nDim] = n * cellLength + m_centerOfGravity[dim];
        }
      }
    }
  }
}


/** \brief creates the initial grid
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * 4.1 Creates an initial cube around the geometry.
 * 4.2 Refines the initial cube to the initialGeometricalRfnLvl.
 *     In each iteration the following is performed:
 *     4.2.1 update the level offsets:
 *           this gets rid of the delettion of lower levels by using the following structure.
 *           Based on the m_initialGeometricalRfnLevel the first initial cube is located at the
 *           beginning or at the end of the collector:
 *             - m_initialGeometricalRfnLevel is even: locate at the beginning
 *             - m_initialGeometricalRfnLevel is odd: locate at the end
 *           The next following level is always written to the other end of the collector, i.e.,
 *           if the first cube is located at the end of the collector, the next higher level starts
 *           at the beginning of the collector. This method iterates, cells are written to the front
 *           back of the collector continuously. This way it is not necessary to delete the lower levels
 *           since they are just overwritten by one of the upcoming levels. This way, at the end of
 *           the refinement the initial level is located at the beginning of the collector. Before refining
 *           and writing to a level at the end of the collector the offsets are determined by estimating the
 *           possible number of cells on the upcoming level.
 *     4.2.2 check if there is enough memory
 *     4.2.3 refine the grid one further level
 *     4.2.4 update the neighborhood information on the newly created cells
 *     4.2.5 update solver affiliation and cells outside the geometry or without solver affiliation are deleted
 * 4.3 Remove the link to the parents on the initial level
 * 4.4 Update the number of generated cells and print information
 * 4.5 Reorder cells after Hilbert id
 * 4.6 If the code runs parallel and no further refinement is requested, finish
 *     the process otherwise nothing is required at this stage. Finishing calls
 *     updateInterRankNeighbors().
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::createInitialGrid() {
  TRACE();

  RECORD_TIMER_START(m_t_createInitialGrid);

  m_log << "  (4) creating initial grid" << endl;
  outStream << "  (4) creating initial grid" << endl;

  m_log << "     + refining from: 0 to " << m_minLevel << endl;
  outStream << "     + refining from: 0 to " << m_minLevel << endl;

  m_log << "     + initial cube size: " << m_lengthOnLevel[0] << endl;
  outStream << "     + initial cube size: " << m_lengthOnLevel[0] << endl;

  MInt in_id = m_levelOffsets[0][0];

  // 4.1 Create initial cube
  for(MInt d = 0; d < nDim; d++)
    a_coordinate(in_id, d) = m_centerOfGravity[d];

  for(MInt d = 0; d < m_noNeighbors; d++)
    a_neighborId(in_id, d) = -1;

  for(MInt c = 0; c < m_maxNoChildren; c++)
    a_childId(in_id, c) = -1;

  a_parentId(in_id) = -1;
  a_level(in_id) = 0;
  a_globalId(in_id) = (MLong)in_id;
  a_noChildren(in_id) = 0;
  for(MInt s = 0; s < m_noSolvers; s++) {
    a_noSolidLayer(in_id, s) = -1;
  }

  for(MInt i = 0; i < 8; i++)
    a_hasProperty(in_id, i) = 0;

  // we need to say that the initial cell has a cut so that the children are checked for a cut
  a_hasProperty(in_id, 1) = 1;

  // all Solvers are cut and inside and to be refined
  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    a_isSolverBoundary(in_id, solver) = 1;
    a_isInSolver(in_id, solver) = 1;
    a_isToRefineForSolver(in_id, solver) = 1;
  }

  m_noCells = 1;

  // 4.2 Refine to the initial refinement level
  for(MInt l = 0; l < m_minLevel; l++) {
    // 4.2.1 Update the level offsets
    if(m_levelOffsets[l][0] == 0) {
      m_levelOffsets[l + 1][0] = m_maxNoCells - (m_levelOffsets[l][1] - m_levelOffsets[l][0]) * m_maxNoChildren;
      m_levelOffsets[l + 1][1] = m_maxNoCells;
    } else {
      m_levelOffsets[l + 1][0] = 0;
      m_levelOffsets[l + 1][1] = (m_levelOffsets[l][1] - m_levelOffsets[l][0]) * m_maxNoChildren;
    }

    // 4.2.2 check if there is enough memory
    checkMemoryAvailability(0, l + 1);

    for(MInt solver = 0; solver < m_noSolvers; solver++)
      markSolverForRefinement(l, solver);

    // 4.2.3 refine
    refineGrid(m_levelOffsets, l, 0);

    // 4.2.4 update neighbors
    m_log << "       * finding the neighbors for the new level" << endl;
    outStream << "       * finding the neighbors for the new level" << endl;
    findChildLevelNeighbors(m_levelOffsets, l);

    // 4.2.5 delete outside cells
    deleteOutsideCellsSerial(l + 1);
  }
  // 4.3 remove parent links on base cell level
  for(MInt i = m_levelOffsets[m_minLevel][0]; i < m_levelOffsets[m_minLevel][1]; i++)
    a_parentId(i) = -1;

  // 4.4 update number of cells generated
  m_noCells = m_levelOffsets[m_minLevel][1] - m_levelOffsets[m_minLevel][0];
  m_log << "     + created " << m_noCells << " cells for level " << m_minLevel << endl;
  outStream << "     + created " << m_noCells << " cells for level " << m_minLevel << endl;

  // 4.5 reorder cells after Hilbert curve
  reorderCellsHilbert();

  RECORD_TIMER_STOP(m_t_createInitialGrid);
}

/** \brief excludes obvious non solver affiliated cells from the inside outside flooding
 *
 * \author Thomas Schilden
 * \date 03.04.2018
 *
 * Before the inside outside check of a specific solver is performed, cells
 * not belonging to the solver are marked as visited, that reduces the number
 * of cells to flood in the following step: markInsideOutside( , ,solver)
 * While doing so, the visited property is reset for cells that might be inside
 *
 * \param[in] level_ the level to check
 * \param[in] solver specific solver
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::excludeInsideOutside(MInt** offsets, MInt level_, MInt solver) {
  TRACE();

  // reset "visited" property
  for(MInt i = offsets[level_][1] - 1; i >= offsets[level_][0]; i--) {
    if(a_isInSolver(i, solver)) {
      a_hasProperty(i, 6) = 0;
    } else {
      a_hasProperty(i, 6) = 1;
    }
  }
}

/** \brief starts the inside-outside flooding for solvers
 *
 * \author Thomas Schilden
 * \date 03.04.2018
 *
 * similar to markInsideOutside(,) but runs only for a specific solver
 *
 * \param[in] offsets: the range of cells to run on
 * \param[in] level_: the level to check
 * \param[in] solver: specific solver
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markInsideOutside(MInt** offsets, MInt level_, MInt solver) {
  TRACE();

  for(MInt i = offsets[level_][1] - 1; i >= offsets[level_][0]; i--) {
    // this cell was already visited, skip
    if(a_hasProperty(i, 6)) continue;

    // this is a solver boundary cell, skip
    if(a_isSolverBoundary(i, solver)) {
      a_isInSolver(i, solver) = 1;
      a_hasProperty(i, 6) = 1;
      continue;
    }

    // cell is inside solver, do inside-solver-flooding
    if(pointIsInsideSolver(&a_coordinate(i, 0), solver)) {
      stack<MInt> fillStack;
      a_isInSolver(i, solver) = 1;
      a_hasProperty(i, 6) = 1;
      fillStack.push(i);
      floodCells(&fillStack, 1, solver);
    } // cell is outside, do outside-flooding
    else {
      stack<MInt> fillStack;
      a_isInSolver(i, solver) = 0;
      a_hasProperty(i, 6) = 1;
      fillStack.push(i);
      floodCells(&fillStack, 0, solver);
    }
  }
}

/** \brief
 *
 * \author
 * \date 16.03.2018
 *
 * \param[in] cell offsets to run on (maybe not needed)
 * \param[in] level to run the check on
 * \param[in] solver to run the check for
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::checkMemoryAvailability(MInt stage, MInt level_) {
  TRACE();

  switch(stage) {
    case 0:
      // the current level is at the beginning of the collector
      if((m_levelOffsets[level_ - 1][0] > m_levelOffsets[level_][0]
          && m_levelOffsets[level_][1] > m_levelOffsets[level_ - 1][0])
         || (m_levelOffsets[level_ - 1][0] < m_levelOffsets[level_][0]
             && m_levelOffsets[level_][0] < m_levelOffsets[level_ - 1][1])) {
        stringstream errorMsg;
        errorMsg << "Not enough memory - normal cells overlap:\n"
                 << "  - cell offsets current level: " << m_levelOffsets[level_][0] << " " << m_levelOffsets[level_][1]
                 << "\n"
                 << "  - cell offset last level: " << m_levelOffsets[level_ - 1][0] << " "
                 << m_levelOffsets[level_ - 1][1] << "\n"
                 << "  - max. no of available cells: " << m_maxNoCells << endl;
        m_log << errorMsg.str();
        mTerm(1, AT_, errorMsg.str());
      }

      break;
    case 1:
      if(noDomains() > 1) {
        if(m_levelOffsets[level_][1] > m_haloCellOffsetsLevel[level_][0]) {
          stringstream errorMsg;
          errorMsg << "Not enough memory - normal and halo cells overlap:\n"
                   << "  - upper normal cell offset: " << m_levelOffsets[level_][1] << "\n"
                   << "  - lower halo cell offset: " << m_haloCellOffsetsLevel[level_][0] << "\n"
                   << "  - max. no of available cells: " << m_maxNoCells << endl;
          m_log << errorMsg.str();
          mTerm(1, AT_, errorMsg.str());
        }
      } else {
        if(m_levelOffsets[level_][1] > m_maxNoCells - 1) {
          stringstream errorMsg;
          errorMsg << "Not enough memory:\n"
                   << "  - upper cell offset: " << m_levelOffsets[level_][1] << "\n"
                   << "  - max. no of available cells: " << m_maxNoCells << endl;
          m_log << errorMsg.str();
          mTerm(1, AT_, errorMsg.str());
        }
      }

      break;
    default:
      break;
  }
}


/** \brief refines the grid on a given level
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * Refines the grid on the level provided by the parameter by running over
 * all cells of the current level. The ids of the cells to be considered are
 * obtained from the arrays m_levelOffsets[level]. Also checks during the run
 * if enough memory is still available.
 * Calls refineCell and provides the id of the cell to refine and the first id to
 * write the newly created cells to.
 *
 * \param[in] offsets the offsets to use in the array of the cells
 * \param[in] level_ the level to be refined
 * \param[in] halo is this a halo refinement?
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::refineGrid(MInt** offsets, MInt level_, MBool halo) {
  TRACE();

  MInt t_refineGrid = 0;
  if(level_ < m_minLevel) {
    NEW_SUB_TIMER_STATIC(t_rfnGrid, "refine Grid", m_t_createInitialGrid);
    t_refineGrid = t_rfnGrid;
  } else if(level_ < m_maxUniformRefinementLevel) {
    NEW_SUB_TIMER_STATIC(t_rfnGrid, "refine Grid", m_t_createStartGrid);
    t_refineGrid = t_rfnGrid;
  }

  RECORD_TIMER_START(t_refineGrid);

  if(halo) {
    m_log << "     + refining halo grid on level " << level_ << ":     " << endl;
    outStream << "     + refining halo grid on level " << level_ << ":     " << endl;
  } else {
    m_log << "     + refining grid on level " << level_ << ":     " << endl;
    outStream << "     + refining grid on level " << level_ << ":     " << endl;
  }
  outStream.flush();

  // this is the start where we want to start writing from
  MInt offsetLevelStart = offsets[level_ + 1][0];

  MInt diff = offsets[level_][1] - offsets[level_][0];

  // Refine each cell on the given level
  for(MInt i = offsets[level_][0]; i < offsets[level_][1]; i++) {
    MInt current = i - offsets[level_][0];

    // display progress
    if(diff > 100 && current % (diff / 100) == 0) {
      if(100 * current / diff < 11)
        outStream << "\b\b\b" << 100 * current / diff << "% ";
      else
        outStream << "\b\b\b\b" << 100 * current / diff << "% ";
      outStream.flush();
    }

    // refine cell
    if(level_ > m_minLevel && noDomains() > 1) {
      refineCell(m_cellIdLUT[i], &offsetLevelStart);
    } else {
      refineCell(i, &offsetLevelStart);
    }
  }

  outStream << "\b\b\b\b\b 100% done." << endl;

  RECORD_TIMER_STOP(t_refineGrid);
}

/** \brief refines the grid on a given level for a provided patch
 *
 * \author Andreas Lintermann
 * \date 23.10.2013
 *
 * Does the same as refineGrid(...) except that only those cells are
 * refined that apply to the patch constraints (b_properties[5] = 1).
 *
 * \param[in] offsets the offsets to use in the array of the cells
 * \param[in] level the level to be refined
 * \param[in] halo is this a halo refinement?
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::refineGridPatch(MInt** offsets, MInt level_, MBool halo) {
  TRACE();

  NEW_SUB_TIMER_STATIC(t_refineGridPatch, "refine Grid", m_t_createComputationalGrid);
  RECORD_TIMER_START(t_refineGridPatch);

  if(halo) {
    m_log << "     + refining halo grid on level " << level_ << ":     ";
    outStream << "     + refining halo grid on level " << level_ << ":     ";
  } else {
    m_log << "     + refining grid on level " << level_ << ":     ";
    outStream << "     + refining grid on level " << level_ << ":     ";
  }
  outStream.flush();

  // this is the start where we want to start writing from
  MInt startchildId = offsets[level_ + 1][0];

  size_t diff = offsets[level_][1] - offsets[level_][0];

  // Refine each cell on the given level
  for(MInt i = offsets[level_][0]; i < offsets[level_][1]; i++) {
    MInt cellId = i;
    if(level_ > m_minLevel && noDomains() > 1) {
      cellId = m_cellIdLUT[i];
    }
    MInt current = i - offsets[level_][0];
    if(diff > 100 && current % (diff / 100) == 0) {
      if(100 * (size_t)current / diff < 11)
        outStream << "\b\b\b" << 100 * (size_t)current / diff << "% ";
      else
        outStream << "\b\b\b\b" << 100 * (size_t)current / diff << "% ";
      outStream.flush();
    }

    if(a_hasProperty(cellId, 5)) {
      MInt currentNoCells = (offsets[level_][1] - offsets[level_][0]) + m_maxNoChildren;
      if(currentNoCells > m_maxNoCells) {
        stringstream errorMsg;
        errorMsg << "Max. no. cells reached: " << m_maxNoCells << endl;
        m_log << errorMsg.str();
        mTerm(1, AT_, errorMsg.str());
      }
      // refine cell
      refineCell(cellId, &startchildId);
    }
  }

  outStream << "\b\b\b\b\b 100% done." << endl;

  RECORD_TIMER_STOP(t_refineGridPatch);
}

/** \brief refinement > minLevel is processed in createComputationalMultisolverGrid
 *
 * \author Thomas Schilden
 * \date 03.04.2018
 *
 * First, the uniform grid up to the minimal maxUniformRefinementLevel of the solvers
 * is created in createStartGrid.
 * Second, the remaining levels are looped up and the local solver refinement methods
 * are called. They mark cells as to be refined for a specific solver. In concludeSolver
 * refinement, the solver refinement request is concluded to a single property no 5.
 * Then, the offsets are updated and the grid is refined.
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::createComputationalMultisolverGrid() {
  TRACE();

  // could be part of the loop below without any problems
  createStartGrid();

  RECORD_TIMER_START(m_t_createComputationalGrid);

  m_log << " createComputationalMultisolverGrid" << endl;

  // l is the level i want to refine
  for(MInt l = m_maxUniformRefinementLevel; l < m_maxRfnmntLvl; l++) {
    // let all solvers mark the cells they want to refine
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      markLocalSolverRefinement(l, solver);
    }

    // cells that are marked for refinement by at least one solver will be refined
    concludeSolverRefinement(l);

    // update the offsets
    updateOffsets(l);

    if(noDomains() > 1) {
      updateHaloOffsets(l, m_rfnCountHalos, m_rfnCountHalosDom);
    }

    // refine marked cells
    refineComputationalGrid(l);
  }

  // deleting useless solid cells
  if(m_keepOutsideBndryCellChildren == 3) {
    if(noDomains() > 1) {
      mTerm(1, AT_, "ERROR: keepOutsideBndryCellChildren in mode 3 is not implemented for parallel usage yet!!!");
      deleteCoarseSolidCellsParallel();
    } else {
      deleteCoarseSolidCellsSerial();
    }
  }

  RECORD_TIMER_STOP(m_t_createComputationalGrid);
}

template <MInt nDim>
void GridgenPar<nDim>::concludeSolverRefinement(MInt gridLvl) {
  TRACE();

  // reset counts of refined cells
  m_rfnCount = 0;
  m_rfnCountHalos = 0;
  if(noDomains() > 1) {
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      m_rfnCountHalosDom[dom] = 0;
    }
  }

  // set internal cells that need to be refined
  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      if(a_isToRefineForSolver(k, solver) || a_noSolidLayer(k, solver) > 0) {
        m_rfnCount++;
        a_hasProperty(k, 5) = 1;
        break;
      }
    }
  }

  // set halo cells that need to be refined
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(a_isToRefineForSolver(k, solver)) {
            m_rfnCountHalos++;
            m_rfnCountHalosDom[dom]++;
            a_hasProperty(k, 5) = 1;
            break;
          }
        }
      }
    }
  }
}

template <MInt nDim>
void GridgenPar<nDim>::markLocalSolverRefinement(MInt level_, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  if(level_ < bp->maxUniformRefinementLevel) {
    markSolverForRefinement(level_, solver);
  } else {
    MInt refLevel = level_ - bp->maxUniformRefinementLevel;
    if(refLevel < bp->noLocalPatchRfnLvls()) markPatchForSolverRefinement(level_, solver);
    if(refLevel < bp->noLocalBndRfnLvls) markBndForSolverRefinement(level_, solver);
  }
}

template <MInt nDim>
void GridgenPar<nDim>::markPatchForSolverRefinement(MInt level_, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  MInt refLevel = level_ - bp->maxUniformRefinementLevel;

  for(MInt patch = 0; patch < bp->noPatchesPerLevel(refLevel); patch++) {
    MString patchStr = bp->localRfnLevelMethods[refLevel].substr(patch, 1);
    if(patchStr == "B") {
      markLocalBox(refLevel, patch, solver);
    } else if(patchStr == "R") {
      markLocalRadius(refLevel, patch, solver);
    } else if(patchStr == "C" || patchStr == "T") {
      markLocalCylinder(refLevel, patch, solver, patchStr);
    } else if(patchStr == "O") {
      markLocalCone(refLevel, patch, solver);
    } else if(patchStr == "H") {
      markLocalHat(refLevel, patch, solver);
    } else if(patchStr == "S") {
      markLocalRectangleAngled(refLevel, patch, solver);
    } else if(patchStr == "W") {
      markLocalCartesianWedge(refLevel, patch, solver);
    } else if(patchStr == "A" || patchStr == "N") {
      markLocalSlicedCone(refLevel, patch, solver, patchStr);
    } else {
      TERMM(1, "Unknown patch type: '" + patchStr + "'");
    }
  }
}

template <MInt nDim>
void GridgenPar<nDim>::markBndForSolverRefinement(MInt level_, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  MInt refLevel = level_ - bp->maxUniformRefinementLevel;

  m_log << "     + running boundary refinement for solver " << solver << endl;

  // 6.2b.1 initialize the distances of all levels
  MIntScratchSpace dists(bp->noLocalBndRfnLvls, bp->noLocalRfnBoundaryIds, AT_, "dists");
  dists.fill(0);

  switch(bp->localBndRfnMethod) {
    case 0: {
      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        dists(bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1, i) = bp->localMinBoundaryThreshold[i];

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        for(MInt lvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 2; lvl >= 0; lvl--)
          dists(lvl, i) = dists(lvl + 1, i) + bp->smoothDistance[i];

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        for(MInt lvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1; lvl >= 0; lvl--)
          dists(lvl, i) = ((MInt)pow(2, refLevel)) * dists(lvl, i);
    } break;
    case 1: {
      MFloatScratchSpace fists(bp->noLocalBndRfnLvls, bp->noLocalRfnBoundaryIds, AT_, "fists");
      fists.fill(F0);

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        fists(bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1, i) =
            m_lengthOnLevel[0] / m_reductionFactor / bp->localMinBoundaryThreshold[i];

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        for(MInt lvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 2; lvl >= 0; lvl--)
          fists(lvl, i) =
              fists(lvl + 1, i) + bp->smoothDistance[i] * m_lengthOnLevel[bp->maxUniformRefinementLevel + lvl + 1];

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        for(MInt lvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1; lvl >= 0; lvl--)
          dists(lvl, i) = (MInt)(fists(lvl, i) / m_lengthOnLevel[bp->maxUniformRefinementLevel + lvl] + 0.5);
    } break;
    case 2: {
      // notes: - localBndRfnDistance : distance in stl units for first layer to be refined
      MFloatScratchSpace distsStlUnits(bp->noLocalBndRfnLvls, bp->noLocalRfnBoundaryIds, AT_, "distsStlUnits");
      distsStlUnits.fill(F0);

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++)
        distsStlUnits(bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1, i) = bp->localBndRfnDistance[i];

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++) {
        const MInt localRfnMaxLvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1;
        const MInt localRfnMinLvl = bp->noLocalBndRfnLvls - bp->localBndRfnMinLvlDiff[i] - 1;
        for(MInt lvl = localRfnMaxLvl - 1; lvl >= localRfnMinLvl; lvl--) {
          distsStlUnits(lvl, i) = distsStlUnits(lvl + 1, i)
                                  + bp->smoothDistance[i] * m_lengthOnLevel[bp->maxUniformRefinementLevel + lvl + 1];
        }
        // ensure a default SmoothDistance to always have a valid grid, even in case of non-meaningful user input
        constexpr MInt defaultSmoothDistance = 2;
        for(MInt lvl = localRfnMinLvl - 1; lvl > -1; lvl--) {
          distsStlUnits(lvl, i) = distsStlUnits(lvl + 1, i)
                                  + defaultSmoothDistance * m_lengthOnLevel[bp->maxUniformRefinementLevel + lvl + 1];
        }
      }

      for(MInt i = 0; i < bp->noLocalRfnBoundaryIds; i++) {
        const MInt localRfnMaxLvl = bp->noLocalBndRfnLvls - bp->localRfnLvlDiff[i] - 1;
        for(MInt lvl = localRfnMaxLvl; lvl >= 0; lvl--)
          dists(lvl, i) = (MInt)(distsStlUnits(lvl, i) / m_lengthOnLevel[bp->maxUniformRefinementLevel + lvl] + 0.5);
      }
    } break;
    default:
      mTerm(1, "no valid localBndRfnMethod");
  }

  // group the Ids with the same distance
  set<MInt> rfnBoundaryGroupDist;
  for(MInt bId = 0; bId < bp->noLocalRfnBoundaryIds; bId++)
    rfnBoundaryGroupDist.insert(dists(refLevel, bId));

  vector<vector<MInt>> rfnBoundaryGroupMemberIds(rfnBoundaryGroupDist.size(), vector<MInt>(0));

  for(MInt bId = 0; bId < bp->noLocalRfnBoundaryIds; bId++) {
    set<MInt>::iterator setIt = rfnBoundaryGroupDist.find(dists(refLevel, bId));
    rfnBoundaryGroupMemberIds[distance(rfnBoundaryGroupDist.begin(), setIt)].push_back(bp->localRfnBoundaryIds[bId]);
  }

  m_log << "       * refinement groups (distance: Ids)" << endl;

  for(MInt gId = 0; gId < (MInt)rfnBoundaryGroupDist.size(); gId++) {
    set<MInt>::iterator setIt = rfnBoundaryGroupDist.begin();
    advance(setIt, gId);
    if(*setIt == 0) continue;

    m_log << "         - group " << *setIt << " :";
    for(MInt mId = 0; mId < (MInt)rfnBoundaryGroupMemberIds[gId].size(); mId++) {
      m_log << " " << rfnBoundaryGroupMemberIds[gId][mId];
    }
    m_log << endl;

    // 6.2b.3.1 do the boundary propagation
    // * we propagate untill we reach cell (*setIt - 1). So this will be the farest marked distance
    propagateDistance(level_, *setIt, rfnBoundaryGroupMemberIds[gId], solver);

    // 6.2b.3.2 mark the according cells
    // * we refine cells with a distance < *setIt, therefore all marked cells
    markBndDistance(level_, solver);
  }
}

template <MInt nDim>
void GridgenPar<nDim>::markSolverForRefinement(MInt level_, MInt solver) {
  TRACE();

  for(MInt k = m_levelOffsets[level_][0]; k < m_levelOffsets[level_][1]; k++) {
    if(a_isInSolver(k, solver)) {
      a_isToRefineForSolver(k, solver) = 1;
    }
  }

  if(noDomains() > 1) {
    MInt lev_pos = 2 * (level_ - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          a_isToRefineForSolver(k, solver) = 1;
        }
      }
    }
  }
}

template <MInt nDim>
void GridgenPar<nDim>::finalizeGrid() {
  TRACE();

  RECORD_TIMER_START(m_t_finalizeGrid);

  // 6.3 Check the refinement validity for LB.
  // Note: can be skipped via property for non-LB grids since this may take a really long time for large grids
  if(m_checkGridLbValidity) {
    checkLBRefinementValidity();
  } else {
    outStream << "     + NOTE: Skipping mesh validity check for LB" << endl;
    m_log << "     + NOTE: Skipping mesh validity check for LB" << endl;
  }

  // 6.4 we are finished, run the final setup, TODO labels:GRIDGEN,totest check the level, what is it needed for
  if(noDomains() > 1) {
    if(m_hasBeenLoadBalanced) {
      communicateHaloGlobalIds(m_maxRfnmntLvl);
    } else {
      updateInterRankNeighbors();
    }

    m_noTotalCells = 0;
    for(MInt d = 0; d < noDomains(); d++) {
      m_noTotalCells += (MLong)m_noCellsPerDomain[d];
    }
  } else {
    m_noTotalCells = (MLong)m_noCells;
    reorderGlobalIdsDF();
    updateGlobalIdsReferences();
  }

  RECORD_TIMER_STOP(m_t_finalizeGrid);
}

/** \brief refines a single cell
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * Refines a single cell by the following algorithm:
 * For 8 child cells do:
 *   a. Create coordinates
 *   b. Init child
 *   c. Check if parent has a cut with the geometry,
 *      if so check this cell for a cut as well
 *      (this calls checkCellForCut(...))
 *   d. Update the parent
 *
 * \param[in] id the id of the cell to refine
 * \param[in] currentChildId the id, where the first child will be written to
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::refineCell(MInt id, MInt* currentChildId) {
  TRACE();

  MInt level_ = a_level(id);
  const MInt childLevel = level_ + 1;
  const MFloat childCellLength = m_lengthOnLevel[childLevel];

  static const MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                           {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};

  a_noChildren(id) = 0;

  for(MInt c = 0; c < m_maxNoChildren; c++) {
    // a. Create coordinates
    for(MInt i = 0; i < nDim; i++) {
      a_coordinate(*currentChildId, i) = a_coordinate(id, i) + F1B2 * signStencil[c][i] * childCellLength;
    }

    // b. Init child
    a_level(*currentChildId) = childLevel;
    a_parentId(*currentChildId) = (MLong)id;
    a_noChildren(*currentChildId) = 0;
    a_globalId(*currentChildId) = (MLong)(*currentChildId);
    for(MInt s = 0; s < m_noSolvers; s++) {
      a_noSolidLayer(*currentChildId, s) = -1;
    }

    // cell is considered outside upon creation and non-boundary
    for(MInt i = 0; i < 8; i++)
      a_hasProperty(*currentChildId, i) = 0;

    for(MInt solver = 0; solver < m_noSolvers; solver++)
      a_isSolverBoundary(*currentChildId, solver) = 0;

    for(MInt solver = 0; solver < m_noSolvers; solver++)
      a_isToRefineForSolver(*currentChildId, solver) = 0;

    // The affiliation is taken from the parent,
    // because there are solver-outside cells that are refined
    // at this stage a_isInSolver = 1 is more a "isMaybeInSolver"
    // a_isInSolver = 0 is isNotInSolver
    for(MInt solver = 0; solver < m_noSolvers; solver++)
      a_isInSolver(*currentChildId, solver) = a_isToRefineForSolver(id, solver) ? 1 : 0;

    // init childs childIds
    for(MInt k = 0; k < m_maxNoChildren; k++)
      a_childId(*currentChildId, k) = -1;

    // init childs neighborIds
    for(MInt i = 0; i < nDim; i++) {
      a_neighborId(*currentChildId, i) = -1;
      a_neighborId(*currentChildId, i + nDim) = -1;
    }

    // c. Check if parent has a cut with the geometry, if so check this cell for a cut as well
    if(a_hasProperty(id, 1)) {
      a_hasProperty(*currentChildId, 1) = checkCellForCut(*currentChildId);
      if(a_hasProperty(*currentChildId, 1)) {
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(!a_isToRefineForSolver(id, solver)) continue;
          for(MInt srfc = 0; srfc < m_noBndIdsPerSolver[solver]; srfc++) {
            if(m_bndCutInfo[solver][srfc]) {
              a_isSolverBoundary(*currentChildId, solver) = 1;
              break;
            }
          }
        }
      }
    }

    // d. Update parent
    a_childId(id, c) = *currentChildId;
    a_noChildren(id) = a_noChildren(id) + 1;


    (*currentChildId)++;
  }
}

/** \brief checks if a cell has a cut with the geometry
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * This does the following:
 *   a. Create a target of the cell for the check
 *   b. Do the intersection test
 *      this calls getIntersectionElements from Geometry
 *   c. Returns the result. If more than one cut has appeared the
 *      cell has a cut
 *
 * \param[in] id the id of the cell to check
 *
 * \return MBool value if the cell has a cut or not
 **/
template <MInt nDim>
MBool GridgenPar<nDim>::checkCellForCut(MInt id) {
  TRACE();
  return m_geometry->getCellIntersectingSurfaces(&a_coordinate(id, 0), m_lengthOnLevel[a_level(id) + 1], m_bndCutInfo);
}

/** \brief updates the children of a given level
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * This algorithm runs over all children of a given level and updates the neighborhood. This
 * is simply be done by checking if a missing connection is considered to be an inside or an
 * outside connection:
 *
 *  - inside connection: the neighbor is one of my siblings
 *  - outside connection: the neighbor is one of the sibling of a neighboring parent
 *
 * The reconstruction is based on hardcoded neighboring connections, i.e., it is known for a given direction
 * which internal sibling is a neighbor due to the relative position of a child beneath a parent.
 * This is also given across neighboring parents. Therefore an outside neighbor can be found by going to the
 * parent and then to the neighboring parent in the desired direction and then inspecting a certain child of
 * this cell.
 *
 * \param[in] level the level of the parent, for which the neighborhood of the children need to be updated
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::findChildLevelNeighbors(MInt** offsets, MInt level_) {
  TRACE();

  MInt t_refineGrid = 0;
  if(level_ < m_minLevel) {
    NEW_SUB_TIMER_STATIC(t_rfnGrid, "find new neighborhood", m_t_createInitialGrid);
    t_refineGrid = t_rfnGrid;
  } else if(level_ < m_maxUniformRefinementLevel) {
    NEW_SUB_TIMER_STATIC(t_rfnGrid, "find new neighborhood", m_t_createStartGrid);
    t_refineGrid = t_rfnGrid;
  } else if(level_ < m_maxRfnmntLvl) {
    NEW_SUB_TIMER_STATIC(t_rfnGrid, "find new neighborhood", m_t_createComputationalGrid);
    t_refineGrid = t_rfnGrid;
  }

  RECORD_TIMER_START(t_refineGrid);

  // loop over all parents
  for(MInt i = offsets[level_][0]; i < offsets[level_][1]; i++) {
    MInt parent = i;
    MLong* childIds = &a_childId(parent, 0);

    // loop over all children
    for(MInt c = 0; c < m_maxNoChildren; c++) {
      if(childIds[c] < 0) continue;

      MInt child = (MInt)childIds[c];
      // loop over all directions of a child
      for(MInt dir = 0; dir < m_noNeighbors; dir++) {
        if(a_neighborId(child, dir) == -1) {
          // this is an inner connection (nghInside3D is a 3D search array but should also work for 2D)
          if(nghInside3D[c][dir] >= 0) {
            a_neighborId(child, dir) = childIds[nghInside3D[c][dir]];
          } else { // this is an outside connection
            // 2D/3D Switch for nghAcrossCell
            const MInt chdir = (nDim == 3) ? nghAcrossCell3D[c][dir] : nghAcrossCell2D[c][dir];

            if(a_neighborId(parent, dir) >= 0 && a_childId((MInt)a_neighborId(parent, dir), chdir) >= 0)
              a_neighborId(child, dir) = a_childId((MInt)a_neighborId(parent, dir), chdir);
          }
        }
      }
    }
  }

  RECORD_TIMER_STOP(t_refineGrid);
}

/** \brief marks inside and outside cells
 *
 * \author Andreas Lintermann
 * \date 29.10.2013
 *
 * This algorithm runs over all cells in the provided offset range and checks for each
 * cell if it has beend marked (b_properties[6]). If not, this cell is checked if is a
 * boundary cell (b_properties[1]). If not, this cell is checked for its position (inside or
 * outside). If it is inside a inside-flooding is performed by calling floodCells
 * with the argument "1" for inside determination. Otherwise an outside-flooding
 * is performed, which call floodCells with the argument "0" for outside
 * determination.
 *
 * \param[in] offsets a pointer to the offsets array, i.e., the normal cells or the
 *                    halo cells
 * \param[in] level_ the level to check
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markInsideOutside(MInt** offsets, MInt level_) {
  TRACE();

  // reset "visited" property, added since markInsideOutsideSolver also uses it
  for(MInt i = offsets[level_][1] - 1; i >= offsets[level_][0]; i--)
    a_hasProperty(i, 6) = 0;

  for(MInt i = offsets[level_][1] - 1; i >= offsets[level_][0]; i--) {
    // this cell was already visited, skip
    if(a_hasProperty(i, 6)) continue;

    // this is a boundary cell, skip
    if(a_hasProperty(i, 1)) {
      a_hasProperty(i, 0) = 1;
      a_hasProperty(i, 6) = 1;
      continue;
    }

    // cell is inside, do inside-flooding
    if(pointIsInside(&a_coordinate(i, 0)) && !(a_hasProperty(i, 1))) {
      stack<MInt> fillStack;
      a_hasProperty(i, 0) = 1;
      a_hasProperty(i, 6) = 1;
      fillStack.push(i);
      floodCells(&fillStack, 1);
    } // cell is outside, do outside-flooding
    else {
      stack<MInt> fillStack;
      a_hasProperty(i, 0) = 0;
      a_hasProperty(i, 6) = 1;
      fillStack.push(i);
      floodCells(&fillStack, 0);
    }
  }
}

/** \brief floods cells
 *
 * \author Andreas Lintermann
 * \date 22.10.2013
 *
 * Having found a cell do a flooding on the inside cells by checking the neighborhood
 * recursively by making use of the stack. Only add cells that are not boundary cells.
 *
 * \param[in] fillStack the stack carrying the initial cells
 * \param[in] marker marks if the cells are marked inside (1) or outside (0)
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::floodCells(stack<MInt>* fillStack, MChar marker) {
  TRACE();

  while(!fillStack->empty()) {
    MInt currentId = fillStack->top();
    fillStack->pop();

    for(MInt n = 0; n < m_noNeighbors; n++)
      if(a_neighborId(currentId, n) >= 0) {
        MInt nghbr = (MInt)a_neighborId(currentId, n);

        // cell is not visited and is not a boundary cell
        if(!(a_hasProperty(nghbr, 6))) {
          if(!(a_hasProperty(nghbr, 1))) {
            fillStack->push(nghbr);
            a_hasProperty(nghbr, 0) = marker;
            a_hasProperty(nghbr, 6) = 1;
          } else {
            a_hasProperty(nghbr, 0) = 1;
            a_hasProperty(nghbr, 6) = 1;
          }
        }
      }
  }
}

template <MInt nDim>
void GridgenPar<nDim>::floodCells(stack<MInt>* fillStack, MChar marker, MInt solver) {
  TRACE();

  while(!fillStack->empty()) {
    MInt currentId = fillStack->top();
    fillStack->pop();

    for(MInt n = 0; n < m_noNeighbors; n++)
      if(a_neighborId(currentId, n) >= 0) {
        MInt nghbr = (MInt)a_neighborId(currentId, n);

        // cell is not visited and is not a boundary cell
        if(!(a_hasProperty(nghbr, 6))) {
          if(!(a_isSolverBoundary(nghbr, solver))) {
            fillStack->push(nghbr);
            a_isInSolver(nghbr, solver) = marker;
            a_hasProperty(nghbr, 6) = 1;
          } else {
            a_isInSolver(nghbr, solver) = 1;
            a_hasProperty(nghbr, 6) = 1;
          }
        }
      }
  }
}

/** \brief deletes the cells outside of the geometry
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * Outside cells are detected in the following way:
 *
 *  a. Mark all inside and outside cells. This calls markInsideOutside(...).
 *  b. perform cutOff
 *  c. marks solver affiliation and declares cells for deletion that are inside but without
 *     solver affiliation
 *  d. keeps some outside cells needed for grid refinement steps along walls in the solver run
 *  e. delete all outside cells. These are the cells which were not marked as inside and
 *     are not a boundary cell. The deletion is performed from the back of the offsetrange.
 *     The deletion moves a cell at the end of the offsetrange to the location of the cell to
 *     be deleted. This is done by calling the copy cell function.
 *
 *         All neighbors are updated to have no reference to the deleted cell anymore.
 *         The size of the offsetrange is decremented by 1 upon deletion.
 *
 * \param[in] level_ the level to check outside cells for
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::deleteOutsideCellsSerial(MInt level_) {
  TRACE();

  MInt t_deleteOutsideCellsSerial = 0;
  if(level_ <= m_minLevel) {
    NEW_SUB_TIMER_STATIC(t_delSer, "delete outside cells serial", m_t_createInitialGrid);
    t_deleteOutsideCellsSerial = t_delSer;
  } else if(level_ <= m_maxUniformRefinementLevel) {
    NEW_SUB_TIMER_STATIC(t_delSer, "delete outside cells serial", m_t_createStartGrid);
    t_deleteOutsideCellsSerial = t_delSer;
  } else if(level_ <= m_maxRfnmntLvl) {
    NEW_SUB_TIMER_STATIC(t_delSer, "delete outside cells serial", m_t_createComputationalGrid);
    t_deleteOutsideCellsSerial = t_delSer;
  }

  RECORD_TIMER_START(t_deleteOutsideCellsSerial);

  MInt num = m_levelOffsets[level_][1] - m_levelOffsets[level_][0];

  m_log << "       * detecting outside cells to delete:" << endl;
  outStream << "       * detecting outside cells to delete:" << endl;

  // a. mark all inside and outside cells
  markInsideOutside(m_levelOffsets, level_);

  // b. perform cut-off
  if(m_cutOff && level_ >= m_minLevel) performCutOff(m_levelOffsets, level_);

  // c. mark solver affiliation and set cells without solver affiliation as outside
  markSolverAffiliation(level_);

  // d keep some of the outside cells
  if(m_keepOutsideBndryCellChildren) {
    keepOutsideBndryCellChildrenSerial(m_levelOffsets[level_], level_);
    // perform another cut-off and mark cell outside the cutOff so that they are surely
    // deleted in the next step!
    if(m_keepOutsideBndryCellChildren == 3 && m_cutOff && level_ >= m_minLevel) {
      performCutOff(m_levelOffsets, level_, true);
    }
  }

  // e. delete all outside cells
  for(MInt i = m_levelOffsets[level_][1] - 1; i >= m_levelOffsets[level_][0]; i--) {
    MInt cell = i;

    // labels:GRID NOTE: hack to keep all refined cells
    // usefull for periodicBoundaries with leveljumps of the different solvers!
    // if(level_ > m_minLevel ) {
    //   a_hasProperty(cell, 0) = true;
    //   a_isInSolver(cell, 0) = true;
    //}

    // this is required, it is an outside  and not an inside or boundary cell
    MBool inSolidLayer = false;
    for(MInt s = 0; s < m_noSolvers; s++) {
      if(a_noSolidLayer(cell, s) >= 0) {
        inSolidLayer = true;
        break;
      }
    }

    if(!(a_hasProperty(cell, 0)) && !(a_hasProperty(cell, 1)) && !inSolidLayer) {
      // update parent
      if(a_parentId(cell) > -1) {
        for(MInt c = 0; c < m_maxNoChildren; c++)
          if(a_childId((MInt)a_parentId(cell), c) == i) {
            a_childId((MInt)a_parentId(cell), c) = -1;
            break;
          }
        a_noChildren((MInt)a_parentId(cell)) = a_noChildren((MInt)a_parentId(cell)) - 1;
      }

      // update all directions
      for(MInt dir = 0; dir < m_noNeighbors; dir++) {
        MInt dirNeighborId = (MInt)a_neighborId(cell, dir);

        // if we have a neighbor
        if(dirNeighborId >= 0) {
          // delete direct reference
          a_neighborId(dirNeighborId, oppositeDirGrid[dir]) = -1;
        }
      }

      // delete the cell from the collector
      if(i != m_levelOffsets[level_][1] - 1) copyCell(m_levelOffsets[level_][1] - 1, i);
      m_levelOffsets[level_][1]--;
    }
  }

  m_noCells = m_levelOffsets[level_][1];

  m_log << "         - outside cells deleted: " << num - (m_levelOffsets[level_][1] - m_levelOffsets[level_][0])
        << endl;
  m_log << "         - new offsets: " << m_levelOffsets[level_][0] << " " << m_levelOffsets[level_][1] << endl;
  outStream << "         - outside cells deleted: " << num - (m_levelOffsets[level_][1] - m_levelOffsets[level_][0])
            << endl;
  outStream << "         - new offsets: " << m_levelOffsets[level_][0] << " " << m_levelOffsets[level_][1] << endl;

  RECORD_TIMER_STOP(t_deleteOutsideCellsSerial);
}

/** \brief performs deleting solid cells lower than maxRefinementLevel
 *
 * \author Jie Ruan, Moritz Waldmann
 * \date 30.09.2018
 *
 * The function call when we set ggp_keepOutsideBndryCellChildren = 3 in property file, for creating
 * solid cells. After reaching maxRefinementLevel when need delete all solid cells which lower than
 * maxRfnmntLvl, thus in this function run through maxRefinementLevel-1 to minlevel to delete all solid
 * cells which havw no child.
 *
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::deleteCoarseSolidCellsSerial() {
  TRACE();

  m_log << "       * detecting coarse solid cells to delete:" << endl;
  cout << "       * detecting coarse solid cells to delete:" << endl;

  MInt noDeletedCellsTotal = 0;

  for(MInt level = m_maxRfnmntLvl - 1; level >= m_minLevel; level--) {
    MInt noDeletedCellsPerLevel = 0;
    cout << "m_levelOffsets[level][0] " << m_levelOffsets[level][0] << " m_levelOffsets[level][1] "
         << m_levelOffsets[level][1] << endl;
    for(MInt cell = m_levelOffsets[level][1] - 1; cell >= m_levelOffsets[level][0]; cell--) {
      if(!(a_hasProperty(cell, 0)) && !(a_hasProperty(cell, 1)) && a_noChildren(cell) == 0) {
        // update parent
        if(a_parentId(cell) > -1) {
          for(MInt c = 0; c < m_maxNoChildren; c++)
            if(a_childId((MInt)a_parentId(cell), c) == cell) {
              a_childId((MInt)a_parentId(cell), c) = -1;
              break;
            }
          a_noChildren((MInt)a_parentId(cell)) = a_noChildren((MInt)a_parentId(cell)) - 1;
        }
        // update all directions
        for(MInt dir = 0; dir < m_noNeighbors; dir++) {
          MInt dirNeighborId = (MInt)a_neighborId(cell, dir);
          // if we have a neighbor
          if(dirNeighborId >= 0) {
            // delete direct reference
            a_neighborId(dirNeighborId, oppositeDirGrid[dir]) = -1;
          }
        }
        // delete the cell from the collector
        if(cell != m_levelOffsets[level][1] - 1) copyCell(m_levelOffsets[level][1] - 1, cell);
        m_levelOffsets[level][1]--;
        noDeletedCellsPerLevel++;
      }
    }
    noDeletedCellsTotal += noDeletedCellsPerLevel;
    for(MInt re = level + 1; re <= m_maxRfnmntLvl; re++) {
      for(MInt i = 0; i < noDeletedCellsPerLevel; i++) {
        copyCell(m_levelOffsets[re][1] - 1 - i, m_levelOffsets[re][0] - 1 - i);
      }
      m_levelOffsets[re][0] = m_levelOffsets[re][0] - noDeletedCellsPerLevel;
      m_levelOffsets[re][1] = m_levelOffsets[re][1] - noDeletedCellsPerLevel;
    }
  }
  m_noCells = m_levelOffsets[m_maxRfnmntLvl][1];

  m_log << "         - coarse solid cells deleted: " << noDeletedCellsTotal << endl;
  cout << "         - coarse solid cells deleted: " << noDeletedCellsTotal << endl;

  for(MInt j = m_minLevel; j <= m_maxRfnmntLvl; j++) {
    m_log << "         Level: " << j << "- new offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1]
          << endl;
    cout << "         Level: " << j << "- new offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1] << endl;
  }
}

/** \brief performs deleting solid cells lower than maxRefinementLevel
 *
 * \author Jie Ruan, Moritz Waldmann
 * \date 30.09.2018
 *
 * The function call when we set ggp_keepOutsideBndryCellChildren = 3 in property file, for creating
 * solid cells. After reaching maxRefinementLevel when need delete all solid cells which lower than
 * maxRfnmntLvl, thus in this function run through maxRefinementLevel-1 to minlevel to delete all solid
 * cells which havw no child.
 *
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::deleteCoarseSolidCellsParallel() {
  TRACE();

  // declare lambda to check for valid neighbors
  auto neighborExists = [&](MInt haloCellId) {
    for(MInt n = 0; n < m_noNeighbors; n++) {
      if(a_neighborId(haloCellId, n) >= 0 && a_neighborId(haloCellId, n) < m_noCells) {
        return true;
      }
    }
    return false;
  };

  for(MInt j = m_minLevel; j <= m_maxRfnmntLvl; j++) {
    m_log << "         Level: " << j << "- old offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1]
          << endl;
    cout << "         Level: " << j << "- old offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1] << endl;
  }

  for(MInt j = m_minLevel; j <= m_maxRfnmntLvl; j++) {
    MInt levPos = 2 * (j - m_minLevel);
    for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
      m_log << "         Level: " << j << "- old offsets: " << m_haloCellOffsets[dom][levPos + 1] << " "
            << m_haloCellOffsets[dom][levPos] << endl;
      cout << "         Level: " << j << "- old offsets: " << m_haloCellOffsets[dom][levPos + 1] << " "
           << m_haloCellOffsets[dom][levPos] << endl;
    }
  }
  m_log << "       * detecting coarse solid cells to delete:" << endl;
  cout << "       * detecting coarse solid cells to delete:" << endl;

  MInt noDeletedCellsTotal = 0;
  MInt noDeletedHalosTotal = 0;
  for(MInt level = m_maxRfnmntLvl - 1; level >= m_minLevel; level--) {
    MInt noDeletedCellsPerLevel = 0;

    for(MInt cell = m_levelOffsets[level][1] - 1; cell >= m_levelOffsets[level][0]; cell--) {
      if(!(a_hasProperty(cell, 0)) && !(a_hasProperty(cell, 1)) && a_noChildren(cell) == 0) {
        // update parent
        if(a_parentId(cell) > -1) {
          for(MInt c = 0; c < m_maxNoChildren; c++)
            if(a_childId((MInt)a_parentId(cell), c) == cell) {
              a_childId((MInt)a_parentId(cell), c) = -1;
              break;
            }
          a_noChildren((MInt)a_parentId(cell)) = a_noChildren((MInt)a_parentId(cell)) - 1;
        }
        // update all directions
        for(MInt dir = 0; dir < m_noNeighbors; dir++) {
          MInt dirNeighborId = (MInt)a_neighborId(cell, dir);
          // if we have a neighbor
          if(dirNeighborId >= 0) {
            // delete direct reference
            a_neighborId(dirNeighborId, oppositeDirGrid[dir]) = -1;
          }
        }
        // delete the cell from the collector
        if(cell != m_levelOffsets[level][1] - 1) copyCell(m_levelOffsets[level][1] - 1, cell);
        m_levelOffsets[level][1]--;
        noDeletedCellsPerLevel++;
      }
    }
    noDeletedCellsTotal += noDeletedCellsPerLevel;
    for(MInt re = level + 1; re <= m_maxRfnmntLvl; re++) {
      for(MInt i = 0; i < noDeletedCellsPerLevel; i++) {
        copyCell(m_levelOffsets[re][1] - 1 - i, m_levelOffsets[re][0] - 1 - i);
      }
      m_levelOffsets[re][0] = m_levelOffsets[re][0] - noDeletedCellsPerLevel;
      m_levelOffsets[re][1] = m_levelOffsets[re][1] - noDeletedCellsPerLevel;
    }

    const MInt levelPos = 2 * (level - m_minLevel);
    MInt noDeletedHalosPerLevel = 0;

    for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
      for(MInt haloCellId = m_haloCellOffsets[dom][levelPos + 1] - 1; haloCellId >= m_haloCellOffsets[dom][levelPos];
          haloCellId--) {
        // if cell hasn't been moved check it
        if(!a_hasProperty(haloCellId, 7)) {
          // 1. for boundary or inside cells check if neighbors still exist
          // 2. if they do skip those cells since they are valid
          if(a_hasProperty(haloCellId, 0) || a_hasProperty(haloCellId, 1) || (a_noSolidLayer(haloCellId, 0) > 0)) {
            if(neighborExists(haloCellId)) {
              continue;
            }
          }

          // delete cell
          cout << " x-dir " << a_coordinate(haloCellId, 0) << " y-dir " << a_coordinate(haloCellId, 1) << " z-dir "
               << a_coordinate(haloCellId, 2) << endl;
          a_hasProperty(haloCellId, 7) = 1;
          deleteCellReferences(haloCellId, haloCellId);
          noDeletedHalosPerLevel++;
        }

        // skip cells marked as moved
        while(m_haloCellOffsets[dom][levelPos] < haloCellId && a_hasProperty(m_haloCellOffsets[dom][levelPos], 7)) {
          m_haloCellOffsets[dom][levelPos]++;
        }

        // copy last cell to current position
        if(haloCellId > m_haloCellOffsets[dom][levelPos]) {
          copyCell(m_haloCellOffsets[dom][levelPos], haloCellId);

          // mark as moved
          a_hasProperty(m_haloCellOffsets[dom][levelPos], 7) = 1;

          // check copied cell
          haloCellId++;
        }

        // increase offset since cell has either been copied or was last cell and invalid
        m_haloCellOffsets[dom][levelPos]++;
      }

      // set new offset
      if(dom > 0) {
        m_haloCellOffsets[dom - 1][levelPos + 1] = m_haloCellOffsets[dom][levelPos];
      }
    }
    noDeletedHalosTotal += noDeletedHalosPerLevel;

    for(MInt re = level + 1; re <= m_maxRfnmntLvl; re++) {
      const MInt levPos = 2 * (re - m_minLevel);
      for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
        for(MInt i = m_haloCellOffsets[dom][levPos + 1] - 1; i >= m_haloCellOffsets[dom][levPos]; i--) {
          copyCell(i, i + noDeletedHalosPerLevel);
        }
        m_haloCellOffsets[dom][levPos] = m_haloCellOffsets[dom][levPos] + noDeletedCellsPerLevel;
        m_haloCellOffsets[dom][levPos + 1] = m_haloCellOffsets[dom][levPos + 1] + noDeletedCellsPerLevel;
      }
    }
  }
  m_noCells = m_levelOffsets[m_maxRfnmntLvl][1];

  m_log << "         - coarse solid cells deleted: " << noDeletedCellsTotal << endl;
  cout << "         - coarse solid cells deleted: " << noDeletedCellsTotal << endl;

  for(MInt j = m_minLevel; j <= m_maxRfnmntLvl; j++) {
    m_log << "         Level: " << j << "- new offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1]
          << endl;
    cout << "         Level: " << j << "- new offsets: " << m_levelOffsets[j][0] << " " << m_levelOffsets[j][1] << endl;
  }
  m_log << "         - coarse solid halos deleted: " << noDeletedHalosTotal << endl;
  cout << "         - coarse solid halos deleted: " << noDeletedHalosTotal << endl;

  for(MInt j = m_minLevel; j <= m_maxRfnmntLvl; j++) {
    MInt levPos = 2 * (j - m_minLevel);
    for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
      m_log << "         Level: " << j << "- new offsets: " << m_haloCellOffsets[dom][levPos + 1] << " "
            << m_haloCellOffsets[dom][levPos] << endl;
      cout << "         Level: " << j << "- new offsets: " << m_haloCellOffsets[dom][levPos + 1] << " "
           << m_haloCellOffsets[dom][levPos] << endl;
    }
  }
}


template <MInt nDim>
void GridgenPar<nDim>::keepOutsideBndryCellChildrenSerial(MInt* offsets, MInt level_) {
  if(m_keepOutsideBndryCellChildren == 1) {
    // check for refinement along a curved wall: if the outside-cell's parent
    // has a boundary-cell-neighbor without children, then the cell is kept
    for(MInt cellId = offsets[0]; cellId < offsets[1]; cellId++) {
      MInt pId = (MInt)a_parentId(cellId);
      for(MInt n = 0; n < m_noNeighbors; n++) {
        MInt nId = (MInt)a_neighborId(pId, n);
        if(nId < 0) continue;
        if(!(a_hasProperty(nId, 1))) continue;
        if(a_noChildren(nId) > 0) continue;
        a_hasProperty(cellId, 0) = 1;
        break;
      }
    }
  } else if(m_keepOutsideBndryCellChildren == 2) {
    // check for refinement along a curved wall: if the outside-cell's parent
    // is a boundary-cell, then the cell is kept
    for(MInt cellId = offsets[0]; cellId < offsets[1]; cellId++) {
      MInt pId = (MInt)a_parentId(cellId);
      if(a_hasProperty(pId, 1)) a_hasProperty(cellId, 0) = 1;
    }
  } else if(m_keepOutsideBndryCellChildren == 3) {
#ifdef MAIA_EXTRA_DEBUG
    MInt SOLID0;
    MInt SOLID1;
    MInt FLUID;
    MInt DELETE;
#endif

    // loop over all solvers:
    for(MInt s = 0; s < m_noSolvers; s++) {
#ifdef MAIA_EXTRA_DEBUG
      FLUID = 0;
      SOLID1 = 0;
      SOLID0 = 0;
#endif
      // initialize all cells;
      // fluid cells have a_noSolidLayer=-1: other cells get a hugh value
      for(MInt cellId = offsets[0]; cellId < offsets[1]; cellId++) {
        if(a_hasProperty(cellId, 0) && !a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = -1;
#ifdef MAIA_EXTRA_DEBUG
          FLUID++;
#endif
        }
        if(a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = m_noSolidLayer[s] * 2;
#ifdef MAIA_EXTRA_DEBUG
          SOLID0++;
#endif
        }
        if(!a_hasProperty(cellId, 0) && !a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = m_noSolidLayer[s] * 2;
#ifdef MAIA_EXTRA_DEBUG
          SOLID1++;
#endif
        }
      }
#ifdef MAIA_EXTRA_DEBUG
      cout << s << " FLUID " << FLUID << " bndry " << SOLID0 << " outside " << SOLID1 << endl;
#endif
      // propagate the solid layers starting from the boundary cells
      for(MInt cellId = offsets[0]; cellId < offsets[1]; cellId++) {
        if(!a_hasProperty(cellId, 1)) continue;
        createSolidCellLayer(cellId, 0, s);
      }
    }

    // loop over all solvers:
    for(MInt s = 0; s < m_noSolvers; s++) {
#ifdef MAIA_EXTRA_DEBUG
      SOLID0 = 0;
      SOLID1 = 0;
      FLUID = 0;
      DELETE = 0;
#endif
      // only outside cells with a_noSolidLayer = -1 will be deleted
      // so reset all cells which
      // have a_noSolidLayer > m_noSolidLayer[s]
      for(MInt cellId = offsets[0]; cellId < offsets[1]; cellId++) {
        if(a_noSolidLayer(cellId, s) > m_noSolidLayer[s]) {
          a_noSolidLayer(cellId, s) = -1;
        }

#ifdef MAIA_EXTRA_DEBUG
        if(a_noSolidLayer(cellId, s) > m_noSolidLayer[s]) DELETE++;
        if(a_noSolidLayer(cellId, s) == -1) FLUID++;
        if(a_noSolidLayer(cellId, s) == 0) SOLID0++;
        if(a_noSolidLayer(cellId, s) == 1) SOLID1++;
#endif
      }
#ifdef MAIA_EXTRA_DEBUG
      cout << s << " FLUID " << FLUID << " bndry " << SOLID0 << " SOLID1 " << SOLID1 << " TODELETE " << DELETE << endl;
#endif
    }
    for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
      for(MInt s = 0; s < m_noSolvers; s++) {
        if(a_noSolidLayer(cellId, s) >= 0) {
          a_isToRefineForSolver(cellId, s) = 1;
          a_isInSolver(cellId, s) = 1;
        }
      }
    }
  } else {
    mTerm(1, AT_, "unknown keepOutsideBndryCellChildren option");
  }
}

template <MInt nDim>
MInt GridgenPar<nDim>::nghborStencil(const MInt neighbor, const MInt dir) {
  static constexpr MInt nghborStencil3d[26][3] = {
      {0, -1, -1}, {1, -1, -1}, {2, -1, -1}, {3, -1, -1}, {4, -1, -1}, {5, -1, -1}, {0, 2, -1}, {0, 3, -1}, {0, 4, -1},
      {0, 5, -1},  {1, 2, -1},  {1, 3, -1},  {1, 4, -1},  {1, 5, -1},  {2, 4, -1},  {2, 5, -1}, {3, 4, -1}, {3, 5, -1},
      {0, 2, 4},   {0, 2, 5},   {0, 3, 4},   {0, 3, 5},   {1, 2, 4},   {1, 2, 5},   {1, 3, 4},  {1, 3, 5}};

  static constexpr MInt nghborStencil2d[26][3] = {
      {0, -1, -1},  {1, -1, -1},  {2, -1, -1},  {3, -1, -1},  {0, 2, -1},   {1, 3, -1},   {2, 1, -1},
      {3, 0, -1},   {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1},
      {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1}};

  IF_CONSTEXPR(nDim == 3) return nghborStencil3d[neighbor][dir];
  return nghborStencil2d[neighbor][dir];
}

template <MInt nDim>
void GridgenPar<nDim>::createSolidCellLayer(MInt cellId, MInt layer, MInt solver) {
  TRACE();

  MInt totalNoNeighbors = (nDim == 3) ? 26 : 8;
  if(a_noSolidLayer(cellId, solver) > layer) {
    a_noSolidLayer(cellId, solver) = layer;
    // Definition: We do not refine at final distance!
    if(layer < m_noSolidLayer[solver]) {
      for(MInt n = 0; n < totalNoNeighbors; n++) {
        //        for (MInt n = 0; n < m_noNeighbors; n++){
        MInt nghborId = cellId;
        for(MInt d = 0; d < nDim; d++) {
          if(nghborStencil(n, d) == -1) break;
          nghborId = (MInt)a_neighborId(nghborId, nghborStencil(n, d));
          if(nghborId == -1) break;
        }
        if(nghborId == -1) continue;
        //        if (a_neighborId(cellId, n) > -1)
        //          createSolidCellLayer((MInt)a_neighborId(cellId, n), layer + 1, solver);
        createSolidCellLayer(nghborId, layer + 1, solver);
      }
    }
  }
}


/** \brief performs a cut-off defined by the properties
 *
 * \author Andreas Lintermann
 * \date 13.12.2013
 *
 * Runs over all cells and checks if cut-offs appear:
 *  - B: a box is used, all inside are kept
 *  - P: a plane is used, all cells that have a positive scalar product with the normal
 *       are kept
 *
 * \param[in] level_ the level to check
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::performCutOff(MInt** offsets, MInt level_, MBool deleteMode) {
  TRACE();

  m_log << "         - performing cut-off on level " << level_ << " with the following options: " << endl;

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    SolverRefinement* bp = m_solverRefinement + solver;
    if(((bp->cutOff == 1) && (level_ == m_minLevel)) || ((bp->cutOff == 2) && (level_ <= m_minLevel))
       || ((bp->cutOff == 3) && (level_ >= m_minLevel)) || (bp->cutOff == 4)) {
      m_log << "         - solver " << solver << endl;

      for(MInt c = 0; c < (signed)bp->cutOffMethods.size(); c++) {
        m_log << "           . type: " << bp->cutOffMethods[c] << endl;
        if(bp->cutOffMethods[c] == "P") {
          m_log << "             : plane point: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : normal: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d + nDim] << " ";
          }
          m_log << endl;
          m_log << "             : layers: ";
          m_log << bp->cutOffNmbrLayers[c][0] << " ";
          m_log << endl;
        } else if(bp->cutOffMethods[c] == "B") {
          m_log << "             : minimal point: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : layers: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffNmbrLayers[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : maximal point: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d + nDim] << " ";
          }
          m_log << endl;
          m_log << "             : layers: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffNmbrLayers[c][d + nDim] << " ";
          }
          m_log << endl;
        } else if(bp->cutOffMethods[c] == "iB") {
          m_log << "             : minimal point: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : layers: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffNmbrLayers[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : maximal point: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffCoordinates[c][d + nDim] << " ";
          }
          m_log << endl;
          m_log << "             : layers: ";
          for(MInt d = 0; d < nDim; d++) {
            m_log << bp->cutOffNmbrLayers[c][d + nDim] << " ";
          }
          m_log << endl;
        } else if(bp->cutOffMethods[c] == "C") {
          m_log << "             : rot axis: ";
          for(MInt d = 0; d < 2; d++) {
            m_log << bp->cutOffCoordinates[c][d] << " ";
          }
          m_log << endl;
          m_log << "             : rot angle: ";
          m_log << bp->cutOffCoordinates[c][2] << " ";
          m_log << endl;
          m_log << "             : rot axis case: ";
          m_log << bp->cutOffCoordinates[c][3] << " ";
          m_log << endl;
          m_log << "             : layers: ";
          m_log << bp->cutOffNmbrLayers[c][0] << " ";
          m_log << endl;
        }
      }

      // check the cut-offs
      for(MInt i = offsets[level_][0]; i < offsets[level_][1]; i++) {
        // skip if the cell is not in the solver
        if(!a_isInSolver(i, solver)) continue;
        // skip if cell is already outside
        // only skip outside cells in the first round!
        if(!deleteMode && a_hasProperty(i, 0) == 0) continue;

        MInt preCut = 0;
        MFloat* coords = &a_coordinate(i, 0);
        if(level_ > m_minLevel) coords = &a_coordinate((MInt)a_parentId(i), 0);
        if(level_ < m_minLevel) preCut = 1;

        MBool keep = true;
        for(MInt c = 0; c < (signed)bp->cutOffMethods.size(); c++) {
          if(bp->cutOffMethods[c] == "P") {
            // scalar product between the difference vector of the plane point and the normal
            MFloat s = 0.0;
            MFloat n = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              s += (coords[d] - bp->cutOffCoordinates[c][d]) * bp->cutOffCoordinates[c][d + nDim];
              n += POW2(bp->cutOffCoordinates[c][d + nDim]);
            }

            if((s / sqrt(n)) < -(bp->cutOffNmbrLayers[c][0] + preCut) * m_lengthOnLevel[level_]) {
              keep = keep && false;
              break;
            }
          } else if(bp->cutOffMethods[c] == "B") {
            for(MInt d = 0; d < nDim; d++)
              if(coords[d]
                     < bp->cutOffCoordinates[c][d] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][d] + preCut)
                 || coords[d] > bp->cutOffCoordinates[c][d + nDim]
                                    + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][d + nDim] + preCut))
                keep = keep && false;
          } else if(bp->cutOffMethods[c] == "iB") {
            if(nDim == 2) {
              if(coords[0]
                     >= bp->cutOffCoordinates[c][0] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][0] + preCut)
                 && coords[0] <= bp->cutOffCoordinates[c][0 + nDim]
                                     + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][0 + nDim] + preCut)
                 && coords[1]
                        >= bp->cutOffCoordinates[c][1] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][1] + preCut)
                 && coords[1] <= bp->cutOffCoordinates[c][1 + nDim]
                                     + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][1 + nDim] + preCut)) {
                keep = keep && false;
              }
            } else if(nDim == 3) {
              if(coords[0]
                     >= bp->cutOffCoordinates[c][0] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][0] + preCut)
                 && coords[0] <= bp->cutOffCoordinates[c][0 + nDim]
                                     + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][0 + nDim] + preCut)
                 && coords[1]
                        >= bp->cutOffCoordinates[c][1] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][1] + preCut)
                 && coords[1] <= bp->cutOffCoordinates[c][1 + nDim]
                                     + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][1 + nDim] + preCut)
                 && coords[2]
                        >= bp->cutOffCoordinates[c][2] - m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][2] + preCut)
                 && coords[2] <= bp->cutOffCoordinates[c][2 + nDim]
                                     + m_lengthOnLevel[level_] * (bp->cutOffNmbrLayers[c][2 + nDim] + preCut)) {
                keep = keep && false;
              }
            }
          } else if(bp->cutOffMethods[c] == "C") {
            MFloat cornerStencil[8][3] = {{1, 1, 1},  {1, -1, 1},  {-1, 1, 1},  {-1, -1, 1},
                                          {1, 1, -1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};

            MFloat cellLength = m_lengthOnLevel[level_ + 1];
            if(level_ > m_minLevel) continue;

            MFloat center[3];
            center[0] = coords[0];
            center[1] = bp->cutOffCoordinates[c][0];
            center[2] = bp->cutOffCoordinates[c][1];
            MFloat dAlpha = (PI / 180.0) * bp->cutOffCoordinates[c][2];

            MFloat phi = (PI / 180.0) * bp->cutOffCoordinates[c][3];

            MFloat normal[3];
            normal[0] = F0;
            normal[1] = sin(phi);
            normal[2] = cos(phi);
            MFloat s = 0.0;
            MFloat n = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              s += (coords[d] - center[d]) * normal[d];
              n += POW2(normal[d]);
            }
            // Cells on opposite site of rot axis
            if(s / sqrt(n) < -(bp->cutOffNmbrLayers[c][0] + preCut) * m_lengthOnLevel[level_]) {
              keep = keep && false;
            }
            if(keep == false) break;


            // Cells outside of cut-offs
            MFloat sgn[2] = {-1.0, 1.0};
            for(MInt p = 0; p < 2; p++) {
              MFloat angle = phi + (sgn[p] * dAlpha / 2.0);
              if(angle < F0) angle = 2 * PI + angle;
              if(angle > 2 * PI) angle = angle - (2 * PI);

              normal[1] = -F1 * sgn[p] * cos(angle);
              normal[2] = sgn[p] * sin(angle);

              MInt cnt = 0;
              for(MInt corn = 0; corn < IPOW2(nDim); corn++) {
                MFloat dummyCoords[3];
                for(MInt d = 0; d < nDim; d++) {
                  dummyCoords[d] = coords[d] + cornerStencil[corn][d] * cellLength;
                }
                s = 0.0;
                n = 0.0;
                for(MInt d = 0; d < nDim; d++) {
                  s += (dummyCoords[d] - center[d]) * normal[d];
                  n += POW2(normal[d]);
                }
                if((s / sqrt(n)) < -(bp->cutOffNmbrLayers[c][0] + preCut) * m_lengthOnLevel[level_]) {
                  cnt++;
                }
              }
              if(cnt >= IPOW2(nDim)) {
                keep = keep && false;
                break;
              }
            }
          }
        }

        if(!keep) {
          a_isInSolver(i, solver) = 0;
          if(deleteMode) {
            a_noSolidLayer(i, solver) = -1;
            // if the last solver: check if removed from all others as well
            if(solver == m_noSolvers - 1) {
              MBool del = true;
              for(MInt s = 0; s < m_noSolvers; s++) {
                if(a_noSolidLayer(i, s) > 0) {
                  del = false;
                  break;
                }
              }
              if(del) {
                a_hasProperty(i, 0) = false;
                a_hasProperty(i, 1) = false;
              }
            }
          }
        }
      }
    }
  }
}

/** \brief deletes the references of a given cell
 *
 * \author Andreas Lintermann
 * \date 06.11.2013
 *
 * 08.05.2014: Removed pointer references (Jerry Grimmen)
 *
 * The algorithm does the following:
 *
 *  a. Run over all neighbors of the given cell and delete the reference
 *     to this cell.
 *  b. Go to the parent, and update the according child-id.
 *
 * \param[in] cell id of the cell to delete
 * \param[in] pos the id of the cell
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::deleteCellReferences(MInt cell, MInt pos) {
  TRACE();

  // a. update all directions
  for(MInt dir = 0; dir < m_noNeighbors; dir++) {
    MInt dirNeighborId = (MInt)a_neighborId(cell, dir);

    // if we have a neighbor
    if(dirNeighborId >= 0) a_neighborId(dirNeighborId, oppositeDirGrid[dir]) = -1;
  }

  // b. update parent
  if(a_parentId(cell) >= 0) {
    MInt parent = (MInt)a_parentId(cell);
    for(MInt c = 0; c < m_maxNoChildren; c++)
      if((MInt)a_childId((MInt)parent, c) == pos) a_childId((MInt)parent, c) = -1;

    a_noChildren((MInt)parent) = a_noChildren((MInt)parent) - 1;
  }
}

/** \brief deletes outside cells in parallel mode
 *
 * \author Andreas Lintermann
 * \date 21.10.2013
 *
 * Similar to deleteOutsideCellsSerial(...) this function does the following:
 *
 *  a. Mark all inside and outside cells: calls markInsideOutside(...) with
 *     the normal cell and halo-cell offsets.
 *  b. perform cutOff
 *  c. Runs over all cells on the given level and deletes those cells which are
 *     marked as outside (b_properties[0] = 0 && b_properties[1] = 0).
 *  d. Runs over all halo cells on the given level and deletes those halo cells which are
 *     marked as outside (b_properties[0] = 0) or not referenced anymore (non-existing neighborhood).
 *      This algorithm runs over all neighbor domains and then over the according halo-cell offsets.
 *  e. Generates look-up table to allow access to cells in globalId order, if required ( i.e.
 *     important for the exchange in updateInterRankNeighbors()).
 *
 * \param[in] level_ the level to check
 *
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::deleteOutsideCellsParallel(MInt level_) {
  TRACE();

  MInt t_deleteOutsideCellsParallel = 0;
  if(level_ <= m_maxUniformRefinementLevel) {
    NEW_SUB_TIMER_STATIC(t_delPar, "delete outside cells parallel", m_t_createStartGrid);
    t_deleteOutsideCellsParallel = t_delPar;
  } else if(level_ <= m_maxRfnmntLvl) {
    NEW_SUB_TIMER_STATIC(t_delPar, "delete outside cells parallel", m_t_createComputationalGrid);
    t_deleteOutsideCellsParallel = t_delPar;
  }


  RECORD_TIMER_START(t_deleteOutsideCellsParallel);

  m_log << "       * detecting outside cells to delete:" << endl;
  outStream << "       * detecting outside cells to delete:" << endl;

  // a. mark all inside and outside cells
  markInsideOutside(m_levelOffsets, level_);
  markInsideOutside(m_haloCellOffsetsLevel, level_);

  // b. perform cut-off
  if(m_cutOff && level_ >= m_minLevel) {
    performCutOff(m_levelOffsets, level_);
    performCutOff(m_haloCellOffsetsLevel, level_);
  }

  // mark solver affiliation set cells without solver affiliation as outside
  markSolverAffiliation(level_);

  // c. delete all outside cells
  const MInt levelPos = 2 * (level_ - m_minLevel);
  const MInt num = m_levelOffsets[level_][1] - m_levelOffsets[level_][0];
  MInt noDeletedCells = 0;

  // b.2 keep some of the outside cells
  if(m_keepOutsideBndryCellChildren) keepOutsideBndryCellChildrenParallel(level_);

  const MInt numHalo = m_haloCellOffsetsLevel[level_][1] - m_haloCellOffsetsLevel[level_][0];

  for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
    MBool anySolver = false;
    for(MInt s = 0; s < m_noSolvers; s++) {
      if(a_noSolidLayer(cellId, s) > 0) {
        anySolver = true;
        break;
      }
    }

    // is neither an inside nor boundary cell -> delete
    if(!a_hasProperty(cellId, 0) && !a_hasProperty(cellId, 1) && !anySolver) {
      // erase entry if not already erased
      if(!a_hasProperty(cellId, 7)) {
        a_hasProperty(cellId, 7) = 1;
        deleteCellReferences(cellId, cellId);
        noDeletedCells++;
      }

      // skip cells marked as moved
      do {
        m_levelOffsets[level_][1]--;
      } while(m_levelOffsets[level_][1] > cellId && a_hasProperty(m_levelOffsets[level_][1], 7));

      // copy last cell to empty position
      if(cellId < m_levelOffsets[level_][1]) {
        copyCell(m_levelOffsets[level_][1], cellId);

        // mark as moved
        a_hasProperty(m_levelOffsets[level_][1], 7) = 1;

        // copied cell needs to be checked
        cellId--;
      }
    }
  }
  m_noCells -= noDeletedCells;

  m_log << "           =deleted " << noDeletedCells << " cells" << endl;
  outStream << "           =deleted " << noDeletedCells << " cells" << endl;

  // d. delete outside as well as now unreferenced halocells
  m_log << "         - halo cell deletion:" << endl;
  outStream << "         - halo cell deletion:" << endl;
  m_log << "           = outside halo cells per domain deleted [new offsets]:" << endl;
  outStream << "           = outside halo cells per domain deleted [new offsets]:" << endl;

  MInt noDeletedHalos = 0;

  // declare lambda to check for valid neighbors
  auto neighborExists = [&](MInt haloCellId) {
    for(MInt n = 0; n < m_noNeighbors; n++) {
      if(a_neighborId(haloCellId, n) >= 0 && a_neighborId(haloCellId, n) < m_noCells) {
        return true;
      }
    }
    return false;
  };

  for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
    for(MInt haloCellId = m_haloCellOffsets[dom][levelPos + 1] - 1; haloCellId >= m_haloCellOffsets[dom][levelPos];
        haloCellId--) {
      // if cell hasn't been moved check it
      if(!a_hasProperty(haloCellId, 7)) {
        // 1. for boundary or inside cells check if neighbors still exist
        // 2. if they do skip those cells since they are valid
        if(a_hasProperty(haloCellId, 0) || a_hasProperty(haloCellId, 1) || (a_noSolidLayer(haloCellId, 0) > 0)) {
          if(neighborExists(haloCellId)) {
            continue;
          }
        }

        // delete cell
        a_hasProperty(haloCellId, 7) = 1;
        deleteCellReferences(haloCellId, haloCellId);
        noDeletedHalos++;
      }

      // skip cells marked as moved
      while(m_haloCellOffsets[dom][levelPos] < haloCellId && a_hasProperty(m_haloCellOffsets[dom][levelPos], 7)) {
        m_haloCellOffsets[dom][levelPos]++;
      }

      // copy last cell to current position
      if(haloCellId > m_haloCellOffsets[dom][levelPos]) {
        copyCell(m_haloCellOffsets[dom][levelPos], haloCellId);

        // mark as moved
        a_hasProperty(m_haloCellOffsets[dom][levelPos], 7) = 1;

        // check copied cell
        haloCellId++;
      }

      // increase offset since cell has either been copied or was last cell and invalid
      m_haloCellOffsets[dom][levelPos]++;
    }

    // set new offset
    if(dom > 0) {
      m_haloCellOffsets[dom - 1][levelPos + 1] = m_haloCellOffsets[dom][levelPos];
    }

    m_log << "             . " << m_neighborDomains[dom] << ": "
          << " [" << m_haloCellOffsets[dom][levelPos] << " " << m_haloCellOffsets[dom][levelPos + 1] << "]" << endl;
    outStream << "             . " << m_neighborDomains[dom] << ": "
              << " [" << m_haloCellOffsets[dom][levelPos] << " " << m_haloCellOffsets[dom][levelPos + 1] << "]" << endl;
  }


  m_haloCellOffsetsLevel[level_][0] += noDeletedHalos;
  m_noTotalHaloCells -= noDeletedHalos;
  m_noHaloCellsOnLevel[level_] -= noDeletedHalos;

  // d. generate look-up table

  // sort index by globalId
  map<MLong, MInt> lut;
  for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
    lut.insert(make_pair(a_globalId(cellId), cellId));
  }
  for(MInt cellId = m_haloCellOffsets[0][levelPos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][levelPos + 1];
      cellId++) {
    lut.insert(make_pair(a_globalId(cellId), cellId));
  }

  // generate lookup table
  auto lutIt = lut.begin();
  for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
    m_cellIdLUT.insert(make_pair(cellId, lutIt->second));
    lutIt++;
  }
  for(MInt cellId = m_haloCellOffsets[0][levelPos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][levelPos + 1];
      cellId++) {
    m_cellIdLUT.insert(make_pair(cellId, lutIt->second));
    lutIt++;
  }

  m_log << "           = outside and unreferenced halo cells deleted [deleted]: "
        << numHalo - (m_haloCellOffsetsLevel[level_][1] - m_haloCellOffsetsLevel[level_][0]) << " [" << noDeletedHalos
        << "]" << endl;
  m_log << "           = new halo offsets: " << m_haloCellOffsetsLevel[level_][0] << " "
        << m_haloCellOffsetsLevel[level_][1] << endl;
  m_log << "           = new halo cells per domain on this level " << level_ << ":" << endl;

  outStream << "           = outside and unreferenced halo cells deleted [deleted]: "
            << numHalo - (m_haloCellOffsetsLevel[level_][1] - m_haloCellOffsetsLevel[level_][0]) << " ["
            << noDeletedHalos << "]" << endl;
  outStream << "           = new halo offsets: " << m_haloCellOffsetsLevel[level_][0] << " "
            << m_haloCellOffsetsLevel[level_][1] << endl;
  outStream << "           = new halo cells per domain on this level " << level_ << ":" << endl;

  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    m_log << "             . " << m_neighborDomains[dom] << " "
          << m_haloCellOffsets[dom][levelPos + 1] - m_haloCellOffsets[dom][levelPos] << endl;

    outStream << "             . " << m_neighborDomains[dom] << " "
              << m_haloCellOffsets[dom][levelPos + 1] - m_haloCellOffsets[dom][levelPos] << endl;
  }

  m_log << "         - outside cells deleted: " << num - (m_levelOffsets[level_][1] - m_levelOffsets[level_][0])
        << endl;
  m_log << "           = new offsets: " << m_levelOffsets[level_][0] << " " << m_levelOffsets[level_][1] << endl;

  outStream << "         - outside cells deleted: " << num - (m_levelOffsets[level_][1] - m_levelOffsets[level_][0])
            << endl;
  outStream << "           = new offsets: " << m_levelOffsets[level_][0] << " " << m_levelOffsets[level_][1] << endl;

  RECORD_TIMER_STOP(t_deleteOutsideCellsParallel);
}

/** \brief marks cells without solver affiliation as outside cells
 *
 * \author Thomas Schilden
 * \date 03.04.2018
 *
 * In multisolver grids, a refined grid might contain cells that are inside the set
 * of solver geometries but are not affiliated to a solver. Furthermore, these cells
 * might be cut by any geometry, i.e., they are boundary cells. Cells without a solver
 * are marked outside and lose their boundary cell status to be deleted.
 *
 * \param[in] level_ the level to check
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markSolverAffiliation(MInt level_) {
  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    excludeInsideOutside(m_levelOffsets, level_, solver);
    if(noDomains() > 1) {
      excludeInsideOutside(m_haloCellOffsetsLevel, level_, solver);
    }
    markInsideOutside(m_levelOffsets, level_, solver);
    if(noDomains() > 1) {
      markInsideOutside(m_haloCellOffsetsLevel, level_, solver);
    }
  }

  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    MChar inside = 0;
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      if(a_isInSolver(i, solver)) {
        inside = 1;
        break;
      }
    }
    if(!inside) {
      a_hasProperty(i, 0) = 0;
      a_hasProperty(i, 1) = 0; // they also lose their status as boundary cells to be deleted
    }
  }

  if(noDomains() > 1) { // might be irrelevant since m_noNeighborDomains = 0
    MInt lev_pos = 2 * (level_ - m_minLevel);

    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
        MChar inside = 0;
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(a_isInSolver(i, solver)) {
            inside = 1;
            break;
          }
        }
        if(!inside) {
          a_hasProperty(i, 0) = 0;
          a_hasProperty(i, 1) = 0;
        }
      }
    }
  }
}

template <MInt nDim>
void GridgenPar<nDim>::keepOutsideBndryCellChildrenParallel(MInt level_) {
  // parallel part
  if(m_keepOutsideBndryCellChildren == 1) {
    keepOutsideBndryCellChildrenSerial(m_levelOffsets[level_], level_);

    const MInt levelPos = 2 * (level_ - m_minLevel);

    // declare lambda to check for valid neighbors
    auto neighborExists = [&](MInt haloCellId) {
      for(MInt n = 0; n < m_noNeighbors; n++) {
        if(a_neighborId(haloCellId, n) >= 0 && a_neighborId(haloCellId, n) < m_noCells) {
          return true;
        }
      }
      return false;
    };

    // delete unreferenced halo cells
    MInt noDeletedHalos = 0;

    for(MInt dom = m_noNeighborDomains - 1; dom >= 0; dom--) {
      for(MInt haloCellId = m_haloCellOffsets[dom][levelPos + 1] - 1; haloCellId >= m_haloCellOffsets[dom][levelPos];
          haloCellId--) {
        // if cell hasn't been moved check it
        if(!a_hasProperty(haloCellId, 7)) {
          // find dereferenced cells
          if(neighborExists(haloCellId)) {
            continue;
          }

          // delete cell
          a_hasProperty(haloCellId, 7) = 1;
          deleteCellReferences(haloCellId, haloCellId);
          noDeletedHalos++;
        }

        // skip cells marked as moved
        while(m_haloCellOffsets[dom][levelPos] < haloCellId && a_hasProperty(m_haloCellOffsets[dom][levelPos], 7)) {
          m_haloCellOffsets[dom][levelPos]++;
        }

        // copy last cell to current position
        if(haloCellId > m_haloCellOffsets[dom][levelPos]) {
          copyCell(m_haloCellOffsets[dom][levelPos], haloCellId);

          // mark as moved
          a_hasProperty(m_haloCellOffsets[dom][levelPos], 7) = 1;

          // check copied cell
          haloCellId++;
        }

        // increase offset since cell has either been copied or was last cell and invalid
        m_haloCellOffsets[dom][levelPos]++;
      }

      // set new offset
      if(dom > 0) {
        m_haloCellOffsets[dom - 1][levelPos + 1] = m_haloCellOffsets[dom][levelPos];
      }
    }

    // update halo offsets
    m_haloCellOffsetsLevel[level_][0] += noDeletedHalos;
    m_noTotalHaloCells -= noDeletedHalos;
    m_noHaloCellsOnLevel[level_] -= noDeletedHalos;

    // create a small LUT for the halos
    map<MInt, MInt> lut;
    map<MInt, MInt> LUT;
    for(MInt cellId = m_haloCellOffsets[0][levelPos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][levelPos + 1];
        cellId++) {
      lut.insert(make_pair(static_cast<MInt>(a_globalId(cellId)), cellId));
    }
    auto lutIt = lut.begin();
    for(MInt cellId = m_haloCellOffsets[0][levelPos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][levelPos + 1];
        cellId++) {
      LUT.insert(make_pair(cellId, lutIt->second));
      lutIt++;
    }

    // identify window cells: we dont use the m_ Look Up Table, so we duplicate some code
    vector<vector<MInt>> winCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));
    vector<vector<MInt>> haloCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));

    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
        a_hasProperty(i, 3) = 0;
      }

      for(MInt p = m_haloCellOffsets[dom][levelPos]; p < m_haloCellOffsets[dom][levelPos + 1]; p++) {
        const MInt cellId = LUT[p];
        haloCellIdsPerDomain[dom].push_back(cellId);
        for(MInt n = 0; n < m_noNeighbors; n++) {
          if(a_neighborId(cellId, n) < m_noCells && a_neighborId(cellId, n) > -1) {
            a_hasProperty((MInt)a_neighborId(cellId, n), 3) = 1;
          }
        }
      }
      for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
        if(a_hasProperty(i, 3)) {
          winCellIdsPerDomain[dom].push_back(i);
        }
      }
    }

    // prepare communication
    MIntScratchSpace noSendWindowPerDomain(m_noNeighborDomains, AT_, "noSendWindowPerDomain");
    MIntScratchSpace noReceiveHaloPerDomain(m_noNeighborDomains, AT_, "noReceiveHaloPerDomain");
    MInt allSend = 0;
    MInt allReceive = 0;

    for(MInt d = 0; d < m_noNeighborDomains; ++d) {
      noReceiveHaloPerDomain[d] = (MInt)haloCellIdsPerDomain[d].size();
      allReceive += noReceiveHaloPerDomain[d];
      noSendWindowPerDomain[d] = (MInt)winCellIdsPerDomain[d].size();
      allSend += noSendWindowPerDomain[d];
    }

    MIntScratchSpace myWindowInt(allSend, AT_, "myWindowInt");
    MIntScratchSpace myHaloInt(allReceive, AT_, "myHaloInt");

    MInt offset = 0;
    for(MInt d = 0; d < m_noNeighborDomains; ++d) {
      for(MInt c = 0; c < noSendWindowPerDomain[d]; c++)
        myWindowInt(offset + c) = a_hasProperty(winCellIdsPerDomain[d][c], 0);
      offset += noSendWindowPerDomain[d];
    }

    communicateIntToNeighbors(myHaloInt, noReceiveHaloPerDomain, myWindowInt, noSendWindowPerDomain, 1);

    // distribute inside/outside to halo cells
    offset = 0;
    for(MInt d = 0; d < m_noNeighborDomains; d++) {
      for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
        a_hasProperty(haloCellIdsPerDomain[d][c], 0) = (MBool)myHaloInt(offset + c);
      }
      offset += noReceiveHaloPerDomain[d];
    }

    // reset the "moved" property
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++)
      for(MInt p = m_haloCellOffsets[dom][levelPos]; p < m_haloCellOffsets[dom][levelPos + 1]; p++)
        a_hasProperty(p, 7) = 0;
  } else if(m_keepOutsideBndryCellChildren == 2) {
    keepOutsideBndryCellChildrenSerial(m_levelOffsets[level_], level_);
    keepOutsideBndryCellChildrenSerial(m_haloCellOffsets[level_], level_);
  } else if(m_keepOutsideBndryCellChildren == 3) {
    mTerm(1, AT_, "ERROR: keepOutsideBndryCellChildren in mode 3 is not implemented for parallel usage yet!!!");

    MInt SOLID0 = 0;
    MInt SOLID1 = 0;
    MInt FLUID = 0;
    MInt DELETE = 0;

    // loop over all solvers!
    for(MInt s = 0; s < m_noSolvers; s++) {
      // initialize all cells; fluid cells have a_noSolidLayer=-1: other cells get a hugh value
      for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
        // skip cell outside the solver and also cells outside the cutOff!
        if(a_isInSolver(cellId, s) == 0) continue;
        if(a_hasProperty(cellId, 0) && !a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = -1;
          FLUID++;
        }
        if(a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = m_noSolidLayer[s] * 2;
          SOLID0++;
        }
        if(!a_hasProperty(cellId, 0) && !a_hasProperty(cellId, 1)) {
          a_noSolidLayer(cellId, s) = m_noSolidLayer[s] * 2;
          SOLID1++;
        }
      }
      cout << "## FLUID " << FLUID << " SOLID0 " << SOLID0 << " SOLID1 " << SOLID1 << " TODELETE " << DELETE << endl;

      // initialize all halos; fluid cells have a_noSolidLayer=-1: other cells get a hugh value
      const MInt lev_pos = 2 * (level_ - m_minLevel);
      for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
        for(MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
          if(a_hasProperty(i, 0)) a_noSolidLayer(i, s) = -1;
          if(a_hasProperty(i, 1)) a_noSolidLayer(i, s) = m_noSolidLayer[s] * 2;
          if(!a_hasProperty(i, 0) && !a_hasProperty(i, 1)) a_noSolidLayer(i, s) = m_noSolidLayer[s] * 2;
        }
      }

      // propagate the solid layers starting from the boundary cells do this
      // only for normal cells, halos will be updated later
      for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
        if(a_isInSolver(cellId, s) == 0) continue;
        if(!a_hasProperty(cellId, 1)) continue;
        createSolidCellLayer(cellId, 0, m_noSolidLayer[s]);
      }

      //############################################################################
      // PREPARE COMMUNICATION
      //############################################################################

      // create a small LUT for the halos
      map<MInt, MInt> lut;
      map<MInt, MInt> LUT;
      for(MInt cellId = m_haloCellOffsets[0][lev_pos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][lev_pos + 1];
          cellId++) {
        lut.insert(make_pair(static_cast<MInt>(a_globalId(cellId)), cellId));
      }
      auto lutIt = lut.begin();
      for(MInt cellId = m_haloCellOffsets[0][lev_pos]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][lev_pos + 1];
          cellId++) {
        LUT.insert(make_pair(cellId, lutIt->second));
        lutIt++;
      }

      // create lists of window and halo cells
      vector<vector<MInt>> winCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));
      vector<vector<MInt>> haloCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));


      for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          a_hasProperty(i, 3) = 0;
        }

        for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++) {
          const MInt cellId = LUT[p];
          haloCellIdsPerDomain[dom].push_back(cellId);
          for(MInt n = 0; n < m_noNeighbors; n++) {
            if(a_neighborId(cellId, n) < m_noCells && a_neighborId(cellId, n) > -1) {
              a_hasProperty((MInt)a_neighborId(cellId, n), 3) = 1;
            }
          }
        }
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          if(a_hasProperty(i, 3)) {
            winCellIdsPerDomain[dom].push_back(i);
          }
        }
      }

      // d. prepare the communication
      MIntScratchSpace noSendWindowPerDomain(m_noNeighborDomains, AT_, "noSendWindowPerDomain");
      MIntScratchSpace noReceiveHaloPerDomain(m_noNeighborDomains, AT_, "noReceiveHaloPerDomain");

      m_log << "         -  process " << domainId() << " needs to transfer to domain [no. cells to transfer]: ";
      cout << "         - process " << domainId() << " needs to transfer to domain [no. cells to transfer]: ";
      for(MInt d = 0; d < m_noNeighborDomains; d++) {
        noSendWindowPerDomain[d] = (MInt)winCellIdsPerDomain[d].size();
        m_log << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
        cout << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
      }
      m_log << endl;
      cout << endl;

      m_log << "         - process " << domainId() << " needs to receive from domain [no. cells to receive]: ";
      cout << "         - process " << domainId() << " needs to receive from domain [no. cells to receive]: ";
      for(MInt d = 0; d < m_noNeighborDomains; d++) {
        noReceiveHaloPerDomain[d] = (MInt)haloCellIdsPerDomain[d].size();
        m_log << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
        cout << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
      }
      m_log << endl;
      cout << endl;

      MInt allSend = 0;
      MInt allReceive = 0;
      vector<MInt> offsetsSend;
      vector<MInt> offsetsReceive;
      for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
        offsetsSend.push_back(allSend);
        offsetsReceive.push_back(allReceive);
        allSend += noSendWindowPerDomain[dom];
        allReceive += noReceiveHaloPerDomain[dom];
      }

      MIntScratchSpace sndBufWin(allSend, AT_, "sndBufWin");
      MIntScratchSpace rcvBufHalo(allReceive, AT_, "rcvBufHalo");

      MUshort finished = 0;

      // e. determine whether we have changes by comparing a_refinementDistance
      //   on halo and window cells
      //   this has to be repeated until we are finished, exchange information as long
      //   as we have local changes
      MInt rounds = 0;
      while(!finished) {
        finished = 1;

        MPI_Request* mpi_request_;
        mAlloc(mpi_request_, m_noNeighborDomains, "mpi_request_", AT_);
        for(MInt d = 0; d < m_noNeighborDomains; d++) {
          for(MInt c = 0; c < noSendWindowPerDomain[d]; c++) {
            sndBufWin[offsetsSend[d] + c] = a_noSolidLayer(winCellIdsPerDomain[d][c], 0);
          }
          MPI_Issend(&(sndBufWin[offsetsSend[d]]), noSendWindowPerDomain[d], MPI_INT, m_neighborDomains[d], 0,
                     mpiComm(), &mpi_request_[d], AT_, "(sndBufWin[offsetsSend[d]])");
        }

        MPI_Status status_;
        for(MInt d = 0; d < m_noNeighborDomains; d++)
          MPI_Recv(&(rcvBufHalo[offsetsReceive[d]]), noReceiveHaloPerDomain[d], MPI_INT, m_neighborDomains[d], 0,
                   mpiComm(), &status_, AT_, "(rcvBufHalo[offsetsReceive[d]])");

        for(MInt d = 0; d < m_noNeighborDomains; d++)
          MPI_Wait(&mpi_request_[d], &status_, AT_);

        for(MInt d = 0; d < m_noNeighborDomains; d++) {
          for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
            const MInt halo = haloCellIdsPerDomain[d][c];
            if(a_isInSolver(halo, s) == 0) continue;
            if(a_noSolidLayer(halo, 0) == -1) {
              for(MInt n = 0; n < m_noNeighbors; n++) {
                MInt nghborId = (MInt)a_neighborId(halo, n);
                if((nghborId < 0) || (nghborId > m_noCells)) continue;
                if((a_noSolidLayer(nghborId, 0) > (rcvBufHalo[c + offsetsReceive[d]] + 1))
                   && (rcvBufHalo[c + offsetsReceive[d]] > -1)) {
                  if(rcvBufHalo[c + offsetsReceive[d]] < m_noSolidLayer[s]) {
                    createSolidCellLayer(nghborId, (rcvBufHalo[c + offsetsReceive[d]] + 1), m_noSolidLayer[s]);
                    finished = 0;
                  }
                }
              }
            }
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &finished, 1, MPI_UNSIGNED_SHORT, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                      "finished");

        rounds++;
      }
      cout << "         - communication rounds required: " << rounds << endl;
    }

    // loop over all solvers:
    for(MInt s = 0; s < m_noSolvers; s++) {
      SOLID0 = 0;
      SOLID1 = 0;
      FLUID = 0;
      DELETE = 0;
      // only outside cells with a_noSolidLayer = -1 will be deleted so reset all cells which
      // have a_noSolidLayer > m_noSolidLayer
      for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
        if(a_isInSolver(cellId, s) == 0) continue;
        if(a_noSolidLayer(cellId, s) == -1) FLUID++;

        if(a_noSolidLayer(cellId, s) > m_noSolidLayer[s]) {
          a_noSolidLayer(cellId, s) = -1;
          DELETE++;
        }

        if(a_noSolidLayer(cellId, s) == 0) SOLID0++;

        if(a_noSolidLayer(cellId, s) == 1) SOLID1++;
      }

      cout << " FLUID " << FLUID << " SOLID0 " << SOLID0 << " SOLID1 " << SOLID1 << " TODELETE " << DELETE << endl;
      const MInt lev_pos = 2 * (level_ - m_minLevel);
      for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
        for(MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
          if(a_noSolidLayer(i, s) > m_noSolidLayer[s]) a_noSolidLayer(i, s) = -1;
        }
      }
    }

    for(MInt cellId = m_levelOffsets[level_][0]; cellId < m_levelOffsets[level_][1]; cellId++) {
      for(MInt s = 0; s < m_noSolvers; s++) {
        if(a_noSolidLayer(cellId, s) >= 0) {
          a_isToRefineForSolver(cellId, s) = 1;
          a_isInSolver(cellId, s) = 1;
        }
      }
    }

    //    for (MInt dom = 0; dom < m_noNeighborDomains; dom++){
    //      for (MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
    //        if(a_noSolidLayer(i) >= 0){
    //          for(MInt block = 0; block < m_noSolvers; block++){
    //            a_isToRefineForSolver(i, block) = 1;
    //            a_isInSolver(i, block) = 1;
    //          }
    //        }
    //      }
    //    }
  } else {
    mTerm(1, AT_, "unknown keepOutsideBndryCellChildren option");
  }
}


/** \brief reorders the cells after Hilbert id
 *
 * \author Andreas Lintermann, Jerry Grimmen
 * \date 16.10.2013, 15.04.2014 (2D Update)
 *
 * Reorders the cells after the Hilbert curve by the following algorithm
 *
 *   a. Find Hilbert ids of all cells. This changes the local coordinates to
 *      be in the unit cube, then calls Hilbert index method to obtain the
 *      Hilbert id. This is only done on a temporary array.
 *   b. The temporary array containing the Hilbert ids is sorted. Additionally
 *      a second array used for lookup is sorted in the same way.
 *   c. The cells are swaped according to the sorting of the Hilbert curve.
 *   d. The global ids of the cells are updated.
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::reorderCellsHilbert() {
  TRACE();

  NEW_SUB_TIMER(t_reorderCellsHilbert, "reordering cells after Hilbert id", m_t_createInitialGrid);
  RECORD_TIMER_START(t_reorderCellsHilbert);

  m_log << "     + creating Hilbert curve" << endl;
  outStream << "     + creating Hilbert curve" << endl;

  // Store center of gravity, length level-0, and Hilbert level locally
  std::array<MFloat, MAX_SPACE_DIMENSIONS> centerOfGravity;
  copy_n(m_centerOfGravity, nDim, &centerOfGravity[0]);
  MFloat lengthLevel0 = m_lengthOnLevel[0];
  MInt hilbertLevel = m_minLevel;

  // Load center of gravity, length level-0, and Hilbert level from target grid
  // file if such a file was specified
  if(m_targetGridFileName != "") {
    TERMM(1, "deprecated");
    m_log << "       * reading grid information from target grid (" << m_targetGridFileName << ")" << endl;
    outStream << "       * reading grid information from target grid (" << m_targetGridFileName << ")" << endl;

    // Read information from target grid file
    using namespace maia::parallel_io;
    ParallelIo file(m_targetGridFileName, PIO_READ, mpiComm());
    file.getAttribute(&centerOfGravity[0], "centerOfGravity", nDim);
    file.getAttribute(&lengthLevel0, "lengthLevel0");
    file.getAttribute(&hilbertLevel, "minLevel");
  } else if(m_multiSolverMinLevel > -1) {
    hilbertLevel = m_multiSolverMinLevel;

    // Determine center of gravity and length level-0 from the given multisolver bounding box
    lengthLevel0 = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      centerOfGravity[i] = 0.5 * (m_multiSolverBoundingBox[nDim + i] + m_multiSolverBoundingBox[i]);
      // Note: computed similar to m_lengthOnLevel
      const MFloat dist =
          (F1 + F1 / FPOW2(30)) * fabs(m_multiSolverBoundingBox[nDim + i] - m_multiSolverBoundingBox[i]);
      lengthLevel0 = std::max(lengthLevel0, dist);
    }

    // Store for later use
    m_multiSolverLengthLevel0 = lengthLevel0;
    m_multiSolverCenterOfGravity.resize(nDim);
    std::copy_n(&centerOfGravity[0], nDim, &m_multiSolverCenterOfGravity[0]);

    // Perform sanity checks to ensure min-level cells match and the same Hilbert order is obtained
    checkMultiSolverGridExtents(nDim, &m_centerOfGravity[0], m_lengthOnLevel[0], m_minLevel, &centerOfGravity[0],
                                lengthLevel0, m_multiSolverMinLevel);

    stringstream msg;
    msg << "       * using multisolver grid information from property file:" << endl
        << "         - minLevel = " << hilbertLevel << endl
        << "         - boundingBox = ";
    for(MInt i = 0; i < nDim; i++) {
      msg << "[" << m_multiSolverBoundingBox[i] << ", " << m_multiSolverBoundingBox[nDim + i] << "]";
      if(i < nDim - 1) {
        msg << " x ";
      }
    }
    msg << endl << "         - centerOfGravity = [";
    for(MInt i = 0; i < nDim; i++) {
      msg << centerOfGravity[i];
      if(i < nDim - 1) {
        msg << ", ";
      }
    }
    msg << "]" << endl << "         - lengthLevel0 = " << lengthLevel0 << endl;

    m_log << msg.str();
    if(domainId() == 0) {
      outStream << msg.str();
    }
  }

  ScratchSpace<MLong> hilbertIds(m_noCells, AT_, "hilbertIds");
  MIntScratchSpace hilbert_lookup(m_noCells, AT_, "hilbert_lookup");
  MIntScratchSpace posTracker(m_noCells, AT_, "posTracker");
  MIntScratchSpace revPosTracker(m_noCells, AT_, "posTracker");
  // a. find Hilbert ids of all cells
  for(MInt i = 0; i < m_noCells; ++i) {
    MFloat* c = &a_coordinate(i, 0);

    // Normalize to unit cube
    array<MFloat, 3> x;
    x[0] = (c[0] - centerOfGravity[0] + lengthLevel0 * 0.5) / lengthLevel0;
    x[1] = (c[1] - centerOfGravity[1] + lengthLevel0 * 0.5) / lengthLevel0;
    IF_CONSTEXPR(nDim == 3) { x[2] = (c[2] - centerOfGravity[2] + lengthLevel0 * 0.5) / lengthLevel0; }

    IF_CONSTEXPR(nDim == 2) { hilbertIds[i] = maia::grid::hilbert::index<2>(&x[0], static_cast<MLong>(hilbertLevel)); }
    else IF_CONSTEXPR(nDim == 3) {
      hilbertIds[i] = maia::grid::hilbert::index<3>(&x[0], static_cast<MLong>(hilbertLevel));
    }
    else {
      TERMM(1, "Bad number of dimensions: " + to_string(nDim));
    }

    hilbert_lookup[i] = i;
    posTracker[i] = i;
    revPosTracker[i] = i;
  }

  m_log << "       * resorting cells after Hilbert id" << endl;
  outStream << "       * resorting cells after Hilbert id" << endl;

  // b. sort temp. array after Hilbert id
  quickSort(hilbertIds.getPointer(), hilbert_lookup.getPointer(), 0, m_noCells - 1);

  // c. swap cells in collector so that the cells are soreted by the Hilbert Id
  for(MInt i = 0; i < m_noCells; ++i) {
    if(i != hilbert_lookup[i]) {
      if(i != posTracker[hilbert_lookup[i]]) swapCells(i, posTracker[hilbert_lookup[i]]);
      MInt iTmp = revPosTracker[i];
      revPosTracker[posTracker[hilbert_lookup[i]]] = iTmp;
      posTracker[iTmp] = posTracker[hilbert_lookup[i]];
    }
  }

  // d. Reset global ids after reordering
  for(MInt i = 0; i < m_noCells; ++i) {
    a_globalId(i) = i;
  }
  RECORD_TIMER_STOP(t_reorderCellsHilbert);
}

/** \brief sorts a list of integers and updates a second one
 *
 * \author Andreas Lintermann
 * \date 16.10.2013
 *
 * This function uses the well known quick-sort algorithm. In addition a second
 * array can be provided to be used as a lookup table. This function is for example to
 * sort the array of Hilbert ids. Additionally, a second array containing a linear
 * increase of cell ids starting from cell 0 is provided. In case the order of the ids
 * is changed in the array to sort the second array is changed accordingly. This can
 * then be used to do a pre-sorting of the Hilbert ids and the second array is
 * incorporated into the swapping of the memory.
 *
 * \tparam[in] T Data type of array to sort.
 * \param[in] globalIdArray the array to sort
 * \param[in] lookup the lookup table to be swapped according to the changes in the globalIdArray
 * \param[in] startIndex the index of the beginning of the fraction to sort
 * \param[in] endIndex the index of the end of the fraction to sort
 *
 **/
template <MInt nDim>
template <class T>
void GridgenPar<nDim>::quickSort(T* globalIdArray, MInt* lookup, MInt startIndex, MInt endIndex) {
  TRACE();

  T pivot = globalIdArray[startIndex];
  MInt splitPoint = 0;

  if(endIndex > startIndex) {
    MInt leftBoundary = startIndex;
    MInt rightBoundary = endIndex;

    while(leftBoundary < rightBoundary) {
      while(pivot < globalIdArray[rightBoundary] && rightBoundary > leftBoundary)
        rightBoundary--;

      swap(globalIdArray[leftBoundary], globalIdArray[rightBoundary]);
      swap(lookup[leftBoundary], lookup[rightBoundary]);

      while(pivot >= globalIdArray[leftBoundary] && leftBoundary < rightBoundary)
        leftBoundary++;

      swap(globalIdArray[rightBoundary], globalIdArray[leftBoundary]);
      swap(lookup[rightBoundary], lookup[leftBoundary]);
    }

    splitPoint = leftBoundary;

    globalIdArray[splitPoint] = pivot;

    quickSort(globalIdArray, lookup, startIndex, splitPoint - 1);
    quickSort(globalIdArray, lookup, splitPoint + 1, endIndex);
  }
}

/** \brief moves a cell from one location in the collector to another
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * Within this function all the information is copied to the new location.
 *
 * \param[in] from the id of the cell to be moved
 * \param[in] to the id of the cell that should be overwritten
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::copyCell(MInt from, MInt to) {
  TRACE();

  // don't copy cells that have been marked moved (messes with neighbors)
  ASSERT(a_hasProperty(from, 7) != 1, "" << globalDomainId() << ": Invalid cell " << from << " copied!");

  for(MInt i = 0; i < 8; i++) {
    a_hasProperty(to, i) = a_hasProperty(from, i);
    a_isInSolver(to, i) = a_isInSolver(from, i);
    a_isSolverBoundary(to, i) = a_isSolverBoundary(from, i);
    a_isToRefineForSolver(to, i) = a_isToRefineForSolver(from, i);
  }

  a_level(to) = a_level(from);
  a_noChildren(to) = a_noChildren(from);
  a_globalId(to) = a_globalId(from);
  a_parentId(to) = a_parentId(from);
  for(MInt s = 0; s < 8; s++) {
    a_noSolidLayer(to, s) = a_noSolidLayer(from, s);
  }

  copy(&a_coordinate(from, 0), &a_coordinate(from, nDim), &a_coordinate(to, 0));
  copy(&a_neighborId(from, 0), &a_neighborId(from, m_noNeighbors), &a_neighborId(to, 0));
  copy(&a_childId(from, 0), &a_childId(from, m_maxNoChildren), &a_childId(to, 0));

  // update the neighbors
  for(MInt i = 0; i < m_noNeighbors; i++)
    if(a_neighborId(to, i) >= 0) a_neighborId((MInt)a_neighborId(to, i), oppositeDirGrid[i]) = (MLong)to;

  // update the parent
  if(a_parentId(to) >= 0)
    for(MInt c = 0; c < m_maxNoChildren; c++)
      if(a_childId((MInt)a_parentId(to), c) == (MLong)from) {
        a_childId((MInt)a_parentId(to), c) = (MLong)to;
        break;
      }

  // update the children
  for(MInt c = 0; c < m_maxNoChildren; c++)
    if(a_childId(to, c) >= 0) a_parentId((MInt)a_childId(to, c)) = (MLong)to;
}

/** \brief swaps two cells in memory
 *
 * \author Andreas Lintermann
 * \date 16.10.2013
 *
 * The swapping is performed as follows:
 *
 * Copy cell1 to the temp, copy cell2 to cell1, copy temp to cell2.
 *
 * \param[in] cellId1 the id of the first cell
 * \param[in] cellId2 the id of the second cell
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::swapCells(MInt cellId1, MInt cellId2) {
  TRACE();

  copyCell(cellId1, m_maxNoCells - 1);
  copyCell(cellId2, cellId1);
  copyCell(m_maxNoCells - 1, cellId2);
}

/** \brief checks if a given point is inside the geometry
 *
 * \author Andreas Lintermann
 * \date 11.10.2013
 *
 * For each direction a ray is created from the given point towards the
 * outer geometry extent. This ray is then used for a cut-test using
 * the getLineIntersectionElements function of Geometry. Two cases can
 * appear:
 *
 *  - the number of cuts is even for the rays: the point is outside
 *  - the number of cuts is odd for the rays: the point is inside
 *
 * \param[in] coordinates the coordinates to check
 *
 * \return if the point is inside or not
 *
 **/
template <MInt nDim>
MBool GridgenPar<nDim>::pointIsInside(MFloat* coordinates) {
  TRACE();
  return m_geometry->isPointInside(coordinates);
}

template <MInt nDim>
MBool GridgenPar<nDim>::pointIsInsideSolver(MFloat* coordinates, MInt solver) {
  TRACE();
  return m_geometry->isPointInsideNode(coordinates, solver);
}

/** \brief creates the start grid
 *
 * \author Andreas Lintermann
 * \date 23.10.2013
 *
 * The algorithm does the following:
 *
 *   5.1. Check if we are running in parallel, if so parallelize the grid.
 *        This calls parallelizeGrid().
 *   5.2. Run over all levels that are still to be refined and do the following:
 *       5.2.1 Update the m_levelOffsets of the next level. This generates a new
 *             offset range for the new level to create under the assumption
 *             that all cells on the current level are refined.
 *       5.2.2 Update the halo cell offsets. This call the updateHaloCellOffsets(...)
 *             function.
 *       5.2.3 Check memory availability.
 *       5.2.4 Refine the normal cells and the halo cells. This calls refineGrid(...)
 *             similar to the serial case. To refine the halo cells the new halo cell
 *             offsets are provided.
 *       5.2.5 Update the number of cells on this domain, since it has changed.
 *       5.2.6 Find the neighbors for all newly generated normal and halo cells.
 *             This calls findChildLevelNeighbors(...) and provides the offsets for
 *             the normal cells and the halo cells of the current level.
 *       5.2.7 Delete all normal outside cells. This has still to be done for the
 *             halo cells.
 *    5.3 Again update the number of cells.
 *    5.4 If the code runs parallel and no further refinement is requested, finish
 *        the process otherwise nothing is required at this stage. Finishing calls
 *        updateInterRankNeighbors().
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::createStartGrid() {
  TRACE();

  RECORD_TIMER_START(m_t_createStartGrid);

  // 5.2. run over all levels that are still to be refined uniformly
  for(MInt l = m_minLevel; l < m_maxUniformRefinementLevel; l++) {
    // 5.2.1 update the m_levelOffsets of the next level
    m_levelOffsets[l + 1][0] = m_noCells;
    m_levelOffsets[l + 1][1] = m_noCells + (m_levelOffsets[l][1] - m_levelOffsets[l][0]) * m_maxNoChildren;

    // 5.2.2 update the halo offset information and size
    if(noDomains() > 1) updateHaloOffsets(l, m_noHaloCellsOnLevel[l], nullptr);

    for(MInt solver = 0; solver < m_noSolvers; solver++)
      markSolverForRefinement(l, solver);

    // 5.2.3 check memory availability
    checkMemoryAvailability(1, l + 1);

    // 5.2.4 refine the normal cells and the halo cells
    refineGrid(m_levelOffsets, l, 0);
    if(noDomains() > 1) {
      refineGrid(m_haloCellOffsetsLevel, l, 1);
    }

    // 5.2.5 update the number of cells on this domain
    m_noCells = m_levelOffsets[l + 1][1];

    // 5.2.6 find the neighbors for all newly generated normal and halo cells
    findChildLevelNeighbors(m_levelOffsets, l);
    if(noDomains() > 1) {
      findChildLevelNeighbors(m_haloCellOffsetsLevel, l);
    }

    // 5.2.7 delete all outside cells
    if(noDomains() > 1) {
      deleteOutsideCellsParallel(l + 1);
    } else {
      deleteOutsideCellsSerial(l + 1);
    }

    // 5.3 update the number of cells
    m_noCells = m_levelOffsets[l + 1][1];
  }

  RECORD_TIMER_STOP(m_t_createStartGrid);
}

/** \brief updates the references of all cells to use the global-ids
 *
 * \author Andreas Lintermann
 * \date 04.12.2013
 *
 * Runs over all cells, checks the referenced neighbors, parents, and children
 * and updates the relation to use the global-id of the according cell.
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::updateGlobalIdsReferences() {
  TRACE();

  for(MInt i = 0; i < m_noCells; i++) {
    for(MInt n = 0; n < m_noNeighbors; n++) {
      if(a_neighborId(i, n) >= 0 && a_neighborId(i, n) < m_noCells) {
        a_neighborId(i, n) = a_globalId((MInt)a_neighborId(i, n));
      }
    }

    for(MInt c = 0; c < m_maxNoChildren; c++)
      if(a_childId(i, c) >= 0) a_childId(i, c) = a_globalId((MInt)a_childId(i, c));

    if(a_parentId(i) >= 0) a_parentId(i) = a_globalId((MInt)a_parentId(i));
  }
}

/** \brief updates the offsets
 *
 * \author Andreas Lintermann
 * \date 02.06.2014
 *
 * \param[in] gridLevel the level for which to update the offsets
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::updateOffsets(MInt gridLevel) {
  TRACE();

  m_levelOffsets[gridLevel + 1][0] = m_noCells;
  m_levelOffsets[gridLevel + 1][1] = m_noCells + m_rfnCount * m_maxNoChildren;

  m_log << "       * cells on level " << gridLevel << " to refine: " << m_rfnCount << endl;
  m_log << "         - new cells to create: " << m_rfnCount * m_maxNoChildren << endl;
  m_log << "         - new no. of cells: " << m_noCells + m_rfnCount * m_maxNoChildren << endl;
  outStream << "       * cells on level " << gridLevel << " to refine: " << m_rfnCount << endl;
  outStream << "         - new cells to create: " << m_rfnCount * m_maxNoChildren << endl;
  outStream << "         - new no. of cells: " << m_noCells + m_rfnCount * m_maxNoChildren << endl;
}

/** \brief updates the offsets of the halos
 *
 * \author Andreas Lintermann
 * \date 12.11.2013
 *
 * Update the halo offset information and size. This generates a new
 * offset range for the new halo cells on the new level under the
 * assumption that all current halo cells are refined. The offset range
 * is created from the back of the collector, without using the very last
 * element (which is reserved as a temporary item for a swapping operation),
 * i.e., the lower the level, the closer we are to the end of the collector.
 * The information is also stored domain-wise.
 *
 * \param[in] l the level to consider
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::updateHaloOffsets(MInt l, MInt noHalos, MInt* rfnCountHalosDom) {
  TRACE();

  m_log << "     + updating the halo offsets for the new level for " << noHalos << " halos" << endl;
  outStream << "     + updating the halo offsets for the new level for " << noHalos << " halos" << endl;

  m_haloCellOffsetsLevel[l + 1][0] = m_haloCellOffsetsLevel[l][0] - noHalos * m_maxNoChildren;
  m_haloCellOffsetsLevel[l + 1][1] = m_haloCellOffsetsLevel[l][0];
  m_noHaloCellsOnLevel[l + 1] = m_haloCellOffsetsLevel[l + 1][1] - m_haloCellOffsetsLevel[l + 1][0];
  m_noTotalHaloCells += m_noHaloCellsOnLevel[l + 1];

  MInt lev_pos = 2 * (l - m_minLevel);
  MInt lev_pos_newlevel = 2 * ((l + 1) - m_minLevel);

  m_log << "       * halo cells on this level for processes [offsets]: " << endl;
  outStream << "       * halo cells on this level for processes [offsets]: " << endl;
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    MInt diff_on_level = 0;
    if(noHalos == m_noHaloCellsOnLevel[l])
      diff_on_level = m_haloCellOffsets[dom][lev_pos + 1] - m_haloCellOffsets[dom][lev_pos];
    else
      diff_on_level = rfnCountHalosDom[dom];

    if(dom == 0)
      m_haloCellOffsets[dom][lev_pos_newlevel] = m_haloCellOffsetsLevel[l + 1][0];
    else
      m_haloCellOffsets[dom][lev_pos_newlevel] = m_haloCellOffsets[dom - 1][lev_pos_newlevel + 1];

    m_haloCellOffsets[dom][lev_pos_newlevel + 1] =
        m_haloCellOffsets[dom][lev_pos_newlevel] + diff_on_level * m_maxNoChildren;
    m_log << "         - " << m_neighborDomains[dom] << ": "
          << " " << m_haloCellOffsets[dom][lev_pos_newlevel + 1] - m_haloCellOffsets[dom][lev_pos_newlevel] << " ["
          << m_haloCellOffsets[dom][lev_pos_newlevel] << " " << m_haloCellOffsets[dom][lev_pos_newlevel + 1] << "]  "
          << endl;
    outStream << "         - " << m_neighborDomains[dom] << ": "
              << " " << m_haloCellOffsets[dom][lev_pos_newlevel + 1] - m_haloCellOffsets[dom][lev_pos_newlevel] << " ["
              << m_haloCellOffsets[dom][lev_pos_newlevel] << " " << m_haloCellOffsets[dom][lev_pos_newlevel + 1]
              << "]  " << endl;
  }
}

/** \brief parallize the present grid
 *
 * \author Andreas Lintermann
 * \date 17.10.2013
 *
 * This function parallelizes the grid with the following algorithm:
 *
 *  a. Generate the rank offsets in the current list of all cells. This
 *     continuously devides the number of the remaining cells by the
 *     number of the remaining ranks.
 *  b. Identify the neighbor domains. This runs over the offset for my
 *     domain and checks if a neighbor belongs to another offset. If so,
 *     a list of neighboring domains is updated by the index of the offset
 *     in which the neighbor is located. The number of halo cells for this
 *     domain is then incremented. The cell in my domain offset is marked
 *     as window cell, while the neighbor on the other domain is marked as
 *     halo cell.
 *  c. Find halo cells in other domains and move them to the end of the
 *     collector. This runs over the offsets of the domains found as a
 *     neighbor in the previous step and moves the cells to the end od the
 *     collector. Therefore the offsets for the halo cells at the end of the
 *     collector is calculated per neighboring domain. The movement also
 *     updates the neighbors which are in my domain and deletes the neighbors
 *     which are outside my domain. The update of local neighbors already
 *     includes the shift performed in the next step.
 *  d. Shift all cells to the beginning of the collector. This completely
 *     moves the solver of my local cells to the beginning of the collector.
 *     Then all neighbor ids are simply updated by shifting the neighbor id
 *     by the size of the shift.
 *  e. Update number of cells and the level offsets. This finally removes
 *     the knowledge about any other cell which does not belong to my domain.
 *
 * \param[in] level_ indicator for the timers
 **/
template <MInt nDim>
void GridgenPar<nDim>::parallelizeGrid() {
  TRACE();

  RECORD_TIMER_START(m_t_parallelizeGrid);

  m_log << "     + parallelizing grid" << endl;
  outStream << "     + parallelizing grid" << endl;

  // a. Generate the rank offsets in the current list of all cells
  m_log << "       * calculating rank offsets" << endl;
  outStream << "       * calculating rank offsets" << endl;
  MIntScratchSpace rank_offsets(noDomains() + 1, AT_, "rank_offsets");

  determineRankOffsets(rank_offsets);

  m_log << "         - rank offsets are : ";
  outStream << "         - rank offsets are : ";
  for(MInt dom = 0; dom < noDomains() + 1; dom++) {
    m_log << rank_offsets[dom] << " ";
    outStream << rank_offsets[dom] << " ";
  }
  m_log << endl;
  outStream << endl;


  // b. Identify the neighbor domains
  m_log << "       * finding communication partners" << endl;
  outStream << "       * finding communication partners" << endl;

  const MInt my_lower_off = rank_offsets[globalDomainId()];
  const MInt my_upper_off = rank_offsets[globalDomainId() + 1];

  // set to save neighbor domains
  set<MInt> neighbor_domains_set;

  // space for saving number of halo cells per domain
  MIntScratchSpace no_haloPerDomain(noDomains(), AT_, "no_haloPerDomain");
  no_haloPerDomain.fill(0);

  // iterate over all assigned cells
  for(MInt i = my_lower_off; i < my_upper_off; i++) {
    for(MInt dir = 0; dir < m_noNeighbors; dir++) {
      const MInt nghbrId = (MInt)a_neighborId(i, dir);

      // check if neighbor exists and is not on own domain
      if(nghbrId >= 0 && (nghbrId < my_lower_off || nghbrId >= my_upper_off)) {
        // mark as window cell
        a_hasProperty(i, 3) = 1;

        // determine to which domain neighbor cell has been assigned
        MInt l = 0;
        while(nghbrId >= rank_offsets[l + 1]) {
          l++;
        }

        neighbor_domains_set.insert(l);

        // mark as halo cell
        if(!a_hasProperty(nghbrId, 4)) {
          no_haloPerDomain[l]++;
          a_hasProperty(nghbrId, 4) = 1;
        }
      }
    }
  }

  m_noNeighborDomains = (MInt)neighbor_domains_set.size();
  mAlloc(m_neighborDomains, m_noNeighborDomains, "m_neighborDomains", -1, AT_);
  mAlloc(m_rfnCountHalosDom, m_noNeighborDomains, "m_rfnCountHalosDom", 0, AT_);

  m_log << "         - domain " << globalDomainId() << " has " << m_noNeighborDomains
        << " neighbor(s) [number of halos]: ";
  outStream << "         - domain " << globalDomainId() << " has " << m_noNeighborDomains
            << " neighbor(s) [number of halos]: ";
  set<MInt>::iterator it = neighbor_domains_set.begin();
  for(MInt j = 0; j < m_noNeighborDomains; j++) {
    m_neighborDomains[j] = *it;
    it++;

    m_log << m_neighborDomains[j] << " [" << no_haloPerDomain[m_neighborDomains[j]] << "]   ";
    outStream << m_neighborDomains[j] << " [" << no_haloPerDomain[m_neighborDomains[j]] << "]   ";
  }
  m_log << endl;
  outStream << endl;


  // c. find halo cells in other domains and move them to the end of the collector
  m_log << "       * finding and moving halo cells that we want to keep" << endl;
  outStream << "       * finding and moving halo cells that we want to keep" << endl;

  // new ordering, remember the offsets for each domains and for each level
  mAlloc(m_haloCellOffsets, m_noNeighborDomains, 2 * (m_maxRfnmntLvl - m_minLevel + 1), "m_haloCellOffsets", 0, AT_);
  const MInt lastNeighborDomPos = m_noNeighborDomains - 1;
  const MInt lastNeighborDomId = m_neighborDomains[lastNeighborDomPos];

  m_haloCellOffsets[lastNeighborDomPos][1] = m_maxNoCells - 1;
  m_haloCellOffsets[lastNeighborDomPos][0] =
      m_haloCellOffsets[lastNeighborDomPos][1] - no_haloPerDomain[lastNeighborDomId];

  for(MInt j = m_noNeighborDomains - 2; j >= 0; j--) {
    m_haloCellOffsets[j][1] = m_haloCellOffsets[j + 1][0];
    m_haloCellOffsets[j][0] = m_haloCellOffsets[j][1] - no_haloPerDomain[m_neighborDomains[j]];
  }

  m_log << "         - haloCellOffsets on this level for processes [offsets]: " << endl;
  outStream << "         - haloCellOffsets on this level for processes [offsets]: " << endl;
  for(MInt j = 0; j < m_noNeighborDomains; j++) {
    m_log << "           = " << m_neighborDomains[j] << ": " << m_haloCellOffsets[j][1] - m_haloCellOffsets[j][0]
          << " [" << m_haloCellOffsets[j][0] << " " << m_haloCellOffsets[j][1] << "]" << endl;
    outStream << "           = " << m_neighborDomains[j] << ": " << m_haloCellOffsets[j][1] - m_haloCellOffsets[j][0]
              << " [" << m_haloCellOffsets[j][0] << " " << m_haloCellOffsets[j][1] << "]" << endl;
  }

  MInt shift = my_lower_off;

  // copy halo cells to the end of the collector
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    MInt domain = m_neighborDomains[dom];
    MInt haloPos = m_haloCellOffsets[dom][0];

    // run over all cells of my neighbor domain
    for(MInt i = rank_offsets[domain]; i < rank_offsets[domain + 1]; i++) {
      // if we found a halo cell that we want to keep
      if(a_hasProperty(i, 4)) {
        copyCell(i, haloPos);
        haloPos++;
      }
    }
  }

  // update the number of halo cells and their offsets
  m_haloCellOffsetsLevel[m_minLevel][0] = m_haloCellOffsets[0][0];
  m_haloCellOffsetsLevel[m_minLevel][1] = m_maxNoCells - 1;

  m_noTotalHaloCells = m_haloCellOffsetsLevel[m_minLevel][1] - m_haloCellOffsetsLevel[m_minLevel][0];
  m_noHaloCellsOnLevel[m_minLevel] = m_noTotalHaloCells;


  // update the halo cell such that it only points to cells that are on my domain
  for(MInt h = m_haloCellOffsetsLevel[m_minLevel][0]; h < m_haloCellOffsetsLevel[m_minLevel][1]; h++) {
    MLong* neigh = &a_neighborId(h, 0);
    for(MInt n = 0; n < m_noNeighbors; n++) {
      if(neigh[n] < my_lower_off || (neigh[n] >= my_upper_off && neigh[n] < m_haloCellOffsetsLevel[m_minLevel][0])) {
        neigh[n] = -1;
      } else if(neigh[n] >= my_lower_off && neigh[n] < my_upper_off) {
        neigh[n] -= shift;
      }
    }
  }

  // d. shift all cells to the beginning of the collector
  MInt noCellsInMyDomain = my_upper_off - my_lower_off;

  m_log << "       * shifting my cell array of size " << noCellsInMyDomain << " by " << my_lower_off << " units"
        << endl;
  outStream << "       * shifting my cell array of size " << noCellsInMyDomain << " by " << my_lower_off << " units"
            << endl;

  m_log << "         - moving memory for pointers in cells" << endl;
  outStream << "         - moving memory for pointers in cells" << endl;

  // perform the cell shift by the shift variable, why not use copyCell(...) to help future us's
  copy(&a_parentId(my_lower_off), &a_parentId(my_lower_off + noCellsInMyDomain), &a_parentId(0));
  copy(&a_globalId(my_lower_off), &a_globalId(my_lower_off + noCellsInMyDomain), &a_globalId(0));
  copy(&a_level(my_lower_off), &a_level(my_lower_off + noCellsInMyDomain), &a_level(0));
  copy(&a_noChildren(my_lower_off), &a_noChildren(my_lower_off + noCellsInMyDomain), &a_noChildren(0));
  copy(&a_neighborId(my_lower_off, 0), &a_neighborId(my_lower_off + noCellsInMyDomain, m_noNeighbors),
       &a_neighborId(0, 0));
  copy(&a_childId(my_lower_off, 0), &a_childId(my_lower_off + noCellsInMyDomain, m_maxNoChildren), &a_childId(0, 0));
  copy(&a_coordinate(my_lower_off, 0), &a_coordinate(my_lower_off + noCellsInMyDomain, nDim), &a_coordinate(0, 0));

  m_log << "         - moving bitsets and updating the neighbors" << endl;
  outStream << "         - moving bitsets and updating the neighbors" << endl;

  // update bitsets and shift the neighbors
  for(MInt i = 0; i < noCellsInMyDomain; i++) {
    MInt cell = i;
    MInt former_cell = i + my_lower_off;
    for(MInt b = 0; b < 8; b++) {
      a_hasProperty(cell, b) = a_hasProperty(former_cell, b);
      a_isInSolver(cell, b) = a_isInSolver(former_cell, b);
      a_isSolverBoundary(cell, b) = a_isSolverBoundary(former_cell, b);
      a_isToRefineForSolver(cell, b) = a_isToRefineForSolver(former_cell, b);
    }
    for(MInt s = 0; s < m_noSolvers; s++) {
      a_noSolidLayer(cell, s) = a_noSolidLayer(former_cell, s);
    }


    MLong* neigh = &a_neighborId(cell, 0);
    for(MInt n = 0; n < m_noNeighbors; n++)
      if(neigh[n] >= 0 && neigh[n] < m_haloCellOffsets[0][0]) neigh[n] -= shift;
  }

  // e. update number of cells and the level offsets
  m_noTotalCells = (MLong)m_noCells;
  m_noCells = noCellsInMyDomain;

  for(MInt l = 0; l < m_minLevel; l++) {
    m_levelOffsets[l][0] = -1;
    m_levelOffsets[l][1] = -1;
  }

  m_levelOffsets[m_minLevel][0] = 0;
  m_levelOffsets[m_minLevel][1] = m_noCells;

  // f. initialize the LUT, needed for updateInterRankNeighbors, (findHaloWindowCells within)
  // add initial halo cells to LUT
  for(MInt cellId = m_haloCellOffsets[0][0]; cellId < m_haloCellOffsets[m_noNeighborDomains - 1][1]; cellId++) {
    m_cellIdLUT.insert(make_pair(cellId, cellId));
  }

  // add intial cells to LUT
  for(MInt cellId = 0; cellId < m_levelOffsets[m_minLevel][1]; cellId++) {
    m_cellIdLUT.insert(make_pair(cellId, cellId));
  }

  RECORD_TIMER_STOP(m_t_parallelizeGrid);
}

template <MInt nDim>
void GridgenPar<nDim>::determineRankOffsets(MIntScratchSpace& rank_offsets) {
  TRACE();

  rank_offsets.fill(0);
  rank_offsets[noDomains()] = m_noCells;
  if(m_weightBndCells || m_weightPatchCells || m_weightSolverUniformLevel) {
    // calculate the cell weights and the sum of the weights
    m_log << "         * using further refinement information" << endl;
    MInt sumOfWeights = 0;
    MIntScratchSpace cellWeights(m_noCells, AT_, "cellWeights");
    cellWeights.fill((MInt)pow(IPOW2(m_maxUniformRefinementLevel - m_minLevel), nDim));
    for(MInt cId = 0; cId < m_noCells; cId++) {
      // solver uniform check
      if(m_weightSolverUniformLevel) {
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          SolverRefinement* bp = m_solverRefinement + solver;
          if(a_isInSolver(cId, solver)) {
            MInt tmpWeight = (MInt)pow(IPOW2(bp->maxUniformRefinementLevel - m_minLevel), nDim);
            cellWeights[cId] = mMax(cellWeights[cId], tmpWeight);
          }
        }
      }
      // Bnd check
      if(m_weightBndCells) {
        // check cell for cut
        if(checkCellForCut(cId)) {
          // get the max refinement level of all cuts
          MInt maxCutLevel = 0;
          for(MInt solver = 0; solver < m_noSolvers; solver++) {
            SolverRefinement* bp = m_solverRefinement + solver;
            if(bp->noLocalBndRfnLvls) { // only check if boundary refinement is enabled for this solver
              for(MInt rfnBndId = 0; rfnBndId < bp->noLocalRfnBoundaryIds; rfnBndId++) {
                if(m_bndCutInfo[solver][bp->localRfnBoundaryIds[rfnBndId]]) {
                  maxCutLevel = mMax(bp->maxBoundaryRfnLvl - bp->localRfnLvlDiff[rfnBndId], maxCutLevel);
                }
              }
            }
          }
          // calculate weight
          MInt cutWeight = (MInt)pow(IPOW2(maxCutLevel - m_minLevel), nDim);
          cellWeights[cId] = mMax(cellWeights(cId), cutWeight);
        }
      }
      // Patch check
      if(m_weightPatchCells) {
        // cell bounding box
        MFloat bbox[2 * nDim];
        for(MInt dim = 0; dim < nDim; dim++) {
          bbox[dim] = a_coordinate(cId, dim) - m_lengthOnLevel[m_minLevel] * 0.5;
          bbox[dim + nDim] = a_coordinate(cId, dim) + m_lengthOnLevel[m_minLevel] * 0.5;
        }
        MInt maxCutLevel = 0;
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(a_isInSolver(cId, solver)) {
            // 1. do the regions overlap?
            SolverRefinement* bp = m_solverRefinement + solver;
            MInt patchLvl = (bp->noLocalPatchRfnLvls() - 1);
            for(; patchLvl >= 0; patchLvl--) { // no patch no enter (patchLvl = -1)
              MBool intersect = false;
              for(MInt patch = 0; patch < (MInt)bp->noPatchesPerLevel(patchLvl); patch++) {
                const MString patchStr = bp->localRfnLevelMethods[patchLvl].substr(patch, 1);
                const MInt pos = bp->localRfnLevelPropertiesOffset[patchLvl] + patch;
                if(patchStr == "B") {
                  intersect = maia::geom::doBoxesOverlap<MFloat, nDim>(bbox, bp->localRfnPatchProperties[pos]);
                } else if(patchStr == "R") {
                  intersect = maia::geom::doBoxAndSphereOverlap<nDim>(bbox, bp->localRfnPatchProperties[pos],
                                                                      bp->localRfnPatchProperties[pos][nDim]);
                } else {
                  m_log << "weighting not implemented for patches other than the BOX and sphere, sorry" << endl;
                  intersect = false;
                }
                if(intersect) break;
              }
              if(intersect) {
                patchLvl += bp->maxUniformRefinementLevel;
                maxCutLevel = mMax(maxCutLevel, patchLvl);
                break;
              }
            }
          }
        }
        MInt patchWeight = (MInt)pow(IPOW2(maxCutLevel - m_minLevel), nDim);
        cellWeights[cId] = mMax(cellWeights(cId), patchWeight);
      }
      // sum of weights for distribution
      sumOfWeights += cellWeights[cId];
    }
    // set offsets
    // average Load
    MInt avgLoad = (MInt)ceil(sumOfWeights / noDomains());
    MInt tmpLoad = 0;
    MInt tmpSumLoad = 0;
    MInt dom = 1;
    for(MInt cId = 0; cId < m_noCells && dom < noDomains(); cId++) {
      if(tmpLoad < avgLoad) {
        // still fits
        tmpLoad += cellWeights[cId];
      } else {
        // too much
        rank_offsets[dom] = cId;
        tmpSumLoad += tmpLoad;
        // update average
        avgLoad = (MInt)ceil((sumOfWeights - tmpSumLoad) / (noDomains() - dom));
        tmpLoad = cellWeights[cId];
        dom++;
      }
    }
  } else {
    rank_offsets[1] = (MInt)ceil(m_noCells / noDomains());
    for(MInt dom = 2; dom < noDomains(); dom++)
      rank_offsets[dom] =
          rank_offsets[dom - 1] + (MInt)ceil((m_noCells - rank_offsets[dom - 1]) / (noDomains() - dom + 1));
  }
}


/** \brief updates the neighbors on the neighboring ranks
 *
 * \author Andreas Lintermann
 * \date 18.10.2013
 *
 * The algorithm does the following:
 *
 *  a. Get the number of cells created on the other domains.
 *     This fills an array m_noCellsPerDomain with the number of
 *     cells created on each domain.
 *
 *  b. Determine my own offset. This calulates the offset from where we
 *     need to write in the I/O-routine. The value is stored in the
 *     variable m_cellOffsetPar.
 *
 *  c. Create lists of window and halo cells. This creates a vector
 *     of vectors of window and halo ids (winCellIdsPerDomain, haloCellIdsPerDomain).
 *     The first dimensions represents the neighbor domain and the second the
 *     according window and halo cell id.
 *
 *  d. Find the halo and window cells for all domains. This runs over all neighbor domains
 *     First, all normal cells are cleared and then remarked as window cell domain-wise.
 *     They are then collected in winCellIdsPerDomain. At the end the global-ids are updated.
 *
 *  e. Reorder the global ids to be sorted depth-first. This calls reorderGlobalIdsDF(), which
 *     reorders only the global-ids. The cells are not moved in moved in memory.
 *
 *  f. Wite parallel geometry if requested.
 *
 *  g. Prepare the communication. This allocates the buffers for the cells to send and
 *     receive.
 *
 *  h. Fill the sendbuffer and send it (domain-wise). The buffer is filled by the
 *     accroding global-id.
 *
 *  i. Receive from neighbors (domain-wise).
 *
 *  j. Update the neighbors, parents and children to use the global-id. Runs over
 *     all cells and checks the global-id of the neighbors, parent, and children
 *     and updates the local information accordingly. This is only done for associated
 *     cells that are not halo cells.
 *
 *  k. Update the global ids of the halos and the acoording window neighbors. The receive
 *     buffer contains all the global-ids from the windows of the neighbors. The halos
 *     are filled by these global-ids. Then the according neighbor, which is a window cell
 *     in the normal cell collector is updated for its neighborhood, to point to the
 *     global-id that we have received.
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::updateInterRankNeighbors() {
  TRACE();

  NEW_SUB_TIMER(t_updateInterRankNeighbors, "update inter rank neighbors", m_t_finalizeGrid);
  m_t_updateInterRankNeighbors = t_updateInterRankNeighbors;
  RECORD_TIMER_START(m_t_updateInterRankNeighbors);

  m_log << "     + updating the inter rank neighbors" << endl;
  outStream << "     + updating the inter rank neighbors" << endl;

  // a. get the number of cells created on the other domains
  MInt* sndBuf = &m_noCells;
  MPI_Allgather(sndBuf, 1, MPI_INT, m_noCellsPerDomain, 1, MPI_INT, mpiComm(), AT_, "sndBuf", "m_noCellsPerDomain");

  // b. determine my own offset
  for(MInt d = 0; d < globalDomainId(); d++)
    m_cellOffsetPar += (MLong)m_noCellsPerDomain[d];

  // c. create lists of window and halo cells
  vector<vector<MInt>> winCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));
  vector<vector<MInt>> haloCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));

  // d. find the halo and window cells for all domains
  if(m_hasBeenLoadBalanced)
    findHaloAndWindowCellsKD(winCellIdsPerDomain, haloCellIdsPerDomain);
  else
    findHaloAndWindowCells(winCellIdsPerDomain, haloCellIdsPerDomain);

  // e. reorder the global ids to be sorted depth-first
  reorderGlobalIdsDF();

  // f. write parallel geometry if requested
  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    if(m_geometry->nodeSurfaceType(solver) == STL) {
      m_geometry->setGeometryPointerToNode(m_STLgeometry, solver);
      if(m_STLgeometry->m_parallelGeometry) writeParallelGeometry();
    }
  }

  // g. prepare the communication
  MIntScratchSpace noSendWindowPerDomain(m_noNeighborDomains, AT_, "noSendWindowPerDomain");
  MIntScratchSpace noReceiveHaloPerDomain(m_noNeighborDomains, AT_, "noReceiveHaloPerDomain");

  m_log << "       * we need to transfer to domain [no. cells to transfer]: ";
  outStream << "       * we need to transfer to domain [no. cells to transfer]: ";
  for(MInt d = 0; d < m_noNeighborDomains; d++) {
    noSendWindowPerDomain[d] = (signed)winCellIdsPerDomain[d].size();
    m_log << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
    outStream << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
  }
  m_log << endl;
  outStream << endl;

  m_log << "       * we need to receive from domain [no. cells to receive]: ";
  outStream << "       * we need to receive from domain [no. cells to receive]: ";
  for(MInt d = 0; d < m_noNeighborDomains; d++) {
    noReceiveHaloPerDomain[d] = (signed)haloCellIdsPerDomain[d].size();
    m_log << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
    outStream << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
  }
  m_log << endl;
  outStream << endl;

  // count all
  MInt allSend = 0;
  MInt allReceive = 0;
  vector<MInt> offsetsSend;
  vector<MInt> offsetsReceive;
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    offsetsSend.push_back(allSend);
    offsetsReceive.push_back(allReceive);
    allSend += noSendWindowPerDomain[dom];
    allReceive += noReceiveHaloPerDomain[dom];
  }

  MLongScratchSpace sndBufWin(allSend, AT_, "sndBufWin");
  MLongScratchSpace rcvBufHalo(allReceive, AT_, "rcvBufHalo");

  // h. fill the sendbuffer and send it
  MPI_Request* mpi_request = nullptr;
  mAlloc(mpi_request, m_noNeighborDomains, "mpi_request", AT_);
  for(MInt d = 0; d < m_noNeighborDomains; d++) {
    for(MInt c = 0; c < noSendWindowPerDomain[d]; c++)
      sndBufWin[offsetsSend[d] + c] = a_globalId(winCellIdsPerDomain[d][c]);

    MPI_Issend(&(sndBufWin[offsetsSend[d]]), noSendWindowPerDomain[d], MPI_LONG, m_neighborDomains[d], 0, mpiComm(),
               &mpi_request[d], AT_, "(sndBufWin[offsetsSend[d]])");
  }

  // i. receive from neighbors
  MPI_Status status;
  for(MInt d = 0; d < m_noNeighborDomains; d++)
    MPI_Recv(&(rcvBufHalo[offsetsReceive[d]]), noReceiveHaloPerDomain[d], MPI_LONG, m_neighborDomains[d], 0, mpiComm(),
             &status, AT_, "(rcvBufHalo[offsetsReceive[d]])");


  for(MInt d = 0; d < m_noNeighborDomains; d++)
    MPI_Wait(&mpi_request[d], &status, AT_);

  // j. update the neighbors, parents and children to use the global id
  updateGlobalIdsReferences();

  // k. update the global ids of the halos and the acoording window neighbors
  for(MInt d = 0; d < m_noNeighborDomains; d++)
    for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
      MInt halo = haloCellIdsPerDomain[d][c];
      a_globalId(halo) = rcvBufHalo[c + offsetsReceive[d]];
      for(MInt n = 0; n < m_noNeighbors; n++) {
        if(a_neighborId(halo, n) >= 0 && a_neighborId(halo, n) < m_noCells) {
          a_neighborId((MInt)a_neighborId(halo, n), oppositeDirGrid[n]) = a_globalId(halo);
        }
      }
    }

  RECORD_TIMER_STOP(m_t_updateInterRankNeighbors);
}

/** \brief recursively traverse the tree depth-first and collect halo cells
 *
 * \author Andreas Lintermann
 * \date 20.10.2013
 *
 * This function calls itself as long as children are avalailable. It traverses
 * the tree starting at given halo node and collects the halos in pre-order
 * depth-first manner in a vector which is passed down as a pointer.
 *
 * \param[in] parentId the id of the current parrent to traverse
 * \param[in] cellIdsPerDomain the collection of the halo cells passed so far
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::collectHaloChildren(MInt parentId_, vector<MInt>* cellIdsPerDomain) {
  TRACE();

  for(MInt c = 0; c < m_maxNoChildren; c++) {
    if(a_childId(parentId_, c) >= 0 && a_hasProperty((MInt)a_childId(parentId_, c), 4)) {
      cellIdsPerDomain->push_back((MInt)a_childId(parentId_, c));
      collectHaloChildren((MInt)a_childId(parentId_, c), cellIdsPerDomain);
    }
  }
}

/** \brief recursively traverse the tree depth-first and collect window cells
 *
 * \author Andreas Lintermann
 * \date 23.10.2013
 *
 * This function calls itself as long as children are avalailable. It traverses
 * the tree starting at given window node and collects the windows in pre-order
 * depth-first manner in a vector which is passed down as a pointer.
 *
 * \param[in] parentId the id of the current parrent to traverse
 * \param[in] cellIdsPerDomain the collection of the window cells passed so far
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::collectWindowChildren(MInt parentId_, vector<MInt>* cellIdsPerDomain) {
  TRACE();

  for(MInt c = 0; c < m_maxNoChildren; c++) {
    if(a_childId(parentId_, c) >= 0 && a_hasProperty((MInt)a_childId(parentId_, c), 3)) {
      cellIdsPerDomain->push_back((MInt)a_childId(parentId_, c));
      collectWindowChildren((MInt)a_childId(parentId_, c), cellIdsPerDomain);
    }
  }
}

// TODO labels:GRIDGEN the grid generator should not know about triangles etc.
template <MInt nDim>
void GridgenPar<nDim>::writeParallelGeometry() {
  TRACE();

  NEW_SUB_TIMER(t_parGeom, "parallel geometry", m_t_updateInterRankNeighbors);
  RECORD_TIMER_START(t_parGeom);

  vector<MInt> cutting_tris;
  MInt noBndCells = 0;
  MIntScratchSpace noTrisPerPartitionCell((MInt)m_partitionCellList.size(), AT_, "noTrisPerPartitionCell");

  MIntScratchSpace dummyInt(1, AT_, "dummyInt");
  MFloatScratchSpace dummyFloat(1, AT_, "dummyFloat");

  // loop over all cells (which are already in Hilbert order)
  for(MInt i = 0; i < (signed)m_partitionCellList.size(); i++) {
    MInt partitionCellId = get<0>(m_partitionCellList[i]);
    noTrisPerPartitionCell[i] = 0;

    // only use the boundary cells
    if(a_hasProperty(partitionCellId, 1)) {
      // a. Create target for check
      MFloat target[6];
      MFloat cellHalfLength = m_lengthOnLevel[a_level(partitionCellId) + 1];

      // add 0.5%
      cellHalfLength += cellHalfLength * 0.005;

      std::vector<MInt> nodeList;

      for(MInt j = 0; j < nDim; j++) {
        target[j] = a_coordinate(partitionCellId, j) - cellHalfLength;
        target[j + nDim] = a_coordinate(partitionCellId, j) + cellHalfLength;
      }

      // b. Do the intersection test
      m_STLgeometry->getIntersectionElements(target, nodeList, cellHalfLength, &a_coordinate(partitionCellId, 0));
      const MInt noNodes = (signed)nodeList.size();

      noTrisPerPartitionCell[i] = noNodes;

      for(MInt t = 0; t < noNodes; t++)
        cutting_tris.push_back(nodeList[t]);

      noBndCells++;
    }
  }

  MInt sumTris = (MInt)cutting_tris.size();

  MIntScratchSpace noTrianglesPerDomain(noDomains(), AT_, "noTrianglesPerDomain");
  MPI_Allgather(&sumTris, 1, MPI_INT, noTrianglesPerDomain.getPointer(), 1, MPI_INT, mpiComm(), AT_, "sumTris",
                "noTrianglesPerDomain.getPointer()");

  MInt noLocalPartitionCells = (MInt)m_partitionCellList.size();
  MIntScratchSpace noPartitionCellsPerDomain(noDomains(), AT_, "noPartitionCellsPerDomain");
  MPI_Allgather(&noLocalPartitionCells, 1, MPI_INT, noPartitionCellsPerDomain.getPointer(), 1, MPI_INT, mpiComm(), AT_,
                "noLocalPartitionCells", "noPartitionCellsPerDomain.getPointer()");

  MInt sumAllTris = 0;
  MInt sumAllPartitionCells = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    sumAllTris += noTrianglesPerDomain[d];
    sumAllPartitionCells += noPartitionCellsPerDomain[d];
  }
  MInt myTriOffset = 0;
  MInt myPartitionCellOffset = 0;
  for(MInt d = 0; d < domainId(); d++) {
    myTriOffset += noTrianglesPerDomain[d];
    myPartitionCellOffset += noPartitionCellsPerDomain[d];
  }

  using namespace maia::parallel_io;
  ParallelIo geomIO(m_STLgeometry->m_parallelGeomFileName, PIO_REPLACE, mpiComm());

  geomIO.defineScalar(PIO_INT, "noTriangles");
  geomIO.defineScalar(PIO_INT, "noRealTriangles");
  geomIO.defineScalar(PIO_FLOAT, "memIncreaseFactor");
  geomIO.defineArray(PIO_INT, "noTrisPerPartitionCell", sumAllPartitionCells);
  geomIO.defineArray(PIO_INT, "originalTriId", sumAllTris);
  geomIO.defineArray(PIO_INT, "segmentId", sumAllTris);

  MString normalNames[3] = {"normals0", "normals1", "normals2"};
  for(MInt d = 0; d < nDim; d++)
    geomIO.defineArray(PIO_FLOAT, normalNames[d], sumAllTris);

  MString vertexNames[3][3] = {{"vertices00", "vertices01", "vertices02"},
                               {"vertices10", "vertices11", "vertices12"},
                               {"vertices20", "vertices21", "vertices22"}};
  for(MInt v = 0; v < nDim; v++)
    for(MInt d = 0; d < nDim; d++)
      geomIO.defineArray(PIO_FLOAT, vertexNames[v][d], sumAllTris);

  geomIO.writeScalar(sumAllTris, "noTriangles");
  geomIO.writeScalar(m_STLgeometry->GetNoElements(), "noRealTriangles");
  geomIO.writeScalar((MFloat)(sumAllTris) / ((MFloat)(m_STLgeometry->GetNoElements())), "memIncreaseFactor");

  geomIO.setOffset(noLocalPartitionCells, myPartitionCellOffset);
  if(noLocalPartitionCells == 0)
    geomIO.writeArray(dummyInt.getPointer(), "noTrisPerPartitionCell");
  else
    geomIO.writeArray(noTrisPerPartitionCell.getPointer(), "noTrisPerPartitionCell");

  geomIO.setOffset(sumTris, myTriOffset);
  if(sumTris == 0)
    geomIO.writeArray(dummyInt.getPointer(), "originalTriId");
  else
    geomIO.writeArray(cutting_tris.data(), "originalTriId");

  {
    MIntScratchSpace segmentIdEntry(sumTris, AT_, "segmentIdEntry");
    for(MInt i = 0; i < sumTris; i++) {
      MInt tri = cutting_tris[i];
      segmentIdEntry[i] = m_STLgeometry->elements[tri].m_segmentId;
    }
    if(sumTris == 0)
      geomIO.writeArray(dummyInt.getPointer(), "segmentId");
    else
      geomIO.writeArray(segmentIdEntry.getPointer(), "segmentId");
  }


  for(MInt d = 0; d < nDim; d++) {
    MFloatScratchSpace normalEntry(sumTris, AT_, "normalEntry");
    for(MInt i = 0; i < sumTris; i++) {
      MInt tri = cutting_tris[i];
      normalEntry[i] = m_STLgeometry->elements[tri].m_normal[d];
    }
    if(sumTris == 0)
      geomIO.writeArray(dummyFloat.getPointer(), normalNames[d]);
    else
      geomIO.writeArray(normalEntry.getPointer(), normalNames[d]);
  }

  for(MInt v = 0; v < nDim; v++)
    for(MInt d = 0; d < nDim; d++) {
      MFloatScratchSpace vertexEntry(sumTris, AT_, "vertexEntry");
      for(MInt i = 0; i < sumTris; i++) {
        MInt tri = cutting_tris[i];
        vertexEntry[i] = m_STLgeometry->elements[tri].m_vertices[v][d];
      }
      if(sumTris == 0)
        geomIO.writeArray(dummyFloat.getPointer(), vertexNames[v][d]);
      else
        geomIO.writeArray(vertexEntry.getPointer(), vertexNames[v][d]);
    }

  RECORD_TIMER_STOP(t_parGeom);
}

/** \brief reorders the globalIds depth-first
 *
 * \author Andreas Lintermann
 * \date 11.12.2013
 *
 * Takes the cells on the coarsest level and traverses them by calling
 * traverseDFGlobalId(...). The code does the following:
 *
 * a. Traverse the cells by running over the cells on the coarsest level and descending
 *    in depth-first order. While doing so update the global ids and store all the ids and the
 *    according offsprings in the 2D scratch partitionCellList.
 * b. Since in the last step all cells were inserted, pick only those, which are really
 *    min-cells, i.e., those cells that have a less or equal number of offspring than
 *    m_partitionCellOffspringThreshold.
 * c. Update the number of min-cells.
 *
 * \param[in] level_ indicator for the timers
 **/
template <MInt nDim>
void GridgenPar<nDim>::reorderGlobalIdsDF() {
  TRACE();

  NEW_SUB_TIMER_STATIC(t_ro, "reorder global ids DF", m_t_finalizeGrid);
  RECORD_TIMER_START(t_ro);

  m_log << "     + reordering global ids" << endl;
  outStream << "     + reordering global ids" << endl;


  m_log << "       * traversing all cells depth-first" << endl;
  outStream << "       * traversing all cells depth-first" << endl;

  MLong currentGlobalId = m_cellOffsetPar - 1;

  // partitionCellList(ids,0) -> localIds
  // partitionCellList(ids,1) -> noOffsprings of cells
  MLongScratchSpace partitionCellList(m_noCells, 2, AT_, "partitionCellList");
  // workload of cells
  MFloatScratchSpace workloadPerCell(m_noCells, 1, AT_, "workloadPerCell");
  // total work load for recursively determine the workload of each cell (=partitionCellList(ids,2))
  MFloatScratchSpace workload(m_noCells, AT_, "workload");
  // set weights if necessary
  MFloatScratchSpace weight(m_noCells, AT_, "weight");

  if(m_weightMethod > 0)
    setCellWeights(weight);
  else
    weight.fill(1);

  // a. traverse
  MInt j = 0;

  // After load balancing, it is possible that some domains have so called partition level shift.
  // This means we have at least one cell which has a halo cell as parent.
  // The real first cell is always the one at m_levelOffsets[m_minLevel + m_noMissingParents][0].
  // The next cell is then either m_levelOffsets[m_minLevel + m_noMissingParents][0] + 1 (The cell 'beside' the first
  // one) Or m_levelOffsets[m_minLevel + m_noMissingParents - 1][0] (The cell 'one level up' the first one) If we got up
  // in level such that we are back on m_minLevel, we continue like always.

  MInt noMissingParents = m_noMissingParents;
  while(noMissingParents > 0) {
    MInt k = 0;
    MInt cellId = m_levelOffsets[m_minLevel + noMissingParents][0];
    while((a_hasProperty((MInt)a_parentId(cellId), 4) == 1) && (a_level(cellId) == (m_minLevel + noMissingParents))) {
      currentGlobalId++;
      a_globalId(cellId) = currentGlobalId;
      partitionCellList(j, 0) = cellId;
      // save the current weight of the initial level cells
      MFloat currentWorkload = weight(cellId);

      // save last pos
      MInt last = j;

      traverseDFGlobalId(cellId, &currentGlobalId, partitionCellList, workloadPerCell, &currentWorkload, weight,
                         workload, &j);

      MLong noOffsprings = 0;
      if(cellId == 0)
        noOffsprings = currentGlobalId - m_cellOffsetPar;
      else
        noOffsprings = currentGlobalId - a_globalId(cellId);

      // update
      partitionCellList(last, 1) = noOffsprings;
      workloadPerCell(last) = currentWorkload;

      ++j;
      ++k;
      cellId = m_levelOffsets[m_minLevel + noMissingParents][0] + k;

    } // end of : while ( (a_hasProperty(a_parentId(cellId), 4) == 1) && (a_level(cellId) == (m_minLevel +
      // noMissingParents)) )

    --noMissingParents;

  } // end of : while ( noMissingParents > 0 )

  for(MInt i = m_levelOffsets[m_minLevel][0]; i < m_levelOffsets[m_minLevel][1]; i++, j++) {
    currentGlobalId++;
    a_globalId(i) = currentGlobalId;
    partitionCellList(j, 0) = i;
    // save the current weight of the initial level cells
    MFloat currentWorkload = weight(i);

    // save last pos
    MInt last = j;

    traverseDFGlobalId(i, &currentGlobalId, partitionCellList, workloadPerCell, &currentWorkload, weight, workload, &j);

    MLong noOffsprings = 0;
    if((i == 0) && (0 == m_noMissingParents))
      noOffsprings = currentGlobalId - m_cellOffsetPar;
    else
      noOffsprings = currentGlobalId - a_globalId(i);

    // update
    partitionCellList(last, 1) = noOffsprings;
    workloadPerCell(last) = currentWorkload;
  }

  m_log << "       * finding partition cells" << endl;
  outStream << "       * finding partition cells" << endl;

  // b. use only those elements that are min-cells
  for(MInt i = 0; i < m_noCells; i++) {
    MLong noOffsprings = partitionCellList(i, 1);

    // Added by Jerry for load balanced grid. Skipping cells with halo childs.
    MBool hasHaloChilds = false;
    for(MInt child = 0; child < m_maxNoChildren; ++child) {
      MInt t_child = (MInt)a_childId(i, child);
      if(t_child == -1) {
        continue;
      }
      if(a_hasProperty(t_child, 4) == 1) {
        hasHaloChilds = true;
        break;
      }
    }

    if(hasHaloChilds && g_dynamicLoadBalancing) {
      continue;
    }

    const MFloat work = workloadPerCell(i);

    MBool isSolverLeafCell = false;
    // Multisolver grid: assert that child cells of a cell that is a leaf cell for one solver are
    // not selected as partition cells to keep the entire subtree of the grid starting at that
    // solver-leaf-cell on one domain for coupling reasons
    if(m_noSolvers > 1) {
      const MLong tmpCellId = partitionCellList(i, 0);
      for(MInt solver = 0; solver < m_noSolvers; solver++) {
        if(!a_isInSolver((MInt)tmpCellId, solver)) { // Cell not relevant for solver
          continue;
        }

        // Check if there are child cells for this solver
        MBool hasSolverChild = false;
        for(MInt child = 0; child < m_maxNoChildren; ++child) {
          const MInt t_child = (MInt)a_childId((MInt)tmpCellId, child);
          if(t_child == -1) {
            continue;
          }
          if(a_isInSolver(t_child, solver)) {
            hasSolverChild = true;
            break;
          }
        }
        // If there is no child cell for this solver the cell is a solver leaf cell and is added to
        // the list of partition cells below
        if(!hasSolverChild) {
          m_log << "found solver leaf cell " << i << " id" << tmpCellId << " b" << solver << " l" << a_level(i) << " l"
                << a_level((MInt)tmpCellId) << std::endl;
          isSolverLeafCell = true;
          break;
        }
      }
    }

    if((noOffsprings <= m_partitionCellOffspringThreshold && work <= (MFloat)m_partitionCellWorkloadThreshold)
       || isSolverLeafCell) {
      if(isSolverLeafCell) {
        m_log << "partition solver leaf cell " << i << " " << partitionCellList(i, 0) << " l" << a_level(i) << " o"
              << noOffsprings << std::endl;
      }
      m_partitionCellList.push_back(make_tuple(partitionCellList(i, 0), noOffsprings, work));
      i += (MInt)noOffsprings;
    }
  }

  // c. update number of min-cells
  m_noPartitionCells = (MInt)m_partitionCellList.size();

  m_log << "           - found " << m_noPartitionCells << endl;
  outStream << "           - found " << m_noPartitionCells << endl;

  RECORD_TIMER_STOP(t_ro);
}

/** \brief recursively traverses the octree
 *
 * \author Andreas Lintermann
 * \date 03.12.2013
 *
 * Recursively calls itself and changes the global-id of the current cell.
 * It also fills the list partitionCellList with the according cell-id and the
 * number of offsprings.
 *
 * \param[in] parentId the id of the cell to traverse
 * \param[in] globalId the new global-id to use next
 * \param[in] partitionCellList a reference to the list containing all cell information
 * \param[in] workloadPerCell stores the work load per cell
 * \param[in] currentWorkload the current incremental work load
 * \param[in] weight the weighting
 * \param[in] workload base work load
 * \param[in] j pointer to the incrementer
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::traverseDFGlobalId(MInt parentId_, MLong* globalId_, MLongScratchSpace& partitionCellList,
                                          MFloatScratchSpace& workloadPerCell, MFloat* currentWorkload,
                                          MFloatScratchSpace& weight, MFloatScratchSpace& workload, MInt* j) {
  TRACE();

  MLong* children = &a_childId(parentId_, 0);

  for(MInt c = 0; c < m_maxNoChildren; c++) {
    if(children[c] > -1) {
      if(g_dynamicLoadBalancing
         && (a_hasProperty((MInt)children[c], 4)
             == 1)) { // Added by Jerry for load balanced grid. Skipping halo childs.
        continue;
      }

      // increase the weight by the weight of the children
      *currentWorkload = *currentWorkload + weight((MInt)children[c]);
      // increase global id by one
      *globalId_ = *globalId_ + 1;
      // save global id of child
      a_globalId((MInt)children[c]) = *globalId_;
      // set workload of child to current total work load or parent analog to the global Id
      workload((MInt)children[c]) = *currentWorkload;

      *j = *j + 1;
      partitionCellList(*j, 0) = children[c];

      // save last pos
      MInt last = *j;

      traverseDFGlobalId((MInt)children[c], globalId_, partitionCellList, workloadPerCell, currentWorkload, weight,
                         workload, j);

      MLong noOffsprings = 0;
      noOffsprings = *(globalId_)-a_globalId((MInt)children[c]);
      // subtract the total workload of the child from the current totally work load to compute the work load of the
      // child
      MFloat work = *currentWorkload - workload((MInt)children[c]) + weight((MInt)children[c]); // add own weight

      // update offsprings + workload of parent
      partitionCellList(last, 1) = noOffsprings;
      // set current work load of the child
      workloadPerCell(last) = work;
    }
  }
}

/**
 * \fn void GridgenPar<nDim>::setCellWeights(MFloatScratchSpace& weight)
 * \brief sets the cell weights according to the box system
 *
 * \author Stephan Schlimpert
 * \date 03.06.2014
 *
 * \todo labels:GRIDGEN if the weight method is suitable then use the box system like Andi's refinement method
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::setCellWeights(MFloatScratchSpace& weight) {
  TRACE();
  // chrs test the m_weight cell property, should be replaced by a meaningful function
  weight.fill(1);

  m_log << "applying weight method " << m_weightMethod << endl;
  switch(m_weightMethod) {
    case 1: {
      MInt noWeightCoordinates = Context::propertyLength("weightCoordinates");

      MFloatScratchSpace weightCoordinates(noWeightCoordinates, AT_, "weightCoordinates");

      for(MInt i = 0; i < noWeightCoordinates; i++) {
        if((i + 1) % 2 == 0)
          weightCoordinates[i] = 10000;
        else
          weightCoordinates[i] = -10000;
      }

      for(MInt i = 0; i < noWeightCoordinates; i++) {
        /*! \page propertyPage1
          \section weightCoordinates
          <code>MFloat CartesianGrid::m_weightCoordinates </code>\n
          default = <code>0</code>\n \n
          defines the weighht coordinates -x,+x,-y,+y,-z,+z direction:
          <ul>
          <li> Any non-negative floating point value smaller than the extent of the domain in the respective direction.
          </li>
          </ul>
          Keywords: <i>WEIGHT, COORDINATES</i>
         */
        weightCoordinates[i] = Context::getBasicProperty<MFloat>("weightCoordinates", AT_, &weightCoordinates.p[i], i);
        m_log << "weight coordinate [" << i << "] = " << weightCoordinates[i] << endl;
      }

      MInt noWeightValues = Context::propertyLength("setWeightValue");
      if(noWeightValues != (m_maxRfnmntLvl - m_minLevel + 1)) {
        stringstream errorMsg;
        errorMsg << "ERROR: noWeightValues not correct = " << noWeightValues;
        m_log << errorMsg.str() << endl;
        mTerm(1, AT_, errorMsg.str());
      }
      MFloatScratchSpace setWeightValue(noWeightValues, AT_, "setWeightValue");
      setWeightValue.fill(1);
      MInt levelCnt = m_maxRfnmntLvl;
      for(MInt i = 0; i < noWeightValues; i++) {
        /*! \page propertyPage1
          \section setWeightValue
          <code>MFloat GridGenpar::setWeightValue </code>\n
          default = <code>1</code>\n \n
          defines the weight values for the levels:
          <ul>
          <li> Any non-negative floating point value. </li>
          </ul>
          Keywords: <i>WEIGHT, VALUE</i>
         */
        setWeightValue(i) = Context::getBasicProperty<MFloat>("setWeightValue", AT_, &setWeightValue(i), i);
        m_log << "weight values [" << i << "] = " << setWeightValue(i) << " on level " << levelCnt << endl;
        levelCnt--;
      }

      for(MInt level_ = m_maxRfnmntLvl; level_ >= m_minLevel; --level_) {
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          if(a_coordinate(i, 0) > weightCoordinates[0] && a_coordinate(i, 0) < weightCoordinates[1]
             && a_coordinate(i, 1) > weightCoordinates[2] && a_coordinate(i, 1) < weightCoordinates[3]
             && a_coordinate(i, 2) > weightCoordinates[4] && a_coordinate(i, 2) < weightCoordinates[5]) {
            weight(i) = setWeightValue(m_maxRfnmntLvl - level_);
          }
        }
      }
      break;
    }
    case 2: {
      MInt noWeightCoordinates = Context::propertyLength("weightCoordinates");

      MFloatScratchSpace weightCoordinates(noWeightCoordinates, AT_, "weightCoordinates");

      for(MInt i = 0; i < noWeightCoordinates; i++) {
        if((i + 1) % 2 == 0)
          weightCoordinates[i] = 10000;
        else
          weightCoordinates[i] = -10000;
      }

      for(MInt i = 0; i < noWeightCoordinates; i++) {
        /*! \page propertyPage1
          \section weightCoordinates
          <code>MFloat CartesianGrid::m_weightCoordinates </code>\n
          default = <code>0</code>\n \n
          defines the weighht coordinates -x,+x,-y,+y,-z,+z direction:
          <ul>
          <li> Any non-negative floating point value smaller than the extent of the domain in the respective direction.
          </li>
          </ul>
          Keywords: <i>WEIGHT, COORDINATES</i>
         */
        weightCoordinates[i] = Context::getBasicProperty<MFloat>("weightCoordinates", AT_, &weightCoordinates.p[i], i);
        m_log << "weight coordinate [" << i << "] = " << weightCoordinates[i] << endl;
      }

      MInt noWeightValues = Context::propertyLength("setWeightValue");
      if(noWeightValues != (m_maxRfnmntLvl - m_minLevel + 1)) {
        stringstream errorMsg;
        errorMsg << "ERROR: noWeightValues not correct = " << noWeightValues;
        m_log << errorMsg.str() << endl;
        mTerm(1, AT_, errorMsg.str());
      }
      MFloatScratchSpace setWeightValue(noWeightValues, AT_, "setWeightValue");
      setWeightValue.fill(1);
      MInt levelCnt = m_maxRfnmntLvl;
      for(MInt i = 0; i < noWeightValues; i++) {
        /*! \page propertyPage1
          \section setWeightValue
          <code>MFloat GridGenpar::setWeightValue </code>\n
          default = <code>1</code>\n \n
          defines the weight values for the levels:
          <ul>
          <li> Any non-negative floating point value. </li>
          </ul>
          Keywords: <i>WEIGHT, VALUE</i>
         */
        setWeightValue(i) = Context::getBasicProperty<MFloat>("setWeightValue", AT_, &setWeightValue(i), i);
        m_log << "weight values [" << i << "] = " << setWeightValue(i) << " on level " << levelCnt << endl;
        levelCnt--;
      }

      MInt level_ = m_maxRfnmntLvl;
      for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
        if(a_coordinate(i, 0) > weightCoordinates[0] && a_coordinate(i, 0) < weightCoordinates[1]
           && a_coordinate(i, 1) > weightCoordinates[2] && a_coordinate(i, 1) < weightCoordinates[3]
           && a_coordinate(i, 2) > weightCoordinates[4] && a_coordinate(i, 2) < weightCoordinates[5]) {
          weight(i) = setWeightValue(m_maxRfnmntLvl - level_);
        }
      }
      break;
    }
    case 3: {
      MInt noWeightValues = Context::propertyLength("setWeightValue");
      if(noWeightValues != (m_maxRfnmntLvl - m_minLevel + 1)) {
        stringstream errorMsg;
        errorMsg << "ERROR: noWeightValues not correct = " << noWeightValues;
        m_log << errorMsg.str() << endl;
        mTerm(1, AT_, errorMsg.str());
      }
      MFloatScratchSpace setWeightValue(noWeightValues, AT_, "setWeightValue");
      setWeightValue.fill(1);
      MInt levelCnt = m_maxRfnmntLvl;
      for(MInt i = 0; i < noWeightValues; i++) {
        /*! \page propertyPage1
          \section setWeightValue
          <code>MFloat GridGenpar::setWeightValue </code>\n
          default = <code>1</code>\n \n
          defines the weight values for the levels:
          <ul>
          <li> Any non-negative floating point value. </li>
          </ul>
          Keywords: <i>WEIGHT, VALUE</i>
         */
        setWeightValue(i) = Context::getBasicProperty<MFloat>("setWeightValue", AT_, &setWeightValue(i), i);
        m_log << "weight values [" << i << "] = " << setWeightValue(i) << " on level " << levelCnt << endl;
        levelCnt--;
      }

      for(MInt level_ = m_maxRfnmntLvl; level_ >= m_minLevel; --level_) {
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          weight(i) = setWeightValue(m_maxRfnmntLvl - level_);
        }
      }
      break;
    }
    // weights dependent
    case 4: {
      MInt noWeightValues = Context::propertyLength("setWeightValue");
      if(noWeightValues != (m_maxRfnmntLvl - m_minLevel + 1)) {
        stringstream errorMsg;
        errorMsg << "ERROR: noWeightValues not correct = " << noWeightValues;
        m_log << errorMsg.str() << endl;
        mTerm(1, AT_, errorMsg.str());
      }
      MFloatScratchSpace setWeightValue(noWeightValues, AT_, "setWeightValue");
      setWeightValue.fill(1);
      MInt levelCnt = m_maxRfnmntLvl;
      for(MInt i = 0; i < noWeightValues; i++) {
        /*! \page propertyPage1
          \section setWeightValue
          <code>MFloat GridGenpar::setWeightValue </code>\n
          default = <code>1</code>\n \n
          defines the weight values for the levels:
          <ul>
          <li> Any non-negative floating point value. </li>
          </ul>
          Keywords: <i>WEIGHT, VALUE</i>
         */
        setWeightValue(i) = Context::getBasicProperty<MFloat>("setWeightValue", AT_, &setWeightValue(i), i);
        m_log << "weight values [" << i << "] = " << setWeightValue(i) << " on level " << levelCnt << endl;
        levelCnt--;
      }

      for(MInt level_ = m_maxRfnmntLvl; level_ >= m_minLevel; --level_) {
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          if(a_noChildren(i) == 0 && a_hasProperty(i, 1) == 1)
            weight(i) = setWeightValue(m_maxRfnmntLvl - level_);
          else if(a_noChildren(i) == 0)
            weight(i) = setWeightValue(m_maxRfnmntLvl - level_) + 0.3;
          else
            weight(i) = 1.0;
        }
      }
      break;
    }
    case 5: {
      // set a weight for all boundary cells

      for(MInt level_ = m_maxRfnmntLvl; level_ >= m_minLevel; --level_) {
        for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
          if(a_noChildren(i) > 0) {
            continue;
          }

          if(checkCellForCut(i)) {
            weight(i) = 5;
          }
        }
      }
      break;
    }
    case 6: {
      IF_CONSTEXPR(nDim == 3) {
        MFloat spawnCoord[nDim];
        // spray simulation weights
        if(Context::propertyLength("spawnCoordinates") != nDim) {
          TERMM(1, "Need to give a Coordinate for every dimension");
        }

        for(MInt i = 0; i < nDim; i++) {
          spawnCoord[i] = Context::getBasicProperty<MFloat>("spawnCoordinates", AT_, i);
        }

        for(MInt i = m_levelOffsets[m_maxRfnmntLvl][0]; i < m_levelOffsets[m_maxRfnmntLvl][1]; i++) {
          if(a_coordinate(i, 0) > spawnCoord[0] - 0.001 && a_coordinate(i, 0) < spawnCoord[0] + 0.001
             && a_coordinate(i, 1) > spawnCoord[1] - 0.001 && a_coordinate(i, 1) < spawnCoord[1] + 0.001
             && a_coordinate(i, 2) > spawnCoord[2] && a_coordinate(i, 2) < spawnCoord[2] + 0.02) {
            weight(i) = 100;
          }
        }
      }
      else {
        TERM(-1);
      }
      break;
    }
    default:
      stringstream errorMessage;
      errorMessage << "ERROR: weight method does not exist! ";
      mTerm(1, AT_, errorMessage.str());
      break;
  }
}

/**
 * \brief Retrieves all direct and diagonal neighboring cells of the given cell
 *
 * \author Lennart Schneiders, Andreas Lintermann
 * \date 27.11.2017
 *
 * \param[in] cellId the id of the cell to return the neighbors from
 * \param[in] adjacentCells pointer to an array that will be filled with the results
 *
 **/
template <MInt nDim>
MInt GridgenPar<nDim>::getAdjacentGridCells(MInt cellId, MInt* adjacentCells) {
  MInt cnt = 0;
  set<MLong> nghbrs;
  for(MInt dir0 = 0; dir0 < 2 * nDim; dir0++) {
    MLong nghbrId0 = -1;
    if(a_neighborId(cellId, dir0) > -1) nghbrId0 = a_neighborId(cellId, dir0);
    if(nghbrId0 < 0) continue;
    nghbrs.insert(nghbrId0);
    for(MInt dir1 = 0; dir1 < 2 * nDim; dir1++) {
      if((dir1 / 2) == (dir0 / 2)) continue;
      MLong nghbrId1 = -1;
      if(a_neighborId((MInt)nghbrId0, dir1) > -1) nghbrId1 = a_neighborId((MInt)nghbrId0, dir1);
      if(nghbrId1 < 0) continue;
      nghbrs.insert(nghbrId1);
      IF_CONSTEXPR(nDim == 3) {
        for(MInt dir2 = 0; dir2 < 2 * nDim; dir2++) {
          if(((dir2 / 2) == (dir0 / 2)) || ((dir2 / 2) == (dir1 / 2))) continue;
          MLong nghbrId2 = -1;
          if(a_neighborId((MInt)nghbrId1, dir2) > -1) nghbrId2 = a_neighborId((MInt)nghbrId1, dir2);
          if(nghbrId2 < 0) continue;
          nghbrs.insert(nghbrId2);
        }
      }
    }
  }
  set<MLong>::iterator it = nghbrs.begin();
  for(it = nghbrs.begin(); it != nghbrs.end(); it++) {
    ASSERT(cnt < 27, "");
    adjacentCells[cnt] = (MInt)*it;
    cnt++;
  }
  return cnt;
}


/** \brief checks if this is a valid mesh also for LB computations
 *
 * \author Andreas Lintermann
 * \date 23.11.2017
 *
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::checkLBRefinementValidity() {
  TRACE();

  NEW_SUB_TIMER_STATIC(t_checkLBRefinementValidity, "check refinement valdity for LB", m_t_finalizeGrid);
  RECORD_TIMER_START(t_checkLBRefinementValidity);

  outStream << "     + checking mesh validity for LB" << endl;
  m_log << "     + checking mesh validity for LB" << endl;

  // run over all cells on the lowest tree level
  // find out if the parent has all children
  // if not, check if the parent's neighbor is a leaf cell

  set<MInt> lberrcells;


  for(MInt i = 0; i < m_noCells; ++i) {
    // find leave cells
    if(a_noChildren(i) == 0) {
      MInt parent = (MInt)a_parentId(i);

      // check if the parent has a missing child
      if(parent >= 0 && a_noChildren(parent) != m_maxNoChildren) {
        MIntScratchSpace nghbrList(27, AT_, "nghbrList");
        for(MInt l = 0; l < 27; l++)
          nghbrList[l] = -1;

        const MInt counter = getAdjacentGridCells(parent, nghbrList.begin());
        for(MInt n = 0; n < counter; n++) {
          if(nghbrList[n] > -1 && a_noChildren(nghbrList[n]) == 0 && !a_hasProperty(nghbrList[n], 5)) {
            lberrcells.insert(parent);
            break;
          }
        }
      }
    }
  }


  MInt errorcells = 0;
  MInt l_lberrcells = (MInt)lberrcells.size();

  MPI_Allreduce(&l_lberrcells, &errorcells, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "l_lberrcells", "errorcells");

  if(errorcells > 0) {
    outStream << "       * \033[1;31mthere are " << errorcells << " errorneous cells\033[0m" << endl;
    outStream
        << "         - \033[1;31mEXPLANATION: there are cells that will be interface cells in the LB computation and "
           "which have missing children\033[0m"
        << endl;
    outStream << "         - \033[1;31mHINTS      : try to refine in- and outlets\033[0m" << endl;
    outStream
        << "                        \033[1;31mmake sure you chose a good combination of localMinBoundaryThreshold and "
           "smoothDistance\033[0m"
        << endl;
    m_log << "       * there are " << errorcells << " errorneous cells" << endl;
    m_log << "         - EXPLANATION: there are cells that will be interface cells in the LB computation and which "
             "have missing children"
          << endl;
    m_log << "         - HINTS      : try to refine in- and outlets" << endl;
    m_log << "                        make sure you chose a good combination of localMinBoundaryThreshold and "
             "smoothDistance"
          << endl;
  } else {
    outStream << "       * mesh is \033[1;32mOK\033[0m" << endl;
    m_log << "       * mesh is OK" << endl;
  }


  RECORD_TIMER_STOP(t_checkLBRefinementValidity);
}


/** \brief refines cells of the computational grid that are marked for refinement
 *
 * \author Andreas Lintermann
 * \date 6.5.2014
 *
 *     6.2a.2.2 check memory availability
 *     6.2a.2.3 Refine the cells that have previously been marked. This calls the function
 *              refineGridPatch(...) which only refines cells that have previously been
 *              marked. If we run in parallel mode, do this also for the halo cells.
 *     6.2b.2.4 Update the number of cells.
 *     6.2b.2.5 Find the neighbors for the new cells. Calls findChildLevelNeighbors(...).
 *              If we run in parallel mode, this is also done for the halo cells.
 *     6.2b.2.6 Delete all outside cells, do this in serial or in parallel. This calls
 *              either deleteOutsideCellsSerial(...) or deleteOutsideCellsParallel(...).
 **/
template <MInt nDim>
void GridgenPar<nDim>::refineComputationalGrid(MInt compRfnLvl) {
  TRACE();

  // 6.2a.2.2 check memory availability
  checkMemoryAvailability(1, compRfnLvl + 1);

  // 6.2a.2.3 refine the cells and the halo cells if we are parallel
  refineGridPatch(m_levelOffsets, compRfnLvl, 0);
  if(noDomains() > 1) {
    refineGridPatch(m_haloCellOffsetsLevel, compRfnLvl, 1);
  }

  // 6.2a.2.4 update the number of cells
  m_noCells = m_levelOffsets[compRfnLvl + 1][1];

  // 6.2a.2.5 find the neighbors on this level, do this also for the halos if we run in parallel
  findChildLevelNeighbors(m_levelOffsets, compRfnLvl);
  if(noDomains() > 1) {
    findChildLevelNeighbors(m_haloCellOffsetsLevel, compRfnLvl);
  }

  // 6.2a.2.6 delete all outside cells, do this in serial or in parallel
  if(noDomains() > 1) {
    deleteOutsideCellsParallel(compRfnLvl + 1);
  } else {
    deleteOutsideCellsSerial(compRfnLvl + 1);
  }

  // update number of cells
  m_noCells = m_levelOffsets[compRfnLvl + 1][1];

  // 6.2a.2.9 in parallel, check if there is a load imbalance in the grid before continuing
  if((noDomains() > 1) && (g_dynamicLoadBalancing)) { // TODO labels:GRIDGEN,DLB skip for last level?
    checkLoadBalance(compRfnLvl);
  }
}

/** \brief propagates the distance away from the boundary
 *
 * \author Andreas Lintermann
 * \date 06.11.2013
 *
 * This function propagates distance information starting at the boundaries and
 * stores them in the local variable m_rfnDistance:
 *
 *  a. Start by collecting all boundary cells. Unfortunately, we have to do an
 *     intersection test again, since we probably want to refine only certain
 *     boundaries (defined by the property localRfnBoundaryIds). The collected
 *     cells are stored in the vector boundaryCells.
 *  b. Now start to recursivley mark the cells by their distance. This calls
 *     progagationStep(...), which recursively runs over all cells and checks for
 *     new unmarked cells. For the serial case, the algorithm is finished at
 *     this point.
 *  c. Otherwise, if we have a parallel calculation, we need to check if we have
 *     to cross a domain boundary with our propagation.
 *  d. Exchange information as long as we have local changes.
 *     determine whether we have changes by comparing a_refinementDistance
 *     on halo and window cells
 *     this has to be repeated until we are finished, exchange information as long
 *     as we have local changes
 *
 * \param[in] level_, the level to run the propagation on
 * \param[in] distance, the final distance for this level
 * \param[in] rfnBoundaryGroup, group of boundaries to refine
 * \param[in] solver, the solver to check for refinement
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::propagateDistance(MInt level_, MInt distance, std::vector<MInt>& rfnBoundaryGroup, MInt solver) {
  TRACE();

  NEW_SUB_TIMER_STATIC(t_boundaryPropagation, "boundary propagation", m_t_createComputationalGrid);
  RECORD_TIMER_START(t_boundaryPropagation);

  m_log << "           marking boundary cells: ";

  // a. collect all boundary cells on this domain an init the refinementDistance
  vector<MInt> boundaryCells;
  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    a_refinementDistance(i) = numeric_limits<MInt>::max();
    if(a_isSolverBoundary(i, solver)) {
      checkCellForCut(i);

      for(MInt b = 0; b < (MInt)rfnBoundaryGroup.size(); b++) {
        if(m_bndCutInfo[solver][rfnBoundaryGroup[b]]) {
          boundaryCells.push_back(i);
          break;
        }
      }
    }
  }

  if(noDomains() > 1) { // i do not really need this, infos would propagate from other domains
    MInt lev_pos = 2 * (level_ - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
        a_refinementDistance(i) = numeric_limits<MInt>::max();
        if(a_isSolverBoundary(i, solver)) {
          checkCellForCut(i);

          for(MInt b = 0; b < (MInt)rfnBoundaryGroup.size(); b++) {
            if(m_bndCutInfo[solver][rfnBoundaryGroup[b]]) {
              boundaryCells.push_back(i);
              break;
            }
          }
        }
      }
    }
  }

  // b. do the recursive marking
  for(MInt i = 0; i < (MInt)boundaryCells.size(); i++)
    propagationStep(boundaryCells[i], 0, distance, solver);

  // if we are parallel, we have to transmit information
  if(noDomains() > 1) { // TODO labels:GRIDGEN the communication setup should be done once per level or only account for
                        // solver window/halos

    // c. create lists of window and halo cells
    vector<vector<MInt>> winCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));
    vector<vector<MInt>> haloCellIdsPerDomain(m_noNeighborDomains, vector<MInt>(0));

    if(m_hasBeenLoadBalanced)
      findHaloAndWindowCellsKD(winCellIdsPerDomain, haloCellIdsPerDomain);
    else
      findHaloAndWindowCells(winCellIdsPerDomain, haloCellIdsPerDomain);

    // d. prepare the communication
    MIntScratchSpace noSendWindowPerDomain(m_noNeighborDomains, AT_, "noSendWindowPerDomain");
    MIntScratchSpace noReceiveHaloPerDomain(m_noNeighborDomains, AT_, "noReceiveHaloPerDomain");

    m_log << "         - we need to transfer to domain [no. cells to transfer]: ";
    outStream << "         - we need to transfer to domain [no. cells to transfer]: ";
    for(MInt d = 0; d < m_noNeighborDomains; d++) {
      noSendWindowPerDomain[d] = (MInt)winCellIdsPerDomain[d].size();
      m_log << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
      outStream << m_neighborDomains[d] << " [" << noSendWindowPerDomain[d] << "]   ";
    }
    m_log << endl;
    outStream << endl;

    m_log << "         - we need to receive from domain [no. cells to receive]: ";
    outStream << "         - we need to receive from domain [no. cells to receive]: ";
    for(MInt d = 0; d < m_noNeighborDomains; d++) {
      noReceiveHaloPerDomain[d] = (MInt)haloCellIdsPerDomain[d].size();
      m_log << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
      outStream << m_neighborDomains[d] << " [" << noReceiveHaloPerDomain[d] << "]   ";
    }
    m_log << endl;
    outStream << endl;

    MInt allSend = 0;
    MInt allReceive = 0;
    vector<MInt> offsetsSend;
    vector<MInt> offsetsReceive;
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      offsetsSend.push_back(allSend);
      offsetsReceive.push_back(allReceive);
      allSend += noSendWindowPerDomain[dom];
      allReceive += noReceiveHaloPerDomain[dom];
    }

    MIntScratchSpace sndBufWin(allSend, AT_, "sndBufWin");
    MIntScratchSpace rcvBufHalo(allReceive, AT_, "rcvBufHalo");

    MUshort finished = 0;

    // e. determine whether we have changes by comparing a_refinementDistance
    //   on halo and window cells
    //   this has to be repeated until we are finished, exchange information as long
    //   as we have local changes
    MInt rounds = 0;
    while(!finished) {
      finished = 1;

      MPI_Request* mpi_request_ = nullptr;
      mAlloc(mpi_request_, m_noNeighborDomains, "mpi_request_", AT_);
      for(MInt d = 0; d < m_noNeighborDomains; d++) {
        for(MInt c = 0; c < noSendWindowPerDomain[d]; c++) {
          sndBufWin[offsetsSend[d] + c] = a_refinementDistance(winCellIdsPerDomain[d][c]);
        }
        MPI_Issend(&(sndBufWin[offsetsSend[d]]), noSendWindowPerDomain[d], MPI_INT, m_neighborDomains[d], 0, mpiComm(),
                   &mpi_request_[d], AT_, "(sndBufWin[offsetsSend[d]])");
      }

      MPI_Status status_;
      for(MInt d = 0; d < m_noNeighborDomains; d++)
        MPI_Recv(&(rcvBufHalo[offsetsReceive[d]]), noReceiveHaloPerDomain[d], MPI_INT, m_neighborDomains[d], 0,
                 mpiComm(), &status_, AT_, "(rcvBufHalo[offsetsReceive[d]])");

      for(MInt d = 0; d < m_noNeighborDomains; d++)
        MPI_Wait(&mpi_request_[d], &status_, AT_);

      for(MInt d = 0; d < m_noNeighborDomains; d++) {
        for(MInt c = 0; c < noReceiveHaloPerDomain[d]; c++) {
          const MInt halo = haloCellIdsPerDomain[d][c];
          if(rcvBufHalo[c + offsetsReceive[d]] < a_refinementDistance(halo)) {
            propagationStep(halo, rcvBufHalo[c + offsetsReceive[d]], distance, solver);
            finished = 0;
          }
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &finished, 1, MPI_UNSIGNED_SHORT, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE",
                    "finished");

      rounds++;
    }
    m_log << "         - communication rounds required: " << rounds << endl;
  }

  RECORD_TIMER_STOP(t_boundaryPropagation);
}

/** \brief recursivley marks the cells with a distance
 *
 * \author Andreas Lintermann
 * \date 01.11.2013
 *
 * This is a recursive function, which is initially called with a boundary cell-id.
 * For this cell, the neighborhood is then recursivley investigated. If a neighbor cell
 * has not been visited, yet (m_rfnDistance < 0) it is traversed and set to an incremental
 * distance. If a neighbor has a distance smaller than the incremented distance it is not
 * traversed.
 *
 * \param[in] cellId the current cell-id to check
 * \param[in] rfnDistance the current distance to use
 * \param[in] finalDistance the final distance to use
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::propagationStep(MInt cellId, MInt rfnDistance, MInt finalDistance, MInt solver) {
  TRACE();

  NEW_SUB_TIMER_STATIC(t_markRfnDistance, "mark refinement distance", m_t_createComputationalGrid);
  if(rfnDistance == 0) {
    RECORD_TIMER_START(t_markRfnDistance);
  }

  if(a_refinementDistance(cellId) > rfnDistance) {
    a_refinementDistance(cellId) = rfnDistance;
    // Definition: We do not refine at final distance!
    if(rfnDistance < (finalDistance - 1)) {
      for(MInt n = 0; n < m_noNeighbors; n++) {
        if(a_neighborId(cellId, n) > -1) {
          if(a_isInSolver((MInt)a_neighborId(cellId, n), solver)) {
            propagationStep((MInt)a_neighborId(cellId, n), rfnDistance + 1, finalDistance, solver);
          }
        }
      }
    }
  }

  if(rfnDistance == 0) {
    RECORD_TIMER_STOP(t_markRfnDistance);
  }
}

/** \brief marks cells that lie in a box-type patch
 *
 * \author Andreas Lintermann
 * \date 04.12.2013
 *
 * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over alls
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given box) and marks these cells
 *     (b_properties[5] = 1).
 *  b. the same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level_ the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalBox(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  // a. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local box for solver " << solver << endl;
  m_log << "         - corner1: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl << "         - corner2: ";
  for(MInt d = nDim; d < 2 * nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl;

  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      if(maia::geom::isPointInsideBox<MFloat, nDim>(&a_coordinate(k, 0), bp->localRfnPatchProperties[pos])) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // b. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          if(maia::geom::isPointInsideBox<MFloat, nDim>(&a_coordinate(k, 0), bp->localRfnPatchProperties[pos])) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}

/** \brief marks cells that lie in a sphere-type patch
 *
 * \author Andreas Lintermann
 * \date 04.12.2013
 *
 * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over all
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given sphere) and marks these cells
 *     (b_properties[5] = 1).
 *  b. the same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalRadius(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  // a. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  // Get coordinates
  const MFloat* sphereCoord = bp->localRfnPatchProperties[pos];

  // Get radius
  const MFloat radius = bp->localRfnPatchProperties[pos][nDim];

  m_log << "       * marking local sphere for solver " << solver << endl;
  m_log << "         - center: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << sphereCoord[d] << " ";
  }
  m_log << endl << "         - radius: " << radius << endl;

  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      if(maia::geom::isPointInsideSphere<nDim>(&a_coordinate(k, 0), sphereCoord, radius)) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // b. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          if(maia::geom::isPointInsideSphere<nDim>(&a_coordinate(k, 0), sphereCoord, radius)) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}

/** \brief marks cells that lie in a cylinder-type patch
 *
 * \author Andreas Lintermann
 * \date 04.12.2013
 *
 * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over alls
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given cylinder) and marks these cells
 *     (b_properties[5] = 1).
 *     The according formula for the calculation of the distance \f$d\f$
 *     of a given point \f$\vec{p}\f$ (cell coordinates) is given by
 *
 *     \f[d = \left|\left(\vec{p}-\vec{o}_1\right)\times\vec{r}\right|,\f]
 *
 *     where \f$\vec{o}_1\f$ is the first center-point provided in the property,
 *     \f$\vec{r}\f$ is the normalized vector between \f$\vec{o}_1\f$ and
 *     \f$\vec{o}_2\f$ (as provieded in thr property file).
 *
 *  b. The same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalCylinder(MInt compRfnLvl, MInt patch, MInt solver, MString patchStr) {
  TRACE();

  if(patchStr != "C" && patchStr != "T") {
    TERMM(1, "This function doesn't support the selected patch type!");
  }

  const MBool isTube = (patchStr == "T") ? true : false;

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  if(isTube) {
    m_log << "       * marking local tube for solver " << solver << endl;
  } else {
    m_log << "       * marking local cylinder for solver " << solver << endl;
  }

  // Get coordinates
  std::array<MFloat, nDim> leftCoord;
  std::array<MFloat, nDim> rightCoord;
  for(MInt d = 0; d < nDim; d++) {
    leftCoord[d] = bp->localRfnPatchProperties[pos][d];
    rightCoord[d] = bp->localRfnPatchProperties[pos][d + nDim];
  }

  m_log << "         - left center: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << leftCoord[d] << " ";
  }
  m_log << "         - right center: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << rightCoord[d] << " ";
  }
  m_log << endl;
  MFloat radius = bp->localRfnPatchProperties[pos][2 * nDim];
  m_log << "         - radius: " << radius << endl;
  MFloat innerRadius = -1.0; // for cylinder, distance is always > innerRadius
  if(isTube) {
    m_log << "         - inner radius: " << bp->localRfnPatchProperties[pos][2 * nDim + 1] << endl;
    innerRadius = bp->localRfnPatchProperties[pos][2 * nDim + 1];
  }

  // a. precalculate the required vector
  std::array<MFloat, nDim> normalDiffAB;
  std::array<MFloat, nDim> normalDiffBA;
  MFloat length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normalDiffAB[d] = rightCoord[d] - leftCoord[d];
    normalDiffBA[d] = leftCoord[d] - rightCoord[d];
    length += normalDiffAB[d] * normalDiffAB[d];
  }
  length = sqrt(length);

  for(MInt d = 0; d < nDim; d++) {
    normalDiffAB[d] /= length;
  }

  // b. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      if(isInsideCylinder(&a_coordinate(k, 0), leftCoord.data(), rightCoord.data(), normalDiffAB.data(),
                          normalDiffBA.data(), radius, innerRadius)) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // c. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          if(isInsideCylinder(&a_coordinate(k, 0), leftCoord.data(), rightCoord.data(), normalDiffAB.data(),
                              normalDiffBA.data(), radius, innerRadius)) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}


/** \brief checks if point (cell) is inside a cylinder
 *
 * \author  Andreas Lintermann
 * \date 04.12.2013
 *
 * \param[in] pointCoord cell coordinate
 * \param[in] (left/right)Coord coordinates of the left and right circles
 * \param[in] normalDiff(AB/BA) cylinder axis
 * \param[in] radius radius of the cylinder
 * \param[in] innerRadius inner radius of the cylinder (tube)
 *
 */
template <MInt nDim>
MBool GridgenPar<nDim>::isInsideCylinder(const MFloat* const pointCoord, const MFloat* const leftCoord,
                                         const MFloat* const rightCoord, const MFloat* const normalDiffAB,
                                         const MFloat* const normalDiffBA, const MFloat radius,
                                         const MFloat innerRadius) {
  TRACE();

  std::array<MFloat, nDim> diffA;
  std::array<MFloat, nDim> diffB;

  MFloat s1 = 0.0;
  MFloat s2 = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    diffA[d] = pointCoord[d] - leftCoord[d];
    diffB[d] = pointCoord[d] - rightCoord[d];
    s1 += diffA[d] * normalDiffAB[d];
    s2 += diffB[d] * normalDiffBA[d];
  }

  std::array<MFloat, nDim> cross;
  MFloat distance = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    MInt p1 = (d + 1) % (nDim);
    MInt p2 = (d + 2) % (nDim);

    cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
    distance += cross[d] * cross[d];
  }
  distance = sqrt(distance);

  if(distance <= radius && distance >= innerRadius && s1 >= 0.0 && s2 >= 0.0) {
    return 1;
  }
  return 0;
}

template <MInt nDim>
void GridgenPar<nDim>::markLocalCartesianWedge(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local cylinder for solver " << solver << endl;
  m_log << "         - min center: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl;
  m_log << "         - radius: " << bp->localRfnPatchProperties[pos][nDim] << endl;
  m_log << "         - length: " << bp->localRfnPatchProperties[pos][nDim + 1] << endl;
  m_log << "         - axis: " << bp->localRfnPatchProperties[pos][nDim + 2] << endl;
  m_log << "         - start angle: " << bp->localRfnPatchProperties[pos][nDim + 3] << endl;
  m_log << "         - end angle: " << bp->localRfnPatchProperties[pos][nDim + 4] << endl;

  MFloat length = bp->localRfnPatchProperties[pos][nDim + 1];
  MInt axis = (MInt)(bp->localRfnPatchProperties[pos][nDim + 2]);

  if(axis >= nDim) mTerm(1, AT_, "cylinder wedge refinement: axis out of bounds! Must be < nDim");

  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      MFloat* coords = &a_coordinate(k, 0);

      // check if cell is within cylinder axis range
      if(coords[axis] < bp->localRfnPatchProperties[pos][axis]
         || coords[axis] > bp->localRfnPatchProperties[pos][axis] + length)
        continue;

      // check if cell is within radius
      std::array<MFloat, nDim> normalVector;
      MFloat distance = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        if(d == axis)
          normalVector[d] = 0.0;
        else
          normalVector[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
        distance += normalVector[d] * normalVector[d];
      }
      distance = sqrt(distance);
      if(distance > bp->localRfnPatchProperties[pos][nDim]) continue;

      // determine angle
      // MInt relativeDirection = (axis + 1) % nDim;
      // MFloat angle = asin(normalVector[relativeDirection])*180.0/PI;
      MFloat angle = atan2(normalVector[(axis + 2) % nDim], normalVector[(axis + 1) % nDim]) * 180.0 / PI;

      if(angle < bp->localRfnPatchProperties[pos][nDim + 3] || angle > bp->localRfnPatchProperties[pos][nDim + 4])
        continue;

      a_isToRefineForSolver(k, solver) = 1;
    }
  }

  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          MFloat* coords = &a_coordinate(k, 0);

          if(coords[axis] < bp->localRfnPatchProperties[pos][axis]
             || coords[axis] > bp->localRfnPatchProperties[pos][axis] + length)
            continue; // outside of cylinder length

          std::array<MFloat, nDim> normalVector;
          MFloat distance = 0.0;
          for(MInt d = 0; d < nDim; d++) {
            if(d == axis)
              normalVector[d] = 0.0;
            else
              normalVector[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
            distance += normalVector[d] * normalVector[d];
          }
          distance = sqrt(distance);
          if(distance > bp->localRfnPatchProperties[pos][nDim]) continue; // outside of radius

          // MInt relativeDirection = (axis + 1) % nDim;
          // MFloat angle = asin(normalVector[relativeDirection])*180.0/PI;
          MFloat angle = atan2(normalVector[(axis + 2) % nDim], normalVector[(axis + 1) % nDim]) * 180.0 / PI;

          if(angle < bp->localRfnPatchProperties[pos][nDim + 3] || angle > bp->localRfnPatchProperties[pos][nDim + 4])
            continue; // outside of angular section

          a_isToRefineForSolver(k, solver) = 1;
        }
      }
    }
  }
}


// Volkan
/** \brief marks cells that lie in a rectangular-angled-type patch
 * * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over alls
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given cylinder) and marks these cells
 *     (b_properties[5] = 1).
 *     The according formula for the calculation of the distance \f$d\f$
 *     of a given point \f$\vec{p}\f$ (cell coordinates) is given by
 *
 *     \f[d = \left|\left(\vec{p}-\vec{o}_1\right)\times\vec{r}\right|,\f]
 *
 *     where \f$\vec{o}_1\f$ is the first center-point provided in the property,
 *     \f$\vec{r}\f$ is the normalized vector between \f$\vec{o}_1\f$ and
 *     \f$\vec{o}_2\f$ (as provieded in thr property file).
 *
 *  b. The same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalRectangleAngled(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local angled rectangle for solver " << solver << endl;
  m_log << "         - left center: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << "         - right center: ";
  for(MInt d = nDim; d < 2 * nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl;
  m_log << "         - height: " << bp->localRfnPatchProperties[pos][2 * nDim] << endl;
  m_log << "         - width: " << bp->localRfnPatchProperties[pos][2 * nDim + 1] << endl;


  // a. precalculate the required vector

  std::array<MFloat, nDim> normalDiffAB;
  std::array<MFloat, nDim> normalDiffBA;
  MFloat length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normalDiffAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - bp->localRfnPatchProperties[pos][d];
    normalDiffBA[d] = bp->localRfnPatchProperties[pos][d] - bp->localRfnPatchProperties[pos][d + nDim];
    length += normalDiffAB[d] * normalDiffAB[d];
  }

  length = sqrt(length);

  for(MInt d = 0; d < nDim; d++) {
    normalDiffAB[d] /= length;
  }

  // b. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      MFloat* coords = &a_coordinate(k, 0);

      std::array<MFloat, nDim> diffA;
      std::array<MFloat, nDim> diffB;

      MFloat s1 = 0.0;
      MFloat s2 = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        diffA[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
        diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
        s1 += diffA[d] * normalDiffAB[d];
        s2 += diffB[d] * normalDiffBA[d];
      }

      std::array<MFloat, nDim> cross;
      MFloat distanceH = 0.0;
      MFloat distanceW = 0.0;


      for(MInt d = 0; d < nDim; d++) {
        MInt p1 = (d + 1) % (nDim);
        MInt p2 = (d + 2) % (nDim);

        cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
        if(d == 0) {
          distanceH += cross[d] * cross[d];
        }
        if(d == 1) {
          distanceH += cross[d] * cross[d];
        }
        if(d == 2) {
          distanceW += cross[d] * cross[d];
        }
      }
      distanceH = sqrt(distanceH);
      distanceW = sqrt(distanceW);

      if(distanceH <= bp->localRfnPatchProperties[pos][2 * nDim]
         && distanceW <= bp->localRfnPatchProperties[pos][2 * nDim + 1] && s1 >= 0.0 && s2 >= 0.0) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // c. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          MFloat* coords = &a_coordinate(k, 0);

          std::array<MFloat, nDim> diffA;
          std::array<MFloat, nDim> diffB;

          MFloat s1 = 0.0;
          MFloat s2 = 0.0;
          for(MInt d = 0; d < nDim; d++) {
            diffA[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
            diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
            s1 += diffA[d] * normalDiffAB[d];
            s2 += diffB[d] * normalDiffBA[d];
          }

          std::array<MFloat, nDim> cross;
          MFloat distanceH = 0.0;
          MFloat distanceW = 0.0;


          for(MInt d = 0; d < nDim; d++) {
            MInt p1 = (d + 1) % (nDim);
            MInt p2 = (d + 2) % (nDim);

            cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
            if(d == 0) {
              distanceH += cross[d] * cross[d];
            }
            if(d == 1) {
              distanceH += cross[d] * cross[d];
            }
            if(d == 2) {
              distanceW += cross[d] * cross[d];
            }
          }
          distanceH = sqrt(distanceH);
          distanceW = sqrt(distanceW);

          if(distanceH <= bp->localRfnPatchProperties[pos][2 * nDim]
             && distanceW <= bp->localRfnPatchProperties[pos][2 * nDim + 1] && s1 >= 0.0 && s2 >= 0.0) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}

/** \brief marks cells that lie in a cone-type patch with a smooth hat
 *
 * \author Andreas Lintermann
 * \date 12.12.2013
 *
 * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over alls
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given cone and sphere) and marks
 *     these cells (b_properties[5] = 1).
 *     First the center of the sphere and the new center point of the
 *     cone are calculated. Then, the distance of a coordinates of a
 *     given cell is calculated to obatain the enclosed angle between
 *     this the vector of the coordinates and the first provided point.
 *     In case the angle is larger than the one provided in the
 *     properties, this cell is marked as outside, otherwise as inside.
 *     The distance of the coordinates to the line is calculated
 *     as described in markLocalCylinder.
 *  b. The same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalCone(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local rounded cone for solver " << solver << endl;
  m_log << "         - left point: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << "         - right point: ";
  for(MInt d = nDim; d < 2 * nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl;
  m_log << "         - opening angle: " << bp->localRfnPatchProperties[pos][2 * nDim] << endl;
  m_log << "         - radius: " << bp->localRfnPatchProperties[pos][2 * nDim + 1] << endl;

  // a. precalculate the required vector

  MFloat Radius = bp->localRfnPatchProperties[pos][2 * nDim + 1];
  MFloat angle = bp->localRfnPatchProperties[pos][2 * nDim];
  MFloat angle_rad = angle * PI / 180.0;
  //  MFloat angle_rad90 = (90.0-angle) * PI / 180.0;
  MFloat radius = Radius * cos(angle_rad);

  std::array<MFloat, nDim> normalDifforigAB;
  MFloat length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normalDifforigAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - bp->localRfnPatchProperties[pos][d];
    length += normalDifforigAB[d] * normalDifforigAB[d];
  }

  length = sqrt(length);

  for(MInt d = 0; d < nDim; d++)
    normalDifforigAB[d] /= length;

  // calucate the new position of A
  std::array<MFloat, nDim> newPosA;
  std::array<MFloat, nDim> sphereCenter;
  std::array<MFloat, nDim> normalDiffAB;
  std::array<MFloat, nDim> normalDiffBA;
  MFloat lengthAB = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    newPosA[d] = bp->localRfnPatchProperties[pos][d] + radius / tan(angle_rad) * normalDifforigAB[d];
    sphereCenter[d] =
        bp->localRfnPatchProperties[pos][d] + radius * (1.0 / tan(angle_rad) + tan(angle_rad)) * normalDifforigAB[d];

    normalDiffAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - newPosA[d];
    normalDiffBA[d] = newPosA[d] - bp->localRfnPatchProperties[pos][d + nDim];

    lengthAB += normalDiffAB[d] * normalDiffAB[d];
  }

  lengthAB = sqrt(lengthAB);

  for(MInt d = 0; d < nDim; d++) {
    normalDiffAB[d] /= lengthAB;
  }

  // b. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      MFloat* coords = &a_coordinate(k, 0);

      std::array<MFloat, nDim> diffA;
      std::array<MFloat, nDim> difforigA;
      std::array<MFloat, nDim> diffB;
      std::array<MFloat, nDim> diffCenter;

      MFloat s1 = 0.0;
      MFloat s2 = 0.0;
      MFloat len_diffOrig = 0.0;
      MFloat len_diffCenter = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        diffA[d] = coords[d] - newPosA[d];
        diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
        difforigA[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
        diffCenter[d] = coords[d] - sphereCenter[d];

        s1 += diffA[d] * normalDiffAB[d];
        s2 += diffB[d] * normalDiffBA[d];
        len_diffOrig += difforigA[d] * difforigA[d];
        len_diffCenter += diffCenter[d] * diffCenter[d];
      }

      len_diffOrig = sqrt(len_diffOrig);
      len_diffCenter = sqrt(len_diffCenter);

      std::array<MFloat, nDim> cross;
      MFloat distance = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        MInt p1 = (d + 1) % (nDim);
        MInt p2 = (d + 2) % (nDim);

        cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
        distance += cross[d] * cross[d];
      }

      distance = sqrt(distance);

      MFloat ang = asin(distance / len_diffOrig) * 180 / PI;

      if((ang <= angle && s1 >= 0.0 && s2 >= 0.0) || len_diffCenter <= Radius) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // b. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);

    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          MFloat* coords = &a_coordinate(k, 0);

          std::array<MFloat, nDim> diffA;
          std::array<MFloat, nDim> difforigA;
          std::array<MFloat, nDim> diffB;
          std::array<MFloat, nDim> diffCenter;

          MFloat s1 = 0.0;
          MFloat s2 = 0.0;
          MFloat len_diffOrig = 0.0;
          MFloat len_diffCenter = 0.0;
          for(MInt d = 0; d < nDim; d++) {
            diffA[d] = coords[d] - newPosA[d];
            diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
            difforigA[d] = coords[d] - bp->localRfnPatchProperties[pos][d];
            diffCenter[d] = coords[d] - sphereCenter[d];

            s1 += diffA[d] * normalDiffAB[d];
            s2 += diffB[d] * normalDiffBA[d];
            len_diffOrig += difforigA[d] * difforigA[d];
            len_diffCenter += diffCenter[d] * diffCenter[d];
          }

          len_diffOrig = sqrt(len_diffOrig);
          len_diffCenter = sqrt(len_diffCenter);

          std::array<MFloat, nDim> cross;
          MFloat distance = 0.0;
          for(MInt d = 0; d < nDim; d++) {
            MInt p1 = (d + 1) % (nDim);
            MInt p2 = (d + 2) % (nDim);

            cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
            distance += cross[d] * cross[d];
          }

          distance = sqrt(distance);

          MFloat ang = asin(distance / len_diffOrig) * 180 / PI;

          if((ang <= angle && s1 >= 0.0 && s2 >= 0.0) || len_diffCenter <= Radius) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}


/** \brief marks cells that lie in a sliced cone-type patch
 *
 * \author Rodrigo Miguez (rodrigo) rodrigo.miguez@rwth-aachen.de
 * \date 12.05.2022
 *
 * This algorithm does the following:
 *
 *  a. Precalculate the required values such as the sliced cone axis, circle normals according to
 *     the type of frustum cone being used (aligned or not aligned circle normals).
 *     For the non-aligned case, a check of the provided normal vectors is done.
 *
 *  b. A dry-run is performed to mark the cells that must be refined in the current level and the
 *     same is done for the halo cells if we run in parallel mode.
 *
 * \param[in] level the level to run over
 * \param[in] patchStr the type of frustum cone patch being used
 *
 */
template <MInt nDim>
void GridgenPar<nDim>::markLocalSlicedCone(MInt compRfnLvl, MInt patch, MInt solver, MString patchStr) {
  TRACE();

  if(patchStr != "A" && patchStr != "N") {
    TERMM(1, "This function doesn't support the selected patch type!");
  }

  const MBool isAligned = (patchStr == "A") ? true : false;

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local sliced cone for solver " << solver << endl;
  // Get coordinates
  std::array<MFloat, nDim> leftCoord;
  std::array<MFloat, nDim> rightCoord;
  for(MInt d = 0; d < nDim; d++) {
    leftCoord[d] = bp->localRfnPatchProperties[pos][d];
    rightCoord[d] = bp->localRfnPatchProperties[pos][d + nDim];
  }

  m_log << "         - left coord: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << leftCoord[d] << " ";
  }
  m_log << endl << "         - right coord: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << rightCoord[d] << " ";
  }

  // Get radius
  const MInt col = (isAligned) ? 2 : 4;
  const MFloat leftR = bp->localRfnPatchProperties[pos][col * nDim];
  const MFloat rightR = bp->localRfnPatchProperties[pos][col * nDim + 1];

  m_log << endl << "         - left radius: " << leftR << endl;
  m_log << "         - right radius: " << rightR << endl;

  // a. precalculate the required values
  std::array<MFloat, nDim> normalDifforigAB;
  MFloat length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normalDifforigAB[d] = rightCoord[d] - leftCoord[d];
    length += normalDifforigAB[d] * normalDifforigAB[d];
  }
  length = sqrt(length);

  // Get unit normals
  std::array<MFloat, nDim> leftNormal;
  std::array<MFloat, nDim> rightNormal;
  for(MInt d = 0, i = 0; d < nDim; d++, i++) {
    if(isAligned) {
      rightNormal[d] = normalDifforigAB[d] / length;
      leftNormal[d] = -rightNormal[d];
    } else {
      leftNormal[d] = bp->localRfnPatchProperties[pos][2 * nDim + d];
      rightNormal[d] = bp->localRfnPatchProperties[pos][2 * nDim + d + nDim];
    }
  }

  m_log << "         - left normal: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << leftNormal[d] << " ";
  }
  m_log << endl << "         - right normal: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << rightNormal[d] << " ";
  }
  m_log << endl;

  // Check if the direction of the circle normals are valid (only necessary when using non
  // axis-aligned normals)
  if(!isAligned) {
    MFloat leftNormalDir = 0.0;
    MFloat rightNormalDir = 0.0;
    MBool isAborting = false;
    stringstream ss;
    for(MInt d = 0; d < nDim; d++) {
      leftNormalDir += (rightCoord[d] - leftCoord[d]) * leftNormal[d];
      rightNormalDir += (leftCoord[d] - rightCoord[d]) * rightNormal[d];
      // Normals are at a plane not crossed by the cone axis
      if(approx(normalDifforigAB[d], 0.0, MFloatEps) && !approx(leftNormal[d], 0.0, MFloatEps)) {
        ss << "Warning: Normal values for the left circle are not valid at dir = " << d << "!" << endl;
        isAborting = true;
      } else if(approx(normalDifforigAB[d], 0.0, MFloatEps) && !approx(rightNormal[d], 0.0, MFloatEps)) {
        ss << "Warning: Normal values for the right circle are not valid at dir = " << d << "!" << endl;
        isAborting = true;
      }
      // Normals are pointing in the same direction
      if((d == nDim - 1) && leftNormalDir >= 0.0 && rightNormalDir >= 0.0) {
        ss << "Warning: Normals are not consistent!" << endl;
        isAborting = true;
      }
    }
    if(isAborting) {
      cerr << ss.str() << "Aborting...";
      m_log << ss.str() << "Aborting...";
      TERMM(1, ss.str() + "Aborting...");
    }
  }

  // b. do a dry run so that we know how many cells we can expect,
  //   also mark these cells
  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(a_isInSolver(k, solver)) {
      if(isInsideSlicedCone(&a_coordinate(k, 0), leftCoord.data(), rightCoord.data(), leftNormal.data(),
                            rightNormal.data(), normalDifforigAB.data(), leftR, rightR)) {
        a_isToRefineForSolver(k, solver) = 1;
      }
    }
  }

  // b. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);

    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(a_isInSolver(k, solver)) {
          if(isInsideSlicedCone(&a_coordinate(k, 0), leftCoord.data(), rightCoord.data(), leftNormal.data(),
                                rightNormal.data(), normalDifforigAB.data(), leftR, rightR)) {
            a_isToRefineForSolver(k, solver) = 1;
          }
        }
      }
    }
  }
}


/** \brief checks if point (cell) is inside a sliced cone
 *
 * \author Rodrigo Miguez (rodrigo) rodrigo.miguez@rwth-aachen.de
 * \date 25.05.2022
 *
 * This algorithm does the following:
 *
 *  a. For any given point we find the perpendicular projection of vector p (the point we are
 *     testing):
 *     lambda = \frac{(\vec{p} - \vec{a}) \cdot (\vec{b} - \vec{a})}
 *                    {(\vec{b} - \vec{a}) \cdot (\vec{b} - \vec{a})}
 *     where, \vec{a} and/or \vec{b}: are the center points of the circles that form the sliced cone
 *
 *  b. Calculate the point of projection of p on the line:
 *     \vec{r_p} = \vec{a} + lambda * (\vec{b} - \vec{a})
 *
 *  c. Calculate the distance of the point to its projection on the line:
 *     lenght = |\vec{p} - \vec{r_p}|
 *
 *  d. Check if point is inside the sliced cone radius:
 *     length <= Ra + lambda * (Rb - Ra), where Ra/Rb: are the radius of the circles that form the
 *                                                     sliced cone
 *
 *  e. Finally, check if the point lies between the circles that form the sliced cone by
 *     calculating the normals of the circles and comparing if the point is in the opposite
 *     direction of both normals
 *
 * \param[in] pointCoord cell coordinate
 * \param[in] (left/right)Coord coordinates of the left and right circles
 * \param[in] (left/right)Normal normals of the left and right circles
 * \param[in] normalDifforigAB sliced cone axis
 * \param[in] (left/right)R radius of the left and right circles
 *
 */
template <MInt nDim>
MBool GridgenPar<nDim>::isInsideSlicedCone(const MFloat* const pointCoord, const MFloat* const leftCoord,
                                           const MFloat* const rightCoord, const MFloat* const leftNormal,
                                           const MFloat* const rightNormal, const MFloat* const normalDifforigAB,
                                           const MFloat leftR, const MFloat rightR) {
  TRACE();

  // Calculate lambda
  MFloat lambda = 0.0;
  MFloat temp0 = 0.0;
  MFloat temp1 = 0.0;
  // dot product
  for(MInt d = 0; d < nDim; d++) {
    temp0 += (pointCoord[d] - leftCoord[d]) * normalDifforigAB[d];
    temp1 += normalDifforigAB[d] * normalDifforigAB[d];
  }
  lambda = temp0 / temp1;

  // Store the length between pointCoord and projectionCoord
  MFloat lengthPP = 0.0;

  std::array<MFloat, nDim> projectionCoord;

  for(MInt d = 0; d < nDim; d++) {
    projectionCoord[d] = leftCoord[d] + lambda * normalDifforigAB[d];
    lengthPP += pow(pointCoord[d] - projectionCoord[d], 2);
  }
  lengthPP = sqrt(lengthPP);
  if(lengthPP <= (leftR + lambda * (rightR - leftR))) {
    MFloat isPointDirLeftNormal = 0.0;
    MFloat isPointDirRightNormal = 0.0;
    // dot product
    for(MInt d = 0; d < nDim; d++) {
      isPointDirLeftNormal += leftNormal[d] * (pointCoord[d] - leftCoord[d]);
      isPointDirRightNormal += rightNormal[d] * (pointCoord[d] - rightCoord[d]);
    }
    // Only mark if pointCoord is in the oposite direction of both circle normals
    if(isPointDirLeftNormal <= 0.0 && isPointDirRightNormal <= 0.0) {
      return true;
    }
  }

  return false;
}


/** \brief marks cells that lie in a cone-type patch with a smooth hat
 *
 * \author Andreas Lintermann
 * \date 12.12.2013
 *
 * This algorithm does the following:
 *
 *  a. It performs a dry run of the refinement, i.e., runs over alls
 *     cells on the current level, counts the number of cells to refine
 *     (the cells which are inside a given cone and sphere) and marks
 *     these cells (b_properties[5] = 1).
 *     First the center of the sphere and the new center point of the
 *     cone are calculated. Then, the distance of a coordinates of a
 *     given cell is calculated to obtain the enclosed angle between
 *     this the vector of the coordinates and the first provided point.
 *     In case the angle is larger than the one provided in the
 *     properties, this cell is marked as outside, otherwise as inside.
 *     The distance of the coordinates to the line is calculated
 *     as described in markLocalCylinder.
 *  b. The same is done for the halo cells if we run in parallel mode
 *
 * \param[in] level the level to run over
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markLocalHat(MInt compRfnLvl, MInt patch, MInt solver) {
  TRACE();

  SolverRefinement* bp = m_solverRefinement + solver;

  const MInt pos = bp->localRfnLevelPropertiesOffset[compRfnLvl] + patch;
  const MInt gridLvl = compRfnLvl + bp->maxUniformRefinementLevel;

  m_log << "       * marking local hat for solver " << solver << endl;
  m_log << "         - left point: ";
  for(MInt d = 0; d < nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << "         - right point: ";
  for(MInt d = nDim; d < 2 * nDim; d++) {
    m_log << bp->localRfnPatchProperties[pos][d] << " ";
  }
  m_log << endl;
  m_log << "         - opening angle: " << bp->localRfnPatchProperties[pos][2 * nDim] << endl;
  m_log << "         - radius: " << bp->localRfnPatchProperties[pos][2 * nDim + 1] << endl;
  m_log << "         - thickness: " << bp->localRfnPatchProperties[pos][2 * nDim + 2] << endl;


  // SolverRefinement* bp = m_solverRefinement + 0;

  MFloat Radius = bp->localRfnPatchProperties[pos][2 * nDim + 1];
  MFloat thickness = bp->localRfnPatchProperties[pos][2 * nDim + 2];
  MFloat angle = bp->localRfnPatchProperties[pos][2 * nDim];
  MFloat angle_rad = angle * PI / 180.0;
  MFloat deltaA = thickness / (F2 * sin(angle_rad));
  // original property
  std::array<MFloat, nDim> normalDifforigAB;
  MFloat length = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normalDifforigAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - bp->localRfnPatchProperties[pos][d];
    length += normalDifforigAB[d] * normalDifforigAB[d];
  }

  length = sqrt(length);

  for(MInt d = 0; d < nDim; d++) {
    normalDifforigAB[d] /= length;
  }

  std::array<MFloat, nDim> newPosA;
  std::array<MFloat, nDim> sphereCenter;
  std::array<MFloat, nDim> normalDiffAB;
  std::array<MFloat, nDim> normalDiffBA;
  MBoolScratchSpace insideLargerCone((m_levelOffsets[gridLvl][1] - m_levelOffsets[gridLvl][0]), AT_,
                                     "insideLargerCone");

  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    insideLargerCone(k - m_levelOffsets[gridLvl][0]) = false;
  }

  for(MInt dir = 1; dir > -2; dir = dir - 2) {
    // a. precalculate the required vector for the inner/outer cone
    MFloat radius = (Radius + dir * thickness / F2) * cos(angle_rad);
    // calucate the new position of A
    MFloat lengthAB = 0.0;
    for(MInt d = 0; d < nDim; d++) {
      // splitpoint between cone and sphere
      newPosA[d] = (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d])
                   + radius / tan(angle_rad) * normalDifforigAB[d];

      sphereCenter[d] = (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d])
                        + radius * (1.0 / tan(angle_rad) + tan(angle_rad)) * normalDifforigAB[d];

      normalDiffAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - newPosA[d];
      normalDiffBA[d] = newPosA[d] - bp->localRfnPatchProperties[pos][d + nDim];

      lengthAB += normalDiffAB[d] * normalDiffAB[d];
    }

    lengthAB = sqrt(lengthAB);


    for(MInt d = 0; d < nDim; d++) {
      normalDiffAB[d] /= lengthAB;
    }

    // b. do a dry run so that we know how many cells we can expect,
    //   also mark these cells
    // MInt rfnCount = 0;
    for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
      if(a_isInSolver(k, solver)) {
        MFloat* coords = &a_coordinate(k, 0);

        std::array<MFloat, nDim> diffA;
        std::array<MFloat, nDim> difforigA;
        std::array<MFloat, nDim> diffB;
        std::array<MFloat, nDim> diffCenter;

        MFloat s1 = 0.0;
        MFloat s2 = 0.0;
        MFloat len_diffOrig = 0.0;
        MFloat len_diffCenter = 0.0;
        for(MInt d = 0; d < nDim; d++) {
          diffA[d] = coords[d] - newPosA[d];
          diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
          difforigA[d] = coords[d] - (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d]);
          diffCenter[d] = coords[d] - sphereCenter[d];

          s1 += diffA[d] * normalDiffAB[d];
          s2 += diffB[d] * normalDiffBA[d];
          len_diffOrig += difforigA[d] * difforigA[d];
          len_diffCenter += diffCenter[d] * diffCenter[d];
        }

        len_diffOrig = sqrt(len_diffOrig);
        len_diffCenter = sqrt(len_diffCenter);

        std::array<MFloat, nDim> cross;
        MFloat distance = 0.0;
        for(MInt d = 0; d < nDim; d++) {
          MInt p1 = (d + 1) % (nDim);
          MInt p2 = (d + 2) % (nDim);

          cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
          distance += cross[d] * cross[d];
        }

        distance = sqrt(distance);

        MFloat ang = asin(distance / len_diffOrig) * 180 / PI;

        if((ang <= angle && s1 >= 0.0 && s2 >= 0.0) || len_diffCenter <= Radius + dir * thickness / F2) {
          insideLargerCone(k - m_levelOffsets[gridLvl][0]) = true;
          if(dir == -1) {
            insideLargerCone(k - m_levelOffsets[gridLvl][0]) = false;
          }
        }
      }
    }
  }

  for(MInt k = m_levelOffsets[gridLvl][0]; k < m_levelOffsets[gridLvl][1]; k++) {
    if(insideLargerCone(k - m_levelOffsets[gridLvl][0])) {
      a_isToRefineForSolver(k, solver) = 1;
    }
  }


  // b. do the same for the halos
  if(noDomains() > 1) {
    MInt lev_pos = 2 * (gridLvl - m_minLevel);
    MBoolScratchSpace insideLargerConeHalo(
        (m_haloCellOffsets[m_noNeighborDomains - 1][lev_pos + 1] - m_haloCellOffsets[0][lev_pos]), AT_,
        "insideLargerCone");
    for(MInt k = m_haloCellOffsets[0][lev_pos]; k < m_haloCellOffsets[m_noNeighborDomains - 1][lev_pos + 1]; k++) {
      insideLargerConeHalo(k - m_haloCellOffsets[0][lev_pos]) = false;
    }

    for(MInt dir = 1; dir > -2; dir = dir - 2) {
      // a. precalculate the required vector for the inner/outer cone
      MFloat radius = (Radius + dir * thickness / F2) * cos(angle_rad);
      // calucate the new position of A
      MFloat lengthAB = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        // splitpoint between cone and sphere
        newPosA[d] = (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d])
                     + radius / tan(angle_rad) * normalDifforigAB[d];
        sphereCenter[d] = (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d])
                          + radius * (1.0 / tan(angle_rad) + tan(angle_rad)) * normalDifforigAB[d];
        normalDiffAB[d] = bp->localRfnPatchProperties[pos][d + nDim] - newPosA[d];
        normalDiffBA[d] = newPosA[d] - bp->localRfnPatchProperties[pos][d + nDim];

        lengthAB += normalDiffAB[d] * normalDiffAB[d];
      }

      lengthAB = sqrt(lengthAB);


      for(MInt d = 0; d < nDim; d++)
        normalDiffAB[d] /= lengthAB;

      for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
        for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
          if(a_isInSolver(k, solver)) {
            MFloat* coords = &a_coordinate(k, 0);

            std::array<MFloat, nDim> diffA;
            std::array<MFloat, nDim> difforigA;
            std::array<MFloat, nDim> diffB;
            std::array<MFloat, nDim> diffCenter;

            MFloat s1 = 0.0;
            MFloat s2 = 0.0;
            MFloat len_diffOrig = 0.0;
            MFloat len_diffCenter = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              diffA[d] = coords[d] - newPosA[d];
              diffB[d] = coords[d] - bp->localRfnPatchProperties[pos][d + nDim];
              difforigA[d] = coords[d] - (bp->localRfnPatchProperties[pos][d] - dir * deltaA * normalDifforigAB[d]);
              diffCenter[d] = coords[d] - sphereCenter[d];

              s1 += diffA[d] * normalDiffAB[d];
              s2 += diffB[d] * normalDiffBA[d];
              len_diffOrig += difforigA[d] * difforigA[d];
              len_diffCenter += diffCenter[d] * diffCenter[d];
            }

            len_diffOrig = sqrt(len_diffOrig);
            len_diffCenter = sqrt(len_diffCenter);

            std::array<MFloat, nDim> cross;
            MFloat distance = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              MInt p1 = (d + 1) % (nDim);
              MInt p2 = (d + 2) % (nDim);

              cross[d] = diffA[p1] * normalDiffAB[p2] - normalDiffAB[p1] * diffA[p2];
              distance += cross[d] * cross[d];
            }

            distance = sqrt(distance);

            MFloat ang = asin(distance / len_diffOrig) * 180 / PI;

            if((ang <= angle && s1 >= 0.0 && s2 >= 0.0) || len_diffCenter <= Radius + dir * thickness / F2) {
              insideLargerConeHalo(k - m_haloCellOffsets[0][lev_pos]) = true;
              if(dir == -1) {
                insideLargerConeHalo(k - m_haloCellOffsets[0][lev_pos]) = false;
              }
            }
          }
        }
      } // for dom k
    }   // for dir

    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt k = m_haloCellOffsets[dom][lev_pos]; k < m_haloCellOffsets[dom][lev_pos + 1]; k++) {
        if(insideLargerConeHalo(k - m_haloCellOffsets[0][lev_pos])) {
          a_isToRefineForSolver(k, solver) = 1;
        }
      }
    }
  } // if domain id
}

/** \brief marks cells that lie in a certain distance to the boundary
 * \author Andreas Lintermann
 * \date 02.06.2014
 *
 * \param[in] level_ the level to run over
 * \param[in] dists the distances to use
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::markBndDistance(MInt level_, MInt solver) {
  TRACE();

  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    if(a_refinementDistance(i) < numeric_limits<MInt>::max()) {
      a_isToRefineForSolver(i, solver) = 1;
    }
  }

  // parallel part
  if(noDomains() > 1) {
    // 6.2b.2.3 dry run for the halo cells, how many new cells will we generate?
    MInt lev_pos = 2 * (level_ - m_minLevel);
    for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
      for(MInt i = m_haloCellOffsets[dom][lev_pos]; i < m_haloCellOffsets[dom][lev_pos + 1]; i++) {
        if(a_refinementDistance(i) < numeric_limits<MInt>::max()) {
          a_isToRefineForSolver(i, solver) = 1;
        }
      }
    }
  }
}

/** \brief checks if a dynamic load balancing is needed.
 *
 * \author Jerry Grimmen
 * \date 02.04.2014
 *
 * This funtions checks the current state of the grid generation.
 *  1. Get the overall number of cells in the whole grid.
 *  2. Calculate the average number of cells and allow a small overhead.
 *  3. Checks if any domain has more than the average.
 *  4. If it's the case, do a load balance.
 *
 */

template <MInt nDim>
void GridgenPar<nDim>::checkLoadBalance(MInt in_level) {
  TRACE();

  // 0. Debug: Use to skip spezified levels.
  if(in_level == 42) {
    return;
  }

  // 1. Getting the global number of cells in the grid.
  MLong globalNoCells = 0;
  MLong localNoCells = (MLong)m_noCells;
  MPI_Allreduce(&localNoCells, &globalNoCells, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "localNoCells", "globalNoCells");

  // 2. Calculate the average number of cells each domain should have and adds a small amount to it.
  MLong averageNoCells = globalNoCells / ((MLong)noDomains());
  // averageNoCells *= 1.5; // Should be fine tuned!
  averageNoCells += (averageNoCells / 2); // Should be fine tuned!

  // 3. Looks if any domain has more cells than the average.
  MInt needsLoadBalancing = 0;
  if(m_noCells > averageNoCells) {
    needsLoadBalancing = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, &needsLoadBalancing, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                "needsLoadBalancing");

  // 4. If any domain has more cells than to the load balancing.
  if(needsLoadBalancing > 0) {
    m_log << "Imbalance found, starting load balancing." << endl;
    dynamicLoadBalancing();
    m_hasBeenLoadBalanced = true;

// Debug: Neighborhood complete check activated only in extra debug mode.
#ifdef MAIA_EXTRA_DEBUG
    for(MInt l = m_minLevel; l < in_level + 1; ++l) {
      m_log << "Checking Neighborhood for level " << l << endl;
      outStream << "Checking Neighborhood for level " << l << endl;
      checkNeighborhood(l);
      checkNeighborhoodDistance(l);
      checkNeighborhoodIntegrity(l);
    }
#endif

    m_log << "****************************************" << endl;
    m_log << "*                                      *" << endl;
    m_log << "*   Finished, the load balancing is    *" << endl;
    m_log << "*  Grid generation, with, we continue  *" << endl;
    m_log << "*                                      *" << endl;
    m_log << "*      May the force be with you       *" << endl;
    m_log << "*                                      *" << endl;
    m_log << "****************************************" << endl;

    //    m_log << "Load balancing finished, continuing with grid generation." << endl;

  } else {
    m_log << "No imbalance in load was found, continuing with grid generation." << endl;
  }

  return;
}

/** \brief Rebalances the grid dynamicaly while grid building.
 *
 * \author Jerry Grimmen
 * \date 26.02.2014
 *
 * This function can be splitted in the following big parts:
 *  - Gathering informations:
 *  -- Calculate current offsets.
 *  -- Calculate new global ids for local cells.
 *  -- Calculate new offsets.
 *  -- Calculate the domain to which the local cells will belong.
 *  -- Calculate the number of cells we will send and receive.
 *
 *  - Updating halo ids:
 *  -- Send and receive new halo global ids.
 *  -- Changes all local ids from local cell to global ids.
 *
 *  - Sending and receiving new cells:
 *  -- Send and receive all cell datas into scratchspaces.
 *
 *  - Recreation of the cell array:
 *  -- Moves all cells we keep to the end of the cell array.
 *  -- In domain and level order, add the cells back to the cell array.
 *
 *  - Regenerate halo cells:
 *  -- Gather the halo ids of all the parents of the first cell in the domain.
 *  -- Tags window cells and collect the halo ids.
 *  -- Communicate halo global ids to the other domains.
 *  -- Send and receive all halo cell datas into scratchspaces.
 *  -- Put the new halo cells into the cell array.
 *
 *  - Global to local:
 *  -- Create a global to local map.
 *  -- Add the local cells to the map.
 *  -- Add the halo cells to the map.
 *  -- Changes all global ids of all cells to local ids.
 */

template <MInt nDim>
void GridgenPar<nDim>::dynamicLoadBalancing() {
  TRACE();

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 0 - Setting timer" << endl;
#endif

  // Setting timer
  MInt t_dynamicLoadBalancing = 0;
  NEW_SUB_TIMER_STATIC(t_dynLoadBal, "dynamic load balancing", m_t_createComputationalGrid);
  t_dynamicLoadBalancing = t_dynLoadBal;

  MInt t_gatherInformation = 0;
  MInt t_updateHaloIds = 0;
  MInt t_communicateCells = 0;
  MInt t_recreateCells = 0;
  MInt t_createHaloCells = 0;
  MInt t_globalToLocal = 0;

  NEW_SUB_TIMER_STATIC(t_dynLoadBal1, "Gathering informations", t_dynamicLoadBalancing);
  NEW_SUB_TIMER_STATIC(t_dynLoadBal2, "Update halo cells ids", t_dynamicLoadBalancing);
  NEW_SUB_TIMER_STATIC(t_dynLoadBal3, "Sending and receiving new cells", t_dynamicLoadBalancing);
  NEW_SUB_TIMER_STATIC(t_dynLoadBal4, "Recreation of cell array", t_dynamicLoadBalancing);
  NEW_SUB_TIMER_STATIC(t_dynLoadBal5, "Regenerate halo cells", t_dynamicLoadBalancing);
  NEW_SUB_TIMER_STATIC(t_dynLoadBal6, "Global to local", t_dynamicLoadBalancing);

  t_gatherInformation = t_dynLoadBal1;
  t_updateHaloIds = t_dynLoadBal2;
  t_communicateCells = t_dynLoadBal3;
  t_recreateCells = t_dynLoadBal4;
  t_createHaloCells = t_dynLoadBal5;
  t_globalToLocal = t_dynLoadBal6;

  RECORD_TIMER_START(t_dynamicLoadBalancing);

  // ==============================
  // #=- Gathering informations -=#
  // ==============================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 1 - Gathering informations" << endl;
#endif

  RECORD_TIMER_START(t_gatherInformation);

  // 1. Get current number of cells on each domain.
  MLong noCells = (MLong)m_noCells;
  MLongScratchSpace globalNoCellsPerDomain(noDomains() + 1, AT_, "globalNoCellsPerDomain");
  MPI_Allgather(&noCells, 1, MPI_LONG, globalNoCellsPerDomain.getPointer(), 1, MPI_LONG, mpiComm(), AT_, "noCells",
                "globalNoCellsPerDomain.getPointer()");

  // 2. Calculate the current offset.
  MLongScratchSpace globalOffset(noDomains() + 1, AT_, "globalOffset");
  MLong globalNoCells = globalNoCellsPerDomain[0];
  globalOffset[0] = 0;

  for(MInt i = 1; i < noDomains(); ++i) {
    globalOffset[i] = globalOffset[i - 1] + globalNoCellsPerDomain[i - 1];
    globalNoCells += globalNoCellsPerDomain[i];
  }

  /* 3. Now we create an information array localCellInfo.
   * [cellId][0] : local id of the cell in its own domain.
   * [cellId][1] : number of offsprings of the cell.
   * [cellId][2] : global id of the cell.
   * [cellId][3] : id of the domain to which the cell will belong.
   */

  MLongScratchSpace localCellInfo(m_noCells, 4, AT_, "localCellInfo");
  MLong currentGlobalId = globalOffset[globalDomainId()] - 1;

  // 3.1. Collecting information [0] and [1] (local id and number of offsprings)
  MInt currentCellId = 0;

  // 3.1.1. Gather all cells which are before the cell at m_levelOffsets[m_minLevel][0] (may happen after load balacing)
  MInt noMissingParents = m_noMissingParents;
#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Gathering informations:: m_noMisingParents is " << m_noMissingParents
        << endl;
#endif

  // Short Fix to match the new traverseDFGlobalId
  MFloatScratchSpace minWorkloadPerCell(m_noCells, AT_, "minWorkloadPerCell");
  minWorkloadPerCell.fill(0.0);
  MFloatScratchSpace minWeightList(m_noCells, AT_, "minWeightList");
  minWeightList.fill(0.0);
  MFloatScratchSpace minWorkloadList(m_noCells, AT_, "minWorkloadList");
  minWorkloadList.fill(0.0);
  MFloat currentWorkload = 0.0;

  while(noMissingParents > 0) {
    MInt k = 0;
    MInt cellId = m_levelOffsets[m_minLevel + noMissingParents][0];
    while((a_hasProperty((MInt)a_parentId(cellId), 4) == 1) && (a_level(cellId) == (m_minLevel + noMissingParents))) {
      currentGlobalId++;
      a_globalId(cellId) = currentGlobalId;
      localCellInfo(currentCellId, 0) = (MLong)cellId;

      MInt last = currentCellId;

      traverseDFGlobalId(cellId, &currentGlobalId, localCellInfo, minWorkloadPerCell, &currentWorkload, minWeightList,
                         minWorkloadList, &currentCellId);

      MLong noOffsprings = 0;

      if(cellId == 0)
        noOffsprings = currentGlobalId - globalOffset[globalDomainId()] + 1;
      else
        noOffsprings = currentGlobalId - a_globalId(cellId);

      localCellInfo(last, 1) = noOffsprings;

      ++currentCellId;
      ++k;
      cellId = m_levelOffsets[m_minLevel + noMissingParents][0] + k;

    } // end of : while ( (a_hasProperty(a_parentId(cellId), 4) == 1) && (a_level(cellId) == (m_minLevel +
      // noMissingParents)) )

    --noMissingParents;

  } // end of : while ( noMissingParents > 0 )

  // 3.1.2. Gather all cells from m_levelOffsets[m_minLevel][0] and after.
  for(MInt cellId = m_levelOffsets[m_minLevel][0]; cellId < m_levelOffsets[m_minLevel][1]; cellId++, currentCellId++) {
    ++currentGlobalId;
    a_globalId(cellId) = currentGlobalId;
    localCellInfo(currentCellId, 0) = (MLong)cellId;
    MInt lastCellId = currentCellId;

    traverseDFGlobalId(cellId, &currentGlobalId, localCellInfo, minWorkloadPerCell, &currentWorkload, minWeightList,
                       minWorkloadList, &currentCellId);

    MLong noOffsprings = 0;
    if(cellId == 0) {
      noOffsprings = currentGlobalId - globalOffset[globalDomainId()] + 1;
    } else {
      noOffsprings = currentGlobalId - a_globalId(cellId) + 1;
    }

    localCellInfo(lastCellId, 1) = noOffsprings;
  } // end of : for (MInt cellId = m_levelOffsets[m_minLevel][0]; cellId < m_levelOffsets[m_minLevel][1]; cellId++,
    // currentCellId++)

  // 3.2. Collecting information [2] (global id)
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    localCellInfo(cellId, 2) = a_globalId((MInt)localCellInfo(cellId, 0));
  }

  // 3.3. Calculate new fixed offsets such that on all domains there is no more than 1 cell difference.
  globalOffset.fill(-1);
  globalNoCellsPerDomain.fill(0);

  globalOffset[0] = 0;
  for(MLong d = 1; d < noDomains(); ++d) {
    globalOffset[d] = globalOffset[d - 1] + ((globalNoCells - globalOffset[d - 1]) / (((MLong)noDomains()) - (d - 1)));
    globalNoCellsPerDomain[d - 1] = globalOffset[d] - globalOffset[d - 1];
  }
  globalNoCellsPerDomain[noDomains() - 1] = globalNoCells - globalOffset[noDomains() - 1];
  globalNoCellsPerDomain[noDomains()] = globalNoCells;
  globalOffset[noDomains()] = globalNoCells;

  // 3.4. Collecting information [3] (id of the new domain the cell will belong)
  MInt currentCpu = 0;
  for(MInt i = 1; i < (noDomains() + 1); ++i) {
    if(globalOffset[i] > localCellInfo(0, 2)) {
      currentCpu = i;
      break;
    }
  }

  for(MInt i = 0; i < m_noCells; ++i) {
    if(localCellInfo(i, 2) < globalOffset[currentCpu]) {
      localCellInfo(i, 3) = (MLong)(currentCpu - 1);
      continue;
    }

    if(localCellInfo(i, 2) == globalOffset[currentCpu]) {
      localCellInfo(i, 3) = (MLong)currentCpu;
      ++currentCpu;
      continue;
    }

  } // end of : for (MInt i = 0; i < m_noCells; ++i)

  // 4. Calculating the number of cells to send and to receive.
  MIntScratchSpace noCellsToReceive(noDomains() + 1, AT_, "noCellsToReceive");
  MIntScratchSpace noCellsToSend(noDomains() + 1, AT_, "noCellsToSend");
  noCellsToReceive.fill(0);
  noCellsToSend.fill(0);

  for(MInt i = 0; i < m_noCells; ++i) {
    if(localCellInfo(i, 3) == (MLong)globalDomainId()) {
      continue;
    }

    noCellsToSend[localCellInfo(i, 3)] += 1;
  }

  for(MInt i = 0; i < noDomains(); ++i) {
    MPI_Scatter(noCellsToSend.getPointer(), 1, MPI_INT, &noCellsToReceive[i], 1, MPI_INT, i, mpiComm(), AT_,
                "noCellsToSend.getPointer()", "noCellsToReceive[i]");
  }

  for(MInt i = 0; i < noDomains(); ++i) {
    noCellsToSend[noDomains()] += noCellsToSend[i];
    noCellsToReceive[noDomains()] += noCellsToReceive[i];
#ifdef MAIA_EXTRA_DEBUG
    m_log << "DynamicLoadBalancing Debug: To domain " << i << ", we will send " << noCellsToSend[i] << " and receive "
          << noCellsToReceive[i] << " cells." << endl;
#endif
  }

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Debug: Total number of cells, we will send " << noCellsToSend[noDomains()]
        << " and receive " << noCellsToReceive[noDomains()] << " cells." << endl;
#endif

  RECORD_TIMER_STOP(t_gatherInformation);

  // ================================
  // #=- Updating halo global ids -=#
  // ================================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 2 - Update halo global ids" << endl;
#endif

  RECORD_TIMER_START(t_updateHaloIds);

  communicateHaloGlobalIds(-42);

  RECORD_TIMER_STOP(t_updateHaloIds);

  // =======================================
  // #=- Sending and receiving new cells -=#
  // =======================================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 3 - Sending and receiving new cells" << endl;
#endif

  RECORD_TIMER_START(t_communicateCells);

  { // extra '{}' span for sending / receiving the new cells.
    MLongScratchSpace newParentId(noCellsToReceive[noDomains()], AT_, "newParentId");
    MLongScratchSpace newGlobalId(noCellsToReceive[noDomains()], AT_, "newGlobalId");
    MIntScratchSpace newLevel(noCellsToReceive[noDomains()], AT_, "newLevel");
    MIntScratchSpace newNoChildren(noCellsToReceive[noDomains()], AT_, "newNoChildren");

    MLongScratchSpace newChildId(noCellsToReceive[noDomains()], m_maxNoChildren, AT_, "newChildId");
    MLongScratchSpace newNeighborId(noCellsToReceive[noDomains()], m_noNeighbors, AT_, "newNeighborId");
    MIntScratchSpace newProperty(noCellsToReceive[noDomains()], 8, AT_, "newProperty");
    MIntScratchSpace newIsInSolver(noCellsToReceive[noDomains()], m_noSolvers, AT_, "newIsInSolver");
    MIntScratchSpace newIsSolverBoundary(noCellsToReceive[noDomains()], m_noSolvers, AT_, "newIsSolverBoundary");
    MIntScratchSpace newIsToRefineForSolver(noCellsToReceive[noDomains()], m_noSolvers, AT_, "newIsToRefineForSolver");
    MFloatScratchSpace newCoordinate(noCellsToReceive[noDomains()], nDim, AT_, "newCoordinate");

    newParentId.fill(-1);
    newGlobalId.fill(-1);
    newLevel.fill(-1);
    newNoChildren.fill(-1);
    newChildId.fill(-1);
    newNeighborId.fill(-1);
    newProperty.fill(-1);
    newCoordinate.fill(-1.0);

    { // sending / receiving parentId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_parentId((MInt)localCellInfo(i, 0));
            ++currentId;
          }
        }
      }

      communicateLong(newParentId, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving parentId

    { // sending / receiving globalId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_globalId((MInt)localCellInfo(i, 0));
            ++currentId;
          }
        }
      }

      communicateLong(newGlobalId, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving globalId

    { // sending / receiving level
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_level((MInt)localCellInfo(i, 0));
            ++currentId;
          }
        }
      }

      communicateInt(newLevel, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving level

    { // sending / receiving noChildren
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_noChildren((MInt)localCellInfo(i, 0));
            ++currentId;
          }
        }
      }

      communicateInt(newNoChildren, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving noChildren

    // sending / receiving childId
    for(MInt childId = 0; childId < m_maxNoChildren; ++childId) {
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MLongScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_childId((MInt)localCellInfo(i, 0), childId);
            ++currentId;
          }
        }
      }

      communicateLong(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newChildId(i, childId) = recvBuffer(i);
      }

    } // end of : sending / receiving childId

    // sending / receiving neighborId
    for(MInt neighborId = 0; neighborId < m_noNeighbors; ++neighborId) {
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MLongScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_neighborId((MInt)localCellInfo(i, 0), neighborId);
            ++currentId;
          }
        }
      }

      communicateLong(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newNeighborId(i, neighborId) = recvBuffer(i);
      }

    } // end of : sending / receiving neighborId

    // sending / receiving property
    for(MInt property = 0; property < 8; ++property) {
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_hasProperty((MInt)localCellInfo(i, 0), property);
            ++currentId;
          }
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newProperty(i, property) = recvBuffer(i);
      }

    } // end of : sending / receiving property

    // sending / receiving solver affiliation
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_isInSolver((MInt)localCellInfo(i, 0), solver);
            ++currentId;
          }
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsInSolver(i, solver) = recvBuffer(i);
      }
    }

    // sending / receiving solver boundary
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_isSolverBoundary((MInt)localCellInfo(i, 0), solver);
            ++currentId;
          }
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsSolverBoundary(i, solver) = recvBuffer(i);
      }
    }

    // sending / receiving solver is to refine
    for(MInt solver = 0; solver < m_noSolvers; solver++) {
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_isToRefineForSolver((MInt)localCellInfo(i, 0), solver);
            ++currentId;
          }
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsToRefineForSolver(i, solver) = recvBuffer(i);
      }
    }

    // sending / receiving coordinate
    for(MInt dim = 0; dim < nDim; ++dim) {
      MFloatScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MFloatScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      MInt currentId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(noCellsToSend[dom] == 0) {
          continue;
        }

        for(MInt i = 0; i < m_noCells; ++i) {
          if(localCellInfo(i, 3) == dom) {
            sendBuffer(currentId) = a_coordinate((MInt)localCellInfo(i, 0), dim);
            ++currentId;
          }
        }
      }

      communicateDouble(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newCoordinate(i, dim) = recvBuffer(i);
      }

    } // end of : sending / receiving coordinate

    RECORD_TIMER_STOP(t_communicateCells);

    // ====================================
    // #=- Recreation of the cell array -=#
    // ====================================

#ifdef MAIA_EXTRA_DEBUG
    m_log << "DynamicLoadBalancing Extra Debug: Step 4 - Recreation of the cell array" << endl;
#endif

    RECORD_TIMER_START(t_recreateCells);

    // 1. Move own cells out of the way (moves the cells that we will keep to the end of the cell array)
    currentCellId = m_maxNoCells - 2;
    for(MInt i = m_noCells - 1; i > -1; --i) {
      if(localCellInfo(i, 3) == globalDomainId()) {
        a_parentId(currentCellId) = a_parentId((MInt)localCellInfo(i, 0));
        a_globalId(currentCellId) = a_globalId((MInt)localCellInfo(i, 0));
        a_noChildren(currentCellId) = a_noChildren((MInt)localCellInfo(i, 0));
        a_level(currentCellId) = a_level((MInt)localCellInfo(i, 0));

        for(MInt child = 0; child < m_maxNoChildren; ++child) {
          a_childId(currentCellId, child) = a_childId((MInt)localCellInfo(i, 0), child);
        }

        for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
          a_neighborId(currentCellId, ngh) = a_neighborId((MInt)localCellInfo(i, 0), ngh);
        }

        for(MInt dim = 0; dim < nDim; ++dim) {
          a_coordinate(currentCellId, dim) = a_coordinate((MInt)localCellInfo(i, 0), dim);
        }

        for(MInt property = 0; property < 8; ++property) {
          a_hasProperty(currentCellId, property) = a_hasProperty((MInt)localCellInfo(i, 0), property);
        }

        for(MInt solver = 0; solver < 8; solver++) {
          a_isInSolver(currentCellId, solver) = a_isInSolver((MInt)localCellInfo(i, 0), solver);
          a_isSolverBoundary(currentCellId, solver) = a_isSolverBoundary((MInt)localCellInfo(i, 0), solver);
          a_isToRefineForSolver(currentCellId, solver) = a_isToRefineForSolver((MInt)localCellInfo(i, 0), solver);
          // not working in parallel yet anyway
          a_noSolidLayer(currentCellId, solver) = -1;
        }

        localCellInfo(i, 0) = currentCellId;
        --currentCellId;
      } // end of : if (localCellInfo(i, 3) == globalDomainId())
    }   // end of : for (MInt i = m_noCells - 1; i > -1; --i)

    // 2. Copy all the cells in the right order in the cell array.
    MInt currentLocalCellId = 0;
    for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
      m_levelOffsets[level][0] = currentLocalCellId;
      currentCellId = 0;
      for(MInt dom = 0; dom < noDomains(); ++dom) {
        if(dom == globalDomainId()) {
          for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
            if((a_level((MInt)localCellInfo(cellId, 0)) == level) && (localCellInfo(cellId, 3) == globalDomainId())) {
              a_parentId(currentLocalCellId) = a_parentId((MInt)localCellInfo(cellId, 0));
              a_globalId(currentLocalCellId) = a_globalId((MInt)localCellInfo(cellId, 0));
              a_noChildren(currentLocalCellId) = a_noChildren((MInt)localCellInfo(cellId, 0));
              a_level(currentLocalCellId) = a_level((MInt)localCellInfo(cellId, 0));

              for(MInt child = 0; child < m_maxNoChildren; ++child) {
                a_childId(currentLocalCellId, child) = a_childId((MInt)localCellInfo(cellId, 0), child);
              }

              for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
                a_neighborId(currentLocalCellId, ngh) = a_neighborId((MInt)localCellInfo(cellId, 0), ngh);
              }

              for(MInt dim = 0; dim < nDim; ++dim) {
                a_coordinate(currentLocalCellId, dim) = a_coordinate((MInt)localCellInfo(cellId, 0), dim);
              }

              for(MInt property = 0; property < 8; ++property) {
                a_hasProperty(currentLocalCellId, property) = a_hasProperty((MInt)localCellInfo(cellId, 0), property);
              }

              for(MInt solver = 0; solver < 8; solver++) {
                a_isInSolver(currentLocalCellId, solver) = a_isInSolver((MInt)localCellInfo(cellId, 0), solver);
                a_isSolverBoundary(currentLocalCellId, solver) =
                    a_isSolverBoundary((MInt)localCellInfo(cellId, 0), solver);
                a_isToRefineForSolver(currentLocalCellId, solver) =
                    a_isToRefineForSolver((MInt)localCellInfo(cellId, 0), solver);
                // not working in parallel yet anyway
                a_noSolidLayer(currentLocalCellId, solver) = -1;
              }

              ++currentLocalCellId;

            } // end of : if ( (a_level((MInt)localCellInfo(cellId, 0)) == level) && (localCellInfo(cellId, 3) ==
              // globalDomainId()) )

          } // end of : for (MInt cellId = 0; cellId < m_noCells; ++cellId)
          continue;

        } // end of : if (dom == globalDomainId())

        if(noCellsToReceive[dom] == 0) {
          continue;
        } // end of : if (noCellsToReceive[dom] == 0)

        if(noCellsToReceive[dom] > 0) {
          for(MInt cellId = currentCellId; cellId < (currentCellId + noCellsToReceive[dom]); ++cellId) {
            if(newLevel(cellId) == level) {
              a_parentId(currentLocalCellId) = newParentId(cellId);
              a_globalId(currentLocalCellId) = newGlobalId(cellId);
              a_noChildren(currentLocalCellId) = newNoChildren(cellId);
              a_level(currentLocalCellId) = newLevel(cellId);

              for(MInt child = 0; child < m_maxNoChildren; ++child) {
                a_childId(currentLocalCellId, child) = newChildId(cellId, child);
              }

              for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
                a_neighborId(currentLocalCellId, ngh) = newNeighborId(cellId, ngh);
              }

              for(MInt dim = 0; dim < nDim; ++dim) {
                a_coordinate(currentLocalCellId, dim) = newCoordinate(cellId, dim);
              }

              for(MInt property = 0; property < 8; ++property) {
                a_hasProperty(currentLocalCellId, property) = (MBool)newProperty(cellId, property);
              }

              for(MInt solver = 0; solver < m_noSolvers; solver++) {
                a_isInSolver(currentLocalCellId, solver) = (MBool)newIsInSolver(cellId, solver);
                a_isSolverBoundary(currentLocalCellId, solver) = (MBool)newIsSolverBoundary(cellId, solver);
                a_isToRefineForSolver(currentLocalCellId, solver) = (MBool)newIsToRefineForSolver(cellId, solver);
                // solid-data is not communicated as its not working for
                // multiple-ranks yet!
                a_noSolidLayer(currentLocalCellId, solver) = -1;
              }

              // who knows what is there
              for(MInt solver = m_noSolvers; solver < 8; solver++) {
                a_isInSolver(currentLocalCellId, solver) = 0;
                a_isSolverBoundary(currentLocalCellId, solver) = 0;
                a_isToRefineForSolver(currentLocalCellId, solver) = 0;
                a_noSolidLayer(currentLocalCellId, solver) = -1;
              }

              ++currentLocalCellId;
            } // end of : if (cells->a[cellId].m_level[0] == level)
          }   // end of : for (MInt cellId = currentCellId; cellId < (currentCellId + cellToReceive[dom]); ++cellId)
          currentCellId += noCellsToReceive[dom];
        } // end of : if (cellToReceive[dom] > 0)


      } // end of :  for (MInt dom = 0; dom < noDomains(); ++dom)
      m_levelOffsets[level][1] = currentLocalCellId;
    } // end of : for (MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level)

    m_noCells = currentLocalCellId;

  } // end of : extra '{}' span for sending / receiving the new cells.

  RECORD_TIMER_STOP(t_recreateCells);

  // =============================
  // #=- Regenerate halo cells -=#
  // =============================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 5 - Regenerate halo cells" << endl;
#endif

  RECORD_TIMER_START(t_createHaloCells);

  // Resetting neighborDomainInformation.
  mDeallocate(m_neighborDomains);
  mDeallocate(m_rfnCountHalosDom);
  m_noNeighborDomains = 0;
  noCellsToReceive.fill(0);

  // Map to collect all missing halo cells
  map<MLong, MInt> haloCellsMap;

  // 1. Get the number of missing parents up to the m_minLevel
  MInt parentDepth = 0;
  MLong parentId = -1;
  for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
    if(globalOffset[globalDomainId()] == a_globalId(m_levelOffsets[level][0])) {
      parentDepth = a_level(m_levelOffsets[level][0]) - m_minLevel;
      parentId = a_parentId(m_levelOffsets[level][0]);
      break;
    }
  }

  m_noMissingParents = parentDepth;
  MInt currentId = 0;

  // 2. Add the first missing parent to the map.
  if(parentId > -1) {
    haloCellsMap.insert(make_pair(parentId, currentId));
    ++currentId;
    --parentDepth;
  }

  // 3. Look for more parents.
  MBool stillLookingForParents = true;
  while(stillLookingForParents) {
    // 3.1. Gather from all domains, the number of parents currently missing.
    MIntScratchSpace parentDepthGathered(noDomains(), AT_, "parentDepthGathered");
    parentDepthGathered.fill(0);
    MPI_Allgather(&parentDepth, 1, MPI_INT, parentDepthGathered.getPointer(), 1, MPI_INT, mpiComm(), AT_, "parentDepth",
                  "parentDepthGathered.getPointer()");

    MLong parentIdToLookFor = -1;
    MPI_Status status;

    // 3.2. Calculate the domain to which our current missing parent belongs.
    MInt parentDomain = -1;
    for(MInt i = 0; i < noDomains(); ++i) {
      if(parentId < globalOffset[i]) {
        parentDomain = i - 1;
        break;
      }
    }

    // 3.3. Gather from all domains, which domain will send us a parent.
    MIntScratchSpace parentDomainGathered(noDomains(), AT_, "parentDomainGathered");
    parentDomainGathered.fill(-1);
    MPI_Allgather(&parentDomain, 1, MPI_INT, parentDomainGathered.getPointer(), 1, MPI_INT, mpiComm(), AT_,
                  "parentDomain", "parentDomainGathered.getPointer()");


    // 3.4. Communication...
    for(MInt dom = 0; dom < noDomains(); ++dom) {
      // 3.4.a) If its our turn:
      if((dom == globalDomainId()) && (parentDepthGathered[globalDomainId()] > 0)) {
        // Send the parent id we have.
        MPI_Send(&parentId, 1, MPI_LONG, parentDomain, parentDomain, mpiComm(), AT_, "parentId");
        // Wait for the parent of the parent.
        MPI_Recv(&parentId, 1, MPI_LONG, parentDomain, parentDomain, mpiComm(), &status, AT_, "parentId");
        // Add the new halo parent to the map.
        //        #ifdef MAIA_EXTRA_DEBUG
        //          m_log << "DynamicLoadBalancing Debug: Adding halo partition level ancestor with parentId " <<
        //          parentId << " to the haloCellsMap." << endl;
        //        #endif
        haloCellsMap.insert(make_pair(parentId, currentId));
        ++currentId;
      } // end of : if ( (dom == globalDomainId()) && (parentDepthGathered[globalDomainId()] > 0) )

      // 3.5.b) If it's not our turn and the other domain sends us a parent:
      if((parentDomainGathered[dom] == globalDomainId()) && (parentDepthGathered[dom] > 0)) {
        // Receive the parent.
        MPI_Recv(&parentIdToLookFor, 1, MPI_LONG, dom, globalDomainId(), mpiComm(), &status, AT_, "parentIdToLookFor");
        // Search for it and save the parent of the parent.
        for(MInt i = 0; i < m_noCells; ++i) {
          if(a_globalId(i) == parentIdToLookFor) {
            parentIdToLookFor = a_parentId(i);
            break;
          }
        }
        // Sends back the parent of the parent.
        MPI_Send(&parentIdToLookFor, 1, MPI_LONG, dom, globalDomainId(), mpiComm(), AT_, "parentIdToLookFor");
      } // end of : if ( (parentDomainGathered[dom] == globalDomainId()) && (parentDepthGathered[dom] > 0) )

    } // end of : for (MInt dom = 0; dom < noDomains(); ++dom)

    // 3.5. Decrease the number of parents we need to find.
    if(parentDepth > 0) {
      --parentDepth;
    }

    // 3.6. Look if any domain still search for more parents.
    MInt continueSearch = 0;
    MPI_Allreduce(&parentDepth, &continueSearch, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "parentDepth", "continueSearch");

    if(continueSearch == 0) {
      stillLookingForParents = false;
    }

  } // end of : while (stillLookingForParents)

  // 4. Taging window cells and adding halo cells to the map.
  for(MInt i = 0; i < m_noCells; ++i) {
    a_hasProperty(i, 3) = 0;
    a_hasProperty(i, 4) = 0;

    // 4.1. Looks for halo neighbors.
    for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
      if(a_neighborId(i, ngh) == -1) {
        continue;
      }
      if((a_neighborId(i, ngh) < globalOffset[globalDomainId()])
         || (a_neighborId(i, ngh) >= globalOffset[globalDomainId()] + m_noCells)) {
        if(haloCellsMap.count(a_neighborId(i, ngh)) == 0) {
          //          #ifdef MAIA_EXTRA_DEBUG
          //            m_log << "DynamicLoadBalancing Debug: Adding halo neighborId " << a_neighborId(i, ngh) << "
          //            to the haloCellsMap." << endl;
          //          #endif
          haloCellsMap.insert(make_pair(a_neighborId(i, ngh), currentId));
          ++currentId;
        }
        a_hasProperty(i, 3) = 1;
      }
    } // end of : for (MInt ngh = 0; ngh < m_noNeighbors; ++ngh)

    // 4.2. Looks for halo childs.
    for(MInt child = 0; child < m_maxNoChildren; ++child) {
      if(a_childId(i, child) == -1) {
        continue;
      }

      if((a_childId(i, child) < globalOffset[globalDomainId()])
         || (a_childId(i, child) >= globalOffset[globalDomainId()] + m_noCells)) {
        if(haloCellsMap.count(a_childId(i, child)) == 0) {
          //          #ifdef MAIA_EXTRA_DEBUG
          //            m_log << "DynamicLoadBalancing Debug: Adding halo childId " << a_childId(i, child) << " to
          //            the haloCellsMap." << endl;
          //          #endif
          haloCellsMap.insert(make_pair(a_childId(i, child), currentId));
          ++currentId;
        }
        a_hasProperty(i, 3) = 1;
      }
    } // end of : for (MInt child = 0; child < m_maxNoChildren; ++child)

    // 4.3. Looks for halo parents.
    if((a_parentId(i) < globalOffset[globalDomainId()])
       || (a_parentId(i) >= globalOffset[globalDomainId()] + m_noCells)) {
      if(a_parentId(i) == -1) {
        continue;
      }

      if(haloCellsMap.count(a_parentId(i)) == 0) {
        //        #ifdef MAIA_EXTRA_DEBUG
        //          m_log << "DynamicLoadBalancing Debug: Adding halo parentId " << a_parentId(i) << " to the
        //          haloCellsMap." << endl;
        //        #endif
        haloCellsMap.insert(make_pair(a_parentId(i), currentId));
        ++currentId;
      }
      a_hasProperty(i, 3) = 1;
    } // end of : if ( (a_parentId(i) < globalOffset[globalDomainId()]) || (a_parentId(i) >=
      // globalOffset[globalDomainId()] + m_noCells) )

  } // end of : for (MInt i = 0; i < m_noCells; ++i)

  // 5. Creating haloInformation which holds:
  //  - [0] == An identifier
  //  - [1] == The global id of the halo cell.
  //  - [2] == The domain to which the corresponding window cell belongs.
  MLongScratchSpace haloInformation((MInt)haloCellsMap.size(), 3, AT_, "haloInformation");
  haloInformation.fill(-2);

  // 5.1. Looks for the domain id of the halo cells
  currentId = 0;
  //  #ifdef MAIA_EXTRA_DEBUG
  //    m_log << "DynamicLoadBalancing Extra Debug: Creating haloInformation array... " << endl;
  //  #endif
  for(map<MLong, MInt>::iterator i = haloCellsMap.begin(); i != haloCellsMap.end(); ++i) {
    for(MLong dom = 0; dom < (MLong)noDomains(); ++dom) {
      if((globalOffset[dom] <= i->first) && (globalOffset[dom + 1] > i->first)) {
        haloInformation(currentId, 0) = (MInt)i->second;
        haloInformation(currentId, 1) = i->first;
        haloInformation(currentId, 2) = dom;
        //        #ifdef MAIA_EXTRA_DEBUG
        //          m_log << "haloInformation(currentId, 0) is " << haloInformation(currentId, 0) << " ";
        //          m_log << "haloInformation(currentId, 1) is " << haloInformation(currentId, 1) << " ";
        //          m_log << "haloInformation(currentId, 2) is " << haloInformation(currentId, 2) << endl;
        //        #endif
        ++currentId;
        ++noCellsToReceive[dom];
        break;
      } // end of : if ( (globalOffset[dom] <= i->first) && (globalOffset[dom+1] > i->first) )
    }   // end of : for (MInt dom = 0; dom < noDomains(); ++dom)

  } // end of : for (map<MInt, MInt>::iterator i = haloCellsMap.begin(); i != haloCellsMap.end(); ++i)

  // 6. Communicating the whole informations of the halo cells.
  // 6.1. Preparation...
  for(MInt i = 0; i < noDomains(); ++i) {
    if(noCellsToReceive[i] > 0) {
      noCellsToReceive[noDomains()] += noCellsToReceive[i];
    }
  }

  noCellsToSend.fill(0);
  for(MInt i = 0; i < noDomains(); ++i) {
    MPI_Scatter(noCellsToReceive.getPointer(), 1, MPI_INT, &noCellsToSend[i], 1, MPI_INT, i, mpiComm(), AT_,
                "noCellsToReceive.getPointer()", "noCellsToSend[i]");
  }

  for(MInt i = 0; i < noDomains(); ++i) {
    if(noCellsToSend[i] > 0) {
      noCellsToSend[noDomains()] += noCellsToSend[i];
    }
  }

  for(MInt i = 0; i < noDomains(); ++i) {
    if((noCellsToSend[i] > 0) || (noCellsToReceive[i] > 0)) {
      ++m_noNeighborDomains;
    }
  }

  mAlloc(m_neighborDomains, m_noNeighborDomains, "m_neighborDomains", -1, AT_);
  mAlloc(m_rfnCountHalosDom, m_noNeighborDomains, "m_rfnCountHalosDom", 0, AT_);

  currentId = 0;
  for(MInt i = 0; i < noDomains(); ++i) {
    if((noCellsToSend[i] > 0) || (noCellsToReceive[i] > 0)) {
      m_neighborDomains[currentId] = i;
      ++currentId;
    }
  }

  // 6.2. Sending first which halo cells are needed to the other domain.
  MLongScratchSpace cellsToSendGlobalId(noCellsToSend[noDomains()], AT_, "cellsToSendGlobalId");

  { // extra '{}' span for communicating globalIds of the window cell we need to receive
    MLongScratchSpace dataToSend(noCellsToReceive[noDomains()], AT_, "dataToSend");
    MLongScratchSpace dataToReceive(noCellsToSend[noDomains()], AT_, "dataToReceive");

    MInt offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      currentId = 0;
      for(MInt i = 0; i < (signed)haloCellsMap.size(); ++i) {
        if(haloInformation(i, 2) == m_neighborDomains[dom]) {
          dataToSend(currentId + offset) = haloInformation(i, 1);
          ++currentId;
        }
      }
      offset += currentId;
    } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)

    communicateLong(dataToReceive, noCellsToSend, dataToSend, noCellsToReceive);

    offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noCellsToSend[m_neighborDomains[dom]]; ++i) {
        //        m_log << "Receiving from Domain " << m_neighborDomains[dom] << " globalId " << dataToReceive(i +
        //        offset) << endl;
        cellsToSendGlobalId(i + offset) = dataToReceive(i + offset);
        //        m_log << "Saved value is " << cellsToSendGlobalId(i + offset) << endl;
      }
      offset += noCellsToSend[m_neighborDomains[dom]];
    } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)
  }   // end of : extra '{}' span for communicating globalIds of the window cell we need to receive

  // 6.3. Send / Receive the new halo cells.
  { // extra '{}' span for sending / receiving the new halo cells.

    MLongScratchSpace newParentId((MInt)haloCellsMap.size(), AT_, "newParentId");
    MLongScratchSpace newGlobalId((MInt)haloCellsMap.size(), AT_, "newGlobalId");
    MIntScratchSpace newLevel((MInt)haloCellsMap.size(), AT_, "newLevel");
    MIntScratchSpace newNoChildren((MInt)haloCellsMap.size(), AT_, "newNoChildren");

    MLongScratchSpace newChildId((MInt)haloCellsMap.size(), m_maxNoChildren, AT_, "newChildId");
    MLongScratchSpace newNeighborId((MInt)haloCellsMap.size(), m_noNeighbors, AT_, "newNeighborId");
    MIntScratchSpace newProperty((MInt)haloCellsMap.size(), 8, AT_, "newProperty");
    MIntScratchSpace newIsInSolver((MInt)haloCellsMap.size(), m_noSolvers, AT_, "newIsInSolver");
    MIntScratchSpace newIsSolverBoundary((MInt)haloCellsMap.size(), m_noSolvers, AT_, "newIsSolverBoundary");
    MIntScratchSpace newIsToRefineForSolver((MInt)haloCellsMap.size(), m_noSolvers, AT_, "newIsToRefineForSolver");
    MFloatScratchSpace newCoordinate((MInt)haloCellsMap.size(), nDim, AT_, "newCoordinate");

    MIntScratchSpace windowCellLocalId(noCellsToSend[noDomains()], AT_, "windowCellLocalId");
    MIntScratchSpace windowCellLocalIdOffset(m_noNeighborDomains, AT_, "windowCellLocalIdOffset");

    // 6.3.1. Getting window cell local ids.
    currentId = 0;
    MInt offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      MInt domain = m_neighborDomains[dom];

      windowCellLocalIdOffset[dom] = currentId;
      //      m_log << "domain is " << domain << endl;
      //      m_log << "noCellsToSend[" << domain<< "] is " << noCellsToSend[domain] << endl;
      //      m_log << "offset is " << offset << endl;

      for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
        for(MInt i = 0; i < noCellsToSend[domain]; ++i) {
          if(cellsToSendGlobalId(i + offset) == a_globalId(cellId)) {
            //            m_log << "Both are Equal!!! ";
            //            m_log << "noCells left to find are " << noCellsToSend[domain] - currentId + offset << endl;
            //            m_log << cellsToSendGlobalId(i + offset) << " == " << a_globalId(cellId) << endl;
            windowCellLocalId[currentId] = cellId;
            ++currentId;
            break;
          } // end of : if (cellsToSendGlobalId(i + offset) == m_cells->a[cellId].m_globalId[0])
        }   // end of : for (MInt i = 0; i < noCellsToSend[domain]; ++i)

      } // end of : for (MInt cellId = 0; cellId < m_noCells; ++cellId)
      offset += noCellsToSend[domain];

    } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)

    MIntScratchSpace haloOffset(m_noNeighborDomains, AT_, "haloOffset");
    if(m_noNeighborDomains > 0) {
      haloOffset[0] = 0;
    }

    for(MInt dom = 1; dom < m_noNeighborDomains; ++dom) {
      haloOffset[dom] = haloOffset[dom - 1] + noCellsToReceive[m_neighborDomains[dom - 1]];
    }

    // 6.3.2. Communicate halo/window data...
    { // sending / receiving parentId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_parentId(windowCellLocalId(i));
          ++currentId;
        }
      }

      communicateLong(newParentId, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving parentId

    { // sending / receiving globalId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_globalId(windowCellLocalId(i));
          ++currentId;
        }
      }

      communicateLong(newGlobalId, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving GlobalId

    { // sending / receiving level
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_level(windowCellLocalId(i));
          ++currentId;
        }
      }

      communicateInt(newLevel, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving level

    { // sending / receiving noChildren
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_noChildren(windowCellLocalId(i));
          ++currentId;
        }
      }

      communicateInt(newNoChildren, noCellsToReceive, sendBuffer, noCellsToSend);
    } // end of : sending / receiving noChildren

    for(MInt childId = 0; childId < m_maxNoChildren; ++childId) { // sending / receiving childId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MLongScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_childId(windowCellLocalId(i), childId);
          ++currentId;
        }
      }

      communicateLong(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newChildId(i, childId) = recvBuffer(i);
      }

    } // end of : sending / receiving childId

    for(MInt neighborId = 0; neighborId < m_noNeighbors; ++neighborId) { // sending / receiving neighborId
      MLongScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MLongScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_neighborId(windowCellLocalId(i), neighborId);
          ++currentId;
        }
      }

      communicateLong(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newNeighborId(i, neighborId) = recvBuffer(i);
      }

    } // end of : sending / receiving neighborId

    for(MInt property = 0; property < 8; ++property) { // sending / receiving property
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_hasProperty(windowCellLocalId(i), property);
          ++currentId;
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newProperty(i, property) = recvBuffer(i);
      }

    } // end of : sending / receiving property

    for(MInt solver = 0; solver < m_noSolvers; solver++) { // sending / receiving of solver affiliation
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_isInSolver(windowCellLocalId(i), solver);
          ++currentId;
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsInSolver(i, solver) = recvBuffer(i);
      }
    }

    for(MInt solver = 0; solver < m_noSolvers; solver++) { // sending / receiving of solver boundary
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_isSolverBoundary(windowCellLocalId(i), solver);
          ++currentId;
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsSolverBoundary(i, solver) = recvBuffer(i);
      }
    }

    for(MInt solver = 0; solver < m_noSolvers; solver++) { // sending / receiving of solver boundary
      MIntScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_isSolverBoundary(windowCellLocalId(i), solver);
          ++currentId;
        }
      }

      communicateInt(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newIsToRefineForSolver(i, solver) = recvBuffer(i);
      }
    }

    for(MInt dim = 0; dim < nDim; ++dim) { // sending / receiving coordinate
      MFloatScratchSpace sendBuffer(noCellsToSend[noDomains()], AT_, "sendBuffer");
      MFloatScratchSpace recvBuffer(noCellsToReceive[noDomains()], AT_, "recvBuffer");

      currentId = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        MInt domain = m_neighborDomains[dom];

        for(MInt i = windowCellLocalIdOffset[dom]; i < (windowCellLocalIdOffset[dom] + noCellsToSend[domain]); ++i) {
          sendBuffer(currentId) = a_coordinate(windowCellLocalId(i), dim);
          ++currentId;
        }
      }

      communicateDouble(recvBuffer, noCellsToReceive, sendBuffer, noCellsToSend);

      for(MInt i = 0; i < noCellsToReceive[noDomains()]; ++i) {
        newCoordinate(i, dim) = recvBuffer(i);
      }

    } // end of : sending / receiving coordinate

    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = haloOffset[dom]; i < (haloOffset[dom] + noCellsToReceive[m_neighborDomains[dom]]); ++i) {
        newProperty(i, 3) = 0;
        newProperty(i, 4) = 1;
      }
    }

    // 7. Recreate halo cells from the communication data.
    mDeallocate(m_haloCellOffsets);
    mAlloc(m_haloCellOffsets, m_noNeighborDomains, 2 * (m_maxRfnmntLvl - m_minLevel + 1), "m_haloCellOffsets", 0, AT_);
    m_noTotalHaloCells = 0;

    MInt currentLocalCellId = m_maxNoCells - 2;
    for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
      m_haloCellOffsetsLevel[level][1] = currentLocalCellId + 1;

      for(MInt dom = m_noNeighborDomains - 1; dom > -1; --dom) {
        m_haloCellOffsets[dom][(2 * (level - m_minLevel)) + 1] = currentLocalCellId + 1;

        for(MInt cellId = (haloOffset[dom] + noCellsToReceive[m_neighborDomains[dom]]) - 1; cellId >= haloOffset[dom];
            --cellId) {
          if(newLevel(cellId) == level) {
            a_parentId(currentLocalCellId) = newParentId(cellId);
            a_globalId(currentLocalCellId) = newGlobalId(cellId);
            a_noChildren(currentLocalCellId) = newNoChildren(cellId);
            a_level(currentLocalCellId) = level;

            for(MInt child = 0; child < m_maxNoChildren; ++child) {
              a_childId(currentLocalCellId, child) = newChildId(cellId, child);
            }

            for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
              a_neighborId(currentLocalCellId, ngh) = newNeighborId(cellId, ngh);
            }

            for(MInt dim = 0; dim < nDim; ++dim) {
              a_coordinate(currentLocalCellId, dim) = newCoordinate(cellId, dim);
            }

            for(MInt property = 0; property < 8; ++property) {
              a_hasProperty(currentLocalCellId, property) = (MBool)newProperty(cellId, property);
            }

            for(MInt solver = 0; solver < m_noSolvers; solver++) {
              a_isInSolver(currentLocalCellId, solver) = (MBool)newIsInSolver(cellId, solver);
              a_isSolverBoundary(currentLocalCellId, solver) = (MBool)newIsSolverBoundary(cellId, solver);
              a_isToRefineForSolver(currentLocalCellId, solver) = (MBool)newIsToRefineForSolver(cellId, solver);
              // not working in parallel yet anyway
              a_noSolidLayer(currentLocalCellId, solver) = -1;
            }

            for(MInt solver = m_noSolvers; solver < 8; solver++) {
              a_isInSolver(currentLocalCellId, solver) = 0;
              a_isSolverBoundary(currentLocalCellId, solver) = 0;
              a_isToRefineForSolver(currentLocalCellId, solver) = 0;
              a_noSolidLayer(currentLocalCellId, solver) = -1;
            }
            --currentLocalCellId;
          } // end of : if (cells->a[cellId].m_level[0] == level)

        } // end of : for (MInt i = (haloOffset[dom] + noCellsToReceive[m_neighborDomains[dom]]) - 1; i >=
          // haloOffset[dom]; --i)

        m_haloCellOffsets[dom][(2 * (level - m_minLevel))] = currentLocalCellId + 1;

      } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)

      m_haloCellOffsetsLevel[level][0] = currentLocalCellId + 1;
      m_noHaloCellsOnLevel[level] = m_haloCellOffsetsLevel[level][1] - m_haloCellOffsetsLevel[level][0];
      m_noTotalHaloCells += m_noHaloCellsOnLevel[level];
    } // end of : for (MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level)

  } // end of : extra '{}' span for sending / receiving the new halo cells.

  RECORD_TIMER_STOP(t_createHaloCells);

  // =======================
  // #=- Global to local -=#
  // =======================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 6 - Global to local" << endl;
#endif

  RECORD_TIMER_START(t_globalToLocal);

  map<MLong, MInt> globalToLocal;

  // 1. Add local cells to the map
  for(MInt i = 0; i < m_noCells; ++i) {
    globalToLocal.insert(make_pair(a_globalId(i), i));
  }

  // 2. Add halo cells to the map
  for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
    if(m_noHaloCellsOnLevel[level] > 0) {
      for(MInt i = m_haloCellOffsetsLevel[level][0]; i < m_haloCellOffsetsLevel[level][1]; ++i) {
        globalToLocal.insert(make_pair(a_globalId(i), i));
      }
    }
  }

  // 3. Convert local cells
  for(MInt i = 0; i < m_noCells; ++i) {
    if(globalToLocal.count(a_parentId(i)) == 0) {
#ifdef MAIA_EXTRA_DEBUG
      if(a_parentId(i) > -1) {
        m_log << "DynamicLoadBalancing Debug: global parenId " << a_parentId(i) << " does not exist" << endl;
      }
#endif
      a_parentId(i) = -1;
    } else {
      a_parentId(i) = globalToLocal[a_parentId(i)];
    }

    for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
      if(globalToLocal.count(a_neighborId(i, ngh)) == 0) {
#ifdef MAIA_EXTRA_DEBUG
        if(a_neighborId(i, ngh) > -1) {
          m_log << "DynamicLoadBalancing Debug: global nghbrId " << a_neighborId(i, ngh) << " does not exist" << endl;
        }
#endif
        a_neighborId(i, ngh) = -1;
      } else {
        a_neighborId(i, ngh) = globalToLocal[a_neighborId(i, ngh)];
      }
    }

    a_noChildren(i) = m_maxNoChildren;
    for(MInt child = 0; child < m_maxNoChildren; ++child) {
      if(globalToLocal.count(a_childId(i, child)) == 0) {
#ifdef MAIA_EXTRA_DEBUG
        if(a_childId(i, child) > -1) {
          m_log << "DynamicLoadBalancing Debug: global childId " << a_childId(i, child) << " does not exist" << endl;
        }
#endif
        a_childId(i, child) = -1;
        --a_noChildren(i);
      } else {
        a_childId(i, child) = globalToLocal[a_childId(i, child)];
      }
    }

  } // end of : for (MInt i = 0; i < m_noCells; ++i)

  // 4. Convert halo cells
  for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
    if(m_noHaloCellsOnLevel[level] > 0) {
      for(MInt i = m_haloCellOffsetsLevel[level][0]; i < m_haloCellOffsetsLevel[level][1]; ++i) {
        if(globalToLocal.count(a_parentId(i)) == 0) {
          a_parentId(i) = -1;
        } else {
          a_parentId(i) = globalToLocal[a_parentId(i)];
        }

        for(MInt ngh = 0; ngh < m_noNeighbors; ++ngh) {
          if(globalToLocal.count(a_neighborId(i, ngh)) == 0) {
            a_neighborId(i, ngh) = -1;
          } else {
            a_neighborId(i, ngh) = globalToLocal[a_neighborId(i, ngh)];
          }
        }

        a_noChildren(i) = 0;
        for(MInt child = 0; child < m_maxNoChildren; ++child) {
          if(globalToLocal.count(a_childId(i, child)) == 0) {
            a_childId(i, child) = -1;
          } else {
            a_childId(i, child) = globalToLocal[a_childId(i, child)];
            ++a_noChildren(i);
          }
        }

      } // end of : for (MInt i = m_haloCellOffsetsLevel[level][0]; i < m_haloCellOffsetsLevel[level][1]; ++i)

    } // end of : if (m_noHaloCellsOnLevel[level] > 0)

  } // end of : for (MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level)

  RECORD_TIMER_STOP(t_globalToLocal);

  // ==========================
  // #=- Reset lookup table -=#
  // ==========================

#ifdef MAIA_EXTRA_DEBUG
  m_log << "DynamicLoadBalancing Extra Debug: Step 6 - Reset lookup table" << endl;
#endif

  m_cellIdLUT.clear();

  for(MInt level = m_minLevel; level <= m_maxRfnmntLvl; ++level) {
    for(MInt cellId = m_levelOffsets[level][0]; cellId < m_levelOffsets[level][1]; cellId++) {
      m_cellIdLUT.insert(m_cellIdLUT.end(), make_pair(cellId, cellId));
    }

    for(MInt cellId = m_haloCellOffsetsLevel[level][0]; cellId < m_haloCellOffsetsLevel[level][1]; cellId++) {
      m_cellIdLUT.insert(m_cellIdLUT.end(), make_pair(cellId, cellId));
    }
  }

  RECORD_TIMER_STOP(t_dynamicLoadBalancing);
  return;
} // end of : void GridgenPar::dynamicLoadBalancing()

/** \brief communicate the global ids of the halo cells.
 *
 * \author Jerry Grimmen
 * \date 26.02.2014
 *
 * The function communicates the new global ids of the halo cells.
 *  1. If we are already finished with building the grid, then call reorderGlobalIdsDF
 *  2. Collect all local halo cell ids into a vector
 *  3. Get the new global ids of the halo cell.
 *   In this part, we send the coordinates of our halo cells to the other domain.
 *   The neighbor looks through it's window cell for the right one and sends back the new global id.
 *  4. Changes the local ids to global ids of internal cells.
 *
 *  The coordinate check is currently the only way i found out to get the right new global id.
 *  Sending the level of the halo cell would reduce the amoung of cells to look through,
 *  but it will also add another communication.
 *
 */

template <MInt nDim>
void GridgenPar<nDim>::communicateHaloGlobalIds(MInt level) {
  TRACE();

  // 1. If we are at the end of the grid creation use reoderGlobalIdsDF
  if(level > -1) {
    // 1.1. Get the number of cells created on the other domains.
    MInt* sndBuf = &m_noCells;
    MPI_Allgather(sndBuf, 1, MPI_INT, m_noCellsPerDomain, 1, MPI_INT, mpiComm(), AT_, "sndBuf", "m_noCellsPerDomain");

    // 1.2. Determine my own offset.
    for(MInt d = 0; d < globalDomainId(); d++)
      m_cellOffsetPar += (MLong)m_noCellsPerDomain[d];

    // 1.3. Reorder the global ids to be sorted depth-first.
    reorderGlobalIdsDF();
  }

  // 2. Collect all halo cells local ids.
  vector<vector<MLong>> winCellIdsPerDomain;
  vector<vector<MLong>> haloCellIdsPerDomain;

  // 2.1. Create vectors to hold the halo cells local id.
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    vector<MLong> w;
    winCellIdsPerDomain.push_back(w);

    vector<MLong> h;
    haloCellIdsPerDomain.push_back(h);

    for(MInt i = 0; i < m_noCells; ++i) {
      a_hasProperty(i, 3) = 0;
      a_hasProperty(i, 4) = 0;
    }

    // 2.2. Collect the halo cells local ids
    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      MInt lev_pos = 2 * (lev - m_minLevel);
      for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++) {
        haloCellIdsPerDomain[dom].push_back(m_cellIdLUT[p]);
      }
    }
  }

  // 3. Get the halo cells global ids.
  { // extra '{ }' span for communication of the globalId of the halo cells.

    // 3.1. Scatter to other domains the number of halo cells we have.
    MIntScratchSpace noDataToSend(noDomains() + 1, AT_, "noDataToSend");
    MIntScratchSpace noDataToReceive(noDomains() + 1, AT_, "noDataToReceive");
    noDataToSend.fill(0);
    noDataToReceive.fill(0);

    for(MInt i = 0; i < m_noNeighborDomains; ++i) {
      noDataToReceive[m_neighborDomains[i]] = (MInt)haloCellIdsPerDomain[i].size();
      noDataToReceive[noDomains()] += (MInt)haloCellIdsPerDomain[i].size();
    }

    for(MInt d = 0; d < noDomains(); ++d) {
      MPI_Scatter(noDataToReceive.getPointer(), 1, MPI_INT, &noDataToSend[d], 1, MPI_INT, d, mpiComm(), AT_,
                  "noDataToReceive.getPointer()", "noDataToSend[d]");
    }

    for(MInt i = 0; i < noDomains(); ++i) {
      noDataToSend[noDomains()] += noDataToSend[i];
    }

    // 3.2. Prepare an array to hold all halo cells coordinates and level we will receive.
    MFloatScratchSpace haloCoordinates(noDataToSend[noDomains()], nDim, AT_, "haloCoordinates");
    MIntScratchSpace haloLevel(noDataToSend[noDomains()], AT_, "haloLevel");

    { // extra '{}' for sending / receiving halo cells coordinates.
      for(MInt dim = 0; dim < nDim; ++dim) {
        MFloatScratchSpace sendBuffer(noDataToSend[noDomains()], AT_, "sendBuffer");
        MFloatScratchSpace recvBuffer(noDataToReceive[noDomains()], AT_, "recvBuffer");

        MInt offset = 0;
        for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
          for(MInt i = 0; i < noDataToReceive[m_neighborDomains[dom]]; ++i) {
            recvBuffer(i + offset) = a_coordinate((MInt)haloCellIdsPerDomain[dom][i], dim);
          }
          offset += noDataToReceive[m_neighborDomains[dom]];
        }

        communicateDouble(sendBuffer, noDataToSend, recvBuffer, noDataToReceive);

        offset = 0;
        for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
          for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
            haloCoordinates(i + offset, dim) = sendBuffer(i + offset);
          }
          offset += noDataToSend[m_neighborDomains[dom]];
        }
      }
    } // end of : extra '{}' for sending / receiving halo cells coordinates.

    { // extra '{}' for sending / receiving halo cells level.
      MIntScratchSpace sendBuffer(noDataToSend[noDomains()], AT_, "sendBuffer");
      MIntScratchSpace recvBuffer(noDataToReceive[noDomains()], AT_, "recvBuffer");

      MInt offset = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        for(MInt i = 0; i < noDataToReceive[m_neighborDomains[dom]]; ++i) {
          recvBuffer(i + offset) = a_level((MInt)haloCellIdsPerDomain[dom][i]);
        }
        offset += noDataToReceive[m_neighborDomains[dom]];
      }

      communicateInt(sendBuffer, noDataToSend, recvBuffer, noDataToReceive);

      offset = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
          haloLevel(i + offset) = sendBuffer(i + offset);
        }
        offset += noDataToSend[m_neighborDomains[dom]];
      }
    } // end of : extra '{}' for sending / receiving halo cells level.

    // 3.3. Get the correct window cell using coordinate check and add it to the vector.

    vector<Point<3>> pts;
    for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
      IF_CONSTEXPR(nDim == 2) {
        Point<3> pt(a_coordinate(cellId, 0), a_coordinate(cellId, 1), 0.0, cellId);
        pts.push_back(pt);
      }
      else {
        Point<3> pt(a_coordinate(cellId, 0), a_coordinate(cellId, 1), a_coordinate(cellId, 2), cellId);
        pts.push_back(pt);
      }
    }

    KDtree<3> tree(pts);

    MInt offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
        MInt cellId = -2;
        MFloat distance = -2.0;
        IF_CONSTEXPR(nDim == 2) {
          Point<3> pt(haloCoordinates(i + offset, 0), haloCoordinates(i + offset, 1), 0.0, cellId);
          cellId = tree.nearest(pt, distance);
        }
        else {
          Point<3> pt(haloCoordinates(i + offset, 0), haloCoordinates(i + offset, 1), haloCoordinates(i + offset, 2),
                      cellId);
          cellId = tree.nearest(pt, distance);
        }

        winCellIdsPerDomain[dom].push_back(a_globalId(cellId));

      } // end of : for (MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i)
      offset += noDataToSend[m_neighborDomains[dom]];
    } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)

    // 3.4. Send and receive the halo cells global ids.
    { // extra '{}' span for sending / receiving halo cells global ids.
      MLongScratchSpace sendBuffer(noDataToSend[noDomains()], AT_, "sendBuffer");
      MLongScratchSpace recvBuffer(noDataToReceive[noDomains()], AT_, "recvBuffer");
      sendBuffer.fill(-2);
      recvBuffer.fill(-2);

      offset = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        for(MInt i = 0; i < (signed)winCellIdsPerDomain[dom].size(); ++i) {
          sendBuffer(i + offset) = winCellIdsPerDomain[dom][i];
        }
        offset += (MInt)winCellIdsPerDomain[dom].size();
      }

      communicateLong(recvBuffer, noDataToReceive, sendBuffer, noDataToSend);

      offset = 0;
      for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
        for(MInt i = 0; i < (signed)haloCellIdsPerDomain[dom].size(); ++i) {
          a_globalId((MInt)haloCellIdsPerDomain[dom][i]) = recvBuffer(i + offset);
        }
        offset += (signed)haloCellIdsPerDomain[dom].size();
      }

    } // end of : extra '{}' span for sending / receiving halo cells global ids.

  } // end of : extra '{ }' span for communication of the globalId of the halo cells.

  // 4. For all internal cells, change all local ids to global ids.
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    if(a_parentId(cellId) > -1) {
      a_parentId(cellId) = a_globalId((MInt)a_parentId(cellId));
    }

    for(MInt i = 0; i < m_maxNoChildren; ++i) {
      if(a_childId(cellId, i) > -1) {
        a_childId(cellId, i) = a_globalId((MInt)a_childId(cellId, i));
      }
    }

    for(MInt i = 0; i < m_noNeighbors; ++i) {
      if(a_neighborId(cellId, i) > -1) {
        a_neighborId(cellId, i) = a_globalId((MInt)a_neighborId(cellId, i));
      }
    }
  } // end of : for (MInt cellId = 0; cellId < m_noCells; ++cellId)
}

/** \brief find window cells using the known halo cells
 *
 * \author Thomas Schilden
 * \date 3.2.2016
 *
 * \param[in] winCellIdsPerDomain : vector of window cells to be filled
 * \param[in] haloCellIdsPerDomain: vector of halo cells to be filled
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::findHaloAndWindowCells(vector<vector<MInt>>& winCellIdsPerDomain,
                                              vector<vector<MInt>>& haloCellIdsPerDomain) {
  TRACE();

  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    for(MInt i = 0; i < m_noCells; i++) {
      a_hasProperty(i, 3) = 0;
    }

    // collect the halo cells and mark the window cells
    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      MInt cnt = 0;
      MInt lev_pos = 2 * (lev - m_minLevel);
      for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++) {
        const MInt cellId = m_cellIdLUT[p];
        haloCellIdsPerDomain[dom].push_back(cellId);
        for(MInt n = 0; n < m_noNeighbors; n++) {
          if(a_neighborId(cellId, n) < m_noCells && a_neighborId(cellId, n) > -1) {
            a_hasProperty((MInt)a_neighborId(cellId, n), 3) = 1;
          }
        }
        cnt++;
      }
    }

    // collect the window cells
    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      MInt cnt = 0;
      for(MInt i = m_levelOffsets[lev][0]; i < m_levelOffsets[lev][1]; i++) {
        if(a_hasProperty(m_cellIdLUT[i], 3)) {
          winCellIdsPerDomain[dom].push_back(m_cellIdLUT[i]);
          cnt++;
        }
      }
    }
  }
}


/** \brief finds halo and window cells using a kd tree
 *
 * \author Thomas Schilden
 * \date 25.08.2015
 *
 *  does the same as findHaloAndWindowCells
 *
 *  Jerry: The coordinate check is currently the only way i found out to get the right ids.
 *  Sending the level of the halo cell reduces the amoung of cells to look through,
 *  but it will also add another communication.
 *
 */
template <MInt nDim>
void GridgenPar<nDim>::findHaloAndWindowCellsKD(vector<vector<MInt>>& winCellIdsPerDomain,
                                                vector<vector<MInt>>& haloCellIdsPerDomain) {
  TRACE();
  // 1. fill vectors that hold the halo cells local id.
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    for(MInt i = 0; i < m_noCells; ++i) {
      a_hasProperty(i, 3) = 0;
      a_hasProperty(i, 4) = 0;
    }

    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      MInt lev_pos = 2 * (lev - m_minLevel);
      for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++) {
        haloCellIdsPerDomain[dom].push_back(m_cellIdLUT[p]);
      }
    }
  }

  // 2. Scatter the number of halo cells between the domains
  MIntScratchSpace noDataToSend(noDomains() + 1, AT_, "noDataToSend");
  MIntScratchSpace noDataToReceive(noDomains() + 1, AT_, "noDataToReceive");
  noDataToSend.fill(0);
  noDataToReceive.fill(0);

  for(MInt i = 0; i < m_noNeighborDomains; ++i) {
    noDataToReceive[m_neighborDomains[i]] = (MInt)haloCellIdsPerDomain[i].size();
    noDataToReceive[noDomains()] += (MInt)haloCellIdsPerDomain[i].size();
  }

  MPI_Alltoall(noDataToReceive.getPointer(), 1, MPI_INT, noDataToSend.getPointer(), 1, MPI_INT, mpiComm(), AT_,
               "noDataToReceive.getPointer()", "noDataToSend.getPointer()");

  for(MInt i = 0; i < noDomains(); ++i) {
    noDataToSend[noDomains()] += noDataToSend[i];
  }

  // 3. distribute the halo cells coordinates among the domains that have to communicate
  MFloatScratchSpace haloCoordinates(noDataToSend[noDomains()], nDim, AT_, "haloCoordinates");
  for(MInt dim = 0; dim < nDim; ++dim) {
    MFloatScratchSpace sendBuffer(noDataToSend[noDomains()], AT_, "sendBuffer");
    MFloatScratchSpace recvBuffer(noDataToReceive[noDomains()], AT_, "recvBuffer");

    MInt offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noDataToReceive[m_neighborDomains[dom]]; ++i) {
        recvBuffer(i + offset) = a_coordinate(haloCellIdsPerDomain[dom][i], dim);
      }
      offset += noDataToReceive[m_neighborDomains[dom]];
    }

    communicateDouble(sendBuffer, noDataToSend, recvBuffer, noDataToReceive);

    offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
        haloCoordinates(i + offset, dim) = sendBuffer(i + offset);
      }
      offset += noDataToSend[m_neighborDomains[dom]];
    }
  }

  // 4. distribute the halo cell levels among the domains that have to communicate
  MIntScratchSpace haloLevel(noDataToSend[noDomains()], AT_, "haloLevel");
  { // extra '{}' for sending / receiving halo cells level.
    MIntScratchSpace sendBuffer(noDataToSend[noDomains()], AT_, "sendBuffer");
    MIntScratchSpace recvBuffer(noDataToReceive[noDomains()], AT_, "recvBuffer");

    MInt offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noDataToReceive[m_neighborDomains[dom]]; ++i) {
        recvBuffer(i + offset) = a_level(haloCellIdsPerDomain[dom][i]);
      }
      offset += noDataToReceive[m_neighborDomains[dom]];
    }

    communicateInt(sendBuffer, noDataToSend, recvBuffer, noDataToReceive);

    offset = 0;
    for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
      for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
        haloLevel(i + offset) = sendBuffer(i + offset);
      }
      offset += noDataToSend[m_neighborDomains[dom]];
    }
  } // end of : extra '{}' for sending / receiving halo cells level.

  // 5. Get the correct window cell using coordinate check and add it to the vector.
  vector<Point<3>> pts;
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    IF_CONSTEXPR(nDim == 2) {
      Point<3> pt(a_coordinate(cellId, 0), a_coordinate(cellId, 1), 0.0, cellId);
      pts.push_back(pt);
    }
    else {
      Point<3> pt(a_coordinate(cellId, 0), a_coordinate(cellId, 1), a_coordinate(cellId, 2), cellId);
      pts.push_back(pt);
    }
  }

  KDtree<3> tree(pts);

  MInt offset = 0;
  for(MInt dom = 0; dom < m_noNeighborDomains; ++dom) {
    for(MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i) {
      MInt cellId = -2;
      MFloat distance = -2.0;
      IF_CONSTEXPR(nDim == 2) {
        Point<3> pt(haloCoordinates(i + offset, 0), haloCoordinates(i + offset, 1), 0.0, cellId);
        cellId = tree.nearest(pt, distance);
      }
      else {
        Point<3> pt(haloCoordinates(i + offset, 0), haloCoordinates(i + offset, 1), haloCoordinates(i + offset, 2),
                    cellId);
        cellId = tree.nearest(pt, distance);
      }

      // contains the global Id of the cell that i have to send
      winCellIdsPerDomain[dom].push_back(cellId);

    } // end of : for (MInt i = 0; i < noDataToSend[m_neighborDomains[dom]]; ++i)
    offset += noDataToSend[m_neighborDomains[dom]];
  } // end of : for (MInt dom = 0; dom < m_noNeighborDomains; ++dom)
}


/** \brief Communicates an array of int values.
 *
 * \author Jerry Grimmen
 * \date 10.02.2014
 *
 * The function sends the parts of sendBuffer to the other domains.
 *  1. Calculate the offsets using the noCellsToSend/Receive arrays.
 *  2. Loop over all domains.
 *  2.1. If the domain is lower or greater than mine and I need to receive cells, open a MPI_Receive.
 *  2.2. If the domain is my domain, loop over all domains.
 *  2.2.1. If I need to send something to the other domain, open a MPI_Send.
 *
 * Even through it looks like a complete blocking communication, it is used mostly during the load balancing.
 * At this stage, each domain only communicate with the it's directly neighboring domains.
 * So that most blocking communication takes places in parallel.
 * I.e. while domain 3 is receiving from 4 and 5 and sending to 1 and 2,
 * domain 1223 is receiving from 1224 and 1225 and sending to 1221 and 1222.
 *
 **/

template <MInt nDim>
void GridgenPar<nDim>::communicateInt(MIntScratchSpace& recvBuffer, MIntScratchSpace& noCellsToReceive,
                                      MIntScratchSpace& sendBuffer, MIntScratchSpace& noCellsToSend) {
  MPI_Status statusRecv;

  MIntScratchSpace sendOffset(noDomains(), AT_, "sendOffset");
  MIntScratchSpace recvOffset(noDomains(), AT_, "recvOffset");
  sendOffset.fill(0);
  recvOffset.fill(0);

  for(MInt i = 1; i < noDomains(); ++i) {
    sendOffset[i] = sendOffset[i - 1] + noCellsToSend[i - 1];
    recvOffset[i] = recvOffset[i - 1] + noCellsToReceive[i - 1];
  }

  for(MInt dom = 0; dom < noDomains(); ++dom) {
    if((dom < globalDomainId()) || (dom > globalDomainId())) {
      if(noCellsToReceive[dom] == 0) {
        continue;
      }

      MIntScratchSpace recvData(noCellsToReceive[dom], AT_, "recvData");
      recvData.fill(-2);
      MPI_Recv(recvData.getPointer(), noCellsToReceive[dom], MPI_INT, dom, dom, mpiComm(), &statusRecv, AT_,
               "recvData.getPointer()");
      for(MInt i = 0; i < noCellsToReceive[dom]; ++i) {
        recvBuffer[i + recvOffset[dom]] = recvData[i];
      }
    } // end of : if ( (dom < globalDomainId()) || (dom > globalDomainId()) )

    if(dom == globalDomainId()) {
      for(MInt toDom = 0; toDom < noDomains(); ++toDom) {
        if((noCellsToSend[toDom] == 0) || (globalDomainId() == toDom)) {
          continue;
        }

        MIntScratchSpace sendData(noCellsToSend[toDom], AT_, "sendData");
        sendData.fill(-2);
        for(MInt i = 0; i < noCellsToSend[toDom]; ++i) {
          sendData[i] = sendBuffer[i + sendOffset[toDom]];
        }

        MPI_Send(sendData.getPointer(), noCellsToSend[toDom], MPI_INT, toDom, globalDomainId(), mpiComm(), AT_,
                 "sendData.getPointer()");
      } // end of : for (MInt toDom = 0; toDom < noDomains(); ++toDom)
    }   // end of : if (dom == globalDomainId())

  } // end of : for (MInt dom = 0; dom < noDomains(); ++dom)

} // end of : void GridgenPar::communicateInt(MIntScratchSpace & recvBuffer, MIntScratchSpace & noCellsToReceive,
  // MIntScratchSpace & sendBuffer, MIntScratchSpace & noCellsToSend)

/** \brief Communicates an array of int values.
 *
 * \author Jerry Grimmen
 * \date 10.02.2014
 *
 * The function sends the parts of sendBuffer to the other domains.
 *  1. Calculate the offsets using the noCellsToSend/Receive arrays.
 *  2. Loop over all domains.
 *  2.1. If the domain is lower or greater than mine and I need to receive cells, open a MPI_Receive.
 *  2.2. If the domain is my domain, loop over all domains.
 *  2.2.1. If I need to send something to the other domain, open a MPI_Send.
 *
 * Even through it looks like a complete blocking communication, it is used mostly during the load balancing.
 * At this stage, each domain only communicate with the it's directly neighboring domains.
 * So that most blocking communication takes places in parallel.
 * I.e. while domain 3 is receiving from 4 and 5 and sending to 1 and 2,
 * domain 1223 is receiving from 1224 and 1225 and sending to 1221 and 1222.
 *
 **/

template <MInt nDim>
void GridgenPar<nDim>::communicateLong(MLongScratchSpace& recvBuffer, MIntScratchSpace& noCellsToReceive,
                                       MLongScratchSpace& sendBuffer, MIntScratchSpace& noCellsToSend) {
  MPI_Status statusRecv;

  MIntScratchSpace sendOffset(noDomains(), AT_, "sendOffset");
  MIntScratchSpace recvOffset(noDomains(), AT_, "recvOffset");
  sendOffset.fill(0);
  recvOffset.fill(0);

  for(MInt i = 1; i < noDomains(); ++i) {
    sendOffset[i] = sendOffset[i - 1] + noCellsToSend[i - 1];
    recvOffset[i] = recvOffset[i - 1] + noCellsToReceive[i - 1];
  }

  for(MInt dom = 0; dom < noDomains(); ++dom) {
    if((dom < globalDomainId()) || (dom > globalDomainId())) {
      if(noCellsToReceive[dom] == 0) {
        continue;
      }

      MLongScratchSpace recvData(noCellsToReceive[dom], AT_, "recvData");
      recvData.fill(-2);
      MPI_Recv(recvData.getPointer(), noCellsToReceive[dom], MPI_LONG, dom, dom, mpiComm(), &statusRecv, AT_,
               "recvData.getPointer()");
      for(MInt i = 0; i < noCellsToReceive[dom]; ++i) {
        recvBuffer[i + recvOffset[dom]] = recvData[i];
      }
    } // end of : if ( (dom < globalDomainId()) || (dom > globalDomainId()) )

    if(dom == globalDomainId()) {
      for(MInt toDom = 0; toDom < noDomains(); ++toDom) {
        if((noCellsToSend[toDom] == 0) || (globalDomainId() == toDom)) {
          continue;
        }

        MLongScratchSpace sendData(noCellsToSend[toDom], AT_, "sendData");
        sendData.fill(-2);
        for(MInt i = 0; i < noCellsToSend[toDom]; ++i) {
          sendData[i] = sendBuffer[i + sendOffset[toDom]];
        }

        MPI_Send(sendData.getPointer(), noCellsToSend[toDom], MPI_LONG, toDom, globalDomainId(), mpiComm(), AT_,
                 "sendData.getPointer()");
      } // end of : for (MInt toDom = 0; toDom < noDomains(); ++toDom)
    }   // end of : if (dom == globalDomainId())

  } // end of : for (MInt dom = 0; dom < noDomains(); ++dom)

} // end of : void GridgenPar::communicateInt(MIntScratchSpace & recvBuffer, MIntScratchSpace & noCellsToReceive,
  // MIntScratchSpace & sendBuffer, MIntScratchSpace & noCellsToSend)


/** \brief Communicates an array of double values.
 *
 * \author Jerry Grimmen
 * \date 10.02.2014
 *
 * The function sends the parts of sendBuffer to the other domains.
 *  1. Calculate the offsets using the noCellsToSend/Receive arrays.
 *  2. Loop over all domains.
 *  2.1. If the domain is lower or greater than mine and I need to receive cells, open a MPI_Receive.
 *  2.2. If the domain is my domain, loop over all domains.
 *  2.2.1. If I need to send something to the other domain, open a MPI_Send.
 *
 * Even through it looks like a complete blocking communication, it is used mostly during the load balancing.
 * At this stage, each domain only communicate with the it's directly neighboring domains.
 * So that most blocking communication takes places in parallel.
 * I.e. while domain 3 is receiving from 4 and 5 and sending to 1 and 2,
 * domain 1223 is receiving from 1224 and 1225 and sending to 1221 and 1222.
 *
 **/

template <MInt nDim>
void GridgenPar<nDim>::communicateDouble(MFloatScratchSpace& recvBuffer, MIntScratchSpace& noCellsToReceive,
                                         MFloatScratchSpace& sendBuffer, MIntScratchSpace& noCellsToSend) {
  MPI_Status statusRecv;

  MIntScratchSpace sendOffset(noDomains(), AT_, "sendOffset");
  MIntScratchSpace recvOffset(noDomains(), AT_, "recvOffset");
  sendOffset.fill(0);
  recvOffset.fill(0);

  for(MInt i = 1; i < noDomains(); ++i) {
    sendOffset[i] = sendOffset[i - 1] + noCellsToSend[i - 1];
    recvOffset[i] = recvOffset[i - 1] + noCellsToReceive[i - 1];
  }

  for(MInt dom = 0; dom < noDomains(); ++dom) {
    if((dom < globalDomainId()) || (dom > globalDomainId())) {
      if(noCellsToReceive[dom] == 0) {
        continue;
      }

      MFloatScratchSpace recvData(noCellsToReceive[dom], AT_, "recvData");
      recvData.fill(-2.0);
      MPI_Recv(recvData.getPointer(), noCellsToReceive[dom], MPI_DOUBLE, dom, dom, mpiComm(), &statusRecv, AT_,
               "recvData.getPointer()");
      for(MInt i = 0; i < noCellsToReceive[dom]; ++i) {
        recvBuffer[i + recvOffset[dom]] = recvData[i];
      }
    } // end of : if ( (dom < globalDomainId()) || (dom > globalDomainId()) )

    if(dom == globalDomainId()) {
      for(MInt toDom = 0; toDom < noDomains(); ++toDom) {
        if((noCellsToSend[toDom] == 0) || (globalDomainId() == toDom)) {
          continue;
        }

        MFloatScratchSpace sendData(noCellsToSend[toDom], AT_, "sendData");
        sendData.fill(-2.0);
        for(MInt i = 0; i < noCellsToSend[toDom]; ++i) {
          sendData[i] = sendBuffer[i + sendOffset[toDom]];
        }

        MPI_Send(sendData.getPointer(), noCellsToSend[toDom], MPI_DOUBLE, toDom, globalDomainId(), mpiComm(), AT_,
                 "sendData.getPointer()");
      } // end of : for (MInt toDom = 0; toDom < noDomains(); ++toDom)
    }   // end of : if (dom == globalDomainId())

  } // end of : for (MInt dom = 0; dom < noDomains(); ++dom)

} // end of : void GridgenPar::communicateDouble(MFloatScratchSpace & recvBuffer, MIntScratchSpace &
  // noCellsToReceive, MFloatScratchSpace & sendBuffer, MIntScratchSpace & noCellsToSend)

/** \brief Communicates Int values, with a variable nmbr of variables to the Neighbors
 *
 * \author Thomas Schilden
 * \date 26.08.2015
 *
 * a check whether a send-recv mismatch exist is not performed, mpi will tell you
 **/
template <MInt nDim>
void GridgenPar<nDim>::communicateIntToNeighbors(MIntScratchSpace& recvMem, MIntScratchSpace& noCellsToReceive,
                                                 MIntScratchSpace& sendMem, MIntScratchSpace& noCellsToSend,
                                                 MInt noVar) {
  MIntScratchSpace sendOffset(m_noNeighborDomains, AT_, "sendOffset");
  MIntScratchSpace recvOffset(m_noNeighborDomains, AT_, "recvOffset");
  sendOffset.fill(0);
  recvOffset.fill(0);

  for(MInt i = 1; i < m_noNeighborDomains; i++) {
    sendOffset[i] = sendOffset[i - 1] + noCellsToSend[i - 1];
    recvOffset[i] = recvOffset[i - 1] + noCellsToReceive[i - 1];
  }

  MInt noRecv = 0;
  for(MInt i = 0; i < m_noNeighborDomains; i++)
    if(noCellsToReceive[i] > 0) noRecv++;

  ScratchSpace<MPI_Request> request(noRecv, AT_, "request");
  ScratchSpace<MPI_Status> status(noRecv, AT_, "status");

  for(MInt i = 0, rCnt = 0; i < m_noNeighborDomains; i++)
    if(noCellsToReceive[i] > 0) {
      MPI_Irecv(recvMem.getPointer() + noVar * recvOffset[i], noVar * noCellsToReceive[i], MPI_INT,
                m_neighborDomains[i], 0, MPI_COMM_WORLD, request.getPointer() + rCnt, AT_,
                "recvMem.getPointer() + noVar * recvOffset[i]");
      rCnt++;
    }

  for(MInt i = 0; i < m_noNeighborDomains; i++)
    if(noCellsToSend[i] > 0) {
      MPI_Send(sendMem.getPointer() + noVar * sendOffset[i], noVar * noCellsToSend[i], MPI_INT, m_neighborDomains[i], 0,
               MPI_COMM_WORLD, AT_, "sendMem.getPointer() + noVar * sendOffset[i]");
    }

  MPI_Waitall(noRecv, request.getPointer(), status.getPointer(), AT_);
}


//##################################################################################################################
//##################################################################################################################
//####                                               IO ROUTINES                                                ####
//##################################################################################################################
//##################################################################################################################

/** \brief writes the grid to file in parallel
 *
 * \author Lennart Schneiders
 * \date July 2017
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::saveGrid() {
  TRACE();

  RECORD_TIMER_START(m_t_saveGrid);

  m_log << "  (7) writing grid to file" << endl;
  outStream << "  (7) writing grid to file" << endl;

  if(m_noTotalCells > numeric_limits<MInt>::max()) {
    outStream << "Exceeding 32bit boundary." << endl;
  }

  m_noPartitionCells = (MInt)m_partitionCellList.size();

  // a. file creation and preprocessing
  using namespace maia::parallel_io;
  MString filename = m_outputDir;
  filename.append(m_gridOutputFileName);

  MInt maxLevel = 0;
  MInt minLevel = m_maxLevels;
  for(MInt j = 0; j < m_noCells; ++j) {
    maxLevel = mMax(a_level(j), maxLevel);
    minLevel = mMin(a_level(j), minLevel);
  }

  MLong maxNoCells = m_noCells;
  MPI_Allreduce(MPI_IN_PLACE, &maxNoCells, 1, type_traits<MLong>::mpiType(), MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxNoCells");
  m_log << "Save grid: maximum number of cells per rank: " << maxNoCells << std::endl;

  MLong maxNoHalos = m_noTotalHaloCells;
  MPI_Allreduce(MPI_IN_PLACE, &maxNoHalos, 1, type_traits<MLong>::mpiType(), MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxNoHalos");
  m_log << "Save grid: maximum number of halo cells per rank: " << maxNoHalos << std::endl;

  MPI_Allreduce(MPI_IN_PLACE, &maxLevel, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "maxLevel");
  MPI_Allreduce(MPI_IN_PLACE, &minLevel, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "minLevel");

  MFloat totalWorkload = F0;
  for(MInt i = 0; i < m_noPartitionCells; ++i) {
    totalWorkload += get<2>(m_partitionCellList[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, &totalWorkload, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "totalWorkload");

  // e. determine maxNoCPUs
  MLong maxNoOffsprings = 0;
  MFloat maxWorkload = F0;
  for(MInt i = 0; i < m_noPartitionCells; ++i) {
    if(get<1>(m_partitionCellList[i]) + 1 > maxNoOffsprings) maxNoOffsprings = get<1>(m_partitionCellList[i]) + 1;
    maxWorkload = mMax(maxWorkload, get<2>(m_partitionCellList[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE, &maxNoOffsprings, 1, MPI_LONG, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "maxNoOffsprings");
  MPI_Allreduce(MPI_IN_PLACE, &maxWorkload, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "maxWorkload");
  // MFloat maxNoCPUs = 0.2 * (MFloat) m_noTotalCells / (MFloat) maxNoOffsprings;
  // m_log << "     + maximum number of possible CPUs with an average deviation of max. 20%: " << maxNoCPUs << endl;
  // cout << "     + maximum number of possible CPUs with an average deviation of max. 20%: " << maxNoCPUs << endl;
  // m_log << "       * obtainable average number of cells: " << m_noTotalCells / (MLong)mMax(F1,maxNoCPUs) << "\n
  // (tunable by reducing the property partitionCellOffspringThreshold)" << endl; cout << "       * obtainable average
  // number of cells: " << m_noTotalCells / (MLong)mMax(F1,maxNoCPUs) << "\n         (tunable by reducing the
  // property partitionCellOffspringThreshold)" << endl;
  MFloat maxNoCPUs = totalWorkload / maxWorkload;

  MLong noTotalMinLevelCells = 0;
  MLong minLevelCellOffset = 0;
  MLong partitionCellOffset = 0;
  MLong noLeafCells = 0;
  MInt maxPartitionLevelShift = 0;
  MInt noMinLevelCells = 0;
  vector<pair<MLong, MInt>> minLevelCells;
  MIntScratchSpace isMinLevelCell(m_noCells, AT_, "isMinLevelCell");
  isMinLevelCell.fill(0);
  for(MInt i = 0; i < m_noCells; ++i) {
    if(a_level(i) == minLevel) {
      minLevelCells.push_back(make_pair(a_globalId(i), i));
      isMinLevelCell(i) = 1;
      noMinLevelCells++;
    }
    if(a_noChildren(i) == 0) noLeafCells++;
  }
  MPI_Allreduce(MPI_IN_PLACE, &noLeafCells, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "noLeafCells");
  sort(minLevelCells.begin(), minLevelCells.end());

  if(noDomains() == 1) {
    m_noTotalCells = (MLong)m_noCells;
    m_noTotalPartitionCells = (MLong)m_noPartitionCells;
    noTotalMinLevelCells = (MLong)noMinLevelCells;
  } else {
    mAlloc(m_noPartitionCellsPerDomain, noDomains(), "m_noPartitionCellsPerDomain", AT_);
    MIntScratchSpace noMinLevelCellsPerDomain(noDomains(), AT_, "noMinLevelCellsPerDomain");
    MPI_Allgather(&m_noPartitionCells, 1, MPI_INT, m_noPartitionCellsPerDomain, 1, MPI_INT, mpiComm(), AT_,
                  "m_noPartitionCells", "m_noPartitionCellsPerDomain");
    MPI_Allgather(&noMinLevelCells, 1, MPI_INT, noMinLevelCellsPerDomain.getPointer(), 1, MPI_INT, mpiComm(), AT_,
                  "noMinLevelCells", "noMinLevelCellsPerDomain.getPointer()");

    for(MInt d = 0; d < globalDomainId(); d++) {
      partitionCellOffset += (MLong)m_noPartitionCellsPerDomain[d];
      minLevelCellOffset += (MLong)noMinLevelCellsPerDomain[d];
    }

    m_noTotalPartitionCells = 0;
    noTotalMinLevelCells = 0;
    for(MInt d = 0; d < noDomains(); d++) {
      m_noTotalPartitionCells += (MLong)m_noPartitionCellsPerDomain[d];
      noTotalMinLevelCells += (MLong)noMinLevelCellsPerDomain[d];
    }
  }
  MFloat avgWorkload = totalWorkload / ((MFloat)m_noTotalPartitionCells);
  MFloat avgOffspring = ((MFloat)m_noTotalCells) / ((MFloat)m_noTotalPartitionCells);

  map<MLong, MInt> globalToLocal;
  for(MInt i = 0; i < m_noCells; ++i) {
    globalToLocal.insert(pair<MLong, MInt>(a_globalId(i), i));
  }

  set<MLong> partitionLevelAncestorIds;
  MInt partitionLevelShift = 0;
  for(MInt i = 0; i < m_noPartitionCells; ++i) {
    MInt levelDiff = a_level(get<0>(m_partitionCellList[i])) - m_minLevel;
    partitionLevelShift = mMax(levelDiff, partitionLevelShift);
  }
  if(partitionLevelShift) {
    for(MInt i = 0; i < m_noPartitionCells; ++i) {
      MInt cellId = get<0>(m_partitionCellList[i]);
      // if ( a_parentId( globalToLocal[a_globalId(cellId)] ) < 0 ) continue;
      MLong parentId = a_parentId(globalToLocal[a_globalId(cellId)]);
      // MLong parentId = a_globalId(cellId);
      while(parentId > -1) {
        if(parentId >= m_cellOffsetPar && parentId < m_cellOffsetPar + m_noCells) {
          partitionLevelAncestorIds.insert(parentId);
          parentId = a_parentId(globalToLocal[parentId]);
        } else {
          parentId = -1;
        }
      }
    }
  }

  MLong totalNoPartitionLevelAncestors = 0;
  MLong localPartitionLevelAncestorCount = (signed)partitionLevelAncestorIds.size();
  MPI_Allreduce(&localPartitionLevelAncestorCount, &totalNoPartitionLevelAncestors, 1, MPI_LONG, MPI_SUM, mpiComm(),
                AT_, "localPartitionLevelAncestorCount", "totalNoPartitionLevelAncestors");
  MPI_Allreduce(&partitionLevelShift, &maxPartitionLevelShift, 1, MPI_INT, MPI_MAX, mpiComm(), AT_,
                "partitionLevelShift", "maxPartitionLevelShift");

  ParallelIo grid(filename, PIO_REPLACE, mpiComm());

  grid.defineArray(PIO_LONG, "partitionCellsGlobalId", m_noTotalPartitionCells);
  grid.defineArray(PIO_FLOAT, "partitionCellsWorkload", m_noTotalPartitionCells);

  grid.defineArray(PIO_LONG, "minLevelCellsTreeId", noTotalMinLevelCells);
  grid.defineArray(PIO_LONG, "minLevelCellsNghbrIds", 2 * nDim * noTotalMinLevelCells);

  grid.defineArray(PIO_UCHAR, "cellInfo", m_noTotalCells);

  if(m_noSolvers > 1 || g_multiSolverGrid) grid.defineArray(PIO_UCHAR, "solver", m_noTotalCells);

  if(m_writeCoordinatesToGridFile) {
    m_log << "NOTE: writing coordinates to grid file" << std::endl;
    for(MInt i = 0; i < nDim; i++) {
      grid.defineArray(PIO_FLOAT, "coordinates_" + std::to_string(i), m_noTotalCells);
    }
  }

  // c. prepare I/O

  MPI_Offset start = -1;
  MPI_Offset startMinLevelCells = -1;
  MPI_Offset startPartitionCells = -1;
  MPI_Offset count = -1;
  MPI_Offset countMinLevelCells = -1;
  MPI_Offset countPartitionCells = -1;

  MInt nodims = nDim;
  MInt tstep = 0;

  outStream << "g_multiSolverGrid: " << g_multiSolverGrid << endl;
  if(g_multiSolverGrid) {
    for(MInt b = 0; b < m_noSolvers; b++) {
      MLong solverCount = 0;
      // for (MInt i=0;i<m_noTotalCells; i++){
      for(MInt i = 0; i < m_noCells; i++) {
        if(a_isInSolver(i, b) == true) {
          solverCount++;
        }
      }
      cerr << "counted: " << solverCount << " for solver: " << b << endl;
      stringstream attributeName;
      attributeName << "noCells_" << b;

      MPI_Allreduce(MPI_IN_PLACE, &solverCount, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "solverCount");

      grid.setAttributes(&solverCount, attributeName.str(), 1);
    }
  }


  grid.setAttributes(&nodims, "nDim", 1);
  grid.setAttributes(&m_noSolvers, "noSolvers", 1);
  grid.setAttributes(&tstep, "globalTimeStep", 1);
  grid.setAttributes(&m_noTotalCells, "noCells", 1);
  grid.setAttributes(&noLeafCells, "noLeafCells", 1);
  grid.setAttributes(&noTotalMinLevelCells, "noMinLevelCells", 1);
  grid.setAttributes(&m_noTotalPartitionCells, "noPartitionCells", 1);
  grid.setAttributes(&totalNoPartitionLevelAncestors, "noPartitionLevelAncestors", 1);
  grid.setAttributes(&minLevel, "minLevel", 1);
  grid.setAttributes(&maxLevel, "maxLevel", 1);
  grid.setAttributes(&m_maxUniformRefinementLevel, "maxUniformRefinementLevel", 1);
  grid.setAttributes(&maxPartitionLevelShift, "maxPartitionLevelShift", 1);
  grid.setAttributes(&m_lengthOnLevel[0], "lengthLevel0", 1);
  grid.setAttributes(&m_centerOfGravity[0], "centerOfGravity", nDim);
  grid.setAttributes(&m_boundingBox[0], "boundingBox", 2 * nDim);

  // Add additional multisolver information if grid is reordered by a different Hilbert curve
  if(m_hasMultiSolverBoundingBox) {
    grid.setAttributes(&m_multiSolverLengthLevel0, "multiSolverLengthLevel0", 1);
    grid.setAttributes(&m_multiSolverMinLevel, "multiSolverMinLevel", 1);
    grid.setAttributes(&m_multiSolverCenterOfGravity[0], "multiSolverCenterOfGravity", nDim);
    grid.setAttributes(&m_multiSolverBoundingBox[0], "multiSolverBoundingBox", 2 * nDim);
  }

  grid.setAttributes(&m_reductionFactor, "reductionFactor", 1);
  grid.setAttributes(&m_decisiveDirection, "decisiveDirection", 1);
  grid.setAttributes(&totalWorkload, "totalWorkload", 1);
  grid.setAttributes(&maxWorkload, "partitionCellMaxWorkload", 1);
  grid.setAttributes(&avgWorkload, "partitionCellAverageWorkload", 1);
  grid.setAttributes(&maxNoOffsprings, "partitionCellMaxNoOffspring", 1);
  grid.setAttributes(&avgOffspring, "partitionCellAverageNoOffspring", 1);
  grid.setAttributes(&m_partitionCellWorkloadThreshold, "partitionCellWorkloadThreshold", 1);
  grid.setAttributes(&m_partitionCellOffspringThreshold, "partitionCellOffspringThreshold", 1);
  grid.setAttributes(&maxNoCPUs, "maxNoBalancedCPUs", 1);

  // Update start and count for parallel writing
  start = m_cellOffsetPar;
  startMinLevelCells = minLevelCellOffset;
  startPartitionCells = partitionCellOffset;
  count = m_noCells;
  countMinLevelCells = noMinLevelCells;
  countPartitionCells = m_noPartitionCells;


  //---
  grid.setOffset(countPartitionCells, startPartitionCells);
  //---

  { // partitionCellsGlobalIds
    MLongScratchSpace tmp(m_noPartitionCells, AT_, "tmp");
    for(MInt i = 0; i < m_noPartitionCells; ++i)
      tmp[i] = a_globalId(get<0>(m_partitionCellList[i]));
    grid.writeArray(tmp.getPointer(), "partitionCellsGlobalId");
  }

  { // partitionCellsWorkload
    MFloatScratchSpace tmp(m_noPartitionCells, AT_, "tmp");
    for(MInt i = 0; i < m_noPartitionCells; ++i)
      tmp[i] = get<2>(m_partitionCellList[i]);
    grid.writeArray(tmp.getPointer(), "partitionCellsWorkload");
  }


  //---
  grid.setOffset(countMinLevelCells, startMinLevelCells);
  //---

  if(minLevel > 21)
    mTerm(1, AT_,
          "Tree id can not be stored using 64 bits. Change minLevelCellsTreeId to 128 bits or reduce minLevel to 21 "
          "or less.");

  { // minLevelCellsTreeId
    MLongScratchSpace tmp(noMinLevelCells, AT_, "tmp");

    // Use multiSolver information if given
    const MLong tmpMinLevel = (m_multiSolverMinLevel > -1) ? m_multiSolverMinLevel : minLevel;
    const MFloat tmpLength0 = (m_multiSolverMinLevel > -1) ? m_multiSolverLengthLevel0 : m_lengthOnLevel[0];
    const MFloat* tmpCog = (m_multiSolverMinLevel > -1) ? &m_multiSolverCenterOfGravity[0] : m_centerOfGravity;

    for(MInt i = 0; i < noMinLevelCells; ++i) {
      maia::grid::hilbert::coordinatesToTreeId<nDim>(tmp[i], &a_coordinate(minLevelCells[i].second, 0), tmpMinLevel,
                                                     tmpCog, tmpLength0);
    }
    grid.writeArray(tmp.getPointer(), "minLevelCellsTreeId");
  }

  { // minLevelCellsNghbrIds
    MLongScratchSpace tmp(noMinLevelCells, 2 * nDim, AT_, "tmp");
    for(MInt i = 0; i < noMinLevelCells; ++i) {
      for(MInt j = 0; j < m_noNeighbors; j++) {
        tmp(i, j) = a_neighborId(minLevelCells[i].second, j);
      }
    }
    grid.setOffset(2 * nDim * countMinLevelCells, 2 * nDim * startMinLevelCells);
    grid.writeArray(tmp.getPointer(), "minLevelCellsNghbrIds");
  }

  //---
  grid.setOffset(count, start);
  //---

  {
    // cellInfo
    ScratchSpace<MUchar> tmp(m_noCells, AT_, "tmp");
    for(MInt i = 0; i < m_noCells; ++i) {
      MLong sortedLocalId = a_globalId(i) - m_cellOffsetPar;
      MUint noChilds = (MUint)a_noChildren(i);
      MUint isMinLvl = (MUint)isMinLevelCell(i);
      MUint position = 0;
      if(a_parentId(i) > -1) {
        MInt parentId = globalToLocal[a_parentId(i)];
        for(MUint j = 0; j < (unsigned)m_maxNoChildren; j++) {
          if(a_childId(parentId, j) == a_globalId(i)) position = j;
        }
      }
      MUint tmpBit = noChilds | (position << 4) | (isMinLvl << 7);
      tmp[sortedLocalId] = static_cast<MUchar>(tmpBit);
    }
    grid.writeArray(tmp.begin(), "cellInfo");

    // solverAffiliation
    if(m_noSolvers > 1 || g_multiSolverGrid) {
      for(MInt i = 0; i < m_noCells; ++i) {
        MLong sortedLocalId = a_globalId(i) - m_cellOffsetPar;
        MUint tmpBit = 0;
        for(MInt solver = 0; solver < m_noSolvers; solver++) {
          if(a_isInSolver(i, solver)) {
            tmpBit |= (1 << solver);
          }
        }
        tmp[sortedLocalId] = static_cast<MUchar>(tmpBit);
      }
      grid.writeArray(tmp.begin(), "solver");
    }
  }

  if(m_writeCoordinatesToGridFile) {
    ScratchSpace<MFloat> tmp(m_noCells, AT_, "tmp");
    for(MInt d = 0; d < nDim; d++) {
      for(MInt i = 0; i < m_noCells; ++i) {
        const MLong sortedLocalId = a_globalId(i) - m_cellOffsetPar;
        tmp[sortedLocalId] = a_coordinate(i, d);
      }
      grid.writeArray(tmp.begin(), "coordinates_" + std::to_string(d));
    }
  }

  m_log << "     + grid file written to file '" << filename << "'" << endl;
  outStream << "     + grid file written to '" << filename << "'" << endl;

  RECORD_TIMER_STOP(m_t_saveGrid);
}


#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic pop
#endif


//##################################################################################################################
//##################################################################################################################
//####                                          DEBUGGING FUNCTIONS                                             ####
//##################################################################################################################
//##################################################################################################################

/** \brief degugging function writing the grid in serial
 *
 * \author Andreas Lintermann
 * \date 05.12.2013
 *
 * Writes the grid in serial per domain.
 *
 * \param[in] level_ the level to write
 * \param[in] tag a tag which is appended to the file name
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::saveGridDomain(MInt level_, MInt tag) {
  TRACE();

  using namespace maia::parallel_io;

  MChar buf0[10];
  MChar buf[10];
  MChar buf1[10];
  MString preName;

  sprintf(buf0, "%d", tag);
  sprintf(buf, "%d", globalDomainId());
  sprintf(buf1, "%d", level_);
  preName = buf0;
  preName.append("G");
  preName.append(buf);
  preName.append("_L");
  preName.append(buf1);
  preName.append("_");
  preName.append(m_gridOutputFileName);

  size_t size = m_noCells + m_noTotalHaloCells;

  size_t sndBuf = size;
  MIntScratchSpace noCellsPerDomain(noDomains(), AT_, "noCellsPerDomain");
  MPI_Allgather(&sndBuf, 1, MPI_INT, noCellsPerDomain.getPointer(), 1, MPI_INT, mpiComm(), AT_, "sndBuf",
                "noCellsPerDomain.getPointer()");

  size_t myOffset = 0;
  for(MInt d = 0; d < globalDomainId(); d++)
    myOffset += noCellsPerDomain[d];

  // WARNING: untested switch from NetCDF/Parallel netCDF to ParallelIo
  // The method previously used direct I/O calls, which were replaced by
  // ParallelIo methods in summer 2015. However, since the method was not
  // used by any of the testcases, this code is still *untested*. Thus,
  // if your code uses this part of the code, please make sure that the
  // I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  TERMM(1, "untested I/O method, please see comment for how to proceed");
  ParallelIo parallelIo(preName.c_str(), maia::parallel_io::PIO_REPLACE, MPI_COMM_SELF);

  CartesianNetcdf CN;

  MIntScratchSpace variablesVar(33, AT_, "variablesVar");

  MString varnames[33] = {"partitionCellsId",
                          "partitionCellsNoOffsprings",
                          "partitionCellsWorkLoad",
                          "partitionCellsLvlDiff",
                          "parentId",
                          "noChildIds",
                          "childIds_0",
                          "childIds_1",
                          "childIds_2",
                          "childIds_3",
                          "childIds_4",
                          "childIds_5",
                          "childIds_6",
                          "childIds_7",
                          "level_0",
                          "level_1",
                          "level_2",
                          "noNghbrIds_0",
                          "noNghbrIds_1",
                          "noNghbrIds_2",
                          "noNghbrIds_3",
                          "noNghbrIds_4",
                          "noNghbrIds_5",
                          "nghbrIds_0",
                          "nghbrIds_1",
                          "nghbrIds_2",
                          "nghbrIds_3",
                          "nghbrIds_4",
                          "nghbrIds_5",
                          "weight",
                          "coordinates_0",
                          "coordinates_1",
                          "coordinates_2"};

  MInt type[33] = {PIO_INT, PIO_INT, PIO_INT, PIO_INT,   PIO_INT,   PIO_INT,  PIO_INT, PIO_INT, PIO_INT,
                   PIO_INT, PIO_INT, PIO_INT, PIO_INT,   PIO_INT,   PIO_INT,  PIO_INT, PIO_INT, PIO_INT,
                   PIO_INT, PIO_INT, PIO_INT, PIO_INT,   PIO_INT,   PIO_INT,  PIO_INT, PIO_INT, PIO_INT,
                   PIO_INT, PIO_INT, PIO_INT, PIO_FLOAT, PIO_FLOAT, PIO_FLOAT};

  parallelIo.defineScalar(PIO_INT, "noCells");

  for(MInt i = 0; i < 33; i++) {
    parallelIo.defineArray(type[i], varnames[i].c_str(), size);
  }
  parallelIo.setOffset(size, 0);

  MInt s = size;
  parallelIo.writeScalar(s, "noCells");

  MInt shift = m_haloCellOffsetsLevel[level_][0] - m_noCells;

  MInt cnt = 0;
  // partitionCellsIdDim
  {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = i + myOffset;
    if(noDomains() > 1)
      for(MInt j = m_haloCellOffsetsLevel[level_][0]; j < m_maxNoCells - 1; j++, i++)
        tmpVars[i] = j - shift + myOffset;
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // partitionCellsIdOffsprings,partitionCellsWorkload,partitionCellsLvlDiff
  for(MInt j = 0; j < 3; j++) {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = 1;
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = 1;
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // parentId
  {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++) {
      if(a_parentId(i) >= 0)
        tmpVars[i] = (MInt)a_parentId(i) + myOffset;
      else
        tmpVars[i] = (MInt)a_parentId(i);
    }
    if(noDomains() > 1) {
      for(MInt j = m_haloCellOffsetsLevel[level_][0]; j < m_maxNoCells - 1; j++, i++) {
        if(a_parentId(j) >= 0)
          tmpVars[i] = (MInt)a_parentId(j) - shift + myOffset;
        else
          tmpVars[i] = (MInt)a_parentId(j);
      }
    }
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // noChildIds
  {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = 0;
    if(noDomains() > 1)
      for(MInt j = m_haloCellOffsetsLevel[level_][0]; j < m_maxNoCells - 1; j++, i++)
        tmpVars[i] = 0;
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // childIds
  for(MInt j = 0; j < 8; j++) {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++) {
      if(a_childId(i, j) >= 0)
        tmpVars[i] = (MInt)a_childId(i, j) + myOffset;
      else
        tmpVars[i] = (MInt)a_childId(i, j);
    }
    if(noDomains() > 1) {
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++) {
        if(a_childId(k, j) >= 0)
          tmpVars[i] = (MInt)a_childId(k, j) - shift + myOffset;
        else
          tmpVars[i] = (MInt)a_childId(k, j);
      }
    }
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // level
  for(MInt j = 0; j < 3; j++) {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = a_level(i);
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = a_level(k);
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // noNghbrIds
  for(MInt j = 0; j < 6; j++) {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = (a_neighborId(i, j)) >= 0 ? 1 : 0;
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = (a_neighborId(k, j)) >= 0 ? 1 : 0;
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // nghbrIds
  for(MInt j = 0; j < 6; j++) {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    if(noDomains() > 1) {
      for(; i < m_noCells; i++) {
        if(a_neighborId(i, j) >= 0)
          if(a_neighborId(i, j) >= m_noCells)
            tmpVars[i] = (MInt)a_neighborId(i, j) - shift; // + myOffset;
          else
            tmpVars[i] = (MInt)a_neighborId(i, j); // + myOffset;
        else
          tmpVars[i] = (MInt)a_neighborId(i, j);
      }
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++) {
        if(a_neighborId(k, j) >= 0)
          if(a_neighborId(k, j) >= m_noCells)
            tmpVars[i] = (MInt)a_neighborId(k, j) - shift; // + myOffset;
          else
            tmpVars[i] = (MInt)a_neighborId(k, j); // + myOffset;
        else
          tmpVars[i] = (MInt)a_neighborId(k, j);
      }
    } else {
      for(; i < m_noCells; i++)
        if(a_neighborId(i, j) >= m_noCells)
          tmpVars[i] = (MInt)a_neighborId(i, j) + myOffset;
        else
          tmpVars[i] = (MInt)a_neighborId(i, j);
    }
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // weight
  {
    MIntScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = 1;
    if(noDomains() > 1)
      for(MInt j = m_haloCellOffsetsLevel[level_][0]; j < m_maxNoCells - 1; j++, i++)
        tmpVars[i] = 1;
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // coordinates
  for(MInt j = 0; j < 3; j++) {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = a_coordinate(i, j);
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = a_coordinate(k, j);
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }
}

/** \brief debugging function writing debug info to a file
 *
 * \author Andreas Lintermann
 * \date 05.12.2013
 *
 * Writes certain information to files in serial.
 *
 * \param[in] level_ the level to write
 * \param[in] tag a tag which is appended to the file name
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::writeGridInformation(MInt level_, MInt tag) {
  TRACE();

  using namespace maia::parallel_io;

  outStream << "Writing grid info" << endl;
  MChar buf0[10];
  MChar buf[10];
  MChar buf1[10];
  MString preName;
  MString preName1;

  sprintf(buf0, "%d", tag);
  sprintf(buf, "%d", globalDomainId());
  sprintf(buf1, "%d", level_);
  preName = buf0;
  preName.append("G");
  preName.append(buf);
  preName.append("_L");
  preName.append(buf1);
  preName.append("_");
  preName.append(m_gridOutputFileName);

  preName1 = buf0;
  preName1.append("Info");
  preName1.append(buf);
  preName1.append("_L");
  preName1.append(buf1);
  preName1.append("_");
  preName1.append(m_gridOutputFileName);

  size_t size = m_noCells + m_noTotalHaloCells;

  // WARNING: untested switch from NetCDF/Parallel netCDF to ParallelIo
  // The method previously used direct I/O calls, which were replaced by
  // ParallelIo methods in summer 2015. However, since the method was not
  // used by any of the testcases, this code is still *untested*. Thus,
  // if your code uses this part of the code, please make sure that the
  // I/O still works as expected and then remove this warning as well as
  // the subsequent TERMM().
  TERMM(1, "untested I/O method, please see comment for how to proceed");
  ParallelIo parallelIo(preName1.c_str(), maia::parallel_io::PIO_REPLACE, MPI_COMM_SELF);

  parallelIo.setAttribute(preName.c_str(), "gridFile");

  CartesianNetcdf CN;

  MString varnames[5] = {"domainId", "haloCode", "windowCode", "bndCut", "inside"};
  MInt type[5] = {PIO_INT, PIO_INT, PIO_INT, PIO_INT, PIO_INT};

  parallelIo.defineScalar(PIO_INT, "noCells");

  for(MInt i = 0; i < 5; i++) {
    parallelIo.defineArray(type[i], varnames[i].c_str(), size);
  }

  MInt cnt = 0;
  MInt s = size;
  parallelIo.writeScalar(s, "noCells");

  parallelIo.setOffset(size, 0);

  vector<vector<MInt>> winCellIdsPerDomain;
  vector<vector<MInt>> haloCellIdsPerDomain;
  for(MInt d = 0; d < m_noNeighborDomains; d++) {
    vector<MInt> c;
    winCellIdsPerDomain.push_back(c);
    vector<MInt> e;
    haloCellIdsPerDomain.push_back(e);
  }

  // find the halo and window cells for all domains
  for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
    // clear all cells
    for(MInt i = 0; i < m_noCells; i++)
      a_hasProperty(i, 3) = 0;

    // collect the halo cells and mark the window cells
    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      MInt lev_pos = 2 * (lev - m_minLevel);
      for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++) {
        haloCellIdsPerDomain[dom].push_back(p);
        MLong* neigh = &a_neighborId(p, 0);
        for(MInt n = 0; n < m_noNeighbors; n++)
          if(neigh[n] < m_noCells && neigh[n] > -1) {
            a_hasProperty(neigh[n], 3) = 1;
          }
      }
    }

    // collect the window cells
    for(MInt lev = m_minLevel; lev <= m_maxRfnmntLvl; lev++) {
      for(MInt i = m_levelOffsets[lev][0]; i < m_levelOffsets[lev][1]; i++)
        if(a_hasProperty(i, 3)) {
          winCellIdsPerDomain[dom].push_back(i);
        }
    }
  }

  // domainId
  {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = globalDomainId();
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = globalDomainId();
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // haloCells
  {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = -1;
    if(noDomains() > 1)
      for(MInt lev = level_; lev >= m_minLevel; lev--)
        for(MInt dom = 0; dom < m_noNeighborDomains; dom++) {
          MInt lev_pos = 2 * (lev - m_minLevel);
          for(MInt p = m_haloCellOffsets[dom][lev_pos]; p < m_haloCellOffsets[dom][lev_pos + 1]; p++, i++)
            tmpVars[i] = m_neighborDomains[dom];
        }
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // windowCells
  {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++) {
      vector<MInt> doms;
      for(MInt dom = 0; dom < m_noNeighborDomains; dom++)
        for(MInt k = 0; k < (signed)winCellIdsPerDomain[dom].size(); k++)
          if(i == winCellIdsPerDomain[dom][k]) doms.push_back(m_neighborDomains[dom]);

      if(doms.size() > 0) {
        MInt code = 0;
        for(MInt l = 0; l < (signed)doms.size(); l++)
          code += doms[l] * pow(10, l);
        tmpVars[i] = code;
      } else
        tmpVars[i] = -1;
    }

    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = -2;

    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // bndCut
  {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = a_hasProperty(i, 1);
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = a_hasProperty(k, 1);
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }

  // inside
  {
    MFloatScratchSpace tmpVars(size, AT_, "tmpVars");
    MInt i = 0;
    for(; i < m_noCells; i++)
      tmpVars[i] = a_hasProperty(i, 0);
    if(noDomains() > 1)
      for(MInt k = m_haloCellOffsetsLevel[level_][0]; k < m_maxNoCells - 1; k++, i++)
        tmpVars[i] = a_hasProperty(k, 0);
    parallelIo.writeArray(tmpVars.getPointer(), varnames[cnt].c_str());
    cnt++;
  }
}

/** \brief writes a pseudo solution file with some grid information
 *
 * \author Andreas Lintermann
 * \date 15.10.2013
 *
 * This simply puts some data into the cell data, so that the result can be shown in ParaView.
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::writeGridInformationPar() {
  TRACE();

  outStream << "WRITING INFO FILE" << endl;

  // the grid is already finalized for output. The cells are somehow
  //   mixed up and need to be maped ( smart c&p from saveGrid )

  // File creation.
  using namespace maia::parallel_io;
  //  ParallelIo grid ("Info.HDF5", PIO_REPLACE, mpiComm());
  ParallelIo grid("out/Info.Netcdf", PIO_REPLACE, mpiComm());

  // Creating file header.
  grid.defineScalar(PIO_INT, "noCells");
  grid.defineArray(PIO_INT, "rfnDistance", m_noTotalCells);
  grid.defineArray(PIO_INT, "bndCut", m_noTotalCells);
  grid.defineArray(PIO_INT, "domainId", m_noTotalCells);

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "solver_" << solver;
    grid.defineArray(PIO_INT, ss.str(), m_noTotalCells);
  }

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "solverBoundary_" << solver;
    grid.defineArray(PIO_INT, ss.str(), m_noTotalCells);
  }

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "toRefine_" << solver;
    grid.defineArray(PIO_INT, ss.str(), m_noTotalCells);
  }
  /*
    for(MInt dim = 0; dim < nDim; dim++){
      stringstream ss;
      ss << "coord_" << dim;
      grid.defineArray(PIO_FLOAT, ss.str(),m_noTotalCells);
    }
  */

  grid.setAttribute(m_gridOutputFileName, "gridFile");

  MPI_Offset start = -1;
  MPI_Offset count = -1;

  MInt* tmpScratchInt = 0;

  // noCells.
  start = 0;
  count = 1;
  grid.setOffset(count, start);

  grid.writeScalar(m_noTotalCells, "noCells");

  // Update start and count for parallel writing
  if(noDomains() == 1)
    start = 0;
  else
    start = m_cellOffsetPar;
  count = m_noCells;

  grid.setOffset(count, start);

  MIntScratchSpace tmp(m_noCells, AT_, "tmp");
  tmpScratchInt = tmp.getPointer();

  for(MInt j = 0; j < m_noCells; ++j) {
    MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
    tmpScratchInt[sortedLocalId] = a_refinementDistance(j);
  }
  grid.writeArray(tmpScratchInt, "rfnDistance");

  for(MInt j = 0; j < m_noCells; ++j) {
    MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
    tmpScratchInt[sortedLocalId] = a_hasProperty(j, 1);
  }
  grid.writeArray(tmpScratchInt, "bndCut");

  for(MInt j = 0; j < m_noCells; ++j) {
    MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
    tmpScratchInt[sortedLocalId] = globalDomainId();
  }
  grid.writeArray(tmpScratchInt, "domainId");

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "solver_" << solver;
    for(MInt j = 0; j < m_noCells; ++j) {
      MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
      tmpScratchInt[sortedLocalId] = a_isInSolver(j, solver);
    }
    grid.writeArray(tmpScratchInt, ss.str());
  }

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "solverBoundary_" << solver;
    for(MInt j = 0; j < m_noCells; ++j) {
      MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
      tmpScratchInt[sortedLocalId] = a_isSolverBoundary(j, solver);
    }
    grid.writeArray(tmpScratchInt, ss.str());
  }

  for(MInt solver = 0; solver < m_noSolvers; solver++) {
    stringstream ss;
    ss << "toRefine_" << solver;
    for(MInt j = 0; j < m_noCells; ++j) {
      MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
      tmpScratchInt[sortedLocalId] = a_isToRefineForSolver(j, solver);
    }
    grid.writeArray(tmpScratchInt, ss.str());
  }
  /*
    for(MInt dim = 0; dim < nDim; dim++){
      MFloatScratchSpace tmpF(m_noCells, AT_, "tmpF");
      MFloat* tmpScratchF = tmpF.getPointer();
      stringstream ss;
      ss << "coord_" << dim;
      for (MInt j = 0; j < m_noCells; ++j){
        MLong sortedLocalId = a_globalId(j) - m_cellOffsetPar;
        tmpScratchF[sortedLocalId] = a_coordinate(j,dim);
      }
      grid.writeArray(tmpScratchF, ss.str());
    }
  */
}

/** \brief checks if the neighborhood is correct
 *
 * \author Andreas Lintermann
 * \date 15.10.2013
 *
 * Finds minimal and maximal cells and walks in each direction
 * and checks the neighborhood by comparing the coordinates
 *
 * \param[in] level_ the level to apply this algorithm to
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::checkNeighborhood(MInt level_) {
  TRACE();

  outStream << "    * performing consistency check" << endl;

  std::array<MFloat, nDim> min;
  std::array<MFloat, nDim> max;
  for(MInt d = 0; d < nDim; d++) {
    min[d] = m_centerOfGravity[d] + 0.5 * (m_reductionFactor * m_geometryExtents[m_decisiveDirection]);
    max[d] = m_centerOfGravity[d] - 0.5 * (m_reductionFactor * m_geometryExtents[m_decisiveDirection]);
  }

  // find minimum coordinates
  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    MFloat* c = &a_coordinate(i, 0);

    for(MInt d = 0; d < nDim; d++) {
      if(c[d] < min[d]) min[d] = c[d];
      if(c[d] > max[d]) max[d] = c[d];
    }
  }


  vector<vector<MInt>> min_cells;
  vector<vector<MInt>> max_cells;
  // find all cells
  for(MInt d = 0; d < nDim; d++) {
    vector<MInt> min_d;
    vector<MInt> max_d;
    for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
      MFloat* c = &a_coordinate(i, 0);

      if(approx(c[d], min[d], MFloatEps)) min_d.push_back(i);

      if(approx(c[d], max[d], MFloatEps)) max_d.push_back(i);
    }

    min_cells.push_back(min_d);
    max_cells.push_back(max_d);
  }

  MFloat eps = 0.00001;

  // do the check
  for(MInt d = 0; d < 1; d++) {
    outStream << d << endl;
    MInt n_dir = 2 * d + 1;
    for(MInt i = 0; i < (MInt)min_cells[d].size(); i++) {
      MInt start_id = min_cells[d][i];

      MFloat* last_c = &a_coordinate(start_id, 0);

      outStream << start_id << ": " << last_c[0] << " " << last_c[1] << " " << last_c[2] << endl;

      if(a_neighborId(start_id, n_dir) >= 0) {
        MInt next_id = a_neighborId(start_id, n_dir);
        MFloat* next_c = &a_coordinate(next_id, 0);

        while(next_id != start_id) {
          if(a_neighborId(next_id, n_dir) >= 0) {
            outStream << next_id << ": " << next_c[0] << " " << next_c[1] << " " << next_c[2] << endl;
            outStream << a_hasProperty(next_id, 1) << " ";
            // test
            MInt od1 = (d + 1) % nDim;
            MInt od2 = (d + 2) % nDim;

            if(!(fabs(last_c[od1] - next_c[od1]) < eps) || !(fabs(last_c[od2] - next_c[od2]) < eps))
              outStream << "ERROR1 " << endl;

            if(!((fabs(last_c[d] - next_c[d]) - m_lengthOnLevel[level_]) < eps))
              outStream << "ERROR2 " << fabs(last_c[d] - next_c[d] - m_lengthOnLevel[level_]) << " "
                        << m_lengthOnLevel[level_] << endl;

            last_c = &a_coordinate(next_id, 0);
            next_id = a_neighborId(next_id, n_dir);
            next_c = &a_coordinate(next_id, 0);
          } else {
            break;
          }
        }
        outStream << endl;
      }
      outStream << endl;
    }
  }


  // do the check
  for(MInt d = 0; d < nDim; d++) {
    // cout << d << endl;
    MInt n_dir = 2 * d;
    for(MInt i = 0; i < (MInt)max_cells[d].size(); i++) {
      MInt start_id = max_cells[d][i];

      MFloat* last_c = &a_coordinate(start_id, 0);

      // cout << start_id << ": " << last_c[0] << " " << last_c[1] << " " << last_c[2] << endl;

      if(a_neighborId(start_id, n_dir) >= 0) {
        MInt next_id = a_neighborId(start_id, n_dir);
        MFloat* next_c = &a_coordinate(next_id, 0);

        while(next_id != start_id) {
          if(a_neighborId(next_id, n_dir) >= 0) {
            // cout << next_id << ": " << next_c[0] << " " << next_c[1] << " " << next_c[2] << endl;
            // test
            MInt od1 = (d + 1) % nDim;
            MInt od2 = (d + 2) % nDim;

            if(!(fabs(last_c[od1] - next_c[od1]) < eps) || !(fabs(last_c[od2] - next_c[od2]) < eps))
              outStream << "ERROR1 " << endl;

            if(!((fabs(last_c[d] - next_c[d]) - m_lengthOnLevel[level_]) < eps))
              outStream << "ERROR2 " << fabs(last_c[d] - next_c[d] - m_lengthOnLevel[level_]) << " "
                        << m_lengthOnLevel[level_] << endl;

            last_c = &a_coordinate(next_id, 0);
            next_id = a_neighborId(next_id, n_dir);
            next_c = &a_coordinate(next_id, 0);
          } else {
            break;
          }
        }
      }
      // cout << endl;
    }
  }
}

/** \brief checks if the distance between neighboring cells is alright
 *
 * \author Andreas Lintermann
 * \date 02.12.2013
 *
 * \param[in] level_ the level to check
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::checkNeighborhoodDistance(MInt level_) {
  TRACE();

  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    MLong* neigh = &a_neighborId(i, 0);
    for(MInt n = 0; n < m_noNeighbors; n++)
      if(neigh[n] >= 0) {
        MFloat dist = 0.0;
        for(MInt dim = 0; dim < nDim; dim++)
          dist += (a_coordinate(neigh[n], dim) - a_coordinate(i, dim))
                  * (a_coordinate(neigh[n], dim) - a_coordinate(i, dim));

        dist = sqrt(dist);
        if(!approx(dist, m_lengthOnLevel[level_], MFloatEps)) outStream << dist << endl;
      }
  }

  for(MInt i = m_haloCellOffsetsLevel[level_][0]; i < m_haloCellOffsetsLevel[level_][1]; i++) {
    MLong* neigh = &a_neighborId(i, 0);
    for(MInt n = 0; n < m_noNeighbors; n++)
      if(neigh[n] >= 0) {
        MFloat dist = 0.0;
        for(MInt dim = 0; dim < nDim; dim++)
          dist += (a_coordinate(neigh[n], dim) - a_coordinate(i, dim))
                  * (a_coordinate(neigh[n], dim) - a_coordinate(i, dim));

        dist = sqrt(dist);
        if(!approx(dist, m_lengthOnLevel[level_], MFloatEps)) {
          outStream << i << " " << n << " " << neigh[n] << " " << dist << " " << m_lengthOnLevel[level_]
                    << a_coordinate(i, 0) << " " << a_coordinate(i, 1) << " " << a_coordinate(i, 2) << " "
                    << a_coordinate(neigh[n], 0) << " " << a_coordinate(neigh[n], 1) << " " << a_coordinate(neigh[n], 2)
                    << " " << a_level(i) << " " << a_level(neigh[n]) << endl;
        }
      }
  }
}

/** \brief checks if the neighborhood relation of the cells is alright
 *
 * \author Andreas Lintermann
 * \date 02.12.2013
 *
 * \param[in] level_ the level to check
 *
 **/
template <MInt nDim>
void GridgenPar<nDim>::checkNeighborhoodIntegrity(MInt level_) {
  TRACE();

  outStream << "Checking neighborhood integrity for level: " << level_ << endl;
  outStream << "Normal:" << endl;
  for(MInt i = m_levelOffsets[level_][0]; i < m_levelOffsets[level_][1]; i++) {
    MLong* neigh = &a_neighborId(i, 0);
    for(MInt n = 0; n < m_noNeighbors; n++)
      if(neigh[n] >= 0)
        if(a_neighborId(neigh[n], oppositeDirGrid[n]) != i)
          // cout << "ERROR: " << i << " is not a neighbor of " << neigh[n] << " " << n << " " << oppositeDirGrid[n]
          outStream << "ERROR: " << i << " is not a neighbor of " << neigh[n] << " " << n << " " << oppositeDirGrid[n]
                    << endl;
  }

  // cout << "Halos" << endl;
  outStream << "Halos" << endl;
  for(MInt i = m_haloCellOffsetsLevel[level_][0]; i < m_haloCellOffsetsLevel[level_][1]; i++) {
    MLong* neigh = &a_neighborId(i, 0);
    for(MInt n = 0; n < m_noNeighbors; n++)
      if(neigh[n] >= 0)
        if(a_neighborId(neigh[n], oppositeDirGrid[n]) != i)
          // cout << "ERROR: " << i << " is not a neighbor of " << neigh[n] << " " << n << " " << oppositeDirGrid[n]
          outStream << "ERROR: " << i << " is not a neighbor of " << neigh[n] << " " << n << " " << oppositeDirGrid[n]
                    << endl;
  }
}
