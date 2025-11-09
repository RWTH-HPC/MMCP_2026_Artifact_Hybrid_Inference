// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSLICES_H_
#define DGSLICES_H_

#include <numeric>
#include "COMM/mpioverride.h"
#include "typetraits.h"

// Predeclaration of auxiliary claas
class DgSliceSeries;
template <MInt nDim, class SysEqn>
class DgCartesianSolver;

/** \brief Determine all coordinates and alloc buffer size create a valid 2D grid which represents a
 *         slice from a 3D grid.
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 06.01.2016
 *
 **/
template <MInt nDim, class SysEqn>
class DgSlices {
 public:
  explicit DgSlices(DgCartesianSolver<nDim, SysEqn>& solver) : m_solver(solver) {}
  void setProperties();
  void init();
  void save(const MBool isFinalTimeStep);
  MBool enabled() const { return m_enabled; }

 private:
  MPI_Comm mpiComm() const { return m_solver.mpiComm(); }
  MInt solverId() const { return m_solver.m_solverId; }
  MInt domainId() const { return m_solver.domainId(); }

  void updateGridFile(DgSliceSeries& sliceSeries, const MString& prefix = "");

  void init(DgSliceSeries& sliceSeries);
  void save(DgSliceSeries& sliceSeries, const MBool isFinalTimeStep);

  // Dimension var to ensure that slice is created only from 3D to 2D
  static constexpr const MInt m_dim = 3;

  // Store access to DG solver
  DgCartesianSolver<nDim, SysEqn>& m_solver;
  // Enables functionality of this class
  MBool m_enabled = false;
  // Holds all properties and buffers for each slice
  std::vector<DgSliceSeries> m_sliceSeries{};

  /// List of slice variables
  std::vector<MInt> m_sliceVarIds{};
  /// Number of variables for each slice variable
  std::vector<MInt> m_noSliceVars{};
  /// List of variable names for each slice variable
  std::vector<std::vector<MString>> m_sliceVarNames{};
};


// auxiliary class contains to store slice information
class DgSliceSeries {
 public:
  DgSliceSeries(const MString& axis, const MFloat intercept, const MInt writeInterval, const MString& gridFileName,
                const MInt fileNo)
    : m_axis(axis),
      m_intercept(intercept),
      m_writeInterval(writeInterval),
      m_gridFileName(gridFileName),
      m_fileNo(fileNo) {}

  // Checks for initialization
  MBool m_isInit = false;

  // Settings of slice
  // Axis which represents the normal vector of the slice layer
  MString m_axis = "";
  // Coordinate of slice
  MFloat m_intercept = 0.0;
  // Number of how many files will be created
  MInt m_writeInterval = -1;
  // Name of the 2D slice grid file
  MString m_gridFileName = "";
  // Number of slice
  const MInt m_fileNo = -1;

  // MPI Comm for current slice
  MPI_Comm m_mpiComm = MPI_COMM_NULL;

  // Buffer
  // Number of hilbertIds in current slice
  MInt m_noHilbertIds = 0;
  // Number of local nodes (also include nodes for dummy elements)
  MInt m_localNoNodes = 0;
  // Number of cells belong to slice
  MInt m_noCellInSlice = -1;
  // Holds coordinates of elements which are relevant to the slice
  std::vector<MFloat> m_coordinates{};
  // Store polyDeg of every slice cell
  std::vector<MInt> m_polyDegs{};
  // Number of cells per HilbertId (for file writing)
  std::vector<MInt> m_cellIdsPerHilbertId{};
  // Offset for polyDeg computed by slice hilbertId
  std::vector<MInt> m_cellIdsOffset{};
  // elementIds in the slice are needed to determine data points
  std::vector<MInt> m_elementIds{};
  // holds position in internal array for each element
  std::vector<MInt> m_elementOffset{};
  // holds number of nodes for each element
  std::vector<MInt> m_elementNodes{};
  // Number of elements per slice hilbertId (for file writing)
  std::vector<MInt> m_elementIdsPerHilbertId{};
  // Number of nodes per slice hilbetrtId (number of local data)
  std::vector<MInt> m_elementNodesPerHilbertIds{};
  // Offset for local nodes computed by slice hilbertId
  std::vector<MInt> m_elementNodesHilbertOffset{};
};


/** \fn void DgSlices<nDim, SysEqn>::setProperties()
 * \brief Read properties for all slices and create corresponding slice objects.
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 06.01.2016
 *
 **/
template <MInt nDim, class SysEqn>
void DgSlices<nDim, SysEqn>::setProperties() {
  TRACE();

  using namespace std;

  if(Context::propertyExists("sliceEnabled", solverId())) {
    // Read properties and save them
    /*! \property
      \page propertyPageDG DG
      \section sliceEnabled
      <code>MString DgSlice::m_enabled</code>\n
      default = false
      This enables the dg slice feature
      Possible values are:
      <ul>
        <li>true, false</li>
      </ul>
      Keywords:
      <i>DISCONTINUOUS_GALERKIN, I/O</i>
    */
    m_enabled = Context::getSolverProperty<MBool>("sliceEnabled", solverId(), AT_);

    // Leave when not enabled
    if(!m_enabled) {
      return;
    }
  }

  // TODO labels:DG,PP unify properties for DG-slices and PP-slices!
  m_log << "Set properties for DG Slice class...";

  // Check for default properties

  /*! \property
    \page propertyPageDG DG
    \section sliceGridFileName
    <code>MString DgSlice::m_gridFileName</code>\n
    default = grid.slice
    Name of the 2D slice grid file, which will be created automatically
    Possible values are:
    <ul>
      <li>any String</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  MString gridFileNameDefault = "grid.slice";
  gridFileNameDefault = Context::getSolverProperty<MString>("sliceGridFileName", solverId(), AT_, &gridFileNameDefault);

  /*! \property
    \page propertyPageDG DG
    \section sliceAxis
    <code>MString DgSlice::m_axis</code>\n
    Specifies the axis which represents the normal vector of the slice layer.\n
    Possible values are:
    <ul>
      <li>'x', 'y' or 'z'</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  MString axisDefault;
  if(Context::propertyExists("sliceAxis", solverId())) {
    axisDefault = Context::getSolverProperty<MString>("sliceAxis", solverId(), AT_, &axisDefault);
  }

  /*! \property
    \page propertyPageDG DG
    \section sliceIntercept
    <code>MFloat DgSlice::m_intercept</code>\n
    Gives the coordinate where the slice layer cuts m_axis.\n
    Possible values are:
    <ul>
      <li>any float</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  MFloat interceptDefault = 0.0;
  interceptDefault = Context::getSolverProperty<MFloat>("sliceIntercept", solverId(), AT_, &interceptDefault);

  /*! \property
    \page propertyPageDG DG
    \section sliceWriteInterval
    <code>MString DgSlice::m_writeInterval</code>\n
    default = <code>restartInterval</code> if specified\n
    Specifies at which interval a slice file should be saved.\n
    Possible values are:
    <ul>
      <li>any positive integer</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  MInt writeIntervalDefault = -1;
  if(m_solver.m_restartInterval > 0) {
    writeIntervalDefault =
        Context::getSolverProperty<MInt>("sliceWriteInterval", solverId(), AT_, &m_solver.m_restartInterval);
  } else {
    writeIntervalDefault =
        Context::getSolverProperty<MInt>("sliceWriteInterval", solverId(), AT_, &writeIntervalDefault);
  }

  // read properties for all slice numbers
  MInt noFiles = 0;
  MBool sliceFound = false;

  /*! \property
    \page propertyPageDG DG
    \section sliceAxis_
    <code>-</code>\n
    Sets the slice `fileNumber` to be written out: `sliceAxis_fileNumber`.
    Keywords: <i>DISCONTINUOUS_GALERKIN, I/O</i>
  */
  // Check if sliceAxis or sliceIntercept property for next slice exists
  if(Context::propertyExists("sliceAxis_" + to_string(noFiles), solverId())
     || Context::propertyExists("sliceIntercept_" + to_string(noFiles), solverId())) {
    sliceFound = true;
  }

  // Check if a full set of default properties are set, if true write one slice file
  if(Context::propertyExists("sliceAxis", solverId()) && writeIntervalDefault > 0) {
    sliceFound = true;
  }

  // Leave method if not all properties are set
  if(!sliceFound) {
    m_enabled = false;
    m_log << "not all necessary slice properties are set!" << std::endl;
    return;
  }

  // Search for slices in properties and create slice objects
  while(sliceFound) {
    // Check for non default values of properties
    /*! \property
      \page propertyPageDG DG
      \section sliceGridFileName_N
      Corresponds to property `sliceGridFileName` when using multiple slices.
    */
    MString gridFileName = Context::getSolverProperty<MString>("sliceGridFileName_" + to_string(noFiles), solverId(),
                                                               AT_, &gridFileNameDefault);

    MString axis =
        Context::getSolverProperty<MString>("sliceAxis_" + to_string(noFiles), solverId(), AT_, &axisDefault);

    /*! \property
      \page propertyPageDG DG
      \section sliceIntercept_N
      Corresponds to property `sliceIntercept` when using multiple slices.
    */
    MFloat intercept =
        Context::getSolverProperty<MFloat>("sliceIntercept_" + to_string(noFiles), solverId(), AT_, &interceptDefault);

    /*! \property
      \page propertyPageDG DG
      \section sliceWriteInterval_N
      Corresponds to property `sliceWriteInterval` when using multiple slices.
    */
    MInt writeInterval = Context::getSolverProperty<MInt>("sliceWriteInterval_" + to_string(noFiles), solverId(), AT_,
                                                          &writeIntervalDefault);

    // Check for incorrect properties
    if(writeInterval < 1) {
      TERMM(1, "WriteInterval for Slice no. " + to_string(noFiles) + " is smaller than 1.");
    }


    if(axis != "x" && axis != "y" && axis != "z") {
      TERMM(1, "Illegal axis for Slice no. " + to_string(noFiles) + ". Only x, y or z is a valid axis");
    }

    // Add slice number to grid file name
    gridFileName += '_' + to_string(noFiles);
    // Create slice objects
    m_sliceSeries.emplace_back(axis, intercept, writeInterval, gridFileName, noFiles);
    noFiles++;

    // Check for noext slice
    sliceFound = Context::propertyExists("sliceAxis_" + to_string(noFiles), solverId())
                 || Context::propertyExists("sliceIntercept_" + to_string(noFiles), solverId());
  }

  m_sliceVarIds.clear();
  m_noSliceVars.clear();
  m_sliceVarNames.clear();
  m_solver.getSolverSamplingProperties(m_sliceVarIds, m_noSliceVars, m_sliceVarNames, "Slice");

  // @ansgar_slice TODO labels:DG not required at the moment
  // Allocate additional storage for the sampling variables if required
  // solver().initSolverSamplingVariables(m_solverSamplingVarIds, m_noSolverSamplingVars);

  // Initalization finished
  m_enabled = true;
  m_log << "done!" << std::endl;

  // Abort to avoid usage of dg slice in 2D case
  IF_CONSTEXPR(nDim == 2) {
    if(m_enabled) {
      TERMM(1, "Can not use DG Slice in 2D calculation.");
    }
  }
}


/** \brief Initialize all slice objects
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 20.06.2016
 *
 **/
template <MInt nDim, class SysEqn>
void DgSlices<nDim, SysEqn>::init() {
  TRACE();
  // Return if dg slice is deactivated
  if(!m_enabled) {
    return;
  }

  // Call init method for all slice objects
  for(MUint i = 0; i < m_sliceSeries.size(); i++) {
    init(m_sliceSeries[i]);
  }
}


/** \brief Create 2D slice grid file and get all elements which belong to specific grid.
 *
 *         First the DG Elements, Number of Nodes, PolyDegree and Coorinates, which belong the
 *         current slice are determined from the slice cells. This is done in data chunks with the
 *         same hilbertId.
 *
 *         After that a MPI Communicator is created for all domains of the current slice.
 *
 *         In the last step, the HilbertId Information is used, to determine global offsets for file
 *         creation. Since the data in the grids and solution files is sorted by their hilbertId, it
 *         is necessary to determine this offset information. It is important to understand, that
 *         the number of hilbertIds is different on every domain and that they can be distributed in
 *         in any order, since a slice of a 3D hilbert curve is performed. To determine the offset
 *         correctly the local information of hilbertIds and number of nodes is exchanged  globally
 *         (globally in the context of slices), that every domain knows how many hilbertIds and
 *         nodes (with same hilbertId) are on each domain and from this the global offsets are
 *         calculated for the local data.
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 06.01.2016
 *
 **/
template <MInt nDim, class SysEqn>
void DgSlices<nDim, SysEqn>::init(DgSliceSeries& sliceSeries) {
  TRACE();

  // Reset variables, needed when DgSlices::init() is called multiple times, e.g., DLB enabled
  // Reset noLocalNodes
  sliceSeries.m_localNoNodes = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 1: Determine DG Elements and number of Nodes
  //////////////////////////////////////////////////////////////////////////////////////////////////

  // Create grid file and get slice cells from grid
  MIntScratchSpace cellIds(m_solver.grid().noInternalCells(), AT_, "cellIds");
  MIntScratchSpace hilbertInfo(m_solver.grid().noInternalCells() * 3, AT_, "hilbertIds");
  // Get cellIds for current slice and hilbertId informaiton. The cellIds are used to determine the
  // DG elements for the slice and the hilbertId information is needed for file creation.
  const MString gridFilePath = m_solver.outputDir() + sliceSeries.m_gridFileName;
  m_solver.grid().raw().createGridSlice(sliceSeries.m_axis, sliceSeries.m_intercept, gridFilePath, solverId(),
                                        &sliceSeries.m_noCellInSlice, &cellIds[0], &sliceSeries.m_noHilbertIds,
                                        &hilbertInfo[0]);
  sliceSeries.m_elementNodesPerHilbertIds.resize(sliceSeries.m_noHilbertIds);
  sliceSeries.m_elementIdsPerHilbertId.resize(sliceSeries.m_noHilbertIds);

  m_log << "Initializing DG Slice " << sliceSeries.m_fileNo << " ...";
  // Search for slice elements if slice cells are found on current domain
  if(sliceSeries.m_noCellInSlice > 0) {
    // Convert grid cell-ids into solver cell-ids (only required if there are multiple solvers)
    if(m_solver.grid().raw().treeb().noSolvers() > 1) {
      for(MInt i = 0; i < sliceSeries.m_noCellInSlice; i++) {
        cellIds[i] = m_solver.grid().tree().grid2solver(cellIds[i]);
      }
    }

    // Initalize polyDegree vector
    sliceSeries.m_polyDegs.resize(sliceSeries.m_noCellInSlice);
    std::fill_n(&sliceSeries.m_polyDegs[0], sliceSeries.m_noCellInSlice, m_solver.m_initPolyDeg);

    // Assume size of element variables
    const MInt noMaxElements = m_solver.m_elements.size();
    // Collectors for element data
    MIntScratchSpace elementIds(noMaxElements, AT_, "elementIds");
    MIntScratchSpace elementNodes(noMaxElements, AT_, "elementNodes");
    MIntScratchSpace elementOffset(noMaxElements, AT_, "elementOffset");

    // Total number of coordinates for all elements
    MInt noCoords = 0;
    // Assume max number of nodes
    const MInt noMaxNodes = pow(
        *std::max_element(&m_solver.m_elements.polyDeg(0), &m_solver.m_elements.polyDeg(0) + m_solver.m_elements.size())
            + 1,
        2);
    MFloatScratchSpace coordinates(noMaxElements * noMaxNodes * m_dim, AT_, "coordinates");

    // Holds coordinates for one element
    std::array<MFloat, m_dim> nodeCoord;
    // This loop searches for elements and node coordinates of the slice. The cells are checked in
    // chunks the same sliceHilbertId to determine the element node offset at the end of this part.
    MInt noElementsSlice = 0;
    for(MInt h = 0, cellSum = 0, noNodes = 0, prevNoNodes = 0, sumElements = 0; h < sliceSeries.m_noHilbertIds; h++) {
      for(MInt c = cellSum, e = 0; (c - cellSum) < hilbertInfo[h * 3 + 1]; e++) {
        // Check to identify cells which are not an element; These cells are
        // counted as dummy nodes in the solution file and have the value 0

        // skip slice which do not have a corresponding element
        if(cellIds[c] < m_solver.m_elements.cellId(e)) {
          // count dummy nodes to cells which are not elements
          noNodes += ipow(m_solver.m_initPolyDeg + 1, m_dim - 1);
          // Go to next cell after counting nodes
          c++;
          // Hold element until all dummy elements are counted
          e--;
          continue;
        }
        // save node coordinates if element is in slice
        if(cellIds[c] == m_solver.m_elements.cellId(e)) {
          // Determine number of nodes
          const MInt noNodes1D = m_solver.m_elements.polyDeg(e) + 1;

          // Tensor to hold node coordinates
          const MFloatTensor nodeCoordinates(&m_solver.m_elements.nodeCoordinates(e), noNodes1D, noNodes1D, noNodes1D,
                                             m_dim);
          // Loop over all nodes to save the coordinates
          for(MInt i = 0; i < noNodes1D; i++) {
            for(MInt j = 0; j < noNodes1D; j++) {
              if(sliceSeries.m_axis == "x") {
                nodeCoord[0] = sliceSeries.m_intercept;
                nodeCoord[1] = nodeCoordinates(0, i, j, 1);
                nodeCoord[2] = nodeCoordinates(0, i, j, 2);
              } else if(sliceSeries.m_axis == "y") {
                nodeCoord[0] = nodeCoordinates(i, 0, j, 0);
                nodeCoord[1] = sliceSeries.m_intercept;
                nodeCoord[2] = nodeCoordinates(i, 0, j, 2);
              } else if(sliceSeries.m_axis == "z") {
                nodeCoord[0] = nodeCoordinates(i, j, 0, 0);
                nodeCoord[1] = nodeCoordinates(i, j, 0, 1);
                nodeCoord[2] = sliceSeries.m_intercept;
              }
              // Save coordinate tuple
              std::copy_n(&nodeCoord[0], m_dim, &coordinates[noCoords]);
              noCoords += m_dim;
            }
          }
          elementIds[noElementsSlice] = e;

          // Otherwise use current offset for the current element
          elementOffset[noElementsSlice] = noNodes - prevNoNodes; // sliceSeries.m_localNoNodes;

          // Increase offset based on element polynomial degree
          const MInt noNodesXD = ipow(m_solver.m_elements.polyDeg(e) + 1, m_dim - 1);
          noNodes += noNodesXD;

          // Set number of nodes to use later and proceed with next element
          elementNodes[noElementsSlice] = noNodesXD;
          noElementsSlice++;

          // save polyDegree (it is easier to collect this data here)
          sliceSeries.m_polyDegs[c] = m_solver.m_elements.polyDeg(e);

          c++;
        }
      }
      // Save number of elements which have the same hilbertId
      sliceSeries.m_elementIdsPerHilbertId[h] = noElementsSlice - sumElements;

      // Save number of nodes per hilbert id and total number
      sliceSeries.m_elementNodesPerHilbertIds[h] = noNodes - prevNoNodes;
      sliceSeries.m_localNoNodes += sliceSeries.m_elementNodesPerHilbertIds[h];
      // Save number of evaluated cells, elements and nodes
      cellSum += hilbertInfo[h * 3 + 1];
      prevNoNodes = noNodes;
      sumElements = noElementsSlice;
    }
    // Resize variables
    sliceSeries.m_elementIds.resize(noElementsSlice);
    sliceSeries.m_elementNodes.resize(noElementsSlice);
    sliceSeries.m_elementOffset.resize(noElementsSlice);
    sliceSeries.m_coordinates.resize(noCoords);
    // Get values from scratches and save them
    std::copy_n(&elementIds[0], noElementsSlice, &sliceSeries.m_elementIds[0]);
    std::copy_n(&elementNodes[0], noElementsSlice, &sliceSeries.m_elementNodes[0]);
    std::copy_n(&elementOffset[0], noElementsSlice, &sliceSeries.m_elementOffset[0]);
    std::copy_n(&coordinates[0], noCoords, &sliceSeries.m_coordinates[0]);
  }

  // If only non element cells are found, make m_coordinates non zero length
  // to create a dummy states scratch
  if(sliceSeries.m_coordinates.empty()) {
    sliceSeries.m_coordinates.resize(1);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 2: Create MPI Communicator
  //////////////////////////////////////////////////////////////////////////////////////////////////

  using namespace maia;
  // Set MPI Comm for ranks that are part of slice
  // Create a new MPI Communicator
  ScratchSpace<MInt> sliceScratch(globalNoDomains(), AT_, "sliceScratch");

  // Get Rank if domain holds elements in slice
  MInt rank = -1;
  if(sliceSeries.m_noCellInSlice > 0) {
    MPI_Comm_rank(mpiComm(), &rank);
  }

  // Combine all ranks to get relevant ranks
  MPI_Allgather(&rank, 1, type_traits<MInt>::mpiType(), &sliceScratch[0], 1, type_traits<MInt>::mpiType(), mpiComm(),
                AT_, "rank", "sliceScratch[0]");

  // Check for relevant ranks and save them to create the new communicator
  const MInt noSliceDomains =
      std::count_if(sliceScratch.begin(), sliceScratch.end(), [](const MInt a) { return a != -1; });
  MIntScratchSpace sliceDomains(noSliceDomains, AT_, "sliceDomains");
  MInt position = 0;
  for(MInt i : sliceScratch) {
    if(i != -1) {
      sliceDomains[position] = i;
      position++;
    }
  }

  // Create new point data mpi group
  MPI_Group globalGroup, localGroup;
  MPI_Comm_group(mpiComm(), &globalGroup, AT_, "globalGroup");
  MPI_Group_incl(globalGroup, noSliceDomains, &sliceDomains[0], &localGroup, AT_);

  // Create new communicator and clean up
  MPI_Comm_create(mpiComm(), localGroup, &sliceSeries.m_mpiComm, AT_, "sliceSeries.m_mpiComm");

  MPI_Group_free(&globalGroup, AT_);
  MPI_Group_free(&localGroup, AT_);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Step 3: Exchange hilbertId data and determine offsets
  //////////////////////////////////////////////////////////////////////////////////////////////////

  if(sliceSeries.m_noCellInSlice > 0) {
    // Get slice domain
    MInt sliceDomain = -1;
    MPI_Comm_rank(sliceSeries.m_mpiComm, &sliceDomain);

    // Determine total no of slice hilbertIds
    // Cells and elements are sorted by slice hilbertIds and also by domain. So if we have hilbertId
    // 2 on domain 3 and 5, then the order is determined by domaindId. This information is
    // distributed globally and the purpose of the next lines.
    MIntScratchSpace noLocalHilbertIds(noSliceDomains, AT_, "noLocalHilbertIds");
    noLocalHilbertIds[sliceDomain] = sliceSeries.m_noHilbertIds;

    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, &noLocalHilbertIds[0], 1, MPI_INT, sliceSeries.m_mpiComm, AT_,
                  "MPI_IN_PLACE", "noLocalHilbertIds[0]");

    // Compute offset for recieve data
    MIntScratchSpace offsetsRecvData(noSliceDomains, AT_, "offsetsRecvData");
    offsetsRecvData[0] = 0;
    for(MInt i = 1; i < noSliceDomains; i++) {
      offsetsRecvData[i] = offsetsRecvData[i - 1] + noLocalHilbertIds[i - 1] * 3;
    }
    MInt noTotalRecvData = offsetsRecvData[noSliceDomains - 1] + noLocalHilbertIds[noSliceDomains - 1] * 3;

    // Create data array for sending
    MIntScratchSpace sendHilbertData(noLocalHilbertIds[sliceDomain] * 3, AT_, "sendHilbertData");
    // Send hilbertId, domainId and number of Nodes per HilbertId
    for(MInt i = 0; i < sliceSeries.m_noHilbertIds; i++) {
      sendHilbertData[i * 3] = hilbertInfo[i * 3];
      sendHilbertData[i * 3 + 1] = domainId();
      sendHilbertData[i * 3 + 2] = sliceSeries.m_elementNodesPerHilbertIds[i];
    }

    // Create data array for recieving
    MIntScratchSpace recvHilbertData(noTotalRecvData, AT_, "recvHilbertData");
    MIntScratchSpace noRecvData(noSliceDomains, AT_, "noRecvData");
    for(MInt i = 0; i < noSliceDomains; i++) {
      noRecvData[i] = noLocalHilbertIds[i] * 3;
    }

    // Exhange elementNoHilbertNodes
    MPI_Allgatherv(&sendHilbertData[0], sendHilbertData.size(), MPI_INT, &recvHilbertData[0], &noRecvData[0],
                   &offsetsRecvData[0], MPI_INT, sliceSeries.m_mpiComm, AT_, "sendHilbertData[0]",
                   "recvHilbertData[0]");

    // Create sliceGlobalHilbertInfo
    ScratchSpace<std::array<MInt, 3>> sliceGlobalHilbertInfo(noTotalRecvData / 3, AT_, "sliceGlobalHilbertInfo");
    for(MInt i = 0; i < noTotalRecvData / 3; i++) {
      // Contains number of cells sorted by hilbertId (first), domaindId (second) and number of
      // nodes per hilbertId
      sliceGlobalHilbertInfo[i][0] = recvHilbertData[i * 3];
      sliceGlobalHilbertInfo[i][1] = recvHilbertData[i * 3 + 1];
      sliceGlobalHilbertInfo[i][2] = recvHilbertData[i * 3 + 2];
    }
    // Sort sliceGlobalHilbertInfo by domainId...
    std::stable_sort(sliceGlobalHilbertInfo.begin(), sliceGlobalHilbertInfo.end(),
                     [](const std::array<MInt, 3>& a, const std::array<MInt, 3>& b) { return a[1] < b[1]; });
    // ... and then sort stable by hilbertId -> list ordered by hilbertId and then nested for each hilbertId by domainId
    std::stable_sort(sliceGlobalHilbertInfo.begin(), sliceGlobalHilbertInfo.end(),
                     [](const std::array<MInt, 3>& a, const std::array<MInt, 3>& b) { return a[0] < b[0]; });

    // Calulate elementNodeOffset by sliceHilbertIds for file writing
    MInt maxNoHilbertIds = *std::max_element(noLocalHilbertIds.begin(), noLocalHilbertIds.end());
    sliceSeries.m_elementNodesHilbertOffset.resize(maxNoHilbertIds);
    sliceSeries.m_cellIdsPerHilbertId.resize(maxNoHilbertIds);
    sliceSeries.m_cellIdsOffset.resize(maxNoHilbertIds);
    // The data in sliceGlobalHilbertInfo is sorted by their hilbertId on the first level. On the
    // second level by domainId. Example:
    //                             [hilbertId, domaindId, noNodes]
    // sliceGlobalHilbertInfo[0] = [0, 0, 20]
    // sliceGlobalHilbertInfo[1] = [0, 2, 10]
    // sliceGlobalHilbertInfo[2] = [1, 1, 6]
    //
    // To determine the global offset for the local data (for each local hilbertId) count the number
    // of nodes of all hilbertIds below the current one. If this hilbertId is distributed on
    // multiple domains, then count also the next number of nodes till the local domainid is reached
    // in sliceGlobalHilbertInfo. In the example above the offset for hilbertId 0 on domain 1 would
    // be 20.
    for(MInt i = 0, j = 0; i < noLocalHilbertIds[sliceDomain]; i++) {
      MInt offset = 0;
      j = 0;
      // Count until current hilbertId is reached
      while(sliceGlobalHilbertInfo[j][0] < hilbertInfo[i * 3]) {
        offset += sliceGlobalHilbertInfo[j][2];
        j++;
      }
      // Count to current domain
      while(sliceGlobalHilbertInfo[j][1] < domainId()) {
        offset += sliceGlobalHilbertInfo[j][2];
        j++;
      }
      // Set local node offset
      sliceSeries.m_elementNodesHilbertOffset[i] = offset;
      // save noCells with same hilbertId  for writing poly degree to the grid
      sliceSeries.m_cellIdsPerHilbertId[i] = hilbertInfo[i * 3 + 1];
      // save offset for writing poly degree to the grid
      sliceSeries.m_cellIdsOffset[i] = hilbertInfo[i * 3 + 2];
    }
  }

  // TODO labels:DG,DOC
  MBool optimizedSliceIo = true;
  optimizedSliceIo = Context::getBasicProperty<MBool>("optimizedSliceIo", AT_, &optimizedSliceIo);

  // Reorganize data before writing
  if(optimizedSliceIo && sliceSeries.m_noCellInSlice > 0) {
    // Create buffers to organize m_elementOffset
    std::vector<MInt> bufferElementIdsPerHilbertId;
    std::vector<MInt> bufferElementNodesPerHilbertIds;
    // Fill buffers
    std::copy(sliceSeries.m_elementIdsPerHilbertId.begin(),
              sliceSeries.m_elementIdsPerHilbertId.end(),
              back_inserter(bufferElementIdsPerHilbertId));
    std::copy(sliceSeries.m_elementNodesPerHilbertIds.begin(),
              sliceSeries.m_elementNodesPerHilbertIds.end(),
              back_inserter(bufferElementNodesPerHilbertIds));

    MInt noHilbertIdChunks = 0;

    // Reconstruct data for writing
    for(MInt h = 1, i = 0, offsetSum = 0; h <= (MInt)sliceSeries.m_elementNodesHilbertOffset.size(); h++) {
      // Reorganize data
      if(h < sliceSeries.m_noHilbertIds) {
        if((sliceSeries.m_elementNodesPerHilbertIds[i] + sliceSeries.m_elementNodesHilbertOffset[i])
           == sliceSeries.m_elementNodesHilbertOffset[h]) {
          sliceSeries.m_elementNodesPerHilbertIds[i] += sliceSeries.m_elementNodesPerHilbertIds[h];
          sliceSeries.m_elementIdsPerHilbertId[i] += sliceSeries.m_elementIdsPerHilbertId[h];

          sliceSeries.m_cellIdsPerHilbertId[i] += sliceSeries.m_cellIdsPerHilbertId[h];

          // Update element offsets of the elements belonging to the added hilbert id
          MInt startIdx = std::accumulate(&bufferElementIdsPerHilbertId[0], &bufferElementIdsPerHilbertId[h], 0);
          MInt endIdx = std::accumulate(&bufferElementIdsPerHilbertId[0], &bufferElementIdsPerHilbertId[h + 1], 0);
          offsetSum += bufferElementNodesPerHilbertIds[h - 1];
          for(MInt j = startIdx; j < endIdx; j++) {
            sliceSeries.m_elementOffset[j] += offsetSum;
          }
        } else {
          i++;
          sliceSeries.m_elementNodesPerHilbertIds[i] = sliceSeries.m_elementNodesPerHilbertIds[h];
          sliceSeries.m_elementNodesHilbertOffset[i] = sliceSeries.m_elementNodesHilbertOffset[h];
          sliceSeries.m_elementIdsPerHilbertId[i] = sliceSeries.m_elementIdsPerHilbertId[h];

          sliceSeries.m_cellIdsPerHilbertId[i] = sliceSeries.m_cellIdsPerHilbertId[h];
          sliceSeries.m_cellIdsOffset[i] = sliceSeries.m_cellIdsOffset[h];

          // Reset offsetSum in case of multiple concatenations
          offsetSum = 0;
        }
      }

      // Set remaining elements to zero at last iteration
      if(h == (MInt)sliceSeries.m_elementNodesHilbertOffset.size()) {
        i++;
        noHilbertIdChunks = i;

        std::fill(sliceSeries.m_elementNodesPerHilbertIds.begin() + i, sliceSeries.m_elementNodesPerHilbertIds.end(),
                  0);
        std::fill(sliceSeries.m_elementNodesHilbertOffset.begin() + i, sliceSeries.m_elementNodesHilbertOffset.end(),
                  0);
        std::fill(sliceSeries.m_elementIdsPerHilbertId.begin() + i, sliceSeries.m_elementIdsPerHilbertId.end(), 0);

        std::fill(sliceSeries.m_cellIdsPerHilbertId.begin() + i, sliceSeries.m_cellIdsPerHilbertId.end(), 0);
        std::fill(sliceSeries.m_cellIdsOffset.begin() + i, sliceSeries.m_cellIdsOffset.end(), 0);
      }
    }

    MInt maxNoHilbertIdChunks = -1;
    MPI_Allreduce(&noHilbertIdChunks, &maxNoHilbertIdChunks, 1, maia::type_traits<MInt>::mpiType(), MPI_MAX,
                  sliceSeries.m_mpiComm, AT_, "noHilbertIdChunks", "maxNoHilbertIdChunks");
    sliceSeries.m_elementNodesHilbertOffset.resize(maxNoHilbertIdChunks);

    sliceSeries.m_noHilbertIds = noHilbertIdChunks;
    sliceSeries.m_elementNodesPerHilbertIds.resize(noHilbertIdChunks);
    sliceSeries.m_elementIdsPerHilbertId.resize(noHilbertIdChunks);
    sliceSeries.m_cellIdsPerHilbertId.resize(noHilbertIdChunks);
    sliceSeries.m_cellIdsOffset.resize(noHilbertIdChunks);
  }

  sliceSeries.m_isInit = true;
  m_log << "done!" << std::endl;
}

/** \brief Call save method for all slice objects
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 20.06.2016
 *
 **/
template <MInt nDim, class SysEqn>
void DgSlices<nDim, SysEqn>::save(const MBool isFinalTimeStep) {
  TRACE();
  // Return if dg slice is deactivated
  if(!m_enabled) {
    return;
  }

  // Call save method for all slice objects
  for(MUint i = 0; i < m_sliceSeries.size(); i++) {
    save(m_sliceSeries[i], isFinalTimeStep);
  }
}


/** \brief Collect data for all slice nodes and create solution files. Data is sorted by their slice
 *         hilbertId and therefore also written in chunks of hilbertIds.
 *
 *  \author Marcus Wiens (marcus) marcus.wiens@rwth-aachen.de
 *  \date 06.01.2016
 *
 **/
template <MInt nDim, class SysEqn>
void DgSlices<nDim, SysEqn>::save(DgSliceSeries& sliceSeries, const MBool isFinalTimeStep) {
  TRACE();

  // Ensure that initalization was done
  if(!sliceSeries.m_isInit) {
    TERMM(1, "DG Slice class was not properly initialized.");
  }

  // Leave when domain is not part of slice
  if(sliceSeries.m_localNoNodes < 1) {
    return;
  }

  // Check for write Interval
  if(!(m_solver.m_timeStep % sliceSeries.m_writeInterval == 0 || isFinalTimeStep)) {
    return;
  }

  // TODO labels:DG @ansgar_slice not required at the moment
  /* // Compute sampling variables if not done already for another time series in this time step */
  /* if(m_solverCalcSamplingVars) { */
  /*   solver().calcSamplingVariables(m_solverSamplingVarIds); */
  /*   m_solverCalcSamplingVars = false; */
  /* } */

  // Write all slice variables to the same file
  const MInt noVars = std::accumulate(m_noSliceVars.begin(), m_noSliceVars.end(), 0);
  const MUint noVarIds = m_noSliceVars.size();

  // Get state at slice coordinates
  MFloatTensor coordinatesTensor(&sliceSeries.m_coordinates[0],
                                 std::ceil((MFloat)sliceSeries.m_coordinates.size() / m_dim), m_dim);

  // Save states of nodes
  MFloatScratchSpace states(std::ceil((MFloat)sliceSeries.m_coordinates.size() / m_dim) * noVars, AT_, "states");
  MFloatTensor stateTensor(&states[0], std::ceil((MFloat)sliceSeries.m_coordinates.size() / m_dim), noVars);
  // TODO labels:DG,totest check performance of state calculation (timer)
  // Calculate states of all nodes
  if(sliceSeries.m_coordinates.size() > 1) {
    for(MUint e = 0, offset = 0; e < sliceSeries.m_elementNodes.size(); e++) {
      const MUint noNodes = sliceSeries.m_elementNodes[e];
      for(MUint n = 0; n < noNodes; n++) {
        MInt varOffset = 0;
        for(MUint v = 0; v < noVarIds; v++) {
          m_solver.calcSamplingVarAtPoint(&coordinatesTensor(offset + n, 0), sliceSeries.m_elementIds[e],
                                          m_sliceVarIds[v], &stateTensor(offset + n, varOffset), true);
          varOffset += m_noSliceVars[v];
        }
      }
      offset += noNodes;
    }
  }

  // Write File
  using namespace maia::parallel_io;

  // TODO labels:DG slice naming with slice axis and intercept? same as in postprocessing
  // Create output file name with slices axis and number
  // TODO labels:DG change bX (block X) into e.g. sX (solver X)? need to change other output files as well
  const MString s = (m_solver.grid().raw().treeb().noSolvers() > 1) ? "_b" + std::to_string(solverId()) : "";
  std::ostringstream fileName;
  fileName << m_solver.outputDir() << "slice" << s << "_" << sliceSeries.m_fileNo << "_" << std::setw(8)
           << std::setfill('0') << m_solver.m_timeStep << ParallelIo::fileExt();
  ParallelIo file(fileName.str(), PIO_REPLACE, sliceSeries.m_mpiComm);

  // Determine offset and global number of nodes
  ParallelIo::size_type nodesOffset, globalNoNodes, offset, totalCount;
  ParallelIo::calcOffset(sliceSeries.m_localNoNodes, &nodesOffset, &globalNoNodes, sliceSeries.m_mpiComm);
  ParallelIo::calcOffset(sliceSeries.m_noCellInSlice, &offset, &totalCount, sliceSeries.m_mpiComm);

  // Set Attributes
  // Grid file name and solver type
  file.setAttribute(sliceSeries.m_gridFileName + ParallelIo::fileExt(), "gridFile");
  file.setAttribute(solverId(), "solverId");
  file.setAttribute("DG", "solverType");
  file.setAttribute(m_solver.m_timeStep, "timeStep");
  file.setAttribute(m_solver.m_time, "time");
  // Save slice axis and intercept to identify the grid file
  file.setAttribute(sliceSeries.m_axis, "sliceAxis");
  file.setAttribute(sliceSeries.m_intercept, "sliceIntercept");

  // Define arrays in file
  // Polyomial degree
  file.defineArray(PIO_INT, "polyDegs", totalCount);

  // Get information about integration method and polynomial type
  MString dgIntegrationMethod = "DG_INTEGRATE_GAUSS";
  dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", solverId(), AT_, &dgIntegrationMethod);
  MString dgPolynomialType = "DG_POLY_LEGENDRE";
  dgPolynomialType = Context::getSolverProperty<MString>("dgPolynomialType", solverId(), AT_, &dgPolynomialType);

  // Add integration method & polynomial type to grid file as attributes
  file.setAttribute(dgIntegrationMethod, "dgIntegrationMethod", "polyDegs");
  file.setAttribute(dgPolynomialType, "dgPolynomialType", "polyDegs");

  // Set Variables
  for(MUint v = 0, totalNoVars = 0; v < noVarIds; v++) {
    for(MInt i = 0; i < m_noSliceVars[v]; i++) {
      const MString name = "variables" + std::to_string(totalNoVars);
      file.defineArray(PIO_FLOAT, name, globalNoNodes);
      file.setAttribute(m_sliceVarNames[v][i], "name", name);
      totalNoVars++;
    }
  }

  // Write poly degree to file (sorted by hilbert id)
  for(MInt h = 0, pOffset = 0; h < (MInt)sliceSeries.m_elementNodesHilbertOffset.size(); h++) {
    if(h < sliceSeries.m_noHilbertIds && sliceSeries.m_cellIdsPerHilbertId[h] > 0) {
      file.setOffset(sliceSeries.m_cellIdsPerHilbertId[h], sliceSeries.m_cellIdsOffset[h]);
      file.writeArray(&sliceSeries.m_polyDegs[0 + pOffset], "polyDegs");
      pOffset += sliceSeries.m_cellIdsPerHilbertId[h];
    } else {
      file.setOffset(0, 0);
      file.writeArray(&sliceSeries.m_polyDegs[0], "polyDegs");
    }
  }

  // Number of nodes for elements per HilbertId
  const MInt tmp = std::accumulate(sliceSeries.m_elementNodes.begin(), sliceSeries.m_elementNodes.end(), 0);
  // create statetensor variable for easier accessing
  MFloatTensor stateTensorFinal(&states[0], std::max(tmp, 1), noVars);

  // Write data
  for(MInt i = 0; i < noVars; i++) {
    for(MInt h = 0, elementOffset = 0, offsetElementNodes = 0; h < (MInt)sliceSeries.m_elementNodesHilbertOffset.size();
        h++) {
      if(h < sliceSeries.m_noHilbertIds && sliceSeries.m_elementNodesPerHilbertIds[h] > 0) {
        // Write data for every (chunk of continuous) hilbertId(s) on this domain

        // Set amount and offset for current hilbertId
        file.setOffset(sliceSeries.m_elementNodesPerHilbertIds[h], sliceSeries.m_elementNodesHilbertOffset[h]);
        // Create buffer to collect data in correct order for saving
        MFloatScratchSpace buffer(sliceSeries.m_elementNodesPerHilbertIds[h], AT_, "buffer");
        std::fill(buffer.begin(), buffer.end(), 0.0);
        const MInt noElementsSlice = sliceSeries.m_elementIdsPerHilbertId[h] + elementOffset;

        // Fill the buffer at positions where the element data belongs (reminder: elements only
        // exist on highest refinement level)
        for(MInt e = elementOffset; e < noElementsSlice; e++) {
          MFloat* const b = &buffer[sliceSeries.m_elementOffset[e]];
          for(MInt j = 0; j < sliceSeries.m_elementNodes[e]; j++) {
            b[j] = stateTensorFinal(offsetElementNodes + j, i);
          }
          offsetElementNodes += sliceSeries.m_elementNodes[e];
        }
        elementOffset = noElementsSlice;

        const MString name = "variables" + std::to_string(i);
        // write data
        file.writeArray(&buffer[0], name);
      } else {
        // Dummy calls for data writing if this domain is finished, but other domains not
        file.setOffset(0, 0);
        const MFloat tmpBuf = 0;
        const MString name = "variables" + std::to_string(i);
        file.writeArray(&tmpBuf, name);
      }
    }
  }
}
#endif
