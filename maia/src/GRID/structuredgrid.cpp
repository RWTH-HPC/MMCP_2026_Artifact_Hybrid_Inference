// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "structuredgrid.h"
#include <numeric>
#include <vector>
#include "COMM/mpioverride.h"
#include "FV/fvstructuredcell.h"
#include "FV/fvstructuredsolver.h"
#include "GRID/structuredpartition.h"
#include "INCLUDE/maiaconstants.h"
#include "IO/parallelio.h"
#include "UTIL/parallelfor.h"
#include "globals.h"
#include "globalvariables.h"

using namespace std;

/**
 * \brief Constructor for structured grids
 * \param[in] solverId_ id of the solver
 * \param[in] comm MPI Communicator
 */
template <MInt nDim>
StructuredGrid<nDim>::StructuredGrid(const MInt solverId_, const MPI_Comm comm)
  : m_solverId(solverId_), m_mpiComm(comm) {
  MPI_Comm_rank(mpiComm(), &m_domainId);
  MPI_Comm_size(mpiComm(), &m_noDomains);

  /*! \property
    \page propertiesFVSTRCTRD
    \section gridInputFileName
    <code>MString StructuredGrid::m_gridInputFileName </code>\n
    default = <code>""</code>\n \n
    Name of the grid file.\n
    Keywords: <i>GRID, STRUCTURED</i>
  */
  m_gridInputFileName = Context::getSolverProperty<MString>("gridInputFileName", m_solverId, AT_);
  m_noGhostLayers = Context::getSolverProperty<MInt>("noGhostLayers", m_solverId, AT_);

  mAlloc(m_nPoints, nDim, "m_nPoints", -1, AT_);
  mAlloc(m_nActivePoints, nDim, "m_nActivePoints", -1, AT_);
  mAlloc(m_nCells, nDim, "m_nCells", -1, AT_);
  mAlloc(m_nActiveCells, nDim, "m_nActiveCells", -1, AT_);
  mAlloc(m_nOffsetCells, nDim, "m_nOffsetCells", -1, AT_);
  mAlloc(m_nOffsetPoints, nDim, "m_nOffsetPoints", -1, AT_);
  mAlloc(m_nBlockCells, nDim, "m_nBlockCells", -1, AT_);

  mAlloc(m_periodicDisplacements, nDim * nDim, "m_periodicDisplacement", F0, AT_);

  m_totalNoCells = 0;
  m_noPoints = 1;
  m_noActiveCells = 1;
  m_noCells = 1;

  mAlloc(m_singularity, 30, "m_singularity", AT_);
}


/**
 * \brief Create decomposition of the grid into partitions for MPI parallelization
 * \param[in] readFromFile Set to true if decomposition should be read from external file
 */
template <MInt nDim>
void StructuredGrid<nDim>::gridDecomposition(MBool readFromFile) {
  if(m_domainId == 0) {
    cout << "Doing  block decomposition..." << endl;
  }
  /*! \property
    \page propertiesFVSTRCTRD
    \section readPartitionFromFile
    <code>MBool StructuredGrid::m_readDecompositionFromFile </code>\n
    default = <code>0</code>\n \n
    Trigger the use to read the MPI partitioning from a file (faster).\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>PARALLEL, PARTITIONING, STRUCTURED</i>
  */
  m_readDecompositionFromFile = false;
  m_readDecompositionFromFile =
      Context::getSolverProperty<MBool>("readPartitionFromFile", m_solverId, AT_, &m_readDecompositionFromFile);

  ParallelIoHdf5 pio(m_gridInputFileName, maia::parallel_io::PIO_READ, m_mpiComm);
  // unique grid identifier to associate grid and solution file on restart
  MInt noBlocksType = pio.getAttributeType("noBlocks", "");
  if(noBlocksType == 1) {
    pio.getAttribute<MInt>(&m_noBlocks, "noBlocks", "");
  } else if(noBlocksType == 0) {
    MFloat noBlocksFloat = -1.0;
    pio.getAttribute<MFloat>(&noBlocksFloat, "noBlocks", "");
    m_noBlocks = (MInt)noBlocksFloat;
  }

  pio.getAttribute(&m_uID, "UID", "");

  m_partition = make_unique<StructuredDecomposition<nDim>>(m_noBlocks, m_gridInputFileName, m_mpiComm, m_noDomains);

  if(readFromFile) {
    m_log << "reading repartition from file ....." << endl;
    MBool success = m_partition->readFromFile();
    if(!success) {
      m_log << "..... Reading from partition file FAILED --> new decomposition activated" << endl;
      m_partition->decompose();
    }
    m_log << "..... Reading in successful" << endl;
  } else {
    m_partition->decompose();
  }

  m_totalNoCells = 0;
  for(MInt i = 0; i < m_noBlocks; i++) {
    MLong temp = 1;
    for(MInt dim = 0; dim < nDim; dim++) {
      temp *= getBlockNoCells(i, dim);
    }
    if(temp != 1) m_totalNoCells += temp;
  }
  if(m_domainId == 0) {
    cout << "Doing  block decomposition... SUCCESSFUL!" << endl;
  }
}

/**
 * \brief Prepares the arrays containing the size of the grid (points/cells)
 *        before the grid coordinates are actually read in
 */
template <MInt nDim>
void StructuredGrid<nDim>::prepareReadGrid() {
  m_blockId = getMyBlockId();
  for(MInt i = 0; i < nDim; i++) {
    m_nOffsetCells[i] = getMyOffset(i);
    if(getMyOffset(i) > 0) {
      m_nOffsetPoints[i] = getMyOffset(i) + 1;
    } else {
      m_nOffsetPoints[i] = 0;
    }
    m_nActivePoints[i] = getMyActivePoints(i);
    m_nPoints[i] = m_nActivePoints[i] + 2 * m_noGhostLayers;

    m_nActiveCells[i] = m_nActivePoints[i] - 1;
    m_nCells[i] = m_nPoints[i] - 1;
    m_noPoints *= m_nPoints[i];
    m_noActiveCells *= m_nActiveCells[i];
    m_noCells *= m_nCells[i];
  }

  mAlloc(m_cells->coordinates, nDim, m_noCells, "m_cells->coordinates", AT_);
  mAlloc(m_coordinates, nDim, m_noPoints, "m_coordinates", -1.01010101, AT_);

  /*! \property
    \page propertiesFVSTRCTRD
    \section movingGrid
    <code>MInt FvStructuredSolver::m_movingGrid </code>\n
    default = <code> 0 </code>\n \n
    Trigger to use moving grid methods\n
    Possible values are:\n
    <ul>
    <li>0: off</li>
    <li>1: on</li>
    </ul>
    Keywords: <i>MOVING, STRUCTURED</i>
  */

  m_movingGrid = false;
  if(Context::propertyExists("movingGrid", m_solverId)) {
    m_movingGrid = Context::getSolverProperty<MBool>("movingGrid", m_solverId, AT_, &m_movingGrid);
  }

  if(m_movingGrid) {
    mAlloc(m_velocity, nDim, m_noPoints, "m_velocity", F0, AT_);
    mAlloc(m_acceleration, nDim, m_noPoints, "m_acceleration", F0, AT_);
    mAlloc(m_oldCoordinates, nDim, m_noPoints, "m_mgOldCoordinates", F0, AT_);
    mAlloc(m_initCoordinates, nDim, m_noPoints, "m_mgInitCoordinates", F0, AT_);
  }
}

/**
 * \brief Allocates memory for the metrics and the Jacobians
 */
template <MInt nDim>
void StructuredGrid<nDim>::allocateMetricsAndJacobians() {
  mAlloc(m_cells->cornerJac, m_noCells, "m_cells->cornerJac", 123.123, AT_);
  mAlloc(m_cells->cellJac, m_noCells, "m_cells->cellJac", 123.123, AT_);
  mAlloc(m_cells->oldCellJac, m_noCells, "m_cells->oldCellJac", 123.123, AT_);
  IF_CONSTEXPR(nDim == 2) mAlloc(m_cells->surfJac, m_noCells * 2, "m_cells->surfJac", AT_);
  mAlloc(m_cells->cellMetrics, nDim * nDim, m_noCells, "cellMetrics", 123.123, AT_);
  mAlloc(m_cells->cornerMetrics, nDim * nDim, m_noCells, "cornerMetrics", 123.123, AT_);
  mAlloc(m_cells->surfaceMetrics, nDim * nDim, m_noCells, "surfaceMetrics", 123.123, AT_);

  // By now only available in 2D
  IF_CONSTEXPR(nDim == 2) {
    if(m_hasSingularity) {
      MInt no_cells = 0;
      for(MInt i = 0; i < m_hasSingularity; ++i) {
        // only correct for bc 6000 not for bc 4000-5000
        if(m_singularity[i].BC == -6000) {
          for(MInt j = 0; j < nDim; j++) {
            const MInt len = m_singularity[i].end[j] - m_singularity[i].start[j];
            ASSERT(len == 1, "");
          }
          no_cells += 1;
        }
      }

      const MInt no_metrics = 4;

      mAlloc(m_cells->surfaceMetricsSingularity, no_metrics, no_cells, "surfaceMetrics", 123.123, AT_);
    }
  }
}


/**
 * \brief Reads in the coordinates (x,y,z) from the grid file
 */
template <MInt nDim>
void StructuredGrid<nDim>::readGrid() {
  ParallelIoHdf5 pio(m_gridInputFileName, maia::parallel_io::PIO_READ, m_mpiComm);

  // create the string to contain the datasetname in the file
  MString sBlockName = "/block";
  stringstream dummy1;
  dummy1 << m_blockId << "/";
  sBlockName += dummy1.str();

  MString varNames[] = {"x", "y", "z"};
  ParallelIo::size_type offset[3];
  ParallelIo::size_type size[3];
  for(MInt i = 0; i < nDim; ++i) {
    offset[i] = m_nOffsetCells[i];
    size[i] = m_nActivePoints[i];
  }

  for(MInt dim = 0; dim < nDim; dim++) {
    pio.readArray(m_coordinates[dim], sBlockName, varNames[dim], nDim, offset, size);
  }

  moveCellPoints(); // shifts points and cells in the respective array to account for the ghost cells
}


/**
 * \brief Saves coordinates for partitioned grid with ghost points.
 *        Useful for debugging.
 * \author Pascal Meysonnat
 */
template <MInt nDim>
void StructuredGrid<nDim>::writePartitionedGrid() {
  TRACE();
  // first every process needs to create the datasets and the structure where the
  // data will be stored!
  m_log << "writing the partitionedGrid.hdf5 File" << endl;
  cout << "writing the partitionedGrid.hdf5 File" << endl;
  const char* fileName = "partitionedGrid.hdf5";
  ParallelIoHdf5 pio(fileName, maia::parallel_io::PIO_REPLACE, m_mpiComm);
  MInt noDomains_ = noDomains();
  pio.setAttribute(noDomains_, "noBlocks", "");
  MString gridVarNames[3] = {"x", "y", "z"};
  for(MInt i = 0; i < noDomains(); i++) {
    // create datasets for the io library
    ParallelIo::size_type noPoints[nDim] = {};
    for(MInt j = 0; j < nDim; j++) {
      noPoints[j] = getActivePoints(i, j) + 2 * m_noGhostLayers;
    }

    stringstream path;
    path << i;
    MString partitionPathStr = "block";
    partitionPathStr += path.str();
    const char* partitionPath = partitionPathStr.c_str();
    for(MInt dim = 0; dim < nDim; ++dim) {
      pio.defineArray(maia::parallel_io::PIO_FLOAT, partitionPath, gridVarNames[dim], nDim, noPoints);
    }
  }

  // write the values into the array so that we can visualize it
  ParallelIo::size_type offset[3] = {0, 0, 0};
  ParallelIo::size_type size[3] = {m_nPoints[0], m_nPoints[1], m_nPoints[2]};
  stringstream path;
  path << domainId();
  MString partitionPathStr = "block";
  partitionPathStr += path.str();
  for(MInt dim = 0; dim < nDim; ++dim) {
    pio.writeArray(&m_coordinates[dim][0], partitionPathStr, gridVarNames[dim], nDim, offset, size);
  }
}

/**
 * \brief Moves the coordinates by the number of ghost-layers  (2D version)
          away from the boundary, this needs to be done once
          after the coordinates are read in
 */
template <>
void StructuredGrid<2>::moveCellPoints() {
  constexpr MInt nDim = 2;

  std::array<MInt, nDim> begin{0, 0};
  std::array<MInt, nDim> end{m_nActivePoints[1], m_nActivePoints[0]};

  // two dimensional case
  maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j) {
    const MInt i_org = end[0] - 1 - i;
    const MInt j_org = end[1] - 1 - j;
    const MInt i_new = i_org + m_noGhostLayers;
    const MInt j_new = j_org + m_noGhostLayers;
    const MInt pointId_org = i_org + (j_org * m_nActivePoints[1]); // position in Array
    const MInt pointId_new = i_new + (j_new * m_nPoints[1]);       // new position in Array
    for(MInt dim = 0; dim < nDim; ++dim) {
      m_coordinates[dim][pointId_new] = m_coordinates[dim][pointId_org]; // copy value to the right place
      m_coordinates[dim][pointId_org] = -1000.0;                         // for test purposes only
    }
  });
}

/**
 * \brief Moves the coordinates by the number of ghost-layers  (3D version)
          away from the boundary, this needs to be done once
          after the coordinates are read in
 */
template <>
void StructuredGrid<3>::moveCellPoints() {
  constexpr MInt nDim = 3;

  std::array<MInt, nDim> begin{0, 0};
  std::array<MInt, nDim> end{m_nActivePoints[2], m_nActivePoints[1], m_nActivePoints[0]};

  // three dimensional case
  maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j, const MInt& k) {
    const MInt i_org = end[0] - 1 - i;
    const MInt j_org = end[1] - 1 - j;
    const MInt k_org = end[2] - 1 - k;
    const MInt i_new = i_org + m_noGhostLayers;
    const MInt j_new = j_org + m_noGhostLayers;
    const MInt k_new = k_org + m_noGhostLayers;

    const MInt pointId_org = i_org + (j_org + k_org * m_nActivePoints[1]) * m_nActivePoints[2];
    const MInt pointId_new = i_new + (j_new + k_new * m_nPoints[1]) * m_nPoints[2];

    for(MInt dim = 0; dim < nDim; ++dim) {
      m_coordinates[dim][pointId_new] = m_coordinates[dim][pointId_org]; // copy value to the right place
      m_coordinates[dim][pointId_org] = F0;
    }
  });
}

/**
 * \brief Writes the current grid (including deformations for moving grids)
          to a file
 * \param[in] solutionOutput Directory to which the output is written to
 * \param[in] outputFormat File ending
 */
template <MInt nDim>
void StructuredGrid<nDim>::writeGrid(MString solutionOutput, MString outputFormat) {
  stringstream fileName;
  fileName << solutionOutput << "grid" << globalTimeStep << outputFormat;
  ParallelIoHdf5 pio(fileName.str(), maia::parallel_io::PIO_REPLACE, m_mpiComm);

  pio.setAttribute(m_noBlocks, "noBlocks", "");
  MString fileTypeName = "grid";
  pio.setAttribute(fileTypeName, "filetype", "");
  MString gridTypeName = "structured";
  pio.setAttribute(gridTypeName, "gridType", "");
  pio.setAttribute(m_uID, "UID", "");
  pio.setAttribute(globalTimeStep, "globalTimeStep", "");

  ParallelIo::size_type allPoints[3]{-1, -1, -1};
  for(MInt i = 0; i < m_noBlocks; ++i) {
    for(MInt j = 0; j < nDim; ++j) {
      allPoints[j] = getBlockNoPoints(i, j);
    }
    MString blockPathStr = "block" + std::to_string(i);
    pio.defineArray(maia::parallel_io::PIO_FLOAT, blockPathStr, "x", nDim, allPoints);
    pio.defineArray(maia::parallel_io::PIO_FLOAT, blockPathStr, "y", nDim, allPoints);
    IF_CONSTEXPR(nDim == 3) { pio.defineArray(maia::parallel_io::PIO_FLOAT, blockPathStr, "z", nDim, allPoints); }
  }

  MString blockPathStr = "block" + std::to_string(m_blockId);
  ParallelIo::size_type offset[nDim]{};
  ParallelIo::size_type size[nDim]{};
  ParallelIo::size_type ghostArray[nDim]{};
  for(MInt dim = 0; dim < nDim; ++dim) {
    offset[dim] = m_nOffsetCells[dim];
    size[dim] = m_nActivePoints[dim];
    ghostArray[dim] = m_noGhostLayers;
  }

  pio.writeArray(&m_coordinates[0][0], blockPathStr, "x", nDim, offset, size, ghostArray);
  pio.writeArray(&m_coordinates[1][0], blockPathStr, "y", nDim, offset, size, ghostArray);
  IF_CONSTEXPR(nDim == 3) { pio.writeArray(&m_coordinates[2][0], blockPathStr, "z", nDim, offset, size, ghostArray); }
}


/**
 * \brief Exchanges the boundary grid points between MPI partitions
 *
 * \author Marian Albers
 *
 * \param[in] sndComm MPI sending maps
 * \param[in] rcvComm MPI receiving maps
 * \param[in] commType Type of the exchange boundary (block boundary, periodic boundary, etc.)
 */
template <MInt nDim>
void StructuredGrid<nDim>::exchangePoints(std::vector<unique_ptr<StructuredComm<nDim>>>& sndComm,
                                          std::vector<unique_ptr<StructuredComm<nDim>>>& rcvComm,
                                          StructuredCommType commType) {
  std::vector<MPI_Request> sndRequests;
  std::vector<MPI_Request> rcvRequests;
  std::vector<MPI_Status> sndStatus;
  std::vector<MPI_Status> rcvStatus;
  sndRequests.reserve(sndComm.size());
  rcvRequests.reserve(rcvComm.size());

  gatherPoints(sndComm, commType);
  sendPoints(sndComm, commType, sndRequests);
  receivePoints(rcvComm, commType, rcvRequests);

  sndStatus.resize(sndRequests.size());
  MPI_Waitall(sndRequests.size(), &sndRequests[0], &sndStatus[0], AT_);

  rcvStatus.resize(rcvRequests.size());
  MPI_Waitall(rcvRequests.size(), &rcvRequests[0], &rcvStatus[0], AT_);

  scatterPoints(rcvComm, commType);
}

/**
 * \brief Gathers the coordinates of the points for all given sending maps
 *        and copies them to a sending buffer
 *
 * \author Marian Albers
 *
 * \param[in] sndComm MPI sending maps
 * \param[in] commType Type of the exchange boundary (block boundary, periodic boundary, etc.)
 * \return
 */
template <MInt nDim>
void StructuredGrid<nDim>::gatherPoints(std::vector<unique_ptr<StructuredComm<nDim>>>& sndComm,
                                        StructuredCommType commType) {
  TRACE();
  for(auto& snd : sndComm) {
    if(commType != snd->commType) continue;


    MFloat** coordinates = m_coordinates;
    std::array<MInt, nDim> INC{};
    for(MInt dim = 0; dim < nDim; ++dim) {
      INC[dim] = m_nPoints[dim];
    }

    MInt plusOne = 1;
    if(snd->commType == SINGULAR || snd->commType == PERIODIC_BC_SINGULAR) {
      coordinates = m_cells->coordinates;
      plusOne = 0;

      for(MInt dim = 0; dim < nDim; ++dim) {
        INC[dim] = m_nCells[dim];
      }
    }

    MBool isPeriodic = false;
    if(snd->commType == PERIODIC_BC || snd->commType == PERIODIC_BC_SINGULAR) {
      isPeriodic = true;
    }

    std::array<MInt, nDim> begin{};
    std::array<MInt, nDim> end{};
    std::array<MInt, nDim> size{};

    for(MInt dim = 0; dim < nDim; ++dim) {
      begin[dim] = snd->startInfoPoints[dim];
      end[dim] = snd->endInfoPoints[dim] + plusOne;
      size[dim] = end[dim] - begin[dim];
    }
    const MInt totalSize = std::accumulate(size.begin(), size.end(), 1, std::multiplies<double>());

    // Aliasing the unique_pointer in snd to a raw point is needed for PSTL on NVHPC
    const MInt bcId = snd->bcId;
    auto* pointBuffer = snd->pointBuffer.get();
    if constexpr(nDim == 3) {
      for(MInt dim = 0; dim < nDim; dim++) {
        maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j, const MInt& k) {
          const MInt pointId = i + (j + k * INC[1]) * INC[2];
          const MInt bufferId =
              totalSize * dim + (i - begin[0]) + ((j - begin[1]) + (k - begin[2]) * size[1]) * size[0];

          if(isPeriodic) {
            MFloat tmppoint = coordinates[dim][pointId];
            periodicPointsChange(tmppoint, bcId, dim);
            pointBuffer[bufferId] = tmppoint;
          } else {
            pointBuffer[bufferId] = coordinates[dim][pointId];
          }
        });
      }
    } else if constexpr(nDim == 2) {
      for(MInt dim = 0; dim < nDim; dim++) {
        maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j) {
          const MInt pointId = i + j * INC[1];
          const MInt bufferId = totalSize * dim + (i - begin[0]) + (j - begin[1]) * size[0];

          if(isPeriodic) {
            MFloat tmppoint = coordinates[dim][pointId];
            periodicPointsChange(tmppoint, bcId, dim);
            pointBuffer[bufferId] = tmppoint;
          } else {
            pointBuffer[bufferId] = coordinates[dim][pointId];
          }
        });
      }
    }
  }
}

/**
 * \brief Displaces the points on periodic boundaries by the distance
 *        between the two periodic boundaries
 *
 * \author Marian Albers
 *
 * \param[in] pt Coordinate point to be displaced
 * \param[in] type Periodic boundary type (4401-4406), depending on the type the displacement is added or subtracted
 * \param[in] dim Dimension (0->x,1->y,2->z) of the coordinate in which the displacement is effective
 */
template <MInt nDim>
inline void StructuredGrid<nDim>::periodicPointsChange(MFloat& pt, const MInt type, const MInt dim) {
  /*
    case 4401 4402: first periodic direction
    case 4403 4404: second periodic direction
    case 4405 4405: third periodic direction
   */

  // SND Map BC is reversed, therefore we need
  // to subtract instead of adding and vice versa compared
  // to the periodicPointsChange in fvstructuredsolverwindowinfo.cpp
  switch(type) {
    case 4401:
    case 4403:
    case 4405: {
      const MInt displacementId = (MFloat)(type - 4400 + 1) / 2.0 - 1;
      pt = pt - m_periodicDisplacements[dim * nDim + displacementId];
      break;
    }

    case 4402:
    case 4404:
    case 4406: {
      const MInt displacementId = (MFloat)(type - 4400 + 1) / 2.0 - 1;
      pt = pt + m_periodicDisplacements[dim * nDim + displacementId];
      break;
    }

    default: {
#ifndef WAR_NVHPC_PSTL
      std::cout << "ERROR!!! periodic type is wrong!!! in BC call BC: " << type << endl;
#endif
    }
  }
}

/**
 * \brief Send the coordinates between partitions to other partitions
 *
 * \author Marian Albers
 *
 * \param[in] sndComm MPI sending maps
 * \param[in] commType Type of the exchange boundary (block boundary, periodic boundary, etc.)
 * \param[in] sndRequests Stores the sending requests
 */
template <MInt nDim>
void StructuredGrid<nDim>::sendPoints(std::vector<unique_ptr<StructuredComm<nDim>>>& sndComm,
                                      StructuredCommType commType,
                                      std::vector<MPI_Request>& sndRequests) {
  TRACE();
  for(auto& snd : sndComm) {
    if(commType != snd->commType) continue;
    MPI_Request request{};
    const MInt tag = m_domainId + (snd->tagHelper) * m_noDomains;
    MInt err = MPI_Isend((void*)&snd->pointBuffer[0],
                         snd->pointBufferSize,
                         MPI_DOUBLE,
                         snd->nghbrId,
                         tag,
                         m_mpiComm,
                         &request,
                         AT_,
                         "snd->pointBuffer");
    sndRequests.push_back(request);
    if(err) cout << "rank" << m_domainId << "sending throws errors" << endl;
  }
}

/**
 * \brief Receives the coordinates between partitions to other partitions
 *
 * \author Marian Albers
 *
 * \param[in] rcvComm MPI receiving maps
 * \param[in] commType Type of the exchange boundary (block boundary, periodic boundary, etc.)
 * \param[in] rcvRequests Stores the receiving requests
 */
template <MInt nDim>
void StructuredGrid<nDim>::receivePoints(std::vector<unique_ptr<StructuredComm<nDim>>>& rcvComm,
                                         StructuredCommType commType,
                                         std::vector<MPI_Request>& rcvRequests) {
  TRACE();
  for(auto& rcv : rcvComm) {
    if(commType != rcv->commType) continue;
    MPI_Request request{};
    const MInt tag = rcv->nghbrId + (rcv->tagHelper) * m_noDomains;
    MInt err = MPI_Irecv((void*)&rcv->pointBuffer[0],
                         rcv->pointBufferSize,
                         MPI_DOUBLE,
                         rcv->nghbrId,
                         tag,
                         m_mpiComm,
                         &request,
                         AT_,
                         "rcv->pointBuffer");
    rcvRequests.push_back(request);
    if(err) cout << "rank" << m_domainId << " receiving throws errors" << endl;
  }
}

/**
 * \brief Distributes the exchanged points from the receiving buffers to the actual coordinates of the grid
 *
 * \author Marian Albers
 *
 * \param[in] rcvComm MPI receiving maps
 * \param[in] commType Type of the exchange boundary (block boundary, periodic boundary, etc.)
 */
template <MInt nDim>
void StructuredGrid<nDim>::scatterPoints(std::vector<unique_ptr<StructuredComm<nDim>>>& rcvComm,
                                         StructuredCommType commType) {
  TRACE();

  // the ordering of the grid points can be different from
  // sending instance ==> reorder it and copy it to the
  // right place
  for(auto& rcv : rcvComm) {
    if(commType != rcv->commType) continue;

    MFloat** coordinates = m_coordinates;
    std::array<MInt, nDim> INC{};
    MInt plusOne = 1;

    for(MInt dim = 0; dim < nDim; ++dim) {
      INC[dim] = m_nPoints[dim];
    }

    if(rcv->commType == SINGULAR || rcv->commType == PERIODIC_BC_SINGULAR) {
      coordinates = m_cells->coordinates;
      plusOne = 0;
      for(MInt dim = 0; dim < nDim; ++dim) {
        INC[dim] = m_nCells[dim];
      }
    }

    std::array<MInt, nDim> begin{};
    std::array<MInt, nDim> end{};
    std::array<MInt, nDim> size{};

    for(MInt dim = 0; dim < nDim; ++dim) {
      begin[dim] = rcv->startInfoPoints[dim];
      end[dim] = rcv->endInfoPoints[dim] + plusOne;
      size[dim] = end[dim] - begin[dim];
    }
    const MInt totalSize = std::accumulate(size.begin(), size.end(), 1, std::multiplies<double>());

    std::array<MInt, nDim> stepBuffer{};
    std::array<MInt, nDim> startBuffer{};
    std::array<MInt, nDim> endBuffer{};
    std::array<MInt, nDim> sizeBuffer{};

    for(MInt j = 0; j < nDim; j++) {
      stepBuffer[rcv->orderInfo[j]] = rcv->stepInfo[j];
    }

    for(MInt j = 0; j < nDim; j++) {
      endBuffer[j] = size[j] - 1;
      sizeBuffer[rcv->orderInfo[j]] = size[j];
      if(stepBuffer[j] < 0) {
        std::swap(startBuffer[j], endBuffer[j]);
      }
    }

    // Aliasing the unique_pointer in rcv to a raw point is needed for PSTL on NVHPC
    auto* pointBuffer = rcv->pointBuffer.get();
    auto* orderInfo = rcv->orderInfo.data();
    auto* startInfoCells = rcv->startInfoCells.data();
    if constexpr(nDim == 3) {
      for(MInt dim = 0; dim < nDim; dim++) {
        maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j, const MInt& k) {
          std::array<MInt, nDim> start{};
          start[orderInfo[0]] = startBuffer[0] + (i - startInfoCells[0]) * stepBuffer[0];
          start[orderInfo[1]] = startBuffer[1] + (j - startInfoCells[1]) * stepBuffer[1];
          start[orderInfo[2]] = startBuffer[2] + (k - startInfoCells[2]) * stepBuffer[2];

          const MInt bufferId = dim * totalSize + start[0] + (start[1] + start[2] * sizeBuffer[1]) * sizeBuffer[0];
          const MInt pointId = i + (j + k * INC[1]) * INC[2];

          coordinates[dim][pointId] = pointBuffer[bufferId];
        });
      }
    } else if constexpr(nDim == 2) {
      for(MInt dim = 0; dim < nDim; dim++) {
        maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j) {
          std::array<MInt, nDim> start{};
          start[orderInfo[0]] = startBuffer[0] + (i - startInfoCells[0]) * stepBuffer[0];
          start[orderInfo[1]] = startBuffer[1] + (j - startInfoCells[1]) * stepBuffer[1];

          const MInt bufferId = dim * totalSize + start[0] + start[1] * sizeBuffer[0];
          const MInt pointId = i + j * INC[1];

          coordinates[dim][pointId] = pointBuffer[bufferId];
        });
      }
    }
  }
}

/**
 * \brief Sets the reference to the cell object
 * \param[in] structuredCell Reference to the structured cell object
 */
template <MInt nDim>
void StructuredGrid<nDim>::setCellReference(StructuredCell* structuredCell) {
  m_cells = structuredCell;
}

/*=================================================
  METRICS/JACOBIAN
  =================================================*/

/**
 * \brief Computes the cell center coordinates of each cell from the surrounding grid points (2D)
 */
template <>
void StructuredGrid<2>::computeCellCenterCoordinates() {
  constexpr MInt nDim = 2;

  maia::parallelFor<true, nDim>(cellBegin(0), cellEnd(0), [=](const MInt& i, const MInt& j) {
    const MInt IJ = getPointIdFromCell(i, j);
    const MInt IP1J = getPointIdFromPoint(IJ, 1, 0);
    const MInt IJP1 = getPointIdFromPoint(IJ, 0, 1);
    const MInt IP1JP1 = getPointIdFromPoint(IJ, 1, 1);
    const MInt cellId = cellIndex(i, j);

    for(MInt dim = 0; dim < nDim; dim++) {
      // average the coordinates for cell centre data
      m_cells->coordinates[dim][cellId] =
          F1B4
          * (m_coordinates[dim][IJ] + m_coordinates[dim][IP1J] + m_coordinates[dim][IJP1] + m_coordinates[dim][IP1JP1]);
    }
  });
}

/**
 * \brief Computes the cell center coordinates of each cell from the surrounding grid points (3D)
 */
template <>
void StructuredGrid<3>::computeCellCenterCoordinates() {
  constexpr MInt nDim = 3;

  maia::parallelFor<true, nDim>(cellBegin(0), cellEnd(0), [=](const MInt& i, const MInt& j, const MInt& k) {
    const MInt pointId = i + (j + k * m_nPoints[1]) * m_nPoints[2];
    const MInt IJK = pointId;
    const MInt IP1JK = pointId + 1;
    const MInt IJP1K = pointId + m_nPoints[2];
    const MInt IP1JP1K = IJP1K + 1;
    const MInt IJKP1 = pointId + m_nPoints[2] * m_nPoints[1];
    const MInt IP1JKP1 = IJKP1 + 1;
    const MInt IJP1KP1 = pointId + m_nPoints[2] + m_nPoints[2] * m_nPoints[1];
    const MInt IP1JP1KP1 = IJP1KP1 + 1;
    const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];

    // average the coordinates for cell centre data
    for(MInt dim = 0; dim < nDim; dim++) {
      m_cells->coordinates[dim][cellId] =
          F1B8
          * (m_coordinates[dim][IJK] + m_coordinates[dim][IP1JK] + m_coordinates[dim][IJP1K]
             + m_coordinates[dim][IP1JP1K] + m_coordinates[dim][IJKP1] + m_coordinates[dim][IP1JKP1]
             + m_coordinates[dim][IJP1KP1] + m_coordinates[dim][IP1JP1KP1]);
    }
  });
}

/**
 * \brief Computes all metrics by calling the functions for each type of metric computation (cell, corner, surface)
 */
template <MInt nDim>
void StructuredGrid<nDim>::computeMetrics() {
  TRACE();
  computeSurfaceMetrics();
  computeCornerMetrics();
  computeCellMetrics();
  IF_CONSTEXPR(nDim == 2) {
    if(m_hasSingularity) computeSurfaceMetricsSingularity();
  }
}

/**
 * \brief Computes the Jacobians by calling the functions for each type of Jacobian computation (corner, cell)
 */
template <MInt nDim>
void StructuredGrid<nDim>::computeJacobian() {
  TRACE();
  computeCornerJacobian();
  computeCellJacobian();
  IF_CONSTEXPR(nDim == 2) computeSurfaceJacobian();
}

/**
 * \brief Computes the corner Jacobian, i.e., the Jacobian of the volume composed by 4 cell centers (2D)
 */
template <>
void StructuredGrid<2>::computeCornerJacobian() {
  TRACE();
  // jacobian in the physical space is the inverse of the jacobian in computational space
  constexpr MInt nDim = 2;

  maia::parallelFor<true, nDim>(
      cellBegin(m_noGhostLayers - 1), cellEnd(m_noGhostLayers), [=](const MInt& i, const MInt& j) {
        const MInt cellId = cellIndex(i, j);
        const MFloat cornerJac =
            m_cells->cornerMetrics[xsd * 2 + xsd][cellId] * m_cells->cornerMetrics[ysd * 2 + ysd][cellId]
            - m_cells->cornerMetrics[ysd * 2 + xsd][cellId] * m_cells->cornerMetrics[xsd * 2 + ysd][cellId];

        // since metric terms are with omitted jacobian
        // there is factor of J^3; multiplied with J^-1 (invJac) we get J^2
        // --> take square root to get J
        m_cells->cornerJac[cellId] = cornerJac;
      });
}

/**
 * \brief Computes the corner Jacobian, i.e., the Jacobian of the volume composed by 8 cell centers (3D)
 */
template <>
void StructuredGrid<3>::computeCornerJacobian() {
  TRACE();
  // jacobian in the physical space is the inverse of the jacobian in computational space
  constexpr MInt nDim = 3;
  maia::parallelFor<true, nDim>(
      cellBegin(m_noGhostLayers - 1), cellEnd(m_noGhostLayers), [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt cellId = cellIndex(i, j, k);

        const MFloat invJac =
            m_cells->cornerMetrics[xsd * nDim + xsd][cellId]
                * (m_cells->cornerMetrics[ysd * nDim + ysd][cellId] * m_cells->cornerMetrics[zsd * nDim + zsd][cellId]
                   - m_cells->cornerMetrics[ysd * nDim + zsd][cellId]
                         * m_cells->cornerMetrics[zsd * nDim + ysd][cellId])
            - m_cells->cornerMetrics[ysd * nDim + xsd][cellId]
                  * (m_cells->cornerMetrics[xsd * nDim + ysd][cellId] * m_cells->cornerMetrics[zsd * nDim + zsd][cellId]
                     - m_cells->cornerMetrics[xsd * nDim + zsd][cellId]
                           * m_cells->cornerMetrics[zsd * nDim + ysd][cellId])
            + m_cells->cornerMetrics[zsd * nDim + xsd][cellId]
                  * (m_cells->cornerMetrics[xsd * nDim + ysd][cellId] * m_cells->cornerMetrics[ysd * nDim + zsd][cellId]
                     - m_cells->cornerMetrics[xsd * nDim + zsd][cellId]
                           * m_cells->cornerMetrics[ysd * nDim + ysd][cellId]);

        // since metric terms are with omitted jacobian
        // there is factor of J^3; multiplied with J^-1 (invJac) we get J^2
        // --> take square root to get J
        this->m_cells->cornerJac[cellId] = sqrt(invJac);
      });
}


/**
 * \brief Computes the modified corner Jacobian (3D) - Not required anymore
 */
template <>
void StructuredGrid<3>::computeModCornerJacobian() {
  TRACE();
  constexpr MInt nDim = 3;

  MFloat** __restrict coords = m_coordinates;

  constexpr MFloat fsqrttwo = F1 / SQRT2;

  constexpr MInt noSubJ = 8;
  MFloatScratchSpace subJ_(m_noCells, noSubJ, AT_, "subJ");
  MFloatScratchSpace subJtmp_(m_noCells, noSubJ, AT_, "subJtmp");
  // Aliasing the scratch space to a raw point is needed for PSTL on NVHPC
  MFloat* subJ = subJ_.data();
  MFloat* subJtmp = subJ_.data();

  maia::parallelFor<true, nDim>(
      cellBegin(m_noGhostLayers), cellEnd(m_noGhostLayers), [=](const MInt& i, const MInt& j, const MInt& k) {
        // tmp storage for vectors
        MFloat CP1[3];
        MFloat CP2[3];
        MFloat CP3[3];
        MFloat CP4[3];
        MFloat CP5[3];
        MFloat CP6[3];
        MFloat CP7[3];
        MFloat CP8[3];
        MFloat CP9[3];
        MFloat CP10[3];
        MFloat CP11[3];
        MFloat CP12[3];
        // vectors for metrics
        MFloat tmpX1[3];
        MFloat tmpX2[3];
        MFloat tmpX3[3];
        MFloat tmpX4[3];
        MFloat tmpX5[3];
        MFloat tmpX6[3];
        // tmp metric storage dxi
        MFloat DX1[3];
        // tmp metric storage deta
        MFloat DX2[3];
        // tmp metric storage dzeta
        MFloat DX3[3];

        // auxilliary variables for surface values
        MFloat S1[3];
        MFloat S2[3];
        MFloat S3[3];
        MFloat S1P[3];
        MFloat S2P[3];
        MFloat S3P[3];

        const MInt cellId = cellIndex(i, j, k);
        const MInt centCellId = cellIndex(i + 1, j + 1, k + 1);
        const MInt tmpId = getPointIdFromCell(i, j, k);

        const MInt ijk = getPointIdFromPoint(tmpId, 1, 1, 1);
        const MInt ijpk = getPointIdFromPoint(ijk, 0, 1, 0);
        const MInt ijkp = getPointIdFromPoint(ijk, 0, 0, 1);
        const MInt ijpkp = getPointIdFromPoint(ijk, 0, 1, 1);
        const MInt ipjk = getPointIdFromPoint(ijk, 1, 0, 0);
        const MInt ipjpk = getPointIdFromPoint(ijk, 1, 1, 0);
        const MInt ipjkp = getPointIdFromPoint(ijk, 1, 0, 1);
        const MInt ipjpkp = getPointIdFromPoint(ijk, 1, 1, 1);

        for(MInt isd = xsd; isd < nDim; isd++) {
          // averaging of the grid points of surface 1 (j+1/2,k+1/2) around corner point
          S1[isd] = F1B4 * (coords[isd][ijpkp] + coords[isd][ijk] + coords[isd][ijkp] + coords[isd][ijpk]);
          // averaging of the grid points of surface 2 (i+1/2,k+1/2) around corner point
          S2[isd] = F1B4 * (coords[isd][ipjk] + coords[isd][ijk] + coords[isd][ipjkp] + coords[isd][ijkp]);
          // averaging of the grid points of surface 3 (i+1/2,j+1/2) around corner point
          S3[isd] = F1B4 * (coords[isd][ipjk] + coords[isd][ijk] + coords[isd][ipjpk] + coords[isd][ijpk]);
          // averaging of the grid points of surface 1p (j+1/2,k+1/2) around corner point
          S1P[isd] = F1B4 * (coords[isd][ipjpkp] + coords[isd][ipjk] + coords[isd][ipjkp] + coords[isd][ipjpk]);
          // averaging of the grid points of surface 2p (i+1/2,k+1/2)
          S2P[isd] = F1B4 * (coords[isd][ipjpk] + coords[isd][ijpk] + coords[isd][ipjpkp] + coords[isd][ijpkp]);
          // averaging of the grid oints of surface 3p (i+1/2,j+1/2)
          S3P[isd] = F1B4 * (coords[isd][ipjkp] + coords[isd][ijkp] + coords[isd][ipjpkp] + coords[isd][ijpkp]);
        }


        ///////////////////////////////////////////
        ////////// subjacobian 1 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          // averaging corner points
          CP1[isd] = F1B2 * (coords[isd][ipjk] + coords[isd][ijk]);
          CP2[isd] = F1B2 * (coords[isd][ijpk] + coords[isd][ijk]);
          CP3[isd] = F1B2 * (coords[isd][ijkp] + coords[isd][ijk]);

          // setting up vectors for new metric terms
          tmpX1[isd] = (CP2[isd] - CP3[isd]) * fsqrttwo;
          tmpX2[isd] = (S1[isd] - coords[isd][ijk]) * fsqrttwo;

          tmpX3[isd] = (CP3[isd] - CP1[isd]) * fsqrttwo;
          tmpX4[isd] = (S2[isd] - coords[isd][ijk]) * fsqrttwo;

          tmpX5[isd] = (S3[isd] - coords[isd][ijk]) * fsqrttwo;
          tmpX6[isd] = (CP2[isd] - CP1[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 0] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 0] +=
              (m_cells->coordinates[isd][centCellId] - coords[isd][ijk]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 2 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP4[isd] = F1B2 * (coords[isd][ipjk] + coords[isd][ipjpk]);
          CP5[isd] = F1B2 * (coords[isd][ipjk] + coords[isd][ipjkp]);

          tmpX1[isd] = (S3[isd] - S2[isd]) * fsqrttwo;
          tmpX2[isd] = (m_cells->coordinates[isd][centCellId] - CP1[isd]) * fsqrttwo;

          tmpX3[isd] = (S2[isd] - coords[isd][ipjk]) * fsqrttwo;
          tmpX4[isd] = (CP5[isd] - CP1[isd]) * fsqrttwo;

          tmpX5[isd] = (CP4[isd] - CP1[isd]) * fsqrttwo;
          tmpX6[isd] = (S3[isd] - coords[isd][ipjk]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 1] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 1] += (S1P[isd] - CP1[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 3 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP6[isd] = F1B2 * (coords[isd][ipjpk] + coords[isd][ijpk]);
          CP7[isd] = F1B2 * (coords[isd][ijpkp] + coords[isd][ijpk]);

          tmpX1[isd] = (coords[isd][ijpk] - S1[isd]) * fsqrttwo;
          tmpX2[isd] = (CP7[isd] - CP2[isd]) * fsqrttwo;

          tmpX3[isd] = (S1[isd] - S3[isd]) * fsqrttwo;
          tmpX4[isd] = (m_cells->coordinates[isd][centCellId] - CP2[isd]) * fsqrttwo;

          tmpX5[isd] = (CP6[isd] - CP2[isd]) * fsqrttwo;
          tmpX6[isd] = (coords[isd][ijpk] - S3[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 2] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 2] += (S2P[isd] - CP2[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 4 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP8[isd] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ipjpk]);

          tmpX1[isd] = (CP6[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;
          tmpX2[isd] = (S2P[isd] - S3[isd]) * fsqrttwo;

          tmpX3[isd] = (m_cells->coordinates[isd][centCellId] - CP4[isd]) * fsqrttwo;
          tmpX4[isd] = (S1P[isd] - S3[isd]) * fsqrttwo;

          tmpX5[isd] = (coords[isd][ipjpk] - S3[isd]) * fsqrttwo;
          tmpX6[isd] = (CP6[isd] - CP4[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 3] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 3] += (CP8[isd] - S3[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 5 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP9[isd] = F1B2 * (coords[isd][ipjkp] + coords[isd][ijkp]);
          CP10[isd] = F1B2 * (coords[isd][ijpkp] + coords[isd][ijkp]);

          tmpX1[isd] = (S1[isd] - coords[isd][ijkp]) * fsqrttwo;
          tmpX2[isd] = (CP10[isd] - CP3[isd]) * fsqrttwo;

          tmpX3[isd] = (coords[isd][ijkp] - S2[isd]) * fsqrttwo;
          tmpX4[isd] = (CP9[isd] - CP3[isd]) * fsqrttwo;

          tmpX5[isd] = (m_cells->coordinates[isd][centCellId] - CP3[isd]) * fsqrttwo;
          tmpX6[isd] = (S1[isd] - S2[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 4] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 4] += (S3P[isd] - CP3[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 6 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP11[isd] = F1B2 * (coords[isd][ipjkp] + coords[isd][ipjpkp]);

          tmpX1[isd] = (m_cells->coordinates[isd][centCellId] - CP9[isd]) * fsqrttwo;
          tmpX2[isd] = (S3P[isd] - S2[isd]) * fsqrttwo;

          tmpX3[isd] = (CP9[isd] - CP5[isd]) * fsqrttwo;
          tmpX4[isd] = (coords[isd][ipjkp] - S2[isd]) * fsqrttwo;

          tmpX5[isd] = (S1P[isd] - S2[isd]) * fsqrttwo;
          tmpX6[isd] = (m_cells->coordinates[isd][centCellId] - CP5[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 5] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 5] += (CP11[isd] - S2[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 7 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          CP12[isd] = F1B2 * (coords[isd][ipjpkp] + coords[isd][ijpkp]);

          tmpX1[isd] = (CP7[isd] - CP10[isd]) * fsqrttwo;
          tmpX2[isd] = (coords[isd][ijpkp] - S1[isd]) * fsqrttwo;

          tmpX3[isd] = (CP10[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;
          tmpX4[isd] = (S3P[isd] - S1[isd]) * fsqrttwo;

          tmpX5[isd] = (S2P[isd] - S1[isd]) * fsqrttwo;
          tmpX6[isd] = (CP7[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 6] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 6] += (CP12[isd] - S1[isd]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }

        ///////////////////////////////////////////
        ////////// subjacobian 8 //////////////////
        ///////////////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          tmpX1[isd] = (S2P[isd] - S3P[isd]) * fsqrttwo;
          tmpX2[isd] = (CP12[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;

          tmpX3[isd] = (S3P[isd] - S1P[isd]) * fsqrttwo;
          tmpX4[isd] = (CP11[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;

          tmpX5[isd] = (CP8[isd] - m_cells->coordinates[isd][centCellId]) * fsqrttwo;
          tmpX6[isd] = (S2P[isd] - S1P[isd]) * fsqrttwo;
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        subJ[noSubJ * cellId + 7] = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          subJ[noSubJ * cellId + 7] +=
              (coords[isd][ipjpkp] - m_cells->coordinates[isd][centCellId]) * (DX1[isd] + DX2[isd] + DX3[isd]);
        }
      });

  //////////////////////////////////////////////
  ///// assemble subjacobians //////////////////
  //////////////////////////////////////////////

  // copy into dummy array
  maia::parallelFor<true>(0, m_noCells, [=](const MInt& cellId) {
    for(MInt j = 0; j < 8; j++) {
      subJtmp[noSubJ * cellId + j] = subJ[noSubJ * cellId + j];
    }
  });

  // shift subjacobians
  maia::parallelFor<true, nDim>(cellBegin(m_noGhostLayers), cellEnd(m_noGhostLayers),
                                [=](const MInt& i, const MInt& j, const MInt& k) {
                                  const MInt cellId = cellIndex(i, j, k);

                                  subJ[noSubJ * cellId + 0] = subJ[noSubJ * cellIndex(i - 1, j - 1, k - 1) + 7];
                                  subJ[noSubJ * cellId + 1] = subJ[noSubJ * cellIndex(i, j - 1, k - 1) + 6];
                                  subJ[noSubJ * cellId + 2] = subJ[noSubJ * cellIndex(i - 1, j, k - 1) + 5];
                                  subJ[noSubJ * cellId + 3] = subJ[noSubJ * cellIndex(i, j, k - 1) + 4];
                                });

  maia::parallelFor<true, nDim>(cellBegin(m_noGhostLayers - 1), cellEnd(m_noGhostLayers),
                                [=](const MInt& i, const MInt& j, const MInt& k) {
                                  const MInt cellId = cellIndex(i, j, k);
                                  subJ[noSubJ * cellId + 4] = subJtmp[noSubJ * cellIndex(i - 1, j - 1, k) + 3];
                                  subJ[noSubJ * cellId + 5] = subJtmp[noSubJ * cellIndex(i, j - 1, k) + 2];
                                  subJ[noSubJ * cellId + 6] = subJtmp[noSubJ * cellIndex(i - 1, j, k) + 1];
                                  subJ[noSubJ * cellId + 7] = subJtmp[noSubJ * cellId + 0];
                                });

  // finally jacobian at corner point!
  maia::parallelFor<true, nDim>(cellBegin(m_noGhostLayers - 1), cellEnd(m_noGhostLayers),
                                [=](const MInt& i, const MInt& j, const MInt& k) {
                                  const MInt cellId = cellIndex(i, j, k);

                                  m_cells->cornerJac[cellId] = F0;
                                  for(MInt jacId = 0; jacId < 8; jacId++) {
                                    m_cells->cornerJac[cellId] += subJ[noSubJ * cellId + jacId];
                                  }

                                  m_cells->cornerJac[cellId] = F1B3 * m_cells->cornerJac[cellId];
                                });
}

/**
 * \brief Computes the modified corner Jacobian (2D) - Not required anymore
 */
template <>
void StructuredGrid<2>::computeModCornerJacobian() {
  TRACE();
  constexpr MInt nDim = 2;
  MFloatScratchSpace subJ(this->m_noCells, 4, AT_, "subJ");
  MFloatScratchSpace subJtmp(this->m_noCells, 4, AT_, "subJtmp");

  MFloat** __restrict coords = m_coordinates;

  subJ.fill(1234.56);
  subJtmp.fill(5678.9);

  for(MInt j = m_noGhostLayers - 2; j < this->m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers - 2; i < this->m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt cellId = cellIndex(i, j);
      const MInt centCellId = cellIndex(i + 1, j + 1);
      const MInt tmpId = getPointIdFromCell(i, j);

      const MInt ij = getPointIdFromPoint(tmpId, 1, 1);
      const MInt ipj = getPointIdFromPoint(ij, 1, 0);
      const MInt ijp = getPointIdFromPoint(ij, 0, 1);
      const MInt ipjp = getPointIdFromPoint(ij, 1, 1);

      // auxilliary variables for surface values
      MFloat S1[2];
      MFloat S2[2];
      MFloat S1P[2];
      MFloat S2P[2];


      for(MInt isd = xsd; isd < nDim; isd++) {
        S1[isd] = F1B2 * (coords[isd][ijp] + coords[isd][ij]);
        S2[isd] = F1B2 * (coords[isd][ipj] + coords[isd][ij]);
        S1P[isd] = F1B2 * (coords[isd][ipjp] + coords[isd][ipj]);
        S2P[isd] = F1B2 * (coords[isd][ipjp] + coords[isd][ijp]);
      }


      // vectors for metrics
      MFloat tmpX1[2];
      MFloat tmpX2[2];

      ///////////////////////////////////////////
      ////////// subjacobian 1 //////////////////
      ///////////////////////////////////////////

      for(MInt isd = xsd; isd < nDim; isd++) {
        // averaging corner points is in 2D not necessary, cause of calculated surfaces mid points above

        // setting up vectors for new metric terms
        tmpX1[isd] = (m_cells->coordinates[isd][centCellId] - coords[isd][ij]) / sqrt(2);
        tmpX2[isd] = (S1[isd] - S2[isd]) / sqrt(2);
      }

      subJ(cellId, 0) = F0;
      subJ(cellId, 0) = crossProduct(tmpX1, tmpX2);

      ///////////////////////////////////////////
      ////////// subjacobian 2 //////////////////
      ///////////////////////////////////////////

      for(MInt isd = xsd; isd < nDim; isd++) {
        tmpX1[isd] = (S1P[isd] - S2[isd]) / sqrt(2);
        tmpX2[isd] = (m_cells->coordinates[isd][centCellId] - m_coordinates[isd][ipj]) / sqrt(2);
      }

      subJ(cellId, 1) = F0;
      subJ(cellId, 1) = crossProduct(tmpX1, tmpX2);

      ///////////////////////////////////////////
      ////////// subjacobian 3 //////////////////
      ///////////////////////////////////////////

      for(MInt isd = xsd; isd < nDim; isd++) {
        tmpX1[isd] = (S2P[isd] - S1[isd]) / sqrt(2);
        tmpX2[isd] = (m_coordinates[isd][ijp] - m_cells->coordinates[isd][centCellId]) / sqrt(2);
      }

      subJ(cellId, 2) = F0;
      subJ(cellId, 2) = crossProduct(tmpX1, tmpX2);

      ///////////////////////////////////////////
      ////////// subjacobian 4 //////////////////
      ///////////////////////////////////////////

      for(MInt isd = xsd; isd < nDim; isd++) {
        tmpX1[isd] = (m_coordinates[isd][ipjp] - m_cells->coordinates[isd][centCellId]) / sqrt(2);
        tmpX2[isd] = (S2P[isd] - S1P[isd]) / sqrt(2);
      }

      subJ(cellId, 3) = F0;
      subJ(cellId, 3) = crossProduct(tmpX1, tmpX2);
    }
  }


  //////////////////////////////////////////////
  ///// assemble subjacobians //////////////////
  //////////////////////////////////////////////

  // copy into dummy array
  for(MInt i = 0; i < m_noCells; i++) {
    for(MInt j = 0; j < 4; j++) {
      subJtmp(i, j) = subJ(i, j);
    }
  }

  // shift subjacobians

  for(MInt j = m_noGhostLayers - 1; j < this->m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < this->m_nCells[1] - m_noGhostLayers; i++) {
      MInt cellId = cellIndex(i, j);

      subJ(cellId, 0) = subJtmp(cellIndex(i - 1, j - 1), 3);
      subJ(cellId, 1) = subJtmp(cellIndex(i, j - 1), 2);
      subJ(cellId, 2) = subJtmp(cellIndex(i - 1, j), 1);
      subJ(cellId, 3) = subJtmp(cellIndex(i, j), 0);
    }
  }


  // finally jacobian at corner point!
  for(MInt j = m_noGhostLayers - 1; j < this->m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < this->m_nCells[1] - m_noGhostLayers; i++) {
      MInt cellId = cellIndex(i, j);

      m_cells->cornerJac[cellId] = F0;
      for(MInt jacId = 0; jacId < 4; jacId++) {
        m_cells->cornerJac[cellId] += subJ(cellId, jacId);
      }


      m_cells->cornerJac[cellId] = m_cells->cornerJac[cellId];
    }
  }
}

/**
 * \brief Computes the Jacobians of all cells (3D)
 */
template <>
void StructuredGrid<3>::computeCellJacobian() {
  TRACE();
  constexpr MInt nDim = 3;
  MFloat** __restrict coords = m_coordinates;

  maia::parallelFor<true, nDim>(
      cellBegin(m_noGhostLayers - 1), cellEnd(1), [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt cellId = cellIndex(i, j, k);

        const MInt ijk = getPointIdFromCell(i, j, k);
        const MInt ijpk = getPointIdFromPoint(ijk, 0, 1, 0);
        const MInt ijkp = getPointIdFromPoint(ijk, 0, 0, 1);
        const MInt ijpkp = getPointIdFromPoint(ijk, 0, 1, 1);
        const MInt ipjk = getPointIdFromPoint(ijk, 1, 0, 0);
        const MInt ipjpk = getPointIdFromPoint(ijk, 1, 1, 0);
        const MInt ipjkp = getPointIdFromPoint(ijk, 1, 0, 1);
        const MInt ipjpkp = getPointIdFromPoint(ijk, 1, 1, 1);

        MFloat DX1[3] = {F0, F0, F0};
        MFloat DX2[3] = {F0, F0, F0};
        MFloat DX3[3] = {F0, F0, F0};

        MFloat tmpX1[3] = {F0, F0, F0};
        MFloat tmpX2[3] = {F0, F0, F0};
        MFloat tmpX3[3] = {F0, F0, F0};
        MFloat tmpX4[3] = {F0, F0, F0};
        MFloat tmpX5[3] = {F0, F0, F0};
        MFloat tmpX6[3] = {F0, F0, F0};

        for(MInt isd = xsd; isd < nDim; isd++) {
          // setting up vectors for new metric terms
          tmpX1[isd] = (coords[isd][ijpk] - coords[isd][ijkp]);
          tmpX2[isd] = (coords[isd][ijpkp] - coords[isd][ijk]);

          tmpX3[isd] = (coords[isd][ijkp] - coords[isd][ipjk]);
          tmpX4[isd] = (coords[isd][ipjkp] - coords[isd][ijk]);

          tmpX5[isd] = (coords[isd][ipjpk] - coords[isd][ijk]);
          tmpX6[isd] = (coords[isd][ijpk] - coords[isd][ipjk]);
        }

        this->crossProduct(DX1, tmpX1, tmpX2);
        this->crossProduct(DX2, tmpX3, tmpX4);
        this->crossProduct(DX3, tmpX5, tmpX6);

        MFloat jac = F0;
        for(MInt isd = xsd; isd < nDim; isd++) {
          jac += (coords[isd][ipjpkp] - coords[isd][ijk]) * F1B2 * (DX1[isd] + DX2[isd] + DX3[isd]);
        }


        m_cells->cellJac[cellId] = F1B3 * fabs(jac);
      });
}

/**
 * \brief Computes the Jacobians of all cells (2D)
 */
template <>
void StructuredGrid<2>::computeCellJacobian() {
  TRACE();
  constexpr MInt nDim = 2;
  const MFloat* const* const RESTRICT coords = m_coordinates;

  maia::parallelFor<true, nDim>(cellBegin(m_noGhostLayers - 1), cellEnd(1), [=](const MInt& i, const MInt& j) {
    const MInt cellId = cellIndex(i, j);
    const MInt IJ = getPointIdFromCell(i, j);
    const MInt IPJ = getPointIdFromPoint(IJ, 1, 0);
    const MInt IJP = getPointIdFromPoint(IJ, 0, 1);
    const MInt IPJP = getPointIdFromPoint(IJ, 1, 1);

    MFloat diag1[2] = {F0, F0};
    MFloat diag2[2] = {F0, F0};

    for(MInt isd = xsd; isd < nDim; isd++) {
      diag1[isd] = (coords[isd][IPJP] - coords[isd][IJ]);
      diag2[isd] = (coords[isd][IJP] - coords[isd][IPJ]);
    }

    const MFloat area = F1B2 * this->crossProduct(diag1, diag2);

    m_cells->cellJac[cellId] = area;
  });
}


/**
 * \brief Computes the surface Jacobian (2D)
 * \return
 */
template <>
void StructuredGrid<2>::computeSurfaceJacobian() {
  TRACE();

  // Jacobian: Dxi/Dx * Deta/Dy - Deta/Dx * Dxi/Dy
  constexpr MInt nDim = 2;
  maia::parallelFor<true, nDim>(cellBegin(0), cellEnd(1), [=](const MInt& i, const MInt& j) {
    const MInt IJ = cellIndex(i, j);
    const MInt IPJ = cellIndex(i + 1, j);
    const MInt IJP = cellIndex(i, j + 1);
    const MInt ipjp = getPointIdFromPoint(getPointIdFromCell(i, j), 1, 1);
    const MInt ipj = getPointIdFromPoint(getPointIdFromCell(i, j), 1, 0);
    const MInt ijp = getPointIdFromPoint(getPointIdFromCell(i, j), 0, 1);

    MFloat DcoordDxi[nDim];
    MFloat DcoordDeta[nDim];

    for(MInt dim = 0; dim < nDim; ++dim) {
      // compute d(x,y,z)/dxi
      DcoordDxi[dim] = m_cells->coordinates[dim][IPJ] - m_cells->coordinates[dim][IJ];

      // compute d(x,y,z)/deta
      DcoordDeta[dim] = m_coordinates[dim][ipjp] - m_coordinates[dim][ipj];
    }

    this->m_cells->surfJac[IJ] = DcoordDeta[1] * DcoordDxi[0] - DcoordDxi[1] * DcoordDeta[0];

    for(MInt dim = 0; dim < nDim; ++dim) {
      // compute d(x,y,z)/dxi
      DcoordDxi[dim] = m_coordinates[dim][ipjp] - m_coordinates[dim][ijp];

      // compute d(x,y,z)/deta
      DcoordDeta[dim] = m_cells->coordinates[dim][IJP] - m_cells->coordinates[dim][IJ];
    }

    this->m_cells->surfJac[m_noCells + IJ] = DcoordDeta[1] * DcoordDxi[0] - DcoordDxi[1] * DcoordDeta[0];
  });
}


/**
 * \brief Computes the surface Jacobian (3D)
 * \return
 */
template <>
void StructuredGrid<3>::computeSurfaceJacobian() {
  TRACE();
  mTerm(1, AT_, "Not implemented for 3D");
}

/**
 * \brief Computes the metrics of the cell at the cell center (2D)
 */
template <>
void StructuredGrid<2>::computeCellMetrics() {
  TRACE();
  constexpr MInt nDim = 2;
  maia::parallelFor<true, nDim>(cellBegin(m_noGhostLayers - 1), cellEnd(1), [=](const MInt& i, const MInt& j) {
    const MInt cellId = cellIndex(i, j);
    // auxilliary variables
    MFloat DcoordDxi[2];
    MFloat DcoordDeta[2];

    for(MInt isd = xsd; isd < nDim; isd++) {
      DcoordDxi[isd] =
          (m_cells->coordinates[isd][cellIndex(i + 1, j)] - m_cells->coordinates[isd][cellIndex(i - 1, j)]) * F1B2;

      DcoordDeta[isd] =
          (m_cells->coordinates[isd][cellIndex(i, j + 1)] - m_cells->coordinates[isd][cellIndex(i, j - 1)]) * F1B2;
    }

    m_cells->cellMetrics[0][cellId] = DcoordDeta[1];  // DxiDx
    m_cells->cellMetrics[1][cellId] = -DcoordDeta[0]; // DxiDy
    m_cells->cellMetrics[2][cellId] = -DcoordDxi[1];  // DetaDx
    m_cells->cellMetrics[3][cellId] = DcoordDxi[0];   // DetaDy
  });
}

/**
 * \brief Computes the metrics of the cell at the cell center (3D)
 */
template <>
void StructuredGrid<3>::computeCellMetrics() {
  TRACE();
  constexpr MInt nDim = 3;
  const MFloat* const* const RESTRICT coords = m_cells->coordinates;

  maia::parallelFor<true, nDim>(
      cellBegin(m_noGhostLayers - 1), cellEnd(1), [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt cellId = cellIndex(i, j, k);
        // auxilliary variables
        MFloat DcoordDxi[3];
        MFloat DcoordDeta[3];
        MFloat DcoordDzeta[3];

        MFloat metricTmp[3];

        for(MInt isd = xsd; isd < nDim; isd++) {
          DcoordDxi[isd] = (coords[isd][cellIndex(i + 1, j, k)] - coords[isd][cellIndex(i - 1, j, k)]) * F1B2;

          DcoordDeta[isd] = (coords[isd][cellIndex(i, j + 1, k)] - coords[isd][cellIndex(i, j - 1, k)]) * F1B2;

          DcoordDzeta[isd] = (coords[isd][cellIndex(i, j, k + 1)] - coords[isd][cellIndex(i, j, k - 1)]) * F1B2;
        }

        // compute metric terms and store them

        // dxi
        this->crossProduct(metricTmp, DcoordDeta, DcoordDzeta);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cellMetrics[xsd * nDim + isd][cellId] = metricTmp[isd];
        }
        // deta
        this->crossProduct(metricTmp, DcoordDzeta, DcoordDxi);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cellMetrics[ysd * nDim + isd][cellId] = metricTmp[isd];
        }
        // dzeta
        this->crossProduct(metricTmp, DcoordDxi, DcoordDeta);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cellMetrics[zsd * nDim + isd][cellId] = metricTmp[isd];
        }
      });
}

/**
 * \brief Computes the surface metrics for the cell surfaces at the surface centroids (3D)
 */
template <>
void StructuredGrid<3>::computeSurfaceMetrics() {
  TRACE();
  constexpr MInt nDim = 3;
  const MFloat* const* const RESTRICT coords = m_coordinates;

  maia::parallelFor<true, nDim>(cellBegin(0), cellEnd(0), [=](const MInt& i, const MInt& j, const MInt& k) {
    // determine global cell ID
    const MInt cellId = cellIndex(i, j, k);
    // determine global point ID for local cell IDs
    const MInt ijk = getPointIdFromCell(i, j, k);
    const MInt ipjk = getPointIdFromPoint(ijk, 1, 0, 0);
    const MInt ipjpk = getPointIdFromPoint(ijk, 1, 1, 0);
    const MInt ipjkp = getPointIdFromPoint(ijk, 1, 0, 1);
    const MInt ipjpkp = getPointIdFromPoint(ijk, 1, 1, 1);
    const MInt ijpk = getPointIdFromPoint(ijk, 0, 1, 0);
    const MInt ijpkp = getPointIdFromPoint(ijk, 0, 1, 1);
    const MInt ijkp = getPointIdFromPoint(ijk, 0, 0, 1);

    // auxilliary variables
    MFloat DcoordDxi[3] = {F0, F0, F0};
    MFloat DcoordDeta[3] = {F0, F0, F0};
    MFloat DcoordDzeta[3] = {F0, F0, F0};

    MFloat metricTmp[3] = {F0, F0, F0};

    //////////////////////////////////////////////////////
    ////////////////////// FACE I ////////////////////////
    //////////////////////////////////////////////////////

    for(MInt isd = xsd; isd < nDim; isd++) {
      DcoordDeta[isd] = ((coords[isd][ipjpkp] + coords[isd][ipjpk]) - (coords[isd][ipjkp] + coords[isd][ipjk])) * F1B2;
      DcoordDzeta[isd] = ((coords[isd][ipjpkp] + coords[isd][ipjkp]) - (coords[isd][ipjpk] + coords[isd][ipjk])) * F1B2;
    }

    // compute Dxi and store
    crossProduct(metricTmp, DcoordDeta, DcoordDzeta);

    for(MInt isd = xsd; isd < nDim; isd++) {
      m_cells->surfaceMetrics[xsd * nDim + isd][cellId] = metricTmp[isd];
    }

    ///////////////////////////////////////////////////////
    ////////////////////// FACE J /////////////////////////
    ///////////////////////////////////////////////////////

    for(MInt isd = xsd; isd < nDim; isd++) {
      DcoordDxi[isd] = ((coords[isd][ipjpkp] + coords[isd][ipjpk]) - (coords[isd][ijpkp] + coords[isd][ijpk])) * F1B2;

      DcoordDzeta[isd] = ((coords[isd][ijpkp] + coords[isd][ipjpkp]) - (coords[isd][ipjpk] + coords[isd][ijpk])) * F1B2;
    }

    // compute Deta and store
    crossProduct(metricTmp, DcoordDzeta, DcoordDxi);

    for(MInt isd = xsd; isd < nDim; isd++) {
      m_cells->surfaceMetrics[ysd * nDim + isd][cellId] = metricTmp[isd];
    }

    ///////////////////////////////////////////////////////
    ////////////////////// FACE K /////////////////////////
    ///////////////////////////////////////////////////////

    for(MInt isd = xsd; isd < nDim; isd++) {
      DcoordDxi[isd] = ((coords[isd][ipjpkp] + coords[isd][ipjkp]) - (coords[isd][ijpkp] + coords[isd][ijkp])) * F1B2;

      DcoordDeta[isd] = ((coords[isd][ipjpkp] + coords[isd][ijpkp]) - (coords[isd][ipjkp] + coords[isd][ijkp])) * F1B2;
    }

    // compute Dzeta and store
    crossProduct(metricTmp, DcoordDxi, DcoordDeta);

    for(MInt isd = xsd; isd < nDim; isd++) {
      m_cells->surfaceMetrics[zsd * nDim + isd][cellId] = metricTmp[isd];
    }
  });
}

/**
 * \brief Computes the surface metrics for the cell surfaces at the surface centroids (2D)
 */
template <>
void StructuredGrid<2>::computeSurfaceMetrics() {
  TRACE();
  constexpr MInt nDim = 2;
  maia::parallelFor<true, nDim>(cellBegin(0), cellEnd(1), [=](const MInt& i, const MInt& j) {
    // determine global cellID
    const MInt cellId = this->cellIndex(i, j);
    // determine global point ID for local cell IDs
    const MInt IJ = getPointIdFromCell(i, j);
    const MInt IPJ = getPointIdFromPoint(IJ, 1, 0);
    const MInt IPJP = getPointIdFromPoint(IJ, 1, 1);
    const MInt IJP = getPointIdFromPoint(IJ, 0, 1);

    // auxiliary variables
    MFloat DcoordDxi[2];

    // Face I //
    for(MInt isd = xsd; isd < nDim; ++isd) {
      DcoordDxi[isd] = m_coordinates[isd][IPJP] - m_coordinates[isd][IPJ];
    }

    // compute Dxi
    m_cells->surfaceMetrics[0][cellId] = DcoordDxi[1];
    m_cells->surfaceMetrics[1][cellId] = -DcoordDxi[0];
    // store Dxi

    // Face I //
    for(MInt isd = xsd; isd < nDim; ++isd) {
      DcoordDxi[isd] = m_coordinates[isd][IPJP] - m_coordinates[isd][IJP];
    }

    m_cells->surfaceMetrics[2][cellId] = -DcoordDxi[1];
    m_cells->surfaceMetrics[3][cellId] = DcoordDxi[0];
  });
}

/**
 * \brief Computes the surface metrics for the cell surfaces at the surface centroids (2D)
          Special version for cells at grid singularities
 */
template <MInt nDim>
void StructuredGrid<nDim>::computeSurfaceMetricsSingularity() {
  TRACE();

  IF_CONSTEXPR(nDim == 3) TERMM(1, "Not implemented in 3D!");

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC == -6000) {
      // Sanity check
      for(MInt j = 0; j < nDim; j++) {
        ASSERT(m_singularity[i].end[j] - m_singularity[i].start[j] == 1, "");
      }

      for(MInt jj = m_singularity[i].start[1]; jj < m_singularity[i].end[1]; ++jj) {
        for(MInt ii = m_singularity[i].start[0]; ii < m_singularity[i].end[0]; ++ii) {
          // Cell Id of singularity cell
          const MInt IJ = cellIndex(ii, jj);

          const MInt sign_xi = 2 * m_singularity[i].Viscous[0] + 1;
          const MInt sign_eta = 2 * m_singularity[i].Viscous[1] + 1;

          const MInt IPMJ = getCellIdFromCell(IJ, sign_xi, 0);
          const MInt IJPM = getCellIdFromCell(IJ, 0, sign_eta);

          // auxiliary variables
          MFloat DcoordD[2];

          // Face I //
          for(MInt isd = xsd; isd < nDim; ++isd) {
            DcoordD[isd] = sign_xi * (m_cells->coordinates[isd][IPMJ] - m_cells->coordinates[isd][IJ]);
          }

          // compute Deta
          m_cells->surfaceMetricsSingularity[2][i] = -DcoordD[1];
          m_cells->surfaceMetricsSingularity[3][i] = DcoordD[0];
          // store Deta

          // Face I //
          for(MInt isd = xsd; isd < nDim; ++isd) {
            DcoordD[isd] = sign_eta * (m_cells->coordinates[isd][IJPM] - m_cells->coordinates[isd][IJ]);
          }

          m_cells->surfaceMetricsSingularity[0][i] = DcoordD[1];
          m_cells->surfaceMetricsSingularity[1][i] = -DcoordD[0];

          //          cout << "dom=" << domainId() << " x|y=" << setprecision(10) << m_cells->coordinates[0][IJ] << "|"
          //          << m_cells->coordinates[1][IJ]
          //            << " " << m_cells->surfaceMetricsSingularity[i][0] << "|" <<
          //            m_cells->surfaceMetricsSingularity[i][1]
          //            << " " << m_cells->surfaceMetricsSingularity[i][2] << "|" <<
          //            m_cells->surfaceMetricsSingularity[i][3]
          //            << " IPMJ=" << m_cells->coordinates[0][IPMJ] << "|" << m_cells->coordinates[1][IPMJ]
          //            << " IJPM=" << m_cells->coordinates[0][IJPM] << "|" << m_cells->coordinates[1][IJPM] << endl;
        }
      }
    }
  }
}

/**
 * \brief Computes the corner metrics for the volume composed the surrounding 8 cell centers (3D)
 */
template <>
void StructuredGrid<3>::computeCornerMetrics() {
  TRACE();
  constexpr MInt nDim = 3;

  const MFloat* const* const RESTRICT coords = m_coordinates;

  for(MInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < this->m_nCells[2] - m_noGhostLayers; i++) {
        // determine global cell ID
        MInt cellId = this->cellIndex(i, j, k); // i + ( k * m_nCells[1] + j ) * m_nCells[2];

        // auxilliary variables
        MFloat DcoordDxi[3] = {F0, F0, F0};
        MFloat DcoordDeta[3] = {F0, F0, F0};
        MFloat DcoordDzeta[3] = {F0, F0, F0};

        MFloat metricTmp[3];

        for(MInt isd = xsd; isd < nDim; isd++) {
          // compute d(x,y,z)/dxi
          DcoordDxi[isd] = F1B2
                           * (coords[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 1, 1, 1)]
                              - coords[isd][getPointIdFromPoint(getPointIdFromCell(i - 1, j, k), 1, 1, 1)]);

          // compute d(x,y,z)/deta
          DcoordDeta[isd] = F1B2
                            * (coords[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 1, 1)]
                               - coords[isd][getPointIdFromPoint(getPointIdFromCell(i, j - 1, k), 1, 1, 1)]);

          // compute d(x,y,z)/dzeta
          DcoordDzeta[isd] = F1B2
                             * (coords[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 1, 1)]
                                - coords[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k - 1), 1, 1, 1)]);
        }

        // compute metric terms and store them

        // dxi
        this->crossProduct(metricTmp, DcoordDeta, DcoordDzeta);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[xsd * nDim + isd][cellId] = metricTmp[isd];
        }
        // deta
        this->crossProduct(metricTmp, DcoordDzeta, DcoordDxi);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[ysd * nDim + isd][cellId] = metricTmp[isd];
        }
        // dzeta
        this->crossProduct(metricTmp, DcoordDxi, DcoordDeta);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[zsd * nDim + isd][cellId] = metricTmp[isd];
        }
      }
    }
  }
}

/**
 * \brief Computes the corner metrics for the volume composed the surrounding 4 cell centers (2D)
 */
template <>
void StructuredGrid<2>::computeCornerMetrics() {
  TRACE();
  constexpr MInt nDim = 2;
  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
      // determine global cell ID
      const MInt cellId = cellIndex(i, j);

      // auxilliary variables
      MFloat DcoordDxi[2];
      MFloat DcoordDeta[2];

      for(MInt isd = xsd; isd < nDim; ++isd) {
        // looks complicated, but what happens is that we always catch the point Id of ipjp
        // from the neighboring cell and build the centered difference

        // compute d(x,y,z)/dxi
        DcoordDxi[isd] = F1B2
                         * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j), 1, 1)]
                            - m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i - 1, j), 1, 1)]);

        // compute d(x,y,z)/deta
        DcoordDeta[isd] = F1B2
                          * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1), 1, 1)]
                             - m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j - 1), 1, 1)]);
      }
      m_cells->cornerMetrics[0][cellId] = DcoordDeta[1];  // DxiDx
      m_cells->cornerMetrics[1][cellId] = -DcoordDeta[0]; // DxiDy
      m_cells->cornerMetrics[2][cellId] = -DcoordDxi[1];  // DetaDx
      m_cells->cornerMetrics[3][cellId] = DcoordDxi[0];   // DetaDy
    }
  }
}

/**
 * \brief Computes the modified corner metrics, call is redirected to the standard corner metrics function (2D)
 */
template <>
void StructuredGrid<2>::computeModCornerMetrics() {
  TRACE();
  computeCornerMetrics();
}

/**
 * \brief More accurate version of the corner metric computation (3D)
 */
template <>
void StructuredGrid<3>::computeModCornerMetrics() {
  TRACE();
  m_log << "computing corner metrics ..." << endl;
  constexpr MInt nDim = 3;

  for(MInt k = m_noGhostLayers - 1; k < this->m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers - 1; j < this->m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < this->m_nCells[2] - m_noGhostLayers; i++) {
        // determine global cell ID
        MInt cellId = this->cellIndex(i, j, k); // i + ( k * m_nCells[1] + j ) * m_nCells[2];
        MFloat metricTmp[3];

        MFloat p1[3];
        MFloat p2[3];
        MFloat p3[3];
        MFloat p4[3];

        MFloat diag1[3];
        MFloat diag2[3];


        ////////////////////////////////
        ////////// DXI /////////////////
        ////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          // looks complicated, but what happens is that we always catch the point Id of ipjpkp
          // from the neighboring cell and build the centered difference

          p1[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k + 1), 1, 0, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k + 1), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k + 1), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k + 1), 1, 1, 1)]);


          p2[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 0, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 1, 1)]);


          p3[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 0, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 1, 1)]);

          p4[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 0, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 1, 1)]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct(metricTmp, diag1, diag2);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[xsd * nDim + isd][cellId] = F1B2 * metricTmp[isd];
        }


        ////////////////////////////////
        ////////// DETA ////////////////
        ////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          // looks complicated, but what happens is that we always catch the point Id of ipjpkp
          // from the neighboring cell and build the centered difference

          p1[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k + 1), 0, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k + 1), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k + 1), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k + 1), 1, 1, 1)]);


          p2[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 0, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 1, 1)]);

          p3[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 0, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 1, 1, 1)]);

          p4[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 0, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 1, 0)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k + 1), 1, 1, 1)]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct(metricTmp, diag1, diag2);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[ysd * nDim + isd][cellId] = F1B2 * metricTmp[isd];
        }

        ////////////////////////////////
        ////////// DZETA ///////////////
        ////////////////////////////////

        for(MInt isd = xsd; isd < nDim; isd++) {
          p1[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j + 1, k), 1, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j + 1, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j + 1, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j + 1, k), 0, 0, 1)]);


          p2[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j, k), 0, 0, 1)]);

          p3[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i, j + 1, k), 0, 0, 1)]);

          p4[isd] = F1B4
                    * (m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 1, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 0, 1, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 1, 0, 1)]
                       + m_coordinates[isd][getPointIdFromPoint(getPointIdFromCell(i + 1, j, k), 0, 0, 1)]);

          diag1[isd] = p1[isd] - p2[isd];
          diag2[isd] = p3[isd] - p4[isd];
        }

        this->crossProduct(metricTmp, diag1, diag2);
        for(MInt isd = xsd; isd < nDim; isd++) {
          m_cells->cornerMetrics[zsd * nDim + isd][cellId] = F1B2 * metricTmp[isd];
        }
      }
    }
  }
}

/**
 * \brief Computes the volume flux for three surfaces (positive i, positive j, positive k)
          for moving grids (3D)
 *
 * \author Marian Albers
 *
 * \param[in] timeStep Current time step
 * \param[in] RKalpha Array with Runge-Kutta alpha coefficients
 * \param[in] RKStep Current Runge-Kutta substeps
 */
template <>
void StructuredGrid<3>::computeDxt(MFloat timeStep, MFloat* RKalpha, MInt RKStep) {
  TRACE();
  const MFloat frk = F1 / (timeStep * RKalpha[RKStep]);
  constexpr MInt nDim = 3;

  const MInt IJK[3] = {m_nCells[2], m_nCells[1], m_nCells[0]};

  const MFloat* const* const RESTRICT coords = m_coordinates;
  const MFloat* const* const RESTRICT oldCoords = m_oldCoordinates;

  const MInt inc[12] = {getPointIdFromPoint(0, 1, 0, 0), getPointIdFromPoint(0, 1, 1, 0),
                        getPointIdFromPoint(0, 1, 0, 1), getPointIdFromPoint(0, 1, 1, 1),

                        getPointIdFromPoint(0, 0, 1, 0), getPointIdFromPoint(0, 0, 1, 1),
                        getPointIdFromPoint(0, 1, 1, 0), getPointIdFromPoint(0, 1, 1, 1),

                        getPointIdFromPoint(0, 0, 0, 1), getPointIdFromPoint(0, 1, 0, 1),
                        getPointIdFromPoint(0, 0, 1, 1), getPointIdFromPoint(0, 1, 1, 1)};


  for(MInt dim = 0; dim < nDim; ++dim) {
    for(MInt k = 0; k < IJK[2]; ++k) {
      for(MInt j = 0; j < IJK[1]; ++j) {
        for(MInt i = 0; i < IJK[0]; ++i) {
          // determine global cell ID
          const MInt cellId = cellIndex(i, j, k);

          // determine global point ID for local cell IDs
          const MInt ijk = getPointIdFromCell(i, j, k);

          const MInt p0 = ijk + inc[dim * 4 + 0];
          const MInt p1 = ijk + inc[dim * 4 + 1];
          const MInt p2 = ijk + inc[dim * 4 + 2];
          const MInt p3 = ijk + inc[dim * 4 + 3];

          // reset the volume flux
          m_cells->dxt[dim][cellId] = F0;

          MFloat diag1[3] = {F0, F0, F0};
          MFloat diag2[3] = {F0, F0, F0};
          MFloat oldNormal[3] = {F0, F0, F0};
          MFloat newNormal1[3] = {F0, F0, F0};
          MFloat newNormal2[3] = {F0, F0, F0};

          //(i+1/2,j-1/4,k-1/4)

          for(MInt isd = xsd; isd < nDim; isd++) {
            diag1[isd] = (oldCoords[isd][p3] - oldCoords[isd][p0]);
            diag2[isd] = (oldCoords[isd][p2] - oldCoords[isd][p1]);
          }

          this->crossProduct(oldNormal, diag1, diag2);

          for(MInt isd = xsd; isd < nDim; isd++) {
            diag1[isd] = (coords[isd][p2] - oldCoords[isd][p0]);
            diag2[isd] = (coords[isd][p0] - oldCoords[isd][p2]);
          }

          this->crossProduct(newNormal1, diag1, diag2);

          for(MInt isd = xsd; isd < nDim; isd++) {
            diag1[isd] = (coords[isd][p0] - oldCoords[isd][p1]);
            diag2[isd] = (coords[isd][p1] - oldCoords[isd][p0]);
          }

          this->crossProduct(newNormal2, diag1, diag2);

          const MFloat jac =
              F1B3
              * ((coords[0][p3] - oldCoords[0][p0]) * (oldNormal[0] + newNormal1[0] + newNormal2[0]) * F1B2
                 + (coords[1][p3] - oldCoords[1][p0]) * (oldNormal[1] + newNormal1[1] + newNormal2[1]) * F1B2
                 + (coords[2][p3] - oldCoords[2][p0]) * (oldNormal[2] + newNormal1[2] + newNormal2[2]) * F1B2);

          m_cells->dxt[dim][cellId] = jac * frk;
        }
      }
    }
  }
}

/**
 * \brief Computes the volume flux for three surfaces (positive i, positive j)
          for moving grids (2D)
 *
 * \author Marian Albers
 *
 * \param[in] timeStep Current time step
 * \param[in] RKalpha Array with Runge-Kutta alpha coefficients
 * \param[in] RKStep Current Runge-Kutta substeps
 */
template <>
void StructuredGrid<2>::computeDxt(MFloat timeStep, MFloat* RKalpha, MInt RKStep) {
  const MFloat frk = F1 / (timeStep * RKalpha[RKStep]);
  constexpr MInt nDim = 2;

  const MInt IJ[2] = {m_nCells[0], m_nCells[1]};

  const MFloat* const* const RESTRICT coords = m_coordinates;
  const MFloat* const* const RESTRICT oldCoords = m_oldCoordinates;

  const MInt inc[4] = {
      getPointIdFromPoint(0, 1, 1),
      getPointIdFromPoint(0, 1, 0),
      getPointIdFromPoint(0, 0, 1),
      getPointIdFromPoint(0, 1, 1),
  };

  for(MInt dim = 0; dim < nDim; ++dim) {
    for(MInt j = 0; j < IJ[0]; ++j) {
      for(MInt i = 0; i < IJ[1]; ++i) {
        // determine global cell ID
        const MInt cellId = cellIndex(i, j);

        // determine global point ID for local cell IDs
        const MInt ij = getPointIdFromCell(i, j);

        const MInt p0 = ij + inc[dim * 2 + 0];
        const MInt p1 = ij + inc[dim * 2 + 1];

        // reset the volume flux
        m_cells->dxt[dim][cellId] = F0;

        MFloat diag1[2] = {F0, F0};
        MFloat diag2[2] = {F0, F0};

        for(MInt isd = xsd; isd < nDim; isd++) {
          diag1[isd] = (oldCoords[isd][p0] - coords[isd][p1]);
          diag2[isd] = (oldCoords[isd][p1] - coords[isd][p0]);
        }

        const MFloat area = F1B2 * (diag1[0] * diag2[1] - diag1[1] * diag2[0]);

        m_cells->dxt[dim][cellId] = area * frk;
      }
    }
  }
}

/**
 * \brief Copies the current state of the grid coordinates to m_oldCoordinates (3D)
 */
template <>
void StructuredGrid<3>::saveGrid() {
  TRACE();
  constexpr MInt nDim = 3;
  for(MInt k = 0; k < m_nPoints[0]; ++k) {
    for(MInt j = 0; j < m_nPoints[1]; ++j) {
      for(MInt i = 0; i < m_nPoints[2]; ++i) {
        const MInt pointId = pointIndex(i, j, k);
        for(MInt isd = xsd; isd < nDim; ++isd) {
          m_oldCoordinates[isd][pointId] = m_coordinates[isd][pointId];
        }
      }
    }
  }
}

/**
 * \brief Copies the current state of the grid coordinates to m_oldCoordinates (2D)
 */
template <>
void StructuredGrid<2>::saveGrid() {
  TRACE();
  constexpr MInt nDim = 2;
  for(MInt j = 0; j < m_nPoints[0]; ++j) {
    for(MInt i = 0; i < m_nPoints[1]; ++i) {
      const MInt pointId = pointIndex(i, j);
      for(MInt isd = xsd; isd < nDim; ++isd) {
        m_oldCoordinates[isd][pointId] = m_coordinates[isd][pointId];
      }
    }
  }
}

/**
 * \brief Copies the current state of the cell Jacobians to m_cells->oldCellJacobian
 */
template <MInt nDim>
void StructuredGrid<nDim>::saveCellJacobian() // saves the old cell Jacobian
{
  TRACE();
  MFloat* const RESTRICT oldCellJac = ALIGNED_MF(m_cells->oldCellJac);
  const MFloat* const RESTRICT cellJac = ALIGNED_MF(m_cells->cellJac);
  for(MInt cellId = 0; cellId < m_noCells; cellId++) {
    oldCellJac[cellId] = cellJac[cellId];
  }
}

/**
 * \brief Simple extrapolation of the grid coordinates onto the ghost point coordinates at the grid boundaries (3D)
 */
template <>
void StructuredGrid<3>::extrapolateGhostPointCoordinates() {
  TRACE();
  // This function mirrors the grid points on the faces
  constexpr MInt nDim = 3;
  // i-direction
  MInt pointId, FixPointId, MirrorPointId;

  for(MInt k = m_noGhostLayers; k < (m_nPoints[0] - m_noGhostLayers); k++) {
    for(MInt j = m_noGhostLayers; j < (m_nPoints[1] - m_noGhostLayers); j++) {
      for(MInt i = 0; i < m_noGhostLayers; i++) {
        pointId = (m_noGhostLayers - 1 - i) + (j + k * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            (m_noGhostLayers - i) + (j + k * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId = (m_noGhostLayers + 1 - i) + (j + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
        // coordinates at the other end!!

        pointId = (m_nPoints[2] - i - 1) + (j + k * m_nPoints[1]) * m_nPoints[2];
        FixPointId = (m_nPoints[2] - m_noGhostLayers - 1) + (j + k * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = (m_nPoints[2] - 1 - (2 * m_noGhostLayers - i)) + (j + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  // j-direction

  for(MInt k = m_noGhostLayers; k < (m_nPoints[0] - m_noGhostLayers); k++) {
    for(MInt j = 0; j < m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < (m_nPoints[2] - m_noGhostLayers); i++) {
        pointId = i + ((m_noGhostLayers - j - 1) + k * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            i + ((m_noGhostLayers - j) + k * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId = i + ((m_noGhostLayers + 1 - j) + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
        // coordinates at the other end!!
        pointId = i + ((m_nPoints[1] - j - 1) + k * m_nPoints[1]) * m_nPoints[2];
        FixPointId = i + ((m_nPoints[1] - m_noGhostLayers - 1) + k * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = i + ((m_nPoints[1] - 1 - (2 * m_noGhostLayers - j)) + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  // k-direction
  for(MInt k = 0; k < m_noGhostLayers; k++) {
    for(MInt j = 0; j < (m_nPoints[1] - m_noGhostLayers); j++) {
      for(MInt i = m_noGhostLayers; i < (m_nPoints[2] - m_noGhostLayers); i++) {
        pointId = i + (j + (m_noGhostLayers - 1 - k) * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            i + (j + (m_noGhostLayers - k) * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId =
            i + (j + (m_noGhostLayers + 1 - k) * m_nPoints[1]) * m_nPoints[2]; // m_noGhostLayers+(m_noGhostLayers-i)
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }


        // coordinates at the other end!!
        pointId = i + (j + (m_nPoints[0] - k - 1) * m_nPoints[1]) * m_nPoints[2];
        FixPointId = i + (j + (m_nPoints[0] - m_noGhostLayers - 1) * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = i + (j + (m_nPoints[0] - 1 - (2 * m_noGhostLayers - k)) * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  // corner points missing yet!! They only need to be calculated if a Visualisation tool is used to
  // show the ghost points and the grid


  // in i-direction

  for(MInt k = 0; k < (m_nPoints[0]); k++) {
    for(MInt j = 0; j < (m_nPoints[1]); j++) {
      for(MInt i = 0; i < m_noGhostLayers; i++) {
        pointId = (m_noGhostLayers - 1 - i) + (j + k * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            (m_noGhostLayers - i) + (j + k * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId =
            (m_noGhostLayers + 1 - i) + (j + k * m_nPoints[1]) * m_nPoints[2]; // m_noGhostLayers+(m_noGhostLayers-i)
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
        // coordinates at the other end!!

        pointId = (m_nPoints[2] - i - 1) + (j + k * m_nPoints[1]) * m_nPoints[2];
        FixPointId = (m_nPoints[2] - m_noGhostLayers - 1) + (j + k * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = (m_nPoints[2] - 1 - (2 * m_noGhostLayers - i)) + (j + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  for(MInt k = 0; k < (m_nPoints[0]); k++) {
    for(MInt j = 0; j < m_noGhostLayers; j++) {
      for(MInt i = 0; i < (m_nPoints[2]); i++) {
        pointId = i + ((m_noGhostLayers - j - 1) + k * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            i + ((m_noGhostLayers - j) + k * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId = i + ((m_noGhostLayers + 1 - j) + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
        // coordinates at the other end!!
        pointId = i + ((m_nPoints[1] - j - 1) + k * m_nPoints[1]) * m_nPoints[2];
        FixPointId = i + ((m_nPoints[1] - m_noGhostLayers - 1) + k * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = i + ((m_nPoints[1] - 1 - (2 * m_noGhostLayers - j)) + k * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }

  for(MInt k = 0; k < m_noGhostLayers; k++) {
    for(MInt j = 0; j < (m_nPoints[1]); j++) {
      for(MInt i = 0; i < (m_nPoints[2]); i++) {
        pointId = i + (j + (m_noGhostLayers - 1 - k) * m_nPoints[1]) * m_nPoints[2]; // pointId in Array
        FixPointId =
            i + (j + (m_noGhostLayers - k) * m_nPoints[1]) * m_nPoints[2]; // point about which everything is mirrored
        MirrorPointId = i + (j + (m_noGhostLayers + 1 - k) * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }

        // coordinates at the other end!!
        pointId = i + (j + (m_nPoints[0] - k - 1) * m_nPoints[1]) * m_nPoints[2];
        FixPointId = i + (j + (m_nPoints[0] - m_noGhostLayers - 1) * m_nPoints[1]) * m_nPoints[2];
        MirrorPointId = i + (j + (m_nPoints[0] - 1 - (2 * m_noGhostLayers - k)) * m_nPoints[1]) * m_nPoints[2];
        for(MInt dim = 0; dim < nDim; dim++) {
          m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
        }
      }
    }
  }
}

/**
 * \brief Simple extrapolation of the grid coordinates onto the ghost point coordinates at the grid boundaries (2D)
 */
template <>
void StructuredGrid<2>::extrapolateGhostPointCoordinates() {
  TRACE();
  // This function mirrors the grid points on the faces
  constexpr MInt nDim = 2;

  // i-direction
  MInt pointId, FixPointId, MirrorPointId;

  for(MInt j = m_noGhostLayers; j < (m_nPoints[0] - m_noGhostLayers); ++j) {
    for(MInt i = 0; i < m_noGhostLayers; ++i) {
      pointId = (m_noGhostLayers - 1 - i) + (j * m_nPoints[1]); // pointId in Array
      FixPointId = (m_noGhostLayers - i) + (j * m_nPoints[1]);  // point about which everything is mirrored
      MirrorPointId = (m_noGhostLayers + 1 - i) + (j * m_nPoints[1]);

      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }

      // coordinates at the other end!!
      pointId = (m_nPoints[1] - i - 1) + (j * m_nPoints[1]);
      FixPointId = (m_nPoints[1] - m_noGhostLayers - 1) + (j * m_nPoints[1]);
      MirrorPointId = (m_nPoints[1] - 1 - (2 * m_noGhostLayers - i)) + (j * m_nPoints[1]);
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }
    }
  }

  // j-direction

  for(MInt j = 0; j < m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers; i < (m_nPoints[1] - m_noGhostLayers); ++i) {
      pointId = i + (m_noGhostLayers - j - 1) * m_nPoints[1]; // pointId in Array
      FixPointId = i + (m_noGhostLayers - j) * m_nPoints[1];  // point about which everything is mirrored
      MirrorPointId = i + (m_noGhostLayers + 1 - j) * m_nPoints[1];
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }

      // coordinates at the other end!!
      pointId = i + (m_nPoints[0] - j - 1) * m_nPoints[1];
      FixPointId = i + (m_nPoints[0] - m_noGhostLayers - 1) * m_nPoints[1];
      MirrorPointId = i + (m_nPoints[0] - 1 - (2 * m_noGhostLayers - j)) * m_nPoints[1];
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }
    }
  }


  // corner points missing yet!! They only need to be calculated if a Visualisation tool is used to
  // show the ghost points and the grid


  // in i-direction
  for(MInt j = 0; j < (m_nPoints[0]); ++j) {
    for(MInt i = 0; i < m_noGhostLayers; ++i) {
      pointId = (m_noGhostLayers - 1 - i) + (j * m_nPoints[1]);       // pointId in Array
      FixPointId = (m_noGhostLayers - i) + (j * m_nPoints[1]);        // point about which everything is mirrored
      MirrorPointId = (m_noGhostLayers + 1 - i) + (j * m_nPoints[1]); // m_noGhostLayers+(m_noGhostLayers-i)
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }
      // coordinates at the other end!!

      pointId = (m_nPoints[1] - i - 1) + (j * m_nPoints[1]);
      FixPointId = (m_nPoints[1] - m_noGhostLayers - 1) + (j * m_nPoints[1]);
      MirrorPointId = (m_nPoints[1] - 1 - (2 * m_noGhostLayers - i)) + (j * m_nPoints[1]);
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }
    }
  }

  // in j-direction

  for(MInt j = 0; j < m_noGhostLayers; ++j) {
    for(MInt i = 0; i < (m_nPoints[1]); ++i) {
      pointId = i + (m_noGhostLayers - j - 1) * m_nPoints[1]; // pointId in Array
      FixPointId = i + (m_noGhostLayers - j) * m_nPoints[1];  // point about which everything is mirrored
      MirrorPointId = i + (m_noGhostLayers + 1 - j) * m_nPoints[1];
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }

      // coordinates at the other end!!
      pointId = i + (m_nPoints[0] - j - 1) * m_nPoints[1];
      FixPointId = i + (m_nPoints[0] - m_noGhostLayers - 1) * m_nPoints[1];
      MirrorPointId = i + (m_nPoints[0] - 1 - (2 * m_noGhostLayers - j)) * m_nPoints[1];
      for(MInt dim = 0; dim < nDim; ++dim) {
        m_coordinates[dim][pointId] = (2 * m_coordinates[dim][FixPointId] - m_coordinates[dim][MirrorPointId]);
      }
    }
  }
}

/**
 * \brief Destructor of the structured grid class
 * \return
 */
template <MInt nDim>
StructuredGrid<nDim>::~StructuredGrid() {}

// Explicit instantations for 2D and 3D
template class StructuredGrid<2>;
template class StructuredGrid<3>;
