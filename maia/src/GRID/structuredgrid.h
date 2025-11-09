// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDGRID_H
#define STRUCTUREDGRID_H

#include <vector>
#include "COMM/mpioverride.h"
#include "FV/fvstructuredcell.h"
#include "FV/fvstructuredcomm.h"
#include "FV/fvstructuredtimers.h"
#include "GRID/structuredpartition.h"

struct SingularInformation {
  MInt start[3];
  MInt end[3];
  MInt Nstar;
  MInt displacement[5][3];
  MInt count;
  MInt totalPoints;
  MInt totalCells;
  MFloat** coordinates = nullptr;
  MFloat** variables = nullptr;
  MFloat** ReconstructionConstants = nullptr;
  MInt Viscous[3];
  MInt BC;
  MInt BCsingular[6]{};
  MInt SingularBlockId[4];
};

/** \brief Structured grid class
 *
 * This class contains the functions for loading the structured grid from the grid file
 * and offers functions for parallel handling of the grid. Furthermore it contains
 * functions to compute grid-dependent quantities (metrics, Jacobians).
 *
 */
template <MInt nDim>
class StructuredGrid {
 public:
  StructuredGrid(const MInt, const MPI_Comm);
  ~StructuredGrid();

  MInt getGridMovingMethod() { return m_gridMovingMethod; };
  MBool isMovingGrid() { return m_movingGrid; };
  void moveCellPoints();
  void readGrid();
  void prepareReadGrid();
  void allocateMetricsAndJacobians();
  void gridDecomposition(MBool);
  void writeGrid(MString, MString);
  void writePartitionedGrid();
  void setCells(StructuredCell* cells) { m_cells = cells; };

  void exchangePoints(std::vector<std::unique_ptr<StructuredComm<nDim>>>&,
                      std::vector<std::unique_ptr<StructuredComm<nDim>>>&,
                      StructuredCommType);
  void gatherPoints(std::vector<std::unique_ptr<StructuredComm<nDim>>>&, StructuredCommType);
  void sendPoints(std::vector<std::unique_ptr<StructuredComm<nDim>>>&, StructuredCommType, std::vector<MPI_Request>&);
  void
  receivePoints(std::vector<std::unique_ptr<StructuredComm<nDim>>>&, StructuredCommType, std::vector<MPI_Request>&);
  void scatterPoints(std::vector<std::unique_ptr<StructuredComm<nDim>>>&, StructuredCommType);
  void periodicPointsChange(MFloat&, const MInt, const MInt);

  void setCellReference(StructuredCell*);

  void saveGrid();
  void saveCellJacobian();
  void computeCellCenterCoordinates();
  void computeMetrics();
  void computeJacobian();
  void computeSurfaceMetrics();
  void computeCornerMetrics();
  void computeModCornerMetrics();
  void computeCellMetrics();
  void computeCellJacobian();
  void computeCornerJacobian();
  void computeModCornerJacobian();
  void computeSurfaceJacobian();
  void computeSurfaceMetricsSingularity();
  void computeDxt(MFloat, MFloat*, MInt);

  void extrapolateGhostPointCoordinates();

  StructuredCell* m_cells;
  std::unique_ptr<StructuredDecomposition<nDim>> m_partition;

  // Singularities
  // class SingularInformation;
  SingularInformation* m_singularity = nullptr;
  MInt m_hasSingularity;

  MString m_uID;
  MString m_gridInputFileName; // is copy of m_gridInputFileName in the structuredBlck
  MInt m_gridFileId;
  MInt m_blockId;

  MInt* m_nPoints = nullptr;       // stores the maximum dimension of the partition with ghost points
  MInt* m_nActivePoints = nullptr; // stores the maximum dimension of the partition without ghost points
  MInt* m_nCells = nullptr;        // cell array dimension with ghost layer
  MInt* m_nActiveCells = nullptr;  // cell array dimension without ghost layer
  MInt* m_nOffsetCells = nullptr;
  MInt* m_nOffsetPoints = nullptr;
  MInt* m_nBlockCells = nullptr;
  MFloat** m_coordinates = nullptr;
  MFloat** m_oldCoordinates = nullptr;
  MFloat** m_initCoordinates = nullptr;
  MFloat** m_velocity = nullptr;
  MFloat** m_acceleration = nullptr;

  MInt m_noBlocks = 1;
  MInt m_noGhostLayers; // is a copy of m_noGhostLayers in the structuredBlck
  MInt m_totalNoCells;
  MInt m_noPoints;
  MInt m_noActiveCells;
  MInt m_noCells;

  MInt m_hasConnectionInfo; // not used yet, but can replace m_hasConnectionInfo in
                            // fvstructuredsolverwindowinfo
  MFloat* m_periodicDisplacements = nullptr;

  /**
   * \brief Returns the block id of the block in which the given domain is located
   * \param[in] domainId Domain ID
   * \return Block ID
   */
  MInt getBlockId(MInt domainId_) { return m_partition->getBlockIdFromPartition(domainId_); };

  /**
   * \brief Returns the block id of the block in which the own domain
   * \return Block ID
   */
  MInt getMyBlockId() { return m_partition->getBlockIdFromPartition(domainId()); };

  /**
   * \brief Returns the offset in the given dimension of the own domain inside the block
   * \param[in] dim Dimension
   * \return offset
   */
  MInt getMyOffset(MInt dim) { return m_partition->getPartitionOffset(domainId(), dim); };

  /**
   * \brief Returns the offset in the given dimension inside the block of the own domain
   * \param[in] domainId_ Domain ID
   * \param[in] dim Dimension
   * \return offset
   */
  MInt getOffset(MInt domainId_, MInt dim) { return m_partition->getPartitionOffset(domainId_, dim); };

  /**
   * \brief Returns the number of active points in the given dimension (without ghost-cells) of the own domain
   * \param[in] dim Dimension
   * \return number of active points
   */
  MInt getMyActivePoints(MInt dim) { return m_partition->getPartitionSize(domainId(), dim); };

  /**
   * \brief Returns the number of active points in the given dimension (without ghost-cells) of the given domain id
   * \param[in] domainId Domain ID
   * \param[in] dim Dimension
   * \return number of active points
   */
  MInt getActivePoints(MInt domainId_, MInt dim) { return m_partition->getPartitionSize(domainId_, dim); };

  /**
   * \brief Returns the number of total block cells in the given dimension for the given block id
   * \param[in] blockId Block ID
   * \param[in] dim Dimension
   * \return number of block cells
   */
  MInt getBlockNoCells(MInt blockId_, MInt dim) { return m_partition->getBlockSize(blockId_, dim); }

  /**
   * \brief Returns the number of total block points in the given dimension for the given block id
   * \param[in] blockId Block ID
   * \param[in] dim Dimension
   * \return number of block points
   */
  MInt getBlockNoPoints(MInt blockId_, MInt dim) { return (m_partition->getBlockSize(blockId_, dim) + 1); }

  /**
   * \brief Returns the number of total block cells in the given dimension for the own block
   * \param[in] dim Dimension
   * \return number of block cells
   */
  MInt getMyBlockNoCells(MInt dim) { return m_partition->getBlockSize(getMyBlockId(), dim); }

  /**
   * \brief Returns the number of total block points in the given dimension for the own block
   * \param[in] dim Dimension
   * \return number of block points
   */
  MInt getMyBlockNoPoints(MInt dim) { return m_partition->getBlockSize(getMyBlockId(), dim) + 1; }

  /**
   * \brief Returns the total number of blocks
   * \return number of blocks
   */
  MInt getNoBlocks() { return m_noBlocks; };

  /**
   * \brief Return the MPI communicator used by this grid
   */
  constexpr MPI_Comm mpiComm() const { return m_mpiComm; }

  /**
   * \brief Return the solver id to which this grid belongs
   */
  constexpr MInt solverId() const { return m_solverId; }

  /**
   * \brief Return the total number of domains (total number of ranks in current MPI communicator)
   */
  MInt noDomains() const { return m_noDomains; }

  /**
   * \brief Return the domainId (rank)
   */
  MInt domainId() const { return m_domainId; }

 protected:
  /**
   * \brief Computes the 3D cross product
   * \param[in] result Pointer to result array
   * \param[in] result Pointer to first input array
   * \param[in] result Pointer to second input array
   */
  inline void crossProduct(MFloat* result, MFloat* vec1, MFloat* vec2) {
    result[xsd] = vec1[ysd] * vec2[zsd] - vec1[zsd] * vec2[ysd];
    result[ysd] = vec1[zsd] * vec2[xsd] - vec1[xsd] * vec2[zsd];
    result[zsd] = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
  }

  /**
   * \brief Computes the 2D cross product
   * \param[in] result Pointer to result array
   * \param[in] result Pointer to first input array
   * \param[in] result Pointer to second input array
   */
  inline MFloat crossProduct(MFloat vec1[2], MFloat vec2[2]) {
    MFloat result = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
    return result;
  }

  /**
   * \brief Compute the point id of the point that has the offset (incI, inJ, incK) to the given point origin
   * \param[in] origin Origin point
   * \param[in] incI Increment in i-direction
   * \param[in] incJ Increment in j-direction
   * \param[in] incK Increment in k-direction
   * \return Point ID
   */
  inline MInt getPointIdFromPoint(const MInt origin, const MInt incI, const MInt incJ, const MInt incK) {
    return origin + incI + incJ * m_nPoints[2] + incK * m_nPoints[2] * m_nPoints[1];
  }

  /**
   * \brief Compute the point id of the point that has the offset (incI, inJ) to the given point origin
   * \param[in] origin Origin point
   * \param[in] incI Increment in i-direction
   * \param[in] incJ Increment in j-direction
   * \return Point ID
   */
  inline MInt getPointIdFromPoint(const MInt origin, const MInt incI, const MInt incJ) {
    return origin + incI + incJ * m_nPoints[1];
  }

  /**
   * \brief Compute the cell id of the cell that has the offset (incI, inJ, incK) to the given cell origin
   * \param[in] origin Origin point
   * \param[in] incI Increment in i-direction
   * \param[in] incJ Increment in j-direction
   * \param[in] incK Increment in k-direction
   * \return Cell ID
   */
  inline MInt getCellIdFromCell(const MInt origin, const MInt incI, const MInt incJ, const MInt incK) {
    return origin + incI + incJ * m_nCells[2] + incK * m_nCells[2] * m_nCells[1];
  }

  /**
   * \brief Compute the cell id of the cell that has the offset (incI, inJ) to the given cell origin
   * \param[in] origin Origin point
   * \param[in] incI Increment in i-direction
   * \param[in] incJ Increment in j-direction
   * \return Cell ID
   */
  inline MInt getCellIdFromCell(MInt origin, MInt incI, MInt incJ) { return origin + incI + incJ * m_nCells[1]; }

  /**
   * \brief Compute the lower point id for cell (i,j,k)
   * \param[in] i Cell i-coordinate
   * \param[in] j Cell j-coordinate
   * \param[in] k Cell k-coordinate
   * \return Point ID
   */
  inline MInt getPointIdFromCell(const MInt i, const MInt j, const MInt k) {
    return i + (k * (m_nCells[1] + 1) + j) * (m_nCells[2] + 1);
  }

  /**
   * \brief Compute the lower point id for cell (i,j)
   * \param[in] i Cell i-coordinate
   * \param[in] j Cell j-coordinate
   * \return Point ID
   */
  inline MInt getPointIdFromCell(const MInt i, const MInt j) { return i + (j * (m_nCells[1] + 1)); }

  /**
   * \brief Compute cell ID for given (i,j,k)
   * \param[in] i Cell i-coordinate
   * \param[in] j Cell j-coordinate
   * \param[in] k Cell k-coordinate
   * \return Cell ID
   */
  inline MInt cellIndex(const MInt i, const MInt j, const MInt k) { return i + (j + k * m_nCells[1]) * m_nCells[2]; }

  /**
   * \brief Compute cell ID for given (i,j)
   * \param[in] i Cell i-coordinate
   * \param[in] j Cell j-coordinate
   * \return Cell ID
   */
  inline MInt cellIndex(const MInt i, const MInt j) { return i + (j * m_nCells[1]); }

  /**
   * \brief Compute point ID for given (i,j,k)
   * \param[in] i Point i-coordinate
   * \param[in] j Point j-coordinate
   * \param[in] k Point k-coordinate
   * \return Point ID
   */
  inline MInt pointIndex(const MInt i, const MInt j, const MInt k) { return i + (j + k * m_nPoints[1]) * m_nPoints[2]; }

  /**
   * \brief Compute point ID for given (i,j)
   * \param[in] i Point i-coordinate
   * \param[in] j Point j-coordinate
   * \return Point ID
   */
  inline MInt pointIndex(const MInt i, const MInt j) { return i + (j * m_nPoints[1]); }

  inline MInt surfId(const MInt point, const MInt isd, const MInt dim) { return point + (isd + 3 * dim) * 9; }

  static const MInt xsd = 0;
  static const MInt ysd = 1;
  static const MInt zsd = 2;

 private:
  const MInt m_solverId;
  const MPI_Comm m_mpiComm;
  MInt m_domainId{};
  MInt m_noDomains{};

  MBool m_movingGrid;
  MInt m_gridMovingMethod;

  // repartitionFile
  MBool m_readDecompositionFromFile;


  inline std::array<MInt, nDim> pointBegin(const MInt plus) {
    std::array<MInt, nDim> result{};
    std::fill_n(result.begin(), nDim, plus);
    return result;
  }

  inline std::array<MInt, nDim> cellBegin(const MInt plus) {
    std::array<MInt, nDim> result{};
    std::fill_n(result.begin(), nDim, plus);
    return result;
  }

  inline std::array<MInt, nDim> cellEnd(const MInt minus) {
    std::array<MInt, nDim> result{};
    for(MInt i = 0; i < nDim; ++i) {
      result[i] = m_nCells[nDim - 1 - i] - minus;
    }
    return result;
  }

  inline std::array<MInt, nDim> pointEnd(const MInt minus) {
    std::array<MInt, nDim> result{};
    for(MInt i = 0; i < nDim; ++i) {
      result[i] = m_nPoints[nDim - 1 - i] - minus;
    }
    return result;
  }
};


#endif
