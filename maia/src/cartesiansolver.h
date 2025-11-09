// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CARTESIANSOLVER_H_
#define CARTESIANSOLVER_H_

#include <algorithm>

#include "COMM/mpioverride.h"
#include "FV/fvcartesiancellcollector.h"
#include "GEOM/geometry.h"
#include "GEOM/geometryanalytic.h"
#include "GEOM/geometryelement.h"
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridpoint.h"
#include "GRID/cartesiangridproxy.h"
#include "GRID/cartesianpointbasedcell.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "IO/parallelio.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "property.h"
#include "solver.h"

// Forward declarations
template <MInt nDim>
class CartesianGrid;

namespace maia {

struct PatchRefinement {
  // refinement patch
  MInt noLocalPatchRfnLvls() { return (MInt)localRfnLevelMethods.size(); };
  MInt noPatchesPerLevel(MInt addLevel) { return (MInt)localRfnLevelMethods[addLevel].size(); };

  std::vector<MString> localRfnLevelMethods;

  MInt* localRfnLevelPropertiesOffset = nullptr;
  MInt* noLocalRfnPatchProperties = nullptr;
  MFloat** localRfnPatchProperties = nullptr;
};

template <MInt nDim, class SolverType>
class CartesianSolver : public Solver {
 public:
  // Types
  using Grid = CartesianGrid<nDim>;
  using GridProxy = typename maia::grid::Proxy<nDim>;
  using Geom = Geometry<nDim>;
  using TreeProxy = maia::grid::tree::TreeProxy<nDim>;
  // using SolverCell = FvCell;
  using Cell = maia::grid::tree::Cell;

  // Methods
  // Constructor takes solverId, reference to grid
  CartesianSolver(const MInt solverId, GridProxy& gridProxy_, const MPI_Comm comm, const MBool checkActive = false);

  /// Read-only accessors for grid data
  MInt minLevel() const { return grid().minLevel(); }
  MInt maxLevel() const { return grid().maxLevel(); }
  MInt maxNoGridCells() const { return grid().maxNoCells(); }
  MInt maxRefinementLevel() const { return grid().maxRefinementLevel(); }
  MInt maxUniformRefinementLevel() const { return grid().maxUniformRefinementLevel(); }
  MInt noNeighborDomains() const { return grid().noNeighborDomains(); }
  MInt neighborDomain(const MInt id) const { return grid().neighborDomain(id); }
  MLong domainOffset(const MInt id) const { return grid().domainOffset(id); }
  MInt noHaloLayers() const { return grid().noHaloLayers(); }
  MInt noHaloCells(const MInt domainId) const { return grid().noHaloCells(domainId); }
  MInt haloCellId(const MInt domainId, const MInt cellId) const { return grid().haloCell(domainId, cellId); }
  MInt noWindowCells(const MInt domainId) const { return grid().noWindowCells(domainId); }
  MInt windowCellId(const MInt domainId, const MInt cellId) const { return grid().windowCell(domainId, cellId); }
  MString gridInputFileName() const { return grid().gridInputFileName(); }
  MFloat reductionFactor() const { return grid().reductionFactor(); }
  MFloat centerOfGravity(const MInt dir) const { return grid().centerOfGravity(dir); }
  MInt neighborList(const MInt cellId, const MInt dir) const { return grid().neighborList(cellId, dir); }
  const MLong& localPartitionCellGlobalIds(const MInt cellId) const {
    return grid().localPartitionCellGlobalIds(cellId);
  }
  MLong localPartitionCellOffsets(const MInt index) const { return grid().localPartitionCellOffsets(index); }
  MInt noMinCells() const { return grid().noMinCells(); }
  MInt minCell(const MInt id) const { return grid().minCell(id); }
  const MInt& haloCell(const MInt domainId, const MInt cellId) const { return grid().haloCell(domainId, cellId); }
  const MInt& windowCell(const MInt domainId, const MInt cellId) const { return grid().windowCell(domainId, cellId); }

  MBool isActive() const override { return grid().isActive(); }

  // Grid and tree proxy accessors
  constexpr GridProxy& grid() const { return m_gridProxy; }
  GridProxy& grid() { return m_gridProxy; }

  virtual void sensorDerivative(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                                std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorDivergence(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                                std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorTotalPressure(std::vector<std::vector<MFloat>>& /*unused*/,
                                   std::vector<std::bitset<64>>& /*unused*/, std::vector<MFloat>& /*unused*/,
                                   MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorEntropyGrad(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                                 std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorEntropyQuot(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                                 std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorVorticity(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                               std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorInterface(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                               std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  void sensorLimit(std::vector<std::vector<MFloat>>& /*sensors*/, std::vector<std::bitset<64>>& /*sensorCellFlag*/,
                   std::vector<MFloat>& /*sensorWeight*/, MInt /*sensorOffset*/, MInt /*sen*/,
                   std::function<MFloat(MInt)> /*value*/, const MFloat /*limit*/, const MInt* /*range*/,
                   const MBool /*refineDiagonals*/, const MBool allowCoarsening = true);
  void sensorSmooth(std::vector<std::vector<MFloat>>& /*sensors*/, std::vector<std::bitset<64>>& /*sensorCellFlag*/,
                    std::vector<MFloat>& /*sensorWeight*/, MInt /*sensorOffset*/, MInt /*sen*/);
  void sensorBand(std::vector<std::vector<MFloat>>& /*sensors*/, std::vector<std::bitset<64>>& /*sensorCellFlag*/,
                  std::vector<MFloat>& /*sensorWeight*/, MInt /*sensorOffset*/, MInt /*sen*/);
  virtual void sensorMeanStress(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                                std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorParticle(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                              std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorSpecies(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                             std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/) {
    TERMM(1, "Not implemented for this solver");
  };
  virtual void sensorPatch(std::vector<std::vector<MFloat>>& sensor, std::vector<std::bitset<64>>& sensorCellFlag,
                           std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) {
    patchRefinement(sensor, sensorCellFlag, sensorWeight, sensorOffset, sen);
  };
  virtual void sensorCutOff(std::vector<std::vector<MFloat>>& /*unused*/, std::vector<std::bitset<64>>& /*unused*/,
                            std::vector<MFloat>& /*unused*/, MInt /*unused*/, MInt /*unused*/)
  {
    TERMM(1, "Not implemented for this solver");
  };
  void saveSensorData(const std::vector<std::vector<MFloat>>& sensors, const MInt& level, const MString& gridFileName,
                      const MInt* const recalcIds) override;

  // Asserts that cellId belongs to a valid grid cell
  void assertValidGridCellId(const MInt) const {}
  /// \brief Returns the grid parent id of the cell \p cellId
  MLong c_parentId(const MInt cellId) const { return solver().grid().tree().parent(cellId); }
  /// \brief Returns the grid neighbor id of the grid cell \p cellId \p dir
  MLong c_neighborId(const MInt cellId, const MInt dir) const { return solver().grid().tree().neighbor(cellId, dir); }
  MInt c_noCells() const { return grid().tree().size(); }
  MInt c_level(const MInt cellId) const { return grid().tree().level(cellId); }

  // \brief Return global grid id
  MLong c_globalGridId(const MInt cellId) { return grid().raw().a_globalId(grid().tree().solver2grid(cellId)); }

  // General Exchange Functions
  template <typename T>
  void exchangeData(T* data, const MInt dataBlockSize = 1);
  template <typename T>
  void exchangeLeafData(std::function<T&(MInt, MInt)> data, const MInt noDat = 1);


  template <class G, class S, class M>
  void exchangeSparseLeafValues(G getData, S setData, const MInt dataSize, M cellMapping);

  template <typename T>
  void exchangeAzimuthalPer(T* data, MInt dataBlockSize = 1, MInt firstBlock = 0);

  template <typename T>
  void collectVariables(T* variablesIn, ScratchSpace<T>& variablesOut, const std::vector<MString>& variablesNameIn,
                        std::vector<MString>& variablesNameOut, const MInt noVars, const MInt noCells,
                        const MBool reverseOrder = false);
  template <typename T>
  void collectVariables(T** variablesIn, ScratchSpace<T>& variablesOut, const std::vector<MString>& variablesNameIn,
                        std::vector<MString>& variablesNameOut, const MInt noVars, const MInt noCells);

  void saveGridFlowVars(const MChar* fileName, const MChar* gridFileName, const MInt noTotalCells,
                        const MInt noInternal, MFloatScratchSpace& dbVariables, std::vector<MString>& dbVariablesName,
                        MInt noDbVars, MIntScratchSpace& idVariables, std::vector<MString>& idVariablesName,
                        MInt noIdVars, MFloatScratchSpace& dbParameters, std::vector<MString>& dbParametersName,
                        MIntScratchSpace& idParameters, std::vector<MString>& idParametersName, MInt* recalcIds,
                        MFloat time);

  template <typename T>
  void collectParameters(T, ScratchSpace<T>&, const MChar*, std::vector<MString>&);

  void calcRecalcCellIdsSolver(const MInt* const recalcIdsTree, MInt& noCells, MInt& noInternalCellIds,
                               std::vector<MInt>& recalcCellIdsSolver, std::vector<MInt>& reorderedCellIds);

 protected:
  // Grid/geometry-related methods
  // stl-boundary related:
  void identifyBoundaryCells(MBool* const isInterface, const std::vector<MInt>& bndCndIds = std::vector<MInt>());
  void identifyBoundaryCells();
  MBool cellIsOnGeometry(MInt cellId, Geometry<nDim>* geom);
  void setBoundaryDistance(const MBool* const interfaceCell, const MFloat* const outerBandWidth,
                           MFloatScratchSpace& distance);
  void markSurrndCells(MIntScratchSpace& inList, const MInt bandWidth, const MInt level,
                       const MBool refineDiagonals = true);
  void receiveWindowTriangles();

  // adaptation related
  void compactCells();
  MInt createCellId(const MInt gridCellId);
  void removeCellId(const MInt cellId);

  MInt inCell(const MInt cellId, MFloat* point, MFloat fac = F1);

  MInt setUpInterpolationStencil(const MInt cellId, MInt*, const MFloat*, std::function<MBool(MInt, MInt)>,
                                 MBool allowIncompleteStencil);

  template <MBool cartesianInterpolation>
  MFloat interpolateFieldData(MInt*, MFloat*, MInt varId, std::function<MFloat(MInt, MInt)> scalarField,
                              std::function<MFloat(MInt, MInt)> coordinate);

  MFloat leastSquaresInterpolation(MInt*, MFloat*, MInt varId, std::function<MFloat(MInt, MInt)> scalarField,
                                   std::function<MFloat(MInt, MInt)> coordinate);

  void checkNoHaloLayers();

  void mapInterpolationCells(std::map<MInt, MInt>& cellMap);
  void setHaloCellsOnInactiveRanks();

  void patchRefinement(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen);

  void reOrderCellIds(std::vector<MInt>& reOrderedCells);
  void recomputeGlobalIds(std::vector<MInt>&, std::vector<MLong>&);

  void extractPointIdsFromGrid(Collector<PointBasedCell<nDim>>*&,
                               Collector<CartesianGridPoint<nDim>>*&,
                               const MBool,
                               const std::map<MInt, MInt>&,
                               MInt levelThreshold = 999999,
                               MFloat* bBox = nullptr,
                               MBool levelSetMb = false) const;

 private:
  // Methods
  // CRTP accessors
  SolverType& solver() { return static_cast<SolverType&>(*this); }
  constexpr const SolverType& solver() const { return static_cast<const SolverType&>(*this); }

 protected:
  MInt* m_rfnBandWidth;
  MInt m_noSensors;
  MInt m_adaptationInterval;
  MInt m_adaptationStep = F0;
  std::vector<MInt> m_maxSensorRefinementLevel;
  std::vector<MFloat> m_sensorWeight;
  std::vector<MFloat> m_sensorDerivativeVariables;
  MBool m_adaptation;
  MBool m_adapts;
  MInt m_noInitialSensors;
  MBool m_resTriggeredAdapt = false;
  MInt m_noSmoothingLayers = -1;
  MInt m_sensorBandAdditionalLayers;

  // adaptation
  MBool m_sensorInterface;
  MBool m_sensorParticle;
  std::vector<MString> m_sensorType;
  MInt* m_recalcIds = nullptr;

  using fun = void (CartesianSolver<nDim, SolverType>::*)(std::vector<std::vector<MFloat>>&,
                                                          std::vector<std::bitset<64>>&, std::vector<MFloat>&, MInt,
                                                          MInt);
  std::vector<fun> m_sensorFnPtr;

  // necessary for generalised interpolation function!
  MInt m_maxNoSets = -1;

  // Azimuthal communication
  std::vector<MFloat> m_azimuthalCartRecCoord;

  // Holds the reverse-direction-Information
  const MInt m_revDir[6] = {1, 0, 3, 2, 5, 4};
  const MInt m_noDirs = 2 * nDim;

 private:
  void readPatchProperties();
  void initAdaptation();
  void addChildren(std::vector<MInt>& reOrderedIds, const MInt parentId);


  // Data
  MString m_gridCutTest{};
  GridProxy& m_gridProxy;

  PatchRefinement* m_patchRefinement = nullptr;
  MBool m_testPatch = false;
  MInt m_patchStartTimeStep = -1;
  MInt m_patchStopTimeStep = -1;
};

template <MInt nDim, class SolverType>
CartesianSolver<nDim, SolverType>::CartesianSolver(const MInt solverId_, GridProxy& gridProxy_, const MPI_Comm comm,
                                                   const MBool checkActive)
  : Solver(solverId_, comm, checkActive ? gridProxy_.isActive() : true), m_gridProxy(gridProxy_) {
  // Get necessary properties
  m_gridCutTest = "SAT";
  m_gridCutTest = Context::getSolverProperty<MString>("gridCutTest", solverId_, AT_, &m_gridCutTest);

  initAdaptation();
}

/** Identifies boundary cells by determining intersections betw. geometry and grid
 * Argument `isInterface` must point to an array with a size of at least
 * "solver().grid().tree().size()" elements
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::identifyBoundaryCells(MBool* const isInterface,
                                                              const std::vector<MInt>& bndCndIds) {
  // TODO labels:GEOM,DLB add mode to tag only halos during DLB if isInterface property is communicated?
  TRACE();

  MFloat target[6] = {0, 0, 0, 0, 0, 0};

  // Check all cells for intersection with geometry
  const MInt noCells = solver().grid().tree().size();

  for(MInt i = 0; i < noCells; i++) {
    const MFloat cellHalfLength = solver().grid().cellLengthAtLevel(solver().grid().tree().level(i) + 1);
    for(MInt j = 0; j < nDim; j++) {
      target[j] = solver().a_coordinate(i, j) - cellHalfLength;
      target[j + nDim] = solver().a_coordinate(i, j) + cellHalfLength;
    }

    std::vector<MInt> nodeList;
    // TODO labels:GEOM why not read gridCutTest once in the geometry and decide there which getIntersectionElements
    // function to use?
    if(m_gridCutTest == "SAT") {
      solver().geometry().getIntersectionElements(target, nodeList, cellHalfLength, &solver().a_coordinate(i, 0));
    } else {
      solver().geometry().getIntersectionElements(target, nodeList);
    }

    isInterface[i] = false;

    const MInt noNodes = nodeList.size();
    if(noNodes > 0 && bndCndIds.empty()) {
      isInterface[i] = true;
    } else if(noNodes > 0) {
      for(MInt n = 0; n < noNodes; n++) {
        for(const auto& bnd : bndCndIds) {
          if(bnd == solver().geometry().elements[nodeList[n]].m_bndCndId) {
            isInterface[i] = true;
          }
        }
      }
    }
  }
}


///  Directly sets "a_isInterface" in solver
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-16
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::identifyBoundaryCells() {
  TRACE();

  const MInt noCells = solver().grid().tree().size();
  MBoolScratchSpace isInterface(noCells, AT_, "isInterface");
  identifyBoundaryCells(&isInterface[0]);
  for(MInt i = 0; i < noCells; i++) {
    solver().a_isInterface(i) = isInterface[i];
  }
}


/** \brief checks whether a cell lies on a certain geometry
 *         copied the essential part from identifyBoundaryCells
 *
 * \author Thomas Schilden
 * \date 18.10.2016
 *
 * \param[in]  cellId
 * \param[in]  geometry
 * \param[out] true/false
 *
 **/
template <MInt nDim, class SolverType>
MBool CartesianSolver<nDim, SolverType>::cellIsOnGeometry(MInt cellId, Geometry<nDim>* geom) {
  TRACE();

  MFloat target[6] = {0, 0, 0, 0, 0, 0};
  const MFloat cellHalfLength = solver().grid().cellLengthAtLevel(solver().a_level(cellId) + 1);
  std::vector<MInt> nodeList;

  for(MInt dim = 0; dim < nDim; dim++) {
    target[dim] = solver().a_coordinate(cellId, dim) - cellHalfLength;
    target[dim + nDim] = solver().a_coordinate(cellId, dim) + cellHalfLength;
  }

  if(m_gridCutTest == "SAT") {
    geom->getIntersectionElements(target, nodeList, cellHalfLength, &solver().a_coordinate(cellId, 0));
  } else {
    geom->getIntersectionElements(target, nodeList);
  }

  return nodeList.size() > 0;
}


/** \brief Receives triangles from neighbors contained in their window cells and inserts them locally
 *
 * \author Andreas Lintermann
 * \date 25.09.2015
 *
 * The algorithm does the following:
 *
 *  1. Collect the traingles contained in my window cells by running over
 *     all window cells provided by the input array windowCells. Store each
 *     triangle in a set of pairs with the triangle id and the original id.
 *  2. Send and receive the unique triangle ids, i.e., the original ids for
 *     comparison, which triangles to keep in the array created in 1. Remove
 *     those not needed, i.e., remove those which are already available on
 *     the neighbor.
 *  3. Inform all neighbors how many triangles will be sent and the send
 *     an receive the according triangles.
 *  4. Insert the received triangles in the geometry.
 *  5. Update the tree and the bounding box of the geometry.
 *
 * \param[in] neighbors array containing all neighbor domain ids
 * \param[in] noNeighbors number of the domain neighbors
 * \param[in] windowCells the ids of the window cells per neighbor domain
 * \param[in] noWindowCells the number of window cells per neighbor domain
 *
 **/
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::receiveWindowTriangles() {
  TRACE();

  using namespace std;

  m_log << "    - exchanging triangles " << std::endl;
  vector<set<pair<MInt, MInt>>> triangleIdsPerDomain;

  MPI_Request* mpi_request = nullptr;
  MPI_Status status;
  mAlloc(mpi_request, solver().grid().noNeighborDomains(), "mpi_request", AT_);

  m_log << "      * collecting window triangles" << std::endl;
  // 1. collect unique triangles per neighbor domain
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    set<pair<MInt, MInt>> triangles;
    for(MInt w = 0; w < (signed)solver().grid().noWindowCells(n); w++) {
      MInt currentId = solver().grid().windowCell(n, w);

      // on coarsest level
      if(solver().grid().tree().parent(currentId) < 0) {
        // a. Create target for check
        MFloat target[6];
        MFloat cellHalfLength = solver().grid().cellLengthAtLevel(solver().grid().tree().level(currentId) + 1);
        cellHalfLength += 0.005 * cellHalfLength;

        std::vector<MInt> nodeList;

        for(MInt j = 0; j < nDim; j++) {
          target[j] = solver().a_coordinate(currentId, j) - cellHalfLength;
          target[j + nDim] = solver().a_coordinate(currentId, j) + cellHalfLength;
        }

        // b. Do the intersection test
        solver().geometry().getIntersectionElements(target, nodeList, cellHalfLength,
                                                    &solver().a_coordinate(currentId, 0));

        for(MInt t = 0; t < (signed)nodeList.size(); t++) {
          triangles.insert(pair<MInt, MInt>(nodeList[t], solver().geometry().elements[nodeList[t]].m_originalId));
        }
      }
    }
    triangleIdsPerDomain.push_back(triangles);
  }

  // 2. let me know about the uniqueIds of my neighbors
  MIntScratchSpace numUniques(solver().grid().noNeighborDomains(), AT_, "numUniques");
  MInt myNumUniques = solver().geometry().m_uniqueOriginalTriId.size();
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MPI_Issend(&myNumUniques, 1, MPI_INT, solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &mpi_request[n], AT_,
               "myNumUniques");
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MPI_Recv(&numUniques[n], 1, MPI_INT, solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &status, AT_,
             "numUniques[n]");
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MPI_Wait(&mpi_request[n], &status, AT_);
  }

  MInt sumuniques = 0;
  MIntScratchSpace uniquesDomOff(solver().grid().noNeighborDomains(), AT_, "uniquesDomOff");
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    uniquesDomOff[n] = sumuniques;
    sumuniques += numUniques[n];
  }

  // send and receive the uniques
  MIntScratchSpace myUniques(myNumUniques, AT_, "myUniques");
  MInt ua = 0;
  for(auto u = solver().geometry().m_uniqueOriginalTriId.begin(); u != solver().geometry().m_uniqueOriginalTriId.end();
      ++u, ua++) {
    myUniques[ua] = *u;
  }
  MIntScratchSpace alluniques(sumuniques, AT_, "alluniques");

  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(!myUniques.empty())
      MPI_Issend(myUniques.getPointer(), myUniques.size(), MPI_INT, solver().grid().neighborDomain(n), 0,
                 MPI_COMM_WORLD, &mpi_request[n], AT_, "myUniques.getPointer()");
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(numUniques[n] > 0) {
      MPI_Recv(&alluniques[uniquesDomOff[n]], numUniques[n], MPI_INT, solver().grid().neighborDomain(n), 0,
               MPI_COMM_WORLD, &status, AT_, "alluniques[uniquesDomOff[n]]");
    }
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(numUniques[n] > 0) {
      MPI_Wait(&mpi_request[n], &status, AT_);
    }
  }

  m_log << "      * sending and receiving unique originalIds" << std::endl;
  m_log << "        + sending " << myUniques.size() << " originalIds to: ";
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    m_log << solver().grid().neighborDomain(n) << " ";
  }
  m_log << std::endl;
  m_log << "      * receiving unique originalIds" << std::endl;
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    m_log << "        + receiving from " << solver().grid().neighborDomain(n) << ": " << numUniques[n] << std::endl;
  }

  m_log << "      * removing doubles from sender list" << std::endl;

  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    // for searching
    set<MInt> un;
    for(MInt off = uniquesDomOff[n]; off < uniquesDomOff[n] + numUniques[n]; off++) {
      un.insert(alluniques[off]);
    }

    for(auto it = triangleIdsPerDomain[n].begin(); it != triangleIdsPerDomain[n].end();) {
      auto iun = un.find(get<1>(*it));
      if(iun != un.end()) {
        triangleIdsPerDomain[n].erase(it);
        it = triangleIdsPerDomain[n].begin();
      } else {
        ++it;
      }
    }
  }


  MIntScratchSpace toReceive(solver().grid().noNeighborDomains(), AT_, "toReceive");


  m_log << "      * sending numbers of triangles to send to other domains" << std::endl;
  // 3. send information how many triangles will be send
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MInt numtris = triangleIdsPerDomain[n].size();
    MPI_Issend(&numtris, 1, MPI_INT, solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &mpi_request[n], AT_,
               "numtris");
  }

  // receive information
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MPI_Recv(&toReceive[n], 1, MPI_INT, solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &status, AT_,
             "toReceive[n]");
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    MPI_Wait(&mpi_request[n], &status, AT_);
  }

  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(toReceive[n] > 0) {
      m_log << "        + receive from domain " << solver().grid().neighborDomain(n) << " : " << toReceive[n]
            << std::endl;
    }
    if(!triangleIdsPerDomain[n].empty()) {
      m_log << "        + send to domain " << solver().grid().neighborDomain(n) << "      : "
            << triangleIdsPerDomain[n].size() << std::endl;
    }
  }


  MInt sumofreceive = 0;
  MInt sumofsend = 0;
  MIntScratchSpace offsetreceive(solver().grid().noNeighborDomains(), AT_, "offsetreceive");
  MIntScratchSpace offsetsend(solver().grid().noNeighborDomains(), AT_, "offsetsend");
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    offsetsend[n] = sumofsend;
    sumofsend += triangleIdsPerDomain[n].size();
    offsetreceive[n] = sumofreceive;
    sumofreceive += toReceive[n];
  }

  m_log << "      * receive offsets:" << std::endl;
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(toReceive[n] > 0) {
      m_log << "        + from domain " << solver().grid().neighborDomain(n) << ": " << offsetreceive[n] << std::endl;
    }
  }

  m_log << "      * send offsets:" << endl;
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(!triangleIdsPerDomain[n].empty()) {
      m_log << "        + from domain " << solver().grid().neighborDomain(n) << ": " << offsetsend[n] << std::endl;
    }
  }


  // create array holding triangles to receive
  MInt trisize = nDim            // normal
                 + (nDim * nDim) // vertcies
                 + 3;            // segmentId, bndCndId, orininalId

  MFloatScratchSpace recTris(sumofreceive * trisize, AT_, "recTris");
  MFloatScratchSpace sndTris(sumofsend * trisize, AT_, "sndTris");

  m_log << "      * sending triangles" << std::endl;
  // send the triangles
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(triangleIdsPerDomain[n].empty()) {
      continue;
    }
    // create buffer

    MInt j = offsetsend[n] * trisize;
    for(const auto& tri : triangleIdsPerDomain[n]) {
      sndTris[j++] = (MFloat)solver().geometry().elements[get<0>(tri)].m_originalId;
      sndTris[j++] = (MFloat)solver().geometry().elements[get<0>(tri)].m_segmentId;
      sndTris[j++] = (MFloat)solver().geometry().elements[get<0>(tri)].m_bndCndId;

      // normal
      for(MInt d = 0; d < nDim; d++, j++) {
        sndTris[j] = solver().geometry().elements[get<0>(tri)].m_normal[d];
      }

      // vertcies
      for(MInt d1 = 0; d1 < nDim; d1++) {
        for(MInt d2 = 0; d2 < nDim; d2++, j++) {
          sndTris[j] = solver().geometry().elements[get<0>(tri)].m_vertices[d1][d2];
        }
      }
    }
    MPI_Issend(sndTris.getPointer() + (offsetsend[n] * trisize), trisize * triangleIdsPerDomain[n].size(), MPI_DOUBLE,
               solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &mpi_request[n], AT_,
               "sndTris.getPointer()+(offsetsend[n]*trisize)");
  }

  m_log << "      * receiving triangles" << endl;
  // receive triangles
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(toReceive[n] > 0) {
      MPI_Recv(recTris.getPointer() + (offsetreceive[n] * trisize), toReceive[n] * trisize, MPI_DOUBLE,
               solver().grid().neighborDomain(n), 0, MPI_COMM_WORLD, &status, AT_,
               "recTris.getPointer()+(offsetreceive[n]*trisize)");
    }
  }
  for(MInt n = 0; n < solver().grid().noNeighborDomains(); n++) {
    if(toReceive[n] > 0) {
      MPI_Wait(&mpi_request[n], &status, AT_);
    }
  }

  m_log << "      * inserting received triangles" << endl;

  // if(solver().geometry().GetNoElements() == 0)
  if(sumofreceive > 0) {
    solver().geometry().resizeCollector(sumofreceive);
  }


  // 4. insert the triangles
  for(MInt tri = 0; tri < sumofreceive * trisize; tri += trisize) {
    auto originalId = (MInt)recTris[tri];
    // check if triangle already in my domain
    // not found, great, insert the triangle
    if(solver().geometry().m_uniqueOriginalTriId.find(originalId) == solver().geometry().m_uniqueOriginalTriId.end()) {
      solver().geometry().addElement(&recTris[tri]);
    }
  }

  // 5. update tree
  solver().geometry().calculateBoundingBox();

  if(solver().geometry().m_debugParGeom && solver().geometry().GetNoElements() > 0) {
    solver().geometry().writeParallelGeometryVTK("allpluswindow");
  }

  if(solver().geometry().GetNoElements() > 0 && sumofreceive > 0) {
    solver().geometry().rebuildAdtTree();
  }
}


// --------------------------------------------------------------------------------------


/**
 * \brief Removes all holes in the cell collector and moves halo cells to the back of the collector
 * \author Lennart Schneiders
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::compactCells() {
  MIntScratchSpace oldCellId(m_gridProxy.maxNoCells(), AT_, "oldCellId");
  MIntScratchSpace isToDelete(m_gridProxy.maxNoCells(), AT_, "isToDelete");
  oldCellId.fill(-1);
  isToDelete.fill(0);
  for(auto& i : solver().m_freeIndices) {
    isToDelete[i] = 1;
  }
  solver().m_freeIndices.clear();

  if(grid().azimuthalPeriodicity()) {
    m_gridProxy.correctAzimuthalHaloCells();
  }

  m_gridProxy.resizeGridMap(solver().m_cells.size());

  // 0. some sanity checks
  ASSERT(solver().m_cells.size() == m_gridProxy.tree().size(),
         "m_cells: " + std::to_string(solver().m_cells.size()) + " tree: " + std::to_string(m_gridProxy.tree().size()));
  for(MInt cellId = 0; cellId < solver().m_cells.size(); cellId++) {
    ASSERT(isToDelete(cellId)
               || solver().a_isHalo(cellId) == solver().grid().raw().a_isHalo(grid().tree().solver2grid(cellId)),
           "");
  }
  for(MInt gridCellId = 0; gridCellId < solver().grid().raw().treeb().size(); gridCellId++) {
    if(solver().grid().raw().treeb().solver(gridCellId, solver().solverId())) {
      ASSERT(
          grid().tree().grid2solver(gridCellId) > -1 && grid().tree().grid2solver(gridCellId) < solver().m_cells.size(),
          std::to_string(gridCellId) + " " + std::to_string(grid().tree().grid2solver(gridCellId)) + " "
              + std::to_string(solver().m_cells.size()) + " " + std::to_string(solver().grid().raw().treeb().size()));
    } else {
      ASSERT(grid().tree().grid2solver(gridCellId) < 0,
             std::to_string(gridCellId) + " " + std::to_string(grid().tree().grid2solver(gridCellId)) + " "
                 + std::to_string(solver().m_cells.size()) + " "
                 + std::to_string(solver().grid().raw().treeb().size()));
    }
  }

  // 1. determine number of cells and internal cells
  MInt noCells = 0;
  MInt noInternalCells = 0;
  oldCellId.fill(-1);
  for(MInt cellId = 0; cellId < solver().m_cells.size(); cellId++) {
    if(isToDelete(cellId) != 0) {
      continue;
    }
    oldCellId(cellId) = cellId;
    if(!solver().a_isHalo(cellId)) {
      noInternalCells++;
    }
    noCells++;
  }

  if(solver().grid().raw().treeb().noSolvers() == 1) {
    ASSERT(noCells == solver().grid().raw().treeb().size(), "");
    ASSERT(noInternalCells == solver().grid().raw().m_noInternalCells, "");
  }

  // 2. remove holes created by previously deleted cells and move halo cells to the back
  MInt otherId = solver().m_cells.size() - 1;
  for(MInt cellId = 0; cellId < noInternalCells; cellId++) {
    if(isToDelete(cellId) || solver().a_isHalo(cellId)) {
      while(isToDelete(otherId) || solver().a_isHalo(otherId)) {
        otherId--;
      }
      ASSERT(cellId < otherId, "");

      solver().swapCells(cellId, otherId);
      grid().swapSolverIds(cellId, otherId);
      std::swap(oldCellId(cellId), oldCellId(otherId));
      std::swap(isToDelete(cellId), isToDelete(otherId));

      ASSERT(grid().tree().solver2grid(cellId) > -1, "");
      ASSERT(grid().tree().solver2grid(otherId) < 0 || solver().a_isHalo(otherId), "");
      ASSERT(!solver().a_isHalo(cellId), "");
      ASSERT(isToDelete(otherId) || solver().a_isHalo(otherId), "");
    }
    ASSERT(!solver().a_isHalo(cellId) && !isToDelete(cellId), "");
  }

  // 3. remove holes in the range of halo cells
  otherId = solver().m_cells.size() - 1;
  for(MInt cellId = noInternalCells; cellId < noCells; cellId++) {
    if(isToDelete(cellId) != 0) {
      while(isToDelete(otherId) != 0) {
        otherId--;
      }
      ASSERT(cellId < otherId, "");
      ASSERT(otherId >= noCells, "");
      ASSERT(solver().a_isHalo(otherId) && !isToDelete(otherId), "");
      solver().swapCells(cellId, otherId);
      grid().swapSolverIds(cellId, otherId);
      std::swap(oldCellId(cellId), oldCellId(otherId));
      std::swap(isToDelete(cellId), isToDelete(otherId));

      ASSERT(grid().tree().solver2grid(cellId) > -1, "");
      ASSERT(grid().tree().solver2grid(otherId) < 0, "");
    }
    ASSERT(solver().a_isHalo(cellId) && !isToDelete(cellId), "");
  }

  for(MInt cellId = noInternalCells; cellId < noCells; cellId++) {
    ASSERT(grid().tree().solver2grid(cellId) > -1
               && grid().tree().solver2grid(cellId) < solver().grid().raw().treeb().size(),
           "");
    ASSERT(solver().a_isHalo(cellId) == solver().grid().raw().a_isHalo(grid().tree().solver2grid(cellId)), "");
  }

  solver().m_cells.size(noCells);
  m_gridProxy.resizeGridMap(solver().m_cells.size());
  ASSERT(solver().m_cells.size() == m_gridProxy.tree().size(), "");


  if(solver().grid().raw().treeb().noSolvers() == 1) {
    for(MInt gridCellId = 0; gridCellId < noCells; gridCellId++) {
      ASSERT(grid().tree().solver2grid(gridCellId) == gridCellId, "");
      ASSERT(grid().tree().grid2solver(gridCellId) == gridCellId, "");
    }
  }


  /*
    MIntScratchSpace newCellId(m_gridProxy.maxNoCells(), AT_, "newCellId");
    newCellId.fill(-1);
    for( MInt cellId = 0; cellId < noCells; cellId++ ) {
      if ( oldCellId( cellId ) < 0 ) continue;
      newCellId( oldCellId( cellId ) ) = cellId;
    }
    for(MInt i=0; i < noNeighborDomains(); i++) {
      for ( MInt j = 0; j < (signed)m_gridProxy.m_haloCells[i].size(); j++ ) {
        MInt cellId = m_gridProxy.m_haloCells[i][j];
        if ( newCellId(cellId) > -1 ) {
          m_gridProxy.m_haloCells[i][j] = newCellId(cellId);
        }
      }
      for ( MInt j = 0; j < (signed)m_gridProxy.m_windowCells[i].size(); j++ ) {
        MInt cellId = m_gridProxy.m_windowCells[i][j];
        if ( newCellId(cellId) > -1 ) {
          m_gridProxy.m_windowCells[i][j] = newCellId(cellId);
        }
      }
    }
  */

  //  m_gridProxy.tree().size(noCells);
  //  m_gridProxy.m_noInternalCells = noInternalCells;
}


// --------------------------------------------------------------------------------------


/**
 * \brief
 * \author Lennart Schneiders
 */
template <MInt nDim, class SolverType>
MInt CartesianSolver<nDim, SolverType>::createCellId(const MInt gridCellId) {
  MInt solverCellId = -1;
  if(solver().m_freeIndices.size() > 0) {
    auto it = solver().m_freeIndices.begin();
    solverCellId = *(it);
    solver().m_freeIndices.erase(it);
    m_gridProxy.resizeGridMap(solver().m_cells.size());
  } else {
    solverCellId = solver().m_cells.size();
    solver().m_cells.append();
    m_gridProxy.resizeGridMap(solver().m_cells.size());
  }
  ASSERT(solverCellId > -1 && solverCellId < solver().m_cells.size(), "");
  if(!g_multiSolverGrid) {
    ASSERT(solverCellId == gridCellId, std::to_string(solverCellId) + " " + std::to_string(gridCellId));
  }

  solver().a_resetPropertiesSolver(solverCellId);
  solver().a_isHalo(solverCellId) = solver().grid().raw().a_isHalo(gridCellId);

  grid().setSolver2grid(solverCellId, gridCellId);
  grid().setGrid2solver(gridCellId, solverCellId);

  return solverCellId;
}


/**
 * \brief
 * \author Lennart Schneiders
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::removeCellId(const MInt cellId) {
  const MInt gridCellId = grid().tree().solver2grid(cellId);
  ASSERT(gridCellId > -1 && gridCellId < solver().grid().raw().treeb().size(), "");

  solver().a_resetPropertiesSolver(cellId);

  grid().setSolver2grid(cellId, std::numeric_limits<MInt>::min());
  grid().setGrid2solver(gridCellId, std::numeric_limits<MInt>::min());
  solver().grid().raw().treeb().solver(gridCellId, solver().solverId()) = false;

  if(cellId == (solver().m_cells.size() - 1)) {
    solver().m_cells.size(solver().m_cells.size() - 1);
    m_gridProxy.resizeGridMap(solver().m_cells.size());
  } else {
    solver().m_freeIndices.insert(cellId);
  }

  ASSERT(g_multiSolverGrid || cellId == gridCellId, "");
}

///  transverses over all neighboring cells for a specified length
//   and sets the approximate distance to the interfaceCells!
///  For refinemnt (or other purposes!)
///
/// \author Tim Wegmann
/// \date 2018-11-14

template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::setBoundaryDistance(const MBool* const interfaceCell,
                                                            const MFloat* const outerBandWidth,
                                                            MFloatScratchSpace& distance) {
  TRACE();

  distance.fill(std::numeric_limits<MFloat>::max());
  std::vector<MInt> currentLayer;
  MIntScratchSpace onCurrentLayer(solver().a_noCells(), AT_, "onCurrentLayer");
  MIntScratchSpace accessed(solver().a_noCells(), AT_, "accessed");
  onCurrentLayer.fill(0);
  accessed.fill(0);

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(interfaceCell[cellId]) {
      currentLayer.push_back(cellId);
      onCurrentLayer(cellId) = 1;
      accessed(cellId) = 1;
      distance(cellId) = F0;
    }
  }
  MInt proceed = 1;
  while(proceed != 0) {
    proceed = 0;
    solver().exchangeData(distance.data());
    solver().exchangeData(onCurrentLayer.data());
    for(MInt i = 0; i < solver().noNeighborDomains(); i++) {
      for(MInt c = 0; c < solver().noHaloCells(i); c++) {
        if(onCurrentLayer(solver().haloCellId(i, c))) {
          currentLayer.push_back(solver().haloCellId(i, c));
          accessed(solver().haloCellId(i, c)) = 1;
        }
      }
    }
    onCurrentLayer.fill(0);
    std::vector<MInt> currentLayerBak(currentLayer);
    for(auto& cellId : currentLayerBak) {
      for(MInt dir = 0; dir < solver().m_noDirs; dir++) {
        if(solver().a_hasNeighbor(cellId, dir)) {
          MInt nghbrId = solver().c_neighborId(cellId, dir);
          if(solver().c_isLeafCell(cellId)) {
            MFloat dx = F1B2 * (solver().c_cellLengthAtCell(cellId) + solver().c_cellLengthAtCell(nghbrId));
            distance(nghbrId) = mMin(distance(nghbrId), distance(cellId) + dx);
            if(accessed(nghbrId) == 0) {
              currentLayer.push_back(nghbrId);
              accessed(nghbrId) = 1;
              onCurrentLayer(nghbrId) = 1;
              if(distance(nghbrId) < outerBandWidth[mMax(solver().minLevel(), solver().a_level(nghbrId) - 1)]
                                         + F2 * solver().c_cellLengthAtCell(nghbrId)) {
                proceed = 1;
              }
            }
          } else {
            for(MInt c = 0; c < solver().grid().m_maxNoChilds; c++) {
              if(!childCode[dir][c]) {
                continue;
              }
              MInt childId = solver().c_childId(nghbrId, c);
              if(childId < 0) {
                continue;
              }
              MFloat dx = F1B2 * (solver().c_cellLengthAtCell(cellId) + solver().c_cellLengthAtCell(childId));
              distance(childId) = mMin(distance(childId), distance(cellId) + dx);
              if(accessed(childId) == 0) {
                currentLayer.push_back(childId);
                accessed(childId) = 1;
                onCurrentLayer(childId) = 1;
                if(distance(childId) < outerBandWidth[mMax(solver().minLevel(), solver().a_level(childId) - 1)]
                                           + F2 * solver().c_cellLengthAtCell(childId)) {
                  proceed = 1;
                }
              }
            }
          }
        } else {
          MInt parentId = solver().c_parentId(cellId);
          while(parentId > -1 && !solver().a_hasNeighbor(parentId, dir)) {
            parentId = solver().c_parentId(parentId);
          }
          if(parentId > -1 && solver().a_hasNeighbor(parentId, dir)) {
            MInt nghbrId = solver().c_neighborId(parentId, dir);
            if(!solver().c_isLeafCell(nghbrId)) {
              continue;
            }
            MFloat dx = F1B2 * (solver().c_cellLengthAtCell(cellId) + solver().c_cellLengthAtCell(nghbrId));
            distance(nghbrId) = mMin(distance(nghbrId), distance(cellId) + dx);
            if(accessed(nghbrId) == 0) {
              currentLayer.push_back(nghbrId);
              accessed(nghbrId) = 1;
              onCurrentLayer(nghbrId) = 1;
              if(distance(nghbrId) < outerBandWidth[mMax(solver().minLevel(), solver().a_level(nghbrId) - 1)]
                                         + F2 * solver().c_cellLengthAtCell(nghbrId)) {
                proceed = 1;
              }
            }
          }
        }
      }
    }
    if(currentLayer.empty()) {
      proceed = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &proceed, 1, MPI_INT, MPI_SUM, solver().mpiComm(), AT_, "MPI_IN_PLACE", "proceed");
  }
}


///  Marks all cells, which are near (within bandWidth) to the given Cells (inList[cellId] = 1) with inList > 0.
///  Useful for refinement (or other purposes!) in all solvers!
///
/// \param[in] bandWidth Number of surrounding cells to be marked
/// \param[in] inList List with origin cells marked with 1
/// \param[in] refineDiagonals Mark diagonal cells
/// \param[in] level Level on which cells are marked
/// \author Tim Wegmann
/// \date 2018-11-14

template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::markSurrndCells(MIntScratchSpace& inList, const MInt bandWidth,
                                                        const MInt level, const MBool refineDiagonals) {
  TRACE();
  static constexpr MInt noDirs = 2 * nDim;
  static constexpr std::array<MInt, 4> diagDirs{3, 2, 0, 1};


  // Loop over the number of refinement-grid-cells and add those to the list:
  for(MInt loopMarker = 1; loopMarker < bandWidth; loopMarker++) {
    for(MInt cellId = 0; cellId < solver().a_noCells(); cellId++) {
      if(solver().grid().tree().level(cellId) != level) {
        continue;
      }
      if(inList(cellId) != loopMarker) {
        continue;
      }

      const MInt cellDist = loopMarker + 1;

      // direct neighbors +x-x+y-y+z-z
      for(MInt mainDir = 0; mainDir < noDirs; mainDir++) {
        const MInt nghbrId = solver().c_neighborId(cellId, mainDir);
        if(nghbrId < 0) {
          continue;
        }
        if(inList(nghbrId) == 0) {
          inList(nghbrId) = cellDist;
        }

        if(refineDiagonals) {
          if(mainDir < 4) {
            // add diagonal neighbors in x-y plane
            const MInt nghbrId2 = solver().c_neighborId(nghbrId, diagDirs.at(mainDir));
            if(nghbrId2 >= 0 && inList(nghbrId2) == 0) {
              inList(nghbrId2) = cellDist;
            }
          }

          IF_CONSTEXPR(nDim == 3) {
            for(MInt dirZ = 4; dirZ < 6; dirZ++) {
              const MInt nghbrIdZ = solver().c_neighborId(cellId, dirZ);
              if(nghbrIdZ < 0) {
                continue;
              }
              if(nghbrIdZ >= 0 && inList(nghbrIdZ) == 0) {
                inList(nghbrIdZ) = cellDist;
              }
              for(MInt dir = 0; dir < 4; dir++) {
                const MInt nghbrId2 = solver().c_neighborId(nghbrIdZ, dir);
                if(nghbrId2 < 0) {
                  continue;
                }
                if(inList(nghbrId2) == 0) {
                  inList(nghbrId2) = cellDist;
                }
                // tridiagonal neighbors
                const MInt triDiagNghbrId = solver().c_neighborId(nghbrId2, diagDirs.at(dir));
                if(triDiagNghbrId >= 0 && inList(triDiagNghbrId) == 0) {
                  inList(triDiagNghbrId) = cellDist;
                }
              }
            }
          }
        }
      }
    }
    exchangeData(inList.data());
    if(grid().azimuthalPeriodicity()) {
      exchangeAzimuthalPer(&inList[0]);
    }
  }
}


/** \brief check that each solver has the required number of haloLayers for leaf cells!!!
 **        TODO labels:toenhance under production, needs to be cleaned up!
 ** @author: Tim Wegmann
 ** @date: 22.05.2020
 **/
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::checkNoHaloLayers() {
#ifndef NDEBUG

  if(Context::propertyExists("periodicCells")) return;
  // NOTE: not working for Alexejs periodicExchange! As periodic cells there are only created
  //      on solver level! The functionality there should be moved into the cartesiangrid!
  //      and once this is done there, the additional exchange Cells are part of the proxy!


  // NOTE: searching only on the same level is not correct if the
  //      domainBoundary is at a lower level and the first haloLayer is refined,
  //      the correct number of leafCell-halo-layers needs to be ensured!

  m_log << "check HaloLayer-Count for solver " << solver().solverId() << " ... ";

  MBool uncorrectLayerCount = false;
  static constexpr MInt cornerIndices[8][3] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
                                               {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};

  /*
  std::multimap< MInt, MInt > firstHaloLayerOld;
  firstHaloLayerOld.clear();

  //loop over all haloCells and find cells which have a window cell as neighbor
  for(MInt cellId =  solver().noInternalCells(); cellId <  solver().c_noCells(); cellId++) {
    ASSERT( solver().a_isHalo(cellId), "");
    //only check for solver leaf cells
    if(!grid().tree().isLeafCell(cellId)) continue;
    for(MInt dir = 0; dir < 2*nDim; dir++) {
      MInt nghbrId = -1;
      if(grid().tree().hasNeighbor(cellId, dir)) {
        nghbrId = grid().tree().neighbor(cellId,dir);
      } else if ( grid().tree().hasParent(cellId) && grid().tree().parent(cellId) < grid().tree().size() ) {
        nghbrId = grid().tree().neighbor(grid().tree().parent(cellId),dir);
      }

      if ( nghbrId < 0 ) continue;


      //if(!grid().tree().isLeafCell(nghbrId) ){
        //get the adjacent Childs
      //  for ( MInt child = 0; child < ipow(2, nDim); child++) {
      //    if ( !childCode[ dir ][ child ] ) continue;
      //    if ( grid().tree().child( nghbrId ,  child ) > -1 ) {
      //      const MInt childId = grid().tree().child( nghbrId ,  child );
      //      if(grid().tree().hasProperty(childId, Cell::IsWindow)){
      //        firstHaloLayerOld.insert( std::make_pair( cellId, m_revDir[dir]));
      //      }
      //    }
      //  }
      //} else
      if(grid().tree().hasProperty(nghbrId, Cell::IsWindow)){
        firstHaloLayerOld.insert(std::make_pair(cellId, m_revDir[dir]));
        ASSERT(cellId > -1, "");
        ASSERT(m_revDir[dir] > -1 && m_revDir[dir] < 2*nDim, std::to_string(m_revDir[dir]));
      }
    }
  }

  for ( std::multimap<MInt,MInt>::iterator it = firstHaloLayerOld.begin();
       it != firstHaloLayerOld.end(); it++ ) {
    MInt cellId = it->first;
    const MInt dir = it->second;
    //const MInt level = grid().tree().level(cellId);
    ASSERT(dir > -1 && dir < 2*nDim, std::to_string(dir));

    if(!grid().tree().isLeafCell(cellId) ) continue;

    for(MInt layer = 1; layer < grid().noHaloLayers(); layer++){
      MInt nghbrId = -1;

      if(grid().tree().hasNeighbor(cellId, dir)) {
        nghbrId = grid().tree().neighbor(cellId,dir);
      } else if (grid().tree().hasParent(cellId) && grid().tree().parent(cellId) < grid().tree().size()) {
        nghbrId = grid().tree().neighbor(grid().tree().parent(cellId),dir);
      }

      if(nghbrId > -1 && !grid().tree().isLeafCell(nghbrId) ){
        //check that the two childs are present!
        for ( MInt child = 0; child < ipow(2, nDim); child++) {
          if ( !childCode[ dir ][ child ] ) continue;
          if ( grid().tree().child( nghbrId ,  child ) < 0 ) {
            std::cerr << "Incorrect halo-Layer count: 1 " << cellId << "  " << solver().domainId() << " "
             << grid().tree().globalId(cellId) << " " <<  solver().solverId() << " " << dir << std::endl;
             uncorrectLayerCount = true;
            break;
          } else {
            //continue the leayer check from one of the childs onward!
            cellId = grid().tree().child( nghbrId ,  child );
          }
        }
      } else if(nghbrId < 0 || nghbrId > grid().tree().size()) {
        //check if the cell is at the domain-bndry
        //meaning that the potential neighbor would have to be outside the geometry!

        //prepare to call pointIsInside in the geometry class:
        MBool outside = true;
        const MFloat cellHalfLength =  F1B2*grid().cellLengthAtCell(cellId);
        MFloat corner[ 3 ] = { 0,0,0 };

        //computing the coordinates the neighbor should be having!
        MFloat coords[nDim];
        for ( MInt k = 0; k < nDim; k++ ) {
          coords[k] = grid().tree().coordinate( cellId, k );
        }
        coords[dir/2] += ((dir%2==0)?-F1:F1) * grid().cellLengthAtCell(cellId);

        for(MInt node = 0; node < ipow(2, nDim); node++){
          for( MInt dim = 0; dim < nDim; dim++ ){
            corner[ dim ] = coords[dim] + cornerIndices[ node ][dim] * cellHalfLength;
          }
          IF_CONSTEXPR(nDim == 2) {
            if(! solver().geometry().pointIsInside( corner )) outside = false;
          } else {
            if ( ! solver().geometry().pointIsInside2( corner )) outside = false;
          }
        }

        if(outside) {
          break;
        } else {

          //check if the cell is at the cutOff
          if(grid().hasCutOff() && grid().tree().cutOff(cellId) > -1) {
            break;
          }


          std::cerr << "Incorrect halo-Layer count: 2 " << cellId << "  " << solver().domainId() << " "
               << grid().tree().globalId(cellId) << " " <<  solver().solverId() << " " << dir << " "
               << layer << std::endl;
               uncorrectLayerCount = true;
          std::cerr << "Coord: " << coords[0] << " " << coords[1] << " " << coords[nDim-1] << std::endl;
          break;
        }
      } else {
        //check passed, continue with next layer
        cellId = nghbrId;
      }




    }
  }

  if(uncorrectLayerCount) {
    mTerm(1,AT_, "Incorrect number of halo-layers!");
  }
*/

  std::set<MInt> haloCells;

  std::multimap<MInt, MInt> firstHaloLayer;
  firstHaloLayer.clear();

  // Create a list of all halo cells to be checked.
  // Add regular halo cells
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noHaloCells(i); j++) {
      const MInt cellId = grid().haloCell(i, j);
      haloCells.insert(cellId);
    }
  }
  // Add azimuthal periodic halo cells
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      const MInt cellId = grid().azimuthalHaloCell(i, j);
      haloCells.insert(cellId);
    }
  }
  for(MInt i = 0; i < grid().noAzimuthalUnmappedHaloCells(); i++) {
    const MInt cellId = grid().azimuthalUnmappedHaloCell(i);
    haloCells.insert(cellId);
  }

  // loop over all haloCells and find cells which have a window cell as neighbor
  // these cells are the first layer of haloCells.
  for(auto it = haloCells.begin(); it != haloCells.end(); ++it) {
    MInt cellId = *it;
    if(!grid().tree().isLeafCell(cellId)) continue;
    ASSERT(solver().a_isHalo(cellId), "");

    // skip partition levell ancestors!
    if(grid().tree().hasProperty(cellId, Cell::IsPartLvlAncestor)) continue;

    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      MInt nghbrId = -1; // adjacent neighbor
      // MInt nghbrId1 = -1; //adjacent child1
      // MInt nghbrId2 = -1; //adjacent child2
      if(grid().tree().hasNeighbor(cellId, dir)) {
        // direct neighbor
        nghbrId = grid().tree().neighbor(cellId, dir);
        // if the direct neighbor has children, get the two closest childs
        /*
        if(!grid().tree().isLeafCell(nghbrId) ){
          for ( MInt child = 0; child < ipow(2, nDim); child++) {
          if( !childCode[ dir ][ child ] ) continue;
          if( grid().tree().child( nghbrId ,  child ) < 0 ) continue;
          if(nghbrId1 < 0) {
          nghbrId1 = grid().tree().child( nghbrId ,  child );
        } else {
          ASSERT(nghbrId2 < 0, "");
          nghbrId2 = grid().tree().child( nghbrId ,  child );
        }
          }
        }
        */
      } else if(grid().tree().hasParent(cellId) && grid().tree().parent(cellId) < grid().tree().size()) {
        nghbrId = grid().tree().neighbor(grid().tree().parent(cellId), dir);
      }
      if(nghbrId < 0) continue;

      // use childs
      // if(nghbrId1 > 0 && grid().tree().hasProperty(nghbrId1, Cell::IsWindow) &&
      //   grid().tree().isLeafCell(cellId)){
      //  firstHaloLayer.insert( std::make_pair( std::make_pair(cellId, i), m_revDir[dir]));
      //}
      // if(nghbrId2 > 0 && grid().tree().hasProperty(nghbrId2, Cell::IsWindow) &&
      //   grid().tree().isLeafCell(cellId)){
      //  firstHaloLayer.insert( std::make_pair( std::make_pair(cellId, i), m_revDir[dir]));
      //}
      // if(/*nghbrId1 < 0 && */solver().a_isWindow(nghbrId) /*&&
      //  grid().tree().isLeafCell(cellId)*/){
      if(!grid().tree().hasProperty(nghbrId, Cell::IsHalo)) {
        firstHaloLayer.insert(std::make_pair(cellId, m_revDir[dir]));
      }
    }
  }


  for(std::multimap<MInt, MInt>::iterator it = firstHaloLayer.begin(); it != firstHaloLayer.end(); it++) {
    MInt cellId = it->first;
    // const MInt nghbrDomain = it->first.second;
    const MInt dir = it->second;
    ASSERT(dir > -1 && dir < 2 * nDim, std::to_string(dir));
    ASSERT(grid().tree().isLeafCell(cellId), "");
    ASSERT(solver().a_isHalo(cellId), "");

    for(MInt layer = 1; layer < grid().noHaloLayers(); layer++) {
      MInt nghbrId = -1;  // adjacent neighbor
      MInt nghbrId1 = -1; // adjacent child1
      // MInt nghbrId2 = -1; //adjacent child2
      if(grid().tree().hasNeighbor(cellId, dir)) {
        nghbrId = grid().tree().neighbor(cellId, dir);
        if(!grid().tree().isLeafCell(nghbrId)) {
          // check the two childs!
          for(MInt child = 0; child < ipow(2, nDim); child++) {
            if(!childCode[dir][child]) continue;
            if(grid().tree().child(nghbrId, child) < 0) {
              if(grid().noHaloLayers() % 2 == 0) {
                // check if the neighboring child could be outside:
                MBool outside = true;
                const MFloat cellHalfLength = F1B4 * grid().cellLengthAtCell(cellId);
                MFloat corner[3] = {0, 0, 0};

                // computing the coordinates the neighbor should be having!
                MFloat coords[nDim];
                for(MInt k = 0; k < nDim; k++) {
                  coords[k] = grid().tree().coordinate(nghbrId, k);
                }
                for(MInt k = 0; k < nDim; k++) {
                  coords[k] += cornerIndices[child][k] * cellHalfLength;
                }

                for(MInt node = 0; node < ipow(2, nDim); node++) {
                  for(MInt dim = 0; dim < nDim; dim++) {
                    corner[dim] = coords[dim] + cornerIndices[node][dim] * F1B2 * cellHalfLength;
                  }
                  IF_CONSTEXPR(nDim == 2) {
                    if(!solver().geometry().pointIsInside(corner)) outside = false;
                  }
                  else {
                    if(!solver().geometry().pointIsInside2(corner)) outside = false;
                  }
                }

                if(outside) continue;

                uncorrectLayerCount = true;
                std::cerr << "Incorrect halo-Layer count: 1 " << cellId << "  " << solver().domainId() << " "
                          << grid().tree().globalId(cellId) << " " << solver().solverId() << " " << dir << std::endl;
                std::cerr << grid().tree().coordinate(cellId, 0) << " " << grid().tree().coordinate(cellId, 1) << " "
                          << grid().tree().coordinate(cellId, nDim - 1) << std::endl;
                break;
              }
            } else {
              if(nghbrId1 < 0) {
                nghbrId1 = grid().tree().child(nghbrId, child);
              } /*else {
                ASSERT(nghbrId2 < 0, "");
                nghbrId2 = grid().tree().child( nghbrId ,  child );
              }
              */
            }
          }
        }
      } else if(grid().tree().hasParent(cellId) && grid().tree().parent(cellId) < grid().tree().size()) {
        nghbrId = grid().tree().neighbor(grid().tree().parent(cellId), dir);
      }

      // check if the cell is at the domain-bndry
      // meaning that the potential neighbor would have to be outside the geometry!
      if(nghbrId < 0 || nghbrId > grid().tree().size()) {
        // prepare to call pointIsInside in the geometry class:
        MBool outside = true;
        const MFloat cellHalfLength = F1B2 * grid().cellLengthAtCell(cellId);
        MFloat corner[3] = {0, 0, 0};

        // computing the coordinates the neighbor should be having!
        MFloat coords[nDim];
        for(MInt k = 0; k < nDim; k++) {
          coords[k] = grid().tree().coordinate(cellId, k);
        }
        coords[dir / 2] += ((dir % 2 == 0) ? -F1 : F1) * grid().cellLengthAtCell(cellId);

        for(MInt node = 0; node < ipow(2, nDim); node++) {
          for(MInt dim = 0; dim < nDim; dim++) {
            corner[dim] = coords[dim] + cornerIndices[node][dim] * cellHalfLength;
          }
          IF_CONSTEXPR(nDim == 2) {
            if(!solver().geometry().pointIsInside(corner)) outside = false;
          }
          else {
            if(!solver().geometry().pointIsInside2(corner)) outside = false;
          }
        }

        if(outside) {
          break;
        } else {
          // check if the cell is at the cutOff
          if(grid().hasCutOff() && grid().tree().cutOff(cellId) > -1) {
            break;
          }
        }
      }

      // if the neighbor is not present at all
      if(nghbrId < 0) {
        uncorrectLayerCount = true;
        std::cerr << "Incorrect halo-Layer count: 2 " << solver().domainId() << " " << solver().solverId() << " "
                  << layer << " " << cellId << "/" << grid().tree().globalId(cellId) << " " << dir << " "
                  << grid().tree().coordinate(cellId, 0) << " " << grid().tree().coordinate(cellId, 1) << " "
                  << grid().tree().coordinate(cellId, nDim - 1) << " " << grid().isPeriodic(cellId) << std::endl;
        break;
      }

      /*
      //check if the corresponding neighbor is a window cell,
      //meaning that the direction changed due to a domain-corner!
      if(nghbrId1 < 0  ) {
        if(!solver().a_isHalo(nghbrId)) {
          ASSERT(grid().tree().hasProperty(nghbrId, Cell::IsWindow), "");
          break;
        }
      } else if ( nghbrId1 > -1 ) {
        if(!solver().a_isHalo(nghbrId1)) {
          //neighbor appears to be a window cell again, meaning, that we went back into the domain.
          ASSERT(grid().tree().hasProperty(nghbrId1, Cell::IsWindow), "");
          break;
        }
      }
      if(nghbrId2 > -1 ) {
        //check ngbrId1
        if(!solver().a_isHalo(nghbrId2)) {
          ASSERT(grid().tree().hasProperty(nghbrId2, Cell::IsWindow), "");
          break;
        }
      }

      // check that the ngbrId is also a haloCell and is a halo for the correct domain!
      MBool cellFound = false;
      MBool cellFound2 = false;
      if(nghbrId1 < 0  ) {
        //check ngbrId
        for(MInt j = 0; j < grid().noHaloCells(nghbrDomain); j++){
          const MInt haloCell = grid().haloCell(nghbrDomain,j);
          if(haloCell == nghbrId) {
            cellFound = true;
            break;
          }
        }
      } else if ( nghbrId1 > -1 ) {
        //check ngbrId1
        for(MInt j = 0; j < grid().noHaloCells(nghbrDomain); j++){
          const MInt haloCell = grid().haloCell(nghbrDomain,j);
          if(haloCell == nghbrId1) {
            cellFound = true;
            break;
          }
        }
      }
      if(nghbrId2 > -1 ) {
        //check ngbrId1
        for(MInt j = 0; j < grid().noHaloCells(nghbrDomain); j++){
          const MInt haloCell = grid().haloCell(nghbrDomain,j);
          if(haloCell == nghbrId2) {
            cellFound2 = true;
            break;
          }
        }
      }

      if(!cellFound) {
        std::cerr << "Incorrect halo-Layer count: 3 " << cellId << "  " << solver().domainId() << " "
                  << grid().tree().globalId(cellId) << " " <<  solver().solverId() << " " << dir << " "
                  << layer << " " << nghbrId << " " << solver().a_isHalo(nghbrId) << std::endl;
        uncorrectLayerCount = true;
        break;
      }

      if(nghbrId2 > -1 && !cellFound2) {
        std::cerr << "Incorrect halo-Layer count: 4 " << cellId << "  " << solver().domainId() << " "
                  << grid().tree().globalId(cellId) << " " <<  solver().solverId() << " " << dir << " "
                  << layer << " " << nghbrId << std::endl;
        uncorrectLayerCount = true;
        break;
      }
      */

      // continue towards the next layer
      if(nghbrId1 < 0) {
        cellId = nghbrId;
      } else {
        // two remaining childs to iterate forward
        // for now just the first child is used...
        cellId = nghbrId1;
      }
    }
  }


  if(uncorrectLayerCount) {
    mTerm(1, AT_, "Incorrect number of halo-layers!");
  }


  m_log << "done" << std::endl;

#endif
}

// Check if point is located in cell
template <MInt nDim, class SolverType>
MInt CartesianSolver<nDim, SolverType>::inCell(const MInt cellId, MFloat* point, MFloat fac) {
  TRACE();

  MFloat xmin, xmax;
  MFloat ymin, ymax;
  MFloat zmin, zmax;
  MFloat cellHalfLength = F1B2 * grid().cellLengthAtCell(cellId);

  xmin = grid().tree().coordinate(cellId, 0) - (cellHalfLength * fac);
  xmax = grid().tree().coordinate(cellId, 0) + (cellHalfLength * fac);
  ymin = grid().tree().coordinate(cellId, 1) - (cellHalfLength * fac);
  ymax = grid().tree().coordinate(cellId, 1) + (cellHalfLength * fac);

  if(nDim == 2) {
    if(point[0] < xmax && point[0] >= xmin && point[1] < ymax && point[1] >= ymin) {
      return true;
    } else {
      return false;
    }
  } else {
    zmin = grid().tree().coordinate(cellId, 2) - (cellHalfLength * fac);
    zmax = grid().tree().coordinate(cellId, 2) + (cellHalfLength * fac);

    if(point[0] < xmax && point[0] >= xmin && point[1] < ymax && point[1] >= ymin && point[2] < zmax
       && point[2] >= zmin) {
      return true;
    } else {
      return false;
    }
  }
}


///  works together with interpolateFieldData()
///  cellId: cell on the highes refinement level which contains the point
///  interpolationCells: collection of stencil cellIds
///  considering the cartesian cell stencil
/// \author Claudia Guenther, Tim Wegmann
/// \date 03/2012, 2020-04-01

template <MInt nDim, class SolverType>
MInt CartesianSolver<nDim, SolverType>::setUpInterpolationStencil(const MInt cellId,
                                                                  MInt* interpolationCells,
                                                                  const MFloat* point,
                                                                  std::function<MBool(MInt, MInt)>
                                                                      neighborCheck,
                                                                  MBool allowIncompleteStencil) {
  TRACE();

  MInt position = 0;
  const MInt maxPosition = IPOW2(nDim) - 1;
  MBool invalidStencil = false;

  //------------------------------

  if(point[0] > grid().tree().coordinate(cellId, 0)) {
    position += 1;
  }
  if(point[1] > grid().tree().coordinate(cellId, 1)) {
    position += 2;
  }

  IF_CONSTEXPR(nDim == 3) {
    if(point[2] > grid().tree().coordinate(cellId, 2)) {
      position += 4;
    }
  }

  position = maxPosition - position;

  // check if cell is boundary cell and located at position at the boundary in the stencil. if yes, take other stencil.

  MInt xNghbrDir = (position + 1) % 2;
  MInt posIncrementX = -1 + 2 * xNghbrDir;
  MInt yNghbrDir = (position / 2 + 1) % 2;
  MInt posIncrementY = -2 + 4 * yNghbrDir;
  yNghbrDir += 2;
  MInt zNghbrDir = -1;
  MInt posIncrementZ = -1;


  if(!neighborCheck(cellId, xNghbrDir) || !grid().tree().hasNeighbor(cellId, xNghbrDir)) {
    position += posIncrementX;
  }
  if(!neighborCheck(cellId, yNghbrDir) || !grid().tree().hasNeighbor(cellId, yNghbrDir)) {
    position += posIncrementY;
  }

  IF_CONSTEXPR(nDim == 3) {
    zNghbrDir = (position / 4 + 1) % 2;
    posIncrementZ = -4 + 8 * zNghbrDir;
    zNghbrDir += 4;

    if(!neighborCheck(cellId, zNghbrDir) || !grid().tree().hasNeighbor(cellId, zNghbrDir)) {
      position += posIncrementZ;
    }
  }

  interpolationCells[position] = cellId;

  IF_CONSTEXPR(nDim == 2) {
    MInt nghbrX = -1;
    MInt nghbrY = -1;
    if(position % 2 == 0) {
      xNghbrDir = 1;
      nghbrX = neighborCheck(cellId, xNghbrDir) ? grid().tree().neighbor(cellId, xNghbrDir) : -1;
      posIncrementX = 1;
    } else {
      xNghbrDir = 0;
      nghbrX = neighborCheck(cellId, xNghbrDir) ? grid().tree().neighbor(cellId, xNghbrDir) : -1;
      posIncrementX = -1;
    }
    if(((MInt)(position / 2)) % 2 == 0) {
      yNghbrDir = 3;
      nghbrY = neighborCheck(cellId, yNghbrDir) ? grid().tree().neighbor(cellId, yNghbrDir) : -1;
      posIncrementY = 2;
    } else {
      yNghbrDir = 2;
      nghbrY = neighborCheck(cellId, yNghbrDir) ? grid().tree().neighbor(cellId, yNghbrDir) : -1;
      posIncrementY = -2;
    }
    interpolationCells[position + posIncrementX] = nghbrX;
    interpolationCells[position + posIncrementY] = nghbrY;
    if(nghbrX != -1)
      interpolationCells[position + posIncrementX + posIncrementY] =
          neighborCheck(nghbrX, yNghbrDir) ? grid().tree().neighbor(nghbrX, yNghbrDir) : -1;
    else if(nghbrY != -1)
      interpolationCells[position + posIncrementX + posIncrementY] =
          neighborCheck(nghbrY, xNghbrDir) ? grid().tree().neighbor(nghbrY, xNghbrDir) : -1;
    else
      interpolationCells[position + posIncrementX + posIncrementY] = -1;
  }
  else {
    MInt nghbrX = -1;
    MInt nghbrY = -1;
    MInt nghbrZ = -1;
    if(position % 2 == 0) {
      xNghbrDir = 1;
      nghbrX = neighborCheck(cellId, xNghbrDir) ? grid().tree().neighbor(cellId, xNghbrDir) : -1;
      posIncrementX = 1;
    } else {
      xNghbrDir = 0;
      nghbrX = neighborCheck(cellId, xNghbrDir) ? grid().tree().neighbor(cellId, xNghbrDir) : -1;
      posIncrementX = -1;
    }
    if(((MInt)(position / 2)) % 2 == 0) {
      yNghbrDir = 3;
      nghbrY = neighborCheck(cellId, yNghbrDir) ? grid().tree().neighbor(cellId, yNghbrDir) : -1;
      posIncrementY = 2;
    } else {
      yNghbrDir = 2;
      nghbrY = neighborCheck(cellId, yNghbrDir) ? grid().tree().neighbor(cellId, yNghbrDir) : -1;
      posIncrementY = -2;
    }
    if((MInt)(position / 4) == 0) {
      zNghbrDir = 5;
      nghbrZ = neighborCheck(cellId, zNghbrDir) ? grid().tree().neighbor(cellId, zNghbrDir) : -1;
      posIncrementZ = 4;
    } else {
      zNghbrDir = 4;
      nghbrZ = neighborCheck(cellId, zNghbrDir) ? grid().tree().neighbor(cellId, zNghbrDir) : -1;
      posIncrementZ = -4;
    }
    interpolationCells[position + posIncrementX] = nghbrX;
    interpolationCells[position + posIncrementY] = nghbrY;
    interpolationCells[position + posIncrementZ] = nghbrZ;
    if(nghbrX > -1) {
      interpolationCells[position + posIncrementX + posIncrementZ] =
          neighborCheck(nghbrX, zNghbrDir) ? grid().tree().neighbor(nghbrX, zNghbrDir) : -1;
      interpolationCells[position + posIncrementX + posIncrementY] =
          neighborCheck(nghbrX, yNghbrDir) ? grid().tree().neighbor(nghbrX, yNghbrDir) : -1;
    } else {
      interpolationCells[position + posIncrementX + posIncrementZ] = -1;
      interpolationCells[position + posIncrementX + posIncrementY] = -1;
    }
    if(nghbrY != -1) {
      interpolationCells[position + posIncrementY + posIncrementZ] =
          neighborCheck(nghbrY, zNghbrDir) ? grid().tree().neighbor(nghbrY, zNghbrDir) : -1;
    } else {
      interpolationCells[position + posIncrementY + posIncrementZ] = -1;
    }
    const MInt nghbrXY = interpolationCells[position + posIncrementX + posIncrementY];
    if(nghbrX > -1 && nghbrXY > -1)
      interpolationCells[position + posIncrementX + posIncrementY + posIncrementZ] =
          neighborCheck(nghbrXY, zNghbrDir) ? grid().tree().neighbor(nghbrXY, zNghbrDir) : -1;
    else
      interpolationCells[position + posIncrementX + posIncrementY + posIncrementZ] = -1;
  }

  for(MInt p = 0; p <= maxPosition; p++) {
    if(interpolationCells[p] == -1) {
      invalidStencil = true;
      break;
    }
  }
  if(invalidStencil && !allowIncompleteStencil) {
    for(MInt p = 0; p <= maxPosition; p++) {
      interpolationCells[p] = cellId;
    }
    position = -1;
  }

  return position;
}


///  bilinear interpolation of cell-Centered scalar(float) field data onto a higher refined mesh
///  works together with setUpInterpolationStencil()
//   mode: cartesian interpolation = true:
//         -> cartesian coordinates (c_coordinate/a_coordinate of non-boundary-cells!)
//            based on pure cartesian cooordinates (i.e. known reconstruction constants!)
//            NOTE:  requires ordering of the cells in the cartesin-directions
//         cartesian interpolation = false:
/// \author Claudia Guenther, Tim Wegmann
/// \date 03/2012, 2020-04-01
template <MInt nDim, class SolverType>
template <MBool cartesianInterpolation>
MFloat CartesianSolver<nDim, SolverType>::interpolateFieldData(MInt* interpolationCells, MFloat* point, MInt varId,
                                                               std::function<MFloat(MInt, MInt)> scalarField,
                                                               std::function<MFloat(MInt, MInt)> coordinate) {
  TRACE();

  IF_CONSTEXPR(cartesianInterpolation) {
    const MFloat xMin = coordinate(interpolationCells[0], 0);
    const MFloat xMax = coordinate(interpolationCells[1], 0);
    const MFloat yMin = coordinate(interpolationCells[0], 1);
    const MFloat yMax = coordinate(interpolationCells[2], 1);
    const MFloat deltaX_Minus = point[0] - xMin;
    const MFloat deltaX_Plus = xMax - point[0];
    const MFloat deltaY_Minus = point[1] - yMin;
    const MFloat deltaY_Plus = yMax - point[1];
    const MFloat Delta = deltaX_Minus + deltaX_Plus;

    IF_CONSTEXPR(nDim == 2) {
      const MFloat phi = (scalarField(interpolationCells[0], varId) * deltaX_Plus * deltaY_Plus
                          + scalarField(interpolationCells[1], varId) * deltaX_Minus * deltaY_Plus
                          + scalarField(interpolationCells[2], varId) * deltaX_Plus * deltaY_Minus
                          + scalarField(interpolationCells[3], varId) * deltaX_Minus * deltaY_Minus)
                         / (Delta * Delta);
      return phi;
    }
    else { // nDim == 3
      const MFloat zMin = coordinate(interpolationCells[0], 2);
      const MFloat zMax = coordinate(interpolationCells[4], 2);
      const MFloat deltaZ_Minus = point[2] - zMin;
      const MFloat deltaZ_Plus = zMax - point[2];
      const MFloat phi = (scalarField(interpolationCells[0], varId) * deltaX_Plus * deltaY_Plus * deltaZ_Plus
                          + scalarField(interpolationCells[1], varId) * deltaX_Minus * deltaY_Plus * deltaZ_Plus
                          + scalarField(interpolationCells[2], varId) * deltaX_Plus * deltaY_Minus * deltaZ_Plus
                          + scalarField(interpolationCells[3], varId) * deltaX_Minus * deltaY_Minus * deltaZ_Plus
                          + scalarField(interpolationCells[4], varId) * deltaX_Plus * deltaY_Plus * deltaZ_Minus
                          + scalarField(interpolationCells[5], varId) * deltaX_Minus * deltaY_Plus * deltaZ_Minus
                          + scalarField(interpolationCells[6], varId) * deltaX_Plus * deltaY_Minus * deltaZ_Minus
                          + scalarField(interpolationCells[7], varId) * deltaX_Minus * deltaY_Minus * deltaZ_Minus)
                         / (Delta * Delta * Delta);
      return phi;
    }
  }
  else {
    MFloat phi = 0;
    MFloat sumDist = 0;
    MFloat dist = -1;
    for(MInt i = 0; i < IPOW2(nDim); i++) {
      if(interpolationCells[i] < 0) continue;
      IF_CONSTEXPR(nDim == 2) {
        dist = std::sqrt(POW2(point[0] - coordinate(interpolationCells[i], 0))
                         + POW2(point[1] - coordinate(interpolationCells[i], 1)));
      }
      else {
        dist = std::sqrt(POW2(point[0] - coordinate(interpolationCells[i], 0))
                         + POW2(point[1] - coordinate(interpolationCells[i], 1))
                         + POW2(point[2] - coordinate(interpolationCells[i], 2)));
      }
      sumDist += dist;
      phi += scalarField(interpolationCells[i], varId) * dist;
    }
    return phi / sumDist;
  }
}

/// least squares interpolation of cell-Centered scalar(float) field data
/// to the given coordinate
/// works together with setUpInterpolationStencil()
/// \author Tim Wegmann
/// \date 2021-01-04
template <MInt nDim, class SolverType>
MFloat CartesianSolver<nDim, SolverType>::leastSquaresInterpolation(MInt* interpolationCells, MFloat* point, MInt varId,
                                                                    std::function<MFloat(MInt, MInt)> scalarField,
                                                                    std::function<MFloat(MInt, MInt)> coordinate) {
  std::array<MFloat, nDim + 1> b;
  std::array<std::array<MFloat, (nDim + 1)>, (nDim + 1)> A;
  std::array<MFloat, nDim + 1> adjA;
  b.fill(0.0);
  for(MInt i = 0; i < nDim + 1; i++) {
    for(MInt j = 0; j < nDim + 1; j++) {
      A[i][j] = 0.0;
    }
  }
  for(MInt n = 0; n < IPOW2(nDim); n++) {
    const MInt cellId = interpolationCells[n];
    b[0] += scalarField(cellId, varId);
    MFloat deltax = point[0] - coordinate(cellId, 0);
    MFloat deltay = point[1] - coordinate(cellId, 1);
    b[1] += scalarField(cellId, varId) * deltax;
    b[2] += scalarField(cellId, varId) * deltay;

    A[0][0]++;
    A[1][1] += POW2(deltax);
    A[2][2] += POW2(deltay);
    A[1][0] += deltax;
    A[2][0] += deltay;
    A[2][1] += deltax * deltay;
    IF_CONSTEXPR(nDim == 3) {
      MFloat deltaz = point[2] - coordinate(cellId, 2);
      b[3] += scalarField(cellId, varId) * deltaz;
      A[3][0] += deltaz;
      A[3][1] += deltaz * deltax;
      A[3][2] += deltaz * deltay;
      A[3][3] += POW2(deltaz);
    }
  }

  A[0][1] = A[1][0];
  A[0][2] = A[2][0];
  A[1][2] = A[2][1];
  IF_CONSTEXPR(nDim == 3) {
    A[0][3] = A[3][0];
    A[1][3] = A[3][1];
    A[2][3] = A[3][2];
  }

  MFloat det = maia::math::determinant(A);
  ASSERT(fabs(det) > MFloatEps, "Poor least-squares stencil!");

  maia::math::adjoint1stRow(A, adjA);

  MFloat sum = 0;
  for(MInt i = 0; i < nDim + 1; i++) {
    sum += adjA[i] * b[i];
  }

  return sum / det;
}


///  maps the cells in the solver to the reference grid file, for Interpolation!
/// \author Tim Wegmann
/// \date 2018-11-14

template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::mapInterpolationCells(std::map<MInt, MInt>& cellMap) {
  TRACE();

  ParallelIo srcGrid("src/grid.Netcdf", maia::parallel_io::PIO_READ, solver().mpiComm());

  MInt noMinCells = solver().m_localMinCellsOffsets[1] - solver().m_localMinCellsOffsets[0] + 1; // see cartesiangrid
  srcGrid.setOffset(noMinCells, solver().m_localMinCellsOffsets[0]);

  // memory
  MIntScratchSpace srcMinCellsGlobalIds(noMinCells, AT_, "srcMinCellsGlobalIds");
  MIntScratchSpace srcMinCellsNoOffsprings(noMinCells, AT_, "srcMinCellsNoOffsprings");
  // read
  srcGrid.readArray(srcMinCellsGlobalIds.getPointer(), "minCellsId");
  srcGrid.readArray(srcMinCellsNoOffsprings.getPointer(), "minCellsNoOffsprings");

  // sum offsprings
  MInt srcFirstMinCellGlobalId = srcMinCellsGlobalIds(0);
  MInt srcNoCells = 0;
  for(MInt i = 0; i < noMinCells; i++) {
    srcNoCells += srcMinCellsNoOffsprings(i);
    srcMinCellsGlobalIds(i) -= srcFirstMinCellGlobalId;
  }
  // new offsets
  srcGrid.setOffset(srcNoCells, srcFirstMinCellGlobalId);


  // more memory
  //   grid
  MFloatScratchSpace srcCoords(srcNoCells, nDim, AT_, "srcCoords");
  MIntScratchSpace srcNoChildIds(srcNoCells, AT_, "srcNoChildIds");
  MIntScratchSpace srcLevel(srcNoCells, AT_, "srcLevel");
  MIntScratchSpace srcChildIds(srcNoCells, IPOW2(nDim), AT_, "srcChildIds");

  // load coords
  std::vector<MString> varNames;
  for(MInt j = 0; j < nDim; ++j) {
    std::stringstream ss;
    ss << "coordinates_" << j;
    varNames.push_back(ss.str());
    ss.clear();
    ss.str("");
  }
  srcGrid.readArray(&srcCoords(0, 0), varNames, nDim);
  varNames.clear();
  // load nmbr child
  srcGrid.readArray(srcNoChildIds.begin(), "noChildIds");
  // level
  srcGrid.readArray(srcLevel.begin(), "level_0");
  // load childIds
  for(MInt j = 0; j < IPOW2(nDim); ++j) {
    std::stringstream ss;
    ss << "childIds_" << j;
    varNames.push_back(ss.str());
    ss.clear();
    ss.str("");
  }
  srcGrid.readArray(&srcChildIds(0, 0), varNames, IPOW2(nDim));
  varNames.clear();
  // correct Ids
  for(MInt i = 0; i < srcNoCells; i++) {
    for(MInt j = 0; j < IPOW2(nDim); j++) {
      srcChildIds(i, j) -= srcFirstMinCellGlobalId;
    }
  }

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(solver().c_noChildren(cellId)) {
      continue;
    }
    auto* coord = (MFloat*)(&(solver().a_coordinate(cellId, 0)));
    MInt srcId = -1;
    // find min Cell
    for(MInt srcMinCellId = 0; srcMinCellId < noMinCells; srcMinCellId++) {
      MInt srcMinCellGlobalId = srcMinCellsGlobalIds(srcMinCellId);
      MFloat* srcCoord = &srcCoords(srcMinCellGlobalId, 0);
      MFloat hdc = solver().c_cellLengthAtLevel(srcLevel(srcMinCellGlobalId) + 1);
      MInt insideCnt = 0;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(coord[dim] < srcCoord[dim] + hdc && coord[dim] >= srcCoord[dim] - hdc) {
          insideCnt++;
        }
      }

      if(insideCnt == nDim) {
        srcId = srcMinCellGlobalId;
        break;
      }
      if(srcMinCellId == (noMinCells - 1)) {
        mTerm(1, " should not happen, loc cell not in src region");
      }
    }
    // loop up (or down)
    while(srcNoChildIds(srcId) != 0) {
      for(MInt cc = 0; cc < IPOW2(nDim); cc++) {
        MInt tcid = srcChildIds(srcId, cc);
        if(tcid >= 0) {
          MFloat* srcCoord = &srcCoords(tcid, 0);
          MFloat hdc = solver().c_cellLengthAtLevel(srcLevel(tcid) + 1);
          // if the cell is inside the child break loop
          MInt insideCnt = 0;
          for(MInt dim = 0; dim < nDim; dim++) {
            if(coord[dim] < srcCoord[dim] + hdc && coord[dim] >= srcCoord[dim] - hdc) {
              insideCnt++;
            }
          }
          if(insideCnt == nDim) {
            srcId = tcid;
            break;
          }
        }
        if(cc == (IPOW2(nDim) - 1)) {
          mTerm(1, " should not happen, loc cell not in any child");
        }
      }
    }

    cellMap[cellId] = srcId;
  }
}


template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::patchRefinement(std::vector<std::vector<MFloat>>& sensors,
                                                        std::vector<std::bitset<64>>& sensorCellFlag,
                                                        std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                        MInt sen) {
  if(solver().m_testPatch && globalTimeStep < solver().m_patchStartTimeStep
     && globalTimeStep > solver().m_patchStopTimeStep) {
    return;
  }

  m_log << "   - Sensor preparation for the patch sensor" << std::endl;

  MFloat coordinates[nDim]{};
  MFloat patchDimensions[nDim * 2]{};

  for(MInt lvl = 0; lvl < solver().m_patchRefinement->noLocalPatchRfnLvls(); lvl++) {
    const MInt level = lvl + solver().maxUniformRefinementLevel();
    for(MInt patch = 0; patch < solver().m_patchRefinement->noPatchesPerLevel(lvl); patch++) {
      MString patchStr = solver().m_patchRefinement->localRfnLevelMethods[lvl].substr(patch, 1);
      const MInt pos = solver().m_patchRefinement->localRfnLevelPropertiesOffset[lvl] + patch;
      for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
        if(solver().grid().tree().level(cellId) > level) continue;
        for(MInt i = 0; i < nDim; i++) {
          coordinates[i] = solver().grid().tree().coordinate(cellId, i);
          patchDimensions[i] = solver().m_patchRefinement->localRfnPatchProperties[pos][i];
          if(patchStr == "B") {
            patchDimensions[nDim + i] = solver().m_patchRefinement->localRfnPatchProperties[pos][nDim + i];
          }
        }
        MBool keepCell = false;

        if(patchStr == "B") {
          if(maia::geom::isPointInsideBox<MFloat, nDim>(&coordinates[0], &patchDimensions[0])) {
            keepCell = true;
          }
        } else if(patchStr == "R") {
          if(maia::geom::isPointInsideSphere<nDim>(&coordinates[0], &patchDimensions[0],
                                                   solver().m_patchRefinement->localRfnPatchProperties[pos][nDim])) {
            keepCell = true;
          }
        } else if(patchStr == "H") {
          // checks if point lies within a ring segment, requires center, r_min, r_max, phi_min, phi_max, axis

          std::array<MFloat, nDim> center{};
          for(MInt dim = 0; dim < nDim; dim++) {
            center[dim] = solver().m_patchRefinement->localRfnPatchProperties[pos][dim];
          }

          std::array<MFloat, 2> R{};
          R[0] = solver().m_patchRefinement->localRfnPatchProperties[pos][nDim];
          R[1] = solver().m_patchRefinement->localRfnPatchProperties[pos][nDim + 1];

          std::array<MFloat, 2> phi{};
          phi[0] = solver().m_patchRefinement->localRfnPatchProperties[pos][nDim + 2];
          phi[1] = solver().m_patchRefinement->localRfnPatchProperties[pos][nDim + 3];

          constexpr MInt axis = 2; // currently hardcoded

          if(maia::geom::isPointInsideRingSegment<nDim>(&coordinates[0], &center[0], R.data(), phi.data(), axis)) {
            keepCell = true;
          }
        }

        MBool refineCell = false;
        if(keepCell && solver().grid().tree().level(cellId) <= level) {
          refineCell = true;
          keepCell = false;
        }

        if(refineCell) {
          const MInt gridCellId = grid().tree().solver2grid(cellId);
          sensors[sensorOffset + sen][gridCellId] = 1.0;
          sensorCellFlag[gridCellId][sensorOffset + sen] = true;
        } else if(keepCell) {
          const MInt gridCellId = grid().tree().solver2grid(cellId);
          sensors[sensorOffset + sen][gridCellId] = 0.0;
          sensorCellFlag[gridCellId][sensorOffset + sen] = false;
        }
      }
    }
  }

  sensorWeight[sensorOffset + sen] = solver().m_sensorWeight[sen];
}

/// \author Thomas Hoesgen
/// \date 2020-02-19
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::setHaloCellsOnInactiveRanks() {
  TRACE();

  ASSERT(solver().m_cells.size() == 0, "");
  ASSERT(!solver().isActive(), "");

  for(MInt cellId = 0; cellId < solver().grid().tree().size(); cellId++) {
    MInt solverCellId = solver().m_cells.size();
    solver().m_cells.append();
    solver().a_resetPropertiesSolver(solverCellId);
    solver().a_isHalo(solverCellId) = true;

    ASSERT(solver().grid().tree().grid2solver(solver().grid().tree().solver2grid(solverCellId)) == solverCellId, "");
  }
}

///  read Information for patch refinement (the same way as for the grid generation!)
/// \author Tim Wegmann
/// \date 2020-04-17
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::readPatchProperties() {
  m_log << "      + reading patch refinement properties for solver " << solver().solverId() << std::endl;

  mAlloc(m_patchRefinement, 1, "m_patchRefinement", AT_);

  MString localRfnLevelMethods = Context::getSolverProperty<MString>("localRfnLvlMethods", solver().solverId(), AT_);

  // sanity check
  if(localRfnLevelMethods.front() == '-' || localRfnLevelMethods.back() == '-')
    mTerm(1, AT_, "ERROR: localRfnLvlMethods begins or ends with hyphen!");

  if(localRfnLevelMethods.find(" ") != MString::npos) mTerm(1, AT_, "ERROR: localRfnLvlMethods contains space!");

  MString::size_type prev_pos = MString::npos;
  MString::size_type nextMarker = 0;

  // save substrings marked with '-' as the patch-refinement methods
  while(nextMarker != MString::npos) {
    prev_pos++;
    nextMarker = localRfnLevelMethods.find("-", prev_pos);
    MString substring(localRfnLevelMethods.substr(prev_pos, nextMarker - prev_pos));
    m_patchRefinement->localRfnLevelMethods.push_back(substring);
    prev_pos = nextMarker;
  }

  // how many patches do we have
  MInt numPatches = 0;
  for(MInt l = 0; l < m_patchRefinement->noLocalPatchRfnLvls(); l++) {
    numPatches += m_patchRefinement->noPatchesPerLevel(l);
  }

  m_log << "      - " << numPatches << " Total-Patches. Ranging from level " << solver().maxUniformRefinementLevel()
        << " - " << solver().maxUniformRefinementLevel() + m_patchRefinement->noLocalPatchRfnLvls() << std::endl;

  // allocate vector containing the offset to each rfnLvl in the patchProperties matrix
  mAlloc(m_patchRefinement->localRfnLevelPropertiesOffset, m_patchRefinement->noLocalPatchRfnLvls() + 1,
         "localRfnLevelPropertiesOffset", 0, AT_);

  // allocate vector containing the number of properties for each patch
  mAlloc(m_patchRefinement->noLocalRfnPatchProperties, numPatches, "noLocalRfnPatchProperties", 0, AT_);

  // calculate number of properties required for each patchRfnLvl and set offsets
  m_patchRefinement->localRfnLevelPropertiesOffset[0] = 0;
  MInt count_req = 0;
  MInt s = 0;
  for(MInt i = 0; i < m_patchRefinement->noLocalPatchRfnLvls(); i++) {
    const MString lvlStr = m_patchRefinement->localRfnLevelMethods[i];
    for(MString::size_type j = 0; j < lvlStr.size(); j++) {
      const MString patchStr = lvlStr.substr(j, 1);
      if(patchStr == "B") {
        m_patchRefinement->noLocalRfnPatchProperties[s] = 2 * nDim;
      } else if(patchStr == "R") {
        m_patchRefinement->noLocalRfnPatchProperties[s] = nDim + 1;
      } else if(patchStr == "H") {
        m_patchRefinement->noLocalRfnPatchProperties[s] = nDim + 4;
      } else {
        mTerm(1, AT_, "Not yet implemented, please do so...");
      }
      count_req += m_patchRefinement->noLocalRfnPatchProperties[s];
      s++;
    }
    m_patchRefinement->localRfnLevelPropertiesOffset[i + 1] = s;
  }

  // allocate matrix containing the properties for each patch
  mAlloc(m_patchRefinement->localRfnPatchProperties, numPatches, m_patchRefinement->noLocalRfnPatchProperties,
         "localRfnPatchProperties", AT_);

  MInt count_prov = Context::propertyLength("localRfnLevelProperties", solver().solverId());

  // sanity check
  if(count_prov < count_req)
    mTerm(1, AT_, "ERROR: number of localRfnLevelProperties does not match the requested value!");

  // load properties of patch refinement methods
  MInt lvl = 0;
  MInt j = 0;
  for(MInt i = 0; i < count_req; i++) {
    m_patchRefinement->localRfnPatchProperties[lvl][j] =
        Context::getSolverProperty<MFloat>("localRfnLevelProperties", solver().solverId(), AT_, i);
    j++;
    if(j == m_patchRefinement->noLocalRfnPatchProperties[lvl]) {
      j = 0;
      lvl++;
    }
  }

  // properties to test patch build-up and reduction
  m_testPatch = Context::getSolverProperty<MBool>("testPatchRefinement", solver().solverId(), AT_, &m_testPatch);

  if(m_testPatch) {
    m_patchStartTimeStep =
        Context::getSolverProperty<MInt>("patchStart", solver().solverId(), AT_, &m_patchStartTimeStep);
    m_patchStopTimeStep = Context::getSolverProperty<MInt>("patchStop", solver().solverId(), AT_, &m_patchStopTimeStep);
  }
}

template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::initAdaptation() {
  /*! \page propertiesAMR
    \section adaptation
    <code>MBool FvMbSolver::m_adaptation</code>\n
    default = <code>false</code>\n \n
    Enable the use of adaptively refined grids. \n \n
    Possible values are:
    <ul>
      <li><code>0</code> (off)</li>
      <li><code>1</code> (on)</li>
    </ul>
    Keywords: <i>GENERAL, ADAPTATION, GRID</i>
  */
  m_adaptation = false;
  m_adaptation = Context::getSolverProperty<MBool>("adaptation", solver().solverId(), AT_, &m_adaptation);

  m_adapts = true;
  m_adapts = Context::getSolverProperty<MBool>("solverAdapts", solver().solverId(), AT_, &m_adapts);

  m_resTriggeredAdapt = false;
  m_resTriggeredAdapt =
      Context::getSolverProperty<MBool>("resTriggeredAdaptation", solver().solverId(), AT_, &m_resTriggeredAdapt);

  m_noSensors = 0;
  m_noInitialSensors = 0;
  m_rfnBandWidth = nullptr;
  m_sensorInterface = false;
  m_sensorParticle = false;

  if(!m_adaptation) return;

  const MInt maxLevel = Context::getSolverProperty<MInt>("maxRfnmntLvl", solver().solverId(), AT_);

  m_log << "##################################################################" << std::endl;
  m_log << "##################### Adaptation is active #######################" << std::endl;
  m_log << "##################################################################" << std::endl << std::endl;

  solver().m_singleAdaptation = false;
  solver().m_singleAdaptation =
      Context::getSolverProperty<MBool>("singleAdaptationLoop", solver().solverId(), AT_, &solver().m_singleAdaptation);

  ASSERT(Context::propertyLength("sensorType") == Context::propertyLength("sensorWeight"),
         "Lenght of sensorType and sensorWeight not equal.");

  if(Context::propertyExists("sensorType", solver().solverId())) {
    m_noSensors = Context::propertyLength("sensorType", solver().solverId());
  }

  if(m_noSensors > 0 && !Context::propertyExists("sensorWeight", solver().solverId())) {
    TERMM(-1, "sensorWeight property not set!");
  }

  ASSERT(m_noSensors == Context::propertyLength("sensorWeight", solver().solverId()),
         "Length of sensorType " + std::to_string(m_noSensors) + " and sensorWeight "
             + std::to_string(Context::propertyLength("sensorWeight", solver().solverId())) + " not equal.");

  m_saveSensorData = false;
  if(Context::propertyExists("saveSensorData", solver().solverId())) {
    m_saveSensorData = Context::getSolverProperty<MBool>("saveSensorData", solver().solverId(), AT_);
  }

  m_log << "  * Sensors for adaptive mesh refinement for solver " << solver().solverId() << " are active:" << std::endl;

  MInt der = 0;
  for(MInt s = 0; s < m_noSensors; s++) {
    const MString sensor = Context::getSolverProperty<MString>("sensorType", solver().solverId(), AT_, s);
    const MFloat weight = Context::getSolverProperty<MFloat>("sensorWeight", solver().solverId(), AT_, s);
    const MInt level = Context::getSolverProperty<MInt>("maxSensorLevel", solver().solverId(), AT_, &maxLevel, s);

    m_sensorType.push_back(sensor);
    m_sensorWeight.push_back(weight);
    m_maxSensorRefinementLevel.push_back(level);

    m_log << "    - sensor: " << sensor << std::endl;
    m_log << "    - Weight: " << weight << std::endl;
    m_log << "    - maxLevel: " << level << std::endl;

    switch(string2enum(sensor)) {
      case DERIVATIVE: {
        ASSERT(Context::propertyLength("sensorType", solver().solverId()) >= der,
               "sensorDerivativeVariables not stated for this sensor");
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorDerivative);
        MInt derivative = Context::getSolverProperty<MInt>("sensorDerivativeVariables", solver().solverId(), AT_, der);
        m_sensorDerivativeVariables.push_back(derivative);
        der++;
        break;
      }
      case DIVERGENCE:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorDivergence);
        break;
      case TOTALPRESSURE:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorTotalPressure);
        break;
      case ENTROPY_GRADIENT:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorEntropyGrad);
        break;
      case ENTROPY_QUOTIENT:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorEntropyQuot);
        break;
      case MEANSTRESS:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorMeanStress);
        break;
      case VORTICITY:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorVorticity);
        break;
      case INTERFACE:
        ASSERT(s == 0, "Interface sensor has to be the first sensor");
        m_sensorInterface = true;
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorInterface);
        ASSERT((MInt)m_sensorWeight[s] == -1, "Must be discrete sensor!");
        m_noInitialSensors++;
        break;
      case PARTICLE:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorParticle);
        m_sensorParticle = true;
        ASSERT((MInt)m_sensorWeight[s] == -1, "Must be discrete sensor!");
        m_noInitialSensors++;
        break;
      case SPECIES:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorSpecies);
        break;
      case PATCH:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorPatch);
        readPatchProperties();
        // Patch sensor must be second or first, as its also not solution dependand!
        ASSERT(s == 0 || (m_sensorInterface && s == 1), "");
        ASSERT((MInt)m_sensorWeight[s] == -1, "Must be discrete sensor!");
        m_noInitialSensors++;
        break;
      case CUTOFF:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorCutOff);
        m_noInitialSensors++;
        break;
      case SMOOTH:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorSmooth);
        m_noInitialSensors++;
        m_noSmoothingLayers = Context::getSolverProperty<MInt>("smoothingLayers", solver().solverId(), AT_);
        break;
      case BAND:
        m_sensorFnPtr.push_back(&CartesianSolver<nDim, SolverType>::sensorBand);
        m_sensorBandAdditionalLayers = 10;
        m_sensorBandAdditionalLayers = Context::getSolverProperty<MInt>("sensorBandAdditionalLayers", m_solverId, AT_,
                                                                        &m_sensorBandAdditionalLayers);
        ASSERT(s == (m_noSensors - 1), "Please ensure this is the last sensor to be called.");
        break;
      default:
        m_sensorDerivativeVariables.push_back(-1);
        mTerm(1, AT_, "Sensor not found! Check property sensorType!");
    }
  }

  // todo labels:toenhance @Moritz: Switch to sensors and uncomment the Assert below!
  // ASSERT(m_noSensors > 0, "Adaptation without sensors is pointless...");

  m_log << std::endl;
}

/**
 * @brief This sensor generates a max refinement band around the cells with max refinement level. In order for it to
 * work, the property addedAdaptationSteps has to be equal to /maxRefinementLevel() - minLevel()/. This sensor also
 * ensures a smooth transition between levels. Do not use together with sensorSmooth.
 * @author Borja Pedro Beltran
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::sensorBand(std::vector<std::vector<MFloat>>& sensors,
                                                   std::vector<std::bitset<64>>& sensorCellFlag,
                                                   std::vector<MFloat>& sensorWeight,
                                                   MInt sensorOffset,
                                                   MInt sen) {
  TRACE();

  if(this->m_adaptationStep < (maxRefinementLevel() - minLevel())) return;

  MIntScratchSpace markedCells(solver().a_noCells(), AT_, "markedCells");
  MIntScratchSpace bandWidth(this->m_maxSensorRefinementLevel[sen], AT_, "bandWidth");
  markedCells.fill(0);
  bandWidth.fill(0);

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(c_level(cellId) != m_maxSensorRefinementLevel[sen]) continue;
    markedCells[cellId] = 1;
    MInt parentId = c_parentId(cellId);
    while(parentId > -1 && parentId < solver().c_noCells()) {
      markedCells[parentId] = 1;
      parentId = c_parentId(parentId);
    }
  }

  bandWidth[this->m_maxSensorRefinementLevel[sen] - 1] =
      m_bandWidth[this->m_maxSensorRefinementLevel[sen] - 1] + this->m_sensorBandAdditionalLayers;

  for(MInt i = this->m_maxSensorRefinementLevel[sen] - 2; i >= 0; i--) {
    bandWidth[i] = (bandWidth[i + 1] / 2) + this->m_sensorBandAdditionalLayers;
  }

  for(MInt level = minLevel(); level < this->m_maxSensorRefinementLevel[sen]; level++) {
    this->markSurrndCells(markedCells, bandWidth[level], level, true);
  }

  MInt refinedCells = 0;
  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(markedCells(cellId) > 0) {
      const MInt gridCellId = grid().tree().solver2grid(cellId);
      refinedCells++;
      sensors[sensorOffset + sen][gridCellId] = 1;
      sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      MInt parent = c_parentId(cellId);
      if(parent > -1 && parent < solver().a_noCells()) {
        MInt parentGridCellId = grid().tree().solver2grid(parent);
        if(parentGridCellId > -1 && parentGridCellId < grid().raw().m_noInternalCells) {
          sensors[sensorOffset + sen][parentGridCellId] = 1;
          sensorCellFlag[parentGridCellId][sensorOffset + sen] = true;
        }
      }
    }
  }
  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
}

/**\brief simple sensor to apply a limit for a value
 * \author Sven Berger
 *
 * @tparam nDim space dimensions
 * @param sensors array of sensors
 * @param sensorCellFlag sensor set?
 * @param sensorWeight sensor weight
 * @param sensorOffset offset to sensors
 * @param sen sensorid of current sensor
 * @param value acessor function using cellid to value
 * @param limit limit used for value at which to apply adaptation
 */

template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::sensorLimit(std::vector<std::vector<MFloat>>& sensors,
                                                    std::vector<std::bitset<64>>& sensorCellFlag,
                                                    std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen,
                                                    std::function<MFloat(MInt)> value, const MFloat limit,
                                                    const MInt* bandWidth, const MBool refineDiagonals,
                                                    const MBool allowCoarsening) {
  TRACE();

  std::ignore = sensorWeight[sensorOffset];

  MIntScratchSpace markedCells(solver().a_noCells(), AT_, "markedCells");
  markedCells.fill(0);

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(value(cellId) >= limit) {
      markedCells[cellId] = 1;
      MInt parentId = c_parentId(cellId);
      while(parentId > -1 && parentId < solver().a_noCells()) {
        markedCells[parentId] = 1;
        parentId = c_parentId(parentId);
      }
    }
  }

  exchangeData(&markedCells[0], 1);

  for(MInt level = minLevel(); level < m_maxSensorRefinementLevel[sen]; level++) {
    this->markSurrndCells(markedCells, bandWidth[level], level, refineDiagonals);
  }

  for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
    if(markedCells(cellId) == 0) { // coarsen
      if(solver().grid().tree().level(cellId) == minLevel()) continue;
      if(!solver().c_isLeafCell(cellId)) continue;
      if(markedCells[solver().c_parentId(cellId)]) continue;
      if(!allowCoarsening) continue;

      const MInt gridCellId = grid().tree().solver2grid(cellId);
      sensors[sensorOffset + sen][gridCellId] = -1.0;
      sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      //}
    } else {
      if(solver().c_level(cellId) < m_maxSensorRefinementLevel[sen]) { // refine cell
        // if(solver().c_noChildren(cellId) > 0) continue;

        const MInt gridCellId = grid().tree().solver2grid(cellId);
        sensors[sensorOffset + sen][gridCellId] = 1.0;
        sensorCellFlag[gridCellId][sensorOffset + sen] = true;

        MInt parent = solver().c_parentId(cellId);
        if(parent > -1 && parent < solver().c_noCells()) {
          MInt parentGridCellId = grid().tree().solver2grid(parent);
          if(parentGridCellId > -1 && parentGridCellId < solver().grid().raw().m_noInternalCells) {
            sensors[sensorOffset + sen][parentGridCellId] = 1.0;
            sensorCellFlag[parentGridCellId][sensorOffset + sen] = true;
          }
        }
      }
    }
  }
}

/** \brief Saves all sensor values for debug/tunig purposes
 *  \param[in] sensors        Sensor values, must be filled via 'setSensors' before
 *  \param[in] level          Iteration interval of adaption
 *  \param[in] gridFileName   Name of the corresponding grid file
 *  \param[in] recalcIdsTree  Mapping of the recalculated cell ids
 *  If 'saveSensorData' is set, for every adaption a 'sensorData_' file with
 *  corresponding grid file is written. This might serve for debugging various
 *  sensor type or helps tuning specific sensors especially if used in
 *  combination. For each solver an individual file with values of each sensor
 *  is written.
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::saveSensorData(const std::vector<std::vector<MFloat>>& sensors,
                                                       const MInt& level, const MString& gridFileName,
                                                       const MInt* const recalcIdsTree) {
  using namespace maia::parallel_io;
  ASSERT(nDim == 2 || nDim == 3, "wrong number of dimensions!");

  // Update recalc ids
  std::vector<MInt> recalcCellIdsSolver(0);
  MInt noCells;
  MInt noInternalCellIds;
  std::vector<MInt> reorderedCellIds(0);
  this->calcRecalcCellIdsSolver(recalcIdsTree, noCells, noInternalCellIds, recalcCellIdsSolver, reorderedCellIds);

  MLong totalNoInternalCells = -1;
  const MLong longNoInternalCells = noInternalCellIds;
  MPI_Allreduce(&longNoInternalCells, &totalNoInternalCells, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_,
                "longNoInternalCells", "totalNoInternalCells");

  // Open file
  std::stringstream fileName;
  fileName << outputDir() << "sensorData_" << solverId() << "_" << std::to_string(level) << "_" << globalTimeStep
           << ParallelIo::fileExt();
  ParallelIo parallelIo(fileName.str(), PIO_REPLACE, mpiComm());
  // Write header information
  parallelIo.defineScalar(PIO_INT, "noCells");
  parallelIo.setAttribute(solverId(), "solverId");
  for(MUint counter = 0; counter < sensors.size(); counter++) {
    const MString arrayName = "variables" + std::to_string(counter);
    parallelIo.defineArray(PIO_FLOAT, arrayName, totalNoInternalCells);
    const MString varName = "sensor_" + m_sensorType[counter];
    parallelIo.setAttribute(varName, "name", arrayName);
  }
  // Write global attributes
  parallelIo.setAttribute(gridFileName, "gridFile", "");
  parallelIo.setAttribute(solver().time(), "time");
  parallelIo.setAttribute(globalTimeStep, "globalTimeStep");
  // Write scalars
  parallelIo.writeScalar(totalNoInternalCells, "noCells");
  // Set file offsets for field data
  ParallelIo::size_type offset = 0;
  MPI_Exscan(&longNoInternalCells, &offset, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "longNoInternalCells", "offset");
  const ParallelIo::size_type noInternalCells = longNoInternalCells;
  parallelIo.setOffset(noInternalCells, offset);
  // Write sensor values
  const MString name = "variables";
  MInt suffix = 0;
  for(auto&& sensor : sensors) {
    ScratchSpace<MFloat> buffer(noInternalCells, AT_, "buffer");
    for(MInt cellId = 0; cellId < noInternalCells; cellId++) {
      const MInt cellIdRecalc = (recalcIdsTree == nullptr) ? cellId : recalcCellIdsSolver[cellId];
      const MInt cellIdGrid = grid().tree().solver2grid(cellIdRecalc);
      buffer[cellId] = sensor[cellIdGrid];
    }
    parallelIo.writeArray(buffer.data(), name + std::to_string(suffix));
    suffix++;
  }
}


/** \brief sensor to smooth level jumps
 *         NOTE: only refines additional cells to ensure a smooth level transition
 *               this requires that all other sensors are frozen i.e. no refine/coarse sensors set!
 * \author Tim Wegmann
 * \date 2022-08-08
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::sensorSmooth(std::vector<std::vector<MFloat>>& sensors,
                                                     std::vector<std::bitset<64>>& sensorCellFlag,
                                                     std::vector<MFloat>& sensorWeight,
                                                     MInt sensorOffset,
                                                     MInt sen) {
  TRACE();

  static constexpr MInt noDirs = 2 * nDim;
  ASSERT(m_noSmoothingLayers > 0, "No-smoothing layers not specified!");
  static constexpr std::array<MInt, 4> diagDirs{3, 2, 0, 1};

  MFloatScratchSpace markedCells(solver().c_noCells(), AT_, "markedCells");
  MIntScratchSpace markedLevel(solver().c_noCells(), AT_, "markedCells");

  // 1) loop over all refinement levels and set sensor for smooth level transition
  //   top-down:
  for(MInt lvl = m_maxSensorRefinementLevel[sen]; lvl > maxUniformRefinementLevel(); lvl--) {
    markedCells.fill(-1.0);

    for(MInt cellId = 0; cellId < c_noCells(); cellId++) {
      markedLevel[cellId] = c_level(cellId);
    }

    // 2) mark cells at matching level jump
    // const MInt parentLvl = lvl - 1;
    for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
      if(solver().c_level(cellId) != lvl) continue;
      MBool atLvlJump = false;
      for(MInt dir = 0; dir < noDirs; dir++) {
        if(solver().a_hasNeighbor(cellId, dir)) continue;
        // note: cell must have a valid parent as its level is above the minLevel!
        if(solver().a_hasNeighbor(c_parentId(cellId), dir)) {
          atLvlJump = true;
        }
      }
      if(atLvlJump) {
        markedCells(c_parentId(cellId)) = 1.0;
        // markedLevel(c_parentId(cellId)) = parentLvl;
      }
    }

    exchangeData(&markedCells[0], 1);

    // 3) mark band around the level jump on the parent level to ensure smooth transition
    // note: similar to markSurrndCells cells, but different:
    //      in this case the refinement is done top-down instead of bottom-up level-wise!
    //      thus marking across level-jumps must be considered
    for(MInt loopMarker = 1; loopMarker <= m_noSmoothingLayers + 1; loopMarker++) {
      for(MInt cellId = 0; cellId < solver().c_noCells(); cellId++) {
        if((MInt)floor(markedCells(cellId)) != loopMarker) continue;
        const MInt orgLvl = markedLevel[cellId];

        MFloat cellDist = -1;
        if(solver().c_level(cellId) == orgLvl) {
          cellDist = loopMarker + 1;
        } else {
          // note: only mark half as many cells when going down-wards
          ASSERT(solver().c_level(cellId) < orgLvl, "");
          cellDist = loopMarker + 2 * (orgLvl - solver().c_level(cellId));
        }

        // direct neighbors +x-x+y-y+z-z
        for(MInt mainDir = 0; mainDir < noDirs; mainDir++) {
          MInt nghbrId = solver().c_neighborId(cellId, mainDir);
          if(nghbrId < 0) {
            if(c_parentId(cellId) > -1 && solver().a_hasNeighbor(c_parentId(cellId), mainDir)) {
              nghbrId = solver().c_neighborId(c_parentId(cellId), mainDir);
            } else {
              continue;
            }
          }
          const MFloat nghbDist = markedCells(nghbrId);
          if(nghbDist < cellDist) {
            markedCells(nghbrId) = cellDist;
            markedLevel(nghbrId) = orgLvl;
          }
          if(mainDir < 4 && loopMarker < m_noSmoothingLayers) {
            // add diagonal neighbors in x-y plane
            MInt nghbrId2 = solver().c_neighborId(nghbrId, diagDirs.at(mainDir));
            if(nghbrId2 < 0) {
              if(c_parentId(nghbrId) > -1 && solver().a_hasNeighbor(c_parentId(nghbrId), diagDirs.at(mainDir))) {
                nghbrId2 = solver().c_neighborId(c_parentId(nghbrId), diagDirs.at(mainDir));
              }
            }
            if(nghbrId2 > -1) {
              const MFloat nghbDist2 = markedCells(nghbrId2);
              if(nghbDist2 < cellDist) {
                markedCells(nghbrId2) = cellDist;
                markedLevel(nghbrId2) = orgLvl;
              }
            }
            IF_CONSTEXPR(nDim == 3) {
              // add both diagonal neighbors for the direct neighbors
              for(MInt dirZ = 4; dirZ < 6; dirZ++) {
                MInt nghbrIdZ = solver().c_neighborId(nghbrId, dirZ);
                if(nghbrIdZ < 0) {
                  if(c_parentId(nghbrId) > -1 && solver().a_hasNeighbor(c_parentId(nghbrId), dirZ)) {
                    nghbrIdZ = solver().c_neighborId(c_parentId(nghbrId), dirZ);
                  }
                }
                if(nghbrIdZ > -1) {
                  const MFloat nghbDistZ = markedCells(nghbrIdZ);
                  if(nghbDistZ < cellDist) {
                    markedCells(nghbrIdZ) = cellDist;
                    markedLevel(nghbrIdZ) = orgLvl;
                  }
                }
              }
            }
          }
        }
      }
      exchangeData(&markedCells[0], 1);
    }

    // 4) marks cells for refinement
    for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
      if(markedCells(cellId) < 0) continue;
      if(solver().c_level(cellId) < m_maxSensorRefinementLevel[sen]) {
        const MInt gridCellId = grid().tree().solver2grid(cellId);

        // refine
        if(markedLevel[cellId] > solver().c_level(cellId)) {
          sensors[sensorOffset + sen][gridCellId] = 1.0;
          sensorCellFlag[gridCellId][sensorOffset + sen] = true;

          MInt parent = solver().c_parentId(cellId);
          if(parent > -1 && parent < solver().c_noCells()) {
            MInt parentGridCellId = grid().tree().solver2grid(parent);
            if(parentGridCellId > -1 && parentGridCellId < solver().grid().raw().m_noInternalCells) {
              sensors[sensorOffset + sen][parentGridCellId] = 1.0;
              sensorCellFlag[parentGridCellId][sensorOffset + sen] = true;
            }
          }
        }
      }
    }
  }
  sensorWeight[sensorOffset + sen] = m_sensorWeight[sen];
}

/**\brief reOrder cellIds before writing the restart file!
 *        This is necessary for example if the minLevel shall be raised at the new restart!
 * \author Tim Wegmann
 * \date 2020-10-30
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::reOrderCellIds(std::vector<MInt>& reOrderedCells) {
  TRACE();

  if(grid().newMinLevel() < 0) {
    return;
  }

  // see cartesiangrid storeMinLevelCells
  std::map<MLong, MInt> minLevelCells;

  for(MInt cellId = 0; cellId < solver().a_noCells(); cellId++) {
    // allow partition-levelAnchestor minLevel halo-cell!
    if(grid().tree().hasProperty(cellId, Cell::IsHalo)
       && !grid().raw().a_hasProperty(grid().tree().solver2grid(cellId), Cell::IsPartLvlAncestor)) {
      continue;
    }
    if(c_level(cellId) == grid().newMinLevel()) {
      const MLong hilbertId = grid().generateHilbertIndex(cellId, grid().newMinLevel());
      if(minLevelCells.count(hilbertId) > 0) {
        mTerm(1, AT_, "duplicate hilbertId.");
      }
      minLevelCells.insert(std::make_pair(hilbertId, cellId));
    }
  }
  // the new order of mincell levels is stored as the block ids of the new grid block
  const MInt size = minLevelCells.size();
  std::vector<MInt> newMinCellOrder;
  newMinCellOrder.reserve(size);
  for(auto& minLevelCell : minLevelCells) {
    newMinCellOrder.push_back(minLevelCell.second);
  }

  // now generate the full new cell order for writeout
  reOrderedCells.reserve(size);
  for(auto it = newMinCellOrder.begin(); it != newMinCellOrder.end(); it++) {
    reOrderedCells.push_back(*it);
    addChildren(reOrderedCells, *it);
  }
}


/**\brief add childs to reOrdered cellIds
 *        This is necessary for example if the minLevel shall be raised at the new restart!
 * \author Tim Wegmann
 * \date 2020-10-30
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::addChildren(std::vector<MInt>& reOrderedIds, const MInt parentId) {
  for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
    const MInt childId = grid().tree().child(parentId, child);
    // skip if there is no child
    if(childId < 0) continue;
    // add the child to the list
    reOrderedIds.push_back(childId);
    // if the cell has even more children repeat this for those children
    if(grid().tree().noChildren(childId) > 0) {
      addChildren(reOrderedIds, childId);
    }
  }
}
/**\brief reOrder cellIds before writing the restart file!
 *        This is necessary for example if the minLevel shall be raised at the new restart!
 * \author Tim Wegmann
 * \date 2020-10-30
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::recomputeGlobalIds(std::vector<MInt>& reOrderedCells,
                                                           std::vector<MLong>& newGlobalIds) {
  std::vector<MLong> domainOffsets;

  MInt countInternal = 0;
  for(MUint i = 0; i < reOrderedCells.size(); i++) {
    if(solver().a_isHalo(reOrderedCells[i])) continue;
    countInternal++;
  }

  MIntScratchSpace newInternalCells(solver().grid().noDomains(), AT_, "newInternalCells");
  newInternalCells[solver().domainId()] = countInternal;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &newInternalCells[0], 1, type_traits<MInt>::mpiType(),
                solver().mpiComm(), AT_, "MPI_IN_PLACE", "newInternalCells");

  domainOffsets.assign(solver().grid().noDomains() + 1, -1);
  domainOffsets[0] = solver().grid().bitOffset();
  for(MInt d = 1; d < solver().grid().noDomains() + 1; d++) {
    domainOffsets[d] = domainOffsets[d - 1] + (MLong)newInternalCells[d - 1];
  }

  newGlobalIds.clear();
  for(MUint i = 0; i < reOrderedCells.size(); i++) {
    if(solver().a_isHalo(reOrderedCells[i])) continue;
    newGlobalIds.push_back(domainOffsets[solver().domainId()] + i);
  }
}
/**\brief generalised helper function for writing restart files!
 *        This is necessary for example if the minLevel shall be raised at the new restart!
 * \author Tim Wegmann
 * \date 2020-10-30
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::collectVariables(T* variablesIn, ScratchSpace<T>& variablesOut,
                                                         const std::vector<MString>& variablesNameIn,
                                                         std::vector<MString>& variablesNameOut, const MInt noVars,
                                                         const MInt noCells, const MBool reverseOrder) {
  TRACE();

  const MInt variablesOffset = variablesNameOut.size();
  const MInt noTotalVars = noVars + variablesOffset;
  for(MInt j = variablesOffset; j < noTotalVars; j++) {
    const MInt k = j - variablesOffset;
    for(MInt i = 0; i < noCells; i++) {
      if(reverseOrder) {
        variablesOut.p[j * noCells + i] = variablesIn[k * noCells + i];
      } else {
        variablesOut.p[j * noCells + i] = variablesIn[i * noVars + k];
      }
    }
    variablesNameOut.push_back(variablesNameIn[k]);
  }
}

/**\brief generalised helper function for writing restart files!
 *        This is necessary for example if the minLevel shall be raised at the new restart!
 * \author Tim Wegmann
 * \date 2020-10-30
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::collectVariables(T** variablesIn, ScratchSpace<T>& variablesOut,
                                                         const std::vector<MString>& variablesNameIn,
                                                         std::vector<MString>& variablesNameOut, const MInt noVars,
                                                         const MInt noCells) {
  TRACE();

  const MInt variablesOffset = variablesNameOut.size();
  const MInt noTotalVars = noVars + variablesOffset;
  for(MInt j = variablesOffset; j < noTotalVars; j++) {
    const MInt k = j - variablesOffset;
    for(MInt i = 0; i < noCells; i++) {
      variablesOut.p[j * noCells + i] = variablesIn[i][k];
    }
    variablesNameOut.push_back(variablesNameIn[k]);
  }
}

/**
 * \brief This function collects a single parameters for the massivley parallel IO functions
 * \author Tim Wegmann
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::collectParameters(T parametersIn, ScratchSpace<T>& parametersOut,
                                                          const MChar* parametersNameIn,
                                                          std::vector<MString>& parametersNameOut) {
  TRACE();

  const MInt parsOffset = parametersNameOut.size();
  parametersOut.p[parsOffset] = parametersIn;
  parametersNameOut.push_back(parametersNameIn);
}

/// \brief  Derive recalc cell ids of the solver and reordered cell ids
/// \author Unkown (refactoring washing my hands in innocence)
/// \date   30.10.2023
/// \param[in]  recalcIdsTree         Recalculated cell id mapping provided by grid if adpated
/// \param[out] noCells               Number of internal cells including halo partition-levelAnchestor
/// \param[out] noInternalCellIds     Number of internal cells
/// \param[out] recalcCellIdsSolver   Number of internal cells
/// \param[out] reorderedCellIds      Map of reorderd cell ids
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::calcRecalcCellIdsSolver(const MInt* const recalcIdsTree, MInt& noCells,
                                                                MInt& noInternalCellIds,
                                                                std::vector<MInt>& recalcCellIdsSolver,
                                                                std::vector<MInt>& reorderedCellIds) {
  MBool needRecalcCellIdsSolver = (recalcIdsTree != nullptr);
  noCells = noInternalCells();
  noInternalCellIds = noInternalCells();

  if(grid().newMinLevel() > 0) {
    if(domainId() == 0) {
      std::cerr << "Increasing minLevel for solver " << m_solverId << std::endl;
    }
    this->reOrderCellIds(reorderedCellIds);
    MInt countInternal = 0;
    for(MUint i = 0; i < reorderedCellIds.size(); i++) {
      if(grid().tree().hasProperty(reorderedCellIds[i], Cell::IsHalo)) {
        continue;
      }
      countInternal++;
    }
    needRecalcCellIdsSolver = false;
    noCells = reorderedCellIds.size();
    noInternalCellIds = countInternal;
  }

  if(needRecalcCellIdsSolver) {
    recalcCellIdsSolver.resize(grid().tree().size());
    // for multisolver recalc size needs to be raw grid size
    MInt recalcCounter = 0;
    for(MInt cellIdGrid = 0; cellIdGrid < grid().raw().noInternalCells(); cellIdGrid++) {
      if(grid().raw().a_hasProperty(cellIdGrid, Cell::IsHalo)) {
        continue;
      }
      const MInt cellIdSolver = grid().tree().grid2solver(recalcIdsTree[cellIdGrid]);
      if(cellIdSolver > -1) {
        recalcCellIdsSolver[recalcCounter] = cellIdSolver;
        recalcCounter++;
        ASSERT(grid().solverFlag(recalcIdsTree[cellIdGrid], m_solverId), "");
      }
    }
    ASSERT(recalcCounter == grid().noInternalCells(), "recalc ids size is wrong");
  }

#ifdef LS_DEBUG
  MInt internalGCells = 0;
  for(MInt k = 0; k < a_noCells(); ++k) {
    if(grid().tree().hasProperty(k, Cell::IsHalo)) continue;
    ++internalGCells;
  }
  ASSERT(internalGCells == noInternalCells(), "");
#endif
}

/**
 * \brief This function writes the parallel Netcdf cartesian grid cell based solution/restart file
 *        currently used in PostData, LPT and LS solvers!
 * \author Tim Wegmann
 */
template <MInt nDim, class SolverType>
void CartesianSolver<nDim, SolverType>::saveGridFlowVars(
    const MChar* fileName, const MChar* gridFileName, const MInt noTotalCells, const MInt noInternal,
    MFloatScratchSpace& dbVariables, std::vector<MString>& dbVariablesName, MInt noDbVars,
    MIntScratchSpace& idVariables, std::vector<MString>& idVariablesName, MInt noIdVars,
    MFloatScratchSpace& dbParameters, std::vector<MString>& dbParametersName, MIntScratchSpace& idParameters,
    std::vector<MString>& idParametersName, MInt* recalcIds, MFloat time) {
  TRACE();

  MString tmpNcVariablesName = "";
  noDbVars = dbVariablesName.size();
  noIdVars = idVariablesName.size();
  MInt noIdPars = idParametersName.size();
  MInt noDbPars = dbParametersName.size();
  std::vector<MString> ncVariablesName;
  MFloat* tmpDbVariables;
  MInt* tmpIdVariables;

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_REPLACE, solver().mpiComm());

  MLong num = -1;
  MLong longNoInternalCells = (MLong)noInternal;
  MPI_Allreduce(&longNoInternalCells, &num, 1, MPI_LONG, MPI_SUM, solver().mpiComm(), AT_, "noInternalCells", "num");

  for(MInt j = 0; j < noDbVars; j++) {
    tmpNcVariablesName = "variables" + std::to_string(j);
    ncVariablesName.push_back(tmpNcVariablesName);
    parallelIo.defineArray(PIO_FLOAT, ncVariablesName[j], num);
    m_log << dbVariablesName[j].c_str() << std::endl;
    parallelIo.setAttribute(dbVariablesName[j], "name", ncVariablesName[j]);
  }
  MInt k;
  for(MInt j = 0; j < noIdVars; j++) {
    k = j + noDbVars;
    tmpNcVariablesName = "variables" + std::to_string(k);
    ncVariablesName.push_back(tmpNcVariablesName);
    parallelIo.defineArray(PIO_INT, ncVariablesName[k], num);
    parallelIo.setAttribute(idVariablesName[j], "name", ncVariablesName[k]);
  }

  for(MInt j = 0; j < noIdPars; j++) {
    parallelIo.defineScalar(PIO_INT, idParametersName[j]);
  }

  for(MInt j = 0; j < noDbPars; j++) {
    parallelIo.defineScalar(PIO_FLOAT, dbParametersName[j]);
  }

  parallelIo.setAttribute(gridFileName, "gridFile");
  if(g_multiSolverGrid) {
    parallelIo.setAttribute(solver().solverId(), "solverId");
  }
  parallelIo.setAttribute(num, "noCells");
  // TODO labels:IO,toenhance make universal by setting as attribute and update header files in testcases!
  if(time > -1) {
    parallelIo.setAttribute(time, "levelSetTime");
  }

  ParallelIo::size_type start = 0;
  ParallelIo::size_type count = 1;

  MInt noInternalCellIds = noInternal;
  MPI_Exscan(&noInternalCellIds, &start, 1, MPI_INT, MPI_SUM, solver().mpiComm(), AT_, "noInternalCellIds", "start");
  count = noInternalCellIds;

  parallelIo.setOffset(count, start);

  for(MInt j = 0; j < noDbVars; j++) {
    tmpDbVariables = &(dbVariables.p[j * noTotalCells]);
    MFloatScratchSpace tmpDouble(count, AT_, "tmpDouble");
    if(recalcIds != nullptr) {
      for(MInt l = 0; l < count; ++l) {
        tmpDouble[l] = tmpDbVariables[recalcIds[l]];
      }
    } else {
      for(MInt l = 0; l < count; ++l) {
        tmpDouble[l] = tmpDbVariables[l];
      }
    }
    parallelIo.writeArray(tmpDouble.begin(), ncVariablesName[j]);
  }

  for(MInt j = 0; j < noIdVars; j++) {
    tmpIdVariables = &(idVariables.p[j * noTotalCells]);
    MIntScratchSpace tmpInt(count, AT_, "tmpInt");
    if(recalcIds != nullptr) {
      for(MInt l = 0; l < count; ++l) {
        tmpInt[l] = tmpIdVariables[recalcIds[l]];
      }
    } else {
      for(MInt l = 0; l < count; ++l) {
        tmpInt[l] = tmpIdVariables[l];
      }
    }
    parallelIo.writeArray(tmpInt.begin(), ncVariablesName[j + noDbVars]);
  }

  for(MInt j = 0; j < noIdPars; j++) {
    parallelIo.writeScalar(idParameters.p[j], idParametersName[j]);
  }

  for(MInt j = 0; j < noDbPars; j++) {
    parallelIo.writeScalar(dbParameters.p[j], dbParametersName[j]);
  }

  m_log << Scratch::printSelfReport();
  m_log << Scratch::printSelfReport();
}

/**
 * \brief Exchange memory in 'data' assuming a solver size of 'dataBlockSize' per cell
 * \author Lennart Schneiders
 *
 * \tparam        T             Data type
 * \param[in,out] data          Data to exchange
 * \param[in]     dataBlockSize Number of variables of type T per cell, default is 1
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::exchangeData(T* data, const MInt dataBlockSize) {
  TRACE();
  if(noNeighborDomains() == 0) {
    return;
  }

  maia::mpi::exchangeData(solver().grid().neighborDomains(), solver().grid().haloCells(), solver().grid().windowCells(),
                          solver().mpiComm(), data, dataBlockSize);
}

/**
 * \brief Blocking exchange memory in 'data' assuming a solver size of 'dataBlockSize' per cell
 *        NOTE: exchange is only performed on leaf-cells and leaf-NeighborDomains
 *              Assumes, that updateLeafCellExchange has been called in the proxy previously!
 * \author                      Tim Wegmann
 * \tparam        T             Data type
 * \param[in,out] data          Data to exchange
 * \param[in]     noDat         Number of variables of type T per cell, default is 1
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::exchangeLeafData(std::function<T&(MInt, MInt)> data, const MInt noDat) {
  TRACE();

  if(grid().noDomains() < 2) {
    return;
  }

  const MInt tag = 613;
  auto DTYPE = type_traits<T>::mpiType();
  ScratchSpace<T> receiveBuffer(noDat * grid().leafRecSize(), AT_, "windowBuffer");
  ScratchSpace<T> sendBuffer(noDat * grid().leafSendSize(), AT_, "windowBuffer");

  ScratchSpace<MPI_Request> recvRequests(grid().noLeafRecvNeighborDomains(), AT_, "recvRequests");
  ScratchSpace<MPI_Request> sendRequests(grid().noLeafSendNeighborDomains(), AT_, "sendRequests");

  // 1) start receiving
  MInt receiveCount = 0;
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    MPI_Irecv(&(receiveBuffer[receiveCount]), noDat * grid().noLeafHaloCells(d), DTYPE, grid().neighborDomain(d), tag,
              grid().mpiComm(), &recvRequests[n], AT_, "receiveBuffer");
    receiveCount += noDat * grid().noLeafHaloCells(d);
  }

  // 2) fill send buffer
  MInt sendCount = 0;
  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    for(MInt j = 0; j < grid().noLeafWindowCells(d); j++) {
      const MInt cellId = grid().leafWindowCell(d, j);
      for(MInt k = 0; k < noDat; k++) {
        sendBuffer[sendCount] = data(cellId, k);
        sendCount++;
      }
    }
  }

  // 3) start sending
  sendCount = 0;
  for(MInt n = 0; n < grid().noLeafSendNeighborDomains(); n++) {
    const MInt d = grid().leafSendNeighborDomain(n);
    MPI_Isend(&sendBuffer[sendCount], noDat * grid().noLeafWindowCells(d), DTYPE, grid().neighborDomain(d), tag,
              grid().mpiComm(), &sendRequests[n], AT_, "&sendBuffer[sendCount]");
    sendCount += noDat * grid().noLeafWindowCells(d);
  }

  // 4) wait for all send and receive requests to finish
  MPI_Waitall(grid().noLeafRecvNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);

  // 5) scatter date from receive buffers
  MInt recvCount = 0;
  for(MInt n = 0; n < grid().noLeafRecvNeighborDomains(); n++) {
    const MInt d = grid().leafRecvNeighborDomain(n);
    for(MInt j = 0; j < grid().noLeafHaloCells(d); j++) {
      const MInt cellId = grid().leafHaloCell(d, j);
      for(MInt k = 0; k < noDat; k++) {
        data(cellId, k) = receiveBuffer(recvCount);
        recvCount++;
      }
    }
  }

  MPI_Waitall(grid().noLeafSendNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);
}

/**
 * \brief Exchange of sparse data structures on max Level.
 *
 * Given a sparse data structure with getData and setData,
 * and a mapping between the cell id and the data structure index,
 * data is exchange between max level window and halo cells.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param getData
 * \param setData
 * \param dataSize
 * \param cellMapping
 */
template <MInt nDim, class SolverType>
template <class G, class S, class M>
void CartesianSolver<nDim, SolverType>::exchangeSparseLeafValues(G getData, S setData, const MInt dataSize,
                                                                 M cellMapping) {
  TRACE();

  auto* mpi_send_req = new MPI_Request[grid().noNeighborDomains()];
  auto* mpi_recv_req = new MPI_Request[grid().noNeighborDomains()];

  MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");

  MIntScratchSpace sendBuffersInt(grid().noNeighborDomains(), AT_, "sendBuffersInt");
  MIntScratchSpace receiveBuffersInt(grid().noNeighborDomains(), AT_, "receiveBuffersInt");

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    mpi_send_req[i] = MPI_REQUEST_NULL;
    mpi_recv_req[i] = MPI_REQUEST_NULL;
  }

  // 0. gather information about the number of values to be exchanged:
  MInt sendBufferCounter;
  MInt sendBufferCounterOverall = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferCounter = 0;
    sendBufferCounter++;
    for(MInt j = 0; j < grid().noLeafWindowCells(i); j++) {
      const MInt cellId = grid().leafWindowCell(i, j);
      const MInt index = cellMapping(cellId, 0);
      if(index < 0) continue;
      for(MInt d = 0; d < dataSize; d++) {
        ASSERT(!std::isnan(getData(index, d)), grid().tree().globalId(cellId));
      }

      //+check
      sendBufferCounter++;

      //+data
      sendBufferCounter += dataSize;
    }
    sendBufferSize[i] = sendBufferCounter;
    sendBuffersInt[i] = sendBufferCounter;
    sendBufferCounterOverall += sendBufferCounter;
  }

  // 0.a. exchange number data solvers to be exchanged:
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Irecv(&receiveBuffersInt[i], 1, MPI_INT, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_recv_req[i], AT_,
              "receiveBuffers[i]");
  }
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Isend(&sendBuffersInt[i], 1, MPI_INT, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_send_req[i], AT_,
              "sendBuffersInt[i]");
  }
  MPI_Waitall(grid().noNeighborDomains(), mpi_recv_req, MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(grid().noNeighborDomains(), mpi_send_req, MPI_STATUSES_IGNORE, AT_);

  // 0.b setup real exchange framework:
  MInt receiveBufferCounterOverall = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    receiveBufferCounterOverall += receiveBuffersInt[i];
    receiveBufferSize[i] = receiveBuffersInt[i];
  }

  MFloatScratchSpace sendBuffersOverall(sendBufferCounterOverall, AT_, "sendBuffersOverall");
  MFloatScratchSpace receiveBuffersOverall(receiveBufferCounterOverall, AT_, "receiveBuffersOverall");
  MFloatPointerScratchSpace sendBuffers(grid().noNeighborDomains(), AT_, "sendBuffers");
  MFloatPointerScratchSpace receiveBuffers(grid().noNeighborDomains(), AT_, "receiveBuffers");
  MInt counterSend = 0;
  MInt counterReceive = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBuffers.p[i] = &sendBuffersOverall.p[counterSend];
    counterSend += sendBufferSize[i];
    receiveBuffers.p[i] = &receiveBuffersOverall.p[counterReceive];
    counterReceive += receiveBufferSize[i];
  }

  // sicherheitshalber zuruecksetzen
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    mpi_send_req[i] = MPI_REQUEST_NULL;
    mpi_recv_req[i] = MPI_REQUEST_NULL;
  }

  // 1. gather relevant information
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferCounter = 0;
    sendBuffers[i][0] = (MFloat)(-1);
    sendBufferCounter++;
    for(MInt j = 0; j < grid().noLeafWindowCells(i); j++) {
      const MInt cellId = grid().leafWindowCell(i, j);
      const MInt index = cellMapping(cellId, 0);
      if(index < 0) continue;
      // check
      sendBuffers[i][sendBufferCounter] = (MFloat)j;
      sendBufferCounter++;
      // data
      for(MInt d = 0; d < dataSize; d++) {
        ASSERT(!std::isnan(getData(index, d)), "");
        sendBuffers[i][sendBufferCounter] = getData(index, d);
        sendBufferCounter++;
      }
    }
    sendBufferSize[i] = sendBufferCounter;
    sendBuffers[i][0] = (MFloat)(sendBufferCounter);
  }


  // 2. receive data
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MInt bufSize = receiveBufferSize[i];
    MPI_Irecv(receiveBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_recv_req[i],
              AT_, "receiveBuffers[i]");
  }

  // 3. send data
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MInt bufSize = sendBufferSize[i];
    MPI_Isend(sendBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_send_req[i], AT_,
              "sendBuffers[i]");
  }
  MPI_Waitall(grid().noNeighborDomains(), mpi_recv_req, MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(grid().noNeighborDomains(), mpi_send_req, MPI_STATUSES_IGNORE, AT_);

  // 4. store recieved data
  MInt receiveBufferCounter = 0;
  MInt j = -1;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    receiveBufferCounter = 0;
    if(receiveBufferSize[i] != receiveBuffersInt[i]) {
      mTerm(1, AT_, " receiveBufferSize doesn't match expected size from previous communication! ");
    }
    if(receiveBufferSize[i] == 0) {
      m_log << "Warning: empty message from rank " << grid().neighborDomain(i) << std::endl;
    }
    receiveBufferCounter++;
    while(receiveBufferCounter < receiveBufferSize[i]) {
      j = (MInt)receiveBuffers[i][receiveBufferCounter];
      receiveBufferCounter++;
      const MInt cellId = grid().leafHaloCell(i, j);
      MBool skip = false;

      if(cellId > grid().tree().size()) skip = true;
      if(!skip) {
        // add halo-cutCandidate
        // const MInt candId = candidates.size();
        // candidates.emplace_back();
        // candidates[candId].cellId = cellId;

        const MInt index = cellMapping(cellId, 0);

        ASSERT(index > -1, "No corresponding halo cell found!");

        // add infomation from computeNodalvalues:
        for(MInt d = 0; d < dataSize; d++) {
          setData(index, d) = receiveBuffers[i][receiveBufferCounter];
          receiveBufferCounter++;
        }

      } else {
        receiveBufferCounter += dataSize;
      }
    }
  }

  delete[] mpi_send_req;
  delete[] mpi_recv_req;
}

/**
 * \brief Exchange of sparse data structures on max Level.
 *
 * \author Thomas Hoesgen <t.hoesgen@aia.rwth-aachen.de>
 *
 * \param dData
 * \param firstBlock
 * \param dataBlockSize
 * \param mode == 0: Nearest neighbor
 * \param mode == 1: Linear interpolation
 */
template <MInt nDim, class SolverType>
template <typename T>
void CartesianSolver<nDim, SolverType>::exchangeAzimuthalPer(T* data, MInt dataBlockSize, MInt firstBlock) {
  TRACE();

  if(grid().noAzimuthalNeighborDomains() == 0) {
    return;
  }

  MUint sndSize = maia::mpi::getBufferSize(grid().azimuthalWindowCells());
  ScratchSpace<T> windowData(sndSize * dataBlockSize, AT_, "windowData");
  windowData.fill(0);
  MUint rcvSize = maia::mpi::getBufferSize(grid().azimuthalHaloCells());
  ScratchSpace<T> haloData(rcvSize * dataBlockSize, AT_, "haloData");
  haloData.fill(0);

  MInt sndCnt = 0;
  if(std::is_same<T, MInt>::value || std::is_same<T, MBool>::value || std::is_same<T, MLong>::value) {
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        for(MInt b = firstBlock; b < firstBlock + dataBlockSize; b++) {
          windowData[sndCnt] = data[windowId * dataBlockSize + b];
          sndCnt++;
        }
      }
    }
  } else if(std::is_same<T, MFloat>::value) {
    std::function<MBool(const MInt, const MInt)> neighborCheck = [&](const MInt cell, const MInt id) {
      return static_cast<MBool>(grid().tree().hasNeighbor(cell, id));
    };
    std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cell, const MInt id) {
      return static_cast<MFloat>(grid().tree().coordinate(cell, id));
    };
    std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cell, const MInt id) {
      return data[cell * dataBlockSize + id];
    };
    MInt cellCnt = 0;
    for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
      for(MInt j = 0; j < grid().noAzimuthalWindowCells(i); j++) {
        MInt windowId = grid().azimuthalWindowCell(i, j);
        std::array<MFloat, nDim> recCoord;
        std::copy_n(&(solver().m_azimuthalCartRecCoord[cellCnt * nDim]), nDim, &recCoord[0]);

        MInt position = -1;
        const MInt magic_number = 8; // pow(2, nDim);
        std::array<MInt, magic_number> interpolationCells = {0, 0, 0, 0, 0, 0, 0, 0};
        position =
            setUpInterpolationStencil(windowId, interpolationCells.data(), recCoord.data(), neighborCheck, false);
        for(MInt b = firstBlock; b < firstBlock + dataBlockSize; b++) {
          if(position > -1) {
            windowData[sndCnt] =
                interpolateFieldData<true>(&interpolationCells[0], &recCoord[0], b, scalarField, coordinate);
          } else {
            windowData[sndCnt] = scalarField(windowId, b);
          }
          sndCnt++;
        }
        cellCnt++;
      }
    }
  } else {
    mTerm(1, AT_, "Non implemented data type.");
  }

  // Exchange
  maia::mpi::exchangeBuffer(grid().azimuthalNeighborDomains(), grid().azimuthalHaloCells(),
                            grid().azimuthalWindowCells(), mpiComm(), windowData.getPointer(), haloData.getPointer(),
                            dataBlockSize);

  // Extract
  MInt rcvCnt = 0;
  for(MInt i = 0; i < grid().noAzimuthalNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().noAzimuthalHaloCells(i); j++) {
      MInt haloId = grid().azimuthalHaloCell(i, j);

      for(MInt b = firstBlock; b < firstBlock + dataBlockSize; b++) {
        data[haloId * dataBlockSize + b] = haloData[rcvCnt * dataBlockSize + b];
      }
      rcvCnt++;
    }
  }

  // Unmapped Halos
  MBool valueSet;
  const MInt magic_number = 126; // Should be the highest possible count
  MIntScratchSpace nghbrList(magic_number, AT_, "nghbrList");
  for(MInt i = 0; i < grid().noAzimuthalUnmappedHaloCells(); i++) {
    MInt haloId = grid().azimuthalUnmappedHaloCell(i);
    valueSet = false;
    MInt counter = grid().getAdjacentGridCells(haloId, 2, nghbrList, grid().tree().level(haloId), 0);
    for(MInt n = 0; n < counter; n++) {
      MInt nghbrId = nghbrList[n];
      if(nghbrId < 0) {
        continue;
      }
      for(MInt b = firstBlock; b < firstBlock + dataBlockSize; b++) {
        data[haloId * dataBlockSize + b] = data[nghbrId * dataBlockSize + b];
      }
      valueSet = true;
      break;
    }
    if(!valueSet) {
      std::cerr << "Unmapped not set:" << domainId() << " " << haloId << " " << grid().tree().coordinate(haloId, 0)
                << " " << grid().tree().coordinate(haloId, 1) << " " << grid().tree().coordinate(haloId, 2)
                << std::endl;
    }
    ASSERT(valueSet, "No value set!");
  }
}

/** \brief Creates a list of unique corner points for all cells using a hash map
 * levelThreshold optionally specifies the maximum cell level to be extracted
 * bBox optionally specifies a bounding to box to which the extracted domain shall be truncated
 * \author Lennart Schneiders
 * \date 11.01.2013
 */
template <MInt nDim, class Solver>
void CartesianSolver<nDim, Solver>::extractPointIdsFromGrid(Collector<PointBasedCell<nDim>>*& extractedCells,
                                                            Collector<CartesianGridPoint<nDim>>*& gridPoints,
                                                            const MBool extractHaloCells,
                                                            const std::map<MInt, MInt>& splitChildToSplitCell,
                                                            MInt levelThreshold, MFloat* bBox, MBool levelSetMb) const {
  TRACE();

  const MBool isFV =
      (string2enum(solver().solverType()) == MAIA_FINITE_VOLUME) || (string2enum(solver().solverType()) == MAIA_FV_MB);

#ifndef NDEBUG
  if(isFV) {
    for(MInt cellId = 0; cellId < solver().c_noCells(); cellId++) {
      ASSERT(solver().a_hasProperty(cellId, FvCell::IsSplitChild)
                 ? splitChildToSplitCell.count(cellId) == 1
                       && splitChildToSplitCell.find(cellId)->second < solver().c_noCells()
                       && splitChildToSplitCell.find(cellId)->second > -1
                 : true,
             "associated BndryCell for splitChild is missing.");
    }
  }
#endif


  const MInt noCells = solver().c_noCells();
  const MInt maxNoGridPoints = noCells + (noCells / 2);
  const MInt noPoints = IPOW2(nDim);
  const MInt sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
                           {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
  const MFloat DX = solver().c_cellLengthAtLevel(maxRefinementLevel());
  const MFloat EPS = 0.1 * DX;
  ScratchSpace<MFloat> cellLength(maxRefinementLevel() + 1, AT_, "cellLength");
  ScratchSpace<MBool> cellIsActive(noCells, AT_, "isActive");
  for(MInt l = 0; l <= maxRefinementLevel(); l++) {
    cellLength.p[l] = solver().c_cellLengthAtLevel(l);
  }
  MFloat bbox_mem[6];
  if(bBox == nullptr) {
    bBox = &bbox_mem[0];
    solver().geometry().getBoundingBox(bBox);
    for(MInt i = 0; i < nDim; i++) {
      bBox[i] = bBox[i] - F1B2 * cellLength(maxUniformRefinementLevel());
      bBox[nDim + i] = bBox[nDim + i] + F1B2 * cellLength(maxUniformRefinementLevel());
    }
  }
  MFloat bbox2[6];
  for(MInt i = 0; i < nDim; i++) {
    bBox[i] -= EPS;
    bBox[nDim + i] += EPS;
    // bbox2[i] = bBox[i] - cellLength(maxRefinementLevel());// - EPS;
    // bbox2[nDim+i] = bBox[nDim+i] + cellLength(maxRefinementLevel());// + EPS;
    bbox2[i] = bBox[i] - cellLength(maxUniformRefinementLevel());
    bbox2[nDim + i] = bBox[nDim + i] + cellLength(maxUniformRefinementLevel());
    // bBox[i] = bbox2[i];
    // bBox[nDim+i] = bbox2[nDim+i];
  }
  if(levelThreshold < minLevel()) {
    if(domainId() == 0)
      std::cerr << "level threshold reset to minimum refinement level, since it was below." << std::endl;
    levelThreshold = minLevel();
  }
  if(levelSetMb) {
    cellIsActive.fill(false);
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      if(solver().a_isBndryGhostCell(cellId)) continue;
      if(!solver().a_hasProperty(cellId, FvCell::IsOnCurrentMGLevel)) continue;
      if(!solver().a_hasProperty(cellId, FvCell::IsSplitChild) && solver().c_noChildren(cellId) > 0) continue;
      cellIsActive(cellId) = true;
      MInt parentId = solver().a_hasProperty(cellId, FvCell::IsSplitChild) ? splitChildToSplitCell.find(cellId)->second
                                                                           : c_parentId(cellId);
      while(parentId > -1) {
        cellIsActive(parentId) = true;
        parentId = c_parentId(parentId);
      }
    }
  }
  MUlong N[3] = {1, 1, 1};
  for(MInt i = 0; i < nDim; i++)
    N[i] = 1 + (size_t)((bbox2[nDim + i] - bbox2[i]) / DX);
  if(N[0] * N[1] * N[2] > std::numeric_limits<MUlong>::max()) {
    std::cerr << "Warning: MUlong exceeded by hash function!" << std::endl;
  }
  std::unordered_map<size_t, MInt> table;

  //------
  // 0. create collectors
  if(extractedCells) {
    std::cerr << "Warning: extractedCells is not a nullptr pointer as expected." << std::endl;
    delete extractedCells;
    extractedCells = nullptr;
  }
  if(gridPoints) {
    std::cerr << "Warning: gridPoints is not a nullptr pointer as expected." << std::endl;
    delete gridPoints;
    gridPoints = nullptr;
  }
  extractedCells = new Collector<PointBasedCell<nDim>>(noCells, nDim, 0, 0);
  gridPoints = new Collector<CartesianGridPoint<nDim>>(maxNoGridPoints, nDim, 0);
  if(!extractedCells) {
    mTerm(1, AT_, "Allocation of extractedCells failed.");
  }
  if(!gridPoints) {
    mTerm(1, AT_, "Allocation of gridPoints failed.");
  }

  //------
  // 1. determine extracted cells
  MInt noExtractedCells = 0;
  extractedCells->resetSize(0);
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(a_isBndryGhostCell(cellId)) continue;
    if(isFV) {
      if((solver().a_hasProperty(cellId, FvCell::IsSplitChild) || solver().c_noChildren(cellId) == 0)
         && !solver().a_hasProperty(cellId, FvCell::IsOnCurrentMGLevel))
        continue;
      if(!solver().a_hasProperty(cellId, FvCell::IsSplitChild) && solver().c_noChildren(cellId) > 0
         && solver().a_level(cellId) < levelThreshold)
        continue;
      if(solver().a_level(
             solver().a_hasProperty(cellId, FvCell::IsSplitChild) ? splitChildToSplitCell.find(cellId)->second : cellId)
         > levelThreshold)
        continue;
      if(levelSetMb && !cellIsActive(cellId)) continue;
      if(!(/*(grid().azimuthalPeriodicity() && solver().a_isPeriodic(cellId)) ||*/ extractHaloCells)
         && solver().a_isHalo(cellId)) {
        continue;
      }
      if(solver().a_level(
             solver().a_hasProperty(cellId, FvCell::IsSplitChild) ? splitChildToSplitCell.find(cellId)->second : cellId)
         < mMin(maxUniformRefinementLevel(), levelThreshold))
        continue;
      //    if ( a_hasProperty( cellId ,  SolverCell::IsSplitChild ) ) continue;
      if(solver().a_hasProperty(cellId, FvCell::IsSplitCell)) continue;
    }
    //    if ( a_isPeriodic(cellId) ) continue;
    MBool outside = false;
    for(MInt i = 0; i < nDim; i++) {
      if(solver().a_coordinate(cellId, i) < bBox[i] || solver().a_coordinate(cellId, i) > bBox[nDim + i])
        outside = true;
    }
    if(isFV) {
      if(grid().azimuthalPeriodicity() && solver().a_isPeriodic(cellId)) {
        outside = false;
      }
    }
    if(outside) continue;
    extractedCells->append();
    extractedCells->a[noExtractedCells].m_cellId = cellId;
    for(MInt p = 0; p < noPoints; p++) {
      extractedCells->a[noExtractedCells].m_pointIds[p] = -1;
    }
    noExtractedCells++;
  }

  //------
  // 2. determine grid points
  gridPoints->resetSize(0);
  for(MInt c = 0; c < noExtractedCells; c++) {
    const MInt cellId = extractedCells->a[c].m_cellId;
    for(MInt p = 0; p < noPoints; p++) {
      MInt gridPointId = -1;
      MFloat coords[3] = {std::numeric_limits<MFloat>::max(), std::numeric_limits<MFloat>::max(),
                          std::numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized
      size_t index[3] = {0, 0, 0};
      if(isFV) {
        for(MInt i = 0; i < nDim; i++) {
          coords[i] = solver().a_coordinate(cellId, i)
                      + sign[p][i] * F1B2
                            * cellLength.p[solver().a_level(solver().a_hasProperty(cellId, FvCell::IsSplitChild)
                                                                ? splitChildToSplitCell.find(cellId)->second
                                                                : cellId)];
          ASSERT(coords[i] + EPS > bbox2[i] && coords[i] - EPS < bbox2[nDim + i],
                 std::to_string(cellId) << " "
                                        << std::to_string(
                                               solver().a_level(solver().a_hasProperty(cellId, FvCell::IsSplitChild)
                                                                    ? splitChildToSplitCell.find(cellId)->second
                                                                    : cellId))
                                        << " " << std::to_string(i) << " " << std::to_string(EPS) << " "
                                        << std::to_string(coords[i]) << " " << std::to_string(bbox2[i]) << " "
                                        << std::to_string(bbox2[nDim + i]));
          index[i] = (size_t)((coords[i] - bbox2[i] + EPS) / DX);
        }
      } else {
        for(MInt i = 0; i < nDim; i++) {
          coords[i] = solver().a_coordinate(cellId, i) + sign[p][i] * F1B2 * cellLength.p[solver().a_level(cellId)];
          ASSERT(coords[i] + EPS > bbox2[i] && coords[i] - EPS < bbox2[nDim + i],
                 std::to_string(cellId) << " " << std::to_string(solver().a_level(cellId)) << " " << std::to_string(i)
                                        << " " << std::to_string(EPS) << " " << std::to_string(coords[i]) << " "
                                        << std::to_string(bbox2[i]) << " " << std::to_string(bbox2[nDim + i]));
          index[i] = (size_t)((coords[i] - bbox2[i] + EPS) / DX);
        }
      }
      size_t key = index[0] + N[0] * index[1] + N[0] * N[1] * index[2];
      std::pair<typename std::unordered_map<size_t, MInt>::iterator, MBool> ret = table.insert(std::make_pair(key, -1));
      if(ret.second) {
        gridPointId = gridPoints->size();
        ASSERT(gridPointId < maxNoGridPoints, "");
        gridPoints->append();
        gridPoints->a[gridPointId].m_noAdjacentCells = 0;
        for(MInt i = 0; i < noPoints; i++)
          gridPoints->a[gridPointId].m_cellIds[i] = -1;
        for(MInt i = 0; i < nDim; i++)
          gridPoints->a[gridPointId].m_coordinates[i] = coords[i];
        ret.first->second = gridPointId;
      } else {
        ASSERT(gridPointId < maxNoGridPoints, "");
        gridPointId = ret.first->second;
      }
      ASSERT(gridPointId > -1 && gridPointId < maxNoGridPoints, "");
      extractedCells->a[c].m_pointIds[p] = gridPointId;
      MInt rootId =
          (solver().a_hasProperty(cellId, FvCell::IsSplitChild)) ? splitChildToSplitCell.find(cellId)->second : cellId;
      MBool found = false;
      for(MInt n = 0; n < gridPoints->a[gridPointId].m_noAdjacentCells; n++) {
        if(gridPoints->a[gridPointId].m_cellIds[n] == rootId) found = true;
      }
      if(!found) {
        // ASSERT( gridPoints->a[ gridPointId ].m_noAdjacentCells > -1 && gridPoints->a[ gridPointId ].m_noAdjacentCells
        // < IPOW2(nDim), "");
        if(gridPoints->a[gridPointId].m_noAdjacentCells < IPOW2(nDim)) {
          gridPoints->a[gridPointId].m_cellIds[gridPoints->a[gridPointId].m_noAdjacentCells] = rootId;
          gridPoints->a[gridPointId].m_noAdjacentCells++;
        } else
          std::cerr << "Warning: grid point with more than " << IPOW2(nDim)
                    << " neighbor cells: " << gridPoints->a[gridPointId].m_noAdjacentCells << std::endl;
      }
    }
    for(MInt p = 0; p < noPoints; p++) {
      if(extractedCells->a[c].m_pointIds[p] < 0) {
        std::cerr << "Warning: no point for cell " << cellId << " " << p << std::endl;
      }
    }
  }

  //------
  // 3. determine grid point connectivity at domain interfaces
  if(!extractHaloCells) {
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      for(MInt c = 0; c < noHaloCells(d); c++) {
        MInt cellId = haloCellId(d, c);
        if(isFV) {
          if((solver().a_hasProperty(cellId, FvCell::IsSplitChild) || solver().c_noChildren(cellId) == 0)
             && !solver().a_hasProperty(cellId, FvCell::IsOnCurrentMGLevel))
            continue;
          if((!solver().a_hasProperty(cellId, FvCell::IsSplitChild) || solver().c_noChildren(cellId) > 0)
             && solver().a_level(cellId) < levelThreshold)
            continue;
          if(levelSetMb && !cellIsActive(cellId)) continue;
          if(solver().a_hasProperty(cellId, FvCell::IsNotGradient)) continue;
          //      if ( solver().a_hasProperty( cellId ,  FvCell::IsSplitChild ) ) continue;
          if(solver().a_hasProperty(cellId, FvCell::IsSplitCell)) continue;
          if(solver().a_level(solver().a_hasProperty(cellId, FvCell::IsSplitChild)
                                  ? splitChildToSplitCell.find(cellId)->second
                                  : cellId)
             > levelThreshold)
            continue;
        }
        MBool outside = false;
        for(MInt i = 0; i < nDim; i++)
          if(solver().a_coordinate(cellId, i) < bBox[i] || solver().a_coordinate(cellId, i) > bBox[nDim + i])
            outside = true;
        if(outside) continue;
        for(MInt p = 0; p < noPoints; p++) {
          MInt gridPointId = -1;
          MFloat coords[3];
          size_t index[3] = {0, 0, 0};
          if(isFV) {
            for(MInt i = 0; i < nDim; i++) {
              coords[i] = solver().a_coordinate(cellId, i)
                          + sign[p][i] * F1B2
                                * cellLength.p[solver().a_level(solver().a_hasProperty(cellId, FvCell::IsSplitChild)
                                                                    ? splitChildToSplitCell.find(cellId)->second
                                                                    : cellId)];
              ASSERT(coords[i] + EPS > bbox2[i] && coords[i] - EPS < bbox2[nDim + i], "");
              index[i] = (size_t)((coords[i] - bbox2[i] + EPS) / DX);
            }
          } else {
            for(MInt i = 0; i < nDim; i++) {
              coords[i] = solver().a_coordinate(cellId, i) + sign[p][i] * F1B2 * cellLength.p[solver().a_level(cellId)];
              ASSERT(coords[i] + EPS > bbox2[i] && coords[i] - EPS < bbox2[nDim + i], "");
              index[i] = (size_t)((coords[i] - bbox2[i] + EPS) / DX);
            }
          }
          size_t key = index[0] + N[0] * index[1] + N[0] * N[1] * index[2];
          std::pair<typename std::unordered_map<size_t, MInt>::iterator, MBool> ret =
              table.insert(std::make_pair(key, -1));
          if(ret.second)
            continue;
          else {
            ASSERT(gridPointId < maxNoGridPoints, "");
            gridPointId = ret.first->second;
          }
          if(gridPointId < 0) continue;
          ASSERT(gridPointId > -1 && gridPointId < maxNoGridPoints, std::to_string(gridPointId));
          ASSERT(gridPoints->a[gridPointId].m_noAdjacentCells > -1
                     && gridPoints->a[gridPointId].m_noAdjacentCells < IPOW2(nDim),
                 std::to_string(gridPoints->a[gridPointId].m_noAdjacentCells));
          gridPoints->a[gridPointId].m_cellIds[gridPoints->a[gridPointId].m_noAdjacentCells] = cellId;
          gridPoints->a[gridPointId].m_noAdjacentCells++;
        }
      }
    }
  }
}


} // namespace maia

#endif // ifndef CARTESIANSOLVER_H_
