// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CARTESIANGRIDIO_H
#define CARTESIANGRIDIO_H

// needed by pv-plugin
#include "INCLUDE/maiamacro.h"
// needed by pv-plugin
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "GEOM/geometry.h"
#include "INCLUDE/maiaconstants.h"
#include "IO/context.h"
#include "MEMORY/alloc.h"
#include "MEMORY/scratch.h"
#include "UTIL/hilbert.h"
#include "cartesiangridtree.h"
#include "partition.h"

namespace maia {
namespace grid {

template <typename Grid>
class IO {
  Grid& grid_;

  using MInt_dyn_array = typename Grid::MInt_dyn_array;
  using Tree = typename Grid::Tree;
  using Cell = typename Tree::Cell;

  // constexpr variables:
  static constexpr MInt nDim = Grid::nDim_();
  static constexpr MInt m_noDirs = 2 * nDim;
  static constexpr MInt m_maxNoChilds = IPOW2(nDim);

  // references:
  MLong& m_noPartitionCellsGlobal;
  MLong (&m_localPartitionCellOffsets)[3];
  MInt& m_noPartitionCells;
  MLong*& m_localPartitionCellGlobalIds;
  MInt*& m_localPartitionCellLocalIds;
  MInt& m_noMinLevelCellsGlobal;
  MLong& m_noCellsGlobal;
  MInt& m_noInternalCells;
  std::vector<MInt>& m_nghbrDomains;
  MLong*& m_domainOffsets;
  MInt& m_maxLevel;
  MInt& m_minLevel;
  MInt& m_maxRfnmntLvl;
  MFloat& m_lengthLevel0;
  MFloat& m_reductionFactor;
  MFloat (&m_boundingBox)[6];
  MFloat (&m_centerOfGravity)[3];
  MInt& m_maxUniformRefinementLevel;
  MInt& m_decisiveDirection;
  MInt& m_partitionCellOffspringThreshold;
  MFloat& m_partitionCellWorkloadThreshold;
  Tree& m_tree;
  std::map<MLong, MInt>& m_globalToLocalId;
  MInt& m_noPartitionLevelAncestors;
  MLong& m_noPartitionLevelAncestorsGlobal;
  MInt& m_noHaloPartitionLevelAncestors;
  MInt& m_maxPartitionLevelShift;
  MBool& m_loadGridPartition;
  MBool& m_restart;
  MLong& m_32BitOffset;
  std::vector<MInt>& m_partitionLevelAncestorIds;
  std::vector<MLong>& m_partitionLevelAncestorChildIds;
  std::vector<MInt>& m_partitionLevelAncestorNghbrDomains;
  std::vector<std::vector<MInt>>& m_partitionLevelAncestorHaloCells;
  std::vector<std::vector<MInt>>& m_partitionLevelAncestorWindowCells;

  // methods:
  MInt noDomains() const { return grid_.noDomains(); }
  MInt noNeighborDomains() const { return grid_.noNeighborDomains(); }
  MInt noAzimuthalNeighborDomains() const { return grid_.noAzimuthalNeighborDomains(); }
  MPI_Comm mpiComm() const { return grid_.mpiComm(); }
  MInt domainId() const { return grid_.domainId(); }
  MInt solverId() const { return grid_.solverId(); }
  void loadDonorGridPartition(const MLong* const partitionCellsId, const MInt noPartitionCells) {
    grid_.loadDonorGridPartition(partitionCellsId, noPartitionCells);
  }
  MLong generateHilbertIndex(const MInt cellId, const MInt minLevel) {
    return grid_.generateHilbertIndex(cellId, minLevel);
  }
  void resetCell(const MInt& cellId) { grid_.resetCell(cellId); }
  MLong& a_globalId(const MInt cellId) { return grid_.a_globalId(cellId); }
  MLong& a_parentId(const MInt cellId) { return grid_.a_parentId(cellId); }
  MInt& a_level(const MInt cellId) { return grid_.a_level(cellId); }

  MInt a_hasNeighbor(const MInt cellId, const MInt dir) const { return grid_.a_hasNeighbor(cellId, dir); }
  MLong& a_neighborId(const MInt cellId, const MInt dir) { return grid_.a_neighborId(cellId, dir); }
  MFloat& a_coordinate(const MInt cellId, const MInt dir) { return grid_.a_coordinate(cellId, dir); }
  maia::grid::cell::BitsetType::reference a_hasProperty(const MInt cellId, const Cell p) {
    return grid_.a_hasProperty(cellId, p);
  }
  MLong& a_childId(const MInt cellId, const MInt pos) { return grid_.a_childId(cellId, pos); }
  MFloat cellLengthAtLevel(const MInt level) const { return grid_.cellLengthAtLevel(level); }
  MFloat cellLengthAtCell(const MInt cellId) const { return grid_.cellLengthAtCell(cellId); }
  MInt a_noChildren(const MInt cellId) { return grid_.a_noChildren(cellId); }
  MInt& a_noOffsprings(const MInt cellId) { return grid_.a_noOffsprings(cellId); }
  template <class U>
  maia::grid::cell::BitsetType::reference a_hasProperty(Collector<U>* NotUsed(cells_), const MInt cellId,
                                                        const Cell p) {
    return grid_.a_hasProperty(cellId, p);
  }

  template <class U>
  MInt& a_noOffsprings(Collector<U>* NotUsed(cells_), const MInt cellId) {
    return grid_.a_noOffsprings(cellId);
  }
  template <class U>
  MFloat& a_workload(Collector<U>* NotUsed(cells_), const MInt cellId) {
    return grid_.a_workload(cellId);
  }
  MFloat a_weight(const MInt cellId) { return grid_.a_weight(cellId); }
  template <class U>
  MInt a_noCells(Collector<U>* NotUsed(cells_)) {
    return m_tree.size();
  }
  // solvers not having cells in the old collector need their own accessor

  IO(Grid& g)
    : grid_(g),
      m_noPartitionCellsGlobal(grid_.m_noPartitionCellsGlobal),
      m_localPartitionCellOffsets(grid_.m_localPartitionCellOffsets),
      m_noPartitionCells(grid_.m_noPartitionCells),
      m_localPartitionCellGlobalIds(grid_.m_localPartitionCellGlobalIds),
      m_localPartitionCellLocalIds(grid_.m_localPartitionCellLocalIds),
      m_noMinLevelCellsGlobal(grid_.m_noMinLevelCellsGlobal),
      m_noCellsGlobal(grid_.m_noCellsGlobal),
      m_noInternalCells(grid_.m_noInternalCells),
      m_nghbrDomains(grid_.m_nghbrDomains),
      m_domainOffsets(grid_.m_domainOffsets),
      m_maxLevel(grid_.m_maxLevel),
      m_minLevel(grid_.m_minLevel),
      m_maxRfnmntLvl(grid_.m_maxRfnmntLvl),
      m_lengthLevel0(grid_.m_lengthLevel0),
      m_reductionFactor(grid_.m_reductionFactor),
      m_boundingBox(grid_.m_boundingBox),
      m_centerOfGravity(grid_.m_centerOfGravity),
      m_maxUniformRefinementLevel(grid_.m_maxUniformRefinementLevel),
      m_decisiveDirection(grid_.m_decisiveDirection),
      m_partitionCellOffspringThreshold(grid_.m_partitionCellOffspringThreshold),
      m_partitionCellWorkloadThreshold(grid_.m_partitionCellWorkloadThreshold),
      m_tree(grid_.m_tree),
      m_globalToLocalId(grid_.m_globalToLocalId),
      m_noPartitionLevelAncestors(grid_.m_noPartitionLevelAncestors),
      m_noPartitionLevelAncestorsGlobal(grid_.m_noPartitionLevelAncestorsGlobal),
      m_noHaloPartitionLevelAncestors(grid_.m_noHaloPartitionLevelAncestors),
      m_maxPartitionLevelShift(grid_.m_maxPartitionLevelShift),
      m_loadGridPartition(grid_.m_loadGridPartition),
      m_restart(grid_.m_restart),
      m_32BitOffset(grid_.m_32BitOffset),
      m_partitionLevelAncestorIds(grid_.m_partitionLevelAncestorIds),
      m_partitionLevelAncestorChildIds(grid_.m_partitionLevelAncestorChildIds),
      m_partitionLevelAncestorNghbrDomains(grid_.m_partitionLevelAncestorNghbrDomains),
      m_partitionLevelAncestorHaloCells(grid_.m_partitionLevelAncestorHaloCells),
      m_partitionLevelAncestorWindowCells(grid_.m_partitionLevelAncestorWindowCells) {}

 public:
  static void load(Grid& g, MString const& fileName) { IO(g).loadGrid(fileName); }
  template <class CELLTYPE>
  static void save(Grid& g, const MChar* fileName, Collector<CELLTYPE>* cpu_cells,
                   const std::vector<std::vector<MInt>>& cpu_haloCells,
                   const std::vector<std::vector<MInt>>& cpu_windowCells,
                   const std::vector<std::vector<MInt>>& cpu_azimuthalHaloCells,
                   const std::vector<MInt>& cpu_azimuthalUnmappedHaloCells,
                   const std::vector<std::vector<MInt>>& cpu_azimuthalWindowCells, MInt* const recalcIdTree) {
    IO(g).saveGrid(fileName, cpu_cells, cpu_haloCells, cpu_windowCells, cpu_azimuthalHaloCells,
                   cpu_azimuthalUnmappedHaloCells, cpu_azimuthalWindowCells, recalcIdTree);
  }
  template <class CELLTYPE>
  static void calculateNoOffspringsAndWorkload(Grid& g, Collector<CELLTYPE>* cpu_cells, MInt cpu_noCells,
                                               const std::vector<std::vector<MInt>>& cpu_haloCells,
                                               const std::vector<std::vector<MInt>>& cpu_windowCells,
                                               const std::vector<std::vector<MInt>>& partitionLevelAncestorWindowCells,
                                               const std::vector<std::vector<MInt>>& partitionLevelAncestorHaloCells) {
    IO(g).calculateNoOffspringsAndWorkload(cpu_cells, cpu_noCells, cpu_haloCells, cpu_windowCells,
                                           partitionLevelAncestorWindowCells, partitionLevelAncestorHaloCells);
  }
  static void partitionParallel(Grid& g, const MInt tmpCount, const MLong tmpOffset,
                                const MFloat* const partitionCellsWorkload, const MLong* const partitionCellsGlobalId,
                                const MFloat totalWorkload, MLong* const partitionCellOffsets,
                                MLong* const globalIdOffsets, const MBool computeOnlyPartition = false) {
    IO(g).domainPartitioningParallel(tmpCount, tmpOffset, partitionCellsWorkload, partitionCellsGlobalId, totalWorkload,
                                     partitionCellOffsets, globalIdOffsets, computeOnlyPartition);
  }

  static void communicatePartitionLevelAncestorData(Grid& g, const MInt noLocalMinLevelCells,
                                                    const MInt* const localMinLevelCells,
                                                    const MInt* const localMinLevelId,
                                                    const MLong* const minLevelCellsTreeId,
                                                    const MUchar* const cellInfo) {
    IO(g).communicatePartitionLevelAncestorData(noLocalMinLevelCells, localMinLevelCells, localMinLevelId,
                                                minLevelCellsTreeId, cellInfo);
  }
  static void propagateNeighborsFromMinLevelToMaxLevel(Grid& g) { IO(g).propagateNeighborsFromMinLevelToMaxLevel(); }

 private:
  /** \brief Parallel domain partitioning
      \author Lennart Schneiders
      \date July 2017
      // Parallel partitioning based on the scheme: 'Scalable high-quality 1D partitioning', by M. Lieber and W.E. Nagel
     (HPCS 2014)
      // ---
      // The hierarchical algorithm consists of the following steps:
      // 1. The globalIds and workloads of the partition cells are read in parallel, equally distributed, by all domains
      // 2. A coarse partitioning into <noGroups> groups is carried out based on a heuristic algorithm (fully parallel)
     and returns the groupOffsets
      // 3. The local partition cell data is redistributed from the individual domains to the group masters, as
     indicated by the groupOffsets
      // 4. Within each group, the master performs a serial algorithm to find the (locally) optimal partitioning
      // 5. All masters gather the final domain offsets and compute the global workload imbalance
      // 6. All masters broadcast the final domain offsets to their group slaves
      // 7. All masters send the local partition cell global ids to their slaves
      //
      // Note: with the parameter computeOnlyPartition set to true step 7 will be skipped, which will
      // prevent any grid variables to be changed. Thus, a new partitioning is only stored in the
      // passed data arrays.
      // ---
    */
  void domainPartitioningParallel(const MInt tmpCount, const MLong tmpOffset,
                                  const MFloat* const partitionCellsWorkload, const MLong* const partitionCellsGlobalId,
                                  const MFloat totalWorkload, MLong* const partitionCellOffsets,
                                  MLong* const globalIdOffsets, const MBool computeOnlyPartition = false) {
    TRACE();

    using namespace std;

    const MFloat safetyFac = 1.5;
    const MInt maxLocalCells =
        (MInt)(((MFloat)Scratch::getAvailableMemory()) / (safetyFac * ((MFloat)(sizeof(MInt) + sizeof(MFloat)))));
    const MInt maxNoPartitionCellsPerGroup =
        mMin(maxLocalCells, (MInt)IPOW2(16)); // this is a rough estimate, due to unbalanced workloads the actual
                                              // distribution can vary significantly
    const MInt intermediateNoGroups_ = 1 + (m_noPartitionCellsGlobal - 1) / maxNoPartitionCellsPerGroup;
    const MInt noGroups =
        mMax(1, mMin(noDomains() / 24,
                     intermediateNoGroups_)); // make sure there are at least twentyfour ranks per group
    const MInt noRanksPerGroup = noDomains() / noGroups;

    // 1. The partition cell globalIds and workloads are read in parallel, equally distributed, by all domains
    MLongScratchSpace groupOffsets(noGroups + 1, AT_, "groupOffsets");
    const MFloat idealWorkload = totalWorkload / ((MFloat)noDomains());
    MInt groupId = 0;
    MInt groupLocalId = domainId();
    groupOffsets(0) = 0;
    groupOffsets(noGroups) = m_noPartitionCellsGlobal;
    MPI_Comm mpiCommGroup = mpiComm();
    MPI_Comm mpiCommMasters = MPI_COMM_SELF;
    if(noGroups > 1) {
      // 2. A coarse partitioning into noGroups groups is carried out based on a heuristic algorithm (fully parallel)
      // and returns the groupOffsets
      maia::grid::heuristicPartitioningParallel(
          &partitionCellsWorkload[0], tmpCount, tmpOffset, m_noPartitionCellsGlobal, totalWorkload, noDomains(),
          domainId(), noGroups, mpiComm(), mpiCommGroup, mpiCommMasters, groupId, groupLocalId, &groupOffsets[0]);
    }
    MLong groupPartitionOffset0 = groupOffsets[groupId];
    MLong groupPartitionOffset1 = groupOffsets[groupId + 1];
    MBool isGroupMaster = (groupLocalId == 0);
    MInt noGroupRanks = 0;
    MPI_Comm_size(mpiCommGroup, &noGroupRanks);

    MInt groupPartitionCellsCnt = isGroupMaster ? groupPartitionOffset1 - groupPartitionOffset0 : 0;
    ASSERT(groupPartitionCellsCnt >= 0, "");
    MFloatScratchSpace partitionCellsWorkloadLocal(mMax(1, groupPartitionCellsCnt), AT_, "partitionCellsWorkloadLocal");
    MLongScratchSpace partitionCellsIdLocal(mMax(1, groupPartitionCellsCnt), AT_, "partitionCellsIdLocal");
    // 3. The local partition cell data is redistributed from the individual domains to the group masters, as indicated
    // by the groupOffsets
    collectDistributedGroupData(tmpCount, tmpOffset, &partitionCellsGlobalId[0], &partitionCellsWorkload[0], noGroups,
                                &groupOffsets[0], noRanksPerGroup, m_noPartitionCellsGlobal, groupPartitionCellsCnt,
                                &partitionCellsIdLocal[0], &partitionCellsWorkloadLocal[0], mpiComm(), domainId());

    MLongScratchSpace partitionCellOffsetsLocal(noGroupRanks + 1, AT_, "partitionCellOffsetsLocal");
    MLongScratchSpace globalIdOffsetsLocal(noGroupRanks + 1, AT_, "globalIdOffsetsLocal");
    if(isGroupMaster) {
      // 4. The master performs a serial algorithm to find the (locally) optimal partitioning
      MFloat maxWorkloadGroup = maia::grid::optimalPartitioningSerial(
          &partitionCellsWorkloadLocal[0], static_cast<MLong>(groupPartitionCellsCnt), static_cast<MLong>(noGroupRanks),
          &partitionCellOffsetsLocal[0]);
      MFloat maxDeviation = maxWorkloadGroup / idealWorkload;
      MPI_Allreduce(MPI_IN_PLACE, &maxDeviation, 1, MPI_DOUBLE, MPI_MAX, mpiCommMasters, AT_, "MPI_IN_PLACE",
                    "maxDeviation");
      if(noDomains() > 1 && groupId == 0) {
        std::cerr << "global workload imbalance is " << (maxDeviation - 1.0) * 100.0 << "% ";
        m_log << "global workload imbalance is " << (maxDeviation - 1.0) * 100.0 << "%  (using " << noGroups
              << " coarse partition groups, each with at least " << noRanksPerGroup << " ranks);"
              << " ideal workload is " << idealWorkload << " per rank." << std::endl;
      }

      MIntScratchSpace recvCnt(noGroups, AT_, "recvCnt");
      MIntScratchSpace recvDispl(noGroups, AT_, "recvDispl");
      for(MLong i = 0; i < noGroupRanks; i++) {
        globalIdOffsetsLocal[i] = partitionCellsIdLocal[partitionCellOffsetsLocal[i]];
        partitionCellOffsetsLocal[i] += groupPartitionOffset0;
      }
      for(MLong i = 0; i < noGroups; i++) {
        recvDispl[i] = i * noRanksPerGroup;
        recvCnt[i] = (i == noGroups - 1) ? noDomains() - i * noRanksPerGroup : noRanksPerGroup;
      }
      // 5. Group masters assemble the offsets
      MPI_Allgatherv(&globalIdOffsetsLocal[0], noGroupRanks, MPI_LONG, &globalIdOffsets[0], &recvCnt[0], &recvDispl[0],
                     MPI_LONG, mpiCommMasters, AT_, "globalIdOffsetsLocal[0]", "globalIdOffsets[0]");
      MPI_Allgatherv(&partitionCellOffsetsLocal[0], noGroupRanks, MPI_LONG, &partitionCellOffsets[0], &recvCnt[0],
                     &recvDispl[0], MPI_LONG, mpiCommMasters, AT_, "partitionCellOffsetsLocal[0]",
                     "partitionCellOffsets[0]");
      globalIdOffsets[noDomains()] = m_noCellsGlobal + m_32BitOffset;
      partitionCellOffsets[noDomains()] = m_noPartitionCellsGlobal;
      for(MInt i = 0; i < noGroupRanks; i++) {
        partitionCellOffsetsLocal[i] -= groupPartitionOffset0;
      }
    }

    // 6. Group masters send the final offsets to their slaves
    globalIdOffsets[0] = m_32BitOffset;
    partitionCellOffsets[0] = 0;
    MPI_Bcast(&globalIdOffsets[0], noDomains() + 1, MPI_LONG, 0, mpiCommGroup, AT_, "globalIdOffsets[0]");
    MPI_Bcast(&partitionCellOffsets[0], noDomains() + 1, MPI_LONG, 0, mpiCommGroup, AT_, "partitionCellOffsets[0]");

    // Grid variables are changes hereafter, skip if enabled and just return a new partitioning in
    // the passed data arrays to be used e.g. for load balancing without enforcing it directly.
    if(computeOnlyPartition) {
      return;
    }

    // 7. Group masters send the local partition cell global ids to their slaves
    m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];
    if(m_noPartitionCells <= 0) {
      mTerm(1, AT_, "Cannot allocate array with " + std::to_string(m_noPartitionCells) + " elements.");
    }
    mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
           AT_);
    mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, AT_);

    if(isGroupMaster) {
      for(MInt i = 1; i < noGroupRanks; i++) {
        MLong sendOffset = partitionCellOffsetsLocal[i];
        MLong sendCount = partitionCellOffsetsLocal[i + 1] - partitionCellOffsetsLocal[i];
        MPI_Send(&partitionCellsIdLocal[sendOffset], sendCount, MPI_LONG, i, 66, mpiCommGroup, AT_,
                 "partitionCellsIdLocal[sendOffset]");
      }
      MInt sendCount = partitionCellOffsetsLocal[1] - partitionCellOffsetsLocal[0];
      std::copy(&partitionCellsIdLocal[0], &partitionCellsIdLocal[0] + sendCount, &m_localPartitionCellGlobalIds[0]);
    } else {
      MPI_Recv(&m_localPartitionCellGlobalIds[0], m_noPartitionCells, MPI_LONG, 0, 66, mpiCommGroup, MPI_STATUS_IGNORE,
               AT_, "m_localPartitionCellGlobalIds[0]");
    }

    // save the offsets of the local partitionCells
    m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];
    m_localPartitionCellOffsets[0] =
        partitionCellOffsets[domainId()]; // begin of the local partitionCells in the gobal partitionCell array
    m_localPartitionCellOffsets[1] = partitionCellOffsets[domainId() + 1]; // end index of the gobal partitionCell array
    m_localPartitionCellOffsets[2] =
        m_noPartitionCellsGlobal; // end of the local partitionCells in the gobal partitionCell array
    MInt noCells = globalIdOffsets[domainId() + 1] - globalIdOffsets[domainId()];
    m_noInternalCells = noCells;

    std::copy(&globalIdOffsets[0], &globalIdOffsets[0] + noDomains() + 1, &m_domainOffsets[0]);
    for(MInt i = 0; i <= noDomains(); ++i) {
      if(m_domainOffsets[i] < 0 || m_domainOffsets[i] < m_32BitOffset
         || m_domainOffsets[i] > m_noCellsGlobal + m_32BitOffset)
        mTerm(1, AT_, "Invalid domain offsets A.");
      if(i < noDomains() && m_domainOffsets[i] >= m_domainOffsets[i + 1]) {
        cout << "m_domainOffsets[i]: " << m_domainOffsets[i] << std::endl;
        cout << "m_domainOffsets[i+1]: " << m_domainOffsets[i + 1] << std::endl;
        mTerm(1, AT_, "Invalid domain offsets B.");
      }
    }
  }


  //-----------------------------------------------------------------------


  /** \brief Load a grid file
      \author Lennart Schneiders
      \date October 2017
    */
  void loadGrid(const MString& fileName) {
    TRACE();

    using namespace std;

    auto logDuration = [this](const MFloat timeStart, const MString comment) {
      logDuration_(timeStart, "GRID", comment, mpiComm(), domainId(), noDomains());
    };
    const MFloat gridTimeStart = wallTime();

    // Log grid information
    auto logGridInfo = [](const MInt minLevel, const MFloat lengthLevel0, const MFloat* const boundingBox,
                          const MFloat* const centerOfGravity, const MString message) {
      stringstream msg;
      msg << endl
          << "  * " << message << ":" << endl
          << "    - minLevel = " << minLevel << endl
          << "    - boundingBox = ";
      for(MInt i = 0; i < nDim; i++) {
        msg << setprecision(15) << "[" << boundingBox[i] << ", " << boundingBox[nDim + i] << "]";
        if(i < nDim - 1) {
          msg << " x ";
        }
      }
      msg << endl << "         - centerOfGravity = [";
      for(MInt i = 0; i < nDim; i++) {
        msg << setprecision(15) << centerOfGravity[i];
        if(i < nDim - 1) {
          msg << ", ";
        }
      }
      msg << "]" << endl << "    - lengthLevel0 = " << setprecision(15) << lengthLevel0;
      m_log << msg.str() << std::endl;
    };


    // 0. Open the file.
    const MFloat openGridTimeStart = wallTime();
    if(domainId() == 0) std::cerr << "Loading grid file " << fileName << std::endl;
    m_log << "Loading grid file..." << std::endl;
    using namespace maia::parallel_io;
    ParallelIo grid(fileName, PIO_READ, mpiComm());
    logDuration(openGridTimeStart, "Open grid file");


    // 1. Read global attributes
    const MFloat attTimeStart = wallTime();
    MFloat totalWorkload = F0;
    MInt gridDim = -1;
    grid.getAttribute(&gridDim, "nDim");
    grid.getAttribute(&m_noCellsGlobal, "noCells");
    grid.getAttribute(&m_noPartitionCellsGlobal, "noPartitionCells");
    grid.getAttribute(&m_noMinLevelCellsGlobal, "noMinLevelCells");
    grid.getAttribute(&m_minLevel, "minLevel");
    grid.getAttribute(&m_maxLevel, "maxLevel");
    if(grid.hasAttribute("maxUniformRefinementLevel"))
      grid.getAttribute(&m_maxUniformRefinementLevel, "maxUniformRefinementLevel", 1);
    grid.getAttribute(&totalWorkload, "totalWorkload");
    grid.getAttribute(&m_noPartitionLevelAncestorsGlobal, "noPartitionLevelAncestors");
    grid.getAttribute(&m_maxPartitionLevelShift, "maxPartitionLevelShift");
    grid.getAttribute(&m_lengthLevel0, "lengthLevel0");
    grid.getAttribute(&m_centerOfGravity[0], "centerOfGravity", nDim);
    // grid.getAttribute(&m_boundingBox[0], "boundingBox", 2*nDim);
    // grid.getAttribute(&m_reductionFactor, "reductionFactor", 1);
    // grid.getAttribute(&m_decisiveDirection, "decisiveDirection", 1);
    if(grid.hasAttribute("boundingBox")) grid.getAttribute(&m_boundingBox[0], "boundingBox", 2 * nDim);
    if(grid.hasAttribute("reductionFactor")) grid.getAttribute(&m_reductionFactor, "reductionFactor", 1);
    if(grid.hasAttribute("decisiveDirection")) grid.getAttribute(&m_decisiveDirection, "decisiveDirection", 1);
    if(nDim != gridDim) mTerm(1, AT_, "Number of dimensions mismatch with grid file.");
    if(m_maxRfnmntLvl < m_maxLevel) {
      if(grid_.paraViewPlugin()) {
        // No property file read for plugin, thus maxRfnmntLvl is not set
        m_maxRfnmntLvl = m_maxLevel;
      } else {
        TERMM(1, "Error: specified default property maxRfnmntLvl < maxLevel! (" + std::to_string(m_maxRfnmntLvl) + " < "
                     + std::to_string(m_maxLevel) + ")");
      }
    }
    if(m_maxUniformRefinementLevel < m_minLevel || m_maxUniformRefinementLevel > m_maxLevel) {
      mTerm(1, AT_, "Warning: maxUniformRefinementLevel invalid.");
    }

    logGridInfo(m_minLevel, m_lengthLevel0, m_boundingBox, m_centerOfGravity, "grid information read from grid file");

    grid_.m_hasMultiSolverBoundingBox = grid.hasAttribute("multiSolverBoundingBox");
    // Read multisolver grid information if given in the grid file
    if(grid_.m_hasMultiSolverBoundingBox) {
      std::vector<MFloat> multiSolverBoundingBox(2 * nDim);
      grid.getAttribute(&multiSolverBoundingBox[0], "multiSolverBoundingBox", 2 * nDim);
      std::copy_n(&multiSolverBoundingBox[0], 2 * nDim, &(grid_.m_targetGridBoundingBox[0]));

      std::vector<MFloat> multiSolverCoG(nDim);
      grid.getAttribute(&multiSolverCoG[0], "multiSolverCenterOfGravity", nDim);

      MInt multiSolverMinLevel = -1;
      grid.getAttribute(&multiSolverMinLevel, "multiSolverMinLevel");
      TERMM_IF_COND(multiSolverMinLevel < 0,
                    "ERROR: invalid multiSolverMinLevel " + std::to_string(multiSolverMinLevel));

      MFloat multiSolverLength0 = -1.0;
      grid.getAttribute(&multiSolverLength0, "multiSolverLengthLevel0");

      // Determine center of gravity and length level-0 from the given multisolver bounding box
      MFloat lengthLevel0 = 0.0;
      for(MInt i = 0; i < nDim; i++) {
        grid_.m_targetGridCenterOfGravity[i] = 0.5 * (multiSolverBoundingBox[nDim + i] + multiSolverBoundingBox[i]);
        TERMM_IF_NOT_COND(approx(grid_.m_targetGridCenterOfGravity[i], multiSolverCoG[i], MFloatEps),
                          "center of gravity mismatch");

        // Note: same as during grid generation
        const MFloat dist =
            (F1 + F1 / FPOW2(30)) * std::fabs(multiSolverBoundingBox[nDim + i] - multiSolverBoundingBox[i]);
        // const MFloat dist = std::fabs(multiSolverBoundingBox[nDim + i] -
        // multiSolverBoundingBox[i]);
        lengthLevel0 = std::max(lengthLevel0, dist);
      }
      // TODO labels:GRID,IO slice grid multisolver bounding box can be smaller than lengthLevel0!
      TERMM_IF_NOT_COND(approx(lengthLevel0, multiSolverLength0, MFloatEps),
                        "length level 0 mismatch: " + to_string(lengthLevel0) + " != " + to_string(multiSolverLength0));

      grid_.m_targetGridLengthLevel0 = lengthLevel0;
      grid_.m_targetGridMinLevel = multiSolverMinLevel;

      // Perform sanity checks to ensure min-level cells match and the same Hilbert order is
      // obtained
      checkMultiSolverGridExtents(nDim, &m_centerOfGravity[0], m_lengthLevel0, m_minLevel,
                                  &grid_.m_targetGridCenterOfGravity[0], grid_.m_targetGridLengthLevel0,
                                  grid_.m_targetGridMinLevel);

      logGridInfo(grid_.m_targetGridMinLevel, grid_.m_targetGridLengthLevel0, &multiSolverBoundingBox[0],
                  &(grid_.m_targetGridCenterOfGravity[0]), "multisolver grid information from grid file");
    } else {
      // Note: for converted old grid files with bad bounding box/center of gravity information, use
      // the center of gravity read from the grid file for treeIdToCoordinates, it is overwritten
      // later with the values computed from the geometry bounding box
      copy(m_centerOfGravity, m_centerOfGravity + nDim, grid_.m_targetGridCenterOfGravity);

      if(grid.hasAttribute("boundingBox")) {
        copy_n(&m_boundingBox[0], 2 * nDim, &(grid_.m_targetGridBoundingBox[0]));
      }

      grid_.m_targetGridLengthLevel0 = m_lengthLevel0;
      grid_.m_targetGridMinLevel = m_minLevel;
    }

    // Ensure that number of solver in grid file matches number of solvers in grid class
    MInt noSolvers = -1;
    grid.getAttribute(&noSolvers, "noSolvers");

    if(noSolvers < m_tree.noSolvers()) {
      ASSERT(grid_.m_addSolverToGrid, "");
    }

    MLong _32BitOffsetBuffer = m_32BitOffset;
    if(grid.hasAttribute("bitOffset")) {
      grid.getAttribute(&m_32BitOffset, "bitOffset");
    } else {
      m_32BitOffset = 0;
    }
    // If no restartGrid is loaded, all globalIds in the gridFile start with 0. The offset is therefore introduced in
    // the solver as soon as the globalIds assigned for the first time. If a restartGird is loaded, all globalIds start
    // with the 32BitOffset and partitioning takes place with offset.

    logDuration(attTimeStart, "Read attributes etc");


    // 2. partition the grid
    const MFloat partitionTimeStart = wallTime();
    if(domainId() == 0) std::cerr << "  * partition grid on " << noDomains() << " domains... ";
    m_log << "  * partition grid on " << noDomains() << " domains; ";
    mAlloc(m_domainOffsets, noDomains() + 1, "m_domainOffsets", static_cast<MLong>(0), AT_);
    if(noDomains() == 1) { // Serial run
      m_noPartitionCells = m_noPartitionCellsGlobal;
      mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
             AT_);
      mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, AT_);
      grid.setOffset(m_noPartitionCellsGlobal, 0);
      grid.readArray(m_localPartitionCellGlobalIds, "partitionCellsGlobalId");
      m_noInternalCells = m_noCellsGlobal;
      m_domainOffsets[0] = 0;
      m_domainOffsets[1] = m_noCellsGlobal;
      m_domainOffsets[0] += m_32BitOffset;
      m_domainOffsets[1] += m_32BitOffset;

      m_localPartitionCellOffsets[0] = 0;
      m_localPartitionCellOffsets[1] = m_noPartitionCellsGlobal;
      m_localPartitionCellOffsets[2] = m_noPartitionCellsGlobal;
    } else if(m_loadGridPartition) { // Partition partition cells according to target grid
      MLong noDataToRead = (domainId() == 0) ? m_noPartitionCellsGlobal : 0;
      MLongScratchSpace partitionCellsGlobalId(mMax(1L, noDataToRead), AT_, "partitionCellsGlobalId");
      grid.setOffset(noDataToRead, 0);
      grid.readArray(partitionCellsGlobalId.data(), "partitionCellsGlobalId");
      loadDonorGridPartition(&partitionCellsGlobalId[0], m_noPartitionCellsGlobal);
      mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
             AT_);
      mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, AT_);
      grid.setOffset(m_noPartitionCells, m_localPartitionCellOffsets[0]);
      grid.readArray(m_localPartitionCellGlobalIds, "partitionCellsGlobalId");
    } else if(grid_.m_loadPartition) {
      // Load partition from file
      MString gridPartitionFileName = "";
      gridPartitionFileName = Context::getBasicProperty<MString>("partitionFileName", AT_);

      MLongScratchSpace partitionCellOffsets(noDomains() + 1, FUN_, "partitionCellOffsets");
      grid_.loadPartitionFile(gridPartitionFileName, &partitionCellOffsets[0]);
      partitionCellOffsets[noDomains()] = m_noPartitionCellsGlobal;

      m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];
      ASSERT(m_noPartitionCells > 0, "Number of local partition cells needs to be > 0");

      m_localPartitionCellOffsets[0] = partitionCellOffsets[domainId()];
      m_localPartitionCellOffsets[1] = partitionCellOffsets[domainId() + 1];
      m_localPartitionCellOffsets[2] = m_noPartitionCellsGlobal;

      mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
             FUN_);
      mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, FUN_);

      grid.setOffset(m_noPartitionCells, m_localPartitionCellOffsets[0]);
      grid.readArray(m_localPartitionCellGlobalIds, "partitionCellsGlobalId");

      const MLong globalIdOffsetsLocal = (domainId() > 0) ? m_localPartitionCellGlobalIds[0] : 0;
      MPI_Allgather(&globalIdOffsetsLocal, 1, MPI_LONG, m_domainOffsets, 1, MPI_LONG, mpiComm(), AT_,
                    "globalIdOffsetsLocal", "m_domainOffsets");
      m_domainOffsets[noDomains()] = m_noCellsGlobal;

      m_noInternalCells = m_domainOffsets[domainId() + 1] - m_domainOffsets[domainId()];
    } else if(grid_.m_partitionParallelSplit) {
      // Evenly distribute all partition-cells among domains for partitioning the grid
      MLong noLocalPartitionCells = std::floor(m_noPartitionCellsGlobal / noDomains());
      MLong offsetPartitionCells = domainId() * noLocalPartitionCells;
      const MLong missingPartitionCells = m_noPartitionCellsGlobal - noLocalPartitionCells * noDomains();
      if(domainId() < missingPartitionCells) {
        noLocalPartitionCells++;
        offsetPartitionCells += domainId();
      } else {
        offsetPartitionCells += missingPartitionCells;
      }

      // Storage for local partition-cell workloads
      MFloatScratchSpace partitionCellsWorkLoad(noLocalPartitionCells, AT_, "partitionCellsWorkload");

      // Read local part of min-cell workloads from grid
      grid.setOffset(noLocalPartitionCells, offsetPartitionCells);
      grid.readArray(&partitionCellsWorkLoad[0], "partitionCellsWorkload");

      MLongScratchSpace partitionCellOffsets(noDomains() + 1, FUN_, "partitionCellOffsets");
      // Partition min-cells according to workload -> partition cell offsets
      maia::grid::partitionParallelSplit(&partitionCellsWorkLoad[0], noLocalPartitionCells, offsetPartitionCells,
                                         static_cast<MLong>(noDomains()), static_cast<MLong>(domainId()), mpiComm(),
                                         &partitionCellOffsets[0]);
      partitionCellOffsets[noDomains()] = m_noPartitionCellsGlobal;

      m_noPartitionCells = partitionCellOffsets[domainId() + 1] - partitionCellOffsets[domainId()];

      m_localPartitionCellOffsets[0] = partitionCellOffsets[domainId()];
      m_localPartitionCellOffsets[1] = partitionCellOffsets[domainId() + 1];
      m_localPartitionCellOffsets[2] = m_noPartitionCellsGlobal;

      mAlloc(m_localPartitionCellGlobalIds, m_noPartitionCells, "m_localPartitionCellGlobalIds", static_cast<MLong>(-1),
             FUN_);
      mAlloc(m_localPartitionCellLocalIds, m_noPartitionCells, "m_localPartitionCellLocalIds", -1, FUN_);

      grid.setOffset(m_noPartitionCells, m_localPartitionCellOffsets[0]);
      grid.readArray(m_localPartitionCellGlobalIds, "partitionCellsGlobalId");

      const MLong globalIdOffsetsLocal = (domainId() > 0) ? m_localPartitionCellGlobalIds[0] : 0;
      MPI_Allgather(&globalIdOffsetsLocal, 1, MPI_LONG, m_domainOffsets, 1, MPI_LONG, mpiComm(), AT_,
                    "globalIdOffsetsLocal", "m_domainOffsets");
      m_domainOffsets[noDomains()] = m_noCellsGlobal;

      m_noInternalCells = m_domainOffsets[domainId() + 1] - m_domainOffsets[domainId()];
    } else { // Parallel partitioning
      const MInt tmpCount =
          (domainId() == noDomains() - 1)
              ? m_noPartitionCellsGlobal - (noDomains() - 1) * (m_noPartitionCellsGlobal / noDomains())
              : m_noPartitionCellsGlobal / noDomains();
      const MLong tmpOffset = domainId() * (m_noPartitionCellsGlobal / noDomains());
      MFloatScratchSpace partitionCellsWorkload(tmpCount, AT_, "partitionCellsWorkload");
      MLongScratchSpace partitionCellsGlobalId(tmpCount, AT_, "partitionCellsGlobalId");
      MLongScratchSpace partitionCellOffsets(noDomains() + 1, AT_, "partitionCellOffsets");
      MLongScratchSpace globalIdOffsets(noDomains() + 1, AT_, "globalIdOffsets");
      grid.setOffset(tmpCount, tmpOffset);
      grid.readArray(partitionCellsWorkload.data(), "partitionCellsWorkload");
      grid.readArray(partitionCellsGlobalId.data(), "partitionCellsGlobalId");
      domainPartitioningParallel(tmpCount, tmpOffset, &partitionCellsWorkload[0], &partitionCellsGlobalId[0],
                                 totalWorkload, &partitionCellOffsets[0], &globalIdOffsets[0]);
    }
    ASSERT(m_noInternalCells == m_domainOffsets[domainId() + 1] - m_domainOffsets[domainId()], "");
    if(domainId() == 0) std::cerr << std::endl;
    logDuration(partitionTimeStart, "Partition");

    if(domainId() == 0) std::cerr << "  * create data structure for " << m_noCellsGlobal << " cells" << std::endl;
    m_log << "  * create data structure for " << m_noCellsGlobal << " cells" << std::endl;

    // 3. Read cellInfo (8bit unsigned integer, bits 0-3: noChildIds, bits 4-6: position in childIds of parent cell, bit
    // 7: cell is on minLevel)
    const MFloat readCellInfoTimeStart = wallTime();
    std::vector<MUchar> cellInfo(m_noInternalCells);
    grid.setOffset(m_noInternalCells, m_domainOffsets[domainId()] - m_32BitOffset);
    grid.readArray(cellInfo.data(), "cellInfo");
    logDuration(readCellInfoTimeStart, "Read cell info");


    // 4. Check cellInfo for a partition level shift and in this case correct the domain offsets
    const MFloat correctPlsTimeStart = wallTime();
    const MInt noPartitionLevelAncestorsGlobal = correctDomainOffsetsAtPartitionLevelShifts(cellInfo);
    TERMM_IF_COND(noPartitionLevelAncestorsGlobal != m_noPartitionLevelAncestorsGlobal,
                  "Wrong number of partition level ancestors: " + std::to_string(noPartitionLevelAncestorsGlobal)
                      + " != " + std::to_string(m_noPartitionLevelAncestorsGlobal));
    logDuration(correctPlsTimeStart, "Correct offsets PLS");


    // 5. Reset data structure
    m_tree.clear();
    m_tree.append(m_noInternalCells);
    for(MInt i = 0; i < m_noInternalCells; ++i) {
      resetCell(i);
    }

    // Read solver use info from grid
    const MFloat readSolverInfoTimeStart = wallTime();
    {
      MUcharScratchSpace solver(m_noInternalCells, AT_, "solver");
      if(noSolvers == 1 && !g_multiSolverGrid) {
        // If number of solvers is one, all cells belong to it
        std::fill(solver.begin(), solver.end(), 1);
      } else {
        // Otherwise read info from grid
        // Note: the number of internal cells might have changed due to
        // correctDomainOffsetsAtPartitionLevelShifts() and thus we have to set the offset again for
        // reading the solver bits!
        grid.setOffset(m_noInternalCells, m_domainOffsets[domainId()] - m_32BitOffset);
        grid.readArray(&solver[0], "solver");
      }

      for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
        m_tree.solverFromBits(cellId, solver[cellId]);
      }

      // set new solver flag for all cells belonging to the reference solver!
      if(grid_.m_addSolverToGrid) {
        const MInt newSolverId = noSolvers;
        noSolvers++;
        // loop over all cells and set flag to cells which also have the referenceSolverId
        for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
          if(m_tree.solver(cellId, grid_.m_referenceSolver)) {
            m_tree.solver(cellId, newSolverId) = true;
          }
        }
      }
    }
    logDuration(readSolverInfoTimeStart, "Read solver info");

    if(domainId() == 0) {
      if(m_minLevel == m_maxLevel)
        std::cerr << "  * setup grid connectivity on level " << m_minLevel << std::endl;
      else
        std::cerr << "  * setup grid connectivity from level " << m_minLevel << " to " << m_maxLevel << std::endl;
    }
    if(m_minLevel == m_maxLevel)
      m_log << "  * setup grid connectivity on level " << m_minLevel << std::endl;
    else
      m_log << "  * setup grid connectivity from level " << m_minLevel << " to " << m_maxLevel << std::endl;


    // If no restartGrid is loaded, introduce the offset for the first time. Otherwise no further shift is neccessary.
    if(m_32BitOffset > 0)
      m_32BitOffset = 0;
    else
      m_32BitOffset = _32BitOffsetBuffer;

    // 6. Set globalId and mapping
    for(MInt i = 0; i < noDomains() + 1; ++i) {
      m_domainOffsets[i] += m_32BitOffset;
    }
    m_globalToLocalId.clear();
    for(MInt i = 0; i < m_noInternalCells; ++i) {
      a_globalId(i) = i + m_domainOffsets[domainId()];
      m_globalToLocalId[a_globalId(i)] = i;
    }


    // 7. Determine min level cells on this domain (identified by the last bit of cellInfo)
    MInt noMinLevelCells = 0;
    std::vector<MInt> localMinLevelCells;
    ScratchSpace<MInt> localMinLevelId(m_noInternalCells, AT_, "localMinLevelId");
    localMinLevelId.fill(-1);
    for(MInt i = 0; i < m_noInternalCells; ++i) {
      MUint tmpBit = static_cast<MUint>(cellInfo[i]);
      MUint isMinLevel = (tmpBit >> 7) & 1;
      if(isMinLevel) {
        localMinLevelId(i) = noMinLevelCells;
        a_level(i) = m_minLevel;
        a_parentId(i) = -1;
        localMinLevelCells.push_back(i);
        noMinLevelCells++;
      }
    }


    // 8. Read treeId and nghbrIds of local min level cells
    const MFloat readMinCellInfoTimeStart = wallTime();
    MLongScratchSpace minLevelCellsTreeId(mMax(1, noMinLevelCells), AT_, "minLevelCellsTreeId");
    {
      MLongScratchSpace minLevelCellsNghbrIds(noMinLevelCells, m_noDirs, AT_, "minLevelCellsNghbrIds");
      MInt minLevelCellOffset = 0;
      MPI_Exscan(&noMinLevelCells, &minLevelCellOffset, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "noMinLevelCells",
                 "minLevelCellOffset");

      grid.setOffset(noMinLevelCells, minLevelCellOffset);
      grid.readArray(minLevelCellsTreeId.data(), "minLevelCellsTreeId");
      grid.setOffset(m_noDirs * noMinLevelCells, m_noDirs * minLevelCellOffset);
      grid.readArray(minLevelCellsNghbrIds.data(), "minLevelCellsNghbrIds");

      for(MInt i = 0; i < noMinLevelCells; i++) {
        MInt cellId = localMinLevelCells[i];
        for(MInt j = 0; j < m_noDirs; j++) {
          a_neighborId(cellId, j) = (minLevelCellsNghbrIds(i, j) > -1) ? minLevelCellsNghbrIds(i, j) + m_32BitOffset
                                                                       : minLevelCellsNghbrIds(i, j);
        }
      }
    }
    logDuration(readMinCellInfoTimeStart, "Read min cell info");


    // 9. Compute coordinates of local min level cells
    for(MInt i = 0; i < noMinLevelCells; i++) {
      maia::grid::hilbert::treeIdToCoordinates<nDim>(&a_coordinate(localMinLevelCells[i], 0), minLevelCellsTreeId[i],
                                                     static_cast<MLong>(grid_.m_targetGridMinLevel),
                                                     &(grid_.m_targetGridCenterOfGravity[0]),
                                                     grid_.m_targetGridLengthLevel0);
    }
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      m_localPartitionCellGlobalIds[i] += m_32BitOffset;
    }

    // Reintroduced the offset since it is used in communicatePartitionLevelAncestorData.
    if(m_restart) m_32BitOffset = _32BitOffsetBuffer;


    // 10. Exchange partition level ancestors, if any, determine their connectivity, and append them at the back of the
    // cell collector
    const MFloat exchangePlaTimeStart = wallTime();
    {
      // exchange and assemble data of all partition level ancestors
      // afterwards, the tree between the minLevel and the partition level is set (level, childId, noChildren, parentId,
      // coordinate, and noOffsprings)
      communicatePartitionLevelAncestorData(noMinLevelCells, localMinLevelCells.data(), &localMinLevelId[0],
                                            &minLevelCellsTreeId[0], &cellInfo[0]);

      for(MInt i = m_noInternalCells; i < m_tree.size(); ++i) {
        ASSERT(a_hasProperty(i, Cell::IsPartLvlAncestor), "Error: Cell not marked as partition level ancestor.");
        a_hasProperty(i, Cell::IsHalo) = true;
      }
    }
    logDuration(exchangePlaTimeStart, "Exchange partition level ancestors");

    // Set the offset to 0 in case a restartGrid is used. traverseAndSetupSubtree has already the offset of the
    // restartGrid.
    if(m_restart) m_32BitOffset = 0;

    // 11. Recursively setup the local subtrees (level, childIds, noChildren, parentId, coordinates, and noOffspring)
    // starting from each partition cell
    const MFloat setupSubtreesTimeStart = wallTime();
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      MInt cellId = m_localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
      if(a_level(cellId) < m_minLevel || a_level(cellId) > m_minLevel + m_maxPartitionLevelShift) {
        TERMM(1, "Wrong level for partition cell: level " + std::to_string(a_level(cellId)) + ", partitionCellGlobalId "
                     + std::to_string(m_localPartitionCellGlobalIds[i]) + " i" + std::to_string(i));
      }
      traverseAndSetupSubtree(cellId, cellInfo.data());
      a_hasProperty(cellId, Cell::IsPartitionCell) = true;
    }
    logDuration(setupSubtreesTimeStart, "Setup local grid subtrees");

    // Reintroduced the offset since it is used in propagateNeighborsFromMinLevelToMaxLevel()
    if(m_restart) m_32BitOffset = _32BitOffsetBuffer;

    // 12. Set neighborIds on each level > minLevel
    const MFloat neighborsTimeStart = wallTime();
    propagateNeighborsFromMinLevelToMaxLevel();
    logDuration(neighborsTimeStart, "Propagate neighbors");

    // Finally, set the offset to the original value.
    m_32BitOffset = _32BitOffsetBuffer;

    // 13. Convert global to local ids
    if(noDomains() > 1) {
      grid_.createGlobalToLocalIdMapping();
      grid_.changeGlobalToLocalIds();
    }

    for(MInt i = 0; i < m_noPartitionCells; i++) {
      m_localPartitionCellLocalIds[i] = globalIdToLocalId(m_localPartitionCellGlobalIds[i]);
    }


    // Fix for old converted grids with bad bounding box/center of gravity information
    if(!grid_.m_hasMultiSolverBoundingBox && !g_multiSolverGrid) {
      MBool boundingBoxDiff = false;
      for(MInt i = 0; i < 2 * nDim; i++) {
        if(!approx(grid_.m_boundingBoxBackup[i], m_boundingBox[i], MFloatEps)) {
          boundingBoxDiff = true;
        }
      }

      // this part may be removed if transition to new grid format is complete and no grids
      // written with the grid-converter are used anymore.  Remove the following line once
      // labels:GRID transition to new grid format is complete. This hack is necessary as the grid files
      // obtained from the converter have a bad bounding box information
      copy_n(&grid_.m_boundingBoxBackup[0], 2 * nDim, &m_boundingBox[0]);

      // Recompute center of gravity
      MBool centerOfGravityDiff = false;
      for(MInt dir = 0; dir < nDim; dir++) {
        const MFloat tmpCenter = m_boundingBox[dir] + 0.5 * (m_boundingBox[dir + nDim] - m_boundingBox[dir]);

        if(!approx(tmpCenter, m_centerOfGravity[dir], MFloatEps)) {
          centerOfGravityDiff = true;
          m_centerOfGravity[dir] = tmpCenter;
        }
      }
      copy(m_centerOfGravity, m_centerOfGravity + nDim, grid_.m_targetGridCenterOfGravity);

      if(boundingBoxDiff || centerOfGravityDiff) {
        cerr0 << std::endl
              << "WARNING: the grid information from the grid file have been updated/corrected, "
                 "which should only be necessary in case of a converted old grid file with bad "
                 "bounding box/center of gravity information (see m_log for the changes)."
              << std::endl
              << std::endl;
        logGridInfo(m_minLevel, m_lengthLevel0, m_boundingBox, m_centerOfGravity,
                    "updated/corrected grid information with geometry bounding box");
      }
    }

    if(grid_.m_addSolverToGrid && !g_multiSolverGrid) {
      g_multiSolverGrid = true;
    }

    // done
    if(domainId() == 0) std::cerr << "done." << std::endl;
    m_log << "done." << std::endl;
    logDuration(gridTimeStart, "Load grid total");
  }


  //-----------------------------------------------------------------------


  /** \brief Given neighbor information on the minLevel, create neighbor information on all higher levels
      \author Lennart Schneiders
      \date October 2017
  */
  void propagateNeighborsFromMinLevelToMaxLevel() {
    TRACE();

    constexpr MInt noInternalConnections = (nDim == 2) ? 4 : 12;
    constexpr MInt connectionDirs[12] = {1, 1, 3, 3, 1, 1, 3, 3, 5, 5, 5, 5};
    constexpr MInt childs0[12] = {0, 2, 0, 1, 4, 6, 4, 5, 0, 1, 2, 3};
    constexpr MInt childs1[12] = {1, 3, 2, 3, 5, 7, 6, 7, 4, 5, 6, 7};
    constexpr MInt dirStencil[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {2, 2, 3, 3, 2, 2, 3, 3}, {4, 4, 4, 4, 5, 5, 5, 5}};
    constexpr MInt sideIds[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {0, 0, 1, 1, 0, 0, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}};
    constexpr MInt otherSide[2] = {1, 0};
    constexpr MInt revDir[6] = {1, 0, 3, 2, 5, 4};

    std::map<MLong, MInt> exchangeCells;
    for(MInt level = m_minLevel; level < m_maxLevel; ++level) {
      exchangeCells.clear();
      MInt noExchangeCells = 0;
      for(MInt i = 0; i < m_tree.size(); ++i) {
        if(level == a_level(i)) {
          for(MInt d = 0; d < m_noDirs; d++) {
            if(a_hasNeighbor(i, d) > 0 && a_neighborId(i, d) > -1) {
              MLong nghbrId = a_neighborId(i, d);
              if(nghbrId < m_domainOffsets[domainId()] || nghbrId >= m_domainOffsets[domainId() + 1]) {
                if(exchangeCells.count(nghbrId) == 0) {
                  MInt cpu = -1;
                  for(MInt c = 0; c < noDomains(); c++) {
                    if(nghbrId >= m_domainOffsets[c] && nghbrId < m_domainOffsets[c + 1]) {
                      cpu = c;
                      c = noDomains();
                    }
                  }
                  if(cpu < 0 || cpu == domainId()) {
                    mTerm(1, AT_, "Neighbor domain not found " + std::to_string(cpu) + " " + std::to_string(nghbrId));
                  }
                  exchangeCells[nghbrId] = cpu;
                  noExchangeCells++;
                }
              }
            }
          }
        }
      }
      MInt noCellsToQuery = exchangeCells.size();
      MLongScratchSpace queryGlobalIds(mMax(1, noCellsToQuery), AT_, "queryGlobalIds");
      MInt cnt = 0;
      for(auto it = exchangeCells.begin(); it != exchangeCells.end(); it++) {
        queryGlobalIds[cnt] = it->first;
        it->second = cnt;
        cnt++;
      }
      MLongScratchSpace recvChildIds(mMax(1, noCellsToQuery), m_maxNoChilds, AT_, "recvChildIds");
      queryGlobalData(noCellsToQuery, &queryGlobalIds[0], &recvChildIds[0], &a_childId(0, 0), m_maxNoChilds);

      for(MInt i = 0; i < m_tree.size(); ++i) {
        if(level == a_level(i)) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            if(a_childId(i, child) < 0) continue;
            MInt childId = globalIdToLocalId(a_childId(i, child));
            if(childId < 0) continue;
            ASSERT(a_level(childId) == level + 1,
                   "wrong child level cell" + std::to_string(i) + " child" + std::to_string(childId));
            MInt nodeId[3];
            MInt revNodeId[3];
            for(MInt j = 0; j < nDim; j++) {
              nodeId[j] = sideIds[j][child] * IPOW2(j);
              revNodeId[j] = otherSide[sideIds[j][child]] * IPOW2(j);
            }
            for(MInt d = 0; d < nDim; d++) {
              const MInt dir = dirStencil[d][child];
              if(a_hasNeighbor(i, dir) == 0) {
                a_neighborId(childId, dir) = -1;
              } else {
                const MInt childNode = child - nodeId[d] + revNodeId[d];
                const MLong nghbr0 = a_neighborId(i, dir);
                MInt nghbrId = globalIdToLocalId(nghbr0);
                MLong childNghbrId = -1;
                if(nghbrId > -1) {
                  childNghbrId = a_childId(nghbrId, childNode);
                } else {
                  if(exchangeCells.count(nghbr0) == 0) {
                    mTerm(1, AT_, "Exchange cell not found.");
                  }
                  MInt idx = exchangeCells[nghbr0];
                  if(queryGlobalIds[idx] != nghbr0) mTerm(1, AT_, "Cell inconsistency.");
                  childNghbrId = recvChildIds(idx, childNode);
                }
                if(childNghbrId < 0) { // may happen due to cells deleted at the boundaries
                  a_neighborId(childId, dir) = -1;
                } else {
                  a_neighborId(childId, dir) = childNghbrId;
                  MInt childNghbrLocalId = globalIdToLocalId(childNghbrId);
                  if(childNghbrLocalId > -1) {
                    a_neighborId(childNghbrLocalId, revDir[dir]) = a_childId(i, child);
                    if(a_neighborId(childNghbrLocalId, revDir[dir]) >= m_noCellsGlobal + m_32BitOffset)
                      mTerm(1, AT_, "Nghbr inconsistency.");
                  }
                }
              }
            }
          }
          for(MInt c = 0; c < noInternalConnections; c++) {
            MInt dir = connectionDirs[c];
            MLong child0 = a_childId(i, childs0[c]);
            MLong child1 = a_childId(i, childs1[c]);
            MInt localId0 = globalIdToLocalId(child0);
            MInt localId1 = globalIdToLocalId(child1);
            if(localId0 > -1) {
              a_neighborId(localId0, dir) = child1;
              if(a_neighborId(localId0, dir) >= m_noCellsGlobal + m_32BitOffset) mTerm(1, AT_, "Nghbr inconsistency.");
            }
            if(localId1 > -1) {
              a_neighborId(localId1, revDir[dir]) = child0;
              if(a_neighborId(localId1, revDir[dir]) >= m_noCellsGlobal + m_32BitOffset)
                mTerm(1, AT_, "Nghbr inconsistency.");
            }
          }
        }
      }
    }
  }


  //-----------------------------------------------------------------------


  inline MInt globalIdToLocalId(const MLong globalId) { return grid_.globalIdToLocalId(globalId); }


  //-----------------------------------------------------------------------


  MInt findNeighborDomainId(MLong globalId) { return grid_.findNeighborDomainId(globalId); }


  //-----------------------------------------------------------------------


  template <class ITERATOR, typename U>
  MInt findIndex(ITERATOR first, ITERATOR last, const U& val) {
    return grid_.findIndex(first, last, val);
  }


  //-----------------------------------------------------------------------


  /** \brief Collect distributed data located in intData and fltData on the master rank (0) of each communicator/group
      \author Lennart Schneiders
      \date October 2017
    */
  void collectDistributedGroupData(const MInt noLocalData, const MLong localDataOffset, const MLong* const intData,
                                   const MFloat* const fltData, const MInt noGroups, const MLong* const groupOffsets,
                                   const MInt noRanksPerGroup, const MLong noGlobalData, const MInt noRecvData,
                                   MLong* const recvInt, MFloat* const recvFlt, MPI_Comm comm, MInt rank) {
    MInt recvCount = 0;
    MInt locCnt = 0;
    MInt noReqs = 0;
    ScratchSpace<MPI_Request> requests(2 * mMin(noRecvData, noGroups * noRanksPerGroup) + 3 * noGroups, AT_,
                                       "requests");
    MLongScratchSpace sendBuf(noGroups, 2, AT_, "sendBuf");
    requests.fill(MPI_REQUEST_NULL);

    // 1. all ranks send their local data to the corresponding group masters as determined by the groupOffsets
    while(locCnt < noLocalData) {
      MInt sendCnt = 0;
      MInt group = 0;
      while(!(localDataOffset + locCnt >= groupOffsets[group] && localDataOffset + locCnt < groupOffsets[group + 1]))
        group++;
      ASSERT(group < noGroups, "");
      while(locCnt + sendCnt < noLocalData && localDataOffset + locCnt + sendCnt >= groupOffsets[group]
            && localDataOffset + locCnt + sendCnt < groupOffsets[group + 1])
        sendCnt++;
      MInt locOffset = localDataOffset + locCnt - groupOffsets[group];
      MInt ndom = group * noRanksPerGroup;
      if(ndom == rank) {
        ASSERT(locOffset + sendCnt <= noRecvData, "");
        std::copy(&fltData[locCnt], &fltData[locCnt] + sendCnt, &recvFlt[locOffset]);
        std::copy(&intData[locCnt], &intData[locCnt] + sendCnt, &recvInt[locOffset]);
        recvCount += sendCnt;
      } else {
        ASSERT(groupOffsets[group] + locOffset + sendCnt <= noGlobalData, "");
        sendBuf(group, 0) = sendCnt;
        sendBuf(group, 1) = locOffset;
        MPI_Isend(&sendBuf(group, 0), 2, MPI_LONG, ndom, 0, comm, &requests[noReqs++], AT_,
                  "sendBuf(group"); // send data count and offset in local group array
#if defined(HOST_Klogin)
        MPI_Isend(const_cast<MLong*>(&intData[locCnt]), sendCnt, MPI_LONG, ndom, 1, comm, &requests[noReqs++], AT_,
                  "const_cast<MLong*>(&intData[locCnt])"); // send partition cell ids
        MPI_Isend(const_cast<MFloat*>(&fltData[locCnt]), sendCnt, MPI_DOUBLE, ndom, 2, comm, &requests[noReqs++], AT_,
                  "const_cast<MFloat*>(&fltData[locCnt])"); // send partition cell workloads
#else
        MPI_Isend(&intData[locCnt], sendCnt, MPI_LONG, ndom, 1, comm, &requests[noReqs++], AT_,
                  "intData[locCnt]"); // send partition cell ids
        MPI_Isend(&fltData[locCnt], sendCnt, MPI_DOUBLE, ndom, 2, comm, &requests[noReqs++], AT_,
                  "fltData[locCnt]"); // send partition cell workloads
#endif
      }
      locCnt += sendCnt;
    }

    // 2. the group masters receive all partition data residing in their group
    const MFloat time0 = MPI_Wtime();
    MFloat time1 = MPI_Wtime();
    while(recvCount < noRecvData) {
      MPI_Status status;
      MInt flag;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, comm, &flag, &status); // wait for messages (asynchronously)
      if(flag) {
        MInt source = status.MPI_SOURCE;
        MLong recvBuf[2];
        MPI_Recv(recvBuf, 2, MPI_LONG, source, 0, comm, &status, AT_,
                 "recvBuf"); // receive data count and offset in local group array
        MInt recvSize = recvBuf[0];
        MLong locOffset = recvBuf[1];
        ASSERT(locOffset + recvSize <= noRecvData, "");
        MPI_Irecv(&recvInt[locOffset], recvSize, MPI_LONG, source, 1, comm, &requests[noReqs++], AT_,
                  "recvInt[locOffset]"); // receive partition cell ids
        MPI_Irecv(&recvFlt[locOffset], recvSize, MPI_DOUBLE, source, 2, comm, &requests[noReqs++], AT_,
                  "recvFlt[locOffset]"); // receive partition cell workloads
        recvCount += recvSize;
      }
      if((MInt)((MPI_Wtime() - time0) / 10) > (MInt)((time1 - time0) / 10)) {
        std::cerr << "Rank " << rank << " already waiting " << MPI_Wtime() - time0
                  << " seconds for incoming partition data..." << std::endl;
      }
      time1 = MPI_Wtime();
    }

    if(noReqs) MPI_Waitall(noReqs, &requests[0], MPI_STATUSES_IGNORE, AT_); // finish asynchronous communication
  }


  //-----------------------------------------------------------------------


  /** \brief Check for partition level shifts and in this case make sure the first partition level ancestors are on the
     domain with the first corresponding partition cell \author Lennart Schneiders \date October 2017
    */
  MInt correctDomainOffsetsAtPartitionLevelShifts(std::vector<MUchar>& cellInfo) {
    TRACE();

    MIntScratchSpace isPartitionLevelAncestor(m_noInternalCells, AT_, "isPartitionLevelAncestor");
    isPartitionLevelAncestor.fill(1);
    MInt noPartitionLevelAncestorsGlobal = 0;
    MInt noPartitionLevelAncestorsLocal = m_noInternalCells;

    // Loop over all partition cells and flag all offspring cells (flag 0), only partition level
    // ancestor cells will not be traversed (flag 1), the local number of partion level ancestors is
    // then the number of internal cells minus the number of traversed cells
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      MInt cellId = m_localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
      noPartitionLevelAncestorsLocal -=
          traverseAndFlagSubtree(cellId, cellInfo.data(), m_noInternalCells,
                                 &isPartitionLevelAncestor[0]); // flag all parent cells of partition level
    }

    // Determine global number of partition level ancestor cells
    m_noPartitionLevelAncestors = noPartitionLevelAncestorsLocal;
    MPI_Allreduce(&noPartitionLevelAncestorsLocal, &noPartitionLevelAncestorsGlobal, 1, MPI_INT, MPI_SUM, mpiComm(),
                  AT_, "noPartitionLevelAncestorsLocal", "noPartitionLevelAncestorsGlobal");

    if(noPartitionLevelAncestorsGlobal > 0) {
      MInt domainShift = 0;
      MInt lastId = m_noInternalCells - 1;

      // Increase the domain shift as long as the last cell is still a partition level ancestor (no
      // descendant partition cell on this domain!)
      while(isPartitionLevelAncestor[lastId]) {
        domainShift++;
        lastId--;
      }

      if(domainId() == noDomains() - 1 && domainShift) {
        TERMM(1, "The last domain cannot have a domain shift since there is no following domain.");
      }

      // Gather the domain shifts
      MIntScratchSpace domainOffsetsDelta(noDomains(), AT_, "domainOffsetsDelta");
      MPI_Allgather(&domainShift, 1, MPI_INT, domainOffsetsDelta.data(), 1, MPI_INT, mpiComm(), AT_, "domainShift",
                    "domainOffsetsDelta.data()");

      MPI_Request req = MPI_REQUEST_NULL;
      ScratchSpace<MUchar> sendBuf(domainShift, AT_, "sendBuf");

      if(domainShift) {
        // Send the partition level ancestors (without a descendant partition cell on this domain)
        // to the next domain and delete these cells on this domain
        std::copy(&cellInfo[m_noInternalCells - domainShift], &cellInfo[m_noInternalCells - domainShift] + domainShift,
                  &sendBuf[0]);
        MPI_Issend(&sendBuf[0], domainShift, MPI_UNSIGNED_CHAR, domainId() + 1, 657, mpiComm(), &req, AT_,
                   "sendBuf[0]");
        m_noInternalCells -= domainShift;
        cellInfo.resize(m_noInternalCells);
      }
      if(domainId() > 0 && domainOffsetsDelta[domainId() - 1]) {
        // Receive the partition level ancestors from the previous domain and store at the beginning
        // of the cellInfo array
        const MInt delta = domainOffsetsDelta[domainId() - 1];
        cellInfo.resize(m_noInternalCells + delta);
        memmove(&cellInfo[delta], &cellInfo[0], m_noInternalCells * sizeof(MUchar));
        MPI_Recv(&cellInfo[0], delta, MPI_UNSIGNED_CHAR, domainId() - 1, 657, mpiComm(), MPI_STATUS_IGNORE, AT_,
                 "cellInfo[0]");
        m_noInternalCells += delta;
      }

      // Correct global domain offsets
      for(MInt i = 0; i < noDomains(); ++i) {
        m_domainOffsets[i + 1] -= domainOffsetsDelta[i];
      }
      MPI_Wait(&req, MPI_STATUS_IGNORE, AT_);
    }

    return noPartitionLevelAncestorsGlobal;
  }


  //-----------------------------------------------------------------------


  /** \brief The subtrees between the local minLevel and the partitionLevel is gathered on the domain which holds the
     corresponding minLevel cell. This domains reconstructs the subtree connectivity, given the collected data and sends
     the connectivity informaton back to the various domains. \author Lennart Schneiders \date October 2017
    */
  void communicatePartitionLevelAncestorData(const MInt noLocalMinLevelCells, const MInt* const localMinLevelCells,
                                             const MInt* const localMinLevelId, const MLong* const minLevelCellsTreeId,
                                             const MUchar* const cellInfo) {
    TRACE();

    using namespace std;

    m_partitionLevelAncestorIds.clear();
    m_partitionLevelAncestorChildIds.clear();
    m_partitionLevelAncestorNghbrDomains.clear();
    m_partitionLevelAncestorHaloCells.clear();
    m_partitionLevelAncestorWindowCells.clear();
    m_noHaloPartitionLevelAncestors = 0;

    if(m_noPartitionLevelAncestorsGlobal == 0) {
      return;
    }

    // Mark partiton cells
    MIntScratchSpace isPartitionCell(m_noInternalCells, AT_, "isPartitionCell");
    isPartitionCell.fill(0);
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      const MInt cellId = m_localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
      isPartitionCell(cellId) = 1;
    }

    // Determine global id of first min level cell on this domain
    MLong firstMinLevelCell = std::numeric_limits<MLong>::max();
    for(MInt i = 0; i < noLocalMinLevelCells; i++) {
      firstMinLevelCell = mMin(firstMinLevelCell, a_globalId(localMinLevelCells[i]));
    }

    // Gather first min level cell ids of all domains
    MLongScratchSpace minLevelCellGlobalIdOffsets(noDomains() + 1, AT_, "minLevelCellGlobalIdOffsets");
    MPI_Allgather(&firstMinLevelCell, 1, MPI_LONG, minLevelCellGlobalIdOffsets.data(), 1, MPI_LONG, mpiComm(), AT_,
                  "firstMinLevelCell", "minLevelCellGlobalIdOffsets.data()");
    minLevelCellGlobalIdOffsets[noDomains()] = m_noCellsGlobal + m_32BitOffset;

    for(MInt cpu = noDomains() - 1; cpu >= 0; cpu--) { // some domains may not have any min level cells at all
      if(minLevelCellGlobalIdOffsets[cpu] > m_noCellsGlobal + m_32BitOffset) {
        minLevelCellGlobalIdOffsets[cpu] = minLevelCellGlobalIdOffsets[cpu + 1];
      }
    }

    std::vector<std::tuple<MLong, MInt, MInt>> partitionLevelAncestorsOnMinLevel;
    std::vector<std::tuple<MLong, MInt, MInt>> partitionLevelAncestorsOther;
    MIntScratchSpace isPartitionLevelAncestor(m_noInternalCells, AT_, "isPartitionLevelAncestor");
    isPartitionLevelAncestor.fill(1);

    MInt noPartitionLevelAncestorsGlobal = 0;
    MInt noPartitionLevelAncestorsLocal = m_noInternalCells;
    TERMM_IF_NOT_COND(m_noInternalCells == m_tree.size(),
                      "Error: tree should only contain internal cells at this moment.");

    // Determine partition level ancestor cells and set the corresponding cell property
    for(MInt i = 0; i < m_noPartitionCells; i++) {
      const MInt cellId = m_localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
      noPartitionLevelAncestorsLocal -=
          traverseAndFlagSubtree(cellId, cellInfo, m_tree.size(),
                                 &isPartitionLevelAncestor[0]); // flag all parent cells of partition level
    }

    m_noPartitionLevelAncestors = noPartitionLevelAncestorsLocal;
    // Determine global number of partition level ancestor cells
    MPI_Allreduce(&noPartitionLevelAncestorsLocal, &noPartitionLevelAncestorsGlobal, 1, MPI_INT, MPI_SUM, mpiComm(),
                  AT_, "noPartitionLevelAncestorsLocal", "noPartitionLevelAncestorsGlobal");
    m_noPartitionLevelAncestorsGlobal = noPartitionLevelAncestorsGlobal;

    for(MInt i = 0; i < m_noInternalCells; i++) {
      a_hasProperty(i, Cell::IsPartLvlAncestor) = isPartitionLevelAncestor(i);
    }

    for(MInt i = 0; i < m_noPartitionCells; i++) {
      const MInt cellId = m_localPartitionCellGlobalIds[i] - m_domainOffsets[domainId()];
      if(localMinLevelId[cellId] < 0) {
        isPartitionLevelAncestor[cellId] = 1; // also flag partition cells which are not on min level
      }
    }

    for(MInt cellId = 0; cellId < m_noInternalCells; cellId++) {
      const MLong globalId = a_globalId(cellId);
      if(isPartitionLevelAncestor(cellId)) {
        if(localMinLevelId[cellId] > -1) {
          // Partition level ancestor is a min cell
          partitionLevelAncestorsOnMinLevel.push_back(std::make_tuple(globalId, cellId, domainId()));
        } else {
          // Partition level ancestor is not a min-cell, find the domain which has the corresponding
          // parent min-cell
          MInt ndom = -1;
          for(MInt cpu = 0; cpu < noDomains(); cpu++) {
            if(globalId >= minLevelCellGlobalIdOffsets[cpu] && globalId < minLevelCellGlobalIdOffsets[cpu + 1]) {
              ndom = cpu;
              break;
            }
          }
          TERMM_IF_COND(ndom < 0 || ndom >= noDomains(), "domain for globalId not found " + std::to_string(ndom));
          partitionLevelAncestorsOther.push_back(std::make_tuple(globalId, cellId, ndom));
        }
      }
    }

    MIntScratchSpace noQuery(noDomains(), AT_, "noQuery");
    MIntScratchSpace queryOffsets(noDomains() + 1, AT_, "queryOffsets");
    noQuery.fill(0);

    // Count the number of queries per domain
    for(MUint i = 0; i < partitionLevelAncestorsOther.size(); i++) {
      noQuery[get<2>(partitionLevelAncestorsOther[i])]++;
    }
    queryOffsets(0) = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      queryOffsets(c + 1) = queryOffsets(c) + noQuery[c];
    }

    MLongScratchSpace queryGlobalId(queryOffsets[noDomains()], AT_, "queryGlobalId");
    ScratchSpace<MUchar> queryCellInfo(queryOffsets[noDomains()], AT_, "queryCellInfo");
    MIntScratchSpace queryIsPartitionCell(queryOffsets[noDomains()], AT_, "queryIsPartitionCell");
    noQuery.fill(0);

    // Gather the local partition level ancestor data that needs to be communicated
    for(MUint i = 0; i < partitionLevelAncestorsOther.size(); i++) {
      MLong globalId = get<0>(partitionLevelAncestorsOther[i]);
      MInt cellId = get<1>(partitionLevelAncestorsOther[i]);
      MInt cpu = get<2>(partitionLevelAncestorsOther[i]);
      queryGlobalId(queryOffsets(cpu) + noQuery(cpu)) = globalId;
      queryCellInfo(queryOffsets(cpu) + noQuery(cpu)) = cellInfo[cellId];
      queryIsPartitionCell(queryOffsets(cpu) + noQuery(cpu)) = isPartitionCell[cellId];
      // it->second = queryOffsets(cpu) + noQuery(cpu);
      noQuery[cpu]++;
    }

    MIntScratchSpace noCollect(noDomains(), AT_, "noCollect");
    MIntScratchSpace collectOffsets(noDomains() + 1, AT_, "collectOffsets");

    // Communicate the number of queries of all domains with each other
    MPI_Alltoall(&noQuery[0], 1, MPI_INT, &noCollect[0], 1, MPI_INT, mpiComm(), AT_, "noQuery[0]", "noCollect[0]");
    // Compute offsets
    collectOffsets(0) = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      collectOffsets(c + 1) = collectOffsets(c) + noCollect[c];
    }

    // Buffers for received data
    MLongScratchSpace collectGlobalId(collectOffsets[noDomains()], AT_, "collectGlobalId");
    ScratchSpace<MUchar> collectCellInfo(collectOffsets[noDomains()], AT_, "collectCellInfo");
    MIntScratchSpace collectIsPartitionCell(collectOffsets[noDomains()], AT_, "collectIsPartitionCell");
    ScratchSpace<MPI_Request> queryReq(3 * noDomains(), AT_, "queryReq");
    queryReq.fill(MPI_REQUEST_NULL);

    // Send the local partition level ancestor data to the domain with the corresponding min-cells
    MInt queryCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(noQuery[c] == 0) continue;
      MPI_Issend(&queryGlobalId[queryOffsets[c]], noQuery[c], MPI_LONG, c, 456, mpiComm(), &queryReq[queryCnt++], AT_,
                 "queryGlobalId[queryOffsets[c]]");
      MPI_Issend(&queryCellInfo[queryOffsets[c]], noQuery[c], MPI_UNSIGNED_CHAR, c, 457, mpiComm(),
                 &queryReq[queryCnt++], AT_, "queryCellInfo[queryOffsets[c]]");
      MPI_Issend(&queryIsPartitionCell[queryOffsets[c]], noQuery[c], MPI_INT, c, 458, mpiComm(), &queryReq[queryCnt++],
                 AT_, "queryIsPartitionCell[queryOffsets[c]]");
    }

    // Receive the data if there is any
    MInt collectCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(noCollect[c] == 0) continue;
      MPI_Recv(&collectGlobalId[collectOffsets[c]], noCollect[c], MPI_LONG, c, 456, mpiComm(), MPI_STATUS_IGNORE, AT_,
               "collectGlobalId[collectOffsets[c]]");
      MPI_Recv(&collectCellInfo[collectOffsets[c]], noCollect[c], MPI_UNSIGNED_CHAR, c, 457, mpiComm(),
               MPI_STATUS_IGNORE, AT_, "collectCellInfo[collectOffsets[c]]");
      MPI_Recv(&collectIsPartitionCell[collectOffsets[c]], noCollect[c], MPI_INT, c, 458, mpiComm(), MPI_STATUS_IGNORE,
               AT_, "collectIsPartitionCell[collectOffsets[c]]");
      collectCnt++;
    }

    if(queryCnt > 0) MPI_Waitall(queryCnt, &queryReq[0], MPI_STATUSES_IGNORE, AT_);

    // Collect the local min-level partition level ancestor data
    std::vector<std::tuple<MLong, MUchar, MInt>> collectData;
    for(MUint i = 0; i < partitionLevelAncestorsOnMinLevel.size(); i++) {
      MLong globalId = get<0>(partitionLevelAncestorsOnMinLevel[i]);
      MInt localId = get<1>(partitionLevelAncestorsOnMinLevel[i]);
      ASSERT(domainId() == get<2>(partitionLevelAncestorsOnMinLevel[i]), "");
      collectData.push_back(std::make_tuple(globalId, cellInfo[localId], isPartitionCell[localId]));
    }

    // ... and that of all descendent partition level ancestors (local and on other domains)
    for(MInt c = 0; c < noDomains(); c++) {
      for(MInt d = 0; d < noCollect[c]; d++) {
        MInt id = collectOffsets[c] + d;
        collectData.push_back(std::make_tuple(collectGlobalId[id], collectCellInfo[id], collectIsPartitionCell[id]));
      }
    }

    sort(collectData.begin(), collectData.end()); // sort by globalId to retain depth-first ordering starting from min
                                                  // level, truncated at partition level

    constexpr MInt noData = 5 + m_maxNoChilds + m_noDirs;
    MLongScratchSpace collectMinLevelTreeId(collectData.size(), AT_, "collectMinLevelTreeId");
    MLongScratchSpace collectParentId(collectData.size(), AT_, "collectParentId");
    MIntScratchSpace collectLevel(collectData.size(), AT_, "collectLevel");
    MIntScratchSpace collectNoChildren(collectData.size(), AT_, "collectNoChildren");
    MIntScratchSpace collectNoOffspring(collectData.size(), AT_, "collectNoOffspring");
    MLongScratchSpace collectChildIds(collectData.size(), m_maxNoChilds, AT_, "collectChildIds");
    MLongScratchSpace collectNghbrIds(collectData.size(), m_noDirs, AT_, "collectNghbrIds");
    std::vector<MInt> newCells;
    MInt newCellsMaxLevel = m_minLevel;
    std::map<MLong, MInt> collectMap;
    collectMinLevelTreeId.fill(-1);
    collectParentId.fill(-1);
    collectLevel.fill(-1);
    collectNoChildren.fill(-1);
    collectNoOffspring.fill(-1);
    collectChildIds.fill(-1);
    collectNghbrIds.fill(-1);

    // Create mapping: part. lvl ancestor globalId -> local position in collectData
    for(MUint i = 0; i < collectData.size(); i++) {
      collectMap[get<0>(collectData[i])] = i;
    }

    // Collect data of min-level partition level ancestors
    for(MUint i = 0; i < partitionLevelAncestorsOnMinLevel.size(); i++) {
      MLong globalId = get<0>(partitionLevelAncestorsOnMinLevel[i]);
      MInt localId = get<1>(partitionLevelAncestorsOnMinLevel[i]);
      MInt id = collectMap[globalId];
      collectLevel(id) = m_minLevel;
      collectParentId(id) = -1;
      for(MInt j = 0; j < m_noDirs; j++) {
        collectNghbrIds(id, j) = a_neighborId(localId, j);
      }
      MInt minId = localMinLevelId[localId];
      ASSERT(minId > -1, "");
      collectMinLevelTreeId(id) = minLevelCellsTreeId[minId];
      collectNoOffspring(id) = (minId == noLocalMinLevelCells - 1)
                                   ? minLevelCellGlobalIdOffsets[domainId() + 1] - globalId
                                   : a_globalId(localMinLevelCells[minId + 1]) - globalId;
    }

    // Setup the truncated subtrees between the min-level partition level ancestors and the
    // partition cells (i.e. fill the collect* data buffers)
    for(MUint i = 0; i < partitionLevelAncestorsOnMinLevel.size(); i++) {
      MInt counter = collectMap[get<0>(partitionLevelAncestorsOnMinLevel[i])];
      traverseAndSetupTruncatedSubtree(counter, collectData, &collectParentId[0], &collectLevel[0],
                                       &collectNoChildren[0], &collectChildIds[0], &collectNoOffspring[0]);
    }

    // Set the child ids and the number of offsprings for the min-level partition level ancestors
    for(MUint i = 0; i < partitionLevelAncestorsOnMinLevel.size(); i++) {
      MLong globalId = get<0>(partitionLevelAncestorsOnMinLevel[i]);
      MInt localId = get<1>(partitionLevelAncestorsOnMinLevel[i]);
      MInt id = collectMap[globalId];
      a_noOffsprings(localId) = collectNoOffspring(id);
      for(MInt j = 0; j < m_maxNoChilds; j++) {
        a_childId(localId, j) = collectChildIds(id, j);
      }
    }


    MIntScratchSpace returnOffsets(collectCnt + 1, AT_, "returnOffsets");
    std::vector<MLong> returnData;

    if(collectCnt > 0) {
      returnOffsets(0) = 0;
      collectCnt = 0;
      for(MInt c = 0; c < noDomains(); c++) {
        if(noCollect[c] > 0) {
          // Create map of cells to return to this domain
          std::map<MLong, MInt> returnChain;

          for(MInt d = 0; d < noCollect[c]; d++) {
            const MInt id0 = collectOffsets[c] + d; // position in collected data buffers
            MLong globalId = collectGlobalId[id0];
            MInt id1 = collectMap[globalId];  // local position in collectData
            if(collectIsPartitionCell[id0]) { // do not send back partition cells to save some communicaton overhead
              id1 = collectMap[collectParentId[id1]];
              globalId = get<0>(collectData[id1]);
            }
            while(globalId > -1) {
              if(returnChain.count(globalId) == 0) { // cell is not already flagged for return
                returnChain[globalId] = id1;         // add to map and continue with parent cell if existing
                if(collectParentId[id1] > -1) {
                  id1 = collectMap[collectParentId[id1]];
                  globalId = get<0>(collectData[id1]);
                } else {
                  globalId = -1;
                }
              } else {
                globalId = -1;
              }
            }
          }

          for(auto it = returnChain.begin(); it != returnChain.end(); it++) {
            MLong globalId = it->first;
            MInt id1 = it->second;

            // Setup the (missing) cells on the current domain
            if(c == domainId()) {
              MInt localId = -1;
              auto it1 = m_globalToLocalId.find(globalId);
              if(it1 == m_globalToLocalId.end()) {
                // Cell is not an internal cells (not existing), append to tree
                localId = m_tree.size();
                m_tree.append();
                resetCell(localId);
                a_globalId(localId) = globalId;
                a_hasProperty(localId, Cell::IsPartLvlAncestor) = true;
                m_globalToLocalId[globalId] = localId;
              } else {
                // Cell is an internal cells and exists, perform sanity checks
                localId = it1->second;
                ASSERT(a_globalId(localId) == globalId, "");
                ASSERT(a_hasProperty(localId, Cell::IsPartLvlAncestor), "");
              }
              a_parentId(localId) = collectParentId(id1);
              a_noOffsprings(localId) = collectNoOffspring(id1);
              a_level(localId) = collectLevel(id1);

              for(MInt i = 0; i < m_maxNoChilds; i++) {
                a_childId(localId, i) = collectChildIds(id1, i);
                auto it2 = m_globalToLocalId.find(a_childId(localId, i));
                // Set the parent id if the child cell exists
                if(it2 != m_globalToLocalId.end()) {
                  a_parentId(it2->second) = globalId;
                }
              }

              newCells.push_back(localId);
              newCellsMaxLevel = mMax(newCellsMaxLevel, a_level(localId));
            } else {
              // Assemble the cell information to return to another domain
              MInt offs = returnData.size();
              returnData.resize(offs + noData);
              returnData[offs++] = globalId;
              returnData[offs++] = collectMinLevelTreeId(id1);
              returnData[offs++] = collectLevel(id1);
              returnData[offs++] = collectParentId(id1);
              returnData[offs++] = collectNoOffspring(id1);
              for(MInt i = 0; i < m_maxNoChilds; i++) {
                returnData[offs++] = collectChildIds(id1, i);
              }
              for(MInt i = 0; i < m_noDirs; i++) {
                returnData[offs++] = collectNghbrIds(id1, i);
              }
            }
          }
          if(c != domainId()) {
            returnOffsets(collectCnt + 1) = returnData.size();
            collectCnt++;
          }
        }
      }
    }

    ScratchSpace<MPI_Request> returnReq(collectCnt, AT_, "returnReq");
    ScratchSpace<MInt> noReturnDataBuffer(collectCnt, AT_, "noReturnDataBuffer");
    returnReq.fill(MPI_REQUEST_NULL);
    MInt returnCnt = 0;

    // Send the amount of cell information which will be transfered back to the other domains
    if(collectCnt > 0) {
      for(MInt c = 0; c < noDomains(); c++) {
        if(c == domainId()) continue;
        if(noCollect[c] == 0) continue;
        // Note: use persistent buffer which does not go out of scope until the send is finished
        noReturnDataBuffer[returnCnt] = returnOffsets(returnCnt + 1) - returnOffsets(returnCnt);
        MPI_Issend(&noReturnDataBuffer[returnCnt], 1, MPI_INT, c, 459, mpiComm(), &returnReq[returnCnt], AT_,
                   "&noReturnDataBuffer[returnCnt]");
        returnCnt++;
      }
    }

    MIntScratchSpace noReceiveData(std::max(queryCnt, 1), AT_, "noReceiveData");
    MIntScratchSpace receiveOffs(queryCnt + 1, AT_, "receiveOffs");
    queryCnt = 0;
    receiveOffs(0) = 0;
    // Receive the cell information counts
    for(MInt c = 0; c < noDomains(); c++) {
      if(c == domainId()) continue;
      if(noQuery[c] == 0) continue;
      MPI_Recv(&noReceiveData[queryCnt], 1, MPI_INT, c, 459, mpiComm(), MPI_STATUS_IGNORE, AT_,
               "noReceiveData[queryCnt]");
      receiveOffs(queryCnt + 1) = receiveOffs(queryCnt) + noReceiveData[queryCnt];
      queryCnt++;
    }

    if(returnCnt > 0) MPI_Waitall(returnCnt, &returnReq[0], MPI_STATUSES_IGNORE, AT_);

    MLongScratchSpace receiveData(std::max(receiveOffs(queryCnt), 1), AT_, "receiveData");

    returnReq.fill(MPI_REQUEST_NULL);
    returnCnt = 0;
    // Send the cell information and receive on corresponding domains
    if(collectCnt > 0) {
      for(MInt c = 0; c < noDomains(); c++) {
        if(c == domainId()) continue;
        if(noCollect[c] == 0) continue;
        MInt noReturnData = returnOffsets(returnCnt + 1) - returnOffsets(returnCnt);
        MInt offs = returnOffsets(returnCnt);
        MPI_Issend(&returnData[offs], noReturnData, MPI_LONG, c, 460, mpiComm(), &returnReq[returnCnt++], AT_,
                   "returnData[offs]");
      }
    }
    queryCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(c == domainId()) continue;
      if(noQuery[c] == 0) continue;
      MPI_Recv(&receiveData(receiveOffs(queryCnt)), noReceiveData(queryCnt), MPI_LONG, c, 460, mpiComm(),
               MPI_STATUS_IGNORE, AT_, "receiveData(receiveOffs(queryCnt))");
      queryCnt++;
    }

    if(returnCnt > 0) MPI_Waitall(returnCnt, &returnReq[0], MPI_STATUSES_IGNORE, AT_);

    queryCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(c == domainId()) continue;
      if(noQuery[c] == 0) continue;
      MInt offs = receiveOffs(queryCnt);
      // Loop over all received cells from this domain and add to tree or just set its information
      while(offs < receiveOffs(queryCnt + 1)) {
        MLong globalId = receiveData(offs++);
        MInt localId = -1;
        auto it = m_globalToLocalId.find(globalId);
        if(it == m_globalToLocalId.end()) {
          // Cell is not already existing, append to tree
          localId = m_tree.size();
          m_tree.append();
          resetCell(localId);
          a_globalId(localId) = globalId;
          a_hasProperty(localId, Cell::IsPartLvlAncestor) = true;
          m_globalToLocalId[globalId] = localId;
        } else {
          // Cell exists
          localId = it->second;
          ASSERT(a_globalId(localId) == globalId, "");
          ASSERT(a_hasProperty(localId, Cell::IsPartLvlAncestor), "");
        }
        MLong treeId = receiveData(offs++);
        a_level(localId) = receiveData(offs++);
        a_parentId(localId) = receiveData(offs++);
        a_noOffsprings(localId) = receiveData(offs++);

        newCells.push_back(localId);
        newCellsMaxLevel = mMax(newCellsMaxLevel, a_level(localId));

        for(MInt i = 0; i < m_maxNoChilds; i++) {
          a_childId(localId, i) = receiveData(offs++);
          auto it1 = m_globalToLocalId.find(a_childId(localId, i));
          // Set the parent id if the child cell exists
          if(it1 != m_globalToLocalId.end()) {
            a_parentId(it1->second) = globalId;
          }
        }
        for(MInt i = 0; i < m_noDirs; i++) {
          a_neighborId(localId, i) = receiveData(offs++);
        }
        // Compute the cell coordinates of the min-level cells
        if(a_level(localId) == m_minLevel) {
          maia::grid::hilbert::treeIdToCoordinates<nDim>(
              &a_coordinate(localId, 0), treeId, (MLong)grid_.m_targetGridMinLevel,
              &grid_.m_targetGridCenterOfGravity[0], grid_.m_targetGridLengthLevel0);
        }
      }
      queryCnt++;
    }

    m_noHaloPartitionLevelAncestors = m_tree.size() - m_noInternalCells;

    for(MInt i = 0; i < m_noPartitionCells; i++) {
      MLong globalId = m_localPartitionCellGlobalIds[i];
      MInt cellId = globalId - m_domainOffsets[domainId()];
      if(localMinLevelId[cellId] < 0) { // partition cell is not on min-level
        MInt parentId = m_globalToLocalId[a_parentId(cellId)];
        a_level(cellId) = a_level(parentId) + 1; // set the level of the partition cell
      }
    }

    static const MFloat childStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                              {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
    // Loop up starting from the min-level and compute the cell coordinates of the child cells
    for(MInt lvl = m_minLevel; lvl <= newCellsMaxLevel; lvl++) {
      for(MUint i = 0; i < newCells.size(); i++) {
        if(a_level(newCells[i]) == lvl) {
          for(MInt child = 0; child < m_maxNoChilds; child++) {
            if(a_childId(newCells[i], child) > -1) {
              auto it = m_globalToLocalId.find(a_childId(newCells[i], child));
              if(it != m_globalToLocalId.end()) {
                for(MInt j = 0; j < nDim; ++j) {
                  a_coordinate(it->second, j) =
                      a_coordinate(newCells[i], j) + childStencil[child][j] * F1B2 * cellLengthAtLevel(lvl + 1);
                }
              }
            }
          }
        }
      }
    }

    for(MInt cellId = m_noInternalCells; cellId < m_tree.size(); cellId++) {
      a_hasProperty(cellId, Cell::IsHalo) = true;
    }

    for(MUint i = 0; i < newCells.size(); i++) {
      MInt cellId = newCells[i];
      ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor), "");
      // Store partition level ancestor cell ids and their child ids
      m_partitionLevelAncestorIds.push_back(cellId);
      for(MInt child = 0; child < m_maxNoChilds; child++) {
        m_partitionLevelAncestorChildIds.push_back(a_childId(cellId, child));
      }

      MInt ndom = findNeighborDomainId(a_globalId(cellId));
      if(ndom == domainId()) {
        MLong lastId = a_globalId(cellId) + a_noOffsprings(cellId) - 1;
        MInt firstDom = domainId() + 1;
        MInt lastDom = findNeighborDomainId(lastId);
        for(MInt d = firstDom; d <= lastDom; d++) {
          MInt idx =
              findIndex(m_partitionLevelAncestorNghbrDomains.begin(), m_partitionLevelAncestorNghbrDomains.end(), d);
          if(idx < 0) {
            idx = m_partitionLevelAncestorNghbrDomains.size();
            m_partitionLevelAncestorNghbrDomains.push_back(d);
            m_partitionLevelAncestorHaloCells.push_back(std::vector<MInt>());
            m_partitionLevelAncestorWindowCells.push_back(std::vector<MInt>());
          }
        }
      } else {
        MInt idx =
            findIndex(m_partitionLevelAncestorNghbrDomains.begin(), m_partitionLevelAncestorNghbrDomains.end(), ndom);
        if(idx < 0) {
          idx = (MInt)m_partitionLevelAncestorNghbrDomains.size();
          m_partitionLevelAncestorNghbrDomains.push_back(ndom);
          m_partitionLevelAncestorHaloCells.push_back(std::vector<MInt>());
          m_partitionLevelAncestorWindowCells.push_back(std::vector<MInt>());
        }
      }
    }

    sort(m_partitionLevelAncestorNghbrDomains.begin(),
         m_partitionLevelAncestorNghbrDomains
             .end()); // sort by ascending neighbor domain id to avoid send/recv deadlocks

    for(MUint i = 0; i < newCells.size(); i++) {
      MInt cellId = newCells[i];
      ASSERT(a_hasProperty(cellId, Cell::IsPartLvlAncestor), "");
      MInt ndom = findNeighborDomainId(a_globalId(cellId));
      if(ndom == domainId()) {
        MLong lastId = a_globalId(cellId) + a_noOffsprings(cellId) - 1;
        MInt firstDom = domainId() + 1;
        MInt lastDom = findNeighborDomainId(lastId);
        for(MInt d = firstDom; d <= lastDom; d++) {
          MInt idx =
              findIndex(m_partitionLevelAncestorNghbrDomains.begin(), m_partitionLevelAncestorNghbrDomains.end(), d);
          ASSERT(idx > -1, "");
          m_partitionLevelAncestorWindowCells[idx].push_back(cellId);
        }
      } else {
        MInt idx =
            findIndex(m_partitionLevelAncestorNghbrDomains.begin(), m_partitionLevelAncestorNghbrDomains.end(), ndom);
        ASSERT(idx > -1, "");
        m_partitionLevelAncestorHaloCells[idx].push_back(cellId);
      }
    }

    // note that window/halo cells are already sorted by ascending globalIds and ascending neighborDomainIds!
    for(MUint i = 1; i < m_partitionLevelAncestorNghbrDomains.size(); i++) {
      if(m_partitionLevelAncestorNghbrDomains[i] <= m_partitionLevelAncestorNghbrDomains[i - 1]) {
        mTerm(1, AT_, "Invalid sorting of neighbor domains.");
      }
    }
    for(MUint i = 0; i < m_partitionLevelAncestorNghbrDomains.size(); i++) {
      for(MUint j = 1; j < m_partitionLevelAncestorHaloCells[i].size(); j++) {
        if(m_partitionLevelAncestorHaloCells[i][j] <= m_partitionLevelAncestorHaloCells[i][j - 1]) {
          mTerm(1, AT_, "Invalid sorting of halo cells.");
        }
      }
      for(MUint j = 1; j < m_partitionLevelAncestorWindowCells[i].size(); j++) {
        if(m_partitionLevelAncestorWindowCells[i][j] <= m_partitionLevelAncestorWindowCells[i][j - 1]) {
          mTerm(1, AT_, "Invalid sorting of window cells.");
        }
      }
    }
  }


  //-----------------------------------------------------------------------


  /** \brief Generic communication of data given by dataSource, matched by queryGlobalIds, stored in receiveData
      \author Lennart Schneiders
      \date October 2017
    */
  void queryGlobalData(const MInt noCellsToQuery, const MLong* const queryGlobalIds, MLong* const receiveData,
                       const MLong* const dataSource, const MInt dataWidth) {
    TRACE();

    MIntScratchSpace noRecv(noDomains(), AT_, "noRecv");
    MIntScratchSpace noSend(noDomains(), AT_, "noSend");
    MIntScratchSpace recvOffsets(noDomains() + 1, AT_, "recvOffsets");
    MIntScratchSpace sendOffsets(noDomains() + 1, AT_, "sendOffsets");
    MIntScratchSpace queryDomains(noCellsToQuery, AT_, "queryDomains");
    MIntScratchSpace originalId(noCellsToQuery, AT_, "originalId");
    noRecv.fill(0);
    queryDomains.fill(-1);
    originalId.fill(-1);

    for(MInt i = 0; i < noCellsToQuery; i++) {
      MLong cellId = queryGlobalIds[i];
      MInt cpu = -1;
      for(MInt c = 0; c < noDomains(); c++) {
        if(cellId >= m_domainOffsets[c] && cellId < m_domainOffsets[c + 1]) {
          cpu = c;
          c = noDomains();
        }
      }
      // if ( cpu < 0 || cpu == domainId() ) {
      if(cpu < 0) {
        mTerm(1, AT_, "Neighbor domain not found.");
      }
      queryDomains[i] = cpu;
      noRecv[cpu]++;
    }

    recvOffsets(0) = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      recvOffsets(c + 1) = recvOffsets(c) + noRecv[c];
    }
    MLongScratchSpace recvIds(recvOffsets[noDomains()], AT_, "recvIds");
    noRecv.fill(0);

    for(MInt i = 0; i < noCellsToQuery; i++) {
      MLong cellId = queryGlobalIds[i];
      MInt cpu = queryDomains[i];
      recvIds(recvOffsets(cpu) + noRecv(cpu)) = cellId;
      originalId(recvOffsets(cpu) + noRecv(cpu)) = i;
      noRecv[cpu]++;
    }
    MPI_Alltoall(&noRecv[0], 1, MPI_INT, &noSend[0], 1, MPI_INT, mpiComm(), AT_, "noRecv[0]", "noSend[0]");
    sendOffsets(0) = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      sendOffsets(c + 1) = sendOffsets(c) + noSend[c];
    }

    MLongScratchSpace recvData(recvOffsets[noDomains()], dataWidth, AT_, "recvChildIds");
    MLongScratchSpace sendData(sendOffsets[noDomains()], dataWidth, AT_, "sendChildIds");
    MLongScratchSpace sendIds(sendOffsets[noDomains()], AT_, "sendIds");
    ScratchSpace<MPI_Request> sendReq(noDomains(), AT_, "sendReq");

    sendReq.fill(MPI_REQUEST_NULL);
    MInt recvCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(noRecv[c] == 0) continue;
      MPI_Issend(&recvIds[recvOffsets[c]], noRecv[c], MPI_LONG, c, 123, mpiComm(), &sendReq[recvCnt], AT_,
                 "recvIds[recvOffsets[c]]");
      recvCnt++;
    }
    for(MInt c = 0; c < noDomains(); c++) {
      if(noSend[c] == 0) continue;
      MPI_Recv(&sendIds[sendOffsets[c]], noSend[c], MPI_LONG, c, 123, mpiComm(), MPI_STATUS_IGNORE, AT_,
               "sendIds[sendOffsets[c]]");
    }
    if(recvCnt > 0) MPI_Waitall(recvCnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);


    for(MInt k = 0; k < sendOffsets[noDomains()]; k++) {
      if(sendIds[k] < m_domainOffsets[domainId()] || sendIds[k] >= m_domainOffsets[domainId() + 1])
        mTerm(1, AT_, "Invalid exchange.");
      auto it = m_globalToLocalId.find(sendIds[k]);
      if(it == m_globalToLocalId.end()) mTerm(1, AT_, "Invalid exchange2.");
      MInt cellId = it->second;
      for(MInt c = 0; c < dataWidth; c++) {
        sendData(k, c) = dataSource[cellId * dataWidth + c];
      }
    }

    sendReq.fill(MPI_REQUEST_NULL);
    MInt sendCnt = 0;
    for(MInt c = 0; c < noDomains(); c++) {
      if(noSend[c] == 0) continue;
      MPI_Issend(&sendData(sendOffsets[c], 0), noSend[c] * dataWidth, MPI_LONG, c, 124, mpiComm(), &sendReq[sendCnt],
                 AT_, "sendData(sendOffsets[c]");
      sendCnt++;
    }
    for(MInt c = 0; c < noDomains(); c++) {
      if(noRecv[c] == 0) continue;
      MPI_Recv(&recvData(recvOffsets[c], 0), noRecv[c] * dataWidth, MPI_LONG, c, 124, mpiComm(), MPI_STATUS_IGNORE, AT_,
               "recvData(recvOffsets[c]");
    }
    if(sendCnt > 0) MPI_Waitall(sendCnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

    for(MInt k = 0; k < recvOffsets[noDomains()]; k++) {
      MInt cellId = originalId[k];
      for(MInt c = 0; c < dataWidth; c++) {
        receiveData[cellId * dataWidth + c] = recvData[k * dataWidth + c];
      }
    }
  }


  //-----------------------------------------------------------------------


  /** \brief Recursively setup subtree connection from partition level to the leaf level
      \author Lennart Schneiders
      \date October 2017
  */
  void traverseAndSetupSubtree(const MInt& cellId, const MUchar* const cellInfo) {
    static constexpr MFloat childStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                                  {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
    const MUint childCnt =
        static_cast<MUint>(cellInfo[cellId]) & 15u; // bitwise operations on char evoke integer promotion
    a_noOffsprings(cellId) = 1;
    std::fill(&a_childId(cellId, 0), &a_childId(cellId, 0) + m_maxNoChilds, -1);
    //---
    for(MUint child = 0; child < childCnt; child++) {
      MInt childId = cellId + a_noOffsprings(cellId);
      ASSERT(childId < m_tree.size(),
             "child out of range: " + std::to_string(childId) + " " + std::to_string(m_tree.size()));
      MUint position = (static_cast<MUint>(cellInfo[childId]) >> 4u) & 7u;
      a_level(childId) = a_level(cellId) + 1;
      a_childId(cellId, (MInt)position) = a_globalId(cellId) + a_noOffsprings(cellId);
      a_parentId(childId) = a_globalId(cellId);
      for(MInt j = 0; j < nDim; ++j) {
        a_coordinate(childId, j) =
            a_coordinate(cellId, j) + childStencil[(MUint)position][j] * F1B2 * cellLengthAtCell(childId);
      }
      traverseAndSetupSubtree(childId, cellInfo);
      a_noOffsprings(cellId) += a_noOffsprings(childId);
    }
  }


  //-----------------------------------------------------------------------


  /** \brief Recursively flag subtree cells from partition level to the leaf level
      \author Lennart Schneiders
      \date October 2017
    */
  MInt traverseAndFlagSubtree(const MUint& cellId, const MUchar* const cellInfo, const MInt N, MInt* const flag) {
    const MUint childCnt =
        static_cast<MUint>(cellInfo[cellId]) & 15u; // bitwise operations on char evoke integer promotion
    flag[cellId] = 0;
    MInt noFlagged = 1;
    MUint childId = cellId;
    MUint totalNoChilds = childCnt;
    //---
    while(totalNoChilds) {
      childId++;
      ASSERT((MInt)childId < N, "child out of range: " + std::to_string(childId) + " " + std::to_string(N));
      MUint isMinLevel = (static_cast<MUint>(cellInfo[childId]) >> 7) & 1;
      ASSERT(!isMinLevel, "");
      flag[childId] = 0;
      noFlagged++;
      totalNoChilds--;
      totalNoChilds += static_cast<MUint>(cellInfo[childId]) & 15u;
    }
    return noFlagged;
  }


  //-----------------------------------------------------------------------


  /** \brief Recursively setup truncated subtree between minLevel and partitionLevel (used to setup partition level
     ancestors) \author Lennart Schneiders \date October 2017
  */
  void traverseAndSetupTruncatedSubtree(MInt& counter, const std::vector<std::tuple<MLong, MUchar, MInt>>& collectData,
                                        MLong* const collectParentId, MInt* const collectLevel,
                                        MInt* const collectNoChildren, MLong* const collectChildIds,
                                        MInt* const collectNoOffspring) {
    using namespace std;

    ASSERT(counter < (signed)collectData.size(), "");
    const MLong globalId = get<0>(collectData[counter]);
    const MInt counterParent = counter;
    const MInt isPartitionCell = get<2>(collectData[counterParent]);
    const MUint childCnt = static_cast<MUint>(get<1>(collectData[counterParent]))
                           & 15u; // bitwise operations on char evoke integer promotion
    if(isPartitionCell) return;
    std::fill(&collectChildIds[m_maxNoChilds * counterParent],
              &collectChildIds[m_maxNoChilds * counterParent] + m_maxNoChilds, -1);
    //---
    for(MUint child = 0; child < childCnt; child++) {
      counter++;
      ASSERT(counter < (signed)collectData.size(), "");
      const MLong childId = get<0>(collectData[counter]);
      const MUint position = (static_cast<MUint>(get<1>(collectData[counter])) >> 4u) & 7u;
      collectLevel[counter] = collectLevel[counterParent] + 1;
      collectChildIds[m_maxNoChilds * counterParent + (MInt)position] = childId;
      collectParentId[counter] = globalId;
      MLong nextId = -1;
      if(child < childCnt - 1) {
        MUint totalNoChilds = 1;
        MInt tmpCnt = 0;
        while(totalNoChilds) {
          ASSERT(counter + tmpCnt < (signed)collectData.size(), "");
          if(!get<2>(collectData[counter + tmpCnt])) {
            totalNoChilds += static_cast<MUint>(get<1>(collectData[counter + tmpCnt])) & 15u;
          }
          tmpCnt++;
          totalNoChilds--;
          ASSERT(counter + tmpCnt < (signed)collectData.size(), "");
          nextId = get<0>(collectData[counter + tmpCnt]);
        }
      } else {
        nextId = globalId + collectNoOffspring[counterParent];
      }
      collectNoOffspring[counter] = nextId - childId;
      traverseAndSetupTruncatedSubtree(counter, collectData, collectParentId, collectLevel, collectNoChildren,
                                       collectChildIds, collectNoOffspring);
    }
    collectNoChildren[counterParent] = childCnt;
  }


  //-----------------------------------------------------------------------


  /** \brief Save the grid
      \author
      \date October 2017
  */
  template <class CELLTYPE>
  void saveGrid(const MChar* fileName, Collector<CELLTYPE>* cpu_cells,
                const std::vector<std::vector<MInt>>& cpu_haloCells,
                const std::vector<std::vector<MInt>>& cpu_windowCells,
                const std::vector<std::vector<MInt>>& cpu_azimuthalHaloCells,
                const std::vector<MInt>& cpu_azimuthalUnmappedHaloCells,
                const std::vector<std::vector<MInt>>& cpu_azimuthalWindowCells, MInt* const recalcIdTree) {
    grid_.setLevel();

    const MInt cpu_noCells = grid_.m_tree.size();

    MLong tmpCnt = 0; // Variable used for various counting activities.

    /* Mark all halo cells and communicate the no. of real cells for each domain
     * -------------------------------------------------------------------------
     *
     * MInt noHaloCellsTotal == The number of halo cells this domain has.
     * MIntScratchSpace noCellsPar[cpuCnt] == Contains for each domain the number of cell without halos.
     */

    std::vector<std::vector<MInt>> partitionLevelAncestorHaloCells;
    std::vector<std::vector<MInt>> partitionLevelAncestorWindowCells;
    partitionLevelAncestorHaloCells.resize(noNeighborDomains());
    partitionLevelAncestorWindowCells.resize(noNeighborDomains());

    MIntScratchSpace noCellsPar(noDomains(), AT_, "noCellsPar");
    MIntScratchSpace offspringPrefix(noDomains(), AT_, "offspringPrefix");
    MInt noHaloCellsTotal = 0;
    MInt noWindowCellsTotal = 0;
    MInt noHalo0 = cpu_noCells - m_noInternalCells;
    if(noNeighborDomains() > 0) {
      for(MInt i = 0; i < noNeighborDomains(); ++i) {
        noHaloCellsTotal += (signed)cpu_haloCells[i].size();
        noWindowCellsTotal += (signed)cpu_windowCells[i].size();
        for(MInt j = 0; j < (signed)cpu_haloCells[i].size(); ++j) {
          if(a_hasProperty(cpu_haloCells[i][j], Cell::IsPartLvlAncestor)) {
            partitionLevelAncestorHaloCells[i].push_back(cpu_haloCells[i][j]);
          }
        }
        for(MInt j = 0; j < (signed)cpu_windowCells[i].size(); ++j) {
          if(a_hasProperty(cpu_windowCells[i][j], Cell::IsPartLvlAncestor)) {
            partitionLevelAncestorWindowCells[i].push_back(cpu_windowCells[i][j]);
          }
        }
      }
    }
    // Azimuthal periodicity
    if(grid_.m_azimuthalPer) {
      if(noAzimuthalNeighborDomains() > 0) {
        for(MInt i = 0; i < noAzimuthalNeighborDomains(); ++i) {
          noHaloCellsTotal += (signed)cpu_azimuthalHaloCells[i].size();
          noWindowCellsTotal += (signed)cpu_azimuthalWindowCells[i].size();
        }
        noHaloCellsTotal += (signed)cpu_azimuthalUnmappedHaloCells.size();
      }
    }
    ASSERT(noHalo0 == noHaloCellsTotal, std::to_string(noHalo0) + " != " + std::to_string(noHaloCellsTotal));

    // Note: the following code was originally used. However, since right now only the FV solver uses
    // this method, it was made the default behavior.
    // // only valid for gcell collector !!!
    // noCellsPar.p[domainId()] = cpu_noCells - noHaloCellsTotal;
    // MInt noInternalCells = cpu_noCells;

    // if (std::is_same<CELLTYPE, FvCell>::value) {
    //   noCellsPar.p[domainId()] = m_noInternalCells ;
    //   noInternalCells = m_noInternalCells;
    // }

    // NOTE: for minLevelIncrease skip all cells below the m_newMinLevel!
    //      the default is -1
    MInt changeSize = grid_.m_newMinLevel > 0 ? m_noInternalCells : 1;
    MIntScratchSpace changeId2(changeSize, AT_, "changeId2");
    changeId2.fill(-1);

    MInt noInternalCells = m_noInternalCells;
    if(grid_.m_newMinLevel > 0) { // re-count noInternalCells
      noInternalCells = 0;
      for(MInt cellId = 0; cellId < cpu_noCells; cellId++) {
        if(a_level(cellId) < grid_.m_newMinLevel) continue;
        if(a_hasProperty(cpu_cells, cellId, Cell::IsHalo)) continue;
        changeId2[cellId] = noInternalCells;
        noInternalCells++;
      }
    }
    noCellsPar.p[domainId()] = noInternalCells;

    if(grid_.m_newMinLevel < 0) {
      ASSERT(!a_hasProperty(cpu_cells, m_noInternalCells - 1, Cell::IsHalo), "");
      if(noNeighborDomains() > 0) {
        ASSERT(a_hasProperty(cpu_cells, m_noInternalCells, Cell::IsHalo), "");
      }
    }

    if(1 < noDomains()) {
      tmpCnt = noCellsPar.p[domainId()];
      MInt sndBuf = static_cast<MInt>(tmpCnt);
      MInt* rcvBuf = noCellsPar.getPointer();
      MPI_Allgather(&sndBuf, 1, MPI_INT, rcvBuf, 1, MPI_INT, mpiComm(), AT_, "sndBuf", "rcvBuf");
      tmpCnt = sndBuf;
    }

    // Calculate noOffsprings and workload
    calculateNoOffspringsAndWorkload(cpu_cells, cpu_noCells, cpu_haloCells, cpu_windowCells,
                                     partitionLevelAncestorWindowCells, partitionLevelAncestorHaloCells);

    // Create minLevelCells std::vector of all most coarse cells
    std::vector<MInt> minLevelCells;
    for(MInt i = 0; i < m_noInternalCells; ++i) {
      if(a_hasProperty(cpu_cells, i, Cell::IsHalo)) { // If we have a halo cell, continue with the next cell.
        continue;
      }
      // increase the minLevel by one when writing the restartGrid
      if(a_level(i) == grid_.m_newMinLevel) {
        minLevelCells.push_back(i);
      }
      if(grid_.m_newMinLevel > 0) continue;

      if(a_parentId(i) == -1) {
        ASSERT(a_level(i) == m_minLevel, std::to_string(a_level(i)) + " " + std::to_string(m_minLevel));
        minLevelCells.push_back(i);
      }
    }

    // Calculating hilbert id
    const MInt minLevelCellsCnt = minLevelCells.size();
    const MInt minLevel = grid_.m_newMinLevel > 0 ? grid_.m_newMinLevel : m_minLevel;


    MIntScratchSpace minLevelCells_id(mMax(1, minLevelCellsCnt), AT_, "minLevelCells_id");
    MLongScratchSpace minLevelCells_hId(mMax(1, minLevelCellsCnt), AT_, "minLevelCells_hId");
    minLevelCells_hId.fill(0); //!< reset the hilbert ids before calling hilbertIndex !!!
    for(MInt i = 0; i < minLevelCellsCnt; ++i) {
      minLevelCells_id[i] = minLevelCells[i];
      minLevelCells_hId[i] = generateHilbertIndex(minLevelCells[i], minLevel);
    }

    // sorting minLevel cells after hilbert id
    sortAfterHilbertIds(minLevelCells_hId.getPointer(), minLevelCells_id.getPointer(), minLevelCellsCnt);

    // count offSprings for minLevel cells!
    MInt offspringSum = 0;
    for(MInt i = 0; i < minLevelCellsCnt; ++i) {
      ASSERT(a_noOffsprings(cpu_cells, minLevelCells_id.p[i]) > -1, "minLevelCells_id.p[i]) " << minLevelCells_id.p[i]);
      offspringSum += a_noOffsprings(cpu_cells, minLevelCells_id.p[i]);
    }
    MPI_Allgather(&offspringSum, 1, MPI_INT, offspringPrefix.getPointer(), 1, MPI_INT, mpiComm(), AT_, "offspringSum",
                  "offspringPrefix.getPointer()");


    // Calculating Offset
    tmpCnt = m_32BitOffset;
    for(MInt i = 0; i < domainId(); ++i) {
      // tmpCnt += noCellsPar.p[i];
      tmpCnt += offspringPrefix.p[i];
    }


    // Calculating the new id (sorted grid-cellId) for all minCells
    MLongScratchSpace minLevelCells_newId(mMax(1, minLevelCellsCnt), AT_, "minLevelCells_newId");
    for(MInt i = 0; i < minLevelCellsCnt; ++i) {
      minLevelCells_newId.p[i] = tmpCnt;
      tmpCnt += a_noOffsprings(cpu_cells, minLevelCells_id.p[i]);
    }


    if(grid_.m_newMinLevel < 0) {
      for(MInt i = 0; i < minLevelCellsCnt; ++i) {
        ASSERT(minLevelCells_newId.p[i] == a_globalId(minLevelCells_id.p[i]),
               "minLevelCells_newId.p[i] " << minLevelCells_newId.p[i] << "a_globalId(minLevelCells_id.p[i] "
                                           << a_globalId(minLevelCells_id.p[i]));
      }
    }

    // Creating new order of ids for own cells

    // changeId from unsorted to sorted grid cells after globalId!
    MLongScratchSpace changeId(cpu_noCells, AT_, "changeId");
    // old unsorted ordering
    MIntScratchSpace oldId(noCellsPar[domainId()], AT_, "oldId");
    // new sorted ordering
    MLongScratchSpace newId(noCellsPar[domainId()], AT_, "newId");

    tmpCnt = -1;
    tmpCnt = createTreeOrderingOfIds(cpu_cells, m_noInternalCells, oldId.getPointer(), newId.getPointer(),
                                     minLevelCells_id.getPointer(), minLevelCells_newId.getPointer(), minLevelCellsCnt,
                                     changeId.getPointer(), partitionLevelAncestorWindowCells,
                                     partitionLevelAncestorHaloCells);

    ASSERT(tmpCnt == noCellsPar[domainId()], "tmpCnt: " << tmpCnt << " noCellsPar: " << noCellsPar[domainId()]);

    // Fill this so that it can be used in the solvers to write their restart-File in the
    // same ordering!
    if(recalcIdTree) {
      for(MInt i = 0; i < tmpCnt; i++) {
        ASSERT(i < tmpCnt, "index out of range");
        recalcIdTree[i] = oldId[i];
        ASSERT(recalcIdTree[i] > -1 && recalcIdTree[i] < m_noInternalCells, "recalcId out of range");
      }
    }

    // Create the changeId array in wich: changeId.p[oldId] = newId.
    // changeId is -2 for halo-cells!
    for(MInt i = 0; i < cpu_noCells; ++i) {
      changeId[i] = -2;
    }

    for(MInt i = 0; i < noCellsPar[domainId()]; ++i) {
      changeId[oldId[i]] = newId[i];
    }

    // Communicate new id of the halo cells.
    if(1 < noDomains() && (noWindowCellsTotal > 0) && (noHaloCellsTotal > 0)) {
      // Creating ScratchSpace for the communication
      MLongScratchSpace sndBuf(noWindowCellsTotal, AT_, "sndBuf");
      MLongScratchSpace rcvBuf(noHaloCellsTotal, AT_, "rcvBuf");
      MIntScratchSpace sndDsp(noNeighborDomains(), AT_, "sndDsp");
      MIntScratchSpace rcvDsp(noNeighborDomains(), AT_, "rcvDsp");
      MIntScratchSpace sndCnt(noNeighborDomains(), AT_, "sndCnt");
      MIntScratchSpace rcvCnt(noNeighborDomains(), AT_, "rcvCnt");

      tmpCnt = 0;
      for(MInt i = 0; i < noNeighborDomains(); ++i) {
        // Fill the send buffer with the new ids of the own window cells.
        for(MInt j = 0; j < (signed)cpu_windowCells[i].size(); ++j) {
          sndBuf.p[tmpCnt] = changeId.p[cpu_windowCells[i][j]];
          ++tmpCnt;
        }

        // Calculate the offsets for sending and receiving datas.
        if(0 == i) {
          sndDsp.p[0] = 0;
          rcvDsp.p[0] = 0;
        } else {
          sndDsp.p[i] = sndDsp.p[i - 1] + ((signed)cpu_windowCells[i - 1].size());
          rcvDsp.p[i] = rcvDsp.p[i - 1] + ((signed)cpu_haloCells[i - 1].size());
        }

        sndCnt.p[i] = (signed)cpu_windowCells[i].size();
        rcvCnt.p[i] = (signed)cpu_haloCells[i].size();
      }

      // Note: it might be required to use a unique mpi-tag for the communication below to avoid any
      // conflicting messages and possibly deadlocks in case there are still other open MPI requests
      const MInt mpiTag = 0;
      if(0 < noNeighborDomains()) {
        MPI_Request* mpiReq;
        mpiReq = new MPI_Request[noNeighborDomains()];
        MInt cntS = 0;
        for(MInt i = 0; i < noNeighborDomains(); i++) {
          mpiReq[i] = MPI_REQUEST_NULL;
          MInt bufSize = (signed)cpu_windowCells[i].size();
          if(bufSize == 0) continue;
          MPI_Issend(&(sndBuf[cntS]), bufSize, MPI_LONG, m_nghbrDomains[i], mpiTag, mpiComm(), &mpiReq[i], AT_,
                     "(sndBuf[cntS])");
          cntS += bufSize;
        }
        MPI_Status status;
        MInt cntR = 0;
        for(MInt i = 0; i < noNeighborDomains(); i++) {
          MInt bufSize = (signed)cpu_haloCells[i].size();
          if(bufSize == 0) continue;
          MPI_Recv(&(rcvBuf[cntR]), bufSize, MPI_LONG, m_nghbrDomains[i], mpiTag, mpiComm(), &status, AT_,
                   "(rcvBuf[cntR])");
          cntR += bufSize;
        }
        for(MInt i = 0; i < noNeighborDomains(); i++) {
          if((signed)cpu_windowCells[i].size() == 0) continue;
          MPI_Wait(&mpiReq[i], &status, AT_);
        }
      }


      // Use the new datas of the receive buffer to complete the missing ids in changeId.
      tmpCnt = 0;
      for(MInt i = 0; i < noNeighborDomains(); ++i) {
        for(MInt j = 0; j < (signed)cpu_haloCells[i].size(); ++j) {
          changeId.p[cpu_haloCells[i][j]] = rcvBuf.p[tmpCnt];
          ++tmpCnt;
        }
      }
    }

    /* Filtering partitionCells using m_partitionCellOffspringThreshold
     * ------------------------------------------
     * To ensure that later we have an 'optimal' distribution of cells on all domains, we now filters out cells that are
     * too big (too many offsprings) or too heavy (too much workload) using a threshold value defined through the
     * property partitionCellOffspringThreshold.
     *
     *  Choosing a right number for partitionCellOffspringThreshold is a user choice... the lower the value, the more
     * domains can be used for calcution but also more cells gets filtered and the amount of partitionCells used for the
     * distribution increases, which can result in increases in time while loading the grid and a bigger grid file.
     * Choosing a value too high will reduce the amount of domains usable for calculation on the grid but could increase
     * the speed at which the grid gets loaded.
     *
     * 1. Calculate the maximum number of all cells our grid has and the threshold for workload and noOffsprings.
     * 2. Set a vector containing all minLevelCells so far.
     * 3. First filter step: noOffsprings
     *   - If a cell has more offsprings than the offspringsThreshold or a cell has more Offspring than the max number
     * of cells allowed on a single domain. Delete it and replace it through her childrens.
     *   - Re-filter each child.
     *   - Continue with next cell.
     * 4. Second filter step: workload
     *   - If a cell has more workload than the workload threshold. Delete it and replace it through her childrens.
     *   - Re-filter each child.
     *   - Continue with next cell.
     * 5. Compute the level difference each cell has to the minLevelCells from the beginning. This value is used in the
     * distribution to check if we have a real partitionCell or a splitted one.
     * 6. Communicate partitionCellsFiltered.size() to all domains and calculate the number of all filtered cells,
     * needed later for offset calculation on file writing.
     *
     * ZfSInt noCellsParMax == addition of all cells on all domains.
     * MInt offSpringThreshold == How much offsprings a cell can have before becoming filtered and splitted into her
     * childrens. MFloat workloadThreshold == How much weight a cell can have before becoming filtered and splitted
     * into her childrens. MInt maxCells == maximum number of cells a single domains can hold. iNTscratchSpace
     * partitionCellLevelDiff == holds the level difference between the cells and the partitionCells. MIntScratchSpace
     * partitionCellsFilteredSizePar == holds the amount of filtered cells for each domains. MInt
     * partitionCellsFilteredSizeParMax == addition of all partitionCellsFilteredSize on all domains.
     *
     */

    // Calculating Offset
    tmpCnt = m_32BitOffset;
    for(MInt i = 0; i < domainId(); ++i) {
      tmpCnt += noCellsPar.p[i];
    }

    for(MInt i = 0; i < noCellsPar[domainId()]; ++i) {
      ASSERT(changeId[oldId[i]] >= tmpCnt && changeId[oldId[i]] < tmpCnt + noCellsPar.p[domainId()], "");
    }

    // Setting thresholds
    const MInt offspringsThreshold = m_partitionCellOffspringThreshold;
    const MFloat workloadThreshold = m_partitionCellWorkloadThreshold;

    std::vector<MLong> partitionCellsFiltered;

    if(!grid_.m_updatedPartitionCells && grid_.m_updatePartitionCellsOnRestart) {
      std::vector<MInt> partitionCellsGlobalIds;

      // 2. Vector containing all minLevelCells
      for(MInt i = 0; i < minLevelCellsCnt; ++i) {
        ASSERT(minLevelCells_newId.p[i] >= tmpCnt && minLevelCells_newId.p[i] < tmpCnt + noCellsPar.p[domainId()],
               "minLevelCells_newId:" << minLevelCells_newId.p[i] << " tmpCnt: " << tmpCnt
                                      << " noCellsPar: " << noCellsPar.p[domainId()]);
        ASSERT(minLevelCells_newId.p[i] == changeId[oldId.p[minLevelCells_newId.p[i] - tmpCnt]], "");
      }

      // 3.Filter steps: noOffsprings/workload
      std::vector<MInt> tmpPartitionCells;
      for(MInt cellId = 0; cellId < cpu_noCells; ++cellId) {
        // if (a_parentId(cpu_cells, cellId) == -1) {
        if(a_level(cellId) == minLevel && !a_hasProperty(cpu_cells, cellId, Cell::IsPeriodic)) {
          //   ASSERT( a_level(cpu_cells, cellId) == m_minLevel, "" );
          if(!grid_.m_newMinLevel) {
            ASSERT(a_parentId(cellId) == -1, "");
          }
          tmpPartitionCells.push_back(cellId);
        }
      }
      // TODO labels:GRID,ADAPTATION @ansgar add multisolver check for leaf cells!
      while(!tmpPartitionCells.empty()) {
        MInt cellId = tmpPartitionCells.front();
        tmpPartitionCells.erase(tmpPartitionCells.begin());
        if(a_noOffsprings(cpu_cells, cellId) > offspringsThreshold
           || a_workload(cpu_cells, cellId) > workloadThreshold) {
          MInt cnt = 0;
          for(MUint j = 0; j < IPOW2(nDim); ++j) {
            if(a_childId(cellId, j) > -1) {
              tmpPartitionCells.insert(tmpPartitionCells.begin() + cnt, a_childId(cellId, j));
              cnt++;
            }
          }
        } else if(!a_hasProperty(cpu_cells, cellId, Cell::IsHalo)) {
          partitionCellsFiltered.push_back(changeId.p[cellId]);
          partitionCellsGlobalIds.push_back(grid_.a_globalId(cellId));
        }
      }

      // update partition cell offsets and partition-cell localIds which are required when writing LPT-restart files!
      MLong noPartitionCellsGlobal = partitionCellsFiltered.size();
      MPI_Allreduce(MPI_IN_PLACE, &noPartitionCellsGlobal, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_,
                    "MPI_IN_PLACE", "m_noPartitionCellsGlobal");
      if(noPartitionCellsGlobal != grid_.m_noPartitionCellsGlobal) {
        cerr0 << "Partition cell update in restart-file! " << noPartitionCellsGlobal << " "
              << grid_.m_noPartitionCellsGlobal << std::endl;
      }

      MInt noLocalPartitionCells = partitionCellsFiltered.size();

      if(grid_.m_localPartitionCellLocalIdsRestart != nullptr) {
        mDeallocate(grid_.m_localPartitionCellLocalIdsRestart);
      }
      mAlloc(grid_.m_localPartitionCellLocalIdsRestart, noLocalPartitionCells, "m_localPartitionCellLocalIdsRestart",
             -1, AT_);

      // sort by globalIds
      sort(partitionCellsGlobalIds.begin(), partitionCellsGlobalIds.end());

      MIntScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
      // determine offset within GLOBAL!!! partition cells
      MPI_Allgather(&noLocalPartitionCells, 1, MPI_INT, &localPartitionCellCounts[0], 1, MPI_INT, mpiComm(), AT_,
                    "noLocalPartitionCells", "localPartitionCellCounts[0]");

      // sum up all previous domains partition cells
      MInt offset = 0;
      for(MInt dId = 0; dId < domainId(); dId++) {
        offset += localPartitionCellCounts[dId];
      }

      // Set local partition cell offset, the next offset and the number of global partition cells
      grid_.m_localPartitionCellOffsetsRestart[0] = offset;
      grid_.m_localPartitionCellOffsetsRestart[1] = offset + noLocalPartitionCells;
      grid_.m_localPartitionCellOffsetsRestart[2] = noPartitionCellsGlobal;

      for(MInt i = 0; i < noLocalPartitionCells; i++) {
        grid_.m_localPartitionCellLocalIdsRestart[i] = grid_.globalIdToLocalId(partitionCellsGlobalIds[i]);
      }

    } else {
      for(MInt cellId = 0; cellId < cpu_noCells; ++cellId) {
        if(a_hasProperty(cellId, Cell::IsPeriodic) || a_hasProperty(cellId, Cell::IsHalo)) {
          continue;
        }
        if(a_level(cellId) < grid_.m_newMinLevel) continue;

        if(a_hasProperty(cellId, Cell::IsPartitionCell)) {
          partitionCellsFiltered.push_back(changeId.p[cellId]);
        }
      }
    }

    sort(partitionCellsFiltered.begin(), partitionCellsFiltered.end());

    // 5. Compute level difference
    MIntScratchSpace partitionCellLevelDiff(partitionCellsFiltered.size(), AT_, "partitionCellLevelDiff");
    for(MInt i = 0; i < (signed)partitionCellsFiltered.size(); ++i) {
      partitionCellLevelDiff.p[i] = a_level(oldId.p[partitionCellsFiltered[i] - tmpCnt]) - minLevel;
    }

    // 6. Communicate partitionCellsFiltered.size() into partitionCellsFilteredSizePar
    MIntScratchSpace partitionCellsFilteredSizePar(noDomains(), AT_, "partitionCellsFilteredSizePar");
    partitionCellsFilteredSizePar.p[domainId()] = partitionCellsFiltered.size();

    if(1 < noDomains()) {
      tmpCnt = partitionCellsFilteredSizePar.p[domainId()];
      MInt sendVar = static_cast<MInt>(tmpCnt);
      MInt* rcvBuf = partitionCellsFilteredSizePar.getPointer();
      MPI_Allgather(&sendVar, 1, MPI_INT, rcvBuf, 1, MPI_INT, mpiComm(), AT_, "sendVar", "rcvBuf");
      tmpCnt = sendVar;
    }

    /* Write to parallel file
     * ----------------------
     */

    writeParallelGridFile(fileName, cpu_cells, partitionCellsFilteredSizePar.data(), noCellsPar.data(), oldId.data(),
                          changeId.data(), partitionCellsFiltered, partitionCellLevelDiff.data(), changeId2.data());
  }


  /** \brief Sorting after hilbert id
   *
   * 1. Sort the array minLevelCells_hId using combSort and uses the same
   *    sorting on the array minLevelCells_id.
   * 2. Using the no. of Offsprings, create a new id for all minLevelCells
   *    (cell 0 get id 0, cell 1 gets id (id 0 + noOffspring of 0),... )
   *
   * MIntScratchSpace minLevelCells_newId == array containing the new id of the minLevelCells.
   *
   * \tparam[in] T Cell type
   * \param[in] array2Sort the first array to sort
   * \param[in] followUp1 the second array to sort in the same order as the first
   * \param[in] cell2SortCnt the size of the arrays
   **/
  void sortAfterHilbertIds(MLong* array2Sort, MInt* followUp1, MInt cell2SortCnt) {
    TRACE();
    // Comb Sort... (perhaps later change to a better sort)
    MInt gap = cell2SortCnt;
    MBool swapped = true;
    while((gap > 1) || swapped) {
      gap = MInt(gap / 1.247330950103979);
      gap = mMax(gap, 1);
      MInt i = 0;
      swapped = false;
      while(i + gap < cell2SortCnt) {
        if(array2Sort[i] > array2Sort[i + gap]) {
          std::swap(array2Sort[i], array2Sort[i + gap]);
          std::swap(followUp1[i], followUp1[i + gap]);
          swapped = true;
        }
        ++i;
      }
    }
  }


  // -----------------------------------------------------------------------------


  /** \brief Caluclate the number of offsprings and the workload for each cell
   *  \author Lennart Schneiders
   *  \date October 2017
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
  template <typename CELLTYPE>
  void calculateNoOffspringsAndWorkload(Collector<CELLTYPE>* input_cells, MInt input_noCells,
                                        const std::vector<std::vector<MInt>>& cpu_haloCells,
                                        const std::vector<std::vector<MInt>>& cpu_windowCells,
                                        const std::vector<std::vector<MInt>>& partitionLevelAncestorWindowCells,
                                        const std::vector<std::vector<MInt>>& partitionLevelAncestorHaloCells) {
    TRACE();

    MInt recvSize = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      recvSize += (signed)partitionLevelAncestorWindowCells[d].size();
    }
    MIntScratchSpace recvBuffer(mMax(1, recvSize), AT_, "recvBuffer");
    MFloatScratchSpace recvBuffer2(mMax(1, recvSize), AT_, "recvBuffer2");
    MIntScratchSpace tmpNoOffs(input_noCells, AT_, "tmpNoOffs");
    MFloatScratchSpace tmpWorkload(input_noCells, AT_, "tmpWorkload");

    MInt noChilds = IPOW2(nDim);

    for(MInt i = 0; i < input_noCells; ++i) {
      if(a_level(i) < grid_.m_newMinLevel) {
        a_noOffsprings(input_cells, i) = 0;
        a_workload(input_cells, i) = 0.0;
        continue;
      }
      if(!a_hasProperty(input_cells, i, Cell::IsHalo)) {
        a_noOffsprings(input_cells, i) = 1;
        a_workload(input_cells, i) = a_weight(i);
        ASSERT(!std::isnan(a_weight(i)), "cellId " + std::to_string(i) + " g" + std::to_string(a_globalId(i)));
      } else {
        a_noOffsprings(input_cells, i) = 0;
        a_workload(input_cells, i) = 0.0;
      }
    }

    const MInt minLevel = (grid_.m_newMinLevel > 0) ? grid_.m_newMinLevel : m_minLevel;

    for(MInt level_ = m_maxLevel; level_ >= minLevel; --level_) {
      for(MInt i = 0; i < input_noCells; ++i) {
        if(level_ == a_level(i)) {
          if(0 != a_noChildren(i)) {
            for(MInt j = 0; j < noChilds; ++j) {
              if(-1 != a_childId(i, j)) {
                a_noOffsprings(input_cells, i) += a_noOffsprings(input_cells, a_childId(i, j));
                a_workload(input_cells, i) += a_workload(input_cells, a_childId(i, j));
              }
            }
          }
        }
      }
    }


    if(m_maxPartitionLevelShift > 0 && noNeighborDomains() > 0) {
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorHaloCells[d].size(); j++) {
          tmpNoOffs[partitionLevelAncestorHaloCells[d][j]] =
              a_noOffsprings(input_cells, partitionLevelAncestorHaloCells[d][j]);
          tmpWorkload[partitionLevelAncestorHaloCells[d][j]] =
              a_workload(input_cells, partitionLevelAncestorHaloCells[d][j]);
        }
      }
      maia::mpi::reverseExchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                                     mpiComm(), &tmpNoOffs[0], &recvBuffer[0]);
      maia::mpi::reverseExchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                                     mpiComm(), &tmpWorkload[0], &recvBuffer2[0]);
      MInt cnt = 0;
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorWindowCells[d].size(); j++) {
          MInt cellId = partitionLevelAncestorWindowCells[d][j];
          a_noOffsprings(input_cells, cellId) += recvBuffer[cnt];
          a_workload(input_cells, cellId) += recvBuffer2[cnt];
          cnt++;
        }
      }

      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorWindowCells[d].size(); j++) {
          tmpNoOffs[partitionLevelAncestorWindowCells[d][j]] =
              a_noOffsprings(input_cells, partitionLevelAncestorWindowCells[d][j]);
          tmpWorkload[partitionLevelAncestorWindowCells[d][j]] =
              a_workload(input_cells, partitionLevelAncestorWindowCells[d][j]);
        }
      }
      maia::mpi::exchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                              mpiComm(), &tmpNoOffs[0]);
      maia::mpi::exchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                              mpiComm(), &tmpWorkload[0]);
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorHaloCells[d].size(); j++) {
          a_noOffsprings(input_cells, partitionLevelAncestorHaloCells[d][j]) =
              tmpNoOffs[partitionLevelAncestorHaloCells[d][j]];
          a_workload(input_cells, partitionLevelAncestorHaloCells[d][j]) =
              tmpWorkload[partitionLevelAncestorHaloCells[d][j]];
        }
      }
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)cpu_windowCells[i].size(); j++) {
        tmpNoOffs[cpu_windowCells[i][j]] = a_noOffsprings(input_cells, cpu_windowCells[i][j]);
        tmpWorkload[cpu_windowCells[i][j]] = a_workload(input_cells, cpu_windowCells[i][j]);
      }
    }
    if(noNeighborDomains() > 0) {
      maia::mpi::exchangeData(m_nghbrDomains, cpu_haloCells, cpu_windowCells, mpiComm(), &tmpNoOffs[0]);
      maia::mpi::exchangeData(m_nghbrDomains, cpu_haloCells, cpu_windowCells, mpiComm(), &tmpWorkload[0]);
    }
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      for(MInt j = 0; j < (signed)cpu_haloCells[i].size(); j++) {
        a_noOffsprings(input_cells, cpu_haloCells[i][j]) = tmpNoOffs[cpu_haloCells[i][j]];
        a_workload(input_cells, cpu_haloCells[i][j]) = tmpWorkload[cpu_haloCells[i][j]];
      }
    }
  }


  // -----------------------------------------------------------------------------


  /** \brief Create a tree ordering of Ids
   *  \author Lennart Schneiders
   *  \date October 2017
   *
   *  Choose a minLevelCells then sort her childs so that each cell is followed by her childs:
   *
   * PartitionCells - Child0 {(Child0 of Child0) [[(Child0 of Child0 of Child0) ...
   * (Child3 of Child0 of Child0)]] (Child1 of Child0) ... (Child3 of Child0)} Child1 ... Child3
   *
   * \tparam[in] T Cell type
   * \param[in] input_cells collector of the cells
   * \param[in] oldId the array which should hold the tree-like ordering after processing
   * \param[in] newId the array which should hold the tree-like ordering after processing
   * \param[in] minLevelCells_id the minLevelCells
   * \param[in] minLevelCellsCnt the size of the array minLevelCells_id
   **/
  template <typename CELLTYPE>
  MInt createTreeOrderingOfIds(Collector<CELLTYPE>* input_cells, MInt noInternalCells, MInt* oldId, MLong* newId,
                               MInt* minLevelCells_id, MLong* minLevelCells_newId, MInt minLevelCellsCnt, MLong* tmp_id,
                               std::vector<std::vector<MInt>>& partitionLevelAncestorWindowCells,
                               std::vector<std::vector<MInt>>& partitionLevelAncestorHaloCells) {
    TRACE();

    std::fill(&tmp_id[0], &tmp_id[0] + a_noCells(input_cells), -1);
    MIntScratchSpace childOffsprings(a_noCells(input_cells), IPOW2(nDim), AT_, "childOffsprings");
    childOffsprings.fill(0);
    for(MInt cellId = 0; cellId < grid_.m_tree.size(); cellId++) {
      if(a_level(cellId) < grid_.m_newMinLevel) continue;
      if(!a_hasProperty(input_cells, cellId, Cell::IsHalo)) {
        for(MInt k = 0; k < IPOW2(nDim); ++k) {
          if(a_childId(cellId, k) > -1) {
            childOffsprings(cellId, k) = a_noOffsprings(input_cells, a_childId(cellId, k));
          }
        }
      }
    }

    if(m_maxPartitionLevelShift > 0) {
      // gather child no offsprings here...
      MIntScratchSpace noSend(noNeighborDomains(), AT_, "noSend");
      MIntScratchSpace noRecv(noNeighborDomains(), AT_, "noRecv");
      noSend.fill(0);
      noRecv.fill(0);
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorHaloCells[d].size(); j++) {
          MInt haloId = partitionLevelAncestorHaloCells[d][j];
          for(MUint k = 0; k < IPOW2(nDim); ++k) {
            MInt childId = a_childId(haloId, k);
            if(childId > -1) {
              if(!a_hasProperty(input_cells, childId, Cell::IsHalo)) {
                noSend(d)++;
              }
            }
          }
        }
      }
      MIntScratchSpace sendOffsets(noNeighborDomains() + 1, AT_, "sendOffsets");
      sendOffsets(0) = 0;
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        sendOffsets(d + 1) = sendOffsets(d) + noSend(d);
      }
      MIntScratchSpace sendData(sendOffsets(noNeighborDomains()), 3, AT_, "sendData");
      noSend.fill(0);
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt j = 0; j < (signed)partitionLevelAncestorHaloCells[d].size(); j++) {
          MInt haloId = partitionLevelAncestorHaloCells[d][j];
          for(MUint k = 0; k < IPOW2(nDim); ++k) {
            MInt childId = a_childId(haloId, k);
            if(childId > -1) {
              if(!a_hasProperty(input_cells, childId, Cell::IsHalo)) {
                sendData(sendOffsets(d) + noSend(d), 0) = a_noOffsprings(input_cells, childId);
                sendData(sendOffsets(d) + noSend(d), 1) = j;
                sendData(sendOffsets(d) + noSend(d), 2) = k;
                noSend(d)++;
              }
            }
          }
        }
      }

      ScratchSpace<MPI_Request> sendReq(noNeighborDomains(), AT_, "sendReq");
      sendReq.fill(MPI_REQUEST_NULL);
      MInt sendCnt = 0;
      for(MInt c = 0; c < noNeighborDomains(); c++) {
        MPI_Issend(&noSend[c], 1, MPI_INT, m_nghbrDomains[c], 12343, mpiComm(), &sendReq[sendCnt], AT_, "noSend[c]");
        sendCnt++;
      }
      for(MInt c = 0; c < noNeighborDomains(); c++) {
        MPI_Recv(&noRecv[c], 1, MPI_INT, m_nghbrDomains[c], 12343, mpiComm(), MPI_STATUS_IGNORE, AT_, "noRecv[c]");
      }
      if(sendCnt > 0) MPI_Waitall(sendCnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

      MIntScratchSpace recvOffsets(noNeighborDomains() + 1, AT_, "recvOffsets");
      recvOffsets(0) = 0;
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        recvOffsets(d + 1) = recvOffsets(d) + noRecv(d);
      }
      MIntScratchSpace recvData(recvOffsets(noNeighborDomains()), 3, AT_, "recvData");

      sendReq.fill(MPI_REQUEST_NULL);
      sendCnt = 0;
      for(MInt c = 0; c < noNeighborDomains(); c++) {
        if(noSend[c] == 0) continue;
        MPI_Issend(&sendData(sendOffsets[c], 0), 3 * noSend[c], MPI_INT, m_nghbrDomains[c], 12344, mpiComm(),
                   &sendReq[sendCnt], AT_, "sendData(sendOffsets[c]");
        sendCnt++;
      }
      for(MInt c = 0; c < noNeighborDomains(); c++) {
        if(noRecv[c] == 0) continue;
        MPI_Recv(&recvData(recvOffsets[c], 0), 3 * noRecv[c], MPI_INT, m_nghbrDomains[c], 12344, mpiComm(),
                 MPI_STATUS_IGNORE, AT_, "recvData(recvOffsets[c]");
      }
      if(sendCnt > 0) MPI_Waitall(sendCnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt k = 0; k < noRecv[d]; k++) {
          MInt idx = recvData(recvOffsets[d] + k, 1);
          MInt child = recvData(recvOffsets[d] + k, 2);
          ASSERT(child > -1 && child < IPOW2(nDim), "");
          MInt windowId = partitionLevelAncestorWindowCells[d][idx];
          childOffsprings(windowId, child) = recvData(recvOffsets[d] + k, 0);
        }
      }


      if(noNeighborDomains() > 0) {
        maia::mpi::exchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                                mpiComm(), &childOffsprings[0], IPOW2(nDim));
      }
    }

    for(MInt i = 0; i < minLevelCellsCnt; ++i) {
      tmp_id[minLevelCells_id[i]] = minLevelCells_newId[i];
      ASSERT(tmp_id[minLevelCells_id[i]] > -1, "");
    }

    if(m_maxPartitionLevelShift > 0 && noDomains() > 1) {
      maia::mpi::exchangeData(m_nghbrDomains, partitionLevelAncestorHaloCells, partitionLevelAncestorWindowCells,
                              mpiComm(), &tmp_id[0]);
    }

    const MInt minLevel = (grid_.m_newMinLevel > 0) ? grid_.m_newMinLevel : m_minLevel;

    for(MInt level = minLevel; level < m_maxLevel; level++) {
      for(MInt cellId = 0; cellId < grid_.m_tree.size(); cellId++) {
        if(a_level(cellId) == level) {
          // skipping all halo cells except for partition-level-ancestors
          if(tmp_id[cellId] < 0) continue;
          MLong cntId = tmp_id[cellId] + 1;
          for(MInt k = 0; k < IPOW2(nDim); ++k) {
            MInt childId = a_childId(cellId, k);
            if(childId > -1) {
              tmp_id[childId] = cntId;
            }
            cntId += childOffsprings(cellId, k);
          }
        }
      }
    }
    std::map<MLong, MInt> cellMap;
    for(MInt cellId = 0; cellId < noInternalCells; cellId++) {
      if(a_level(cellId) < grid_.m_newMinLevel) continue;
      if(!a_hasProperty(input_cells, cellId, Cell::IsHalo)) {
        ASSERT(tmp_id[cellId] > -1, "");
        cellMap.insert(std::make_pair(tmp_id[cellId], cellId));
      }
    }

    MInt tmpCnt = 0;
    for(auto it = cellMap.begin(); it != cellMap.end(); ++it) {
      oldId[tmpCnt] = it->second;
      newId[tmpCnt] = it->first;
      tmpCnt++;
    }

    return tmpCnt;
  }


  // -----------------------------------------------------------------------------


  /** \brief Save a grid file
      \author Lennart Schneiders
      \date October 2017
    */
  template <typename CELLTYPE>
  void writeParallelGridFile(const MChar* fileName, Collector<CELLTYPE>* input_cells,
                             MInt* partitionCellsFilteredSizePar, MInt* noCellsPar, MInt* oldId, MLong* changeId,
                             const std::vector<MLong>& partitionCellsFiltered, MInt* partitionCellLevelDiff,
                             MInt* changeId2) {
    TRACE();

    // Create File.
    using namespace maia::parallel_io;
    ParallelIo grid(fileName, PIO_REPLACE, mpiComm());

    // Calculate Offsets
    MLongScratchSpace partitionOffset(noDomains(), AT_, "partitionOffset");
    MLongScratchSpace allOffset(noDomains(), AT_, "allOffset");

    MLong partitionCellsFilteredSizeParMax = 0;
    MLong noCellsParMax = 0;
    MLong tmp = 0;

    partitionOffset[0] = 0;
    allOffset[0] = m_32BitOffset;

    for(MInt i = 1; i < noDomains(); ++i) {
      partitionOffset[i] = partitionOffset[i - 1] + partitionCellsFilteredSizePar[i - 1];
      allOffset[i] = allOffset[i - 1] + noCellsPar[i - 1];
    }

    for(MInt i = 0; i < noDomains(); ++i) {
      partitionCellsFilteredSizeParMax += (MLong)partitionCellsFilteredSizePar[i];
      noCellsParMax += (MLong)noCellsPar[i];
      tmp += (MLong)noCellsPar[i];
    }
    // test whether grid can be written out
    m_log << "total number of cells to be written out " << tmp << std::endl;

    MInt minLevel = std::numeric_limits<MInt>::max();
    MInt maxLevel = -1;


    const MInt noCells = grid_.m_newMinLevel > 0 ? m_noInternalCells : noCellsPar[domainId()];

    MInt newCellCount = 0;
    for(MInt i = 0; i < noCells; ++i) {
      MInt cellId = i;
      if(grid_.m_newMinLevel > 0) cellId = changeId2[i];
      if(a_level(oldId[cellId]) < grid_.m_newMinLevel) continue;
      minLevel = mMin(minLevel, a_level(oldId[cellId]));
      maxLevel = mMax(maxLevel, a_level(oldId[cellId]));
      newCellCount++;
    }

    if(grid_.m_newMinLevel > 0) {
      std::cerr << "new-cout " << newCellCount << " old " << noCellsPar[domainId()] << std::endl;
      ASSERT(newCellCount == noCellsPar[domainId()], "");
    }

    MPI_Allreduce(&m_minLevel, &minLevel, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "m_minLevel", "minLevel");
    MPI_Allreduce(&m_maxLevel, &maxLevel, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "m_maxLevel", "maxLevel");


    MInt noMinLevelCells = 0;
    MLong noLeafCells = 0;

    if(grid_.m_newMinLevel > 0) minLevel = grid_.m_newMinLevel;

    std::vector<std::pair<MLong, MInt>> minLevelCells;
    MIntScratchSpace isMinLevelCell(noCellsPar[domainId()], AT_, "isMinLevelCell");
    isMinLevelCell.fill(0);
    for(MInt i = 0; i < noCells; ++i) {
      if(grid_.m_newMinLevel > 0 && changeId2[i] < 0) continue;
      MInt cellId = i;
      if(grid_.m_newMinLevel > 0) cellId = changeId2[i];

      if(a_level(oldId[cellId]) == minLevel) {
        if(grid_.m_newMinLevel < 0) {
          ASSERT(a_parentId(oldId[cellId]) == -1, "");
        }
        minLevelCells.push_back(std::make_pair(changeId[oldId[cellId]], oldId[cellId]));
        isMinLevelCell(cellId) = 1;
        noMinLevelCells++;
      }
      if(a_noChildren(oldId[cellId]) == 0) noLeafCells++;
    }
    MPI_Allreduce(MPI_IN_PLACE, &noLeafCells, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "noLeafCells");
    // sort( minLevelCells.begin(), minLevelCells.end() );

    MIntScratchSpace noMinLevelCellsPerDomain(noDomains(), AT_, "noMinLevelCellsPerDomain");
    MPI_Allgather(&noMinLevelCells, 1, MPI_INT, noMinLevelCellsPerDomain.getPointer(), 1, MPI_INT, mpiComm(), AT_,
                  "noMinLevelCells", "noMinLevelCellsPerDomain.getPointer()");
    MLong minLevelCellOffset = 0;
    for(MInt d = 0; d < domainId(); d++) {
      minLevelCellOffset += (MLong)noMinLevelCellsPerDomain[d];
    }
    MLong noTotalMinLevelCells = 0;
    for(MInt d = 0; d < noDomains(); d++) {
      noTotalMinLevelCells += (MLong)noMinLevelCellsPerDomain[d];
    }

    std::set<MLong> partitionLevelAncestorIds;
    MInt partitionLevelShift = 0;
    MInt maxPartitionLevelShift = 0;
    for(MInt i = 0; i < partitionCellsFilteredSizePar[domainId()]; ++i) {
      MInt levelDiff = a_level(oldId[partitionCellsFiltered[i] - allOffset[domainId()]]) - minLevel;
      ASSERT(levelDiff == partitionCellLevelDiff[i], "");
      partitionLevelShift = mMax(levelDiff, partitionLevelShift);
    }
    if(partitionLevelShift) {
      for(MInt i = 0; i < partitionCellsFilteredSizePar[domainId()]; ++i) {
        MInt cellId = oldId[partitionCellsFiltered[i] - allOffset[domainId()]];
        MInt parentId = a_parentId(cellId);
        while(parentId > -1 && a_level(cellId) >= grid_.m_newMinLevel) {
          if(!a_hasProperty(input_cells, parentId, Cell::IsHalo)) {
            partitionLevelAncestorIds.insert(changeId[parentId]);
            parentId = a_parentId(parentId);
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


    MLong maxNoOffsprings = 0;
    MFloat maxNoCPUs = F0;
    MFloat totalWorkload = F0;
    MFloat maxWorkload = F0;
    MLong partitionCellOffspringThreshold = (MLong)m_partitionCellOffspringThreshold;
    for(MInt i = 0; i < partitionCellsFilteredSizePar[domainId()]; ++i) {
      maxNoOffsprings =
          mMax((MLong)a_noOffsprings(input_cells, oldId[partitionCellsFiltered[i] - allOffset[domainId()]]),
               maxNoOffsprings);
      totalWorkload += a_workload(input_cells, oldId[partitionCellsFiltered[i] - allOffset[domainId()]]);
      maxWorkload =
          mMax(maxWorkload, a_workload(input_cells, oldId[partitionCellsFiltered[i] - allOffset[domainId()]]));
    }
    MPI_Allreduce(MPI_IN_PLACE, &maxNoOffsprings, 1, MPI_LONG, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "maxNoOffsprings");
    MPI_Allreduce(MPI_IN_PLACE, &maxWorkload, 1, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE", "maxWorkload");
    MPI_Allreduce(MPI_IN_PLACE, &totalWorkload, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                  "totalWorkload");
    // maxNoCPUs = noCellsParMax / maxNoOffsprings;
    maxNoCPUs = totalWorkload / maxWorkload;
    MFloat avgWorkload = totalWorkload / ((MFloat)partitionCellsFilteredSizeParMax);
    MFloat avgOffspring = ((MFloat)noCellsParMax) / ((MFloat)partitionCellsFilteredSizeParMax);

    MInt noSolvers = m_tree.noSolvers();

    if(g_multiSolverGrid) {
      MLong solverCount;
      for(MInt b = 0; b < noSolvers; b++) {
        solverCount = 0;
        for(MInt i = 0; i < m_noInternalCells; i++) {
          if(a_level(i) < grid_.m_newMinLevel) continue;
          if(m_tree.solver(i, b) == true) {
            solverCount++;
          }
        }
        MPI_Allreduce(MPI_IN_PLACE, &solverCount, 1, MPI_LONG, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "solverCount");
        grid.setAttributes(&solverCount, "noCells_" + std::to_string(b), 1);
      }
    }

    MInt nodims = nDim;
    grid.setAttributes(&nodims, "nDim", 1);
    grid.setAttributes(&noSolvers, "noSolvers", 1);
    grid.setAttributes(&globalTimeStep, "globalTimeStep", 1);
    grid.setAttributes(&noCellsParMax, "noCells", 1);
    grid.setAttributes(&noLeafCells, "noLeafCells", 1);
    grid.setAttributes(&noTotalMinLevelCells, "noMinLevelCells", 1);
    grid.setAttributes(&partitionCellsFilteredSizeParMax, "noPartitionCells", 1);
    grid.setAttributes(&totalNoPartitionLevelAncestors, "noPartitionLevelAncestors", 1);
    grid.setAttributes(&minLevel, "minLevel", 1);
    grid.setAttributes(&maxLevel, "maxLevel", 1);
    grid.setAttributes(&m_maxUniformRefinementLevel, "maxUniformRefinementLevel", 1);
    grid.setAttributes(&maxPartitionLevelShift, "maxPartitionLevelShift", 1);
    grid.setAttributes(&m_lengthLevel0, "lengthLevel0", 1);
    grid.setAttributes(&m_centerOfGravity[0], "centerOfGravity", nDim);
    grid.setAttributes(&m_boundingBox[0], "boundingBox", 2 * nDim);

    // Add additional multisolver information if grid is reordered by a different Hilbert curve
    if(grid_.m_hasMultiSolverBoundingBox) {
      grid.setAttributes(&grid_.m_targetGridLengthLevel0, "multiSolverLengthLevel0", 1);
      grid.setAttributes(&grid_.m_targetGridMinLevel, "multiSolverMinLevel", 1);
      grid.setAttributes(&(grid_.m_targetGridCenterOfGravity[0]), "multiSolverCenterOfGravity", nDim);
      grid.setAttributes(&(grid_.m_targetGridBoundingBox[0]), "multiSolverBoundingBox", 2 * nDim);
    }

    grid.setAttributes(&m_reductionFactor, "reductionFactor", 1);
    grid.setAttributes(&m_decisiveDirection, "decisiveDirection", 1);
    grid.setAttributes(&totalWorkload, "totalWorkload", 1);
    grid.setAttributes(&maxWorkload, "partitionCellMaxWorkload", 1);
    grid.setAttributes(&avgWorkload, "partitionCellAverageWorkload", 1);
    grid.setAttributes(&maxNoOffsprings, "partitionCellMaxNoOffspring", 1);
    grid.setAttributes(&avgOffspring, "partitionCellAverageNoOffspring", 1);
    grid.setAttributes(&m_partitionCellWorkloadThreshold, "partitionCellWorkloadThreshold", 1);
    grid.setAttributes(&partitionCellOffspringThreshold, "partitionCellOffspringThreshold", 1);
    grid.setAttributes(&maxNoCPUs, "maxNoBalancedCPUs", 1);
    if(m_32BitOffset > 0) {
      grid.setAttributes(&m_32BitOffset, "bitOffset", 1);
    }

    grid.defineArray(PIO_LONG, "partitionCellsGlobalId", partitionCellsFilteredSizeParMax);
    grid.defineArray(PIO_FLOAT, "partitionCellsWorkload", partitionCellsFilteredSizeParMax);

    grid.defineArray(PIO_LONG, "minLevelCellsTreeId", noTotalMinLevelCells);
    grid.defineArray(PIO_LONG, "minLevelCellsNghbrIds", 2 * nDim * noTotalMinLevelCells);

    grid.defineArray(PIO_UCHAR, "cellInfo", noCellsParMax);
    if(m_tree.noSolvers() > 1 || g_multiSolverGrid) {
      grid.defineArray(PIO_UCHAR, "solver", noCellsParMax);
    }

    m_log << "header definition finished" << std::endl;

    // Writing Data.

    // settings Offset to write partitionCells data.
    grid.setOffset(partitionCellsFilteredSizePar[domainId()], partitionOffset[domainId()]);

    // partitionCellsId
    grid.writeArray(&partitionCellsFiltered[0], "partitionCellsGlobalId");

    // partitionCellsWorkLoad.
    {
      MFloatScratchSpace tmp_partitionCellsWorkLoad(partitionCellsFilteredSizePar[domainId()], AT_,
                                                    "tmp_partitionCellsWorkLoad");
      for(MInt i = 0; i < partitionCellsFilteredSizePar[domainId()]; ++i) {
        tmp_partitionCellsWorkLoad[i] =
            a_workload(input_cells, oldId[partitionCellsFiltered[i] - allOffset[domainId()]]);
      }
      grid.writeArray(tmp_partitionCellsWorkLoad.data(), "partitionCellsWorkload");
    }


    // partitionCellsNghbrIds.
    {
      grid.setOffset(m_noDirs * noMinLevelCells, m_noDirs * minLevelCellOffset);
      MLongScratchSpace tmp_nghbrIds(noMinLevelCells, m_noDirs, AT_, "tmp_nghbrIds");
      for(MInt i = 0; i < (signed)minLevelCells.size(); i++) {
        for(MInt j = 0; j < m_noDirs; j++) {
          if(a_neighborId(minLevelCells[i].second, j) > -1) {
            tmp_nghbrIds(i, j) = changeId[a_neighborId(minLevelCells[i].second, j)];
          } else {
            tmp_nghbrIds(i, j) = -1;
          }
        }
      }
      grid.writeArray(tmp_nghbrIds.data(), "minLevelCellsNghbrIds");
    }

    // partitionCellsCoordinates.
    {
      grid.setOffset(noMinLevelCells, minLevelCellOffset);
      MLongScratchSpace tmp_treeId(noMinLevelCells, AT_, "tmp_treeId");
      for(MUint i = 0; i < minLevelCells.size(); i++) {
        maia::grid::hilbert::coordinatesToTreeId<nDim>(
            tmp_treeId[i], &a_coordinate(minLevelCells[i].second, 0), (MLong)grid_.m_targetGridMinLevel,
            &(grid_.m_targetGridCenterOfGravity[0]), grid_.m_targetGridLengthLevel0);
      }
      grid.writeArray(tmp_treeId.data(), "minLevelCellsTreeId");
    }

    // cellInfo
    {
      ScratchSpace<MUchar> tmp_cellInfo(noCellsPar[domainId()], AT_, "tmp_cellInfo");
      MInt cellCount = 0;
      for(MInt i = 0; i < noCells; ++i) {
        MInt cellId = i;
        if(grid_.m_newMinLevel > 0) cellId = changeId2[i];
        if(a_level(oldId[cellId]) < grid_.m_newMinLevel) continue;
        MUint noChilds = (MUint)a_noChildren(oldId[cellId]);
        MUint isMinLevel = (MUint)isMinLevelCell(cellId);
        MUint position = 0;
        MInt parentId = a_parentId(oldId[cellId]);
        if(parentId > -1 && a_level(oldId[cellId]) > grid_.m_newMinLevel) {
          for(MUint j = 0; j < (unsigned)m_maxNoChilds; j++) {
            if(a_childId(parentId, j) == oldId[cellId]) position = j;
          }
        }
        MUint tmpBit = noChilds | (position << 4) | (isMinLevel << 7);
        tmp_cellInfo[cellCount] = static_cast<MUchar>(tmpBit);
        cellCount++;
      }
      grid.setOffset(noCellsPar[domainId()], allOffset[domainId()] - m_32BitOffset);
      grid.writeArray(tmp_cellInfo.data(), "cellInfo");
    }


    // solver
    if(m_tree.noSolvers() > 1 || g_multiSolverGrid) {
      ScratchSpace<MUchar> tmptmp(m_tree.size(), AT_, "tmptmp");
      for(MInt i = 0; i < noCells; ++i) {
        MInt cellId = i;
        if(grid_.m_newMinLevel > 0) cellId = changeId2[i];
        if(a_level(oldId[cellId]) < grid_.m_newMinLevel) continue;
        MUint tmpBit = 0;
        for(MInt solver = 0; solver < m_tree.noSolvers(); solver++) {
          if(m_tree.solver(oldId[cellId], solver)) {
            tmpBit |= (1 << solver);
          }
        }
        tmptmp[cellId] = static_cast<MUchar>(tmpBit);
      }
      grid.setOffset(noCellsPar[domainId()], allOffset[domainId()] - m_32BitOffset);
      grid.writeArray(tmptmp.data(), "solver");
    }
  }


  // -----------------------------------------------------------------------------
};

} // namespace grid
} // namespace maia

//#if defined(MAIA_GCC_COMPILER)
//#pragma GCC diagnostic pop
//#endif

#endif // CARTESIANGRIDIO_H
