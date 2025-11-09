// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSOLVERSTRUCTPARTITION_H
#define FVSOLVERSTRUCTPARTITION_H

#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include <memory>

class ParallelIoHdf5;

class treeNode {
 public:
  treeNode(){};
  ~treeNode(){};
  MInt level = 0;
  MInt noChilds = 1;
  MInt noLeaves = 1;
  treeNode** childs = nullptr;
};

template <MInt nDim>
class PartitionInfo {
 public:
  PartitionInfo() {
    size = std::make_unique<MInt[]>(nDim);
    offset = std::make_unique<MInt[]>(nDim);
  };
  ~PartitionInfo(){};

  MInt blockId = -1;
  MInt partitionId = -1;
  MInt cpu = -1;
  std::unique_ptr<MInt[]> size{};
  std::unique_ptr<MInt[]> offset{};
  MLong totalSize = 0;
  MFloat weight = 1.0;
};

/**
 * \brief Class for the decomposition (partition) of structured grids
 *
 * The class will decompose a given multi-block grid into a number of partitions
 * using the balanced cut trees method described in
 *
 * G. Geiser, W. Schröder,
 * Structured multi-block grid partitioning using balanced cut trees,
 * J. Parallel Distrib. Comput. 138 (2020)
 *
 * Note that throughout the class the term 'block' and the corresponding 'blockId'
 * refer to input block, i.e., the blocks of the grid file, whereas the
 * term 'partition' refers to the partitions created by the decomposition. Finally,
 * the number of partitions should be equal to the given number of MPI domains.
 */
template <MInt nDim>
class StructuredDecomposition {
 public:
  StructuredDecomposition(const MInt, const MString, MPI_Comm, const MInt);
  StructuredDecomposition(const MInt, ParallelIoHdf5* pio, MPI_Comm, const MInt);
  ~StructuredDecomposition(){};
  void decompose();
  MBool readFromFile();
  MInt getBlockIdFromPartition(const MInt domainId_) { return m_partitionInfo[domainId_]->blockId; };
  MInt getPartitionOffset(const MInt domainId_, const MInt dim) { return m_partitionInfo[domainId_]->offset[dim]; };
  MInt getPartitionSize(const MInt domainId_, const MInt dim) { return m_partitionInfo[domainId_]->size[dim]; };
  MInt getBlockOffset(const MInt domainId_, const MInt dim) { return m_blockInfo[domainId_]->offset[dim]; };
  MInt getBlockSize(const MInt domainId_, const MInt dim) { return m_blockInfo[domainId_]->size[dim]; };

 private:
  MInt rounding(MFloat x);
  void initializeChilds(treeNode*& treePointer);
  void setPartitionInfoHelper(treeNode*&, const MInt, const MInt (&)[nDim], MInt&, MInt (&)[nDim], MInt (&)[nDim],
                              MInt (&)[nDim]);
  void setPartitionInfo(treeNode**&, MIntScratchSpace&);
  void traverseForInsertionNode(treeNode*& treePointer, const MInt, const MInt (&)[nDim], treeNode*&);
  void destroyChilds(treeNode*& treePointer);
  MInt sumLeaves(treeNode*& treePointer);
  void insertChildAtNode(treeNode*&, const MInt, const MInt (&)[nDim]);
  void addLeaf(treeNode**& treeRoot, MIntScratchSpace& level2dimension);
  MInt noDomains() const { return m_noDomains; }

  const MInt m_noBlocks;
  const MPI_Comm m_mpiComm;
  const MInt m_noDomains;
  MInt m_noPartitions;
  const MFloat m_eps;

  std::vector<std::unique_ptr<PartitionInfo<nDim>>> m_blockInfo{};
  std::vector<std::unique_ptr<PartitionInfo<nDim>>> m_partitionInfo{};
};

#endif
