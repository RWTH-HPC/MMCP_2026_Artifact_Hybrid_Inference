// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "GRID/structuredpartition.h"
#include <cstring>
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "globals.h"

using namespace std;


/**
 * \brief Constructor for HDF5 grid file name
 *        Reads the size of the coordinate datasets for all blocks
 *
 * \param[in] noBlock Number of blocks in the grid file
 * \param[in] fileName File name of the grid file
 * \param[in] comm MPI communicator
 * \param[in] noDomains_ Number of MPI domains to be used for decomposition
 */
template <MInt nDim>
StructuredDecomposition<nDim>::StructuredDecomposition(const MInt noBlocks_,
                                                       const MString fileName_,
                                                       MPI_Comm comm_,
                                                       const MInt noDomains_)
  : m_noBlocks(noBlocks_),
    m_mpiComm(comm_),
    m_noDomains(noDomains_),
    m_noPartitions(1),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  for(MInt j = 0; j < m_noBlocks; j++) {
    m_blockInfo.push_back(make_unique<PartitionInfo<nDim>>());
  }

  MIntScratchSpace dataSetSize(m_noBlocks * nDim, AT_, "dataSetSize");
  MFloatScratchSpace weight(m_noBlocks, AT_, "weights");

  dataSetSize.fill(0);

  if(globalDomainId() == 0) {
    ParallelIoHdf5 pio(fileName_, maia::parallel_io::PIO_APPEND, MPI_COMM_SELF);
    for(MInt i = 0; i < m_noBlocks; i++) {
      MString sBlockName = "/block" + std::to_string(i);
      std::vector<ParallelIo::size_type> tmp(nDim, 0);
      pio.getArraySize("x", sBlockName, &tmp[0]);
      std::copy(tmp.begin(), tmp.end(), &dataSetSize[i * nDim]);
      weight[i] = 1.0;
      if(pio.hasAttribute("weight", sBlockName)) { // get the weighting factor
        pio.getAttribute<MFloat>(&weight[i], "weight", sBlockName);
      }
    }

    MPI_Bcast(&dataSetSize[0], m_noBlocks * nDim, MPI_INT, 0, m_mpiComm, AT_, "dataSetSize.getPointer()");
    MPI_Bcast(weight.getPointer(), m_noBlocks, MPI_DOUBLE, 0, m_mpiComm, AT_, "weight.getPointer()");
  } else {
    MPI_Bcast(&dataSetSize[0], m_noBlocks * nDim, MPI_INT, 0, m_mpiComm, AT_, "dataSetSize.getPointer()");
    MPI_Bcast(weight.getPointer(), m_noBlocks, MPI_DOUBLE, 0, m_mpiComm, AT_, "weight.getPointer()");
  }

  for(MInt i = 0; i < m_noBlocks; i++) {
    m_blockInfo[i]->blockId = i;
    for(MInt j = 0; j < nDim; j++) {
      m_blockInfo[i]->size[j] = dataSetSize[i * nDim + j] - 1;
      m_blockInfo[i]->weight = weight[i];
      m_blockInfo[i]->offset[j] = 0;
    }
  }
}


/**
 * \brief Constructor for HDF5 grid file handle
 *        Reads the size of the coordinate datasets for all blocks
 *
 * \param[in] noBlock Number of blocks in the grid file
 * \param[in] pio ParallelIO file handle of the grid file
 * \param[in] noDomains_ Number of MPI domains to be used for decomposition
 * \param[in] comm MPI communicator
 */
template <MInt nDim>
StructuredDecomposition<nDim>::StructuredDecomposition(const MInt noBlocks_,
                                                       ParallelIoHdf5* pio,
                                                       MPI_Comm comm_,
                                                       const MInt noDomains_)
  : m_noBlocks(noBlocks_),
    m_mpiComm(comm_),
    m_noDomains(noDomains_),
    m_noPartitions(1),
    m_eps(std::numeric_limits<MFloat>::epsilon()) {
  for(MInt j = 0; j < m_noBlocks; j++) {
    m_blockInfo.push_back(make_unique<PartitionInfo<nDim>>());
  }

  ParallelIo::size_type ijk_max[3] = {0};
  for(MInt i = 0; i < m_noBlocks; i++) {
    m_blockInfo[i]->blockId = i;
    MString sBlockName = "/block" + std::to_string(i);
    pio->getArraySize("x", sBlockName, &ijk_max[0]);

    // assign ijk_max to m_blockInfo and also assign the offset
    for(MInt j = 0; j < nDim; j++) {
      m_blockInfo[i]->size[j] = ijk_max[j] - 1;
      m_blockInfo[i]->offset[j] = 0;
    }
  }
}

/**
 * \brief Helper function for setPartionInfo
 *
 * Goes trough the tree recursively and computes the partition information
 *
 * \param[in] treePointer Pointer to a tree node
 * \param[in] blockId block id of the subtree
 * \param[in] level2dimension1D Array of level to dimension mapping
 * \param[in] partitionCounter Counter of the number of partitions
 * \param[in] beginPartitionCounter
 * \param[in] endPartitionCounter
 * \param[in] divisor
 * \return
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::setPartitionInfoHelper(treeNode*& treePointer,
                                                           const MInt blockId,
                                                           const MInt (&level2dimension1D)[nDim],
                                                           MInt& partitionCounter,
                                                           MInt (&beginPartitionHelper)[nDim],
                                                           MInt (&endPartitionHelper)[nDim],
                                                           MInt (&divisor)[nDim]) {
  treeNode* myTreePointer;
  MInt beginBlock[nDim]{0};
  MInt endBlock[nDim]{0};

  divisor[level2dimension1D[treePointer->level]] = treePointer->noLeaves;

  if(treePointer->level < (nDim - 1)) {
    beginPartitionHelper[level2dimension1D[treePointer->level]] = 0;
    for(MInt i = 0; i < treePointer->noChilds; i++) {
      if(i > 0) {
        beginPartitionHelper[level2dimension1D[treePointer->level]] =
            beginPartitionHelper[level2dimension1D[treePointer->level]] + treePointer->childs[i - 1]->noLeaves;
      }
      endPartitionHelper[level2dimension1D[treePointer->level]] =
          beginPartitionHelper[level2dimension1D[treePointer->level]] + treePointer->childs[i]->noLeaves;
      myTreePointer = treePointer->childs[i];
      setPartitionInfoHelper(myTreePointer, blockId, level2dimension1D, partitionCounter, beginPartitionHelper,
                             endPartitionHelper, divisor);
    }
  } else {
    for(MInt i = 0; i < treePointer->noChilds; i++) {
      beginPartitionHelper[level2dimension1D[treePointer->level]] = i;
      endPartitionHelper[level2dimension1D[treePointer->level]] = i + 1;
      for(MInt j = 0; j < nDim; j++) {
        beginBlock[j] = rounding(((MFloat)m_blockInfo[blockId]->size[j]) * ((MFloat)beginPartitionHelper[j])
                                 / ((MFloat)divisor[j]));
        endBlock[j] =
            rounding(((MFloat)m_blockInfo[blockId]->size[j]) * ((MFloat)endPartitionHelper[j]) / ((MFloat)divisor[j]));
      }
      m_partitionInfo[partitionCounter]->weight = m_blockInfo[blockId]->weight;
      m_partitionInfo[partitionCounter]->totalSize = 1;

      for(MInt j = 0; j < nDim; j++) {
        m_partitionInfo[partitionCounter]->size[j] = endBlock[j] - beginBlock[j] + 1;
        m_partitionInfo[partitionCounter]->totalSize *= (m_partitionInfo[partitionCounter]->size[j]);
        m_partitionInfo[partitionCounter]->offset[j] = beginBlock[j];
      }

      m_partitionInfo[partitionCounter]->blockId = blockId;
      m_partitionInfo[partitionCounter]->partitionId = partitionCounter;
      partitionCounter++;
    }
  }
}

/**
 * \brief Computes partition information from given BCT
 * \param[in] treeRoot Node of a (sub)tree
 * \param[in] level2dimension Array of level to space dimension mapping
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::setPartitionInfo(treeNode**& treeRoot, MIntScratchSpace& level2dimension) {
  treeNode* myTreeRootPointer{};
  MInt beginPartitionHelper[nDim]{0};
  MInt endPartitionHelper[nDim]{0};
  MInt divisor[nDim]{0};

  m_partitionInfo.clear();
  for(MInt j = 0; j < m_noPartitions; j++) {
    m_partitionInfo.push_back(make_unique<PartitionInfo<nDim>>());
  }

  MInt partitionCounter = 0;
  for(MInt i = 0; i < m_noBlocks; i++) {
    myTreeRootPointer = treeRoot[i];
    MInt level2dimension1D[nDim]{0};

    for(MInt j = 0; j < nDim; j++) {
      level2dimension1D[j] = level2dimension(j, i);
    }

    setPartitionInfoHelper(myTreeRootPointer, i, level2dimension1D, partitionCounter, beginPartitionHelper,
                           endPartitionHelper, divisor);
  }
}

/**
 * \brief Recursively traverses the given tree to find the next insertion position
 * \param[in] level2dimension1D Array of level to space dimension mapping
 * \param[in] insertTreePointer Pointer to the node where new node is inserted
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::traverseForInsertionNode(treeNode*& treePointer,
                                                             const MInt blockId,
                                                             const MInt (&level2dimension1D)[nDim],
                                                             treeNode*& insertTreePointer) {
  treeNode* myTreePointer{};

  // compare visited nodes with current candidate insertion node
  // by their normalized number of children
  if(((((MFloat)treePointer->noChilds)
       * ((MFloat)m_blockInfo[blockId]->size[level2dimension1D[insertTreePointer->level]] + 1.0))
      < (((MFloat)insertTreePointer->noChilds)
         * ((MFloat)m_blockInfo[blockId]->size[level2dimension1D[treePointer->level]] + 1.0)))
     || (approx(((MFloat)(treePointer->noChilds)
                 * ((MFloat)(m_blockInfo[blockId]->size[level2dimension1D[insertTreePointer->level]]) + 1.0)),
                ((MFloat)(insertTreePointer->noChilds)
                 * ((MFloat)(m_blockInfo[blockId]->size[level2dimension1D[treePointer->level]]) + 1.0)),
                m_eps)
         && (treePointer->level < insertTreePointer->level))) {
    insertTreePointer = treePointer;
  }

  if(treePointer->level < (nDim - 1)) {
    // look for leftmost child with least leaves
    myTreePointer = treePointer->childs[0];
    for(MInt i = 1; i < treePointer->noChilds; i++) {
      if(treePointer->childs[i]->noLeaves < myTreePointer->noLeaves) {
        myTreePointer = treePointer->childs[i];
      }
    }

    // now traverse in the leftmost child with least leaves
    traverseForInsertionNode(myTreePointer, blockId, level2dimension1D, insertTreePointer);
  }
}


/**
 * \brief Initializes the childs of the given tree node
 * \param[in] treePointer Pointer to the tree
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::initializeChilds(treeNode*& treePointer) {
  if(treePointer->level < (nDim - 1)) {
    treePointer->childs = new treeNode*[treePointer->noChilds];
    for(MInt i = 0; i < treePointer->noChilds; i++) {
      treePointer->childs[i] = new treeNode;
    }
    for(MInt countChild = 0; countChild < treePointer->noChilds; countChild++) {
      treePointer->childs[countChild]->level = (treePointer->level) + 1;
      treeNode* myTreePointer = treePointer->childs[countChild];
      initializeChilds(myTreePointer);
    }
  }
}

/**
 * \brief Recursively destroy all childs of a tree node
 * \param[in] treePointer Pointer to a tree node
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::destroyChilds(treeNode*& treePointer) {
  if(treePointer->level < (nDim - 1)) {
    for(MInt i = 0; i < treePointer->noChilds; i++) {
      treeNode* myTreePointer = treePointer->childs[i];
      destroyChilds(myTreePointer);
    }
    delete[] treePointer->childs;
  }
}

/**
 * \brief Compute the number of leaves of a tree nodes
 * \param[in] treePointer Pointer to a tree node
 * \return Number of leaves
 */
template <MInt nDim>
MInt StructuredDecomposition<nDim>::sumLeaves(treeNode*& treePointer) {
  treeNode* myTreePointer{};
  MInt sumLeaf{};

  if(treePointer->level < (nDim - 1)) {
    sumLeaf = 0;
    for(MInt i = 0; i < treePointer->noChilds; i++) {
      myTreePointer = treePointer->childs[i];
      sumLeaf += sumLeaves(myTreePointer);
    }
  } else {
    sumLeaf = treePointer->noChilds;
  }
  treePointer->noLeaves = sumLeaf;
  return sumLeaf;
}

/**
 * \brief Inserts new child node below given tree node
 * \param[in] insertTreePointer Pointer to the tree node below which new child is inserted
 * \param[in] blockId block id of the current subtree
 * \param[in] level2dimension1D Array of level to space dimension mapping
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::insertChildAtNode(treeNode*& insertTreePointer,
                                                      const MInt blockId,
                                                      const MInt (&level2dimension1D)[nDim]) {
  destroyChilds(insertTreePointer);
  insertTreePointer->noChilds++;
  MInt myNoLeaves = insertTreePointer->noLeaves + 1;
  initializeChilds(insertTreePointer);

  if(insertTreePointer->level < (nDim - 1)) {
    for(MInt i = insertTreePointer->noChilds + 1; i <= myNoLeaves; i++) {
      treeNode* myInsertTreePointer = insertTreePointer;
      traverseForInsertionNode(insertTreePointer, blockId, level2dimension1D, myInsertTreePointer);
      insertChildAtNode(myInsertTreePointer, blockId, level2dimension1D);
      insertTreePointer->noLeaves = sumLeaves(insertTreePointer);
    }
  }
}

/**
 * \brief Add a leaf into the tree
 * \param[in] treeRoot Roots of the trees, same number as number of blocks
 * \param[in] level2dimension Array to map the current level to the space dimension
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::addLeaf(treeNode**& treeRoot, MIntScratchSpace& level2dimension) {
  // select the block id which currently contains the largest partition
  MInt blockId = 0;
  MLongFloat maxWeightedSize = 0;
  for(MInt i = 0; i < m_noPartitions; i++) {
    const MLongFloat weightedSize = m_partitionInfo[i]->totalSize * m_partitionInfo[i]->weight;
    if(weightedSize > maxWeightedSize) {
      blockId = m_partitionInfo[i]->blockId;
      maxWeightedSize = weightedSize;
    }
  }

  // set this tree root for further processing
  treeNode* myTreeRootPointer = treeRoot[blockId];
  treeNode* myInsertTreePointer = treeRoot[blockId];

  // pick corresponding level->dimension map
  MInt level2dimension1D[nDim];
  for(MInt i = 0; i < nDim; i++) {
    level2dimension1D[i] = level2dimension(i, blockId);
  }

  // now find insertion point
  traverseForInsertionNode(myTreeRootPointer, blockId, level2dimension1D, myInsertTreePointer);
  insertChildAtNode(myInsertTreePointer, blockId, level2dimension1D);
  myTreeRootPointer->noLeaves = sumLeaves(myTreeRootPointer);

  // update partition information
  m_noPartitions++;
  setPartitionInfo(treeRoot, level2dimension);
}

/**
 * \brief Decompose the grid into partitions for the set number of domains
 */
template <MInt nDim>
void StructuredDecomposition<nDim>::decompose() {
  treeNode* myTreeRootPointer;
  treeNode** myTreeRoot;
  MInt countBlock;
  m_noPartitions = m_noBlocks;
  myTreeRoot = new treeNode*[m_noBlocks];

  // Create one tree root for each block
  for(MInt i = 0; i < m_noBlocks; i++) {
    myTreeRoot[i] = new treeNode;
  }

  for(countBlock = 0; countBlock < m_noBlocks; countBlock++) {
    myTreeRootPointer = myTreeRoot[countBlock];
    initializeChilds(myTreeRootPointer);
  }

  // prepare by filling level2dimension
  MIntScratchSpace level2dimensionA(nDim, m_noBlocks, AT_, "level2dimensionA");
  level2dimensionA.fill(0);

  for(MInt k = 0; k < m_noBlocks; k++) {
    std::vector<std::pair<MInt, MInt>> dummy{};
    for(MInt i = 0; i < nDim; i++) {
      dummy.push_back({m_blockInfo[k]->size[i], i});
    }
    // sort pairs by first value in descending order
    std::sort(dummy.begin(), dummy.end(), [](auto& left, auto& right) { return left.first > right.first; });
    for(MInt i = 0; i < nDim; i++) {
      level2dimensionA(i, k) = dummy[i].second;
    }
  }

  setPartitionInfo(myTreeRoot, level2dimensionA);

  // now add leaves successively
  while(m_noPartitions < m_noDomains) {
    addLeaf(myTreeRoot, level2dimensionA);
  }

  // Clean up tree
  for(MInt i = 0; i < m_noBlocks; i++) {
    myTreeRootPointer = myTreeRoot[i];
    destroyChilds(myTreeRootPointer);
  }

  for(MInt j = 0; j < m_noBlocks; j++) {
    delete myTreeRoot[j];
  }

  delete[] myTreeRoot;
}

/**
 * \brief Read a precomputed partition info from an external file
 * \return true or false depending on the success
 */
template <MInt nDim>
MBool StructuredDecomposition<nDim>::readFromFile() {
  // temporary scratch to send/receive the partition Info
  MIntScratchSpace filePartitionInfo((3 + 3 + 4) * m_noDomains, "filePartitionInfo", AT_);
  filePartitionInfo.fill(0);
  MInt noCPUs = 0;

  if(globalDomainId() == 0) { // only root reads in the information
    ParallelIoHdf5 pio("partitionFile.hdf5", maia::parallel_io::PIO_APPEND, MPI_COMM_SELF);
    pio.getAttribute<MInt>(&noCPUs, "noCPUs", "");
  }

  // broadcast the number of cpus
  MPI_Bcast(&noCPUs, 1, MPI_INT, 0, m_mpiComm, AT_, "noCPUs");
  if(noCPUs != m_noDomains) {
    m_log << "WARNING no CPUs is not equal to ones given in partitionFile" << endl;
    return false;
  }

  if(globalDomainId() == 0) {
    // read the information out of the file
    ParallelIoHdf5 pio("partitionFile.hdf5", maia::parallel_io::PIO_APPEND, MPI_COMM_SELF);
    ParallelIo::size_type asize[2] = {noCPUs, 10};
    pio.readArray(filePartitionInfo.getPointer(), "", "partitionData", 2, asize);
  }

  MPI_Bcast(filePartitionInfo.getPointer(), noCPUs * (10), MPI_INT, 0, m_mpiComm, AT_,
            "filePartitionInfo.getPointer()");

  m_noPartitions = noCPUs;
  for(MInt j = 0; j < m_noPartitions; j++) {
    m_partitionInfo.push_back(make_unique<PartitionInfo<nDim>>());
  }

  // put the information into the right place
  for(MInt j = 0; j < m_noPartitions; j++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      m_partitionInfo[j]->size[dim] = filePartitionInfo[j * 10 + dim];
      m_partitionInfo[j]->offset[dim] = filePartitionInfo[j * 10 + 3 + dim];
    }

    m_partitionInfo[j]->blockId = filePartitionInfo[j * 10 + 6];
    m_partitionInfo[j]->partitionId = filePartitionInfo[j * 10 + 7];
    m_partitionInfo[j]->cpu = filePartitionInfo[j * 10 + 8];
    m_partitionInfo[j]->totalSize = filePartitionInfo[j * 10 + 9];
  }

  return true;
}

template <MInt nDim>
MInt StructuredDecomposition<nDim>::rounding(MFloat x) {
  MFloat a;
  ((x - floor(x)) >= 0.5) ? a = ceil(x) : a = floor(x);
  return (MInt)a;
}

// Explicit instantiations for 2D and 3D
template class StructuredDecomposition<2>;
template class StructuredDecomposition<3>;
