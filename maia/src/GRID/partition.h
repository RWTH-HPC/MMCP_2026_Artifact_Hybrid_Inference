// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PARTITION_H_
#define PARTITION_H_

#include <algorithm>
#include <numeric>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "typetraits.h"

namespace maia {
namespace grid {

/// \brief Translates a list of partition cell offsets to global offsets.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Jerry Grimmen
/// \note changed, since unnecessarily complicated, Lennart
/// \date 2014-12-02
///
/// \tparam WeightType Numeric type for weight information.
/// \tparam IdType Integer type to store counters and size information.
/// \param[in] partitionCellOffsets List of partition cell offsets (access by domain id).
/// \param[in] partitionCellToGlobalIds List of global ids for each partition cell (access
///                               by partition cell id).
/// \param[in] noDomains Number of domains to which the cells are distributed.
/// \param[out] globalIdOffsets Pointer to list that will hold the global id of
///             the first cell that goes into the respective domain. Must hold
///             at least `noDomains` elements.
///
/// This algorithm is only useful for partition cells and is generally not applicable
/// to other problems.
template <typename IdType>
void partitionCellToGlobalOffsets(const IdType* const partitionCellOffsets,
                                  const IdType* const partitionCellToGlobalIds,
                                  const IdType noDomains,
                                  IdType* const globalIdOffsets) {
  // Set global offset to zero for first domain
  globalIdOffsets[0] = 0;

  // Calculate global offsets (domain 0 is skipped since it is already zero)
  for(IdType domainId = 1; domainId < noDomains; domainId++) {
    globalIdOffsets[domainId] = partitionCellToGlobalIds[partitionCellOffsets[domainId]];
  }
}


/** \brief Serial algorithm to find the optimal (not ideal) partitioning with given workloads based on a probing
   algorithm with bisection \author Lennart Schneiders \date 07.2017 The algorithm guarantees to find the optimal
   workload distribution. To speed up, the termination criterion <eps> is used to stop the bisection iterations after
   the partitioning has converged to a workload imbalance not more that 0.01% above its optimal value. For details on
   the overall scheme see 'Scalable high-quality 1D partitioning', by M. Lieber and W.E. Nagel (HPCS 2014)
  */
template <typename WeightType, typename IdType>
WeightType optimalPartitioningSerial(const WeightType* const weights, const IdType noWeights, const IdType noPartitions,
                                     IdType* const offsets) {
  WeightType totalWorkload = F0;
  WeightType maxWorkload = F0;
  ScratchSpace<WeightType> prefixSum(noWeights + 1, AT_, "prefixSum");
  prefixSum(0) = F0;
  for(IdType i = 0; i < noWeights; i++) {
    prefixSum(i + 1) = prefixSum(i) + weights[i];
    maxWorkload = mMax(maxWorkload, weights[i]);
  }
  totalWorkload = prefixSum(noWeights);
  WeightType optimalWorkload = totalWorkload / ((WeightType)noPartitions);
  ScratchSpace<IdType> offsetTmp(noPartitions + 1, AT_, "offsetTmp");
  WeightType eps = mMin(0.0001 * optimalWorkload, 0.5 * maxWorkload / optimalWorkload);
  // WeightType eps = 0.0001*optimalWorkload;
  IdType noPartitionsOpt = 0;
  WeightType maxPartitionWorkload = F0;
  for(IdType i = 0; i <= noPartitions; i++) {
    offsets[i] = -1;
  }
  WeightType lowerVal = mMax(maxWorkload, optimalWorkload);
  WeightType upperVal = lowerVal + mMax(maxWorkload, F2 * eps);
  WeightType probeWorkload = lowerVal; // try the optimal workload first
  WeightType workloadUsed = -F1;
  MInt itCnt = 0;
  if(noWeights <= noPartitions) {
    std::cerr << "Warning: less or equal number of weights than number of partitions!" << std::endl;
    for(IdType i = 0; i <= noPartitions; i++) {
      offsets[i] = mMin(noPartitions, i);
    }
    for(IdType i = 0; i < noPartitions; i++) {
      maxPartitionWorkload = mMax(maxPartitionWorkload, weights[i]);
    }
    return maxPartitionWorkload;
  }
  while(upperVal - lowerVal > eps) {
    offsetTmp.fill(-1);
    offsetTmp(0) = 0;
    IdType id = noWeights / noPartitions; // initial offset guess
    IdType lastCpu = 1;
    for(IdType cpu = 1; cpu < noPartitions; cpu++) {
      MBool found = false;
      while(!found) {
        if(prefixSum(id) - prefixSum(offsetTmp(cpu - 1)) > probeWorkload)
          id--;
        else if(!(prefixSum(id + 1) - prefixSum(offsetTmp(cpu - 1)) > probeWorkload))
          id++;
        else
          found = true;
        if(id >= noWeights) {
          cpu = noPartitions;
          id = noWeights;
          found = true;
        }
      }
      if(cpu < noPartitions) {
        while((noWeights - id) < (noPartitions - cpu))
          id--;
        id = mMax(offsetTmp[cpu - 1] + 1, id);
        offsetTmp(cpu) = id;
        lastCpu = cpu;
        id = mMin(noWeights - 1, 2 * offsetTmp(cpu) - offsetTmp(cpu - 1)); // initial offset guess
      }
    }
    if(prefixSum(noWeights) - prefixSum(offsetTmp(lastCpu)) > probeWorkload) {
      lowerVal = probeWorkload;
    } else {
      std::copy(offsetTmp.begin(), offsetTmp.begin() + noPartitions + 1, offsets);
      workloadUsed = probeWorkload;
      noPartitionsOpt = lastCpu + 1;
      upperVal = probeWorkload;
    }
    if(itCnt > 0 && offsets[0] == -1) {
      upperVal += F2 * maxWorkload; // enlarge the interval if no valid distribution was found previously
    }
    probeWorkload = F1B2 * (lowerVal + upperVal);
    itCnt++;
  }
  probeWorkload = workloadUsed;

  //--
  for(IdType i = noPartitionsOpt; i <= noPartitions; i++) {
    offsets[i] = noWeights;
  }

  // fine tuning, avoid domains with tiny workload, however, this is not affecting the overall imbalance
  for(IdType i = 1; i < noPartitionsOpt; i++) {
    while((prefixSum[offsets[i + 1]] - prefixSum[offsets[i]] + weights[offsets[i] - 1])
              < (prefixSum[offsets[i]] - prefixSum[offsets[i - 1]] - weights[offsets[i] - 1])
          && offsets[i] > offsets[i - 1] + 1) {
      offsets[i]--;
    }
  }
  // This may occur for strong local refinement and oddly leaves one or more last ranks empty, i.e., they are not even
  // required for the optimal partitioning. Note, however, that the optimal partitioning is usually far away from the
  // ideal partitioning (totalWorkload/noRanks) for these cases, in the first place! In this case, try to redistribute
  // the workload with a certain penalty, so they are not idling.
  if(noPartitionsOpt < noPartitions) {
    // try to split existing partitions approximately in half
    WeightType probeSplit = 0.5;
    IdType diff = noPartitions - noPartitionsOpt;
    while(diff > 0) {
      probeSplit -= 0.05;
      const IdType diff0 = diff;
      for(IdType k = 0; k < diff0; k++) {
        for(IdType i = 0; i < noPartitionsOpt; i++) {
          if(prefixSum[offsets[i + 1]] - prefixSum[offsets[i]] > F2 * probeSplit) {
            IdType id = (offsets[i + 1] + offsets[i]) / 2;
            MBool found = false;
            while(!found && id > offsets[i] && id < offsets[i + 1]) {
              if(prefixSum(id) - prefixSum(offsets[i]) < probeSplit * probeWorkload)
                id++;
              else if(!(prefixSum(id - 1) - prefixSum(offsets[i]) < probeSplit * probeWorkload))
                id--;
              else
                found = true;
            }
            //            std::cout << offsets[i] << ":" << prefixSum(offsets[i]) << " " << id << ":"
            //                      << prefixSum(id) << " " << offsets[i + 1] << ":" << prefixSum(offsets[i + 1])
            //                      << " " << probeSplit << " " << probeWorkload << std::endl;
            if(id > offsets[i] && id < offsets[i + 1]
               && prefixSum(id) - prefixSum(offsets[i]) > probeSplit * probeWorkload
               && prefixSum(offsets[i + 1]) - prefixSum(id) > probeSplit * probeWorkload) {
              //              std::cout << "DEBUG1" << std::endl;
              for(IdType j = noPartitions; j > i + 1; j--) {
                offsets[j] = offsets[j - 1]; // shift right
              }
              offsets[i + 1] = id; // set intermediate offset
              i = noPartitionsOpt;
              diff--;
            }
          }
        }
      }
      if(probeSplit <= 0) {
        std::cout << "DIFF " << diff << std::endl;
        mTerm(1, AT_, "partition algorithm failed due to redistribution error");
      }
    }
  }

  offsets[noPartitions] = noWeights;

  // Return the maximum local workload such that the quality of the final partitioning can be determined later on.
  maxPartitionWorkload = F0;
  for(IdType i = 0; i < noPartitions; i++) {
    maxPartitionWorkload = mMax(maxPartitionWorkload, (prefixSum[offsets[i + 1]] - prefixSum[offsets[i]]));
  }
  return maxPartitionWorkload;
}


/** \brief Fully parallel partitioning into <noGroups> partitions based on a heuristical approach
    \author Lennart Schneiders
    \date 07.2017

    Creates a rough estimate of the domain boundaries based on heuristics.
    This is used to create a coarse partitioning into several groups.
    Within the groups, a more accurate partition algorithm will be executed at a later stage.
    Note that for noGroups=1 one obtains a domain decompositioning similar to the one we had previously.
    For details on the overall scheme see 'Scalable high-quality 1D partitioning', by M. Lieber and W.E. Nagel (HPCS
   2014)
  */
template <typename WeightType, typename IdType>
void heuristicPartitioningParallel(const WeightType* const localWeights, const IdType noLocalWeights,
                                   const MLong localOffset, const MLong noGlobalWeights, const WeightType totalWorkload,
                                   const IdType noGlobalPartitions, const IdType globalPartitionId,
                                   const IdType noGroups, const MPI_Comm comm, MPI_Comm& commGroup,
                                   MPI_Comm& commMasters, IdType& groupId, IdType& groupLocalId,
                                   MLong* const groupOffsets) {
  // create some sub-communicators for later local data processing
  const IdType noRanksPerGroup = noGlobalPartitions / noGroups;
  groupId = (globalPartitionId >= noGroups * noRanksPerGroup) ? noGroups - 1 : globalPartitionId / noRanksPerGroup;
  groupLocalId = globalPartitionId - groupId * noRanksPerGroup;
  MBool isGroupMaster = (groupLocalId == 0);
  IdType masterTag = isGroupMaster ? 0 : MPI_UNDEFINED;

  MPI_Comm_split(comm, groupId, globalPartitionId, &commGroup, AT_, "commGroup");
  MPI_Comm_split(comm, masterTag, globalPartitionId, &commMasters, AT_, "commMaster");

  // compute global workload offsets/prefixes and the weights of the various groups according to their number of ranks
  WeightType localOffset0 = F0;
  WeightType localWorkload = F0;
  for(IdType i = 0; i < noLocalWeights; i++) {
    localWorkload += localWeights[i];
  }
  MPI_Exscan(&localWorkload, &localOffset0, 1, MPI_DOUBLE, MPI_SUM, comm, AT_, "localWorkload", "localOffset0");
  ScratchSpace<WeightType> groupWeight(noGroups, AT_, "groupWeight");
  for(IdType i = 0; i < noGroups; i++)
    groupWeight[i] = ((WeightType)(i * noRanksPerGroup)) / ((WeightType)noGlobalPartitions);

  // check heuristically for partition boundaries within my local data range and store it
  std::vector<std::pair<IdType, MLong>> localPartitionCellOffsets;
  localWorkload = localOffset0;
  for(IdType i = 0; i < noLocalWeights; i++) {
    WeightType f0 = localWorkload / totalWorkload;
    WeightType f1 = (localWorkload + localWeights[i]) / totalWorkload;
    for(IdType g = 1; g < noGroups; g++) {
      if(f0 <= groupWeight(g) && f1 > groupWeight(g)) {
        // localPartitionCellOffsets.push_back( make_pair(g,localOffset+i+1) );
        if(fabs(groupWeight(g) - f0) < fabs(f1 - groupWeight(g)))
          localPartitionCellOffsets.push_back(std::make_pair(g, localOffset + i));
        else
          localPartitionCellOffsets.push_back(std::make_pair(g, localOffset + i + 1));
      }
    }
    localWorkload += localWeights[i];
  }

  // send all local partition offsets to rank 0
  IdType recvCount = 0;
  IdType noOffsetsFound = (signed)localPartitionCellOffsets.size();
  if(noOffsetsFound > 0) {
    MLongScratchSpace sendBuf(noOffsetsFound, 2, AT_, "sendBuf");
    for(IdType i = 0; i < noOffsetsFound; i++) {
      sendBuf(i, 0) = localPartitionCellOffsets[i].first;
      sendBuf(i, 1) = localPartitionCellOffsets[i].second;
    }
    if(globalPartitionId == 0) {
      for(IdType i = 0; i < noOffsetsFound; i++) {
        groupOffsets[sendBuf(i, 0)] = sendBuf(i, 1);
        recvCount++;
      }
    } else {
      MPI_Send(&sendBuf[0], 2 * noOffsetsFound, MPI_LONG, 0, 100, comm, AT_, "sendBuf[0]");
    }
  }

  // rank 0 receives all global offsets
  if(globalPartitionId == 0) {
    const WeightType time0 = MPI_Wtime();
    WeightType time1 = MPI_Wtime();
    groupOffsets[0] = 0;
    groupOffsets[noGroups] = noGlobalWeights;
    while(recvCount < noGroups - 1) {
      MPI_Status status;
      MInt recvSize;
      MInt flag;
      MPI_Iprobe(MPI_ANY_SOURCE, 100, comm, &flag, &status);
      if(flag) {
        MPI_Get_count(&status, MPI_LONG, &recvSize, AT_);
        MInt source = status.MPI_SOURCE;
        MInt recvOffsets = recvSize / 2;
        MLongScratchSpace recvBuf(recvOffsets, 2, AT_, "recvBuf");
        MPI_Recv(&recvBuf[0], recvSize, MPI_LONG, source, 100, comm, &status, AT_, "recvBuf[0]");
        for(IdType i = 0; i < recvOffsets; i++) {
          groupOffsets[recvBuf(i, 0)] = recvBuf(i, 1);
        }
        recvCount += recvOffsets;
      }
      if((MInt)((MPI_Wtime() - time0) / 10) > (MInt)((time1 - time0) / 10)) {
        std::cerr << "Rank " << globalPartitionId << " already waiting " << MPI_Wtime() - time0
                  << " seconds for incoming group offsets..." << std::endl;
      }
      time1 = MPI_Wtime();
    }
  }

  // distribute the offsets
  MPI_Bcast(&groupOffsets[0], noGroups + 1, MPI_LONG, 0, comm, AT_, "groupOffsets[0]");
}


template <typename WeightType, typename IdType>
void partitionParallelSplit(const WeightType* const localWeights,
                            const IdType noLocalWeights,
                            const IdType localWeightOffset,
                            const IdType noDomains,
                            const IdType domainId,
                            const MPI_Comm mpiComm,
                            IdType* const offsets) {
  // Compute sum of local weights
  WeightType localTotalWeight = WeightType();
  for(IdType i = 0; i < noLocalWeights; i++) {
    localTotalWeight += localWeights[i];
  }

  // Compute global prefix sum offsets
  WeightType prefixSumOffset = WeightType();
  MPI_Exscan(&localTotalWeight, &prefixSumOffset, 1, type_traits<WeightType>::mpiType(), MPI_SUM, mpiComm, AT_,
             "localTotalWeight", "prefixSumOffset");

  // Exchange total global weight from last process
  WeightType globalWeight = prefixSumOffset + localTotalWeight;
  MPI_Bcast(&globalWeight, 1, type_traits<WeightType>::mpiType(), noDomains - 1, mpiComm, AT_, "globalWeight");

  // 'Optimal' weight on a single domain
  const WeightType optimalWeight = globalWeight / static_cast<WeightType>(noDomains);

  // Set all offsets to zero
  std::fill_n(offsets, noDomains, IdType());

  // Compute the current domain
  IdType currentDomainId = std::floor(prefixSumOffset / optimalWeight);

  ScratchSpace<IdType> foundDomain(noDomains, AT_, "foundDomain");
  std::fill(foundDomain.begin(), foundDomain.end(), IdType());

  // Last prefix sum value, start with global offset
  WeightType lastPrefixSum = prefixSumOffset;
  for(IdType i = 1; i < noLocalWeights + 1; i++) {
    // Increase prefix sum by last weight (exclusive prefix sum)
    const WeightType prefixSum = localWeights[i - 1] + lastPrefixSum;
    // Determine optimal splitter, i.e. the next domain offset
    const WeightType splitter = (currentDomainId + 1) * optimalWeight;
    // Check if splitter position reached
    if(!(prefixSum < splitter)) {
      // Set offset of next domain
      currentDomainId++;
      offsets[currentDomainId] = localWeightOffset + i;

      // Compare distances between splitter and prefix sum, and splitter and the last prefix sum
      const WeightType distance = ABS(prefixSum - splitter);
      const WeightType lastDistance = ABS(lastPrefixSum - splitter);
      // If the previous distance is closer to the optimal splitter decrease the offset by one
      // This limits the maximum domain weight deviation by w_max/w_opt, with w_max the maximum
      // weight of a single element/object and w_opt the optimal weight for a domain.
      // However, only do this if the previous domain is not empty afterwards, i.e. has the same
      // offset.
      if(lastDistance < distance && offsets[currentDomainId] > offsets[currentDomainId - 1] + 1) {
        offsets[currentDomainId]--;
      }
      // Store the domain id on which the offset was found
      foundDomain[currentDomainId] = domainId;
    }
    lastPrefixSum = prefixSum;

    if(currentDomainId == noDomains - 1) {
      break;
    }
  }

  if(lastPrefixSum + 0.0001 * optimalWeight > (currentDomainId + 1) * optimalWeight) {
    offsets[currentDomainId + 1] = localWeightOffset + noLocalWeights;
    foundDomain[currentDomainId + 1] = domainId;
    std::cerr << "DEBUG partition: set offset " << currentDomainId + 1 << " on domain " << domainId
              << " with last local partition cell " << localWeightOffset + noLocalWeights << " " << lastPrefixSum << " "
              << (currentDomainId + 1) * optimalWeight << std::endl;
  }

  // Take max of domains in case an offset is found twice which can happen in special situations
  MPI_Allreduce(MPI_IN_PLACE, &foundDomain[0], noDomains, type_traits<IdType>::mpiType(), MPI_MAX, mpiComm, AT_,
                "MPI_IN_PLACE", "foundDomain");
  for(IdType i = 0; i < noDomains; i++) {
    if(offsets[i] > 0 && foundDomain[i] != domainId) {
      std::cerr << "DEBUG partition: reset offset " << i << " on domain " << domainId << ", also found on domain "
                << foundDomain[i] << std::endl;
      offsets[i] = 0;
    }
  }

  // Combine offsets of all processes
  MPI_Allreduce(MPI_IN_PLACE, &offsets[0], noDomains, type_traits<IdType>::mpiType(), MPI_MAX, mpiComm, AT_,
                "MPI_IN_PLACE", "offsets");

  // Check for invalid offsets!
  for(IdType i = 1; i < noDomains; i++) {
    if(offsets[i] <= offsets[i - 1]) {
      m_log << "Partition: invalid offsets " << i - 1 << " " << offsets[i - 1] << " and " << i << " " << offsets[i]
            << "; setting offset " << i << " to " << offsets[i - 1] + 1 << std::endl;
      offsets[i] = offsets[i - 1] + 1;
    }
  }

  // Storage for total domain weights
  ScratchSpace<WeightType> domainWeights(noDomains, AT_, "domainWeights");
  std::fill(domainWeights.begin(), domainWeights.end(), WeightType());

  // Find local weight offset in determined domain offsets
  IdType* domainIt = std::lower_bound(&offsets[0], &offsets[0] + noDomains, localWeightOffset);
  // Target domain of first local weight on this domain
  IdType targetDomain = (*domainIt == localWeightOffset) ? std::distance(&offsets[0], domainIt)
                                                         : std::distance(&offsets[0], domainIt) - 1;

  // Add weights to corresponding domains
  for(IdType i = 0; i < noLocalWeights; i++) {
    domainWeights[targetDomain] += localWeights[i];

    // Increase domain id if the next global weight index corresponds to the next domain offset
    if(targetDomain < noDomains - 1 && localWeightOffset + i + 1 == offsets[targetDomain + 1]) {
      targetDomain++;
    }
  }

  // Combine locally determined domain weights
  MPI_Allreduce(MPI_IN_PLACE, &domainWeights[0], noDomains, type_traits<WeightType>::mpiType(), MPI_SUM, mpiComm, AT_,
                "MPI_IN_PLACE", "domainWeights");

  WeightType maxWeight = WeightType();
  WeightType minWeight = domainWeights[0];

  // DEBUG output
  for(IdType i = 0; i < noDomains; i++) {
    /* m_log << "domain weight " << i << " " << domainWeights[i] << endl; */
    maxWeight = std::max(maxWeight, domainWeights[i]);
    minWeight = std::min(minWeight, domainWeights[i]);
  }
  m_log << "maximum domain weight " << maxWeight << " (max/opt = " << maxWeight / optimalWeight << ")" << std::endl;
  m_log << "minimum domain weight " << minWeight << " (min/opt = " << minWeight / optimalWeight << ")" << std::endl;
}

} // namespace grid
} // namespace maia

#endif // PARTITION_H_
