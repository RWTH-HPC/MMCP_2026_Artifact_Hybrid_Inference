// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangridcontroller.h"
#include "UTIL/maiamath.h"
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wduplicated-branches"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-dtor"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-result"
#endif
#ifdef MAIA_CLANG_COMPILER
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wextra-semi-stmt"
#pragma clang diagnostic ignored "-Wfloat-equal"
#endif
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic pop
#endif
#ifdef MAIA_CLANG_COMPILER
#pragma clang diagnostic pop
#endif

/// \brief Compute computational weights for different components in the
///        simulation based on the current distribution and the domain loads.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// \param[in] loads Domain loads.
/// \param[in] domainWeight Local domain weight or processing capacity.
/// \param[out] weights Vector of computed weights.
template <MInt nDim>
void maia::grid::Controller<nDim>::computeWeights(const MFloat* loads, const MFloat domainWeight,
                                                  std::vector<MFloat>& weights) {
  TRACE();
  using namespace maia;

  // Number of load quantities (should be the same for every rank)
  const MInt noLoadTypes = globalNoLoadTypes();
  if(noLoadTypes < 1) {
    TERMM(1,
          "There are no load quantities present in the solver(s), but needed to determine a new "
          "partition for dynamic load balancing!");
  }
  weights.resize(noLoadTypes);
  std::fill(weights.begin(), weights.end(), 1.0);

  MIntScratchSpace loadQuantities(noLoadTypes, AT_, "loadQuantities");
  MFloatScratchSpace loadQuantitiesScaled(noLoadTypes, AT_, "loadQuantitiesScaled");
  MFloatScratchSpace globalLoadQuantities(noDomains(), noLoadTypes, AT_, "globalLoadQuantities");

  // Get the solver specific load quantities
  getLoadQuantities(&loadQuantities[0]);

  // Scale the load quantities with the domain weight
  for(MInt i = 0; i < noLoadTypes; i++) {
    loadQuantitiesScaled[i] = domainWeight * loadQuantities[i];
  }

  // Gather (scaled) load quantities on rank 0
  MPI_Gather(&loadQuantitiesScaled[0], noLoadTypes, type_traits<MFloat>::mpiType(), &globalLoadQuantities[0],
             noLoadTypes, type_traits<MFloat>::mpiType(), 0, gridb().mpiComm(), AT_, "loadQuantitiesScaled[0]",
             "globalLoadQuantities[0]");

  if(domainId() == 0) {
    // Estimate the parameters, i.e. the new weights
    estimateParameters(noDomains(), noLoadTypes, &globalLoadQuantities[0], &loads[0], &weights[0]);
    // limit weights for certain types
    if(!m_testDynamicLoadBalancing) {
      MInt offset = 0;
      for(MInt i = 0; i < noSolvers(); i++) {
        solver(i).limitWeights(&weights[offset]);
        offset += solver(i).noLoadTypes();
      }
      // scale weights
      const MFloat norm = std::accumulate(&weights[0], &weights[0] + noLoadTypes, 0.0);
      for(MInt i = 0; i < noLoadTypes; i++) {
        weights[i] /= norm;
      }
    }
  }

  // Broadcast the estimated parameters/weights
  MPI_Bcast(&weights[0], noLoadTypes, type_traits<MFloat>::mpiType(), 0, gridb().mpiComm(), AT_, "weights[0]");

  /* storeLoadsAndWeights(&loads[0], noLoadTypes, &loadQuantities[0], domainWeight, &weights[0]); */

  { // Output to m_log
    m_log << " * computed weights:";
    MInt count = 0;
    for(MInt i = 0; i < noSolvers(); i++) {
      const MInt solverCount = solver(i).noLoadTypes();
      std::vector<MString> names(solverCount);
      std::vector<MFloat> weightsTmp(solverCount);

      solver(i).getDefaultWeights(&weightsTmp[0], names);
      for(MInt j = 0; j < solverCount; j++) {
        m_log << " s" << i << "_" << names[j] << ":" << count;
        count++;
      }
    }

    m_log << std::endl << " * computed weights:";
    for(MInt i = 0; i < noLoadTypes; i++) {
      TERMM_IF_COND(std::isnan(weights[i]), "computed weight is NaN");
      m_log << " " << std::scientific << weights[i];
    }
    m_log << std::endl;
  }

  // use static pre-defined weights instead
  if(m_dlbStaticWeights != nullptr && m_dlbStaticWeightMode > -1) {
    if(m_dlbStaticWeightMode == 0
       || (m_dlbStaticWeightMode > 0 && m_dlbStep < m_dlbStaticWeightMode && m_dlbStaticWeightMode != 98)
       || (m_dlbStaticWeightMode == 98 && m_dlbStep % 2 == 0)) {
      m_log << "Using static weights at dlb-step " << m_dlbStep << " :" << std::endl;
      for(MInt i = 0; i < noLoadTypes; i++) {
        weights[i] = m_dlbStaticWeights[i];
        m_log << " " << std::scientific << weights[i];
      }
      m_log << std::endl;
    }
    if(m_dlbStaticWeightMode > 0 && m_dlbStaticWeightMode < 90 && m_dlbStep == m_dlbStaticWeightMode) {
      m_dlbStaticWeightMode = -1;
      m_log << "Switching to computed weights from now on!" << std::endl;
    }
  }


  // store weights and use last weights at restart
  if(m_dlbLastWeights == nullptr) {
    mAlloc(m_dlbLastWeights, noLoadTypes, "dlbLastWeights", AT_);
  }
  for(MInt i = 0; i < noLoadTypes; i++) {
    m_dlbLastWeights[i] = weights[i];
  }
}
/// \brief Solve the parameter estimation problem A*x=b.
///
/// The linear least squares problem A*x=b is regularized with a Tikhonov
/// regularization and solved using Eigen. Parameters x
/// are assumed to be positive. The regularization parameter is increased as long as
/// the residual decreases. Estimated parameters are normalized in the 1-norm and any negative
/// parameter is set to zero.
///
/// \param[in] m Number of rows of the matrix A.
/// \param[in] n Number of columns of the matrix A (x=1 if n==1).
/// \param[in] A Matrix (size m x n).
/// \param[in] b Right hand side vector (size m).
/// \param[out] x Estimated parameter vector (size n).
template <MInt nDim>
void maia::grid::Controller<nDim>::estimateParameters(MInt m, MInt n, const MFloat* const A, const MFloat* const b,
                                                      MFloat* const x) {
  TRACE();

  // Always return the same parameters for testing
  if(m_testDynamicLoadBalancing) {
    // For testing: switch weights in consecutive steps to force a partitioning change
    std::cerr << "estimateParameters: n = " << n << ", setting weights to 1.0" << std::endl;
    std::fill_n(x, n, 1.0);
    if(n > 1 && m_dlbStep % 2 == 0) {
      if(Context::propertyExists("solverWeights_0")) {
        std::vector<MFloat> weights;
        getSpecifiedSolverWeights(weights);
        std::copy(weights.begin(), weights.end(), &x[0]);
        std::cerr << "estimateParameters: using specified solver weights" << std::endl;
      } else {
        std::cerr << "estimateParameters: n = " << n << ", setting w[n-1] = 5.0" << std::endl;
        x[n - 1] = 5.0;
      }
    }
    return;
  }

  if(n == 1) {
    // Single parameter, return 1
    x[0] = 1.0;
  } else {
    // Solve the (overdetermined) system of linear equations (linear least squares)

    // Regularization parameter
    MFloat alpha = 0.01;
    // Maximum number of iterations
    const MInt maxNoIt = 30;

    // Determine column sums
    MFloatScratchSpace sum_j(n, AT_, "sum_j");
    for(MInt j = 0; j < n; j++) {
      sum_j[j] = 0.0;
      for(MInt i = 0; i < m; i++) {
        sum_j[j] += A[i * n + j];
      }
      if(m_debugBalance) {
        m_log << " * estimateParameters regularize " << j << " proportional to " << sum_j[j] << std::endl;
      }
    }

    MFloat oldResidual = 0.0;
    MFloatScratchSpace oldParam(n, AT_, "oldParam");

    MInt it = 0;
    while(it < maxNoIt) {
      // Work on copies of A and b since they get overwritten and additional
      // rows are added to allow a regularization of the overdetermined system
      MFloatScratchSpace A_work((m + n) * n, AT_, "A_work");
      MFloatScratchSpace b_work(m + n, AT_, "b_work");

      // Copy matrix A and fill additional rows with zeros
      std::copy_n(&A[0], m * n, &A_work[0]);
      std::fill_n(&A_work[m * n], n * n, 0.0);

      // Set 'alpha' on diagonal of additional (n x n) solver (Tikhonov regular.)
      for(MInt i = 0; i < n; i++) {
        // Set proportional to average column value or to 1 if column sum is zero
        const MFloat reg = (sum_j[i] > 0) ? alpha * sum_j[i] / noDomains() : 1.0;
        A_work[m * n + i * (n + 1)] = reg;
      }

      // Copy vector b and add zeros for regularization
      std::copy_n(&b[0], m, &b_work[0]);
      std::fill_n(&b_work[m], n, 0.0);

      // Compute the minimum-norm solution to the linear least squares problem:
      // minimize 2-norm(|b - A*x|) using singular value decomposition (SVD) of
      // A. A is an M-by-N matrix which can be rank deficient.

      // TODO labels:CONTROLLER,totest check if this works properly!

      std::vector<MFloat> singularValues = maia::math::svd(A_work.data(), b_work.data(), m + n, n, x);

      // Compute condition number of A (cond_2-norm = S(0)/S(min(m,n)-1))
      // TODO labels:CONTROLLER this is not correct it should be max(S)/min(S)
      const MFloat cond_2 = singularValues[0] / singularValues[std::min(m, n) - 1];

      // Compute residual
      // TODO labels:CONTROLLER could be replaced by
      // maia::math::MFloatVectorXd residual = A_eigen*x_eigen - b_eigen
      MFloat residual = 0.0;
      MFloat maxDiff = 0.0;
      MFloatScratchSpace mxv(m, AT_, "mxv");
      for(MInt i = 0; i < m; i++) {
        mxv[i] = 0.0;
        for(MInt j = 0; j < n; j++) {
          mxv[i] += A[i * n + j] * x[j];
        }
        residual += POW2(mxv[i] - b[i]);
        maxDiff = std::max(maxDiff, mxv[i] - b[i]);
      }

      if(m_debugBalance) {
        m_log << " * estimateParameters iteration " << it << ":";
        for(MInt i = 0; i < n; i++) {
          m_log << std::scientific << " " << x[i];
        }
        m_log << "; residual " << residual << "; maxDiff " << maxDiff << std::endl;
      }

      /* if (std::any_of(&x[0], &x[0] + n, [](MFloat p) { return p < 0.0; }) || cond_2 > maxCond)
       * { */
      /* const MBool hasNegativeParam = std::any_of(&x[0], &x[0] + n, [](MFloat p) { return p <
       * 0.0; }); */

      // If the residual decreased increase regularization constant and repeat
      if(it == 0 || residual <= oldResidual) {
        alpha *= 1.5;
        it++;

        oldResidual = residual;
        std::copy_n(&x[0], n, &oldParam[0]);

        if(m_debugBalance) {
          if(it == maxNoIt) {
            m_log << " * estimateParameters max iterations reached " << it << ", alpha " << alpha
                  << ", condition number " << cond_2 << std::endl;
          } else {
            m_log << " * estimateParameters increase alpha to " << alpha << ", condition number " << cond_2
                  << ", residual " << residual << std::endl;
          }
        }
      } else {
        // Take values of previous iteration
        std::copy_n(&oldParam[0], n, &x[0]);

        if(m_debugBalance) {
          m_log << " * estimateParameters found solution! iteration " << it << ", condition number " << cond_2
                << ", residual " << residual << ", previous " << oldResidual << ", maxDiff " << maxDiff;
          for(MInt i = 0; i < n; i++) {
            m_log << std::scientific << " " << x[i];
          }
          m_log << std::endl;
        }
        break;
      }
    } // alpha loop
  }

  // Prevent negative weights and Nans, i.e. set to zero
  for(MInt i = 0; i < n; i++) {
    if(x[i] < 0.0 || std::isnan(x[i])) {
      x[i] = 0.0;
    }
  }
  // Renormalize parameter vector using 1-norm, such that parameters are better suited for
  // comparison and will add up to one.
  const MFloat norm = std::accumulate(&x[0], &x[0] + n, 0.0);
  for(MInt i = 0; i < n; i++) {
    x[i] /= norm;
  }
}


/// \brief Determine new partitioning for dynamic load balancing.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
///
/// Based on the domain loads and current domain decomposition determine a new
/// partitioning that is expected to reduce imbalances among domains.
///
/// \param[in] loads Vector of domain loads.
/// \param[in] imbalance Load imbalance among domains.
/// \param[out] partitionCellOffsets New partition-cell offsets (size noDomains()+1).
/// \param[out] globalIdOffsets New global domain offsets.
///
/// Available DLB partition methods:
/// 1. DLB_PARTITION_WEIGHT
///   - compute new weights for the different load types
///   - determine the new workload of all min-cells and solve 1D partitioning problem
/// 2. DLB_PARTITION_SHIFT_OFFSETS
///   - compute new weights and domain weights (quantifies processing capacities)
///   - determine summed load errors for each domain offset and shift offsets until the imbalance is
///     predicted to be counterbalanced
/// 3. DLB_PARTITION_DEFAULT
///   - combination of methods 1 and 2: 3 steps with pure weighting, the remaining steps use the
///     offset shifting method
/// 4. DLB_PARTITION_TEST
///   - for testing, use old partitioning to check if solver still works after reinitialization
template <MInt nDim>
MBool maia::grid::Controller<nDim>::loadBalancingPartition(const MFloat* loads, const MFloat imbalance,
                                                           MLong* const partitionCellOffsets,
                                                           MLong* const globalIdOffsets) {
  TRACE();

  using namespace maia;

  m_log << " * determine new partition" << std::endl;

  MLongScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
  MLongScratchSpace localPartitionCellOffsets(noDomains() + 1, AT_, "localPartitionCellOffsets");
  // Determine the number of partition cells on all domains and the partition cell offsets
  gridb().determineNoPartitionCellsAndOffsets(&localPartitionCellCounts[0], &localPartitionCellOffsets[0]);

  const MInt noLocalPartitionCells = localPartitionCellCounts[domainId()];

  // Set the partition method
  MInt partitionMethod = m_dlbPartitionMethod;

  // Change the actual used partition method depending on the DLB step
  switch(m_dlbPartitionMethod) {
    case DLB_PARTITION_DEFAULT: {
      // Default partition method: combination of 3 steps with pure weighting,
      // the remaining steps use the offset shifting method
      const MInt noWeightPartitionSteps = 3;
      if(m_dlbStep < noWeightPartitionSteps) {
        partitionMethod = DLB_PARTITION_WEIGHT;
      } else {
        partitionMethod = DLB_PARTITION_SHIFT_OFFSETS;
      }
      break;
    }
    // Else: dont change the partition method
    case DLB_PARTITION_WEIGHT:
    case DLB_PARTITION_SHIFT_OFFSETS:
    case DLB_PARTITION_TEST:
      break;
    default:
      TERMM(1, "Unknown DLB partition method.");
      break;
  }

  // Initialize domain weights in first DLB step
  if(m_dlbStep == 0) {
    m_domainWeights.assign(noDomains(), 1.0);
    m_lastOffsetShiftDirection.assign(noDomains() + 1, 0);
  }

  // Increase DLB step
  m_dlbStep++;

  // Compute new partitioning with the activated partition method
  switch(partitionMethod) {
    case DLB_PARTITION_WEIGHT: {
      m_log << "Partition method #0 (compute weights, determine min cells workload and use the "
               "partition() algorithm)"
            << std::endl;
      // Use load quantities, i.e. number of DOF or number of active cells, and loads to estimate
      // new weights. The new domain distribution is then determined by the partition() algorithm.
      std::vector<MFloat> weights;
      computeWeights(loads, 1.0, weights);

      updateWeightsAndWorkloads(weights, false);

      // Note: in case of a partition level shift the global id offsets might not be correct here,
      // i.e., they are not corrected. This is done later in loadBalancingCalcNewGlobalOffsets().
      partition(&partitionCellOffsets[0], &globalIdOffsets[0], true);
      break;
    }
    case DLB_PARTITION_TEST: {
      // TESTING: use old offsets to check if solver still works after reinit
      for(MInt i = 0; i < noDomains() + 1; i++) {
        partitionCellOffsets[i] = localPartitionCellOffsets[i];
      }
      if(domainId() == 0) {
        std::cerr << "Using original partition!" << std::endl;
      }
      break;
    }
    case DLB_PARTITION_SHIFT_OFFSETS: {
      m_log << "Partition method #3 (compute weights and domain weights iteratively, determine "
               "summed load error at each domain offset and vary min cell offset such that "
               "predicted error is zero)"
            << std::endl;

      std::vector<MFloat> weights;
      std::vector<MFloat> oldWeights;
      std::vector<MFloat> oldDomainWeights;
      oldDomainWeights = m_domainWeights;
      std::fill(oldDomainWeights.begin(), oldDomainWeights.end(),
                1.0); // TODO labels:CONTROLLER,DLB reset in each step?

      MFloatScratchSpace domainWorkLoad(noDomains(), AT_, "domainWorkLoad");
      MFloat averageWorkLoad = -1.0;

      MFloatScratchSpace localPartitionCellsWorkload(noLocalPartitionCells, AT_, "localPartitionCellsWorkload");
      MFloatScratchSpace partitionCellsWorkload(gridb().m_noPartitionCellsGlobal, AT_, "partitionCellsWorkload");

      const MInt maxNoIt = 20;
      // Iterate until new weights and domain weights are found, i.e. converged
      for(MInt it = 0; it < maxNoIt; it++) {
        // Compute weights
        computeWeights(loads, m_domainWeights[domainId()], weights);

        // @ansgar_mb TODO labels:CONTROLLER,DLB test this
        updateWeightsAndWorkloads(weights, false);

        // Assemble local partition-cell workloads
        for(MInt i = 0; i < noLocalPartitionCells; i++) {
          const MLong globalPartitionCellId = gridb().m_localPartitionCellGlobalIds[i];
          const MInt partitionCellId = gridb().globalIdToLocalId(globalPartitionCellId, true);
          TERMM_IF_NOT_COND(gridb().a_hasProperty(partitionCellId, Cell::IsPartitionCell),
                            "Error: cell is not a partition cell.");

          const MFloat workload = gridb().a_workload(partitionCellId);
          TERMM_IF_NOT_COND(workload > 0.0, "Error: partition cell workload needs to be > 0.0");
          localPartitionCellsWorkload[i] = workload;
        }

        // @ansgar TODO labels:DLB Temporary fix; use MLong!? only an issue for more than 2billion partition cells...
        MIntScratchSpace localPartitionCellCounts_(noDomains(), AT_, "localPartitionCellCounts_");
        MIntScratchSpace localPartitionCellOffsets_(noDomains() + 1, AT_, "localPartitionCellOffsets_");
        for(MInt i = 0; i < noDomains(); i++) {
          localPartitionCellCounts_[i] = (MInt)localPartitionCellCounts[i];
          localPartitionCellOffsets_[i] = (MInt)localPartitionCellOffsets[i];
        }
        localPartitionCellOffsets_[noDomains()] = localPartitionCellOffsets[noDomains()];

        // Gather partition-cell workloads on root
        MPI_Gatherv(&localPartitionCellsWorkload[0], noLocalPartitionCells, type_traits<MFloat>::mpiType(),
                    &partitionCellsWorkload[0], &localPartitionCellCounts_[0], &localPartitionCellOffsets_[0],
                    type_traits<MFloat>::mpiType(), 0, gridb().mpiComm(), AT_, "localPartitionCellsWorkload[0]",
                    "partitionCellsWorkload[0]");
        MInt error = 0;
        if(domainId() == 0) {
          // Output some information on partition cell workloads
          const MFloat maxPartitionWorkload =
              *std::max_element(partitionCellsWorkload.begin(), partitionCellsWorkload.end());
          const MFloat minPartitionWorkload =
              *std::min_element(partitionCellsWorkload.begin(), partitionCellsWorkload.end());
          const MFloat avgPartitionWorkload =
              std::accumulate(partitionCellsWorkload.begin(), partitionCellsWorkload.end(), 0.0)
              / static_cast<MFloat>(gridb().m_noPartitionCellsGlobal);
          m_log << " * Maximum/minimum/average partition cell workload: " << maxPartitionWorkload << ", "
                << minPartitionWorkload << ", " << avgPartitionWorkload << std::endl;
          if(minPartitionWorkload <= 0) {
            m_log << globalTimeStep << " ERROR: minimum partition cell workload is " << minPartitionWorkload
                  << std::endl;
            cerr0 << globalTimeStep << " ERROR: minimum partition cell workload is " << minPartitionWorkload
                  << std::endl;
            error = 1;
          }

          // Compute domain workloads
          for(MInt i = 0; i < noDomains(); i++) {
            const MInt offset = localPartitionCellOffsets[i];
            const MInt count = localPartitionCellCounts[i];
            domainWorkLoad[i] =
                std::accumulate(&partitionCellsWorkload[offset], &partitionCellsWorkload[offset] + count, 0.0);
          }
          // Compute average workload
          averageWorkLoad = std::accumulate(domainWorkLoad.begin(), domainWorkLoad.end(), 0.0) / noDomains();

          m_log << " * Average domain workload: " << averageWorkLoad << std::endl;

          // Update domain weights
          for(MInt i = 0; i < noDomains(); i++) {
            m_domainWeights[i] = (m_useDomainFactor) ? loads[i] * averageWorkLoad / domainWorkLoad[i] : 1.0;
          }

          const MFloat domainWeightMinLoadThreshold = 0.2;

          MFloat averageDomainWeight = 0.0;
          MInt countDomainWeights = 0;
          for(MInt i = 0; i < noDomains(); i++) {
            if(loads[i] >= domainWeightMinLoadThreshold) {
              averageDomainWeight += m_domainWeights[i];
              countDomainWeights++;
            }
          }
          averageDomainWeight /= (MFloat)countDomainWeights;

          for(MInt i = 0; i < noDomains(); i++) {
            // For domains with low load reset domain weight since estimation might fail
            if(loads[i] < domainWeightMinLoadThreshold) {
              m_domainWeights[i] = 1.0;
            } else {
              // Scale domain weights with average such that the mean domain weight is one
              m_domainWeights[i] /= averageDomainWeight;

              // Weighting with previous values
              m_domainWeights[i] = 0.5 * (m_domainWeights[i] + oldDomainWeights[i]);
            }
          }

          const MFloat maxDomainWeight = *std::max_element(m_domainWeights.begin(), m_domainWeights.end());
          const MFloat minDomainWeight = *std::min_element(m_domainWeights.begin(), m_domainWeights.end());
          m_log << " * Domain weights: max = " << maxDomainWeight << ", min = " << minDomainWeight << std::endl;
        }

        // Broadcast error flag
        MPI_Bcast(&error, 1, type_traits<MInt>::mpiType(), 0, gridb().mpiComm(), AT_, "error");
        if(error > 0) {
          return false;
        }

        // Broadcast new domain weights
        MPI_Bcast(&m_domainWeights[0], noDomains(), type_traits<MFloat>::mpiType(), 0, gridb().mpiComm(), AT_,
                  "m_domainWeights[0]");

        // Compute residuals and check for convergence
        if(it > 0) {
          const MInt noWeights = weights.size();
          MFloat weightResidual = 0.0;
          for(MInt i = 0; i < noWeights; i++) {
            weightResidual += ABS(weights[i] - oldWeights[i]);
          }
          weightResidual /= noWeights;

          MFloat domainWeightResidual = 0.0;
          for(MInt i = 0; i < noDomains(); i++) {
            domainWeightResidual += ABS(m_domainWeights[i] - oldDomainWeights[i]);
          }
          domainWeightResidual /= noDomains();

          const MFloat weightResThreshold = 1e-3;
          const MFloat domainWeightResThreshold = 1e-2;
          if(weightResidual < weightResThreshold && domainWeightResidual < domainWeightResThreshold) {
            m_log << " * weights and domain weights converged, residuals: " << weightResidual << ", "
                  << domainWeightResidual << std::endl;
            break;
          } else if(it == maxNoIt - 1) {
            m_log << " * WARNING: Weights and domain weights not converged, but "
                     "maximum number of iterations reached. "
                  << weightResidual << " " << domainWeightResidual << std::endl;
          } else {
            m_log << " * Weights and domain weights iteration " << it << " " << weightResidual << " "
                  << domainWeightResidual << std::endl;
          }
        }
        oldWeights = weights;
        oldDomainWeights = m_domainWeights;
      }

      if(m_debugBalance) {
        const MInt noLoadTypes = globalNoLoadTypes();
        MIntScratchSpace loadQuantities(noLoadTypes, AT_, "loadQuantities");
        getLoadQuantities(&loadQuantities[0]);
        storeLoadsAndWeights(&loads[0], noLoadTypes, &loadQuantities[0], m_domainWeights[domainId()], &weights[0]);
      }

      // New Weights and domain weights are computed, rank 0 has already the partition
      // cell workload information, now determine the new partition cell offsets
      if(domainId() == 0) {
        MFloatScratchSpace errors(noDomains(), AT_, "errors");
        MFloatScratchSpace summedErrors(noDomains() + 1, AT_, "summedErrors");
        MFloatScratchSpace summedErrAbs(noDomains() + 1, AT_, "summedErrAbs");
        MIntScratchSpace summedErrDir(noDomains() + 1, AT_, "summedErrDir");

        // Compute load error for each domain ...
        for(MInt i = 0; i < noDomains(); i++) {
          errors[i] = loads[i] - 1.0;
        }
        // ... and summed load errors (and their absolut value) at each domain offset
        summedErrors[0] = 0.0;
        summedErrAbs[0] = 0.0;
        summedErrDir[0] = 0;
        for(MInt i = 1; i < noDomains() + 1; i++) {
          summedErrors[i] = summedErrors[i - 1] + errors[i - 1];
          summedErrAbs[i] = ABS(summedErrors[i]);
          summedErrDir[i] = -signum(summedErrors[i]);
        }

        const MFloat maxSummedErr = *std::max_element(summedErrors.begin(), summedErrors.end());
        const MFloat minSummedErr = *std::min_element(summedErrors.begin(), summedErrors.end());
        m_log << globalTimeStep << " * Summed load deviations: max = " << maxSummedErr << ", min = " << minSummedErr
              << std::endl;

        // Threshold below which an offset will not be changed
        const MFloat errorThreshold = m_dlbImbalanceThreshold;

        // Restrict the amount of imbalance to be shifted in a single step for higher degrees of
        // parallelism
        const MFloat restrictionFactor = 1.0 - std::min(0.5, 0.005 * noDomains());

        MFloat loadVariance = 0.0;
        for(MInt i = 0; i < noDomains(); i++) {
          loadVariance += POW2(loads[i] - 1.0);
        }
        loadVariance /= noDomains();
        const MFloat loadStd = std::sqrt(loadVariance);
        m_log << globalTimeStep << " * load distribution: variance = " << loadVariance << ", stdev = " << loadStd
              << std::endl;

        const MFloat imbalancePenaltyThreshold = 2.0 * errorThreshold;
        const MFloat imbPenaltyFactor =
            (imbalance < imbalancePenaltyThreshold)
                ? (1.0 + (imbalancePenaltyThreshold - imbalance) / imbalancePenaltyThreshold)
                : 1.0;
        m_log << globalTimeStep << " * Imbalance penalty factor " << imbPenaltyFactor << " (imbalance " << imbalance
              << "); restriction factor " << restrictionFactor << std::endl;

        // Final local shift for refining the best found partitioning
        const MBool finalLocalShift =
            m_dlbNoFinalLocalShifts > 0
            && globalTimeStep >= m_loadBalancingStopTimeStep - m_loadBalancingInterval * (m_dlbNoFinalLocalShifts + 1);
        // Intermediate local shift
        const MBool intermediateLocalShift =
            !finalLocalShift && (m_dlbNoLocalShifts > 0 && (m_dlbStep % (1 + m_dlbNoLocalShifts) != 1));
        // Check for global shift step
        const MBool globalShiftStep = !finalLocalShift && !intermediateLocalShift;

        m_log << globalTimeStep << " * Shift step: global " << globalShiftStep << "; intermediate local "
              << intermediateLocalShift << "; final local " << finalLocalShift << std::endl;

        // Flag to indicate if a global shift is required for each offset
        ScratchSpace<MInt> globalShiftFlag(noDomains() + 1, AT_, "globalShiftFlag");
        globalShiftFlag.fill(0);
        // Minimum distance of an offset to the next offset with a global shift
        ScratchSpace<MInt> minGlobalDist(noDomains() + 1, AT_, "minGlobalDist");
        minGlobalDist.fill(noDomains());

        { // Determine offsets with a global shift
          ScratchSpace<MInt> globalShiftFlagTmp(noDomains() + 1, AT_, "globalShiftFlagTmp");
          globalShiftFlagTmp.fill(0);

          MInt noGlobalShifts = 0;
          // 1. Determine offsets which require a global shift, i.e., the load is too imbalanced
          // among the neighboring domains
          for(MInt i = 1; i < noDomains(); i++) {
            const MInt dist = 4; // consider offsets in this distance to both sides
            const MInt firstOffset = std::max(0, i - dist);
            const MInt lastOffset = std::min(noDomains(), i + dist);
            // Difference in summed imbalance
            const MFloat summedErrorDiff = ABS(summedErrors[lastOffset] - summedErrors[firstOffset]);
            const MBool sameSign = (summedErrDir[lastOffset] == summedErrDir[firstOffset]);
            const MBool absErrCondition = (summedErrAbs[i] > errorThreshold);

            const MFloat smoothShiftThreshold = (lastOffset - firstOffset) * errorThreshold; // 0.5
            const MBool smoothShiftCondition =
                (absErrCondition && (sameSign && summedErrorDiff > smoothShiftThreshold));
            // Global shift required if the difference in the summed imbalance is too high
            MBool globalShiftCondition = false;
            if(globalShiftStep) {
              globalShiftCondition = (m_dlbSmoothGlobalShifts) ? smoothShiftCondition : absErrCondition;
            }

            globalShiftFlagTmp[i] = (globalShiftCondition) ? summedErrDir[i] : 0;
            if(globalShiftFlagTmp[i] != 0) {
              noGlobalShifts++;
              if(m_debugBalance) {
                m_log << globalTimeStep << " * global shift " << i << " diff=" << summedErrorDiff
                      << " tr=" << smoothShiftThreshold << " s_first=" << summedErrors[firstOffset]
                      << " s_last=" << summedErrors[lastOffset] << " s_abs=" << summedErrAbs[i] << std::endl;
              }
            }
          }
          m_log << globalTimeStep << " * number of global shifts " << noGlobalShifts << std::endl;

          // 2. For each offset marked with a global shift, mark also the nearby offsets up to a
          // certain distance while reducing the size of the shift such that global shifts are
          // smoothed out over some neighborhood. For this the minimum distance to a global shift is
          // stored, which later reduces the summed imbalance that needs to be counterbalanced by
          // the shift.
          for(MInt i = 1; i < noDomains(); i++) {
            // Check for offset with global shift
            if(globalShiftFlagTmp[i] != 0) {
              minGlobalDist[i] = 0;
              globalShiftFlag[i] = globalShiftFlagTmp[i];

              // Distance over which to smooth out the global shift
              const MInt dist = std::ceil(summedErrAbs[i] / errorThreshold);
              if(m_debugBalance) {
                m_log << globalTimeStep << " * global shift " << i << " marking distance " << dist << " "
                      << summedErrAbs[i] << std::endl;
              }

              // Loop over lower offsets
              for(MInt nId = i - 1; nId >= std::max(i - dist, 1); nId--) {
                const MInt distance = ABS(nId - i);

                // Stop marking offsets if shift direction is not the same anymore or if the
                // summed imbalance is too small for this distance, which would result in a summed
                // imbalance with opposite sign
                if(summedErrDir[nId] != summedErrDir[i] || errorThreshold * distance > summedErrAbs[nId]) {
                  break;
                }

                minGlobalDist[nId] = std::min(minGlobalDist[nId], distance);
                globalShiftFlag[nId] = summedErrDir[nId];
              }
              // Loop over higher offsets
              for(MInt nId = i + 1; nId <= std::min(i + dist, noDomains() - 1); nId++) {
                const MInt distance = ABS(nId - i);

                if(summedErrDir[nId] != summedErrDir[i] || errorThreshold * distance > summedErrAbs[nId]) {
                  break;
                }

                minGlobalDist[nId] = std::min(minGlobalDist[nId], distance);
                globalShiftFlag[nId] = summedErrDir[nId];
              }
            }
          }

          MInt noGlobalShiftsTotal = 0;
          for(MInt i = 1; i < noDomains(); i++) {
            if(globalShiftFlag[i] != 0) {
              noGlobalShiftsTotal++;
            }
          }
          m_log << globalTimeStep << " * number of total global shifts " << noGlobalShiftsTotal << std::endl;
        }

        // Variables for some user output about the shifts
        MInt noShifts = 0;
        MInt noShiftsWithoutChange = 0;
        MFloat maxWeight = 0.0;
        MFloat maxDiff = 0.0;
        // Storage for shift information that is written to file in debug mode
        ScratchSpace<MInt> shifts(noDomains() + 1, AT_, "shifts");
        shifts.fill(0);
        ScratchSpace<MFloat> shiftedWorkload(noDomains() + 1, AT_, "shiftedWorkload");
        shiftedWorkload.fill(0.0);

        partitionCellOffsets[0] = 0;
        // vary min cell offsets individually
        for(MInt offsetId = 1; offsetId < noDomains(); offsetId++) {
          MFloat summedWorkload = 0.0;
          // Determine current workload of the domain left to this offset (with new previous offset)
          for(MInt i = partitionCellOffsets[offsetId - 1]; i < localPartitionCellOffsets[offsetId]; i++) {
            summedWorkload += (m_useDomainFactor) ? m_domainWeights[offsetId - 1] * partitionCellsWorkload[i]
                                                  : partitionCellsWorkload[i];
          }

          // Check for a global or local shift
          const MBool globalShiftCondition = (globalShiftFlag[offsetId] != 0);
          MBool localShiftCondition = false;
          MInt localShiftDir = 0;
          MFloat localShiftDiff = 0.0;

          // No global shift for this offset and the neighboring ones, use a local shift to
          // distribute load among neighboring domains
          if(!globalShiftCondition && globalShiftFlag[offsetId - 1] == 0 && globalShiftFlag[offsetId + 1] == 0) {
            const MFloat errLeft = errors[offsetId - 1];
            const MFloat errRight = errors[offsetId];
            localShiftDiff = 0.5 * (errLeft - errRight);

            if(errLeft > errorThreshold && errRight < errLeft) {
              // overload left
              localShiftDir = -1;
            } else if(errRight > errorThreshold && errRight > errLeft) {
              // overload right
              localShiftDir = 1;
            } else if(errLeft < -errorThreshold && errLeft < errRight) {
              // underload left (l/r no overload)
              localShiftDir = 1;
            } else if(errRight < -errorThreshold && errRight < errLeft) {
              // underload right (l/r no overload)
              localShiftDir = -1;
            }
            localShiftCondition = (localShiftDir != 0);
          }

          const MBool shift = globalShiftCondition || localShiftCondition;

          if(m_debugBalance) {
            m_log << globalTimeStep << " * shift conditions " << offsetId << " shift=" << shift
                  << " g=" << globalShiftCondition << " l=" << localShiftCondition << " err=" << summedErrors[offsetId]
                  << " minDist=" << minGlobalDist[offsetId] << " lDir=" << localShiftDir << " diff=" << localShiftDiff
                  << " eL=" << errors[offsetId - 1] << " eR=" << errors[offsetId] << std::endl;
          }

          if(shift) {
            // Set shift direction depending on global or local shift
            const MInt dir = (globalShiftCondition) ? summedErrDir[offsetId] : localShiftDir;
            MFloat penaltyFactor = 1.25;

            const MInt lastDir = (globalShiftCondition) ? m_lastOffsetShiftDirection[offsetId] : 0;
            // Check for current shift in opposite direction of last shift (only if shift is global)
            if(lastDir != dir && lastDir != 0) {
              penaltyFactor = 2.0;
              if(m_debugBalance) {
                m_log << "Opposite shift " << offsetId << ", lastDir " << lastDir << ", dir " << dir << std::endl;
              }
            }
            penaltyFactor *= imbPenaltyFactor;

            // Determine 'new' and 'old' domain id depending on the direction in
            // which the partition cell offset will be shifted
            const MInt domainOld = (dir < 0) ? offsetId - 1 : offsetId;
            const MInt domainNew = (dir < 0) ? offsetId : offsetId - 1;

            // Initialize estimated summed load error
            // To smooth out the globally enforced shifts, the error threshold times
            // the minimum distance to the next domain which requires a global shift is subtracted,
            // i.e., the shift of an offset further away from a global shift is reduced
            const MFloat summedErrorTmp =
                restrictionFactor * (summedErrors[offsetId] + dir * errorThreshold * minGlobalDist[offsetId]);
            //(dir > 0) ? std::max(-2.0, summedErrorTmp) : std::min(2.0, summedErrorTmp);
            const MFloat summedError = summedErrorTmp;
            MFloat oldEstimate = (globalShiftCondition) ? summedError : localShiftDiff;
            const MFloat initialEstimate = oldEstimate;

            MFloat penalty = 1.0;

            // Domain id of the currently considered partition cell
            MInt partitionCellDomain = domainOld;

            if(m_debugBalance) {
              m_log << "DLB_DEBUG: Start offset shift " << offsetId << ", domainOld " << domainOld << ", domainNew "
                    << domainNew << ", oldEstimate " << oldEstimate << ", summedError " << summedErrors[offsetId]
                    << ", penaltyFactor " << penaltyFactor << std::endl;
            }

            // Shift the partition cell offset such that the new estimated summed load error is
            // minimized
            for(MInt i = 1;; i++) {
              // TODO labels:DLB check if first/last partition cell is reached and stop!
              // Currently considered partition cell id
              const MInt newPartitionCellId = localPartitionCellOffsets[offsetId] + dir * i;

              // Update domain id on which the current partition cell belongs
              if(dir > 0 && newPartitionCellId == localPartitionCellOffsets[partitionCellDomain + 1]) {
                // Next partition cell offsets reached, increase domain id
                partitionCellDomain++;
              } else if(dir < 0 && newPartitionCellId == localPartitionCellOffsets[partitionCellDomain] - 1) {
                // Previous partition cell offset exceeded, decrease domain id
                partitionCellDomain--;
              }

              const MBool domainFactorShift = true; //(oldEstimate < 1.0); // TODO labels:CONTROLLER,DLB
              // Performance factor of domains
              const MFloat domainFactor = (m_useDomainFactor && domainFactorShift)
                                              ? (m_domainWeights[domainNew] / m_domainWeights[partitionCellDomain])
                                              : 1.0;

              // Compute new estimate of summed load error
              const MFloat estimate = oldEstimate
                                      + dir * penaltyFactor * penalty * domainFactor * loads[partitionCellDomain]
                                            * partitionCellsWorkload[newPartitionCellId]
                                            / domainWorkLoad[partitionCellDomain];
              const MFloat diff = estimate - oldEstimate;

              // Reset summed workload if the previous new offset is still larger than the current
              // one
              if(newPartitionCellId <= partitionCellOffsets[offsetId - 1]) {
                summedWorkload = 0;
              }
              const MFloat workloadIncrement =
                  (m_useDomainFactor) ? m_domainWeights[offsetId - 1] * partitionCellsWorkload[newPartitionCellId]
                                      : partitionCellsWorkload[newPartitionCellId];
              const MFloat summedWorkloadNew = summedWorkload + dir * workloadIncrement;
              const MFloat summedWorkloadRel = summedWorkload / averageWorkLoad;

#ifndef NDEBUG
              // Only useful for debugging small cases
              if(m_debugBalance) {
                m_log << "DLB_DEBUG: offset " << offsetId << ", partCellId " << newPartitionCellId
                      << ", partCellDomain " << partitionCellDomain << ", domainFactor " << domainFactor
                      << ", oldEstimate " << oldEstimate << ", estimate " << estimate << ", diff " << diff
                      << ", penalty " << penalty << ", partCellWorkload " << partitionCellsWorkload[newPartitionCellId]
                      << ", domainWorkload " << domainWorkLoad[partitionCellDomain] << ", load "
                      << loads[partitionCellDomain] << std::endl;
              }
#endif

              // Accept the previous new partition cell offset if the predicted summed load error
              // reached its minimum or if the previous domain has only a single parititon cell left
              // (when decreasing the current partition cell offset).
              const MBool prevOffsetReached = (dir < 0 && newPartitionCellId == partitionCellOffsets[offsetId - 1]);
              const MBool lastCellReached = (newPartitionCellId == gridb().m_noPartitionCellsGlobal - 1);

              // Continue shift if the workload is too large
              const MBool workloadCheck = (m_dlbMaxWorkloadLimit > 1.0 && globalShiftCondition)
                                              ? (summedWorkloadRel < m_dlbMaxWorkloadLimit)
                                              : true;

              if(prevOffsetReached || lastCellReached || (ABS(estimate) >= ABS(oldEstimate) && workloadCheck)
                 || (dir == 1 && !workloadCheck)) {
                // New domain offset
                partitionCellOffsets[offsetId] = localPartitionCellOffsets[offsetId] + dir * (i - 1);

                if(m_debugBalance) {
                  m_log << globalTimeStep << " * shifted offset " << offsetId << " " << dir * (i - 1) << " "
                        << localPartitionCellOffsets[offsetId] << " " << partitionCellOffsets[offsetId] << " "
                        << summedErrors[offsetId] << " " << estimate << " " << oldEstimate << " " << penalty << " "
                        << summedWorkloadRel << std::endl;
                }
                shifts[offsetId] = dir * (i - 1);

                // Set direction of shift (if it was a global shift)
                m_lastOffsetShiftDirection[offsetId] = (i > 1 && globalShiftCondition) ? dir : 0;

                // Make sure the current domain has at least a single min cell
                // TODO labels:DLB handle last domain offset!
                if(partitionCellOffsets[offsetId] <= partitionCellOffsets[offsetId - 1]) {
                  partitionCellOffsets[offsetId] = partitionCellOffsets[offsetId - 1] + 1;
                  shifts[offsetId] = partitionCellOffsets[offsetId] - localPartitionCellOffsets[offsetId];
                  m_log << "WARNING: setting domain offset #" << offsetId
                        << " with just a single partition cell on this domain. " << partitionCellOffsets[offsetId - 1]
                        << " " << partitionCellOffsets[offsetId] << std::endl;
                } else if(partitionCellOffsets[offsetId] == partitionCellOffsets[offsetId - 1] + 1) {
                  m_log << "WARNING: domain #" << offsetId - 1 << " has just a single partition cell. "
                        << partitionCellOffsets[offsetId - 1] << " " << partitionCellOffsets[offsetId] << std::endl;
                }

                if(shifts[offsetId] != 0) {
                  noShifts++;
                } else {
                  noShiftsWithoutChange++;
                }

                // Continue with next domain offset
                break;
              } else {
                // Update estimate and continue with next partition cell
                oldEstimate = estimate;
                // Update summed workload for the domain left to the offset
                summedWorkload = summedWorkloadNew;

                // Determine max difference and weight for some user output
                maxDiff = std::max(maxDiff, ABS(diff));
                maxWeight = std::max(maxWeight, partitionCellsWorkload[newPartitionCellId]);
                // Increase shifted workload
                shiftedWorkload[offsetId] += dir * partitionCellsWorkload[newPartitionCellId];

                // Compute new penalty factor: 1.0/cos(pi/2*((s-e)/s))
                // Goes to infinity as the estimate goes to zero and is intended to prevent
                // 'overshooting' of the domain offsets, i.e. shifting the domain offsets too far
                const MFloat s = initialEstimate;
                const MFloat x = ABS((ABS(s) - ABS(estimate)) / ABS(s));
                penalty = 1.0 / cos(M_PI / 2.0 * x);
              }
            }
          } else {
            // Accept old partition cell offset
            partitionCellOffsets[offsetId] = localPartitionCellOffsets[offsetId];
            m_lastOffsetShiftDirection[offsetId] = 0;
            shifts[offsetId] = 0;

            if(m_debugBalance) {
              m_log << globalTimeStep << " * keep offset " << offsetId << " " << summedErrors[offsetId] << std::endl;
            }
          }
        }

        // Correct offsets which were not shifted to prevent domains with an invalid number of cells!
        for(MInt offsetId = 1; offsetId < noDomains(); offsetId++) {
          if(partitionCellOffsets[offsetId] <= partitionCellOffsets[offsetId - 1]) {
            partitionCellOffsets[offsetId] = partitionCellOffsets[offsetId - 1] + 1;
            shifts[offsetId] = partitionCellOffsets[offsetId] - localPartitionCellOffsets[offsetId];
            m_log << "WARNING: correcting domain offset #" << offsetId
                  << " with now just a single partition cell on this domain. " << partitionCellOffsets[offsetId - 1]
                  << " " << partitionCellOffsets[offsetId] << std::endl;
          }
        }

        m_log << " * number of shifted offsets: " << noShifts
              << "; offsets to shift without change: " << noShiftsWithoutChange << std::endl;
        m_log << " * max weight during shifting: " << maxWeight << "; max difference: " << maxDiff << std::endl;

        // Debug: store partition offset shifts to file
        if(m_debugBalance) {
          FILE* shiftsFile = nullptr;
          const MString fileName = "shifts_" + std::to_string(globalTimeStep) + ".dat";
          shiftsFile = fopen(fileName.c_str(), "w");

          fprintf(shiftsFile, "# 1:offsetId 2:partitionCellShift 3:shiftedWorkload "
                              "4:summedImbalance 5:load 6:domainWorkload 7:domainWeight\n");
          for(MInt i = 0; i < noDomains(); i++) {
            fprintf(shiftsFile, "%d %d %.8e %.8e %.8e %.8e %.8e\n", i, shifts[i], shiftedWorkload[i], summedErrors[i],
                    loads[i], domainWorkLoad[i], m_domainWeights[i]);
          }
          fclose(shiftsFile);
        }
      }

      // Broadcast new partition cell offsets
      MPI_Bcast(&partitionCellOffsets[0], noDomains(), type_traits<MLong>::mpiType(), 0, gridb().mpiComm(), AT_,
                "partitionCellOffsets[0]");
      partitionCellOffsets[noDomains()] = gridb().m_noPartitionCellsGlobal;

      break;
    }
    default: {
      TERMM(1, "Partition method not known.");
      break;
    }
  }

  // Check if new distribution is different to the old one
  const MBool newPartition =
      !std::equal(&localPartitionCellOffsets[0], &localPartitionCellOffsets[0] + noDomains(), &partitionCellOffsets[0]);

  m_log << " * checking new partition" << std::endl;
  MBool validPartition =
      (partitionCellOffsets[0] == 0) && (partitionCellOffsets[noDomains() - 1] < gridb().m_noPartitionCellsGlobal);
  if(!validPartition) {
    m_log << "   * invalid first/last domain offset:" << partitionCellOffsets[0] << " "
          << partitionCellOffsets[noDomains() - 1] << "; noPartitionCellsGlobal = " << gridb().m_noPartitionCellsGlobal
          << std::endl;
  }
  // Check if new partitioning is valid, i.e. the domain offsets are in ascending order without duplicates
  for(MInt i = 1; i < noDomains(); i++) {
    if(partitionCellOffsets[i] <= partitionCellOffsets[i - 1]) {
      validPartition = false;
      m_log << "   * invalid offsets " << i - 1 << " " << partitionCellOffsets[i - 1] << " and " << i << " "
            << partitionCellOffsets[i] << std::endl;
    }
  }

  // Return if partitioning is not valid
  // Note: this may not work if m_forceLoadBalancing is active
  if(!validPartition) {
    m_log << " * new partition is not valid!" << std::endl;
    return false;
  }

  // Return if partitioning did not change and DLB is not forced
  if(!newPartition && !m_forceLoadBalancing) {
    m_log << " * partitioning did not change, return." << std::endl;
    return false;
  } else {
    m_log << " * new partition (" << newPartition << "), forced balance (" << m_forceLoadBalancing << ")" << std::endl;
  }

  // Calculate new global domain offsets from the new partition cell offsets
  loadBalancingCalcNewGlobalOffsets(&localPartitionCellOffsets[0], &partitionCellOffsets[0], &globalIdOffsets[0]);

  return true;
}

template class maia::grid::Controller<2>;
template class maia::grid::Controller<3>;
