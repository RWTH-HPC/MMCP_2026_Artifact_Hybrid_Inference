// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "couplerfvmultilevel.h"

#include <cmath>
#include "FV/fvcartesiansolverxd.h"
#include "MEMORY/scratch.h"

using namespace std;

/// Constructor only stores arguments. All non-trivial initialization is performed in init().
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-10-01
template <MInt nDim, class SysEqn>
CouplerFvMultilevel<nDim, SysEqn>::CouplerFvMultilevel(
    const MInt couplingId, std::vector<FvCartesianSolverXD<nDim, SysEqn>*> solvers)
  : Base(couplingId), BaseFv(couplingId, solvers, solvers.size()) {
  TRACE();

  // Initialize timers
  initTimers();
}


template <MInt nDim, class SysEqn>
CouplerFvMultilevel<nDim, SysEqn>::~CouplerFvMultilevel() {
  TRACE();

  // Stop coupler timer
  RECORD_TIMER_STOP(timer("Total"));
}


/// \brief Initialize the FV-multilevel coupler
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::init() {
  TRACE();

  /*! \page propertyPage1
    \section prolongationMethod
    <code>MInt CouplerFvMultilevel::m_prolongationMethod</code> \n
    possible values: \n
    0: simple injection (first-order accuracy)
    1: bilinear interpolation (second-order accuracy)
    2: slope interpolation (second-order accuracy)
    default = <code>0 </code> \n \n
    Possible values are:
    <ul>
    <li><code>false</code> (off)</li>
    <li><code>true</code> (on)</li>
    </ul>
    Keywords: <i>FV, MULTILEVEL, COUPLING</i>
  */
  m_prolongationMethod = 0;
  m_prolongationMethod = Context::getBasicProperty("prolongationMethod", AT_, &m_prolongationMethod);

  m_log << "Multilevel calculate gradients: " << m_prolongationMethod << std::endl;

  /*! \page propertyPage1
    \section correctCoarseBndry
    <code>MBool CouplerFvMultilevel::m_correctCoarseBndry</code> \n
    default = <code>false</code> \n \n
    Possible values are:
    <ul>
    <li><code>false</code> (off)</li>
    <li><code>true</code> (on)</li>
    </ul>
    Keywords: <i>FV, MULTILEVEL, COUPLING</i>
  */
  m_correctCoarseBndry = false;
  m_correctCoarseBndry = Context::getBasicProperty("correctCoarseBndry", AT_, &m_correctCoarseBndry);

  // Mark solver 0 as primary multilevel solver
  fvSolver(0).setMultilevelPrimary();

  initRestrictedCells();

  initMapping();

  initLeafCells();

  resetMultiLevelVariables();

  if(!fvSolver(0).m_adaptation || (globalTimeStep > 0 && !fvSolver(0).m_levelSetMb)) {
    resetSecondaryFlowVariables();
  }

  // Perform checks
  sanityCheck(0);
}

/// \brief call after the initial adaptation when all cells are refined correctly
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::finalizeCouplerInit() {
  // for initial adaptation, check again
  if(fvSolver(0).m_adaptation && globalTimeStep == 0) {
    sanityCheck(1);
  }

  if(m_correctCoarseBndry && !fvSolver(0).m_levelSetMb) {
    for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
      correctCoarseBndryCells(solverId);
    }
  }

  const auto nan = std::numeric_limits<MFloat>::quiet_NaN();
  auto& fine = fvSolver(0);
  for(MInt cellId = 0; cellId < fine.a_noCells(); cellId++) {
    if(fine.a_isInactive(cellId)) {
      for(MInt varId = 0; varId < fine.CV->noVariables; varId++) {
        fine.a_variable(cellId, varId) = nan;
      }
    }
  }

  initSplitMapping();
}

/// \brief call after the initial adaptation when all cells are refined correctly
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::finalizeAdaptation(MInt solverId) {
  if(solverId == fvSolver(noSolvers() - 1).solverId()) {
    // after all solvers finished their finalization

    m_restrictedCells.clear();
    m_coarse2fine.clear();
    m_fine2coarse.clear();
    m_leafCells.clear();
    m_splitCellMapping.clear();

    initRestrictedCells();

    initMapping();

    initLeafCells();

    if(globalTimeStep < 0) {
      resetMultiLevelVariables();

      resetSecondaryFlowVariables();
    }
  }
}


template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::preCouple(MInt solverId) {
  // apply coarse boundary correction to all secondary solvers after the primary solver pre-time step
  if(solverId == fvSolver(0).solverId()) {
    // correct every time step only after it has been updated in the primary solver
    if(m_correctCoarseBndry && fvSolver(0).m_levelSetMb
       && (fvSolver(0).m_reComputedBndry || globalTimeStep == fvSolver(0).m_restartTimeStep)) {
      for(MInt solver = 1; solver < noSolvers(); solver++) {
        correctCoarseBndryCells(solver);
      }
    }
  }

  if(fvSolver(0).m_reComputedBndry || globalTimeStep == fvSolver(0).m_restartTimeStep) {
    initSplitMapping();
  }
}

template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::postCouple(MInt) {
  TRACE();
  startTimer("Main");

  // Determine currently active solver/level via the solvers status
  MInt activeLevel = -1;
  for(MInt solverId = 0; solverId < noSolvers(); solverId++) {
    if(fvSolver(solverId).getSolverStatus()) {
      activeLevel = solverId;
      break;
    }
  }
  if(activeLevel == -1) {
    TERMM(1, "Multilevel postCouple(): no active solver/level found!");
  }

  // Restrict solution unless if this is the coarsest level (sawtooth cycle)
  if(activeLevel != noSolvers() - 1) {
    restriction(activeLevel);
  } else {
    for(MInt i = noSolvers() - 1; i > 0; i--) {
      prolongation(i);
    }
  }

  stopTimer("Main");
}


/// \brief Restrict the solution on the given fine level onto the next coarser level
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::restriction(const MInt level) {
  TRACE();

  startTimer("Restriction");

  // Sanity check for argument
  if(level < 0 || level > noSolvers() - 2) {
    TERMM(1,
          "Bad restriction level argument. Level designates the level which will be restricted, "
          "starting from 0 (finest level) to noSolvers() - 2 (next-to-coarsest level). Current value: "
              + to_string(level));
  }

  const MInt fineLevel = level;
  const MInt coarseLevel = level + 1;
  auto& fine = fvSolver(fineLevel);
  auto& coarse = fvSolver(coarseLevel);

  // re-Compute the right-hand side at fine level for restriction
  startTimer("ComputeFineRHS");
  fine.rhs();
  fine.rhsBnd();
  stopTimer("ComputeFineRHS");

  // Restrict data from fine to coarse level and store restricted RHS
  restrictData(fineLevel);

  coarse.exchangeData(&coarse.a_variable(0, 0), coarse.CV->noVariables);
  coarse.exchangeData(&coarse.a_oldVariable(0, 0), coarse.CV->noVariables);

  // Run lhsBnd to apply small cell correction and exchange data after restriction
  startTimer("ComputeCoarseLhsBnd");
  coarse.lhsBnd();

  stopTimer("ComputeCoarseLhsBnd");

  if(m_prolongationMethod == 2) {
    // Calculate and store the gradients on the coarse grid based on the new restricted data
    startTimer("CalcGradientsRestriction");
    coarse.LSReconstructCellCenterCons(0);
    stopTimer("CalcGradientsRestriction");
  }

  // Re-calculate the residual of the coarse solver
  coarse.rhs();
  coarse.rhsBnd();

  // Calculate coarse grid correction factor
  // (i.e. recompute the residual and compute the difference between the new residual
  // and the residual
  computeCoarseLevelCorrection(coarseLevel);

  // Store restricted variables
  // (i.e. store variables in the coarse grid after the restriction before the coarse time step!)
  storeRestrictedVariables(coarseLevel);

  stopTimer("Restriction");
}


/// \brief Reset the coarse level correction tau
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::resetTau(const MInt level) {
  TRACE();

  auto& coarse = fvSolver(level);
  const MInt noVars = coarse.CV->noVariables;

  // Reset coarse level correction "tau"
  const MInt noInternalCells = coarse.grid().noInternalCells();
  for(MInt cellId = 0; cellId < noInternalCells; cellId++) {
    for(MInt varId = 0; varId < noVars; varId++) {
      coarse.a_tau(cellId, varId) = 0.0;
    }
  }
}


/// \brief Restrict the data (fine-coarse) on the given level
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::restrictData(const MInt level) {
  TRACE();

  startTimer("RestrictData");

  // Store references and variables for convenience
  auto& fine = fvSolver(level);
  auto& coarse = fvSolver(level + 1);
  const MInt noVars = fine.CV->noVariables;

  ScratchSpace<MFloat> sumVolumes(fine.noInternalCells(), AT_, "sumVolumes");
  std::fill_n(&sumVolumes[0], fine.noInternalCells(), std::numeric_limits<MFloat>::quiet_NaN());

  // Initialize cell volumes on fine grid with the coarse grid cell volumes, if a matching fine grid
  // cell has children this value will be overwritten below by the sum of volumes of the child cells
  for(MInt cellId = 0; cellId < coarse.noInternalCells(); cellId++) {
    const MInt fineCellId = m_coarse2fine[level + 1][cellId];
    sumVolumes[fineCellId] = coarse.a_cellVolume(cellId);
  }

  // First, loop over all cells of the fine solver that need to be restricted to and restrict their
  // child cells' data
  vector<MInt> childLessCells;
  const MInt noRestrictedCells = m_restrictedCells[level].size();
  for(MInt restrictedCellId = 0; restrictedCellId < noRestrictedCells; restrictedCellId++) {
    const MInt cellId = m_restrictedCells[level][restrictedCellId];

    // Reset all relevant variables
    for(MInt varId = 0; varId < noVars; varId++) {
      fine.a_variable(cellId, varId) = 0.0;
      fine.a_oldVariable(cellId, varId) = 0.0;
      if(fine.m_dualTimeStepping) {
        fine.a_dt1Variable(cellId, varId) = 0.0;
        fine.a_dt2Variable(cellId, varId) = 0.0;
      }
      fine.a_rightHandSide(cellId, varId) = 0.0;
      fine.a_tau(cellId, varId) = 0.0;
      if(fine.m_levelSetMb) {
        fine.m_rhs0[cellId][varId] = 0.0;
      }
    }

    sumVolumes[cellId] = 0;

    // Loop over child cells
    MFloat volume = 0.0;
    for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
      // Skip if child does not exist
      const MInt child = fine.c_childId(cellId, childId);
      if(child == -1) {
        continue;
      }
      if(fine.a_isInactive(child)) {
        continue;
      }

      // Update total volume of restricted cells
      const MFloat childVolume = fine.a_cellVolume(child);
      volume += childVolume;
      sumVolumes[cellId] += childVolume;

      // Restrict data by volume average of the child data
      for(MInt varId = 0; varId < noVars; varId++) {
        fine.a_variable(cellId, varId) += fine.a_variable(child, varId) * childVolume;
        fine.a_oldVariable(cellId, varId) += fine.a_oldVariable(child, varId) * childVolume;
        if(fine.m_dualTimeStepping) {
          fine.a_dt1Variable(cellId, varId) += fine.a_dt1Variable(child, varId) * childVolume;
          fine.a_dt2Variable(cellId, varId) += fine.a_dt2Variable(child, varId) * childVolume;
        }
        fine.a_rightHandSide(cellId, varId) += fine.a_rightHandSide(child, varId);
        fine.a_tau(cellId, varId) += fine.a_tau(child, varId);
        if(fine.m_levelSetMb) {
          fine.m_rhs0[cellId][varId] += fine.m_rhs0[child][varId];
        }
      }
    }

    // Volume-average variables
    if(volume > 0) {
      for(MInt varId = 0; varId < noVars; varId++) {
        fine.a_variable(cellId, varId) = fine.a_variable(cellId, varId) / volume;
        fine.a_oldVariable(cellId, varId) = fine.a_oldVariable(cellId, varId) / volume;
        if(fine.m_dualTimeStepping) {
          fine.a_dt1Variable(cellId, varId) = fine.a_dt1Variable(cellId, varId) / volume;
          fine.a_dt2Variable(cellId, varId) = fine.a_dt2Variable(cellId, varId) / volume;
        }
      }
    } else {
      ASSERT(volume > -MFloatEps, "");
      const MInt coarseCellId = m_fine2coarse[level][restrictedCellId];
      if(coarseCellId > -1 && !coarse.a_isInactive(coarseCellId)) {
        childLessCells.push_back(restrictedCellId);
      }
    }
  }

  // parent/coarse cells which do not have
  ASSERT(childLessCells.empty(), "Warning, not matching children have been found for a multi-level solver");

  // Next, copy all data from fine to coarse solver and store restricted RHS
  for(MInt cellId = 0; cellId < coarse.noInternalCells(); cellId++) {
    const MInt fineCellId = m_coarse2fine[level + 1][cellId];
    for(MInt varId = 0; varId < noVars; varId++) {
      coarse.a_variable(cellId, varId) = fine.a_variable(fineCellId, varId);
      coarse.a_oldVariable(cellId, varId) = fine.a_oldVariable(fineCellId, varId);
      if(fine.m_dualTimeStepping) {
        coarse.a_dt1Variable(cellId, varId) = fine.a_dt1Variable(fineCellId, varId);
        coarse.a_dt2Variable(cellId, varId) = fine.a_dt2Variable(fineCellId, varId);
      }

      // Volume weighting factor for restricting the RHS and tau
      const MFloat volumeFactor = coarse.a_cellVolume(cellId) / sumVolumes[fineCellId];
      coarse.a_rightHandSide(cellId, varId) =
          volumeFactor * (fine.a_rightHandSide(fineCellId, varId) - fine.a_tau(fineCellId, varId));
      coarse.a_restrictedRHS(cellId, varId) = coarse.a_rightHandSide(cellId, varId);
      coarse.a_tau(cellId, varId) = volumeFactor * fine.a_tau(fineCellId, varId);
      if(coarse.m_levelSetMb) {
        coarse.m_rhs0[cellId][varId] = volumeFactor * (fine.m_rhs0[fineCellId][varId] - fine.a_tau(fineCellId, varId));
      }
    }
  }

  // update coarse split children based on closes finer child
  for(MInt sc = 0; sc < coarse.a_noSplitCells(); sc++) {
    MInt totalSplitId = 0;
    for(MInt scc = 0; scc < coarse.a_noSplitChilds(sc); scc++) {
      const MInt splitChild = coarse.a_splitChildId(sc, scc);
      const MInt finceChild = m_splitCellMapping[level + 1][totalSplitId++];
      for(MInt varId = 0; varId < noVars; varId++) {
        coarse.a_variable(splitChild, varId) = fine.a_variable(finceChild, varId);
        coarse.a_oldVariable(splitChild, varId) = fine.a_oldVariable(finceChild, varId);
        if(fine.m_dualTimeStepping) {
          coarse.a_dt1Variable(splitChild, varId) = fine.a_dt1Variable(finceChild, varId);
          coarse.a_dt2Variable(splitChild, varId) = fine.a_dt2Variable(finceChild, varId);
        }
        // Volume weighting factor for restricting the RHS and tau
        const MFloat volumeFactor = coarse.a_cellVolume(splitChild) / sumVolumes[finceChild];
        coarse.a_rightHandSide(splitChild, varId) =
            volumeFactor * (fine.a_rightHandSide(finceChild, varId) - fine.a_tau(finceChild, varId));
        coarse.a_restrictedRHS(splitChild, varId) = coarse.a_rightHandSide(splitChild, varId);
        coarse.a_tau(splitChild, varId) = volumeFactor * fine.a_tau(finceChild, varId);
        if(coarse.m_levelSetMb) {
          coarse.m_rhs0[splitChild][varId] =
              volumeFactor * (fine.m_rhs0[finceChild][varId] - fine.a_tau(finceChild, varId));
        }
      }
    }
  }


  stopTimer("RestrictData");
}


/// \brief Compute the coarse level correction tau
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::computeCoarseLevelCorrection(const MInt level) {
  TRACE();

  startTimer("ComputeCoarseCorrection");

  // Store references and variables for convenience
  auto& coarse = fvSolver(level);

  // Calculate coarse level correction "tau"
  // only required for internal cells, as the runge-Kutta step is only performed for those
  for(MInt cellId = 0; cellId < coarse.noInternalCells(); cellId++) {
    for(MInt varId = 0; varId < coarse.CV->noVariables; varId++) {
      coarse.a_tau(cellId, varId) = coarse.a_rightHandSide(cellId, varId) - coarse.a_restrictedRHS(cellId, varId);
    }
  }

  stopTimer("ComputeCoarseCorrection");
}


/// \brief Store the restricted variables (required for the prolongation)
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::storeRestrictedVariables(const MInt level) {
  TRACE();

  startTimer("StoreRestrictedVars");

  // Store references and variables for convenience
  auto& coarse = fvSolver(level);

  // compute conservative variables on halo cells
  // before conservative variables have been transferred from the fine grid for all internal cells,
  // primitive variables been calculated and primitive variables been exchanged
  //! There is an exception in the Enthalpy case where the conservative variables are exchanged!

  for(MInt cellId = coarse.noInternalCells(); cellId < coarse.a_noCells(); cellId++) {
    coarse.setConservativeVariables(cellId);
  }

  // Store restricted variables for all cells, as halo information might bee needed for interpolation
  for(MInt cellId = 0; cellId < coarse.a_noCells(); cellId++) {
    for(MInt varId = 0; varId < coarse.CV->noVariables; varId++) {
      ASSERT(!std::isnan(coarse.a_variable(cellId, varId)),
             "Restricted variables of cellId " + std::to_string(cellId) + " of domainId "
                 + std::to_string(coarse.domainId()) + "are nan");
      coarse.a_restrictedVar(cellId, varId) = coarse.a_variable(cellId, varId);
    }
  }

  stopTimer("StoreRestrictedVars");
}


/// \brief Prolong the solution on the given coarse level onto the next finer level
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::prolongation(const MInt level) {
  TRACE();

  startTimer("Prolongation");

  // Store references and variables for convenience
  auto& coarse = static_cast<SolverType&>(fvSolver(level));
  auto& fine = static_cast<SolverType&>(fvSolver(level - 1));

  // only halo cells only primary variables are  correct at this point!
  if(m_prolongationMethod == 2) {
    // Temporarily store gradients
    startTimer("StoreGradients");
    const MInt noCells = coarse.noInternalCells();
    const MInt noVars = coarse.CV->noVariables;
    MFloatScratchSpace slopes(noCells * noVars * nDim, AT_, "slopes");
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noVars; varId++) {
        for(MInt i = 0; i < nDim; i++) {
          slopes[cellId * noVars * nDim + varId * nDim + i] = coarse.a_storedSlope(cellId, varId, i);
        }
      }
    }
    stopTimer("StoreGradients");

    // Calculate new gradients after the completed coarse time step
    startTimer("CalcGradientsProlongation");
    coarse.LSReconstructCellCenterCons(0);
    stopTimer("CalcGradientsProlongation");

    // Compute delta in gradients
    startTimer("CalcGradientDelta");
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noVars; varId++) {
        for(MInt i = 0; i < nDim; i++) {
          coarse.a_storedSlope(cellId, varId, i) -= slopes[cellId * noVars * nDim + varId * nDim + i];
        }
      }
    }
    stopTimer("CalcGradientDelta");
  }

  prolongData(level);

  // Re-compute boundary conditions etc. based on new/corrected variable state in the fine cell
  startTimer("FineLhsBnd");
  fine.exchangeData(&fine.a_variable(0, 0), fine.CV->noVariables);

  fine.lhsBnd();

  stopTimer("FineLhsBnd");

  stopTimer("Prolongation");
}


/// \brief Prolong the data of a coarse level
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::prolongData(const MInt level) {
  TRACE();

  startTimer("ProlongData");

  // Store references and variables for convenience
  const MInt coarseLevel = level;
  const MInt fineLevel = level - 1;
  auto& coarse = fvSolver(coarseLevel);
  auto& fine = fvSolver(fineLevel);
  const MInt noVars = coarse.CV->noVariables;
  const MInt noInternalCells = coarse.noInternalCells();

  // compute difference to the restricted variables on the coarse grid
  // necessary on halo and ghost-cells as their change might be required for interpolation!
  for(MInt cellId = 0; cellId < coarse.a_noCells(); cellId++) {
    if(cellId < coarse.c_noCells() && !coarse.c_isLeafCell(cellId)) {
      continue;
    }
    for(MInt varId = 0; varId < noVars; varId++) {
      coarse.a_variable(cellId, varId) -= coarse.a_restrictedVar(cellId, varId);
    }
  }


  // Reset restricted cell variables in fine solver
  const MInt noRestrictedCells = m_restrictedCells[fineLevel].size();
  for(MInt restrictedCellId = 0; restrictedCellId < noRestrictedCells; restrictedCellId++) {
    const MInt cellId = m_restrictedCells[fineLevel][restrictedCellId];
    for(MInt varId = 0; varId < noVars; varId++) {
      fine.a_variable(cellId, varId) = 0.0;
    }
  }

  // Add coarse-level correction to cells on fine level and copy slopes
  for(MInt cellId = 0; cellId < noInternalCells; cellId++) {
    const MInt fineCellId = m_coarse2fine[coarseLevel][cellId];
    switch(m_prolongationMethod) {
      case 0:
      case 1: {
        for(MInt varId = 0; varId < noVars; varId++) {
          fine.a_variable(fineCellId, varId) += coarse.a_variable(cellId, varId);
        }
        break;
      }
      case 2: {
        for(MInt varId = 0; varId < noVars; varId++) {
          fine.a_variable(fineCellId, varId) += coarse.a_variable(cellId, varId);
          for(MInt i = 0; i < nDim; i++) {
            fine.a_storedSlope(fineCellId, varId, i) = coarse.a_storedSlope(cellId, varId, i);
          }
        }
        break;
      }
      default: {
        mTerm(1, AT_, "Unknown prolongation method");
      }
    }
  }

  std::function<MFloat(const MInt, const MInt)> consVar = [&](const MInt cellId, const MInt varId) {
    return static_cast<MFloat>(coarse.a_variable(cellId, varId));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(coarse.a_coordinate(cellId, id));
  };


  // interpolate correction from restricted cells to their children
  for(MInt restrictedId = 0; restrictedId < noRestrictedCells; restrictedId++) {
    const MInt cellId = m_restrictedCells[fineLevel][restrictedId];
    const MInt coarseCellId = m_fine2coarse[fineLevel][restrictedId];

    // Loop over child cells to prolong correction
    for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
      // Skip if child does not exist
      const MInt child = fine.c_childId(cellId, childId);
      if(child == -1) {
        continue;
      }
      if(fine.a_isInactive(child)) {
        continue;
      }

      switch(m_prolongationMethod) {
        case 0: {
          for(MInt varId = 0; varId < noVars; varId++) {
            fine.a_variable(child, varId) += fine.a_variable(cellId, varId);
          }
          break;
        }
        case 1: {
          MInt interpolationCells[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
          MFloat point[nDim];
          for(MInt i = 0; i < nDim; i++) {
            point[i] = fine.a_coordinate(child, i);
          }
          if(!fine.a_isBndryCell(child)) {
            // 1) find regular Cartesian stencil

            // NOTE: prevent access to inactive neighbors here
            MBool backup = coarse.m_deleteNeighbour;
            coarse.m_deleteNeighbour = false;
            std::function<MBool(const MInt, const MInt)> alwaysTrue = [&](const MInt, const MInt) { return true; };

            const MInt position =
                coarse.setUpInterpolationStencil(coarseCellId, interpolationCells, point, alwaysTrue, false);
            coarse.m_deleteNeighbour = backup;

            ASSERT(position > -1, "Return value of setUpInterpolationStencil is > -1");

            // 2) stencil check:
            //   - replace inactive cell with bndry-ghost cell of the coarse cell
            //   - replace with bndry ghost cell if a small-cell
            //  NOTE: possible alternative would be the ghost cell of the largest cell the stencil!
            for(MInt id = 0; id < IPOW2(nDim); id++) {
              const MInt interpolationCell = interpolationCells[id];
              if(interpolationCell < 0) {
                mTerm(1, AT_, "Incorrect stencil!");
              }

              if(coarse.a_isInactive(interpolationCell)) {
                ASSERT(coarse.a_isBndryCell(coarseCellId), "");
                const MInt bndryId = coarse.a_bndryId(coarseCellId);
                const MInt ghostCellId = coarse.a_bndryGhostCellId(bndryId, 0);
                if(ghostCellId > -1) {
                  // meaning not a split-cell
                  interpolationCells[id] = ghostCellId;
                }
              } else if(coarse.a_cellVolume(interpolationCell)
                            / coarse.c_cellVolumeAtLevel(coarse.a_level(interpolationCell))
                        < coarse.m_fvBndryCnd->m_volumeLimitWall) {
                const MInt bndryId = coarse.a_bndryId(interpolationCell);
                const MInt ghostCellId = coarse.a_bndryGhostCellId(bndryId, 0);
                if(ghostCellId > -1) {
                  interpolationCells[id] = ghostCellId;
                }
              }
            }

            // use least-squares interpolation
            // NOTE: other alternatives could be a cartesian-interpolation
            for(MInt varId = 0; varId < noVars; varId++) {
              fine.a_variable(child, varId) +=
                  coarse.leastSquaresInterpolation(&interpolationCells[0], &point[0], varId, consVar, coordinate);
            }
          } else {
            // change stencil when the cell is a bndry-Cell
            const MInt position = coarse.setUpBndryInterpolationStencil(coarseCellId, interpolationCells, point);

            for(MInt id = 0; id < IPOW2(nDim); id++) {
              const MInt intCellId = interpolationCells[id];
              if(intCellId < 0) {
                cerr << position << " " << coarseCellId << " " << coarse.domainId() << endl;
                mTerm(1, AT_, "Incorrect stencil at bndry!");
              } else {
                ASSERT(!coarse.a_isInactive(intCellId), "");
              }
            }

            ASSERT(position == IPOW2(nDim), "");
            for(MInt varId = 0; varId < noVars; varId++) {
              fine.a_variable(child, varId) +=
                  coarse.leastSquaresInterpolation(&interpolationCells[0], &point[0], varId, consVar, coordinate);
            }
          }
          break;
        }
        case 2: {
          // Compute distance between cell centers based on the volumetric coordinates
          array<MFloat, nDim> dx{};
          for(MInt i = 0; i < nDim; i++) {
            dx.at(i) = fine.a_coordinate(child, i) - fine.a_coordinate(cellId, i);
          }
          for(MInt varId = 0; varId < noVars; varId++) {
            fine.a_variable(child, varId) += fine.a_variable(cellId, varId);
            for(MInt i = 0; i < nDim; i++) {
              fine.a_variable(child, varId) += fine.a_storedSlope(cellId, varId, i) * dx.at(i);
            }
          }
          break;
        }
        default: {
          mTerm(1, AT_, "Unknown prolongation method");
        }
      }
    }
  }

  stopTimer("ProlongData");
}


/// \brief for each solver that will be restricted (all except the coarsest), store restricted cells
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::initRestrictedCells() {
  TRACE();

  m_restrictedCells.resize(noSolvers());
  for(MInt solverId = 0; solverId < noSolvers() - 1; solverId++) {
    auto& fineGrid = fvSolver(solverId).grid();
    const MInt noInternalCells = fineGrid.noInternalCells();
    const MInt restrictedLevel = fineGrid.maxLevel() - 1;

    // First, count restricted cells and store them temporarily
    MInt noRestrictedCells = 0;
    MIntScratchSpace restrictedCells(noInternalCells, AT_, "restrictedCells");
    for(MInt i = 0; i < noInternalCells; i++) {
      if(fineGrid.tree().level(i) == restrictedLevel && fineGrid.tree().hasChildren(i)) {
        restrictedCells[noRestrictedCells] = i;
        noRestrictedCells++;
      }
    }

    // Second, move recorded cells to permanent storage
    m_restrictedCells[solverId].resize(noRestrictedCells);
    copy_n(restrictedCells.data(), noRestrictedCells, m_restrictedCells[solverId].data());
  }
}


/// \brief for each coarse solver, store cell map to next finer solver
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::initMapping() {
  TRACE();

  m_coarse2fine.resize(noSolvers());
  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    auto& fineGrid = fvSolver(solverId - 1).grid();
    auto& coarseGrid = fvSolver(solverId).grid();
    const MInt noInternalCells = coarseGrid.noInternalCells();
    m_coarse2fine[solverId].resize(noInternalCells);
    for(MInt i = 0; i < noInternalCells; i++) {
      m_coarse2fine[solverId][i] = fineGrid.tree().grid2solver(coarseGrid.tree().solver2grid(i));
    }
  }

  // based on restricted cells/restricted Id, so only for the number of restricted cells!
  m_fine2coarse.resize(noSolvers());
  for(MInt solverId = 0; solverId < noSolvers() - 1; solverId++) {
    auto& fineGrid = fvSolver(solverId).grid();
    auto& coarseGrid = fvSolver(solverId + 1).grid();
    const MInt noRestrictedCells = m_restrictedCells[solverId].size();
    m_fine2coarse[solverId].resize(noRestrictedCells);
    for(MInt i = 0; i < noRestrictedCells; i++) {
      const MInt fineCellId = m_restrictedCells[solverId][i];
      m_fine2coarse[solverId][i] = coarseGrid.tree().grid2solver(fineGrid.tree().solver2grid(fineCellId));
    }
  }
}

/// \brief for each coarse solver, store cell map of split-cells to
///  next finer solver childs
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::initSplitMapping() {
  TRACE();

  m_splitCellMapping.clear();
  m_splitCellMapping.resize(noSolvers());
  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    auto& fine = fvSolver(solverId - 1);
    auto& coarse = fvSolver(solverId);
    if(coarse.m_totalnosplitchilds == 0) {
      continue;
    }
    m_splitCellMapping[solverId].resize(coarse.m_totalnosplitchilds);
    MInt totalSplitId = 0;
    for(MInt sc = 0; sc < coarse.a_noSplitCells(); sc++) {
      const MInt splitCell = coarse.a_splitCellId(sc);
      const MInt fineCell = m_coarse2fine[solverId][splitCell];
      for(MInt scc = 0; scc < coarse.a_noSplitChilds(sc); scc++) {
        const MInt splitChild = coarse.a_splitChildId(sc, scc);
        std::map<MFloat, MInt> distances;
        // loop over all fine-cell childs and find child with lowest distance to the split-child
        for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
          const MInt child = fine.c_childId(fineCell, childId);
          if(child == -1) {
            continue;
          }
          MFloat dist = 0;
          if(fine.a_isSplitCell(fineCell)) {
            for(MInt fsc = 0; fsc < fine.a_noSplitCells(); fsc++) {
              const MInt fineSplitCell = fine.a_splitCellId(fsc);
              if(fineSplitCell == fineCell) {
                for(MInt fscc = 0; fscc < fine.a_noSplitChilds(fsc); fscc++) {
                  dist = 0;
                  const MInt fineSplitChild = fine.a_splitChildId(fsc, fscc);
                  for(MInt i = 0; i < nDim; i++) {
                    dist += POW2(coarse.a_coordinate(splitChild, i) - fine.a_coordinate(fineSplitChild, i));
                  }
                  dist = sqrt(dist);
                  distances.insert(make_pair(dist, splitChild));
                }
              }
            }
            continue;
          }
          if(fine.a_isInactive(child)) {
            continue;
          }
          for(MInt i = 0; i < nDim; i++) {
            dist += POW2(coarse.a_coordinate(splitChild, i) - fine.a_coordinate(child, i));
          }
          dist = sqrt(dist);
          distances.insert(make_pair(dist, child));
        }

        if(distances.empty()) {
          cerr << "Neighbor search necessary!" << distances.size() << " " << distances.begin()->second << " "
               << coarse.a_isHalo(splitCell) << endl;
        }

        const MInt childId = distances.begin()->second;
        m_splitCellMapping[solverId][totalSplitId++] = childId;
      }
    }
  }
}
/// \brief store internal leaf cells for each solver except for finest solver, for which they are not needed
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::initLeafCells() {
  TRACE();

  // Count internal leaf cells in each solver except for finest solver, for which they are not needed
  std::vector<MInt> noLeafCells(noSolvers(), 0);
  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    const MInt noInternalCells = fvSolver(solverId).grid().noInternalCells();
    for(MInt i = 0; i < noInternalCells; i++) {
      if(fvSolver(solverId).grid().tree().isLeafCell(i)) {
        noLeafCells[solverId]++;
      }
    }
  }

  m_leafCells.resize(noSolvers());
  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    m_leafCells[solverId].resize(noLeafCells[solverId]);
    const MInt noInternalCells = fvSolver(solverId).grid().noInternalCells();
    MInt index = 0;
    for(MInt i = 0; i < noInternalCells; i++) {
      if(fvSolver(solverId).grid().tree().isLeafCell(i)) {
        m_leafCells[solverId][index] = i;
        index++;
      }
    }
  }
}

/// \brief for each solver, reset multilevel variables
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::resetMultiLevelVariables() {
  TRACE();

  const auto nan = std::numeric_limits<MFloat>::quiet_NaN();
  for(MInt solverId = 0; solverId < noSolvers(); solverId++) {
    auto& current = fvSolver(solverId);
    const MInt noVars = current.CV->noVariables;
    const MInt noCells = current.a_noCells();
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noVars; varId++) {
        current.a_tau(cellId, varId) = 0.0;
        current.a_restrictedRHS(cellId, varId) = nan;
        current.a_restrictedVar(cellId, varId) = nan;
        for(MInt i = 0; i < nDim; i++) {
          current.a_storedSlope(cellId, varId, i) = nan;
        }
      }
    }
  }
}

/// \brief for each solver except the finest, reset variables to NaN to ensure that only restricted
/// information is used
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::resetSecondaryFlowVariables() {
  TRACE();

  const auto nan = std::numeric_limits<MFloat>::quiet_NaN();

  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    auto& current = fvSolver(solverId);
    const MInt noVars = current.CV->noVariables;
    const MInt noCells = current.noInternalCells();
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noVars; varId++) {
        current.a_variable(cellId, varId) = nan;
        current.a_rightHandSide(cellId, varId) = nan;
        if(current.m_levelSetMb) {
          current.m_rhs0[cellId][varId] = nan;
        }
      }

      for(MInt varId = 0; varId < noVars; varId++) {
        current.a_pvariable(cellId, varId) = nan;
      }
    }
  }
}

/**
 * corrects the volume and the center of gravity of the boundary cells
 * using fine cell information
 * creates body surfaces using fine cell information
 * does not work correctly with multipleGhostCell formulation (claudia)!
 * See D.Hartman et al Computers & Fluids 37 (2008) 1103-1125: concept of averaged control volumes
 * it is required that small and master cells are not yet merged!!!
 *
 * @author Daniel Hartmann, March 10, 2006, Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::correctCoarseBndryCells(const MInt solverId) {
  TRACE();

  auto& fine = fvSolver(solverId - 1);
  auto& coarse = fvSolver(solverId);

  if(coarse.domainId() == 0 && globalTimeStep < 1) {
    cerr << "Correxting secondary bndry-cell volumes and centers " << endl;
  }

  const MInt noCoarseCells = coarse.a_noCells();
  ScratchSpace<MFloat> bndryCoordinates(noCoarseCells, nDim, AT_, "bndryCoordinates");
  ScratchSpace<MFloat> cellCoordinates(noCoarseCells, nDim, AT_, "cellCoordinates");
  ScratchSpace<MFloat> srfcCoordinates(noCoarseCells, nDim, AT_, "srfcCoordinates");
  ScratchSpace<MFloat> volumes(noCoarseCells, AT_, "volumes");

  const auto nan = std::numeric_limits<MFloat>::quiet_NaN();
  fill(bndryCoordinates.begin(), bndryCoordinates.end(), nan);
  fill(cellCoordinates.begin(), cellCoordinates.end(), nan);
  fill(srfcCoordinates.begin(), srfcCoordinates.end(), nan);
  fill(volumes.begin(), volumes.end(), nan);

  // loop over all boundary cells on coarse solver
  for(MInt bndryId = 0; bndryId < coarse.m_bndryCells->size(); bndryId++) {
    // choose only those boundary cells at the current level
    const MInt coarseCellId = coarse.m_bndryCells->a[bndryId].m_cellId;
    ASSERT(coarseCellId > -1, "");
    if(!coarse.c_isLeafCell(coarseCellId)) {
      continue;
    }

    if(coarse.a_isHalo(coarseCellId)) {
      continue;
    }
    // NOTE: data is exchanged for halo-cells below as the different solvers might not have matching
    //      lower level halo-layers
    const MInt fineCellId = m_coarse2fine[solverId][coarseCellId];
    if(fineCellId < 0) {
      continue;
    }

    if(fine.grid().tree().hasChildren(fineCellId)) {
      const MInt noChildren = fine.grid().tree().noChildren(fineCellId);
      if(noChildren == 0) {
        continue;
      }

      // reset temp. variables
      MFloat volume = 0;
      MFloat bndryCellVolumes = 0;
      MFloat srfcArea = 0;

      array<MFloat, nDim> bndryXyz{};
      array<MFloat, nDim> cellXyz{};
      array<MFloat, nDim> srfcXyz{};

      bndryXyz.fill(0);
      cellXyz.fill(0);
      srfcXyz.fill(0);

      // -------- average children cell coordinates and volumes ---------
      for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
        const MInt childCellId = fine.c_childId(fineCellId, childId);
        if(childCellId == -1) {
          continue;
        }
        if(fine.a_isInactive(childCellId)) {
          continue;
        }

        ASSERT(fine.c_isLeafCell(childCellId), "");

        const MInt bndryCellId = fine.a_bndryId(childCellId);
        if(bndryCellId > -1) {
          volume += fine.m_bndryCells->a[bndryCellId].m_volume;
          bndryCellVolumes += fine.m_bndryCells->a[bndryCellId].m_volume;

          for(MInt i = 0; i < nDim; i++) {
            bndryXyz.at(i) += fine.a_coordinate(childCellId, i) * fine.m_bndryCells->a[bndryCellId].m_volume;
          }

          // average body surface of children
          srfcArea += fine.m_bndryCells->a[bndryCellId].m_srfcs[0]->m_area;

          for(MInt i = 0; i < nDim; i++) {
            srfcXyz.at(i) += fine.m_bndryCells->a[bndryCellId].m_srfcs[0]->m_coordinates[i]
                             * fine.m_bndryCells->a[bndryCellId].m_srfcs[0]->m_area;
          }

        } else {
          const MFloat childVolume = fine.c_cellVolumeAtLevel(fine.a_level(childCellId));
          volume += childVolume;
          for(MInt i = 0; i < nDim; i++) {
            cellXyz.at(i) += fine.a_coordinate(childCellId, i) * childVolume;
          }
        }
      }

      for(MInt i = 0; i < nDim; i++) {
        cellXyz[i] += bndryXyz[i];
        cellXyz[i] = cellXyz[i] / volume;
        bndryXyz[i] = bndryXyz[i] / bndryCellVolumes;
        srfcXyz[i] = srfcXyz[i] / srfcArea;
      }
      // ----------------------------------------------------------------

      // correct the coordinate shift of the boundary cell
      // set the coordinates of the boundary cell
      for(MInt i = 0; i < nDim; i++) {
        coarse.m_bndryCells->a[bndryId].m_coordinates[i] += cellXyz[i] - coarse.a_coordinate(coarseCellId, i);
        coarse.a_coordinate(coarseCellId, i) = cellXyz[i];
        coarse.m_bndryCells->a[bndryId].m_srfcs[0]->m_coordinates[i] = srfcXyz[i];


        bndryCoordinates(coarseCellId, i) = coarse.m_bndryCells->a[bndryId].m_coordinates[i];
        cellCoordinates(coarseCellId, i) = cellXyz[i];
        srfcCoordinates(coarseCellId, i) = srfcXyz[i];
      }

      // update coarse cell volume
      coarse.m_bndryCells->a[bndryId].m_volume = volume;

      volumes(coarseCellId) = volume;
    }
  }

  // exchange data
  coarse.exchangeData(&volumes(0, 0), 1);
  coarse.exchangeData(&bndryCoordinates(0, 0), nDim);
  coarse.exchangeData(&cellCoordinates(0, 0), nDim);
  coarse.exchangeData(&srfcCoordinates(0, 0), nDim);

  // set data for halo cells
  for(MInt bndryId = 0; bndryId < coarse.m_bndryCells->size(); bndryId++) {
    const MInt coarseCellId = coarse.m_bndryCells->a[bndryId].m_cellId;

    if(coarseCellId < 0) {
      continue;
    }
    if(!coarse.a_isHalo(coarseCellId)) {
      continue;
    }
    if(coarse.a_level(coarseCellId) != coarse.maxLevel()) {
      continue;
    }

    for(MInt i = 0; i < nDim; i++) {
      coarse.m_bndryCells->a[bndryId].m_coordinates[i] = bndryCoordinates(coarseCellId, i);
      coarse.a_coordinate(coarseCellId, i) = cellCoordinates(coarseCellId, i);
      coarse.m_bndryCells->a[bndryId].m_srfcs[0]->m_coordinates[i] = srfcCoordinates(coarseCellId, i);
    }

    ASSERT(volumes(coarseCellId) > -MFloatEps, "coarse multi-level solver has negative-cell volume after correction!");
    coarse.m_bndryCells->a[bndryId].m_volume = volumes(coarseCellId);
  }

  coarse.computeCellVolumes();
}


/// Perform different checks to ensure that everything is properly initialized
/// and configured.
/// mode 0: skip checks if the solver cells are created during the initial adaptation
/// mode 1: second call after the initial adaptation (=> don't skip checks!)
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-10-01
template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::sanityCheck(const MInt mode) {
  TRACE();

  // Check #1: Ensure that no small boundary cells are detected
  //           (they are not yet handled properly)
  for(MInt i = 0; i < noSolvers(); i++) {
    auto& b = fvSolver(i);
    if(b.m_fvBndryCnd->m_smallBndryCells->size() > 0) {
      TERMM(1, "Old concept of merged cells not (yet) supported");
    }
    if(!b.m_fvBndryCnd->m_smallCutCells.empty()) {
      cerr << "Multilevel: small cells might not be handled properly yet, continue anyway..." << endl;
    }
  }

  // Check #2: Ensure that we have at least two solvers
  if(noSolvers() < 2) {
    TERMM(1, "Multilevel computation does not make sense with less than two solvers");
  }

  // if the solver cells are created during an initial adaptation the checks below need
  // to be done at a later point!
  if(mode == 0 && fvSolver(0).m_adaptation && globalTimeStep == 0) {
    return;
  }

  // Check #3: Ensure that the maximum level is reduced by one with each solver and that the fv-solvers
  //           are of the same type
  const MInt overallMaxLevel = fvSolver(0).grid().maxLevel();
  for(MInt solverId = 1; solverId < noSolvers(); solverId++) {
    if(fvSolver(solverId).grid().maxLevel() != overallMaxLevel - solverId) {
      TERMM(1,
            "Maximum refinement level not reduced by one between solvers " + to_string(solverId - 1) + " and "
                + to_string(solverId));
    }
    if(fvSolver(solverId).m_levelSetMb != fvSolver(0).m_levelSetMb) {
      TERMM(1, "Multilevel coupling of different fv-solver types not allowed yet!");
    }
  }

  // Check #4: Make sure that for each solver that is restricted, the cells at the restricted level
  // match the leaf cells of the next coarser solver
  for(MInt solverId = 0; solverId < noSolvers() - 1; solverId++) {
    const MInt noInternalCells = fvSolver(solverId).grid().noInternalCells();
    const MInt restrictedLevel = fvSolver(solverId).grid().maxLevel() - 1;
    MInt count = 0;
    for(MInt i = 0; i < noInternalCells; i++) {
      // Check all cells that are either at the restricted level or that are coarser and leaf cells,
      // as they should match the leaf cells of the next coarser solver
      if(fvSolver(solverId).grid().tree().level(i) == restrictedLevel
         || (fvSolver(solverId).grid().tree().level(i) < restrictedLevel
             && fvSolver(solverId).grid().tree().isLeafCell(i))) {
        // The grid cell id of the restricted cells should match the grid cell id of the leaf cell
        if(fvSolver(solverId).grid().tree().solver2grid(i)
           == fvSolver(solverId + 1).grid().tree().solver2grid(m_leafCells[solverId + 1][count])) {
          count++;
        } else {
          TERMM(1,
                "Mismatch between restricted cells of solver " + to_string(solverId) + " and leaf cells of solver "
                    + to_string(solverId + 1));
        }
      }
    }

    // Also check if count matches the total number of leaf cells of the coarse solver
    if(count != static_cast<MInt>(m_leafCells[solverId + 1].size())) {
      TERMM(1, "Number of restricted cells and number of leaf cells does not match");
    }
  }


  if(mode == 1) {
    // Check #N: Ensure that "restricted" cells at level n are identical to leaf cells at level n + 1
    for(MInt solverId = 0; solverId < noSolvers() - 1; solverId++) {
      auto& fine = fvSolver(solverId);
      auto& coarse = fvSolver(solverId + 1);
      for(MInt restrictedId = 0; restrictedId < (signed)m_restrictedCells[solverId].size(); restrictedId++) {
        const MInt fineCellId = m_restrictedCells[solverId][restrictedId];
        const MInt coarseCellId = m_fine2coarse[solverId][restrictedId];
        ASSERT(fineCellId == m_coarse2fine[solverId + 1][coarseCellId],
               "fineCellId is eqal to coarseCellId on the n+1 solver");
        ASSERT(fine.grid().tree().solver2grid(fineCellId) == coarse.grid().tree().solver2grid(coarseCellId),
               "solver2grid of the finecellId is equal to solver2grid of the coarseCellId ");
        ASSERT(fine.a_level(fineCellId) == coarse.a_level(coarseCellId) && fine.c_noChildren(fineCellId) > 0,
               "fineCellId is on the same level as coarseCellId and children exists");
      }
    }
  }
}


/// Return existing timer id for given timer name or die.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-11-26
///
/// \param[in] name Name of the timer. If it does not exist, MAIA aborts.
template <MInt nDim, class SysEqn>
MInt CouplerFvMultilevel<nDim, SysEqn>::timer(const MString& name) {
  if(m_timers.count(name) == 0) {
    TERMM(1, "Timer '" + name + "' does not exist.");
  }
  return m_timers.at(name);
}


template <MInt nDim, class SysEqn>
MInt& CouplerFvMultilevel<nDim, SysEqn>::createTimer(const MString& name) {
  TRACE();

  // First ensure that there are no name conflicts
  if(m_timers.count(name) > 0) {
    TERMM(1, "Timer '" + name + "' already exists.");
  }

  // Create new timer and initialize with dummy id
  m_timers[name] = -1;

  // Return reference to timer id
  return m_timers[name];
}


template <MInt nDim, class SysEqn>
MBool CouplerFvMultilevel<nDim, SysEqn>::startTimer(const MString& name) {
  if(name.empty()) {
    TERMM(1, "Empty timer name");
  }
  RECORD_TIMER_START(timer(name));
  return true;
}


template <MInt nDim, class SysEqn>
MBool CouplerFvMultilevel<nDim, SysEqn>::stopTimer(const MString& name) {
  if(name.empty()) {
    TERMM(1, "Empty timer name");
  }
  RECORD_TIMER_STOP(timer(name));
  return true;
}


template <MInt nDim, class SysEqn>
void CouplerFvMultilevel<nDim, SysEqn>::initTimers() {
  TRACE();

  // Create timer group & timer for coupler, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "CouplerFvMultilevel");
  NEW_TIMER_NOCREATE(createTimer("Total"), "Total object lifetime", m_timerGroup);
  RECORD_TIMER_START(timer("Total"));

  // Timer for initialization
  NEW_SUB_TIMER_NOCREATE(createTimer("Init"), "Initialization", timer("Total"));

  // Create timer for performance-relevant run information
  NEW_SUB_TIMER_NOCREATE(createTimer("Main"), "Post couple", timer("Total"));
  NEW_SUB_TIMER_NOCREATE(createTimer("Restriction"), "Restriction", timer("Main"));
  NEW_SUB_TIMER_NOCREATE(createTimer("Prolongation"), "Prolongation", timer("Main"));
  NEW_SUB_TIMER_NOCREATE(createTimer("I/O"), "I/O", timer("Main"));

  // Restriction
  NEW_SUB_TIMER_NOCREATE(createTimer("ComputeFineRHS"), "fine.rhs(),fine.rhsBnd()", timer("Restriction"));
  NEW_SUB_TIMER_NOCREATE(createTimer("RestrictData"), "restrictData()", timer("Restriction"));
  NEW_SUB_TIMER_NOCREATE(createTimer("ComputeCoarseLhsBnd"), "coarse.lhsBnd()", timer("Restriction"));
  NEW_SUB_TIMER_NOCREATE(
      createTimer("CalcGradientsRestriction"), "coarse.LSReconstructCellCenterCons()", timer("Restriction"));
  NEW_SUB_TIMER_NOCREATE(
      createTimer("ComputeCoarseCorrection"), "computeCoarseLevelCorrection()", timer("Restriction"));
  NEW_SUB_TIMER_NOCREATE(createTimer("StoreRestrictedVars"), "storeRestrictedVariables()", timer("Restriction"));

  // Prolongation
  NEW_SUB_TIMER_NOCREATE(createTimer("StoreGradients"), "store gradients", timer("Prolongation"));
  NEW_SUB_TIMER_NOCREATE(
      createTimer("CalcGradientsProlongation"), "coarse.LSReconstructCellCenterCons()", timer("Prolongation"));
  NEW_SUB_TIMER_NOCREATE(createTimer("CalcGradientDelta"), "calc gradient delta", timer("Prolongation"));
  NEW_SUB_TIMER_NOCREATE(createTimer("ProlongData"), "prolongData()", timer("Prolongation"));
  NEW_SUB_TIMER_NOCREATE(createTimer("FineLhsBnd"), "fine.lhsBnd()", timer("Prolongation"));
}


/// \brief Constructor for multilevel interpolation coupler
template <MInt nDim, class SysEqn>
CouplerFvMultilevelInterpolation<nDim, SysEqn>::CouplerFvMultilevelInterpolation(
    const MInt couplingId, std::vector<FvCartesianSolverXD<nDim, SysEqn>*> solvers)
  : Base(couplingId), BaseCoupler(couplingId, solvers) {
  TRACE();
}


/// \brief Coupling routine to transfer coarse level restart data to the finest level/solver
template <MInt nDim, class SysEqn>
void CouplerFvMultilevelInterpolation<nDim, SysEqn>::postCouple(MInt) {
  TRACE();
  startTimer("Main");

  for(MInt solverId = 0; solverId < noSolvers(); solverId++) {
    if(fvSolver(solverId).getSolverStatus()) {
      TERMM(1, "Multilevel interpolation postCouple(): found active solver/level - all solvers should be inactive!");
    }
  }

  // Reset solution unless for the coarsest level with the actual restart data
  for(MInt i = 0; i < noSolvers() - 1; i++) {
    reset(i); // Resets also the restricted variables for the next coarser level
  }

  // Prolong the coarse level restart data to the finest level (since the restricted vars and the stored slopes are set
  // to 0 above the full solution is transferred to the next finer level)
  for(MInt i = noSolvers() - 1; i > 0; i--) {
    prolongation(i);
  }

  stopTimer("Main");
}


/// \brief Reset variables/slopes of a level/solver
template <MInt nDim, class SysEqn>
void CouplerFvMultilevelInterpolation<nDim, SysEqn>::reset(const MInt level) {
  TRACE();

  auto& fine = fvSolver(level);
  auto& coarse = fvSolver(level + 1);
  const MInt noVars = fine.CV->noVariables;

  for(MInt cellId = 0; cellId < fine.noInternalCells(); cellId++) {
    for(MInt varId = 0; varId < noVars; varId++) {
      // Reset variables on this level
      fine.a_variable(cellId, varId) = 0.0;
      for(MInt dim = 0; dim < nDim; dim++) {
        // ... and slopes (of conservative variables)
        fine.a_storedSlope(cellId, varId, dim) = 0.0;
      }
    }
  }

  for(MInt cellId = 0; cellId < coarse.noInternalCells(); cellId++) {
    for(MInt varId = 0; varId < noVars; varId++) {
      // Reset restricted variables on next coarser level
      coarse.a_restrictedVar(cellId, varId) = 0.0;
      for(MInt dim = 0; dim < nDim; dim++) {
        // ... and slopes (of conservative variables)
        coarse.a_storedSlope(cellId, varId, dim) = 0.0;
      }
    }
  }
}


// Explicit template instantiations
template class CouplerFvMultilevel<2, FvSysEqnNS<2>>;
template class CouplerFvMultilevel<3, FvSysEqnNS<3>>;
template class CouplerFvMultilevel<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class CouplerFvMultilevel<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplerFvMultilevel<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class CouplerFvMultilevel<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplerFvMultilevel<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplerFvMultilevel<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;

template class CouplerFvMultilevelInterpolation<2, FvSysEqnNS<2>>;
template class CouplerFvMultilevelInterpolation<3, FvSysEqnNS<3>>;
