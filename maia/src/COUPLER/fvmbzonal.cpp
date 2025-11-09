// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvmbzonal.h"

#include <algorithm>
#include <cmath>
#include <stack>
#include <vector>
#include "FV/fvcartesiansolverxd.h"
#include "MEMORY/alloc.h"
#include "MEMORY/scratch.h"
#include "globals.h"
#include "globalvariables.h"

using namespace std;


template <MInt nDim, class SysEqn>
CouplerFvMbZonal<nDim, SysEqn>::CouplerFvMbZonal(const MInt couplingId, FvMbSolver* upStream_, FvMbSolver* downStream_)
  : Coupling(couplingId), CouplingFvMb<nDim, SysEqn>(couplingId, upStream_) {
  TRACE();

  // set solver pointer and ids
  m_solverUp = upStream_;
  m_solverDown = downStream_;
  m_couplingId = couplingId;

  m_upStreamId = upStream().m_solverId;
  m_downStreamId = downStream().m_solverId;

  readProperties();
}

/** \brief reads lsfvmb-coupling-specific data
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::readProperties() {
  m_noRKSteps = 5;
  m_noRKSteps = Context::getSolverProperty<MInt>("noRKSteps", upStream().solverId(), AT_, &m_noRKSteps);

  ASSERT(m_noRKSteps == Context::getSolverProperty<MInt>("noRKSteps", downStream().solverId(), AT_, &m_noRKSteps), "");
  m_zonalMethod = "";
  /*! \page propertyPage1
    \section cutOffMethod
    <code>MString zonalMethod </code>\n
    default = <code>""</code>\n \n
    Specifies the zonal method to be used.\n
    Any method which is also specified as a cutOff-method is possible.
    Currently however only the following are implemented yet:
    <ul>
    <li> B - box cut-off  </li>
    <li> P - plane cut-off  </li>
    </ul>
    Keywords: <i>ZONAL</i>
  */
  m_zonalMethod = Context::getSolverProperty<MString>("zonalMethod", couplerId(), AT_, &m_zonalMethod);
  /*! \page propertyPage1
    \section cutOffMethod
    <code>MBool zonalDualTimeStepping </code>\n
    default = <code>"false"</code>\n \n
    Allow for different time-Steps in the different Fv-Mb zones.
    Keywords: <i>ZONAL</i>
  */
  m_zonalDualTimeStepping =
      Context::getSolverProperty<MBool>("zonalDualTimeStepping", couplerId(), AT_, &m_zonalDualTimeStepping);

  if(m_zonalMethod == "P") {
    /*! \page propertyPage1
    \section outsideDefault
    <code>MFloat* CouplerFvMbZonal::zonalCoordinate </code>\n
    Set the zonal coordinate
    <ul>
    <li>any float number </li>
    </ul>
    Keywords: <i>ZONAL, FVMB-COUPLING</i>
    */
    m_zonalCoordinate = Context::getSolverProperty<MFloat>("zonalCoordinate", couplerId(), AT_);

    /*! \page propertyPage1
   \section outsideDefault
   <code>MInt* CouplerFvMbZonal::zonalDir </code>\n
   Set the zonal direction, must be a cartesian direction!
   <ul>
   <li>0 : x-dir </li>
   <li>1 : y-dir </li>
   <li>2 : z-dir </li>
   </ul>
   Keywords: <i>ZONAL, FVMB-COUPLING</i>
   */
    m_zonalDir = Context::getSolverProperty<MInt>("zonalDir", couplerId(), AT_);

    if(upStream().domainId() == 0) {
      cerr << "Using cartesian zonal Method with plane-cut at " << m_zonalCoordinate << " in direction " << m_zonalDir
           << endl;
    }

  } else if(m_zonalMethod == "B") {
  } else {
    mTerm(1, AT_, "ZonalMethod not implemented yet!");
  }
}


/** \brief create first zonal mapping for the zonal exchange
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::finalizeCouplerInit() {
  createZonalMapping();

  ASSERT(upStream().a_timeStepComputationInterval() == downStream().a_timeStepComputationInterval(), "");
}

/** \brief update zonal mapping after adaptation
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::finalizeAdaptation(MInt solverLoopId) {
  if(solverLoopId >= m_upStreamId && solverLoopId >= m_downStreamId) {
    createZonalMapping();
    unifyTimeStep();
  }
}

/** \brief update zonal mapping after balance
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::finalizeBalance(MInt solverLoopId) {
  if(solverLoopId >= m_upStreamId && solverLoopId >= m_downStreamId) {
    createZonalMapping();
  }
}

/** \brief exchange zonal variables each runge-Kutta Step
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::subCouple(MInt recipeStep, MInt solverId,
                                               std::vector<MBool>& /*solverCompleted*/) {
  // only during FVMb computation
  if(recipeStep == 1 && (solverId == m_upStreamId || solverId == m_downStreamId)) {
    m_RKStep++;

    if(m_RKStep == 2 * m_noRKSteps) {
      // complete exchange during the last runge-Kutta Step
      exchangeZonalValues(-1);
    } else {
      // only update from the solver who just made its timeStep!
      exchangeZonalValues(solverId);
    }

    // check the timeStep contraint at the last-RK-step of both solvers
    // if the timeStep might have changed!
    if(m_RKStep == 2 * m_noRKSteps && upStream().a_timeStepComputationInterval() > 0
       && globalTimeStep % upStream().a_timeStepComputationInterval() == 0) {
      unifyTimeStep();
    }
  }
}

/** \brief exchange zonal variables before each runge-Kutta Step
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::preCouple(MInt solverStep) {
  // past the LS-solver
  // before the FvMb-solver
  if(solverStep == 1) {
    // updateZonalMapping();
    exchangeZonalValues(-1);
    m_RKStep = 0;
  }
}

/** \brief unify the time-Step
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::postCouple(MInt solverStep) {
  if(solverStep == 1) {
    if(upStream().isActive() && downStream().isActive()
       && fabs(upStream().timeStep() - downStream().timeStep()) > upStream().m_eps) {
      cerr << setprecision(16) << "UpStream-TimeStep " << upStream().timeStep() << " DownStream-TimeStep "
           << upStream().timeStep() << endl;
      mTerm(1, AT_, "Time in Fv-Mb-Solvers differs!.");
    }
  }
}

/** \brief create the cell mapping which allows for a faster exchange during each
 *         RK-Step. This requires little additional memory.
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::createZonalMapping() {
  TRACE();

  m_upDown.clear();
  m_downUp.clear();

  m_upStreamOffset = 0;
  m_downStreamOffset = 0;

  if(upStream().isActive()) {
    if(m_zonalMethod == "P") {
      // first create mapping for valid grid cells
      // this mapping will only change duing adaptation or balance
      for(MInt upId = 0; upId < upStream().c_noCells(); upId++) {
        if(upStream().a_coordinate(upId, m_zonalDir) > m_zonalCoordinate) {
          MBool lastLayer = false;
          if(!upStream().a_hasNeighbor(upId, m_zonalDir * 2 + 1, false) && !upStream().a_isHalo(upId)) {
            const MInt parentId = upStream().c_parentId(upId);
            if(parentId > -1 && !upStream().a_hasNeighbor(parentId, m_zonalDir * 2 + 1, false)) {
              lastLayer = true;
            }
          }
          const MInt downId = up2downId(upId, lastLayer);
          if(downId > -1) {
            m_downUp.push_back(make_pair(downId, upId));
          }
        }
      }
    } else if(m_zonalMethod == "B") {
    }
    m_downStreamOffset = m_downUp.size();
  }

  if(downStream().isActive()) {
    if(m_zonalMethod == "P") {
      // first create mapping for valid grid cells
      // this mapping will only change duing adaptation or balance
      for(MInt downId = 0; downId < downStream().c_noCells(); downId++) {
        if(downStream().a_coordinate(downId, m_zonalDir) < m_zonalCoordinate) {
          MBool lastLayer = false;
          if(!downStream().a_hasNeighbor(downId, m_zonalDir * 2, false) && !downStream().a_isHalo(downId)) {
            const MInt parentId = downStream().c_parentId(downId);
            if(parentId > -1 && !downStream().a_hasNeighbor(parentId, m_zonalDir * 2, false)) {
              lastLayer = true;
            }
          }
          const MInt upId = down2upId(downId, lastLayer);
          if(upId > -1) {
            m_upDown.push_back(make_pair(upId, downId));
          }
        }
      }
    } else if(m_zonalMethod == "B") {
    }
    m_upStreamOffset = m_upDown.size();
  }
}

/** \brief update the zonal mapping before the FVMb-timeStep as the
 *         bndryCell order might have changed!
 *   NOTE: currently unused, as its not necessary, the bndry-ghostCell values are
 *          updated inside the solutionStep.
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::updateZonalMapping() {
  TRACE();

  // resize to grid-cell offset
  m_upDown.resize(m_upStreamOffset);
  m_downUp.resize(m_downStreamOffset);

  // re-build bndryCell mapping
  if(upStream().isActive()) {
    if(m_zonalMethod == "P") {
      for(MInt upId = upStream().c_noCells(); upId < upStream().a_noCells(); upId++) {
        if(upStream().a_coordinate(upId, m_zonalDir) > m_zonalCoordinate) {
          MBool lastLayer = false;
          const MInt associatedId = upStream().getAssociatedInternalCell(upId);
          if(!upStream().a_hasNeighbor(associatedId, m_zonalDir * 2 + 1, false) && !upStream().a_isHalo(associatedId)) {
            const MInt parentId = upStream().c_parentId(associatedId);
            if(parentId > -1 && !upStream().a_hasNeighbor(parentId, m_zonalDir * 2 + 1, false)) {
              lastLayer = true;
            }
          }
          const MInt downId = up2downId(upId, lastLayer);
          if(downId > -1) {
            m_upDown.push_back(make_pair(upId, downId));
          }
        }
      }
    } else if(m_zonalMethod == "B") {
    }
  }


  if(downStream().isActive()) {
    if(m_zonalMethod == "P") {
      for(MInt downId = downStream().c_noCells(); downId < downStream().a_noCells(); downId++) {
        if(downStream().a_coordinate(downId, m_zonalDir) < m_zonalCoordinate) {
          MBool lastLayer = false;
          const MInt associatedId = downStream().getAssociatedInternalCell(downId);
          if(!downStream().a_hasNeighbor(associatedId, m_zonalDir * 2, false) && !downStream().a_isHalo(associatedId)) {
            const MInt parentId = downStream().c_parentId(associatedId);
            if(parentId > -1 && !downStream().a_hasNeighbor(parentId, m_zonalDir * 2, false)) {
              lastLayer = true;
            }
          }
          const MInt upId = down2upId(downId, lastLayer);
          if(upId > -1) {
            m_upDown.push_back(make_pair(downId, upId));
          }
        }
      }
    } else if(m_zonalMethod == "B") {
    }
  }
}


/** \brief exchange zonal variables around the zonal coordinate
 *          from upstream to downstream solver and the other way around!
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::exchangeZonalValues(const MInt solverId) {
  TRACE();

  if(m_upDown.empty()) {
    ASSERT(m_downUp.empty(), "");
    return;
  }

  // transfer data from up-stream to down-stream
  if(solverId == m_upStreamId || solverId == -1) {
    for(auto it = m_upDown.begin(); it != m_upDown.end(); it++) {
      const MInt upId = it->first;
      const MInt downId = it->second;

      for(MInt var = 0; var < upStream().m_sysEqn.CV->noVariables; var++) {
        downStream().a_variable(downId, var) = upStream().a_variable(upId, var);
        downStream().a_pvariable(downId, var) = upStream().a_pvariable(upId, var);
      }
    }
  }

  // transfer data from down-stream to up-stream
  if(solverId == m_downStreamId || solverId == -1) {
    for(auto it = m_downUp.begin(); it != m_downUp.end(); it++) {
      const MInt downId = it->first;
      const MInt upId = it->second;

      for(MInt var = 0; var < downStream().m_sysEqn.CV->noVariables; var++) {
        upStream().a_variable(upId, var) = downStream().a_variable(downId, var);
        upStream().a_pvariable(upId, var) = downStream().a_pvariable(downId, var);
      }
    }
  }
}


/** \brief conversion from upStream solverId to the downStream solverId
 *         NOTE: also handles bndry-ghost cell conversion!
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
MInt CouplerFvMbZonal<nDim, SysEqn>::up2downId(MInt upCell, MBool lastLayer) {
  TRACE();

  if(!downStream().isActive()) return -1;

  // to hanlde bndry ghostCells and others as well
  //   1: regular grid cell
  // < 1: bndry-ghost cell (index represents matching surface)
  MInt cellType = 1;
  if(upCell >= upStream().c_noCells()) {
    const MInt gridCell = upStream().getAssociatedInternalCell(upCell);
    if(upStream().a_isBndryGhostCell(upCell)) {
      // find the matching surface of the grid cell to which the upCell is the bndry-ghostcell
      const MInt bndryId = upStream().a_bndryId(gridCell);
      MInt srfc = -1;
      for(MInt surf = 0; surf < upStream().m_bndryCells->a[bndryId].m_noSrfcs; surf++) {
        if(upStream().m_bndryCells->a[bndryId].m_srfcVariables[surf]->m_ghostCellId == upCell) {
          srfc = surf;
          break;
        }
      }
      ASSERT(srfc >= 0, "");
      cellType = -srfc;
    } else if(upStream().a_hasProperty(upCell, FvCell::IsSplitChild)
              || upStream().a_hasProperty(upCell, FvCell::IsSplitClone)) {
      cellType = 2;
    }
    upCell = gridCell;
  }

  upStream().assertValidGridCellId(upCell);

  const MInt gridId = upStream().grid().tree().solver2grid(upCell);
  ASSERT(upStream().grid().solverFlag(gridId, m_upStreamId), "");

  if(!downStream().grid().solverFlag(gridId, m_downStreamId)) return -1;

  MInt downId = downStream().grid().tree().grid2solver(gridId);

  if(downId > 0) {
    ASSERT(downStream().a_level(downId) == upStream().a_level(upCell), "");
  }

  if(cellType == 1 || downId < 0) return downId;

  if(!lastLayer) {
    ASSERT(downStream().a_isBndryCell(upCell), "");
  } else {
    return -1;
  }

  if(cellType < 1) {
    // bndry-ghost cell
    const MInt srfc = -cellType;
    const MInt bndryId = downStream().a_bndryId(downId);
    downId = downStream().m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_ghostCellId;
    return downId;

  } else {
    // split child
    mTerm(1, AT_, "Not yet implemented for splitchilds");
  }

  return -1;
}

/** \brief conversion from downStream to upStream cellId
 *         NOTE: also handles bndry-ghost cells
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
MInt CouplerFvMbZonal<nDim, SysEqn>::down2upId(MInt downCell, MBool lastLayer) {
  TRACE();

  if(!upStream().isActive()) return -1;

  // to hanlde bndry ghostCells and others as well
  //   1: regular grid cell
  // < 1: bndry-ghost cell (index represents matching surface)
  MInt cellType = 1;
  if(downCell >= downStream().c_noCells()) {
    const MInt gridCell = downStream().getAssociatedInternalCell(downCell);
    if(downStream().a_isBndryGhostCell(downCell)) {
      // find the matching surface of the grid cell to which the upCell is the bndry-ghostcell
      const MInt bndryId = downStream().a_bndryId(gridCell);
      MInt srfc = -1;
      for(MInt surf = 0; surf < downStream().m_bndryCells->a[bndryId].m_noSrfcs; surf++) {
        if(downStream().m_bndryCells->a[bndryId].m_srfcVariables[surf]->m_ghostCellId == downCell) {
          srfc = surf;
          break;
        }
      }
      ASSERT(srfc >= 0, "");
      cellType = -srfc;
    } else if(downStream().a_hasProperty(downCell, FvCell::IsSplitChild)
              || downStream().a_hasProperty(downCell, FvCell::IsSplitClone)) {
      cellType = 2;
    }
    downCell = gridCell;
  }

  downStream().assertValidGridCellId(downCell);

  const MInt gridId = downStream().grid().tree().solver2grid(downCell);
  ASSERT(downStream().grid().solverFlag(gridId, m_downStreamId), "");

  if(!upStream().grid().solverFlag(gridId, m_upStreamId)) return -1;

  MInt upCell = upStream().grid().tree().grid2solver(gridId);

  if(upCell > 0) {
    ASSERT(downStream().a_level(downCell) == upStream().a_level(upCell), "");
  }

  if(cellType == 1 || upCell < 0) return upCell;

  if(!lastLayer) {
    ASSERT(upStream().a_isBndryCell(upCell), "");
  } else {
    return -1;
  }

  if(cellType < 1) {
    // bndry-ghost cell
    const MInt srfc = -cellType;
    const MInt bndryId = upStream().a_bndryId(upCell);
    upCell = upStream().m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_ghostCellId;
    return upCell;

  } else {
    // split child
    mTerm(1, AT_, "Not yet implemented for splitchilds");
  }

  return -1;
}


/** \brief set the same timeStep both FvMb-solvers
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvMbZonal<nDim, SysEqn>::unifyTimeStep() {
  const MFloat invalidTimeStep = std::numeric_limits<MFloat>::max();
  MFloat timeUp = invalidTimeStep;
  MFloat timeDown = invalidTimeStep;
  MInt maxLevelUp = -1;
  MInt maxLevelDown = -1;
  MFloat newTimeUp = invalidTimeStep;
  MFloat newTimeDown = invalidTimeStep;

  // determine the new timeStep on ranks which have both up-Stream and downStream solvers!
  if(upStream().isActive() && downStream().isActive()) {
    timeUp = upStream().timeStep(true);
    maxLevelUp = upStream().maxLevel();
    timeDown = downStream().timeStep(true);
    maxLevelDown = downStream().maxLevel();

    MInt levelDiff = maxLevelDown - maxLevelUp;
    if(levelDiff == 0) {
      newTimeUp = mMin(timeUp, timeDown);
      newTimeDown = newTimeUp;
    } else {
      ASSERT(m_zonalDualTimeStepping, "");
      if(levelDiff == 1) {
        // downStream has higher level => lower timeStep!
        newTimeDown = mMin(newTimeDown, newTimeUp / 2.0);
        newTimeUp = 2.0 * newTimeDown;
      } else if(levelDiff == -1) {
        // upStream has the higher level => lower timeStep!
        newTimeUp = mMin(newTimeUp, newTimeDown / 2.0);
        newTimeDown = 2.0 * newTimeUp;
      } else {
        mTerm(1, AT_, "Zonal dual-timeStepping currently only implemented for 1 level difference!");
      }
    }
  }

  if(downStream().isActive()) {
    MPI_Allreduce(MPI_IN_PLACE, &newTimeDown, 1, MPI_DOUBLE, MPI_MIN, downStream().mpiComm(), AT_, "MPI_IN_PLACE",
                  "newTimeDown");
  }
  if(upStream().isActive()) {
    MPI_Allreduce(MPI_IN_PLACE, &newTimeUp, 1, MPI_DOUBLE, MPI_MIN, upStream().mpiComm(), AT_, "MPI_IN_PLACE",
                  "newTimeUp");
  }

  if(upStream().isActive()) {
    upStream().forceTimeStep(newTimeUp);
  }
  if(downStream().isActive()) {
    downStream().forceTimeStep(newTimeDown);
  }
}

// Explicit template instantiations
template class CouplerFvMbZonal<2, FvSysEqnNS<2>>;
template class CouplerFvMbZonal<3, FvSysEqnNS<3>>;
template class CouplerFvMbZonal<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class CouplerFvMbZonal<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplerFvMbZonal<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class CouplerFvMbZonal<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplerFvMbZonal<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplerFvMbZonal<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
