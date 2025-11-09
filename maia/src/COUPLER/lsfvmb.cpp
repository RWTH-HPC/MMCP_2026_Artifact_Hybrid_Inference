// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

//#define COUPLING_DEBUG_

#include "lsfvmb.h"

#include <algorithm>
#include <stack>
#include <vector>
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "UTIL/kdtree.h"
#include "coupling.h"
#include "globals.h"
#include "globalvariables.h"

using namespace std;

template <MInt nDim, class SysEqn>
LsFvMb<nDim, SysEqn>::LsFvMb(const MInt couplingId, LsSolver* ls, FvMbSolver* fvMb)
  : Coupling(couplingId), CouplingLS<nDim>(couplingId, ls), CouplingFvMb<nDim, SysEqn>(couplingId, fvMb) {
  TRACE();

  initData();

  readProperties();
}


/** \brief Initialize coupling-class-specific Data
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::initData() {
  TRACE();

  // solver-specific data:
  m_fvSolverId = fvMbSolver().m_solverId;
  m_lsSolverId = lsSolver().m_solverId;

  if(lsSolver().m_closeGaps) {
    mAlloc(m_hadGapCells, lsSolver().m_noGapRegions, "m_hadGapCells", 0, AT_);
    mAlloc(m_hasGapCells, lsSolver().m_noGapRegions, "m_hasGapCells", 0, AT_);
  }

  m_noRfJumps = 0;
  for(MInt set = 0; set < lsSolver().m_noSets; set++) {
    const MInt levelDif = fvMbSolver().m_lsCutCellLevel[set] - lsSolver().a_maxGCellLevel(set);
    m_noRfJumps = mMax(m_noRfJumps, levelDif);
  }

  const MInt maxlevelDif = fvMbSolver().maxRefinementLevel() - lsSolver().maxRefinementLevel();
  m_noRfJumps = mMax(m_noRfJumps, maxlevelDif);

  if(lsSolver().m_maxLevelChange) m_noRfJumps = 0;

  m_log << " Ls - FvMb Coupler detected " << m_noRfJumps << " level-jumps" << endl;


  fvMbSolver().m_maxLsValue = lsSolver().m_outsideGValue;

  if(fvMbSolver().m_levelSetRans) {
    mAlloc(fvMbSolver().m_levelSetValues, lsSolver().m_noSets, "fvMbSolver().m_levelSetValues", AT_);
  }
}

/** \brief Checks property-data which is read in by both ls-and Fv-Solver
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::checkProperties() {
  TRACE();

  ASSERT(lsSolver().m_noEmbeddedBodies == fvMbSolver().m_noEmbeddedBodies, "");
  ASSERT(lsSolver().m_constructGField == fvMbSolver().m_constructGField, "");
  ASSERT(lsSolver().m_noGapRegions == fvMbSolver().m_noGapRegions, "");

  if(!m_allowLsInterpolation) {
    ASSERT(lsSolver().a_maxGCellLevel(0) == fvMbSolver().maxRefinementLevel() || lsSolver().m_maxLevelChange, "");
    ASSERT(m_noRfJumps == 0,
           to_string(lsSolver().a_maxGCellLevel(0)) + " " + to_string(fvMbSolver().maxRefinementLevel()));
  } else {
    ASSERT(m_noRfJumps > 0, "");
    if(lsSolver().domainId() == 0) {
      cerr << "Allowing ls interpolation with " << m_noRfJumps << " level jumps!" << endl;
    }
    ASSERT(lsSolver().m_gShadowWidth <= fvMbSolver().m_bandWidth[fvMbSolver().maxRefinementLevel() - 1], "");
  }

  if(!lsSolver().m_maxLevelChange) {
    ASSERT(lsSolver().a_maxGCellLevel(0) == fvMbSolver().m_lsCutCellBaseLevel, "");
  }

  ASSERT(lsSolver().m_closeGaps == fvMbSolver().m_closeGaps, "");
  ASSERT(lsSolver().m_noBodiesInSet == fvMbSolver().m_noBodiesInSet, "");
  ASSERT(lsSolver().m_bodyToSetTable == fvMbSolver().m_bodyToSetTable, "");
  ASSERT(lsSolver().m_setToBodiesTable == fvMbSolver().m_setToBodiesTable, "");
  ASSERT(lsSolver().m_startSet == fvMbSolver().m_startSet, "");
  ASSERT(lsSolver().m_maxNoSets == fvMbSolver().m_noLevelSetsUsedForMb, "");
  ASSERT(lsSolver().m_noSets == fvMbSolver().m_noSets, "");

  // If fvMbSolver is inactive it does not have any cells
  if(fvMbSolver().isActive()) {
    // Check that the fv-solver-cell-count is correct!
    ASSERT(fvMbSolver().a_noCells(), "");
  }

  if(!fvMbSolver().m_constructGField) ASSERT(g_multiSolverGrid, "");

  ASSERT(lsSolver().m_engineSetup == fvMbSolver().m_engineSetup, "");
  if(lsSolver().m_engineSetup && lsSolver().isActive() && fvMbSolver().isActive()) {
    ASSERT(fabs(lsSolver().crankAngle(lsSolver().m_time, 0) - fvMbSolver().crankAngle(fvMbSolver().m_physicalTime, 0))
               < fvMbSolver().m_eps,
           "");
  }
}

/** \fn void LsFvMb<nDim, SysEqn>::readProperties()
 * \brief reads lsfvmb-coupling-specific data
 * \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::readProperties() {
  TRACE();

  /*! \page propertiesCoupling
  \section outsideDefault
  <code>MBool* LsFvMb::outsideDefault </code>\n
  default = <code>true</code>\n \n
  A trigger which determines the default levelset values of a fv-cell if no connection a ls-cell
  can be found! This is the case if the ls-domain is smaller than the fv-domain!
  Possible values are:
  <ul>
  <li>true: possitive outsideGValue </li>
  <li>false: negative outsideGValue </li>
  </ul>
  Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
  */

  m_outsideDefault = true;
  if(Context::propertyExists("outsideDefault", m_fvSolverId)) {
    m_outsideDefault = Context::getSolverProperty<MBool>("outsideDefault", m_fvSolverId, AT_, &m_outsideDefault);
  }

  /*! \page propertiesCoupling
  \section allowLsInterpolation
  <code>MBool* LsFvMb::allowLsInterpolation </code>\n
  default = <code>false</code>\n \n
  Trigger whether the level set values are interpolated between LS and FV-MB grids.
  Interpolation is required when the FvMb Solver has a higher resolution at the boundary,
  otherwise the version without interpolation is faster!
  <ul>
  <li>true: enable grid interpolation </li>
  <li>false: disable grid interpolation </li>
  </ul>
  Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
  */
  m_allowLsInterpolation = false;
  m_allowLsInterpolation =
      Context::getSolverProperty<MBool>("allowLsInterpolation", m_fvSolverId, AT_, &m_allowLsInterpolation);
}


/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *    \author  Tim Wegmann
 */

template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::testCoupling() {
  TRACE();

  if(fvMbSolver().m_constructGField) return;

  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    ASSERT(fvMbSolver().grid().tree().solver2grid(fc) >= 0, "");
    ASSERT(fvMbSolver().grid().solverFlag(fvMbSolver().grid().tree().solver2grid(fc), m_fvSolverId), "");
#ifdef COUPLING_DEBUG_
    if(lsSolver().grid().solverFlag(fvMbSolver().grid().tree().solver2grid(fc), m_lsSolverId)) {
      ASSERT(ls2fvId(fv2lsId(fc)) == fc,
             to_string(fc) + " " + to_string(fv2lsId(fc)) + " " + to_string(ls2fvId(fv2lsId(fc))));
    }
#endif
  }
  for(MInt cellId = 0; cellId < a_noLsCells(); cellId++) {
    ASSERT(lsSolver().grid().tree().solver2grid(cellId) >= 0, "");
    ASSERT(lsSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_lsSolverId), "");
#ifdef COUPLING_DEBUG_
    if(fvMbSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_fvSolverId)) {
      ASSERT(fv2lsId(ls2fvId(cellId)) == cellId,
             to_string(cellId) + " " + to_string(ls2fvId(cellId)) + " " + to_string(fv2lsId(ls2fvId(cellId))));
    }
#endif
  }

  if(!lsSolver().m_levelSetMb) return;


#ifdef COUPLING_DEBUG_
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    for(MInt dir = 0; dir < nDim; dir++) {
      if(lsSolver().grid().solverFlag(fvMbSolver().grid().tree().solver2grid(fc), m_lsSolverId)) {
        ASSERT(abs(fvMbSolver().c_coordinate(fc, dir) - lsSolver().c_coordinate(fv2lsId(fc), dir)) < 0.00000001, "");
      }
    }
  }
#endif
}

/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *    \author  Tim Wegmann
 */

template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::testLsValues() {
  TRACE();

  if(fvMbSolver().m_constructGField) return;

  if(!fvMbSolver().m_levelSetMb) return;

#ifdef COUPLING_DEBUG_
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
      if(lsSolver().grid().solverFlag(fvMbSolver().grid().tree().solver2grid(fc), m_lsSolverId)) {
        ASSERT(abs(a_levelSetFunctionG(fv2lsId(fc), set) - fvMbSolver().a_levelSetValuesMb(fc, set)) < 0.00000001, "");
      }
    }
  }
#endif
}

/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *    \author  Tim Wegmann
 */

template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::testGapProperty() {
  TRACE();

  if(!lsSolver().m_closeGaps) return;
  if(lsSolver().m_levelSetMb && fvMbSolver().m_constructGField) return;

#ifdef COUPLING_DEBUG_
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    if(fvMbSolver().a_isGapCell(fc) || fvMbSolver().a_wasGapCell(fc)) {
      MInt gc = lsSolver().grid().tree().grid2solver(fvMbSolver().grid().tree().solver2grid(fc));
      ASSERT(gc > -1, "");
      ASSERT(lsSolver().a_potentialGapCellClose(gc) > 0
                 && lsSolver().a_potentialGapCellClose(gc) <= fvMbSolver().m_noEmbeddedBodies,
             "");
    }
  }
#endif
}


/** \brief Sets the Levelset-Values in fvSolver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferLevelSetValues() {
  TRACE();

  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    fvMbSolver().m_levelSetValues[s].resize(fvMbSolver().a_noCells());
  }

  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    for(MInt fc = 0; fc < a_noFvCells(); fc++) {
      MInt lsId = fv2lsIdParent(fc);
      fvMbSolver().a_levelSetFunction(fc, s) = a_outsideGValue();
      if(lsId > -1) {
        if(fvMbSolver().a_hasProperty(fc, FvCell::IsMovingBnd)) {
          // Interplate levelset for cut-cells
          // direct interpolation from ls to fv for only 1 level difference
          MFloat point[3] = {fvMbSolver().a_coordinate(fc, 0), fvMbSolver().a_coordinate(fc, 1),
                             fvMbSolver().a_coordinate(fc, 2)};
          fvMbSolver().a_levelSetFunction(fc, s) = interpolateLevelSet(lsId, point, s);
        } else {
          fvMbSolver().a_levelSetFunction(fc, s) = lsSolver().a_levelSetFunctionG(fv2lsIdParent(fc), s);
        }
      }
    }
  }
}


/** \brief transfers the LevelSetValues for all cells from the levelset to the moving boundary Part
 *    \author Christoph Siewert, Tim Wegmann
 */

template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferLevelSetFieldValues(MBool exchangeLVS) {
  TRACE();
  // Nothing to transfer if fvMb is inactive!
  if(!fvMbSolver().isActive()) return;

  if(exchangeLVS) {
    // lsSolver().startLoadTimer(AT_);
    // lsSolver().checkHaloCells();
    // lsSolver().stopLoadTimer(AT_);
  }

#ifdef COUPLING_DEBUG_
  testCoupling();
#endif

  std::map<MInt, MInt> interpolationParents;
  interpolationParents.clear();

  // Set outside default if lssolver is not active otherwise map
  if(lsSolver().isActive()) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt cellId = 0; cellId < a_noFvCells(); cellId++) {
      // reset values with invalid/defaults
      for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
        if(m_outsideDefault) {
          fvMbSolver().a_levelSetValuesMb(cellId, set) = a_outsideGValue();
        } else {
          fvMbSolver().a_levelSetValuesMb(cellId, set) = -a_outsideGValue();
        }
        fvMbSolver().a_associatedBodyIds(cellId, set) = -1;
      }
      const MInt gCellId = fv2lsId(cellId);
      // direct transfer of matching cells
      if(gCellId > -1) {
        for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
          fvMbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellId, set);
          fvMbSolver().a_levelSetValuesMb(cellId, set) = a_levelSetFunctionG(gCellId, set);
        }
      } else {
        // ls and fv grid differs

        if(!g_multiSolverGrid) {
          ASSERT(cellId >= a_noFvGridCells(), "");
        } else {                            // multiSolverGrid, where the level and the grid-extension might differ!
          if(cellId <= a_noFvGridCells()) { // grid-cell!
            const MInt gCellParent = fv2lsIdParent(cellId);

            if(gCellParent > -1) { // just level difference

              if(m_noRfJumps == 0) {
                // same maxRefinementLevel, meaning just a difference due to different refinement-widths
                for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                  fvMbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                  // ASSERT(!lsSolver().a_inBandG( gCellParent ,  set ), "");
                  fvMbSolver().a_levelSetValuesMb(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                }
              } else {
                // ls-interpolation towards higher fv-level
                ASSERT(m_allowLsInterpolation, "");
                if(fvMbSolver().a_isHalo(cellId) || lsSolver().a_isHalo(gCellParent)) continue;

                // simple transfer if not a band-cell
                if(lsSolver().m_buildCollectedLevelSetFunction && !lsSolver().a_inBandG(gCellParent, 0)
                   && lsSolver().a_level(gCellParent) <= lsSolver().a_maxGCellLevel(0)) {
                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    fvMbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                    fvMbSolver().a_levelSetValuesMb(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                  }
                  continue;
                }

                const MInt levelDifference = fvMbSolver().c_level(cellId) - lsSolver().a_level(gCellParent);
                if(levelDifference == 1) {
                  // direct interpolation from ls to fv for only 1 level difference
                  MFloat point[3] = {fvMbSolver().c_coordinate(cellId, 0), fvMbSolver().c_coordinate(cellId, 1),
                                     fvMbSolver().c_coordinate(cellId, 2)};
                  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                  MInt position = lsSolver().setUpLevelSetInterpolationStencil(gCellParent, interpolationCells, point);

                  // interpolate levelset data
                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    fvMbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                    if(position < 0) {
                      // no valid interpolation stencil found
                      fvMbSolver().a_levelSetValuesMb(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                    } else { // interpolation for all sets
                      const MFloat phi = lsSolver().interpolateLevelSet(interpolationCells, point, set);
                      fvMbSolver().a_levelSetValuesMb(cellId, set) = phi;
                    }
                  }
                  // build collected levelset data
                  buildCollectedLevelSet(cellId);

                } else {
                  // add cells based on which the fv-interpolation needs to start
                  // that is cells for which the level difference is above 1
                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    fvMbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                  }

                  MInt parent = cellId;
                  for(MInt i = 0; i < levelDifference - 1; i++) {
                    parent = fvMbSolver().c_parentId(parent);
                  }
                  ASSERT(parent > -1, "");
                  // parent of this parent has a direct link, meaning that the values for the parent
                  // are set above!
                  ASSERT(fv2lsId(fvMbSolver().c_parentId(parent)) == gCellParent, "");
                  interpolationParents.insert(make_pair(parent, fvMbSolver().a_level(parent)));
                }
              }
            }
          }
        }
      }
    }
  } else {
// lssolver is inactive
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(MInt cellId = 0; cellId < a_noFvCells(); cellId++) {
      for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
        if(m_outsideDefault) {
          fvMbSolver().a_levelSetValuesMb(cellId, set) = a_outsideGValue();
        } else {
          fvMbSolver().a_levelSetValuesMb(cellId, set) = -a_outsideGValue();
        }
        fvMbSolver().a_associatedBodyIds(cellId, set) = -1;
      }
    }
  }

  // now interpolate through all childs
  // iterate downwards level by level and exchange values in between!
  if(m_allowLsInterpolation) {
    std::list<MInt> newInterpolationParents;

    for(MInt lvl = lsSolver().maxUniformRefinementLevel() + 1; lvl < fvMbSolver().maxRefinementLevel(); lvl++) {
      fvMbSolver().exchangeLevelSetData();

      MInt noInterpolationCells = (MInt)interpolationParents.size();
      MPI_Allreduce(MPI_IN_PLACE, &noInterpolationCells, 1, MPI_INT, MPI_MAX, fvMbSolver().mpiComm(), AT_, "INPLACE",
                    "noInterpolationCells");
      if(noInterpolationCells == 0) break;

      if(interpolationParents.empty()) continue;

      ASSERT(m_noRfJumps > 1, "");

      newInterpolationParents.clear();

      for(auto it = interpolationParents.begin(); it != interpolationParents.end(); it++) {
        const MInt parent = it->first;
        const MInt level = it->second;
        if(level > lvl) {
          newInterpolationParents.push_back(parent);
        } else {
          for(MInt c = 0; c < fvMbSolver().grid().m_maxNoChilds; c++) {
            const MInt childId = fvMbSolver().c_childId(parent, c);
            if(childId < 0) {
              continue;
            }
            interpolateLsFV(parent, childId);
            if(!fvMbSolver().c_isLeafCell(childId)) {
              newInterpolationParents.push_back(childId);
            }
          }
        }
      }
      interpolationParents.clear();
      for(auto it = newInterpolationParents.begin(); it != newInterpolationParents.end(); it++) {
        const MInt cellId = (*it);
        interpolationParents.insert(make_pair(cellId, fvMbSolver().a_level(cellId)));
      }
    }
  }

  if(exchangeLVS) fvMbSolver().exchangeLevelSetData();

  for(MUint sc = 0; sc < fvMbSolver().m_splitCells.size(); sc++) {
    const MInt cellId = fvMbSolver().m_splitCells[sc];
    for(MUint ssc = 0; ssc < fvMbSolver().m_splitChilds[sc].size(); ssc++) {
      const MInt splitChildId = fvMbSolver().m_splitChilds[sc][ssc];
      for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
        fvMbSolver().a_associatedBodyIds(splitChildId, set) = fvMbSolver().a_associatedBodyIds(cellId, set);
        fvMbSolver().a_levelSetValuesMb(splitChildId, set) = fvMbSolver().a_levelSetValuesMb(cellId, set);
      }
    }
  }

  testLsValues();
}

/** \brief interpolate levelset values on the fv-grid
 *
 * \author Tim Wegmann
 * \date March 2020
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::interpolateLsFV(const MInt from, const MInt to) {
  TRACE();

  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MFloat point[3] = {fvMbSolver().c_coordinate(to, 0), fvMbSolver().c_coordinate(to, 1),
                     fvMbSolver().c_coordinate(to, 2)};

  std::function<MBool(const MInt, const MInt)> alwaysTrue = [&](const MInt, const MInt) { return true; };

  MBool backup = fvMbSolver().m_deleteNeighbour;
  fvMbSolver().m_deleteNeighbour = false;
  const MInt position = fvMbSolver().setUpInterpolationStencil(from, interpolationCells, point, alwaysTrue, false);
  fvMbSolver().m_deleteNeighbour = backup;

  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
    if(position < 0) {
      fvMbSolver().a_levelSetValuesMb(to, set) = fvMbSolver().a_levelSetValuesMb(from, set);
    } else {
      MFloat phi = interpolateLevelSetMb(interpolationCells, point, set);
      fvMbSolver().a_levelSetValuesMb(to, set) = phi;
    }
  }
  buildCollectedLevelSet(to);
}

///  interpolates the levelSet value further on fv-cells
/// \author Jannik Borgelt
/// \date 2020-04-01
template <MInt nDim, class SysEqn>
MFloat LsFvMb<nDim, SysEqn>::interpolateLevelSet(MInt cellId, MFloat* point, MInt set) {
  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MInt position = 0;

  position = lsSolver().setUpLevelSetInterpolationStencil(cellId, interpolationCells, point);

  // Interpolate level set
  if(position > -1) {
    return lsSolver().interpolateLevelSet(interpolationCells, point, set);
  } else {
    return lsSolver().a_levelSetFunctionG(cellId, set);
  }
}

///  interpolates the levelSet value further on fv-cells
/// \author Tim Wegmann
/// \date 2020-04-01
template <MInt nDim, class SysEqn>
MFloat LsFvMb<nDim, SysEqn>::interpolateLevelSetMb(MInt* interpolationCells, MFloat* point, const MInt set) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt refSet) {
    return static_cast<MFloat>(fvMbSolver().a_levelSetValuesMb(cellId, refSet));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(fvMbSolver().c_coordinate(cellId, id));
  };

  return fvMbSolver().template interpolateFieldData<true>(&interpolationCells[0], &point[0], set, scalarField,
                                                          coordinate);
}

//---------------------------------------------------------------------------

/** \brief Sets the gap-cell-property
 *         mode 1 : regular call each timeStep (advance is/was-GapCell)
 *         mode 0 : initialisation and after balance (wasGapCell = isGapCell)!
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferGapCellProperty(MInt mode) {
  TRACE();

  if(globalTimeStep < 0) {
    mode = 0;
  }

  // return if lssolver is not active
  // if(!(lsSolver().isActive())) return;

  if(!lsSolver().m_closeGaps) return;
  if(lsSolver().m_noGapRegions <= 0) return;
  if(lsSolver().m_levelSetMb && fvMbSolver().m_constructGField) return;

#ifdef COUPLING_DEBUG_
  testCoupling();
#endif

  // 1) reset has/had Gap-Cells
  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    m_hadGapCells[region] = 0;
    m_hasGapCells[region] = 0;
  }

  // 2) set a_wasGapCell, m_hadGapCells and initialize a_isGapCell
  for(MInt fc = 0; fc < fvMbSolver().noInternalCells(); fc++) {
    if(fvMbSolver().a_isGapCell(fc)) {
      MInt lsId = fv2lsId(fc);
      if(lsId < 0) {
        if(fc > fvMbSolver().c_noCells()) continue;
        lsId = fv2lsIdParent(fc);
      }
      ASSERT(lsId > -1, "");
      // skip G0-regions
      if(lsSolver().m_G0regionId > -1 && a_potentialGapCellClose(lsId) == lsSolver().m_G0regionId) continue;
      MInt regionId = a_potentialGapCellClose(lsId) - 2;
      ASSERT(a_potentialGapCellClose(lsId) > 0 && a_potentialGapCellClose(lsId) <= lsSolver().m_noEmbeddedBodies, "");
      if(regionId > -1 && regionId < lsSolver().m_noGapRegions && globalTimeStep > 0) {
        // don't set hadGapCells during initialisation, so that initGapClosure is still called
        // at the first timeStep!
        m_hadGapCells[regionId]++;
      }
    }
    if(mode == 1) {
      fvMbSolver().a_wasGapCell(fc) = fvMbSolver().a_isGapCell(fc);
    } else {
      fvMbSolver().a_wasGapCell(fc) = false;
    }
    fvMbSolver().a_isGapCell(fc) = false;
  }


  for(MInt fc = a_noFvGridCells(); fc < fvMbSolver().noInternalCells(); fc++) {
    if(mode == 1) {
      fvMbSolver().a_wasGapCell(fc) = fvMbSolver().a_isGapCell(fc);
    } else {
      fvMbSolver().a_wasGapCell(fc) = false;
    }
    fvMbSolver().a_isGapCell(fc) = false;
  }

  // 4) set a_isGapCell and hasGapCells
  if(lsSolver().isActive()) {
    for(MInt gCellId = 0; gCellId < lsSolver().a_noCells(); gCellId++) {
      if(a_nearGapG(gCellId)) {
        const MInt fvId = ls2fvId(gCellId);
        if(fvId < 0) {
          ASSERT(abs(abs(a_levelSetFunctionG(gCellId, 0)) - a_outsideGValue()) < 0.000001,
                 "ERROR, no fv-Cell found for relevant Gap-Cells!");
        }
        if(a_potentialGapCellClose(gCellId) > 0 && a_potentialGapCellClose(gCellId) <= lsSolver().m_noEmbeddedBodies) {
          if(fvMbSolver().c_isLeafCell(fvId)) {
            fvMbSolver().a_isGapCell(fvId) = true;
            const MInt regionId = a_potentialGapCellClose(gCellId) - 2;
            if(regionId > -1 && regionId < lsSolver().m_noGapRegions) {
              m_hasGapCells[regionId]++;
            }
          } else if(!fvMbSolver().c_isLeafCell(fvId) && !fvMbSolver().a_isHalo(fvId) && !lsSolver().a_isHalo(gCellId)) {
            ASSERT((m_allowLsInterpolation && m_noRfJumps == 1) || fvMbSolver().m_maxLevelChange, "");
            // if the fv-cell is not a leaf cell, find leaf-cell parent
            // and set property for all corresponding childs
            for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
              const MInt child = fvMbSolver().c_childId(fvId, childId);
              if(child < 0) continue;
              fvMbSolver().a_isGapCell(child) = true;
              const MInt regionId = a_potentialGapCellClose(gCellId) - 2;
              if(regionId > -1 && regionId < lsSolver().m_noGapRegions) {
                m_hasGapCells[regionId]++;
              }
            }
          }
        }
      }
    }
  }

  testGapProperty();

  // set isGapCell for g0regionId:
  if(lsSolver().m_G0regionId > -1) {
    MInt noG0regionCells = 0;
    for(MInt gCellId = 0; gCellId < lsSolver().a_noCells(); gCellId++) {
      if(lsSolver().a_potentialGapCellClose(gCellId) == lsSolver().m_G0regionId && lsSolver().a_gapWidth(gCellId) > 0) {
        MInt fvId = ls2fvId(gCellId);
        if(fvId < 0) {
          ASSERT(abs(abs(a_levelSetFunctionG(gCellId, 0)) - a_outsideGValue()) < 0.000001,
                 "ERROR, no fv-Cell found for relevant Gap-Cells!");
        }
        noG0regionCells++;
        fvMbSolver().a_isGapCell(fvId) = true;
        ASSERT(!lsSolver().a_nearGapG(gCellId), "");
      }
    }
#if defined COUPLING_DEBUG_ || !defined NDEBUG
    MPI_Allreduce(MPI_IN_PLACE, &noG0regionCells, 1, MPI_INT, MPI_SUM, lsSolver().mpiComm(), AT_, "INPLACE",
                  "noG0regionCells");
    if(lsSolver().domainId() == 0) {
      cerr << " No of g0-region Cells " << noG0regionCells << endl;
    }
#endif
  }


  // at the beginning of the restart
  if(lsSolver().m_restart && globalTimeStep == fvMbSolver().m_restartTimeStep) {
    for(MInt fc = a_noFvGridCells(); fc < a_noFvCells(); fc++) {
      fvMbSolver().a_wasGapCell(fc) = fvMbSolver().a_isGapCell(fc);
    }
    for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
      m_hadGapCells[region] = m_hasGapCells[region];
    }
  }

  // 3) exchange hadGapCells and hasGapCells for all regions
  MPI_Allreduce(MPI_IN_PLACE, &m_hadGapCells[0], lsSolver().m_noGapRegions, MPI_INT, MPI_SUM, fvMbSolver().mpiComm(),
                AT_, "INPLACE", "m_hadGapCells");

  MPI_Allreduce(MPI_IN_PLACE, &m_hasGapCells[0], lsSolver().m_noGapRegions, MPI_INT, MPI_SUM, fvMbSolver().mpiComm(),
                AT_, "INPLACE", "m_hasGapCells");

  fvMbSolver().exchangeGapInfo();

#if defined COUPLING_DEBUG_ || !defined NDEBUG
  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    if(lsSolver().domainId() == 0) {
      cerr << globalTimeStep << " region " << region << " has " << m_hasGapCells[region] << " had "
           << m_hadGapCells[region] << endl;
    }
  }
#endif
}


/**
 * \brief returns a specific property of the specifuec body
 *        used to provide a unique function for both level-set and moving boundary code
 *        return mode:
 *        1: body venter
 *        2: body velocity
 *        3: body acceleration
 *        4: body temperature
 *
 *
 * \author Claudia Guenther, Tim Wegmann
 * \date 03/2011
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::computeBodyProperties(MInt returnMode, MFloat* bodyData, MInt body, MFloat time) {
  TRACE();

  // TODO labels:COUPLER,LS,toremove remove this pass through function!
  lsSolver().computeBodyPropertiesForced(returnMode, &bodyData[0], body, time);
}

/** \brief Updates the fv-mb-solver flow solver
 *         (after a completed levelSet TimeStep and finalizeLevelSet() )
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::updateFlowSolver() {
  TRACE();

  if(fvMbSolver().m_constructGField) return;

  fvMbSolver().m_structureStep = 0;

  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    fvMbSolver().m_gapState[region] = 0;
  }

  MBool& firstRun = m_static_updateLevelSetFlowSolver_firstRun;

  if(firstRun
     || (fvMbSolver().m_trackMovingBndry && globalTimeStep >= fvMbSolver().m_trackMbStart
         && globalTimeStep < fvMbSolver().m_trackMbEnd)) {
    testCoupling();

    // 1) transfer Gap-Cell property
    if(fvMbSolver().m_levelSetRans) {
      transferLevelSetValues();
    }
    transferGapCellProperty(1);

    // 2) transfer the levelSet Data to the flow solver
    transferLevelSetFieldValues(true);

    // 3) check gap cell-status
    setGapState();

    // 4) mark gap cells with status and region
    initFvGapCells();

    // 5) transfer body position/velocity/acceleration
    transferBodyProperties();
  }

  firstRun = false;
}

/** \brief Sets the Gap-State in the fvmb-solver (stati are the same for all domains!)
 * possible Gap-States:
 *
 * 0: default value   : should be overwritten!
 *-2: init Gap-Closure: the Gap switches from open state in the previous timeStep to a closed state
 * 2: init Gap-Opening: the Gap switches from closed state in the previous timeStep to an open state
 * 1: fully open Gap  : the Gap is wide open and no special Gap-Handling is necessary
 *-1: shrinking Gap   : the Gap is closed and is closing further (body-distance is decreasing)
 * 3: widening Gap    : the Gap is closed but the body-distance is increasing (moving appart)
 *
 *
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::setGapState() {
  TRACE();

  if(!lsSolver().m_closeGaps || lsSolver().m_gapInitMethod == 0) return;

  MBool(&earlyOpened)[m_maxNoGapRegions] = m_static_setGapState_earlyOpened;
  MBool& first = m_static_setGapState_first;

  if(first) {
    first = false;
    for(MInt region = 0; region < m_maxNoGapRegions; region++) {
      earlyOpened[region] = false;
    }
  }

  MBoolScratchSpace forceNoGaps(lsSolver().m_noGapRegions, AT_, "forceNoGaps");
  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    forceNoGaps[region] = false;
    if(lsSolver().m_forceNoGaps > 0) {
      const MFloat cad = lsSolver().crankAngle(lsSolver().m_time, 0);
      if(lsSolver().m_forceNoGaps == 1
         && ((lsSolver().m_gapSign[region] > F0 && cad > lsSolver().m_gapAngleOpen[region]
              && cad < lsSolver().m_gapAngleClose[region])
             || (lsSolver().m_gapSign[region] < F0
                 && (cad > lsSolver().m_gapAngleOpen[region] || cad < lsSolver().m_gapAngleClose[region])))) {
        forceNoGaps[region] = true;
      } else if(lsSolver().m_forceNoGaps == 2
                && ((cad > lsSolver().m_gapAngleOpen[region] && cad < lsSolver().m_gapAngleClose[region])
                    || (cad > lsSolver().m_gapAngleOpen[region + 1] && cad < lsSolver().m_gapAngleClose[region + 1]))) {
        forceNoGaps[region] = true;
      }
    }
  }

  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    if(forceNoGaps[region]) {
      if(!earlyOpened[region] && globalTimeStep > 0 && globalTimeStep != lsSolver().m_restartTimeStep) {
        fvMbSolver().m_gapState[region] = 2;
        earlyOpened[region] = true;
        if(fvMbSolver().domainId() == 0) {
          cerr << "--> Ls Gap Forced Opening Initialized for region " << region << " at timestep " << globalTimeStep
               << endl;
        }

      } else if(earlyOpened[region] || globalTimeStep == 0 || globalTimeStep == lsSolver().m_restartTimeStep) {
        fvMbSolver().m_gapState[region] = 1;
        earlyOpened[region] = true;
        if(fvMbSolver().domainId() == 0) {
          cerr << "--> Ls Gap Forced open for region " << region << " at timestep " << globalTimeStep << endl;
        }
      }
    }
    if(!forceNoGaps[region] && earlyOpened[region]) {
      lsSolver().m_minGapWidthDt1[region] = 11 * lsSolver().m_outsideGValue;
      earlyOpened[region] = false;
    }
  }

  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    if(forceNoGaps[region]) continue;

    if(earlyOpened[region]) {
      MFloat eps2 =
          lsSolver().c_cellLengthAtLevel(lsSolver().a_maxGCellLevel(0)) * sqrt(nDim) * lsSolver().m_gapDeltaMin;
      MFloat deviation = (fabs(lsSolver().m_minGapWidth[region] - eps2) / eps2) * 100;
      if(fvMbSolver().domainId() == 0) {
        cerr << "--> Ls Gap is Open early for region " << region << " at timestep " << globalTimeStep
             << " with deviation " << deviation << endl;
      }
      fvMbSolver().m_gapState[region] = 1;
      ASSERT(!m_hasGapCells[region] && !m_hadGapCells[region], "");
      if(lsSolver().m_minGapWidthDt1[region] < lsSolver().m_minGapWidth[region]
         && lsSolver().m_minGapWidth[region] > (10 * lsSolver().m_outsideGValue)) {
        earlyOpened[region] = false;
      }
      continue;
    }

    if(lsSolver().m_minGapWidthDt1[region] > (10 * lsSolver().m_outsideGValue)
       && lsSolver().m_minGapWidth[region] < lsSolver().m_minGapWidthDt1[region]) {
      if(fvMbSolver().domainId() == 0) {
        cerr << "--> Ls Gap Closure Initialized for region " << region << " at timestep " << globalTimeStep << " "
             << endl;
      }
      ASSERT(m_hasGapCells[region] && !m_hadGapCells[region],
             to_string(m_hasGapCells[region]) + " " + to_string(m_hadGapCells[region]));

      fvMbSolver().m_gapState[region] = -2;

    } else if(lsSolver().m_minGapWidthDt1[region] < lsSolver().m_minGapWidth[region]
              && lsSolver().m_minGapWidth[region] > (10 * lsSolver().m_outsideGValue)) {
      if(fvMbSolver().domainId() == 0) {
        cerr << "--> Ls Gap Opening Initialized for region " << region << " at timestep " << globalTimeStep << endl;
      }
      fvMbSolver().m_gapState[region] = 2;
      ASSERT(m_hadGapCells[region] && !m_hasGapCells[region], "");

    } else if(lsSolver().m_minGapWidth[region] < (10 * lsSolver().m_outsideGValue)
              && lsSolver().m_minGapWidthDt1[region] < (10 * lsSolver().m_outsideGValue)) {
      if(lsSolver().m_minGapWidthDt1[region] >= lsSolver().m_minGapWidth[region]) {
        // gap-shrinking at restart for >=
        if(fvMbSolver().domainId() == 0) {
          cerr << "--> Ls Gap is Shrinking for region " << region << " at timestep " << globalTimeStep << endl;
        }
        fvMbSolver().m_gapState[region] = -1;
        ASSERT(m_hasGapCells[region] && m_hadGapCells[region],
               to_string(m_hasGapCells[region]) + " " + to_string(m_hadGapCells[region]));

      } else {
        MFloat eps2 =
            lsSolver().c_cellLengthAtLevel(lsSolver().a_maxGCellLevel(0)) * sqrt(nDim) * lsSolver().m_gapDeltaMin;
        MFloat deviation = (fabs(lsSolver().m_minGapWidth[region] - eps2) / eps2) * 100;
        // for small deviations from the gap-Closing distance, check if any gap-Cells are remaining!
        // if non are remaining, initialize gap-Opening
        if(deviation < 1.0) {
          if(m_hasGapCells[region] == 0 && m_hadGapCells[region] > 0) {
            if(fvMbSolver().domainId() == 0) {
              cerr << "--> Ls Gap Opening Initialized early for region " << region << " at timestep " << globalTimeStep
                   << " with deviation " << deviation << endl;
            }
            fvMbSolver().m_gapState[region] = 2;
            ASSERT(m_hadGapCells[region] && !m_hasGapCells[region], "");
            earlyOpened[region] = true;

            continue;
          }
        }
        if(fvMbSolver().domainId() == 0) {
          cerr << "--> Ls Gap is Widening for region " << region << " at timestep " << globalTimeStep << endl;
        }
        fvMbSolver().m_gapState[region] = 3;
        ASSERT(lsSolver().m_minGapWidthDt1[region] < lsSolver().m_minGapWidth[region], "");
        ASSERT(m_hasGapCells[region] && m_hadGapCells[region],
               to_string(m_hasGapCells[region]) + " " + to_string(m_hadGapCells[region]));
      }

    } else if(lsSolver().m_minGapWidthDt1[region] > (10 * lsSolver().m_outsideGValue)
              && lsSolver().m_minGapWidth[region] > (10 * lsSolver().m_outsideGValue)) {
      if(fvMbSolver().domainId() == 0) {
        cerr << "--> Ls Gap is Open for region " << region << " at timestep " << globalTimeStep << endl;
      }
      fvMbSolver().m_gapState[region] = 1;
      ASSERT(!m_hasGapCells[region] && !m_hadGapCells[region], "");
    }
  }


  // possible trigger adaptation due to gap Opening/shrinking
  /*
  if(globalTimeStep != lsSolver().m_restartTimeStep && lsSolver().m_engineSetup) {
    for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
      const MInt state = fvMbSolver().m_gapState[region];
      if(fvMbSolver().m_gapState[region] == 2) { // opening always forces adaptation!
        if(fvMbSolver().domainId() == 0) {
          cerr << "Forcing Adaptation due to gap-state " << state << " of region " << region << endl;
        }
        lsSolver().m_forceAdaptation = true;
      } else if(fvMbSolver().m_gapState[region] == -1 && fvMbSolver().maxLevel() != fvMbSolver().maxRefinementLevel()
                && lsSolver().m_minGapWidth[region]
                       < (lsSolver().c_cellLengthAtLevel(lsSolver().a_maxGCellLevel(0)) * lsSolver().m_gapDeltaMin)) {
        if(fvMbSolver().domainId() == 0) {
          cerr << "Forcing Adaptation due to gap-state " << state << " of region " << region << endl;
        }
        lsSolver().m_forceAdaptation = true;
      }
    }
  }
  */
}


/** \brief Initialises fv-Gap Cells for gapHandling
 *         (must be called each timeStep before the gapHandling call in the fvmb-solver and
 *          after the gapHandling call in the lsSolver!)
 *
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::initFvGapCells() {
  fvMbSolver().m_gapCells.clear();
  // tuple of cellId, regionId, GapCell-Type (for isGapCell > 0, for wasGapCell < 0)

  if(!lsSolver().m_closeGaps) return;

  if(globalTimeStep < 0) return;

  if(lsSolver().m_gapInitMethod > 0) {
    fvMbSolver().m_noneGapRegions = true;
    for(MInt region = 0; region < fvMbSolver().m_noGapRegions; region++) {
      if(fvMbSolver().m_gapState[region] != 1) fvMbSolver().m_noneGapRegions = false;
      ASSERT(fvMbSolver().m_gapState[region] != 0 && fvMbSolver().m_gapState[region] >= -2
                 && fvMbSolver().m_gapState[region] <= 3,
             "Invalid GapState");
    }

    if(lsSolver().m_G0regionId > -1) fvMbSolver().m_noneGapRegions = false;
    if(fvMbSolver().m_noneGapRegions) return;
  }

  MIntScratchSpace lsInfo(fvMbSolver().c_noCells(), 3, AT_, "lsInfo");

  // a) list all Gap-Cells and previous Gap-Cells
  //   Gap-Cells have default Types:
  //   0: Gap-Cells that were not a Gap-Cell before
  //  -1: Gap-Cells that were a Gap-Cell before
  //   1: Cells that are not a Gap-Cell anymore but used to be a GapCell
  for(MInt cellId = 0; cellId < fvMbSolver().noInternalCells(); cellId++) {
    if(fvMbSolver().a_hasProperty(cellId, FvCell::IsSplitChild)) continue;
    fvMbSolver().assertValidGridCellId(cellId);
    if(fvMbSolver().a_isGapCell(cellId)) {
      lsGapInfo(cellId, &lsInfo(cellId, 0));
      const MInt status = fvMbSolver().a_wasGapCell(cellId) ? -1 : 0;
      fvMbSolver().m_gapCells.emplace_back(cellId, lsInfo(cellId, 0), status, lsInfo(cellId, 1), lsInfo(cellId, 2));
    } else if(fvMbSolver().a_wasGapCell(cellId)) {
      lsGapInfo(cellId, &lsInfo(cellId, 0));
      fvMbSolver().m_gapCells.emplace_back(cellId, lsInfo(cellId, 0), 1, lsInfo(cellId, 1), lsInfo(cellId, 2));
    }
  }

  fvMbSolver().exchangeData(&lsInfo[0], 3);

  for(MInt cellId = fvMbSolver().noInternalCells(); cellId < fvMbSolver().c_noCells(); cellId++) {
    if(fvMbSolver().a_hasProperty(cellId, FvCell::IsSplitChild)) continue;
    if(fvMbSolver().a_isBndryGhostCell(cellId)) continue;

    if(fvMbSolver().a_isGapCell(cellId)) {
      const MInt status = fvMbSolver().a_wasGapCell(cellId) ? -1 : 0;
      fvMbSolver().m_gapCells.emplace_back(cellId, lsInfo(cellId, 0), status, lsInfo(cellId, 1), lsInfo(cellId, 2));
    } else if(fvMbSolver().a_wasGapCell(cellId)) {
      fvMbSolver().m_gapCells.emplace_back(cellId, lsInfo(cellId, 0), 1, lsInfo(cellId, 1), lsInfo(cellId, 2));
    }
  }

  fvMbSolver().m_initGapCell = true;
}

/** \brief preCoupler
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::preCouple(MInt recepiStep) {
  TRACE();

  // PRE LS
  if(recepiStep == 0) {
    // check that the timeStep matches!
    if(fvMbSolver().isActive() && lsSolver().isActive()) {
      ASSERT(fabs(lsSolver().m_time - (fvMbSolver().m_physicalTime + lsSolver().m_timeStep)) < fvMbSolver().m_eps,
             "Different times in lSolver and fv-Solver! Expected times: " + to_string(lsSolver().m_time) + " "
                 + to_string(fvMbSolver().m_physicalTime) + " "
                 + to_string(fvMbSolver().m_physicalTime + lsSolver().m_timeStep));
    }
  }

  if(fvMbSolver().grid().wasAdapted()) {
    if(fvMbSolver().m_levelSetRans) {
      transferLevelSetValues();
    }
  }
}

/** \brief postCoupler
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::postCouple(MInt recepiStep) {
  TRACE();

  // POST LS
  if(recepiStep == 0) {
    updateFlowSolver();
  } else { // POST FVMB

    // update time and timeStep in the levelSet-solver
    // the timeStep is then increased in preTimeStep by the individual solver!
    transferTimeStep();

    fvMbSolver().m_initGapCell = false;
  }
}


/** \brief called after each solver-balance
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::finalizeBalance(const MInt solverId) {
  TRACE();

  // POST LS
  if(solverId == 0) {
    // 1) is/was Gap cell property is not exchanged in the fv-solver balance
    if(fvMbSolver().m_levelSetRans) {
      transferLevelSetValues();
    }
    transferGapCellProperty(0);

    // 2) ls-values are not exchanged in the fv-solver balance
    transferLevelSetFieldValues(true);

    // 3) gapState information is available on all ranks

    // 4) set fv-gap cells again
    initFvGapCells();
  }
}

/** \brief finalizeAdaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::finalizeAdaptation(const MInt solverId) {
  TRACE();

  if(solverId > 0) return; // call only once after the ls-solver!

  if(globalTimeStep > -1) { // after the lsSolver!


    // 1) transfer levelset values again after
    //   they have been set completely in the ls-solver!
    //   they are not interpolated in the fv-solver interpolation!
    transferLevelSetFieldValues(true);
    if(fvMbSolver().m_levelSetRans) {
      transferLevelSetValues();
    }

    // 2) is/was Gap cell property is handeled in the fv-solver during the adaptation
    //   => transferGapCellProperty is not necessary!

    // 3) gapState information does not change during adaptation!

    // 4) set fv-gap cells again, as they are not swapped in the fv-solver-adaptation!
    initFvGapCells();

    // 5) correct the time if the timeStep changes during adaptation!
    if(fvMbSolver().maxLevel() != fvMbSolver().m_maxLevelBeforeAdaptation) {
      ASSERT(fvMbSolver().isActive(), to_string(fvMbSolver().m_solverId) + " " + to_string(fvMbSolver().domainId()));
      const MFloat oldTimeStep = fvMbSolver().timeStep(true);
      if(fvMbSolver().domainId() == 0) {
        cerr << "Max-Level change from " << fvMbSolver().m_maxLevelBeforeAdaptation << " to " << fvMbSolver().maxLevel()
             << " during adaptation => timeStep ";
      }

      fvMbSolver().setTimeStep();
      // NOTE: the sweptVolume is updated according to the timeStep change in setTimeStep!


      const MFloat newTimeStep = fvMbSolver().timeStep(true);
      if(fvMbSolver().domainId() == 0) {
        cerr << "change from " << oldTimeStep << " to " << newTimeStep << endl;
      }

      // revert time and transfer the new timeStep!
      lsSolver().m_time = lsSolver().m_time - oldTimeStep;

      // reverse time in the fvmb-solver!
      fvMbSolver().m_time = fvMbSolver().m_time - oldTimeStep * fvMbSolver().m_timeRef;
      fvMbSolver().m_physicalTime = fvMbSolver().m_physicalTime - oldTimeStep;
      fvMbSolver().m_physicalTimeDt1 = fvMbSolver().m_physicalTimeDt1 - oldTimeStep;

      transferTimeStep();
      lsSolver().m_time = lsSolver().m_time + newTimeStep;

      fvMbSolver().advanceTimeStep();

      // compute the levelSet movement to the new TimeStep
      lsSolver().solutionStep();
      lsSolver().postTimeStep();

      // update levelsetfield!
      transferLevelSetFieldValues(true);
      if(fvMbSolver().m_levelSetRans) {
        transferLevelSetValues();
      }
    }


  } else if(globalTimeStep < 0 && fvMbSolver().m_geometryChange != nullptr) {
    // caution: lsSolver().m_geometryChange has already been reset to default!

    transferLevelSetFieldValues(true);
    if(fvMbSolver().m_levelSetRans) {
      transferLevelSetValues();
    }

    // init BndryLayer again to set wasInactive again
    fvMbSolver().initBndryLayer();

    // transfer oldG0 cells to the fv-block to apply possible damping in this area!
    /*
    if (!lsSolver().m_oldG0Cells.empty()) {
      for(auto it = lsSm_oldG0Cells.begin(); it != m_oldG0Cells.end(); it++) {
      const MInt cellId = it->first;
      const MInt set = it->second;

    }
    */
  }
}

/** \brief finalizeAdaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::postAdaptation() {
  TRACE();

  // if(globalTimeStep < 0 || m_allowLsInterpolation) {
  // handels the transfer during the initial adaptation and
  // also the interpolation of the levelset if the fv-solver is refined further!
  // TODO labels:COUPLER,FVMB first simple interpolation approximation in the fvmb-solver refineCell function!
  //                    then this call will be unncessary!
  //                    the ls-value will be updated in finalizeAdaptation, this is just for the sensors!
  transferLevelSetFieldValues(true);
  if(fvMbSolver().m_levelSetRans) {
    transferLevelSetValues();
  }
  //}

  updateLevelSet();
}

/** \brief finalizeAdaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::prepareAdaptation() {
  TRACE();

  updateLevelSet();
}

/** \brief update the levelset outside of the band
 *         default should be case 2 -> don't do anything!
 *  \author old function from Claudia, without warenty! Please test before use!
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::updateLevelSet() {
  ASSERT(!fvMbSolver().m_constructGField, "ERROR: invalid function call!");

  switch(fvMbSolver().m_levelSetAdaptationScheme) {
    case 0:
      updateLevelSetOutsideBandPar();
      break;
    case 2:
      return; // don't update OutsideCells at all!
      break;
    case 1:
    default:
      fvMbSolver().updateLevelSetOutsideBand();
      break;
  }
}


/** \brief computes an approximate level set value for cells outside the level-set computing band
 *         should not be used, as it is very slow,
 *         just for performance measurements of the parallel fast marching method!
 *    \author Claudia Guenther
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::updateLevelSetOutsideBandPar() {
  MFloatScratchSpace points(lsSolver().a_noG0Cells(0), nDim, AT_, "points");
  MFloatScratchSpace pointsLVS(lsSolver().a_noG0Cells(0), AT_, "pointsLVS");

  //--------------------
  MInt tempCnt = 0;
  for(MInt id = 0; id < lsSolver().a_noG0Cells(0); id++) {
    MInt cellId = lsSolver().a_G0CellId(id, 0);
    pointsLVS[tempCnt] = a_levelSetFunctionG(cellId, 0);
    for(MInt d = 0; d < nDim; d++) {
      points(tempCnt, d) = a_coordinateG(cellId, d);
    }
    tempCnt++;
  }

  // exchange points with all other processors:
  MIntScratchSpace noPoints(1, AT_, "noPoints");
  MIntScratchSpace globalNoPoints(fvMbSolver().noDomains(), AT_, "globalNoPoints");
  noPoints[0] = a_noG0Cells(0);
  MPI_Allgather(noPoints.getPointer(), 1, MPI_INT, globalNoPoints.getPointer(), 1, MPI_INT, fvMbSolver().mpiComm(), AT_,
                "noPoints.getPointer()", "globalNoPoints.getPointer()");

  MIntScratchSpace displs(fvMbSolver().noDomains(), AT_, "displs");
  MIntScratchSpace recCounts(fvMbSolver().noDomains(), AT_, "recCounts");
  MInt dataBlockSize = nDim + 1;
  noPoints[0] = 0;
  for(MInt i = 0; i < fvMbSolver().noDomains(); i++) {
    displs[i] = noPoints[0] * dataBlockSize;
    recCounts[i] = globalNoPoints[i] * dataBlockSize;
    noPoints[0] += globalNoPoints[i];
  }

  // gather all points from all other processors:
  MFloatScratchSpace globalPoints(noPoints[0], nDim, AT_, "globalPoints");
  MFloatScratchSpace globalPointsLVS(noPoints[0], AT_, "globalPointsLVS");
  MFloatScratchSpace tmpsend(globalNoPoints[fvMbSolver().domainId()] * dataBlockSize, AT_, "tmpsend");
  MFloatScratchSpace tmpreceive(noPoints[0] * dataBlockSize, AT_, "tmpreceive");
  for(MInt p = 0; p < globalNoPoints[fvMbSolver().domainId()]; p++) {
    for(MInt d = 0; d < nDim; d++)
      tmpsend[p * dataBlockSize + d] = points(p, d);
    tmpsend[p * dataBlockSize + dataBlockSize - 1] = pointsLVS[p];
  }

  MPI_Allgatherv(tmpsend.getPointer(), globalNoPoints[fvMbSolver().domainId()] * dataBlockSize, MPI_DOUBLE,
                 tmpreceive.getPointer(), recCounts.getPointer(), displs.getPointer(), MPI_DOUBLE,
                 fvMbSolver().mpiComm(), AT_, "tmpsend.getPointer()", "tmpreceive.getPointer()");

  for(MInt p = 0; p < noPoints[0]; p++) {
    for(MInt d = 0; d < nDim; d++)
      globalPoints(p, d) = tmpreceive[p * dataBlockSize + d];
    globalPointsLVS[p] = tmpreceive[p * dataBlockSize + dataBlockSize - 1];
  }

  for(MInt cell = 0; cell < fvMbSolver().a_noCells(); cell++) {
    const MInt gCellId = fv2lsId(cell);
    if(gCellId > -1 && a_inBandG(gCellId, 0)) continue;
    if(fvMbSolver().a_levelSetValuesMb(cell, 0) > F0) {
      fvMbSolver().a_levelSetValuesMb(cell, 0) = fvMbSolver().c_cellLengthAtLevel(0);
    } else {
      fvMbSolver().a_levelSetValuesMb(cell, 0) = -fvMbSolver().c_cellLengthAtLevel(0);
    }
    for(MInt p = 0; p < noPoints[0]; p++) {
      MFloat distTmp = F0;
      MFloat lvsValTmp = globalPointsLVS[p];
      for(MInt d = 0; d < nDim; d++) {
        distTmp += (fvMbSolver().a_coordinate(cell, d) - globalPoints(p, d))
                   * (fvMbSolver().a_coordinate(cell, d) - globalPoints(p, d));
      }
      distTmp = sqrt(distTmp);
      if(fvMbSolver().a_levelSetValuesMb(cell, 0) >= F0) {
        lvsValTmp += distTmp;
      } else {
        lvsValTmp -= distTmp;
      }
      if(abs(lvsValTmp) < abs(fvMbSolver().a_levelSetValuesMb(cell, 0))) {
        fvMbSolver().a_levelSetValuesMb(cell, 0) = lvsValTmp;
      }
    }
  }
}


/** \brief return restart time for lsSolver, when the fvSolver is not initialised yet!
 *  \author Tim Wegmann
 */

template <MInt nDim, class SysEqn>
MFloat LsFvMb<nDim, SysEqn>::restartTime() {
  TRACE();

  MInt timeStep = 0;
  MFloat time = 0;
  MFloat physicalTime = -std::numeric_limits<MFloat>::max();

  if(fvMbSolver().isActive() && fvMbSolver().m_restart) {
    stringstream varFileName;

    varFileName << fvMbSolver().restartDir() << "restartVariables";
    if(!fvMbSolver().m_useNonSpecifiedRestartFile) {
      if(!fvMbSolver().m_multipleFvSolver) {
        varFileName << "_" << fvMbSolver().m_restartTimeStep;
      } else {
        varFileName << fvMbSolver().solverId() << "_" << fvMbSolver().m_restartTimeStep;
      }
    }
    varFileName << ParallelIo::fileExt();

    if(fvMbSolver().domainId() == 0) {
      fvMbSolver().loadRestartTime((varFileName.str()).c_str(), timeStep, time, physicalTime);
    } else {
      physicalTime = 0;
    }
  } else if(fvMbSolver().isActive()) {
    physicalTime = 0;
  }

  if(fvMbSolver().isActive()) {
    MPI_Allreduce(MPI_IN_PLACE, &physicalTime, 1, MPI_DOUBLE, MPI_MAX, fvMbSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "physicalTime");
  }

  if(fvMbSolver().grid().hasInactiveRanks()) {
    if(lsSolver().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &physicalTime, 1, MPI_DOUBLE, MPI_MAX, lsSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "physicalTime");
    }
  }

  ASSERT(physicalTime > -fvMbSolver().m_eps, "");
  return physicalTime;
}


/** \brief returns the gap region for a given fvCell
 *  \author Tim Wegmann
 */

template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::lsGapInfo(const MInt fvCellId, MInt* info) {
  TRACE();

  MInt gCellId = fv2lsId(fvCellId);

  if(gCellId < 0) {
    ASSERT(m_allowLsInterpolation, "");
    gCellId = fv2lsIdParent(fvCellId);
  }
  ASSERT(gCellId > -1, "");
  MInt regionId = a_potentialGapCellClose(gCellId) - 2;
  ASSERT(regionId < fvMbSolver().m_noGapRegions, "");
  // labels:COUPLER G0region-hack
  if(lsSolver().m_G0regionId > -1 && a_potentialGapCellClose(gCellId) == lsSolver().m_G0regionId) {
    regionId = fvMbSolver().m_noGapRegions;
  }
  if(lsSolver().m_gapInitMethod > 0) {
    ASSERT(regionId > -1, "");
  }

  info[0] = regionId;
  info[1] = lsSolver().a_bodyIdG(gCellId, 0);
  info[2] = lsSolver().a_secondBodyId(gCellId);
}

/** \brief performs the coupling after solver initialization
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::init() {
  TRACE();

  // copies the bodyToSet information from the levelset solver to the fvmb-solver
  fvMbSolver().m_startSet = lsSolver().m_startSet;
  fvMbSolver().m_noSets = lsSolver().m_noSets;
  fvMbSolver().m_bodyToSetTable = lsSolver().m_bodyToSetTable;
  fvMbSolver().m_noBodiesInSet = lsSolver().m_noBodiesInSet;
  fvMbSolver().m_setToBodiesTable = lsSolver().m_setToBodiesTable;

  fvMbSolver().updateGeometry();

  // set the time in the lsSolver:
  lsSolver().m_time = -99;
  lsSolver().m_time = restartTime();

  if(lsSolver().isActive() && fvMbSolver().isActive()) {
    ASSERT(fabs(lsSolver().m_time - fvMbSolver().m_physicalTime) < fvMbSolver().m_eps, "");
  }

  checkProperties();

  // check for any in-active ranks
  MInt noInactiveFv = 0;
  MInt noInactiveLs = 0;
  if(!fvMbSolver().isActive()) {
    noInactiveFv++;
  }
  if(!lsSolver().isActive()) {
    noInactiveLs++;
  }

  MPI_Allreduce(MPI_IN_PLACE, &noInactiveFv, 1, MPI_INT, MPI_SUM, globalMaiaCommWorld(), AT_, "MPI_IN_PLACE", "noInactiveFv");
  MPI_Allreduce(MPI_IN_PLACE, &noInactiveLs, 1, MPI_INT, MPI_SUM, globalMaiaCommWorld(), AT_, "MPI_IN_PLACE", "noInactiveLs");
  if(noInactiveFv > 0) {
    cerr0 << "FvMb solver has " << noInactiveFv << " inactive ranks!" << endl;
  }
  if(noInactiveLs > 0) {
    cerr0 << "Ls solver has " << noInactiveLs << " inactive ranks!" << endl;
  }

  // initialise the g-Field in the fv-solver!
  transferLevelSetFieldValues(true);
  if(fvMbSolver().m_levelSetRans) {
    transferLevelSetValues();
  }

  // can only be called after the levelset values have been set!
  fvMbSolver().initBndryLayer();

  // TODO labels:COUPLER,FVMB,LS extend own shorter/simplified/easier/better coupler initBodyProperties version!
  //      don't call initBodyProperties() in the fvmb-solver at all for constructedGField with
  //      levelset motion function!
  if(fvMbSolver().m_LsMovement) {
    initBodyProperties();
  }
}

/** \brief performs the final coupling after finalization
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::finalizeCouplerInit() {
  TRACE();

  // ls-Solver time step is chosen corresponding to the fv-time step
  transferTimeStep();

  if(fvMbSolver().m_levelSetRans) {
    transferLevelSetValues();

    fvMbSolver().rhs();
    fvMbSolver().rhsBnd();
    fvMbSolver().advanceSolution();
    fvMbSolver().advanceBodies();
  }
}

/** \brief performs the coupling in solver finalization
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::finalizeSubCoupleInit(MInt currentSolver) {
  TRACE();

  if(lsSolver().isActive() && fvMbSolver().isActive()) {
    ASSERT(fabs(lsSolver().m_time - fvMbSolver().m_physicalTime) < fvMbSolver().m_eps, "");
  }

  if(currentSolver == m_lsSolverId) {
    updateFlowSolver();
    if(lsSolver().m_LsRotate) {
      transferBodyRadius();
    }
  }
}

/** \brief transfers the gcell time step from the fvMbSolver!
 *         time should already be matching!
 *
 * \author  Tim Wegmann
 * \date Jan 2020
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferTimeStep() {
  TRACE();

  // If fvMbSolver is inactive on one rank we need to take the time step from another rank!
  MFloat fvTimeStep = std::numeric_limits<MFloat>::max();

  // Get local fv time step
  if(fvMbSolver().isActive()) {
    fvTimeStep = fvMbSolver().timeStep(true);
  }

  // Get minimum fv time step on all active ls ranks
  if(fvMbSolver().grid().hasInactiveRanks()) {
    if(lsSolver().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &fvTimeStep, 1, MPI_DOUBLE, MPI_MIN, lsSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "fvTimeStep");
    }
  }

  if(lsSolver().isActive() && fvMbSolver().isActive()
     && fabs(lsSolver().m_time - fvMbSolver().m_physicalTime) > fvMbSolver().m_eps) {
    cerr << "Fv-Solver " << fvMbSolver().isActive() << " ls-Solver " << lsSolver().isActive() << endl;
    cerr << "Fv-Time " << setprecision(12) << fvMbSolver().m_physicalTime << " ls-Time " << lsSolver().m_time << endl;
    cerr << "Fv-TS " << setprecision(12) << fvTimeStep << " ls-TS " << lsSolver().m_timeStep << endl;
    mTerm(1, AT_, "Time in Solvers differs!.");
  }
  lsSolver().m_timeStep = fvTimeStep;
}

/** \brief Copies the bodyRadius from the LS-solver into the FV-solver
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferBodyRadius() {
  TRACE();

  const MInt noBodies = lsSolver().m_noEmbeddedBodies;
  for(MInt b = 0; b < noBodies; b++) {
    fvMbSolver().m_bodyRadius[b] = lsSolver().m_bodyRadius[b];
  }
}

/** \brief transfer current body properties from the ls-Solver to the fv-Solver
 *          bodyPosition
 *          bodyVelocity
 *          bodyAcceleration
 *          bodyTemperature
 *          bodyAngularVelocity
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferBodyProperties() {
  TRACE();

  if(!fvMbSolver().m_LsMovement) return;

  // the fv-time has not been increased yet, but will be first thing in preTimeStep!
  const MFloat nextFvTime = fvMbSolver().m_physicalTime + fvMbSolver().timeStep();

  for(MInt body = 0; body < fvMbSolver().m_noEmbeddedBodies; body++) {
    computeBodyProperties(1, &fvMbSolver().m_bodyCenter[body * nDim], body, nextFvTime);
    computeBodyProperties(2, &fvMbSolver().m_bodyVelocity[body * nDim], body, nextFvTime);
    computeBodyProperties(3, &fvMbSolver().m_bodyAcceleration[body * nDim], body, nextFvTime);

    if(fvMbSolver().m_movingBndryCndId == 3008 || fvMbSolver().m_movingBndryCndId == 3010) {
      computeBodyProperties(4, &fvMbSolver().m_bodyTemperature[body], body, nextFvTime);
    }
  }

  if(lsSolver().m_LsRotate) {
    transferAngularVelocity();
  }
}


/** \brief Writes the angular velocity and the angular acceleration from the
 *         LS-solver into the FV-solver
 *  \author Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::transferAngularVelocity() {
  TRACE();

  for(MInt b = 0; b < lsSolver().m_noEmbeddedBodies; b++) {
    for(MInt i = 0; i < nDim; i++) {
      if(lsSolver().isActive()) {
        fvMbSolver().m_bodyAngularVelocity[b * nDim + i] = lsSolver().m_bodyAngularVelocity[b * nDim + i];
        fvMbSolver().m_bodyAngularAcceleration[b * nDim + i] = lsSolver().m_bodyAngularAcceleration[b * nDim + i];
      } else {
        fvMbSolver().m_bodyAngularVelocity[b * nDim + i] = NAN;
        fvMbSolver().m_bodyAngularAcceleration[b * nDim + i] = NAN;
      }
    }
  }
}

/** \brief build the combined levelSet for the given cellId
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::buildCollectedLevelSet(const MInt cellId) {
  TRACE();

  if(!fvMbSolver().a_isGapCell(cellId)) {
    fvMbSolver().a_levelSetValuesMb(cellId, 0) = fvMbSolver().a_levelSetValuesMb(cellId, 1);
    fvMbSolver().a_associatedBodyIds(cellId, 0) = fvMbSolver().a_associatedBodyIds(cellId, 1);
    for(MInt set = 2; set < a_noLevelSetsMb(); set++) {
      MFloat phi0 = fvMbSolver().a_levelSetValuesMb(cellId, 0);
      MFloat phi1 = fvMbSolver().a_levelSetValuesMb(cellId, set);
      MInt body0 = fvMbSolver().a_associatedBodyIds(cellId, 0);
      MInt body1 = fvMbSolver().a_associatedBodyIds(cellId, set);

      //
      if(phi0 >= F0 && phi1 >= F0) {
        if(phi1 < phi0) {
          body0 = body1;
        }
        phi0 = mMin(phi0, phi1);
      } else if(phi0 <= F0 && phi1 <= F0) {
        if(abs(phi1) > abs(phi0)) {
          body0 = body1;
        }
        phi0 = -sqrt(phi0 * phi0 + phi1 * phi1);
      } else if(phi0 * phi1 <= F0) {
        if(phi0 < F0) {
          // phi0 = phi0;
        } else {
          phi0 = phi1;
          body0 = body1;
        }
      }

      fvMbSolver().a_levelSetValuesMb(cellId, 0) = phi0;
      fvMbSolver().a_associatedBodyIds(cellId, 0) = body0;
    }
  }
}

/** \brief initialise the body properties in the fvmb-solver based
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void LsFvMb<nDim, SysEqn>::initBodyProperties() {
  TRACE();

  const MFloat time = fvMbSolver().m_physicalTime;

  for(MInt body = 0; body < fvMbSolver().m_noEmbeddedBodies; body++) {
    computeBodyProperties(1, &fvMbSolver().m_bodyCenter[body * nDim], body, time);
    computeBodyProperties(2, &fvMbSolver().m_bodyVelocity[body * nDim], body, time);
    computeBodyProperties(3, &fvMbSolver().m_bodyAcceleration[body * nDim], body, time);

    if(fvMbSolver().m_movingBndryCndId == 3008 || fvMbSolver().m_movingBndryCndId == 3010) {
      computeBodyProperties(4, &fvMbSolver().m_bodyTemperature[body], body, time);
    }
  }
}

template class LsFvMb<2, FvSysEqnNS<2>>;
template class LsFvMb<3, FvSysEqnNS<3>>;
template class LsFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class LsFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class LsFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class LsFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class LsFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class LsFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
