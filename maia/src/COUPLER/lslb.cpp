// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

//#define COUPLING_DEBUG_
#include "lslb.h"
#include <algorithm>
#include <stack>
#include <vector>
#include "GRID/cartesiannetcdf.h"
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "globals.h"

#include "globalvariables.h"

using namespace std;

template <MInt nDim, MInt nDist, class SysEqn>
LsLb<nDim, nDist, SysEqn>::LsLb(MInt couplingId, LsSolver* ls, LbSolver* lb)
  : Coupling(couplingId), CouplingLS<nDim>(couplingId, ls), CouplingLB<nDim, nDist, SysEqn>(couplingId, lb) {
  TRACE();

  initData();
  readProperties();
  checkProperties();
  updateGeometry();
}

/** \brief Initialize coupling-class-specific Data
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::initData() {
  TRACE();

  // solver-specific data:
  m_lbSolverId = lbSolver().solverId();
  m_lsSolverId = lsSolver().solverId();

  /*! \page propertyPage1
    \section transferRegion
    <code>MFloat LsLb::m_transferRegion</code>\n
    default = <code></code>\n\n
    Set the region in which the LS values are transferred to the LB solver.\n
    The array is filled by the following convention: -x, -y, -z, x, y, z\n
    Keywords: <i>LATTICE_BOLTZMANN, LEVEL-SET</i>
  */
  if(Context::propertyLength("transferRegion", m_lbSolverId) > 0) {
    if(Context::propertyLength("transferRegion", m_lbSolverId) == 2 * nDim) {
      mAlloc(m_transferBoundingBox, 2 * nDim, "m_transferBoundingBox", AT_);
      for(MInt d = 0; d < 2 * nDim; d++) {
        m_transferBoundingBox[d] = Context::getSolverProperty<MFloat>("transferRegion", m_lbSolverId, AT_, d);
      }
    } else {
      stringstream errorMessage;
      errorMessage << "ERROR: Wrong number of array entries in transferRegion. Must be 2 * nDim!!";
      mTerm(1, AT_, errorMessage.str());
    }
  }
}

/** \brief Checks property-data which is read in by both ls-and Lb-Solver
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::checkProperties() {
  TRACE();

  lbSolver().m_noEmbeddedBodies = lsSolver().m_noEmbeddedBodies;

  if(lsSolver().m_maxNoSets > 0) lbSolver().m_noLevelSetsUsedForMb = lsSolver().m_maxNoSets;
  if(lbSolver().m_useOnlyCollectedLS) lbSolver().m_noLevelSetsUsedForMb = 1;

  lbSolver().m_maxNoSets = lbSolver().m_noLevelSetsUsedForMb;

  lbSolver().m_constructGField = lsSolver().m_constructGField;
}

/** \brief reads lslb-coupling-specific data
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::readProperties() {
  TRACE();

  /*! \page propertyPage1
  \section outsideDefault
  <code>MBool* LsLb::outsideDefault </code>\n
  default = <code>true</code>\n \n
  A trigger which determines the default levelset values of a lb-cell if no connection a ls-cell
  can be found! This is the case if the ls-domain is smaller than the lb-domain!
  Possible values are:
  <ul>
  <li>true: possitive outsideGValue </li>
  <li>false: negative outsideGValue </li>
  </ul>
  Keywords: <i>LEVELSET, MULTIPLE LEVEL SET FUNCTIONS</i>
  */

  m_outsideDefault = true;
  if(Context::propertyExists("outsideDefault", m_lbSolverId)) {
    m_outsideDefault = Context::getSolverProperty<MBool>("outsideDefault", m_lbSolverId, AT_, &m_outsideDefault);
  }
}

/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::testCoupling() {
  TRACE();

  if(lsSolver().m_constructGField) return;

  for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
    ASSERT(lbSolver().grid().tree().solver2grid(lbCellId) >= 0, "");
    ASSERT(lbSolver().grid().solverFlag(lbSolver().grid().tree().solver2grid(lbCellId), m_lbSolverId), "");
#ifdef COUPLING_DEBUG_
    if(lsSolver().grid().solverFlag(lbSolver().grid().tree().solver2grid(lbCellId), m_lsSolverId)) {
      ASSERT(ls2lbId(lb2lsId(lbCellId)) == lbCellId,
             to_string(lbCellId) + " " + to_string(lb2lsId(lbCellId)) + " " + to_string(ls2lbId(lb2lsId(lbCellId))));
    }
#endif
  }

  for(MInt cellId = 0; cellId < a_noLsCells(); cellId++) {
    ASSERT(lsSolver().grid().tree().solver2grid(cellId) >= 0, "");
    ASSERT(lsSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_lsSolverId), "");
#ifdef COUPLING_DEBUG_
    if(lbSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_lbSolverId)) {
      ASSERT(lb2lsId(ls2lbId(cellId)) == cellId,
             to_string(cellId) + " " + to_string(ls2lbId(cellId)) + " " + to_string(lb2lsId(ls2lbId(cellId))));
    }
#endif
  }

  if(!lsSolver().m_levelSetLb) return;

#ifdef COUPLING_DEBUG_
  for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
    for(MInt dir = 0; dir < nDim; dir++) {
      if(lsSolver().grid().solverFlag(lbSolver().grid().tree().solver2grid(lbCellId), m_lsSolverId)) {
        ASSERT(abs(lbSolver().c_coordinate(lbCellId, dir) - a_coordinateG(lb2lsId(lbCellId), dir)) < 0.00000001, "");
      }
    }
  }
#endif
}

/** \brief transfers the LevelSetValues for all cells from the levelset to the moving boundary Part
 *    \author Christoph Siewert, Tim Wegmann
 */

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::transferLevelSetFieldValues(MBool exchangeLVS) {
  TRACE();

  const MInt noRfJumps = mMax(0, lbSolver().maxLevel() - lsSolver().maxRefinementLevel());

  // Nothing to transfer if lbSolver is inactive!
  if(!lbSolver().isActive()) return;

  if(exchangeLVS) {
    if(lsSolver().isActive()) {
      lsSolver().exchangeAllLevelSetData();
      lsSolver().checkHaloCells();
    }
  }

  if(lsSolver().isActive()) {
    std::list<MInt> interpolationParents;
    interpolationParents.clear();

    for(MInt cellId = 0; cellId < a_noLbCells(); cellId++) {
      // reset values with invalid/defaults
      for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
        if(m_outsideDefault) {
          lbSolver().a_levelSetFunctionMB(cellId, set) = a_outsideGValue();
        } else {
          lbSolver().a_levelSetFunctionMB(cellId, set) = -a_outsideGValue();
        }
        lbSolver().a_associatedBodyIds(cellId, set) = -2;
      }
      if(m_transferBoundingBox != nullptr) {
        if(m_transferBoundingBox[0] > lbSolver().c_coordinate(cellId, 0)
           || m_transferBoundingBox[0 + nDim] < lbSolver().c_coordinate(cellId, 0)) {
          continue;
        }
        if(m_transferBoundingBox[1] > lbSolver().c_coordinate(cellId, 1)
           || m_transferBoundingBox[1 + nDim] < lbSolver().c_coordinate(cellId, 1)) {
          continue;
        }
        IF_CONSTEXPR(nDim == 3) {
          if(m_transferBoundingBox[2] > lbSolver().c_coordinate(cellId, 2)
             || m_transferBoundingBox[2 + nDim] < lbSolver().c_coordinate(cellId, 2)) {
            continue;
          }
        }
      }
      const MInt gCellId = lb2lsId(cellId);
      // direct transfer of matching cells
      if(gCellId > -1) {
        for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
          lbSolver().a_levelSetFunctionMB(cellId, set) = a_levelSetFunctionG(gCellId, set);
          lbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellId, set);
        }
      } else {
        if(!g_multiSolverGrid) {
          ASSERT(cellId >= a_noLbCells(), "");
        } else {
          // multiSolverGrid, where the level and the grid-extension might differ!
          if(cellId <= a_noLbCells()) { // grid-cell!
            const MInt gCellParent = lb2lsIdParent(cellId);

            if(gCellParent > -1) { // just a different level

              if(noRfJumps == 0) {
                // same maxRefinementLevel, meaning just a difference due to different refinement-widths
                for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                  lbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                  // ASSERT(!lsSolver().a_inBandG( gCellParent ,  set ), "");
                  lbSolver().a_levelSetFunctionMB(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                }
              } else {
                // ls-interpolation towards higher lb-level
                if(lbSolver().a_isHalo(cellId) || lsSolver().a_isHalo(gCellParent)) continue;
                // simple transfer if not a band-cell
                if(lsSolver().m_buildCollectedLevelSetFunction && !lsSolver().a_inBandG(gCellParent, 0)) {
                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    lbSolver().a_levelSetFunctionMB(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                    lbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                  }
                  continue;
                }
                const MInt levelDifference = lbSolver().c_level(cellId) - lsSolver().a_level(gCellParent);
                if(levelDifference == 1) {
                  // direct interpolation from ls to fv for only 1 level difference
                  MFloat point[3] = {lbSolver().a_coordinate(cellId, 0), lbSolver().a_coordinate(cellId, 1),
                                     lbSolver().a_coordinate(cellId, 2)};
                  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
                  MInt position = lsSolver().setUpLevelSetInterpolationStencil(gCellParent, interpolationCells, point);

                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    lbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                    if(position < 0) {
                      // no valid interpolation stencil found
                      lbSolver().a_levelSetFunctionMB(cellId, set) = a_levelSetFunctionG(gCellParent, set);
                    } else { // interpolation for all sets
                      const MFloat phi = lsSolver().interpolateLevelSet(interpolationCells, point, set);

                      lbSolver().a_levelSetFunctionMB(cellId, set) = phi;
                    }
                  }
                  // build collected levelset data
                  buildCollectedLevelSet(cellId);

                } else {
                  // add cells based on which the fv-interpolation needs to start
                  // that is cells for which the level difference is above 1
                  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
                    lbSolver().a_associatedBodyIds(cellId, set) = a_bodyIdG(gCellParent, set);
                  }
                  MInt parent = cellId;
                  for(MInt i = 0; i < levelDifference - 1; i++) {
                    parent = lbSolver().c_parentId(parent);
                  }
                  ASSERT(parent > -1, "");
                  // parent of this parent has a direct link, meaning that the values for the parent
                  // are set above!
                  ASSERT(lb2lsId(lbSolver().c_parentId(parent)) == gCellParent, "");
                  interpolationParents.push_back(parent);
                }
              }
            }
          }
        }
      }
    }

    // now interpolate through all childs
    // iterate downwards level by level and exchange values in between!
    if(!interpolationParents.empty()) {
      cerr << "This part is not tested yet!!!" << endl;
      std::list<MInt> newInterpolationParents;

      for(MInt lvl = lsSolver().maxUniformRefinementLevel(); lvl < lbSolver().maxRefinementLevel(); lvl++) {
        lbSolver().exchangeData(&(lbSolver().a_levelSetFunctionMB(0, 0)), lbSolver().m_maxNoSets);
        lbSolver().exchangeData(&(lbSolver().a_associatedBodyIds(0, 0)), lbSolver().m_maxNoSets);

        MInt noInterpolationCells = (MInt)interpolationParents.size();
        MPI_Allreduce(MPI_IN_PLACE, &noInterpolationCells, 1, MPI_INT, MPI_MAX, lbSolver().mpiComm(), AT_, "INPLACE",
                      "noInterpolationCells");
        if(noInterpolationCells == 0) break;

        if(interpolationParents.empty()) continue;

        ASSERT(noRfJumps > 1, "");

        interpolationParents.sort();
        interpolationParents.unique();
        newInterpolationParents.clear();

        for(auto it = interpolationParents.begin(); it != interpolationParents.end(); it++) {
          const MInt parent = (*it);
          for(MInt c = 0; c < lbSolver().grid().m_maxNoChilds; c++) {
            const MInt childId = lbSolver().c_childId(parent, c);
            if(childId < 0) continue;
            // interpolateLsLb(parent, childId);
            if(!lbSolver().c_isLeafCell(childId)) {
              newInterpolationParents.push_back(childId);
            }
          }
        }
        interpolationParents.clear();
        for(auto it = newInterpolationParents.begin(); it != newInterpolationParents.end(); it++) {
          const MInt cellId = (*it);
          interpolationParents.push_back(cellId);
        }
      }
    }
  } else {
    for(MInt cellId = 0; cellId < a_noLbCells(); cellId++) {
      for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
        if(m_outsideDefault) {
          lbSolver().a_levelSetFunctionMB(cellId, set) = a_outsideGValue();
        } else {
          lbSolver().a_levelSetFunctionMB(cellId, set) = -a_outsideGValue();
        }
        lbSolver().a_associatedBodyIds(cellId, set) = -2;
      }
    }
  }
}

/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *
 */

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::createBodyTree() {
  TRACE();
  IF_CONSTEXPR(nDim == 3) { updateGeometry(); }
}

//---------------------------------------------------------------------------

template <MInt nDim, MInt nDist, class SysEqn>
MInt LsLb<nDim, nDist, SysEqn>::noLevelSetFieldData() {
  MInt noLevelSetFieldData = 0;

  if(lsSolver().m_writeOutAllLevelSetFunctions) {
    noLevelSetFieldData += noLevelSetFieldData + 2 * a_noSets();
  } else {
    noLevelSetFieldData += 1;
  }
  if(lsSolver().m_writeOutAllCurvatures) {
    noLevelSetFieldData += a_noSets();
  } else {
    noLevelSetFieldData += 1;
  }

  return noLevelSetFieldData;
}


/** \brief Updates the member-variables in the geometry-intersection class
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::updateGeometry() {
  TRACE();

  lbSolver().m_geometryIntersection->m_noEmbeddedBodies = lbSolver().m_noEmbeddedBodies;
  lbSolver().m_geometryIntersection->m_noLevelSetsUsedForMb = lbSolver().m_noLevelSetsUsedForMb;
  lbSolver().m_geometryIntersection->m_bodyToSetTable = lsSolver().m_bodyToSetTable;
  lbSolver().m_geometryIntersection->m_setToBodiesTable = lsSolver().m_setToBodiesTable;
  lbSolver().m_geometryIntersection->m_noBodiesInSet = lsSolver().m_noBodiesInSet;
}

/** \brief Updates the fv-mb-solver flow solver
 *         (after a completed levelSet TimeStep and finalizeLevelSet() )
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::updateLevelSetFlowSolver() {
  TRACE();

  if(lbSolver().m_constructGField) return;

  MBool& firstRun = m_static_updateLevelSetFlowSolver_firstRun;

  if(firstRun
     || (lbSolver().m_trackMovingBndry && globalTimeStep >= lbSolver().m_trackMbStart
         && globalTimeStep < lbSolver().m_trackMbEnd)) {
    testCoupling();

    // 1) transfer the levelSet Data to the flow solver
    transferLevelSetFieldValues(true);
  }

  firstRun = false;
}

//-----

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::init() {
  TRACE();

  if(lbSolver().isActive()) {
    lbSolver().initializeMovingBoundaries();
    lbBndCnd().initializeBndMovingBoundaries();
  }
  updateLevelSetFlowSolver();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::finalizeCouplerInit() {
  TRACE();
  lsSolver().m_timeStep = F1;
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::finalizeSubCoupleInit(MInt couplingStep) {
  TRACE();
  std::ignore = couplingStep;
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::preCouple(MInt step) {
  TRACE();

  switch(step) {
    case 0: {
      break;
    }
    case 1: {
      updateLevelSetFlowSolver();
      std::vector<MInt> maxGCellLevels(lsSolver().m_maxNoSets);
      for(MInt set = 0; set < lsSolver().m_maxNoSets; set++) {
        maxGCellLevels[set] = lsSolver().a_maxGCellLevel(set);
      }

      lbSolver().preCoupleLs(maxGCellLevels);

      lbSolver().createBndryToBodyMapping(bndryToBodyMapping, bodyToBndryMapping);

      lbBndCnd().createMBComm();

      initializeSolidDomain();

      MFloatScratchSpace bodyVelocities(m_maxNoEmbeddedBodies, nDim, AT_, "bodyVelocities");
      for(MInt body = 0; body < m_maxNoEmbeddedBodies; body++) {
        lsSolver().computeBodyPropertiesForced(2, &bodyVelocities(body, 0), body, globalTimeStep);
      }

      for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
        for(auto body : bndryToBodyMapping[mbCell]) {
          for(MInt n = 0; n < nDim; n++) {
            a_mbCell().velocity(mbCell, n) = bodyVelocities(body, n);
          }
        }
      }

      break;
    }
    default: {
      break;
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::postCouple(MInt step) {
  TRACE();

  switch(step) {
    case 0: {
      break;
    }
    case 1: {
      lbBndCnd().postCouple();
      lsSolver().m_timeStep = F1;
      break;
    }
    default: {
      break;
    }
  }
}

/** \brief finalizeAdaptation
 *  \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::postAdaptation() {
  TRACE();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::finalizeAdaptation(const MInt solverId) {
  TRACE();

  if(solverId == 1) {
    if(lbSolver().isActive()) {
      lbSolver().initializeMovingBoundaries();
      lbBndCnd().initializeBndMovingBoundaries();
    }
    updateLevelSetFlowSolver();
  }
}

// Attempt to move initSolidDomain to coupler ...
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::initializeSolidDomain() {
  TRACE();

  if(!lbSolver().isActive()) {
    return;
  }

  MFloatScratchSpace bodyVelocities(a_noEmbeddedBodies(), nDim, AT_, "bodyVelocities");

  // Get the body Velocity for each embedded body
  for(MInt body = 0; body < a_noEmbeddedBodies(); body++) {
    lsSolver().computeBodyPropertiesForced(2, &bodyVelocities(body, 0), body, globalTimeStep);
  }

  for(MInt i = 0; i < lbSolver().a_noCells(); i++) {
    // Regular fluid cell
    if(a_isActive(i) && a_wasActive(i)) {
      continue;
    }

    // Regular Solid cell
    if(!a_isActive(i)) {
      // TODO labels:COUPLER,toremove Skip non-leaf cells ...
      /*if(!lbSolver().c_isLeafCell(i)){
        continue;
      }*/

      // determine the LS-Body to which the cell belongs, to set the right body velocity
      MInt bodyId = -1;
      MInt setOfBody = 0;
      for(MInt set = lbSolver().m_levelSetId; set < lbSolver().m_maxNoSets; set++) {
        if(a_associatedBodyIdsMb(i, set) >= 0) {
          bodyId = a_associatedBodyIdsMb(i, set);
          setOfBody = set;
          break;
        }
      }

      // the Velocity of the deactivated cell is set to the body velocity
      // the Density is set to 1.0
      if((bodyId >= 0) && (bodyId < a_noEmbeddedBodies()) && (a_levelSetFunctionMb(i, setOfBody) < 0)) {
        for(MInt j = 0; j < nDim; j++) {
          a_variable(i, j) = bodyVelocities(bodyId, j);
          a_oldVariable(i, j) = bodyVelocities(bodyId, j);
        }
        MFloat squaredVelocity = 0.0;
        for(MInt n = 0; n < nDim; n++) {
          squaredVelocity += bodyVelocities(bodyId, n) * bodyVelocities(bodyId, n);
        }
        lbSolver().setEqDists(i, 1.0, &bodyVelocities(bodyId, 0));
        if(a_isThermal()) {
          a_variable(i, a_pvt()) = a_initTemperatureKelvin();
          lbSolver().setEqDistsThermal(i, a_initTemperatureKelvin(), F1, squaredVelocity, &bodyVelocities(bodyId, 0));
        }
      } else {
        for(MInt j = 0; j < a_noVariables(); j++) {
          a_variable(i, j) = F0;
          a_oldVariable(i, j) = F0;
        }
        MFloat squaredVelocity = 0.0;
        MFloat bodyVelocitiesOutsideLs[nDim] = {F0};
        lbSolver().setEqDists(i, 1.0, bodyVelocitiesOutsideLs);
        if(a_isThermal()) {
          a_variable(i, a_pvt()) = a_initTemperatureKelvin();
          lbSolver().setEqDistsThermal(i, a_initTemperatureKelvin(), F1, squaredVelocity, bodyVelocitiesOutsideLs);
        }
      }
      a_variable(i, a_pvrho()) = 1.0;
      a_oldVariable(i, a_pvrho()) = 1.0;
    } else {
      lbBndCnd().refillEmergedCell(i);
    }
  }
}

/** \brief build the combined levelSet for the given cellId
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::buildCollectedLevelSet(const MInt cellId) {
  TRACE();

  lbSolver().a_levelSetFunctionMB(cellId, 0) = lbSolver().a_levelSetFunctionMB(cellId, 1);
  lbSolver().a_associatedBodyIds(cellId, 0) = lbSolver().a_associatedBodyIds(cellId, 1);
  for(MInt set = 2; set < a_noLevelSetsMb(); set++) {
    MFloat phi0 = lbSolver().a_levelSetFunctionMB(cellId, 0);
    MFloat phi1 = lbSolver().a_levelSetFunctionMB(cellId, set);
    MInt body0 = lbSolver().a_associatedBodyIds(cellId, 0);
    MInt body1 = lbSolver().a_associatedBodyIds(cellId, set);

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

    lbSolver().a_levelSetFunctionMB(cellId, 0) = phi0;
    lbSolver().a_associatedBodyIds(cellId, 0) = body0;
  }
}

/** \brief interpolate levelset values on the lb-grid
 *
 * \author Tim Wegmann
 * \date March 2020
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLb<nDim, nDist, SysEqn>::interpolateLsLb(const MInt from, const MInt to) {
  TRACE();

  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MFloat point[3] = {lbSolver().a_coordinate(to, 0), lbSolver().a_coordinate(to, 1), lbSolver().a_coordinate(to, 2)};

  const MInt position = lbSolver().setUpLbInterpolationStencil(from, interpolationCells, point);

  for(MInt set = 0; set < a_noLevelSetsMb(); set++) {
    if(position < 0) {
      lbSolver().a_levelSetFunctionMB(to, set) = lbSolver().a_levelSetFunctionMB(from, set);
    } else {
      MFloat phi = interpolateLevelSet(interpolationCells, point, set);
      lbSolver().a_levelSetFunctionMB(to, set) = phi;
    }
  }
  buildCollectedLevelSet(to);
}

///  interpolates the levelSet value further on fv-cells
/// \author Tim Wegmann
/// \date 2020-04-01
template <MInt nDim, MInt nDist, class SysEqn>
MFloat LsLb<nDim, nDist, SysEqn>::interpolateLevelSet(MInt* interpolationCells, MFloat* point, const MInt set) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt refSet) {
    return static_cast<MFloat>(lbSolver().a_levelSetFunctionMB(cellId, refSet));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(lbSolver().a_coordinate(cellId, id));
  };

  return lbSolver().interpolateFieldDataLb(&interpolationCells[0], &point[0], set, scalarField, coordinate);
}


template class LsLb<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class LsLb<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LsLb<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
