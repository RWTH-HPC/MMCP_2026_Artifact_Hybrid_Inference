// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lsfv.h"

#include <algorithm>
#include <stack>
#include <vector>
#include "FV/fvcartesiansolverxd.h"
#include "LS/lscartesiansolver.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "UTIL/kdtree.h"
#include "globals.h"

#include "globalvariables.h"

#include "coupling.h"

using namespace std;

template <MInt nDim, class SysEqn>
CouplingLsFv<nDim, SysEqn>::CouplingLsFv(const MInt couplingId, LsSolver* ls, FvCartesianSolver* fv)
  : Coupling(couplingId), CouplingLS<nDim>(couplingId, ls), CouplingFv<nDim, SysEqn>(couplingId, fv) {
  TRACE();

  initData();
  readProperties();
  checkProperties();

  for(MInt dir = 0; dir < nDim; dir++) {
    lsSolver().a_meanCoord(dir) = a_meanCoord(dir);
  }
}

/** \initializes the coupler in context of the unified run loop
 *    \author Thomas Lürkens
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::init() {
  mAlloc(fvSolver().m_levelSetValues, lsSolver().m_noSets, "fvSolver().m_levelSetValues", AT_);
  mAlloc(fvSolver().m_curvatureG, lsSolver().m_noSets, "fvSolver().m_curvatureG", AT_);

  // set the time in the lsSolver:
  lsSolver().m_time = -99;
  lsSolver().m_time = a_time();
  // transferGapCellProperty();
  // computeGCellTimeStep();
  return;
}

/** \initializes the coupler in context of the unified run loop
 *    \author Thomas Lürkens
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::finalizeCouplerInit() {
  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    fvSolver().m_levelSetValues[s].resize(fvSolver().a_noCells());
    fvSolver().m_curvatureG[s].resize(fvSolver().a_noCells());
  }

  transferGapCellProperty();
  transferLevelSetValues();
  computeGCellTimeStep();
}

/** \ preCouple routine in context of the unified run loop
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::preCouple(MInt /*recipeStep*/) {
  if(!lsSolver().m_combustion && !lsSolver().m_levelSetMb) {
    if(lsSolver().m_semiLagrange) {
      returnStep_semiLagrange();
    }
  }
}

/** \ postCouple routine in context of the unified run loop
 *    \author Thomas Lürkens
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::postCouple(MInt /*recipeStep*/) {
  lsSolver().m_time = a_time();

  transferLevelSetValues();
  transferGapCellProperty();
  testCoupling();

  if(!lsSolver().m_combustion && !lsSolver().m_levelSetMb && !lsSolver().m_semiLagrange && !lsSolver().m_LSSolver) {
    returnStep();
  }
}

/** \brief finalizeAdaptation
 *  \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::postAdaptation() {
  TRACE();

  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    fvSolver().m_levelSetValues[s].clear();
    fvSolver().m_levelSetValues[s].resize(fvSolver().a_noCells());
    fvSolver().m_curvatureG[s].clear();
    fvSolver().m_curvatureG[s].resize(fvSolver().a_noCells());
  }

  transferLevelSetValues();
}


/** \brief Initialize coupling-class-specific Data
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::initData() {
  TRACE();

  // solver-specific data:
  m_fvSolverId = fvSolver().m_solverId;
  m_lsSolverId = lsSolver().m_solverId;

  if(lsSolver().m_closeGaps) {
    mAlloc(m_hadGapCells, lsSolver().m_noGapRegions, "m_hadGapCells", 0, AT_);
    mAlloc(m_hasGapCells, lsSolver().m_noGapRegions, "m_hasGapCells", 0, AT_);
  }
}

/** \brief Checks property-data which is read in by both ls-and Fv-Solver
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::checkProperties() {
  TRACE();

  ASSERT(lsSolver().a_maxGCellLevel() == fvSolver().maxRefinementLevel(), "");
}

/** \brief Checks Propertty-Data which is read in by both ls-and Fv-Solver
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::readProperties() {
  TRACE();

  // claudia: ugly hardcoded patch for fast checkup of level-set boundary motion - sorry...
  /*! \page propertyPage1
    \section maxLevelsetVelocityForTimestep
    <code>MFloat LsCartesianSolver::computeGCellTimeStep()::maxVelocity</code>\n
    default = <code>14.5</code>\n \n
    Set maximum velocity used to compute m_timeStep if m_timeStepMethod == 7. \n \n
    Possible values are:
    <ul>
      <li>any positive floating point value</li>
    </ul>
    Keywords: <i>LEVELSET, TIMESTEP</i>
  */
  m_maxVelocity = 14.5;
  if(Context::propertyExists("maxLevelsetVelocityForTimestep", m_lsSolverId)) {
    m_maxVelocity =
        Context::getSolverProperty<MFloat>("maxLevelsetVelocityForTimestep", m_lsSolverId, AT_, &m_maxVelocity);
  }

  m_cfl = -1;
  if(Context::propertyExists("cfl", m_lsSolverId)) {
    m_cfl = Context::getSolverProperty<MFloat>("cfl", m_lsSolverId, AT_);
  }

  m_G0regionId = -1;
  if(Context::propertyExists("G0regionId", m_fvSolverId)) {
    m_G0regionId = Context::getSolverProperty<MInt>("G0regionId", m_fvSolverId, AT_, &m_G0regionId);
  }

  m_initialCrankAngle = F0;
  if(Context::propertyExists("initialCrankAngle", m_fvSolverId)) {
    m_initialCrankAngle =
        Context::getSolverProperty<MFloat>("initialCrankAngle", m_fvSolverId, AT_, &m_initialCrankAngle);
  }

  m_timeStepMethod = 8;
  if(Context::propertyExists("timeStepMethod", m_lsSolverId)) {
    m_timeStepMethod = Context::getSolverProperty<MInt>("timeStepMethod", m_lsSolverId, AT_);
  }

  m_solverMethod = Context::getSolverProperty<MString>("solverMethod", m_lsSolverId, AT_);

  /*! \page propertyPage1
    \section bandWidthRef bandWidthRefMax
    <code>MFloat CouplingLsFv::bandWidthRef</code>\n
    default = <code> 4/5 </code>\n \n
    Defines the number of cells (on all levels and on the maxRefinementLevel) which will be refined around the G0-cells
    <ul>
    <li>any positive floating point value</li>
    </ul>
    Keywords: <i> LEVELSET, Refinement </i>
  */
  m_bandWidthRef = 4;
  m_bandWidthRefMax = 5;
  if(Context::propertyExists("bandWidthRef", m_fvSolverId)) {
    m_bandWidthRef = Context::getSolverProperty<MInt>("bandWidthRef", m_fvSolverId, AT_);
    m_bandWidthRefMax = Context::getSolverProperty<MInt>("bandWidthRefMax", m_fvSolverId, AT_);
  }
}


/** \brief computes the gcell time step
 *    \author Daniel Hartmann, Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::computeGCellTimeStep() {
  TRACE();

  if(m_timeStepMethod < 6 || m_timeStepMethod > 8) {
    if(string2enum(m_solverMethod) == MAIA_SEMI_LAGRANGE_LEVELSET
       || string2enum(m_solverMethod) == MAIA_RUNGE_KUTTA_LEVELSET) {
      mTerm(1, AT_, "Computation with pure LVS not possible with timeStepMethod other than 6 ,7 or 8! Please check!");
    }
  }

  ASSERT(m_cfl > 0, "Couldn't read cfl property for the ls-solver timestepping!");

  if(m_timeStepMethod == 8) {
    lsSolver().m_timeStep = m_cfl * lsSolver().m_gCellDistance;

  } else if(m_timeStepMethod == 6) {
    lsSolver().m_timeStep = m_cfl * lsSolver().m_gCellDistance;
    fvSolver().forceTimeStep(lsSolver().m_timeStep / a_timeRef());

  } else if(m_timeStepMethod == 7) {
    lsSolver().m_timeStep = m_cfl * lsSolver().m_gCellDistance / m_maxVelocity;
    fvSolver().forceTimeStep(lsSolver().m_timeStep / a_timeRef());

  } else {
    lsSolver().m_timeStep = lsTimeStep();
  }
}

/** \brief mimics the behaviour of the rungeKuttaStep() methods with respect to increasing time
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
MBool CouplingLsFv<nDim, SysEqn>::returnStep() {
  TRACE();

  fvSolver().m_time += lsTimeStep();
  fvSolver().m_physicalTime += lsTimeStep() * a_timeRef();

  return true;
}

/** \brief mimics the behaviour of the rungeKuttaStep() methods with respect to increasing time
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::returnStep_semiLagrange() {
  TRACE();

  fvSolver().m_time += lsTimeStep() * a_timeRef();
  fvSolver().m_physicalTime += lsTimeStep() * a_timeRef();
}


template <MInt nDim, class SysEqn>
MFloat CouplingLsFv<nDim, SysEqn>::interpolateLevelSet(MInt cellId, MFloat* point, MInt set) {
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


template <MInt nDim, class SysEqn>
MInt CouplingLsFv<nDim, SysEqn>::noLevelSetFieldData() {
  MInt noLevelSetFieldData = 0;
  ASSERT(fvSolver().m_levelSet, "");

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


template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::setLsInList(MIntScratchSpace&) {
  TRACE();

  cerr0 << "Setting Sensors for Fv-Solver adaptation!" << endl;

  MInt bandWidthRef = m_bandWidthRef;
  MInt bandWidthRefmax = m_bandWidthRefMax;

  MIntScratchSpace sendBufferSize(fvSolver().grid().noNeighborDomains(), AT_, "sendBufferSize");
  MIntScratchSpace receiveBufferSize(fvSolver().grid().noNeighborDomains(), AT_, "receiveBufferSize");
  MInt listCount = 0;
  MIntScratchSpace inList(a_noFvGridCells(), AT_, "inList");
  for(MInt c = 0; c < a_noFvGridCells(); c++)
    inList[c] = 0;

  // 2) add the g0-Cells and all parents to the list:
  MInt endSet = a_noSets();
  if(lsSolver().m_buildCollectedLevelSetFunction) {
    endSet = 1;
  }

  for(MInt set = 0; set < endSet; set++) {
    for(MInt id = 0; id < lsSolver().a_noG0Cells(set); id++) {
      MInt fvCellId = ls2fvIdParent(a_G0CellId(id, set));
      if(fvCellId < 0) continue;
      if(fvSolver().a_isHalo(fvCellId)) continue;
      inList[fvCellId] = 1;
      listCount++;
      MInt parentId = fvSolver().c_parentId(fvCellId);
      while(parentId > -1) {
        if(parentId < a_noFvGridCells()) inList[parentId] = 1;
        parentId = fvSolver().c_parentId(parentId);
      }
    }
  }

  // Exchange the listCount on all Domains
  MPI_Allreduce(&listCount, &listCount, 1, MPI_INT, MPI_SUM, fvSolver().mpiComm(), AT_, "listCount", "listCount");

  if(listCount == 0) {
    if(fvSolver().domainId() == 0) cerr << "No G0-Cells found!" << endl;
  }

  fvSolver().exchangeData(&inList[0], 1);

  for(MInt level = fvSolver().minLevel(); level < fvSolver().maxRefinementLevel(); level++) {
    if(level == fvSolver().maxRefinementLevel() - 1) bandWidthRef = bandWidthRefmax;

    // 4) Loop over the number of revinement-grid-cells and add those to the list:
    for(MInt loopMarker = 1; loopMarker < bandWidthRef; loopMarker++) {
      for(MInt cellId = 0; cellId < a_noFvGridCells(); cellId++) {
        if(fvSolver().a_level(cellId) != level) continue;
        if(inList[cellId] != loopMarker) continue;

        // direct neighbor
        for(MInt n = 0; n < fvSolver().m_noDirs; n++) {
          MInt nghbrId = fvSolver().c_neighborId(cellId, n, false);
          if(nghbrId < 0) continue;
          if(inList[nghbrId] == 0) inList[nghbrId] = loopMarker + 1;

          // diagonal neighbors
          if(n == 0 || n == 1) {
            for(MInt nn = 2; nn < fvSolver().m_noDirs; nn++) {
              MInt nghbrId2 = fvSolver().c_neighborId(nghbrId, nn, false);
              if(nghbrId2 < 0) continue;
              if(inList[nghbrId2] == 0) inList[nghbrId2] = loopMarker + 1;
            }
          }
          if(n == 4 || n == 5) {
            for(MInt nn = 2; nn < 4; nn++) {
              MInt nghbrId2 = fvSolver().c_neighborId(nghbrId, nn, false);
              if(nghbrId2 < 0) continue;
              if(inList[nghbrId2] == 0) inList[nghbrId2] = loopMarker + 1;

              for(MInt nnn = 0; nnn < 1; nnn++) {
                MInt nghbrId3 = fvSolver().c_neighborId(nghbrId2, nnn, false);
                if(nghbrId3 < 0) continue;
                if(inList[nghbrId3] == 0) inList[nghbrId3] = loopMarker + 1;
              }
            }
          }
        }
      }

      fvSolver().exchangeData(&inList[0], 1);
    }
  }
}


/** \brief Sets the Levelset-Values in fvSolver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::transferLevelSetValues() {
  TRACE();

  fvSolver().a_noSets() = a_noSets();
  fvSolver().a_noLevelSetFieldData() = noLevelSetFieldData();

  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    fvSolver().m_levelSetValues[s].resize(fvSolver().a_noCells());
  }

  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    for(MInt fc = 0; fc < a_noFvCells(); fc++) {
      MInt lsId = fv2lsIdParent(fc);
      fvSolver().a_levelSetFunction(fc, s) = a_outsideGValue();
      if(lsId > -1) {
        fvSolver().a_curvatureG(fc, s) = lsSolver().a_curvatureG(lsId, 0);
        if(fvSolver().a_isInterface(fc)) {
          // Interplate levelset for cut-cells
          MFloat point[nDim];
          for(MInt d = 0; d < nDim; d++) {
            point[d] = fvSolver().a_coordinate(fc, d);
          }
          fvSolver().a_levelSetFunction(fc, s) = interpolateLevelSet(fv2lsIdParent(fc), point, s);
        } else {
          fvSolver().a_levelSetFunction(fc, s) = lsSolver().a_levelSetFunctionG(fv2lsIdParent(fc), s);
        }
      }
    }
  }
}

/** \brief Sets the gapcell-property
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::transferGapCellProperty() {
  TRACE();

  if(!lsSolver().m_closeGaps) return;
  if(lsSolver().m_noGapRegions <= 0) return;
  if((lsSolver().m_levelSetMb) && fvSolver().m_constructGField) return;

#ifdef COUPLING_DEBUG_
  testCoupling();
#endif

  // 1) reset has/had Gap-Cells
  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    m_hadGapCells[region] = 0;
    m_hasGapCells[region] = 0;
  }

  // 2) set a_wasGapCell, m_hadGapCells and initialize a_isGapCell
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    if(fvSolver().a_isGapCell(fc)) {
      MInt lsId = fv2lsId(fc);
      ASSERT(lsId > -1, "");
      // skip G0-regions
      if(m_G0regionId > -1 && a_potentialGapCellClose(lsId) == m_G0regionId) continue;
      MInt regionId = a_potentialGapCellClose(lsId) - 2;
      ASSERT(a_potentialGapCellClose(lsId) > 0 && a_potentialGapCellClose(lsId) <= lsSolver().m_noEmbeddedBodies, "");
      if(regionId > -1 && regionId < lsSolver().m_noGapRegions && globalTimeStep > 0) {
        // don't set hadGapCells during initialisation, so that initGapClosure is still called
        // at the first timeStep!
        m_hadGapCells[regionId]++;
      }
    }
    if(globalTimeStep > 0) {
      fvSolver().a_wasGapCell(fc) = fvSolver().a_isGapCell(fc);
    } else {
      fvSolver().a_wasGapCell(fc) = false;
    }
    fvSolver().a_isGapCell(fc) = false;
  }

  for(MInt fc = a_noFvGridCells(); fc < a_noFvCells(); fc++) {
    if(globalTimeStep > 0) {
      fvSolver().a_wasGapCell(fc) = fvSolver().a_isGapCell(fc);
    } else {
      fvSolver().a_wasGapCell(fc) = false;
    }
    fvSolver().a_isGapCell(fc) = false;
  }

  // 4) set a_isGapCell and hasGapCells
  for(MInt gCellId = 0; gCellId < lsSolver().a_noCells(); gCellId++) {
    if(a_nearGapG(gCellId)) {
      MInt fvId = ls2fvId(gCellId);
      if(fvId < 0) {
        ASSERT(abs(abs(a_levelSetFunctionG(gCellId, 0)) - a_outsideGValue()) < 0.000001,
               "ERROR, no fv-Cell found for relevant Gap-Cells!");
      }
      if(a_potentialGapCellClose(gCellId) > 0 && a_potentialGapCellClose(gCellId) <= lsSolver().m_noEmbeddedBodies) {
        fvSolver().a_isGapCell(fvId) = true;
        MInt regionId = a_potentialGapCellClose(gCellId) - 2;
        if(regionId > -1 && regionId < lsSolver().m_noGapRegions) {
          m_hasGapCells[regionId]++;
        }
      }
    }
  }
  testGapProperty();

  // set isGapCell for g0regionId:
  if(m_G0regionId > -1) {
    MInt noG0regionCells = 0;
    for(MInt gCellId = 0; gCellId < lsSolver().a_noCells(); gCellId++) {
      if(lsSolver().a_potentialGapCellClose(gCellId) == m_G0regionId && lsSolver().a_gapWidth(gCellId) > 0) {
        MInt fvId = ls2fvId(gCellId);
        if(fvId < 0) {
          ASSERT(abs(abs(a_levelSetFunctionG(gCellId, 0)) - a_outsideGValue()) < 0.000001,
                 "ERROR, no fv-Cell found for relevant Gap-Cells!");
        }
        noG0regionCells++;
        fvSolver().a_isGapCell(fvId) = true;
        ASSERT(!lsSolver().a_nearGapG(gCellId), "");
      }
    }
#if defined COUPLING_DEBUG_ || !defined NDEBUG
    MPI_Allreduce(MPI_IN_PLACE, &noG0regionCells, 1, MPI_INT, MPI_SUM, fvSolver().mpiComm(), AT_, "INPLACE",
                  "noG0regionCells");
    if(lsSolver().domainId() == 0) {
      cerr << " No of g0-region Cells " << noG0regionCells << endl;
    }
#endif
  }


  // at the beginning of the restart
  if(lsSolver().m_restart && globalTimeStep == fvSolver().m_restartTimeStep) {
    for(MInt fc = a_noFvGridCells(); fc < a_noFvCells(); fc++) {
      fvSolver().a_wasGapCell(fc) = fvSolver().a_isGapCell(fc);
    }
    for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
      m_hadGapCells[region] = m_hasGapCells[region];
    }
  }

  // 3) exchange hadGapCells and hasGapCells for all regions
  MPI_Allreduce(MPI_IN_PLACE, &m_hadGapCells[0], lsSolver().m_noGapRegions, MPI_INT, MPI_SUM, fvSolver().mpiComm(), AT_,
                "INPLACE", "m_hadGapCells");

  MPI_Allreduce(MPI_IN_PLACE, &m_hasGapCells[0], lsSolver().m_noGapRegions, MPI_INT, MPI_SUM, fvSolver().mpiComm(), AT_,
                "INPLACE", "m_hasGapCells");

  fvSolver().exchangeGapInfo();

#if defined COUPLING_DEBUG_ || !defined NDEBUG
  for(MInt region = 0; region < lsSolver().m_noGapRegions; region++) {
    if(lsSolver().domainId() == 0) {
      cerr << globalTimeStep << " region " << region << " has " << m_hasGapCells[region] << " had "
           << m_hadGapCells[region] << endl;
    }
  }
#endif
}


/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *    \author  Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::testGapProperty() {
  TRACE();

  if(!lsSolver().m_closeGaps) return;
  if((lsSolver().m_levelSetMb) && fvSolver().m_constructGField) return;

#ifdef COUPLING_DEBUG_
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    if(fvSolver().a_isGapCell(fc) || fvSolver().a_wasGapCell(fc)) {
      MInt gc = lsSolver().grid().tree().grid2solver(fvSolver().grid().tree().solver2grid(fc));
      ASSERT(gc > -1, "");
      ASSERT(lsSolver().a_potentialGapCellClose(gc) > 0
                 && lsSolver().a_potentialGapCellClose(gc) <= fvSolver().m_noEmbeddedBodies,
             "");
    }
  }
#endif
}

/** \brief transfers the LevelSetValues from the levelset to the moving boundary Part
 *    \author  Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::testCoupling() {
  TRACE();

  if(lsSolver().m_constructGField) return;

  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    ASSERT(fvSolver().grid().tree().solver2grid(fc) >= 0, "");
    ASSERT(fvSolver().grid().solverFlag(fvSolver().grid().tree().solver2grid(fc), m_fvSolverId), "");
#ifdef COUPLING_DEBUG_
    if(lsSolver().grid().solverFlag(fvSolver().grid().tree().solver2grid(fc), m_lsSolverId)) {
      ASSERT(ls2fvId(fv2lsId(fc)) == fc,
             to_string(fc) + " " + to_string(fv2lsId(fc)) + " " + to_string(ls2fvId(fv2lsId(fc))));
    }
#endif
  }
  for(MInt cellId = 0; cellId < a_noLsCells(); cellId++) {
    ASSERT(lsSolver().grid().tree().solver2grid(cellId) >= 0, "");
    ASSERT(lsSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_lsSolverId), "");
#ifdef COUPLING_DEBUG_
    if(fvSolver().grid().solverFlag(lsSolver().grid().tree().solver2grid(cellId), m_fvSolverId)) {
      ASSERT(fv2lsId(ls2fvId(cellId)) == cellId,
             to_string(cellId) + " " + to_string(ls2fvId(cellId)) + " " + to_string(fv2lsId(ls2fvId(cellId))));
    }
#endif
  }

  if(!lsSolver().m_levelSetMb) return;


#ifdef COUPLING_DEBUG_
  for(MInt fc = 0; fc < a_noFvGridCells(); fc++) {
    for(MInt dir = 0; dir < nDim; dir++) {
      if(lsSolver().grid().solverFlag(fvSolver().grid().tree().solver2grid(fc), m_lsSolverId)) {
        ASSERT(abs(fvSolver().c_coordinate(fc, dir) - a_coordinateG(fv2lsId(fc), dir)) < 0.00000001, "");
      }
    }
  }
#endif
}


/**
 * \brief returns the current position of the body center(s) in bodyCenter
 * used to provide a unique function for both level-set and moving boundary code
 * if invoked with moving boundary code (m_levelSetMb is true), the respective mb function is called!
 * \author Claudia Guenther
 * \date 03/2011
 */
template <MInt nDim, class SysEqn>
void CouplingLsFv<nDim, SysEqn>::computeBodyProperties(MInt returnMode, MFloat* bodyData, MInt body, MFloat time) {
  TRACE();

  MFloat elapsedTime = time;
  MFloat angle = F0;
  MBool& first = m_static_computeBodyProperties_first;
  MFloat(&amplitude)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_amplitude;
  MFloat(&freqFactor)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_freqFactor;
  MFloat(&initialBodyCenter)[m_maxNoEmbeddedBodies * 3] = m_static_computeBodyProperties_initialBodyCenter;
  MFloat& Strouhal = m_static_computeBodyProperties_Strouhal;
  MFloat(&mu)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_mu;
  MFloat(&mu2)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_mu2;
  MFloat(&liftStartAngle1)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftStartAngle1;
  MFloat(&liftEndAngle1)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftEndAngle1;
  MFloat(&liftStartAngle2)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftStartAngle2;
  MFloat(&liftEndAngle2)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_liftEndAngle2;
  MFloat(&circleStartAngle)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_circleStartAngle;
  MFloat(&normal)[m_maxNoEmbeddedBodies * 3] = m_static_computeBodyProperties_normal;
  MInt(&bodyToFunction)[m_maxNoEmbeddedBodies] = m_static_computeBodyProperties_bodyToFunction;
  MFloat& omega = m_static_computeBodyProperties_omega;
  MFloat& rotAngle = m_static_computeBodyProperties_rotAngle;

  if(first) {
    MInt noEmbeddedBodies = lsSolver().m_noEmbeddedBodies;

    if(noEmbeddedBodies > m_maxNoEmbeddedBodies) {
      mTerm(1, AT_, "Error in computeBodyProperties: too many embedded Bodies!");
    }

    // 1: set default values:
    Strouhal = 0.2;
    for(MInt k = 0; k < m_maxNoEmbeddedBodies; k++) {
      amplitude[k] = 0.1;
      freqFactor[k] = 1.0;
      bodyToFunction[k] = 1;
      for(MInt i = 0; i < nDim; i++) {
        initialBodyCenter[k * nDim + i] = F0;
        normal[k * nDim + i] = F0;
      }
      normal[k * nDim + 0] = 1.0;
      liftStartAngle1[k] = F0;
      liftEndAngle1[k] = PI;
      liftStartAngle2[k] = 3.0 * PI;
      liftEndAngle2[k] = 4.0 * PI;
      circleStartAngle[k] = F0;
    }

    // 2: read Properties

    /*! \page propertyPage1
      \section amplitudes
      <code>MFloat* LsCartesianSolver::m_static_computeBodyProperties_amplitude</code>\n
      default = <code>none</code>\n \n
      Amplitude for body motion of embedded bodies. \n
      NOTE: also used in FV-MB solver for some special cases. \n \n
      Possible values are:
      <ul>
        <li>list of floating point numbers</li>
      </ul>
      Keywords: <i>LEVELSET, MOVING, BODY, BODY_MOTION</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      amplitude[i] = Context::getSolverProperty<MFloat>("amplitudes", m_lsSolverId, AT_, &amplitude[i], i);

    /*! \page propertyPage1
      \section freqFactors
      <code>MFloat* LsCartesianSolver::m_static_computeBodyProperties_freqFactor</code>\n
      default = <code>none</code>\n \n
      Set the frequency factors for prescribing body motion for all embedded bodies. \n
      NOTE: also used in FV-MB solver for some special cases. \n \n
      Possible values are:
      <ul>
        <li>list of positive floating point numbers</li>
      </ul>
      Keywords: <i>LEVELSET, MOVING, BODY, BODY_MOTION</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      freqFactor[i] = Context::getSolverProperty<MFloat>("freqFactors", m_lsSolverId, AT_, &freqFactor[i], i);

    /*! \page propertyPage1
      \section bodyMovementFunctions
      <code>MInt* LsCartesianSolver::m_static_computeBodyProperties_bodyToFunction</code>\n
      default = <code>1</code>\n \n
      Prescribes the functions for the body movement. Check the switch case in
      computeBodyProperties() for what each case actually does. \n \n
      Possible values are:
      <ul>
        <li>integer from 1 to 7</li>
      </ul>
      Keywords: <i>LEVELSET, BODY, BODY_MOTION, MOVING</i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      bodyToFunction[i] =
          Context::getSolverProperty<MInt>("bodyMovementFunctions", m_lsSolverId, AT_, &bodyToFunction[i], i);

    Strouhal = Context::getSolverProperty<MFloat>("Strouhal", m_lsSolverId, AT_, &Strouhal);

    /*! \page propertyPage1
\section initialBodyCenters
<code>MFloat LsCartesianSolver::initialBodyCenter</code>\n
      For each body, nDim Float values.
      With this property, one can move an STL around to its initial position.
      In initialInsidePoints is the "real" center of the stl files.
      This property, is the movement from the "real" center to the initial center used during calculation.

Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      for(MInt j = 0; j < nDim; j++)
        initialBodyCenter[i * nDim + j] = Context::getSolverProperty<MFloat>(
            "initialBodyCenters", m_lsSolverId, AT_, &initialBodyCenter[i * nDim + j], i * nDim + j);

    for(MInt i = 0; i < noEmbeddedBodies; i++)
      for(MInt j = 0; j < nDim; j++)
        normal[i * nDim + j] = Context::getSolverProperty<MFloat>("bodyMotionNormals", m_lsSolverId, AT_,
                                                                  &normal[i * nDim + j], i * nDim + j);

    /*! \page propertyPage1
\section liftStartAngles1
<code>MFloat LsCartesianSolver::liftStartAngle1</code>\n
default = <code>0.0</code>\n \n
      For each body, sets the start angle of the translation.
      The translation is described by a:
        - bodyToFunction case 1: cosine function
        - bodyToFunction case 2: valve lift shifted quadratic sine function (first angle)
<ul>
<li>Between 0 and 4.0 * PI </li>
</ul>
Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      liftStartAngle1[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles1", m_lsSolverId, AT_, &liftStartAngle1[i], i);

    /*! \page propertyPage1
\section liftStartAngles2
<code>MFloat LsCartesianSolver::liftStartAngle2</code>\n
default = <code>3.0 * PI</code>\n \n
      For each body, sets the start angle of the translation.
      The translation is described by a:
        - bodyToFunction case 2: valve lift shifted quadratic sine function (second angle)
<ul>
<li>Between 0 and 4.0 * PI </li>
</ul>
Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      liftStartAngle2[i] =
          Context::getSolverProperty<MFloat>("liftStartAngles2", m_lsSolverId, AT_, &liftStartAngle2[i], i);

    /*! \page propertyPage1
    \section liftEndAngles1
    <code>MFloat LsCartesianSolver::liftEndAngle1</code>\n
    default = <code>PI</code>\n \n
          For each body, sets the end angle of the translation.
    <ul>
    <li>Between 0 and 4.0 * PI </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      liftEndAngle1[i] = Context::getSolverProperty<MFloat>("liftEndAngles1", m_lsSolverId, AT_, &liftEndAngle1[i], i);

    /*! \page propertyPage1
    \section liftEndAngles2
    <code>MFloat LsCartesianSolver::liftEndAngle2</code>\n
    default = <code>4.0 * PI</code>\n \n
          For each body, sets the end angle of the translation.
    <ul>
    <li>Between 0 and 4.0 * PI </li>
    </ul>
    Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      liftEndAngle2[i] = Context::getSolverProperty<MFloat>("liftEndAngles2", m_lsSolverId, AT_, &liftEndAngle2[i], i);

    /*! \page propertyPage1
\section circleStartAngles
<code>MFloat LsCartesianSolver::circleStartAngle</code>\n
default = <code>0.0</code>\n \n
      For each body, sets the start angle for the circular motion case.
      (i.e. bodyToFunction case 5).
<ul>
<li>Between 0 and 4.0 * PI </li>
</ul>
Keywords: <i>LEVELSET, MULTILEVELSET, MB </i>
    */
    for(MInt i = 0; i < noEmbeddedBodies; i++)
      circleStartAngle[i] =
          Context::getSolverProperty<MFloat>("circleStartAngles", m_lsSolverId, AT_, &circleStartAngle[i], i);

    /*! \page propertyPage1
    \section rotAngle
    <code> MFloat LsPar::rotAngle </code>  \n
    default = <code>"0.0"</code>\n \n
    Used to rotate the primary direction of the body movement.
    information in the property-file is expected in degree!
    <ul>
    <li>Any positive Float</li>
    </ul>
    Keywords: <i>LEVELSET EMBEDED BOUNDARY, MOVEMENT FUNCTIONS </i>
   */

    rotAngle = 0.0;
    rotAngle = Context::getSolverProperty<MFloat>("rotAngle", m_lsSolverId, AT_, &rotAngle);
    rotAngle *= -PI / 180;

    // 3: compute relevant values:
    const MFloat freq0 = Strouhal * a_UInfinity();
    const MFloat freq02 = Strouhal;
    for(MInt k = 0; k < noEmbeddedBodies; k++) {
      // when using mu  : has a dimension!
      // when using mu2 : dimensionless!
      mu[k] = freqFactor[k] * freq0 * F2 * PI;
      mu2[k] = freqFactor[k] * freq02 * F2 * PI;
    }

    // if bodyMovementFunction is 6 or 7, adjust start and end angles:

    for(MInt i = 0; i < noEmbeddedBodies; i++) {
      if(bodyToFunction[i] == 6 || bodyToFunction[i] == 7) {
        liftStartAngle1[i] = liftStartAngle1[i] * PI;
        liftEndAngle1[i] = liftEndAngle1[i] * PI - liftStartAngle1[i];
      }
    }

    omega = freqFactor[body] * sqrt(a_TInfinity()) * a_Ma() * PI / (2.0 * amplitude[body]);

    first = false;
  }


  //--------------------------------

  switch(bodyToFunction[body]) {
    case 1: // cosine function

      angle = mu[body] * elapsedTime - liftStartAngle1[body];
      //     cerr << " time: " << elapsedTime << " angle: " << angle << endl;
      switch(returnMode) {
        case 1: // return body center
          if(angle > liftEndAngle1[body]) angle = liftEndAngle1[body];
          if(angle > 0) {
            bodyData[0] = -amplitude[body] * cos(angle);
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = -amplitude[body];
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 2: // return body velocity
          if((angle > 0 && angle < liftEndAngle1[body])) {
            bodyData[0] = mu[body] * amplitude[body] * sin(angle);
            // cerr << " time pos vel body "<< body << " is: " << bodyData[0] << " " << bodyData[1]  << endl;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = F0;
            // cerr<< " time neg vel body " << body << " is: " << bodyData[0] << " " << bodyData[1]  << endl;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if((angle > 0 && angle < liftEndAngle1[body])) {
            bodyData[0] = mu[body] * mu[body] * amplitude[body] * cos(angle);
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }
      break;

    case 2: { // valve lift shifted quadratic sine

      // cerr << " body position valve. body: " << body << endl;
      angle = mu2[body] * elapsedTime;

      switch(returnMode) {
        case 1: // return body center
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = amplitude[body] * POW2(sin(mu2[body] * elapsedTime)) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 2: // return body velocity
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = amplitude[body] * mu2[body] * sin(2 * mu2[body] * elapsedTime) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                          * normal[body * nDim + 0];
            bodyData[1] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                          * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = 2 * amplitude[body] * mu2[body] * mu2[body] * cos(2 * mu2[body] * elapsedTime)
                          * normal[body * nDim + 2];
          } else {
            bodyData[0] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 0];
            bodyData[1] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = 2 * amplitude[body] * mu2[body] * mu2[body] * normal[body * nDim + 2];
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }

      break;
    }


    case 3: // piston movement

      // cerr << " body position piston. body: " << body << endl;
      angle = mu2[body] * elapsedTime;


      switch(returnMode) {
        case 1: // return body center
          if(elapsedTime > F0) {
            bodyData[0] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 0];
            bodyData[1] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 2];
          } else {
            bodyData[0] = -amplitude[body] * normal[body * nDim + 0];
            bodyData[1] = -amplitude[body] * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = -amplitude[body] * normal[body * nDim + 2];
          }
          break;
        case 2: // return body velocity
          if(elapsedTime > F0) {
            bodyData[0] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 0];
            bodyData[1] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = mu2[body] * amplitude[body] * sin(mu2[body] * elapsedTime) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if(elapsedTime > F0) {
            bodyData[0] =
                mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 0];
            bodyData[1] =
                mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] =
                mu2[body] * mu2[body] * amplitude[body] * cos(mu2[body] * elapsedTime) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }

      break;


    case 4: // cosine function with normal

      angle = mu2[body] * elapsedTime;

      switch(returnMode) {
        case 1: // return body center
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = amplitude[body] * cos(angle) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * cos(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * cos(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = amplitude[body] * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * normal[body * nDim + 2];
          }
          break;
        case 2: // return body velocity
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 0];
            bodyData[1] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = -mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if((angle > liftStartAngle1[body] && angle <= liftEndAngle1[body])
             || (angle > liftStartAngle2[body] && angle <= liftEndAngle2[body])) {
            bodyData[0] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 0];
            bodyData[1] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -mu2[body] * mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }
      break;


    case 5:
      // circular motion:
      // initialBodyCenter: center of rotation
      // amplitude        : radius
      // only 2d rotation implemented until now!
      angle = mu2[body] * elapsedTime + circleStartAngle[body];

      switch(returnMode) {
        case 1: // return body center
          if(elapsedTime > F0) {
            bodyData[0] = amplitude[body] * cos(angle);
            bodyData[1] = amplitude[body] * sin(angle);
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = amplitude[body] * cos(circleStartAngle[body]);
            bodyData[1] = amplitude[body] * sin(circleStartAngle[body]);
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 2: // return body velocity
          if(elapsedTime > F0) {
            bodyData[0] = -amplitude[body] * sin(angle);
            bodyData[1] = amplitude[body] * cos(angle);
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if(elapsedTime > F0) {
            bodyData[0] = -amplitude[body] * cos(angle);
            bodyData[1] = -amplitude[body] * sin(angle);
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }
      break;

    case 6: // simplified piston motion
    {
      angle = omega * elapsedTime - liftStartAngle1[body];
      if(angle > liftEndAngle1[body]) angle = F0;
      //  if( domainId() == 0 ) cerr << "angle for piston:" <<angle<< " elapsedTime: " << elapsedTime << endl;

      switch(returnMode) {
        case 1: // return body center
          if(elapsedTime > F0) {
            bodyData[0] = -amplitude[body] * cos(angle) * normal[body * nDim + 0];
            bodyData[1] = -amplitude[body] * cos(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = -amplitude[body] * cos(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = -amplitude[body] * normal[body * nDim + 0];
            bodyData[1] = -amplitude[body] * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = -amplitude[body] * normal[body * nDim + 2];
          }
          break;
        case 2: // return body velocity
          if(elapsedTime > F0) {
            bodyData[0] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 0];
            bodyData[1] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = omega * amplitude[body] * sin(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if(elapsedTime > F0) {
            bodyData[0] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 0];
            bodyData[1] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = omega * omega * amplitude[body] * cos(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }


      break;
    }

    case 7: // valve lift shifted quadratic sine matching piston motion 6

      //       cerr << " body position valve. body: " << body << endl;
      angle = omega * elapsedTime - liftStartAngle1[body];
      if(angle > liftEndAngle1[body]) angle = F0;

      switch(returnMode) {
        case 1: // return body center
          if(angle > F0) {
            bodyData[0] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * POW2(sin(angle)) * normal[body * nDim + 2];
          } else {
            bodyData[0] = 0;
            bodyData[1] = 0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
          }
          break;
        case 2: // return body velocity
          if(angle > F0) {
            bodyData[0] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * omega * sin(2 * angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if(angle > F0) {
            bodyData[0] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 0];
            bodyData[1] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = 2 * amplitude[body] * omega * omega * cos(2 * angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 0];
            bodyData[1] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = 2 * amplitude[body] * omega * omega * normal[body * nDim + 2];
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }

      break;

    case 8:
      // translational movement in normal-direction with constant Ma-number (Ma_trans)
      // Ma_trans is directly given by the amplitude-factor for each body!
      // The body velocity is then based on the free-stream speed of sound.
      angle = sqrt(a_TInfinity());

      if(std::isnan(angle)) {
        cerr << "ERROR in the initialisation of the ls-Solver coupling class!" << endl;
        angle = 0;
      }

      switch(returnMode) {
        case 1: // return body center
          bodyData[0] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 0];
          bodyData[1] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 1];
          IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * angle * elapsedTime * normal[body * nDim + 2];
          break;
        case 2: // return body velocity
          bodyData[0] = amplitude[body] * angle * normal[body * nDim + 0];
          bodyData[1] = amplitude[body] * angle * normal[body * nDim + 1];
          IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * angle * normal[body * nDim + 2];
          break;
        case 3: // return body acceleration
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;

        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }

      break;

    case 9: // sine function with normel (=> periodic motion around the initialBodyCenter!)

      angle = mu2[body] * elapsedTime;

      switch(returnMode) {
        case 1: // return body center
          if(angle >= liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
            bodyData[0] = amplitude[body] * sin(angle) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * sin(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * sin(angle) * normal[body * nDim + 2];
          } else if(angle < liftStartAngle1[body]) {
            bodyData[0] = 0;
            bodyData[1] = 0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = 0;
          } else if(angle > liftEndAngle1[body]) {
            bodyData[0] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 0];
            bodyData[1] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = amplitude[body] * sin(liftEndAngle1[body]) * normal[body * nDim + 2];
          }
          break;
        case 2: // return body velocity
          if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
            bodyData[0] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 0];
            bodyData[1] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3) bodyData[2] = mu2[body] * amplitude[body] * cos(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        case 3: // return body acceleration
          if(angle > liftStartAngle1[body] && angle <= liftEndAngle1[body]) {
            bodyData[0] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 0];
            bodyData[1] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 1];
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -mu2[body] * mu2[body] * amplitude[body] * sin(angle) * normal[body * nDim + 2];
          } else {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          break;
      }
      break;

    case 10: {
      // piston motion equation:
      // reference coordincate system is the cylindeer head!
      // l            : rod length
      // r            : crank radius
      // TDC          : distance at TDC from the cylinder head
      // normal       : motion normal
      // phi          : crank angle in radian

      const MFloat l = liftStartAngle1[body];
      const MFloat TDC = liftEndAngle1[body];
      const MFloat r = amplitude[body];

      MFloat phi = mu2[body] * elapsedTime;

      // consider initial crank-angle:
      phi = phi + m_initialCrankAngle * PI / 180;

      switch(returnMode) {
        case 1: // return body center
          bodyData[0] =
              normal[body * nDim + 0] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));
          bodyData[1] =
              normal[body * nDim + 1] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));
          IF_CONSTEXPR(nDim == 3)
          bodyData[2] =
              normal[body * nDim + 2] * (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))));

          if(lsSolver().domainId() == 0 && fabs(elapsedTime - a_physicalTime()) < 0.0000000001) {
            const MFloat cad = crankAngle(a_physicalTime());
            cerr << "Crank-angle-degree : " << cad << endl;
            if((cad / 45) > 0.989 && (cad / 45) < 1.01) {
              cerr << "Physical-Time : " << a_physicalTime() << endl;
              cerr << "                " << a_physicalTime() * 0.002898783653689 << endl;
            }
            cerr << "Piston-position    : " << (l + r + TDC - (r * cos(phi) + sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                 << " in m " << bodyData[1] * 0.075 << endl;
            cerr << "Piston-velocity    : "
                 << mu2[body]
                        * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))))
                 << endl;
          }

          break;
        case 2: // return body velocity
          bodyData[0] = normal[body * nDim + 0] * mu2[body]
                        * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));

          bodyData[1] = normal[body * nDim + 1] * mu2[body]
                        * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
          IF_CONSTEXPR(nDim == 3)
          bodyData[2] = normal[body * nDim + 2] * mu2[body]
                        * (r * sin(phi) + (POW2(r) * sin(phi) * cos(phi) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));

          break;
        case 3: // return body acceleration
          bodyData[0] =
              normal[body * nDim + 0] * mu2[body] * mu2[body]
              * (r * cos(phi)
                 - (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                 + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / POW3((sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
          bodyData[1] =
              normal[body * nDim + 1] * mu2[body] * mu2[body]
              * (r * cos(phi)
                 - (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                 + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / (POW3(sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
          IF_CONSTEXPR(nDim == 3)
          bodyData[2] =
              normal[body * nDim + 2] * mu2[body] * mu2[body]
              * (r * cos(phi)
                 - (POW2(r) * (POW2(cos(phi)) - POW2(sin(phi))) / (sqrt(POW2(l) - POW2(r) * POW2(sin(phi)))))
                 + (POW4(r) * POW2(sin(phi)) * POW2(cos(phi))) / (POW3(sqrt(POW2(l) - POW2(r) * POW2(sin(phi))))));
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
      }
      break;
    }
    case 11: // valve motion for single-valve engine:
    {
      // reference coordincate system is the cylindeer head!
      // L            : maximum valve lift (2*amplitude[body])
      // normal       : motion normal
      // phi          : crank angle in radian
      // cad          : crank angle degree
      // maxNoCycles  : maximum number of considered cycles


      const MFloat cad = crankAngle(elapsedTime);

      MFloat IO = liftStartAngle1[body];
      MFloat IC = liftEndAngle1[body];
      MFloat EO = liftStartAngle2[body];
      MFloat EC = liftEndAngle2[body];

      // possible values are:
      // IO: -3          // CAD BTDC 3
      // IC: 180+47      // CAD ABDC 47
      // EO: 3*180-47;   // CAD BBDC 47
      // EC: 720+3;      // CAD ATDC 3

      switch(returnMode) {
        case 1: // return body center
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= 0 && cad < IC) {
            MFloat phi_i = (cad - IO) * PI / 180;
            MFloat delta_i = (IC - IO) * PI / 180;
            MFloat freq_i = 1 / (delta_i / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq_i * phi_i - phase));

          } else if(cad > IC && cad < EO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            MFloat phi_e = (cad - EO) * PI / 180;
            MFloat delta_e = (EC - EO) * PI / 180;
            MFloat freq_e = 1 / (delta_e / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq_e * phi_e - phase));
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          if(lsSolver().domainId() == 0 && fabs(elapsedTime - a_physicalTime()) < 0.0000000001) {
            cerr << "Crank-angle-degree : " << cad << endl;
            cerr << "Valve-position     : " << bodyData[0] << " in mm" << bodyData[0] * 75 << endl;
          }

          break;
        case 2: // return body velocity
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= 0 && cad < IC) {
            MFloat phi_i = (cad - IO) * PI / 180;
            MFloat delta_i = (IC - IO) * PI / 180;
            MFloat freq_i = 1 / (delta_i / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq_i * cos(freq_i * phi_i - phase);

          } else if(cad > IC && cad < EO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            MFloat phi_e = (cad - EO) * PI / 180;
            MFloat delta_e = (EC - EO) * PI / 180;
            MFloat freq_e = 1 / (delta_e / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq_e * cos(freq_e * phi_e - phase);
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          if(lsSolver().domainId() == 0 && fabs(elapsedTime - a_physicalTime()) < 0.0000000001) {
            cerr << "Valve-velocity     : " << bodyData[0] << endl;
          }

          break;
        case 3: // return body acceleration
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= 0 && cad < IC) {
            MFloat phi_i = (cad - IO) * PI / 180;
            MFloat delta_i = (IC - IO) * PI / 180;
            MFloat freq_i = 1 / (delta_i / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                          * sin(freq_i * phi_i - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                          * sin(freq_i * phi_i - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_i)
                          * sin(freq_i * phi_i - phase);

          } else if(cad > IC && cad < EO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            MFloat phi_e = (cad - EO) * PI / 180;
            MFloat delta_e = (EC - EO) * PI / 180;
            MFloat freq_e = 1 / (delta_e / 2 / PI);
            MFloat phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                          * sin(freq_e * phi_e - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                          * sin(freq_e * phi_e - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq_e)
                          * sin(freq_e * phi_e - phase);
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
      }
      break;
    }
    case 12: // valve motion for inlet valve:
    {
      // reference coordincate system is the cylindeer head!
      // amplitude    : dimensionless maximum  valve lift
      // normal       : normalised motion normal
      // phi          : crank angle in radian
      // cad          : crank angle degree
      // maxNoCycles  : maximum number of considered cycles

      const MFloat cad = crankAngle(elapsedTime);

      MFloat scaleToMeter = 0.075;

      // valve-timing in crank-angle-degree

      MFloat IO = liftStartAngle1[body];
      MFloat IC = liftEndAngle1[body];

      // possible values are:
      // IO: -24         // 24 CAD BTDC
      // IC: 232         // CAD ATDC

      MFloat phi = 0;
      MFloat delta = -1;
      MFloat freq = 0;
      MFloat phase = -PI / 2;

      switch(returnMode) {
        case 1: // return body center
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= IO && cad <= IC) {
            phi = (cad - IO) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));

          } else if(cad >= 720 + IO) {
            phi = (cad - (720 + IO)) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));

          } else if(cad >= IC && cad < 720 + IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          if(lsSolver().domainId() == 0 && fabs(elapsedTime - a_physicalTime()) < 0.0000000001) {
            cerr << "Crank-angle-degree           : " << cad << endl;
            cerr << "dimensionless-Inlet-Valve-position : " << amplitude[body] * (1 - sin(freq * phi - phase)) << endl;
            cerr << "dimensional-Inlet-Valve-position   : "
                 << amplitude[body] * (1 - sin(freq * phi - phase)) * scaleToMeter << " in m " << endl;
            cerr << "dimensionless-Inlet-Valve-velocity : "
                 << -amplitude[body] * mu2[body] * freq * cos(freq * phi - phase) << endl;
          }

          break;
        case 2: // return body velocity
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= IO && cad <= IC) {
            phi = (cad - IO) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);

          } else if(cad >= 720 + IO) {
            phi = (cad - (720 + IO)) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);

          } else if(cad >= IC && cad < 720 + IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          break;
        case 3: // return body acceleration
          if(cad < IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= IO && cad <= IC) {
            phi = (cad - IO) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);

          } else if(cad >= 720 + IO) {
            phi = (cad - (720 + IO)) * PI / 180;
            delta = (IC - IO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);

          } else if(cad >= IC && cad < 720 + IO) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
      }
      break;
    }
    case 13: // valve motion for outlet valve:
    {
      // reference coordincate system is the cylindeer head!
      // amplitude    : dimensionless maximum  valve lift
      // normal       : motion normal
      // phi          : crank angle in radian
      // cad          : crank angle degree
      // maxNoCycles  : maximum number of considered cycles
      /*
      MFloat cad = mu2[body] * elapsedTime;
      MInt maxNoCycles = 20;

      for (MInt cycle = maxNoCycles; cycle > 0; cycle-- ) {
        if(cad >= 4*PI*cycle) cad = cad - 4*PI*cycle ;
      }
      cad = cad *180/PI;
      */
      const MFloat cad = crankAngle(elapsedTime);

      // valve-timing in crank-angle-degree

      MFloat EO = liftStartAngle1[body];
      MFloat EC = liftEndAngle1[body];
      MFloat scaleToMeter = 0.075;

      // possible values are:
      // EO: 480;    // CAD ATDC
      // EC: 16;     // CAD ATDC

      MFloat phi = 0;
      MFloat delta = -1;
      MFloat freq = 0;
      MFloat phase = -PI / 2;

      switch(returnMode) {
        case 1: // return body center
          if(cad < EO && cad >= EC) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            phi = (cad - EO) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));

          } else if(cad < EC) {
            phi = -(EC - cad) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * (1 - sin(freq * phi - phase));
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * (1 - sin(freq * phi - phase));
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * (1 - sin(freq * phi - phase));

          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          if(lsSolver().domainId() == 0 && fabs(elapsedTime - a_physicalTime()) < 0.0000000001) {
            cerr << "Crank-angle-degree           : " << cad << endl;
            cerr << "dimensionless-Outlet-Valve-position : " << amplitude[body] * (1 - sin(freq * phi - phase)) << endl;
            cerr << "dimensional-Outlet-Valve-position   : "
                 << amplitude[body] * (1 - sin(freq * phi - phase)) * scaleToMeter << " in m " << endl;
            cerr << "dimensionless-Outlet-Valve-velocity : "
                 << -amplitude[body] * mu2[body] * freq * cos(freq * phi - phase) << endl;
          }

          break;
        case 2: // return body velocity
          if(cad < EO && cad >= EC) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            phi = (cad - EO) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);


          } else if(cad < EC) {
            phi = -(EC - cad) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = -normal[body * nDim + 0] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            bodyData[1] = -normal[body * nDim + 1] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = -normal[body * nDim + 2] * amplitude[body] * mu2[body] * freq * cos(freq * phi - phase);
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }

          break;
        case 3: // return body acceleration
          if(cad < EO && cad >= EC) {
            bodyData[0] = F0;
            bodyData[1] = F0;
            IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
          } else if(cad >= EO && cad < 720) {
            phi = (cad - EO) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);

          } else if(cad < EC) {
            phi = -(EC - cad) * PI / 180;
            delta = (720 + EC - EO) * PI / 180;
            freq = 1 / (delta / 2 / PI);
            phase = -PI / 2;

            bodyData[0] = normal[body * nDim + 0] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            bodyData[1] = normal[body * nDim + 1] * amplitude[body] * mu2[body] * mu2[body] * POW2(freq)
                          * sin(freq * phi - phase);
            IF_CONSTEXPR(nDim == 3)
            bodyData[2] = normal[body * nDim + 2] * amplitude[body] * mu2[body] * mu2[body]
                          * (1 + POW2(freq) * sin(freq * phi - phase));
          } else {
            mTerm(1, AT_, "Unexpected crank-angle-degree!");
          }
          break;
        default:
          bodyData[0] = F0;
          bodyData[1] = F0;
          IF_CONSTEXPR(nDim == 3) bodyData[2] = F0;
      }
      break;
    }
    default:
      mTerm(1, AT_, "function type not implemented. Please check bodyMovementFunctions property!");
  }

  // add the initialBodyCenter to the bodyDato for the body-positioning:
  if(returnMode == 1) {
    for(MInt dir = 0; dir < nDim; dir++) {
      bodyData[dir] += initialBodyCenter[body * nDim + dir];
    }
  }

  // rotate the final result around the z-axis

  if(rotAngle > 0 || rotAngle < 0) {
    MFloat tmp0 = bodyData[0] * cos(rotAngle) + bodyData[1] * sin(rotAngle);
    MFloat tmp1 = bodyData[1] * cos(rotAngle) - bodyData[0] * sin(rotAngle);
    bodyData[0] = tmp0;
    bodyData[1] = tmp1;
    bodyData[2] = bodyData[2];
  }
}

/** \brief help-function for engine calculations which returns the crank-angle
 *         for a given time
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
MFloat CouplingLsFv<nDim, SysEqn>::crankAngle(MFloat elapsedTime) {
  TRACE();

  MFloat& Strouhal = m_static_crankAngle_Strouhal;
  MBool& first = m_static_crankAngle_first;

  if(first) {
    Strouhal = Context::getSolverProperty<MFloat>("Strouhal", m_lsSolverId, AT_, &Strouhal);
    first = false;
  }

  const MFloat mu2 = Strouhal * F2 * PI;
  MFloat cad = mu2 * elapsedTime;
  const MInt maxNoCycles = 20;

  for(MInt cycle = maxNoCycles; cycle > 0; cycle--) {
    if(cad >= 4 * PI * cycle) cad = cad - 4 * PI * cycle;
  }
  cad = cad * 180 / PI;

  // consider initial crank angle
  cad = cad + m_initialCrankAngle;

  return cad;
}

template class CouplingLsFv<2, FvSysEqnNS<2>>;
template class CouplingLsFv<3, FvSysEqnNS<3>>;
template class CouplingLsFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingLsFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingLsFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class CouplingLsFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplingLsFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplingLsFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
