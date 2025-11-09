// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbrb.h"

#include "GRID/cartesiannetcdf.h"
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "globals.h"
#include "globalvariables.h"

#include <algorithm>
#include <stack>
#include <vector>

/**
 * \brief C'tor for the lattice Boltzmann rigid bodies coupler
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param couplingId[in] Unique id to identify the couplre
 * \param lb[in]         Pointer to the lattice Botzmann solver
 * \param rb[in]         Pointer to the rigid bodies solver
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbRb<nDim, nDist, SysEqn>::LbRb(MInt couplingId, LbSolver* lb, RBodies* rb)
  : Coupling(couplingId), CouplingLb(couplingId, lb), CouplingRb(couplingId, rb) {
  TRACE();

  // Init timers as the first action
  initTimers();

  initData();
  readProperties();
  checkProperties();
  updateGeometry();

  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

/**
 * \brief D'tor for the lattice Boltzmann rigid bodies coupler
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbRb<nDim, nDist, SysEqn>::~LbRb() {
  TRACE();

  RECORD_TIMER_STOP(m_timers[Timers::Class]);

  averageTimer();
}

/**
 * \brief Creates all timers and subtimers
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::initTimers() {
  TRACE();

  // Create timer group and coupler timer, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timers[Timers::TimerGroup], "LbRb");
  NEW_TIMER_NOCREATE(m_timers[Timers::Class], "Total object lifetime", m_timers[Timers::TimerGroup]);
  RECORD_TIMER_START(m_timers[Timers::Class]);

  // Create and start constructor timer
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Constructor], "Constructor", m_timers[Timers::Class]);
  RECORD_TIMER_START(m_timers[Timers::Constructor]);

  // Create Couple timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CouplePostRb], "Couple post RB", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ConstructGField], "ConstructGField", m_timers[Timers::CouplePostRb]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Preparation], "Prepare construct", m_timers[Timers::ConstructGField]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FindBoundaryCells], "FindBoundaryCells", m_timers[Timers::CouplePostRb]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PreCouple], "Pre Couple", m_timers[Timers::FindBoundaryCells]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FindBoundaryMapping], "Find Boundary Mapping",
                         m_timers[Timers::FindBoundaryCells]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitSolidDomain], "InitSolidDomain", m_timers[Timers::CouplePostRb]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SetBoundaryVelocity], "SetBoundaryVelocity", m_timers[Timers::CouplePostRb]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CouplePostLb], "Couple post LB", m_timers[Timers::Class]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CreateComm], "CreateComm", m_timers[Timers::CouplePostLb]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ApplyBC], "ApplyBC", m_timers[Timers::CouplePostLb]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ComputeBodyForces], "ComputeBodyForces", m_timers[Timers::CouplePostLb]);
}


/**
 * \brief Initializes coupler specific data
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::initData() {
  TRACE();

  // update RB grid meta data
  bodies().updateMaxLevel(lbSolver().maxLevel());

  // Conversion factors
  const MFloat dx = a_cellLengthAtLevel(lbSolver().maxLevel());

  conversionLbRb.length = dx;
  conversionRbLb.length = 1 / conversionLbRb.length;

  conversionLbRb.velocity = sqrt(3) / a_Ma();
  conversionRbLb.velocity = 1 / conversionLbRb.velocity;

  conversionLbRb.time = conversionLbRb.length / conversionLbRb.velocity;
  conversionRbLb.time = 1 / conversionLbRb.time;

  conversionLbRb.force = POW2(conversionLbRb.velocity) * POW2(conversionLbRb.length);
  conversionRbLb.force = 1 / conversionLbRb.force;

  conversionLbRb.torque = conversionLbRb.force * conversionLbRb.length;
  conversionRbLb.torque = 1 / conversionLbRb.torque;

  bodies().setTimestep(conversionLbRb.time);

#ifndef NDEBUG
  if(!lbSolver().domainId()) {
    std::cout << "CONVERSION" << std::endl
              << "LENGTH LB RB " << conversionLbRb.length << std::endl
              << "LENGTH RB LB " << conversionRbLb.length << std::endl
              << "VEL LB RB " << conversionLbRb.velocity << std::endl
              << "VEL RB LB " << conversionRbLb.velocity << std::endl
              << "TIME LB RB " << conversionLbRb.time << std::endl
              << "TIME RB LB " << conversionRbLb.time << std::endl
              << "FORCE LB RB " << conversionLbRb.force << std::endl
              << "FORCE RB LB " << conversionRbLb.force << std::endl;
  }
#endif
}

/**
 * \brief Checks property-data which is read in by both ls-and Lb-Solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::checkProperties() {
  TRACE();

  lbSolver().m_noLevelSetsUsedForMb = 1;
  lbSolver().m_maxNoSets = 1;
}

/**
 * \brief Read all relevant properties
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::readProperties() {}

/**
 * \brief Initialize solver data needed for this coupling
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::init() {
  TRACE();

  // LB solver
  lbSolver().initializeMovingBoundaries();
  lbBndCnd().initializeBndMovingBoundaries();

  // RB solver
  bodies().setTimestep(conversionLbRb.time);

  // Construct level set for initial adaptation
  RECORD_TIMER_START(m_timers[Timers::CouplePostRb]);
  RECORD_TIMER_START(m_timers[Timers::ConstructGField]);
  constructGField();
  RECORD_TIMER_STOP(m_timers[Timers::ConstructGField]);
  RECORD_TIMER_STOP(m_timers[Timers::CouplePostRb]);
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::postAdaptation() {
  // Update conversion factors
  initData();
  // Construct level set to initialize newly refined cells
  RECORD_TIMER_START(m_timers[Timers::CouplePostRb]);
  RECORD_TIMER_START(m_timers[Timers::ConstructGField]);
  constructGField();
  RECORD_TIMER_STOP(m_timers[Timers::ConstructGField]);
  RECORD_TIMER_STOP(m_timers[Timers::CouplePostRb]);
}

/**
 * \brief Coupling between solver substeps
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] rs              Current recipe step
 * \param[in] sid             Current solver id after which this function is called
 * \param[in] solverCompleted Current completion status of all solvers
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::finalizeAdaptation(const MInt solverId) {
  // LB bnd cnd is restarted after adaptation and needs to be initialized again
  if(solverId == lbSolver().solverId()) {
    lbSolver().initializeMovingBoundaries();
    lbBndCnd().initializeBndMovingBoundaries();
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::subCouple(MInt /*rs*/, MInt sid, std::vector<MBool>& solverCompleted) {
  // Do nothing if both solver already finished
  if(solverCompleted[lbSolver().solverId()] && solverCompleted[bodies().solverId()]) {
    return;
  }

  if(sid == bodies().solverId()) {
    // POST RB
    RECORD_TIMER_START(m_timers[Timers::CouplePostRb]);
    RECORD_TIMER_START(m_timers[Timers::ConstructGField]);
    constructGField();
    RECORD_TIMER_STOP(m_timers[Timers::ConstructGField]);

    std::vector<MInt> maxGCellLevels(lbSolver().m_maxNoSets);
    for(MInt set = 0; set < lbSolver().m_maxNoSets; set++) {
      maxGCellLevels[set] = lbSolver().maxLevel();
    }

    RECORD_TIMER_START(m_timers[Timers::FindBoundaryCells]);
    RECORD_TIMER_START(m_timers[Timers::PreCouple]);
    lbSolver().preCoupleLs(maxGCellLevels);
    RECORD_TIMER_STOP(m_timers[Timers::PreCouple]);

    RECORD_TIMER_START(m_timers[Timers::FindBoundaryMapping]);
    lbSolver().createBndryToBodyMapping(bndryToBodyMapping, bodyToBndryMapping);
    RECORD_TIMER_STOP(m_timers[Timers::FindBoundaryMapping]);
    RECORD_TIMER_STOP(m_timers[Timers::FindBoundaryCells]);

    RECORD_TIMER_START(m_timers[Timers::InitSolidDomain]);
    initializeSolidDomain();
    RECORD_TIMER_STOP(m_timers[Timers::InitSolidDomain]);

    RECORD_TIMER_START(m_timers[Timers::SetBoundaryVelocity]);

    struct {
      MFloat velocity;
      MFloat angularVelocity;
      MFloat length;
      MFloat temperature;
    } conversion{conversionRbLb.velocity, conversionLbRb.time, 1.0 / a_cellLengthAtLevel(lbSolver().maxLevel()), 1.0};

    maia::coupling::setBoundaryVelocity<nDim>(bodies(), a_mbCell(), bodyToBndryMapping, conversion);

    RECORD_TIMER_STOP(m_timers[Timers::SetBoundaryVelocity]);
    RECORD_TIMER_STOP(m_timers[Timers::CouplePostRb]);

  } else if(sid == lbSolver().solverId()) {
    // POST LB
    RECORD_TIMER_START(m_timers[Timers::CouplePostLb]);
    RECORD_TIMER_START(m_timers[Timers::CreateComm]);
    lbBndCnd().createMBComm();
    RECORD_TIMER_STOP(m_timers[Timers::CreateComm]);

    RECORD_TIMER_START(m_timers[Timers::ApplyBC]);
    lbBndCnd().postCouple();
    RECORD_TIMER_STOP(m_timers[Timers::ApplyBC]);

    RECORD_TIMER_START(m_timers[Timers::ComputeBodyForces]);

    struct {
      MFloat force;
      MFloat length;
    } const conversion{conversionLbRb.force, 1.0};

    if(a_noCollectorBodies() > 0) {
      maia::coupling::setBoundaryForceAndTorque<nDim>(a_mbCell(), bodies(), bndryToBodyMapping, conversion);
    }

#ifndef NDEBUG
    for(MInt b = 0; b < a_noCollectorBodies(); b++) {
      for(MInt n = 0; n < nDim; n++) {
        if(std::isnan(bodies().a_bodyForce(b, n))) {
          std::cout << "SUM FORCE IS NAN AFTER SETBOUNDANDFORCE AT TIMESETEP: " << lbSolver().getCurrentTimeStep()
                    << std::endl;
          // TERMM(1, "SUM FORCE AFTER IS NAN");
        }
      }
    }
#endif

    RECORD_TIMER_STOP(m_timers[Timers::ComputeBodyForces]);
    RECORD_TIMER_STOP(m_timers[Timers::CouplePostLb]);
  }
}

/**
 * \brief Get body velocity converted to LB units
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] body         Body id
 * \param[in] bodyVelocity Body velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::getBodyVelocity(const MInt body, MFloat* const bodyVelocity) {
  for(MInt n = 0; n < nDim; n++) {
    bodyVelocity[n] = bodies().a_bodyVelocity(body, n) * conversionRbLb.velocity;
  }
}

/**
 * \brief Get angular body velocity converted to LB units
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] body         Body id
 * \param[in] bodyVelocity Angular body velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::getBodyAngularVelocity(const MInt body, MFloat* const angularVelocity) {
  for(MInt n = 0; n < nRot; n++) {
    angularVelocity[n] = bodies().a_angularVelocity(body, n) / conversionRbLb.time;
  }
}

/**
 * \brief Updates the member-variables in the geometry-intersection class
 *
 * \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::updateGeometry() {
  TRACE();

  lbSolver().m_noEmbeddedBodies = 1;
  lbSolver().m_geometryIntersection->m_noLevelSetsUsedForMb = 1;
  lbSolver().m_geometryIntersection->m_noEmbeddedBodies = bodies().size();
}

/**
 * \brief Initialize cells which are inisde the solid domain or just entered the fluid domain
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::initializeSolidDomain() {
  TRACE();

  if(lbSolver().noNeighborDomains() > 0) {
    lbSolver().exchangeOldDistributions();
  }

  MFloatScratchSpace bodyVelocities(a_noCollectorBodies(), nDim, AT_, "bodyVelocities");

  // Get the body Velocity for each embedded body
  for(MInt body = 0; body < a_noCollectorBodies(); body++) {
    getBodyVelocity(body, &bodyVelocities(body, 0));
  }

  for(MInt i = 0; i < a_noCells(); i++) {
    // Regular fluid cell
    if(a_isActive(i) && a_wasActive(i)) {
      continue;
    }

    // Regular Solid cell
    if(!a_isActive(i)) {
      if(!lbSolver().c_isLeafCell(i)) {
        continue;
      }

      // determine the body to which the cell belongs, to set the right body velocity
      MInt bodyId = -1;
      MInt setOfBody = 0;
      for(MInt set = lbSolver().m_levelSetId; set < lbSolver().m_maxNoSets; set++) {
        if(a_associatedBodyIdsMb(i, set) >= 0) {
          bodyId = a_associatedBodyIdsMb(i, set);
          setOfBody = set;
          break;
        }
      }

      // if body was deleted, the halo cells that were occupied by the body in the previous timeStep, must be refilled
      // now
      if(bodyId == -1 && lbSolver().a_isHalo(i) && bodies().m_bodyWasDeleted) {
        lbBndCnd().refillEmergedCell(i);
        lbSolver().a_isActive(i) = 1;
        continue;
      }

      ASSERT(bodyId > -1, "No valid bodyId for solid cell! (bodyId=" << bodyId << ") (" << lbSolver().a_coordinate(i, 0)
                                                                     << " " << lbSolver().a_coordinate(i, 1) << ") ("
                                                                     << i << " " << lbSolver().c_globalId(i) << ")"
                                                                     << " isHalo " << lbSolver().a_isHalo(i));

      // the Velocity of the deactivated cell is set to the body velocity
      // the Density is set to 1.0
      if((bodyId >= 0) && (bodyId < a_noCollectorBodies()) && (a_levelSetFunctionMb(i, setOfBody) < 0)) {
        for(MInt j = 0; j < nDim; j++) {
          a_variable(i, j) = bodyVelocities(bodyId, j);
          a_oldVariable(i, j) = bodyVelocities(bodyId, j);
        }
      } else {
        for(MInt j = 0; j < a_noVariables(); j++) {
          a_variable(i, j) = F0;
          a_oldVariable(i, j) = F0;
        }
      }
      a_variable(i, a_pvrho()) = 1.0;
      a_oldVariable(i, a_pvrho()) = 1.0;

      if(a_isThermal()) {
        a_variable(i, a_pvt()) = a_initTemperatureKelvin();
      }

      // Distributions are not set since they are not used for solid nodes
    } else {
      // New fluid cell
      lbBndCnd().refillEmergedCell(i);
    }
  }
  bodies().m_bodyWasDeleted = false;
}

/**
 * \brief Dispatch function for the body specific construction of the level set
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::constructGField() {
  if(bodies().a_bodyType() == 1)
    constructGField_<1>();
  else if(bodies().a_bodyType() == 2)
    constructGField_<2>();
  else if(bodies().a_bodyType() == 3)
    constructGField_<3>();
  else if(bodies().a_bodyType() == 4)
    constructGField_<4>();
  else if(bodies().a_bodyType() == 7)
    constructGField_<7>();
  else
    mTerm(1, AT_, "Body type not implemented!");
}

/**
 * \brief Constructs the level-set field after each time step
 *
 * \author Lennart Schneiders, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \tparam bodyType Body type for which the level set is constructed
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt bodyType>
void LbRb<nDim, nDist, SysEqn>::constructGField_() {
  TRACE();

  RECORD_TIMER_START(m_timers[Timers::Preparation]);

  /* bandwidth at "level" has the no. of cells at "level",
  whose total width equates to all summed up bands up to "level" */
  MInt level = lbSolver().maxUniformRefinementLevel() - 1;
  MFloat m_bodyDistThreshold = 0.0;
  if(lbSolver().m_adaptation) {
    m_bodyDistThreshold = lbSolver().m_bandWidth[level] * a_cellLengthAtLevel(level);
    // m_bodyDistThreshold =  4 * a_cellLengthAtLevel(lbSolver().maxLevel());
  } else {
    m_bodyDistThreshold = 2 * a_cellLengthAtLevel(lbSolver().maxLevel());
  }

  // To be changed later ~jv
  const MInt m_noLevelSetsUsedForMb = 1;

  // Reset all sets
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
      a_associatedBodyIdsMb(cellId, set) = -1;
      a_levelSetFunctionMb(cellId, set) = m_bodyDistThreshold + 1e-14;
    }
  }

  const MInt noRelevantBodies = a_noCollectorBodies();

  if(noRelevantBodies == 0) {
    RECORD_TIMER_STOP(m_timers[Timers::Preparation]);
    return;
  }

  RECORD_TIMER_STOP(m_timers[Timers::Preparation]);

  // construct vector with connecting Bodies
  std::vector<MInt> collectorBodyIds(a_noCollectorBodies());
  for(MInt i = 0; i < a_noCollectorBodies(); i++) {
    collectorBodyIds[i] = i;
  }

  for(MInt i = 0; i < noMinCells(); i++) {
    const MInt cellId = minCell(i);
    const MFloat minDist = m_bodyDistThreshold + 1e-14;

    for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
      a_associatedBodyIdsMb(cellId, set) = -1;
      a_levelSetFunctionMb(cellId, set) = minDist;
    }

    if(a_noCollectorBodies() > 0) {
      descendLevelSetValue<bodyType>(cellId, collectorBodyIds.data(), noRelevantBodies);
    }

    MInt parentId = a_parentId(cellId);
    while(parentId > -1) {
      for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
        a_levelSetFunctionMb(parentId, set) = a_levelSetFunctionMb(cellId, set);
        a_associatedBodyIdsMb(parentId, set) = a_associatedBodyIdsMb(cellId, set);
      }
      parentId = a_parentId(parentId);
    }
  }
}

/**
 * \brief Descend the level set value to the cells children
 *
 * This function is used recursively starting from the min level cells.
 * The descend is skipped if the distance on the lower level is greater than a given threshold.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] cellId  Cell id of the coarse cell
 * \param[in] bodyIds List of all relevant body ids
 * \param[in] bodyCnt Number of bodies
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt bodyType>
void LbRb<nDim, nDist, SysEqn>::descendLevelSetValue(const MInt cellId, const MInt* const bodyIds, const MInt bodyCnt) {
  const MFloat minLevelThreshold = a_cellLengthAtLevel(lbSolver().a_level(cellId));

  MBool skipDescend = true;

  for(MInt b = 0; b < bodyCnt; b++) {
    const MInt k = bodyIds[b];
    const MInt set = 0;

    MFloat dist = bodies().template getDistance<bodyType>(&lbSolver().a_coordinate(cellId, 0), k);

    if(dist < minLevelThreshold) {
      skipDescend = false;
    }

    if(fabs(dist) < fabs(a_levelSetFunctionMb(cellId, set)) || (dist < F0)) {
      a_levelSetFunctionMb(cellId, set) = dist;
      a_associatedBodyIdsMb(cellId, set) = k;
    }
  }

  if(skipDescend || !lbSolver().grid().tree().hasChildren(cellId)) {
    return;
  }

  // Recursively descend via children
  constexpr MInt maxNoChildren = nDim == 3 ? 8 : 4;
  for(MInt child = 0; child < maxNoChildren; child++) {
    if(a_childId(cellId, child) < 0) continue;
    descendLevelSetValue<bodyType>(a_childId(cellId, child), bodyIds, bodyCnt);
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbRb<nDim, nDist, SysEqn>::averageTimer() {
  TRACE();
  if(!lbSolver().grid().isActive()) return;

  // Get timer operation
  m_timerType = "max";
  m_timerType = Context::getSolverProperty<MString>("timerType", lbSolver().solverId(), AT_, &m_timerType);

  // 0) map timer ids for safety
  const MInt noTimers = m_timers.size();

  // 1) fill buffer with local timer values
  std::vector<MFloat> timerValues_;
  timerValues_.reserve(noTimers);
  for(MInt i = 0; i < noTimers; i++) {
    timerValues_.emplace_back(RETURN_TIMER_TIME(m_timers[i]));
  }

  // 2) collect values from all ranks
  if(m_timerType == "average") {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_SUM,
                  lbSolver().mpiComm(), AT_, "MPI_IN_PLACE", "timerValues_");
  } else {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_MAX,
                  lbSolver().mpiComm(), AT_, "MPI_IN_PLACE", "timerValues_");
  }

  // 3) perform averaging on timer and4) set new timer values
  if(m_timerType == "average") {
    const MInt noDomains_ = lbSolver().noDomains();
    for(MInt i = 0; i < noTimers; i++) {
      const MFloat meanValue = timerValues_[i] / noDomains_;
      SET_RECORD(m_timers[i], meanValue);
    }
  } else {
    for(MInt i = 0; i < noTimers; i++) {
      SET_RECORD(m_timers[i], timerValues_[i]);
    }
  }
}
template class LbRb<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class LbRb<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LbRb<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
