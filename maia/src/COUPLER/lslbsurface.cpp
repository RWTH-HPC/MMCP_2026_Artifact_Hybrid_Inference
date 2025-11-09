// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

//#define COUPLING_DEBUG_
#include "lslbsurface.h"
#include "GRID/cartesiannetcdf.h"
#include "IO/parallelio.h"
#include "LB/lblatticedescriptor.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "globals.h"
#include "globalvariables.h"

#include <algorithm>
#include <stack>
#include <vector>

using namespace std;

template <MInt nDim, MInt nDist, class SysEqn>
LsLbSurface<nDim, nDist, SysEqn>::LsLbSurface(MInt couplingId, LsSolver* ls, LbSolver* lb)
  : Coupling(couplingId), CouplingLS<nDim>(couplingId, ls), CouplingLB<nDim, nDist, SysEqn>(couplingId, lb) {
  TRACE();

  // Init timers as the first action
  initTimers();

  readProperties();

  initData();

  checkProperties();

  updateGeometry();

  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

template <MInt nDim, MInt nDist, class SysEqn>
LsLbSurface<nDim, nDist, SysEqn>::~LsLbSurface() {
  TRACE();

  RECORD_TIMER_STOP(m_timers[Timers::Class]);
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::initTimers() {
  TRACE();

  // Create timer group and coupler timer, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timers[Timers::TimerGroup], "Coupler LS LB Surface");
  NEW_TIMER_NOCREATE(m_timers[Timers::Class], "Total object lifetime", m_timers[Timers::TimerGroup]);
  RECORD_TIMER_START(m_timers[Timers::Class]);

  // Create and start constructor timer
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Constructor], "Constructor", m_timers[Timers::Class]);
  RECORD_TIMER_START(m_timers[Timers::Constructor]);

  // Create PreCouple timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PreCouple], "PreCouple", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TransferToLevelSet], "TransferToLevelSet", m_timers[Timers::PreCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SetExtensionVelocity], "SetExtensionVelocity",
                         m_timers[Timers::TransferToLevelSet]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ExtendVelocity], "ExtendVelocity", m_timers[Timers::TransferToLevelSet]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::FindBoundaryCells], "FindBoundaryCells", m_timers[Timers::PreCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitOutsideDomain], "InitOutsideDomain", m_timers[Timers::PreCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SetBoundaryCondition], "SetBoundaryCondition", m_timers[Timers::PreCouple]);

  // Create PostCouple timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PostCouple], "PostCouple", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ComputeBoundaryValues], "ComputeBoundaryValues",
                         m_timers[Timers::PostCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CreateComm], "CreateComm", m_timers[Timers::PostCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ApplyBoundaryCondition], "ApplyBoundaryCondition",
                         m_timers[Timers::PostCouple]);
}


/** \brief Initialize coupling-class-specific Data
 *    \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::initData() {
  TRACE();

  // solver-specific data:
  m_lbSolverId = lbSolver().solverId();
  m_lsSolverId = lsSolver().solverId();

  const MFloat dx = a_cellLengthAtLevel(lbSolver().maxLevel());

  conversionLbLs.velocity = sqrt(3) / a_Ma();
  conversionLsLb.velocity = 1 / conversionLbLs.velocity;

  conversionLbLs.length = dx;
  conversionLsLb.length = 1 / conversionLbLs.length;

  conversionLbLs.time = conversionLbLs.length / conversionLbLs.velocity;
  conversionLsLb.time = 1 / conversionLbLs.time;

  conversionLbLs.pressure = POW2(conversionLbLs.velocity);
  conversionLsLb.pressure = 1 / conversionLbLs.pressure;

  if(!lbSolver().domainId()) {
    std::cout << "CONVERSION" << std::endl
              << "VEL LB LS " << conversionLbLs.velocity << std::endl
              << "VEL LS LB " << conversionLsLb.velocity << std::endl
              << "LEN LB LS " << conversionLbLs.length << std::endl
              << "LEN LS LB " << conversionLsLb.length << std::endl
              << "TIME LB LS " << conversionLbLs.time << std::endl
              << "TIME LS LB " << conversionLsLb.time << std::endl
              << "PRESSURE LB LS " << conversionLbLs.pressure << std::endl
              << "PRESSURE LS LB " << conversionLsLb.pressure << std::endl;
  }

  lsSolver().m_timeStep = conversionLbLs.time;

  const MFloat nu = a_Ma() * LBCS / a_Re() * lbSolver().m_referenceLength;

  m_gravity = POW2(m_Ga) * POW2(nu) / POW3(lbSolver().m_referenceLength);

  m_surfaceTension = 1.0 * m_gravity * POW2(lbSolver().m_referenceLength) / m_Eo;

  // Some quality of life output

  const MFloat tau = (1 + 6 * nu) / 2.0;

  stringstream ss;
  const MInt maxLineLength = 256;
  MChar b[maxLineLength];
  const MString divider = "-------------------------------------------------------------------------";

  ss << "\n";
  ss << "FREE SURFACE PROBLEM SUMARRY"
     << "\n";
  ss << divider << "\n";

  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "Mach", a_Ma());
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "Reynolds", a_Re());
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "Galileo", m_Ga);
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "Eotvos", m_Eo);
  ss << b;

  ss << "\n";
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "viscosity nu", nu);
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "relaxation tau", tau);
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "gravity g", m_gravity);
  ss << b;
  snprintf(b, maxLineLength, " | %-45s | %-20E | \n", "surface tension sigma", m_surfaceTension);
  ss << b;

  if(!lbSolver().domainId()) {
    cout << ss.str() << std::endl;
  }
}

/** \brief Checks property-data which is read in by both ls-and Lb-Solver
 *    \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::checkProperties() {
  TRACE();

  lbSolver().m_noLevelSetsUsedForMb = 1;
  if(lsSolver().m_maxNoSets > 1) {
    lbSolver().m_noLevelSetsUsedForMb = lsSolver().m_maxNoSets;
  }

  lbSolver().m_maxNoSets = lsSolver().m_maxNoSets;
  lbSolver().m_constructGField = lsSolver().m_constructGField;

  // Check that the lb-solver-cell-count is correct!
  ASSERT(lbSolver().a_noCells(), "");
}

/** \brief reads lsfvmb-coupling-specific data
 *    \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::readProperties() {
  TRACE();

  m_Ga = 0.0;
  m_Ga = Context::getBasicProperty<MFloat>("Ga", AT_, &m_Ga);


  m_Eo = 0.0;
  m_Eo = Context::getBasicProperty<MFloat>("Eo", AT_, &m_Eo);

  m_initCurvature = 0.0;
  m_initCurvature = Context::getBasicProperty<MFloat>("initCurvature", AT_, &m_initCurvature);

  m_initHeight = 0.0;
  m_initHeight = Context::getBasicProperty<MFloat>("initHeight", AT_, &m_initHeight);
  m_initHeight /= lbSolver().c_cellLengthAtLevel(lbSolver().maxLevel());
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::init() {
  TRACE();

  lbSolver().initializeMovingBoundaries();
  lbBndCnd().initializeBndMovingBoundaries();
}


template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::preCouple(const MInt recipeStep) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::PreCouple]);

  // Pre LS
  if(recipeStep == 0) {
    if(!lbSolver().domainId()) {
      std::cout << "Pre LS couple ..." << std::endl;
    }

    RECORD_TIMER_START(m_timers[Timers::TransferToLevelSet]);
    // TODO labels:COUPLER Decide based on number of bnd cells
    // FIXME labels:COUPLER julian
    if(lbSolver().getCurrentTimeStep() > 1) {
      // Reset, set and extend transport velocity for the LS
      // FIXME labels:COUPLER julian
      // lsSolver().resetExtensionVelocity();

      RECORD_TIMER_START(m_timers[Timers::SetExtensionVelocity]);
      setExtensionVelocity();
      RECORD_TIMER_STOP(m_timers[Timers::SetExtensionVelocity]);

      RECORD_TIMER_START(m_timers[Timers::ExtendVelocity]);
      // FIXME labels:COUPLER julian
      // lsSolver().extendVelocity(0);
      RECORD_TIMER_STOP(m_timers[Timers::ExtendVelocity]);
    }
    RECORD_TIMER_STOP(m_timers[Timers::TransferToLevelSet]);
  }

  // Pre LB
  else {
    if(!lbSolver().domainId()) {
      std::cout << "Pre LB couple ..." << std::endl;
    }

    std::vector<MInt> maxGCellLevels(lsSolver().m_maxNoSets);
    for(MInt set = 0; set < lsSolver().m_maxNoSets; set++) {
      maxGCellLevels[set] = lsSolver().a_maxGCellLevel(set);
    }
    RECORD_TIMER_START(m_timers[Timers::FindBoundaryCells]);
    lbSolver().preCoupleLs(maxGCellLevels);
    RECORD_TIMER_STOP(m_timers[Timers::FindBoundaryCells]);

    RECORD_TIMER_START(m_timers[Timers::InitOutsideDomain]);
    refillEmergedCells();
    RECORD_TIMER_STOP(m_timers[Timers::InitOutsideDomain]);

    // FIXME labels:COUPLER Curvature semms wrong for the first time step
    for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
      a_mbCell().density(mbCell) = 1.0;
    }
    if(lbSolver().getCurrentTimeStep() == 1) {
      return;
    }

    RECORD_TIMER_START(m_timers[Timers::SetBoundaryCondition]);
    MFloatScratchSpace curvatureInterp(a_mbCell().size(), AT_, "curvature");

    interpolateCurvature(curvatureInterp);
    interpolateNormal();

    if(a_mbCell().size() > 0) {
      std::cout << "curv " << curvatureInterp[0] << std::endl;
    }

    // Compute and transfer pressure from curvature (Laplace)
    MFloat l_curvature_sum = 0;
    MFloat l_rho_sum = 0;
    MFloat l_rho_hyd = 0;
    MFloat l_drho_curv = 0;

    for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
      l_curvature_sum += curvatureInterp[mbCell];
    }

    const MInt l_count = a_mbCell().size();
    MInt g_count = 0;
    MPI_Allreduce(&l_count, &g_count, 1, MPI_INT, MPI_SUM, lbSolver().mpiComm(), AT_, "count", "count");

    MFloat g_curvature_sum = 0;
    MPI_Allreduce(&l_curvature_sum, &g_curvature_sum, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "curv",
                  "curv");

    const MFloat dpCorr = 0.1 * (POW3(g_curvature_sum / g_count / m_initCurvature) - 1.0);

    for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
      const MFloat curvature = curvatureInterp[mbCell];

      const MFloat dx = a_cellLengthAtLevel(lbSolver().maxLevel());

      MFloat dpLB = -m_initHeight * m_gravity + dpCorr;

      l_rho_hyd += 1 + 3.0 * dpLB;

      // TODO labels:COUPLER Check Factor convention 2D/3D ...
      IF_CONSTEXPR(nDim == 2) { dpLB += 2 * (curvature - m_initCurvature) * dx * m_surfaceTension; }
      else {
        dpLB += 1 * (curvature - m_initCurvature) * dx * m_surfaceTension;
        l_drho_curv += 3 * (curvature - m_initCurvature) * dx * m_surfaceTension;
      }

      const MFloat rhoLB = 1 + 3.0 * dpLB;

      a_mbCell().density(mbCell) = rhoLB;

      l_rho_sum += rhoLB;
    }

    MFloat g_rho_sum = 0;
    MPI_Allreduce(&l_rho_sum, &g_rho_sum, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "curv", "curv");

    MFloat g_rho_hyd = 0;
    MPI_Allreduce(&l_rho_hyd, &g_rho_hyd, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "curv", "curv");

    MFloat g_drho_curv = 0;
    MPI_Allreduce(&l_drho_curv, &g_drho_curv, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "curv", "curv");

    if(!lbSolver().domainId()) {
      std::cout << " AVG curv " << g_curvature_sum / g_count << " AVG rho " << g_rho_sum / g_count << " AVG hyd rho "
                << g_rho_hyd / g_count << " AVG curv drho " << g_drho_curv / g_count << std::endl;
    }

    RECORD_TIMER_STOP(m_timers[Timers::SetBoundaryCondition]);
  }
  RECORD_TIMER_STOP(m_timers[Timers::PreCouple]);
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::postCouple(const MInt recipeStep) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::PostCouple]);

  // Post LS
  if(recipeStep == 0) {
    if(!lbSolver().domainId()) {
      std::cout << "Post LS couple ..." << std::endl;
    }

    // TODO labels:COUPLER,LS,LB Transfer LS FROM LS to LB in a more generalized manner
    // TODO labels:COUPLER Only band cells ?
    const MInt setId = 0;
    for(MInt cellId = 0; cellId < a_noLsCells(); cellId++) {
      // TODO labels:COUPLER Use id conversion
      a_levelSetFunctionMb(cellId, setId) = a_levelSetFunctionG(cellId, setId);
    }

    evaluateContour();
  }

  // Post LB
  else if(recipeStep == 1) {
    if(!lbSolver().domainId()) {
      std::cout << "Post LB couple ..." << std::endl;
    }

    RECORD_TIMER_START(m_timers[Timers::ComputeBoundaryValues]);
    lbBndCnd().extrapolateVelocitiesMb();

    // Exchange velocity between boundary cells
    auto setVel = static_cast<MFloat& (CouplingLb::MbCellCollector::*)(const MInt, const MInt)>(
        &CouplingLb::MbCellCollector::velocity);
    auto getVel = static_cast<MFloat (CouplingLb::MbCellCollector::*)(const MInt, const MInt) const>(
        &CouplingLb::MbCellCollector::velocity);
    auto setData = std::bind(setVel, &a_mbCell(), std::placeholders::_1, std::placeholders::_2);
    auto getData = std::bind(getVel, &a_mbCell(), std::placeholders::_1, std::placeholders::_2);

    auto cellMapping = std::bind(&CouplingLb::a_boundaryCellMb, this, std::placeholders::_1, std::placeholders::_2);

    lbSolver().exchangeSparseLeafValues(getData, setData, nDim, cellMapping);

    // Check all vels in mbcellcollector for nan
    for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
      const MFloat u = a_mbCell().velocity(mbCell, 0);

      if(std::isnan(u)) {
        mTerm(1, "Should not happen");
      }
    }
    RECORD_TIMER_STOP(m_timers[Timers::ComputeBoundaryValues]);

    RECORD_TIMER_START(m_timers[Timers::CreateComm]);
    lbBndCnd().createMBComm();
    RECORD_TIMER_STOP(m_timers[Timers::CreateComm]);

    RECORD_TIMER_START(m_timers[Timers::ApplyBoundaryCondition]);
    lbBndCnd().postCouple();
    RECORD_TIMER_STOP(m_timers[Timers::ApplyBoundaryCondition]);
  }

  RECORD_TIMER_STOP(m_timers[Timers::PostCouple]);
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::setExtensionVelocity() {
  // Revert the algorithm:
  // Try to find a reasonable velocity value for every G0 cell...
  for(MInt id = 0; id < a_noG0Cells(0); id++) {
    const MInt cellId = a_G0CellId(id, 0);

    if(lsSolver().a_isHalo(cellId)) {
      continue;
    }

    const MInt mbCell = a_boundaryCellMb(cellId);

    if(mbCell > -1) {
      // Boundary cell found

      for(MInt n = 0; n < nDim; n++) {
        a_extensionVelocityG(cellId, n, 0) = a_mbCell().velocity(mbCell, n) * conversionLbLs.velocity;
      }

      if(std::isnan(a_extensionVelocityG(cellId, 0, 0))) {
        std::cout << "Velocity nan in boundary cell..." << std::endl;
      }

    } else {
      // No boundary cell found

      // Find the nearest boundary cell by maximizing e * n
      MFloat minBndDist = std::numeric_limits<MFloat>::min();
      MInt minDistDir = -1;
      for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
        MFloat bndDist = 0;
        for(MInt n = 0; n < nDim; n++) {
          bndDist += a_normalVectorG(cellId, n, 0) * LbLatticeDescriptor<nDim, nDist>::ppdfDir(dist, n);
        }
        IF_CONSTEXPR(nDim == 2) {
          if(dist >= 4) {
            bndDist /= SQRT2;
          }
        }
        else {
          if(dist >= 6 && dist < 18) {
            bndDist /= SQRT2;
          } else if(dist >= 18) {
            bndDist /= SQRT3;
          }
        }
        bndDist *= a_levelSetFunctionG(cellId, 0);

        if(minBndDist < bndDist) {
          minBndDist = bndDist;
          minDistDir = dist;
        }
      }
      const MInt neighbor = lbSolver().c_neighborId(cellId, minDistDir);

      const MInt nextMbCell = a_boundaryCellMb(neighbor);

      if(nextMbCell > -1) {
        // Found the next possible boundary cell
        for(MInt n = 0; n < nDim; n++) {
          a_extensionVelocityG(cellId, n, 0) = a_mbCell().velocity(nextMbCell, n) * conversionLbLs.velocity;
        }

        if(std::isnan(a_extensionVelocityG(cellId, 0, 0))) {
          std::cout << "Velocity nan in boundary cell..." << std::endl;
        }

      } else {
        for(MInt n = 0; n < nDim; n++) {
          a_extensionVelocityG(cellId, n, 0) = 0.0;
        }
        std::cout << "G0 cell " << lbSolver().c_globalId(cellId) << " neighbor " << lbSolver().c_globalId(neighbor)
                  << std::endl;

        for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
          MFloat bndDist = 0;
            for(MInt n = 0; n < nDim; n++) {
              bndDist += a_normalVectorG(cellId, n, 0) * LbLatticeDescriptor<nDim, nDist>::ppdfDir(dist, n);
            }
            IF_CONSTEXPR(nDim == 2) {
              if(dist >= 4) {
                bndDist /= SQRT2;
              }
            }
            else {
              if(dist >= 6 && dist < 18) {
                bndDist /= SQRT2;
              } else if(dist >= 18) {
                bndDist /= SQRT3;
              }
            }
          bndDist *= a_levelSetFunctionG(cellId, 0);
          std::cout << "Dist " << dist << " bndDist " << bndDist << " normal_y " << a_normalVectorG(cellId, 1, 0)
                    << " ppdfdir " << LbLatticeDescriptor<nDim, nDist>::ppdfDir(dist, 1) << std::endl;
        }

        mTerm(1, "Cant find corresponding LB bnd cell for given G0 cell. Should not happen.");
      }
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::setExtensionVelocityB() {
  // Collector is created after this, such that these velocities are from the time step n-1
  for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
    const MInt cellId = a_mbCell().cellId(mbCell);

    const MFloat uLB = a_mbCell().velocity(mbCell, 0);
    const MFloat vLB = a_mbCell().velocity(mbCell, 1);

    const MFloat uLS = uLB * conversionLbLs.velocity;
    const MFloat vLS = vLB * conversionLbLs.velocity;

    if(!std::isnan(uLS) && !std::isnan(vLS)) {
      a_extensionVelocityG(cellId, 0, 0) = uLS;
      a_extensionVelocityG(cellId, 1, 0) = vLS;


      MInt minDistDir = -1;
      MFloat minWallDistance = F2;
      for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
        if(minWallDistance > a_mbCell().distance(mbCell, dist)) {
          minWallDistance = a_mbCell().distance(mbCell, dist);
          minDistDir = dist;
        }
      }

      ASSERT(minDistDir > -1, "minDistDir not found!");

      const MInt neighbor = lbSolver().c_neighborId(cellId, minDistDir);

      if(neighbor == -1 || lsSolver().a_isHalo(neighbor)) {
        continue;
      }

      // a_extensionVelocityG(cellId, 0, 0) = 0.0;
      // a_extensionVelocityG(cellId, 1, 0) = v;
      a_extensionVelocityG(neighbor, 0, 0) = uLS;
      a_extensionVelocityG(neighbor, 1, 0) = vLS;
    }

    if(!lsSolver().a_isHalo(cellId)) {
      // Set value to LS cell correspondig to the LB interface cell
      if(std::isnan(uLS) || std::isnan(vLS)) {
        std::cout << "attempt to set nan at cell " << cellId << " which is an halo " << lsSolver().a_isHalo(cellId)
                  << " which is in " << (a_levelSetFunctionMb(cellId, 0) > 0) << std::endl;
        mTerm(1, "nan set");
      }

      a_extensionVelocityG(cellId, 0, 0) = uLS;
      a_extensionVelocityG(cellId, 1, 0) = vLS;
    }

    // The LS solver has two layers of interface cells

    // Find direction to the missing cell
    MInt minDistDir = -1;
    MFloat minWallDistance = F2;
    for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
      if(minWallDistance > a_mbCell().distance(mbCell, dist)) {
        minWallDistance = a_mbCell().distance(mbCell, dist);
        minDistDir = dist;
      }
    }


    ASSERT(minDistDir > -1, "minDistDir not found!");

    const MInt neighbor = lbSolver().c_neighborId(cellId, minDistDir);

    const MInt opposite = LbLatticeDescriptor<nDim, nDist>::oppositeDist(minDistDir);
    const MInt oppositeNeighbor = lbSolver().c_neighborId(cellId, opposite);

    if(neighbor == -1 || lsSolver().a_isHalo(neighbor) || lsSolver().a_isHalo(oppositeNeighbor)) {
      continue;
    }

    if(std::isnan(uLS) || std::isnan(vLS)) {
      std::cout << "attempt to set nan at neighbor " << neighbor << std::endl;
      std::cout << "cell itself is halo " << lsSolver().a_isHalo(cellId) << std::endl;
      std::cout << "extrap neighbor is halo " << lsSolver().a_isHalo(oppositeNeighbor) << std::endl;
      mTerm(1, "nan set");
    }

    a_extensionVelocityG(neighbor, 0, 0) = uLS;
    a_extensionVelocityG(neighbor, 1, 0) = vLS;
  }
}


// FIXME labels:COUPLER Only valid for congruent grids
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::interpolateCurvature(MFloatScratchSpace& curvature) {
  for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
    const MInt pCellId = a_mbCell().cellId(mbCell);

    // Find minWallDistDir
    MInt minDistDir = -1;
    MFloat minWallDistance = F2;
    for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
      if(minWallDistance > a_mbCell().distance(mbCell, dist)) {
        minWallDistance = a_mbCell().distance(mbCell, dist);
        minDistDir = dist;
      }
    }

    ASSERT(minDistDir > -1, "minDistDir not found!");

    const MFloat q = lbBndCnd().getDistanceMb(pCellId, mbCell, minDistDir);

    const MInt neighbor = lbSolver().c_neighborId(pCellId, minDistDir);

    if(neighbor == -1) {
      continue;
    }

    // const MFloat lin = a_curvatureG(pCellId, 0) * (1-q) + a_curvatureG(neighbor, 0) * q;

    /*const MFloat rad1 = maia::math::sgn(a_curvatureG(pCellId, 0))
                        * sqrt(1/abs(a_curvatureG(pCellId, 0)));
    const MFloat rad2 = maia::math::sgn(a_curvatureG(neighbor,0))
                        * sqrt(1/abs(a_curvatureG(neighbor,0)));*/

    const MFloat rad1 = 1 / a_curvatureG(pCellId, 0);
    const MFloat rad2 = 1 / a_curvatureG(neighbor, 0);

    const MFloat rad = rad1 * (1 - q) + rad2 * q;

    // const MFloat lin = maia::math::sgn(rad)*1/POW2(rad);
    const MFloat lin = 1 / rad;

    curvature[mbCell] = lin;
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::interpolateNormal() {
  for(MInt mbCell = 0; mbCell < a_mbCell().size(); mbCell++) {
    const MInt pCellId = a_mbCell().cellId(mbCell);

    // Find minWallDistDir
    MInt minDistDir = -1;
    MFloat minWallDistance = F2;
    for(MInt dist = 0; dist < a_noDistributions() - 1; dist++) {
      if(minWallDistance > a_mbCell().distance(mbCell, dist)) {
        minWallDistance = a_mbCell().distance(mbCell, dist);
        minDistDir = dist;
      }
    }

    ASSERT(minDistDir > -1, "minDistDir not found!");

    const MFloat q = lbBndCnd().getDistanceMb(pCellId, mbCell, minDistDir);

    ASSERT(minDistDir > -1, "No min dist found.");

    const MInt neighbor = lbSolver().c_neighborId(pCellId, minDistDir);

    if(neighbor == -1) {
      continue;
    }

    for(MInt n = 0; n < nDim; n++) {
      const MFloat lin =
          lsSolver().a_normalVectorG(pCellId, n, 0) * (1 - q) + lsSolver().a_normalVectorG(neighbor, n, 0) * q;

      a_mbCell().normal(mbCell, n) = lin;
    }
  }
}


/** \brief Updates the member-variables in the geometry-intersection class
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::updateGeometry() {
  TRACE();

  lbSolver().m_geometryIntersection->m_noEmbeddedBodies = lbSolver().m_noEmbeddedBodies;
  lbSolver().m_geometryIntersection->m_noLevelSetsUsedForMb = lbSolver().m_noLevelSetsUsedForMb;
  lbSolver().m_geometryIntersection->m_bodyToSetTable = lsSolver().m_bodyToSetTable;
  lbSolver().m_geometryIntersection->m_setToBodiesTable = lsSolver().m_setToBodiesTable;
  lbSolver().m_geometryIntersection->m_noBodiesInSet = lsSolver().m_noBodiesInSet;
}

// TODO labels:COUPLER Move up
template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::refillEmergedCells() {
  TRACE();

  for(MInt i = 0; i < a_noCells(); i++) {
    if(lbSolver().a_isHalo(i)) {
      continue;
    }

    // Regular fluid cell
    if(a_isActive(i) && a_wasActive(i)) {
      continue;
    }

    // Regular outside cell
    // TODO labels:COUPLER Access via PV
    if(!a_isActive(i)) {
      for(MInt n = 0; n < nDim; n++) {
        lbSolver().a_variable(i, n) = 0.0;
      }
      lbSolver().a_variable(i, nDim) = 1.0;
    }

    if(a_isActive(i) && !a_wasActive(i)) {
      lbBndCnd().refillEmergedCell(i);
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LsLbSurface<nDim, nDist, SysEqn>::evaluateContour() {
  MFloat COG[nDim]{0.0};
  MFloat mass = 0.0;

  for(MInt cellId = 0; cellId < a_noLsCells(); cellId++) {
    if(lsSolver().a_isHalo(cellId)) {
      continue;
    }

    if(a_levelSetFunctionG(cellId, 0) > 0) {
      continue;
    }

    const MFloat dm = 1.0;

    for(MInt n = 0; n < nDim; n++) {
      COG[n] += dm * lsSolver().c_coordinate(cellId, n);
    }
    mass += dm;
  }

  if(lsSolver().noDomains() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, COG, 2, MPI_DOUBLE, MPI_SUM, lsSolver().mpiComm(), AT_, "COG", "COG");
    MPI_Allreduce(MPI_IN_PLACE, &mass, 1, MPI_DOUBLE, MPI_SUM, lsSolver().mpiComm(), AT_, "mass", "mass");
  }

  for(MInt n = 0; n < nDim; n++) {
    COG[n] /= mass;
  }

  const MFloat dx = a_cellLengthAtLevel(lbSolver().maxLevel());

  if(!lsSolver().domainId()) {
    std::cout << "Bubble COG " << COG[0] << " " << COG[1] << " mass " << mass << " hydrostatic pressure rho "
              << 1.0 - 3.0 * COG[1] / dx * m_gravity << std::endl;
  }

  if(!lsSolver().domainId()) {
    FILE* log;
    log = fopen("bubble.log", "a+");
    fprintf(log, "%d ", globalTimeStep);
    fprintf(log, "%f ", COG[1]);
    fprintf(log, "\n");
    fclose(log);
  }
}

// Explicit instantiations
template class LsLbSurface<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class LsLbSurface<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LsLbSurface<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
