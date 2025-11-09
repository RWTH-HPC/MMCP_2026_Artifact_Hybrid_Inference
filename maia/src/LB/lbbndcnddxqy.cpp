// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbbndcnddxqy.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_set>
#include "COMM/mpioverride.h"
#include "GEOM/geometryelement.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "UTIL/parallelfor.h"
#include "lbfunctions.h"
#include "lbsyseqn.h"

using namespace lbconstants;

/**
 * \brief Contructor for the nDim/nDist specific boundary conditions
 *
 * \param solver LB solver to be manipulated by the boundary condition
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbBndCndDxQy<nDim, nDist, SysEqn>::LbBndCndDxQy(LbSolver<nDim>* solver) : LbBndCnd<nDim>(solver) {
  TRACE();

  NEW_TIMER_GROUP(tgrp_BC3D, "Boundary Condition 3D");
  NEW_TIMER(t_BC3DAll, "complete BC setup", tgrp_BC3D);
  RECORD_TIMER_START(t_BC3DAll);

  m_log << std::endl;
  m_log << "#########################################################################################################"
           "#############"
        << std::endl;
  m_log << "##                                             " << nDim
        << "D Boundary Conditions                                    "
           "           ##"
        << std::endl;
  m_log << "#########################################################################################################"
           "#############"
        << std::endl
        << std::endl;

  NEW_SUB_TIMER(t_initMembers, "init members", t_BC3DAll);
  RECORD_TIMER_START(t_initMembers);

  m_solver = static_cast<LbSolverDxQy<nDim, nDist, SysEqn>*>(solver);

  // for prescibed rho at inflow / outflow boundaries
  m_rho1 = m_solver->m_rho1;
  m_rho2 = m_solver->m_rho2;

  /*! \page propertyPage1
    \section lbControlInflow
    <code>MInt LbBndCndDxQy::m_lbControlInflow</code>\n
    default = <code>0</code>\n\n
    Controls, if a non-rectangular profile is prescribed at an in/ouflow BC.
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_lbControlInflow = 0;
  m_lbControlInflow = Context::getSolverProperty<MInt>("lbControlInflow", m_solverId, AT_, &m_lbControlInflow);

  /*! \page propertyPage1
    \section lbZeroInflowVelocity
    <code>MFloat LbBndCndDxQy::m_lbZeroInflowVelocity</code>\n
    default = <code>1.0</code>\n\n
    Sets the inflow velocity to 0.0
    <ul>
    <li><code>0.0</code> </li>
    <li><code>1.0</code> </li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_zeroInflowVelocity = 1.0;
  m_zeroInflowVelocity =
      Context::getSolverProperty<MFloat>("zeroInflowVelocity", m_solverId, AT_, &m_zeroInflowVelocity);

  RECORD_TIMER_STOP(t_initMembers);

  if(m_solver->isActive()) {
    NEW_SUB_TIMER(t_calcWDist, "calculate wall distances", t_BC3DAll);
    g_tc_geometry.push_back(std::pair<MString, MInt>("calculate wall distances", t_calcWDist));
    RECORD_TIMER_START(t_calcWDist);
    if(m_calcSublayerDist) {
      mAlloc(m_distIntersectionElementId, m_solver->a_noCells(), nDist, "distIntersectionElementId", -1, AT_);
    }
    calculateWallDistances(); // needed for bounce back with inclined walls
    RECORD_TIMER_STOP(t_calcWDist);

    NEW_SUB_TIMER(t_parInflow, "apply parabolic inflow", t_BC3DAll);
    g_tc_geometry.push_back(std::pair<MString, MInt>("apply parabolic inflow", t_parInflow));
    RECORD_TIMER_START(t_parInflow);

    if(m_lbControlInflow) {
      switch(m_lbControlInflow) {
        case 1:
          parabolicInflow(1);
          break;
        case 2:
          parabolicInflow(2);
          break;
        default:
          m_log << "LbBndCndDxQy::() unknown index for parabolic inflow " << std::endl;
          std::stringstream errorMessage;
          errorMessage << "LbBndCndDxQy::() unknown index for parabolic inflow";
          DEBUG("lbbndcnddxqy::parabolicInflow() unknown index for parabolic inflow" << std::endl,
                MAIA_DEBUG_ASSERTION);
          TERMM(1, errorMessage.str());
      }
    }
    RECORD_TIMER_STOP(t_parInflow);
  }

  m_maxDeltaRho = 1.0 - m_rho1;
  m_rhoLast = m_rho1;
  m_lRho = 0.0;

  m_currentTimeStep = globalTimeStep;

  /*! \page propertyPage1
    \section bounceBackSchemeMb
    <code>MFloat LbBndCndDxQy::bounceBackSchemeMb</code>\n
    default = <code>BOUZIDI_QUADRATIC</code>\n\n
    Controls the used bounce back scheme for moving boundaries.
    <ul>
    <li><code>BOUZIDI_LINEAR</code> </li>
    <li><code>BOUZIDI_QUADRATIC</code> </li>
    <li><code>YU_QUADRATIC</code> </li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  const MString bounceBackSchemeMb_default = "BOUZIDI_QUADRATIC";
  bounceBackSchemeMb = string2enum(
      Context::getSolverProperty<MString>("bounceBackSchemeMb", m_solverId, AT_, &bounceBackSchemeMb_default));

  switch(bounceBackSchemeMb) {
    case(BOUZIDI_LINEAR):
#ifdef WAR_NVHPC_PSTL
      m_bounceBackFunctionMbCase = 0;
#else
      bounceBackFunctionMb = &LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_lin;
#endif
      break;
    case(BOUZIDI_QUADRATIC):
#ifdef WAR_NVHPC_PSTL
      m_bounceBackFunctionMbCase = 1;
#else
      bounceBackFunctionMb = &LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_qua;
#endif
      break;
    case(YU_QUADRATIC):
#ifdef WAR_NVHPC_PSTL
      m_bounceBackFunctionMbCase = 2;
#else
      bounceBackFunctionMb = &LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Yu_qua;
#endif
      break;
    default:
      TERMM(1, "Bounceback Scheme Mb not set!");
      break;
  }

  /*! \page propertyPage1
    \section lbRefillMethodOrder
    <code>MInt LbBndCndDxQy::m_lbRefillMethodOrder</code>\n
    default = <code>2</code>\n\n
    Determines order of normal extrapolation for refill method of emerged cells (moving boundaries)
    <ul>
    <li><code>1</code> (linear)</li>
    <li><code>2</code> (quadratic)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_refillMethodOrder = 2;
  m_refillMethodOrder = Context::getSolverProperty<MInt>("refillMethodOrder", m_solverId, AT_, &m_refillMethodOrder);

  RECORD_TIMER_STOP(t_BC3DAll);

  if(m_solver->m_geometry->m_debugParGeom) {
    DISPLAY_ALL_GROUP_TIMERS(tgrp_BC3D);
    m_log << std::endl;
  }

  m_log << std::endl << std::endl;
}

/**
 * \brief Destructor for the nDim/nDist specific boundary condition
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbBndCndDxQy<nDim, nDist, SysEqn>::~LbBndCndDxQy() {
  if(m_calcSublayerDist) mDeallocate(m_distIntersectionElementId);
  if(m_oldWallTemp != nullptr) mDeallocate(m_oldWallTemp);
  if(m_mucousDist != nullptr) mDeallocate(m_mucousDist);
  if(m_fluidDist != nullptr) mDeallocate(m_fluidDist);
}

/** \LBBC{secLBBC_bc0, bc0, 0}
 * Dummy function for exception handling<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc0(MInt)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc0(MInt NotUsed(index)) {}

/** \brief Writes the rim of a boundary as VTP to file
 *
 * \author Andreas Lintermann
 * \date 21.01.2013
 *
 * The file rim_domainId.vtp, where domainId is the current MPI-rank
 * is written per domain participating in a segment with id segmentId.
 *
 * \param[in] vertices the vertices as pointer to a pointer
 * \param[in] num the number of vertices
 * \param[in] segmentId the id of the segment
 *
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::writeBoundaryVTK(MFloat** vertices, MInt num, MInt segmentId) {
  TRACE();

  std::ofstream ofl;
  std::stringstream fname;
  fname << "rim_s" << segmentId << "_" << m_solver->domainId() << ".vtp";
  ofl.open(fname.str().c_str());

  //================== VTKFile =================
  ofl << "<?xml version=\"1.0\"?>" << std::endl;
  ofl << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  ofl << "<PolyData>" << std::endl;
  ofl << "<Piece NumberOfPoints=\"" << num << "\" NumberOfVerts=\"" << num << "\" NumberOfLines=\"" << 1 << "\">"
      << std::endl;

  ofl << "<PointData TCoords=\"Texture Coordinates\">" << std::endl;
  ofl << "<DataArray type=\"Float32\" Name=\"Array_Position\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  for(MInt i = 0; i < num; i++)
    ofl << i << " ";
  ofl << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "</PointData>" << std::endl;

  //================== Points =================
  ofl << "<Points>" << std::endl;
  ofl << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  for(MInt i = 0; i < num; i++) {
    for(MInt d = 0; d < nDim; d++) {
      ofl << vertices[i][d] << " ";
    }
  }
  ofl << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "</Points>" << std::endl;
  //================== /Points =================

  /*
  //================== Verts ===================
  ofl << "<Verts>" << std::endl;
  ofl << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  for(MInt i = 0; i < num; i++)
    ofl << i << " ";
  ofl << "0" << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  for(MInt i = 1; i <= num; i++)
    ofl << i << " ";
  ofl << num+1 << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "</Verts>" << std::endl;
  //================== /Verts ==================
  */
  //================== Lines ===================
  ofl << "<Lines>" << std::endl;
  ofl << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  for(MInt i = 0; i < num; i++)
    ofl << i << " ";
  ofl << "0" << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  // for(MInt i = 1; i <= num; i++)
  //  ofl << i << " ";
  ofl << num + 1 << std::endl;
  ofl << "</DataArray>" << std::endl;
  ofl << "</Lines>" << std::endl;
  //================== /Lines ==================

  ofl << "</Piece>" << std::endl;
  ofl << "</PolyData>" << std::endl;
  ofl << "</VTKFile>" << std::endl;
  ofl.close();
}

/** \brief calculates the intersection point between the boundary wall and the
 *        trajectories of the distribution functions.
 *
 * This function is essential for the implementation of an improved bounce
 * back scheme which can deal with arbitrary boundaries (i.e. the wall
 * needn't be located halfway between the cell centers as with simple bounce
 * back rule)
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateWallDistances() {
  TRACE();

  m_log << "  + Calculating wall distances..." << std::endl;
  m_log << "    - method:          " << m_interpolationDistMethod << std::endl;
  m_log << "    - grid cut method: " << m_gridCutTest << std::endl;

  // init the volumes
  for(auto& i : m_bndCells) {
    i.m_isFluid = false;
  }

  if(m_interpolationDistMethod == "sphere" || m_interpolationDistMethod == "circle") {
    const MInt nDim_ = (m_interpolationDistMethod == "sphere") ? 3 : 2;
    for(MInt i = 0; i < (MInt)m_bndCells.size(); i++) {
      const MInt currentId = m_bndCells[i].m_cellId;

      const MFloat cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(currentId) + 1);
      const MFloat cellRadiusSq =
          std::inner_product(&m_solver->a_coordinate(currentId, 0), &m_solver->a_coordinate(currentId, nDim_),
                             &m_solver->a_coordinate(currentId, 0), .0);

      //--Determine whether cell center is inside of fluid or not---------------
      // const MFloat radius = 0.5;
      const MFloat radiusSq = 0.25;
      m_bndCells[i].m_isFluid = (cellRadiusSq > radiusSq);

      //--Calculate distance for each distribution------------------------------
      // init distances
      m_bndCells[i].m_distances.clear();
      m_bndCells[i].m_distances.resize(IPOW3[nDim] - 1, 0.0);
      for(MInt dist = 0; dist < nDist - 1; dist++) {
        const MFloat dxSq = Ld::distType(dist) * 4.0 * cellHalfLength * cellHalfLength;

        // Set distance to a large value...
        m_bndCells[i].m_distances[dist] = std::numeric_limits<MFloat>::max();

        // solve quadratic equation for q: R^2 = (x_cell + q*c_dist)^2
        // if cell belongs to the wall q>1/2
        // otherwise q<=1/2

        MFloat qValue = 0.0; // normalized wall distance: q = distance between cell center and wall / dx
        for(MInt j = 0; j < nDim_; j++) {
          qValue -=
              m_solver->a_coordinate(currentId, j) * ((MFloat)Ld::idFld(dist, j) - 1.0) * 2.0 * cellHalfLength / dxSq;
        }

        const MFloat tmpQValue = sqrt(qValue * qValue + (radiusSq - cellRadiusSq) / dxSq);

        if(qValue + tmpQValue >= 0.0 && fabs(qValue + tmpQValue) <= 0.5) {
          qValue += tmpQValue;
          m_bndCells[i].m_distances[dist] = qValue;
        }
        if(qValue - tmpQValue >= 0.0 && fabs(qValue - tmpQValue) <= 0.5) {
          qValue -= tmpQValue;
          m_bndCells[i].m_distances[dist] = qValue;
        }
      }
    }
  }

  else if(m_interpolationDistMethod == "perpOp" || m_interpolationDistMethod == "cylinder") {
    // TODO labels:LB This method differs from the 2D and all other 3D variants, such
    // that all inOutSegments boundaries are completely set to m_isFluid = false.
    // In the other cases this is calculated, too, which may lead to confusion!

    // find ids of non-periodic and non-inOutSegments (until now this is wall)
    std::vector<MInt> wallIds;
    for(MInt i = 0; i < (MInt)(m_bndCndSegIds.size()); i++) {
      MBool is_periodic = false;
      if(m_noPeriodicSegments != 0)
        for(MInt j = 0; j < m_noPeriodicSegments; j++)
          if(m_bndCndSegIds[i] == m_periodicSegmentsIds[j]) {
            is_periodic = true;
            break;
          }

      MBool is_inout = false;
      for(MInt j = 0; j < m_noInOutSegments; j++)
        if(m_bndCndSegIds[i] == m_inOutSegmentsIds[j]) {
          is_inout = true;
          break;
        }

      if(!is_periodic & !is_inout) wallIds.push_back(i);
    }

    if(wallIds.size() > 0) {
      std::unordered_set<MInt> reconsider;

      for(auto wallId : wallIds) {
        m_log << "    - BC " << m_bndCndIds[wallId] << std::endl;
        m_log << "      * number or cells: " << (m_bndCndOffsets[wallId + 1] - m_bndCndOffsets[wallId]) << std::endl;
        m_log << "      * segment id:      " << m_bndCndSegIds[wallId] << std::endl;

        for(MInt i = m_bndCndOffsets[wallId]; i < m_bndCndOffsets[wallId + 1]; i++) {
          const MInt currentId = m_bndCells[i].m_cellId;

          //--Calculate distance for each distribution--------------------------
          const MFloat cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(currentId) + 1);

          MFloat d[nDim], e[nDim];
          MFloat target[2 * nDim];
          for(MInt j = 0; j < nDim; j++) {
            target[j] = m_solver->a_coordinate(currentId, j) - cellHalfLength;
            target[j + nDim] = m_solver->a_coordinate(currentId, j) + cellHalfLength;

            // Take cell center as first point of discrete velocity trajectories
            d[j] = m_solver->a_coordinate(currentId, j);
          }

          // get all triangles in currentCell (target)
          std::vector<MInt> nodeList;
          m_solver->m_geometry->getIntersectionElements(target, nodeList, cellHalfLength,
                                                        &m_solver->a_coordinate(currentId, 0));

          // init distances
          m_bndCells[i].m_distances.clear();
          m_bndCells[i].m_distances.resize(IPOW3[nDim] - 1, 0.0);

          // run over all distributions and calculate distances
          for(MInt dist = 0; dist < IPOW3[nDim] - 1; dist++) {
            const MFloat dx = F2 * SQRT[Ld::distType(dist)] * cellHalfLength;

            // Set distance to a large value...
            m_bndCells[i].m_distances[dist] = std::numeric_limits<MFloat>::max();

            // Construct second point of discrete velocity trajectories
            for(MInt j = 0; j < nDim; j++)
              e[j] = d[j] + ((MFloat)Ld::idFld(dist, j) - 1.0) * cellHalfLength;

            MFloat targetsmall[2 * nDim];
            for(MInt k = 0; k < nDim; k++) {
              targetsmall[k] = d[k];
              targetsmall[k + nDim] = e[k];
            }

            MFloat min_dist;
            MBool has_cut = m_solver->m_geometry->getClosestLineIntersectionLength(m_bndCndIds[wallId], nodeList,
                                                                                   targetsmall, &min_dist);
            if(has_cut)
              m_bndCells[i].m_distances[dist] =
                  min_dist / dx; // normalized wall distance: q = distance between cell center and wall / dx
            else
              m_bndCells[i].m_distances[dist] = 1.0;
          }

          //--Determine whether cell center is inside of fluid or not-----------
          // In the following the state (fluid/solid) of each boundary cell is
          // detected based on the state of the neighboring cells being located
          // in direction of a cut.
          m_bndCells[i].m_isFluid = false;
          MBool hasBndryCutNghbr = false;
          MBool hasNoCutNghbr = false;
          MBool hasFluidCutNghbr = false;
          for(MInt dist = 0; dist < IPOW3[nDim] - 1; dist++) {
            if(m_bndCells[i].m_distances[dist] < 1.0) {
              if(m_solver->a_hasNeighbor(currentId, dist)) {
                const MInt nghId = m_solver->c_neighborId(currentId, dist);
                if(m_solver->a_isBndryCell(nghId)) {
                  hasBndryCutNghbr = true;
                } else {
                  // In case of any non-boundary cut neighbor, the state is
                  // directly definable (it is non-fluid) -> break
                  hasFluidCutNghbr = true;
                  break;
                }
              } else {
                // Having no cut neighbor does not directly specify the state.
                // In case of intersecting boundaries a neighbor can be missing
                // for fluid as well as for non-fluid state.
                hasNoCutNghbr = true;
              }
            }
          }
          if(hasFluidCutNghbr)
            m_bndCells[i].m_isFluid = false;
          else if(hasBndryCutNghbr)
            reconsider.insert(i);
          else if(hasNoCutNghbr)
            m_bndCells[i].m_isFluid = true;
          // else :
          // In this case something must have gone wrong previously. A cut with
          // an existing nghbr, which is neither a boundary nor simple fluid
          // cells is new to my knowledge.
        }
      }
      // Boundary cells that only have cut-neighbors which are boundary cells,
      // too, need to be reconsidered. This might has to be done multiple times
      // until the state of at least one neighbor is known for each reconsidered
      // cell.
      MUint reconsiderSize = 0;
      while(reconsiderSize != reconsider.size()) {
        reconsiderSize = reconsider.size();
        m_log << "    - number of cells to reconsider: " << reconsiderSize << std::endl;
        auto r = reconsider.begin();
        while(r != reconsider.end()) {
          const MInt cellId = m_bndCells[*r].m_cellId;
          MBool chk = false;
          for(MInt dist = 0; dist < IPOW3[nDim] - 1; dist++) {
            if(m_bndCells[*r].m_distances[dist] < 1.0 && m_solver->a_hasNeighbor(cellId, dist)) {
              const MInt nghId = m_solver->c_neighborId(cellId, dist);
              if(m_solver->a_isBndryCell(nghId)) {
                const MInt bndIdNeigh = m_solver->a_bndId(nghId);
                if(reconsider.find(bndIdNeigh) == reconsider.end()) {
                  // only one cut between cellId and nghId
                  if(approx(m_bndCells[bndIdNeigh].m_distances[Ld::oppositeDist(dist)], 1.0, MFloatEps)) {
                    m_bndCells[*r].m_isFluid = !m_bndCells[bndIdNeigh].m_isFluid; // cellId has opposite state of nghId
                    chk = true;
                    break;
                  }
                }
              }
            }
          }
          if(chk)
            reconsider.erase(r++);
          else
            r++;
        }
      }
    }

    // analytical correction of wall distance for unit cylinder aligned with origin z-axis
    if(m_interpolationDistMethod == "cylinder") {
      // get correct ID
      const MInt cylinderBcId = Context::getSolverProperty<MInt>("cylinderBcId", m_solverId, AT_);
      MInt cylinderId = -1;
      for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
        if(m_bndCndIds[i] == cylinderBcId) {
          cylinderId = i;
          break;
        }
      }
      if(cylinderId != -1) {
        const MFloat radius = 0.5;
        const MFloat radiusSqr = radius * radius;

        for(MInt i = m_bndCndOffsets[cylinderId]; i < m_bndCndOffsets[cylinderId + 1]; i++) {
          const MInt currentId = m_bndCells[i].m_cellId;

          const MFloat cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(currentId) + 1);

          const MFloat cellRadiusSq = m_solver->a_coordinate(currentId, 0) * m_solver->a_coordinate(currentId, 0)
                                      + m_solver->a_coordinate(currentId, 1) * m_solver->a_coordinate(currentId, 1);

          //--Determine whether cell center is inside of fluid or not---------------
          m_bndCells[i].m_isFluid = cellRadiusSq > radiusSqr;

          //--Calculate distance for each distribution------------------------------
          // init distances
          m_bndCells[i].m_distances.clear();
          m_bndCells[i].m_distances.resize(IPOW3[nDim] - 1, 0.0);
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            MFloat dx;
            if(dist < 6)
              dx = F2 * cellHalfLength;
            else if(dist < 18)
              dx = SQRT2 * F2 * cellHalfLength;
            else
              dx = SQRT3 * F2 * cellHalfLength;
            const MFloat dxSqr = dx * dx;

            // Set distance to a large value...
            m_bndCells[i].m_distances[dist] = std::numeric_limits<MFloat>::max();

            // solve quadratic equation for q: R^2 = (x_cell + q*c_dist)^2
            // if cell belongs to the wall q>1/2
            // otherwise q<=1/2

            MFloat qValue = 0.0;     // normalized wall distance: q = distance between cell center and wall / dx
            const MInt nDim_cyl = 2; // only consider x and y coordinates
            for(MInt j = 0; j < nDim_cyl; j++) {
              qValue -= m_solver->a_coordinate(currentId, j) * ((MFloat)Ld::idFld(dist, j) - 1.0) * 2.0 * cellHalfLength
                        / dxSqr;
            }

            const MFloat tmpQValue = sqrt(qValue * qValue - (cellRadiusSq - radiusSqr) / dxSqr);

            if(qValue + tmpQValue >= 0.0 && fabs(qValue + tmpQValue) <= 0.5) {
              m_bndCells[i].m_distances[dist] = qValue + tmpQValue;
            } else if(qValue - tmpQValue >= 0.0 && fabs(qValue - tmpQValue) <= 0.5) {
              m_bndCells[i].m_distances[dist] = qValue - tmpQValue;
            }
          }
        }
      }
    }
  }

  // distance determination for x-oriented pipe
  else if(m_interpolationDistMethod == "pipe") {
  }

  else if(m_interpolationDistMethod == "STD") {
    //--I./II. find cut distance in bounding box of 2*dx---------------------------
    constexpr MFloat relEps = 1e-16;
    const MFloat eps = relEps * m_solver->c_cellLengthAtLevel(m_solver->maxLevel());
    auto getInsideOutsideAndDistances = [&](const MInt i, const MInt index) {
      auto& bndCell = m_bndCells[i];
      const MInt currentId = bndCell.m_cellId;
      //--I. in-/outside check--------------------------------------------------
      bndCell.m_isFluid = !(m_solver->m_geometry->pointIsInside2(&m_solver->a_coordinate(currentId, 0)));
      //--determine distances for each bndry cell's distribution--------------
      // get all triangles that intersect currentId
      const MFloat cellLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(currentId));
      // A fluid cell searches for the closest cut, whereas a solid one searches
      // largest cut distance.
      // This assures that in case of multiple cuts between solid-fluid bndCells
      // both find the same cut closest to the fluid region. If q > 0.5 for the
      // fluid cell and a nghbr solid cell is existing, the cut from the fluid
      // is removed to avoid multiple execution.
      // In case of multiple cuts between solid-solid bndCells (p.e.thin gap)
      // this results in q > 0.5 for both cells.
      const MFloat bbHalfLenght = cellLength;
      MFloat target[2 * nDim]; // bounding box of the currentId
      MFloat d[nDim];          // starting point of each ray
      for(MInt j = 0; j < nDim; j++) {
        target[j] = m_solver->a_coordinate(currentId, j) - bbHalfLenght;
        target[j + nDim] = m_solver->a_coordinate(currentId, j) + bbHalfLenght;
        d[j] = m_solver->a_coordinate(currentId, j);
      }
      std::vector<MInt> nodeList;
      if(m_gridCutTest == "SAT")
        m_solver->m_geometry->getIntersectionElements(target, nodeList, bbHalfLenght,
                                                      &m_solver->a_coordinate(currentId, 0));
      else
        m_solver->m_geometry->getIntersectionElements(target, nodeList);
      // init distances
      bndCell.m_distances.clear();
      bndCell.m_distances.resize(IPOW3[nDim] - 1, 0.0);
      // Loop over each distribution
      for(MInt dist = 0; dist < nDist - 1; dist++) {
        // smallest cut distance:
        MFloat e[nDim]; // end point of current ray
        for(MInt j = 0; j < nDim; j++) {
          e[j] = d[j] + ((MFloat)Ld::idFld(dist, j) - 1.0) * bbHalfLenght;
        }
        MFloat trgDistance =
            (bndCell.m_isFluid) ? std::numeric_limits<MFloat>::max() : std::numeric_limits<MFloat>::lowest();
        MInt trgNodeId = -1;
        for(MInt n = 0; n < (signed)nodeList.size(); n++) {
          MFloat distance = 0.0;
          MFloat** const v = m_solver->m_geometry->elements[nodeList[n]].m_vertices;
          // TODO labels:LB D2Qx: use getLineTriangleIntersectionSimpleDistance
          // alternativ that is working for 2D.
          const MBool hasCut =
              m_solver->m_geometry->getLineTriangleIntersectionSimpleDistance(d, e, v[0], v[1], v[2], &distance);
          if(hasCut) {
            distance = (distance < 0.0) ? 0.0 : distance;
            // With the following it shall be assured, that the closest
            // distance is found. Within a tolerance eps the  bndry Id
            // belonging to the currentId is preferred.
            if(fabs(distance - trgDistance) < eps) {
              if(m_solver->m_geometry->elements[nodeList[n]].m_bndCndId == m_bndCndIds[index]) {
                trgDistance = distance;
                trgNodeId = n;
              }
            } else if((bndCell.m_isFluid && distance < trgDistance) || // Choose smallest distance for fluid cell
                      (!bndCell.m_isFluid && distance > trgDistance)   // Choose largest distance for solid cell
            ) {
              trgDistance = distance;
              trgNodeId = n;
            }
          }
        }
        if(trgNodeId > -1 && m_solver->m_geometry->elements[nodeList[trgNodeId]].m_bndCndId == m_bndCndIds[index]) {
          const MFloat dx = SQRT[Ld::distType(dist)] * cellLength;
          bndCell.m_distances[dist] = trgDistance / dx;
          if(m_calcSublayerDist) {
            m_distIntersectionElementId[currentId][dist] = nodeList[trgNodeId];
          }
        } else {
          bndCell.m_distances[dist] = std::numeric_limits<MFloat>::max();
        }
      }
      //      // << DBG
      //      std::stringstream ss;
      //      constexpr MInt dbgGlobalId = 132150;
      //      //constexpr MInt dbgGlobalId2 = 131454; // in dir 22
      //      constexpr MInt dbgGlobalId2 = 37446; // in dir 21
      //      if( m_solver->c_globalId(currentId) == dbgGlobalId ||
      //          m_solver->c_globalId(currentId) == dbgGlobalId2
      //          ) {
      //        MInt dbgId = currentId;
      //        ss << "DBG: " << m_solver->c_globalId(dbgId) << ", " << bndCell.m_isFluid << ", "
      //           << m_solver->a_isHalo(dbgId) << std::endl;
      //        ss << "DBG: nghbrStates:" << std::endl;
      //        for(MInt z = 0; z < nDist; z++) {
      //          ss << "     " << z << " : ";
      //          ss << bndCell.m_distances[z]  << ", ";
      //          const MBool chk = m_solver->a_hasNeighbor(dbgId, z);
      //          const MInt nghbrId = chk ? m_solver->c_neighborId(dbgId, z) : -1;
      //          ss << m_solver->c_globalId(nghbrId) << ", ";
      //          if(chk)
      //            ss << m_solver->a_isBndryCell(nghbrId) << ", " << m_solver->a_isHalo(nghbrId);
      //          else
      //            ss << "nan"
      //               << " nan";
      //          ss << std::endl;
      //        }
      //        std::cout << ss.str();
      //      }
      //      // DBG >>
    };

    for(MInt index = 0; index < (MInt)m_bndCndSegIds.size(); index++) {
      for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
        getInsideOutsideAndDistances(i, index);
      }
    }
    //--III. check for fluid split cells: appending new m_bndCells--------------
    const MInt oldNumberOfBndCells = m_bndCells.size();
    MBool needResorting = false;
    for(MInt i = 0; i < oldNumberOfBndCells; i++) {
      if(m_bndCells[i].m_isFluid) {
        const MInt currentId = m_bndCells[i].m_cellId;
        for(MInt j = 0; j < nDist - 1; j++) {
          if(m_bndCells[i].m_distances[j] > 0.5) continue;
          if(!m_solver->a_hasNeighbor(currentId, j)) continue;
          const MInt nghbrId = m_solver->c_neighborId(currentId, j);
          // TODO labels:LB miro: halo cell is inactive
          if(m_solver->a_isBndryCell(nghbrId) || !m_solver->a_isActive(nghbrId) || m_solver->a_isHalo(nghbrId))
            continue;
          // add new boundary cell
          needResorting = true;
          m_bndCells.emplace_back();
          auto& newBndCell = m_bndCells.back();
          newBndCell.m_cellId = nghbrId;
          newBndCell.m_isFluid = true;
          newBndCell.m_segmentId.push_back(m_bndCells[i].m_segmentId[0]);
          newBndCell.m_bndCndId.push_back(m_bndCells[i].m_bndCndId[0]);
          m_solver->a_isBndryCell(nghbrId) = true;
          // Determine in-/outside and cuts for new cells (All cuts will be
          // q > 0.5, but at least some will be q < 1.0 )
          const MInt index = this->m_mapBndCndSegId2Index[m_bndCells[i].m_segmentId[0]];
          const MInt newBndCellId = m_bndCells.size() - 1;
          getInsideOutsideAndDistances(newBndCellId, index);
          // Remove cuts if neighbor is solid boundary cell to avoid double execution
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            if(newBndCell.m_distances[dist] < 1.0) {
              if(m_solver->a_hasNeighbor(nghbrId, dist)) {
                const MInt nghbrId2 = m_solver->c_neighborId(nghbrId, dist);
                if(m_solver->a_isBndryCell(nghbrId2)) {
                  const MInt nghbrBndId2 = m_solver->a_bndId(nghbrId2);
                  if(m_bndCells[nghbrBndId2].m_isFluid == false) {
                    newBndCell.m_distances[dist] = std::numeric_limits<MFloat>::max();
                  }
                }
              }
            }
          }
        }
      }
    }
    //--IV. sort boundary cells if split cells are found------------------------
    if(needResorting) {
      m_log << "    ## redoing sortBoundaryCells" << std::endl;
      this->sortBoundaryCells();
      m_log << "    redoing sortBoundaryCells ## " << std::endl;
    }
    //--V. Correct double execution by solid-fluid------------------------------
    const MInt newNumberOfBndCells = m_bndCells.size();
    maia::parallelFor(0, newNumberOfBndCells, [&](MInt i) {
      auto& bndCell = m_bndCells[i];
      const MInt currentId = bndCell.m_cellId;
      if(bndCell.m_isFluid) { // fluid
        for(MInt j = 0; j < nDist - 1; j++) {
          const MInt opp = Ld::oppositeDist(j);
          const MBool cutj = bndCell.m_distances[j] < 1.0;
          const MBool cutopp = bndCell.m_distances[opp] < 1.0;
          // j
          if(cutj) {
            // perform simple bounce back if cut in opposite directions
            if(bndCell.m_distances[j] <= 0.5) {
              if(cutopp) bndCell.m_distances[j] = 0.5;
            }
            // remove cut, if solid neighbor to avoid multiple execution
            else { // 0.5 < q < 1.0
              if(m_solver->a_hasNeighbor(currentId, j)) {
                const MInt nghbrId = m_solver->c_neighborId(currentId, j);
                if(m_solver->a_isBndryCell(nghbrId)) {
                  if(m_bndCells[m_solver->a_bndId(nghbrId)].m_isFluid == false) {
                    bndCell.m_distances[j] = std::numeric_limits<MFloat>::max();
                  }
                }
              }
            }
          }
          // opp
          if(cutopp) {
            // perform simple bounce back if cut in opposite directions
            if(bndCell.m_distances[opp] <= 0.5) {
              if(cutj) bndCell.m_distances[opp] = 0.5;
            }
            // remove cut, if solid neighbor to avoid multiple execution
            else { // 0.5 < q < 1.0
              if(m_solver->a_hasNeighbor(currentId, opp)) {
                const MInt nghbrId = m_solver->c_neighborId(currentId, opp);
                if(m_solver->a_isBndryCell(nghbrId)) {
                  if(m_bndCells[m_solver->a_bndId(nghbrId)].m_isFluid == false) {
                    bndCell.m_distances[opp] = std::numeric_limits<MFloat>::max();
                  }
                }
              }
            }
          }
        }
      } else { // solid
        // remove cuts in direction of other solid boundary cell
        for(MInt j = 0; j < nDist - 1; j++) {
          if(bndCell.m_distances[j] < 0.5) {
            if(m_solver->a_hasNeighbor(currentId, j)) {
              const MInt nghbrId = m_solver->c_neighborId(currentId, j);
              if(m_solver->a_isBndryCell(nghbrId)) {
                if(m_bndCells[m_solver->a_bndId(nghbrId)].m_isFluid == false) {
                  bndCell.m_distances[j] = std::numeric_limits<MFloat>::max();
                }
              }
            }
          }
        }
      }
    });
  }
  // Dump out distance field and inside outside state
  if(m_outputWallDistanceField) {
    std::stringstream fileName;
    fileName << m_solver->restartDir() << "lbDistanceField" << m_solver->getIdentifier() << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, m_solver->mpiComm());
    //--define global attributes
    parallelIo.setAttribute(m_solver->solverId(), "solverId");
    parallelIo.setAttribute(m_solver->gridInputFileName(), "gridFile", "");
    //--define arrays
    const MPI_Offset firstGlobalId = m_solver->domainOffset(m_solver->domainId());
    const MPI_Offset localNoCells = m_solver->noInternalCells();
    parallelIo.setOffset(localNoCells, firstGlobalId);
    addWallDistanceFieldToOutputFile(parallelIo, true, false);
    addWallDistanceFieldToOutputFile(parallelIo, false, true);
    m_outputWallDistanceField = false;
  }
}

// TODO labels:LB D2Qx: Make this unified with 3D
// TODO labels:LB Does not need SysEqn, in fact. Is resolved as soon as 2D 3D are unified.
template <MInt nDim_, MInt nDist_, class SysEqn>
inline void LbBndCndDxQy<nDim_, nDist_, SysEqn>::calculateWallDistances2D() {
  TRACE();

  constexpr MInt nDim = 2;
  constexpr MInt nDist = 9;
  // init the volumes
  for(MInt i = 0; i < (MInt)m_bndCells.size(); i++) {
    m_bndCells[i].m_isFluid = false;
  }

  if(m_interpolationDistMethod == "perpOp") {
    MBool triangleIntersection;
    MFloat target[6] = {0, 0, 0, 2, 2, 4};
    MFloat cellHalfLength = 0.0;
    MFloat currentDistance = 0;
    //    MFloat dxB2;

    MFloat a[2];      // point in plane
    MFloat b[2];      // point in plane
    MFloat c[2];      // start of piercing edge
    MFloat d[2];      // end of piercing edge
    MFloat s2, gamma; // s1; For pierce point calculation
    // MFloat pP[2];  uncomment for piercePoint calculation

    // find ids of non-periodic and non-inOutSegments (until now this is wall)
    std::vector<MInt> wallIds;
    for(MInt i = 0; i < (MInt)(m_bndCndSegIds.size()); i++) {
      MBool is_periodic = false;
      if(m_noPeriodicSegments != 0)
        for(MInt j = 0; j < m_noPeriodicSegments; j++)
          if(m_bndCndSegIds[i] == m_periodicSegmentsIds[j]) {
            is_periodic = true;
            break;
          }

      MBool is_inout = false;
      for(MInt j = 0; j < m_noInOutSegments; j++)
        if(m_bndCndSegIds[i] == m_inOutSegmentsIds[j]) {
          is_inout = true;
          break;
        }

      if(!is_periodic & !is_inout) wallIds.push_back(i);
    }

    for(auto wallId : wallIds) {
      for(MInt i = m_bndCndOffsets[wallId]; i < m_bndCndOffsets[wallId + 1]; i++) {
        const MInt currentId = m_bndCells[i].m_cellId;

        //--Determine whether cell center is inside of fluid or not-------------
        // TODO labels:LB,toenhance dxqy: this does not conform with 3D yet ! Should be adapted
        if(m_solver->m_geometry->pointIsInside(&m_solver->a_coordinate(currentId, 0))) {
          //       m_log << " BndCell " << i << " [" << currentId << "] is inside geometry. " << std::endl;
          m_bndCells[i].m_isFluid = false;
        } else {
          //       m_log << " BndCell " << i << " [" << currentId << "] is outside geometry. " << std::endl;
          m_bndCells[i].m_isFluid = true;
        }

        //--Calculate distance for each distribution---------------------------
        // init distances
        m_bndCells[i].m_distances.clear();
        m_bndCells[i].m_distances.resize(IPOW3[nDim] - 1, 0.0);

        // Define corners of current cell in target
        for(MInt j = 0; j < nDim; j++) {
          cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(currentId) + 1);
          target[j] = m_solver->a_coordinate(currentId, j) - cellHalfLength;
          target[j + nDim] = m_solver->a_coordinate(currentId, j) + cellHalfLength;

          // Take cell center as first point of discrete velocity trajectories
          c[j] = m_solver->a_coordinate(currentId, j);
        }

        // get all triangles in currentCell (target)
        std::vector<MInt> nodeList;
        if(m_gridCutTest == "SAT")
          m_solver->m_geometry->getIntersectionElements(target, nodeList, cellHalfLength,
                                                        &m_solver->a_coordinate(currentId, 0));
        else
          m_solver->m_geometry->getIntersectionElements(target, nodeList);

        triangleIntersection = false;
        for(MInt dist = 0; dist < 8; dist++) {
          // Set distance to a large value...
          m_bndCells[i].m_distances[dist] = m_solver->c_cellLengthAtLevel(0);
          //          if(dist < 4)
          //            dxB2 = cellHalfLength;
          //          else
          //            dxB2 = SQRT2 * cellHalfLength;

          // Construct second point of discrete velocity trajectories
          for(MInt j = 0; j < nDim; j++) {
            d[j] = c[j] + ((MFloat)Ld::idFld(dist, j) - 1.0) * cellHalfLength;
          }
          currentDistance = 1;
          m_bndCells[i].m_distances[dist] = currentDistance;
          // Check for intersection with discrete velocity trajectories
          for(MInt n = 0; n < (signed)nodeList.size(); n++) {
            if(m_solver->m_geometry->elements[nodeList[n]].m_bndCndId != m_bndCndIds[wallId]) continue;
            if(m_solver->m_geometry->edgeTriangleIntersection(m_solver->m_geometry->elements[nodeList[n]].m_vertices[0],
                                                              m_solver->m_geometry->elements[nodeList[n]].m_vertices[1],
                                                              0, c, d)) {
              triangleIntersection = true;
              // Calculate Distance
              for(MInt k = 0; k < nDim; k++) {
                a[k] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[0][k];
                b[k] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[1][k];
                // d and e are already set
              }
              gamma = (b[0] - a[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - a[1]);

              // s1 = ((c[0]-d[0]) * (a[1]-c[1]) - (c[1]-d[1])*(a[0]-c[0])) / gamma;


              s2 = ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / gamma;

              // 1. b Pierce point pP in plane j:
              /*for (MInt k = 0; k < nDim; k++){
            if(s1*s1 < s2*s2)
                pP[k] = c[k] + s2 * ( d[k] - c[k] );
            else
                pP[k] = a[k] + s1 * ( b[k] - a[k] );
                }*/
              // Take only the magnitude !
              if(s2 < 0.0) s2 = -s2;

              if(s2 > 1.0) continue;

              currentDistance = s2 * F1B2;
              // Take shortest distance
              m_bndCells[i].m_distances[dist] = (currentDistance < m_bndCells[i].m_distances[dist])
                                                    ? currentDistance
                                                    : m_bndCells[i].m_distances[dist];

              // see pierce point calculation in geometry.cpp
            }
          }
        }
        if(!triangleIntersection) {
          m_log << " BndCell[" << i << "] has no triangle intersection!" << std::endl;
        }
      }
    }
  }

  else if(m_interpolationDistMethod == "circle") {
    MFloat radius = 0.5;
    MFloat x = 0.0;
    MFloat y = 0.0;
    std::cout << "Prescribe analytic circle with radius: " << radius << ", x:" << x << ", y:" << y << std::endl;

    MFloat cellHalfLength, cellRadiusSq, dxB2, dxSq;
    MInt currentId;

    // Loop over m_bndCells
    for(MInt i = 0; i < (MInt)m_bndCells.size(); i++) {
      currentId = m_bndCells[i].m_cellId;
      cellHalfLength = m_solver->grid().cellLengthAtLevel(m_solver->a_level(currentId) + 1);
      cellRadiusSq = (m_solver->a_coordinate(currentId, 0) - x) * (m_solver->a_coordinate(currentId, 0) - x)
                     + (m_solver->a_coordinate(currentId, 1) - y) * (m_solver->a_coordinate(currentId, 1) - y);

      // Mark cells as in- or outside the wall.
      if(cellRadiusSq <= radius * radius) {
        m_bndCells[i].m_isFluid = false; // cell is not inside fluid
      } else {
        m_bndCells[i].m_isFluid = true; // cell belongs to fluid
      }

      // Prescribe distance to wall for every PPDF
      m_bndCells[i].m_distances.clear();
      m_bndCells[i].m_distances.resize(IPOW3[nDim] - 1, 0.0);
      for(MInt dist = 0; dist < nDist - 1; dist++) {
        // Set distance to a large valule
        m_bndCells[i].m_distances[dist] = m_solver->grid().cellLengthAtLevel(0);

        // Compute squared length and length for each direction
        if(dist < 4) {
          dxSq = 4.0 * cellHalfLength * cellHalfLength;
          dxB2 = cellHalfLength;
        } else {
          dxSq = 2.0 * 2.0 * cellHalfLength * cellHalfLength; // Diagonal directions
          dxB2 = SQRT2 * cellHalfLength;
        }

        // Intersection circle directions delivers quadratic equation:
        // d: lenght cell center -- cut
        // a: x-component of d
        // b: y-component of d
        // cx: x-coord of cell center
        // cy: y-coord of cell center
        // rc: radius of cell center
        // 3 cases: 1) diagonal: a=b 2) horizontal a=d,b=0 3) vertikal: b=d,a=0
        // E.g. for diagonal case
        // d = 1/sqrt(2)(cx+cy)+-sqrt((1/sqrt(2)(cx+cy))^2-rc^2+r^2)
        // with p1 = 1/sqrt(2)(cx+cy) this can be generalized for all 3 cases
        // d = p1+-sqrt(p1^2+r^2-rc^2) = p1+-p2 e.g. p=cx for dir 0.
        // Finally, d is computed normalized as q=d/cellLength_of_direction
        MFloat p1 = 0;
        MFloat p2 = 0;
        for(MInt j = 0; j < nDim; j++) {
          p1 -= m_solver->a_coordinate(currentId, j) * (Ld::idFld(dist, j) - 1.0) * 2.0 * cellHalfLength / dxSq;
        }
        p2 = sqrt(p1 * p1 + (radius * radius - cellRadiusSq) / dxSq);

        // Quadratic equation has two solutions. Either q>1/2 or q<=1/2 is correct.
        if(p1 + p2 >= 0.0 && fabs(p1 + p2) <= 0.5) {
          // In d2q9, the not normalized wall distance is stored
          m_bndCells[i].m_distances[dist] = (p1 + p2) * 2.0 * dxB2;
        }
        if(p1 - p2 >= 0.0 && fabs(p1 - p2) <= 0.5) {
          m_bndCells[i].m_distances[dist] = (p1 - p2) * 2.0 * dxB2;
        }
      }
    }
  }

  else {
    std::stringstream ss;
    ss << "ERROR: Unknown interpolationDistMethod: " << m_interpolationDistMethod << std::endl;
    m_log << ss.str();
    TERMM(1, ss.str());
  }
}

template <>
inline void LbBndCndDxQy<2, 9, maia::lb::LbSysEqnCompressible<2, 9>>::calculateWallDistances() {
  calculateWallDistances2D();
}

template <>
inline void LbBndCndDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::calculateWallDistances() {
  calculateWallDistances2D();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::addWallDistanceFieldToOutputFile(ParallelIo& parallelIo,
                                                                         const MBool writeHeader,
                                                                         const MBool writeData) {
  // total #of cells in the grid == global cell index of the last cell in the last process:
  const MInt totalNoCells = m_solver->domainOffset(m_solver->noDomains());
  const MInt internalCells = m_solver->noInternalCells();
  if(writeHeader) {
    // write header
    {
      const MString varName = "isFluid";
      parallelIo.defineArray(maia::parallel_io::PIO_INT, varName, totalNoCells);
      parallelIo.setAttribute(varName, "name", varName);
    }
    for(MInt i = 0; i < nDist - 1; i++) {
      const MString varName = "distance_" + std::to_string(i);
      parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, varName, totalNoCells);
      parallelIo.setAttribute(varName, "name", varName);
    }
  }
  if(writeData) {
    { // write isFluid data
      MIntScratchSpace tmp(internalCells, AT_, "tmp");
      const MString varName = "isFluid";
      for(MInt i = 0; i < internalCells; i++) {
        tmp[i] = -1;
      }
      for(auto& bndCell : m_bndCells) {
        if(m_solver->a_isHalo(bndCell.m_cellId)) continue;
        tmp[bndCell.m_cellId] = (bndCell.m_isFluid) ? 1 : 0;
      }
      parallelIo.writeArray(tmp.getPointer(), varName);
    }
    MFloatScratchSpace tmp(internalCells, AT_, "tmp");
    for(MInt j = 0; j < nDist - 1; j++) {
      const MString varName = "distance_" + std::to_string(j);
      for(MInt i = 0; i < internalCells; i++) {
        tmp[i] = -1.0;
      }
      for(auto& bndCell : m_bndCells) {
        if(m_solver->a_isHalo(bndCell.m_cellId)) continue;
        tmp[bndCell.m_cellId] = (bndCell.m_distances[j] < 2.0) ? bndCell.m_distances[j] : -1.0; // set -1.0 for no cut
      }
      parallelIo.writeArray(tmp.getPointer(), varName);
    }
  }
}

/** \brief Special initialization for boundary conditions only for 3D!
 *
 *  For each cell in the desired segment(s) the smalles distance
 *  to the rim of the boundary is determined.
 *
 *  \param[in] index Boundary index
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::parabolicInflow(MInt index) {
  TRACE();

  m_log << "  + Generating parabolic inflow..." << std::endl;

  MInt noInflowSegments = 0;
  MInt* inflowSegmentIds = nullptr;

  /*! \page propertyPage1
    \section inflowSegmentIds
    <code>MFloat LbBndCndDxQy::m_inflowSegmentIds</code>\n
    This defines the inflow segments to control the flow for (lbControlInflow)
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  if(Context::propertyExists("inflowSegmentIds", m_solverId)) {
    noInflowSegments = Context::propertyLength("inflowSegmentIds", m_solverId);
    mAlloc(inflowSegmentIds, noInflowSegments, "inflowSegmentIds", AT_);
    for(MInt i = 0; i < noInflowSegments; i++) {
      inflowSegmentIds[i] = Context::getSolverProperty<MInt>("inflowSegmentIds", m_solverId, AT_, i);
    }
  } else {
    m_log << "    - none found" << std::endl << std::endl;
    return;
  }

  MInt num, currentId, nearestNode, neighborNode0, neighborNode1;
  MFloat v1[nDim] = {0.0};
  MFloat v2[nDim] = {0.0};
  MFloat v3[nDim] = {0.0};
  MFloat dmin, dtemp, maxRadius;

  MFloat v1_length, v2_length, v3_length, angle0, angle1, proj0, proj1;


  // loop over all Segments
  for(MInt i = 0; i < noInflowSegments; i++) {
    MInt own = m_solver->m_geometry->m_ownSegmentId[inflowSegmentIds[i]];
    MFloat** a = nullptr;

    if(m_solver->m_geometry->m_parallelGeometry && nDim == 3)
      // in this case all need to participate
      a = allOwnersGetBoundaryVertices(inflowSegmentIds[i], &num);
    else if(own)
      a = m_solver->m_geometry->GetBoundaryVertices(inflowSegmentIds[i], nullptr, nullptr, 0, &num);

    // do the following only if this segment Id exists in this domain!
    if(own) {
      MInt id = -1;

      for(MInt j = 0; j < (MInt)(m_bndCndSegIds.size()); j++)
        if(m_bndCndSegIds[j] == inflowSegmentIds[i]) id = j;

      if(id == -1) continue;
      m_log << "    - segment " << inflowSegmentIds[i] << std::endl;

      // for debugging
      // writeBoundaryVTK(a,num,inflowSegmentIds[i]);

      maxRadius = m_solver->m_geometry->getBndMaxRadius(a, num);
      m_log << "      * max. radius:                 " << maxRadius << std::endl;

      // loop over all m_bndCells in Segment
      for(MInt j = m_bndCndOffsets[id]; j < m_bndCndOffsets[id + 1]; j++) {
        currentId = m_bndCells[j].m_cellId;

        nearestNode = 0;
        dmin = 0;
        for(MInt k = 0; k < nDim; k++) {
          dmin += (m_solver->a_coordinate(currentId, k) - a[0][k]) * (m_solver->a_coordinate(currentId, k) - a[0][k]);
        }

        // run over all rim-nodes and look for nearest node
        for(MInt k = 1; k < num; k++) {
          dtemp = 0;
          for(MInt l = 0; l < nDim; l++) {
            dtemp +=
                (m_solver->a_coordinate(currentId, l) - a[k][l]) * (m_solver->a_coordinate(currentId, l) - a[k][l]);
          }
          if(dtemp < dmin) {
            dmin = dtemp;
            nearestNode = k;
          }
        }

        // vector from nearest node to cell
        v1_length = 0;
        for(MInt k = 0; k < nDim; k++) {
          v1[k] = m_solver->a_coordinate(currentId, k) - a[nearestNode][k];
          v1_length += v1[k] * v1[k];
        }
        v1_length = sqrt(v1_length);

        dmin = v1_length;


        // vector from nearest node to neighbor node
        neighborNode0 = (nearestNode + 1) % num;

        if(nearestNode == 0)
          neighborNode1 = num - 1;
        else
          neighborNode1 = nearestNode - 1;

        v2_length = 0;
        v3_length = 0;
        for(MInt k = 0; k < nDim; k++) {
          v2[k] = a[neighborNode0][k] - a[nearestNode][k];
          v3[k] = a[neighborNode1][k] - a[nearestNode][k];
          v2_length += v2[k] * v2[k];
          v3_length += v3[k] * v3[k];
        }
        v2_length = sqrt(v2_length);
        v3_length = sqrt(v3_length);

        // check if projections are both negative, then we are done
        proj0 = 0.0;
        proj1 = 0.0;
        for(MInt d = 0; d < nDim; d++) {
          proj0 += (v1[d] * v2[d]);
          proj1 += (v1[d] * v3[d]);
        }

        if(proj0 <= 0 && proj1 <= 0)
          dmin = v1_length;
        else {
          // Get angles:
          angle0 = acos(proj0 / (v1_length * v2_length));
          angle1 = acos(proj1 / (v1_length * v3_length));

          if(angle0 <= angle1)
            dmin = sin(angle0) * v1_length;
          else
            dmin = sin(angle1) * v1_length;
        }

        m_bndCells[j].m_multiplier = dmin;
      }

      switch(index) {
        case 1: // poiseuille profile for pipe flow

          m_log << "      * poiseuille inflow:           |v|=Ma*c and v ~ r^2" << std::endl;

          // normalize distance and calculate parabolic profile
          for(MInt j = m_bndCndOffsets[id]; j < m_bndCndOffsets[id + 1]; j++) {
            m_bndCells[j].m_multiplier = m_bndCells[j].m_multiplier / maxRadius;
            m_bndCells[j].m_multiplier = 1.0 - m_bndCells[j].m_multiplier;
            m_bndCells[j].m_multiplier = 2.0 * (1.0 - (m_bndCells[j].m_multiplier * m_bndCells[j].m_multiplier));
          }
          break;

        case 2: // r^4 profile (broader than poiseuille)

          m_log << "      * poiseille inflow (broader): |v|=Ma*c and v ~ r^4" << std::endl;

          // normalize distance and calculate parabolic profile
          for(MInt j = m_bndCndOffsets[id]; j < m_bndCndOffsets[id + 1]; j++) {
            m_bndCells[j].m_multiplier = m_bndCells[j].m_multiplier / maxRadius;
            m_bndCells[j].m_multiplier = 1.0 - m_bndCells[j].m_multiplier;
            m_bndCells[j].m_multiplier = 1.5
                                         * (1.0
                                            - (m_bndCells[j].m_multiplier * m_bndCells[j].m_multiplier
                                               * m_bndCells[j].m_multiplier * m_bndCells[j].m_multiplier));
          }
          break;

        default:
          m_log << "lbbndcnddxqy::parabolicInflow() unknown index for parabolic inflow " << std::endl;
          std::stringstream errorMessage;
          errorMessage << "lbbndcnddxqy::parabolicInflow() unknown index for parabolic inflow";
          DEBUG("lbbndcnddxqy::parabolicInflow() unknown index for parabolic inflow" << std::endl,
                MAIA_DEBUG_ASSERTION);
          TERMM(1, errorMessage.str());
      }
    }
  }
  mDeallocate(inflowSegmentIds);
}

/**
 * \brief Obtain all boundary vertices for a given segmentId
 *
 * \param[in] segmentId Segment id for which the vertices are returned
 * \param[out] num Number of vertices found for the given segment id
 *
 * /returns Boundary vertices
 */
template <MInt nDim, MInt nDist, class SysEqn>
MFloat** LbBndCndDxQy<nDim, nDist, SysEqn>::allOwnersGetBoundaryVertices(MInt segmentId, MInt* num) {
  TRACE();

  MInt own;
  MInt sumowners;
  MInt firstOwner;
  MIntScratchSpace owners(m_solver->noDomains(), AT_, "owners");

  m_log << "      * segment owned by:          ";
  m_solver->m_geometry->determineSegmentOwnership(segmentId, &own, &sumowners, &firstOwner, owners.getPointer());
  for(MInt d = 0; d < m_solver->noDomains(); d++)
    if(owners[d] > 0) m_log << d << " ";

  m_log << std::endl;
  m_log << "      * sum of owners:             " << sumowners << std::endl;
  m_log << "      * root of communication:     " << firstOwner << std::endl;


  // my domain owns the whole segment
  if(sumowners == 1)
    return m_solver->m_geometry->GetBoundaryVertices(segmentId, nullptr, nullptr, 0, num);
  else {
    // build communicator
    MPI_Comm charComm;
    m_solver->createMPIComm(owners.getPointer(), sumowners, &charComm);

    if(own) {
      // collect the triangles for testing first

      MInt offStart = 0;
      MInt offEnd = 0;

      if(m_solver->m_geometry->m_parallelGeometry) {
        offStart = m_solver->m_geometry->m_segmentOffsets[segmentId];
        offEnd = m_solver->m_geometry->m_segmentOffsets[segmentId + 1];
      } else {
        offStart = m_solver->m_geometry->m_segmentOffsetsWithoutMB[segmentId];
        offEnd = m_solver->m_geometry->m_segmentOffsetsWithoutMB[segmentId + 1];
      }
      MInt numElements = offEnd - offStart;

      m_log << "      * number of local triangles: " << numElements << std::endl;
      m_log << "      * segment offsets:           " << offStart << " - " << offEnd << std::endl;

      // (normals + vertices)
      MInt noTriInfo = nDim * nDim;
      MIntScratchSpace myOriginalIds(numElements, AT_, "myOriginalIds");
      MFloatScratchSpace segTriangles(noTriInfo * numElements, AT_, "segTriangles");

      for(MInt t = offStart, i = 0, j = 0; t < offEnd; t++, i++) {
        myOriginalIds[i] = m_solver->m_geometry->elements[t].m_originalId;
        for(MInt v = 0; v < nDim; v++)
          for(MInt d = 0; d < nDim; d++, j++)
            segTriangles[j] = m_solver->m_geometry->elements[t].m_vertices[v][d];
      }

      MIntScratchSpace numElemPerCPU(sumowners, AT_, "numElemPerCPU");
      MPI_Allgather(&numElements, 1, MPI_INT, numElemPerCPU.getPointer(), 1, MPI_INT, charComm, AT_, "numElements",
                    "numElemPerCPU.getPointer()");
      MIntScratchSpace numTriInfoPerCPU(sumowners, AT_, "numTriInfoPerCPU");

      m_log << "      * triangles per domain:      ";
      MInt sumallelem = 0;
      for(MInt i = 0; i < sumowners; i++) {
        m_log << numElemPerCPU[i] << " ";
        sumallelem += numElemPerCPU[i];
        numTriInfoPerCPU[i] = numElemPerCPU[i] * noTriInfo;
      }
      m_log << std::endl;
      m_log << "      * sum of global triangles:   " << sumallelem << std::endl;

      MIntScratchSpace displOrig(sumowners, AT_, "displOrig");
      MIntScratchSpace displTris(sumowners, AT_, "displTris");
      displOrig[0] = 0;
      displTris[0] = 0;
      for(MInt d = 1; d < sumowners; d++) {
        displOrig[d] = displOrig[d - 1] + numElemPerCPU[d - 1];
        displTris[d] = displTris[d - 1] + numTriInfoPerCPU[d - 1];
      }

      MIntScratchSpace allOriginalIds(sumallelem, AT_, "allOriginalIds");
      MPI_Allgatherv(myOriginalIds.getPointer(), numElements, MPI_INT, allOriginalIds.getPointer(),
                     numElemPerCPU.getPointer(), displOrig.getPointer(), MPI_INT, charComm, AT_,
                     "myOriginalIds.getPointer()", "allOriginalIds.getPointer()");

      MFloatScratchSpace allSegTriangles(noTriInfo * sumallelem, AT_, "allSegTriangles");
      MPI_Allgatherv(segTriangles.getPointer(), noTriInfo * numElements, MPI_DOUBLE, allSegTriangles.getPointer(),
                     numTriInfoPerCPU.getPointer(), displTris.getPointer(), MPI_DOUBLE, charComm, AT_,
                     "segTriangles.getPointer()", "allSegTriangles.getPointer()");

      std::set<MInt> uniqueTriangles;
      for(MInt i = 0; i < sumallelem; i++)
        uniqueTriangles.insert(allOriginalIds[i]);

      MInt noUniqueTris = uniqueTriangles.size();
      m_log << "      * sum of unique triangles:   " << noUniqueTris << std::endl;

      MIntScratchSpace dbl(noUniqueTris, AT_, "dbl");
      for(MInt i = 0; i < noUniqueTris; i++)
        dbl[i] = 0;

      // this contains the offsets in the list of triangles which we want to use for the calculation
      MIntScratchSpace keepOffsets(noUniqueTris, AT_, "keepOffsets");
      for(MInt i = 0, j = 0; i < sumallelem; i++) {
        MInt dist = distance(uniqueTriangles.begin(), uniqueTriangles.find(allOriginalIds[i]));
        if(dist != noUniqueTris && dbl[dist] == 0) {
          keepOffsets[j] = i * noTriInfo;
          dbl[dist] = 1;
          j++;
        }
      }

      return m_solver->m_geometry->GetBoundaryVertices(segmentId, allSegTriangles.getPointer(),
                                                       keepOffsets.getPointer(), noUniqueTris, num);
    } else
      return nullptr;
  }
}

/**
 * \brief Interpolated bounce back for moving walls using the linear BFL scheme
 *
 * Performs a bounce back for a single boundary cell using the linear Bouzidi-Firdaouss-Lallemand (BFL) scheme.
 * If the boundary cell is inside, it is updated. If the boundary cell is outside, its inside neighbours are updated.
 * For reference see https://doi.org/10.1063/1.1399290
 *
 * \author Moritz Waldmann, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] cellIndex Boundary cell id from where the bounce back is performed.
 * \param[in] set Level set for which the bounce back is performed
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_lin(const MInt cellIndex, const MInt set) {
  TRACE();

  const MInt pCellId = m_solver->m_G0CellList[cellIndex];
  MFloat rho = m_solver->a_variable(pCellId, PV->RHO);
  std::array<MFloat, nDim> uW{};
  getBoundaryVelocityMb(cellIndex, uW.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          const MFloat F2q = F2 * q;
          m_solver->a_oldDistribution(pCellId, opposite) = (F2q * m_solver->a_distribution(pCellId, j))
                                                           + ((F1 - F2q) * m_solver->a_oldDistribution(pCellId, j))
                                                           + firstMomentSourceTerm(uW.data(), rho, opposite);
#ifndef WAR_NVHPC_PSTL
          incrementForces(pCellId, cellIndex, j, uW.data());
#endif
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else {
          m_solver->a_oldDistribution(pCellId, opposite) =
              m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
        }
      }

      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      // else if(m_bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      // else if( m_wallBoundaryCellList[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistribution(pCellId, opposite) =
            m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
      }
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif
      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);
        rho = m_solver->a_variable(neighborId, PV->RHO);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(q < 0.5) {
            const MFloat F2q = F2 * (F1 - q);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
              m_solver->a_oldDistribution(neighborId, j) =
                  (m_solver->a_distribution(neighborId, opposite) / F2q)
                  + (m_solver->a_distribution(neighborId, j) * (F2q - F1) / F2q)
                  + (F1 / F2q) * firstMomentSourceTerm(uW.data(), rho, j);

#ifndef WAR_NVHPC_PSTL
              incrementForces(neighborId, cellIndex, opposite, uW.data());
#endif
            }
          }
        }
      }
    } // end of the loop over all PPDF directions in case 2
    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDists(pCellId, F1, uW.data());
  }
}


/**
 * \brief Interpolated bounce back for moving walls using the quadratic BFL scheme
 *
 * Performs a bounce back for a single boundary cell using the quadratic Bouzidi-Firdaouss-Lallemand (BFL) scheme.
 * If the boundary cell is inside, it is updated. If the boundary cell is outside, its inside neighbours are updated.
 * For reference see https://doi.org/10.1063/1.1399290
 *
 * \author Moritz Waldmann, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] cellIndex Boundary cell id from where the bounce back is performed.
 * \param[in] set Level set for which the bounce back is performed
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_qua(const MInt cellIndex, const MInt set) {
  TRACE();

  const MInt pCellId = m_solver->m_G0CellList[cellIndex];
  MFloat rho = m_solver->a_variable(pCellId, PV->RHO);
  std::array<MFloat, nDim> uW{};
  getBoundaryVelocityMb(cellIndex, uW.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite)) {
          const MInt neighborId = m_solver->c_neighborId(pCellId, opposite);
          // const MInt nneighborId = m_solver->c_neighborId(neighborId, opposite);
          // std::cout << "nn" << nneighborId << std::endl;
          const MFloat F2q = F2 * q;
          const MFloat Fq = q;

          MFloat quadratic = Fq * (F2q + 1) * m_solver->a_distribution(pCellId, j)
                             + (1 + F2q) * (1 - F2q) * m_solver->a_distribution(neighborId, j)
                             - Fq * (1 - F2q) * m_solver->a_oldDistribution(neighborId, j)
                             + firstMomentSourceTerm(uW.data(), rho, opposite);

          m_solver->a_oldDistribution(pCellId, opposite) = quadratic;
#ifndef WAR_NVHPC_PSTL
          incrementForces(pCellId, cellIndex, j, uW.data());
#endif
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else {
          m_solver->a_oldDistribution(pCellId, opposite) =
              m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
        }
      }

      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistribution(pCellId, opposite) =
            m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
      }
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);
        rho = m_solver->a_variable(neighborId, PV->RHO);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(q < 0.5) {
            const MFloat F2q = F2 * (F1 - q);
            const MFloat Fq = (1 - q);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
              if(m_solver->a_hasNeighbor(neighborId, j)) {
                const MInt nneighborId = m_solver->c_neighborId(neighborId, j);

                MFloat quadratic =
                    1 / Fq / (F2q + 1)
                        * (m_solver->a_distribution(neighborId, opposite) + firstMomentSourceTerm(uW.data(), rho, j))
                    + (F2q - 1) / Fq * m_solver->a_distribution(neighborId, j)
                    - (F2q - 1) / (F2q + 1) * m_solver->a_distribution(nneighborId, j);

                m_solver->a_oldDistribution(neighborId, j) = quadratic;
              } else {
                // Linear fallback for confined direction
                MFloat linear = (m_solver->a_distribution(neighborId, opposite) / F2q)
                                + (m_solver->a_distribution(neighborId, j) * (F2q - F1) / F2q)
                                + (F1 / F2q) * firstMomentSourceTerm(uW.data(), rho, j);

                m_solver->a_oldDistribution(neighborId, j) = linear;
              }
#ifndef WAR_NVHPC_PSTL
              incrementForces(neighborId, cellIndex, opposite, uW.data());
#endif
            }
          }
        }
      }
    } // end of the loop over all PPDF directions in case 2
    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDists(pCellId, F1, uW.data());
  }
}


/**
 * \brief Interpolated bounce back for moving walls using the quadratic scheme by Yu et al.
 *
 * Performs a bounce back for a single boundary cell using the quadratic scheme by Yu et al.
 * If the boundary cell is inside, it is updated. If the boundary cell is outside, its inside neighbours are updated.
 * For reference see https://doi.org/10.1016/S0376-0421(03)00003-4
 *
 * \author Moritz Waldmann, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] cellIndex Boundary cell id from where the bounce back is performed.
 * \param[in] set Level set for which the bounce back is performed
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Yu_qua(const MInt cellIndex, const MInt set) {
  TRACE();

  const MInt pCellId = m_solver->m_G0CellList[cellIndex];
  MFloat rho = m_solver->a_variable(pCellId, PV->RHO);
  std::array<MFloat, nDim> uW{};
  getBoundaryVelocityMb(cellIndex, uW.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          const MInt neighborId = m_solver->c_neighborId(pCellId, opposite);
          const MInt nneighborId = m_solver->c_neighborId(neighborId, opposite);

          const MFloat F2q = F2 * q;
          const MFloat Fq = q;

          // Calculate temporary distribution at position which will propagrate exactly to the wall
          const MFloat wall = Fq * (Fq + 1) / 2 * m_solver->a_distribution(pCellId, j)
                              + (1 - Fq) * (1 + Fq) * m_solver->a_distribution(neighborId, j)
                              - Fq * (1 - Fq) / 2 * m_solver->a_distribution(nneighborId, j);

          // Instantaneous bounce-back at the wall
          const MFloat wall_bb = wall + firstMomentSourceTerm(uW.data(), rho, opposite);

          // Interpolate bounce-back value
          const MFloat interp = 2 / (1 + Fq) / (2 + Fq) * wall_bb
                                + F2q / (1 + Fq) * m_solver->a_distribution(pCellId, opposite)
                                - Fq / (2 + Fq) * m_solver->a_distribution(neighborId, opposite);

          m_solver->a_oldDistribution(pCellId, opposite) = interp;
#ifndef WAR_NVHPC_PSTL
          incrementForces(pCellId, cellIndex, j, uW.data());
#endif
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else {
          m_solver->a_oldDistribution(pCellId, opposite) =
              m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
        }
      }

      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistribution(pCellId, opposite) =
            m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW.data(), rho, opposite);
      }
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);
        rho = m_solver->a_variable(neighborId, PV->RHO);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(q < 0.5) {
            const MFloat F2q = F2 * (F1 - q);
            const MFloat Fq = (1 - q);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
              const MInt nneighborId = m_solver->c_neighborId(neighborId, j);
              const MInt nnneighborId = m_solver->c_neighborId(nneighborId, j);

              // Calculate temporary distribution at position which will propagrate exactly to the wall
              const MFloat wall = Fq * (Fq + 1) / 2 * m_solver->a_distribution(neighborId, opposite)
                                  + (1 - Fq) * (1 + Fq) * m_solver->a_distribution(nneighborId, opposite)
                                  - Fq * (1 - Fq) / 2 * m_solver->a_distribution(nnneighborId, opposite);

              // Instantaneous bounce-back at the wall
              const MFloat wall_bb = wall + firstMomentSourceTerm(uW.data(), rho, j);


              // Interpolate bounce-back value
              const MFloat interp = 2 / (1 + Fq) / (2 + Fq) * wall_bb
                                    + F2q / (1 + Fq) * m_solver->a_distribution(neighborId, j)
                                    - Fq / (2 + Fq) * m_solver->a_distribution(nneighborId, j);


              m_solver->a_oldDistribution(neighborId, j) = interp;
#ifndef WAR_NVHPC_PSTL
              incrementForces(neighborId, cellIndex, opposite, uW.data());
#endif
            }
          }
        }
      }
    } // end of the loop over all PPDF directions in case 2
    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDists(pCellId, F1, uW.data());
  }
}


/**
 * \brief Calculates the forces exerted by a single distribution during bounce back
 *
 * The force is calculated by the momentum exchange method (MEA) in its Galilei invariant formulation by Caiazzo et al.
 * The result is stored per distribution for each boundary cell.
 * This allows a more accurate torque calculation later on.
 * For reference see https://doi.org/10.1016/j.camwa.2007.08.004
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] cellId Cell id of the boundary cell
 * \param[in] mbCellId Moving boundary cell id
 * \param[in] j Distribution which is intersecting the boundary
 * \param[in] uW Boundary velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::incrementForces(const MInt cellId, const MInt mbCellId, const MInt j,
                                                               const MFloat* uW) {
  const MInt opposite = Ld::oppositeDist(j);

  const MFloat& Fin = m_solver->a_distribution(cellId, j);
  const MFloat& Fout = m_solver->a_oldDistribution(cellId, opposite);

  // Galilei invariant momentum exchange method
  for(MInt d = 0; d < nDim; d++) {
    const MFloat Vin = Ld::ppdfDir(j, d) - uW[d];
    const MFloat Vout = Ld::ppdfDir(opposite, d) - uW[d];

    m_boundaryCellsMb.force(mbCellId, j, d) = (Fin * Vin - Fout * Vout);
  }
}


/**
 * \brief Reads boundary velocity from the moving boundary cell collector
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  mbCellId Moving boundary cell id
 * \param[out] uW Boundary velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::getBoundaryVelocityMb(const MInt mbCellId, MFloat* uW) {
  for(MInt n = 0; n < nDim; n++) {
    uW[n] = m_boundaryCellsMb.velocity(mbCellId, n);
  }
}

/**
 * \brief Reads boundary velocity from properties
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  mbCellId Moving boundary cell id
 * \param[out] uW Boundary velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::getBoundaryVelocity(const MInt index, MFloat* uW) {
  for(MInt d = 0; d < nDim; d++)
    uW[d] = F0;
  for(MInt i = 0; i < m_lbNoMovingWalls; i++) {
    if(m_segIdMovingWalls[i] == m_bndCndSegIds[index]) {
      for(MInt d = 0; d < nDim; d++) {
        uW[d] = m_lbWallVelocity[i * nDim + d];
      }
      break;
    }
  }
}

/**
 * \brief Reads boundary temperature from properties
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  mbCellId Moving boundary cell id
 * \param[out] uW Boundary velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline MFloat LbBndCndDxQy<nDim, nDist, SysEqn>::getBoundaryTemperature(const MInt index) {
  MFloat wT = m_solver->m_initTemperatureKelvin;
  for(MInt i = 0; i < m_lbNoHeatedWalls; i++) {
    if(m_segIdHeatedWalls[i] == m_bndCndSegIds[index]) {
      wT = m_lbWallTemperature[i];
      break;
    }
  }
  return wT;
}

/**
 * \brief Momentum source term used in bounce back schemes for moving boundaries
 *
 * Attention: This term corresponds to the 'compressible' formulation of the equilibirum distribution
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] uW Boundary velocity
 * \param[in] rho Fluid density at boundary
 * \param[in] dist Distribution to be set by the bounce back scheme
 *
 * /returns Momentum source term
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline MFloat LbBndCndDxQy<nDim, nDist, SysEqn>::firstMomentSourceTerm(const MFloat* const uW, const MFloat rho,
                                                                       const MInt dist) {
#ifdef WAR_NVHPC_PSTL
  MFloat ppdfDir[nDim] = {F0};
  for(MInt d = 0; d < nDim; d++) {
    ppdfDir[d] = m_solver->m_ppdfDir[dist * nDim + d];
  }
  const MFloat scalarProduct = std::inner_product(uW, uW + nDim, &ppdfDir[0], .0);
  const MFloat weight = m_solver->m_tp[m_solver->m_distType[dist]];
#else
  const MFloat scalarProduct = std::inner_product(uW, uW + nDim, lbDescriptor::ppdfDir<nDim>[dist], .0);
  const MFloat weight = Ld::tp(Ld::distType(dist));
#endif

  return 2.0 * F1BCSsq * weight * rho * scalarProduct;
}


/** \LBBC{secLBBC_bc66666, bc66666, MB}
 *
 *  \author Moritz Waldmann, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date 01.2020
 *
 * Applies no-slip boundary condition to moving boundaries
 *  The bounce back scheme can be chosen via property.
 *
 *  \param[in] set The boundary condition is applied to all cells which belong to this level set.<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc66666(MInt set)
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc66666(MInt set) {
  TRACE();

  const MInt noMbCells = m_boundaryCellsMb.size();

#ifdef WAR_NVHPC_PSTL
  mAlloc(m_distances, noMbCells, nDist - 1, "distances", F1, AT_);
  for(MInt i = 0; i < noMbCells; i++) {
    const MInt pCellId = m_solver->m_G0CellList[i];
    for(MInt d = 0; d < nDist - 1; d++) {
      m_distances[i][d] = getDistanceMb(pCellId, i, d);
    }
  }
#endif

  maia::parallelFor<true>(0, noMbCells, [=](MInt mbId) {

#ifdef WAR_NVHPC_PSTL
    // Perform bounce back set by properties
    if(m_bounceBackFunctionMbCase == 0) {
      interpolatedBounceBackMb_Bouzidi_lin(mbId, set);
    } else if(m_bounceBackFunctionMbCase == 1) {
      interpolatedBounceBackMb_Bouzidi_qua(mbId, set);
    } else {
      interpolatedBounceBackMb_Yu_qua(mbId, set);
    }
#else
    // Perform bounce back set by properties
    (this->*bounceBackFunctionMb)(mbId, set);
#endif

    if(m_solver->m_isThermal) {
      interpolatedBounceBackMb_Bouzidi_lin_thermal(mbId, set);
    }
    if(m_solver->m_isTransport) {
      interpolatedBounceBackMb_Bouzidi_lin_transport(mbId, set);
    }
  });

#ifdef WAR_NVHPC_PSTL
  mDeallocate(m_distances);
#endif
  if(m_calcWallForces && m_solver->m_currentNoG0Cells > 0) {
    calculateWallForcesMb(set);
  }
}

/** \LBBC{secLBBC_bc66668, bc66668, MB}
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date 01.2020
 *
 * Applies pressure boundary condition to moving boundaries
 *
 *  The quadratic Bouzidi anti bounce back is used.
 *
 *  \param[in] set The boundary condition is applied to all cells which belong to this level set.<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc66668(MInt set)
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc66668(MInt set) {
  TRACE();

  MInt noMbCells = m_boundaryCellsMb.size();
  if(m_solver->noDomains() > 1) {
    MPI_Allreduce(&noMbCells, &noMbCells, 1, MPI_INT, MPI_SUM, m_solver->mpiComm(), AT_, "noMbCells", "noMbCells");
  }
  if(!m_solver->domainId()) {
    std::cout << "Apply free surface bc to " << noMbCells << " no mbCells" << std::endl;
  }

  for(MInt mbCellId = 0; mbCellId < m_boundaryCellsMb.size(); mbCellId++) {
    const MInt pCellId = m_boundaryCellsMb.cellId(mbCellId);

    ASSERT(m_solver->a_isG0CandidateOfSet(pCellId, (set - m_solver->m_levelSetId)), "Inconsistent collector!");

    interpolatedAntiBounceBackMb_Bouzidi_qua(mbCellId, set);
  }
}


/** \LBBC{secLBBC_bc20000, bc20000, 2000}
 *
 * Lattice Boltzmann enhanced no slip (bounce back) condition for
 * inclined walls. Bouzidi 2001 (aka "BFL rule")
 *
 * There are "dry" nodes and "wet" nodes (cell is center inside/outside).
 * Dry nodes do not take part in the collision/streaming process.
 * An interpolated distribution is reflected at a certain point
 * between the cell centers. If the wall is located exactly halfway between
 * the nodes this condition reduces to simple halfway bounce-back.
 * Must be performed after propagation.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20000(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20000(MInt index) {
  TRACE();

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    interpolatedBounceBackSingleSpecies(i, uW.data());
  });

  if(m_calcWallForces) {
    calculateWallForces(index);
  }
}


/** \LBBC{secLBBC_bc20000, bc20001, 2001}
 *
 * Lattice Boltzmann no slip (halfway bounce-back) condition
 *
 * If no neihbor exists in outgoing direction, the outgoing
 * distribution is reflected.
 * -> The wall is located at the outer cell edge.
 * Must be performed after propagation.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20001(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20001(MInt index) {
  TRACE();

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep;
  const MInt begin = m_bndCndOffsets[index];
  const MInt end = m_bndCndOffsets[index + 1];
  const MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    const MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    if(m_solver->a_isHalo(m_bndCells[i].m_cellId)) return;

    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MInt opposite = Ld::oppositeDist(j);
#endif
      if(m_solver->a_hasNeighbor(m_bndCells[i].m_cellId, opposite) == 0) {
        // leave out interface cells
        if(m_solver->c_parentId(m_bndCells[i].m_cellId) > -1) {
          if(m_solver->a_hasNeighbor(m_solver->c_parentId(m_bndCells[i].m_cellId), opposite)) continue;
        }
        m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) =
            m_solver->a_distribution(m_bndCells[i].m_cellId, opposite);
      }
    }
  });
}


/** \LBBC{secLBBC_bc20002, bc20002, 2002}
 *
 * Lattice Boltzmann no slip (bounce back) condition
 *
 * If no neihbor exists in incoming direction, the incoming
 * distribution is reflected.
 * -> The wall is located at the cell center!
 * Must be performed after collision and BEFORE propagation.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20002(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20002(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndCells[i].m_cellId;
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // leave out halo cells
    if(m_solver->a_isHalo(cellId)) continue;

    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(cellId, Ld::oppositeDist(j)) == 0) {
        m_solver->a_oldDistribution(cellId, j) = m_solver->a_oldDistribution(cellId, Ld::oppositeDist(j));
      }
    }
  }
}


/** \LBBC{secLBBC_bc20003, bc20003, 2003}
 *
 * Lattice Boltzmann no slip (bounce back) condition
 *
 * no slip condition for periodic channel
 * pressure is extrapolated
 * velocity is set to zero
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20003(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20003(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if( (globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    const MInt cellId = m_bndCells[i].m_cellId;

    // extrapolate inner density values
    MInt neighborId{};
    if(m_solver->a_hasNeighbor(cellId, 2) == 0) {
      neighborId = m_solver->c_neighborId(cellId, 3);
    } else {
      neighborId = m_solver->c_neighborId(cellId, 2);
    }

    MFloat rho = m_solver->a_variable(neighborId, PV->RHO);

    const MFloat old_rho = m_solver->a_oldVariable(cellId, PV->RHO);

    // non-reflecting bc (Finck, Haenel 2008)
    // rho = ( old_rho + F1BCS * ( sqrt(tmp2) - sqrt(old_tmp2) ) + m_rho1 ) / 2.0 ;
    rho = (old_rho + rho) / 2.0;

    // set velocity to zero
    const std::array<MFloat, nDim> u{};
    const MFloat squaredVelocity{};

    m_solver->setEqDists(cellId, rho, squaredVelocity, u.data());

    m_solver->a_variable(cellId, PV->RHO) = rho;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(cellId, PV->VV[n]) = u[n];
    }
  }
}

/** \LBBC{secLBBC_bc20004, bc20004, 2004}
 *
 * Lattice Boltzmann enhanced no slip (bounce back) condition for
 * inclined walls. Haenel 1998, improved by Mei 1999/2000
 *
 * There are "dry" nodes and "wet" nodes (cell is center inside/outside).
 * Dry nodes do not take part in the collision/streaming process.
 * An interpolated distribution is reflected at a certain point
 * between the cell centers. If the wall is located exactly halfway between
 * the nodes this condition reduces to simple halfway bounce-back.
 * Must be performed after propagation.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20004(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20004(MInt index) {
  TRACE();

  MFloat F2q, rhoF, uF[nDim], uBF[nDim], tmpUF, tmpUBF, tmp2UF, b[2 * nDim];

  const MInt dist1 = Ld::distFld(0);
  const MInt dist2 = Ld::distFld(1) + Ld::distFld(0);

  MInt id, k, tpIndex;
  MFloat fEq, xi;
  MInt nghbrId, tmpDistId;

  // field for force evaluation
  static std::ofstream ofl;

  ScratchSpace<MFloat> distributions(nDist - 1, AT_, "distributions");
  for(MInt i = 0; i < nDist - 1; i++) {
    distributions[i] = 0.0;
  }


  //---------------------------------------------------------------------------------
  // go through m_bndcells with wallbndcnd;
  // the halo cells have to be included, since they may be of relevance in case 2
  //---------------------------------------------------------------------------------
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    id = m_bndCells[i].m_cellId;

    // If there are several levels of boundary cells use only the highest
    // nearest neighbors, must be of the same level
    if(m_solver->a_level(id) < m_solver->maxLevel()) continue;

    m_omega = 2.0 / (1.0 + 6.0 * m_solver->a_nu(id) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(id)));

    rhoF = 1.0; // standard value which is usually overwritten in the following. Only in exceptional case it remains 1

    //----------------------------------------------------------------------------------------------------
    // perform bounce back for distributions piercing the wall - distinquish between inner and outer cells
    //----------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // 1. boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      rhoF = m_solver->a_variable(id, PV->RHO);
      for(MInt n = 0; n < nDim; n++) {
        uF[n] = m_solver->a_variable(id, PV->VV[n]);
        b[2 * n] = -uF[n];
        b[2 * n + 1] = uF[n];
      }

      tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

      for(MInt j = 0; j < nDist - 1; j++) {
        tmpUF = F0;
        tpIndex = 0;
        if(j < dist1) {
          // dxB2 = cellHalfLength;
          tmpUF = b[j];
          tpIndex = 1;
        } else {
          if(j < dist2) {
            // dxB2 = SQRT2 * cellHalfLength;
            k = j - Ld::distFld(0);
            tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          } else {
            // dxB2 = SQRT3 * cellHalfLength;
            if(nDim == 3) {
              k = j - (Ld::distFld(0) + Ld::distFld(1));
              tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
              tpIndex = 3;
            }
          }
        }

        //---------------------------------------------------------------------------------
        // if q < 0.5, perform bounce back
        //---------------------------------------------------------------------------------

        if(m_bndCells[i].m_distances[j] < 0.5) {
          F2q = F2 * m_bndCells[i].m_distances[j]; // q = distance / dx; distance: distance between cell center and wall

          if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(j)) > 0) {
            nghbrId = m_solver->c_neighborId(id, Ld::oppositeDist(j));

            // extrapolate velocity from next inner cell
            tmpUBF = F0;
            for(MInt n = 0; n < nDim; n++) {
              uBF[n] = m_solver->a_variable(nghbrId, n);
              tmpUBF += (Ld::idFld(j, n) - 1) * uBF[n];
            }
          } else {
            tmpUBF = tmpUF;
          }

          fEq = Ld::tp(tpIndex) * rhoF
                * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
          xi = m_omega * (F2q - 1.0) / (1.0 - 2.0 * m_omega);

          // perform interpolated bounce back
          m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) =
              (1.0 - xi) * m_solver->a_distribution(id, j) + xi * fEq;

          // add to momentum sum
          distributions[Ld::oppositeDist(j)] +=
              m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) + m_solver->a_distribution(id, j);

        }

        // if the cell is located at a corner, there might be distributions missing
        else {
          // In this case simple bounce back is used. This could decrease the overall accuracy.
          if(m_solver->a_hasNeighbor(id, j) == 0) {
            // do simple bounce back
            m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) = m_solver->a_distribution(id, j);

            // add to momentum sum
            distributions[Ld::oppositeDist(j)] +=
                m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) + m_solver->a_distribution(id, j);
          }
        }

      } // end of loop over all directions

    }

    //---------------------------------------------------------------------------------
    // 2. boundary cell is outside fluid
    //---------------------------------------------------------------------------------

    // The outer cell does not belong to the flow field, thus the macroscopic values are held constant.
    // It is possible that there are cells without any cutting velocity! Those have to be controlled, too.
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        // if ( j < 6 )
        // 	dxB2 = cellHalfLength;
        // else {
        // 	if ( j < 18)
        // 	    dxB2 = SQRT2 * cellHalfLength;
        // 	else
        // 	    dxB2 = SQRT3 * cellHalfLength;
        // }

        //----------------------------------------------------------------------------------------------------------------
        // if q_innerNeighbor = 1 - q_bndCell >= 0.5, perform bounceback on inner neighbor
        // if the inner neighbor is a halo cell, it can be skipped
        //----------------------------------------------------------------------------------------------------------------
        //   	    if ( m_solver->a_hasNeighbor(id, j) > 0 && m_solver->c_neighborId(id, j) < m_solver->noInternalCells()
        //   ) {
        if(m_solver->a_hasNeighbor(id, j) > 0 && !m_solver->a_isHalo(m_solver->c_neighborId(id, j))) {
          if(m_bndCells[i].m_distances[j] <= 0.5) {
            F2q = F2 - F2 * m_bndCells[i].m_distances[j]; // q = 1 - 0.5*(distance/cellHalfLength);

            nghbrId = m_solver->c_neighborId(id, j);

            rhoF = m_solver->a_variable(nghbrId, PV->RHO);
            for(MInt n = 0; n < nDim; n++) {
              uF[n] = m_solver->a_variable(nghbrId, PV->VV[n]);
              b[2 * n] = -uF[n];
              b[2 * n + 1] = uF[n];
            }

            tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

            tmpUF = F0;
            tpIndex = 0;
            if(j < dist1) {
              tmpUF = b[Ld::oppositeDist(j)];
              tpIndex = 1;
            } else {
              if(j < dist2) {
                k = Ld::oppositeDist(j) - Ld::distFld(0);
                tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
                tpIndex = 2;
              } else {
                if(nDim == 3) {
                  k = Ld::oppositeDist(j) - (Ld::distFld(0) + Ld::distFld(1));
                  tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
                  tpIndex = 3;
                }
              }
            }

            tmpUBF = 0.0;
            for(MInt dim = 0; dim < nDim; dim++) {
              uBF[dim] = (1.0 - 3.0 / F2q) * uF[dim];
              tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
            }

            fEq = Ld::tp(tpIndex) * rhoF
                  * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
            xi = m_omega * (F2q - 1.0) / (1 + 0.5 * m_omega);
            // xi = m_omega * (F2q - 1.0);

            // perform interpolated bounce back
            m_solver->a_oldDistribution(nghbrId, j) =
                (1.0 - xi) * m_solver->a_distribution(nghbrId, Ld::oppositeDist(j)) + xi * fEq;

            // add to momentum sum
            distributions[j] +=
                m_solver->a_oldDistribution(nghbrId, j) + m_solver->a_distribution(nghbrId, Ld::oppositeDist(j));
            //		    distributions[j] += 2.0 * m_solver->a_distribution(neighborId, Ld::oppositeDist(j));
          }
        }

      } // end of loop over all directions

      // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
      // It is possible that there are cells without any cutting velocity! These are considered too.
      m_solver->a_variable(id, PV->RHO) = 1.0;
      m_solver->a_oldVariable(id, PV->RHO) = 1.0;

      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_variable(id, n) = 0.0;
        m_solver->a_oldVariable(id, n) = 0.0;
      }

      // Calculation of eq-distributions in bnd cell
      //--------------------------------------------
      // Calculation of distributions for directions with only one component
      for(MInt j = 0; j < Ld::distFld(0); j++) {
        m_solver->a_oldDistribution(id, j) = Ld::tp(1);
      }

      // Calculation of distributions for directions with two components
      tmpDistId = Ld::distFld(0);
      for(MInt j = 0; j < Ld::distFld(1); j++) {
        m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(2);
      }

      // Calculation of distributions for directions with three components
      tmpDistId = Ld::distFld(0) + Ld::distFld(1);
      for(MInt j = 0; j < Ld::distFld(2); j++) {
        m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(3);
      }

      // Calculation of distribution for rest particle distribution (center)
      m_solver->a_oldDistribution(id, Ld::lastId()) = Ld::tp(0);
    }

  } // end of loop over all m_bndcells
}

/** LBBC{secLBBC_bc20005, bc20005, 2005}
 *
 * Wall bc with extrapolation of non-eq (Guo)
 *
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20005(MInt index) {
  TRACE();

  MFloat q, rhoF, uF[nDim], tmpUF, tmpUBF, tmp2UF, tmp2UBF = 0.0, b[2 * nDim];
  MFloat uBF[nDim] = {F0};
  for(MInt d = 0; d < nDim; d++) {
    uBF[d] = std::numeric_limits<MFloat>::max();
  }

  MInt id = 0, k, tpIndex;
  MFloat fEqW, fEqF;
  MInt nghbrId, tmpDistId;

  // field for force evaluation
  static std::ofstream ofl;

  ScratchSpace<MFloat> distributions(nDist - 1, AT_, "distributions");
  for(MInt i = 0; i < nDist - 1; i++) {
    distributions[i] = 0.0;
  }


  // go through m_bndcells with wallbndcnd
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // If there are several levels of boundary cells use only the highest
    // nearest neighbors, must be of the same level
    if(m_solver->a_level(id) < m_solver->maxLevel()) continue;

    id = m_bndCells[i].m_cellId;

    m_omega = 2.0 / (1.0 + 6.0 * m_solver->a_nu(id) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(id)));

    rhoF = 1.0; // standard value which is usually overwritten in the following. Only in exceptional case it remains 1

    //----------------------------------------------------------------------------------------------------
    // perform bounce back for distributions piercing the wall - distinquish between inner and outer cells
    //----------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // 1. boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      rhoF = m_solver->a_variable(id, PV->RHO);
      tmp2UF = 0.0;
      for(MInt d = 0; d < nDim; d++) {
        uF[d] = m_solver->a_variable(id, d);
        tmp2UF += uF[d] * uF[d];
        b[2 * d] = -uF[d];
        b[2 * d + 1] = uF[d];
      }

      for(MInt j = 0; j < nDist - 1; j++) {
        if(j < Ld::distFld(0)) {
          // dxB2 = cellHalfLength;
          tmpUF = b[j];
          tpIndex = 1;
        } else {
          if(j < (Ld::distFld(0) + Ld::distFld(1))) {
            // dxB2 = SQRT2 * cellHalfLength;
            k = j - Ld::distFld(0);
            tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          } else {
            // dxB2 = SQRT3 * cellHalfLength;
            k = j - (Ld::distFld(0) + Ld::distFld(1));
            tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
            tpIndex = 3;
          }
        }

        //---------------------------------------------------------------------------------
        // if q < 0.5, perform bounce back
        //---------------------------------------------------------------------------------

        if(m_bndCells[i].m_distances[j] < 0.5) {
          q = m_bndCells[i].m_distances[j]; // q = distance / dx; distance: distance between cell center and wall

          if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(j)) > 0) {
            nghbrId = m_solver->c_neighborId(id, Ld::oppositeDist(j));

            // extrapolate velocity from next inner cell
            uBF[0] = (q - 1.0) * m_solver->a_variable(id, PV->U)
                     + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->U) / (q + 1.0);
            uBF[1] = (q - 1.0) * m_solver->a_variable(id, PV->V)
                     + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->V) / (q + 1.0);

            IF_CONSTEXPR(nDim == 3)
            uBF[2] = (q - 1.0) * m_solver->a_variable(id, PV->W)
                     + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->W) / (q + 1.0);

            tmpUBF = 0.0;
            tmp2UBF = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              tmpUBF += (Ld::idFld(j, d) - 1) * uBF[d];
              tmp2UBF += uBF[d] * uBF[d];
            }
          } else {
            tmpUBF = 0.0;
            for(MInt dim = 0; dim < nDim; dim++) {
              uBF[dim] = (q - 1.0) * uF[dim] / q;
              tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
              tmp2UBF += uBF[dim] * uBF[dim];
            }
          }

          fEqW = Ld::tp(tpIndex) * rhoF
                 * (1.0 + tmpUBF * F1BCSsq + tmpUBF * tmpUBF * F1BCSsq * F1BCSsq * F1B2 - tmp2UBF * F1BCSsq * F1B2);
          fEqF = Ld::tp(tpIndex) * rhoF
                 * (1.0 + tmpUF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);

          // perform relaxation to extrapolated values
          m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) =
              (1.0 - m_omega) * (m_solver->a_distribution(id, Ld::oppositeDist(j)) - fEqF) + fEqW;

          // add to momentum sum
          distributions[Ld::oppositeDist(j)] +=
              m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) + m_solver->a_distribution(id, j);


        }

        // if the cell is located at a corner, there might be distributions missing
        else {
          // In this case simple bounce back is used. This could decrease the overall accuracy.
          if(m_solver->a_hasNeighbor(id, j) == 0) {
            m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) = m_solver->a_distribution(id, j);

            // add to momentum sum
            distributions[Ld::oppositeDist(j)] +=
                m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) + m_solver->a_distribution(id, j);
          }
        }

      } // end of loop over all directions

    }

    //---------------------------------------------------------------------------------
    // 2. boundary cell is outside fluid
    //---------------------------------------------------------------------------------

    // The outer cell does not belong to the flow field, thus the macroscopic values are held constant.
    // It is possible that there are cells without any cutting velocity! Those have to be controlled, too.
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        // if ( j < 6 )
        // 	dxB2 = cellHalfLength;
        // else {
        // 	if ( j < 18)
        // 	    dxB2 = SQRT2 * cellHalfLength;
        // 	else
        // 	    dxB2 = SQRT3 * cellHalfLength;
        // }

        //----------------------------------------------------------------------------------------------------------------
        // if q_innerNeighbor = 1 - q_bndCell > 0.5, perform bounceback on inner neighbor
        // if the inner neighbor is a halo cell, it can be skipped
        //----------------------------------------------------------------------------------------------------------------

        if(m_solver->a_hasNeighbor(id, j) > 0 && !m_solver->a_isHalo(m_solver->c_neighborId(id, j))) {
          if(m_bndCells[i].m_distances[j] <= 0.5) {
            q = F1 - m_bndCells[i].m_distances[j]; // q = 1 - 0.5*(distance/cellHalfLength);

            nghbrId = m_solver->c_neighborId(id, j);

            rhoF = m_solver->a_variable(nghbrId, PV->RHO);
            tmp2UF = 0.0;
            for(MInt d = 0; d < nDim; d++) {
              uF[d] = m_solver->a_variable(id, d);
              tmp2UF += uF[d] * uF[d];
              b[2 * d] = -uF[d];
              b[2 * d + 1] = uF[d];
            }

            if(j < Ld::distFld(0)) {
              tmpUF = b[Ld::oppositeDist(j)];
              tpIndex = 1;
            } else {
              if(j < (Ld::distFld(0) + Ld::distFld(1))) {
                k = Ld::oppositeDist(j) - Ld::distFld(0);
                tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
                tpIndex = 2;
              } else {
                k = Ld::oppositeDist(j) - (Ld::distFld(0) + Ld::distFld(1));
                tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
                tpIndex = 3;
              }
            }

            // 0.75 <= q
            if(q >= F3B4) {
              tmpUBF = 0.0;
              for(MInt dim = 0; dim < nDim; dim++) {
                uBF[dim] = (q - 1.0) * uF[dim] / q;
                tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
                tmp2UBF += uBF[dim] * uBF[dim];
              }

            }

            // 0.5 < q < 0.75
            else {
              if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(j)) > 0) {
                nghbrId = m_solver->c_neighborId(id, Ld::oppositeDist(j));

                // extrapolate velocity from next inner cell
                uBF[0] = (q - 1.0) * m_solver->a_variable(id, PV->U)
                         + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->U) / (q + 1.0);
                uBF[1] = (q - 1.0) * m_solver->a_variable(id, PV->V)
                         + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->V) / (q + 1.0);
                IF_CONSTEXPR(nDim == 3)
                uBF[2] = (q - 1.0) * m_solver->a_variable(id, PV->W)
                         + (1.0 - q) * (q - 1.0) * m_solver->a_variable(nghbrId, PV->W) / (q + 1.0);
                tmpUBF = 0.0;
                tmp2UBF = 0.0;
                for(MInt d = 0; d < nDim; d++) {
                  tmpUBF += (Ld::idFld(j, d) - 1) * uBF[d];
                  tmp2UBF += uBF[d] * uBF[d];
                }
              } else {
                tmpUBF = 0.0;
                for(MInt dim = 0; dim < nDim; dim++) {
                  uBF[dim] = (q - 1.0) * uF[dim] / q;
                  tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
                  tmp2UBF += uBF[dim] * uBF[dim];
                }
              }
            }

            fEqW = Ld::tp(tpIndex) * rhoF
                   * (1.0 + tmpUBF * F1BCSsq + tmpUBF * tmpUBF * F1BCSsq * F1BCSsq * F1B2 - tmp2UBF * F1BCSsq * F1B2);
            fEqF = Ld::tp(tpIndex) * rhoF
                   * (1.0 + tmpUF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);

            // perform relaxation to extrapolated values
            m_solver->a_oldDistribution(nghbrId, j) =
                (1.0 - m_omega) * (m_solver->a_distribution(nghbrId, j) - fEqF) + fEqW;

            // add to momentum sum
            distributions[j] +=
                m_solver->a_oldDistribution(nghbrId, j) + m_solver->a_distribution(nghbrId, Ld::oppositeDist(j));
            //		    distributions[j] += 2.0 * m_solver->a_distribution(neighborId, Ld::oppositeDist(j));
          }
        }

      } // end of loop over all directions

      // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
      // It is possible that there are cells without any cutting velocity! These are considered too.
      m_solver->a_variable(id, PV->RHO) = 1.0;
      for(MInt d = 0; d < nDim; d++) {
        m_solver->a_variable(id, d) = 0.0;
      }

      m_solver->a_oldVariable(id, PV->RHO) = 1.0;
      for(MInt d = 0; d < nDim; d++) {
        m_solver->a_oldVariable(id, d) = 0.0;
      }

      // Calculation of eq-distributions in bnd cell
      //--------------------------------------------
      // Calculation of distributions for directions with only one component
      for(MInt j = 0; j < Ld::distFld(0); j++) {
        m_solver->a_oldDistribution(id, j) = Ld::tp(1);
      }

      // Calculation of distributions for directions with two components
      tmpDistId = Ld::distFld(0);
      for(MInt j = 0; j < Ld::distFld(1); j++) {
        m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(2);
      }

      // Calculation of distributions for directions with three components
      tmpDistId = Ld::distFld(0) + Ld::distFld(1);
      for(MInt j = 0; j < Ld::distFld(2); j++) {
        m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(3);
      }

      // Calculation of distribution for rest particle distribution (center)
      m_solver->a_oldDistribution(id, Ld::lastId()) = Ld::tp(0);
    }

  } // end of loop over all m_bndcells
}

/** \LBBC{secLBBC_bc20020, bc20020, 2020}
 *
 * \author Andreas Lintermann
 * \date 24.02.2011
 *
 * Lattice Boltzmann enhanced no slip (bounce back) condition for
 * inclined walls. Bouzidi 2001 (aka "BFL rule" - see BC20000) including Thermal LBGK.
 *
 *
 * At the moment only the equilibrium distribution functions are evaluated for all cells having a cut. No BFL-rule is
 * applied to Thermal LBGK.
 *
 * \todo labels:LB,toenhance Consider using BFL-rule for Thermal LBGK for higher order wall-approxiomation
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20020(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20020(MInt index) {
  TRACE();

  const MFloat wT = getBoundaryTemperature(index);

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //      if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    interpolatedBounceBackSingleSpecies(i, uW.data());
    const MInt pCellId = m_bndCells[i].m_cellId;

    // Set constant temperature
    calculateEqDistsWallSingleSpeciesThermal(pCellId, wT);
  }
}

/** \LBBC{secLBBC_bc20023, bc20023, 2023}
 * \author Andreas Lintermann
 * \date 13.04.2011
 *
 * Lattice Boltzmann no slip (bounce back) condition
 *
 * If no neihbor exists in incoming direction, the incoming
 * distribution is reflected.
 * -> The wall is located at the cell center!
 * Must be performed after collision and BEFORE propagation.
 * This is an extension of BC 2002 for Thermal LBGK.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20023(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20023(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // leave out halo cells
    if(m_solver->a_isHalo(currentId)) continue;

    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(j)) == 0) {
        m_solver->a_oldDistribution(currentId, j) = m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j));
        m_solver->a_oldDistributionThermal(currentId, j) =
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j));
      }
    }
  }
}

/** \LBBC{secLBBC_bc20022, bc20022, 2022}
 *
 * \author Andreas Lintermann
 * \date 14.03.2011
 *
 * Lattice Boltzmann enhanced no slip (bounce back) condition for
 * inclined walls. Bouzidi 2001 (aka "BFL rule" - see BC20000) including Thermal LBGK.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20022(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20022(MInt index) {
  TRACE();

  //---------------------------------------------------------------------------------
  // go through m_bndcells with wallbndcnd;
  // the halo cells have to be included, since they may be of relevance in case 2
  //---------------------------------------------------------------------------------
  MInt p = 0;
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++, p++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // Set constant temperature
    //    m_solver->a_variable(currentId, PV->T) = m_solver->m_initTemperatureKelvin;
    // m_solver->a_oldVariable(currentId, PV->T) = m_solver->m_initTemperatureKelvin;

    // If there are several levels of boundary cells use only the highest
    // nearest neighbors, must be of the same level
    if(m_solver->a_level(currentId) < m_solver->maxLevel()) continue;

    //----------------------------------------------------------------------------------------------------
    // perform bounce back for distributions piercing the wall - distinquish between inner and outer cells
    //----------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // case 1: boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      for(MInt j = 0; j < nDist - 1; j++) {
        // if ( j < 6 )
        // 	dxB2 = cellHalfLength;
        // else {
        // 	if ( j < 18)
        // 	    dxB2 = SQRT2 * cellHalfLength;
        // 	else
        // 	    dxB2 = SQRT3 * cellHalfLength;
        // }

        //---------------------------------------------------------------------------------
        // if q <= 0.5, perform bounce back
        //---------------------------------------------------------------------------------

        if(m_bndCells[i].m_distances[j] <= 0.5) {
          const MFloat F2q =
              F2 * m_bndCells[i].m_distances[j]; // q = distance / dx; distance: distance between cell center and wall
          if(m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(j)) > 0) {
            m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j)) =
                (F2q * m_solver->a_distribution(currentId, j))
                + ((1 - F2q) * m_solver->a_oldDistribution(currentId, j));
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j)) =
                (F2q * m_solver->a_distributionThermal(currentId, j))
                + ((1 - F2q) * m_solver->a_oldDistributionThermal(currentId, j));
          }
          // If opposite neighbors are missing, an incoming distribution has to be extrapolated from another cell inside
          // the fluid This could decrease the overall accuracy
          else {
            // do simple bounce back
            m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j)) = m_solver->a_distribution(currentId, j);
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j)) =
                m_solver->a_distributionThermal(currentId, j);
          }
        }
        // if no cut was found and there is no neighbor
        // then cell is located at a corner.
        else {
          // In this case simple bounce back is used. This could decrease the overall accuracy.
          if(m_solver->a_hasNeighbor(currentId, j) == 0) {
            m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j)) = m_solver->a_distribution(currentId, j);
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j)) =
                m_solver->a_distributionThermal(currentId, j);
          }
        }
      } // end of loop over all directions
    }

    //---------------------------------------------------------------------------------
    // case 2: boundary cell is outside fluid
    //---------------------------------------------------------------------------------

    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        // if ( j < 6 )
        // 	dxB2 = cellHalfLength;
        // else {
        // 	if ( j < 18)
        // 	    dxB2 = SQRT2 * cellHalfLength;
        // 	else
        // 	    dxB2 = SQRT3 * cellHalfLength;
        // }

        //----------------------------------------------------------------------------------------------------------------
        // if q_innerNeighbor = 1 - q_bndCell > 0.5, perform bounceback on inner neighbor
        // if the inner neighbor is a halo cell, it can be skipped
        //----------------------------------------------------------------------------------------------------------------
        //   	    if ( m_solver->a_hasNeighbor(currentId, j) > 0 && m_solver->c_neighborId(currentId, j) <
        //   m_solver->noInternalCells()
        //   )
        if(m_solver->a_hasNeighbor(currentId, j) > 0 && !m_solver->a_isHalo(m_solver->c_neighborId(currentId, j))) {
          const MInt neighborId = m_solver->c_neighborId(currentId, j);

          const MFloat F2q = F2 - F2 * m_bndCells[i].m_distances[j]; // q = 1 - 0.5*(distance/cellHalfLength);

          if(m_bndCells[i].m_distances[j] < 0.5) {
            m_solver->a_oldDistribution(neighborId, j) =
                (m_solver->a_distribution(neighborId, Ld::oppositeDist(j)) / F2q)
                + (m_solver->a_distribution(neighborId, j) * (F2q - 1.0) / F2q);
            m_solver->a_oldDistributionThermal(neighborId, j) =
                (m_solver->a_distributionThermal(neighborId, Ld::oppositeDist(j)) / F2q)
                + (m_solver->a_distributionThermal(neighborId, j) * (F2q - 1.0) / F2q);
          }
        }

      } // end of loop over all directions

      // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
      // It is possible that there are cells without any cutting velocity! These are considered too.

      m_solver->a_variable(currentId, PV->RHO) = 1.0;
      m_solver->a_oldVariable(currentId, PV->RHO) = 1.0;
      for(MInt dir = 0; dir < nDim; dir++) {
        m_solver->a_variable(currentId, PV->U + dir) = 0.0;
        m_solver->a_oldVariable(currentId, PV->U + dir) = 0.0;
      }
      m_solver->a_variable(currentId, PV->T) = m_solver->m_initTemperatureKelvin;
      m_solver->a_oldVariable(currentId, PV->T) = m_solver->m_initTemperatureKelvin;

      const MFloat T = m_solver->a_variable(currentId, PV->T);

      // Calculation of eq-distributions in bnd cell
      //--------------------------------------------
      // Calculation of distributions for directions with only one component
      for(MInt j = 0; j < Ld::distFld(0); j++) {
        m_solver->a_oldDistribution(currentId, j) = Ld::tp(1);
        // new
        m_solver->a_oldDistributionThermal(currentId, j) = Ld::tp(1) * F1BCSsq * T * CSsq;
      }

      // Calculation of distributions for directions with two components
      MInt tmpDistId = Ld::distFld(0);
      for(MInt j = 0; j < Ld::distFld(1); j++) {
        m_solver->a_oldDistribution(currentId, tmpDistId + j) = Ld::tp(2);
        // m_solver->a_oldDistributionThermal(currentId, j) = Ld::tp(1) * F1BCSsq * T * CSsq;
        m_solver->a_oldDistributionThermal(currentId, tmpDistId + j) = Ld::tp(2) * F1BCSsq * T * CSsq;
      }

      // Calculation of distributions for directions with three components
      tmpDistId = Ld::distFld(0) + Ld::distFld(1);
      for(MInt j = 0; j < Ld::distFld(2); j++) {
        m_solver->a_oldDistribution(currentId, tmpDistId + j) = Ld::tp(3);
        // m_solver->a_oldDistributionThermal(currentId, tmpDistId + j) = Ld::tp(2) * F1BCSsq * T * CSsq;
        m_solver->a_oldDistributionThermal(currentId, tmpDistId + j) = Ld::tp(3) * F1BCSsq * T * CSsq;
      }

      // Calculation of distribution for rest particle distribution (center)
      m_solver->a_oldDistribution(currentId, Ld::lastId()) = Ld::tp(0);
      m_solver->a_oldDistributionThermal(currentId, Ld::lastId()) = Ld::tp(0) * T;
    }
  } // end of loop over all m_bndcells
}

/** \LBBC{secLBBC_bc20024, bc20024, 2024}
 * \author Andreas Lintermann
 * \date 13.04.2011
 *
 * Lattice Boltzmann no slip (bounce back) condition
 *
 * no slip condition for periodic channel
 * pressure is extrapolated
 * velocity is set to zero
 * Similar to BC2003 including Thermal LBGK
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20024(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20024(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    // if( (globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // extrapolate inner density values
    const MInt neighborDir = (m_solver->a_hasNeighbor(currentId, 2) == 0) ? 3 : 2;
    const MInt neighborId = m_solver->c_neighborId(currentId, neighborDir);

    // get old values
    const MFloat old_rho = m_solver->a_oldVariable(currentId, PV->RHO);
    const MFloat old_T = m_solver->m_initTemperatureKelvin;

    // non-reflecting bc (Finck, Haenel 2008)
    // rho = ( old_rho + F1BCS * ( sqrt(tmp2) - sqrt(old_tmp2) ) + m_rho1 ) / 2.0 ;
    const MFloat rho = F1B2 * (m_solver->a_variable(neighborId, PV->RHO) + old_rho);
    const MFloat T = F1B2 * (m_solver->a_variable(neighborId, PV->T) + old_T);


    // Calculation of eq-distributions in bnd cell
    //--------------------------------------------

    // Calculation of distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      m_solver->a_oldDistribution(currentId, j) = Ld::tp(1) * F1BCSsq * rho * CSsq;
      m_solver->a_distribution(currentId, j) = m_solver->a_oldDistribution(currentId, j);

      m_solver->a_oldDistributionThermal(currentId, j) = Ld::tp(1) * F1BCSsq * T * CSsq;
      m_solver->a_distributionThermal(currentId, j) = m_solver->a_oldDistributionThermal(currentId, j);
    }

    // Calculation of distributions for directions with two components
    MInt tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      m_solver->a_oldDistribution(currentId, tmpDistId + j) = Ld::tp(2) * F1BCSsq * rho * CSsq;
      m_solver->a_distribution(currentId, tmpDistId + j) = m_solver->a_oldDistribution(currentId, tmpDistId + j);

      m_solver->a_oldDistributionThermal(currentId, tmpDistId + j) = Ld::tp(2) * F1BCSsq * T * CSsq;
      m_solver->a_distributionThermal(currentId, tmpDistId + j) =
          m_solver->a_oldDistributionThermal(currentId, tmpDistId + j);
    }

    // Calculation of distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      m_solver->a_oldDistribution(currentId, tmpDistId + j) = Ld::tp(3) * F1BCSsq * rho * CSsq;
      m_solver->a_distribution(currentId, tmpDistId + j) = m_solver->a_oldDistribution(currentId, tmpDistId + j);

      m_solver->a_oldDistributionThermal(currentId, tmpDistId + j) = Ld::tp(3) * F1BCSsq * T * CSsq;
      m_solver->a_distributionThermal(currentId, tmpDistId + j) =
          m_solver->a_oldDistributionThermal(currentId, tmpDistId + j);
    }

    // Calculation of distribution for rest particle distribution (center)
    m_solver->a_oldDistribution(currentId, Ld::lastId()) = Ld::tp(0) * rho;
    m_solver->a_distribution(currentId, Ld::lastId()) = m_solver->a_oldDistribution(currentId, Ld::lastId());

    m_solver->a_oldDistributionThermal(currentId, Ld::lastId()) = Ld::tp(0) * T;
    m_solver->a_distributionThermal(currentId, Ld::lastId()) =
        m_solver->a_oldDistributionThermal(currentId, Ld::lastId());

    m_solver->a_variable(currentId, PV->RHO) = rho;
    for(MInt dir = 0; dir < nDim; dir++) {
      m_solver->a_variable(currentId, PV->U + dir) = 0.0;
    }
    m_solver->a_variable(currentId, PV->T) = T;
  }
}

/** \LBBC{secLBBC_bc20025, bc20025, 2025}
 * \author Andreas Lintermann
 * \date 13.04.2011
 *
 * Wall BC for thermal flows
 *
 * Sets the equilibrum distribution functions for given macroscopic variables.<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20025(MInt index)
 */

template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20025(MInt index) {
  TRACE();

  MInt id;
  MInt tmpDistId;

  //---------------------------------------------------------------------------------
  // go through m_bndcells with wallbndcnd;
  // the halo cells have to be included, since they may be of relevance in case 2
  //---------------------------------------------------------------------------------
  MInt p = 0;
  // MFloat min_x=16.0;
  // MFloat min_z=7.0;
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++, p++) {
    //    if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    id = m_bndCells[i].m_cellId;

    // Set constant temperature
    m_solver->a_variable(id, PV->T) = m_solver->m_initTemperatureKelvin;
    m_solver->a_oldVariable(id, PV->T) = m_solver->m_initTemperatureKelvin;

    for(MInt d = 0; d < nDim; d++) {
      m_solver->a_variable(id, d) = 0.0;
    }

    MFloat T = m_solver->a_variable(id, PV->T);
    MFloat rho = m_solver->a_variable(id, PV->RHO);

    for(MInt j = 0; j < Ld::distFld(0); j++) {
      m_solver->a_oldDistribution(id, j) = Ld::tp(1) * F1BCSsq * rho * CSsq;
      m_solver->a_oldDistributionThermal(id, j) = Ld::tp(1) * F1BCSsq * T * CSsq;
    }

    // Calculation of eq-distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(2) * F1BCSsq * rho * CSsq;
      m_solver->a_oldDistributionThermal(id, tmpDistId + j) = Ld::tp(2) * F1BCSsq * T * CSsq;
    }

    // Calculation of eq-distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      m_solver->a_oldDistribution(id, tmpDistId + j) = Ld::tp(3) * F1BCSsq * rho * CSsq;
      m_solver->a_oldDistributionThermal(id, tmpDistId + j) = Ld::tp(3) * F1BCSsq * T * CSsq;
    }

    m_solver->a_oldDistribution(id, Ld::lastId()) = Ld::tp(0) * rho;
    m_solver->a_oldDistributionThermal(id, Ld::lastId()) = Ld::tp(0) * T;
  }
}


/** \LBBC{secLBBC_bc20026, bc20026, 2026}
 * \author Andreas Lintermann
 * \date 28.04.2011
 *
 * Lattice Boltzmann no slip condition also for Thermal LBGK
 *
 * If no neihbor exists in outgoing direction, the outgoing
 * distribution is reflected. The wall is located at the outer cell edge. Must be performed after propagation.
 *
 * The temperture of the wall sided cells is evaluated based on the current temperature and the wall temperature.
 * A ghost layer is initialized holding the inner wall temperatures. The eq-dist functions are calculated and a bounce
 * back is performed for the temperature.
 *
 * \todo labels:LB Halo layer for all boundary conditions.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20026(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20026(MInt index) {
  TRACE();

  MInt tmpDistId;

  constexpr MInt lb_proj_zplane_neigh[27] = {-1, -1, -1, -1, 26, -1, -1, -1, -1, -1, 0, -1, 1, -1,
                                             2,  -1, 3,  -1, 6,  -1, 7,  -1, 8,  -1, 9, -1, -1};

  MInt nocells = m_bndCndOffsets[index + 1] - m_bndCndOffsets[index];
  ScratchSpace<MFloat> ghosts(nocells, nDist, AT_, "ghosts");

  MInt bndId = 0;
  MFloat T_inner = 0.0, T_outer = 0.0, T_wall = m_solver->m_initTemperatureKelvin;

  for(MInt num = 0; num < nocells; num++) {
    bndId = num + m_bndCndOffsets[index];
    T_inner = m_solver->a_variable(m_bndCells[bndId].m_cellId, PV->T);

    T_outer = 2 * T_wall - T_inner;

    for(MInt j = 0; j < Ld::distFld(0); j++)
      ghosts(num, j) = Ld::tp(1) * F1BCSsq * T_outer * CSsq;

    // Calculation of eq-distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++)
      ghosts(num, tmpDistId + j) = Ld::tp(2) * F1BCSsq * T_outer * CSsq;

    // Calculation of eq-distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++)
      ghosts(num, tmpDistId + j) = Ld::tp(3) * F1BCSsq * T_outer * CSsq;

    ghosts(num, Ld::lastId()) = Ld::tp(0) * T_outer;
  }

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // leave out halo cells
    //    if (m_bndCells[i].m_cellId >= m_solver->noInternalCells())
    //	continue;

    if(m_solver->a_isHalo(m_bndCells[i].m_cellId)) continue;


    MInt proj_dir = 0;
    MInt proj_neigh = 0;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(m_bndCells[i].m_cellId, Ld::oppositeDist(j)) == 0) {
        m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) =
            m_solver->a_distribution(m_bndCells[i].m_cellId, Ld::oppositeDist(j));

        proj_dir = lb_proj_zplane_neigh[Ld::oppositeDist(j)];

        if(proj_dir < 0) continue;

        if(proj_dir == 26)
          proj_neigh = i - m_bndCndOffsets[index];
        else {
          for(MInt l = m_bndCndOffsets[index]; l < m_bndCndOffsets[index + 1]; l++)
            if(m_bndCells[l].m_cellId == m_solver->c_neighborId(m_bndCells[i].m_cellId, proj_dir)) {
              proj_neigh = l - m_bndCndOffsets[index];
              break;
            }
        }


        if(proj_neigh >= 0 && proj_neigh < nocells)
          m_solver->a_oldDistributionThermal(m_bndCells[i].m_cellId, j) = ghosts(proj_neigh, j);
      }
    }
  }


  //  for(MInt i = 0; i < nocells; i++) delete[] ghosts[i];
  //  delete[] ghosts;
}


/** \LBBC{secLBBC_bc20027, bc20027, 2027}
 * \author Andreas Lintermann
 * \date 21.05.2013
 *
 * Lattice Boltzmann no slip interpolated bounce-back condition also for Thermal LBGK
 *
 * This BC applies BC 2000 for the velocity and the density. In contrast to BC 2026 this BC allows
 * the definition of an interpolated temperature dependent on the distance to the wall. The
 * algorithm is divided into two parts:
 *
 * 1. do an interpolated bounce back as done in BC 2000 for the velocity and the density distributions
 * 2. do an interpolated prescription of the PPDFs for the temperature:
 *    - case 1: cell center is inside fluid:
 *      - if a neighbor does not exist overwrite the opposite direction by a distribution, which
 *        is calculated from the eq-distributions of a temperature previously extrapolated based
 *        on the distance to the wall
 *    - case 2: cell center is outside fluid:
 *      - calculate extrapolated temperature of outer cell based on wall-temperature and
 *        inner temperature
 *      - calculate eq-distributions based on this temperature and propagate to neighboring
 *        fluid cells
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20027(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20027(MInt index) {
  TRACE();

  MInt start = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];

  MFloat T_wall = m_solver->m_initTemperatureKelvin;
  MFloat eps = 0.0000000000001;

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

  for(MInt i = start; i < end; i++) {
    //      if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    // 1. do an interpolated bounce back as done in BC 2000 for the velocity and the
    //    density distributions
    interpolatedBounceBackSingleSpecies(i, uW.data());

    MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    // 2. do an interpolated prescription of the PPDFs for the temperature

    //---------------------------------------------------------------------------------
    // case 1: boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      for(MInt j = 0; j < nDist - 1; j++) {
        // if we do not have a neighbor in the current direction
        if(m_solver->a_hasNeighbor(pCellId, j) == 0) {
          MFloat T_inner = m_solver->a_variable(currentId, PV->T);
          MFloat T_outer = 0.0;
          if(m_bndCells[i].m_distances[j] > eps) {
            T_outer = (T_wall - T_inner) / (m_bndCells[i].m_distances[j]) + T_inner;
            // Calculate single eq-distribution in direction j and overwrite neighbor
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j)) =
                Ld::tp(Ld::distType(j)) * F1BCSsq * T_outer * CSsq;

          } else {
            // Calculate single eq-distribution in direction j and overwrite neighbor
            m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(j)) =
                Ld::tp(Ld::distType(j)) * F1BCSsq * T_wall * CSsq;
          }
        }
      }
    }

    //---------------------------------------------------------------------------------
    // case 2: boundary cell is outside fluid
    //---------------------------------------------------------------------------------
    else {
      for(MInt j = 0; j < nDist; j++) {
        // if we have a neighbor
        if(m_solver->a_hasNeighbor(pCellId, j) != 0)
          if(!m_solver->a_isBndryCell(m_solver->c_neighborId(pCellId, j))) {
            MInt neighborId = m_solver->c_neighborId(pCellId, j);

            MFloat T_inner = m_solver->a_variable(neighborId, PV->T);
            MFloat T_outer = (T_wall - T_inner) / (1.0 - m_bndCells[i].m_distances[j]) + T_inner;

            // Calculate single eq-distribution in direction j and overwrite neighbor
            m_solver->a_oldDistributionThermal(neighborId, j) = Ld::tp(Ld::distType(j)) * F1BCSsq * T_outer * CSsq;
          }
      }
    }
  }
}

/** \LBBC{secLBBC_slidingWall, slidingWall, 205x}
 *
 * Haenel wall bc with non zero velocity
 * planar wall moving at v=Ma*c_s along axis-direction
 *
 * \tparam direction Wall orientation in Cartesian direction
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::slidingWall(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt direction>
void LbBndCndDxQy<nDim, nDist, SysEqn>::slidingWall(MInt index) {
  TRACE();

  // constexpr sanity check
  if constexpr(direction > nDim - 1) {
    std::stringstream ss;
    ss << "Function not defined for nDim=" << nDim << " and direction =" << direction << std::endl;
    TERMM(1, ss.str());
    return;
  }

  // set wall velocity
  std::array<MFloat, nDim> uW;
  std::array<MFloat, 2 * nDim> c;
  uW.fill(0.0);
  uW[direction / 2] = m_Ma * LBCS;
  for(MInt d = 0; d < nDim; d++) {
    c[2 * d] = -uW[d];
    c[2 * d + 1] = uW[d];
  }

  // go through m_bndcells with wallbndcnd
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt id = m_bndCells[i].m_cellId;

    // leave out halo cells
    if(m_solver->a_isHalo(id)) continue;

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    //  neighborId = m_solver->c_neighborId(id,  Ld::oppositeDist(direction) );

    //     uW[0] = F2B3 * m_solver->a_variable(id, PV->U) - F1B2 * m_solver->a_variable(neighborId, PV->U);
    //     uW[1] = F2B3 * m_solver->a_variable(id, PV->V) - F1B2 * m_solver->a_variable(neighborId, PV->V);
    //     uW[2] = F2B3 * m_solver->a_variable(id, PV->W) - F1B2 * m_solver->a_variable(neighborId, PV->W);

    const MFloat rhoF = m_solver->a_variable(id, PV->RHO);

    /*
        for( MInt j = 0; j < Ld::dxQyFld(); j++){

      currentDist = Ld::componentFld( Ld::oppositeDist(direction) , j);

      if ( currentDist < 6 ) {
          tmpUW = c[currentDist];
          tpIndex = 1;
      }
      else {
          if ( currentDist < 18){
        k = currentDist - Ld::distFld(0);
        tmpUW = (c[Ld::mFld1(2*k)] + c[Ld::mFld1(2*k+1)]);
        tpIndex = 2;
          }
          else {
        k = currentDist - (Ld::distFld(0) + Ld::distFld(1));
        tmpUW = (c[Ld::mFld2(3*k)] + c[Ld::mFld2(3*k+1)] + c[Ld::mFld2(3*k+2)] );
        tpIndex = 3;
          }
      }

      // 1. bounce back part
    //	if(m_solver->a_hasNeighbor(id, currentDist) > 0){
          m_solver->a_oldDistribution(id, currentDist) = m_solver->a_distribution(id, Ld::oppositeDist(currentDist));
    //	}
    // 	else{//if cell is at corner, extrapolate distribution from inner cell
    // 	    for (MInt k=0; k < nDist-1; k++){
    // 		if(m_solver->a_hasNeighbor(id, k) > 0 && !m_solver->a_onlyBoundary(m_solver->c_neighborId(id, k)) )
    // 		    m_solver->a_oldDistribution(id, currentDist) = m_solver->a_oldDistribution(m_solver->c_neighborId(id, k),
    Ld::oppositeDist(currentDist));
    // 	    }
    // 	}

      // 2. add momentum
      m_solver->a_oldDistribution(id, currentDist) += 2.0 * Ld::tp(tpIndex) * rhoF * F1BCSsq * tmpUW;// the
    minus corresponds to the opposite direction

        }
    */
    for(MInt currentDist = 0; currentDist < nDist - 1; currentDist++) {
      // 1. bounce back part
      if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(currentDist)) == 0) {
        m_solver->a_oldDistribution(id, currentDist) = m_solver->a_distribution(id, Ld::oppositeDist(currentDist));
      }
      // 	else{//if cell is at corner, extrapolate distribution from inner cell
      // 	    for (MInt k=0; k < nDist-1; k++){
      // 		if(m_solver->a_hasNeighbor(id, k) > 0 && !m_solver->a_onlyBoundary(m_solver->c_neighborId(id, k)) )
      // 		    m_solver->a_oldDistribution(id, currentDist) = m_solver->a_oldDistribution(m_solver->c_neighborId(id,
      // k), Ld::oppositeDist(currentDist));
      // 	    }
      // 	}

      // 2. add momentum
      MFloat tmpUW;
      if(currentDist < Ld::distFld(0)) {
        tmpUW = c[currentDist];
      } else if(currentDist < Ld::distFld(0) + Ld::distFld(1)) {
        const MInt k = currentDist - Ld::distFld(0);
        tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
      } else {
        const MInt k = currentDist - (Ld::distFld(0) + Ld::distFld(1));
        tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
      }
      m_solver->a_oldDistribution(id, currentDist) += 2.0 * Ld::tp(Ld::distType(currentDist)) * rhoF * F1BCSsq * tmpUW;
    }


  } // end of loop over all m_bndcells
}


/** \LBBC{secLBBC_bc20220, bc20220, 2220}
 * \author Moritz Waldmann
 * \date 10.12.2019
 *
 * Lattice Boltzmann enhanced no slip (bounce back) condition for
 * inclined walls. Bouzidi 2001 (aka "BFL rule" - see BC20000) including Thermal LBGK.
 *
 * At the moment only the equilibrium distribution functions are evaluated for all cells having a cut. No BFL-rule is
 * applied to Thermal LBGK.
 *
 * \todo labels:LB,toenhance Consider using BFL-rule for Thermal LBGK for higher order wall-approxiomation
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20220(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20220(MInt index) {
  TRACE();

  const MFloat wT = getBoundaryTemperature(index);

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    interpolatedBounceBackSingleSpecies(i, uW.data());

    interpolatedBounceBackSingleSpeciesThermal(i, wT, uW.data());
  });
}

/** \brief Lattice Boltzmann enhanced no slip (bounce back) condition for
 *        inclined walls with the Dirichlet condition for the thermal LBGK.
 *        Mode = 0 is a 2D/3D cylinder flow.
 *        Mode = 1 is a 2D/3D channel flow.
 *        Mode = 2 is a 2D channel flow.
 *
 * \author Shota Ito
 * \date 11.09.2022
 **
 * \param[in] index Boundary index
 *            bcMode type of bc (see above)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::heatFluxBc(MInt index, MInt bcMode) {
  TRACE();

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MFloat wT = F0;

    const MInt plugFlowDir = (bcMode == 0 && nDim == 3) ? 2 : 0;
    const MBool restingFluid = (bcMode == 0 && nDim == 2);
    const MInt pCellId = m_bndCells[i].m_cellId;

    if(bcMode == 0) {
      if(nDim == 3) {
        // wave number and heat flux at the wall
        MFloat k = F2 * M_PI / (m_referenceLength * 2.0);

        // since it is a rhsBnd it should be performed after propagation
        MFloat dz = m_solver->c_cellLengthAtLevel(m_solver->a_level(m_bndCells[i].m_cellId));
        MFloat z = m_solver->a_coordinate(m_bndCells[i].m_cellId, 2) / dz;
        z = z - m_referenceLength;

        // 3D Pipe flow
        wT = cos(k * z);
      } else {
        // since it is a rhsBnd it should be performed after propagation
        MFloat x = m_solver->a_coordinate(m_bndCells[i].m_cellId, 0);
        MFloat y = m_solver->a_coordinate(m_bndCells[i].m_cellId, 1);

        MFloat phi = atan(y / x);
        MFloat n = 2.0;
        wT = cos(n * phi);
      }
      // perform bounce back
      interpolatedBounceBackSingleSpeciesThermal(i, wT, uW.data());
    } else if(bcMode == 1) {
      // wave number and heat flux at the wall
      MFloat k = F2 * M_PI / (m_referenceLength);

      // since it is a rhsBnd it should be performed after propagation
      MFloat dx = m_solver->c_cellLengthAtLevel(m_solver->a_level(m_bndCells[i].m_cellId));
      MFloat x = m_solver->a_coordinate(m_bndCells[i].m_cellId, 0) / dx;

      // channel flow
      // heat flux in y direction
      wT = cos(k * x);

      // perform bounce back
      interpolatedBounceBackSingleSpeciesThermal(i, wT, uW.data());
    } else if(bcMode == 2) {
      // wave number and heat flux at the wall
      MFloat k = F2 * M_PI / m_referenceLength;
      // since it is a rhsBnd it should be performed after propagation
      MFloat dx = m_solver->c_cellLengthAtLevel(m_solver->a_level(m_bndCells[i].m_cellId));
      MFloat x = m_solver->a_coordinate(m_bndCells[i].m_cellId, 0) / dx;

      // channel flow
      // heat flux in y direction
      MFloat wq = (m_solver->m_kappa / m_referenceLength) * cos(k * x);

      // perform bounce back
      // TODO: So far only valid for 2D simulations
      ASSERT(nDim == 2, "Not implemented for 3D channels yet.");
      interpolatedBounceBackSingleSpeciesThermalFlux(i, wq, index);
    }

    // Set the density to 1.0 and the velocity to obtain a plug flow
    MFloat l_rho = 1.0;
    MFloat l_uu[nDim];
    for(MInt d = 0; d < nDim; d++)
      l_uu[d] = F0;

    if(!restingFluid) l_uu[plugFlowDir] = m_Ma * LBCS;

    MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    for(MInt n = 0; n < nDim; n++)
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
  }
}

/** \brief Lattice Boltzmann boundary condition for modelling nasal mucosa layer.
 * Two ideas are implemented:
 * 1. Additional artificial sublayer is created which models the heat conduction
 * inside the mucous layer with a constant organ-side temperature.
 * 2. Latent heat of evaporation is introduced to model cooling/heating effects
 * depending on the humidity flux across the boundary.
 *
 * 	This is a coupled Dirichlet-type boundary condition for the scalar transport
 * 	of the humidity concentration and temperature. The Dirichlet condition is
 * 	modelled by the idea of the interpolated bounce back for scalar transport
 * 	problems as introduced in the paper Li et al. (2013).
 *
 * 	This version of the boundary condition calculates the latent heat term with the
 * 	temperature of the previous time step.
 *
 * \author Shota Ito
 * \date 28.09.2022
 **
 * \param[in] index Boundary index
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20501(MInt index) {
  TRACE();

  std::array<MFloat, nDim> uW{};
  getBoundaryVelocity(index, uW.data());

  // loop over cells and call coupled boundary condition
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    // -----------------------------------------------------
    // 1. Calculate wall concentration and wall temperature
    // -----------------------------------------------------
    MFloat Cw[nDist] = {F0};
    MFloat Tw[nDist] = {F0};

    calculateWallInterface(i, Cw, Tw);

    // ------------------------------------
    // 2. Perform interpolated bounce back
    // ------------------------------------

    if(!m_mucosaModel.plugFlow) {
      // perform interpolated bounce back
      interpolatedBounceBackSingleSpecies(i, uW.data()); // no-slip velocity condition
    } else {
      // Set the density to 1.0 and the velocity to obtain a plug flow
      const MInt pCellId = m_bndCells[i].m_cellId;
      MFloat l_rho = 1.0;
      MFloat l_uu[nDim] = {0.0};
      if(nDim == 3) l_uu[2] = m_Ma * LBCS;
      MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);
      m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);

      m_solver->a_variable(pCellId, PV->RHO) = l_rho;
      for(MInt n = 0; n < nDim; n++)
        m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }

    interpolatedBounceBackSingleSpeciesTransport(i, Cw, uW.data());
    interpolatedBounceBackSingleSpeciesThermal(i, Tw, uW.data());
  }
}

/** \brief Initialize bc20501.
 *
 * Initialize bc20501. The parameters of the nasal mucosa model can be defined/adjusted here.
 * They are used to transform the simulation variables into dimensionalized variables in the
 * calculateWallInterface function. The fluid-side and mucosa-side distances are also
 * calculated by the calculateSublayerDistances function.
 *
 * \author Shota Ito
 * \date 07.09.2022
 **
 * \param[in] index Boundary index
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20501_init(MInt index) {
  TRACE();

  /*! \page propertyPage1
    \section m_plugFlow
    <code>MInt LbBndCndDxQy::plugFlow</code>\n
    Initiates a plugFlow. Used for testing and validation
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_mucosaModel.plugFlow = Context::getSolverProperty<MInt>("plugFlow", m_solverId, AT_);

  /*! \page propertyPage1
    \section referenceLengthSegId
    <code>MInt LbBndCndDxQy::segmentId</code>\n
    This property is mandatory if property referenceLengthLB is not given.
    Defines the segment id to be used for the calculation of the characteristic
    length in cell units.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MInt segmentId = Context::getSolverProperty<MInt>("charSegId", m_solverId, AT_);

  // Get the boundary points of the segment
  MInt num = 0;
  MFloat** bndVs = m_solver->m_geometry->GetBoundaryVertices(segmentId, nullptr, nullptr, 0, &num);

  // Calculate circumference
  MFloat circ = m_solver->m_geometry->calcCircumference(bndVs, num);

  // Get size of surface
  MFloat size = m_solver->m_geometry->GetBoundarySize(segmentId);

  // Get the hydraulic diameter
  MFloat tracheaDiameter = 4 * size / circ;

  // define reference variables for conversion between dimensioned and dimensionless
  // length of one cell in relation to the referenceLength in meter
  m_mucosaModel.refDx = tracheaDiameter / m_referenceLength;

  /*! \page propertyPage1
    \section conductivity
    <code>MInt LbBndCndDxQy::conductivity</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional conductivity
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat conductivity = 0.025;
  conductivity = Context::getSolverProperty<MFloat>("conductivity", m_solverId, AT_, &conductivity);

  /*! \page propertyPage1
    \section density
    <code>MInt LbBndCndDxQy::density</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional density
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat density = 1.225;
  density = Context::getSolverProperty<MFloat>("density", m_solverId, AT_, &density);

  /*! \page propertyPage1
    \section viscosity
    <code>MInt LbBndCndDxQy::viscosity</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional viscosity
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat viscosity = 16.45 * 1e-6;
  viscosity = Context::getSolverProperty<MFloat>("viscosity", m_solverId, AT_, &viscosity);

  /*! \page propertyPage1
    \section mucosaConductivity
    <code>MInt LbBndCndDxQy::mucosaConductivity</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional mucosaConductivity
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat mucosaCond = 0.6;
  mucosaCond = Context::getSolverProperty<MFloat>("mucosaCond", m_solverId, AT_, &mucosaCond);

  /*! \page propertyPage1
    \section organSideTemp
    <code>MInt LbBndCndDxQy::organSideTemp</code>\n
    default = <code>0</code>\n\n
    Provides the temperature of the interior of the mucosa in Celsius
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat organSideTemp = 37.0;
  organSideTemp = Context::getSolverProperty<MFloat>("organSideTemp", m_solverId, AT_, &organSideTemp);

  /*! \page propertyPage1
    \section mucosaThickness
    <code>MInt LbBndCndDxQy::mucosaThickness</code>\n
    default = <code>0</code>\n\n
    Provides the thickness of the mucosa layer in meter
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat mucosaThickness = 0.001;
  mucosaThickness = Context::getSolverProperty<MFloat>("mucosaThickness", m_solverId, AT_, &mucosaThickness);

  /*! \page propertyPage1
    \section mucosaDiffusivity
    <code>MInt LbBndCndDxQy::mucosaDiff</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional diffusivity of the mucosa layer
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat mucosaDiff = 2.6 * 1e-5;
  mucosaDiff = Context::getSolverProperty<MFloat>("mucosaDiffusivity", m_solverId, AT_, &mucosaDiff);

  /*! \page propertyPage1
    \section bndryDiffusivity
    <code>MInt LbBndCndDxQy::bndryDiff</code>\n
    default = <code>0</code>\n\n
    Provides the dimensional diffusivity of the bndry layer
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  MFloat bndryDiff = 3.0 * 1e-5;
  bndryDiff = Context::getSolverProperty<MFloat>("bndryDiffusivity", m_solverId, AT_, &bndryDiff);

  // reference velocity in m/s
  MFloat refU = (m_solver->m_Re * viscosity) / tracheaDiameter;

  // reference temperature in celsius
  m_mucosaModel.refT = organSideTemp;

  // reference water fraction
  // Inthavong et al (2018) "Wet surface wall model for latent heat exchange during evaporation"
  MFloat T = m_mucosaModel.refT;
  MFloat refWaterFrac = (2.027 + 0.6036 * T - 0.010972 * T * T + 0.0006312 * T * T * T) * 1e-3;

  // reference concentration
  m_mucosaModel.refC = refWaterFrac * density;

  // reference diffusivity
  m_mucosaModel.refDiff = (refU * m_mucosaModel.refDx) / (m_Ma * LBCS);

  // reference conductivity
  m_mucosaModel.refCondF = conductivity;

  // concentration in the interior of the mucosa
  MFloat organSideConc = m_mucosaModel.refC;

  // non-dimensional parameters of the nasal cavity flow
  m_mucosaModel.mucosaThickness = mucosaThickness / m_mucosaModel.refDx;
  m_mucosaModel.T_o = organSideTemp / m_mucosaModel.refT;
  m_mucosaModel.C_o = organSideConc / m_mucosaModel.refC;
  m_mucosaModel.condRatio = mucosaCond / conductivity;
  m_mucosaModel.diffRatio = mucosaDiff / bndryDiff;

  // initialize an array for the old wall temperatures
  mAlloc(m_oldWallTemp, m_solver->a_noCells(), nDist, "oldWallTemp", m_mucosaModel.T_o, AT_);

  // -----------------------------------------------
  //  Calculate distances inside the mucousLayer
  // -----------------------------------------------

  mAlloc(m_mucousDist, m_solver->a_noCells(), nDist, "mucousDist", -1.0, AT_);
  mAlloc(m_fluidDist, m_solver->a_noCells(), nDist, "fluidDist", -1.0, AT_);
  calculateSublayerDistances(index);
}

/** \brief Initialize bc20222 and bc20223. To apply the interpolated bounce back scheme with
 * a prescribed thermal flux as a Neumann condition we need to know the number of missing
 * distributions to calculate the factor to split the wall normal flux.
 *
 * \author Shota Ito
 * \date 24.07.2022
 **
 * \param[in] index Boundary index
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bcIBBNeumannInit(MInt index) {
  TRACE();
  // obtain a map where the number of missing distributions per cell are stored
  ibbBndIds.push_back(index);

  MIntScratchSpace noMissDist(m_solver->a_noCells(), AT_, "noMissDist");
  noMissDist.fill(0);

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt bndCellId = i;
    const MInt pCellId = m_bndCells[bndCellId].m_cellId;
    // case 1: boundary cell is fluid cell
    if(m_bndCells[bndCellId].m_isFluid) {
      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
        noMissDist(pCellId) = -1;
        continue;
      }
      for(MInt j = 0; j < nDist - 1; j++) {
        const MInt opposite = Ld::oppositeDist(j);
        if(m_bndCells[bndCellId].m_distances[j] <= 0.5 && m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          MInt neighborId = m_solver->c_neighborId(pCellId, opposite);
          MInt nbndid = m_solver->a_bndId(neighborId);
          if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) noMissDist(pCellId) = noMissDist(pCellId) + 1;
        }
      }
    }
    // case 2: boundary cell is solid cell
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        if(m_solver->a_hasNeighbor(pCellId, j) && m_bndCells[bndCellId].m_distances[j] < 0.5) {
          const MInt neighborId = m_solver->c_neighborId(pCellId, j);
          if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
            if(m_bndCells[bndCellId].m_distances[j] < 0.5) {
              MInt nbndid = m_solver->a_bndId(neighborId);
              if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
                if(m_solver->a_hasNeighbor(neighborId, j)) {
                  const MInt neighbor2Id = m_solver->c_neighborId(neighborId, j);
                  MInt n2bndid = m_solver->a_bndId(neighbor2Id);
                  if(n2bndid == -1 || m_bndCells[n2bndid].m_isFluid)
                    noMissDist(neighborId) = noMissDist(neighborId) + 1;
                }
              }
            }
          }
        }
      }
    }
  }

  for(MInt i = 0; i < m_solver->a_noCells(); i++) {
    MInt no = noMissDist(i);
    noMissDistBnd.push_back(no);
  }
}

/** \LBBC{secLBBC_bc20220, bc20220, 2220}
 * \author Moritz Waldmann
 * \date 10.12.2019
 *
 * Haenel wall bc with non zero velocity and thermal treatment
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20230(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20230(MInt index) {
  TRACE();

  MFloat F2q, rhoF, uF[nDim], uBF[nDim], tmpUF, tmpUBF, tmp2UF, b[2 * nDim], tmpUW, c[2 * nDim];
  MFloat uW[nDim] = {F0};
  getBoundaryVelocity(index, uW);

  MInt id, k, tpIndex;
  MFloat fEq, xi;
  MInt nghbrId;

  for(MInt d = 0; d < nDim; d++) {
    c[2 * d] = -uW[d];
    c[2 * d + 1] = uW[d];
  }

  const MFloat wT = getBoundaryTemperature(index);

  // go through m_bndcells with wallbndcnd
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    id = m_bndCells[i].m_cellId;

    m_omega = 2.0 / (1.0 + 6.0 * m_solver->a_nu(id) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(id)));

    rhoF = 1.0; // standard value which is usually overwritten in the following. Only in exceptional case it remains 1

    //----------------------------------------------------------------------------------------------------
    // perform bounce back for distributions piercing the wall - distinquish between inner and outer cells
    //----------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // 1. boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      if(m_densityFluctuations)
        rhoF = 1.0 + m_solver->a_variable(id, PV->RHO);
      else
        rhoF = m_solver->a_variable(id, PV->RHO);

      for(MInt n = 0; n < nDim; n++) {
        uF[n] = m_solver->a_variable(id, PV->VV[n]);
        b[2 * n] = -uF[n];
        b[2 * n + 1] = uF[n];
      }

      tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

      for(MInt j = 0; j < nDist - 1; j++) {
        if(j < 6) {
          tmpUF = b[j];
          tmpUW = c[j];
          tpIndex = 1;
        } else {
          if(j < 18) {
            k = j - Ld::distFld(0);
            tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
            tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          } else {
            k = j - (Ld::distFld(0) + Ld::distFld(1));
            tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
            tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
            tpIndex = 3;
          }
        }

        //---------------------------------------------------------------------------------
        // if q < 0.5, perform bounce back
        //---------------------------------------------------------------------------------

        if(m_bndCells[i].m_distances[j] < 0.5) {
          F2q = F2 * m_bndCells[i].m_distances[j];

          if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(j)) > 0) {
            nghbrId = m_solver->c_neighborId(id, Ld::oppositeDist(j));

            tmpUBF = 0.0;
            for(MInt n = 0; n < nDim; n++) {
              uBF[n] = m_solver->a_variable(nghbrId, PV->VV[n]);
              tmpUBF += (Ld::idFld(j, n) - 1) * uBF[n];
            }
          } else {
            tmpUBF = tmpUF;
          }

          fEq = Ld::tp(tpIndex) * rhoF
                * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
          xi = m_omega * (F2q - 1.0) / (1.0 - 2.0 * m_omega);
          // xi = m_omega * (F2q - 1.0) / (1.0 - m_omega);

          // perform interpolated bounce back
          m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) =
              (1.0 - xi) * m_solver->a_distribution(id, j) + xi * fEq
              - 2.0 * Ld::tp(tpIndex) * rhoF * F1BCSsq * tmpUW; // the minus corresponds to the opposite direction

        }

        // if the cell is located at a corner, there might be distributions missing
        else {
          // In this case simple bounce back is used. This could decrease the overall accuracy.
          if(m_solver->a_hasNeighbor(id, j) == 0) {
            m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) = m_solver->a_distribution(id, j);
          }
        }

      } // end of loop over all directions

    }

    //---------------------------------------------------------------------------------
    // 2. boundary cell is outside fluid
    // if the inner neighbor is a halo cell, it can be skipped
    //---------------------------------------------------------------------------------

    // The outer cell does not belong to the flow field, thus the macroscopic values are held constant.
    // It is possible that there are cells without any cutting velocity! Those have to be controlled, too.
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        //----------------------------------------------------------------------------------------------------------------
        // if q_innerNeighbor = 1 - q_bndCell >= 0.5, perform bounceback on inner neighbor
        //----------------------------------------------------------------------------------------------------------------
        if(m_solver->a_hasNeighbor(id, j) > 0 && !m_solver->a_isHalo(m_solver->c_neighborId(id, j))) {
          if(m_bndCells[i].m_distances[j] <= 0.5) {
            F2q = F2 - F2 * m_bndCells[i].m_distances[j]; // q = 1 - 0.5*(distance/cellHalfLength);

            nghbrId = m_solver->c_neighborId(id, j);

            if(m_densityFluctuations)
              rhoF = 1.0 + m_solver->a_variable(nghbrId, PV->RHO);
            else
              rhoF = m_solver->a_variable(nghbrId, PV->RHO);

            for(MInt n = 0; n < nDim; n++) {
              uF[n] = m_solver->a_variable(nghbrId, PV->VV[n]);
              b[2 * n] = -uF[n];
              b[2 * n + 1] = uF[n];
            }

            tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

            if(j < 6) {
              tmpUF = b[Ld::oppositeDist(j)];
              tmpUW = c[Ld::oppositeDist(j)];
              tpIndex = 1;
            } else {
              if(j < 18) {
                k = Ld::oppositeDist(j) - Ld::distFld(0);
                tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
                tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
                tpIndex = 2;
              } else {
                k = Ld::oppositeDist(j) - (Ld::distFld(0) + Ld::distFld(1));
                tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
                tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
                tpIndex = 3;
              }
            }

            tmpUBF = 0.0;
            for(MInt dim = 0; dim < nDim; dim++) {
              uBF[dim] = (1.0 - 3.0 / F2q) * uF[dim] + (3.0 / F2q) * uW[dim];
              tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
            }

            fEq = Ld::tp(tpIndex) * rhoF
                  * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
            xi = m_omega * (F2q - 1.0) / (1 + 0.5 * m_omega);
            // xi = m_omega * (F2q - 1.0);

            // perform interpolated bounce back
            m_solver->a_oldDistribution(nghbrId, j) =
                (1.0 - xi) * m_solver->a_distribution(nghbrId, Ld::oppositeDist(j)) + xi * fEq
                - 2.0 * Ld::tp(tpIndex) * rhoF * F1BCSsq * tmpUW; // the minus corresponds to the opposite direction
          }
        }

      } // end of loop over all directions

      calculateEqDistsWallSingleSpecies(id);
    }

    interpolatedBounceBackSingleSpeciesThermal(i, wT, uW);

  } // end of loop over all m_bndcells
}

/** \LBBC{secLBBC_bc30000, bc30000, 3000}
 * Lattice Boltzmann outflow boundary condition
 *
 * extrapolation of incoming distributions in arbitrary direction (rhs)
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc30000(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc30000(MInt index) {
  TRACE();

  MInt currentId, neighborId1, neighborId2;

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //      if( (globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    currentId = m_bndCells[i].m_cellId;

    // extrapolate inner distributions
    if(m_exDirs[index][0] >= 0 && m_solver->a_hasNeighbor(currentId, m_exDirs[ind][0]) > 0) {
      neighborId1 = m_solver->c_neighborId(currentId, m_exDirs[ind][0]);
      if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(currentId, m_exDirs[ind][1]) > 0) {
        neighborId2 = m_solver->c_neighborId(currentId, m_exDirs[ind][1]);
        for(MInt j = 0; j < nDist - 1; j++) {
          m_solver->a_oldDistribution(currentId, j) =
              m_exWeights[ind][0] * m_solver->a_oldDistribution(neighborId1, j)
              + m_exWeights[ind][1] * m_solver->a_oldDistribution(neighborId2, j);
        }
      } else {
        for(MInt j = 0; j < nDist - 1; j++) {
          m_solver->a_oldDistribution(currentId, j) = m_solver->a_oldDistribution(neighborId1, j);
        }
      }
    } else {
      if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(currentId, m_exDirs[ind][1]) > 0) {
        neighborId2 = m_solver->c_neighborId(currentId, m_exDirs[ind][1]);
        for(MInt j = 0; j < nDist - 1; j++) {
          m_solver->a_oldDistribution(currentId, j) = m_solver->a_oldDistribution(neighborId2, j);
        }
      }
    }
  }
}

/** \brief performes an interpolated bounce back as proposed by Bouzidi et al. 2001
 *
 * \author Andreas Lintermann
 * \date 04.10.2012
 *
 * In contrast to interpolatedBounceBackSingleSpecies, this additionally checks for neighboring cells that lie
 * on the other side of two walls and yields constitency if two neighbors want to modify the same cell.
 *
 * \param[in] cellId the ID of the current wall-cell
 *
 **/
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackSingleSpecies(const MInt cellId,
                                                                                   const MFloat* const uW) {
  TRACE();

  const MInt pCellId = m_bndCells[cellId].m_cellId;
  MFloat l_rho = m_solver->a_variable(pCellId, PV->RHO);

  // case 1: boundary cell is inside fluid
  if(m_bndCells[cellId].m_isFluid) {
    for(MInt j = 0; j < nDist - 1; j++) {
      // if q < 1.0, perform bounce back
      if(m_bndCells[cellId].m_distances[j] < 1.0) {
#ifdef WAR_NVHPC_PSTL
        const MInt opposite = m_solver->m_oppositeDist[j];
#else
        const MInt opposite = Ld::oppositeDist(j);
#endif
        if(m_bndCells[cellId].m_distances[opposite] < 1.0) {
          // if cut in opposite dir -> perform simple bounce back to avoid
          // multiple definition
          m_solver->a_oldDistribution(pCellId, opposite) =
              m_solver->a_distribution(pCellId, j) + firstMomentSourceTerm(uW, l_rho, opposite);
        } else {
          if(m_bndCells[cellId].m_distances[j] <= 0.5) {
            // TODO labels:LB it might be that no neighbor exist in opposite dir. For this
            // case a special treatment setting oldDist previously is demanded.
            // p.e.: intersecting of this and another BC ..
            // - ..inlet/outlet: if using calcluateEqDist(..) a_oldDistribution as
            //    well as a_distribution are set and hence no problem appears
            // - ..slipFlow: has to be performed before this BC
            const MFloat F2q = F2 * m_bndCells[cellId].m_distances[j];
            m_solver->a_oldDistribution(pCellId, opposite) = F2q * m_solver->a_distribution(pCellId, j)
                                                             + (F1 - F2q) * m_solver->a_oldDistribution(pCellId, j)
                                                             + firstMomentSourceTerm(uW, l_rho, opposite);
          } else {
            // if 0.5 < q < 1.0, perform bounce back, similar to bounce back from
            // a solid boundary cell's point of view
            const MFloat F2q = F2 * m_bndCells[cellId].m_distances[j];
            m_solver->a_oldDistribution(pCellId, opposite) =
                (m_solver->a_distribution(pCellId, j) / F2q)
                + (m_solver->a_distribution(pCellId, opposite) * (F2q - F1) / F2q)
                + firstMomentSourceTerm(uW, l_rho, opposite);
          }
        }
      }
    }
  }

  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
      // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
      if(m_bndCells[cellId].m_distances[j] < 0.5) {
        // do we have a non-halo fluid neighbor?
        if(m_solver->a_hasNeighbor(pCellId, j)) {
          const MInt neighborId = m_solver->c_neighborId(pCellId, j);
          if(!m_solver->a_isHalo(neighborId) && m_solver->a_isActive(neighborId)) {
            l_rho = m_solver->a_variable(neighborId, PV->RHO);
#ifdef WAR_NVHPC_PSTL
            const MInt opposite = m_solver->m_oppositeDist[j];
#else
            const MInt opposite = Ld::oppositeDist(j);
#endif
            const MFloat F2q = F2 * (F1 - m_bndCells[cellId].m_distances[j]);
            m_solver->a_oldDistribution(neighborId, j) = (m_solver->a_distribution(neighborId, opposite) / F2q)
                                                         + (m_solver->a_distribution(neighborId, j) * (F2q - F1) / F2q)
                                                         + (F1 / F2q) * firstMomentSourceTerm(uW, l_rho, j);
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDists(pCellId, 1.0, uW);
    for(MInt d = 0; d < nDim; d++) {
      m_solver->a_variable(pCellId, d) = uW[d];
      m_solver->a_oldVariable(pCellId, d) = uW[d];
    }
    m_solver->a_variable(pCellId, nDim) = F1;
    m_solver->a_oldVariable(pCellId, nDim) = F1;
  }
}


/** \brief performes an interpolated bounce back for thermal as proposed by Like Li et al. 2012
 *
 * \author Moritz Waldmann, Jie Ruan
 * \date 24.08.2018
 *
 * This function implement bounce back idea to thermal boundary condition.
 *
 * \param[in]  cellId Cell id of current boundary cell
 * \param[in]  wT     Wall temperature at boundary cell
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackSingleSpeciesThermal(const MInt cellId,
                                                                                          const MFloat* const wTPtr,
                                                                                          const MFloat* const uW) {
  TRACE();

  const MInt pCellId = m_bndCells[cellId].m_cellId;

  std::array<MBool, nDist - 1> ibb{};
  ibb.fill(false);

  // case 1: boundary cell is inside fluid
  if(m_bndCells[cellId].m_isFluid) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MInt opposite = Ld::oppositeDist(j);
#endif
      MFloat s = F0;
      if(m_solver->m_innerEnergy) {
        s = zerothMomentSourceTerm<1>(pCellId, opposite, wTPtr[j], uW);
      } else if(m_solver->m_totalEnergy) {
        s = zerothMomentSourceTerm<2>(pCellId, opposite, wTPtr[j], uW);
      } else {
        s = zerothMomentSourceTerm<0>(pCellId, opposite, wTPtr[j], uW);
      }
      const MFloat sourceTerm = s;
      // if q <= 0.5, perform bounce back
      if(m_bndCells[cellId].m_distances[j] <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          MFloat F2q = F2 * m_bndCells[cellId].m_distances[j];
          m_solver->a_oldDistributionThermal(pCellId, opposite) =
              (-F2q * m_solver->a_distributionThermal(pCellId, j))
              + ((F2q - F1) * m_solver->a_oldDistributionThermal(pCellId, j)) + sourceTerm;
          ibb[j] = true;
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else
          m_solver->a_oldDistributionThermal(pCellId, opposite) = m_solver->a_distributionThermal(pCellId, j);
      }

      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      else if(m_bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
        m_solver->a_oldDistributionThermal(pCellId, opposite) = m_solver->a_distributionThermal(pCellId, j);
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MInt opposite = Ld::oppositeDist(j);
#endif
      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        MInt neighborId = m_solver->c_neighborId(pCellId, j);

        MFloat s = F0;
        if(m_solver->m_innerEnergy) {
          s = zerothMomentSourceTerm<1>(neighborId, j, wTPtr[j], uW);
        } else if(m_solver->m_totalEnergy) {
          s = zerothMomentSourceTerm<2>(neighborId, j, wTPtr[j], uW);
        } else {
          s = zerothMomentSourceTerm<0>(neighborId, j, wTPtr[j], uW);
        }
        const MFloat sourceTerm = s;

        if(m_bndCells[cellId].m_distances[j] < 0.5) {
          ibb[j] = true; // must be set here otherwise the wT is computed wrongly
        }
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          MInt nbndid = m_solver->a_bndId(neighborId);

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(m_bndCells[cellId].m_distances[j] < 0.5) {
            MFloat F2q = F2 * (F1 - m_bndCells[cellId].m_distances[j]);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
              m_solver->a_oldDistributionThermal(neighborId, j) =
                  -(m_solver->a_distributionThermal(neighborId, opposite) / F2q)
                  + (m_solver->a_distributionThermal(neighborId, j) * (F2q - F1) / F2q) + sourceTerm / F2q;
            }
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    MFloat wTMean = F0;
    MInt counterT = 0;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(ibb[j]) {
        wTMean += wTPtr[j];
        counterT++;
      }
    }
    if(counterT != 0) {
      wTMean /= counterT;
    } else {
      wTMean = 1.0; // this case should not happen;
    }
    m_solver->setEqDistsThermal(pCellId, wTMean, F1, uW);
    m_solver->a_variable(pCellId, PV->T) = wTMean;
    m_solver->a_oldVariable(pCellId, PV->T) = wTMean;
  }
}

/** \brief performes an interpolated bounce back for transport phenomena as proposed by Like Li et al. 2012
 *
 * \author Shota Ito
 * \date 08.06.2022
 *
 * This function implement bounce back idea to transport boundary condition.
 *
 * \param[in]  cellId Cell id of current boundary cell
 * \param[in]  wC     Wall concentration at boundary cell
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackSingleSpeciesTransport(const MInt cellId,
                                                                                            const MFloat* const wCPtr,
                                                                                            const MFloat* const uW) {
  TRACE();

  const MInt pCellId = m_bndCells[cellId].m_cellId;

  std::array<MBool, nDist - 1> ibb{};
  ibb.fill(false);

  // case 1: boundary cell is inside fluid
  if(m_bndCells[cellId].m_isFluid) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MInt opposite = Ld::oppositeDist(j);
#endif
      const MFloat sourceTerm = zerothMomentSourceTerm<0>(pCellId, opposite, wCPtr[j], uW);
      // if q <= 0.5, perform bounce back
      if(m_bndCells[cellId].m_distances[j] <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          MFloat F2q = F2 * m_bndCells[cellId].m_distances[j];
          m_solver->a_oldDistributionTransport(pCellId, opposite) =
              (-F2q * m_solver->a_distributionTransport(pCellId, j))
              + ((F2q - F1) * m_solver->a_oldDistributionTransport(pCellId, j)) + sourceTerm;
          ibb[j] = true;
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else
          m_solver->a_oldDistributionTransport(pCellId, opposite) = m_solver->a_distributionTransport(pCellId, j);
      }

      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      else if(m_bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
        m_solver->a_oldDistributionTransport(pCellId, opposite) = m_solver->a_distributionTransport(pCellId, j);
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MInt opposite = Ld::oppositeDist(j);
#endif
      const MFloat sourceTerm = zerothMomentSourceTerm<0>(cellId, j, wCPtr[j], uW);
      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        MInt neighborId = m_solver->c_neighborId(pCellId, j);
        if(m_bndCells[cellId].m_distances[j] < 0.5) {
          ibb[j] = true; // this must be set here otherwise wCMean is computed wrong
        }
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          MInt nbndid = m_solver->a_bndId(neighborId);

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(m_bndCells[cellId].m_distances[j] < 0.5) {
            MFloat F2q = F2 * (F1 - m_bndCells[cellId].m_distances[j]);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
              m_solver->a_oldDistributionTransport(neighborId, j) =
                  -(m_solver->a_distributionTransport(neighborId, opposite) / F2q)
                  + (m_solver->a_distributionTransport(neighborId, j) * (F2q - F1) / F2q) + sourceTerm / F2q;
            }
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    MFloat wCMean = F0;
    MInt counterC = 0;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(ibb[j]) {
        wCMean += wCPtr[j];
        counterC++;
      }
    }
    if(counterC != 0) {
      wCMean /= counterC;
    } else {
      wCMean = 1.0; // this case should not happen;
    }
    m_solver->setEqDistsTransport(pCellId, wCMean, uW);
    m_solver->a_variable(pCellId, nDim + 2) = wCMean;
    m_solver->a_oldVariable(pCellId, nDim + 2) = wCMean;
  }
}

/** \brief Calculates the source term for the zeroth moment with is required for thermal bounce back schemes
 * (Li et al. 2012)
 *
 * \author Moritz Waldmann
 * \date 17.11.2023
 *
 * This function calculates the source term for the zeroth moment for bounce back schemes. The template
 * mode switches between the different types of equilibrium functions, e.g., mode = 0 corresponds to the standard
 * Eq. function, mode = 1 corresponds to the inner energy Eq., and mode = 2 corresponds to the total energy Eq.
 *
 * \tparam     mode     EQ distribution type: 0=standard, 1=inner energy, 2=total energy
 * \param[in]  pCellId 	Cell id of current boundary cell
 * \param[in]  dist     Current direction of the PPDF
 * \param[in]  var      Value of the variable used for the source term
 * \param[in]  uW       wall velocity
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt mode>
MFloat LbBndCndDxQy<nDim, nDist, SysEqn>::zerothMomentSourceTerm(const MInt pCellId, const MInt dist, const MFloat var,
                                                                 const MFloat* uW) {
  TRACE();

  MFloat sourceTerm = F0;
  const MFloat rho = m_solver->a_variable(pCellId, PV->RHO);
  const MFloat D = (MFloat)nDim;
#ifdef WAR_NVHPC_PSTL
  const MFloat distType = (MFloat)m_solver->m_distType[dist];
  const MFloat weight = m_solver->m_tp[m_solver->m_distType[dist]];
  MFloat sP = F0;
  MFloat sV = F0;
  for(MInt d = 0; d < nDim; d++) {
    sP += m_solver->m_ppdfDir[dist * nDim + d] * uW[d];
    sV += uW[d] * uW[d];
  }
  const MFloat tmp = sP * sP * F1B2mulF1BCSsq * F1BCSsq;
  IF_CONSTEXPR(mode == 0) { sourceTerm = F2 * weight * var * (F1 + tmp - sV * F1B2mulF1BCSsq); }
  IF_CONSTEXPR(mode == 1) { sourceTerm = weight * var * D * rho * (distType + tmp - sV * F1B2mulF1BCSsq); }
  IF_CONSTEXPR(mode == 2) {
    const MFloat E = D * var * F1B2;
    sourceTerm = F2 * weight * rho
                 * (tmp - sV + F1B2 * (distType - D * CSsq) + E * (F1 + tmp * F1B2mulF1BCSsq - sV * F1B2mulF1BCSsq));
  }
#else
  const MFloat weight = Ld::tp(Ld::distType(dist));
  const MFloat distType = (MFloat)Ld::distType(dist);
  const MFloat sP = std::inner_product(uW, uW + nDim, lbDescriptor::ppdfDir<nDim>[dist], .0);
  const MFloat sV = std::inner_product(uW, uW + nDim, uW, .0);
  const MFloat tmp = sP * sP * F1B2mulF1BCSsq * F1BCSsq;
  IF_CONSTEXPR(mode == 0) { sourceTerm = F2 * weight * var * (F1 + tmp - sV * F1B2mulF1BCSsq); }
  IF_CONSTEXPR(mode == 1) { sourceTerm = weight * var * D * rho * (distType + tmp - sV * F1B2mulF1BCSsq); }
  IF_CONSTEXPR(mode == 2) {
    const MFloat E = D * var * F1B2;
    sourceTerm = F2 * weight * rho
                 * (tmp - sV + F1B2 * (distType - D * CSsq) + E * (F1 + tmp * F1B2mulF1BCSsq - sV * F1B2mulF1BCSsq));
  }
#endif
  return sourceTerm;
}

/** \brief performes an interpolated bounce back for thermal LBGK with a heat flux to the wall normal proposed by Like
 * Li et al. 2012
 *
 * \author Shota Ito
 * \date 02.07.2022
 *
 * This function implement bounce back idea to thermal boundary condition with a heat flux.
 *
 * \param[in]  cellId 	Boundary cell id of current boundary cell
 * \param[in]  qT     	Magnitude of the heat flux normal to the boundary wall
 * \param[in]  bcIndex  Index of the boundary condition
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackSingleSpeciesThermalFlux(MInt cellId, MFloat qT,
                                                                                              MInt bcIndex) {
  TRACE();

  const MInt pCellId = m_bndCells[cellId].m_cellId;
  const MFloat l_dt = FPOW2(m_solver->maxLevel() - m_solver->a_level(pCellId)); // local time step
  const MFloat LB_dx = 1.0;                                                     // local grid space in LB units
  const MInt noCells = m_solver->a_noCells();

  MInt noMissDist = 0;
  MFloat l_q = qT; // Heat flux normal to wall

  // case 1: boundary cell is inside fluid
  if(m_bndCells[cellId].m_isFluid) {
    // Calculate the factor to split the flux
    MInt bndId = 0;
    for(MInt n = 0; n < (MInt)ibbBndIds.size(); n++) {
      if(ibbBndIds.at(n) == bcIndex) {
        bndId = n;
        break;
      }
    }
    noMissDist = noMissDistBnd.at(bndId * noCells + pCellId);
    if(m_solver->a_hasProperty(pCellId, Cell::IsHalo))
      return;
    else if(noMissDist == 0) {
      std::cout << "ERROR this should not happen!" << std::endl;
    }
    const MFloat splittingFactor = 1.0 / (1.0 * noMissDist); // in a fluid cell this is constant!
    MFloat fluxInDistDir = splittingFactor * l_q;            // split flux in the distribution directions

    // now loop over the distributions
    for(MInt j = 0; j < nDist - 1; j++) {
      const MInt opposite = Ld::oppositeDist(j);

      if(m_bndCells[cellId].m_distances[j] <= 0.5) {
        // Now perform bounce back
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          MInt neighborId = m_solver->c_neighborId(pCellId, opposite);
          MInt nbndid = m_solver->a_bndId(neighborId);

          // check if my neighbor is an inside cell or a boundary cell with the cell center inside the fluid
          if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
            MFloat FFq = F2 * m_bndCells[cellId].m_distances[j] + F1;
            MFloat Fq = F2 * m_bndCells[cellId].m_distances[j] - F1;
            m_solver->a_oldDistributionThermal(pCellId, opposite) =
                m_solver->a_distributionThermal(pCellId, j)
                - (Fq / FFq) * m_solver->a_distributionThermal(neighborId, j)
                + (Fq / FFq) * m_solver->a_distributionThermal(pCellId, opposite)
                + (F2 / FFq) * (l_dt / LB_dx) * fluxInDistDir;
          }
          // this is a strange case where we don't have valid neighbors; do simple bounce back
          else {
            m_solver->a_oldDistributionThermal(pCellId, opposite) = m_solver->a_distributionThermal(pCellId, j);
          }
        }
        // this is a strange case where we don't have any neighbors in the direction where we want to set the
        // distribution; do simple bounce back
        else {
          m_solver->a_oldDistributionThermal(pCellId, opposite) = m_solver->a_distributionThermal(pCellId, j);
        }
      }
      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      else if(m_bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistributionThermal(pCellId, opposite) = m_solver->a_distributionThermal(pCellId, j);
      }
    }
  }

  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
      const MInt opposite = Ld::oppositeDist(j);

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        MInt neighborId = m_solver->c_neighborId(pCellId, j);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          MInt nbndid = m_solver->a_bndId(neighborId);
          // if q < 0.5, perform bounce back for the neighbor fluid cell, where the fluid cell has q > 0.5
          if(m_bndCells[cellId].m_distances[j] < 0.5) {
            // Calculate q for the neighbor cell
            MFloat nei_q = F1 - m_bndCells[cellId].m_distances[j];
            // check if my neighbor is an inside cell or a boundary cell with cell center inside the fluid
            if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
              // check if there is a neighbor of the neighbor cell since we need that cell too for the interpolation
              if(m_solver->a_hasNeighbor(neighborId, j)) {
                MInt neighbor2Id = m_solver->c_neighborId(neighborId, j);
                MInt n2bndid = m_solver->a_bndId(neighbor2Id);
                // check if the neighbors neighbor cell is a valid cell
                if(n2bndid == -1 || m_bndCells[n2bndid].m_isFluid) {
                  // Calculate factor to split the flux
                  // Here the splitting factors could differ for each distributions, thus the calculation is performed
                  // inside the loop
                  MInt bndId = 0;
                  for(MInt n = 0; n < (MInt)ibbBndIds.size(); n++) {
                    if(ibbBndIds.at(n) == bcIndex) {
                      bndId = n;
                      break;
                    }
                  }
                  noMissDist = noMissDistBnd.at(bndId * noCells + neighborId);
                  if(noMissDist == 0) {
                    std::cout << "ERROR this should not happen!" << std::endl;
                  }
                  MFloat splittingFactor = 1.0 / (1.0 * noMissDist);
                  MFloat fluxInDistDir = splittingFactor * l_q;

                  // perform interpolated bounce back
                  MFloat FFq = F2 * nei_q + F1;
                  MFloat Fq = F2 * nei_q - F1;
                  m_solver->a_oldDistributionThermal(neighborId, j) =
                      m_solver->a_distributionThermal(neighborId, opposite)
                      - (Fq / FFq) * m_solver->a_distributionThermal(neighbor2Id, opposite)
                      + (Fq / FFq) * m_solver->a_distributionThermal(neighborId, j)
                      + (F2 / FFq) * (l_dt / LB_dx) * fluxInDistDir;
                }
                // if there is no second neighbor, then do simple bounce back
                else {
                  m_solver->a_oldDistributionThermal(neighborId, j) =
                      m_solver->a_distributionThermal(neighborId, opposite);
                }
              }
              // if there is no second neighbor, then do simple bounce back
              else {
                m_solver->a_oldDistributionThermal(neighborId, j) =
                    m_solver->a_distributionThermal(neighborId, opposite);
              }
            }
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incomming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! There are considered too.
    std::vector<MFloat> wallTemp;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_bndCells[cellId].m_distances[j] < 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, j)) {
          MFloat temp = m_solver->a_variable(m_solver->c_neighborId(pCellId, j), PV->T);
          wallTemp.push_back(temp);
        }
      }
    }
    const MInt noCuts = (MInt)wallTemp.size();
    MFloat sumTemp = F0;
    for(MInt n = 0; n < (signed)wallTemp.size(); n++) {
      sumTemp += wallTemp.at(n);
    }
    MFloat wT = sumTemp / ((MInt)noCuts);
    calculateEqDistsWallSingleSpeciesThermal(pCellId, wT);
  }
}

/** \brief Calculates the distances through a specified sublayer in each distribution direction.
 *
 * \author Shota Ito
 * \date 27.09.2022
 *
 * \param[in]  cellId Cell id of current boundary cell
 * \param[in]  index  index of the current BC
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateSublayerDistances(MInt index) {
  TRACE();

  for(MInt cellId = m_bndCndOffsets[index]; cellId < m_bndCndOffsets[index + 1]; cellId++) {
    const MInt pCellId = m_bndCells[cellId].m_cellId;
    // const MFloat l_dx = m_solver->c_cellLengthAtLevel(m_solver->a_level(pCellId));
    for(MInt j = 0; j < nDist - 1; j++) {
      const MInt opposite = Ld::oppositeDist(j);
      // skip the cases where no bounce back is performed
      if(m_bndCells[cellId].m_isFluid) {
        if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;
        if(m_bndCells[cellId].m_distances[j] > 0.5) {
          continue;
        } else {
          if(m_solver->a_hasNeighbor(pCellId, opposite) == 0) continue;
        }
      } else {
        if(!m_solver->a_hasNeighbor(pCellId, j)) {
          continue;
        } else {
          if(m_bndCells[cellId].m_distances[j] > 0.5) {
            continue;
          } else {
            MInt nbndid = m_solver->a_bndId(m_solver->c_neighborId(pCellId, j));
            if(!(nbndid == -1 || m_bndCells[nbndid].m_isFluid)) continue;
          }
        }
      }

      // in order to compute the wall temperature we first need the distance in the sublayer from the wall cut point
      // along each distribution. For this we need the normal vector in each wall cut point which points into the fluid
      // domain
      MFloat wallNormal[nDim] = {F0};

      MFloat lDir = F0;
      for(MInt d = 0; d < nDim; d++)
        lDir += Ld::ppdfDir(j, d) * Ld::ppdfDir(j, d);
      lDir = sqrt(lDir); // this is to correct the different length for the diagonal and normal directions
      MFloat dist2Wall;  // distance from the cell center to the wall in cartesian coordinates
      if(m_bndCells[cellId].m_isFluid) {
        dist2Wall = m_bndCells[cellId].m_distances[j] * lDir;
      } else {
        dist2Wall = (1.0 - m_bndCells[cellId].m_distances[j]) * lDir;
      }

      if(m_calcSublayerDist) {
        MFloat** const v = m_solver->m_geometry->elements[m_distIntersectionElementId[pCellId][j]].m_vertices;
        MFloat edge[nDim];
        MFloat edge2[nDim];
        if(nDim == 2) {
          for(MInt dim = 0; dim < nDim; dim++)
            edge[dim] = v[0][dim] - v[1][dim];
          wallNormal[0] = edge[1];
          wallNormal[1] = -edge[0];
        } else if(nDim == 3) {
          for(MInt dim = 0; dim < nDim; dim++) {
            edge[dim] = v[0][dim] - v[1][dim];
            edge2[dim] = v[0][dim] - v[2][dim];
          }
          wallNormal[0] = edge[1] * edge2[2] - edge[2] * edge2[1];
          wallNormal[1] = edge[2] * edge2[0] - edge[0] * edge2[2];
          wallNormal[2] = edge[0] * edge2[1] - edge[1] * edge2[0];
        }
      }

      // direction where the cut was found
      MFloat bounceBackDir[nDim];
      for(MInt d = 0; d < nDim; d++) {
        bounceBackDir[d] = Ld::ppdfDir(j, d);
      }

      // dot product and calculate length of normal vector and direction vector
      MFloat dp = F0;
      MFloat lv = F0;
      MFloat lb = F0;
      for(MInt d = 0; d < nDim; d++) {
        dp += bounceBackDir[d] * wallNormal[d];
        lv += wallNormal[d] * wallNormal[d];
        lb += bounceBackDir[d] * bounceBackDir[d];
      }
      lv = sqrt(lv);
      lb = sqrt(lb);

      // invert normal vector if needed
      if(dp < 0) {
        for(MInt d = 0; d < nDim; d++) {
          wallNormal[d] = -wallNormal[d];
        }
        dp = -dp;
      }

      // Check if the normal vector is valid
      if(approx(lv, 0.0, MFloatEps)) continue;

      // Now calculate the angle between the normal vector and the distribution
      MFloat phi = acos(dp / (lv * lb));

      // The distance in the sublayer is given as
      m_mucousDist[pCellId][j] = m_mucosaModel.mucosaThickness / cos(phi);
      m_fluidDist[pCellId][j] = dist2Wall;
    }
  }
}

/** \brief Calculates the interface concentration and temperature for the interpolated bounce back
 * in each distribution direction.
 *
 * \author Shota Ito
 * \date 27.09.2022
 *
 * \param[in]  cellId Cell id of current boundary cell
 * \param[in]  index  index of the current BC
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateWallInterface(MInt cellId, MFloat* wallConc, MFloat* wallTemp) {
  TRACE();
  const MInt pCellId = m_bndCells[cellId].m_cellId;
  MInt fluidCellId = pCellId;

  for(MInt j = 0; j < nDist - 1; j++) {
    const MInt opposite = Ld::oppositeDist(j);
    // skip the cases where no bounce back is performed
    if(m_bndCells[cellId].m_isFluid) {
      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;
      if(m_bndCells[cellId].m_distances[j] > 0.5) {
        continue;
      } else {
        if(!m_solver->a_hasNeighbor(pCellId, opposite)) continue;
      }
    } else {
      if(!m_solver->a_hasNeighbor(pCellId, j)) {
        continue;
      } else {
        if(m_bndCells[cellId].m_distances[j] > 0.5) {
          continue;
        } else {
          MInt nbndid = m_solver->a_bndId(m_solver->c_neighborId(pCellId, j));
          if(!(nbndid == -1 || m_bndCells[nbndid].m_isFluid)) {
            continue;
          } else {
            fluidCellId = m_solver->c_neighborId(pCellId, j);
          }
        }
      }
    }

    // ----------------------------------------------
    // 1. Calculate Cw
    // ----------------------------------------------

    // Get the necessary dimensioned variables
    MFloat d_Cf = m_solver->a_variable(fluidCellId, PV->C) * m_mucosaModel.refC;
    MFloat d_C_o = m_mucosaModel.C_o * m_mucosaModel.refC;
    MFloat fD_mD = m_fluidDist[pCellId][j] / m_mucousDist[pCellId][j];

    MFloat d_Cw = (d_Cf + m_mucosaModel.diffRatio * fD_mD * d_C_o) / (1.0 + m_mucosaModel.diffRatio * fD_mD);
    wallConc[j] = d_Cw / m_mucosaModel.refC;

    // ----------------------------------------------
    // 2. Calculate latent heat (if flag is on)
    // ----------------------------------------------

    MFloat d_latentSourceTerm = F0;
    if(m_latentHeat) {
      // First, get the old wall temperature
      MFloat Tw_old = m_oldWallTemp[pCellId][j];

      // Convert the dimensionless variables
      MFloat d_Tw_old = Tw_old * m_mucosaModel.refT;
      MFloat d_dist2Wall = m_fluidDist[pCellId][j] * m_mucosaModel.refDx;
      MFloat d_diff = m_solver->m_diffusivity * m_mucosaModel.refDiff;
      MFloat d_cond_f = m_mucosaModel.refCondF;

      // Calculate the specific latent heat
      // Inthavong et al (2018) "Wet surface wall model for latent heat exchange during evaporation"
      MFloat d_h_lh = (2500.8 - 2.3641 * d_Tw_old + 1.5893 * 1e-3 * d_Tw_old * d_Tw_old
                       - 6.1434 * 1e-6 * d_Tw_old * d_Tw_old * d_Tw_old)
                      * 1000;

      // Compute the water flux across the interface
      MFloat d_j_flux = (d_Cw - d_Cf) / d_dist2Wall * d_diff;

      // Compute the latent heat along this distribution
      MFloat d_q_latent = -d_j_flux * d_h_lh;

      d_latentSourceTerm = (d_q_latent * d_dist2Wall / d_cond_f) / (1.0 + m_mucosaModel.condRatio * fD_mD);
    }

    // ----------------------------------------------
    // 3. Calculate Tw, Cw
    // ----------------------------------------------

    // Finally, compute the wall temperature
    // Again first dimensinon the necessary variables
    MFloat d_Tf = m_solver->a_variable(fluidCellId, PV->T) * m_mucosaModel.refT;
    MFloat d_T_o = m_mucosaModel.T_o * m_mucosaModel.refT;

    // Calculate the wall temperature in celsius
    MFloat d_Tw =
        (d_Tf + m_mucosaModel.condRatio * fD_mD * d_T_o) / (1.0 + m_mucosaModel.condRatio * fD_mD) + d_latentSourceTerm;

    // Calculate the dimensionless wall temperature
    wallTemp[j] = d_Tw / m_mucosaModel.refT;

    // save the old wall temperature
    if(m_latentHeat) {
      m_oldWallTemp[pCellId][j] = wallTemp[j];
    }
  }
}

/** \brief Applies interpolated anti bounce back to given moving boundary cell
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date 01.2020
 *
 *  For reference, see Dubois et al., On anti bounce back boundary condition for lattice Boltzmann schemes
 *  https://arxiv.org/abs/1812.04305
 *
 *  \param[in] cellIndex Moving boundary cell the BC is applied to
 *  \param[in] set The boundary condition is applied to all cells which belong to this level set.
 **/
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedAntiBounceBackMb_Bouzidi_qua(const MInt cellIndex,
                                                                                        const MInt set) {
  TRACE();

  const MInt pCellId = m_boundaryCellsMb.cellId(cellIndex);

  const MFloat rho = m_boundaryCellsMb.density(cellIndex);

  std::array<MFloat, nDim> uBoundary{};
  getBoundaryVelocityMb(cellIndex, uBoundary.data());
  const MFloat squaredVelocity = std::inner_product(uBoundary.begin(), uBoundary.end(), uBoundary.begin(), .0);

  std::array<MFloat, nDist> eqDist{};
  sysEqn().calcEqDists(rho, squaredVelocity, uBoundary.data(), eqDist.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
        continue;
      }

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite)) {
          const MInt neighborId = m_solver->c_neighborId(pCellId, opposite);

          const MFloat F2q = F2 * q;
          const MFloat Fq = q;

          const MFloat quadratic = -Fq * (F2q + 1) * m_solver->a_distribution(pCellId, j)
                                   - (1 + F2q) * (1 - F2q) * m_solver->a_oldDistribution(pCellId, j)
                                   + Fq * (1 - F2q) * m_solver->a_oldDistribution(neighborId, j) + eqDist[j]
                                   + eqDist[opposite];

          m_solver->a_oldDistribution(pCellId, opposite) = quadratic;
        }
        // this is the case in which we have no neighbor in the direction we want to set (strange case)
        // can in my opinion only appear if neighbors were set wrong or at an interface
        else {
          m_solver->a_oldDistribution(pCellId, opposite) =
              -m_solver->a_distribution(pCellId, j) + eqDist[j] + eqDist[opposite];
        }
      }
      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      // else if(m_bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      // else if( m_wallBoundaryCellList[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistribution(pCellId, opposite) =
            -m_solver->a_distribution(pCellId, j) + eqDist[j] + eqDist[opposite];
      }
    }
  }
  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);

      const MInt opposite = Ld::oppositeDist(j);

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);

        // if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)){

        const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

        // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
        if(q < 0.5) {
          const MFloat F2q = F2 * (F1 - q);
          const MFloat Fq = (1 - q);

          // check if my neighbor is an inside cell or a boundary cell with the cell center inside
          if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
            const MInt nneighborId = m_solver->c_neighborId(neighborId, j);

            const MFloat quadratic =
                1 / Fq / (F2q + 1) * (-m_solver->a_distribution(neighborId, opposite) + eqDist[j] + eqDist[opposite])
                + (F2q - 1) / Fq * m_solver->a_distribution(neighborId, j)
                - (F2q - 1) / (F2q + 1) * m_solver->a_distribution(nneighborId, j);

            m_solver->a_oldDistribution(neighborId, j) = quadratic;
          }
        }
        //}
      }
    } // end of the loop over all PPDF directions in case 2

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDists(pCellId, 1.0, squaredVelocity, uBoundary.data());
  }
}


/** \brief Calculates the equilibrium PPDFs for walls
 *
 * \author Andreas Lintermann
 * \date 04.10.2012
 *
 * \param[in] pCellId Current cell id
 */
template <MInt nDim, MInt nDist, class SysEqn>
ATTRIBUTES1(ATTRIBUTE_NO_AUTOVEC)
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateEqDistsWallSingleSpecies(const MInt pCellId) {
  const MFloat rho = 1.0;
  const std::array<MFloat, nDim> vel{0.0};

  // Set macroscopic variables
  m_solver->a_variable(pCellId, PV->RHO) = rho;
  m_solver->a_oldVariable(pCellId, PV->RHO) = rho;
  for(MInt n = 0; n < nDim; n++) {
    m_solver->a_variable(pCellId, n) = vel[n];
    m_solver->a_oldVariable(pCellId, n) = vel[n];
  }

  m_solver->setEqDists(pCellId, rho, vel.data());
}

/** \brief Calculates the equilibrium PPDFs for walls for thermal LBGK
 *
 * \author Andreas Lintermann
 * \date 04.10.2012
 *
 * \param[in] pCellId Cell id of current boundary cell
 * \param[in] wT      Wall temperature
 **/
template <MInt nDim, MInt nDist, class SysEqn>
ATTRIBUTES1(ATTRIBUTE_NO_AUTOVEC)
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateEqDistsWallSingleSpeciesThermal(const MInt pCellId,
                                                                                        const MFloat wT) {
  TRACE();

  m_solver->a_variable(pCellId, PV->T) = wT;
  m_solver->a_oldVariable(pCellId, PV->T) = wT;

  constexpr std::array<MFloat, nDim> u{};
  m_solver->setEqDistsThermal(pCellId, wT, F1, F0, u.data());
}

/** \brief Calculates the equilibrium PPDFs for walls for transport LBGK
 *
 * \author Shota Ito
 * \date 08.06.2022
 *
 * \param[in] pCellId Cell id of current boundary cell
 * \param[in] wC      Wall concentration
 **/
template <MInt nDim, MInt nDist, class SysEqn>
ATTRIBUTES1(ATTRIBUTE_NO_AUTOVEC)
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateEqDistsWallSingleSpeciesTransport(const MInt pCellId,
                                                                                          const MFloat wC) {
  TRACE();

  m_solver->a_variable(pCellId, PV->C) = wC;
  m_solver->a_oldVariable(pCellId, PV->C) = wC;

  constexpr std::array<MFloat, nDim> u{};
  m_solver->setEqDistsTransport(pCellId, wC, 0.0, u.data());
}

/** \brief Extrapolates inner variables to a cell
 *
 * \author Andreas Lintermann
 * \date 02.10.2012
 *
 * Note that if no neighbors in the extrapolation directions exist, then no value is assigned to p_var.
 *
 * \param[in]  index   Index of the BC
 * \param[in]  pCellId Cell id
 * \param[out] p_var   Pointer to the density
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::extrapolateVariable(const MInt index, const MInt pCellId, const MInt var,
                                                                   MFloat* const p_var) {
  TRACE();

  if(m_exDirs[index][0] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][0]) > 0) {
    MInt neighborId1 = m_solver->c_neighborId(pCellId, m_exDirs[index][0]);

    if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][1]) > 0) {
      MInt neighborId2 = m_solver->c_neighborId(pCellId, m_exDirs[index][1]);

      MFloat ex0 = m_exWeights[index][0];
      MFloat ex1 = m_exWeights[index][1];

      *p_var = ex0 * m_solver->a_variable(neighborId1, var) + ex1 * m_solver->a_variable(neighborId2, var);
    } else
      *p_var = m_solver->a_variable(neighborId1, var);
  } else {
    if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][1]) > 0) {
      MInt neighborId2 = m_solver->c_neighborId(pCellId, m_exDirs[index][1]);
      *p_var = m_solver->a_variable(neighborId2, var);
    }
  }
}

/** \brief Extrapolates inner velocities to a cell
 *
 * \author Andreas Lintermann
 * \date 28.09.2012
 *
 * Note that if no neighbors in the extrapolation directions exist, then no value is assigned to the velocities.
 *
 * \param[in]  index   The index of the BC
 * \param[in]  pCellId Cell id
 * \param[out] l_uu    Pointer to the velocity
 **/
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::extrapolateVelocities(MInt index, const MInt pCellId, MFloat* l_uu) {
  TRACE();

  if(m_exDirs[index][0] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][0]) > 0) {
    const MInt neighborId1 = m_solver->c_neighborId(pCellId, m_exDirs[index][0]);

    if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][1]) > 0) {
      const MInt neighborId2 = m_solver->c_neighborId(pCellId, m_exDirs[index][1]);

      const MFloat ex0 = m_exWeights[index][0];
      const MFloat ex1 = m_exWeights[index][1];

      for(MInt n = 0; n < nDim; n++) {
        l_uu[n] = ex0 * m_solver->a_variable(neighborId1, n) + ex1 * m_solver->a_variable(neighborId2, n);
      }

    } else {
      for(MInt n = 0; n < nDim; n++) {
        l_uu[n] = m_solver->a_variable(neighborId1, n);
      }
    }
  } else if(m_exDirs[index][1] >= 0 && m_solver->a_hasNeighbor(pCellId, m_exDirs[index][1]) > 0) {
    const MInt neighborId2 = m_solver->c_neighborId(pCellId, m_exDirs[index][1]);

    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_solver->a_variable(neighborId2, n);
    }
  } else {
#ifndef WAR_NVHPC_PSTL
    TERMM(1, "No value is extrapolated !");
#endif
  }
}

/** \LBBC{sec_LBBC_bc10000, bc10000, 1000}
 * \author Georg Eitel-Amor, Andreas Lintermann
 * \date 02.10.2012
 *
 * Lattice Boltzmann inflow boundary condition
 *
 * Diriclet condition with prescribed eq-distributions.
 * Velocity vector is read from properties and prescribed on boundary.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10000(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10000(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat ma_x_F1BCS = m_Ma * LBCS;
  const MInt pvrho = PV->RHO;

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    const MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    MFloat l_rho = 1.0;
    MFloat l_uu[nDim] = {0.0};

    // extrapolate density
    extrapolateVariable(ind, pCellId, pvrho, &l_rho);

    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * ma_x_F1BCS * (*bndCell).m_multiplier * m_zeroInflowVelocity;
    }

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, l_uu);

    m_solver->a_variable(pCellId, pvrho) = l_rho;
    for(MInt n = 0; n < nDim; n++) {
#ifdef WAR_NVHPC_PSTL
      m_solver->a_variable(pCellId, n) = l_uu[n];
#else
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
#endif
    }
  });
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc10001, bc10001, 1001}
 * \author Georg Eitel-Amor, Andreas Lintermann
 * \date 02.10.2012
 *
 * Lattice Boltzmann inflow boundary condition
 *
 * Diriclet condition with prescribed eq-distributions.
 * Velocity vector is read from properties and prescribed on boundary.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10001(MInt index)
 *
 */
// TODO labels:LB Remove this and update Testcases to bc10000 (2D_cylinder...)
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10001(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat ma_x_F1BCS = m_Ma * LBCS;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    const MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    MFloat l_rho = 1.0;
    MFloat l_uu[nDim] = {0.0};

    // extrapolate density
    extrapolateVariable(ind, pCellId, PV->RHO, &l_rho);

    if(nDim == 2) {
      l_rho = (m_solver->a_variable(pCellId, PV->RHO) + l_rho) * F1B2;
    }

    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * ma_x_F1BCS * (*bndCell).m_multiplier * m_zeroInflowVelocity;
    }

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, l_uu);

    //#####################################################################
    // WARNING!!!
    // Why is it not used in 2D
    IF_CONSTEXPR(nDim == 3) {
      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
      }
      m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    }
    //#####################################################################
  }
  //  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc10002, bc10002, 1002}
 * \author Andreas Lintermann
 * \date 21.04.2010
 *
 * Lattice Boltzmann inflow boundary condition
 *
 * Velocities are prescribed, rho is extrapolated from inside;
 * only valid for channel with center \f$y=0\f$.
 * Density \f$\rho\f$ is extrapolated from neighbors inside. The prescirbed profile
 * is the parabolic one for channel flows (see script Schroeder).
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10002(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10002(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm
  MFloat rho = 1.0, u[nDim];
  MFloat tmpWidth;
  MInt ind = m_mapSegIdsInOutCnd[index];
  ScratchSpace<MFloat> bBox(2 * nDim, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];
  m_solver->m_geometry->getBoundingBox(bBoxPtr);

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MInt currentId = m_bndCells[i].m_cellId;

    extrapolateVariable(ind, currentId, PV->RHO, &rho);

    // This is different to bc1000
    tmpWidth = 0.5 * fabs(bBox[nDim + 1] - bBox[1]);

    const MFloat relPos = m_solver->a_coordinate(currentId, 1) / tmpWidth;
    const MFloat linear = relPos;
    const MFloat parabola = (1.0 - relPos * relPos);

    u[0] = m_Ma * LBCS
           * ((1.0 - m_solver->m_CouettePoiseuilleRatio) * linear + m_solver->m_CouettePoiseuilleRatio * parabola);
    u[1] = 0.0;
    IF_CONSTEXPR(nDim == 3) u[2] = 0.0;

    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);
    m_solver->setEqDists(currentId, rho, squaredVelocity, u);
  }
}

/** \LBBC{sec_LBBC_bc10003, bc10003, 1003}
 * Lattice Boltzmann inflow boundary condition using the nonequilibrium PPDF
 *
 * Simple Dirichlet boundary condition u=mach number. Macroscopic values
 * are prescribed and a pseudo collision step is performed, which
 * overwrites the distributions computed by the real collision step.<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10003(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10003(MInt index) {
  TRACE();

  MFloat rho = F1;
  MFloat l_uu[nDim] = {F0};
  MInt currentId;
  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    currentId = m_bndCells[i].m_cellId;

    MFloat nu = m_solver->a_nu(currentId);
    m_omega = 2.0 / (1.0 + 6.0 * nu * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)));

    // 1. Calculate eq-Distributions from actual macroscopic values
    //--------------------------------------------

    // macroscopic values
    rho = m_solver->a_variable(currentId, PV->RHO);
    MFloat squaredVelocity = F0;
    for(MInt d = 0; d < nDim; d++) {
      l_uu[d] = m_solver->a_variable(currentId, d);
      squaredVelocity += l_uu[d] * l_uu[d];
    }

    std::array<MFloat, nDist> eDistributions{};
    sysEqn().calcEqDists(rho, l_uu, eDistributions.data());

    // 2. Calculation of non-eq-parts3
    //--------------------------------------------
    std::array<MFloat, nDist> nePart{};
    for(MInt j = 0; j < 8; j++) {
      nePart[j] = m_solver->a_distribution(currentId, j) - eDistributions[j];
    }

    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * m_Ma * F1BCS * m_bndCells[i].m_multiplier * m_zeroInflowVelocity;
    }

    // 3. Calculation of incoming distributions in bnd cell
    //--------------------------------------------

    std::array<MFloat, nDist> newEqDists{};
    sysEqn().calcEqDists(rho, l_uu, newEqDists.data());

    for(MInt j = 0; j < nDist - 1; j++) {
      m_solver->a_oldDistribution(currentId, j) = newEqDists[j] + (F1 - m_omega) * nePart[j];
    }
  }
}

/** \brief Lattice Boltzmann inflow boundary condition
 *
 * \author Andreas Lintermann, Gabriel Faustini
 * \date 26.04.2023
 *
 * Velocities are prescribed, rho is extrapolated from inside;
 * only valid for channel with center \f$y=0\f$.
 * Density \f$\rho\f$ is extrapolated from neighbors inside. The prescirbed profile
 * is the parabolic one for channel flows (see script Schroeder).
 *
 * Edited from bc1002 to accept non-Newtonian fluid
 *
 * \param[in] index the index of the BC
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10004(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  std::array<MFloat, 2 * nDim> bBox{};
  m_solver->m_geometry->getBoundingBox(bBox.data());
  const MFloat channelHalfHeight = fabs(bBox[nDim + 1] - bBox[1]) * F1B2;

  const MFloat exp = F1 + F1 / m_solver->m_n;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MInt currentId = m_bndCells[i].m_cellId;

    MFloat rho = 1.0;
    extrapolateVariable(ind, currentId, PV->RHO, &rho);

    const MFloat relPos = fabs(m_solver->a_coordinate(currentId, 1) / channelHalfHeight - F1);
    const MFloat parabola = F1 - pow(relPos, exp);

    std::array<MFloat, nDim> u{m_Ma * LBCS * parabola};

    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);
    m_solver->setEqDists(currentId, rho, squaredVelocity, u.data());

    // Saving bc results in m_solver
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(currentId, PV->VV[n]) = u[n];
    }
    m_solver->a_variable(currentId, PV->RHO) = rho;
  }
}

/** \LBBC{sec_LBBC_bc10022, bc10022, 1022}
 * \author Moritz Waldmann
 * \date 10.12.2019
 *
 * Lattice Boltzmann inflow boundary condition similar to 10002 but with thermal treatment
 *
 * Velocities are prescribed, rho is extrapolated from inside;
 * only valid for channel with center \f$y=0\f$.
 * Density \f$\rho\f$ is extrapolated from neighbors inside. The prescirbed profile
 * is a combination of Coutte and Poiseuille flow. The temperature is set to ambient temperature.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10022(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10022(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm
  const MFloat T = 1.0;
  MInt ind = m_mapSegIdsInOutCnd[index];
  std::array<MFloat, 2 * nDim> bBox{};
  m_solver->m_geometry->getBoundingBox(bBox.data());

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MInt currentId = m_bndCells[i].m_cellId;

    MFloat rho = 1.0;
    extrapolateVariable(ind, currentId, PV->RHO, &rho);

    // This is different to bc1000
    // Channel height
    const MFloat tmpWidth = fabs(bBox[nDim + 1] - bBox[1]);
    // WARNING, I dont understand this, but without is the testcase does not work
    const MFloat bottom = bBox[1];

    const MFloat relPos = (m_solver->a_coordinate(currentId, 1) - bottom) / tmpWidth;
    const MFloat parabola = (relPos - relPos * relPos);

    const MFloat linear = relPos;
    std::array<MFloat, nDim> l_uu{};
    l_uu[0] = m_Ma * LBCS
              * ((1.0 - m_solver->m_CouettePoiseuilleRatio) * linear + m_solver->m_CouettePoiseuilleRatio * parabola)
              * rho;
    l_uu[1] = 0.0;
    IF_CONSTEXPR(nDim == 3) { l_uu[2] = 0.0; }

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    std::array<MFloat, nDim> u{};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = l_uu[d] / rho;
    }
    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // TODO: This is not correct. It should be changed in a separat branch, since the modification of testcases is
    // required
    m_solver->setEqDists(currentId, rho, squaredVelocity, u.data());
    m_solver->setEqDistsThermal(currentId, T, rho, squaredVelocity, u.data());

    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(currentId, PV->VV[n]) = l_uu[n];
    }
    m_solver->a_variable(currentId, PV->RHO) = rho;
    m_solver->a_variable(currentId, PV->T) = T;
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc20030, bc20030, 2030}
 *
 * Haenel wall bc with non zero velocity<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc20030(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc20030(MInt index) {
  TRACE();

  MFloat F2q, rhoF, uF[nDim], uBF[nDim], tmpUF, tmpUBF, tmp2UF, b[2 * nDim], tmpUW, c[2 * nDim];
  MFloat uW[nDim] = {F0};
  getBoundaryVelocity(index, uW);

  MInt id, k, tpIndex;
  MFloat fEq, xi;
  MInt nghbrId;

  for(MInt d = 0; d < nDim; d++) {
    c[2 * d] = -uW[d];
    c[2 * d + 1] = uW[d];
  }

  // go through m_bndCells with wallbndcnd
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    //    if((globalTimeStep - 1 ) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0 )
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    id = m_bndCells[i].m_cellId;

    m_omega = 2.0 / (1.0 + 6.0 * m_solver->a_nu(id) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(id)));

    rhoF = 1.0; // standard value which is usually overwritten in the following. Only in exceptional case it remains 1

    //----------------------------------------------------------------------------------------------------
    // perform bounce back for distributions piercing the wall - distinquish between inner and outer cells
    //----------------------------------------------------------------------------------------------------

    //---------------------------------------------------------------------------------
    // 1. boundary cell is inside fluid
    //---------------------------------------------------------------------------------
    if(m_bndCells[i].m_isFluid) {
      if(m_densityFluctuations)
        rhoF = 1.0 + m_solver->a_variable(id, PV->RHO);
      else
        rhoF = m_solver->a_variable(id, PV->RHO);

      for(MInt n = 0; n < nDim; n++) {
        uF[n] = m_solver->a_variable(id, PV->VV[n]);
        b[2 * n] = -uF[n];
        b[2 * n + 1] = uF[n];
      }

      tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

      for(MInt j = 0; j < nDist - 1; j++) {
        if(j < Ld::distFld(0)) {
          tmpUF = b[j];
          tmpUW = c[j];
          tpIndex = 1;
        } else {
          if(j < (Ld::distFld(0) + Ld::distFld(1))) {
            k = j - Ld::distFld(0);
            tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
            tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          } else {
            k = j - (Ld::distFld(0) + Ld::distFld(1));
            tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
            tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
            tpIndex = 3;
          }
        }

        //---------------------------------------------------------------------------------
        // if q < 0.5, perform bounce back
        //---------------------------------------------------------------------------------

        if(m_bndCells[i].m_distances[j] < 0.5) {
          F2q = F2 * m_bndCells[i].m_distances[j];

          if(m_solver->a_hasNeighbor(id, Ld::oppositeDist(j)) > 0) {
            nghbrId = m_solver->c_neighborId(id, Ld::oppositeDist(j));

            tmpUBF = F0;
            for(MInt d = 0; d < nDim; d++) {
              uBF[d] = m_solver->a_variable(nghbrId, d);
              tmpUBF += (Ld::idFld(j, d) - 1) * uBF[d];
            }
          } else {
            tmpUBF = tmpUF;
          }

          fEq = Ld::tp(tpIndex) * rhoF
                * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
          xi = m_omega * (F2q - 1.0) / (1.0 - 2.0 * m_omega);
          // xi = m_omega * (F2q - 1.0) / (1.0 - m_omega);

          // perform interpolated bounce back
          m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) =
              (1.0 - xi) * m_solver->a_distribution(id, j) + xi * fEq
              - 2.0 * Ld::tp(tpIndex) * rhoF * F1BCSsq * tmpUW; // the minus corresponds to the opposite direction

        }

        // if the cell is located at a corner, there might be distributions missing
        else {
          // In this case simple bounce back is used. This could decrease the overall accuracy.
          if(m_solver->a_hasNeighbor(id, j) == 0) {
            m_solver->a_oldDistribution(id, Ld::oppositeDist(j)) = m_solver->a_distribution(id, j);
          }
        }

      } // end of loop over all directions

    }

    //---------------------------------------------------------------------------------
    // 2. boundary cell is outside fluid
    // if the inner neighbor is a halo cell, it can be skipped
    //---------------------------------------------------------------------------------

    // The outer cell does not belong to the flow field, thus the macroscopic values are held constant.
    // It is possible that there are cells without any cutting velocity! Those have to be controlled, too.
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        //----------------------------------------------------------------------------------------------------------------
        // if q_innerNeighbor = 1 - q_bndCell >= 0.5, perform bounceback on inner neighbor
        //----------------------------------------------------------------------------------------------------------------
        if(m_solver->a_hasNeighbor(id, j) > 0 && !m_solver->a_isHalo(m_solver->c_neighborId(id, j))) {
          if(m_bndCells[i].m_distances[j] <= 0.5) {
            F2q = F2 - F2 * m_bndCells[i].m_distances[j]; // q = 1 - 0.5*(distance/cellHalfLength);

            nghbrId = m_solver->c_neighborId(id, j);

            if(m_densityFluctuations)
              rhoF = 1.0 + m_solver->a_variable(nghbrId, PV->RHO);
            else
              rhoF = m_solver->a_variable(nghbrId, PV->RHO);

            for(MInt n = 0; n < nDim; n++) {
              uF[n] = m_solver->a_variable(nghbrId, PV->VV[n]);
              b[2 * n] = -uF[n];
              b[2 * n + 1] = uF[n];
            }

            tmp2UF = std::inner_product(&uF[0], &uF[nDim], &uF[0], .0);

            if(j < Ld::distFld(0)) {
              tmpUF = b[Ld::oppositeDist(j)];
              tmpUW = c[Ld::oppositeDist(j)];
              tpIndex = 1;
            } else {
              if(j < (Ld::distFld(0) + Ld::distFld(1))) {
                k = Ld::oppositeDist(j) - Ld::distFld(0);
                tmpUF = (b[Ld::mFld1(2 * k)] + b[Ld::mFld1(2 * k + 1)]);
                tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
                tpIndex = 2;
              } else {
                k = Ld::oppositeDist(j) - (Ld::distFld(0) + Ld::distFld(1));
                tmpUF = (b[Ld::mFld2(3 * k)] + b[Ld::mFld2(3 * k + 1)] + b[Ld::mFld2(3 * k + 2)]);
                tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
                tpIndex = 3;
              }
            }

            tmpUBF = 0.0;
            for(MInt dim = 0; dim < nDim; dim++) {
              uBF[dim] = (1.0 - 3.0 / F2q) * uF[dim] + (3.0 / F2q) * uW[dim];
              tmpUBF += (Ld::idFld(Ld::oppositeDist(j), dim) - 1) * uBF[dim];
            }

            fEq = Ld::tp(tpIndex) * rhoF
                  * (1.0 + tmpUBF * F1BCSsq + tmpUF * tmpUF * F1BCSsq * F1BCSsq * F1B2 - tmp2UF * F1BCSsq * F1B2);
            xi = m_omega * (F2q - 1.0) / (1.0 + 0.5 * m_omega);
            // xi = m_omega * (F2q - 1.0);

            // perform interpolated bounce back
            m_solver->a_oldDistribution(nghbrId, j) =
                (1.0 - xi) * m_solver->a_distribution(nghbrId, Ld::oppositeDist(j)) + xi * fEq
                - 2.0 * Ld::tp(tpIndex) * rhoF * F1BCSsq * tmpUW; // the minus corresponds to the opposite direction
          }
        }

      } // end of loop over all directions
      m_solver->setEqDists(id, F1, uW);
    }
  } // end of loop over all m_bndCells
  if(m_calcWallForces) calculateWallForces(index);
}

/** \brief Calculation of the wall forces for walls (single solver)
 *
 * \author Moritz Waldmann
 * \date 01.05.2020
 *
 * Based on the momentum exchange of the PPDFs and the wall
 *
 * \param[in] index the index of the BC
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateWallForces(MInt index) {
  TRACE();

  if(globalTimeStep % m_calcWallForcesInterval != 0) return;

  std::array<MFloat, nDim> forceWall_;
  forceWall_.fill(0.0);

  MFloat uW[nDim] = {F0};
  getBoundaryVelocity(index, uW);

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt pCellId = m_bndCells[i].m_cellId;

    if(m_solver->a_level(pCellId) != m_solver->maxLevel()) continue;

    // case 1: boundary cell is inside fluid
    if(m_bndCells[i].m_isFluid) {
      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
        continue;
      }

      for(MInt j = 0; j < nDist - 1; j++) {
        const MInt opposite = Ld::oppositeDist(j);

        MFloat Vin[nDim], Vout[nDim];
        for(MInt d = 0; d < nDim; d++) {
          Vin[d] = Ld::ppdfDir(j, d) - uW[d];
          Vout[d] = Ld::ppdfDir(opposite, d) - uW[d];
        }

        // if q < 1.0, bounce-back was performed
        if(m_bndCells[i].m_distances[j] < 1.0) {
          if(m_solver->a_hasNeighbor(pCellId, opposite)) {
            for(MInt d = 0; d < nDim; d++) {
              forceWall_[d] += m_solver->a_distribution(pCellId, j) * Vin[d]
                               - m_solver->a_oldDistribution(pCellId, opposite) * Vout[d];
            }
          }
        }
      }
    }
    // case 2: boundary cell is inside solid
    else {
      for(MInt j = 0; j < nDist - 1; j++) {
        // do we have a neighbor at all?
        if(m_solver->a_hasNeighbor(pCellId, j)) {
          // make sure that the neighbor is not a halo cell
          const MInt neighborId = m_solver->c_neighborId(pCellId, j);
          if(m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
            continue;
          }
          const MInt opposite = Ld::oppositeDist(j);

          MFloat Vin[nDim], Vout[nDim];
          for(MInt d = 0; d < nDim; d++) {
            Vin[d] = Ld::ppdfDir(opposite, d) - uW[d];
            Vout[d] = Ld::ppdfDir(j, d) - uW[d];
          }

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          const MFloat q = m_bndCells[i].m_distances[j];
          if(q < 0.5) {
            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            const MInt nbndid = m_solver->a_bndId(neighborId);
            if(nbndid == -1 || m_bndCells[nbndid].m_isFluid) {
              // Add differences of PPDFS before and after collision to momentum sum
              for(MInt d = 0; d < nDim; d++) {
                forceWall_[d] += m_solver->a_distribution(neighborId, opposite) * Vin[d]
                                 - m_solver->a_oldDistribution(neighborId, j) * Vout[d];
              }
            }
          }
        }
      }
    }
  }

  std::array<MFloat, nDim> forceWall;
  forceWall.fill(0.0);

  const MInt segId = this->m_mapIndex2BndCndSegId[index];
  const auto& cwfc = this->m_mapWallForceContainer[segId];
  // communication if required
  if(cwfc.noComm > 1) {
    MPI_Reduce(forceWall_.data(), forceWall.data(), nDim, MPI_DOUBLE, MPI_SUM, 0, cwfc.comm, AT_, "forceWall_",
               "forceWall");
  } else {
    forceWall = forceWall_;
  }

  // output
  if(cwfc.isRoot) {
    std::FILE* forceFile;
    forceFile = fopen(cwfc.fileName.c_str(), "a+");
    fprintf(forceFile, "%10d\t", globalTimeStep);
    for(MInt d = 0; d < nDim; d++) {
      fprintf(forceFile, "%15e\t", forceWall[d]);
    }
    fprintf(forceFile, "\n");
    fclose(forceFile);
  }
}

/** \brief Calculation of the wall forces for moving walls described by a level set
 *
 * Based on the momentum exchange of the PPDFs and the wall
 *
 * \author Moritz Waldmann
 * \date 01.05.2020
 *
 * \param[in] set Index of the level set
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::calculateWallForcesMb(MInt set) {
  TRACE();
  ASSERT(set < m_solver->m_maxNoSets, "ERROR: wrong set is choosen");

  if(globalTimeStep % m_calcWallForcesInterval != 0) return;

  std::array<MFloat, nDim> forceWall_;
  forceWall_.fill(0.0);

  for(MInt cellIndex = 0; cellIndex < m_boundaryCellsMb.size(); cellIndex++) {
    const MInt pCellId = m_boundaryCellsMb.cellId(cellIndex);
    MFloat uW[nDim] = {F0};
    getBoundaryVelocityMb(cellIndex, uW);

    // case 1: boundary cell is inside fluid
    if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
        continue;
      }
      for(MInt j = 0; j < nDist - 1; j++) {
        const MInt opposite = Ld::oppositeDist(j);

        MFloat Vin[nDim], Vout[nDim];
        for(MInt d = 0; d < nDim; d++) {
          Vin[d] = Ld::ppdfDir(j, d) - uW[d];
          Vout[d] = Ld::ppdfDir(opposite, d) - uW[d];
        }

        const MFloat q = getDistanceMb(pCellId, cellIndex, j);

        if(q <= 0.5) {
          if(m_solver->a_hasNeighbor(pCellId, opposite)) {
            for(MInt d = 0; d < nDim; d++) {
              forceWall_[d] += m_solver->a_distribution(pCellId, j) * Vin[d]
                               - m_solver->a_oldDistribution(pCellId, opposite) * Vout[d];
            }
          }
        }
      }
    } else {
      for(MInt j = 0; j < nDist - 1; j++) {
        // do we have a neighbor at all?
        if(m_solver->a_hasNeighbor(pCellId, j)) {
          // make sure that the neighbor is not a halo cell
          const MInt neighborId = m_solver->c_neighborId(pCellId, j);
          if(m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
            continue;
          }
          const MInt opposite = Ld::oppositeDist(j);

          MFloat Vin[nDim], Vout[nDim];
          for(MInt d = 0; d < nDim; d++) {
            Vin[d] = Ld::ppdfDir(opposite, d) - uW[d];
            Vout[d] = Ld::ppdfDir(j, d) - uW[d];
          }

          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          const MFloat q = getDistanceMb(pCellId, cellIndex, j);
          if(q < 0.5) {
            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));
            if((!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0)) {
              // Add differences of PPDFS before and after collision to momentum sum
              for(MInt d = 0; d < nDim; d++) {
                forceWall_[d] += m_solver->a_distribution(neighborId, opposite) * Vin[d]
                                 - m_solver->a_oldDistribution(neighborId, j) * Vout[d];
              }
            }
          }
        }
      }
    }
  }

  std::array<MFloat, nDim> forceWall;
  forceWall.fill(0.0);

  // communication if required
  MBool isRoot = true;
  if(m_solver->noDomains() > 1) {
    MPI_Reduce(forceWall_.data(), forceWall.data(), nDim, MPI_DOUBLE, MPI_SUM, 0, m_BCWallMBComm, AT_, "forceWall_",
               "forceWall");
    isRoot = m_solver->domainId() == m_BCWallMBNeighbors[0];
  } else {
    forceWall = forceWall_;
  }

  // output
  if(isRoot) {
    std::FILE* forceFile;
    forceFile = fopen(m_forceFile.c_str(), "a+");
    fprintf(forceFile, "%10d\t", globalTimeStep);
    for(MInt d = 0; d < nDim; d++) {
      fprintf(forceFile, "%15e\t", forceWall[d]);
    }
    // drag coefficient for a sphere
    // const MFloat factor = 1 / (0.5 * m_Ma * m_Ma * F1B3 * PI * m_referenceLength * m_referenceLength * 0.25);
    // for(MInt d = 0; d < nDim; d++) {
    //  fprintf(forceFile, "%e\t", forceWall[d]*factor);
    // }
    fprintf(forceFile, "\n");
    fclose(forceFile);
  }
}

/** \LBBC{sec_LBBC_bc10050, bc10050, 1050}
 * \author Andreas Lintermann
 * \date 28.04.2015
 *
 * Lattice Boltzmann inflow boundary condition
 *
 * Diriclet condition with prescribed eq-distributions.
 * As a velocity vector the reference velocity is used and prescribed on boundary for direction X.
 * Density is set to rho = 1.0.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10050(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10050(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    const MFloat l_rho = 1.0;
    MFloat l_uu[nDim] = {0.0};
    l_uu[0] = m_Ma * LBCS;

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    // calulate the equilibrium distribution functions
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
}

/** \LBBC{sec_LBBC_bc10060, bc10060, 1060}
 * \author Andreas Lintermann, Moritz Waldmann
 * \date 23.02.2011
 *
 * Inlet boundary condidtion similar to bc10000 including Thermal boundary treatment (also available for inner
 * energy distribution function)
 *
 * Temperature is prescribed on boundary and eq-PPFDs are evaluated.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10060(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10060(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  const MFloat ma_x_F1BCS = m_Ma * LBCS;
  const MInt pvrho = PV->RHO;
  const MInt pvt = PV->T;

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  const MInt begin = m_bndCndOffsets[index];
  const MInt end = m_bndCndOffsets[index + 1];
  const MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    const MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    MFloat l_rho = 1.0;

    // extrapolate density
    extrapolateVariable(ind, pCellId, pvrho, &l_rho);

    std::array<MFloat, nDim> l_uu{};
    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * ma_x_F1BCS * (*bndCell).m_multiplier * m_zeroInflowVelocity * l_rho;
    }

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    std::array<MFloat, nDim> u{};
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, u.data());
    MFloat l_t = 1.0;
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, squaredVelocity, u.data());

    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, n) = l_uu[n];
    }
    m_solver->a_variable(pCellId, pvrho) = l_rho;
    m_solver->a_variable(pCellId, pvt) = l_t;
#ifdef WAR_NVHPC_PSTL
  });
#else
  });
#endif
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc10061, bc10061, 1061}
 * \author Andreas Lintermann
 * \date 02.10.2012
 *
 * Inlet boundary condidtion similar to bc10000 including Thermal boundary treatment. Blasius profile is
 * prescribed.
 *
 * A Blasius solution is prescribed for the velocity and temperature.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10061(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10061(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  const MInt x_pos = m_solver->m_referenceLength * m_solver->m_blasiusPos;
  const MInt pvrho = PV->RHO;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    const MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    // Use Blasius
    const MFloat eta = (*bndCell).m_eta;
    const MInt pos = (eta / m_blasius_delta);

    MFloat l_rho = 1.0;

    // extrapolate density
    extrapolateVariable(ind, pCellId, pvrho, &l_rho);

    std::array<MFloat, nDim> l_u{};
    std::array<MFloat, nDim> u{};
    const MFloat l_t = m_solver->m_initTemperatureKelvin - m_blasius[pos][2];
    l_u[0] = (m_Ma * LBCS) * m_blasius[pos][2] * l_rho;
    IF_CONSTEXPR(nDim == 3) l_u[1] = 0.0 * l_rho;
    l_u[nDim - 1] =
        0.5 * sqrt((m_Ma * LBCS) * m_solver->m_nu / x_pos) * (eta * m_blasius[pos][2] - m_blasius[pos][1]) * l_rho;

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_u[n] / l_rho;
    }
    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, squaredVelocity, u.data());

    m_solver->a_variable(pCellId, PV->U) = l_u[0];
    m_solver->a_variable(pCellId, PV->V) = l_u[1];
    IF_CONSTEXPR(nDim == 3) m_solver->a_variable(pCellId, PV->W) = l_u[2];
    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    m_solver->a_variable(pCellId, PV->T) = l_t;
  }
}

/** \brief Lattice Boltzmann inflow boundary condition similar to BC1060
 * but for concentration transport. Inflow concentration is set to 0.
 *
 * \author Shota Ito
 * \date 08.06.22
 *
 * \param[in] index the index of the BC
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10090(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  const MFloat ma_x_F1BCS = m_Ma * LBCS;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    const MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    MFloat l_rho = 1.0;
    const MFloat l_t = 1.0;

    // extrapolate density
    extrapolateVariable(ind, pCellId, PV->RHO, &l_rho);

    MFloat l_c = 0.0;
    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> u{};
    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * ma_x_F1BCS * (*bndCell).m_multiplier * m_zeroInflowVelocity * l_rho;
    }

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // calculate the equlibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, squaredVelocity, u.data());
    m_solver->setEqDistsTransport(pCellId, l_c, squaredVelocity, u.data());

    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    if(m_solver->m_isThermal) m_solver->a_variable(pCellId, PV->T) = l_t;
    m_solver->a_variable(pCellId, PV->C) = l_c;
  }
  writeBCOutput(index);
}

/**\LBBC{sec_LBBC_bc10070, bc10070, 1070}
 * \author Miro Gondrum
 * \date 22.12.2022
 *
 * Neumann density BC adjusting the local mass flux to be idential to rho_0*Ma*Cs
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10070(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10070(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  constexpr MFloat rho0 = 1.0;
  const MFloat trgRhoU = rho0 * m_Ma * LBCS;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    auto& bndCell = m_bndCells[i];
    const MInt cellId = bndCell.m_cellId;
    MFloat rho = 1.0;
    extrapolateVariable(ind, cellId, PV->RHO, &rho);
    const MFloat trgU = trgRhoU / rho;
    MFloat u[nDim];
    for(MInt n = 0; n < nDim; n++) {
      u[n] = trgU * m_initialVelocityVecs[ind][n];
    }
    m_solver->setEqDists(cellId, rho, u);
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(cellId, PV->U + n) = u[n];
    }
    m_solver->a_variable(cellId, PV->RHO) = rho;
  }
  // writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc10080, bc10080, 1080}
 * \author Andreas Lintermann
 * \date 30.05.2012, 28.04.2015
 *
 * Lattice Boltzmann inflow boundary condition for pulsatile flow
 *
 * Diriclet condition with prescribed eq-distributions, the
 * velocity vector is read from properties and prescribed on boundary.
 * Volume flux varies over time sinusodially.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10080(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10080(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat pulsatile = -1.0 * sin(2.0 * PI * m_pulsatileFrequency * globalTimeStep);
  MFloat ma_x_pulsatile_F1BCS = m_Ma * pulsatile * LBCS;
  const MInt pvrho = PV->RHO;

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    LbGridBoundaryCell<nDim>* bndCell = &m_bndCells[i];
    const MInt currentId = (*bndCell).m_cellId;
    const MInt pCellId = currentId;

    MFloat l_rho = 1.0;
    MFloat l_uu[nDim] = {0.0};

    // extrapolate density
    extrapolateVariable(ind, pCellId, pvrho, &l_rho);

    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_initialVelocityVecs[ind][n] * (*bndCell).m_multiplier * ma_x_pulsatile_F1BCS;
    }

    // TODO: This is not correct. It should be changed in a separat branch, since the modification of testcases is
    // required
    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] *= l_rho;
    }
    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);

#ifdef WAR_NVHPC_PSTL
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, n) = l_uu[n];
    }
#else
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
#endif
    m_solver->a_variable(pCellId, pvrho) = l_rho;
#ifdef WAR_NVHPC_PSTL
  });
#else
  });
#endif
  writeBCOutput(index);
}


/** \LBBC{sec_LBBC_bc10111, bc10111, 1111}
 * Povitsky cavity flow
 *
 * Boundary cells are set to equilibrium with initialVelocityVecs.
 *
 *  \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10111(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10111(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm

  MInt ind = m_mapSegIdsInOutCnd[index];

  const MFloat velMag = sqrt(std::inner_product(&m_initialVelocityVecs[ind][0], &m_initialVelocityVecs[ind][nDim],
                                                &m_initialVelocityVecs[ind][0], 0.0));
  const MFloat rho = 1.0;

  std::array<MFloat, nDim> u{};
  for(MInt n = 0; n < nDim; n++) {
    u[n] = m_Ma * F1BCS * m_initialVelocityVecs[ind][n] / velMag;
  }

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt cellId = m_bndCells[i].m_cellId;

    m_solver->setEqDists(cellId, rho, u.data());

    m_solver->a_oldVariable(cellId, PV->RHO) = rho;
    m_solver->a_variable(cellId, PV->RHO) = rho;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_oldVariable(cellId, PV->VV[n]) = u[n];
      m_solver->a_variable(cellId, PV->VV[n]) = u[n];
    }
  }
}


/** \LBBC{sec_LBBC_bc40000, bc40000, 4000}
 * \author Rainhill Freitas, Andreas Lintermann
 * \date 28.04.2015
 *
 * Lattice Boltzmann inflow boundary condition
 *
 * Diriclet condition for the density rho, which is read from the properties file. This implements
 * the non-reflecting boundary condition from
 * Finck, Haenel, (2007), https://doi.org/10.1016/j.compbiomed.2006.06.013,
 * according to either
 *
 * \f[ \rho = \frac{1}{2} \rho_{old} + \frac{1}{c_{s}} * ||\vec{v}|| - ||\vec{v}_{old}|| + (\rho_{prop}-1) \f]
 *
 * or
 *
 *  \f[ \rho = \frac{1}{2} \rho_{old} + \frac{1}{c_{s}} * ||\vec{v}|| - ||\vec{v}_{old}|| + \rho_{prop} \f]
 *
 * depending on the property m_densityFluctuations. Additionally a von Neumann condition is applied for the
 * velcoity. Note that in the rare case that the extrapolation directions have not been found (-1) or the
 * according neighbors in all extrapolation directions are missing, the velocity is set to 0. This was a wrong
 * implemnetation in the previous version of this algorithm, which set the velocity to some arbitrary value
 * still stored in memory from a previous calculation.
 *
 * \param[in] index the index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40000(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40000(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];
  const MInt pvrho = PV->RHO;

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  const MInt begin = m_bndCndOffsets[index];
  const MInt end = m_bndCndOffsets[index + 1];
  const MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    MFloat l_uu[nDim] = {0.0};

    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    // non-reflecting bc (Finck, Haenel 2008)
    MFloat l_rho = 1.0;
    // TODO labels:LB Remove the if constexpr(nDim) and update Testcases (2D_channel...)
    if constexpr(nDim == 3) {
      MFloat old_squared_velocity = F0;
      MFloat squaredVelocity = F0;
      for(MInt n = 0; n < nDim; n++) {
        squaredVelocity += (l_uu[n] * l_uu[n]);
        const MFloat l_old_u = m_solver->a_oldVariable(pCellId, n);
        old_squared_velocity += (l_old_u * l_old_u);
      }
      const MFloat old_rho = m_solver->a_oldVariable(pCellId, pvrho);
      if(m_densityFluctuations)
        l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + (m_rho1 - 1.0)) / 2.0;
      else
        l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + m_rho1) / 2.0;
    } else {
      l_rho = (1.0 + m_solver->a_variable(pCellId, pvrho)) * F1B2;
    }

    m_solver->setEqDists(pCellId, l_rho, l_uu);

    m_solver->a_variable(pCellId, pvrho) = l_rho;
#ifdef WAR_NVHPC_PSTL
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, n) = l_uu[n];
    }
#else
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
#endif
  });
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40020, bc40020, 4020}
 *
 * Enforce density at boundary (Dirichlet condition)
 *
 * The mean density at inflow is relaxed against the fixed
 * density value.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40020(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40020(MInt index) {
  TRACE();

  MInt cellsWeigthed = 0;
  MInt weight = 0;
  // For the calculation of the mean density the different sizes of the
  // cells must be considered.
  MInt maxLevel = m_solver->a_level(m_bndCells[m_bndCndOffsets[index]].m_cellId);
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    maxLevel =
        (maxLevel > m_solver->a_level(m_bndCells[i].m_cellId)) ? maxLevel : m_solver->a_level(m_bndCells[i].m_cellId);
  }

  MFloat meanDensity = 0.0;
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    weight = pow(IPOW2(maxLevel - m_solver->a_level(m_bndCells[i].m_cellId)), nDim);
    if(m_solver->c_noChildren(m_bndCells[i].m_cellId) == 0) {
      for(MInt j = 0; j < nDist; j++) {
        meanDensity += weight * m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j);
      }
      cellsWeigthed += weight;
    }
  }
  meanDensity = 1.0 - meanDensity / cellsWeigthed;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    for(MInt j = 0; j < nDist; j++) {
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) += Ld::tp(Ld::distType(j)) * meanDensity;
    }
  }
}

/** \brief Lattice Boltzmann inflow boundary condition
 * TODO labels:LB called nowhere
 * Do nothing condition in axis direction: only for LBGK!!! (rhs)
 *
 * \param[in] index     Boundary index
 * \param[in] direction Boundary orientation in Cartesian direction
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::dnt(MInt index, MInt direction) {
  TRACE();

  if(direction > 2 * nDim) {
    TERMM(1, "invalid direction for dnt!");
  }

  MFloat rho, u[nDim], uInner[nDim];
  MInt currentId, currentDist;

  MFloatScratchSpace eDistributions(nDist, AT_, "eDistributions");
  MFloatScratchSpace eDistributions1(nDist, AT_, "eDistributions1");
  MFloat M;
  MInt neighborId;

  const MInt d = std::floor(direction / 2.0);

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    // since it is a rhsBnd it should be performed after propagation
    // TODO labels:LB,totest Not tested yet for a locally refined mesh

    currentId = m_bndCells[i].m_cellId;

    // leave out halo cells
    if(m_solver->a_isHalo(currentId)) continue;

    neighborId = m_solver->c_neighborId(currentId, Ld::oppositeDist(direction));

    m_omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)));

    // extrapolate inner velocity
    for(MInt dim = 0; dim < nDim; dim++) {
      uInner[dim] = m_solver->a_variable(neighborId, PV->VV[dim]);
    }

    MFloat zeroVel[nDim] = {0.0};

    sysEqn().calcEqDists(1.0, 0.0, zeroVel, &eDistributions1[0]);

    // prescribe velocity (extrapolated) via BFL rule
    for(MInt j = 0; j < Ld::dxQyFld(); j++) {
      currentDist = Ld::componentFld(Ld::oppositeDist(direction), j);

      m_solver->a_oldDistribution(currentId, currentDist) =
          m_solver->a_oldDistribution(currentId, Ld::oppositeDist(currentDist));

      for(MInt k = 0; k < nDim; k++) {
        m_solver->a_oldDistribution(currentId, currentDist) +=
            6 * eDistributions1[currentDist] * (Ld::idFld(currentDist, k) - 1.0) * uInner[k];
      }
    }

    // calculate eq-distributions from actual variables

    rho = m_solver->a_variable(m_bndCells[i].m_cellId, PV->RHO);
    for(MInt dim = 0; dim < nDim; dim++) {
      u[dim] = m_solver->a_variable(m_bndCells[i].m_cellId, PV->VV[dim]);
    }

    sysEqn().calcEqDists(rho, u, &eDistributions[0]);

    // ----------------------------------
    // calculate diagonal Matrix M from eq-distributions
    // by now only for BGK-collisions!!!
    // ----------------------------------

    M = (Ld::idFld(direction, d) - 1.0) * (-m_solver->a_nu(m_bndCells[i].m_cellId) * m_omega)
        * (m_solver->a_oldDistribution(m_bndCells[i].m_cellId, direction) - eDistributions[direction]);

    // calculate eq-distributions for rho = 1 and actual velocity
    sysEqn().calcEqDists(1.0, u, &eDistributions1[0]);

    // ---------------------------------
    // overwrite incoming distribution along coordinate axes
    // ---------------------------------
    if(m_solver->a_hasNeighbor(currentId, direction) == 0) {
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, Ld::oppositeDist(direction)) =
          eDistributions1[Ld::oppositeDist(direction)] + M
          + (m_solver->a_oldDistribution(m_bndCells[i].m_cellId, direction) - eDistributions[direction]);
    }
  }
}

/** \LBBC{sec_LBBC_bc40030, bc40030, 4030}
 * Guo, Extrapolation of velocity and non-eq parts
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40030(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40030(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0) continue;

    // leave out halo cells
    if(m_solver->a_isHalo(currentId)) continue;

    // macroscopic values
    const MFloat rho = m_solver->a_variable(currentId, PV->RHO);
    std::array<MFloat, nDim> u;
    std::array<MFloat, nDim> oldU;
    for(MInt d = 0; d < nDim; d++) {
      u[d] = m_solver->a_variable(currentId, PV->VV[d]);
      oldU[d] = m_solver->a_oldVariable(currentId, PV->VV[d]);
    }
    const MFloat squaredVelocity = std::inner_product(u.begin(), u.end(), u.begin(), 0.0);

    // WARNING!!! The density should be calculated as in 3d
    // However, result is correct
    // TODO labels:LB Remove the if and update Testcases (2D_cylinder...)
    MFloat rhoNonRefl;
    if constexpr(nDim == 3) {
      // non-reflecting bc (Finck, Haenel 2007, Eq.(12-13))
      const MFloat oldSquaredVelocity = std::inner_product(oldU.begin(), oldU.end(), oldU.begin(), 0.0);
      const MFloat rho_out = (m_densityFluctuations) ? 0.0 : 1.0;
      rhoNonRefl = (rho + F1BCS * (sqrt(squaredVelocity) - sqrt(oldSquaredVelocity)) + rho_out) / 2.0;
    } else {
      rhoNonRefl = (rho + 1.0) / 2.0;
    }

    // Calculate current and new equilibrium distribution
    std::array<MFloat, nDist> eqDist, eqNew;
    if(m_solver->isCompressible()) {
      sysEqn().calcEqDists(rho, squaredVelocity, u.data(), eqDist.data());
      sysEqn().calcEqDists(rhoNonRefl, squaredVelocity, u.data(), eqNew.data());
    } else {
      sysEqn().calcEqDists(rho, squaredVelocity, u.data(), eqDist.data());
      sysEqn().calcEqDists(rhoNonRefl, squaredVelocity, u.data(), eqNew.data());
    }

    // Add new equilibrium and old non-equilibrium
    const MInt lvlDiff = m_solver->maxLevel() - m_solver->a_level(currentId);
    const MFloat omega = 1.0 / (0.5 + F1BCSsq * m_solver->a_nu(currentId) * FFPOW2(lvlDiff));
    for(MInt j = 0; j < nDist - 1; j++) {
      const MFloat neqDist = m_solver->a_distribution(currentId, j) - eqDist[j];
      m_solver->a_oldDistribution(currentId, j) = eqNew[j] + (1 - omega) * neqDist;
    }
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40070, bc40070, 4070}
 * \author Andreas Lintermann
 * \date 25.08.2012
 *
 * Lattice Boltzmann inflow boundary condition as proposed by Ingold Hoerschler
 *
 * Using the equation of Saint Vernant and Wanzel, one can obtain the pressure by using
 *
 * \f[p_1 = p_0\left(1 - \frac{\gamma - 1}{2\gamma}\cdot\frac{\rho_0}{p_0}\cdot v_0^2\right)^{\frac{\gamma}{\gamma
 * -1}}\f]
 *
 * A non-dimensionalization with
 *
 * \f[\rho^\ast = \frac{\rho}{\rho_0},\quad c_{s}^{\ast} = \frac{c_s}{\xi_0},\quad v^{\ast} = \frac{v}{\xi_0}\f]
 *
 * and using \f[p=\rho c_s^2\f], \f[c_s=\frac{1}{\sqrt{3}}\xi_0\f] and extending with
 * \f[\frac{\rho^{\ast^2}}{\rho^{\ast^2}}\f] leads to
 *
 * \f[\rho_{i+1}=\left(1 - \frac{\gamma - 1}{2\gamma}\cdot \frac{3}{\rho^{\ast^2}_{i}}\left(\rho^{\ast}
 * v^\ast\right)^2\right)^{\frac{\gamma}{\gamma -1}}\f]
 *
 * So, first all the velocities are extrapolated to the boundary cells and the new density ratio is calculated
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40070(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40070(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  const MInt pvrho = PV->RHO;

  const MFloat gamma = lb_gamma2;
#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    MFloat l_uu[nDim] = {0.0};

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    // extrapolate inner velocity values
    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    // okay, now recalculate density
    const MFloat l_old_rho = m_solver->a_oldVariable(pCellId, pvrho);

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    // this formula comes from Ingolf Hoerschlers dissitation
    const MFloat l_rho =
        pow((1.0 - (gamma - 1.0) * (1.0 / (2.0 * gamma)) * (1.0 / (l_old_rho * l_old_rho)) * 3.0 * squaredVelocity),
            (gamma / (gamma - 1.0)));

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);

#ifdef WAR_NVHPC_PSTL
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, n) = l_uu[n];
    }
#else
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
#endif
    m_solver->a_variable(pCellId, pvrho) = l_rho;
#ifdef WAR_NVHPC_PSTL
  });
#else
  });
#endif
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40071, bc40071, 4071}
 * \author Andreas Lintermann
 * \date 25.08.2012
 *
 * Lattice Boltzmann outflow boundary condition with prescribed pressure at outlet
 *
 * This BC holds the density ratio at a level provided by m_rho1 and extrapolates the velocities.
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40071(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40071(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MFloat l_uu[nDim] = {0.0};

    const MInt currentId = m_bndCells[i].m_cellId;

    const MInt pCellId = currentId;

    // extrapolate inner velocity values
    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    // hold density
    const MFloat l_rho = m_rho1;

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, l_uu);

    m_solver->a_variable(pCellId, PV->U) = l_uu[0];
    m_solver->a_variable(pCellId, PV->V) = l_uu[1];
    if(nDim == 3) m_solver->a_variable(pCellId, PV->W) = l_uu[2];
    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40072, bc40072, 4072}
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * Lattice Boltzmann outflow boundary condition with adjusting pressure at outlet
 *
 * This BC adjusts the density ratio to obtain a volumeflux given by the Reynolds number.
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40072(MInt index)
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40072(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat mean_velocity = 0.0;
  MInt numberOfCellsPerDomain = 0;

  // first interpolate velocities and find out mean velocity
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MFloat l_uu[nDim] = {0.0};
    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    if(m_solver->c_noChildren(pCellId) == 0 && !m_solver->a_hasProperty(currentId, Cell::IsHalo)) {
      // extrapolate inner velocity values
      extrapolateVelocities(ind, pCellId, &l_uu[0]);

      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
        l_uu[n] /= m_rhoLast;
      }
      for(MInt n = 0; n < nDim; n++) {
        mean_velocity += (m_bndNormals[ind][n] * l_uu[n]);
      }

      numberOfCellsPerDomain++;
    }
  }

  MFloat mean_velocity_all = 0.0;

  MInt totalNumberOfCells = 0;
  if(m_noBCNeighbors[m_mapBndCndIdSegId[index]] > 1) {
    MPI_Allreduce(&numberOfCellsPerDomain, &totalNumberOfCells, 1, MPI_INT, MPI_SUM,
                  m_BCComm[m_mapBndCndIdSegId[index]], AT_, "numberOfCellsPerDomain",
                  "totalNumberOfCells"); // without Halo cells

    // Get the sum by communication in my group
    MPI_Allreduce(&mean_velocity, &mean_velocity_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_velocity", "mean_velocity_all");
    mean_velocity_all /= totalNumberOfCells;
  } else {
    totalNumberOfCells = numberOfCellsPerDomain;
    mean_velocity_all = mean_velocity / totalNumberOfCells;
  }

  // okay, no we can calculate the new local Reynolds number:
  MFloat l_Re = mean_velocity_all * m_referenceLength / m_solver->m_nu;
  MFloat Re_diff = m_solver->m_Re - l_Re;
  MFloat eps = 0.01;

  // how do we have to set the local density to obtain the real Reynolds number?
  MFloat l_rho = 0.0;
  if(globalTimeStep == 1) {
    l_rho = m_rho1;
    // this defines the initial step-size
    m_deltaRho = 1.0 - l_rho;
    m_maxDeltaRho = m_deltaRho;
  } else {
    l_rho = m_rhoLast;
  }

  // get new density
  if(!(fabs(Re_diff) < eps)) {
    // we still need no increase the Reynolds number, but how far?
    // this means we have to decrease the density more
    if(Re_diff > 0) {
      // this is the value by which we have to scale the decrease of the density
      MFloat scaleRe = (m_solver->m_Re - l_Re) / (m_solver->m_Re - m_ReLast);

      m_deltaRho = fabs(scaleRe * m_deltaRho);
      if(m_deltaRho > m_maxDeltaRho) m_deltaRho = m_maxDeltaRho;
      //    else if(m_deltaRho < 0.00000001)
      //      m_deltaRho = 0.00000001;

      l_rho -= m_deltaRho;

    }
    // we have to decrease the Reynolds number again, but how far?
    // this means we have to increase the density again
    else if(Re_diff < 0) {
      // this is the value by which we have to scale the decrease of the density
      MFloat scaleRe = (l_Re - m_solver->m_Re) / (l_Re - m_ReLast);

      m_deltaRho = fabs(scaleRe * m_deltaRho);
      if(m_deltaRho > m_maxDeltaRho) m_deltaRho = m_maxDeltaRho;

      l_rho += m_deltaRho;
    }
  }

  if(this->m_calcBcResidual && m_solver->domainId() == m_BCneighbors[m_mapBndCndIdSegId[index]][0]) {
    m_BCResidualStream[m_mapBndCndIdSegId[index]]
        << globalTimeStep << " final Re: " << m_solver->m_Re << " local Re: " << l_Re
        << " delta Re: " << (m_solver->m_Re - l_Re) << " local rho: " << l_rho << " delta rho: " << m_deltaRho
        << " max delta rho: " << m_maxDeltaRho << " nu " << m_solver->m_nu << " #Cells " << totalNumberOfCells
        << std::endl;
  }
  m_rhoLast = l_rho;
  m_ReLast = l_Re;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;

    const MInt pCellId = currentId;

    MFloat l_uu[nDim] = {0.0};
    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_solver->a_variable(pCellId, PV->VV[n]);
    }
    m_solver->a_variable(pCellId, PV->RHO) = l_rho;

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, l_uu);
  }
  writeBCOutput(index);
}

/** \brief recalculates the local density for BC 4073
 *
 * \author Andreas Lintermann
 * \date 10.07.2015
 *
 * This function is only executed onece per timestep in case the interval for adaption is met or
 * a report has to be written to disk.
 *
 * \param[in] index Boundary index
 *
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::recalcRho(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];

  MFloat l_Re = m_ReLast;
  // do only adaption if we meet the interval criterion or we need to report to file
  if(globalTimeStep % m_localReCutInterval == 0 || globalTimeStep % m_localReCutReportInterval == 0) {
    // calculate local mean velocity if we have a cut with the calculation plane
    MFloat mean_velocity = 0.0;
    if(m_hasLocalReCut) {
      for(MInt i = 0; i < (MInt)m_localReCutCells.size(); i++) {
        MInt pCellId = m_localReCutCells[i];
        MFloat density = m_solver->a_variable(pCellId, PV->RHO);

        MFloat l_uu[nDim] = {0.0};

        for(MInt n = 0; n < nDim; n++) {
          m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
          l_uu[n] /= density;
        }
        for(MInt n = 0; n < nDim; n++) {
          mean_velocity += (m_bndNormals[ind][n] * l_uu[n]);
        }
      }
    }

    // Get the sum by communication in my group
    MFloat mean_velocity_all;
    MPI_Allreduce(&mean_velocity, &mean_velocity_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_velocity", "mean_velocity_all");

    mean_velocity_all /= m_totalNoBcCells[m_mapBndCndIdSegId[index]];

    // okay, now we can calculate the new local Reynolds number:
    l_Re = mean_velocity_all * m_localReCutDiameter / m_solver->m_nu;
    MFloat Re_diff = m_localReCutRe - l_Re;

    MFloat eps = 0.01;

    // update rho only if we meet the interval
    if(globalTimeStep % m_localReCutInterval == 0) {
      // shall we update the density?
      if(!(fabs(Re_diff) < eps)) {
        MFloat s = Re_diff / m_localReCutRe;

        // we still need no increase the Reynolds number, but how far?
        // this means we have to decrease the density more
        if(Re_diff > 0) {
          // this is the value by which we have to scale the decrease of the density
          MFloat scaleRe = Re_diff / (m_localReCutRe - m_ReLast);

          if(fabs(s) < m_localReCutAdpPerc) {
            m_deltaRho = fabs(scaleRe * m_deltaRho);
            if(m_deltaRho > m_maxDeltaRho) m_deltaRho = m_maxDeltaRho;
          }
          m_lRho -= m_deltaRho;

        }
        // we have to decrease the Reynolds number again, but how far?
        // this means we have to increase the density again

        else if(Re_diff < 0) {
          // this is the value by which we have to scale the decrease of the density
          MFloat scaleRe = Re_diff / (l_Re - m_ReLast);

          if(fabs(s) < m_localReCutAdpPerc) {
            m_deltaRho = fabs(scaleRe * m_deltaRho);
            if(m_deltaRho > m_maxDeltaRho) m_deltaRho = m_maxDeltaRho;
          }
          m_lRho += m_deltaRho;
        }
      }
      m_rhoLast = m_lRho;
      m_ReLast = l_Re;
    }

    if(this->m_calcBcResidual && m_solver->domainId() == m_firstBCinComm
       && globalTimeStep % m_localReCutReportInterval == 0)
      m_BCResidualStream[m_mapBndCndIdSegId[index]] << globalTimeStep << "\t" << m_localReCutRe << "\t" << l_Re << "\t"
                                                    << m_ReLast << "\t" << Re_diff << "\t" << m_lRho << "\t"
                                                    << m_rhoLast << "\t" << m_deltaRho << std::endl;
  }
}

/** \LBBC{sec_LBBC_bc40073, bc40073, 4073}
 * \author Andreas Lintermann
 * \date 28.04.2015
 *
 * Lattice Boltzmann outflow boundary condition with adjusting pressure at outlet.
 *
 * This BC adjusts the density ratio to obtain a volumeflux given by the Reynolds number.
 * In contrast the local Reynolds number is not measured at the BC location but somewhere
 * else in the flow field, i.e., at a location defined in the geometry file.
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40073(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40073(MInt index) {
  TRACE();

  if(globalTimeStep == 1) {
    m_lRho = m_rho1;
    m_deltaRho = 1.0 - m_lRho;

    /*! \page propertyPage1
    \section deltaRho
    <code>MFloat LbBndCndDxQy::m_deltaRho</code>\n
    default = <code>1.0 - m_lRho</code>\n\n
    This is the inital step size for the adaptive BC 4073.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_deltaRho = Context::getSolverProperty<MFloat>("deltaRho", m_solverId, AT_, &m_deltaRho);

    m_maxDeltaRho = m_deltaRho;
  } else
    m_lRho = m_rhoLast;


  // rho is only recalculated if the BC is called for the first time on this domain at this time step
  if(globalTimeStep != m_currentTimeStep) {
    recalcRho(index);
    m_currentTimeStep = globalTimeStep;
  }

  MInt ind = m_mapSegIdsInOutCnd[index];
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    MFloat l_uu[nDim] = {0.0};

    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    m_solver->setEqDists(pCellId, m_lRho, l_uu);

    m_solver->a_variable(pCellId, PV->RHO) = m_lRho;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40080, bc40080, 4080}
 * \author Andreas Lintermann
 * \date 01.10.2012
 *
 * This is the same as bc40070, except that it additionally prescribes an inflow temperature
 *
 * For more details refer to bc40070. The temperture is set to the reference temperature. In general, the
 * wall gets the different temperature.
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40080(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40080(MInt index) {
  TRACE();

  const MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> u{};
    // the wall-BC gets the higher temerpature, the inflow temperature is the reference temperature
    const MFloat l_t = 1.0;

    const MInt currentId = m_bndCells[i].m_cellId;

    const MInt pCellId = currentId;

    // extrapolate inner velocity values
    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    // okay, now recalculate density
    const MFloat l_old_rho = m_solver->a_oldVariable(pCellId, PV->RHO);

    const MFloat squaredVelocityRho = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    // this formula comes from Ingolf Hoerschlers dissitation
    MFloat l_rho = pow(
        (1.0
         - (lb_gamma2 - 1.0) * (1.0 / (2.0 * lb_gamma2)) * (1.0 / (l_old_rho * l_old_rho)) * 3.0 * squaredVelocityRho),
        (lb_gamma2 / (lb_gamma2 - 1.0)));

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, u.data());

    if(m_solver->m_isTransport) {
      MFloat l_c = F1;
      m_solver->setEqDistsTransport(pCellId, l_c, u.data());
      m_solver->a_variable(pCellId, PV->C) = l_c;
    }

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    m_solver->a_variable(pCellId, PV->T) = l_t;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40081, bc40081, 4081}
 * \author Andreas Lintermann
 * \date 25.08.2012
 *
 * Lattice Boltzmann outflow boundary condition with prescribed pressure at outlet
 *
 * This BC holds the density ratio at a level provided by m_rho1 and extrapolates the velocities and the temperature.
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40081(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40081(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MFloat l_t = 0.0;
    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> u{};

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    // extrapolate inner velocity values
    extrapolateVelocities(ind, pCellId, &l_uu[0]);
    extrapolateVariable(ind, pCellId, PV->T, &l_t);

    // this formula comes from Ingolf Hoerschlers dissitation
    MFloat l_rho = m_rho1;

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, squaredVelocity, u.data());

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    m_solver->a_variable(pCellId, PV->T) = l_t;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
  writeBCOutput(index);
}

/** \LBBC{sec_LBBC_bc40082, bc40082, 4082}
 * \author Moritz Waldmann
 * \date 03.05.2018
 *
 * Lattice Boltzmann outflow boundary condition with adjusting pressure and extrapolating temperature at the
 * outlet
 *
 * Similar to BC 40072, but with additional extrapolation of the temperature
 *
 * \param[in] index the index of the segment of this BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40082(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40082(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat mean_velocity = 0.0;
  MInt numberOfCellsPerDomain = 0;

  // first interpolate velocities and find out mean velocity
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    if(m_solver->c_noChildren(pCellId) == 0 && !m_solver->a_hasProperty(currentId, Cell::IsHalo)) {
      MFloat l_uu[nDim] = {0.0};

      // extrapolate inner velocity values
      extrapolateVelocities(ind, pCellId, &l_uu[0]);

      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
        l_uu[n] /= m_rhoLast;
      }

      for(MInt n = 0; n < nDim; n++) {
        mean_velocity += (m_bndNormals[ind][n] * l_uu[n]);
      }
      numberOfCellsPerDomain++;
    }
  }

  MFloat mean_velocity_all = 0.0;

  //  MInt totalNumberOfCells = m_totalNoBcCells[m_mapBndCndIdSegId[index]]; Including Halo cells
  MInt totalNumberOfCells = 0;
  if(m_noBCNeighbors[m_mapBndCndIdSegId[index]] > 1) {
    MPI_Allreduce(&numberOfCellsPerDomain, &totalNumberOfCells, 1, MPI_INT, MPI_SUM,
                  m_BCComm[m_mapBndCndIdSegId[index]], AT_, "numberOfCellsPerDomain",
                  "totalNumberOfCells"); // without Halo cells

    // Get the sum by communication in my group
    MPI_Allreduce(&mean_velocity, &mean_velocity_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_velocity", "mean_velocity_all");
    mean_velocity_all /= totalNumberOfCells;
  } else {
    totalNumberOfCells = numberOfCellsPerDomain;
    mean_velocity_all = mean_velocity / totalNumberOfCells;
  }

  // okay, no we can calculate the new local Reynolds number:
  MFloat l_Re = mean_velocity_all * m_referenceLength / m_solver->m_nu;
  MFloat Re_diff = m_solver->m_Re - l_Re;
  MFloat eps = 0.01;

  // how do we have to set the local density to obtain the real Reynolds number?
  MFloat l_rho = 0.0;
  if(globalTimeStep == 1) {
    l_rho = m_rho1;
    // this defines the initial step-size
    m_deltaRho = 1.0 - l_rho;
    m_maxDeltaRho = m_deltaRho;
  } else {
    l_rho = m_rhoLast;
  }

  // get new density
  if(!(fabs(Re_diff) < eps)) {
    // we still need no increase the Reynolds number, but how far?
    // this means we have to decrease the density more
    if(Re_diff > 0) {
      // this is the value by which we have to scale the decrease of the density
      MFloat scaleRe = (m_solver->m_Re - l_Re) / (m_solver->m_Re - m_ReLast);

      m_deltaRho = fabs(scaleRe * m_deltaRho);
      if(m_deltaRho > m_maxDeltaRho)
        m_deltaRho = m_maxDeltaRho;
      else if(m_deltaRho < 0.00000001)
        m_deltaRho = 0.00000001;

      l_rho -= m_deltaRho;

    }
    // we have to decrease the Reynolds number again, but how far?
    // this means we have to increase the density again
    else if(Re_diff < 0) {
      // this is the value by which we have to scale the decrease of the density
      MFloat scaleRe = (l_Re - m_solver->m_Re) / (l_Re - m_ReLast);

      m_deltaRho = fabs(scaleRe * m_deltaRho);
      if(m_deltaRho > m_maxDeltaRho)
        m_deltaRho = m_maxDeltaRho;
      else if(m_deltaRho < 0.00000001)
        m_deltaRho = 0.00000001;

      l_rho += m_deltaRho;
    }
  }

  if(this->m_calcBcResidual && m_solver->domainId() == m_BCneighbors[m_mapBndCndIdSegId[index]][0]) {
    m_BCResidualStream[m_mapBndCndIdSegId[index]]
        << globalTimeStep << " final Re: " << m_solver->m_Re << " local Re: " << l_Re
        << " delta Re: " << (m_solver->m_Re - l_Re) << " local rho: " << l_rho << " delta rho: " << m_deltaRho
        << " max delta rho: " << m_maxDeltaRho << " nu " << m_solver->m_nu << " kappa " << m_solver->m_kappa
        << " #Cells " << totalNumberOfCells << std::endl;
  }
  m_rhoLast = l_rho;
  m_ReLast = l_Re;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MInt currentId = m_bndCells[i].m_cellId;

    const MInt pCellId = currentId;

    MFloat l_t = 0.0;

    extrapolateVariable(ind, currentId, PV->T, &l_t);

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    m_solver->a_variable(pCellId, PV->T) = l_t;
    MFloat l_uu[nDim] = {0.0};
    for(MInt n = 0; n < nDim; n++) {
      l_uu[n] = m_solver->a_variable(pCellId, PV->VV[n]) / l_rho;
    }

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, squaredVelocity, l_uu);
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, squaredVelocity, l_uu);

    if(m_solver->m_isTransport) {
      MFloat l_c = 0.0;
      extrapolateVariable(ind, currentId, PV->C, &l_c);
      m_solver->setEqDistsTransport(pCellId, l_c, squaredVelocity, l_uu);
    }
  }
  writeBCOutput(index);
}

/** \brief Initialization of BC 4072 or 4082 after restart to determine
 * deltaRho, ReLast and rhoLast
 *
 * \author Mario Rüttgers
 * \date 15.02.2021
 *
 * \param[in] index the index of the segment of this BC
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40072_40082_init(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  MFloat mean_velocity = 0.0;
  MFloat mean_rho = 0.0;
  MFloat mean_oldRho = 0.0;
  MInt numberOfCellsPerDomain = 0;

  // first calculate rhoLast and oldRhoLast
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt pCellId = m_bndCells[i].m_cellId;

    if(m_solver->c_noChildren(pCellId) == 0 && !m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
      MFloat l_uu[nDim] = {0.0};
      MFloat l_rho = 0.0;
      MFloat l_oldRho = 0.0;

      for(MInt n = 0; n < nDim; n++) {
        l_uu[n] = m_solver->a_variable(pCellId, PV->VV[n]);
      }

      l_rho = m_solver->a_variable(pCellId, PV->RHO);
      l_oldRho = m_solver->a_oldVariable(pCellId, PV->RHO);

      mean_rho += l_rho;
      mean_oldRho += l_oldRho;

      for(MInt n = 0; n < nDim; n++) {
        mean_velocity += (m_bndNormals[ind][n] * l_uu[n]);
      }

      numberOfCellsPerDomain++;
    }
  }

  MFloat mean_rho_all = 0.0;
  MFloat mean_oldRho_all = 0.0;
  MFloat mean_velocity_all = 0.0;

  //  MInt totalNumberOfCells = m_totalNoBcCells[m_mapBndCndIdSegId[index]]; Including Halo cells
  MInt totalNumberOfCells = 0;
  if(m_noBCNeighbors[m_mapBndCndIdSegId[index]] > 1) {
    MPI_Allreduce(&numberOfCellsPerDomain, &totalNumberOfCells, 1, MPI_INT, MPI_SUM,
                  m_BCComm[m_mapBndCndIdSegId[index]], AT_, "numberOfCellsPerDomain",
                  "totalNumberOfCells"); // without Halo cells

    // Get the sum by communication in my group
    MPI_Allreduce(&mean_rho, &mean_rho_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_rho", "mean_rho_all");
    mean_rho_all /= totalNumberOfCells;

    // Get the sum by communication in my group
    MPI_Allreduce(&mean_oldRho, &mean_oldRho_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_oldRho", "mean_oldRho_all");
    mean_oldRho_all /= totalNumberOfCells;

    // Get the sum by communication in my group
    MPI_Allreduce(&mean_velocity, &mean_velocity_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                  "mean_velocity", "mean_velocity_all");
    mean_velocity_all /= totalNumberOfCells;
  } else {
    totalNumberOfCells = numberOfCellsPerDomain;
    mean_rho_all = mean_rho / totalNumberOfCells;
    mean_oldRho_all = mean_oldRho / totalNumberOfCells;
    mean_velocity_all = mean_velocity / totalNumberOfCells;
  }

  // Set m_rhoLast and m_deltaRho
  m_rhoLast = mean_rho_all;
  m_deltaRho = fabs(mean_rho_all - mean_oldRho_all);

  // Calculate the local Reynolds number:
  MFloat l_Re = mean_velocity_all / m_rhoLast * m_referenceLength / m_solver->m_nu;

  // Set m_ReLast
  m_ReLast = l_Re;
}

/** \LBBC{sec_LBBC_bc40100, bc40100, 4100}
 * \author Andreas Lintermann
 * \date 02.10.2012
 *
 * Outlet boundary condidtion similar to bc40000 including Thermal boundary treatment.
 *
 * \param[in] index index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40100(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40100(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> l_old_uu{};
    std::array<MFloat, nDim> u{};

    extrapolateVelocities(ind, pCellId, &l_uu[0]);

    MFloat old_rho = m_solver->a_oldVariable(pCellId, PV->RHO);
    for(MInt n = 0; n < nDim; n++) {
      l_old_uu[n] = m_solver->a_oldVariable(pCellId, PV->VV[n]);
    }

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    const MFloat old_squared_velocity = std::inner_product(&l_old_uu[0], &l_old_uu[nDim], &l_old_uu[0], .0);

    // non-reflecting bc (Finck, Haenel 2008)
    MFloat l_rho = 1.0;
    if(m_densityFluctuations)
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + (m_rho1 - 1.0)) / 2.0;
    else
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + m_rho1) / 2.0;

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    MFloat l_t = 1.0;
    m_solver->setEqDists(pCellId, l_rho, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, u.data());

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    m_solver->a_variable(pCellId, PV->T) = l_t;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
}

/** \LBBC{sec_LBBC_bc40110, bc40110, 4110}
 * \author Andreas Lintermann, Moritz Waldmann
 * \date 02.10.2012
 *
 * Outlet boundary condition similar to bc40110 including Thermal boundary treatment.
 *
 * This BC is similar to BC40100, except that the temperature is extrapolated from the inside.
 *
 * \param[in] index index of the BC<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40110(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40110(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];
  const MInt pvrho = PV->RHO;
  const MInt pvt = PV->T;

#ifdef WAR_NVHPC_PSTL
  const MInt globalTimeStep_ = globalTimeStep - 1;
  MInt begin = m_bndCndOffsets[index];
  MInt end = m_bndCndOffsets[index + 1];
  MInt offset = end - begin;

  maia::parallelFor<true>(0, offset, [=](MInt id) {
    MInt i = begin + id;
    if((globalTimeStep_) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#else
  maia::parallelFor(m_bndCndOffsets[index], m_bndCndOffsets[index + 1], [=](MInt i) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) return;
#endif

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> l_old_uu{};
    std::array<MFloat, nDim> u{};
    MFloat l_t = 1.0;

    extrapolateVelocities(ind, pCellId, &l_uu[0]);
    extrapolateVariable(ind, pCellId, PV->T, &l_t);

    MFloat old_rho = m_solver->a_oldVariable(pCellId, pvrho);
    for(MInt n = 0; n < nDim; n++) {
      l_old_uu[n] = m_solver->a_oldVariable(pCellId, n);
    }

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    const MFloat old_squared_velocity = std::inner_product(&l_old_uu[0], &l_old_uu[nDim], &l_old_uu[0], .0);

    // non-reflecting bc (Finck, Haenel 2008)
    MFloat l_rho = 1.0;
    if(m_densityFluctuations)
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + (m_rho1 - 1.0)) / 2.0;
    else
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + m_rho1) / 2.0;

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    m_solver->setEqDists(pCellId, l_rho, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, u.data());

    m_solver->a_variable(pCellId, pvrho) = l_rho;
    m_solver->a_variable(pCellId, pvt) = l_t;
#ifdef WAR_NVHPC_PSTL
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, n) = l_uu[n];
    }
  });
#else
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  });
#endif
  writeBCOutput(index);
}

/** \brief Outlet boundary condition similar to bc40110 including Transport boundary treatment.
 * \author Shota Ito
 * \date 08.06.2022
 *
 * This BC is similar to BC40100, except that the concentration is extrapolated from the inside.
 *
 * \param[in] index index of the BC
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40120(MInt index) {
  TRACE();

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;
    const MInt pCellId = currentId;

    std::array<MFloat, nDim> l_uu{};
    std::array<MFloat, nDim> l_old_uu{};
    std::array<MFloat, nDim> u{};
    MFloat l_c = 1.0;
    MFloat l_t = 1.0;

    extrapolateVelocities(ind, pCellId, &l_uu[0]);
    if(m_solver->m_isThermal) {
      extrapolateVariable(ind, pCellId, PV->T, &l_t);
    }
    extrapolateVariable(ind, pCellId, PV->C, &l_c);

    MFloat old_rho = m_solver->a_oldVariable(pCellId, PV->RHO);
    for(MInt n = 0; n < nDim; n++) {
      l_old_uu[n] = m_solver->a_oldVariable(pCellId, PV->VV[n]);
    }

    const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    const MFloat old_squared_velocity = std::inner_product(&l_old_uu[0], &l_old_uu[nDim], &l_old_uu[0], .0);

    // non-reflecting bc (Finck, Haenel 2008)
    MFloat l_rho = 1.0;
    if(m_densityFluctuations)
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + (m_rho1 - 1.0)) / 2.0;
    else
      l_rho = (old_rho + F1BCS * (sqrt(squaredVelocity) - sqrt(old_squared_velocity)) + m_rho1) / 2.0;

    // TODO: As long as rho * u is saved as a primitive variable in the thermal collision steps,
    // rho * u should be set in the BCs as well. The Eq dists, however, require only u.
    for(MInt n = 0; n < nDim; n++) {
      u[n] = l_uu[n] / l_rho;
    }
    // calulate the equilibrium distribution functions for the new density
    m_solver->setEqDists(pCellId, l_rho, u.data());
    m_solver->setEqDistsThermal(pCellId, l_t, l_rho, u.data());
    m_solver->setEqDistsTransport(pCellId, l_c, u.data());

    m_solver->a_variable(pCellId, PV->RHO) = l_rho;
    if(m_solver->m_isThermal) m_solver->a_variable(pCellId, PV->T) = l_t;
    m_solver->a_variable(pCellId, PV->C) = l_c;
    for(MInt n = 0; n < nDim; n++) {
      m_solver->a_variable(pCellId, PV->VV[n]) = l_uu[n];
    }
  }
}

/** \LBBC{sec_LBBC_bc10010, bc10010, 1010}
 * Lattice Boltzmann inflow boundary condition
 *
 * extrapolation Chen, Martinez 1996, prescribed velocity
 * for arbitrary directions (lhs)<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10010(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10010(MInt index) {
  TRACE();

  MFloat rho = F0;
  std::array<MFloat, nDim> u{};
  for(MInt d = 0; d < nDim; d++) {
    u[d] = std::numeric_limits<MFloat>::max();
  }

  MFloat tmp;
  MFloat b[2 * nDim] = {0.0};
  MInt tmpDistId, currentId;

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    currentId = m_bndCells[i].m_cellId;

    m_omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)));
    // --------------------------------------------------------------------
    // perform collision step once more for boundary cell (with new eq-distributions)

    rho = m_solver->a_variable(currentId, PV->RHO);

    for(MInt d = 0; d < nDim; d++) {
      u[d] = m_initialVelocityVecs[ind][d] * m_Ma * LBCS * m_bndCells[i].m_multiplier;
      b[2 * d] = -u[d];
      b[2 * d + 1] = u[d];
    }

    const MFloat tmp2 = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // Calculation of distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      m_solver->a_distribution(currentId, j) =
          m_solver->a_oldDistribution(currentId, j)
          + m_omega
                * (Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2)
                   - m_solver->a_oldDistribution(currentId, j));
    }

    // Calculation of distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      m_solver->a_distribution(currentId, tmpDistId + j) =
          m_solver->a_oldDistribution(currentId, tmpDistId + j)
          + m_omega
                * (Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
                   - m_solver->a_oldDistribution(currentId, tmpDistId + j));
    }

    // Calculation of distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      m_solver->a_distribution(currentId, tmpDistId + j) =
          m_solver->a_oldDistribution(currentId, tmpDistId + j)
          + m_omega
                * (Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
                   - m_solver->a_oldDistribution(currentId, tmpDistId + j));
    }
    // Calculation of distribution for rest particle distribution (center)
    m_solver->a_distribution(currentId, Ld::lastId()) =
        m_solver->a_oldDistribution(currentId, Ld::lastId())
        + m_omega * (Ld::tp(0) * (rho - F1B2 * F1BCSsq * tmp2) - m_solver->a_oldDistribution(currentId, Ld::lastId()));


    // --------------------------------------------------------------------
    // extrapolation of incoming distributions

    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(j)) == 0) {
        m_solver->a_oldDistribution(currentId, j) = m_solver->a_distribution(currentId, j);
      }
    }

    m_solver->a_variable(currentId, PV->RHO) = rho;
    for(MInt d = 0; d < nDim; d++) {
      m_solver->a_variable(currentId, d) = u[d];
    }
  }
}

/** \LBBC{sec_LBBC_bc10020, bc10020, 1020}
 * non-eq, prescribed velocity<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10020(MInt index)
 */

template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10020(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm
  MFloat rho;
  std::array<MFloat, nDim> u{};
  MFloat tmp, tmp2;
  std::array<MFloat, 2 * nDim> b{};
  MInt tmpDistId, currentId;

  MFloatScratchSpace eDistributions(nDist, AT_, "eDistributions");
  MFloatScratchSpace nePart(nDist, AT_, "nePart");

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    currentId = m_bndCells[i].m_cellId;

    m_omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)));

    // 1. Calculate eq-Distributions from actual macroscopic values
    //--------------------------------------------

    // Calculate macroscopic values


    //     for( MInt j = 0; j < nDist - 1; j ++){
    // 	if (m_solver->a_hasNeighbor(currentId, j) > 0 && !m_solver->a_onlyBoundary(m_solver->c_neighborId(currentId,
    // j))){ 	    if (m_solver->a_onlyBoundary(m_solver->c_neighborId(currentId, j)) == true) cerr << "extrapolation
    // from boundary cell!" << endl; 	    neighborId = m_solver->c_neighborId(currentId, j); 	    break;
    // 	}
    //     }

    //     // extrapolate inner density values
    //     if (m_solver->a_hasNeighbor(currentId, 2) == 0) {
    // 	neighborId = m_solver->c_neighborId(currentId, 3);
    //     }
    //     else {
    // 	neighborId = m_solver->c_neighborId(currentId, 2);
    //     }

    //     rho = m_solver->a_variable(neighborId, PV->RHO);
    //     u[0] = m_solver->a_variable(neighborId, PV->U);
    //     u[1] = m_solver->a_variable(neighborId, PV->V);
    //     u[2] = m_solver->a_variable(neighborId, PV->W);

    // macroscopic values
    rho = m_solver->a_variable(currentId, PV->RHO);
    for(MInt n = 0; n < nDim; n++) {
      u[n] = m_solver->a_variable(currentId, PV->VV[n]);
      b[2 * n] = -u[n];
      b[2 * n + 1] = u[n];
    }

    tmp2 = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // Calculation of eq-distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      eDistributions[j] = Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }

    // 2. Calculation of non-eq-parts3
    //--------------------------------------------
    for(MInt j = 0; j < nDist; j++) {
      nePart[j] = m_solver->a_oldDistribution(currentId, j) - eDistributions[j];
    }

    // 3. Calculation of incoming distributions in bnd cell
    //--------------------------------------------

    //     if (m_lbControlInflow == 0) {
    //       for(MInt d=0; d < nDim; d++){
    // 	u[d] = m_initialVelocityVecs[ind][d] * m_Ma * LBCS ;
    //       }
    //     }
    //     else {
    //       for(MInt d=0; d < nDim; d++){
    // 	u[d] = m_initialVelocityVecs[ind][d] * m_Ma * LBCS * m_bndCells[i].m_multiplier ;
    //       }
    //     }

    // set velocity to zero
    for(MInt d = 0; d < nDim; d++) {
      u[d] = m_initialVelocityVecs[ind][d] * m_Ma * LBCS * m_bndCells[i].m_multiplier;
      b[2 * d] = -u[d];
      b[2 * d + 1] = u[d];
    }

    tmp2 = std::inner_product(&u[0], &u[nDim], &u[0], .0);

    // Calculation of distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) =
          Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[j]);
      // m_solver->a_distribution(m_bndCells[i].m_cellId, j) = Ld::tp(1) * F1BCSsq *  ( rho * CSsq + b[j] +
      // b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2) + ((1-m_omega) * nePart[j]);
    }

    // Calculation of distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, tmpDistId + j) =
          Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[tmpDistId + j]);
      // m_solver->a_distribution(m_bndCells[i].m_cellId, tmpDistId + j) = Ld::tp(2) * F1BCSsq * ( rho * CSsq +
      // tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2) + ((1-m_omega) * nePart[tmpDistId + j]);
    }

    // Calculation of distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, tmpDistId + j) =
          Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[tmpDistId + j]);
      // m_solver->a_distribution(m_bndCells[i].m_cellId, tmpDistId + j) = Ld::tp(3) * F1BCSsq * ( rho * CSsq +
      // tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2) + ((1-m_omega) * nePart[tmpDistId + j]);
    }
  }
}

/** \LBBC{sec_LBBC_bc40130, bc40130, 4130}
 *
 * Guo, Extrapolation of velocity and non-eq parts
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40130(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40130(MInt index) {
  TRACE();

  // For now testing only the D3Q19 algorithm
  MFloat rho, u[nDim], old_u[nDim];
  MFloat tmp = F0;
  MFloat tmp2 = F0;
  MFloat old_tmp2 = F0;
  MFloat b[2 * nDim];
  MInt tmpDistId, currentId;
  MFloatScratchSpace eDistributions(nDist, AT_, "eDistributions");
  MFloatScratchSpace nePart(nDist, AT_, "nePart");

  MInt ind = m_mapSegIdsInOutCnd[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    currentId = m_bndCells[i].m_cellId;

    // leave out halo cells
    //  if (currentId >= m_solver->noInternalCells())
    //	continue;
    if(m_solver->a_isHalo(currentId)) continue;

    MFloat l_t = 1.0;

    extrapolateVariable(ind, currentId, PV->T, &l_t);

    MFloat nu =
        m_Ma * LBCS / m_solver->m_Re * m_referenceLength * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId));
    m_omega = 2.0 / (1.0 + 6.0 * nu);
    // 1. Calculate eq-Distributions from actual macroscopic values
    //--------------------------------------------

    // macroscopic values
    rho = m_solver->a_variable(currentId, PV->RHO);
    u[0] = m_solver->a_variable(currentId, PV->U);
    u[1] = m_solver->a_variable(currentId, PV->V);
    IF_CONSTEXPR(nDim == 3) u[2] = m_solver->a_variable(currentId, PV->W);

    old_u[0] = m_solver->a_oldVariable(currentId, PV->U);
    old_u[1] = m_solver->a_oldVariable(currentId, PV->V);
    IF_CONSTEXPR(nDim == 3) old_u[2] = m_solver->a_oldVariable(currentId, PV->W);

    tmp2 = F0;
    old_tmp2 = F0;
    for(MInt d = 0; d < nDim; d++) {
      tmp2 += (u[d] * u[d]);
      old_tmp2 += (old_u[d] * old_u[d]);
      b[2 * d] = -u[d];
      b[2 * d + 1] = u[d];
    }

    // Calculation of eq-distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      eDistributions[j] = Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }

    // 2. Calculation of non-eq-parts3
    //--------------------------------------------
    for(MInt j = 0; j < nDist - 1; j++) {
      nePart[j] = m_solver->a_distribution(currentId, j);
    }
    for(MInt j = 0; j < nDist - 1; j++) {
      nePart[j] -= eDistributions[j];
    }
    // had to be split leads to compiler error with optimization

    // WARNING!!! The density should be calculated as in 3d
    // However, result is correct
    IF_CONSTEXPR(nDim == 3) {
      if(m_densityFluctuations) {
        rho = (rho + F1BCS * (sqrt(tmp2) - sqrt(old_tmp2))) / 2.0;
      } else {
        // rho = 1.0;
        // rho = ( m_solver->a_variable(m_bndCells[i].m_cellId, PV->RHO) + rho ) / 2.0 ;
        // WARNING!!! The density should be calculated as in 3d
        // However, result is correct
        rho = (rho + F1BCS * (sqrt(tmp2) - sqrt(old_tmp2)) + 1.0) / 2.0;
      }
    }
    else {
      rho = (rho + 1.0) / 2.0;
    }

    // 3. Calculation of incoming distributions in bnd cell
    //--------------------------------------------

    // Calculation of distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      // if (!m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(j))){
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) =
          Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[j]);
      //}
    }

    // Calculation of distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      // if (!m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(tmpDistId + j))){
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, tmpDistId + j) =
          Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[tmpDistId + j]);
      //}
    }

    // Calculation of distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      // if (!m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(tmpDistId + j))){
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      m_solver->a_oldDistribution(m_bndCells[i].m_cellId, tmpDistId + j) =
          Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
          + ((1 - m_omega) * nePart[tmpDistId + j]);
      //      }
    }

    m_solver->setEqDistsThermal(currentId, l_t, rho, u);
  }
  writeBCOutput(index);
}


/** \LBBC{sec_LBBC_outflow, outflow, 301x}
 * Lattice Boltzmann inflow boundary condition
 *
 * Simple condition which sets the incoming distributions
 * to the last inner value (rhs) in a given axis direction
 *
 * \tparam    direction Boundary orientation in Cartesian direction
 * \param[in] index     Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::outflow(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt direction>
void LbBndCndDxQy<nDim, nDist, SysEqn>::outflow(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(m_bndCells[i].m_cellId, Ld::oppositeDist(j)) == 0) {
        const MInt neighbor1 = m_solver->c_neighborId(m_bndCells[i].m_cellId, Ld::oppositeDist(direction));
        m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) = m_solver->a_oldDistribution(neighbor1, j);
      }
    }
  }
}


/** \LBBC{sec_LBBC_slipFlow, slipFlow, 3020}
 * \author Miro Gondrum
 * \date   01.09.2020
 *
 * Lattice Boltzmann slip wall boundary condition for Cartesian planes
 *
 * Slip wall condition for Cartesian walls
 *
 * \note If this BC and a no-slip wall are intersecting,
 * multiBCTreatment='multiple' should be used and no-slip BC must be performed
 * after this BC.
 *
 * TODO labels:LB,totest THERMAL: Untested for thermal.
 *
 * \param[in] index Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::slipFlow(MInt index)
 **/
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool thermal>
void LbBndCndDxQy<nDim, nDist, SysEqn>::slipFlow(MInt index) {
  TRACE();
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndCells[i].m_cellId;
    const MInt maxLvlDiff = m_solver->maxLevel() - m_solver->a_level(cellId);
    //--check: is BC required for this cell now?--------------------------------
    if(globalTimeStep % IPOW2(maxLvlDiff) != 0 || m_solver->a_isHalo(cellId) || !m_solver->a_isActive(cellId)) {
      continue;
    }
    for(MInt trgDist = 0; trgDist < nDist - 1; trgDist++) {
      //--check: is BC required for this distribution?--------------------------
      const MInt srcDir1 = Ld::oppositeDist(trgDist); // init with original propagation src direction
      if(!(m_bndCells[i].m_distances[srcDir1] <= 0.5)) {
        continue; // No cut in opposite dir, hence continue with next dist.
      }
      //--get srcDir1 vector with values {0,1,2}---------------------------------
      std::array<MInt, nDim> srcDirVec, srcDistVec;
      for(MInt d = 0; d < nDim; d++) {
        srcDirVec[d] = Ld::idFld(srcDir1, d);  // opposite trgDistVec
        srcDistVec[d] = Ld::idFld(trgDist, d); // init with value of trgDistVec
      }
      //--get srcDistVec--------------------------------------------------------
      // check for neighbors in those Cartesian directions the srcDirVec is composed of
      MBool normalPropagation = true;
      for(MInt d = 0; d < nDim; d++) {
        if(srcDirVec[d] == 1) continue;
        // else: mirror if no neighbor exist in this direction
        const MInt tmpDir = (srcDirVec[d] / 2) + (d * 2); // convert {0,2} into id of Cartesian dir. That is
                                                          // {0,1}, {2,3}, or {4,5} for x,y, or z respectively
        if(m_bndCells[i].m_distances[tmpDir] <= 0.5) {
          srcDistVec[d] = 2 - srcDistVec[d]; // mirror
          normalPropagation = false;
        }
      }
      if(normalPropagation) continue;
      //--get srcDirVec---------------------------------------------------------
      // srcDirVector is in the opposite direction of the vector which is
      // located between srcDistVec and trgDistVec
      for(MInt d = 0; d < nDim; d++) {
        srcDirVec[d] = 2 - (Ld::idFld(trgDist, d) + srcDistVec[d]) / 2;
      }
      //--convert from vec to distribution id-----------------------------------
      const MInt srcDist = Ld::dirFld(srcDistVec[0], srcDistVec[1], srcDistVec[2]);
      const MInt srcDir2 = Ld::dirFld(srcDirVec[0], srcDirVec[1], srcDirVec[2]);
      //--set value-------------------------------------------------------------
      // I) determine cell state
      enum CellState { NORMAL, AT_INTERFACE, WALL, UNDEFINED };
      CellState state;
      MInt srcCell;
      if(srcDir2 == Ld::lastId()) {
        srcCell = cellId;
        state = NORMAL;
      } else if(m_solver->a_hasNeighbor(cellId, srcDir2)) {
        srcCell = m_solver->c_neighborId(cellId, srcDir2);
        if(m_solver->a_isActive(srcCell)) {
          state = NORMAL;
        } else {
          if(m_solver->m_isRefined && m_solver->a_isInterfaceParent(cellId)) {
            state = AT_INTERFACE;
          } else {
            state = WALL;
          }
        }
      } else if(m_solver->m_isRefined && m_solver->a_isInterfaceChild(cellId)) {
        const MInt parentId = m_solver->c_parentId(cellId);
        if(m_solver->a_hasNeighbor(parentId, srcDir2)) {
          state = AT_INTERFACE;
        } else {
          state = UNDEFINED;
        }
      } else {
        // In case of missing neighbor (even for interface parent cell) no
        // behavior can be defined. This must come from the BC defined from
        // the missing direction.
        state = UNDEFINED;
      }
      // TODO labels:LB miro: Until here, everthing should be precomputed and stored in
      //            advance (state, srcDist, trgDist, srcCell).

      // II) perfrom calculation depending on cell state
      switch(state) {
        default:
        case NORMAL: {
          m_solver->a_oldDistribution(cellId, trgDist) = m_solver->a_distribution(srcCell, srcDist);
          if constexpr(thermal) {
            m_solver->a_oldDistributionThermal(cellId, trgDist) = m_solver->a_distributionThermal(srcCell, srcDist);
          }
          break;
        }
        case AT_INTERFACE: {
          m_solver->a_oldDistribution(cellId, trgDist) = m_solver->a_oldDistribution(cellId, srcDist);
          if constexpr(thermal) {
            m_solver->a_oldDistributionThermal(cellId, trgDist) = m_solver->a_oldDistributionThermal(cellId, srcDist);
          }
          break;
        }
        case WALL: {
          // interpolated bounce back: here only for solid nghbr cell
          // TODO labels:LB problem if nghbr belongs to more than 1 boundary. Indeed,
          // this is the case here. So wall type has to be set as last BC in
          // input file -.-
          const MInt solidNghbrBndId = m_solver->a_bndId(srcCell);
          if(m_bndCells[solidNghbrBndId].m_isFluid == false && m_bndCells[solidNghbrBndId].m_distances[srcDist] < 0.5) {
            const MFloat F2q = F2 * (F1 - m_bndCells[solidNghbrBndId].m_distances[srcDist]);
            const MInt oppositeDir = Ld::oppositeDist(trgDist);
            m_solver->a_oldDistribution(cellId, trgDist) =
                (m_solver->a_distribution(cellId, oppositeDir) / F2q)
                + (m_solver->a_distribution(cellId, trgDist) * (F2q - F1) / F2q);
            if constexpr(thermal) {
              m_solver->a_oldDistributionThermal(cellId, trgDist) =
                  (m_solver->a_distributionThermal(cellId, oppositeDir) / F2q)
                  + (m_solver->a_distributionThermal(cellId, trgDist) * (F2q - F1) / F2q);
            }
          } else {
            // SOMETHING has to be went wrong
          }
          break;
        }
        case UNDEFINED: {
          break;
        }
      }
    }
  }
}

/** \LBBC{sec_LBBC_outflowLinear, outflowLinear, 303x}
 * linear extrapolation
 *
 * \tparam    direction Boundary orientation in Cartesian direction
 * \param[in] index     Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::outflowLinear(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt direction>
void LbBndCndDxQy<nDim, nDist, SysEqn>::outflowLinear(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;
    for(MInt j = 0; j < nDist - 1; j++) {
      if(m_solver->a_hasNeighbor(m_bndCells[i].m_cellId, Ld::oppositeDist(j)) == 0) {
        const MInt neighbor1 = m_solver->c_neighborId(m_bndCells[i].m_cellId, Ld::oppositeDist(direction));
        const MInt neighbor2 = m_solver->c_neighborId(neighbor1, Ld::oppositeDist(direction));
        m_solver->a_oldDistribution(m_bndCells[i].m_cellId, j) =
            2.0 * m_solver->a_oldDistribution(neighbor1, j) - m_solver->a_oldDistribution(neighbor2, j);
      }
    }
  }
}

/** \LBBC{sec_LBBC_pab, pab, 4060}
 * Lattice Boltzmann outflow boundary condition
 *
 * Pressure anti bounce back method (I. Ginzburg et al. 2008)
 * see also: Izquierdo et al. 2008
 *
 * constant pressure is prescribed for a boundary located at q=1/2
 * velocity is extrapolated
 *
 * \param[in] index     Boundary index<br>
 *
 * LbBndCndDxQy<nDim, nDist, SysEqn>::pab(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::pab(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    const MInt currentId = m_bndCells[i].m_cellId;

    const MFloat rho = m_solver->a_variable(currentId, PV->RHO);
#ifdef WAR_NVHPC_PSTL
    MFloat vel[nDim] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      vel[d] = m_solver->a_variable(currentId, d);
    }
    const MFloat* const u = vel;
#else
    const MFloat* const u = &m_solver->a_variable(currentId, PV->U);
#endif

    // Calculate symmetric part of eq-distributions at the border
    //---------------------------------------------------
    const MFloat tmp2 = std::inner_product(&u[0], &u[nDim], &u[0], 0.0);
    MFloat b[2 * nDim];
    for(MInt n = 0; n < nDim; n++) {
      b[2 * n] = -u[n];
      b[2 * n + 1] = u[n];
    }

    // Calculation of eq. distributions for directions with only one component
    MFloat eqSym[nDist - 1];
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      eqSym[j] = Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2);
    }

    // Calculation of eq. distributions for directions with two components
    MInt tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      const MInt tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      eqSym[tmpDistId + j] = Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }

    // Calculation of eq. distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      const MInt tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      eqSym[tmpDistId + j] = Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }

    for(MInt j = 0; j < nDist - 1; j++) {
      if(!m_solver->a_hasNeighbor(currentId, j))
        m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j)) =
            -m_solver->a_distribution(currentId, j) + 2.0 * eqSym[j];
    }
  }
}

/** \brief calculates the averaged flow values at a certain outlet/inlet
 *
 * Attention: This function is assuming the compressible LBM formulation, i.e. a_variable storing rho*u
 *
 * \author Moritz Waldmann
 * \date 12.03.2019
 *
 * \param[in] index Boundary index (as specified in the geometry file) for wich is output is written
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::writeBCOutput(MInt index) {
  TRACE();

  if(!this->m_calcBcResidual) return;
  MInt ind = m_mapSegIdsInOutCnd[index];

  MFloat mean_velocity = 0.0;
  MFloat mean_rho = 0.0;
  MFloat mean_u[nDim] = {0.0};
  MInt numberOfCellsPerDomain = 0;
  MFloat mean_p_d = 0.0;
  MFloat mean_p_stat = 0.0;
  MFloat mean_t = 0.0;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(m_bndCells[i].m_cellId)) != 0) continue;

    MInt pCellId = m_bndCells[i].m_cellId;

    if(m_solver->c_noChildren(pCellId) == 0 && !m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
      MFloat l_rho = m_solver->a_variable(pCellId, PV->RHO);
      MFloat l_uu[nDim] = {0.0};
      for(MInt n = 0; n < nDim; n++) {
        l_uu[n] = m_solver->a_variable(pCellId, PV->VV[n]) / l_rho;
      }

      MFloat l_t = 1.0;
      if(m_solver->m_isThermal) {
        l_t = m_solver->a_variable(pCellId, PV->T);
      }

      for(MInt n = 0; n < nDim; n++) {
        mean_velocity += (m_bndNormals[ind][n] * l_uu[n]);
      }
      const MFloat squaredVelocity = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);
      mean_rho += l_rho;
      numberOfCellsPerDomain++;

      for(MInt n = 0; n < nDim; n++) {
        mean_u[n] += l_uu[n];
      }

      mean_p_d += F1B2 * l_rho * squaredVelocity;
      mean_p_stat += F1B3 * l_rho;
      mean_t += l_t;
    }
  }

  MFloat mean_velocity_all = 0.0;
  MFloat mean_rho_all = 0.0;
  MFloat mean_u_all[nDim] = {0.0};
  MFloat mean_p_d_all = 0.0;
  MFloat mean_p_stat_all = 0.0;
  MFloat mean_t_all = 0.0;
  // MInt totalNumberOfCells = m_totalNoBcCells[m_mapBndCndIdSegId[index]]; Including Halo cells
  MInt totalNumberOfCells = 0;
  switch(m_noBCNeighbors[m_mapBndCndIdSegId[index]]) {
    case 0: {
      break;
    }

    case 1: {
      totalNumberOfCells = numberOfCellsPerDomain;
      mean_velocity_all = mean_velocity / totalNumberOfCells;
      mean_rho_all = mean_rho / totalNumberOfCells;

      for(MInt n = 0; n < nDim; n++) {
        mean_u_all[n] = mean_u[n] / totalNumberOfCells;
      }

      mean_p_d_all = mean_p_d / totalNumberOfCells;
      mean_p_stat_all = mean_p_stat / totalNumberOfCells;
      if(m_solver->m_isThermal) {
        mean_t_all = mean_t / totalNumberOfCells;
      } else {
        mean_t_all = 1.0;
      }
      MFloat ReynoldsNr = m_referenceLength * mean_velocity_all / m_solver->m_nu;

      if(m_solver->domainId() == m_BCneighbors[m_mapBndCndIdSegId[index]][0]) {
        m_BCResidualStream[m_mapBndCndIdSegId[index]].open(m_BCOutputFileName[m_mapBndCndIdSegId[index]],
                                                           std::ios_base::app);
        m_BCResidualStream[m_mapBndCndIdSegId[index]] << globalTimeStep << " p_dyn " << mean_p_d_all << " p_stat "
                                                      << mean_p_stat_all << " Mag " << mean_velocity_all << " v_x "
                                                      << mean_u_all[0] << " v_y " << mean_u_all[1];
        IF_CONSTEXPR(nDim == 3) m_BCResidualStream[m_mapBndCndIdSegId[index]] << " v_z " << mean_u_all[2];

        m_BCResidualStream[m_mapBndCndIdSegId[index]] << " rho " << mean_rho_all << " nu " << m_solver->m_nu << " Re "
                                                      << ReynoldsNr;

        if(m_solver->m_isThermal) {
          m_BCResidualStream[m_mapBndCndIdSegId[index]] << " t " << mean_t_all << " kappa " << m_solver->m_kappa;
        }
        m_BCResidualStream[m_mapBndCndIdSegId[index]] << " #Cells " << totalNumberOfCells << std::endl;
        m_BCResidualStream[m_mapBndCndIdSegId[index]].close();
      }
      break;
    }

    default: {
      MPI_Allreduce(&numberOfCellsPerDomain, &totalNumberOfCells, 1, MPI_INT, MPI_SUM,
                    m_BCComm[m_mapBndCndIdSegId[index]], AT_, "numberOfCellsPerDomain",
                    "totalNumberOfCells"); // without Halo cells
      // Get the sum by communication in my group
      MPI_Allreduce(&mean_velocity, &mean_velocity_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]],
                    AT_, "mean_velocity", "mean_velocity_all");
      MPI_Allreduce(&mean_rho, &mean_rho_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                    "mean_rho", "mean_rho_all");
      MPI_Allreduce(mean_u, mean_u_all, nDim, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_, "mean_u",
                    "mean_u_all");
      MPI_Allreduce(&mean_p_d, &mean_p_d_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                    "mean_p_d", "mean_p_d_all");
      MPI_Allreduce(&mean_p_stat, &mean_p_stat_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_,
                    "mean_p_stat", "mean_p_stat_all");

      if(m_solver->m_isThermal) {
        MPI_Allreduce(&mean_t, &mean_t_all, 1, MPI_DOUBLE, MPI_SUM, m_BCComm[m_mapBndCndIdSegId[index]], AT_, "mean_t",
                      "mean_t_all");
        mean_t_all /= totalNumberOfCells;
      } else {
        mean_t_all = 1.0;
      }

      mean_velocity_all /= totalNumberOfCells;
      mean_rho_all /= totalNumberOfCells;

      for(MInt n = 0; n < nDim; n++) {
        mean_u_all[n] /= totalNumberOfCells;
      }

      mean_p_d_all /= totalNumberOfCells;
      mean_p_stat_all /= totalNumberOfCells;

      MFloat ReynoldsNr = m_referenceLength * mean_velocity_all / m_solver->m_nu;

      if(m_solver->domainId() == m_BCneighbors[m_mapBndCndIdSegId[index]][0]) {
        m_BCResidualStream[m_mapBndCndIdSegId[index]].open(m_BCOutputFileName[m_mapBndCndIdSegId[index]],
                                                           std::ios_base::app);
        m_BCResidualStream[m_mapBndCndIdSegId[index]]
            << globalTimeStep << " p_dyn " << mean_p_d_all
            << " p_stat "
            //                                                      << mean_p_stat_all << " v_n,BC " <<
            //                                                      mean_velocity_all << " v_x "
            << mean_p_stat_all << " Mag " << mean_velocity_all << " v_x " << mean_u_all[0] << " v_y " << mean_u_all[1];
        IF_CONSTEXPR(nDim == 3) m_BCResidualStream[m_mapBndCndIdSegId[index]] << " v_z " << mean_u_all[2];

        m_BCResidualStream[m_mapBndCndIdSegId[index]] << " rho " << mean_rho_all << " nu " << m_solver->m_nu << " Re "
                                                      << ReynoldsNr;
        if(m_solver->m_isThermal) {
          m_BCResidualStream[m_mapBndCndIdSegId[index]] << " t " << mean_t_all << " kappa " << m_solver->m_kappa;
        }
        m_BCResidualStream[m_mapBndCndIdSegId[index]] << " #Cells " << totalNumberOfCells << std::endl;
        m_BCResidualStream[m_mapBndCndIdSegId[index]].close();
      }
    }
  }
}

/**
 * \brief calls refill-function corresponding to the refillMethodOrder, defined by property "refillMethodOrder"
 *
 * \author Johannes Grafen <johannes.grafen@rwth-aachen.de
 *
 * \param[in] pCellId Cell id of the newly emerged cell
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::refillEmergedCell(const MInt pCellId) {
  TRACE();
  if(m_refillMethodOrder == 1) {
    (this->LbBndCndDxQy<nDim, nDist, SysEqn>::refillEmergedCellNormalExtrapolationLinear)(pCellId);
  } else {
    (this->LbBndCndDxQy<nDim, nDist, SysEqn>::refillEmergedCellNormalExtrapolationQuadratic)(pCellId);
  }
}

/**
 * \brief Initialize emerged cell by performing linear extrapolation in boundary normal direction
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] pCellId Cell id of the newly emerged cell
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::refillEmergedCellNormalExtrapolationLinear(const MInt pCellId) {
  TRACE();

  ASSERT(m_solver->noHaloLayers() > 1,
         "linear extrapolation for refillEmergedCell, at least 2 halo layers are required!");

  const MInt mbCell = m_boundaryCellMappingMb[pCellId];
  ASSERT(mbCell > -1, "Cell " << pCellId << " is no G0 Cell, thus cant be refilled!");

  // Find direction with shortest path to boundary
  MInt minDistDir = -1;
  MFloat minWallDistance = F2;
  for(MInt dist = 0; dist < nDist - 1; dist++) {
    if(minWallDistance > m_boundaryCellsMb.distance(mbCell, dist)) {
      minWallDistance = m_boundaryCellsMb.distance(mbCell, dist);
      minDistDir = Ld::oppositeDist(dist);
    }
  }

  MFloat l_rho = F0;
  std::array<MFloat, nDim> l_vv{};

  if(m_solver->a_hasNeighbor(pCellId, minDistDir)
     && m_solver->a_isActive(m_solver->c_neighborId(pCellId, minDistDir))) {
    const MInt nghbor = m_solver->c_neighborId(pCellId, minDistDir);
    if(m_solver->a_hasNeighbor(nghbor, minDistDir)
       && m_solver->a_isActive(m_solver->c_neighborId(nghbor, minDistDir))) {
      const MInt nnghbor = m_solver->c_neighborId(nghbor, minDistDir);

      // Two neighbors found
      for(MInt j = 0; j < nDist; j++) {
        m_solver->a_oldDistribution(pCellId, j) =
            F2 * m_solver->a_oldDistribution(nghbor, j) - F1 * m_solver->a_oldDistribution(nnghbor, j);

        m_solver->a_distribution(pCellId, j) =
            F2 * m_solver->a_distribution(nghbor, j) - F1 * m_solver->a_distribution(nnghbor, j);
      }
    } else {
      // One neighbor found
      for(MInt j = 0; j < nDist; j++) {
        m_solver->a_oldDistribution(pCellId, j) = m_solver->a_oldDistribution(nghbor, j);
        m_solver->a_distribution(pCellId, j) = m_solver->a_distribution(nghbor, j);
      }
    }
  }

  m_solver->calculateMacroscopicVariables(pCellId, l_rho, l_vv.data());

  m_solver->a_variable(pCellId, PV->RHO) = l_rho;
  for(MInt n = 0; n < nDim; n++) {
    m_solver->a_variable(pCellId, PV->VV[n]) = l_vv[n];
  }
}

/**
 * \brief Initialize emerged cell by performing quadratic extrapolation in boundary normal direction
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in] pCellId Cell id of the newly emerged cell
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::refillEmergedCellNormalExtrapolationQuadratic(const MInt pCellId) {
  TRACE();

  ASSERT(m_solver->noHaloLayers() > 2,
         "quadratic extrapolation for refillEmergedCell, at least 3 halo layers are required!");

  const MInt mbCell = m_boundaryCellMappingMb[pCellId];
  ASSERT(mbCell > -1, "Cell " << pCellId << " is no G0 Cell, thus cant be refilled!");

  // Find direction with shortest path to boundary
  MInt minDistDir = -1;
  MFloat minWallDistance = F2;
  for(MInt dist = 0; dist < nDist - 1; dist++) {
    if(minWallDistance > m_boundaryCellsMb.distance(mbCell, dist)) {
      minWallDistance = m_boundaryCellsMb.distance(mbCell, dist);
      minDistDir = Ld::oppositeDist(dist);
    }
  }

  MFloat l_rho = F0;
  std::array<MFloat, nDim> l_vv{};

  if(m_solver->a_hasNeighbor(pCellId, minDistDir)
     && m_solver->a_isActive(m_solver->c_neighborId(pCellId, minDistDir))) {
    const MInt nghbor = m_solver->c_neighborId(pCellId, minDistDir);
    if(m_solver->a_hasNeighbor(nghbor, minDistDir)
       && m_solver->a_isActive(m_solver->c_neighborId(nghbor, minDistDir))) {
      const MInt nnghbor = m_solver->c_neighborId(nghbor, minDistDir);
      if(m_solver->a_hasNeighbor(nnghbor, minDistDir)
         && m_solver->a_isActive(m_solver->c_neighborId(nnghbor, minDistDir))) {
        const MInt nnnghbor = m_solver->c_neighborId(nnghbor, minDistDir);

        // Three neighbors found
        for(MInt j = 0; j < nDist; j++) {
          m_solver->a_oldDistribution(pCellId, j) = F3 * m_solver->a_oldDistribution(nghbor, j)
                                                    - F3 * m_solver->a_oldDistribution(nnghbor, j)
                                                    + F1 * m_solver->a_oldDistribution(nnnghbor, j);

          m_solver->a_distribution(pCellId, j) = F3 * m_solver->a_distribution(nghbor, j)
                                                 - F3 * m_solver->a_distribution(nnghbor, j)
                                                 + F1 * m_solver->a_distribution(nnnghbor, j);
        }

      } else {
        // Two neighbors found
        for(MInt j = 0; j < nDist; j++) {
          m_solver->a_oldDistribution(pCellId, j) =
              F2 * m_solver->a_oldDistribution(nghbor, j) - F1 * m_solver->a_oldDistribution(nnghbor, j);

          m_solver->a_distribution(pCellId, j) =
              F2 * m_solver->a_distribution(nghbor, j) - F1 * m_solver->a_distribution(nnghbor, j);
        }
      }

    } else {
      // One neighbor found
      for(MInt j = 0; j < nDist; j++) {
        m_solver->a_oldDistribution(pCellId, j) = m_solver->a_oldDistribution(nghbor, j);
        m_solver->a_distribution(pCellId, j) = m_solver->a_distribution(nghbor, j);
      }
    }
  }

  m_solver->calculateMacroscopicVariables(pCellId, l_rho, l_vv.data());

  m_solver->a_variable(pCellId, PV->RHO) = l_rho;
  for(MInt n = 0; n < nDim; n++) {
    m_solver->a_variable(pCellId, PV->VV[n]) = l_vv[n];
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_lin_thermal(const MInt cellIndex,
                                                                                     const MInt set) {
  TRACE();

  const MInt pCellId = m_solver->m_G0CellList[cellIndex];
  MFloat l_t = m_solver->m_initTemperatureKelvin;
  std::array<MFloat, nDim> uW{};
  getBoundaryVelocityMb(cellIndex, uW.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif
      MFloat s = F0;
      if(m_solver->m_innerEnergy) {
        s = zerothMomentSourceTerm<1>(pCellId, opposite, l_t, uW.data());
      } else if(m_solver->m_totalEnergy) {
        s = zerothMomentSourceTerm<2>(pCellId, opposite, l_t, uW.data());
      } else {
        s = zerothMomentSourceTerm<0>(pCellId, opposite, l_t, uW.data());
      }
      const MFloat sourceTerm = s;

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          const MFloat F2q = F2 * q;
          m_solver->a_oldDistributionThermal(pCellId, opposite) =
              (-F2q * m_solver->a_distributionThermal(pCellId, j))
              + ((F2q - F1) * m_solver->a_oldDistributionThermal(pCellId, j)) + sourceTerm;

        } else {
          // this is the case in which we have no neighbor in the direction we want to set (strange case)
          // can in my opinion only appear if neighbors were set wrong or at an interface
          m_solver->a_oldDistributionThermal(pCellId, opposite) =
              -m_solver->a_distributionThermal(pCellId, j) + sourceTerm;
        }
      }
      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      // else if(bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      // else if( m_wallBoundaryCellList[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistributionThermal(pCellId, opposite) =
            m_solver->a_distributionThermal(pCellId, j) + sourceTerm;
      }
    }
  }

  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

          MFloat s = F0;
          if(m_solver->m_innerEnergy) {
            s = zerothMomentSourceTerm<1>(neighborId, j, l_t, uW.data());
          } else if(m_solver->m_totalEnergy) {
            s = zerothMomentSourceTerm<2>(neighborId, j, l_t, uW.data());
          } else {
            s = zerothMomentSourceTerm<0>(neighborId, j, l_t, uW.data());
          }
          const MFloat sourceTerm = s;
          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(q < 0.5) {
            const MFloat F2q = F2 * (F1 - q);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
              m_solver->a_oldDistributionThermal(neighborId, j) =
                  -(m_solver->a_distributionThermal(neighborId, opposite) / F2q)
                  + (m_solver->a_distributionThermal(neighborId, j) * (F2q - F1) / F2q) + (F1 / F2q) * sourceTerm;
            }
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDistsThermal(pCellId, l_t, F1, uW.data());

    m_solver->a_variable(pCellId, PV->T) = l_t;
    m_solver->a_oldVariable(pCellId, PV->T) = l_t;
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::interpolatedBounceBackMb_Bouzidi_lin_transport(const MInt cellIndex,
                                                                                       const MInt set) {
  TRACE();

  const MInt pCellId = m_solver->m_G0CellList[cellIndex];
  MFloat l_c = m_solver->m_initCon;
  std::array<MFloat, nDim> uW{};
  getBoundaryVelocityMb(cellIndex, uW.data());

  // case 1: boundary cell is inside fluid
  if(m_solver->a_levelSetFunctionMB(pCellId, set) > 0) {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif
      const MFloat sourceTerm = zerothMomentSourceTerm<0>(pCellId, opposite, l_c, uW.data());

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) continue;

      // if q <= 0.5, perform bounce back
      if(q <= 0.5) {
        if(m_solver->a_hasNeighbor(pCellId, opposite) != 0) {
          const MFloat F2q = F2 * q;
          m_solver->a_oldDistributionTransport(pCellId, opposite) =
              (-F2q * m_solver->a_distributionTransport(pCellId, j))
              + ((F2q - F1) * m_solver->a_oldDistributionTransport(pCellId, j)) + sourceTerm;

        } else {
          // this is the case in which we have no neighbor in the direction we want to set (strange case)
          // can in my opinion only appear if neighbors were set wrong or at an interface
          m_solver->a_oldDistributionTransport(pCellId, opposite) =
              -m_solver->a_distributionTransport(pCellId, j) + sourceTerm;
        }
      }
      // this is the case in which we do not have a neighbor and no cut, do simple bounce back (strange case)
      // can in my opinoin only appear if something has gone wrong either with the grid or with the distance calculation
      // else if(bndCells[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      // else if( m_wallBoundaryCellList[cellId].m_distances[j] > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0)
      else if(q > 0.5 && m_solver->a_hasNeighbor(pCellId, j) == 0) {
        m_solver->a_oldDistributionTransport(pCellId, opposite) =
            m_solver->a_distributionTransport(pCellId, j) + sourceTerm;
      }
    }
  }

  // case 2: boundary cell is outside fluid
  else {
    for(MInt j = 0; j < nDist - 1; j++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat q = m_distances[cellIndex][j];
      const MInt opposite = m_solver->m_oppositeDist[j];
#else
      const MFloat q = getDistanceMb(pCellId, cellIndex, j);
      const MInt opposite = Ld::oppositeDist(j);
#endif

      // do we have a neighbor at all?
      if(m_solver->a_hasNeighbor(pCellId, j)) {
        // make sure that the neighbor is not a halo cell
        const MInt neighborId = m_solver->c_neighborId(pCellId, j);
        if(!m_solver->a_hasProperty(neighborId, Cell::IsHalo)) {
          const MBool nbndid = m_solver->a_isG0CandidateOfSet(neighborId, (set - m_solver->m_levelSetId));

          const MFloat sourceTerm = zerothMomentSourceTerm<0>(neighborId, j, l_c, uW.data());
          // if q < 0.5, perform bounce back (this is meant from the outer cell, the inner cell then has q >= 0.5)
          if(q < 0.5) {
            const MFloat F2q = F2 * (F1 - q);

            // check if my neighbor is an inside cell or a boundary cell with the cell center inside
            if(!nbndid || m_solver->a_levelSetFunctionMB(neighborId, set) > 0) {
              m_solver->a_oldDistributionTransport(neighborId, j) =
                  -(m_solver->a_distributionTransport(neighborId, opposite) / F2q)
                  + (m_solver->a_distributionTransport(neighborId, j) * (F2q - F1) / F2q) + (F1 / F2q) * sourceTerm;
            }
          }
        }
      }
    }

    // The outer cell does not belong to the flow field, thus its incoming distributions are overwritten.
    // It is possible that there are cells without any cutting velocity! These are considered too.
    m_solver->setEqDistsTransport(pCellId, l_c, uW.data());
    m_solver->a_variable(pCellId, PV->C) = l_c;
    m_solver->a_oldVariable(pCellId, PV->C) = l_c;
  }
}

/**
 * \brief Extrapolates velocity to all moving boundary surfaces
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::extrapolateVelocitiesMb() {
  // TODO labels:LB Calculate velocities from oldDistribution, BUT BE CAREFUL, LAST HALO LAYER NOT SYNCHED!
  // since a_variable is from the beginning of the timestep (pre coll, pre prop)
  const MInt setId = 0;

  for(MInt cellIndex = 0; cellIndex < m_boundaryCellsMb.size(); cellIndex++) {
    for(MInt dist = 0; dist < nDist - 1; dist++) {
      ASSERT(m_boundaryCellsMb.distance(cellIndex, dist) >= 0.0, "boundary cell collector: no valid dist");
    }
  }

  for(MInt cellIndex = 0; cellIndex < m_boundaryCellsMb.size(); cellIndex++) {
    const MInt pCellId = m_boundaryCellsMb.cellId(cellIndex);

    // Find minWallDistDir
    MInt minDistDir = -1;
    MFloat minWallDistance = F2;
    for(MInt dist = 0; dist < nDist - 1; dist++) {
      if(minWallDistance > m_boundaryCellsMb.distance(cellIndex, dist)) {
        minWallDistance = m_boundaryCellsMb.distance(cellIndex, dist);
        minDistDir = dist;
      }
    }

    ASSERT(minDistDir > -1, "minDistDir not found!");

    minWallDistance = getDistanceMb(pCellId, cellIndex, minDistDir);

    ASSERT(minDistDir > -1, "No min dist found.");

    if(m_solver->a_levelSetFunctionMB(pCellId, setId) > 0) {
      // Interface cell is inside the fluid

      if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
        continue;
      }

      const MInt opposite = Ld::oppositeDist(minDistDir);
      const MFloat q = minWallDistance;
      const MInt neighbor = m_solver->c_neighborId(pCellId, opposite);

      if(neighbor > -1) {
        // Neighbor exists: linear extrapolation to surface
        for(MInt n = 0; n < nDim; n++) {
          m_boundaryCellsMb.velocity(cellIndex, n) =
              m_solver->a_variable(pCellId, n) * (1 + q) - m_solver->a_variable(neighbor, n) * q;
        }

      } else {
        // Neighbor doesnt exists: nearest neighbor value set to surface
        std::cerr << "Velocity extrapolation to surface: zero order fallback. (inside)"
                  << "(gid " << m_solver->c_globalId(pCellId) << " )" << std::endl;
        for(MInt n = 0; n < nDim; n++) {
          m_boundaryCellsMb.velocity(cellIndex, n) = m_solver->a_variable(pCellId, n);
        }
      }
    } else {
      // Interface cell is outside the fluid

      const MFloat q = 1 - minWallDistance;
      const MInt neighbor = m_solver->c_neighborId(pCellId, minDistDir);

      const MInt nneighbor = m_solver->c_neighborId(neighbor, minDistDir);

      if(nneighbor > -1) {
        // Next neighbor exists: linear extrapolation to surface
        for(MInt n = 0; n < nDim; n++) {
          m_boundaryCellsMb.velocity(cellIndex, n) =
              m_solver->a_variable(neighbor, n) * (2 - q) - m_solver->a_variable(nneighbor, n) * (1 - q);
        }

      } else {
        for(MInt n = 0; n < nDim; n++) {
          m_boundaryCellsMb.velocity(cellIndex, n) = m_solver->a_variable(neighbor, n);
        }
      }
    }
  }

  // SANITY CHECKS
  for(MInt mbCell = 0; mbCell < m_boundaryCellsMb.size(); mbCell++) {
    const MInt pCellId = m_boundaryCellsMb.cellId(mbCell);
    if(m_solver->a_hasProperty(pCellId, Cell::IsHalo)) {
      continue;
    }

    if(std::isnan(m_boundaryCellsMb.velocity(mbCell, 0)) || std::isnan(m_boundaryCellsMb.velocity(mbCell, 1))) {
      std::cout << "vel nan for "
                << "(gid " << m_solver->c_globalId(pCellId) << " )" << std::endl;
      TERMM(1, "Velocity extrapolation nan");
    }
  }
}

/** \brief  Calculate characetristic values based on LODI-equations at (x_b, t+1)
 *  \author Miro Gondrum
 *  \date   26.03.2021
 *  Izquierdo et al. 2008: https://doi.org/10.1103/PhysRevE.78.046707
 *  But in contrast value is calculated in the cell-center such that equilibrium
 *  can be forced. Also, no relaxation against target rho or u is performed.
 *
 *  \param[in]  cellId  cellId of the boundary cell
 *  \param[out] rho_b   Characteristic density at boundary cell
 *  \param[out] u_b     Charecteristic velocity at boundary cell
 *  \note template parameter 'type': 0-incoming wave zeros, 1-velocity, 2-pressure
 **/
template <MInt nDim, MInt nDist, class SysEqn>
template <MUint type>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calcCharValues(const MInt index, const MInt cellId, MFloat& rho_b,
                                                              MFloat* u_b) {
  TRACE();
#ifdef WAR_NVHPC_PSTL
  MFloat var[nDim + 1] = {F0};
  for(MInt d = 0; d < nDim + 1; d++) {
    var[d] = m_solver->a_variable(cellId, d);
  }
  MFloat* const ptr0 = var;
#else
  MFloat* const ptr0 = &m_solver->a_variable(cellId, 0);
#endif
  const MFloat F1Brho0 = 1.0 / ptr0[PV->RHO];
  MFloat rho_t{};
  MFloat u_t[nDim]{};
  //--calculate time derivative-----------------------------------------------
  // Heubes et al. 2012: http://dx.doi.org/10.1016/j.cam.2013.09.019
  //  Defined switching factor gamma (here factorThomson) which yield pure
  //  Thomson for gamma=1.0 or pure LODI for gamma=0.0 (as proposed by Izquierdo
  //  et al.). By Heubes et al. it was found out that the convex combination of
  //  both approaches using gamma=0.75 is superior.
  constexpr MFloat factorThomson = 1.0;

  constexpr MFloat sigma = 0.59;
  constexpr MFloat rho_trg = 1.0;
  const MFloat k_relax =
      sigma * (1 - m_Ma * m_Ma) * LBCS / (m_domainLength * FFPOW2(m_solver->maxLevel() - m_solver->a_level(cellId)));
  MFloat L_approx[nDim]{};
  switch(type) {
    case 1: { // relax against velocity
      const MInt ind = m_mapSegIdsInOutCnd[index];
      for(MInt i = 0; i < nDim; i++) {
        const MFloat u_trg = m_initialVelocityVecs[ind][i] * m_Ma * LBCS; //* m_bndCells[].m_multiplier;
        L_approx[i] = k_relax * (ptr0[PV->VV[i]] - u_trg);
      }
      break;
    }
    case 2: { // relax against pressure
      const MFloat L_approx_ = k_relax * CSsq * (ptr0[PV->RHO] - rho_trg);
      for(MInt i = 0; i < nDim; i++) {
        L_approx[i] = L_approx_;
      }
      break;
    }
    default: { // otherwise zero incoming wave
      break;
    }
  }
  const MInt pvVv[3] = {PV->VV[0], PV->VV[1], PV->VV[2]};
  for(MInt c = 0; c < nDim; c++) {
    const MInt dir1 = 2 * c;     // in negative axis dir
    const MInt dir2 = 2 * c + 1; // in positive axis dir
    const MBool isTangentialDir = (m_solver->a_hasNeighbor(cellId, dir1) && m_solver->a_hasNeighbor(cellId, dir2));
    if(isTangentialDir) {
      continue;
      const MInt nghbrId1 = m_solver->c_neighborId(cellId, dir1);
      const MInt nghbrId2 = m_solver->c_neighborId(cellId, dir2);
#ifdef WAR_NVHPC_PSTL
      MFloat var_n1[nDim + 1] = {F0};
      for(MInt d = 0; d < nDim + 1; d++) {
        var_n1[d] = m_solver->a_variable(nghbrId1, d);
      }
      MFloat* const ptr1 = var_n1;
      MFloat var_n2[nDim + 1] = {F0};
      for(MInt d = 0; d < nDim + 1; d++) {
        var_n2[d] = m_solver->a_variable(nghbrId2, d);
      }
      MFloat* const ptr2 = var_n2;
#else
      MFloat* const ptr1 = &m_solver->a_variable(nghbrId1, 0);
      MFloat* const ptr2 = &m_solver->a_variable(nghbrId2, 0);
#endif
      const MFloat F1Brho1 = 1.0 / ptr1[PV->RHO];
      const MFloat F1Brho2 = 1.0 / ptr2[PV->RHO];
      // calculate derivatives based on central differences
      // ...(3)--1--0--2--(4)...
      const MFloat rho_x = 0.5 * (ptr2[PV->RHO] - ptr1[PV->RHO]);
      MFloat u_x[nDim];
      for(MInt j = 0; j < nDim; j++) {
        u_x[j] = 0.5 * (ptr2[PV->VV[j]] * F1Brho2 - ptr1[PV->VV[j]] * F1Brho1);
      }
      // time derivative part for this direction
      rho_t -= factorThomson * (ptr0[pvVv[c]] * F1Brho0 * rho_x + u_x[c] * ptr0[PV->RHO]);
      u_t[c] -= factorThomson * (F1Brho0 * CSsq * rho_x);
      for(MInt j = 0; j < nDim; j++) {
        u_t[j] -= factorThomson * (F1Brho0 * ptr0[pvVv[c]] * u_x[j]);
      }
      // Add artificial damping to otherwise undamped central difference
      /*
      constexpr MFloat eps = 1.0e-2;
      if(m_solver->a_hasNeighbor(nghbrId1, dir1) && m_solver->a_hasNeighbor(nghbrId2, dir2)) {
        const MInt nghbrId3 = m_solver->c_neighborId(nghbrId1, dir1);
        const MInt nghbrId4 = m_solver->c_neighborId(nghbrId2, dir2);
        MFloat* const ptr3 = &m_solver->a_variable(nghbrId3, 0);
        MFloat* const ptr4 = &m_solver->a_variable(nghbrId4, 0);
        const MFloat F1Brho3 = 1.0 / ptr3[PV->RHO];
        const MFloat F1Brho4 = 1.0 / ptr4[PV->RHO];
        rho_t -=
            eps * (ptr3[PV->RHO] - 4.0 * ptr1[PV->RHO] + 6.0 * ptr0[PV->RHO] - 4.0 * ptr2[PV->RHO] + ptr4[PV->RHO]);
        for(MInt j = 0; j < nDim; j++) {
          u_t[j] -= eps
                    * (ptr3[PV->VV[j]] * F1Brho3 - 4.0 * ptr1[PV->VV[j]] * F1Brho1 + 6.0 * ptr0[PV->VV[j]] * F1Brho2
                       - 4.0 * ptr2[PV->VV[j]] * F1Brho2 + ptr4[PV->VV[j]] * F1Brho4);
        }
      } else {
        rho_t -= eps * (ptr1[PV->RHO] - F2 * ptr0[PV->RHO] + ptr2[PV->RHO]);
        for(MInt j = 0; j < nDim; j++) {
          u_t[j] -= eps * (ptr1[PV->VV[j]] * F1Brho1 - F2 * ptr0[PV->VV[j]] * F1Brho0 + ptr2[PV->VV[j]] * F1Brho2);
        }
      }
      */
    } else {
      MInt dir;
      if(m_solver->a_hasNeighbor(cellId, dir1)) {
        dir = dir1;
      } else if(m_solver->a_hasNeighbor(cellId, dir2)) {
        dir = dir2;
      } else {
        TERMM(1, "Boundary cell with no neighbor in normal direction?");
      }
      const MInt nghbrId1 = m_solver->c_neighborId(cellId, dir);
      const MInt nghbrId2 = m_solver->c_neighborId(nghbrId1, dir);
      const MInt sign = (dir % 2 == 0) ? 1 : -1; // negativ for -x bndry
#ifdef WAR_NVHPC_PSTL
      MFloat var_n1[nDim + 1];
      for(MInt d = 0; d < nDim + 1; d++) {
        var_n1[d] = m_solver->a_variable(nghbrId1, d);
      }
      MFloat* const ptr1 = var_n1;
      MFloat var_n2[nDim + 1];
      for(MInt d = 0; d < nDim + 1; d++) {
        var_n2[d] = m_solver->a_variable(nghbrId2, d);
      }
      MFloat* const ptr2 = var_n2;
#else
      MFloat* const ptr1 = &m_solver->a_variable(nghbrId1, 0);
      MFloat* const ptr2 = &m_solver->a_variable(nghbrId2, 0);
#endif
      const MFloat F1Brho1 = 1.0 / ptr1[PV->RHO];
      const MFloat F1Brho2 = 1.0 / ptr2[PV->RHO];
      // calculate derivatives based on one-sided differences
      // |--0--1--2--...
      const MFloat rho_x = sign * (F3B2 * ptr0[PV->RHO] - F2 * ptr1[PV->RHO] + F1B2 * ptr2[PV->RHO]);
      MFloat u_x[nDim];
      for(MInt j = 0; j < nDim; j++) {
        // if dir is pointing in negative direction (dir%2==0) -> backward difference, otherwise forward
        u_x[j] =
            sign
            * (F3B2 * ptr0[PV->VV[j]] * F1Brho0 - F2 * ptr1[PV->VV[j]] * F1Brho1 + F1B2 * ptr2[PV->VV[j]] * F1Brho2);
      }
      // assemble L - wave amplitude variation, while vanishing incoming waves
      MFloat L[nDim + 1]{}; // for 3D corresponding to {u_n-c, u_t1, u_t2, u_n+c}
      const MBool flowPointsInwards = (ptr0[pvVv[c]] * sign < 0);
      if(flowPointsInwards) {
        MInt id = 1;
        for(MInt j = 0; j < nDim; j++) {
          if(j != c) {
            L[id] = L_approx[j];
            id++;
          }
        }
      } else {
        MInt id = 1;
        for(MInt j = 0; j < nDim; j++) {
          if(j != c) {
            L[id] = ptr0[pvVv[c]] * F1Brho0 * u_x[j];
            id++;
          }
        }
      }
      if(sign > 0) {
        // boundary outward normal points in positive axis direction
        L[0] = L_approx[c];
        L[nDim] = (ptr0[pvVv[c]] * F1Brho0 + LBCS) * (CSsq * rho_x + LBCS * ptr0[PV->RHO] * u_x[c]);
      } else {
        // boundary outward normal points in negative axis direction
        L[0] = (ptr0[pvVv[c]] * F1Brho0 - LBCS) * (CSsq * rho_x - LBCS * ptr0[PV->RHO] * u_x[c]);
        L[nDim] = L_approx[c];
      }
      // time derivative part for this direction
      rho_t -= F1B2 * F1BCSsq * (L[nDim] + L[0]);
      MInt id = 1;
      for(MInt j = 0; j < nDim; j++) {
        if(j != c) {
          u_t[j] -= L[id];
          id++;
        } else {
          u_t[c] -= F1B2 * F1BCS * F1Brho0 * (L[nDim] - L[0]);
        }
      }
    }
  }
  //--calc characteristic values----------------------------------------------
  rho_b = m_solver->a_variable(cellId, PV->RHO) + rho_t;
  for(MInt j = 0; j < nDim; j++) {
    u_b[j] = m_solver->a_variable(cellId, PV->VV[j]) * F1Brho0 + u_t[j];
  }
}

/** \brief  Calculate characetristic values based on LODI-equations at (x_b+1/2, t-1/2)
 *  \author Miro Gondrum
 *  \date   01.03.2021
 *  Izquierdo et al. 2008: https://doi.org/10.1103/PhysRevE.78.046707
 *
 *  \param[in]  index     Index of the boundary condition
 *  \param[in]  direction Direction the boundary outward normal is pointing to
 *  \param[in]  bndCellId cellId of the cell in boundary context
 *  \param[out] rho_b     Characteristic density at boundary cell + dx/2 for t-1/2
 *  \param[out] u_b       Charecteristic velocity at boundary cell + dx/2 for t-1/2
 *
 *  \note This function does not provide valid values for edges and corner, yet.
 **/
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::calcCharValuesOnFace(const MInt index, const MInt /*direction*/,
                                                                    const MInt bndCellId, MFloat& rho_b, MFloat* u_b) {
  constexpr MFloat deltaT = 1;
  constexpr MFloat sigma = 0.59;
  constexpr MFloat kappa = 1.0; // in original kappa=1.2 is used in combination with corrected LBE
  constexpr MFloat sqrtKappa = 1.0;
  // speed of sound c_s = sqrt(kappa*p/rho) where kappa = 1.2 is used at the boundary
  // const MFloat sqrtKappa = sqrt(kappa);

  MFloat* const formerVariables = this->m_bndCndData[index].data;
  const MInt noBcVars = this->m_bndCndData[index].noVars;
  const MInt offset = m_bndCndOffsets[index];

  const MInt currentId = m_bndCells[bndCellId].m_cellId;
  const MInt iMo = (bndCellId - offset) * noBcVars;

  const MFloat k_relax =
      sigma * (1 - m_Ma * m_Ma) * sqrtKappa * LBCS
      / (m_domainLength
         * FFPOW2(m_solver->maxLevel()
                  - m_solver->a_level(currentId))); // m_domainLength represents the distance between the inlet
                                                    // and the outlet on the highest level of refinement

  const MFloat rho_trg = (m_densityFluctuations) ? 0.0 : 1.0; // target pressure at boundary
  const MFloat rho_offset = (m_densityFluctuations) ? 1.0 : 0.0;

  // ...
  auto subTimeDerivatives = [&](MInt direction2, MFloat& rho_t, MFloat* u_t) {
    const MInt nghbrId = m_solver->c_neighborId(currentId, Ld::oppositeDist(direction2));
    // x_b is the location between fluid and boundary center. x_b-i is i-th fluid cell
    // compute derivatives:
    // (x_b,   t-1) : formerVariables
    // (x_b-1, t-1) : [a_variable(currentId,..) + a_oldVariable(currentId,..)]/2
    // (x_b-2, t-1) : [a_variable(nghbrId,..)   + a_oldVariable(nghbrId,..)]  /2
    const MFloat der_rho =
        F1B3
        * (8 * formerVariables[iMo + PV->RHO]
           - 9 * F1B2 * (m_solver->a_variable(currentId, PV->RHO) + m_solver->a_oldVariable(currentId, PV->RHO))
           + F1B2 * (m_solver->a_variable(nghbrId, PV->RHO) + m_solver->a_oldVariable(nghbrId, PV->RHO)));
    const MFloat der_p = kappa * CSsq * der_rho;
    const MFloat der_u =
        F1B3
        * (8 * formerVariables[iMo + PV->U]
           - 9 * F1B2 * (m_solver->a_variable(currentId, PV->U) + m_solver->a_oldVariable(currentId, PV->U))
           + F1B2 * (m_solver->a_variable(nghbrId, PV->U) + m_solver->a_oldVariable(nghbrId, PV->U)));
    const MFloat der_v =
        F1B3
        * (8 * formerVariables[iMo + PV->V]
           - 9 * F1B2 * (m_solver->a_variable(currentId, PV->V) + m_solver->a_oldVariable(currentId, PV->V))
           + F1B2 * (m_solver->a_variable(nghbrId, PV->V) + m_solver->a_oldVariable(nghbrId, PV->V)));
    const MFloat der_w =
        F1B3
        * (8 * formerVariables[iMo + PV->W]
           - 9 * F1B2 * (m_solver->a_variable(currentId, PV->W) + m_solver->a_oldVariable(currentId, PV->W))
           + F1B2 * (m_solver->a_variable(nghbrId, PV->W) + m_solver->a_oldVariable(nghbrId, PV->W)));

    // L: wave amplitude variations determined for (x_b, t-1/2)
    // incoming ones are approximated
    //                           | BC+:       BC-:      // in Izquierdo2008
    MFloat L1 = F0; // u_n - c_s | incoming   outgoing  // L1
    MFloat L2 = F0; // energy    | outgoing   incoming  // L3 (=0 for isothermal)
    MFloat L3 = F0; // u_t1      | outgoing   incoming  // L2.1
    MFloat L4 = F0; // u_t2      | outgoing   incoming  // L2.2 (since it's 3D)
    MFloat L5 = F0; // u_n + c_s | outgoing   incoming  // L4

    switch(direction2) {
      case 0:
        // negative derivatives due to negative direction!
        L1 = (formerVariables[iMo + PV->U] - sqrtKappa * LBCS)
             * (-der_p - (formerVariables[iMo + PV->RHO] + rho_offset) * (-der_u) * sqrtKappa * LBCS);
        L2 = 0.0;
        L3 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        L4 = L3;
        L5 = L3;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        u_t[1] += L3;
        u_t[2] += L4;
        break;
      case 1:
        L1 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        // L2 = formerVariables[iMo + PV->U] * (POW2(LBCS) * der_rho - der_p); // = 0 for isothermal LB
        L2 = 0.0;
        L3 = formerVariables[iMo + PV->U] * der_v;
        L4 = formerVariables[iMo + PV->U] * der_w;
        L5 = (formerVariables[iMo + PV->U] + sqrtKappa * LBCS)
             * (der_p + (formerVariables[iMo + PV->RHO] + rho_offset) * der_u * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        u_t[1] += L3;
        u_t[2] += L4;
        break;
      case 2:
        // negative derivatives due to negative direction!
        L1 = (formerVariables[iMo + PV->V] - sqrtKappa * LBCS)
             * (-der_p - (formerVariables[iMo + PV->RHO] + rho_offset) * (-der_v) * sqrtKappa * LBCS);
        L2 = 0.0;
        L3 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        L4 = L3;
        L5 = L3;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += L3;
        u_t[1] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        u_t[2] += L4;
        break;
      case 3:
        L1 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        L2 = 0.0;
        L3 = formerVariables[iMo + PV->V] * der_u;
        L4 = formerVariables[iMo + PV->V] * der_w;
        L5 = (formerVariables[iMo + PV->V] + sqrtKappa * LBCS)
             * (der_p + (formerVariables[iMo + PV->RHO] + rho_offset) * der_v * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += L3;
        u_t[1] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        u_t[2] += L4;
        break;
      case 4:
        // negative derivatives due to negative direction!
        L1 = (formerVariables[iMo + PV->W] - sqrtKappa * LBCS)
             * (-der_p - (formerVariables[iMo + PV->RHO] + rho_offset) * (-der_w) * sqrtKappa * LBCS);
        L2 = 0.0;
        L3 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        L4 = L3;
        L5 = L3;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += L3;
        u_t[1] += L4;
        u_t[2] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        break;
      case 5:
        L1 = k_relax * kappa * CSsq * (formerVariables[iMo + PV->RHO] - rho_trg);
        L2 = 0.0;
        L3 = formerVariables[iMo + PV->W] * der_u;
        L4 = formerVariables[iMo + PV->W] * der_v;
        L5 = (formerVariables[iMo + PV->W] + sqrtKappa * LBCS)
             * (der_p + (formerVariables[iMo + PV->RHO] + rho_offset) * der_w * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u_t[0] += L3;
        u_t[1] += L4;
        u_t[2] += (L5 - L1) * F1B2 * (F1BCS / sqrtKappa) / (formerVariables[iMo + PV->RHO] + rho_offset);
        break;
      default:
        TERMM(1, "Unrealistic case, please check your Boundary Conditions!!!");
        break;
    }
    // compute rho using the characteristic wave amplitudes
    rho_t = (F1BCSsq / kappa) * (F1B2 * (L5 + L1) + L2);
  };
  MFloat rho_t = 0;
  MFloat u_t[3]{}; // TODO labels:LB dxqy: only for 3D defined, yet
  for(MInt d = 0; d < nDim * 2; d++) {
    if(m_solver->a_hasNeighbor(currentId, d) == 0) {
      subTimeDerivatives(d, rho_t, u_t);
    }
  }
  // forward euler step
  rho_b = formerVariables[iMo + PV->RHO] - deltaT * rho_t;
  u_b[0] = formerVariables[iMo + PV->U] - deltaT * u_t[0];
  u_b[1] = formerVariables[iMo + PV->V] - deltaT * u_t[1];
  u_b[2] = formerVariables[iMo + PV->W] - deltaT * u_t[2];
  // save variables in virtual boundary cell for next timestep
  formerVariables[iMo + PV->RHO] = rho_b;
  formerVariables[iMo + PV->U] = u_b[0];
  formerVariables[iMo + PV->V] = u_b[1];
  formerVariables[iMo + PV->W] = u_b[2];
}

/** \LBBC{sec_LBBC_charVelocity, charVelocity, 104x}
 *  in-/outflow bc
 *
 *  Non-reflecting characteristic BC
 *  Izquierdo et al. 2008: https://doi.org/10.1103/PhysRevE.78.046707
 *  The velocity is set to Ma*c_s, in opposite boundary normal direction<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::charVelocity(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::charVelocity(MInt index, MInt direction) {
  TRACE();

  constexpr MFloat kappa = 1.0; // in original kappa=1.2 is used in combination with corrected LBE

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0
       || m_solver->a_isHalo(currentId) || !m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(direction))) {
      continue;
    }

    //--determine macroscopic values based on LODI------------------------------
    MFloat rho_b, u_b[3];
    calcCharValuesOnFace(index, direction, i, rho_b, &u_b[0]);

    //--perform bounce back with velocity term----------------------------------
    const MFloat rho_offset = (m_densityFluctuations) ? 1.0 : 0.0;
    MFloat c[2 * nDim];
    for(MInt n = 0; n < nDim; n++) {
      c[2 * n] = -u_b[n];
      c[2 * n + 1] = u_b[n];
    }

    for(MInt j = 0; j < Ld::dxQyFld(); j++) {
      const MInt dir = Ld::componentFld(direction, j);
      if(m_solver->a_hasNeighbor(currentId, dir) == 0) {
        MInt tpIndex;
        MFloat tmpUW;
        if(dir < Ld::distFld(0)) {
          tmpUW = c[dir];
          tpIndex = 1;
        } else {
          if(dir < (Ld::distFld(0) + Ld::distFld(1))) {
            const MInt k = dir - Ld::distFld(0);
            tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          } else {
            const MInt k = dir - (Ld::distFld(0) + Ld::distFld(1));
            tmpUW = (c[Ld::mFld2(3 * k)] + c[Ld::mFld2(3 * k + 1)] + c[Ld::mFld2(3 * k + 2)]);
            tpIndex = 3;
          }
        }

        m_solver->a_oldDistribution(currentId, Ld::oppositeDist(dir)) =
            m_solver->a_distribution(currentId, dir)
            - 2.0 * Ld::tp(tpIndex) * (rho_b + rho_offset) * (F1BCSsq / kappa) * tmpUW;
      }
    }
  }
}

template <>
inline void LbBndCndDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::charVelocity(MInt index, MInt direction) {
  TRACE();

  constexpr MInt nDim = 2;

  MFloat L1, L2, L3, L4;
  L1 = L2 = L3 = L4 = F0;
  MFloat rho, u1, u2, ub;
  MFloat deltaT = 1;
  MFloat sigma = 1;

  MFloat kappa = 1.2;
  MFloat sqrtKappa = sqrt(1.2); // speed of sound c_s = sqrt(kappa*p/rho) where kappa = 1.2 is used at the boundary
  MInt currentID, nghbrID, tpIndex, k;
  // MFloat tmpWidth = m_referenceLength * (m_solver->c_cellLengthAtLevel(m_solver->m_solver->maxLevel()]);
  MFloat k_relax = sigma * (1 - m_Ma * m_Ma) * sqrtKappa * LBCS / m_solver->c_cellLengthAtLevel(0);
  MFloat der_u, der_v, der_p;
  MFloat c[4], tmpUW;

  MFloat* const formerVariables = this->m_bndCndData[index].data;
  const MInt offset = m_bndCndOffsets[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    currentID = m_bndCells[i].m_cellId;
    nghbrID = m_solver->c_neighborId(currentID, Ld::oppositeDist(direction));
    constexpr MInt noBcVars = nDim + 1;
    const MInt iMo = (i - offset) * noBcVars;

    // parabolic inflow
    // ub =  1.5 * m_Ma * sqrtKappa * LBCS * (1.0 - (2*m_solver->a_coordinate(currentID,
    // 1)/tmpWidth)*(2*m_solver->a_coordinate(currentID, 1)/tmpWidth));

    // uniform inflow
    ub = m_Ma * LBCS;

    // calculate incoming wave amplitudes. to obtain L1 spatial derivatives der_p and der_u are needed.
    // t^= t + 1/2 => t^ - 1 = t - 1/2 => Averaging of m_variables and m_oldVariables

    if(direction == 0) {
      der_p = F1B3 * CSsq
              * (-8 * formerVariables[iMo + PV->RHO]
                 + 9 * F1B2 * (m_solver->a_variable(currentID, PV->RHO) + m_solver->a_oldVariable(currentID, PV->RHO))
                 - F1B2 * (m_solver->a_oldVariable(nghbrID, PV->RHO) + m_solver->a_variable(nghbrID, PV->RHO)));


      der_u = F1B3
              * (-8 * formerVariables[iMo + PV->U]
                 + 9 * F1B2 * (m_solver->a_variable(currentID, PV->U) + m_solver->a_oldVariable(currentID, PV->U))
                 - F1B2 * (m_solver->a_variable(nghbrID, PV->U) + m_solver->a_oldVariable(nghbrID, PV->U)));


      der_v = F1B3
              * (-8 * formerVariables[iMo + PV->V]
                 + 9 * F1B2 * (m_solver->a_variable(currentID, PV->V) + m_solver->a_oldVariable(currentID, PV->V))
                 - F1B2 * (m_solver->a_variable(nghbrID, PV->V) + m_solver->a_oldVariable(nghbrID, PV->V)));


      L1 = (formerVariables[iMo + PV->U] - sqrtKappa * LBCS)
           * (der_p - formerVariables[iMo + PV->RHO] * der_u * sqrtKappa * LBCS);
      L2 = k_relax * F1B3 * (formerVariables[iMo + PV->U] - 0); // set u_y=0
      L3 = L2;
      L4 = L2;
    } else if(direction == 1) {
      der_p = CSsq * F1B3
              * (8 * formerVariables[iMo + PV->RHO]
                 - 9 * F1B2 * (m_solver->a_variable(currentID, PV->RHO) + m_solver->a_oldVariable(currentID, PV->RHO))
                 + F1B2 * (m_solver->a_oldVariable(nghbrID, PV->RHO) + m_solver->a_variable(nghbrID, PV->RHO)));


      der_u = F1B3
              * (8 * formerVariables[iMo + PV->U]
                 - 9 * F1B2 * (m_solver->a_variable(currentID, PV->U) + m_solver->a_oldVariable(currentID, PV->U))
                 + F1B2 * (m_solver->a_variable(nghbrID, PV->U) + m_solver->a_oldVariable(nghbrID, PV->U)));


      der_v = F1B3
              * (8 * formerVariables[iMo + PV->V]
                 - 9 * F1B2 * (m_solver->a_variable(currentID, PV->V) + m_solver->a_oldVariable(currentID, PV->V))
                 + F1B2 * (m_solver->a_variable(nghbrID, PV->V) + m_solver->a_oldVariable(nghbrID, PV->V)));

      L1 = k_relax * F1B3 * (formerVariables[iMo + PV->U] - ub);
      L2 = formerVariables[iMo + PV->U] * der_v;
      L3 = 0;
      L4 = (formerVariables[iMo + PV->U] + sqrtKappa * LBCS)
           * (der_p + formerVariables[iMo + PV->RHO] * der_u * sqrtKappa * LBCS);
    } else {
      TERMM(1, "Direction can only be 0 or 1 in D2Q9!");
    }

    // compute u1,u2 and rho by using the characteristic wave amplitudes
    u1 = formerVariables[iMo + PV->U] - deltaT * F1B2 / formerVariables[iMo + PV->RHO] * F1BCS / sqrtKappa * (L4 - L1);
    u2 = formerVariables[iMo + PV->V] - deltaT * L2;
    rho = formerVariables[iMo + PV->RHO] - deltaT * F1B2 / kappa * F1BCSsq * (L4 + L1)
          - F1BCSsq * L3 / kappa; // where 1.2 is kappa²

    c[0] = -ub;
    c[1] = ub;
    c[2] = F0;
    c[3] = F0;

    for(MInt j = 0; j < 8; j++) {
      if(m_solver->a_hasNeighbor(currentID, j) == 0) {
        if(j < 4) {
          tmpUW = c[j];
          tpIndex = 1;
        } else {
          if(j < 8) {
            k = j - Ld::distFld(0);
            tmpUW = (c[Ld::mFld1(2 * k)] + c[Ld::mFld1(2 * k + 1)]);
            tpIndex = 2;
          }
        }

        m_solver->a_oldDistribution(currentID, Ld::oppositeDist(j)) =
            m_solver->a_distribution(currentID, j)
            - 2.0 * Ld::tp(tpIndex) * rho * F1BCSsq * tmpUW; // the minus corresponds to the opposite direction
      }
    }

    // // only for bounce back in corners!!! (upper and lower boundary need to be walls)
    // if (m_solver->a_hasNeighbor(currentID, 3) == 0) {
    // 	m_solver->a_oldDistribution(currentID, 6) = m_solver->a_distribution(currentID, 4);
    // 	m_solver->a_oldDistribution(currentID, 2) = m_solver->a_distribution(currentID, 3);
    // 	m_solver->a_oldDistribution(currentID, 5) = m_solver->a_distribution(currentID, 7);
    // }

    // if (m_solver->a_hasNeighbor(currentID, 2) == 0) {
    // 	m_solver->a_oldDistribution(currentID, 4) = m_solver->a_distribution(currentID, 6);
    // 	m_solver->a_oldDistribution(currentID, 3) = m_solver->a_distribution(currentID, 2);
    // 	m_solver->a_oldDistribution(currentID, 7) = m_solver->a_distribution(currentID, 5);
    // }

    // set new boundary variables
    formerVariables[iMo + 0] = u1;
    formerVariables[iMo + 1] = u2;
    formerVariables[iMo + 2] = rho;
  }
}

/** \LBBC{sec_LBBC_bc10046, bc10046, 1046}
 *  \author Miro Gondrum
 *
 *  Non-reflecting BC based on LODI - setting equilibrium
 *
 *  Izquierdo et al. 2008: https://doi.org/10.1103/PhysRevE.78.046707<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc10046(MInt index)
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc10046(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0
       || m_solver->a_isHalo(currentId)) {
      continue;
    }
    //--LODI calculation--------------------------------------------------------
    MFloat rho_b, u_b[nDim];
    calcCharValues<1>(index, currentId, rho_b, &u_b[0]);
    //--set equilibrium state---------------------------------------------------
#ifdef WAR_NVHPC_PSTL
    MFloat oldDist[nDist] = {F0};
    for(MInt d = 0; d < nDist; d++) {
      oldDist[d] = m_solver->a_oldDistribution(currentId, d);
    }
    sysEqn().calcEqDists(rho_b, &u_b[0], oldDist);
    for(MInt d = 0; d < nDist; d++) {
      m_solver->a_oldDistribution(currentId, d) = oldDist[d];
    }
#else
    sysEqn().calcEqDists(rho_b, &u_b[0], &m_solver->a_oldDistribution(currentId, 0));
#endif
  }
}

/** \LBBC{sec_LBBC_charPressure, charPressure, 404x}
 *  in-/outflow bc
 *
 *  Non-reflecting characteristic BC
 *  Izquierdo et al. 2008: https://doi.org/10.1103/PhysRevE.78.046707
 *  The density is relaxed to 1.0 at the in- or outflow.
 *  The characteristics are solved at x_b = x+dq/2 and the solution is saved in m_formerVariables for the next
 * timestep.<br> LbBndCndDxQy<nDim, nDist, SysEqn>::charPressure(MInt index)
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::charPressure(MInt index, const MInt direction) {
  TRACE();

  // This boundary condition constitutes the following steps:
  // 1) Determining the macroscopic variables at the outlet based on LODI
  //    (linearized one-dimensional invsicd) equation : rho_w, u_w
  // 2) Performing pressure antibounce-back with previosly calculated rho_w, u_w

  constexpr MFloat kappa = 1.0; // in original kappa=1.2 is used in combination with corrected LBE

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0
       || m_solver->a_isHalo(currentId) || !m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(direction))) {
      continue;
    }

    //--determine macroscopic values based on LODI------------------------------
    MFloat rho_b, u_b[3];
    calcCharValuesOnFace(index, direction, i, rho_b, &u_b[0]);

    //--PAB---------------------------------------------------------------------
    //-Calculate symmetric part of eq-distributions at the boundary
    const MFloat squaredU_b = std::inner_product(&u_b[0], &u_b[nDim], &u_b[0], .0);

    const MInt tmpDistId0 = Ld::distFld(0);
    const MInt tmpDistId1 = Ld::distFld(1);
    const MInt tmpDistId2 = Ld::distFld(2);

    MFloat b[2 * nDim];
    for(MInt n = 0; n < nDim; n++) {
      b[2 * n] = -u_b[n];
      b[2 * n + 1] = u_b[n];
    }

    MFloat eqSym[nDist - 1];
    // Calculation of eq. distributions for directions with only one component
    for(MInt j = 0; j < tmpDistId0; j++) {
      eqSym[j] = Ld::tp(1) * (F1BCSsq / kappa)
                 * (rho_b * (CSsq * kappa) + b[j] * b[j] * (F1BCSsq / kappa) * F1B2 - squaredU_b * F1B2);
    }
    // Calculation of eq. distributions for directions with two components
    for(MInt j = 0; j < tmpDistId1; j++) {
      const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      eqSym[tmpDistId0 + j] = Ld::tp(2) * (F1BCSsq / kappa)
                              * (rho_b * (CSsq * kappa) + tmp * tmp * (F1BCSsq / kappa) * F1B2 - squaredU_b * F1B2);
    }
    // Calculation of eq. distributions for directions with three components
    {
      const MInt tmpDistId = tmpDistId0 + tmpDistId1;
      for(MInt j = 0; j < tmpDistId2; j++) {
        const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
        eqSym[tmpDistId + j] = Ld::tp(3) * (F1BCSsq / kappa)
                               * (rho_b * (CSsq * kappa) + tmp * tmp * (F1BCSsq / kappa) * F1B2 - squaredU_b * F1B2);
      }
    }

    // Calculate eq-distributions at boundary cell center (for correction-term)
    //------------------------------------------------------------
#ifdef WAR_NVHPC_PSTL
    MFloat u[nDim] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = m_solver->a_variable(currentId, d);
    }
    const MFloat* const l_uu = u;
#else
    const MFloat* const l_uu = &m_solver->a_variable(currentId, PV->VV[0]);
#endif
    const MFloat l_rho = m_solver->a_variable(currentId, PV->RHO);
    for(MInt n = 0; n < nDim; n++) {
      b[2 * n] = -l_uu[n];
      b[2 * n + 1] = l_uu[n];
    }
    const MFloat squaredU = std::inner_product(&l_uu[0], &l_uu[nDim], &l_uu[0], .0);

    MFloat eqSymC[nDist - 1];
    // Calculation of eq. distributions for directions with only one component
    for(MInt j = 0; j < tmpDistId0; j++) {
      eqSymC[j] =
          Ld::tp(1) * (F1BCSsq / kappa) * (l_rho * (CSsq * kappa) + b[j] * b[j] * F1BCSsq * F1B2 - squaredU * F1B2);
    }

    // Calculation of eq. distributions for directions with two components
    {
      const MInt tmpDistId = tmpDistId0;
      for(MInt j = 0; j < tmpDistId1; j++) {
        const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
        eqSymC[tmpDistId + j] = Ld::tp(2) * (F1BCSsq / kappa)
                                * (l_rho * (CSsq * kappa) + tmp * tmp * (F1BCSsq / kappa) * F1B2 - squaredU * F1B2);
      }
    }

    // Calculation of eq. distributions for directions with three components
    {
      const MInt tmpDistId = tmpDistId0 + tmpDistId1;
      for(MInt j = 0; j < tmpDistId2; j++) {
        const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
        eqSymC[tmpDistId + j] = Ld::tp(3) * (F1BCSsq / kappa)
                                * (l_rho * (CSsq * kappa) + tmp * tmp * (F1BCSsq / kappa) * F1B2 - squaredU * F1B2);
      }
    }

    // symmetric distribution function for 2nd order correction term
    //---------------------------------------------------
    MFloat meanDist[nDist - 1];
    // Calculate average for opposite directions
    for(MInt j = 0; j < nDist - 1; j++) {
      meanDist[j] =
          F1B2
          * (m_solver->a_oldDistribution(currentId, j) + m_solver->a_oldDistribution(currentId, Ld::oppositeDist(j)));
    }

    // set missing incoming distributions (normal to the boundary)
    //--------------------------------------
    const MFloat omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)));

    for(MInt j = 0; j < Ld::dxQyFld(); j++) {
      const MInt dir = Ld::componentFld(direction, j);
      m_solver->a_oldDistribution(currentId, Ld::oppositeDist(dir)) =
          -m_solver->a_distribution(currentId, dir) + 2 * eqSym[dir] + (2 - omega) * (meanDist[dir] - eqSymC[dir]);
    }
  }
}

template <>
inline void LbBndCndDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::charPressure(MInt index, MInt direction) {
  TRACE();

  constexpr MInt nDim = 2;

  MFloat L1, L2, L3, L4;
  L1 = L2 = L3 = L4 = F0;
  MFloat rho, u1, u2, rho_b, u1c, u2c, rhoc;
  MFloat deltaT = 1;
  MFloat sigma = 1;

  MFloat F9B2 = 9.0 / 2.0;

  MFloat kappa = 1.2;
  MFloat sqrtKappa = sqrt(kappa); // speed of sound c_s = sqrt(kappa*p/rho) where kappa = 1.2 is used at the boundary
  MInt currentID, nghbrID;
  MFloat k_relax;
  MFloat der_u, der_p, der_v;
  MFloat eqD[9];
  MFloat eqDc[8];
  MFloat meanDist[8];

  MFloat* const formerVariables = this->m_bndCndData[index].data;
  const MInt offset = m_bndCndOffsets[index];

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    currentID = m_bndCells[i].m_cellId;
    constexpr MInt noBcVars = nDim + 1;
    const MInt iMo = (i - offset) * noBcVars;

    k_relax = sigma * (1 - m_Ma * m_Ma) * sqrtKappa * LBCS
              / (m_domainLength
                 * FFPOW2(m_solver->maxLevel()
                          - m_solver->a_level(currentID))); // m_domainLength represents the distance between the inlet
                                                            // and the outlet on the highest level of refinement

    m_omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentID) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentID)));

    // if t<5 set formerVariables equal to variables
    if(globalTimeStep < 5) {
      formerVariables[iMo + 0] = m_solver->a_variable(currentID, PV->U);
      formerVariables[iMo + 1] = m_solver->a_variable(currentID, PV->V);
      formerVariables[iMo + 2] = m_solver->a_variable(currentID, PV->RHO);
    } else {
      nghbrID = m_solver->c_neighborId(currentID, Ld::oppositeDist(direction));

      if(m_solver->a_hasNeighbor(currentID, direction) == 0) {
        rho_b = 1;

        if(direction == 1) {
          der_p =
              CSsq * F1B3
              * (8 * formerVariables[iMo + PV->RHO]
                 - 9 * F1B2 * (m_solver->a_variable(currentID, PV->RHO) + m_solver->a_oldVariable(currentID, PV->RHO))
                 + F1B2 * (m_solver->a_oldVariable(nghbrID, PV->RHO) + m_solver->a_variable(nghbrID, PV->RHO)));


          der_u = F1B3
                  * (8 * formerVariables[iMo + PV->U]
                     - 9 * F1B2 * (m_solver->a_variable(currentID, PV->U) + m_solver->a_oldVariable(currentID, PV->U))
                     + F1B2 * (m_solver->a_variable(nghbrID, PV->U) + m_solver->a_oldVariable(nghbrID, PV->U)));


          der_v = F1B3
                  * (8 * formerVariables[iMo + PV->V]
                     - 9 * F1B2 * (m_solver->a_variable(currentID, PV->V) + m_solver->a_oldVariable(currentID, PV->V))
                     + F1B2 * (m_solver->a_variable(nghbrID, PV->V) + m_solver->a_oldVariable(nghbrID, PV->V)));

          L1 = k_relax * F1B3 * (formerVariables[iMo + PV->RHO] - rho_b);
          L2 = formerVariables[iMo + PV->U] * der_v;
          L3 = 0;
          L4 = (formerVariables[iMo + PV->U] + sqrtKappa * LBCS)
               * (der_p + formerVariables[iMo + PV->RHO] * der_u * sqrtKappa * LBCS);
        } else if(direction == 0) {
          der_p =
              F1B3 * CSsq
              * (-8 * formerVariables[iMo + PV->RHO]
                 + 9 * F1B2 * (m_solver->a_variable(currentID, PV->RHO) + m_solver->a_oldVariable(currentID, PV->RHO))
                 - F1B2 * (m_solver->a_oldVariable(nghbrID, PV->RHO) + m_solver->a_variable(nghbrID, PV->RHO)));


          der_u = F1B3
                  * (-8 * formerVariables[iMo + PV->U]
                     + 9 * F1B2 * (m_solver->a_variable(currentID, PV->U) + m_solver->a_oldVariable(currentID, PV->U))
                     - F1B2 * (m_solver->a_variable(nghbrID, PV->U) + m_solver->a_oldVariable(nghbrID, PV->U)));


          der_v = F1B3
                  * (-8 * formerVariables[iMo + PV->V]
                     + 9 * F1B2 * (m_solver->a_variable(currentID, PV->V) + m_solver->a_oldVariable(currentID, PV->V))
                     - F1B2 * (m_solver->a_variable(nghbrID, PV->V) + m_solver->a_oldVariable(nghbrID, PV->V)));


          L1 = (formerVariables[iMo + PV->U] - sqrtKappa * LBCS)
               * (der_p - formerVariables[iMo + PV->RHO] * der_u * sqrtKappa * LBCS);
          L2 = k_relax * F1B3 * (formerVariables[iMo + PV->RHO] - rho_b);
          L3 = L2;
          L4 = L2;
        } else {
          TERMM(1, "Direction can only be 0 or 1 in D2Q9!");
        }

        // compute u1,u2 and rho using the characteristic wave amplitudes
        u1 = formerVariables[iMo + PV->U]
             - deltaT * F1B2 / formerVariables[iMo + PV->RHO] * F1BCS / sqrtKappa * (L4 - L1);
        u2 = formerVariables[iMo + PV->V] - deltaT * L2;
        rho = formerVariables[iMo + PV->RHO] - deltaT * F1B2 / kappa * F1BCSsq * (L4 + L1)
              - F1BCSsq * L3 / kappa; // where 1.2 is kappa²

        // symmetric equilibriumfunctions rho_0=1 for the boundary
        eqD[0] = F1B9 * (rho + F9B2 * ((-u1) * (-u1) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[1] = F1B9 * (rho + F9B2 * ((u1) * (u1)-F1B3 * (u1 * u1 + u2 * u2)));
        eqD[2] = F1B9 * (rho + F9B2 * ((-u2) * (-u2) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[3] = F1B9 * (rho + F9B2 * ((u2) * (u2)-F1B3 * (u1 * u1 + u2 * u2)));
        eqD[4] = F1B36 * (rho + F9B2 * ((u1 + u2) * (u1 + u2) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[5] = F1B36 * (rho + F9B2 * ((u1 - u2) * (u1 - u2) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[6] = F1B36 * (rho + F9B2 * ((-u1 - u2) * (-u1 - u2) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[7] = F1B36 * (rho + F9B2 * ((-u1 + u2) * (-u1 + u2) - F1B3 * (u1 * u1 + u2 * u2)));
        eqD[8] = 0;

        u1c = m_solver->a_variable(currentID, PV->U);
        u2c = m_solver->a_variable(currentID, PV->V);
        rhoc = m_solver->a_variable(currentID, PV->RHO);

        // symmetric equilibriumfunctions rho_0=1 for the boundary
        eqDc[0] = F1B9 * (rhoc + F9B2 * ((-u1c) * (-u1c) - F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[1] = F1B9 * (rhoc + F9B2 * ((u1c) * (u1c)-F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[2] = F1B9 * (rhoc + F9B2 * ((-u2c) * (-u2c) - F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[3] = F1B9 * (rhoc + F9B2 * ((u2c) * (u2c)-F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[4] = F1B36 * (rhoc + F9B2 * ((u1c + u2c) * (u1c + u2c) - F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[5] = F1B36 * (rhoc + F9B2 * ((u1c - u2c) * (u1c - u2c) - F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[6] = F1B36 * (rhoc + F9B2 * ((-u1c - u2c) * (-u1c - u2c) - F1B3 * (u1c * u1c + u2c * u2c)));
        eqDc[7] = F1B36 * (rhoc + F9B2 * ((-u1c + u2c) * (-u1c + u2c) - F1B3 * (u1c * u1c + u2c * u2c)));

        // symmetric distribution function
        // Calculate average for opposite directions
        for(MInt j = 0; j < 8; j++) {
          meanDist[j] = F1B2
                        * (m_solver->a_oldDistribution(currentID, j)
                           + m_solver->a_oldDistribution(currentID, Ld::oppositeDist(j)));
        }

        for(MInt j = 0; j < 8; j++) {
          if(m_solver->a_hasNeighbor(currentID, j) == 0)
            m_solver->a_oldDistribution(currentID, Ld::oppositeDist(j)) =
                -m_solver->a_distribution(currentID, j) + 2 * eqD[j] + (2 - m_omega) * (meanDist[j] - eqDc[j]);
        }

        // only for bounce back in corners!!! (upper and lower boundary need to be walls)
        if(m_solver->a_hasNeighbor(currentID, 3) == 0) {
          m_solver->a_oldDistribution(currentID, 6) = m_solver->a_distribution(currentID, 4);
          m_solver->a_oldDistribution(currentID, 2) = m_solver->a_distribution(currentID, 3);
          m_solver->a_oldDistribution(currentID, 5) = m_solver->a_distribution(currentID, 7);
        }

        if(m_solver->a_hasNeighbor(currentID, 2) == 0) {
          m_solver->a_oldDistribution(currentID, 4) = m_solver->a_distribution(currentID, 6);
          m_solver->a_oldDistribution(currentID, 3) = m_solver->a_distribution(currentID, 2);
          m_solver->a_oldDistribution(currentID, 7) = m_solver->a_distribution(currentID, 5);
        }

        // set new boundary variables
        formerVariables[iMo + 0] = u1;
        formerVariables[iMo + 1] = u2;
        formerVariables[iMo + 2] = rho;
      }
    }
  }
}

/** \LBBC{sec_LBBC_bc40046, bc40046, 4046}
 *  \author Miro Gondrum
 *
 * Non-reflecting BC based on LODI - setting equilibrium<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::bc40046(MInt index)
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbBndCndDxQy<nDim, nDist, SysEqn>::bc40046(MInt index) {
  TRACE();

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    const MInt currentId = m_bndCells[i].m_cellId;
    if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->a_level(currentId)) != 0
       || m_solver->a_isHalo(currentId)) {
      continue;
    }
    //--LODI calculation--------------------------------------------------------
    MFloat rho_b, u_b[nDim];
    calcCharValues<2>(index, currentId, rho_b, &u_b[0]);
    //--set equilibrium state---------------------------------------------------
#ifdef WAR_NVHPC_PSTL
    MFloat oldDist[nDist] = {F0};
    for(MInt d = 0; d < nDist; d++) {
      oldDist[d] = m_solver->a_oldDistribution(currentId, d);
    }
    sysEqn().calcEqDists(rho_b, &u_b[0], oldDist);
    for(MInt d = 0; d < nDist; d++) {
      m_solver->a_oldDistribution(currentId, d) = oldDist[d];
    }
#else
    sysEqn().calcEqDists(rho_b, &u_b[0], &m_solver->a_oldDistribution(currentId, 0));
#endif
  }
}

/** \LBBC{sec_LBBC_charPressure2, charPressure2, 405x}
 *  in-/outflow bc
 *
 *  Non-reflecting characteristic bc according to Izquerdo et al. 2008 - combined with bc40030
 *  The density is set to 1.0 at the in- or outflow.<br>
 * LbBndCndDxQy<nDim, nDist, SysEqn>::charPressure2(MInt index, MInt direction)
 */
template <MInt nDim, MInt nDist, class SysEqn>
inline void LbBndCndDxQy<nDim, nDist, SysEqn>::charPressure2(MInt index, MInt direction) {
  TRACE();

  MFloat L1, L2, L3, L4, L5;
  MFloat rho, u[3] = {0.0, 0.0, 0.0}, rho_b;
  MFloat nghbr_rho, nghbr_u[3] = {0.0, 0.0, 0.0};

  MFloat deltaT = 1;
  MFloat sigma = 0.59;

  MFloat tmp, tmp2, b[6];
  MInt tmpDistId;

  MFloat kappa = 1.2;
  MFloat sqrtKappa = sqrt(1.2); // speed of sound c_s = sqrt(kappa*p/rho) where kappa = 1.2 is used at the boundary

  MInt currentID, nghbrID;

  MFloat k_relax; // relaxation factor

  MFloat der_u, der_p, der_v, der_w;

  MFloat rho_offset;

  MFloatScratchSpace nePart(nDist - 1, AT_, "nePart");
  MFloatScratchSpace eDistributions(nDist - 1, AT_, "eDistributions");

  // MInt ind = m_mapSegIdsInOutCnd[index];

  L1 = F0;
  L2 = F0;
  L3 = F0;
  L4 = F0;
  L5 = F0;

  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    currentID = m_bndCells[i].m_cellId;
    nghbrID = m_solver->c_neighborId(currentID, Ld::oppositeDist(direction));

    k_relax = sigma * (1 - m_Ma * m_Ma) * sqrtKappa * LBCS
              / (m_domainLength
                 * FFPOW2(m_solver->maxLevel()
                          - m_solver->a_level(currentID))); // m_domainLength represents the distance between the inlet
                                                            // and the outlet on the highest level of refinement

    m_omega =
        2.0 / (1.0 + 6.0 * m_solver->a_nu(currentID) * FFPOW2(m_solver->maxLevel() - m_solver->a_level(currentID)));

    if(m_densityFluctuations) {
      rho_b = 0.0;
      rho_offset = 1.0;
    } else {
      rho_b = 1.0;
      rho_offset = 0.0;
    }

    nghbr_rho = m_solver->a_variable(nghbrID, PV->RHO);
    nghbr_u[0] = m_solver->a_variable(nghbrID, PV->U);
    nghbr_u[1] = m_solver->a_variable(nghbrID, PV->V);
    nghbr_u[2] = m_solver->a_variable(nghbrID, PV->W);

    // determine non-equilibrium components in neighbor cell
    // -------------------------------------------------------
    tmp2 = nghbr_u[0] * nghbr_u[0] + nghbr_u[1] * nghbr_u[1] + nghbr_u[2] * nghbr_u[2];

    b[0] = -u[0];
    b[1] = u[0];
    b[2] = -u[1];
    b[3] = u[1];
    b[4] = -u[2];
    b[5] = u[2];

    // Calculation of eq-distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      eDistributions[j] = Ld::tp(1) * F1BCSsq * (nghbr_rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(2) * F1BCSsq * (nghbr_rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Calculation of eq-distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
      eDistributions[tmpDistId + j] =
          Ld::tp(3) * F1BCSsq * (nghbr_rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2);
    }
    // Determine non-eq parts
    for(MInt j = 0; j < nDist - 1; j++) {
      nePart[j] = m_solver->a_distribution(currentID, j) - eDistributions[j];
    }

    // Introduce characteristics:
    // ------------------------------------------
    // compute derivatives:
    der_p = CSsq * (m_solver->a_oldVariable(currentID, PV->RHO) - m_solver->a_oldVariable(nghbrID, PV->RHO));
    der_u = (m_solver->a_oldVariable(currentID, PV->U) - m_solver->a_oldVariable(nghbrID, PV->U));
    der_v = (m_solver->a_oldVariable(currentID, PV->V) - m_solver->a_oldVariable(nghbrID, PV->V));
    der_w = (m_solver->a_oldVariable(currentID, PV->W) - m_solver->a_oldVariable(nghbrID, PV->W));

    switch(direction) {
      case 0:
        // negative derivatives due to negative direction!
        L1 = (m_solver->a_oldVariable(currentID, PV->U) - sqrtKappa * LBCS)
             * (-der_p - (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * (-der_u) * sqrtKappa * LBCS);
        L2 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L3 = L2;
        L4 = L2;
        L5 = L2;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        u[1] = m_solver->a_oldVariable(currentID, PV->V) - deltaT * L3;
        u[2] = m_solver->a_oldVariable(currentID, PV->W) - deltaT * L4;
        break;

      case 1:
        L1 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L2 = 0;
        L3 = m_solver->a_oldVariable(currentID, PV->U) * der_v;
        L4 = m_solver->a_oldVariable(currentID, PV->U) * der_w;
        L5 = (m_solver->a_oldVariable(currentID, PV->U) + sqrtKappa * LBCS)
             * (der_p + (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * der_u * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        u[1] = m_solver->a_oldVariable(currentID, PV->V) - deltaT * L3;
        u[2] = m_solver->a_oldVariable(currentID, PV->W) - deltaT * L4;
        break;

      case 2:
        // negative derivatives due to negative direction!
        L1 = (m_solver->a_oldVariable(currentID, PV->V) - sqrtKappa * LBCS)
             * (-der_p - (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * (-der_v) * sqrtKappa * LBCS);
        L2 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L3 = L2;
        L4 = L2;
        L5 = L2;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U) - deltaT * L3;
        u[1] = m_solver->a_oldVariable(currentID, PV->V)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        u[2] = m_solver->a_oldVariable(currentID, PV->W) - deltaT * L4;
        break;

      case 3:
        L1 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L2 = 0;
        L3 = m_solver->a_oldVariable(currentID, PV->V) * der_u;
        L4 = m_solver->a_oldVariable(currentID, PV->V) * der_w;
        L5 = (m_solver->a_oldVariable(currentID, PV->V) + sqrtKappa * LBCS)
             * (der_p + (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * der_v * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U) - deltaT * L3;
        u[1] = m_solver->a_oldVariable(currentID, PV->V)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        u[2] = m_solver->a_oldVariable(currentID, PV->W) - deltaT * L4;
        break;

      case 4:
        // negative derivatives due to negative direction!
        L1 = (m_solver->a_oldVariable(currentID, PV->W) - sqrtKappa * LBCS)
             * (-der_p - (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * (-der_w) * sqrtKappa * LBCS);
        L2 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L3 = L2;
        L4 = L2;
        L5 = L2;

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U) - deltaT * L3;
        u[1] = m_solver->a_oldVariable(currentID, PV->V) - deltaT * L4;
        u[2] = m_solver->a_oldVariable(currentID, PV->W)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        break;

      case 5:
        L1 = k_relax * F1B3 * (m_solver->a_oldVariable(currentID, PV->RHO) - rho_b);
        L2 = 0;
        L3 = m_solver->a_oldVariable(currentID, PV->W) * der_u;
        L4 = m_solver->a_oldVariable(currentID, PV->W) * der_v;
        L5 = (m_solver->a_oldVariable(currentID, PV->W) + sqrtKappa * LBCS)
             * (der_p + (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * der_w * sqrtKappa * LBCS);

        // compute u1,u2,u3 using the characteristic wave amplitudes
        u[0] = m_solver->a_oldVariable(currentID, PV->U) - deltaT * L3;
        u[1] = m_solver->a_oldVariable(currentID, PV->V) - deltaT * L4;
        u[2] = m_solver->a_oldVariable(currentID, PV->W)
               - deltaT * F1B2 / (m_solver->a_oldVariable(currentID, PV->RHO) + rho_offset) * F1BCS / sqrtKappa
                     * (L5 - L1);
        break;

      default:
        TERMM(1, "Unrealistic case, please check your Boundary Conditions!!!");
    }

    // compute rho using the characteristic wave amplitudes
    rho = m_solver->a_oldVariable(currentID, PV->RHO) - deltaT * F1B2 / kappa * F1BCSsq * (L5 + L1)
          - F1BCSsq * L2 / kappa;

    // set new variables in boundary cell
    m_solver->a_variable(currentID, PV->U) = u[0];
    m_solver->a_variable(currentID, PV->V) = u[1];
    m_solver->a_variable(currentID, PV->W) = u[2];
    m_solver->a_variable(currentID, PV->RHO) = rho;

    // Calculation of missing incoming distributions in bnd cell
    //--------------------------------------------

    tmp2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    b[0] = -u[0];
    b[1] = u[0];
    b[2] = -u[1];
    b[3] = u[1];
    b[4] = -u[2];
    b[5] = u[2];

    // Calculation of distributions for directions with only one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      if(!m_solver->a_hasNeighbor(currentID, Ld::oppositeDist(j))) {
        m_solver->a_oldDistribution(currentID, j) =
            Ld::tp(1) * F1BCSsq * (rho * CSsq + b[j] + b[j] * b[j] * F1BCSsq * F1B2 - tmp2 * F1B2)
            + ((1 - m_omega) * nePart[j]);
      }
    }

    // Calculation of distributions for directions with two components
    tmpDistId = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      if(!m_solver->a_hasNeighbor(currentID, Ld::oppositeDist(tmpDistId + j))) {
        tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
        m_solver->a_oldDistribution(currentID, tmpDistId + j) =
            Ld::tp(2) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
            + ((1 - m_omega) * nePart[tmpDistId + j]);
      }
    }

    // Calculation of distributions for directions with three components
    tmpDistId = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      if(!m_solver->a_hasNeighbor(currentID, Ld::oppositeDist(tmpDistId + j))) {
        tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
        m_solver->a_oldDistribution(currentID, tmpDistId + j) =
            Ld::tp(3) * F1BCSsq * (rho * CSsq + tmp + tmp * tmp * F1BCSsq * F1B2 - tmp2 * F1B2)
            + ((1 - m_omega) * nePart[tmpDistId + j]);
      }
    }
  }
}

template <>
inline void LbBndCndDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::charPressure2(MInt /*index*/,
                                                                                      MInt /*direction*/) {
  TERMM(1, "This Boundary Condition is not implemented for 2D");
}
