// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbbndcnd.h"

#include <algorithm>
#include <set>
#include "COMM/mpioverride.h"
#include "GEOM/geometrycontext.h"
#include "IO/context.h"
#include "IO/infoout.h"
#include "lbconstants.h"
#include "lbgridboundarycell.h"
#include "lblatticedescriptor.h"
#include "lbsolver.h"
#include "property.h"

using namespace std;
using namespace lbconstants;

template <MInt nDim>
LbBndCnd<nDim>::LbBndCnd(LbSolver<nDim>* solver)
  : m_noInternalCells(solver->grid().noInternalCells()),
    m_solverId(solver->m_solverId),
    PV(solver->PV),
    m_noDistributions(solver->m_noDistributions),
    m_methodId(solver->m_methodId),
    m_Ma(solver->m_Ma),
    m_referenceLength(solver->m_referenceLength),
    m_domainLength(solver->m_domainLength),
    m_densityFluctuations(solver->m_densityFluctuations),
    m_pulsatileFrequency(solver->m_pulsatileFrequency) {
  TRACE();

  NEW_TIMER_GROUP(tgrp_BC, "Boundary Condition LB (solverId=" + std::to_string(m_solverId) + ")");
  NEW_TIMER(t_BCAll, "complete BC setup", tgrp_BC);
  m_t_BCAll = t_BCAll;
  RECORD_TIMER_START(t_BCAll);

  m_solver = solver;

  m_log << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl;
  m_log << "##                                             Boundary Conditions                                       "
           "           ##"
        << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl
        << endl;

  // external forcing terms
  m_Fext = solver->m_Fext;

  m_gridCutTest = "SAT";
  m_gridCutTest = Context::getSolverProperty<MString>("gridCutTest", m_solverId, AT_, &m_gridCutTest);

  /*! \page propertyPage1
    \section interpolationDistMethod
    <code>MString LbBndCnd::m_interpolationDistMethod</code>\n
    default = <code>perpOp</code>\n\n
    This property defines the way the distances to the STL are calculated.
    <ul>
    <li><code>sphere</code> analytically defined sphere</li>
    <li><code>perpOp</code> uses the perpendicular operator</li>
    <li><code>pipe</code> analystically defined cylinder</li>
    <li><code>STD</code> old standard method (replaced by perpOp)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_interpolationDistMethod = "perpOp";
  m_interpolationDistMethod =
      Context::getSolverProperty<MString>("interpolationDistMethod", m_solverId, AT_, &m_interpolationDistMethod);

  /*! \page propertyPage1
    \section outputWallDistanceField
    <code>MString LbBndCnd::m_outputWallDistanceField</code>\n
    default = <code>false</code>\n\n
    This property defines, whether the isFluid state and the wall distance fields are dumped in a file.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_outputWallDistanceField =
      Context::getSolverProperty<MBool>("outputWallDistanceField", m_solverId, AT_, &m_outputWallDistanceField);

  /*! \page propertyPage1
    \section multiBCTreatment
    <code>MString LbBndCnd::multiBCTreatment</code>\n
    default = <code>I-W-P</code>\n\n
    This property defines the way BCs are treated for cells being part of multiple BCs
    <ul>
    <li><code>I-W-P</code> order: in/oultet - wall - periodic</li>
    <li><code>I-P-W</code> order: in/outlet - periodic - wall</li>
    <li><code>W-P-I</code> order: wall - periodic - in/outlet</li>
    <li><code>W-I-P</code> order: wall - in/outlet - periodic</li>
    <li><code>P-W-I</code> order: periodic - wall - in/outlet</li>
    <li><code>P-I-W</code> order: periodic - in/outlet</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_multiBCTreatment = "I-W-P";
  m_multiBCTreatment = Context::getSolverProperty<MString>("multiBCTreatment", m_solverId, AT_, &m_multiBCTreatment);

  /*! \page propertyPage1
    \section lbNoMovingWalls
    <code>MInt LbBndCndDxQy::m_lbNoMovingWalls</code>\n
    default = <code>0</code>\n\n
    Number of walls with a wall velocity abs(u_w) > 0.
    <ul>
    <li><code>any integer >= 0</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_lbNoMovingWalls = 0;
  m_lbNoMovingWalls = Context::getSolverProperty<MInt>("lbNoMovingWalls", m_solverId, AT_, &m_lbNoMovingWalls);

  if(m_lbNoMovingWalls > 0) {
    mAlloc(m_segIdMovingWalls, m_lbNoMovingWalls, "m_segIdMovingWalls", 0, AT_);

    /*! \page propertyPage1
      \section segIdMovingWalls
      <code>MInt LbBndCndDxQy::m_segIdMovingWalls</code>\n
      default = <code>0</code>\n\n
      Segment id of walls with a wall velocity abs(u_w) > 0.
      <ul>
      <li><code>integer >= 0</code> (on)</li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    for(MInt i = 0; i < m_lbNoMovingWalls; i++) {
      m_segIdMovingWalls[i] = Context::getSolverProperty<MInt>("segIdMovingWalls", m_solverId, AT_, i);
    }

    mAlloc(m_lbWallVelocity, m_lbNoMovingWalls * nDim, "m_lbWallVelocity", F0, AT_);

    /*! \page propertyPage1
      \section lbWallVelocity
      <code>MFloat LbBndCndDxQy::m_lbWallVelocity</code>\n
      default = <code>0.0</code>\n\n
      Velocity components of each wall with a velocity abs(u_w) > 0.
      <ul>
      <li><code>any float </code> (on)</li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    for(MInt i = 0; i < m_lbNoMovingWalls; i++) {
      for(MInt d = 0; d < nDim; d++) {
        m_lbWallVelocity[i * nDim + d] =
            Context::getSolverProperty<MFloat>("lbWallVelocity", m_solverId, AT_, (i * nDim + d));
        m_lbWallVelocity[i * nDim + d] *= m_Ma * LBCS;
      }
    }
  }

  /*! \page propertyPage1
    \section calcWallForces
    <code>MBool LbBndCnd::m_calcWallForces</code>\n
    default = <code>false</code>\n\n
    Activates the calculation of wall forces
    <ul>
    <li><code>false</code> </li>
    <li><code>true</code> </li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_calcWallForces = false;
  m_calcWallForces = Context::getSolverProperty<MBool>("calculateWallForces", m_solverId, AT_, &m_calcWallForces);

  if(m_calcWallForces) {
    /*! \page propertyPage1
      \section calcWallForcesInterval
      <code>MInt LbBndCnd::m_calcWallForcesInterval</code>\n
      default = <code>1</code>\n\n
      Determine the update interval for wall force calculation
      <ul>
      <li><code>Any integer number > 0</code> </li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_calcWallForcesInterval = 1;
    m_calcWallForcesInterval =
        Context::getSolverProperty<MInt>("calculateWallForcesInterval", m_solverId, AT_, &m_calcWallForcesInterval);

    /*! \page propertyPage1
      \section forceFileName
      <code>MInt LbBndCnd::m_forceFile</code>\n
      default = <code>1</code>\n\n
      Determine the output file name for wall force calculation
      <ul>
      <li><code>Any valid file name</code> </li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_forceFile = "forces.log";
    m_forceFile = Context::getSolverProperty<MString>("forceFileName", m_solverId, AT_, &m_forceFile);
  }

  /*! \page propertyPage1
    \section calcBcResidual
    <code>MInt LbBndCnd::m_calcBcResidual</code>\n
    default = <code>false</code>\n\n
    Determine whether residuals for BC are written out
    <ul>
    <li><code>false</code> off </li>
    <li><code>true</code> on</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_calcBcResidual = false;
  m_calcBcResidual = Context::getSolverProperty<MBool>("calcBcResidual", m_solverId, AT_, &m_calcBcResidual);

  /*! \page propertyPage1
    \section lbNoHeatedWalls
    <code>MInt LbBndCndDxQy::m_lbNoHeatedWalls</code>\n
    default = <code>0</code>\n\n
    Number of walls with a wall temperature abs(T) > 0.
    <ul>
    <li><code>any integer >= 0</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_lbNoHeatedWalls = 0;
  m_lbNoHeatedWalls = Context::getSolverProperty<MInt>("lbNoHeatedWalls", m_solverId, AT_, &m_lbNoHeatedWalls);

  if(m_lbNoHeatedWalls > 0) {
    mAlloc(m_segIdHeatedWalls, m_lbNoHeatedWalls, "m_segIdHeatedWalls", 0, AT_);

    /*! \page propertyPage1
      \section segIdHeatedWalls
      <code>MInt LbBndCndDxQy::m_segIdHeatedWalls</code>\n
      default = <code>0</code>\n\n
      Segment id of walls with a wall temperature abs(T) > 0.
      <ul>
      <li><code>integer >= 0</code> (on)</li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    for(MInt i = 0; i < m_lbNoHeatedWalls; i++) {
      m_segIdHeatedWalls[i] = Context::getSolverProperty<MInt>("segIdHeatedWalls", m_solverId, AT_, i);
    }

    mAlloc(m_lbWallTemperature, m_lbNoHeatedWalls, "m_lbWallTemperature", F0, AT_);

    /*! \page propertyPage1
      \section lbWallTemperature
      <code>MFloat LbBndCndDxQy::m_lbWallTemperature</code>\n
      default = <code>0.0</code>\n\n
      Velocity components of each wall with a temperature abs(T) > 0.
      <ul>
      <li><code>any float </code> (on)</li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    for(MInt i = 0; i < m_lbNoHeatedWalls; i++) {
      m_lbWallTemperature[i] = Context::getSolverProperty<MFloat>("lbWallTemperature", m_solverId, AT_, i);
    }
  }

  // the number of segments that carry periodic conditions
  m_noPeriodicSegments = 0;
  m_periodicSegmentsIds = nullptr;
  if(m_solver->m_geometry->geometryContext().propertyExists("periodicSegmentsIds", 0)) {
    m_noPeriodicSegments = m_solver->m_geometry->geometryContext().getProperty("periodicSegmentsIds", 0)->count();
    // the segment ids that carry periodic conditions
    mAlloc(m_periodicSegmentsIds, m_noPeriodicSegments, "m_periodicSegmentsIds", 0, AT_);
    // m_periodicSegmentsIds = new MInt[m_noPeriodicSegments];
    for(MInt i = 0; i < m_noPeriodicSegments; i++)
      m_periodicSegmentsIds[i] =
          *m_solver->m_geometry->geometryContext().getProperty("periodicSegmentsIds", 0)->asInt(i);
  }


  if(m_solver->grid().isActive()) {
    // This initializes the initial velocity vectors and the normal vectors for each boundary segment
    initMembers();

    createBoundaryCells();
    sortBoundaryCells();

    calculateVectors();
    processAngles();
    printBndVelInfo();
  }

  // for certain BCs it is important to have a neighbor communicator
  if(m_solver->grid().isActive()) {
    setBCNeighborCommunicator();
    setBCWallNeighborCommunicator();
  }

  /*! \page propertyPage1
    \section latentHeat
    <code>MBool LbBndCnd::latentHeat</code>\n
    default = <code>1</code>\n\n
    This property activates latentHeat calculation
    <ul>
    <li><code>0</code> off </li>
    <li><code>1</code> on</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, THERMAL, TRANSPORT</i>
  */
  m_latentHeat = false;
  m_latentHeat = Context::getSolverProperty<MBool>("latentHeat", m_solverId, AT_, &m_latentHeat);

  if(m_latentHeat) {
    /*! \page propertyPage1
      \section calcSubLayerDist
      <code>MBool LbBndCnd::calcSublayerDist</code>\n
      default = <code>1</code>\n\n
      This property activates the calculation of wall distances to a second wall,
      not given by an STL, along the direction of the PPDFs. This is required if
      sublayers outside of the STL are used. For example, it is used for the
      calculation of heat conduction in the solid domain.
      <ul>
      <li><code>false</code> off </li>
      <li><code>true</code> on</li>
      </ul>\n
      Keywords: <i>LATTICE_BOLTZMANN, THERMAL, TRANSPORT</i>
    */
    m_calcSublayerDist = false;
    m_calcSublayerDist = Context::getSolverProperty<MBool>("calcSublayerDist", m_solverId, AT_, &m_calcSublayerDist);
  }

  if(m_solver->isActive()) {
    setBndCndHandler();
  }

  RECORD_TIMER_STOP(t_BCAll);

  if(m_solver->m_geometry->m_debugParGeom) {
    DISPLAY_ALL_GROUP_TIMERS(tgrp_BC);
    m_log << endl;
  }

  m_log << endl;
}


//! The destructor
/*!
 */
template <MInt nDim>
LbBndCnd<nDim>::~LbBndCnd() {
  TRACE();

  if(m_numberOfCommBCs > 0) {
    mDeallocate(m_allDomainsHaveBC);
    mDeallocate(m_noBCNeighbors);
    mDeallocate(m_totalNoBcCells);
    mDeallocate(tmp_group);
    mDeallocate(BCGroup);
    mDeallocate(m_BCComm);
    if(m_calcBcResidual) mDeallocate(m_BCResidualStream);
    mDeallocate(m_BCOutputFileName);
    if((MInt)(m_bndCndIds.size()) > 0) mDeallocate(m_mapBndCndIdSegId);
  }

  mDeallocate(m_inOutSegmentsIds);
  mDeallocate(m_bndNormals);
  mDeallocate(m_initialVelocityVecs);
  mDeallocate(m_bndNormalDirection);
  mDeallocate(m_exDirs);
  mDeallocate(m_exWeights);
  mDeallocate(m_phi);
  mDeallocate(m_theta);

  // Delete boundary condition handlers
  bndCndHandlerVariables.clear();
  bndCndHandlerRHS.clear();
  bndCndInitiator.clear();

  // Delete boundary cells
  m_bndCells.clear();
  m_boundaryCellsMb.clear();
}

template <MInt nDim>
void LbBndCnd<nDim>::processAngles() {
  TRACE();

  NEW_SUB_TIMER(t_processAngles, "process angles", m_t_BCAll);
  RECORD_TIMER_START(t_processAngles);


  m_log << "  + Processing angles..." << endl;

  // -------------------------------------------------------------------------------
  // Determine data for arbitrary outlet directions
  // (polar angles & directions are derived from the normal vector of the outlet)
  // -------------------------------------------------------------------------------

  MFloat projection, angle1, angle2;
  MFloat length;

  MIntScratchSpace axesDirs(nDim, AT_, "axesDirs");
  MFloat projections[3];

  for(MInt id = 0; id < m_noInOutSegments; id++) {
    MInt seg_id = -1;
    for(MInt o = 0; o < (MInt)m_bndCndSegIds.size(); o++)
      if(m_inOutSegmentsIds[id] == m_bndCndSegIds[o]) {
        seg_id = o;
        break;
      }

    if(seg_id < 0) continue;

    if(m_bndNormals[id][0] > 0) {
      m_phi[id] = atan(m_bndNormals[id][1] / m_bndNormals[id][0]);

      axesDirs[0] = 0;

      if(m_bndNormals[id][1] > 0)
        axesDirs[1] = 2;
      else
        axesDirs[1] = 3;
    } else {
      axesDirs[0] = 1;

      if(m_bndNormals[id][1] > 0) {
        if(approx(m_bndNormals[id][0], 0.0, MFloatEps)) {
          m_phi[id] = PI / 2;
        } else {
          m_phi[id] = PI + atan(m_bndNormals[id][1] / m_bndNormals[id][0]);
        }

        axesDirs[1] = 2;

      } else {
        if(approx(m_bndNormals[id][0], 0.0, MFloatEps)) {
          m_phi[id] = -PI / 2;
        } else {
          m_phi[id] = -PI + atan(m_bndNormals[id][1] / m_bndNormals[id][0]);
        }
        axesDirs[1] = 3;
      }
    }

    IF_CONSTEXPR(nDim == 3) {
      m_theta[id] = asin(m_bndNormals[id][2]);

      if(m_bndNormals[id][2] > 0)
        axesDirs[2] = 4;
      else
        axesDirs[2] = 5;
    }
    else m_theta[id] = 0;

    for(MInt i = 0; i < nDim; i++) {
      projections[i] = 0.0;
      for(MInt j = 0; j < nDim; j++)
        projections[i] += (LbLatticeDescriptorBase<3>::idFld(axesDirs[i], j) - 1.0)
                          * m_bndNormals[id][j]; // TODO labels:LB dxqy: 3 replaceable by nDim ?
    }

    IF_CONSTEXPR(nDim == 3) {
      m_log << "    - InOutSegment " << id << ":" << endl;
      m_log << "      * segment id:                                          " << m_bndCndSegIds[seg_id] << endl;
      m_log << "      * distributions along axes which point into the fluid: " << axesDirs[0] << " " << axesDirs[1]
            << " " << axesDirs[0] << endl;
      m_log << "      * scalar product with bndNormal:                       " << projections[0] << " "
            << projections[1] << " " << projections[2] << endl;
    }
    else {
      m_log << "    - InOutSegment " << id << ":" << endl;
      m_log << "      * angle phi:                                           " << 180 * m_phi[id] / PI << endl;
      m_log << "      * distributions along axes which point into the fluid: " << axesDirs[0] << " " << axesDirs[1]
            << endl;
      m_log << "      * scalar product with bndNormal:                       " << projections[0] << " "
            << projections[1] << endl;
    }


    // -----------------------------------------------------------------------------------------------------------
    // Determine directions and weights for extrapolation, i.e. go through all
    // directions and take those two which have the highest projection on the inverse normal vector.
    // The weights are calculated from the ratio of enclosed angles.
    // -----------------------------------------------------------------------------------------------------------

    angle1 = 0;
    angle2 = 0;
    m_exDirs[id][0] = -1;
    m_exDirs[id][1] = -1;

    // number of neighbors was raised to maximum to improve extrapolation
    for(MInt i = 0; i < pow(3.0, nDim) - 1; i++) {
      //	for (MInt i=0; i < m_noDistributions - 1; i++){
      projection = 0.0;
      length = 0.0;
      for(MInt j = 0; j < nDim; j++) {
        projection += -1.0 * (LbLatticeDescriptorBase<3>::idFld(i, j) - 1.0)
                      * m_bndNormals[id][j]; // TODO labels:LB dxqy: 3 replaceable by nDim ?
        length += (LbLatticeDescriptorBase<3>::idFld(i, j) - 1.0)
                  * (LbLatticeDescriptorBase<3>::idFld(i, j) - 1.0); // TODO labels:LB dxqy: 3 replaceable by nDim ?
      }
      length = sqrt(length);
      projection = projection / length;

      if(projection > angle1) {
        if(angle1 > angle2) {
          angle2 = angle1;
          m_exDirs[id][1] = m_exDirs[id][0];
        }
        m_exDirs[id][0] = i;
        angle1 = projection;
      } else {
        if(projection > angle2) {
          m_exDirs[id][1] = i;
          angle2 = projection;
        }
      }
    }

    angle1 = abs(acos(angle1));
    angle2 = abs(acos(angle2));

    m_exWeights[id][0] = angle2 / (angle1 + angle2);
    m_exWeights[id][1] = angle1 / (angle1 + angle2);

    m_log << "      * extrapolation directions:                            " << m_exDirs[id][0] << " "
          << m_exDirs[id][1] << endl;
    m_log << "      * angle enclosed with bndNormal:                       " << 180 * angle1 / PI << " "
          << 180 * angle2 / PI << endl;
    m_log << "      * extrapolation weights:                               " << m_exWeights[id][0] << " "
          << m_exWeights[id][1] << endl;

    // Oha, there seems to be something wrong, no extrapolation neighbors found
    if(m_exDirs[id][0] == -1 && m_exDirs[id][1] == -1) {
      stringstream errorMsg;
      errorMsg << "ERROR: no neighbors found for extrapolation for segment id " << m_inOutSegmentsIds[id] << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }
  }

  m_log << endl;

  RECORD_TIMER_STOP(t_processAngles);
}

/** \brief This function prints the gathered information on the boundary vectors and the initial velocity vector into
 *the log file. \author Andreas Lintermann \date 28.01.2010, 27.04.2015
 *
 * Prints all information about the boundary conditions.
 **/
template <MInt nDim>
void LbBndCnd<nDim>::printBndVelInfo() {
  TRACE();

  NEW_SUB_TIMER(t_printBndVelInfo, "print boundary velcoity info", m_t_BCAll);
  RECORD_TIMER_START(t_printBndVelInfo);


  m_log << "  + We have the following inflow / outflow boundary information: " << endl;

  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
    MInt mapId = m_mapSegIdsInOutCnd[i];
    if(mapId != -1) {
      m_log << "    - information for segment " << m_bndCndSegIds[i] << endl;
      m_log << "      * BC:         " << m_bndCndIds[i] << endl;
      m_log << "      * no cells:   " << m_noBndCellsPerSegment[m_bndCndSegIds[i]] << endl;
      m_log << "      * normal id:  " << m_mapSegIdsInOutCnd[i] << endl;
      m_log << "      * normal dir: " << m_bndNormalDirection[mapId] << " ("
            << (m_bndNormalDirection[mapId] == -1 ? "outside" : "inside") << ")" << endl;
      m_log << "      * normal:     ";
      for(MInt j = 0; j < nDim; j++)
        m_log << m_bndNormals[mapId][j] << " ";
      m_log << endl;
      m_log << "      * init vel:   ";
      for(MInt j = 0; j < nDim; j++)
        m_log << m_initialVelocityVecs[mapId][j] << " ";
      m_log << endl;
      m_log << "      * ext dirs:   " << m_exDirs[mapId][0] << " " << m_exDirs[mapId][1] << endl;
      m_log << "      * ext wght:   " << m_exWeights[mapId][0] << " " << m_exWeights[mapId][1] << endl;
    }
  }

  m_log << endl;

  RECORD_TIMER_STOP(t_printBndVelInfo);
}

/** \brief This function initializes member vectors.
 *
 * \author Andreas Lintermann
 * \date 17.07.2015
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::initMembers() {
  TRACE();

  NEW_SUB_TIMER(t_initMembers, "init members", m_t_BCAll);
  RECORD_TIMER_START(t_initMembers);

  m_log << "  + Initializing some members..." << endl;

  // the number of all segments
  m_noSegments = m_solver->m_geometry->geometryContext().getNoSegments();

  if(m_solver->m_geometry->geometryContext().propertyExists("inOutSegmentsIds", 0)) {
    // the number of segments that carry in-/outflow conditions
    m_noInOutSegments = m_solver->m_geometry->geometryContext().getProperty("inOutSegmentsIds", 0)->count();

    // the segment ids that carry in-/outflow conditions
    mAlloc(m_inOutSegmentsIds, m_noInOutSegments, "m_inOutSegmentsIds", 0, AT_);
  } else
    m_noInOutSegments = 0;

  // init the velocity vectors and the segment normals, the direction array and the array holding the ids of the
  // in-/outflow conditions
  if(m_noInOutSegments > 0) {
    mAlloc(m_bndNormals, m_noInOutSegments, nDim, "m_bndNormals", F0, AT_);
    mAlloc(m_initialVelocityVecs, m_noInOutSegments, nDim, "m_initialVelocityVecs", F0, AT_);
    mAlloc(m_bndNormalDirection, m_noInOutSegments, "m_bndNormalDirection", -1, AT_);

    mAlloc(m_exDirs, m_noInOutSegments, 2, "m_exDirs", 0, AT_);
    mAlloc(m_exWeights, m_noInOutSegments, 2, "m_exWeights", F0, AT_);
    mAlloc(m_phi, m_noInOutSegments, "m_phi", F0, AT_);
    mAlloc(m_theta, m_noInOutSegments, "m_theta", F0, AT_);
  }
  for(MInt i = 0; i < m_noInOutSegments; i++) {
    m_inOutSegmentsIds[i] = *m_solver->m_geometry->geometryContext().getProperty("inOutSegmentsIds", 0)->asInt(i);
  }

  /*! \page propertyPage1
    \section bndNormalMethod
    <code>MString LbBndCnd::bndNormalMethod</code>\n
    default = <code>calcNormal</code>\n\n
    This property defines the way normals for in- and outlets are treated.
    <ul>
    <li><code>read</code> reads the normals from the property file (bndNormalVectors)</li>
    <li><code>calcNormal</code> calculates the normals based on triangle information and averages them</li>
    <li><code>fromSTL</code> reads the normals from STL and averages them</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_bndNormalMethod = "calcNormal";
  m_bndNormalMethod = Context::getSolverProperty<MString>("bndNormalMethod", m_solverId, AT_, &m_bndNormalMethod);

  /*! \page propertyPage1
    \section initVelocityMethod
    <code>MString LbBndCnd::initVelocityMethod</code>\n
    default = <code>calcNormal</code>\n\n
    This property defines the way initial velocity vectors for in- and outlets are treated. In general, if the
    nmormals point outward, such that the according velocity vectors should point inside.
    <ul>
    <li><code>read</code> reads the vectors from the property file (initialVelocityVectors)</li>
    <li><code>calcNormal</code> vectors are set to bndNormalVectors</li>
    <li><code>fromSTL</code> vectors are set to -1*bndNormalVectors</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_initVelocityMethod = "calcNormal";
  m_initVelocityMethod =
      Context::getSolverProperty<MString>("initVelocityMethod", m_solverId, AT_, &m_initVelocityMethod);

  /*! \page propertyPage1
    \section fastParallelGeomNormals
    <code>MInt LbBndCnd::fastParallelGeomNormals</code>\n
    default = <code>0</code>\n\n
    This property defines how normals are calculated for a parallel geometry.
    <ul>
    <li><code>0</code> normals are calucated for each processor participating in a BC and then the average is
    calculated; involves inside/outsie detection on all processes.</li> <li><code>1</code> normals are calculates on
    only one process in the BC communicator, only one inside-outside detetction is performed.</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_fastParallelGeomNormals = 0;
  m_fastParallelGeomNormals =
      Context::getSolverProperty<MInt>("fastParallelGeomNormals", m_solverId, AT_, &m_fastParallelGeomNormals);

  m_log << endl;

  RECORD_TIMER_STOP(t_initMembers);
}

/** \brief This function calculates member vectors.
 *
 * \author Andreas Lintermann
 * \date 17.07.2015
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::calculateVectors() {
  TRACE();

  // Vector calculation is only performed for in/out segments.
  if(m_noInOutSegments < 1) return;

  NEW_SUB_TIMER(t_calcVecs, "calculate vectors", m_t_BCAll);
  g_tc_geometry.push_back(pair<MString, MInt>("calculate vectors", t_calcVecs));
  RECORD_TIMER_START(t_calcVecs);

  m_log << "  + Calculating vectors..." << endl;

  // read from file
  if(m_bndNormalMethod == "read") {
    m_log << "    - reading normals in 3D (bndNormalMethod is set to '" << m_bndNormalMethod << "')" << endl;
    for(MInt i = 0, k = 0; i < m_noInOutSegments; i++) {
      MBool own = false;
      for(MInt o = 0; o < (MInt)m_bndCndSegIds.size(); o++) {
        if(m_inOutSegmentsIds[i] == m_bndCndSegIds[o]) { // m_inOutSegmentsIds is from geometry-file
          own = true;
          break;
        }
      }

      for(MInt j = 0; j < nDim; j++) {
        /*! \page propertyPage1
         \section bndNormalVectors
         <code>MFloat** LbBndCnd::bndNormalVecs</code>\n
         This property defines the boundary normal vectors per in/outlet BC in the format x,y,z.\n
         Keywords: <i>LATTICE_BOLTZMANN</i>
        */
        m_bndNormals[i][j] = Context::getSolverProperty<MFloat>("bndNormalVectors", m_solverId, AT_, k);

        // normal vector points out of the body
        m_bndNormalDirection[i] = -1;

        k++;
      }

      if(own) {
        m_log << "      * segment " << m_inOutSegmentsIds[i] << " has normal: (" << m_bndNormals[i][0] << " "
              << m_bndNormals[i][1];
        IF_CONSTEXPR(nDim > 2)
        m_log << " " << m_bndNormals[i][2] << ")" << endl;
        else m_log << ")" << endl;
      }
    }
  }
  // Get the normals from the STL
  else {
    if(m_bndNormalMethod == "calcNormal")
      retrieveNormal = &LbBndCnd::calculateNormalFromTriangle;
    else if(m_bndNormalMethod == "fromSTL")
      retrieveNormal = &LbBndCnd::getNormalFromSTL;

    calculateBndNormals();
  }

  // fill initial velocity vectors
  for(MInt i = 0, k = 0; i < m_noInOutSegments; i++) {
    for(MInt j = 0; j < nDim; j++) {
      // read from file
      if(m_initVelocityMethod == "read")
        /*! \page propertyPage1
          \section initialVelocityVectors
          <code>MFloat** LbBndCnd::initialVelocityVecs</code>\n
          This property defines the velocity vectors per in/outlet BC in the format x,y,z.\n
          Keywords: <i>LATTICE_BOLTZMANN</i>
        */
        m_initialVelocityVecs[i][j] = Context::getSolverProperty<MFloat>("initialVelocityVectors", m_solverId, AT_, k);

      // just use the segment normals
      else if(m_initVelocityMethod == "calcNormal")
        m_initialVelocityVecs[i][j] = m_bndNormals[i][j];

      // use inverted segment normals
      else if(m_initVelocityMethod == "fromSTL")
        m_initialVelocityVecs[i][j] = -1 * m_bndNormals[i][j];

      k++;
    }

    // Normalize velocity vectors
    MFloat normalLength = 0.0;
    for(MInt j = 0; j < nDim; j++) {
      normalLength += m_initialVelocityVecs[i][j] * m_initialVelocityVecs[i][j];
    }
    normalLength = sqrt(normalLength);
    for(MInt j = 0; j < nDim; j++) {
      m_initialVelocityVecs[i][j] = m_initialVelocityVecs[i][j] / normalLength;
    }
  }

  m_log << endl;

  RECORD_TIMER_STOP(t_calcVecs);
}

/** \brief Calculate the normal based on a triangle
 *
 * \author Andreas Lintermann
 * \date 27.04.2015
 *
 * Calculate the normal from a traingle using the cross product. Just in case the normal has a length of 0
 * false is returned, otherwise true. If invalid, the triangle is skipped for the calculation of
 * the segment normal in calculateBndNormals().
 *
 * \param[in] ge the GeometryElement (a triangle)
 * \param[in] normal address of the normal
 * \return the validity of the triangle normal
 *
 **/
template <MInt nDim>
inline MBool LbBndCnd<nDim>::calculateNormalFromTriangle(GeometryElement<nDim> ge, MFloat* normal) {
  TRACE();

  std::array<MFloat, nDim> edge1;
  std::array<MFloat, nDim> edge2;

  // the edges
  for(MInt d = 0; d < nDim; d++) {
    edge1[d] = ge.m_vertices[1][d] - ge.m_vertices[0][d];
    edge2[d] = ge.m_vertices[2][d] - ge.m_vertices[0][d];
  }

  MFloat normal_len = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    MInt next = (d + 1) % nDim;
    MInt nnext = (d + 2) % nDim;
    normal[d] = edge1[next] * edge2[nnext] - edge1[nnext] * edge2[next];
    normal_len += normal[d] * normal[d];
  }

  normal_len = sqrt(normal_len);

  for(MInt d = 0; d < nDim; d++)
    normal[d] /= normal_len;

  if(approx(normal_len, 0.0, MFloatEps)) {
    return false;
  } else {
    return true;
  }
}

/** \brief Get the normal from the STL
 *
 * \author Andreas Lintermann
 * \date 27.04.2015
 *
 * Gets the normal directly from the STL information. Just in case the normal has a length of 0
 * false is returned, otherwise true. If invalid, the triangle is skipped for the calculation of
 * the segment normal in calculateBndNormals().
 *
 * \param[in] ge the GeometryElement (a triangle)
 * \param[in] normal address of the normal
 *
 * \return the validity of the triangle normal
 *
 **/
template <MInt nDim>
inline MBool LbBndCnd<nDim>::getNormalFromSTL(GeometryElement<nDim> ge, MFloat* normal) {
  TRACE();

  MFloat normal_len = 0.0;
  for(MInt d = 0; d < nDim; d++) {
    normal[d] = ge.m_normal[d];
    normal_len += normal[d] * normal[d];
  }

  normal_len = sqrt(normal_len);

  for(MInt d = 0; d < nDim; d++)
    normal[d] /= normal_len;
  if(approx(normal_len, 0.0, MFloatEps))
    return false;
  else
    return true;
}

/** returns the center, the normal and the count of a segment
 *
 * \author Andreas Lintermann
 * \date 14.01.2016
 *
 * \param[in] segmentId the id of the segment to test
 * \param[in] normal a MFloatScratchSpace to be filled with the avaregaed normal
 * \param[in] center a MFloatScratchSpace to be filled with the center
 * \param[in] count pointer to a counter
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::calculateAveragedNormalAndCenter(MInt segmentId, MFloat* const normal, MFloat* const center) {
  TRACE();

  for(MInt d = 0; d < nDim; d++) {
    normal[d] = 0.0;
    center[d] = 0.0;
  }

  MInt offStart = 0;
  MInt offEnd = 0;

  if(m_solver->m_geometry->m_parallelGeometry) {
    offStart = m_solver->m_geometry->m_segmentOffsets[segmentId];
    offEnd = m_solver->m_geometry->m_segmentOffsets[segmentId + 1];
  } else {
    offStart = m_solver->m_geometry->m_segmentOffsetsWithoutMB[segmentId];
    offEnd = m_solver->m_geometry->m_segmentOffsetsWithoutMB[segmentId + 1];
  }
  MInt count = 0;

  for(MInt e = offStart; e < offEnd; e++) {
    GeometryElement<nDim> triangle = m_solver->m_geometry->elements[e];
    std::array<MFloat, nDim> normal_t;

    // this retrieves the normal either from the STL or recalculates it from the triangle
    MBool is_valid = (this->*retrieveNormal)(triangle, normal_t.data());

    // make sure that this is a valid triangle, otherwise, skip
    if(is_valid) {
      std::array<MFloat, nDim> tmp;
      for(MInt d = 0; d < nDim; d++)
        tmp[d] = normal[d] + normal_t[d];

      // this is the dangerous case in which we might go back directly onto the surface again
      MBool flip = true;
      for(MInt d = 0; d < nDim; d++)
        flip = flip && approx(tmp[d], 0.0, MFloatEps);

      if(flip)
        for(MInt d = 0; d < nDim; d++)
          normal_t[d] *= -1;

      for(MInt d = 0; d < nDim; d++) {
        normal[d] += normal_t[d];
        for(MInt v = 0; v < nDim; v++)
          center[d] += triangle.m_vertices[v][d];
      }
      count++;
    }
  }

  // Oha, there seems to be something wrong with the geometry, no normal found!
  if(count == 0) {
    stringstream errorMsg;
    errorMsg << "ERROR: no normal found for segment id " << segmentId << endl;
    m_log << errorMsg.str();
    TERMM(1, errorMsg.str());
  }

  // normalize
  MFloat avg_normal_len = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

  for(MInt d = 0; d < nDim; d++) {
    normal[d] /= avg_normal_len;
    center[d] /= nDim * count;
  }
}

/** Caluates the normal for one participatiung process for an in/outlet BC and distributes
 *
 * \author Andreas Lintermann
 * \date 12.06.2016
 *
 * Only one process in the BC communicator performs a normal calculation and uses this to generate a point
 * in the center of gravity. It is then checked if this point is inside/outside the geometry to define the
 * boundary normal.
 *
 * param[in] own_segments is a vector of tuples of positions and segmentId
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::fastParallelGeomNormals3D(vector<pair<MInt, MInt>> own_segments) {
  TRACE();

  MIntScratchSpace numAllOwners(m_solver->noDomains(), AT_, "numAllOwners");
  for(MInt d = 0; d < m_solver->noDomains(); d++)
    numAllOwners[d] = 0;
  numAllOwners[m_solver->domainId()] = own_segments.size();

  MPI_Allreduce(MPI_IN_PLACE, numAllOwners.getPointer(), m_solver->noDomains(), MPI_INT, MPI_SUM, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "numAllOwners.getPointer()");

  MInt sumOwners = 0;
  MInt myOffsetStart = 0;
  for(MInt d = 0; d < m_solver->noDomains(); d++) {
    if(d == m_solver->domainId()) myOffsetStart = sumOwners;
    sumOwners += numAllOwners[d];
  }

  // calculate my local normal and the test point (used later) and print to log
  MFloatScratchSpace test_pts(sumOwners * nDim, AT_, "test_pts");
  MIntScratchSpace io(sumOwners * nDim, AT_, "io");
  for(MInt i = 0; i < sumOwners * nDim; i++) {
    test_pts[i] = 0.0;
    io[i] = 0;
  }

  for(MInt s = 0; s < (MInt)own_segments.size(); s++) {
    MInt pos = own_segments[s].first;
    MInt segId = own_segments[s].second;
    MInt posInTestPts = (myOffsetStart + s) * nDim;

    std::array<MFloat, nDim> normal;
    std::array<MFloat, nDim> center;
    calculateAveragedNormalAndCenter(segId, normal.data(), center.data());

    MFloat max_diag = (m_solver->c_cellLengthAtLevel(m_solver->grid().maxUniformRefinementLevel())) * SQRT3;

    for(MInt d = 0; d < nDim; d++) {
      test_pts[posInTestPts + d] = center[d] + normal[d] * max_diag;
      m_bndNormals[pos][d] = normal[d];
    }

    m_log << "      * local values for segment " << segId << endl;
    m_log << "        = local normal: (" << normal[0] << " " << normal[1];
    if constexpr(nDim == 3) m_log << " " << normal[2];
    m_log << ")" << endl;
    m_log << "        = local center: (" << center[0] << " " << center[1];
    if constexpr(nDim == 3) m_log << " " << center[2];
    m_log << ")" << endl;
    m_log << "        = local test point: (" << test_pts[posInTestPts] << " " << test_pts[posInTestPts + 1] << " "
          << test_pts[posInTestPts + 2] << endl;
  }

  MPI_Allreduce(MPI_IN_PLACE, test_pts.getPointer(), sumOwners * nDim, MPI_DOUBLE, MPI_SUM, m_solver->mpiComm(), AT_,
                "MPI_IN_PLACE", "test_pts.getPointer()");

  // this runs bascially over all test points and checks for each domain if there exist intersections with the geometry
  for(MInt i = 0; i < sumOwners; i++) {
    MInt pos = i * nDim;
    vector<vector<MInt>> int_nodes;
    m_solver->m_geometry->determineRayIntersectedElements(&test_pts[pos], &int_nodes);

    // now for each direction we need to know if the original ids are unique
    for(MInt d = 0; d < nDim; d++) {
      // this will contain the number of intersections per domain
      MIntScratchSpace doms_with_ints(m_solver->noDomains(), AT_, "doms_with_ints");
      for(MInt dom = 0; dom < m_solver->noDomains(); dom++)
        doms_with_ints[dom] = 0;

      doms_with_ints[m_solver->domainId()] = int_nodes[d].size();

      // after the following doms_with_ints contains for each process the number of intersection for direcetion d
      MPI_Allreduce(MPI_IN_PLACE, doms_with_ints.getPointer(), m_solver->noDomains(), MPI_INT, MPI_SUM,
                    m_solver->mpiComm(), AT_, "MPI_IN_PLACE", "doms_with_ints.getPointer()");

      MInt no_allCuts = 0;
      MInt myOff = 0;
      for(MInt dom = 0; dom < m_solver->noDomains(); dom++)
        no_allCuts += doms_with_ints[dom];

      for(MInt dom = 0; dom < m_solver->domainId(); dom++)
        myOff += doms_with_ints[dom];

      MIntScratchSpace origCmps(no_allCuts, AT_, "origCmps");
      for(MInt l = 0; l < no_allCuts; l++)
        origCmps[l] = 0;

      for(MInt l = 0; l < (MInt)int_nodes[d].size(); l++)
        origCmps[myOff + l] = m_solver->m_geometry->elements[int_nodes[d][l]].m_originalId;

      MPI_Allreduce(MPI_IN_PLACE, origCmps.getPointer(), no_allCuts, MPI_INT, MPI_SUM, m_solver->mpiComm(), AT_,
                    "MPI_IN_PLACE", "origCmps.getPointer()");

      set<MInt> uniques;
      for(MInt l = 0; l < no_allCuts; l++)
        uniques.insert(origCmps[l]);

      io[i * nDim + d] = uniques.size();
    }
  }

  // summarize the number of cuts to dermine the global inside/outside
  MPI_Allreduce(MPI_IN_PLACE, io.getPointer(), sumOwners * nDim, MPI_INT, MPI_SUM, m_solver->mpiComm(), AT_,
                "MPI_IN_PLACE", "io.getPointer()");

  for(MInt s = 0; s < (MInt)own_segments.size(); s++) {
    MInt p = own_segments[s].first;
    MInt posInTestPts = (myOffsetStart + s) * nDim;
    MInt sum = 0;

    for(MInt d = 0; d < nDim; d++)
      sum += io[posInTestPts + d];

    MBool is_inside = sum % 2;

    updateBndNormals(p, is_inside, m_bndNormals[p]);
    m_log << "      * segment " << m_inOutSegmentsIds[p] << endl;
    m_log << "        = normal: (" << m_bndNormals[p][0] << " " << m_bndNormals[p][1] << " " << m_bndNormals[p][2]
          << ")" << endl;
  }
}

/** Caluates the normal for all participatiung process for an in/outlet BC and avarages
 *
 * \author Andreas Lintermann
 * \date 12.06.2016
 *
 * Each process in the BC communicator generates a normal from the triangles and then avarages it together with
 * the other participating processes. All of them furthermore run the inside-outside determination.
 *
 * param[in] own_segments is a vector of tuples of positions and segmentId
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::normalParallelGeomNormals3D(vector<pair<MInt, MInt>> own_segments) {
  TRACE();

  MIntScratchSpace allDoOwn(m_noInOutSegments, AT_, "allDoOwn");
  MIntScratchSpace IDoOwn(m_noInOutSegments, AT_, "IDoOwn");
  MFloatScratchSpace test_pts(m_noInOutSegments, nDim, AT_, "test_pts");

  for(MInt i = 0; i < m_noInOutSegments; i++)
    for(MInt d = 0; d < nDim; d++)
      test_pts(i, d) = 0.0;

  for(MInt i = 0; i < m_noInOutSegments; i++)
    IDoOwn[i] = 0;

  for(MInt s = 0; s < (MInt)own_segments.size(); s++)
    IDoOwn[own_segments[s].first] = 1;


  MPI_Allreduce(IDoOwn.getPointer(), allDoOwn.getPointer(), m_noInOutSegments, MPI_INT, MPI_SUM, m_solver->mpiComm(),
                AT_, "IDoOwn.getPointer()", "allDoOwn.getPointer()");

  // calculate my local normal and the test point (used later) and print to log
  for(MInt s = 0; s < (MInt)own_segments.size(); s++) {
    MInt pos = own_segments[s].first;
    MInt segId = own_segments[s].second;

    std::array<MFloat, nDim> normal;
    std::array<MFloat, nDim> center;
    calculateAveragedNormalAndCenter(segId, normal.data(), center.data());

    MFloat max_diag = (m_solver->c_cellLengthAtLevel(m_solver->grid().maxUniformRefinementLevel())) * SQRT3;

    for(MInt d = 0; d < nDim; d++) {
      test_pts(pos, d) = center[d] + normal[d] * max_diag;
      m_bndNormals[pos][d] = normal[d];
    }

    m_log << "      * local values for segment " << segId << endl;
    m_log << "        = local normal: (" << normal[0] << " " << normal[1];
    if constexpr(nDim == 3) m_log << " " << normal[2];
    m_log << ")" << endl;
    m_log << "        = local center: (" << center[0] << " " << center[1];
    if constexpr(nDim == 3) m_log << " " << center[2];
    m_log << ")" << endl;
    m_log << "        = local test point: (" << test_pts(pos, 0) << " " << test_pts(pos, 1) << " " << test_pts(pos, 2)
          << endl;
  }

  // this is now for all the others

  // count the number of non-single entries in allDoOwn
  // also store the roots of segment ids for domains that own the segment completely
  MIntScratchSpace owner_roots(m_solver->noDomains(), AT_, "owner_roots");
  MInt numComm = 0;
  vector<MInt> posToComm;
  for(MInt i = 0; i < m_noInOutSegments; i++) {
    if(allDoOwn[i] == 1 && IDoOwn[i] == 1)
      owner_roots[m_solver->domainId()] = i;
    else
      owner_roots[m_solver->domainId()] = -1;

    if(allDoOwn[i] > 1) {
      posToComm.push_back(i);
      numComm++;
    }
  }

  // create the communicators
  vector<MPI_Comm> comms;
  for(MInt i = 0; i < numComm; i++) {
    MInt pos = posToComm[i];

    // place holders for all owners
    MInt sumowners = 0;
    MInt firstOwner = -1;
    MIntScratchSpace owners(m_solver->noDomains(), AT_, "owners");
    MPI_Comm tmpComm;

    // first determine all owners
    m_solver->m_geometry->determineSegmentOwnership(m_inOutSegmentsIds[pos], &IDoOwn[pos], &sumowners, &firstOwner,
                                                    owners.getPointer());

    if(m_solver->domainId() == firstOwner) owner_roots[m_solver->domainId()] = pos;
    m_log << "      * creating MPI communicators for segment " << m_inOutSegmentsIds[pos] << endl;
    m_log << "        = sum of owners:             " << sumowners << endl;
    m_log << "        = root of communication:     " << firstOwner << endl;
    m_log << "        = owners:                    ";
    for(MInt d = 0; d < m_solver->noDomains(); d++)
      if(owners[d] > 0) m_log << d << " ";
    m_log << endl;

    // build communicator for the subgroup
    m_solver->createMPIComm(owners.getPointer(), sumowners, &tmpComm);
    comms.push_back(tmpComm);
  }
  MPI_Allreduce(MPI_IN_PLACE, owner_roots.getPointer(), m_solver->noDomains(), MPI_INT, MPI_SUM, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "owner_roots.getPointer()");

  // do the exchange
  for(MInt i = 0; i < numComm; i++) {
    MInt p = posToComm[i];
    // I am an owner, I need to participate in the communication
    if(IDoOwn[p]) {
      // sum the normals and test points and average
      MFloatScratchSpace testenv(2 * nDim, AT_, "testenv");
      MInt j = 0;
      for(MInt d = 0; d < nDim; d++, j++)
        testenv[j] = m_bndNormals[p][d];
      for(MInt d = 0; d < nDim; d++, j++)
        testenv[j] = test_pts(p, d);

      MPI_Allreduce(MPI_IN_PLACE, testenv.getPointer(), 2 * nDim, MPI_DOUBLE, MPI_SUM, comms[i], AT_, "MPI_IN_PLACE",
                    "testenv.getPointer()");

      MInt size;
      MPI_Comm_size(comms[i], &size);

      // average over the number of domains
      for(MInt d = 0; d < 2 * nDim; d++)
        testenv[d] /= size;

      // now testenv contains the averaged normal and the averaged test point
      // copy back to original position
      j = 0;
      for(MInt d = 0; d < nDim; d++, j++)
        m_bndNormals[i][d] = testenv[j];
      for(MInt d = 0; d < nDim; d++, j++)
        test_pts(p, d) = testenv[j];
    }
  }

  // at this position all domains owning (a piece) of a segment have their normal and test point available

  // now test all the test points
  // allocate space on all domains for the testing points
  MFloatScratchSpace rec_test_pts(nDim * m_noInOutSegments, AT_, "rec_test_pts");
  for(MInt i = 0; i < nDim * m_noInOutSegments; i++)
    rec_test_pts[i] = 0.0;

  for(MInt d = 0; d < m_solver->noDomains(); d++) {
    MInt segInOut = owner_roots[m_solver->domainId()];
    if(segInOut >= 0) {
      MInt p = segInOut * nDim;
      for(MInt k = 0; k < nDim; k++)
        rec_test_pts[p + k] = test_pts(segInOut, k);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, rec_test_pts.getPointer(), nDim * m_noInOutSegments, MPI_DOUBLE, MPI_SUM,
                m_solver->mpiComm(), AT_, "MPI_IN_PLACE", "rec_test_pts.getPointer()");

  // now perform inside/outside determination
  MIntScratchSpace num_cuts(m_noInOutSegments * nDim, AT_, "num_cuts");

  for(MInt i = 0; i < m_noInOutSegments; i++) {
    MInt pos = i * nDim;
    m_solver->m_geometry->pointIsInside2(&rec_test_pts[pos], &num_cuts[pos]);
  }

  // summarize the number of cuts to dermine the global inside/outside
  MPI_Allreduce(MPI_IN_PLACE, num_cuts.getPointer(), nDim * m_noInOutSegments, MPI_INT, MPI_SUM, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "num_cuts.getPointer()");

  for(MInt s = 0; s < (MInt)own_segments.size(); s++) {
    MInt p = own_segments[s].first;
    MInt pos = p * nDim;
    MInt sum = 0;

    for(MInt d = 0; d < nDim; d++)
      sum += num_cuts[pos + d];

    MBool is_inside = sum % 2;

    updateBndNormals(p, is_inside, m_bndNormals[p]);
    m_log << "      * segment " << m_inOutSegmentsIds[p] << endl;
    m_log << "        = normal: (" << m_bndNormals[p][0] << " " << m_bndNormals[p][1] << " " << m_bndNormals[p][2]
          << ")" << endl;
  }
}


template <MInt nDim>
void LbBndCnd<nDim>::calculateBndNormals() {
  TRACE();
  // TODO labels:LB dxqy: enable also for 2D
  std::stringstream ss;
  ss << "Calculate boundary normal not implemented for " << nDim << "D, yet !";
  TERMM(1, ss.str());
}

/** \brief Calculates the averaged normal on a boundary segment.
 *
 * \author Andreas Lintermann
 * \date 23.01.2013, 27.04.2015, 14.01.2016
 *
 * For each segment the averaged normal \f$n\f$ is calculated by running over
 * all triangles and adding either the cross-product of the spanning edges or by reading the
 * normals stored in the STL.
 * Additionally the center \f$c\f$ is evaluated by averaging the vertex locations
 * over the number of vertices. Then it is tested if
 *
 * \f[\vec{p}=\vec{c}+d\cdot\vec{n}\f],
 *
 * where $d$ is the length of the diagonal of the smallest cell, is inside the
 * geometry (using pointIsInside2 from CartesianGrid).
 *
 **/
template <>
void LbBndCnd<3>::calculateBndNormals() {
  constexpr MInt nDim = 3;
  TRACE();

  m_log << "    - calculating normals in 3D (bndNormalMethod is set to " << m_bndNormalMethod << ")" << endl;

  // 1. find out what segments we own and print to log
  vector<pair<MInt, MInt>> own_segments;
  for(MInt i = 0; i < m_noInOutSegments; i++) {
    for(MInt o = 0; o < (MInt)m_bndCndSegIds.size(); o++)
      if(m_inOutSegmentsIds[i] == m_bndCndSegIds[o]) {
        own_segments.push_back(pair<MInt, MInt>(i, m_inOutSegmentsIds[i]));
        break;
      }

    // initialize all by (0,0,0)
    for(MInt d = 0; d < nDim; d++)
      m_bndNormals[i][d] = 0.0;
  }

  if(m_solver->m_geometry->m_parallelGeometry) {
    if(m_fastParallelGeomNormals)
      fastParallelGeomNormals3D(own_segments);
    else
      normalParallelGeomNormals3D(own_segments);
  } else {
    for(MInt s = 0; s < (MInt)own_segments.size(); s++) {
      MInt pos = own_segments[s].first;
      MInt segId = own_segments[s].second;
      std::array<MFloat, nDim> normal;
      std::array<MFloat, nDim> center;
      calculateAveragedNormalAndCenter(segId, normal.data(), center.data());

      MFloat max_diag = (m_solver->c_cellLengthAtLevel(m_solver->grid().maxUniformRefinementLevel())) * SQRT3;
      std::array<MFloat, nDim> test_pt;

      for(MInt d = 0; d < nDim; d++)
        test_pt[d] = center[d] + normal[d] * max_diag;

      updateBndNormals(pos, !m_solver->m_geometry->pointIsInside2(test_pt.data()), normal.data());
      m_log << "      * segment " << segId << endl;
      m_log << "        = normal: (" << m_bndNormals[pos][0] << " " << m_bndNormals[pos][1] << " "
            << m_bndNormals[pos][2] << ")" << endl;
      m_log << "        = center: (" << center[0] << " " << center[1] << " " << center[2] << ")" << endl;
    }
  }
}

/** \brief Updates the normals of an inlet/outlet based on inside / outside detection
 *
 * \author Andreas Lintermann
 * \date 21.09.2015
 *
 * \param[in] segId the id of the segment
 * \param[in] is_inside the result of the inside / outside determination
 * \param[in] svg_normal the averaged normal used to override
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::updateBndNormals(MInt segId, MBool is_inside, MFloat* avg_normal) {
  // init the right direction, they should always point outside
  m_bndNormalDirection[segId] = -1;

  // outside
  if(!is_inside) {
    for(MInt d = 0; d < nDim; d++)
      m_bndNormals[segId][d] = avg_normal[d];
  }
  // inside
  else {
    for(MInt d = 0; d < nDim; d++)
      m_bndNormals[segId][d] = -1.0 * avg_normal[d];
  }
}

/** \brief This function returns the index in the array of the boundary conditions that are inflow/outflow conditions
   for a given index in all boundary cells. \author Andreas Lintermann \date 27.01.2010

    The function simply runs over the the array containing the segment ids of the boundaries that are inflow/outflow and
   checks if in the original array of all boundary cell at a given offset this segment id is found.

  */
template <MInt nDim>
MInt LbBndCnd<nDim>::findBndCnd(MInt index) {
  TRACE();

  MInt ind = -1;
  for(MInt i = 0; i < m_noInOutSegments; i++)
    if(m_bndCells[m_bndCndOffsets[index]].m_segmentId[0] == m_inOutSegmentsIds[i]) ind = i;

  return ind;
}

/** \brief Allocate data for given boundary index
 *  \author Miro Gondrum
 *  \date 09.02.2021
 *  \param[in]  index the index of the BC
 */
template <MInt nDim>
void LbBndCnd<nDim>::bcDataAllocate(MInt index, MInt noVars) {
  TRACE();
  m_bndCndData[index]; // creates default entry of bndCndData
  LbBndCndData& p = m_bndCndData[index];
  p.noVars = noVars;
  p.noBndCells = 0;
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    if(!m_solver->a_isHalo(m_bndCells[i].m_cellId)) {
      p.noBndCells++;
    }
  }
  p.noBndCellsWithHalos = m_bndCndOffsets[index + 1] - m_bndCndOffsets[index];
  const MInt noEntries = p.noBndCellsWithHalos * p.noVars;
  mAlloc(p.data, noEntries, "m_bndCndData[" + std::to_string(index) + "].data", F0, AT_);
  m_segIdUseBcData[m_bndCndSegIds[index]] = 1;
}

/** \brief Write bndCndData info in given ParallelIo file's header
 *  \author Miro Gondrum
 *  \date 10.02.2021
 *  \param[in]  parallelIo  Parallel file handler opened in define mode
 *  \note Needs to be called before bcDataWriteRestartData
 */
template <MInt nDim>
void LbBndCnd<nDim>::bcDataWriteRestartHeader(ParallelIo& parallelIo) {
  TRACE();
  using namespace maia::parallel_io;
  // Here it is looped over all global segment ids. Since most values are only
  // set for local segment ids a mapping is performed. In case of local missing
  // segment the number of cells per segment is zeroed.
  for(MInt i = 0; i < m_noSegments; i++) {
    if(m_segIdUseBcData[i]) {
      const MInt localSegIndex = m_mapBndCndSegId2Index[i];
      if(localSegIndex != -1) {
        // local domain hold this segment
        LbBndCndData& p = m_bndCndData[localSegIndex];
        const ParallelIo::size_type noEntries = p.noBndCells * p.noVars;
        ParallelIo::size_type noEntriesGlobal;
        ParallelIo::calcOffset(noEntries, &p.globalOffset, &noEntriesGlobal, m_solver->mpiComm());
        if(noEntriesGlobal > 0) {
          const MString name = "bcData_segment" + to_string(i);
          parallelIo.defineArray(PIO_FLOAT, name, noEntriesGlobal);
          parallelIo.setAttribute(name, "name", name);
        }
      } else {
        // local domain does NOT hold this segment
        const ParallelIo::size_type noEntries = 0;
        ParallelIo::size_type noEntriesGlobal, offset;
        ParallelIo::calcOffset(noEntries, &offset, &noEntriesGlobal, m_solver->mpiComm());
        if(noEntriesGlobal > 0) {
          const MString name = "bcData_segment" + to_string(i);
          parallelIo.defineArray(PIO_FLOAT, name, noEntriesGlobal);
          parallelIo.setAttribute(name, "name", name);
        }
      }
    }
  }
}

/** \brief  Write bndCndData data in given ParallelIo file
 *  \author Miro Gondrum
 *  \date 11.02.2021
 *  \param[in]  parallelIo  Parallel file handler opened in write mode
 *  \note Needs to be called after bcDataWriteRestartHeader
 */
template <MInt nDim>
void LbBndCnd<nDim>::bcDataWriteRestartData(ParallelIo& parallelIo) {
  TRACE();
  using namespace maia::parallel_io;
  for(MInt i = 0; i < m_noSegments; i++) {
    if(m_segIdUseBcData[i]) {
      const MInt localSegIndex = m_mapBndCndSegId2Index[i];
      const MString name = "bcData_segment" + to_string(i);
      if(localSegIndex != -1) {
        // local domain hold this segment
        LbBndCndData& p = m_bndCndData[localSegIndex];
        const MInt noEntries = p.noBndCells * p.noVars;
        // create tmp array only containing non-halo data that are written to target file
        MFloatScratchSpace tmp(max(noEntries, 1), AT_, "tmp");
        MInt l = 0;
        const MInt offset = m_bndCndOffsets[localSegIndex];
        for(MInt j = 0; j < p.noBndCellsWithHalos; j++) {
          if(!m_solver->a_isHalo(m_bndCells[j + offset].m_cellId)) {
            for(MInt k = 0; k < p.noVars; k++) {
              tmp[l * p.noVars + k] = p.data[j * p.noVars + k];
            }
            l++;
          }
        }
        parallelIo.setOffset(noEntries, p.globalOffset);
        parallelIo.writeArray(&tmp[0], name);
      } else {
        // local domain does NOT hold this segment -> write no data
        const MFloat* dummyData = nullptr;
        parallelIo.setOffset(0, 0);
        parallelIo.writeArray(&dummyData[0], name);
      }
    }
  }
}

/** \brief  Read bndCndData data in given ParallelIo file
 *  \author Miro Gondrum
 *  \date 11.02.2021
 *  \param[in]  parallelIo  Parallel file handler opened in write mode
 */
template <MInt nDim>
void LbBndCnd<nDim>::bcDataReadRestartData(ParallelIo& parallelIo) {
  TRACE();
  using namespace maia::parallel_io;
  for(MInt i = 0; i < m_noSegments; i++) {
    if(m_segIdUseBcData[i]) {
      const MInt localSegIndex = m_mapBndCndSegId2Index[i];
      if(localSegIndex != -1) {
        // local domain hold this segment
        LbBndCndData& p = m_bndCndData[localSegIndex];
        const MInt noEntries = p.noBndCells * p.noVars;
        ParallelIo::size_type noEntriesGlobal;
        ParallelIo::calcOffset(noEntries, &p.globalOffset, &noEntriesGlobal, m_solver->mpiComm());
        const MString name = "bcData_segment" + to_string(i);
        // create tmp array only containing non-halo data that are written to target file
        MFloatScratchSpace tmp(max(noEntries, 1), AT_, "tmp");
        parallelIo.setOffset(noEntries, p.globalOffset);
        parallelIo.readArray(&tmp[0], name);
        // fill tmp array into data (which also contains halo data)
        MInt l = 0;
        const MInt offset = m_bndCndOffsets[localSegIndex];
        for(MInt j = 0; j < p.noBndCellsWithHalos; j++) {
          if(!m_solver->a_isHalo(m_bndCells[j + offset].m_cellId)) {
            for(MInt k = 0; k < p.noVars; k++) {
              p.data[j * p.noVars + k] = tmp[l * p.noVars + k];
            }
            l++;
          }
        }
      } else {
        // local domain does NOT hold this segment -> only participate in offset calculation
        const ParallelIo::size_type noEntries = 0;
        ParallelIo::size_type noEntriesGlobal, offset;
        ParallelIo::calcOffset(noEntries, &offset, &noEntriesGlobal, m_solver->mpiComm());
        const MString name = "bcData_segment" + to_string(i);
        MFloat* dummyData = nullptr;
        parallelIo.setOffset(noEntries, offset);
        parallelIo.readArray(&dummyData[0], name);
      }
    }
  }
}

/** \brief This function checks if for an inflow boundary the normal points into the according direction and changes it
 *depending on it. \author Andreas Lintermann \date 27.01.2010, 27.04.2015
 *
 * A previously calculated direction array is used to check if the normal points into the right direction. If not, the
 *normal is multiplicated with -1 and the direction indicator is switched. Note that bounday normals alway point
 *outside.
 *
 * \param[in] index the index of the BC
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::applyDirectionChangeInflow(MInt index) {
  TRACE();

  if(m_initVelocityMethod != "read" && m_initVelocityMethod != "fromSTL") {
    MInt ind = m_mapSegIdsInOutCnd[index];

    // wrong direction, mutiplicate with -1
    if(m_bndNormalDirection[ind] == -1)
      for(MInt i = 0; i < nDim; i++)
        m_initialVelocityVecs[ind][i] *= -1.0;
  }
}

/** \brief This function checks if for an outflow boundary the normal points into the according direction and changes it
 *depending on it. \author Andreas Lintermann \date 27.01.2010
 *
 * A previously calculated direction array is used to check if the normal points into the right direction. If not, the
 *normal is multiplicated with -1 and the direction indicator is switched. Note that boundary normals always point
 *outside.
 *
 * \param[in] index the index of the BC
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::applyDirectionChangeOutflow(MInt index) {
  TRACE();

  if(m_initVelocityMethod != "read" && m_initVelocityMethod != "fromSTL") {
    MInt ind = m_mapSegIdsInOutCnd[index];

    // wrong direction, mutiplicate with -1
    if(m_bndNormalDirection[ind] == 1)
      for(MInt i = 0; i < nDim; i++)
        m_initialVelocityVecs[ind][i] *= -1.0;
  }
}

/** \brief Solves the Blasius equation for f,f',f".
 *
 * \author Andreas Lintermann
 * \date 24.05.2011
 *
 * Solves the Blasius equation for cells in positive z-direction and stores the result in an array m_blasius.
 * The values are calculated beginning with a starting f" = 0.332051914927913096446.The similarity parameter eta of the
 * boundary cells is stored in the boundary cell variable m_eta.
 *
 * \param index The index ot the boundary to be treated
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::solveBlasiusZ(MInt index) {
  TRACE();
  MInt currentId;
  MFloat minmax[2] = {10000000.0, -10000000.0};

  MInt x_pos = m_solver->m_referenceLength * m_solver->m_blasiusPos;
  MFloat z_pos = 0.0;

  set<MFloat> etas;
  MFloat eta;

  MFloat etamax;
  MFloat fppwall = 0.332051914927913096446;

  // This is fine enough!
  MFloat deta = 0.001;
  MInt steps;

  // Find min and max Z
  for(MInt i = m_bndCndOffsets[index], j = 0; i < m_bndCndOffsets[index + 1]; i++, j++) {
    if(m_solver->c_noChildren(m_bndCells[i].m_cellId) > 0) continue;
    if(m_solver->a_coordinate(m_bndCells[i].m_cellId, nDim - 1) < minmax[0])
      minmax[0] = m_solver->a_coordinate(m_bndCells[i].m_cellId, nDim - 1);
    if(m_solver->a_coordinate(m_bndCells[i].m_cellId, nDim - 1) > minmax[1])
      minmax[1] = m_solver->a_coordinate(m_bndCells[i].m_cellId, nDim - 1);
  }

  // Find all etas and store in m_eta
  for(MInt i = m_bndCndOffsets[index]; i < m_bndCndOffsets[index + 1]; i++) {
    currentId = m_bndCells[i].m_cellId;
    if(m_solver->c_noChildren(currentId) > 0) {
      m_bndCells[i].m_eta = 0.0;
      continue;
    }
    z_pos = ((m_solver->a_coordinate(currentId, nDim - 1)) - minmax[0]
             + ((m_solver->c_cellLengthAtLevel(m_solver->maxLevel())) / 2.0))
            / (m_solver->c_cellLengthAtLevel(m_solver->maxLevel()));
    eta = z_pos * sqrt((m_Ma * LBCS) / (x_pos * m_solver->m_nu));
    etas.insert(eta);

    m_bndCells[i].m_eta = eta;
  }

  // Sort entries
  list<MFloat> etas_sorted;
  for(set<MFloat>::iterator i = etas.begin(); i != etas.end(); i++)
    etas_sorted.push_back(*i);

  etas_sorted.sort();

  list<MFloat>::iterator it1 = etas_sorted.end();
  it1--;
  list<MFloat>::iterator it2 = etas_sorted.begin();
  it2++;


  // Calculate Blasius solution
  etamax = *(it1);
  it1 = etas_sorted.begin();

  steps = etamax / deta;

  // This is accurate enough
  m_blasius_delta = deta;

  mAlloc(m_blasius, steps + 1, 5, "m_blasius", F0, AT_);
  m_blasius[0][3] = fppwall;

  // Lets use the shooting method for getting the Blasius solution
  for(MInt i = 0; i < steps; i++) {
    m_blasius[i + 1][0] = m_blasius[i][0] + deta;
    m_blasius[i + 1][1] = m_blasius[i][1] + m_blasius[i][2] * deta;
    m_blasius[i + 1][2] = m_blasius[i][2] + m_blasius[i][3] * deta;
    m_blasius[i + 1][3] = m_blasius[i][3] + (-0.5 * m_blasius[i][1] * m_blasius[i][3]) * deta;
  }
}

/** \brief This function sorts the boundary cells according to the BC id.
    \author Andreas Lintermann
    \date 27.01.2010

    The sorting of the boundary cells is necessary because all cells are
    stored in one single collector. The boundary condition functions get
    the position of a type of boundary cells in the collector (offset) and
    the number of boundary cells of equal type ( i.e. with the same
    boundary id) as parameter.
  */
template <MInt nDim>
void LbBndCnd<nDim>::sortBoundaryCells() {
  TRACE();

  // TODO labels:LB,TIMERS fix timers when sortBoundaryCells is called from lbbndcnddxqy.cpp
  /* NEW_SUB_TIMER(t_sortBCCells, "sort boundary cells", m_t_BCAll); */
  /* RECORD_TIMER_START(t_sortBCCells); */

  // Sort the boundary cells
  const MInt tmpCellSize = m_bndCells.size();

  m_log << "  + Sorting boundary cells..." << endl;
  m_log << "    - no. boundary cells: " << tmpCellSize << endl;

  // using stable_sort to preserve the order by cellId within segmentIds
  stable_sort(m_bndCells.begin(), m_bndCells.end(),
              [](auto& a, auto& b) { return a.m_segmentId[0] < b.m_segmentId[0]; });

  MInt tmpSegmentId = -1; // holds the current id
  MInt tmpBndCndId = -1;  // holds the current boundary cells bndCndId
  MInt counter = 0;       // Counts the boundary cells

  // setting some default states, which is important in case of redoing this
  // method call
  m_bndCndIds.clear();
  m_bndCndOffsets.clear();
  m_bndCndSegIds.clear();
  m_mapBndCndSegId2Index.resize(m_noSegments, -1);
  m_mapIndex2BndCndSegId.resize(m_noSegments, -1);
  m_noBndCellsPerSegment.resize(m_noSegments, 0);
  // TODO labels:LB maybe also useful to reset: a_isBndryCell, a_bndId

  for(MInt i = 0; i < tmpCellSize; i++) {
    if(tmpSegmentId != m_bndCells[i].m_segmentId[0]) {
      // Since m_bndCells has been sorted previously and new segmentid differs
      // from tmpSegmentId, we sort now for a new tmpSegmentId and co
      tmpSegmentId = m_bndCells[i].m_segmentId[0];
      tmpBndCndId = m_bndCells[i].m_bndCndId[0];

      m_bndCndIds.push_back(tmpBndCndId);
      m_bndCndOffsets.push_back(counter);
      m_bndCndSegIds.push_back(tmpSegmentId);
      m_mapBndCndSegId2Index[tmpSegmentId] = m_bndCndSegIds.size() - 1;
      m_mapIndex2BndCndSegId[m_bndCndSegIds.size() - 1] = tmpSegmentId;
    }
    m_noBndCellsPerSegment[tmpSegmentId]++;
    m_solver->a_isBndryCell(m_bndCells[counter].m_cellId) = true;
    m_solver->a_bndId(m_bndCells[counter].m_cellId) = counter;

    counter++;
  }

  m_bndCndOffsets.push_back(counter);

  // Make the mapping from all segments ids to the one that are just in/outflow
  m_mapSegIdsInOutCnd.resize((MInt)m_bndCndSegIds.size());
  for(MInt i = 0; i < (MInt)(m_bndCndSegIds.size()); i++) {
    MInt seg_id = m_bndCndSegIds[i];

    m_log << "    - BC " << m_bndCndIds[i] << endl;
    m_log << "      * index:      " << i << endl;
    m_log << "      * segment id: " << seg_id << endl;
    m_log << "      * no cells:   " << m_noBndCellsPerSegment[seg_id] << endl;
    m_log << "      * offsets:    " << m_bndCndOffsets[i] << " - " << m_bndCndOffsets[i + 1] - 1 << endl;

    MInt pos = -1;
    for(MInt j = 0; j < m_noInOutSegments; j++)
      if(seg_id == m_inOutSegmentsIds[j]) pos = j;

    m_mapSegIdsInOutCnd[i] = pos;

    m_log << "      * mapping:    " << m_mapSegIdsInOutCnd[i] << " "
          << ((m_mapSegIdsInOutCnd[i] < 0) ? "(this is not an inflow / outflow BC)" : "") << endl;
  }
  m_log << endl;

  /* RECORD_TIMER_STOP(t_sortBCCells); */
}

//! Sets the BndCndHandler objects at solver setup
/*! setBndCndHandler() needs to be called after loadBndCells(). It needs the
    boundary condition to set the boundary condition handler.
    Each new boundary condition needs to be implemented here, too.
*/
template <MInt nDim>
void LbBndCnd<nDim>::setBndCndHandler() {
  TRACE();
  if(!m_solver->isActive()) return;

  NEW_SUB_TIMER(t_setBCHandler, "set BC handler", m_t_BCAll);
  RECORD_TIMER_START(t_setBCHandler);

  m_log << "  + Setting the boundary condition handler..." << endl << endl;

  // Allocate space for pointer lists which store the boundary functions
  bndCndHandlerVariables.resize(m_bndCndIds.size());
  bndCndHandlerRHS.resize(m_bndCndIds.size());
  bndCndInitiator.resize(m_bndCndIds.size());

  m_segIdUseBcData.resize(m_noSegments, 0);

  // Fill the function pointer lists with correct bc functions
  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
    switch(m_bndCndIds[i]) {
      case 0:
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        DEBUG("ERROR: bndCndHandlers only point to a dummy function", MAIA_DEBUG_LEVEL1);
        break;

        //----------------------------------------------------
        // velocity boundary conditions

      case 1000: // velocity (eq)
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10000;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1001: // velocity (eq)
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10001;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1002: // Inflow condition for periodic channel (z periodic, x/y is inflow)
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10002;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1004: // velocity (eq)
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10004;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1022: // Inflow condition for periodic channel (z periodic, x/y is inflow)
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10022;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1010: // Inflow condition
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10010;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1020: // Inflow condition
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc10020;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1040: // velocity NRCBC Izquerdo et al. 2008
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10040;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1041:
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10041;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1042:
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10042;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1043:
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10043;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1044:
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10044;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1045:
        bcDataAllocate(i, nDim + 1);
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10045;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1046: // velocity CBC LODI (force equilibrium)
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10046;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1050:
        applyDirectionChangeInflow(i);
        bndCndHandlerRHS[i] = &LbBndCnd::bc10050;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1060: // velocity (eq), thermal
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10060;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1061: // velocity (eq), thermal, Blasius for velocity
        applyDirectionChangeInflow(i);
        solveBlasiusZ(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10061;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 1070: // keeping local inflow mass flux to rho_0*Ma*Cs
        if(m_densityFluctuations) TERMM(1, "BC 1070 not working for densityFluctuations");
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10070;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 1090: // velocity (eq), transport, concentration = 0.0
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10090;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

        //----------------------------------------------------
        // Pulsatile BCs

      case 1080: // velocity (eq) - sinus
        applyDirectionChangeInflow(i);
        bndCndHandlerVariables[i] = &LbBndCnd::bc10080;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;


        //----------------------------------------------------
        // Povitsky

      case 1111: // Povitsky cavity condition (equilibirum dist)
        bndCndHandlerVariables[i] = &LbBndCnd::bc10111;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

        //----------------------------------------------------
        // wall boundary conditions

      case 2000: // no slip, BFL
        bndCndHandlerRHS[i] = &LbBndCnd::bc20000;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2001: // no slip, simple bounce back, halfway between nodes
        bndCndHandlerRHS[i] = &LbBndCnd::bc20001;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2002: // no slip, simple bounce back, on nodes
        bndCndHandlerRHS[i] = &LbBndCnd::bc20002;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2003: // no slip, equilibrium functions
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc20003;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2004: // no slip, Haenel interpolated
        bndCndHandlerRHS[i] = &LbBndCnd::bc20004;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2005: // no slip, Guo interpolated
        bndCndHandlerRHS[i] = &LbBndCnd::bc20005;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2006: // free slip, only for 2d
        bndCndHandlerRHS[i] = &LbBndCnd::bc20006;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2010: // no slip, BFL, with force calculation
        bndCndHandlerRHS[i] = &LbBndCnd::bc20010;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2020: // no slip, BFL, with force calculation
        bndCndHandlerRHS[i] = &LbBndCnd::bc20020;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2220: // no slip, BFL , Thermal
        bndCndHandlerRHS[i] = &LbBndCnd::bc20220;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 2226: // no slip, BFL, Thermal (Dirichlet condition), cylinder flow
        bndCndHandlerRHS[i] = &LbBndCnd::bc20226;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2227: // no slip, BFL, Thermal (Dirichlet condition), channel flow
        bndCndHandlerRHS[i] = &LbBndCnd::bc20227;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2228: // no slip, BFL, Thermal flux (Neumann condition), channel flow
        bndCndHandlerRHS[i] = &LbBndCnd::bc20228;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bcIBBNeumannInit;
        break;

      case 2022: // no slip, BFL , also for Thermal
        bndCndHandlerRHS[i] = &LbBndCnd::bc20022;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2023: // no slip, simple bounce back, on nodes for Thermal
        bndCndHandlerRHS[i] = &LbBndCnd::bc20023;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2024: // no slip, Thermal bc2004
        bndCndHandlerRHS[i] = &LbBndCnd::bc20024;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2025: // no slip, Thermal, sets eq. dist. fnc. for given mac. vars.
        bndCndHandlerRHS[i] = &LbBndCnd::bc20025;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2026: // no slip, simple bounce back, halfway between nodes, same as 2001, but also for Thermal
        bndCndHandlerRHS[i] = &LbBndCnd::bc20026;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2027: // no slip, interpolated bounce back, thermal interpolated
        bndCndHandlerRHS[i] = &LbBndCnd::bc20027;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2030: // no slip, Haenel interpolated with velocity
        bndCndHandlerRHS[i] = &LbBndCnd::bc20030;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2230: // no slip, Haenel interpolated with velocity
        bndCndHandlerRHS[i] = &LbBndCnd::bc20230;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2050: // sliding wall in axis direction
        bndCndHandlerVariables[i] = &LbBndCnd::bc20050;
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2051:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20051;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2052:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20052;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2053:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20053;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2054:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20054;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 2055:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20055;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 2501: // nasal mucosa model, latent heat
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc20501;
        bndCndInitiator[i] = &LbBndCnd::bc20501_init;
        break;

        //----------------------------------------------------
        // zero gradient boundary conditions

        //----------------------------------------------------
      case 3000: // extrapolation in normal vector direction
        bndCndHandlerRHS[i] = &LbBndCnd::bc30000;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 3010: // extrapolation in axis direction
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30010;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3011:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30011;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3012:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30012;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3013:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30013;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3014:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30014;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3015:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30015;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 3020:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30020;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3021:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30021;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3022:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30022;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3023:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30023;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3024:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30024;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3025:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30025;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 3030:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30030;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3031:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30031;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3032:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30032;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3033:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30033;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3034:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30034;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3035:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30035;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 3040:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30040;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3041:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30041;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3042:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30042;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3043:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30043;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3044:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30044;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 3045:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30045;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 3050:
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndHandlerRHS[i] = &LbBndCnd::bc30050;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;


        //----------------------------------------------------
        // pressure bc's

      case 4000: // prescribe equilibrium distributions
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40000;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4001: // prescribe equilibrium distributions
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40001;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4010: // Extrapolation Chen
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40010;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4020: // relaxed pressure
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40020;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4030: // Extrapolation with non-eq
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40030;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4040: // pressure NRCBC Izquerdo et al. 2008
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40040;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4041:
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40041;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4042:
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40042;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4043:
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40043;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4044:
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40044;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4045:
        bcDataAllocate(i, nDim + 1);
        bndCndHandlerRHS[i] = &LbBndCnd::bc40045;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 4046: // non-reflective based on LODI by prescribing equilibrium
        bndCndHandlerRHS[i] = &LbBndCnd::bc40046;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 4060:
        bndCndHandlerRHS[i] = &LbBndCnd::bc40060;
        bndCndHandlerVariables[i] = &LbBndCnd::bc0;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 4070: // relax pressure based on formulations by Hoerschler
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40070;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4071: // hold pressure constant at a given value
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40071;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4072: // adjusts pressure to yield a specified local Reynolds number
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40072;
        bndCndInitiator[i] = &LbBndCnd::bc40072_40082_init;
        break;
      case 4073: // adjusts pressure to yield a specified local Reynolds number (measured at a specified location)
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40073;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4080: // relax pressure based on formulations by Hoerschler, additionally set temperature
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40080;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4081: // hold pressure constant at a given value, extrapolate temperature and velcoities
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40081;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4082: // adjusts pressure to yield a specified local Reynolds number
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40082;
        bndCndInitiator[i] = &LbBndCnd::bc40072_40082_init;
        break;

      case 4100: // prescribe equilibrium distributions
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40100;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      case 4110: // prescribe equilibrium distributions (similar to bc4030)
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40110;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 4120: // prescribe equilibrium distributions (similar to bc4030)
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40120;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;

      case 4130: // Extrapolation with non-eq
        bndCndHandlerRHS[i] = &LbBndCnd::bc0;
        bndCndHandlerVariables[i] = &LbBndCnd::bc40130;
        bndCndInitiator[i] = &LbBndCnd::bc0;
        break;
      default:
        stringstream errorMessage;
        errorMessage << " LbBncCndDxQy::setBndCndHandler : Unknown boundary condition " << m_bndCndIds[i]
                     << " exiting program.";
        TERMM(1, errorMessage.str());
    }
  }

  // Communications
  MPI_Allreduce(MPI_IN_PLACE, m_segIdUseBcData.data(), m_segIdUseBcData.size(), MPI_INT, MPI_SUM, m_solver->mpiComm(),
                AT_, "MPI_IN_PLACE", "m_SegIdUseBcData.data()");

  RECORD_TIMER_STOP(t_setBCHandler);
}


/** \brief Checks if a BC exists that requires communication
 *
 * \author Andreas Lintermann
 * \date 29.05.2015
 *
 * \return a MBool that defines if such a BC exists or not
 *
 **/
template <MInt nDim>
MBool LbBndCnd<nDim>::checkForCommForce() {
  TRACE();

  m_noAllBoundaryIds = 0;
  m_allBoundaryIds = m_solver->m_geometry->GetBoundaryIds(&m_noAllBoundaryIds);

  for(MInt i = 0; i < m_noAllBoundaryIds; i++) {
    switch(m_allBoundaryIds[i]) {
      case 2000: {
        return true;
        // The following statement was commented since it is not reachable. Remove if not needed anymore.
        // break;
      }
      default: {
      }
    }
  }
  return false;
}

/** \brief Sets up a neighbor-communicator for certain BCs
 *
 * \author Andreas Lintermann
 * \date 20.12.2019
 *
 * This function checks if a BC is used which requires internal communication.
 * If so, a new communicator group is created holding only the processes having a
 * part of the BC.
 **/
template <MInt nDim>
void LbBndCnd<nDim>::setBCWallNeighborCommunicator() {
  TRACE();

  NEW_SUB_TIMER(t_setBCNC, "set BC neighbor communicators", m_t_BCAll);
  RECORD_TIMER_START(t_setBCNC);

  // first of all check if we have a BC which requries communication
  m_hasCommForce = checkForCommForce();
  if(!m_hasCommForce) {
    RECORD_TIMER_STOP(t_setBCNC);
    return;
  }

  m_log << "  + Found a BC which requires communication:" << endl;

  // currently working only for one BC
  for(MInt i = 0; i < m_noAllBoundaryIds; i++) {
    MBool found = false;
    switch(m_allBoundaryIds[i]) {
      case 2000: {
        prepareBC2000();
        found = true;
        break;
      }
      default: {
      }
    }
    if(found) {
      break;
    }
  }

  RECORD_TIMER_STOP(t_setBCNC);
}

/** \brief Prepares the BC 4072
 *
 * \author Andreas Lintermann
 * \date 29.05.2015
 *
 * Determines the number of cells that belong to this BC and then communicates this
 * to all other domains. Those domains having such a BC are then grouped into a new
 * communicator m_BCComm and a file is opened to record the BC residual (the density
 * and the local Reynolds number.
 *
 * \param[in] bndCndIdIndex the index of the BC 4072 in the BC array
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::prepareBC2000() {
  TRACE();

  if(!m_calcWallForces) return;

  MInt noAllBoundaryIds = 0;
  const MInt* const allBoundaryIds = m_solver->m_geometry->GetBoundaryIds(&noAllBoundaryIds);
  TERMM_IF_COND(noAllBoundaryIds != m_noSegments, "noAllBoundaryIds != m_noSegments");

  m_mapWallForceContainer.clear();
  for(MInt i = 0; i < m_noSegments; i++) {
    if(allBoundaryIds[i] == 2000) {
      m_mapWallForceContainer[i];
    }
  }

  for(auto& mapWallForceIterator : m_mapWallForceContainer) {
    const MInt segId = mapWallForceIterator.first;
    const MInt index = m_mapBndCndSegId2Index[segId];
    auto& cwfc = mapWallForceIterator.second;

    const MInt noCalcForceCells = (index > -1) ? (m_bndCndOffsets[index + 1] - m_bndCndOffsets[index]) : 0;

    // unfortunately a global communication is requried
    std::vector<MInt> noCalcForceCellsPerDomain(m_solver->noDomains(), 0);
    MPI_Allgather(&noCalcForceCells, 1, MPI_INT, noCalcForceCellsPerDomain.data(), 1, MPI_INT, m_solver->mpiComm(), AT_,
                  "noCalcForceCells", "noCalcForceCellsPerDomain");
    std::vector<MInt> bcWallNeighbors;
    for(MInt i = 0; i < m_solver->noDomains(); i++) {
      if(noCalcForceCellsPerDomain[i] > 0) {
        bcWallNeighbors.push_back(i);
      }
    }

    cwfc.noComm = bcWallNeighbors.size();

    // build a new communicator
    if(cwfc.noComm != 0) {
      MPI_Group tmp_bcGroupWall;
      MPI_Group bcGroupWall;
      MPI_Comm_group(m_solver->mpiComm(), &tmp_bcGroupWall, AT_, "tmp_group");
      MPI_Group_incl(tmp_bcGroupWall, cwfc.noComm, bcWallNeighbors.data(), &bcGroupWall, AT_);
      MPI_Comm_create(m_solver->mpiComm(), bcGroupWall, &cwfc.comm, AT_, "cwfc.comm");

      cwfc.isRoot = (m_solver->domainId() == bcWallNeighbors[0]);
      if(cwfc.isRoot) {
        if(m_mapWallForceContainer.size() > 1) {
          std::stringstream fileNameS;
          fileNameS << m_forceFile << m_bndCndSegIds[index];
          cwfc.fileName = fileNameS.str();
        } else {
          cwfc.fileName = m_forceFile;
        }
        // Write header
        std::FILE* forceFile;
        forceFile = fopen(cwfc.fileName.c_str(), "a+");
        fprintf(forceFile, "# Re=%f, Ma=%f, refLength_LB=%f, dx(maxLvl)=%e\n", m_solver->m_Re, m_solver->m_Ma,
                m_solver->m_referenceLength, m_solver->c_cellLengthAtLevel(m_solver->maxLevel()));
        fprintf(forceFile, "#");
        fprintf(forceFile, "%s\t", "1:TS");
        constexpr MInt columnLength = 15;
        fprintf(forceFile, "%*s\t", columnLength, "2:F_x");
        fprintf(forceFile, "%*s\t", columnLength, "3:F_y");
        if constexpr(nDim == 3) {
          fprintf(forceFile, "%*s\t", columnLength, "4:F_z");
        }
        fprintf(forceFile, "\n");
        fclose(forceFile);
      }
    }
  }
}


/** \brief Checks if a BC exists that requires communication
 *
 * \author Andreas Lintermann
 * \date 29.05.2015
 *
 * \return a MBool that defines if such a BC exists or not
 *
 **/
template <MInt nDim>
MInt LbBndCnd<nDim>::checkForCommBC() {
  TRACE();

  if(m_solver->m_geometry->m_parallelGeometry) {
    set<MInt> tmplist;
    for(MInt i = 0; i < m_noSegments; i++)
      tmplist.insert(*(m_solver->m_geometry->geometryContext().getProperty("BC", i)->asInt()));
    m_noAllBoundaryIds = tmplist.size();
    mAlloc(m_allBoundaryIds, m_noAllBoundaryIds, "m_allBoundaryIds", 0, AT_);
    MInt pos = 0;
    for(set<MInt>::iterator it = tmplist.begin(); it != tmplist.end(); ++it, pos++)
      m_allBoundaryIds[pos] = *it;
  } else {
    m_noAllBoundaryIds = 0;
    m_allBoundaryIds = m_solver->m_geometry->GetBoundaryIds(&m_noAllBoundaryIds);
  }

  for(MInt i = 0; i < m_noAllBoundaryIds; i++) {
    switch(m_allBoundaryIds[i]) {
      case 1000: {
        m_numberOfCommBCs++;
        break;
      }
      case 1022: {
        m_numberOfCommBCs++;
        break;
      }
      case 1060: {
        m_numberOfCommBCs++;
        break;
      }
      case 1080: {
        m_numberOfCommBCs++;
        break;
      }
      case 4000: {
        m_numberOfCommBCs++;
        break;
      }
      case 4030: {
        m_numberOfCommBCs++;
        break;
      }
      case 4130: {
        m_numberOfCommBCs++;
        break;
      }
      case 4070: {
        m_numberOfCommBCs++;
        break;
      }
      case 4071: {
        m_numberOfCommBCs++;
        break;
      }
      case 4072: {
        m_numberOfCommBCs++;
        break;
      }
      case 4073: {
        m_numberOfCommBCs++;
        break;
      }
      case 4080: {
        m_numberOfCommBCs++;
        break;
      }
      case 4081: {
        m_numberOfCommBCs++;
        break;
      }
      case 4082: {
        m_numberOfCommBCs++;
        break;
      }
      case 4110: {
        m_numberOfCommBCs++;
        break;
      }
      default: {
      }
    }
  }
  return m_numberOfCommBCs;
}

/** \brief Prepares the BC 4070, 4071, 4072, 4080, 4081, and 4082
 *
 * \author Andreas Lintermann
 * \date 29.05.2015
 *
 * Determines the number of cells that belong to this BC and then communicates this
 * to all other domains. Those domains having such a BC are then grouped into a new
 * communicator m_BCComm and a file is opened to record the BC residual (the density
 * and the local Reynolds number.
 *
 * \param[in] index the index of the BC, counter the counter of BCs requiring communication
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::prepareBC(MInt index, MInt BCCounter, MInt segId) {
  TRACE();

  for(MInt count = 0; count < (MInt)m_bndCndSegIds.size(); count++) {
    if(m_bndCndIds[count] == index && m_bndCndSegIds[count] == segId) {
      const MInt bndCndIdIndex = count;
      m_allDomainsHaveBC[BCCounter][m_solver->domainId()] =
          m_bndCndOffsets[bndCndIdIndex + 1] - m_bndCndOffsets[bndCndIdIndex];
      m_mapBndCndIdSegId[count] = BCCounter;
      break;
    } else
      m_allDomainsHaveBC[BCCounter][m_solver->domainId()] = 0;
  }

  if(m_solver->noDomains() > 1) {
    // unfortunately a global communication is requried
    const MInt sndBuf = m_allDomainsHaveBC[BCCounter][m_solver->domainId()];
    MIntScratchSpace domainsHaveBC(m_solver->noDomains(), AT_, "domainsHaveBC");
    MPI_Allgather(&sndBuf, 1, MPI_INT, domainsHaveBC.getPointer(), 1, MPI_INT, m_solver->mpiComm(), AT_, "sndBuf",
                  "m_allDomainsHaveBC");

    for(MInt count = 0; count < m_solver->noDomains(); count++) {
      m_allDomainsHaveBC[BCCounter][count] = domainsHaveBC(count);
    }

    std::vector<MInt> BCneighborsPerSeg;
    m_totalNoBcCells[BCCounter] = 0;
    for(MInt i = 0; i < m_solver->noDomains(); i++) {
      if(m_allDomainsHaveBC[BCCounter][i] > 0) {
        BCneighborsPerSeg.push_back(i);
        m_totalNoBcCells[BCCounter] += m_allDomainsHaveBC[BCCounter][i];
      }
    }

    m_BCneighbors.push_back(BCneighborsPerSeg);

    m_noBCNeighbors[BCCounter] = (MInt)m_BCneighbors[BCCounter].size();

    if(m_noBCNeighbors[BCCounter] > 1) {
      MIntScratchSpace ptBCneighbors(m_noBCNeighbors[BCCounter], AT_, "ptBCneighbors");
      for(MInt i = 0; i < m_noBCNeighbors[BCCounter]; i++)
        ptBCneighbors.p[i] = m_BCneighbors[BCCounter][i];

      MPI_Comm_group(m_solver->mpiComm(), &tmp_group[BCCounter], AT_, "tmp_group");

      MPI_Group_incl(tmp_group[BCCounter], m_noBCNeighbors[BCCounter], ptBCneighbors.getPointer(), &BCGroup[BCCounter],
                     AT_);
      MPI_Comm_create(m_solver->mpiComm(), BCGroup[BCCounter], &m_BCComm[BCCounter], AT_, "m_BCComm");
    }
  } else {
    std::vector<MInt> BCneighborsPerSeg;
    BCneighborsPerSeg.push_back(m_solver->domainId());
    m_BCneighbors.push_back(BCneighborsPerSeg);
    m_totalNoBcCells[BCCounter] = m_allDomainsHaveBC[BCCounter][0];
    m_noBCNeighbors[BCCounter] = 1;
  }
}

template <MInt nDim>
void LbBndCnd<nDim>::prepareBC4073(MInt BCCounter, MInt segId) {
  m_log << "    - BC 4073" << endl;
  m_log << "      * reading properties from file" << endl;
  // get the plane from the geomtry file that specifies the plane that should be used for
  // the calculation of the local Reynolds number
  if(m_solver->m_geometry->geometryContext().propertyExists("localReCut", 0)) {
    MInt num = m_solver->m_geometry->geometryContext().getProperty("localReCut", 0)->count();
    if((nDim == 2 && num != 4) || (nDim == 3 && num != 6)) {
      stringstream errorMsg;
      errorMsg << "ERROR: no wrong number of entries in the geometry property 'localReCut' for BC 4073" << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    // this holds the cut information for the local Reynolds number calculation
    mAlloc(m_localReCutPoint, nDim, "m_localReCutPoint", 0.0, AT_);
    mAlloc(m_localReCutNormal, nDim, "m_localReCutNormal", 0.0, AT_);

    MFloat len = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      m_localReCutPoint[i] = *m_solver->m_geometry->geometryContext().getProperty("localReCut", 0)->asFloat(i);
      m_localReCutNormal[i] = *m_solver->m_geometry->geometryContext().getProperty("localReCut", 0)->asFloat(i + nDim);
      len += m_localReCutNormal[i] * m_localReCutNormal[i];
    }
    len = sqrt(len);
    if(approx(len, 0.0, MFloatEps)) {
      stringstream errorMsg;
      errorMsg << "ERROR: the normal defined by the geometry property 'localReCut' seems to be wrong" << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    // normalize
    for(MInt i = 0; i < nDim; i++)
      m_localReCutNormal[i] /= len;

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCutDistance", 0))
      m_localReCutDistance = *m_solver->m_geometry->geometryContext().getProperty("localReCutDistance", 0)->asFloat();
    else
      m_localReCutDistance = m_solver->c_cellLengthAtLevel(0);

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCutInterval", 0))
      m_localReCutInterval = *m_solver->m_geometry->geometryContext().getProperty("localReCutInterval", 0)->asInt();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCutInterval' defined for BC 4073" << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }
  } else {
    stringstream errorMsg;
    errorMsg << "ERROR: no geometry property 'localReCut' defined for BC 4073" << endl;
    m_log << errorMsg.str();
    TERMM(1, errorMsg.str());
  }

  m_localReCutReportInterval = 1000;
  if(m_solver->m_geometry->geometryContext().propertyExists("localReCutReportInterval", 0))
    m_localReCutReportInterval =
        *m_solver->m_geometry->geometryContext().getProperty("localReCutReportInterval", 0)->asInt();

  m_localReCutRe = m_solver->m_Re;
  if(m_solver->m_geometry->geometryContext().propertyExists("localReCutRe", 0))
    m_localReCutRe = *m_solver->m_geometry->geometryContext().getProperty("localReCutRe", 0)->asFloat();

  m_localReCutDiameter = m_referenceLength;
  if(m_solver->m_geometry->geometryContext().propertyExists("localReCutDiameter", 0))
    m_localReCutDiameter = *m_solver->m_geometry->geometryContext().getProperty("localReCutDiameter", 0)->asFloat();

  m_localReCutAdpPerc = 0.2;
  if(m_solver->m_geometry->geometryContext().propertyExists("localReCutAdpPerc", 0))
    m_localReCutAdpPerc = *m_solver->m_geometry->geometryContext().getProperty("localReCutAdpPerc", 0)->asFloat();

  // load from restart? -> overwrite local data
  m_log << "      * reading information for restart" << endl;
  if(m_solver->m_restartFile) {
    if(m_solver->m_geometry->geometryContext().propertyExists("localReCut_rho1", 0))
      m_solver->m_rho1 = *m_solver->m_geometry->geometryContext().getProperty("localReCut_rho1", 0)->asFloat();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCut_rho1' defined for BC 4073 (as required for restart)" << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCut_lRho", 0))
      m_lRho = *m_solver->m_geometry->geometryContext().getProperty("localReCut_lRho", 0)->asFloat();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCut_lRho' defined for BC 4073 (as required for restart)" << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCut_rhoLast", 0))
      m_rhoLast = *m_solver->m_geometry->geometryContext().getProperty("localReCut_rhoLast", 0)->asFloat();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCut_rhoLast' defined for BC 4073 (as required for restart)"
               << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCut_deltaRho", 0))
      m_deltaRho = *m_solver->m_geometry->geometryContext().getProperty("localReCut_deltaRho", 0)->asFloat();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCut_deltaRho' defined for BC 4073 (as required for restart)"
               << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }

    if(m_solver->m_geometry->geometryContext().propertyExists("localReCut_ReLast", 0))
      m_ReLast = *m_solver->m_geometry->geometryContext().getProperty("localReCut_ReLast", 0)->asFloat();
    else {
      stringstream errorMsg;
      errorMsg << "ERROR: no geometry property 'localReCut_ReLast' defined for BC 4073 (as required for restart)"
               << endl;
      m_log << errorMsg.str();
      TERMM(1, errorMsg.str());
    }
  }

  // find all cells that have a cut with the line/plane and the according distance
  m_log << "      * finding cells that have a cut with the reference plane" << endl;
  for(MInt i = 0; i < m_solver->m_cells.size(); i++) {
    // only leaf cells
    if(m_solver->c_noChildren(i) > 0 || m_solver->a_isHalo(i)) continue;

    MFloat cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(i) + 1);
    const MFloat* const coordinates = &(m_solver->a_coordinate(i, 0));

    std::array<MFloat, nDim> edge;

    // fill edge and projection dist = edge * normal (which is normalized)
    MFloat prodist = 0.0;
    MFloat dist = 0.0;
    for(MInt d = 0; d < nDim; d++) {
      edge[d] = coordinates[d] - m_localReCutPoint[d];
      prodist += edge[d] * m_localReCutNormal[d];
      dist += edge[d] * edge[d];
    }

    dist = sqrt(dist);

    if(fabs(prodist) < cellHalfLength && dist <= m_localReCutDistance) m_localReCutCells.push_back(i);
  }
  m_log << "        = no of cells for this domain: " << m_localReCutCells.size() << endl;
  m_hasLocalReCut = (m_localReCutCells.size() > 0);

  // does this domain have the according BC?
  MBool has_bc4073 = false;
  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++)
    if(m_bndCndIds[i] == 4073 && m_mapSegIdsInOutCnd[i] != -1 && m_bndCndSegIds[i] == segId) {
      has_bc4073 = true;
      break;
    }

  m_allDomainsHaveBC[BCCounter][m_solver->domainId()] = m_hasLocalReCut ? m_localReCutCells.size() : -1;

  if(m_allDomainsHaveBC[BCCounter][m_solver->domainId()] == -1 && has_bc4073)
    m_allDomainsHaveBC[BCCounter][m_solver->domainId()] = 0;

  // unfortunately a global communication is requried
  MInt* const sndBuf = &(m_allDomainsHaveBC[BCCounter][m_solver->domainId()]);
  MPI_Allgather(sndBuf, 1, MPI_INT, m_allDomainsHaveBC[BCCounter], 1, MPI_INT, m_solver->mpiComm(), AT_, "sndBuf",
                "m_allDomainsHaveBC");

  std::vector<MInt> BCneighborsPerSeg;

  m_totalNoDomainsReCut = 0;
  m_totalNoBcCells[BCCounter] = 0;
  for(MInt i = 0; i < m_solver->noDomains(); i++)
    if(m_allDomainsHaveBC[BCCounter][i] != -1) {
      BCneighborsPerSeg.push_back(i);
      if(m_allDomainsHaveBC[BCCounter][i] > 0) {
        m_totalNoDomainsReCut++;
        m_totalNoBcCells[BCCounter] += m_allDomainsHaveBC[BCCounter][i];
      }
    }

  m_BCneighbors.push_back(BCneighborsPerSeg);
  m_noBCNeighbors[BCCounter] = (MInt)m_BCneighbors[BCCounter].size();

  m_firstBCinComm = 0;
  if(m_noBCNeighbors[BCCounter] > 0) m_firstBCinComm = m_BCneighbors[BCCounter][0];

  MIntScratchSpace ptBCneighbors(m_noBCNeighbors[BCCounter], AT_, "ptBCneighbors");
  for(MInt i = 0; i < m_noBCNeighbors[BCCounter]; i++)
    ptBCneighbors.p[i] = m_BCneighbors[BCCounter][i];

  // build a new communicator
  m_log << "      * building MPI communicator" << endl;
  if(m_noBCNeighbors[BCCounter] != 0 && m_solver->noDomains() > 1) {
    //      MPI_Group tmp_group, BCGroup;

    MPI_Comm_group(m_solver->mpiComm(), &tmp_group[BCCounter], AT_, "tmp_group");

    MPI_Group_incl(tmp_group[BCCounter], m_noBCNeighbors[BCCounter], ptBCneighbors.getPointer(), &BCGroup[BCCounter],
                   AT_);
#ifdef MAIA_MPI_DEBUG
    MPI_Comm_create(m_solver->mpiComm(), BCGroup[BCCounter], &m_BCComm[BCCounter], AT_, "m_BCComm");
#else
    MPI_Comm_create(m_solver->mpiComm(), BCGroup[BCCounter], &m_BCComm[BCCounter], AT_, "m_BCComm");
#endif
  }
  if(m_calcBcResidual && m_solver->domainId() == m_firstBCinComm) {
    m_BCResidualStream[BCCounter].open(m_BCOutputFileName[BCCounter], ios_base::app);

    if(!m_solver->m_restartFile)
      m_BCResidualStream[BCCounter] << "############################\n"
                                    << "# Order of appearance:\n"
                                    << "#   1: globalTimeStep\n"
                                    << "#   2: m_localReCutRe\n"
                                    << "#   3: l_Re\n"
                                    << "#   4: m_ReLast\n"
                                    << "#   5: Re_diff\n"
                                    << "#   6: m_lRho\n"
                                    << "#   7: m_rhoLast\n"
                                    << "#   8: m_deltaRho\n"
                                    << "############################\n"
                                    << endl;
  }

  // special treatment for those domains that have a cut with the plane/line but do not contain any BC 4073 cells
  m_log << "      * updating boundary condition list" << endl;
  if(!has_bc4073 && m_hasLocalReCut) {
    // we need to add bc4073
    m_bndCndIds.push_back(4073);
    m_bndCndOffsets.push_back(m_bndCndOffsets[m_bndCndOffsets.size() - 1]);
    // empty treatment
    m_bndCndSegIds.push_back(-1);
    m_mapSegIdsInOutCnd.push_back(-1);

    m_log << "   + created new BC 4073" << endl;
  }

  m_log << "      * update interval:     " << m_localReCutInterval << endl;
  m_log << "      * update percentage:   " << m_localReCutAdpPerc << endl;
  m_log << "      * report interval:     " << m_localReCutReportInterval << endl;
  m_log << "      * point:               ";
  for(MInt i = 0; i < nDim; i++)
    m_log << m_localReCutPoint[i] << " ";
  m_log << endl;
  m_log << "      * distance:            " << m_localReCutDistance << endl;
  m_log << "      * normal:              ";
  for(MInt i = 0; i < nDim; i++)
    m_log << m_localReCutNormal[i] << " ";
  m_log << endl;
  m_log << "      * target Re:           " << m_localReCutRe
        << (approx(m_localReCutRe, m_solver->m_Re, MFloatEps) ? " (same as global Re)" : " (different than global Re)")
        << endl;
  m_log << "      * reference length     " << m_localReCutDiameter
        << (approx(m_localReCutDiameter, m_referenceLength, MFloatEps) ? " (same as global ref. length)"
                                                                       : " (different than global ref. length)")
        << endl;
  m_log << "      * local cut:           " << (m_hasLocalReCut ? "yes" : "no") << endl;
  m_log << "      * no. cut cells:       " << m_localReCutCells.size() << endl;
  m_log << "      * no. total cut cells: " << m_totalNoBcCells[BCCounter] << endl;
  m_log << "      * domain owns BC:      " << (has_bc4073 ? "yes" : "no") << endl;
  m_log << "      * no. domains cut:     " << m_totalNoDomainsReCut << endl;
  m_log << "      * no. comm. domains:   " << m_noBCNeighbors[BCCounter] << endl;
  m_log << "      * first comm. domain:  " << m_firstBCinComm << endl;
  m_log << "      * comm. domains:       ";
  for(MInt i = 0; i < m_noBCNeighbors[BCCounter]; i++)
    m_log << m_BCneighbors[BCCounter][i] << " ";
  m_log << endl;
  if(m_solver->m_restartFile) {
    m_log << "      * restart details:     " << endl;
    m_log << "        > rho1:     " << m_solver->m_rho1 << endl;
    m_log << "        > lRho:     " << m_lRho << endl;
    m_log << "        > rhoLast:  " << m_rhoLast << endl;
    m_log << "        > deltaRho: " << m_deltaRho << endl;
    m_log << "        > ReLast:   " << m_ReLast << endl;
  }
  m_log << endl;
}


/** \brief Sets up a neighbor-communicator for certain BCs
 *
 * \author Andreas Lintermann
 * \date 27.09.2012
 *
 * This function checks if a BC is used which requires internal communication.
 * If so, a new communicator group is created holding only the processes having a
 * part of the BC. A file is opened, which will hold the residual for this BC.
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::setBCNeighborCommunicator() {
  TRACE();

  NEW_SUB_TIMER(t_setBCNC, "set BC neighbor communicators", m_t_BCAll);
  RECORD_TIMER_START(t_setBCNC);

  // first of all check if we have a BC which requries communication
  m_numberOfCommBCs = checkForCommBC();
  if(m_numberOfCommBCs == 0) {
    RECORD_TIMER_STOP(t_setBCNC);
    return;
  }

  mAlloc(m_allDomainsHaveBC, m_numberOfCommBCs, m_solver->noDomains(), "m_allDomainsHaveBC", 0, AT_);
  mAlloc(m_noBCNeighbors, m_numberOfCommBCs, "m_noBCNeighbors", 0, AT_);
  mAlloc(m_totalNoBcCells, m_numberOfCommBCs, "m_totalNoBcCells", 0, AT_);

  mAlloc(tmp_group, m_numberOfCommBCs, "tmp_group", AT_);
  mAlloc(BCGroup, m_numberOfCommBCs, "BCGroup", AT_);
  mAlloc(m_BCComm, m_numberOfCommBCs, "m_BCComm", AT_);

  if(m_calcBcResidual) mAlloc(m_BCResidualStream, m_numberOfCommBCs, "m_BCResidualStream", AT_);
  mAlloc(m_BCOutputFileName, m_numberOfCommBCs, "m_BCOutputFileName", AT_);

  if((MInt)(m_bndCndIds.size()) > 0)
    mAlloc(m_mapBndCndIdSegId, (MInt)(m_bndCndIds.size()), "m_mapBndCndIdSegId", 0, AT_);

  m_log << "  + Found a BC which requires communication:" << endl;

  MInt counterCommBC = 0;

  for(MInt segId = 0; segId < m_noAllBoundaryIds; segId++) {
    stringstream s;
    MBool found = false;
    switch(m_allBoundaryIds[segId]) {
      case 1000: {
        s << "Output_SegNo_" << segId << "_BC_" << 1000 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(1000, counterCommBC, segId);
        found = true;
        break;
      }
      case 1022: {
        s << "Output_SegNo_" << segId << "_BC_" << 1022 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(1022, counterCommBC, segId);
        found = true;
        break;
      }
      case 1060: {
        s << "Output_SegNo_" << segId << "_BC_" << 1060 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(1060, counterCommBC, segId);
        found = true;
        break;
      }
      case 1080: {
        s << "Output_SegNo_" << segId << "_BC_" << 1080 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(1080, counterCommBC, segId);
        found = true;
        break;
      }
      case 4000: {
        s << "Output_SegNo_" << segId << "_BC_" << 4000 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4000, counterCommBC, segId);
        found = true;
        break;
      }
      case 4030: {
        s << "Output_SegNo_" << segId << "_BC_" << 4030 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4030, counterCommBC, segId);
        found = true;
        break;
      }
      case 4130: {
        s << "Output_SegNo_" << segId << "_BC_" << 4130 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4130, counterCommBC, segId);
        found = true;
        break;
      }
      case 4070: {
        s << "Output_SegNo_" << segId << "_BC_" << 4070 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4070, counterCommBC, segId);
        found = true;
        break;
      }
      case 4071: {
        s << "Output_SegNo_" << segId << "_BC_" << 4071 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4071, counterCommBC, segId);
        found = true;
        break;
      }
      case 4072: {
        s << "Output_SegNo_" << segId << "_BC_" << 4072 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4072, counterCommBC, segId);
        found = true;
        break;
      }
      case 4073: {
        s << "Output_SegNo_" << segId << "_BC_" << 4073 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC4073(counterCommBC, segId);
        found = true;
        break;
      }
      case 4080: {
        s << "Output_SegNo_" << segId << "_BC_" << 4080 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4080, counterCommBC, segId);
        found = true;
        break;
      }
      case 4081: {
        s << "Output_SegNo_" << segId << "_BC_" << 4081 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4081, counterCommBC, segId);
        found = true;
        break;
      }
      case 4082: {
        s << "Output_SegNo_" << segId << "_BC_" << 4082 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4082, counterCommBC, segId);
        found = true;
        break;
      }
      case 4110: {
        s << "Output_SegNo_" << segId << "_BC_" << 4110 << ".dat";
        m_BCOutputFileName[counterCommBC] = s.str();
        prepareBC(4110, counterCommBC, segId);
        found = true;
        break;
      }
      default: {
      }
    }

    if(found) {
      counterCommBC++;
    }
  }

  RECORD_TIMER_STOP(t_setBCNC);
}

//! Dereferences bndCndHandlerVariables
template <MInt nDim>
void LbBndCnd<nDim>::updateVariables() {
  TRACE();

  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
    (this->*bndCndHandlerVariables[i])(i);
  }
}


//! Dereferences bndCndHandlerRHS
template <MInt nDim>
void LbBndCnd<nDim>::updateRHS() {
  TRACE();
  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
    (this->*bndCndHandlerRHS[i])(i);
  }
}

/** \brief Creates boundary cells according to the geometry information
    \author Andreas Lintermann
    \date 29.01.2010

    This function runs over all cells and makes an intersection test with the triangles of nearby segments.
    All candidates carry a segment id and a boundary condition. Hence, a corner cell can have multiple properties.
    Depending on the allowance of multiple BCs, only one single (no multiple BCs allowed) boundary cell is added to
    a temporary vector. Anyhow, the boundary cell can carry multiple information on segment id and boundary condition.
    In this case the wall boundary condition is chosen by setting the first element of both the boundary condition array
    and the segment id to this wall condition. In the other case of multiple BCs, the boundary cell is added multiple
    times to the vector of temporary boundary cells, each time with its according different conditions. After that,
    the collector of boundary cells is initialized. The size is now known by the size of the temporary vector. Finally,
    the information in the vector is copied to the collector.
  */
template <MInt nDim>
void LbBndCnd<nDim>::createBoundaryCells() {
  TRACE();

  NEW_SUB_TIMER(t_createBCCells, "create boundary cells", m_t_BCAll);
  RECORD_TIMER_START(t_createBCCells);

  MInt bndCellsIt = 0, bndCellsIt2 = 0;
  MFloat cellHalfLength = 0.0;
  ScratchSpace<MFloat> target(2 * nDim, AT_, "target");

  m_log << "  + Creating boundary cells..." << endl;
  m_log << "    - Multiple BC treatment: " << m_multiBCTreatment << endl;

  // This temporary, has to be filled into m_bndCells later on
  for(MInt i = 0; i < m_solver->m_cells.size(); i++) {
    // skip inActive cells (neither leaf nor interfaceParent)
    if(!m_solver->c_isLeafCell(i) && !m_solver->a_isInterfaceParent(i)) continue;

    // Define corners of current cell in target
    for(MInt j = 0; j < nDim; j++) {
      cellHalfLength = m_solver->c_cellLengthAtLevel(m_solver->a_level(i) + 1);
      target[j] = m_solver->a_coordinate(i, j) - cellHalfLength;
      target[j + nDim] = m_solver->a_coordinate(i, j) + cellHalfLength;
    }

    // Check for intersection with geometry elements
    MFloat* targetPtr = &target[0];
    std::vector<MInt> nodeList;
    if(m_gridCutTest == "SAT")
      m_solver->m_geometry->getIntersectionElements(targetPtr, nodeList, cellHalfLength, &m_solver->a_coordinate(i, 0));
    else
      m_solver->m_geometry->getIntersectionElements(targetPtr, nodeList);
    const MInt noNodes = nodeList.size();

    if(noNodes > 0) {
      m_solver->a_onlyBoundary(i) = true;
      m_solver->a_isBndryCell(i) = true;

      // Create a new boundary cell for our found cuts
      m_bndCells.emplace_back();

      bndCellsIt = m_bndCells.size() - 1;
      m_bndCells[bndCellsIt].m_cellId = i;

      // Fill the segmentIds of the boundary cell
      for(MInt j = 0; j < noNodes; j++) {
        MBool already_in = false;
        for(MInt k = 0; k < (MInt)(m_bndCells[bndCellsIt].m_segmentId.size()); k++)
          if(m_bndCells[bndCellsIt].m_segmentId[k] == m_solver->m_geometry->elements[nodeList[j]].m_segmentId) {
            already_in = true;
            break;
          }
        if(!already_in) {
          m_bndCells[bndCellsIt].m_segmentId.push_back(m_solver->m_geometry->elements[nodeList[j]].m_segmentId);
          m_bndCells[bndCellsIt].m_bndCndId.push_back(m_solver->m_geometry->elements[nodeList[j]].m_bndCndId);
        }
      }

      if(m_bndCells[bndCellsIt].m_segmentId.size() > 1) {
        // exists inout / wall / periodic ?
        MInt positions[3] = {-1, -1, -1};

        for(MInt j = 0; j < (MInt)(m_bndCells[bndCellsIt].m_segmentId.size()); j++) {
          MBool is_inout = false;
          MBool is_periodic = false;

          for(MInt k = 0; k < m_noInOutSegments; k++)
            if(m_bndCells[bndCellsIt].m_segmentId[j] == m_inOutSegmentsIds[k]) {
              is_inout = true;
              break;
            }
          if(m_noPeriodicSegments != 0)
            for(MInt k = 0; k < m_noPeriodicSegments; k++)
              if(m_bndCells[bndCellsIt].m_segmentId[j] == m_periodicSegmentsIds[k]) {
                is_periodic = true;
                break;
              }
          if(is_inout)
            positions[0] = j;
          else if(!is_inout && !is_periodic)
            positions[1] = j;
          else if(is_periodic)
            positions[2] = j;
        }

        MInt swap_pos = 0;
        if(m_multiBCTreatment == "W-P-I") {
          // Make wall first (if exists)
          if(positions[1] != -1) swap_pos = positions[1];
          // Make periodic first (if exists)
          else if(positions[2] != -1)
            swap_pos = positions[2];

        } else if(m_multiBCTreatment == "W-I-P") {
          // Make wall first (if exists)
          if(positions[1] != -1) swap_pos = positions[1];
          // Make inout first (if exists)
          else if(positions[0] != -1)
            swap_pos = positions[0];
        } else if(m_multiBCTreatment == "I-W-P") {
          // Make inout first (if exists)
          if(positions[0] != -1) swap_pos = positions[0];
          // Make wall first (if exists)
          else if(positions[1] != -1)
            swap_pos = positions[1];
        } else if(m_multiBCTreatment == "I-P-W") {
          // Make inout first (if exists)
          if(positions[0] != -1) swap_pos = positions[0];
          // Make periodic first (if exists)
          else if(positions[2] != -1)
            swap_pos = positions[2];
        } else if(m_multiBCTreatment == "P-W-I") {
          // Make periodic first (if exists)
          if(positions[2] != -1) swap_pos = positions[2];
          // Make wall first (if exists)
          else if(positions[1] != -1)
            swap_pos = positions[1];
        } else if(m_multiBCTreatment == "P-I-W") {
          // Make periodic first (if exists)
          if(positions[2] != -1) swap_pos = positions[2];
          // Make inout first (if exists)
          else if(positions[0] != -1)
            swap_pos = positions[0];
        }

        /********************************/
        /* Georg 24.06.2010             */
        /*                              */
        /********************************/
        else if(m_multiBCTreatment == "multiple") {
          // Add the bndCell n more times according to the number of different bc's.
          // No sorting of segmentIds is applied.

          bndCellsIt2 = bndCellsIt;
          for(MInt j = 1; j < (MInt)(m_bndCells[bndCellsIt].m_segmentId.size()); j++) {
            // make sure that m_bndCells are only added if there is another bc
            MBool already_in = false;
            for(MInt k = 0; k <= (bndCellsIt2 - bndCellsIt); k++) {
              if(m_bndCells[bndCellsIt].m_bndCndId[bndCellsIt2 - bndCellsIt] == m_bndCells[bndCellsIt].m_bndCndId[j]) {
                already_in = true;
                break;
              }
            }

            if(!already_in) {
              bndCellsIt2++;
              m_bndCells.emplace_back();            // add new m_bndCells
              m_bndCells[bndCellsIt2].m_cellId = i; // new m_bndCells has the same cellId as the previous one

              m_bndCells[bndCellsIt2].m_segmentId.push_back(m_bndCells[bndCellsIt].m_segmentId[j]);
              m_bndCells[bndCellsIt2].m_bndCndId.push_back(m_bndCells[bndCellsIt].m_bndCndId[j]);

              // 			    m_bndCells[bndCellsIt2].m_segmentId[0] = m_bndCells[bndCellsIt].m_segmentId[j];
              // 			    m_bndCells[bndCellsIt2].m_bndCndId[0] = m_bndCells[bndCellsIt].m_bndCndId[j];
            }
          }
        }

        MInt tmp = m_bndCells[bndCellsIt].m_segmentId[0];
        m_bndCells[bndCellsIt].m_segmentId[0] = m_bndCells[bndCellsIt].m_segmentId[swap_pos];
        m_bndCells[bndCellsIt].m_segmentId[swap_pos] = tmp;

        tmp = m_bndCells[bndCellsIt].m_bndCndId[0];
        m_bndCells[bndCellsIt].m_bndCndId[0] = m_bndCells[bndCellsIt].m_bndCndId[swap_pos];
        m_bndCells[bndCellsIt].m_bndCndId[swap_pos] = tmp;
      }
    }
  }

  m_log << "    - noBndCells:            " << m_bndCells.size() << endl << endl;

  for(auto& bndCell : m_bndCells) {
    bndCell.m_multiplier = 1.0;
  }

  RECORD_TIMER_STOP(t_createBCCells);
}

/** \brief Initialize BC data variables
 *  \author Miro Gondrum
 *  \date 09.02.2021
 */
template <MInt nDim>
void LbBndCnd<nDim>::initializeBcData() {
  // Initialize BC
  for(MInt i = 0; i < (MInt)(m_bndCndIds.size()); i++) {
    (this->*bndCndInitiator[i])(i);
  }
  if(!m_solver->m_restartFile) {
    for(auto p : m_bndCndData) {
      const MInt& p_index = p.first;
      const LbBndCndData& p_data = p.second;
      const MInt bndCndOffset = m_bndCndOffsets[p_index];
      for(MInt i = 0; i < p_data.noBndCellsWithHalos; i++) {
        for(MInt j = 0; j < p_data.noVars; j++) {
          p_data.data[i * p_data.noVars + j] = m_solver->a_oldVariable(m_bndCells[i + bndCndOffset].m_cellId, j);
        }
      }
    }
  }
}

/** \brief This function creates the communicator for calculating the wall forces of the
 * level-set boundaries
 * \author Moritz Waldmann
 * \date 21.10.2019
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::createMBComm() {
  TRACE();

  if(!m_calcWallForces) return;

  if(!m_solver->grid().isActive()) {
    return;
  }

  if(m_solver->noDomains() > 1) {
    if(/*!firstInit &&*/ (m_BCWallMBComm != MPI_COMM_NULL)) {
      MPI_Comm_free(&m_BCWallMBComm, AT_, "m_BCWallMBComm");
    }

    m_allDomainsCalcForceMB[m_solver->domainId()] = m_solver->m_currentNoG0Cells;

    if(!m_BCWallMBNeighbors.empty()) m_BCWallMBNeighbors.clear();

    // unfortunately a global communication is requried
    MInt* const sndBuf = &(m_allDomainsCalcForceMB[m_solver->domainId()]);
    MPI_Allgather(sndBuf, 1, MPI_INT, m_allDomainsCalcForceMB, 1, MPI_INT, m_solver->mpiComm(), AT_, "sendBuf",
                  "m_allDomainsCalcForceMB");

    //    for(MInt count = 0; count < m_solver->noDomains(); count++) {
    //      if(m_allDomainsCalcForceMB[count] != 0)
    //        m_log << "    Domain[" << count << "] carries " << m_allDomainsCalcForceMB[count] << " MB cells. " <<
    //        endl;
    //    }
    for(MInt i = 0; i < m_solver->noDomains(); i++) {
      if(m_allDomainsCalcForceMB[i] > 0) {
        m_BCWallMBNeighbors.push_back(i);
      }
    }

    m_noBCWallMBNeighbors = (MInt)m_BCWallMBNeighbors.size();

    MIntScratchSpace ptBCWallMBNeighbors(m_noBCWallMBNeighbors, AT_, "ptBCWallMBNeighbors");
    for(MInt i = 0; i < m_noBCWallMBNeighbors; i++)
      ptBCWallMBNeighbors[i] = m_BCWallMBNeighbors[i];

    // build a new communicator
    if(m_noBCWallMBNeighbors != 0) {
      MPI_Group tmp_groupMB, BCGroupMB;
      MPI_Comm_group(m_solver->mpiComm(), &tmp_groupMB, AT_, "tmp_groupMB");
      MPI_Group_incl(tmp_groupMB, m_noBCWallMBNeighbors, ptBCWallMBNeighbors.getPointer(), &BCGroupMB, AT_);

      MPI_Comm_create(m_solver->mpiComm(), BCGroupMB, &m_BCWallMBComm, AT_, "m_BCWallMBComm");
      MPI_Group_free(&tmp_groupMB, AT_);
      MPI_Group_free(&BCGroupMB, AT_);
      if(m_solver->domainId() == m_BCWallMBNeighbors[0]) {
        m_forceFile = "forcesMB.log";
        // Write header
        std::FILE* forceFile;
        forceFile = fopen(m_forceFile.c_str(), "a+");
        fprintf(forceFile, "# Re=%f, Ma=%f, refLength_LB=%f, dx(maxLvl)=%e\n", m_solver->m_Re, m_solver->m_Ma,
                m_solver->m_referenceLength, m_solver->c_cellLengthAtLevel(m_solver->maxLevel()));
        fprintf(forceFile, "#");
        fprintf(forceFile, "%s\t", "1:timeStep");
        constexpr MInt columnLength = 15;
        fprintf(forceFile, "%*s\t", columnLength, "2:F_x");
        fprintf(forceFile, "%*s\t", columnLength, "3:F_y");
        if(nDim == 3) {
          fprintf(forceFile, "%*s\t", columnLength, "4:F_z");
        }
        fprintf(forceFile, "\n");
        fclose(forceFile);
      }
    }

  } else {
    m_noBCWallMBNeighbors = 1;
    m_BCWallMBNeighbors.push_back(m_solver->domainId());
  }
}


/** \brief This function initializes the LbBndCnd for coupled simulations
 * \author Moritz Waldmann
 * \date 21.10.2019
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::initializeBndMovingBoundaries() {
  TRACE();

  //#####################################################################
  // Allocate space for pointer lists which store the boundary functions
  //#####################################################################

  m_log << "  + Setting the boundary condition handler for LB_LS..." << endl << endl;
  bndCndHandlerVariables_MB.resize(1);
  bndCndHandlerVariables_MB[0] = &LbBndCnd::bc0;
  bndCndHandlerRHS_MB.resize(1);
  // TODO labels:LB,toenhance Generalize the configuration of mb bnd cnd
  // Use real boundary conditions for InitGFieldFromSTL stls ??
  bndCndHandlerRHS_MB[0] = &LbBndCnd::bc66666;
  /*if(m_solver->m_bernoulliBeam)
    bndCndHandlerRHS_MB[0] = &LbBndCnd::bc66667;*/
  // Free surface
  // bndCndHandlerRHS_MB[0] = &LbBndCnd::bc66668;

  //#####################################################################
  // Allocation of all lists needed
  // Setting the required constants
  //#####################################################################

  // Create Communicator
  mAlloc(m_allDomainsCalcForceMB, m_solver->noDomains(), "m_allDomainsCalcForceMB", 0, AT_);

  //#####################################################################
  // Create Communicator
  //#####################################################################

  m_BCWallMBComm = MPI_COMM_NULL;
}

/** \brief This function does the sub coupling step called from the coupling class
 * \author Moritz Waldmann
 * \date 21.10.2019
 *
 **/
template <MInt nDim>
void LbBndCnd<nDim>::postCouple() {
  if(!m_solver->grid().isActive()) {
    return;
  }

  const MInt startSet = m_solver->m_levelSetId;
  for(MInt set = startSet; set < m_solver->m_maxNoSets; set++) {
    (this->*bndCndHandlerVariables_MB[0])(set);
    (this->*bndCndHandlerRHS_MB[0])(set);
  }
}

// Explicit instantiations for 2D and 3D
template class LbBndCnd<2>;
template class LbBndCnd<3>;
