// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptcollision.h"
#include "COMM/mpioverride.h"
#include "IO/parallelio.h"
#include "UTIL/maiamath.h"
#include "lpt.h"
#include "lptellipsoidal.h"
#include "lptellipsoiddistance.h"
#include "lptspherical.h"

using namespace std;
using namespace maia::math;
using namespace maia::lpt;

/// \brief Constructor of the ParticleCollision class.
///
/// \author Sven Berger, Tim Wegmann
/// \date   June 2016
/// \param[in] collisionModel Collision model to use.
/// \param[in] lptSolver Pointer to LPT solver.
/// \param[in] fvSolver Pointer to matching FV solver.
template <MInt nDim>
ParticleCollision<nDim>::ParticleCollision(LPT<nDim>* lptSolver, MInt collisionModel, const MInt solverId,
                                           const MInt domainId, const MInt noDomains, const MPI_Comm comm)
  : m_collisionModel(collisionModel),
    m_lpt(lptSolver),
    m_solverId(solverId),
    m_domainId(domainId),
    m_noDomains(noDomains),
    m_mpiComm(comm) {
  init();
}

/// \fn void ParticleCollision<nDim>::init()
/// \brief Initialisation of collision model.
///
/// \author Sven Berger
/// \date   June 2016
template <MInt nDim>
void ParticleCollision<nDim>::init() {
  /*! \page propertyPageLPT
  \section particleCollOutputStep
  <code>MInt LPT::m_outputStep</code>\n
  default = <code>50</code> \n \n
  Number of time steps between output of collision data. \n
  Keywords: <i>PARTICLE</i>
   */
  m_outputStep = Context::getSolverProperty<MInt>("particleCollOutputStep", solverId(), AT_, &m_outputStep);

  /*! \page propertyPageLPT
      \section particleCollOffset
      <code>MFloat LPT::m_offset</code>\n
      default = <code>-1000.0</code> \n \n
      Specify a lower limit for the x coordinate of a collision to be \n
      detected. This property is useful to allow for an "relaxation \n
      distance" for newly introduced particles. \n
      Note that at the default value no limit is used. \n
      Keywords: <i>PARTICLE</i>
   */
  m_offset = Context::getSolverProperty<MFloat>("particleCollOffset", solverId(), AT_, &m_offset);

  /*! \page propertyPageLPT
  \section particleIncludeEllipsoids
  <code>MBool LPT::m_particleIncludeEllipsoids</code>\n
  default = false \n \n
  Exclude (false) or include (true) ellipsoidal particles \n
  in addition to spheres.\n
  Keywords: <i>PARTICLE</i>
  */
  MBool ellipsoids = false;
  ellipsoids = Context::getSolverProperty<MBool>("particleIncludeEllipsoids", solverId(), AT_, &ellipsoids);

  if(ellipsoids > 0) {
    m_includeEllipsoids = true;

    /*! \page propertyPageLPT
         \section particleEllipsoidCCD
         <code>MFloat LPT::m_particleEllipsoidCCD</code>\n
         default = <code>1</code> \n \n
         Specifies if ContinousCollisionDetction (Choi) or DistanceOfClosedApproach is used \n
         Keywords: <i>PARTICLE</i>
      */
    m_ellipsoidCCD = Context::getSolverProperty<MInt>("particleEllipsoidCCD", solverId(), AT_, &m_ellipsoidCCD);
  }

  /*! \page propertyPageLPT
      \section particleSubDomainFactor
      default = <code>1</code> \n \n
      Number of subdomains per coordinate direction for collision detection. Note that each solver is divided
      separately into the specified number of subdomains per direction, and the
      subdomains are rectangular divisions of the solver bounding box. \n
      Keywords: <i>PARTICLE</i>
   */
  MInt particleSubDomainFactor[MAX_SPACE_DIMENSIONS];
  for(MInt i = 0; i < nDim; i++) {
    particleSubDomainFactor[i] = 1;
    particleSubDomainFactor[i] =
        Context::getSolverProperty<MInt>("particleSubDomainFactor", solverId(), AT_, &particleSubDomainFactor[i], i);
  }

  // Analyze grid for subdomain division
  //------------------------------------
  // Note: It is to be noted that the current method for the division into
  // subdomains is intended for convex domains. For non-convex domains
  // a lot of memory is wasted by declaring subdomains outside of the
  // fluid...

  // find minimum and maximum coordinate values
  MFloat maxCoord[MAX_SPACE_DIMENSIONS], halfLength;
  for(MInt i = 0; i < nDim; i++) {
    halfLength = m_lpt->c_cellLengthAtLevel(m_lpt->c_level(0) + 1);
    m_subDomainCoordOffset[i] = m_lpt->c_coordinate(0, i) - halfLength;
    maxCoord[i] = m_lpt->c_coordinate(0, i) + halfLength;
  }

  for(MInt j = 1; j < m_lpt->a_noCells(); j++) {
    for(MInt i = 0; i < nDim; i++) {
      halfLength = m_lpt->c_cellLengthAtLevel(m_lpt->c_level(j) + 1);
      if(m_subDomainCoordOffset[i] > m_lpt->c_coordinate(j, i) - halfLength) {
        m_subDomainCoordOffset[i] = m_lpt->c_coordinate(j, i) - halfLength;
      }
      if(maxCoord[i] < m_lpt->c_coordinate(j, i) + halfLength) {
        maxCoord[i] = m_lpt->c_coordinate(j, i) + halfLength;
      }
    }
  }

  // find the appropriate number of subdomains in each direction
  for(MInt i = 0; i < nDim; i++) {
    m_subDomainSize[i] = (maxCoord[i] - m_subDomainCoordOffset[i]) / ((MFloat)particleSubDomainFactor[i]);
    m_noOfSubDomains[i] = particleSubDomainFactor[i];
    m_totalSubDomains *= m_noOfSubDomains[i];
  }
  m_log << m_totalSubDomains << " subdomains for collision detection in this solver." << endl;

  subDomainContent = new subDomainCollector<nDim>[m_totalSubDomains];
  if(m_includeEllipsoids) {
    subDomainContentEllipsoid = new subDomainCollectorEllipsoid<nDim>[m_totalSubDomains];
  }
}




/// Disables the following warnings with gcc: -Warray-bounds
#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

/** void LPT::detectPartColl()
 *
 *  \brief Particle collision detection & execution
 *
 *  version: March-09
 *  @author Rudie Kunnen, Sven Berger
 */
template <MInt nDim>
void ParticleCollision<nDim>::detectPartColl(std::vector<LPTSpherical<nDim>>& partList,
                                             std::vector<LPTEllipsoidal<nDim>>& partListEllipsoid,
                                             const MFloat timeStep,
                                             const MFloat time) {
  TRACE();

  collStruct<nDim> thisColl;

  m_timeStep = timeStep;

  MFloat currCollTime = 0.0, outtime;
  MInt subIndex[3] = {numeric_limits<MInt>::max(), numeric_limits<MInt>::max(), numeric_limits<MInt>::max()};

  MInt index = 1, index2 = 0;
  MInt factor[MAX_SPACE_DIMENSIONS] = {numeric_limits<MInt>::max(), numeric_limits<MInt>::max(),
                                       numeric_limits<MInt>::max()};
  MFloat thisCollTime;

  // factor[i] is used in the positioning in the storage array:
  // index = factor[0] * (x_steps) + factor[1] * (y_steps)
  //       + factor[2] * (z_steps)
  factor[0] = 1;
  factor[1] = m_noOfSubDomains[0];
  IF_CONSTEXPR(nDim == 3) { factor[2] = m_noOfSubDomains[0] * m_noOfSubDomains[1]; }

  // *** construct the lists of particles in the subdomains ***
  partListIterator<nDim> i1;
  for(i1 = partList.begin(); i1 != partList.end(); i1++) {
    // skip particles that have moved out of the current solver
    if(!(*i1).toBeDeleted()) {
      index = 0;
      index2 = 0;
      MBool error = false;
      for(MInt i = 0; i < nDim; i++) {
        index2 = MInt(floor(((*i1).m_position[i] - m_subDomainCoordOffset[i]) / m_subDomainSize[i]));
        if((index2 < 0) || (index2 >= m_noOfSubDomains[i])) // IMPORTANT added = to > !!!
        {
          error = true;
          break;
        } else {
          index += factor[i] * index2;
        }
      }

      if(error) {
        if(!(*i1).isInvalid()) {
          stringstream errorMessage;
          errorMessage << " ERROR: timestep " << globalTimeStep << " proc " << domainId() << " - particle "
                       << (*i1).m_partId << " is lost!" // << endl;
                       << " coord " << (*i1).m_position[0] << " " << (*i1).m_position[1] << " " << (*i1).m_position[2]
                       << " oldcoord " << (*i1).m_oldPos[0] << " " << (*i1).m_oldPos[1] << " " << (*i1).m_oldPos[2]
                       << " cellId " << (*i1).m_cellId << " coords " << m_lpt->c_coordinate((*i1).m_cellId, 0) << " "
                       << m_lpt->c_coordinate((*i1).m_cellId, 1) << " " << m_lpt->c_coordinate((*i1).m_cellId, 2)
                       << " halfLength " << m_lpt->c_cellLengthAtLevel(m_lpt->c_level((*i1).m_cellId) + 1) << endl;
          m_log << errorMessage.str();
        }
        if(!(*i1).wasSend() && !(*i1).toBeDeleted()) {
          (*i1).toBeDeleted() = true; // delete the lost particle
          m_log << "particle was lost " << endl;
        }
      } else {
        subDomainContent[index].subDomain.push_back(i1);
      }
    }
  }

  if(m_includeEllipsoids) {
    // *** construct the lists of particles in the subdomains ***
    ellipsListIterator<nDim> i2 = partListEllipsoid.begin();
    while(i2 != partListEllipsoid.end()) {
      // skip particles that have moved out of the current solver
      if(!(*i2).toBeDeleted()) {
        index = 0;
        MBool error = false;
        for(MInt i = 0; i < nDim; i++) {
          index2 = MInt(floor(((*i2).m_position[i] - m_subDomainCoordOffset[i]) / m_subDomainSize[i]));
          if((index2 < 0) || (index2 >= m_noOfSubDomains[i])) {
            error = true;
            break;
          } else {
            index += factor[i] * index2;
          }
        }

        if(error) {
          if(!(*i2).isInvalid()) {
            stringstream errorMessage;
            errorMessage << " ERROR: timestep " << globalTimeStep << " proc " << domainId() << " - particle "
                         << (*i2).m_partId << " is lost!" // << endl;
                         << " coord " << (*i2).m_position[0] << " " << (*i2).m_position[1] << " " << (*i2).m_position[2]
                         << " cellId " << (*i1).m_cellId << " coords " << m_lpt->c_coordinate((*i2).m_cellId, 0) << " "
                         << m_lpt->c_coordinate((*i2).m_cellId, 1) << " " << m_lpt->c_coordinate((*i2).m_cellId, 2)
                         << " halfLength " << m_lpt->c_cellLengthAtLevel(m_lpt->c_level((*i2).m_cellId) + 1) << endl;
            m_log << errorMessage.str();
          }
          if(!(*i2).toBeDeleted() && !(*i2).wasSend()) {
            (*i2).toBeDeleted() = true; // delete the lost particle
            (*i2).toBeRespawn() = true;
            (*i2).isInvalid() = true;
            m_log << "particle was lost " << endl;
          }
        } else {
          subDomainContentEllipsoid[index].subDomain.push_back(i2);
        }
      }
      i2++;
    }
  }


  // Construct the list of collisions for the current time step.
  // Take each subdomain and its direct "forward" neighbors,
  // then check for collisions between particle pairs in these restricted
  // regions.
  // array indices:
  // i - loops over all subdomains
  // j - loops over particles in subdomains to be considered as particle 1
  //     in collision detection
  // k - loops over neighboring subdomains: (k % 2) steps in x-dir,
  //     (k / 2) % 2 steps in y-dir, (k % 4) steps in z-dir
  // l - loops over particles in subdomains to be considered as particle 2
  //     in collision detection
  for(MInt i = 0; i < m_totalSubDomains; i++) {
    subIndex[0] = i % m_noOfSubDomains[0];
    subIndex[1] = (i / m_noOfSubDomains[0]) % m_noOfSubDomains[1];
    subIndex[2] = i / (m_noOfSubDomains[0] * m_noOfSubDomains[1]);
    auto& currDomain = subDomainContent[i];

    for(MInt j = 0; j < (MInt)currDomain.subDomain.size(); j++) {
      // search for collisions in the current subdomain i:
      // (only considering particle pairs that have not been considered yet)
      for(MInt l = j + 1; l < (MInt)currDomain.subDomain.size(); l++) {
        thisCollTime =
            collisionCheck(currDomain.subDomain.at((MUlong)j), currDomain.subDomain.at((MUlong)l), currCollTime);
        if(thisCollTime >= 0.0) // when no collision within the current
        {                       // time step a negative time is returned
          // add new collision event to the list
          thisColl.particle0 = currDomain.subDomain.at((MUlong)j);
          thisColl.particle1 = currDomain.subDomain.at((MUlong)l);
          thisColl.collTime = thisCollTime;
          thisColl.bndryId = -1;
          collList.push_back(thisColl);
        }
      }

      // search for collisions in neighboring subdomains
      // only the "forward" neighboring subdomains and the three "backward
      // diagonal" neighbors are included, the others have already been
      // investigated
      for(MInt k = 1; k < IPOW2(nDim); k++) { // these are the "forward" neighbors (7 total for 3D, 3 for 2D)

        // edge detection: no particle subdomains outside of domain
        if(((subIndex[0] + (k % 2)) < m_noOfSubDomains[0]) && ((subIndex[1] + ((k / 2) % 2)) < m_noOfSubDomains[1])) {
          IF_CONSTEXPR(nDim == 3) {
            if(!((subIndex[2] + (k / 4)) < m_noOfSubDomains[2])) {
              continue;
            }
          }
            index =
                i + (k % 2) + ((k / 2) % 2) * m_noOfSubDomains[0] + (k / 4) * m_noOfSubDomains[0] * m_noOfSubDomains[1];
            for(MInt l = 0; l < (MInt)subDomainContent[index].subDomain.size(); l++) {
              thisCollTime = collisionCheck(subDomainContent[i].subDomain.at((MUlong)j),
                                            subDomainContent[index].subDomain.at((MUlong)l),
                                            currCollTime);
              if(thisCollTime >= 0.0) // when no collision within the
              {                       // current time step a negative time is returned
                // add new collision event to the list
                thisColl.particle0 = subDomainContent[i].subDomain.at((MUlong)j);
                thisColl.particle1 = subDomainContent[index].subDomain.at((MUlong)l);
                thisColl.collTime = thisCollTime;
                thisColl.bndryId = -1; // RUDIE new since particle vector -> list
                collList.push_back(thisColl);
              }
            }
        }
      }
      // now the "backward diagonal" neighbors (3 for 3D, 1 for 2D)
      if((subIndex[0] > 0) && (subIndex[1] < (m_noOfSubDomains[1] - 1))) {
        index = i + m_noOfSubDomains[0] - 1;
        for(MInt l = 0; l < (MInt)subDomainContent[index].subDomain.size(); l++) {
          thisCollTime = collisionCheck(subDomainContent[i].subDomain.at((MUlong)j),
                                        subDomainContent[index].subDomain.at((MUlong)l),
                                        currCollTime);
          if(thisCollTime >= 0.0) // when no collision within the current
          {                       // time step a negative time is returned
            // add new collision event to the list
            thisColl.particle0 = subDomainContent[i].subDomain.at((MUlong)j);
            thisColl.particle1 = subDomainContent[index].subDomain.at((MUlong)l);
            thisColl.collTime = thisCollTime;
            thisColl.bndryId = -1;
            collList.push_back(thisColl);
          }
        }
      }
      IF_CONSTEXPR(nDim == 3) {
        if((subIndex[0] > 0) && (subIndex[2] < (m_noOfSubDomains[2] - 1))) {
          index = i + m_noOfSubDomains[0] * m_noOfSubDomains[1] - 1;
          for(MInt l = 0; l < (MInt)subDomainContent[index].subDomain.size(); l++) {
            thisCollTime = collisionCheck(subDomainContent[i].subDomain.at((MUlong)j),
                                          subDomainContent[index].subDomain.at((MUlong)l),
                                          currCollTime);
            if(thisCollTime >= 0.0) // when no collision within the current
            {                       // time step a negative time is returned
              // add new collision event to the list
              thisColl.particle0 = subDomainContent[i].subDomain.at((MUlong)j);
              thisColl.particle1 = subDomainContent[index].subDomain.at((MUlong)l);
              thisColl.collTime = thisCollTime;
              thisColl.bndryId = -1;
              collList.push_back(thisColl);
            }
          }
        }
      }
      IF_CONSTEXPR(nDim == 3) {
        if((subIndex[0] > 0) && (subIndex[1] < (m_noOfSubDomains[1] - 1))
           && (subIndex[2] < (m_noOfSubDomains[2] - 1))) {
          index = i + m_noOfSubDomains[0] * (1 + m_noOfSubDomains[1]) - 1;
          for(MInt l = 0; l < (MInt)subDomainContent[index].subDomain.size(); l++) {
            thisCollTime = collisionCheck(subDomainContent[i].subDomain.at((MUlong)j),
                                          subDomainContent[index].subDomain.at((MUlong)l),
                                          currCollTime);
            if(thisCollTime >= 0.0) // when no collision within the current
            {                       // time step a negative time is returned
              // add new collision event to the list
              thisColl.particle0 = subDomainContent[i].subDomain.at((MUlong)j);
              thisColl.particle1 = subDomainContent[index].subDomain.at((MUlong)l);
              thisColl.collTime = thisCollTime;
              thisColl.bndryId = -1; // RUDIE new since particle vector -> list
              collList.push_back(thisColl);
            }
          }
        }
      }

      if(m_includeEllipsoids) { // search for collisions of spheres with ellipsoids...
        // note that these collisions will only be detected & recorded;
        // no attempt at elastic collisions etc. is made!
        // NOTE: This is placed here since it is untested!
        for(MInt k = subIndex[0] - 1; k < subIndex[0] + 2; k++) {
          if((k < 0) || (k >= m_noOfSubDomains[0])) {
            continue;
          }
          for(MInt l = subIndex[1] - 1; l < subIndex[1] + 2; l++) {
            if((l < 0) || (l >= m_noOfSubDomains[1])) {
              continue;
            }
            for(MInt m = subIndex[2] - 1; m < subIndex[2] + 2; m++) {
              if((m < 0) || (m >= m_noOfSubDomains[2])) {
                continue;
              }
              index = k + m_noOfSubDomains[0] * l + m * m_noOfSubDomains[0] * m_noOfSubDomains[1];
              for(MInt n = 0; n < (MInt)subDomainContentEllipsoid[index].subDomain.size(); n++) {
                collisionCheckSphereEllipsoid(subDomainContent[i].subDomain.at((MUlong)j),
                                              subDomainContentEllipsoid[index].subDomain.at((MUlong)n));
              }
            }
          }
        }
      }
    }
  }

  if(m_includeEllipsoids) {
    index = 0;
    index2 = 0;
    // Detect overlap of ellipsoids.
    // Take each subdomain and its direct "forward" neighbors,
    // then check for collisions between particle pairs in these restricted
    // regions.
    // array indices:
    // i - loops over all subdomains
    // j - loops over particles in subdomains to be considered as particle 1
    //     in collision detection
    // k - loops over neighboring subdomains: (k % 2) steps in x-dir,
    //     (k / 2) % 2 steps in y-dir, (k % 4) steps in z-dir
    // l - loops over particles in subdomains to be considered as particle 2
    //     in collision detection
    for(MInt i = 0; i < m_totalSubDomains; i++) {
      subIndex[0] = i % m_noOfSubDomains[0];
      subIndex[1] = (i / m_noOfSubDomains[0]) % m_noOfSubDomains[1];
      subIndex[2] = i / (m_noOfSubDomains[0] * m_noOfSubDomains[1]);

      for(MInt j = 0; j < (MInt)subDomainContentEllipsoid[i].subDomain.size(); j++) {
        // search for collisions in the current subdomain i:
        // (only considering particle pairs that have not been considered yet)
        for(MInt l = j + 1; l < (MInt)subDomainContentEllipsoid[i].subDomain.size(); l++) {
          collisionCheckEllipsoidEllipsoid(subDomainContentEllipsoid[i].subDomain.at((MUlong)j),
                                           subDomainContentEllipsoid[i].subDomain.at((MUlong)l));
        }


        // search for collisions in neighboring subdomains
        // only the "forward" neighboring subdomains and the three "backward
        // diagonal" neighbors are included, the others have already been
        // investigated
        for(MInt k = 1; k < 8; k++) { // these are the "forward" neighbors (7 total for 3D)

          // edge detection: no particle subdomains outside of domain
          if(((subIndex[0] + (k % 2)) < m_noOfSubDomains[0]) && ((subIndex[1] + ((k / 2) % 2)) < m_noOfSubDomains[1])) {
            IF_CONSTEXPR(nDim == 3) {
              if(!((subIndex[2] + (k / 4)) < m_noOfSubDomains[2])) {
                continue;
              }
            }
            index =
                i + (k % 2) + ((k / 2) % 2) * m_noOfSubDomains[0] + (k / 4) * m_noOfSubDomains[0] * m_noOfSubDomains[1];
            for(MInt l = 0; l < (MInt)subDomainContentEllipsoid[index].subDomain.size(); l++) {
              collisionCheckEllipsoidEllipsoid(subDomainContentEllipsoid[i].subDomain.at((MUlong)j),
                                               subDomainContentEllipsoid[index].subDomain.at((MUlong)l));
            }
          }
        }
        // now the "backward diagonal" neighbors (3 for 3D)
        if((subIndex[0] > 0) && (subIndex[1] < (m_noOfSubDomains[1] - 1))) {
          index = i + m_noOfSubDomains[0] - 1;
          for(MInt l = 0; l < (MInt)subDomainContentEllipsoid[index].subDomain.size(); l++) {
            collisionCheckEllipsoidEllipsoid(subDomainContentEllipsoid[i].subDomain.at((MUlong)j),
                                             subDomainContentEllipsoid[index].subDomain.at((MUlong)l));
          }
        }
        if((subIndex[0] > 0) && (subIndex[2] < (m_noOfSubDomains[2] - 1))) {
          index = i + m_noOfSubDomains[0] * m_noOfSubDomains[1] - 1;
          for(MInt l = 0; l < (MInt)subDomainContentEllipsoid[index].subDomain.size(); l++) {
            collisionCheckEllipsoidEllipsoid(subDomainContentEllipsoid[i].subDomain.at((MUlong)j),
                                             subDomainContentEllipsoid[index].subDomain.at((MUlong)l));
          }
        }
        if((subIndex[0] > 0) && (subIndex[1] < (m_noOfSubDomains[1] - 1))
           && (subIndex[2] < (m_noOfSubDomains[2] - 1))) {
          index = i + m_noOfSubDomains[0] * (1 + m_noOfSubDomains[1]) - 1;
          for(MInt l = 0; l < (MInt)subDomainContentEllipsoid[index].subDomain.size(); l++) {
            collisionCheckEllipsoidEllipsoid(subDomainContentEllipsoid[i].subDomain.at((MUlong)j),
                                             subDomainContentEllipsoid[index].subDomain.at((MUlong)l));
          }
        }
      }
    }
  }

  MInt bndryId;

  // loop until all collisions are accounted for
  while(!collList.empty()) {
    // sort the collision list and take the first element
    collList.sort();
    thisColl = collList.front();

    partListIterator<nDim> collPair[2];
    bndryId = thisColl.bndryId;
    outtime = thisColl.collTime;
    collPair[0] = thisColl.particle0;
    collPair[1] = thisColl.particle1;

    if(bndryId > -1) // wall collision
    {
      if(!(*collPair[0]).toBeDeleted()) {
        switch(m_lpt->a_bndryCellId(bndryId)) {
          case 1012:
          case 1002: {
            // OUTLET BOUNDARY
            if(!(*collPair[0]).wasSend() && !(*collPair[0]).toBeDeleted()) {
              m_log << "particle left domain through outlet boundary" << endl;
              (*collPair[0]).toBeDeleted() = true;
              (*collPair[0]).toBeRespawn() = true;
              (*collPair[0]).isInvalid() = true;
            }
            break;
          }

          case 1021:
          case 1001: {
            // INLET BOUNDARY
            // check if particle really left the computational domain
            MFloat l = checkBndryCross(bndryId, collPair[0], 0.0);
            if(l >= 0.0 && l <= 1.1) {
              // particle has left the area through inlet boundary
              if(!(*collPair[0]).wasSend() && !(*collPair[0]).toBeDeleted()) {
                m_log << "particle left domain through inlet boundary" << endl;
                (*collPair[0]).toBeDeleted() = true;
                (*collPair[0]).toBeRespawn() = true;
                (*collPair[0]).isInvalid() = true;
              }
            }
            break;
          }

          case 3002:
          case 3003: // (mgc)
          case 3803: // driven cavity lid
          {
            // WALL COLLISION
            // particle has left the area through wall boundary
            // => wall collision is executed
            MFloat normVel = 0.0;
            MFloat bndryBasis[MAX_SPACE_DIMENSIONS][MAX_SPACE_DIMENSIONS];
            MFloat normN[MAX_SPACE_DIMENSIONS - 1];
            MFloat reboundVel[MAX_SPACE_DIMENSIONS]; //, pointOfImpact[MAX_SPACE_DIMENSIONS];
            MFloat oldCoordinates[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                                           numeric_limits<MFloat>::max()};
            std::array<std::array<MFloat, MAX_SPACE_DIMENSIONS>, MAX_SPACE_DIMENSIONS> A;
            std::array<MFloat, MAX_SPACE_DIMENSIONS> b;
            fill_n(A[0].data(), MAX_SPACE_DIMENSIONS * MAX_SPACE_DIMENSIONS, 0.0);
            fill_n(b.data(), MAX_SPACE_DIMENSIONS, 0.0);
            //            MFloat Ab[MAX_SPACE_DIMENSIONS * (MAX_SPACE_DIMENSIONS + 1)];
            MFloat lambda = outtime / timeStep;

            for(MInt i = 0; i < nDim; i++) {
              oldCoordinates[i] =
                  (*collPair[0]).m_oldPos[i] + lambda * ((*collPair[0]).m_position[i] - (*collPair[0]).m_oldPos[i]);
              normVel += POW2(m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i]
                              * ((*collPair[0]).m_position[i] - (*collPair[0]).m_oldPos[i]) / timeStep);
              // RUDIE the minus sign in front of the 0.5 is there assuming
              // that the surface normal always points into the fluid!
              //               pointOfImpact[i] = -0.5 * (*collPair[0]).m_diameter
              //                       * m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normalVector[i]
              //                       + lambda  * (*collPair[0]).m_position[i]
              //                + (1.0 - lambda) * (*collPair[0]).m_oldPos[i];
            }

            // wall-normal velocity at collision
            // normVel = sqrt(normVel);

            // compute the basis w.r.t. the boundary surface
            for(MInt i = 0; i < nDim - 1; i++) {
              normN[i] = 0.0;
            }

            for(MInt i = 0; i < nDim - 1; i++) {
              for(MInt j = 0; j < nDim; j++) {
                normN[i] += POW2(m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[i + 1][j]
                                 - m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[0][j]);
              }
            }

            for(MInt i = 0; i < nDim - 1; i++) {
              normN[i] = sqrt(normN[i]);
            }

            for(MInt i = 0; i < nDim; i++) {
              bndryBasis[i][0] = m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i];
              for(MInt j = 1; j < nDim; j++) {
                bndryBasis[i][j] = (m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[j][i]
                                    - m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[0][i])
                                   / normN[j - 1];
              }
            }

            // perform Gaussian elimination to transform the velocity into the
            // boundary coordinate system
            for(MInt i = 0; i < nDim; i++) {
              A[0][i] = m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i];
              for(MInt j = 1; j < nDim; j++) {
                A[j][i] = (m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[j][i]
                           - m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[0][i])
                          / normN[j - 1];
              }
              b[i] = (*collPair[0]).m_oldVel[i];
            }

            maia::math::solveQR<3>(A, b);

            // Last column of Ab now contains the velocity components w.r.t. the
            // boundary coordinate system. Compute rebound velocity, set old
            // velocity to zero, to compute transformed velocity.
            (*collPair[0]).m_oldVel[0] = 0.0;
            reboundVel[0] = -b[0];

            for(MInt i = 1; i < nDim; i++) {
              (*collPair[0]).m_oldVel[i] = 0.0;
              (*collPair[0]).m_velocity[i] = 0.0;
              reboundVel[i] = b[i];
            }

            // Transform the rebound velocity into the x,y,z coordinate system
            // and save as old velocity. Compute new coordinates.
            for(MInt i = 0; i < nDim; i++) {
              for(MInt j = 0; j < nDim; j++) {
                (*collPair[0]).m_oldVel[i] += bndryBasis[i][j] * reboundVel[j];
                (*collPair[0]).m_velocity[i] = (*collPair[0]).m_oldVel[i];
              }
              (*collPair[0]).m_position[i] = oldCoordinates[i] + (1.0 - lambda) * timeStep * (*collPair[0]).m_oldVel[i];
            }
            // update cellId
            (*collPair[0]).checkCellChange(&(*collPair[0]).m_oldPos[0]);

            // compute old coordinates
            for(MInt i = 0; i < nDim; i++) {
              (*collPair[0]).m_oldPos[i] = oldCoordinates[i] + (lambda - 1.0) * timeStep * (*collPair[0]).m_oldVel[i];
            }

            (*collPair[0]).hasCollided() = true;
            (*collPair[0]).firstStep() = true;

            // MFloat collTime = solverPtr->m_time - m_particleTimeStep + outtime;

            // write collision event if particle is real on this proc.
            // RUDIE writing of wall collisions switched off for now...
            //	    if (writeThisColl)
            //            writeCollEvent(collTime, (*collPair[0]).m_partId, -1, (*collPair[0]).m_diameter,
            //		  	   0.0, (*collPair[0]).m_densityRatio, 0.0, normVel, -1.0,
            //			  &pointOfImpact[0]);

            break;
          }
          default: {
            // boundary condition id not implemented for particles!
            cerr << "WARNING: This boundary condition not implemented for particles!" << endl;
            cerr << "Treated as a solid wall" << endl;
          }
        }
      }

      // remove collisions still in the list that involve collPair[0]
      // collPair[1] is not relevant since the current collision is with a wall
      collList.remove_if(compareParticleIds<nDim>(collPair[0]));

      // If particle is still active, check for new collisions of this particle
      // with others. If yes, add new collision to the list and re-sort it
      if(!(*collPair[0]).toBeDeleted()) {
        MInt subCollIndex[3] = {numeric_limits<MInt>::max(), numeric_limits<MInt>::max(),
                                numeric_limits<MInt>::max()}; // gcc 4.8.2 maybe uninitialized
        MInt _index2;
        for(MInt i = 0; i < nDim; i++) {
          subCollIndex[i] = (MInt)(((*collPair[0]).m_position[i] - m_subDomainCoordOffset[i]) / m_subDomainSize[i]);
        }

        IF_CONSTEXPR(nDim == 2) {
          // 2D -> lots of for loops in the following:
          // i - loops over -1, 0, +1 subdomain in x direction
          // j - loops over -1, 0, +1 subdomain in y direction
          // l - loops over particles in subdomain
          for(MInt i = (subCollIndex[0] - 1); i < (subCollIndex[0] + 2); i++) {
            for(MInt j = (subCollIndex[1] - 1); j < (subCollIndex[1] + 2); j++) {
              if((i >= 0) && (i < m_noOfSubDomains[0]) && (j >= 0) && (j < m_noOfSubDomains[1])) {
                _index2 = i + j * m_noOfSubDomains[0];
                for(MInt l = 0; l < (MInt)subDomainContent[_index2].subDomain.size(); l++) {
                  if((collPair[0] != subDomainContent[_index2].subDomain.at((MUlong)l))
                     && !((*subDomainContent[_index2].subDomain.at((MUlong)l)).toBeDeleted())) {
                    MFloat _thisCollTime =
                        collisionCheck(collPair[0], subDomainContent[_index2].subDomain.at((MUlong)l), currCollTime);
                    if(_thisCollTime >= 0.0) // when no collision within the
                    {                        // current time step a negative time is returned
                      // add new collision event to the list
                      thisColl.particle0 = collPair[0];
                      thisColl.particle1 = subDomainContent[_index2].subDomain.at((MUlong)l);
                      thisColl.collTime = _thisCollTime;
                      thisColl.bndryId = -1;
                      collList.push_back(thisColl);
                    }
                  }
                }
              }
            }
          }
        }
        else {
          // 3D -> lots of for loops in the following:
          // i - loops over -1, 0, +1 subdomain in x direction
          // j - loops over -1, 0, +1 subdomain in y direction
          // k - loops over -1, 0, +1 subdomain in z direction
          // l - loops over particles in subdomain
          for(MInt i = (subCollIndex[0] - 1); i < (subCollIndex[0] + 2); i++) {
            for(MInt j = (subCollIndex[1] - 1); j < (subCollIndex[1] + 2); j++) {
              for(MInt k = (subCollIndex[2] - 1); k < (subCollIndex[2] + 2); k++) {
                if((i >= 0) && (i < m_noOfSubDomains[0]) && (j >= 0) && (j < m_noOfSubDomains[1]) && (k >= 0)
                   && (k < m_noOfSubDomains[2])) {
                  _index2 = i + j * m_noOfSubDomains[0] + k * m_noOfSubDomains[0] * m_noOfSubDomains[1];
                  for(MInt l = 0; l < (MInt)subDomainContent[_index2].subDomain.size(); l++) {
                    if((collPair[0] != subDomainContent[_index2].subDomain.at((MUlong)l))
                       && !((*subDomainContent[_index2].subDomain.at((MUlong)l)).toBeDeleted())) {
                      MFloat __thisCollTime =
                          collisionCheck(collPair[0], subDomainContent[_index2].subDomain.at((MUlong)l), currCollTime);
                      if(__thisCollTime >= 0.0) // when no collision within the
                      {                         // current time step a negative time is returned
                        // add new collision event to the list
                        thisColl.particle0 = collPair[0];
                        thisColl.particle1 = subDomainContent[_index2].subDomain.at((MUlong)l);
                        thisColl.collTime = __thisCollTime;
                        thisColl.bndryId = -1;
                        collList.push_back(thisColl);
                      }
                    }
                  }
                }
              }
            }
          }
        }

        // check for new wall collisions of this particle
        geometryInteraction();
      }

      currCollTime = outtime;

    } else { // current collision is a particle-particle collision
      MFloat wdotr = 0.0, relVel = 0.0, meanVel = 0.0, deltap = 0.0;
      MFloat vi0[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                          numeric_limits<MFloat>::max()};
      MFloat vi1[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                          numeric_limits<MFloat>::max()};
      MFloat r[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                        numeric_limits<MFloat>::max()};
      MFloat x0old[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                            numeric_limits<MFloat>::max()};
      MFloat x1old[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                            numeric_limits<MFloat>::max()};
      MFloat partMass0, partMass1;
      MFloat pointOfImpact[MAX_SPACE_DIMENSIONS] = {numeric_limits<MFloat>::max(), numeric_limits<MFloat>::max(),
                                                    numeric_limits<MFloat>::max()};
      MFloat rLength = (*collPair[1]).m_diameter / ((*collPair[0]).m_diameter + (*collPair[1]).m_diameter);
      MInt qmax = 0;

      // do some MBool algebra to determine whether collision is written by this processor
      // in order to prevent multiple outputs of collisions on domain boundaries
      // 1) collisions between 2 ghost particles are not written by this proc.
      // 2) collisions between 2 real particles are written by this proc. by default
      // 3) a collision between a real and a ghost particle is written by this proc. only
      //    if the real particle has the lowest m_partId
      MBool writeThisColl = true;
      MBool haloPart0 = (*collPair[0]).wasSend();
      MBool haloPart1 = (*collPair[1]).wasSend();
      MBool larger = ((*collPair[0]).m_partId > (*collPair[1]).m_partId);
      if(haloPart0 & haloPart1) {
        writeThisColl = false;
      }
      if(haloPart0 ^ haloPart1) {
        if(haloPart1 & larger) {
          writeThisColl = false;
        }
        if(haloPart0 & !larger) {
          writeThisColl = false;
        }
      }

      // particle velocities before collision, and separation vector at instant
      // of collision:
      for(MInt k = 0; k < nDim; k++) {
        vi0[k] = ((*collPair[0]).m_position[k] - (*collPair[0]).m_oldPos[k]) / timeStep;
        vi1[k] = ((*collPair[1]).m_position[k] - (*collPair[1]).m_oldPos[k]) / timeStep;
        r[k] = ((*collPair[0]).m_oldPos[k] + outtime * vi0[k]) - ((*collPair[1]).m_oldPos[k] + outtime * vi1[k]);
        x0old[k] = (*collPair[0]).m_oldPos[k];
        x1old[k] = (*collPair[1]).m_oldPos[k];
        wdotr += r[k] * (vi0[k] - vi1[k]);
        relVel += POW2(vi0[k] - vi1[k]);
        meanVel += POW2(vi0[k] + vi1[k]);
        // subtract the mean flow velocity for the calculation of the meanVel
        pointOfImpact[k] = (*collPair[1]).m_oldPos[k] + outtime * vi1[k] + r[k] * rLength;
      }

      // relVel is relative velocity at instance of collision
      // meanVel is mean velocity at instance of collision
      // pointOfImpact is position where the particles touch
      relVel = sqrt(relVel);
      meanVel = 0.5 * sqrt(meanVel);

      partMass0 = F1B6 * PI * (*collPair[0]).m_densityRatio * POW3((*collPair[0]).m_diameter);
      partMass1 = F1B6 * PI * (*collPair[1]).m_densityRatio * POW3((*collPair[1]).m_diameter);

      if((m_collisionModel == 1) || (m_collisionModel == 2)) { // hard-sphere collision (no coagulation)
        // particle velocities and positions after collision:
        qmax = 2;                        // after collision, check for new collisions for both particles
        for(MInt k = 0; k < nDim; k++) { // deltap is change in momentum without m1*m2 in the numerator
          deltap = -2.0 * wdotr * r[k]
                   / ((partMass0 + partMass1) * 0.25 * POW2((*collPair[0]).m_diameter + (*collPair[1]).m_diameter));
          (*collPair[0]).m_velocity[k] = vi0[k] + deltap * partMass1;
          (*collPair[1]).m_velocity[k] = vi1[k] - deltap * partMass0;
          (*collPair[0]).m_position[k] =
              (*collPair[0]).m_oldPos[k] + vi0[k] * outtime + (*collPair[0]).m_velocity[k] * (timeStep - outtime);
          (*collPair[1]).m_position[k] =
              (*collPair[1]).m_oldPos[k] + vi1[k] * outtime + (*collPair[1]).m_velocity[k] * (timeStep - outtime);
        }

        // check cell change before updating the old particle coordinates
        (*collPair[0]).checkCellChange(x0old);
        (*collPair[1]).checkCellChange(x1old);

        for(MInt k = 0; k < nDim; k++) { // updating particle oldCoords -> necessary in case of multiple
          // collisions in current time step!
          (*collPair[0]).m_oldPos[k] = (*collPair[0]).m_position[k] - (*collPair[0]).m_velocity[k] * timeStep;
          (*collPair[1]).m_oldPos[k] = (*collPair[1]).m_position[k] - (*collPair[1]).m_velocity[k] * timeStep;
        }

        (*collPair[0]).hasCollided() = true;
        (*collPair[0]).firstStep() = true;
        (*collPair[1]).hasCollided() = true;
        (*collPair[1]).firstStep() = true;

        if(m_collisionModel == 1) {
          MFloat part2CFL = 0.0;
          for(MInt k = 0; k < nDim; k++) {
            part2CFL += POW2((*collPair[1]).m_velocity[k]);
          }
          part2CFL *= POW2(timeStep / (*collPair[1]).m_diameter);
          if(part2CFL > 0.04) {
            cout << "Caution: a particle CFL number of " << sqrt(part2CFL) << " has been found!\n";
          }
        }

        // remove future collisions involving these two particles from list
        collList.remove_if(compareParticleIds<nDim>(collPair[0]));
        collList.remove_if(compareParticleIds<nDim>(collPair[1]));
      }

      if((m_collisionModel == 3) || (m_collisionModel == 4)) { // instantaneous coagulation
        // particle 0 is the single particle after collision, with partId the
        // lesser of the two original partIds
        if((*collPair[0]).m_partId > (*collPair[1]).m_partId) { // swap partIds if necessary
          qmax = (*collPair[0]).m_partId;
          (*collPair[0]).m_partId = (*collPair[1]).m_partId;
          (*collPair[1]).m_partId = qmax;
        }
        qmax = 1; // after collision, only collPair[0] needs to be checked for
        // new collisions

        MFloat FSumMass = 1.0 / (partMass0 + partMass1);

        for(MInt k = 0; k < nDim; k++) {
          (*collPair[0]).m_velocity[k] =
              FSumMass * (partMass0 * (*collPair[0]).m_velocity[k] + partMass1 * (*collPair[1]).m_velocity[k]);
          (*collPair[0]).m_oldVel[k] = (*collPair[1]).m_velocity[k];
          (*collPair[0]).m_position[k] = FSumMass
                                             * (partMass0 * ((*collPair[0]).m_oldPos[k] + outtime * vi0[k])
                                                + partMass1 * ((*collPair[1]).m_oldPos[k] + outtime * vi1[k]))
                                         + (timeStep - outtime) * (*collPair[0]).m_velocity[k];
        }

        // check wheter particle changed cells before updating old coordinates
        (*collPair[0]).checkCellChange(x0old);

        for(MInt k = 0; k < nDim; k++) {
          (*collPair[0]).m_oldPos[k] = (*collPair[0]).m_position[k] - timeStep * (*collPair[0]).m_velocity[k];
        }

        MFloat diam0 = (*collPair[0]).m_diameter;
        MFloat diam1 = (*collPair[1]).m_diameter;
        MFloat diam03 = POW3(diam0);
        MFloat diam13 = POW3(diam1);
        (*collPair[0]).m_diameter = pow(diam03 + diam13, F1B3);
        (*collPair[0]).m_densityRatio =
            (diam03 * (*collPair[0]).m_densityRatio + diam13 * (*collPair[1]).m_densityRatio) / (diam03 + diam13);

        // change the m_firstStep flag to true: particle acceleration is not
        // known anymore and thus simpler time integration required
        (*collPair[0]).hasCollided() = true;
        (*collPair[0]).firstStep() = true;

        // remove future collisions involving these two particles from list
        collList.remove_if(compareParticleIds<nDim>(collPair[0]));
        collList.remove_if(compareParticleIds<nDim>(collPair[1]));

        // particle collPair[1] of coagulation event is effectively removed
        (*collPair[1]).toBeDeleted() = true;
        (*collPair[0]).isInvalid() = true;
      }

      if((m_collisionModel == 5) || (m_collisionModel == 6)) { // just a registration of the collision event;
        // particles move through each other without interaction

        if(m_collisionModel == 5) {
          MFloat part2CFL = 0.0;
          for(MInt k = 0; k < nDim; k++) {
            part2CFL += POW2((*collPair[1]).m_velocity[k]);
          }
          part2CFL *= POW2(timeStep / (*collPair[1]).m_diameter);
          if(part2CFL > 0.04) {
            cout << "Caution: a particle CFL number of " << sqrt(part2CFL) << " has been found!\n";
          }
        }
        // remove the current collision from the list
        collList.pop_front();
      }

      // particle CFL is only relevant for retroactive particle collision
      // detection!
      if((m_collisionModel == 1) || (m_collisionModel == 3) || (m_collisionModel == 5)) {
        MFloat part0CFL = 0.0;
        for(MInt k = 0; k < nDim; k++) {
          part0CFL += POW2((*collPair[0]).m_velocity[k]);
        }
        part0CFL *= POW2(timeStep / (*collPair[0]).m_diameter);
        if(part0CFL > 0.04) {
          cout << "Caution: a particle CFL number of " << sqrt(part0CFL) << " has been found!\n";
        }
      }

      // update partial time in this time step, to be able to exclude unphysical
      // earlier collisions that may be found after updates of the "old"
      // particle velocities
      currCollTime = outtime;

      if((m_collisionModel != 5) && (m_collisionModel != 6)) {
        // check for new collisions of these particles with others
        // if yes, add new collision to the list and re-sort it
        IF_CONSTEXPR(nDim == 2) {
          // 2D -> lots of for loops in the following:
          // q - loops collPair[q], 0 and 1
          // i - loops over -1, 0, +1 subdomain in x direction
          // j - loops over -1, 0, +1 subdomain in y direction
          // l - loops over particles in subdomain
          for(MInt q = 0; q < qmax; q++) {
            if((*collPair[q]).toBeDeleted()) {
              continue;
            }
            for(MInt i = 0; i < nDim; i++) {
              subIndex[i] = (MInt)(((*collPair[q]).m_position[i] - m_subDomainCoordOffset[i]) / m_subDomainSize[i]);
            }

            for(MInt i = (subIndex[0] - 1); i < (subIndex[0] + 2); i++) {
              for(MInt j = (subIndex[1] - 1); j < (subIndex[1] + 2); j++) {
                if((i >= 0) && (i < m_noOfSubDomains[0]) && (j >= 0) && (j < m_noOfSubDomains[1])) {
                  MInt _index2 = i + j * m_noOfSubDomains[0];
                  for(MInt l = 0; l < (MInt)subDomainContent[_index2].subDomain.size(); l++) {
                    if((collPair[q] != subDomainContent[_index2].subDomain.at((MUlong)l))
                       && !((*subDomainContent[_index2].subDomain.at((MUlong)l)).toBeDeleted())) {
                      MFloat _thisCollTime =
                          collisionCheck(collPair[q], subDomainContent[_index2].subDomain.at((MUlong)l), currCollTime);
                      if(_thisCollTime >= 0.0) // when no collision within the
                      {                        // current time step a negative time is returned
                        // add new collision event to the list
                        thisColl.particle0 = collPair[q];
                        thisColl.particle1 = subDomainContent[_index2].subDomain.at((MUlong)l);
                        thisColl.collTime = _thisCollTime;
                        thisColl.bndryId = -1;
                        collList.push_back(thisColl);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        else {
          // 3D -> lots of for loops in the following:
          // q - loops collPair[q], 0 and 1
          // i - loops over -1, 0, +1 subdomain in x direction
          // j - loops over -1, 0, +1 subdomain in y direction
          // k - loops over -1, 0, +1 subdomain in z direction
          // l - loops over particles in subdomain
          for(MInt q = 0; q < qmax; q++) {
            if((*collPair[q]).toBeDeleted()) {
              continue;
            }
            for(MInt i = (subIndex[0] - 1); i < (subIndex[0] + 2); i++) {
              for(MInt j = (subIndex[1] - 1); j < (subIndex[1] + 2); j++) {
                for(MInt k = (subIndex[2] - 1); k < (subIndex[2] + 2); k++) {
                  if((i >= 0) && (i < m_noOfSubDomains[0]) && (j >= 0) && (j < m_noOfSubDomains[1]) && (k >= 0)
                     && (k < m_noOfSubDomains[2])) {
                    MInt __index2 = i + j * m_noOfSubDomains[0] + k * m_noOfSubDomains[0] * m_noOfSubDomains[1];
                    for(MInt l = 0; l < (MInt)subDomainContent[__index2].subDomain.size(); l++) {
                      if((collPair[q] != subDomainContent[__index2].subDomain.at((MUlong)l))
                         && !((*subDomainContent[__index2].subDomain.at((MUlong)l)).toBeDeleted())) {
                        MFloat __thisCollTime = collisionCheck(
                            collPair[q], subDomainContent[__index2].subDomain.at((MUlong)l), currCollTime);
                        if(__thisCollTime >= 0.0) // when no collision within the
                        {                         // current time step a negative time is returned
                          // add new collision event to the list
                          thisColl.particle0 = collPair[q];
                          thisColl.particle1 = subDomainContent[__index2].subDomain.at((MUlong)l);
                          thisColl.collTime = __thisCollTime;
                          thisColl.bndryId = -1;
                          collList.push_back(thisColl);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }

        if(m_lpt->m_wallCollisions) {
          // check for new wall collisions of these particles
          for(MInt q = 0; q < qmax; q++) {
            if((*collPair[q]).toBeDeleted()) {
              continue;
            }
            geometryInteraction();
          }
        }
      }

      // write collision event
      if(writeThisColl) {
        // only record collisions after x = m_particleCollOffset so that the particles have had time to
        // adapt to the flow
        if(pointOfImpact[0] > m_offset) {
          writeCollEvent(time, (*collPair[0]).m_partId, (*collPair[1]).m_partId, (*collPair[0]).m_diameter,
                         (*collPair[1]).m_diameter, (*collPair[0]).m_densityRatio, (*collPair[1]).m_densityRatio,
                         relVel, meanVel, pointOfImpact);
        }
      }
    }
  }

  for(MInt i = 0; i < m_totalSubDomains; i++) {
    subDomainContent[i].subDomain.clear();
  }

  if(m_includeEllipsoids) {
    for(MInt i = 0; i < m_totalSubDomains; i++) {
      subDomainContentEllipsoid[i].subDomain.clear();
    }
  }
}

#ifdef MAIA_GCC_COMPILER
#pragma GCC diagnostic pop
#endif

template <MInt nDim>
void ParticleCollision<nDim>::geometryInteraction() {
  TRACE();


  //  collStruct<nDim> thisColl;
  //
  //  // in case of a boundary cell, always check for wall collisions
  //  if(m_fvSolver->a_boundaryId((*partVectorElement).m_cellId) > -1) {
  //
  //    MFloat lambda = checkBndryCross(
  //        m_fvSolver->a_boundaryId((*partVectorElement).m_cellId), partVectorElement,
  //        (*partVectorElement).m_diameter);
  //
  //    if(lambda >= 0.0 && lambda <= 1.0) {
  //      thisColl.particle0 = partVectorElement;
  //      thisColl.particle1 = partVectorElement;
  //      thisColl.bndryId = m_fvSolver->a_boundaryId((*partVectorElement).m_cellId);
  //      thisColl.collTime = lambda * m_timeStep;
  //      theCollList.push_back(thisColl);
  //    }
  //  }
  //
  //  // The particle may encounter a wall which lies inside a neighboring cell.
  //  // The check is only carried out at the highest grid refinement level, as this
  //  // level is found near the solid walls and the check requires quite some time

  //    MFloat lambda;
  //
  //    MFloat centerCoords[MAX_SPACE_DIMENSIONS];
  //    MInt curCellId = (*partVectorElement).m_cellId;
  //    MFloat halfLength = m_lpt->c_cellLengthAtLevel(m_lpt->c_level(curCellId) + 1);
  //    MFloat diagCoords[MAX_SPACE_DIMENSIONS];
  //
  //    //now just c_coordinate!
  //
  //    // diagCoords temporarily holds the difference between particle and
  //    // cell-center coordinates
  //    for(MInt m = 0; m < nDim; m++) {
  //      diagCoords[m] = (*partVectorElement).m_position[m] - centerCoords[m];
  //    }
  //
  //    // direct neighbors in the direction of the relative displacement between
  //    // particle and cell center
  //    for(MInt d = 0; d < nDim; d++) {
  //      for(MInt n = m_fvSolver->a_identNghbrId((*partVectorElement).m_cellId * 2 * nDim + d);
  //          n < m_fvSolver->a_identNghbrId((*partVectorElement).m_cellId * 2 * nDim + d + 1);
  //          n++) {
  //        if(m_fvSolver->a_boundaryId(m_fvSolver->a_storeNghbrId(n)) > -1) {
  //          lambda = checkBndryCross(m_fvSolver->a_boundaryId(m_fvSolver->a_storeNghbrId(n)),
  //                                   partVectorElement, (*partVectorElement).m_diameter);
  //          if(lambda >= 0.0 && lambda <= 1.0) {
  //            thisColl.particle0 = partVectorElement;
  //            thisColl.particle1 = partVectorElement;
  //            thisColl.bndryId = m_fvSolver->a_boundaryId(m_fvSolver->a_storeNghbrId(n));
  //            thisColl.collTime = lambda * m_timeStep;
  //            theCollList.push_back(thisColl);
  //          }
  //        }
  //      }
  //
  //    // check the diagonal neighbor(s)
  //    for(MInt d = 0; d < nDim; d++) {
  //      diagCoords[d] /= fabs(diagCoords[d]);
  //      diagCoords[d] *= halfLength * 1.2;
  //    }
  //
  //
  //    MFloat tempDiagCoords[3];
  //
  //    for(MInt d = 0; d < 4; d++) {
  //      curCellId = (*partVectorElement).m_cellId;
  //      // cycle through directions xy (d=0), xz (1), yz (2), xyz (3)
  //      tempDiagCoords[0] = centerCoords[0] + float((d - 2) != 0) * diagCoords[0];
  //      tempDiagCoords[1] = centerCoords[1] + float((d - 1) != 0) * diagCoords[1];
  //      IF_CONSTEXPR(nDim == 2) {
  //        tempDiagCoords[2] = centerCoords[2] + float(d != 0) * diagCoords[2];
  //      }
  //
  //      checkCoordCellChange(m_solver, tempDiagCoords, centerCoords, curCellId);
  //
  //      if((curCellId != (*partVectorElement).m_cellId) &&
  //         (m_fvSolver->a_boundaryId(curCellId) > -1)) {
  //        lambda = checkBndryCross(m_fvSolver->a_boundaryId(curCellId),
  //                                 partVectorElement, (*partVectorElement).m_diameter);
  //        if(lambda >= 0.0 && lambda <= 1.0) {
  //          thisColl.particle0 = partVectorElement;
  //          thisColl.particle1 = partVectorElement;
  //          thisColl.bndryId = m_fvSolver->a_boundaryId(curCellId);
  //          thisColl.collTime = lambda * m_timeStep;
  //          theCollList.push_back(thisColl);
  //        }
  //      }
  //    }
  //  }
}

/** MFloat LPT::collisionCheck(MInt   i,
 *                                         MInt   j,
 *                                         MFloat currCollTime)
 *
 *  \brief Checks whether particles i and j collide
 *  after currCollTime, and returns collision time
 *
 *  version: December-09
 *  @author Rudie Kunnen
 */
template <MInt nDim>
MFloat ParticleCollision<nDim>::collisionCheck(partListIteratorConst<nDim> i, partListIteratorConst<nDim> j,
                                               MFloat currCollTime) {
  //   TRACE();

  MFloat collTime = -1.0;

  switch(m_collisionModel) {
    case 1:
    case 3:
    case 5: { // retroactive collision detection: search for overlap between particles
      MFloat currDist, collDist, squaredw = 0.0, squaredr0 = 0.0;
      MFloat wdotr0 = 0.0, tempa, tempb;
      // currDist = squared distance between centres of particles i and j
      // collDist = squared distance between centres of particles i and j
      //            at collision, i.e. the square of the sum of radii
      currDist = distance(i, j);
      collDist = POW2(0.5 * ((*i).m_diameter + (*j).m_diameter));

      // collision condition: particles overlap
      if(currDist < collDist) {
        // calculation of collision time, uses square of vector w,
        // square of vector r0, and inner product (w.r0)
        for(MInt k = 0; k < nDim; k++) {
          tempa = (((*i).m_position[k] - (*i).m_oldPos[k]) - ((*j).m_position[k] - (*j).m_oldPos[k])) / m_timeStep;
          tempb = (*i).m_oldPos[k] - (*j).m_oldPos[k];
          squaredw += POW2(tempa);
          squaredr0 += POW2(tempb);
          wdotr0 += tempa * tempb;
        }
        collTime = -(wdotr0 / squaredw) * (1.0 - sqrt(1.0 - squaredw * (squaredr0 - collDist) / POW2(wdotr0)));
      }
    } break;

    case 2:
    case 4:
    case 6: { // proactive collision detection: calculate time of collision and
      // compare with the particle time step
      MFloat bigK, squaredw = 0.0, squaredr0 = 0.0, wdotr0 = 0.0;
      MFloat tempa, tempb;

      for(MInt k = 0; k < nDim; k++) {
        tempa = (((*i).m_position[k] - (*i).m_oldPos[k]) - ((*j).m_position[k] - (*j).m_oldPos[k])) / m_timeStep;
        tempb = (*i).m_oldPos[k] - (*j).m_oldPos[k];
        squaredw += POW2(tempa);
        squaredr0 += POW2(tempb);
        wdotr0 += tempa * tempb;
      }

      // tempb holds the square of the separation distance upon collision
      // (sum of radii, squared)
      tempb = POW2(0.5 * ((*i).m_diameter + (*j).m_diameter));
      bigK = squaredw * (squaredr0 - tempb) / POW2(wdotr0);

      if(bigK <= 1.0) {
        // There are real solutions to the equation that gives the particle
        // collision times, i.e. there will be a collision given the current
        // particle paths.
        // tempa holds the time of minimal separation
        // squaredw and squaredr0 hold the two collision times
        bigK = sqrt(1.0 - bigK);
        tempa = -wdotr0 / squaredw;
        squaredw = tempa * (1.0 - bigK);
        squaredr0 = tempa * (1.0 + bigK);
        tempb = (squaredw < squaredr0) ? squaredw : squaredr0;
        // tempb holds the smaller of the two times
        if((tempb < m_timeStep) && (tempb >= currCollTime)) {
          collTime = tempb;
        }
      } // when bigK > 1.0 there are no real solutions: particles move without
      // (ever) touching (based on current velocity & position)
    } break;

    default: {
      mTerm(1, AT_, "Unknown particle collision type");
    }
  }

  // schedule the collision only when both particles are still valid and present within the domain
  if(!(*i).toBeDeleted() && !(*j).toBeDeleted()) {
    return collTime;
  } else {
    return -1.0;
  }
}

/** void LPT::collisionCheckSphereEllipsoid(ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 * \brief  Collision check for ellipsoids i2 and i3
 *         calls the corresponding proactive or retroactive method
 *
 * @author Christoph Siewert May-2013
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckSphereEllipsoid(
    partListIteratorConst<nDim> i1, ellipsListIteratorConst<nDim> i2) {
  TRACE();

  // check offset
  // particles before the plane x = m_particleCollOffset are not considered
  if((*i1).m_position[0] < m_offset) {
    return;
  }
  if((*i2).m_position[0] < m_offset) {
    return;
  }

  switch(m_collisionModel) {
    case 1:
    case 3:
    case 5: { // retroactive collision detection: search for overlap between particles
      //       collisionCheckSERetroActive(i1,i2);
      collisionCheckSEP(i1, i2);
    } break;

    case 2:
    case 4:
    case 6: { // proactive collision detection: calculate time of collision and
      // compare with the particle time step
      if(m_ellipsoidCCD) {
        collisionCheckSECCD(i1, i2);
      } else {
        collisionCheckSEP(i1, i2);
      }
    } break;

    default: {
      mTerm(1, AT_, "Unknown particle collision type");
    }
  }
}

/** void LPT::collisionCheckEllipsoidEllipsoid(ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 * \brief  Collision check for ellipsoids i2 and i3
 *         calls the corresponding proactive or retroactive method
 *
 * @author Christoph Siewert May-2013
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckEllipsoidEllipsoid(
    ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3) {
  TRACE();

  // check offset
  // particles before the plane x = m_particleCollOffset are not considered
  if((*i2).m_position[0] < m_offset) {
    return;
  }
  if((*i3).m_position[0] < m_offset) {
    return;
  }

  switch(m_collisionModel) {
    case 1:
    case 3:
    case 5: { // retroactive collision detection: search for overlap between particles
      //       collisionCheckEERetroActive(i2,i3);
      collisionCheckEEProActive(i2, i3);
    } break;

    case 2:
    case 4:
    case 6: { // proactive collision detection: calculate time of collision and
      // compare with the particle time step
      if(m_ellipsoidCCD) {
        collisionCheckEECCD(i2, i3);
      } else {
        collisionCheckEEProActive(i2, i3);
      }
    } break;

    default: {
      mTerm(1, AT_, "Unknown particle collision type");
    }
  }
}

/** void LPT::collisionCheckSEP((partListIteratorConst<nDim> i1, ellipsListIteratorConst<nDim> i2)
 * \brief  Collision check for sphere i1 and ellipsoid i2, specialized version of ellips ellips
 *         checks and stores collision event
 *         this method works proactive, but assumes that the closet distance between the centroids
 *         corresponds to the closest distance betweeen the two particles
 *
 * @author Christoph Siewert May-2013
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckSEP(partListIteratorConst<nDim> i1, ellipsListIteratorConst<nDim> i2) {
  TRACE();

  MFloat deltaXSq = F0;
  MFloat deltaVSq = F0;
  MFloat deltaVdeltaX = F0;
  MFloat dist[nDim];
  for(MInt d = 0; d < nDim; ++d) {
    dist[d] = i1->m_oldPos[d] - i2->m_oldPos[d];
    MFloat deltaV = ((i1->m_position[d] - i1->m_oldPos[d]) - (i2->m_position[d] - i2->m_oldPos[d])) / m_timeStep;
    deltaXSq += POW2(dist[d]);
    deltaVSq += POW2(deltaV);
    deltaVdeltaX += dist[d] * deltaV;
  }
  if(approx(deltaVSq, F0, MFloatEps)) { // no relative velocity, check for overlap only
    collisionCheckSERetroActive(i1, i2);
    return;
  }

  MFloat i1minR = i1->m_diameter / F2;
  MFloat& i1maxR = i1minR;
  MFloat i2minR = min(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat i2maxR = max(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat RminSq = POW2(i1minR + i2minR);
  MFloat RmaxSq = POW2(i1maxR + i2maxR);

  // as a first rough criterion, check the distance between the centroids
  // if larger than the sum of major axes, then
  // for sure there is no collision between the two
  MFloat bigKmax = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RmaxSq);
  if(bigKmax < F0) {
    return; // only complex solutions, no collision on current route
  }
  bigKmax = sqrt(bigKmax);
  MFloat tempt1 = (-deltaVdeltaX - bigKmax) / deltaVSq;
  MFloat tempt2 = (-deltaVdeltaX + bigKmax) / deltaVSq;
  MFloat t1max = min(tempt1, tempt2);
  MFloat t2max = max(tempt1, tempt2);
  if(t1max > m_timeStep) {
    return; // collision begin in the future
  }
  if(t2max < F0) {
    return; // collision end in the past
  }

  // if smaller than the sum of the minor axis, then
  //  a collision between the two is certain
  MFloat bigKmin = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RminSq);
  if(bigKmin >= F0) {
    bigKmin = sqrt(bigKmin);
    tempt1 = (-deltaVdeltaX - bigKmin) / deltaVSq;
    tempt2 = (-deltaVdeltaX + bigKmin) / deltaVSq;
    MFloat t1min = min(tempt1, tempt2);
    MFloat t2min = max(tempt1, tempt2);
    if((t1min <= m_timeStep) && (t2min >= F0)) {
      collQueueElemEllipsoid thisColl;
      thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + t1min;
      thisColl.part0 = (*i2).m_partId;
      thisColl.part1 = (*i1).m_partId;
      thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
      thisColl.semiMinorAxis1 = i1minR;
      thisColl.aspectRatio0 = (*i2).m_aspectRatio;
      thisColl.aspectRatio1 = -1.0;
      thisColl.dens0 = (*i2).m_densityRatio;
      thisColl.dens1 = (*i1).m_densityRatio;
      thisColl.collPosX = i1->m_oldPos[0] + (i1->m_position[0] - i1->m_oldPos[0]) * t1min / m_timeStep;
      thisColl.collPosY = i1->m_oldPos[1] + (i1->m_position[1] - i1->m_oldPos[1]) * t1min / m_timeStep;
      thisColl.collPosZ = i1->m_oldPos[2] + (i1->m_position[2] - i1->m_oldPos[2]) * t1min / m_timeStep;
      collQueueEllipsoid.push(thisColl);
      return;
    }
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur
  // consider distance between the bodies with the orientation at the closet distance
  // between the two centroids (might not be a good choose with the ellipsoids rotate very fast
  // while tranlate very slowly
  MFloat distTminSq = F0;
  MFloat a[3], b[3], c[3];
  MFloat l2[3] = {1.0, 0.0, 0.0}, m2[3] = {0.0, 1.0, 0.0}, n2[3] = {0.0, 0.0, 1.0};
  MFloat tmin = (-deltaVdeltaX / deltaVSq) / m_timeStep;
  if((tmin < F0) || (tmin > F1)) {
    MFloat partQuatA[4];
    if(tmin < F0) {
      for(MInt i = 0; i < 4; ++i) {
        partQuatA[i] = i2->m_oldQuaternion[i];
      }
      distTminSq = deltaXSq;
    } else {
      for(MInt i = 0; i < 4; ++i) {
        partQuatA[i] = i2->m_quaternion[i];
      }
      for(MInt i = 0; i < nDim; ++i) {
        distTminSq += POW2(i1->m_position[i] - i2->m_position[i]);
      }
    }
    // Calculate the transformation matrix A of orientation of ellipsoid i2
    a[0] = 1.0 - 2.0 * (partQuatA[2] * partQuatA[2] + partQuatA[3] * partQuatA[3]);
    a[1] = 2.0 * (partQuatA[1] * partQuatA[2] + partQuatA[0] * partQuatA[3]);
    a[2] = 2.0 * (partQuatA[1] * partQuatA[3] - partQuatA[0] * partQuatA[2]);
    b[0] = 2.0 * (partQuatA[1] * partQuatA[2] - partQuatA[0] * partQuatA[3]);
    b[1] = 1.0 - 2.0 * (partQuatA[1] * partQuatA[1] + partQuatA[3] * partQuatA[3]);
    b[2] = 2.0 * (partQuatA[2] * partQuatA[3] + partQuatA[0] * partQuatA[1]);
    c[0] = 2.0 * (partQuatA[1] * partQuatA[3] + partQuatA[0] * partQuatA[2]);
    c[1] = 2.0 * (partQuatA[2] * partQuatA[3] - partQuatA[0] * partQuatA[1]);
    c[2] = 1.0 - 2.0 * (partQuatA[1] * partQuatA[1] + partQuatA[2] * partQuatA[2]);
  } else {
    //         slerp(i2->m_oldQuaternion, i2->m_quaternion, tmin, partQuatA, 4);
    // slerp between old and new orientation;
    const MFloat* e;
    e = i2->m_oldQuaternion.data();
    a[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    a[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    a[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    e = i2->m_quaternion.data();
    b[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    b[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    b[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    slerp(a, b, tmin, c, 3);

    // generate orthonormal basis using Gram-Schmidt and cross-product
    for(MInt dim = 0; dim < 3; ++dim) {
      a[dim] = c[dim];
    }
    a[0] -= F1;
    normalize(a, 3);
    MFloat innerProd = F0;
    for(MInt dim = 0; dim < 3; ++dim) {
      innerProd += c[dim] * a[dim];
    }
    for(MInt dim = 0; dim < 3; ++dim) {
      b[dim] = a[dim] - innerProd * c[dim];
    }
    normalize(b, 3);
    a[0] = b[1] * c[2] - b[2] * c[1];
    a[1] = b[2] * c[0] - b[0] * c[2];
    a[2] = b[0] * c[1] - b[1] * c[0];
    normalize(a, 3);

    for(MInt dim = 0; dim < nDim; ++dim) {
      distTminSq += POW2((i1->m_position[dim] - i2->m_position[dim]) * tmin
                         + (i1->m_oldPos[dim] - i2->m_oldPos[dim]) * (1 - tmin));
    }
  }

  // next, apply the method proposed by Zheng et al., Phys. Rev. E 79, 057702 (2009)
  // to determine the minimal distance between the two bodies
  // at the current orientation/separation

  // determine minimal distance at contact between the particles
  // distance = separation vector of midpoints
  // l1, m1, n1 = principal orientation vectors of ellipsoid i2
  // l2, m2, n2 = principal orientation vectors of ellipsoid i3
  // store the minimal distance in largestCollDistSquared
  EllipsoidDistance theDistance(dist, a, b, c, l2, m2, n2, (*i2).m_semiMinorAxis, (*i2).m_semiMinorAxis,
                                (*i2).m_semiMinorAxis * (*i2).m_aspectRatio, i1minR, i1minR, i1minR);
  MFloat largestCollDistSquared = theDistance.ellipsoids();
  largestCollDistSquared = POW2(largestCollDistSquared);

  if(distTminSq < largestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + tmin*m_particleTimeStep;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i1).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = i1minR;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = -1.0;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i1).m_densityRatio;
    thisColl.collPosX = i1->m_oldPos[0] + (i1->m_position[0] - i1->m_oldPos[0]) * tmin;
    thisColl.collPosY = i1->m_oldPos[1] + (i1->m_position[1] - i1->m_oldPos[1]) * tmin;
    thisColl.collPosZ = i1->m_oldPos[2] + (i1->m_position[2] - i1->m_oldPos[2]) * tmin;
    collQueueEllipsoid.push(thisColl);
  }
}

/** void LPT::collisionCheckSERetroActive(partListIteratorConst<nDim> i1,
 *                              ellipsListIteratorConst<nDim> i2)
 * \brief  Collision check for sphere i1 and ellipsoid i2
 *         based on overlap criterion, and store collision event when found
 *
 * @author Rudie Kunnen, Feb-2011
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckSERetroActive(
    partListIteratorConst<nDim> i1, ellipsListIteratorConst<nDim> i2) {
  TRACE();

  MFloat dist[nDim];
  MFloat totalDistanceSquared = 0.0;

  for(MInt i = 0; i < nDim; i++) {
    dist[i] = (*i2).m_position[i] - (*i1).m_position[i];
    totalDistanceSquared += POW2(dist[i]);
  }

  MFloat Length1 = 0.5 * (*i1).m_diameter;
  MFloat minLength2, maxLength2;
  minLength2 = (*i2).m_semiMinorAxis * (*i2).m_aspectRatio;
  if(minLength2 < (*i2).m_semiMinorAxis) {
    maxLength2 = (*i2).m_semiMinorAxis;
  } else {
    maxLength2 = minLength2;
    minLength2 = (*i2).m_semiMinorAxis;
  }

  // as a first rough criterion, check the distance between the centroids
  // if larger than radius_sphere + semi_major_axis_ellipsoid, then
  // for sure there is no collision between the two
  MFloat largestCollDistSquared = POW2(Length1 + maxLength2);
  if(totalDistanceSquared > largestCollDistSquared) {
    return;
  }
  MFloat smallestCollDistSquared = POW2(Length1 + minLength2);
  if(totalDistanceSquared < smallestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i1).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = Length1;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = -1.0;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i1).m_densityRatio;
    thisColl.collPosX = i1->m_position[0];
    thisColl.collPosY = i1->m_position[1];
    thisColl.collPosZ = i1->m_position[2];
    collQueueEllipsoid.push(thisColl);
    return;
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur
  // consider distace between the bodies

  MFloat A[3][3];
  //   MFloat invA[3][3];

  // Calculate the transformation matrix A of the ellipsoid orientation
  A[0][0] = 1.0 - 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[2] + (*i2).m_quaternion[3] * (*i2).m_quaternion[3]);
  A[0][1] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[2] + (*i2).m_quaternion[0] * (*i2).m_quaternion[3]);
  A[0][2] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[3] - (*i2).m_quaternion[0] * (*i2).m_quaternion[2]);
  A[1][0] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[2] - (*i2).m_quaternion[0] * (*i2).m_quaternion[3]);
  A[1][1] = 1.0 - 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[1] + (*i2).m_quaternion[3] * (*i2).m_quaternion[3]);
  A[1][2] = 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[3] + (*i2).m_quaternion[0] * (*i2).m_quaternion[1]);
  A[2][0] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[3] + (*i2).m_quaternion[0] * (*i2).m_quaternion[2]);
  A[2][1] = 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[3] - (*i2).m_quaternion[0] * (*i2).m_quaternion[1]);
  A[2][2] = 1.0 - 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[1] + (*i2).m_quaternion[2] * (*i2).m_quaternion[2]);
  // and its inverse (inverse == transpose for rotation matrix!)
  //   for (MInt i = 0; i < 3; i++)
  //     for (MInt j = 0; j < 3; j++)
  //       invA[i][j] = A[j][i];
  //
  //   // calculate the matrix product invA * X_A * A, is stored in A
  //   // [X_A is a diagonal matrix with elements: semi_axis_length^2]
  //   MFloat minorAxis2 = POW2(minLength2);
  //   MFloat majorAxis2 = POW2(maxLength2);
  //   for (MInt i = 0; i < 3; i++)
  //   { // this is the right-hand side of the product: X_A * A
  //     // stored in A
  //     A[0][i] *= minorAxis2;
  //     A[1][i] *= minorAxis2;
  //     A[2][i] *= majorAxis2;
  //   }
  //   // left-hand-side part, again stored in A
  //   matrixMultiplyLeft(invA, A);
  //
  //   // prepare input vectors
  //   MFloat l1[3], m1[3], n1[3];
  //   MFloat l2[3] = {1.0, 0.0, 0.0}, m2[3] = {0.0, 1.0, 0.0}, n2[3] = {0.0, 0.0, 1.0};
  //   for (MInt i = 0; i < 3; i++)
  //   {
  //     l1[i] = A[i][0];
  //     m1[i] = A[i][1];
  //     n1[i] = A[i][2];
  //   }

  // next, apply the method proposed by Zheng et al., Phys. Rev. E 79, 057702 (2009)
  // to determine the minimal distance between the two bodies
  // at the current orientation/separation

  // prepare input vectors
  MFloat l1[3], m1[3], n1[3];
  MFloat l2[3] = {1.0, 0.0, 0.0}, m2[3] = {0.0, 1.0, 0.0}, n2[3] = {0.0, 0.0, 1.0};
  for(MInt i = 0; i < 3; i++) {
    l1[i] = A[0][i];
    m1[i] = A[1][i];
    n1[i] = A[2][i];
  }

  // determine minimal distance at contact between the particles
  // distance = separation vector of midpoints
  // l1, m1, n1 = principal orientation vectors of ellipsoid
  // l2, m2, n2 = orientation vectors of sphere (consider sphere as an ellipsoid with AR = 1)
  // store the minimal distance in largestCollDistSquared
  EllipsoidDistance theDistance(dist, l1, m1, n1, l2, m2, n2, (*i2).m_semiMinorAxis, (*i2).m_semiMinorAxis,
                                (*i2).m_semiMinorAxis * (*i2).m_aspectRatio, Length1, Length1, Length1);
  largestCollDistSquared = theDistance.ellipsoids();
  largestCollDistSquared = POW2(largestCollDistSquared);

  if(totalDistanceSquared < largestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i1).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = Length1;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = -1.0;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i1).m_densityRatio;
    thisColl.collPosX = i1->m_position[0];
    thisColl.collPosY = i1->m_position[1];
    thisColl.collPosZ = i1->m_position[2];
    collQueueEllipsoid.push(thisColl);
  }
}

/** void LPT::collisionCheckEEProActive(ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 * \brief  Collision check for ellipsoids i2 and i3
 *         checks and stores collision event
 *         this method works proactive, but assumes that the closet distance between the centroids
 *         corresponds to the closest distance betweeen the two ellipsoids
 *
 * @author Christoph Siewert May-2013
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckEEProActive(
    ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3) {
  TRACE();


  MFloat deltaXSq = F0;
  MFloat deltaVSq = F0;
  MFloat deltaVdeltaX = F0;
  MFloat dist[nDim];
  for(MInt d = 0; d < nDim; ++d) {
    dist[d] = i2->m_oldPos[d] - i3->m_oldPos[d];
    MFloat deltaV = ((i2->m_position[d] - i2->m_oldPos[d]) - (i3->m_position[d] - i3->m_oldPos[d])) / m_timeStep;
    deltaXSq += POW2(dist[d]);
    deltaVSq += POW2(deltaV);
    deltaVdeltaX += dist[d] * deltaV;
  }
  if(approx(deltaVSq, F0, MFloatEps)) { // no relative velocity, check for overlap only
    collisionCheckEERetroActive(i2, i3);
    return;
  }

  MFloat i2minR = min(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat i2maxR = max(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat i3minR = min(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat i3maxR = max(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat RminSq = POW2(i2minR + i3minR);
  MFloat RmaxSq = POW2(i2maxR + i3maxR);

  // as a first rough criterion, check the distance between the centroids
  // if larger than the sum of major axes, then
  // for sure there is no collision between the two
  MFloat bigKmax = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RmaxSq);
  if(bigKmax < F0) {
    return; // only complex solutions, no collision on current route
  }
  bigKmax = sqrt(bigKmax);
  MFloat tempt1 = (-deltaVdeltaX - bigKmax) / deltaVSq;
  MFloat tempt2 = (-deltaVdeltaX + bigKmax) / deltaVSq;
  MFloat t1max = min(tempt1, tempt2);
  MFloat t2max = max(tempt1, tempt2);
  if(t1max > m_timeStep) {
    return; // collision begin in the future
  }
  if(t2max < F0) {
    return; // collision end in the past
  }

  // if smaller than the sum of the minor axis, then
  //  a collision between the two is certain
  MFloat bigKmin = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RminSq);
  if(bigKmin >= F0) {
    bigKmin = sqrt(bigKmin);
    tempt1 = (-deltaVdeltaX - bigKmin) / deltaVSq;
    tempt2 = (-deltaVdeltaX + bigKmin) / deltaVSq;
    MFloat t1min = min(tempt1, tempt2);
    MFloat t2min = max(tempt1, tempt2);
    if((t1min <= m_timeStep) && (t2min >= F0)) {
      collQueueElemEllipsoid thisColl;
      thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + t1min;
      thisColl.part0 = (*i2).m_partId;
      thisColl.part1 = (*i3).m_partId;
      thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
      thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
      thisColl.aspectRatio0 = (*i2).m_aspectRatio;
      thisColl.aspectRatio1 = (*i3).m_aspectRatio;
      thisColl.dens0 = (*i2).m_densityRatio;
      thisColl.dens1 = (*i3).m_densityRatio;
      thisColl.collPosX = i2->m_oldPos[0] + (i2->m_position[0] - i2->m_oldPos[0]) * t1min / m_timeStep;
      thisColl.collPosY = i2->m_oldPos[1] + (i2->m_position[1] - i2->m_oldPos[1]) * t1min / m_timeStep;
      thisColl.collPosZ = i2->m_oldPos[2] + (i2->m_position[2] - i2->m_oldPos[2]) * t1min / m_timeStep;
      collQueueEllipsoid.push(thisColl);
      return;
    }
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur
  // consider distance between the bodies with the orientation at the closet distance
  // between the two centroids (might not be a good choose with the ellipsoids rotate very fast
  // while tranlate very slowly
  MFloat a1[3], b1[3], c1[3];
  MFloat a2[3], b2[3], c2[3];
  MFloat distTminSq;
  MFloat tmin = (-deltaVdeltaX / deltaVSq) / m_timeStep;
  if((tmin < F0) || (tmin > F1)) {
    MFloat partQuatA[4], partQuatB[4];
    if(tmin < F0) {
      for(MInt i = 0; i < 4; ++i) {
        partQuatA[i] = i2->m_oldQuaternion[i];
        partQuatB[i] = i3->m_oldQuaternion[i];
      }
      distTminSq = deltaXSq;
    } else {
      for(MInt i = 0; i < 4; ++i) {
        partQuatA[i] = i2->m_quaternion[i];
        partQuatB[i] = i3->m_quaternion[i];
      }
      distTminSq = F0;
      for(MInt i = 0; i < nDim; ++i) {
        distTminSq += POW2(i2->m_position[i] - i3->m_position[i]);
      }
    }
    // Calculate the transformation matrix A of orientation of ellipsoid i2
    a1[0] = 1.0 - 2.0 * (partQuatA[2] * partQuatA[2] + partQuatA[3] * partQuatA[3]);
    a1[1] = 2.0 * (partQuatA[1] * partQuatA[2] + partQuatA[0] * partQuatA[3]);
    a1[2] = 2.0 * (partQuatA[1] * partQuatA[3] - partQuatA[0] * partQuatA[2]);
    b1[0] = 2.0 * (partQuatA[1] * partQuatA[2] - partQuatA[0] * partQuatA[3]);
    b1[1] = 1.0 - 2.0 * (partQuatA[1] * partQuatA[1] + partQuatA[3] * partQuatA[3]);
    b1[2] = 2.0 * (partQuatA[2] * partQuatA[3] + partQuatA[0] * partQuatA[1]);
    c1[0] = 2.0 * (partQuatA[1] * partQuatA[3] + partQuatA[0] * partQuatA[2]);
    c1[1] = 2.0 * (partQuatA[2] * partQuatA[3] - partQuatA[0] * partQuatA[1]);
    c1[2] = 1.0 - 2.0 * (partQuatA[1] * partQuatA[1] + partQuatA[2] * partQuatA[2]);

    // Calculate the transformation matrix B of orientation of ellipsoid i3
    a2[0] = 1.0 - 2.0 * (partQuatB[2] * partQuatB[2] + partQuatB[3] * partQuatB[3]);
    a2[1] = 2.0 * (partQuatB[1] * partQuatB[2] + partQuatB[0] * partQuatB[3]);
    a2[2] = 2.0 * (partQuatB[1] * partQuatB[3] - partQuatB[0] * partQuatB[2]);
    b2[0] = 2.0 * (partQuatB[1] * partQuatB[2] - partQuatB[0] * partQuatB[3]);
    b2[1] = 1.0 - 2.0 * (partQuatB[1] * partQuatB[1] + partQuatB[3] * partQuatB[3]);
    b2[2] = 2.0 * (partQuatB[2] * partQuatB[3] + partQuatB[0] * partQuatB[1]);
    c2[0] = 2.0 * (partQuatB[1] * partQuatB[3] + partQuatB[0] * partQuatB[2]);
    c2[1] = 2.0 * (partQuatB[2] * partQuatB[3] - partQuatB[0] * partQuatB[1]);
    c2[2] = 1.0 - 2.0 * (partQuatB[1] * partQuatB[1] + partQuatB[2] * partQuatB[2]);
  } else {
    //         slerp(i2->m_oldQuaternion, i2->m_quaternion, tmin, partQuatA, 4);
    //         slerp(i3->m_oldQuaternion, i3->m_quaternion, tmin, partQuatB, 4);
    const MFloat* e;
    e = i2->m_oldQuaternion.data();
    a1[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    a1[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    a1[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    e = i2->m_quaternion.data();
    b1[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    b1[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    b1[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    slerp(a1, b1, tmin, c1, 3);

    // generate orthonormal basis using Gram-Schmidt and cross-product
    for(MInt dim = 0; dim < 3; ++dim) {
      a1[dim] = c1[dim];
    }
    a1[0] -= F1;
    normalize(a1, 3);
    MFloat innerProd = F0;
    for(MInt dim = 0; dim < 3; ++dim) {
      innerProd += c1[dim] * a1[dim];
    }
    for(MInt dim = 0; dim < 3; ++dim) {
      b1[dim] = a1[dim] - innerProd * c1[dim];
    }
    normalize(b1, 3);
    a1[0] = b1[1] * c1[2] - b1[2] * c1[1];
    a1[1] = b1[2] * c1[0] - b1[0] * c1[2];
    a1[2] = b1[0] * c1[1] - b1[1] * c1[0];
    normalize(a1, 3);
    // second ellipsoid
    e = i3->m_oldQuaternion.data();
    a2[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    a2[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    a2[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    e = i3->m_quaternion.data();
    b2[0] = 2.0 * (e[1] * e[3] + e[0] * e[2]);
    b2[1] = 2.0 * (e[2] * e[3] - e[0] * e[1]);
    b2[2] = 1.0 - 2.0 * (e[1] * e[1] + e[2] * e[2]);
    slerp(a2, b2, tmin, c2, 3);

    // generate orthonormal basis using Gram-Schmidt and cross-product
    for(MInt dim = 0; dim < nDim; ++dim) {
      a2[dim] = c2[dim];
    }
    a2[0] -= F1;
    normalize(a2, 3);
    innerProd = F0;
    for(MInt dim = 0; dim < nDim; ++dim) {
      innerProd += c2[dim] * a2[dim];
    }
    for(MInt dim = 0; dim < nDim; ++dim) {
      b2[dim] = a2[dim] - innerProd * c2[dim];
    }
    normalize(b2, 3);
    a2[0] = b2[1] * c2[2] - b2[2] * c2[1];
    a2[1] = b2[2] * c2[0] - b2[0] * c2[2];
    a2[2] = b2[0] * c2[1] - b2[1] * c2[0];
    normalize(a2, 3);

    distTminSq = F0;
    for(MInt dim = 0; dim < nDim; ++dim) {
      distTminSq += POW2((i2->m_position[dim] - i3->m_position[dim]) * tmin
                         + (i2->m_oldPos[dim] - i3->m_oldPos[dim]) * (1 - tmin));
    }
  }

  // next, apply the method proposed by Zheng et al., Phys. Rev. E 79, 057702 (2009)
  // to determine the minimal distance between the two bodies
  // at the current orientation/separation

  // determine minimal distance at contact between the particles
  // distance = separation vector of midpoints
  // l1, m1, n1 = principal orientation vectors of ellipsoid i2
  // l2, m2, n2 = principal orientation vectors of ellipsoid i3
  // store the minimal distance in largestCollDistSquared
  EllipsoidDistance theDistance(dist, a1, b1, c1, a2, b2, c2, (*i2).m_semiMinorAxis, (*i2).m_semiMinorAxis,
                                (*i2).m_semiMinorAxis * (*i2).m_aspectRatio, (*i3).m_semiMinorAxis,
                                (*i3).m_semiMinorAxis, (*i3).m_semiMinorAxis * (*i3).m_aspectRatio);
  MFloat largestCollDistSquared = theDistance.ellipsoids();
  largestCollDistSquared = POW2(largestCollDistSquared);

  if(distTminSq < largestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + tmin*m_particleTimeStep;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i3).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = (*i3).m_aspectRatio;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i3).m_densityRatio;
    thisColl.collPosX = i2->m_oldPos[0] + (i2->m_position[0] - i2->m_oldPos[0]) * tmin;
    thisColl.collPosY = i2->m_oldPos[1] + (i2->m_position[1] - i2->m_oldPos[1]) * tmin;
    thisColl.collPosZ = i2->m_oldPos[2] + (i2->m_position[2] - i2->m_oldPos[2]) * tmin;
    collQueueEllipsoid.push(thisColl);
  }
}

/** void LPT::collisionCheckEERetroActive(
 *   ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 * \brief  Collision check for ellipsoids i2 and i3
 *         based on overlap criterion, and store collision event when found
 *
 * @author Rudie Kunnen, Feb-2011
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckEERetroActive(
    ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3) {
  TRACE();

  MFloat dist[nDim];
  MFloat totalDistanceSquared = 0.0;

  // check offset
  // particles before the plane x = m_particleCollOffset are not considered
  if((*i2).m_position[0] < m_offset) {
    return;
  }
  if((*i3).m_position[0] < m_offset) {
    return;
  }

  for(MInt i = 0; i < nDim; i++) {
    dist[i] = (*i3).m_position[i] - (*i2).m_position[i];
    totalDistanceSquared += POW2(dist[i]);
  }

  MFloat minLength2, maxLength2, minLength3, maxLength3;
  minLength2 = (*i2).m_semiMinorAxis * (*i2).m_aspectRatio;
  if(minLength2 < (*i2).m_semiMinorAxis) {
    maxLength2 = (*i2).m_semiMinorAxis;
  } else {
    maxLength2 = minLength2;
    minLength2 = (*i2).m_semiMinorAxis;
  }
  minLength3 = (*i3).m_semiMinorAxis * (*i3).m_aspectRatio;
  if(minLength3 < (*i3).m_semiMinorAxis) {
    maxLength3 = (*i3).m_semiMinorAxis;
  } else {
    maxLength3 = minLength3;
    minLength3 = (*i3).m_semiMinorAxis;
  }

  // as a first rough criterion, check the distance between the centroids
  // if larger than the sum of major axes, then
  // for sure there is no collision between the two
  MFloat largestCollDistSquared = POW2(maxLength2 + maxLength3);
  if(totalDistanceSquared > largestCollDistSquared) {
    return;
  }
  // if smaller than the sum of the minor axis, then
  // for sure there is a collision between the two
  MFloat smallestCollDistSquared = POW2(minLength2 + minLength3);
  if(totalDistanceSquared < smallestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i3).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = (*i3).m_aspectRatio;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i3).m_densityRatio;
    thisColl.collPosX = i2->m_position[0];
    thisColl.collPosY = i2->m_position[1];
    thisColl.collPosZ = i2->m_position[2];
    collQueueEllipsoid.push(thisColl);
    return;
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur
  // consider distace between the bodies
  MFloat A[3][3], B[3][3];
  //   MFloat invA[3][3], invB[3][3];

  // Calculate the transformation matrix A of orientation of ellipsoid i2
  A[0][0] = 1.0 - 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[2] + (*i2).m_quaternion[3] * (*i2).m_quaternion[3]);
  A[0][1] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[2] + (*i2).m_quaternion[0] * (*i2).m_quaternion[3]);
  A[0][2] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[3] - (*i2).m_quaternion[0] * (*i2).m_quaternion[2]);
  A[1][0] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[2] - (*i2).m_quaternion[0] * (*i2).m_quaternion[3]);
  A[1][1] = 1.0 - 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[1] + (*i2).m_quaternion[3] * (*i2).m_quaternion[3]);
  A[1][2] = 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[3] + (*i2).m_quaternion[0] * (*i2).m_quaternion[1]);
  A[2][0] = 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[3] + (*i2).m_quaternion[0] * (*i2).m_quaternion[2]);
  A[2][1] = 2.0 * ((*i2).m_quaternion[2] * (*i2).m_quaternion[3] - (*i2).m_quaternion[0] * (*i2).m_quaternion[1]);
  A[2][2] = 1.0 - 2.0 * ((*i2).m_quaternion[1] * (*i2).m_quaternion[1] + (*i2).m_quaternion[2] * (*i2).m_quaternion[2]);
  // and its inverse (inverse == transpose for rotation matrix!)
  //   for (MInt i = 0; i < 3; i++)
  //     for (MInt j = 0; j < 3; j++)
  //       invA[i][j] = A[j][i];

  // Calculate the transformation matrix B of orientation of ellipsoid i3
  B[0][0] = 1.0 - 2.0 * ((*i3).m_quaternion[2] * (*i3).m_quaternion[2] + (*i3).m_quaternion[3] * (*i3).m_quaternion[3]);
  B[0][1] = 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[2] + (*i3).m_quaternion[0] * (*i3).m_quaternion[3]);
  B[0][2] = 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[3] - (*i3).m_quaternion[0] * (*i3).m_quaternion[2]);
  B[1][0] = 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[2] - (*i3).m_quaternion[0] * (*i3).m_quaternion[3]);
  B[1][1] = 1.0 - 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[1] + (*i3).m_quaternion[3] * (*i3).m_quaternion[3]);
  B[1][2] = 2.0 * ((*i3).m_quaternion[2] * (*i3).m_quaternion[3] + (*i3).m_quaternion[0] * (*i3).m_quaternion[1]);
  B[2][0] = 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[3] + (*i3).m_quaternion[0] * (*i3).m_quaternion[2]);
  B[2][1] = 2.0 * ((*i3).m_quaternion[2] * (*i3).m_quaternion[3] - (*i3).m_quaternion[0] * (*i3).m_quaternion[1]);
  B[2][2] = 1.0 - 2.0 * ((*i3).m_quaternion[1] * (*i3).m_quaternion[1] + (*i3).m_quaternion[2] * (*i3).m_quaternion[2]);
  // and its inverse (inverse == transpose for rotation matrix!)
  //   for (MInt i = 0; i < 3; i++)
  //     for (MInt j = 0; j < 3; j++)
  //       invB[i][j] = B[j][i];

  //   // calculate the matrix product invA * X_A * A, is stored in A
  //   // [X_A is a diagonal matrix with elements: semi_axis_length^2]
  //   MFloat minorAxis2 = POW2(minLength2);
  //   MFloat majorAxis2 = POW2(maxLength3);
  //   for (MInt i = 0; i < 3; i++)
  //   { // this is the right-hand side of the product: X_A * A
  //     // stored in A
  //     A[0][i] *= minorAxis2;
  //     A[1][i] *= minorAxis2;
  //     A[2][i] *= majorAxis2;
  //   }
  //   // left-hand-side part, again stored in A
  //   matrixMultiplyLeft(invA, A);
  //
  //   // calculate the matrix product invB * X_B * B, is stored in B
  //   // [X_B is a diagonal matrix with elements: semi_axis_length^2]
  //   minorAxis2 = POW2(minLength3);
  //   majorAxis2 = POW2(maxLength3);
  //   for (MInt i = 0; i < 3; i++)
  //   { // this is the right-hand side of the product: X_B * B
  //     // stored in B
  //     B[0][i] *= minorAxis2;
  //     B[1][i] *= minorAxis2;
  //     B[2][i] *= majorAxis2;
  //   }
  //   // left-hand-side part, again stored in B
  //   matrixMultiplyLeft(invB, B);
  //
  //   // next, apply the method proposed by Zheng et al., Phys. Rev. E 79, 057702 (2009)
  //   // to determine the minimal distance between the two bodies
  //   // at the current orientation/separation
  //
  //   // prepare input vectors
  //   MFloat l1[3], m1[3], n1[3];
  //   MFloat l2[3], m2[3], n2[3];
  //   for (MInt i = 0; i < 3; i++)
  //   {
  //     l1[i] = A[i][0];
  //     m1[i] = A[i][1];
  //     n1[i] = A[i][2];
  //     l2[i] = B[i][0];
  //     m2[i] = B[i][1];
  //     n2[i] = B[i][2];
  //   }

  // next, apply the method proposed by Zheng et al., Phys. Rev. E 79, 057702 (2009)
  // to determine the minimal distance between the two bodies
  // at the current orientation/separation
  MFloat l1[3], m1[3], n1[3];
  MFloat l2[3], m2[3], n2[3];
  for(MInt i = 0; i < 3; i++) {
    l1[i] = A[0][i];
    m1[i] = A[1][i];
    n1[i] = A[2][i];
    l2[i] = B[0][i];
    m2[i] = B[1][i];
    n2[i] = B[2][i];
  }

  // determine minimal distance at contact between the particles
  // distance = separation vector of midpoints
  // l1, m1, n1 = principal orientation vectors of ellipsoid i2
  // l2, m2, n2 = principal orientation vectors of ellipsoid i3
  // store the minimal distance in largestCollDistSquared
  EllipsoidDistance theDistance(dist, l1, m1, n1, l2, m2, n2, (*i2).m_semiMinorAxis, (*i2).m_semiMinorAxis,
                                (*i2).m_semiMinorAxis * (*i2).m_aspectRatio, (*i3).m_semiMinorAxis,
                                (*i3).m_semiMinorAxis, (*i3).m_semiMinorAxis * (*i3).m_aspectRatio);
  largestCollDistSquared = theDistance.ellipsoids();
  largestCollDistSquared = POW2(largestCollDistSquared);

  if(totalDistanceSquared < largestCollDistSquared) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i3).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = (*i3).m_aspectRatio;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i3).m_densityRatio;
    thisColl.collPosX = i2->m_position[0];
    thisColl.collPosY = i2->m_position[1];
    thisColl.collPosZ = i2->m_position[2];
    collQueueEllipsoid.push(thisColl);
  }
}

/**  void LPT::collisionCheckSECCD(partListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 *   \brief  checks continuously for the overlaps of the two ellipsoids i2 and i3 assuming that they do not rotate
 * during the timestep after Choi et al., IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 15, NO. 2,
 * MARCH/APRIL 2009, 311 Continuous Collision Detection for Ellipsoids
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckEECCD(ellipsListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3) {
  TRACE();

  const MInt dim3 = 3;

  MFloat deltaXSq = F0;
  MFloat deltaVSq = F0;
  MFloat deltaVdeltaX = F0;
  MFloat dist[dim3];
  for(MInt d = 0; d < dim3; ++d) {
    dist[d] = i2->m_oldPos[d] - i3->m_oldPos[d];
    MFloat deltaV = ((i2->m_position[d] - i2->m_oldPos[d]) - (i3->m_position[d] - i3->m_oldPos[d])) / m_timeStep;
    deltaXSq += POW2(dist[d]);
    deltaVSq += POW2(deltaV);
    deltaVdeltaX += dist[d] * deltaV;
  }
  if(approx(deltaVSq, F0, MFloatEps)) { // no relative velocity, check for overlap only
    collisionCheckEERetroActive(i2, i3);
    return;
  }

  MFloat i2minR = min(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat i2maxR = max(i2->m_semiMinorAxis, i2->m_semiMinorAxis * i2->m_aspectRatio);
  MFloat i3minR = min(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat i3maxR = max(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat RminSq = POW2(i2minR + i3minR);
  MFloat RmaxSq = POW2(i2maxR + i3maxR);

  // as a first rough criterion, check the distance between the centroids
  // if larger than the sum of major axes, then
  // for sure there is no collision between the two
  MFloat bigKmax = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RmaxSq);
  if(bigKmax < F0) {
    return; // only complex solutions, no collision on current route
  }
  bigKmax = sqrt(bigKmax);
  MFloat tempt1 = (-deltaVdeltaX - bigKmax) / deltaVSq;
  MFloat tempt2 = (-deltaVdeltaX + bigKmax) / deltaVSq;
  MFloat t1max = min(tempt1, tempt2);
  MFloat t2max = max(tempt1, tempt2);
  if(t1max > m_timeStep) {
    return; // collision begin in the future
  }
  if(t2max < F0) {
    return; // collision end in the past
  }

  // if smaller than the sum of the minor axis, then
  //  a collision between the two is certain
  MFloat bigKmin = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RminSq);
  if(bigKmin >= F0) {
    bigKmin = sqrt(bigKmin);
    tempt1 = (-deltaVdeltaX - bigKmin) / deltaVSq;
    tempt2 = (-deltaVdeltaX + bigKmin) / deltaVSq;
    MFloat t1min = min(tempt1, tempt2);
    MFloat t2min = max(tempt1, tempt2);
    if((t1min <= m_timeStep) && (t2min >= F0)) {
      collQueueElemEllipsoid thisColl;
      thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + t1min;
      thisColl.part0 = (*i2).m_partId;
      thisColl.part1 = (*i3).m_partId;
      thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
      thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
      thisColl.aspectRatio0 = (*i2).m_aspectRatio;
      thisColl.aspectRatio1 = (*i3).m_aspectRatio;
      thisColl.dens0 = (*i2).m_densityRatio;
      thisColl.dens1 = (*i3).m_densityRatio;
      thisColl.collPosX = i2->m_oldPos[0] + (i2->m_position[0] - i2->m_oldPos[0]) * t1min / m_timeStep;
      thisColl.collPosY = i2->m_oldPos[1] + (i2->m_position[1] - i2->m_oldPos[1]) * t1min / m_timeStep;
      thisColl.collPosZ = i2->m_oldPos[2] + (i2->m_position[2] - i2->m_oldPos[2]) * t1min / m_timeStep;
      collQueueEllipsoid.push(thisColl);
      return;
    }
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur


  MFloat e0_1 = i2->m_quaternion[0];
  MFloat e1_1 = i2->m_quaternion[1];
  MFloat e2_1 = i2->m_quaternion[2];
  MFloat e3_1 = i2->m_quaternion[3];

  MFloat e0_2 = i3->m_quaternion[0];
  MFloat e1_2 = i3->m_quaternion[1];
  MFloat e2_2 = i3->m_quaternion[2];
  MFloat e3_2 = i3->m_quaternion[3];

  MFloat MA[16];
  MA[0] = ((e0_1 * e0_1 + e1_1 * e1_1) - e2_1 * e2_1) - e3_1 * e3_1;
  MA[4] = 2.0 * e1_1 * e2_1 - 2.0 * e0_1 * e3_1;
  MA[8] = 2.0 * e0_1 * e2_1 + 2.0 * e1_1 * e3_1;
  MA[12] = i2->m_oldPos[0];
  MA[1] = 2.0 * e0_1 * e3_1 + 2.0 * e1_1 * e2_1;
  MA[5] = ((e0_1 * e0_1 - e1_1 * e1_1) + e2_1 * e2_1) - e3_1 * e3_1;
  MA[9] = 2.0 * e2_1 * e3_1 - 2.0 * e0_1 * e1_1;
  MA[13] = i2->m_oldPos[1];
  MA[2] = 2.0 * e1_1 * e3_1 - 2.0 * e0_1 * e2_1;
  MA[6] = 2.0 * e0_1 * e1_1 + 2.0 * e2_1 * e3_1;
  MA[10] = ((e0_1 * e0_1 - e1_1 * e1_1) - e2_1 * e2_1) + e3_1 * e3_1;
  MA[14] = i2->m_oldPos[2];
  MA[3] = F0;
  MA[7] = F0;
  MA[11] = F0;
  MA[15] = F1;

  MFloat invMB[16];
  invMB[0] = ((e0_2 * e0_2 + e1_2 * e1_2) - e2_2 * e2_2) - e3_2 * e3_2;
  invMB[4] = 2.0 * e0_2 * e3_2 + 2.0 * e1_2 * e2_2;
  invMB[8] = 2.0 * e1_2 * e3_2 - 2.0 * e0_2 * e2_2;
  invMB[12] = -(invMB[0] * i3->m_oldPos[0] + invMB[4] * i3->m_oldPos[1] + invMB[8] * i3->m_oldPos[2]);
  invMB[1] = 2.0 * e1_2 * e2_2 - 2.0 * e0_2 * e3_2;
  invMB[5] = ((e0_2 * e0_2 - e1_2 * e1_2) + e2_2 * e2_2) - e3_2 * e3_2;
  invMB[9] = 2.0 * e0_2 * e1_2 + 2.0 * e2_2 * e3_2;
  invMB[13] = -(invMB[1] * i3->m_oldPos[0] + invMB[5] * i3->m_oldPos[1] + invMB[9] * i3->m_oldPos[2]);
  invMB[2] = 2.0 * e0_2 * e2_2 + 2.0 * e1_2 * e3_2;
  invMB[6] = 2.0 * e2_2 * e3_2 - 2.0 * e0_2 * e1_2;
  invMB[10] = ((e0_2 * e0_2 - e1_2 * e1_2) - e2_2 * e2_2) + e3_2 * e3_2;
  invMB[14] = -(invMB[2] * i3->m_oldPos[0] + invMB[6] * i3->m_oldPos[1] + invMB[10] * i3->m_oldPos[2]);
  invMB[3] = F0;
  invMB[7] = F0;
  invMB[11] = F0;
  invMB[15] = F1;

  MFloat v_1[3];
  v_1[0] = (i2->m_position[0] - i2->m_oldPos[0]); /// m_particleTimeStep;
  v_1[1] = (i2->m_position[1] - i2->m_oldPos[1]); /// m_particleTimeStep;
  v_1[2] = (i2->m_position[2] - i2->m_oldPos[2]); /// m_particleTimeStep;

  MFloat v_2[3];
  v_2[0] = (i3->m_position[0] - i3->m_oldPos[0]); /// m_particleTimeStep;
  v_2[1] = (i3->m_position[1] - i3->m_oldPos[1]); /// m_particleTimeStep;
  v_2[2] = (i3->m_position[2] - i3->m_oldPos[2]); /// m_particleTimeStep;

  MFloat A[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0};
  A[0] = 1 / (i2->m_semiMinorAxis * i2->m_semiMinorAxis);
  A[5] = A[0];
  A[10] = 1 / (i2->m_semiMinorAxis * i2->m_semiMinorAxis * i2->m_aspectRatio * i2->m_aspectRatio);

  MFloat detMinusB = -1
                     / (i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_semiMinorAxis
                        * i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_aspectRatio * i3->m_aspectRatio);

  MFloat m[5];
  MFloat b_Matrix[16];
  createMatrix_CCD(MA, v_1, invMB, v_2, i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio, b_Matrix, m);
  MInt iteration = 0;
  MFloat t1 = 0.0;
  MFloat t2 = 1.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  setup until here */

  MBool collision;
  MFloat positive1 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t1);
  if(positive1 < 0.0) {
    collision = true;
  } else {
    MFloat positive2 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t2);
    if(positive2 < 0.0) {
      collision = true;
    } else {
      collision = checkCSI_CCD(A, b_Matrix, m, detMinusB, t1, t2, positive1, positive2, iteration);
    }
  }
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  collision detection until here */
  if(collision) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i3).m_partId;
    thisColl.semiMinorAxis0 = (*i2).m_semiMinorAxis;
    thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
    thisColl.aspectRatio0 = (*i2).m_aspectRatio;
    thisColl.aspectRatio1 = (*i3).m_aspectRatio;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i3).m_densityRatio;
    thisColl.collPosX = i2->m_position[0];
    thisColl.collPosY = i2->m_position[1];
    thisColl.collPosZ = i2->m_position[2];
    collQueueEllipsoid.push(thisColl);
  }
}

/**  void LPT::collisionCheckSECCD(partListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3)
 *   \brief  checks continuously for the overlaps of the sphere i2 and the ellipsoid i3 assuming that they do not rotate
 * during two timesteps after Choi et al., IEEE TRANSACTIONS ON VISUALIZATION AND COMPUTER GRAPHICS, VOL. 15, NO. 2,
 * MARCH/APRIL 2009, 311 Continuous Collision Detection for Ellipsoids
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
void ParticleCollision<nDim>::collisionCheckSECCD(partListIteratorConst<nDim> i2, ellipsListIteratorConst<nDim> i3) {
  TRACE();

  const MInt dim3 = 3;

  MFloat deltaXSq = F0;
  MFloat deltaVSq = F0;
  MFloat deltaVdeltaX = F0;
  MFloat dist[dim3];
  for(MInt d = 0; d < dim3; ++d) {
    dist[d] = i2->m_oldPos[d] - i3->m_oldPos[d];
    MFloat deltaV = ((i2->m_position[d] - i2->m_oldPos[d]) - (i3->m_position[d] - i3->m_oldPos[d])) / m_timeStep;
    deltaXSq += POW2(dist[d]);
    deltaVSq += POW2(deltaV);
    deltaVdeltaX += dist[d] * deltaV;
  }
  if(approx(deltaVSq, F0, MFloatEps)) { // no relative velocity, check for overlap only
    collisionCheckSERetroActive(i2, i3);
    return;
  }

  MFloat i1minR = i2->m_diameter / F2;
  MFloat& i1maxR = i1minR;
  MFloat i2minR = min(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat i2maxR = max(i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio);
  MFloat RminSq = POW2(i1minR + i2minR);
  MFloat RmaxSq = POW2(i1maxR + i2maxR);

  // as a first rough criterion, check the distance between the centroids
  // if larger than the sum of major axes, then
  // for sure there is no collision between the two
  MFloat bigKmax = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RmaxSq);
  if(bigKmax < F0) {
    return; // only complex solutions, no collision on current route
  }
  bigKmax = sqrt(bigKmax);
  MFloat tempt1 = (-deltaVdeltaX - bigKmax) / deltaVSq;
  MFloat tempt2 = (-deltaVdeltaX + bigKmax) / deltaVSq;
  MFloat t1max = min(tempt1, tempt2);
  MFloat t2max = max(tempt1, tempt2);
  if(t1max > m_timeStep) {
    return; // collision begin in the future
  }
  if(t2max < F0) {
    return; // collision end in the past
  }

  // if smaller than the sum of the minor axis, then
  //  a collision between the two is certain
  MFloat bigKmin = POW2(deltaVdeltaX) - deltaVSq * (deltaXSq - RminSq);
  if(bigKmin >= F0) {
    bigKmin = sqrt(bigKmin);
    tempt1 = (-deltaVdeltaX - bigKmin) / deltaVSq;
    tempt2 = (-deltaVdeltaX + bigKmin) / deltaVSq;
    MFloat t1min = min(tempt1, tempt2);
    MFloat t2min = max(tempt1, tempt2);
    if((t1min <= m_timeStep) && (t2min >= F0)) {
      collQueueElemEllipsoid thisColl;
      thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + t1min;
      thisColl.part0 = (*i3).m_partId;
      thisColl.part1 = (*i2).m_partId;
      thisColl.semiMinorAxis0 = (*i3).m_semiMinorAxis;
      thisColl.semiMinorAxis1 = i1minR;
      thisColl.aspectRatio0 = (*i3).m_aspectRatio;
      thisColl.aspectRatio1 = -1.0;
      thisColl.dens0 = (*i3).m_densityRatio;
      thisColl.dens1 = (*i2).m_densityRatio;
      thisColl.collPosX = i2->m_oldPos[0] + (i2->m_position[0] - i2->m_oldPos[0]) * t1min / m_timeStep;
      thisColl.collPosY = i2->m_oldPos[1] + (i2->m_position[1] - i2->m_oldPos[1]) * t1min / m_timeStep;
      thisColl.collPosZ = i2->m_oldPos[2] + (i2->m_position[2] - i2->m_oldPos[2]) * t1min / m_timeStep;
      collQueueEllipsoid.push(thisColl);
      return;
    }
  }

  // distance is smaller than maximal distance but larger than minimal distance
  // so collision may occur


  //  MFloat e0_1 = 0.5;
  //  MFloat e1_1 = 0.5;
  //  MFloat e2_1 = 0.5;
  //  MFloat e3_1 = 0.5;

  MFloat e0_2 = i3->m_quaternion[0];
  MFloat e1_2 = i3->m_quaternion[1];
  MFloat e2_2 = i3->m_quaternion[2];
  MFloat e3_2 = i3->m_quaternion[3];

  MFloat MA[16];
  MA[0] = 1.0;
  MA[4] = 0.0;
  MA[8] = 0.0;
  MA[12] = i2->m_oldPos[0];
  MA[1] = 0.0;
  MA[5] = 1.0;
  MA[9] = 0.0;
  MA[13] = i2->m_oldPos[1];
  MA[2] = 0.0;
  MA[6] = 0.0;
  MA[10] = 1.0;
  MA[14] = i2->m_oldPos[2];
  MA[3] = F0;
  MA[7] = F0;
  MA[11] = F0;
  MA[15] = F1;

  MFloat invMB[16];
  invMB[0] = ((e0_2 * e0_2 + e1_2 * e1_2) - e2_2 * e2_2) - e3_2 * e3_2;
  invMB[4] = 2.0 * e0_2 * e3_2 + 2.0 * e1_2 * e2_2;
  invMB[8] = 2.0 * e1_2 * e3_2 - 2.0 * e0_2 * e2_2;
  invMB[12] = -(invMB[0] * i3->m_oldPos[0] + invMB[4] * i3->m_oldPos[1] + invMB[8] * i3->m_oldPos[2]);
  invMB[1] = 2.0 * e1_2 * e2_2 - 2.0 * e0_2 * e3_2;
  invMB[5] = ((e0_2 * e0_2 - e1_2 * e1_2) + e2_2 * e2_2) - e3_2 * e3_2;
  invMB[9] = 2.0 * e0_2 * e1_2 + 2.0 * e2_2 * e3_2;
  invMB[13] = -(invMB[1] * i3->m_oldPos[0] + invMB[5] * i3->m_oldPos[1] + invMB[9] * i3->m_oldPos[2]);
  invMB[2] = 2.0 * e0_2 * e2_2 + 2.0 * e1_2 * e3_2;
  invMB[6] = 2.0 * e2_2 * e3_2 - 2.0 * e0_2 * e1_2;
  invMB[10] = ((e0_2 * e0_2 - e1_2 * e1_2) - e2_2 * e2_2) + e3_2 * e3_2;
  invMB[14] = -(invMB[2] * i3->m_oldPos[0] + invMB[6] * i3->m_oldPos[1] + invMB[10] * i3->m_oldPos[2]);
  invMB[3] = F0;
  invMB[7] = F0;
  invMB[11] = F0;
  invMB[15] = F1;

  MFloat v_1[3];
  v_1[0] = (i2->m_position[0] - i2->m_oldPos[0]); /// m_particleTimeStep;
  v_1[1] = (i2->m_position[1] - i2->m_oldPos[1]); /// m_particleTimeStep;
  v_1[2] = (i2->m_position[2] - i2->m_oldPos[2]); /// m_particleTimeStep;

  MFloat v_2[3];
  v_2[0] = (i3->m_position[0] - i3->m_oldPos[0]); /// m_particleTimeStep;
  v_2[1] = (i3->m_position[1] - i3->m_oldPos[1]); /// m_particleTimeStep;
  v_2[2] = (i3->m_position[2] - i3->m_oldPos[2]); /// m_particleTimeStep;

  MFloat A[16] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0};
  A[0] = 4.0 / (i2->m_diameter * i2->m_diameter);
  A[5] = A[0];
  A[10] = A[0];

  MFloat detMinusB = -1
                     / (i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_semiMinorAxis
                        * i3->m_semiMinorAxis * i3->m_semiMinorAxis * i3->m_aspectRatio * i3->m_aspectRatio);

  MFloat m[5];
  MFloat b_Matrix[16];
  createMatrix_CCD(MA, v_1, invMB, v_2, i3->m_semiMinorAxis, i3->m_semiMinorAxis * i3->m_aspectRatio, b_Matrix, m);
  MInt iteration = 0;
  MFloat t1 = 0.0;
  MFloat t2 = 1.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  setup until here */

  MBool collision;
  MFloat positive1 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t1);
  if(positive1 < 0.0) {
    collision = true;
  } else {
    MFloat positive2 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t2);
    if(positive2 < 0.0) {
      collision = true;
    } else {
      collision = checkCSI_CCD(A, b_Matrix, m, detMinusB, t1, t2, positive1, positive2, iteration);
    }
  }

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /*  collision detection until here */
  if(collision) {
    collQueueElemEllipsoid thisColl;
    thisColl.collTimeStep = globalTimeStep; // solverPtr->m_time - m_particleTimeStep + tmin*m_particleTimeStep;
    thisColl.part0 = (*i2).m_partId;
    thisColl.part1 = (*i3).m_partId;
    thisColl.semiMinorAxis1 = (*i3).m_semiMinorAxis;
    thisColl.semiMinorAxis0 = i2->m_diameter / 2.0;
    thisColl.aspectRatio1 = (*i3).m_aspectRatio;
    thisColl.aspectRatio0 = -1.0;
    thisColl.dens0 = (*i2).m_densityRatio;
    thisColl.dens1 = (*i3).m_densityRatio;
    thisColl.collPosX = i2->m_oldPos[0];
    thisColl.collPosY = i2->m_oldPos[1];
    thisColl.collPosZ = i2->m_oldPos[2];
    collQueueEllipsoid.push(thisColl);
  }
}

/** MFloat LPT::checkBndryCross(MInt    bndryId,
 *                                          MInt    vectorId,
 *                                          MFloat particleDiameter)
 *  \brief Checks crossing of boundary surface of the boundary cell
 *
 *  The current cell in which the particle resides, is a boundary cell.
 *  It is checked whether the particle crossed the boundary surface of that
 *  cell in the last time step using the new and old particle positions.
 *
 *  @author Martin Brinks
 */
template <MInt nDim>
MFloat ParticleCollision<nDim>::checkBndryCross(MInt bndryId, partListIteratorConst<nDim> vectorId,
                                                MFloat particleDiameter) {
  TRACE();

  MFloat nx0 = 0.0, na = 0.0, nb = 0.0;

  for(MInt i = 0; i < nDim; i++) {
    nx0 += m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i]
           * m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_planeCoordinates[0][i];
    na += m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i] * (*vectorId).m_oldPos[i];
    nb +=
        m_lpt->m_bndryCells->a[bndryId].m_srfcs[0]->m_normal[i] * ((*vectorId).m_position[i] - (*vectorId).m_oldPos[i]);
  }

  return (nx0 - na + (0.5 * particleDiameter)) / nb;
}


/** void LPT::writeCollEvent(MFloat collisionTime,
 *                                    MInt   particle0,
 *                                    MInt   particle1,
 *                                    MFloat diam0,
 *                                    MFloat diam1,
 *                                    MFloat dens0,
 *                                    MFloat dens1,
 *                                    MFloat relVel,
 *                                    MFloat meanVel,
 *                                    MFloat *x)
 *  \brief  Add collision data to queue which is to be written
 *
 * store all data concerning particle-particle collisions:
 * collisionTime = exact time of impact
 * particle0     = index of particle 0
 * particle1     = index of particle 1
 * relVel        = relative velocity  |w| = |v_0 - v_1|
 * meanVel       = mean velocity      |W| = 0.5 * |v_0 + v_1|
 * x[nDim]   = point of impact
 *
 *  @author Rudie Kunnen, modified by Jerry Grimmen (29-4-2010)
 */
template <MInt nDim>
void ParticleCollision<nDim>::writeCollEvent(MFloat collisionTime, MInt particle0, MInt particle1, MFloat diam0,
                                             MFloat diam1, MFloat dens0, MFloat dens1, MFloat relVel, MFloat meanVel,
                                             MFloat* x) {
  TRACE();

  // Create a collision element for the collision queue.
  // Needed to write collision datas into a collective Netcdf file.
  // m_log << std::setprecision(12) << "collTime " << collisionTime << " partId1 " << particle0 << " partId2 " <<
  //        particle1 << " timestep " << m_timeStep
  //        << " [" << m_time << ", " << m_time + m_timeStep << "]" << endl;
  collQueueElem collData = {collisionTime, diam0, diam1, dens0, dens1,     relVel,
                            meanVel,       x[0],  x[1],  x[2],  particle0, particle1};
  collQueue.push(collData);
}

/**  void LPT::createMatrix_CCD(const MFloat (&MA)[16], const MFloat (&u1)[3], const MFloat (&invMB)[16], const
 * MFloat (&u2)[3], const MFloat& a_2, const MFloat& c_2, MFloat (&n)[16], MFloat (&m)[5]) \brief  calculates
 * M_A^T*(M_B^(-1))^T*B*(M_B^(-1))*M_A, write it the part independent of t in n and the other values in m
 *
 *           M_A is the matrix of Ellipsoid2
 *           u1 are the time-dependent values of M_A
 *           invM_B is the matrix of Ellipsoid3
 *           u2 are the time-dependent values of M_B
 *           a_2 is the semiAxis of Ellipsoid3
 *           c_2 is the semiAxis*aspectRatio of Ellipsoid3
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
void ParticleCollision<nDim>::createMatrix_CCD(const MFloat (&MA)[16], const MFloat (&u1)[3], const MFloat (&invMB)[16],
                                               const MFloat (&u2)[3], const MFloat& a_2, const MFloat& c_2,
                                               MFloat (&n)[16], MFloat (&m)[5]) {
  TRACE();

  MFloat b_j1;
  MFloat k1, l1;
  MFloat d1, d2, d3;
  MFloat bla1, bla2, bla3;
  MFloat y;
  MFloat b_y, c_y, d_y, e_y, f_y;
  b_j1 = (-invMB[0] * u2[0] - invMB[4] * u2[1]) - invMB[8] * u2[2];
  k1 = (-invMB[1] * u2[0] - invMB[5] * u2[1]) - invMB[9] * u2[2];
  l1 = (-invMB[2] * u2[0] - invMB[6] * u2[1]) - invMB[10] * u2[2];
  d1 = (MA[0] * invMB[2] + MA[1] * invMB[6]) + MA[2] * invMB[10];
  d2 = (MA[0] * invMB[1] + MA[1] * invMB[5]) + MA[2] * invMB[9];
  d3 = (MA[0] * invMB[0] + MA[1] * invMB[4]) + MA[2] * invMB[8];
  bla1 = (invMB[8] * d3 / (a_2 * a_2) + invMB[9] * d2 / (a_2 * a_2)) + invMB[10] * d1 / (c_2 * c_2);
  bla2 = (invMB[4] * d3 / (a_2 * a_2) + invMB[5] * d2 / (a_2 * a_2)) + invMB[6] * d1 / (c_2 * c_2);
  bla3 = (invMB[0] * d3 / (a_2 * a_2) + invMB[1] * d2 / (a_2 * a_2)) + invMB[2] * d1 / (c_2 * c_2);
  n[0] = (MA[2] * bla1 + MA[1] * bla2) + MA[0] * bla3;
  n[4] = (MA[6] * bla1 + MA[5] * bla2) + MA[4] * bla3;
  n[8] = (MA[10] * bla1 + MA[9] * bla2) + MA[8] * bla3;

  d1 = (MA[4] * invMB[2] + MA[5] * invMB[6]) + MA[6] * invMB[10];
  d2 = (MA[4] * invMB[1] + MA[5] * invMB[5]) + MA[6] * invMB[9];
  d3 = (MA[4] * invMB[0] + MA[5] * invMB[4]) + MA[6] * invMB[8];
  bla1 = (invMB[8] * d3 / (a_2 * a_2) + invMB[9] * d2 / (a_2 * a_2)) + invMB[10] * d1 / (c_2 * c_2);
  bla2 = (invMB[4] * d3 / (a_2 * a_2) + invMB[5] * d2 / (a_2 * a_2)) + invMB[6] * d1 / (c_2 * c_2);
  bla3 = (invMB[0] * d3 / (a_2 * a_2) + invMB[1] * d2 / (a_2 * a_2)) + invMB[2] * d1 / (c_2 * c_2);
  n[1] = (MA[2] * bla1 + MA[1] * bla2) + MA[0] * bla3;
  n[5] = (MA[6] * bla1 + MA[5] * bla2) + MA[4] * bla3;
  n[9] = (MA[10] * bla1 + MA[9] * bla2) + MA[8] * bla3;

  d1 = (MA[8] * invMB[2] + MA[9] * invMB[6]) + MA[10] * invMB[10];
  d2 = (MA[8] * invMB[1] + MA[9] * invMB[5]) + MA[10] * invMB[9];
  d3 = (MA[8] * invMB[0] + MA[9] * invMB[4]) + MA[10] * invMB[8];
  bla1 = (invMB[8] * d3 / (a_2 * a_2) + invMB[9] * d2 / (a_2 * a_2)) + invMB[10] * d1 / (c_2 * c_2);
  bla2 = (invMB[4] * d3 / (a_2 * a_2) + invMB[5] * d2 / (a_2 * a_2)) + invMB[6] * d1 / (c_2 * c_2);
  bla3 = (invMB[0] * d3 / (a_2 * a_2) + invMB[1] * d2 / (a_2 * a_2)) + invMB[2] * d1 / (c_2 * c_2);
  n[2] = (MA[2] * bla1 + MA[1] * bla2) + MA[0] * bla3;
  n[6] = (MA[6] * bla1 + MA[5] * bla2) + MA[4] * bla3;
  n[10] = (MA[10] * bla1 + MA[9] * bla2) + MA[8] * bla3;

  d1 = ((invMB[14] + invMB[2] * MA[12]) + invMB[6] * MA[13]) + invMB[10] * MA[14];
  d2 = ((invMB[13] + invMB[1] * MA[12]) + invMB[5] * MA[13]) + invMB[9] * MA[14];
  d3 = ((invMB[12] + invMB[0] * MA[12]) + invMB[4] * MA[13]) + invMB[8] * MA[14];
  bla1 = (invMB[8] * d3 / (a_2 * a_2) + invMB[9] * d2 / (a_2 * a_2)) + invMB[10] * d1 / (c_2 * c_2);
  bla2 = (invMB[4] * d3 / (a_2 * a_2) + invMB[5] * d2 / (a_2 * a_2)) + invMB[6] * d1 / (c_2 * c_2);
  bla3 = (invMB[0] * d3 / (a_2 * a_2) + invMB[1] * d2 / (a_2 * a_2)) + invMB[2] * d1 / (c_2 * c_2);
  n[3] = (MA[2] * bla1 + MA[1] * bla2) + MA[0] * bla3;
  n[7] = (MA[6] * bla1 + MA[5] * bla2) + MA[4] * bla3;
  n[11] = (MA[10] * bla1 + MA[9] * bla2) + MA[8] * bla3;
  n[15] = (((((-1.0 + MA[12] * bla3) + MA[13] * bla2) + MA[14] * bla1) + invMB[12] * d3 / (a_2 * a_2))
           + invMB[13] * d2 / (a_2 * a_2))
          + invMB[14] * d1 / (c_2 * c_2);
  y = u1[0] * bla3;
  b_y = u1[1] * bla2;
  c_y = u1[2] * bla1;
  d_y = b_j1 * d3 / (a_2 * a_2);
  e_y = k1 * d2 / (a_2 * a_2);
  f_y = l1 * d1 / (c_2 * c_2);
  n[12] = n[3];
  n[13] = n[7];
  n[14] = n[11];
  d1 = ((l1 + invMB[2] * u1[0]) + invMB[6] * u1[1]) + invMB[10] * u1[2];
  d2 = ((k1 + invMB[1] * u1[0]) + invMB[5] * u1[1]) + invMB[9] * u1[2];
  d3 = ((b_j1 + invMB[0] * u1[0]) + invMB[4] * u1[1]) + invMB[8] * u1[2];
  bla1 = (invMB[8] * d3 / (a_2 * a_2) + invMB[9] * d2 / (a_2 * a_2)) + invMB[10] * d1 / (c_2 * c_2);
  bla2 = (invMB[4] * d3 / (a_2 * a_2) + invMB[5] * d2 / (a_2 * a_2)) + invMB[6] * d1 / (c_2 * c_2);
  bla3 = (invMB[0] * d3 / (a_2 * a_2) + invMB[1] * d2 / (a_2 * a_2)) + invMB[2] * d1 / (c_2 * c_2);
  m[0] = (MA[2] * bla1 + MA[1] * bla2) + MA[0] * bla3;
  m[1] = (MA[6] * bla1 + MA[5] * bla2) + MA[4] * bla3;
  m[2] = (MA[10] * bla1 + MA[9] * bla2) + MA[8] * bla3;
  m[3] = ((((((((((y + b_y) + c_y) + d_y) + e_y) + f_y) + MA[12] * bla3) + MA[13] * bla2) + MA[14] * bla1)
           + invMB[12] * d3 / (a_2 * a_2))
          + invMB[13] * d2 / (a_2 * a_2))
         + invMB[14] * d1 / (c_2 * c_2);
  m[4] = ((((u1[0] * bla3 + u1[1] * bla2) + u1[2] * bla1) + b_j1 * d3 / (a_2 * a_2)) + k1 * d2 / (a_2 * a_2))
         + l1 * d1 / (c_2 * c_2);
}

/**  MFloat LPT::checkStaticOverlap_CCD(const MFloat (&b_Matrix)[16], const MFloat (&A)[16], const MFloat
 * (&m)[5], const MFloat& detMinusB, MFloat& t) \brief  calculates if there is a local maximum of det(A*(u-1)-u*B)
 * in u for a fixed value of t and if the value there is above zero since this necessary for separation returns this
 * value otherwise -1. b_Matrix is the build matrix by createMatrix_CCD independent of t m are the values build by
 * createMatrix_CCD dependent of t A is the matrix of Ellpsoid2 detMinusB is det(-B) where B is the matrix of Ellipsoid3
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
MFloat ParticleCollision<nDim>::checkStaticOverlap_CCD(const MFloat (&b_Matrix)[16], const MFloat (&A)[16],
                                                       const MFloat (&m)[5], const MFloat& detMinusB, MFloat& t) {
  TRACE();


  const MFloat eps = 1.0E-13;

  MFloat b[16];
  copy(b_Matrix, b_Matrix + 16, b);
  b[12] = b[12] + m[0] * t;
  b[13] = b[13] + m[1] * t;
  b[14] = b[14] + m[2] * t;
  b[15] = b[15] + m[3] * t + m[4] * t * t;

  MFloat b12s = b[4] * b[4];
  MFloat b13s = b[8] * b[8];
  MFloat b14s = b[12] * b[12];
  MFloat b23s = b[9] * b[9];
  MFloat b24s = b[13] * b[13];
  MFloat b34s = b[14] * b[14];
  MFloat b2233 = b[5] * b[10];
  MFloat termA = (b[0] * A[5] * A[10] + b[5] * A[0] * A[10]) + b[10] * A[0] * A[5];
  MFloat termB = ((b2233 - b23s) * A[0] + (b[0] * b[10] - b13s) * A[5]) + (b[0] * b[5] - b12s) * A[10];

  MFloat T4 = -A[0] * A[5] * A[10];
  MFloat T3 = termA + b[15] * T4;
  MFloat T2 = (((termA * b[15] - termB) - b34s * A[0] * A[5]) - b14s * A[5] * A[10]) - b24s * A[0] * A[10];
  MFloat tmp1 = termB * b[15];
  MFloat tmp2 = b[0] * (b2233 + A[5] * b34s + A[10] * b24s - b23s);
  MFloat tmp3 = b[5] * (A[0] * b34s + A[10] * b14s - b13s);
  MFloat tmp4 = b[10] * (A[0] * b24s + A[5] * b14s - b12s);
  MFloat tmp5 = b[14] * (A[0] * b[9] * b[13] + A[5] * b[8] * b[12]) + b[4] * (A[10] * b[12] * b[13] - b[8] * b[9]);
  tmp5 = tmp5 + tmp5;
  MFloat T1 = -tmp1 + tmp2 + tmp3 + tmp4 - tmp5;
  MFloat T0 = detMinusB;

  MFloat t4 = (T0 + T1 + T2 + T3 + T4);
  MFloat t3 = (-T1 - 2 * T2 - 3 * T3 - 4 * T4);
  MFloat t2 = (T2 + 3 * T3 + 6 * T4);
  MFloat t1 = (-T3 - 4 * T4);
  MFloat t0 = T4;

  MFloat a = 4 * t4;
  MFloat bb = 3 * t3;
  MFloat c = 2 * t2;
  MFloat d = 1 * t1;
  if(fabs(a) > eps) {
    MFloat p = (3 * a * c - bb * bb) / (3 * a * a);
    MFloat q = (2 * bb * bb * bb - 9 * a * bb * c + 27 * a * a * d) / (27 * a * a * a);
    if(fabs(p) > eps) {
      if(p > 0) {
        MFloat term = sqrt(p / 3.0);
        MFloat u = -2 * term * sinh((1.0 / 3.0) * asinh((3 * q / (2 * p)) / term)) - bb / (3 * a);
        return testU_CCD(u, t4, t3, t2, t1, t0);
      } else {
        MFloat test = 4 * p * p * p + 27 * q * q;
        MFloat term = sqrt(-p / 3.0);
        if(test > 0) {
          MFloat u =
              -2 * fabs(q) / q * term * cosh((1.0 / 3.0) * acosh((-3 * fabs(q) / (2 * p)) / term)) - bb / (3 * a);
          return testU_CCD(u, t4, t3, t2, t1, t0);
        } else {
          for(MInt that = 0; that < 3; that++) {
            MFloat u =
                2 * term * cos((1.0 / 3.0) * acos((3 * q / (2 * p)) / term) - that * 2 * PI / 3.0) - bb / (3 * a);
            MFloat positive = testU_CCD(u, t4, t3, t2, t1, t0);
            if(positive > -1) {
              return positive;
            }
          }
        }
      }
    } else {
      if(fabs(q) > eps) {
        MFloat u = pow(-q, (1.0 / 3.0)) - bb / (3 * a);
        return testU_CCD(u, t4, t3, t2, t1, t0);
      } else {
        MFloat u = -bb / (3 * a);
        return testU_CCD(u, t4, t3, t2, t1, t0);
      }
    }
  } else {
    if(fabs(bb) > eps) {
      MFloat p_2 = c / bb / 2.0;
      MFloat q = d / bb;
      MFloat temp = p_2 * p_2 - q;
      if(temp >= 0) {
        temp = sqrt(temp);
        MFloat u = p_2 + temp;
        MFloat positive = testU_CCD(u, t4, t3, t2, t1, t0);
        if(positive > -1) {
          return positive;
        }
        u = p_2 - temp;
        return testU_CCD(u, t4, t3, t2, t1, t0);
      } else {
        return -F1;
      }
    } else {
      if(fabs(c) > eps) {
        MFloat u = d / c;
        return testU_CCD(u, t4, t3, t2, t1, t0);
      } else {
        return -F1;
      }
    }
  }

  return -F1;
}

/**  MBool LPT::checkCSI_CCD(const MFloat (&A)[16], const MFloat (&b_Matrix)[16], const MFloat (&m)[5], const
 * MFloat& detMinusB, MFloat& t1, MFloat& t2, MFloat& positive1, MFloat& positive2, MInt& iteration) \brief
 * checks recursively if the ellipsoids are separated by subdividing the time intervals in smaller CSIs. returns false
 * if separated, otherwise returns true
 *
 *           A is the matrix of Ellpsoid2
 *           b_Matrix is the build matrix by createMatrix_CCD independent of t
 *           m are the values build by createMatrix_CCD dependent of t
 *           detMinusB is det(-B) where B is the matrix of Ellipsoid3
 *           t1 is the beginning time
 *           t2 is the ending time
 *           positive1 is the u value of the local maximum for t1
 *           positive2 is the u value of the local maximum for t2
 *           iteration is the number of recursive calls
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
MBool ParticleCollision<nDim>::checkCSI_CCD(const MFloat (&A)[16], const MFloat (&b_Matrix)[16], const MFloat (&m)[5],
                                            const MFloat& detMinusB, MFloat& t1, MFloat& t2, MFloat& positive1,
                                            MFloat& positive2, MInt& iteration) {
  TRACE();

  if(t2 - t1 < 1.0E-13) {
    return false;
  } else {
    iteration++;

    MFloat b_A[16];
    for(MInt q0 = 0; q0 < 16; q0++) {
      b_A[q0] = A[q0] * (positive1 - 1.0) - positive1 * b_Matrix[q0];
    }

    MFloat b_positive1[5];
    for(MInt q0 = 0; q0 < 5; q0++) {
      b_positive1[q0] = -positive1 * m[q0];
    }

    t1 = checkDynamicOverlap_CCD(b_A, b_positive1, t1, t2);
    if(t1 >= t2) {
      return false;
    } else {
      for(MInt q0 = 0; q0 < 16; q0++) {
        b_A[q0] = A[q0] * (positive2 - 1.0) - positive2 * b_Matrix[q0];
      }

      for(MInt q0 = 0; q0 < 5; q0++) {
        b_positive1[q0] = -positive2 * m[q0];
      }

      t2 = checkDynamicOverlap_CCD(b_A, b_positive1, t2, t1);
      if(t1 >= t2) {
        return false;
      } else {
        MFloat t_neu = (t1 + t2) / 2.0;
        MFloat positive_neu = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t_neu);
        if(positive_neu < 0.0) {
          return true;
        } else {
          positive1 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t1);
          MBool value = checkCSI_CCD(A, b_Matrix, m, detMinusB, t1, t_neu, positive1, positive_neu, iteration);
          if(value) {
            return true;
          } else {
            positive2 = checkStaticOverlap_CCD(b_Matrix, A, m, detMinusB, t2);
            if(positive2 < 0.0) {
              return true;
            } else {
              return checkCSI_CCD(A, b_Matrix, m, detMinusB, t_neu, t2, positive_neu, positive2, iteration);
            }
          }
        }
      }
    }
  }
}

/**  MFloat LPT::testU_CCD(const MFloat& u ,const MFloat& t4, const MFloat& t3, const MFloat& t2, const
 * MFloat& t1, const MFloat& t0) \brief  calculates the value of the forth order polynomial and returns it if
 * positive other returns -1
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
MFloat ParticleCollision<nDim>::testU_CCD(const MFloat& u, const MFloat& t4, const MFloat& t3, const MFloat& t2,
                                          const MFloat& t1, const MFloat& t0) {
  TRACE();

  MFloat positive = -1;
  if((u > 0) && (u <= 1)) {
    MFloat result = t4 * u * u * u * u + t3 * u * u * u + t2 * u * u + t1 * u + t0;
    if(result >= 0) {
      positive = u;
    }
  }

  return positive;
}

/**  MFloat LPT::checkDynamicOverlap_CCD(const MFloat (&n)[16], const MFloat (&m)[5], const MInt& side)
 *   \brief  Performes BezierShot from "side", returns the first zero-crossing
 *           results in finding the roots of the determinate of the Matrix for fixed "u" value, which is a quadratic
 *           equation in t.
 *           returns the t value of the root if one exists between 0 and 1 otherwise returns side.
 *           "n" is the symmetric matrix containing the values independent of t, "m" are the coefficients t is
 * multiplied to "side" is -1 or 1 depending on the direction of the BezierShot
 *
 *
 *   Aug-2013
 *   @author Christoph Siewert
 */
template <MInt nDim>
MFloat ParticleCollision<nDim>::checkDynamicOverlap_CCD(const MFloat (&n)[16], const MFloat (&m)[5], MFloat t1,
                                                        MFloat t2) {
  TRACE();

  MInt side = 1;
  if(t2 < t1) {
    side = -1;
    MInt temp = (MInt)t2;
    t2 = t1;
    t1 = temp;
  }

  MFloat t;
  MFloat AA;
  MFloat BB;
  MFloat CC;
  MFloat p_2;
  AA = (((((((((((((((m[0] * m[0] * (n[9] * n[9]) - n[5] * n[10] * (m[0] * m[0])) + 2.0 * n[10] * m[0] * m[1] * n[4])
                    - 2.0 * m[0] * m[1] * n[8] * n[9])
                   - 2.0 * m[0] * m[2] * n[4] * n[9])
                  + 2.0 * n[5] * m[0] * m[2] * n[8])
                 + m[1] * m[1] * (n[8] * n[8]))
                - n[0] * n[10] * (m[1] * m[1]))
               - 2.0 * m[1] * m[2] * n[4] * n[8])
              + 2.0 * n[0] * m[1] * m[2] * n[9])
             + m[2] * m[2] * (n[4] * n[4]))
            - n[0] * n[5] * (m[2] * m[2]))
           - m[4] * n[10] * (n[4] * n[4]))
          + 2.0 * m[4] * n[4] * n[8] * n[9])
         - m[4] * n[5] * (n[8] * n[8]))
        - m[4] * n[0] * (n[9] * n[9]))
       + m[4] * n[0] * n[5] * n[10];
  BB = (((((((((((((((((((((2.0 * n[13] * m[1] * (n[8] * n[8]) + 2.0 * n[14] * m[2] * (n[4] * n[4]))
                           + 2.0 * n[12] * m[0] * (n[9] * n[9]))
                          - m[3] * n[0] * (n[9] * n[9]))
                         - m[3] * (n[8] * n[8]) * n[5])
                        - m[3] * (n[4] * n[4]) * n[10])
                       - 2.0 * n[13] * m[2] * n[4] * n[8])
                      - 2.0 * m[1] * n[14] * n[4] * n[8])
                     - 2.0 * n[12] * m[1] * n[8] * n[9])
                    - 2.0 * m[0] * n[13] * n[8] * n[9])
                   - 2.0 * n[12] * m[2] * n[4] * n[9])
                  + 2.0 * n[12] * m[2] * n[8] * n[5])
                 - 2.0 * m[0] * n[14] * n[4] * n[9])
                + 2.0 * m[0] * n[14] * n[8] * n[5])
               + 2.0 * n[13] * m[2] * n[0] * n[9])
              + 2.0 * m[1] * n[14] * n[0] * n[9])
             - 2.0 * n[14] * m[2] * n[0] * n[5])
            + 2.0 * n[12] * m[1] * n[4] * n[10])
           + 2.0 * m[0] * n[13] * n[4] * n[10])
          - 2.0 * n[13] * m[1] * n[0] * n[10])
         - 2.0 * n[12] * m[0] * n[5] * n[10])
        + 2.0 * m[3] * n[4] * n[8] * n[9])
       + m[3] * n[0] * n[5] * n[10];
  CC = (((((((((((((((n[12] * n[12] * (n[9] * n[9]) - n[5] * n[10] * (n[12] * n[12]))
                     + 2.0 * n[10] * n[12] * n[13] * n[4])
                    - 2.0 * n[12] * n[13] * n[8] * n[9])
                   - 2.0 * n[12] * n[14] * n[4] * n[9])
                  + 2.0 * n[5] * n[12] * n[14] * n[8])
                 + n[13] * n[13] * (n[8] * n[8]))
                - n[0] * n[10] * (n[13] * n[13]))
               - 2.0 * n[13] * n[14] * n[4] * n[8])
              + 2.0 * n[0] * n[13] * n[14] * n[9])
             + n[14] * n[14] * (n[4] * n[4]))
            - n[0] * n[5] * (n[14] * n[14]))
           - n[15] * n[10] * (n[4] * n[4]))
          + 2.0 * n[15] * n[4] * n[8] * n[9])
         - n[15] * n[5] * (n[8] * n[8]))
        - n[15] * n[0] * (n[9] * n[9]))
       + n[15] * n[0] * n[5] * n[10];
  if(fabs(AA) > 1.0E-13) {
    p_2 = BB / (2.0 * AA);
    BB = p_2 * p_2 - CC / AA;
    if(BB > 0.0) {
      BB = side * sqrt(BB);
      t = -p_2 - BB;
      if((t >= t1) && (t <= t2)) {
      } else {
        t = -p_2 + BB;
        if((t >= t1) && (t <= t2)) {
        } else {
          t = side;
        }
      }
    } else {
      t = side;
    }
  } else if(fabs(BB) > 1.0E-13) {
    t = CC / BB;
    if((t >= t1) && (t <= t2)) {
    } else {
      t = side;
    }
  } else {
    t = side;
  }

  return t;
}

/**
 *  \brief write collision data queue to Netcdf file (using parallel output)
 *
 *  @author Jerry Grimmen, Apr-10
 */
template <MInt nDim>
void ParticleCollision<nDim>::writeCollData() {
  TRACE();

  // Read the size of the collision queue and put it into the Netcdf file.
  MInt ncmpiCollQueueSize = (MInt)collQueue.size();
  // Read the amount of collisions of each domain.
  MIntScratchSpace ncmpiCollCount(noDomains(), AT_, "ncmpiCollCount");
  MPI_Allgather(&ncmpiCollQueueSize, 1, MPI_INT, &ncmpiCollCount[0], 1, MPI_INT, mpiComm(), AT_, "ncmpiCollQueueSize",
                "ncmpiCollCount[0]");

  MInt ncmpiCollCountMax = 0;
  for(MInt i = 0; i < noDomains(); ++i) {
    ncmpiCollCountMax += ncmpiCollCount[i];
  }

  // Only write collision data if there is any collision at all.
  if(ncmpiCollCountMax != 0) {
    // Change file name for collision data and create a new file. (Leaves file
    // open and in define mode).
    stringstream ncmpistream;
    string ncmpiFileName;
    ncmpistream << "out/collData." << globalTimeStep << ParallelIo::fileExt();
    ncmpistream >> ncmpiFileName;
    ncmpistream.clear();

    using namespace maia::parallel_io;
    ParallelIo parallelIo(ncmpiFileName, PIO_REPLACE, mpiComm());

    // Define Attribute timestep (double).
    parallelIo.setAttribute(m_lpt->m_timeStep, "particleTimestep");

    // Define Dimension collCount [noDomains] & Define Variable collCount
    // (int) [collCount].
    parallelIo.defineArray(PIO_INT, "collCount", noDomains());

    // Define dimensions and variables for all collision data.
    parallelIo.defineArray(PIO_FLOAT, "collTime", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_INT, "part0", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_INT, "part1", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "p0diam", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "p1diam", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "p0dens", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "p1dens", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "relVel", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "meanVel", ncmpiCollCountMax);
    parallelIo.defineArray(PIO_FLOAT, "collPos", 3 * ncmpiCollCountMax);

    // EndDef (Go into data mode).

    ParallelIo::size_type ncmpiStart = domainId();
    ParallelIo::size_type ncmpiCount = 1;

    parallelIo.setOffset(ncmpiCount, ncmpiStart);
    parallelIo.writeArray(&ncmpiCollQueueSize, "collCount");

    ncmpiStart = 0;
    for(MInt i = 0; i < domainId(); ++i) {
      ncmpiStart += ncmpiCollCount[i];
    }
    ncmpiCount = ncmpiCollCount[domainId()];
    if(ncmpiStart >= ncmpiCollCountMax) {
      if(ncmpiCount == 0) {
        ncmpiStart = 0;
      } else {
        mTerm(1, AT_,
              "Error in LPT::writeCollData!! ncmpiStart >= "
              "ncmpiCollCountMax but ncmpiCount != 0");
      }
    }

    // Arrays to hold collisions data.
    MFloatScratchSpace ncmpiCollTime(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MIntScratchSpace ncmpiPart0(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MIntScratchSpace ncmpiPart1(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiP0Diam(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiP1Diam(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiP0Dens(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiP1Dens(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiRelVel(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiMeanVel(ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MFloatScratchSpace ncmpiCollPos(3 * ncmpiCollCount[domainId()], AT_, "ncmpiCollTime");
    MInt ncmpiCountId = 0;

    // Write the data from the queue into the array.
    while(!collQueue.empty()) {
      ncmpiCollTime[ncmpiCountId] = collQueue.front().collTime;
      ncmpiPart0[ncmpiCountId] = collQueue.front().part0;
      ncmpiPart1[ncmpiCountId] = collQueue.front().part1;
      ncmpiP0Diam[ncmpiCountId] = collQueue.front().diam0;
      ncmpiP1Diam[ncmpiCountId] = collQueue.front().diam1;
      ncmpiP0Dens[ncmpiCountId] = collQueue.front().dens0;
      ncmpiP1Dens[ncmpiCountId] = collQueue.front().dens1;
      ncmpiRelVel[ncmpiCountId] = collQueue.front().relVel;
      ncmpiMeanVel[ncmpiCountId] = collQueue.front().meanVel;
      ncmpiCollPos[(3 * ncmpiCountId)] = collQueue.front().collPosX;
      ncmpiCollPos[(3 * ncmpiCountId) + 1] = collQueue.front().collPosY;
      ncmpiCollPos[(3 * ncmpiCountId) + 2] = collQueue.front().collPosZ;

      collQueue.pop();
      ++ncmpiCountId;
    }

    // Put the arrays into the Netcdf file.
    parallelIo.setOffset(ncmpiCount, ncmpiStart);
    parallelIo.writeArray(&ncmpiCollTime[0], "collTime");
    parallelIo.writeArray(&ncmpiPart0[0], "part0");
    parallelIo.writeArray(&ncmpiPart1[0], "part1");
    parallelIo.writeArray(&ncmpiP0Diam[0], "p0diam");
    parallelIo.writeArray(&ncmpiP1Diam[0], "p1diam");
    parallelIo.writeArray(&ncmpiP0Dens[0], "p0dens");
    parallelIo.writeArray(&ncmpiP1Dens[0], "p1dens");
    parallelIo.writeArray(&ncmpiRelVel[0], "relVel");
    parallelIo.writeArray(&ncmpiMeanVel[0], "meanVel");

    ncmpiCount *= 3;
    ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;
    parallelIo.setOffset(ncmpiCount, ncmpi3Start);
    parallelIo.writeArray(&ncmpiCollPos[0], "collPos");

    // Free memory.
    while(!collQueue.empty()) {
      collQueue.pop();
    }
  }

  if(m_includeEllipsoids) {
    // Read the size of the collision queue and put it into the Netcdf file.
    ncmpiCollQueueSize = (MInt)collQueueEllipsoid.size();
    // Read the amount of collisions of each domain.
    MIntScratchSpace ncmpiEllipsoidCollCount(noDomains(), AT_, "ncmpiEllipsoidCollCount");
    MPI_Allgather(&ncmpiCollQueueSize, 1, MPI_INT, &ncmpiEllipsoidCollCount[0], 1, MPI_INT, mpiComm(), AT_,
                  "ncmpiCollQueueSize", "ncmpiEllipsoidCollCount[0]");

    ncmpiCollCountMax = 0;
    for(MInt i = 0; i < noDomains(); ++i) {
      ncmpiCollCountMax += ncmpiEllipsoidCollCount[i];
    }

    // Only write collision data if there is any collision at all.
    if(ncmpiCollCountMax != 0) {
      // Change file name for collision data and create a new file. (Leaves
      // file open and in define mode).
      stringstream ncmpistream;
      string ncmpiFileName;
      ncmpistream << "out/collEllipsoid." << globalTimeStep << ParallelIo::fileExt();
      ncmpistream >> ncmpiFileName;
      ncmpistream.clear();

      // WARNING: untested switch from NetCDF/Parallel netCDF to ParallelIo
      // The method previously used direct I/O calls, which were replaced by
      // ParallelIo methods in summer 2015. However, since the method was not
      // used by any of the testcases, this code is still *untested*. Thus,
      // if your code uses this part of the code, please make sure that the
      // I/O still works as expected and then remove this warning as well as
      // the subsequent TERMM().
      TERMM(1, "untested I/O method, please see comment for how to proceed");
      using namespace maia::parallel_io;
      ParallelIo parallelIo(ncmpiFileName, PIO_REPLACE, mpiComm());
      ncmpiFileName.clear();

      // Define Attribute timestep (double).
      MFloat ncmpiTimeStep = globalTimeStep;
      parallelIo.setAttribute(ncmpiTimeStep, "timestep");

      // Define Dimension overlapCount [noDomains] & Define Variable collCount
      // (int) [collCount].
      parallelIo.defineArray(PIO_INT, "overlapCount", noDomains());

      // Define dimensions and variables for all collision datas.
      parallelIo.defineArray(PIO_FLOAT, "collTimeStep", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_INT, "part0", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_INT, "part1", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p0semiMinorAxis", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p1semiMinorAxis", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p0aspectRatio", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p1aspectRatio", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p0dens", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "p1dens", ncmpiCollCountMax);
      parallelIo.defineArray(PIO_FLOAT, "collPos", 3 * ncmpiCollCountMax);

      // EndDef (Go into data mode).

      ParallelIo::size_type ncmpiStart = domainId();
      ParallelIo::size_type ncmpiCount = 1;

      parallelIo.setOffset(ncmpiCount, ncmpiStart);
      parallelIo.writeArray(&ncmpiCollQueueSize, "overlapCount");

      ncmpiStart = 0;
      for(MInt i = 0; i < domainId(); ++i) {
        ncmpiStart += ncmpiEllipsoidCollCount[i];
      }
      ncmpiCount = ncmpiEllipsoidCollCount[domainId()];
      if(ncmpiStart >= ncmpiCollCountMax) {
        if(ncmpiCount == 0) {
          ncmpiStart = 0;
        } else {
          mTerm(1, AT_,
                "Error in LPT::writeCollData "
                "(m_includeEllipsoids)!! ncmpiStart >= "
                "ncmpiCollCountMax but ncmpiCount != 0");
        }
      }

      // Arrays to hold collisions data.
      MFloatScratchSpace ncmpiCollTimeStep(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiCollTimeStep");
      MIntScratchSpace ncmpiPart0(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiPart0");
      MIntScratchSpace ncmpiPart1(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiPart1");
      MFloatScratchSpace ncmpiP0SemiMinorAxis(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP0SemiMinorAxis");
      MFloatScratchSpace ncmpiP1SemiMinorAxis(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP1SemiMinorAxis");
      MFloatScratchSpace ncmpiP0AspectRatio(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP0AspectRatio");
      MFloatScratchSpace ncmpiP1AspectRatio(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP1AspectRatio");
      MFloatScratchSpace ncmpiP0Dens(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP0Dens");
      MFloatScratchSpace ncmpiP1Dens(ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiP1Dens");
      MFloatScratchSpace ncmpiCollPos(3 * ncmpiEllipsoidCollCount[domainId()], AT_, "ncmpiCollPos");
      MInt ncmpiCountId = 0;


      // Write the data from the queue into the array.
      while(!collQueueEllipsoid.empty()) {
        ncmpiCollTimeStep[ncmpiCountId] = collQueueEllipsoid.front().collTimeStep;
        ncmpiPart0[ncmpiCountId] = collQueueEllipsoid.front().part0;
        ncmpiPart1[ncmpiCountId] = collQueueEllipsoid.front().part1;
        ncmpiP0SemiMinorAxis[ncmpiCountId] = collQueueEllipsoid.front().semiMinorAxis0;
        ncmpiP1SemiMinorAxis[ncmpiCountId] = collQueueEllipsoid.front().semiMinorAxis1;
        ncmpiP0AspectRatio[ncmpiCountId] = collQueueEllipsoid.front().aspectRatio0;
        ncmpiP1AspectRatio[ncmpiCountId] = collQueueEllipsoid.front().aspectRatio1;
        ncmpiP0Dens[ncmpiCountId] = collQueueEllipsoid.front().dens0;
        ncmpiP1Dens[ncmpiCountId] = collQueueEllipsoid.front().dens1;
        ncmpiCollPos[(3 * ncmpiCountId)] = collQueueEllipsoid.front().collPosX;
        ncmpiCollPos[(3 * ncmpiCountId) + 1] = collQueueEllipsoid.front().collPosY;
        ncmpiCollPos[(3 * ncmpiCountId) + 2] = collQueueEllipsoid.front().collPosZ;

        collQueueEllipsoid.pop();
        ++ncmpiCountId;
      }

      // Put the arrays into the Netcdf file.
      parallelIo.setOffset(ncmpiCount, ncmpiStart);
      parallelIo.writeArray(&ncmpiCollTimeStep[0], "collTimeStep");
      parallelIo.writeArray(&ncmpiPart0[0], "part0");
      parallelIo.writeArray(&ncmpiPart1[0], "part1");
      parallelIo.writeArray(&ncmpiP0SemiMinorAxis[0], "p0semiMinorAxis");
      parallelIo.writeArray(&ncmpiP1SemiMinorAxis[0], "p1semiMinorAxis");
      parallelIo.writeArray(&ncmpiP0AspectRatio[0], "p0aspectRatio");
      parallelIo.writeArray(&ncmpiP1AspectRatio[0], "p1aspectRatio");
      parallelIo.writeArray(&ncmpiP0Dens[0], "p0dens");
      parallelIo.writeArray(&ncmpiP1Dens[0], "p1dens");

      ncmpiCount *= 3;
      ParallelIo::size_type ncmpi3Start = 3 * ncmpiStart;

      parallelIo.setOffset(ncmpiCount, ncmpi3Start);
      parallelIo.writeArray(&ncmpiCollPos[0], "collPos");

      // Free memory.
      while(!collQueueEllipsoid.empty()) {
        collQueueEllipsoid.pop();
      }
    }
  }
}

// Explicit instantiations for 2D and 3D
template class ParticleCollision<3>;
