// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <cmath>
#include "UTIL/functions.h"
#include "lpt.h"

/**
 * \fn void LPT::readModelProps()
 *
 * \brief Read model configuration properties from configuration file.
 *
 * \author Sven Berger, Tim Wegmann
 * \date   August 2015
 */
template <MInt nDim>
void LPT<nDim>::readModelProps() {
  TRACE();

  m_log << "reading model properties...";

  m_restart = false;
  m_restart = Context::getSolverProperty<MBool>("restartFile", solverId(), AT_, &m_restart);

  m_restartTimeStep = Context::getSolverProperty<MInt>("restartTimeStep", solverId(), AT_);

  m_maxNoParticles = 1E6;
  m_maxNoParticles = Context::getSolverProperty<MInt>("maxNoParticles", solverId(), AT_, &m_maxNoParticles);

  /*! \page propertyPageLPT
    \section particleExchangeBufferSize
    <code>MInt LPT::m_particleExchangeBufferSize</code>\n
    default = <code>1000</code> \n \n
    Specifies the size of the buffer used for transfer between cpus \n
    of particles.\n
    Keywords: <i>PARTICLE</i>
 */
  m_exchangeBufferSize =
      Context::getSolverProperty<MInt>("particleExchangeBufferSize", solverId(), AT_, &m_exchangeBufferSize);


  /*! \page propertyPageLPT
    \section particleSizeLimit
    <code>MFloat LPT::m_sizeLimit</code>\n
    default = 1E-12 \n \n
    Remove particles smaller than this size \n
    Keywords: <i>PARTICLE</i>
    */
  m_sizeLimit = 1E-12;
  m_sizeLimit = Context::getSolverProperty<MFloat>("particleSizeLimit", solverId(), AT_, &m_sizeLimit);

  /*! \page propertyPageLPT
    \section spawnParticles
    <code>MBool LPT::m_spawnParticles</code>\n
    default = false \n \n
    Spawn particles \n
    Keywords: <i>PARTICLE</i>
    */
  m_spawnParticles = false;
  m_spawnParticles = Context::getSolverProperty<MBool>("spawnParticles", solverId(), AT_, &m_spawnParticles);

  if(m_spawnParticles) {
    m_log << "Particle spawning has been activated." << std::endl;
  }

  /*! \page propertyPageLPT
    \section sprayPrimaryBUp
    <code>MBool LPT::m_activePrimaryBUp</code>\n
    default = <code>false</code> \n \n
    Activate primary break-up \n
    Keywords: <i>PARTICLE</i>
  */
  m_activePrimaryBUp = false;
  m_activePrimaryBUp = Context::getSolverProperty<MBool>("sprayPrimaryBUp", solverId(), AT_, &m_activePrimaryBUp);
  /*! \page propertyPageLPT
    \section spraySecondaryBUp
    <code>MBool LPT::m_activeSecondaryBUp</code>\n
    default = <code>false</code> \n \n
    Activate secondary break-up \n
    Keywords: <i>PARTICLE</i>
  */
  m_activeSecondaryBUp = false;
  m_activeSecondaryBUp = Context::getSolverProperty<MBool>("spraySecondaryBUp", solverId(), AT_, &m_activeSecondaryBUp);
  /*! \page propertyPageLPT
    \section particleIncludeEllipsoids
    <code>MBool LPT::m_particleIncludeEllipsoids</code>\n
    default = false \n \n
    Exclude (false) or include (true) ellipsoidal particles \n
    in addition to spheres.\n
    Keywords: <i>PARTICLE</i>
    */
  m_ellipsoids = false;
  m_ellipsoids = Context::getSolverProperty<MBool>("particleIncludeEllipsoids", solverId(), AT_, &m_ellipsoids);

  /*! \page propertyPageLPT
      \section particleCollisions
      <code>MInt LPT::m_collisions</code>\n
      default = <code>0</code> \n \n
      Specify the method for collision detection and execution. \n
      (Note that for ellipsoids the collision detection is different, \n
      it is currently hard-coded as option 6) \n \n
      Possible values: \n
      <ul>
      <li><code>0</code> No collision detection</li>
      <li><code>1</code> Retroactive detection; elastic collisions</li>
      <li><code>2</code> Proactive detection; elastic collision</li>
      <li><code>3</code> Retroactive detection; coagulation</li>
      <li><code>4</code> Proactive detection; coagulation</li>
      <li><code>5</code> Retroactive detection; no particle interaction</li>
      <li><code>6</code> Proactive detection; no particle interaction</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_collisions = Context::getSolverProperty<MInt>("particleCollisions", solverId(), AT_, &m_collisions);

  if(m_collisions > 0) {
    m_collisionModel = std::unique_ptr<ParticleCollision<nDim>>(
        new ParticleCollision<nDim>(this, m_collisions, solverId(), domainId(), noDomains(), mpiComm()));
  }

  /*! \page propertyPageLPT
  \section particleWallCollisions
  <code>MInt LPT::m_particleWallCollisions</code>\n
  default = <code>true</code> \n \n
  Specify if wall-collisions should be taken into account. \n \n \n
  Possible values: true/false \n
  Keywords: <i>PARTICLE</i>
  */
  m_wallCollisions = Context::getSolverProperty<MBool>("particleWallCollisions", solverId(), AT_, &m_wallCollisions);

  /*! \page propertyPageLPT
      \section particleInitializationMethod
      <code>MInt LPT::m_initializationMethod</code>\n
      default = <code>0</code> \n \n
      Choose the method for the particle initialization: \n
      <ul>
        <li><code>0</code> Initializes particleInitNoPartperCell particles per cell</li>
        <li><code>1</code> Initializes particles in the cell center of all inflow cells</li>
        <li>
          <code>2/3</code> Initializes particles according to the input file part.txt.
          Specify particle diameter, density and position
        </li>
        <li><code>4</code> Same as 2/3 but now specify also particle velocity</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_initializationMethod =
      Context::getSolverProperty<MInt>("particleInitializationMethod", solverId(), AT_, &m_initializationMethod);

  /*! \page propertyPageLPT
      \section particleDrag
      <code>MInt LPT::m_dragModelType</code>\n
      default = <code>1</code> \n \n
      Choose drag formulation (spheres): \n
      <ul>
      <li><code>0</code> No drag force</li>
      <li><code>1</code> Linear Stokes drag</li>
      <li><code>2</code> Schiller-Naumann correlation \f$ f_D = 1 + 0.15 Re^{0.687} \f$</li>
      <li><code>3</code> Pinsky-Khain correlation \f$ f_D = 1 + 0.17 Re^{0.632} + 10^{-6}Re^{2.25} \f$</li>
      <li><code>4</code> Baumgarten, 2005</li>
      <li><code>4</code> Rowe, 1961</li>
      <li><code>5</code> Drag limited to 10% of fluid fraction</li>
      </ul>
      Choose drag formulation (ellipsoids): \n
      <ul>
      <li><code>0</code> No drag force</li>
      <li><code>1</code> Linear</li>
      <li><code>2</code> Schiller-Naumann correlation \f$ f_D = 1 + 0.15 Re^{0.687} \f$</li>
      <li><code>1</code> Drag correlation for prolates by Konstantin et al. (2020)</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_dragModelType = Context::getSolverProperty<MInt>("particleDrag", solverId(), AT_, &m_dragModelType);

  /*! \page propertyPageLPT
      \section particleLift
      <code>MInt LPT::m_liftModelType</code>\n
      default = <code>1</code> \n \n
      Choose lift formulation: \n
      <ul>
      <li><code>0</code> No lift force</li>
      <li><code>2</code> Linear</li>
      <li><code>1</code> Lift correlation for prolates by Konstantin et al. (2020)</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_liftModelType = Context::getSolverProperty<MInt>("particleLift", solverId(), AT_, &m_liftModelType);

  /*! \page propertyPageLPT
      \section particleTorque
      <code>MInt LPT::m_torqueModelType</code>\n
      default = <code>1</code> \n \n
      Choose torque formulation: \n
      <ul>
      <li><code>0</code> No torque</li>
      <li><code>1</code> Linear</li>
      <li><code>1</code> Torque correlation for prolates by Konstantin et al. (2020)</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_torqueModelType = Context::getSolverProperty<MInt>("particleTorque", solverId(), AT_, &m_torqueModelType);

  /*! \page propertyPageLPT
      \section motionEquation
      <code>MInt LPT::m_motionEquationType</code>\n
      default = <code>0</code> \n \n
      Choose correctorStep formulation: \n
      <ul>
        <li>
          <code>0</code> Motion equation based on the analytical solution of the simplified motion equation under the
            assumption of constant drag. See Crowe et al., Multiphase Flows with Droplets and Particles, 1998.
        </li>
        <li>
          <code>1</code> Non dimensional version based on C. Baumgarten, Mixture formation in internal combustion
          engines, 2005
        </li>
        <li><code>2</code> Particles move with fluid velocity</li>
      </ul>

      Keywords: <i>PARTICLE</i>
   */
  m_motionEquationType = Context::getSolverProperty<MInt>("motionEquation", solverId(), AT_, &m_motionEquationType);

  /*! \page propertyPageLPT
      \section particleRespawn
      <code>MInt LPT::m_particleRespawn</code>\n
      default = <code>false</code> \n \n
      Specifies whether particles that leave the domain are reintroduced. \n
      <ul>
      <li><code>0</code> Particles are not reintroduced</li>
      <li><code> >0 </code> Particles are reintroduced in a plane perpendicular
      to direction (m_particleRespawn - 1)</li>
      </ul>
      Keywords: <i>PARTICLE</i>
   */
  m_respawn = Context::getSolverProperty<MBool>("particleRespawn", solverId(), AT_, &m_respawn);

  if(m_respawn) {
    /*! \page propertyPageLPT
        \section particleRespawnPlane
        <code>MFloat LPT::m_respawnPlane</code>\n
        default = nullptr \n \n
        The coordinate plane in which particles are respawned. \n
        Keywords: <i>PARTICLE</i>
     */
    m_respawnPlane = Context::getSolverProperty<MFloat>("particleRespawnPlane", solverId(), AT_);

    // create a vector of cells which are used for the "respawning" of particles
    MFloat halfLength = NAN;
    for(MInt i = 0; i < noInternalCells(); i++) {
      if(!c_isLeafCell(i)) continue;
      MInt level = c_level(i);
      halfLength = c_cellLengthAtLevel(level + 1);
      if(fabs(c_coordinate(i, 0) - m_respawnPlane) < halfLength) {
        auto weightedTimes = static_cast<MInt>(pow(pow(2, nDim), (maxRefinementLevel() - level)));
        for(MInt count = 0; count < weightedTimes; ++count) {
          m_respawnCells.push_back(i);
        }
      }
    }
    m_log << m_respawnCells.size() << " respawn cells found in this solver" << std::endl;
  }

  /*! \page propertyPageLPT
      \section particlePartOffset
      <code>MFloat LPT::m_xCutOff</code>\n
      default = <code>-1000.0</code> \n \n
      Specify a lower limit for the x coordinate of particle to be \n
      written. This property is useful to allow for an "relaxation \n
      distance" for newly introduced particles. \n
      Note that at the default value no limit is used. \n
      Keywords: <i>PARTICLE</i>
   */
  m_xCutOff = Context::getSolverProperty<MFloat>("particlePartOffset", solverId(), AT_, &m_xCutOff);

  /*! \page propertyPageLPT
      \section particleMomentumCoupling
      <code>MBool LPT::m_momentumCoupling</code>\n
      default = <code>false</code> \n \n
      Activate momentum coupling. \n
      Keywords: <i>PARTICLE</i>
   */
  m_momentumCoupling = false;
  m_momentumCoupling =
      Context::getSolverProperty<MBool>("particleMomentumCoupling", solverId(), AT_, &m_momentumCoupling);
  if(m_momentumCoupling) {
    m_log << "Momentum coupling has been activated." << std::endl;
  }

  /*! \page propertyPageLPT
      \section particleHeatCoupling
      <code>MBool LPT::m_heatCoupling</code>\n
      default = <code>false</code> \n \n
      Activate heat coupling. \n
      Keywords: <i>PARTICLE</i>
   */
  m_heatCoupling = false;
  m_heatCoupling = Context::getSolverProperty<MBool>("particleHeatCoupling", solverId(), AT_, &m_heatCoupling);

  if(m_heatCoupling) {
    m_log << "Heat coupling has been activated." << std::endl;
  }

  /*! \page propertyPageLPT
    \section particleMassCoupling
    <code>MBool LPT::m_massCoupling</code>\n
    default = <code>false</code> \n \n
    Activate mass coupling. \n
    Keywords: <i>PARTICLE</i>
 */
  m_massCoupling = false;
  m_massCoupling = Context::getSolverProperty<MBool>("particleMassCoupling", solverId(), AT_, &m_massCoupling);
  if(m_massCoupling) {
    m_log << "Mass coupling has been activated." << std::endl;
  }

  /*! \page propertyPageLPT
      \section particleEvaporation
      <code>MBool LPT::m_evaporation</code>\n
      default = <code>false</code> \n \n
      Activate evaporation modeling. \n
      Keywords: <i>PARTICLE</i>
   */
  m_evaporation = false;
  m_evaporation = Context::getSolverProperty<MBool>("particleEvaporation", solverId(), AT_, &m_evaporation);

  if(m_evaporation) {
    m_log << "Evaporation modeling has been activated." << std::endl;
  }

  /*! \page propertyPageLPT
    \section maxNoBndryCells
    <code> MInt LPT::m_maxNoBndryCells</code>  \n
    default = <code>""</code>\n \n
    Maximum number of boundary cells to be allocated. (Used for wall-collision)
    <ul> <li>Any positiv integer</li>
    </ul>
    Keywords: <i>LPT, ALLOCATION, BNDRY</i>
   */
  m_maxNoBndryCells = Context::getSolverProperty<MInt>("maxNoBndryCells", solverId(), AT_);


  /*! \page propertyPageLPT
  \section particleAdapRange
  <code>MFloat cell-range</code>\n
  default = {2.0, 0.0}\n
  Range around particles, which should be refined during adaptation.
  First value is on maxRefinementLevel, second value is on all lower levels\n
  Keywords: <i>PARTICLE, ADAPTATION</i>
  */
  MInt range[2] = {2, 0};
  for(MInt i = 0; i < 2; i++) {
    range[i] = Context::getSolverProperty<MInt>("particleAdapRange", solverId(), AT_, &range[i], i);
  }

  mAlloc(m_bandWidth, maxRefinementLevel() + 1, "m_bandWidth", 0, AT_);
  m_bandWidth[maxRefinementLevel() - 1] = range[0];
  for(MInt i = maxRefinementLevel() - 2; i >= 0; i--) {
    m_bandWidth[i] = (m_bandWidth[i + 1] / 2) + 1 + range[1];
  }

  m_innerBound = Context::getSolverProperty<MInt>("innerBound", solverId(), AT_, &m_innerBound);

  MFloat referenceTemperature = 273.15;
  /*! \page propertyPageLPT
    \section referenceTemperature
    <code>MFloat LPT::m_referenceTemperature </code>\n
    default = <code>273.15</code>\n \n
   Reference temperature \f$ T_{\mathrm{ref}}\f$
   Used to scale the Sutherland's constant as follows: \f$ S/T_{\mathrm{ref}} \f$
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>LPT, VARIABLES</i>
  */
  referenceTemperature =
      Context::getSolverProperty<MFloat>("referenceTemperature", m_solverId, AT_, &referenceTemperature);

  m_sutherlandConstant = 110.4;
  /*! \page propertyPageLPT
    \section sutherlandConstant
    <code>MFloat LPT::m_sutherlandConstant </code>\n
    default = <code>110.4 K</code>\n \n
    Sutherland's constant. Used by Sutherland's law.
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>LPT, VARIABLES</i>
  */
  m_sutherlandConstant =
      Context::getSolverProperty<MFloat>("sutherlandConstant", solverId(), AT_, &m_sutherlandConstant);

  m_sutherlandConstant /= referenceTemperature;
  m_sutherlandPlusOne = m_sutherlandConstant + F1;

  /*! \page propertyPageLPT
      \section cfl
    <code>MFloat LPT::m_cfl </code>\n
    default = <code>1.0</code>\n \n
    Courant number C - Factor of the CFL condition \n \n
    possible values are:
    <ul>
      <li>positive floating point values</li>
    </ul>
    Keywords: <i> LPT, TIME_INTEGRATION </i>
  */

  m_cfl = 1.0;
  m_cfl = Context::getSolverProperty<MFloat>("cfl", solverId(), AT_, &m_cfl);

  m_timeStepComputationInterval = m_restart ? -1 : 0;
  /*! \page propertyPageLPT
  \section timeStepVar
  <code>MInt LPT::m_timeStepComputationInterval </code>\n
  default = <code> -1 </code>\n
  Specifies on which interval the time-step will be recomputed.
  <ul>
  <li> -1 -> default (when restart): never - read from restart file, requires m_restart == true. </li>
  <li> 0  -> default (when no restart): once at the beginning </li>
  <li> 1  -> every time-step </li>
  <li> 2  -> every two time-steps </li>
  <li> N  -> every Ntime-steps </li>
  </ul>
  Keywords: <i>LPT, TIME_INTEGRATION</i>
*/
  m_timeStepComputationInterval =
      Context::getSolverProperty<MInt>("timeStepComputationInterval", solverId(), AT_, &m_timeStepComputationInterval);


  m_nonBlockingComm = Context::getSolverProperty<MBool>("nonBlockingComm", solverId(), AT_, &m_nonBlockingComm);

  if(m_nonBlockingComm) {
    m_nonBlockingStage = 0;
    if(m_ellipsoids) {
      cerr0 << "LPT-Warning: nonBlocking currently only works with either spherical or ellipsoidal particles "
            << "--> disabling exchange of spherical particles" << std::endl;
    }
  }

  /*! \page propertyPageLPT
    \section weightBaseCell
    <code>MFloat LPT::m_weightBaseCell</code>\n
    default = <code>0.0</code>\n \n
    Weight applied for any lpt-cell during static weight computation for domain decomposition during
    balance.
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */
  m_weightBaseCell = 0.0;
  m_weightBaseCell = Context::getSolverProperty<MFloat>("weightBaseCell", solverId(), AT_, &m_weightBaseCell);
  /*! \page propertyPageLPT
    \section weightBaseCell
    <code>MFloat LPT::m_weightBaseCell</code>\n
    default = <code>0.0</code>\n \n
    Weight applied for any lpt-cell during static weight computation for domain decomposition during
    balance, good value could be 0.05.
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */

  m_weightLeafCell = 0.0;
  m_weightLeafCell = Context::getSolverProperty<MFloat>("weightLeafCell", solverId(), AT_, &m_weightLeafCell);
  /*! \page propertyPageLPT
    \section weightBaseCell
    <code>MFloat LPT::m_weightBaseCell</code>\n
    default = <code>0.0</code>\n \n
    Weight applied for any lpt-cell during static weight computation for domain decomposition during
    balance, good value could be 0.1
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */

  m_weightParticleCell = 0.1;
  m_weightParticleCell =
      Context::getSolverProperty<MFloat>("weightParticleCell", solverId(), AT_, &m_weightParticleCell);
  /*! \page propertyPageLPT
    \section weightParticle
    <code>MFloat LPT::m_weightParticle</code>\n
    default = <code>0.0</code>\n \n
    Weight applied for each particle in a lpt-leaf cell during static weight computation
    for domain decomposition during balance, good value could be 0.1
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */
  m_weightParticle = 1.0;
  m_weightParticle = Context::getSolverProperty<MFloat>("weightParticle", solverId(), AT_, &m_weightParticle);
  /*! \page propertyPageLPT
    \section weightSpawnCell
    <code>MFloat LPT::m_weightSpawnCell</code>\n
    default = <code>0.0</code>\n \n
    Weight applied for the cell holding the injector during static weight computation
    for domain decomposition during balance, good value could be 5.0
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */
  m_weightSpawnCell = 0.0;
  m_weightSpawnCell = Context::getSolverProperty<MFloat>("weightSpawnCell", solverId(), AT_, &m_weightSpawnCell);
  /*! \page propertyPageLPT
    \section weightMulitSolverFactor
    <code>MFloat LPT::m_weightMulitSolverFactor</code>\n
    default = <code>1.0</code>\n \n
    Multi-solver factor for lpt weight in comparison with the computational load of the other
    numerical schems -> 1.0 for single solver
    Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */
  m_weightMulitSolverFactor = 1.0;
  m_weightMulitSolverFactor =
      Context::getSolverProperty<MFloat>("weightMulitSolverFactor", solverId(), AT_, &m_weightMulitSolverFactor);
  /*! \page propertyPageLPT
  \section limitWeights
  <code>MBool LPT::m_limitWeights</code>\n
  default = <code>false</code>\n \n
  Weight applied for any lpt-cell during static weight computation for domain decomposition during
  balance.
  Keywords: <i> LPT, WEIGHTING, BALANCE</i>
  */
  m_limitWeights = false;
  m_limitWeights = Context::getSolverProperty<MBool>("limitDLBWeights", solverId(), AT_, &m_limitWeights);

  m_weightSourceCells = false;
  m_weightSourceCells =
      Context::getSolverProperty<MBool>("weightLPTSourceCells", solverId(), AT_, &m_weightSourceCells);

  m_domainIdOutput = false;
  m_domainIdOutput = Context::getSolverProperty<MBool>("domainIdOutput", m_solverId, AT_, &m_domainIdOutput);

  // LPT sleep time in milli-seconds
  m_sleepLPT = -1;
  m_sleepLPT = Context::getSolverProperty<MInt>("sleepTimeLPT", m_solverId, AT_, &m_sleepLPT);

  m_skipLPT = false;
  m_skipLPT = Context::getSolverProperty<MBool>("skipLPT", m_solverId, AT_, &m_skipLPT);

  m_log << "finished" << std::endl;
}

/**
 * \fn void LPT::readSpawnProps()
 *
 *  \brief Read properties related to creating particles from configuration file.
 *
 *  \author Sven Berger
 *  \date   March 2016
 */
template <MInt nDim>
void LPT<nDim>::readSpawnProps() {
  /*! \page propertyPageLPT
      \section spawnSeed
      <code>MLong LPT::m_spawnSeed</code>\n
      default = Default Seed (5489u)\n \n
      Initialize PRNG with given seed. \n
      Keywords: <i>PARTICLE</i>
      */
  m_spawnSeed = 5489U;
  if(Context::propertyExists("spawnSeed", solverId())) {
    m_spawnSeed = (MLong)Context::getSolverProperty<MInt>("spawnSeed", solverId(), AT_);
  }
  // initialise the psydo-randon number generator with the given seed
  // this is done identically on all ranks and for all random numbers!
  m_PRNGSpawn.seed(m_spawnSeed);
  m_PRNGRespawn.seed(m_spawnSeed);

  /*! \page propertyPageLPT
      \section particleDiameter
      <code>MFloat LPT::m_spawnDiameter</code>\n
      default = 0.00001\n
      REQUIRED\n
      Spawn particles with the given diameter \n
      Keywords: <i>PARTICLE</i>
      */
  m_spawnDiameter = 0.00001;
  m_spawnDiameter = Context::getSolverProperty<MFloat>("spawnDiameter", solverId(), AT_, &m_spawnDiameter);

  /*! \page propertyPageLPT
      \section spawnParticlesInitVelo
      <code>MFloat LPT::m_spawnParticlesInitVelo</code>\n
      default = 0\n
      REQUIRED\n
      Spawn particles with the given speed \n
      Keywords: <i>PARTICLE</i>
       */
  m_spawnVelocity = Context::getSolverProperty<MFloat>("spawnParticlesInitVelo", solverId(), AT_, &m_spawnVelocity);

  /*! \page propertyPageLPT
      \section spawnParticlesCount
      <code>MInt LPT::m_spawnParticlesCount</code>\n
      default = 1 \n
      REQUIRED\n
      Spawn the given number of particles per time unit \n
      Keywords: <i>PARTICLE</i>
       */
  m_spawnParticlesCount =
      Context::getSolverProperty<MInt>("spawnParticlesCount", solverId(), AT_, &m_spawnParticlesCount);

  /*! \page propertyPageLPT
      \section spawnDistSigmaCoeff
      <code>MFloat LPT::m_spawnDistSigmaCoeff</code>\n
      default = 1\n
      Sigma coefficient for the normal distribution. \n
      Keywords: <i>PARTICLE</i>
       */
  m_spawnDistSigmaCoeff =
      Context::getSolverProperty<MFloat>("spawnDistSigmaCoeff", solverId(), AT_, &m_spawnDistSigmaCoeff);

  m_spawnParticlesConeAngle =
      Context::getSolverProperty<MFloat>("spawnParticlesConeAngle", solverId(), AT_, &m_spawnParticlesConeAngle);
  /*! \page propertyPageLPT
    \section spawnDir
    <code>MFloat* LPT::m_spawnDir</code>\n
    default = {1, 0, 0} \n \n
    Normal vector of the particle spawn direction. \n
    Keywords: <i>PARTICLE</i>
  */

  for(MInt i = 0; i < nDim; i++) {
    m_spawnDir[i] = 0;
    if(i == 0) m_spawnDir[i] = 1.0;
  }

  if(Context::propertyExists("spawnDir", solverId())) {
    if(Context::propertyLength("spawnDir", solverId()) != nDim) {
      TERMM(1, "Need to give a Coordinate for every dimension");
    }

    for(MInt i = 0; i < nDim; i++) {
      m_spawnDir[i] = Context::getSolverProperty<MFloat>("spawnDir", solverId(), AT_, i);
    }
  }

  /*! \page propertyPage1
    \section engineSetup
    <code>MBool LPT::m_engineSetup</code>\n
    default = false\n
    Trigger specific engine features! \n
    Keywords: <i>PARTICLE</i>
  */
  m_engineSetup = false;
  m_engineSetup = Context::getSolverProperty<MBool>("engineSetup", solverId(), AT_, &m_engineSetup);
}

/**
 * \fn void LPT::readEllipsoidProps()
 *
 *  \brief Read properties specific to ellipsoids.
 *
 *  \author Sven Berger
 *  \date   September 2016
 */
template <MInt nDim>
void LPT<nDim>::readEllipsoidProps() {
  /*! \page propertyPageLPT
      \section particleEllipsoidRandomOrientation
      <code>MInt LPT::m_ellipsoidRandomOrientation</code>\n
      default = <code>1</code> \n \n
      if any number unequal 0, ellipoids are initialized with random orientation,
      if equal 0 b-axis is orientated in x-direction \n
      Keywords: <i>PARTICLE</i>
   */
  m_ellipsoidRandomOrientation = Context::getSolverProperty<MInt>("particleEllipsoidRandomOrientation", solverId(), AT_,
                                                                  &m_ellipsoidRandomOrientation);
  /*! \page propertyPage1
      \section particleEllipsoidRandomOrientationSeed
      <code>MLong LPT::m_ellipsoidRandomOrientationSeed</code>\n
      default = Default Seed (4865U)\n \n
      Keywords: <i>PARTICLE</i>
      */
  m_ellipsoidRandomOrientationSeed = 4865U;
  if(Context::propertyExists("particleEllipsoidRandomOrientationSeed", solverId())) {
    m_ellipsoidRandomOrientationSeed =
        (MLong)Context::getSolverProperty<MInt>("particleEllipsoidRandomOrientationSeed", solverId(), AT_);
  }
}

/**
 * \fn void LPT::readMomentumCouplingProps()
 *
 *  \brief Read properties specific to the momentum coupling.
 *
 *  \author Sven Berger
 *  \date   September 2016
 */
template <MInt nDim>
void LPT<nDim>::readMomentumCouplingProps() {
  /*! \page propertyPageLPT
      \section particleCouplingRedist
      <code>MBool LPT::m_couplingRedist</code>\n
      default = <code>false</code> \n \n
      Activate coupling redistribution. \n
      Keywords: <i>PARTICLE</i>
   */
  m_couplingRedist = false;
  m_couplingRedist = Context::getSolverProperty<MBool>("particleCouplingRedist", solverId(), AT_, &m_couplingRedist);

  if(m_couplingRedist) {
    /*! \page propertyPageLPT
    \section particleNoRedistLayer
    <code>MBool LPT::m_noRedistLayer</code>\n
    default = <code>1</code> \n \n
    Activate momentum coupling redistribution area. \n
    Keywords: <i>PARTICLE</i>
 */
    m_noRedistLayer = Context::getSolverProperty<MInt>("particleNoRedistLayer", solverId(), AT_, &m_noRedistLayer);

    if(m_noRedistLayer > grid().noHaloLayers() && noDomains() > 1) {
      mTerm(-1, AT_, "ERROR: Number of halo layers needs to be higher than the redistribution area!");
    }

    m_log << "Momentum redistribution is active with a redistribution area of " << m_noRedistLayer << std::endl;
  }
}

/**
 *  \fn void LPT::findSpawnCellId()
 *
 *  \brief Find spawn cellId
 *
 *  \author Sven Berger
 *  \date   December 2019
 */
template <MInt nDim>
void LPT<nDim>::findSpawnCellId() {
  const MBool redetermination = m_spawnCellId >= 0;
  const MInt lastSpawnCellId = m_spawnCellId;
  const MInt lastSpawnDomainId = m_spawnDomainId;

  m_spawnCellId = -1;
  /*! \page propertyPageLPT
      \section spawnCoordinates
      <code>MFloat* spawnCoord</code>\n
      default = N/A \n \n
      Spawn particles in the given coordinate \n
      Keywords: <i>PARTICLE</i>
       */
  if(Context::propertyExists("spawnCoordinates", solverId())) {
    if(Context::propertyLength("spawnCoordinates", solverId()) != nDim) {
      TERMM(1, "Need to give a Coordinate for every dimension");
    }

    for(MInt i = 0; i < nDim; i++) {
      m_spawnCoord[i] = Context::getSolverProperty<MFloat>("spawnCoordinates", solverId(), AT_, i);
    }

    m_spawnCellId = grid().findContainingLeafCell(&m_spawnCoord[0]);

    // If we are running parallel we need to find the cell which is
    // closest to the given coordinates on all domains.
    if(noDomains() > 1) {
      MInt spawnDomain = domainId();
      MInt minDomain = -1;
      if(m_spawnCellId < 0 || a_isHalo(m_spawnCellId)) {
        spawnDomain = std::numeric_limits<MInt>::max();
      }

      // If more than one domain has cells with the same distance only
      // spawn on the domain with the lowest Id
      MPI_Allreduce(&spawnDomain, &minDomain, 1, MPI_INT, MPI_MIN, mpiComm(), AT_, "spawnDomain", "minDomain");

      if(minDomain != domainId()) {
        m_spawnCellId = -1;
      }

      m_spawnDomainId = minDomain;

      m_log << " LPT selected spawn-cell " << m_spawnCellId << " on domain " << m_spawnDomainId << std::endl;
      if(domainId() == m_spawnDomainId) {
        std::cerr << " LPT selected spawn-cell " << m_spawnCellId << " on domain " << m_spawnDomainId
                  << " with coordinates " << c_coordinate(m_spawnCellId, 0) << " " << c_coordinate(m_spawnCellId, 1)
                  << " " << c_coordinate(m_spawnCellId, nDim - 1) << std::endl;
      }
    }

    if(redetermination) {
      std::cerr << " LPT spawn-cellId was redetermined";
      if(m_spawnCellId < 0) {
        std::cerr << " and spawnDomain changed from " << domainId() << " to " << m_spawnDomainId;
        std::cerr << " left-over mass was : " << m_particleResiduum << std::endl;
      } else {
        std::cerr << std::endl;
      }
    }
    if(lastSpawnCellId < 0 && m_spawnCellId > 0 && lastSpawnDomainId > 0) {
      std::cerr << "Revieved left-over mass is : " << m_particleResiduum << std::endl;
    }
    if((m_spawnDomainId + 1 > noDomains() && noDomains() > 1) || (m_spawnCellId < 0 && noDomains() == 1)) {
      TERMM(-1, "Invalid spawn coordinates, values are probably outside the fluid-domain/box?");
    }
  }
}
