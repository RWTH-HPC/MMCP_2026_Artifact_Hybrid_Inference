// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptspray.h"
#include <cmath>
#include "IO/tomlutils.h"
#include "lpt.h"
#include "lptlib.h"

using namespace maia::lpt;
using namespace std;
using namespace maia;

template <MInt nDim>
MFloat SprayModel<nDim>::s_lengthFactor{};
template <MInt nDim>
LPT<nDim>* SprayModel<nDim>::s_backPtr = nullptr;

/// \fn void SprayModel<nDim>::init(const MInt spawnCellId, LPT<nDim>* _particle)
/// \brief Spray model initialisation
///
/// \author Sven Berger
/// \date   October 2015
/// \param[in] spawnCellId Location at which particles are created.
template <MInt nDim>
void SprayModel<nDim>::init(LPT<nDim>* _particle) {
  m_log << "       Initialising spray model" << endl;

  s_backPtr = _particle;
  s_lengthFactor = s_backPtr->m_partList[0].s_lengthFactor;

  m_activePrimaryBUp = s_backPtr->m_activePrimaryBUp;
  m_activeSecondaryBUp = s_backPtr->m_activeSecondaryBUp;

  if(!s_backPtr->m_restart) {
    s_backPtr->m_particleResiduum = 0.0;
  }

  // allocate spray data
  mAlloc(m_needleLiftTime, 2, "needleLiftTime", AT_);
  mAlloc(m_sprayAngle, 2, "m_sprayAngle", F0, AT_);
  mAlloc(m_injData, m_injDataSize, "m_injData", AT_);

  readSprayProperties();

  logDistributionFunction();
}


/**
 *  \brief Read spray model and injector properties
 *  \author Tim Wegmann
 *  \date   January 2023
 */
template <MInt nDim>
void SprayModel<nDim>::readSprayProperties() {
  // debug parameter to validate injection/conservation
  m_injStopTimeStep = -1;
  m_injStopTimeStep = Context::getSolverProperty<MInt>("injStopTimeStep", solverId(), AT_, &m_injStopTimeStep);

  m_injStartTimeStep = -1;
  m_injStartTimeStep = Context::getSolverProperty<MInt>("injStartTimeStep", solverId(), AT_, &m_injStartTimeStep);

  if(m_activePrimaryBUp) {
    m_needleLiftTime[0] = 0.01;
    m_needleLiftTime[1] = 0.02;
    /*! \page propertyPageLPT
    \section injectorNeedleLiftTime
    <code>MFloat SprayModel::m_needleLiftTime</code>\n
    default = 0 \n \n
    Dimensionless time for the needle to reach maximum lift (e.g. maximum flow-rate) \n
    Non-dimensionalised by the injection duration(which can also be non-dimensional!)
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    if(Context::propertyExists("injectorNeedleLiftTime", solverId())) {
      m_needleLiftTime[0] = Context::getSolverProperty<MFloat>("injectorNeedleLiftTime", solverId(), AT_, 0);
      m_needleLiftTime[1] = Context::getSolverProperty<MFloat>("injectorNeedleLiftTime", solverId(), AT_, 1);
      m_useNeedleLiftTime = true;
    }

    /*! \page propertyPageLPT
    \section injectorInjectionTime
    <code>MFloat SprayModel::m_injDuration</code>\n
    default = 0.00055 \n \n
    Time during which fuel is injected.\n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_injDuration = 0.00055;
    m_injDuration = Context::getSolverProperty<MFloat>("injectorInjectionTime", solverId(), AT_, &m_injDuration);

    /*! \page propertyPageLPT
    \section injectionCrankAngle
    <code>MFloat SprayModel::m_injectionCrankAngle</code>\n
    default = -1 \n \n
    Meaning, injection crank angle is not applied and the injection starts at timeStep zero!
    Otherwise specify an injection crank-angle, valid values are in the range of 0-720. \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_injectionCrankAngle = -1;
    m_injectionCrankAngle =
        Context::getSolverProperty<MFloat>("injectionCrankAngle", solverId(), AT_, &m_injectionCrankAngle);

    // NOTE: the Strouhal number must be different from the FV Strouhal-number
    //        if a different non-dimensionalisation is used in the LPT-solver!
    m_Strouhal = Context::getSolverProperty<MFloat>("Strouhal", solverId(), AT_, &m_Strouhal);

    m_initialCad = Context::getSolverProperty<MFloat>("initialCrankAngle", solverId(), AT_, &m_initialCad);

    /*! \page propertyPageLPT
      \section sprayPrimRosinRammler
      <code>MInt SprayModel::m_primRosinRammler</code>\n
      default = <code>true</code> \n \n
      Use Rosin-Rammler distribution to determine initial droplet size. \n
      Keywords: <i>PARTICLE</i>
    */
    m_primRosinRammler = true;
    m_primRosinRammler =
        Context::getSolverProperty<MBool>("sprayPrimRosinRammler", solverId(), AT_, &m_primRosinRammler);

    /*! \page propertyPageLPT
    \section injectorType
    <code>MFloat SprayModel::m_injectorType</code>\n
    default = single-hole fullCone injector \n \n
    Specific typ of injector to be used for the simulation, falls into the following categories:
    single-hole, multi-hole,hollow-cone \n
    Options are: "FULLCONE", "HOLLOWCONE", "MULTIHOLE", "MULTIHOLE_OPT"
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    MString injectorType = "FULLCONE";
    injectorType = Context::getSolverProperty<MString>("injectorType", solverId(), AT_, &injectorType);
    m_injectorType = string2enum(injectorType);

    // injector categories:
    m_multiHoleInjector = false;
    m_singleHoleInjector = false;
    m_hollowConeInjector = false;

    switch(m_injectorType) {
      case MULTIHOLE:
      case MULTIHOLE_OPT:
      case MULTIHOLE_TME: {
        m_multiHoleInjector = true;
        break;
      }
      case FULLCONE: {
        m_singleHoleInjector = true;
        break;
      }
      case HOLLOWCONE: {
        m_hollowConeInjector = true;
        break;
      }
      default: {
        mTerm(1, AT_, "Unknown injector-type!");
      }
    }

    /*! \page propertyPageLPT
    \section injectorNozzleDiameter
    <code>MFloat SprayModel::m_injectorNozzleDiameter</code>\n
    default = NA \n \n
    Injector exit diameter \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_injectorNozzleDiameter = Context::getSolverProperty<MFloat>("injectorNozzleDiameter", solverId(), AT_);

    // read injector type specific properties
    if(m_multiHoleInjector) {
      /*! \page propertyPageLPT
      \section injectorDiameter
        <code>MFloat SprayModel::m_injectorDiameter</code>\n
        default = NA \n \n
        Diameter of multi-hole injector.
        If specified as an array, the radius to the central injector point can be specified
        for individually for each hole..\n
        Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      if(Context::propertyLength("injectorDiameter") == 1) {
        mAlloc(m_injectorDiameter, 1, "injectorDiameter", AT_);
        m_injectorDiameter[0] = Context::getSolverProperty<MFloat>("injectorDiameter", solverId(), AT_);
      } else {
        const MInt length = Context::propertyLength("injectorDiameter");
        mAlloc(m_injectorDiameter, length + 1, "injectorDiameter", AT_);
        m_injectorDiameter[0] = 0.0;
        for(MInt i = 0; i < length; i++) {
          m_injectorDiameter[i + 1] = Context::getSolverProperty<MFloat>("injectorDiameter", solverId(), AT_, i);
          m_injectorDiameter[0] = mMax(m_injectorDiameter[0], m_injectorDiameter[i + 1]);
        }
      }

      /*! \page propertyPageLPT
      \section injectorOrificeAngle
        <code>MFloat SprayModel::m_injectorOrificeAngle</code>\n
        default = NA \n \n
        Orientation angle of the individual orifices in degree.
        default = 37 degree (spray-G injector)
        Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      if(!Context::propertyExists("injectorOrificeAngle") || Context::propertyLength("injectorOrificeAngle") == 1) {
        mAlloc(m_injectorOrificeAngle, 1, "injectorOrificeAngle", AT_);
        m_injectorOrificeAngle[0] = 37.0;
        m_injectorOrificeAngle[0] =
            Context::getSolverProperty<MFloat>("injectorOrificeAngle", solverId(), AT_, &m_injectorOrificeAngle[0]);
      } else {
        const MInt length = Context::propertyLength("injectorOrificeAngle");
        mAlloc(m_injectorOrificeAngle, length, "injectorOrificeAngle", AT_);
        for(MInt i = 0; i < length; i++) {
          m_injectorOrificeAngle[i] = Context::getSolverProperty<MFloat>("injectorOrificeAngle", solverId(), AT_, i);
        }
      }
      /*! \page propertyPageLPT
      \section orificePositionAngle
        <code>MFloat SprayModel::m_orificePositionAngle</code>\n
        default = NA \n \n
        Orientation angle of the individual orifices in degree.
        default = number-holes/360 degree (spray-G injector)
        Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      if(Context::propertyExists("orificePositionAngle")) {
        const MInt length = Context::propertyLength("orificePositionAngle");
        mAlloc(m_orificePositionAngle, length, "orificePositionAngle", AT_);
        for(MInt i = 0; i < length; i++) {
          m_orificePositionAngle[i] = Context::getSolverProperty<MFloat>("orificePositionAngle", solverId(), AT_, i);
        }
      }
    }

    /*! \page propertyPageLPT
    \section injectorOrficeSize
    <code>MFloat SprayModel::m_injectorOrficeSize</code>\n
    default = NA \n \n
    Injector outer orfices size of the injector\n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_injectorOrficeSize = Context::getSolverProperty<MFloat>("injectorOrficeSize", solverId(), AT_);


    if(m_primRosinRammler && m_singleHoleInjector) {
      // determine diameter factors for the rosin rammler initial droplet size distriution (IDSD)
      // default values are fraction of the nominal injector nozzle diameter for spray-A as in:
      // LARGE EDDY SIMULATION OF HIGH-VELOCITY FUEL SPRAYS:
      // STUDYING MESH RESOLUTION AND BREAKUP MODEL EFFECTS FOR SPRAY A
      // A. Wehrfritz, V. Vuorinen, O. Kaario, & M. Larmi
      // Atomization and Sprays, 23 (5): 419–442 (2013)
      m_RosinRammlerMin = 90;
      m_RosinRammlerMean = 15;
      m_RosinRammlerMax = 5;
      m_RosinRammlerSpread = 3;

    } else if(m_primRosinRammler && m_hollowConeInjector) {
      // determine diameter factors for the rosin rammler initial droplet size distriution (IDSD)
      // default values are fraction of the liquid sheet length for hollow-cone injectors
      m_RosinRammlerMin = 3.33333333;
      m_RosinRammlerMean = 1;
      m_RosinRammlerMax = 0.33333333;
      m_RosinRammlerSpread = 3;
    } else if(m_primRosinRammler && m_multiHoleInjector) {
      m_RosinRammlerMin = 3.33333333;
      m_RosinRammlerMean = 1;
      m_RosinRammlerMax = 0.33333333;
      m_RosinRammlerSpread = 3;
    }

    // read specified values from the property file
    if(m_primRosinRammler) {
      /*! \page propertyPageLPT
      \section RosinRammlerSpread
      <code>MFloat SprayModel::m_RosinRammlerSpread</code>\n
      default = 3 \n \n
      Spread of Rosin-Rammler distribution function\n
      Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      m_RosinRammlerSpread =
          Context::getSolverProperty<MFloat>("RosinRammlerSpread", solverId(), AT_, &m_RosinRammlerSpread);

      /*! \page propertyPageLPT
      \section RosinRammlerMin
      <code>MFloat SprayModel::m_RosinRammlerMin</code>\n
      default = 90 \n \n
      Factor of nominal nozzle diameter for Rosin-Rammler min distribution parameter\n
      Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      m_RosinRammlerMin = Context::getSolverProperty<MFloat>("RosinRammlerMin", solverId(), AT_, &m_RosinRammlerMin);

      /*! \page propertyPageLPT
      \section injectorNozzleDiameter
      <code>MFloat SprayModel::m_RosinRammlerMax</code>\n
      default = 5 \n \n
      Factor of nominal nozzle diameter for Rosin-Rammler max distribution parameter\n
      Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      m_RosinRammlerMax = Context::getSolverProperty<MFloat>("RosinRammlerMax", solverId(), AT_, &m_RosinRammlerMax);

      /*! \page propertyPageLPT
      \section injectorNozzleDiameter
      <code>MFloat SprayModel::m_RosinRammlerMean</code>\n
      default = 15 \n \n
      Factor of nominal nozzle diameter for Rosin-Rammler mean distribution parameter \n
      Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      m_RosinRammlerMean = Context::getSolverProperty<MFloat>("RosinRammlerMean", solverId(), AT_, &m_RosinRammlerMean);
    }

    /*! \page propertyPageLPT
    \section spawnDistSigmaCoeff
    <code>MFloat SprayModel::m_injectionDistSigmaCoeff</code>\n
    default = 1\n
    Sigma coefficient for the normal distribution. \n
    Keywords: <i>PARTICLE</i>
    */
    m_injectionDistSigmaCoeff =
        Context::getSolverProperty<MFloat>("spawnDistSigmaCoeff", solverId(), AT_, &m_injectionDistSigmaCoeff);

    /*! \page propertyPageLPT
        \section partEmittDist
        <code>MInt SprayModel::m_partEmittDist</code>\n
        default = None \n \n
        Distributes particles according to a given distribution\n
        Keywords: <i>PARTICLE</i>
         */
    MString partDist = "PART_EMITT_DIST_NONE";
    m_partEmittDist = string2enum(Context::getSolverProperty<MString>("partEmittDist", solverId(), AT_, &partDist));

    MFloat diameter = m_injectorNozzleDiameter;
    if(m_multiHoleInjector) {
      diameter = m_injectorDiameter[0];
    } else if(m_singleHoleInjector) {
      diameter = m_injectorOrficeSize;
    }

    if(diameter < (2 * s_backPtr->c_cellLengthAtLevel(s_backPtr->minLevel()))) {
      m_broadcastInjected = false;
      if(s_backPtr->domainId() == 0) {
        cerr << "Injection broadcast not necessary!" << endl;
      }
    }

    /*! \page propertyPageLPT
    \section injectionSpeed
    <code>MFloat SprayModel::m_injectionSpeed</code>\n
    default = NA \n \n
    Injection speed \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_injectionSpeed = Context::getSolverProperty<MFloat>("injectionSpeed", solverId(), AT_, &m_injectionSpeed);

    ASSERT(m_injectionSpeed > 0, "");

    if(Context::propertyExists("sprayInjectionRate", solverId())) {
      const MInt numInjections = Context::propertyLength("sprayInjectionRate", solverId());

      if(numInjections == 1) {
        m_currentInjectionRate = Context::getSolverProperty<MFloat>("sprayInjectionRate", solverId(), AT_);

      } else {
        m_injectionRateList.resize(numInjections);
        m_injectionTimings.resize(numInjections);

        for(MInt i = 0; i < numInjections; i++) {
          /*! \page propertyPageLPT
          \section sprayInjectionRate
          <code>MFloatTensor SprayModel::m_injectionRateList</code>\n
          default = NA \n \n
          Time variable injection rate \n
          Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
          */
          m_injectionRateList[i] = Context::getSolverProperty<MFloat>("sprayInjectionRate", solverId(), AT_, i);

          /*! \page propertyPageLPT
          \section sprayInjectionTiming
          <code>MFloatTensor SprayModel::m_injectionTimings</code>\n
          default = NA \n \n
          Time of provided injection rate \n
          Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
          */
          m_injectionTimings[i] = Context::getSolverProperty<MFloat>("sprayInjectionTiming", solverId(), AT_, i);
        }
      }
    } else {
      ASSERT(!m_activePrimaryBUp, "");
    }

    /*! \page propertyPageLPT
        \section injectorDir
        <code>MFloat* SprayModel::injectorDir</code>\n
        default = {0, 0, 1} \n \n
        Normal vector of the injector direction. \n
        Keywords: <i>PARTICLE</i>
    */

    for(MInt i = 0; i < nDim; i++) {
      m_injectorDir[i] = 0;
      if(i == nDim - 1) m_injectorDir[i] = 1.0;
    }

    if(Context::propertyExists("injectorDir", solverId())) {
      if(Context::propertyLength("injectorDir", solverId()) != nDim) {
        TERMM(1, "Need to give a Coordinate for every dimension");
      }

      for(MInt i = 0; i < nDim; i++) {
        m_injectorDir[i] = Context::getSolverProperty<MFloat>("injectorDir", solverId(), AT_, i);
      }
    }

    if(m_hollowConeInjector) {
      /*! \page propertyPageLPT
        \section sprayPrimaryMaxPrimParcels
        <code>MInt SprayModel::m_maxNoPrimParcels</code>\n
        default = <code>1</code> \n \n
        Set the maximum number (and spawn locations) of parcels generated per time step. \n
        Keywords: <i>PARTICLE</i>
      */
      m_maxNoPrimParcels = 1;
      m_maxNoPrimParcels =
          Context::getSolverProperty<MInt>("sprayPrimaryMaxPrimParcels", solverId(), AT_, &m_maxNoPrimParcels);

      /*! \page propertyPageLPT
         \section sprayPrimaryParcelsPerInj
         <code>MInt SprayModel::m_primParcelsPerInj</code>\n
         default = <code>1000000</code> \n \n
         Number of parcels to be generated per injection \n
         Keywords: <i>PARTICLE</i>
      */
      m_primParcelsPerInj = 1000000;
      m_primParcelsPerInj =
          Context::getSolverProperty<MInt>("sprayPrimaryParcelsPerInj", solverId(), AT_, &m_primParcelsPerInj);

      /*! \page propertyPageLPT
        \section injectorAngularGap
        default = -1 \n \n
        Injector angular gap which can be used as a condition for the primary breakup. \n
        -1 Meaning that the angularGap is not used!
        Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
      */
      m_angularGap = Context::getSolverProperty<MFloat>("injectorAngularGap", solverId(), AT_, &m_angularGap);
    }

    if(m_singleHoleInjector || m_multiHoleInjector) {
      /*! \page propertyPageLPT
      \section sprayPrimaryBUpParcelSize
      <code>MInt SprayModel::m_primBrkUpParcelSize</code>\n
      default = <code>1</code> \n
      Activate primary break-up. \n
      Valid values: \n
       Any positive integer values. \n
      Keywords: <i>PARTICLE</i>
      */
      m_primBrkUpParcelSize = 1;
      m_primBrkUpParcelSize =
          Context::getSolverProperty<MInt>("sprayPrimaryBUpParcelSize", solverId(), AT_, &m_primBrkUpParcelSize);
      ASSERT(m_primBrkUpParcelSize > 0, "ERROR: Invalid parcel size.");
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  /// model configuration properties
  //////////////////////////////////////////////////////////////////////////////////////////////////


  if(Context::propertyLength("sprayAngle", solverId()) == 1) {
    m_sprayAngle[0] = Context::getSolverProperty<MFloat>("sprayAngle", solverId(), AT_, &m_sprayAngle[0]);
    m_sprayAngle[1] = Context::getSolverProperty<MFloat>("sprayAngle", solverId(), AT_, &m_sprayAngle[0]);
  } else {
    m_sprayAngle[0] = Context::getSolverProperty<MFloat>("sprayAngle", solverId(), AT_, 0);
    m_sprayAngle[1] = Context::getSolverProperty<MFloat>("sprayAngle", solverId(), AT_, 1);
  }


  m_sprayConeAngle = Context::getSolverProperty<MFloat>("sprayConeAngle", solverId(), AT_, &m_sprayConeAngle);

  /*! \page propertyPageLPT
    \section spawnSeed
    <code>MLong LPT::m_spawnSeed</code>\n
    default = Default Seed (5489u)\n \n
    Initialize PRNG with given seed. \n
    Keywords: <i>PARTICLE</i>
    */
  m_spraySeed = 5489U;
  if(Context::propertyExists("spawnSeed", solverId())) {
    m_spraySeed = (MLong)Context::getSolverProperty<MInt>("spawnSeed", solverId(), AT_);
  }

  if(m_activePrimaryBUp) m_PRNGPrimBreakUp.seed(m_spraySeed);
  if(m_activeSecondaryBUp) m_PRNGSecBU.seed(m_spraySeed);

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // secondary breakup properties
  //////////////////////////////////////////////////////////////////////////////////////////////////
  if(m_activeSecondaryBUp) {
    /*! \page propertyPageLPT
    \section sprayBUPEventOutput
    <code>MInt SprayModel::m_secBUDisplayMessage</code>\n
    default = <code>false</code> \n \n
    Activate break-up event output\n
    Keywords: <i>PARTICLE</i>
    */
    m_secBUDisplayMessage = false;
    m_secBUDisplayMessage =
        Context::getSolverProperty<MBool>("spraysecBUPOutput", solverId(), AT_, &m_secBUDisplayMessage);

    m_RTsecBreakUp = true;
    m_RTsecBreakUp = Context::getSolverProperty<MBool>("RayleighTaylorBreakUp", solverId(), AT_, &m_RTsecBreakUp);

    m_KHsecBreakUp = true;
    m_KHsecBreakUp = Context::getSolverProperty<MBool>("KelvinHelmholtzBreakUp", solverId(), AT_, &m_KHsecBreakUp);

    m_RTsecBreakUpLength = true;
    m_RTsecBreakUpLength =
        Context::getSolverProperty<MBool>("RayleighTaylorBreakUpLength", solverId(), AT_, &m_RTsecBreakUpLength);


    /*! \page propertyPageLPT
    \section predictivePRNG
    <code>MBool SprayModel::m_predictivePRNG</code>\n
    default = true \n \n
    Trigger a predictive PRNG which is independend of the number of used Ranks and particle ordering!
    However, in this case the PRNG needs to be initialised each time before generating a number for
    the secondary Breakup. This is computationally very expensive and thus not recommended for large
    production runs!
    \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
     */
    m_predictivePRNG = true;
    m_predictivePRNG = Context::getSolverProperty<MBool>("predictablePRNGSecBU", solverId(), AT_, &m_predictivePRNG);

    // initialise PRNG
    m_PRNGSecBU.discard(globalTimeStep);

    /*! \page propertyPageLPT
    \section secondaryBreakup
    <code>MFloat SprayModel::m_sprayAngleKH</code>\n
    default = 1 \n \n
    Defines spray angle descripting for KH-secondary break-up. Options:
    -1: const. property based spray angle
     2: const. property based spray angle with near nozzle correction
    0-1: calculated spray angle with A1 value given by secSprayAngleMode!
    \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_sprayAngleKH = 2;
    m_sprayAngleKH = Context::getSolverProperty<MFloat>("sprayAngleKH", solverId(), AT_, &m_sprayAngleKH);

    // default values for secondary break-up are based on:
    // LARGE EDDY SIMULATION OF HIGH-VELOCITY FUEL SPRAYS:
    // STUDYING MESH RESOLUTION AND BREAKUP MODEL EFFECTS FOR SPRAY A
    // A. Wehrfritz, V. Vuorinen, O. Kaario, & M. Larmi
    // Atomization and Sprays, 23 (5): 419–442 (2013)


    /*! \page propertyPageLPT
    \section sprayWeberLimitRT
    <code>MFloat SprayModel::m_weberLimitRT</code>\n
    default = 300 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_weberLimitRT = Context::getSolverProperty<MFloat>("sprayWeberLimitRT", solverId(), AT_, &m_weberLimitRT);

    /*! \page propertyPageLPT
    \section sprayB0
    <code>MFloat SprayModel::m_B0</code>\n
    default = 0.61 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_B0 = Context::getSolverProperty<MFloat>("sprayB0", solverId(), AT_, &m_B0);

    /*! \page propertyPageLPT
    \section sprayB1
    <code>MFloat SprayModel::m_B1</code>\n
    default = 40 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // valid range of B1 is 10-60 according to Reitz1999
    m_B1 = Context::getSolverProperty<MFloat>("sprayB1", solverId(), AT_, &m_B1);

    /*! \page propertyPageLPT
    \section sprayC3
    <code>MFloat SprayModel::m_c3</code>\n
    default = 0.1 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // other values used by Sven Berger:
    // 0.15 AIA-Nozzle MPa 10
    // 0.20 AIA-Nozzle MPa 20
    m_CRT = Context::getSolverProperty<MFloat>("sprayCRT", solverId(), AT_, &m_CRT);

    /*! \page propertyPageLPT
    \section sprayCT
    <code>MFloat SprayModel::m_CT</code>\n
    default = 1.0 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_CT = Context::getSolverProperty<MFloat>("sprayCT", solverId(), AT_, &m_CT);

    /*! \page propertyPageLPT
    \section sprayWeberLimitKH
    <code>MFloat SprayModel::m_weberLimitKH</code>\n
    default = 6 \n \n
    KHRT-model property \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_weberLimitKH = Context::getSolverProperty<MFloat>("sprayWeberLimitKH", solverId(), AT_, &m_weberLimitKH);

    /*! \page propertyPageLPT
    \section voumeKHLimit
    <code>MFloat SprayModel::m_massLimit</code>\n
    default = 0.03 \n \n
    KHRT-model property \n
    Dimensionless mass fraction of shedded particle to particle \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // In literature sometimes denoted as B2
    m_massLimit = Context::getSolverProperty<MFloat>("sprayMassLimit", solverId(), AT_, &m_massLimit);

    /*! \page propertyPageLPT
    \section sprayMaxRTChilds
    <code>MInt SprayModel::m_maxRTChildParcels</code>\n
    default = 3 \n \n
    Maximum number of parcels generated by RT \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    m_maxRTChildParcels = Context::getSolverProperty<MInt>("sprayMaxRTChilds", solverId(), AT_, &m_maxRTChildParcels);

    /*! \page propertyPageLPT
    \section sprayRTDiameterMode
    <code>MInt SprayModel::m_</code>\n
    default = 0 \n \n
    Choosen version of RT child diameter \n
    0: use RT-diameter
    1: use d_old^(2/3)* d_RT^(1/3)
    2: use Rosin-Rammer distribution function for the diameter
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // NOTE: version 1 is following: "Implementation and validation of a Lagrangian spray model using
    //      experimental data of the EVN Spray G injector" Horacia J. Aguerre, Norberto M. Nigro
    //      Computers and Fluids 190(2019) 30-48
    m_RTDiameterMode = Context::getSolverProperty<MInt>("sprayRTDiameterMode", solverId(), AT_, &m_RTDiameterMode);

    /*! \page propertyPageLPT
    \section sprayCbu
    <code>MFloat SprayModel::m_Cbu</code>\n
    default = m_B1/2 \n \n
    KHRT-model property \n
    Model parameter for liquid breakup length
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // default based on Reitz et al.
    //"Modeling Spray Atomization with the Kelvin-Helmholtz/Rayleigh-Taylor
    // Hybrid Model" Atomization and Sprays 1999
    m_Cbu = m_B1 / 2;
    m_Cbu = Context::getSolverProperty<MFloat>("sprayCbu", solverId(), AT_, &m_Cbu);

    MFloat breakUpLength = m_Cbu * m_injectorNozzleDiameter * sqrt(material().ambientDensityRatio());

    cerr0 << "Liquid break-up length for KH-RT model is " << breakUpLength << endl;
  }
}


/// \brief Perform a spray injection
/// \author Sven Berger
/// \date   November 2015
template <MInt nDim>
void SprayModel<nDim>::injection(const MFloat dt) {
  // injection has not started yet, this also means that zero particles should be present!
  if(m_injectionCrankAngle > maia::math::crankAngle(s_backPtr->m_time, m_Strouhal, m_initialCad, 0)) {
    ASSERT(s_backPtr->a_noParticles() == 0, "");
    if(domainId() == 0) {
      cerr << "Before injection" << m_injectionCrankAngle << " "
           << maia::math::crankAngle(s_backPtr->m_time, m_Strouhal, m_initialCad, 0) << endl;
    }
    for(MInt i = 0; i < m_injDataSize; i++) {
      m_injData[i] = 0;
    }
    return;
  }

  m_timeSinceSOI += dt;


  MBool endOfInjection = false;
  if(m_injStopTimeStep > -1 && globalTimeStep > m_injStopTimeStep) endOfInjection = true;
  if(m_timeSinceSOI > m_injDuration) endOfInjection = true;
  if(m_injStartTimeStep > -1 && globalTimeStep < m_injStartTimeStep) endOfInjection = true;

  // injection has already ended
  if(endOfInjection) {
    cerr0 << "After injection " << m_timeSinceSOI << " " << m_injDuration << endl;

    for(MInt i = 0; i < m_injDataSize; i++) {
      m_injData[i] = 0;
    }
    if(s_backPtr->m_spawnCellId > -1) {
      m_injData[0] = m_timeSinceSOI;
    }
    return;
  }


  // injection location is on another domain
  if(s_backPtr->m_spawnCellId < 0) {
    ASSERT(s_backPtr->noDomains() > 1, "");

    // receive injected particles instead
    if(m_broadcastInjected) {
      s_backPtr->recvInjected();
    }
    for(MInt i = 0; i < m_injDataSize; i++) {
      m_injData[i] = 0;
    }

  } else {
    // injection on this domain
    ASSERT(domainId() == s_backPtr->m_spawnDomainId, "");

    // update injection rate
    updateInjectionRate();

    // current number of particles
    const MInt prevNo = s_backPtr->a_noParticles();

    // inject new particles
    primaryBreakUp(dt);

    if(m_broadcastInjected) {
      s_backPtr->broadcastInjected(prevNo);
    }
  }
}

/// \brief Inject new particles from the injector exit.
///
/// \author Sven Berger
/// \date   November 2015
template <MInt nDim>
void SprayModel<nDim>::primaryBreakUp(const MFloat dt) {
  TRACE();

  ASSERT(domainId() == s_backPtr->m_spawnDomainId, "");
  ASSERT(s_backPtr->m_spawnCellId > -1, "");

  // invalidate injection data:
  MInt noInjParticles = -1;
  MInt noInjDroplets = -1;
  MFloat injMass = NAN;
  MFloat injParticleDiameter = NAN;
  static MFloat lastD = -1;

  // 1) ramp-up and other general injection parameters:

  MFloat injectionProgress = m_timeSinceSOI / m_injDuration;


  // the rampFactor is used to control/ramp the massflowrate and velocity
  // during the injector opening and closing
  // which linearly increases/decreases to the full injection velocity massflowrate
  MFloat rampFactor = 1.0;
  if(injectionProgress < m_needleLiftTime[0]) {
    rampFactor = injectionProgress / m_needleLiftTime[0];
  } else if(injectionProgress > (1 - m_needleLiftTime[1])) {
    rampFactor = 1 - (injectionProgress - (1 - m_needleLiftTime[1])) / m_needleLiftTime[1];
  }

  // theoretically injected mass during this timestep
  // current injection rate + left over mass from previous timestep
  injMass = m_currentInjectionRate * rampFactor * dt + s_backPtr->m_particleResiduum;

  // considering ramp-factor for parcel count
  MInt noDroplets = mMax(1, (MInt)(m_primBrkUpParcelSize * rampFactor));

  // 2) injector dependant injection parameters
  //   and injection of new particles
  if(m_singleHoleInjector) { //(Spray A)

    if(m_primRosinRammler) {
      // droplet diameter distribution function

      const MFloat diameterMin = m_injectorNozzleDiameter / m_RosinRammlerMin;
      const MFloat diameterMean = m_injectorNozzleDiameter / m_RosinRammlerMean;
      const MFloat diameterMax = m_injectorNozzleDiameter / m_RosinRammlerMax;

      // mass of the smallest possible droplet
      const MFloat minMass = sphericalMass(diameterMin, material().T());

      noInjParticles = 0;
      noInjDroplets = 0;

      MFloat remainingMass = injMass;
      injMass = 0;

      while(remainingMass > minMass) {
        if(lastD > 0) {
          injParticleDiameter = lastD;
        } else {
          injParticleDiameter =
              rosinRammler(diameterMin, diameterMean, diameterMax, m_RosinRammlerSpread, randomPrimBreakUp(1));
        }

        const MFloat particleMass = noDroplets * sphericalMass(injParticleDiameter, material().T());
        const MFloat tempMass = remainingMass - particleMass;

        // limit mass to be positive
        if(tempMass < 0) {
          lastD = injParticleDiameter;
          s_backPtr->m_particleResiduum = remainingMass;
          break;
        }

        // the initial cone diameter needs to be reduced by half of the droplet diameter
        // since the droplet is supposed to be fully within the cone
        const MFloat coneDiameter = m_injectorOrficeSize - 0.5 * injParticleDiameter;
        lastD = -1;
        injectParticles(1, s_backPtr->m_spawnCoord, injParticleDiameter, material().densityRatio(),
                        m_injectionSpeed * rampFactor, m_sprayConeAngle, coneDiameter, false, 0.0, noDroplets, 1);

        noInjParticles++;
        noInjDroplets += noDroplets;
        injMass += particleMass;
        remainingMass -= particleMass;
      }

    } else {
      // blob-injection with all particles at the injector center

      injParticleDiameter = m_injectorNozzleDiameter;
      noInjDroplets = 1;

      const MFloat initParticleMass = sphericalMass(injParticleDiameter, material().T());


      MFloat noNewParticles = (m_currentInjectionRate * dt) / initParticleMass + s_backPtr->m_particleResiduum;
      noInjParticles = floor(noNewParticles);

      s_backPtr->m_particleResiduum = noNewParticles - noInjParticles;

      injectParticles(noInjParticles, s_backPtr->m_spawnCoord, injParticleDiameter, material().densityRatio(),
                      m_injectionSpeed, m_sprayConeAngle, 0.0, false, 0.0, 1, 1);
    }
  } else if(m_hollowConeInjector) {
    // hollow cone injector, such as the TMFB Injector
    // This implementation of the hollow cone injector model of:
    // Spray Modeling for Outwardly-Opening Hollow-Cone Injectors, Sim, Badra et al., SAE, 2016

    // maximum number of parcels to be generated during this time step
    const MInt parcelTargetNo = m_primParcelsPerInj / m_injDuration * dt;

    // limit the number of parcels to be injected to maxNoPrimParcels
    noInjParticles = std::min(parcelTargetNo, m_maxNoPrimParcels);

    // injection velocity is scaled by this factor
    const MFloat scaledInjectionV = rampFactor * m_injectionSpeed;

    // thickness of the fluid flow within the angular gap
    MFloat initialLiquidSheetThickness = NAN;

    if(m_angularGap > -1) {
      initialLiquidSheetThickness = max(m_angularGap * rampFactor, m_primMinDropletSize);
    } else {
      initialLiquidSheetThickness =
          max(0.5 * m_injectorNozzleDiameter
                  - sqrt(M_PI * POW2(m_injectorNozzleDiameter) * material().density() * scaledInjectionV
                         - 4.0 * m_currentInjectionRate * rampFactor)
                        / (2.0 * sqrt(M_PI) * sqrt(material().density() * scaledInjectionV)),
              m_primMinDropletSize);
    }

    // initial particle size when considered as string like coneshaped
    // starting from the injector angular gap
    injParticleDiameter = 4.0 / M_PI * initialLiquidSheetThickness;
    if(m_primRosinRammler) {
      const MFloat diameterMin = injParticleDiameter / m_RosinRammlerMin;
      const MFloat diameterMean = injParticleDiameter / m_RosinRammlerMean;
      const MFloat diameterMax = injParticleDiameter / m_RosinRammlerMax;
      injParticleDiameter =
          rosinRammler(diameterMin, diameterMean, diameterMax, m_RosinRammlerSpread, randomPrimBreakUp(1));
    }

    // mass of a particle to be added
    const MFloat initialDropletMass = sphericalMass(injParticleDiameter, material().T());

    const MFloat newDroplets = injMass / initialDropletMass;

    noInjDroplets = floor(newDroplets / noInjParticles);

    if(noInjDroplets > 0) {
      const MFloat injectorGapDiameter = m_injectorNozzleDiameter - 0.5 * injParticleDiameter;
      s_backPtr->m_particleResiduum = (newDroplets - noInjParticles * noInjDroplets) * initialDropletMass;

      injectParticles(noInjParticles, s_backPtr->m_spawnCoord, injParticleDiameter, material().densityRatio(),
                      scaledInjectionV, m_sprayConeAngle, injectorGapDiameter * s_lengthFactor, true, m_angularGap,
                      noInjDroplets, 1.0);

      injMass = noInjParticles * noInjDroplets * initialDropletMass;

    } else {
      injMass = 0;
      noInjParticles = 0;
    }

  } else if(m_multiHoleInjector) {
    // injectors with multiple injection holes
    // MULTIHOLE     : ECN Spray-G
    // MULTIHOLE_Opt : optimized 7-hole Spray-G
    // MULTIHOLE_TME : 5-hole Injection of the TME

    MInt noHoles = -1;
    switch(m_injectorType) {
      case MULTIHOLE: {
        noHoles = 8;
        break;
      }
      case MULTIHOLE_OPT: {
        noHoles = 7;
        break;
      }
      case MULTIHOLE_TME: {
        noHoles = 5;
        break;
      }
      default:
        mTerm(1, AT_, "Unknown multi-hole injector-type!");
    }

    static constexpr MFloat C_a = 0.65;
    // NOTE: Wehrfritz et. al. propose to use the nominal injector nozzle diameter!
    //      otherwise the area-coefficient of ~0.65 can be used

    const MFloat initialConeAngle = m_sprayConeAngle;
    // as in:
    //"Validation of a comprehensive computational fluid dynamics methodology to predict the direct
    // injection process of gasoline sprays using Spray G experimental data"
    // Davide Paredi , Tommaso Lucchini , Gianluca D’Errico, Angelo Onorati,
    // Lyle Pickett and Joshua Lacey
    // International J of Engine Research 2020, Vol. 21(1) 199–216

    ASSERT(m_primRosinRammler, "");

    const MFloat effectiveDiam = sqrt(C_a) * m_injectorNozzleDiameter;
    const MFloat diameterMin = effectiveDiam / m_RosinRammlerMin;
    const MFloat diameterMean = effectiveDiam / m_RosinRammlerMean;
    const MFloat diameterMax = effectiveDiam / m_RosinRammlerMax;

    MFloat remainingMass = injMass;
    const MFloat minMass = sphericalMass(diameterMin, material().T());

    injMass = 0;
    noInjParticles = 0;
    noInjDroplets = 0;

    std::function<MFloat(MInt)> holePosition = [&](MInt id) {
      const MFloat angleDelta = 2 * M_PI / 8;
      switch(m_injectorType) {
        case MULTIHOLE: {
          return (id * angleDelta);
        }
        case MULTIHOLE_OPT: {
          // 0: x=0,  z=-1
          // 2: x=-1, z=0
          // 4: x=0,  z=1
          // 6: x=1,  z=0
          if(id > 1) {
            id = id + 1;
          }
          return id * angleDelta;
        }
        case MULTIHOLE_TME: {
          return m_orificePositionAngle[id] / 180 * M_PI;
        }
        default:
          mTerm(1, AT_, "Unknown multi-hole injector-type!");
      }
    };

    std::function<MFloat(MInt)> injectorDiameter = [&](MInt id) {
      switch(m_injectorType) {
        case MULTIHOLE:
        case MULTIHOLE_OPT: {
          return m_injectorDiameter[0];
        }
        case MULTIHOLE_TME: {
          return m_injectorDiameter[id + 1];
        }
        default:
          mTerm(1, AT_, "Unknown multi-hole injector-type!");
      }
    };

    // return angle [degree] of the hole orientation
    std::function<MFloat(MInt)> holeAngle = [&](MInt id) {
      switch(m_injectorType) {
        case MULTIHOLE:
        case MULTIHOLE_OPT: {
          return 37.0;
        }
        case MULTIHOLE_TME: {
          return m_injectorOrificeAngle[id];
        }
        default:
          mTerm(1, AT_, "Unknown multi-hole injector-type!");
      }
    };

    while(remainingMass > minMass) {
      if(lastD > 0) {
        injParticleDiameter = lastD;
      } else {
        injParticleDiameter =
            rosinRammler(diameterMin, diameterMean, diameterMax, m_RosinRammlerSpread, randomPrimBreakUp(1));
      }

      const MFloat particleMass = noHoles * noDroplets * sphericalMass(injParticleDiameter, material().T());

      const MFloat tempMass = remainingMass - particleMass;

      if(tempMass < 0) {
        lastD = injParticleDiameter;
        s_backPtr->m_particleResiduum = remainingMass;
        break;
      }

      lastD = -1;
      injectParticles(noHoles, s_backPtr->m_spawnCoord, injParticleDiameter, material().densityRatio(),
                      m_injectionSpeed, initialConeAngle, m_injectorOrficeSize, holePosition, injectorDiameter,
                      holeAngle, false, noDroplets, noHoles);

      injMass += particleMass;
      noInjParticles += noHoles;
      noInjDroplets += noHoles * noDroplets;
      remainingMass -= particleMass;
    }
  }

  // 3) store injection data for possible post-processing:
  m_injData[0] = m_timeSinceSOI;
  m_injData[1] = injectionProgress;
  m_injData[2] = rampFactor;
  m_injData[3] = noInjParticles;
  m_injData[4] = noInjDroplets;
  m_injData[5] = injParticleDiameter;
  m_injData[6] = injMass;

  array<MFloat, 3> injMomentum = {};
  for(MInt i = 0; i < 3; i++) {
    injMomentum[i] = 0;
  }
  MFloat injEnergy = 0;
  for(MInt i = 0; i < noInjParticles; i++) {
    const MInt id = s_backPtr->a_noParticles() - 1 - i;
    const MFloat mass = sphericalMass(s_backPtr->m_partList[id].m_diameter, s_backPtr->m_partList[id].m_temperature)
                        * s_backPtr->m_partList[id].m_noParticles;
    MFloat velMagSquared = 0;
    for(MInt j = 0; j < nDim; j++) {
      injMomentum[j] += mass * s_backPtr->m_partList[id].m_velocity[j];
      velMagSquared += POW2(s_backPtr->m_partList[id].m_velocity[j]);
    }
    injEnergy += mass * s_backPtr->m_material->cp(s_backPtr->m_partList[id].m_temperature)
                 * s_backPtr->m_partList[id].m_temperature * 1 / s_backPtr->m_material->gammaMinusOne();
    injEnergy += 0.5 * mass * velMagSquared;
  }
  m_injData[7] = injMomentum[0];
  m_injData[8] = injMomentum[1];
  m_injData[9] = injMomentum[2];
  m_injData[10] = injEnergy;

  cerr << globalTimeStep << " " << m_timeSinceSOI << " injMass " << injMass << " injParticles " << noInjParticles
       << endl;
}

/// \brief Atomization of particles.
/// \author Sven Berger
/// \date   November 2015
template <MInt nDim>
void SprayModel<nDim>::secondaryBreakUp(const MFloat dt) {
  TRACE();

  m_noRTsecBreakUp = 0;
  m_noKHsecBreakUp = 0;

  if(!m_RTsecBreakUp && !m_KHsecBreakUp) return;

  // use const to avoid breakup of particles which have just been created
  const MInt noPart = s_backPtr->a_noParticles();

  // increase size of particle vector if necessary before joining next loop
  if(0.5 * s_backPtr->m_partList.capacity() < s_backPtr->a_noParticles()) {
    s_backPtr->m_partList.reserve(10 * s_backPtr->m_partList.capacity());
  }

  // 0) Compute break-up length:
  // (usually used when no initial droplet size distribution (IDS) is available and
  // the blob method is used
  // Break-up length can be calculated by Levich theory according to Levich, 1962

  const MFloat breakupLength =
      m_RTsecBreakUpLength ? m_Cbu * m_injectorNozzleDiameter * sqrt(material().ambientDensityRatio()) : 0;

  // if(droplet.m_breakUpTime > 0.43839) {
  // reset for particles which have been affected significantly by tumble motion!
  //  breakupLength = 0.0;
  //}


  ASSERT(m_injectorNozzleDiameter > 0, "ERROR: No valid nozzle diameter set.");
  ASSERT(breakupLength > -MFloatEps, "ERROR: Invalid breakup length.");

  // 1) Loop over all particles and check for break-up
  for(MInt i = 0; i < noPart; i++) {
    auto& droplet = s_backPtr->m_partList.at(i);

    if(droplet.firstStep()) continue;
    if(droplet.isInvalid()) continue;
    if(droplet.fullyEvaporated()) continue;
    ASSERT(!isnan(droplet.m_diameter), "Invalid diameter! ");
    if(droplet.m_diameter < s_backPtr->m_sizeLimit) continue;
    if(droplet.hadWallColl()) continue;

    ASSERT(droplet.m_cellId > 0, "ERROR: Invalid cellId");

    // increase breakup time
    droplet.m_breakUpTime += dt;

    ASSERT(!isnan(droplet.m_shedDiam) && droplet.m_shedDiam > -MFloatEps, "");

    // 2) Computations for Kelvin-Helmholtz Break up Model mainly based on:
    //   "Modeling Spray Atomization with the Kelvin-Helmholtz/Rayleigh-Taylor Hibrid Model"
    //    by Beale and Reitz, Atomization and Sprays (1999)
    // NOTE: in this case the radius and not the diameter is used as reference length!

    MFloat liquidLength = MFloatMax;
    if(m_RTsecBreakUpLength) {
      liquidLength = maia::math::distance(s_backPtr->m_spawnCoord, droplet.m_position);
      /*
      if(m_injectorType == MULTIHOLE) {
        const MFloat dist1 = maia::math::distance(s_backPtr->m_spawnCoord, droplet.m_position);
        const MFloat diff = POW2(dist1) - POW2(m_injectorDiameter);
        if(diff > 0) {
          liquidLength = sqrt(diff);
        } else {
          liquidLength = 0;
        }
      }
      */
    }

    const MFloat dropRadius = 0.5 * droplet.m_diameter;
    const MFloat relV = droplet.magRelVel(&droplet.m_fluidVel[0], &droplet.m_velocity[0]);

    //#ifdef LPT_DEBUG
    if(isnan(relV)) {
      cerr << droplet.m_fluidVel[0] << " " << droplet.m_fluidVel[1] << " " << droplet.m_fluidVel[nDim - 1] << endl;
      mTerm(1, AT_, "");
    }
    //#endif

    // radius based particle Weber-number
    const MFloat We_l =
        droplet.WeberNumber(s_backPtr->m_material->density(droplet.m_temperature), pow(relV, 2), droplet.m_temperature)
        * droplet.s_We;

    // if the weber number is to small the further calculated values are nan or inf
    if(We_l < m_weberLimitKH) continue;

    // particle Reynolds number
    // factor 0.5 to convert from diameter to radius based
    const MFloat Re_l = 0.5
                        * droplet.particleRe(material().density(droplet.m_temperature), relV,
                                             material().dynamicViscosity(droplet.m_temperature))
                        * droplet.s_Re;

    // gas Weber number
    const MFloat We_g = droplet.WeberNumber(droplet.m_fluidDensity, pow(relV, 2), droplet.m_temperature) * droplet.s_We;
    // Ohnesorge number
    const MFloat Z = sqrt(We_l) / Re_l;
    // Taylor number
    const MFloat T = Z * sqrt(We_g);

    // Lamda_KH
    const MFloat KH_waveL = 9.02 * dropRadius * (1.0 + 0.45 * sqrt(Z)) * (1.0 + 0.4 * pow(T, 0.7))
                            / pow(1.0 + 0.865 * pow(We_g, 1.67), 0.6);


    // Kelvin-Helmholtz reference diameter (thus factor 2.0)
    const MFloat KH_diameter = 2.0 * m_B0 * KH_waveL;

    // ASSERT(KH_diameter > MFloatEps,
    //       to_string(KH_diameter) + " surfWL " + to_string(KH_waveL));
    if(KH_diameter < 2 * MFloatEps) {
      cerr << KH_diameter << " d_KH " << dropRadius << " " << Z << " " << T << " " << We_g << " " << We_l << endl;
    }

    // Omega_KH
    const MFloat KH_growthR = ((0.34 + 0.38 * pow(We_g, 1.5)) / ((1.0 + Z) * (1.0 + 1.4 * pow(T, 0.6))))
                              * sqrt(material().spraySurfaceTension(droplet.m_temperature)
                                     / (material().density(droplet.m_temperature) * pow(dropRadius, 3)))
                              * sqrt(1 / droplet.s_We);


    // tau_kH
    const MFloat KH_time = 3.726 * m_B1 * dropRadius / (KH_waveL * KH_growthR);

    // const MFloat KH_speed = (droplet.m_diameter - KH_diameter) / KH_time;

    // 3) prepare initialisation of new droplets
    // parent velocity magnitude
    const MFloat magVel = droplet.magVel();
    array<MFloat, nDim> parentTrajectory{};
    // get parent droplet direction
    for(MInt v = 0; v < nDim; v++) {
      parentTrajectory.at(v) = droplet.m_velocity.at(v) / magVel;
    }

    // in case parent parcel is not moving
    if(magVel < MFloatEps) {
      // TODO labels:LPT,totest check if this exception is necessary!
      continue;
      /*
      // get direction from flow direction
      MFloat fluidVel = 0;
      for(MInt j = 0; j < nDim; j++) {
        fluidVel += POW2(s_backPtr->a_fluidVelocity(droplet.m_cellId, j));
      }
      fluidVel = sqrt(fluidVel);

      for(MInt v = 0; v < nDim; v++) {
        parentTrajectory[v] = s_backPtr->a_fluidVelocity(droplet.m_cellId, v) / fluidVel;
      }
      */
    }

    // magnitude of droplet acceleration
    MFloat drop_accel = 0.0;
    for(MInt dir = 0; dir < nDim; dir++) {
      drop_accel += POW2(droplet.m_accel.at(dir));
    }
    drop_accel = sqrt(drop_accel);

    // storing some old droplet properties
    const MFloat oldMass = droplet.sphericalMass() * droplet.m_noParticles;
    vector<MFloat> oldMom(nDim);
    for(MInt n = 0; n < nDim; n++) {
      oldMom[n] = oldMass * droplet.m_velocity[n];
    }
    const MFloat oldDiameter = droplet.m_diameter;

    // 4) Computations for Rayleigh-Taylor Break up Model mainly based on:
    //   "Modeling Spray Atomization with the Kelvin-Helmholtz/Rayleigh-Taylor Hibrid Model"
    //    by Beale and Reitz, Atomization and Sprays (1999)
    //    or Baumgarten

    // NOTE wavelength = 2 * PI * 1/ waveNumber(=K_RT)
    const MFloat RT_waveL =
        2.0 * PI
        * sqrt(3.0 * material().spraySurfaceTension(droplet.m_temperature)
               / (drop_accel * (material().density(droplet.m_temperature) - droplet.m_fluidDensity)))
        / sqrt(droplet.s_We);

    // Omega_RT
    const MFloat RT_frequency =
        sqrt(2.0 / (3.0 * sqrt(3.0 * material().spraySurfaceTension(droplet.m_temperature)))
             * pow(drop_accel * (material().density(droplet.m_temperature) - droplet.m_fluidDensity), 3.0 / 2.0)
             / (material().density(droplet.m_temperature) + droplet.m_fluidDensity))
        * pow(droplet.s_We, 1 / 4);

    const MFloat RT_time = m_CT / RT_frequency;

    // NOTE: in this case its the diameter!
    const MFloat RT_diameter = m_CRT * RT_waveL;

    ASSERT(RT_diameter > 0, "");

    // const MFloat RT_speed = (droplet.m_diameter - RT_diameter) / RT_time;


    // 6) Rayleigh-Taylor break-up
    if(m_RTsecBreakUp && liquidLength >= breakupLength && droplet.m_breakUpTime >= RT_time
       && droplet.m_diameter > RT_diameter && We_l < m_weberLimitRT) {
      // choose RT child diameter version
      MFloat childDiameter = RT_diameter;
      if(m_RTDiameterMode == 2) {
        initPRNG(globalTimeStep, s_backPtr->c_globalId(droplet.m_cellId));
        const MFloat diameterMin = RT_diameter / m_RosinRammlerMin;
        const MFloat diameterMean = RT_diameter / m_RosinRammlerMean;
        const MFloat diameterMax = RT_diameter / m_RosinRammlerMax;
        childDiameter = rosinRammler(diameterMin, diameterMean, diameterMax, m_RosinRammlerSpread, randomSecBreakUp());
        if(droplet.m_diameter < childDiameter) {
          // no breakup
          continue;
        }
      } else if(m_RTDiameterMode == 1) {
        childDiameter = pow(oldDiameter, 2 / 3) * pow(RT_diameter, 1 / 3);
      }

      if(m_maxRTChildParcels == 0) {
        // no child parcles are added, only the diameter and no-particles is changed
        auto noDroplets = static_cast<MInt>(droplet.m_noParticles * POW3(oldDiameter / childDiameter));
        // NOTE: casting obove means a runding down for all positive values
        //      version below applies additiona round up for values above 0.75

        // NOTE: applying rounding to avoid mass-losses!
        if(noDroplets > 0) {
          childDiameter = pow(6.0 / PI, 1.0 / 3.0)
                          * pow(oldMass / (noDroplets * material().density(droplet.m_temperature)), 1.0 / 3.0);
          droplet.m_diameter = childDiameter;
          droplet.m_shedDiam = childDiameter;
          droplet.m_noParticles = noDroplets;
          droplet.m_breakUpTime = 0.0;
          m_noRTsecBreakUp++;

          const MFloat newMass = sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles;
          if(fabs(newMass - oldMass) > MFloatEps) {
            cerr << "Missing mass is RT-breakup " << oldMass << " " << newMass << endl;
          }
        }

      } else {
        if(m_RTDiameterMode != 2) {
          initPRNG(globalTimeStep, s_backPtr->c_globalId(droplet.m_cellId));
        }

        // calculate new number of droplets
        const MFloat dropletVolume = pow(droplet.m_diameter, 3.0);
        const MFloat dropletRTVolume = pow(childDiameter, 3.0);

        MInt noBrokenUpDroplets = dropletVolume / dropletRTVolume;

        // don't allow inconsequential breakup
        if(noBrokenUpDroplets == 1) {
          continue;
        }
        noBrokenUpDroplets *= droplet.m_noParticles;

        // number of drops per new parcel
        auto dropsPerParcel = static_cast<MInt>(
            noBrokenUpDroplets > m_maxRTChildParcels ? floor(noBrokenUpDroplets / m_maxRTChildParcels) : 1);

        // limit the number of new parcels
        MInt noNewParcels = noBrokenUpDroplets;
        if(dropsPerParcel > 1) {
          noNewParcels = m_maxRTChildParcels - 1;
        } else {
          --noNewParcels;
        }

        const MFloat RT_dropletMass = sphericalMass(childDiameter, droplet.m_temperature);
        const MFloat parentMass = oldMass - dropsPerParcel * RT_dropletMass * noNewParcels;
        const MFloat parentDiam =
            pow(6.0 / PI, 1.0 / 3.0)
            * pow(parentMass / (dropsPerParcel * material().density(droplet.m_temperature)), 1.0 / 3.0);

        ASSERT(parentMass > 0, "Invalid mass!");
        ASSERT(parentDiam > 0 && !isnan(parentDiam), "");
        ASSERT(dropsPerParcel > 0, "");

        // Reset values
        droplet.m_diameter = parentDiam;
        droplet.m_shedDiam = parentDiam;
        droplet.m_noParticles = dropsPerParcel;
        droplet.m_breakUpTime = 0.0;
        m_noRTsecBreakUp++;

#ifdef _OPENMP
#pragma omp critical
#endif
        if(m_secBUDisplayMessage) {
          cerr << "#################################################" << endl;
          cerr << "Break-Up type: RT" << endl;
          cerr << "parent diameter: " << droplet.m_diameter << endl;
          cerr << "Time: " << m_timeSinceSOI << endl;
          cerr << "old number of particles: " << droplet.m_noParticles << endl;
          cerr << "new particles/parcels: " << noNewParcels << endl;
          cerr << "parcel size: " << dropsPerParcel << endl;
          cerr << "Weber number: " << We_l << endl;
          cerr << "#################################################" << endl;
        }

        MFloat sprayAngle = m_sprayAngle[0];
        if(m_hollowConeInjector) {
          // hollow cone injectors have a large spray angle because of the hollow cone
          // geometry the actual relevant spray angle is much smaller
          static constexpr MFloat angularGapAngle = 45.0;
          sprayAngle = m_sprayAngle[0] - 2 * angularGapAngle;
        }

#ifdef LPT_DEBUG
        vector<MFloat> mom(nDim);
#endif

        // create new droplets
        array<MFloat, nDim> childVelocity{};
        for(MInt s = 0; s < noNewParcels; s++) {
          randomVectorInCone(&childVelocity[0], &parentTrajectory[0], magVel, sprayAngle, m_partEmittDist,
                             randomSecBreakUp(), m_injectionDistSigmaCoeff);

          const MInt childId = s_backPtr->addParticle(droplet.m_cellId, childDiameter, droplet.m_densityRatio, 0, 3,
                                                      &childVelocity[0], &droplet.m_position[0], dropsPerParcel);

          // set child temperature based on parent value
          s_backPtr->m_partList[childId].m_temperature = droplet.m_temperature;
          s_backPtr->m_partList[childId].m_heatFlux = 0.0;
          s_backPtr->m_partList[childId].m_dM = 0.0;


#ifdef LPT_DEBUG
          // check kin. energy conservation
          MFloat magVelChild = 0;
          for(MInt n = 0; n < nDim; n++) {
            magVelChild += POW2(childVelocity[n]);
          }
          magVelChild = sqrt(magVelChild);
          if(fabs(magVel - magVelChild) > MFloatEps) {
            cerr << "RT kin. energy loss: " << magVel << " " << magVelChild << endl;
          }
          for(MInt n = 0; n < nDim; n++) {
            mom[n] += sphericalMass(childDiameter, droplet.m_temperature) * dropsPerParcel * childVelocity[n];
          }
#endif
        }

#ifdef LPT_DEBUG
        // check momentum conservation:
        for(MInt n = 0; n < nDim; n++) {
          mom[n] +=
              droplet.m_velocity[n] * sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles;

          if(fabs(oldMom[n] - mom[n]) > MFloatEps) {
            cerr << "RT momentum change: " << oldMom[n] << " " << mom[n] << endl;
          }
        }

        // check conservation:
        const MFloat newMass = sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles;
        MFloat childMass = 0;
        for(MInt n = 0; n < noNewParcels; n++) {
          const MInt id = s_backPtr->a_noParticles() - 1 - n;
          childMass += (sphericalMass(s_backPtr->m_partList[id].m_diameter, s_backPtr->m_partList[id].m_temperature)
                        * s_backPtr->m_partList[id].m_noParticles);
        }
        const MFloat massLoss_RT = oldMass - newMass - childMass;
        if(fabs(massLoss_RT) > MFloatEps) {
          cerr << "RT Mass loss! Old " << oldMass << " parent " << newMass << " childs " << childMass << " diff "
               << massLoss_RT << endl;
        }
#endif
      }
      continue;
    }

    // Kelvin-Helmholtz induced break-up
    // if((KH_speed >= RT_speed || liquidLength < breakupLength) && droplet.m_diameter >= KH_diameter && m_KHsecBreakUp)
    // {
    if(m_KHsecBreakUp && droplet.m_diameter > KH_diameter) {
      const MFloat KH_mass = sphericalMass(KH_diameter, droplet.m_temperature);
      if(KH_mass < MFloatEps) continue;
      // ASSERT(KH_mass > MFloatEps, to_string(KH_mass) + " KH_diameter:" + to_string(KH_diameter));

      // Kelvin-Helmholtz induced mass stripping
      // NOTE: shedDiam is the diameter of the parent particle after undergoing KH-breakup!
      // i.e. the remaining mass in the parent!
      droplet.m_shedDiam -= (droplet.m_shedDiam - KH_diameter) / KH_time * dt;

      // if the shed diameter is smaller than the KH diameter
      // no child parcel is added, the parent diameter is set to the KH-diameter
      // and the noParticles is set by mass conservation
      if(droplet.m_shedDiam < KH_diameter) {
        MInt noDroplets = floor(oldMass / KH_mass);
        if(noDroplets > 0) {
          const MFloat diameter = pow(6.0 / PI, 1.0 / 3.0)
                                  * pow(oldMass / (noDroplets * material().density(droplet.m_temperature)), 1.0 / 3.0);
          // fully transversion to KH-diameter
          // after applying rounding to avoid mass-losses!
          // droplet.m_breakUpTime = -dt;
          droplet.m_diameter = diameter;
          droplet.m_shedDiam = diameter;
          droplet.m_noParticles = noDroplets;
          m_noKHsecBreakUp++;

          const MFloat newMass = sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles;
          if(fabs(newMass - oldMass) > MFloatEps) {
            cerr << "Missing mass at KH-breakup 1 is " << oldMass << " " << newMass << endl;
          }
          continue;
        }
      }

      // mass going into the particle = oldMass - mass remaining in the parant parcel
      const MFloat sheddedMass =
          oldMass - sphericalMass(droplet.m_shedDiam, droplet.m_temperature) * droplet.m_noParticles;

      if(m_massLimit > sheddedMass / oldMass) {
        continue;
      }

      MInt noSheddedOffDrops = floor(sheddedMass / KH_mass);
      if(noSheddedOffDrops <= 0) continue;
      const MFloat parentMass = oldMass - KH_mass * noSheddedOffDrops;
      ASSERT(parentMass > 0, "");
      // when reformulating the equations above, parentDiam is equal to the shedDiam!
      const MFloat parentDiam =
          pow(6.0 / PI * parentMass / (droplet.m_noParticles * material().density(droplet.m_temperature)), 1.0 / 3.0);

      // adding additional child parcel
      // child parcel properties set as prescribed in
      //"Modeling Atomization Processes in High-Pressure Vaporizing Sprays"
      // R. Reitz , Atom.& Sprayz 3 (1987) 309-337
      initPRNG(globalTimeStep, s_backPtr->c_globalId(droplet.m_cellId));

      // droplet.m_breakUpTime = -dt;
      droplet.m_diameter = parentDiam;
      droplet.m_shedDiam = parentDiam;

      // different spray-angle versions
      MFloat angle = m_sprayAngle[0];
      if(m_sprayAngleKH > 1 && m_sprayAngleKH < 2.9) {
        angle = m_sprayAngle[0];
        if(liquidLength < breakupLength) {
          angle = m_sprayConeAngle;
        }
      } else if(m_sprayAngleKH > 0 && m_sprayAngleKH < 1) {
        angle = 2 * atan(m_sprayAngleKH * KH_growthR * KH_waveL / m_injectionSpeed);
        angle = min(m_sprayAngle[0], angle);
      } else if(m_sprayAngleKH > 2.9) {
        if(liquidLength < breakupLength) {
          angle = m_sprayAngle[0];
        } else {
          angle = m_sprayAngle[1];
        }
      }

      if(angle < 0 || angle > 90) {
        cerr << "Large KH-spray angle : " << angle << endl;
      }

      array<MFloat, nDim> childVelocity{};
      randomVectorInCone(&childVelocity[0], &parentTrajectory[0], magVel, angle, string2enum("PART_EMITT_DIST_NONE"),
                         randomSecBreakUp(), m_injectionDistSigmaCoeff);

      // adding a single additional parcel with KH_diameter
      const MInt childId = s_backPtr->addParticle(droplet.m_cellId, KH_diameter, droplet.m_densityRatio, 0, 3,
                                                  &childVelocity[0], &droplet.m_position[0], noSheddedOffDrops);

      // set child temperature based on parent value
      s_backPtr->m_partList[childId].m_temperature = droplet.m_temperature;
      s_backPtr->m_partList[childId].m_heatFlux = 0.0;
      s_backPtr->m_partList[childId].m_dM = 0.0;


      m_noKHsecBreakUp++;

#ifdef _OPENMP
#pragma omp critical
#endif
      if(m_secBUDisplayMessage) {
        cerr << "#################################################" << endl;
        cerr << "Break-Up type: KH" << endl;
        cerr << "parent diameter: " << droplet.m_diameter << " " << parentDiam << endl;
        cerr << "KH size: " << KH_diameter << endl;
        cerr << "KH time: " << KH_time << endl;
        cerr << "Time: " << m_timeSinceSOI << endl;
        cerr << "old #droplets: " << droplet.m_noParticles << endl;
        cerr << "# new droplets: " << noSheddedOffDrops << endl;
        cerr << "Weber number: " << We_l << endl;
        cerr << "#################################################" << endl;
      }

      // apply velocity change to parent particle
      // under momentum and kin. energy conservation!
      /*
      vector<MFloat> childMom(nDim);
      for(MInt n = 0; n < nDim; n++){
        childMom[n] = sphericalMass(KH_diameter) * noSheddedOffDrops * childVelocity[n];
        droplet.m_velocity[n] = (oldMom[n] - childMom[n])
                              / (sphericalMass(droplet.m_diameter) * droplet.m_noParticles);
      }

      maia::math::normalize(&droplet.m_velocity[0], 3);
      for(MInt n = 0; i < nDim; n++){
        droplet.m_velocity[i] *= magVel;
      }
      */

#ifdef LPT_DEBUG
      // check momentum conservation
      vector<MFloat> mom(nDim);
      for(MInt n = 0; n < nDim; n++) {
        mom[n] =
            droplet.m_velocity[n] * sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles
            + sphericalMass(KH_diameter, droplet.m_temperature) * noSheddedOffDrops * childVelocity[n];
        if(fabs(oldMom[n] - mom[n]) > MFloatEps) {
          cerr << "KH momentum change: " << oldMom[n] << " " << mom[n] << " "
               << droplet.m_velocity[n] * sphericalMass(droplet.m_diameter, droplet.m_temperature)
                      * droplet.m_noParticles
               << endl;
        }
      }

      // check kin. energy conservation
      MFloat magVelChild = 0;
      for(MInt n = 0; n < nDim; n++) {
        magVelChild += POW2(childVelocity[n]);
      }
      magVelChild = sqrt(magVelChild);
      if(fabs(magVel - magVelChild) > MFloatEps) {
        cerr << "KH kin. energy loss: " << magVel << " " << magVelChild << endl;
      }
      MFloat magVelNew = droplet.magVel();
      if(fabs(magVel - magVelNew) > MFloatEps) {
        cerr << "KH kin. energy loss: " << magVel << " " << magVelNew << endl;
      }
      // check mass conservation
      const MFloat newMass = sphericalMass(droplet.m_diameter, droplet.m_temperature) * droplet.m_noParticles;
      const MInt childId = s_backPtr->a_noParticles() - 1;
      const MFloat childMass =
          sphericalMass(s_backPtr->m_partList[childId].m_diameter, s_backPtr->m_partList[childId].m_temperature)
          * s_backPtr->m_partList[childId].m_noParticles;
      const MFloat massLoss_KH = oldMass - newMass - childMass;
      if(fabs(massLoss_KH) > MFloatEps) {
        cerr << "KH Mass loss! Old " << oldMass << " parent " << newMass << " childs " << childMass << " diff "
             << massLoss_KH << endl;
      }
#endif
    } /*else if( m_KHsecBreakUp && KH_diameter > droplet.m_diameter &&
               droplet.m_breakUpTime > KH_time ) {
      //version by Reitz (1987)
      MFloat d_KH = 2 * mMin(
         pow(3 * PI * POW2(droplet.m_diameter) * relV / ( 2 * 4 * KH_growthR), 0.33),
         pow(3 * POW2(droplet.m_diameter) * KH_waveL / 16, 0.33));

      auto noDroplets =  static_cast<MInt>(droplet.m_noParticles *
                                           POW3(droplet.m_diameter) / POW3(d_KH));
      if(noDroplets > 0 ) {
        const MFloat diameter = pow(6.0 / PI, 1.0 / 3.0) *
                                pow(oldMass / (noDroplets * material().density(droplet.m_temperature)), 1.0/ 3.0);
        //change droplet diameter
        droplet.m_breakUpTime = 0.0;
        droplet.m_diameter = diameter;
        droplet.m_shedDiam = diameter;
        droplet.m_noParticles = noDroplets;
        m_noKHsecBreakUp++;
      }
    }
    */
  }
}

/// \brief Update injection rate
///
/// \author Piotr Duda, Sven Berger
/// \date   June 2016
template <MInt nDim>
void SprayModel<nDim>::updateInjectionRate() {
  const MInt numInjections = m_injectionRateList.dim0();
  if(numInjections == 0) return;

  if(m_injectionTimings[0] > m_timeSinceSOI) {
    m_currentInjectionRate = (m_timeSinceSOI + 1) / (m_injectionTimings[0] + 1) * m_injectionRateList[0];
  } else {
    for(MInt i = 0; i <= numInjections; i++) {
      if(m_injectionTimings[i] > m_timeSinceSOI) {
        const MFloat t =
            (m_timeSinceSOI - m_injectionTimings[i - 1]) / (m_injectionTimings[i] - m_injectionTimings[i - 1]);
        m_currentInjectionRate = m_injectionRateList[i - 1] + t * (m_injectionRateList[i] - m_injectionRateList[i - 1]);
        break;
      }
    }
  }
}


/// \brief Create new random particles using the provided options.
///
/// \author Sven Berger, Tim Wegmann
/// \date   November 2015
/// \param[in] spawnParticlesCount Number of particles created in the given location.
/// \param[in] spawnCoord Coordinates of the location to spawn the particles in. (Note: Check if
/// this is on the calling procs domain!!!)
/// \param[in] particleDiameter Size of the particles that
/// get created.
/// \param[in] particleDensityRatio Density ratio of the newly created particles.
/// \param[in] particleVelocity initial velocity of the created particles.
/// \param[in] sprayConeAngle Angle of the cone in which particles are created
/// \param[in] coneDiameter intial diameter of the spray cone (0.0 for blob-injection)
/// \param[in] hull Only create particles on the hull of the cone
/// \param[in] holeNozzleAngle Opening angle of the nozzle
/// \param[in] parcelSize Parcel size that is used for the newly created particles
/// \param[in] noInjectionHoles Number of holes of the injector
/// \param[in] injectorDiameter Diameter of the injector.
template <MInt nDim>
void SprayModel<nDim>::injectParticles(const MInt spawnParticlesCount, const MFloat* spawnCoord,
                                       const MFloat particleDiameter, const MFloat particleDensityRatio,
                                       const MFloat particleVelocity, const MFloat sprayConeAngle,
                                       const MFloat coneDiameter, std::function<MFloat(MInt)> holePositionAngle,
                                       std::function<MFloat(MInt)> injectorDiameter,
                                       std::function<MFloat(MInt)> holeNozzleAngle, const MBool hull,
                                       const MInt parcelSize, const MInt noInjectionHoles) {
  std::array<MFloat, nDim> spawnParticlesInitVelo;
  array<MFloat, 3> injectionLocation{};
  array<MFloat, 3> holeLocation{};

  if(spawnParticlesCount > 0) {
    m_injStep++;
  }

  if(particleDiameter < s_backPtr->m_sizeLimit) {
    cerr << "Inj. small particle " << particleDiameter << " " << s_backPtr->m_sizeLimit << endl;
  }

  for(MInt i = 0; i < spawnParticlesCount; i++) {
    MInt randomShift = 1;
    if(m_partEmittDist != PART_EMITT_DIST_NONE) {
      randomShift = 0;
      m_PRNGPrimBreakUpCount +=
          randomVectorInCone(&spawnParticlesInitVelo[0], &m_injectorDir[0], particleVelocity, sprayConeAngle,
                             m_partEmittDist, randomPrimBreakUp(0), m_injectionDistSigmaCoeff, 0);

      if(coneDiameter > numeric_limits<MFloat>::epsilon() && !hull) {
        // determine a random point within the injector orfice
        randomPointInCircle(&holeLocation[0], &m_injectorDir[0], coneDiameter, randomPrimBreakUp(2));

        if(noInjectionHoles > 1) {
          // determine position of the hole
          pointOnCircle(&injectionLocation[0], &m_injectorDir[0], injectorDiameter(i), holePositionAngle(i));
        }
      } else if(hull) {
        m_PRNGPrimBreakUpCount +=
            randomVectorInCone(&spawnParticlesInitVelo[0], &m_injectorDir[0], particleVelocity, sprayConeAngle,
                               m_partEmittDist, randomPrimBreakUp(0), m_injectionDistSigmaCoeff, holeNozzleAngle(i));

        // determine a point just on the hull of the cone (like for a hollow-cone injector)
        // TODO labels:LPT make settable
        static constexpr MBool randomDist = false;
        if(randomDist) {
          randomPointOnCircle(&injectionLocation[0], &m_injectorDir[0], coneDiameter, randomPrimBreakUp(1),
                              spawnParticlesCount, i);
        } else {
          // generate one particle per 1 deg
          const static MFloat angleBetween = 360.0 / spawnParticlesCount;
          const MFloat phi = (i * angleBetween + ((m_injStep - 1) % static_cast<MInt>(angleBetween))) / 180.0 * M_PI;
          pointOnCircle(&injectionLocation[0], &m_injectorDir[0], coneDiameter, phi);
        }
      }
    }

    // rotate injection by the angle of the nozzle
    if(fabs(holeNozzleAngle(i)) > numeric_limits<MFloat>::epsilon()) {
      const MFloat holeAngleRad = -holeNozzleAngle(i) / 180 * M_PI;
      const MFloat absDist = maia::math::norm(injectionLocation);

      array<MFloat, nDim> crossP{};
      array<MFloat, nDim> revInjLoc{};
      for(MInt j = 0; j < nDim; j++) {
        revInjLoc[j] = -injectionLocation[j] / absDist;
      }

      maia::math::cross(&m_injectorDir[0], &revInjLoc[0], &crossP[0]);

      // boundary basis at the injection location
      MFloatScratchSpace injLocBasis(nDim, nDim, FUN_, "R");
      for(MInt j = 0; j < nDim; j++) {
        injLocBasis(j, 0) = m_injectorDir[j];
        injLocBasis(j, 1) = revInjLoc[j];
        injLocBasis(j, 2) = crossP[j];
      }

      // find inverse for reverse-tranformation
      MFloatScratchSpace inverse(nDim, nDim, AT_, "inverse");
      for(MInt j = 0; j < nDim; j++) {
        for(MInt k = 0; k < nDim; k++) {
          inverse(j, k) = injLocBasis(j, k);
        }
      }

      maia::math::invert(&inverse(0, 0), 3, 3);
      // NOTE: inverse might not be normalized correctly, meaning that the magnitude is not conserved!

      // rotation matrix around the z-axis in the injection location basis
      MFloatScratchSpace R(nDim, nDim, FUN_, "R");
      R.fill(0.0);
      R(2, 2) = 1.0;
      R(0, 0) = cos(holeAngleRad);
      R(1, 0) = sin(holeAngleRad);
      R(0, 1) = -sin(holeAngleRad);
      R(1, 1) = cos(holeAngleRad);

      // global coordinate transformation is then found be B*R*B^-1
      MFloatScratchSpace result1(nDim, nDim, FUN_, "result1");
      MFloatScratchSpace result(nDim, nDim, FUN_, "result");
      MFloatScratchSpace tempV(nDim, FUN_, "tempV");
      for(MInt j = 0; j < nDim; j++) {
        tempV[j] = spawnParticlesInitVelo[j];
      }

      maia::math::multiplyMatricesSq(injLocBasis, R, result1, nDim);
      maia::math::multiplyMatricesSq(result1, inverse, result, nDim);
      for(MInt j = 0; j < nDim; j++) {
        spawnParticlesInitVelo[j] = 0;
        for(MInt k = 0; k < nDim; k++) {
          spawnParticlesInitVelo[j] += result(j, k) * tempV[k];
        }
      }

      // reset vector length and direction
      maia::math::normalize(&spawnParticlesInitVelo[0], nDim);
      // reverse sign if velocity has opposite direction to the injector direction
      if(scalarProduct(&spawnParticlesInitVelo[0], &m_injectorDir[0], nDim) < 0.0) {
        for(MInt j = 0; j < nDim; j++) {
          spawnParticlesInitVelo[j] *= -1;
        }
      }
      for(MInt j = 0; j < nDim; j++) {
        spawnParticlesInitVelo[j] *= particleVelocity;
      }
    }
    // final injection position = hole position + orifice position + injector center
    for(MInt j = 0; j < nDim; j++) {
      injectionLocation[j] += spawnCoord[j] + holeLocation[j];
    }

    s_backPtr->addParticle(s_backPtr->m_spawnCellId, particleDiameter, particleDensityRatio, randomShift, 2,
                           &spawnParticlesInitVelo[0], &injectionLocation[0], parcelSize);
  }
}

/**
 *  \brief create log of distribution functions for plots and debugging
 *  \author Tim Wegmann
 *  \date   January 2023
 */
template <MInt nDim>
void SprayModel<nDim>::logDistributionFunction() {
  if(domainId() != 0) return;

  const MInt noDraws = 1000000;
  MFloat minValue = m_injectorNozzleDiameter / m_RosinRammlerMin;
  MFloat maxValue = m_injectorNozzleDiameter / m_RosinRammlerMax;
  MFloat meanValue = m_injectorNozzleDiameter / m_RosinRammlerMean;

  std::mt19937_64 randomNumberGenerator;
  randomNumberGenerator.seed(m_spraySeed);

  const MInt noBinns = 1000;
  MFloat binWidth = (maxValue - minValue) / noBinns;
  MIntScratchSpace binValue(noBinns, AT_, "binValue");
  binValue.fill(0);

  for(MInt i = 0; i < noDraws; i++) {
    const MFloat drawValue = rosinRammler(minValue, meanValue, maxValue, m_RosinRammlerSpread, randomNumberGenerator);
    const MInt binId = std::floor((drawValue - minValue) / binWidth);
    binValue[binId]++;
  }

  m_log << "IDSD Binning: " << endl;
  m_log << "Rosin-Rammler distribution: " << endl;
  m_log << "Min-value : " << minValue << endl;
  m_log << "Max-value : " << maxValue << endl;
  m_log << "Mean-value: " << meanValue << endl;
  m_log << "Spread    : " << m_RosinRammlerSpread << endl;
  m_log << "Number-draws: " << noDraws << endl;

  for(MInt i = 0; i < noBinns; i++) {
    m_log << minValue + i * binWidth << " " << binValue[i] << endl;
  }

  minValue = 0.0;
  maxValue = 5 * meanValue;
  binWidth = (maxValue - minValue) / noBinns;

  binValue.fill(0);

  for(MInt i = 0; i < noDraws; i++) {
    const MFloat drawValue = NTDistribution(meanValue, randomNumberGenerator);
    const MInt binId = std::floor((drawValue - minValue) / binWidth);
    binValue[binId]++;
  }

  m_log << "Nukiyama-Tanasawa distribution: " << endl;
  m_log << "Mean-value: " << meanValue << endl;
  for(MInt i = 0; i < noBinns; i++) {
    m_log << minValue + i * binWidth << " " << binValue[i] << endl;
  }
}
// Explicit instantiations for 2D and 3D
template class SprayModel<3>;
