// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_SPRAY_H
#define MAIA_SPRAY_H

#include <random>
#include "UTIL/maiamath.h"
#include "UTIL/materialstate.h"
#include "UTIL/tensor.h"
#include "globals.h"
#include "lptspherical.h"

template <MInt nDim>
class LPTSpherical;

template <MInt nDim>
class SprayModel {

 public:
  SprayModel(){};

  void init(LPT<nDim>* _particle);

  void injection(const MFloat dt);

  void secondaryBreakUp(const MFloat dt);

  MFloat timeSinceSOI() const { return m_timeSinceSOI; }
  MFloat& timeSinceSOI() { return m_timeSinceSOI; }

  MFloat injectionDuration() const { return m_injDuration; }
  MFloat injectionSpeed() const { return m_injectionSpeed; }

  static constexpr MInt m_injDataSize = 11;
  MFloat* m_injData = nullptr;

  MBool m_broadcastInjected = true;

  MInt m_injStep = 0;
  MInt m_injStopTimeStep = -1;
  MInt m_injStartTimeStep = -1;

  MBool soonInjection(const MFloat time) {
    if(m_injStartTimeStep > -1 && globalTimeStep + 100 < m_injStartTimeStep) {
      return true;
    }
    if(m_injectionCrankAngle < 0) {
      return true;
    } else if(m_injectionCrankAngle > -1
              && maia::math::crankAngle(time, m_Strouhal, m_initialCad, 0) > m_injectionCrankAngle) {
      return true;
    }
    return false;
  }

  MFloat m_injectionCrankAngle = -1;
  MFloat m_Strouhal = -99;
  MFloat m_initialCad = 0;

 private:
  // Point back to particle container
  static LPT<nDim>* s_backPtr;
  static MFloat s_lengthFactor;

  MInt domainId() const { return s_backPtr->domainId(); }
  MInt solverId() const { return s_backPtr->solverId(); }

  MaterialState<nDim>& material() const { return *s_backPtr->m_material; }

  MFloat sphericalMass(const MFloat diameter, const MFloat temperature) {
    return 1.0 / 6.0 * PI * material().density(temperature) * POW3(diameter);
  }

  void readSprayProperties();
  void logDistributionFunction();

  // ---------- properties which must be matching the used non-dimensionalisation ------------

  // current injection rate [kg/sec] ( rho_ref * L_ref^2 * u_ref )
  MFloat m_currentInjectionRate = 0.0;

  // list of injections rates [kg/sec] rho_ref * L_ref^2 * u_ref
  MFloatTensor m_injectionRateList;

  // list of injection timings [sec]  L_ref/u_ref
  MFloatTensor m_injectionTimings;

  // minimal size of droplets during primary breakup [m] / L_ref
  MFloat m_primMinDropletSize = 5e-7;

  // time until injector closes [sec]  L_ref/u_ref
  MFloat m_injDuration = 0.00055; // 0.00055 20MPa 0.000625 10MPa

  // Injector design parameter (diameter at the end of the injector nozzle)
  // The size of the Blobs in the Blob-method is identical to this value! [m] / L_ref
  MFloat m_injectorNozzleDiameter = 0.0;

  // injection speed/velocity [m/sec] / u_ref
  MFloat m_injectionSpeed = -1.0;

  // diameter of multi-hole injector [m] / L_ref
  MFloat* m_injectorDiameter = nullptr;

  // orifice angle of multi-hole injector [degree]
  MFloat* m_injectorOrificeAngle = nullptr;

  // diameter of multi-hole injector [m] / L_ref
  MFloat* m_orificePositionAngle = nullptr;

  // diameter of multi-hole injector orfice  [m] / L_ref
  MFloat m_injectorOrficeSize{};

  // min diameter for Rosin-Rammler distribution [m] / L_ref
  MFloat m_RosinRammlerMin = -1;

  // mean diameter for Rosin-Rammler distribution [m] / L_ref
  MFloat m_RosinRammlerMean = -1;

  // max diameter for Rosin-Rammler distribution [m] / L_ref
  MFloat m_RosinRammlerMax = -1;

  //--------------------------------------------------------------------------------------------

  // spread for Rosin-Rammler distribution [-]
  MFloat m_RosinRammlerSpread = -1;

  // time until the injector is fully open [ percentage of the injection duration]
  MFloat* m_needleLiftTime = nullptr;

  // injection direction
  MFloat m_injectorDir[nDim]{};

  MFloat m_injectionDistSigmaCoeff = 1.0;

  // using string2enum
  MInt m_partEmittDist = 0;

  // cone angle of the spray
  MFloat* m_sprayAngle = nullptr;
  MFloat m_sprayConeAngle = 0.0;

  // Properties for primary breakup Rosin-Rammler distribution
  MBool m_primRosinRammler = true;

  MFloat m_timeSinceSOI = 0.0;

  // Model properties of secondary break-up
  MFloat m_B0 = 0.61;
  MFloat m_B1 = 40;
  MFloat m_CRT = 0.1;
  MFloat m_CT = 1.0;
  MFloat m_weberLimitKH = 6;
  MFloat m_massLimit = 0.03;
  MFloat m_weberLimitRT = 300;
  MBool m_secBUDisplayMessage = false;
  MInt m_maxRTChildParcels = 3;
  MFloat m_Cbu = 20;
  MInt m_RTDiameterMode = 0;
  MFloat m_sprayAngleKH = 2;

  MBool m_activePrimaryBUp = true;
  MBool m_activeSecondaryBUp = false;

  MInt m_primBrkUpParcelSize = 1;

  MBool m_useNeedleLiftTime = false;

  MFloat m_angularGap = -1;

  // using string2enum!
  MInt m_injectorType = 0;
  MBool m_multiHoleInjector = false;
  MBool m_singleHoleInjector = false;
  MBool m_hollowConeInjector = false;

  MBool m_predictivePRNG = true;
  MBool m_RTsecBreakUp = true;
  MBool m_KHsecBreakUp = true;
  MBool m_RTsecBreakUpLength = true;

  MPI_Comm m_primBUp;

  MInt m_maxNoPrimParcels = 1;
  MInt m_primParcelsPerInj = 1000000;

  // random number for secondary breakup
  std::mt19937_64 m_PRNGSecBU;

  void updateInjectionRate();

  void primaryBreakUp(const MFloat dt);

  void injectParticles(const MInt spawnParticlesCount, const MFloat* spawnCoord, const MFloat particleDiameter,
                       const MFloat particleDensityRatio, const MFloat particleVelocity, const MFloat sprayConeAngle,
                       const MFloat coneDiameter, std::function<MFloat(MInt)> holePositionAngle,
                       std::function<MFloat(MInt)> injectorDiameter, std::function<MFloat(MInt)> holeNozzleAngle,
                       const MBool hull, const MInt parcelSize, const MInt noInjectionHoles);

  void injectParticles(const MInt spawnParticlesCount, const MFloat* spawnCoord, const MFloat particleDiameter,
                       const MFloat particleDensityRatio, const MFloat particleVelocity, const MFloat sprayConeAngle,
                       const MFloat coneDiameter, const MBool hull, const MFloat nozzleAngle, const MInt parcelSize,
                       const MInt noInjectionHoles) {
    std::function<MFloat(MInt)> dummyDefault = [&](MInt) { return -1.0; };
    std::function<MFloat(MInt)> holeNozzleAngle = [&](MInt) { return nozzleAngle; };
    ASSERT(noInjectionHoles == 1, "Incorrect function call!");

    injectParticles(spawnParticlesCount, spawnCoord, particleDiameter, particleDensityRatio, particleVelocity,
                    sprayConeAngle, coneDiameter, dummyDefault, dummyDefault, holeNozzleAngle, hull, parcelSize,
                    noInjectionHoles);
  }


  // randon number generator for the secondary breakup
  // set seed by particle-position and discard by globalTimeStep!
  // only reset if the PRNG is required to be very predictable!
  void initPRNG(const MInt seed, const MInt discard) {
    if(m_predictivePRNG) {
      m_PRNGSecBU.seed(seed);
      m_PRNGSecBU.discard(discard);
    }
  }

  inline std::mt19937_64& randomSecBreakUp() { return m_PRNGSecBU; }

 public:
  // randon number for primary breakup at spwanDomainId
  std::mt19937_64 m_PRNGPrimBreakUp;
  MInt m_PRNGPrimBreakUpCount = 0;

  inline std::mt19937_64& randomPrimBreakUp(const MInt calls) {
    ASSERT(domainId() == s_backPtr->m_spawnDomainId, "");
    m_PRNGPrimBreakUpCount += calls;
    return m_PRNGPrimBreakUp;
  }

  MLong m_spraySeed{};

  MInt m_noRTsecBreakUp = 0;
  MInt m_noKHsecBreakUp = 0;
};

#endif // MAIA_SPRAY_H
