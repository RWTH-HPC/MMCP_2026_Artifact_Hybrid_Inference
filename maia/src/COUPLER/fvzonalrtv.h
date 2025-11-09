// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVZONALRTV_H_
#define FVZONALRTV_H_

#include <vector>
#include "FV/fvcartesiansyseqnns.h"
#include "FV/fvcartesiansyseqnrans.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "UTIL/functions.h"
#include "coupling.h"
#include "fvzonal.h"

template <MInt nDim, class SysEqn>
class FvZonal;

/*author Jannik: December 2019
Coupling class for zonal RANS-LES method
*/
template <MInt nDim, class SysEqn>
class FvZonalRTV : public FvZonal<nDim, SysEqn> {
 private:
  friend class FvZonal<nDim, SysEqn>;

 public:
  using RANS = FvCartesianSolverXD<nDim, SysEqn>;
  using LES = FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>;

  // Constructor
  FvZonalRTV(const MInt couplingId, RANS* R, LES* L);

  using Base = FvZonal<nDim, SysEqn>;
  using Base::a_noFvCellsLES;
  using Base::a_noFvCellsRANS;
  using Base::a_noFvGridCellsLES;
  using Base::a_noFvGridCellsRANS;
  using Base::getAveragingFactor;
  using Base::initLESValues;
  using Base::initRANSValues;
  using Base::LESSolver;
  using Base::m_LESNoVarAverage;
  using Base::m_LESSolverId;
  using Base::m_noReconstructNutVars;
  using Base::m_RANSSolverId;
  using Base::m_restartLESAverage;
  using Base::m_zonalAveragingTimeStep;
  using Base::m_zonalTransferInterval;
  using Base::noExchangeVariables;
  using Base::noLESVariables;
  using Base::noRANSVariables;
  using Base::RANSSolver;

  void init();
  void finalizeCouplerInit();
  void preCouple(MInt);
  void postAdaptation();

 private:
  // suppressed functions
  void postCouple(MInt){};
  void finalizeSubCoupleInit(MInt){};
  void finalizeAdaptation(const MInt){};
  void prepareAdaptation(){};
  void subCouple(MInt, MInt, std::vector<MBool>&){};
  void cleanUp(){};
  void initData(){};
  void checkProperties(){};
  void readProperties(){};

  void determineNutReconstructionCells();
  void determineZonalPositions();
  void reconstructNutilde();
  void transferSolverData();
  void finalizeRANSSolverAfterReset();
  void reconstructAverageFromNut();

  MInt m_rntStartTimeStep;
  MBool m_reconstructNut = false;
  MBool m_reconstructAverageFromNut = false;

  MFloat m_nuTildeInfinity;
  MInt m_averageTimeSteps;

  MFloat m_turbulentIntensity;
  MFloat m_tuLengthScaleCorrection;

  std::vector<MInt> m_rntBcCells;

  MInt m_noRntBcCells;

  MInt m_bcId7902;
  MFloat m_7902Position;
  MInt m_7902faceNormalDir;
  MInt m_7902wallDir;
  MInt m_7902periodicDir;

  MFloat m_averagePos;
  MInt m_averageDir;
  MInt m_averageFaceDir;

  const MFloat eps = 1e-16;
  const MFloat epss = 1e-8;

 protected:
  struct LESVarAverageData {
    static constexpr const MInt noAvgVars = 12;
    static constexpr const MInt UM = 0;
    static constexpr const MInt VM = 1;
    static constexpr const MInt WM = 2;
    static constexpr const MInt RHOM = 3;
    static constexpr const MInt PM = 4;
    static constexpr const MInt NUT = 5;
    static constexpr const MInt UU = 6;
    static constexpr const MInt UV = 7;
    static constexpr const MInt UW = 8;
    static constexpr const MInt VV = 9;
    static constexpr const MInt VW = 10;
    static constexpr const MInt WW = 11;
  };
};


#endif // ifndef ZONALRTV_H_
