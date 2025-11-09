// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVZONALSTG_H_
#define FVZONALSTG_H_

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
class FvZonalSTG : public FvZonal<nDim, SysEqn> {
 private:
  friend class FvZonal<nDim, SysEqn>;
  // friend class FvZonal<nDim, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>;

 public:
  using RANS = FvCartesianSolverXD<nDim, SysEqn>;
  using LES = FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>;

  // Constructor
  FvZonalSTG(const MInt couplingId, RANS* R, LES* L);

  using Base = FvZonal<nDim, SysEqn>;
  using Base::a_noFvCellsLES;
  using Base::a_noFvCellsRANS;
  using Base::a_noFvGridCellsLES;
  using Base::a_noFvGridCellsRANS;
  using Base::getAveragingFactor;
  using Base::initLESValues;
  using Base::initRANSValues;
  using Base::LESSolver;
  using Base::m_averageTimeSteps;
  using Base::m_azimuthalAngle;
  using Base::m_cylindricCommunication;
  using Base::m_LESNoVarAverage;
  using Base::m_LESSolverId;
  using Base::m_noReconstructNutVars;
  using Base::m_RANSSolverId;
  using Base::m_restartLESAverage;
  using Base::m_STGSponge;
  using Base::m_zonalAveragingTimeStep;
  using Base::m_zonalTransferInterval;
  using Base::noExchangeVariables;
  using Base::noLESVariables;
  using Base::noRANSVariables;
  using Base::RANSSolver;

  void init();
  void finalizeCouplerInit();
  void preCouple(MInt);

  void finalizeAdaptation(const MInt solverId) override;
  void finalizeBalance(const MInt) override{};
  void balancePre() override{};
  void balancePost() override;

 private:
  // suppressed functions
  void postCouple(MInt){};
  void finalizeSubCoupleInit(MInt){};
  void subCouple(MInt, MInt, std::vector<MBool>&){};
  void cleanUp(){};
  void initData(){};
  void checkProperties(){};
  void readProperties(){};

  // void finalizeLESAverage();
  void determineZonalPositions();
  void calcPeriodicSpongeAverage();
  void calcLESSectorAverage();
  void transferSolverData();
  void transferSpongeData();
  void resetRANSSolver();
  void finalizeRANSSolverAfterReset();
  void saveSpongeData();
  void loadSpongeData();

  void initSpongeExchange();
  void resetSTGSpongeAfterAdaptation(){};
  void initCylinderExchange();
  void cylinderExchange();
  void calcRANSSectorValues();

  const MFloat eps = 1e-16;
  const MFloat epss = 1e-8;

  // MBool m_wasAdapted = false;

  MFloat m_nuTildeInfinity;
  MFloat m_7901Position;
  MInt m_7901faceNormalDir;
  MInt m_7901wallDir;
  MInt m_7901periodicDir;
  MFloat m_7909Position;
  MInt m_7909faceNormalDir;
  MInt m_7909wallDir;
  MInt m_7909periodicDir;
  MInt m_bcId7909;
  MFloat m_averagePos;
  MInt m_averageDir;

  std::vector<MInt> m_periodicSpongeCylinderExchangeIndex;
  std::vector<MInt>* m_periodicSpongeInterpolationIndex = nullptr;

  MFloat* m_uvErr = nullptr;
  MFloat* m_uvRans = nullptr;
  MFloat* m_uvInt = nullptr;
  const MFloat alpha = 10.0;
  const MFloat beta = 2.0;

  // scale uv of RANS solution
  MFloat m_uvRANSFactor;

  MFloat* m_RANSSectors;
  MFloat* m_RANSSectorLimits;

  MBool m_cylCommActive = false;
  MInt m_commSizeCylExchange = 0;
  MInt m_cylRoot = -1;
  MPI_Comm m_commCyl;
  MPI_Comm m_commStg;
  MInt m_noRANSCylinderExchangeVariables;
  MInt m_noCylindricalGlobalExchangeLocations;
  MInt m_noCylindricalGlobalRANSExchangeValues;
  MInt m_noCylindricalGlobalLESExchangeValues;
  MInt m_noCylindricalGlobalExchangeIds;
  MInt m_noRANSExchangeCells;
  MInt m_noGlobalRANSExchangeCells;
  MFloat* m_globalCylinderExchangeLocations = nullptr;
  MFloat* m_globalCylinderRANSExchangeValues = nullptr;
  MFloat* m_globalCylinderLESExchangeValues = nullptr;
  MInt* m_globalCylinderExchangeIds = nullptr;
  MInt m_cylinderExchangeIdsOffset;
  MFloat* m_cylinderInterpolationAngle = nullptr;
  std::vector<MInt>* m_globalCylinderInterpolationIndex = nullptr;
  std::vector<std::pair<MInt, MFloat>>* m_globalCylinderInterpolationCell = nullptr;
  std::vector<MInt> m_globalCylinderInterpolationNumber;


 protected:
  struct LESVarAverageData {
    static constexpr const MInt noAvgVars = 11;
    static constexpr const MInt UM = 0;
    static constexpr const MInt VM = 1;
    static constexpr const MInt WM = 2;
    static constexpr const MInt RHOM = 3;
    static constexpr const MInt PM = 4;
    static constexpr const MInt UU = 5;
    static constexpr const MInt UV = 6;
    static constexpr const MInt UW = 7;
    static constexpr const MInt VV = 8;
    static constexpr const MInt VW = 9;
    static constexpr const MInt WW = 10;
  };
};


#endif // ifndef ZONALSTG_H_
