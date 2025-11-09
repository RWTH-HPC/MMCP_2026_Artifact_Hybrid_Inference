// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CouplerLbFvEEMultiphase_H_
#define CouplerLbFvEEMultiphase_H_

#include <algorithm>
#include "FV/fvcartesiansolverxd.h"
#include "LB/lbbndcnd.h"
#include "LB/lbbndcnddxqy.h"
#include "LB/lbcellcollector.h"
#include "LB/lbcellproperties.h"
#include "LB/lbgridboundarycell.h"
#include "LB/lbsolver.h"
#include "LB/lbsolverdxqy.h"
#include "UTIL/pointbox.h"
#include "enums.h"
#include "solver.h"

#include "couplerlbfv.h"
#include "coupling.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;
template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
class CouplerLbFv;

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
class CouplerLbFvEEMultiphase : public CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv> {
 private:
  void initData();
  void initAlpha();
  void readProperties();
  void cleanUp() override{};

  void revertLbVariables();
  void revertLbDistributions();
  void revertLbOldVariables();
  void revertFvVariables();
  MBool checkAlphaConverged();

  using Base = CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>;

  using Base::a_fvSolverId;
  using Base::a_lbSolverId;
  using Base::a_noLbCells;
  using Base::conversionFvLb;
  using Base::conversionLbFv;
  using Base::fv2lbId;
  using Base::fvSolver;
  using Base::lb2fvId;
  using Base::lbBndCnd;
  using Base::lbSolver;

  using CouplingFv<nDim, SysEqnFv>::a_noFvCells;

 public:
  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqnLb>;
  using FvCartesianSolver = FvCartesianSolverXD<nDim, SysEqnFv>;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  CouplerLbFvEEMultiphase(const MInt couplerId, LbSolver* lb, FvCartesianSolver* fv);

  void initDepthcorrection();
  void correctInvalidAlpha();
  std::vector<MInt> findRedistCells(const MInt cellId, const MBool searchUp, const MFloat limit);
  void transferPressureLb2Fv(const MFloat rkAlpha, const MBool update);
  template <MBool updateGradients>
  void transferULb2Fv(const MFloat rkAlpha);
  void transferUFv2Lb();
  void transferNuLb2Fv(const MFloat rkAlpha);
  void transferAlphaFv2Lb();

  void preCouple(MInt) override;
  void postcouple(MInt){};
  void subCouple(MInt, MInt, std::vector<MBool>& solverCompleted) override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt) override;
  void init() override;
  void balancePre() override;
  void balancePost() override{};

  void revertLb();
  void revertFv();

  MInt m_alphaConvergenceCheck{};
  MInt m_maxNoAlphaIterations{};
  MFloat m_epsAlpha{};
  MInt m_initAlphaMethod{};
  MFloat m_initialAlpha{};
  MFloat m_alphaInf{};
  MBool m_redistributeAlpha{};
  MBool m_disableSubstepAlphaRedist{};
  std::array<MFloat, nDim> m_gravityRefCoords{};
  MFloat* m_depthCorrectionValues = nullptr;
  std::array<MFloat, nDim> m_depthCorrectionCoefficients{};
  MBool m_updateAfterPropagation{};
  MFloat m_interpolationFactor{};
  MBool m_updateFVBC{};
  MFloat m_alphaCeil{};
  MFloat m_alphaFloor{};
};

#endif
