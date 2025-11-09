// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSFVMB_H_
#define LSFVMB_H_

#include "FV/fvmbcartesiansolverxd.h"
#include "LS/lscartesiansolver.h"
#include "solver.h"

#include "coupling.h"
#include "couplingutils.h"

template <MInt nDim>
class CouplingLS;

template <MInt nDim, class SysEqn>
class CouplingFvMb;

template <MInt nDim_, class SysEqn>
class LsFvMb : public CouplingLS<nDim_>, public CouplingFvMb<nDim_, SysEqn> {
 public:
  static constexpr MInt nDim = nDim_;

 private:
  friend class FvMbCartesianSolverXD<nDim, SysEqn>;
  friend class LsCartesianSolver<nDim>;

 public:
  // Simplifications
  using LsSolver = LsCartesianSolver<nDim>;
  using FvMbSolver = FvMbCartesianSolverXD<nDim, SysEqn>;
  using FVBndryCnd = FvBndryCndXD<nDim, SysEqn>;

  // Type for cell properties
  using LsCell = typename LsSolver::Cell;
  // FvCell is already set in one of the headers and will be used as such here!
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  // Constructor
  LsFvMb<nDim_, SysEqn>(const MInt couplingId, LsSolver* ls, FvMbSolver* fvMb);

 private:
  using Base = CouplingFvMb<nDim_, SysEqn>;
  using Base::startLoadTimer;
  using Base::stopLoadTimer;

  using Base::returnIdleRecord;
  using Base::returnLoadRecord;

  MInt m_fvSolverId{};
  MInt m_lsSolverId{};

  MBool m_outsideDefault{};

  MBool m_allowLsInterpolation = false;
  MInt m_noRfJumps = 0;

  MInt* m_hadGapCells{};
  MInt* m_hasGapCells{};

  static constexpr const MInt m_noCorners = (nDim == 2) ? 4 : 8;
  static constexpr const MInt m_maxNoGapRegions = 5;

  // updateLevelSetFlowSolver
  MBool m_static_updateLevelSetFlowSolver_firstRun = true;

  // setGapState
  MBool m_static_setGapState_first = true;
  MBool m_static_setGapState_earlyOpened[m_maxNoGapRegions];

  // Call to solvers and base class functions
  using CouplingFvMb<nDim_, SysEqn>::fvMbSolver;
  using CouplingFvMb<nDim_, SysEqn>::a_noFvCells;
  using CouplingFvMb<nDim_, SysEqn>::a_noFvGridCells;
  using CouplingFvMb<nDim_, SysEqn>::a_noLevelSetsMb;

  using CouplingLS<nDim_>::lsSolver;
  using CouplingLS<nDim_>::a_noLsCells;
  using CouplingLS<nDim_>::a_outsideGValue;
  using CouplingLS<nDim_>::a_levelSetFunctionG;
  using CouplingLS<nDim_>::a_bodyIdG;
  using CouplingLS<nDim_>::a_potentialGapCellClose;
  using CouplingLS<nDim_>::a_inBandG;
  using CouplingLS<nDim_>::a_noG0Cells;
  using CouplingLS<nDim_>::a_coordinateG;
  using CouplingLS<nDim_>::a_nearGapG;

  //--------------------------------- functions ------------------------------------------------------

 private:
  void initData();
  void checkProperties();
  void readProperties();

  void interpolateLsFV(const MInt, const MInt);
  MFloat interpolateLevelSetMb(MInt* interpolationCells, MFloat* point, const MInt set);
  MFloat interpolateLevelSet(MInt cellId, MFloat* point, MInt set);
  void buildCollectedLevelSet(const MInt);

  void updateFlowSolver();

  void transferTimeStep();
  void transferAngularVelocity();
  void transferBodyRadius();
  void transferGapCellProperty(MInt mode);
  void transferBodyProperties();
  void transferLevelSetValues();
  void transferLevelSetFieldValues(MBool);

  MFloat restartTime();

  void setGapState();
  void setGapGhostCellVariables(MInt);

  void updateLevelSet();
  void updateLevelSetOutsideBandPar();

  void initFvGapCells();
  void lsGapInfo(const MInt, MInt*);

  void initBodyProperties();

  MInt ls2fvId(const MInt lsId) { return convertId(lsSolver(), fvMbSolver(), lsId); };
  MInt ls2fvIdParent(const MInt lsId) { return convertIdParent(lsSolver(), fvMbSolver(), lsId); };
  MInt fv2lsId(const MInt fvId) {
    // TODO labels:COUPLER,FV,LS,totest Is this check needed? This functions should never be called for these modes i
    // guess..
    if(fvMbSolver().m_levelSetMb && fvMbSolver().m_constructGField) {
      return -1;
    }
    return convertId(fvMbSolver(), lsSolver(), fvId);
  };
  MInt fv2lsIdParent(const MInt fvId) {
    // TODO labels:COUPLER,FV,LS,totest Is this check needed? This functions should never be called for these modes i
    // guess..
    if(fvMbSolver().m_levelSetMb && fvMbSolver().m_constructGField) {
      return -1;
    }
    return convertIdParent(fvMbSolver(), lsSolver(), fvId);
  };

  void testCoupling();
  void testLsValues();
  void testGapProperty();
  void computeBodyProperties(MInt returnMode, MFloat* bodyData, MInt body, MFloat time);

 public:
  void init() override;

  void preCouple(MInt recepiStep) override;
  void subCouple(MInt /* recipeStep = 0*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override{};
  void postCouple(MInt recipeStep = 0) override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt) override;

  void postAdaptation() override;
  void finalizeAdaptation(const MInt) override;
  void prepareAdaptation() override;

  void finalizeBalance(const MInt) override;
  void balancePre() override{};
  void balancePost() override{};
  MInt noLevelSetFieldData();

  MInt noCouplingTimers(const MBool NotUsed(allTimings)) const override { return 2; }

  void getCouplingTimings(std::vector<std::pair<MString, MFloat>>& timings, const MBool NotUsed(allTimings)) override {
    const MString namePrefix = "c" + std::to_string(this->couplerId()) + "_";

    timings.emplace_back(namePrefix + "loadCouplerLsFvMb", returnLoadRecord());
    timings.emplace_back(namePrefix + "idleCouplerLsFvMb", returnIdleRecord());
  };
};


#endif // ifndef LSFVMB_H_
