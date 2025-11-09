// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBLPT_H_
#define LBLPT_H_

#include "LPT/lpt.h"
#include "coupling.h"
#include "couplinglpt.h"
#include "solver.h"

// Forward declarations
template <MInt nDim>
class CouplingParticle;

class LbSolverD2Q9;
template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;
template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;

namespace maia {
namespace lblpt {} // namespace lblpt
} // namespace maia

template <MInt nDim, MInt nDist, class SysEqn>
class LbLpt : public CouplingLpt<nDim, CouplingLB<nDim, nDist, SysEqn>> {
 public:
  using BaseLb = CouplingLB<nDim, nDist, SysEqn>;
  using LbSolver = typename BaseLb::solverType;

  using BaseLpt = CouplingParticle<nDim>;
  using BaseLptX = CouplingLpt<nDim, CouplingLB<nDim, nDist, SysEqn>>;
  using LptSolver = typename BaseLpt::solverType;

  // Constructor
  LbLpt(const MInt couplingId, LPT<nDim>* particle, LbSolver* lb);

 protected:
  LPT<nDim>& lpt() const override { return *m_particle; }

  using BaseLb::a_bndCellId;
  using BaseLb::a_cellLengthAtLevel;
  using BaseLb::a_Ma;
  using BaseLb::a_noBndCells;
  using BaseLb::a_noCells;
  using BaseLb::lbBndCnd;
  using BaseLb::lbSolver;

  using typename BaseLptX::ConversionFactor;

  ConversionFactor& conversionLbLpt = BaseLptX::conversionFlowLpt;
  ConversionFactor& conversionLptLb = BaseLptX::conversionLptFlow;

 public:
  void init() override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt /*unused*/) override{};

  void preCouple(MInt) override;
  void postCouple(MInt) override;

  void finalizeBalance(const MInt /*unused*/) override;
  void balancePre() override{};

  void cleanUp() override{};

  // pass-through calls in the LPT-solver
  // TODO: remove these coupler calls in the next step!
  MFloat a_lbVariable(const MInt lptCellId, const MInt varId, const MInt mode) {
    MInt lbCellId = lpt2lbIdParent(lptCellId);
    if(lbCellId < 0 || lbCellId > lbSolver().a_noCells()) return -99;
    if(mode == 0) {
      return lbSolver().a_variable(lbCellId, varId);
    } else {
      return lbSolver().a_oldVariable(lbCellId, varId);
    }
  }

  MFloat a_lbCellCoordinate(const MInt lptCellId, const MInt dir) {
    MInt lbCellId = lpt2lbIdParent(lptCellId);
    if(lbCellId < 0 || lbCellId > lbSolver().a_noCells()) return -99;
    return lbSolver().a_coordinate(lbCellId, dir);
  }

 private:
  void checkProperties() override;
  void readProperties();
  void initData();
  void initConversion();
  void calculateGridBoundaries();
  void initParticleVelocity();

  void updateLbSolver();

  // void updateFractions();
  // void updateForcing();
  void updateLPTBndry();
  void transferFlowField();
  void transferCellStatus();
  void transferTimeStep();

  void computeParticleInterphaseExchangeRate();

  MInt lpt2lbId(const MInt lptId) { return convertId(lpt(), lbSolver(), lptId); };
  MInt lpt2lbIdParent(const MInt lptId) { return convertIdParent(lpt(), lbSolver(), lptId); };
  MInt lb2lptId(const MInt lbId) { return convertId(lbSolver(), lpt(), lbId); };
  MInt lb2lptIdParent(const MInt lbId) { return convertIdParent(lbSolver(), lpt(), lbId); };

  LPT<nDim>* m_particle = nullptr;
  LbSolver* m_lbSolver = nullptr;

  MInt m_lptSolverId{};
  MInt m_lbSolverId{};
  MInt m_lptSolverOrder{};
  MInt m_noSolverSteps{};
  std::array<MFloat, 2 * nDim> m_gridBoundaries{};

  // Properties
  MBool m_lptlbInterpolation = false;
  MBool m_forceLbTimeStep = true;
  MBool m_CalcSurface = false;
};

#endif // ifndef LBLPT_H_
