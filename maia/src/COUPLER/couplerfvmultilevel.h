// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLERFVMULTILEVEL_H_
#define COUPLERFVMULTILEVEL_H_

#include <map>
#include <vector>
#include "INCLUDE/maiatypes.h"
#include "coupling.h"

// Forward declarations
class Solver;
template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;
namespace maia {
namespace grid {
template <MInt nDim>
class Proxy;
} // namespace grid
} // namespace maia
template <MInt nDim>
class CartesianGrid;


template <MInt nDim, class SysEqn>
class CouplerFvMultilevel : public CouplingFv<nDim, SysEqn> {
 public:
  using Base = Coupling;
  using BaseFv = CouplingFv<nDim, SysEqn>;
  using SolverType = typename BaseFv::solverType;
  using Grid = CartesianGrid<nDim>;
  using GridProxy = typename SolverType::GridProxy;

  using BaseFv::fvSolver;
  using BaseFv::noSolvers;

  CouplerFvMultilevel(const MInt couplingId, std::vector<FvCartesianSolverXD<nDim, SysEqn>*> solvers);
  ~CouplerFvMultilevel();

  // Main coupling functions
  void init() override;
  void finalizeSubCoupleInit(MInt){};
  void finalizeCouplerInit();
  void preCouple(MInt) override;
  void subCouple(MInt, MInt, std::vector<MBool>&) override{};
  void postCouple(MInt) override;
  void cleanUp(){};
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& NotUsed(domainInfo)) override{};
  void finalizeAdaptation(MInt solverId);

  void readProperties() override{};
  void checkProperties() override{};

  void restriction(const MInt level);
  void resetTau(const MInt level); // TODO labels:COUPLER,toremove not used
  void prolongation(const MInt level);
  void sanityCheck();
  MBool startTimer(const MString& name);
  MBool stopTimer(const MString& name);

 private:

  // Multilevel methods
  void restrictData(const MInt level);
  void computeCoarseLevelCorrection(const MInt level);
  void storeRestrictedVariables(const MInt level);
  void prolongData(const MInt level);

  void correctCoarseBndryCells(const MInt solverId);

  void resetSecondaryFlowVariables();
  void resetMultiLevelVariables();
  void initLeafCells();
  void initMapping();
  void initSplitMapping();
  void initRestrictedCells();

  void sanityCheck(const MInt mode = 0);

  // Timer management
  MInt timer(const MString& name);
  MInt& createTimer(const MString& name);
  std::map<MString, MInt> m_timers{};
  MInt m_timerGroup = -1;
  void initTimers();

  // Multilevel data:
  std::vector<std::vector<MInt>> m_leafCells{};
  std::vector<std::vector<MInt>> m_coarse2fine{};
  std::vector<std::vector<MInt>> m_fine2coarse{};
  std::vector<std::vector<MInt>> m_restrictedCells{};
  std::vector<std::vector<MInt>> m_splitCellMapping{};

  // Multilevel properties:
  MInt m_prolongationMethod = false;
  MBool m_correctCoarseBndry = false;
};


/// \brief FV multilevel interpolation coupler to transfer solution data of a coarse grid onto a finer grid
template <MInt nDim, class SysEqn>
class CouplerFvMultilevelInterpolation final : public CouplerFvMultilevel<nDim, SysEqn> {
 public:
  using BaseCoupler = CouplerFvMultilevel<nDim, SysEqn>;
  using Base = typename BaseCoupler::Base;
  using BaseFv = typename BaseCoupler::BaseFv;
  using BaseCoupler::fvSolver;
  using BaseCoupler::noSolvers;
  using BaseCoupler::prolongation;
  using BaseCoupler::startTimer;
  using BaseCoupler::stopTimer;

  CouplerFvMultilevelInterpolation(const MInt couplingId,
                                   // std::vector<std::unique_ptr<Solver>>* const solvers
                                   std::vector<FvCartesianSolverXD<nDim, SysEqn>*>
                                       solvers);
  ~CouplerFvMultilevelInterpolation() { TRACE(); };

  // Override the coupling routine of the FvMultilevel coupler
  void postCouple(MInt) override;

 private:
  void reset(const MInt level);
};

#endif // #ifndef COUPLERFVMULTILEVEL_H_
