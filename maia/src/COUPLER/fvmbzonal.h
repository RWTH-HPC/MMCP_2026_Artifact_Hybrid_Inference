// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLERFVMBZONAL_H_
#define COUPLERFVMBZONAL_H_

#include <map>
#include <vector>
#include "FV/fvmbcartesiansolverxd.h"
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
class CouplerFvMbZonal final : public CouplingFvMb<nDim, SysEqn> {
 public:
  using FvMbSolver = FvMbCartesianSolverXD<nDim, SysEqn>;
  using FVBndryCnd = FvBndryCndXD<nDim, SysEqn>;


  CouplerFvMbZonal(const MInt couplingId, FvMbSolver* upStream, FvMbSolver* downStream);
  ~CouplerFvMbZonal(){};

  // empty main functions
  void init() override{};
  void finalizeSubCoupleInit(MInt){};
  void postCouple(MInt) override;
  void cleanUp(){};
  void checkProperties() override{};

  // initialisation functions
  void readProperties() override;
  void finalizeCouplerInit();
  void finalizeBalance(const MInt) override;
  void finalizeAdaptation(const MInt) override;

  // update functions
  void subCouple(MInt, MInt, std::vector<MBool>&);
  void preCouple(MInt) override;

 private:
  FvMbSolver& upStream() const { return *m_solverUp; }
  FvMbSolver& downStream() const { return *m_solverDown; }
  MInt couplerId() const { return m_couplingId; }

  FvMbSolver* m_solverUp;
  FvMbSolver* m_solverDown;

  MInt m_upStreamId{};
  MInt m_downStreamId{};
  MInt m_couplingId{};

  MString m_zonalMethod{};
  MInt m_zonalDir{};
  MFloat m_zonalCoordinate{};

  MBool m_zonalDualTimeStepping = false;

  std::vector<std::pair<MInt, MInt>> m_upDown;
  std::vector<std::pair<MInt, MInt>> m_downUp;

  MInt m_upStreamOffset = 0;
  MInt m_downStreamOffset = 0;
  MInt m_RKStep = 0;
  MInt m_noRKSteps = -1;

  void createZonalMapping();
  void updateZonalMapping();
  void exchangeZonalValues(const MInt);
  void unifyTimeStep();

  MInt up2downId(MInt, const MBool);
  MInt down2upId(MInt, const MBool);
};

#endif // #ifndef COUPLERFVMULTILEVEL_H_
