// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVINTERPOLATION_H_
#define FVINTERPOLATION_H_

#include <algorithm>
#include "FV/fvcartesiancellcollector.h"
#include "FV/fvcartesiancellproperties.h"
#include "GRID/cartesiangrid.h"
#include "coupling.h"
#include "couplingutils.h"
#include "enums.h"

using namespace std;

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

class Coupling;

/** Coupling class for nonZonalRestart
    (can be expended to a general FV-FV coupling class)
*   \author Jannik Borgelt
*/
template <MInt nDim_, class SysEqnOld, class SysEqnNew>
class FvCartesianInterpolation : virtual public Coupling {
 public:
  static constexpr MInt nDim = nDim_;

  // using CartesianSolver = maia::CartesianSolver<nDim, Solver>;
  using OldFvSolver = FvCartesianSolverXD<nDim, SysEqnOld>;
  using NewFvSolver = FvCartesianSolverXD<nDim, SysEqnNew>;

  FvCartesianInterpolation(const MInt couplingId, OldFvSolver* oldS, NewFvSolver* newS);
  ~FvCartesianInterpolation(){};

  void init() override;
  void preCouple(MInt) override;
  void postCouple(MInt) override;

 private:
  OldFvSolver* m_oldSolver;
  NewFvSolver* m_newSolver;

  OldFvSolver& oldSolver() const { return *m_oldSolver; }
  NewFvSolver& newSolver() const { return *m_newSolver; }

  MInt noExchangeVariables() { return nDim + 2 /*+ (MInt)m_reconstructNut*/; };
  // suppressed functions
  void finalizeSubCoupleInit(MInt){};
  void finalizeCouplerInit(){};
  void subCouple(MInt, MInt, std::vector<MBool>&){};
  void cleanUp(){};
  void initData(){};
  void checkProperties(){};
  void readProperties(){};

  void transferSolverData();

  MInt a_noFvCellsOld() const { return oldSolver().a_noCells(); }
  MInt a_noFvGridCellsOld() const { return oldSolver().c_noCells(); }
  MInt a_noFvCellsNew() const { return newSolver().a_noCells(); }
  MInt a_noFvGridCellsNew() const { return newSolver().c_noCells(); }

  const MFloat eps = 1e-16;
  const MFloat epss = 1e-8;

  MInt m_oldSolverId;
  MInt m_newSolverId;

  MBool m_nonZonalRestart;
};

#endif // ifndef FVINTERPOLATION_H_
