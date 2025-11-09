// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSFVCOMBUSTION_H_
#define LSFVCOMBUSTION_H_

#include "FV/fvcartesiancellcollector.h"
#include "LS/lscartesiansolver.h"
#include "FV/fvcartesiancellproperties.h"
#include "coupling.h"
#include "couplingutils.h"


template <MInt nDim>
class CouplingLS;

template <MInt nDim, class SysEqn>
class CouplingFv;

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim_, class SysEqn>
class LsFvCombustion : public Coupling {
 public:
  static constexpr MInt nDim = nDim_;

 private:
  // template <MInt nDim, SysEqn> friend class FvCartesianSolverXD;
  // template <SysEqn> friend class FvMbSolver2D;
  // template <SysEqn> friend class FvMbSolver3D;
  friend class LsCartesianSolver<nDim>;
  friend class FvCartesianSolverXD<nDim, SysEqn>;

 public:
  // Type for cell properties
  using Cell = typename LsCartesianSolver<nDim_>::Cell;
  using SolverCell = FvCell;

  using FvCartesianSolver = FvCartesianSolverXD<nDim_, SysEqn>;
  using LsSolver = LsCartesianSolver<nDim_>;

  LsSolver* m_lsSolver;
  LsSolver& lsSolver() const { return *m_lsSolver; }

  FvCartesianSolver* m_fvSolver;
  FvCartesianSolver& fvSolver() const { return *m_fvSolver; }

  LsFvCombustion<nDim_, SysEqn>(const MInt couplingId, LsSolver* ls, FvCartesianSolver* fv);

  MInt a_noSets() const { return lsSolver().m_noSets; }
  MFloat a_levelSetFunctionG(MInt gcellId, MInt set) const { return lsSolver().a_levelSetFunctionG(gcellId, set); }
  MFloat a_outsideGValue() const { return lsSolver().m_outsideGValue; }

  void init();
  void finalizeSubCoupleInit(MInt){};
  void finalizeCouplerInit() override;
  void preCouple(MInt) override;
  void subCouple(MInt, MInt, std::vector<MBool>&){};
  void postCouple(MInt);
  void cleanUp(){};
  void checkProperties(){};
  void readProperties(){};

  void computeSourceTerms();
  MFloat collectGFromCouplingClass(MInt);
  MFloat collectCurvFromCouplingClass(MInt);
  MInt noLevelSetFieldData();
  void saveOutputLS();
  void setRhoFlameTubeInLs();
  void setRhoInfinityInLs();
  MFloat collectFvDataForCollectGEquationModelDataOpt(MInt, MInt);
  MInt m_maxNoSets;
  void computeGCellTimeStep();
  void setLsTimeStep(MFloat);

  // Id conversion
  MInt ls2fvId(const MInt lsId) { return convertId(lsSolver(), fvSolver(), lsId); };
  MInt fv2lsId(const MInt fvId) { return convertId(fvSolver(), lsSolver(), fvId); };

  void collectGEquationModelDataOptInterpolate(MFloat* fluidDensity, MInt set);
  void fastInterfaceExtensionVelocity();
  void exchangeCouplingData();
  void constructExtensionVelocity();
  void collectGEquationModelDataOpt(MFloat* fluidDensity, MInt set);
};

#endif // ifndef LSFVCOMBUSTION_H_
