// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLERLBFV_H_
#define COUPLERLBFV_H_

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

#include "coupling.h"
#include "couplingutils.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;
template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy;

template <class SysEqn>
class FvCartesianSolver2D;
template <class SysEqn>
class FvCartesianSolver3D;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;
template <MInt nDim, class SysEqn>
class CouplingFv;


template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
class CouplerLbFv : public CouplingFv<nDim, SysEqnFv>, public CouplingLB<nDim, nDist, SysEqnLb> {
 private:
  friend class LbSolverDxQy<nDim, nDist, SysEqnLb>;
  friend class LbBndCndDxQy<nDim, nDist, SysEqnLb>;
  friend class FvCartesianSolverXD<nDim, SysEqnFv>;

 public:
  // Simplifications
  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqnLb>;
  using FvCartesianSolver = FvCartesianSolverXD<nDim, SysEqnFv>;

  using BaseLb = CouplingLB<nDim, nDist, SysEqnLb>;
  using BaseFv = CouplingFv<nDim, SysEqnFv>;

  using BaseLb::lbBndCnd;

  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  // Constructor
  CouplerLbFv(const MInt couplingId, LbSolver* lb, FvCartesianSolver* fv);

 private:
  MInt m_lbSolverId{};
  MInt m_fvSolverId{};

  MFloat m_cfl{};
  MFloat m_maxVelocity{};
  MInt m_timeStepMethod{};
  MString m_solverMethod;

 protected:
  using BaseLb::a_cellLengthAtLevel;
  using BaseLb::a_childId;
  using BaseLb::a_noCells;
  using BaseLb::a_noLbCells;
  using BaseLb::a_parentId;
  using BaseLb::lbSolver;
  using BaseLb::minCell;
  using BaseLb::noMinCells;

  using BaseFv::a_noFvGridCells;
  using BaseFv::fvSolver;

  //--------------------------------- functions ------------------------------------------------------

 private:
  void initData();

  void cleanUp() override{};

  // Conversion-Factors
  struct ConversionFactors {
    MFloat length{};
    MFloat velocity{};
    MFloat pressure{};
    MFloat viscosity{};
  };


 public:
  ConversionFactors conversionLbFv;
  ConversionFactors conversionFvLb;

  virtual void init(){};
  virtual void initConversionFactors();
  void postCouple(MInt){};
  virtual void finalizeCouplerInit(){};
  virtual void finalizeSubCoupleInit(MInt){};
  virtual void subCouple(MInt, MInt, std::vector<MBool>&){};
  virtual void checkProperties();
  virtual void readProperties();

  MInt fv2lbId(const MInt fvId) const { return convertId(fvSolver(), lbSolver(), fvId); };
  MInt lb2fvId(const MInt lbId) const { return convertId(lbSolver(), fvSolver(), lbId); };

  void testCoupling();

  MInt a_fvSolverId() const { return m_fvSolverId; }
  MInt a_lbSolverId() const { return m_lbSolverId; }
};

#endif // ifndef COUPLERLBFV_H_
