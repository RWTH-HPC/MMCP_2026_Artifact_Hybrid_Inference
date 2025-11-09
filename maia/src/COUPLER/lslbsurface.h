// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSLBSURFACE_H_
#define LSLBSURFACE_H_

#include "UTIL/pointbox.h"
#include "UTIL/timer.h"
#include "cartesiansolver.h"
#include "enums.h"

#include "coupling.h"
#include "surfacecoupling.h"

#include <algorithm>

template <MInt nDim>
class CouplingLS;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;

template <MInt nDim, MInt nDist, class SysEqn>
class LsLbSurface : public CouplingLS<nDim>, public CouplingLB<nDim, nDist, SysEqn> {
 public:
  // Simplifications
  using LsSolver = LsCartesianSolver<nDim>;
  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqn>;
  using CouplingLs = CouplingLS<nDim>;
  using CouplingLb = CouplingLB<nDim, nDist, SysEqn>;

  // Constructor
  LsLbSurface<nDim, nDist, SysEqn>(MInt couplingId, LsSolver* ls, LbSolver* lb);
  ~LsLbSurface<nDim, nDist, SysEqn>();

  MInt m_lbSolverId{};
  MInt m_lsSolverId{};

  static constexpr MInt m_noCorners = (nDim == 2) ? 4 : 8;

  MBool m_static_updateLevelSetFlowSolver_firstRun = true;

  using CouplingLs::a_bodyIdG;
  using CouplingLs::a_bodyToSet;
  using CouplingLs::a_coordinateG;
  using CouplingLs::a_curvatureG;
  using CouplingLs::a_extensionVelocityG;
  using CouplingLs::a_G0CellId;
  using CouplingLs::a_levelSetFunctionG;
  using CouplingLs::a_noEmbeddedBodies;
  using CouplingLs::a_noG0Cells;
  using CouplingLs::a_noLsCells;
  using CouplingLs::a_normalVectorG;
  using CouplingLs::a_noSets;
  using CouplingLs::a_outsideGValue;
  using CouplingLs::lsSolver;

  using CouplingLb::a_associatedBodyIdsMb;
  using CouplingLb::a_boundaryCellMb;
  using CouplingLb::a_cellLengthAtLevel;
  using CouplingLb::a_childId;
  using CouplingLb::a_isActive;
  using CouplingLb::a_levelSetFunctionMb;
  using CouplingLb::a_Ma;
  using CouplingLb::a_mbCell;
  using CouplingLb::a_noCells;
  using CouplingLb::a_noDistributions;
  using CouplingLb::a_noEmbeddedBodiesLB;
  using CouplingLb::a_noLbCells;
  using CouplingLb::a_noLevelSetsMb;
  using CouplingLb::a_noVariables;
  using CouplingLb::a_oldVariable;
  using CouplingLb::a_parentId;
  using CouplingLb::a_Re;
  using CouplingLb::a_variable;
  using CouplingLb::a_wasActive;
  using CouplingLb::lbBndCnd;
  using CouplingLb::lbSolver;
  using CouplingLb::minCell;
  using CouplingLb::noMinCells;

 private:
  // Properties
  MBool m_calcWallForces{};

  // Conversion-Factors
  struct ConversionFactors {
    MFloat velocity{};
    MFloat time{};
    MFloat length{};
    MFloat force{};
    MFloat torque{};
    MFloat pressure{};
  };

  ConversionFactors conversionLsLb;
  ConversionFactors conversionLbLs;

  MFloat m_surfaceTension = 0.0;
  MFloat m_gravity{};
  MFloat m_Ga = 0.0;
  MFloat m_Eo = 1.0;
  MFloat m_initCurvature = 0.0;
  MFloat m_initHeight = 0.0;

  struct Timers {
    enum {
      TimerGroup,
      Class,

      Constructor,

      PreCouple,
      TransferToLevelSet,
      SetExtensionVelocity,
      ExtendVelocity,
      FindBoundaryCells,
      InitOutsideDomain,
      SetBoundaryCondition,

      PostCouple,
      ComputeBoundaryValues,
      CreateComm,
      ApplyBoundaryCondition,

      _count
    };
  };

  std::array<MInt, Timers::_count> m_timers;

  // Data
  maia::coupling::Mapping bndryToVolumeMap{};
  maia::coupling::Mapping volumeToBndryMap{};

  //--------------------------------- functions
 private:
  void initData();
  void initTimers();
  void checkProperties();
  void readProperties();

  void updateGeometry();

 public:
  void init() override;
  void finalizeSubCoupleInit(MInt /*solverId*/) override{};
  void finalizeCouplerInit() override{};

  void preCouple(MInt /*recipeStep = 0*/) override;
  void subCouple(MInt /*recipeStep = 0*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override{};
  void postCouple(MInt recipeStep = 0) override;

  void cleanUp() override{};
  void reinitAfterBalance(){};

  MInt ls2lbId(MInt);
  MInt lb2lsId(MInt);
  MInt ls2lbIdParent(MInt);
  MInt lb2lsIdParent(MInt);

  void updateBoundaryCellsFromGField();

  // LB stuff
  void refillEmergedCells();

  void interpolateCurvature(MFloatScratchSpace& curvature);
  void interpolateNormal();
  void evaluateContour();

  // Body stuff

  void bc3060(MInt);

  void transferLevelSetFieldValues(MBool);

  // LS stuff
  void setExtensionVelocity();
  void setExtensionVelocityB();
  void computeGCellTimeStep();
};

#endif // ifndef LSLBSURFACE_H_
