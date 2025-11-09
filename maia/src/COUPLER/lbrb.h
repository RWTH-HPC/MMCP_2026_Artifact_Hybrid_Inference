// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBRB_H_
#define LBRB_H_

#include "UTIL/pointbox.h"
#include "UTIL/timer.h"
#include "cartesiansolver.h"
#include "enums.h"

#include "coupling.h"
#include "lslb.h"
#include "surfacecoupling.h"

#include <algorithm>

template <MInt nDim>
class CouplingRigidBodies;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;

template <MInt nDim, MInt nDist, class SysEqn>
class LbRb : public CouplingLB<nDim, nDist, SysEqn>, public CouplingRigidBodies<nDim> {
 public:
  // Simplifications
  using CouplingLb = CouplingLB<nDim, nDist, SysEqn>;
  using CouplingRb = CouplingRigidBodies<nDim>;

  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqn>;
  using RBodies = typename CouplingRb::RBodies;

  // Constructor
  LbRb<nDim, nDist, SysEqn>(MInt couplingId, LbSolver* lb, RBodies* rb);
  ~LbRb<nDim, nDist, SysEqn>();

  MInt m_lbSolverId{};
  MString m_timerType;

  maia::coupling::Mapping bndryToBodyMapping;
  maia::coupling::Mapping bodyToBndryMapping;
  MFloat** forces = nullptr;

  static constexpr MInt m_noCorners = (nDim == 2) ? 4 : 8;
  static constexpr MInt nRot = (nDim == 3) ? 3 : 1;

  MBool m_static_updateLevelSetFlowSolver_firstRun = true;

  using CouplingLb::a_associatedBodyIdsMb;
  using CouplingLb::a_cellLengthAtLevel;
  using CouplingLb::a_childId;
  using CouplingLb::a_initTemperatureKelvin;
  using CouplingLb::a_isActive;
  using CouplingLb::a_isThermal;
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
  using CouplingLb::a_pvrho;
  using CouplingLb::a_pvt;
  using CouplingLb::a_pvu;
  using CouplingLb::a_pvv;
  using CouplingLb::a_pvw;
  using CouplingLb::a_Re;
  using CouplingLb::a_variable;
  using CouplingLb::a_wasActive;
  using CouplingLb::lbBndCnd;
  using CouplingLb::lbSolver;
  using CouplingLb::minCell;
  using CouplingLb::noMinCells;

  using CouplingRb::a_noCollectorBodies;
  using CouplingRb::bodies;

 private:
  // Conversion-Factors
  struct ConversionFactors {
    MFloat velocity{};
    MFloat time{};
    MFloat length{};
    MFloat force{};
    MFloat torque{};
    MFloat pressure{};
  };

  ConversionFactors conversionRbLb;
  ConversionFactors conversionLbRb;

  struct Timers {
    enum {
      TimerGroup,
      Class,

      Constructor,
      CouplePostLb,
      CouplePostRb,

      ConstructGField,
      Preparation,

      FindBoundaryCells,
      InitSolidDomain,
      SetBoundaryVelocity,
      CreateComm,
      ApplyBC,
      ComputeBodyForces,

      PreCouple,
      FindBoundaryMapping,

      _count
    };
  };

  std::array<MInt, Timers::_count> m_timers;

  //--------------------------------- functions
  //------------------------------------------------------

 private:
  void initData();
  void initTimers();
  void averageTimer();
  void checkProperties();
  void readProperties();

  void updateGeometry();

 public:
  void init() override;
  void finalizeSubCoupleInit(MInt /*solverId*/) override{};
  void finalizeCouplerInit() override{};

  void preCouple(MInt /*recipeStep = 0*/) override{};
  void subCouple(MInt /*recipeStep = 0*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override;
  void postCouple(MInt /*recipeStep = 0*/) override{};

  void postAdaptation() override;
  void finalizeAdaptation(const MInt solverId) override;

  void cleanUp() override{};
  void reinitAfterBalance(){};
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& /*domainInfo*/) override{};

  // LB stuff
  void initializeSolidDomain();

  // Body stuff
  void getBodyVelocity(const MInt body, MFloat* const velocity);
  void getBodyAngularVelocity(const MInt body, MFloat* const angularVelocity);
  void createBodyTree();

  void bc3060(MInt);

  void transferLevelSetFieldValues(MBool);

  // LS stuff
  void computeGCellTimeStep();
  MInt returnNoActiveCorners(MInt);
  void returnLevelSetSignForFluidCellCorners(MInt, MIntScratchSpace* levelSetCornerSigns, MInt set = 0);
  MInt noLevelSetFieldData();
  MInt returnLevelSetSignForFluidFaceCentroid(MInt, MInt, MInt set = 0);

  void constructGField();
  template <MInt bodyType>
  void constructGField_();
  template <MInt bodyType>
  void descendLevelSetValue(const MInt cellId, const MInt* bodyId, const MInt bodyCnt);
};

#endif // ifndef LBRB_H_
