// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSLB_H_
#define LSLB_H_

#include <algorithm>
#include "UTIL/pointbox.h"
#include "cartesiansolver.h"
#include "enums.h"

#include "coupling.h"
#include "couplingutils.h"
#include "surfacecoupling.h"

template <MInt nDim>
class CouplingLS;
template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;

template <MInt nDim_, MInt nDist, class SysEqn>
class LsLb : public CouplingLS<nDim_>, public CouplingLB<nDim_, nDist, SysEqn> {
 public:
  static constexpr MInt nDim = nDim_;

 public:
  // Simplifications
  using LsSolver = LsCartesianSolver<nDim>;
  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqn>;

  // Constructor
  LsLb(MInt couplingId, LsSolver* ls, LbSolver* lb);

 private:
  MInt m_lbSolverId{};
  MInt m_lsSolverId{};

  MBool m_outsideDefault{};

  maia::coupling::Mapping bndryToBodyMapping;
  maia::coupling::Mapping bodyToBndryMapping;

  static constexpr const MInt m_maxNoEmbeddedBodies = 1;
  static constexpr const MInt m_noCorners = (nDim == 2) ? 4 : 8;

  // ... computeBodyProperties
  MBool m_static_computeBodyProperties_first = true;
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_amplitude{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_freqFactor{};
  std::array<MFloat, m_maxNoEmbeddedBodies * 3> m_static_computeBodyProperties_initialBodyCenter{};
  MFloat m_static_computeBodyProperties_Strouhal{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_mu{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_mu2{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_liftStartAngle1{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_liftEndAngle1{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_liftStartAngle2{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_liftEndAngle2{};
  std::array<MFloat, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_circleStartAngle{};
  std::array<MFloat, m_maxNoEmbeddedBodies * 3> m_static_computeBodyProperties_normal{};
  std::array<MInt, m_maxNoEmbeddedBodies> m_static_computeBodyProperties_bodyToFunction{};
  MFloat m_static_computeBodyProperties_omega{};
  MFloat m_static_computeBodyProperties_rotAngle{};


  MBool m_static_updateLevelSetFlowSolver_firstRun = true;

 public:
  using CouplingLS<nDim_>::lsSolver;
  using CouplingLS<nDim_>::a_levelSetFunctionG;
  using CouplingLS<nDim_>::a_noSets;
  using CouplingLS<nDim_>::a_noLsCells;
  using CouplingLS<nDim_>::a_outsideGValue;
  using CouplingLS<nDim_>::a_bodyIdG;
  using CouplingLS<nDim_>::a_bodyToSet;
  using CouplingLS<nDim_>::a_noEmbeddedBodies;

  using CouplingLB<nDim, nDist, SysEqn>::lbSolver;
  using CouplingLB<nDim, nDist, SysEqn>::lbBndCnd;
  using CouplingLB<nDim, nDist, SysEqn>::a_noLbCells;
  using CouplingLB<nDim, nDist, SysEqn>::a_noLevelSetsMb;
  using CouplingLB<nDim, nDist, SysEqn>::a_levelSetFunctionMb;
  using CouplingLB<nDim, nDist, SysEqn>::a_associatedBodyIdsMb;
  using CouplingLB<nDim, nDist, SysEqn>::a_parentId;
  using CouplingLB<nDim, nDist, SysEqn>::a_childId;
  using CouplingLB<nDim, nDist, SysEqn>::minCell;
  using CouplingLB<nDim, nDist, SysEqn>::noMinCells;
  using CouplingLB<nDim, nDist, SysEqn>::a_noCells;
  using CouplingLB<nDim, nDist, SysEqn>::a_cellLengthAtLevel;
  using CouplingLB<nDim, nDist, SysEqn>::a_noEmbeddedBodiesLB;
  using CouplingLB<nDim, nDist, SysEqn>::a_noVariables;
  using CouplingLB<nDim, nDist, SysEqn>::a_variable;
  using CouplingLB<nDim, nDist, SysEqn>::a_oldVariable;
  using CouplingLB<nDim, nDist, SysEqn>::a_isActive;
  using CouplingLB<nDim, nDist, SysEqn>::a_wasActive;
  using CouplingLB<nDim, nDist, SysEqn>::a_Ma;
  using CouplingLB<nDim, nDist, SysEqn>::a_Re;
  using CouplingLB<nDim, nDist, SysEqn>::a_noDistributions;
  using CouplingLB<nDim, nDist, SysEqn>::a_initTemperatureKelvin;
  using CouplingLB<nDim, nDist, SysEqn>::a_pvu;
  using CouplingLB<nDim, nDist, SysEqn>::a_pvv;
  using CouplingLB<nDim, nDist, SysEqn>::a_pvw;
  using CouplingLB<nDim, nDist, SysEqn>::a_pvrho;
  using CouplingLB<nDim, nDist, SysEqn>::a_isThermal;
  using CouplingLB<nDim, nDist, SysEqn>::a_pvt;
  using CouplingLB<nDim, nDist, SysEqn>::a_mbCell;

  static constexpr const MBool m_constructGField = false;
  MFloat* m_transferBoundingBox{};
  //--------------------------------- functions ------------------------------------------------------

 private:
  void initData();
  void checkProperties();
  void readProperties();

  void updateGeometry();

 public:
  void init() override;
  void preCouple(MInt step) override;
  void postCouple(MInt step) override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt couplingStep) override;
  void subCouple(MInt /*recipeStep*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override{};
  void cleanUp() override{};

  void postAdaptation() override;
  void finalizeAdaptation(const MInt solverId) override;

  // Id conversion

  MInt ls2lbId(const MInt lsId) { return convertId(lsSolver(), lbSolver(), lsId); };
  MInt ls2lbIdParent(const MInt lsId) { return convertIdParent(lsSolver(), lbSolver(), lsId); };
  MInt lb2lsId(const MInt lbId) { return convertId(lbSolver(), lsSolver(), lbId); };
  MInt lb2lsIdParent(const MInt lbId) { return convertIdParent(lbSolver(), lsSolver(), lbId); };

  void testCoupling();

  void updateBoundaryCellsFromGField();

  void createBodyTree();

  void updateLevelSetFlowSolver();
  void updateFlowSolverLevelSet();
  void transferLevelSetFieldValues(MBool);
  void buildCollectedLevelSet(const MInt cellId);
  MFloat interpolateLevelSet(MInt* interpolationCells, MFloat* point, const MInt set);
  void interpolateLsLb(const MInt from, const MInt to);

  MInt noLevelSetFieldData();

  void initializeSolidDomain();

  void constructGField();
  template <MInt bodyType>
  void constructGField_();
  template <MInt bodyType>
  void descendLevelSetValue(const MInt cellId, const MInt* bodyId, const MInt bodyCnt);
};

#endif // ifndef LSLB_H_
