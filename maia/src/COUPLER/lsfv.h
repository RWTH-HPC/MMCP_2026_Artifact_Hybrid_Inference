// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSFV_H_
#define LSFV_H_

#include "LS/lscartesiancellproperties.h"
#include "LS/lscartesiansolver.h"
#include "cartesiansolver.h"
#include "coupling.h"
#include "couplingutils.h"
#include "enums.h"

// Workaround to avoid adding this->func() for every function call to the base class
template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

// template <MInt nDim> class CouplingLS;
template <MInt nDim, class SysEqn>
class CouplingFv;


template <MInt nDim_, class SysEqn>
class CouplingLsFv : public CouplingLS<nDim_>, CouplingFv<nDim_, SysEqn> {
 public:
  static constexpr MInt nDim = nDim_;

 private:
  friend class FvCartesianSolverXD<nDim, SysEqn>;
  friend class LsCartesianSolver<nDim>;

 public:
  // Simplifications
  using LsSolver = LsCartesianSolver<nDim>;
  // using FvCartesianSolver = typename maia::lsfv::FvCartesianSolverXD<nDim, SysEqn>::type;
  using FvCartesianSolver = FvCartesianSolverXD<nDim, SysEqn>;

  // Type for cell properties
  // FvCell is already set in one of the headers and will be used as such here!
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  // Constructor
  CouplingLsFv<nDim_, SysEqn>(const MInt couplingId, LsSolver* ls, FvCartesianSolver* fv);

  using CouplingLS<nDim_>::lsSolver;
  using CouplingLS<nDim_>::a_levelSetFunctionG;
  using CouplingLS<nDim_>::a_noSets;
  using CouplingLS<nDim_>::a_noLsCells;
  using CouplingLS<nDim_>::a_outsideGValue;
  using CouplingLS<nDim_>::a_bodyIdG;
  using CouplingLS<nDim_>::a_bodyToSet;
  using CouplingLS<nDim_>::a_noEmbeddedBodies;
  using CouplingLS<nDim_>::a_potentialGapCellClose;
  using CouplingLS<nDim_>::a_nearGapG;

  using CouplingFv<nDim, SysEqn>::fvSolver;
  using CouplingFv<nDim, SysEqn>::noSolvers;

 private:
  MInt m_bandWidthRef;
  MInt m_bandWidthRefMax;

  MInt m_fvSolverId;
  MInt m_lsSolverId;

  MInt m_G0regionId;
  MFloat m_initialCrankAngle;

  MInt* m_hadGapCells{};
  MInt* m_hasGapCells{};

  MFloat m_cfl;
  MFloat m_maxVelocity;
  MInt m_timeStepMethod;
  MString m_solverMethod;

  static constexpr const MInt m_maxNoEmbeddedBodies = 20;

  // ... computeBodyProperties
  MBool m_static_computeBodyProperties_first = true;
  MFloat m_static_computeBodyProperties_amplitude[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_freqFactor[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_initialBodyCenter[m_maxNoEmbeddedBodies * 3]{};
  MFloat m_static_computeBodyProperties_Strouhal{};
  MFloat m_static_computeBodyProperties_mu[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_mu2[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftStartAngle1[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftEndAngle1[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftStartAngle2[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftEndAngle2[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_circleStartAngle[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_normal[m_maxNoEmbeddedBodies * 3]{};
  MInt m_static_computeBodyProperties_bodyToFunction[m_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_omega{};
  MFloat m_static_computeBodyProperties_rotAngle{};

  // crankAngle
  MBool m_static_crankAngle_first = true;
  MFloat m_static_crankAngle_Strouhal;


  //----------------------------------------- functions --------------------------------------------

  void initData();
  void checkProperties() override;
  void readProperties() override;

 public:
  void init() override;
  void preCouple(MInt) override;
  void subCouple(MInt /*recipeStep = 0*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override{};
  void postCouple(MInt recipeStep = 0) override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt /*solverId*/) override{};
  void postAdaptation() override;
  void cleanUp() override{};

 private:
  void returnStep_semiLagrange();
  void transferGapCellProperty();
  void testGapProperty();
  void computeGCellTimeStep();
  void computeBodyProperties(MInt returnMode, MFloat* bodyData, MInt body, MFloat time);
  void testCoupling();
  void setLsInList(MIntScratchSpace&);
  void transferLevelSetValues();

  // Id conversion
  MInt ls2fvId(const MInt lsId) { return convertId(lsSolver(), fvSolver(), lsId); };
  MInt ls2fvIdParent(const MInt lsId) { return convertIdParent(lsSolver(), fvSolver(), lsId); };
  MInt fv2lsId(const MInt fvId) {
    // TODO labels:COUPLER,FVMB Is this check needed? This functions should never be called for these modes i guess..
    if(fvSolver().m_levelSetMb && fvSolver().m_constructGField) {
      return -1;
    }
    return convertId(fvSolver(), lsSolver(), fvId);
  };
  MInt fv2lsIdParent(const MInt fvId) {
    // TODO labels:COUPLER,FVMB Is this check needed? This functions should never be called for these modes i guess..
    if(fvSolver().m_levelSetMb && fvSolver().m_constructGField) {
      return -1;
    }
    return convertIdParent(fvSolver(), lsSolver(), fvId);
  };

  MInt noLevelSetFieldData();
  MFloat interpolateLevelSet(MInt cellId, MFloat* point, MInt set);

  MBool returnStep();
  MFloat lsTimeStep() const { return fvSolver().m_timeStep; }
  MFloat crankAngle(MFloat);

  MFloat a_meanCoord(MInt dir) const { return fvSolver().m_meanCoord[dir]; }

  MFloat a_UInfinity() const { return fvSolver().m_UInfinity; }
  MFloat a_TInfinity() const { return fvSolver().m_TInfinity; }
  MFloat a_Ma() const { return fvSolver().a_Ma(); }

  MFloat a_time() const { return fvSolver().m_time; }
  MFloat a_timeRef() const { return fvSolver().m_timeRef; }
  MFloat a_physicalTime() const { return fvSolver().m_physicalTime; }
  MInt a_noFvCells() const { return fvSolver().a_noCells(); }
  MInt a_noFvGridCells() const { return fvSolver().c_noCells(); }
  MInt a_G0CellId(MInt id, MInt set) const { return lsSolver().a_G0CellId(id, set); }
  MInt a_RKStep() const { return fvSolver().m_RKStep; }
  MInt a_noRKSteps() const { return fvSolver().m_noRKSteps; }
};

#endif // ifndef LSFV_H_
