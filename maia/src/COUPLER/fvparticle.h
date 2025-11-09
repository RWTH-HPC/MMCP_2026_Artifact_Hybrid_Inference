// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVPARTICLE_H_
#define FVPARTICLE_H_

#include "LPT/lpt.h"
#include "coupling.h"
#include "couplinglpt.h"
#include "couplingutils.h"
#include "solver.h"

// Forward declarations
template <MInt nDim>
class CouplingParticle;

template <MInt nDim, class SysEqn>
class CouplerFvParticle : public CouplingLpt<nDim, CouplingFv<nDim, SysEqn>> {
 private:
  friend class LPT<nDim>;
  friend class FvCartesianSolverXD<nDim, SysEqn>;

  using FvCartesianSolver = FvCartesianSolverXD<nDim, SysEqn>;

  using BaseFv = CouplingFv<nDim, SysEqn>;
  using BaseLpt = CouplingParticle<nDim>;
  using BaseLptX = CouplingLpt<nDim, BaseFv>;

  using BaseFv::returnIdleRecord;
  using BaseFv::returnLoadRecord;
  using BaseFv::startLoadTimer;
  using BaseFv::stopLoadTimer;

  using typename BaseLptX::ConversionFactor;

  ConversionFactor& conversionLptFv = BaseLptX::conversionLptFlow;
  ConversionFactor& conversionFvLpt = BaseLptX::conversionFlowLpt;

 public:
  // Constructor
  CouplerFvParticle<nDim, SysEqn>(const MInt couplingId, LPT<nDim>* particle, FvCartesianSolver* fv);

 protected:
  LPT<nDim>& lpt() const override { return *m_particle; }
  FvCartesianSolverXD<nDim, SysEqn>& fvSolver() const { return *m_fvSolver; }

 public:
  void init() override;
  void finalizeCouplerInit() override;
  void finalizeSubCoupleInit(MInt /*solverId*/) override{};

  void preCouple(MInt) override;
  void postCouple(MInt) override;
  void subCouple(MInt, MInt, std::vector<MBool>&) override;

  void postAdaptation() override{};
  void finalizeAdaptation(const MInt) override;
  void prepareAdaptation() override;

  void finalizeBalance(const MInt unused) override;
  void balancePre() override{};
  void balancePost() override{};

  void cleanUp() override{};

  void writeRestartFile(const MInt) override{};

  MInt noCouplingTimers(const MBool NotUsed(allTimings)) const override { return 2; }

  void getCouplingTimings(std::vector<std::pair<MString, MFloat>>& timings, const MBool NotUsed(allTimings)) override {
    const MString namePrefix = "c" + std::to_string(this->couplerId()) + "_";

    timings.emplace_back(namePrefix + "loadCouplerFvParticle", returnLoadRecord());
    timings.emplace_back(namePrefix + "idleCouplerFvParticle", returnIdleRecord());
  };


 private:
  void checkProperties() override;
  void readProperties();
  void initData();
  void initConversion();
  void initParticleVelocity();

  void transferNoParticlesInCell();
  void writeParticleCellStats();

  void transferExternalSources();
  void updateLPTBndry();
  void transferFlowField();
  void transferCellStatus();
  void transferVelocitySlopes();
  void transferTimeStep();

  void setExternalSourceInCell(const MInt, const MInt, const MFloat);

  void unifyTimeStep();

  MInt childLoop(MInt cellId);

  // Id conversion
  MInt lpt2fvId(const MInt lptId) { return convertId(lpt(), fvSolver(), lptId); };
  MInt lpt2fvIdParent(const MInt lptId) { return convertIdParent(lpt(), fvSolver(), lptId); };
  MInt fv2lptId(const MInt fvId) { return convertId(fvSolver(), lpt(), fvId); };
  MInt fv2lptIdParent(const MInt fvId) { return convertIdParent(fvSolver(), lpt(), fvId); };

  void interpolateFVLPT(const MInt from, const MInt to);
  void interpolateVelocitySlopesFVLPT(const MInt from, const MInt to);
  MFloat interpolateVariable(MInt*, MFloat*, MInt);
  MFloat interpolateSlope(MInt*, MFloat*, MInt, MInt);


  LPT<nDim>* m_particle = nullptr;
  FvCartesianSolverXD<nDim, SysEqn>* m_fvSolver = nullptr;

  MInt m_lptSolverId{};
  MInt m_fvSolverId{};
  MInt m_lptSolverOrder{};
  MInt m_fvSolverOrder{};
  MInt m_noSolverSteps{};

  MBool m_lptFvInterpolation = false;
  MInt m_fvLPTSpeciesId{};

  SysEqn m_sysEqn;

  MInt m_solutionStep;
  MInt m_noSolutionSteps;
  MBool m_forceFvTimeStep;
  MBool m_interLeafed;
};

#endif // ifndef FVPARTICLE_H_
