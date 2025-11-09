// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGFVLPT_H_
#define POSTPROCESSINGFVLPT_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "FV/fvransmodelconstants.h"
#include "INCLUDE/maiatypes.h"
#include "LPT/lpt.h"
#include "globals.h"
#include "postdata.h"
#include "postprocessing.h"
#include "postprocessingfv.h"
#include "postprocessinglpt.h"


template <MInt nDim, class SysEqn>
class PostProcessingFvLPT : public PostProcessingLPT<nDim>, PostProcessingFv<nDim, SysEqn> {
 private:
  template <MInt nDim_, class ppType>
  friend class PostProcessing;
  template <MInt nDim_>
  friend class LPT;

 public:
  using lptSolver = LPT<nDim>;
  using fv = FvCartesianSolverXD<nDim, SysEqn>;
  using BaseFv = PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>;
  using BaseLpt = PostProcessing<nDim, PostProcessingLPT<nDim>>;
  using ppFv = PostProcessingFv<nDim, SysEqn>;
  using ppLpt = PostProcessingLPT<nDim>;

  // Constructor
  PostProcessingFvLPT(MInt postprocessingId_, PostData<nDim>* data, fv* ppFvSolver_, lptSolver* ppLptSolver_);

  virtual ~PostProcessingFvLPT(){};

  Solver* mSolver() const override { return static_cast<Solver*>(m_ppSolverFv); };

  static MString s_ppType;

  // void postprocessPreInit() override { BaseFv::postprocessPreInit(); };
  void initPostProcessing() override { BaseFv::initPostProcessing(); };
  void postprocessPreSolve() override { BaseFv::postprocessPreSolve(); };
  void postprocessPostSolve() override { BaseFv::postprocessPostSolve(); };
  void postprocessInSolve(const MBool finalTimeStep) override { BaseFv::postprocessInSolve(finalTimeStep); };
  void postprocessSolution() override { BaseFv::postprocessSolution(); };
  // NOTE: all initalization functions need to be called in both base solvers
  //      if they are required there!
  void initSprayData() override {
    ppFv::initSprayData();
    ppLpt::initSprayData();
  }

 protected:
  void computeSprayData() override;
  void writeSprayData() override;
  void initLPTSolutionFile() override { ppLpt::initLPTSolutionFile(); }
  void writeLPTSolutionFile() override { ppLpt::writeLPTSolutionFile(); };

 private:
  lptSolver* m_ppSolverLpt;
  fv* m_ppSolverFv;

  lptSolver& lpt() const { return *m_ppSolverLpt; }
  fv& fvSolver() const { return *m_ppSolverFv; }

  MInt m_fvLPTSpeciesId = 0;

  void updateData();

  using BaseFv::m_sprayComputeInterval;
  using BaseFv::m_sprayDataSize;
  using BaseFv::m_sprayWriteInterval;

  using BaseFv::m_vapourCV;
  using BaseFv::m_vapourPen;

  using BaseLpt::m_injectionData;
  using BaseLpt::m_LPTSolutionInterval;
  using BaseLpt::m_particleCV;
  using BaseLpt::m_particlePen;
  using BaseLpt::m_sprayStat;
  using BaseLpt::m_sprayTimes;

  using ppLpt::m_cInjData;

  MInt m_sprayDataStep = 0;

  void advanceDataStep() {
    ppFv::advanceDataStep();
    ppLpt::advanceDataStep();
    m_sprayDataStep++;
  };
  void resetDataStep() {
    ppFv::resetDataStep();
    ppLpt::resetDataStep();
    m_sprayDataStep = 0;
  };
};

#endif // POSTPROCESSINGLPT_H_
