// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGLBLPT_H_
#define POSTPROCESSINGLBLPT_H_


#include "postdata.h"
#include "postprocessing.h"
#include "postprocessinglb.h"
#include "postprocessinglpt.h"

template <MInt nDim>
class PostProcessingLbLPT : public PostProcessingLPT<nDim>, PostProcessingLb<nDim> {
 private:
  template <MInt nDim_, class ppType>
  friend class PostProcessing;
  template <MInt nDim_>
  friend class LPT;

 public:
  using lpt = LPT<nDim>;
  using lb = LbSolver<nDim>;
  using BaseLb = PostProcessing<nDim, PostProcessingLb<nDim>>;
  using BaseLpt = PostProcessing<nDim, PostProcessingLPT<nDim>>;
  using ppLb = PostProcessingLb<nDim>;
  using ppLpt = PostProcessingLPT<nDim>;

  // Constructor
  PostProcessingLbLPT(MInt postprocessingId_, PostData<nDim>* data, lb* ppLbSolver_, lpt* ppLptSolver_);

  virtual ~PostProcessingLbLPT(){};

  Solver* mSolver() const override { return static_cast<Solver*>(m_ppSolverLb); };

  void initPostProcessing() override { BaseLb::initPostProcessing(); };
  void postprocessPreSolve() override { BaseLb::postprocessPreSolve(); };
  void postprocessPostSolve() override { BaseLb::postprocessPostSolve(); };
  void postprocessInSolve(const MBool finalTimeStep) override { BaseLb::postprocessInSolve(finalTimeStep); };
  void postprocessSolution() override { BaseLb::postprocessSolution(); };

 protected:
  void initLPTSolutionFile() override { ppLpt::initLPTSolutionFile(); }
  void writeLPTSolutionFile() override { ppLpt::writeLPTSolutionFile(); };

  void initParticleStatistics() override { ppLpt::initParticleStatistics(); };
  void computeParticleStatistics() override { ppLpt::computeParticleStatistics(); };

  void initPLIsoTurbulenceStatistics() override;
  void computePLIsoTurbulenceStatistics() override;

 private:
  lpt* m_ppSolverLpt;
  lb* m_ppSolverLb;

  MFloat m_conversionLbLptLength;
  MFloat m_conversionLptLbLength;

  lpt& lptSolver() const { return *m_ppSolverLpt; }
  lb& lbSolver() const { return *m_ppSolverLb; }

  using BaseLb::m_finalTimeStep;
  using BaseLb::m_postprocessingId;
  using BaseLpt::m_LPTSolutionInterval;
  using ppLb::m_tau_eta;

  void getConversionFactors();
};

#endif // POSTPROCESSINGLBLPT_H_
