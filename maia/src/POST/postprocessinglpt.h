// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGLPT_H_
#define POSTPROCESSINGLPT_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "INCLUDE/maiatypes.h"
#include "LPT/lpt.h"
#include "globals.h"
#include "postdata.h"
#include "postprocessing.h"


template <MInt nDim>
class PostProcessingLPT : public PostProcessing<nDim, PostProcessingLPT<nDim>> {
 private:
  template <MInt nDim_, class ppType>
  friend class PostProcessing;
  template <MInt nDim_>
  friend class LPT;

 public:
  using SolverType = LPT<nDim>;
  using Base = PostProcessing<nDim, PostProcessingLPT<nDim>>;
  using Base::m_injectionData;
  using Base::m_LPTSolutionInterval;
  using Base::m_particleCV;
  using Base::m_particlePen;
  using Base::m_postprocessingId;
  using Base::m_sprayDataSize;
  using Base::m_sprayDataStep;
  using Base::m_sprayStat;
  using Base::m_sprayTimes;

  // Constructor
  PostProcessingLPT(MInt postprocessingId_, PostData<nDim>* data, SolverType* ppSolver_);

  virtual ~PostProcessingLPT(){};

  SolverType& solver() const { return *m_ppSolver; }

  void initSprayData() override;

 protected:
  void initLPTSolutionFile() override;
  void writeLPTSolutionFile() override;

  void initParticleStatistics() override;
  void computeParticleStatistics() override;
  SolverType* m_ppSolver;

  void initParticleLog();
  void writeParticleLog();

  void particleMass();
  void parcelStatistics();
  void particlePenetration();
  MInt getInjectionData();

  void advanceDataStep() { m_sprayDataStep++; };
  void resetDataStep() { m_sprayDataStep = 0; };

  MFloat** m_cInjData = nullptr;
  MFloat m_AvgRe_p = F0;
  MFloat m_VpRMS = F0;
  MFloat m_EkinP = F0;

  MBool m_writeParLog = false;
  MInt m_parLogInterval = 50;
  MBool m_parLogApp = false;
};

#endif // POSTPROCESSINGLPT_H_
