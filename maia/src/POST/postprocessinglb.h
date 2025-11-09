// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGLB_H_
#define POSTPROCESSINGLB_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "LB/lbsolver.h"
#include "globals.h"
#include "postprocessing.h"
#include "samplingdata.h"

template <MInt nDim>
class PostProcessingLb : public PostProcessing<nDim, PostProcessingLb<nDim>> {
 private:
  template <MInt nDim_, class ppType>
  friend class PostProcessing;

 public:
  using SolverType = LbSolver<nDim>;
  using Base = PostProcessing<nDim, PostProcessingLb<nDim>>;

  // Constructor
  PostProcessingLb(MInt postprocessingId_, PostData<nDim>* data, SolverType* ppSolver_);

  virtual ~PostProcessingLb();

  SolverType& solver() const { return *m_ppSolver; }

  using Base::initPostProcessing;
  using Base::m_finalTimeStep;
  using Base::m_postData;
  using Base::postData;
  using Base::postprocessInSolve;
  using Base::postprocessPostSolve;
  // using Base::postprocessPreInit;
  using Base::getSampleVariables;
  using Base::m_averageInterval;
  using Base::m_averageStartTimestep;
  using Base::m_averageStopTimestep;
  using Base::m_localVars;
  using Base::m_noVariables;
  using Base::postprocessPreSolve;

 private:
  SolverType* m_ppSolver;

  // Average_inSolve
  std::vector<MString> m_varNames;
  MBool m_needVelocityGradient = false;

 protected:
  void initPointSamplingData() override;
  void savePointSamplingData() override;

  void initSurfaceSamplingData() override;
  void saveSurfaceSamplingData() override;

  void initVolumeSamplingData() override;
  void saveVolumeSamplingData() override;

  void initIsoTurbulenceStatistics() override;
  void computeIsoTurbulenceStatistics() override;

  std::unique_ptr<PointData<nDim, SolverType>> m_pointData;
  std::unique_ptr<SurfaceData<nDim, SolverType>> m_surfaceData;
  std::unique_ptr<VolumeData<nDim, SolverType>> m_volumeData;

  MFloat m_tau_eta;

  MBool getSampleVarsDerivatives(const MInt cellId, std::vector<MFloat>& vars) {
    return solver().getSampleVarsDerivatives(cellId, vars);
  };
};

#endif // POSTPROCESSINGLB_H_
