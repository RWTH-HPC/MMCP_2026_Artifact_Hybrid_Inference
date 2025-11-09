// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGDG_H_
#define POSTPROCESSINGDG_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "DG/dgcartesiansolver.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "DG/dgcartesiansyseqnlinearscalaradv.h"
#include "globals.h"
#include "postdata.h"
#include "postprocessing.h"

template <MInt nDim, class SysEqn>
class PostProcessingDg : public PostProcessing<nDim, PostProcessingDg<nDim, SysEqn>> {
 public:
  using SolverType = DgCartesianSolver<nDim, SysEqn>;
  using Base = PostProcessing<nDim, PostProcessingDg<nDim, SysEqn>>;

  // Constructor
  PostProcessingDg(MInt postprocessingId_, PostData<nDim>* data, SolverType* ppSolver_);

  virtual ~PostProcessingDg(){};

  using Base::initPostProcessing;
  using Base::m_postData;
  using Base::postprocessInSolve;
  using Base::postprocessPostSolve;
  // using Base::postprocessPreInit;
  using Base::postprocessPreSolve;

  SolverType& solver() const { return *m_ppSolver; }

 private:
  SolverType* m_ppSolver;
};

#endif // POSTPROCESSINGDG_H_
