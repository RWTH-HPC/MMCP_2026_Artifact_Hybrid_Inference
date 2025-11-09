// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGCONTROLLER_H_
#define POSTPROCESSINGCONTROLLER_H_

#include "postprocessing.h"

template <MInt nDim>
class PostProcessingController {
 public:
  PostProcessingController(std::vector<PostProcessingInterface*> pp);

  void init();
  void preSolve();
  void inSolve(MBool finalTimeStep);
  void ppSolution();
  void postSolve();

  void setStep(const MInt step) {
    m_step = (step < 0) ? 0 : step;
  }

 private:
  std::vector<PostProcessingInterface*> m_PP;

  MInt m_step = -1;
  MInt m_maxNoSteps;
  std::vector<std::vector<MBool>> m_ppOrder{};

  MInt noPP() { return (signed)m_PP.size(); };
};

template <MInt nDim>
PostProcessingController<nDim>::PostProcessingController(std::vector<PostProcessingInterface*> pp) : m_PP(pp) {
  m_maxNoSteps = 1;
  if(Context::propertyExists("recipeMaxNoSteps")) {
    m_maxNoSteps = Context::getBasicProperty<MInt>("recipeMaxNoSteps", AT_, &m_maxNoSteps);
  }
}

template <MInt nDim>
void PostProcessingController<nDim>::init() {
  TRACE();

  m_log << "##################################################################" << std::endl;
  m_log << "##                        Postprocessing                        ##" << std::endl;
  m_log << "##################################################################" << std::endl << std::endl;

  // read callOrder
  m_ppOrder.resize(m_maxNoSteps, std::vector<MBool>(noPP(), true));
  if(m_maxNoSteps > 1) {
    for(MInt p = 0; p < noPP(); p++) {
      const MString propName = "postprocessingOrder_" + std::to_string(p);
      for(MInt step = 0; step < m_maxNoSteps; step++) {
        m_ppOrder[step][p] = (MBool)Context::getBasicProperty<MInt>(propName, AT_, step);
      }
    }
  }

  // init Postprocessing instances
  for(auto&& pp : m_PP) {
    pp->mSolver()->startLoadTimer(AT_);
    pp->initPostProcessing();
    pp->mSolver()->stopLoadTimer(AT_);
  }
}

template <MInt nDim>
void PostProcessingController<nDim>::preSolve() {
  TRACE();

  for(auto&& pp : m_PP) {
    pp->mSolver()->startLoadTimer(AT_);
    pp->postprocessPreSolve();
    pp->mSolver()->stopLoadTimer(AT_);
  }
}


template <MInt nDim>
void PostProcessingController<nDim>::inSolve(MBool finalTimeStep) {
  TRACE();

  for(auto&& pp : m_PP) {
    if(m_ppOrder[m_step][pp->a_postprocessingId()]) {
      pp->mSolver()->startLoadTimer(AT_);
      pp->postprocessInSolve(finalTimeStep);
      pp->mSolver()->stopLoadTimer(AT_);
    }
  }
}

template <MInt nDim>
void PostProcessingController<nDim>::ppSolution() {
  TRACE();

  for(auto&& pp : m_PP) {
    pp->mSolver()->disableDlbTimers();
    pp->postprocessSolution();
    pp->mSolver()->enableDlbTimers();
  }
}

template <MInt nDim>
void PostProcessingController<nDim>::postSolve() {
  TRACE();

  for(auto&& pp : m_PP) {
    pp->mSolver()->startLoadTimer(AT_);
    pp->postprocessPostSolve();
    pp->mSolver()->stopLoadTimer(AT_);
  }
}


#endif // POSTPROCESSINGCONTROLLER_H_
