// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessingdg.h"
//#include "postprocessing.h"
#include "postdata.h"
#include "postprocessing.cpp"

using namespace std;

template <MInt nDim, class SysEqn>
PostProcessingDg<nDim, SysEqn>::PostProcessingDg(MInt postprocessingId_,
                                                 PostData<nDim>* data,
                                                 DgCartesianSolver<nDim, SysEqn>* ppSolver_)
  : PostProcessingInterface(postprocessingId_),
    PostProcessing<nDim, PostProcessingDg<nDim, SysEqn>>(postprocessingId_, data) {
  m_ppSolver = ppSolver_;
}

template class PostProcessingDg<2, DgSysEqnAcousticPerturb<2>>;
template class PostProcessingDg<3, DgSysEqnAcousticPerturb<3>>;
template class PostProcessingDg<2, DgSysEqnLinearScalarAdv<2>>;
template class PostProcessingDg<3, DgSysEqnLinearScalarAdv<3>>;
