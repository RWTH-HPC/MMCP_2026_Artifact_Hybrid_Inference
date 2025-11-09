// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianwmsurface.h"

template <MInt nDim>
FvWMSurface<nDim>::~FvWMSurface() {
  return;
}

template <MInt nDim>
void FvWMSurface<nDim>::init(MInt bndryCellId, MInt bndrySrfcId, MFloat utau) {
  m_bndryCellId = bndryCellId;
  m_bndrySrfcId = bndrySrfcId;

  for(MInt v = 0; v < 5; v++) {
    m_wmImgVars[v] = F0;
  }
  m_wmTauW = F0;
  m_wmMUEWM = F0;
  m_wmUII = F0;
  m_wmUTAU = utau;
}

template class FvWMSurface<2>;
template class FvWMSurface<3>;
