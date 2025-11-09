// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVWMSURFACE_H
#define FVWMSURFACE_H

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class FvWMSurface {
 public:
  FvWMSurface(){};
  ~FvWMSurface();

  void init(MInt, MInt, MFloat);

  MBool m_wmHasImgCell = false;

  MInt m_bndryCellId = -1;
  MInt m_bndrySrfcId = -1;

  MFloat m_wmTauW;
  MFloat m_wmMUEWM;
  MFloat m_wmImgVars[5];
  MFloat m_wmUII;
  MFloat m_wmUTAU;
};

#endif
