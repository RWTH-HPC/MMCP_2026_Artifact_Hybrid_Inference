// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbparentcell.h"

MInt LbParentCell::m_maxNoCells;

void LbParentCell::allocateElements(void* /*memPointer*/, void* /*basePointer*/, MInt /*cellId*/) {
  //  MFloat * tmpPointer = (MFloat *)memPointer;
  // MFloat * tmpPointerDataBase = (MFloat *)basePointer;
  // MFloat * tmpPointer;

  //   tmpPointer = (MFloat*)((MInt*) tmpPointerDataBase +  m_noInterpolationNeighbors * cellId);
  //   m_interpolationNeighbors = (MInt*) tmpPointer;
  //   tmpPointerDataBase = (MFloat*)((MInt*)tmpPointerDataBase + m_noInterpolationNeighbors * m_maxNoCells);

  //   tmpPointer = (MFloat*)((MFloat*) tmpPointerDataBase +  m_noInterpolationNeighbors * cellId);
  //   m_interpolationCoefficients = tmpPointer;
  //   tmpPointerDataBase = (MFloat*)((MFloat*) tmpPointerDataBase +  m_noInterpolationNeighbors * m_maxNoCells);
}
