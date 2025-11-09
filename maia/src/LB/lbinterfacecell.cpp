// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbinterfacecell.h"

MInt LbInterfaceCell::m_noInterpolationNeighbors;
MInt LbInterfaceCell::m_maxNoCells;

#if defined(MAIA_CLANG_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
#endif

void LbInterfaceCell::allocateElements(void* /*memPointer*/, void* basePointer, MInt cellId) {
  //  TRACE();
  //  MFloat * tmpPointer = (MFloat *)memPointer;
  MFloat* tmpPointerDataBase = (MFloat*)basePointer;
  MFloat* tmpPointer;

  tmpPointer = (MFloat*)((MInt*)tmpPointerDataBase + m_noInterpolationNeighbors * cellId);
  m_interpolationNeighbors = (MInt*)tmpPointer;
  tmpPointerDataBase = (MFloat*)((MInt*)tmpPointerDataBase + m_noInterpolationNeighbors * m_maxNoCells);

  tmpPointer = tmpPointerDataBase + m_noInterpolationNeighbors * cellId;
  m_interpolationCoefficients = tmpPointer;
  tmpPointerDataBase = tmpPointerDataBase + m_noInterpolationNeighbors * m_maxNoCells;
}

#if defined(MAIA_CLANG_COMPILER)
#pragma GCC diagnostic pop
#endif
