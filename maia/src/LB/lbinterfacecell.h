// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBINTERFACECELL_H
#define LBINTERFACECELL_H

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/debug.h"

class LbInterfaceCell {
 public:
  LbInterfaceCell() {}
  static void init(MInt dimension, MInt /*noDistributions*/, MInt maxNoCells) {
    //  TRACE();
    m_noInterpolationNeighbors = IPOW2(dimension);
    m_maxNoCells = maxNoCells;
  };

  virtual ~LbInterfaceCell(){};


 public:
  MInt m_cellId;
  MInt m_position;
  // Should not be static for complex interpolation rules (with different neighbor relations)
  static MInt m_noInterpolationNeighbors;
  static MInt m_maxNoCells;

  MInt* m_interpolationNeighbors;
  MFloat* m_interpolationCoefficients;

  static MInt staticElementSize() {
    return (m_noInterpolationNeighbors) * sizeof(MFloat) + (m_noInterpolationNeighbors) * sizeof(MInt);
  }

  void allocateElements(void*, void*, MInt cellId);
};

#endif
