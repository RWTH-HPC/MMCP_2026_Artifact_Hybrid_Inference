// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBPARENTCELL_H
#define LBPARENTCELL_H

#include "INCLUDE/maiatypes.h"

class LbParentCell {
 public:
  MInt m_cellId;
  static MInt m_maxNoCells;

  //  LbParentCell()
  //  {
  //  }
  static void init(MInt /*dimension*/, MInt /*noDistributions*/, MInt maxNoCells) {
    m_maxNoCells = maxNoCells;
    // m_noChildIds = noChildIds;
  };

  //  virtual ~LbParentCell(){};

  static MInt staticElementSize() { return 0; }

  void allocateElements(void*, void*, MInt cellId);
};

#endif
