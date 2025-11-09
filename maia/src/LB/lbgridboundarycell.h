// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBGRIDBOUNDARYCELL_H
#define LBGRIDBOUNDARYCELL_H

#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "UTIL/debug.h"

/** \brief This class contains the necessary data to define a boundary cell for
 *  the LB method
 */
template <MInt nDim>
struct LbGridBoundaryCell {
  MInt m_cellId;
  std::vector<MFloat> m_distances;
  MFloat m_multiplier;
  MFloat m_eta;
  MBool m_isFluid;
  std::vector<MInt> m_segmentId;
  std::vector<MInt> m_bndCndId;
};

#endif
