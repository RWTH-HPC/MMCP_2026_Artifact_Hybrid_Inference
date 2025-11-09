// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FCGRIDBNDRYCELL_H
#define FCGRIDBNDRYCELL_H

#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "UTIL/debug.h"

/** /brief This class contains the necessary data to define
 *        a boundary cell for the FC method
 */

template <MInt nDim>
class FcGridBndryCell {
 public:
  MInt m_cellId;

  std::vector<std::vector<MFloat>> m_cutFaces;
  std::vector<std::vector<MFloat>> m_faceNormals;
  std::vector<MInt> m_segmentIdOfCutFace;

  MFloat m_avgFaceNormal[nDim];
  MFloat m_directionOfAction[nDim];
  MFloat m_displacement[nDim];
  MFloat m_load[nDim];

  std::vector<MInt> m_segmentId;
  std::vector<MInt> m_bndryCndId;
};


#endif
