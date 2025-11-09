// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYADTNODE_H_
#define GEOMETRYADTNODE_H_

#include "INCLUDE/maiatypes.h"

/** Part of the ADT implementation
 *
 *  For further details see geometryadt.h
 */
class GeometryAdtNode {
 public:
  MInt m_parent = -1;
  MInt m_leftSubtree = 0;
  MInt m_rightSubtree = 0;
  MInt m_depth = -1;

  MFloat m_partition = MFloatNaN;

  /// Holds the minimum value
  MFloat m_a = MFloatNaN;

  /// Holds the maximum value
  MFloat m_b = MFloatNaN;

  /// Holds the id of the connected element
  MInt m_element = -1;

  MInt getStaticElementSize() { return 5 * sizeof(MInt) + 3 * sizeof(MFloat); };
};

#endif
