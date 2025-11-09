// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POINTBOX_DEFINED
#define POINTBOX_DEFINED

#include <cmath>
#include "INCLUDE/maiatypes.h"
/**
 *   Needed for building a K-D-tree
 *     Numercial Recipies in C: The Art of Scientific Computing Thrid Edition
 *     Authors:
 *
 *
 */
template <MInt DIM>
struct Point {
  MFloat x[DIM]{};
  MInt cellId;
  Point(const Point& p) {
    for(MInt i = 0; i < DIM; i++)
      x[i] = p.x[i];
    cellId = p.cellId;
  }
  Point& operator=(const Point& p) {
    for(MInt i = 0; i < DIM; i++)
      x[i] = p.x[i];
    cellId = p.cellId;
    return *this;
  }
  MBool operator==(const Point& p) const {
    for(MInt i = 0; i < DIM; i++)
      if(x[i] != p.x[i]) return false;
    return true;
  }
  explicit Point(MFloat x0 = 0.0, MFloat x1 = 0.0, MFloat x2 = 0.0, MInt id = -1) {
    x[0] = x0;
    cellId = id;
    if(DIM > 1) x[1] = x1;
    if(DIM > 2) x[2] = x2;
    //	if (DIM > 3) throw("Point not implemented for DIM > 3");
  }
};
template <MInt DIM>
struct Box {
  Point<DIM> lo, hi;
  Box() = default;
  Box(const Point<DIM>& mylo, const Point<DIM>& myhi) : lo(mylo), hi(myhi) {}
};
template <MInt DIM>
MFloat dist(const Point<DIM>& p, const Point<DIM>& q) {
  MFloat dd = 0.0;
  for(MInt j = 0; j < DIM; j++) {
    dd += std::pow(q.x[j] - p.x[j], 2);
  }
  return std::sqrt(dd);
}
template <MInt DIM>
MFloat dist(const Box<DIM>& b, const Point<DIM>& p) {
  MFloat dd = 0;
  for(MInt i = 0; i < DIM; i++) {
    if(p.x[i] < b.lo.x[i]) dd += std::pow(p.x[i] - b.lo.x[i], 2);
    if(p.x[i] > b.hi.x[i]) dd += std::pow(p.x[i] - b.hi.x[i], 2);
  }
  return std::sqrt(dd);
}

#endif
