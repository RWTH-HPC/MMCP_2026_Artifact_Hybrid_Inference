// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ELLIPSOIDS_H_
#define ELLIPSOIDS_H_

#include <complex>
#include "INCLUDE/maiatypes.h"

class EllipsoidDistance {
  template <MInt nDim>
  friend class LPT;

 public:
  EllipsoidDistance(MFloat di[3], MFloat l1i[3], MFloat m1i[3], MFloat n1i[3], MFloat l2i[3], MFloat m2i[3],
                    MFloat n2i[3], MFloat a1, MFloat b1, MFloat c1, MFloat a2, MFloat b2, MFloat c2);
  void crossP(MFloat* x, MFloat* y, MFloat* z);                      //
  void norm(MFloat* vec, MFloat* nvec);                              //
  MFloat mag(MFloat* V);                                             //
  MFloat dotP(MFloat* Vec1, MFloat* Vec2);                           // Functions
  MFloat plane_int(MFloat);                                          //
  MFloat distance2d(MFloat, MFloat, MFloat, MFloat, MFloat, MFloat); //
  std::complex<MFloat> c_cbrt(std::complex<MFloat>);                 //
  MFloat ellipsoids(void);
  MFloat l1i[3], l2i[3], m1i[3], m2i[3], n1i[3], n2i[3], di[3];                      // Input vectors
  MFloat a1, a2, b1, b2, c1, c2;                                                     // Input semiaxes lenghts
  MFloat d[3], p0[3], p[3], s[3], l1[3], l2[3], m1[3], m2[3], n1[3], n2[3], dxp0[3]; // Normalized vectors
  MFloat pi;
  MFloat gratio;
  MFloat o1g;
  MFloat tolerance;
  MFloat delt;
};

#endif
