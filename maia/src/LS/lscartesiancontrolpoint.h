// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSCONTROLPOINT_H
#define LSCONTROLPOINT_H

#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class Geometry;

template <MInt nDim>
class LsControlPoint {
 private:
  // Control Point Method 1, Lagrangian point
  MFloat tr[3][3];
  MFloat trinv[3][3];
  // Control Point Method 2, Moving STL
  Geometry<nDim>* geometry;
  MInt noElements;

 protected:
  // Auxilary funtions
  MBool MATXINVT(MFloat T[3][3], MFloat** R);
  MBool LUDECMP(MFloat** a, MInt* indx, MFloat& d);
  void LUBKSB(MFloat** a, MInt* indx, MFloat* b);
  MFloat norm(MFloat* r);
  // Control Point Method 1, Lagrangian point
  void CtrlPnt1_UpdateTR();

 public:
  LsControlPoint();
  ~LsControlPoint();

  // Control Point Method 1, Lagrangian point
  // CtrlPnt1_LagrangPnt[i][j]
  // i=0, mid points; i=1, head in u; i=2, head in v; i=3, head in w
  // Always consider everything as 3D. In 2D, set the 3rd dimension to zero
  MFloat CtrlPnt1_LagrangPnt[4][3];
  void CtrlPnt1_Initialize(Geometry<nDim>* g, MInt dimensions);
  void CtrlPnt1_InitPosition(MFloat* p);
  // Warning!!! orientation about current reference origin
  void CtrlPnt1_InitOrientation(MFloat* u, MFloat* v, MFloat* w);
  // Waring!!! vector q will be changed
  void CtrlPnt1_quvw(MFloat* q);
  // Warning!!!! No orientation implemented, only moving origin
  void CtrlPnt1_CtrlPntToSTL(const MChar* FILNAME);


  // Control Point Method 2, Movable STL
  // pointer to geometry elements
  void CtrlPnt2_shiftSTL(MInt bcId, MFloat* dx);
  void CtrlPnt2_rotateSTL(MInt bcId, MFloat* dphi, MFloat* ori);
  void CtrlPnt2_Initialize(Geometry<nDim>* g, MInt dimensions);
  void CtrlPnt2_InitPosition(MFloat* p);
  // Warning!!! orientation about current reference origin
  void CtrlPnt2_InitOrientation(MFloat* u, MFloat* v, MFloat* w);
  void CtrlPnt2_UpdateAllNormalVector();
  void CtrlPnt2_UpdateNormalVector(MInt e);
  void CtrlPnt2_Update();
  void CtrlPnt2_InitPosition(MFloat*, MInt);
  void CtrlPnt2_InitOrientation(MFloat*, MFloat*, MFloat*, MInt);
  void CtrlPnt2_CtrlPntToSTL(const MChar* FILENAME, MInt bcId);
  void CtrlPnt2_CtrlPntToGNUPlot(const MChar* FILENAME);

 private:
  void CtrlPnt2_MoveMBElementVertex(MInt e, MInt v, MFloat* dx);
  void CtrlPnt2_MoveMBElementVertex(MInt, MInt, MFloat*, MInt);
  MInt CtrlPnt2_getNoElements() { return noElements; }
  void CtrlPnt2_getMBElementVertex(MInt e, MInt v, MFloat* r);
};

#endif
