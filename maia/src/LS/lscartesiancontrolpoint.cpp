// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lscartesiancontrolpoint.h"

#include <fstream>
#include <iomanip>
#include "GEOM/geometry.h"
#include "GEOM/geometryelement.h"
#include "INCLUDE/maiaconstants.h"
#include "MEMORY/alloc.h"
#include "MEMORY/collector.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/timer.h"

using namespace std;

template <MInt nDim>
LsControlPoint<nDim>::LsControlPoint() {
  TRACE();
  geometry = nullptr;
}

template <MInt nDim>
LsControlPoint<nDim>::~LsControlPoint() {
  TRACE();
  geometry = nullptr;
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_Initialize(Geometry<nDim>* g, MInt NotUsed(dimensions)) {
  TRACE();

  geometry = g;
  noElements = geometry->m_mbelements->size();
  // initialize geometry at ref. geometry mid-point
  for(MInt j = 0; j < nDim; j++)
    CtrlPnt1_LagrangPnt[0][j] = 0.0;
  // the other 3 points at bounding box maximum in the direction
  for(MInt i = 0; i < 4; i++)
    for(MInt j = 0; j < 3; j++)
      CtrlPnt1_LagrangPnt[i][j] = 0.0;
  for(MInt i = 1; i <= nDim; i++) {
    for(MInt j = 0; j < nDim; j++) {
      if(j == (i - 1)) CtrlPnt1_LagrangPnt[i][j] = geometry->m_mbminMax[j + nDim];
    }
  }
  // Special pre-defined for 2D case
  IF_CONSTEXPR(nDim == 2) CtrlPnt1_LagrangPnt[3][2] = 1.0;
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_InitPosition(MFloat* p) {
  TRACE();

  for(MInt i = 0; i <= 3; i++)
    for(MInt j = 0; j < 3; j++) {
      CtrlPnt1_LagrangPnt[i][j] += p[j];
    }
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_InitOrientation(MFloat* u, MFloat* v, MFloat* w) {
  TRACE();

  // Transpose already
  for(MInt i = 0; i < 3; i++) {
    tr[0][i] = u[i];
    tr[1][i] = v[i];
    tr[2][i] = w[i];
  }

  MFloat** R = nullptr;
  mAlloc(R, 3, 3, "R", F0, AT_);
  if(!MATXINVT(tr, R)) {
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++)
        R[i][j] = tr[i][j];
  }
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      trinv[i][j] = R[i][j];
  for(MInt i = 0; i < 4; i++) {
    MFloat quvw[3] = {0, 0, 0}, qxyz[3] = {0, 0, 0};
    for(MInt j = 0; j < 3; j++)
      quvw[j] = CtrlPnt1_LagrangPnt[i][j];
    for(MInt k = 0; k < 3; k++)
      for(MInt j = 0; j < 3; j++) {
        qxyz[k] += trinv[k][j] * quvw[j];
      }
    for(MInt j = 0; j < 3; j++)
      CtrlPnt1_LagrangPnt[i][j] = qxyz[j];
  }

  mDeallocate(R);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_UpdateTR() {
  TRACE();

  MFloat u[3], v[3], w[3];
  for(MInt i = 0; i < 3; i++) {
    u[i] = CtrlPnt1_LagrangPnt[1][i] - CtrlPnt1_LagrangPnt[0][i];
    v[i] = CtrlPnt1_LagrangPnt[2][i] - CtrlPnt1_LagrangPnt[0][i];
    w[i] = CtrlPnt1_LagrangPnt[3][i] - CtrlPnt1_LagrangPnt[0][i];
  }
  MFloat norm_u = norm(u);
  MFloat norm_v = norm(v);
  MFloat norm_w = norm(w);
  for(MInt i = 0; i < 3; i++) {
    u[i] = u[i] / norm_u;
    v[i] = v[i] / norm_v;
    w[i] = w[i] / norm_w;
  }
  for(MInt j = 0; j < 3; j++) {
    tr[0][j] = u[j];
    tr[1][j] = v[j];
    tr[2][j] = w[j];
  }
}

// Warning vector q will be changed
template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_quvw(MFloat* q) {
  TRACE();

  CtrlPnt1_UpdateTR();
  MFloat quvw[3] = {0, 0, 0};
  for(MInt i = 0; i < 3; i++)
    q[i] = q[i] - CtrlPnt1_LagrangPnt[0][i];
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++) {
      quvw[i] += tr[i][j] * q[j];
    }
  for(MInt i = 0; i < 3; i++)
    q[i] = quvw[i];
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt1_CtrlPntToSTL(const MChar* FILENAME) {
  TRACE();

  IF_CONSTEXPR(nDim == 2) {
    cerr << "CtrlPntToSTL is not valid for 2D." << endl;
    return;
  }
  // Update TR
  CtrlPnt1_UpdateTR();
  // Find invert TR
  MFloat** R = nullptr;
  mAlloc(R, 3, 3, "R", F0, AT_);

  if(!MATXINVT(tr, R)) {
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++)
        R[i][j] = tr[i][j];
  }
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      trinv[i][j] = R[i][j];

  ofstream fo(FILENAME);
  fo << "solid CTRL-PNT" << endl;
  fo << scientific;
  for(MInt e = 0; e < noElements; e++) {
    fo << setw(15) << "facet normal";
    MFloat nuvw[3] = {0, 0, 0}, nxyz[3] = {0, 0, 0};
    for(MInt i = 0; i < 3; i++) {
      nuvw[i] = geometry->mbelements[e].m_normal[i];
    }
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++) {
        nxyz[i] += trinv[i][j] * nuvw[j];
      }
    for(MInt i = 0; i < 3; i++) {
      fo << setw(20) << nxyz[i] << '\t';
    }
    fo << endl;
    fo << setw(15) << "outer loop" << endl;
    for(MInt v = 0; v < 3; v++) {
      MFloat quvw[3] = {0, 0, 0}, qxyz[3] = {0, 0, 0};
      fo << setw(15) << "vertex";
      for(MInt i = 0; i < 3; i++)
        quvw[i] = geometry->mbelements[e].m_vertices[v][i];
      for(MInt i = 0; i < 3; i++)
        for(MInt j = 0; j < 3; j++) {
          qxyz[i] += trinv[i][j] * quvw[j];
        }
      for(MInt i = 0; i < 3; i++)
        qxyz[i] += CtrlPnt1_LagrangPnt[0][i];
      for(MInt i = 0; i < 3; i++) {
        fo << setw(20) << qxyz[i];
      }
      fo << endl;
    }
    fo << setw(15) << "endloop" << endl;
    fo << setw(15) << "endfacet" << endl;
  }
  fo << "endsolid" << endl;
  fo.close();

  mDeallocate(R);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_Initialize(Geometry<nDim>* g, MInt NotUsed(dimensions)) {
  TRACE();

  geometry = g;
  noElements = geometry->m_mbelements->size();
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_InitPosition(MFloat* p) {
  TRACE();

  geometry->MoveAllMBElementVertex(p);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_InitOrientation(MFloat* u, MFloat* v, MFloat* w) {
  TRACE();

  // Transpose already
  for(MInt i = 0; i < 3; i++) {
    tr[0][i] = u[i];
    tr[1][i] = v[i];
    tr[2][i] = w[i];
  }
  MFloat** R = nullptr;
  mAlloc(R, 3, 3, "R", F0, AT_);

  if(!MATXINVT(tr, R)) {
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++)
        R[i][j] = tr[i][j];
  }
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      trinv[i][j] = R[i][j];
  mDeallocate(R);
  for(MInt e = 0; e < noElements; e++) {
    for(MInt vt = 0; vt < nDim; vt++) {
      MFloat quvw[3] = {0, 0, 0}, qxyz[3] = {0, 0, 0};
      for(MInt i = 0; i < nDim; i++)
        quvw[i] = geometry->mbelements[e].m_vertices[vt][i];
      for(MInt i = 0; i < 3; i++)
        for(MInt j = 0; j < 3; j++) {
          qxyz[i] += trinv[i][j] * quvw[j];
        }
      geometry->ReplaceMBElementVertex(e, vt, qxyz);
    }
  }
  CtrlPnt2_UpdateAllNormalVector();
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_getMBElementVertex(MInt e, MInt v, MFloat* r) {
  TRACE();

  for(MInt i = 0; i < nDim; i++)
    r[i] = geometry->mbelements[e].m_vertices[v][i];
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_MoveMBElementVertex(MInt e, MInt v, MFloat* dx) {
  TRACE();

  geometry->MoveMBElementVertex(e, v, dx);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_Update() {
  TRACE();

  geometry->UpdateMBBoundingBox();
  geometry->UpdateADT();
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_UpdateAllNormalVector() {
  TRACE();

  for(MInt e = 0; e < noElements; e++)
    geometry->UpdateMBNormalVector(e);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_UpdateNormalVector(MInt e) {
  TRACE();

  geometry->UpdateMBNormalVector(e);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_InitPosition(MFloat* p, MInt bcId) {
  TRACE();

  for(MInt e = 0; e < noElements; e++)
    if(geometry->mbelements[e].m_bndCndId == bcId)
      for(MInt v = 0; v < nDim; v++)
        geometry->MoveMBElementVertex(e, v, p);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_InitOrientation(MFloat* u, MFloat* v, MFloat* w, MInt bcId) {
  TRACE();

  // Transpose already
  for(MInt i = 0; i < 3; i++) {
    tr[0][i] = u[i];
    tr[1][i] = v[i];
    tr[2][i] = w[i];
  }
  MFloat** R = nullptr;
  mAlloc(R, 3, 3, "R", F0, AT_);

  if(!MATXINVT(tr, R)) {
    for(MInt i = 0; i < 3; i++)
      for(MInt j = 0; j < 3; j++)
        R[i][j] = tr[i][j];
  }
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      trinv[i][j] = R[i][j];
  // for(MInt i=0;i<3;i++) {delete [] R[i];} delete [] R; R=nullptr;
  mDeallocate(R);
  for(MInt e = 0; e < noElements; e++) {
    if(geometry->mbelements[e].m_bndCndId == bcId) {
      for(MInt vt = 0; vt < nDim; vt++) {
        MFloat quvw[3] = {0, 0, 0}, qxyz[3] = {0, 0, 0};
        for(MInt i = 0; i < nDim; i++)
          quvw[i] = geometry->mbelements[e].m_vertices[vt][i];
        for(MInt i = 0; i < 3; i++)
          for(MInt j = 0; j < 3; j++) {
            qxyz[i] += trinv[i][j] * quvw[j];
          }
        geometry->ReplaceMBElementVertex(e, vt, qxyz);
      }
    }
  }
  CtrlPnt2_UpdateAllNormalVector();
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_MoveMBElementVertex(MInt e, MInt v, MFloat* dx, MInt bcId) {
  TRACE();

  if(geometry->mbelements[e].m_bndCndId == bcId) geometry->MoveMBElementVertex(e, v, dx);
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_shiftSTL(MInt bcId, MFloat* dx) {
  TRACE();

  // Movable STL
  for(MInt e = 0; e < CtrlPnt2_getNoElements(); e++) {
    for(MInt v = 0; v < nDim; v++) {
      // Move control points to the new time
      CtrlPnt2_MoveMBElementVertex(e, v, dx, bcId);
    }
  }
  CtrlPnt2_Update();
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_rotateSTL(MInt bcId, MFloat* dphi, MFloat* ori) {
  TRACE();

  IF_CONSTEXPR(nDim != 3) {
    mTerm(1, AT_, "Untested in 2D!");
    // nPos[0] = ori[0] - cos(dphi[0]) * (oPos[0] - ori[0]) + sin(dphi[0]) * (oPos[0] - ori[0]);
    // nPos[1] = ori[1] - sin(dphi[0]) * (oPos[0] - ori[0]) - cos(dphi[0]) * (oPos[0] - ori[0]);
  }

  MFloat oPos[nDim];
  MFloat nPos[nDim];
  MFloat dx[nDim];

  MFloat rotMat[3][3];
  rotMat[0][0] = (cos(dphi[1]) * cos(dphi[0]));
  rotMat[0][1] = (cos(dphi[1]) * sin(dphi[0]));
  rotMat[0][2] = (-sin(dphi[1]));
  rotMat[1][0] = (sin(dphi[2]) * sin(dphi[1]) * cos(dphi[0]) - cos(dphi[2]) * sin(dphi[0]));
  rotMat[1][1] = (sin(dphi[2]) * sin(dphi[1]) * sin(dphi[0]) + cos(dphi[2]) * cos(dphi[0]));
  rotMat[1][2] = (sin(dphi[2]) * cos(dphi[1]));
  rotMat[2][0] = (cos(dphi[2]) * sin(dphi[1]) * cos(dphi[0]) + sin(dphi[2]) * sin(dphi[0]));
  rotMat[2][1] = (cos(dphi[2]) * sin(dphi[1]) * sin(dphi[0]) - sin(dphi[2]) * cos(dphi[0]));
  rotMat[2][2] = (cos(dphi[2]) * cos(dphi[1]));

  for(MInt e = 0; e < CtrlPnt2_getNoElements(); e++) {
    if(geometry->mbelements[e].m_bndCndId != bcId) continue;

    for(MInt v = 0; v < nDim; v++) {
      for(MInt n = 0; n < nDim; n++) {
        oPos[n] = geometry->mbelements[e].m_vertices[v][n];
      }

      nPos[0] = ori[0] + rotMat[0][0] * (oPos[0] - ori[0]) + rotMat[0][1] * (oPos[1] - ori[1])
                + rotMat[0][2] * (oPos[2] - ori[2]);
      nPos[1] = ori[1] + rotMat[1][0] * (oPos[0] - ori[0]) + rotMat[1][1] * (oPos[1] - ori[1])
                + rotMat[1][2] * (oPos[2] - ori[2]);
      nPos[2] = ori[2] + rotMat[2][0] * (oPos[0] - ori[0]) + rotMat[2][1] * (oPos[1] - ori[1])
                + rotMat[2][2] * (oPos[2] - ori[2]);

      for(MInt n = 0; n < nDim; n++) {
        dx[n] = nPos[n] - oPos[n];
      }

      // Move control points to the new time
      CtrlPnt2_MoveMBElementVertex(e, v, dx, bcId);
    }
    geometry->UpdateMBNormalVector(e);
  }
  CtrlPnt2_Update();
}

template <MInt nDim>
MBool LsControlPoint<nDim>::MATXINVT(MFloat T[3][3], MFloat** R) {
  TRACE();

  MFloat** a = nullptr;
  MFloat** y = nullptr;
  MFloat col[3];
  MInt indx[3];
  MFloat d;

  mAlloc(a, 3, 3, "a", F0, AT_);
  mAlloc(y, 3, 3, "y", F0, AT_);

  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      a[i][j] = T[i][j];
  if(!LUDECMP(a, indx, d)) {
    return false;
  }
  for(MInt j = 0; j < 3; j++) {
    for(MInt i = 0; i < 3; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    LUBKSB(a, indx, col);
    for(MInt i = 0; i < 3; i++)
      y[i][j] = col[i];
  }
  for(MInt i = 0; i < 3; i++)
    for(MInt j = 0; j < 3; j++)
      R[i][j] = y[i][j];
  mDeallocate(a);
  mDeallocate(y);

  return true;
}

template <MInt nDim>
MBool LsControlPoint<nDim>::LUDECMP(MFloat** a, MInt* indx, MFloat& d) {
  TRACE();

  const MFloat TINY = 1.0e-20;
  MInt imax = std::numeric_limits<MInt>::min();
  // matrix size
  static const MInt n = 3;
  // scale factor
  MFloat vv[3];
  for(MInt i = 0; i < n; i++) {
    MFloat big = 0.0;
    for(MInt j = 0; j < n; j++) {
      MFloat temp = fabs(a[i][j]);
      if(temp > big) big = temp;
    }
    if(approx(big, 0.0, MFloatEps)) {
      return false;
    }
    vv[i] = 1.0 / big;
  }
  for(MInt j = 0; j < n; j++) {
    for(MInt i = 0; i < j; i++) {
      MFloat sum = a[i][j];
      for(MInt k = 0; k < i; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    MFloat big = 0.0;
    for(MInt i = j; i < n; i++) {
      MFloat sum = a[i][j];
      for(MInt k = 0; k < j; k++)
        sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      MFloat dum = vv[i] * fabs(sum);
      if(dum >= big) {
        big = dum;
        imax = i;
      }
    }
    if(j != imax) {
      for(MInt k = 0; k < n; k++) {
        MFloat dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(approx(a[j][j], 0.0, MFloatEps)) a[j][j] = TINY;
    if(j != n - 1) {
      MFloat dum = 1.0 / (a[j][j]);
      for(MInt i = j + 1; i < n; i++)
        a[i][j] *= dum;
    }
  }

  return true;
}

template <MInt nDim>
void LsControlPoint<nDim>::LUBKSB(MFloat** a, MInt* indx, MFloat* b) {
  TRACE();

  MInt i, ii = 0, ip, j;
  MFloat sum;
  MInt n = 3;
  for(i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if(ii != 0)
      for(j = ii - 1; j < i; j++)
        sum -= a[i][j] * b[j];
    else if(!approx(sum, 0.0, MFloatEps))
      ii = i + 1;
    b[i] = sum;
  }
  for(i = n - 1; i >= 0; i--) {
    sum = b[i];
    for(j = i + 1; j < n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_CtrlPntToSTL(const MChar* FILENAME, MInt bcId) {
  TRACE();

  IF_CONSTEXPR(nDim == 2) {
    MFloat middle[3] = {0.0, 0.0, 0.0};
    ofstream fo(FILENAME);
    MInt noElem = 0;
    // compute middle
    for(MInt i = 0; i < 2; i++) {
      middle[i] = F0;
    }
    for(MInt e = 0; e < noElements; e++) {
      if(geometry->mbelements[e].m_bndCndId == bcId) {
        noElem++;
        for(MInt v = 0; v < 2; v++) {
          for(MInt i = 0; i < 2; i++) {
            middle[i] += geometry->mbelements[e].m_vertices[0][i];
          }
        }
      }
    }
    for(MInt i = 0; i < 2; i++) {
      middle[i] = middle[i] / (noElem * 2);
    }

    fo << "solid CTRL-PNT" << endl;
    fo << scientific;
    for(MInt e = 0; e < noElements; e++) {
      if(geometry->mbelements[e].m_bndCndId == bcId) {
        fo << setw(15) << "facet normal";
        for(MInt i = 0; i < 2; i++) {
          fo << setw(20) << 0.0;
        }
        fo << setw(20) << 1.0;
        fo << endl;
        fo << setw(15) << "outer loop" << endl;
        for(MInt v = 0; v < 2; v++) {
          fo << setw(15) << "vertex";
          for(MInt i = 0; i < 2; i++) {
            fo << setw(20) << geometry->mbelements[e].m_vertices[v][i];
          }
          fo << setw(20) << F0;
          fo << endl;
        }
        fo << setw(15) << "vertex";
        for(MInt i = 0; i < 3; i++) {
          fo << setw(20) << middle[i];
        }
        fo << endl;
        fo << setw(15) << "endloop" << endl;
        fo << setw(15) << "endfacet" << endl;
      }
    }
    fo << "endsolid" << endl;
    fo.close();
  }
  else {
    ofstream fo(FILENAME);
    fo << "solid CTRL-PNT" << endl;
    fo << scientific;
    for(MInt e = 0; e < noElements; e++) {
      if(geometry->mbelements[e].m_bndCndId == bcId) {
        fo << setw(15) << "facet normal";
        for(MInt i = 0; i < 3; i++) {
          fo << setw(20) << geometry->mbelements[e].m_normal[i];
        }
        fo << endl;
        fo << setw(15) << "outer loop" << endl;
        for(MInt v = 0; v < 3; v++) {
          fo << setw(15) << "vertex";
          for(MInt i = 0; i < 3; i++) {
            fo << setw(20) << geometry->mbelements[e].m_vertices[v][i];
          }
          fo << endl;
        }
        fo << setw(15) << "endloop" << endl;
        fo << setw(15) << "endfacet" << endl;
      }
    }
    fo << "endsolid" << endl;
    fo.close();
  }
}

template <MInt nDim>
void LsControlPoint<nDim>::CtrlPnt2_CtrlPntToGNUPlot(const MChar* FILENAME) {
  TRACE();

  ofstream fo(FILENAME);
  for(MInt e = 0; e < noElements; e++) {
    for(MInt v = 0; v < nDim; v++) {
      for(MInt i = 0; i < nDim; i++) {
        fo << geometry->mbelements[e].m_vertices[v][i] << '\t';
      }
      fo << endl;
    }
  }
  fo.close();
}

template <MInt nDim>
MFloat LsControlPoint<nDim>::norm(MFloat* r) {
  return sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));
}

// Explicit instantiations for 2D and 3D
template class LsControlPoint<2>;
template class LsControlPoint<3>;
