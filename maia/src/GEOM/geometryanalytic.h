// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ANALYTICGEOMETRY_H
#define ANALYTICGEOMETRY_H
#include <iostream>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"
#include "UTIL/maiamath.h"

namespace maia {
namespace geom {

template <typename T, MInt nDim>
MBool doBoxesOverlap(const T* const a, const T* const b) {
  for(MInt dim = 0; dim < nDim; dim++) {
    if((a[dim] > b[dim + nDim]) || b[dim] > a[dim + nDim]) {
      return false;
    }
  }
  return true;
}

template <typename T, MInt nDim>
MBool isBoxInsideBox(const T* const a, const T* const b) {
  for(MInt dim = 0; dim < nDim; dim++) {
    if((a[dim] < b[dim]) || a[dim + nDim] > b[dim + nDim]) {
      return false;
    }
  }
  return true;
}

template <typename T, MInt nDim>
MBool isPointInsideBox(const T* const p, const T* const b) {
  for(MInt dim = 0; dim < nDim; dim++) {
    if((p[dim] < b[dim]) || p[dim] > b[dim + nDim]) {
      return false;
    }
  }
  return true;
}

template <MInt nDim>
MBool isPointInsideSphere(const MFloat* const p, const MFloat* const c, const MFloat R) {
  MFloat r = F0;
  for(MInt dim = 0; dim < nDim; dim++) {
    r += POW2(p[dim] - c[dim]);
  }
  return std::sqrt(r) < R;
}


// checks if point p is within a ring segment centered around c. The ring segment is given by r_min and r_max (R),
// phi_min and phi_max (phi) and its rotational axis. Note, that R[0] < R[1] and phi[0] < phi[1] for this to work as
// intended.
template <MInt nDim>
MBool isPointInsideRingSegment(const MFloat* const p, const MFloat* const c, const MFloat* const R,
                               const MFloat* const phi, const MInt axis) {
  IF_CONSTEXPR(nDim == 2) {
    std::cerr << "H-Patch not implemented for 2D, not refining any cells! Take care with the angle calculation in 2D "
                 "when fixing this!"
              << std::endl;
    return false;
  }

  std::array<MFloat, nDim> normalVector;
  MFloat distance = F0;

  // compute normal vector from axis to current point
  for(MInt d = 0; d < nDim; d++) {
    if(d == axis) {
      normalVector[d] = F0;
    } else {
      normalVector[d] = p[d] - c[d];
      distance += normalVector[d] * normalVector[d];
    }
  }

  // check if cell is within r_min and r_max
  distance = sqrt(distance);
  if(distance < R[0] || distance > R[1]) return false;

  // check if cell is within phi_min and phi_max
  MFloat angle = atan2(normalVector[(axis + 2) % nDim], normalVector[(axis + 1) % nDim]) * 180.0 / PI;
  if(angle < phi[0] || angle > phi[1]) return false;

  return true;
}

template <MInt nDim>
MBool isBoxInsideSphere(const MFloat* const p, const MFloat* const c, const MFloat R) {
  MFloat pp[nDim];
  const MInt cornerCode[8][3] = {{0, 0, 0},    {nDim, 0, 0},    {0, nDim, 0},    {nDim, nDim, 0},
                                 {0, 0, nDim}, {nDim, 0, nDim}, {0, nDim, nDim}, {nDim, nDim, nDim}};
  for(MInt corner = 0; corner < IPOW2(nDim); corner++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      pp[dim] = p[dim + cornerCode[corner][dim]];
    }
    if(!(maia::geom::isPointInsideSphere<nDim>(pp, c, R))) {
      return false;
    }
  }
  return true;
}

template <typename T, MInt nDim>
MBool isSphereInsideBox(const T* const c, const T R, const T* const p) {
  T bbox[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    bbox[dim] = c[dim] - R;
    bbox[dim + nDim] = c[dim] + R;
  }
  return maia::geom::isBoxInsideBox<T, nDim>(bbox, p);
}


template <MInt nDim>
MBool doBoxAndSphereOverlap(const MFloat* const b, const MFloat* const c, const MFloat R) {
  MFloat cp[nDim];
  const MInt cornerCode[8][3] = {{0, 0, 0},    {nDim, 0, 0},    {0, nDim, 0},    {nDim, nDim, 0},
                                 {0, 0, nDim}, {nDim, 0, nDim}, {0, nDim, nDim}, {nDim, nDim, nDim}};

  // easy part
  // is sphere center inside box?
  if(maia::geom::isPointInsideBox<MFloat, nDim>(c, b)) {
    return true;
  }

  // is a box corner in sphere?
  for(MInt corner = 0; corner < IPOW2(nDim); corner++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      cp[dim] = b[dim + cornerCode[corner][dim]];
    }
    if(maia::geom::isPointInsideSphere<nDim>(cp, c, R)) {
      return true;
    }
  }

  // 'tricky part' reducing the dimensions
  // what region relative to the cube is the sphere center located in
  MInt insideCnt = 0;
  MFloat p[nDim];
  MFloat lc[nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    if(c[dim] > b[dim] && c[dim] < b[dim + nDim]) {
      p[dim] = F0;
      lc[dim] = F0;
      insideCnt++;
    } else {
      lc[dim] = c[dim];
      p[dim] = (c[dim] <= b[dim]) ? b[dim] : b[dim + nDim];
    }
  }

  // the closest point is a corner, we already checked them
  if(insideCnt == 0) {
    return false;
  }

  // the cosest point is found checking the outside dimensions
  return maia::geom::isPointInsideSphere<nDim>(p, lc, R);
}

template <MInt nDim>
MBool doesLinePenetrateBox(const MFloat* const line, const MFloat* const bbox) {
  IF_CONSTEXPR(nDim > 3 || nDim < 2) { mTerm(1, AT_, "neither 2 nor 3 dimensional"); }

  // create a line bounding box
  MFloat bbox_line[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    if(line[dim] < line[dim + nDim]) {
      bbox_line[dim] = line[dim];
      bbox_line[dim + nDim] = line[dim + nDim];
    } else {
      bbox_line[dim] = line[dim + nDim];
      bbox_line[dim + nDim] = line[dim];
    }
  }

  // check bounding boxes for overlap
  // a separation axis can be found in the cartesian axis
  if(!maia::geom::doBoxesOverlap<MFloat, nDim>(bbox_line, bbox)) {
    return false;
  }

  // translate box to center and scale to unit cube
  MFloat line_c[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    const MFloat halfLength = (bbox[dim + nDim] - bbox[dim]) / F2;
    line_c[dim] = (line[dim] - bbox[dim]) / halfLength - 0.5;
    line_c[dim + nDim] = (line[dim + nDim] - bbox[dim]) / halfLength - 0.5;
  }

  MFloat line_span[3] = {F0, F0, F0};
  MFloat lineLength = F0;
  for(MInt dim = 0; dim < nDim; dim++) {
    line_span[dim] = line_c[dim + nDim] - line_c[dim];
    lineLength += POW2(line_span[dim]);
  }
  lineLength = std::sqrt(lineLength);
  for(MInt dim = 0; dim < nDim; dim++) {
    line_span[dim] /= lineLength;
  }

  const MFloat radius = nDim * F1B4;
  for(MInt e = 2; e > ((nDim == 2) ? 1 : -1); e--) {
    // generate separation axis
    MFloat axis[3];
    MFloat e_vec[3] = {F0, F0, F0};
    e_vec[e] = F1;
    maia::math::cross(line_span, e_vec, axis);

    // project one of the points on the axis
    MFloat line_projection = F0;
    for(MInt dim = 0; dim < nDim; dim++) {
      line_projection += axis[dim] * line_c[dim];
    }

    // larger equal radius means no penetration but maybe a contactpoint
    if(POW2(line_projection) >= radius) {
      return false;
    }

    // project the corner of the cube ...
    // which corner to project
    MFloat corner_projection = F0;
    for(MInt dim = 0; dim < nDim; dim++) {
      corner_projection += axis[dim] * ((axis[dim] > 0) ? 0.5 : -0.5);
    }

    // larger equal radius means no penetration but maybe a contactpoint
    if(std::fabs(line_projection) >= std::fabs(corner_projection)) {
      return false;
    }
  }

  // if no separation axis can be found, there has to be a penetration
  return true;
}

template <MInt nDim>
MBool doesLinePenetrateSphere(const MFloat* const line, const MFloat* const center, const MFloat radius) {
  const MBool p1_inside = maia::geom::isPointInsideSphere<nDim>(line, center, radius);
  const MBool p2_inside = maia::geom::isPointInsideSphere<nDim>(line + nDim, center, radius);

  if(p1_inside && p2_inside) {
    return false;
  }

  if(p1_inside != p2_inside) {
    return true;
  }

  // at this point both end points of the line are outside the sphere

  // create a line and sphere bounding box
  MFloat bbox_line[2 * nDim];
  MFloat bbox_sphere[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    bbox_sphere[dim] = center[dim] - radius;
    bbox_sphere[dim + nDim] = center[dim] + radius;
    if(line[dim] < line[dim + nDim]) {
      bbox_line[dim] = line[dim];
      bbox_line[dim + nDim] = line[dim + nDim];
    } else {
      bbox_line[dim] = line[dim + nDim];
      bbox_line[dim + nDim] = line[dim];
    }
  }

  // check bounding boxes for overlap
  // a separation axis can be found in the cartesian axis
  if(!maia::geom::doBoxesOverlap<MFloat, nDim>(bbox_line, bbox_sphere)) {
    return false;
  }

  MFloat a = F0;
  MFloat b = F0;
  MFloat c = -POW2(radius);
  for(MInt dim = 0; dim < nDim; dim++) {
    const MFloat dl = line[dim + nDim] - line[dim];
    const MFloat dlc = line[dim] - center[dim];
    a += POW2(dl);
    b += dl * dlc;
    c += POW2(dlc);
  }

  const MFloat D = POW2(b / a) - c / a;
  if(D <= F0) {
    return false;
  }

  // at this point there have to be two valid intersections or none
  // it should be enough to check one of them

  const MFloat t_p = b / a + sqrt(D);

  // check the intersection if they lie on the line segment
  return (t_p <= F1) && (t_p >= F0);
}

} // namespace geom
} // namespace maia
#endif // ANALYTICGEOMETRY_H
