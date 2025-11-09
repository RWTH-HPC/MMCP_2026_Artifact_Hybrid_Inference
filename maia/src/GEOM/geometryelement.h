// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYELEMENT_H
#define GEOMETRYELEMENT_H

/** Representation of the geometrical concept of an element
 *
 * Collector class!
 * In 2D an element consists of two vertices and is part of a segment.
 * In 3D an element consists of at least 3 vertices and is part of a segment.
 */
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class GeometryElement {
 public:
  MInt m_segmentId;
  MInt m_bndCndId;
  MInt m_originalId;
  MFloat** m_vertices;
  MFloat* m_normal;
  // Lower nDim values hold minimum
  // Upper nDim values hold maximu
  MFloat* m_minMax;

 public:
  void boundingBox();
  void writeElement() const;
  void calcNormal(const MFloat* const vertices, MFloat* normal) const;
  void calcCentroid(const MFloat* const vertices, MFloat* centroid) const;
  void getVertices(MFloat* vertices) const;
  void allocateElements(void*, void*, MInt&);
  static void init(MInt NotUsed(dimensions), MInt /*noDistributions*/, MInt /*maxSize*/){
      //    TRACE();
  };

  static MInt staticElementSize() { return (nDim + 3) * nDim * sizeof(MFloat) + sizeof(MFloat*) * nDim; }
};

#endif
