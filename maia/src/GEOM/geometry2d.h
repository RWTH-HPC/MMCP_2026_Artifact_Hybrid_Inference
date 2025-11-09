// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRY2D_H
#define GEOMETRY2D_H

#include "geometry.h"

/** Specialized 2D implementation
 *
 */
class Geometry2D : public Geometry<2> {
 public:
  Geometry2D(const MInt solverId_, const MPI_Comm comm);
  Geometry2D(const MInt solverId_, const MString filename, const MPI_Comm comm);
  ~Geometry2D();

  virtual MInt getLineIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList);
  virtual MInt getLineIntersectionElementsOld1(MFloat* targetRegion, std::vector<MInt>& nodeList);
  virtual MInt getLineIntersectionElementsOld2(MFloat* targetRegion, MInt* spaceDirection, std::vector<MInt>& nodeList);
  virtual MInt getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList);
  virtual MInt getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList, MFloat cellHalfLength,
                                       const MFloat* const cell_coords);
  virtual MBool edgeTriangleIntersection(MFloat* trianglePoint1, MFloat* trianglePoint2, MFloat* trianglePoint3,
                                         MFloat* edgePoint1, MFloat* edgePoint2);

  virtual MInt getIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList);
  virtual MInt getLineIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList);
  virtual MInt getSphereIntersectionMBElements(MFloat* P, MFloat radius, std::vector<MInt>& nodeList);
  virtual void MoveAllMBElementVertex(MFloat* dx);
  virtual void MoveMBElementVertex(MInt e, MInt v, MFloat* dx);
  virtual void ReplaceMBElementVertex(MInt e, MInt v, MFloat* np);
  virtual void UpdateMBNormalVector(MInt /*e*/){};
  virtual void UpdateMBBoundingBox();
  virtual void UpdateADT();

 protected:
  virtual void readSegments();
  inline void countSegmentLinesASCII(const MString& fileName, MInt* noElements);
  inline void readSegmentLinesASCII(MString fileName, Collector<element<2>>* elemCollector, MInt bndCndId,
                                    MInt segmentId, MInt* offset);
  void calculateBoundingBox();

  MBool m_GFieldInitFromSTL;
  MInt m_levelSetIntfBndId{};
  MInt* m_levelSetIntfBndIds{};
  MInt m_noLevelSetIntfBndIds{};

  static constexpr const MInt nDim = 2;
};

template <>
struct GeometryXD<2> {
  using type = Geometry2D;
};

#endif
