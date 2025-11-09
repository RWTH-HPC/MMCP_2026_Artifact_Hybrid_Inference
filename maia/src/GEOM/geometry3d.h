// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRY3D_H
#define GEOMETRY3D_H

#include "geometry.h"

/** Specialized 3D implementation
 *
 */
class Geometry3D : public Geometry<3> {
 public:
  Geometry3D(const MInt solverId_, const MPI_Comm comm);
  Geometry3D(const MInt solverId_, const MString& filename, const MPI_Comm comm);
  ~Geometry3D();
  MInt getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList, MFloat cellHalfLength,
                               const MFloat* const cell_coords) override;
  MInt getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) override;
  virtual MInt getIntersectionElementsTetraeder(MFloat* targetRegion, std::vector<MInt>& nodeList);
  MInt getLineIntersectionElementsOld2(MFloat* targetRegion, MInt* spaceDirection,
                                       std::vector<MInt>& nodeList) override;
  MInt getLineIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) override;
  MInt getLineIntersectionElements(MFloat* targetRegion) override;
  MBool getClosestLineIntersectionLength(MInt bndCndId, const std::vector<MInt>& nodeList, MFloat* targetRegion,
                                         MFloat* dist) override;

  MBool edgeTriangleIntersection(MFloat* trianglePoint1,
                                 MFloat* trianglePoint2,
                                 MFloat* trianglePoint3,
                                 MFloat* edgePoint1,
                                 MFloat* edgePoint2) override;
  MBool edgeTriangleIntersectionLB(MFloat* trianglePoint1,
                                   MFloat* trianglePoint2,
                                   MFloat* trianglePoint3,
                                   MFloat* edgePoint1,
                                   MFloat* edgePoint2) override;

  MBool getLineTriangleIntersectionSimple(MFloat* p1, MFloat* p2, MFloat* v1, MFloat* v2, MFloat* v3) override;
  MBool getLineTriangleIntersectionSimpleDistance(MFloat* p1, MFloat* p2, MFloat* v1, MFloat* v2, MFloat* v3,
                                                  MFloat* dist) override;
  MBool getLineTriangleIntersection(const MFloat* const p1, const MFloat* const p2, const MFloat radius,
                                    const MFloat* const v1, const MFloat* const v2, const MFloat* const v3,
                                    MFloat* intersection, MFloat* normal, MFloat* lambda2, MFloat* dist) override;

  // moving boundary
  MInt getIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) override;
  MInt getLineIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) override;
  MInt getLineIntersectionMBElements2(MFloat* targetRegion, MInt* spaceDirection, std::vector<MInt>& nodeList,
                                      MInt bcIc) override;
  MInt getSphereIntersectionMBElements(MFloat* P, MFloat radius, std::vector<MInt>& nodeList) override;
  void MoveAllMBElementVertex(MFloat* dx) override;
  void MoveMBElementVertex(MInt e, MInt v, MFloat* dx) override;
  void ReplaceMBElementVertex(MInt e, MInt v, MFloat* np) override;
  void UpdateMBNormalVector(MInt e) override;
  void UpdateMBBoundingBox() override;
  void UpdateADT() override;

  void writeSTL(const char* fileName) override;
  void writeSTLMB(const char* fileName, MInt& noNodes, MInt*& nodeList) override;
  void writeADTAndSTLToNetCDF(const char* fileName) override;
  void readSTLNetCDF(const char* fileName) override;

  void logStatistics() override;

  MFloat** GetBoundaryVertices(MInt segmentId, MFloat* tri_vx, MInt* keepOffsets, MInt size, MInt* num) override;
  virtual inline std::vector<std::pair<MFloat*, MFloat*>> GetUniqueSegmentEdgesParGeom(MFloat* tri_vx,
                                                                                       MInt* keepOffsets, MInt size);
  virtual inline std::vector<std::pair<MFloat*, MFloat*>> GetUniqueSegmentEdges(MInt segmentId);
  inline MBool isEdgeAlreadyInCollection(std::vector<std::pair<MFloat*, MFloat*>> tmp_edges, MFloat* p1, MFloat* p2,
                                         MInt* pos) override;
  MFloat GetBoundarySize(MInt segmentId) override;
  MFloat GetBoundarySize(MFloat* vertices, MInt* keepOffset, MInt size) override;
  void determineSegmentOwnership(MInt segmentId, MInt* own, MInt* sumowners, MInt* firstOwner, MInt* owners) override;
  MFloat getBndMaxRadius(MFloat** vertices, MInt num) override;

 protected:
  void readSegments() override;
  virtual void correctVertexCoordinates();
  virtual inline void swap4BytesToBE(char* buf);
  virtual MInt is_big_endian();

  virtual inline Collector<element<3>>* readSegmentsSerial();
  virtual inline Collector<element<3>>* readSegmentsParallel();
  virtual inline void countSegmentTrianglesASCII(MString fileName, MInt* noElements);
  virtual inline void countSegmentTrianglesBINARY(MString fileName, MInt* noElements);
  virtual inline void countSegmentTrianglesNETCDF(MString fileName, MInt* noElements, const MPI_Comm comm);
  virtual inline void readSegmentTrianglesASCII(MString fileName, Collector<element<3>>* elemCollector, MInt bndCndId,
                                                MInt segmentId, MInt* offset);
  virtual inline void readSegmentTrianglesBINARY_BE(MString fileName, Collector<element<3>>* elemCollector,
                                                    MInt bndCndId, MInt segmentId, MInt* offset);
  virtual inline void readSegmentTrianglesBINARY_LE(MString fileName, Collector<element<3>>* elemCollector,
                                                    MInt bndCndId, MInt segmentId, MInt* offset);
  virtual inline void readSegmentTrianglesNETCDF(MString fileName, Collector<element<3>>* elemCollector, MInt bndCndId,
                                                 MInt segmentId, MInt* offset, const MPI_Comm comm);

  void rebuildAdtTree() override;
  void calculateBoundingBox() override;
  void writeParallelGeometryVTK(MString filename) override;

  void printMemoryUsage();
  void collectGlobalMemoryUsage() override;
  void addElement(MFloat* tri) override;
  void copyElement(MInt from, MInt to) override;
  void copyElement(MInt from, MInt to, element<3>* fromPtr, element<3>* toPtr);
  void resizeCollector(MInt new_size) override;

  MInt otherCalls = 0;

  MBool m_GFieldInitFromSTL = false;
  MInt m_levelSetIntfBndId = 0;
  MInt* m_levelSetIntfBndIds{};
  MInt m_noLevelSetIntfBndIds = 0;
  MBool m_forceBoundingBox = false;

  MString m_gridCutTest;

  MString m_geomFileType;
  MString m_gridFileName;
  MInt m_noAllTriangles{};

  MInt m_getLIE2CallCounter = 0;
  MInt getLIE2CommCounter = 0;
  MInt getIECallCounter = 0;
  MInt getIECommCounter = 0;
  MInt getIETCallCounter = 0;
  MInt edgeTICallCounter = 0;
  static constexpr const MInt nDim = 3;

  MBool m_communicateSegmentsSerial = true;

  // TIMERS
 private:
  MInt m_tg_geometry{};
  MInt m_t_geometryAll{};
  MInt m_t_readGeometry{};
};

template <>
struct GeometryXD<3> {
  using type = Geometry3D;
};

#endif
