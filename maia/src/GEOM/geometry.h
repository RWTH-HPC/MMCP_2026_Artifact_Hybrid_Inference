// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <array>
#include <set>
#include <vector>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "geometrycontext.h"
#include "geometrycontexttypes.h"

template <MInt nDim>
class GeometryAdt;
template <MInt nDim>
class GeometryElement;
template <class T>
class Collector;

// Auxiliary type traits to support selecting the correct geometry class, instantiated in
// geometry2/3d.h
template <MInt nDim>
struct GeometryXD {};

template <MInt nDim>
using element = GeometryElement<nDim>;
/** Basic implementation
 *
 */
template <MInt nDim>
class Geometry {
 public:
  template <MInt nDim_>
  friend class GeometryAdt;
  Geometry(const MInt solverId_, const MPI_Comm comm);

  virtual ~Geometry() = default;

  MInt solverId() const { return m_solverId; }
  MPI_Comm mpiComm() const { return m_mpiComm; }
  MInt domainId() const { return m_domainId; }
  MInt noDomains() const { return m_noDomains; }
  GeometryContext& geometryContext() { return m_geometryContext; }
  void getBoundingBox(MFloat* const bBox) const;
  const MFloat* boundingBox() const { return &m_minMax[0]; }
  void setHaloElementOffset(MInt off) { m_haloElementOffset = off; }
  MInt getHaloElementOffset() const { return m_haloElementOffset; }
  MBool isOnGeometry(const MFloat, const MFloat*, MString);

  // for all dimensions
  MBool vectorsEqual(MFloat* a, MFloat* b);
  MFloat calcCircumference(MFloat** bndVs, MInt num);

  virtual MInt boundaryCheck(MFloat* /*targetRegion*/, MFloat /*cellHalfLength*/, MFloat* /*cell_coords*/,
                             MInt* /*bndIds*/) {
    return false;
  };
  virtual MInt getIntersectionElements(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/) { return 0; };
  virtual MInt getIntersectionElements(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/,
                                       MFloat /*cellHalfLength*/, const MFloat* const /*cell_coords*/) {
    return 0;
  };
  virtual MInt getLineIntersectionElementsOld1(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/) { return 0; };
  virtual MInt getLineIntersectionElementsOld2(MFloat* /*targetRegion*/, MInt* /*spaceDirection*/,
                                               std::vector<MInt>& /*nodeList*/) {
    return 0;
  };
  virtual MInt getLineIntersectionElements(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/) { return 0; };
  virtual MInt getLineIntersectionElements(MFloat* /*targetRegion*/) { return 0; };
  virtual MBool getClosestLineIntersectionLength(MInt /*bndCndId*/, const std::vector<MInt>& /*nodeList*/,
                                                 MFloat* /*targetRegion*/, MFloat* /*dist*/) {
    return false;
  };

  void getLineIntersectingElementsBcIds(const MFloat* const line, std::set<MInt>& bcIds);

  virtual MBool edgeTriangleIntersection(MFloat* /*trianglePoint1*/,
                                         MFloat* /*trianglePoint2*/,
                                         MFloat* /*trianglePoint3*/,
                                         MFloat* /*edgePoint1*/,
                                         MFloat* /*edgePoint2*/) {
    return 0;
  };
  virtual MBool edgeTriangleIntersectionLB(MFloat* /*trianglePoint1*/,
                                           MFloat* /*trianglePoint2*/,
                                           MFloat* /*trianglePoint3*/,
                                           MFloat* /*edgePoint1*/,
                                           MFloat* /*edgePoint2*/) {
    return 0;
  };

  virtual MBool getLineTriangleIntersectionSimple(MFloat* /*p1*/, MFloat* /*p2*/, MFloat* /*v1*/, MFloat* /*v2*/,
                                                  MFloat* /*v3*/) {
    return false;
  };
  virtual MBool getLineTriangleIntersectionSimpleDistance(MFloat* /*p1*/, MFloat* /*p2*/, MFloat* /*v1*/,
                                                          MFloat* /*v2*/, MFloat* /*v3*/, MFloat* /*dist*/) {
    return false;
  };
  virtual MBool getLineTriangleIntersection(const MFloat* const /*p1*/,
                                            const MFloat* const /*p2*/,
                                            const MFloat /*radius*/,
                                            const MFloat* const /*v1*/,
                                            const MFloat* const /*v2*/,
                                            const MFloat* const /*v3*/,
                                            MFloat* /*intersection*/,
                                            MFloat* /*normal*/,
                                            MFloat* /*lambda2*/,
                                            MFloat* /*dist*/) {
    return false;
  };
  // moving boundary
  void getBoundingBoxMB(MFloat* const bBox) const;
  virtual MInt getIntersectionMBElements(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/) { return 0; }
  virtual MInt getLineIntersectionMBElements(MFloat* /*targetRegion*/, std::vector<MInt>& /*nodeList*/) { return 0; }
  virtual MInt getLineIntersectionMBElements2(MFloat* /*targetRegion*/, MInt* /*spaceDirection*/,
                                              std::vector<MInt>& /*nodeList*/, MInt /*bcId*/) {
    return 0;
  }
  virtual MInt getSphereIntersectionMBElements(MFloat* /* P*/, MFloat /*radius*/, std::vector<MInt>& /*nodeList*/) {
    return 0;
  };
  virtual void MoveAllMBElementVertex(MFloat* /*dx*/){};
  virtual void MoveMBElementVertex(MInt /*e*/, MInt /*v*/, MFloat* /*dx*/){};
  virtual void ReplaceMBElementVertex(MInt /*e*/, MInt /*v*/, MFloat* /*np*/){};
  virtual void UpdateMBNormalVector(MInt /*e*/){};
  virtual void UpdateMBBoundingBox(){};
  virtual void UpdateADT(){};
  virtual void collectGlobalMemoryUsage(){};


  virtual void writeSTL(const MChar* /*fileName*/){};
  virtual void writeADTAndSTLToNetCDF(const MChar* /*fileName*/){};
  virtual void writeSTLMB(const MChar* /*fileName*/, MInt& /*noNodes*/, MInt*& /*nodeList*/){};
  virtual void readSTLNetCDF(const MChar* /*fileName*/){};

  virtual void logStatistics(){};

  virtual MInt GetNoElements() { return m_noElements; };
  virtual MInt GetNoSegments() { return m_noSegments; };
  virtual MInt* GetBoundaryIds(MInt* noAllBcs) {
    *noAllBcs = m_noAllBCs;
    return m_allBCs;
  };

  virtual MFloat** GetBoundaryVertices(MInt /*segmentId*/, MFloat* /*tri_vx*/, MInt* /*keepOffsets*/, MInt /*size*/,
                                       MInt* /*num*/) {
    return nullptr;
  };
  virtual inline MBool isEdgeAlreadyInCollection(std::vector<std::pair<MFloat*, MFloat*>> /*tmp_edges*/, MFloat* /*p1*/,
                                                 MFloat* /*p2*/, MInt* /*pos*/) {
    return false;
  };
  virtual MFloat GetBoundarySize(MInt /*segmentId*/) { return F0; };
  virtual MFloat GetBoundarySize(MFloat* /*vertices*/, MInt* /*keepOffset*/, MInt /*size*/) { return F0; };
  virtual void determineSegmentOwnership(MInt /*segmentId*/, MInt* /*own*/, MInt* /*sumowners*/, MInt* /*firstOwner*/,
                                         MInt* /*owners*/){};
  virtual MFloat getBndMaxRadius(MFloat** /*vertices*/, MInt /*num*/) { return F0; };

  virtual void rebuildAdtTree(){};
  virtual void calculateBoundingBox(){};
  virtual void writeParallelGeometryVTK(MString /*filename*/){};
  virtual void addElement(MFloat* /*tri*/){};
  virtual void copyElement(MInt /*from*/, MInt /*to*/){};
  virtual void resizeCollector(MInt /*new_size*/){};

  // Methods migrated from CartesianGrid
  MBool pointIsInside(const MFloat* const coordinates);
  MBool pointIsInside(const MFloat* const coordinates, MInt* numcutsperdir);
  MBool pointIsInside2(const MFloat* const coordinates, MInt* numcutsperdir = nullptr);
  MBool pointIsInsideMBElements(const MFloat* const coordinates, MInt*, MInt*, MInt);
  MBool pointIsInsideMBElements2(const MFloat* const coordinates, MInt*, MInt*, MInt);
  void determineRayIntersectedElements(const MFloat* const coordinates, std::vector<std::vector<MInt>>* resultnodes);

 protected:
  virtual void readSegments(){};
  virtual void writeSegmentsToDX();

 protected:
  MInt m_noSegments;

  std::array<MFloat, 2 * nDim> m_minMax{};

  MInt m_noElements;

  MString m_segmentBaseName;
  bodyMap m_bodyMap;
  bodyIterator m_bodyIt;

  MInt m_noMBElements;
  MInt m_noBoundaryIds;
  MInt* m_boundaryIds;
  MInt* m_allBCs = nullptr;
  MInt m_noAllBCs;

  MBool m_flowSolver;

 private:
  const MInt m_solverId;
  const MPI_Comm m_mpiComm;
  MInt m_domainId;
  MInt m_noDomains;
  GeometryContext m_geometryContext;
  MInt m_haloElementOffset;

 public:
  std::vector<MInt> m_segmentOffsets;
  std::vector<MInt> m_segmentOffsetsWithoutMB;
  Collector<element<nDim>>* m_elements;
  element<nDim>* elements;
  GeometryAdt<nDim>* m_adt;
  Collector<element<nDim>>* m_mbelements;
  element<nDim>* mbelements;
  std::array<MFloat, 2 * nDim> m_mbminMax{};
  MFloat m_mbMidPnt[3];
  MInt noBoundaryIds() { return m_noBoundaryIds; }
  MBool* m_ownSegmentId = nullptr;
  std::set<MInt> m_uniqueOriginalTriId;
  MFloat m_parGeomMemFactor;
  MString m_inOutTest;

  // Parallel geometry
  MBool m_parallelGeometry = false;
  MBool m_debugParGeom = false;
  MString m_parallelGeomFileName;
};


#endif
