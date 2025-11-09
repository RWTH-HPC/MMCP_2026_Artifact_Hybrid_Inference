// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYCONCEPT_H
#define GEOMETRYCONCEPT_H

#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"

template <MInt nDim>
class Geometry;

class SolverSurface {
 public:
  SolverSurface(){};

  virtual MBool getCellIntersectingElements(const MFloat* const /*coords*/, const MFloat /*cellHalfLength*/,
                                            MBool* const /* cutInfo */) {
    mTerm(1, AT_, "only virt");
    return 0;
  };
  virtual MInt countLineIntersectingElements(const MFloat* const /* line */) {
    mTerm(1, AT_, "only virt");
    return 0;
  };
  virtual void boundingBox(MFloat* const /* bbox */) { mTerm(1, AT_, "only virt"); };
  virtual MInt noSegments() {
    mTerm(1, AT_, "only virt");
    return 0;
  };

  // temporary to make it work with the demand of stl elements
  virtual void setGeometryPointer(Geometry<2>*& /* geometryPointer */) { mTerm(1, AT_, "only virt"); };
  virtual void setGeometryPointer(Geometry<3>*& /* geometryPointer */) { mTerm(1, AT_, "only virt"); };

 private:
};

template <MInt nDim>
class SolverSTLSurface : public SolverSurface {
 public:
  SolverSTLSurface(const MPI_Comm comm, const MInt solver);
  MBool getCellIntersectingElements(const MFloat* const coords, const MFloat cellHalfLength, MBool* const cutInfo);
  MInt countLineIntersectingElements(const MFloat* const line);
  void boundingBox(MFloat* const bbox);
  MInt noSegments();
  void setGeometryPointer(Geometry<nDim>*& geometryPointer) { geometryPointer = m_stlGeometry; }; // temporary

 private:
  Geometry<nDim>* m_stlGeometry;
};

template <MInt nDim>
class SolverAnalyticBoxSurface : public SolverSurface {
 public:
  SolverAnalyticBoxSurface(const MFloat* const bbox);
  MBool getCellIntersectingElements(const MFloat* const cell_coords, const MFloat cellHalfLength, MBool* const cutInfo);
  MInt countLineIntersectingElements(const MFloat* const line);
  void boundingBox(MFloat* const bbox);
  MInt noSegments();

 private:
  // just to have something here
  MFloat m_bbox[2 * nDim];
};

template <MInt nDim>
class SolverAnalyticSphereSurface : public SolverSurface {
 public:
  SolverAnalyticSphereSurface(const MFloat* const c, const MFloat R);
  MBool getCellIntersectingElements(const MFloat* const cell_coords, const MFloat cellHalfLength, MBool* const cutInfo);
  MInt countLineIntersectingElements(const MFloat* const line);
  void boundingBox(MFloat* const bbox);
  MInt noSegments();

 private:
  // just to have something here
  MFloat m_center[nDim];
  MFloat m_radius;
};

// hierarchic geometry concept

struct GeometryDistribution {
  MInt noNodes;
  MInt* noSurfacesPerNode = nullptr;
  MInt** nodeSurfaceIds = nullptr;
};

class GeometryBase {
 public:
  GeometryBase(){};
  GeometryBase(const MPI_Comm comm) : m_mpiComm(comm){};
  virtual MBool isPointInside(const MFloat* const /* point */) {
    mTerm(1, AT_, "only virt");
    return false;
  };

  virtual MBool isPointInsideNode(const MFloat* const /* point */, const MInt /* node */) {
    mTerm(1, AT_, "only virt");
    return false;
  };
  virtual MBool getCellIntersectingSurfaces(const MFloat* const /* coords */, const MFloat /* cellHalfLength */,
                                            MBool* const* const /* cutInfo */) {
    mTerm(1, AT_, "only virt");
    return 0;
  };
  virtual MBool getCellIntersectingSurfacesOfNode(const MFloat* const /* coords */, const MFloat /* cellHalfLength */,
                                                  MBool* const /* cutInfo */, const MInt /* node */) {
    mTerm(1, AT_, "only virt");
    return 0;
  };


  GeometryDistribution m_distribution;

  // The MPI communicator to be used by a node (or the root)
  MPI_Comm m_mpiComm;
};

class GeometryRoot : public GeometryBase {
 public:
  GeometryRoot(){};
  GeometryRoot(const MInt noSolvers, const MInt nDim_, const MPI_Comm comm);

  template <MInt nDim>
  void initGeometry(const MInt noSolvers);

  MBool isPointInside(const MFloat* const point);
  MBool isPointInsideNode(const MFloat* const point, const MInt node);
  MBool getCellIntersectingSurfaces(const MFloat* const coords, const MFloat cellHalfLength,
                                    MBool* const* const cutInfo);
  MBool getCellIntersectingSurfacesOfNode(const MFloat* const coords, const MFloat cellHalfLength, MBool* const cutInfo,
                                          const MInt node);

  void boundingBox(MFloat* const bBox);
  void boundingBoxOfNode(MFloat* const bBox, const MInt node);
  MInt noSegmentsOfNode(const MInt node);
  MInt noNodes() { return m_distribution.noNodes; };
  MInt nodeSurfaceType(const MInt node) { return m_solverSurfaceType[node]; };

  // temporary
  template <MInt nDim>
  void setGeometryPointerToNode(Geometry<nDim>*& geometryPointer, const MInt node);

 private:
  MInt m_nDim;

  MInt m_noSolverSurfaces;
  MInt* m_solverSurfaceType = nullptr;
  SolverSurface** m_solverSurface;
};

class GeometryNode : public GeometryBase {
 public:
  GeometryNode(GeometryBase* geometryRoot, MInt node, MPI_Comm comm);

 private:
  // pointer to lower node (my trees dont grow in australia)
  GeometryBase* m_geometry;

  // node number (identifier) at lower node
  MInt m_node;
};

#endif // GEOMETRYCONCEPT_H
