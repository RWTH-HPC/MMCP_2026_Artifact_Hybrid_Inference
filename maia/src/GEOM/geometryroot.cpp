// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryroot.h"
#include <climits>
#include <iostream>
#include "INCLUDE/maiaconstants.h"
#include "IO/context.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "enums.h"
#include "geometry2d.h"
#include "geometry3d.h"
#include "geometryadt.h"
#include "geometryanalytic.h"
#include "property.h"

using namespace std;


// implementations of various SolverXXXSurface

template <MInt nDim>
SolverSTLSurface<nDim>::SolverSTLSurface(const MPI_Comm comm, const MInt solver) {
  m_stlGeometry = new typename GeometryXD<nDim>::type(solver, comm);
}

template <MInt nDim>
MInt SolverSTLSurface<nDim>::noSegments() {
  return m_stlGeometry->geometryContext().getNoSegments();
}

template <MInt nDim>
void SolverSTLSurface<nDim>::boundingBox(MFloat* const bbox) {
  m_stlGeometry->getBoundingBox(bbox);
}

template <MInt nDim>
MBool SolverSTLSurface<nDim>::getCellIntersectingElements(const MFloat* const coords, const MFloat cellHalfLength,
                                                          MBool* const cutInfo) {
  std::vector<MInt> nodeList;
  MFloat targetRegion[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    targetRegion[dim] = coords[dim] - cellHalfLength;
    targetRegion[dim + nDim] = coords[dim] + cellHalfLength;
  }

  // reset cutInfo
  for(MInt seg = 0; seg < noSegments(); seg++) {
    cutInfo[seg] = false;
  }

  // call geometry
  m_stlGeometry->getIntersectionElements(targetRegion, nodeList, cellHalfLength, coords);
  const MInt noNodes = nodeList.size();

  // translate from nodes to segments/elements/partialsurfaces (no unique nomenclature)
  for(MInt n = 0; n < noNodes; n++) {
    cutInfo[m_stlGeometry->elements[nodeList[n]].m_segmentId] = true;
  }
  return (MBool)noNodes;
}

template <MInt nDim>
MInt SolverSTLSurface<nDim>::countLineIntersectingElements(const MFloat* const line) {
  std::vector<MInt> nodeList;
  MFloat tmp_line[2 * nDim];

  for(MInt dim = 0; dim < nDim; dim++) {
    tmp_line[dim] = line[dim];
    tmp_line[dim + nDim] = line[dim + nDim];
  }

  // call geometry
  m_stlGeometry->getLineIntersectionElements(tmp_line, nodeList);

  return nodeList.size();
}

// analytic box

template <MInt nDim>
SolverAnalyticBoxSurface<nDim>::SolverAnalyticBoxSurface(const MFloat* const bbox) {
  for(MInt i = 0; i < nDim; i++) {
    if(bbox[i + nDim] < bbox[i]) mTerm(1, AT_, "invalid bounding box");
    m_bbox[i] = bbox[i];
    m_bbox[i + nDim] = bbox[i + nDim];
  }
}

template <MInt nDim>
MInt SolverAnalyticBoxSurface<nDim>::noSegments() {
  return 1;
}

template <MInt nDim>
void SolverAnalyticBoxSurface<nDim>::boundingBox(MFloat* const bbox) {
  for(MInt dim = 0; dim < nDim; dim++) {
    bbox[dim] = m_bbox[dim];
    bbox[dim + nDim] = m_bbox[dim + nDim];
  }
}

template <MInt nDim>
MBool SolverAnalyticBoxSurface<nDim>::getCellIntersectingElements(const MFloat* const cell_coords,
                                                                  const MFloat cellHalfLength,
                                                                  MBool* const cutInfo) {
  MFloat cell_bbox[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    cell_bbox[dim] = cell_coords[dim] - cellHalfLength;
    cell_bbox[dim + nDim] = cell_coords[dim] + cellHalfLength;
  }

  if(maia::geom::doBoxesOverlap<MFloat, nDim>(cell_bbox, m_bbox)
     && !maia::geom::isBoxInsideBox<MFloat, nDim>(cell_bbox, m_bbox)) {
    cutInfo[0] = true;
  } else {
    cutInfo[0] = false;
  }
  return cutInfo[0];
}

template <MInt nDim>
MInt SolverAnalyticBoxSurface<nDim>::countLineIntersectingElements(const MFloat* const line) {
  // check both ends of the line
  const MBool p1_inside = maia::geom::isPointInsideBox<MFloat, nDim>(line, m_bbox);
  const MBool p2_inside = maia::geom::isPointInsideBox<MFloat, nDim>(line + nDim, m_bbox);

  if(p1_inside && p2_inside) return 0;

  if(p1_inside != p2_inside) return 1;

  // in case both end are outside we need to take a look
  if(maia::geom::doesLinePenetrateBox<nDim>(line, m_bbox))
    return 2;
  else
    return 0;
}

// analytic sphere

template <MInt nDim>
SolverAnalyticSphereSurface<nDim>::SolverAnalyticSphereSurface(const MFloat* const c, const MFloat R) {
  for(MInt i = 0; i < nDim; i++) {
    m_center[i] = c[i];
  }
  m_radius = R;
}

template <MInt nDim>
MInt SolverAnalyticSphereSurface<nDim>::noSegments() {
  return 1;
}

template <MInt nDim>
void SolverAnalyticSphereSurface<nDim>::boundingBox(MFloat* const bbox) {
  for(MInt dim = 0; dim < nDim; dim++) {
    bbox[dim] = m_center[dim] - m_radius;
    bbox[dim + nDim] = m_center[dim] + m_radius;
  }
}

template <MInt nDim>
MBool SolverAnalyticSphereSurface<nDim>::getCellIntersectingElements(const MFloat* const cell_coords,
                                                                     const MFloat cellHalfLength,
                                                                     MBool* const cutInfo) {
  MFloat cell_bbox[2 * nDim];
  for(MInt dim = 0; dim < nDim; dim++) {
    cell_bbox[dim] = cell_coords[dim] - cellHalfLength;
    cell_bbox[dim + nDim] = cell_coords[dim] + cellHalfLength;
  }

  if(maia::geom::doBoxAndSphereOverlap<nDim>(cell_bbox, m_center, m_radius)
     && !maia::geom::isBoxInsideSphere<nDim>(cell_bbox, m_center, m_radius)) {
    cutInfo[0] = true;
  } else {
    cutInfo[0] = false;
  }
  return cutInfo[0];
}

template <MInt nDim>
MInt SolverAnalyticSphereSurface<nDim>::countLineIntersectingElements(const MFloat* const line) {
  if(!maia::geom::doesLinePenetrateSphere<nDim>(line, m_center, m_radius)) return 0;

  // check start of line
  const MBool p1_inside = maia::geom::isPointInsideSphere<nDim>(line, m_center, m_radius);

  if(p1_inside) return 1;

  // check end of line
  const MBool p2_inside = maia::geom::isPointInsideSphere<nDim>(line + nDim, m_center, m_radius);

  if(p2_inside)
    return 1;
  else
    return 2;
}

// GeometryRoot functions

GeometryRoot::GeometryRoot(const MInt noSolvers, const MInt nDim_, const MPI_Comm comm)
  : GeometryBase(comm), m_nDim(nDim_) {
  if(nDim_ == 2)
    initGeometry<2>(noSolvers);
  else
    initGeometry<3>(noSolvers);
}

template <MInt nDim>
void GeometryRoot::initGeometry(const MInt noSolvers) {
  // how are the surfaces distributed among the nodes
  // right now, every node has one surface
  m_distribution.noNodes = noSolvers;
  mAlloc(m_distribution.noSurfacesPerNode, m_distribution.noNodes, AT_, 1, "noSurfacesPerSolver");
  mAlloc(m_distribution.nodeSurfaceIds, m_distribution.noNodes, m_distribution.noSurfacesPerNode, AT_,
         "noSurfacesPerNode");
  m_noSolverSurfaces = 0;
  for(MInt node = 0; node < m_distribution.noNodes; node++) {
    for(MInt srfc = 0; srfc < m_distribution.noSurfacesPerNode[node]; srfc++) {
      m_distribution.nodeSurfaceIds[node][srfc] = node;
      m_noSolverSurfaces++;
    }
  }

  // what kind of geometry/surface to expect
  mAlloc(m_solverSurfaceType, m_noSolverSurfaces, AT_, 0, "solverSurfaceType");
  for(MInt node = 0; node < m_distribution.noNodes; node++) {
    for(MInt srfc = 0; srfc < m_distribution.noSurfacesPerNode[node]; srfc++) {
      const MInt srfcId = m_distribution.nodeSurfaceIds[node][srfc];
      MString surfaceType = "STL";
      surfaceType = Context::getSolverProperty<MString>("solverSurfaceType", node, AT_, &surfaceType, srfc);
      m_solverSurfaceType[srfcId] = string2enum(surfaceType);
    }
  }

  // init
  m_solverSurface = new SolverSurface*[m_noSolverSurfaces];
  for(MInt i = 0; i < m_noSolverSurfaces; i++) {
    switch(m_solverSurfaceType[i]) {
      case STL: {
        m_solverSurface[i] = new SolverSTLSurface<nDim>(m_mpiComm, i);
        break;
      }
      case ANALYTIC_BOX: {
        // here we hardcode the box, just a usefull example
        MFloat bbox[2 * nDim];
        for(MInt dir = 0; dir < 2 * nDim; dir++) {
          bbox[dir] = Context::getSolverProperty<MFloat>("analyticSolverBoundingBox", i, AT_, dir);
        }
        m_solverSurface[i] = new SolverAnalyticBoxSurface<nDim>(bbox);
        break;
      }
      case ANALYTIC_SPHERE: {
        // here we hardcode the sphere, just a usefull example
        MFloat center[nDim];
        for(MInt dir = 0; dir < nDim; dir++) {
          center[dir] = Context::getSolverProperty<MFloat>("analyticSolverBoundingSphere", i, AT_, dir);
        }
        const MFloat radius = Context::getSolverProperty<MFloat>("analyticSolverBoundingSphere", i, AT_, nDim);
        m_solverSurface[i] = new SolverAnalyticSphereSurface<nDim>(center, radius);
        break;
      }
      default: {
        mTerm(1, AT_, "unknown surface type");
      }
    }
  }
}

MInt GeometryRoot::noSegmentsOfNode(const MInt node) {
  return m_solverSurface[node]->noSegments(); // temporary since there wil be several surfaces per node
}

void GeometryRoot::boundingBox(MFloat* const bBox) {
  // init temporary bounding box
  MFloat tmpBB[6];
  // init bounding box
  for(MInt dim = 0; dim < m_nDim; dim++) {
    bBox[dim] = numeric_limits<MFloat>::max();
    bBox[dim + m_nDim] = numeric_limits<MFloat>::lowest();
  }
  // check all surfaces for their extend
  for(MInt i = 0; i < m_noSolverSurfaces; i++) {
    // init temporary bounding box
    for(MInt dim = 0; dim < m_nDim; dim++) {
      tmpBB[dim] = numeric_limits<MFloat>::max();
      tmpBB[dim + m_nDim] = numeric_limits<MFloat>::lowest();
    }
    // get bounding box of surface
    m_solverSurface[i]->boundingBox(tmpBB);
    // get extrem values
    for(MInt dim = 0; dim < m_nDim; dim++) {
      bBox[dim] = mMin(bBox[dim], tmpBB[dim]);
      bBox[dim + m_nDim] = mMax(bBox[dim + m_nDim], tmpBB[dim + m_nDim]);
    }
  }
}

void GeometryRoot::boundingBoxOfNode(MFloat* const bBox, const MInt node) {
  // init temporary bounding box
  MFloat tmpBB[6];
  // init bounding box
  for(MInt dim = 0; dim < m_nDim; dim++) {
    bBox[dim] = numeric_limits<MFloat>::max();
    bBox[dim + m_nDim] = numeric_limits<MFloat>::lowest();
  }
  // check all surfaces for their extend
  for(MInt i = 0; i < m_distribution.noSurfacesPerNode[node]; i++) {
    const MInt srfcId = m_distribution.nodeSurfaceIds[node][i];
    // init temporary bounding box
    for(MInt dim = 0; dim < m_nDim; dim++) {
      tmpBB[dim] = numeric_limits<MFloat>::max();
      tmpBB[dim + m_nDim] = numeric_limits<MFloat>::lowest();
    }
    // get bounding box of surface
    m_solverSurface[srfcId]->boundingBox(tmpBB);
    // get extrem values
    for(MInt dim = 0; dim < m_nDim; dim++) {
      bBox[dim] = mMin(bBox[dim], tmpBB[dim]);
      bBox[dim + m_nDim] = mMax(bBox[dim + m_nDim], tmpBB[dim + m_nDim]);
    }
  }
}

MBool GeometryRoot::isPointInside(const MFloat* const point) {
  for(MInt node = 0; node < m_distribution.noNodes; node++) {
    if(isPointInsideNode(point, node)) {
      return true;
    }
  }
  return false;
}

MBool GeometryRoot::isPointInsideNode(const MFloat* const point, const MInt node) {
  MFloat bbox[2 * MAX_SPACE_DIMENSIONS];
  MFloat line[2 * MAX_SPACE_DIMENSIONS];
  MFloat bbExtend[MAX_SPACE_DIMENSIONS];

  boundingBoxOfNode(bbox, node);

  for(MInt dim = 0; dim < m_nDim; dim++) {
    line[dim] = point[dim];
    bbExtend[dim] = bbox[dim + m_nDim] - bbox[dim];
  }

  MBool retVal = true;
  for(MInt rayDim = 0; rayDim < m_nDim; rayDim++) { // why in all dimensions?
    for(MInt dim = 0; dim < m_nDim; dim++) {
      line[dim + m_nDim] = point[dim];
    }
    line[rayDim + m_nDim] += 2 * bbExtend[rayDim];
    MInt noIntersec = 0;
    for(MInt i = 0; i < m_distribution.noSurfacesPerNode[node]; i++) {
      MInt srfcId = m_distribution.nodeSurfaceIds[node][i];
      noIntersec += m_solverSurface[srfcId]->countLineIntersectingElements(line);
    }
    retVal = retVal && ((noIntersec % 2) == 1);
    if(!retVal) return false;
  }
  return true;
}

MBool GeometryRoot::getCellIntersectingSurfaces(const MFloat* const coords, const MFloat cellHalfLength,
                                                MBool* const* const cutInfo) {
  MBool retVal = false;
  for(MInt node = 0; node < m_distribution.noNodes; node++) {
    const MBool newCut = getCellIntersectingSurfacesOfNode(coords, cellHalfLength, cutInfo[node], node);
    retVal = retVal || newCut;
  }

  return retVal;
}

MBool GeometryRoot::getCellIntersectingSurfacesOfNode(const MFloat* const coords, const MFloat cellHalfLength,
                                                      MBool* const cutInfo, const MInt node) {
  MBool retVal = false;
  for(MInt i = 0; i < m_distribution.noSurfacesPerNode[node]; i++) {
    const MInt srfcId = m_distribution.nodeSurfaceIds[node][i];
    const MInt cutOffset = 0; // needs to be set for the m_solverSurface to know where to set intersections
                              // but this depends on the concept: How many elements, i.e., different STL surfaces,
                              // can an adt contain?
    retVal =
        retVal || m_solverSurface[srfcId]->getCellIntersectingElements(coords, cellHalfLength, cutInfo + cutOffset);
  }

  return retVal;
}

template <MInt nDim>
void GeometryRoot::setGeometryPointerToNode(Geometry<nDim>*& geometryPointer, const MInt node) {
  m_solverSurface[node]->setGeometryPointer(geometryPointer);
}

// geometry node

GeometryNode::GeometryNode(GeometryBase* geometryManager, MInt node, MPI_Comm comm) {
  cout << "SolverGeometry constructor" << endl;
  m_geometry = geometryManager;
  m_node = node;
  m_mpiComm = comm;
}

// Explicit instantiations for 2D and 3D
template class SolverSTLSurface<2>;
template class SolverSTLSurface<3>;
template class SolverAnalyticBoxSurface<2>;
template class SolverAnalyticBoxSurface<3>;
template void GeometryRoot::initGeometry<2>(MInt noSolvers);
template void GeometryRoot::initGeometry<3>(MInt noSolvers);
template void GeometryRoot::setGeometryPointerToNode<2>(Geometry<2>*& geometryPointer, MInt noSolvers);
template void GeometryRoot::setGeometryPointerToNode<3>(Geometry<3>*& geometryPointer, MInt noSolvers);
