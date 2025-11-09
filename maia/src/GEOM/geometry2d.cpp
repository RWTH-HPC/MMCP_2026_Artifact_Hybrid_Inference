// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometry2d.h"

#include <bitset>
#include "geometryadt.h"
#include "geometrycontext.h"
#include "geometryelement.h"
#include "globals.h"

using namespace std;

Geometry2D::Geometry2D(const MInt solverId_, const MPI_Comm comm) : Geometry<2>(solverId_, comm) {
  TRACE();

  m_adt = 0;
  m_elements = 0;
  m_noElements = 0;

  // moving boundary
  m_mbelements = 0;
  m_noMBElements = 0;
  m_GFieldInitFromSTL = 0;
  m_levelSetIntfBndId = 0;
  m_noLevelSetIntfBndIds = 0;
  m_GFieldInitFromSTL = Context::getSolverProperty<MBool>("GFieldInitFromSTL", solverId(), AT_, &m_GFieldInitFromSTL);

  if(m_GFieldInitFromSTL) {
    m_levelSetIntfBndId = Context::getSolverProperty<MInt>("levelSetIntfBndId", solverId(), AT_, &m_levelSetIntfBndId);
    if(Context::propertyExists("bodyBndryCndIds", solverId())
       || Context::propertyExists("GFieldInitFromSTLBndCndIds", solverId())) {
      const MInt noBodyBndIds = Context::propertyLength("bodyBndryCndIds", solverId());
      const MInt noBndIds = Context::propertyLength("GFieldInitFromSTLBndCndIds", solverId());

      m_noLevelSetIntfBndIds = noBodyBndIds + noBndIds;
      m_levelSetIntfBndIds = new MInt[m_noLevelSetIntfBndIds];

      for(MInt i = 0; i < noBodyBndIds; i++) {
        m_levelSetIntfBndIds[i] = Context::getSolverProperty<MInt>("bodyBndryCndIds", solverId(), AT_, i);
        if(domainId() == 0) {
          cerr << " levelSetIntfBndIds: " << m_levelSetIntfBndIds[i] << endl;
        }
      }

      for(MInt i = noBodyBndIds; i < m_noLevelSetIntfBndIds; i++) {
        m_levelSetIntfBndIds[i] =
            Context::getSolverProperty<MInt>("GFieldInitFromSTLBndCndIds", solverId(), AT_, i - noBodyBndIds);
        if(domainId() == 0) {
          cerr << " levelSetIntfBndIds: " << m_levelSetIntfBndIds[i] << endl;
        }
      }

      if(domainId() == 0) {
        cerr << " noLevelSetIntfBndIds: " << m_noLevelSetIntfBndIds << endl;
      }

    } else {
      m_levelSetIntfBndIds = nullptr;
      //       for(MInt i=0; i < m_noLevelSetIntfBndIds; i++){
      //         m_levelSetIntfBndIds[i] = m_levelSetIntfBndId;
      //       }
    }
  }

  MString testcaseDir = "./";
  testcaseDir = Context::getSolverProperty<MString>("testcaseDir", solverId(), AT_, &testcaseDir);
  MString inputDir = "./";
  inputDir = Context::getSolverProperty<MString>("inputDir", solverId(), AT_, &inputDir);

  MString tmpFileName =
      testcaseDir + inputDir + Context::getSolverProperty<MString>("geometryInputFileName", solverId(), AT_);

  if(tmpFileName.find(".toml") == MString::npos) {
    geometryContext().readPropertyFile(NETCDF, tmpFileName.c_str());
  } else {
    geometryContext().readPropertyFile(TOML, tmpFileName.c_str());
  }
  m_bodyMap = geometryContext().getBodies();
  readSegments();
  m_adt = new GeometryAdt<2>(this);
  m_adt->buildTree();

  // moving boundary
  if(m_GFieldInitFromSTL) {
    m_adt->buildTreeMB();
  }
}

/*
 * brief: generates an auxiliary geometry from an ASCII file that can be used for whatever
 *
 * Thomas Schilden, December 2016
 */
Geometry2D::Geometry2D(const MInt solverId_, const MString filename, const MPI_Comm comm) : Geometry(solverId_, comm) {
  TRACE();

  m_log << "reading the auxiliary stl " << filename << endl;

  m_GFieldInitFromSTL = 0;

  m_noElements = 0;
  m_noMBElements = 0;

  countSegmentLinesASCII(filename, &m_noElements);

  m_log << "  number of Elements: " << m_noElements << endl;

  m_elements = new Collector<element<nDim>>(m_noElements, nDim, 0);

  MInt bla = 0;
  readSegmentLinesASCII(filename, m_elements, 0, 0, &bla);

  elements = &(m_elements->a[0]);

  Geometry2D::calculateBoundingBox();

  m_adt = new GeometryAdt<2>(this);
  m_adt->buildTree();
}

/** \brief
 *
 *
 */
Geometry2D::~Geometry2D() {
  TRACE();

  delete m_adt;
  delete m_elements;

  if(m_mbelements) delete m_mbelements;
}
/** \brief
 *
 *
 */
void Geometry2D::readSegments() {
  TRACE();

  MInt counter = 0;
  MInt segmentId = 0;
  MString fileName;
  // This loop only counts all elements of all segments
  for(m_bodyIt = m_bodyMap.begin(); m_bodyIt != m_bodyMap.end(); m_bodyIt++) {
    // Do not create default body!
    if(m_bodyIt->second->name != "default") {
      for(MInt i = 0; i < m_bodyIt->second->noSegments; i++) {
        fileName = *geometryContext().getProperty("filename", segmentId)->asString();
        if(geometryContext().noPropertySegments("filename") == 1
           && (m_bodyIt->second->noSegments > 1 || m_bodyMap.size() > 2)) {
          fileName += "." + to_string(segmentId);
        }
        countSegmentLinesASCII(fileName, &counter);
        m_noElements += counter;
        segmentId++;
      }
    }
  }

  auto* m_allelements = new Collector<element<nDim>>(m_noElements, nDim, 0);
  element<nDim>* allelements = m_allelements->a;
  vector<MInt> allBCs;

  segmentId = 0;
  counter = 0;
  // This loop reads all elements from all segments
  for(m_bodyIt = m_bodyMap.begin(); m_bodyIt != m_bodyMap.end(); m_bodyIt++) {
    // Do not create default body!
    if(m_bodyIt->second->name != "default") {
      for(MInt i = 0; i < m_bodyIt->second->noSegments; i++) {
        m_log << " Reading Ascii coordinates from file " << endl;
        fileName = *geometryContext().getProperty("filename", segmentId)->asString();
        if(geometryContext().noPropertySegments("filename") == 1
           && (m_bodyIt->second->noSegments > 1 || m_bodyMap.size() > 2)) {
          fileName += "." + to_string(segmentId);
        }
        MInt bndCndId = *geometryContext().getProperty("BC", segmentId)->asInt();
        allBCs.push_back(bndCndId);
        readSegmentLinesASCII(fileName, m_allelements, bndCndId, segmentId, &counter);
        segmentId++;
      }
    }
  }
  m_log << m_noElements << " line elements read. " << endl;

  if(m_GFieldInitFromSTL) {
    m_mbelements = new Collector<element<nDim>>(m_noMBElements, nDim, 0);
    m_noElements -= m_noMBElements;
    mbelements = m_mbelements->a;
  }

  m_elements = new Collector<element<nDim>>(m_noElements, nDim, 0);
  elements = m_elements->a;

  MInt mbelem_counter = 0, elem_counter = 0;
  MBool mbElem = false;

  for(MInt allelem_counter = 0; allelem_counter < m_allelements->size(); allelem_counter++) {
    if(m_GFieldInitFromSTL && m_noLevelSetIntfBndIds > 0) {
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++) {
        if(allelements[allelem_counter].m_bndCndId == m_levelSetIntfBndIds[id]) {
          mbElem = true;
          break;
        }
      }
      if(mbElem) {
        // LevelSet Moving boundary elements
        m_mbelements->append();
        for(MInt j = 0; j < 2; j++) {
          mbelements[mbelem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 2; i++) {
            mbelements[mbelem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        mbelements[mbelem_counter].boundingBox();
        mbelements[mbelem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        mbelements[mbelem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        mbelem_counter++;
        mbElem = false;
      } else {
        // Non movable elements
        m_elements->append();
        for(MInt j = 0; j < 2; j++) {
          elements[elem_counter].m_normal[j] = allelements[allelem_counter].m_normal[j];
          for(MInt i = 0; i < 2; i++) {
            elements[elem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        elements[elem_counter].boundingBox();
        elements[elem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        elements[elem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        elem_counter++;
      }
    } else {
      if(m_GFieldInitFromSTL && (m_levelSetIntfBndId == allelements[allelem_counter].m_bndCndId)) {
        // LevelSet Moving boundary elements
        m_mbelements->append();
        for(MInt i = 0; i < 2; i++) {
          mbelements[mbelem_counter].m_normal[i] = allelements[allelem_counter].m_normal[i];
        }
        for(MInt j = 0; j < 2; j++) {
          for(MInt i = 0; i < 2; i++) {
            mbelements[mbelem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        mbelements[mbelem_counter].boundingBox();
        mbelements[mbelem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        mbelements[mbelem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        mbelem_counter++;
      } else {
        // Non movable elements
        m_elements->append();
        for(MInt i = 0; i < 2; i++) {
          elements[elem_counter].m_normal[i] = allelements[allelem_counter].m_normal[i];
        }
        for(MInt j = 0; j < 2; j++) {
          for(MInt i = 0; i < 2; i++) {
            elements[elem_counter].m_vertices[j][i] = allelements[allelem_counter].m_vertices[j][i];
          }
        }
        elements[elem_counter].boundingBox();
        elements[elem_counter].m_bndCndId = allelements[allelem_counter].m_bndCndId;
        elements[elem_counter].m_segmentId = allelements[allelem_counter].m_segmentId;
        elem_counter++;
      }
    }
  }

  delete m_allelements;

  calculateBoundingBox();

  m_noAllBCs = allBCs.size();
  mAlloc(m_allBCs, m_noAllBCs, "m_allBCs", 0, AT_);
  for(MInt i = 0; i < m_noAllBCs; i++)
    m_allBCs[i] = allBCs[i];
}

/** \brief Detemines whether a triangle and a gridcell intersect with the separating axis theorem (SAT)
    \author Andreas Lintermann
    \date 09.11.2009

    The separating axis theorem (SAT) states, that two convex polyhedra, that are to be tested against intersection are
    disjoint, if they can be separated along either an axis parallel to a normal of a face of the first or the second
   polyhedra, or along an axis formed from the cross product of an edge from them.

    <b>Algorithm:</b>
    <ul>
    <li>
    <b>1)</b> Test projection onto \f$x,y\f$-axis
    Project the triangle and the gridcell onto both axes and check overlapping
    </li>

    <li>
    <b>2)</b> Test if we have parallel lines. If this is the case we are done
    </li>

    <li>
    <b>3)</b> Otherwise, rotate everything so that the axis to be tested becomes axis-aligned. Then test again agains
   projection overlapping.
    </li>
    </ul>

    This function replaces the other Geometry2D::getIntersectionElements, which cannot be used on IBM_BLUE_GENE.
 */
MInt Geometry2D::getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList, MFloat /*cellHalfLength*/,
                                         const MFloat* const /*cell_coords*/) {
  //  TRACE();


  MInt noReallyIntersectingNodes = 0;
  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat cos_alpha, sin_alpha, length_diff;
  MFloat diff_vec[nDim];

  for(MInt i = 0; i < nDim; i++) {
    diff_vec[i] = numeric_limits<MFloat>::max(); // gcc 4.8.2 maybe uninitialized
  }

  MFloat rot_point;
  MFloat rot_corners[nDim];

  MBool overlap[2];
  MBool parallel;
  MInt parallel_dim;

  for(MInt n = 0; n < noNodes; n++) {
    length_diff = 0.0;
    parallel = false;
    overlap[0] = false;
    overlap[1] = false;

    parallel_dim = 0;

    // coordinate-axis tests
    for(MInt i = 0; i < nDim; i++) // run over dimensions (x and y)
    {
      // is first point in projection of quad?
      if(elements[nodeList[n]].m_vertices[0][i] >= targetRegion[i]
         && elements[nodeList[n]].m_vertices[0][i] <= targetRegion[i + nDim]) {
        overlap[i] = true;
        continue;
      }

      // is second point in projection of quad?
      if(elements[nodeList[n]].m_vertices[1][i] >= targetRegion[i]
         && elements[nodeList[n]].m_vertices[1][i] <= targetRegion[i + nDim]) {
        overlap[i] = true;
        continue;
      }

      // do both points clasp the projection of quad?
      if(elements[nodeList[n]].m_vertices[0][i] <= targetRegion[i]
         && elements[nodeList[n]].m_vertices[1][i] >= targetRegion[i + nDim]) {
        overlap[i] = true;
        continue;
      }
      if(elements[nodeList[n]].m_vertices[1][i] <= targetRegion[i]
         && elements[nodeList[n]].m_vertices[0][i] >= targetRegion[i + nDim]) {
        overlap[i] = true;
        continue;
      }
      if(!overlap[i]) break;
    }

    if(!(overlap[0] && overlap[1])) continue;

    // is the line parallel to axis? Then we are done if we don't
    for(parallel_dim = 0; parallel_dim < nDim; parallel_dim++) {
      parallel = parallel
                 || (approx(elements[nodeList[n]].m_vertices[0][parallel_dim],
                            elements[nodeList[n]].m_vertices[1][parallel_dim], MFloatEps));
      if(parallel) break;
    }

    if(parallel) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    } else {
      // For the sake of simplicity rotate everything so that line becomes parallel to x-axis. Check overlapping in
      // y-direction.
      for(MInt i = 0; i < nDim; i++) {
        diff_vec[i] = elements[nodeList[n]].m_vertices[1][i] - elements[nodeList[n]].m_vertices[0][i];
        length_diff += diff_vec[i] * diff_vec[i];
      }
      length_diff = sqrt(length_diff);

      cos_alpha = fabs(diff_vec[0]) / length_diff;
      sin_alpha = fabs(diff_vec[1]) / length_diff;

      if((diff_vec[0] * diff_vec[1]) < 0) // negative slope
      {
        rot_point =
            sin_alpha * elements[nodeList[n]].m_vertices[0][0] + cos_alpha * elements[nodeList[n]].m_vertices[0][1];
        rot_corners[0] = sin_alpha * (targetRegion[2]) + cos_alpha * (targetRegion[3]);
        rot_corners[1] = sin_alpha * (targetRegion[0]) + cos_alpha * (targetRegion[1]);
      } else // positive slope
      {
        rot_point = -1 * sin_alpha * elements[nodeList[n]].m_vertices[0][0]
                    + cos_alpha * elements[nodeList[n]].m_vertices[0][1];
        rot_corners[0] = -1 * sin_alpha * (targetRegion[0]) + cos_alpha * (targetRegion[3]);
        rot_corners[1] = -1 * sin_alpha * (targetRegion[2]) + cos_alpha * (targetRegion[1]);
      }

      // final cut test
      if(rot_corners[0] >= rot_point && rot_corners[1] <= rot_point) {
        nodeList[noReallyIntersectingNodes] = nodeList[n];
        noReallyIntersectingNodes++;
      }
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}

/** \brief Determines all elements that are inside or intersect the target region
    \author Rainhill Freitas
    \date unknown

    \bug labels:GEOM This function cannot be used on IBM_BLUE_GENE. It uses Cramer's rule to solve
   a 2 by 2 linear equation system to find intersection points. The determinat in the denominator
   of Cramer's rule can be become 0 which cannot be handled by IBM_BLUE_GENE.
   Use the other Geometry2D::getIntersectionElements instead.
  */
MInt Geometry2D::getIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  //  TRACE();

  MInt noReallyIntersectingNodes = 0;
  // get all candidates for an intersection
  // Check for intersection...
  bitset<4> points[2];
  bitset<4> faceCodes[4];
  bitset<4> pCode;

  // Edges of the targetRegion (using targetPoints)
  MInt rejection;
  MBool piercePointInside;
  MBool triviallyAccepted;

  // Each of the following arrays holds one different point for
  // all of the 6 planes. Points are built with the targetRegion
  // i.e. a = targetRegion[pointAInPlane[0]],
  //          targetRegion[pointAInPlane[1]],
  // etc.
  MInt pointAInPlane[4][2] = {{0, 1}, {2, 1}, {0, 1}, {0, 3}};

  MInt pointBInPlane[4][2] = {{0, 3}, {2, 3}, {2, 1}, {2, 3}};

  MFloat a[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat b[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat c[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // start of piercing edge
  MFloat d[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // end of piercing edge
  MFloat s1, s2, gamma;                          // For pierce point calculation
  MFloat pP[2] = {numeric_limits<MFloat>::max(),
                  numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // piercePoint
  faceCodes[0] = IPOW2(0);
  faceCodes[1] = IPOW2(1);
  faceCodes[2] = IPOW2(2);
  faceCodes[3] = IPOW2(3);
  bitset<4> result;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  //  return noNodes;
  noReallyIntersectingNodes = 0;

  for(MInt n = 0; n < noNodes; n++) {
    // create edges (in 3D) AB, AC, BC and point A B C
    // Determine outcode (see Aftosmis, Solution Adaptive Cartesian Grid Methods for Aerodynamic flows ...)

    // Loop over all points of an element<2>
    for(MInt p = 0; p < nDim; p++) {
      points[p] = 0;
      // Calculate outcode for point
      for(MInt j = 0; j < nDim; j++) {
        if(elements[nodeList[n]].m_vertices[p][j] < targetRegion[j]) {
          points[p] |= faceCodes[2 * j];
        }
        if(elements[nodeList[n]].m_vertices[p][j] > targetRegion[j + nDim]) {
          points[p] |= faceCodes[2 * j + 1];
        }
      }
      //      m_log << points[p] << endl;
    }
    rejection = 0;
    // check outcode combinations for edges for trivial rejection
    for(MInt i = 0; i < nDim; i++) {
      if((points[0] & points[1]) != 0) {
        rejection++;
        break;
      } else {
        // If one point is inside region the element<2> is trivially accepted
        triviallyAccepted = false;
        for(MInt k = 0; k < nDim; k++) {
          if(points[k] == 0) {
            triviallyAccepted = true;
            rejection = 0;
            break;
          }
        }
        if(triviallyAccepted) {
          break;
        }
        // No trivial rejection, check for rejection of subsegment:
        // For all pierce points!
        // 1. Calculate pierce point:
        //    a - determine plane for pierce point calculation
        //    b - calculate pierce point
        // 2. Check for rejection of new segment
        //    a - calculate new outcode
        //    b - check for containment in pierce planes face
        // 3. If all(!) pierce points are rejected -> reject edge
        // TODO labels:GEOM This algorith might get a problem if a triangle lies completely on
        //    a face.

        // 1.a
        result = (points[0] | points[1]);
        piercePointInside = false;
        for(MInt j = 0; j < 2 * nDim; j++) {
          if(result[j] == 1) {
            // pierce plane found
            for(MInt k = 0; k < nDim; k++) {
              a[k] = targetRegion[pointAInPlane[j][k]];
              b[k] = targetRegion[pointBInPlane[j][k]];
              c[k] = elements[nodeList[n]].m_vertices[0][k];
              d[k] = elements[nodeList[n]].m_vertices[1][k];
            }
            gamma = (b[0] - a[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - a[1]);

            s1 = ((c[0] - d[0]) * (a[1] - c[1]) - (c[1] - d[1]) * (a[0] - c[0])) / gamma;


            s2 = ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / gamma;
            // 1. b Pierce point pP in plane j:
            for(MInt k = 0; k < nDim; k++) {
              if(s1 * s1 < s2 * s2)
                pP[k] = c[k] + s2 * (d[k] - c[k]);
              else
                pP[k] = a[k] + s1 * (b[k] - a[k]);
            }
            pCode = 0;
            // 2. a Calculate outcode for pierce point
            for(MInt k = 0; k < nDim; k++) {
              if(pP[k] < targetRegion[k]) {
                pCode |= faceCodes[2 * k];
              }
              if(pP[k] > targetRegion[k + nDim]) {
                pCode |= faceCodes[2 * k + 1];
              }
            }

          } else {
            continue;
          }
          // 2. b
          result = faceCodes[j];
          result.flip();
          result = (result & pCode);
          if(result == 0) { // -> is contained
            piercePointInside = true;
            break;
          }
        }
        // reject if all pierce points are off coresponding face
        // else accept
        if(!piercePointInside) {
          rejection++;
        } else {
          rejection = 0;
          break;
        }
      }
    }
    // If not all edges are rejected a cutting element<2> has been found
    //    m_log << rejection << endl;
    if(!rejection) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }

    // write nodes in reallyIntersectingNodes
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}

/** \brief Determine intersection between an edge and a triangle (in 2D
 *        between 2 edges)
 *
 */
inline MBool Geometry2D::edgeTriangleIntersection(MFloat* trianglePoint1,
                                                  MFloat* trianglePoint2,
                                                  MFloat* /*trianglePoint3*/,
                                                  MFloat* edgePoint1,
                                                  MFloat* edgePoint2) {
  //  TRACE();

  MFloat a[2]; // point in plane
  MFloat b[2]; // point in plane
  MFloat c[2]; // start of piercing edge
  MFloat d[2]; // end of piercing edge
  MFloat s1, s2, gamma;


  // Calculate the intersection between edge and triangle and check if
  // the intersection point lies on the edge.

  // pierce plane found
  // TODO labels:GEOM why is a,b,c,d used here?
  for(MInt k = 0; k < nDim; k++) {
    a[k] = edgePoint1[k];
    b[k] = edgePoint2[k];
    c[k] = trianglePoint1[k];
    d[k] = trianglePoint2[k];
  }

  gamma = (b[0] - a[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - a[1]);

  // TODO labels:GEOM,toenhance the small number should be replaced by an eps value defined in maia.h
#ifdef IBM_BLUE_GENE
  if(gamma < 0.00000000001) {
    return false;
  }
#endif

  s1 = ((c[0] - d[0]) * (a[1] - c[1]) - (c[1] - d[1]) * (a[0] - c[0])) / gamma;

  s2 = ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / gamma;

  if(s2 < 0 || s2 > 1.0 || s1 < 0 || s1 > 1.0) {
    return false;
  } else {
    return true;
  }
}

/** \brief Returns the ids of all elements, cut by a orthogonal line (or rectangular region)
 *
 *
 */
MInt Geometry2D::getLineIntersectionElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  //  TRACE();

  MInt noReallyIntersectingNodes = 0;

  m_adt->retrieveNodes(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat e1[2], e2[2];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
  }

  for(MInt n = 0; n < noNodes; n++) {
    if(edgeTriangleIntersection(elements[nodeList[n]].m_vertices[0], elements[nodeList[n]].m_vertices[1], 0, e1, e2)) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}

MInt Geometry2D::getLineIntersectionElementsOld1(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  return getLineIntersectionElements(targetRegion, nodeList);
}

MInt Geometry2D::getLineIntersectionElementsOld2(MFloat* targetRegion, MInt* /*spaceDirection*/,
                                                 std::vector<MInt>& nodeList) {
  return getLineIntersectionElements(targetRegion, nodeList);
}

MInt Geometry2D::getIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  //   TRACE();

  MInt noReallyIntersectingNodes = 0;
  // get all candidates for an intersection
  // Check for intersection...
  bitset<4> points[2];
  bitset<4> faceCodes[4];
  bitset<4> pCode;

  // Edges of the targetRegion (using targetPoints)
  MInt rejection;
  MBool piercePointInside;
  MBool triviallyAccepted;

  // Each of the following arrays holds one different point for
  // all of the 6 planes. Points are built with the targetRegion
  // i.e. a = targetRegion[pointAInPlane[0]],
  //          targetRegion[pointAInPlane[1]],
  // etc.
  MInt pointAInPlane[4][2] = {{0, 1}, {2, 1}, {0, 1}, {0, 3}};

  MInt pointBInPlane[4][2] = {{0, 3}, {2, 3}, {2, 1}, {2, 3}};

  MFloat a[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat b[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // point in plane
  MFloat c[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // start of piercing edge
  MFloat d[2] = {numeric_limits<MFloat>::max(),
                 numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // end of piercing edge
  MFloat s1, s2, gamma;                          // For pierce point calculation
  MFloat pP[2] = {numeric_limits<MFloat>::max(),
                  numeric_limits<MFloat>::max()}; // gcc 4.8.2 maybe uninitialized // piercePoint
  faceCodes[0] = IPOW2(0);
  faceCodes[1] = IPOW2(1);
  faceCodes[2] = IPOW2(2);
  faceCodes[3] = IPOW2(3);
  bitset<4> result;

  m_adt->retrieveNodesMBElements(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();
  //  return noNodes;
  noReallyIntersectingNodes = 0;

  for(MInt n = 0; n < noNodes; n++) {
    // create edges (in 3D) AB, AC, BC and point A B C
    // Determine outcode (see Aftosmis, Solution Adaptive Cartesian Grid Methods for Aerodynamic flows ...)

    // Loop over all points of an element<2>
    for(MInt p = 0; p < nDim; p++) {
      points[p] = 0;
      // Calculate outcode for point
      for(MInt j = 0; j < nDim; j++) {
        if(mbelements[nodeList[n]].m_vertices[p][j] < targetRegion[j]) {
          points[p] |= faceCodes[2 * j];
        }
        if(mbelements[nodeList[n]].m_vertices[p][j] > targetRegion[j + nDim]) {
          points[p] |= faceCodes[2 * j + 1];
        }
      }
      //      m_log << points[p] << endl;
    }
    rejection = 0;
    // check outcode combinations for edges for trivial rejection
    for(MInt i = 0; i < nDim; i++) {
      if((points[0] & points[1]) != 0) {
        rejection++;
        break;
      } else {
        // If one point is inside region the element<2> is trivially accepted
        triviallyAccepted = false;
        for(MInt k = 0; k < nDim; k++) {
          if(points[k] == 0) {
            triviallyAccepted = true;
            rejection = 0;
            break;
          }
        }
        if(triviallyAccepted) {
          break;
        }
        // No trivial rejection, check for rejection of subsegment:
        // For all pierce points!
        // 1. Calculate pierce point:
        //    a - determine plane for pierce point calculation
        //    b - calculate pierce point
        // 2. Check for rejection of new segment
        //    a - calculate new outcode
        //    b - check for containment in pierce planes face
        // 3. If all(!) pierce points are rejected -> reject edge
        // TODO labels:GEOM This algorith might get a problem if a triangle lies completely on
        //    a face.

        // 1.a
        result = (points[0] | points[1]);
        piercePointInside = false;
        for(MInt j = 0; j < 2 * nDim; j++) {
          if(result[j] == 1) {
            // pierce plane found
            for(MInt k = 0; k < nDim; k++) {
              a[k] = targetRegion[pointAInPlane[j][k]];
              b[k] = targetRegion[pointBInPlane[j][k]];
              c[k] = mbelements[nodeList[n]].m_vertices[0][k];
              d[k] = mbelements[nodeList[n]].m_vertices[1][k];
            }
            gamma = (b[0] - a[0]) * (c[1] - d[1]) - (c[0] - d[0]) * (b[1] - a[1]);

            s1 = ((c[0] - d[0]) * (a[1] - c[1]) - (c[1] - d[1]) * (a[0] - c[0])) / gamma;


            s2 = ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / gamma;
            // 1. b Pierce point pP in plane j:
            for(MInt k = 0; k < nDim; k++) {
              if(s1 * s1 < s2 * s2)
                pP[k] = c[k] + s2 * (d[k] - c[k]);
              else
                pP[k] = a[k] + s1 * (b[k] - a[k]);
            }
            pCode = 0;
            // 2. a Calculate outcode for pierce point
            for(MInt k = 0; k < nDim; k++) {
              if(pP[k] < targetRegion[k]) {
                pCode |= faceCodes[2 * k];
              }
              if(pP[k] > targetRegion[k + nDim]) {
                pCode |= faceCodes[2 * k + 1];
              }
            }

          } else {
            continue;
          }
          // 2. b
          result = faceCodes[j];
          result.flip();
          result = (result & pCode);
          if(result == 0) { // -> is contained
            piercePointInside = true;
            break;
          }
        }
        // reject if all pierce points are off coresponding face
        // else accept
        if(!piercePointInside) {
          rejection++;
        } else {
          rejection = 0;
          break;
        }
      }
    }
    // If not all edges are rejected a cutting element<2> has been found
    //    m_log << rejection << endl;
    if(!rejection) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
      continue;
    }

    // write nodes in reallyIntersectingNodes
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}

MInt Geometry2D::getLineIntersectionMBElements(MFloat* targetRegion, std::vector<MInt>& nodeList) {
  //   TRACE();

  MInt noReallyIntersectingNodes = 0;

  m_adt->retrieveNodesMBElements(targetRegion, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat e1[2], e2[2];
  for(MInt i = 0; i < nDim; i++) {
    e1[i] = targetRegion[i];
    e2[i] = targetRegion[i + nDim];
  }

  for(MInt n = 0; n < noNodes; n++) {
    if(edgeTriangleIntersection(mbelements[nodeList[n]].m_vertices[0], mbelements[nodeList[n]].m_vertices[1], 0, e1,
                                e2)) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}


MInt Geometry2D::getSphereIntersectionMBElements(MFloat* P, MFloat radius, std::vector<MInt>& nodeList) {
  // TRACE();

  MInt noReallyIntersectingNodes = 0;

  // compute minimum target region that contains the sphere:
  MFloat target[4] = {0, 0, 0, 0};
  MFloat enlargeFactor = 1.05;
  for(MInt i = 0; i < nDim; i++) {
    target[i] = P[i] - radius * enlargeFactor;
    target[i + nDim] = P[i] + radius * enlargeFactor;
  }

  // fetch all possibly intersecting triangles in this region:
  m_adt->retrieveNodesMBElements(target, nodeList);
  const MInt noNodes = nodeList.size();

  MFloat A[2], B[2], AB[2], Q1[2];

  // loop over all elements
  for(MInt n = 0; n < noNodes; n++) {
    // store the vertices of the current element<2> in A, B, C:
    for(MInt i = 0; i < nDim; i++) {
      A[i] = mbelements[nodeList[n]].m_vertices[0][i];
      B[i] = mbelements[nodeList[n]].m_vertices[1][i];
    }

    // compute separation algorithm from http://realtimecollisiondetection.net/blog/?p=103:
    for(MInt i = 0; i < nDim; i++) {
      A[i] = A[i] - P[i];
      B[i] = B[i] - P[i];
      AB[i] = B[i] - A[i];
    }
    MFloat rr = radius * radius;
    MFloat aa = F0;
    MFloat ab = F0;
    MFloat bb = F0;
    MFloat e1 = F0;
    for(MInt i = 0; i < nDim; i++) {
      aa += A[i] * A[i];
      ab += A[i] * B[i];
      bb += B[i] * B[i];
      e1 += AB[i] * AB[i];
    }
    MBool sep2 = (aa > rr) && (ab > aa);
    MBool sep3 = (bb > rr) && (ab > bb);
    MFloat d1 = ab - aa;
    MFloat qq1 = F0;
    for(MInt i = 0; i < nDim; i++) {
      Q1[i] = A[i] * e1 - d1 * AB[i];
      qq1 += Q1[i] * Q1[i];
    }
    MBool sep5 = (qq1 > rr * e1 * e1);

    MBool separated = sep2 || sep3 || sep5;

    if(!separated) {
      nodeList[noReallyIntersectingNodes] = nodeList[n];
      noReallyIntersectingNodes++;
    }
  }
  nodeList.resize(noReallyIntersectingNodes);

  return noReallyIntersectingNodes;
}


void Geometry2D::MoveAllMBElementVertex(MFloat* dx) {
  TRACE();

  for(MInt e = 0; e < m_noMBElements; e++) {
    for(MInt j = 0; j < nDim; j++) {
      for(MInt i = 0; i < nDim; i++) {
        mbelements[e].m_vertices[j][i] += dx[i];
      }
    }
  }
}

void Geometry2D::MoveMBElementVertex(MInt e, MInt v, MFloat* dx) {
  TRACE();

  for(MInt i = 0; i < nDim; i++) {
    mbelements[e].m_vertices[v][i] += dx[i];
  }
}

void Geometry2D::ReplaceMBElementVertex(MInt e, MInt v, MFloat* np) {
  TRACE();

  for(MInt i = 0; i < nDim; i++)
    mbelements[e].m_vertices[v][i] = np[i];
}

void Geometry2D::UpdateMBBoundingBox() {
  TRACE();

  for(MInt e = 0; e < m_mbelements->size(); e++) {
    mbelements[e].boundingBox();
  }
  for(MInt j = 0; j < nDim; j++) {
    m_mbminMax[j] = mbelements[0].m_minMax[j];
    m_mbminMax[j + nDim] = mbelements[0].m_minMax[j + nDim];
  }
  for(MInt i = 0; i < m_mbelements->size(); i++) {
    for(MInt j = 0; j < nDim; j++) {
      m_mbminMax[j + nDim] = (m_mbminMax[j + nDim] < mbelements[i].m_minMax[j + nDim])
                                 ? mbelements[i].m_minMax[j + nDim]
                                 : m_mbminMax[j + nDim];
      m_mbminMax[j] = (m_mbminMax[j] > mbelements[i].m_minMax[j]) ? mbelements[i].m_minMax[j] : m_mbminMax[j];
    }
  }
}

void Geometry2D::UpdateADT() {
  TRACE();

  m_adt->buildTreeMB();
}

/** \brief counts the number of lines in an ASCII file, should be the in the first line
 *
 * \author Thomas Schilden
 * \date 27.12.2016
 *
 * \param[in] fileName the name of the file to open
 * \param[in] the number of entries to be returned (gets nulled in function again)
 *
 **/
inline void Geometry2D::countSegmentLinesASCII(const MString& fileName, MInt* noElements) {
  TRACE();

  *noElements = 0;

  m_log << " Counting segment file: " << fileName << endl;

  ifstream ifl(fileName);
  if(!ifl) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment (counting loop), couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }
  // If number of elements unknown, count them...
  ifl >> (*noElements);
  ifl.close();
  (*noElements)--; // The ascii files contain the points(!) not the elements
  m_log << " number of vertices = " << *noElements << endl;
}

/** \brief reads the lines in an ASCII file
 *
 * \author Thomas Schilden
 * \date 27.12.2016
 *
 * \param[in] fileName the name of the file to open
 * \param[in] etc..
 *
 **/
inline void Geometry2D::readSegmentLinesASCII(MString fileName, Collector<element<2>>* elemCollector, MInt bndCndId,
                                              MInt segmentId, MInt* offset) {
  TRACE();
  m_log << " Reading segment file: " << fileName << endl;
  m_log << " BC : " << bndCndId << endl;


  ifstream ifl(fileName);
  if(!ifl) {
    stringstream errorMessage;
    errorMessage << " ERROR in segment::readSegment (reading loop), couldn't find file : " << fileName << "!";
    mTerm(1, AT_, errorMessage.str());
  }

  element<2>* allelements = elemCollector->a;

  // Read no of coordinates
  MInt tmp;
  ifl >> tmp;
  elemCollector->append();
  for(MInt j = 0; j < 2; j++) {
    ifl >> allelements[*offset].m_vertices[0][j];
    allelements[*offset].m_normal[j] = F0;
  }
  MInt fileElementCounter = 0;
  while(!ifl.eof()) {
    // Read coordinates
    for(MInt j = 0; j < 2; j++) {
      ifl >> allelements[*offset].m_vertices[1][j];
      allelements[*offset].m_normal[j] = F0;
    }
    allelements[*offset].boundingBox();
    allelements[*offset].m_bndCndId = bndCndId;
    allelements[*offset].m_segmentId = segmentId;
    if(m_GFieldInitFromSTL) {
      for(MInt id = 0; id < m_noLevelSetIntfBndIds; id++) {
        if(bndCndId == m_levelSetIntfBndIds[id]) {
          m_noMBElements++;
        }
      }
    }

    if(m_noElements > 1000 && !((*offset) % (m_noElements / 10))) {
      m_log << (MInt)(((*offset) - 1) / ((MFloat)m_noElements / 100)) + 1 << " % read." << endl;
    }
    fileElementCounter++;
    (*offset)++;
    // Check if already all points read
    if(fileElementCounter < tmp - 1) { // set start point of next line
      elemCollector->append();
      for(MInt j = 0; j < 2; j++) {
        allelements[*offset].m_vertices[0][j] = allelements[(*offset) - 1].m_vertices[1][j];
        allelements[*offset].m_normal[j] = F0;
      }
    } else { // If all points read leave loop
      break;
    }
  }
  ifl.close();
}

/** \brief Calculates the global bounding box of the geometry
 *
 * \author Thomas Schilden
 * \date 27.12.2016
 *
 **/
void Geometry2D::calculateBoundingBox() {
  TRACE();

  // Calculate bounding box of segment
  for(MInt j = 0; j < nDim; j++) {
    m_minMax[j] = elements[0].m_minMax[j];
    m_minMax[j + nDim] = elements[0].m_minMax[j + nDim];
  }
  for(MInt i = 0; i < m_noElements; i++) {
    for(MInt j = 0; j < nDim; j++) {
      m_minMax[j + nDim] =
          (m_minMax[j + nDim] < elements[i].m_minMax[j + nDim]) ? elements[i].m_minMax[j + nDim] : m_minMax[j + nDim];
      m_minMax[j] = (m_minMax[j] > elements[i].m_minMax[j]) ? elements[i].m_minMax[j] : m_minMax[j];
    }
  }
  m_log << " Bounding box for segment :" << endl;
  m_log << " max = ";
  for(MInt i = 0; i < nDim; i++) {
    m_log << m_minMax[i + nDim] << " ";
  }
  m_log << endl;

  m_log << " min = ";
  for(MInt i = 0; i < nDim; i++) {
    m_log << m_minMax[i] << " ";
  }
  m_log << endl;

  if(m_GFieldInitFromSTL) {
    // Calculate bounding box of Moving Boundary segment
    UpdateMBBoundingBox();
    m_log << " Bounding box for Moving Boundary segment:" << endl;
    m_log << " max = ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_mbminMax[i + nDim] << " ";
    }
    m_log << endl;
    m_log << " min = ";
    for(MInt i = 0; i < nDim; i++) {
      m_log << m_mbminMax[i] << " ";
    }
    m_log << endl;
  }
}
