// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometry.h"
#include "IO/context.h"
#include "MEMORY/collector.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "geometryadt.h"
#include "property.h"

using namespace std;

template <MInt nDim>
Geometry<nDim>::Geometry(const MInt solverId_, const MPI_Comm comm)
  : m_solverId(solverId_), m_mpiComm(comm), m_geometryContext(comm), m_haloElementOffset(0) {
  // Determine domain id and number of domains
  MPI_Comm_rank(mpiComm(), &m_domainId);
  MPI_Comm_size(mpiComm(), &m_noDomains);


  /*! \page propertyPage1
    \section parallelGeometry
    <code>MInt Geometry::m_parallelGeometry </code>\n
    default = <code>0</code>\n \n
    Test if point is inside or outside of stl.\n
    possible values are:
    <ul>
    <li>STD: Standard ?</li>
    <li>SAT: Separation Access Theorem</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, MASSIVE, TRIGGER, GEOMETRY</i>
  */
  // TODO labels:GEOM,DOC
  m_inOutTest = "perpOp";
  m_inOutTest = Context::getSolverProperty<MString>("inOutTest", m_solverId, AT_, &m_inOutTest);

  m_parallelGeometry = false;
  /*! \page propertyPage1
    \section parallelGeometry
    <code>MBool Geometry::m_parallelGeometry </code>\n
    default = <code>0</code>\n \n
    Trigger the use of parallel geometry.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, MASSIVE, TRIGGER, GEOMETRY</i>
  */
  m_parallelGeometry = Context::getSolverProperty<MBool>("parallelGeometry", m_solverId, AT_, &m_parallelGeometry);
  if(m_noDomains == 1) m_parallelGeometry = 0;

  m_debugParGeom = 0;
  /*! \page propertyPage1
    \section debugParGeom
    <code>MBool Geometry::m_debugParGeom</code>\n
    default = <code>0</code>\n \n
    Trigger the debug mode of the parallel geometry.\n
    possible values are:
    <ul>
    <li>0 : deactivated</li>
    <li>1 : activated</li>
    </ul>
    Keywords: <i>GRID, GENERATOR, PARALLEL, MASSIVE, TRIGGER, GEOMETRY</i>
  */
  m_debugParGeom = Context::getSolverProperty<MBool>("debugParGeom", m_solverId, AT_, &m_debugParGeom);

  // Store in if flow solver is enabled to differentiate the cases grid generator of flow solver for
  // parallel geometry
  m_flowSolver = false;
  m_flowSolver = Context::getBasicProperty<MBool>("flowSolver", AT_, &m_flowSolver);
}


/** \brief Returns the bounding box for the geometry
 *
 */
template <MInt nDim>
void Geometry<nDim>::getBoundingBox(MFloat* const bBox) const {
  for(MInt i = 0; i < nDim * 2; i++) {
    bBox[i] = m_minMax[i];
  }
}


template <MInt nDim>
void Geometry<nDim>::getBoundingBoxMB(MFloat* const bBox) const {
  for(MInt i = 0; i < nDim * 2; i++) {
    bBox[i] = m_mbminMax[i];
  }
}

template <MInt nDim>
void Geometry<nDim>::writeSegmentsToDX() {
  m_adt->writeTreeToDx();
}

/** \brief Compares two vectors entry by entry
 *
 * \author Andreas Lintermann
 * \date 18.09.2015
 *
 * \param[in] a vector a
 * \param[in] b vector b
 * return if both input vectors are equal
 *
 **/
template <MInt nDim>
MBool Geometry<nDim>::vectorsEqual(MFloat* a, MFloat* b) {
  MBool ret = true;
  for(MInt d = 0; d < nDim; d++)
    ret = ret && approx(a[d], b[d], MFloatEps);

  return ret;
}

/** \brief Returns the circumference of a segment
 *
 * \author Andreas Lintermann
 * \date 17.09.2015
 *
 * \param[in] bndVs holds the edges of the boundary
 * \param[in] num size of the input array
 * \return the circumference
 *
 **/
template <MInt nDim>
MFloat Geometry<nDim>::calcCircumference(MFloat** bndVs, MInt num) {
  MFloat circ = 0.0;
  MFloat tmp_length = 0.0;
  std::array<MFloat, nDim> edge;
  for(MInt i = 0; i < num; i++) {
    tmp_length = 0.0;
    for(MInt j = 0; j < nDim; j++) {
      edge[j] = bndVs[(i + 1) % num][j] - bndVs[i][j];
      tmp_length += edge[j] * edge[j];
    }
    circ += sqrt(tmp_length);
  }

  return circ;
}

/** \brief  Determines if a point is inside of the geometry
 *  \author (Miro Gondrum refactored; original by Lintermann/Schlottke)
 *  \date   20.10.2020
 *
 *  \param[in]  coordinates   the coordinates to test for inside/outside
 *  \param[out] numcutsperdir the number of cuts per direction
 *  \return if the point is inside the geometry (inside=true; outside=false)
 */
template <MInt nDim>
MBool Geometry<nDim>::pointIsInside(const MFloat* const coordinates) {
  std::array<MInt, nDim> numcutsperdir{};
  return pointIsInside(coordinates, numcutsperdir.data());
}

/** \brief  Determines if a point is inside of the geometry
 *  \author (Miro Gondrum refactored; original by Lintermann/Schlottke)
 *  \date   20.10.2020
 *
 *  \param[in]  coordinates   the coordinates to test for inside/outside
 *  \param[out] numcutsperdir the number of cuts per direction
 *  \return if the point is inside the geometry (inside=true; outside=false)
 */
template <MInt nDim>
MBool Geometry<nDim>::pointIsInside(const MFloat* const coordinates, MInt* numcutsperdir) {
  // TODO labels:GEOM,toenhance This function is not considering the distance from the next triangle.
  // If coordinates is located exactly on the surface the result
  // (inside/outside) is depended on its orientation.
  // By using distance to next element and one single ray the state
  // (inside/outside) should be defined better.
  TRACE();

  const MFloat* const box = boundingBox();
  MBool isInside = true;
  for(MInt dir = 0; dir < nDim; dir++) {
    // cast ray between point in question and a point out side the domain
    array<MFloat, 2 * nDim> ray;
    for(MInt i = 0; i < nDim; i++) {
      ray[i] = coordinates[i];
      ray[i + nDim] = coordinates[i]; // point in question
    }
    ray[dir] = box[dir] - (ray[dir] - box[dir]); // outside point

    // get list of intersecting geometry elements
    std::vector<MInt> nodeList;
    getLineIntersectionElements(&ray[0], nodeList);
    numcutsperdir[dir] = nodeList.size();
    isInside = isInside && (numcutsperdir[dir] % 2 == 0);
  }
  return isInside;
}


/** \brief Determines if a point is in or outside the geometry
 *
 * \author Andreas Lintermann
 * \date 21.09.2015
 *
 * \param[in] coordinates the coordinates to test for inside/outside
 * \param[in] the number of cuts per direction to be retunred
 * \return if the point is inside (false) or outside (true)
 *
 **/
// TODO labels:GEOM,DLB optimize this further to speedup reinitialization during DLB!
template <MInt nDim>
MBool Geometry<nDim>::pointIsInside2(const MFloat* const coordinates, MInt* numcutsperdir) {
  TRACE();
  MInt noNodes = 0;
  std::vector<MInt> nodeList;

  const MBool onlyCheck = (numcutsperdir == nullptr);

  MInt spaceDirectionId[3];
  MFloat target[6];
  const MFloat* tmp = boundingBox();

  // 1. Take an arbitrary unmarked cell
  // Cast rays in x,y and z direction and check for intersections
  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }

  target[0] = tmp[0] - (target[0] - tmp[0]); // x-ray to bounding box

  if(m_inOutTest == "perpOp") {
    getLineIntersectionElements(target, nodeList);
  } else {
    spaceDirectionId[0] = 1;
    spaceDirectionId[1] = 2;
    spaceDirectionId[2] = 0;
    getLineIntersectionElementsOld2(target, spaceDirectionId, nodeList);
  }
  const MInt noNodesX = nodeList.size();

  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }

  target[1] = tmp[1] - (target[1] - tmp[1]); // y-ray to bounding box
  nodeList.clear();
  if(m_inOutTest == "perpOp") {
    getLineIntersectionElements(target, nodeList);
  } else {
    spaceDirectionId[0] = 2;
    spaceDirectionId[1] = 0;
    spaceDirectionId[2] = 1;
    getLineIntersectionElementsOld2(target, spaceDirectionId, nodeList);
  }
  const MInt noNodesY = nodeList.size();

  // if modulo for x and y do not match return early if the number of cuts per direction is not needed
  if(onlyCheck) {
    if(noNodesX % 2 != noNodesY % 2) {
      return false;
    }
  }

  // only for 3D
  MInt noNodesZ = 0;
  IF_CONSTEXPR(nDim == 3) {
    for(MInt j = 0; j < nDim; j++) {
      target[j] = coordinates[j];
      target[j + nDim] = coordinates[j];
    }

    target[2] = tmp[2] - (target[2] - tmp[2]); // z-ray to bounding box
    // 2.b and c
    nodeList.clear();
    if(m_inOutTest == "perpOp") {
      getLineIntersectionElements(target, nodeList);
    } else {
      spaceDirectionId[0] = 0;
      spaceDirectionId[1] = 1;
      spaceDirectionId[2] = 2;
      getLineIntersectionElementsOld2(target, spaceDirectionId, nodeList);
    }
    noNodesZ = nodeList.size();

    if(noNodesX % 2 == noNodesY % 2 && noNodesX % 2 == noNodesZ % 2) {
      noNodes = noNodesX;
    } else {
      // cerr << " Different results for ray intersection, casting further rays. " << endl;
      // 	if(noNodesX)
      // 	  noNodes = noNodesX;
      // 	else if(noNodesY)
      // 	  noNodes = noNodesY;
      // 	else
      // 	  noNodes = noNodesZ;
      noNodes = 1;
    }
  }
  else {
    if(noNodesX % 2 == noNodesY % 2) {
      noNodes = noNodesX;
    } else {
      // cerr << " Different results for ray intersection, casting further rays. " << endl;
      noNodes = 1;
    }
  }

  if(!onlyCheck) {
    numcutsperdir[0] = noNodesX;
    numcutsperdir[1] = noNodesY;
    numcutsperdir[2] = noNodesZ;
  }

  if(noNodes % 2) {
    return false;
  } else {
    return true;
  }
}

// not very efficient - replace by some more sophisticated algorithm for multiple bodies
// 2D-version! does not allow for overlapping mb-Elements!
template <MInt nDim>
MBool Geometry<nDim>::pointIsInsideMBElements(const MFloat* const coordinates,
                                              MInt* bodyBndryCndIds,
                                              MInt* setToBodiesTable,
                                              MInt noBodiesInSet) {
  TRACE();

  MInt noNodes;
  MInt noNodesX_set = 0;
  MInt noNodesY_set = 0;
  MInt noNodesZ_set = 0;
  MInt body;
  MInt bcId;
  std::vector<MInt> nodeList;

  MFloat target[2 * nDim]{};
  MFloat tmp[2 * nDim]{};
  getBoundingBoxMB(&tmp[0]);

  // 1. Take an arbitrary unmarked cell
  // Cast rays in x,y and z direction and check for intersections
  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }
  target[0] = tmp[0] - (target[0] - tmp[0]); // x-ray to bounding box
  getLineIntersectionMBElements(target, nodeList);
  const MInt noNodesX = nodeList.size();
  for(MInt node = 0; node < noNodesX; node++) {
    for(MInt b = 0; b < noBodiesInSet; b++) {
      body = setToBodiesTable[b];
      bcId = bodyBndryCndIds[body];
      if(mbelements[nodeList[node]].m_bndCndId == bcId) {
        noNodesX_set++;
        break;
      }
    }
  }

  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }
  target[1] = tmp[1] - (target[1] - tmp[1]); // y-ray to bounding box
  nodeList.clear();
  getLineIntersectionMBElements(target, nodeList);
  const MInt noNodesY = nodeList.size();
  for(MInt node = 0; node < noNodesY; node++) {
    for(MInt b = 0; b < noBodiesInSet; b++) {
      body = setToBodiesTable[b];
      bcId = bodyBndryCndIds[body];
      if(mbelements[nodeList[node]].m_bndCndId == bcId) {
        noNodesY_set++;
        break;
      }
    }
  }

  // only for 3D
  IF_CONSTEXPR(nDim == 3) {
    for(MInt j = 0; j < nDim; j++) {
      target[j] = coordinates[j];
      target[j + nDim] = coordinates[j];
    }

    target[2] = tmp[2] - (target[2] - tmp[2]); // z-ray to bounding box
    // 2.b and c
    nodeList.clear();
    getLineIntersectionMBElements(target, nodeList);
    const MInt noNodesZ = nodeList.size();
    for(MInt node = 0; node < noNodesZ; node++) {
      for(MInt b = 0; b < noBodiesInSet; b++) {
        body = setToBodiesTable[b];
        bcId = bodyBndryCndIds[body];
        if(mbelements[nodeList[node]].m_bndCndId == bcId) {
          noNodesZ_set++;
          break;
        }
      }
    }

    if(noNodesX_set % 2 == noNodesY_set % 2 && noNodesX_set % 2 == noNodesZ_set % 2) {
      noNodes = noNodesX_set;
    } else {
      // cerr << " Different results for ray intersection, casting further rays. " << endl;
      //  if(noNodesX)
      //    noNodes = noNodesX;
      //  else if(noNodesY)
      //    noNodes = noNodesY;
      //  else
      //    noNodes = noNodesZ;
      noNodes = 1;
    }
  }
  else {
    if(noNodesX_set % 2 == noNodesY_set % 2) {
      noNodes = noNodesX_set;
    } else {
      // cerr << " Different results for ray intersection, casting further rays. " << endl;
      noNodes = 1;
    }
  }

  // delete [] tmp; tmp=nullptr;

  if(noNodes % 2) {
    // 3.
    return true;
  } else {
    return false;
  }
}

// not very efficient - replace by some more sophisticated algorithm for multiple bodies
// only intended for 3D, getLineIntersectionMBElements2 has no 2D implementation!
// allows for multiple overlapping geometries (however untested)
template <MInt nDim>
MBool Geometry<nDim>::pointIsInsideMBElements2(const MFloat* const coordinates,
                                               MInt* bodyBndryCndIds,
                                               MInt* setToBodiesTable,
                                               MInt noBodiesInSet) {
  TRACE();

  MInt noNodes;
  MInt noNodesX = 0;
  MInt noNodesY = 0;
  MInt noNodesZ = 0;
  MInt noNodesX_set = 0;
  MInt noNodesY_set = 0;
  MInt noNodesZ_set = 0;
  MInt body;
  MInt bcId;
  MBool debug = false;
  const MInt maxNoNodesDebug = 1000;
  MIntScratchSpace nodesXBackup(maxNoNodesDebug, AT_, "nodesXBackup");
  MIntScratchSpace nodesYBackup(maxNoNodesDebug, AT_, "nodesYBackup");
  MIntScratchSpace nodesZBackup(maxNoNodesDebug, AT_, "nodesZBackup");

  MInt spaceDirectionId[3];
  MFloat target[6];

  MFloat tmp[6];
  MFloat* tmpPtr = &tmp[0];
  getBoundingBoxMB(tmpPtr);

  // 1. Take an arbitrary unmarked cell
  // Cast rays in x,y and z direction and check for intersections
  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }
  spaceDirectionId[0] = 1;
  spaceDirectionId[1] = 2;
  spaceDirectionId[2] = 0;


  target[0] = tmp[0] - (target[0] - tmp[0]); // x-ray to bounding box
  MInt count = 0;
  for(MInt b = 0; b < noBodiesInSet; b++) {
    body = setToBodiesTable[b];
    bcId = bodyBndryCndIds[body];
    std::vector<MInt> nodeList;
    getLineIntersectionMBElements2(target, spaceDirectionId, nodeList, bcId);
    noNodesX = nodeList.size();
    noNodesX_set += noNodesX;
    if(debug) {
      if(noNodesX > maxNoNodesDebug)
        mTerm(1, AT_, " too many intersecting nodes to debug... please check or increase maxNoNodesDebug! ");
      for(MInt n = 0; n < noNodesX; n++) {
        nodesXBackup[count] = nodeList[n];
        count++;
      }
    }
  }

  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }
  spaceDirectionId[0] = 2;
  spaceDirectionId[1] = 0;
  spaceDirectionId[2] = 1;

  target[1] = tmp[1] - (target[1] - tmp[1]); // y-ray to bounding box
  count = 0;
  for(MInt b = 0; b < noBodiesInSet; b++) {
    body = setToBodiesTable[b];
    bcId = bodyBndryCndIds[body];
    std::vector<MInt> nodeList;
    getLineIntersectionMBElements2(target, spaceDirectionId, nodeList, bcId);
    noNodesY = nodeList.size();
    noNodesY_set += noNodesY;
    if(debug) {
      if(noNodesY > maxNoNodesDebug)
        mTerm(1, AT_, " too many intersecting nodes to debug... please check or increase maxNoNodesDebug! ");
      for(MInt n = 0; n < noNodesY; n++) {
        nodesYBackup[count] = nodeList[n];
        count++;
      }
    }
  }

  // only for 3D
  IF_CONSTEXPR(nDim == 3) {
    for(MInt j = 0; j < nDim; j++) {
      target[j] = coordinates[j];
      target[j + nDim] = coordinates[j];
    }
    spaceDirectionId[0] = 0;
    spaceDirectionId[1] = 1;
    spaceDirectionId[2] = 2;


    target[2] = tmp[2] - (target[2] - tmp[2]); // z-ray to bounding box
    count = 0;
    // 2.b and c
    for(MInt b = 0; b < noBodiesInSet; b++) {
      body = setToBodiesTable[b];
      bcId = bodyBndryCndIds[body];
      std::vector<MInt> nodeList;
      getLineIntersectionMBElements2(target, spaceDirectionId, nodeList, bcId);
      noNodesZ = nodeList.size();
      noNodesZ_set += noNodesZ;
      if(debug) {
        if(noNodesZ > maxNoNodesDebug)
          mTerm(1, AT_, " too many intersecting nodes to debug... please check or increase maxNoNodesDebug! ");
        for(MInt n = 0; n < noNodesZ; n++) {
          nodesZBackup[count] = nodeList[n];
          count++;
        }
      }
    }

    if(noNodesX_set % 2 == noNodesY_set % 2 && noNodesX_set % 2 == noNodesZ_set % 2) {
      noNodes = noNodesX_set;
    } else {
      cerr << " Different results for ray intersection, continue with most probable result. Please check your "
              "geometry! "
           << endl;
      cerr << " nodesX: " << noNodesX << " " << noNodesX_set << " nodesY: " << noNodesY << " " << noNodesY_set
           << " nodesZ: " << noNodesZ << " " << noNodesZ_set << endl;
      cerr << " write stl of all intersecting nodes, point is " << coordinates[0] << " " << coordinates[1];
      IF_CONSTEXPR(nDim == 3) cerr << " " << coordinates[2];
      cerr << endl;
      cerr << " bounding box is " << tmp[0] << " " << tmp[1] << " " << tmp[2] << " " << tmp[3] << " " << tmp[4] << " "
           << tmp[5] << endl;
      for(MInt j = 0; j < nDim; j++) {
        target[j] = coordinates[j];
        target[j + nDim] = coordinates[j];
      }

      ofstream ofl;
      stringstream filename;
      filename << "IntersectingNodes_D" << domainId() << ".stl";
      ofl.open((filename.str()).c_str());
      ofl.precision(12);

      if(ofl) {
        if(noNodesX_set) {
          ofl << "solid intersecting triangles X " << domainId() << endl;

          for(MInt i = 0; i < noNodesX_set; i++) {
            ofl << "  facet normal";
            for(MInt j = 0; j < 3; j++)
              ofl << " " << mbelements[nodesXBackup[i]].m_normal[j];
            ofl << endl << "  outer loop" << endl;

            for(MInt k = 0; k < 3; k++) {
              ofl << "   vertex";
              for(MInt l = 0; l < 3; l++) {
                ofl << " " << mbelements[nodesXBackup[i]].m_vertices[k][l];
              }
              ofl << endl;
            }

            ofl << "  endloop" << endl << " endfacet" << endl;
          }
          ofl << "endsolid intersecting triangles X " << domainId() << endl;
        }
        if(noNodesY_set) {
          ofl << "solid intersecting triangles Y " << domainId() << endl;

          for(MInt i = 0; i < noNodesY_set; i++) {
            ofl << "  facet normal";
            for(MInt j = 0; j < 3; j++)
              ofl << " " << mbelements[nodesYBackup[i]].m_normal[j];
            ofl << endl << "  outer loop" << endl;

            for(MInt k = 0; k < 3; k++) {
              ofl << "   vertex";
              for(MInt l = 0; l < 3; l++) {
                ofl << " " << mbelements[nodesYBackup[i]].m_vertices[k][l];
              }
              ofl << endl;
            }

            ofl << "  endloop" << endl << " endfacet" << endl;
          }
          ofl << "endsolid intersecting triangles Y " << domainId() << endl;
        }
        if(noNodesZ_set) {
          ofl << "solid intersecting triangles Z " << domainId() << endl;

          for(MInt i = 0; i < noNodesZ_set; i++) {
            ofl << "  facet normal";
            for(MInt j = 0; j < 3; j++)
              ofl << " " << mbelements[nodesZBackup[i]].m_normal[j];
            ofl << endl << "  outer loop" << endl;

            for(MInt k = 0; k < 3; k++) {
              ofl << "   vertex";
              for(MInt l = 0; l < 3; l++) {
                ofl << " " << mbelements[nodesZBackup[i]].m_vertices[k][l];
              }
              ofl << endl;
            }

            ofl << "  endloop" << endl << " endfacet" << endl;
          }
          ofl << "endsolid intersecting triangles Z " << domainId() << endl;
        }
      }


      if(noNodesX % 2 == noNodesY % 2)
        noNodes = noNodesX;
      else if(noNodesX % 2 == noNodesZ % 2)
        noNodes = noNodesX;
      else if(noNodesY % 2 == noNodesZ % 2)
        noNodes = noNodesY;
      else
        noNodes = 1;
    }
  }
  else {
    if(noNodesX_set % 2 == noNodesY_set % 2) {
      noNodes = noNodesX_set;
    } else {
      cerr << " Different results for ray intersection, continue with most probable result. Please check your "
              "geometry! "
           << endl;
      noNodes = 1;
    }
  }

  // delete [] tmp; tmp=nullptr;

  if(noNodes % 2) {
    return true;
  } else {
    return false;
  }
}


/** \brief returns the geometry elements that have a cut with rays originating in the provided coordinates
 *
 * \author Andreas Lintermann
 * \date 04.05.2016
 *
 * The algorithm does the following
 *   1. create the orginal target from the coordinates
 *   2. do the check for all directions
 *   2.1 create a ray from current target to bounding box in direction dim
 *   2.2 get the intersecting elements
 *   2.3 create the resulting nodes vector
 *
 * \param[in] coordinates the original coordinates to cast rays from
 * \param[in] resultnodes contains the intersected elements per direction
 **/
template <MInt nDim>
void Geometry<nDim>::determineRayIntersectedElements(const MFloat* const coordinates,
                                                     vector<vector<MInt>>* resultnodes) {
  TRACE();

  std::vector<MInt> nodeList;
  const MFloat* bb = boundingBox();

  // 1. create the orginal target from the coordinates
  std::array<MFloat, nDim * 2> target;
  for(MInt j = 0; j < nDim; j++) {
    target[j] = coordinates[j];
    target[j + nDim] = coordinates[j];
  }

  // 2. do the check for all directions
  for(MInt dim = 0; dim < nDim; dim++) {
    // 2.1 create a ray from current target to bounding box in direction dim
    std::array<MFloat, nDim * 2> raybb;
    for(MInt j = 0; j < 2 * nDim; j++)
      raybb[j] = target[j];
    raybb[dim] = bb[dim] - (target[dim] - bb[dim]);

    // 2.2 get the intersecting elements
    getLineIntersectionElements(raybb.data(), nodeList);

    // 2.3 create the resulting nodes vector
    vector<MInt> v;
    for(MInt i = 0; i < (signed)nodeList.size(); i++)
      v.push_back(nodeList[i]);

    resultnodes->push_back(v);
  }
}


/// \brief Return the set of boundary condition ids of the elements cut by the given line
template <MInt nDim>
void Geometry<nDim>::getLineIntersectingElementsBcIds(const MFloat* const line, std::set<MInt>& bcIds) {
  std::vector<MInt> nodeList;
  MFloat tmp_line[2 * nDim];
  std::copy_n(line, 2 * nDim, tmp_line);

  getLineIntersectionElements(tmp_line, nodeList);

  bcIds.clear();
  for(MInt i = 0; i < (signed)nodeList.size(); i++) {
    const MInt bndCndId = elements[nodeList[i]].m_bndCndId;
    bcIds.insert(bndCndId);
  }
}

template <MInt nDim>
MBool Geometry<nDim>::isOnGeometry(const MFloat cellHalfLength, const MFloat* coordinates, MString gridCutTest) {
  MFloat target[2 * nDim] = {0.0};
  std::vector<MInt> nodeList;

  for(MInt dim = 0; dim < nDim; dim++) {
    target[dim] = coordinates[dim] - cellHalfLength;
    target[dim + nDim] = coordinates[dim] + cellHalfLength;
  }

  if(gridCutTest == "SAT") {
    this->getIntersectionElements(target, nodeList, cellHalfLength, coordinates);
  } else {
    this->getIntersectionElements(target, nodeList);
  }

  return nodeList.size() > 0;
}

// Explicit instantiations for 2D and 3D
template class Geometry<2>;
template class Geometry<3>;
