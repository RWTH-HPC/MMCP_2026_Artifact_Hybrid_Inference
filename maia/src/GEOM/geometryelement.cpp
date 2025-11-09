// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryelement.h"

#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/maiamath.h"
#include "UTIL/timer.h"

using namespace std;

//----------------------------------------------------------
// Definitions of class element
//----------------------------------------------------------

template <MInt nDim>
using element = GeometryElement<nDim>;

template <MInt nDim>
void GeometryElement<nDim>::allocateElements(void* tmpPointer, void* /*dummy*/, MInt& /*dummy2*/) {
  //  TRACE();

  m_vertices = (MFloat**)tmpPointer;
  // Jump memory for pointers to pointers
  for(MInt i = 0; i < nDim; i++) {
    tmpPointer = (void*)((MFloat**)tmpPointer + 1);
  }

  for(MInt i = 0; i < nDim; i++) {
    // Set pointers to pointers
    m_vertices[i] = (MFloat*)tmpPointer;
    // Jump memory for values
    tmpPointer = (void*)((MFloat*)tmpPointer + nDim);
  }
  m_normal = (MFloat*)tmpPointer;
  tmpPointer = (void*)((MFloat*)tmpPointer + nDim);
  m_minMax = (MFloat*)tmpPointer;
}

template <MInt nDim>
void GeometryElement<nDim>::boundingBox() {
  //  TRACE();

  for(MInt j = 0; j < nDim; j++) {
    m_minMax[j] = m_vertices[0][j];
    m_minMax[j + nDim] = m_vertices[0][j];
  }
  for(MInt i = 0; i < nDim; i++) {
    for(MInt j = 0; j < nDim; j++) {
      // Find maximum
      m_minMax[j + nDim] = (m_minMax[j + nDim] < m_vertices[i][j]) ? m_vertices[i][j] : m_minMax[j + nDim];
      // Find minimum
      m_minMax[j] = (m_minMax[j] > m_vertices[i][j]) ? m_vertices[i][j] : m_minMax[j];
    }
  }
}

template <MInt nDim>
void GeometryElement<nDim>::writeElement() const {
  TRACE();

  cerr.precision(10);
  cerr << "------------------------------" << endl;
  // M_Normal vector only in 3D (for now)
  IF_CONSTEXPR(nDim == 3) { cerr << "M_Normal : " << m_normal[0] << " " << m_normal[1] << " " << m_normal[2] << endl; }

  for(MInt i = 0; i < nDim; i++) {
    cerr << "Vertex " << i << ": ";
    for(MInt j = 0; j < nDim; j++) {
      cerr << m_vertices[i][j] << " ";
    }
    cerr << endl;
  }
  cerr << "Bounding box: " << endl << "m_min: ";

  for(MInt i = 0; i < nDim; i++) {
    cerr << m_minMax[i] << " ";
  }
  cerr << endl;
  cerr << "Bounding box: " << endl << "m_max: ";

  for(MInt i = 0; i < nDim; i++) {
    cerr << m_minMax[i + nDim] << " ";
  }
  cerr << endl;
}


/// \brief Calculates the normal vector from the  geometry element  vertices
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Rodrigo Miguez
/// \date 2014-03-18, 2019-10-31
///
/// \param[in]  vertices Geometry element vertices
/// \param[out] normal After a call to calcNormal(), `normal` holds the
///                    normalized normal vector. Must point to storage of at
///                    least size nDim.
template <MInt nDim>
void GeometryElement<nDim>::calcNormal(const MFloat* const vertices, MFloat* normal) const {
  TRACE();

  using Vec2D = array<MFloat, 2>;
  using Vec3D = array<MFloat, 3>;

  // Since the used algorithms depend on the dimension, we need to branch here
  IF_CONSTEXPR(nDim == 2) {
    Vec2D edge, norm;

    // Calculate edge vector
    edge[0] = vertices[2] - vertices[0];
    edge[1] = vertices[3] - vertices[1];

    // Calculate vector normal to edge
    norm[0] = edge[1];
    norm[1] = -edge[0];

    // Normalize normal vector and copy to result
    maia::math::normalize(norm);
    copy(norm.begin(), norm.end(), normal);
  }
  else IF_CONSTEXPR(nDim == 3) {
    Vec3D edge0, edge1, norm;

    // Calculate edge vectors
    edge0[0] = vertices[3] - vertices[0];
    edge0[1] = vertices[4] - vertices[1];
    edge0[2] = vertices[5] - vertices[2];
    edge1[0] = vertices[6] - vertices[0];
    edge1[1] = vertices[7] - vertices[1];
    edge1[2] = vertices[8] - vertices[2];

    // Calculate vector normal to plane
    norm = maia::math::cross(edge0, edge1);

    // Normalize normal vector and copy to result
    maia::math::normalize(norm);
    copy(norm.begin(), norm.end(), normal);
  }
  else {
    TERMM(1, "Bad number of dimensions.");
  }
}


/// \brief Calculate the centroid of the geometry element  vertices
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>, Rodrigo Miguez
/// \date 2014-03-19, 2019-10-31
///
/// \param[in]  vertices Geometry element vertices
/// \param[out] centroid After a call to calcCentroid(), `centroid` holds the
///                      centroid of the element. Must point to storage of at
///                      least size nDim.
template <MInt nDim>
void GeometryElement<nDim>::calcCentroid(const MFloat* const vertices, MFloat* centroid) const {
  TRACE();

  // Calculate sums over vertices and divide by dimensions (works for 2D and 3D)
  // Example 2D:
  // vertices: (x0, y0), (x1, y1)
  // centroid: ((x0 + x1)/2, (y0 + y1)/2)
  for(MInt dim = 0; dim < nDim; dim++) {
    centroid[dim] = 0.0;
    for(MInt vertexId = 0; vertexId < nDim; vertexId++) {
      centroid[dim] += vertices[(vertexId * nDim) + dim];
    }
    centroid[dim] /= static_cast<MFloat>(nDim);
  }
}


/// \brief Return the vertices of the geometry element
///
/// \author Rodrigo Barros Miguez (rodrigo) <rodrigo.miguez@rwth-aachen.de>
/// \date June 2019
///
/// \param[out] vertices After a call to getVertices(), 'vertices' holds the
///                      vertices of the geometry element.
template <MInt nDim>
void GeometryElement<nDim>::getVertices(MFloat* vertices) const {
  TRACE();

  IF_CONSTEXPR(nDim == 2) {
    // Copy vertices information to vectors
    vertices[0] = m_vertices[0][0];
    vertices[1] = m_vertices[0][1];
    vertices[2] = m_vertices[1][0];
    vertices[3] = m_vertices[1][1];
  }
  else IF_CONSTEXPR(nDim == 3) {
    // Copy vertices information to vectors
    vertices[0] = m_vertices[0][0];
    vertices[1] = m_vertices[0][1];
    vertices[2] = m_vertices[0][2];
    vertices[3] = m_vertices[1][0];
    vertices[4] = m_vertices[1][1];
    vertices[5] = m_vertices[1][2];
    vertices[6] = m_vertices[2][0];
    vertices[7] = m_vertices[2][1];
    vertices[8] = m_vertices[2][2];
  }
}


// Explicit instantiations for 2D and 3D
template class GeometryElement<2>;
template class GeometryElement<3>;
