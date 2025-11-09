// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYADT_H
#define GEOMETRYADT_H

#include <array>
#include <stack>
#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "config.h"
#include "geometryelement.h"

template <MInt nDim>
class GeometryElement;

/** This class implements an alternating digital tree
 *
 *  The ADT stores geometrical elements in a hierarchical 2-d dimensional
 *  order, allowing fast proximity queries.
 *  This implementation follows the structure outlined by Aftosmis:
 *
 *  von Karman Institute for Fluid Dynamics
 *  Lecture Series 1997-02
 *  28th computational fluid dynamics
 *  Solution Adapative Cartesian Grid Methods for Aerodynamic Flows
 *  with Complex Geometries
 */

/** This namespace holds the function objects for sorting
 *  classes.
 *
 */
namespace sortFunctions {
// Function objects for sorting the elements
template <MInt nDim>
struct lessMinMax {
  using element = GeometryElement<nDim>;
  lessMinMax(const MInt dir, const element* const elements) : m_dir(dir), m_elements(elements) {}

  MBool operator()(const MInt a, const MInt b) { return m_elements[a].m_minMax[m_dir] < m_elements[b].m_minMax[m_dir]; }

 private:
  const MInt m_dir = -1;
  const element* const m_elements = nullptr;
};
} // namespace sortFunctions

template <MInt nDim>
class Geometry;
class GeometryAdtNode;

template <MInt nDim>
class GeometryAdt {
 public:
  GeometryAdtNode* m_nodes;
  std::array<MFloat, 2 * nDim> m_minMax{};
  GeometryAdtNode* m_mbnodes;
  std::array<MFloat, 2 * nDim> m_mbminMax{};

 protected:
  MInt m_root;
  const Geometry<nDim>* m_geometry;

 public:
  explicit GeometryAdt(const Geometry<nDim>* geometry);

  ~GeometryAdt();
  void buildTree();
  void buildTreeOld();
  void buildTreeNew();
  void writeTreeToDx();
  void retrieveNodes(MFloat* targetRegion, std::vector<MInt>& nodeList) const;

  void splitTree(MInt noSubTrees);

  void buildTreeMB();
  void retrieveNodesMBElements(const MFloat* targetRegion, std::vector<MInt>& nodeList) const;

  MInt get_root();

  MLong memoryUsage();

 private:
  void writeNode(MInt node);
};

#endif
