// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FCDESCRIPTOR_H
#define FCDESCRIPTOR_H

#include "INCLUDE/maiaconstants.h"
#include "UTIL/functions.h"
#include "UTIL/maiamath.h"

// Description: This file contains arrays depending on number of space
// dimensions D. All these template arrays are wrapped in
// the struct FcDescriptor.
// Usage: Defining a template specific typedef once in a class all arrays
// can be easily accesed, p.e.:
//    using Fd = fcDescriptor<nDim>;

//--declarations of lattice descriptor------------------------------------------
namespace fcDescriptor {

// Holds the lobatto points for the first 4 poly degrees
// The array is filled with -1
inline constexpr MFloat lobattoPoints[4][3] = {
    {-F1, -F1, -F1}, {F0, -F1, -F1}, {-F1 / SQRT5, F1 / SQRT5, -F1}, {-SQRT3 / SQRT7, F0, SQRT3 / SQRT7}};

// Holds the position of the neighboring cells in the neighbor list that share
// vertex number v of the current cell with the current cell
template <MInt D>
inline constexpr MInt vertexNghbrs[IPOW2(D)][IPOW2(D) - 1] = {};

// Holds the position of the neighboring cells in the neighbor list that share
// edge number e of the current cell with the current cell
inline constexpr MInt edgeNghbrs3d[12][3] = {{0, 2, 6},  {0, 4, 10}, {0, 5, 11}, {0, 3, 7},  {2, 4, 14}, {2, 5, 15},
                                             {3, 4, 16}, {3, 5, 17}, {1, 2, 8},  {1, 4, 12}, {1, 5, 13}, {1, 3, 9}};

// Holds the position of the neighboring cell in the neighbor list that shares
// surface number s of the current cell with the current cell
template <MInt D>
inline constexpr MInt surfaceNghbrs[D * 2] = {};

// Holds the directions in which the vertices are located for a cell
// The order is based on the order of the legendre points in the cell
template <MInt D>
inline constexpr MFloat vertexPos[IPOW2(D)][D] = {};

// Holds the directions in which the edges are located for a 3D cell
// The order is based on the order of the legendre points in the cell
// Edges of surface 0: (-x,-y) -> (-x,-z) -> (-x,+z) -> (-x,-y)
// Edges of surface 1: (-y,-z) -> (-x,+z)
// Edges of surface 4: (+y,-z) -> (+y,+z)
// Edges of surface 5: (+x,-y) -> (+x,-z) -> (+x,+z) -> (+x,+y)
inline constexpr MFloat edgePos3d[12][3] = {{-F1, -F1, F0}, {-F1, F0, -F1}, {-F1, F0, F1}, {-F1, F1, F0},
                                            {F0, -F1, -F1}, {F0, -F1, F1},  {F0, F1, -F1}, {F0, F1, F1},
                                            {F1, -F1, F0},  {F1, F0, -F1},  {F1, F0, F1},  {F1, F1, F0}};

template <MInt D>
inline constexpr MFloat surfacePos[D * 2][D] = {};

template <MInt D>
inline constexpr MInt oppositeVertex[IPOW2(D)][IPOW2(D) - 1] = {};

// Holds the opposite edges for each edge of a 3D cell for each neighboring cell.
// Works as the oppositeVertex3d array
inline constexpr MInt oppositeEdge3d[12][3] = {{8, 3, 11}, {9, 2, 10}, {10, 1, 9}, {11, 0, 8}, {6, 5, 7}, {7, 4, 6},
                                               {4, 7, 5},  {5, 6, 4},  {0, 11, 3}, {1, 10, 2}, {2, 9, 1}, {3, 8, 0}};

template <MInt D>
inline constexpr MInt oppositeSurface[2 * D] = {};
} // namespace fcDescriptor

//--definitions of descriptor for specific d ------------------------
namespace fcDescriptor {

// Holds the position of the neighboring cells in the neighbor list that share
// vertex number v of the current cell with the current cell
template <>
inline constexpr MInt vertexNghbrs<2>[4][3] = {{0, 2, 6}, {0, 3, 7}, {1, 2, 5}, {1, 3, 4}};

// Holds the position of the neighboring cells in the neighbor list that share
// vertex number v of the current cell with the current cell
template <>
inline constexpr MInt vertexNghbrs<3>[8][7] = {
    {0, 2, 4, 6, 10, 14, 18}, {0, 2, 5, 6, 11, 15, 19}, {0, 3, 4, 7, 10, 16, 20}, {0, 3, 5, 7, 11, 17, 21},
    {1, 2, 4, 8, 12, 14, 22}, {1, 2, 5, 8, 13, 15, 23}, {1, 3, 4, 9, 12, 16, 24}, {1, 3, 5, 9, 13, 17, 25}};

// Holds the position of the neighboring cell in the neighbor list that shares
// surface number s of the current cell with the current cell
template <>
inline constexpr MInt surfaceNghbrs<2>[4] = {0, 2, 3, 1};

// Holds the position of the neighboring cell in the neighbor list that shares
// surface number s of the current cell with the current cell
template <>
inline constexpr MInt surfaceNghbrs<3>[6] = {0, 2, 4, 5, 3, 1};

// Holds the directions in which the vertices are located for a 2D cell
// The order is based on the order of the legendre points in the cell
template <>
inline constexpr MFloat vertexPos<2>[4][2] = {{-F1, -F1}, {-F1, F1}, {F1, -F1}, {F1, F1}};

template <>
// Holds the directions in which the vertices are located for a 3D cell
// The order is based on the order of the legendre points in the cell
inline constexpr MFloat vertexPos<3>[8][3] = {{-F1, -F1, -F1}, {-F1, -F1, F1}, {-F1, F1, -F1}, {-F1, F1, F1},
                                              {F1, -F1, -F1},  {F1, -F1, F1},  {F1, F1, -F1},  {F1, F1, F1}};

// Holds the directions in which the edges are located for a 2D cell
//     __2__
//    |     |
//   0|     |3
//    |__ __|
//       1
template <>
inline constexpr MFloat surfacePos<2>[4][2] = {{-F1, F0}, {F0, -F1}, {F0, F1}, {F1, F0}};

// Holds the direction in which the surfaces are located for a 3D cell
// The order is based on the order of the legendre points in the cell
// -x -> -y -> -z -> +z -> +y -> +x
template <>
inline constexpr MFloat surfacePos<3>[6][3] = {{-F1, F0, F0}, {F0, -F1, F0}, {F0, F0, -F1},
                                               {F0, F0, F1},  {F0, F1, F0},  {F1, F0, F0}};

// Holds the opposite vertices for each vertex of a 3D cell for each neighboring cell.
// Based on vertexPos3d and vertexNghbrs3d
// For example, vertex 0 of cell a is shared with 7 other cells. Their location is
// given in the arry vertexNghbrs3d. In the neighbors vertex order the corresponding
// vertex has the number as given in this array
//->For vertex 0 in cell A the 3. neighbor is located at position 4 of the neighbor list
// the corresponding vertex position is 1
template <>
inline constexpr MInt oppositeVertex<3>[8][7] = {{4, 2, 1, 6, 5, 3, 7}, {5, 3, 0, 7, 4, 2, 6}, {6, 0, 3, 4, 7, 1, 5},
                                                 {7, 1, 2, 5, 6, 0, 4}, {0, 6, 5, 2, 1, 7, 3}, {1, 7, 4, 3, 0, 6, 2},
                                                 {2, 4, 7, 0, 3, 5, 1}, {3, 5, 6, 1, 2, 4, 0}};

// Holds the opposite vertices for each vertex of a 2D cell for each neighboring cell.
// Based on vertexPos2d and vertexNghbrs2d
// For example, vertex 0 of cell a is shared with 3 other cells. Their location is
// given in the arry vertexNghbrs2d. In the neighbors vertex order the corresponding
// vertex has the number as given in this array
//->For vertex 0 in cell A the 3. neighbor is located at position 6 of the neighbor list
// the corresponding vertex position is 1
template <>
inline constexpr MInt oppositeVertex<2>[4][3] = {{2, 1, 3}, {3, 0, 2}, {0, 3, 1}, {1, 2, 0}};

// Holds the opposite edge for each edge of a 2D cell for each neighboring cell.
// Works as the oppositeVertex2d array
template <>
inline constexpr MInt oppositeSurface<2>[4] = {3, 2, 1, 0};

// Holds the opposite surface for each edge of a 3D cell for each neighboring cell.
// Works as the oppositeVertex3d array
template <>
inline constexpr MInt oppositeSurface<3>[6] = {5, 4, 3, 2, 1, 0};

} // namespace fcDescriptor

//--wrapper struct for lattice descriptor---------------------------------------
/** \brief  LB lattice descriptor for arrays depending on D
 *  \author Miro Gondrum
 *  \date   27.05.2021
 *  Wrapper struct to easily access LB arrays depending on number of space
 *  dimension D.
 **/
template <MInt D>
struct FcDescriptor {
  FcDescriptor() = delete;

  static constexpr MFloat legendreFunction(const MInt order, const MFloat pos) {
    MFloat p = pos;
    MFloat p_deriv = F1;
    maia::math::calculateLegendrePolyAndDeriv(order, pos, &p, &p_deriv);
    return p;
  }

  static constexpr MFloat legendreDerivFunction(const MInt order, const MFloat pos) {
    MFloat p = pos;
    MFloat p_deriv = F1;
    maia::math::calculateLegendrePolyAndDeriv(order, pos, &p, &p_deriv);
    return p_deriv;
  }

  static constexpr void gaussPoint(const MInt order, MFloat* nodes, MFloat* weights) {
    maia::math::calculateLegendreGaussNodesAndWeights(order, nodes, weights);
  }

  static constexpr void nodePosition(const MInt node, const MInt nodesPerDir, MInt* nodePos) {
    MInt remainNoNodes = node;
    for(MInt d = D - 1; d >= 0; d--) {
      nodePos[d] = remainNoNodes % nodesPerDir;
      remainNoNodes = remainNoNodes / nodesPerDir;
    }
  }

  static constexpr MFloat nodePosEquidist(const MInt order, const MInt node) {
    const MFloat p = (MFloat)order;
    const MFloat base = F2 / (p + F1);
    const MFloat exp = (MFloat)node;

    MFloat nodePos = -F1 + exp * base;
    return nodePos;
  }

  static constexpr MFloat nodePosLobattoPoints(MInt i, MInt j) { return fcDescriptor::lobattoPoints[i][j]; }

  static constexpr MInt nghbrCellOfVertex(MInt i, MInt j) { return fcDescriptor::vertexNghbrs<D>[i][j]; }

  static constexpr MInt nghbrCellOfEdge(MInt i, MInt j) { return fcDescriptor::edgeNghbrs3d[i][j]; }

  static constexpr MInt nghbrCellOfSurface(MInt i) { return fcDescriptor::surfaceNghbrs<D>[i]; }

  static constexpr MFloat vertexPosition(MInt i, MInt j) { return fcDescriptor::vertexPos<D>[i][j]; }

  static constexpr MFloat edgePosition(MInt i, MInt j) { return fcDescriptor::edgePos3d[i][j]; }

  static constexpr MFloat surfacePosition(MInt i, MInt j) { return fcDescriptor::surfacePos<D>[i][j]; }

  static constexpr MInt vertexIdOfOppCell(MInt i, MInt j) { return fcDescriptor::oppositeVertex<D>[i][j]; }

  static constexpr MInt edgeIdOfOppCell(MInt i, MInt j) { return fcDescriptor::oppositeEdge3d[i][j]; }

  static constexpr MInt surfaceIdOfOppCell(MInt i) { return fcDescriptor::oppositeSurface<D>[i]; }
};

#endif
