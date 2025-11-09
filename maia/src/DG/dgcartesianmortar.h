// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGMORTAR_H_
#define DGMORTAR_H_

#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"
#include "dgcartesianinterpolation.h"

#include <string>

namespace maia {
namespace dg {
namespace mortar {


const MBool forward = true;
const MBool reverse = false;
const MInt lower = 0;
const MInt upper = 1;


/// Calculate p-refinement mortar projection matrix from source to target
/// polynomial degree.
///
/// \author Sven Berger, Michael Schlottke
/// \date 2015-02-27
///
/// \param[in] sourcePolyDeg Source polynomial degree.
/// \param[in] sourceNodes Node locations in [-1,1] of source polynomial degree.
/// \param[in] sourceBaryWeights Barycentric weights of source polynomial
///                              degree.
/// \param[in] targetPolyDeg Target polynomial degree.
/// \param[in] targetNodes Node locations in [-1,1] of target polynomial degree.
/// \param[out] matrix Projection matrix of size
///                    (targetPolyDeg + 1) x (sourcePolyDeg + 1).
///
/// This method can be used to calculate the projection matrix between any two
/// polynomial degrees. For each forward projection, the corresponding reverse
/// projection can be found by reversing the source and target polynomial
/// degrees.
///
/// When projecting from polynomial degree N_s to N_t, then sourcePolyDeg = N_s
/// and targetPolyDeg = N_t. The reverse projection can then be obtained by
/// setting sourcePolyDeg = N_t and targetPolyDeg = N_s.
///
/// The content of the projection matrix P_ij is from polynomial degree N_s to
/// N_t is as follows:
///
///    P_ij = l_j(x_i)
///
/// where l_j are the N_s + 1 Lagrange polynomials of degree N_s while x_i are
/// the N_t + 1 Gauss nodes of degree N_t.
///
/// See also p. 22, Eq. 3.33 and Eq. 3.34 in
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
inline void calcMortarProjectionMatrixP(const MInt sourcePolyDeg,
                                        const MFloat* const sourceNodes,
                                        const MFloat* const sourceBaryWeights,
                                        const MInt targetPolyDeg,
                                        const MFloat* const targetNodes,
                                        MFloat* const matrix) {
  const MInt sourceNoNodes = sourcePolyDeg + 1;
  const MInt targetNoNodes = targetPolyDeg + 1;

  // Matrix is Lagrange polynomials of source polynomial degree evaluated at
  // nodes of target polynomial degree
  for(MInt i = 0; i < targetNoNodes; i++) {
    interpolation::calcLagrangeInterpolatingPolynomials(
        targetNodes[i], sourcePolyDeg, sourceNodes, sourceBaryWeights, &matrix[i * sourceNoNodes]);
  }
}


/// Calculate h-refinement forward projection matrix from a coarse element to
/// the mortar element (2D).
///
/// \author Sven Berger, Michael Schlottke
/// \date 2015-02-27
///
/// \param[in] polyDeg Polynomial degree of element/mortar.
/// \param[in] nodes Node locations of Lagrange polynomial in [-1,1].
/// \param[in] wBary Barycentric weights of Lagrange polynomial.
/// \param[in] position Determines lower or upper element.
/// \param[out] matrix Projection matrix of size (polyDeg + 1) x (polyDeg + 1).
///
/// Projection:
/// lower:       upper:
///       |         _  |
/// |     |     |   /| |
/// |           |  /
/// | \         |
/// | _\| |     |      |
///       |            |
///
/// See also pp. 25-26 in
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
inline void calcMortarProjectionMatrixHForward(const MInt polyDeg,
                                               const MFloat* const nodes,
                                               const MFloat* const wBary,
                                               const MInt position,
                                               MFloat* const matrix) {
  const MInt noNodes = polyDeg + 1;

  // Forward matrix is the Lagrange polynomials evaluated at other nodes
  const MFloat shift = (position == dg::mortar::lower) ? -0.5 : 0.5;
  for(MInt i = 0; i < noNodes; i++) {
    maia::dg::interpolation::calcLagrangeInterpolatingPolynomials(
        0.5 * nodes[i] + shift, polyDeg, nodes, wBary, &matrix[i * noNodes]);
  }
}


/// Calculate h-refinement reverse projection matrix from a mortar element to
/// the coarse element (2D).
///
/// \author Sven Berger, Michael Schlottke
/// \date 2015-02-27
///
/// \param[in] polyDeg Polynomial degree of element/mortar.
/// \param[in] nodes Node locations of Lagrange polynomial in [-1,1].
/// \param[in] wBary Barycentric weights of Lagrange polynomial.
/// \param[in] position Determines lower or upper element.
/// \param[out] matrix Projection matrix of size (polyDeg + 1) x (polyDeg + 1).
///
/// Projection:
/// lower:       upper:
///       |           |
/// |     |     |   / |
/// |  _        | |/_
/// | |\        |
/// |   \ |     |     |
///       |           |
///
/// See also pp. 25-26 in
///     Sven Berger: Implementation and validation of an adaptive hp-refinement
///     method for the discontinuous Galerkin spectral element method. Master
///     thesis, RWTH Aachen University, 2014.
/// and pp. 11-14 in
///     Patrick Antony: Development of a coupled discontinuous Galerkin
///     -finite volume scheme, Bachelor thesis, RWTH Aachen University, 2018.
///
///  Contrary to the method described as in the first reference by Berger, the reverse projection
///  operator performs the following steps (see second reference by Antony):
///  1) Transform input from any quadrature node types to Gauss nodes
///  2) Apply mortar projection on Gauss nodes (identical to the original formulation in Berger)
///  3) Transform result from Gauss nodes to original quadrature node types
///
///  Therefore, the projection operator contructed here is a matrix multiplication of three
///  matrices. If the integration nodes are Gauss nodes as well, this is equivalent to the method in
///  the first reference.
inline void calcMortarProjectionMatrixHReverse(const MInt polyDeg,
                                               const MFloat* const nodes,
                                               const MFloat* const wBary,
                                               const MInt position,
                                               MFloat* const matrix) {
  const MInt noNodes = polyDeg + 1;

  // Calculate nodes and quadrature weights for projection (always Gauss)
  MFloatVector projectionNodes;
  MFloatVector projectionWInt;
  projectionNodes.resize(noNodes);
  projectionWInt.resize(noNodes);
  maia::dg::interpolation::calcLegendreGaussNodesAndWeights(polyDeg, &projectionNodes[0], &projectionWInt[0]);

  // Calculate barycentric weights for projection nodes
  MFloatVector projectionWBary;
  projectionWBary.resize(noNodes);
  maia::dg::interpolation::calcBarycentricWeights(polyDeg, &projectionNodes[0], &projectionWBary[0]);

  // Calculate projection matrix from mortar element to coarse element in projection nodes
  const MFloat shift = (position == dg::mortar::lower) ? -0.5 : 0.5;
  MFloatTensor mortarToCoarseElement(noNodes, noNodes);
  for(MInt i = 0; i < noNodes; i++) {
    maia::dg::interpolation::calcLagrangeInterpolatingPolynomials(0.5 * projectionNodes[i] + shift,
                                                                  polyDeg,
                                                                  &projectionNodes[0],
                                                                  &projectionWBary[0],
                                                                  &mortarToCoarseElement(i, 0));
  }

  // Transpose and multiply
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < i + 1; j++) {
      const MFloat temp = mortarToCoarseElement(i, j);
      mortarToCoarseElement(i, j) = 0.5 * mortarToCoarseElement(j, i) * projectionWInt[j] / projectionWInt[i];
      mortarToCoarseElement(j, i) = 0.5 * temp * projectionWInt[i] / projectionWInt[j];
    }
  }


  // Transformation matrix from grid to projection nodes (this is the identity when using Gauss
  // nodes anyways)
  MFloatTensor p2dg(noNodes, noNodes);
  calcMortarProjectionMatrixP(polyDeg, &projectionNodes[0], &projectionWBary[0], polyDeg, nodes, &p2dg[0]);

  // Transformation matrix from projection to grid nodes (this is the identity when using Gauss
  // nodes anyways)
  MFloatTensor dg2p(noNodes, noNodes);
  calcMortarProjectionMatrixP(polyDeg, nodes, wBary, polyDeg, &projectionNodes[0], &dg2p[0]);

  // The total projection matrix m is a product of the transformation matrices to(dg2p) and
  // from(p2dg) quadrature nodes and the projection matrix in quadrature nodes, in this order: m =
  // p2dg * mortarToCoarseElement * dg2p (where "*" denotes a matrix-matrix product)
  MFloatTensor m(matrix, noNodes, noNodes);
  m.set(0.0);
  for(MInt q = 0; q < noNodes; q++) {
    for(MInt j = 0; j < noNodes; j++) {
      for(MInt i = 0; i < noNodes; i++) {
        for(MInt n = 0; n < noNodes; n++) {
          m(q, j) += p2dg(q, i) * mortarToCoarseElement(i, n) * dg2p(n, j);
        }
      }
    }
  }
}

} // namespace mortar
} // namespace dg
} // namespace maia

#endif // define DGMORTAR_H_
